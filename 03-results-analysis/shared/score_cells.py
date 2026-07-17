import os
import re
import csv
import json
import argparse
from pathlib import Path
from openai import OpenAI

# Shared by both the database-marker-tools (Phase 4, was P7a) and llm-ontology-scoring
# (Phase 4, was P7b) phases -- the two original per-phase copies were identical apart from
# these paths, so pass --raw-dir/--out-dir/--checkpoint-prefix to point at your own phase's
# raw_data/with_score directories (not included in this archive).
SKIP_FILES = {"runtime_dataset.csv"}
BATCH_SIZE = 25

SYSTEM_PROMPT = """You are a cell biology expert scoring predicted vs ground-truth cell types.
Score each pair as exactly one of: 1, 0.5, or 0
  1   = same type (synonyms, plurals, minor naming differences, abbreviations)
  0.5 = same broad lineage/tissue but different specificity
  0   = different cell type or biologically unrelated

Key rules:
- Know common scRNA-seq abbreviations: OPC/COP/NFOL/MFOL/MOL=oligodendrocyte lineage,
  IN-CTX=cortex interneuron, IPC=intermediate progenitor, NK=natural killer, RG=radial glia
- Kupffer cells = liver macrophages -> score 1 against macrophage predictions
- Keratinocytes = epithelial -> score 1 against epithelial predictions
- "CD13+ sorted cells" style labels are FACS gates, score on biological plausibility
- Subtype vs parent: "Macrophage C1QB" vs "Macrophages" = 1
- Wrong lineage entirely = 0; related but distinct = 0.5

Return ONLY a JSON array of numbers, one per pair, in the same order as the input.
Example for 3 pairs: [1, 0.5, 0]"""


def strip_numbering(s):
    return re.sub(r"^\d+\.\s*", "", s).strip()


def load_checkpoint(path):
    if os.path.exists(path):
        with open(path) as f:
            raw = json.load(f)
        return {tuple(json.loads(k)): v for k, v in raw.items()}
    return {}


def save_checkpoint(scores, path):
    serializable = {json.dumps(list(k)): v for k, v in scores.items()}
    with open(path, "w") as f:
        json.dump(serializable, f)


def collect_all_pairs(raw_dir):
    unique_pairs = set()
    file_rows = {}
    for fname in sorted(os.listdir(raw_dir)):
        if not fname.endswith(".csv") or fname in SKIP_FILES:
            continue
        fpath = os.path.join(raw_dir, fname)
        with open(fpath, encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f, delimiter=",")
            rows_data = []
            for row in reader:
                gt = row.get("ground_truth_celltype", "").strip()
                pred = row.get("raw_pred", "").strip()
                rows_data.append((gt, pred))
                unique_pairs.add((gt, pred))
        file_rows[fpath] = rows_data
    return unique_pairs, file_rows


def _call_api(client, model, pairs_list, service_tier=None):
    """Call the API for a list of pairs and return a scores dict. Raises on failure."""
    n = len(pairs_list)
    pairs_json = json.dumps(
        [{"i": i, "gt": gt, "pred": pred} for i, (gt, pred) in enumerate(pairs_list)],
        ensure_ascii=False,
    )
    user_msg = (
        f"{pairs_json}\n\n"
        f"Score exactly {n} pairs. Return exactly {n} numbers as a JSON array."
    )
    token_kwarg = (
        {"max_completion_tokens": 512}
        if model.startswith("gpt-5")
        else {"max_tokens": 512}
    )
    extra_kwargs = {}
    if service_tier is not None:
        extra_kwargs["service_tier"] = service_tier
    response = client.chat.completions.create(
        model=model,
        **token_kwarg,
        **extra_kwargs,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": user_msg},
        ],
    )
    raw_text = response.choices[0].message.content.strip()
    raw_text = re.sub(r"^```(?:json)?\s*", "", raw_text)
    raw_text = re.sub(r"\s*```$", "", raw_text)
    scores_list = json.loads(raw_text)
    if len(scores_list) != len(pairs_list):
        raise ValueError(f"Expected {len(pairs_list)} scores, got {len(scores_list)}")
    result = {}
    for (gt, pred), raw_val in zip(pairs_list, scores_list):
        val = float(raw_val)
        val = min([0, 0.5, 1], key=lambda x: abs(x - val))
        result[(gt, pred)] = val
    return result


def score_batch(client, model, pairs_list, service_tier=None):
    """Score a batch with one retry: on failure split in half and retry each half."""
    try:
        return _call_api(client, model, pairs_list, service_tier=service_tier)
    except Exception as e:
        if len(pairs_list) == 1:
            raise
        print(f"    Batch failed ({e}), splitting and retrying...")
        mid = len(pairs_list) // 2
        result = {}
        result.update(_call_api(client, model, pairs_list[:mid], service_tier=service_tier))
        result.update(_call_api(client, model, pairs_list[mid:], service_tier=service_tier))
        return result


def run_scoring_pass(client, model, pairs_list, checkpoint_path, label, service_tier=None):
    all_scores = load_checkpoint(checkpoint_path)
    if all_scores:
        print(f"  [{label}] Resumed from checkpoint: {len(all_scores)} pairs already scored")

    remaining = [p for p in pairs_list if p not in all_scores]
    total_batches = (len(remaining) + BATCH_SIZE - 1) // BATCH_SIZE

    for i in range(0, len(remaining), BATCH_SIZE):
        batch = remaining[i : i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        print(f"  [{label}] Batch {batch_num}/{total_batches} ({len(batch)} pairs)...")
        try:
            batch_scores = score_batch(client, model, batch, service_tier=service_tier)
            all_scores.update(batch_scores)
            save_checkpoint(all_scores, checkpoint_path)
        except Exception as e:
            print(f"  [{label}] ERROR on batch {batch_num}: {e}")
            print("  Saving progress and continuing...")
            save_checkpoint(all_scores, checkpoint_path)

    return all_scores


def apply_scores_to_files(file_rows, scores1, scores2, scores3, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for fpath, _ in file_rows.items():
        with open(fpath, encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f, delimiter=",")
            fieldnames = list(reader.fieldnames)
            rows = list(reader)

        new_fieldnames = list(fieldnames)
        for col in ("score_1", "score_2", "score_3"):
            if col not in new_fieldnames:
                new_fieldnames.append(col)

        for row in rows:
            gt = row.get("ground_truth_celltype", "").strip()
            pred = row.get("raw_pred", "").strip()

            if not pred or pred.upper() == "NA":
                row["score_1"] = 0
                row["score_2"] = 0
                row["score_3"] = 0
                continue

            key = (gt, strip_numbering(pred))
            key_orig = (gt, pred)

            s1 = scores1.get(key, scores1.get(key_orig, ""))
            s2 = scores2.get(key, scores2.get(key_orig, ""))
            s3 = scores3.get(key, scores3.get(key_orig, ""))

            if s1 == "":
                print(f"WARNING: no score_1 for {key} in {os.path.basename(fpath)}")
            if s2 == "":
                print(f"WARNING: no score_2 for {key} in {os.path.basename(fpath)}")
            if s3 == "":
                print(f"WARNING: no score_3 for {key} in {os.path.basename(fpath)}")

            row["score_1"] = s1
            row["score_2"] = s2
            row["score_3"] = s3

        stem = Path(fpath).stem
        out_path = os.path.join(out_dir, stem + ".tsv")
        with open(out_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=new_fieldnames,
                delimiter="\t",
                quoting=csv.QUOTE_MINIMAL,
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)

    print(f"  Written {len(file_rows)} files.")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", required=True, help="Directory of raw *_pred.csv files to score")
    parser.add_argument("--out-dir", required=True, help="Directory to write *_score.tsv files to")
    parser.add_argument(
        "--checkpoint-prefix",
        default="/tmp/score_cells",
        help="Prefix for the three resumable checkpoint JSON files (default: /tmp/score_cells)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    checkpoint_1 = f"{args.checkpoint_prefix}_1_checkpoint.json"
    checkpoint_2 = f"{args.checkpoint_prefix}_2_checkpoint.json"
    checkpoint_3 = f"{args.checkpoint_prefix}_3_checkpoint.json"

    openrouter_client = OpenAI(
        base_url="https://openrouter.ai/api/v1",
        api_key=os.environ["OPENROUTER_API_KEY"],
    )
    openai_client = OpenAI(
        api_key=os.environ["OPENAI_API_KEY"],
    )

    print("Step 1: Collecting all unique pairs...")
    unique_pairs, file_rows = collect_all_pairs(args.raw_dir)
    print(f"  Total files: {len(file_rows)}")
    print(f"  Total unique pairs: {len(unique_pairs)}")

    na_pairs = {(gt, pred) for gt, pred in unique_pairs if not pred or pred.upper() == "NA"}
    api_pairs = {
        (gt, strip_numbering(pred))
        for gt, pred in unique_pairs
        if pred and pred.upper() != "NA"
    }
    print(f"  NA predictions (auto-score 0): {len(na_pairs)}")
    print(f"  Pairs needing API scoring: {len(api_pairs)}")

    api_pairs_list = list(api_pairs)

    print("\nStep 2: Scoring via gemini-3-flash-preview (score_1)...")
    scores1 = run_scoring_pass(
        openrouter_client,
        "google/gemini-3-flash-preview",
        api_pairs_list,
        checkpoint_1,
        "score_1",
        service_tier="flex",
    )

    print("\nStep 3: Scoring via claude-haiku-4-5 (score_2)...")
    scores2 = run_scoring_pass(
        openrouter_client,
        "anthropic/claude-haiku-4-5",
        api_pairs_list,
        checkpoint_2,
        "score_2",
        service_tier="flex",
    )

    print("\nStep 4: Scoring via gpt-5.4-nano-2026-03-17 (score_3)...")
    scores3 = run_scoring_pass(
        openai_client,
        "gpt-5.4-nano-2026-03-17",
        api_pairs_list,
        checkpoint_3,
        "score_3",
    )

    print(f"\nStep 5: Writing score_1/score_2/score_3 to {len(file_rows)} files...")
    apply_scores_to_files(file_rows, scores1, scores2, scores3, args.out_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
