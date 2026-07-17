# Tool-specific libraries (load as needed)
install.packages('SCINA') #SCINA
install.packages('scSorter') #scSorter

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R") #scType
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R") #scType
install.packages("Matrix") 
install.packages("openxlsx")
install.packages("HGNChelper")

system("git clone https://github.com/bioinfo-ibms-pumc/SCSA.git") #SCSA
system("cd SCSA")
system('pip install "pandas<2.0.0"')
system("pip install --upgrade statsmodels")
system(pip install --upgrade scikit-learn
pip install --upgrade xgboost)

install.packages("scCATCH") #scCATCH

devtools::install_github("Irrationone/cellassign") #cellassign
reticulate::py_install("tensorflow-probability[tf]", pip = TRUE)
reticulate::py_install("tensorflow", pip = TRUE)

devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3") #garnett
devtools::install_github('cole-trapnell-lab/monocle3')
BiocManager::install(c('DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))
devtools::install_github('satijalab/seurat-wrappers')

BiocManager::install("SingleR") #SingleR

BiocManager::install("scmap") #scmap

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("zwj-tina/scibetR") #SciBetR

BiocManager::install("rnabioco/clustifyr") #clustifyr

devtools::install_github("immunogenomics/harmony") #scpred
devtools::install_github(repo="powellgenomicslab/scPred",  ref="9f407b7436f40d44224a5976a94cc6815c6e837f")

install.packages("xgboost") 

devtools::install_github("souravc83/fastAdaboost")


BiocManager::install('CHETAH') #CHETAH

BiocManager::install(c("S4Vectors", "hopach", "limma")) #scClassify 
devtools::install_github("SydneyBioX/scClassify")

BiocManager::install("scAnnotatR") #scAnnotatR

devtools::install_github("pcahan1/singleCellNet") #singleCellNet

# setup_celltypist.R is not included in this archive; it built a conda env for CellTypist.
# See classic-ML-based/celltypist_helper.py for how the tool itself invokes that environment.
# setup_celltypist_environment()
# verify_celltypist_installation()

devtools::install_github("pcahan1/singleCellNet") #singleCellNet

install.packages("scAnnotate") #scAnnotate
BiocManager::install('glmGamPoi')

devtools::install_github("BatadaLab/scID", force = TRUE) #scID
BiocManager::install("MAST")

install.packages("mLLMCelltype") #mLLMCelltype

devtools::install_github("ElliotXie/CASSIA/CASSIA_R", ref = "install", force = TRUE)
library(CASSIA)

# Automatically set up the Python environment if needed
setup_cassia_env()

install.packages("openai")
remotes::install_github("Winnie09/GPTCelltype") #GPTCelltype

devtools::install_github("liuhong-jia/scAnno")  #scAnno

library(devtools)
library(SingleCellExperiment)
library(M3Drop)
install_github("bm2-lab/scLearn") #scLearn

remotes::install_github("haoharryfeng/NeuCA") #NeuCA


install_github("ziyili20/CAMLU", build_vignettes = FALSE) #CAMLU
install.packages("keras")


devtools::install_github("atakanekiz/CIPR-Package", build_vignettes = FALSE) #CIPR

devtools::install_github("swainasish/ScInfeR") #ScInfeR


BiocManager::install(c(
  "ontologyIndex", "ontologySimilarity", "ontoProc", # Bioconductor
  "stringdist", "stringr", "dplyr", "irr", "caret", "MLmetrics" # CRAN
))


SCLEARN 

install_github("bm2-lab/scLearn")

