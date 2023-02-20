##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages #####
source("FUN_Package_InstLoad.R")
PKG_Basic.set <- c("tidyverse","circlize")
PKG_BiocManager.set <- c("ComplexHeatmap")

FUN_Package_InstLoad(Basic.set = PKG_Basic.set, BiocManager.set = PKG_BiocManager.set)

##### Function setting #####
source("FUN_Beautify_ggplot.R")
source("FUN_ggPlot_vline.R")


##### Import setting* #####
SetImportPath_FOL <- "Input_TCGA_PAAD"  # Input Folder Name
SetImport_RNAFileName <- "TCGA_PAAD_HiSeqV2" # RNAFileName <- "TCGA_PAAD_HiSeqV2"
# SetImport_CNVFileName <- "TCGA_Gistic2_CopyNumber_Gistic2_all_data_by_genes"
SetImport_PhenoFileName <- "TCGA-PAAD.GDC_phenotype.tsv"


##### Conditions setting* #####

Set_Target_geneset <-


##### Export setting* #####
Set_Export_Name1 <- "TCGA"
Set_Export_Name2 <- "PAAD"
Set_Export_Name3 <- "Test"

Result_Folder_Name <- paste0("Result_",Sys.Date(),"_",Set_Export_Name1,"_",Set_Export_Name2,"_",Set_Export_Name3) ## Generate output folder automatically
dir.create(Result_Folder_Name)


##### Import Data #####
## Import RNA expression data
GeneExp.df <- read.table(paste0(SetImportPath_FOL,"/",SetImport_RNAFileName),
                         header=T, row.names = 1, sep="\t")
colnames(GeneExp.df) <-  gsub("\\.", "-", colnames(GeneExp.df))

## Import Metadata
Pheno.df <- read.delim(paste0(SetImportPath_FOL,"/",PhenoFileName), header=T,  sep="\t") # row.names = 1,


##### Data preprocessing #####




##### Correlation analysis #####





##### Plot Heatmap #####



#### Export Result ####
## PDF


## TSV



#### Save RData ####

