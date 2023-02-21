##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages #####
source("FUN_Package_InstLoad.R")
PKG_Basic.set <- c("tidyverse","circlize")
PKG_BiocManager.set <- c("ComplexHeatmap")

FUN_Package_InstLoad(Basic.set = PKG_Basic.set, BiocManager.set = PKG_BiocManager.set)

if(!require(devtools)) install.packages("devtools")
if(!require(ggpubr)) devtools::install_github("kassambara/ggpubr")
library("ggpubr")

##### Function setting #####
source("FUN_Beautify_ggplot.R")
source("FUN_ggPlot_vline.R")


##### Import setting* #####
SetImportPath_FOL <- "Input_TCGA_PAAD"  # Input Folder Name
SetImport_RNAFileName <- "TCGA_PAAD_HiSeqV2" # RNAFileName <- "TCGA_PAAD_HiSeqV2"
# SetImport_CNVFileName <- "TCGA_Gistic2_CopyNumber_Gistic2_all_data_by_genes"
SetImport_PhenoFileName <- "TCGA-PAAD.GDC_phenotype.tsv"


##### Conditions setting* #####
Set_Target_geneset1 <- c("ARHGEF10L","RNF10","RNF11","RNF13","GTF2IP1","REM1","TSKS","ASS1")
Set_Target_geneset2 <- c("NCBP2","DISC1","RNF115","RNF112","SPN","DHX8","TCOF1","LRRTM3","NUP98")

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
Pheno.df <- read.delim(paste0(SetImportPath_FOL,"/",SetImport_PhenoFileName), header=T,  sep="\t") # row.names = 1,


##### Data preprocessing #####
## Extract Primary tumor



##### Correlation analysis #####
## Ref: http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
## Ref: https://www.scribbr.com/statistics/pearson-correlation-coefficient/#:~:text=Revised%20on%20December%205%2C%202022,the%20relationship%20between%20two%20variables.


## For loop
COR_Rvalue.df <- as.data.frame(matrix(data = NA, ncol = length(Set_Target_geneset1), nrow = length(Set_Target_geneset2)))
colnames(COR_Rvalue.df) <- Set_Target_geneset1
row.names(COR_Rvalue.df) <- Set_Target_geneset2

COR_Pvalue.df <- COR_Rvalue.df

TTT <- cor.test(GeneExp.df[1,] %>% as.numeric(), GeneExp.df[2,]%>% as.numeric(), method = c("pearson"))

for (i in 1:length(Set_Target_geneset1)) {


  for (j in 1:length(Set_Target_geneset2)) {


  }

}



## Apply function



##### Plot Heatmap #####



#### Export Result ####
## PDF


## TSV



#### Save RData ####

