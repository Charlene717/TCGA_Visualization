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

## Test Correlation analysis
library("ggpubr")
Test <- cor.test(GeneExp.df[1,] %>% as.numeric(), GeneExp.df[2,]%>% as.numeric(), method = c("pearson"))
Test_RValue <- Test[["estimate"]][["cor"]]
Test_PValue <- Test[["p.value"]]


ggscatter(GeneExp.df %>% t() %>% as.data.frame(), x = row.names(GeneExp.df)[1], y = row.names(GeneExp.df)[2],
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = row.names(GeneExp.df)[1], ylab = row.names(GeneExp.df)[2])


for (j in 1:length(Set_Target_geneset2)) {
  for (i in 1:length(Set_Target_geneset1)) {
    Temp <- cor.test(GeneExp.df[Set_Target_geneset2[j],] %>% as.numeric(),
                     GeneExp.df[Set_Target_geneset1[i],] %>% as.numeric(),
                     method = c("pearson"))

    COR_Rvalue.df[j,i] <- Temp[["estimate"]][["cor"]]
    COR_Pvalue.df[j,i] <- Temp[["p.value"]]

  }

}
rm(Temp,i,j)

## Check Result
ggscatter(GeneExp.df %>% t() %>% as.data.frame(), x = Set_Target_geneset2[1], y = Set_Target_geneset1[1],
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = Set_Target_geneset2[1], ylab = Set_Target_geneset1[1])



## Apply function



##### Plot Heatmap #####
## Ref: https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/#:~:text=35%20mins-,Hierarchical%20Clustering%20in%20R%3A%20The%20Essentials,clusters%20of%20samples%20and%20features.

library(ComplexHeatmap)
Heatmap(COR_Rvalue.df,
        name = "Correlation", #title of legend
        column_title = "Variables", row_title = "Samples",
        # row_names_gp = gpar(fontsize = 7) # Text size for row names
)



#### Export Result ####
## PDF


## TSV



#### Save RData ####

