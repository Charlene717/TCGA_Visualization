##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages #####
source("FUN_Package_InstLoad.R")
PKG_Basic.set <- c("tidyverse","circlize")
PKG_BiocManager.set <- c("ComplexHeatmap","limma")

FUN_Package_InstLoad(Basic.set = PKG_Basic.set, BiocManager.set = PKG_BiocManager.set)

if(!require(devtools)) install.packages("devtools")
if(!require(ggpubr)) devtools::install_github("kassambara/ggpubr")
library("ggpubr")

##### Function setting #####
source("FUN_Beautify_ggplot.R")

##### Import setting* #####
SetImportPath_FOL <- "Input_TCGA_PAAD"  # Input Folder Name
SetImport_RNAFileName <- "TCGA_PAAD_HiSeqV2" # RNAFileName <- "TCGA_PAAD_HiSeqV2"
# SetImport_CNVFileName <- "TCGA_Gistic2_CopyNumber_Gistic2_all_data_by_genes"
SetImport_PhenoFileName <- "TCGA.PAAD.sampleMap_PAAD_clinicalMatrix"


##### Conditions setting* #####
# Set_Target_geneset1 <- c("ARHGEF10L","RNF10","RNF11","RNF13","GTF2IP1","REM1","TSKS","ASS1")
Set_Target_geneset1 <- c("FN1", "MTDH", "CDH17", "VTN", "SPP1", "FGA", "FGB", "FGG","COL1A1", "COL1A2",
                         "COL4A1-A6", "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6",
                         "COL5A1", "COL5A2", "COL5A3",
                         "LAMA1-A5", "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5",
                         "LAMB1-B4", "LAMB1", "LAMB2", "LAMB3", "LAMB4",
                         "LAMC1-C3", "LAMC1", "LAMC2", "LAMC3", "AATK"
                         )
Set_Target_geneset1 <- alias2Symbol(Set_Target_geneset1, species = "Hs", expand.symbols = FALSE)


# Set_Target_geneset2 <- c("NCBP2","DISC1","RNF115","RNF112","SPN","DHX8","TCOF1","LRRTM3","NUP98")
Set_Target_geneset2 <- c("ITGA5", "ITGAV", "ITGA8", "ITGB1", "ITGB3", "ITGB5", "ITGB6", "ITGB8", "ITGA2B",
                         "ITGA1", "ITGA2", "ITGA10", "ITGA11",
                         "ITGA9", "ITGA4", "ITGB7", "ITGAE", "ITGAM", "ITGB2", "ITGAL", "ITGAD", "ITGAX",
                         "ITGA3", "ITGA6", "ITGA7", "ITGB4")
Set_Target_geneset2 <- alias2Symbol(Set_Target_geneset2, species = "Hs", expand.symbols = FALSE)

Set_col_Range  <- c(-1,-0.5, 0, 0.5,1)
Set_col  <- c("#1f5294", "#366cb3", "white", "#c44d75", "#ad2653")
Set_col_fun = colorRamp2(Set_col_Range, Set_col)


##### Export setting* #####
Set_Export_Name1 <- "TCGA"
Set_Export_Name2 <- "PAAD"
Set_Export_Name3 <- "PT_Grp1"

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

Plt.Barplot <- ggplot(Pheno.df, aes(x=as.factor(Pheno.df[,"sample_type"]), fill=as.factor(Pheno.df[,"sample_type"]))) + geom_bar()
Plt.Barplot

Plt.Barplot + labs(fill="sample_type", x="sample_type", y = "count")+
  theme_classic() %>% FUN_BeautifyggPlot(AxisTitleSize=2,LegPos = c(0.82, 0.85), OL_Thick = 1.5)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.65, vjust = 0.7),) -> Plt.Barplot1
Plt.Barplot1


## Extract Primary tumor
PT.set <- Pheno.df[Pheno.df$sample_type == "Primary Tumor" ,]$sampleID
GeneExp_Ori.df <- GeneExp.df
GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% PT.set]

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


##### Visualization #####
#### Plot Heatmap ####
## Heatmap in R: Static and Interactive Visualization
## Ref: https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/#:~:text=35%20mins-,Hierarchical%20Clustering%20in%20R%3A%20The%20Essentials,clusters%20of%20samples%20and%20features.
## Ref(ComplexHeatmap): https://jokergoo.github.io/ComplexHeatmap-reference/book/

library(ComplexHeatmap)
Heatmap(COR_Rvalue.df,
        name = "Correlation", #title of legend
        # column_title = "Variables", row_title = "Samples",
        col = Set_col_fun,
        # row_names_gp = gpar(fontsize = 7) # Text size for row names
) -> Plt.Heatmap

Plt.Heatmap

Heatmap(COR_Rvalue.df,
        name = "Correlation", #title of legend
        # column_title = "Variables", row_title = "Samples",
        col = Set_col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE
        # row_names_gp = gpar(fontsize = 7) # Text size for row names
) -> Plt.Heatmap_NonClu

Plt.Heatmap_NonClu

#### Correlation plot ####
## Ref: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

#### Dot plot ####
## By ChatGPT: https://openai.com/blog/chatgpt/
library(ggplot2)

# 將P值的dataframe轉換成長格式，並為每一列增加一個列名稱
# P_df_long <- reshape2::melt(COR_Pvalue.df, varnames = c("row", "col"), value.name = "P")

P_df_long <- reshape2::melt(data.frame(Gene= rownames(COR_Pvalue.df),-log10(COR_Pvalue.df)) ,
                            varnames = c("row", "col"), value.name = "log10P")

colnames(P_df_long)[1:2] <-c("Gene1","Gene2")

# 將R值的dataframe轉換成長格式，並為每一列增加一個列名稱
# R_df_long <- reshape2::melt(COR_Rvalue.df, varnames = c("row", "col"), value.name = "R")
R_df_long <- reshape2::melt(data.frame(Gene= rownames(COR_Rvalue.df),COR_Rvalue.df),
                            varnames = c("row", "col"), value.name = "R")

colnames(R_df_long)[1:2] <-c("Gene1","Gene2")

# 合併P值和R值的dataframe
df <- merge(P_df_long, R_df_long)

df$R <- df$R %>% as.numeric()

df$Gene1  <- factor(df$Gene1, levels=Set_Target_geneset2)
df$Gene2  <- factor(df$Gene2, levels=Set_Target_geneset1)

# 繪製氣泡圖
ggplot(df, aes(x = Gene2, y = Gene1, size = log10P, fill = R)) +
  geom_point(shape = 21) +
  scale_fill_gradientn(colours = Set_col, # c("#1f5294", "#366cb3", "white", "#c44d75", "#ad2653"),
                       #values = c(-1,-0.5, 0, 0.5, 1), # values = Set_col_Range
                       limits = c(-1, 1)) +
  scale_size_continuous(range = c(2, 10)) +
  labs(size = "-log10Pvalue", fill = "R value")+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9)) -> Plt.Dot
  # labs(title = "Bubble Plot", x = "Column", y = "Row", size = "P value", fill = "R value")
Plt.Dot+labs(x=NULL, y=NULL)


#### Export Result ####
## PDF
pdf(
  file = paste0(Result_Folder_Name,"/",Result_Folder_Name,"_Heatmap.pdf"),
  width = 12,  height = 10
)
Plt.Heatmap_NonClu
Plt.Heatmap
Plt.Barplot1
dev.off()

pdf(
  file = paste0(Result_Folder_Name,"/",Result_Folder_Name,"_DotPlot.pdf"),
  width = 12,  height = 10
)
Plt.Dot+labs(x=NULL, y=NULL)
Plt.Barplot1
dev.off()
## TSV
write.table(data.frame(Rvalue=row.names(COR_Rvalue.df),COR_Rvalue.df),
            file = paste0(Result_Folder_Name,"/",Result_Folder_Name,"_COR_Rvalue.tsv"),
            sep="\t", row.names= F, quote = FALSE)
write.table(data.frame(Pvalue=row.names(COR_Pvalue.df), COR_Pvalue.df),
            file = paste0(Result_Folder_Name,"/",Result_Folder_Name,"_COR_Pvalue.tsv"),
            sep="\t", row.names= F, quote = FALSE)


#### Save RData ####
save.image(paste0(Result_Folder_Name,"/",Result_Folder_Name,".RData"))
