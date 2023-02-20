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


##### Conditions setting* #####


##### Current path and new folder setting* #####




##### Import Data #####



##### Data preprocessing #####




##### Correlation analysis #####





##### Plot Heatmap #####



#### Export Result ####
## PDF


## TSV



#### Save RData ####

