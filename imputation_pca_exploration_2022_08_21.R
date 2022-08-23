# load packages
library(terra)
library(tidyverse)
library(yaImpute)
library(randomForest)
library(viridis)
library(factoextra)

###################################################
###LOAD POLYGON DATASET AND ADD LIDAR ATTRIBUTES###
###################################################

# load extracted data frame
# dat_fri_100 is loaded with LiDAR arributes derived from 
# masking pixels that do not fall 100 % inside polygons
load('D:/ontario_inventory/dat/dat_fri_100.RData')

# # save full original dataset
# dat_orig <- dat
# 
# # only keep forested polygons
# dat_orig <- filter(dat_orig, POLYTYPE == 'FOR')

# load photo interpreted polygons
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <- read.csv('D:/ontario_inventory/imputation/fri_data_screen_1bc_2a_2p95_2cc_10perc_2022_07_06.csv')

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID,]

# create subset for PCA analysis with only ALS p95, ba, and cc
dat_pca <- dat %>% select(p95, ba, cc) %>% na.omit

# run pca
pca <- prcomp(dat_pca, scale. = T)

# plot results
fviz_eig(pca)

fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# calculate cumulative proportion of variance explained
std <- pca$sdev
var <- std^2
varex <- var/sum(var)

plot(cumsum(varex), xlab = "Principal Component",ylab = "Cumulative Proportion of Variance Explained",type = "b")
abline(h=0.975,col='red',v=30)
