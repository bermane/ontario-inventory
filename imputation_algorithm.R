# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

# load packages
library(terra)
library(tidyverse)
library(yaImpute)
library(randomForest)

# example 1 from yaImpute user guide
data('iris')
set.seed(1)
refs <- sample(rownames(iris), 50)
x <- iris[, 1:3]
y <- iris[refs, 4:5]

ex <- yai(x = x, y = y, method = 'mahalanobis')
plot(ex)

head(impute(ex))

# example 2 from yaImpute user guide
data("MoscowMtStJoe")
x <- MoscowMtStJoe[, c("EASTING", "NORTHING", "ELEVMEAN",
                       "SLPMEAN", "ASPMEAN", "INTMEAN", "HTMEAN", "CCMEAN")]
x[, 5] <- (1 - cos((x[, 5] - 30) * pi/180))/2
names(x)[5] = "TrASP"
y <- MoscowMtStJoe[, c(1, 9, 12, 14, 18)]
mal <- yai(x = x, y = y, method = "mahalanobis")
msn <- yai(x = x, y = y, method = "msn")
gnn <- yai(x = x, y = y, method = "gnn")
ica <- yai(x = x, y = y, method = "ica")
d
y2 <- cbind(whatsMax(y[, 1:4]), y[, 5])
names(y2) <- c("MajorSpecies", "BasalAreaMajorSp", "TotalBA")
rf <- yai(x = x, y = y2, method = "randomForest")
head(y2)

plot(rf, vars = yvars(rf))

rfImp <- impute(rf)
rmsd <- compare.yai(mal, msn, gnn, rfImp, ica)
apply(rmsd, 2, mean, na.rm = TRUE)
plot(rmsd)

# load polygon attribute dataset (post screening)
dat <- read.csv('D:/ontario_inventory/imputation/dat_screen.csv')
  
# load polygon object for CRS purposes
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp')

# load LiDAR datasets we need to extract over polygons
# create named vector with variable names and data links

lidar <- c('cc' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_2m_cov.tif',
           'avg' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_avg.tif',
           'max' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_max.tif',
           'qav' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_qav.tif',
           'ske' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_ske.tif',
           'kur' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_kur.tif',
           'cv' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/SPL100 metrics/RMF_20m_T130cm_cv.tif',
           'lor' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_lor.tif',
           'ba' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_ba_ha.tif',
           'qmdbh' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_qmdbh.tif',
           'dens' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_dens.tif',
           'agb' = 'D:/ontario_inventory/romeo/RMF_EFI_layers/ABA layers SPL 2018/RMF_20m_T130cm_AGB_ha.tif')
