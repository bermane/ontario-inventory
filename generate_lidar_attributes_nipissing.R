# this code generates lidar gridded attributes from point cloud
# data in the nipissing forest
# data are already normalized to height above ground

# load packages
library(lidR)
library(terra)
library(tidyverse)
library(gdalUtils)

#Create lascatalog
las_infolder <- "D:/ontario_inventory/nipissing/als/NPS_B2J_norm"
ctg <- readLAScatalog(las_infolder, filter = "-drop_withheld -keep_random_fraction 0.25")

#Run in para
library(future)
plan(multisession, workers = 4L)

# #Produce rumple
# rumple <- rast(grid_metrics(ctg, ~lidRmetrics::metrics_rumple(X,Y,Z, pixel_size = 1), res = 20))
# writeRaster(rumple,'rumple.tif')

#Produce Canopy cover
#You can just keep changing the height for 5 and 15 m :
# Compute Canopy Cover > 2m ####
myMetrics = function(z,rn,height){
  first  = rn == 1L
  zfirst = z[first]
  nfirst = length(zfirst)
  firstabove2 = sum(zfirst > height)
  x = (firstabove2/nfirst)*100
  metrics = list(
    above2aboven1st = x, # Num of returns above 'height' divided by num of 1st returns
    zsqmean = sqrt(mean(z^2))  # Quadratic mean of z
  )
  metrics
  # Combined with standard metrics
  return( c(stdmetrics_z(z),metrics))
}

# option to write each chunk to disk
opt_output_files(ctg) <- "D:/ontario_inventory/nipissing/als/chunks/alsmetrics_{ORIGINALFILENAME}_{ID}"

# calc metrics
canCover = grid_metrics(ctg, ~myMetrics(Z, rn=ReturnNumber,2))

# check which chunks didn't process
# load vector of processed chunk ids
chunks <- list.files('D:/ontario_inventory/nipissing/als/chunks', pattern = '*.tif') %>%
  str_split('_')
chunks <- sapply(chunks, function(x) x[3])
chunks <- str_replace(chunks, '.tif', '') %>% as.numeric

# create seq from 1:length of ctg
seq <- 1:length(ctg$filename)

# find which files didn't get processed
ids <- which(!(seq %in% chunks))

# find which files didn't get processed
files <- ctg$filename[ids]

# set ctg to only process missing files
ctg$processed <- FALSE
ctg$processed[ids] <- TRUE

# rerun calc metrics
canCover = grid_metrics(ctg, ~myMetrics(Z, rn=ReturnNumber,2))

# build vrt
vrt <- vrt(list.files('D:/ontario_inventory/nipissing/als/chunks', pattern = '*.tif', full.names = T),
                    filename = 'D:/ontario_inventory/nipissing/als/chunks/vrt.vrt')

# check layer names
ras <- rast(list.files('D:/ontario_inventory/nipissing/als/chunks', pattern = '*.tif', full.names = T)[1])

# change vrt layer names
names(vrt) <- names(ras)

# Check files and generate vrt
plot(canCover$above2aboven1st, main = 'Haliburton Canopy Coverage at > 2m')
writeRaster(canCover$above2aboven1st,'CanopyCoverage2m.tif')

# # load las file
# las <- readLAS('D:/ontario_inventory/nipissing/als/NPS_B2J_norm/1kmZ175400515002020L.laz')
# 
# las <- readLAS('D:/ontario_inventory/nipissing/als/NPS_B2J_norm/1kmZ175400515002020L.laz', 
#                filter = '-keep_first')
# 
# # check the data
# las_check(las)
# 
# # plot the data
# plot(las)
# 
# # create a DTM with TIN method
# dtm <- rasterize_terrain(las, res = 1, algorithm = tin())
# 
# # DTM normalization
# nlas <- las - dtm # pretty fast?
# 
# hist(filter_ground(nlas)$Z, main = "", xlab = "Elevation")
# 
# # Point cloud normalization
# nlas2 <- normalize_height(las, tin())
#   
# hist(filter_ground(nlas2)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")
# 
# # define pixel metrics we want to calculate
# met <- function(x){
#   metrics = list(cv = sd(x) / mean(x))
#   return(c(metrics, stdmetrics))
# }
# 
# # calculate pixel metrics
# pmet <- pixel_metrics(nlas2, .stdmetrics, 20)
# 
# plot(hmean, col = height.colors(50))
