# this code generates lidar gridded attributes from point cloud
# data in the nipissing forest

# load packages
library(lidR)


# load las file
las <- readLAS('D:/ontario_inventory/nipissing/als/NPS_B2J_norm/1kmZ175400515002020L.laz')

las <- readLAS('D:/ontario_inventory/nipissing/als/NPS_B2J_norm/1kmZ175400515002020L.laz', 
               filter = '-keep_first')

# check the data
las_check(las)

# plot the data
plot(las)

# create a DTM with TIN method
dtm <- rasterize_terrain(las, res = 1, algorithm = tin())

# DTM normalization
nlas <- las - dtm # pretty fast?

hist(filter_ground(nlas)$Z, main = "", xlab = "Elevation")

# Point cloud normalization
nlas2 <- normalize_height(las, tin())
  
hist(filter_ground(nlas2)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")

# define pixel metrics we want to calculate
met <- function(x){
  metrics = list(cv = sd(x) / mean(x))
  return(c(metrics, stdmetrics))
}

# calculate pixel metrics
pmet <- pixel_metrics(nlas2, .stdmetrics, 20)

plot(hmean, col = height.colors(50))
