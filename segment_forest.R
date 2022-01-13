# This function applies mean shift image segmentation to 3-band LiDAR raster
# in order to create forest stand polygons

raster_in <- 'D:/ontario_inventory/segmentation/spl_stack.tif'
otb_dir <- 'C:/OTB/bin'
temp_dir <- 'D:/temp'
out_path <- 'D:/'
out_name <- 'ms_10_10_100'
spatial_r <- 10
range_r <- 10
min_size <- 100

segment_forest <- function(
  raster_in, # location of input 3-band LiDAR raster
  otb_dir, # location of otb directory
  temp_dir, # directory where OTB will output temporary files
  out_path, # location where mean shift vector file will be written
  out_name, # name of file to be written
  spatial_r = 10, # mean shift spatial radius
  range_r = 10, # mean shift range radius
  min_size = 100 #mean shift minimum polygon size
 ){
  
  # set working directory where temp files will be output
  setwd(temp_dir)
  
  # create function to run mean shift
  meanshift_otb <- function(otb_path = "", raster_in = "", out_path = "", name ="", spatialr = "10",
                            ranger = "10", minsize = "100", tilesizex = "500", tilesizey = "500",
                            outmode = "vector", cleanup = "true", ram = "256"){
    # Set configuration
    conf <- paste("-in", raster_in, "-spatialr", spatialr, "-ranger", ranger,
                  "-minsize", minsize, "-tilesizex", tilesizex, "-tilesizey", tilesizey,
                  "-mode", outmode, "-mode.vector.out", paste(out_path, "/", name, ".shp", sep=""),
                  "-cleanup", cleanup,"-ram", ram)
    
    # apply function in command line
    system(paste(otb_path, "/otbcli_LargeScaleMeanShift", " ", conf, sep=""))
    
    # save configuration for further use
    write.table(x = conf,file = paste(out_path,"/",name,"_conf.txt",sep=""),row.names = F, col.names = F)
  }
  
  # run mean shift
  meanshift_otb(otb_path = otb_dir,
                raster_in = raster_in,
                out_path = out_path,
                name = out_name,
                spatialr = as.character(spatial_r),
                ranger = as.character(range_r),
                minsize = as.character(min_size),
                ram = "1024")
  
}