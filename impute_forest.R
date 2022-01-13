# This function runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

polygons <- 'D:/ontario_inventory/imputation/test_polygon_attributes/fri_lidar_vars.shp'
interp_vars <- c('HT', 'CC', 'BA')
lidar_vars <- c('lor', 'cc', 'ba')
percentile <- 0.3

impute_forest <- function(
  polygons, # FRI polygons shapefile with lidar attributes
  interp_vars, # character vector of column names of interpreter derived variables for data screening. 1 min/3 max. Order corresponds to lidar_vars.
  lidar_vars, # character vector of column names of lidar derived variables for data screening. 1 min/3 max. Order corresponds to interp_vars.
  percentile # fraction of data to screen from top and bottom of linear models
){

  # load polygon dataset
  poly <- terra::vect(polygons)
  
  # extract as dataframe
  dat <- terra::as.data.frame(poly)
  
  ######################
  ### DATA SCREENING ###
  ######################
  
  # first var
  # remove rows with missing values
  dat_a <- dat[is.na(dat[,colnames(dat) %in% interp_vars[1]]) == F & is.na(dat[,colnames(dat) %in% lidar_vars[1]]) == F,]
  
  # only keep land polygons
  dat_a <- dat_a[dat_a$POLYTYPE != 'WAT',]
  
  # run simple linear model
  lm_a <- lm(dat_a[,colnames(dat_a) %in% interp_vars[1]] ~ dat_a[,colnames(dat_a) %in% lidar_vars[1]])
  
  # bind residuals to data
  dat_a <- tibble::add_column(dat_a, resid = resid(lm_a))
  
  # find percentile
  lm_a_perc <- quantile(dat_a$resid, probs = c(percentile, 1-percentile))
  
  # remove percentile of data
  dat_a <- dat_a[dat_a$resid > lm_a_perc[1] & dat_a$resid < lm_a_perc[2],]
  
  # second var
  if(length(interp_vars) > 1 & length(lidar_vars) > 1){
    # remove rows with missing values
    dat_b <- dat[is.na(dat[,colnames(dat) %in% interp_vars[2]]) == F & is.na(dat[,colnames(dat) %in% lidar_vars[2]]) == F,]
    
    # only keep land polygons
    dat_b <- dat_b[dat_b$POLYTYPE != 'WAT',]
    
    # run simple linear model
    lm_b <- lm(dat_b[,colnames(dat_b) %in% interp_vars[2]] ~ dat_b[,colnames(dat_b) %in% lidar_vars[2]])
    
    # bind residuals to data
    dat_b <- tibble::add_column(dat_b, resid = resid(lm_b))
    
    # find percentile
    lm_b_perc <- quantile(dat_b$resid, probs = c(percentile, 1-percentile))
    
    # remove percentile of data
    dat_b <- dat_b[dat_b$resid > lm_b_perc[1] & dat_b$resid < lm_b_perc[2],]
  }
  
  # third var
  if(length(interp_vars) == 3 & length(lidar_vars) == 3){
    # remove rows with missing values
    dat_c <- dat[is.na(dat[,colnames(dat) %in% interp_vars[3]]) == F & is.na(dat[,colnames(dat) %in% lidar_vars[3]]) == F,]
    
    # only keep land polygons
    dat_c <- dat_c[dat_c$POLYTYPE != 'WAT',]
    
    # run simple linear model
    lm_c <- lm(dat_c[,colnames(dat_c) %in% interp_vars[3]] ~ dat_c[,colnames(dat_c) %in% lidar_vars[3]])
    
    # bind residuals to data
    dat_c <- tibble::add_column(dat_c, resid = resid(lm_c))
    
    # find percentile
    lm_c_perc <- quantile(dat_c$resid, probs = c(percentile, 1-percentile))
    
    # remove percentile of data
    dat_c <- dat_c[dat_c$resid > lm_c_perc[1] & dat_c$resid < lm_c_perc[2],]
  }
  
  # combine results
  dat_screen <- dat_a[,'POLYID']
  if(length(interp_vars) > 1 & length(lidar_vars) > 1){
    dat_screen <- intersect(dat_screen, dat_b[,'POLYID'])
  }
  if(length(interp_vars) == 3 & length(lidar_vars) == 3){
    dat_screen <- intersect(dat_screen, dat_c[,'POLYID'])
  }
  
}