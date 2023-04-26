library(sen2r)

sen2r(preprocess = T,
      s2_levels = c('11c', '12a'),
      timewindow = c("2018-06-01", "2018-09-30"),
      extent = 'D:/ontario_inventory/FSF/sentinel/extent_fsf.geojson',
      list_prods = 'BOA',
      list_indices = c('EVI', 'NBR', 'NDVI', 'SAVI'),
      index_source = 'BOA',
      mask_type = 'clear_sky',
      max_mask = 100,
      reference_path = "D:/ontario_inventory/FSF/ALS/zq95.img",
      res = c(20, 20),
      proj = "+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs",
      rgb_compression = 'DEFLATE',
      overwrite = T,
      path_l1c = 'D:/ontario_inventory/FSF/sentinel/l1c',
      path_l2a = 'D:/ontario_inventory/FSF/sentinel/l2a',
      path_out = "D:/ontario_inventory/FSF/sentinel/Processed",
      parallel = T)

sen2r(param_list = 'D:/ontario_inventory/FSF/sentinel/sen2r_parameters_fsf.json')
