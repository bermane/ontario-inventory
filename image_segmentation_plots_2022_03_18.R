# This code generates plots of polygon distribution from segmentation

# load packages
library(terra)
library(tidyverse)
library(exactextractr)
library(sf)
library(janitor)


# dataset parameters
name_o <- 'D:/ontario_inventory/segmentation/mask_wat_ucl/ms_10_10_100_add_pt.shp'
name_out <- 'ms_10_10_100'
name <- 'LiDAR Polygons'

# extract polygon stats from qgis meanshift test runs
# load polygon and remove polygons with 1 pixel (eventually need to merge these maybe?)
# need a spatially contiguous file
p <- vect(name_o) %>% as.data.frame

# subset reliable forested polygons
p2 <- p %>% filter(dom_for == 'Yes', is.na(POLYTYPE), p95 >= 5)

# subset non masked WAT and UCL polygons
p3 <- p %>% filter(is.na(POLYTYPE))

# load interpreter derived polygons to extract statistics
poly <- vect('D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp') %>%
  as.data.frame

# subset forested polygons only
poly_for <- poly[poly$POLYTYPE=='FOR',]

# subset all non water/ucl polygons
poly_land <- poly[!(poly$POLYTYPE %in% c('WAT', 'UCL')),]

# combine datasets
poly_combine <- rbind(data.frame(ha = (poly_land$AREA)/10000,
                           lab = 'FRI Polygons'),
                      data.frame(ha = (p3$nbPixels)*0.04,
                                 lab = 'LiDAR Polygons'))

# calculate mean
mu <- poly_combine %>% 
  group_by(lab) %>%
  summarise(mean = mean(ha))

# plot density
ggplot(data.frame(nbPixels = p3$nbPixels), aes(x = nbPixels*0.04)) +
  geom_density() +
  xlim(c(0,100)) +
  ylim(c(0, 0.2)) +
  geom_vline(aes(xintercept = median(nbPixels*0.04, na.rm = T)), 
             linetype = "dashed", size = 0.6) +
  theme_bw() +
  xlab('Number of Hectares') +
  ylab('Density') +
  ggtitle(str_c(name, ' (Non WAT/UCL)')) +
  theme(text = element_text(size = 20))

# save plot
ggsave(str_c('D:/ontario_inventory/segmentation/mask_wat_ucl/plots/', name_out, '_mask_wat_ucl.png'),
       width = 2100, height = 2100, units = 'px')

# plot histogram
ggplot(poly_combine, aes(x = ha, color = lab)) +
  geom_histogram(binwidth = 1, fill = 'white', alpha = 0.3, position = 'identity') +
  xlim(c(0,50)) +
  theme_bw() +
  xlab('Number of Hectares') +
  ylab('Number of Polygons') +
  ggtitle('Polygon Distribution (non wat/ucl)') +
  theme(text = element_text(size = 16)) +
  scale_y_continuous(breaks = c(0, 2500, 5000, 7500, 10000)) +
  geom_vline(data = mu, aes(xintercept = mean, color = lab),
             linetype="dashed") +
  scale_color_manual(values = c('FRI Polygons' = '#EE7733',
                                'LiDAR Polygons' = '#0077BB'),
                     name = 'Dataset') +
  annotate(geom = 'text', x = 14, y = 10500, label = 'mean lines')
  
# save plot
ggsave(str_c('D:/ontario_inventory/segmentation/mask_wat_ucl/plots/fri_lidar_polygons_mask_wat_ucl.png'),
       width = 2100, height = 1369, units = 'px')

# plot histogram
ggplot(poly_combine, aes(x = ha, color = lab)) +
  geom_histogram(binwidth = 1, fill = 'white', position = 'dodge') +
  xlim(c(0,50)) +
  theme_bw() +
  xlab('Number of Hectares') +
  ylab('Number of Polygons') +
  ggtitle('Polygon Distribution (non wat/ucl)') +
  theme(text = element_text(size = 20)) +
  geom_vline(data = mu, aes(xintercept = mean, color = lab),
             linetype="dashed")

 

# # plot histogram
# ggplot(data.frame(nbPixels = p3$nbPixels), aes(x = nbPixels*0.04)) +
#   geom_histogram(binwidth = 1) +
#   xlim(c(0,100))+
#   geom_vline(aes(xintercept = median(nbPixels*0.04, na.rm = T)), 
#              linetype = "dashed", size = 0.6) +
#   theme_bw() +
#   xlab('Number of Hectares') +
#   ylab('Number of Polygons') +
#   ggtitle(str_c(name, ' (Non WAT/UCL)')) +
#   theme(text = element_text(size = 20))

