# This code runs kNN imputation between FRI polygons and GRM output
# using best variable combination in order to generate results
# to use in imputation paper

# load packages
library(terra)
library(tidyverse)
library(magrittr)
library(RANN)
library(gtools)
library(janitor)
library(viridis)
library(gridExtra)
library(reshape2)

########################
### LOAD FRI DATASET ###
########################

# load extracted data frame
load('D:/ontario_inventory/dat/dat_fri_extr.RData')

# load photo interpreted polygons
poly_fri <-
  vect(
    'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp'
  )

# cbind centroids to dat
dat <- cbind(dat, centroids(poly_fri) %>% crds)

################################################
### SUBSET DATA BASED ON SCREENING PROCEDURE ###
################################################

# load final dataset after screening
# dom for, p95 > 5, cc > 50%
dat_screen <-
  read.csv(
    'D:/ontario_inventory/imputation/dat_screen_domfor_p95_cc.csv'
  )

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID, ]

# change age values to 2018 value
dat$AGE2018 <- 2018 - dat$YRORG

# select columns we need
dat %<>% select(POLYID, WG, SPCOMP,
                cc, avg, max, p95, 
                qav, ske, kur, cv, lor, 
                ba, qmdbh, dens, agb, 
                top_height, v, v_merch,
                B6, depth_q25, rumple, zentropy,
                zpcum8, zsd, x, y, AGE2018)

# remove any missing values
dat <- na.omit(dat)

############################
### ADD SPCOMP VARIABLES ###
############################

# look at first few rows
head(dat$SPCOMP)

# parse SPCOMP strings
sp <- str_split(dat$SPCOMP, pattern = "\\s{2}")

# add first species to dat
dat$SP1 <- sapply(sp, FUN = function(x){
  str <- x[1]
  str <- str_sub(str, start = 1, end = 2)
  return(str)
})

# add first species percent to dat
dat$SP1P <- sapply(sp, FUN = function(x){
  str <- x[2]
  if(is.na(str)){
    str <- 100
  } else{
    str <- str_sub(str, start = 1, end = 2)
  }
  return(str)
})

# add second species to dat
dat$SP2 <- sapply(sp, FUN = function(x){
  str <- x[2]
  if(is.na(str) == F){
    str <- str_sub(str, start = 3, end = 4)
  }
  return(str)
})

# add second species percent to dat
dat$SP2P <- sapply(sp, FUN = function(x){
  str <- x[3]
  if(is.na(str) == F){
    str <- str_sub(str, start = 1, end = 2)
  }
  return(str)
})

# load sp_group data
sp_group <- read.csv(
  'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest_SPGROUP.shp'
)

# change POLYID to numeric
dat %<>% mutate(POLYID = as.numeric(POLYID))
sp_group %<>% mutate(POLYID = as.numeric(POLYID))

# join to dat
dat <- left_join(dat, 
                 sp_group %>% select(POLYID, SpeciesGroup2, SpeciesGroup3), 
                 by = 'POLYID')

# change name to dat_fri
dat_fri <- dat
rm(dat)

###################################
### LOAD GRM DATASET AND SCREEN ###
###################################

# load GRM polygons
poly_grm <- vect('D:/ontario_inventory/segmentation/grm/shp/grm_10_01_05.shp')

# project to same as FRI polygons to make sure
poly_grm <- project(poly_grm, poly_fri)

# load GRM polygon data frame
load('D:/ontario_inventory/dat/dat_grm_extr.RData')

# cbind centroids to dat_grm
dat_grm <- cbind(dat_grm, centroids(poly_grm) %>% crds)

# screen for forested polygons
dat_grm_scr <- dat_grm %>% filter(dom_for == 'Yes',
                                  p95 >= 5,
                                  cc >= 50)

######################################################
### FUNCTIONS TO RUN K NEAREST NEIGHBOR IMPUTATION ###
######################################################

# create mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# create rmsd function
rmsd <- function(obs, est){
  sqrt(mean((est - obs) ^ 2))
}

# create rrmsd function
rrmsd <- function(obs, est){
  sqrt(mean((est - obs) ^ 2)) / mean(obs) * 100
}

# create md function
md <- function(obs, est){
  mean(est - obs)
}

# create rmd function
rmd <- function(obs, est){
  mean(est - obs) / mean(obs) * 100
}

# create mae function
mae <- function(obs, est){
  mean(abs(est - obs))
}

#################################
### CLEAN DATA FOR IMPUTATION ###
#################################

# create data frame for grm and fri metrics used in imputation
dat_grm_imp <- dat_grm_scr %>% select(id, avg, B6, 
                                      rumple, x,
                                      y, zpcum8, zsd) %>% na.omit
dat_fri_imp <- dat_fri %>% select(avg, B6, 
                                  rumple, x,
                                  y, zpcum8, zsd) %>% na.omit

# need to combine and scale all values together then separate again
dat_comb_scaled <- rbind(dat_grm_imp %>% select(-id), 
                         dat_fri_imp) %>% scale
dat_grm_scaled <- dat_comb_scaled[1:NROW(dat_grm_imp),]
dat_fri_scaled <- dat_comb_scaled[(NROW(dat_grm_imp)+1):(NROW(dat_grm_imp)+NROW(dat_fri_imp)),]

######################
### RUN IMPUTATION ###
######################

# create vector of imputation variables
vars <- c('avg', 'B6', 
          'rumple', 'x',
          'y', 'zpcum8', 'zsd')

# create vector of age and species variables
aux_vars <- c('AGE2018', 'WG', 'SP1', 'SP2',
              'SpeciesGroup3', 'SpeciesGroup2')

# run nearest neighbor imputation k = 5
nn <- nn2(dat_fri_scaled, dat_grm_scaled, k = 5)

# get nn indices
nni <- nn[[1]]

# add imputed model variables to GRM data frame
for(i in seq_along(vars)){
  dat_grm_imp %<>% add_column(
    !!str_c(vars[i], '_imp') := apply(
      nni, 
      MARGIN = 1, 
      FUN = function(x){
        mean(dat_fri_imp[x, vars[i]])
      }
    ))
}

# add age and species variables to GRM data frame
for(i in seq_along(aux_vars)){
  if(i == 'AGE2018'){
    dat_grm_imp %<>% add_column(
      !!str_c(aux_vars[i]) := apply(
        nni, 
        MARGIN = 1, 
        FUN = function(x){
          mean(dat_fri[x, aux_vars[i]])
        }
      ))
  } else{
    dat_grm_imp %<>% add_column(
      !!str_c(aux_vars[i]) := apply(
        nni, 
        MARGIN = 1, 
        FUN = function(x){
          getmode(dat_fri[x, aux_vars[i]])
        }
      ))
  }
}

# update colnames
dat_grm_imp %<>% rename(AGE = AGE2018,
                        CLASS3 = SpeciesGroup3,
                        CLASS5 = SpeciesGroup2)

# calculate performance across imputation attributes
for(i in seq_along(vars)){
  if(i == 1){
    perf <- tibble(!!str_c(vars[i], '_rmsd') := 
                     rmsd(dat_grm_imp[, vars[i]],
                         dat_grm_imp[, str_c(vars[i], '_imp')]),
                   !!str_c(vars[i], '_rrmsd') := 
                     rrmsd(dat_grm_imp[, vars[i]],
                         dat_grm_imp[, str_c(vars[i], '_imp')]),
                   !!str_c(vars[i], '_md') := 
                     md(dat_grm_imp[, vars[i]],
                           dat_grm_imp[, str_c(vars[i], '_imp')]),
                   !!str_c(vars[i], '_rmd') := 
                     rmd(dat_grm_imp[, vars[i]],
                           dat_grm_imp[, str_c(vars[i], '_imp')]),
                   !!str_c(vars[i], '_mae') := 
                     mae(dat_grm_imp[, vars[i]],
                           dat_grm_imp[, str_c(vars[i], '_imp')]))
  }else{
    perf %<>% add_column(!!str_c(vars[i], '_rmsd') := 
                           rmsd(dat_grm_imp[, vars[i]],
                                dat_grm_imp[, str_c(vars[i], '_imp')]),
                         !!str_c(vars[i], '_rrmsd') := 
                           rrmsd(dat_grm_imp[, vars[i]],
                                 dat_grm_imp[, str_c(vars[i], '_imp')]),
                         !!str_c(vars[i], '_md') := 
                           md(dat_grm_imp[, vars[i]],
                              dat_grm_imp[, str_c(vars[i], '_imp')]),
                         !!str_c(vars[i], '_rmd') := 
                           rmd(dat_grm_imp[, vars[i]],
                               dat_grm_imp[, str_c(vars[i], '_imp')]),
                         !!str_c(vars[i], '_mae') := 
                           mae(dat_grm_imp[, vars[i]],
                               dat_grm_imp[, str_c(vars[i], '_imp')]))
  }
}

# melt and round to two decimal places
perf %<>% melt %>% mutate(value = round(value, 2))

# add values back into main GRM data frame (missing values for polygons not
# included in the imputation)
dat_grm <- left_join(dat_grm, dat_grm_imp)

# add new data frame to polygons
values(poly_grm) <- dat_grm

# save grm polygon output
writeVector(poly_grm, 'D:/ontario_inventory/imputation/vanilla/grm_10_01_05_imp.shp', overwrite = T)

###############################################
### PREPARE FRI POLYGONS FOR FIGURE OUTPUTS ###
###############################################

# set age == 0 to NA
# we don't need to do this for AGE we will just
# do it for AGE2018
# poly_fri$AGE[poly_fri$AGE == 0] <- NA

# load whole FRI data frame again
dat <- as.data.frame(poly_fri)

# change POLYID to numeric
dat %<>% mutate(POLYID = as.numeric(POLYID))

# rename species groups in dat_fri
dat_fri %<>% rename(class3 = SpeciesGroup3, class5 = SpeciesGroup2)

# join species classes from dat_fri
# since those are the polygons we actually use in imputation
dat <- left_join(dat, 
                     dat_fri %>% select(POLYID, class3, class5, SP1, AGE2018), 
                     by = 'POLYID')

# re-input attributes into FRI polygons
values(poly_fri) <- dat
rm(dat)

# output shp files
writeVector(poly_fri[,c('AGE2018', 'class3', 'class5')],
            'D:/ontario_inventory/imputation/vanilla/fri_for_fig.shp', 
            overwrite = T)

writeVector(poly_grm[,c('AGE', 'CLASS3', 'CLASS5')],
            'D:/ontario_inventory/imputation/vanilla/grm_for_fig.shp', 
            overwrite = T)

# ##########################
# ### CREATE MAP OUTPUTS ###
# ##########################
# 
# # plot age
# png(file = 'C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/age_imp.png',
#     width = 12, height = 12, units = "in",
#     res = 300)
# 
# par(mfrow=c(2,1))
# 
# plot(poly_grm, 'AGE', col = viridis(14), type = 'interval', breaks = seq(0, 140, 10),
#      plg = list(x = 'topright', cex = 1.5, title = 'Age'),
#      main = '')
# title(main = 'Imputed Age of GRM Segmented Forest Stands', cex.main = 3)
# 
# plot(poly_fri, 'AGE2018', col = viridis(14), type = 'interval', breaks = seq(0, 140, 10),
#      plg = list(x = 'topright', cex = 1.5, title = 'Age'),
#      main = '')
# title(main = 'Age of FRI Forest Stands', cex.main = 3)
# 
# dev.off()
# 
# # plot group of 3 species
# png(file = 'C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/class3_imp.png',
#     width = 18, height = 12, units = 'in',
#     res = 300)
# 
# par(mfrow=c(2,1))
# 
# plot(poly_grm, 'CLASS3', col = viridis(3), type = 'classes',
#      plg = list(x = 'topright', cex = 2, title = 'Species Class'),
#      main = '')
# title(main = 'Imputed Species Class (3) of GRM Segmented Forest Stands', cex.main = 3)
# 
# plot(poly_fri, 'class3', col = viridis(3), type = 'classes',
#      plg = list(x = 'topright', cex = 2, title = 'Species Class'),
#      main = '')
# title(main = 'Species Class (3) FRI Forest Stands', cex.main = 3)
# 
# dev.off()
# 
# # plot group of 5 species
# par(mfrow=c(2,1))
# 
# plot(poly_grm, 'class5', col = viridis(5), type = 'classes',
#      plg = list(x = 'topright', cex = 2, title = 'Species Class'),
#      main = '')
# title(main = 'Imputed Species Class (5) of GRM Segmented Forest Stands', cex.main = 3)
# 
# plot(poly_fri, 'class5', col = viridis(5), type = 'classes',
#      plg = list(x = 'topright', cex = 2, title = 'Species Class'),
#      main = '')
# title(main = 'Species Class (5) FRI Forest Stands', cex.main = 3)
# 
# png(file = 'C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/class5_imp.png',
#     width = 18, height = 12, units = 'in')

#############################
### CREATE FIGURE OUTPUTS ###
#############################

# get median values of age
median(dat_grm$AGE, na.rm = T)
median(poly_fri$AGE2018, na.rm = T)

# density plots of age
p1 <- ggplot(dat_grm, aes(x = AGE)) +
  geom_density(fill = 'grey') +
  geom_vline(aes(xintercept = median(AGE, na.rm = T)),
             linetype = "dashed",
             size = 0.6) +
  xlim(c(0,200)) +
  ylim(c(0, 0.02)) +
  theme_bw() +
  xlab('Age') +
  ylab('Density') +
  ggtitle('Imputed Age of GRM Segmented \nForest Stands') +
  theme(text = element_text(size = 25),
        plot.title = element_text(size=30))

p2 <- ggplot(as.data.frame(poly_fri), aes(x = AGE2018)) +
  geom_density(fill = 'grey') +
  geom_vline(aes(xintercept = median(AGE2018, na.rm = T)),
             linetype = "dashed",
             size = 0.6) +
  xlim(c(0, 200)) +
  ylim(c(0, 0.02)) +
  theme_bw() +
  xlab('Age') +
  ylab('Density') +
  ggtitle('Age of FRI \nForest Stands') +
  theme(text = element_text(size = 25),
        plot.title = element_text(size=30))

g <- arrangeGrob(p1, p2, ncol = 2)

ggsave('C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/age_dens_imp.png',
       g, width = 18, height = 12)

# density plots of age for combined figure
ggplot(dat_grm, aes(x = AGE)) +
  geom_density(fill = 'grey') +
  geom_vline(aes(xintercept = median(AGE, na.rm = T)),
             linetype = "dashed",
             size = 0.6) +
  xlim(c(0,200)) +
  ylim(c(0, 0.02)) +
  theme_bw() +
  xlab('Age') +
  ylab('Density') +
  theme(text = element_text(size = 25),
        plot.title = element_text(size=30))

ggsave('C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/grm_dens_imp.png', 
       width = 6, height = 4)

ggplot(as.data.frame(poly_fri), aes(x = AGE2018)) +
  geom_density(fill = 'grey') +
  geom_vline(aes(xintercept = median(AGE2018, na.rm = T)),
             linetype = "dashed",
             size = 0.6) +
  xlim(c(0, 200)) +
  ylim(c(0, 0.02)) +
  theme_bw() +
  xlab('Age') +
  ylab('Density') +
  theme(text = element_text(size = 25),
        plot.title = element_text(size=30))

ggsave('C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/fri_dens_imp.png', 
       width = 6, height = 4)

# distribution of 3 species classes
# create data frame for GRM 3 classes
dat_grm_c3 <- dat_grm %>% 
  tabyl(CLASS3) %>%  
  filter(is.na(CLASS3) == F) %>% 
  arrange(desc(CLASS3)) %>%
  mutate(prop = n / sum(.$n)*100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop) %>%
  mutate(lbl = round(prop))

# create data frame for FRI 3 classes
dat_fri_c3 <- poly_fri %>% as.data.frame %>%
  tabyl(class3) %>%  
  filter(is.na(class3) == F) %>% 
  arrange(desc(class3)) %>%
  mutate(prop = n / sum(.$n)*100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop) %>%
  mutate(lbl = round(prop))

# plot
p1 <- ggplot(dat_grm_c3, aes(x = "", y = prop, fill = CLASS3)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 12) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30)) +
  scale_fill_manual(values = c('#228833', '#aa3377', '#ccbb44')) +
  labs(fill = "3 Species Classification") +
  ggtitle("Imputed 3 Species Classification \nof GRM Segmented Forest Stands")

p2 <- ggplot(dat_fri_c3, aes(x = "", y = prop, fill = class3)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 12) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30)) +
  scale_fill_manual(values = c('#228833', '#aa3377', '#ccbb44')) +
  labs(fill = "3 Species Classification") +
  ggtitle("3 Species Classification \nof FRI Forest Stands")

g <- arrangeGrob(p1, p2, ncol = 2)

ggsave('C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/class3_pie_imp.png',
       g, width = 18, height = 12)

# plot for combined figure
ggplot(dat_grm_c3, aes(x = "", y = prop, fill = CLASS3)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 10) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30),
        legend.position = 'none') +
  scale_fill_manual(values = c('#228833', '#aa3377', '#ccbb44')) +
  labs(fill = "3 Species Classification")

ggsave('C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/grm_class3_pie_imp.png',
       width = 6, height = 4)

ggplot(dat_fri_c3, aes(x = "", y = prop, fill = class3)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 10) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30),
        legend.position = 'none') +
  scale_fill_manual(values = c('#228833', '#aa3377', '#ccbb44')) +
  labs(fill = "3 Species Classification")

ggsave('C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/fri_class3_pie_imp.png',
       width = 6, height = 4)

# distribution of 5 species classes
# create data frame for GRM 5 classes
dat_grm_c5 <- dat_grm %>% 
  tabyl(CLASS5) %>%  
  filter(is.na(CLASS5) == F) %>% 
  arrange(desc(CLASS5)) %>%
  mutate(prop = n / sum(.$n)*100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop) %>%
  mutate(lbl = round(prop))

# create data frame for FRI 5 classes
dat_fri_c5 <- poly_fri %>% as.data.frame %>%
  tabyl(class5) %>%  
  filter(is.na(class5) == F) %>% 
  arrange(desc(class5)) %>%
  mutate(prop = n / sum(.$n)*100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop) %>%
  mutate(lbl = round(prop))

# plot
p1 <- ggplot(dat_grm_c5, aes(x = "", y = prop, fill = CLASS5)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 12) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30)) +
  scale_fill_manual(values = c('#ccbb44', '#228833', '#4477aa',
                                        '#ee6677', '#aa3377')) +
  labs(fill = "5 Species Classification") +
  ggtitle("Imputed 5 Species Classification \nof GRM Segmented Forest Stands")

p2 <- ggplot(dat_fri_c5, aes(x = "", y = prop, fill = class5)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 12) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30)) +
  scale_fill_manual(values = c('#ccbb44', '#228833', '#4477aa',
                                        '#ee6677', '#aa3377')) +
  labs(fill = "5 Species Classification") +
  ggtitle("5 Species Classification \nof FRI Forest Stands")

g <- arrangeGrob(p1, p2, ncol = 2)

ggsave('C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/class5_pie_imp.png',
       g, width = 18, height = 12)

# plot for combined figure
ggplot(dat_grm_c5, aes(x = "", y = prop, fill = CLASS5)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 10) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30),
        legend.position = 'none') +
  scale_fill_manual(values = c('#ccbb44', '#228833', '#4477aa',
                                        '#ee6677', '#aa3377')) +
                                          labs(fill = "5 Species Classification")

ggsave('C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/grm_class5_pie_imp.png',
       width = 6, height = 4)

ggplot(dat_fri_c5, aes(x = "", y = prop, fill = class5)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(y = ypos, label = str_c(lbl, "%")), size = 10) +
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, 'cm'),
        plot.title = element_text(size=30),
        legend.position = 'none') +
  scale_fill_manual(values = c('#ccbb44', '#228833', '#4477aa',
                                        '#ee6677', '#aa3377')) +
                                          labs(fill = "5 Species Classification")

ggsave('C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/figures/fri_class5_pie_imp.png',
       width = 6, height = 4)

################################################
### MAKE FINAL TABLE OF SPECIES DISTRIBUTION ###
################################################

# distribution of leading species class
# create data frame for GRM
dat_grm_sp1 <- dat_grm %>% 
  tabyl(SP1) %>%  
  filter(is.na(SP1) == F) %>% 
  arrange(SP1) %>%
  mutate(prop = n / sum(.$n)*100) %>%
  mutate(lbl = round(prop))

# create data frame for FRI
dat_fri_sp1 <- poly_fri %>% as.data.frame %>%
  tabyl(SP1) %>%  
  filter(is.na(SP1) == F) %>% 
  arrange(SP1) %>%
  mutate(prop = n / sum(.$n)*100) %>%
  mutate(lbl = round(prop))

# join leading species to make sure match
dat_sp1 <- dat_fri_sp1 %>% 
  left_join(dat_grm_sp1 %>% select(SP1, prop),
                                     by = 'SP1') %>%
  rename(prop_fri = prop.x,
         prop_grm = prop.y)

# rearrange in order
dat_grm_c3 %<>% arrange(CLASS3)
dat_fri_c3 %<>% arrange(class3)
dat_grm_c5 %<>% arrange(CLASS5)
dat_fri_c5 %<>% arrange(class5)

# add all sp groups
sp_dist <- tibble(class = dat_grm_c3$CLASS3, 
                  grm = dat_grm_c3$prop, 
                  fri = dat_fri_c3$prop) %>%
  add_row(class = dat_grm_c5$CLASS5,
          grm = dat_grm_c5$prop,
          fri = dat_fri_c5$prop) %>%
  add_row(class = dat_sp1$SP1,
          grm = dat_sp1$prop_grm,
          fri = dat_sp1$prop_fri)

# write csv
write.csv(sp_dist, 
          'C:/Users/bermane/Team Braintree Dropbox/Ethan Berman/UBC/Ontario Inventory/papers/paper 2/tables/sp_dist.csv',
          row.names = F)
