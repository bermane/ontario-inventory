# check correlation of ALS variables that could have
# been used for segmentation
# we arbitrarily selected p95, cc, cv

# load packages
library(terra)
library(tidyverse)
library(sf)
library(magrittr)
library(corrplot)

# load dataset with extracted als attributes
load('D:/ontario_inventory/dat/dat_fri_extr.RData')

# remove water and unclassified polygons
dat %<>% filter(!(POLYTYPE %in% c('WAT', 'UCL')))

# select ALS variables
vars <- c('cc', 'avg', 'max', 'p95',
          'qav', 'ske', 'kur', 'cv')

# look at correlation of variables
cor <- cor(dat[, vars], use = 'complete.obs')
corrplot(cor, method = 'number')
