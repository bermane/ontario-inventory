# This code runs an imputation algorithm on LiDAR attributes connecting them to interpreter
# derived response attributes resulting in a model that can be applied to forest stands
# segmented from LiDAR rasters

# The purpose is to compare the distributions of various attributes (HT, CC, BA, AGE, POLYTYPE, WG)
# from FRI polygons, LiDAR attributes in FRI polygons, imputed vars into FRI polys, imputed vars into LiDAR polys

# this new version only runs forested pixels since we are only really interested in modelling
# and imputing attributes in forests
# it uses data screened for 10 and 90th percentile of residuals

# This version of the code runs a bootstrap out of bag validation

# load packages
library(terra)
library(tidyverse)
library(yaImpute)
library(randomForest)
library(viridis)
library(janitor)

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
poly <-
  vect(
    'D:/ontario_inventory/romeo/RMF_EFI_layers/Polygons Inventory/RMF_PolygonForest.shp'
  )

##############################################
###SUBSET DATA BASED ON SCREENING PROCEDURE###
##############################################

# load final dataset after screening
dat_screen <-
  read.csv(
    'D:/ontario_inventory/imputation/fri_data_screen_1bc_2a_2p95_2cc_10perc_2022_07_06.csv'
  )

# subset dat based on screening
dat <- dat[dat$FOREST_ID %in% dat_screen$FOREST_ID, ]

# change all non-numeric variables to factor
dat[sapply(dat, is.character)] <-
  lapply(dat[sapply(dat, is.character)],
         as.factor)

##################
###EXPLORE DATA###
##################

# # lets check the most prominent working groups
# # subset data based on top 5
# wg <- tabyl(dat$WG)
# wg <- wg[order(-wg$n),]
# wg <- wg$`dat$WG`[1:5]
# dat_wg <- dat[dat$WG %in% wg,]
# 
# # p95 vs wg
# ggplot(dat_wg, aes(x = WG, y = p95)) +
#   geom_boxplot(fill='#A4A4A4', color="black")+
#   theme_classic()
# 
# # cc vs wg
# ggplot(dat_wg, aes(x = WG, y = cc)) +
#   geom_boxplot(fill='#A4A4A4', color="black")+
#   theme_classic()
# 
# # cv vs wg
# ggplot(dat_wg, aes(x = WG, y = cv)) +
#   geom_boxplot(fill='#A4A4A4', color="black")+
#   theme_classic()
# 
# # max vs wg
# ggplot(dat_wg, aes(x = WG, y = max)) +
#   geom_boxplot(fill='#A4A4A4', color="black")+
#   theme_classic()
# 
# # avg vs wg
# ggplot(dat_wg, aes(x = WG, y = avg)) +
#   geom_boxplot(fill='#A4A4A4', color="black")+
#   theme_classic()
# 
# # qav vs wg
# ggplot(dat_wg, aes(x = WG, y = qav)) +
#   geom_boxplot(fill='#A4A4A4', color="black")+
#   theme_classic()
# 
# # kur vs wg
# ggplot(dat_wg, aes(x = WG, y = kur)) +
#   geom_boxplot(fill='#A4A4A4', color="black")+
#   theme_classic()
# 
# # ske vs wg
# ggplot(dat_wg, aes(x = WG, y = ske)) +
#   geom_boxplot(fill='#A4A4A4', color="black")+
#   theme_classic()
# 
# # p95 vs age

##############################
###RUN IMPUTATION ALGORITHM###
##############################

# set seed so reproducible
set.seed(5)

# als data only
# build dataframe of reference variables
ref_vars <-
  c('cc', 'p95', 'avg', 'qav', 'cv', 'kur', 'max', 'ske')
x <- dat[, ref_vars]

# ref_vars <-
#   c('cc', 'p95', 'cv')
# x <- dat[, ref_vars]

# build dataframe of target variables
tar_vars <- c('AGE')
y <- dat %>% select(all_of(tar_vars))

# build dataframe of ancillary data
anci <-
  dat[, c('HT', 'CC', 'BA', 'WG', 'AGE', 'p95', 'cc', 'cv', 'avg', 'max')]

# run imputation algorithms
imp <-
  yai(
    x = x,
    y = y,
    method = 'randomForest',
    ntree = 500 * NCOL(y),
    bootstrap = TRUE,
    rfMode = 'regression'
  )

# unsupervised
# imp <-
#   yai(
#     x = x,
#     method = 'randomForest',
#     ntree = 100 * NCOL(x)
#   )

# impute predictions over test data
pred <- impute(imp, ancillaryData = anci,
               method.factor = 'dstWeighted')

# subset values not used in RF
pred <- pred[!(rownames(pred) %in% imp$bootstrap),]
pred <- pred[str_detect(rownames(pred), "\\.") == F,]

# # predictions for multiple k
# pred <- predict(
#   object = imp,
#   newdata = test_data,
#   ancillaryData = anci,
#   observed = T,
#   method = 'median',
#   method.factor = 'mean or median'
# )

# # not sure why getting NAs for observed data
# # but can fill from test_data df
# pred$HT.o <- test_data$HT
# pred$CC.o <- test_data$CC
# pred$BA.o <- test_data$BA
# pred$WG.o <- test_data$WG
# pred$AGE.o <- test_data$AGE
# pred$p95.o <- test_data$p95
# pred$cc.o <- test_data$cc
# pred$cv.o <- test_data$cv
# pred$avg.o <- test_data$avg
# pred$max.o <- test_data$max

# calculate rmse different variables
rmse <- function(obs, est){
  sqrt(mean((obs - est) ^ 2))
}


# calculate rmse metrics
perform_df <-
  data.frame(
    AGE_raw_rmse = rmse(pred$AGE.o, pred$AGE),
    HT_raw_rmse = rmse(pred$HT.o, pred$HT),
    BA_raw_rmse = rmse(pred$BA.o, pred$BA),
    CC_raw_rmse = rmse(pred$CC.o, pred$CC),
    p95_raw_rmse = rmse(pred$p95.o, pred$p95),
    cc_raw_rmse = rmse(pred$cc.o, pred$cc),
    cv_raw_rmse = rmse(pred$cv.o, pred$cv),
    avg_raw_rmse = rmse(pred$avg.o, pred$avg),
    max_raw_rmse = rmse(pred$max.o, pred$max),
    AGE_scaled_rmse = rmse(pred$AGE.o, pred$AGE) /
      sd(pred$AGE.o),
    HT_scaled_rmse = rmse(pred$HT.o, pred$HT) /
      sd(pred$HT.o),
    BA_scaled_rmse = rmse(pred$BA.o, pred$BA) /
      sd(pred$BA.o),
    CC_scaled_rmse = rmse(pred$CC.o, pred$CC) /
      sd(pred$CC.o),
    p95_scaled_rmse = rmse(pred$p95.o, pred$p95) /
      sd(pred$p95.o),
    cc_scaled_rmse = rmse(pred$cc.o, pred$cc) /
      sd(pred$cc.o),
    cv_scaled_rmse = rmse(pred$cv.o, pred$cv) /
      sd(pred$cv.o),
    avg_scaled_rmse = rmse(pred$avg.o, pred$avg) /
      sd(pred$avg.o),
    max_scaled_rmse = rmse(pred$max.o, pred$max) /
      sd(pred$max.o)
  )

# calculate wg accuracy
# create df of WG
wg <- data.frame(obs = pred$WG.o,
                 est = pred$WG)

# make estimate levels match obs levels
levels(wg$est) <-
  c(levels(wg$est), levels(wg$obs)[!(levels(wg$obs) %in% levels(wg$est))])
wg$est <- factor(wg$est, levels = levels(wg$obs))

# create column of match or not
wg$match <- wg$obs == wg$est

# add total percent of matching WG to perform_df
perform_df <- cbind(perform_df,
                    data.frame(wg_accuracy = NROW(wg[wg$match == T,]) /
                                 NROW(wg)))

##########################################
###CREATE WORKING GROUP FREQUENCY TABLE###
##########################################

# create df of WG
wg <- data.frame(obs = pred$WG.o,
                 est = pred$WG)

# make estimate levels match obs levels
levels(wg$est) <- c(levels(wg$est), levels(wg$obs)[!(levels(wg$obs) %in% levels(wg$est))])
wg$est <- factor(wg$est, levels = levels(wg$obs))

# create column of match or not
wg$match <- wg$obs == wg$est

# total percent of matching WG
NROW(wg[wg$match == T,])/NROW(wg)

# check count of different working groups
plyr::count(wg, 'obs')
plyr::count(wg, 'est')

# build accuracy table
accmat_wg <- table("pred" = wg$est, "ref" = wg$obs)

# UA
ua_wg <- diag(accmat_wg) / rowSums(accmat_wg) * 100

# PA
pa_wg <- diag(accmat_wg) / colSums(accmat_wg) * 100

# OA
oa_wg <- sum(diag(accmat_wg)) / sum(accmat_wg) * 100

# build confusion matrix
accmat_wg_ext <- addmargins(accmat_wg)
accmat_wg_ext <- rbind(accmat_wg_ext, "Users" = c(pa_wg, NA))
accmat_wg_ext <- cbind(accmat_wg_ext, "Producers" = c(ua_wg, NA, oa_wg))
#colnames(accmat_wg_ext) <- c(levels(as.factor(shp.train$classes)), "Sum", "PA")
#rownames(accmat_wg_ext) <- c(levels(as.factor(shp.train$classes)), "Sum", "UA")
accmat_wg_ext <- round(accmat_wg_ext, digits = 1)
dimnames(accmat_wg_ext) <- list("Prediction" = colnames(accmat_wg_ext),
                                "Reference" = rownames(accmat_wg_ext))
class(accmat_wg_ext) <- "table"
accmat_wg_ext

# write as csv
# write.csv(perform_avg, wfile = 'D:/ontario_inventory/imputation/ten_fold_validation/rmsd_10_fold_age_wg.csv')
