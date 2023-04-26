library(terra)
library(lidR)
library(sf)
library(stringr)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)
library(ggplot2)
library(readr)

assign_common_name <- function(sp_abbrev) {
  sp_abbrev  <- toupper(sp_abbrev)
  
  dict <- data.frame(SB = "black spruce", 
                     LA = "eastern larch", 
                     BW = "white birch", 
                     BF = "balsam fir", 
                     CE = "cedar", 
                     SW = "white spruce", 
                     PT = "trembling aspen", 
                     PJ = "jack pine", 
                     PO = "poplar", 
                     PB = "balsam poplar", 
                     PR = "red pine", 
                     PW = "white pine", 
                     SX = "spruce", 
                     MR = "red maple", 
                     AB = "black ash", 
                     BY = "yellow birch",
                     OR = 'red oak',
                     CW = 'eastern white cedar',
                     MH = 'hard maple',
                     HE = 'eastern hemlock',
                     BD = 'basswood',
                     CB = 'black cherry',
                     BE = 'american beech',
                     AW = 'white ash',
                     PL = 'largetooth aspen',
                     AG = 'red ash',
                     OW = 'white oak',
                     IW = 'ironwood',
                     OB = 'bur oak',
                     EW = 'white elm',
                     MS = 'silver maple',
                     PS = 'scots pine',
                     OH = 'other hardwoods',
                     BG = 'grey birch',
                     AL = 'alder',
                     SR = 'red spruce',
                     BB = 'blue beech',
                     MT = 'mountain maple',
                     MB = 'black maple',
                     OC = 'other conifers',
                     SN = 'norway spruce',
                     PE = 'silver poplar',
                     HI = 'hickory',
                     AX = 'ash') %>%
    pivot_longer(everything(), names_to = "abb", values_to = "common")
  
  dict$common[match(sp_abbrev, dict$abb)]
  
}

assign_type <- function(sp_common) {
  sp_common  <- tolower(sp_common)
  ifelse(stringr::str_detect(sp_common, pattern = "pine|spruce|fir|cedar|larch|conifers|hemlock"), "Coniferous", "Deciduous")
}

# FSF
VRI <- st_read('D:/ontario_inventory/FSF/FRI/FSF_opi_polygon_CSRS_NAD83_17.shp')

VRI_FOR <- VRI %>%
  st_drop_geometry() %>%
  filter(POLYTYPE == "FOR") %>%
  select(OPI_ID, POLYTYPE, OLEADSPC, OSPCOMP) %>%
  mutate(new_SP = str_match_all(OSPCOMP, "[A-Z]{2}[ ]+[0-9]+")) %>%
  unnest(new_SP) %>%
  mutate(new_SP = as.character(new_SP)) %>%
  separate(new_SP, into = c("SP", "PROP")) %>%
  mutate(PROP = as.numeric(PROP),
         Common = assign_common_name(SP),
         sp_type = assign_type(Common))

# Polygon-level species groups ---- 

poly_dom_type <- VRI_FOR %>%
  group_by(OPI_ID) %>%
  summarize(per_conif = sum(PROP[sp_type == "Coniferous"]), 
            per_decid = sum(PROP[sp_type == "Deciduous"]))

poly_dom_sp <- VRI_FOR %>%
  group_by(OPI_ID) %>%
  slice_max(PROP, n = 1, with_ties = FALSE)


poly_dom_sp_group <- inner_join(poly_dom_type, poly_dom_sp, by = "OPI_ID")

poly_dom_sp_group <- poly_dom_sp_group %>%
  mutate(SpeciesGroup1 = ifelse(PROP >= 70, Common, 
                                ifelse(PROP < 70 & per_conif >= 70, "Mixed Coniferous",
                                       ifelse(PROP < 70 & per_decid >= 70, "Mixed Deciduous", "Mixedwoods"))), 
         SpeciesGroup2 = ifelse(Common == "jack pine" & PROP >= 50 & per_conif >= 70, "Jack Pine Dominated", ifelse(
           Common == "black spruce" & PROP >= 50 & per_conif >= 70, "Black Spruce Dominated", ifelse(
             per_decid >= 70, "Hardwood", ifelse(
               per_decid >= 30 & per_decid <= 70 & per_conif >= 30 & per_conif <= 70, "Mixedwood", "Mixed Conifers"
             )))), 
         SpeciesGroup3 = ifelse(per_conif >= 70, "Softwood", 
                                ifelse(per_decid >= 70, "Hardwood", "Mixedwood"))) %>%
  mutate(across(.cols = starts_with("SpeciesGroup"), .fns = as.factor))

abbreviation_dict <- data.frame(
  SpeciesGroup2 = c("Black Spruce Dominated", 
                    "Jack Pine Dominated", 
                    "Hardwood", 
                    "Mixedwood", 
                    "Mixed Conifers"), 
  SpeciesGroup2_short = c("BS", "JP", "HW", "MW", "MC")
)

poly_dom_sp_group <- inner_join(poly_dom_sp_group, 
                                abbreviation_dict)

write.csv(poly_dom_sp_group,
          'D:/ontario_inventory/imputation/fsf/fsf_SPGROUP.csv',
          row.names = F)
