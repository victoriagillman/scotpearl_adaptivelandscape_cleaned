#### Importing RASTERS Wilderness, Soil ph, TOC, Phosphorous and then adding on fw predictors, host use, catchment area

library(dplyr)      # For data manipulation
library(terra)      # For working with spatial data
library(raster)     # For working with raster data
library(sf) # For shp files

# load("code/R_code/backup_workspace_images/compiling_full_predictor_dataset_for_trimming_next.RData")

gc() #garbage cleanup
rm(list = ls())  # Removes all objects from the environment

# Point data
coords_df <- read.csv("./data/confidential/all_freshwater_points2nearestcell_values.csv") %>%
  mutate(
    Longitude = as.numeric(Longitude),
    Latitude = as.numeric(Latitude)
  )
coords_terra <- vect(coords_df, geom = c("Longitude", "Latitude"), crs = "epsg:4326") #WGS84
coords_sf <- st_as_sf(coords_terra) %>% dplyr::select(Sample_ID) # because some are extracted using sf
plot(coords_terra)

#######################################################################################################
# Worldclim/Freshwater sdmpredictors- I preferred to extract these in the sdmpredictors code and attaching later

# terr_layers_cut<-brick("./output/predictor_download/worldclim_downloads_fullrasters/worldclim_layers_cut/worldclim_layers_cut.grd")
# fw_layers_cut<-brick("./output/predictor_download/freshwater_sdmpredictors/layers_cut/layers_cut.grd")
# crs<-crs(terr_layers_cut) #Deprecated Proj.4 representation: +proj=longlat +datum=WGS84 +no_defs 
# 
# ## fw extract 
# fw_orig_values <- terra::extract(rast(fw_layers_cut), coords_terra, na.rm=TRUE, method="bilinear")

#######################################################################################################
##  wildness
wild_rm <- raster("./data/too_large_files/other_predictors/WILDNESS-REM_SCOTLAND_TIFF_27700/WILDNESS-REM_SCOTLAND.tif") #rast() does not work
wild_rm_reproj <- projectRaster(wild_rm, crs = CRS("EPSG:4326")) #12  min, #rast() did not work
plot(wild_rm_reproj, xlim = c(-6.5, -5.5), ylim = c(56.5, 57)) #Check everything lines up
plot(coords_terra, add=T)
plot(fw_layers_cut[[1]], add=T)

coords_terra_w <- terra::extract(rast(wild_rm_reproj), coords_terra, na.rm=TRUE)
wild_points <- cbind(coords_df, wildness=coords_terra_w$OBJECTID)
coords_terra_w <- cbind(coords_terra, coords_terra_w)

# #######################################################################################################
# ## AWC- water capacity
# # AWC <- vect("./output/predictor_download/other_predictors/Hutton_AWC_OSGB_2019/AWC_2019_1.shp") #not the best way
# AWC <- st_read("./output/predictor_download/other_predictors/Hutton_AWC_OSGB_2019") # if I piped it it takes forever so I broke it down
# AWC_transformed <- st_transform(AWC, crs = crs)
# AWC_cleaned <- st_make_valid(AWC_transformed)
# awc_points_full <- st_intersection(coords_sf, AWC_cleaned) 
# awc_points<-awc_points_full %>% dplyr::select(Sample_ID, AWC) %>% st_drop_geometry()
# 
# #######################################################################################################
# ## soil pH
# ph <- st_read("./output/predictor_download/other_predictors/Hutton_pH_OpenData") %>%
#   st_transform(crs = crs) %>%
#   st_make_valid()
# ph_points_full<- st_intersection(coords_sf, ph)
# ph_points<-ph_points_full %>% dplyr::select(Sample_ID, pH_W_Med) %>% st_drop_geometry()
# 
# #######################################################################################################
# ## TOC
# TOC <- st_read("./output/predictor_download/other_predictors/Hutton_TOC_WGS84/") %>%
#   st_transform(crs = crs) %>%
#   st_make_valid()
# TOC_points_full<-st_intersection(coords_sf, TOC) 
# TOC_points<-TOC_points_full %>%  dplyr::select(Sample_ID, organicCar) %>% st_drop_geometry()
# 
# #######################################################################################################
# ## Phos
# Phos <- st_read("./output/predictor_download/other_predictors/Soil_Phosphorus_sorption_capacity/") %>%
#   st_transform(crs = crs) %>%
#   st_make_valid()
# plot(Phos, xlim = c(-6.5, -5.5), ylim = c(56.5, 57))
# Phos_points_full<-st_intersection(coords_sf, Phos)  
# Phos_points<-Phos_points_full %>%  dplyr::select(Sample_ID, NEW_P_INDE) %>% st_drop_geometry()
#######################################################################################################
# Host use from VPritchard 
populations <- coords_df %>%
  distinct(Population) %>%
  pull(Population)
populations #Pull pop names in dataset

# Add the Host column using mutate() and case_when()
host <- coords_df %>%
  mutate(Host = case_when(
    AnonSiteCode == "EI1" ~ "Salmon",
    AnonSiteCode == "OH1" ~ "Unknown",
    AnonSiteCode == "WS4" ~ "Trout",
    AnonSiteCode == "ES1" ~ "Unknown",
    AnonSiteCode == "WR4" ~ "Trout",
    AnonSiteCode == "WS5" ~ "Unknown",
    AnonSiteCode == "WS2" ~ "Only trout available",
    AnonSiteCode == "WS1" ~ "Only trout available",
    AnonSiteCode == "WI4" ~ "Only trout available",
    AnonSiteCode == "EI5" ~ "Salmon",
    AnonSiteCode == "WR3" ~ "Trout",
    AnonSiteCode == "WR2" ~ "Trout",
    AnonSiteCode == "WS3" ~ "Only trout available",
    AnonSiteCode == "WI2" ~ "Only trout available",
    AnonSiteCode == "WI1" ~ "Only trout available",
    AnonSiteCode == "WR1" ~ "Salmon",
    AnonSiteCode == "NE1" ~ "Unknown",
    AnonSiteCode == "WI3" ~ "Salmon",
    AnonSiteCode == "ER2" ~ "Unknown",
    AnonSiteCode == "ER1" ~ "Unknown",
    TRUE ~ NA_character_ # Default case to handle any other values
  ))

#######################################################################################################
## Gather all data above
compiled_data <-
  # awc_points  %>%
  wild_points %>% 
  # left_join(wild_points %>% dplyr::select(Sample_ID, wildness), by = "Sample_ID") %>%
  # left_join(ph_points, by = "Sample_ID") %>%
  # left_join(TOC_points, by = "Sample_ID") %>%
  # left_join(Phos_points, by = "Sample_ID") %>%
  left_join(host %>% dplyr::select(Sample_ID, Host), by = "Sample_ID")

# save.image("code/R_code/backup_workspace_images/heavy_predictor_mapping_hopefullylasttime.RData") # As a backup :) 
# load("code/R_code/backup_workspace_images/heavy_predictor_mapping_hopefullylasttime.RData")

# write.csv(compiled_data, "data/confidential/other_env_data.csv", row.names=F)
# compiled_data<-read.csv("data/confidential/other_env_data.csv")
# 
# #######################################################################################################
# # Adding in all freshwater sdmpredictoris data
# fw_points_moved_df <- read.csv("data/confidential/all_freshwater_points2nearestcell_values.csv")
# 
# mega_env_df <- inner_join(fw_points_moved_df, compiled_data, by = "Sample_ID")
# # write.csv(mega_env_df, "data/confidential/all_fw_and_other_env_data.csv", row.names=F)
# 
# # This shouldnt be needed on run again- as I added the anon names right at the start. However this is a way so I dont have to do all that data extraction again!
# 
# #######################################################################################################
# # Adding anon locations 
# swabLOC_anon<-read.csv("data/confidential/swab_coord_vc_anon_names.csv")
# final_meta <- left_join(swabLOC_anon, mega_env_df, by = "Sample_ID")
# # Identify common columns
# common_columns <- intersect(names(swabLOC_anon), names(mega_env_df))
# common_columns <- setdiff(common_columns, "Sample_ID")  # Exclude the join column
# 
# # Remove common columns from new_meta
# new_meta_filtered <- mega_env_df %>% dplyr::select(-one_of(common_columns))
# 
# # Perform the join on the Sample_ID column
# final_meta <- left_join(swabLOC_anon, new_meta_filtered, by = "Sample_ID") %>%
#   dplyr::select(-X, -X.1, ID)
# 
# #######################################################################################################
# ## Adding catchment area
# catch_col<-read.csv("data/confidential/catchment_area_for_each_sample.csv")
# final_meta_final <- final_meta %>%
#   left_join(catch_col %>% dplyr::select(Sample_ID, catchment_area_km), by = c("Sample_ID"))%>%
#   dplyr::select(-Host, Host)
# write.csv(final_meta_final, "data/confidential/anon_all_fw_and_other_env_data.csv",row.names=F)

## Other altitude can be added in here but that is for me post-thesis probably (obtained from Env_altitude_from_os_terrain50)
# os_alt<-read.csv( "data/os_altitude50_to_sample_id.csv")
# final_meta_final <- final_meta %>%
#   left_join(os_alt, by = c("Sample_ID"))%>%
#   dplyr::select(-Host, Host)
# write.csv(final_meta_final, "data/confidential/anon_all_fw_and_other_env_data.csv",row.names=F)
write.csv(compiled_data, "data/confidential/anon_all_fw_and_other_env_data.csv",row.names=F)

# save.image("code/R_code/backup_workspace_images/compiling_full_predictor_dataset_for_trimming_next.RData") # As a backup :) 
