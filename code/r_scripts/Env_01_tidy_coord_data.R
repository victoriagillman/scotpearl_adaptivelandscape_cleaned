#Victoria Gillman 8th April 2024. 

# This code preprocesses environmental predictor data and links it with species occurrence records.
#     - Reads, cleans, and merges sample IDs with sequencing barcodes.
#     - Adds AnonSiteCodes
#     - Imports, subsets, and downloads environmental predictor data.
#     - Converts occurrence points and visualises them with the modelling region to check same crs.
#     - Converts raster layers, stacks them, and saves the RasterBrick (best method of preserving associated metadata).
#     - Extracts environmental values with occurrence records and saves combined data.
library(dplyr)      # For data manipulation
library(terra)      # For working with spatial data
library(sdmpredictors)  # For accessing environmental predictor datasets
library(raster)     # For working with raster data
library(fuzzySim)   # For working with species distribution modelling data

#### Read in and screen environmental predictors
# At the moment I will just use worldclim factors to keep it simple howver I remember having other features in my SDM and also the new Lancaster paper for river temp? 
# I need to link the sample ID and the barcode IDs. I need to link these with lat/long and I need to link these with extracted worldclim values. First 
sampLOC <- read.csv("data/confidential/pearl_meta_updated_with_WI3andWI4.csv")
barIDs<-read.table("data/PM_Barcode_IDs.txt", header =TRUE)

sampLOC %>%
  group_by(Tissue_Type) %>%
  summarise(count = n())
sampLOC %>%  group_by(Collection_River) %>%  summarise(count = n())
print(sampLOC %>% group_by(VG.River.Names) %>% summarise(count = n()), n = 40)

print(sampLOC %>%  group_by(Latitude, Longitude) %>%  summarise(count = n(), .groups = 'drop'),n=40) #Unique coordinate/sites counts


#Adding ICW763 dublicates 
icx763_row <- filter(sampLOC, Sample.ID == "ICX763")

icx763a_row <- icx763_row %>%
  mutate(Sample.ID = "ICX763a") 
icx763b_row <- icx763_row %>%
  mutate(Sample.ID = "ICX763b")

#Adding the issue of ICW164:171 misnamed as ICW664:671
icw_extended <- bind_rows(
  sampLOC %>%
    filter(Sample.ID %in% paste0("ICW", 164:171)),
  sampLOC %>%
    filter(Sample.ID %in% paste0("ICW", 164:171)) %>%
    mutate(Sample.ID = paste0("ICW", 664:671))
)

# Calling ICW558 ICW558a
icw558a_row <- sampLOC %>%
  filter(Sample.ID == "ICW558") %>%
  mutate(Sample.ID = "ICW558a")


# Binding these all together
sampLOC_extended <- bind_rows(sampLOC, icx763a_row, icx763b_row, icw_extended, icw558a_row)

swabLOC <- full_join(
  dplyr::select(sampLOC_extended, Sample_ID = Sample.ID, Latitude, Longitude, Length.mm., Width.mm., Life_Stage),
  dplyr::select(barIDs, Sample_ID = ID, Barcode, Population),
  by = "Sample_ID"
) %>%
  filter(!is.na(Barcode))

#Add missing moidart locations #already done
# moiLOC <- read.csv("data/confidential/moidart_loc.csv")
# swabLOC <- mutate(swabLOC, Latitude = as.double(Latitude))
# swabLOC <- mutate(swabLOC, Longitude = as.double(Longitude))
# merged_df <- left_join(swabLOC, moiLOC, by = c("Sample_ID"="Sample.ID"))
# swabLOC <- mutate(merged_df, 
#                     Latitude = ifelse(is.na(Latitude.x), Latitude.y, Latitude.x),
#                     Longitude = ifelse(is.na(Longitude.x), Longitude.y, Longitude.x)) %>%
#   dplyr::select(-Latitude.x, -Latitude.y, -Longitude.x, -Longitude.y)


# Manually putting in wrong sample id latitude and longitude for ICW558b
# Manually input 
swabLOC$Latitude[swabLOC$Sample_ID == 'ICW558b'] <- 57.73493
swabLOC$Longitude[swabLOC$Sample_ID == 'ICW558b'] <- -4.890829


#Get the missing lat/long sample ids

sample_ids_missing_latitude <- swabLOC %>%
  filter(!is.na(Barcode) & is.na(Latitude)) %>%
  dplyr::select(Sample_ID, Population)

nrow(sample_ids_missing_latitude)
dim(swabLOC)
remove(icx763a_row, icx763b_row, icw_extended, icw558a_row, icx763_row, sampLOC, barIDs, merged_df,moiLOC)

# Save the cleaned pearl dataset
swabLOC_clean<-swabLOC %>% filter(!is.na(Latitude)) 
write.csv(swabLOC_clean, "./data/confidential/swab_locations_cleaned.csv",row.names=F)
swabLOC_clean <- read.csv("./data/confidential/swab_locations_cleaned.csv")


library(dplyr)
library(sf)

## Importing site locations and getting their vice county information
# Import csv file containing coordinates
coords = read.csv("data/confidential/swab_locations_cleaned.csv")

length(unique(coords$Population))
# coords<-coords%>% group_by(Population) %>%
#   summarise(Latitude = first(Latitude), Longitude = first(Longitude)) # I dont think I need to summarise this

name_mapping <- c( #some of the names in my coords file dont match those in the genomics data
  "Bhiaghlann" = "Briaghlann",
  "Conon_Gharbhrain" = "Conon_Glascarnoch",
  "Sruthan_Braigh" = "Struthan_Bhraigh"
)
coords$Population <- recode(coords$Population, !!!name_mapping)

# Add a vice-county column
coords_sf <- st_as_sf(coords, coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

polys <- st_read("data/scotland_vice_counties/scotland_vicecounties.shp") %>%
  st_transform(st_crs(coords_sf))

coords_sf <- st_join(coords_sf, polys)
coords_vc <- coords_sf %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  bind_cols(coords %>% dplyr::select(Latitude, Longitude))
##
coords_anon_names <- coords_vc %>%
  mutate(
    VCInitials = gsub("(?<=[A-Z])[^A-Z]+", "", VCNAME, perl = TRUE)  # Extract initials from VCNAME correctly
  ) %>%
  group_by(VCNAME) %>%
  mutate(
    AnonSiteCode = paste0(VCInitials, dense_rank(Population))  # Combine initials with dense rank of SiteCode
  ) %>%
  ungroup() %>%
  dplyr::select(-VCInitials)  # Remove the intermediate VCInitials column if not needed
length(unique(coords_anon_names$Population)) #23
length(unique(coords_anon_names$AnonSiteCode)) #23
unique(coords_anon_names$Population)
unique(coords_anon_names$AnonSiteCode)


write.csv(coords_anon_names, "data/confidential/swab_coord_vc_anon_names.csv", row.names = F)

swabLOC_clean<-read.csv("data/confidential/swab_coord_vc_anon_names.csv") ### New anonomised version
plot(polys$geometry)
plot(coords_sf$geometry, add=T, col="deeppink", cex=5,pch=17)
# library(mapview)
# mapview(coords_sf)




#### Importing WorldClim terrestrial environmental data ####
# First we want Scotland to subset climate raster maps so they are easier to handle
# scotland<-vect("./data/scotland_outline/govgeoportalsubset_scot.shp")
# plot(scotland)
# crs(scotland)
my_window<-ext(polys)

#### Commented out section is Selecting and Downloading appropriate predictor variables ####

# install.packages("sdmpredictors")
library(sdmpredictors)

# Access available datasets and layers
vignette("quickstart", package = "sdmpredictors")
pred_datasets <- list_datasets(terrestrial = TRUE, marine = FALSE)
pred_layers <- list_layers(datasets = pred_datasets)

#
# # Save dataset and layer information to a CSV file
if (!file.exists("./data/too_large_files")) {
  dir.create("./data/too_large_files")
}
if (!file.exists("./data/predictor_download")) {
  dir.create("./data/predictor_download")
}
write.csv(pred_layers, "./data/predictor_download/worldclim_predictors_metadata.csv")

# Choose specific dataset and layers
layers_choices <- unique(pred_layers[pred_layers$dataset_code == "WorldClim", c("name", "layer_code")])
layers_choices$row_number <- seq_len(nrow(layers_choices))
layers_choice <- layers_choices[c(1:20),]


write.csv(layers_choice, "./data/predictor_download/worldclim_layerschoice_metadata.csv")

# # Define folder for downloading/fetching the variables' map layers
if (!file.exists("./data/too_large_files/worldclim_downloads_fullrasters/")) {
  dir.create("./data/too_large_files/worldclim_downloads_fullrasters")
}
options(sdmpredictors_datadir = "./data/too_large_files/worldclim_downloads_fullrasters")

# # Download the layers
library(raster)

layers <- load_layers(layers_choice$layer_code, rasterstack = FALSE) #downloads the files, now just need to import them below
# # Convert to terra as faster
 layers <- lapply(layers, rast)

#### Importing the World-size environmental predictors downloaded above
layers <- rast(list.files("./data/too_large_files/worldclim_downloads_fullrasters", pattern = "\\.tif$", full.names = TRUE))

# Crop layers to the extent of Scotland
layers <- lapply(layers, crop, ext(polys))
length(layers)

plot(layers[[1]], main = names(layers)[1]) 
# Check resolution is all the same
unique(pred_layers[pred_layers$dataset_code == "WorldClim", ]$cellsize_lonlat) 


# Making occurrence points a Spatvector
swabLOC_clean$Longitude <- as.numeric(swabLOC_clean$Longitude)
swabLOC_clean$Latitude <- as.numeric(swabLOC_clean$Latitude)
occurrence_points <- vect(swabLOC_clean, geom = c("Longitude", "Latitude"), crs = "epsg:4326") #WGS84
plot(occurrence_points$geometry)
dim(occurrence_points)


#Checking everything overlaps correctly
# mod_region<-scotland
plot(layers[[1]])
plot(my_window, border = "red", add = TRUE)  # it appears only if it's within the plotting region
plot(occurrence_points, col = "darkgreen", add = TRUE)
# plot(mod_region, border = "green", lwd = 2, add = TRUE)


# plot(scotland)
# plot(occurrence_points, bg="#AB1A24", col="#AB1A24", add=T, pch=25)

# Check overlap in more detail (this is where Freshwater points do not overlap- hence using terrestrial)
plot(layers[[1]], xlim = c(-6.5, -5.5), ylim = c(56.5, 57)) 
plot(occurrence_points, col = "black", cex= 5, add=TRUE)
# plot(mod_region, border = "green", lwd = 2, add = TRUE) ###continued issue of points not overlapping with raster


# Okay great everything overlaps so lets export the Raster environment files
# Convert each SpatRaster to Raster* object
raster_layers <- lapply(layers, raster)
# Stack the Raster* objects into a RasterStack
stack_layers <- stack(raster_layers)
names(stack_layers)
plot(stack_layers[[1]])
# Make a new directory "./data/too_large_files/worldclim_downloads_fullrasters
if (!file.exists("./data/too_large_files/worldclim_downloads_fullrasters/worldclim_layers_cut")) dir.create("./data/too_large_files/worldclim_downloads_fullrasters/worldclim_layers_cut")
# Writing as a .grd file preserves the layer names, if have to import as a single-file layered .tif I have the layer names below. 
writeRaster(stack_layers,overwrite=TRUE, filename = "./data/too_large_files/worldclim_downloads_fullrasters/worldclim_layers_cut/worldclim_layers_cut.grd") 

# If must save as a single-file layered .tif the names and metadata goes so must rename before we get confused
# layers_cut <- rast("./data/predictor_download/worldclim_downloads_fullrasters/worldclim_layers_cut/worldclim_layers_cut.tif", full.names = TRUE)
# names(layers_cut)
# original_names <- c(
#   "WC_alt_lonlat", "WC_bio1_lonlat", "WC_bio10_lonlat", "WC_bio11_lonlat",
#   "WC_bio12_lonlat", "WC_bio13_lonlat", "WC_bio14_lonlat", "WC_bio15_lonlat",
#   "WC_bio16_lonlat", "WC_bio17_lonlat", "WC_bio18_lonlat", "WC_bio19_lonlat",
#   "WC_bio2_lonlat", "WC_bio3_lonlat", "WC_bio4_lonlat", "WC_bio5_lonlat",
#   "WC_bio6_lonlat", "WC_bio7_lonlat", "WC_bio8_lonlat", "WC_bio9_lonlat"
# )
# names(layers_cut) <- original_names
# plot(layers_cut)

# layers_cut<-brick("./data/predictor_download/worldclim_downloads_fullrasters/worldclim_layers_cut/worldclim_layers_cut.grd")
# names(layers_cut)
# plot(layers[[1]], xlim = c(-6.5, -5.5), ylim = c(56.5, 57)) 
# plot(occurrence_points, col = "black", cex= 5, add=TRUE)
# # plot(mod_region, border = "green", lwd = 2, add = TRUE) ###continued issue of points not overlapping with raster
# # summary(layers_cut) #a RasterBrick is more efficient way of storing rasterstack with same extent and resolution
# object.size(layers_cut) # 16160 bytes
# object.size(layers) # 26288 bytes

plot(stack_layers[[1]], xlim = c(-6.5, -5.5), ylim = c(56.5, 57)) 
plot(occurrence_points, col = "black", cex= 5, add=TRUE)
# plot(mod_region, border = "green", lwd = 2, add = TRUE) ###continued issue of points not overlapping with raster

head(occurrence_points)

##### DONT THINK I NEED TO DO BELOW GRIDDING THIS AS THIS ISNT AN SDM #####
# # install.packages("fuzzySim") DONT THINK I NEED TO DO THIS AS THIS ISNT AN SDM
# # library(fuzzySim)
# ?gridRecords
# gridded_data <- gridRecords(rst = layers, pres.coords = occurrence_points[ , c("Longitude", "Latitude")])
# head(gridded_data)
# 
# nrow(gridded_data)  # should be the same number as:
# sum(!is.na(values(layers_cut[[1]])))
# 
# names(gridded_data)
# myspecies
#####


#### Convert sf to terra file and extract the env parameters for the swab locations!! 
spat_raster <- rast(stack_layers) #rasterbrick is sf but extract is terra
extr_values<-terra::extract(spat_raster, occurrence_points) #rasterbrick is raster package! 
combined_data_sv <- cbind(occurrence_points, extr_values)
combined_data <- cbind(swabLOC_clean, extr_values)
dim(combined_data_sv)

head(combined_data)
write.csv(combined_data, "data/confidential/swab_enviroment.csv", row.names = FALSE)


### Done! ###
# save.image("code/R_code/backup_workspace_images/extract_env_for_swab_loc.RData") # As a backup :) 
# 
# #### Ggplots
# library(ggplot2)
# library(sf)
# env<-read.csv("data/confidential/swab_enviroment.csv")
# # env<-na.omit(env)
# env <- env[!is.na(env$Latitude),]
# 
# env_sf <- st_as_sf(env, coords = c("Longitude", "Latitude"), crs = st_crs(scotland))
# mapview(env_sf)
# sum(is.na(env)) #163 #17/4 is 67
# 
# ggplot()+
#   geom_sf(scotland)+
#   theme_void()
