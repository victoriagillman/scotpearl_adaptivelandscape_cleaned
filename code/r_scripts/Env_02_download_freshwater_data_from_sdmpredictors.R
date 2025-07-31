###This is me trying to use freshwater variables
# install.packages("sdmpredictors")
library(sdmpredictors)

# load("code/R_code/backup_workspace_images/extract_all_fw_from_swab_loc.RData") # As a backup :) 

#### Installing freshwater raster datasets from sdmpredictors ####===========================================
# Access available datasets and layers
vignette("quickstart", package = "sdmpredictors")
pred_datasets <- list_datasets(terrestrial = TRUE, marine = FALSE)
pred_layers <- list_layers(datasets = pred_datasets)

datadirectory<-"./data/too_large_files/all_freshwater_sdmpredictors/"

# Save dataset and layer information to a CSV file
if (!file.exists(datadirectory)) {
  dir.create(datadirectory)
}
write.csv(pred_layers, file= paste0(datadirectory,"metadata.csv"))

# Choose specific dataset and layers
layers_choices <- unique(pred_layers[pred_layers$dataset_code == "Freshwater", c("name", "layer_code")])
layers_choices$row_number <- seq_len(nrow(layers_choices))
# layers_choice <- layers_choices[c(4, 8, 10, 83:101,157:168,303:312),]


write.csv(layers_choices, file=paste0(datadirectory,"layerschoice_metadata.csv"))

# Define folder for downloading/fetching the variables' map layers
if (!file.exists(paste0(datadirectory,"downloads_fullrasters"))) {
  dir.create(paste0(datadirectory,"downloads_fullrasters"))
}
options(sdmpredictors_datadir = paste0(datadirectory,"downloads_fullrasters"))

# Download or fetch the layers
library(raster)
library(terra)
remove(l)
layers <- load_layers(layers_choices$layer_code, rasterstack = FALSE) #downloads the files, now just need to import them below
# Convert to terra as faster
layers <- lapply(layers, rast)
# datadirectory<-"C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/predictor_download/freshwater_sdmpredictors/" #quick import for the histogram plot :( C:\Users\r01vg21\OneDrive - University of Aberdeen\2.PROJECT\scotpearl_popgen\output\predictor_download\freshwater_sdmpredictors\downloads_fullrasters

layers <- rast(list.files(paste0(datadirectory,"downloads_fullrasters/"), pattern = "\\.tif$", full.names = TRUE)) #imports the global tiffs we just downloaded
names(layers)
layer_names <- names(layers)


# sampLOCo <- read.csv("../../data/confidential/swab_locations_cleaned.csv") #older version
# sampLOC_all_col<-read.csv("data/confidential/swab_coord_vc_anon_names.csv") # with new anonymised dataset
sampLOC<-read.csv("data/confidential/catchment_area_for_each_sample.csv")
sampLOC <- na.omit(sampLOC[, c("Longitude", "Latitude")])
sampLOC$Longitude <- as.numeric(sampLOC$Longitude)
sampLOC$Latitude <- as.numeric(sampLOC$Latitude)
occurrence_points <- vect(sampLOC, geom = c("Longitude", "Latitude"), crs = "epsg:4326") #WGS84
ext(occurrence_points)
plot(occurrence_points$geometry)
plot(layers[[1]])

layers <- lapply(layers, crop, ext(occurrence_points)+1) #somehow removes layer names
names(layers) <- layer_names
# lengthlayers_choicelength(layers)

plot(layers[[1]], main = names(layers)[1]) 
# plot(layers)
# Check resolution is all the same
unique(pred_layers[pred_layers$dataset_code == "Freshwater", ]$cellsize_lonlat) 

#Checking everything overlaps correctly
#  VISUALISING OCCURANCE POINTS AND MODDELING REGION
mod_region<-ext(occurrence_points)+1
plot(mod_region)
plot(occurrence_points, col = "darkgreen", add = TRUE)
plot(mod_region, border = "green", lwd = 2, add = TRUE)
plot(layers[[1]], add=TRUE)
title(paste0(names(layers[[1]])))


plot(layers[[42]], xlim = c(-5.5, -6.5), ylim = c(56.5, 57), main=paste0(names(layers[[42]]))) 
plot(occurrence_points, col = "black", cex= 5, add=TRUE)
plot(mod_region, border = "green", lwd = 2, add = TRUE) ###continued issue of points not overlapping with raster

# Okay great everything overlaps so lets export the Raster environment files
# Convert each SpatRaster to Raster* object
raster_layers <- lapply(layers, raster)
# Stack the Raster* objects into a RasterStack
stack_layers <- stack(raster_layers)
names(stack_layers)
# Make a new directory
# if (!file.exists(paste0(datadirectory,"layers_cut"))) dir.create(paste0(datadirectory,"layers_cut"))
# # Writing as a .grd file preserves the layer names, if have to import as a single-file layered .tif I have the layer names below. 
# writeRaster(stack_layers,overwrite=TRUE, filename = paste0(datadirectory,"layers_cut/layers_cut.grd")) #save
# layers_cut<-brick(paste0(datadirectory,"layers_cut/layers_cut.grd")) #import
# names(layers_cut)
# 
# 
# ## Count rivers which do not overlap ################
# # Extract values from the raster layer at occurrence points
# values_at_points <- extract(layers[[1]], occurrence_points, bind=T) %>% as.data.frame(.) %>%
#   cbind(sampLOC_all_col, .) %>% dplyr::select(Latitude,Longitude, AnonSiteCode, names(layers[[1]]) )
# 
# # Filter occurrence points that overlap (non-NA values) and count by longitude
# values_at_points %>%
#   group_by(AnonSiteCode) %>% summarise(n = n()) %>% print(n=)
# values_at_points %>%
#   filter(is.na(.[[names(layers[[1]])]])) %>%
#   pull(AnonSiteCode) %>%
#   unique() %>%
#   print()
# read.csv("data/confidential/meta_for_filtered_genomes.csv") %>%    pull(AnonSiteCode) %>%
#   unique() 
# 
# # Read in the filtered dataset
# filtered_sites <- read.csv("data/confidential/meta_for_filtered_genomes.csv") %>%
#   pull(AnonSiteCode) %>%
#   unique()
# 
# # Get the sites with NA in the raster extraction
# na_sites <- values_at_points %>%
#   filter(is.na(.[[names(layers[[1]])]])) %>%
#   pull(AnonSiteCode) %>%
#   unique()
# 
# # Find sites from the filtered dataset that did not have a successful extraction (NA values)
# missing_sites <- setdiff(filtered_sites, na_sites)
# 
# # Print missing sites
# print(missing_sites)
# 
# # Count the number of missing sites
# length(missing_sites)



######### EXTRACTING VARIABLES WITHIN THIS CODE WORKS BETTER THAN TRYING TO SAVE AND IMPORT THE LAYERS AGAIN ########
# Read and prepare sample location data
library(rSDM)
sampLOC <- read.csv("./data/confidential/catchment_area_for_each_sample.csv")
sampLOC$Longitude <- as.numeric(sampLOC$Longitude)
sampLOC$Latitude <- as.numeric(sampLOC$Latitude)
crs<-crs(stack_layers[[1]])
plot(stack_layers[[1]])
# Convert sample locations to sf object
occurrence_points_rSDM <- locs2sf(sampLOC, lon.col = "Longitude", lat.col = "Latitude", crs = crs)
plot(stack_layers[[1]])
points(occurrence_points_rSDM)
# Move coordinates to nearest cell
move_coords <- points2nearestcell(
  locs = occurrence_points_rSDM,
  ras = rast(stack_layers[[1]]), #I added a subset here- not sure if this will be good or not? i think I rembmer the program naturally just uses the first layer
  move = TRUE,
  distance = 10050,
  table = TRUE,
  map = "base"
)

# Convert moved coordinates to vect object
vect_movecoords <- vect(move_coords)

rast<-rast(stack_layers)
fw_points_moved <- terra::extract(rast, vect_movecoords, bind = FALSE)

# Combine with original sample location data
fw_points_moved_df <- cbind(sampLOC[!is.na(sampLOC$Latitude), ], fw_points_moved)

# Write to CSV
# write.csv(fw_points_moved_df, "data/confidential/all_freshwater_points2nearestcell_values.csv", row.names = FALSE)

# save.image("code/R_code/backup_workspace_images/extract_all_fw_from_swab_loc.RData") # As a backup :) 



# Histogram of moved points

# Extract original coordinates from occurrence_points_rSDM
orig_coords <- st_coordinates(occurrence_points_rSDM)

# Extract new coordinates from move_coords
new_coords <- st_coordinates(move_coords)

# Calculate the distances (in metres) between original and new points
# Assumes CRS is projected (e.g., UTM) for distances in metres, or it uses st_distance otherwise
distances <- st_distance(occurrence_points_rSDM, move_coords, by_element = TRUE)

# Convert distances to numeric vector (in case it's a matrix)
distances <- as.numeric(distances)

# Create a data frame to store distances
distance_df <- data.frame(distance = distances)

# Plot the histogram of distances moved
library(ggplot2)
hist_pl<-ggplot(distance_df, aes(x = distance)) +
  geom_histogram(binwidth = 500, fill = "skyblue", color = "black") +
  labs(
       x = "Distance Moved (metres)", 
       y = "Number of Points") +
  theme_classic()
  # theme_minimal()
hist_pl <- ggplot(distance_df, aes(x = distance)) +
  geom_histogram(binwidth = 500, fill = "skyblue", color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Distance Moved (metres)", 
    y = "Number of Points") +
  theme_classic()
  # theme(
  #   panel.grid.major = element_line(colour = "grey80", linewidth = 0.4),
  #   panel.grid.minor = element_line(colour = "grey90", linewidth = 0.2)
  # )

hist_pl
ggsave(plot= hist_pl,
       width = 5,
       height = 3,
       units = c( "in"),
       dpi = 600,
       filename = "./output/figures/historgram_dispances_bumped.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/historgram_dispances_bumped.png")

###################





### to use the freshwater variables using a bumping 
rastlayers <- lapply(list.files("./output/predictor_download/downloads_fullrasters/", pattern = "\\.tif$", full.names = TRUE), raster) #imports as rasters
# Convert the list of raster files to a RasterStack
stacklayers <- stack(lapply(list.files("./output/predictor_download/downloads_fullrasters/", pattern = "\\.tif$", full.names = TRUE), raster))
nlayers(stacklayers)
res(stacklayers)
crs(stacklayers)
names(stacklayers)
points <- SpatialPoints(coords = sampLOC[, c("Longitude", "Latitude")], proj4string = CRS("+proj=longlat +datum=WGS84"))


myExpl <- points %>% 
  raster::extract(stacklayers, ., df = TRUE) %>% 
  dplyr::select(-ID)
dim(myExpl) # 1128
myExpl
# Loop through each layer in myExpl
for (layer_name in names(myExpl)) {
  # Check for NA values in the current layer
  na_index <- is.na(myExpl[[layer_name]])
  
  # Replace NA values with the mean value extracted from the corresponding layer in predictors
  myExpl[na_index, layer_name] <- raster::extract(stacklayers[[layer_name]], points[na_index,], buffer = 100000, fun = mean)
}

dim(myExpl)
myExpl



# Function to replace NA values in myExpl with extracted raster values
fill_na_with_extraction <- function(myExpl, stacklayers, points, buffer_size = 100) {
  for (layer_name in names(myExpl)) {
    na_index <- is.na(myExpl[[layer_name]])
    if (sum(na_index) > 0) {
      extracted_values <- raster::extract(stacklayers[[layer_name]], points[na_index,], buffer = buffer_size)
      myExpl[na_index, layer_name] <- extracted_values
    }
  }
  return(myExpl)
}

# Fill missing predictor values in myExpl using extraction with a buffer
myExpl_filled <- fill_na_with_extraction(myExpl, stacklayers, points, buffer_size = 100)

# Check if NA values have been filled
summary(myExpl_filled)




#====================================
# # Enlarge the extent of the raster brick ####
# buffer_extent <- extent(fw_layers_cut)
# buffer_extent <- buffer_extent + 10000  # Increase the extent by 10 units (adjust as needed)
# # Define a matrix for the focal function (e.g., a 3x3 window)
# focal_matrix <- matrix(1, nrow=3, ncol=3)
# 
# # Compute the focal mean for each layer in the brick
# focal_mean_brick <- stack(lapply(1:nlayers(fw_layers_cut), function(i) {
#   focal(fw_layers_cut[[i]], w=focal_matrix, fun=mean, na.rm=TRUE)
# }))
# # Function to fill NoData areas with focal mean values
# fill_nodata <- function(layer, focal_layer) {
#   overlay(layer, focal_layer, fun=function(x, y) {
#     ifelse(is.na(x), y, x)
#   })
# }
# 
# # Apply this function to each layer in the brick
# filled_brick <- stack(mapply(fill_nodata, as.list(fw_layers_cut), as.list(focal_mean_brick)))
# plot(filled_brick[[1]],xlim = c(-6.5, -5.5), ylim = c(56.5, 57))
# plot(occurrence_points, add=T)
# plot(fw_layers_cut[[1]], add=F)
# 
# cellStats(fw_layers_cut[[1]], range)
# cellStats(filled_brick[[1]], range)
#############
extracted_values <- extract(rast(fw_layers_cut), occurrence_points, method='simple')
na_indices <- which(is.na(extracted_values))
point_coords <- geom(occurrence_points)[1, c("x", "y")]
non_na_cells <- which(!is.na(values(fw_layers_cut)))

## Nearest neighbour approach
# Function to find the nearest non-NA cell value
extract_nearest <- function(points, raster_layer) {
  # Get the initial extraction results
  extracted_values <- extract(rast(raster_layer), points, method='simple')
  
  # Find points with NA values
  na_indices <- which(is.na(extracted_values))
  
  if(length(na_indices) > 0) {
    for(i in na_indices) {
      # Ensure that i does not exceed the number of rows in points
      if (i <= nrow(points)) {
        # Get the coordinates of the point
        point_coords <- geom(points)[i, c("x", "y")]
        
        # Find the nearest cell with a non-NA value
        non_na_cells <- which(!is.na(values(raster_layer)))
        cell_coords <- xyFromCell(raster_layer, non_na_cells)
        
        # Compute distances to all non-NA cells
        distances <- spDistsN1(cell_coords, point_coords, longlat = FALSE)
        
        # Get the value of the nearest cell
        nearest_cell <- non_na_cells[which.min(distances)]
        extracted_values[i] <- raster_layer[nearest_cell]
      }
    }
  }
  
  return(extracted_values)
}

# Apply this function to each layer in the brick
extracted_values_list <- lapply(1:nlayers(fw_layers_cut), function(i) {
  extract_nearest(occurrence_points, fw_layers_cut[[i]])
})


# Combine the extracted values into a data frame
extracted_values_df <- as.data.frame(do.call(cbind, extracted_values_list))
names(extracted_values_df) <- paste0("layer_", 1:nlayers(fw_layers_cut))

# Add the extracted values to the points data
occurrence_points@data <- cbind(occurrence_points@data, extracted_values_df)

