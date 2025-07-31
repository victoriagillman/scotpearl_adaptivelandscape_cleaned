### Lets get the catchment info for the pearl points
## 1- By assigning each a catchment they are in

# setwd("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen")

library(dplyr)
library(sf)

scotland_polygon <- st_read("data/scotland_vice_counties/scotland_vicecounties.shp") %>% rmapshaper::ms_simplify(keep = 0.05) %>% st_union() %>% st_cast("POLYGON")

# coords <- read.csv("data/confidential/meta_for_filtered_genomes.csv") #not cleaned names
coords <- read.csv("./data/confidential/swab_enviroment.csv")

coords_sf<-coords %>% st_as_sf( coords = c("Longitude", "Latitude"), crs = 4326) %>% st_transform(st_crs(scotland_polygon))

polys <- st_read("data/too_large_files/catchment/MDRExternalStorage.Recordset_12840.shp") %>%
  st_set_crs(st_crs(scotland_polygon))

# Attach catchment area column and save
polys$catchment_area_km <- st_area(polys)%>% 
  units::set_units("km^2") %>% 
  as.numeric()  # Convert to numeric, removing units

#Join catch polys to pointd
coords_join <- st_join(coords_sf, polys)
cdf<-coords_join %>%
  st_transform(coords_sf, crs = 4326) %>% 
  mutate(Longitude = st_coordinates(.)[,1],
         Latitude = st_coordinates(.)[,2]) %>% st_drop_geometry() %>%  # Remove geometry column
     as.data.frame()  # Convert to data frame
write.csv(cdf, "data/confidential/catchment_area_for_each_sample.csv",row.names = F)
v<-read.csv("data/confidential/catchment_area_for_each_sample.csv")

# # Select and rename the catchment column
# # Bind Latitude and Longitude from original coords
# # coords_catch <- coords_join %>%
# #   st_drop_geometry() %>%  # Remove geometry column
# #   as.data.frame() %>%     # Convert to data frame
# #   select(catchment = S_CatNAME) %>%  # Select and rename the catchment column
# #   bind_cols(coords %>% select(Latitude, Longitude))  # Add Latitude and Longitude columns
# 
# coords_catch <- coords %>%
#   bind_cols(coords_join %>% st_drop_geometry() %>% dplyr::select(catchment = S_CatNAME))
# 
# catch_counts<- coords_catch %>%
#   group_by(catchment) %>%   # Group by the catchment column
#   summarise(count = n(),    # Count the number of entries in each catchment
#             .groups = 'drop') # Drop the grouping structure
# 
# # If not, transform it to EPSG: 27700
# coords_catch_sf <- coords_catch %>%
#   st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
#   st_transform(crs = 27700)
# 
# sw<-scotland_polygon%>%st_transform(crs = 4326)
# cs<-catch_sum %>%  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
# library(ggplot2)
# library(ggrepel)
# ggplot() +
#   geom_sf(data = sw , fill = NA, colour = "black") +  # Plot Scotland outline
#   geom_sf(data = cs, aes(colour = catchment), size = 2, alpha = 0.7) +  # Plot points coloured by catchment
#    geom_label_repel(
#     data = catch_sum, 
#     aes(x = Longitude, y = Latitude, label = catchment),
#     size = 3, 
#     colour = "black", 
#     box.padding = 10,  # Increased padding around the labels
#     point.padding = 0,  # Padding around the points
#     segment.size = 0.5,  # Thicker lines connecting labels
#     segment.color = 'grey50', 
#     max.overlaps = Inf
#   ) +
#   scale_colour_viridis_d() +  # Use a colour scale that's distinct for categorical data
#   theme_minimal() +           # Use a minimal theme for clarity
#   labs(title = "Distribution of Points by Catchment", colour = "Catchment")  # Add a title and colour legend label
# ggplot() +
#   geom_sf(data = sw, fill = NA, colour = "black") +  # Plot Scotland outline
#   geom_sf(data = cs, aes(colour = catchment), size = 2, alpha = 0.7) +  # Plot points coloured by catchment
#   geom_label_repel(
#     data = catch_sum, 
#     aes(x = Longitude, y = Latitude, label = catchment),
#     size = 3, 
#     colour = "black", 
#     box.padding = 5,  # Increased padding around the labels
#     point.padding = 0,  # Padding around the points
#     segment.size = 0.5,  # Thicker lines connecting labels
#     segment.color = 'grey50', 
#     max.overlaps = Inf,
#     nudge_x = 0.2,  # Optional: Fine-tune label positions
#     nudge_y = 0.2,
#     min.segment.length = 0  # Allow lines to be very short
#   ) +
#   scale_colour_viridis_d() +  # Use a colour scale that's distinct for categorical data
#   theme_minimal() +           # Use a minimal theme for clarity
#   guides(colour = "none") + # Removes the legend for colour
#   labs(title = "Distribution of Points by Catchment", colour = "Catchment")  # Add a title and colour legend label
# ggplot() +
#   geom_sf(data = sw, fill = NA, colour = "black") +  # Plot Scotland outline
#   geom_sf(data = cs, aes(colour = catchment), size = 2, alpha = 0.7) +  # Plot points coloured by catchment
#   geom_label_repel(
#     data = catch_sum, 
#     aes(x = Longitude, y = Latitude, label = catchment),
#     size = 3, 
#     colour = "black", 
#     box.padding = 0,  # Increased padding around the labels
#     point.padding = 0,  # Padding around the points
#     segment.size = 0.5,  # Thicker lines connecting labels
#     segment.color = 'grey50', 
#     max.overlaps = Inf,
#     force = 10,
#     nudge_x = 0.5,  # Move labels to the right
#     nudge_y = 0.5,  # Adjust vertical position if needed
#     min.segment.length = 1000  # Allow lines to be very short
#   ) +
#   scale_colour_viridis_d() +  # Use a colour scale that's distinct for categorical data
#   theme_minimal() +           # Use a minimal theme for clarity
#   guides(colour = "none") +  # Removes the legend for colour
#   labs(title = "Distribution of Points by Catchment", colour = "Catchment")  # Add a title and colour legend label
# 
# 
# # Merge polygons with sample counts
# polys_with_counts <- merge(polys, catch_counts, by.x = "S_CatNAME", by.y = "catchment", all.x = TRUE)%>% st_intersection(scotland_polygon)%>% rmapshaper::ms_simplify(keep = 0.005)
# 
# 
# palette1 <- c("#f2f0f7", "#dadaeb", "#bcbddc", "#9e9ac8", "#756bb1", "#54278f")
# palette2 <- c("#bcbddc", "#756bb1", "#54278f")
# theme_clean <- function(base_size = 12) {
#   require(grid)
#   theme_grey(base_size) %+replace%
#     theme(
#       axis.title = element_blank(),
#       axis.text = element_blank(),
#       panel.background = element_blank(),
#       plot.background = element_blank(), 
#       legend.key=element_blank(),
#       panel.grid = element_blank(),
#       axis.ticks.length = unit(0,"cm"), 
#       # axis.ticks.margin = unit(0,"cm"),
#       # panel.margin = unit(0,"lines"),
#       panel.spacing = unit(c(0, 0, 0, 0), "lines"),
#       complete = TRUE
#     )}
# 
# # Create plot
# ggplot() +
#   geom_sf(data = polys_with_counts, aes(fill = count), alpha = 1) +
#   geom_sf_text(data = polys_with_counts, aes(label = count), size = 3, color = "black") +
#   scale_fill_gradientn(
#     colors = palette2,
#     limits = c(0, 90),
#     na.value = "transparent",
#     name = "Number of samples in catchment",
#     guide = guide_legend(
#       keyheight = unit(3, units = "mm"),
#       keywidth = unit(12, units = "mm"),
#       label.position = "bottom",
#       title.position = 'top',
#       nrow = 1
#     )
#   ) +
#   labs(caption = "Number of samples taken for each vice-county") +
#   theme_clean() +
#   theme(
#     text = element_text(color = "#22211d"),
#     plot.title = element_text(size = 22, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
#     plot.subtitle = element_text(size = 17, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.43, l = 2, unit = "cm")),
#     plot.caption = element_text(size = 12, color = "#4e4d47", margin = margin(b = 0.3, r = -99, unit = "cm")),
#     legend.position.inside = c(-0.35, 0.9)
#   )
