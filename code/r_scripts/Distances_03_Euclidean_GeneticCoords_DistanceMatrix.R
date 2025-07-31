### distance between coordinates, modified from Nick Jeff Eelgrass distance script
library(sf)


coords<-read.csv( "data/confidential/meta_for_filtered_genomes.csv")
coords_grouped<-coords%>% group_by(SiteCode) %>%
  summarise(Latitude = first(Latitude), Longitude = first(Longitude))

latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"

# Convert to sf object
site_coords <- coords_grouped %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = latlong, remove = FALSE)
# Calculate "as the bird would fly" distances
geo_dist <- site_coords %>%
  st_distance() %>%
  as.matrix() / 1000 # Convert to kilometres

geo_dist <- geo_dist %>% units::drop_units()

rownames(geo_dist) <- site_coords %>% pull(SiteCode)
colnames(geo_dist) <- site_coords %>% pull(SiteCode)
geo_dist
plot(geo_dist)
saveRDS(geo_dist, file = "data/confidential/distance_matrix.rds")

geo_dist_long <- as.data.frame(as.table(geo_dist))

# Plot the heatmap
 ggplot(geo_dist_long, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(x = "Site", y = "Site", fill = "Distance (km)", title = "Geographic Distance Heatmap") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 
 
 
 total_distances <- geo_dist_long %>%
   group_by(Var1) %>%
   summarise(total_distance = sum(Freq, na.rm = TRUE)) %>%
   arrange(total_distance)
plot(total_distances)
ggplot(total_distances, aes(x = reorder(Var1, total_distance), y = total_distance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(x = "Location", y = "Total Distance (km)", title = "Total Distances by Location") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Plot the heatmap
# Merge total distances into the long format distance data
geo_dist_long <- geo_dist_long %>%
  left_join(total_distances, by = c("Var1" = "Var1")) %>%
  left_join(total_distances, by = c("Var2" = "Var1"), suffix = c(".Var1", ".Var2"))

# Plot the heatmap reordered by total distance
ggplot(geo_dist_long, aes(x = reorder(Var1, total_distance.Var1), y = reorder(Var2, total_distance.Var2), fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(x = "Site", y = "Site", fill = "Distance (km)", title = "Geographic Distance Heatmap") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

