### Here I load in a river linestring dataset, give connected linestrings a common river_id name and subset to just rivers the sampled mussels are in

# Load necessary libraries
library(dplyr)
library(igraph)

# Example data frame
scotland_polygon <- st_read("data/scotland_vice_counties/scotland_vicecounties.shp") %>% rmapshaper::ms_simplify(keep = 0.001) %>% st_union() %>% st_cast("POLYGON")
river_segments <- st_read("data/too_large_files/river_lines/River lines clipped.shp") %>% rmapshaper::ms_simplify(keep = 0.01)%>% st_transform(st_crs(scotland_polygon)) # %>% st_intersection(river_network, scotland_polygon) #simplify to make it actually able to process
summary(river_segments)
points_data <- read.csv("data/confidential/meta_for_filtered_genomes.csv")%>% group_by(SiteCode) %>%
  summarise(Latitude = first(Latitude), Longitude = first(Longitude))%>%
  st_as_sf( coords = c("Longitude", "Latitude"), crs = 4326) %>% st_transform(st_crs(scotland_polygon))
linestr <- river_segments[st_geometry_type(river_segments) == "LINESTRING", ]

#Give rivers a unique code
# install.packages("riverdist")
# library(riverdist)
# MyRivernetwork <- line2network(sf=linestr)
#  
# # Convert to spatial format compatible with riverdist
# river_lines <- as(river_segments, "Spatial")
 
# Step 1: Create an edge list (pairs of startpointid and endpointid)
edges <- linestr %>%  st_drop_geometry()%>%
  dplyr::select(startNode, endNode) %>%
  as.matrix()
dim(edges) #checks whether geom has gone otherwise next wont work

# Step 2: Create a graph from the edge list
g <- graph_from_edgelist(edges, directed = FALSE)

# Step 3: Identify connected components (rivers)
components <- clusters(g)

V(g)$name
# Step 4: Create a data frame with the river IDs
river_ids <- data.frame(
  startNode = V(g)$name,
  river_id = components$membership
)

# Step 5: Merge the river_ids back with the original river_segments
river_segments <- linestr %>%
  left_join(river_ids, by = "startNode")

# View the result
print(river_segments)
# Filter for river_id 1
river_1 <- river_segments %>%
  filter(river_id == 2)

ggplot(data = river_1) +
  geom_sf() +
  ggtitle("Geometry of River with ID 1") +
  theme_minimal()

# Match unique code to point data and subset
# Step 1: Find the nearest river segment for each point
nearest_rivers <- st_nearest_feature(points_data, river_segments)

# Step 2: Subset the river segments using the nearest river indices
nearest_river_segments <- river_segments[nearest_rivers, ]

# Step 3: Extract the unique river IDs of the nearest river segments
nearest_river_ids <- unique(nearest_river_segments$river_id)

# Step 4: Filter river_segments to include only those with the nearest river IDs
my_rivers <- river_segments %>%
  filter(river_id %in% nearest_river_ids)
head(my_rivers)
# Step 5: Plot the nearest river segments
bbox <- st_bbox(my_rivers)
ggplot() +
  geom_sf(data = my_rivers, aes(color = as.character(river_id))) +  # Plot nearest river segments
  geom_sf(data = points_data, color = 'red') +  # Optional: Plot points
  geom_sf(data = scotland_polygon, fill = NA, color = 'black', size = 0.5) +  # Coastline as reference
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +  # Set plot limits
  ggtitle("Nearest River Segments to Points") +
  theme_minimal()

riv_counts<- my_rivers %>%
  group_by(river_id) %>%   # Group by the catchment column
  summarise(count = n(),    # Count the number of entries in each catchment
            .groups = 'drop') # Drop the grouping structure
print(riv_counts)

st_write(my_rivers, "data/too_large_files/river_lines/just_sampled_rivers_linestring_subset.shp", delete_layer = TRUE)

