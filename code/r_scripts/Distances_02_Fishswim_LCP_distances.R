## Using river dist package to measure lenght of river and grab river mouth coordinates

# Install and load required packages
# install.packages(c("sf", "igraph", "units", "geosphere"))
library(sf)
library(igraph)
library(units)
library(geosphere)
library(mapview)
# Load Scotland outline (polygon) and points data
scotland_polygon <- st_read("data/scotland_vice_counties/scotland_vicecounties.shp") %>% rmapshaper::ms_simplify(keep = 0.001) %>% st_union() %>% st_cast("POLYGON")
points_data <- read.csv("data/confidential/meta_for_filtered_genomes.csv")%>% group_by(AnonSiteCode) %>%
  summarise(Latitude = first(Latitude), Longitude = first(Longitude))%>%
  st_as_sf( coords = c("Longitude", "Latitude"), crs = 4326) %>% st_transform(st_crs(scotland_polygon))

catch <- st_read("data/too_large_files/catchment/MDRExternalStorage.Recordset_12840.shp") %>%
  st_set_crs(st_crs(scotland_polygon))
# plot(catch$geometry)
ggplot() +
  geom_sf(data = catch, fill = NA)
# Load or define your river network data
river_network <- st_read("data/too_large_files/river_lines/just_sampled_rivers_linestring_subset.shp")%>% st_transform(st_crs(scotland_polygon)) 
bbox <- st_bbox(river_network)

head(river_network)
ggplot() +
  geom_sf(data = river_network, aes(colour=identifier)) +
  geom_sf(data = scotland_polygon, fill = NA) +  # Add the scotland_polygon layer
  theme_minimal() +
  guides(colour = "none") +  # Removes the legend for colour
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE)  # Zoom to bbox
# ggplot() +
#   geom_sf(data = long_river_line_sf, aes(colour = factor(segment_id))) +  # Plot river lines
#   geom_sf(data = scotland_polygon, fill = NA) +  # Plot Scotland polygon
#   geom_sf(data = points_data %>% filter(SiteCode == "Briaghlann"), size = 3, colour = "blue") +  # Plot the filtered point
#   points_data %>%
#   filter(SiteCode == "Briaghlann") %>%  # Filter for the desired SiteCode
#   st_buffer(dist = 10000) %>%  # Create a buffer around the selected point (adjust distance as needed)
#   st_bbox() %>%  # Extract the bounding box
#   {coord_sf(xlim = c(.$xmin, .$xmax), ylim = c(.$ymin, .$ymax), expand = FALSE)} +  # Zoom to bbox
#   theme_minimal() +
#   guides(colour = "none")  # Removes the legend for colour

# Extract end points of rivers (estuary points)
riv_multi<- river_network %>%
  group_by(river_id) %>%   # Group by the catchment column
  summarise(count = n(),    # Count the number of entries in each catchment
            .groups = 'drop') # Drop the grouping structure
riv_multi %>%
  filter(river_id %in% c(1748, 1767, 1769)) %>%
  { 
    bbox <- st_bbox(.)
    ggplot() +
      geom_sf(data = ., aes(color = as.factor(river_id))) +
      geom_sf(data = points_data, size = 2) +
      geom_sf_text(data = points_data, aes(label = AnonSiteCode), size = 3, nudge_y = 100, nudge_x = 800) +
      coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
      scale_color_viridis_d(name = "River ID") +
      ggtitle("Subset Rivers with Points (cropped and labelled)") +
      theme_minimal()
  }

library(riverdist)
river_ids <- unique(river_network$river_id) #So now there is a very lengthy process I think of manually changing river id to measure each river
river_ids

filt_riv <- river_network %>%
  filter(river_id == 1769)##1769 is WI4
points_data_filt <- points_data %>%
  filter(AnonSiteCode == "WI4")
ggplot() +
  geom_sf(data = filt_riv, colour="blue", alpha=0.5) +
  geom_sf(data = scotland_polygon, fill = NA) +  # Add the scotland_polygon layer
  geom_sf(data = points_data_filt, colour="red")+
  filt_riv %>%
    st_buffer(dist = 100) %>%  # Create a buffer around the selected point (adjust distance as needed)
    st_bbox() %>%  # Extract the bounding box
    {coord_sf(xlim = c(.$xmin, .$xmax), ylim = c(.$ymin, .$ymax), expand = FALSE)} +  # Zoom to bbox
  # coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE)+  # Zoom to bbox
  theme_minimal()
ggplot() +
  geom_sf(data = filt_riv, colour = "blue", alpha = 0.5) +
  geom_sf(data = scotland_polygon, fill = NA) +
  geom_sf(data = points_data_filt, colour = "red") +
  {  # Insert coord_sf() correctly as a layer
    bbox <- filt_riv %>%
      st_buffer(dist = 1000) %>%
      st_bbox()
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE)
  } +
  theme_minimal()

st_write(filt_riv, "data/too_large_files/oneriver_river_network.shp", delete_layer = TRUE)
river_network_riverdist <- line2network(path = "data/too_large_files", layer = "oneriver_river_network")
# # Perform cleanup if necessary (visual gaps/overlaps)
river_network_riverdist <- cleanup(river_network_riverdist) #(y, y (100))
par(mfrow = c(1, 1))

# filtered_river_network_riverdist <- setmouth(seg = 1, vert = 1, rivers = filtered_river_network_riverdist) #hopeuflly did that in cleanup()
coords<-st_coordinates(points_data_filt)
points_on_river <- xy2segvert(x = coords[, "X"], y = coords[, "Y"], river = river_network_riverdist)
mouthdist(seg = points_on_river$seg,
  vert = points_on_river$vert,
  rivers = river_network_riverdist
) #WI4 1201.622
river_network_riverdist$mouth$mouth.seg

# Plot the river network with highlighted points
riverpoints(
  seg = c(points_on_river$seg, river_network_riverdist$mouth$mouth.seg), 
  vert = c(points_on_river$vert, river_network_riverdist$mouth$mouth.vert), 
  col = "green",  # Colour of the points
  pch = 15, # Plotting character (e.g., filled squares)
  rivers = river_network_riverdist  # The river network object
)
riverdistance(
  startseg = points_on_river$seg,
  endseg = river_network_riverdist$mouth$mouth.seg,
  startvert = points_on_river$vert,
  endvert = river_network_riverdist$mouth$mouth.vert,
  rivers = river_network_riverdist,
  map = TRUE
) #### 26538.6 meters (WI4=1201.622)

###################################### LETS TRY A FUNCTION
#Initialize an empty data frame
df_river_data <- data.frame(SiteCode = character(),
                            Mouth_Lat = numeric(),
                            Mouth_Long = numeric(),
                            Distance_to_mouth = numeric(),
                            stringsAsFactors = FALSE)
list_mouthvertsegs<- data.frame(
  river_id=numeric(),
  SiteCode=character(),
  mouth_segment = numeric(),
  mouth_vertex=numeric()
)
# output_dir <- "data/too_large_files/processed_rivers"  # Directory to save the cleaned river networks

river_ids
length(river_ids)
target_river_id <-river_ids[15] ##1769 is WI4
target_river_id
# Apply the function to each river_id
# dir.create(output_dir, showWarnings = FALSE)
# process_river <- function(river_id) {
  # Filter river network
  filt_riv <- river_network %>% filter(river_id == target_river_id)
  points_data_filt<-points_data %>%
    st_transform(st_crs(filt_riv)) %>%  # Ensure CRS matches
    st_intersection(
      st_buffer(
        st_transform(filt_riv, st_crs(points_data)),  # Ensure CRS matches
        dist = 100
      )
    ) %>% 
    group_by(AnonSiteCode) %>%  # Replace 'SiteCode' with your geographic attribute
    summarise(
      geometry = st_union(geometry)  # Optionally merge geometries into a single feature per group
    ) 
  head(points_data_filt)
  ggplot() +
    geom_sf(data = filt_riv, colour="blue", alpha=0.5) +
    geom_sf(data = scotland_polygon, fill = NA) +  # Add the scotland_polygon layer
    geom_sf(data = points_data_filt, colour="red")+
    filt_riv %>%
    st_buffer(dist = 10000) %>%  # Create a buffer around the selected point (adjust distance as needed)
    st_bbox() %>%  # Extract the bounding box
    {coord_sf(xlim = c(.$xmin, .$xmax), ylim = c(.$ymin, .$ymax), expand = FALSE)} +  # Zoom to bbox
    # coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE)+  # Zoom to bbox
    theme_minimal()
  # Create river network object
  river_network_riverdist <- line2network(filt_riv)
  
  # Perform cleanup
  river_network_riverdist_cl<-cleanup(river_network_riverdist) #yes,yes, min dist 500m, id mouth, dont remove additional segments, build segment route
  

  # Extract the SiteCode
  site_code <- points_data_filt$AnonSiteCode[1]
  # Save the cleaned river network using SiteCode in the filename
  output_file <- file.path(output_dir, paste0("river_network_riverdist_", site_code, ".rds"))
  saveRDS(river_network_riverdist_cl, output_file)
  assign(paste0("river_cleaned_riverdist_", site_code), river_network_riverdist_cl)
  
  # river_network_riverdist<-river_network_riverdist %>% removeduplicates() %>% dissolve() %>% removemicrosegs() %>% splitsegments() %>% addverts( mindist = 1000)
  # plot(river_network_riverdist_cl)
  
  # Convert coordinates to segments and vertices
  coords <- st_coordinates(points_data_filt)
  points_on_river <- xy2segvert(x = coords[, "X"], y = coords[, "Y"], river = river_network_riverdist_cl)
  
  # Set river mouth (assuming you need this step)
  # river_network_riverdist <- setmouth(seg = 2, vert = 2, rivers = river_network_riverdist) # Adjust as necessary
  
plot(river_network_riverdist_cl)
riverpoints(points_on_river$seg, points_on_river$vert, river_network_riverdist_cl, pch=16, col="blue",cex=3)
riverpoints(river_network_riverdist_cl$mouth$mouth.seg, river_network_riverdist_cl$mouth$mouth.vert, river_network_riverdist_cl, pch=15, col="green",cex=3) 
  # Calculate distance to mouth
  distances <- mouthdist(seg = points_on_river$seg, vert = points_on_river$vert, rivers = river_network_riverdist_cl)
  # river_mouth_coords <- river_network_riverdist_cl$mouth
  # mouth_coords <- seg2xy(seg = river_network_riverdist_cl$mouth$mouth.seg, vert = river_network_riverdist_cl$mouth$mouth.vert, rivers = river_network_riverdist_cl)
  # river_network_riverdist_cl$
  
    mouth_segment <- river_network_riverdist_cl$mouth$mouth.seg
    mouth_vertex <- river_network_riverdist_cl$mouth$mouth.vert
    mouth_coords <- river_network_riverdist_cl$lines[[mouth_segment]][mouth_vertex, ]
    mouth_coords
    mouth_sf <- st_as_sf(data.frame(x = mouth_coords[1], y = mouth_coords[2]), 
                         coords = c("x", "y"), 
                         crs = st_crs(river_network_riverdist_cl$sf)) # Use the CRS from the original data

    ggplot() +
      geom_sf(data = filt_riv, colour = "blue", alpha = 0.5) +
      geom_sf(data = scotland_polygon, fill = NA) +  # Add the Scotland polygon layer
      geom_sf(data = points_data_filt, colour = "red") +  # Points of interest
      geom_sf(data = mouth_sf, colour = "green", size = 3, shape = 17) +  # River mouth, styled as a green triangle
      filt_riv %>%
      st_buffer(dist = 10000) %>%  # Create a buffer around the selected point (adjust distance as needed)
      st_bbox() %>%  # Extract the bounding box
      {coord_sf(xlim = c(.$xmin, .$xmax), ylim = c(.$ymin, .$ymax), expand = FALSE)} +  # Zoom to bbox
      theme_minimal()  # Minimal theme

    points_on_river$seg[1]
  # Calculate and visualise river distance #will not work if more than one point on river
  river_distance <- riverdistance(
    startseg = points_on_river$seg,
    endseg = river_network_riverdist_cl$mouth$mouth.seg,
    startvert = points_on_river$vert,
    endvert = river_network_riverdist_cl$mouth$mouth.vert,
    rivers = river_network_riverdist_cl,
    map = TRUE
  )
  distances[1]
  points_data_filt$AnonSiteCode[1]
  # Add mouth coordinates and distance to df_river_data
  df_river_data <- rbind(df_river_data, data.frame(
    SiteCode =   "Struthan_Bhraigh",
    # SiteCode =   points_data_filt$AnonSiteCode[1],
    Mouth_Lat = st_coordinates(mouth_sf)[2],  # Y-coordinate (latitude)
    Mouth_Long = st_coordinates(mouth_sf)[1], # X-coordinate (longitude)
    Distance_to_mouth = distances[1] #could do distances
  ))
  list_mouthvertsegs<-rbind(list_mouthvertsegs, data.frame(
    river_id=target_river_id,
    # SiteCode=points_data_filt$AnonSiteCode[1],
    SiteCode =   "Struthan_Bhraigh",
    mouth_segment = 1, #UPDATE THESE FOR EACH RIVER
    mouth_vertex=5
  ))
  list_mouthvertsegs<-na.omit(list_mouthvertsegs)
  df_combined <- df_river_data %>%
    left_join(list_mouthvertsegs, by = "SiteCode") %>% distinct()
  # df_river_data <- rbind(df_river_data, data.frame(
  #   SiteCode =   points_data_filt$SiteCode[2],
  #   Mouth_Lat = st_coordinates(mouth_sf)[2],  # Y-coordinate (latitude)
  #   Mouth_Long = st_coordinates(mouth_sf)[1], # X-coordinate (longitude)
  #   Distance_to_mouth = distances[2] #could do distances
  # ))
  # 
  #   return(list(distances = distances, river_distance = river_distance))
# }
  
  # write.csv(df_combined, "data/confidential/river_mouth_and_distances_to_sample_points.csv", row.names = F)
  # write.csv(list_mouthvertsegs, "data/confidential/river_mouth_rivdist_verticesandsegmentsused.csv", row.names = F)
  df_combined_og<-(read.csv("data/confidential/river_mouth_and_distances_to_sample_points.csv"))
  df_combined
  df_combined_og
  
  # Update with new data
  df_combined_og[df_combined_og$SiteCode == df_combined$SiteCode, ] <- df_combined
  df_combined_og
  # Optional: write updated dataframe to file
  write.csv(df_combined_og, "data/confidential/river_mouth_and_distances_to_sample_points_updateed2025.csv", row.names = FALSE)
  df_combined<-df_combined_og 
  df_combined<-read.csv("data/confidential/river_mouth_and_distances_to_sample_points_updateed2025.csv")
  # I have migrated the working code here and deleted the old commented out code
  ## Least cost path using a mix of Nick Jefferys code on github, Chloe Cargills and my own thrown in! 
  # Create Scotland raster
  # Move points to nearest raster pixel
  # Create transition object
  # Calc least cost path analysis
  # Output that with distance to river and point bump distances
  
  # Load necessary packages
  library(sf)
  library(rmapshaper)
  library(ggplot2)
  library(raster)
  library(fasterize)
  library(rSDM)
  library(gdistance)
  library(reshape2)
  library(dplyr)
  library(viridis)  # Optional: for colour scales
  
  common_crs<- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
  scotland_polygon <- st_read("data/scotland_vice_counties/scotland_vicecounties.shp") %>% rmapshaper::ms_simplify(keep = 0.001) %>% st_union() %>% st_cast("POLYGON")
  scotland_polygon <-st_transform(scotland_polygon,common_crs)
  coords<-read.csv("data/confidential/river_mouth_and_distances_to_sample_points_updateed2025.csv") %>% 
    rename(long = Mouth_Long,
           lat = Mouth_Lat) %>% st_as_sf( coords = c("long", "lat"), crs = 27700) %>% st_transform(st_crs(scotland_polygon))
  river_network <- st_read("data/too_large_files/river_lines/just_sampled_rivers_linestring_subset.shp")%>% st_transform(st_crs(scotland_polygon)) 
  # coords<-points_data %>% left_join(coords,by = "SiteCode")
  ggplot() +
    geom_sf(data = coords) +
    geom_sf(data = scotland_polygon, fill = NA) +  # Add the scotland_polygon layer
    theme_minimal() +
    guides(colour = "none") 
  
  
  bound_box<- st_bbox(c(xmin = 63000 , xmax = 350000 , ymax = 990000, ymin = 740000 ), crs = 27700)%>%
    st_as_sfc()%>%
    st_as_sf()
  
  bound_box<- st_bbox(scotland_polygon, crs = 27700)%>%
    st_as_sfc()%>%
    st_as_sf()
  
  ggplot()+
    geom_sf(data=scotland_polygon)+
    # geom_sf(data=bound_box_small,fill=NA)+
    geom_sf(data=coords)+
    theme_bw()
  library(sf)
  
  # coords = coords%>%st_transform(common_crs)
  plot(coords)
  #create dynamic projection to set the boundaries of the least-cost path analysis
  pts_center <- coords%>%
    summarise()%>%
    st_centroid()%>%
    st_transform(common_crs)%>%
    suppressWarnings()%>%
    suppressMessages()
  pts_center
  aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                  st_coordinates(pts_center)[2], st_coordinates(pts_center)[1])
  
  #set cropping buffer to create a transition object
  
  max_dist <- (max(st_distance(coords)%>%as.numeric()/1000)+25)%>%round() #maximum distance between stations + 25km radius
  max_dist
  buffer_pts <- st_buffer(pts_center%>%st_transform(aeqd),dist=units::set_units(max_dist,"km"))%>%st_bbox()%>%st_as_sfc()%>%st_sf()
  buffer_pts
  plot(st_geometry(buffer_pts))
  plot(st_geometry(coords), add = TRUE, col = 'red')  # Overlay points
  print(st_bbox(buffer_pts))  # Check extent of the buffer
  
  #create raster at 5 km scale that can be used to develop the transition object
  resolution<-1
  r=suppressWarnings(raster(raster::extent(buffer_pts),res=resolution*1000,crs=aeqd)%>%projectRaster(.,crs=common_crs))
  scotland_polygon<-scotland_polygon %>% st_transform(common_crs)
  bound_box<-bound_box %>% st_transform(common_crs)
  coast <- scotland_polygon%>%st_intersection(bound_box)
  plot(coast)
  
  ind <- r%>%
    fasterize(coast,.)%>%
    values(.)%>%
    as.data.frame()%>%
    rename(val=1)%>%
    mutate(val=is.na(val))%>%
    pull(val)
  plot(ind)
  #set raster value for the transition matrix
  r[] <- NA #land - no conductance
  r[ind] <- 1 #water (note set to 1, if you set to 10 you change the costDistance calculation)
  #raster used to ensure the sample sites are in water
  r2 <- r
  r2[] <- 10
  r2[ind] <- NA
  plot(r)
  plot(r2)
  
  
  # Make sure coords in raster
  coords_df<- coords %>%
    st_coordinates() %>%
    as.data.frame() %>%
    bind_cols(st_drop_geometry(coords)) %>%
    rename(longitude = X, latitude = Y)
  
  
  # ## Change raster to fit points BUT leaves a hole in scotland so at least on point will not work!
  # pts <- coords_df %>%
  #   dplyr::select(longitude, latitude) %>%
  #   as.matrix()
  # r3<-r
  # icell <- cellFromXY(r3, pts) # identify cells in which sampling locations occur
  # 
  # # r[icell[1:11]] #check whether any points sit on land: #1, 4, 8, 11 do
  # # r[icell[c(1, 4, 8, 11)]] <- 1 #coerce sampling locations to sea.
  # r3<-r
  # r3[icell]
  # 
  # r3[icell] <- ifelse(is.na(r3[icell]), 1, r3[icell]) #essextially instead of bumping we just make sea anyway
  # r3[icell]
  # plot(r3) #forget about bumping coords, move the sea to points/river mouths
  # points(coords)
  # 
  # plot(r)
  # points(coords)
  # 
  # z.use <- rast(r)
  # z.use[icell[10]] 
  # z.use[icell[10]] <- 2 
  # d <- gridDist(z.use, target = 2, scale=1000) #calculate shortest distance in km to sampling locations, ignoring NA values (land)
  # d[icell] #distances in km from sampling points
  # 
  # plot(d) # check output is logical
  # points(coords)
  
  
  #### Locs2nearestcell
  library(rSDM)
  # Read and prepare sample location data
  coords_df<- coords %>%
    st_coordinates() %>%
    as.data.frame() %>%
    bind_cols(st_drop_geometry(coords)) %>%
    rename(longitude = X, latitude = Y)
  
  # Convert sample locations to sf object
  coords_sf <- locs2sf(coords_df, lon.col = "longitude", lat.col = "latitude", crs = common_crs)
  plot(r)
  points(coords)
  # Move coordinates to nearest cell
  move_coords <- points2nearestcell(
    locs = coords_sf,
    ras = terra::rast(r), #I added a subset here- not sure if this will be good or not? i think I rembmer the program naturally just uses the first layer
    move = TRUE,
    distance = 100000000,
    table = TRUE,
    map = "base"
  )
  
  move_coords_df<- move_coords %>%
    st_coordinates() %>%
    as.data.frame() %>%
    bind_cols(st_drop_geometry(coords)) %>%
    rename(longitude = X, latitude = Y)
  
  plot(r)
  points(coords)
  points(move_coords, col="red")
  
  # Get distances coords were moved
  move_coords$offset=9999 #placeholder for the distances 
  for(i in move_coords$SiteCode){
    
    move_site <- move_coords%>%filter(SiteCode == i)
    coords_site<-coords%>%filter(SiteCode == i)
    
    move_coords[move_coords$SiteCode == i, "offset"] <- st_distance(
      move_site, coords_site
    ) %>% as.numeric() / 1000  # Distance in km
    
    rm(move_site)
    rm(coords_site)
  }
  
  ## create row index that can be used for matching
  move_coords$ind <- 1:nrow(move_coords)
  
  trans <- transition(r,transitionFunction = min,directions=16)%>%
    geoCorrection(.,type="c",multpl = FALSE) #may take time
  
  ## Least-cost path analysis
  # install.packages("gdistance")
  library(gdistance)
  lcp_df <- suppressWarnings(costDistance(trans,as_Spatial(move_coords))/1000)
  
  
  ## Least-cost path lines for plotting
  nb.loc <- nrow(move_coords)
  path <- list()
  comb <- combn(1:nb.loc, 2)
  
  pairwise_lines <- NULL
  pb <- txtProgressBar(min = 0, max = ncol(comb), style = 3)
  for (i in 1:ncol(comb)) {
    
    origin <- move_coords[comb[1, i],]%>%as_Spatial()
    goal <- move_coords[comb[2, i],]%>%as_Spatial()
    temp <- gdistance::shortestPath(trans, origin, goal,
                                    output = "SpatialLines")
    
    #convert to sf
    #if there is more than one point then create a linestring, else create a point
    
    linestring_logic <- length(temp@lines[[1]]@Lines[[1]]@coords)>1
    
    if(linestring_logic){
      
      line_coords <- data.frame(temp@lines[[1]]@Lines[[1]]@coords)
      
      temp_sf <- st_as_sf(line_coords,coords=c("x","y"),crs=common_crs)%>%
        summarise(do_union = FALSE) %>%
        st_cast("LINESTRING") %>%
        mutate(origin = move_coords[comb[1, i],]%>%pull(SiteCode),
               dest = move_coords[comb[2, i],]%>%pull(SiteCode),
               id = paste(origin,dest,sep="-"),
               len = nrow(line_coords))%>%
        dplyr::select(origin,dest,id,len,geometry)
      #}
      
      #combine with i results
      pairwise_lines <- rbind(pairwise_lines,temp_sf)
      
      setTxtProgressBar(pb, i)
      
    } #end for (i in 1:ncol(comb))
  }
  
  # st_write(pairwise_lines, "data/confidential/fishswim_path_lines.shp", driver = "ESRI Shapefile")
  
  # Assuming pairwise_lines is already created and is an sf object
  str(pairwise_lines)
  summary(pairwise_lines)
  st_is_valid(pairwise_lines)
  st_crs(pairwise_lines)
  nrow(pairwise_lines)
  
  # Transform layers to BNG (EPSG:27700)
  scotland_polygon_bng <- st_transform(scotland_polygon, 27700)
  river_network_bng <- st_transform(river_network, 27700)
  pairwise_lines_bng <- st_transform(pairwise_lines, 27700)
  
  # Plot with BNG
  leastcost_plot <- ggplot() +
    geom_sf(data = scotland_polygon_bng, fill = "grey90", colour = "black") +
    geom_sf(data = river_network_bng, colour = "grey70") +
    geom_sf(data = pairwise_lines_bng[pairwise_lines_bng$len > 1, ], aes(colour = origin), linewidth = 1) + 
    scale_colour_viridis_d() + 
    theme_bw() +
    annotation_scale(location = "bl", width_hint = 0.3) + 
    annotation_north_arrow(location = "tl", 
                           which_north = "true",
                           height = unit(0.75, "cm"),
                           width = unit(0.75, "cm")) +
    theme(legend.position = "none") +
    labs(
      x = "Longitude",
      y = "Latitude")
  
  leastcost_plot
ggsave(plot= leastcost_plot,
           width = 6,
           height = 9,
           units = c( "in"),
           dpi = 600,
           filename = "./output/figures/least_cost_plot.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/least_cost_plot.png")
  
  pairwise_lines$id[pairwise_lines$len <= 1] ## these are those that go nowhere because they are in the same river/ catchment so share the same rivermouth to rivermouth dist
  
  ## Lets make a df
  output <- melt(as.matrix(lcp_df), varnames = c("row", "col"))%>%
    mutate(start = move_coords%>%data.frame()%>%.[row,"SiteCode"],
           end = move_coords%>%data.frame()%>%.[col,'SiteCode'])%>%
    rename(dist=value)%>%
    dplyr::select(start,end,dist)
  
  #adjust for the bumps to water
  output_original<-output
  for(i in move_coords$SiteCode){
    output[output$start == i & output$dist!=0 | output$end == i & output$dist!=0,"dist"] <- output[output$start == i & output$dist!=0 | 
                                                                                                     output$end == i & output$dist!=0,"dist"] + move_coords%>%filter(SiteCode == i)%>%pull(offset)
    
    
  }
  output_with_bumptorast<-output
  # 3. Adjust the distances by adding move_coords$Distancetomouth
  move_coords$Distance_to_mouth_km<-move_coords$Distance_to_mouth%>% as.numeric() / 1000  # Distance in km
  for(i in move_coords$SiteCode){
    output[output$start == i & output$dist != 0 | output$end == i & output$dist != 0, "dist"] <- 
      output[output$start == i & output$dist != 0 | output$end == i & output$dist != 0, "dist"] + 
      move_coords %>% filter(SiteCode == i) %>% pull(Distance_to_mouth_km)
  }
  plot(output$dist)
  write.csv(output, "data/confidential/point_distance_matrix_as_fish_swims", row.names=F)
  otest<-read.csv("data/confidential/point_distance_matrix_as_fish_swims_og")
  plot(otest$dist, output$dist)
  identical(otest$dist, output$dist)
  which(otest$dist != output$dist)
  
