### A heatmap of vice counties for genomic analysis
## Better heatmap code- old stuff at the bottom

### A heatmap of vice counties for genomic analysis? Can I whip this up before lunch? 
library(tidyr)
library(sf)
library(dplyr)
library(ggspatial)  # For scale bar and north arrow
library(ggplot2)
library(ggrepel)


polys <- st_read("data/scotland_vice_counties/scotland_vicecounties.shp") %>%
  st_transform(crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% 
  rmapshaper::ms_simplify(keep = 0.002)


# Import the metadata from filtered 
swab_filtered<-read.csv("data/confidential/meta_for_filtered_genomes.csv") %>% group_by(VCNAME) %>% summarise(swabs_filtered = n()) #filt
swab_sequenced<- read.csv("data/confidential/trimmed_env_data.csv") %>%  group_by(VCNAME) %>% summarise(swabs_sequenced = n())#sequenced
swab_all <- read.csv("data/confidential/pearl_meta_updated_with_WI3andWI4.csv") %>%
  mutate(Population = Collection_River) %>%
  mutate(Population = recode(Population,
                             "Bhiaghlann" = "Briaghlann",
                             "Conon_Gharbhrain" = "Conon_Glascarnoch",
                             "Sruthan_Braigh" = "Struthan_Bhraigh")) %>%
  filter(!is.na(Latitude) & Latitude != "na") %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% 
  st_join(., polys) %>% 
  group_by(VCNAME) %>% summarise(swabs_all = n()) %>%
  st_drop_geometry() 


# Temporarily drop geometry to perform joins
swab_counts <- polys %>%
  st_drop_geometry() %>%  # Drop geometry for the joins
  select(VCNAME) %>%      # Keep only the VCNAME column
  full_join(swab_all, by = "VCNAME") %>%
  full_join(swab_filtered, by = "VCNAME") %>%
  full_join(swab_sequenced, by = "VCNAME")

# Reattach geometry using a spatial join
swab_counts <- polys %>%
  select(VCNAME, geometry) %>%  # Select VCNAME and geometry from polys
  left_join(swab_counts, by = "VCNAME") %>%  # Join with the existing swab_counts (with count data)
  st_transform(crs = 27700)  # Transform to the desired CRS (British National Grid)


head(swab_counts)
plot(swab_counts)

ggplot() +
  geom_sf(data = swab_counts[!(swab_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(fill = swabs_all), alpha = 1) +
  geom_sf_text(data = swab_counts[!(swab_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(label = swabs_all), size = 4, color = "white") +
  scale_fill_gradientn(
    colours = c("#ffa2c7", "#d01e79"),  # Gradient from white to deeppink
    na.value = "grey80",
    name = "Number of samples sent for sequencing"
  ) +
  guides(fill = "none") +  # Removes the legend
  theme(
    text = element_text(color = "grey0"),
    plot.title = element_text(size = 22, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
    plot.subtitle = element_text(size = 17, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.43, l = 2, unit = "cm")),
    plot.caption = element_text(size = 12, color = "#4e4d47", margin = margin(b = 0.3, r = -99, unit = "cm")),
    legend.position = "none"  # Removes the legend
  ) +
  labs(x = "Latitude", y = "Longitude") +
  theme_bw() +  # Optional
  annotation_scale(location = "bl", width_hint = 0.3) + 
  annotation_north_arrow(location = "tr", 
                         which_north = "true",
                         height = unit(0.75, "cm"), 
                         width = unit(0.75, "cm"))

## Make nice facet
swab_counts_long <- swab_counts %>%
  pivot_longer(
    cols = c(swabs_all, swabs_filtered, swabs_sequenced),  # Columns to pivot
    names_to = "Swab_Type",  # Name for the new 'type' column
    values_to = "Count"      # Name for the new 'value' column
  ) %>%
  mutate(Swab_Type = factor(Swab_Type,
                            levels = c("swabs_all", "swabs_sequenced", "swabs_filtered"),
                            labels = c("All samples", "Sequenced samples", "Filtered samples")))

# Plot facet
vc_samp_count_facet<-
ggplot(data = swab_counts_long[!(swab_counts_long$VCNAME %in% c("Shetland", "Orkney")),]) +
  geom_sf(data = swab_counts_long[!(swab_counts_long$VCNAME %in% c("Shetland", "Orkney")),], 
          aes(fill = Count), alpha = 1) +
  geom_sf_text(data = swab_counts_long[!(swab_counts_long$VCNAME %in% c("Shetland", "Orkney")),], 
               aes(label = Count), size = 3, color = "white") +
  scale_fill_gradientn(
    colours = c("#ffa2c7", "#d01e79"),  # Gradient from light to deep pink
    na.value = "grey80",
    name = "Individual\nSample Count"
  ) +
  facet_wrap(~Swab_Type) +  # Facet by swab type
  theme(
    text = element_text(color = "grey0"),
    plot.title = element_text(size = 22, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
    plot.subtitle = element_text(size = 17, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.43, l = 2, unit = "cm")),
    plot.caption = element_text(size = 12, color = "#4e4d47", margin = margin(b = 0.3, r = -99, unit = "cm")),
    legend.position = "bottom"  # Legend at the bottom
  ) +
  labs(x = "Latitude", y = "Longitude") +
  theme_bw() +  # Optional for clean background
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tr", 
                         which_north = "true",
                         height = unit(0.75, "cm"), 
                         width = unit(0.75, "cm"))
vc_samp_count_facet # error missing 95, this is because 95 NA values for vice counties
summary(swab_counts_long[!(swab_counts_long$VCNAME %in% c("Shetland", "Orkney")),])

ggsave(plot= vc_samp_count_facet,
       width = 12,
       height = 5,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/map_sample_counts_per_vc_facet.png")

shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/map_sample_counts_per_vc_facet.png")

#### Nice sample vis plot

swab_filtered_sites <- read.csv("data/confidential/meta_for_filtered_genomes.csv") %>%
  group_by(VCNAME) %>%
  summarise(unique_sites = n_distinct(AnonSiteCode), .groups = "drop")

polys_bng<-polys %>% st_transform(crs = 27700) 
labels_sf <- polys_bng %>%
  st_centroid() %>%  # Get centroids
  st_coordinates() %>%  # Extract coordinates
  as.data.frame() %>%  # Convert to data frame
  cbind(polys_bng %>% st_drop_geometry()) %>%   # Add non-geometry data
  filter(VCNAME %in% swab_filtered$VCNAME)  %>%  # Filter relevant rows
  left_join(swab_filtered_sites, by = "VCNAME")%>%
  mutate(Label = paste0(VCNAME, " (N=", unique_sites, ")")) 
# colnames(labels_sf) <- c("Longitude", "Latitude", names(polys %>% st_drop_geometry()))
plot(polys_bng$geometry)
points(labels_sf)

# labels_sf <- labels_sf %>%
  # mutate(Label = paste0(VCNAME, " (N=", unique_sites, ")"))


# # Plot
# ggplot() +
#   geom_sf(data = polys_bng[!(polys_bng$VCNAME %in% c("Shetland", "Orkney")),], 
#           col = "grey1", fill = "grey80") + 
#   geom_point(
#     data = labels_sf,
#     aes(x = X, y = Y),
#     size = 2,  # Adjust point size as needed
#     colour = "black", fill = "deeppink", shape = 21
#   ) +
#   geom_label_repel(
#     data = labels_sf, 
#     aes(x = X, y = Y, label = Label),
#     box.padding = 0.5,       # Padding around label box
#     point.padding = 0.1,     # Padding around point
#     segment.color = "black", # Line segments
#     size = 4, 
#     colour = "black", 
#     max.overlaps = Inf
#   ) +
#   labs(x = "Longitude",
#        y = "Latitude") +
#   theme_bw() +
#   # Scale bar
#   annotation_scale(location = "bl", width_hint = 0.3) + 
#   # North arrow
#   annotation_north_arrow(location = "tr",
#                          which_north = "true",
#                          height = unit(0.75, "cm"), 
#                          width = unit(0.75, "cm"))
# ggplot() +
#   geom_sf(data = swab_counts[!(swab_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(fill = swabs_filtered), alpha = 1) +
#   # geom_sf_text(data = swab_counts[!(swab_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(label = swabs_all), size = 4, color = "white") +
#   scale_fill_gradientn(
#     colours = c("#ffa2c7", "#d01e79"),  # Gradient from white to deeppink
#     na.value = "grey80",
#     name = "Number of filtered samples"
#   ) +
#   geom_point(
#     data = labels_sf,
#     aes(x = X, y = Y),
#     size = 2,  # Adjust point size as needed
#     colour = "black", fill = "deeppink", shape = 21
#   ) +
#   geom_label_repel(
#     data = labels_sf, 
#     aes(x = X, y = Y, label = Label),
#     box.padding = 0.5,       # Padding around label box
#     point.padding = 0.1,     # Padding around point
#     segment.color = "black", # Line segments
#     size = 4, 
#     colour = "black", 
#     max.overlaps = Inf
#   ) +
#   labs(x = "Longitude",
#        y = "Latitude") +
#   theme_bw() +
#   # Scale bar
#   annotation_scale(location = "bl", width_hint = 0.3) + 
#   # North arrow
#   annotation_north_arrow(location = "tr",
#                          which_north = "true",
#                          height = unit(0.75, "cm"), 
#                          width = unit(0.75, "cm"))
# 
# ggplot() +
#   geom_sf(data = swab_counts[!(swab_counts$VCNAME %in% c("Shetland", "Orkney")),], fill="grey80", alpha = 1) +
#   # geom_sf_text(data = swab_counts[!(swab_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(label = swabs_all), size = 4, color = "white") +
#   geom_sf(data = swab_counts[!(swab_counts$VCNAME %in% c("Shetland", "Orkney")) & !is.na(swab_counts$swabs_filtered),],
#           fill = "grey50", alpha = 1, colour = "grey20")+
#   # scale_fill_gradientn(
#   #   colours = c("#ffa2c7", "#d01e79"),  # Gradient from white to deeppink
#   #   na.value = "grey80",
#   #   name = "Number of filtered samples"
#   # ) +
#   geom_point(
#     data = labels_sf,
#     aes(x = X, y = Y),
#     size = 2,  # Adjust point size as needed
#     colour = "black", fill = "deeppink", shape = 21
#   ) +
#   geom_label_repel(
#     data = labels_sf, 
#     aes(x = X, y = Y, label = Label),
#     box.padding = 0.5,       # Padding around label box
#     point.padding = 0.1,     # Padding around point
#     segment.color = "black", # Line segments
#     size = 4, 
#     colour = "black", 
#     max.overlaps = Inf
#   ) +
#   labs(x = "Longitude",
#        y = "Latitude") +
#   theme_bw() +
#   # Scale bar
#   annotation_scale(location = "bl", width_hint = 0.3) + 
#   # North arrow
#   annotation_north_arrow(location = "tr",
#                          which_north = "true",
#                          height = unit(0.75, "cm"), 
#                          width = unit(0.75, "cm"))
# 
# 


###
# Create a label for each VCNAME

labels_with_populations <- read.csv("data/confidential/meta_for_filtered_genomes.csv") %>%
  group_by(VCNAME) %>%
  summarise(
    Populations = paste(sort(unique(AnonSiteCode)), collapse = ", ")  # Sorted populations
  ) %>%
  ungroup() %>%
  arrange(VCNAME) %>%
  mutate(Label = Populations)

# 
# labels_with_populations <- read.csv("data/confidential/meta_for_filtered_genomes.csv")  %>%
#   group_by(VCNAME) %>%
#   summarise(
#     Populations = paste(sort(unique(AnonSiteCode)), collapse = ", ")  # Sort and concatenate unique AnonSiteCodes
#   ) %>%
#   mutate(Label = paste0(VCNAME, ",\nPopulations = ", Populations))  # Format label with line break

# Join the label data to the centroids
labels_sf<- polys_bng %>%
  st_centroid() %>%  # Get centroids
  st_coordinates() %>%  # Extract coordinates
  as.data.frame() %>%  # Convert to data frame
  cbind(polys_bng %>% st_drop_geometry()) %>%   # Add non-geometry data
  left_join(labels_with_populations, by = "VCNAME") %>%  # Join labels
  filter(!is.na(Populations))  # Remove rows with NA Populations

labels_sf <- labels_sf %>% #this is physically placing the labels at coordinates
  mutate(
    Label_X = case_when( #going sideways
      VCNAME == "West Sutherland" ~ 210561.2, # Top middle
      VCNAME =="East Sutherland" ~ 350000, #one clock
      VCNAME =="East Ross" ~ 400000,  # two clock
      VCNAME =="East Inverness-shire"~ 400000, #three clock
      VCNAME =="West Inverness-shire" ~42000, #6 clock
      VCNAME =="North Ebudes" ~ 35000, #7 clock
      VCNAME =="West Ross" ~ 100000, #9 clock
      VCNAME == "Outer Hebrides" ~40000,  # 10 clock
TRUE ~ X  # Default to original X
    ),
    Label_Y = case_when( #going up
      VCNAME == "West Sutherland" ~ 1020000, # Top middle
      VCNAME =="East Sutherland" ~ 1010000, #one o clock
      VCNAME =="East Ross" ~ 920000, #two clock
      VCNAME =="East Inverness-shire"~ 755109.5, #three clock
      VCNAME =="West Inverness-shire" ~ 720000, #6 clock
      VCNAME =="North Ebudes" ~ 834512.9, #7 clock
      VCNAME =="West Ross" ~ 985000, #9 clock
      VCNAME == "Outer Hebrides" ~ 923672.5,  # 10 clock
      TRUE ~ Y  # Default to original Y
    )
  )

# colnames(labels_sf) <- c("X", "Y", names(polys %>% st_drop_geometry()), "Populations", "Label")

# Plot vc labels in speciefied locations
expand<-40000
vc_labels_with_pops<-ggplot() +
  geom_sf(data = polys_bng[!(polys_bng$VCNAME %in% c("Shetland", "Orkney")),], 
          col = "grey1", fill = "grey80") + 
  geom_sf(data = swab_counts[!(swab_counts$VCNAME %in% c("Shetland", "Orkney")) & !is.na(swab_counts$swabs_filtered),],
          fill = "grey50", alpha = 1, colour = "grey1")+
  geom_point(
    data = labels_sf,
    aes(x = X, y = Y),
    size = 2,  # Adjust point size as needed
    colour = "black", fill = "deeppink", shape = 21
  ) +
  geom_segment(data = labels_sf, aes(x = Label_X, y = Label_Y, xend = X, yend = Y), colour = "black") +
  geom_point(data = labels_sf, 
            aes(x = Label_X, y = Label_Y), colour="blue")+
  geom_label(data = labels_sf, 
              aes(x = Label_X, y = Label_Y, label = Label), size=2.5)+
   labs(x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  annotation_scale(location = "bl", width_hint = 0.3) + 
  annotation_north_arrow(location = "tr",
                         which_north = "true",
                         height = unit(0.75, "cm"), 
                         width = unit(0.75, "cm"))+
  coord_sf(xlim = c(7582.734 -expand , 413539.7 +expand ),
           ylim = c(530942, 979277.9+60000))
vc_labels_with_pops
# ggsave(plot= vc_labels_with_pops,
#        width = 6.5,
#        height = 6,
#        units = c( "in"),
#        dpi = 600,
#        filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/map_vc_labels_with_pops.png")
# 
# shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/map_vc_labels_with_pops.png")


vc_labels_with_pops_col_swabs<-ggplot() +
  geom_sf(data = polys_bng[!(polys_bng$VCNAME %in% c("Shetland", "Orkney")),], col="grey80", fill="grey80") +
  geom_sf(data = swab_counts[!is.na(swab_counts$swabs_filtered), ], aes(fill = swabs_filtered), alpha = 1) +
  # geom_sf_text(data = swab_counts[!(swab_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(label = swabs_all), size = 4, color = "white") +
  scale_fill_gradientn(
    colours = c("#ffa2c7", "#d01e79"),  # Gradient from white to deeppink
    na.value = "transparent",
    name = "Number of filtered samples"
  ) +
  geom_segment(data = labels_sf, aes(x = Label_X, y = Label_Y, xend = X, yend = Y), colour = "black", linewidth = 0.3) +
  geom_point(data = labels_sf, 
             aes(x = Label_X, y = Label_Y), colour="blue")+
  geom_label(data = labels_sf, 
             aes(x = Label_X, y = Label_Y, label = Label), size=3, label.r = unit(0, "pt"),label.size = 0,  #fill = NA,
             label.padding = unit(.5, "lines")  # increase padding
)+
  labs(x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  # Scale bar
  annotation_scale(location = "bl", width_hint = 0.3) + 
  # North arrow
  annotation_north_arrow(location = "tr",
                         which_north = "true",
                         height = unit(0.75, "cm"), 
                         width = unit(0.75, "cm"))+
  
  coord_sf(xlim = c(7582.734 -expand , 413539.7 +expand ),
           ylim = c(530942, 979277.9+60000))+
  theme(legend.position  = c(0.72, 0.04),  # x=0.5 (centre), y=0.95 (near top)
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.size = unit(1, "lines"), # smaller or larger to match text
        legend.text = element_text(size = 10)  # adjust number as needed
        # legend.justification = c(0.5, 0.95)  # anchors legend by centre top
  ) 
vc_labels_with_pops_col_swabs
ggsave(plot= vc_labels_with_pops_col_swabs,
       width = 6.5,
       height = 6.0,
       units = c( "in"),
       dpi=600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/map_vc_labels_with_pops_col_swabs.png")

# shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/map_vc_labels_with_pops_col_swabs.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/map_vc_labels_with_pops_col_swabs.png")

library(reshape2)
geo_dist<-readRDS( "data/confidential/distance_matrix.rds")
geo_dist
geo_dist_long <- melt(geo_dist) %>%
  arrange(value) %>%
  filter(as.character(Var1) < as.character(Var2)) 


## For VP for host manuscript: ------------------------------------------------------------------------------------------------
vc_outline_only <- ggplot() +
  geom_sf(
    data = polys_bng[!(polys_bng$VCNAME %in% c("Shetland", "Orkney")),],
    fill = "white",
    colour = "black",
    lwd = 0.3
  ) +
  theme_bw() +  # Clean plot, no axes
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.75, "cm"), width = unit(0.75, "cm"))

# View plot
vc_outline_only
ggsave(plot= vc_outline_only,
       width = 6.5,
       height = 6.0,
       units = c( "in"),
       dpi=600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/vc_just_in_white_no_pop_labels.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/vc_just_in_white_no_pop_labels.png")

labels_sf <- labels_sf %>% #this is physically placing the labels at coordinates
  mutate(
    Label_X = case_when( #going sideways
      VCNAME == "West Sutherland" ~ 210561.2, # Top middle
      VCNAME =="East Sutherland" ~ 350000, #one clock
      VCNAME =="East Ross" ~ 400000,  # two clock
      VCNAME =="East Inverness-shire"~ 442000, #three clock
      VCNAME =="West Inverness-shire" ~42000, #6 clock
      VCNAME =="North Ebudes" ~ 27000, #7 clock
      VCNAME =="West Ross" ~ 100000, #9 clock
      VCNAME == "Outer Hebrides" ~40000,  # 10 clock
      TRUE ~ X  # Default to original X
    ),
    Label_Y = case_when( #going up
      VCNAME == "West Sutherland" ~ 1020000, # Top middle
      VCNAME =="East Sutherland" ~ 1010000, #one o clock
      VCNAME =="East Ross" ~ 920000, #two clock
      VCNAME =="East Inverness-shire"~ 755109.5, #three clock
      VCNAME =="West Inverness-shire" ~ 720000, #6 clock
      VCNAME =="North Ebudes" ~ 834512.9, #7 clock
      VCNAME =="West Ross" ~ 985000, #9 clock
      VCNAME == "Outer Hebrides" ~ 923672.5,  # 10 clock
      TRUE ~ Y  # Default to original Y
    )
  )
expand<-60000

VP_plot<-ggplot() +
  geom_sf(data = polys_bng[!(polys_bng$VCNAME %in% c("Shetland", "Orkney")),], col="grey80", fill="grey80") +
  geom_sf(data = swab_counts[!is.na(swab_counts$swabs_filtered), ], fill = "white", alpha = 1) +
  # geom_sf_text(data = swab_counts[!(swab_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(label = swabs_all), size = 4, color = "white") +
  # scale_fill_gradientn(
  #   colours = c("#ffa2c7", "#d01e79"),  # Gradient from white to deeppink
  #   na.value = "transparent",
  #   name = "Number of filtered samples"
  # ) +
  geom_segment(data = labels_sf, aes(x = Label_X, y = Label_Y, xend = X, yend = Y), colour = "black", linewidth = 0.3) +
  geom_point(data = labels_sf, 
             aes(x = Label_X, y = Label_Y), colour="blue")+
  geom_label(data = labels_sf, 
             aes(x = Label_X, y = Label_Y, label = VCNAME), size=3, label.r = unit(0, "pt"),label.size = 0,  #fill = NA,
             label.padding = unit(.5, "lines")  # increase padding
  )+
  labs(x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  # Scale bar
  annotation_scale(location = "bl", width_hint = 0.3) + 
  # North arrow
  annotation_north_arrow(location = "tr",
                         which_north = "true",
                         height = unit(0.75, "cm"), 
                         width = unit(0.75, "cm"))+
  
  coord_sf(xlim = c(7582.734 -30000 , 413539.7 +70000 ),
           ylim = c(530942, 979277.9+60000))+
  theme(legend.position  = c(0.72, 0.04),  # x=0.5 (centre), y=0.95 (near top)
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.size = unit(1, "lines"), # smaller or larger to match text
        legend.text = element_text(size = 10)  # adjust number as needed
        # legend.justification = c(0.5, 0.95)  # anchors legend by centre top
  ) 
VP_plot
ggsave(plot= VP_plot,
       width = 6.5,
       height = 6.0,
       units = c( "in"),
       dpi=600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/vc_just_in_white_vicecounty_labels.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/vc_just_in_white_vicecounty_labels.png")

scotland <- polys_bng[!(polys_bng$VCNAME %in% c("Shetland", "Orkney")), ] %>%
  st_union() %>%
  st_sf()  # Convert back to sf object
VP_plot_white_land <- ggplot() +
  
  # Land polygons in white with black borders
  geom_sf(data = scotland, 
          fill = "white") +
  
  # Filtered swab polygons on top
  geom_sf(data = swab_counts[!is.na(swab_counts$swabs_filtered), ], fill = "grey90") +
  
  # # Labels
  # geom_segment(data = labels_sf, aes(x = Label_X, y = Label_Y, xend = X, yend = Y), 
  #              colour = "black", linewidth = 0.2) +
  # geom_label(data = labels_sf, 
  #            aes(x = Label_X, y = Label_Y, label = VCNAME), size = 3, 
  #            label.r = unit(0, "pt"), label.size = 0,
  #            label.padding = unit(.5, "lines")) +
  
  # Axes and theme
  labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 10) +
  
  # Scale bar and north arrow
  annotation_scale(location = "bl", width_hint = 0.3) + 
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.75, "cm"), width = unit(0.75, "cm")) +
  
  coord_sf(xlim = c(7582.734 - 30000, 413539.7 + 70000),
           ylim = c(530942, 979277.9 + 60000)) +
  
  theme(
    legend.position = c(0.72, 0.04),
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key.size = unit(1, "lines"),
    legend.text = element_text(size = 10)
  )

VP_plot_white_land
ggsave(plot= VP_plot_white_land,
       width = 6.5,
       height = 6.0,
       units = c( "in"),
       dpi=600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/vc_greyer_vicecounty_nolabels_whitescotland.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/vc_greyer_vicecounty_nolabels_whitescotland.png")
#############   OLD CODE  #####################################################################################################

library(tidyr)
library(sf)
library(dplyr)
library(ggspatial)  # For scale bar and north arrow
library(ggplot2)


# Import the metadata from filtered 
swabs1<-read.csv("data/confidential/meta_for_filtered_genomes.csv") #filt
swabs2<- read.csv("data/confidential/trimmed_env_data.csv") #sequenced
swabs3<- read.csv("data/confidential/pearl_meta3_with_moidart.csv") #all this is full >350 swab dataset will need to do some fiddling on anon_vc_names file to get the file set up essentially ####
          coords<-swabs3
          coords$Population<-coords$Collection_River
          name_mapping <- c( #some of the names in my coords file dont match those in the genomics data
            "Bhiaghlann" = "Briaghlann",
            "Conon_Gharbhrain" = "Conon_Glascarnoch",
            "Sruthan_Braigh" = "Struthan_Bhraigh"
          )
          coords$Population <- recode(coords$Population, !!!name_mapping)
          coords <- coords %>%
                 filter(!is.na(Latitude) & Latitude != "na")
          coords_sf <- st_as_sf(coords, coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
          polys <- st_read("data/scotland_vice_counties/scotland_vicecounties.shp") %>%
             st_transform(st_crs(coords_sf))
          coords_sf <- st_join(coords_sf, polys)
          coords_vc <- coords_sf %>%
               st_drop_geometry() %>%
               as.data.frame() %>%
               bind_cols(coords %>% dplyr::select(Latitude, Longitude))
          swabs2<-coords_vc ####
          swabs2 %>%count(Population)
          
          swabs2 %>%count(AnonSiteCode)
          
          # Extract unique latitudes from swabs2
          latitudes_swabs2 <- swabs2 %>%
            dplyr::select(Latitude) %>%
            distinct()
          
          # Filter swabs3 for these latitudes
          swabs3 %>%
            filter(Latitude %in% latitudes_swabs2$Latitude) %>%
            dplyr::select(Collection_River, Latitude) %>%
            count(Collection_River)
          
# library(viridis)
palette1 <- c("#f2f0f7", "#dadaeb", "#bcbddc", "#9e9ac8", "#756bb1", "#54278f")
palette2 <- c("#bcbddc", "#756bb1", "#54278f")
theme_clean <- function(base_size = 12) {
  require(grid)
  theme_grey(base_size) %+replace%
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(), 
      legend.key=element_blank(),
      panel.grid = element_blank(),
      axis.ticks.length = unit(0,"cm"), 
      # axis.ticks.margin = unit(0,"cm"),
      # panel.margin = unit(0,"lines"),
      panel.spacing = unit(c(0, 0, 0, 0), "lines"),
      complete = TRUE
    )}


# Function to process swabs data
process_swabs_data <- function(swabs_data) {
  
  # Example of your existing code inside the function
  
  # Import shape polygons for vice counties
  polys <- st_read("data/scotland_vice_counties/scotland_vicecounties.shp") %>% rmapshaper::ms_simplify(keep = 0.002)  # Simplify polygons
  
  # Group by VCNAME and summarise counts
  vc_size <- swabs_data %>% group_by(VCNAME) %>% summarise(n = n())
  
  # Merge polygons with sample counts
  polys_with_counts <- merge(polys, vc_size, by = "VCNAME", all.x = TRUE)
  
  # Create plot
  plot <- ggplot() +
    geom_sf(data = polys_with_counts, aes(fill = n), alpha = 1) +
    geom_sf_text(data = polys_with_counts, aes(label = n), size = 6, color = "black") +
    scale_fill_gradientn(
      colors = palette2,
      limits = c(0, 90),
      na.value = "transparent",
      name = "Number of samples",
      guide = guide_legend(
        keyheight = unit(3, units = "mm"),
        keywidth = unit(12, units = "mm"),
        label.position = "bottom",
        title.position = 'top',
        nrow = 1
      )
    ) +
    labs(caption = "Number of samples taken for each vice-county") +
    theme_clean() +
    theme(
      text = element_text(color = "#22211d"),
      plot.title = element_text(size = 22, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
      plot.subtitle = element_text(size = 17, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.43, l = 2, unit = "cm")),
      plot.caption = element_text(size = 12, color = "#4e4d47", margin = margin(b = 0.3, r = -99, unit = "cm")),
      legend.position.inside = c(-0.35, 0.9)
    )
  
  return(list(polys_with_counts = polys_with_counts, plot = plot))
}

x<-process_swabs_data(swabs2)
polys_with_counts<-x$polys_with_counts
# Create plot

plot <- ggplot() +
  geom_sf(data = polys_with_counts[!(polys_with_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(fill = n), alpha = 1) +
  geom_sf_text(data = polys_with_counts[!(polys_with_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(label = n), size = 4, color = "white") +
  scale_fill_gradientn(
    colors = rev(paste0("grey", round(seq(20, 80, length.out = 4)))),  # Round to integers to avoid errors
    limits = c(0, 60),
    na.value = "transparent",
    name = "Number of samples sent for sequencing",
    breaks = seq(0, 60, by = 15),  # Set breaks at intervals of 15
    guide = guide_legend(
      keyheight = unit(3, units = "mm"),
      keywidth = unit(12, units = "mm"),
      label.position = "bottom",
      title.position = 'top',
      nrow = 1
    )
  ) +
  theme(
    text = element_text(color = "grey0"),
    plot.title = element_text(size = 22, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
    plot.subtitle = element_text(size = 17, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.43, l = 2, unit = "cm")),
    plot.caption = element_text(size = 12, color = "#4e4d47", margin = margin(b = 0.3, r = -99, unit = "cm")),
    legend.position = unit(c(0.9, 0.1), "npc"),  # Absolute position using grid units
  ) +
  labs(x="Latitude",
       y="Longitude")+
  theme_bw()+
  annotation_scale(location = "bl", width_hint = 0.3) + 
  annotation_north_arrow(location = "tr", 
                         #style = north_arrow_fancy_orienteering(),
                         which_north = "true" ,
                         height=unit(0.75, "cm"), width=unit(0.75, "cm"))
plot

plot <- ggplot() +
  geom_sf(data = polys_with_counts[!(polys_with_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(fill = as.factor(n)), alpha = 1) +
  geom_sf_text(data = polys_with_counts[!(polys_with_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(label = n), size = 4, color = "white") +
  scale_fill_grey(na.value = "transparent",start = 0.8, end = 0.2)+
   guides(fill = "none") +  # Removes the legend
    theme(
    text = element_text(color = "grey0"),
    plot.title = element_text(size = 22, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
    plot.subtitle = element_text(size = 17, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.43, l = 2, unit = "cm")),
    plot.caption = element_text(size = 12, color = "#4e4d47", margin = margin(b = 0.3, r = -99, unit = "cm")),
    legend.position = "none"  # Removes the legend
  ) +
  labs(x = "Latitude", y = "Longitude") +
  theme_bw() +  # You can remove this if necessary
  annotation_scale(location = "bl", width_hint = 0.3) + 
  annotation_north_arrow(location = "tr", 
                         which_north = "true",
                         height = unit(0.75, "cm"), 
                         width = unit(0.75, "cm"))

plot


plot <- ggplot() +
  geom_sf(data = polys_with_counts[!(polys_with_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(fill = n), alpha = 1) +
  geom_sf_text(data = polys_with_counts[!(polys_with_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(label = n), size = 4, color = "white") +
  # scale_fill_gradientn(
  #   colors = rev(paste0("grey", round(seq(20, 80, length.out = 4)))),  # Round to integers to avoid errors
  #   limits = c(0, 60),
  #   na.value = "transparent",
  #   name = "Number of samples sent for sequencing",
  #   breaks = seq(0, 60, by = 15)  # Set breaks at intervals of 15
  # ) +
  scale_color_distiller(type = "seq",
                        direction = -1,
                        palette = "Greys", 
                        na.value = "transparent")+
  guides(fill = "none") +  # Removes the legend
  theme(
    text = element_text(color = "grey0"),
    plot.title = element_text(size = 22, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
    plot.subtitle = element_text(size = 17, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.43, l = 2, unit = "cm")),
    plot.caption = element_text(size = 12, color = "#4e4d47", margin = margin(b = 0.3, r = -99, unit = "cm")),
    legend.position = "none"  # Removes the legend
  ) +
  labs(x = "Latitude", y = "Longitude") +
  theme_bw() +  # You can remove this if necessary
  annotation_scale(location = "bl", width_hint = 0.3) + 
  annotation_north_arrow(location = "tr", 
                         which_north = "true",
                         height = unit(0.75, "cm"), 
                         width = unit(0.75, "cm"))

plot
plot <- ggplot() +
  geom_sf(data = polys_with_counts[!(polys_with_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(fill = n), alpha = 1) +
  geom_sf_text(data = polys_with_counts[!(polys_with_counts$VCNAME %in% c("Shetland", "Orkney")),], aes(label = n), size = 4, color = "white") +
  scale_fill_gradientn(
    colours = c("#ffa2c7", "#d01e79"),  # Gradient from white to deeppink
    na.value = "grey80",
    name = "Number of samples sent for sequencing"
  ) +
  guides(fill = "none") +  # Removes the legend
  theme(
    text = element_text(color = "grey0"),
    plot.title = element_text(size = 22, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
    plot.subtitle = element_text(size = 17, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.43, l = 2, unit = "cm")),
    plot.caption = element_text(size = 12, color = "#4e4d47", margin = margin(b = 0.3, r = -99, unit = "cm")),
    legend.position = "none"  # Removes the legend
  ) +
  labs(x = "Latitude", y = "Longitude") +
  theme_bw() +  # Optional
  annotation_scale(location = "bl", width_hint = 0.3) + 
  annotation_north_arrow(location = "tr", 
                         which_north = "true",
                         height = unit(0.75, "cm"), 
                         width = unit(0.75, "cm"))

plot

# List of swabs data frames
swabs_list <- list(swabs1, swabs2, swabs3)
# Apply the function to each swabs data frame
plots <- lapply(swabs_list, process_swabs_data)
plots

# Display individual plots
names(plots) <- c("Filtered swabs count Plot", "Sequenced swabs count Plot", "All swabs Plot")

# Save each plot
for (i in seq_along(plots)) {
  ggsave(paste0("output/heatmap_plots/",names[i], ".png"), plots[[i]], width = 9.8, height = 8.7, units = "in", dpi = 600)
}



#Get table vc counts
vc_size <- swabs %>% group_by(VCNAME) %>% summarise(n = n()) %>% print()

# Import shape polygons for vice counties
coords_sf <- st_as_sf(swabs, coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

polys <- st_read("data/scotland_vice_counties/scotland_vicecounties.shp") 
# install.packages("rmapshaper")
polys <- rmapshaper::ms_simplify(polys, keep = 0.005)  # 25% of points retained
coun1 <- merge(polys, vc_size, by = "VCNAME", all.x=TRUE) #####YAYY 
polys <- left_join(polys, vc_size, by = "VCNAME")

###Now make pretty plot ######

# Create the plot
ggplot() +
  geom_sf(data = polys, aes(fill = n), alpha=1) +
  geom_sf_text(data = polys, aes(label = n), size = 6, color = "black") +  # Add text labels was #4e4d47 text colour
  scale_fill_gradientn(
    colors = palette2,
    limits = c(0, 90),
    na.value = "transparent",
    name = "Number of samples",
    guide = guide_legend(
      keyheight = unit(3, units = "mm"),
      keywidth = unit(12, units = "mm"),
      label.position = "bottom",
      title.position = 'top',
      nrow = 1
    )
  ) +
  labs(caption = "Number of samples taken for each vice-county") +
  theme_clean()+
  theme(
    text = element_text(color = "#22211d"),
    plot.title = element_text(size = 22, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
    plot.subtitle = element_text(size = 17, hjust = 0.01, color = "#4e4d47", margin = margin(b = -0.1, t = 0.43, l = 2, unit = "cm")),
    plot.caption = element_text(size = 12, color = "#4e4d47", margin = margin(b = 0.3, r = -99, unit = "cm")),
    legend.position.inside = c(-0.35, 0.9)
  )
all
seq
filt

install.packages("gganimate")
library(gganimate)
plot_list <- list(all, seq, filt)
gganimate::gganimate(
  plot_list,
  interval = 1,  # Adjust as needed (seconds per frame)
  "transition_layers"
)




