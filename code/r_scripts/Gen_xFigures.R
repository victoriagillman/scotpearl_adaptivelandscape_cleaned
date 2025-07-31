## Large composite plot
library(ggplot2)
library(ggrepel)
library(ggforce)
library(reshape2)
library(patchwork)
library(ggthemes) 
library(paletteer) 
library(ggspatial)  # For scale bar and north arrow


#Save all data from above for easy replot
# dir.create("output/figures/plot_data_dir", showWarnings = FALSE)
font<-theme(
  text = element_text(family = "Arial"),
  plot.title = element_text(family = "Arial"),
  axis.title = element_text(family = "Arial"),
  axis.text = element_text(family = "Arial"),
  legend.text = element_text(family = "Arial"),
  legend.title = element_text(family = "Arial")
)

## Map scotland plot with labelled loc
crs(polys)
polys_bng <- st_transform(polys, crs = 27700)
labels_sf <- labels %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%  # Convert to sf with WGS84
  st_transform(crs = 27700) %>%  # Transform to British National Grid
  mutate(Longitude = st_coordinates(.)[,1],  # Extract the transformed X coordinates
         Latitude = st_coordinates(.)[,2]) %>%  # Extract the transformed Y coordinates
  st_drop_geometry()  # Drop the geometry to keep it as a data frame with x and y
# Transform CRS to British National Grid
labvc<- ggplot() +
  geom_sf(data = polys_bng[!(polys_bng$VCNAME %in% c("Shetland", "Orkney")),], col="grey1", fill="grey80") + # Plot the polygons
  geom_label_repel(
    data = labels_sf, 
    aes(x = Longitude, y = Latitude, label = AnonSiteCode),
    box.padding = 0.1,       # Reduce the padding around the label box
    point.padding = 0.1,     # Reduce the padding around the point
    segment.color = NA,      # Remove the line segments by setting color to NA
    size = 4, 
    colour = "black", 
    max.overlaps = Inf
  ) +
  labs(y="Latitude",
       y="Longitude")+
  theme_bw()+
  annotation_scale(location = "bl", width_hint = 0.3) + 
  annotation_north_arrow(location = "tr", 
                         #style = north_arrow_fancy_orienteering(),
                         which_north = "true" ,
                         height=unit(0.75, "cm"), width=unit(0.75, "cm")
                         ) 
  # coord_sf(crs = "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")  # Orthographic projection centered on Scotland

labvc


unique(polys$VCNAME)

## Plot just VC names
# Convert labels to sf, transform CRS, and join with polys to find matching polygons
labels_sf_join <- coords %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = 27700) #%>%
  st_join(polys_bng, join = st_within)  # Find polygons containing label points

# Get centroids of filtered polygons
  polys_centroids <- polys_bng %>%
    filter(VCNAME %in% unique(labels_sf_join$VCNAME)) %>%
    st_centroid() %>%
    st_coordinates() %>%
    as.data.frame() %>%
    mutate(VCNAME = polys_bng$VCNAME[polys_bng$VCNAME %in% unique(labels_sf$VCNAME)]) %>%
    rename(Latitude = Y, Longitude = X)

# Plot
vcs <- ggplot() +
  geom_sf(data = polys_bng[!(polys_bng$VCNAME %in% c("Shetland", "Orkney")),], col="grey80", fill="grey80") +
  geom_sf(data = polys_bng %>% filter(VCNAME %in% polys_centroids$VCNAME), col="grey10", fill="grey60") +
  geom_point(data = polys_centroids, aes(x = Longitude, y = Latitude), color = "grey40") +
  geom_label_repel(
    data = polys_centroids,
    aes(x = Longitude, y = Latitude, label = VCNAME),
       size = 3,  max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.5,
  ) +
  theme_bw() +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering())

vcs
ggsave(plot= vcs,
       width = 6.5,
       height = 6,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/vc_map.png")

# Mini map for Shetland and Orkney
# mini_map <- ggplot() +
#   geom_sf(data = polys[polys$VCNAME %in% c("Shetland", "Orkney"), ]) +
#   # coord_sf(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat)) + # Adjust limits as necessary
#   theme_bw() +
#   theme(plot.margin = margin(0, 0, 0, 0))+
#   annotation_scale(location = "bl", width_hint = 0.8)   # Scale bar at bottom left
#   # annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering()) 
# mini_map
# labvc + mini_map + 
#   plot_layout(ncol = 2, widths = c(4, 1))  # Adjust widths as needed

save(polys, labels, coords, file = "output/figures/plot_data_dir/lab_vc_plot.RData")
load("output/figures/plot_data_dir/map_scotland_data.RData")  # Load polygons and labels

# FST
load(file = "../../data/too_large_files/Gen03_FST.RData")

fst_matrix_pl<-ggplot(anon_fst_matrix_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(colour = "grey89",
            width = 1, height = 1) +
  geom_text(aes(label = ifelse(is.na(value), "", sprintf("%.2f", value))), 
            color = "black", size = 2.5) +
  scale_fill_gradient(low = "grey97", high = "#7b66d2", na.value = "transparent", limits = c(0, 1)) + ##8074a8 or 7b66d2
  labs(x = "Population", y = "Population") +
  guides(fill = guide_colourbar(title = bquote(F[ST]))) +
  theme_classic() +
  theme(
    # plot.margin = margin(0, 0, 0, 0),     # Remove margin around the plot
    axis.line = element_blank(),  # Remove axis lines
    axis.title = element_text(size=10),
    axis.text=element_text(size=8), #colour="grey1"
    # panel.grid.major = element_blank(),  # Remove major grid lines
    # panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),       # Remove plot border
    legend.position = c(0.09, 0.8),        # Adjust legend position
    axis.ticks.length = unit(0, "cm")  # Set length of ticks
    # axis.ticks = element_line(size = 0.5) # Adjust tick line size
  )+
  # geom_rect(aes(xmin = 0.5, xmax = max(as.numeric(Var1)) + 0.5, 
  #               ymin = 0.5, ymax = max(as.numeric(Var2)) + 0.5), 
  #           colour = "grey1", fill = NA) #size = 0.5
  geom_segment(aes(x = 0.5, xend = max(as.numeric(Var1)) + 0.5, y = 0.5, yend = 0.5), 
               colour = "grey1") +  # Bottom (x-axis)
  geom_segment(aes(x = 0.5, xend = 0.5, y = 0.5, yend = max(as.numeric(Var2)) + 0.5), 
               colour = "grey1")  # Left (y-axis)
  # font
fst_matrix_pl
# save(anon_fst_matrix_melted, labels, file = "output/figures/plot_data_dir/anon_fst_matrix_melted.RData")



# FST vs DIST
pl_fish<-ggplot(merged_data_fish, aes(x = fish_dist, y = Linearised_Fst)) +
  geom_point(size=3, alpha=0.6, colour="darkblue") +  # Scatter plot of points
  # geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear regression line
  labs(#title = bquote("Fish Distance vs " ~ F[ST]),
    x = "Least-Cost-Path Distance (km)", y = expression(paste(F[ST], "/ (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 2)) +  # Set y-axis limits+
  theme_classic()+
  theme(text = element_text(size = 10)) #+font

pl_fish
pl_bird<-ggplot(merged_data_fish, aes(x = Distance, y = Linearised_Fst)) +
  geom_point(size=3, alpha=0.6, colour="darkred") +  # Scatter plot of points
  # geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear regression line
  labs(#title = bquote("Bird Distance vs " ~ F[ST]), 
    x = "Euclidean Distance (km)", y = expression(paste(F[ST], "/ (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 2)) +  # Set y-axis limits+
  theme_classic()+
  theme(text = element_text(size = 10))#+font

pl_bird
# save(merged_data_fish, labels, file = "output/figures/plot_data_dir/merged_data_fish.RData")
pl_env<-ggplot(ibe_data, aes(x = Env_Dist, y = Linearised_Fst)) +
  geom_point(size = 3, alpha = 0.6, colour = "darkgreen") +
  # geom_smooth(method = "lm", se = FALSE, colour = "darkgreen") +
  labs(x = paste0("Environmental Distance (PC1 = ", round(pc1_var, 2), "%)"),
    y = expression(paste(F[ST], " / (1 - ", F[ST], ")"))
  )+
  scale_y_continuous(limits = c(0, 2)) +  # Set y-axis limits+
  theme_classic()+
  theme(text = element_text(size = 10)) #+font
pl_env



library(patchwork)
combi_pl <- pl_fish | pl_bird 
combi_pl_s<-combi_pl+plot_annotation(
  tag_levels = 'a',  # Labels start with 'a', 'b', etc.
  tag_suffix = ')'   # Suffix for labels, e.g., a)
)
combi_pl_s

# Admix plot
generate_admix_plot <- function(Knumber, meta) {
  snmf_meta <- snmf_Qscores_meta(leaproj = pearl_lea, Knum = Knumber, meta = meta)
  snmf_meta <-as.data.frame(snmf_meta)
  qlong <- reshape2::melt(snmf_meta, id.vars = c("ID", "AnonSiteCode", "VCNAME", "SiteCode"))
  # cols <- gg_color_hue(length(unique(qlong$variable)))
  
  ggplot(data = qlong, aes(x = ID, y = value, fill = variable)) +
    geom_bar(stat = "identity", width = 1) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(~AnonSiteCode, scales = "free", switch = "x", space = "free") +
    scale_fill_paletteer_d("ggthemes::Classic_Purple_Gray_12") +
    # scale_fill_manual(values = cols) +
    ylab("Ancestry proportion") +
    # font+
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text = element_text(colour = "black", size = 8,
                                angle=90, hjust = 1),  # Smaller facet text
      strip.background = element_rect(fill = "grey100", color = "grey100"),  # Facet background and border
      panel.border = element_rect(color = "grey1", fill = NA),  # Border around each panel #linewidth = 0.5
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.spacing = unit(0.0, "cm"),
      # legend.position = "none", #or "top"
      legend.position = "top", #or "top"
      legend.direction = "horizontal",  # Horizontal legend
      legend.box = "horizontal",        # Arrange legend in one row
      legend.title = element_blank(),
      legend.text = element_text(size = 8),  # Smaller text size
      legend.key.size = unit(0.5, "cm"),      # Smaller legend keys
      plot.title.position = "plot",          # Title aligned with the legend
      plot.title = element_text(hjust = 0.1, vjust = -10)  # Lower title and center it
    ) +
    guides(fill = guide_legend(nrow = 1)) + 
    ggtitle(paste("K =", Knumber))
}

# cols <- gg_color_hue(9)
# cols<-sample(cols)
opt_k_pl <- generate_admix_plot(best_K, meta_reduced)
opt_k_pl

# save(generate_admix_plot,best_K, meta_reduced, pearl_lea, file = "output/figures/plot_data_dir/admixture_plot_white.RData")

# PCA plot
base_plot <- ggplot(data = PC_K4_meta, aes(x = PC1_genetic, y = PC2_genetic)) +
  labs(x = sprintf("PC1 (%.2f%%)", percentage_var[1]), 
       y = sprintf("PC2 (%.2f%%)", percentage_var[2]),
       fill = "Vice County Anon \nPopulations"
  ) +
  theme_classic() +
  theme(text = element_text(size = 10))
base_plot
ellipsebypop_anon<- base_plot +
  geom_point(aes(fill = AnonSiteCode), size = 5, shape = 21, colour = "black", stroke = 0.5) +
  geom_mark_ellipse(aes(group = AnonSiteCode, fill = AnonSiteCode, label = AnonSiteCode), expand = unit(3, "mm"), 
                    label.fontsize = 8,
                    label.buffer = unit(1, 'mm'),
                    label.fontface = "plain",
                    con.cap = 0
)+
  theme(legend.position = "none")
ellipsebypop_anon
base_plot
# ellipsebypop_anon <- base_plot +
#   geom_point(aes(fill = AnonSiteCode), size = 5, shape = 21, colour = "black", stroke = 0.5) +
#   stat_ellipse(aes(group = AnonSiteCode, fill = AnonSiteCode, colour = AnonSiteCode), 
#                level = 0.95, # 95% confidence interval
#                type = "t",   # Use t-distribution
#                size = 0.8) + # Thickness of ellipse border
#   theme(legend.position = "none")
# save(PC_K4_meta, percentage_var, file = "output/figures/plot_data_dir/pca_ld_data.RData")




# ##========    patchwork      =========
# labvc | fst_matrix_pl | pl_fish | pl_bird |opt_k_pl | ellipsebypop_anon
# 
p1<-labvc
p2<-fst_matrix_pl 
p3<-pl_fish 
p4<-pl_bird
p5<-opt_k_pl
p6<-ellipsebypop_anon +  coord_fixed(ratio = 1)  # 1:1 aspect ratio

layout <- "
AAABBBCCDD
AAABBBCCDD
AAABBBCCDD
AAABBBFFFF
AAABBBFFFF
EEEEEEFFFF
EEEEEEFFFF
EEEEEEFFFF
"
layout <- "
AABBCD
AABBFF
EEEEFF
"


free(p1) + p2 + p3 + p4 + p5 + p6 +
  plot_layout(design = layout)+
  plot_annotation(
    tag_levels = 'a',  # Labels start with 'a', 'b', etc.
    tag_suffix = ')'   # Suffix for labels, e.g., a)
  )

# 
# (p1 | p2) /
#   p5
# (p3 | p4) /
#   p6
# 

scaling_factor<-
#Saving so all right dim
dir.create("output/figures/comp_plot_figs", showWarnings = FALSE)
ggsave(plot= labvc,
       width = 6.5,
       height = 6,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned//output/figures/comp_plot_figs/map.png")
ggsave(plot= fst_matrix_pl,
       width = 6.5,
       height = 6,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/fst_matrix.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/fst_matrix.png")

ggsave(plot= pl_fish,
       width = 3,
       height = 3,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/fst_fish.png")
ggsave(plot= pl_bird,
       width = 3,
       height = 3,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/fst_bird.png")
ggsave(plot= pl_env,
       width = 3,
       height = 3,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/fst_env_dist.png")
ggsave(plot= opt_k_pl, 
       width = 12.5,
       height = 3.2,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/optimumK_snmf_white_background.png")
ggsave(plot= ellipsebypop_anon, 
       width = 6.5,
       height = 6,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/pca_ld.png")

#-----------------   Presentatation plot ------
ggsave(plot= opt_k_pl, 
       width = 8.5,
       height = 3.5,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/pres_optimumK_snmf_white_background.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/pres_optimumK_snmf_white_background.png")

pres_pca<- base_plot +
  geom_point(aes(fill = AnonSiteCode), size = 5, shape = 21, colour = "black", stroke = 0.5) +
  geom_mark_ellipse(aes(group = AnonSiteCode, fill = AnonSiteCode, label = AnonSiteCode), expand = unit(0.1, "in"), 
                    label.fontsize = 14,
                    label.buffer = unit(0, 'in'), #6 was best so far
                    label.fontface = "plain",
                    con.cap = 0,
                    # label.width = 0.1,
                    
  )+
  theme(legend.position = "none")
pres_pca
ggsave(plot= pres_pca, 
       width = 9.5,
       height = 7,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/pres_pca.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/pres_pca.png")
base_plot$data$Host <- factor(base_plot$data$Host, levels = c("Salmon", "Trout", "Only trout available", "Unknown"))
host<-base_plot + 
  geom_point(aes(fill = Host), size = 5, shape = 21, colour = "black", stroke = 0.5)+
  # scale_fill_manual(values = c("Salmon" = "#FA8072", "Trout" = "#A0856A", "Only trout available" = "#C5A78E", "Unknown" = "grey")) +
  labs(fill="Host")+
  theme(legend.position = c(.85, .85))  # Position legend at the top-right

host
ggsave(plot= host, 
       width = 9.5,
       height = 7,
       units = c( "in"),
       dpi = 600,
       filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/pres_host_pca.png")
