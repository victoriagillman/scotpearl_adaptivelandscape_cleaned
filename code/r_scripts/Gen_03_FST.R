# Plink for FST 

#Set wd to be where .bed .bim .fam files are
setwd("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/data/plink_VP/")
getwd()
load(file = "../../data/too_large_files/Gen03_FST.RData")

# Use system coding to talk to plink- check if its working
system("cmd.exe /c dir")
system("./plink.exe")

system("plink --file LDpruned_ped --fst --family --allow-extra-chr --out FstResults")
# Total genotyping rate is 0.960745.
# 3456 variants and 156 people pass filters and QC.
# Note: No phenotypes present.
# Writing --fst report (18 populations) to FstResults.fst ... done.
# 3427 markers with valid Fst estimates (29 excluded).
# Mean Fst estimate: 0.187523
# Weighted Fst estimate: 0.218283


# evaluate plink output from
library(ggplot2)
library(dplyr)
library(data.table)
if (!require("qqman")) {
  install.packages("qqman", dependencies = TRUE)
  library(qqman)
}

fst <- fread("FstResults.fst")



data <- read.delim("FstResults.fst") %>%
  tidyr::drop_na() 


# Create a mapping from scaffold names to numeric values
scaffold_names <- unique(data$CHR)
scaffold_mapping <- setNames(seq_along(scaffold_names), scaffold_names)

# Replace scaffold names with their corresponding numeric values
data <- data %>%
  mutate(CHR = scaffold_mapping[CHR])

# Check the head of the modified data
head(data)

manhattan(data, chr = "CHR", bp = "POS", p = "FST", logp = F)

global_fst <- mean(data$FST, na.rm = TRUE)
global_fst
###=================================================================================
#Creating a matrix?
# install.packages("StAMPP")
library(StAMPP)
library(adegenet)
# install.packages("ape")
library(ape)
library(dplyr)
# BiocManager::install("ggtree")
library(ggtree)

system("./plink --file LDpruned_ped --allow-extra-chr --recode A --out LDpruned_recodeA")
popgenli<-read.PLINK(file="LDpruned_recodeA.raw")
        ## system("plink --file pearlmusselfiles --recode A --allow-extra-chr --out RAW")
        ## popgenli<-read.PLINK(file="FstRAW.raw") #using NON LD data


# import of Stampp
stamp<-stamppConvert(popgenli, "genlight")

# #Calc pairwise FST val between each pop (takes time)
# pop_fst<-stamppFst(stamp, nboots = 1000, percent = 95, nclusters = 1)
# pop_fst
# saveRDS(pop_fst, "../confidential/pop_fst_1000boots.rds")
# pop_fst<-stamppFst(stamp, nboots = 10000, percent = 99, nclusters = 1)
# saveRDS(pop_fst, "../confidential/pop_fst_10000boots_99perc.rds")

summary(popgenli)
table(pop(popgenli))
pop_fst$Fsts


# ###dfiff method
# genlight_obj <- gl.read.vcf("ld_pruned_het.vcf")
# head(genlight_obj)
# indNames(genlight_obj)
# 
# pop(genlight_obj) <- sapply(indNames(genlight_obj), function(x) strsplit(x, "_")[[1]][1])
# pop(genlight_obj) <- sapply(indNames(genlight_obj), function(x) {
#   parts <- unlist(strsplit(x, "_"))
#   paste(parts[-length(parts)], collapse = "_")
# }) #split pop names from IDs by last _ in name
# 
# pop(genlight_obj)
# dartFst<-gl.fst.pop(popgenli, nboots = 1000, percent = 95, nclusters = 1, verbose = NULL)
# 
# dartFst$Pvalues


pop_fst<-readRDS( "../confidential/pop_fst_1000boots.rds")
str(pop_fst)

# fst.matrix <- as.data.frame(pop_fst$Fsts) #might be the better way but not changing now
# Extract the FST matrix
fst_matrix <- pop_fst$Fsts
max(fst_matrix, na.rm=T)
min(fst_matrix, na.rm=T)

# Extract p values
pvals <- as.data.frame(pop_fst$Pvalues) #all 0
pvals
pvals_adj <- matrix(
  p.adjust(unlist(pop_fst$Pvalues), method = "bonferroni"),
  nrow = nrow(pop_fst$Pvalues),
  dimnames = dimnames(pop_fst$Pvalues)
)
pvals_adj
summary(pvals_adj)
nrow(pvals)


pop_fst$Pvalues
#####   Plot phylogenies ##################
tree<-nj(as.dist(fst_matrix))

# Basic ways to plot
plot(tree) #some negative values, because algorithm struggling to fit the data
plot.phylo(tree, adj=0.5)

# Better plot
# Read the vicec data and handle special cases 
vicec<-read.csv("../../data/confidential/meta_for_filtered_genomes.csv") %>% dplyr::select(-X) %>% mutate(Population = ifelse(Population == "Conon_Gharbhrain_wrong_sample_ID", "Conon_Glascarnoch", Population)) %>% mutate(AnonSiteCode = ifelse(AnonSiteCode == "ER2", "ER1", AnonSiteCode)) %>% dplyr::select(AnonSiteCode, VCNAME, Population, ID)
length(unique(vicec$AnonSiteCode))
# Create a mapping of population to anonymous site code and VICEC group
name_map <- setNames(vicec$AnonSiteCode, vicec$Population)
vicec_group_map <- setNames(vicec$VCNAME, vicec$Population)
anon_df <- unique(vicec%>% dplyr::select(AnonSiteCode, VCNAME, Population))
# write.csv(anon_df, "../../data/confidential/anon_pop_names_reduced.csv", row.names = F)
anon_df<-read.csv("../../data/confidential/anon_pop_names_reduced.csv")

# Update the tree with anonymous labels
anon_tree <- tree
anon_tree$tip.label <- name_map[tree$tip.label]

# Create annotation data frame for anonymous tree
annot_df_anon <- data.frame(label = anon_tree$tip.label, VCNAME = vicec_group_map[tree$tip.label])

# Plot the anonymous tree with annotations
options(ignore.negative.edge=T)
p_anon<-ggtree(anon_tree, aes(size-1)) %<+% annot_df_anon +
  geom_tiplab(aes(label = label, colour = VCNAME), offset = 0.001) + # Adjust offset for better label positioning
  theme(legend.position = "right") +
  # geom_tippoint(size = 3, alpha = 0.5, aes(colour = VCNAME)) +
  geom_treescale() +
  ggplot2::labs(colour = "Vice County") + # Rename the legend
  ggplot2::guides(colour = guide_legend(override.aes = list(label = ""))) # Show only dots in the legend
p_anon

p_anon<-ggtree(anon_tree, aes(colour = VCNAME), linewidth=1) %<+% annot_df_anon +
  # geom_tree()+
  geom_tiplab(aes(label = label, colour = VCNAME), offset = 0.001) +
  # geom_tippoint(size = 3, alpha = 0.5, aes(colour = VCNAME)) +
  geom_treescale() +
  coord_cartesian(clip = "off")+
  theme(legend.position = "right") +
  ggplot2::labs(colour = "Vice County") +
  theme(legend.position = "none", # Remove legend
        plot.margin = margin(5.5, 40, 5.5, 5.5)) 
  # guides(colour = guide_legend(override.aes = list(label = "")))
p_anon
# ggsave(plot= p_anon, #orininal size
#        width = 700/100,
#        height = 500/100,
#        units = c( "in"),
#        dpi = 600,
#        filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/anon_njtree_truenegativevals.png")
# shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/figures/anon_njtree_truenegativevals.png")

# ggsave(plot= p_anon, #updated for comp plot
#        width = 12,
#        height = 6,
#        units = c( "in"),
#        dpi = 600,
#        filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/anon_njtree_truenegativevals_dim12vs6.png")
# shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/figures/comp_plot_figs/anon_njtree_truenegativevals_dim12vs6.png")
# 
# ggsave(plot= p_anon, #updated for comp plot
#        width = 6.5,
#        height = 6,
#        units = c( "in"),
#        dpi = 600,
#        filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/comp_plot_figs/anon_njtree_truenegativevals_dim6.5vs6.png")
# shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/figures/comp_plot_figs/anon_njtree_truenegativevals_dim6.5vs6.png")


# Create annotation data frame for original tree (non-anonymous)
annot_df_original <- data.frame(label = tree$tip.label, VCNAME = vicec_group_map[tree$tip.label])

# Plot the original tree with annotations
p_original <- ggtree(tree) %<+% annot_df_original +
  geom_tiplab(aes(label = label, colour = VCNAME), offset = 0.001) + # Adjust offset for better label positioning
  theme(legend.position = "right") +
  geom_tippoint(size = 3, alpha = 0.5, aes(colour = VCNAME)) +
  geom_treescale() +
  ggplot2::labs(colour = "Vice County") + # Rename the legend
  ggplot2::guides(colour = guide_legend(override.aes = list(label = ""))) # Show only dots in the legend

# Display both plots
p_anon  # Display anonymous tree plot
p_original  # Display original tree plot

sum(anon_tree$edge.length < 0, na.rm = TRUE)  # Count negative edges
anon_tree$edge.length[anon_tree$edge.length < 0] <- 0  # Set negatives to zero


## FST matrix plot #######################

# Assuming fst_matrix is already prepared and cleaned
# If fst_matrix is not in the correct format, reshape it using melt function
library(reshape2)
vicec<-read.csv("../../data/confidential/trimmed_env_data.csv")%>%
  mutate(Sample_ID = ifelse(Sample_ID == "ICY837","ICY837notICY773", Sample_ID)) %>% mutate(Population = ifelse(Population == "Conon_Gharbhrain_wrong_sample_ID", "Conon_Glascarnoch", Population))%>% mutate(AnonSiteCode = ifelse(AnonSiteCode == "ER2", "ER1", AnonSiteCode))
# vicec<-read.csv("../../data/confidential/trimmed_env_data.csv") %>% dplyr::select(-X) %>% mutate(Population = ifelse(Population == "Conon_Gharbhrain_wrong_sample_ID", "Conon_Glascarnoch", Population))
# Create a mapping of population to anonymous site code and VICEC group
name_map <- setNames(vicec$AnonSiteCode, vicec$Population)
vicec_group_map <- setNames(vicec$VCNAME, vicec$Population)

print(name_map)

# Rename row and column names based on name_map
# fst_matrix <- fst_matrix[order(rownames(fst_matrix)), order(colnames(fst_matrix))]
anon_fst_matrix<-fst_matrix
rownames(anon_fst_matrix) <- name_map[rownames(anon_fst_matrix)]
colnames(anon_fst_matrix) <- name_map[colnames(anon_fst_matrix)]

# Print renamed matrix
length(rownames(anon_fst_matrix))
print(anon_fst_matrix)
summary(anon_fst_matrix)
print(fst_matrix)

# Symmetrize the matrix: Fill NA values by mirroring across the diagonal
anon_fst_matrix[upper.tri(anon_fst_matrix)] <- t(anon_fst_matrix)[upper.tri(anon_fst_matrix)]

# #Make alphabetical 
ordered_anon_names <- sort(rownames(anon_fst_matrix))  # Get the anon names in alphabetical order
anon_fst_matrix <- anon_fst_matrix[ordered_anon_names, ordered_anon_names]

# Remove the upper triangle by setting it to NA
anon_fst_matrix[upper.tri(anon_fst_matrix)] <- NA

# Melt matrices for plotting
fst_matrix_melted <- melt(as.matrix(fst_matrix))
anon_fst_matrix_melted <- melt(as.matrix(anon_fst_matrix))
summary(anon_fst_matrix_melted)
# # Plot non-anonymous heatmap
# ggplot(fst_matrix_melted,  aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile(colour="grey95") +
#   geom_text(aes(label = ifelse(is.na(value), "", round(value, 2))), color = "white", size = 4) +
#   scale_fill_gradient(low = "grey90", high = "darkblue", na.value = "transparent", limits = c(0, 1)) +
#   labs(title = "Heatmap of FST Values", x = "Populations", y = "Populations") +
#   guides(fill = guide_colourbar(title = "FST")) +
#   theme_minimal()+
#   theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     panel.border = element_blank(),       # Remove plot border
#     axis.text.x= element_text(angle = -45, vjust = 0.5, hjust = 0)  # Rotate x-axis labels
#   )



# Plot anonymous heatmap
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
# ggsave(plot= fst_matrix_pl,
#        width = 6.5,
#        height = 6,
#        units = c( "in"),
#        dpi = 600,
#        filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/figures/fst_matrix.png")
# shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/figures/fst_matrix.png")

## FST between host types
vicec_f <- vicec %>%
  filter(AnonSiteCode %in% rownames(anon_fst_matrix)) %>%
  mutate(Host = ifelse(AnonSiteCode %in% c("NE1", "WS5"), 
                       "Suspected salmon", Host)) 
# Exclude NA values and update Host to 'Suspected salmon' where necessary
vicec_f <- vicec %>% 
  filter(AnonSiteCode %in% rownames(anon_fst_matrix)) %>% 
  mutate(Host = ifelse(AnonSiteCode %in% c("NE1", "WS5"), "Suspected salmon", Host))

# Now filter for salmon populations (excluding NA values from the AnonSiteCode column)
salmon <- anon_fst_matrix_melted %>% 
  filter(Var1 %in% vicec_f$AnonSiteCode[!is.na(vicec_f$AnonSiteCode) & vicec_f$Host %in% c("Salmon", "Suspected salmon")] & 
           Var2 %in% vicec_f$AnonSiteCode[!is.na(vicec_f$AnonSiteCode) & vicec_f$Host %in% c("Salmon", "Suspected salmon")]) %>% 
  na.omit()
unique(vicec_f$AnonSiteCode[!is.na(vicec_f$AnonSiteCode) & vicec_f$Host %in% c("Salmon", "Suspected salmon")])
union(salmon$Var1, salmon$Var2)
unique(vicec_f$AnonSiteCode[vicec$Host %in% c("Salmon","Suspected salmon" )])
trout <- anon_fst_matrix_melted %>%
  filter(Var1 %in% vicec_f$AnonSiteCode[!is.na(vicec_f$AnonSiteCode) & vicec_f$Host %in% c("Trout", "Only trout available")] & 
           Var2 %in% vicec_f$AnonSiteCode[!is.na(vicec_f$AnonSiteCode) & vicec_f$Host %in% c("Trout", "Only trout available")])%>% 
  na.omit()
union(trout$Var1, trout$Var2)
unique(vicec_f$AnonSiteCode[vicec$Host %in% c("Trout", "Only trout available","Unknown")])

unknown <- anon_fst_matrix_melted %>%
  filter(Var1 %in% vicec_f$AnonSiteCode[!is.na(vicec_f$AnonSiteCode) & vicec_f$Host %in% c("Unknown")] & 
           Var2 %in% vicec_f$AnonSiteCode[!is.na(vicec_f$AnonSiteCode) & vicec_f$Host %in% c("Unknown")])%>% 
  na.omit()
union(unknown$Var1, unknown$Var2)

# Filter for salmon-trout pairs with correct handling of NA values and Host conditions
salmon_trout_pairs <- anon_fst_matrix_melted %>%
  filter(
    (Var1 %in% vicec_f$AnonSiteCode[!is.na(vicec_f$AnonSiteCode) & vicec_f$Host %in% c("Salmon", "Suspected salmon")] & 
       Var2 %in% vicec_f$AnonSiteCode[!is.na(vicec_f$AnonSiteCode) & vicec_f$Host %in% c("Trout", "Only trout available")]) |
      (Var1 %in% vicec_f$AnonSiteCode[!is.na(vicec_f$AnonSiteCode) & vicec_f$Host %in% c("Trout", "Only trout available")] & 
         Var2 %in% vicec_f$AnonSiteCode[!is.na(vicec_f$AnonSiteCode) & vicec_f$Host %in% c("Salmon", "Suspected salmon")])
  )%>% 
  na.omit()
range(salmon$value)
range(trout$value)
range(salmon_trout_pairs$value)

nrow(na.omit(anon_fst_matrix_melted))
nrow(bind_rows(salmon, trout, salmon_trout_pairs))

all(na.omit(anon_fst_matrix_melted) %in% bind_rows(salmon, trout, salmon_trout_pairs))
salmon_trout_pairs <- anon_fst_matrix_melted %>%
  filter(
    (Var1 %in% vicec$AnonSiteCode[vicec$Host == "Salmon"] & 
       Var2 %in% vicec$AnonSiteCode[vicec$Host %in% c("Trout", "Only trout available")]) | 
      (Var1 %in% vicec$AnonSiteCode[vicec$Host %in% c("Trout", "Only trout available")] & 
         Var2 %in% vicec$AnonSiteCode[vicec$Host == "Salmon"])
  )
head(anon_fst_matrix_melted)
### =========  Plot FST versus distance between populations =================

geo_dist<-readRDS( "../confidential/distance_matrix.rds")
geo_dist
geo_dist_long <- melt(geo_dist)
# Sort geo_dist_long by Var1 and Var2
geo_dist_long <- geo_dist_long %>%
  arrange(Var1, Var2)

# Sort fst_matrix_melted by Var1 and Var2
fst_matrix_melted <- fst_matrix_melted %>%
  arrange(Var1, Var2)
head(geo_dist_long)
head(fst_matrix_melted)


merged_data <- merge(geo_dist_long, fst_matrix_melted, by = c("Var1", "Var2"))%>%
  rename(Distance = value.x,
    FST = value.y)%>% na.omit()
head(merged_data)

# Function to linearise Fst values
linearise_fst <- function(fst) {
  return(fst / (1 - fst))
}

# Apply the function to the Fst column
merged_data$Linearised_Fst <- sapply(merged_data$FST, linearise_fst)


library(ggforce)
ggplot(merged_data, aes(x = Distance, y = FST)) +
  geom_point(size=3, alpha=0.5, colour="darkblue") +  # Scatter plot of points
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear regression line
  labs(title = bquote("Distance vs " ~ F[ST]), x = "Distance (km)", y = bquote(F[ST])) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits+
  theme_classic()
  # theme_classic()+
  # ggforce::geom_mark_hull(
  #   aes(label = "OH1 vs WS1",
  #       filter = FST > 0.6),
  #   # description = "Commuters had to deal with severe disruptions in public transport on July 9 and August 6",
  #   color = "black",
  #   # label.family = "Asap SemiCondensed",
  #   label.fontsize = c(14)
  # )+
  # ggforce::geom_mark_hull(
  #   aes(label = paste0(Var1, " vs ", Var2),
  #       filter = FST > 0.5,
  #   # description = "Commuters had to deal with severe disruptions in public transport on July 9 and August 6",
  #   color = FST),
  #   # label.family = "Asap SemiCondensed",
  #   label.fontsize = c(14)
  # )

head(geo_dist)
head(fst_matrix)

# Perform linear regression
lm_model <- lm(FST ~ Distance, data = merged_data)

# Calculate residuals
merged_data$residuals <- residuals(lm_model)

# Define a threshold for significant divergence (adjust as needed)
residual_threshold <- 0.2  # Example threshold for significant divergence

# Plot Distance vs FST with labels for significantly divergent points
# install.packages("ggrepel")

ggplot(merged_data, aes(x = Distance, y = FST)) +
  geom_point() +  # Scatter plot of points
  geom_text(data = subset(merged_data, abs(residuals) > residual_threshold),
            aes(label = paste0(Var1, " vs ", Var2)),
            hjust = 1, vjust = -1, size = 3) +  # Label points with high residuals
  labs(title = bquote("Distance vs " ~ F[ST]), x = "Distance (km)", y = bquote(F[ST])) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_minimal()
# Plot Distance vs FST with labels for significantly divergent points
ggplot(merged_data, aes(x = Distance, y = Linearised_Fst)) +
  geom_point() +  # Scatter plot of points
  # ggrepel::geom_text_repel(data = subset(merged_data, abs(residuals) > residual_threshold),aes(label = paste0(Var1, " vs ", Var2)), size = 3, box.padding = 0.5, max.overlaps = Inf) +  # Label points with high residuals, repel labels
  labs(#title = bquote("Distance vs " ~ F[ST]),
       x = "Distance (km)", y = expression(paste(F[ST], "/ (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_classic()



#### Fish distance ====================================================================

head(merged_data)
fishswim<-read.csv("../confidential/point_distance_matrix_as_fish_swims")
head(fishswim)

# Assuming merged_data and fishswim are data frames
merged_data_fish <- merged_data %>%
  left_join(fishswim, by = c("Var1" = "start", "Var2" = "end"))

merged_data_fish <- merged_data_fish %>%
  mutate(fish_dist = ifelse(dist == 0, Distance, dist)) %>%  #some pops in same river so I take the eucliden distance here (grig vs briagh and lax vs lax_unballach and polly vs pollystaccburn)
  dplyr::select(-dist)
head(merged_data_fish)

## Plot fish and bird
pl_fish<-ggplot(merged_data_fish, aes(x = fish_dist, y = Linearised_Fst)) +
  geom_point(size=3, alpha=0.6, colour="darkblue") +  # Scatter plot of points
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear regression line
  labs(#title = bquote("Fish Distance vs " ~ F[ST]),
       x = "Least-Cost-Path Distance (km)", y = expression(paste(F[ST], "/ (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 2)) +  # Set y-axis limits+
  theme_classic()
pl_fish
pl_bird<-ggplot(merged_data_fish, aes(x = Distance, y = Linearised_Fst)) +
  geom_point(size=3, alpha=0.6, colour="darkred") +  # Scatter plot of points
   geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear regression line
  labs(#title = bquote("Bird Distance vs " ~ F[ST]), 
       x = "Euclidean Distance (km)", y = expression(paste(F[ST], "/ (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 2)) +  # Set y-axis limits+
  theme_classic()
pl_bird
library(patchwork)
combi_pl <- pl_fish | pl_bird 
  combi_pl_s<-combi_pl+plot_annotation(
    tag_levels = 'a',  # Labels start with 'a', 'b', etc.
    tag_suffix = ')'   # Suffix for labels, e.g., a)
  )
combi_pl_s


head(merged_data_fish)
head(vicec)
ggplot(
  merged_data_fish %>%
    filter(Var1 %in% vicec$Population[vicec$Host == "Salmon"] & 
             Var2 %in% vicec$Population[vicec$Host == "Salmon"]), 
  aes(x = fish_dist, y = Linearised_Fst)
) +
  geom_point(size=3, alpha=0.5, colour="darkblue") +  # Scatter plot of points
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear regression line
  labs(title = bquote("Distance vs " ~ F[ST]),x = "Least-Cost-Path Distance (km)", y = expression(paste(F[ST], "/ (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_classic()
ggplot(
  merged_data_fish %>%
    filter(Var1 %in% vicec$Population[vicec$Host %in% c("Trout", "Only trout available")] & 
             Var2 %in% vicec$Population[vicec$Host %in% c("Trout", "Only trout available")]), 
  aes(x = fish_dist, y = Linearised_Fst)
) +
  geom_point(size=3, alpha=0.5, colour="darkblue") +  # Scatter plot of points
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear regression line
  labs(title = bquote("Distance vs " ~ F[ST]),x = "Least-Cost-Path Distance (km)", y = expression(paste(F[ST], "/ (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_classic()
unique(merged_data_fish %>%
         filter(Var1 %in% vicec$Population[vicec$Host %in% c("Trout", "Only trout available")] & 
                  Var2 %in% vicec$Population[vicec$Host %in% c("Trout", "Only trout available")]))
# Perform linear regression
lm_model_fish <- lm(Linearised_Fst ~ Distance, data = merged_data_fish)

#Grab stats from model for write up 
summary(lm_model_fish)
summary(lm_model)
plot(lm_model_fish$fitted.values, sqrt(abs(lm_model_fish$residuals)),
     xlab = "Fitted Values", ylab = "Square Root of Absolute Residuals",
     main = "Scale-Location")
qqnorm(lm_model_fish$residuals)

cor(merged_data$FST, merged_data$Distance, method = "pearson")
cor(merged_data_fish$FST, merged_data_fish$fish_dist, method = "pearson")
# Calculate residuals
merged_data_fish$residuals <- residuals(lm_model_fish)

# Define a threshold for significant divergence (adjust as needed)
residual_threshold <- 0.2  # Example threshold for significant divergence

ggplot(merged_data_fish, aes(x = fish_dist, y = FST)) +
  geom_point() +  # Scatter plot of points
  geom_text(data = subset(merged_data_fish, abs(residuals) > residual_threshold),
            aes(label = paste0(Var1, " vs ", Var2)),
            hjust = 1, vjust = -1, size = 3) +  # Label points with high residuals
  labs(title = bquote("Distance vs " ~ F[ST]), x = "Distance (km)", y = bquote(F[ST])) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_minimal()
ggplot(merged_data_fish, aes(x = fish_dist, y = FST)) +
  geom_point() +  # Scatter plot of points
  geom_text(data = subset(merged_data_fish, fish_dist == 0),
            aes(label = paste0(Var1, " vs ", Var2)),
            hjust = -1, vjust = -1, size = 3) +  # Label points with high residuals
  labs(title = bquote("Distance vs " ~ F[ST]), x = "Distance (km)", y = bquote(F[ST])) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_minimal()
ggplot(merged_data_fish, aes(x = fish_dist, y = FST)) +
  geom_point() +  # Scatter plot of points
  geom_text(data = merged_data_fish,
            aes(label = paste0(Var1, " vs ", Var2)),
            hjust = 1, vjust = -1, size = 2.5) +  # Label points with high residuals
  labs(title = bquote("Distance vs " ~ F[ST]), x = "Distance (km)", y = bquote(F[ST])) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_minimal()
# Plot Distance vs FST with labels for significantly divergent points
ggplot(merged_data_fish, aes(x = fish_dist, y = Linearised_Fst)) +
  geom_point() +  # Scatter plot of points
  # ggrepel::geom_text_repel(data = subset(merged_data, abs(residuals) > residual_threshold),aes(label = paste0(Var1, " vs ", Var2)), size = 3, box.padding = 0.5, max.overlaps = Inf) +  # Label points with high residuals, repel labels
  labs(#title = bquote("Distance vs " ~ F[ST]),
    x = "As-the-fish-swims Distance (km)", y = expression(paste(F[ST], "/ (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_classic()


ggplot(merged_data_fish, aes(x = fish_dist)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of Distance Measurements", x = "Distance", y = "Density") +
  theme_minimal()
save.image(file = "../../data/too_large_files/Gen03_FST.RData")


### Linearised FST distances using Mantel test to account for the non-inpep pairwise testing =========================================
library(vegan)

# Duplicate for both directions of matrix!!
linearised_full <- merged_data_fish %>%
  # select(Var1, Var2, Linearised_Fst, Distance, fish_dist) %>%
  bind_rows(merged_data_fish %>% rename(Var1 = Var2, Var2 = Var1))
head(linearised_full)
# Make matrix from the square format table
fst_matrix_lin <- acast(linearised_full, Var1 ~ Var2, value.var = "Linearised_Fst")
fst_matrix_lin
# fst_matrix_lin <- acast(linearised_full, Var1 ~ Var2, value.var = "FST") #using non-linear fst gives the same result
# fst_matrix_lin
all.equal(fst_matrix_lin, t(fst_matrix_lin), check.attributes = FALSE)
# diag(fst_matrix_lin) <- 0 #mantel test doesnt like the NA vals?

euclidist_matrix_lin <- acast(linearised_full, Var1 ~ Var2, value.var = "Distance")
euclidist_matrix_lin
all.equal(euclidist_matrix_lin, t(euclidist_matrix_lin), check.attributes = FALSE)
# diag(euclidist_matrix_lin) <- 0 #mantel test doesnt like the NA vals?

fishdist_matrix_lin <- acast(linearised_full, Var1 ~ Var2, value.var = "fish_dist")
fishdist_matrix_lin
all.equal(fishdist_matrix_lin, t(fishdist_matrix_lin), check.attributes = FALSE)
# diag(fishdist_matrix_lin) <- 0 #mantel test doesnt like the NA vals?

# Make the dist() format


# Now run Mantel tests
mantel_geo <- mantel(fst_matrix_lin, euclidist_matrix_lin, method = "spearman", permutations = 9999) #Mantel r = 0.2054, p = 0.0557
mantel_geo
mantel_fish <- mantel(fst_matrix_lin, fishdist_matrix_lin, method = "spearman", permutations = 9999) #Mantel r = 0.0004, p = 0.487
mantel_fish


#---------------------------------------------------------------------------------------------------------------------------------------

## Comparing Fst in slamon and trout
# Count number of populations by host type
vicec %>%
  group_by(Host) %>%
  summarise(n_populations = n_distinct(Population))

#mantel test here old cr test below
# Filter for salmon population from the full linearised data
salmon <- linearised_full %>%
  filter(Var1 %in% vicec$Population[vicec$Host %in% c("Salmon")] & 
           Var2 %in% vicec$Population[vicec$Host %in% c("Salmon")])

# Create Fst distance matrix for salmon
fst_matrix_salmon <- acast(salmon, Var1 ~ Var2, value.var = "Linearised_Fst")
euclidist_matrix_salmon <- acast(salmon, Var1 ~ Var2, value.var = "Distance")
fishdist_matrix_salmon <- acast(salmon, Var1 ~ Var2, value.var = "fish_dist")

# Check for symmetry of matrices (important for Mantel test)
all.equal(fst_matrix_salmon, t(fst_matrix_salmon), check.attributes = FALSE)
all.equal(euclidist_matrix_salmon, t(euclidist_matrix_salmon), check.attributes = FALSE)
all.equal(fishdist_matrix_salmon, t(fishdist_matrix_salmon), check.attributes = FALSE)

# Filter for trout population from the full linearised data
trout <- linearised_full %>%
  filter(Var1 %in% vicec$Population[vicec$Host %in% c("Trout", "Only trout available")] & 
           Var2 %in% vicec$Population[vicec$Host %in% c("Trout", "Only trout available")])

# Create Fst distance matrix for trout
fst_matrix_trout <- acast(trout, Var1 ~ Var2, value.var = "Linearised_Fst")
euclidist_matrix_trout <- acast(trout, Var1 ~ Var2, value.var = "Distance")
fishdist_matrix_trout <- acast(trout, Var1 ~ Var2, value.var = "fish_dist")

# Check for symmetry of matrices (important for Mantel test)
all.equal(fst_matrix_trout, t(fst_matrix_trout), check.attributes = FALSE)
all.equal(euclidist_matrix_trout, t(euclidist_matrix_trout), check.attributes = FALSE)
all.equal(fishdist_matrix_trout, t(fishdist_matrix_trout), check.attributes = FALSE)


# Run Mantel test between Fst and Euclidean distance (Spearman method)
mantel_salmon_geo <- mantel(fst_matrix_salmon, euclidist_matrix_salmon, method = "spearman", permutations = 9999) #No evidence of IBD (r = –0.5, p = 0.833).
print(mantel_salmon_geo)

# Run Mantel test between Fst and Fish distance
mantel_salmon_fish <- mantel(fst_matrix_salmon, fishdist_matrix_salmon, method = "spearman", permutations = 9999) #Weak positive correlation but not significant (r = 0.5, p = 0.5).
print(mantel_salmon_fish)


# Run Mantel test between Fst and Euclidean distance (Spearman method)
mantel_trout_geo <- mantel(fst_matrix_trout, euclidist_matrix_trout, method = "spearman", permutations = 9999) #A significant positive correlation was found (r = 0.3895, p = 0.007)
print(mantel_trout_geo)

# Run Mantel test between Fst and Fish distance 
mantel_trout_fish <- mantel(fst_matrix_trout, fishdist_matrix_trout, method = "spearman", permutations = 9999) #A weaker correlation approaching significance (r = 0.2776, p = 0.0587)
print(mantel_trout_fish)


# Salmon plots
pl_salmon_fish <- ggplot(salmon, aes(x = fish_dist, y = Linearised_Fst)) +
  geom_point(size = 3, alpha = 0.6, colour = "darkblue") +
  geom_smooth(method = "lm", se = FALSE, colour = "blue") +
  labs(x = "Least-Cost-Path Distance (km)", 
       y = expression(paste(F[ST], " / (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 2)) +
  theme_classic()
pl_salmon_fish
pl_salmon_geo <- ggplot(salmon, aes(x = Distance, y = Linearised_Fst)) +
  geom_point(size = 3, alpha = 0.6, colour = "darkred") +
  geom_smooth(method = "lm", se = FALSE, colour = "blue") +
  labs(x = "Euclidean Distance (km)", 
       y = expression(paste(F[ST], " / (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 2)) +
  theme_classic()

# Trout plots
pl_trout_fish <- ggplot(trout, aes(x = fish_dist, y = Linearised_Fst)) +
  geom_point(size = 3, alpha = 0.6, colour = "darkblue") +
  geom_smooth(method = "lm", se = FALSE, colour = "blue") +
  labs(x = "Least-Cost-Path Distance (km)", 
       y = expression(paste(F[ST], " / (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 2)) +
  theme_classic()

pl_trout_geo <- ggplot(trout, aes(x = Distance, y = Linearised_Fst)) +
  geom_point(size = 3, alpha = 0.6, colour = "darkred") +
  geom_smooth(method = "lm", se = FALSE, colour = "blue") +
  labs(x = "Euclidean Distance (km)", 
       y = expression(paste(F[ST], " / (1 - ", F[ST], ")"))) +
  scale_y_continuous(limits = c(0, 2)) +
  theme_classic()

# Combine plots
# library(patchwork)
final_plot <- (
  (pl_salmon_fish | pl_salmon_geo) + 
    plot_annotation(title = "Isolation by Distance in Salmon Populations")
) /
  (
    (pl_trout_fish | pl_trout_geo) + 
      plot_annotation(title = "Isolation by Distance in Trout Populations")
  ) +
  plot_annotation(tag_levels = 'a', tag_suffix = ')')

# View combined plot
final_plot

# Save it
ggsave("../../output/figures/combined_ibd_plot_abis salmon and bc is trout pops_trout mantel is sig but probs because so unequal distrib of dist.png", final_plot, width = 8, height = 8, dpi = 600)
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/combined_ibd_plot_abis salmon and bc is trout pops_trout mantel is sig but probs because so unequal distrib of dist.png")


#============= End of neat code  =========


##======================

merged_data_fish
salmon<-merged_data_fish %>%
  filter(Var1 %in% vicec$Population[vicec$Host %in% c("Salmon")] & 
           Var2 %in% vicec$Population[vicec$Host %in% c("Salmon")])
union(salmon$Var1, salmon$Var2)
unique(vicec$Population[vicec$Host %in% c("Salmon")])
trout<-merged_data_fish %>%
  filter(Var1 %in% vicec$Population[vicec$Host %in% c("Trout", "Only trout available")] & 
           Var2 %in% vicec$Population[vicec$Host %in% c("Trout", "Only trout available")])
union(trout$Var1, trout$Var2)
unique(vicec$Population[vicec$Host %in% c("Trout", "Only trout available")])

range(salmon$FST)
range(salmon$dist)
range(trout$FST)
range(trout$dist)

lm_model_fish <- lm(Linearised_Fst ~ fish_dist, data = salmon)

#Grab stats from model for write up 
summary(lm_model_fish)
plot(lm_model_fish$fitted.values, sqrt(abs(lm_model_fish$residuals)),
     xlab = "Fitted Values", ylab = "Square Root of Absolute Residuals",
     main = "Scale-Location")
qqnorm(lm_model_fish$residuals)

cor(salmon$Linearised_Fst, salmon$Distance, method = "pearson")
cor(merged_data_fish$Linearised_Fst, merged_data_fish$fish_dist, method = "pearson")
# Calculate residuals
salmon$residuals <- residuals(lm_model_fish)

ggplot(salmon, aes(x = fish_dist, y = Linearised_Fst)) +
  geom_point() +  # Scatter plot of points
  labs(title = bquote("Distance vs " ~ F[ST]), x = "Distance (km)", y = bquote(F[ST])) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_minimal()

cor(salmon$FST, salmon$fish_dist, method = "spearman")
######Trout
lm_model_fish <- lm(Linearised_Fst ~ fish_dist, data = trout)

#Grab stats from model for write up 
summary(lm_model_fish)
plot(lm_model_fish$fitted.values, sqrt(abs(lm_model_fish$residuals)),
     xlab = "Fitted Values", ylab = "Square Root of Absolute Residuals",
     main = "Scale-Location")
qqnorm(lm_model_fish$residuals)

cor(trout$Linearised_Fst, trout$Distance, method = "pearson")
cor(merged_data_fish$Linearised_Fst, merged_data_fish$fish_dist, method = "pearson")
# Calculate residuals
trout$residuals <- residuals(lm_model_fish)

ggplot(trout, aes(x = fish_dist, y = Linearised_Fst)) +
  geom_point() +  # Scatter plot of points
  labs(title = bquote("Distance vs " ~ F[ST]), x = "Distance (km)", y = bquote(F[ST])) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_minimal()

cor(trout$FST, trout$fish_dist, method = "spearman")

