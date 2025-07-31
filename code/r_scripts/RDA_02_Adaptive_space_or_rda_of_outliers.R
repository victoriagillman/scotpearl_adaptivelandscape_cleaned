## RDA based on outlier loci
#Leading on from RDA code 
library(ggforce)
library(vegan)
library(dplyr)

load("data/confidential/RDA_host_condition_updated.RData")

gen.imp<-gen.imp
outliers<-Out_qv
nrow(outliers)
pred_prune<-pred_prune
labels<-labels

gen.imp_subset <- gen.imp[, colnames(gen.imp) %in% rownames(Out_qv)]

out.rda<-rda(gen.imp_subset ~ FW_dem_range + wildness + WC_alt + FW_hydro_wavg_03 +FW_hydro_wavg_05+ FW_hydro_wavg_17+ catchment_area_km + Condition(Host), data = pred_prune, scale=T)
out.rda
RsquareAdj(out.rda)

out.signif.full <- anova.cca(out.rda, parallel=getOption("mc.cores"), permutations = how(nperm=999)) #p<0.001
out.signif.axis <- anova.cca(out.rda, by="axis", parallel=getOption("mc.cores"), permutations = how(nperm=999)) #all <0.001
out.signif.term <- anova.cca(out.rda, by="term", parallel=getOption("mc.cores")) #all <0.001

out.signif.full 
out.signif.axis 
out.signif.term
var_table<-data.frame(Predictor= rownames(out.signif.term), Variance=out.signif.term$Variance) %>%
  filter(Predictor != "Residual") %>% 
  mutate(Proportion = Variance / sum(Variance)) %>% 
  left_join(labels, by = c("Predictor" = "layer_code")) %>% 
  arrange(desc(Proportion)) 
var_table$name <- factor(var_table$name, levels = var_table$name[order(var_table$Proportion, decreasing = TRUE)])
var_table
var_table <- var_table %>%
  mutate(
    name = case_when(
      name == "BIO3 = Isothermality (BIO2/BIO7) (x100)" ~ "BIO3",   
      name == "BIO5 = Max Temperature of Warmest Month" ~ "BIO5",   
      name == "BIO17 = Precipitation of Driest Quarter" ~ "BIO17", 
      name == "Catchment area (km)" ~ "Catchment area", 
      name == "Wildness from roads" ~ "Wildness", 
      TRUE ~ as.character(name)         # Keep other names unchanged
    )
  ) %>% 
  arrange(desc(Proportion)) 

var_table <- data.frame(Predictor = rownames(out.signif.term), Variance = out.signif.term$Variance) %>%
  filter(Predictor != "Residual") %>%
  mutate(Proportion = Variance / sum(Variance)) %>%
  left_join(labels, by = c("Predictor" = "layer_code")) %>%
  mutate(
    name = case_when(
      name == "BIO3 = Isothermality (BIO2/BIO7) (x100)" ~ "BIO3",
      name == "BIO5 = Max Temperature of Warmest Month" ~ "BIO5",
      name == "BIO17 = Precipitation of Driest Quarter" ~ "BIO17",
      name == "Catchment area" ~ "Catchment area (km)",
      name == "Wildness from roads" ~ "Wildness",
      TRUE ~ name
    )
  ) %>%
  arrange(desc(Proportion)) %>%
  mutate(name = factor(name, levels = name))

var_imp_plot<-ggplot(var_table, aes(x = name, y = Variance, fill = name)) +
  geom_bar(stat = "identity") +
  labs(x = "Predictor",
       y = "Variance") +
  # scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")
var_imp_plot
# ggsave(var_imp_plot, width = 5,
#        height = 3.5,
#        units = c( "in"),
#        dpi = 600,
#        filename = "./output/figures/rda_on_outqv_varimportance.png")
# shell.exec(paste0(dir,"/output/figures/rda_on_outqv_varimportance.png"))

out.rda[["call"]]
summary(eigenvals(out.rda, model = "constrained"))

# Extract the eigenvalues and calculate the proportion of variance explained
eigenvalues <- eigenvals(out.rda, model = "constrained")
eigenvalues
variance_explained <- eigenvalues / sum(eigenvalues)
variance_explained
varience_for_writeup<-round(variance_explained * 100, 3)  
varience_for_writeup
# Extract the proportions for RDA1 and RDA2
rda1_prop <- round(variance_explained[1] * 100, 1)  
rda2_prop <- round(variance_explained[2] * 100, 1) 
out.rda$terminfo$xlev

num_axis<-6
rda_stats <- extract_rda_stats(rda_obj = out.rda, gen_imp = gen.imp_subset, num_axes = num_axis, scaling = "none")

# Extract the scores for environmental variables (the arrows and lines)
env_scores_df <- as.data.frame(scores(out.rda, display = "bp", choices=c(1,2)))
env_scores_df$variable <- rownames(env_scores_df)
unique_predictors <- unique(env_scores_df$variable)
columns_without_layer_code <- setdiff(unique_predictors, labels$layer_code)
# Add new labels to labels df based on ones missing
labels<-labels%>% dplyr::select(layer_code,name)
labels <- rbind(labels, data.frame(layer_code = columns_without_layer_code, name = columns_without_layer_code))
env_scores_df <- left_join(env_scores_df, labels, by = c("variable" = "layer_code"))
env_scores_df
unique(env_scores_df$name)

env_scores_df %>%
  mutate(RDA1 = abs(RDA1),
         RDA2 = abs(RDA2)) %>%
  arrange(desc(RDA2)) %>%
  dplyr::select(variable, RDA1, RDA2)

site_scores <- as.data.frame(scores(out.rda, display = "sites", choices=1:4, scaling = 3))
site_scores$Sample_ID <- rownames(site_scores)
rda_env<-left_join(env,site_scores, by = "Sample_ID") %>% filter(!is.na(RDA1))

# Plotting with ggplot2
snp_gg_rda<-ggplot() +
  # geom_point(color = "grey90", size = 1.4) +  # Plot all SNPs with normal scaling in grey
  geom_point(data = rda_stats, aes(x = RDA1*3, y = RDA2*3), colour = "deeppink", size = 3, alpha = 0.2)+ # SNP loadings
  labs(x =paste0("RDA1 (",rda1_prop, "%)"), y =paste0("RDA2 (", rda2_prop, "%)")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Horizontal line at y=0
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Vertical line at x=0
  geom_segment(data = env_scores_df, aes(x = 0, y = 0, xend = RDA1 * 1.5, yend = RDA2 * 1.5),                colour = "black", linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +  # Arrows from origin to scaled env scores
  geom_text_repel(data = env_scores_df, aes(x = 1.5 * RDA1, y = 1.5 * RDA2, label = str_wrap(name, width = 35)), 
                  size=3.5,
                  min.segment.length = Inf,
                  position = position_nudge_center(0.2, 0.1, 0, 0),
                  max.overlaps = Inf, # Adjust as necessary for more flexibility
                  box.padding = 0.1, # Adjust padding around labels
                  point.padding = 0.1) +  # Adjust padding around points
  theme_classic()+
  theme(legend.position = "right",  # Position legend at the bottom
        legend.key.size = unit(0.5, "cm"))+  # Adjust legend key size geom_point(data = rda_env, aes(color = Population), size = 3)+  # Overlay candidate SNPs with special scaling
  labs(colour="Anon Population")#  coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) 
snp_gg_rda
# ggsave(snp_gg_rda, width = 8.5,height = 6,units = c( "in"),dpi = 600, filename = "./output/figures/rda_just_outliers.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/rda_just_outliers.png")
head(rda_stats)

summary(rda_stats$RDA1)

pop_ggrda<-ggplot() + #clustered populations suggest they are adapted to similar environmental conditions
  geom_point(data = rda_env, aes(x = RDA1, y = RDA2, colour = AnonSiteCode), size = 3) + 
  ggforce::geom_mark_ellipse(data = rda_env %>% filter(AnonSiteCode %in% c("WS1", "OH1", "ER1", "WS3", "WS2")), 
                    aes(x = RDA1, y = RDA2, group = AnonSiteCode, label = AnonSiteCode), 
                    fill = NA, 
                    expand = unit(2, "mm"), 
                    label.fontsize = 7,
                    label.buffer = unit(-1, 'mm'),
                    label.fontface = "italic",
                    label.fill = NA,
                    con.type = "none") +
  labs(x = paste0("RDA1 (", rda1_prop, "%)"), y = paste0("RDA2 (", rda2_prop, "%)"), colour = "Population"
       ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Horizontal line at y=0
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Vertical line at x=0
  geom_segment(data = env_scores_df, aes(x = 0, y = 0, xend = RDA1 * 10, yend = RDA2 * 10), colour = "black", 
               linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +  # Arrows from origin to scaled env scores
  geom_text_repel(data = env_scores_df, aes(x = 10.5 * RDA1, y = 10.5 * RDA2, label = str_wrap(name, width = 20)), 
                  size = 3.5,
                  min.segment.length = Inf,
                  position = position_nudge_center(0.2, 0.1, 0, 0),
                  max.overlaps = Inf,
                  box.padding = 0.1, 
                  point.padding = 0.1) + 
  theme_classic() +
  theme(#legend.position = "bottom",  
        legend.position = c(0.505, 0), 
        legend.justification = c(0, 0), 
        legend.box.just = "left",  
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_blank(),
        legend.title=element_text(hjust=0.5), 
        legend.direction = "horizontal"  # Set horizontal direction for the legend
  )+
  guides(colour = guide_legend(title.position = "top",nrow = 3,  # Wrap into 3 rows
         byrow = TRUE))

pop_ggrda
# ggsave(pop_ggrda, width = 8.5,height = 6,units = c( "in"),dpi = 600, filename = "./output/figures/rda_just_outliers_pop_in_adaptive_space.png")
# shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/rda_just_outliers_pop_in_adaptive_space.png")


library(patchwork)
comp_plot <- snp_gg_rda / pop_ggrda + 
  plot_annotation(tag_levels = 'a', tag_suffix = ")")
comp_plot
ggsave(comp_plot, width = 8.5,height = 12,units = c( "in"),dpi = 600, filename = "./output/figures/rda_just_outliers_composite_in_adaptive_space.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/rda_just_outliers_composite_in_adaptive_space.png")
# ggplot() +
#   geom_point(data = rda_stats, aes(x = RDA1*8, y = RDA2*8), colour = "grey", size = 3, alpha = 0.2)+ # SNP loadings
#   geom_point(data = rda_env, aes(x = RDA1, y = RDA2, colour = AnonSiteCode), size = 3) + 
#   geom_mark_ellipse(data = rda_env %>% filter(AnonSiteCode %in% c("WS1", "OH1", "ER1", "WS3", "WS2")), 
#                     aes(x = RDA1, y = RDA2, group = AnonSiteCode, label = AnonSiteCode), 
#                     fill = NA, 
#                     expand = unit(2, "mm"), 
#                     label.fontsize = 7,
#                     label.buffer = unit(-1, 'mm'),
#                     label.fontface = "italic",
#                     label.fill = NA,
#                     con.type = "none") +
#   labs(x =paste0("RDA1 (",rda1_prop, "%)"), y =paste0("RDA2 (", rda2_prop, "%)"), title = "Outlier RDA showing all SNPs in adaptive space") +
#   geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Horizontal line at y=0
#   geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Vertical line at x=0
#   geom_segment(data = env_scores_df, aes(x = 0, y = 0, xend = RDA1 * 1.5, yend = RDA2 * 1.5),                colour = "black", linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +  # Arrows from origin to scaled env scores
#   geom_text_repel(data = env_scores_df, aes(x = 1.5 * RDA1, y = 1.5 * RDA2, label = str_wrap(name, width = 35)), 
#                   size=3.5,
#                   min.segment.length = Inf,
#                   position = position_nudge_center(0.2, 0.1, 0, 0),
#                   max.overlaps = Inf, # Adjust as necessary for more flexibility
#                   box.padding = 0.1, # Adjust padding around labels
#                   point.padding = 0.1) +  # Adjust padding around points
#   theme_classic()+
#   theme(legend.position = "right",  # Position legend at the bottom
#         legend.key.size = unit(0.5, "cm"))+  # Adjust legend key size geom_point(data = rda_env, aes(color = Population), size = 3)+  # Overlay candidate SNPs with special scaling
#   labs(colour="Anon Population")#  coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) 




library(tibble)
env_scores_df %>%
  as.data.frame() %>%
  arrange(desc(abs(RDA1))) %>%
  dplyr::select(name, RDA1)
env_scores_df %>%
  as.data.frame() %>%
  arrange(desc(abs(RDA2))) %>%
  select(name, RDA2)
