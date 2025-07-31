#Victoria Gillman 25th June 2024. 
# Script for redundancy analysis for genome-environment association of swabbed freshwater pearl mussels in Scotland
# Using VP filtered data
# Needs trimmed env data and plink conversion to .raw from Gen_01_PlinkFiltering

# This code is based on the tutorial https://popgen.nescent.org/2018-03-27_RDA_GEA.html 
#
setwd("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned")

#### Installing and loading necessary packages ####
# install.packages(c("psych","vegan"), dependencies=TRUE)
# Load packages
# -------------
library(dplyr)      # For data manipulation
library(tidyr)
# Analysis
library(adegenet)   # Used for genetic data preparation
library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RDA
library(qvalue)
library(robust)

# Spatial
library(terra)      # For working with spatial data
library(sdmpredictors)  # For accessing environmental predictor datasets
library(raster)     # For working with raster data
library(fuzzySim)   # For working with species distribution modelling data


# Data visualisation
library(ggplot2) 
library(patchwork) #for loadings histograms
library(ggrepel)
library(ggpp)
library(stringr)

# library(poppr)
# -------------

#### Preparing genetic data ####
genRAW <- read.PLINK("data/plink_VP/pearlmusselfiles_raw.raw",
                     map.file = "data/plink_VP/pearlmusselfiles.map" )
length(genRAW@position)
summary(genRAW)
map <- read.table("data/plink_VP/pearlmusselfiles.map", header = FALSE, stringsAsFactors = FALSE)
map_SNPs <- map$V4
colnames(map) <- c("Chromosome", "SNP_ID", "Genetic_Distance", "Position")
genRAW@position<-map_SNPs
genRAW@position

dim(genRAW) #156 5486
locNames(genRAW)
pop(genRAW)
pop(genRAW)[grep("Conon_Gl", pop(genRAW))]

# locNames(genRAW)=paste("SNP", 1:nLoc(genRAW),Sep="_") #give snps unique ids
sum(is.na(genRAW))  #0
gen<-as.matrix(genRAW)
sum(is.na(gen))  #34632 (so 21647/(167*5173)=0.025=2.5% missing data)
cat("Percentage of missing data:", (((sum(is.na(gen))) / (dim(genRAW)[1] * dim(genRAW)[2])) * 100), "%\n")

gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) #a function across each column (SNP) of the genomic data matrix. Within the function, missing values (is.na(x)) are replaced with the most common genotype (determined by which.max(table(x))) for that SNP.
sum(is.na(gen.imp)) # No NAs
###THERE ARE OTHER IMPUTATION METHODS HIGHLIGHTED BY THE TUTORIAL


#### Read in environmental predictors ####
env<-read.csv("data/confidential/trimmed_env_data.csv")  %>% mutate(Population = ifelse(Population == "Conon_Gharbhrain_wrong_sample_ID", "Conon_Glascarnoch", Population)) %>% mutate(AnonSiteCode = ifelse(AnonSiteCode == "ER2", "ER1", AnonSiteCode))
length(unique(env$AnonSiteCode))
# p_env<-read.csv("data/confidential/SNP_metadata_downsampled_geneticPCs_K4.csv")%>%rename(Sample_ID = ID) %>%dplyr::select(-SiteCode)
# setdiff(colnames(p_env), colnames(env))
# env<-p_env #adding the 4 pca results
env$Sample_ID <- as.character(env$Sample_ID) # Make individual names characters (not factors)

# Reorder so geno and env are in the same order
setdiff(rownames(gen.imp), env$Sample_ID) #id_not_in_env<-
setdiff(env$Sample_ID, rownames(gen.imp)) #id_not_in_gen.imp<-
rownames(gen.imp)[rownames(gen.imp) == "ICY837notICY773"] <- "ICY837"

identical(rownames(gen.imp), env$Sample_ID)  # False
env <- env %>%
  # dplyr::select(-X)%>%
  filter(Sample_ID %in% rownames(gen.imp)) %>%
  arrange(match(Sample_ID, rownames(gen.imp)))
identical(rownames(gen.imp), env$Sample_ID)  # TRUE

## Import PCA val for pop rel RDA
PC_K4_meta<- read.csv( "data/confidential/PCA_values_for_indi_K4.csv") %>% rename(Sample_ID = ID) %>% mutate(Sample_ID = ifelse(Sample_ID == "ICY837notICY773","ICY837", Sample_ID))
PC_env <- env %>%
  left_join(PC_K4_meta %>% dplyr::select(Sample_ID, "PC1_genetic", "PC2_genetic", "PC3_genetic", "PC4_genetic"), by = c("Sample_ID"))%>%
  dplyr::select(-Host, Host)
sum(is.na(PC_env))
identical(rownames(gen.imp), PC_env$Sample_ID)  # TRUE
head(PC_env)

## Colinearity checked in earlier code
env %>%
  group_by(AnonSiteCode) %>%
  summarise(avg_wildness = mean(wildness, na.rm = TRUE),
            range_wildness = max(wildness, na.rm = TRUE) - min(wildness, na.rm = TRUE)
  )%>%
  arrange(desc(avg_wildness)) 

env %>%
  filter(AnonSiteCode == "OH1") %>%
  dplyr::select(Latitude, Longitude)
# Labels for graphing later!
labels<-read.csv("./data/too_large_files/all_freshwater_sdmpredictors/layerschoice_metadata.csv")
legend_labels <- setNames(labels$name, labels$layer_code)

#=========================        THE RDA STARTS HERE            ========================================================

#### Ordi2step model for forward selection #############
# Prep data for RDA- remove useless columns
# pred<-env[, 7:ncol(env), drop = FALSE]
# colnames(pred)
# # Null model
# RDA0 <- rda(gen.imp ~ 1, data = pred, scale=T)
# # Full model
# RDAfull <- rda(gen.imp ~ ., data = pred, scale=T)
# vif.cca(RDAfull)
# RDAfull <- rda(gen.imp ~ FW_dem_range + FW_flow_acc + FW_hydro_wavg_03 + FW_hydro_wavg_05 + 
#                  FW_hydro_wavg_09 + FW_hydro_wavg_10 + FW_hydro_wavg_17 + FW_lc_wavg_03 + 
#                  FW_lc_wavg_04 + FW_lc_wavg_05 + FW_lc_wavg_06 + FW_lc_wavg_07 + FW_lc_wavg_08 + 
#                  FW_lc_wavg_11 + FW_soil_wavg_01 + FW_soil_wavg_02 + FW_soil_wavg_04 + 
#                  FW_soil_wavg_07 + wildness + organicCar + Host, data = pred, scale = TRUE)

# Do forward selection/permutation test
# mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)
# Extract selected variables
# mod$anova
# added_rows <- mod$anova[grep("\\+ ", rownames(mod$anova)), ]
# selected_variables <- gsub("\\+ ", "", rownames(added_rows))
# selected_variables
# selected_formula <- as.formula(paste("gen.imp ~", paste(selected_variables, collapse = " + ")))
# selected_formula
# RDAprune <- rda(selected_formula, data = pred, scale=T)
# self_form<-as.formula(paste0("gen.imp ~ Host + FW_soil_wavg_07 + FW_soil_wavg_01 +
#     FW_dem_range + wildness + WC_alt +
#     FW_hydro_wavg_03  + FW_hydro_wavg_05 + FW_hydro_wavg_10
#     "))
#### End Ordi2step code ###############################

# pred_prune <- env[, c("FW_soil_wavg_07", "FW_soil_wavg_01", ##this has catchment area in!!
#                       "FW_dem_range", "wildness", "WC_alt",
#                       "FW_hydro_wavg_03", "FW_hydro_wavg_05", "FW_hydro_wavg_17","catchment_area_km", "Host")]
unique(pred_prune$Host)
pred_prune <- env[, c( "FW_dem_range", "wildness", "WC_alt", 
                       "FW_hydro_wavg_03", "FW_hydro_wavg_05", "FW_hydro_wavg_17","catchment_area_km", "Host")]
RDAprune <- rda(gen.imp ~ ., data = pred_prune, scale=T)
vif.cca(RDAprune)

host_as_cond <- rda(gen.imp ~ FW_dem_range + wildness + WC_alt +
                    FW_hydro_wavg_03 +FW_hydro_wavg_05+ FW_hydro_wavg_17+ catchment_area_km + Condition(Host), data = pred_prune, scale=T)

# host_as_cond <- rda(gen.imp ~ wildness + WC_alt +
#                       FW_hydro_wavg_03 +FW_hydro_wavg_05+ FW_hydro_wavg_17+ catchment_area_km + Condition(Host), data = pred_prune, scale=T)

vif.cca(host_as_cond)
pearl.rda<-host_as_cond
# pred_prune_with_pop_struct<-PC_env[, c("FW_soil_wavg_07", "FW_soil_wavg_01", 
#                                     "FW_dem_range", "wildness", "WC_alt",
#                                     "FW_hydro_wavg_03", "FW_hydro_wavg_05", "FW_hydro_wavg_17","catchment_area_km", 
#                                     "PC1_genetic","PC2_genetic", "PC3_genetic", "PC4_genetic", "Host")]

# RDA_PCCond <- rda(gen.imp ~ FW_soil_wavg_07 + FW_soil_wavg_01 + 
#                   FW_dem_range + wildness + WC_alt +
#                   FW_hydro_wavg_03 + FW_hydro_wavg_05 + FW_hydro_wavg_17 + catchment_area_km +
#                   Host + Condition(PC1_genetic + PC2_genetic + PC3_genetic + PC4_genetic), data = pred_prune_with_pop_struct, scale=T)
# pearl.rda<-RDA_PCCond

# Run forward selected model

pearl.rda[["call"]]
# pearl.rda <- readRDS("output/rda_with_pca_condition.rds") #rda(gen.imp ~ Host + FW_lc_wavg_03 + FW_lc_wavg_05 + wildness + FW_lc_wavg_11 + FW_dem_range + FW_soil_wavg_01 + FW_soil_wavg_07 + FW_lc_wavg_04 +   FW_lc_wavg_06 + FW_hydro_wavg_05 +FW_hydro_wavg_10 + FW_flow_acc +   FW_lc_wavg_08 + Condition(PC1_genetic + PC2_genetic + PC3_genetic + PC4_genetic), data = pred, scale=T)
summary(pred_prune)

# Check vif
vif.cca(pearl.rda) 
pearl.rda
summary(pearl.rda)
# library(stats)
# alias(pearl.rda) #no aliasing required
screeplot(pearl.rda)

# # Normality of residuals ######
# residuals_rda <- residuals(pearl.rda)
# qqnorm(residuals_rda, main = "Q-Q Plot of RDA Residuals")
# qqline(residuals_rda, col = "blue") #Not great, and taking out soil param made it worse :'(
# hist(residuals_rda, main = "Histogram of RDA Residuals", xlab = "Residuals", breaks = 30)
# 
# # Homogeneity of Variance
# # Plot residuals vs. fitted values
# plot(fitted(pearl.rda), residuals_rda, main = "Residuals vs Fitted",
#      xlab = "Fitted Values", ylab = "Residuals")
# abline(h = 0, col = "red")

# install.packages("nortest")
# library(nortest)
# ad.test(residuals_rda) #########


#Check RDA significance (takes time so why we save the model)
# signif.full <- anova.cca(pearl.rda, parallel=getOption("mc.cores"), permutations = how(nperm=999)) # default is permutation=999
# signif.axis <- anova.cca(pearl.rda, by="axis", parallel=getOption("mc.cores"), permutations = how(nperm=999)) #takes 50
# signif.term <- anova.cca(pearl.rda, by="term", parallel=getOption("mc.cores")) #takes time

signif.full #takes time #full model is significant >0.001
signif.axis #<0.001 for RDA1:5, RDA6 is <0.01
signif.term #<0.001 for RDA1:5, RDA6 is <0.01
RsquareAdj(pearl.rda)


save.image(file = "data/confidential/RDA_host_condition_updated.RData")
load("data/confidential/RDA_host_condition_updated.RData")

# load("code/R_code/backup_workspace_images/HOST_CONDITION.RData")


pearl.rda[["call"]]
summary(eigenvals(pearl.rda, model = "constrained"))

# Extract the eigenvalues and calculate the proportion of variance explained
eigenvalues <- eigenvals(pearl.rda, model = "constrained")
eigenvalues
variance_explained <- eigenvalues / sum(eigenvalues)
varience_for_writeup<-round(variance_explained * 100, 3)  
varience_for_writeup
# Extract the proportions for RDA1 and RDA2
rda1_prop <- round(variance_explained[1] * 100, 1)  
rda2_prop <- round(variance_explained[2] * 100, 1) 
pearl.rda$terminfo$xlev
#Predictor variance percent for writeup
 (signif.term$Variance / sum(signif.term$Variance)) * 100
signif.term %>%
  mutate(Proportion = (Variance / sum(Variance[-nrow(.)])) * 100) %>% 
  mutate(Predictors = rownames(.)) %>%  
  dplyr::select(Predictors, Proportion) %>%  
  pivot_wider(names_from = Predictors, values_from = Proportion)  %>% 
  as.data.frame()



# Just clim RDA
# clim.RDA <- rda(gen.imp ~ FW_hydro_wavg_03 + FW_hydro_wavg_05 +  FW_hydro_wavg_17 , data = pred, scale=T)
# vif.cca(clim.RDA)
# pearl.rda<-clim.RDA

# ##Another way to term all permutations of factors (BUT perhaps less suited to genomic rda)
# spec.RDA<-step(RDAfull, scope=formula(RDAfull), test="perm")
# summary(spec.RDA)
#===========================================
## So back to the original self-selected variables #####

plot(pearl.rda, scaling=3)          # default is axes 1 and 2
plot(pearl.rda, choices = c(1, 3), scaling=3)  # axes 1 and 3
# install.packages("ggord")
# library(ggord)
pearl.rda$CCA
#######I have moved this section up in the code, hopefully it all works here!
#########################

rdadapt <- function(rda,K) ###from Capblanq tutorial, needed for the extract_rda_stats function
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

extract_rda_stats <- function(rda_obj, gen_imp, num_axes = 2, scaling = 3) {
  # Extract species loadings on specified number of RDA axes with specified scaling
  loadings <- scores(rda_obj, choices = 1:num_axes, display = "species", scaling = scaling)
  
  # Extract SNP names from gen.imp (assuming gen.imp is a dataframe with SNP data)
  snp_names <- colnames(gen_imp)
  
  # Extract SNP positions
  snp_positions <- data.frame(
    pos = 1:length(snp_names),
    snp = snp_names
  )
  
  # Prepare a matrix to hold loadings for all axes
  loading_matrix <- matrix(NA_real_, nrow = length(snp_names), ncol = num_axes)
  
  # Fill the loading matrix with loadings
  loading_matrix[, 1:num_axes] <- loadings
  
  # Calculate p-values or q-values using rdadapt function (assuming rdadapt is available)
  rdadapt_results <- rdadapt(rda_obj, num_axes)
  p_values <- rdadapt_results$p.values
  q_values <- rdadapt_results$q.values
  
  # Combine SNP names, positions, loadings, p-values, and q-values into a dataframe
  loadings_df <- data.frame(
    SNP = snp_positions$snp,
    Position = snp_positions$pos,
    loading_matrix,
    p_values = p_values,
    q_values = q_values
  )
  
  # Rename columns to RDA1, RDA2, RDA3, ..., and include p-values and q-values
  colnames(loadings_df)[-(1:2)] <- c(paste0("RDA", 1:num_axes), "p_values", "q_values")
  
  return(loadings_df)
}

## For dynamic plot saving
run_name<-"Host_as_condition"
dir <- getwd() 

# Extract RDA statistics for X axes with scaling = 3
num_axis<-6
rda_stats <- extract_rda_stats(rda_obj = pearl.rda, gen_imp = gen.imp, num_axes = num_axis, scaling = 3)
head(rda_stats)
summary(rda_stats)
# Extract RDA statistics for 2 axes with normal scaling = 1
rda_stats_normal <- extract_rda_stats(rda_obj = pearl.rda, gen_imp = gen.imp, num_axes = num_axis, scaling = 3) #means no scaling done
head(rda_stats_normal)

# Extract the scores for environmental variables (the arrows and lines)
env_scores_df <- as.data.frame(scores(pearl.rda, display = "bp", scaling = 3, choices=c(1,2)))
env_scores_df$variable <- rownames(env_scores_df)
unique_predictors <- unique(env_scores_df$variable)
columns_without_layer_code <- setdiff(unique_predictors, labels$layer_code)
columns_without_layer_code
# Add new labels to labels df based on ones missing
labels<-labels%>% dplyr::select(layer_code,name)
labels <- rbind(labels, data.frame(layer_code = columns_without_layer_code, name = columns_without_layer_code))
# Define label replacements
label_map <- c(
  "wildness" = "Wildness",
  "WC_alt" = "Altitude",
  "catchment_area_km" = "Catchment area (km²)"
)
# Update names in labels
labels <- labels %>%
  mutate(name = ifelse(layer_code %in% names(label_map),
                       label_map[layer_code],
                       name))
labels$name <- ifelse(grepl("FW_hydro", labels$layer_code),
                      toupper(gsub(".*_(\\d+)$", "BIOCLIM\\1", labels$layer_code)),
                      labels$name)
labels$name <- gsub("BIOCLIM0([1-9])", "BIOCLIM\\1", labels$name)

env_scores_df <- left_join(env_scores_df, labels, by = c("variable" = "layer_code"))
head(env_scores_df)
unique(env_scores_df$name)
site_scores <- as.data.frame(scores(pearl.rda, display = "sites", choices=1:4, scaling = 3))
site_scores$Sample_ID <- rownames(site_scores)
rda_env<-left_join(env,site_scores, by = "Sample_ID") %>% filter(!is.na(RDA1))

##### Plots ######
        
# Plotting with ggplot2
base_plot<-ggplot(rda_stats, aes(x = RDA1, y = RDA2)) +
  geom_point(color = "grey90", size = 1.4) +  # Plot all SNPs with normal scaling in grey
  labs(x =paste0("RDA1 (",rda1_prop, "%)"), y =paste0("RDA2 (", rda2_prop, "%)"), title = "RDA1 vs RDA2 for SNP Loadings") +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), linewidth = 0.6) +  # Horizontal line at y=0
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), linewidth = 0.6) +  # Vertical line at x=0
  geom_segment(data = env_scores_df, aes(x = 0, y = 0, xend = RDA1 * 7, yend = RDA2 * 7), 
               colour = "black", linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +  # Arrows from origin to scaled env scores
  # geom_text_repel(data = env_scores_df, aes(x = 5 * RDA1, y = 5 * RDA2, label = str_wrap(name, width = 35)), 
  #                 size=3,
  #                 min.segment.length = Inf,
  #                 position = position_nudge_center(0.2, 0.1, 0, 0),
  #                 max.overlaps = Inf, # Adjust as necessary for more flexibility
  #                 box.padding = 0.1, # Adjust padding around labels
  #                 point.padding = 0.1) +  # Adjust padding around points
   theme_classic() +
  theme(legend.position = "right",  # Position legend at the bottom
        legend.key.size = unit(0.5, "cm")) + # Adjust legend key size
  scale_x_continuous(expand = expansion(mult = 0.5))
base_plot

bigger_rda<-ggplot(rda_stats, aes(x = RDA1, y = RDA2)) +
  geom_point(color = "grey90", size = 1.4) +  # Plot all SNPs with normal scaling in grey
  geom_point(data = rda_env, aes(color = AnonSiteCode), size = 3, alpha=1)+ # Population loadings
  labs(x =paste0("RDA1 (",rda1_prop, "%)"), y =paste0("RDA2 (", rda2_prop, "%)"), title = "RDA1 vs RDA2 for SNP Loadings") +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Horizontal line at y=0
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Vertical line at x=0
  geom_segment(data = env_scores_df, aes(x = 0, y = 0, xend = RDA1 * 10, yend = RDA2 * 10),                colour = "black", linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +  # Arrows from origin to scaled env scores
  geom_text_repel(data = env_scores_df, aes(x = 10.5 * RDA1, y = 10.5 * RDA2, label = str_wrap(name, width = 35)), 
                  size=3.5,
                  min.segment.length = Inf,
                  position = position_nudge_center(0.2, 0.1, 0, 0),
                  max.overlaps = Inf, # Adjust as necessary for more flexibility
                  box.padding = 0.1, # Adjust padding around labels
                  point.padding = 0.1) +  # Adjust padding around points
  theme_classic()+
  theme(legend.position = "right",  # Position legend at the bottom
        legend.key.size = unit(0.5, "cm"))+  # Adjust legend key size geom_point(data = rda_env, aes(color = Population), size = 3)+  # Overlay candidate SNPs with special scaling
  labs(colour="Anon Population")+  coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) 
  # ggforce::geom_mark_ellipse(data = rda_env %>% filter(AnonSiteCode %in% c("WS1", "OH1", "ER1", "WS3", "WS2", "WR4", "WI1", "WI2", "WI4")), 
  #                            aes(x = RDA1, y = RDA2, group = AnonSiteCode, label = AnonSiteCode), 
  #                            fill = NA, 
  #                            expand = unit(2, "mm"), 
  #                            label.fontsize = 7,
  #                            label.buffer = unit(0, 'mm'),
  #                            label.margin = margin(0,0,0,0),
  #                            label.fontface = "italic",
  #                            label.fill = NA,
  #                            con.type = "none") 
bigger_rda
ggsave(plot= bigger_rda,
       width = 9.5,
       height = 7,
       units = c( "in"),
       dpi = 600,
       filename = paste0("./output/figures/", run_name, " bigger base rda.png"))
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/Host_as_condition bigger base rda pca2.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned//output/figures/Host_as_condition bigger base rda pca_labelled.png")
# ggsave(filename = paste0("./output/figures/", run_name, " bigger base rda pca2.png"), plot = bigger_rda)
shell.exec(paste0(dir,"/output/figures/", run_name, " bigger base rda pca_labelled.png"))


base_plot +
  geom_point(data = rda_env, aes(x = RDA1, y = RDA2, color = Population), size = 3) +  # Overlay candidate SNPs with special scaling
  labs(colour = "Population") +  # Legend label
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5))  # Adjust coordinate limits if needed

base_plot +
  geom_point(data = rda_env, aes(x = RDA1, y = RDA2, color = FW_soil_wavg_07), size = 3) +  # Overlay candidate SNPs with special scaling
  labs(colour = "Cation") +  # Legend label
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5))  # Adjust coordinate limits if needed
#######

######     Identify candidate SNPs involved in local adaptation #################################################
################################################################################################################
# Look at distribution of the data
ggplot() +
  geom_line(aes(x=c(1:length(pearl.rda$CCA$eig)), y=as.vector(pearl.rda$CCA$eig)), linetype="dotted",
            size = 1.5, color="darkgrey") +
  geom_point(aes(x=c(1:length(pearl.rda$CCA$eig)), y=as.vector(pearl.rda$CCA$eig)), size = 3,
             color="darkgrey") +
  scale_x_discrete(name = "Ordination axes", limits=c(1:9)) +
  ylab("Inertia") +
  theme_bw()

load.rda <- scores(pearl.rda, choices=c(1:7), display="species")  # Species scores for the first three constrained axes
wrap_plots(lapply(1:7, function(i) {
  ggplot(data.frame(loading = load.rda[, i]), aes(x = loading)) +
    geom_histogram(binwidth = 0.05, fill = "grey80", colour = "black") +
    labs(x = paste("RDA", i, "Loadings"), y = "Frequency")+
    theme_classic()
  # xlim(-0.4, 0.5) +  # Set x-axis limits
  # ylim(0, 2500)  # Set y-axis limits
}), ncol = 3) + plot_annotation(title = paste0("Run name: '", run_name, "'")) 

dim(load.rda)

rda_stats <- extract_rda_stats(rda_obj = pearl.rda, gen_imp = gen.imp, num_axes = 6, scaling = 3)

# Top 1% snps
# Extract the top 1% of SNPs based on RDA1 loadings
# top_one_percent_snps_rda1 <- rda_stats %>%
#   arrange(desc(.)) %>%
#   filter(RDA1 > quantile(RDA1, 0.99))

top_one_percent_snps <- rda_stats %>%
  mutate(across(starts_with("RDA"), ~ . > quantile(., 0.99), .names = "{.col}_is_top1")) %>%
  rowwise() %>%
  mutate(axis = paste0(names(across(ends_with("_is_top1")))[which(c_across(ends_with("_is_top1")))], collapse = ", ")) %>%
  mutate(axis = sub("_is_top1", "", axis)) %>%
  filter(axis != "") %>%
  dplyr::select(-ends_with("_is_top1")) %>%
  ungroup()
nrow(top_one_percent_snps)
# 3SD from the mean for each RDA axis
# This is just for gettig the outliers
# threeSD_snps <- rda_stats %>% 
#   mutate(across(starts_with("RDA"), ~ (. - mean(.)) / sd(.))) %>%  # Z-score calculation
#   filter(if_any(starts_with("RDA"), ~ abs(.) > 3))  # Filter SNPs where Z-scores are greater than 3
# Subset SNPs 3SD from their RDA axis
threeSD_snps <- rda_stats %>%
  mutate(across(starts_with("RDA"), ~ (. - mean(.)) / sd(.), .names = "{.col}")) %>%
  pivot_longer(cols = starts_with("RDA"), 
               names_to = "axis", 
               values_to = "z_score") %>%
  filter(abs(z_score) > 3) %>% # Add the corresponding loading for each SNP and axis
  left_join(rda_stats %>% dplyr::select(SNP, starts_with("RDA")), by = "SNP") %>%
  as.data.frame()
nrow(threeSD_snps)
# # Get predictor correlations for each SNP- function later on
# cor_df <- cor(gen.imp[, threeSD_snps$SNP], pred_prune %>% select(-Host))  %>% as.data.frame()
# pred_Host_expand <- bind_cols(pred_prune %>% select(-Host), model.matrix(~ Host - 1, data = pred_prune)[, -1])
# # cor_df <- cor(gen.imp[, threeSD_snps$SNP], pred_Host_expand)  %>% as.data.frame()
# # Bind correlations with SNPs
# threeSD_snps_bind <- cbind(threeSD_snps, cor_df)
# 
# # Find the most correlated predictor for each SNP
# # Define column range for correlation calculation
# column_window <- which(colnames(threeSD_snps_bind) %in% colnames(pred_Host_expand))
# 
# # Vectorised function to find the most correlated predictor and its value
# results <- apply(threeSD_snps_bind[, column_window], 1, function(row) {
#   max_index <- which.max(abs(row))
#   c(predictor = colnames(threeSD_snps_bind)[column_window][max_index], correlation = abs(row[max_index]))
# })
# threeSD_snps_bind[c("predictor", "correlation")] <- t(results)
# table(threeSD_snps_bind$predictor) #doesnt have to be all predictors-some arnet here because they arent "most corr var"
# 
# threeSD_snps_bind$correlation <- as.numeric(threeSD_snps_bind$correlation)
# 
# threeSD_snps_bind %>%
#   group_by(predictor) %>%
#   summarise(
#     count = n(),
#     average_correlation = mean(correlation, na.rm = TRUE)
#   )
# threeSD_snps_bind %>%
#   group_by(axis) %>%
#   summarise(
#     count = n()
#   )
# 
# 
# head(threeSD_snps_bind)
# head(cand)


## Loci based on p-value
# P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rda_stats$p_values) 
thres_env

outliers_below_pvalue <- rda_stats %>%
  # arrange(desc(p_values)) %>%
  filter(p_values < thres_env)
nrow(outliers_below_pvalue)

# Out_qvalue <- rda_stats %>%
#   # arrange(desc(p_values)) %>%
#   filter(q_values < 0.1)
# nrow(Out_qvalue)
Out_qvalue <- rda_stats %>%
  # arrange(desc(p_values)) %>%
  filter(q_values < 0.05)
nrow(Out_qvalue)
## Get most correlated predictor
find_most_correlated_predictor <- function(gen_imp, snp_df, pred_host_expand) {
  cor_df <- cor(gen_imp[, snp_df$SNP], pred_host_expand) %>% as.data.frame()
  snp_bind <- cbind(snp_df, cor_df)
  column_window <- which(colnames(snp_bind) %in% colnames(pred_host_expand))
  results <- apply(snp_bind[, column_window], 1, function(row) {
    max_index <- which.max(abs(row))
    c(predictor = colnames(snp_bind)[column_window][max_index], correlation = abs(row[max_index]))
  })
  snp_bind[c("predictor", "correlation")] <- t(results)
  snp_bind$correlation <- as.numeric(snp_bind$correlation)
  return(snp_bind)
}
Out_3SD <- find_most_correlated_predictor(gen.imp, threeSD_snps, pred_prune %>% dplyr::select(-Host))
Out_pv<-find_most_correlated_predictor(gen.imp, outliers_below_pvalue, pred_prune %>% dplyr::select(-Host))
Out_top1<-find_most_correlated_predictor(gen.imp, top_one_percent_snps, pred_prune %>% dplyr::select(-Host))
Out_qv<-find_most_correlated_predictor(gen.imp, Out_qvalue, pred_prune %>% dplyr::select(-Host))

write.csv(Out_qv, "output/rda_qvalue_outlier_snps.csv", row.names=F)


head(Out_qv)
forkara<-Out_qv %>% dplyr::select(SNP,FW_dem_range, wildness, WC_alt, FW_hydro_wavg_03, FW_hydro_wavg_05, FW_hydro_wavg_17, catchment_area_km,         predictor, correlation)
Out_3SD %>%
  group_by(predictor) %>%
  summarise(
    count = n(),
    average_correlation = mean(correlation, na.rm = TRUE)
  )
Out_qv %>%
  group_by(predictor) %>%
  summarise(
    count = n(),
    average_correlation = mean(correlation, na.rm = TRUE)
  )
out_qv_summary <- Out_qv %>%
  group_by(predictor) %>%
  summarise(
    count = n(),
    average_correlation = mean(correlation, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = -predictor, names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = predictor, values_from = value)
out_qv_summary
write.csv(out_qv_summary, "output/outlier_qv_countsperabsolutetopcorrelation.csv", row.names = FALSE)
Out_qv %>%
  # filter(predictor == "FW_hydro_wavg_03") %>%
  ggplot(aes(FW_hydro_wavg_03)) +
  geom_histogram(fill = "deeppink", color = "white", binwidth = 0.01) +
  labs(title = "Correlation Distribution (FW_hydro_wavg_03)", x = "Correlation") +
  xlim(c(-.5, .5))+
  theme_classic()
Out_qv %>%
  # filter(predictor == "FW_hydro_wavg_03") %>%
  summarise(
    count_positive = sum(FW_hydro_wavg_03 > 0),
    count_negative = sum(FW_hydro_wavg_03 < 0)
  ) %>%
  print()



# outliers_qvalue <- outliers_function(rda_stats, gen.imp, pred_prune, subset_type = "q_value", q_value_threshold = 0.05)
head(outliers_qvalue)
head(Out_qv)
pie_out3SD<-ggplot(Out_3SD %>%
         group_by(predictor) %>%
         summarise(
           count = n(),
           average_correlation = mean(correlation, na.rm = TRUE)
         ), aes(x = "", y = count, fill = predictor)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "3 SD approach Predictors by Count") +
  theme_void() +
  theme(legend.title = element_blank())
pie_outqv<-ggplot(Out_qv %>%
                    group_by(predictor) %>%
                    summarise(
                      count = n(),
                      average_correlation = mean(correlation, na.rm = TRUE)
                    ), aes(x = "", y = count, fill = predictor)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "q-value approach of Predictors by Count") +
  theme_void() +
  theme(legend.title = element_blank())
pie_outqv
com<-pie_out3SD+pie_outqv
ggsave(plot= pie_outqv,
       width = 9.5,
       height = 7,
       units = c( "in"),
       dpi = 600,
       filename = paste0("./output/figures/", run_name, " q-value approach of Predictors by Count.png"))
ggsave(plot= com,
       width = 12,
       height = 7,
       units = c( "in"),
       dpi = 600,
       filename = paste0("./output/figures/", run_name, " SD and q-value approach of Predictors by Count.png"))
# Pull SNPs significant with pcadapt

## Run PCAdapt script
pcadapt_snps<-readRDS("output/pcadapt_significant_snps_bonferroni_adjusted_under_0.05.rds")
head(pcadapt_snps)
head(rda_stats)
# Filter rda_stats based on SNP names in pcadapt_snps
pcadapt_snps_with_rda_pvals <- rda_stats %>%
  filter(sub("_[ACGT]$", "", SNP) %in% pcadapt_snps$SNP)
head(pcadapt_snps_with_rda_pvals)


#### Lots of manhatten plots!
man_pcad<-ggplot(rda_stats, aes(x = Position, y = -log10(p_values))) +
  geom_point(color = "grey80") +  # Plot all SNPs with normal scaling
  # geom_point(data = Out_3SD, aes(colour =axis ), size = 3) +  # Overlay candidate SNPs with special scaling
  geom_point(data = pcadapt_snps_with_rda_pvals, colour = "#14a8ff", size = 3) +  # Overlay candidate SNPs with special scaling
  labs(x = "Position in Genome", y = expression(-log[10](italic("p")["RDA"])), title= paste("Outliers based on PCADAPT ( N =",nrow(pcadapt_snps_with_rda_pvals),")")) +
  scale_color_discrete(name = "Predictor") +  # Adjust legend title if needed
  theme_classic() +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  # Ensure x-axis starts at 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22))  # Ensure y-axis starts at 0man<-ggplot(rda_stats, aes(x = Position, y = -log10(p_values))) +

man_pcad

manqv<-ggplot(rda_stats, aes(x = Position, y = -log10(p_values))) +
  geom_point(color = "grey80") +  
  geom_point(data = Out_qv, colour = "deeppink", size = 3) +  
  labs(x = "Position in Genome", y = expression(-log[10](italic("p")["RDA"])), title= paste("Outliers based on q-values < 0.05 ( N =", nrow(Out_qv),")")) +
  scale_color_discrete(name = "Predictor") +  # Adjust legend title if needed
  theme_classic() +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22))  
manqv
ggsave(plot= manqv,
       width = 9,
       height = 3,
       units = c( "in"),
       dpi = 600,
       filename = paste0("./output/figures/", run_name, "qval_manhatten.png"))
shell.exec(paste0(dir,"/output/figures/", run_name, "qval_manhatten.png"))

manSD<-ggplot(rda_stats, aes(x = Position, y = -log10(p_values))) +
  geom_point(color = "grey80") +  
  # geom_point(data = Out_3SD, aes(colour =axis ), size = 3) +  
  geom_point(data = Out_3SD, colour = "deeppink", size = 3) +  
  labs(x = "Position in Genome", y = expression(-log[10](italic("p")["RDA"])), title= paste("Outliers based on 3 SD ( N =", nrow(Out_3SD),")")) +
  scale_color_discrete(name = "Predictor") +  
  theme_classic() +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22))  
manSD
manpv<-ggplot(rda_stats, aes(x = Position, y = -log10(p_values))) +
  geom_point(color = "grey80") +  
  geom_point(data = Out_pv, colour = "deeppink", size = 3) +  
  labs(x = "Position in Genome", y = expression(-log[10](italic("p")["RDA"])), title = paste("Outliers based on p-values < 0.05 ( N =", nrow(Out_pv), ")")) +
  scale_color_discrete(name = "Predictor") +  
  theme_classic() +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22))  
manpv
# man1p<-ggplot(rda_stats, aes(x = Position, y = -log10(p_values))) +
#   geom_point(color = "grey80") +  # Plot all SNPs with normal scaling
#   geom_point(data = Out_top1, colour = "deeppink", size = 3) +  # Overlay candidate SNPs with special scaling
#   # geom_point(data = Out_top1, aes(colour = axis ), size = 3) +  # Overlay candidate SNPs with special scaling
#   labs(x = "Position in Genome", y = "-log10(p.values)", title= "Outliers based on top 1% across 6 Axes") +
#   scale_color_discrete(name = "Predictor") +  # Adjust legend title if needed
#   theme_classic() +
#   geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  # Ensure x-axis starts at 0
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 22))  # Ensure y-axis starts at 0
# man1p
comp<-manSD / manpv / manqv / man_pcad +
  plot_annotation(tag_levels = 'a', tag_suffix = ")")
comp

ggsave(plot= comp,
       width = 8,
       height = 10,
       units = c( "in"),
       dpi = 600,
       filename = paste0("./output/figures/", run_name, "comparing_outlier_techniques.png"))
shell.exec(paste0(dir,"/output/figures/", run_name, "comparing_outlier_techniques.png"))

manqv_qv<-ggplot(rda_stats, aes(x = Position, y = -log10(q_values))) +
  geom_point(color = "grey80", size=2, alpha=1) +  
  geom_point(data = Out_qv, colour = "deeppink", size = 2) +  
  labs(x = "Position in Genome", y = expression(-log[10](italic("q")["RDA"])), title= paste("Outliers based on q-values < 0.05 ( N =", nrow(Out_qv),")")) +
  scale_color_discrete(name = "Predictor") +  # Adjust legend title if needed
  theme_classic() +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = gray(.60), size=0.6) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 21), 
                     breaks = c(seq(0, 20, by = 5), -log10(0.05)),  # Add default breaks + threshold
                     labels = c(seq(0, 20, by = 5), "q = 0.05")     # Label the threshold
  ) 
manqv_qv
ggsave(plot= manqv_qv,
       width = 9,
       height = 3,
       units = c( "in"),
       dpi = 600,
       filename = paste0("./output/figures/", run_name, "qval_qval_manhatten.png"))
shell.exec(paste0(dir,"/output/figures/", run_name, "qval_qval_manhatten.png"))


compare_pcadapt_with_rda<-ggplot(rda_stats, aes(x = Position, y = -log10(p_values))) +
  geom_point(aes(colour = "All SNPs"), size = 2) +  # Assign a label for all SNPs
  geom_point(data = Out_pv, aes(colour = "RDA (p-value based)"), size = 4) +
  geom_point(data = Out_3SD, aes(colour = "RDA (3-SD based)"), size = 3, alpha = 1) +
  geom_point(data = pcadapt_snps_with_rda_pvals, aes(colour = "pcadapt"), size = 2, pch = 17) +
  geom_hline(yintercept = -log10(thres_env), linetype = "dashed", aes(colour = "Threshold"), size = 0.6) +
  labs(
    x = "Position in Genome",
    y = "-log10(p.values)",
    colour = "Outlier SNPs"  # Title for the legend
  ) +
  theme_classic() +
  scale_colour_manual(
    values = c(
      "All SNPs" = "grey80",
      "RDA (p-value based)" = "deeppink",
      "RDA (3-SD based)" = "#ff8ac9",
      "pcadapt" = "#111",
      "Threshold" = "grey"
    )
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22))
pcadapt_snps_with_rda_pvals %>% #how many pcaadapt snps are there?
  #filter(-log10(p_values) > -log10(thres_env)) %>%
  nrow()
compare_pcadapt_with_rda
ggsave(compare_pcadapt_with_rda,width = 8.5,
       height = 3.5,
       units = c( "in"),
       dpi = 600,
       filename = paste0("./output/figures/", run_name, "manhatten_with_both_pcadapt_and_RDA_SNPS.png"))

# Combine and pivot data for easier handling
rda_combined <- rda_stats %>%
  mutate(type = "All SNPs") %>%
  bind_rows(
    Out_3SD %>% mutate(type = "Outliers"),
    Out_pv %>% mutate(type = "Top outliers"),
    pcadapt_snps_with_rda_pvals %>% mutate(type = "pcadapt outliers")
    
  ) %>%
  # Ensure duplicates are labelled as "Outliers p-value"
  group_by(Position) %>%
  arrange(desc(type)) %>%  # Prioritise "Outliers p-value" in case of duplicates
  slice(1) %>%  # Keep only the first occurrence based on the priority
  ungroup()
 ggplot(rda_combined, aes(x = Position, y = -log10(p_values), color = type)) +
  geom_point(size = 2, shape=16) +
  labs(x = "Position in Genome", y = "-log10(p-values)", title = "Outliers based on p-values across 6 Axes") +
  scale_color_manual(
    name = "Legend",
    values = c("All SNPs" = "grey80", "Outliers" = "#ff8ac9", "Top outliers" = "deeppink", "pcadapt outliers"="#14a8ff")
  ) +
  theme_classic() +
  geom_hline(yintercept = -log10(thres_env), linetype = "dashed", color = "gray80", size = 0.6) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22))

ggplot() +
  geom_point(data = rda_combined %>% filter(type == "All SNPs"), 
             aes(x = Position, y = -log10(p_values)), 
             color = "grey80", size = 2, shape = 16) +  # Plot all SNPs at the bottom
  geom_point(data = rda_combined %>% filter(type == "Outliers"), 
             aes(x = Position, y = -log10(p_values)), 
             color = "#ff8ac9", size = 2.5, shape = 16) +  # Plot outliers on top of all SNPs
  geom_point(data = rda_combined %>% filter(type == "Top outliers"), 
             aes(x = Position, y = -log10(p_values)), 
             color = "deeppink", size = 3, shape = 16) +  # Plot top outliers at the very top
  labs(x = "Position in Genome", y = "-log10(p-values)", title = "Outliers based on p-values across 6 Axes") +
  theme_classic() +
  geom_hline(yintercept = -log10(thres_env), linetype = "dashed", color = "gray80", size = 0.6) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22))
man_com<-ggplot() +
  geom_point(data = rda_combined %>% filter(type == "All SNPs"), 
             aes(x = Position, y = -log10(p_values), color = "All SNPs"), 
             size = 2, shape = 16) +  # Plot all SNPs at the bottom
  geom_point(data = rda_combined %>% filter(type == "Outliers"), 
             aes(x = Position, y = -log10(p_values), color = "Outliers"), 
             size = 2.5, shape = 16) +  # Plot outliers on top of all SNPs
  geom_point(data = rda_combined %>% filter(type == "Top outliers"), 
             aes(x = Position, y = -log10(p_values), color = "Top outliers"), 
             size = 3, shape = 16) +  # Plot top outliers at the very top
  labs(x = "Position in Genome", y = "-log10(p-values)") +
  scale_color_manual(
    name = " ",
    values = c("All SNPs" = "grey80", "Outliers" = "#ff9dd2", "Top outliers" = "deeppink")
  ) +
  theme_classic() +
  geom_hline(yintercept = -log10(thres_env), linetype = "dashed", color = "gray80", size = 0.6) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22))
man_com
ggsave(man_com,width = 8.5,
       height = 3.5,
       units = c( "in"),
       dpi = 600,
       filename = paste0("./output/figures/", run_name, "manhatten_with_both_pval_and_3SD_SNPS.png"))

# Add a label column for top outliers (e.g., top 5 based on p-values)
top_outliers <- rda_combined %>%
  filter(type != "All SNPs") %>%
  arrange(desc(-log10(p_values))) %>%
  slice_head(n = 5)

ggsave(filename = paste0("./output/figures/", run_name, " comp_outlier_methods.png"), plot = comp)

## Plot rda plot
ggplot(rda_stats_normal, aes(x = RDA1, y = RDA2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Horizontal line at y=0
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Vertical line at x=0
  geom_point(color="grey90") +  # Plot all SNPs with normal scaling
  geom_point(data = Out_3SD, size = 3, colour="deeppink") +  # Overlay candidate SNPs with special scaling
  labs(x =paste0("RDA1 (",rda1_prop, "%)"), y =paste0("RDA2 (", rda2_prop, "%)"), title = "RDA1 vs RDA2 for SNP Loadings") +
  scale_color_discrete(name = "Predictor") +  # Adjust legend title if needed
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  geom_segment(data = env_scores_df, aes(xend = RDA1*6, yend = RDA2*6, x = 0, y = 0), colour = "black", linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_text_repel(data = env_scores_df, aes(x = 6 * RDA1, y = 6 * RDA2, label = str_wrap(name, width =20)), 
                  size=3.5,
                  min.segment.length = Inf,
                  position = position_nudge_center(0.2, 0.1, 0, 0),
                  max.overlaps = Inf, # Adjust as necessary for more flexibility
                  box.padding = 0.1, # Adjust padding around labels
                  point.padding = 0.1) +  # Adjust padding around points
  # geom_text(data = env_scores_df, aes(x = RDA1, y = RDA2, label = variable), 
  #           color = "blue", vjust = 1.5, hjust = 1.5)
  # theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11))
  scale_x_continuous(expand = expansion(mult = 0.5))+
  theme_classic() +
  theme(legend.position = "bottom",  # Position legend at the bottom
        legend.spacing.x = unit(0.1, "cm"),
        legend.key.size = unit(0.1, "cm"))  # Adjust legend key size




##### Investigating overlap of outliers

head(threeSD_snps_bind)
head(top_one_percent_snps_rda1)
head(outliers_below_pvalue)
dim(outliers_below_pvalue)

# Looking at overlaps between methods
snp_lists <- list(
  Out_3SD = Out_3SD$SNP,
  Out_pv = Out_pv$SNP,
  Out_qv = Out_qv$SNP,
  pcadapt = pcadapt_snps_with_rda_pvals$SNP
)
length(Reduce(intersect, snp_lists))       # All lists
length(intersect(snp_lists$pcadapt, snp_lists$Out_qv))  # pcadapt vs Out_qv

# All pairwise overlaps
overlaps <- combn(names(snp_lists), 2, function(x) {
  intersect_snps <- intersect(snp_lists[[x[1]]], snp_lists[[x[2]]])
  data.frame(
    Method1 = x[1],
    Method2 = x[2],
    Overlap = length(intersect_snps),
    SNPs = paste(intersect_snps, collapse = ", ")
  )
}, simplify = FALSE)

overlaps_df <- do.call(rbind, overlaps)
print(overlaps_df)
common_snps <- Reduce(intersect, snp_lists)
common_snps


save.image(file = "data/confidential/RDA_host_condition_updated.RData")

#######################################################
# Base R way of filtering SNPS based on 3SD from mean
load.rda <- scores(pearl.rda, choices=c(1:6), display="species")  # Species scores for the first six constrained axes

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) # 31     # so identifying outliers *3 sd* away from the mean
cand2 <- outliers(load.rda[,2],3) # 7
cand3 <- outliers(load.rda[,3],3) # 55
cand4 <- outliers(load.rda[,4],3) # 72
cand5 <- outliers(load.rda[,5],3) # 17
cand6 <- outliers(load.rda[,6],3) # 57

ncand <- length(cand1) + length(cand2) + length(cand3) + length(cand4) + length(cand5) + length(cand6)
ncand


cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))
cand5 <- cbind.data.frame(rep(5,times=length(cand5)), names(cand5), unname(cand5))
cand6 <- cbind.data.frame(rep(6,times=length(cand6)), names(cand6), unname(cand6))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- colnames(cand4)<- colnames(cand5)<- colnames(cand6)<- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3, cand4, cand5, cand6)
cand$snp <- as.character(cand$snp)
pred2 <- cbind(pred_prune, model.matrix(~ Host - 1, data=pred_prune)[, -1]) # Removing the intercept
pred2$Host <- NULL

foo <- matrix(nrow=(ncand), ncol=ncol(pred2))  # number columns for number predictors
# colnames(foo) <- c("bio2_MDR", "bio_3Isotherm", "bio5_MaxT_WarmMonth", "bio9_Mean_T_Dry_Quarter", "bio10Mean_T_Warm_Quarter", "bio17precip_driestq")
# colnames(foo)<-colnames(pca)
 colnames(foo)<-colnames(pred2)
length(cand$snp)
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2] #grabs snp cand name from cand 
  snp.gen <- gen.imp[,nam] #looks for snp in original gen.imp file
  foo[i,] <- apply(pred2,2,function(x) cor(x,snp.gen)) #adds pred info for snp to empty foo file
}


cand <- cbind.data.frame(cand,foo)  
head(cand) #so now should have snp and info for each predictor

# Investigate the candidates
cand$snp[duplicated(cand$snp)]
length(cand$snp[duplicated(cand$snp)])  # 89 duplicate detections BECAUSE THIS IS LOOKING FOR DUPLICATE SNP NAMES

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # 27 duplicates on axis 1
table(foo[foo[,1]==2,2]) #  7 duplicates on axis 2
table(foo[foo[,1]==3,2]) # 55 duplicates on axis 3
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

# column_window<-paste0("4:",(ncol(cand)))
column_window <- 4:ncol(cand) # Update to match your range
variable <- ncol(cand) +1
correlation <- ncol(cand) + 2
# for (i in 1:length(cand$snp)) {
#   bar <- cand[i,]
#   cand[i,16] <- names(which.max(abs(bar[4:16]))) # gives the variable
#   cand[i,17] <- max(abs(bar[4:16]))              # gives the correlation
# }
# for (i in 1:length(cand$snp)) {
#   bar <- cand[i,]
#   cand[i,(variable)] <- names(which.max(abs(bar[4:(column_window)]))) # gives the variable
#   cand[i,(correlation)] <- max(abs(bar[4:(column_window)]))              # gives the correlation
# }
for (i in 1:nrow(cand)) {
  bar <- cand[i,]
  cand[i, variable] <- names(which.max(abs(bar[column_window]))) # gives the variable
  cand[i, correlation] <- max(abs(bar[column_window]))           # gives the correlation
}

colnames(cand)[variable] <- "predictor"
colnames(cand)[correlation] <- "correlation"

table(cand$predictor) #doesnt have to be all predictors-some arnet here because they arent "most corr var"

cand %>%
  group_by(predictor) %>%
  summarise(
    count = n(),
    average_correlation = mean(correlation, na.rm = TRUE)
  )


# save.image("code/R_code/backup_workspace_images/PC_CONDITION_catchment_vifsub10_RDA.RData") # As a backup :)
# load("code/R_code/backup_workspace_images/catchment_vifsub10_RDA.RData") #this is catchment and vif>10 predictors
# load("code/R_code/backup_workspace_images/PC_CONDITION_catchment_vifsub10_RDA.RData") #this is catchment and vif>10 predictors

  ############ ggplots after cand identified
#CAND IS THE LIST OF IMPORTANT SNPS
base_plot<-ggplot(rda_stats, aes(x = RDA1, y = RDA2)) +
  geom_point(color = "grey90", size = 1.4) +  # Plot all SNPs with normal scaling in grey
  labs(x =paste0("RDA1 (",rda1_prop, "%)"), y =paste0("RDA2 (", rda2_prop, "%)"), title = "RDA1 vs RDA2 for SNP Loadings") +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), linewidth = 0.6) +  # Horizontal line at y=0
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), linewidth = 0.6) +  # Vertical line at x=0
  geom_segment(data = env_scores_df, aes(x = 0, y = 0, xend = RDA1 * 7, yend = RDA2 * 7), 
               colour = "black", linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +  # Arrows from origin to scaled env scores
  geom_text_repel(data = env_scores_df, aes(x = 5 * RDA1, y = 5 * RDA2, label = str_wrap(name, width = 35)),
                  size=3,
                  min.segment.length = Inf,
                  position = position_nudge_center(0.2, 0.1, 0, 0),
                  max.overlaps = Inf, # Adjust as necessary for more flexibility
                  box.padding = 0.1, # Adjust padding around labels
                  point.padding = 0.1) +  # Adjust padding around points
  theme_classic() +
  theme(legend.position = "right",  # Position legend at the bottom
        legend.key.size = unit(0.5, "cm")) + # Adjust legend key size
  scale_x_continuous(expand = expansion(mult = 0.5))
base_plot
# just plot the outliers
unsc_cand_rda_stats <- rda_stats[rda_stats$SNP %in% cand$snp, ]
base_plot +  geom_point(data = unsc_cand_rda_stats, aes(color = "Candidate SNPs"), size = 3) +  # Overlay candidate SNPs
  scale_color_manual(name = "SNP Type", values = c("Candidate SNPs" = "red"))   # Colour legend

  
# Extract RDA statistics for 2 axes with special scaling = 3 for candidate SNPs (cand) *CARE* scaling must stay the same
cand_rda_stats <- rda_stats[rda_stats$SNP %in% cand$snp, ]
head(cand_rda_stats)
remove(combined_df)
combined_df <- left_join(cand_rda_stats, cand %>% dplyr::select(snp, predictor, correlation), by = c("SNP" = "snp"))
head(combined_df)

unique(combined_df$predictor)
ggplot(rda_stats_normal, aes(x = RDA1, y = RDA2)) +
  geom_point() +  # Plot all SNPs with normal scaling
  geom_point(data = combined_df, aes(color = predictor), size = 3) +  # Overlay candidate SNPs with special scaling
  # scale_color_manual(name = "SNP Type", values = c("Candidate SNPs" = "red")) +  # Colour legend
  labs(x = "RDA1", y = "RDA2", title = "RDA1 vs RDA2 for SNP Loadings") +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  theme_minimal()

# Getting better labels
head(combined_df)
combined_df$predictor <- sub("_lonlat$", "", combined_df$predictor) #make sure names match layers$layer_code
unique_predictors <- unique(combined_df$predictor)
columns_without_layer_code <- setdiff(unique_predictors, labels$layer_code)
# Add new labels to labels df based on ones missing
labels<-labels%>% dplyr::select(layer_code,name)
labels <- rbind(labels, data.frame(layer_code = columns_without_layer_code, name = columns_without_layer_code))
combined_df <- left_join(combined_df, labels, by = c("predictor" = "layer_code"))
head(combined_df)
# install.packages("stringr")
library(stringr)


base_plot+  geom_point(data = combined_df, aes(color = name), size = 3)+  # Overlay candidate SNPs with special scaling
  labs(colour="Predictor")


outlier_plot<-ggplot(rda_stats_normal, aes(x = RDA1, y = RDA2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Horizontal line at y=0
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +  # Vertical line at x=0
  geom_point(color="grey90") +  # Plot all SNPs with normal scaling
  geom_point(data = combined_df, aes(color = name), size = 3) +  # Overlay candidate SNPs with special scaling
  labs(x =paste0("RDA1 (",rda1_prop, "%)"), y =paste0("RDA2 (", rda2_prop, "%)"), title = "RDA1 vs RDA2 for SNP Loadings") +
  scale_color_discrete(name = "Predictor") +  # Adjust legend title if needed
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  geom_segment(data = env_scores_df, aes(xend = RDA1*6, yend = RDA2*6, x = 0, y = 0), colour = "black", linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_text_repel(data = env_scores_df, aes(x = 6 * RDA1, y = 6 * RDA2, label = str_wrap(name, width =20)), 
                  size=3.5,
                  min.segment.length = Inf,
                  position = position_nudge_center(0.2, 0.1, 0, 0),
                  max.overlaps = Inf, # Adjust as necessary for more flexibility
                  box.padding = 0.1, # Adjust padding around labels
                  point.padding = 0.1) +  # Adjust padding around points
  # geom_text(data = env_scores_df, aes(x = RDA1, y = RDA2, label = variable), 
  #           color = "blue", vjust = 1.5, hjust = 1.5)
  # theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11))
  scale_x_continuous(expand = expansion(mult = 0.5))+
  theme_classic() +
  theme(legend.position = "bottom",  # Position legend at the bottom
        legend.spacing.x = unit(0.1, "cm"),
        legend.key.size = unit(0.1, "cm"))  # Adjust legend key size
outlier_plot
# ggsave(filename = paste0("./output/figures/", run_name, " outlier rda.png"), plot = outlier_plot)
shell.exec(paste0(dir,"/output/figures/", run_name, " outlier rda.png"))
### My own manhatten plot
# First I need snp posititions- I have adjusted the rda extract function
# Filter rda_stats_normal based on cand1
rda_stats_cand1 <- rda_stats_normal %>%
  filter(SNP %in% cand1$snp)
ggplot(rda_stats_normal, aes(x = Position, y = abs(RDA1))) +
  geom_point(color = "grey80") +  # Plot all SNPs with normal scaling
  geom_point(data = rda_stats_cand1, colour = "red", size = 3) +  # Overlay candidate SNPs with special scaling
  labs(x = "Position in Genome", y = paste0("RDA1 (",rda1_prop, "%)"), title = "RDA1 SNP Loadings and Position in Genome") +
  scale_color_discrete(name = "Predictor") +  # Adjust legend title if needed
  theme_classic() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  # Ensure x-axis starts at 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))  # Ensure y-axis starts at 0
 
rda_stats_cand2 <- rda_stats_normal %>%
  filter(SNP %in% cand2$snp)
ggplot(rda_stats_normal, aes(x = Position, y = abs(RDA1))) +
  geom_point(color = "grey80") +  # Plot all SNPs with normal scaling
  geom_point(data = rda_stats_cand2, colour = "red", size = 3) +  # Overlay candidate SNPs with special scaling
  labs(x = "Position in Genome", y = paste0("RDA1 (",rda2_prop, "%)"), title = "RDA1 SNP Loadings and Position in Genome") +
  scale_color_discrete(name = "Predictor") +  # Adjust legend title if needed
  theme_classic() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  # Ensure x-axis starts at 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))  # Ensure y-axis starts at 0

 
ggplot(rda_stats_normal, aes(x = Position, y = -log10(p_values))) +
  geom_point(color="grey80") +  # Plot all SNPs with normal scaling
  geom_point(data = combined_df, color = "red", size = 3) +  # Overlay candidate SNPs with special scaling
  labs(x = "Position in Genome", y = "-log10(p.values)", title = "P-values from Capblanq") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  scale_color_discrete(name = "Predictor") +  # Adjust legend title if needed
  # theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
  theme_classic()+
  theme(legend.position = "right",  # Position legend at the bottom
        legend.key.size = unit(0.5, "cm"))  # Adjust legend key size

# P-values threshold after Bonferroni correction
length(rda_stats_normal$p_values)
thres_env <- 0.01/length(rda_stats_normal$p_values)
rda_stats_cand <- rda_stats_normal %>%
  filter(SNP %in% cand$snp)
ggplot() +
  geom_point(data = rda_stats_normal, aes(x=Position, y=-log10(p_values)), size=1.4) +
  geom_point(data = rda_stats_cand, aes(x=Position, y=-log10(p_values)), size=1.4, colour="red") + #weird because cand isnt here made by p-values, but on sd from rda axis
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11) +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

head(map)
head(cand)
# Remove suffixes from SNP names in cand
cand_clean <- cand %>%
  mutate(snp_clean = sub("_[A-Z]$", "", snp))

# Filter the map file based on the cleaned SNP names in cand
cand_contigs <- map %>%
  filter(SNP_ID %in% cand_clean$snp_clean)

# Count the number of candidate SNPs in each chromosome
cand_contig_dist<-cand_contigs %>%
  group_by(Chromosome) %>%
  summarise(count = n())
nrow(cand_contig_dist)


# General map file contig number
contig_dist<-map %>%
  group_by(Chromosome) %>%
  summarise(count = n())
nrow(contig_dist)
####################      More code         ##########################################
# capblanq tutorial
#### Function to conduct a RDA based genome scan, # needs packages robust and qvalue
rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

## Running the function 
rdadapt_env<-rdadapt(pearl.rda, 6)
summary(rdadapt_env)


# P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_env$p.values)

## Identifying the loci that are below the p-value threshold
outliers <- data.frame(Loci = colnames(gen.imp)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(gen.imp)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))

## Top hit outlier per contig
outliers <- outliers[order(outliers$contig, outliers$p.value),]

## List of outlier names
outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$contig)])

## Formatting table for ggplot
locus_scores <- scores(pearl.rda, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(pearl.rda, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
TAB_var_lab <- TAB_var %>% 
  tibble::rownames_to_column(var = "variable") %>% left_join(labels, by = c("variable" = "layer_code"))

ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  # geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  geom_text_repel(data = TAB_var_lab, aes(x = 1.1 * RDA1, y = 1.1 * RDA2, label = str_wrap(name, width = 35)),
                  size=2.5,
                  min.segment.length = Inf,
                  # position = position_nudge_center(0.2, 0.1, 0, 0),
                  max.overlaps = Inf, # Adjust as necessary for more flexibility
                  box.padding = 0.1, # Adjust padding around labels
                  point.padding = 0.1) +  # Adjust padding around points
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

## Manhattan plot
Outliers <- rep("Neutral", length(colnames(gen.imp)))
Outliers[colnames(gen.imp)%in%outliers$Loci] <- "All outliers"
Outliers[colnames(gen.imp)%in%outliers_rdadapt_env] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_manhatan <- data.frame(pos = 1:length(colnames(gen.imp)), 
                           pvalues = rdadapt_env$p.values, 
                           Outliers = Outliers)
TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Outliers),]
ggplot(data = TAB_manhatan) +
  geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
  scale_color_manual(values = c("gray80", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  # facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))


# Saving as CSV
# write.csv(outliers_rdadapt_env, file = "data/capblanq_RDA_outliers_rdadapt_env.csv", row.names = FALSE)
# write.csv(cand, file = "data/forester_RDA_cand.csv", row.names = FALSE)
cand_PC<-read.csv("data/forester_PCCondition_RDA_cand.csv")
outliers_rdadapt_env_PC<-read.csv("data/capblanq_PCCondition_RDA_outliers_rdadapt_env.csv")
cand_o<-read.csv("data/forester_RDA_cand.csv")
outliers_rdadapt_env_o<-read.csv("data/capblanq_RDA_outliers_rdadapt_env.csv")

## Overlap in cand/outliers using sd and p-val methods
overlap <- intersect(outliers_rdadapt_env, cand$snp)
overlap <- intersect(PC_outliers_rdadapt_env$x, PC_cand$snp)

# Percent overlap for each method
percent_overlap_rdadapt <- (length(overlap) / length(outliers_rdadapt_env)) * 100
percent_overlap_zscore <- (length(overlap) / length(cand$snp)) * 100

# Output results
cat("Overlap in rdadapt:", percent_overlap_rdadapt, "%\n")
cat("Overlap in Z-score:", percent_overlap_zscore, "%\n")
cat("Total overlapping SNPs:", length(overlap), "\n")
####################################

# Katie E. Lotterhos 
rda_trait_pred <- function(rdaobj, env_row, K){
  #rdaobj is RDA object
  #env_row is the row of the environment in the biplot output
  #K is the number of RDA axes
  scores <- scores(rdaobj, choices=1:K)
  ind.sc <- scores$sites
  pred <- matrix(NA, nrow=nrow(ind.sc), ncol=K)
  for (k in 1:K){
    pred[,k] <- ind.sc[,k]*eigenvals(rdaobj)[k]*summary(rdaobj)$biplot[env_row,k]
  }
  trait_pred <- scale(rowSums(pred))
  return(trait_pred) 
}
pearl.rda
pearl.rda$CCA$biplot
FW_soil_wavg_07_traitpredict <- rda_trait_pred(pearl.rda, 1, 4)
plot(scale(env$FW_soil_wavg_07), FW_soil_wavg_07_traitpredict, xlab="Evolved trait value in simulations",
     ylab="RDA trait prediction")
abline(0,1)
###################################



##### Plot RDA loadings against outlier positions #####################

snp_positions<-data.frame(SNPname= locNames(genRAW), 
                          Position =genRAW@position)

outlier_positions <- snp_positions$Position[match(outlier_SNPs, snp_positions$SNP)]

cand1$position <- snp_positions$Position[match(cand1$snp, snp_positions$SNP)]
cand1$squared_RDA1 <- cand1$loading^2

plot(cand1$position, cand1$loading)
plot(cand1$position, cand1$squared_RDA1)
# Extract SNP names from row names of load.rda
snp_names <- rownames(load.rda)

# Create a data frame with SNP names, positions, and loadings
loadings_data <- data.frame(
  SNP = snp_names,
  Position = snp_positions$Position[match(snp_names, snp_positions$SNPname)],
  RDA1 = load.rda[, "RDA1"],
  RDA2 = load.rda[, "RDA2"],
  RDA3 = load.rda[, "RDA3"],
  RDA4 = load.rda[, "RDA4"],
  RDA5 = load.rda[, "RDA5"],
  RDA6 = load.rda[, "RDA6"]
)
loadings_data$squared_RDA1 <- loadings_data$RDA1^2
plot(loadings_data$Position, loadings_data$squared_RDA1, log="")
points(cand1$position, cand1$squared_RDA1, col="red" )

#################################################
# Manhattan plot
install.packages("qqman")
library(qqman)
library(qvalue)
manhattan(pearl.rda, chr="Chr",bp="Pos", p="Pval", snp="SNP", main="Fv'/Fm' Control at 3 d.a.s.", ylim=c(0,12), col = c("turquoise3", "purple4"), suggestiveline=FALSE, genomewideline=-log10(b))

# Prepare table for graphing:
args <- commandArgs(trailingOnly=TRUE)
results$Pval <-as.numeric(as.character(results$Pval))


# 1. Extract RDA loadings for the outlier SNPs
outlier_SNPs <- unique(c(cand1$snp, cand2$snp, cand3$snp, cand4$snp, cand5$snp, cand6$snp))

outlier_loadings <- data.frame(SNP = outlier_SNPs,
                               RDA1 = rep(NA, length(outlier_SNPs)),
                               RDA2 = rep(NA, length(outlier_SNPs)),
                               RDA3 = rep(NA, length(outlier_SNPs)),
                               RDA4 = rep(NA, length(outlier_SNPs)),
                               RDA5 = rep(NA, length(outlier_SNPs)),
                               RDA6 = rep(NA, length(outlier_SNPs)))

for (i in 1:length(outlier_SNPs)) {
  SNP <- outlier_SNPs[i]
  outlier_loadings[i, 2:7] <- load.rda[SNP, ]
}

# 2. Extract positions for the outlier SNPs
outlier_indices <- match(outlier_SNPs, rownames(genRAW@genotypes))
outlier_positions <- genRAW@position[outlier_indices]
# 3. Plot RDA loadings against locus positions
plot(outlier_positions, outlier_loadings$RDA1, 
     xlab = "Locus Position", ylab = "RDA Loading", 
     main = "RDA Loadings against Locus Positions (Axis 1)")
points(outlier_positions, outlier_loadings$RDA2, col = "red")
points(outlier_positions, outlier_loadings$RDA3, col = "blue")
legend("topright", legend = c("RDA1", "RDA2", "RDA3"), col = c("black", "red", "blue"), pch = 1)



























head(pred)
head(gen.imp)
# # pred
# #### PCA to flatten bioclim variable
# pc <- prcomp(envscaled,
#              center = TRUE,
#              scale. = FALSE)
# attributes(pc)
# pc$center
# pc$scale
# print(pc)
# summary(pc)
# pairs.panels(pc$x, scale=T)
# 
# pc_scores <- predict(pc, newdata = envscaled)[, 1:3]  # Extracting first two principal components
# pca<-as.data.frame(pc_scores)
# pearl.rda <- rda(gen.imp ~ ., data=pca, scale=T)
# pearl.rda
# RsquareAdj(pearl.rda)
# summary(eigenvals(pearl.rda, model = "constrained"))
# screeplot(pearl.rda)
# 
# signif.full <- anova.cca(pearl.rda, parallel=getOption("mc.cores")) # default is permutation=999
# signif.full #takes time
# # 
# vif.cca(pearl.rda) #### most are above 10 so multicoliniarity is a concern....
# 
# plot(pearl.rda, scaling=3)          # default is axes 1 and 2
#  plot(pearl.rda, choices = c(1, 3), scaling=3)  # axes 1 and 3
#  load.rda <- scores(pearl.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
#  hist(load.rda[,1], main="Loadings on RDA1")
#  hist(load.rda[,2], main="Loadings on RDA2")
#  hist(load.rda[,3], main="Loadings on RDA3") 
#  outliers <- function(x,z){
#      lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
#      x[x < lims[1] | x > lims[2]]               # locus names in these tails
#    }
#  cand1 <- outliers(load.rda[,1],3) # 38
#  cand2 <- outliers(load.rda[,2],3) # 69
#  cand3 <- outliers(load.rda[,3],3) # 34
#  ncand <- length(cand1) + length(cand2) + length(cand3)
#  ncand
# 
#  cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
#  cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
#  cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
#  colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
#  cand <- rbind(cand1, cand2, cand3)
#  cand$snp <- as.character(cand$snp)
#  foo <- matrix(nrow=(ncand), ncol=3)  # 2 columns for 2 predictors
# colnames(foo)<-colnames(pca)
# length(cand$snp)
# for (i in 1:length(cand$snp)) {
#   nam <- cand[i,2]
#   snp.gen <- gen.imp[,nam]
#   foo[i,] <- apply(pca,2,function(x) cor(x,snp.gen))
# }
# 
# cand <- cbind.data.frame(cand,foo)  
# head(cand)
# 
# # Investigate the candidates
# 
# length(cand$snp[duplicated(cand$snp)])  # 89 duplicate detections
# 
# foo <- cbind(cand$axis, duplicated(cand$snp)) 
# table(foo[foo[,1]==1,2]) # 27 duplicates on axis 1
# table(foo[foo[,1]==2,2]) #  7 duplicates on axis 2
# table(foo[foo[,1]==3,2]) # 55 duplicates on axis 3
# cand <- cand[!duplicated(cand$snp),] # remove duplicate detections
# 
# 
# for (i in 1:length(cand$snp)) {
#   bar <- cand[i,]
#   cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable
#   cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation
# }
# 
# 
# colnames(cand)[7] <- "predictor"
# colnames(cand)[8] <- "correlation"

#A more complex approach to a colour scheme #####################################################
# Calculate distances
install.packages("geosphere")
library(geosphere)
dist_matrix <- distm(env_unique[, c("Longitude", "Latitude")])

# Scale distances to range [0, 1]
scaled_dist <- (dist_matrix - min(dist_matrix)) / (max(dist_matrix) - min(dist_matrix))*100

# Generate colours using viridis colour scheme
num_colors <- length(eco)
color_palette <- turbo(num_colors)

# Assign colours to locations based on scaled distances
location_colors <- apply(scaled_dist, 1, function(x) color_palette[which.max(x)])
swatch(location_colors)

# Generate a colour palette from viridis
library(viridis)
palette <- turbo(100)  # You can adjust the number to change the resolution of colours
install.packages("hues")
library(hues)
swatch(palette)

scaled_dist <- (dist_matrix - min(dist_matrix)) / (max(dist_matrix) - min(dist_matrix))*100
extracted_hex <- palette[scaled_dist]
swatch(extracted_hex)

# Interpolate the hex codes based on the distribution of x along a 0 to 1 scale
bg <- palette[findInterval(mean_distances, seq(0, 1, length.out = length(palette)))]
swatch(bg)



#########################
# Remove duplicate entries in env data frame
env_unique <- env[!duplicated(env$Population), ]
dim(env_unique)
# Calculate distances
dist_matrix <- distm(env_unique[, c("Longitude", "Latitude")])
dim(dist_matrix)
# Calculate mean distances for each population
mean_distances <- tapply(dist_matrix, INDEX = factor(env_unique$Population), FUN = mean)

# Scale mean distances to range [0, 1]
scaled_distances <- 1 - (mean_distances - min(mean_distances)) / (max(mean_distances) - min(mean_distances))

# Generate colours using viridis colour scheme
num_colors <- length(unique(env_unique$Population))
color_palette <- viridis(num_colors)

# Assign colours to populations based on scaled mean distances
population_colors <- color_palette[findInterval(scaled_distances, seq(0, 1, length.out = num_colors))]

# Plot
plot(env$Longitude, env$Latitude, col = population_colors[match(env$Population, names(population_colors))], pch = 16, cex = 2, main = "Population Colour Palette", xlab = "Longitude", ylab = "Latitude")
legend("topright", legend = unique(env_unique$Population), col = population_colors, pch = 16, cex = 1)

#########################
# Calculate distances
env_unique <- env[!duplicated(env$Population), ]
dist_matrix <- distm(env_unique[, c("Longitude", "Latitude")])

# Scale distances to range [0, 1]
scaled_distances <- 1 - (dist_matrix - min(dist_matrix)) / (max(dist_matrix) - min(dist_matrix))

# Generate colours using viridis colour scheme with 100 colours
num_colors <- 100
color_palette <- viridis(num_colors)

# Extract subset of colours based on scaled distances
population_colors <- color_palette[findInterval(scaled_distances, seq(0, 1, length.out = num_colors))]
swatch(population_colors)

#########################
###########  Capblancq way of plotting
## Formatting table for ggplot
locus_scores <- scores(pearl.rda, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(pearl.rda, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

## Manhattan plot
Outliers <- rep("Neutral", length(colnames(pred)))
Outliers[colnames(pred)%in%outliers$Loci] <- "All outliers"
Outliers[colnames(pred)%in%outliers_rdadapt_env] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_manhatan <- data.frame(pos = 1:length(colnames(pred)), 
                           pvalues = rdadapt_env$p.values, 
                           Outliers = Outliers)
TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Outliers),]
ggplot(data = TAB_manhatan) +
  geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
###########

# # Function to calculate mean distance for each population
# calculate_distance_matrix <- function(lat_long_data) {
#   num_points <- nrow(lat_long_data)
#   distance_matrix <- matrix(NA, nrow = num_points, ncol = num_points)
#   for (i in 1:num_points) {
#     for (j in 1:num_points) {
#       distance_matrix[i, j] <- sqrt((lat_long_data[i, "Latitude"] - lat_long_data[j, "Latitude"])^2 +
#                                       (lat_long_data[i, "Longitude"] - lat_long_data[j, "Longitude"])^2)
#     }
#   }
#   return(distance_matrix)
# }
# # Calculate distance matrix
# distance_matrix <- calculate_distance_matrix(env)
# 
# # Flatten the distance matrix to single numbers for each population
# mean_distances <- calculate_mean_distances(distance_matrix, unique(env$Population))
# plot(mean_distances)
# # Function to normalize distances
# normalize_distances <- function(distances) {
#   min_dist <- min(distances)
#   max_dist <- max(distances)
#   normalized_distances <- (distances - min_dist) / (max_dist - min_dist)
#   return(normalized_distances)
# }
# normalized_distances <- normalize_distances(distance_matrix)
# mean_distances <- calculate_mean_distances(normalized_distances, unique(env$Population))
# plot(mean_distances)








# Useful code for plotting pca biplot
install.packages("devtools")
library(devtools)
#install_github("vqv/ggbiplot")
install.packages("ggbiplot")
library(ggbiplot)
g <- ggbiplot(pc,
              obs.scale = 1,
              var.scale = 1,
              groups = env$Population,
              ellipse = TRUE,
              circle = TRUE,
              ellipse.prob = 0.68)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)


##### Tutorial base R code
# 
# env$Population<-as.factor(env$Population)
# levels(env$Population)
# eco<-env$Population
# levels(eco)
# # bg <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31") #rendom colours
# # library(viridisLite)
# # bg<-turbo(22)
# bg <- c("#008080", "#CC00CC", "#99FF66", "#FF5050", "#FF0000", "#CC66FF", "#0066FF", "#66CCFF", "#00CC99", "#CCCC00", "#B8B400", "#6666FF", "#FF9999", "#00FF99", "#00FF00", "#00CC00", "#FF3399", "#CC3399", "#FF99CC", "#CC99FF", "#CC3300", "#99CCFF") #eyeballed distance
# 
# # 
# # bg <- c("#008b8b",#lax ###vickys colours
# #         "#00cdcd",#lax_un
# #         "#6959cd", #clas
# #         "#4682b4", #bad
# #         "#00bfff",#kir
# #         "#8ee5ee",#rhe
# #         "#27408b", #poll
# #         "#4169e1",#pollst
# #         "#8b0000",#kerry
# #         "#ff4500",#tirr
# #         "#fa8072", #cono
# #         "#d02090",#moris
# #         "#daa520",#hor
# #         "#ffd700",#kilm
# #         "#b3ee3a", #moid
# #         "#32cd32", #bhri
# #         "#b4eeb4", #grig
# #         "#228b22", #struth
# #         "#DB3A07FF", 
# #         "#C22402FF", 
# #         "#A11201FF",
# #         "#7A0403FF")
# # levels(eco)
# # bg <- c("#30123BFF", "#3D358CFF", "#4457C8FF", "#4777EFFF", "#4195FFFF", "#2EB3F3FF", "#1BD0D5FF", "#1AE4B6FF", "#35F393FF", "#62FC6BFF", "#90FF48FF", "#B3F836FF", "#D2E935FF", "#EBD339FF", "#FABA39FF", "#FE9B2DFF", "#F9771EFF", "#EC5410FF", "#DB3A07FF", "#C22402FF", "#A11201FF", "#7A0403FF")
# 
# 
# 
# # axes 1 & 2
# plot(pearl.rda, type="n", scaling=3)
# points(pearl.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
# points(pearl.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the wolves
# text(pearl.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
# legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# 
# # axes 1 & 3
# plot(pearl.rda, type="n", scaling=3, choices=c(1,3))
# points(pearl.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
# points(pearl.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco], choices=c(1,3))
# text(pearl.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
# legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# 
# # Axes 3 vs 4
# plot(pearl.rda, type="n", scaling=3, choices=c(3,4))
# points(pearl.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(3,4))
# points(pearl.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco], choices=c(3,4))
# text(pearl.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(3,4))
# legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# 
# # Axes 4 vs 5
# plot(pearl.rda, type="n", scaling=3, choices=c(4,5))
# points(pearl.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(4,5))
# points(pearl.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco], choices=c(4,5))
# text(pearl.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(4,5))
# legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# 
# # Axes 5 vs 6
# plot(pearl.rda, type="n", scaling=3, choices=c(5,6))
# points(pearl.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(5,6))
# points(pearl.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco], choices=c(5,6))
# text(pearl.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(5,6))
# legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
# 
# 
# 
######## Base R code- after CAND SNPs identified
# unique_predictors<-unique(cand$predictor)
# # Function to generate colours
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# # Generate colours using gg_color_hue function
# cols = gg_color_hue(length(unique(cand$predictor)))
# # Create predictor_colours as a named vector
# predictor_colours <- setNames(cols, unique_predictors)
# names(predictor_colours)
# # Map predictor names to colours in 'tenv'
# tenv <- cand$predictor
# tenv <- predictor_colours[tenv]
# # tenv[tenv=="FW_lc_wavg_03_lonlat"] <- 'deeppink'
# # tenv[tenv=="WC_bio17_lonlat"] <- 'turquoise'
# # tenv[tenv=="WC_bio2_lonlat"] <- 'purple'
# # tenv[tenv=="WC_bio3_lonlat"] <- 'pink'
# # tenv[tenv=="WC_bio5_lonlat"] <- 'yellow'
# # tenv[tenv=="WC_bio9_lonlat"] <- 'green'
# 
# # color by predictor:
# sel <- cand$snp
# col.pred <- rownames(pearl.rda$CCA$v) # pull the SNP names
# head(col.pred)
# 
# # Colouring snps based on pred
# for (i in 1:length(sel)) {           # color code candidate SNPs
#   foo <- match(sel[i],col.pred)
#   col.pred[foo] <- tenv[i]
# }
# 
# col.pred[grep("scaf",col.pred)] <- '#f1eef6' # non-candidate SNPs ****MUST CHANGE "SNP" TO A COMMON STRING OF ALL SNPS
# empty <- col.pred
# empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
# empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
# ou <- cols
# 
# 
# # axes 1 & 2
# plot(pearl.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
# points(pearl.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
# points(pearl.rda, display="species", pch=21, cex=1.5, col=empty.outline, bg=empty, scaling=3)
# text(pearl.rda, scaling=3, display="bp", col="#0868ac", cex=1)
# legend("bottomright", legend=unique_predictors, bty="n", col="gray32", pch=21, cex=1, pt.bg=ou)
# # axes 1 & 3
# plot(pearl.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
# points(pearl.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
# points(pearl.rda, display="species", pch=21, cex=1.5, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
# text(pearl.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
# legend("bottomright", legend=c("bio10Mean_T_Warm_Quarter", "bio17precip_driestq", "bio2_MDiurnalRange", "bio_3Isotherm", "bio5_MaxT_WarmMonth", "bio9_Mean_T_Dry_Quarter"), bty="n", col="gray32", pch=21, cex=1.5, pt.bg=ou)
# 
# # axes 2 & 3
# plot(pearl.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3))
# points(pearl.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(2,3))
# points(pearl.rda, display="species", pch=21, cex=1.5, col=empty.outline, bg=empty, scaling=3, choices=c(2,3))
# text(pearl.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(2,3))
# legend("bottomright", legend=c("bio10Mean_T_Warm_Quarter", "bio17precip_driestq", "bio2_MDiurnalRange", "bio_3Isotherm", "bio5_MaxT_WarmMonth", "bio9_Mean_T_Dry_Quarter"), bty="n", col="gray32", pch=21, cex=1.5, pt.bg=ou)
# 
# # axes 3 & 4
# plot(pearl.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(3,4))
# points(pearl.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(3,4))
# points(pearl.rda, display="species", pch=21, cex=1.5, col=empty.outline, bg=empty, scaling=3, choices=c(3,4))
# text(pearl.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(3,4))
# legend("bottomright", legend=c("bio10Mean_T_Warm_Quarter", "bio17precip_driestq", "bio2_MDiurnalRange", "bio_3Isotherm", "bio5_MaxT_WarmMonth", "bio9_Mean_T_Dry_Quarter"), bty="n", col="gray32", pch=21, cex=1.5, pt.bg=ou)
# 
# # axes 5 & 6
# plot(pearl.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(5,6))
# points(pearl.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(5,6))
# points(pearl.rda, display="species", pch=21, cex=1.5, col=empty.outline, bg=empty, scaling=3, choices=c(5,6))
# text(pearl.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(5,6))
# legend("bottomright", legend=c("bio10Mean_T_Warm_Quarter", "bio17precip_driestq", "bio2_MDiurnalRange", "bio_3Isotherm", "bio5_MaxT_WarmMonth", "bio9_Mean_T_Dry_Quarter"), bty="n", col="gray32", pch=21, cex=1.5, pt.bg=ou)
# ###### Base R code totorial
# 
















