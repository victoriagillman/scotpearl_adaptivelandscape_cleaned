### General population statistics using VP results and editing from Sam Beck code

# Contains
  # Basic descriptive data 
  # PCA analysis using pcadpat
  # Admixture analysis using LEA and snmf
  # Creating admixture pie charts and plotting on a map

# Housekeeping!
rm(list=ls())

# install.packages("adegenet")
# install.packages("snpStats")
# library(adegenet)   # Used for genetic data preparation
# library(snpStats)

#Set wd to be where .bed .bim .fam files are
setwd("./data/plink_VP/")
getwd()


###====================================================================================

library(dplyr)      # For data manipulation
library(data.table) # For reading in gen files
library(pcadapt)    # For PCA data
# install.packages("devtools")
# library("devtools")
# install_github("jdstorey/qvalue")
library(qvalue) #For getting pcadapt scores
# install.packages("tidyverse")
library(tidyverse) #for how many snps do we have bit

#===========================

# Match metadata to individual IDs in order of genotype file
meta <- fread("PM_Dec4.fam") %>% as.data.frame() %>% mutate(SiteCode = V1, ID = V2) %>%
  dplyr::select(SiteCode, ID) %>% mutate(ID = ifelse(ID == "ICY837","ICY837notICY773", ID))
# vicec <- read.csv("../../data/confidential/swab_coord_vc_anon_names.csv")  %>% mutate(Sample_ID = ifelse(Sample_ID == "ICY837","ICY837notICY773", Sample_ID))
vicec<-read.csv("../../data/confidential/trimmed_env_data.csv")%>%
  mutate(Sample_ID = ifelse(Sample_ID == "ICY837","ICY837notICY773", Sample_ID)) %>% mutate(Population = ifelse(Population == "Conon_Gharbhrain_wrong_sample_ID", "Conon_Glascarnoch", Population))%>% mutate(AnonSiteCode = ifelse(AnonSiteCode == "ER2", "ER1", AnonSiteCode))


# Find mismatches between datasets
meta %>% filter(ID %in% setdiff(ID, vicec$Sample_ID)) #Prints those in meta not in vc
vicec %>% filter(Sample_ID %in% setdiff(Sample_ID, meta$ID)) %>%dplyr::select(Sample_ID, Population) #Those in vc not in meta- eg should be those filtered out in plink etc


# Join genetic metadata to main metadata csv
meta.fam <- left_join(meta, vicec, by = c("ID" = "Sample_ID"))
write.csv(meta.fam, "../../data/confidential/meta_for_filtered_genomes.csv")


# Check how even sample sizes are : sample sizes vary between 2-15 
nrow(meta.fam) #156 individuals
sample_size <- meta.fam %>% group_by(AnonSiteCode) %>% summarise(n = n()) %>% print()
min(sample_size$n) #2
max(sample_size$n) #15
vc_size <- meta.fam %>% group_by(VCNAME) %>% summarise(n = n()) %>% print()


# SNP maps of both full and LD datasets
chrom_snp_map <- fread("PM_Dec4.bim") %>%  
  mutate(CHROM =V1, SNP = V2, BP = V4) %>%  
  dplyr::select(CHROM, BP, SNP)
chrom_snp_map_LDpruned <- fread("LDpruned.bim") %>%  
  mutate(CHROM =V1, SNP = V2, BP = V4) %>%  
  dplyr::select(CHROM, BP, SNP)


# PCA functions ------

PCADAPT_import_PCA <- function(filename, Knum) {
  to_pc <- read.pcadapt(paste0(filename, ".bed"), type = "bed")
  PCA <- pcadapt(to_pc, K = Knum, min.maf = 0.01)
  return(PCA)
} # a function to import data from .bed file and run a PCA with specific number of components and minimum maf of 0.01

PCADAPT_scores_meta <- function(PCAobj, filename, Knum, meta) {
  FAM <- fread(paste0(filename, ".fam")) %>% 
    dplyr::select(V1,V2) %>% 
    mutate(SiteCode = V1, ID = V2) %>%  
    dplyr::select(-V1, -V2)
  meta_ordered <- inner_join(FAM, meta)
  PCAscores <- data.frame(PCAobj$scores)
  colnames(PCAscores) <- paste0("PC", rep(1:Knum))
  PCA_meta <- bind_cols(meta_ordered, PCAscores)
  return(PCA_meta)
} ## a function to match and merge PCA scores with metadata- read .fam file, aligns IDs and binds Knum PCA scores to metadata

PCADAPT_var_per_axis <- function(PCAobj) {
  var_per_axis <- (PCAobj$singular.values^2)
  return(var_per_axis)} #a function to get variance explained by each PCA axis by squaring the singular values from the PCA object

PCADAPT_scores_qval_map  <- function(PCAobj, mapobj, Knum) {
  qvalues <- qvalue(PCAobj$pvalues)
  PCAobj_map_pval_qval <- bind_cols(chrom_snp_map_LDpruned, PCAobj$pvalues, qvalues$qvalues)
  colnames(PCAobj_map_pval_qval)[4:5] <- c("pvalues", "qvalues")
  loadings <- data.frame(PCAobj$loadings)
  names(loadings) <- paste0("PC", rep(1:Knum), "_loading")
  PCAobj_map_pval_qval_loadings <- bind_cols(PCAobj_map_pval_qval, loadings)
  return( PCAobj_map_pval_qval_loadings)} #a function to get PCA pvalues, qvalues, loadings and SNP locations

#--------

# Run full PCA and chose how many axis going forward for analysis 
tot_knum<- as.numeric(length(unique(meta.fam$SiteCode))) # 18 is max number of pops
PC_K18 <- PCADAPT_import_PCA(filename = "LDpruned", Knum = tot_knum) # Run full PCA with knum equal to number pops
#screeplot
PC_K18
scree <- plot(PC_K18, option = "screeplot") + theme_classic()  #plateaus at ~K = 4  and possible 10?
#ggsave("figures/scree_K18.png") # K = 3
variance_threshold <- 0.50  # Set desired explained variance for used Knum (e.g., 50%, 90%)
chosen_K <- which(cumsum((PC_K18$singular.values^2) / sum(PC_K18$singular.values^2)) >= variance_threshold)[1]
chosen_K = 4

# Run chosen PCA
PC_K4 <- PCADAPT_import_PCA(filename = "LDpruned", Knum = chosen_K)
summary(PC_K4)
# Now to pull the PCA scores of chose PCA, plus other metrics and bind with map file
# if error with below, re-run chrom_snp_map without using chrom names
PC_K4_scores <- PCADAPT_scores_qval_map(PCAobj = PC_K4, mapobj = chrom_snp_map_LDpruned, Knum = 4)

##Population structure analysis with PCA - bind PCA scores with metadata
PC_K4_meta <-  PCADAPT_scores_meta(PCAobj = PC_K4, filename =  "LDpruned", Knum = 4, meta = meta.fam)
PC_K4_meta <- PC_K4_meta %>%
  dplyr::rename(PC1_genetic = PC1,
         PC2_genetic = PC2,
         PC3_genetic = PC3,
         PC4_genetic = PC4) #PC_genetic is new name
colnames(PC_K4_meta)

# export metadata + PCA scores 
write.csv(PC_K4_meta, "../../data/confidential/PCA_values_for_indi_K4.csv")



###### Plots ######
library(ggforce)# for ellipeses
library(stringr) ##to wrap the legend labels

### PCA of genetic variation: 
labels<-read.csv("../../data/too_large_files/all_freshwater_sdmpredictors/layerschoice_metadata.csv")
legend_labels <- setNames(labels$name, labels$layer_code)

percentage_var <- (PC_K18$singular.values)^2*100 # 10.393347  6.560905  4.909740  3.501487
percentage_var
base_plot <- ggplot(data = PC_K4_meta, aes(x = PC1_genetic, y = PC2_genetic)) +
  labs(x = sprintf("PC1 (%.2f%%)", percentage_var[1]), 
       y = sprintf("PC2 (%.2f%%)", percentage_var[2]),
       fill = "Vice County Anon \nPopulations"
  ) +
  theme_classic() +
  theme(text = element_text(size = 18))
base_plot


# PC1 vs PC2
just_pc1v2 <- base_plot +
  geom_point(aes(fill = AnonSiteCode), size = 5, shape = 21, colour = "black", stroke = 0.5) +
  labs(
    title = "PCA of filtered genetic variation",
    subtitle = "Coloured by Population, 156 individuals, 3456 SNPs",
  )
just_pc1v2


# Ellipse around vice county and pop
ellipsebyvc<- base_plot +
   geom_point(aes(fill = AnonSiteCode), size = 5, shape = 21, colour = "black", stroke = 0.5) +
   labs(
     title = "PCA of filtered genetic variation",
     subtitle = "Coloured by population, 156 individuals, 3456 SNPs",
   ) +
   geom_mark_ellipse(aes(group = VCNAME, label = VCNAME), fill = "darkblue", alpha = 0.15, expand = unit(5, "mm"))
ellipsebyvc

ellipsebypop_anon<- base_plot +
  geom_point(aes(fill = AnonSiteCode), size = 5, shape = 21, colour = "black", stroke = 0.5) +
  labs(
    title = "PCA of filtered genetic variation",
    subtitle = "Coloured by Population, 156 individuals, 3456 SNPs",
  ) +
  geom_mark_ellipse(aes(group = AnonSiteCode, fill = AnonSiteCode, label = AnonSiteCode), expand = unit(5, "mm"))+
  theme(legend.position = "none")
ellipsebypop_anon

ellipsebypop<-base_plot +
  geom_point(aes(fill = SiteCode), size = 5, shape = 21, colour = "black", stroke = 0.5) +
  labs(
    title = "PCA of filtered genetic variation",
    subtitle = "Coloured by Population, 156 individuals, 3456 SNPs",
    fill = "Vice County \nPopulations"
  ) +
  geom_mark_ellipse(aes(group = SiteCode, fill = SiteCode, label = SiteCode), expand = unit(5, "mm"))
ellipsebypop
ggsave(plot= ellipsebypop_anon, 
       width = 1000/100,
       height = 900/100,
       units = c( "in"),
       dpi = 600,
       filename = "../../output/figures/pca_ld.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/figures/pca_ld.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/pca_ld.png")


# Colour by env variation
host<-base_plot + 
  geom_point(aes(fill = Host), size = 5, shape = 21, colour = "black", stroke = 0.5)+
  labs(fill="Host")
host
cation<-base_plot + 
  geom_point(aes(fill = FW_soil_wavg_07), size = 5, shape = 21, colour = "black", stroke = 0.5)+
  labs(fill= str_wrap(paste0(labels$name[labels$layer_code == "FW_soil_wavg_07"]), width=20))
wild<-base_plot + 
  geom_point(aes(fill = wildness), size = 5, shape = 21, colour = "black", stroke = 0.5)+
  labs(fill= str_wrap(paste0(labels$name[labels$layer_code == "wildness"]), width=20))
bio10<-base_plot + 
  geom_point(aes(fill = FW_hydro_wavg_10), size = 5, shape = 21, colour = "black", stroke = 0.5)+
  labs(fill= str_wrap(paste0(labels$name[labels$layer_code == "FW_hydro_wavg_10"]), width=20))

library(patchwork)
host+cation+wild+bio10+
  plot_layout(nrow = 2, byrow = 2) +
  plot_annotation(title = "PCA of filtered genetic variation", subtitle = "156 individuals, 3456 LD-trimmed SNPS", 
                  theme=theme(text = element_text(size = 18)), tag_levels = 'a', tag_suffix=")")
      
      
length(unique(PC_K4_meta$VCNAME))

# #dev.off()
######

#=========================
#######################################
# snmf has been run and output is saved as qlong (proceed from qlong onwards)
### code bvelow only needs to be run if needing new qlong file!

##LEA snmf
###SNMF #ALL  ### # only need to do ped2lfmm once
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("LEA")
library(LEA)
# ped2lfmm(input.file = "LDpruned_ped.ped") #conversion only needs to be done once

# snmf for admixture proportions
pearl_lea = snmf(input.file = "LDpruned_ped.lfmm", K = 1:19, entropy = TRUE, project = "continue", CPU = 10, repetitions=10) #K = 18 for number populations?)
##plot cross-entropy criterion of all runs of the project
plot(pearl_lea, col = "blue", pch = 19, cex = 1.2)  #big declines until 4
#ggsave("figures/snmf_scree.png") 
# summary of the project
summary(pearl_lea)
snmf_summary<- summary(pearl_lea) 
snmf_summary$crossEntropy[,8]
best_K<-as.numeric(which.min(snmf_summary$crossEntropy[1,]))# K = 8 = BEST, but big declines until 4
best_K #9

snmf_Qscores_meta <- function(leaproj, Knum, meta) {
  leaQ <-Q(object = leaproj, K = Knum, run = 1)
  colnames(leaQ) <- paste0("Q", rep(1:Knum))
  leaQ_meta <- bind_cols(meta, leaQ)
  return(leaQ_meta) }  # a function to import Q scores from snmf and match to metadata

Make_admix_table<-  function(Admix_meta_object, Knum){Admix_table <- Admix_meta_object
rownames(Admix_table) <-Admix_table$ID
plot_data <-  Admix_table  %>% 
  gather('pop', 'prob', Q1:paste0("Q", Knum))
return(plot_data)} #a function to make a plottable table for snmf 

meta_reduced <- meta.fam %>%
  dplyr::select(ID, SiteCode, AnonSiteCode, VCNAME) %>%
  droplevels()

# snmf.K4.meta <- snmf_Qscores_meta(leaproj = pearl_lea, Knum = 4, meta = meta.fam)
snmf.Kbest.meta <- snmf_Qscores_meta(leaproj = pearl_lea, Knum = best_K, meta = meta.fam)

#write_csv(snmf.K4.meta, col_names = T, "snmf.K4.meta.csv")
#meta.fam <- fread("snmf.K4.meta.csv")
#head(snmf.K4.meta)
# #### lea documentaion code ####
# project=pearl_lea
# # plot cross-entropy criterion of all runs of the project
# plot(project, cex = 1.2, col = "lightblue", pch = 19)
# 
# # get the cross-entropy of all runs for K = 4
# ce = cross.entropy(project, K = 4)
# 
# # select the run with the lowest cross-entropy for K = 4
# best = which.min(ce)
# 
# # display the Q-matrix
# my.colors <- c("tomato", "lightblue",
#                "olivedrab", "gold", "pink", "purple", "brown")
# barchart(project, K = 4, run = best,
#          border = NA, space = 0,
#          col = my.colors,
#          xlab = "Individuals",
#          ylab = "Ancestry proportions",
#          main = "Ancestry matrix") -> bp
# axis(1, at = 1:length(bp$order),
#      labels = bp$order, las=1,
#      cex.axis = .4)
# ###########

########
snmf.Kbest.admixtable <- Make_admix_table(snmf.Kbest.meta, Knum = best_K)

##########################################
########### admixture plot ###############
##########################################
## reorder levels by region and where fish would hit first
## MT-QC-NL
# snmf.K4.admixtable$SiteCode <- factor(snmf.K4.admixtable$SiteCode, levels = c("NSH","MSW","MUN","UPS", "CMP","CNR","NPR","TNR","WAB","SH","ENG"))

snmf.K4.meta <- as.data.frame(snmf.K4.meta)
myvars <- c("ID","SiteCode", "AnonSiteCode","VCNAME", "Q1","Q2","Q3","Q4")
snmf.K4.meta.reduced <- snmf.K4.meta[myvars]
snmf.K4.meta.reduced <- droplevels(snmf.K4.meta.reduced)
# snmf.K4.meta.reduced$SiteCode <- factor(snmf.K4.meta.reduced$SiteCode, levels = c("NSH","MSW","MUN","UPS", "CMP","CNR","NPR","TNR","WAB","SH","ENG"))

# Convert dataframe to long format
# library(reshape2)
qlong <- reshape2::melt(snmf.K4.meta.reduced, id.vars=c("ID","AnonSiteCode","VCNAME", "SiteCode")) #Make sure have ALL id vars in here!

nameslong <- reshape2::melt(snmf.K4.meta.reduced, id.vars=c("ID","SiteCode","VCNAME"))
# write_csv(qlong, col_names = T, "../../output/anon_qlong_maf02.csv")

# Define colour palette
# pal = colorRampPalette(c("green","blue", "deeppink", "orange"))
# pal = colorRampPalette(c("#ce1256","#df65b0", "#d7b5d8", "#f1eef6"))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(length(unique(qlong$variable)))

library(paletteer)# for the classicpurple colours
# cols = pal(length(unique(qlong$variable)))
# Create admixture plot
admix.bar = ggplot(data=qlong, aes(x=ID, y=value, fill=variable))+
  geom_bar(stat = "identity", width = 1)+
  scale_fill_paletteer_d("ggthemes::Classic_Purple_Gray_12") +
  scale_y_continuous(expand = c(0,0))+
  facet_grid(~AnonSiteCode,scales = "free", switch = "x", space = "free")+
  # scale_fill_manual(values = cols)+
  ylab("Admixture proportion")+
  # xlab("Individual")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(colour="black", size=12),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
admix.bar

##===============================
### Attempt multiple admix plots
# Colour palette function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# install.packages("paletteer")
library(paletteer)

scale_colour_paletteer_d("ggthemes::Classic_Purple_Gray_12") #8074a8 is main purple
# Function to generate admixture plot for a given K
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
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text = element_text(colour = "black", size = 8,
                                angle=90),  # Smaller facet text
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.spacing = unit(0.05, "cm"),
      # panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Border around facets
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
# ggsave(plot= opt_k_pl, 
#        width = 1000/100,
#        height = 360/100,
#        units = c( "in"),
#        dpi = 600,
#        filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/figures/optimumK_snmf.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/figures/optimumK_snmf.png")

# Generate admixture plots for K = 4, 5, 6, 7, 8, 9
# plots <- lapply(4:9, function(k) generate_admix_plot(k, meta_reduced))
# # plots
# # Combine all plots in a single layout
# combined_plot <- wrap_plots(plots, ncol = 2, nrow = 3)+ plot_annotation(tag_levels = 'a')
# combined_plot

# Generate admixture plots for K = 4, 5, 6, 7, 8, 9, 10 ,11
plots <- lapply(5:12, function(k) generate_admix_plot(k, meta_reduced))
# plots
# Combine all plots in a single layout
combined_plot <- wrap_plots(plots, ncol = 2, nrow = 4)+ plot_annotation(tag_levels = 'a',    tag_suffix = ") "
)
combined_plot
ggsave(filename = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/figures/all_snmf_admix_plots.png", plot = combined_plot,
       width = 1741 / 100,  # Convert pixels to inches (assuming 100 dpi)
       height = 1110 / 100,  # Convert pixels to inches (assuming 100 dpi)
       units = "in",  # Specify units as inches
       dpi = 600  # Specify resolution)
)
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/figures/all_snmf_admix_plots.png")

#======================================================================================
#### Off the wall pie charts https://github.com/Tom-Jenkins/admixture_pie_chart_map_tutorial
# ----------------- #
#
# Prepare pie charts
#
# ----------------- #

# Calculate mean admixture proportions for each site
dyn_admix_reduc <- function(Knumber, meta) {
  snmf_meta <- snmf_Qscores_meta(leaproj = pearl_lea, Knum = Knumber, meta = meta)
  snmf_meta_test <-as.data.frame(snmf_meta)
  qlong <- reshape2::melt(snmf_meta_test, id.vars = c("ID", "AnonSiteCode", "VCNAME", "SiteCode"))
  return(snmf_meta_test)
}
snmf.K4.meta.reduced<-dyn_admix_reduc(best_K,meta_reduced) #too lazy to chenge the K4

head(snmf.K4.meta.reduced)
clusters = grep("Q", names(snmf.K4.meta.reduced)) # indexes of cluster columns
avg_admix = aggregate(snmf.K4.meta.reduced[,clusters], list(snmf.K4.meta.reduced$AnonSiteCode), mean)

# Order alphabetically by site
avg_admix = avg_admix[order(as.character(avg_admix$Group.1)), ]
avg_admix_perc<-avg_admix
avg_admix_perc[, 2:ncol(avg_admix_perc)] <- round(avg_admix[, 2:ncol(avg_admix_perc)] * 100, 2)
avg_admix_perc
avg_admix_perc$Ancestry <- apply(avg_admix_perc[, 2:ncol(avg_admix_perc)], 1, function(row) {
  max_value <- max(row)
  if (max_value > 75) {
    return(names(row)[which.max(row)])  # Assign the ancestry with the highest value
  } else {
    return("Mixed")  # Assign "Mixed" if no ancestry is over 75%
  }
})

# Display the updated data frame
avg_admix_perc
write.table(avg_admix_perc, file = "C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/avg_admix_perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/avg_admix_perc.tsv")

# Convert dataframe from wide to long format
avg_admix = melt(avg_admix, id.vars = "Group.1")
head(avg_admix)
avg_admix
cols <- ggthemes::Classic_Purple_Gray_12[1:8]  # Adjust based on available colours

# Define a function to plot pie charts using ggplot for each site
pie_charts = function(admix_df, site, cols){
  # admix_df = dataframe in long format of admixture proportions per site 
  # site = string 
  # cols = vector of colours of length(clusters)
  ggplot(data = subset(admix_df, Group.1 == site),
         aes(x = "", y = value, fill = variable))+
    geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE)+
    coord_polar(theta = "y")+
    scale_fill_paletteer_d("ggthemes::Classic_Purple_Gray_12") +
    # scale_fill_manual(values = cols)+
    theme_void()
}

# Test function on one site
pie_charts(avg_admix, site = "EI1", cols=cols)

# Apply function to all sites using for loop
pies = list()
for (i in avg_admix$Group.1){
  pies[[i]] = pie_charts(admix_df = avg_admix, site = i, cols = cols) 
}

 ###################
#Saving indi pies
# pie_charts = function(admix_df, site, cols){
#   # admix_df = dataframe in long format of admixture proportions per site 
#   # site = string 
#   # cols = vector of colours of length(clusters)
#   p = ggplot(data = subset(admix_df, Group.1 == site),
#              aes(x = "", y = value, fill = variable))+
#     geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE)+
#     coord_polar(theta = "y")+
#     scale_fill_manual(values = cols)+
#     theme_void()+
#     ggtitle(site)  # Set title based on SiteCode
#   return(p)
# }
# # Test function on one site
# pie_charts(avg_admix, site = "WS1", cols=cols)
# 
# # Apply function to all sites and save each pie chart
# output_dir <- "../../output/admix_pies"  # Specify your directory path
# 
# for (site in avg_admix$Group.1) {
#   p <- pie_charts(admix_df = avg_admix, site = site, cols = cols)
#   filename <- paste0(output_dir, "/pie_", gsub(" ", "_", site), ".png")  # Generate filename with SiteCode
#   ggsave(filename, plot = p, width = 6, height = 6)  # Adjust width and height as needed
# }


# ----------------- #
#
# Prepare basemap
#
# ----------------- #

# Obtain coordinates

coords<-meta.fam%>% group_by(SiteCode) %>%
  summarise(Latitude = first(Latitude), Longitude = first(Longitude))
coords <- meta.fam %>%
  group_by(SiteCode) %>%
  mutate(Latitude = first(Latitude), Longitude = first(Longitude)) %>%
  distinct(SiteCode, .keep_all = TRUE) %>%
  select(SiteCode, Population, VCNAME, Latitude, Longitude, AnonSiteCode)

# Order alphabetically by site
coords = coords[order(coords$SiteCode), ] 
coords

# Check order matches coords order
as.character(avg_admix$Group.1) == as.character(coords$SiteCode)

# Set map boundary (xmin, xmax, ymin, ymax)
# install.packages("rworldmap")
# install.packages("rworldxtra")
library(rworldxtra)
library(rworldmap)
library(sf)
library(raster)
boundary = extent(-10, 0, 55, 60)
boundary

# Get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()
# install.packages("ggsn")
# library(ggsn)

# Plot basemap
basemap = ggplot()+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="grey",
               colour="black", size=0.5)+
  coord_quickmap(expand=F)+
  # ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
  # ggsn::scalebar(data = map.outline, dist = 200, dist_unit = "km", height = 0.01,
  #                transform = TRUE, model = "WGS84", 
  #                location = "bottomleft", anchor = c(x = -12.5, y = 45),
  #                st.bottom = FALSE, st.size = 4, st.dist = 0.015)+
  xlab("Longitude")+
  ylab("Latitude")

basemap
polys <- st_read("../../data/scotland_vice_counties/scotland_vicecounties.shp") %>% st_transform(crs = 4326)

ggplot() +
  geom_sf(data = polys, fill = "grey", colour = "black", size = 0.5, alpha = 0.9) +
  coord_sf(expand = FALSE) +
  xlab("Longitude") +
  ylab("Latitude")
ggplot() +
  geom_sf(data = polys, fill = "#F0F0F0", colour = "#CCCCCC", size = 0.2, alpha = 0.7) +  # Soft grey fill and lighter grey lines
  coord_sf(expand = FALSE) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_minimal() +  # Optional: Use a minimal theme
  theme(panel.background = element_rect(fill = "transparent", colour=NA)) 

# ----------------- #
#
# Add pie charts to basemap
#
# ----------------- #

# Extract coordinates for each site
rm(coord.list)
coord.list = list()
for (i in avg_admix$Group.1){
  coord.list[[i]] = c(subset(coords, Population == i)$Longitude, subset(coords, Population == i)$Latitude)
}
coord.list

# Define pie chart sizes
radius = 0.2

# Convert ggplot pie charts to annotation_custom layers
length(avg_admix$Group.1)
length(pies)

pies.ac = list()
for (i in 1:min(length(avg_admix$Group.1), length(pies))) {
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}


# Add layers to basemap
pie.map = basemap + pies.ac
pie.map
 
############# lets randomly move pies and annotate
# Load necessary libraries
library(ggplot2)
library(ggrepel)

# Assuming bs is your original basemap with the polygons plotted

# Define jitter range (adjust as needed)
jitter_range <- 0.05

# Randomly jitter coordinates
jittered_coords <- lapply(coord.list, function(coord) {
  jittered_x <- coord[[1]] + runif(1, -jitter_range, jitter_range)
  jittered_y <- coord[[2]] + runif(1, -jitter_range, jitter_range)
  list(x = jittered_x, y = jittered_y)
})

# Convert ggplot pie charts to annotation_custom layers
pies.ac <- list()
for (i in 1:min(length(avg_admix$Group.1), length(pies))) {
  pies.ac[[i]] <- annotation_custom(grob = ggplotGrob(pies[[i]]),
                                    xmin = jittered_coords[[i]]$x - radius,
                                    xmax = jittered_coords[[i]]$x + radius,
                                    ymin = jittered_coords[[i]]$y - radius,
                                    ymax = jittered_coords[[i]]$y + radius)
}

# Add layers to basemap
bs1 <- basemap +
  pies.ac
bs1
# Create a data frame for labels
labels <- data.frame(
  Longitude = unlist(lapply(jittered_coords, `[[`, "x")),
  Latitude = unlist(lapply(jittered_coords, `[[`, "y")),
  AnonSiteCode = coords$AnonSiteCode[match(avg_admix$Group.1, coords$SiteCode)]
)%>%
  distinct(AnonSiteCode, .keep_all = TRUE)

# Annotate with AnonSiteCode using geom_label_repel
bs2 <- bs1 +
  geom_label_repel(data = labels, aes(x = Longitude, y = Latitude, label = AnonSiteCode),
                   box.padding = 0.5, point.padding = 0.3,
                   segment.color = "black", segment.size = 0.5,
                   size = 3, colour = "black",  force=20, max.overlaps = Inf)

# Display the map with pie charts and annotations
print(bs2)

### JUST TEXT LABELS NO PIE CHARTS
# Annotate with AnonSiteCode using geom_label_repel
lab1<- basemap +
  geom_label_repel(
    data = labels, 
    aes(x = Longitude, y = Latitude, label = AnonSiteCode),
    box.padding = 0.1,       # Reduce the padding around the label box
    point.padding = 0.3,     # Reduce the padding around the point
    # segment.color = NA,      # Remove the line segments by setting color to NA
    size = 3, 
    colour = "black", 
    max.overlaps = Inf
  )

## with vice county data
# install.packages("smoothr")
# library(smoothr)
# smoothed_polys <- smooth(polys, method = "ksmooth", smoothness = 3)
# install.packages("rmapshaper")
library(rmapshaper)
smoothed_polys <- ms_simplify(polys, keep = 0.005)  # 25% of points retained

labvc<- ggplot() +
   geom_sf(data = smoothed_polys) + # Plot the polygons
   geom_label_repel(
     data = labels, 
     aes(x = Longitude, y = Latitude, label = AnonSiteCode),
     box.padding = 0.1,       # Reduce the padding around the label box
     point.padding = 0.1,     # Reduce the padding around the point
     segment.color = NA,      # Remove the line segments by setting color to NA
     size = 3, 
     colour = "black", 
     max.overlaps = Inf
   ) +
   theme_void()
labvc
# library(patchwork)
lab1 | labvc
#############
#Try to make pies not overlap -
# Assuming bs1 is your original basemap with the polygons added
# Also assuming you have defined jittered_coords and pies
# Function to calculate positions with minimum distance
# Function to calculate positions with minimum distance
# Function to adjust coordinates to ensure minimum distance
# Example function to adjust coordinates to prevent overlap

for (group in avg_admix$Group.1) {
  subset_coords <- subset(coords, Population == group)
  coord.list[[as.character(group)]] <- subset_coords[, c("Longitude", "Latitude")]
}

coord.list

adjust_coordinates <- function(coord_list, radius = 0.5) {
  for (i in seq_along(coord_list)) {
    if (i > 1) {
      # Compare with previous coordinates
      while (check_overlap(coord_list[[i]], coord_list[[i-1]], radius)) {
        # Adjust current coordinates (longitude and latitude)
        coord_list[[i]]$Longitude <- coord_list[[i]]$Longitude + 0.1  # Example adjustment
        coord_list[[i]]$Latitude <- coord_list[[i]]$Latitude + 0.1    # Example adjustment
      }
    }
  }
  return(coord_list)
}

# Function to check overlap between two sets of coordinates
check_overlap <- function(coord1, coord2, radius) {
  distance <- sqrt((coord2$Longitude - coord1$Longitude)^2 + (coord2$Latitude - coord1$Latitude)^2)
  return(distance < 2 * radius)  # Check if distance is less than twice the radius
}

# Adjusting coordinates in coord.list to prevent overlap
adjusted_coord_list <- adjust_coordinates(coord.list)

# Output adjusted_coord_list to check the updated coordinates
adjusted_coord_list

coord.list <- lapply(adjusted_coord_list, function(df) {
  c(df$Longitude, df$Latitude)
})
coord.list

############ another go
# Step 1: Define your functions
# Function to check overlap between two sets of coordinates
check_overlap <- function(coord1, coord2, radius) {
  distance <- sqrt((coord2[1] - coord1[1])^2 + (coord2[2] - coord1[2])^2)
  return(distance < 2 * radius)  # Check if distance is less than twice the radius
}

# Function to adjust coordinates to prevent overlap
adjust_coordinates <- function(coord_list, radius = 0.5) {
  for (i in 2:length(coord_list)) {
    while (check_overlap(coord_list[[i]], coord_list[[i-1]], radius)) {
      # Adjust current coordinates (Longitude and Latitude)
      coord_list[[i]][1] <- coord_list[[i]][1] + 0.1  # Example adjustment for Longitude
      coord_list[[i]][2] <- coord_list[[i]][2] + 0.1  # Example adjustment for Latitude
      # You may adjust these increments (0.1) based on your map scale and density
    }
  }
  return(coord_list)
}

# Step 2: Adjust coordinates in your original `coord.list`
adjusted_coord_list <- adjust_coordinates(coord.list)

# Step 3: Create annotation layers for pie charts with adjusted coordinates
library(ggplot2)
library(gridExtra)

# Assuming `pies` contains your pie chart objects and `coord_list` is the adjusted coordinates
pies.ac <- list()
radius <- 0.2  # Define pie chart radius

for (i in 1:min(length(avg_admix$Group.1), length(pies))) {
  pies.ac[[i]] <- annotation_custom(grob = ggplotGrob(pies[[i]]),
                                    xmin = adjusted_coord_list[[i]][1] - radius,
                                    xmax = adjusted_coord_list[[i]][1] + radius,
                                    ymin = adjusted_coord_list[[i]][2] - radius,
                                    ymax = adjusted_coord_list[[i]][2] + radius)
}

# Step 4: Combine pie chart layers with your basemap (`basemap`)
pie.map <- basemap + pies.ac

# Step 5: Visualise the map
print(pie.map)


###########

##Zoomed map
# install.packages("ggspatial")
zoom.map<- ggplot()+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="grey",
               colour="black", size=0.5)+
  coord_quickmap(xlim = c(-5.5, -4.75), ylim = c(58, 58.5),expand=F)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggspatial::annotation_scale(location = "tr", plot_unit = "km")
# Define pie chart sizes
radius = 0.025

# Convert ggplot pie charts to annotation_custom layers
pies.ac = list()
for (i in 1:length(avg_admix$Group.1)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}

# Add layers to basemap
zoomed<-zoom.map+pies.ac
zoomed


##Zoomed map
zoom.map<- ggplot()+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="grey",
               colour="black", size=0.5)+
  coord_quickmap(xlim = c(-6.5, -5.5), ylim = c(56.5, 57),expand=F)+
  xlab("Longitude")+
  ylab("Latitude")
# Define pie chart sizes
radius = 0.025

# Convert ggplot pie charts to annotation_custom layers
pies.ac = list()
for (i in 1:length(avg_admix$Group.1)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}

# Add layers to basemap
zoomed<-zoom.map+pies.ac
zoomed

# Combine ggplots
ggarrange(zoomed + theme(legend.position = "right") + labs(title = "Individual admixture proportions", tag = "A"),
          pie.map + labs(title = "Mean admixture proportions for each site", tag = "B"))
ggsave("3.Admixture_bar_map.png", width = 15, height = 10, dpi = 600)





### # How many SNPs do we have? 
# Turn ped and map file into something more readable
#' Read a ped and map file into R, produce a dataframe for general use so that we can then filter so metadata matches genotype data
readPedMap_tsv_fmt = function(ped_file, map_file = "", headers = FALSE){
  if(headers == TRUE){
    #
    map_data = read_tsv(map_file, skip = 1, col_names = c("chromosome", "snp" , "genetic_distance", "physical_distance"))
    
    header_data = c('#family', 'individual', 'sire', 'dam', 'sex', 'pheno')
    header_data = c(header_data, map_data$snp)
    
    snp_data = read_tsv(ped_file, skip = 1, col_names = header_data)
    return(snp_data)
    
  }else{
    
    map_data = read_tsv(map_file, col_names = c("chromosome", "snp" , "genetic_distance", "physical_distance"))
    
    header_data = c('#family', 'individual', 'sire', 'dam', 'sex', 'pheno')
    header_data = c(header_data, map_data$snp)
    
    snp_data = read_tsv(ped_file, col_names = header_data)
    return(snp_data)
  }
  
}

ped_file <- 'LDpruned_ped.ped'
map_file <- 'LDpruned_ped.map'
snp_data <- readPedMap_tsv_fmt(ped_file, map_file)
gc()
snp_data <- as_tibble(snp_data)
# get the SNP columns
snp_columns <- names(snp_data)[7:length(names(snp_data))]
length(snp_columns) # 3456 SNPs and 156 individuals