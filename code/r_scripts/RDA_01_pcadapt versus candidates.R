### PCAdpat with the candidate values 
library(pcadapt)
library(dplyr)
library(ggplot2)
load(file = "data/confidential/RDA_host_condition_updated.RData")
# These are the 3 types of candidate files
head(threeSD_snps_bind)
head(top_one_percent_snps_rda1)
head(outliers_below_pvalue)

dir<-"C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned"
plinkdir<-"C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/data/plink_VP/"

filename<-paste0(plinkdir, "PM_Dec4")

PCADAPT_import_PCA <- function(filename, Knum) {
  to_pc <- read.pcadapt(paste0(filename, ".bed"), type = "bed")
  PCA <- pcadapt(to_pc, K = Knum, min.maf = 0.02)
  return(PCA)
} # a function to import data from .bed file and run a PCA with specific number of components and minimum maf of 0.01
pcad_import <- PCADAPT_import_PCA(filename = filename, Knum = 18)
plot(pcad_import , option = "manhattan")
plot(pcad_import, option = "qqplot")
hist(pcad_import$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
pcad_import$pvalues
pcad_import$pass

padj <- p.adjust(pcad_import$pvalues, method = "bonferroni")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)

#  Load .bim file with SNP metadata
library(data.table)
bim_data <- fread(paste0(filename, ".bim"), header = FALSE)
colnames(bim_data) <- c("chr", "SNP", "cm", "bp", "A1", "A2")

# Step 2: Combine SNP positions with p-values from pcadapt
pcadapt_snp_info <- tibble(
  chr = bim_data$chr,
  SNP = bim_data$SNP,
  pos = 1:length(bim_data$SNP),  # Use sequential index instead of genomic positions
  pvalue = pcad_import$pvalues,
  padj = p.adjust(pcad_import$pvalues, method = "bonferroni")
)

# Step 2: Filter for significant SNPs (padj < alpha, where alpha is 0.1)
alpha <- 0.05
alpha <- 0.05/length(pcad_import$pvalues) 

pcadapt_significant_snps <- pcadapt_snp_info %>%
  filter(padj < alpha) %>% as.data.frame()
head(pcadapt_significant_snps)

saveRDS(pcadapt_significant_snps, "output/pcadapt_significant_snps_bonferroni_adjusted_under_0.05.rds")

# Step 5: Create Manhattan plot
man<-ggplot(pcadapt_snp_info, aes(x = pos, y = -log10(padj))) +
  geom_point(color = "grey80") +  # Plot all SNPs with normal scaling
  geom_point(data = pcadapt_significant_snps, colour = "#14a8ff", size = 2) +  # Overlay candidate SNPs with special scaling
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = "red") +
  labs(
    x = "Position in Genome", y = "-log10(p.values)" ) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  # Ensure x-axis starts at 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50))
man
ggsave(man, width = 8.5,
       height = 3.5,
       units = c( "in"),
       dpi = 600,
       filename = "./output/figures/manhatten_pcadapt_pval.png")


#comparing to SD and rda p vals
head(Out_pv)
head(Out_3SD)
head(pcadapt_snp_info)
# Filter pcadapt_snp_info for SNPs in Out_3SD
Out_3SD_in_pcad <- pcadapt_snp_info %>%
  filter(SNP %in% sub("_[ACGT]$", "", rownames(Out_3SD))) %>% as.data.frame()
Out_pv_in_pcad <- pcadapt_snp_info %>%
  filter(SNP %in% sub("_[ACGT]$", "", rownames(Out_pv))) %>% as.data.frame()
dim(Out_3SD_in_pcad)
dim(Out_3SD)
head(Out_3SD)
head(Out_3SD_in_pcad)
#for some reason there is a duplicate
man<-ggplot(pcadapt_snp_info, aes(x = pos, y = -log10(padj))) +
  geom_point(color = "grey80") +  # Plot all SNPs with normal scaling
  geom_point(data = Out_pv_in_pcad, colour = "deeppink", size = 4) +  
  geom_point(data = Out_3SD_in_pcad, colour = "#ff8ac9", size = 2) +  
  geom_point(data = pcadapt_significant_snps, colour = "#111", size = 3, alpha=1, pch=3) +  
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = "red") +
  labs(
    x = "Position in Genome", y = "-log10(p.values)" ) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  # Ensure x-axis starts at 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75))
man
ggplot(pcadapt_snp_info, aes(x = pos, y = -log10(padj))) +
  geom_point(aes(colour = "All SNPs"), size = 1) +  # Assign a label for the legend
  geom_point(data = Out_pv_in_pcad, aes(colour = "Outliers (PV)"), size = 4) +
  geom_point(data = Out_3SD_in_pcad, aes(colour = "Outliers (3SD)"), size = 2) +
  geom_point(data = pcadapt_significant_snps, aes(colour = "Significant SNPs"), size = 3, alpha = 1, pch = 3) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = "red", aes(colour = "Threshold")) +
  labs(
    x = "Position in Genome", 
    y = "-log10(p.values)", 
    colour = "Legend"  # Title for the legend
  ) +
  theme_classic() +
  scale_colour_manual(
    values = c(
      "All SNPs" = "grey80",
      "Outliers (PV)" = "deeppink",
      "Outliers (3SD)" = "#ff8ac9",
      "Significant SNPs" = "#111",
      "Threshold" = "red"
    )
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75))

head(Out_3SD_in_pcad)
Out_3SD_below_pcadapt_pthres <- Out_3SD_in_pcad %>%
  filter(padj < alpha) %>% as.data.frame()
Out_pv_below_pcadapt_pthres <- Out_pv_in_pcad %>%
  filter(padj < alpha) %>% as.data.frame()
pv_in_3SD <- Out_3SD_in_pcad %>%
  filter(Out_pv_in_pcad == TRUE) %>%
  as.data.frame()
pv_in_3SD <- Out_3SD_in_pcad %>%
  filter(SNP %in% Out_pv_in_pcad$SNP) %>% as.data.frame()
pv_sd_overlap<-(Out_3SD_below_pcadapt_pthres %>%
  filter(SNP %in% Out_pv_below_pcadapt_pthres$SNP) %>% as.data.frame())
pv_sd_overlap_below_pcadapt_pthres <- pv_sd_overlap %>%
  filter(padj < alpha) %>% as.data.frame()
nrow(pv_sd_overlap_below_pcadapt_pthres)
# P-values threshold after Bonferroni correction
thres_env <- 0.05/length(pcad_import$pvalues) 
thres_env

pcadapt_out_pvalue <- pcad_import$pvalues[pcad_import$pvalues < thres_env]
pcadapt_out_pvalue
pcadapt_out_pvalue <- pcad_import$pvalues %>%
  # arrange(desc(p_values)) %>%
  filter(pvalues < thres_env)
