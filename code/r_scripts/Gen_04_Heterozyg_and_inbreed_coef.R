# https://github.com/ruthrivkin/Polarbear-maladaptation/blob/main/PopGen.R
# Script modified from above to get more basic pop gen stats


#Calculate basic stats
 # install.packages("dartR")
 # BiocManager::install("SNPRelate")
 # gl.install.vanilla.dartR()
library(hierfstat)
library(dartR)
library(StAMPP)
setwd("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/data/plink_VP/")

# packages <- c("hierfstat", "dartR", "StAMPP", "radiator", "tidyr", "readr")
# # Create a directory for citations if it doesn't exist
# if (!dir.exists("../../citations")) {
#   dir.create("citations")
# }
# 
# # Save each package citation as a separate BibTeX file
# for (pkg in packages) {
#   bib_file <- paste0("citations/", pkg, "_citation.bib")
#   write(capture.output(toBibtex(citation(pkg))), 
#         file = bib_file)
# }



#for anon names
# Rename on anon pop names
library(dplyr)
vicec<-read.csv("../../data/confidential/trimmed_env_data.csv")%>%
  mutate(Sample_ID = ifelse(Sample_ID == "ICY837","ICY837notICY773", Sample_ID)) %>% mutate(Population = ifelse(Population == "Conon_Gharbhrain_wrong_sample_ID", "Conon_Glascarnoch", Population))%>% mutate(AnonSiteCode = ifelse(AnonSiteCode == "ER2", "ER1", AnonSiteCode)) %>% dplyr::select(Population, AnonSiteCode)%>%
  dplyr::rename(subpop = Population) #pull in anon site names
head(vicec)


# Read in data
rawfile <- "LDpruned_recodeA.raw"
myData <- read.PLINK(rawfile)
myData.nona <- gl.impute(myData, method = "HW")
unique(myData.nona@pop)
myData


##### Calculate private alleles #####################
myData <- gl.compliance.check(myData)
pa <- gl.report.pa(myData, method = "one2rest")
pa
pa <- pa %>%
  mutate(AnonSiteCode = vicec$AnonSiteCode[match(pop1, vicec$subpop)])#give table anon site names
length(unique(pa$AnonSiteCode)) #should be 18
# Remove columns pop1, p1, and p2, then rename AnonSiteCode and reorder alphabetically
pa_nice <- pa %>%
  dplyr::select(-pop1, -p1, -p2) %>%            # Remove columns pop1, p1, and p2
  rename(pop1 = AnonSiteCode) %>%     # Rename AnonSiteCode to pop1
  dplyr::select(pop1, everything()) %>%      # Move pop1 to be the first column
  arrange(pop1)   
# Display the first few rows to confirm changes
head(pa_nice)

head(pa)

write.csv(pa_nice, "../../data/confidential/privatealleles.ld.csv", row.names = FALSE)

##### Calcuclate nucleotide diversity ##############
# library(BiocManager)
# install.packages("radiator")
# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("thierrygosselin/radiator")
# devtools::package_info(pkgs = "SeqArray") # to verify version
# # If manually installing SeqArray is necessary
# install.packages("BiocManager")
# BiocManager::install("SeqArray")
library(radiator)
library(tidyr)
library(readr)

#Make strata file with ind id and subset info (POP_ID and INDIVIDUALS)
df <- fread("LDpruned.fam") %>% as.data.frame() %>% dplyr::select(V1, V2)
df_strata <- df %>%
  mutate(SampleID = paste(V1, V2, sep = "_")) %>%  # Combine Pop and ID
  dplyr::select(INDIVIDUALS = SampleID , STRATA = V1)                   # Select and rename columns
head(df_strata)
write.table(df_strata, file = "strata.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Make LD pruned VCF files for radiator to work
# system("./plink --file LDpruned_ped  --out LD_VCF  --allow-extra-chr --recode vcf") #--missing
# system("cmd.exe /c dir")

# Check Strata works
radiator::summary_strata("strata.tsv") #care strata matches the col/id names in the vcf

# Read the vcf- if it doesnt work it will create hell files in the wd called read_vcf_DATE@TIME that you will struggle to delete :) If it works you will get vcf in your environment
vcf <- radiator::read_vcf(data = "LD_VCF.vcf", strata = "strata.tsv", verbose = T, parallel.core = 1L) # error in `.DynamicClusterCall()`: needs the paralell core
tidyg <- radiator::tidy_genomic_data(data = "LD_VCF.vcf", strata = "strata.tsv", verbose = T, parallel.core = 1L) # error in `.DynamicClusterCall()`: needs the paralell core

# remove(pi)
vcf$filename
pi <- radiator::pi(vcf, verbose=T)

pi$pi.populations
pi$boxplot.pi
as.data.frame(pi$pi.individuals) %>% group_by(POP_ID) %>%
  summarise(PI = mean(PI), .groups = 'drop') #


pipop <- as.data.frame(pi$pi.populations)  %>%
  mutate(AnonSiteCode = vicec$AnonSiteCode[match(POP_ID, vicec$subpop)]) %>%
  arrange(AnonSiteCode) %>%
  dplyr::select(-POP_ID)
pipop
#give table anon site names
head(pipop)

pidf <- as.data.frame(pi$pi.individuals)
head(pidf)
# Extract the ID portion and add as a new column
library(stringr)
pidf2 <- pidf %>%  mutate(IID = str_extract(INDIVIDUALS, "IC[A-Z]\\d+")) #care code is correct so doenst give NAs
sum(is.na(pidf2$IID))  # Counts NA values in ID column
head(pidf2)

# Add Anon names 
pidf3 <- pidf2 %>%
  mutate(AnonSiteCode = vicec$AnonSiteCode[match(POP_ID, vicec$subpop)])#give table anon site names
length(unique(pidf3$AnonSiteCode)) #should be 18
individual_summary <- pidf3 %>%
  group_by(POP_ID) %>%
  summarise(PI = mean(PI, na.rm = TRUE), .groups = 'drop')

print(individual_summary)

meanpi<- Rmisc::summarySE(pidf3, measurevar = "PI", groupvars = "AnonSiteCode", na.rm = TRUE)

write.csv(meanpi, "../../output/Nucleotide_diversity_ldpruned.csv",row.names = FALSE)
meanpi<-read.csv("../../output/Nucleotide_diversity_ldpruned.csv")
head(meanpi)

####### Heterozygosity and Inbreeding  ########################################
# system("./plink --file LDpruned_ped  --out ld_pruned_het --allow-extra-chr --het --recode vcf") #--missing 
# hetdata <- read.delim("ld_pruned_het.het", sep = "")
# str(hetdata)
# head(hetdata)
# View(hetdata) 
# #Prep het file
# het <- hetdata %>% rename(subpop = FID)
# head(het)
# het <- het %>% #this doesnt work because its using plink individual estimates and then I average them out for population calcs - population calcs should have a different equation which dartR uses
#   mutate(
#     prob.obs.hom = O.HOM./N.NM., 
#     prob.exp.hom = E.HOM./N.NM.,
#     prob.obs.het = 1 - prob.obs.hom,
#     prob.exp.het = 1 - prob.exp.hom
#   )
# head(het)
# str(het)
# range(het$prob.exp.het)
# ggplot(het, aes(x = prob.exp.het, y = prob.obs.het, colour = subpop)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
#   theme_minimal() +
#   labs(x = "Expected heterozygosity", y = "Observed heterozygosity",
#        title = "Observed vs Expected Heterozygosity by Subpopulation")
# 

###another version of heterozyg
# Read VCF as genlight
genlight_obj <- gl.read.vcf("ld_pruned_het.vcf")
head(genlight_obj)
indNames(genlight_obj)

pop(genlight_obj) <- sapply(indNames(genlight_obj), function(x) strsplit(x, "_")[[1]][1])
pop(genlight_obj) <- sapply(indNames(genlight_obj), function(x) {
  parts <- unlist(strsplit(x, "_"))
  paste(parts[-length(parts)], collapse = "_")
}) #split pop names from IDs by last _ in name

pop(genlight_obj)
dartrhet<-gl.report.heterozygosity(genlight_obj, method = "pop")
head(dartrhet)
write.csv(dartrhet, "../../data/confidential/dartR_Geneticdiversity_ld.csv", row.names = FALSE)
dartrhet_ind<-gl.report.heterozygosity(genlight_obj, method = "ind")
head(dartrhet_ind)
# View(dartrhet)
head(dartrhet)
head(vicec)

dartheter<-read.csv("../../data/confidential/dartR_Geneticdiversity_ld.csv")
all.equal(dartheter$He.adj, dartheter$He)
all.equal(dartheter$Ho.adj, dartheter$Ho)


###########################################

# write.csv(het, "../../data/confidential/Geneticdiversity_ld.csv", row.names = FALSE)
# 
het <- read.csv("../../data/confidential/Geneticdiversity_ld.csv")
head(het)
nrow(het)
mean(het$prob.obs.het)

mean(het$prob.exp.het)
# #Remame with anon-names
# het <- het %>%
#   mutate(AnonSiteCode = vicec$AnonSiteCode[match(subpop, vicec$subpop)])#give table anon site names
# length(unique(het$AnonSiteCode)) #should be 18
# head(het)
# # install.packages("Rmisc")
# meanhet <- Rmisc::summarySE(het, measurevar = "prob.obs.het", groupvars = "AnonSiteCode", na.rm = TRUE)
# meanhet
# head(meanhet)
# meanexphet <- Rmisc::summarySE(het, measurevar = "prob.exp.het", groupvars = "AnonSiteCode", na.rm = TRUE)
# meanexphet
# head(meanexphet)
# # write.csv(meanhet, "../../data/confidential/anon_heterozyg_by_pop.csv", row.names = FALSE)
# meanhet<-read.csv("../../data/confidential/anon_heterozyg_by_pop.csv")
# meanF <- Rmisc::summarySE(het, measurevar = "F", groupvars ="AnonSiteCode", na.rm = TRUE)
# meanF
# # write.csv(meanF, "../../data/confidential/anon_inbreeding_cooeff.csv", row.names = FALSE)
# meanF<-read.csv("../../data/confidential/anon_inbreeding_cooeff.csv")


#### Nice composite table ########################
pa_nice<-read.csv("../../output/privatealleles.ld.csv")
dartrhet<-read.csv("../../data/confidential/dartR_Geneticdiversity_ld.csv") %>% 
  mutate(AnonSiteCode = vicec$AnonSiteCode[match(pop, vicec$subpop)])
meanpi<-read.csv("../../output/Nucleotide_diversity_ldpruned.csv")
pipop<-read.csv("../../output/Nucleotide_diversity_ldpruned_pipop.csv")

catcha<-read.csv("../../data/confidential/trimmed_env_data.csv")%>%
  mutate(Sample_ID = ifelse(Sample_ID == "ICY837","ICY837notICY773", Sample_ID)) %>% mutate(Population = ifelse(Population == "Conon_Gharbhrain_wrong_sample_ID", "Conon_Glascarnoch", Population))%>% mutate(AnonSiteCode = ifelse(AnonSiteCode == "ER2", "ER1", AnonSiteCode)) %>% group_by(catchment_area_km, AnonSiteCode, Host) %>%
  summarise(catchment_area_km = mean(catchment_area_km, na.rm = TRUE),
            .groups = 'drop')
#%>% dplyr::select(Population, AnonSiteCode)%>%
  # dplyr::rename(subpop = Population) #pull in anon site names
head(catcha)
head(dartrhet)

# Merge and select relevant columns with meaningful CI names
final_table <- pa_nice %>%
  dplyr::select(-pop2,-N2, -fixed, -priv1, -priv2, -AFD) %>% 
  rename(AnonSiteCode = pop1, Private_Alleles = totalpriv, N=N1) %>%
  left_join(dartrhet %>% dplyr::select(-pop, -nInd, -nLoc, -nLoc.adj, -monoLoc, -polyLoc, -all_NALoc), by = "AnonSiteCode") %>%
  # left_join(meanhet[, c("AnonSiteCode", "prob.obs.het", "ci")]%>% rename(CI_Ho = ci, Ho = prob.obs.het), by = "AnonSiteCode") %>%
  # left_join(meanexphet[, c("AnonSiteCode", "prob.exp.het", "ci")]%>% rename(CI_He = ci, He=prob.exp.het), by = "AnonSiteCode") %>%
  # left_join(meanF[, c("AnonSiteCode", "F", "ci")]%>% rename(CI_F = ci), by = "AnonSiteCode") %>%
  left_join(pipop[, c("AnonSiteCode", "PI_NEI")], by = "AnonSiteCode") %>%  # Joining pipop for PI_NEI
  left_join(catcha[, c("AnonSiteCode", "catchment_area_km", "Host")], by = "AnonSiteCode") %>%  
  mutate(across(where(is.numeric) & !c(N, catchment_area_km), ~ round(.x, 3))) %>%
  mutate(across(catchment_area_km, ~ round(.x, 0))) %>% 
  # rename("Catchment area (km)" = catchment_area_km, "Private Alleles"= Private_Alleles) %>%
  arrange(AnonSiteCode) %>% 
  dplyr::select(AnonSiteCode, catchment_area_km, N, FIS, Ho.adj, Ho.adjSD, He.adj, He.adjSD, PI_NEI, Private_Alleles, Host)
head(final_table)

plot(final_table$FIS, final_table$F)
# plot(final_table$He.x, final_table$He.y) #comparing the plink Het and dartR het estimates
# plot(final_table$Ho.x, final_table$Ho.y)
# View result

write.csv(final_table, "../../output/anon_basic_stats.csv", row.names = FALSE)
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/anon_basic_stats.csv")
final_table<-read.csv("../../output/anon_basic_stats.csv")



summary(final_table)
ggplot(final_table)+
  geom_point(aes(N, Private_Alleles))
final_table$AnonSite[final_table$FIS >= 0.5]
plot(final_table$FIS, final_table$Ho)
final_table[, c("AnonSiteCode", "FIS")]

final_table[order(-final_table$FIS), c("AnonSiteCode", "FIS")]
final_table[order(-final_table$Ho), c("AnonSiteCode", "Ho")]
final_table[order(-final_table$He.adj), c("AnonSiteCode", "He.adj")]
final_table[order(-final_table$PI_NEI), c("", "PI_NEI")]
final_table[order(-final_table$Private_Alleles), c("", "Private_Alleles")]

###### Analysis #########################
library(ggplot2)
library(tidyr)

# Testing inverse relationship between F and Ho
final_table

# Scatter plot to check linearity
ggplot(final_table, aes(x = FIS, y = Ho.adj)) +
  geom_point() +
  geom_smooth(method = "lm", colour = "blue") +
  labs(x = "Fis (Inbreeding Coefficient)", y = "Ho (Observed Heterozygosity)") +
  theme_minimal()

# Normality tests
shapiro.test(final_table$Ho.adj)  # Normality of Ho.adj
shapiro.test(final_table$FIS)   # Normality of F

# Q-Q plots for normality visualization
qqnorm(final_table$Ho.adj); qqline(final_table$Ho.adj, col="red")
qqnorm(final_table$FIS); qqline(final_table$FIS, col="red")

# Boxplots to check for outliers
boxplot(final_table$Ho.adj, main="Boxplot of Ho.adj")
boxplot(final_table$FIS, main="Boxplot of Fis")

# Linear regression to check homoscedasticity
model <- lm(Ho.adj ~ FIS, data = final_table)
plot(model, which = 3)  # Scale-location plot

cor.test(final_table$Ho.adj, final_table$FIS, method = "pearson")

model <- lm(Ho.adj ~ FIS, data = final_table)
summary(model)
ggplot(final_table, aes(x = FIS, y = Ho.adj)) + #no statistical correlation between fis and Ho.adj
  geom_point() +
  geom_smooth(method = "lm", colour = "blue") +
  labs(x = "Inbreeding Coefficient (F)", y = "Observed Heterozygosity (Ho.adj)") +
  theme_minimal()

#for correlations  of catchment area to SNP stats
corr_data <- final_table %>% dplyr::select(, catchment_area_km, N, FIS, Ho.adj, PI_NEI, Private_Alleles) %>% 
  pivot_longer(cols=c(FIS, Ho.adj, PI_NEI, Private_Alleles))
library(broom)

catch_p_values <- corr_data %>%
    group_by(name) %>%
  do({
    model <- lm(value ~ catchment_area_km, data = .)
    tidied <- tidy(model)
    r2 <- summary(model)$r.squared
    tidied %>% filter(term == "catchment_area_km") %>%
      mutate(r.squared = r2)
  }) 
    # do(tidy(lm(value ~ catchment_area_km, data = .)))  %>%  # Fit the model and tidy the results
    # filter(term == "catchment_area_km")# %>%
    # select(name, p.value)  # Select the relevant columns
p_values 

ggplot(corr_data)+
  geom_point(aes(catchment_area_km, value))+
  # geom_smooth(aes(catchment_area_km, value), method = "lm", se = T, colour = "blue") +  # Add correlation line
  facet_wrap(~name, scales="free")+
  geom_text(data = catch_p_values, aes(x = Inf, y = Inf, label = paste("r2 =", round(r.squared, 3), ", p =", round(p.value, 3))), 
            hjust = 1.1, vjust = 1.1, size = 4, colour = "red", 
            position = position_nudge(x = 0.5, y = 0.5)) 

n_p_values <- corr_data %>%
  group_by(name) %>%
  do({
    model <- lm(value ~ N, data = .)
    tidied <- tidy(model)
    r2 <- summary(model)$r.squared
    tidied %>% filter(term == "N") %>%
      mutate(r.squared = r2)
  }) 
# do(tidy(lm(value ~ catchment_area_km, data = .)))  %>%  # Fit the model and tidy the results
# filter(term == "catchment_area_km")# %>%
# select(name, p.value)  # Select the relevant columns
n_p_values 

ggplot(corr_data)+
  geom_point(aes(N, value))+
  # geom_smooth(aes(catchment_area_km, value), method = "lm", se = T, colour = "blue") +  # Add correlation line
  facet_wrap(~name, scales="free")+
  geom_text(data = n_p_values, aes(x = Inf, y = Inf, label = paste("r2 =", round(r.squared, 3), ", p =", round(p.value, 3))), 
            hjust = 1.1, vjust = 1.1, size = 4, colour = "red", 
            position = position_nudge(x = 0.5, y = 0.5)) 

ggplot(final_table, aes(x = , y = PI_NEI)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7)


# Pivot the data to long format for easier plotting with ggplot2
final_table_long <- final_table %>%
  pivot_longer(cols = c(FIS, Ho.adj,He.adj, PI_NEI), names_to = "Metric", values_to = "Value")
  #pivot_longer(cols = c(CI_F, CI_Ho), names_to = "CI_Type", values_to = "CI") %>%
  # filter(Metric == substr(CI_Type, 4, nchar(CI_Type))) %>%  # Ensures matching CI columns
  # select(-CI_Type) # Drop the CI_Type column

# Plot
  library(paletteer)# for the classicpurple colours
statfig<-ggplot(final_table_long, aes(x = AnonSiteCode, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_paletteer_d("tvthemes::Pearl")+
  # scale_fill_paletteer_d("ggthemes::Purple_Pink_Gray")+
  # scale_fill_paletteer_d("ggthemes::Classic_Purple_Gray_12") +
  # geom_errorbar(aes(ymin = Value - CI, ymax = Value + CI), position = position_dodge(width = 0.9), width = 0.2) +
  labs(x = "Population", y = "Metric Value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
statfig
ggsave(plot= statfig,
       width = 6.5,
       height = 3,
       units = c( "in"),
       dpi = 600,
       filename = "../../output/figures/popgenstats.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/popgenstats.png")




# Rescale data for visualisation
library(scales)
data <- final_table %>%
  mutate(N_scaled = rescale(N),
         F_scaled = rescale(FIS),
         HO_scaled = rescale(Ho),
         Pi_scaled = rescale(PI_NEI),
         Private_Alleles_scaled = rescale(Private_Alleles))

  # Convert to long format for ggplot
data_long <- data %>%
  dplyr::select(, N_scaled, F_scaled, HO_scaled, Pi_scaled, Private_Alleles_scaled, CI_F, CI_Ho) %>%
  pivot_longer(cols = c(N_scaled, F_scaled, HO_scaled, Pi_scaled, Private_Alleles_scaled),
               names_to = "Metric", values_to = "Value") %>%
  mutate(CI = ifelse(Metric == "F_scaled", data$CI_F, 
                     ifelse(Metric == "HO_scaled", data$CI_Ho, 0))) # Set CI to 0 where not needed

# Plot
ggplot(data_long, aes(x = , y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_paletteer_d("tvthemes::Pearl")+
  geom_errorbar(aes(ymin = Value - CI, ymax = Value + CI),
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(title = "Comparison of Scaled Genetic Metrics across Populations",
       y = "Scaled Metric Values", x = "Population") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Investigate inbreeding/het
meanF<-read.csv( "../../data/confidential/anon_inbreeding_cooeff.csv")
head(meanF) 
vicec<-read.csv("../../data/confidential/trimmed_env_data.csv")%>%
  mutate(Sample_ID = ifelse(Sample_ID == "ICY837","ICY837notICY773", Sample_ID)) %>% mutate(Population = ifelse(Population == "Conon_Gharbhrain_wrong_sample_ID", "Conon_Glascarnoch", Population))%>% mutate(AnonSiteCode = ifelse(AnonSiteCode == "ER2", "ER1", AnonSiteCode)) 
head(vicec)
nrow(vicec)
length(unique(vicec$catchment_area_km))
merged_data <- het %>%
  inner_join(vicec, by = c("IID" = "Sample_ID"))
# Select specific columns to avoid duplicates
merged_data <- het %>%
  inner_join(vicec, by = c("IID" = "Sample_ID")) %>%
  dplyr::select(IID, O.HOM., E.HOM., N.NM., F, 
         everything(),  # This keeps all other columns
         -c(AnonSiteCode.y, X)) %>%   # Exclude the unwanted columns
  rename(AnonSiteCode = AnonSiteCode.x)  
head(merged_data)
nrow(merged_data)
length(unique(merged_data$catchment_area_km))

plot(merged_data$catchment_area_km, merged_data$F)
cor.test(merged_data$F, merged_data$catchment_area_km)

## Mean data
mean_data <- merged_data %>%
  group_by(catchment_area_km, AnonSiteCode) %>%
  summarise(mean_F = mean(F, na.rm = TRUE),
            se_F = sd(F, na.rm = TRUE) / sqrt(n()),
            .groups = 'drop')

# Update the Host variable to change "Only trout available" to "Trout"
merged_data <- merged_data %>%
  mutate(Host = recode(Host, "Only trout available" = "Trout"))

# Catchment area vs inbreeding coefficient
ggplot(merged_data, aes(x = catchment_area_km, y = F,
                        # colour = Host
                        )) +
  geom_point(size = 2, alpha = 0.4,) +  # Plot individual points
  # geom_errorbar(data = mean_data, aes(x = catchment_area_km, ymin = mean_F - se_F, ymax = mean_F + se_F), width = 0.2, colour = "green") +  
  geom_point(data = mean_data, aes(x = catchment_area_km, y = mean_F), size = 3, colour = "red") +  
  geom_text(data = mean_data, aes(x = catchment_area_km, y = mean_F, label = AnonSiteCode), 
            vjust = -0.5, size = 3.5, colour = "black") +  
  labs(x = "Catchment Area (km²)", y = "Inbreeding Coefficient (F)") +
  theme_classic()


# # Violin plot for F vs AnonSiteCode
ggplot(merged_data, aes(x = AnonSiteCode, y = F)) +
  geom_violin(trim = FALSE, fill = "lightblue") +
  geom_boxplot(width = 0.1, fill = "white", outlier.colour = NA) +  # Optional: Add a boxplot inside the violin
  labs(title = "Violin Plot of F vs AnonSiteCode",
       x = "AnonSite Code",
       y = "F Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Optional: Rotate x-axis labels for better visibility

ggplot(merged_data, aes(x = Host, y = F, fill = Host)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot
  geom_boxplot(width = 0.1, fill = "white", outlier.colour = NA) +  # Add boxplot inside the violin
  labs(x = "Host", y = "Inbreeding Coefficient (F)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility

# Violin plot for Heterozygosity vs AnonSiteCode
ggplot(merged_data, aes(x = AnonSiteCode, y = prob.obs.het)) +
  geom_violin(trim = FALSE, fill = "lightgreen") +
  geom_boxplot(width = 0.1, fill = "white", outlier.colour = NA) +  # Optional: Add a boxplot inside the violin
  labs(title = "Violin Plot of Heterozygosity vs AnonSiteCode",
       x = "AnonSite Code",
       y = expression(Observed~Heterozygosity~(H[obs]))) +  # Correctly formatted y-axis label
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Optional: Rotate x-axis labels for better visibility

mean_het_data <- merged_data %>%
  group_by(catchment_area_km, AnonSiteCode) %>%
  summarise(mean_obs_h = mean(prob.obs.het, na.rm = TRUE),
            se_obs_h = sd(prob.obs.het, na.rm = TRUE) / sqrt(n()),
            .groups = 'drop')
ggplot(merged_data, aes(x = catchment_area_km, y = prob.obs.het, colour = Host)) +
  geom_point(size = 2, alpha = 0.4) +  # Plot individual points
  # geom_errorbar(data = mean_data, aes(x = catchment_area_km, ymin = mean_F - se_F, ymax = mean_F + se_F), width = 0.2, colour = "green") +  
  geom_point(data = mean_het_data, aes(x = catchment_area_km, y = mean_obs_h), size = 3, colour = "red") +  
  geom_text(data = mean_het_data, aes(x = catchment_area_km, y = mean_obs_h, label = AnonSiteCode), 
            vjust = -0.5, size = 3.5, colour = "black") +  
  labs(x = "Catchment Area (km²)",
       y = expression(Observed~Heterozygosity~(H[obs]))) +  # Correctly formatted y-axis label
  theme_classic()
ggplot(merged_data, aes(x = Host, y = prob.obs.het, fill = Host)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot
  geom_boxplot(width = 0.1, fill = "white", outlier.colour = NA) +  # Add boxplot inside the violin
  labs(x = "Host", y = expression(Observed~Heterozygosity~(H[obs]))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility



## Stat test
# Subset data to include only Salmon and Trout
subset_data <- merged_data %>%
  filter(Host %in% c("Salmon", "Trout"))
# Check normality for F
shapiro.test(subset_data$F)
# Check normality for prob.obs.het
shapiro.test(subset_data$prob.obs.het)
# t-test for F
t_test_F <- t.test(F ~ Host, data = subset_data)
print(t_test_F)
# t-test for prob.obs.het
t_test_obs_het <- t.test(prob.obs.het ~ Host, data = subset_data)
print(t_test_obs_het)
# Wilcoxon test for F
wilcox_test_F <- wilcox.test(F ~ Host, data = subset_data)
print(wilcox_test_F)

# Wilcoxon test for prob.obs.het
wilcox_test_obs_het <- wilcox.test(prob.obs.het ~ Host, data = subset_data)
print(wilcox_test_obs_het)
# Load necessary libraries
library(dplyr)

# Define the function
test_significance <- function(data, group_var, var1, var2) {
  # Subset data to include only the specified groups
  subset_data <- data %>%
    filter(!!sym(group_var) %in% c("Salmon", "Trout"))
  
  # Check normality for var1
  shapiro_test_var1 <- shapiro.test(subset_data[[var1]])
  print(paste("Shapiro-Wilk test for", var1, ":", shapiro_test_var1$p.value))
  
  # Check normality for var2
  shapiro_test_var2 <- shapiro.test(subset_data[[var2]])
  print(paste("Shapiro-Wilk test for", var2, ":", shapiro_test_var2$p.value))
  
  # Function to perform the appropriate test
  perform_test <- function(var, group) {
    if (shapiro.test(subset_data[[var]])$p.value > 0.05) {
      # If normally distributed, use t-test
      return(t.test(subset_data[[var]] ~ subset_data[[group]]))
    } else {
      # If not normally distributed, use Wilcoxon test
      return(wilcox.test(subset_data[[var]] ~ subset_data[[group]]))
    }
  }
  
  # Perform tests for both variables
  result_var1 <- perform_test(var1, group_var)
  result_var2 <- perform_test(var2, group_var)
  
  # Return the results
  return(list(
    var1_result = result_var1,
    var2_result = result_var2
  ))
}

# Call the function
results <- test_significance(merged_data, "Host", "F", "prob.obs.het")

# Print results
print(results$var1_result)
print(results$var2_result)


######
ggplot(merged_data, aes(x = Host, y = catchment_area_km, fill = Host)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Alpha controls transparency
  geom_boxplot(width = 0.1, fill = "white", outlier.colour = NA) +  # Optional: boxplot inside violin
  labs(x = "Host", y = "Catchment Area (km²)", title = "Violin Plot of Catchment Area by Host") +
  theme_classic()

ggplot(merged_data, aes(x = AnonSiteCode, y = catchment_area_km, fill = Host)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot, with transparency control
  labs(x = "Site", y = "Catchment Area (km²)", title = "Violin Plot of Catchment Area by Site and Host") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility
ggplot(merged_data, aes(x = AnonSiteCode, y = catchment_area_km, colour = Host)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot points
  labs(x = "Site", y = "Catchment Area (km²)", title = "Catchment Area by Site and Host") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
merged_data %>%
  group_by(catchment_area_km) %>%  # Group by catchment area
  summarise(n_samples = n(),  # Count the number of sites
            sites = paste(unique(AnonSiteCode), collapse = ", "))  # Optional: list of sites in each catchment
catchment_summary <- merged_data %>%
  group_by(catchment_area_km) %>%  # Group by catchment area
  summarise(
    n_samples = n(),  # Count the number of samples
    sites = paste(unique(AnonSiteCode), collapse = ", ")  # List of sites in each catchment
  ) %>%
  ungroup()
# Create a new dataset with the Catchment_Label as the site names
merged_data_r <- merged_data %>%
  left_join(catchment_summary, by = "catchment_area_km") %>%  # Join to bring in sites info
  mutate(Catchment_Label = sites)  # Assign catchment label as the site names

## Create new catchment name column to group tose in the same catchment
# merged_data_r <- merged_data %>%
  # mutate(Catchment_Label = paste0("Catchment_", dense_rank(catchment_area_km)))  # Assign catchment labels without grouping
merged_data_r %>%
  group_by(catchment_area_km) %>%  # Group by catchment area
  summarise(n_samples = n(),  # Count the number of sites
            sites = paste(unique(Catchment_Label), collapse = ", "))  # Optional: list of sites in each catchment
ggplot(merged_data_r, aes(x = Catchment_Label, y = catchment_area_km, colour = Host)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot points
  labs(x = "Site", y = "Catchment Area (km²)", title = "Catchment Area by Site and Host") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
ggplot(merged_data_r, aes(x = Host, y = catchment_area_km, colour = Host)) +
  geom_point(size = 3, alpha = 0.7) +  # Plot points for each catchment area
  labs(x = "Host", y = "Catchment Area (km²)", title = "Catchment Area by Host") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
ggplot(merged_data_r, aes(x = Host, y = catchment_area_km, colour = Host)) +
  geom_point(size = 3, alpha = 0.7) +  # Plot points for each catchment area
  labs(x = "Host", y = "Catchment Area (km²)", title = "Catchment Area by Host") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

site_group <- merged_data_r %>%
  group_by(AnonSiteCode) %>%
  summarise(
    catchment_area_km = mean(catchment_area_km),  # or any other summary statistic
    Host = first(Host),
    Catchment_Label = first(Catchment_Label) # Assuming Host is consistent for each AnonSiteCode
  ) %>%
  ungroup()

head(site_group)  
ggplot(site_group, aes(x = Host, y = catchment_area_km, label = Catchment_Label, colour = Host)) +
  geom_text(size = 3, alpha = 0.7) +  # Use geom_text to plot AnonSiteCode
  labs(x = "Host", y = "Catchment Area (km²)", title = "Catchment Area by Host") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

lm<-lm(merged_data_r$F~merged_data_r$Host+merged_data_r$catchment_area_km+merged_data_r$AnonSiteCode)
summary(lm)
plot(lm)
ggplot(merged_data_r, aes(x = Host, y = F, colour = Host)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add points for individual samples
  labs(title = "Inbreeding Coefficient (F) by Host Type",
       x = "Host Type",
       y = "Inbreeding Coefficient (F)") +
  theme_classic()
