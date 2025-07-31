### Env data trimming and selection is seperate script! Should go straight to the RDA and popgen from here
library(dplyr)
library(ggplot2)
library(tidyverse)
# Labels for graphing later!
labels<-read.csv("data/too_large_files/all_freshwater_sdmpredictors/layerschoice_metadata.csv")
labels<-labels[,2:3]
legend_labels <- setNames(labels$name, labels$layer_code)

# env<-read.csv("data/confidential/fw_and_other_env_data.csv")
env<-read.csv("data/confidential/anon_all_fw_and_other_env_data.csv")
str(env)
colnames(env) <- sub("_lonlat$", "", colnames(env))
colnames(env)
# Remove uneccerasary NA values
sum(is.na(env)) #47 from Host col, lenght and width

env %>%
  filter(is.na(Host)) %>%
  dplyr::select(Sample_ID, Population) #the wrong sampleID is a problem because I gave it lat/long too

sum(is.na(env))
which(is.na(env$Host), arr.ind = TRUE)
unique(env$Host)
env$Host[env$Population == "Conon_Gharbhrain_wrong_sample_ID"] <- "Unknown"
which(is.na(env$Host), arr.ind = TRUE)

#####

##### A bit of tidying up ####
## Lets look at the data and remove useless columns

# Find the names of the columns that are entirely zero
cols_to_remove <- env %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(where(~ sum(.) == 0)) %>%
  colnames()
cols_to_remove <- env %>%
  dplyr::select(where(is.numeric)) %>%              # Select numeric columns
  summarise(across(everything(), ~ all(. == 0))) %>%  # Check if all values in each column are zero
  pivot_longer(everything(), names_to = "column", values_to = "all_zero") %>%  # Convert to long format
  filter(all_zero) %>%                       # Filter columns where all values are zero
  pull(column)                              # Extract column names

# Print the names of the columns to be removed
print(cols_to_remove) 
# "FW_geo_wsum_01_lonlat" "FW_geo_wsum_02_lonlat" "FW_geo_wsum_03_lonlat" "FW_geo_wsum_04_lonlat"
# [5] "FW_geo_wsum_05_lonlat" "FW_geo_wsum_06_lonlat" "FW_geo_wsum_07_lonlat" "FW_geo_wsum_08_lonlat"
# [9] "FW_geo_wsum_09_lonlat" "FW_geo_wsum_11_lonlat" "FW_geo_wsum_12_lonlat" "FW_geo_wsum_13_lonlat"
# [13] "FW_geo_wsum_14_lonlat" "FW_geo_wsum_15_lonlat" "FW_geo_wsum_16_lonlat" "FW_geo_wsum_17_lonlat"
# [17] "FW_geo_wsum_18_lonlat" "FW_geo_wsum_19_lonlat" "FW_geo_wsum_20_lonlat" "FW_geo_wsum_21_lonlat"
# [21] "FW_geo_wsum_22_lonlat" "FW_geo_wsum_23_lonlat" "FW_geo_wsum_24_lonlat" "FW_geo_wsum_25_lonlat"
# [25] "FW_geo_wsum_26_lonlat" "FW_geo_wsum_27_lonlat" "FW_geo_wsum_29_lonlat" "FW_geo_wsum_30_lonlat"
# [29] "FW_geo_wsum_31_lonlat" "FW_geo_wsum_32_lonlat" "FW_geo_wsum_33_lonlat" "FW_geo_wsum_34_lonlat"
# [33] "FW_geo_wsum_35_lonlat" "FW_geo_wsum_36_lonlat" "FW_geo_wsum_37_lonlat" "FW_geo_wsum_38_lonlat"
# [37] "FW_geo_wsum_39_lonlat" "FW_geo_wsum_40_lonlat" "FW_geo_wsum_41_lonlat" "FW_geo_wsum_42_lonlat"
# [41] "FW_geo_wsum_43_lonlat" "FW_geo_wsum_44_lonlat" "FW_geo_wsum_45_lonlat" "FW_geo_wsum_46_lonlat"
# [45] "FW_geo_wsum_47_lonlat" "FW_geo_wsum_48_lonlat" "FW_geo_wsum_49_lonlat" "FW_geo_wsum_50_lonlat"
# [49] "FW_geo_wsum_51_lonlat" "FW_geo_wsum_53_lonlat" "FW_geo_wsum_54_lonlat" "FW_geo_wsum_55_lonlat"
# [53] "FW_geo_wsum_56_lonlat" "FW_geo_wsum_57_lonlat" "FW_geo_wsum_58_lonlat" "FW_geo_wsum_59_lonlat"
# [57] "FW_geo_wsum_60_lonlat" "FW_geo_wsum_61_lonlat" "FW_geo_wsum_62_lonlat" "FW_geo_wsum_63_lonlat"
# [61] "FW_geo_wsum_64_lonlat" "FW_geo_wsum_65_lonlat" "FW_geo_wsum_66_lonlat" "FW_geo_wsum_68_lonlat"
# [65] "FW_geo_wsum_69_lonlat" "FW_geo_wsum_70_lonlat" "FW_geo_wsum_71_lonlat" "FW_geo_wsum_72_lonlat"
# [69] "FW_geo_wsum_73_lonlat" "FW_geo_wsum_74_lonlat" "FW_geo_wsum_75_lonlat" "FW_geo_wsum_77_lonlat"
# [73] "FW_geo_wsum_78_lonlat" "FW_geo_wsum_79_lonlat" "FW_geo_wsum_81_lonlat" "FW_geo_wsum_82_lonlat"
# [77] "FW_geo_wsum_83_lonlat" "FW_geo_wsum_84_lonlat" "FW_geo_wsum_85_lonlat" "FW_geo_wsum_86_lonlat"
# [81] "FW_geo_wsum_87_lonlat" "FW_geo_wsum_88_lonlat" "FW_geo_wsum_89_lonlat" "FW_geo_wsum_90_lonlat"
# [85] "FW_geo_wsum_91_lonlat" "FW_geo_wsum_92_lonlat" "FW_lc_avg_02_lonlat"   "FW_lc_avg_09_lonlat"  
# [89] "FW_lc_avg_10_lonlat"   "FW_lc_max_02_lonlat"   "FW_lc_max_09_lonlat"   "FW_lc_max_10_lonlat"  
# [93] "FW_lc_min_01_lonlat"   "FW_lc_min_02_lonlat"   "FW_lc_min_03_lonlat"   "FW_lc_min_08_lonlat"  
# [97] "FW_lc_min_09_lonlat"   "FW_lc_min_10_lonlat"   "FW_lc_min_11_lonlat"   "FW_lc_min_12_lonlat"  
# [101] "FW_lc_range_02_lonlat" "FW_lc_range_09_lonlat" "FW_lc_range_10_lonlat" "FW_lc_wavg_02_lonlat" 
# [105] "FW_lc_wavg_09_lonlat"  "FW_lc_wavg_10_lonlat" 
summary(env%>%dplyr::select(all_of(cols_to_remove)))

# Remove the columns that are entirely zero# Remove the columns that are entirely zerorow.names = 
env <- env %>%
  dplyr::select((-all_of(cols_to_remove)))

# Next remove unnecasary metadata columns
env <- env%>%   dplyr::select(-Length.mm.,-Width.mm.,-Life_Stage, -Barcode, -ID, -VCNUMBER, -EXTENTS)

# ## Leta attACH THE worldclim dataset
# wenv<-read.csv("data/confidential/swab_enviroment.csv")

# Merge wenv$WC_alt_lonlat into env based on Sample_ID
# env <- env %>% 
#   left_join(select(wenv, Sample_ID, WC_alt_lonlat), by = "Sample_ID")%>% 
#   dplyr::select(-Host, everything(), Host)
# 
# env <- env%>%
#   left_join(wenv %>% dplyr::select(Sample_ID, WC_alt_lonlat), by = "Sample_ID")%>%
#   dplyr::select(-Host, Host)


## Function for correlation testing/vis
## Function to text correlations
library(corrplot)

correlations <- function(data) {
  # Remove specific columns (VCNUMBER, Latitude, Longitude)
  columns_to_remove <- c("VCNUMBER", "Latitude", "Longitude")
  selected_cols <- data[, !names(data) %in% columns_to_remove]
  
  # Select numeric and integer columns
  numeric_cols <- sapply(selected_cols, is.numeric)
  integer_cols <- sapply(selected_cols, is.integer)
  selected_cols <- selected_cols[, numeric_cols | integer_cols]
  
  # Convert selected columns to numeric if needed
  selected_cols <- lapply(selected_cols, as.numeric)
  
  # Create correlation matrix
  data_matrix <- as.matrix(do.call(cbind, selected_cols))
  data_df <- as.data.frame(data_matrix)
  
  #legend labels
  # Set legend labels
  columns_without_layer_code <- setdiff(colnames(data_df), labels$layer_code)
  # Create labels data frame with new entries
  labels <- data.frame(layer_code = c(labels$layer_code, columns_without_layer_code),
                       name = c(labels$name, columns_without_layer_code))
  setdiff(colnames(data_df), labels$layer_code)
  
  colnames(data_df) <-  setNames(labels$name[labels$layer_code %in% colnames(data_df)], labels$layer_code[labels$layer_code %in% colnames(data_df)])
  ###
  
  # Calculate Spearman's correlation
  cor_spear <- cor(x = data_df,
                   use = "pairwise.complete",
                   method = "spearman") #pearson?
  
  # Set diagonal elements to NA
  diag(cor_spear) <- NA
  
  # Plot correlation matrix using corrplot
  corrplot(cor_spear, type = "lower", method = "number", 
           tl.col = 'black', tl.srt = 45, tl.cex = 0.7,number.cex = 0.9)
  apply(abs(cor_spear), MARGIN = 2, FUN = max, na.rm = TRUE) # everything is super correlated other than bio3
  # Identify the maximum absolute correlation for each variable
  max_cor <- apply(abs(cor_spear), MARGIN = 2, FUN = max, na.rm = TRUE)
  
  # Print the maximum correlations
  print(max_cor)
}

 # correlations(env)

# Removing repeated data groups (max, min, av, wavg, range, precipitation and temp variables)
env_red <- env %>% #this is original before I started fiddling for slope
  dplyr::select(
    -matches("WC_bio"), #remove terrestriral params
    -matches("X.x"),
    -matches("FW_geo_wsum"), #underlying geology weighted sum
    -matches("FW_tmin"), #temp min
    -matches("FW_tmax"), #temp max
    # -matches("slope_"), #removes slope
    -matches("dem_av"), #this leaves elevation range
    -matches("dem_m"),
    -matches("flow_len"),
    -matches("FW_lc_m"), #this is land cover- I dont like the model accuracy?
    -matches("FW_lc_range"),
    -matches("FW_lc_av"),
    -matches("FW_lc"), #**** I have removed all landcover stuff
    -matches("hydro_avg"), #leaves bioclim wavg
    -matches("prec_sum"),
    -matches("prec_wsum"),
    -matches("soil_avg"), #leaves soil weighted avg
    -matches("soil_range"),
    -matches("soil_m"),
     -ID.1,
    -COUNT_S_ID, -Shape_Leng, -Shape_Area
  )
# env_red <- env %>% ###aaah I want to put slope in now 17.02.2025 but I probably dont have time
#   dplyr::select(
#     -matches("X.x"),
#     -matches("FW_geo_wsum"), #underlying geology weighted sum
#     -matches("FW_tmin"), #temp min
#     -matches("FW_tmax"), #temp max
#     # -matches("slope_"), #removes slope
#     -matches("dem_av"), #this leaves elevation range
#     -matches("dem_m"),
#     -matches("flow_len"),
#     -matches("FW_lc_m"), #this is land cover- I dont like the model accuracy? 
#     -matches("FW_lc_range"),
#     -matches("FW_lc_av"),
#     -matches("FW_lc"), #**** I have removed all landcover stuff
#     -matches("hydro_avg"), #leaves bioclim wavg
#     -matches("prec_sum"),
#     -matches("prec_wsum"),
#     -matches("soil_avg"), #leaves soil weighted avg
#     -matches("soil_range"),
#     -matches("soil_m")
#   )

colnames(env_red)
correlations(env_red)


# #================= testing #########
# 
# # Remove specific columns (VCNUMBER, Latitude, Longitude)
# columns_to_remove <- c("VCNUMBER", "Latitude", "Longitude")
# selected_cols <- env[, !names(env) %in% columns_to_remove]
# 
# # Select numeric and integer columns
# numeric_cols <- sapply(selected_cols, is.numeric)
# integer_cols <- sapply(selected_cols, is.integer)
# selected_cols <- selected_cols[, numeric_cols | integer_cols]
# 
# # Convert selected columns to numeric if needed
# selected_cols <- lapply(selected_cols, as.numeric)
# 
# # Create correlation matrix
# data_matrix <- as.matrix(do.call(cbind, selected_cols))
# data_df <- as.data.frame(data_matrix)
# 
# # Set legend labels
# setdiff(colnames(data_df), labels$layer_code)
# # Assuming 'data_df' is your data frame and 'labels' is your labels data frame
# columns_without_layer_code <- setdiff(colnames(data_df), labels$layer_code)
# # Create labels2 data frame with new entries
# labels <- data.frame(layer_code = c(labels$layer_code, columns_without_layer_code),
#                       name = c(labels$name, columns_without_layer_code))
# setdiff(colnames(data_df), labels2$layer_code)
# 
# colnames(data_df) <-  setNames(labels$name[labels$layer_code %in% colnames(data_df)], labels$layer_code[labels$layer_code %in% colnames(data_df)])
# colnames(data_df)
# ###
# 
# ###
# # Calculate Spearman's correlation
# cor_spear <- cor(data_df, use = "pairwise.complete.obs", method = "spearman")
# 
# # Set diagonal elements to NA
# diag(cor_spear) <- NA
# 
# # Plot correlation matrix using corrplot
# corrplot(cor_spear, type = "lower", method = "number", 
#          tl.col = 'black', tl.srt = 45, tl.cex = 0.5, number.cex = 0.5)
# 
# # Identify the maximum absolute correlation for each variable
# max_cor <- apply(abs(cor_spear), MARGIN = 2, FUN = max, na.rm = TRUE)
# 
# # Print the maximum correlations
# print(max_cor)
# 
# #================= ############

env_red_biocl <- env_red %>% # trimming bioclim factors
  dplyr::select(
    -FW_hydro_wavg_01, -FW_hydro_wavg_02,-FW_hydro_wavg_09,
    -FW_hydro_wavg_04, -FW_hydro_wavg_06, -FW_hydro_wavg_07,
    -FW_hydro_wavg_08,  -FW_hydro_wavg_12, -FW_hydro_wavg_11,
    -FW_hydro_wavg_13, -FW_hydro_wavg_14, -FW_hydro_wavg_15,
    -FW_hydro_wavg_16, -FW_hydro_wavg_18, -FW_hydro_wavg_19
  )
correlations(env_red_biocl)

env_finaltrim<-env_red_biocl %>%
  dplyr::select((-matches("soil_wavg_03"))) %>% #trim out other high corr env factors
  dplyr::select((-matches("NEW_P_IN"))) %>%
  # dplyr::select((-matches("pH_W_Med"))) %>% #basically exact same as organic carbon
  dplyr::select((-matches("FW_lc_wavg_01"))) %>%
  dplyr::select((-matches("AWC"))) %>%
  dplyr::select((-matches("soil_wavg_05")))%>%
  dplyr::select((-matches("soil_wavg_06")))%>%
  dplyr::select((-matches("soil_wavg_08")))%>%
  dplyr::select((-matches("soil_wavg_09")))%>%
  dplyr::select((-matches("soil_wavg_10")))%>%
  dplyr::select((-matches("FW_lc_wavg_12")))

if (!file.exists("output/figures")) dir.create("output/figures")
png("output/figures/correlation_plot.png", width = 1800, height = 1800, res = 150)
correlations(env_finaltrim)
dev.off()
# Check the result
colnames(env_finaltrim)
labels$name[match(colnames(env_finaltrim), labels$layer_code)]
summary(env_finaltrim)

colnames(env_finaltrim)
png("output/figures/correlation_plot_just_RDA.png", width = 6, height = 6, units = "in", res = 150)
correlations(env_finaltrim %>% dplyr::select(FW_dem_range,FW_hydro_wavg_03,FW_hydro_wavg_05,FW_hydro_wavg_17,wildness, WC_alt, catchment_area_km))
dev.off()

##########################################################################################
write.csv(env_finaltrim, "data/confidential/trimmed_env_data.csv", row.names = F)
# trim<-read.csv("data/confidential/trimmed_env_data.csv")
##########################################################################################




# Consolidate the veg columns
# env_finaltrim$VegetationMeans <- rowMeans(env_finaltrim[, c("FW_lc_wavg_03", "FW_lc_wavg_04", "FW_lc_wavg_05",
#                                      "FW_lc_wavg_06", "FW_lc_wavg_07", "FW_lc_wavg_08", "FW_lc_wavg_11")])

## Lets plot the remaining variables for variation
library(tidyr)
library(ggplot2)
# Normalggplot2# Normalize the numeric columns
just_env <- data.frame(env_finaltrim[, 7:ncol(env_finaltrim)])
just_env<-as.data.frame(just_env)
env_normalized <- as.data.frame(just_env) %>%
  dplyr::select(where(is.numeric)) %>%
  scale() %>%
  as.data.frame()
ncol(just_env)
# Convert to long format
env_long <- env_normalized %>%
  pivot_longer(cols = everything(), names_to = "Column", values_to = "Value")
# Plotting boxplot
ggplot(env_long, aes(x = Column, y = Value)) +
  geom_boxplot(fill = "skyblue") +
  geom_jitter(width = 0.2, alpha = 0.2, colour = "deepskyblue4") +  # Add jittered points
  labs(x = "Column", y = "Normalized Value", title = "Normalized Range of Numeric Columns Across Samples")+
  coord_flip()+
  theme(axis.text.y = element_text(size = 6))  # Adjust the size as needed
env_long <- as.data.frame(just_env) %>%
  pivot_longer(cols = everything(), names_to = "Column", values_to = "Value")
ggplot(env_long, aes(x = Column, y = Value)) +
  geom_boxplot(fill = "skyblue") +
  geom_jitter(width = 0.2, alpha = 0.2, colour = "deepskyblue4") +  # Add jittered points
  labs(x = "Column", y = "Value", title = "Normalized Range of Numeric Columns Across Samples")+
  facet_wrap(~ Column, scales = "free") +  # Facet by column with independent y scales
  coord_flip()+
  theme(axis.text.y = element_text(size = 6))  # Adjust the size as needed
######


# Summarise the data by population (mean, median, or another stat)
env_summary <- env %>%
  group_by(Population) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))  # Replace 'mean' with 'median' if needed

# Convert summarised dataframe to long format for plotting
env_long <- env_summary %>%
  pivot_longer(cols = -Population, names_to = "Column", values_to = "Value")  # Exclude 'pop' from pivoting

# Plot boxplots grouped by population, facet by numeric columns
ggplot(env_long, aes(x = Population, y = Value)) +
  geom_boxplot(fill = "skyblue") +
  geom_jitter(width = 0.2, alpha = 0.6, colour = "darkblue") +  # Add jittered points for visualisation
  labs(x = "Population", y = "Value", title = "Summary of Numeric Columns by Population") +
  facet_wrap(~ Column, scales = "free") +  # Facet by numeric columns with independent y scales
  coord_flip() +
  theme(axis.text.y = element_text(size = 6))



##Base R version ####
# # Confirm that genotypes and environmental data are in the same order
# identical(rownames(gen.imp), env$Sample_ID) # FALSE
# #Get rownames of gen.imp 
# rownames(gen.imp)
# # Get a vector to instigate "matching"
# matching_indices<-match(rownames(gen.imp), env$Sample_ID)
# # Filter out env data where no genomic data (aka the samples that where filtered out in the bioinformatic cleaning step)
# env <- env[matching_indices, ]
# # Are any barcodes redundant?
# unused_samples<-env%>% filter(is.na(Sample_ID))
# # Check if now the genotypes and environmental data are in the same order
# identical(rownames(gen.imp), env$Sample_ID) # FALSE
###########

# Removing correlating env variables as can skew analysis (aka multicolinearity)
# pairs.panels(env[,10:29], scale=T) #library(psych)

# Check for collinearity among predictors #####
# Pairwise correlations - |rho| < 0.70
just_env <- data.frame(env_finaltrim[, 7:(ncol(env_finaltrim))])
# pairs.panels(just_env, scale=T, hist.border = NA, stars = T) #library(psych)

just_env <- lapply(just_env, as.numeric)  # Identify numeric columns
just_env_matrix <- as.matrix(do.call(cbind, just_env))
just_env_df <- as.data.frame(just_env_matrix)

cor_spear <- cor(x = just_env_df,
                 use = "pairwise.complete",
                 method = "spearman") # Spearman's Correlation

diag(cor_spear) <- NA
# install.packages("corrplot")
# library(corrplot)
corrplot(cor_spear, type="lower",method="number", tl.col = 'black',tl.srt = 45, tl.cex=0.5)
apply(abs(cor_spear), MARGIN = 2, FUN = max, na.rm = TRUE) # keeping it belo 0.8
# Identify the maximum absolute correlation for each variable
max_cor <- apply(abs(cor_spear), MARGIN = 2, FUN = max, na.rm = TRUE)

# Print the maximum correlations
print(max_cor)
summary(env)
# Identify variables with maximum correlation less than 0.70
low_cor_vars<-names(max_cor[max_cor < 0.90])
length(low_cor_vars)
length(max_cor)
# just_env <- data.frame(env[, 5:(ncol(env))])
# low_cor_vars <- just_env_df[, low_cor_vars, drop = FALSE]

# filtered_env2 <- filtered_env[, 1:26]
# filtered_env2 <- lapply(filtered_env2, as.numeric)  # Identify numeric columns
# just_env_matrix <- as.matrix(do.call(cbind, filtered_env2))
# just_env_df <- as.data.frame(just_env_matrix)
# 
# cor_spear <- cor(x = just_env_df,
#                  use = "pairwise.complete",
#                  method = "spearman") # Spearman's Correlation
# 
# diag(cor_spear) <- NA
# apply(abs(cor_spear), MARGIN = 2, FUN = max, na.rm = TRUE) # everything is super correlated other than bio3
# # Identify the maximum absolute correlation for each variable
# max_cor <- apply(abs(cor_spear), MARGIN = 2, FUN = max, na.rm = TRUE)
# 
# # Print the maximum correlations
# print(max_cor)






######

#=====================I just want final plot of env in RDA correlations #######################
colnames(env)
corr_cols<-env %>% dplyr::select(FW_dem_range, FW_hydro_wavg_03 , FW_hydro_wavg_05, FW_hydro_wavg_17, wildness, catchment_area_km, WC_alt_lonlat)
correlations(corr_cols)

#========================= A simple printable map of loc and env ========================================
# 
library(mapview)
library(sf)
geo_env <- env[!duplicated(env$Population), ]
geo_env <- geo_env[!is.na(geo_env$Latitude),]
geo_env<-st_as_sf(geo_env, coords = c("Longitude", "Latitude"), crs = "epsg:4326")
# geo_env$Population <- factor(geo_env$Population, levels = sort(unique(geo_env$Population)))

mapview(geo_env, zcol="Population")
mapview(geo_env, zcol = "WC_bio9_lonlat")
mapview(geo_env, zcol = "WC_bio10_lonlat")
mapview(geo_env, zcol = "WC_bio17_lonlat")
mapview(geo_env, zcol = "WC_bio2_lonlat")
mapview(geo_env, zcol = "WC_bio3_lonlat")
mapview(geo_env, zcol = "WC_bio5_lonlat")



plot(env$Population, env$WC_alt_lonlat)
summary(env)