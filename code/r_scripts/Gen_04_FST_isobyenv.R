## Isolation by enviroment
#do FST code first
# Assuming `env` has Population in column 1 and env variables in the rest
trimmed_gen_pops<-read.csv("../../data/confidential/meta_for_filtered_genomes.csv") %>% dplyr::select(-X) %>% mutate(Population = ifelse(Population == "Conon_Gharbhrain_wrong_sample_ID", "Conon_Glascarnoch", Population)) %>% mutate(AnonSiteCode = ifelse(AnonSiteCode == "ER2", "ER1", AnonSiteCode)) %>% dplyr::select(AnonSiteCode, VCNAME, Population, ID)
trimmed_gen_pops
colnames(vicec)
env <- vicec %>%
  filter(Sample_ID %in% trimmed_gen_pops$ID) %>%
  dplyr::select(Sample_ID, AnonSiteCode, Population, FW_dem_range, wildness, WC_alt, 
         FW_hydro_wavg_03, FW_hydro_wavg_05, FW_hydro_wavg_17, catchment_area_km) %>%
  group_by(Population) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

env
env_pca <- prcomp(env %>% 
                    dplyr::select(FW_dem_range, wildness, WC_alt, FW_hydro_wavg_03, FW_hydro_wavg_05, FW_hydro_wavg_17, catchment_area_km),
                  scale. = TRUE)
env$PC1 <- env_pca$x[, 1]  # Extract PC1
colnames(env)
library(proxy)
# Create pairwise matrix
env_dist_matrix <- as.matrix(dist(env$PC1))
rownames(env_dist_matrix) <- env$Population
colnames(env_dist_matrix) <- env$Population
env_dist_matrix
# Reorder env_dist_matrix to match FST matrix
fst_matrix_lin
env_dist_matrix <- env_dist_matrix[rownames(fst_matrix_lin), colnames(fst_matrix_lin)]
mantel_env <- vegan::mantel(fst_matrix_lin, env_dist_matrix, method = "spearman", permutations = 9999)
mantel_env #r=−0.1884 with a significance value of p=0.8819


library(reshape2)
env_dist_long <- melt(env_dist_matrix, varnames = c("Var1", "Var2"), value.name = "Env_Dist")
ibe_data <- merged_data %>%
  left_join(env_dist_long, by = c("Var1", "Var2")) %>%
  na.omit()

pc1_var <- summary(env_pca)$importance[2, 1] * 100
round(pc1_var, 3)

ggplot(ibe_data, aes(x = Env_Dist, y = Linearised_Fst)) +
  geom_point(size = 3, alpha = 0.6, colour = "darkgreen") +
  geom_smooth(method = "lm", se = FALSE, colour = "darkgreen") +
  theme_classic()+
  labs(
  x = paste0("Environmental Distance (PC1 = ", round(pc1_var, 2), "%)"),
  y = expression(paste(F[ST], " / (1 - ", F[ST], ")"))
)

