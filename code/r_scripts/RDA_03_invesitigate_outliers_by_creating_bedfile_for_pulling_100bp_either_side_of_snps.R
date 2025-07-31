# Investigate the outlier SNPs

outliers<-read.csv("output/rda_qvalue_outlier_snps.csv")

# Example: your SNP vector
snp_names <- outliers$SNP
head(outliers)
# Convert to BED format
snp_df <- do.call(rbind, strsplit(snp_names, "_"))
colnames(snp_df) <- c("contig", "position", "allele")
snp_df <- as.data.frame(snp_df)
snp_df$position <- as.numeric(as.character(snp_df$position))

# BED is 0-based, half-open: start = pos - 101, end = pos + 100
snp_df$start <- snp_df$position - 101
snp_df$end <- snp_df$position + 100
snp_df$name <- paste0(snp_df$contig, "_", snp_df$position, "_", snp_df$allele)

# Reorder to BED columns
bed_df <- snp_df[, c("contig", "start", "end", "name")]

# check contig names in the ref genome: grep "^>" GCA_029931535.1_MarmarV2_genomic.fna | cut -c 2- > ref_contigs.txt

# Extract contig number from scaf name
snp_df$contig_num <- sprintf("%04d", as.numeric(gsub("scaf", "", snp_df$contig)))

# Build the FASTA contig name (up to the first space)
snp_df$fasta_contig <- paste0("JAQPZY01000", snp_df$contig_num, ".1")

bed_df <- snp_df[, c("fasta_contig", "start", "end", "name")]

# Write BED file
write.table(bed_df, file = "output/rda_qvalue_outliers_snp_flanks_100bp.bed", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
## 1000 bp
# BED with ±1000 bp flanks
snp_df$start_1000 <- snp_df$position - 1001
snp_df$end_1000 <- snp_df$position + 1000

bed_df_1000 <- data.frame(
  fasta_contig = snp_df$fasta_contig,
  start = snp_df$start_1000,
  end = snp_df$end_1000,
  name = snp_df$name
)

# Write the new BED file
write.table(bed_df_1000, file = "output/rda_qvalue_outliers_snp_flanks_1000bp.bed", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# BED with ±500 bp flanks
snp_df$start_500 <- snp_df$position - 501
snp_df$end_500 <- snp_df$position + 500

bed_df_500 <- data.frame(
  fasta_contig = snp_df$fasta_contig,
  start = snp_df$start_500,
  end = snp_df$end_500,
  name = snp_df$name
)

# Write the new BED file
write.table(bed_df_500, file = "output/rda_qvalue_outliers_snp_flanks_500bp.bed", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

##########################################################################
### Now head onto the HPC and run Unix code 5_investigating_outlier_SNPs
##########################################################################
load("data/confidential/RDA_host_condition_updated.RData")

### Import and parse Blast results
### Load libraries
library(xml2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(purrr)


###3rd attemp
# ---- Helper Function: Parse and Filter BLAST XML ----
parse_blast_hits <- function(xml_path) {
  blast_xml <- read_xml(xml_path)
  iterations <- xml_find_all(blast_xml, ".//Iteration")
  
  extract_hits <- function(iter) {
    query_id <- xml_text(xml_find_first(iter, ".//Iteration_query-def"))
    hits <- xml_find_all(iter, ".//Hit")
    
    if (length(hits) == 0) return(NULL)
    
    do.call(rbind, lapply(hits, function(hit) {
      data.frame(
        Query_ID    = query_id,
        Subject_ID  = xml_text(xml_find_first(hit, ".//Hit_id")),
        Description = xml_text(xml_find_first(hit, ".//Hit_def")),
        E_value     = as.numeric(xml_text(xml_find_first(hit, ".//Hsp_evalue"))),
        stringsAsFactors = FALSE
      )
    }))
  }
  
  xml_hits <- map_dfr(iterations, extract_hits)
  
  xml_hits %>%
    filter(E_value <= 1e-5) %>%
    mutate(SNP = sub("::.*", "", Query_ID)) %>%
    filter(!grepl("genome assembly|microsatellite|chromosome|sequence", Description, ignore.case = TRUE)) #filter out general NCBI uploads
}

regions <- c("100bp", "500bp", "1000bp")
plots <- list()
combined_hits <- list()

for (region in regions) {
  # --- Load BLAST hits ---
  full_hits <- parse_blast_hits(file.path("output/blast_results/", paste0("blast_results_", region, "_full.xml"))) %>%
    mutate(Database = "full", Proximity = region)
  
  mollusc_hits <- parse_blast_hits(file.path("output/blast_results/", paste0("blast_results_", region, "_mollusc.xml"))) %>%
    mutate(Database = "molluscan", Proximity = region)
  
  # --- Combine hits for export ---
  combined_hits[[region]] <- bind_rows(full_hits, mollusc_hits)
  all_blast_hits <- bind_rows(combined_hits)
  
  # 
  # # --- Get top full hits and annotate molluscan ---
  # top_hits <- full_hits %>%
  #   group_by(SNP) %>%
  #   slice_min(order_by = E_value, n = 1, with_ties = FALSE) %>%
  #   ungroup() %>%
  #   mutate(Molluscan = SNP %in% mollusc_hits$SNP)
  # 
  # # --- Match top hits with outliers ---
  # outliers_funct <- Out_qv %>%
  #   inner_join(top_hits, by = "SNP") %>%
  #   distinct(SNP, .keep_all = TRUE)
  # 
  # # --- Manhattan plot ---
  # manqv <- ggplot(rda_stats, aes(x = Position, y = -log10(p_values))) +
  #   geom_point(colour = "grey80") +
  #   geom_point(data = outliers_funct,
  #              aes(colour = Molluscan), size = 2.5, alpha = 0.7) +
  #   geom_text_repel(data = outliers_funct,
  #                   aes(label = str_wrap(Description, width = 30)),
  #                   size = 1.5, box.padding = 0.5,
  #                   segment.colour = "grey50", segment.size = 0.4) +
  #   scale_colour_manual(name = "Molluscan Hit",
  #                       values = c("TRUE" = "#E57373", "FALSE" = "#14FFF6")) +
  #   geom_hline(yintercept = -log10(thres_env), linetype = "dashed", colour = "grey60") +
  #   labs(x = "Position in Genome", y = "-log10(p-value)",
  #        title = paste0("Manhattan Plot with BLAST Hits: ", region)) +
  #   theme_classic()
  # 
  # plots[[region]] <- manqv
}

# --- Output ---
view(all_blast_hits)

# # Optional: show plots
# print(plots[["1000bp"]])


# Remove duplicates (one row per SNP, with molluscan hits in preference to full and the hit from the shortest flank (100, 500,1000))
filtered_hits <- all_blast_hits %>% 
  group_by(SNP) %>%
  # 1) If any molluscan hit exists for this SNP, keep only those; otherwise keep full
  filter(if (any(Database == "molluscan")) Database == "molluscan" else Database == "full") %>%
  # 2) Convert Proximity to numeric (100, 500, 1000)
  mutate(ProxNum = as.numeric(str_remove(Proximity, "bp"))) %>%
  # 3) Sort by proximity (smallest first) and then by best E-value
  arrange(ProxNum, E_value) %>%
  slice(1) %>%       # take the top row per SNP
  ungroup() %>%
  dplyr::select(-ProxNum) %>%   # drop the helper column
  mutate(
    Database = recode(
      Database,
      molluscan = "Molluscan",
      full = "Full"))
########

filtered_hits_with_raw_correlations <- Out_qv %>%
    inner_join(filtered_hits, by = "SNP") %>%
    distinct(SNP, .keep_all = TRUE) %>%
  dplyr::select(SNP, Query_ID, Subject_ID, Description, E_value, p_values, predictor, correlation, Database, Proximity) %>% 
  mutate(
    predictor = recode(
      predictor,
      catchment_area_km = "Catchment Area",
      WC_alt            = "Altitude",
      FW_hydro_wavg_05  = "BIO5",           # updated
      FW_hydro_wavg_03  = "BIO3",
      FW_dem_range      = "Elevation Range",
      FW_hydro_wavg_17  = "BIO17",
      wildness          = "Wildness Index"
    ),
    E_value   = formatC(E_value, format = "e", digits = 3),
    p_values   = signif(p_values, 4),
    correlation = signif(correlation, 3))
unique(filtered_hits_with_raw_correlations$predictor)
head(filtered_hits_with_raw_correlations)

write.csv(filtered_hits_with_raw_correlations, file = "output/rda_qvalue_outliers_blast_results_evalcutoff1e-5_all100_500_1000bp_flanks.csv", row.names = FALSE)


##Plot it up
# Re‑join to get manhatten positionsand also the ncbi database type
annotated_hits <- Out_qv %>%
  inner_join(filtered_hits, by = "SNP") %>%
  distinct(SNP, .keep_all = TRUE)

# Manhattan plot
manqv <- ggplot(rda_stats, aes(x = Position, y = -log10(p_values))) +
  geom_point(colour = "grey80") +
  geom_point(data = annotated_hits, aes(colour = Database),size  = 3, alpha = 0.5) +
  geom_text_repel(data = annotated_hits,aes(label = str_wrap(Description, width = 30)),
    size         = 2,box.padding  = 0.5,segment.colour = "grey50",segment.size = 0.5, max.overlaps = 10000) +
  scale_colour_manual(name   = "Database",values = c(Full = "#14FFF6", Molluscan = "#E57373")  # or colours of your choice
  ) +
  geom_hline(yintercept = -log10(thres_env),linetype   = "dashed",colour     = "grey80") +
  labs(
    x     = "Position in Genome",
    y     = "-log10(p-value)",
    title = paste0(
      "Outlier SNPs with BLAST Hits (n = ",
      nrow(annotated_hits),
      ")"
    )) +
  theme_classic()

manqv
ggsave(manqv, width = 15, height = 7, units = "in", dpi = 600,
       filename = "./output/figures/ncbi_hit_from_qval_snps_manhatten_evalcutoff1e-5_all100_500_1000bp_flank.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_adaptivelandscape_cleaned/output/figures/ncbi_hit_from_qval_snps_manhatten_evalcutoff1e-5_all100_500_1000bp_flank.png")


#======OLD CODE========
### 1. Parse BLAST XML
blast_xml <- read_xml("output/blast_results_1000bp.xml")
# blast_xml <- read_xml("output/blast_results/blast_results_1000bp_full.xml")

iterations <- xml_find_all(blast_xml, ".//Iteration")

extract_hits <- function(iter) {
  query_id <- xml_text(xml_find_first(iter, ".//Iteration_query-def"))
  hits <- xml_find_all(iter, ".//Hit")
  
  if (length(hits) == 0) return(NULL)
  
  do.call(rbind, lapply(hits, function(hit) {
    data.frame(
      Query_ID    = query_id,
      Subject_ID  = xml_text(xml_find_first(hit, ".//Hit_id")),
      Description = xml_text(xml_find_first(hit, ".//Hit_def")),
      E_value     = as.numeric(xml_text(xml_find_first(hit, ".//Hsp_evalue"))),
      stringsAsFactors = FALSE
    )
  }))
}

# blast_hits_all <- lapply(iterations, extract_hits) %>% bind_rows() #takes time, gives 5595 hits

### 2. Filter and extract top hits (5595 hits to 2640 hits)
blast_hits_filtered <- blast_hits_all %>%
  filter(E_value <= 1e-5) %>% #Evalue threshold of 1e-5
  mutate(SNP = sub("::.*", "", Query_ID)) %>%  #get the clean SNP name
  filter(!grepl("genome assembly", Description, ignore.case = TRUE)) %>%  #remove SNP hits with the words "genome assembly" in because its just chromosomes sequencing
  filter(!grepl("microsatellite", Description, ignore.case = TRUE))  #remove SNP hits with the words "microsattelite" it doesnt mean anything functionally

top_hits <- blast_hits_filtered %>% # extract top hits (top e value hit for each individual snp) (so we get SNPs with hits not just 1000 hits for one SNP)
  group_by(SNP) %>%
  slice_min(order_by = E_value, n = 1, with_ties = FALSE) %>% #selects the n=1 best e value (smallest)
  ungroup() # 61 SNPs  had hits

### 4. Tag outlier SNPs with BLAST hit info
outliers_funct <- outliers %>%
  inner_join(top_hits, by = "SNP") %>%
  distinct(SNP, .keep_all = TRUE)

### 5. Plot Manhattan plot

manqv <- ggplot(rda_stats, aes(x = Position, y = -log10(p_values))) +
  geom_point(colour = "grey80") +
  geom_point(data = outliers %>% 
               filter(SNP %in% blast_hits_clean$SNP),
             aes(colour = "BLAST Hits"), size = 3, alpha = 0.5) +
  geom_text_repel(data = outliers_funct,
                  aes(label = str_wrap(Description, width = 30)),
                  size = 1.5, box.padding = 0.5,
                  segment.colour = "grey50", segment.size = 0.5) +
  scale_colour_manual(name = "SNP Category",
                      values = c("BLAST Hits" = "#14FFF6")) +
  geom_hline(yintercept = -log10(thres_env), linetype = "dashed", colour = "grey80") +
  labs(x = "Position in Genome", y = "-log10(p-value)",
       title = paste0("Outliers with NCBI hits (n = ", nrow(outliers_funct), ")")) +
  theme_classic()
manqv
ggsave(manqv, width = 9, height = 3, units = "in", dpi = 600,
       filename = "./output/figures/ncbi_hit_from_qval_snps_manhatten.png")
top_hits %>% print(n=34)








####### Original code #############

# Install necessary package if not already installed
# if (!require("xml2")) install.packages("xml2")

# Load the XML package
library(xml2)

# Function to extract hits from each iteration
extract_hits <- function(iter) {
  query_id <- xml_text(xml_find_first(iter, ".//Iteration_query-def"))
  hits <- xml_find_all(iter, ".//Hit")
  
  if (length(hits) == 0) return(NULL)
  
  do.call(rbind, lapply(hits, function(hit) {
    subject_id <- xml_text(xml_find_first(hit, ".//Hit_id"))
    description <- xml_text(xml_find_first(hit, ".//Hit_def"))
    evalue <- xml_text(xml_find_first(hit, ".//Hsp_evalue"))  # only first HSP
    
    data.frame(
      Query_ID = query_id,
      Subject_ID = subject_id,
      Description = description,
      E_value = as.numeric(evalue),
      stringsAsFactors = FALSE
    )
  }))
}

#### 1000bp
# Load the two BLAST XML files
blast1 <- read_xml("output/blast_results_1000bp.xml")

# Extract <Iteration> nodes (one per query)
iter1 <- xml_find_all(blast1, ".//Iteration")


# Combine all iterations
all_iterations <- c(iter1)

# Apply extraction to all iterations
blast_hits_list <- lapply(all_iterations, extract_hits)

# Combine into one data frame
blast_hits_1000bp <- do.call(rbind, blast_hits_list)

# Preview results
head(blast_hits)
length(unique(blast_hits_1000bp$Query_ID)) #5595 hits

# Optional: write to CSV
# write.csv(blast_hits_1000, "output/blast_hits.csv", row.names = FALSE)

library(dplyr)

# Ensure E_value is numeric
blast_hits_1000bp$E_value <- as.numeric(blast_hits_1000bp$E_value)

# Get top hit (lowest e-value) for each query
top_hits_1000bp_nofilt <- blast_hits_1000bp %>%
  group_by(Query_ID) %>%
  slice_min(order_by = E_value, n = 1, with_ties = FALSE) %>%
  ungroup()

# Filter by E-value cutoff first (e.g., 1e-5)
blast_hits_1000bp_evalcutoff_nofilt<- blast_hits_1000bp %>%
  filter(E_value <= 1e-5)

# Get top hit (lowest e-value) for each query
top_hits_1000bp_evalcutoff_nofilt <- blast_hits_1000bp_evalcutoff_nofilt %>%
  group_by(Query_ID) %>%
  slice_min(order_by = E_value, n = 1, with_ties = FALSE) %>%
  ungroup()
nrow(top_hits_1000bp_evalcutoff_nofilt)

# View result
head(top_hits)
# write.csv(top_hits_1000bp, file = "output/rda_qvalue_outliers_blast_results_tophits_1000bpflanks.csv", row.names = FALSE)
# write.csv(blast_hits_1000bp, file = "output/rda_qvalue_outliers_blast_results_allhits_1000bpflanks.csv", row.names = FALSE)
# write.csv(blast_hits_1000bp_evalcutoff, file = "output/rda_qvalue_outliers_blast_results_evalcutoff1e-5_1000bpflanks.csv", row.names = FALSE)

## Investigating the SNPs
view(blast_hits_1000bp_evalcutoff_nofilt)

blast_hits_1000bp_evalcutoff_nofilt %>%
  # filter(grepl("splicing|tubulin|heat|Se-GPX|SETMAR", Description)) %>% #i dont think the SETMAR is a great match, but these all seem p cool functional genes
  distinct(Query_ID, .keep_all = TRUE, ignore.case = TRUE)%>%
  arrange(E_value)  # Sort by E-value smallest to largers

# Extract SNP names from Query_ID by splitting on '::' and keeping only the SNP portion
func_snp_names_nofilt <- blast_hits_1000bp_evalcutoff_nofilt %>%
  # filter(grepl("splicing|tubulin|heat|Se-GPX|SETMAR", Description)) %>% #i dont think the SETMAR is a great match, but these all seem p cool functional genes
  distinct(Query_ID) %>%
  mutate(SNP_name = sub("::.*", "", Query_ID)) %>%  # Extract SNP portion
  pull(SNP_name)
func_snp_names
# Subset the outliers dataset based on SNP names
outliers_funct_nofilt <- outliers %>%
  filter(SNP %in% func_snp_names_nofilt)
outliers_funct_nofilt
# outliers_funct_tagged_nofilt <- outliers_funct_nofilt %>%
#   left_join(blast_hits_1000bp_evalcutoff_nofilt %>% dplyr::select(SNP, Description), by = "SNP")
# outliers_funct
# head(blast_hits_1000bp_evalcutoff)
# 
# 
# # Extract SNP names from Query_ID by splitting on '::' and keeping only the SNP portion
# blasthit_snpname <- blast_hits_1000bp_evalcutoff %>%
#   # filter(grepl("splicing|tubulin|heat|Se-GPX|SETMAR", Description)) %>% #i dont think the SETMAR is a great match, but these all seem p cool functional genes
#   distinct(Query_ID) %>%
#   mutate(SNP_name = sub("::.*", "", Query_ID)) %>%  # Extract SNP portion
#   pull(SNP_name)
# blasthit_snpname



library(dplyr)

# 1. Create a “clean” BLAST hits table with SNP and Description
blast_hits_clean_nofilt <- blast_hits_1000bp_evalcutoff_nofilt %>%
  mutate(
    SNP = sub("::.*", "", Query_ID)                # drop the ::… suffix
  ) %>%
  dplyr::select(SNP, Description, E_value)

# 2. Pull out only those outliers that have a cool functional hit,
#    and carry the Description along in one join/filter step
outliers_funct_tagged_nofilt <- outliers %>%
  inner_join(blast_hits_clean_nofilt,
             by = "SNP") %>%
  distinct(SNP, .keep_all = TRUE)
outliers_funct_tagged_nofilt
# 3. Now plot — with geom_text_repel() labelling only your functional outliers
library(ggplot2)
library(ggrepel)

manqv <- ggplot(rda_stats, aes(x = Position, y = -log10(p_values))) +
  geom_point(colour = "grey80") +
  geom_point(data = outliers %>% filter(SNP %in% blast_hits_clean_nofilt$SNP[blast_hits_clean_nofilt$E_value <= 1e-5]),
             aes(colour = "Top BLAST Hits"), size = 3, alpha = 0.5) +
  geom_point(data = outliers_funct_tagged_nofilt, aes(colour = "Functional Outliers"), size = 3) +
  geom_text_repel(data = outliers_funct_tagged_nofilt,
                  aes(label = Description), size = 1.5, box.padding = 0.5, max.overlaps = Inf, 
                  segment.colour = "grey50",    # line colour
                  segment.size   = 0.5,         # line thickness
                  segment.alpha  = 0.8          # line transparency
  ) +
  scale_colour_manual(name   = "SNP Category",
                      values = c("Top BLAST Hits"     = "#14FFF6",
                                 "Functional Outliers"= "deeppink")) +
  geom_hline(yintercept = -log10(thres_env),
             linetype   = "dashed",
             colour     = "grey80") +
  labs(x     = "Position in Genome",
       y     = "-log10(p-value)",
       title = paste0("Outliers with interesting hits on NCBI (n = ",
                      nrow(outliers_funct_tagged_nofilt),
                      ")")) +
  theme_classic()

manqv
library(stringr)

manqv <- ggplot(rda_stats, aes(x=Position, y=-log10(p_values))) +
  geom_point(colour="grey80") +
  geom_point(data = outliers %>% filter(SNP %in% blast_hits_clean_nofilt$SNP[blast_hits_clean_nofilt$E_value <= 1e-5]),
             aes(colour = "BLAST Hits"), size = 3, alpha = 0.5) +
  # geom_point(data=outliers_funct_tagged_nofilt, aes(colour="Functional Outliers"), size=3) +
  geom_text_repel(data=outliers_funct_tagged_nofilt,
                  aes(label=str_wrap(Description, width=30)),
                  size=1.5, box.padding=0.5, max.overlaps=Inf,
                  segment.colour="grey50", segment.size=0.5, segment.alpha=0.8) +
  scale_colour_manual(name="SNP Category",
                      values=c("BLAST Hits"="#14FFF6","Functional Outliers"="deeppink")) +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", colour="grey80") +
  labs(x="Position in Genome", y="-log10(p-value)",
       title=paste0("Outliers with interesting hits on NCBI -whole database search (n=",nrow(outliers_funct_tagged_nofilt),")")) +
  theme_classic()

manqv
ggsave(plot= manqv,
       width = 9,
       height = 3,
       units = c( "in"),
       dpi = 600,
       filename = "./output/figures/ncbi_hit_from_qval_snps_manhatten_nomolluscanBLASTfilter.png")
shell.exec("C:/Users/r01vg21/OneDrive - University of Aberdeen/2.PROJECT/scotpearl_popgen/output/figures/ncbi_hit_from_qval_snps_manhatten_nomolluscanBLASTfilter.png")

unique(outliers_funct_tagged$Description)

#How many SNPs had hits?
nrow(outliers_funct_tagged_nofilt)
nrow(top_hits_1000bp_evalcutoff_nofilt)
nrow(top_hits_1000bp_evalcutoff)
nrow(outliers_funct_tagged)
top_hits_1000bp_evalcutoff
view(blast_hits_1000bp_evalcutoff)
view(outliers_funct_tagged_nofilt)
# write.csv(blast_hits_1000bp_evalcutoff, file = "output/rda_qvalue_outliers_blast_results_evalcutoff1e-5_1000bpflanks.csv", row.names = FALSE)
