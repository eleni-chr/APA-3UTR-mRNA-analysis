# Post-APAlyzer analysis - TDP 12m

# Before running this, run the R script "APAlyzer on TDP 12m files"

# Load required libraries
library(multiMiR)
library(dplyr)
library(readr)
library(Biostrings)
library(ggplot2)
library(purrr)
library(biomaRt)
library(tidyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(broom)

# Clear workspace
rm(list = ls())

# ---------------------------------------------------------------------------
# Get miRNA–gene interaction data via multiMiR for dysregulated genes (DESeq2)
# ---------------------------------------------------------------------------

# Read DESeq2 mRNA-seq results (with Ensembl IDs)
mrna_de <- read_csv("DESeq2_mRNA_TDP_12m.csv", col_types = cols())

# Filter for significantly dysregulated genes (using p-adjusted < 0.05)
dysregulated_mrna <- mrna_de %>% filter(`p-adjusted` < 0.05)

# Count upregulated and downregulated genes
num_upregulated <- dysregulated_mrna %>% filter(LFC > 0) %>% nrow()
num_downregulated <- dysregulated_mrna %>% filter(LFC < 0) %>% nrow()

# Remove version numbers from Ensembl IDs if needed
dysregulated_mrna <- dysregulated_mrna %>% mutate(Gene = gsub("\\..*$", "", Gene))

# Map Ensembl IDs to gene symbols using biomaRt
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 filters = "ensembl_gene_id",
                 values = dysregulated_mrna$Gene,
                 mart = mart)
mapping <- mapping %>% mutate(ensembl_gene_id = as.character(ensembl_gene_id))

# Merge mapping with DESeq2 results
dysregulated_mrna <- dysregulated_mrna %>% 
  left_join(mapping, by = c("Gene" = "ensembl_gene_id"))

names(dysregulated_mrna)[names(dysregulated_mrna) == "external_gene_name"] <- "gene_symbol"

dysregulated_mrna <- dysregulated_mrna %>% filter(!is.na(gene_symbol)) %>%
  mutate(gene_symbol = as.character(gene_symbol))

# Get list of dysregulated genes (unique gene symbols)
dysregulated_genes <- unique(dysregulated_mrna$gene_symbol)

# Query multiMiR for only validated interactions for these genes in mouse.
batch_size <- 50
gene_batches <- split(dysregulated_genes, ceiling(seq_along(dysregulated_genes) / batch_size))

query_batch <- function(genes) {
  cat("Querying", length(genes), "genes...\n")
  res <- get_multimir(target = genes, table = "validated", org = "mmu")
  return(res@data)  
}

results_list <- map(gene_batches, query_batch)
combined_results <- bind_rows(results_list) %>% distinct()

# Select gene–miRNA pairs and remove duplicates.
selected_interactions <- combined_results %>%
  dplyr::select(target_symbol, mature_mirna_id) %>% 
  distinct()
colnames(selected_interactions) <- c("gene_symbol", "miRNA")

# Retain only mouse miRNAs (names beginning with "mmu") because MultiMiR can also give hsa miRNAs targeting mmu genes.
selected_interactions <- selected_interactions %>%
  filter(grepl("^mmu", miRNA)) %>%
  filter(!grepl("p\\.", miRNA)) #also discard those like "mmu-miR-21a-5p.1" which are still duplicates.

write.csv(selected_interactions, file = "miRNAs_targeting_dysregulated_genes_TDP_12m.csv", row.names = FALSE)

# ----------------------------------------------------------------
# Get list of dysregulated genes that are also APA-regulated genes
# ----------------------------------------------------------------

# Read APAlyzer results.
apa_results <- read_csv("test_3UTR_TDP_12m.csv", col_types = cols())

# Filter for APA regulated genes (only those with APAreg "DN" or "UP")
apa_regulated_genes <- apa_results %>%
  filter(APAreg %in% c("DN", "UP"))
apa_gene_list <- unique(apa_regulated_genes$gene_symbol)

# Total unique APA-regulated genes
total_APA <- length(apa_gene_list)

# Count unique genes with UP and DN APA regulation
num_UP <- apa_results %>% 
  filter(APAreg == "UP") %>% 
  pull(gene_symbol) %>% 
  unique() %>% 
  length()

num_DN <- apa_results %>% 
  filter(APAreg == "DN") %>% 
  pull(gene_symbol) %>% 
  unique() %>% 
  length()

# Calculate percentages
perc_UP <- (num_UP / total_APA) * 100
perc_DN <- (num_DN / total_APA) * 100

# Get the overlapping genes between APA-regulated and dysregulated gene lists
overlap_genes <- intersect(apa_gene_list, dysregulated_genes)

# Count how many genes are in the overlap
num_overlap <- length(overlap_genes)

# -------------------------------------------------------------------------------------------
# Extract miRNA seeds for miRNAs targeting genes that are BOTH dysregulated and APA-regulated
# -------------------------------------------------------------------------------------------

# The seeds are extracted from miRNAs that were predicted by MultiMiR that are actually available on miRBase AND also detected in our data.

# Get unique miRNAs from the gene-miRNA interactions.
predicted_miRNAs <- unique(selected_interactions$miRNA) # These are predicted from MultiMiR

# Read the miRBase FASTA file containing mature miRNA sequences.
mouse_miRNAs <- readDNAStringSet("mature.fa")  # downloaded from miRBase.

# Filter by names starting with "mmu"
mmu_idx <- startsWith(names(mouse_miRNAs), "mmu")
mouse_mmu_miRNAs <- mouse_miRNAs[mmu_idx]

# Trim names: remove anything after the first whitespace (if needed)
names(mouse_miRNAs) <- sapply(strsplit(names(mouse_miRNAs), " "), `[`, 1)

# Filter the predicted miRNAs to include only those that exist on miRBase.
target_miRNAs <- mouse_miRNAs[names(mouse_miRNAs) %in% predicted_miRNAs] # These are those from MultiMiR after filtering to keep only those from miRBase.

# Further filter to include only those miRNAs that are detected.
# Read DESeq2 miRNA results. 
miRNA_de <- read_csv("DESeq2_miRNA_TDP_12m.csv", col_types = cols())

# Extract the list of miRNAs detected in our data (irrespective of DESeq2 significance).
# This tells us which miRNAs are actually expressed in our samples.
detected_miRNAs <- unique(miRNA_de$miRNA)

# Further filter to include only those miRNAs that are detected.
target_miRNAs_final <- target_miRNAs[names(target_miRNAs) %in% detected_miRNAs] # These are those from MultiMiR after filtering to keep only those from miRBase and then filtering again to include only those detected in our samples

# Extract three seeds per miRNA (6-mer, 7-mer, and 8-mer)
extract_seeds <- function(seq) {
  seq_str <- as.character(seq)
  if(nchar(seq_str) < 8) {
    warning("Sequence is too short to extract all seeds. Returning NA.")
    return(list(seed6 = NA, seed7 = NA, seed8 = NA))
  } else {
    seed6 <- substr(seq_str, 2, 7)   # 6-mer: positions 2-7
    seed7 <- substr(seq_str, 2, 8)   # 7-mer: positions 2-8
    seed8 <- substr(seq_str, 1, 8)   # 8-mer: positions 1-8
    return(list(seed6 = seed6, seed7 = seed7, seed8 = seed8))
  }
}
miR_seeds_list <- lapply(target_miRNAs_final, extract_seeds)

# Remove entries with NA seeds, if necessary.
miR_seeds_list <- miR_seeds_list[!sapply(miR_seeds_list, function(x) any(is.na(unlist(x))))]

# Create a data frame for reference.
seed_df <- data.frame(
  miRNA = names(miR_seeds_list),
  seed6 = sapply(miR_seeds_list, function(x) x$seed6),
  seed7 = sapply(miR_seeds_list, function(x) x$seed7),
  seed8 = sapply(miR_seeds_list, function(x) x$seed8),
  stringsAsFactors = FALSE
)

write.csv(seed_df, file = "seeds_predicted_to_target_dysregulated_genes_TDP_12m.csv", row.names = FALSE)

# Create a named list of seeds per miRNA for downstream analysis.
miR_seeds <- lapply(1:nrow(seed_df), function(i) {
  list(seed6 = seed_df$seed6[i],
       seed7 = seed_df$seed7[i],
       seed8 = seed_df$seed8[i])
})
names(miR_seeds) <- seed_df$miRNA

# Check how many of the predicted miRNAs are actually detected in our samples
predicted_detected <- intersect(predicted_miRNAs, detected_miRNAs)
num_predicted_detected <- length(predicted_detected)

# ---------------------------------------------------------------------------------------------------
# miRNA binding site analysis on 3'UTR isoforms of genes that are BOTH dysregulated and APA-regulated
# ---------------------------------------------------------------------------------------------------

# Extract 3′UTR sequences for isoforms (only for dysregulated genes)

# The sequences themselves are extracted from the reference genome using the 
# genomic coordinates stored in the APAlyzer workspace. In other words, the UTR 
# segments (both constitutive and alternative) are obtained via the 
# BSgenome.Mmusculus.UCSC.mm10 package using coordinates (from UTRdbraw) that 
# were determined during the APA analysis. So although the APA-regulated genes 
# were identified based on the RNA‑seq data, the actual short and long isoform
# sequences come from a reference database (the annotated mm10 genome), not 
# directly from the RNA‑seq reads.So, they represent the complete sequences for
# the 3′UTR isoforms

# Load the APAlyzer workspace
load("APAlyzer_TDP_12m_workspace.RData")

# Split UTRdbraw into cUTR and aUTR parts using the element names
cUTR_list <- UTRdbraw[grep("__cUTR$", names(UTRdbraw))]
aUTR_list <- UTRdbraw[grep("__aUTR$", names(UTRdbraw))]

# Extract gene symbols from the GRanges metadata for cUTR and aUTR
cUTR_genes <- sapply(cUTR_list, function(gr) as.character(mcols(gr)$gene_symbol))
aUTR_genes <- sapply(aUTR_list, function(gr) as.character(mcols(gr)$gene_symbol))

# All genes represented in UTRdbraw (from the cUTR parts)
all_genes <- unique(cUTR_genes)

# Filter the gene list to include only genes that are BOTH dysregulated and APA-regulated
regulated_gene_list <- intersect(all_genes, overlap_genes)

# Create empty lists to store sequences for each gene (only genes that are BOTH dysregulated and APA-regulated)
short_isoform <- list()  # cUTR only
long_isoform  <- list()  # cUTR concatenated with aUTR (if available)

# Load the reference genome for mm10
genome <- BSgenome.Mmusculus.UCSC.mm10

# Loop over each gene (that is BOTH dysregulated and APA-regulated)
for (gene in regulated_gene_list) {
  # Get the cUTR GRanges for this gene
  idx_c <- which(cUTR_genes == gene)
  if (length(idx_c) == 0) next  # Skip if no cUTR found (should not happen)
  
  # Combine all cUTR segments for the gene
  gr_c <- unlist(cUTR_list[idx_c])
  gr_c <- sort(gr_c)  # Order by genomic coordinate
  
  # Extract the sequence for the cUTR segments
  seq_c <- getSeq(genome, gr_c)
  
  # Concatenate segments to form the complete short isoform
  seq_c_concat <- do.call(xscat, as.list(seq_c))
  short_isoform[[gene]] <- seq_c_concat
  
  # Check if any aUTR segments exist for the gene
  idx_a <- which(aUTR_genes == gene)
  if (length(idx_a) > 0) {
    gr_a <- unlist(aUTR_list[idx_a])
    gr_a <- sort(gr_a)
    seq_a <- getSeq(genome, gr_a)
    seq_a_concat <- do.call(xscat, as.list(seq_a))
    
    # The long isoform is the cUTR sequence concatenated with the aUTR sequence.
    long_isoform[[gene]] <- xscat(seq_c_concat, seq_a_concat)
  } else {
    
    # If no aUTR segments are found, the long isoform equals the short isoform.
    long_isoform[[gene]] <- seq_c_concat
  }
}

# Convert the lists to DNAStringSet objects for export
short_isoform_set <- DNAStringSet(short_isoform)
long_isoform_set  <- DNAStringSet(long_isoform)

# Write the sequences to FASTA files
writeXStringSet(short_isoform_set, filepath = "TDP_12m_ShortIsoforms_DYS_plus_APA_reg.fasta")
writeXStringSet(long_isoform_set, filepath = "TDP_12m_LongIsoforms_DYS_plus_APA_reg.fasta")

# Define a helper function to remove overlapping matches (e.g., a 7-mer may be a sub-sequence of an 8-mer).
removeOverlappingMatches <- function(matchList) {
  if (length(matchList) == 0) return(matchList)
  df <- do.call(rbind, lapply(matchList, function(x) {
    data.frame(seedType = x$seedType,
               start = x$start,
               width = x$width,
               stringsAsFactors = FALSE)
  }))
  df$end <- df$start + df$width - 1
  df <- df[order(-df$width, df$start), ]
  keep <- rep(TRUE, nrow(df))
  for (i in seq_len(nrow(df) - 1)) {
    if (keep[i]) {
      for (j in (i + 1):nrow(df)) {
        if (keep[j] && df$start[j] >= df$start[i] && df$end[j] <= df$end[i]) {
          keep[j] <- FALSE
        }
      }
    }
  }
  return(df[keep, ])
}

# Perform binding site search for each predicted miRNA using all three seed types.
results <- list()
for (miRNA in names(miR_seeds)) {
  seeds_for_miRNA <- miR_seeds[[miRNA]]
  # Ensure seeds are uppercase and compute their reverse complements.
  seeds_rc <- lapply(seeds_for_miRNA, function(s) {
    reverseComplement(DNAString(toupper(s)))
  })
  
  # Function to search for matches in a given isoform set (short or long)
  search_matches <- function(isoform_set) {
    lapply(names(isoform_set), function(gene) {
      gene_seq <- isoform_set[[gene]]
      match_list <- list()
      for (seed_type in names(seeds_rc)) {
        pattern <- seeds_rc[[seed_type]]
        matches <- matchPattern(pattern, gene_seq, fixed = TRUE)
        if (length(matches) > 0) {
          for (i in seq_along(matches)) {
            match_list[[length(match_list) + 1]] <- list(
              seedType = seed_type,
              start = start(matches)[i],
              width = width(matches)[i]
            )
          }
        }
      }
      if (length(match_list) > 0) {
        unique_matches <- removeOverlappingMatches(match_list)
      } else {
        unique_matches <- data.frame()  # No matches found.
      }
      return(unique_matches)
    })
  }
  
  matches_short <- search_matches(short_isoform_set)
  names(matches_short) <- names(short_isoform_set)
  
  matches_long <- search_matches(long_isoform_set)
  names(matches_long) <- names(long_isoform_set)
  
  results[[miRNA]] <- list(short = matches_short, long = matches_long)
}

# Helper function to count matches per gene.
count_matches <- function(match_list) {
  sapply(match_list, function(x) {
    if (is.data.frame(x) && nrow(x) > 0) nrow(x) else 0
  })
}

counts_short <- lapply(results, function(x) count_matches(x$short))
counts_long  <- lapply(results, function(x) count_matches(x$long))

# Combine counts into a data frame.
gene_ids <- names(short_isoform_set)
df_counts <- data.frame(Gene = gene_ids, stringsAsFactors = FALSE)
for (miRNA in names(miR_seeds)) {
  df_counts[[paste0(miRNA, "_short")]] <- counts_short[[miRNA]][gene_ids]
  df_counts[[paste0(miRNA, "_long")]]  <- counts_long[[miRNA]][gene_ids]
  df_counts[[paste0(miRNA, "_diff")]]  <- df_counts[[paste0(miRNA, "_long")]] - df_counts[[paste0(miRNA, "_short")]]
}

write.csv(df_counts, file = "miRNA_seed_counts_comparison_TDP_12m.csv", row.names = FALSE)

# ---------------------------------------------------------
# Integrate miRNA count data using DESeq2 results on miRNAs
# ---------------------------------------------------------

# Get the names of the miRNAs predicted to target the genes that are BOTH dysregulated and APA-regulated
# (after they have been filtered to include only those that exist on miRBase and detected in our data).
target_miRNA_names <- names(target_miRNAs_final)

# Count the numbers
num_target_detected <- length(target_miRNA_names)

# Filter miRNAs that are dysregulated (p-value < 0.05)
dysregulated_miRNAs <- miRNA_de %>%
  filter(`p-value` < 0.05) %>%
  pull(miRNA)

# Count common miRNAs
common_dys_target_miRNAs <- length(intersect(target_miRNA_names, dysregulated_miRNAs))

# ---------------------------------------------------------------------------------------
# Count how many gene-miRNA pairs have gained or lost binding sites due to APA regulation
# ---------------------------------------------------------------------------------------

# This analysis involves only the gene-miRNA pairs containing miRNAs that exist 
# in miRBase, which are predicted (according to MultiMiR) to target genes that
# are BOTH dysregulated and APA-regulated.

# Reshape df_counts to long format so that each row is a gene-miRNA pair.
# Assumes that columns are named like "miRNA_short", "miRNA_long", and "miRNA_diff".
# We pivot the data so that we can easily compare short, long, and diff counts.
df_long <- df_counts %>%
  pivot_longer(
    cols = -Gene,
    names_to = c("miRNA", "measure"),
    names_sep = "_",
    values_drop_na = TRUE
  ) %>%
  pivot_wider(
    names_from = measure,
    values_from = value
  )

# Filter out gene-miRNA pairs where both isoforms have zero binding sites.
# Those filtered out are impossible gene-miRNA interactions according to our mRNA-seq data based on our gene sequences.
df_long_filtered <- df_long %>% filter(!(short == 0 & long == 0))

# Further filter gene–miRNA interactions to only include those where the miRNA is actually detected in our samples, according to our miRNA-seq data.
df_long_filtered_detected <- df_long_filtered %>% 
  filter(miRNA %in% detected_miRNAs)

write.csv(df_long_filtered_detected, file = "filtered_gene_miRNA_interactions_detected_TDP_12m.csv", row.names = FALSE)

# -----------------------------------------------------------------------------------------------------------------------
# Get a summary of total gene-miRNA pairs that have gained, lost, or had no change in binding sites due to APA regulation
# -----------------------------------------------------------------------------------------------------------------------

# Count gene-miRNA pairs by net change in binding site counts.
summary_counts <- df_long_filtered_detected %>%
  summarise(
    total_pairs = n(),
    gain = sum(diff > 0, na.rm = TRUE),
    loss = sum(diff < 0, na.rm = TRUE),
    no_change = sum(diff == 0, na.rm = TRUE)
  )

# Also get counts per miRNA
per_miRNA_summary <- df_long_filtered_detected %>%
  group_by(miRNA) %>%
  summarise(
    pairs = n(),
    gain = sum(diff > 0, na.rm = TRUE),
    loss = sum(diff < 0, na.rm = TRUE),
    no_change = sum(diff == 0, na.rm = TRUE)
  )

write.csv(per_miRNA_summary, file = "gene_miRNA_binding_site_summary_TDP_12m.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# Calculate the proportion of genes that have gained/lost binding sites
# ---------------------------------------------------------------------

# Aggregate gene–miRNA pairs by gene using the sum of binding site differences.
gene_summary <- df_long_filtered_detected %>%
  group_by(Gene) %>%
  summarise(
    sum_diff = sum(diff, na.rm = TRUE),
    n_interactions = n()
  ) %>%
  # Classify each gene based on its total binding site difference.
  mutate(change_category = case_when(
    sum_diff > 0 ~ "gain",
    sum_diff < 0 ~ "loss",
    TRUE ~ "no_change"
  ))

# Count the number of unique genes in each category and calculate proportions.
gene_summary_counts <- gene_summary %>%
  summarise(
    total_genes = n(),
    gain = sum(change_category == "gain"),
    loss = sum(change_category == "loss"),
    no_change = sum(change_category == "no_change")
  ) %>%
  mutate(
    gain_prop = gain / total_genes,
    loss_prop = loss / total_genes,
    no_change_prop = no_change / total_genes
  )

# Optionally, write the gene target proportions to a CSV file.
write.csv(gene_summary, file = "gene_target_binding_site_sums_TDP_12m.csv", row.names = FALSE)

# ----------------
# Plot the results
# ----------------

# 1. Bar plot of Dysregulated Genes (Upregulated vs Downregulated)
dysregulated_df <- data.frame(
  Regulation = c("Upregulated", "Downregulated"),
  Count = c(num_upregulated, num_downregulated)
)

p1 <- ggplot(dysregulated_df, aes(x = Regulation, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste(Count, "genes")), vjust = -0.5) +
  labs(title = "Dysregulated Genes in Tg vs. NTg", 
       subtitle = paste("Total =", length(dysregulated_genes)),
       y = "Number of Genes", 
       x = NULL) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none")

print(p1)
ggsave("p1_TDP_12m.svg", plot = p1, width = 8, height = 6, device = "svg")

# 2. Bar plot of APA-regulated Genes (UP vs DN)
APA_df <- data.frame(
  APA_Regulation = c("3'UTR Lengthening", "3'UTR Shortening"),
  Count = c(num_UP, num_DN)
)

p2 <- ggplot(APA_df, aes(x = APA_Regulation, y = Count, fill = APA_Regulation)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste(Count, "genes")), vjust = -0.5) +
  labs(title = "APA-regulated Genes in Tg vs. NTg",
       subtitle = paste("Total =", total_APA),
       y = "Number of Genes", 
       x = "APA Regulation") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none")

print(p2)
ggsave("p2_TDP_12m.svg", plot = p2, width = 8, height = 6, device = "svg")

# 3. Bar plot for miRNA detection from MultiMiR after filtering steps
miRNA_detection_df <- data.frame(
  Category = c("Total Predicted\nby MultiMiR", "Detected and\non miRBase", "Predicted and\non miRBase", "Predicted and\ndetected in\nour Samples"),
  Count = c(length(predicted_miRNAs), length(target_miRNAs_final), length(target_miRNAs), num_predicted_detected)
)

p3 <- ggplot(miRNA_detection_df, aes(x = reorder(Category, -Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste(Count, "miRNAs")), vjust = -0.5) +
  labs(title = "miRNAs Predicted by MultiMiR to Target Genes\nthat are both Dysregulated and APA-regulated",
       subtitle = sprintf("Detected: %.2f%% of predicted miRNAs", (length(target_miRNAs) / length(predicted_miRNAs)) * 100),
       y = "Number of miRNAs", 
       x = NULL) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none")

print(p3)
ggsave("p3_TDP_12m.svg", plot = p3, width = 8, height = 6, device = "svg")

# -----------------
# Prepare some data
# -----------------

# Merge mRNA DESeq2 data with APA results by gene symbol.
merged_gene_data <- dplyr::inner_join(
  dysregulated_mrna %>% dplyr::select(gene_symbol, LFC),
  apa_results %>% dplyr::select(gene_symbol, RED, APAreg),
  by = "gene_symbol"
)

# For the subset of genes that are part of the "actually possible" gene-miRNA pairs,
# extract the unique gene symbols from df_long_filtered_detected.
possible_genes <- unique(df_long_filtered_detected$Gene)

# Filter the merged data to only include these genes.
merged_gene_data_possible <- merged_gene_data %>% filter(gene_symbol %in% possible_genes)

# Create a composite grouping variable combining Direction and APA regulation type.
merged_gene_data_possible <- merged_gene_data_possible %>% 
  mutate(
    Direction = ifelse(LFC > 0, "Upregulated", "Downregulated"),
    APA_type = ifelse(APAreg == "UP", "3'UTR lengthening", "3'UTR shortening"),
    Group = paste(Direction, APA_type, sep = " + ")
  )

# Merge gene-level binding site gain data (gene_summary) with APA/DESeq2 merged data (merged_gene_data_possible)
# Assuming gene_summary has columns: Gene, sum_diff, change_category, and merged_gene_data_possible has gene_symbol, LFC, RED, APAreg, etc.
merged_gene_data_possible2 <- inner_join(
  merged_gene_data_possible,
  gene_summary,
  by = c("gene_symbol" = "Gene")
)

# Recode APA regulation for clarity
merged_gene_data_possible2 <- merged_gene_data_possible2 %>%
  mutate(APA_type = ifelse(APAreg == "UP", "3'UTR lengthening", "3'UTR shortening"))

# -----------------
# Plot more results
# -----------------

# 5. Bar plot for gene-miRNA pairs by binding site change due to 3'UTR lengthening

# Filter gene summary for genes with 3'UTR lengthening
lengthening_gene_summary <- gene_summary[gene_summary$Gene %in% 
                                           merged_gene_data_possible2$gene_symbol[merged_gene_data_possible2$APA_type == "3'UTR lengthening"], ]

# Calculate summary counts for binding site changes for 3'UTR lengthening.
summary_counts_long <- list(
  gain = sum(lengthening_gene_summary$sum_diff > 0),
  no_change = sum(lengthening_gene_summary$sum_diff == 0),
  loss = sum(lengthening_gene_summary$sum_diff < 0)
)

# Create the data frame for the bar plot for lengthening genes
binding_site_df_long <- data.frame(
  Change = c("Gain", "No Change", "Loss"),
  Count = c(summary_counts_long$gain,
            summary_counts_long$no_change,
            summary_counts_long$loss)
)

p4 <- ggplot(binding_site_df_long, aes(x = Change, y = Count, fill = Change)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +
  labs(title = "Gene-miRNA Binding Site Changes\nDue to 3'UTR Lengthening",
       y = "Number of Gene-miRNA Pairs Affected", 
       x = "Binding Site Change") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none")

print(p4)
ggsave("p4_TDP_12m.svg", plot = p4, width = 8, height = 6, device = "svg")


# 6. Bar plot for gene-miRNA pairs by binding site change due to 3'UTR shortening

# Filter the gene summary (or merged data) to include only genes with 3'UTR shortening.
short_gene_summary <- gene_summary[gene_summary$Gene %in% 
                                     merged_gene_data_possible2$gene_symbol[merged_gene_data_possible2$APA_type == "3'UTR shortening"], ]

# Because the 'sum_diff' column in gene_summary is computed as (long - short),
# for SHORTENING genes the binding site loss is the absolute value of that difference.
short_gene_summary <- short_gene_summary %>% 
  mutate(binding_site_loss = abs(sum_diff))

# Recalculate the summary counts for shortening.
summary_counts_short <- list(
  loss = sum(short_gene_summary$sum_diff > 0), # For shortening, a positive value indicates loss, because the difference was calculated as long minus short.
  no_change = sum(short_gene_summary$sum_diff == 0),
  gain = sum(short_gene_summary$sum_diff < 0)
)

# Create the data frame for the bar plot.
binding_site_df_short <- data.frame(
  Change = c("Loss", "No Change", "Gain"),
  Count = c(summary_counts_short$loss,
            summary_counts_short$no_change,
            summary_counts_short$gain)
)

# Plot the bar plot for shortening genes.
p5 <- ggplot(binding_site_df_short, aes(x = Change, y = Count, fill = Change)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +
  labs(title = "Gene-miRNA Binding Site Changes\nDue to 3'UTR Shortening",
       y = "Number of Gene-miRNA Pairs Affected", 
       x = "Binding Site Change") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none")

print(p5)
ggsave("p5_TDP_12m.svg", plot = p5, width = 8, height = 6, device = "svg")

# Create a scatter plot with four trendlines.
p6 <- ggplot(merged_gene_data_possible, aes(x = RED, y = LFC, color = Group, shape = Group)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_shape_manual(values = c("Upregulated + 3'UTR lengthening" = 16,   # filled circle
                                "Upregulated + 3'UTR shortening" = 17,     # filled triangle
                                "Downregulated + 3'UTR lengthening" = 18,    # filled diamond
                                "Downregulated + 3'UTR shortening" = 15)) +  # filled square
  labs(title = "APA regulation vs. LFC for Genes Involved\nin Actually Possible Gene-miRNA Interactions",
       x = "Relative Expression Difference (APA Metric)\n[negative = 3'UTR shortening; positive = 3'UTR lengthening]",
       y = "mRNA log2 Fold Change\n(<0: downregulation     >0: upregulation)",
       color = NULL,   # Remove legend title for color
       shape = NULL) + # Remove legend title for shape
  theme_minimal(base_size = 18) +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE))

print(p6)
ggsave("p6_TDP_12m.svg", plot = p6, width = 8, height = 6, device = "svg")

# ------------------------------------------------------------------------------
# Compute Spearman correlation for each composite group (related to plot p6)
cor_results <- merged_gene_data_possible %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    spearman_r = if (n > 1) cor(RED, LFC, method = "spearman") else NA_real_,
    p_value   = if (n > 1) cor.test(RED, LFC, method = "spearman")$p.value else NA_real_
  )

print(cor_results)
# ------------------------------------------------------------------------------

# Focus on genes with 3'UTR lengthening
lengthening_genes <- merged_gene_data_possible2 %>% 
  filter(APA_type == "3'UTR lengthening")

# Create a scatter plot and compute Spearman correlation between net binding site gain (sum_diff) and mRNA LFC.
p7 <- ggplot(lengthening_genes, aes(x = sum_diff, y = LFC)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Net Binding Site Gain vs mRNA LFC",
       subtitle = "For Genes with 3'UTR Lengthening",
       x = "Number of Binding Sites Gained\nAcross All miRNAs Targeting a Gene",
       y = "mRNA log2 Fold Change\n(<0: downregulation     >0: upregulation)") +
  theme_minimal(base_size = 18)

print(p7)
ggsave("p7_TDP_12m.svg", plot = p7, width = 8, height = 6, device = "svg")

# ------------------------------------------------------------------------------
# Compute Spearman correlation between binding site loss and mRNA LFC (related to plot p7)
cor_lengthening <- cor.test(lengthening_genes$sum_diff, lengthening_genes$LFC, method = "spearman")
print(cor_lengthening)
# ------------------------------------------------------------------------------

# Focus on genes with 3'UTR shortening and compute binding site loss
shortening_genes <- merged_gene_data_possible2 %>% 
  filter(APA_type == "3'UTR shortening") %>%
  mutate(lost_diff = abs(-sum_diff))

# Create a scatter plot for binding site loss vs mRNA LFC
p8 <- ggplot(shortening_genes, aes(x = lost_diff, y = LFC)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Net Binding Site Loss vs mRNA LFC",
       subtitle = "For Genes with 3'UTR Shortening",
       x = "Number of Binding Sites Lost\nAcross All miRNAs Targeting a Gene",
       y = "mRNA log2 Fold Change\n(<0: downregulation     >0: upregulation)") +
  theme_minimal(base_size = 18)

print(p8)
ggsave("p8_TDP_12m.svg", plot = p8, width = 8, height = 6, device = "svg")

# ------------------------------------------------------------------------------
# Compute Spearman correlation between binding site loss and mRNA LFC (related to plot p8)
cor_shortening <- cor.test(shortening_genes$lost_diff, shortening_genes$LFC, method = "spearman")
print(cor_shortening)
# ------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Plot the results for a subset of the data (e.g., most dysregulated genes/miRNAs)
# --------------------------------------------------------------------------------

# Merge miRNA expression data with gene-miRNA interactions
# miRNA_de is the DESeq2 miRNA result data frame and has columns: "miRNA" and "Counts" (aka baseMean)
# df_long_filtered_detected contains gene-miRNA interactions with columns: Gene, miRNA, etc.
df_interactions_expr <- df_long_filtered_detected %>%
  left_join(miRNA_de %>% dplyr::select(miRNA, Counts), by = "miRNA")

# Summarise additional metrics per gene
# Compute the average miRNA expression (baseMean) and count the number of interactions per gene.
gene_miRNA_expr <- df_interactions_expr %>%
  group_by(Gene) %>%
  summarise(
    avg_miRNA_expr = mean(Counts, na.rm = TRUE),
    n_interactions = n()  # Number of gene-miRNA interactions (proxy for number of binding sites)
  )

# Merge with gene-level APA/DESeq2 data
# Assume merged_gene_data_possible2 (from previous steps) has columns: gene_symbol, LFC, RED, APAreg, sum_diff, etc.
# Merge gene_miRNA_expr (which uses "Gene") with merged_gene_data_possible2 (using gene_symbol)
merged_gene_data_final <- merged_gene_data_possible2 %>%
  left_join(gene_miRNA_expr, by = c("gene_symbol" = "Gene"))

# Re-label APA regulation for clarity
merged_gene_data_final <- merged_gene_data_final %>% 
  mutate(
    APA_type = ifelse(APAreg == "UP", "3'UTR lengthening", "3'UTR shortening")
  )

# Select the top 10 most upregulated genes with 3'UTR shortening (since we hypothesise to see a strong correlation with this combination).
top10_upregulated <- merged_gene_data_final %>%
  filter(APA_type == "3'UTR shortening") %>%
  slice_max(order_by = LFC, n = 10) %>%
  mutate(binding_site_loss = abs(sum_diff))  # Convert to positive values

# Create a scatter plot of binding site loss vs. mRNA LFC for these top 10 genes.
p9 <- ggplot(top10_upregulated, aes(x = binding_site_loss, y = LFC)) +
  geom_point(alpha = 0.7, size = 3, color = "darkgreen") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Binding Site Loss vs mRNA LFC",
       subtitle = "Top 10 Most Upregulated Genes (3'UTR Shortening)",
       x = "Binding Sites Lost",
       y = "mRNA log2 Fold Change\n(<0: downregulation     >0: upregulation)") +
  theme_minimal(base_size = 18)

print(p9)
ggsave("p9_TDP_12m.svg", plot = p9, width = 8, height = 6, device = "svg")

# ------------------------------------------------------------------------------
# Compute Spearman correlation for these top 10 genes (related to plot p9)
cor_top10 <- cor.test(top10_upregulated$binding_site_loss, top10_upregulated$LFC, method = "spearman")
print(cor_top10)
# ------------------------------------------------------------------------------

# Select the top 10 most downregulated genes with 3'UTR lengthening (since we hypothesise to see a strong correlation with this combination).
top10_downregulated <- merged_gene_data_final %>%
  filter(APA_type == "3'UTR lengthening") %>%  # downregulated events are expected here
  slice_min(order_by = LFC, n = 10)

# Create a scatter plot of net binding site gain (sum_diff) vs. mRNA LFC for these genes.
p10 <- ggplot(top10_downregulated, aes(x = sum_diff, y = LFC)) +
  geom_point(alpha = 0.7, size = 3, color = "darkred") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Binding Site Gain vs mRNA LFC",
       subtitle = "Top 10 Most Downregulated Genes (3'UTR lengthening)",
       x = "Binding Sites Gained",
       y = "mRNA log2 Fold Change\n(<0: downregulation     >0: upregulation)") +
  theme_minimal(base_size = 18)

print(p10)
ggsave("p10_TDP_12m.svg", plot = p10, width = 8, height = 6, device = "svg")

# ------------------------------------------------------------------------------
# Compute Spearman correlation for these top 10 genes (related to plot p10)
cor_top10_down <- cor.test(top10_downregulated$sum_diff, top10_downregulated$LFC, method = "spearman")
print(cor_top10_down)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Further subset the data to focus only on dysregulated miRNAs.
# ------------------------------------------------------------------------------

# Select the top 10 most upregulated genes with 3'UTR shortening, and further filter for gene-miRNA interactions involving downregulated miRNAs.

# 1. Filter miRNAs for those that are significantly downregulated (p-value < 0.05 and LFC < 0).
dysregulated_miRNAs_down <- miRNA_de %>%
  filter(`p-value` < 0.05, LFC < 0) %>%
  pull(miRNA)

# 2. Filter gene–miRNA interactions to include only those with downregulated miRNAs.
dys_interactions <- df_interactions_expr %>% 
  filter(miRNA %in% dysregulated_miRNAs_down)

# 3. Summarise additional metrics per gene using these dysregulated interactions.
gene_miRNA_expr_dys <- dys_interactions %>%
  group_by(Gene) %>%
  summarise(
    avg_miRNA_expr_dys = mean(Counts, na.rm = TRUE),
    n_interactions_dys = n()
  )

# 4. Merge with the gene-level APA/DESeq2 data.
merged_gene_data_dys <- merged_gene_data_final %>%
  filter(gene_symbol %in% gene_miRNA_expr_dys$Gene)

# 5. Select the top 10 most upregulated genes.
top10_upregulated2 <- merged_gene_data_dys %>%
  filter(APA_type == "3'UTR shortening") %>%
  slice_max(order_by = LFC, n = 10) %>%
  mutate(binding_site_loss = abs(sum_diff))  # Convert to positive values

# 6. Create the scatter plot (using sum_diff as binding site loss).
p11 <- ggplot(top10_upregulated2, aes(x = sum_diff, y = LFC)) +
  geom_point(alpha = 0.7, size = 3, color = "darkgreen") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Binding Site Loss vs. mRNA LFC",
       subtitle = "Top 10 Most Upregulated Genes (targeted by downregulated miRNAs)",
       x = "Binding Sites Lost (for downregulated miRNAs only)",
       y = "mRNA log2 Fold Change\n(<0: downregulation     >0: upregulation)") +
  theme_minimal(base_size = 18)

print(p11)

# ------------------------------------------------------------------------------
# Compute the Spearman correlation (related to plot p11)
cor_top10_up2 <- cor.test(top10_upregulated2$sum_diff, top10_upregulated2$LFC, method = "spearman")
print(cor_top10_up2)
# ------------------------------------------------------------------------------

# Select the top 10 most downregulated genes with 3'UTR shortening, and further filter for gene-miRNA interactions invovling upregulated miRNAs.

# 1. Filter miRNAs for those that are significantly upregulated (p-value < 0.05 and LFC > 0)
dysregulated_miRNAs_up <- miRNA_de %>%
  filter(`p-value` < 0.05, LFC > 0) %>%
  pull(miRNA)

# 2. Filter gene–miRNA interactions to include only those with upregulated miRNAs
dys_interactions <- df_interactions_expr %>% 
  filter(miRNA %in% dysregulated_miRNAs_up)

# 3. Summarize additional metrics per gene using these dysregulated interactions
gene_miRNA_expr_dys <- dys_interactions %>%
  group_by(Gene) %>%
  summarise(
    avg_miRNA_expr_dys = mean(Counts, na.rm = TRUE),
    n_interactions_dys = n()
  )

# 4. Merge with gene-level APA/DESeq2 data (merged_gene_data_final)
merged_gene_data_dys2 <- merged_gene_data_final %>%
  filter(gene_symbol %in% gene_miRNA_expr_dys$Gene)

# 5. For downregulated genes we typically expect 3'UTR lengthening (which yields binding site gain).
#    Select the top 10 most downregulated genes (lowest LFC) from those with 3'UTR lengthening.
top10_downregulated2 <- merged_gene_data_dys2 %>%
  filter(APA_type == "3'UTR lengthening") %>%
  slice_min(order_by = LFC, n = 10)

# 6. Create the scatter plot (using sum_diff as binding site gain) - !!!!!!!!!!!!!!!THIS PLOT DOES NOT WORK BECAUSE THERE ARE 0 ELIGIBLE GENES!!!!!!!!!!!!
p12 <- ggplot(top10_downregulated2, aes(x = sum_diff, y = LFC)) +
  geom_point(alpha = 0.7, size = 3, color = "darkred") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Binding Site Gain vs mRNA LFC",
       subtitle = "Top 10 Most Downregulated Genes (targeted by upregulated miRNAs)",
       x = "Binding Sites Gained (for upregulated miRNAs only)",
       y = "mRNA log2 Fold Change\n(<0: downregulation     >0: upregulation)") +
  theme_minimal(base_size = 18)

print(p12)

# ------------------------------------------------------------------------------
# Compute the Spearman correlation (related to plot p12) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THIS WILL NOT WORK BECAUSE OF 0 ELIGIBLE GENES!!!!!!!!!!!!!
cor_top10_down2 <- cor.test(top10_downregulated2$sum_diff, top10_downregulated2$LFC, method = "spearman")
print(cor_top10_down2)
# ------------------------------------------------------------------------------


# ------------------------
# Print summary of results
# ------------------------

cat("Number of dysregulated genes:", length(dysregulated_genes), "\n")
cat("Number of upregulated dysregulated genes:", num_upregulated, "(", (num_upregulated / length(dysregulated_genes))*100, "%)\n")
cat("Number of downregulated dysregulated genes:", num_downregulated, "(", (num_downregulated / length(dysregulated_genes))*100, "%)\n")

cat("Number of APA-regulated genes:", total_APA, "\n")
cat("Number of APA-regulated genes with 3'UTR lengthening:", num_UP, "(", round(perc_UP, 2), "%)\n")
cat("Number of APA-regulated genes with 3'UTR shortening:", num_DN, "(", round(perc_DN, 2), "%)\n")

cat("Number of dysregulated genes that are also APA-regulated genes:", num_overlap, "\n")

cat("Number of mmu miRNAs available on miRBase:", length(mouse_mmu_miRNAs), "\n")

cat("Number of predicted miRNAs (according to MultiMiR) targeting genes that are BOTH dysregulated and APA-regulated:", length(predicted_miRNAs), "\n")
cat("Number of miRNAs predicted (by MultiMiR) to target genes (that are BOTH dysregulated and APA-regulated) that also exist on miRBase AND detected in our data:", num_target_detected, "(", num_target_detected/length(target_miRNAs)*100, "%)\n")
cat("Number of miRNAs from MultiMiR (BEFORE filtering based on miRBase) that are detected in our samples:", num_predicted_detected, "(", num_predicted_detected/length(predicted_miRNAs)*100, "%)\n")

cat("Number of miRNAs detected in our samples that are also dysregulated:", common_dys_target_miRNAs, "\n")

cat("Total genes in UTRdbraw (cUTR) [i.e., total genes in genome known to be APA-regulated]:", length(all_genes), "\n")

cat("Number of gene-miRNA interactions according to MultiMiR:", nrow(df_long), "\n")
cat("Unique number of genes in gene-miRNA interactions:", length(unique(df_long$Gene)), "\n")
cat("Unique number of miRNAs in gene-miRNA interactions:", length(unique(df_long$miRNA)), "\n")

cat("Number of IMPOSSIBLE gene-miRNA interactions (both short and long isoforms have zero binding sites):", nrow(df_long)-nrow(df_long_filtered), "(", ((nrow(df_long)-nrow(df_long_filtered)) / nrow(df_long))*100, "%)\n")

cat("Number of LIKELY possible gene-miRNA interactions (at least one isoform has at least one binding site):", nrow(df_long_filtered), "(", nrow(df_long_filtered) / nrow(df_long)*100, "%)\n")
cat("Unique number of genes in gene-miRNA interactions:", length(unique(df_long_filtered$Gene)), "\n")
cat("Unique number of miRNAs in gene-miRNA interactions:", length(unique(df_long_filtered$miRNA)), "\n")

cat("Number of ACTUALLY possible gene-miRNA interactions because the miRNA was detected in our samples:", nrow(df_long_filtered_detected), "\n")
cat("Unique number of genes in ACTUALLY possible gene-miRNA interactions:", length(unique(df_long_filtered_detected$Gene)), "\n")
cat("Unique number of miRNAs in ACTUALLY possible gene-miRNA interactions:", length(unique(df_long_filtered_detected$miRNA)), "\n")


# ------------------------------------------------------------------------------
# Save the workspace
save.image(file = "Post-APAlyzer analysis - TDP 12m.RData")
