# APAlyzer on TDP 6m files
# The files were aligned to mm10, so this is the reference genome that is used throughout this script
# ------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())


# Install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("APAlyzer")


# Load libraries
library(APAlyzer)
library(repmis)
library(stats)
library(ggplot2)
library(dplyr)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)

# Define BAM filenames
# These are the three non-transgenic 6-month-old mice and the three M337V transgenic 6-month-old mice.
# These are the ACTUAL names of the files. If they don't match, the code will not work.
bam1 <- "T2.bam"
bam2 <- "T3.bam"
bam3 <- "T4.bam"
bam4 <- "N1.bam"
bam5 <- "N2.bam"
bam6 <- "N3.bam"

flsall <- c(bam1, bam2, bam3, bam4, bam5, bam6)
names(flsall) <- c("T2", "T3", "T4", "N1", "N2", "N3")


# Define reference genome
URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
file="mm10_REF.RData"
source_data(paste0(URL,file,"?raw=True"))

# Building 3’UTR and intronic PAS reference region at once
PASREF=REF4PAS(refUTRraw,dfIPA,dfLE)
UTRdbraw=PASREF$UTRdbraw
dfIPA=PASREF$dfIPA
dfLE=PASREF$dfLE


# -------------------------
# Analysis of APA in 3’UTRs
# -------------------------

# Building aUTR and cUTR references
UTRdbraw=REF3UTR(refUTRraw)

# Calculation of relative expression
DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="forward")

# Export relative expression results
write.table(DFUTRraw, file = "DFUTRraw_TDP_6m.csv", sep = ",", row.names = FALSE, col.names = TRUE)


# --------------------------
# Analysis of APA in introns
# --------------------------

# Calculation of relative expression
IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="forward", nts=1)

# Export relative expression results
write.table(IPA_OUTraw, file = "IPA_OUTraw_TDP_6m.csv", sep = ",", row.names = FALSE, col.names = TRUE)


# -----------------------------------
# Significance analysis of APA events
# -----------------------------------

# Build the sample table with replicates
sampleTable = data.frame(samplename = c(names(flsall)),
                         condition = c(rep("Tg_6m",3),rep("NTg_6m",3)))
sampleTable

# Significantly regulated APA in 3’UTRs
test_3UTR=APAdiff(sampleTable,
                  DFUTRraw, 
                  conKET='NTg_6m',
                  trtKEY='Tg_6m',
                  PAS='3UTR',
                  CUTreads=10,
                  p_adjust_methods="fdr",
                  MultiTest='unpaired t-test')

table(test_3UTR$APAreg)

# Get stats for each group
medians <- aggregate(RED ~ APAreg, data = test_3UTR, FUN = median)
print(medians)

iqrs <- aggregate(RED ~ APAreg, data = test_3UTR, FUN = IQR)
print(iqrs)

# Significantly regulated APA in introns
test_IPA=APAdiff(sampleTable,
                 IPA_OUTraw, 
                 conKET='NTg_6m',
                 trtKEY='Tg_6m',
                 PAS='IPA',
                 CUTreads=10,
                 p_adjust_methods="fdr",
                 MultiTest='unpaired t-test')

table(test_IPA$APAreg)

# Export t-test results
write.table(test_3UTR, file = "test_3UTR_TDP_6m.csv", sep = ",", row.names = FALSE, col.names = TRUE)
write.table(test_IPA, file = "test_IPA_TDP_6m.csv", sep = ",", row.names = FALSE, col.names = TRUE)


# ---------------------------------
# Visualization of analysis results
# ---------------------------------

# 3UTR
# Volcano plot
v <- APAVolcano(test_3UTR,
                PAS='3UTR',
                Pcol = "pvalue",
                xlab = 'Relative Expression Difference (RED)')
v + scale_fill_manual(values = c("blue", "gray", "red")) +
  theme(legend.position = "none")

# Boxplot
test_3UTR$APAreg <- factor(test_3UTR$APAreg,
                           levels = c("DN", "NC", "UP"),
                           labels = c("Negative", "No change", "Positive"))

# Calculate the number of observations (n) for each APA regulation group
counts <- test_3UTR %>% 
  group_by(APAreg) %>% 
  summarise(n = n())

# Create the boxplot using APABox
b <- APABox(test_3UTR, 
            xlab = "APA regulation", 
            ylab = "Relative Expression Difference (RED)", 
            plot_title = NULL)

# Determine a suitable y position for the labels.
# Here, we calculate the maximum RED for each group and then add a small offset.
max_vals <- test_3UTR %>% 
  group_by(APAreg) %>% 
  summarise(ymax = max(RED))

counts <- left_join(counts, max_vals, by = "APAreg")

# Add the text annotations (n-numbers) to the plot
b + scale_fill_manual(values = c("blue", "gray", "red")) +
  theme(legend.position = "none") +
  geom_text(data = counts, 
            aes(x = APAreg, y = ymax + 0.3, label = paste("n =", n)), 
            vjust = 0)

# IPA
APAVolcano(test_IPA, PAS='IPA', Pcol = "pvalue", top=5)
APABox(test_IPA, xlab = "APAreg", ylab = "RED", plot_title = NULL)

# Both
test_3UTR$APA="3'UTR"
test_IPA$APA="IPA"
dfplot=rbind(test_3UTR[,c('RED','APA')],test_IPA[,c('RED','APA')])

###violin
ggplot(dfplot, aes(x = APA, y = RED)) + 
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.2)+ theme_bw() + 
  geom_hline(yintercept=0, linetype="dashed", color = "red")

###CDF
ggplot(dfplot, aes( x = RED, color = APA)) + 
  stat_ecdf(geom = "step") +
  ylab("cumulative fraction")+ 
  geom_vline(xintercept=0, linetype="dashed", color = "gray")+ theme_bw() + 
  geom_hline(yintercept=0.5, linetype="dashed", color = "gray")

# ------------------------------------------------------------------------------
# Save gene sequences to FASTA files
# ------------------------------------------------------------------------------

# Subset significantly APA‐regulated genes (UP or DN)
sigGenes <- subset(test_3UTR, APAreg %in% c("UP", "DN"))
geneIDs <- sigGenes$gene_symbol  # Ensure this column matches your test_3UTR data

# Function to extract isoform sequence using UTRdbraw names
# Assumes UTRdbraw is a GRangesList with names like "GeneSymbol__cUTR" or "GeneSymbol__aUTR"
getIsoformSeq <- function(gene, isoform = c("short", "long")) {
  isoform <- match.arg(isoform)
  # Find all entries in UTRdbraw for this gene (names starting with the gene symbol and '__')
  idx <- grep(paste0("^", gene, "__"), names(UTRdbraw))
  if (length(idx) == 0) {
    message("No UTR entries for gene: ", gene)
    return(DNAStringSet())
  }
  # For short isoform, select only the cUTR; for long, combine cUTR and aUTR if available.
  if (isoform == "short") {
    sub_idx <- idx[grepl("__cUTR$", names(UTRdbraw)[idx])]
  } else {  # long isoform: get both cUTR and aUTR entries
    sub_idx <- idx[grepl("__(cUTR|aUTR)$", names(UTRdbraw)[idx])]
  }
  if (length(sub_idx) == 0) {
    message("No UTR entries for gene: ", gene, " for isoform ", isoform)
    return(DNAStringSet())
  }
  # Combine the GRanges objects into one (reducing overlapping/adjacent ranges)
  gr_combined <- reduce(unlist(UTRdbraw[sub_idx]))
  # Extract the genomic sequence from the mm10 genome
  seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr_combined)
  return(seq)
}

# Create empty DNAStringSet objects to store the sequences
shortSeqs <- DNAStringSet()
longSeqs <- DNAStringSet()

# Loop over each significant gene to extract sequences
for (g in geneIDs) {
  shortSeq <- getIsoformSeq(g, "short")
  longSeq <- getIsoformSeq(g, "long")
  
  # Check that sequences were found for both isoforms
  if (length(shortSeq) == 0 || length(longSeq) == 0) {
    message("Skipping gene ", g, " due to missing sequence data.")
    next
  }
  
  # Name the sequences with the gene symbol
  names(shortSeq) <- g
  names(longSeq) <- g
  
  # Append to the DNAStringSet objects
  shortSeqs <- append(shortSeqs, shortSeq)
  longSeqs <- append(longSeqs, longSeq)
}

# Write the sequences to FASTA files if any sequences were retrieved
if (length(shortSeqs) > 0) {
  writeXStringSet(shortSeqs, filepath = "TDP_6m_Short_isoforms_APA_reg.fasta")
} else {
  message("No short isoform sequences found.")
}

if (length(longSeqs) > 0) {
  writeXStringSet(longSeqs, filepath = "TDP_6m_Long_isoforms_APA_reg.fasta")
} else {
  message("No long isoform sequences found.")
}


# ------------------------------------------------------------------------------
# Save the workspace
save.image(file = "APAlyzer_TDP_6m_workspace.RData")
