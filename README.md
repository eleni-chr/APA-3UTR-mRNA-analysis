# TDP-43 muscle APA analysis code

R scripts supporting the manuscript “Early Metabolic Insufficiency and Alternative Polyadenylation Remodelling in ALS-Linked TDP-43M337V Skeletal Muscle” (submitted for peer-review; link will be provided here upon publication; title may be updated).

Copyright © 2026 Dr Eleni Christoforidou. The code is licensed under the Apache License, Version 2.0; see LICENSE.

## Study context

The manuscript investigates early metabolic insufficiency and alternative polyadenylation (APA) remodelling in skeletal muscle from mice carrying the ALS-linked TDP-43M337V mutation. The R scripts here cover the age-specific APA analyses at 6, 9, and 12 months, including 3′UTR and intronic polyadenylation, and connect those results with differential mRNA and miRNA expression, miRNA-target predictions, and Gene Ontology (GO)-Slim visualisations. Sequencing quality control, trimming, alignment, read counting, differential-expression analysis, and generation of enrichment inputs were performed separately, as described in the manuscript.

## Repository contents

This code repository contains the eight R scripts listed below, the Apache 2.0 LICENSE file, and this README. It does not contain study data, reference sequence files, or analysis outputs. The R scripts expect any required input files to be placed directly in the same directory as the scripts. They write generated output files to that working directory.

## Analysis scripts

| Script | Purpose |
| --- | --- |
| APAlyzer on TDP 6m files.R, APAlyzer on TDP 9m files.R, APAlyzer on TDP 12m files.R | APA analysis of aligned mRNA-seq BAM files at each age, including 3′UTR and intronic polyadenylation analyses. |
| Post-APAlyzer analysis - TDP 6m.R, Post-APAlyzer analysis - TDP 9m.R, Post-APAlyzer analysis - TDP 12m.R | Cross-analysis of differential expression, APA and miRNA target sites. |
| Volcano_plots.R | Volcano plots for mRNA expression, miRNA expression and APA results. |
| GO_results_visualisation.R | Dot plots from PANTHER GO-Slim result spreadsheets. |

The manuscript describes FASTQ quality control, trimming, alignment, counting, differential-expression analysis and some enrichment analyses as work performed outside these R scripts (including on Galaxy or online services). Those upstream steps are not reimplemented here.

## Input dataset manifest

Download the study-specific inputs from the associated dataset publication when it becomes available. Place the required files directly alongside the R scripts, keeping the filenames below. These study files are not included in the code repository; they are intended to become publicly available upon publication of the manuscript through a separate dataset publication. The dataset DOI is forthcoming.

| Input files | Used by |
| --- | --- |
| T2.bam, T3.bam, T4.bam, N1.bam, N2.bam, N3.bam; T7.bam, T8.bam, T15.bam, N8.bam, N9.bam, N10.bam; T11.bam, T13.bam, T21.bam, N15.bam, N16.bam, N17.bam | The three APAlyzer scripts. These are 18 study BAM files in total. |
| DESeq2_mRNA_TDP_6m.csv, DESeq2_mRNA_TDP_9m.csv, DESeq2_mRNA_TDP_12m.csv; DESeq2_miRNA_TDP_6m.csv, DESeq2_miRNA_TDP_9m.csv, DESeq2_miRNA_TDP_12m.csv | Post-APAlyzer and volcano plot scripts. |
| PANTHER_GO-Slim_DESeq2_mRNA_6m.xlsx, PANTHER_GO-Slim_DESeq2_mRNA_9m.xlsx, PANTHER_GO-Slim_miRPath_miRNA_6m.xlsx, PANTHER_GO-Slim_miRPath_miRNA_9m.xlsx, PANTHER_GO-Slim_miRPath_miRNA_12m.xlsx, PANTHER_GO-Slim_APA_mRNA_6m.xlsx | GO_results_visualisation.R. |

The post-APAlyzer scripts also require mature.fa, a miRNA FASTA file. Download it from the miRBase downloads page at https://mirbase.org/download/ and place it alongside the scripts. Follow miRBase's citation guidance when using its sequences.

The APAlyzer scripts retrieve the mm10_REF.RData PAS reference from the external RJWANGbioinfo/PAS_reference_RData GitHub repository at runtime. This reference is not included in this repository. The post-APAlyzer analyses also query online miRNA-target and Ensembl annotation resources, so internet access is required and results may reflect updated external databases.

## Requirements

Use R with internet access for package installation and online reference lookups. Package versions were not recorded with the source scripts, so results may vary as R packages and online annotation services change.

Install the required packages once in R:

~~~r
install.packages(c(
  "BiocManager", "repmis", "ggplot2", "dplyr", "readr", "purrr",
  "tidyr", "broom", "tidyverse", "readxl", "stringr", "janitor",
  "forcats", "scales"
))

BiocManager::install(c(
  "APAlyzer", "Biostrings", "BSgenome.Mmusculus.UCSC.mm10",
  "GenomicRanges", "multiMiR", "biomaRt"
))
~~~

## Running the scripts

Use a working copy containing the code files and the input files listed in the manifest. Set R's working directory to that location, then run the three APAlyzer scripts first:

~~~r
source("APAlyzer on TDP 6m files.R")
source("APAlyzer on TDP 9m files.R")
source("APAlyzer on TDP 12m files.R")
~~~

These scripts generate APA result files and workspace snapshots in the working directory. Then download mature.fa from miRBase and run the age-matched post-APAlyzer scripts:

~~~r
source("Post-APAlyzer analysis - TDP 6m.R")
source("Post-APAlyzer analysis - TDP 9m.R")
source("Post-APAlyzer analysis - TDP 12m.R")
~~~

The APAlyzer scripts generate the test_3UTR result CSVs and APAlyzer workspace snapshots in the working directory. The post-APAlyzer scripts and Volcano_plots.R use the generated APA result CSVs; the post-APAlyzer scripts also load the generated workspace snapshots. These are intermediate outputs, not files included in the code repository. Run the APAlyzer scripts first for a complete workflow. Volcano_plots.R also uses the DESeq2 result tables. GO_results_visualisation.R uses the PANTHER spreadsheets and can run independently.

## Manuscript, dataset and code citations

When using these scripts, cite all three records once their publication details are available:

1. **Manuscript:** _Pending publication; link will be provided here upon publication._
2. **Dataset:** _Pending publication; link will be provided here upon publication._
3. **Code:** Christoforidou, E. (2026). Alternative Polyadenylation Analysis in ALS-Linked TDP-43M337V Mouse Skeletal Muscle: R Scripts (Version v1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24374857

## License and disclaimer

The R code is licensed under the Apache License, Version 2.0. See LICENSE for the full license text. The license applies to the code only; it does not license study data or third-party reference data.

The code is provided “AS IS”, without warranties or guarantees of any kind. To the extent permitted by law, the copyright holder accepts no responsibility or liability for use of the code, its outputs, or decisions made from them. See Sections 7 and 8 of the Apache License, Version 2.0.

Copyright © 2026 Dr Eleni Christoforidou.
