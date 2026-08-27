#!/usr/bin/env Rscript
# novel_variants_analysis_corrected.R

library(VennDiagram)
library(UpSetR)
library(dplyr)
library(readr)
library(ggplot2)

# Set working directory
setwd("/mnt/genomics/pilot_PPR/uniquevariants/extract_csq_tsv/plots_populationdb/R_analysis")

# Create output directory
dir.create(".", recursive = TRUE, showWarnings = FALSE)

# Function to create variant ID
create_variant_id <- function(df) {
  paste(df$CHROM, df$POS, df$REF, df$ALT, sep = ":")
}

# Function to load and process data for one dataset
load_dataset <- function(dataset_name, base_dir) {
  cat("Loading data for", dataset_name, "...\n")
  
  # Load all files
  not_gnomad <- read_csv(file.path(base_dir, paste0(dataset_name, "_NOT_in_gnomAD.csv")), show_col_types = FALSE)
  not_1000g <- read_csv(file.path(base_dir, paste0(dataset_name, "_NOT_in_1000G.csv")), show_col_types = FALSE)
  not_dbsnp <- read_csv(file.path(base_dir, paste0(dataset_name, "_NOT_in_dbSNP.csv")), show_col_types = FALSE)
  
  in_gnomad <- read_csv(file.path(base_dir, paste0(dataset_name, "_in_gnomAD.csv")), show_col_types = FALSE)
  in_1000g <- read_csv(file.path(base_dir, paste0(dataset_name, "_in_1000G.csv")), show_col_types = FALSE)
  in_dbsnp <- read_csv(file.path(base_dir, paste0(dataset_name, "_in_dbSNP.csv")), show_col_types = FALSE)
  
  # Create variant IDs
  not_gnomad_ids <- create_variant_id(not_gnomad)
  not_1000g_ids <- create_variant_id(not_1000g)
  not_dbsnp_ids <- create_variant_id(not_dbsnp)
  
  in_gnomad_ids <- create_variant_id(in_gnomad)
  in_1000g_ids <- create_variant_id(in_1000g)
  in_dbsnp_ids <- create_variant_id(in_dbsnp)
  
  # Calculate total variants
  total_variants <- length(in_gnomad_ids) + length(not_gnomad_ids)
  
  # Find variants absent from all databases
  absent_variants <- intersect(intersect(not_gnomad_ids, not_1000g_ids), not_dbsnp_ids)
  
  # Print statistics
  cat("  Total variants:", format(total_variants, big.mark = ","), "\n")
  cat("  NOT in gnomAD:", format(length(not_gnomad_ids), big.mark = ","), "\n")
  cat("  NOT in 1000G:", format(length(not_1000g_ids), big.mark = ","), "\n")
  cat("  NOT in dbSNP:", format(length(not_dbsnp_ids), big.mark = ","), "\n")
  cat("  Absent from all databases:", format(length(absent_variants), big.mark = ","), "\n")
  
  # VERIFICATION: Check overlaps
  cat("  VERIFICATION:\n")
  cat("    Absent variants in gnomAD:", sum(absent_variants %in% in_gnomad_ids), "(should be 0)\n")
  cat("    Absent variants in 1000G:", sum(absent_variants %in% in_1000g_ids), "(should be 0)\n")
  cat("    Absent variants in dbSNP:", sum(absent_variants %in% in_dbsnp_ids), "(should be 0)\n")
  
  cat("\n")
  
  return(list(
    dataset_name = dataset_name,
    total_variants = total_variants,
    not_gnomad = not_gnomad_ids,
    not_1000g = not_1000g_ids,
    not_dbsnp = not_dbsnp_ids,
    in_gnomad = in_gnomad_ids,
    in_1000g = in_1000g_ids,
    in_dbsnp = in_dbsnp_ids,
    absent_variants = absent_variants
  ))
}

# Main analysis
cat("=== VARIANTS ABSENT FROM POPULATION DATABASES ANALYSIS ===\n\n")

base_dir <- "/mnt/genomics/pilot_PPR/uniquevariants/extract_csq_tsv/plots_frequencies"
datasets <- c("GATK_Illumina", "PPR_Illumina", "GATK_Aviti", "PPR_Aviti")

# Load all datasets
all_data <- list()
for (dataset in datasets) {
  all_data[[dataset]] <- load_dataset(dataset, base_dir)
}

# Create Venn diagrams for each dataset
create_venn_diagram <- function(data, output_file) {
  # Create proper title with "Specific"
  title_map <- list(
    "GATK_Illumina" = "GATK Illumina Specific",
    "PPR_Illumina" = "PPR Illumina Specific", 
    "GATK_Aviti" = "GATK Aviti Specific",
    "PPR_Aviti" = "PPR Aviti Specific"
  )
  
  venn.diagram(
    x = list(
      "NOT in gnomAD" = data$not_gnomad,
      "NOT in 1000G" = data$not_1000g,
      "NOT in dbSNP" = data$not_dbsnp
    ),
    category.names = c("NOT in\ngnomAD", "NOT in\n1000G", "NOT in\ndbSNP"),
    filename = output_file,
    output = TRUE,
    imagetype = "png",
    height = 3000,
    width = 3000,
    resolution = 300,
    compression = "lzw",
    lwd = 4,
    lty = 'blank',
    fill = c('#FF6B6B', '#4ECDC4', '#45B7D1'),
    cex = 2.5,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 2.0,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    rotation = 1,
    main = paste(title_map[[data$dataset_name]], "- Variants Absent from Databases"),
    main.cex = 2.0,
    main.fontface = "bold"
  )
  
  cat("Created Venn diagram:", output_file, "\n")
}

# Create Venn diagrams for each dataset
for (dataset in datasets) {
  output_file <- paste0("venn_", dataset, "_specific.png")
  create_venn_diagram(all_data[[dataset]], output_file)
}

# Create comparison Venn diagrams
cat("\n=== PLATFORM COMPARISON VENN DIAGRAMS ===\n")

gatk_illumina_absent <- all_data[["GATK_Illumina"]]$absent_variants
gatk_aviti_absent <- all_data[["GATK_Aviti"]]$absent_variants

cat("GATK comparison:\n")
cat("  GATK Illumina Specific absent:", format(length(gatk_illumina_absent), big.mark = ","), "\n")
cat("  GATK Aviti Specific absent:", format(length(gatk_aviti_absent), big.mark = ","), "\n")
cat("  Common:", format(length(intersect(gatk_illumina_absent, gatk_aviti_absent)), big.mark = ","), "\n")

venn.diagram(
  x = list(
    "GATK Illumina\nSpecific" = gatk_illumina_absent,
    "GATK Aviti\nSpecific" = gatk_aviti_absent
  ),
  category.names = c("GATK Illumina\nSpecific", "GATK Aviti\nSpecific"),
  filename = "venn_GATK_platform_comparison.png",
  output = TRUE,
  imagetype = "png",
  height = 2500,
  width = 2500,
  resolution = 300,
  compression = "lzw",
  lwd = 4,
  lty = 'blank',
  fill = c('#0173B2', '#029E73'),
  cex = 2.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 2.0,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  main = "GATK: Variants Absent from Population Databases",
  main.cex = 2.0,
  main.fontface = "bold"
)

ppr_illumina_absent <- all_data[["PPR_Illumina"]]$absent_variants
ppr_aviti_absent <- all_data[["PPR_Aviti"]]$absent_variants

cat("PPR comparison:\n")
cat("  PPR Illumina Specific absent:", format(length(ppr_illumina_absent), big.mark = ","), "\n")
cat("  PPR Aviti Specific absent:", format(length(ppr_aviti_absent), big.mark = ","), "\n")
cat("  Common:", format(length(intersect(ppr_illumina_absent, ppr_aviti_absent)), big.mark = ","), "\n")

venn.diagram(
  x = list(
    "PPR Illumina\nSpecific" = ppr_illumina_absent,
    "PPR Aviti\nSpecific" = ppr_aviti_absent
  ),
  category.names = c("PPR Illumina\nSpecific", "PPR Aviti\nSpecific"),
  filename = "venn_PPR_platform_comparison.png",
  output = TRUE,
  imagetype = "png",
  height = 2500,
  width = 2500,
  resolution = 300,
  compression = "lzw",
  lwd = 4,
  lty = 'blank',
  fill = c('#DE8F05', '#CC78BC'),
  cex = 2.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 2.0,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  main = "PPR: Variants Absent from Population Databases",
  main.cex = 2.0,
  main.fontface = "bold"
)

# Save summary statistics
summary_stats <- data.frame(
  Dataset = datasets,
  Total_Variants = sapply(all_data, function(x) x$total_variants),
  NOT_in_gnomAD = sapply(all_data, function(x) length(x$not_gnomad)),
  NOT_in_1000G = sapply(all_data, function(x) length(x$not_1000g)),
  NOT_in_dbSNP = sapply(all_data, function(x) length(x$not_dbsnp)),
  Absent_from_All = sapply(all_data, function(x) length(x$absent_variants)),
  Absent_Percentage = sapply(all_data, function(x) round(length(x$absent_variants) / x$total_variants * 100, 2))
)

write.csv(summary_stats, "R_analysis_summary.csv", row.names = FALSE)
cat("Created summary: R_analysis_summary.csv\n")

# Save absent variants lists
for (dataset in datasets) {
  absent_variants <- all_data[[dataset]]$absent_variants
  if (length(absent_variants) > 0) {
    # Convert back to CHROM, POS, REF, ALT
    variant_parts <- strsplit(absent_variants, ":")
    absent_df <- data.frame(
      CHROM = sapply(variant_parts, function(x) x[1]),
      POS = as.integer(sapply(variant_parts, function(x) x[2])),
      REF = sapply(variant_parts, function(x) x[3]),
      ALT = sapply(variant_parts, function(x) x[4])
    )
    write.csv(absent_df, paste0("R_", dataset, "_absent_variants.csv"), row.names = FALSE)
    cat("Created absent variants list:", paste0("R_", dataset, "_absent_variants.csv"), "\n")
  }
}

# Final verification: Check common absent variants across all datasets
common_all <- Reduce(intersect, list(
  all_data[["GATK_Illumina"]]$absent_variants,
  all_data[["PPR_Illumina"]]$absent_variants,
  all_data[["GATK_Aviti"]]$absent_variants,
  all_data[["PPR_Aviti"]]$absent_variants
))

cat("\n=== FINAL VERIFICATION ===\n")
cat("Common variants absent from all databases across ALL datasets:", format(length(common_all), big.mark = ","), "\n")

if (length(common_all) > 0) {
  variant_parts <- strsplit(common_all, ":")
  common_df <- data.frame(
    CHROM = sapply(variant_parts, function(x) x[1]),
    POS = as.integer(sapply(variant_parts, function(x) x[2])),
    REF = sapply(variant_parts, function(x) x[3]),
    ALT = sapply(variant_parts, function(x) x[4])
  )
  write.csv(common_df, "R_common_absent_all_datasets.csv", row.names = FALSE)
  cat("Created common variants list: R_common_absent_all_datasets.csv\n")
} else {
  cat("No common variants found across all datasets.\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("- Individual Venn diagrams: venn_*_specific.png\n")
cat("- Platform comparison Venn diagrams: venn_*_platform_comparison.png\n")
cat("- Summary statistics: R_analysis_summary.csv\n")
cat("- Absent variants lists: R_*_absent_variants.csv\n")
if (length(common_all) > 0) {
  cat("- Common variants: R_common_absent_all_datasets.csv\n")
}
