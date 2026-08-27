#!/usr/bin/env Rscript

library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {
  cat("Usage: Rscript compare_happy_4pipelines.R <pipeline1_name> <pipeline1_dir> <pipeline2_name> <pipeline2_dir> <pipeline3_name> <pipeline3_dir> <pipeline4_name> <pipeline4_dir>\n")
  quit(status = 1)
}

pipeline1_name <- args[1]; pipeline1_dir <- args[2]
pipeline2_name <- args[3]; pipeline2_dir <- args[4]
pipeline3_name <- args[5]; pipeline3_dir <- args[6]
pipeline4_name <- args[7]; pipeline4_dir <- args[8]

for (i in c(1,3,5,7)) {
  d <- args[i+1]
  if (!dir.exists(d)) { cat("directory not found:", d, "\n"); quit(status = 1) }
}

pipeline1_safe <- gsub("[^A-Za-z0-9]", "_", pipeline1_name)
pipeline2_safe <- gsub("[^A-Za-z0-9]", "_", pipeline2_name)
pipeline3_safe <- gsub("[^A-Za-z0-9]", "_", pipeline3_name)
pipeline4_safe <- gsub("[^A-Za-z0-9]", "_", pipeline4_name)

output_prefix <- paste0(pipeline1_safe, "_vs_", pipeline2_safe, "_vs_",
                        pipeline3_safe, "_vs_", pipeline4_safe)

bed_region_map <- c(
  "15GC85" = "15GC85", "25GC65" = "25GC65", "ACMG" = "ACMG",
  "AD" = "AD", "BB" = "BB", "LM" = "LM", "OD" = "OD",
  "OMIM" = "OMIM", "REFSEQ" = "REFSEQ", "TRHP" = "TRHP"
)

region_order <- c("15GC85", "25GC65", "ACMG", "AD", "BB",
                  "LM", "OD", "OMIM", "REFSEQ", "TRHP")

read_happy_summaries <- function(base_dir, pipeline_name) {
  cat("reading", pipeline_name, "from", base_dir, "\n")
  regions <- list.dirs(path = base_dir, full.names = FALSE, recursive = FALSE)
  regions <- regions[regions != ""]
  summary_list <- list()
  for (region in regions) {
    file_path <- file.path(base_dir, region, paste0(region, ".summary.csv"))
    if (file.exists(file_path)) {
      df <- read.csv(file_path)
      df$Region <- region
      df$Pipeline <- pipeline_name
      summary_list[[region]] <- df
    } else {
      warning("file not found: ", file_path)
    }
  }
  all_data <- bind_rows(summary_list)
  cat(" ", nrow(all_data), "rows from", length(summary_list), "regions\n")
  return(all_data)
}

data1 <- read_happy_summaries(pipeline1_dir, pipeline1_name)
data2 <- read_happy_summaries(pipeline2_dir, pipeline2_name)
data3 <- read_happy_summaries(pipeline3_dir, pipeline3_name)
data4 <- read_happy_summaries(pipeline4_dir, pipeline4_name)

combined_data <- bind_rows(data1, data2, data3, data4)
combined_data$Region_Label <- recode(combined_data$Region, !!!bed_region_map)
combined_data$Region_Label <- factor(combined_data$Region_Label, levels = region_order)

plot_data <- combined_data %>%
  select(Region_Label, Pipeline, Type, Filter,
         Precision = METRIC.Precision,
         Recall = METRIC.Recall,
         F1_Score = METRIC.F1_Score) %>%
  pivot_longer(cols = c(Precision, Recall, F1_Score),
               names_to = "Metric", values_to = "Value")

pipeline_names  <- unique(combined_data$Pipeline)
pipeline_colors <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"), pipeline_names)

plot_metric <- function(variant_type, metric_name, data) {
  df <- data %>% filter(Type == variant_type, Metric == metric_name)

  ggplot(df, aes(x = Region_Label, y = Value,
                 color = Pipeline,
                 group = interaction(Pipeline, Filter),
                 linetype = Filter)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 3.5) +
    scale_color_manual(values = pipeline_colors) +
    scale_linetype_manual(values = c("ALL" = "solid", "PASS" = "dashed")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(expand = expansion(mult = c(0.03, 0.03))) +
    labs(title = metric_name, x = NULL, y = "Score",
         color = "Pipeline", linetype = "Filter") +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1, size = 18,
                                  color = "black", face = "bold"),
      axis.text.y  = element_text(size = 12),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right",
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      plot.margin  = margin(10, 10, 10, 10),
      panel.spacing = unit(0.2, "lines")
    )
}

create_combined_plot <- function(variant_type, data) {
  p1 <- plot_metric(variant_type, "Precision", data)
  p2 <- plot_metric(variant_type, "Recall",    data)
  p3 <- plot_metric(variant_type, "F1_Score",  data)

  (p1 / p2 / p3) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = paste0(variant_type, " Performance Comparison: HG002 vs GIAB v5.0q Truth Set"),
      subtitle = paste0("hap.py vcfeval results — ",
                        pipeline1_name, " \u2022 ", pipeline2_name, " \u2022 ",
                        pipeline3_name, " \u2022 ", pipeline4_name,
                        "\nSolid lines = ALL variants | Dashed lines = PASS variants"),
      theme = theme(
        plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 20))
      )
    )
}

cat("generating SNP plot...\n")
snp_plot <- create_combined_plot("SNP", plot_data)
ggsave(paste0("HG002_SNP_comparison_",   output_prefix, ".pdf"), snp_plot,  width = 12, height = 14)
ggsave(paste0("HG002_SNP_comparison_",   output_prefix, ".png"), snp_plot,  width = 12, height = 14, dpi = 300)

cat("generating INDEL plot...\n")
indel_plot <- create_combined_plot("INDEL", plot_data)
ggsave(paste0("HG002_INDEL_comparison_", output_prefix, ".pdf"), indel_plot, width = 12, height = 14)
ggsave(paste0("HG002_INDEL_comparison_", output_prefix, ".png"), indel_plot, width = 12, height = 14, dpi = 300)

# comparison table
comparison_table <- plot_data %>%
  group_by(Region_Label, Type, Filter, Metric, Pipeline) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Pipeline, values_from = Value) %>%
  mutate(
    Diff_P2_vs_P1 = .data[[pipeline2_name]] - .data[[pipeline1_name]],
    Diff_P3_vs_P1 = .data[[pipeline3_name]] - .data[[pipeline1_name]],
    Diff_P4_vs_P1 = .data[[pipeline4_name]] - .data[[pipeline1_name]],
    Diff_P3_vs_P2 = .data[[pipeline3_name]] - .data[[pipeline2_name]],
    Diff_P4_vs_P2 = .data[[pipeline4_name]] - .data[[pipeline2_name]],
    Diff_P4_vs_P3 = .data[[pipeline4_name]] - .data[[pipeline3_name]]
  ) %>%
  arrange(Type, Filter, Region_Label, Metric)

table_name <- paste0("comparison_", output_prefix, "_happy.csv")
write.csv(comparison_table, table_name, row.names = FALSE)

# summary stats
summary_stats <- comparison_table %>%
  group_by(Type, Metric) %>%
  summarise(
    !!paste0("Mean_", pipeline1_safe) := mean(.data[[pipeline1_name]], na.rm = TRUE),
    !!paste0("Mean_", pipeline2_safe) := mean(.data[[pipeline2_name]], na.rm = TRUE),
    !!paste0("Mean_", pipeline3_safe) := mean(.data[[pipeline3_name]], na.rm = TRUE),
    !!paste0("Mean_", pipeline4_safe) := mean(.data[[pipeline4_name]], na.rm = TRUE),
    Mean_Diff_P2_vs_P1 = mean(Diff_P2_vs_P1, na.rm = TRUE),
    Mean_Diff_P3_vs_P1 = mean(Diff_P3_vs_P1, na.rm = TRUE),
    Mean_Diff_P4_vs_P1 = mean(Diff_P4_vs_P1, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)

# best performer
best_performers <- plot_data %>%
  group_by(Type, Filter, Metric, Region_Label) %>%
  slice_max(Value, n = 1, with_ties = TRUE) %>%
  count(Pipeline, name = "Best_Count") %>%
  arrange(desc(Best_Count))
print(best_performers)

cat("done.\n")
cat("SNP:   HG002_SNP_comparison_",   output_prefix, ".pdf/png\n", sep="")
cat("INDEL: HG002_INDEL_comparison_", output_prefix, ".pdf/png\n", sep="")
cat("table:", table_name, "\n")

