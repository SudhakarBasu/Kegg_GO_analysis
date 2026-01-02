# KEGG and GO Term Frequency Analysis and Visualization
# Simple counting and plotting of terms

library(dplyr)
library(ggplot2)

cat("KEGG and GO Term Frequency Analysis\n")
cat(rep("=", 70), "\n\n", sep="")

# ============================================================================
# Read your data files (single column each)
# ============================================================================
kegg_data <- read.table("Kegg_biotic.txt", header = TRUE, stringsAsFactors = FALSE)
go_data <- read.table("GO_biotic.txt", header = TRUE, stringsAsFactors = FALSE)

cat("Input data summary:\n")
cat("  Total KEGG ID entries:", nrow(kegg_data), "\n")
cat("  Unique KEGG IDs:", length(unique(kegg_data$Kegg_id)), "\n")
cat("  Total GO term entries:", nrow(go_data), "\n")
cat("  Unique GO terms:", length(unique(go_data$GO_id)), "\n\n")

# ============================================================================
# Part 1: Count KEGG ID frequencies
# ============================================================================
cat(rep("=", 70), "\n", sep="")
cat("PART 1: KEGG ID FREQUENCY ANALYSIS\n")
cat(rep("=", 70), "\n\n", sep="")

# Count frequency of each KEGG ID
kegg_counts <- kegg_data %>%
  group_by(Kegg_id) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))

cat("Top 20 most frequent KEGG IDs:\n\n")
print(head(kegg_counts, 20))

# Save results
write.csv(kegg_counts, "KEGG_ID_frequency.csv", row.names = FALSE)
cat("\nSaved to: KEGG_ID_frequency.csv\n")

# ============================================================================
# Part 2: Count GO term frequencies
# ============================================================================
cat("\n\n", rep("=", 70), "\n", sep="")
cat("PART 2: GO TERM FREQUENCY ANALYSIS\n")
cat(rep("=", 70), "\n\n", sep="")

# Count frequency of each GO term
go_counts <- go_data %>%
  group_by(GO_id) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))

cat("Total unique GO terms:", nrow(go_counts), "\n")
cat("Top 20 most frequent GO terms:\n\n")
print(head(go_counts, 20))

# Save results
write.csv(go_counts, "GO_term_frequency.csv", row.names = FALSE)
cat("\nSaved to: GO_term_frequency.csv\n")

# ============================================================================
# Get GO term information and create visualization
# ============================================================================

# Check if GO.db is available
if (requireNamespace("GO.db", quietly = TRUE)) {
  library(GO.db)
  
  cat("\n\nAdding GO term descriptions...\n")
  
  # Get GO term information
  go_counts$Description <- sapply(go_counts$GO_id, function(x) {
    tryCatch(Term(x), error = function(e) x)
  })
  
  go_counts$Ontology <- sapply(go_counts$GO_id, function(x) {
    tryCatch(Ontology(x), error = function(e) NA)
  })
  
  # Remove terms without ontology info
  go_counts <- go_counts[!is.na(go_counts$Ontology), ]
  
  # Map ontology to full names
  go_counts$Ontology_full <- factor(go_counts$Ontology,
                                    levels = c("BP", "CC", "MF"),
                                    labels = c("Biological process",
                                               "Cellular component",
                                               "Molecular function"))
  
  cat("\nGO terms by category:\n")
  cat("  Biological Process (BP):", sum(go_counts$Ontology == "BP"), "\n")
  cat("  Cellular Component (CC):", sum(go_counts$Ontology == "CC"), "\n")
  cat("  Molecular Function (MF):", sum(go_counts$Ontology == "MF"), "\n\n")
  
  # Save full results with descriptions
  write.csv(go_counts, "GO_term_frequency_annotated.csv", row.names = FALSE)
  cat("Saved annotated results to: GO_term_frequency_annotated.csv\n\n")
  
  # ============================================================================
  # Create visualization
  # ============================================================================
  cat("Creating visualization...\n")
  
  # Select top 10 from each category
  plot_data <- go_counts %>%
    filter(!is.na(Ontology)) %>%
    group_by(Ontology) %>%
    arrange(desc(Count)) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # Truncate long descriptions
  plot_data$Label <- ifelse(nchar(plot_data$Description) > 45,
                            paste0(substr(plot_data$Description, 1, 42), "..."),
                            plot_data$Description)
  
  # Order for plotting (by ontology and count)
  plot_data <- plot_data %>%
    arrange(Ontology, desc(Count))
  
  plot_data$Label <- factor(plot_data$Label, levels = rev(plot_data$Label))
  
  # Define colors matching reference image
  ontology_colors <- c("Biological process" = "#8B3A3A",    # Dark red/brown
                       "Cellular component" = "#4682B4",      # Steel blue
                       "Molecular function" = "#3CB371")      # Medium sea green
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Count, y = Label, fill = Ontology_full)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = ontology_colors) +
    labs(title = "GO functional analysis",
         x = "Gene Count",
         y = NULL,
         fill = NULL) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11, color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black")
    ) +
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(0, max(plot_data$Count) * 1.1))
  
  # Save plot
  ggsave("GO_functional_analysis.pdf", p, width = 12, height = 10)
  ggsave("GO_functional_analysis.png", p, width = 12, height = 10, dpi = 300)
  
  cat("\nPlots saved:\n")
  cat("  - GO_functional_analysis.pdf\n")
  cat("  - GO_functional_analysis.png\n\n")
  
  # Display plot
  print(p)
  
  # Print summary of plot data
  cat("\nTerms included in plot:\n")
  for (ont in c("BP", "CC", "MF")) {
    ont_name <- c(BP = "Biological Process", 
                  CC = "Cellular Component", 
                  MF = "Molecular Function")[ont]
    ont_data <- plot_data[plot_data$Ontology == ont, ]
    
    if (nrow(ont_data) > 0) {
      cat("\n", ont_name, " (", nrow(ont_data), " terms):\n", sep = "")
      for (i in 1:nrow(ont_data)) {
        cat("  ", i, ". ", ont_data$Description[i], 
            " (Count: ", ont_data$Count[i], ")\n", sep = "")
      }
    }
  }
  
} else {
  cat("\n\nGO.db package not found.\n")
  cat("Installing GO.db for term descriptions...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("GO.db")
  cat("\nPlease run the script again after installation.\n")
}

# ============================================================================
# Summary
# ============================================================================
cat("\n\n", rep("=", 70), "\n", sep="")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n\n", sep="")

cat("Summary:\n")
cat("  - ", length(unique(kegg_data$Kegg_id)), " unique KEGG IDs analyzed\n", sep = "")
cat("  - ", length(unique(go_data$GO_id)), " unique GO terms analyzed\n", sep = "")
cat("\nOutput files:\n")
cat("  1. KEGG_ID_frequency.csv - KEGG ID counts\n")
cat("  2. GO_term_frequency.csv - GO term counts (basic)\n")
cat("  3. GO_term_frequency_annotated.csv - GO terms with descriptions\n")
cat("  4. GO_functional_analysis.pdf - Publication-quality figure\n")
cat("  5. GO_functional_analysis.png - High-res PNG (300 dpi)\n")
cat("\nWorking directory: ", getwd(), "\n", sep = "")