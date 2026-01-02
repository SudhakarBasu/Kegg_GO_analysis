# KEGG Enrichment Analysis with Publication-Quality Plots
# Analyzes KEGG pathways and creates visualizations

library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(enrichplot)

cat("KEGG Enrichment Analysis\n")
cat(rep("=", 70), "\n\n", sep="")

# ============================================================================
# Read your data
# ============================================================================
cat("Reading KEGG annotation data...\n")
kegg_data <- read.table("Kegg_biotic.txt", header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE)

cat("Total rows in file:", nrow(kegg_data), "\n")
cat("Unique KEGG IDs:", length(unique(kegg_data$Kegg_id)), "\n\n")

# Display sample of data
cat("Sample of your data:\n")
print(head(kegg_data, 10))
cat("\n")

# ============================================================================
# KEGG Enrichment Analysis
# ============================================================================
cat(rep("=", 70), "\n", sep="")
cat("Running KEGG Enrichment Analysis\n")
cat(rep("=", 70), "\n\n")

# Get unique KO IDs
ko_ids <- unique(kegg_data$Kegg_id)

cat("Analyzing", length(ko_ids), "unique KEGG orthologs...\n\n")

# Try with 'ko' organism (KEGG Orthology)
kegg_result_ko <- tryCatch({
  enrichKEGG(
    gene = ko_ids,
    organism = 'ko',
    keyType = 'kegg',
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2
  )
}, error = function(e) {
  cat("Error with 'ko' organism:", e$message, "\n")
  return(NULL)
})

# ============================================================================
# Process and visualize results
# ============================================================================

if (!is.null(kegg_result_ko) && nrow(kegg_result_ko) > 0) {
  cat("\nSUCCESS! Found", nrow(kegg_result_ko), "enriched KEGG pathways\n\n")
  
  # Show top 10 results
  cat("Top 10 enriched pathways:\n")
  result_df <- as.data.frame(kegg_result_ko)
  print(result_df[1:min(10, nrow(result_df)), 
                  c("ID", "Description", "GeneRatio", "pvalue", "p.adjust", "Count")])
  
  # Save full results
  write.csv(result_df, "KEGG_enrichment_results.csv", row.names = FALSE)
  cat("\nFull results saved to: KEGG_enrichment_results.csv\n\n")
  
  # ============================================================================
  # Create visualizations
  # ============================================================================
  cat(rep("=", 70), "\n", sep="")
  cat("Creating Publication-Quality Plots\n")
  cat(rep("=", 70), "\n\n")
  
  # Determine how many pathways to show
  n_show <- min(20, nrow(kegg_result_ko))
  cat("Creating plots for top", n_show, "pathways...\n")
  
  # 1. Dot Plot
  cat("\n1. Creating dot plot...\n")
  p_dot <- dotplot(kegg_result_ko, 
                   showCategory = n_show,
                   title = "KEGG Pathway Enrichment",
                   font.size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = 9))
  
  ggsave("KEGG_dotplot.pdf", p_dot, width = 10, height = 8)
  ggsave("KEGG_dotplot.png", p_dot, width = 10, height = 8, dpi = 300)
  cat("   Saved: KEGG_dotplot.pdf and KEGG_dotplot.png\n")
  
  # 2. Bar Plot
  cat("2. Creating bar plot...\n")
  p_bar <- barplot(kegg_result_ko, 
                   showCategory = n_show,
                   title = "KEGG Pathway Enrichment") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = 9))
  
  ggsave("KEGG_barplot.pdf", p_bar, width = 10, height = 8)
  ggsave("KEGG_barplot.png", p_bar, width = 10, height = 8, dpi = 300)
  cat("   Saved: KEGG_barplot.pdf and KEGG_barplot.png\n")
  
  # 3. Custom horizontal bar plot (like GO plot style)
  cat("3. Creating custom horizontal bar plot...\n")
  
  plot_data <- result_df[1:n_show, ]
  plot_data$log_pvalue <- -log10(plot_data$pvalue)
  
  # Truncate long pathway names
  plot_data$Description_short <- ifelse(nchar(plot_data$Description) > 50,
                                        paste0(substr(plot_data$Description, 1, 47), "..."),
                                        plot_data$Description)
  
  # Order by p-value for plotting
  plot_data$Description_short <- factor(plot_data$Description_short,
                                        levels = rev(plot_data$Description_short))
  
  p_custom <- ggplot(plot_data, aes(x = log_pvalue, y = Description_short)) +
    geom_bar(stat = "identity", fill = "#2E86AB", width = 0.7) +
    labs(title = "KEGG Pathway Enrichment Analysis",
         x = "-log10 (P-value)",
         y = NULL) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11, color = "black"),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black")
    ) +
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(0, max(plot_data$log_pvalue) * 1.1))
  
  ggsave("KEGG_horizontal_barplot.pdf", p_custom, width = 12, height = 10)
  ggsave("KEGG_horizontal_barplot.png", p_custom, width = 12, height = 10, dpi = 300)
  cat("   Saved: KEGG_horizontal_barplot.pdf and KEGG_horizontal_barplot.png\n")
  
  # 4. Enrichment Map (if sufficient pathways)
  if (nrow(kegg_result_ko) >= 5) {
    cat("4. Creating enrichment map...\n")
    tryCatch({
      kegg_result_pairwise <- pairwise_termsim(kegg_result_ko)
      p_emap <- emapplot(kegg_result_pairwise, 
                         showCategory = min(30, nrow(kegg_result_ko)),
                         cex_label_category = 0.6) +
        ggtitle("KEGG Pathway Network") +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      
      ggsave("KEGG_enrichment_map.pdf", p_emap, width = 12, height = 10)
      ggsave("KEGG_enrichment_map.png", p_emap, width = 12, height = 10, dpi = 300)
      cat("   Saved: KEGG_enrichment_map.pdf and KEGG_enrichment_map.png\n")
    }, error = function(e) {
      cat("   Could not create enrichment map:", e$message, "\n")
    })
  }
  
  # 5. Gene-Concept Network (cnetplot)
  if (nrow(kegg_result_ko) >= 3) {
    cat("5. Creating gene-concept network...\n")
    tryCatch({
      p_cnet <- cnetplot(kegg_result_ko, 
                         showCategory = min(5, nrow(kegg_result_ko)),
                         cex_label_gene = 0.5,
                         cex_label_category = 0.8) +
        ggtitle("Gene-Pathway Network") +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      
      ggsave("KEGG_cnetplot.pdf", p_cnet, width = 12, height = 10)
      ggsave("KEGG_cnetplot.png", p_cnet, width = 12, height = 10, dpi = 300)
      cat("   Saved: KEGG_cnetplot.pdf and KEGG_cnetplot.png\n")
    }, error = function(e) {
      cat("   Could not create gene-concept network:", e$message, "\n")
    })
  }
  
  # Display the main plots
  cat("\n\nDisplaying plots...\n")
  print(p_dot)
  
  # ============================================================================
  # Summary
  # ============================================================================
  cat("\n\n", rep("=", 70), "\n", sep="")
  cat("ANALYSIS COMPLETE\n")
  cat(rep("=", 70), "\n\n", sep="")
  
  cat("Summary:\n")
  cat("  Total enriched pathways:", nrow(kegg_result_ko), "\n")
  cat("  Pathways shown in plots:", n_show, "\n")
  cat("  Significance threshold: p < 0.05, q < 0.2\n\n")
  
  cat("Output files:\n")
  cat("  Data:\n")
  cat("    - KEGG_enrichment_results.csv\n\n")
  cat("  Plots:\n")
  cat("    - KEGG_dotplot.pdf/png (dot plot)\n")
  cat("    - KEGG_barplot.pdf/png (bar plot)\n")
  cat("    - KEGG_horizontal_barplot.pdf/png (custom style)\n")
  if (nrow(kegg_result_ko) >= 5) {
    cat("    - KEGG_enrichment_map.pdf/png (network)\n")
  }
  if (nrow(kegg_result_ko) >= 3) {
    cat("    - KEGG_cnetplot.pdf/png (gene-pathway network)\n")
  }
  
  cat("\nTop 5 pathways:\n")
  for (i in 1:min(5, nrow(result_df))) {
    cat("  ", i, ". ", result_df$Description[i], "\n", sep = "")
    cat("      p-value: ", formatC(result_df$pvalue[i], format = "e", digits = 2),
        ", Count: ", result_df$Count[i], "\n", sep = "")
  }
  
} else {
  cat("\n", rep("=", 70), "\n", sep="")
  cat("NO ENRICHMENT FOUND\n")
  cat(rep("=", 70), "\n\n", sep="")
  
  cat("No significantly enriched KEGG pathways were found.\n\n")
  
  cat("Possible reasons:\n")
  cat("1. Your KEGG IDs may not map to specific pathways\n")
  cat("2. The pathways may be too diverse (no clear enrichment)\n")
  cat("3. The p-value/q-value cutoffs may be too stringent\n\n")
  
  cat("Suggestions:\n")
  cat("1. Try relaxing the cutoffs:\n")
  cat("   pvalueCutoff = 0.1 (instead of 0.05)\n")
  cat("   qvalueCutoff = 0.3 (instead of 0.2)\n\n")
  cat("2. Check if your KEGG IDs are correct format (e.g., K00001)\n\n")
  cat("3. View the frequency distribution instead:\n")
  cat("   Run the frequency analysis script for basic counts\n")
}

cat("\nWorking directory:", getwd(), "\n")