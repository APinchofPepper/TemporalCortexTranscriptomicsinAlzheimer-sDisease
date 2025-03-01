#!/usr/bin/env Rscript
# Simplified GSEA script without GO term conversion
# To run: Rscript basic-gsea.R

# Install required packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org")

# Try to load needed packages, install if missing
packages <- c("fgsea", "org.Hs.eg.db", "reactome.db", "dplyr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    if (pkg %in% c("org.Hs.eg.db", "reactome.db", "fgsea")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg, repos="https://cloud.r-project.org")
    }
  }
}

# Load libraries
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(reactome.db))
suppressPackageStartupMessages(library(dplyr))

# Read the ranked gene list 
if (file.exists("ranked_genes_for_gsea.csv")) {
  ranked_genes <- read.csv("ranked_genes_for_gsea.csv")
} else {
  stop("Cannot find ranked_genes_for_gsea.csv file. Run the Python pipeline first.")
}

message("Creating gene ranking...")
# Create a named vector of gene scores
ranks <- setNames(ranked_genes$score, ranked_genes$gene_symbol)

message("Converting gene symbols to Entrez IDs...")
# Convert gene symbols to Entrez IDs for pathway databases
symbols_to_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                          keys = names(ranks),
                                          column = "ENTREZID",
                                          keytype = "SYMBOL",
                                          multiVals = "first")

# Remove NAs
valid_genes <- !is.na(symbols_to_entrez)
entrez_ranks <- ranks[valid_genes]
names(entrez_ranks) <- symbols_to_entrez[valid_genes]

message(paste("Total genes in ranked list:", length(ranks)))
message(paste("Genes successfully mapped to Entrez IDs:", sum(valid_genes)))

message("Getting GO pathways...")
# Get GO Biological Process gene sets (without requiring term names)
go_bp <- as.list(org.Hs.egGO2ALLEGS)

# Filter pathways by size
go_filtered <- go_bp[sapply(go_bp, length) >= 15 & 
                    sapply(go_bp, length) <= 500]
message(paste("GO pathways after filtering:", length(go_filtered)))

message("Running GSEA with GO pathways...")
# Run GSEA with GO BP pathways
gsea_results_go <- fgsea(pathways = go_filtered, 
                       stats = entrez_ranks,
                       minSize = 15,
                       maxSize = 500)

# Just use the pathway IDs as names
gsea_results_go$pathway_name <- gsea_results_go$pathway

# Filter significant results
significant_go <- gsea_results_go %>% 
  filter(padj < 0.1) %>% 
  arrange(padj)

message("Getting Reactome pathways...")
# Get Reactome pathways
reactome_pathways <- as.list(reactomePATHID2EXTID)
reactome_names <- as.list(reactomePATHID2NAME)

# Filter pathways by size
reactome_filtered <- reactome_pathways[sapply(reactome_pathways, length) >= 15 & 
                                      sapply(reactome_pathways, length) <= 500]
message(paste("Reactome pathways after filtering:", length(reactome_filtered)))

message("Running GSEA with Reactome pathways...")
# Run GSEA with Reactome pathways
gsea_results_reactome <- fgsea(pathways = reactome_filtered, 
                             stats = entrez_ranks,
                             minSize = 15,
                             maxSize = 500)

# Add pathway names
gsea_results_reactome$pathway_name <- sapply(gsea_results_reactome$pathway, function(id) {
  if(id %in% names(reactome_names)) {
    return(reactome_names[[id]])
  } else {
    return(id)
  }
})

# Filter significant results
significant_reactome <- gsea_results_reactome %>% 
  filter(padj < 0.1) %>% 
  arrange(padj)

message("Saving results to CSV files...")

# Convert leadingEdge column which contains lists to character strings
prepare_for_csv <- function(results) {
  # Convert the leadingEdge column from list to comma-separated strings
  results$leadingEdge <- sapply(results$leadingEdge, function(x) paste(x, collapse=","))
  return(results)
}

# Prepare results for CSV export
gsea_results_go_csv <- prepare_for_csv(gsea_results_go)
gsea_results_reactome_csv <- prepare_for_csv(gsea_results_reactome)

# Save results to CSV
write.csv(gsea_results_go_csv, "gsea_results_go.csv", row.names = FALSE)
write.csv(gsea_results_reactome_csv, "gsea_results_reactome.csv", row.names = FALSE)

# Try to add GO terms to the results for better interpretability
message("Trying to add GO term descriptions (optional)...")
tryCatch({
  if (requireNamespace("AnnotationDbi", quietly = TRUE) && 
      requireNamespace("GO.db", quietly = TRUE)) {
    
    # Load GO.db and get the TERM mapping for GO IDs
    suppressPackageStartupMessages(library(GO.db))
    go_terms <- select(GO.db, keys=gsea_results_go$pathway, columns="TERM", keytype="GOID")
    
    # Create a lookup table
    term_lookup <- setNames(go_terms$TERM, go_terms$GOID)
    
    # Add terms to the results
    gsea_results_go$description <- term_lookup[gsea_results_go$pathway]
    
    # Save enhanced results
    write.csv(gsea_results_go, "gsea_results_go_with_terms.csv", row.names = FALSE)
    message("Successfully added GO term descriptions!")
  }
}, error = function(e) {
  message("Could not add GO term descriptions: ", e$message)
  message("Results still saved with GO IDs.")
})

# Print summary
cat("\n\n=== GSEA Results Summary ===\n")
cat("GO Biological Process pathways tested:", nrow(gsea_results_go), "\n")
cat("Significant GO pathways (padj < 0.1):", nrow(significant_go), "\n")
cat("Reactome pathways tested:", nrow(gsea_results_reactome), "\n")
cat("Significant Reactome pathways (padj < 0.1):", nrow(significant_reactome), "\n\n")

# Print top significant pathways if available
if(nrow(significant_go) > 0) {
  cat("Top GO pathways:\n")
  print(head(significant_go[, c("pathway_name", "NES", "pval", "padj")], 10))
} else {
  cat("No significant GO pathways found.\n")
}

cat("\n")

if(nrow(significant_reactome) > 0) {
  cat("Top Reactome pathways:\n")
  print(head(significant_reactome[, c("pathway_name", "NES", "pval", "padj")], 10))
} else {
  cat("No significant Reactome pathways found.\n")
}

cat("\nResults saved to gsea_results_go.csv and gsea_results_reactome.csv\n")
if(nrow(significant_go) > 0 || nrow(significant_reactome) > 0) {
  cat("You can use these files for downstream analysis and visualization\n")
}

message("GSEA analysis complete!")