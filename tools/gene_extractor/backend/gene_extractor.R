#!/usr/bin/env Rscript

# ============================================================
# gene Marker and Gene Extractor — headless backend
#
# This file contains the SAME extract_gene_data() logic that
# was inside the original Shiny app.R (server-side function).
# It has not been changed in any way.
#
# What changed is only *how it is invoked*: instead of Shiny's
# fileInput/actionButton/downloadHandler reactivity, this script
# is called from the command line (by the Python/Qt backend) with
# plain file paths, and it prints a single JSON object to stdout
# describing the result. That JSON is what the PySide6 bridge
# reads and forwards to the HTML frontend.
# ============================================================

suppressWarnings(suppressMessages({
  library(readxl)
  library(dplyr)
  library(writexl)
  library(jsonlite)
}))

# ------------------------------------------------------------
# ORIGINAL FUNCTION — unchanged from app.R
# ------------------------------------------------------------

extract_gene_data <- function(extracted_markers_path, gene_markers_path, rice_gene_path, output_dir) {
  # Extract markers with physical positions
  extracted_markers_df <- read_excel(extracted_markers_path, 
                                     col_types = c("text", "text", "text"), 
                                     col_names = c("Markers", "Chromosome", "Position"))
  
  gene_markers_df <- read_excel(gene_markers_path, 
                               col_types = c("text", "text"), 
                               col_names = c("Markers", "Chromosome"))
  
  merged_df <- left_join(gene_markers_df, extracted_markers_df, by = c("Markers", "Chromosome"))
  markers_final <- merged_df %>% select(Markers, Chromosome, Position)
  
  markers_output_path <- file.path(output_dir, 'Markers_with_physical_position.xlsx')
  write_xlsx(markers_final, markers_output_path)
  
  # Extract genes within gene regions
  rice_gene_list <- read_excel(rice_gene_path, skip = 1, 
                               col_names = c('Chromosome', 'Start', 'Stop', 'Strand', 'Gene_id', 'Description'))
  
  # Use the markers with positions to define gene regions
  gene_markers_list <- markers_final %>% 
    filter(!is.na(Position)) %>%
    group_by(Chromosome) %>%
    summarise(start = min(as.numeric(Position)), stop = max(as.numeric(Position))) %>%
    ungroup()
  
  # Convert to numeric
  rice_gene_list$Start <- as.numeric(rice_gene_list$Start)
  rice_gene_list$Stop <- as.numeric(rice_gene_list$Stop)
  
  # Initialize results
  results <- data.frame(Chromosome = character(), Start = numeric(), Stop = numeric(), 
                        Strand = character(), Gene_id = character(), Description = character(), 
                        stringsAsFactors = FALSE)
  
  # Loop through each gene marker region
  if(nrow(gene_markers_list) > 0) {
    for (i in 1:nrow(gene_markers_list)) {
      chromosome <- gene_markers_list$Chromosome[i]
      gene_start_position <- gene_markers_list$start[i]
      gene_stop_position <- gene_markers_list$stop[i]
      
      # Filter genes within the gene region
      filtered_genes <- rice_gene_list %>%
        filter(Chromosome == chromosome &
                 Start >= gene_start_position &
                 Stop <= gene_stop_position)
      
      # Add gene position information
      if(nrow(filtered_genes) > 0) {
        filtered_genes <- filtered_genes %>%
          mutate(gene_start = gene_start_position, gene_stop = gene_stop_position)
        
        results <- rbind(results, filtered_genes)
      }
    }
  }
  
  genes_output_path <- file.path(output_dir, 'Genes_in_gene_regions.xlsx')
  write_xlsx(results, genes_output_path)
  
  return(list(markers = markers_final, genes = results))
}

# ------------------------------------------------------------
# CLI ENTRYPOINT (new — replaces the Shiny UI/server wiring)
#
# Usage:
#   Rscript gene_extractor.R <extracted_markers_path> <gene_markers_path> <rice_gene_path> <output_dir>
#
# Prints one JSON object to stdout:
#   { "status": "success", "markers": [...], "genes": [...],
#     "markers_file": "...", "genes_file": "..." }
#   or
#   { "status": "error", "message": "..." }
# ------------------------------------------------------------

emit_error <- function(message) {
  cat(toJSON(list(status = "error", message = message), auto_unbox = TRUE, na = "null"))
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  emit_error("Usage: Rscript gene_extractor.R <extracted_markers_path> <gene_markers_path> <rice_gene_path> <output_dir>")
  quit(status = 1, save = "no")
}

extracted_markers_path <- args[1]
gene_markers_path       <- args[2]
rice_gene_path         <- args[3]
output_dir             <- args[4]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

result <- tryCatch({

  # Same validation the Shiny server did before calling extract_gene_data()
  if (!file.exists(extracted_markers_path) ||
      !file.exists(gene_markers_path) ||
      !file.exists(rice_gene_path)) {
    stop("Please upload all three required files.")
  }

  res <- extract_gene_data(
    extracted_markers_path,
    gene_markers_path,
    rice_gene_path,
    output_dir
  )

  list(
    status = "success",
    markers = res$markers,
    genes = res$genes,
    markers_file = file.path(output_dir, "Markers_with_physical_position.xlsx"),
    genes_file = file.path(output_dir, "Genes_in_gene_regions.xlsx")
  )

}, error = function(e) {
  list(
    status = "error",
    message = paste("Error processing files:", conditionMessage(e))
  )
})

cat(toJSON(result, auto_unbox = TRUE, na = "null"))
