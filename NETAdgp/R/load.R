# NEURAL EPILEPTOTRANSCRIPTOMIC ADVANCEMENT Proprietary Software
# © 2024 NEURAL EPILEPTOTRANSCRIPTOMIC ADVANCEMENT. All rights reserved.
#
# This software is licensed under the NEURAL EPILEPTOTRANSCRIPTOMIC ADVANCEMENT Software License Agreement.
# You may not use, modify, or distribute this software except in compliance with the license.
# Please refer to the LICENSE file for the full terms of the license.


NETA_dgpRNAseq <- function(count_data, phenotype, p_value_cutoff = 0.05) {
  
  if (!is.factor(phenotype)) {
    phenotype <- factor(phenotype)
  }
  
  if (length(phenotype) != ncol(count_data)) {
    stop("The length of 'phenotype' must match the number of columns in 'count_data'.")
  }
  
  count_dataDGE <- DGEList(counts = count_data, group = phenotype)
  
  count_dataDGE <- calcNormFactors(count_dataDGE)
  
  count_dataDGE <- estimateDisp(count_dataDGE)
  fit <- glmFit(count_dataDGE)
  lrt <- glmLRT(fit)
  
  Top_table <- topTags(lrt, n = Inf)
  degs <- Top_table$table
  
  significant_degs <- degs[degs$PValue < p_value_cutoff, ]
  gene_list <- rownames(significant_degs)
  map_genes_to_go <- function(gene_list) {
    go_results <- gost(query = gene_list, 
                       organism = "hsapiens", 
                       sources = c("GO:BP", "GO:CC", "GO:MF"))
    
    go_terms <- unique(go_results$result$term_name)
    go_matrix <- sapply(go_terms, function(go_term) as.integer(go_results$result$term_name == go_term))
    rownames(go_matrix) <- go_results$result$query
    return(list(go_mappings = go_results$result, go_matrix = go_matrix))
  }
  go_results <- map_genes_to_go(gene_list)
  find_ppi <- function(gene_list) {
    string_db <- STRINGdb$new(version = "12", species = 9606, score_threshold = 400, input_directory = "", protocol = "http",  network_type="full", link_data='combined_only')
    mapped_genes <- string_db$map(data.frame(gene_list), "gene_list", removeUnmappedRows = TRUE)
    ppi_network <- string_db$get_interactions(mapped_genes$STRING_id)
    return(ppi_network)
  }
  ppi_results <- find_ppi(gene_list)
  return(list(degs = significant_degs, go_results = go_results, ppi_results = ppi_results))
}




NETA_dgpOligo <- function(oligo_data, phenotype, p_value_cutoff = 0.05) {
  
  if (!is.factor(phenotype)) {
    phenotype <- factor(phenotype)
  }
  
  if (length(phenotype) != ncol(exprs(oligo_data))) {
    stop("The length of 'phenotype' must match the number of samples in 'oligo_data'.")
  }
  
  norm_data <- rma(oligo_data)  
  exprs_data <- exprs(norm_data)  
  
  design <- model.matrix(~ phenotype)
  fit <- lmFit(exprs_data, design)     
  fit <- eBayes(fit)                   
  top_table <- topTable(fit, number = Inf, adjust.method = "BH")
  significant_degs <- top_table[top_table$adj.P.Val < p_value_cutoff, ]
  
  gene_list <- rownames(significant_degs)
  
  map_genes_to_go <- function(gene_list) {
    go_results <- gost(query = gene_list, 
                       organism = "hsapiens", 
                       sources = c("GO:BP", "GO:CC", "GO:MF"))
    
    go_terms <- unique(go_results$result$term_name)
    go_matrix <- sapply(go_terms, function(go_term) as.integer(go_results$result$term_name == go_term))
    
    rownames(go_matrix) <- go_results$result$query
    
    return(list(go_mappings = go_results$result, go_matrix = go_matrix))
  }
  go_results <- map_genes_to_go(gene_list)
  
  find_ppi <- function(gene_list) {
    string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "", protocol = "http", network_type="full", link_data='combined_only')
    
    mapped_genes <- string_db$map(data.frame(gene_list), "gene_list", removeUnmappedRows = TRUE)
    
    ppi_network <- string_db$get_interactions(mapped_genes$STRING_id)
    
    return(ppi_network)
  }
  
  ppi_results <- find_ppi(gene_list)
  return(list(degs = significant_degs, go_results = go_results, ppi_results = ppi_results))
}








NETA_dgpIlluminaTxt <- function(illumina_data_file, phenotype, p_value_cutoff = 0.05) {
  
  if (!is.factor(phenotype)) {
    phenotype <- factor(phenotype)
  }
  
  raw_data <- lumiR(illumina_data_file) 
  
  norm_data <- lumiExpresso(raw_data, bg.correct = TRUE, variance.stabilize = TRUE, normalize.param = list(method = "rsn"))
  exprs_data <- exprs(norm_data)  
  
  design <- model.matrix(~ phenotype)  
  fit <- lmFit(exprs_data, design)     
  fit <- eBayes(fit)                   
  
  top_table <- topTable(fit, number = Inf, adjust.method = "BH")
  significant_degs <- top_table[top_table$adj.P.Val < p_value_cutoff, ]
  
  gene_list <- rownames(significant_degs)
  
  map_genes_to_go <- function(gene_list) {
    go_results <- gost(query = gene_list, 
                       organism = "hsapiens", 
                       sources = c("GO:BP", "GO:CC", "GO:MF"))
    
    go_terms <- unique(go_results$result$term_name)
    go_matrix <- sapply(go_terms, function(go_term) as.integer(go_results$result$term_name == go_term))
    
    rownames(go_matrix) <- go_results$result$query
    
    return(list(go_mappings = go_results$result, go_matrix = go_matrix))
  }
  
  go_results <- map_genes_to_go(gene_list)
  
  find_ppi <- function(gene_list) {
    string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "", protocol = "http", network_type="full", link_data='combined_only')
    
    mapped_genes <- string_db$map(data.frame(gene_list), "gene_list", removeUnmappedRows = TRUE)
    
    ppi_network <- string_db$get_interactions(mapped_genes$STRING_id)
    
    return(ppi_network)
  }
  
  ppi_results <- find_ppi(gene_list)
  
  return(list(degs = significant_degs, go_results = go_results, ppi_results = ppi_results))
}



NETA_dgpIlluminaIDAT <- function(idat_file_paths, phenotype, p_value_cutoff = 0.05) {
  
  if (!is.factor(phenotype)) {
    phenotype <- factor(phenotype)
  }
  
  raw_data_list <- lapply(idat_file_paths, readIDAT)  # Read all .idat files
  
  raw_data_matrix <- do.call(cbind, lapply(raw_data_list, function(x) x$Quants[, "Mean"]))
  
  lumi_data <- lumiR.idat(idat_file_paths)  # Read the .idat files into a LumiBatch object
  norm_data <- lumiExpresso(lumi_data, bg.correct = TRUE, variance.stabilize = TRUE, normalize.param = list(method = "quantile"))
  
  exprs_data <- exprs(norm_data)
  
  design <- model.matrix(~ phenotype)  
  fit <- lmFit(exprs_data, design)     
  fit <- eBayes(fit)                   
  top_table <- topTable(fit, number = Inf, adjust.method = "BH")
  significant_degs <- top_table[top_table$adj.P.Val < p_value_cutoff, ]
  
  gene_list <- rownames(significant_degs)
  
  map_genes_to_go <- function(gene_list) {
    go_results <- gost(query = gene_list, 
                       organism = "hsapiens", 
                       sources = c("GO:BP", "GO:CC", "GO:MF"))
    
    go_terms <- unique(go_results$result$term_name)
    go_matrix <- sapply(go_terms, function(go_term) as.integer(go_results$result$term_name == go_term))
    
    rownames(go_matrix) <- go_results$result$query
    
    return(list(go_mappings = go_results$result, go_matrix = go_matrix))
  }
  
  go_results <- map_genes_to_go(gene_list)
  
  find_ppi <- function(gene_list) {
    string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "", protocol = "http", network_type="full", link_data='combined_only')
    
    mapped_genes <- string_db$map(data.frame(gene_list), "gene_list", removeUnmappedRows = TRUE)
    
    ppi_network <- string_db$get_interactions(mapped_genes$STRING_id)
    
    return(ppi_network)
  }
  
  # Step 10: Call the PPI function to retrieve protein-protein interactions
  ppi_results <- find_ppi(gene_list)
  
  # Step 11: Return the significant DEGs, GO results, and PPI network
  return(list(degs = significant_degs, go_results = go_results, ppi_results = ppi_results))
}









NETA_dgpAgilent <- function(targets_file, phenotype, p_value_cutoff = 0.05) {
  
  if (!is.factor(phenotype)) {
    phenotype <- factor(phenotype)
  }
  
  targets <- readTargets(targets_file)
  
  raw_data <- read.maimages(targets$FileName, source = "agilent", green.only = TRUE)
  
  bg_corrected_data <- backgroundCorrect(raw_data, method = "normexp", offset = 50)
  
  norm_data <- normalizeBetweenArrays(bg_corrected_data, method = "quantile")
  
  exprs_data <- norm_data$E
  
  design <- model.matrix(~ phenotype)  
  fit <- lmFit(exprs_data, design)     
  fit <- eBayes(fit)                   
  
  top_table <- topTable(fit, number = Inf, adjust.method = "BH")
  significant_degs <- top_table[top_table$adj.P.Val < p_value_cutoff, ]
  
  gene_list <- rownames(significant_degs)
  
  map_genes_to_go <- function(gene_list) {
    go_results <- gost(query = gene_list, 
                       organism = "hsapiens", 
                       sources = c("GO:BP", "GO:CC", "GO:MF"))
    
    go_terms <- unique(go_results$result$term_name)
    go_matrix <- sapply(go_terms, function(go_term) as.integer(go_results$result$term_name == go_term))
    
    rownames(go_matrix) <- go_results$result$query
    
    return(list(go_mappings = go_results$result, go_matrix = go_matrix))
  }
  
  go_results <- map_genes_to_go(gene_list)
  
  find_ppi <- function(gene_list) {
    string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "", protocol = "http", network_type="full", link_data='combined_only')
    
    mapped_genes <- string_db$map(data.frame(gene_list), "gene_list", removeUnmappedRows = TRUE)
    
    ppi_network <- string_db$get_interactions(mapped_genes$STRING_id)
    
    return(ppi_network)
  }
  
  ppi_results <- find_ppi(gene_list)
  
  return(list(degs = significant_degs, go_results = go_results, ppi_results = ppi_results))
}



NETA_dgpAffy <- function(affy_batch, phenotype, p_value_cutoff = 0.05) {
  
  if (!is.factor(phenotype)) {
    phenotype <- factor(phenotype)
  }
  
  if (length(phenotype) != length(affy_batch)) {
    stop("The length of 'phenotype' must match the number of samples in 'affy_batch'.")
  }
  
  norm_data <- mas5(affy_batch)
  exprs_data <- exprs(norm_data)  # Extract normalized expression values
  
  library(limma)  
  
  design <- model.matrix(~ phenotype)  
  fit <- lmFit(exprs_data, design)    
  fit <- eBayes(fit)                   
  
  top_table <- topTable(fit, number = Inf, adjust.method = "BH")
  significant_degs <- top_table[top_table$adj.P.Val < p_value_cutoff, ]
  
  gene_list <- rownames(significant_degs)
  
  map_genes_to_go <- function(gene_list) {
    go_results <- gost(query = gene_list, 
                       organism = "hsapiens", 
                       sources = c("GO:BP", "GO:CC", "GO:MF"))
    
    go_terms <- unique(go_results$result$term_name)
    go_matrix <- sapply(go_terms, function(go_term) as.integer(go_results$result$term_name == go_term))
    
    rownames(go_matrix) <- go_results$result$query
    
    return(list(go_mappings = go_results$result, go_matrix = go_matrix))
  }
  
  go_results <- map_genes_to_go(gene_list)
  
  find_ppi <- function(gene_list) {
    string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "", protocol = "http", network_type="full", link_data='combined_only')
    
    mapped_genes <- string_db$map(data.frame(gene_list), "gene_list", removeUnmappedRows = TRUE)
    
    ppi_network <- string_db$get_interactions(mapped_genes$STRING_id)
    
    return(ppi_network)
  }
  
  ppi_results <- find_ppi(gene_list)
  
  return(list(degs = significant_degs, go_results = go_results, ppi_results = ppi_results))
}





NETA_dgpscRNAseq <- function(seurat_object, p_value_cutoff = 0.05) {
  
  # Step 1: Quality Control on Seurat Object
  seurat_object$sample <- rownames(seurat_object@meta.data)
  
  # Split sample column
  seurat_object@meta.data <- separate(seurat_object@meta.data, col = 'sample', into = c('Patient', 'Barcode', 'Tissue'), sep = '_')
  
  # Calculate mitochondrial percentage
  seurat_object$mitoPercent <- PercentageFeatureSet(seurat_object, pattern='^MT-')
  
  # Filter cells based on quality metrics
  seurat_filtered <- subset(seurat_object, subset = nCount_RNA > 800 &
                              nFeature_RNA > 500 &
                              mitoPercent < 10)
  
  # Step 2: Standard Workflow in Seurat for Batch Effect Detection
  seurat_filtered <- NormalizeData(seurat_filtered)
  seurat_filtered <- FindVariableFeatures(seurat_filtered)
  seurat_filtered <- ScaleData(seurat_filtered)
  seurat_filtered <- RunPCA(seurat_filtered)
  ElbowPlot(seurat_filtered)
  seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:20)
  seurat_filtered <- FindClusters(seurat_filtered)
  seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:20)
  
  # Step 3: Find Differentially Expressed Genes (DEGs) per Cluster
  cluster_markers <- FindAllMarkers(seurat_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 3)
  top_markers <- cluster_markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  
  # Step 4: Annotate Clusters with Clustermole
  cluster_annotations <- clustermole_overlaps(cluster_markers$gene, species = 'hs')  # Automated cell type annotation using clustermole
  
  # Step 5: Extract Gene List for GO and PPI Mapping
  gene_list <- unique(cluster_markers$gene)
  
  # Step 6: Map DEGs to GO Terms and Create a GO Term Binary Matrix
  map_genes_to_go <- function(gene_list) {
    # Perform GO enrichment analysis using gProfiler2
    go_results <- gost(query = gene_list, 
                       organism = "hsapiens", 
                       sources = c("GO:BP", "GO:CC", "GO:MF"))
    
    # Create a binary matrix for GO terms
    go_terms <- unique(go_results$result$term_name)
    go_matrix <- sapply(go_terms, function(go_term) as.integer(go_results$result$term_name == go_term))
    
    # Add gene names as rownames
    rownames(go_matrix) <- go_results$result$query
    
    # Return both GO mappings and the binary matrix
    return(list(go_mappings = go_results$result, go_matrix = go_matrix))
  }
  
  # Step 7: Call the GO Term Mapping Function
  go_results <- map_genes_to_go(gene_list)
  
  # Step 8: Find Protein-Protein Interactions (PPIs) Using STRINGdb
  find_ppi <- function(gene_list) {
    # Initialize the STRINGdb object for human (taxID = 9606)
    
    library(httr)
    
    # URL for the file
    url <- "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
    
    response <- GET(url, write_disk("9606.protein.aliases.v12.0.txt.gz", overwrite = TRUE), config(ssl_verifypeer = FALSE))
    
    # Check if the download was successful
    if (http_status(response)$category != "Success") {
      print("Download failed!")
    } else {
      print("Download successful!")
    }
    string_db <- STRINGdb$new(version = "12", species = 9606, score_threshold = 400, input_directory = "",  protocol="http", network_type="full", link_data='combined_only')
    
    # Map gene names to STRING IDs
    mapped_genes <- string_db$map(data.frame(gene_list), "gene_list", removeUnmappedRows = TRUE)
    
    
    # Retrieve PPI network
    ppi_network <- string_db$get_interactions(mapped_genes$STRING_id)
    
    return(ppi_network)
  }
  
  # Step 9: Call the PPI Function to Retrieve Protein-Protein Interactions
  ppi_results <- find_ppi(gene_list)
  
  # Step 10: Return the Marker Genes, Cluster Annotations, GO Results, and PPI Network
  return(list(cluster_markers = cluster_markers, 
              top_markers = top_markers, 
              cluster_annotations = cluster_annotations, 
              go_results = go_results, 
              ppi_results = ppi_results))
}

