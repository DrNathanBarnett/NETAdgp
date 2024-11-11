\name{NETA_dgpscRNAseq}
\alias{NETA_dgpscRNAseq}

\title{
Automated Differential Gene Expression Analysis and Pathway Mapping for Single-Cell and Single-Nucleus RNA Sequencing
}
\description{
This function automates differential gene expression analysis and pathway mapping for single-cell and single-nucleus RNA sequencing data. The pathway enrichment analysis is performed through Gene Ontology mapping, while protein-protein interactions are identified using the STRING.db connection. Additionally, cluster marker genes are mapped to specific cell types.
}
\usage{
NETA_dgpscRNAseq(seurat_object, p_value_cutoff = 0.05)
}

\arguments{
  \item{seurat_object}{
An object of class \code{Seurat} containing single-cell or single-nucleus RNA sequencing data, including expression matrices, metadata, and any relevant clustering information.
}
  \item{p_value_cutoff}{
A numeric value specifying the significance threshold for differential gene expression analysis. The default value is 0.05.
}
}
\details{
This function handles the differential gene expression analysis of single-cell RNA sequencing data using Seurat objects. It also performs Gene Ontology (GO) enrichment analysis and protein-protein interaction mapping through STRING.db, providing a comprehensive workflow for analyzing single-cell data.
}
\value{
  \item{cluster_markers}{A data frame of differentially expressed genes, including fold changes, p-values, and adjusted p-values.}
  \item{go_results}{A list of enriched Gene Ontology terms with corresponding p-values and enrichment scores.}
  \item{ppi_results}{A data frame or graph object representing the protein-protein interaction network for significant genes.}
  \item{cluster_annotations}{A data frame of mapped cell annotations with the number of genes according to the cell type and the organ of origin.}
  \item{marker_genes}{A list of identified marker genes mapped to their corresponding cell types.}
}
\references{
Barnett, N.J. et al., Neural Epileptotranscriptomic Advancement (NETA) Discovery, 2024.
% Add more references as necessary.
}
\author{
Dr Nathan James Barnett, PhD \\ 
Chief Executive Officer & Managing Director \\ 
Neural Epileptotranscriptomic Advancement (NETA) Discovery \\ 
Prof Craig S McLachlan \\ 
Chairman \\ 
Neural Epileptotranscriptomic Advancement (NETA) Discovery
}
\note{
This function is part of the NETA Discovery platform and is distributed under its license. For more information, visit the NETA Discovery website.
}

\seealso{
\code{\link{Seurat}}, \code{\link{STRINGdb}}, \code{\link{GOstats}}, \code{\link{topTable}}
% Add more relevant links as necessary.
}
\examples{
## Example usage:
library(Seurat)
# Assuming 'seurat_data' is a pre-processed Seurat object
results <- NETA_dgpscRNAseq(seurat_data, p_value_cutoff = 0.05)
top_markers <- results$top_markers
GO_results <- results$GO_enrichment
ppi_results <- results$ppi_results
cluster_markers <- results$cluster_markers
}
\keyword{SingleCellRNAseq}
\keyword{DifferentialExpression}
\keyword{PathwayAnalysis}
\concept{GeneExpression}
\concept{SingleCell}
