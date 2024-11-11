\name{NETA_dgpOligo}
\alias{NETA_dgpOligo}
\title{
Automated Differential Gene Expression Analysis and Pathway Mapping for NimbleGen Oligo Arrays
}
\description{
Automates differential gene expression analysis, Gene Ontology pathway enrichment, and protein-protein interaction mapping for NimbleGen-based oligonucleotide arrays. These high-density arrays can also support comparative genomic hybridization (CGH) and DNA methylation studies. The input data should have Gene Symbols as row names.
}
\usage{
NETA_dgpOligo(oligo_data, phenotype, p_value_cutoff = 0.05)
}
\arguments{
  \item{oligo_data}{
A matrix containing the raw oligonucleotide array data with genes as row names and samples as column names. 
}
  \item{phenotype}{
A vector or factor describing the experimental conditions or sample groupings corresponding to the columns in the data matrix.
}
  \item{p_value_cutoff}{
A numeric value specifying the significance threshold for differential gene expression. Default is 0.05.
}
}
\details{
This function performs preprocessing, normalization, and differential gene expression analysis for NimbleGen oligonucleotide arrays. It also includes Gene Ontology enrichment and protein-protein interaction (PPI) mapping for further functional analysis of significant genes.
}
\value{
  \item{DE_results}{A data frame of differentially expressed genes with corresponding fold changes and p-values.}
  \item{GO_enrichment}{A list of enriched Gene Ontology terms with associated statistics.}
  \item{PPI_network}{A data frame or graph object representing the protein-protein interaction network for significant genes.}
}
\references{
Barnett, N.J. et al., Neural Epileptotranscriptomic Advancement (NETA) Discovery, 2024.
}
\author{
Dr Nathan James Barnett PhD 

Chief Executive Officer & Managing Director

Neural Epileptotranscriptomic Advancement (NETA) Discovery

Prof Craig S McLachlan

Chairman

Neural Epileptotranscriptomic Advancement (NETA) Discovery
}
\note{
This function is part of the NETA Discovery platform and is distributed under its license. For more information, visit the NETA Discovery website.
}
\seealso{
\code{\link{limma}}, \code{\link{oligo}}, \code{\link{topTable}}, \code{\link{GOstats}}
}
\examples{
## Example of usage:
oligo_data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(oligo_data) <- paste0("Gene", 1:100)
phenotype <- factor(c(rep("control", 5), rep("treated", 5)))

results <- NETA_dgpOligo(oligo_data, phenotype)
top_genes <- results$DE_results
GO_results <- results$GO_enrichment
PPI_network <- results$PPI_network
}
\keyword{OligoArray}
\keyword{DifferentialExpression}
\keyword{PathwayAnalysis}
\concept{GeneExpression}
\concept{NimbleGen}
