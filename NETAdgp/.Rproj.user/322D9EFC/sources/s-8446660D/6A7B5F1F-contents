\name{NETA_dgpIlluminaTxt}
\alias{NETA_dgpIlluminaTxt}
\title{
Automated Differential Gene Expression Analysis and Pathway Mapping for Illumina Microarray Platforms without .idat Files
}
\description{
This function automates differential gene expression analysis, Gene Ontology pathway enrichment, and protein-protein interaction mapping for Illumina microarray platforms where raw data in .idat files are not available. The input data should be in a matrix format, with Gene Symbols as row names and samples as column names.
}
\usage{
NETA_dgpIlluminaTxt(illumina_data_file, phenotype, p_value_cutoff = 0.05)
}
\arguments{
  \item{illumina_data_file}{
A matrix or data frame containing the Illumina microarray data, where rows represent genes and columns represent samples. The row names should correspond to Gene Symbols.
}
  \item{phenotype}{
A vector or factor representing the experimental conditions or groupings of the samples. This should match the order of the columns in the input data.
}
  \item{p_value_cutoff}{
A numeric value indicating the significance threshold for differential gene expression analysis. The default value is 0.05.
}
}
\details{
This function is designed for scenarios where Illumina microarray data is provided in text format (e.g., CSV or TSV) and the original .idat files are not available. It performs differential gene expression analysis followed by Gene Ontology (GO) enrichment and protein-protein interaction (PPI) mapping. The function is compatible with various Illumina microarray platforms.
}
\value{
  \item{DE_results}{A data frame containing the differentially expressed genes, including fold changes, p-values, and adjusted p-values.}
  \item{GO_enrichment}{A list of enriched Gene Ontology terms, including p-values and enrichment scores.}
  \item{PPI_network}{A data frame or graph object representing the protein-protein interaction network for the significant genes.}
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
\code{\link{limma}}, \code{\link{lumi}}, \code{\link{topTable}}, \code{\link{GOstats}}
}
\examples{
## Example usage:
illumina_data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(illumina_data) <- paste0("Gene", 1:100)
phenotype <- factor(c(rep("control", 5), rep("treated", 5)))

results <- NETA_dgpIlluminaTxt(illumina_data, phenotype)
top_genes <- results$DE_results
GO_results <- results$GO_enrichment
PPI_network <- results$PPI_network
}
\keyword{Illumina}
\keyword{DifferentialExpression}
\keyword{PathwayAnalysis}
\concept{GeneExpression}
\concept{Microarray}
