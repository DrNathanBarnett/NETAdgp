\name{NETA_dgpAffy}
\alias{NETA_dgpAffy}

\title{
Automated Differential Gene Expression Analysis and Pathway Mapping for Affymetrix Microarray Platforms utilizing AffyBatch Objects from CEL Files
}
\description{
This function automates differential gene expression analysis, Gene Ontology pathway enrichment, and protein-protein interaction mapping for Affymetrix microarray platforms. It processes raw AffyBatch objects generated from CEL files using the `affy` package and performs downstream analysis.
}
\usage{
NETA_dgpAffy(affy_batch, phenotype, p_value_cutoff = 0.05)
}
\arguments{
  \item{affy_batch}{
An \code{AffyBatch} object created from Affymetrix CEL files. The object contains raw microarray data, which will be processed for normalization and analysis.
}
  \item{phenotype}{
A vector or factor describing the experimental conditions or groups that correspond to the samples in the \code{AffyBatch} object.
}
  \item{p_value_cutoff}{
A numeric value specifying the significance threshold for differential gene expression analysis. The default value is 0.05.
}
}
\details{
This function handles the preprocessing (including background correction and normalization) of AffyBatch objects from Affymetrix CEL files, followed by differential gene expression analysis. It also performs Gene Ontology (GO) enrichment analysis and protein-protein interaction (PPI) mapping, providing a complete workflow for Affymetrix microarray data analysis.
}
\value{
  \item{DE_results}{A data frame of differentially expressed genes, including fold changes, p-values, and adjusted p-values.}
  \item{GO_enrichment}{A list of enriched Gene Ontology terms with corresponding p-values and enrichment scores.}
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
\code{\link{affy}}, \code{\link{limma}}, \code{\link{GOstats}}, \code{\link{topTable}}
}
\examples{
## Example usage:
library(affy)
cel_files <- list.celfiles("path/to/cel_files", full.names = TRUE)
affy_data <- ReadAffy(filenames = cel_files)

phenotype <- factor(c("control", "treated", "control", "treated"))

results <- NETA_dgpAffy(affy_data, phenotype)
top_genes <- results$DE_results
GO_results <- results$GO_enrichment
PPI_network <- results$PPI_network
}
\keyword{Affymetrix}
\keyword{DifferentialExpression}
\keyword{PathwayAnalysis}
\concept{GeneExpression}
\concept{Microarray}
