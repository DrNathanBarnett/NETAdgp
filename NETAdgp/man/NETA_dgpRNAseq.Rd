\name{NETA_dgpRNAseq}
\alias{NETA_dgpRNAseq}
\title{
Automated RNA-seq Differential Expression and Pathway Analysis Tool
}
\description{
Automates differential gene expression analysis, Gene Ontology pathway enrichment, and protein-protein interaction mapping for raw bulk RNA-seq count data. The data should have Gene Symbols as row names.
}
\usage{
NETA_dgpRNAseq(count_data, phenotype, p_value_cutoff = 0.05)
}
\arguments{
  \item{count_data}{A matrix of raw RNA-seq count data with genes as row names and samples as column names.}
  \item{phenotype}{A vector or factor describing the experimental conditions or sample groups.}
  \item{p_value_cutoff}{A numeric value for the significance threshold in differential expression analysis. Default is 0.05.}
}
\details{
This function handles RNA-seq differential expression analysis, followed by gene ontology enrichment analysis and protein-protein interaction mapping using publicly available databases.
}
\value{
  \item{DE_results}{A data frame containing the differentially expressed genes, their fold changes, and p-values.}
  \item{GO_enrichment}{A list of enriched Gene Ontology terms and their associated statistics.}
  \item{PPI_network}{A data frame or graph object representing protein-protein interaction networks for significant genes.}
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
Owned and distributed by Neural Epileptotranscriptomic Advancement Discovery.
}
\examples{
## Example of usage:
count_data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(count_data) <- paste0("Gene", 1:100)
phenotype <- factor(c(rep("control", 5), rep("treated", 5)))

results <- NETA_dgpRNAseq(count_data, phenotype)
}
\keyword{RNA-seq}
\keyword{Differential Expression}
\keyword{Pathway Enrichment}
