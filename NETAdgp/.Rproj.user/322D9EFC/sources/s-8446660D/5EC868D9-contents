\name{NETA_dgpIlluminaIDAT}
\alias{NETA_dgpIlluminaIDAT}
\title{
Automated Differential Gene Expression Analysis and Pathway Mapping for Illumina Microarray Platforms with .idat Files
}
\description{
This function automates differential gene expression analysis, Gene Ontology pathway enrichment, and protein-protein interaction mapping for Illumina microarray platforms where raw data in .idat files is available. The input data should be processed and converted from the raw .idat format, with Gene Symbols as row names and samples as column names.
}
\usage{
function(idat_file_paths, phenotype, p_value_cutoff = 0.05) 
}
\arguments{
  \item{idat_file_paths}{
A character vector containing the file paths to the Illumina .idat files.
}
  \item{phenotype}{
A vector or factor representing the experimental conditions or groupings of the samples, matching the order in the `sample_info`.
}
  \item{p_value_cutoff}{
A numeric value specifying the significance threshold for differential gene expression analysis. The default value is 0.05.
}
}
\details{
This function is specifically designed for Illumina microarray platforms where raw .idat files are available. It handles the preprocessing of .idat files, including background correction and normalization, followed by differential gene expression analysis, Gene Ontology (GO) enrichment analysis, and protein-protein interaction (PPI) network mapping. It is compatible with multiple Illumina array types.
}
\value{
  \item{DE_results}{A data frame of differentially expressed genes, including fold changes, p-values, and adjusted p-values.}
  \item{GO_enrichment}{A list of enriched Gene Ontology terms with corresponding p-values and enrichment scores.}
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
\code{\link{limma}}, \code{\link{lumi}}, \code{\link{illuminaio}}, \code{\link{GOstats}}, \code{\link{oligo}}
}
\examples{
## Example usage:
idat_files <- c("path/to/sample1.idat", "path/to/sample2.idat")
sample_info <- data.frame(SampleID = c("Sample1", "Sample2"), Group = c("control", "treated"))
phenotype <- factor(c("control", "treated"))

results <- NETA_dgpIlluminaIDAT(idat_files, sample_info, phenotype)
top_genes <- results$DE_results
GO_results <- results$GO_enrichment
PPI_network <- results$PPI_network
}
\keyword{Illumina}
\keyword{DifferentialExpression}
\keyword{PathwayAnalysis}
\concept{GeneExpression}
\concept{Microarray}
