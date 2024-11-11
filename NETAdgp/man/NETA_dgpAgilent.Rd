\name{NETA_dgpAgilent}
\alias{NETA_dgpAgilent}
%- Add aliases for related topics if necessary
\title{
Automated Differential Gene Expression Analysis and Pathway Mapping for Agilent Microarray Platforms
}
\description{
This function automates the workflow for differential gene expression analysis, Gene Ontology (GO) pathway enrichment, and protein-protein interaction (PPI) mapping for Agilent microarray data. The input data is expected to be in a tab-delimited format with gene symbols as row names and samples as column names.
}
\usage{
NETA_dgpAgilent(targets_file, phenotype, p_value_cutoff = 0.05)
}
\arguments{
  \item{targets_file}{
    The raw Agilent microarray tab-delimited file. The file should contain gene symbols as row names and sample expression values as columns.
  }
  \item{phenotype}{
    A factor vector indicating the phenotype or condition for each sample (e.g., control and treatment groups). The length of the phenotype must match the number of columns (samples) in the data.
  }
  \item{p_value_cutoff}{
    A numeric value specifying the p-value cutoff for identifying differentially expressed genes (DEGs). Default is 0.05.
  }
}
\details{
The function reads raw Agilent microarray data from a tab-delimited file, performs background correction and quantile normalization using the `neqc` function from the \code{limma} package, and fits a linear model to identify DEGs. GO term enrichment is performed using \code{gProfiler2}, and protein-protein interactions are mapped using \code{STRINGdb}. 
}
\value{
A list containing the following components:
  \item{degs}{A data frame of differentially expressed genes (DEGs), including gene symbols and relevant statistics such as log fold changes and adjusted p-values.}
  \item{go_results}{A list containing Gene Ontology (GO) mappings, including biological processes, molecular functions, and cellular components.}
  \item{ppi_results}{A data frame of protein-protein interaction (PPI) networks derived from \code{STRINGdb}.}
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
This function assumes that the input file is properly formatted and contains raw gene expression data from Agilent microarrays.
}
\seealso{
  \code{\link{limma}}, \code{\link{gprofiler2}}, \code{\link{STRINGdb}}, \code{\link{clustermole}}
}
\examples{
# Example usage:
# Read Agilent microarray data and perform differential expression analysis
phenotype <- factor(c("control", "treatment", "control", "treatment"))
results <- NETA_dgpAgilent("path_to_agilent_data.txt", phenotype)

# View results
head(results$degs)
head(results$go_results)
head(results$ppi_results)
}
\keyword{microarray}
\keyword{differential expression}
\keyword{GO enrichment}
\keyword{protein-protein interaction}
