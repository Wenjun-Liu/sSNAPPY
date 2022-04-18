#' logCPM_example: Normalised logCPM of patient-derived explant models obtained from 5 ER-positive primamry breast cancer patients (GSE80098)
#'
#' This data was adopted from a study by Singhal H, et al., which was published as *Genomic agonism and phenotypic antagonism between estrogen and progesterone receptors in breast cancer*
#' in 2016. In this study, 12 primary malignant breast tissues (8PR+ and 4 PR-) were developed into patient-derived explants and treated with Vehicle, E2, E2+R5020, or R5020 for 24 or 48 hrs.
#' Raw data for 48-hr Vehicle-, R5020-treated and E2+R5020-treated samples were retrieved from GEO (GSE80098) and pre-processed into raw count. Filtration was sequentially performed toremove undetectable genes
#' and the filtered counts were normalised using [conditional quantile normalisation](bioconductor.org/packages/devel/bioc/vignettes/cqn/inst/doc/cqn.pdf) to offset effects of systematic artefacts,
#' such as gene length and GC contents. To reduce computing time, we randomly sampled half of the genes after filtration and used their logCPM value as the example data.
#'
#' @format A matrix with 7672 rows and 15 columns
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80098}
#' @usage data(logCPM_example)
"logCPM_example"


#' metadata_example: Sample metadata for malignant breast cancer tumors PDE from 5 ER+ breast cancer patients (GSE80098)
#'
#' @format A data frame with 15 rows and 4 columns
#' \describe{
#'   \item{patient}{patient N2-3, P4-6 }
#'   \item{treatment}{treatment: Vehicle, E2 or E2+R5020}
#'   \item{PR}{progesterone receptor status}
#'   \item{sample}{sample name, corresponding to column names of the logCPM matrix}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4928895/}
#' @usage data(metadata_example)
"metadata_example"
