#' logCPM_example: Normalised logCPM of patient-derived explant models obtained from 12 ER-positive primamry breast cancer patients (GSE80098)
#'
#' This data was adopted from a study by Singhal H, et al., which was published as *Genomic agonism and phenotypic antagonism between estrogen and progesterone receptors in breast cancer*
#' in 2016. In this study, 12 primary malignant breast tissues (8PR+ and 4 PR-) were developed into patient-derived explants and treated with Vehicle, E2, E2+R5020, or R5020 for 24 or 48 hrs.
#' Raw data for Vehicle-, R5020-treated and E2+R5020-treated samples were retrieved from GEO (GSE80098) and pre-processed into raw count. Filtration was sequentially performed toremove undetectable genes
#' and the filtered counts were normalised using [conditional quantile normalisation](bioconductor.org/packages/devel/bioc/vignettes/cqn/inst/doc/cqn.pdf) to offset effects of systematic artefacts,
#' such as gene length and GC contents.
#'
#' @format A matrix with 18481 rows and 16 columns
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80098}
"logCPM_example"


#' metadata_example: Sample metadata for malignant breast cancer tumors PDE from 8 ER+PR+ breast cancer patients (GSE80098)
#'
#' @format A data frame with 16 rows and 6 columns
#' \describe{
#'   \item{patient}{patient p1-p8}
#'   \item{treatment}{treatment: Vehicle or E2}
#'   \item{time}{amount of time the sample was treated for, in hours}
#'   \item{sample}{sample name, corresponding to column names of the logCPM matrix}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4928895/}
"metadata_example"
