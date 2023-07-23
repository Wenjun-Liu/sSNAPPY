#' @title Plot genes' contribution to a pathway's perturbation as a heatmap
#'
#' @description
#' Plot individual genes' contributions to the pathway-level perturbation score
#'
#'
#' @param genePertMatr A matrix of gene-wise perturbation scores corresponding
#' to a pathway. An element of the output generated using function `
#' raw_gene_pert()`
#' @param mapEntrezID Optional. A `data.frame` matching genes' entrez IDs to
#' another identifier with preferred labels. Must contain the columns:
#' `"entrezid"` and `"mapTo"`
#' @param topGene Numeric(1). The number of top genes to plot
#' @param filterBy Filter top genes by the mean, variability (sd), maximum
#' value, or maximum absolute values
#' @param tieMethod Method for handling ties in ranking (i.e. in values returned
#' by `filterBy`, two or many genes share the same value). Default to "min". See
#' `?rank` for other options.
#' @param annotation_df  A `data.frame` for annotating heatmap columns. Must
#' contain a "sample" column with sample names matching to the column names of
#' the `genePertMatr`
#' @param ... Used to pass various potting parameters to
#' [`pheatmap::pheatmap()`]
#'
#'
#' @details The single-sample pathway-level perturbation score for a given
#' pathway is derived from aggregating all the gene-wise perturbation scores of
#' genes in that pathway. This function visualises individual pathway genes'
#' perturbation scores as a heatmap to demonstrate pathway genes' contribution
#' to a pathway perturbation.
#'
#' Plotting of the heatmap is done through [`pheatmap::pheatmap()`] so all
#' plotting parameters accepted by [`pheatmap::pheatmap()`] could also be passed
#' to this function.
#'
#' @references Kolde R (2019). _pheatmap: Pretty Heatmaps_. R package version
#' 1.0.12, <https://CRAN.R-project.org/package=pheatmap>.
#' @examples
#' #compute weighted single sample logFCs
#' data(metadata_example)
#' data(logCPM_example)
#' metadata_example <- dplyr::mutate(metadata_example, treatment = factor(
#'    treatment, levels = c("Vehicle", "E2+R5020", "R5020")))
#' # compute single-sample logFCs for all treated samples
#' ls <- weight_ss_fc(logCPM_example, metadata = metadata_example,
#' groupBy = "patient", treatColumn = "treatment", sampleColumn = "sample")
#'
#' # extract all the KEGG pathways
#' gsTopology <- retrieve_topology(database = "kegg", species = "hsapiens")
#'
#' # compute raw gene-wise perturbation scores
#' genePertScore <- raw_gene_pert(ls$weighted_logFC, gsTopology)

#' # sum gene-wise perturbation scores to derive the pathway-level single-sample perturbation scores
#' pathwayPertScore <- pathway_pert(genePertScore)
#'
#' # Genes with top 10 mean absolute gene-wise perturbation scores in the
#' # Estrogen signaling pathway was visualised.
#' plot_gene_contribution(genePertScore$`kegg.Estrogen signaling pathway`,
#' filterBy = "mean", topGene = 10)
#'
#' # Columns of the heatmap could be annotated by the pathway-level perturbation
#' # and treatments. Firstly, create a `data.frame` with the two annotation
#' # attributes and sample names matching the column names of the perturbation
#' # score matrix.
#' annotation_df <- dplyr::select(metadata_example, sample, treatment)
#' pathwayLevel <- dplyr::filter(pathwayPertScore,
#' gs_name == "kegg.Estrogen signaling pathway")
#' pathwayLevel$`pathway-level` <- ifelse(
#' pathwayLevel$score > 0, "Activated", "Inhibited")
#' annotation_df <- dplyr::left_join(
#' dplyr::select(pathwayLevel, sample, `pathway-level`),
#' annotation_df, unmatched = "drop")
#' # To make the gene labels more informative, also map genes' entrez id
#' # to chosen identifiers.
#' load(system.file("extdata", "entrez2name.rda", package = "sSNAPPY"))
#' plot_gene_contribution(genePertScore$`kegg.Estrogen signaling pathway`,
#' topGene = 10, filterBy = "mean", annotation_df = annotation_df,
#' mapEntrezID = entrez2name)
#'
#' # Plotting parameters accepted by `pheatmap::pheatmap()` could be passed to
#' # this function to customise the plot. For example, changin the color of annotations
#' plot_gene_contribution(genePertScore$`kegg.Estrogen signaling pathway`,
#' topGene = 10, filterBy = "mean", annotation_df = annotation_df,
#' mapEntrezID = entrez2name, annotation_colors = list(
#' treatment = c(R5020 = "black", `E2+R5020` = "white"),
#' `pathway-level` = c(Activated = "darkred", Inhibited = "lightskyblue")))
#' @importFrom tibble column_to_rownames
#' @importFrom pheatmap pheatmap
#' @export
plot_gene_contribution <- function(
        genePertMatr, mapEntrezID = NULL, topGene = 10,
        filterBy = c("mean", "sd", "max.abs"), tieMethod = "min",
        annotation_df = NULL, ...
){
    entrezid <- NULL
    filterBy <- match.arg(filterBy)
    if (filterBy == "max.abs") {
        f <- function(x) {max(abs(x))}
    } else {
        f <- match.fun(filterBy)
    }
    stopifnot(is.numeric(topGene))
    topGene <- min(topGene, nrow(genePertMatr))

    # filter genePertMatrx
    vals <- apply(genePertMatr, 1, f)
    topRanked <- rank(1/abs(vals), ties.method = tieMethod ) <= topGene
    genePertMatr <- genePertMatr[topRanked,]
    ids <- rownames(genePertMatr)

    # match Entrez ID to other gene identifiers
    if (all(c("entrezid","mapTo") %in% colnames(mapEntrezID))) {
        if (any(rownames(genePertMatr) %in% mapEntrezID$entrezid)) {
            mapEntrezID <- dplyr::filter(mapEntrezID, entrezid %in% ids)
            genePertMatr <- genePertMatr[ids %in% mapEntrezID$entrezid,]
            map <- match(rownames(genePertMatr), mapEntrezID$entrezid)
            mapEntrezID <- mapEntrezID[map, ]
            rownames(genePertMatr) <- mapEntrezID$mapTo
        } else {
            warning(
                "None of the EntrezIDs in mapEntrezID mapped to gsTopology.\n",
                "The original EntrezIDs will be retained as rownames."
            )

        }
    }

    if (is.null(annotation_df)) {
        pheatmap(genePertMatr, ...)
    } else{

        if (!"sample" %in% colnames(annotation_df) |
            any(!colnames(genePertMatr) %in% annotation_df$sample)) {
            message(
                "Column names of the perturbation score matrix must match
                the sample column of the annotation_df. Annotation df ignored."
            )
            pheatmap(genePertMatr, ...)
        } else {
            annotation_df <- column_to_rownames(annotation_df, "sample")
            pheatmap(genePertMatr, annotation_col = annotation_df, ...)
        }

    }

}
