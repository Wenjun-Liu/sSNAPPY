#' Plot genes' contribution to a specific pathway's perturbation as heatmap
#' @param genePertScore List of gene-wise perturbation scores generated using function `raw_gene_pert()`
#' @param gsToPlot `character` Name of the pathway to be plotted
#' @param mapEntrezID Optional. A `data.frame` matching genes' entrez IDs to other identifiers, such as gene names.
#' Must contain 2 columns: `"entrezid"` and `"mapTo"`
#' @param metadata  A `data.frame` containing sample metadata for heatmap annotation
#' @param annotation_attribute `character` Vector specifying attributes to draw annotations for. Default to "pathwayPertScore" (
#' ie. pathway-level perturbation scores)
#' @param pathwayPertScore A `data.frame` containing pathway-level perturbation scores for each pathway each treated sample. Output of function `pathway_pert()`
#' @param ... Used to pass various potting parameters to `pheatmap::pheatmap()`
#' @details The single-sample pathway-level perturbation score for a given pathway is derived from aggregating all the gene-wise perturbation scores of genes in that pathway.
#' This function visualises individual pathway genes' perturbation scores as a heatmap to demonstrate genes' contribution to a pathway perturbation.
#'
#' Plotting of the heatmap is done through `pheatmap::pheatmap()` so all plotting parameters accepted by `pheatmap::pheatmap()` could also be passed to this function.
#'
#' It is recommended to provide the pathway-level perturbation scores derived using the `pathway_pert()` function to visualise the directions of changes at pathway-level
#' as a column annotation, which helps the identification of genes driving or antagonizing the perturbation.
#'
#' Additional annotation attributes could be specified through the `annotation_attribute` parameter and the specified attributes must be provided by columns of the sample
#' metadata `data.frame` provided through the `metadata` parameter, otherwise the attributes will be ignored.
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr filter mutate
#' @importFrom pheatmap pheatmap
#' @references Kolde R (2019). _pheatmap: Pretty Heatmaps_. R package version 1.0.12, <https://CRAN.R-project.org/package=pheatmap>.
#' @examples
#' #compute weighted single sample logFCs
#' data(metadata_example)
#' data(logCPM_example)
#'
#' # compute single-sample logFCs for all treated samples
#' ls <- weight_ss_fc(logCPM_example, metadata = metadata_example,
#' factor = "patient", control = "Vehicle")
#'
#' # extract all the KEGG pathways
#' gsTopology <- retrieve_topology(database = "kegg")
#'
#' # compute raw gene-wise perturbation scores
#' genePertScore <- raw_gene_pert(ls$logFC, gsTopology)

#' # sum gene-wise perturbation scores to derive the pathway-level single-sample perturbation scores
#' pathwayPertScore <- pathway_pert(genePertScore)
#'
#' # Genes' contribution to the perturbation of Estrogen signaling pathway was
#' # visuaulised with pathway-level perturbation scores
#' # and treatments as column annotation attributes.
#' plot_gene_contribution(genePertScore = genePertScore, gsToPlot =
#' "Estrogen signaling pathway", metadata = metadata_example,
#' annotation_attribute = c("pathwayPertScore", "treatment"),
#' pathwayPertScore = pathwayPertScore)
#' @export
plot_gene_contribution <- function(genePertScore, gsToPlot, mapEntrezID = NULL, metadata = NULL, annotation_attribute = "pathwayPertScore",  pathwayPertScore = NULL,
                                ...){
    gs_name <- tA <- weight <- entrezid <- NULL
    if (!gsToPlot %in% names(genePertScore)  )
        stop("Gene-wise perturbation scores not provided for the chosen pathway.")

    rankMatrix <- genePertScore[[gsToPlot]]

    if (is.null(annotation_attribute)){
        pheatmap::pheatmap(
            rankMatrix,
            ...
        )
    } else{
        if (!identical(annotation_attribute, "pathwayPertScore")){
            other_attr <- setdiff(annotation_attribute, "pathwayPertScore")

            if (any(other_attr %in% colnames(metadata))){
                    anno_col_df <- dplyr::filter(metadata, sample %in% colnames(rankMatrix))
                    anno_col_df <- column_to_rownames(anno_col_df, "sample")
                    anno_col_df <- anno_col_df[,setdiff(annotation_attribute, "pathwayPertScore"), drop = FALSE]
            }

            if ("pathwayPertScore" %in% annotation_attribute){
                if (is.null(pathwayPertScore)){
                    warning("To use pathway-level perturbation as annotation, pathway-level annotation scores much be provided. Parameter ignored.")
                    if (exists("anno_col_df")){anno_col_df <- anno_col_df} else {anno_col_df <- NULL}
                } else {
                    if (exists("anno_col_df")) {
                        anno_col_df2 <- dplyr::filter(pathwayPertScore, gs_name == gsToPlot)
                        anno_col_df2 <- mutate(anno_col_df2, `Pathway-level Perturbation` = ifelse(tA < 0, "Inhibited", "Activated"))
                        anno_col_df2 <- tibble::column_to_rownames(anno_col_df2[,c("sample", "Pathway-level Perturbation")], "sample")
                        anno_col_df <- cbind(anno_col_df, anno_col_df2)
                    } else {
                        anno_col_df <- dplyr::filter(pathwayPertScore, gs_name == gsToPlot)
                        anno_col_df <- mutate(anno_col_df, `Pathway-level Perturbation` = ifelse(tA < 0, "Inhibited", "Activated"))
                        anno_col_df <- tibble::column_to_rownames(anno_col_df[,c("sample", "Pathway-level Perturbation")], "sample")
                    }
                }
            }

        } else {
            if (is.null(pathwayPertScore) | !gsToPlot %in% pathwayPertScore$gs_name){
                warning("To use pathway-level perturbation as annotation, pathway-level annotation scores much be provided. Parameter ignored.")
                anno_col_df <- NULL
            } else {
                anno_col_df <- dplyr::filter(pathwayPertScore, gs_name == gsToPlot)
                anno_col_df <- mutate(anno_col_df, `Pathway-level Perturbation` = ifelse(tA < 0, "Inhibited", "Activated"))
                anno_col_df <- tibble::column_to_rownames(anno_col_df[,c("sample", "Pathway-level Perturbation")], "sample")
            }

        }

        # if (!is.null(colOrder)){
        #     rankMatrix <- rankMatrix[,colOrder, drop = FALSE]
        # }

        if (all(c("entrezid","mapTo") %in% colnames(mapEntrezID))){
            if (any(rownames(rankMatrix) %in% mapEntrezID$entrezid)){
                mapEntrezID <- dplyr::filter(mapEntrezID, entrezid %in% rownames(rankMatrix))
                rankMatrix <- rankMatrix[rownames(rankMatrix) %in% mapEntrezID$entrezid,]
                mapEntrezID <- mapEntrezID[match(rownames(rankMatrix), mapEntrezID$entrezid), ]
                rownames(rankMatrix) <- mapEntrezID$mapTo
            } else {
                warning("None of the Entrez IDs in mapEntrezID mapped to gsTopology.\n
                        Still using entrez IDs as rownames.")

            }
        }

        pheatmap::pheatmap(
            # tibble::column_to_rownames(geneRank$`Estrogen signaling pathway`[,colOrder], "gene_id"),
            rankMatrix,
            annotation_col = anno_col_df,
            ...
        )
    }


}
