#' Plot genes' contribution to pathway perturbation as heatmap
#' @description
#' @param geneRankList List of gene-wise rank matrices generated using function `raw_gene_pert()`
#' @param gsToPlot `character` Name of pathway to be plotted
#' @param mapRownameTo `character` Rownames of heatmap. Default to `NULL`, which plots the rownames of the gene-wise rank matrix.
#' @param metadata Sample metadata data frame containing information for heatmap annotation
#' @param annotation_attribute `character` Vector specifying attributes to draw annotation for. Default to "pathwayPertScore" (
#' ie. pathway-level perturbation scores)
#' @param pathwayPertScore A dataframe. Output of function `pathway_pert()`
#' @param colOrder
#' @param ... Used to pass various potting parameters to `pheatmap::pheatmap()`
#' @details
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr filter mutate
#' @import pheatmap pheatmap
#' @return
#' @examples
#' #compute weighted single sample logFCs
#' data(metadata_example)
#' data(logCPM_example)
#' ls <- weight_ss_fc(logCPM_example, metadata = metadata_example,
#' factor = "patient", control = "Vehicle")
#' # extract all the KEGG pathways
#' gsTopology <- retrieve_topology(database = "kegg")
#' # compute raw gene-wise perturbation scores
#' genePertScore <- raw_gene_pert(ls$logFC, gsTopology)
#' # rank genes by gene-wise perturbation scores within each sample to compare their contributions to pathway perturbation
#' geneRank <- rank_gene_pert(ls$logFC,genePertScore,  gsTopology)
#' # sum gene-wise perturbation scores to derive the pathway-level single-sample perturbation scores
#' pathwayPertScore <- pathway_pert(ls$logFC, genePertScore)
#' plot_gene_contribution(geneRankList = geneRank, gsToPlot = "Estrogen signaling pathway", metadata = metadata_example,
#' annotation_attribute = c("pathwayPertScore", "treatment"))
#' @export
plot_gene_contribution <- function(geneRankList, gsToPlot, mapRownameTo = NULL, metadata = NULL, annotation_attribute = "pathwayPertScore",  pathwayPertScore = NULL,
                                   colOrder = NULL, annotation_colors = list(`Pathway-level Perturbation` = c("Inhibited" = "blue", "Activated" = "red")), ...){

    if (!gsToPlot %in% names(geneRankList)  )
        stop("Gene-wise perturbation score not provided for the chosen pathway.")

    rankMatrix <- column_to_rownames(geneRankList[[gsToPlot]], "gene_id")

    if (is.null(annotation_attribute)){
        pheatmap::pheatmap(
            rankMatrix[,colOrder, drop = FALSE],
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

        if (!is.null(colOrder)){
            rankMatrix <- rankMatrix[,colOrder, drop = FALSE]
        }

        if (!is.null(mapRownameTo)){
            rownames(rankMatrix) <- mapRownameTo
        }

        pheatmap::pheatmap(
            # tibble::column_to_rownames(geneRank$`Estrogen signaling pathway`[,colOrder], "gene_id"),
            rankMatrix,
            annotation_col = anno_col_df,
            annotation_colors = annotation_colors,

            ...
        )
    }


}
