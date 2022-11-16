

plot_gene_contribution <- function(geneRankList, gsToPlot, rownames, metadata, annotation_attribute = c("pathwayPertScore"), weightedFC, pathwayPertScore,
                                   colOrder){

    rankMatrix <- tibble::column_to_rownames(geneRankList[[gsToPlot]], "gene_id")

    if (!is.null(rownames)){
        rownames(rankMatrix) <- rownames
    }
    anno_col_df <- dplyr::filter(metadata, sample %in% colnames(rankMatrix))
    anno_col_df <- tibble::column_to_rownames(anno_col_df, "sample")
    anno_col_df <- anno_col_df[,setdiff(annotation_attribute, "pathwayPertScore"), drop = FALSE]

    if ("pathwayPertScore" %in% annotation_attribute & !is.null(pathwayPertScore)){
        pathwayPertScore <- dplyr::filter(pathwayPertScore, gs_name == gsToPlot)
        pathwayPertScore <- mutate(pathwayPertScore, dir = ifelse(tA < 0, "Inhibited", "Activated"))
        pathwayPertScore <- tibble::column_to_rownames(pathwayPertScore[,c("sample", "dir")], "sample")
    }

    anno_col_df <- cbind(anno_col_df, pathwayPertScore)
    pheatmap::pheatmap(
        # tibble::column_to_rownames(geneRank$`Estrogen signaling pathway`[,colOrder], "gene_id"),
        rankMatrix[,colOrder, drop = FALSE],
        annotation_col = anno_col_df
    )
}
