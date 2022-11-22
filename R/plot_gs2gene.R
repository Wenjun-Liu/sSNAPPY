plot_gs2gene <- function(normalisedScores, gsTopology, geneFC = NULL, mapEntrezID = NULL, colorGS_By = c("robustZ", "pvalue"), foldGSname = TRUE, foldafter = 2, layout = "fr",
                         edgeAlpha = 0.8,  upGS_col = "brown3", downGS_col = "steelblue3", upGene_col = "pink", downGene_col = "lightblue", GeneNode_size = 3, GeneNode_shape = 17,
                         GsNode_size = 2, GsNode_shape = 16,label_Gene = TRUE, GeneName_size = 3, GsName_size = 6, color_lg = TRUE, gene_lg_title = "Changes in Gene Expression",
                         gs_lg_title = "Pathway Perturbation", arc_strength = 0.5, ...){

    name <- weight <- color <- size <- NULL
    ## check if input has required columns
    stopifnot(colorGS_By %in% c("robustZ", "pvalue"))
    if (!all(c(colorGS_By, "gs_name") %in% colnames(normalisedScores))) stop("Normalised Scores must include gs_name and column for coloring")

    # Make sure the gs topologies are a named list with at least two elements
    stopifnot(length(names(gsTopology)) == length(gsTopology))


    # create igraph object
    g <- make_gs2gene_network(normalisedScores, gsTopology, geneFC, upGene_col, downGene_col, colorGS_By,
                        mapEntrezID)

    if(foldGSname){
        g <- set_vertex_attr(g, "name", value = vapply(V(g)$name, function(x){ifelse(length(strsplit(x, " ")[[1]]) > foldafter,
                                                                                     str_replace_nth(x, " ", "\n", foldafter),
                                                                                     x)}, character(1))) }
    pl <- ggraph(g, layout = layout) +
            ggraph::geom_edge_arc(alpha = edgeAlpha, colour='darkgrey', strength = arc_strength)

    # plot gene-set node
    if (colorGS_By == "robustZ"){
        pl <- pl + geom_node_point(data = .%>% dplyr::filter(type == "GS"), aes(color = color),
                                   stroke = 0.5, size = GeneNode_size, shape = GsNode_shape) +
            scale_color_manual(values = c("Activated" = upGS_col, "Inhibited" = downGS_col), name =  gs_lg_title)
    } else(
        pl <- pl +
            geom_node_point(data = .%>% dplyr::filter(type == "GS"), aes(color = color),
                            stroke = 0.5, size = GeneNode_size) +
            scale_color_continuous(low="red", high="blue", name = gene_lg_title)
        )

    # label gene-set names
    pl <- pl + ggraph::geom_node_label(data = . %>% dplyr::filter(type == "GS"), aes(label = name, color = color), size = GsName_size,
                                       repel = TRUE, fill = NA, ...)

    # plot gene nodes
    pl <- pl + ggnewscale::new_scale_color() +
        geom_node_point(data = .%>% dplyr::filter(type == "GENE"), aes(color = color),
                               stroke = 0.5, size = GsNode_size, shape = GeneNode_shape) +
        scale_color_manual(values = c("Up-regulated" = upGene_col, "Down-regulated" = downGene_col, "NoFC" = "black"), name =  gs_lg_title)

    if (label_Gene){
        pl <- pl + geom_node_text(data = .%>% dplyr::filter(type == "GENE"),aes(label = name), size = GeneName_size, repel = TRUE )
    }
    if(!color_lg){
        pl <- pl + guides(color = "none")
    }

    pl

}

#' @importFrom reshape2 melt
#' @importFrom igraph E V graph.data.frame set_edge_attr set_vertex_attr degree delete_vertices delete_edges
make_gs2gene_network <- function(normalisedScores, gsTopology, geneFC, upGene_col, downGene_col,
                                 colorGS_By = c("robustZ", "pvalue"), mapEntrezID){

    # create dummy variable to pass R CMD CHECK
    from <- to <- E <- robustZ <- NULL
    gsTopology <- gsTopology[names(gsTopology) %in% normalisedScores$gs_name]
    GS2Gene <- get_GSgenelist(gsTopology, mapEntrezID)
    GS2Gene <- left_join(normalisedScores, GS2Gene, by = "gs_name")

    # set color for gene nodes
    if(!is.null(geneFC)){
        if(!is.vector(geneFC) | is.null(names(geneFC)) | length(intersect(GS2Gene$entrezid, names(geneFC))) == 0){
            warning("Genes's logFCs not provided as a named vector. All genes will be colored identically.")
            gene_col <- rep("NoFC", times = nrow(GS2Gene))
            gene_col <- set_names(gene_col, GS2Gene$entrezid)
        } else {
            GS2Gene <- dplyr::filter(GS2Gene, entrezid %in% names(geneFC))
            geneFC <- geneFC[names(geneFC) %in% GS2Gene$entrezid]
            gene_col <- ifelse(geneFC > 0, "Up-regulated", "Down-regulated")
        }

    } else{
        gene_col <- rep("NoFC", times = nrow(GS2Gene))
        gene_col <- set_names(gene_col, GS2Gene$entrezid)
    }



    if (!is.null(mapEntrezID)){
        mapEntrezID <- dplyr::filter(unique(dplyr::select(GS2Gene, entrezid, mapTo)), entrezid %in% names(geneFC))
        mapEntrezID <-  mapEntrezID[match(names(gene_col), mapEntrezID$entrezid),]
        names(gene_col) <- mapEntrezID$mapTo
        g <- graph.data.frame(dplyr::select(GS2Gene, gs_name, mapTo), directed = FALSE)
    } else{
        g <- graph.data.frame(dplyr::select(GS2Gene, gs_name, entrezid), directed = FALSE)
    }

    g <- set_vertex_attr(g, "color", index = names(gene_col), value = gene_col)

    # set color for gene-set nodes

    if (colorGS_By == "robustZ"){
        GScolor <- mutate(GS2Gene, color =  ifelse(robustZ < 0, "Inhibited", "Activated"))
        g <- set_vertex_attr(g, "color", index = GScolor$gs_name, value = GScolor$color)
    }

    if (colorGS_By == "pvalue"){
        GSpvalue <- unique(GS2Gene[,c("gs_name", "pvalue")])
        g <- set_vertex_attr(g, "color", index = GSpvalue$gs_name, value = GSpvalue$pvalue)
    }

    # set types of nodes
    nodeType <- ifelse(V(g)$name %in% GS2Gene$gs_name, "GS", "GENE")
    g <- set_vertex_attr(g, "type", value = nodeType)

    g
}
