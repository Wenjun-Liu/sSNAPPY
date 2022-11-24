#' Plot pathways and genes contained in them as a network
#'
#' @details Taking the perturbation scores of a list of gene-sets derived from `normalise_by_permu()` as input, this function matches gene-set to
#' their associated genes by utilizing information from pathway topology matrices.
#'
#' It's optional to provide genes' logFCs as a named vector, where the names must be genes' entrez IDs in the format of "ENTREZID:XXXX".
#' This is because pathway topology matrices retrieved through `retrieve_topology()` always use entrez ID as identifiers.
#'
#' However, it might not be very informative to label genes with their entrez ID. So users can also choose to proivde a `mapEntrezID`
#' `data.frame` to match entrez IDs to their chosen identifiers. The `data.frame` should contain two columns: `"entrezid"` and `"mapTo"`.
#'
#' If `geneFC` is provided, gene nodes will be colored by genes' directions of changes. Otherwise, all gene nodes will be black.
#'
#' Since some gene-sets could can contain hundreds of genes, it is not recommended to plot all of those genes. If `mapEntrezID` `data.frame` is provided,
#' only genes included in that `data.frame` will be used in the plot.
#'
#' It is strongly recommended to filter for genes with the highest magnitude of changes. If all pathway genes
#' have to be plotted, consider setting `label_Gene` to FALSE to turn off plotting all gene names.
#'
#' @param normalisedScores A `data.frame` derived from the `normalise_by_permu()` function. Only gene-sets of interest should be included
#' @param gsTopology List of pathway topology matrices generated using function `retrieve_topology()`
#' @param geneFC An optional named vector of pathways' fold changes
#' @param mapEntrezID Optional. A `data.frame` matching genes' entrez IDs to other identifier. Must contain 2 columns: `"entrezid"` and `"mapTo"`
#' @param colorGS_By Choose to color nodes by *robustZ* or *pvalue*. A column must exist in the normalisedScores `data.frame`
#' for the chosen parameter
#' @param foldGSname `logical`. Should long gene-set names be folded into two lines
#' @param foldafter The number of words after which gene-set names should be folded. Defaulted to 2
#' @param layout The layout algorithm to apply. Accept all layout supported by `igraph`.
#' @param edgeAlpha Transparency of edges. Default to 0.8
#' @param upGS_col Color for activated gene-sets. Only applicable if `colorGS_By` is set to be "robustZ"
#' @param downGS_col Color for inhibited gene-sets. Only applicable if `colorGS_By` is set to be "robustZ"
#' @param upGene_col Color for up-regulated genes. Only applicable if `geneFC` is not NULL
#' @param downGene_col Color for down-regulated genes. Only applicable if `geneFC` is not NULL
#' @param GeneNode_size Size for gene nodes
#' @param GeneNode_shape Shape for gene nodes
#' @param GsNode_size Size for gene-set nodes
#' @param GsNode_shape Shape for gene nodes
#' @param label_Gene `logical` Should gene name be plotted
#' @param GeneName_size Size of gene name label
#' @param GsName_size Size of gene-set name label
#' @param gene_lg_title `character. Legend for gene nodes color
#' @param gs_lg_title `character` Legend for gene-set nodes color
#' @param arc_strength The bend of edges. 1 approximates a halfcircle while 0 will give a straight line.
#' @importFrom ggraph geom_edge_arc geom_node_label
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' load(system.file("extdata", "gsTopology.rda", package = "sSNAPPY"))
#' load(system.file("extdata", "normalisedScores.rda", package = "sSNAPPY"))
#' #Subset pathways significantly perturbed in sample R5020_N2_48
#' subset <- dplyr::filter(normalisedScores, adjPvalue < 0.05, sample == "R5020_N2_48")
#'
#' # Color gene-sets nodes by robust z-scores.
#' plot_gs2gene(subset, gsTopology, colorGS_By = "robustZ", label_Gene = FALSE,
#' GeneNode_size = 1)
#' # When genes' fold-changes are not provided, gene nodes are colored in black.
#'
#' # To color genes by their directions of changes, firstly compute genes' single-sample logFCs
#' data(logCPM_example)
#' data(metadata_example)
#' ls <- weight_ss_fc(logCPM_example, metadata = metadata_example,
#'  factor = "patient", control = "Vehicle")
#' # Provide fold-changes of sample R5020_N2_48 as a named vector
#' plot_gs2gene(subset, gsTopology, geneFC = ls$logFC[,"R5020_N2_48"], colorGS_By = "robustZ",
#' label_Gene = FALSE)
#'
#' # There are still a large number of genes, making the plot cumbersome. There only fold-changes of
#' # genes with top 500 absolute fold-changes are provide so only pathway genes in that list of 500
#' # genes were plotted.
#' FC <- sort(abs(ls$logFC[,"R5020_N2_48"]), decreasing = TRUE)[1:500]
#' plot_gs2gene(subset, gsTopology, geneFC = FC, colorGS_By = "robustZ")
#'
#' # To make the gene labels more informative, map genes' entrez id to chosen identifiers.
#' load(system.file("extdata", "entrez2name.rda", package = "sSNAPPY"))
#' plot_gs2gene(subset, gsTopology, geneFC = FC, mapEntrezID = entrez2name, colorGS_By = "robustZ")
plot_gs2gene <- function(normalisedScores, gsTopology, geneFC = NULL, mapEntrezID = NULL, colorGS_By = c("robustZ", "pvalue"), foldGSname = TRUE, foldafter = 2, layout = "fr",
                         edgeAlpha = 0.8,  upGS_col = "brown3", downGS_col = "steelblue3", upGene_col = "pink", downGene_col = "lightblue", GeneNode_size = 3, GeneNode_shape = 17,
                         GsNode_size = 2, GsNode_shape = 16,label_Gene = TRUE, GeneName_size = 3, GsName_size = 6,  gene_lg_title = "Changes in Gene Expression",
                         gs_lg_title = "Pathway Perturbation", arc_strength = 0.5){

    name <- weight <- color <- size <- entrezid <- mapTo <- gs_Name <- . <- type <-  NULL
    ## check if input has required columns
    stopifnot(colorGS_By %in% c("robustZ", "pvalue"))
    if (!all(c(colorGS_By, "gs_name") %in% colnames(normalisedScores))) stop("Normalised Scores must include gs_name and column for coloring")

    # Make sure the gs topologies are a named list with at least two elements
    stopifnot(length(names(gsTopology)) == length(gsTopology))
    gsTopology <- gsTopology[names(gsTopology) %in% normalisedScores$gs_name]

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
                                   stroke = 0.5, size = GsNode_size, shape = GsNode_shape) +
            scale_color_manual(values = c("Activated" = upGS_col, "Inhibited" = downGS_col), name =  gs_lg_title) +
            ggraph::geom_node_label(data = . %>% dplyr::filter(type == "GS"), aes(label = name, color = color),
                                    size = GsName_size, repel = TRUE, fill = NA)
    } else(
        pl <- pl +
            geom_node_point(data = .%>% dplyr::filter(type == "GS"), aes(color = as.numeric(color)),
                            stroke = 0.5, size = GsNode_size, shape = GsNode_shape) +
            scale_color_continuous(low="red", high="blue", name = gene_lg_title) +
            ggraph::geom_node_label(data = . %>% dplyr::filter(type == "GS"), aes(label = name,
                                                                                  color = as.numeric(color)),
                                    size = GsName_size, repel = TRUE, fill = NA)
        )

    # plot gene nodes
    pl <- pl + ggnewscale::new_scale_color() +
        geom_node_point(data = .%>% dplyr::filter(type == "GENE"), aes(color = color),
                               stroke = 0.5, size = GeneNode_size, shape = GeneNode_shape) +
        scale_color_manual(values = c("Up-regulated" = upGene_col, "Down-regulated" = downGene_col, "NoFC" = "black"), name =  gene_lg_title)

    if (label_Gene){
        pl <- pl + geom_node_text(data = .%>% dplyr::filter(type == "GENE"),aes(label = name), size = GeneName_size, repel = TRUE )
    }
    pl

}

#' @importFrom reshape2 melt
#' @importFrom purrr set_names
#' @importFrom igraph E V graph.data.frame set_edge_attr set_vertex_attr degree delete_vertices delete_edges
make_gs2gene_network <- function(normalisedScores, gsTopology, geneFC, upGene_col, downGene_col,
                                 colorGS_By = c("robustZ", "pvalue"), mapEntrezID){

    # create dummy variable to pass R CMD CHECK
    from <- to <- E <- robustZ <- entrezid <- mapTo <- gs_name <- NULL
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
        mapEntrezID <- dplyr::filter(unique(dplyr::select(mapEntrezID , entrezid, mapTo)), entrezid %in% names(gene_col))
        mapEntrezID <-  mapEntrezID[match(names(gene_col), mapEntrezID$entrezid),]
        names(gene_col) <- mapEntrezID$mapTo
        g <- graph.data.frame(dplyr::select(GS2Gene, gs_name, mapTo), directed = FALSE)
    } else{
        g <- graph.data.frame(dplyr::select(GS2Gene, gs_name, entrezid), directed = FALSE)
    }

    # set color for gene nodes
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
