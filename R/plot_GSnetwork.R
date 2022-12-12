#' @title Plot significantly perturbed gene-sets as a network
#'
#' @param normalisedScores A `data.frame` of pathway perturbation scores derived from the `normalise_by_permu()` function
#' @param gsTopology List of pathway topology matrices generated using function `retrieve_topology()`
#' @param colorBy Choose to color nodes either by *robustZ* or *pvalue*. A column must exist in the `normalisedScores`
#' `data.frame` for the chosen parameter
#' @param foldGSname `logical` Should long gene-set names be folded into two lines
#' @param foldafter The number of words after which gene-set names should be folded. Defaulted to 2
#' @param layout The layout algorithm to apply. Accept all layout supported by `igraph`
#' @param edgeAlpha `numerical` Transparency of edges. Default to 0.8
#' @param up_col The color used to label activated gene-sets. Only applicable if `colorBy` is set to be "robustZ"
#' @param down_col The color used to label inhibited gene-sets. Only applicable if `colorBy` is set to be "robustZ"
#' @param scale_edgeWidth A numerical vector of length 2 to be provided to `ggraph::scale_edge_width_continuous()` for specifying
#' the minimum and maximum edge widths after transformation. Defaulted to `c(0.5, 3)`
#' @param edgeLegend logical` Should edge weight legend be shown
#' @param scale_nodeSize A numerical vector of length 2 to be provided to `ggplot2::scale_size()` for specifying
#' the minimum and maximum node sizes after transformation. Defaulted to `c(3,6)`
#' @param nodeShape Shape of nodes
#' @param color_lg `logical` Should color legend be shown
#' @param color_lg_title Optional. Title for the color legend
#' @param lb_size Size of node text labels
#' @param lb_color Color of node text labels
#' @param plotIsolated `logical` Should nodes not connected to any other nodes be plotted.  Default to FALSE
#' @importFrom ggraph ggraph geom_node_point geom_node_text geom_edge_link
#' @importFrom ggplot2 scale_color_continuous scale_color_manual aes_ guides aes scale_size
#' @importFrom ggraph scale_edge_width_continuous
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' load(system.file("extdata", "gsTopology.rda", package = "sSNAPPY"))
#' load(system.file("extdata", "normalisedScores.rda", package = "sSNAPPY"))
#'
#' #Subset pathways significantly perturbed in sample R5020_N2_48
#' subset <- dplyr::filter(normalisedScores, adjPvalue < 0.05, sample == "R5020_N2_48")
#'
#' # Color network plot nodes by robust z-score
#' plot_gs_network(subset, gsTopology,
#' colorBy = "robustZ", layout = "dh",
#' color_lg_title = "Direction of pathway Perturbation")
#'
#' # Color network plot nodes by p-values
#' plot_gs_network(subset, gsTopology, layout = "dh",
#' colorBy = "pvalue", color_lg_title = "P-value")
plot_gs_network <- function(normalisedScores, gsTopology, colorBy = c("robustZ", "pvalue"), foldGSname = TRUE, foldafter = 2, layout = "fr",
                           edgeAlpha = 0.8,  up_col = "brown3", down_col = "steelblue3", scale_edgeWidth = c(0.5, 3), edgeLegend = FALSE, scale_nodeSize = c(3,6),
                           nodeShape = 16, color_lg = TRUE, color_lg_title = NULL, lb_size = 3, lb_color = "black", plotIsolated = FALSE){

    name <- weight <- color <- size <- NULL
    ## check if input has required columns
    stopifnot(colorBy %in% c("robustZ", "pvalue"))
    if (!all(c(colorBy, "gs_name") %in% colnames(normalisedScores))) stop("Normalised Scores must include gs_name and column for coloring")

    # Make sure the gs topologies are a named list with at least two elements
    stopifnot(length(names(gsTopology)) == length(gsTopology))
    if(length(unique(normalisedScores$gs_name)) < 2) stop("At least 2 gene-sets are required for a network plot")
    gsTopology <- gsTopology[names(gsTopology) %in% normalisedScores$gs_name]

    # create igraph object
    g <- make_gsNetwork(normalisedScores, gsTopology, colorBy = colorBy, plotIsolated = plotIsolated)

    if(foldGSname){
        g <- set_vertex_attr(g, "name", value = vapply(V(g)$name, function(x){ifelse(length(strsplit(x, " ")[[1]]) > foldafter,
                                                                                     str_replace_nth(x, " ", "\n", foldafter),
                                                                                     x)}, character(1))) }

    # plot network edges
    pl <- ggraph(g, layout = layout) +
        geom_edge_link(alpha = edgeAlpha, aes(width=weight), colour='darkgrey') +
        scale_edge_width_continuous(range = scale_edgeWidth, guide = "none")


    if (colorBy == "robustZ"){
        pl <- pl + geom_node_point(aes(color = color, size = size), shape = nodeShape,stroke = 0.5) +
            scale_color_manual(values = c("Activated" = up_col, "Inhibited" = down_col),name = color_lg_title) +
            scale_size(range =  scale_nodeSize, guide = "none")
    } else(
        pl <- pl +
            geom_node_point(aes(color = color, size = size), shape = nodeShape,stroke = 0.5) +
            scale_color_continuous(low="red", high="blue", name = color_lg_title) +
            scale_size(range =  scale_nodeSize, guide = "none"))

    if(!color_lg){
     pl <- pl + guides(color = "none")
    }


    pl + geom_node_text(aes(label = name), size = lb_size, repel = TRUE, colour = lb_color)

}

#' @importFrom reshape2 melt
#' @importFrom dplyr bind_rows
#' @importFrom igraph E V graph.data.frame set_edge_attr set_vertex_attr degree delete_vertices delete_edges
make_gsNetwork <- function(normalisedScores, gsTopology,  colorBy = c("robustZ", "pvalue", "community"),
                           plotIsolated ){

    # create dummy variable to pass R CMD CHECK
    from <- to <- E <- robustZ <- NULL
    GS2Gene <- get_GSgenelist(gsTopology)
    GS2Gene <- left_join(normalisedScores, GS2Gene, by = "gs_name")

    GSlist <- split(GS2Gene[,c("gs_name", "entrezid")], f = GS2Gene$gs_name)
    nGS <- length(GSlist)
    GSname <- names(GSlist)

    w <- lapply(seq_len(nGS-1), function(x){
        lapply((x+1):nGS, function(y){
            data.frame(from = GSname[x], to = GSname[y], weight = jacIdex_func(GSlist[[x]]$entrezid, GSlist[[y]]$entrezid))
        })
    })

    w <- bind_rows(lapply(w, bind_rows))
    w <- filter(w, from != to)

    g <- graph.data.frame(dplyr::select(w, from, to), directed = FALSE)
    g <- set_edge_attr(g, "weight", value = w$weight)

    GSsize <- melt(lapply(GSlist, nrow))
    colnames(GSsize) <- c("size", "from")
    g <- set_vertex_attr(g, "size", index = GSsize$from, value = GSsize$size)

    if (colorBy == "robustZ"){
        GScolor <- mutate(GS2Gene, color =  ifelse(robustZ < 0, "Inhibited", "Activated"))
        g <- set_vertex_attr(g, "color", index = GScolor$gs_name, value = GScolor$color)
    }

    if (colorBy == "pvalue"){
        GSpvalue <- unique(GS2Gene[,c("gs_name", "pvalue")])
        g <- set_vertex_attr(g, "color", index = GSpvalue$gs_name, value = GSpvalue$pvalue)
    }

    if(!plotIsolated){
        removeEdge <- which(E(g)$weight == 0)
        g <-  delete_edges(g, removeEdge)
        IsolatedNode <- which(degree(g) == 0)
        g <- delete_vertices(g, IsolatedNode)
    }


    g
}




#' @importFrom reshape2 melt
get_GSgenelist <- function(gsTopology, mapEntrezID = NULL){
    GStoGene <- lapply(gsTopology, rownames)
    GStoGene <- reshape2::melt(GStoGene)
    colnames(GStoGene) <- c("entrezid", "gs_name")

    if (all(c("entrezid","mapTo") %in% colnames(mapEntrezID))){
        if (any(GStoGene$entrezid %in% mapEntrezID$entrezid)){
            left_join(GStoGene, mapEntrezID[,c("entrezid","mapTo")], by = "entrezid")
        } else {
            warning("None of the Entrez IDs in mapEntrezID mapped to gsTopology.")
            return(GStoGene)

        }
    } else {
        return(GStoGene)
    }
    }


jacIdex_func <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}

# This function was a modified version of the str_replace_nth function from martinctc/textworks
str_replace_nth <- function(x, pattern, replacement, n) {
    g <- gregexpr(pattern, x)[[1]][n]
    s <- scan(text = gsub("[()]", "", pattern),
              sep = "|",
              what = "",
              quiet = TRUE)
    substr(x, g, g) <- replacement[match(substr(x, g, g), s)]
    x
}
