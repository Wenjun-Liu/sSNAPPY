#' @title Plot gene-set network
#'
#' @param normalisedScores A dataframe as described in the details section
#' @param gsTopology List of pathway topology matrices generated using function `weightedAdjMatrix`
#' @param colorBy Choose to color nodes either by "robustZ" or "pvalue". A column must exist in the normalisedScores for the chosen parameter
#' @param foldGSname logical(1). Should long gene-set names be folded into two lines
#' @param foldafter The number of words after which gene-set names should be folded. Defaulted to 2
#' @param layout The layout algorithm to apply
#' @param edgeAlpha Transparency of edges
#' @param up_col The color to label activated gene-sets. Only applicable if `colorBy` is set to be "robustZ"
#' @param down_col The color to label inhibited gene-sets. Only applicable if `colorBy` is set to be "robustZ"
#' @param scale_edgeWidth Parameter for scaling edge width. Defaulted to 10. Higher numbers reduce all edge width
#' @param scale_nodeSize Parameter for scaling node size. Defaulted to 15. Higher numbers decreases all node sizes
#' @param nodeShape The shape to use for nodes
#' @param color_lg logical(1). Should color legend be shown
#' @param color_lg_title Title for the color legend
#' @param lb_size Size of node text labels
#' @param lb_color Color of node text labels
#' @importFrom ggraph ggraph geom_node_point geom_node_text geom_edge_link
#' @importFrom ggplot2 scale_color_continuous scale_color_manual aes_ guides aes
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' load(system.file("extdata", "gsTopology.rda", package = "sSNAPPY"))
#' load(system.file("extdata", "normalisedScores.rda", package = "sSNAPPY"))
#' # Extract pathways significantly perturbed in a sample
#' normScores_sub <- dplyr::filter(normalisedScores, sample == "E2+R5020_N2_48", adjPvalue < 0.05)
#'
#' # Color network plot nodes by robust z-score
#' plot_gsNetwork(normScores_sub, gsTopology,
#' colorBy = "robustZ", layout = "dh",
#' color_lg_title = "Robust Z-score")
#'
#' # Color network plot nodes by p-values
#' plot_gsNetwork(normScores_sub, gsTopology, layout = "dh",
#' colorBy = "pvalue", color_lg_title = "P-value")
plot_gsNetwork <- function(normalisedScores, gsTopology, colorBy = c("robustZ", "pvalue"), foldGSname = TRUE, foldafter = 2, layout = "fr",
                           edgeAlpha = 0.8,  up_col = "brown3", down_col = "steelblue3", scale_edgeWidth = 10, scale_nodeSize = 15,
                           nodeShape = 16, color_lg = TRUE, color_lg_title = NULL, lb_size = 3, lb_color = "black"){

    name <- NULL
    ## check if input has required columns
    stopifnot(colorBy %in% c("robustZ", "pvalue"))
    if (!all(c(colorBy, "gs_name") %in% colnames(normalisedScores))) stop(paste("Normalised Scores must include gs_name and", colorBy))

    # Make sure the gs topologies are a named list with at least two elements
    stopifnot(length(names(gsTopology)) == length(gsTopology))
    if(length(unique(normalisedScores$gs_name)) < 2) stop("At least 2 gene-sets are required for a network plot")

    # create igraph object
    g <- make_gsNetwork(normalisedScores, gsTopology, colorBy = colorBy, foldGSname = TRUE)

    # plot network edges
    pl <- ggraph(g, layout = layout) +
        geom_edge_link(alpha = edgeAlpha, aes_(width=~I(weight*scale_edgeWidth)), colour='darkgrey')

    if (colorBy == "robustZ"){
        pl <- pl + geom_node_point(aes_(color = ~I(color), size = ~I(size/scale_nodeSize)), shape = nodeShape,
            stroke = 0.5) +
            scale_color_manual(values = c(up_col, down_col),name = color_lg_title)
    } else(
        pl <- pl + geom_node_point(
            aes_(color = ~color, size = ~I(size/scale_nodeSize)), shape = nodeShape, stroke = 0.5) +
            scale_color_continuous(low="red", high="blue", name = color_lg_title))

    if(!color_lg){
     pl <- pl + guides(color = "none")
    }

    pl + geom_node_text(aes(label = name), size = lb_size, repel = TRUE, colour = lb_color)

}

#' Create gene-set network plot object
#'
#' @param normalisedScores A dataframe as described in the details section
#' @param gsTopology List of pathway topology matrices generated using function `weightedAdjMatrix`
#' @param colorBy Choose to color nodes either by "robustZ" or "pvalue". A column must exist in the normalisedScores for the chosen parameter
#' @param foldGSname logical(1). Should long gene-set names be folded into two lines
#' @param foldafter The number of words after which gene-set names should be folded. Defaulted to 2
#' @return igraph object
#' @importFrom reshape2 melt
#' @importFrom igraph E V graph.data.frame set_edge_attr set_vertex_attr
make_gsNetwork <- function(normalisedScores, gsTopology,  colorBy = c("robustZ", "pvalue"), foldGSname = TRUE, foldafter = 2){

    # create dummy variable to pass R CMD CHECK
    from <- to <- E <- robustZ <- NULL
    GS2Gene <- get_GSgenelist(gsTopology)
    GS2Gene <- left_join(normalisedScores, GS2Gene, by = "gs_name")

    GSlist <- split(GS2Gene[,c("gs_name", "gene_id")], f = GS2Gene$gs_name)
    nGS <- length(GSlist)
    GSname <- names(GSlist)

    w <- sapply(seq_len(nGS-1), function(x){
        sapply((x+1):nGS, function(y){
            data.frame(from = GSname[x], to = GSname[y], weight = jacIdex_func(GSlist[[x]], GSlist[[y]]))
        }, simplify = FALSE)
    }, simplify = FALSE)

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

    if(foldGSname){
        g <- set_vertex_attr(g, "name", value = vapply(V(g)$name, function(x){ifelse(length(strsplit(x, " ")[[1]]) > foldafter,
                                                                                     str_replace_nth(x, " ", "\n", foldafter),
                                                                                     x)}, character(1)))
    }
    g
}



#' Get gene-set to gene_id data frame
#'
#' @param gsTopology List of pathway topology matrices generated using function `weightedAdjMatrix`
#' @return dataframe
#' @importFrom reshape2 melt
get_GSgenelist <- function(gsTopology){
    GStoGene <- lapply(gsTopology, rownames)
    GStoGene <- reshape2::melt(GStoGene)
    colnames(GStoGene) <- c("gene_id", "gs_name")
    GStoGene
}

#' Comppute Jaccard index
#'
#' @param x first element
#' @param y second element
#' @return numeric
jacIdex_func <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}

#' Replace the nth occurrence of string with another string
#'
#' @param x string to search for pattern
#' @param pattern pattern to search for
#' @param replacement replacement string
#' @param n the nth orrcurrences of pattern that will be replaced
#' @return string
str_replace_nth <- function(x, pattern, replacement, n) {
    g <- gregexpr(pattern, x)[[1]][n]
    s <- scan(text = gsub("[()]", "", pattern),
              sep = "|",
              what = "",
              quiet = TRUE)
    substr(x, g, g) <- replacement[match(substr(x, g, g), s)]
    x
}
