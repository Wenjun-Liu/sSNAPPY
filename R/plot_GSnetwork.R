#' @title Plot significantly perturbed gene-sets as a network
#'
#' @param normalisedScores A `data.frame` of pathway perturbation scores derived
#' from the `normalise_by_permu()` function
#' @param gsTopology List of pathway topology matrices generated using function 
#' `retrieve_topology()`
#' @param colorBy Choose to color nodes either by *robustZ* or *pvalue*. A 
#' column must exist in the `normalisedScores` `data.frame` for the chosen 
#' parameter
#' @param foldGSname logical(1) Should long gene-set names fold across multiple
#' lines
#' @param foldafter The number of words after which gene-set names should be 
#' folded. Defaults to 2
#' @param layout The layout algorithm to apply. Accept all layout supported by 
#' `igraph`
#' @param edgeAlpha numeric(1) Transparency of edges. Default to 0.8
#' @param scale_edgeWidth A numerical vector of length 2 to be provided to 
#' `ggraph::scale_edge_width_continuous()` for specifying
#' the minimum and maximum edge widths after transformation. Defaults to 
#' `c(0.5, 3)`
#' @param edgeLegend logical(1) Should edge weight legend be shown
#' @param scale_nodeSize A numeric vector of length 2 to be provided to 
#' `ggplot2::scale_size()` for specifying the minimum and maximum node sizes 
#' after transformation. Defaulted to `c(3,6)`
#' @param nodeShape Shape of nodes
#' @param color_lg logical(1) Should color legend be shown
#' @param color_lg_title Optional title for the color legend
#' @param lb_size Size of node text labels
#' @param lb_color Color of node text labels
#' @param plotIsolated logical(1) Should nodes not connected to any other nodes 
#' be plotted.  Defaults to FALSE
#' @param max_overlaps passed to \link[ggraph]{geom_node_text}
#' @param ... Not used
#' 
#' @return A ggplot2 object
#'
#' @examples
#' load(system.file("extdata", "gsTopology.rda", package = "sSNAPPY"))
#' load(system.file("extdata", "normalisedScores.rda", package = "sSNAPPY"))
#'
#' # Subset pathways significantly perturbed in sample R5020_N2_48
#' subset <- dplyr::filter(
#'   normalisedScores, adjPvalue < 0.05, sample == "R5020_N2_48"
#' )
#' subset[["status"]] <- ifelse(subset[["robustZ"]]> 0, "Activated", "Inhibited")
#' 
#'
#' # Color network plot nodes by status
#' plot_gs_network(subset, gsTopology,
#' colorBy = "status", layout = "dh",
#' color_lg_title = "Direction of pathway Perturbation")
#'
#' # Color network plot nodes by p-values
#' plot_gs_network(subset, gsTopology, layout = "dh",
#' colorBy = "pvalue", color_lg_title = "P-value")
#' 
#' @importFrom ggraph ggraph geom_node_point geom_node_text geom_edge_link
#' @import ggplot2 
#' @export
plot_gs_network <- function(
        normalisedScores, gsTopology, colorBy = NULL, 
        foldGSname = TRUE, foldafter = 2, layout = "fr", edgeAlpha = 0.8,  
        scale_edgeWidth = c(0.5, 3), edgeLegend = FALSE, scale_nodeSize = c(3,6),
        nodeShape = 16, color_lg = TRUE, color_lg_title = NULL, lb_size = 3, 
        lb_color = "black", plotIsolated = FALSE, max_overlaps = 10, ...
){
    
    name <- weight <- color <- size <- NULL
    ## check if input has required columns
    cols <- colnames(normalisedScores)
    if (!is.null(colorBy)) colorBy <- match.arg(colorBy, cols)
    
    if (!"gs_name" %in% cols) stop("Normalised Scores must include gs_name")
    
    # Make sure the gs topologies are a named list with at least two elements
    stopifnot(length(names(gsTopology)) == length(gsTopology))
    if(length(unique(normalisedScores$gs_name)) < 2) 
        stop("At least 2 gene-sets are required for a network plot")
    gsTopology <- gsTopology[names(gsTopology) %in% normalisedScores$gs_name]
    
    # create igraph object
    g <- make_gsNetwork(
        normalisedScores, gsTopology, colorBy = colorBy, 
        plotIsolated = plotIsolated
    )
    
    if (foldGSname){
        nm <-  str_replace_nth(
            V(g)$name, pattern = " ", replacement = "\n", n = foldafter
        )
        g <- set_vertex_attr(g, "name", value = nm)
    }
    
    # plot network edges
    ## Given the deprecation of width from line parameters in ggplot2 3.4.0,
    ## this raises an error which needs to be resolved. Hoping for resolution
    ## but ggraph has not been updated since Sep-2022
    pl <- ggraph(g, layout = layout) +
        geom_edge_link(
            aes(edge_width = weight), alpha = edgeAlpha, colour = 'darkgrey'
        ) +
        scale_edge_width_continuous(range = scale_edgeWidth, guide = "none")
    
    
    if (!is.null(colorBy)) {
        pl <- pl + 
            geom_node_point(
                aes(color = color, size = size), shape = nodeShape, stroke = 0.5
            ) +
            scale_size(range =  scale_nodeSize, guide = "none") +
            labs(colour = color_lg_title)
        if(!color_lg) pl <- pl + guides(color = "none")
    } else {
        pl <- pl +
            geom_node_point(aes(size = size), shape = nodeShape,stroke = 0.5) +
            scale_size(range =  scale_nodeSize, guide = "none")
    }
    
    pl + 
        geom_node_text(
            aes(label = name), size = lb_size, repel = TRUE, colour = lb_color,
            max.overlaps = max_overlaps
        ) 
    
}


#' @importFrom reshape2 melt
#' @importFrom dplyr bind_rows
#' @import igraph 
#' @keywords internal
make_gsNetwork <- function(
        normalisedScores, gsTopology, colorBy = NULL, plotIsolated
){
    
    # create dummy variable to pass R CMD CHECK
    from <- to <- NULL
    GS2Gene <- get_GSgenelist(gsTopology)
    GS2Gene <- left_join(
        normalisedScores, GS2Gene, by = "gs_name", multiple = "all"
        ## Setting multiple = "all" requires dplyr > 1.1.0
    )
    
    GSlist <- split(GS2Gene[,c("gs_name", "entrezid")], f = GS2Gene$gs_name)
    nGS <- length(GSlist)
    GSname <- names(GSlist)
    
    w <- lapply(
        seq_len(nGS - 1), 
        function(x){
            lapply(
                seq(x + 1, nGS, by = 1),
                function(y) {
                    w <- jacIdex_func(GSlist[[x]]$entrezid, GSlist[[y]]$entrezid)
                    data.frame(from = GSname[x], to = GSname[y], weight = w)
                }
            )
        }
    )
    
    w <- bind_rows(lapply(w, bind_rows))
    w <- dplyr::filter(w, from != to)
    
    g <- graph.data.frame(dplyr::select(w, from, to), directed = FALSE)
    g <- set_edge_attr(g, "weight", value = w$weight)
    
    GSsize <- melt(lapply(GSlist, nrow))
    colnames(GSsize) <- c("size", "from")
    g <- set_vertex_attr(g, "size", index = GSsize$from, value = GSsize$size)
    
    if (!is.null(colorBy)) {
        GSvalue <- unique(GS2Gene[,c("gs_name", colorBy)])
        g <- set_vertex_attr(
            g, "color", index = GSvalue$gs_name, value = GSvalue[[colorBy]]
        )
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
#' @keywords internal
get_GSgenelist <- function(gsTopology, mapEntrezID = NULL){
    GStoGene <- lapply(gsTopology, rownames)
    GStoGene <- reshape2::melt(GStoGene)
    colnames(GStoGene) <- c("entrezid", "gs_name")
    
    if (all(c("entrezid","mapTo") %in% colnames(mapEntrezID))){
        if (any(GStoGene$entrezid %in% mapEntrezID$entrezid)){
            left_join(
                GStoGene, mapEntrezID[,c("entrezid","mapTo")], by = "entrezid",
                multiple = "all"
            )
        } else {
            warning("None of the Entrez IDs in mapEntrezID mapped to gsTopology.")
            return(GStoGene)
        }
    } else {
        return(GStoGene)
    }
}

#' @keywords internal
jacIdex_func <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    cmn <- intersect(x, y)
    universe <- union(x, y)
    length(cmn) / length(universe)
}

#' @importFrom stringr str_split
#' @keywords internal
str_replace_nth <- function(x, pattern, replacement, n) {
    x_list <- str_split(x, pattern = pattern)
    fold <- vapply(x_list, length, integer(1)) > n
    x_list[fold] <- lapply(
        x_list[fold],
        function(x){
            index <- seq(n, length(x), by = n)
            x[index] <- paste(x[index], replacement, sep = "")
            not_index <- setdiff(seq_along(x), index)
            x[not_index] <- paste(x[not_index], pattern, sep = "")
            x <- paste(x, collapse = "")
        }
    )
    x_list[!fold] <- lapply(x_list[!fold], paste, collapse = pattern)
    unlist(x_list)
}
