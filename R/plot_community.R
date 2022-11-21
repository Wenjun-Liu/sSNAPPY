#' Visualize community structure in the significantly perturbed gene-set network
#'
#' @param normalisedScores A `data.frame` derived from the `normalise_by_permu()` function
#' @param gsTopology List of pathway topology matrices generated using function `retrieve_topology()`
#' @param gsAnnotation  A `data.frame` containing gene-sets categories used for pathway annotation. Must contain at least two columns:
#' `c("gs_name", "category")`, where `gs_name` denotes gene-sets names that are matched to names of pathway topology matrices, and
#' `category` records the categorization of the given pathway. If customized annotation is not provided, it's assumed that the pathways
#' investigated were from KEGG database and the inbuilt KEGG pathway annotation information will be used,
#' @param colorBy  Choose to color nodes either by *"community*, *robustZ* or *pvalue*. To color by *robustZ* or *pvalue*, a
#' column must exist in the normalisedScores for the chosen parameter
#' @param communityMethod A community Detection method supported by `igraph`. See details for all methods available.
#' @param foldGSname `logical`. Should long gene-set names be folded into two lines
#' @param foldafter The number of words after which gene-set names should be folded. Defaulted to 2
#' @param layout The layout algorithm to apply. Accept all layout supported by `igraph`. See details for all layout options available.
#' @param markCommunity `character` A *geom_mark_* method supported by `ggforce` to annotate sets of nodes belonging to the same community.
#' Either *NULL*, *ellipse*, *circle*, *hull*, *rect*
#' @param markAlpha Transparency of annotation areas. Default to 0.2
#' @param edgeAlpha Transparency of edges. Default to 0.8
#' @param up_col The color used to label activated gene-sets. Only applicable if `colorBy` is set to be "robustZ"
#' @param down_col The color used to label inhibited gene-sets. Only applicable if `colorBy` is set to be "robustZ"
#' @param scale_edgeWidth A numerical vector of length 2 to be provided to `ggraph::scale_edge_width_continuous()` for specifying
#' the minimum and maximum edge widths after transformation. Defaulted to c(0.5, 3)
#' @param scale_nodeSize A numerical vector of length 2 to be provided to `ggplot2::scale_size()` for specifying
#' the minimum and maximum node sizes after transformation. Defaulted to c(3,6)
#' @param nodeShape The shape to use for nodes
#' @param color_lg `logical` Should color legend be shown
#' @param color_lg_title Title for the color legend
#' @param edgeLegend logical` Should edge weight legend be shown
#' @param lb_size Size of node text labels
#' @param lb_color Color of node text labels
#' @param plotIsolated `logical`.Should nodes not connected to any other nodes be plotted.  Default to FALSE
#' @param ... Used to pass various potting parameters to `ggforce::geom_mark_*()`
#' @import igraph
#' @import ggforce
#' @importFrom ggnewscale new_scale_color
#' @importFrom utils data
#' @return A ggplot2 object
#' @examples
#' load(system.file("extdata", "gsTopology.rda", package = "sSNAPPY"))
#' load(system.file("extdata", "normalisedScores.rda", package = "sSNAPPY"))
#' #Subset the first 10 rows of the normalisedScores data.frame as an example
#' subset <- normalisedScores[1:15,]
#' # Color network plot nodes by the community they were assigned to and mark nodes belonging
#' # to the same community by ellipses
#' plot_community(subset, gsTopology, colorBy = "community",layout = "kk",
#' color_lg_title = "Community")
#'
#' # Color network plot nodes by pathways' directions of changes and mark nodes belonging
#' # to the same community by ellipses
#' plot_community(subset, gsTopology, colorBy = "robustZ",layout = "kk",
#' color_lg_title = "Direction of pathway Perturbation")
#'
#' # Plotting parameters accepted by `geom_mark_*` could be passed to the function
#' # to adjust the annotation area or the annotation label. See
#' # `?ggforce::geom_mark_*` for more details. For example, to change the linetype
#' # of the connector:
#' plot_community(subset, gsTopology, colorBy = "robustZ",layout = "kk",
#' color_lg_title = "Direction of pathway Perturbation", con.linetype = "dashed")
#'
#' # To change the colour and fill of `geom_mark_*` annotation, use any
#' # `scale_fill_*` and/or `scale_color_*`
#' # functions supported by `ggplot2`. For example:
#' p <- plot_community(subset, gsTopology, colorBy = "robustZ",layout = "kk",
#' markCommunity = "rect",color_lg_title = "Direction of pathway Perturbation")
#' p + ggplot2::scale_color_ordinal() + ggplot2::scale_fill_ordinal()
#' @export
plot_community <- function(normalisedScores, gsTopology, gsAnnotation = NULL, colorBy = c("robustZ", "pvalue", "community"), communityMethod = "cluster_louvain",
                           foldGSname = TRUE, foldafter = 2, layout = "fr", markCommunity = "ellipse", markAlpha = 0.2, edgeAlpha = 0.8,
                           up_col = "brown3", down_col = "steelblue3", scale_edgeWidth = c(0.5, 3), edgeLegend = FALSE, scale_nodeSize = c(3,6),
                            nodeShape = 16, color_lg = TRUE, color_lg_title = NULL, lb_size = 3, lb_color = "black", plotIsolated = FALSE, ...){

    name <- data <- community <- category <- category_n <- weight <- color <-
        size <- Community <- x <- y <- NULL
    if (!communityMethod %in% c("cluster_walktrap", "cluster_spinglass", "cluster_leading_eigen",
                                "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_label_prop",
                                "cluster_leiden", "cluster_louvain"))
        stop("CommunityMethod must be a community detection algorithm specified in the description")

    ## check if input has required columns
    stopifnot(colorBy %in% c("robustZ", "pvalue", "community") | length(colorBy) != 1)
    if (colorBy != "community" & !all(c(colorBy, "gs_name") %in% colnames(normalisedScores)))
        stop("Normalised Scores must include gs_name and column for coloring")

    # Make sure the gs topologies are a named list with at least two elements
    stopifnot(length(names(gsTopology)) == length(gsTopology))
    if(length(unique(normalisedScores$gs_name)) < 2) stop("At least 2 gene-sets are required for a network plot")

    if (is.null(gsAnnotation)){
        # if gene-set annotation info is not provided, use built-in KEGG pathways' annotations
       data(gsAnnotation_df, package = "sSNAPPY")
        } else {
            # if user provided gene-set annotation, gs_name and category column must be in the df
            if (!all(c("gs_name", "category") %in% colnames(gsAnnotation)))
                stop("gsAnnotation dataframe must include gs_name and category columns") else{
                    gsAnnotation_df = gsAnnotation
                }
        }


    if (length(intersect(names(gsTopology), gsAnnotation_df$gs_name)) == 0)
        stop("Gene-set annotation does not match with topology provided")

    # create igraph object
    g <- make_gsNetwork(normalisedScores, gsTopology, colorBy = colorBy,  plotIsolated = plotIsolated)

    # perform community detection
    commuDec_method <- get(communityMethod, envir = rlang::ns_env("igraph"))
    commuDec_result <- commuDec_method(g)
    commuDec_df <- data.frame(
        gs_name = commuDec_result$name,
        community = commuDec_result$membership
    )

    # left_join community detection results to pathway annotation
    commuDec_df <- left_join(commuDec_df, gsAnnotation_df, by = "gs_name")

    # Find the category with highest occurrence for each community
    commuDec_summary <- dplyr::group_by(commuDec_df, community, category)
    commuDec_summary <- dplyr::summarise(commuDec_summary,category_n = dplyr::n())
    commuDec_summary <- dplyr::group_by(commuDec_summary, community)
    commuDec_summary <- dplyr::filter(commuDec_summary, category_n == max(category_n))

    # if there's a tie between two categories for a given community, category names are pasted together
    commuDec_summary <- mutate(commuDec_summary, category = paste(category, collapse  = " & "))
    commuDec_summary <- left_join(commuDec_df[, c("gs_name", "community")], commuDec_summary, by = "community")
    g <- set_vertex_attr(g, "Community", index = commuDec_summary$gs_name, value = commuDec_summary$category)

    if(foldGSname){
        g <- set_vertex_attr(g, "name", value = vapply(V(g)$name, function(x){ifelse(length(strsplit(x, " ")[[1]]) > foldafter,
                                                                                     str_replace_nth(x, " ", "\n", foldafter),
                                                                                     x)}, character(1))) }

    # Find xy coordinators of nodes based on the layout chosen
    layout_method <- get(paste("layout_with_", layout, sep = ""), envir = rlang::ns_env("igraph"))
    xy <- layout_method(g)

    # plot network edges
    pl <- ggraph(g, layout = "manual", x = xy[,1], y = xy[,2]) +
        geom_edge_link(alpha = edgeAlpha, aes(width=weight, colour='darkgrey')) +
        scale_edge_width_continuous(range = scale_edgeWidth, guide = "none")

    # plot node points
    if (colorBy == "robustZ"){
        pl <- pl + geom_node_point(aes(color = color, size = size), shape = nodeShape,stroke = 0.5) +
            scale_color_manual(values = c("Activated" = up_col, "Inhibited" = down_col), name = color_lg_title) +
            scale_size(range =  scale_nodeSize, guide = "none")
    } else if (colorBy == "pvalue"){
        pl <- pl +
            geom_node_point(aes(color = color, size = size), shape = nodeShape,stroke = 0.5) +
            scale_color_continuous(low="red", high="blue", name = color_lg_title) +
            scale_size(range =  scale_nodeSize, guide = "none")
    } else {
        pl <- pl +
            geom_node_point(aes(color = Community, size = size), shape = nodeShape,stroke = 0.5) +
            scale_size(range =  scale_nodeSize, guide = "none")
    }

    if (!is.null(markCommunity)){
        mark_method <- get(paste("geom_mark_", markCommunity, sep = ""), envir = rlang::ns_env("ggforce"))
        pl <- pl +
            ggnewscale::new_scale_color() +
            mark_method(aes(x = x, y = y, fill = Community, color = Community, label = Community),
                        alpha = markAlpha, ...)
    }

    if(!color_lg){
        pl <- pl + guides(color = "none")
    }


    pl + geom_node_text(aes(label = name), size = lb_size, repel = TRUE, colour = lb_color)

}
