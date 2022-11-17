plot_community <- function(normalisedScores, gsTopology, gsAnnotation = NULL, colorBy = c("robustZ", "pvalue", "community"), communityMethod = "cluster_louvain",
                           foldGSname = TRUE, foldafter = 2, layout = "fr", markCommunity = c(NULL, "ellipse", "cirle", "hull", "rect"),edgeAlpha = 0.8,
                           up_col = "brown3", down_col = "steelblue3", scale_edgeWidth = c(0.5, 3), edgeLegend = FALSE, scale_nodeSize = c(3,6),
                            nodeShape = 16, color_lg = TRUE, color_lg_title = NULL, lb_size = 3, lb_color = "black", plotIsolated = FALSE, ...){

    name <- NULL
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

    browser()

    if (is.null(gsAnnotation)){
        # if gene-set annotation info is not provided, use built-in KEGG pathways' annotations
       data(gsAnnotation, package = "sSNAPPY")
        } else {
            # if user provided gene-set annotation, gs_name and category column must be in the df
            if (!all(c("gs_name", "catogery") %in% colnames(gsAnnotation)))
                stop("gsAnnotation dataframe must include gs_name and catogery columns") else{
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
            scale_color_manual(values = c("Activated" = up_col, "Inhibited" = down_col),name = color_lg_title) +
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
            mark_method(aes(x = x, y = y, fill = Community, label = Community), ...)
    }

    if(!color_lg){
        pl <- pl + guides(color = "none")
    }


    pl + geom_node_text(aes(label = name), size = lb_size, repel = TRUE, colour = lb_color)

}
