#' Visualise the community structure in significantly perturbed gene-set network
#'
#' @details A community detection strategy specified by `communityMethod` will
#' be applied to the pathway-pathway network, and communities will be annotated
#' with the pathway category that had the highest number of occurrence,
#' denoting the main biological processes perturbed in that community.
#'
#' At the moment, only KEGG pathway categories are provided with the
#' package, so if the provided `normalisedScores` contains perturbation
#' scores of pathways derived from other databases, annotation of communities
#' will not be performed unless pathway information is provided through
#' the `gsAnnotation` object. The category information needs to be
#' provided in a `data.frame` containing `gs_name` (gene-set names) and
#' `category` (categorising the given pathways).
#'
#' Plotting parameters accepted by `geom_mark_*` could be passed to the
#' function to adjust the annotation area or the annotation label. See
#' \link[ggforce]{geom_mark_ellipse} for more details.
#'
#' @param normalisedScores A `data.frame` derived from \link{normalise_by_permu}
#' @param gsTopology List of pathway topology matrices generated using
#' \link{retrieve_topology}
#' @param gsAnnotation  A `data.frame` containing gene-sets categories for
#' pathway annotation. Must contain the two columns:
#' `c("gs_name", "category")`, where `gs_name` denotes gene-sets names that are
#' matched to names of pathway topology matrices, and `category` records a
#' higher level category for each pathway. If customized annotation is not
#' provided, it will be assumed that the pathways were obtained from the KEGG
#' database and inbuilt KEGG pathway annotation information will be used
#' @param colorBy  Can be any column with in the `normalisedScores` object, or
#' the additional value "community".
#' @param communityMethod A community detection method supported by `igraph`.
#' See details for all methods available.
#' @param foldGSname `logical`. Should long gene-set names be folded into two
#' lines
#' @param foldafter The number of words after which gene-set names should be
#' folded. Defaults to 2
#' @param labelFun function to manipulate or modify gene-set labels. By default,
#' any database will be stripped from the prefix using a regex pattern
#' @param layout The layout algorithm to apply. Accepted layouts are
#' `"fr", "dh", "gem", "graphopt", "kk", "lgl", "mds" and "sugiyama"`
#' @param markCommunity `character` A `geom_mark_*` method supported by
#' `ggforce` to annotate sets of nodes belonging to the same community.
#' Either `*NULL*, *ellipse*, *circle*, *hull*, *rect*`
#' @param markAlpha Transparency of annotation areas.
#' @param edgeAlpha Transparency of edges.
#' @param scale_edgeWidth A numerical vector of length 2 to be provided to
#' `ggraph::scale_edge_width_continuous()` for specifying the minimum and
#' maximum edge widths after transformation.
#' @param scale_nodeSize A numerical vector of length 2 to be provided to
#' `ggplot2::scale_size()` for specifying
#' the minimum and maximum node sizes after transformation.
#' @param nodeShape The shape to use for nodes
#' @param color_lg_title Title for the color legend
#' @param edgeLegend `logical` Should edge weight legend be shown
#' @param lb_size Size of node text labels
#' @param lb_color Color of node text labels
#' @param plotIsolated `logical(1)` Should nodes not connected to any other
#' nodes be plotted. Defaults to FALSE
#' @param ... Used to pass various potting parameters to `ggforce::geom_mark_*()`
#'
#' @return A ggplot2 object
#' @examples
#' load(system.file("extdata", "gsTopology.rda", package = "sSNAPPY"))
#' load(system.file("extdata", "normalisedScores.rda", package = "sSNAPPY"))
#' #Subset the first 10 rows of the normalisedScores data.frame as an example
#' subset <- normalisedScores[1:15,]
#' subset$status <- ifelse(subset$robustZ > 0, "Activated", "Inhibited")
#' # Color network plot nodes by the community they were assigned to and mark
#' # nodes belonging to the same community by ellipses
#' plot_community(subset, gsTopology, colorBy = "community",layout = "kk",
#' color_lg_title = "Community")
#'
#' # Color network plot nodes by pathways' directions of changes and mark nodes
#' # belonging to the same community by ellipses
#' plot_community(subset, gsTopology, colorBy = "status",layout = "kk",
#' color_lg_title = "Direction of pathway perturbation")
#'
#' # To change the colour and fill of `geom_mark_*` annotation, use any
#' # `scale_fill_*` and/or `scale_color_*`
#' # functions supported by `ggplot2`. For example:
#' p <- plot_community(subset, gsTopology, colorBy = "status",layout = "kk",
#' markCommunity = "rect",color_lg_title = "Direction of pathway perturbation")
#' p + ggplot2::scale_color_ordinal() + ggplot2::scale_fill_ordinal()
#'
#' @import igraph
#' @import ggforce
#' @import ggplot2
#' @import ggraph
#' @importFrom utils data
#' @importFrom dplyr left_join group_by summarise mutate ungroup distinct
#' @importFrom tidyr drop_na
#'
#' @export
plot_community <- function(
        normalisedScores, gsTopology, gsAnnotation = NULL, colorBy = "community",
        communityMethod = c(
            "louvain", "walktrap", "spinglass", "leading_eigen",
            "edge_betweenness", "fast_greedy", "label_prop", "leiden"
        ),
        foldGSname = TRUE, foldafter = 2, labelFun = .rm_prefix,
        layout = c(
            "fr", "dh", "gem", "graphopt", "kk", "lgl", "mds", "sugiyama"
        ),
        markCommunity = "ellipse", markAlpha = 0.2, color_lg_title = NULL,
        edgeAlpha = 0.8, scale_edgeWidth = c(0.5, 3), edgeLegend = FALSE,
        scale_nodeSize = c(3,6), nodeShape = 16,  lb_size = 3,
        lb_color = "black", plotIsolated = FALSE, ...
){

    name <- data <- community <- category <- category_n <- weight <- color <-
        size <- Community <- x <- y <- NULL
    communityMethod <- match.arg(communityMethod)
    communityMethod <- paste0("cluster_", communityMethod)
    layout <- match.arg(layout)

    ## check if input has required columns
    cols <- colnames(normalisedScores)
    colorBy <- match.arg(colorBy, c(cols, "community"))
    stopifnot("gs_name" %in% cols)

    # Make sure the gs topologies are a named list with at least two elements
    stopifnot(length(names(gsTopology)) == length(gsTopology))
    if (length(unique(normalisedScores$gs_name)) < 2)
        stop("At least 2 gene-sets are required for a network plot")
    gsTopology <- gsTopology[names(gsTopology) %in% normalisedScores$gs_name]

    if (is.null(gsAnnotation)) {
        ## if gene-set annotation info is not provided, use built-in KEGG
        ## pathway annotations
        gsAnnotation_df <- c()
        data("gsAnnotation_df", envir = environment())
        gsAnnotation <- gsAnnotation_df
    } else {
        ## if user provided gene-set annotation, gs_name and category column
        ## must be in the df
        if (!all(c("gs_name", "category") %in% colnames(gsAnnotation)))
            stop("gsAnnotation must include gs_name and category columns")
    }
    # create igraph object
    if (colorBy != "community") {
        g <- .make_gsNetwork(
            normalisedScores, gsTopology, colorBy = colorBy,
            plotIsolated = plotIsolated, labelFun = NULL
        )
    } else {
        g <- .make_gsNetwork(
            normalisedScores, gsTopology, colorBy = NULL,
            plotIsolated = plotIsolated, NULL
        )
    }

    # perform community detection
    comm_method <- get(communityMethod, envir = rlang::ns_env("igraph"))
    comm_result <- comm_method(g)
    comm_df <- data.frame(
        gs_name = comm_result$names, community = comm_result$membership
    )

    if (length(intersect(names(gsTopology), gsAnnotation$gs_name)) == 0) {
        warning(
            "Gene-set annotation does not match with topology provided. ",
            "Communities won't be annotated"
        )
        g <- set_vertex_attr(
            g, "Community", index = comm_df$gs_name,
            value = as.character(comm_df$community)
        )
        silent_commu <- TRUE
    } else {
        ## left_join community detection results to pathway annotation
        comm_df <- left_join(comm_df, gsAnnotation, by = "gs_name")

        ## Find the category with highest occurrence for each community
        comm_summary <- group_by(comm_df, community, category)
        comm_summary <- summarise(
            comm_summary, category_n = dplyr::n(), .groups = "keep"
        )
        comm_summary <- group_by(comm_summary, community)
        comm_summary <- dplyr::filter(
            comm_summary, category_n == max(category_n)
        )
        comm_summary <- drop_na(comm_summary)

        ## if there's a tie between two categories for a given community,
        ## category names are pasted together
        comm_summary <- mutate(
            comm_summary, category = paste(category, collapse  = " &\n")
        )
        comm_summary <- ungroup(comm_summary)
        comm_summary <- distinct(comm_summary, community, category)
        comm_summary <- left_join(
            comm_df[, c("gs_name", "community")], comm_summary, by = "community"
        )
        g <- set_vertex_attr(
            g, "Community", index = comm_summary$gs_name,
            value = comm_summary$category
        )
        silent_commu <- FALSE
    }

    if (foldGSname) {
        nm <-  .str_replace_nth(
            V(g)$name, pattern = " ", replacement = "\n", n = foldafter
        )
        g <- set_vertex_attr(g, "name", value = as.character(nm) )
    }

    # Find xy coordinators of nodes based on the layout chosen
    layout_method <- get(
        paste("layout_with_", layout, sep = ""),
        envir = rlang::ns_env("igraph")
    )
    xy <- layout_method(g)

    ## Tidy up node labels if a function has been passed to this argument
    if (!is.null(labelFun)) {
        ## This allows for stadard label replacement, but also for users to
        ## provide more complex methods of string manipulation
        stopifnot(is(labelFun, "function"))
        nm <- vertex_attr(g, "name")
        new_nm <- labelFun(nm)
        stopifnot(length(nm) == length(new_nm))
        g <- set_vertex_attr(g, "name", seq_along(nm), new_nm)
    }

    # plot network edges
    pl <- ggraph(g, layout = "manual", x = xy[,1], y = xy[,2]) +
        geom_edge_link(
            alpha = edgeAlpha, aes(width = weight), colour = 'darkgrey'
        ) +
        scale_edge_width_continuous(range = scale_edgeWidth, guide = "none")

    # plot node points
    if (!is.null(colorBy)) {
        color <- ifelse(colorBy == "community", "Community", "color")
        pl <- pl + geom_node_point(
            aes(color = !!sym(color), size = size), shape = nodeShape,
            stroke = 0.5
        ) +
            scale_size(range =  scale_nodeSize, guide = "none")  +
            labs(colour = color_lg_title)
    } else {
        pl <- pl +
            geom_node_point(aes(size = size), shape = nodeShape,stroke = 0.5) +
            scale_size(range =  scale_nodeSize, guide = "none")
    }

    if (!is.null(markCommunity)) {
        mark_method <- get(
            paste("geom_mark_", markCommunity, sep = ""),
            envir = rlang::ns_env("ggforce")
        )
        if (silent_commu) {
            pl <- pl + mark_method(
                aes(x = x, y = y, fill = Community), alpha = markAlpha, ...
            )
        } else {
            pl <- pl + mark_method(
                aes(x = x, y = y, fill = Community, label = Community),
                color = NA, alpha = markAlpha, ...
            )
        }
    }

    pl + geom_node_text(
        aes(label = name), size = lb_size, repel = TRUE, colour = lb_color
    )

}
