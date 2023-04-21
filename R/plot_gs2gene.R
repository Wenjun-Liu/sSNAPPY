#' Plot pathways and genes contained in them as a network
#'
#' @details Taking the perturbation scores of a list of gene-sets derived from
#' `normalise_by_permu()` as input, this function matches gene-sets to
#' their associated genes by utilizing information from pathway topology matrices.
#'
#' If providing logFC values as a named vector, the names must be entrezgene IDs
#' in the format of "ENTREZID:XXXX" for compatibility with the values returned
#' by `retrieve_topology()`. If not providing this vector, only genes associated
#' with two or more pathways will be added to the plot, however, it should be
#' noted that if omitting this vector, network plots can easily become
#' unmanageable.
#'
#' Users can also choose to provide a `mapEntrezID` `data.frame` to match
#' entrezgene IDs to their chosen identifiers. The `data.frame` should contain
#' the columns: `"entrezid"` and `"mapTo"`.
#'
#' If `geneFC` is provided, gene nodes will be colored by values provided,
#' otherwise all gene nodes will drawn in grey.
#'
#' Since some gene-sets could can contain hundreds of genes, it is not
#' recommended to plot all genes. If `mapEntrezID` `data.frame` is provided,
#' only genes included in that `data.frame` will be used in the plot.
#'
#' It is strongly recommended to filter genes using some criteria, such as those
#' with the largest magnitude of change. If all pathway genes are desired,
#' please consider setting `labelGene` to FALSE to remove gene names.
#'
#' @param normalisedScores A `data.frame` derived from the
#' `normalise_by_permu()` function. Only gene-sets of interest should be included
#' @param gsTopology List of pathway topology matrices generated using function
#' `retrieve_topology()`
#' @param geneFC An optional named vector of pathways' fold changes
#' @param mapEntrezID Optional. A `data.frame` matching genes' entrez IDs to
#' another identifier with preferred labels. Must contain the columns:
#' `"entrezid"` and `"mapTo"`
#' @param colorGsBy Column within `normalisedScores` to color gene-set/pathway
#' nodes by
#' @param foldGSname `logical`. Should long gene-set names be folded into two
#' lines
#' @param foldafter The number of words after which gene-set names should be
#' folded.
#' @param labelFun function to manipulate or modify gene-set labels. By default,
#' any database will be stripped from the prefix using a regex pattern
#' @param filterGeneBy Filtration cut-off applied to genes' connectivity (ie.
#' how many pathways was a gene involved in).
#' @param layout The layout algorithm to apply. Accepts all layout supported by
#' `igraph`.
#' @param edgeColor,edgeAlpha Color and transparency of edges
#' @param edgeArc The bend of edges. 1 approximates a semi-circle whilst 0
#' will give a straight line.
#' @param geneNodeSize,geneNodeShape Size and shape for gene nodes
#' @param geneNameSize,geneNameColor,geneNameFace Size, color and fontface
#' to use for gene labels
#' @param labelGene `logical(1)` Should the gene names be included
#' @param gsNodeSize Size for gene-set/pathway nodes
#' @param gsNodeShape Shape for gene-set/pathway nodes. Should be a shape with
#' a fill parameter, such as 21:25
#' @param gsNodeStroke,gsNodeOutline Border thickness and color for
#' gene-set/pathway nodes
#' @param gsNameSize,gsNameColor Size and color of gene-set/pathway labels
#' @param geneLegTitle `character(1)`. Legend title for gene nodes
#' @param gsLegTitle `character(1)` Legend title for gene-set/pathway nodes
#' @param maxOverlaps passed to \link[ggraph]{geom_node_text}
#' @param ... Not used
#'
#' @return A ggplot2 object
#'
#' @examples
#'
#' load(system.file("extdata", "gsTopology.rda", package = "sSNAPPY"))
#' load(system.file("extdata", "normalisedScores.rda", package = "sSNAPPY"))
#'
#' # Subset pathways significantly perturbed in sample R5020_N2_48
#' subset <- dplyr::filter(normalisedScores, adjPvalue < 0.05, sample == "R5020_N2_48")
#' subset$response <- ifelse(subset$robustZ > 0, "Activated", "Inhibited")
#'
#' # Color gene-sets nodes by robust z-scores.
#' plot_gs2gene(
#'   subset, gsTopology, colorGsBy = "robustZ", labelGene = FALSE, geneNodeSize = 1,
#'   gsNodeSize = 4
#' ) + scale_fill_gradient2()
#' # When fold-changes are not provided, gene nodes are colored grey.
#'
#' # To color genes by their direction of change, firstly compute single-sample logFC
#' data(logCPM_example)
#' data(metadata_example)
#' metadata_example <- dplyr::mutate(metadata_example, treatment = factor(
#'    treatment, levels = c("Vehicle", "E2+R5020", "R5020")))
#' ls <- weight_ss_fc(
#'   logCPM_example, metadata = metadata_example,
#'   groupBy = "patient", treatColumn = "treatment",
#'   sampleColumn = "sample"
#' )
#' # Provide fold-changes of sample R5020_N2_48 as a named vector
#' plot_gs2gene(
#'   subset, gsTopology, geneFC = ls$logFC[,"R5020_N2_48"],
#'   colorGsBy = "response", labelGene = FALSE
#' ) + scale_colour_gradient2()
#'
#' # By default, the function only include genes involved in at least 2 pathways,
#' # which can be overwritten by the `filterGeneBy` parameter. But there are still
#' # a large number of genes, making the plot cumbersome. Instead, only include
#' # fold-changes of genes within the top 500 absolute values for fold-change
#' top500 <- rank(1/abs(ls$logFC[,"R5020_N2_48"])) <= 500
#' fcByDir <- ifelse(ls$logFC[top500,"R5020_N2_48"] > 0, "Up-Regulated", "Down-Regulated")
#' plot_gs2gene(subset, gsTopology, geneFC = fcByDir, colorGsBy = "response") +
#'   scale_fill_manual(values = c("darkred", "lightskyblue")) +
#'   scale_colour_manual(values = c("red", "blue"))
#'
#' # To make the gene labels more informative, map genes' entrez id to chosen identifiers.
#' load(system.file("extdata", "entrez2name.rda", package = "sSNAPPY"))
#' plot_gs2gene(
#'   subset, gsTopology, geneFC = fcByDir, mapEntrezID = entrez2name,
#'   colorGsBy = "response", gsNodeSize = 4
#' ) +
#'   scale_fill_manual(values = c("darkred", "lightskyblue"), name = "Pathway") +
#'   scale_colour_manual(values = c("blue", "red"), name = "Gene\nDirection")
#'
#' @importFrom ggraph geom_node_label geom_node_text geom_node_point
#' @importFrom ggraph geom_edge_arc
#' @import ggplot2
#' @export
plot_gs2gene <- function(
        normalisedScores, gsTopology, geneFC = NULL, mapEntrezID = NULL,
        colorGsBy = NULL, foldGSname = TRUE, foldafter = 2,
        labelFun = .rm_prefix, filterGeneBy = 2,
        layout = c(
            "fr", "dh", "gem", "graphopt", "kk", "lgl", "mds", "sugiyama"
        ),
        edgeColor = "darkgrey", edgeAlpha = 0.8, edgeArc = 0.5,
        geneNodeSize = 3, geneNodeShape = 17,
        geneNameFace = c("italic", "plain", "bold", "bold-italic"),
        geneNameColor = "grey30", geneNameSize = 3, labelGene = TRUE,
        gsNodeSize = 2, gsNodeShape = 21, gsNodeStroke = 0.5,
        gsNodeOutline = "white", gsNameSize = 6, gsNameColor = "black",
        geneLegTitle = "Mean logFC", gsLegTitle = colorGsBy,
        maxOverlaps = 10, ...
){

    name <- weight <- color <- size <- fill <- type <- NULL
    ## check if input has required columns
    cols <- colnames(normalisedScores)
    if (!is.null(colorGsBy)) colorGsBy <- match.arg(colorGsBy, cols)
    if (!"gs_name" %in% cols) stop("Normalised Scores must include gs_name")
    layout <- match.arg(layout)
    geneNameFace <- match.arg(geneNameFace)
    if (!gsNodeShape %in% 21:25)
        stop("Nodes require a fill parameter and can only be shapes 21:25")

    # Make sure the gs topologies are a named list with at least two elements
    stopifnot(length(names(gsTopology)) == length(gsTopology))
    if (length(unique(normalisedScores$gs_name)) < 2)
        stop("At least 2 gene-sets are required for a network plot")
    gsTopology <- gsTopology[names(gsTopology) %in% normalisedScores$gs_name]

    # create igraph object
    g <- .make_gs2gene_network(
        normalisedScores, gsTopology, geneFC, colorGsBy, mapEntrezID
    )

    # filter genes by connectivity
    if (filterGeneBy >= 2) 
        g <- delete_vertices(g, degree(g) < filterGeneBy & V(g)$type == "GENE")
    

    if (foldGSname) {
        nm <-  .str_replace_nth(
            V(g)$name, pattern = " ", replacement = "\n", n = foldafter
        )
        g <- set_vertex_attr(g, "name", value = nm)
    }
    
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

    pl <- ggraph(g, layout = layout) +
        geom_edge_arc(
            alpha = edgeAlpha, colour = edgeColor, strength = edgeArc
        )


    ## Add gene-set nodes
    if (!is.null(colorGsBy)) {
        pl <- pl + geom_node_point(
            aes(fill = fill), data = dplyr::filter(pl$data, type == "GS"),
            colour = gsNodeOutline, stroke = gsNodeStroke, size = gsNodeSize,
            shape = gsNodeShape
        )
    } else {
        pl <- pl + geom_node_point(
            data = dplyr::filter(pl$data, type == "GS"), fill = "grey50",
            colour = gsNodeOutline, stroke = gsNodeStroke, size = gsNodeSize,
            shape = gsNodeShape
        )
    }

    ## Add gene nodes
    if (!all(is.na(V(g)$color))) {
        pl <- pl + geom_node_point(
            aes(color = color), data = dplyr::filter(pl$data, type == "GENE"),
            size = geneNodeSize, shape = geneNodeShape
        )
    } else {
        pl <- pl + geom_node_point(
            data = dplyr::filter(pl$data, type == "GENE"),
            size = geneNodeSize, shape = geneNodeShape, colour = "grey50"
        )
    }

    if (labelGene) {
        pl <- pl + geom_node_text(
            aes(label = name),
            data = dplyr::filter(pl$data, type == "GENE"),
            size = geneNameSize, repel = TRUE, show.legend = FALSE,
            fontface = geneNameFace, colour = geneNameColor,
            max.overlaps = maxOverlaps
        )
    }

    ## Add gene-set labels as the final step
    pl <- pl + geom_node_text(
        aes(label = name), data = dplyr::filter(pl$data, type == "GS"),
        color = gsNameColor, size = gsNameSize, repel = TRUE,
        show.legend = FALSE, max.overlaps = maxOverlaps
    ) +
        labs(fill = gsLegTitle, colour = geneLegTitle)

}

#' @importFrom reshape2 melt
#' @importFrom igraph E V graph.data.frame set_edge_attr set_vertex_attr degree
#' @importFrom igraph delete_vertices delete_edges
#' @importFrom dplyr left_join
#' @importFrom rlang "!!" sym
#' @importFrom stats setNames
#' @keywords internal
.make_gs2gene_network <- function(
        normalisedScores, gsTopology, geneFC, colorGsBy, mapEntrezID
){

    # create dummy variable to pass R CMD CHECK
    from <- to <- entrezid <- gs_name <- NULL
    GS2Gene <- .get_GSgenelist(gsTopology, mapEntrezID)
    GS2Gene <- left_join(
        normalisedScores, GS2Gene, by = "gs_name", multiple = "all"
    )
    GS2Gene <- unique(GS2Gene)
    id_col <- "entrezid"
    if (!is.null(mapEntrezID)) {
        id_col <- "mapTo"
        GS2Gene <- dplyr::filter(GS2Gene, !is.na(!!sym(id_col)))
    }

    # set color for gene nodes
    if (
        !is.vector(geneFC) | is.null(names(geneFC)) |
        length(intersect(GS2Gene$entrezid, names(geneFC))) == 0
    ) {
        message(
            "Gene fold-changes were not provided as a named vector. ",
            "All genes will be colored identically."
        )
        gene_col <- rep(NA, times = nrow(GS2Gene))
        gene_col <- setNames(gene_col, GS2Gene[[id_col]])
        # ## In this case, maybe we should restrict the output to genes that
        # ## only appear in more than one pathway. The next 3 lines are experimental...
        # dups <- names(which(table(names(gene_col)) > 1))
        # GS2Gene <- dplyr::filter(GS2Gene, !!sym(id_col) %in% dups)
        # gene_col <- gene_col[dups]
    } else {
        GS2Gene <- dplyr::filter(GS2Gene, entrezid %in% names(geneFC))
        geneFC <- geneFC[names(geneFC) %in% GS2Gene$entrezid]
        id2Name <- setNames(GS2Gene[[id_col]], GS2Gene$entrezid)
        gene_col <- setNames(geneFC, as.character(id2Name[names(geneFC)]))
    }

    g <- graph.data.frame(
        dplyr::select(GS2Gene, gs_name, !!sym(id_col)), directed = FALSE
    )

    # set color for gene nodes
    g <- set_vertex_attr(g, "color", index = names(gene_col), value = gene_col)

    ## Fill for Gene-Set nodes
    if (!is.null(colorGsBy)) {
        GSvalue <- unique(GS2Gene[,c("gs_name", colorGsBy)])
        g <- set_vertex_attr(
            g, "fill", index = GSvalue$gs_name, value = GSvalue[[colorGsBy]]
        )
    } else {
        g <- set_vertex_attr(
            g, "fill", index = unique(GS2Gene$gs_name), value = NA
        )
    }

    ## set types of nodes & return
    nodeType <- ifelse(V(g)$name %in% GS2Gene$gs_name, "GS", "GENE")
    g <- set_vertex_attr(g, "type", value = nodeType)
    g

}
