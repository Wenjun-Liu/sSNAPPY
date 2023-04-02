pathwayDir <- system.file("extdata", "gsTopology.rda", package = "sSNAPPY")
load(pathwayDir)
set.seed(123)
Scores <- data.frame(
    gs_name = c("RNA degradation","Wnt signaling pathway","Histidine metabolism","Ascorbate and aldarate metabolism","Lipid and atherosclerosis"),
    robustZ = runif(5, -1, 1),
    pvalue = runif(5)) %>%
    mutate(color_Z = ifelse(robustZ < 0, "Inhibited", "Activated"))
GS <- Scores$gs_name
g_Zscore <- .make_gsNetwork(Scores, gsTopology, colorBy = "robustZ", plotIsolated = TRUE)
g_pvalue <- .make_gsNetwork(Scores, gsTopology, colorBy = "pvalue", plotIsolated = TRUE)

test_that("make_gsNetwork produces the expected outcome",{
    expect_s3_class(g_Zscore, "igraph")
    expect_equal(dim(igraph::as_data_frame(g_Zscore, "vertices")), c(5,3))
    expect_equal(dim(igraph::as_data_frame(g_Zscore, "edges")), c(10,3))
    # expect_true(is.character(igraph::V(g_Zscore)$color))
    expect_true(is.numeric(igraph::V(g_pvalue)$color))
    # expect_equal(stringr::str_subset(V(g_Zscore)$name, "Histidine"), "Histidine metabolism")
    # expect_equal(stringr::str_subset(V(g_Zscore)$name, "Ascorbate"), "Ascorbate and\naldarate metabolism")
    # g_Zscore_n3 <- .make_gsNetwork(Scores, gsTopology, colorBy = "robustZ",  plotIsolated = TRUE)
    # expect_equal(stringr::str_subset(V(g_Zscore_n3)$name, "Ascorbate"), "Ascorbate and aldarate metabolism")
})

test_that("plot_community returns error when expected", {
    expect_error(plot_community(Scores, gsTopology, colorBy = "robustZ", communityMethod = "cluster_random"), "'arg' should be one of.+")
    expect_error(plot_community(Scores, colorBy = "random"))
    expect_error(plot_community(Scores, colorBy = c("community", "robustZ")))
    gsTopology_noName <- gsTopology
    names(gsTopology_noName) <- NULL
    expect_error(plot_community(Scores, gsTopology_noName, colorBy = "robustZ"))
    expect_error(plot_community(Scores[1, ], gsTopology, colorBy = "robustZ"), "At least 2 gene-sets are required for a network plot")
    expect_error(plot_community(Scores[, -2], gsTopology, colorBy = "robustZ"), "'arg' should be one of .+")
    expect_warning(plot_community(Scores, gsTopology,  colorBy = "community", gsAnnotation = data.frame(gs_name = c(1:5), category = letters[1:5])), "Gene-set annotation does not match with topology provided. Communities won't be annotated")
})

test_that("plot_community produces the expected outcome", {
    expect_s3_class(plot_community(Scores, gsTopology, colorBy = "community"), "ggraph")
    expect_s3_class(plot_community(Scores, gsTopology, colorBy = "pvalue"), "ggraph")
    expect_s3_class(plot_community(Scores, gsTopology, colorBy = "robustZ"), "ggraph")
})
