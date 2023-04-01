pathwayDir <- system.file("extdata", "gsTopology.rda", package = "sSNAPPY")
load(pathwayDir)
set.seed(123)
Scores <- data.frame(
    gs_name = c("RNA degradation","Wnt signaling pathway","Histidine metabolism","Ascorbate and aldarate metabolism","Lipid and atherosclerosis"),
    robustZ = runif(5, -1, 1),
    pvalue = runif(5))
Scores <- dplyr::mutate(Scores, color_Z = ifelse(robustZ < 0, "Inhibited", "Activated"))
GS <- Scores$gs_name
g_Zscore <- make_gsNetwork(Scores, gsTopology, colorBy = "robustZ", plotIsolated = TRUE)
g_pvalue <- make_gsNetwork(Scores, gsTopology, colorBy = "pvalue", plotIsolated = TRUE)
entrez2name_dir <- system.file("extdata", "entrez2name.rda", package = "sSNAPPY")
load(entrez2name_dir)

test_that("get_GSgenelist produces the expected outcome", {
    expect_true(setequal(colnames(get_GSgenelist(gsTopology)), c("entrezid", "gs_name")))
    expect_true(setequal(colnames(get_GSgenelist(gsTopology, mapEntrezID = "random")), c("entrezid", "gs_name")))
    expect_true(setequal(colnames(get_GSgenelist(gsTopology[1:5], mapEntrezID = entrez2name)), c("entrezid", "gs_name", "mapTo")))
})

test_that("make_gsNetwork produces the expected outcome",{
    expect_s3_class(g_Zscore, "igraph")
    expect_equal(dim(igraph::as_data_frame(g_Zscore, "vertices")), c(5,3))
    expect_equal(dim(igraph::as_data_frame(g_Zscore, "edges")), c(10,3))
    # expect_true(is.character(igraph::V(g_Zscore)$color))
    expect_true(is.numeric(igraph::V(g_pvalue)$color))
    # expect_equal(stringr::str_subset(V(g_Zscore)$name, "Histidine"), "Histidine metabolism")
    # expect_equal(stringr::str_subset(V(g_Zscore)$name, "Ascorbate"), "Ascorbate and\naldarate metabolism")
    g_Zscore_n3 <- make_gsNetwork(Scores, gsTopology, colorBy = "robustZ",  plotIsolated = TRUE)
    # expect_equal(stringr::str_subset(V(g_Zscore_n3)$name, "Ascorbate"), "Ascorbate and aldarate metabolism")
})

test_that("plot_gs_network returns error when expected", {
    expect_error(plot_gs_network(Scores[, "gs_name"], gsTopology, colorBy = "robustZ"), "'arg' should be one of")
    expect_error(plot_gs_network(Scores, colorBy = "random"))
    gsTopology_noName <- gsTopology
    names(gsTopology_noName) <- NULL
    expect_error(plot_gs_network(Scores, gsTopology_noName, colorBy = "robustZ"))
    expect_error(plot_gs_network(Scores[1, ], gsTopology, colorBy = "robustZ"), "At least 2 gene-sets are required for a network plot")
})

test_that("plot_gs_network produces the expected outcome", {
    expect_s3_class(plot_gs_network(Scores, gsTopology, colorBy = "robustZ"), "ggraph")
    expect_s3_class(plot_gs_network(Scores, gsTopology, colorBy = "pvalue"), "ggraph")
    g_noLegend <- plot_gs_network(Scores, gsTopology, colorBy = "robustZ", color_lg = FALSE)
    expect_null(cowplot::get_legend(g_noLegend))
})
