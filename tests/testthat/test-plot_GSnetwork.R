pathwayDir <- system.file("testdata", "test_BminsI.rda", package = "SSPT")
load(pathwayDir)
set.seed(123)
Scores <- data.frame(
    gs_name = c("RNA degradation","Wnt signaling pathway","Histidine metabolism","Ascorbate and aldarate metabolism","Lipid and atherosclerosis"),
    robustZ = runif(5, -1, 1),
    pvalue = runif(5)) %>%
    mutate(color_Z = ifelse(robustZ < 0, "Inhibited", "Activated"))
GS <- Scores$gs_name
g_Zscore <- make_gsNetwork(Scores, BminsI, colorBy = "robustZ")
g_pvalue <- make_gsNetwork(Scores, BminsI, colorBy = "pvalue")

test_that("make_gsNetwork produces the expected outcome",{
    expect_s3_class(g_Zscore, "igraph")
    expect_equal(dim(igraph::as_data_frame(g_Zscore, "vertices")), c(5,3))
    expect_equal(dim(igraph::as_data_frame(g_Zscore, "edges")), c(10,3))
    expect_true(is.character(igraph::V(g_Zscore)$color))
    expect_true(is.numeric(igraph::V(g_pvalue)$color))
    expect_equal(stringr::str_subset(V(g_Zscore)$name, "Histidine"), "Histidine metabolism")
    expect_equal(stringr::str_subset(V(g_Zscore)$name, "Ascorbate"), "Ascorbate and\naldarate metabolism")
    g_Zscore_n3 <- make_gsNetwork(Scores, BminsI, colorBy = "robustZ", foldafter = 3)
    expect_equal(stringr::str_subset(V(g_Zscore_n3)$name, "Ascorbate"), "Ascorbate and aldarate\nmetabolism")
})

test_that("plot_gsNetwork returns error when expected", {
    expect_error(plot_gsNetwork(Scores[, "gs_name"], BminsI, colorBy = "robustZ"), "Normalised Scores must include gs_name and robustZ")
    expect_error(plot_gsNetwork(Scores, colorBy = "random"))
    BminsI_noName <- BminsI
    names(BminsI_noName) <- NULL
    expect_error(plot_gsNetwork(Scores, BminsI_noName, colorBy = "robustZ"))
    expect_error(plot_gsNetwork(Scores[1, ], BminsI, colorBy = "robustZ"), "At least 2 gene-sets are required for a network plot")
})

test_that("plot_gsNetwork produces the expected outcome", {
    expect_s3_class(plot_gsNetwork(Scores, BminsI, colorBy = "robustZ"), "ggraph")
    expect_s3_class(plot_gsNetwork(Scores, BminsI, colorBy = "pvalue"), "ggraph")
    g_noLegend <- plot_gsNetwork(Scores, BminsI, colorBy = "robustZ", color_lg = FALSE)
    expect_null(cowplot::get_legend(g_noLegend))
})
