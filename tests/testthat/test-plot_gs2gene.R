pathwayDir <- system.file("extdata", "gsTopology.rda", package = "sSNAPPY")
load(pathwayDir)
set.seed(123)
Scores <- data.frame(
    gs_name = c("RNA degradation","Wnt signaling pathway","Histidine metabolism","Ascorbate and aldarate metabolism","Lipid and atherosclerosis"),
    robustZ = runif(5, -1, 1),
    pvalue = runif(5)) %>%
    mutate(color_Z = ifelse(robustZ < 0, "Inhibited", "Activated"))
GS <- Scores$gs_name
load(system.file("extdata", "entrez2name.rda", package = "sSNAPPY"))

test_that("make_gs2gene_network returns warning message when expected",{
    expect_warning( make_gs2gene_network(Scores, gsTopology[names(gsTopology) %in% Scores$gs_name], colorGS_By = "robustZ", mapEntrezID =entrez2name ,
                                                      geneFC = c(1:5) ), "Genes's logFCs not provided as a named vector. All genes will be colored identically.")
})

test_that("make_gs2gene_network produces the expected outcome",{
    g_Zscore <- make_gs2gene_network(Scores, gsTopology, colorGS_By = "robustZ", mapEntrezID = NULL,
                                     geneFC = NULL, downGene_col = "blue", upGene_col = "red" )
    g_pvalue <- make_gs2gene_network(Scores, gsTopology, colorGS_By = "pvalue", mapEntrezID = NULL,
                                     geneFC = NULL, downGene_col = "blue", upGene_col = "red" )
    expect_s3_class(g_Zscore, "igraph")
    expect_true(setequal(colnames(igraph::as_data_frame(g_Zscore, "vertices")), c("name", "color", "type")))
    expect_equal(unique(V(g_Zscore)$color[-c(1:5)]), "NoFC")
    expect_true(setequal(colnames(igraph::as_data_frame(g_pvalue, "vertices")), c("name", "color", "type")))
    expect_equal(unique(V(g_pvalue)$color[-c(1:5)]), "NoFC")
    expect_true(is.character(igraph::V(g_Zscore)$color))
    expect_true(is.character(igraph::V(g_pvalue)$color))

})

test_that("plot_gs2gene returns error when expected", {

    expect_error(plot_gs2gene(Scores, colorGS_By = "random"))
    expect_error(plot_gs2gene(Scores, colorGS_By = c("community", "robustZ")))
    gsTopology_noName <- gsTopology
    names(gsTopology_noName) <- NULL
    expect_error(plot_gs2gene(Scores, gsTopology_noName, colorGS_By = "robustZ"))
    expect_error(plot_gs2gene(Scores[, -2], gsTopology, colorGS_By = "robustZ", ), "Normalised Scores must include gs_name and column for coloring")
})

test_that("plot_gs2gene produces the expected outcome", {
    expect_s3_class(plot_gs2gene(Scores, gsTopology, colorGS_By = "pvalue"), "ggraph")
    expect_s3_class(plot_gs2gene(Scores, gsTopology, colorGS_By = "robustZ"), "ggraph")
    p_1gene <- plot_gs2gene(Scores, gsTopology, colorGS_By = "pvalue", geneFC = c("ENTREZID:4128" = 0.1))
    expect_true(nrow(dplyr::filter(p_1gene$data, type == "GENE")) == 1)
    p_geneName <- plot_gs2gene(Scores, gsTopology, colorGS_By = "pvalue", mapEntrezID =entrez2name, geneFC = c("ENTREZID:4128" = 0.1))
    expect_true(dplyr::filter(p_geneName$data, type == "GENE")$name == dplyr::filter(entrez2name, entrezid == "ENTREZID:4128")$mapTo)
})
