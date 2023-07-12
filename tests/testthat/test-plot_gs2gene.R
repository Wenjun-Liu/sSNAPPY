pathwayDir <- system.file("extdata", "gsTopology.rda", package = "sSNAPPY")
load(pathwayDir)
set.seed(123)
Scores <- data.frame(
    gs_name = c("RNA degradation","Wnt signaling pathway","Histidine metabolism","Ascorbate and aldarate metabolism","Lipid and atherosclerosis"),
    robustZ = runif(5, -1, 1),
    pvalue = runif(5)) %>%
    mutate(color_Z = ifelse(robustZ < 0, "Inhibited", "Activated"))
Scores <- dplyr::mutate(Scores, gs_name = paste("kegg.", gs_name, sep = ""))
GS <- Scores$gs_name
load(system.file("extdata", "entrez2name.rda", package = "sSNAPPY"))

test_that("make_gs2gene_network returns warning message when expected",{
    expect_message(
        .make_gs2gene_network(
            Scores, gsTopology[names(gsTopology) %in% Scores$gs_name],
            colorGsBy = "robustZ", mapEntrezID = entrez2name, geneFC = c(1:5)
        ),
        "Gene fold-changes were not provided as a named vector. All genes will be colored identically."
    )
})

test_that("make_gs2gene_network produces the expected outcome",{
    g_Zscore <- .make_gs2gene_network(
        Scores, gsTopology, colorGsBy = "color_Z", mapEntrezID = NULL, geneFC = NULL
    )
    g_pvalue <- .make_gs2gene_network(
        Scores, gsTopology, colorGsBy = "pvalue", mapEntrezID = NULL, geneFC = NULL
    )
    expect_s3_class(g_Zscore, "igraph")
    correct_cols <- c("name", "color", "fill", "type")
    expect_true(setequal(colnames(igraph::as_data_frame(g_Zscore, "vertices")), correct_cols))
    expect_equal(unique(V(g_Zscore)$color[-c(1:5)]), NA)
    expect_true(setequal(colnames(igraph::as_data_frame(g_pvalue, "vertices")), correct_cols))
    expect_equal(unique(V(g_pvalue)$color[-c(1:5)]), NA)
    expect_true(is.character(igraph::V(g_Zscore)$fill))
    expect_true(is.numeric(igraph::V(g_pvalue)$fill))

})

test_that("plot_gs2gene returns error when expected", {

    expect_error(plot_gs2gene(Scores, colorGsBy = "random"))
    expect_error(plot_gs2gene(Scores, colorGsBy = c("community", "robustZ")))
    expect_error(plot_gs2gene(Scores[, -c("gs_name")], gsTopology, colorGsBy = "pvalue"))
    gsTopology_noName <- gsTopology
    names(gsTopology_noName) <- NULL
    expect_error(plot_gs2gene(Scores, gsTopology_noName, colorGsBy = "robustZ"))
    expect_error(plot_gs2gene(Scores[, -2], gsTopology, colorGsBy = "robustZ"), "'arg' should be.+")
    expect_error(plot_gs2gene(Scores, gsTopology, colorGsBy = "pvalue",
                              gsNodeShape = 1))
    expect_error(plot_gs2gene(Scores[1,], gsTopology, colorGsBy = "pvalue"))
})

test_that("plot_gs2gene produces the expected outcome", {
    expect_s3_class(plot_gs2gene(Scores, gsTopology, colorGsBy = "pvalue"), "ggraph")
    expect_s3_class(plot_gs2gene(Scores, gsTopology, colorGsBy = NULL), "ggraph")
    expect_s3_class(plot_gs2gene(Scores, gsTopology, colorGsBy = "robustZ"), "ggraph")
    p_1gene <- plot_gs2gene(Scores, gsTopology, colorGsBy = "pvalue",
                            geneFC = c("ENTREZID:4128" = 0.1), filterGeneBy = 0)
    expect_true(nrow(dplyr::filter(p_1gene$data, type == "GENE")) == 1)
    p_NOgene <- plot_gs2gene(Scores, gsTopology, colorGsBy = "pvalue",
                            geneFC = c("ENTREZID:4128" = 0.1), filterGeneBy = 2)
    expect_true(nrow(dplyr::filter(p_NOgene$data, type == "GENE")) == 0)
    p_geneName <- plot_gs2gene(Scores, gsTopology, colorGsBy = "pvalue",
                               mapEntrezID =entrez2name,
                               geneFC = c("ENTREZID:4128" = 0.1),
                               filterGeneBy = 0)
    expect_true(dplyr::filter(p_geneName$data, type == "GENE")$name == dplyr::filter(entrez2name, entrezid == "ENTREZID:4128")$mapTo)
})

