# Stimulate logCPM matrix for 5 genes and 6 sample
# 6 samples are from 2 patients and 3 treatment levels: Control, Treat1, Treat2

y <- matrix(c(c(1:5, 2:6, 3:7), c(1:5, 2:6, 3:7)+ 0.2), 5, 6)
rownames(y) <- c("7105","8813","57147","55732","2268" )
colnames(y) <- c("patient1_control", "patient1_treat1", "patient1_treat2", "patient2_control", "patient2_treat1", "patient2_treat2")
sample <- colnames(y) %>%
    as.data.frame()
colnames(sample) <- c("sample")
sample <- sample %>%
    dplyr::mutate(
        treatment = vapply(.$sample, function(x){
            stringr::str_split(x, "_")[[1]][2]
        }, character(1)),
        patient = vapply(.$sample, function(x){
            stringr::str_split(x, "_")[[1]][1]
        }, character(1)))
ssFC <- weight_ss_fc(y, sample, "patient", "control")
pathwayDir <- system.file("extdata", "gsTopology.rda", package = "sSNAPPY")
load(pathwayDir)
load(system.file("extdata", "entrez2name.rda", package = "sSNAPPY"))
# the number of pathways with at least one of those five genes in it
interesectName <- names(gsTopology[lapply(gsTopology, function(x){length(intersect(rownames(ssFC$logFC),rownames(x)))}) != 0])
# compute raw gene-wise perturbation scores
genePertScore <- raw_gene_pert(ssFC$logFC, gsTopology)
# sum gene-wise perturbation scores to derive the pathway-level single-sample perturbation scores
pathwayPertScore <- pathway_pert(genePertScore)

test_that("plot_gene_contribution returns error when expected", {
    expect_error(plot_gene_contribution( genePertScore, "random name"), "Gene-wise perturbation scores not provided for the chosen pathway.")
})

test_that("plot_gene_contribution returns a pheatmap object as expected", {
    hp <- plot_gene_contribution(genePertScore, gsToPlot = "Chemokine signaling pathway", metadata = sample,
                                 annotation_attribute = c("pathwayPertScore", "treatment"))
    expect_equal(class(hp), "pheatmap")
    hp2 <- plot_gene_contribution(genePertScore, gsToPlot = "Chemokine signaling pathway",
                                 annotation_attribute = NULL)
    expect_equal(class(hp2), "pheatmap")
    hp3 <- plot_gene_contribution(genePertScore, gsToPlot = "Chemokine signaling pathway")
    expect_equal(class(hp3), "pheatmap")

    hp4 <- plot_gene_contribution(genePertScore, gsToPlot = "Chemokine signaling pathway", annotation_attribute = "treatment", metadata = sample)
    expect_equal(class(hp4), "pheatmap")

    hp5 <- plot_gene_contribution(genePertScore, gsToPlot = "Chemokine signaling pathway", pathwayPertScore = pathwayPertScore)
    expect_equal(class(hp5), "pheatmap")

    randomRownames <- sample(1:1000, nrow(genePertScore$`Chemokine signaling pathway`))
    hp6 <- plot_gene_contribution(genePertScore, gsToPlot = "Chemokine signaling pathway", pathwayPertScore = pathwayPertScore, mapRownameTo = randomRownames)
    expect_equal(class(hp6), "pheatmap")

    hp7 <- plot_gene_contribution(genePertScore, gsToPlot = "Chemokine signaling pathway", annotation_attribute = c("pathwayPertScore", "treatment"),
                                  pathwayPertScore = pathwayPertScore)
    expect_equal(class(hp7), "pheatmap")

    hp8 <- plot_gene_contribution(genePertScore, gsToPlot = "Chemokine signaling pathway", annotation_attribute = c("pathwayPertScore", "treatment"),
                                  pathwayPertScore = pathwayPertScore, metadata = sample)
    expect_equal(class(hp8), "pheatmap")

    hp9 <- plot_gene_contribution(genePertScore, gsToPlot = "Chemokine signaling pathway", annotation_attribute = c("pathwayPertScore", "treatment"),
                                   metadata = sample)
    expect_equal(class(hp9), "pheatmap")

    hp10 <- plot_gene_contribution(genePertScore, gsToPlot = "Chemokine signaling pathway", annotation_attribute = c("pathwayPertScore", "treatment"),
                                  mapEntrezID = entrez2name)
    expect_equal(class(hp10), "pheatmap")
})

