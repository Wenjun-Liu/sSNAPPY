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
        treatment = factor(treatment, levels = c("control", "treat1", "treat2")),
        patient = vapply(.$sample, function(x){
            stringr::str_split(x, "_")[[1]][1]
        }, character(1)))
ssFC <- weight_ss_fc(y, sample, groupBy  = "patient", sampleColumn = "sample", treatColumn = "treatment")
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
    expect_error(plot_gene_contribution( genePertScore$`kegg.Chemokine signaling pathway`, topGene = "random"))
    expect_warning(
        plot_gene_contribution( genePertScore$`kegg.Chemokine signaling pathway`, mapEntrezID = entrez2name[1,])
    )
    expect_message(
        plot_gene_contribution(genePertScore$`kegg.Chemokine signaling pathway`,
                               annotation_df = dplyr::select(sample, -sample))
    )
})

test_that("plot_gene_contribution returns a pheatmap object as expected", {
    hp <- plot_gene_contribution(genePertScore$`kegg.Chemokine signaling pathway`)
    expect_equal(class(hp), "pheatmap")
    hp2 <- plot_gene_contribution(genePertScore$`kegg.Chemokine signaling pathway`,
                                  annotation_df = sample)
    expect_equal(class(hp2), "pheatmap")
    hp3 <- plot_gene_contribution(genePertScore$`kegg.Chemokine signaling pathway`,
                                  annotation_df = sample, mapEntrezID = entrez2name)
    expect_equal(class(hp3), "pheatmap")
})

test_that("incorrect function errors",{
    expect_error(
        plot_gene_contribution(genePertScore$`kegg.Chemokine signaling pathway`, filterBy = "")
    )
})