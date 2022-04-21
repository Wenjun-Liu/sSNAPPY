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
# the number of pathways with at least one of those five genes in it
interesectName <- names(gsTopology[lapply(gsTopology, function(x){length(intersect(rownames(ssFC$logFC),rownames(x)))}) != 0])

# create logCPM matrix with gene_id as rownames (instead of entrezID required)
y_wrongIdentifier <- y
rownames(y_wrongIdentifier) <- c("ENSG00000000003","ENSG00000000419","ENSG00000000457","ENSG00000000460","ENSG00000000938")
ssFC_wrongIdentifier <- weight_ss_fc(y_wrongIdentifier, sample, "patient", "control")

test_that("computePerturbationScore returns error when expected", {
    expect_error(computePerturbationScore(ssFC_wrongIdentifier$logFC, gsTopology), "None of the expressed gene was matched to pathways. Check if gene identifiers match")
})

notExpressed <- setdiff(unique(unlist(unname(lapply(gsTopology, rownames)))), rownames(ssFC$logFC))
if (length(notExpressed) != 0){
    # set the FCs of unexpressed pathway genes to 0
    temp <- matrix(0, nrow = length(notExpressed), ncol = ncol(ssFC$logFC))
    rownames(temp) <- notExpressed
    colnames(temp) <- colnames(ssFC$logFC)
    # set the weights of unexpressed pathway genes to 0
    ssFC$logFC <- rbind(ssFC$logFC, temp)}

test_that("ssPertScore_RCPP produces the expected outcome",{
    ls <- ssPertScore_RCPP(gsTopology, ssFC$logFC, rownames(ssFC$logFC), colnames(ssFC$logFC))
    expect_equal(names(ls), names(gsTopology))
    expect_true(is.vector(ls[[1]]))
    expect_equal(names(ls[[1]]), stringr::str_subset(sample$sample, "control", negate = TRUE))
})

test_that("computePerturbationScore produces the expected outcome", {
    output <- computePerturbationScore(ssFC$logFC, gsTopology)
    expect_equal(colnames(output), c("sample", "tA", "gs_name"))
    expect_false(anyNA(output$tA))
    expect_equal(unique(output$sample), stringr::str_subset(sample$sample, "control", negate = TRUE))
    expect_true(length(setdiff(output$gs_name, interesectName)) == 0)
})
