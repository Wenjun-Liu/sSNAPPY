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
ssFC <- weight_ssFC(y, sample, "patient", "control")
pathwayDir <- system.file("testdata", "test_BminsI.rda", package = "SSPT")
load(pathwayDir)
# the number of pathways with at least one of those five genes in it
interesectName <- names(BminsI[lapply(BminsI, function(x){length(intersect(rownames(ssFC$logFC),rownames(x)))}) != 0])

y_withNA <- y
y_withNA[2,2] <- NA
# create logCPM matrix with gene_id as rownames (instead of entrezID required)
y_wrongIdentifier <- y
rownames(y_wrongIdentifier) <- c("ENSG00000000003","ENSG00000000419","ENSG00000000457","ENSG00000000460","ENSG00000000938")
ssFC_wrongIdentifier <- weight_ssFC(y_wrongIdentifier, sample, "patient", "control")

test_that("generate_PermutedScore returns error when expected", {
    expect_error(generate_PermutedScore(y_wrongIdentifier, numOfTreat =  3,
                                        NB = 1000,
                                        gsTopology = BminsI, weight = ssFC_wrongIdentifier$weight), "None of the expressed gene was matched to pathways. Check if gene identifiers match")
    expect_error(generate_PermutedScore(y, numOfTreat = 3,
                                        NB = 1000,
                                        gsTopology = BminsI, weight =ssFC$weight[1:10]), "Gene-wise weights do not match with the dimension of logCPM")
    expect_error(generate_PermutedScore(y, numOfTreat = 4,
                                        NB = 1000,
                                        gsTopology = BminsI, weight =ssFC$weight), "Number of samples must be divisible by the number of treatments")
    expect_error(generate_PermutedScore(y_withNA, numOfTreat = 3,
                                        NB = 1000,
                                        gsTopology = BminsI, weight = ssFC$weight), "NA values not allowed")
})

test_that("generate_PermutedScore produces the expected outcome", {
    results_sub <- generate_PermutedScore(y[, 1:4], numOfTreat =2, NB = 100, gsTopology = BminsI, weight = ssFC$weight)
    expect_equal(length(results_sub), length(BminsI))
    expect_equal(length(results_sub[[1]]), factorial(4)*2)
})

test_that("generate_PermutedScore produces the expected outcome", {
    results <- generate_PermutedScore(y, numOfTreat = 3, NB = 10, gsTopology = BminsI, weight = ssFC$weight)
    expect_equal(length(results[[1]]), 10*(6-3+1))
})

test_that("normaliseByPermutation produces the expected outcome",{
    perS <- list(
        "Chemokine signaling pathway"= rnorm(40, mean = 1, sd = 0.3)
    )
    ssPertScore <- perturbationScore(ssFC$logFC, BminsI)
    output <- normaliseByPermutation(perS, ssPertScore)
    expect_equal(unique(output$sample), c("patient1_treat1", "patient1_treat2", "patient2_treat1", "patient2_treat2"))
    expect_false(anyNA(output$robustZ))
    expect_true(length(intersect(output$MAD, 0)) == 0)
})
