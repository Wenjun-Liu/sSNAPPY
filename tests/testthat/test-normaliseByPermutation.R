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
intersectName <- names(gsTopology[lapply(gsTopology, function(x){length(intersect(rownames(ssFC$logFC),rownames(x)))}) != 0])

y_withNA <- y
y_withNA[2,2] <- NA
# create logCPM matrix with gene_id as rownames (instead of entrezID required)
y_wrongIdentifier <- y
rownames(y_wrongIdentifier) <- c("ENSG00000000003","ENSG00000000419","ENSG00000000457","ENSG00000000460","ENSG00000000938")
ssFC_wrongIdentifier <- weight_ss_fc(y_wrongIdentifier, sample, "patient", "control")

test_that("generate_permuted_scores returns error when expected", {
    expect_error(generate_permuted_scores(y_wrongIdentifier, numOfTreat =  3,
                                        NB = 1000,
                                        gsTopology = gsTopology, weight = ssFC_wrongIdentifier$weight), "None of the expressed gene was matched to pathways. Check if gene identifiers match")
    expect_error(generate_permuted_scores(y, numOfTreat = 3,
                                        NB = 1000,
                                        gsTopology = gsTopology, weight =ssFC$weight[1:10]), "Gene-wise weights do not match with the dimension of expreMatrix")
    expect_error(generate_permuted_scores(y, numOfTreat = 4,
                                        NB = 1000,
                                        gsTopology = gsTopology, weight =ssFC$weight), "Number of samples must be divisible by the number of treatments")
    expect_error(generate_permuted_scores(y_withNA, numOfTreat = 3,
                                        NB = 1000,
                                        gsTopology = gsTopology, weight = ssFC$weight), "NA values not allowed")
})

test_that(".generate_permutedFC produces the expected outcome", {
    temp <- .generate_permutedFC(y, numOfTreat = 3,
                                 NB = 2, weight = ssFC$weight)
    expect_equal(length(temp), 2)
    expect_equal(ncol(temp[[1]]), ncol(y) - (ncol(y)/3))
})


test_that(".generate_permutedFC produces the expected outcome", {
    temp <- .generate_permutedFC(y, numOfTreat = 3,
                                 NB = 2, weight = ssFC$weight)
    expect_equal(length(temp), 2)
    expect_equal(ncol(temp[[1]]), ncol(y) - (ncol(y)/3))
})


test_that("Test data.frame input for generate_permuted_scores", {
    temp <- generate_permuted_scores(as.data.frame(y), numOfTreat = 3,
                                 NB = 2, weight = ssFC$weight, gsTopology = gsTopology[intersectName])
    expect_equal(length(temp), length(intersectName))
    expect_equal(length(temp[[1]]), (ncol(y) - (ncol(y)/3))*2)
})

test_that("Test DGEList input for generate_permuted_scores", {
    dge <- edgeR::DGEList(counts = y)
    temp <- generate_permuted_scores(dge, numOfTreat = 3,
                                   NB = 2, weight = ssFC$weight, gsTopology = gsTopology[intersectName])
    expect_equal(length(temp), length(intersectName))
})


test_that("Test SummarizedExperiment input for generate_permuted_scores", {
    dge <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=y))
    temp <- generate_permuted_scores(dge, numOfTreat = 3,
                                   NB = 2, weight = ssFC$weight, gsTopology = gsTopology[intersectName])
    expect_equal(length(temp), length(intersectName))
})

notExpressed <- setdiff(unique(unlist(unname(lapply(gsTopology, rownames)))), rownames(y))
if (length(notExpressed) != 0){
    temp <- matrix(0, nrow = length(notExpressed), ncol = ncol(y))
    rownames(temp) <- notExpressed
    colnames(temp) <- colnames(y)
    y <- rbind(y, temp)
    ssFC$weight <- c(ssFC$weight, rep(0, length(notExpressed)))

}

permutedFC <- .generate_permutedFC(y, numOfTreat = 3,
                                   NB = 2, weight = ssFC$weight)
test_that("permutedPertScore_RCPP produces the expected outcome", {
    ls <- permutedPertScore_RCPP(X = gsTopology[[1]],
                                 pathwayG = rownames(gsTopology[[1]]),
                                 rownames(y),
                                 permutedFC = permutedFC, ncol(permutedFC[[1]]))
    expect_equal(length(ls), 2)
    expect_equal(length(ls[[1]]), ncol(y) - (ncol(y)/3))
})


test_that("generate_permuted_scores produces the expected outcome", {
    results <- generate_permuted_scores(y, numOfTreat = 3, NB = 3, gsTopology = gsTopology, weight = ssFC$weight)
    expect_equal(length(results[[1]]), 3*(6-3+1))
})

test_that("normalise_by_permu produces the expected outcome",{
    perS <- list(
        "Chemokine signaling pathway"= rnorm(40, mean = 1, sd = 0.3)
    )
    genePertScore <- raw_gene_pert(ssFC$logFC, gsTopology)
    ssPertScore <- pathway_pert( genePertScore)
    output <- normalise_by_permu(perS, ssPertScore)
    expect_equal(unique(output$sample), c("patient1_treat1", "patient1_treat2", "patient2_treat1", "patient2_treat2"))
    expect_false(anyNA(output$robustZ))
    expect_true(length(intersect(output$MAD, 0)) == 0)
})
