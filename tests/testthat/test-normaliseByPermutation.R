# Stimulate logCPM matrix for 5 genes and 6 sample
# 6 samples are from 2 patients and 3 treatment levels: Control, Treat1, Treat2

y <- matrix(c(c(1:5, 2:6, 3:7), c(1:5, 2:6, 3:7)+ 0.2), 5, 6)
rownames(y) <- c("7105","8813","57147","55732","2268" )
colnames(y) <- c("patient1_control", "patient1_treat1", "patient1_treat2",
                 "patient2_control", "patient2_treat1", "patient2_treat2")
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
# the number of pathways with at least one of those five genes in it
intersectName <- names(gsTopology[lapply(gsTopology, function(x){length(intersect(rownames(ssFC$weighted_logFC),rownames(x)))}) != 0])

y_withNA <- y
y_withNA[2,2] <- NA
# create logCPM matrix with gene_id as rownames (instead of entrezID required)
y_wrongIdentifier <- y
rownames(y_wrongIdentifier) <- c("ENSG00000000003","ENSG00000000419","ENSG00000000457","ENSG00000000460","ENSG00000000938")
ssFC_wrongIdentifier <- weight_ss_fc(y_wrongIdentifier, sample, groupBy  = "patient", sampleColumn = "sample", treatColumn = "treatment")

test_that("generate_permuted_scores returns error when expected", {
    expect_error(generate_permuted_scores(y_wrongIdentifier, gsTopology = gsTopology,
                                          weight = ssFC_wrongIdentifier$weight), "None of the expressed gene was matched to pathways. Check if gene identifiers match")
    expect_error(generate_permuted_scores(y, gsTopology = gsTopology,
                                          weight =ssFC$weight[1:10]), "Gene-wise weights do not match with the dimension of expreMatrix")
    expect_error(generate_permuted_scores(y_withNA, gsTopology = gsTopology,
                                          weight = ssFC$weight), "NA values not allowed")
    # create a testScore data.frame containing wrong gs_name
    testScore <- data.frame(
        gs_name = 1:5
    )
    expect_error(generate_permuted_scores(y, gsTopology = gsTopology,
                                          weight = ssFC$weight,
                                          testScore = testScore),
                 "None of the pathways had non-zero test perturbation scores")
})

test_that(".generate_permutedFC produces the expected outcome", {
    # if NB is provided, the number of column in the returned output should
    # equal to the specified NB
    temp <- .generate_permutedFC(y, NB = 2, weight = ssFC$weight)
    expect_equal(dim(temp), c(nrow(y), 2))
})

test_that("generate_permuted_scores produces the expected outcome", {
    genePertScore <- raw_gene_pert(ssFC$weighted_logFC, gsTopology)
    ssPertScore <- pathway_pert( genePertScore)
    temp <- generate_permuted_scores(y, NB = 2, weight = ssFC$weight,
                                     gsTopology = gsTopology,
                                     testScore = ssPertScore)
    expect_equal(length(temp), length(unique(ssPertScore$gs_name)))
    expect_equal(length(temp[[1]]), 2)

    # when NB isn't specified, the number of permuted scores generated for
    # each pathway should equal sample size x (sample size -1)
    temp <- generate_permuted_scores(y, weight = ssFC$weight,
                                     gsTopology = gsTopology,
                                     testScore = ssPertScore)
    expect_equal(length(temp[[1]]), ncol(y)*(ncol(y) -1))

    # When the column number isn't an even number, the maximum permutation pairs
    # should still be sample size x (sample size -1)
    temp <- generate_permuted_scores(y[,-1], weight = ssFC$weight,
                                     gsTopology = gsTopology,
                                     testScore = ssPertScore)
    expect_equal(length(temp[[1]]), (ncol(y) -1)*(ncol(y) -2))

    # If the required NB is bigger than the maxmimum number of permutations
    # possible, exhuast permutation choice
    temp <- generate_permuted_scores(y,  NB = 200,
                                     weight = ssFC$weight,
                                     gsTopology = gsTopology,
                                     testScore = ssPertScore)
    expect_equal(length(temp[[1]]), ncol(y)*(ncol(y) -1))


})

test_that("Test data.frame input for generate_permuted_scores", {
    temp <- generate_permuted_scores(as.data.frame(y),
                                     weight = ssFC$weight,
                                 gsTopology = gsTopology[intersectName])
    expect_equal(length(temp), length(intersectName))
    expect_equal(length(temp[[1]]), ncol(y)*(ncol(y) -1))
})

test_that("Test DGEList input for generate_permuted_scores", {
    dge <- edgeR::DGEList(counts = y)
    temp <- generate_permuted_scores(dge,  NB = 2, weight = ssFC$weight,
                                   gsTopology = gsTopology[intersectName])
    expect_equal(length(temp), length(intersectName))
})


test_that("Test SummarizedExperiment input for generate_permuted_scores", {
    dge <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=y))
    temp <- generate_permuted_scores(dge, NB = 2, weight = ssFC$weight,
                                   gsTopology = gsTopology[intersectName])
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

permutedFC <- .generate_permutedFC(y, NB = 2, weight = ssFC$weight)


test_that("normalise_by_permu produces the expected outcome",{
    perS <- list(
        "kegg.Chemokine signaling pathway"= rnorm(40, mean = 1, sd = 0.3)
    )
    genePertScore <- raw_gene_pert(ssFC$weighted_logFC, gsTopology)
    ssPertScore <- pathway_pert( genePertScore)
    output <- normalise_by_permu(perS, ssPertScore)
    expect_equal(levels(output$sample),
                 c("patient1_treat1", "patient1_treat2", "patient2_treat1", "patient2_treat2"))
    expect_false(anyNA(output$robustZ))
    expect_false(any(output$MAD == 0))
    expect_true(is.factor(output$gs_name))

    output_sorted <- normalise_by_permu(perS, ssPertScore, sortBy = "pvalue")
    expect_equal(which.min(pull(output_sorted, pvalue)), 1)
})

