# Stimulate logCPM matrix for 5 genes and 6 sample
# 6 samples are from 2 patients and 3 treatment levels: Control, Treat1, Treat2

y <- matrix(c(c(1:5, 2:6, 3:7), c(1:5, 2:6, 3:7)+ 0.2), 5, 6)
rownames(y) <- paste("Gene",1:5)
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


logFC <- matrix(rep(c(1, 2, 1,2), each = 5), 5, 4)
rownames(logFC) <- paste("Gene",1:5)
colnames(logFC) <- c("patient1_treat1", "patient1_treat2",  "patient2_treat1", "patient2_treat2")

weight <- rep(0.2, 5)
weighted_FC <- matrix(rep(c(0.2, 0.4 ,0.2, 0.4), each = 5), 5, 4)
rownames(weighted_FC) <- paste("ENTREZID:Gene",1:5)
colnames(weighted_FC) <- c("patient1_treat1", "patient1_treat2",  "patient2_treat1", "patient2_treat2")

# Generate logCPM matrix with random NA values
y_NA <- apply(y, 2, function(x){
    x[sample(c(1:nrow(y)), 1)] <-NA; x
})

test_that(".compute_ssFC returns erros when the treatment column is not a factor",{
    expect_error(
        .compute_ssFC(y, sample, groupBy  = "patient",
                      sampleColumn = "sample", treatColumn = "treatment"),
        "The specified treatment column must be a factor with at least 2 levels")
})

# Change the treatment column to a factor
sample <- sample %>%
    mutate(treatment = factor(treatment, levels = c("control", "treat1", "treat2")))

test_that(".compute_ssFC produces expected output", {
    ssFC <- .compute_ssFC(y, sample, groupBy  = "patient", sampleColumn = "sample", treatColumn = "treatment")
    expect_equal(ssFC, logFC)
})

# # Generate sample metadata df that will produce error
# sample_nofactor <- colnames(y) %>%
#     as.data.frame() %>%
#     magrittr::set_colnames("sample") %>%
#     dplyr::mutate(
#         treatment = vapply(.$sample, function(x){
#             stringr::str_split(x, "_")[[1]][2]
#         }, character(1)))
# sample_notreat <- colnames(y) %>%
#     as.data.frame() %>%
#     magrittr::set_colnames("sample") %>%
#     dplyr::mutate(
#         patient = vapply(.$sample, function(x){
#             stringr::str_split(x, "_")[[1]][1]
#         }, character(1)))
# sample_noCont <- sample %>%
#     dplyr::mutate(treatment = ifelse(treatment == "control", "treatment1", treatment))
# sample_wrongDim <- sample[1:5,]
# sample_onlyContr <- sample %>%
#     dplyr::mutate(treatment = "control")


test_that(".compute_ssFC returns erros when expected",{
    expect_error(
        .compute_ssFC(y, sample, groupBy  = "patient", sampleColumn = "samle", treatColumn = "treatment"),
        "Sample metadata does not contain the sample name column specified")
    expect_error(
        .compute_ssFC(y, sample, groupBy  = "pat", sampleColumn = "sample", treatColumn = "treatment"),
        "Sample metadata must include the columns matching the treatColumn and groupBy parameter")
    expect_error(
        .compute_ssFC(y, sample, groupBy  = "patient", sampleColumn = "patient", treatColumn = "treatment"),
        "Sample names in the metadaata does not match with logCPM's column names")
    expect_error(
        .compute_ssFC(y_NA, sample, groupBy  = "patient", sampleColumn = "sample", treatColumn = "treatment"),
        "NA values not allowed")
})

test_that("weight_ss_fc produces expected output", {
    output <- weight_ss_fc(y, sample, groupBy  = "patient", sampleColumn = "sample", treatColumn = "treatment")
    expect_equal(output$weighted_logFC, weighted_FC)
    expect_equal(output$weight, weight)
})

test_that("weight_ss_fc returns erros when expected", {
    expect_error(weight_ss_fc(y,  groupBy  = "patient", sampleColumn = "sample", treatColumn = "treatment"), "sample metadata must be provided as a data frame")
})

test_that("Test DGEList input for weight_ss_fc", {
    dge <- edgeR::DGEList(counts = y, samples = sample)
    output <- weight_ss_fc(dge,  groupBy  = "patient", sampleColumn = "sample", treatColumn = "treatment")
    expect_equal(dim(output$weighted_logFC), c(5, 4))
})

test_that("Test SummarizedExperiment input for weight_ss_fc", {
    dge <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=y),
                                                      colData = sample)
    output <- weight_ss_fc(dge, groupBy  = "patient", sampleColumn = "sample", treatColumn = "treatment")
    expect_equal(dim(output$weighted_logFC), c(5, 4))
})
