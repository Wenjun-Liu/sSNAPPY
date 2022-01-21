# Stimulate logCPM matrix for 10 genes and 6 sample
# 6 samples are from 2 patients and 3 treatment levels: Control, Treat1, Treat2

# y <- matrix(rnbinom(60,size=1,mu=10),10,6)
# y <- cpm(y)
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

# Generate sample metadata df that will produce error
sample_nofactor <- colnames(y) %>%
    as.data.frame() %>%
    magrittr::set_colnames("sample") %>%
    dplyr::mutate(
        treatment = vapply(.$sample, function(x){
            stringr::str_split(x, "_")[[1]][2]
        }, character(1)))
sample_notreat <- colnames(y) %>%
    as.data.frame() %>%
    magrittr::set_colnames("sample") %>%
    dplyr::mutate(
        patient = vapply(.$sample, function(x){
            stringr::str_split(x, "_")[[1]][1]
        }, character(1)))
sample_noCont <- sample %>%
    dplyr::mutate(treatment = ifelse(treatment == "control", "treatment1", treatment))
sample_wrongDim <- sample[1:5,]
sample_onlyContr <- sample %>%
    dplyr::mutate(treatment = "control")

test_that("compute_ssFC produces expected output", {
    ssFC <- compute_ssFC(y, sample, factor = "patient", control = "control")
    expect_equal(ssFC, logFC)
})

test_that("compute_ssFC returns erros when expected",{
    expect_error(compute_ssFC(y, sample_wrongDim, factor = "patient", control = "control"), "Sample metadata dimesion does not match with logCPM")
    expect_error(compute_ssFC(y, sample_nofactor, factor = "patient", control = "control"), "factor %in% colnames(metadata) is not TRUE ")
    expect_error(compute_ssFC(y, sample_notreat, factor = "patient", control = "control"), "Sample metadata must contain a column named treatment")
    expect_error(compute_ssFC(y, sample_noCont, factor = "patient", control = "control"), "Control level not detected in sample metadata")
    expect_error(compute_ssFC(y, sample_onlyContr, factor = "patient", control = "control"), "At least 2 levels are required treatment")
})
