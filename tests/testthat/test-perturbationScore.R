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
pathwayDir <- "data/BminsI.rds"

.ssPertScore( BminsI, ssFC$logFC)
# create logCPM matrix with gene_id as rownames (instead of entrezID required)
y_wrongIdentifier <- y
rownames(y_wrongIdentifier) <- C("ENSG00000000003","ENSG00000000419","ENSG00000000457","ENSG00000000460","ENSG00000000938")
ssFC_wrongIdentifier <- weight_ssFC(y_wrongIdentifier, sample, "patient", "control")
test_that("perturbationScore returns error when expected", {
    expect_error(perturbationScore(ssFC$logFC, "data/random.rds"), "Pathway topology matrices not detected in the specified file path. Check the file path provided.")
    expect_error(perturbationScore(ssFC$logFC_wrongIdentifier, pathwayDir), "Weighted ssFCs and pathwy topologies must use the same gene identifiers")
})
