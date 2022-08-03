#' @title Compute weighted single sample LogFCs from normalised logCPM
#'
#' @description Compute weighted single sample logFCs for each treated samples
#' using normalised logCPM values. Fit a lowess curve on variance of single
#' sample logFCs ~ mean of logCPM, and use it to predict a gene-wise weight.
#' The weighted single sample logFCs are ready for computing perturbation scores.
#'
#' @details
#'
#' This function computes weighted single sample logFCs from normalised logCPM
#' values, used for computing single sample perturbation scores. Since genes with
#' smaller logCPM turn to have a larger variance among single sample logFCs.
#' A lowess curve is fitted to estimate the relationship between variance of
#' single sample logFCs and mean of logCPM, and the relationship is used to estimate
#' the variance of each mean logCPM value. Gene-wise weights, which are inverse of
#' variances, are then multiplied to single sample logFCs to downweight genes with
#' low counts. It is assumed that the genes with extremely low counts have been
#' removed and the count matrix has been normalised prior to logCPM matrix was
#' derived. Rownames of the matrix must be genes' entrez ID. To convert other gene
#' identifiers to entrz ID, see example.
#'
#' If a S4 object of \code{DGEList or SummarizedExperiment} is provided as input
#' to `expreMatrix`, gene expression matrix will be extracted from it and converted
#' to logCPM matrix. Sample metadata will also be extracted from the same S4 object
#' unless otherwise specified.
#'
#' Provided sample metadata should have the same number of rows as the number of columns in
#' the logCPM matrix. Metadata also must have a column called "sample" storing
#' sample names (column names of logCPM matrix), and a column called "treatment"
#' storing treatment of each sample.The control treatment level specified by `control`
#' parameter must exist in the treatment column.
#'
#' This analysis was designed for experimental designs that include matched pairs of
#' samples, such as when tissues collected from the same patient were treated with
#' different treatments to study different treatment effects. Parameter `factor` tells
#' the function how samples can be put into matching pairs. It must also be included as
#' a column in the metadata.
#'
#' @param expreMatrix matrix and data.frame of logCPM, or DGEList/SummarizedExperiment
#' storing gene expression counts and sample metadata. Feature names need to be
#' gene entrez IDs, and column names need to be sample names
#' @param metadata Sample metadata data frame as described in the details section.
#' @param factor Factor defines how samples can be put into matching pairs (eg. patient).
#' @param control Treatment level that is the control.
#'
#' @importFrom stats approxfun lowess var
#' @return A list with two elements:
#' $weight  gene-wise weights;
#' $logFC weighted single sample logFC matrix
#' @examples
#' # Inspect metadata data frame to make sure it has treatment, sample and patient columns
#' data(metadata_example)
#' data(logCPM_example)
#' length(setdiff(colnames(logCPM_example), metadata_example$sample)) == 0
#' ls <- weight_ss_fc(logCPM_example, metadata = metadata_example,
#'  factor = "patient", control = "Vehicle")
#' @export
setGeneric("weight_ss_fc", function(expreMatrix, metadata = NULL, factor, control)
    standardGeneric("weight_ss_fc"))

#' @rdname weight_ss_fc
setMethod("weight_ss_fc",
          signature = signature(expreMatrix = "matrix"),
          function(expreMatrix, metadata = NULL, factor, control){
              if (is.null(metadata)) stop("sample metadata must be provided")
              ssFC <- .compute_ssFC(expreMatrix, metadata, factor, control)
              varFC <- apply(ssFC, 1, var)
              meanCPM <- apply(expreMatrix, 1, mean)

              # make sure varFC & meanCPM are in correct order
              meanCPM <- meanCPM[match(names(meanCPM),names(varFC))]
              l <- lowess(meanCPM, varFC)
              f <- approxfun(l, rule = 2, ties = list("ordered", mean))
              weight <- 1/f(meanCPM)
              scaled_w <- weight/sum(weight)
              rownames(ssFC) <- paste("ENTREZID:", rownames(ssFC), sep = "")
              list(
                  weight = scaled_w,
                  logFC = ssFC * scaled_w
              )
          })

#' @rdname weight_ss_fc
setMethod("weight_ss_fc",
          signature = signature(expreMatrix = "data.frame"),
          function(expreMatrix, metadata = NULL, factor, control){
              weight_ss_fc(as.matrix(expreMatrix), metadata, factor, control)
          })

#' @rdname weight_ss_fc
setMethod("weight_ss_fc",
          signature = signature(expreMatrix = "DGEList"),
          function(expreMatrix, metadata = NULL, factor, control){
              cpm <- cpm(expreMatrix$counts, log = TRUE)
              if(is.null(metadata)){
                  metadata <- expreMatrix$samples
              }
              weight_ss_fc(cpm, metadata, factor, control)
          })

#' @rdname weight_ss_fc
setMethod("weight_ss_fc",
          signature = signature(expreMatrix = "SummarizedExperiment"),
          function(expreMatrix, metadata = NULL, factor, control){
              cpm <- cpm(SummarizedExperiment::assay(expreMatrix), log = TRUE)
              if(is.null(metadata)){
                  metadata <- SummarizedExperiment::colData(expreMatrix)
              }
              weight_ss_fc(cpm, metadata, factor, control)
          })

#' @title Compute single sample logFCs
#' @param logCPM Matrix of normaslised logCPM
#' @param metadata Sample metadata data frame as described in the details section.
#' @param factor The factor defines how samples can be put into matching pairs.
#' @param control The treatment level that is the control.
#' @importFrom dplyr pull filter
#' @importFrom rlang sym
#' @importFrom magrittr set_colnames
#' @return A matrix of single sample logFC
#' @keywords internal
.compute_ssFC <- function(logCPM, metadata, factor, control){

    metadata <- as.data.frame(metadata)

    # checks
    if (!all(c("treatment", "sample", factor) %in% colnames(metadata))) stop("Sample metadata must include factor, treatment and sample")
    if (any(c(!control %in% unique(metadata$treatment), length(unique(metadata[,"treatment"])) <2))) stop(
        "Treatment needs at least 2 levels where one is the control specified")
    if (!setequal(colnames(logCPM), metadata[,"sample"])) stop("Sample metadaata does not match with logCPM's column names")
    m <- min(logCPM)
    if (is.na(m)) stop("NA values not allowed")

    pairs <- unique(as.character(pull(metadata, sym(factor))))
    ls <- lapply(pairs, function(x){
        contrSample <- dplyr::filter(metadata, metadata$treatment == control, !!sym(factor) == x)
        contrSample <- as.character(pull(contrSample, sample))
        treatedSample <- dplyr::filter(metadata, metadata$treatment != control, !!sym(factor) == x)
        treatedSample <- as.character(pull(treatedSample, sample))

        if (length(unique(metadata[,"treatment"])) == 2){
            set_colnames(as.matrix(logCPM[, treatedSample] - logCPM[, contrSample]), x)
        } else {
            logCPM[, treatedSample] - logCPM[, contrSample]
        }

        })
    do.call(cbind,ls)
}


