#' @title Compute weighted single sample LogFCs from normalised logCPM
#'
#' @description Compute weighted single sample logFCs for each treated samples using normalised logCPM values. Fit a lowess
#' curve on variance of single sample logFCs ~ mean of logCPM, and use this to predict a gene-wise weight. The weighted
#' single sample logFCs are ready for computing perturbation scores.
#'
#' @details
#'
#' This function computes weighted single sample logFCs from normalised logCPM values, used for computing single sample perturbation scores.
#' Since genes with smaller logCPM turn to have a larger variance among single sample logFCs. A lowess curve is fitted to estimate the
#' relationship between variance of single sample logFCs and mean of logCPM, and use this relationship to estimate the variance of each
#' mean logCPM value. Weights, which are inverse variance, are then multiplied to single sample logFCs to downweight genes with low counts.
#'
#' It is assumed that the genes with extremely low counts have been removed and the count matrix has been normalised prior to logCPM
#' matrix was derived. Rownames of the matrix must be genes' entrez ID. To convert ensemble ID to entrz ID, see example.
#'
#' The number of rows in the sample metadata must equal to the number of columns in the logCPM matrix.
#' Metadata also must have a column called "sample" storing sample names, and a column called "treatment" storing
#' treatment of each sample.The control treatment level specified by `control` parameter must exist in the treatment column.
#'
#' This analysis was designed for experimental designs that include matched pairs of samples, such as when the tissues collected from the
#' same patient were treated with different treatments. Parameter `factor` tells the function how samples can be put into pairs. It must also
#' be a column of the metadata.
#'
#' @param logCPM Matrix of normaslised logCPM where rows are genes and columns are samples. Rownames need to be gene entrez ID.
#' @param metadata Sample metadata data frame as described in the details section.
#' @param factor The factor defines how samples can be put into matching pairs.
#' @param control The treatment level that is the control.
#'
#' @importFrom stats approxfun lowess
#' @return A list
#' @export
weight_ssFC <- function(logCPM, metadata, factor, control){

    ssFC <- .compute_ssFC(logCPM, metadata, factor, control)
    varFC <- apply(ssFC, 1, var)
    meanCPM <- apply(logCPM, 1, mean)

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
}


#' @title Compute single sample logFCs
#' @param logCPM Matrix of normaslised logCPM
#' @param metadata Sample metadata data frame as described in the details section.
#' @param factor The factor defines how samples can be put into matching pairs.
#' @param control The treatment level that is the control.
#' @importFrom dplyr pull filter
#' @importFrom rlang sym
#' @return A matrix of single sample logFC
#' @keywords internal
.compute_ssFC <- function(logCPM, metadata, factor, control){

    ## checks
    # if (is.null(logCPM)) stop("LogCPM has to be provided")
    # if (is.null(metadata)) stop("Sample metadata has to be provided")
    if (missing(factor)) stop("Factor defining matching samples must be provided")
    if (missing(control)) stop("Control treatment must be specified")
    if (!"treatment" %in% colnames(metadata)) stop("Sample metadata must contain a column named treatment")
    if (!control %in% unique(metadata$treatment)) stop("Control level not detected in sample metadata")
    if (!"sample" %in% colnames(metadata)) stop ("Sample name must be specific in a column named sample")
    stopifnot(factor %in% colnames(metadata))
    stopifnot(ncol(logCPM) == nrow(metadata))
    m <- min(logCPM)
    if (is.na(m)) stop("NA values not allowed")

    logCPM <- as.matrix(logCPM)
    metadata <- as.data.frame(metadata)
    if(length(unique(metadata[,"treatment"])) <2) stop("At least 2 levels are required treatment")
    pairs <- unique(metadata[,factor])
    sapply(pairs, function(x){
       contrSample <- dplyr::filter(metadata, treatment == control, !!sym(factor) == x)
       contrSample <- pull(contrSample, sample)

       treatedSample <- dplyr::filter(metadata, treatment != control, !!sym(factor) == x)
       treatedSample <- pull(treatedSample, sample)

       logCPM[, treatedSample] - logCPM[, contrSample]
    }, simplify = FALSE) %>%
        do.call(cbind,.)

}


