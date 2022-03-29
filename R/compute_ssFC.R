#' @title Compute weighted single sample LogFCs from normalised logCPM
#'
#' @description Compute weighted single sample logFCs for each treated samples using normalised
#' logCPM values. Fit a lowess curve on variance of single sample logFCs ~ mean of logCPM, and use
#' it to predict a gene-wise weight. The weighted single sample logFCs are ready for computing perturbation scores.
#'
#' @details
#'
#' This function computes weighted single sample logFCs from normalised logCPM values, used for
#' computing single sample perturbation scores. Since genes with smaller logCPM turn to have a larger
#' variance among single sample logFCs. A lowess curve is fitted to estimate the relationship between
#' variance of single sample logFCs and mean of logCPM, and the relationship is used to estimate the
#' variance of each mean logCPM value. Gene-wise weights, which are inverse of variances, are then
#' multiplied to single sample logFCs to downweight genes with low counts. It is assumed that the genes
#' with extremely low counts have been removed and the count matrix has been normalised prior to logCPM
#' matrix was derived. Rownames of the matrix must be genes' entrez ID. To convert other gene identifiers
#' to entrz ID, see example.
#'
#' Sample metadata should have the same number of rows as the number of columns in the logCPM matrix.
#' Metadata also must have a column called "sample" storing sample names (column names of logCPM matrix),
#' and a column called "treatment" storing treatment of each sample.The control treatment level specified
#' by `control` parameter must exist in the treatment column.
#'
#' This analysis was designed for experimental designs that include matched pairs of samples, such as when
#' tissues collected from the same patient were treated with different treatments to study different treatment
#' effects. Parameter `factor` tells the function how samples can be put into matching pairs. It must also be
#' included as a column in the metadata.
#'
#' @param logCPM Matrix of normaslised logCPM where rows are genes and columns are samples. Row names need to
#' be gene entrez ID and column names need to be sample names
#' @param metadata Sample metadata data frame as described in the details section.
#' @param factor Factor defines how samples can be put into matching pairs (eg. patient).
#' @param control Treatment level that is the control.
#'
#' @importFrom stats approxfun lowess var
#' @return A list with two elements:
#' $weight  gene-wise weights;
#' $logFC weighted single sample logFC matrix
#' @examples
#' require(AnnotationHub)
#' require(ensembldb)
#' # convert rownamews of logCPM from gene ids to gene entrez IDs through `AnnotationHub`
#' ah <- AnnotationHub()
#' ah <- subset(ah,genome == "GRCh38" & title == "Ensembl 101 EnsDb for Homo sapiens")
#' ensDb <- ah[[1]]
#' rownames(logCPM_example) <- mapIds(ensDb, rownames(logCPM_example), "ENTREZID", keytype = "GENEID")
#'
#' # Remove genes that couldn't be matched to entrez IDs
#' logCPM_example <- logCPM_example[!is.na(rownames(logCPM_example)),]
#'
#' # Inspect metadata data frame to make sure it has treatment, sample and patient columns
#' head(metadata_example)
#' length(setdiff(colnames(logCPM_example), metadata_example$sample)) == 0
#' ls <- weight_ssFC(logCPM_example, metadata = metadata_example,
#'  factor = "patient", control = "Vehicle")
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

    logCPM <- as.matrix(logCPM)
    metadata <- as.data.frame(metadata)

    # checks
    if (!all(c("treatment", "sample", factor) %in% colnames(metadata))) stop("Sample metadata must include factor, treatment and sample")
    if (any(c(!control %in% unique(metadata$treatment), length(unique(metadata[,"treatment"])) <2))) stop(
        "Treatment needs at least 2 levels where one is the control specified")
    if (!setequal(colnames(logCPM), metadata[,"sample"])) stop("Sample metadaata does not match with logCPM's column names")
    m <- min(logCPM)
    if (is.na(m)) stop("NA values not allowed")


    pairs <- unique(as.character(pull(metadata, sym(factor))))
    ls <- sapply(pairs, function(x){
    contrSample <- dplyr::filter(metadata, metadata$treatment == control, !!sym(factor) == x)
    contrSample <- as.character(pull(contrSample, sample))

    treatedSample <- dplyr::filter(metadata, metadata$treatment != control, !!sym(factor) == x)
    treatedSample <- as.character(pull(treatedSample, sample))

    logCPM[, treatedSample] - logCPM[, contrSample]
    }, simplify = FALSE)
    do.call(cbind,ls)

}


