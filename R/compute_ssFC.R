#' @title Compute single sample logFCs
#'
#' @description Compute single sample logFCs for each treated samples using normalised logCPM values
#'
#' @details The number of rows must equal to the number of columns in the logCPM matrix. It also must have a column called "sample" storing sample names, and a column called "treatment" storing treatment of each sample.
#'
#' @param logCPM Matrix of normaslised logCPM
#' @param metadata Sample metadata data frame as described in the details section.
#' @param factor The factor defines how samples can be put into matching pairs.
#' @param control The treatment level that is the control.
#'
#' @return A matrix of single sample logFC
#'
#' @export
#'
compute_ssFC <- function(logCPM, metadata, factor, control){

    ## checks
    # if (is.null(logCPM)) stop("LogCPM has to be provided")
    # if (is.null(metadata)) stop("Sample metadata has to be provided")
    if (missing(factor)) stop("Factor defining matching samples must be provided")
    if (missing(control)) stop("Control treatment must be specified")
    if (!"treatment" %in% colnames(metadata)) stop("Sample metadata must contain a column named *treatment*")
    if (!control %in% unique(metadata$treatment)) stop("Control level not detected in sample metadata")
    if ("sample" %in% colnames(metadata)){
        if (!setequal(metadata[,"sample"], colnames(logCPM))) stop("Sample metadata does not match with logCPM")
    } else stop ("Sample name must be specific in a column named *sample*")
    stopifnot(factor %in% colnames(metadata))
    stopifnot(ncol(logCPM) == nrow(metadata))

    ##
    logCPM <- as.matrix(logCPM)
    metadata <- as.data.frame(metadata)
    if(length(unique(metadata[,"treatment"])) <2) stop("At least 2 levels are required for the treatment column")
    pairs <- unique(metadata[,factor])
    sapply(pairs, function(x){
       contrSample <- metadata %>%
           dplyr::filter(
               treatment == control,
               !!sym(factor) == x
           ) %>%
           pull(sample)
       treatedSample <- metadata %>%
           dplyr::filter(
               treatment != control,
               !!sym(factor) ==x
           ) %>%
           pull(sample)
       logCPM[, treatedSample] - logCPM[, contrSample] %>%
           head()
    }, simplify = FALSE) %>%
        do.call(cbind,.)

}

