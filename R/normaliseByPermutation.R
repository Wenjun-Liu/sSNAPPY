#' Title
#'
#' @param logCPM
#' @param metadata
#' @param factor
#' @param control
#' @param seed
#' @param NB
#' @param scores
#' @param filPath
#' @param weight
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join filter mutate select
#' @return
#' @export
#'
#' @examples
normaliseByPermutation <- function(logCPM, metadata, factor, control,
                                   seed = sample.int(1e+06, 1), NB = 1000,
                                   scores, filePath, weight){
    # checks
    if (missing(filePath) | !file.exists(filePath)) stop("Pathway topology matrices not detected")
    if (missing(factor)) stop("Factor defining matching samples must be provided")
    if (missing(control)) stop("Control treatment must be specified")
    if (!"treatment" %in% colnames(metadata)) stop("Sample metadata must contain a column named treatment")
    if (!control %in% unique(metadata$treatment)) stop("Control level not detected in sample metadata")
    if (!"sample" %in% colnames(metadata)) stop ("Sample name must be specific in a column named sample")
    stopifnot(factor %in% colnames(metadata))
    stopifnot(ncol(logCPM) == nrow(metadata))
    m <- min(logCPM)
    if (is.na(m)) stop("NA values not allowed")

#     # set BPPARAM
#     if(is.null(BPPARAM)){
#         BPPARAM <- BiocParallel::registered()[[1]]
#     }
#     BPPARAM$workers <- cores

    # load pathway topologies

    BminsI <- readRDS(filePath)

    # if gene-wise weights are not provided, estimate again
    if (is.null(weight))stop("Gene-wise weight must be provided. See details.")

    permutedFC <- .generate_permutedFC(logCPM, metadata, factor, control, weight, NB, seed)

    # Remove pathways with 0 expressed genes in it
    kg2keep <- sapply(names(BminsI), function(x){
        length(intersect(rownames(permutedFC[[1]]),
                         rownames(BminsI[[x]]))) > 0
    })
    BminsI <- BminsI[kg2keep]
    if(length(BminsI) == 0) stop("None of the expressed gene was matched to pathways")

    # compute permuted perturbation scores and remove pathways returned to be NULL

    permutedScore <- lapply(permutedFC, .ssPertScore, BminsI = BminsI)
    permutedScore <- do.call(mapply, c(FUN=c, lapply(permutedScore, `[`, names(BminsI))))
    permutedScore <- permutedScore[!sapply(permutedScore, is.null)]

    summary_func <- function(x){c(MAD = mad(x), MEDIAN = median(x))}
    summaryScore <- as.data.frame(t(sapply(permutedScore, summary_func)))
    summaryScore <- rownames_to_column(summaryScore,"gs_name")
    summaryScore <- filter(summaryScore, MAD != 0)
    summaryScore <- left_join(summaryScore, scores, by = "gs_name")
    mutate(summaryScore, robustZ = (tA - MEDIAN)/MAD )

}



#' Title
#'
#' @param logCPM
#' @param metadata
#' @param factor
#' @param control
#' @param BPPARAM
#' @param seed
#'
#' @return
#'
#' @examples
.generate_permutedFC <- function(logCPM, metadata, factor, control, weight, NB, seed){

    logCPM <- as.matrix(logCPM)
    rownames(logCPM) <- paste("ENTREZID:", rownames(logCPM), sep = "")

    metadata <- as.data.frame(metadata)
    pairs <- unique(metadata[,factor])
    sampleInpairs <- sapply(pairs, function(x){
        contrSample <- dplyr::filter(metadata, treatment == control, !!sym(factor) == x)
        contrSample <- pull(contrSample, sample)

        treatedSample <- dplyr::filter(metadata, treatment != control, !!sym(factor) == x)
        treatedSample <- pull(treatedSample, sample)
        list(contrSample = contrSample, treatedSample = treatedSample)
    }, simplify = FALSE)

    nSample <- nrow(metadata)
    NB <- min(factorial(nSample), NB)
    set.seed(seed)
    sapply(1:NB, function(x){
        # permute sample labels to get permuted logCPM
        colnames(logCPM) <- sample(colnames(logCPM), ncol(logCPM))

        # Built permuted logFCs based on the permuted logCPM
        permutedFC <- sapply(names(sampleInpairs), function(y){
            logCPM[, sampleInpairs[[y]]$treatedSample] - logCPM[, sampleInpairs[[y]]$contrSample]
        }, simplify = FALSE)  %>%
            do.call(cbind,.)

        # Multiply permuted FCs by gene-wise weights
        permutedFC * weight

    }, simplify = FALSE)

}
