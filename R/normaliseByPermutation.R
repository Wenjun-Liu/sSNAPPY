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
    if (!file.exists(filePath)) stop("Pathway topology matrices not detected")

    if (!all(c("treatment", "sample", factor) %in% colnames(metadata))) stop("Sample metadata must include factor, treatment and sample")
    if (any(c(!control %in% unique(metadata$treatment), length(unique(metadata[,"treatment"])) <2))) stop(
        "Treatment needs at least 2 levels where one is the control specified")
    if (ncol(logCPM) != nrow(metadata)) stop("Sample metadaata does not match with logCPM's dimension")
    m <- min(logCPM)
    if (is.na(m)) stop("NA values not allowed")

#     # set BPPARAM
#     if(is.null(BPPARAM)){
#         BPPARAM <- BiocParallel::registered()[[1]]
#     }
#     BPPARAM$workers <- cores

    # load pathway topologies
    BminsI <- readRDS(filePath)
    logCPM <- as.matrix(logCPM)
    rownames(logCPM) <- paste("ENTREZID:", rownames(logCPM), sep = "")
    if (length(intersect(rownames(logCPM), unlist(unname(lapply(BminsI, rownames))))) == 0)
        stop("None of the expressed gene was matched to pathways. Check if gene identifiers match")

    # if gene-wise weights are not provided, estimate again
    if (is.null(weight)) stop("Gene-wise weight must be provided. See details.")

    notExpressed <- setdiff(unique(unlist(unname(lapply(BminsI, rownames)))), rownames(logCPM))
    if (length(notExpressed) != 0){
        temp <- matrix(0, nrow = length(notExpressed), ncol = ncol(logCPM))
        rownames(temp) <- notExpressed
        colnames(temp) <- colnames(logCPM)
        rbind(temp, logCPM)

    }

    # this step will probably benefit from BiocParallel but my laptop freezes with 2 workers
    permutedFC <- .generate_permutedFC(logCPM, metadata, factor, control, weight, NB, seed)

    # Remove pathways with 0 expressed genes in it
    kg2keep <- sapply(names(BminsI), function(x){
        length(intersect(rownames(permutedFC[[1]]),
                         rownames(BminsI[[x]]))) > 0
    })
    BminsI <- BminsI[kg2keep]

    # compute permuted perturbation scores and remove pathways returned to be NULL

    permutedScore <- lapply(permutedFC, ssPertScore_RCPP, BminsI = BminsI)
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
        }, simplify = FALSE)
        permutedFC <- do.call(cbind,permutedFC)

        # Multiply permuted FCs by gene-wise weights
        permutedFC * weight

    }, simplify = FALSE)

}


#' @return
.permutedFC_parallel <- function(logCPM, metadata, factor, control, weight, NB, seed){
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

    #     # set BPPARAM
    #     if(is.null(BPPARAM)){
    BPPARAM <- BiocParallel::registered()[[1]]
    #     }
    #     BPPARAM$workers <- cores
    BiocParallel::bplapply(1:NB, function(x){
        # permute sample labels to get permuted logCPM
        colnames(logCPM) <- sample(colnames(logCPM), ncol(logCPM))

        # Built permuted logFCs based on the permuted logCPM
        permutedFC <- sapply(names(sampleInpairs), function(y){
            logCPM[, sampleInpairs[[y]]$treatedSample] - logCPM[, sampleInpairs[[y]]$contrSample]
        }, simplify = FALSE)
        permutedFC <- do.call(cbind,permutedFC)

        # Multiply permuted FCs by gene-wise weights
        permutedFC * weight

    }, BPPARAM = BPPARAM)

}

