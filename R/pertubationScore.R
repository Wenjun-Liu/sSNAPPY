#' @title Compute Single Sample Perturbation Score
#'
#' @description Propagate weighted single sample logFCs down the pathway topologies to compute single sample perturbation scores for each pathway
#'
#' @details This function use the algorithm adopted from `SPIA` (see citation) to compute a single sample perturbation score per sample per
#' pathway. The rownames of the weighted single sample logFC matrix and the pathway topology matrices must use the same type of gene identifier (ie. entrez ID).
#'
#' @param weightedFC A matrix of weighted single sample logFCs derived from function `weight_ssFC()`
#' @param gsTopology List of pathway topology matrices generated using function `weightedAdjMatrix()`
#'
#' @importFrom purrr set_names
#' @importFrom plyr compact
#' @importFrom dplyr bind_rows mutate
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr set_colnames
#'
#' @references Tarca AL, Draghici S, Khatri P, Hassan SS, Mittal P, Kim JS, Kim CJ, Kusanovic JP, Romero R. A novel signaling pathway impact analysis.
#' Bioinformatics. 2009 Jan 1;25(1):75-82.
#' @return A list where each element is a matrix corresponding to a pathway. Each column of an element corresponds to a sample.
#' @examples
#' #compute weighted single sample logFCs
#' ls <- weight_ssFC(logCPM_example, metadata = metadata_example,
#' factor = "patient", control = "Vehicle")
#'
#' # explore all species and databases supported by graphite
#' \dontrun{
#' graphite::pathwayDatabases()
#' weightedAdjMatrix(species = "hsapiens", database = "kegg",
#' outputDir = "gsTopology.rda")
#' }
#' load(system.file("extdata", "gsTopology.rda", package = "sSNAPPY"))
#' ssPertScore <- perturbationScore(ls$logFC, gsTopology)
#' @export
perturbationScore <- function(weightedFC, gsTopology){
    if (length(intersect(rownames(weightedFC), unlist(unname(lapply(gsTopology, rownames))))) == 0)
        stop("None of the expressed gene was matched to pathways. Check if gene identifiers match")
    # extract all unique pathway genes and find ones that are not expressed
    notExpressed <- setdiff(unique(unlist(unname(lapply(gsTopology, rownames)))), rownames(weightedFC))
    if (length(notExpressed) != 0){
        # set the FCs of unexpressed pathway genes to 0
        temp <- matrix(0, nrow = length(notExpressed), ncol = ncol(weightedFC))
        rownames(temp) <- notExpressed
        colnames(temp) <- colnames(weightedFC)
        # set the weights of unexpressed pathway genes to 0
        weightedFC <- rbind(weightedFC, temp)}

    PF <-  ssPertScore_RCPP(gsTopology, weightedFC, rownames(weightedFC), colnames(weightedFC))

    # Remove list elements that are null or all zeros
    suppressWarnings(PF <- PF[sapply(PF, any)])
    PF <- sapply(names(PF), function(x){
        temp <- as.data.frame(PF[[x]])
        temp <- set_colnames(temp, "tA")
        temp <- rownames_to_column(temp,"sample")
        temp <- mutate(temp, gs_name = x)
    }, simplify = FALSE)
    bind_rows(PF)
}
