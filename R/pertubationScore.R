
#' @title Single Sample Perturbation Score
#'
#' @description Propagate weighted single sample logFCs down the pathway topologies to compute single sample perturbation scores for each pathway
#'
#' @details This function use the algorithm adopted from `SPIA` (see citation) to compute a single sample perturbation score per sample per
#' pathway. The rownames of the weighted single sample logFC matrix and the pathway toplogy matrices must use the same type of gene identifier.
#'
#' @param weightedFC A matrix of weighted single sample logFCs derived from function `weight_ssFC`
#' @param filePath The file path to pathway topology matrices generated using function `weightedAdjMatrix`
#'
#' @importFrom purrr set_names
#' @importFrom plyr compact
#' @importFrom dplyr bind_rows
#'
#' @references Tarca AL, Draghici S, Khatri P, Hassan SS, Mittal P, Kim JS, Kim CJ, Kusanovic JP, Romero R. A novel signaling pathway impact analysis.
#' Bioinformatics. 2009 Jan 1;25(1):75-82.
#' @return A list where each element is a matrix corresponding to a pathway. Each column of an element corresponds to a sample.
#' @export
perturbationScore <- function(weightedFC, filePath){

    if ( !file.exists(filePath)) stop("Pathway topology matrices not detected in the specified file path. Check the file path provided.")

    BminsI <- readRDS(filePath)

    if (length(intersect(rownames(weightedFC), unlist(unname(lapply(BminsI, rownames))))) == 0)
        stop("None of the expressed gene was matched to pathways. Check if gene identifiers match")

    # extract all unique pathway genes and find ones that are not expressed
    notExpressed <- setdiff(unique(unlist(unname(lapply(BminsI, rownames)))), rownames(weightedFC))
    if (length(notExpressed) != 0){
        temp <- matrix(0, nrow = length(notExpressed), ncol = ncol(weightedFC))
        rownames(temp) <- notExpressed
        colnames(temp) <- colnames(weightedFC)
        weightedFC <- rbind(weightedFC, temp)}

    PF <-  ssPertScore_RCPP(BminsI, weightedFC)

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


#' Title
#'
#' @param adjMatrix
#' @param weightedFC
#'
#' @return
.ssPertScore <- function(adjMatrix, weightedFC, tol = 1e-7){

    # if pathway adjacency matrix is not invertible, output NULL
    d <- abs(det(adjMatrix))
    if (d < tol) return(NULL)

    # subset pathway genes' expression
    x <- weightedFC[rownames(adjMatrix), ]

    apply(x, 2, function(y)sum(.Internal(La_solve(adjMatrix, -y, tol)) - y))

}

