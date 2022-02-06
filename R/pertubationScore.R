
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

    if (length(intersect(rownames(weightedFC), rownames(BminsI[[1]]))) == 0)
        stop("Weighted ssFCs and pathwy topologies must use the same gene identifiers.")

    #  remove pathway with 0 expressed gene in it
    kg2keep <- sapply(names(BminsI), function(x){
        length(intersect(rownames(weightedFC),
                         rownames(BminsI[[x]]))) > 0
    })
    BminsI <- BminsI[kg2keep]
    if(length(BminsI) == 0) stop("None of the expressed gene was matched to pathways")

    PF <- .ssPertScore(BminsI, weightedFC)

    # Remove list elements that are null or all zeros
    PF <- PF[!sapply(PF, is.null)]
    PF <- PF[sapply(PF, any)]

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
#' @param BminsI
#' @param weightedFC
#'
#' @return
.ssPertScore <- function(BminsI, weightedFC){
   sapply(names(BminsI), function(x){

        if (abs(det(BminsI[[x]])) > 1e-7){
             sapply(colnames(weightedFC), function(y){
                delE  <- weightedFC[,y]
                delE  <- delE[rownames(BminsI[[x]])]
                # If any of the pathway gene was not expressed, set the ssLogFC to 0
                delE  <- replace(delE, is.na(delE), 0)
                PF <- solve(BminsI[[x]], -delE)
                x <- sum(PF - delE)
            })
        } else {
        # if determinant of the pathway topology is not positive, the equation does not have a unique solution
            x <- NULL
        }

    }, simplify = FALSE)
}



