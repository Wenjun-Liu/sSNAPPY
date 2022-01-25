
#' Title
#'
#' @param weightedFC
#' @param kegg_dir
#'
#' @import purrr set_names
#' @import plyr compact
#' @return
#' @export
#'
#' @examples
perturbationScore <- function(weightedFC, kegg_dir){

    BminsI <- readRDS(here::here(kegg_dir))
    #  remove pathway with 0 expressed gene in it
    kg2keep <- sapply(names(BminsI), function(x){
        length(intersect(rownames(weightedFC),
                         rownames(BminsI[[x]]))) >5
    })
    BminsI <- BminsI[kg2keep]
    if(length(BminsI) == 0) stop("None of the expressed gene was matched to pathways")

    weightedFC <- apply(weightedFC, 2, function(x){
        purrr::set_names(x, rownames(weightedFC))
    })

    PF <- sapply(names(BminsI), function(x){
        # if BminsI is invertible, solve for (B-I)PF = -delE
        if (abs(det(BminsI[[x]])) > 1e-7){
            sapply(colnames(weightedFC), function(y){
                delE  <- weightedFC[,y]
                delE  <- delE[rownames(BminsI[[x]])]
                # If any of the pathway gene was not expressed, set the ssLogFC to 0
                delE  <- replace(delE, is.na(delE), 0)
                solve(BminsI[[x]], -delE)
            })
        } else {
             NULL
        }
    }, simplify = FALSE)

    # Remove list elements that are null
    compact(PF)
}

perturbationScore_alt <- function(weightedFC, kegg_dir){

    BminsI <- readRDS(here::here(kegg_dir))
    #  remove pathway with 0 expressed gene in it
    kg2keep <- sapply(names(BminsI), function(x){
        length(intersect(rownames(weightedFC),
                         rownames(BminsI[[x]]))) >5
    })
    BminsI <- BminsI[kg2keep]
    if(length(BminsI) == 0) stop("None of the expressed gene was matched to pathways")

    weightedFC <- apply(weightedFC, 2, function(x){
        purrr::set_names(x, rownames(weightedFC))
    })

    .ssPertScore(weightedFC, BminsI)

    .ssPertScore <- function()
    PF <- sapply(names(BminsI), function(x){
        # if BminsI is invertible, solve for (B-I)PF = -delE
        if (abs(det(BminsI[[x]])) > 1e-7){
            sapply(colnames(weightedFC), function(y){
                delE  <- weightedFC[,y]
                delE  <- delE[rownames(BminsI[[x]])]
                # If any of the pathway gene was not expressed, set the ssLogFC to 0
                delE  <- replace(delE, is.na(delE), 0)
                solve(BminsI[[x]], -delE)
            })
        } else {
            NULL
        }
    }, simplify = FALSE)

    # Remove list elements that are null
    compact(PF)
}

rel<-c("activation","compound","binding/association","expression","inhibition",
       "activation_phosphorylation","phosphorylation","inhibition_phosphorylation",
       "inhibition_dephosphorylation","dissociation","dephosphorylation",
       "activation_dephosphorylation","state change","activation_indirect effect",
       "inhibition_ubiquination","ubiquination", "expression_indirect effect",
       "inhibition_indirect effect","repression","dissociation_phosphorylation",
       "indirect effect_phosphorylation","activation_binding/association",
       "indirect effect","activation_compound","activation_ubiquination")

#' Title Create weighted adjacent matrix
#'
#' @param list
#' @param beta
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
.weightedAdjMatrix <- function(list, beta, filename){

    if(is.null(beta)){
        beta=c(1,0,0,1,-1,1,0,-1,-1,0,0,1,0,1,-1,0,1,-1,-1,0,0,1,0,1,1)
        names(beta)<-rel
    }else{
        if(!all(names(beta) %in% rel) | length(names(beta))!=length(rel)){
            stop(paste("beta must be a numeric vector of length",length(rel), "with the following names:", "\n", paste(rel,collapse=",")))
        }
    }


    BminsI <- sapply(names(datpT), function(x){
        g2gInteraction <- sapply(rel, function(y){
            datpT[[x]][[y]] * beta[y]

        }, simplify = FALSE)
        g2gInteraction <-   Reduce('+', g2gInteraction)

        numDownstream <- sapply(rel, function(y){
            datpT[[x]][[y]] * abs(sign(beta[y]))
        }, simplify = FALSE)
        numDownstream <- Reduce('+', numDownstream)
        numDownstream <- apply(numDownstream, 2, sum)
        numDownstream[numDownstream==0]<-1
        B <- t(t(g2gInteraction)/numDownstream)
        diag(B) <- diag(B)-1
        B

    }, simplify = FALSE)

    saveRDS(BminsI, filename)


}




