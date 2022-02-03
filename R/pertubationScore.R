
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

    if (missing(filePath) | !file.exists(filePath)) stop("Pathway topology matrices not detected. Check the file path provided.")

    BminsI <- readRDS(filePath)

    if (length(intersect(rownames(weightedFC), rownames(BminsI[[1]]))) == 0)
        stop("Weighted ssFCs and pathwy topologies must use the same gene identifiers")

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

BminsI <- BminsI %>%
    lapply(function(x){})

#' @title Create weighted adjacent matrix
#'
#' @description Convert pathway topology matrices to normalized weighted directed adjacency matrices describing the gene signaling network.
#'
#' @param species
#' @param database
#' @param pathwayName Optional.
#' @param beta Optional. A named numeric vector of weights to be assigned to each type of gene/protein relation type. See details for more information.
#' @param outputDir A file directory specifying where the output should be stored.
#'
#' @details
#'
#' This function takes the pathway topology information retrieved using `graphite` and convert them to normalized weighted directed adjacency
#' matrices describing the gene signaling network. See cited document for more details.
#'
#' The beta parameter speicifies weights to be assigned to each type of gene-gene interaction. It should be a named numeric vector of length 23,
#' whose names must be: c("activation","compound","binding/association","expression","inhibition","activation_phosphorylation","phosphorylation",
#' "indirect","inhibition_phosphorylation","dephosphorylation_inhibition","dissociation","dephosphorylation","activation_dephosphorylation",
#' "state","activation_indirect","inhibition_ubiquination","ubiquination","expression_indirect","indirect_inhibition",
#' "repression","binding/association_phosphorylation","dissociation_phosphorylation","indirect_phosphorylation"). If unspecified, beta will be by default
#' chosen as: c(1,0,0,1,-1,1,0,0,-1,-1,0,0,1,0,1,-1,0,1,-1,-1,0,0,0).
#'
#' @references Tarca AL, Draghici S, Khatri P, Hassan SS, Mittal P, Kim JS, Kim CJ, Kusanovic JP, Romero R. A novel signaling pathway impact analysis.
#' Bioinformatics. 2009 Jan 1;25(1):75-82.
#' @export
#'
#' @examples
weightedAdjMatrix <- function(species, database, pathwayName = NULL, beta = NULL, outputDir){

    supportedDatabase <- graphite::pathwayDatabases()
    if(any(c(!species %in% supportedDatabase$species,!database %in% supportedDatabase$database)))stop(
        "Requested species or database currently not supported by `grahpite`. Run `pathwayDatabses`
        to get databases available."
    )
    if(is.null(outputDir)) stop("Output directory must be specified")

    rel<-c("activation","compound","binding/association","expression","inhibition",
           "activation_phosphorylation","phosphorylation","inhibition_phosphorylation",
           "inhibition_dephosphorylation","dissociation","dephosphorylation",
           "activation_dephosphorylation","state change","activation_indirect effect",
           "inhibition_ubiquination","ubiquination", "expression_indirect effect",
           "inhibition_indirect effect","repression","dissociation_phosphorylation",
           "indirect effect_phosphorylation","activation_binding/association",
           "indirect effect","activation_compound","activation_ubiquination")

    if(is.null(beta)){
        beta=c(1,0,0,1,-1,1,0,-1,-1,0,0,1,0,1,-1,0,1,-1,-1,0,0,1,0,1,1)
        names(beta)<-rel
    }else{
        if(!all(names(beta) %in% rel) | length(names(beta))!=length(rel)){
            stop(paste("beta must be a numeric vector of length",length(rel), "with the following names:", "\n", paste(rel,collapse=",")))
        }
    }

    datpT <- .retrieveTopology(species, database, pathwayName)

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

    saveRDS(BminsI, file = outputDir)


}


#' Title
#'
#' @param species
#' @param database
#' @param pathwayName
#'
#' @importFrom graphite pathways convertIdentifiers prepareSPIA
#' @return
#' @examples
.retrieveTopology <- function(species, database, pathwayName = NULL){
    pys <- pathways(species, database)
    if(!is.null(pathwayName)){
        if (any(pathwayName %in% names(pys))){
            pys[names(pys) %in% pathwayName]
        } else stop("Pathway names provided not detected in retrieved database")
    }
    # always convert pathway nodes identifier to entrez ID
    pys <- suppressMessages(convertIdentifiers(pys, "ENTREZID"))
    # prepare the topologies for SPIA algorithm and store as a temporary file
    outputDir <- tempfile()
    prepareSPIA(pys, outputDir)
    # read the RData into the env
    get(load(paste(outputDir, "SPIA.RData", sep = "")))
}

