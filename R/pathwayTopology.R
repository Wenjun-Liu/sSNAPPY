#' @title Retrieve pathway topology as weighted adjacent matrix
#'
#' @description Retrieve pathway topology matrices and convert them to normalized
#' weighted directed adjacency matrices describing gene signaling networks.
#'
#' @param database A character vector of supported databases.
#' @param species One of the supported species,
#' @param keyword Optional. Areas of interests as a character vector.
#' @param beta Optional. A named numeric vector of weights to be assigned to
#' each type of gene/protein relation type.See details for more information.
#' @import org.Hs.eg.db
#' @return A list where each element is a matrix corresponding to a pathway
#' @details
#' This function takes the pathway topology information retrieved using `graphite`
#' and convert them to normalized weighted directed adjacency matrices describing
#' the gene signaling network, which can be used to compute gene-wise and
#' pathway-level perturbation score through the scoring algorithm derived from
#' the *SPIA* algorithm. See cited document for more details.
#'
#' The `databse` parameter could be a vector containing multiple databases but
#' only one species name could be provided to the `species` parameter.
#'
#' Users can provide areas of interests as keywords, which will be matched to
#' pathway names from chosen databases to filter for pathways. Cases will be
#' ignored in string matchings.
#'
#' The beta parameter specifies weights to be assigned to each type of gene-gene
#' interaction. It should be a named numeric vector of length 23,
#' whose names must be: `c("activation","compound","binding/association",
#' "expression","inhibition","activation_phosphorylation","phosphorylation",
#' "indirect","inhibition_phosphorylation","dephosphorylation_inhibition",
#' "dissociation","dephosphorylation","activation_dephosphorylation",
#' "state","activation_indirect","inhibition_ubiquination","ubiquination",
#' "expression_indirect","indirect_inhibition","repression",
#' "binding/association_phosphorylation","dissociation_phosphorylation",
#' "indirect_phosphorylation")`.
#'
#' If unspecified, beta will be by default chosen as: `c(1,0,0,1,-1,1,0,0,-1,-1,
#' 0,0,1,0,1,-1,0,1,-1,-1,0,0,0)`.
#'
#' The converted weighted adjacent matrices will be stored in a list. We
#' recommend users to store the returned list as a file so this step only
#' needs to be performed once for each database.
#'
#' @references Tarca AL, Draghici S, Khatri P, Hassan SS, Mittal P, Kim JS,
#' Kim CJ, Kusanovic JP, Romero R. A novel signaling pathway impact analysis.
#' Bioinformatics. 2009 Jan 1;25(1):75-82.
#' Sales, G., Calura, E., Cavalieri, D. et al. graphite - a Bioconductor package
#' to convert pathway topology to gene network. BMC Bioinformatics 13, 20 (2012).
#' @export
#' @importFrom stringr str_detect
#' @examples
#' \donttest{
#' # retrieve pathway topology matrices of all KEGG pathway
#' gsTopology <- retrieve_topology(database = "kegg", species = "hsapiens")
#'
#' # retrieve pathway topology matrices from multiple databases
#' gsTopology <- retrieve_topology(database = c("kegg","reactome"),
#' species = "hsapiens")}
#'
#' # If only interested in selected pathways, specify the areas of interest as
#' # keywords
#' gsTopology <- retrieve_topology(database = "kegg",
#' keyword = c("metabolism", "signaling"), species = "hsapiens")
retrieve_topology <-  function(
        database = c("kegg", "pathbank", "wikipathways", "reactome", "panther",
                     "pharmgkb", "smpdb"),
        species = c("athaliana", "btaurus", "celegans", "cfamiliaris",
                    "dmelanogaster", "drerio", "ecoli", "ggallus", "hsapiens",
                    "mmusculus", "rnorvegicus", "scerevisiae", "sscrofa", "xlaevis"),
        keyword = NULL, beta = NULL){

    database <- match.arg(database, several.ok = TRUE)
    species <- match.arg(species)

    rel<-c("activation","compound","binding/association","expression","inhibition",
           "activation_phosphorylation","phosphorylation","inhibition_phosphorylation",
           "inhibition_dephosphorylation","dissociation","dephosphorylation",
           "activation_dephosphorylation","state change","activation_indirect effect",
           "inhibition_ubiquination","ubiquination", "expression_indirect effect",
           "inhibition_indirect effect","repression","dissociation_phosphorylation",
           "indirect effect_phosphorylation","activation_binding/association",
           "indirect effect","activation_compound","activation_ubiquination")

    if(is.null(beta)){
        beta <- c(1,0,0,1,-1,1,0,-1,-1,0,0,1,0,1,-1,0,1,-1,-1,0,0,1,0,1,1)
        names(beta)<-rel
    }else{
        if(!all(names(beta) %in% rel) | length(names(beta))!=length(rel)){
            stop("Beta has wrong length or names. See details for requirements")}}

    datpT <- .retrieveTopology(database, species, keyword)
    int2keep <- names(beta[beta !=0])
    datpT <- lapply(datpT, function(x) x[names(x) %in% int2keep] )

  #sometimes SPIA add a list element called `<graphite_placeholder>`.
  # Remove it if it exists
    datpT <- datpT[str_detect(names(datpT), "placeholder",  negate = TRUE) ]

    gsTopology <- lapply(names(datpT), function(x){
        g2gInteraction <- lapply(int2keep, function(y){
            datpT[[x]][[y]] * beta[y]})
        g2gInteraction <-   Reduce('+', g2gInteraction)

        numDownstream <- Reduce('+', datpT[[x]])
        numDownstream <- apply(numDownstream, 2, sum)
        numDownstream[numDownstream==0]<-1
        B <- g2gInteraction/numDownstream
        diag(B) <- diag(B)-1
        B
    })
    names(gsTopology) <- names(datpT)
    gsTopology

}

#' @importFrom graphite pathways convertIdentifiers prepareSPIA
#' @importFrom stringr str_subset regex
.retrieveTopology <- function(database, species, keyword){
    . <- NULL
    pys <- sapply(database, pathways, species = species, simplify = FALSE)

    # If keywords are provided, use pattern match to extract pathways whose
    # names match with the keywords
    if(!is.null(keyword)){
        pys <- lapply(pys, function(x){
            x %>%
                .[str_detect(
                    names(.), regex(
                        paste(keyword, collapse = "|"), ignore_case = T
                    )
                )]
        })
        if (
            all(
                sapply(pys, function(x)length(x) == 0))
        ) stop("Keywords provided not detected in retrieved database")
    }
    # always convert pathway nodes identifier to entrez ID
    pys <- lapply(pys, convertIdentifiers,  "ENTREZID")

    # prepare the topologies for SPIA algorithm and store as a temporary file
    outputDir <- tempdir()
    lapply(database, function(x){
        prepareSPIA(pys[[x]], paste(outputDir, x, sep = "/"))
    })

    # read the RData into the env
    ls <- sapply(database, function(x){
        get(load(paste(outputDir, "/",x, "SPIA.RData", sep = "")))
    }, simplify = FALSE)
    unlist(ls, recursive = FALSE)

}
