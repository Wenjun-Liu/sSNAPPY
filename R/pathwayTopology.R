#' @title Retrieve pathway topology as weighted adjacency matrix
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
#' This function takes pathway topology information retrieved using `graphite`
#' and converts these to normalized weighted directed adjacency matrices
#' describing the gene signaling network, which can be used to compute gene-wise
#' and pathway-level perturbation score through the scoring algorithm derived
#' from the *SPIA* algorithm. See cited document for more details.
#'
#' The `database` parameter may specify multiple databases but only one species
#' can provided to the `species` parameter.
#'
#' Users can provide areas of interests as keywords, which will be matched to
#' pathway names from chosen databases to subset pathways. Cases will be
#' ignored in string matchings.
#'
#' The beta parameter specifies weights to be assigned to each type of gene-gene
#' interaction. It should be a named numeric vector of length 25,
#' whose names must be: `c("activation","compound","binding/association",
#' "expression","inhibition","activation_phosphorylation","phosphorylation",
#' "indirect","inhibition_phosphorylation","dephosphorylation_inhibition",
#' "dissociation","dephosphorylation","activation_dephosphorylation",
#' "state","activation_indirect","inhibition_ubiquination","ubiquination",
#' "expression_indirect","indirect_inhibition","repression",
#' "binding/association_phosphorylation","dissociation_phosphorylation",
#' "indirect_phosphorylation")`.
#'
#' If unspecified, beta will be set as an integer vector with: a) values of 1
#' for interactions which match 'expression' or 'activation'; b) values of -1
#' for interactions which match 'repression' or 'inhibition'; and c) 0
#' elsewhere.
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
#' # If only interested in selected pathways, specify the areas of interest as
#' # keywords
#' gsTopology <- retrieve_topology(database = "kegg",
#' keyword = c("metabolism", "signaling"), species = "hsapiens")
#' }
retrieve_topology <-  function(
    database = c(
      "kegg", "pathbank", "wikipathways", "reactome", "panther",
      "pharmgkb", "smpdb"
    ),
    species = c(
      "hsapiens", "athaliana", "btaurus", "celegans", "cfamiliaris",
      "dmelanogaster", "drerio", "ecoli", "ggallus", "mmusculus",
      "rnorvegicus", "scerevisiae", "sscrofa", "xlaevis"
    ),
    keyword = NULL, beta = NULL
){
  ## Whatever I have done in this or the below function has removed the prefix
  ## 'database.' from the pathways. This breaks the vignette due to the
  ## pre-loaded nrmalised scores. We need to discuss
  database <- match.arg(database, several.ok = TRUE)
  species <- match.arg(species)

  rel <- c(
    "activation", "compound", "binding/association", "expression",
    "inhibition", "activation_phosphorylation", "phosphorylation",
    "inhibition_phosphorylation", "inhibition_dephosphorylation",
    "dissociation", "dephosphorylation", "activation_dephosphorylation",
    "state change", "activation_indirect effect", "inhibition_ubiquination",
    "ubiquination", "expression_indirect effect",
    "inhibition_indirect effect", "repression","dissociation_phosphorylation",
    "indirect effect_phosphorylation","activation_binding/association",
    "indirect effect","activation_compound","activation_ubiquination"
  )

  if (is.null(beta)) {
    beta <- rep(0, length(rel))
    names(beta) <- rel
    beta[grepl("inhibition|repression", rel)] <- -1
    beta[grepl("activation|expression", rel)] <- 1
  } else {
    if (!all(names(beta) %in% rel) | length(names(beta)) != length(rel)) {
      stop("Beta has wrong length or names. See details for requirements")
    }
  }

  datpT <- .retrieveTopology(database, species, keyword)
  int2keep <- names(beta[beta != 0])
  datpT <- lapply(datpT, function(x) x[names(x) %in% int2keep] )

  ## sometimes SPIA adds a list element called `<graphite_placeholder>`.
  ## Remove it if it exists
  datpT <- datpT[str_detect(names(datpT), "placeholder",  negate = TRUE) ]

  gsTopology <- lapply(
    names(datpT),
    function(x){
      g2gInteraction <- lapply(
        int2keep, function(y) datpT[[x]][[y]] * beta[y]
      )
      g2gInteraction <-   Reduce('+', g2gInteraction)
      numDownstream <- apply(g2gInteraction , 1, function(x){sum(x!=0)})
      numDownstream[numDownstream == 0] <- 1
      B <- g2gInteraction/numDownstream
      diag(B) <- diag(B) - 1
      B
    }
  )
  names(gsTopology) <- names(datpT)
  gsTopology

}

#' @importFrom graphite pathways convertIdentifiers prepareSPIA
.retrieveTopology <- function(database, species, keyword){

  pys <- lapply(database, function(x) pathways(species, x))
  names(pys) <- database

  # If keywords are provided, use pattern match to extract pathways whose
  # names match with the keywords
  if (!is.null(keyword)) {
    pys <- lapply(
      pys,
      function(x){
        pat <- paste(tolower(keyword), collapse = "|")
        ind <- grepl(pat, tolower(names(x)))
        x[ind]
      }
    )
    if (all(vapply(pys, length, integer(1)) == 0))
      stop("Keywords provided not detected in retrieved database")
  }
  # always convert pathway nodes identifier to entrez ID
  pys <- lapply(pys, convertIdentifiers, "ENTREZID")

  # prepare the topologies for SPIA algorithm and store as a temporary file
  outputDir <- tempdir()
  lapply(database, function(x) prepareSPIA(pys[[x]], file.path(outputDir, x)))

  # read the RData into the env
  ls <- lapply(
    database,
    function(x){
      rd_name <- file.path(outputDir, paste0(x, "SPIA.RData"))
      get(load(rd_name))
    }
  )
  names(ls) <- database
  unlist(ls, recursive = FALSE)

}
