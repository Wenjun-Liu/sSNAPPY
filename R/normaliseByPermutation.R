#' @title Permute sample labels to simulate null distribution of perturbation scores
#'
#' @description Simulate null distributions of perturbation scores for each
#' pathway through sample permutation.
#'
#' @details
#' This \code{generate_permuted_scores} function is a generic function that can
#' deal with multiple types of inputs. It firstly randomly permute sample labels
#' to form permuted pairs and generate permuted logFCs, which are then used to
#' compute permuted perturbation scores for each pathway.
#'
#' The function outputs a list that is of the same length as the list storing
#' pathway topology matrices. Each element of the output list corresponds to one pathway
#' and contains a vector of permuted perturbation scores. The permuted perturbation
#' scores will be used to approximate the null distributions of perturbation scores.
#'
#'If the input is S4 object of \code{DGEList or SummarizedExperiment}, gene expression
#'matrix will be extracted and converted to a logCPM matrix.
#'
#' The number of unique permutations when choosing 2 samples from a dataset with
#' N samples without repetitions is N x (N-1). By default, permuted logFCs will
#' be computed for all N x (N-1) permuted pairs of samples. However, this can be
#' overwritten by the `NB` parameter. If the `NB` specified is smaller than N x (N-1),
#' `NB` possibilities will be randomly sampled from all the possible permutation.
#' If the `NB` specified is larger than the maximum number of permutation possible,
#' the parameter will be ignored.
#'
#' @param expreMatrix `matrix` and `data.frame` of logCPM, or `DGEList`/`SummarizedExperiment`
#' storing gene expression counts. Feature names need to be in entrez IDs
#' @param NB Number of permuted sample pairs to compute permuted logFCs for.
#' Default to be the maximum number of possibilities.
#' @param testScore Optional. Users can provide the test perturbation score
#' `data.frame` (ie. output of `pathwayPertScore`) to restrict the permutation
#' step only to pathways with non-zero test scores in at least one sample.
#' @param gsTopology List of pathway topology matrices generated using function `retrieve_topology`
#' @param weight A vector of gene-wise weights derived from function `weight_ss_fc`
#' @import Rcpp
#' @return A list where each element is a vector of perturbation scores for a pathway.
#' @examples
#' #compute weighted single sample logFCs
#' data(metadata_example)
#' data(logCPM_example)
#' metadata_example <- dplyr::mutate(metadata_example, treatment = factor(
#'    treatment, levels = c("Vehicle", "E2+R5020", "R5020")))
#' ls <- weight_ss_fc(logCPM_example, metadata = metadata_example,
#'  groupBy = "patient", treatColumn = "treatment", sampleColumn = "sample")
#' \dontrun{
#' load(system.file("extdata", "gsTopology.rda", package = "sSNAPPY"))
#'
#' # simulate the null distribution of scores through sample permutation
#' permutedScore <- generate_permuted_scores(logCPM_example,
#' gsTopology = gsTopology, weight = ls$weight)
#'
#'  }
#' @export
setGeneric("generate_permuted_scores",
           function(expreMatrix, NB = NULL, testScore = NULL, gsTopology, weight)
               standardGeneric("generate_permuted_scores"))

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "matrix"),
          function(expreMatrix, NB = NULL, testScore = NULL, gsTopology, weight){
              # checks

              m <- min(expreMatrix)
              if (is.na(m)) stop("NA values not allowed")

              if (nrow(expreMatrix) != length(weight))
                  stop("Gene-wise weights do not match with the dimension of expreMatrix")

              rownames(expreMatrix) <- paste("ENTREZID:", rownames(expreMatrix), sep = "")

              if (!is.null(testScore) && "gs_name" %in% colnames(testScore)){
                  gsTopology <- gsTopology[
                      names(gsTopology) %in% unique(testScore$gs_name)]
              }

              if (length(gsTopology) == 0) stop(
                  "None of the pathways had non-zero test perturbation scores"
              )

              if (
                  length(
                      intersect(rownames(expreMatrix),
                                unlist(unname(lapply(gsTopology, rownames)))
                                )) == 0)
                  stop("None of the expressed gene was matched to pathways. Check if gene identifiers match")

              # Exhaust all possible permutation combinations by default
              nSample <- ncol(expreMatrix)
              if (is.null(NB)){NB <- nSample * (nSample -1)}

              # set expression values and weights of unexpressed pathway genes to 0
              notExpressed <- setdiff(
                  unique(unlist(unname(lapply(gsTopology, rownames)))),
                  rownames(expreMatrix))
              if (length(notExpressed) != 0){
                  temp <- matrix(0, nrow = length(notExpressed), ncol = ncol(expreMatrix))
                  rownames(temp) <- notExpressed
                  colnames(temp) <- colnames(expreMatrix)
                  expreMatrix <- rbind(expreMatrix, temp)
                  weight <- c(weight, rep(0, length(notExpressed)))

              }

              # Generate permuted FCs for permuted pairs
              permutedFC <- .generate_permutedFC(expreMatrix, NB, weight)

              # Compute gene-level perturbation scores for each column of permuted
              # FCs
              permute_Gene <- GenePertScore_RCPP(
                  gsTopology, permutedFC, rownames(permutedFC)
                  )

              # Sum gene-wise permuted scores to pathway-wise
              permute_path <- lapply(permute_Gene, function(x){
                  apply(x, 2, sum)
              })

              names(permute_path) <- names(gsTopology)
              permute_path

          })

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "data.frame"),
          function(expreMatrix, NB = NULL, testScore = NULL, gsTopology, weight){
              generate_permuted_scores(as.matrix(expreMatrix), NB,
                                      testScore, gsTopology, weight)
          })

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "DGEList"),
          function(expreMatrix, NB = NULL, testScore = NULL, gsTopology, weight){
              expreMatrix <- cpm(expreMatrix$counts, log = TRUE)
              generate_permuted_scores(expreMatrix, NB,
                                       testScore, gsTopology, weight)
          })

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "SummarizedExperiment"),
          function(expreMatrix, NB = NULL, testScore = NULL, gsTopology, weight){
              expreMatrix <- cpm(SummarizedExperiment::assay(expreMatrix), log = TRUE)
              generate_permuted_scores(expreMatrix, NB,
                                       testScore, gsTopology, weight)
          })

#' @title Normalise test perturbation scores by permutation results
#'
#' @details Normalise the test perturbation scores generated by `weight_ss_fc()`
#' through the permuted perturbation scores derived from the
#' `generate_permuted_scores()` function. The mean absolute deviation(MAD) and
#' median of perturbation scores for each pathway are firstly derived from the
#' permuted perturbation scores. The test perturbation scores are then converted
#' to robust z-scores using MADs and medians calculated.
#' @param permutedScore A list. Output of `generate_permuted_scores`
#' @param testScore A `data.frame`. Output of `pathwayPertScore``
#' @param pAdj_method Method for adjusting p-values for multiple comparisons.
#' See `?p.adjust` for methods available. Default to FDR.
#' @param sortBy Sort the output by p-value, gene-set name or sample names.
#'
#' @importFrom stats mad median p.adjust pnorm p.adjust.methods
#' @importFrom dplyr left_join mutate
#' @return A `data.frame`
#' @export
#' @examples
#' \dontrun{
#' load(system.file("extdata", "gsTopology.rda", package = "sSNAPPY"))
#' data(metadata_example)
#' data(logCPM_example)
#' metadata_example <- dplyr::mutate(metadata_example, treatment = factor(
#'    treatment, levels = c("Vehicle", "E2+R5020", "R5020")))
#' ls <- weight_ss_fc(logCPM_example, metadata = metadata_example,
#' groupBy = "patient", sampleColumn = "sample", treatColumn = "treatment")
#'
#' # compute raw gene-wise perturbation scores
#' genePertScore <- raw_gene_pert(ls$logFC, gsTopology)
#'
#' # sum gene-wise perturbation scores to derive the pathway-level
#' # single-sample perturbation scores
#' pathwayPertScore <- pathway_pert(genePertScore)
#'
#' # simulate the null distribution of scores through sample permutation
#' permutedScore <- generate_permuted_scores(logCPM_example,
#' gsTopology = gsTopology, weight = ls$weight)
#'
#' # normlise the test perturbation scores using the permutation results
#' normalisedScores <- normalise_by_permu(permutedScore, pathwayPertScore,
#' sortBy = "pvalue")
#'  }
normalise_by_permu <- function(permutedScore, testScore,
                               pAdj_method = "fdr",
                               sortBy = c("gs_name", "sample", "pvalue")){
    pvalue <- NULL
    pAdj_method <- match.arg(pAdj_method, p.adjust.methods)
    sortBy <- match.arg(sortBy)

    summary_func <- function(x){c(MAD = mad(x), MEDIAN = median(x))}
    summaryScore <- t(sapply(permutedScore, summary_func))
    gs_name <- rownames(summaryScore)
    summaryScore <- mutate(as.data.frame(summaryScore), gs_name = gs_name)
    summaryScore <- filter(summaryScore, summaryScore$MAD != 0)
    summaryScore <- left_join(summaryScore, testScore, by = "gs_name",
                              multiple = "all")
    summaryScore <- mutate(
        summaryScore, robustZ =
            (summaryScore$score - summaryScore$MEDIAN)/summaryScore$MAD)
    summaryScore <- mutate(
        summaryScore, pvalue = 2*pnorm(-abs(summaryScore$robustZ)))
    summaryScore <- split(summaryScore, f = summaryScore$sample)
    summaryScore <- lapply(summaryScore, mutate,
                           adjPvalue = p.adjust(pvalue, pAdj_method))
    summaryScore <-bind_rows(summaryScore)
    summaryScore <- dplyr::mutate_at(
        summaryScore, vars(c("gs_name", "sample")), as.factor
    )
    summaryScore[
        order(pull(summaryScore, sym(sortBy)), decreasing = FALSE),]

}

#' @importFrom gtools permutations
.generate_permutedFC <- function(expreMatrix, NB, weight){

    colNumber <- seq_len(ncol(expreMatrix))

    # generate all permutation pairs
    Pairs <- gtools::permutations(
        n = length(colNumber), r = 2, v = colNumber)

    if (NB < nrow(Pairs)){
        Pairs <- Pairs[sample(nrow(Pairs), NB, replace = FALSE),]}

    apply(Pairs, 1, function(x){
        (expreMatrix[,x[1]] -
             expreMatrix[,x[2]]) * weight
    })
}



