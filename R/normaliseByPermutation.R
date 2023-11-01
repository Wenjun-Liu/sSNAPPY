#' @title Permute sample labels to simulate null distribution of perturbation scores
#'
#' @description Simulate null distributions of perturbation scores for each
#' pathway through sample permutation.
#'
#' @details
#' This \code{generate_permuted_scores} function firstly randomly permute sample labels
#' to form permuted pairs and generate permuted logFCs, which are then used to
#' compute permuted perturbation scores for each pathway.
#'
#' The function outputs a list that is of the same length as the list storing
#' pathway topology matrices (i.e. `gsTopology`). Each element of the output list
#' corresponds to one pathway and contains a vector of permuted perturbation scores.
#' The permuted perturbation scores will be used to approximate the null distributions
#' of perturbation scores and compute permuted p-values.
#'
#'If the input is S4 object of \code{DGEList or SummarizedExperiment}, gene expression
#'matrix will be extracted and converted to a logCPM matrix.
#'
#' The sample permutation is operated by randomly choosing combination of 2 samples
#' from the column names of the `expreMatrix`. Hence, sample metadata is not a
#' required input. The number of maximum unique permutations in a dataset with N
#' samples is N x (N-1). By default, permuted logFCs will be computed for all N x (N-1)
#' permuted pairs. However, this can be overwritten by the `NB` parameter.
#' If the `NB` specified is smaller than N x (N-1), `NB` possibilities will be randomly
#' sampled from all the possible permutation. If the `NB` specified is larger than
#' the maximum number of permutation possible, the parameter will be ignored.
#' Since the smallest achievable permutation p-value is 1/`NB`, accurate estimation
#' of small p-value requires a large number of permutations that is only feasible
#' in data with a large sample size.
#'
#' @param expreMatrix `matrix` and `data.frame` of logCPM, or `DGEList`/`SummarizedExperiment`
#' storing gene expression counts. Feature names need to be in entrez IDs
#' @param NB Number of permuted sample pairs to compute permuted logFCs for.
#' Default to be the maximum number of possibilities (See details).
#' @param testScore Optional. Output of `pathway_pert()` for restricting the
#' permutation only to pathways with non-zero test scores in at least one sample.
#' @param gsTopology List of pathway topology matrices generated using function `retrieve_topology`
#' @param weight A vector of gene-wise weights derived from function `weight_ss_fc`
#' @param drop logic(1). Whether to drop pathways with all zero scores
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
           function(
        expreMatrix, NB = NULL, testScore = NULL, gsTopology, weight, drop = TRUE)
               standardGeneric("generate_permuted_scores"))

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "matrix"),
          function(
        expreMatrix, NB = NULL, testScore = NULL, gsTopology, weight, drop = TRUE){
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

              # Exhaust all possible permutation combinations by default. If NB
              # is specified and greater than the maximum number of permutations
              # possible, change NB to be the maxNB
              nSample <- ncol(expreMatrix)
              maxNB <- nSample * (nSample -1)
              if (is.null(NB)){
                  NB <- maxNB
              } else {
                  if (NB > maxNB){NB <- maxNB}
              }

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

              # Generate a mtrix of permuted FCs for permuted pairs
              permutedFC <- .generate_permutedFC(expreMatrix, NB, weight)

              # Compute gene-level perturbation scores for each column of permuted
              # FCs

              allGene <- rownames(permutedFC)

              permute_Gene <- lapply(gsTopology, function(x){
                  #gs_sub <- x[rownames(x) %in% allGene, colnames(x) %in% allGene ]

                  # the matrix system only has unique solution when the det
                  # is non-zero. A buffer is set cause det smaller than that
                  # is considered as unsolveable by R too
                  if (abs(det(x))>1e-7){
                      # Loop through each column of the permutedFC matrix
                      netP_ls <- lapply(seq_len(ncol(permutedFC)), function(y){
                          # subset the FC vector to only pathway genes
                          de <- permutedFC[rownames(permutedFC) %in% rownames(x), y]
                          # make sure the row names of the topology matrix matches with
                          # the names of the FC vector
                          de <- de[match(rownames(x), names(de))]
                          # Solve for gene-wise pert score and substract FCs
                          solve(t(x), -de) - de
                      })
                      do.call(cbind,netP_ls)
                  } else {
                      NULL
                  }
              })

              names(permute_Gene) <- names(gsTopology)

              # remove the list element for pathways with linearly dependent
              # topology matrices
              permute_Gene <-  permute_Gene[!sapply( permute_Gene, is.null)]

              # Sum gene-wise permuted scores to pathway-wise permtued scores
              permute_path <- lapply(permute_Gene, function(x){
                  apply(x, 2, sum)
              })

              # if drop parameter is set to TRUE, drop pathways with all zero scores
              if (drop){
                  permute_path[!sapply(permute_path, function(x){all(x == 0)})]
              } else{
                  permute_path
              }


          })

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "data.frame"),
          function(
        expreMatrix, NB = NULL, testScore = NULL, gsTopology, weight, drop = TRUE){
              generate_permuted_scores(as.matrix(expreMatrix),NB, testScore,
                                       gsTopology, weight, drop)
          })

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "DGEList"),
          function(
        expreMatrix, NB = NULL, testScore = NULL, gsTopology, weight, drop = TRUE){
              expreMatrix <- cpm(expreMatrix$counts, log = TRUE)
              generate_permuted_scores(expreMatrix, NB, testScore,
                                       gsTopology, weight, drop)
          })

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "SummarizedExperiment"),
          function(
        expreMatrix, NB = NULL, testScore = NULL, gsTopology, weight, drop = TRUE){
              expreMatrix <- cpm(SummarizedExperiment::assay(expreMatrix), log = TRUE)
              generate_permuted_scores(expreMatrix, NB, testScore,
                                       gsTopology, weight, drop)
          })

#' @title Normalise test perturbation scores by permutation result and compute
#' permutation p-values
#'
#' @details Normalise the test perturbation scores generated by `weight_ss_fc()`
#' through the permuted perturbation scores derived from the
#' `generate_permuted_scores()` function. The mean absolute deviation(MAD) and
#' median of perturbation scores for each pathway are firstly derived from the
#' permuted perturbation scores. The test perturbation scores are then converted
#' to robust z-scores using MADs and medians calculated.
#'
#' Additionally, by assessing the proportion of permuted scores that are more
#' extreme than the test perturbation score within each pathway, the permuted
#' p-value of individual test perturbation scores will be computed.
#'
#' @param permutedScore A list. Output of `generate_permuted_scores`
#' @param testScore A `data.frame`. Output of `pathway_pert`
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
#' genePertScore <- raw_gene_pert(ls$weighted_logFC, gsTopology)
#'
#' # sum gene-wise perturbation scores to derive the pathway-level
#' # single-sample perturbation scores
#' pathwayPertScore <- pathway_pert(genePertScore, ls$weighted_logFC)
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
    tested_gs <- unique(testScore$gs_name)

    # Check if null distribution was simulated for all tested gene-sets
    if (!all(tested_gs %in% names(permutedScore))){
        warning("Permutation scores not available for some tested gene-sets.
                Those gene-sets will be removed.")
        #if not, gene-sets without permuted scores will be removed
        testScore <- dplyr::filter(testScore, gs_name %in% names(permutedScore))
    }

    pAdj_method <- match.arg(pAdj_method, p.adjust.methods)
    sortBy <- match.arg(sortBy)

    # Compute permuted p-values by counting the numbers of permuted scores that
    # are more extreme than the test scores

    pvalue_ls <- lapply(unique(testScore$gs_name), function(x){
        temp <- dplyr::filter(testScore, gs_name == x)
        mutate(
            temp,
            pvalue = vapply(temp$score, function(y){
                # count # of permuted values as or more extreme than test score
                sum(
                    # use >= so the smallest possible pvalue is 1/NB
                    abs(permutedScore[[x]]) >= abs(y)
                )  /
                    # divide by the total number of permuted values
                    length(permutedScore[[x]])
            }, numeric(1))
        )
    })

    pvalues <- bind_rows(pvalue_ls)

    # Compute MAD and MEDIAN from the permuted scores

    summary_func <- function(x){c(MAD = mad(x), MEDIAN = median(x))}
    summaryScore <- t(sapply(permutedScore, summary_func))
    gs_name <- rownames(summaryScore)
    summaryScore <- mutate(as.data.frame(summaryScore), gs_name = gs_name)

    # Remove gene-set with MAD equals 0 as denominator cannot be 0

    summaryScore <- dplyr::filter(summaryScore, summaryScore$MAD != 0)
    summaryScore <- left_join(summaryScore, testScore, by = "gs_name",
                              multiple = "all")

    # compute robust z-scores
    summaryScore <- mutate(
        summaryScore, robustZ =
            (summaryScore$score - summaryScore$MEDIAN)/summaryScore$MAD)
    # summaryScore <- mutate(
    #     summaryScore, z_pvalue = 2*pnorm(-abs(summaryScore$robustZ)))

    summaryScore <- left_join(summaryScore, pvalues,
                              by = c("gs_name", "score", "sample"))
    summaryScore <- split(summaryScore, f = summaryScore$sample)
    summaryScore <- lapply(summaryScore, mutate,
                           adjPvalue = p.adjust(pvalue, pAdj_method),
                           # adjPvalue_z = p.adjust(z_pvalue, pAdj_method)
                           )
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

    # generate all permutation pairs of column numbers
    Pairs <- gtools::permutations(
        n = length(colNumber), r = 2, v = colNumber)

    # if the user specified NB values is smaller than the maximum number of
    # possibilities, randomly sample NB permutations from all permutations
    if (NB < nrow(Pairs)){
        Pairs <- Pairs[sample(nrow(Pairs), NB, replace = FALSE),]}

    # Compute weighted FCs fro all permuted pairs
    apply(Pairs, 1, function(x){
        (expreMatrix[,x[1]] -
             expreMatrix[,x[2]]) * weight
    })
}

