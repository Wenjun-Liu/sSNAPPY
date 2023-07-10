#' @title Permute sample labels to simulate null distribution of perturbation scores
#'
#' @description Simulate null distributions of perturbation scores for each pathway through sample permutation.
#'
#' @details
#' This \code{generate_permuted_scores} function is a generic function that can deal with
#' multiple types of inputs. It firstly randomly permute sample labels `NB` times to generate
#' permuted logFCs, which are then used to compute permuted perturbation scores for each pathway.
#'
#' The function outputs a list that is of the same length as the list storing pathway topology
#' matrices. Each element of the output list is for a pathway and contains a vector of permuted
#' perturbation score of length. The permuted perturbation scores will be used to estimate the
#' null distributions of perturbation scores.
#'
#'If the input is S4 object of \code{DGEList or SummarizedExperiment}, gene expression matrix will
#'be extracted and converted to a logCPM matrix.
#'
#' The default number of permutation (`NB`) is set to 1000. If the requested `NB` is larger than the
#' maximum number of permutations possible, `NB` will be set to the largest number of permutations
#' possible instead.
#'
#' @param expreMatrix `matrix` and `data.frame` of logCPM, or `DGEList`/`SummarizedExperiment`
#' storing gene expression counts. Feature names need to be in entrez IDs
#' @param numOfTreat Number of treatments (including control)
#' @param NB Number of permutations to perform
#' @param testScore Optional. Users can provide the test perturbation score
#' `data.frame` (ie. output of `pathwayPertScore`) to restrict the permutatio
#' step only to pathways with non-zero test scores in at least one sample.
#' @param gsTopology List of pathway topology matrices generated using function `retrieve_topology`
#' @param weight A vector of gene-wise weights derived from function `weight_ss_fc`
#' @param BPPARAM The parallel back-end to uses, if not specified, it is defaulted to the one returned by \code{BiocParallel::bpparam()}.
#' @import Rcpp
#' @importFrom BiocParallel bpparam bplapply
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
#' permutedScore <- generate_permuted_scores(logCPM_example, numOfTreat = 3,
#'  NB = 10, gsTopology = gsTopology, weight = ls$weight)
#'
#' # To see what other parallel back-end can be used:
#'  BiocParallel::registered()
#'  }
#' @export
setGeneric("generate_permuted_scores",
           function(expreMatrix, numOfTreat,NB = 1000, testScore = NULL,
                    gsTopology, weight, BPPARAM = BiocParallel:: bpparam())
               standardGeneric("generate_permuted_scores"))

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "matrix"),
          function(expreMatrix, numOfTreat, NB = 1000, testScore = NULL,
                   gsTopology, weight, BPPARAM = BiocParallel:: bpparam()){
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

              # test if the required number of permutations is bigger than the maximum number of permutations possible
              NB <- min(NB, factorial(ncol(expreMatrix)))

              # if test scores are provided, filter out pathways whose test scores
              # were all zeros across all samples


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

              permutedFC <- .generate_permutedFC(expreMatrix, numOfTreat, NB, weight)

              allG <- rownames(expreMatrix)
              newS <-  ncol(permutedFC[[1]])


              BiocParallel::bplapply(gsTopology, function(x){
                  temp <- permutedPertScore_RCPP(X = x, pathwayG = rownames(x), allG,
                                                 permutedFC = permutedFC, newS)
                  unlist(temp)
              }, BPPARAM = BPPARAM)
          })

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "data.frame"),
          function(expreMatrix, numOfTreat,NB = 1000, testScore = NULL,
                   gsTopology, weight, BPPARAM = BiocParallel:: bpparam()){
              generate_permuted_scores(as.matrix(expreMatrix), numOfTreat,NB,
                                      testScore, gsTopology, weight,
                                      BPPARAM = BiocParallel:: bpparam())
          })

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "DGEList"),
          function(expreMatrix, numOfTreat,NB = 1000, testScore = NULL,
                   gsTopology, weight, BPPARAM = BiocParallel:: bpparam()){
              expreMatrix <- cpm(expreMatrix$counts, log = TRUE)
              generate_permuted_scores(expreMatrix, numOfTreat,NB,
                                       testScore, gsTopology, weight,
                                       BPPARAM = BiocParallel:: bpparam())
          })

#' @rdname generate_permuted_scores
setMethod("generate_permuted_scores",
          signature = signature(expreMatrix = "SummarizedExperiment"),
          function(expreMatrix, numOfTreat,NB = 1000, testScore = NULL,
                   gsTopology, weight, BPPARAM = BiocParallel:: bpparam()){
              expreMatrix <- cpm(SummarizedExperiment::assay(expreMatrix), log = TRUE)
              generate_permuted_scores(expreMatrix, numOfTreat,NB,
                                       testScore, gsTopology, weight,
                                       BPPARAM = BiocParallel:: bpparam())
          })

#' @title Normalise test perturbation scores by permutation results
#'
#' @details Normalise the test perturbation scores generated by `weight_ss_fc()` through the permuted perturbation scores derived from
#' the `generate_permuted_scores()` function. The mean absolute deviation(MAD) and median of perturbation scores for each pathway are firstly
#' derived from the permuted perturbation scores. The test perturbation scores are then converted to robust z-scores using MADs and medians
#' calculated.
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
#' permutedScore <- generate_permuted_scores(logCPM_example, numOfTreat = 3,
#'  NB = 5, gsTopology = gsTopology, weight = ls$weight)
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


.generate_permutedFC <- function(expreMatrix, numOfTreat,
                                 NB, weight){
    nSample <- ncol(expreMatrix)

    # calculate number of possible permutation pairs
    nPairs <- nSample %/% 2

    lapply(seq_len(NB), function(x){
        # generate random index used to put samples into pairs
        index <- sample(1:nSample, size = 2*nPairs, replace = FALSE)
        # indexes to index the random index
        index_pair <- seq(1, length(index), by = 2)
        temp <- lapply(index_pair, function(y){
            (expreMatrix[,index[y+1]] -
                 expreMatrix[,index[[y]]]) * weight
        })
        do.call(cbind, temp)
    })

}



