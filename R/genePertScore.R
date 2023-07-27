#' @title Compute Single-sample Pathway-level Perturbation Score
#'
#' @description Sum gene-wise raw perturbation scores within each sample to derive single-sample perturbation scores for each pathway
#' @param genePertScore List of gene-wise raw perturbation score matrices generated using function `raw_gene_pert()`
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr set_colnames
#'
#' @references Tarca AL, Draghici S, Khatri P, Hassan SS, Mittal P, Kim JS, Kim CJ, Kusanovic JP, Romero R. A novel signaling pathway impact analysis.
#' Bioinformatics. 2009 Jan 1;25(1):75-82.
#' @return A data.frame with 3 columns: score (single-sample pathway-level perturbation score), sample, and gs_name (gene-set name)
#' @examples
#' #compute weighted single sample logFCs
#' data(metadata_example)
#' data(logCPM_example)
#' metadata_example <- dplyr::mutate(metadata_example, treatment = factor(
#'    treatment, levels = c("Vehicle", "E2+R5020", "R5020")))
#' ls <- weight_ss_fc(logCPM_example, metadata = metadata_example,
#' groupBy = "patient", treatColumn = "treatment", sampleColumn = "sample")
#' # extract all the KEGG pathways
#' gsTopology <- retrieve_topology(database = "kegg", species = "hsapiens")
#' # compute raw gene-wise perturbation scores
#' genePertScore <- raw_gene_pert(ls$weighted_logFC, gsTopology)
#' # sum gene-wise perturbation scores to derive the pathway-level single-sample perturbation scores
#' pathwayPertScore <- pathway_pert(genePertScore, ls$weighted_logFC)
#' @export
pathway_pert <- function(genePertScore, weightedFC){
    browser()
    # extract all unique pathway genes and find ones that are not expressed
    notExpressed <- setdiff(unique(unlist(unname(lapply(genePertScore, rownames)))), rownames(weightedFC))
    if (length(notExpressed) != 0){
        # set the FCs of unexpressed pathway genes to 0
        temp <- matrix(0, nrow = length(notExpressed), ncol = ncol(weightedFC))
        rownames(temp) <- notExpressed
        colnames(temp) <- colnames(weightedFC)
        # set the weights of unexpressed pathway genes to 0
        weightedFC <- rbind(weightedFC, temp)
        }

    # sum pathway perturbation scores for each pathway
    PF <- lapply(names(genePertScore), function(x){
        raw_per <- genePertScore[[x]][rownames(genePertScore[[x]]) %in% rownames(weightedFC), ]
        sub_FC <- weightedFC[rownames(weightedFC) %in% rownames(raw_per),]
        sub_FC <- sub_FC[
            match(rownames(raw_per), rownames(sub_FC)),
            match(colnames(raw_per), colnames(sub_FC))]
        net_per <- raw_per - sub_FC
        temp <- as.data.frame(apply(net_per, 2, sum))
        temp <- set_colnames(temp, "score")
        temp <- rownames_to_column(temp, "sample")
        temp <- mutate(temp, gs_name = x)
    })

    bind_rows(PF)
}

#' @title Rank genes by perturbation scores within each sample
#'
#' @description Rank genes by gene-wise raw perturbation scores within each sample to compare
#' genes' contributions to pathway perturbations.
#' @details Ranking is performed within each sample each pathway. If in a given pathway, both positive and negative gene-wise perturbation scores exist, positive
#' and negative scores are ranked separately, where the larger a positive rank, the more the gene contributed to the pathway's activation, and the smaller a negative
#' rank, the more the gene contributed to the pathways' inhibition. When there's a tie in two gene's perturbation score within a sample, the mean of the indices is used.
#' @param gsTopology List of pathway topology matrices generated using function `retrieve_topology()`
#' @param genePertScore List of gene-wise raw perturbation score matrices generated using function `raw_gene_pert()`
#' @importFrom tibble rownames_to_column enframe
#' @importFrom magrittr set_colnames
#' @return A list where each element is a matrix corresponding to a pathway. Each column of an element corresponds to a sample, and each row corresponds to a pathway gene.
#' @examples
#' #compute weighted single sample logFCs
#' data(metadata_example)
#' data(logCPM_example)
#' metadata_example <- dplyr::mutate(metadata_example, treatment = factor(
#'    treatment, levels = c("Vehicle", "E2+R5020", "R5020")))
#' ls <- weight_ss_fc(logCPM_example, metadata = metadata_example,
#' groupBy = "patient", treatColumn = "treatment", sampleColumn = "sample")
#' # extract all the KEGG pathways
#' gsTopology <- retrieve_topology(database = "kegg", species = "hsapiens")
#' # compute raw gene-wise perturbation scores
#' genePertScore <- raw_gene_pert(ls$weighted_logFC, gsTopology)
#' # rank genes by gene-wise perturbation scores within each sample
#' # to compare their contributions to pathway perturbation
#' geneRank <- rank_gene_pert(genePertScore, gsTopology)
#' @export
rank_gene_pert <- function(genePertScore, gsTopology){

    # check if all topology infor is available for pathways in the PertScore list
    if(!all(names(genePertScore) %in% names(gsTopology)))
        stop("Pathway topology information missing for some pathways.")

    # # extract sample names from the FC matrix
    # sampleName <- colnames(genePertScore[[1]])

    output <- lapply(names(genePertScore), function(x){
        temp <- genePertScore[[x]]
        # remove genes whose perturbation scores are 0 acorss all samples
        temp <- temp[apply(temp, 1,function(y){any(y != 0)}), , drop = FALSE]

        if (nrow(temp) > 1){

            if(any(temp > 0) & any(temp < 0)){
                pos_rank <- apply(temp, 2, function(x){rank(x[x > 0])})
                neg_rank  <- apply(temp, 2, function(x){-rank(abs(x[x < 0]))})
                temp <- sapply(names(pos_rank), function(y){
                    c(pos_rank[[y]], neg_rank[[y]])}, simplify = FALSE)
                temp <- lapply(names(temp), function(y){tibble::enframe(temp[[y]], value = y, name = "gene_id")})
                suppressMessages(Reduce(left_join, temp))

            } else {
                temp <-  apply(temp, 2, function(x){
                    sign(x)*rank(abs(x))
                })
                temp
            }

        } else {
            temp <- sign(temp)
            rownames_to_column(as.data.frame(temp), "gene_id")
        }

    })

    names(output) <- names(genePertScore)
    output

}


#' @title Compute Gene-wise Perturbation Score
#'
#' @description Propagate weighted single sample logFCs down the pathway topologies
#' to compute gene-wise perturbation score per gene per sample per pathway
#'
#' @details This function use the algorithm adopted from `SPIA` (see citation) to
#' integrate genes' changes in expression and gene-gene interaction to compute
#' gene-wise perturbation score per gene per sample per pathway. The rownames of
#' the weighted single sample logFC matrix and the pathway topology matrices must
#' use the same type of gene identifier (ie. entrez ID).
#'
#' Pathways with zero perturbation scores across all genes and samples will be
#' dropped from the output.
#'
#' @param weightedFC A matrix of weighted single sample logFCs
#' derived from function `weight_ss_fc()`
#' @param gsTopology List of pathway topology matrices generated using function
#' `retrieve_topology()`
#' @references Tarca AL, Draghici S, Khatri P, Hassan SS, Mittal P, Kim JS,
#' Kim CJ, Kusanovic JP, Romero R. A novel signaling pathway impact analysis.
#' Bioinformatics. 2009 Jan 1;25(1):75-82.
#' @return A list where each element is a matrix corresponding to a pathway.
#' Each column of an element corresponds to a sample, and each row corresponds
#' to a pathway gene.
#' @examples
#' #compute weighted single sample logFCs
#' data(metadata_example)
#' data(logCPM_example)
#' metadata_example <- dplyr::mutate(metadata_example, treatment = factor(
#'    treatment, levels = c("Vehicle", "E2+R5020", "R5020")))
#' ls <- weight_ss_fc(logCPM_example, metadata = metadata_example,
#' groupBy = "patient", treatColumn = "treatment", sampleColumn = "sample")
#' # extract all the KEGG pathways
#' gsTopology <- retrieve_topology(database = "kegg", species = "hsapiens")
#' # compute raw gene-wise perturbation scores
#' genePertScore <- raw_gene_pert(ls$weighted_logFC, gsTopology)
#' @export
raw_gene_pert <- function(weightedFC, gsTopology){

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

    sampleName <- colnames(weightedFC)
    allGene <- rownames(weightedFC)
    GP <- lapply(gsTopology, function(x){
        gs_sub <- x[rownames(x) %in% allGene, colnames(x) %in% allGene ]
        if (abs(det(gs_sub))>1e-7){
            geneP_ls <- lapply(seq_len(ncol(weightedFC)), function(y){
                de <- weightedFC[rownames(weightedFC) %in% rownames(gs_sub), y]
                de <- de[match(rownames(gs_sub), names(de))]
                solve(t(gs_sub), -de)
            })
            geneP_df <- do.call(cbind,geneP_ls)
            colnames(geneP_df) <- sampleName
            geneP_df
        } else {
            NULL
        }
    })

    # Remove pathways that were unsolvable
    GP <- GP[!sapply(GP, is.null)]
    # Remove list elements that are null or all zeros
    GP[sapply(GP, function(x){any(x != 0)})]

}
