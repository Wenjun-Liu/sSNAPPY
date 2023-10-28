#' @title Compute Single-sample Pathway-level Perturbation Score
#'
#' @description Substract ssFC from the raw gene-level perturbation scores within
#' each sample and sum gene-wise raw perturbation scores to derive single-sample
#' perturbation scores for each pathway.
#' @param genePertScore List of gene-wise raw perturbation score matrices
#' generated using function `raw_gene_pert()`
#' @param weightedFC A matrix of weighted ssFC generated using function `weight_ss_fc`
#' @param drop logic(1). Whether to drop pathways with all zero scores
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr set_colnames
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
#' # sum gene-wise perturbation scores to derive the pathway-level single-sample
#' # perturbation scores
#' pathwayPertScore <- pathway_pert(genePertScore, ls$weighted_logFC)
#' @export
pathway_pert <- function(genePertScore, weightedFC, drop = TRUE){
    # extract all unique pathway genes and find ones that are not expressed
    notExpressed <- setdiff(
        unique(
            unlist(unname(lapply(genePertScore, rownames)))),
        rownames(weightedFC))
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
        # extract FC for pathway genes
        sub_FC <- weightedFC[rownames(weightedFC) %in% rownames(raw_per),]
        # match FC matrix with gene pert score gene and sample order
        sub_FC <- sub_FC[
            match(rownames(raw_per), rownames(sub_FC)),
            match(colnames(raw_per), colnames(sub_FC))]
        # extract FC from raw gene pert to get net gene pert
        net_per <- raw_per - sub_FC
        # sum within each column to get pathway-level pert
        temp <- as.data.frame(apply(net_per, 2, sum))
        temp <- set_colnames(temp, "score")
        if(all(temp$score == 0) & drop){
            return(NULL)
        } else{
            temp <- rownames_to_column(temp, "sample")
            mutate(temp, gs_name = x)
        }
    })

    bind_rows(PF)
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

