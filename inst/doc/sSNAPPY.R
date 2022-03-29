## ---- echo=FALSE--------------------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE, crop = NULL)

## ----setup--------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly=TRUE))install.packages("BiocManager")
BiocManager::install("sSNAPPY")
library(sSNAPPY)

## ----otherPackages------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(ggplot2)
library(cowplot)
library(DT)

## ----data---------------------------------------------------------------------
# check if samples included in the logCPM matrix and metadata dataframe are identical
setequal(colnames(logCPM_example), metadata_example$sample)
# View sample metadata
metadata_example %>%
    datatable(
        filter = "top"
    )


## ----convertRownames----------------------------------------------------------
if (!requireNamespace("AnnotationHub", quietly=TRUE))BiocManager::install("AnnotationHub")
if (!requireNamespace("ensembldb", quietly=TRUE))BiocManager::install("ensembldb")
ah <- AnnotationHub::AnnotationHub()
ah <- AnnotationHub::subset(ah,genome == "GRCh38" & title == "Ensembl 101 EnsDb for Homo sapiens")
ensDb <- ah[[1]]
rownames(logCPM_example) <- ensembldb::mapIds(ensDb, rownames(logCPM_example), "ENTREZID", keytype = "GENEID")
# Remove genes that couldn't be matched to entrez IDs
logCPM_example <- logCPM_example[!is.na(rownames(logCPM_example)),]
head(logCPM_example)

## ----ssFC---------------------------------------------------------------------
#compute weighted single sample logFCs
weightedFC <- weight_ssFC(logCPM_example, metadata = metadata_example,
factor = "patient", control = "Vehicle")

## ----lowess, fig.width=10,fig.height=4, echo=FALSE----------------------------
perSample_FC <- sapply(levels(metadata_example$patient), function(x){
    temp <- logCPM_example[1:1000,str_detect(colnames(logCPM_example), x)] 
    ratio <- temp[, str_detect(colnames(temp), "Vehicle", negate = TRUE)] - temp[, str_detect(colnames(temp), "Vehicle")] 
}, simplify = FALSE) %>%
    do.call(cbind,.)
aveCPM <- logCPM_example[1:1000,] %>%
    rowMeans() %>%
    enframe(name = "gene_id", 
            value = "aveCPM")
p1 <- perSample_FC %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(cols = -"gene_id",
                 names_to = "name",
                 values_to = "logFC") %>%
    left_join(aveCPM) %>%
    ggplot(aes(aveCPM, logFC)) +
    geom_point() +
    labs(y = "sslogFC", 
         x = "Average logCPM") +
    theme(
        panel.background = element_blank()
    )
p2 <- data.frame(
    gene_id = rownames(perSample_FC),
    variance = perSample_FC %>%
        apply(1,var)) %>%
    left_join(aveCPM) %>%
    ggplot(aes(aveCPM, variance)) +
    geom_point() +
    geom_smooth(method = "loess") +
    labs(y = "Variance in ssLogFCs", 
         x = "Average logCPM") +
    theme(
        panel.background = element_blank()
    )
plot_grid(p1, 
          p2)


## ----pathwayDatabases---------------------------------------------------------
if (!requireNamespace("graphite", quietly=TRUE))
    install.packages("graphite")
graphite::pathwayDatabases() %>%
  dplyr::filter(species ==  "hsapiens") %>%
  pander::pander()

## ----weightedAdjMatrix--------------------------------------------------------
weightedAdjMatrix(database = "kegg", outputDir = "gsTopology.rda")
load("gsTopology.rda")
head(names(gsTopology))

## -----------------------------------------------------------------------------
# weightedAdjMatrix( database = "kegg", pathwayName = c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)","Pentose phosphate pathway"), outputDir = "gsTopology.rda")
# load("gsTopology.rda")
# names(gsTopology)

## ----ssPertScore--------------------------------------------------------------
ssPertScore <- perturbationScore(weightedFC$logFC, gsTopology)
head(ssPertScore)

## ----permutedScore------------------------------------------------------------
permutedScore <- generate_PermutedScore(logCPM_example, numOfTreat = 3, NB = 1000, gsTopology = gsTopology, weight = weightedFC$weight)

## ----NullDistribution, fig.width=10, fig.height=6-----------------------------
pl <- permutedScore %>%
    keep(~all(.!=0)) %>%
    .[sample(seq_along(.), 6)] %>%
    lapply(function(x){
        ggplot(mapping = aes(x)) + 
            geom_histogram() +
            xlab("Perturbation Score")
    }) 
suppressMessages(plot_grid(plotlist = pl, 
                           nrow = 2))

## ----permutedScore------------------------------------------------------------
permutedScore <- generate_PermutedScore(logCPM_example, numOfTreat = 3, NB = 1000, gsTopology = gsTopology, weight = weightedFC$weight)

## ----NullDistribution, fig.width=10, fig.height=6-----------------------------
pl <- permutedScore %>%
    keep(~all(.!=0)) %>%
    .[sample(seq_along(.), 6)] %>%
    lapply(function(x){
        ggplot(mapping = aes(x)) + 
            geom_histogram() +
            xlab("Perturbation Score")
    }) 
suppressMessages(plot_grid(plotlist = pl, 
                           nrow = 2))

## ----normalisedScores---------------------------------------------------------
normalisedScores <- normaliseByPermutation(permutedScore, ssPertScore)
normalisedScores %>%
    dplyr::filter(adjPvalue < 0.05) %>%
    left_join(metadata_example) %>%
    mutate_at(vars(c("sample", "gs_name")), as.factor) %>%
    mutate_if(is.numeric, sprintf, fmt = '%#.4f') %>%
    mutate(Direction = ifelse(robustZ < 0, "Inhibited", "Activation")) %>%
    dplyr::select(
        sample,patient, Treatment = treatment, `Perturbation Score` = robustZ, Direction,
        `Gene-set name` = gs_name, 
        `P-value` = pvalue, 
        FDR = adjPvalue
    ) %>%
    datatable(
        filter = "top", 
        options = list(
            columnDefs = list(list(targets = "Direction", visible = FALSE))
        )) %>% 
    formatStyle(
        'Perturbation Score', 'Direction',
        backgroundColor = styleEqual(c("Inhibited", "Activation"), c('lightblue', 'indianred'))
    )

## ----sigGS_nt_zscore, fig.width= 15, fig.height=5-----------------------------
pl <- normalisedScores %>%
    dplyr::filter(adjPvalue < 0.05) %>%
    split(f = .$sample) %>%
    lapply(
        plot_gsNetwork, 
        layout = "dh",
        gsTopology = gsTopology, 
        colorBy = "robustZ"
        
    )
plot_grid(
    plotlist = pl, 
    nrow = 1
)

## ----sigGS_nt_pvalue, fig.width= 15, fig.height=5-----------------------------
pl <- normalisedScores %>%
    dplyr::filter(adjPvalue < 0.05) %>%
    split(f = .$sample) %>%
    lapply(
        plot_gsNetwork, 
        layout = "dh",
        gsTopology = gsTopology, 
        colorBy = "pvalue", 
        color_lg_title = "P-value"
    )
plot_grid(
    plotlist = pl, 
    nrow = 1
)

## ----fit----------------------------------------------------------------------
fit <- normalisedScores %>%
    left_join(metadata_example) %>%
    split(f = .$gs_name) %>%
    #.["Estrogen signaling pathway"] %>%
    lapply(function(x)lm(robustZ ~ 0 + treatment + PR, data = x)) %>%
    lapply(summary)
treat_sig <- sapply(names(fit), function(x){
    fit[[x]]$coefficients %>%
        as.data.frame() %>%
        .[1:2,] %>%
        dplyr::select(Estimate, 
                      pvalue = `Pr(>|t|)` ) %>%
        rownames_to_column("Treatment") %>%
        mutate(gs_name = x, 
               FDR = p.adjust(pvalue, "fdr"), 
               Treatment = str_remove_all(Treatment, "treatment")) 
}, simplify = FALSE) %>%
    bind_rows() 

## ----treat_sig_DT-------------------------------------------------------------
treat_sig %>% 
    dplyr::filter(FDR < 0.05) %>%
    mutate_at(vars(c("Treatment", "gs_name")), as.factor) %>%
    mutate_if(is.numeric, sprintf, fmt = '%#.4f') %>%
    mutate(Direction = ifelse(Estimate < 0, "Inhibited", "Activation")) %>%
    dplyr::select(
        Treatment, `Perturbation Score` = Estimate, Direction,
        `Gene-set name` = gs_name, 
        `P-value` = pvalue, 
        FDR
    ) %>%
    datatable(
        filter = "top", 
        options = list(
            columnDefs = list(list(targets = "Direction", visible = FALSE))
        )) %>% 
    formatStyle(
        'Perturbation Score', 'Direction',
        backgroundColor = styleEqual(c("Inhibited", "Activation"), c('lightblue', 'indianred'))
    )

## ----fig.width=6, fig.height=4------------------------------------------------
treat_sig %>% 
    dplyr::filter(FDR < 0.05, Treatment == "R5020") %>%
    dplyr::rename(robustZ = Estimate) %>%
    plot_gsNetwork(
        layout = "stress",
        gsTopology = gsTopology, 
        colorBy = "robustZ"
    ) +
    theme(
        panel.grid = element_blank(), 
        panel.background = element_blank()
    ) 

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

