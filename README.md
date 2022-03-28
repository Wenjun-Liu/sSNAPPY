# sSNAPPY: Singel Sample directioNAl Pathway Perturbation analYsis


`sSNAPPY` is a package for computing directional pathway perturbation scores from logCPM in RNA-seq data. The final output can be used to test the significance of pathway perturbation at both individual-sample and overall treatment levels.
    
The overall workflow of the analysis is:
- Compute weighted single-sample logFCs for each of the treated samples
- Retrieve pathway topologies from chosen database, such as *KEGG* or *wikiPathways*
- Calculate a raw directional pathway perturbation score per pathway per treated sample
- Simulate the null distribution of perturbation scores through sample permutation
- Use the `median` and `MAD` of the null distribution to normalise raw scores and derive significance of pathway perturbation for individual samples
- Fit models to test overall treatment effects on pathway perturbation

# Installation Instructions

To install this package, please use BiocManager.

`install.packages("BiocManager")
BiocManager::install("steveped/extraChIPs")`
