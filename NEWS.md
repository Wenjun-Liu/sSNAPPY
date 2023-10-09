# Changes in version 1.4.4 (2023-10-09)
- Remove databases that are no longer supported: pathbank, panther, pharmgkb, smpdb

# Changes in version 1.4.3 (2023-07-31)
- Updated the perturbation scoring step to account for the orientiaton of topology 
matrices for KEGG pathways (row - Downstream genes; column - Upstream genes)
- Add `prefix` parameter to the `weight_ss_fc` function to allow user-specified 
prefix that are other than "ENTREZID:"

# Changes in version 1.4.2 (2023-07-12)
- Updated the permtuation strategy so that the `generate_permuted_scores` function now construct all possible permuted pairs by default
- Updated the `plot_community` function so KEGG pathways that are not assigned to
categories will be ignored in community labeling

# Changes in version 1.0.3 (2022-11-24)
- Replaced `compute_perturbation_score` with `raw_gene_pert`, `path_gene_per` 
and `rank_gene_pert` to allow the scoring of gene-wise perturbation
- Added 3 new visualisation functions: `plot_gene_contribution`, `plot_community` and `plot_gs2gene`
- Demonstrated the use of new functions in the updated version of vignette

# Changes in version 1.0.2 (2022-08-03)
- Fixed a missing column name problems in `weight_ss_fc`

# Changes in version 1.1.0 (2022-06-27)
- Fixed a few typos

# Changes in version 0.99.0 (2022-03-28)
- Submitted to Bioconductor

