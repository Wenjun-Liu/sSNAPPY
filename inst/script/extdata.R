#' gsTopology: retrieved KEGG pathway topologies
#'
#' Pathway topologies retrieved from KEGG pathway as weightede adjaency matrix using the `retrieve_topology()` function
#' @format A list
"gsTopology"

#' permutedScore: perturbation scores derived from 6 rounds of permutation
#'
#' Initially, 1000 rounds of permutation were performed using the `generate_permuted_scores()` function
#' and six of them were randomly sampled and saved as an example.
#' @format A list
"permutedScore"

#' normalisedScores: test perturbation scores noramlised through permutation
#'
#' Test perturbation scores normalized through 1000 rounds of permutation and associated median, robust z-score and pvalue (derived using `normalisedScores()`)
#' @format A data frame
"normalisedScores"
