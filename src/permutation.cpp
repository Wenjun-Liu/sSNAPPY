#include <Rcpp.h>
using namespace Rcpp;

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

List permutedFC_RCPP(int NB, NumericMatrix logCPM, List sampleInpairs, NumericVector weight){

    //  create all sample label permutations
    CharacterVector s = colnames(logCPM);
    std::sort(s.begin(), s.end());

    List Permu(NB);

    do {
        Permu.push_back(group_sorted);
    } while (std::next_permutation(group_sorted.begin(), group_sorted.end())) ;

    int nTotalPer = Permu.size();
    }

// sapply(1:NB, function(x){
//
// # permute sample labels to get permuted logCPM
//     permutedCPM <- logCPM
//         colnames(permutedCPM) <- sample(colnames(permutedCPM), ncol(permutedCPM))
//
// # Built permuted logFCs based on the permuted logCPM
//         permutedFC <- sapply(names(sampleInpairs), function(y){
//             permutedCPM[, sampleInpairs[[y]]$treatedSample] - permutedCPM[, sampleInpairs[[y]]$contrSample]
//         }, simplify = FALSE)  %>%
//             do.call(cbind,.)
//
// # Multiply permuted FCs by gene-wise weights
//                 permutedFC * weight
//
// }, simplify = FALSE)
