#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List permutedFC_RCPP(arma::mat logCPM, int NB, int sEachp){

    List output(NB);
    int s = logCPM.n_cols;
    int p = s/sEachp;
    arma::uvec indices = arma::regspace<arma::uvec>(0, sEachp, s-1);

    for (int i=0; i<NB; i++){
        arma::uvec permutationI = randperm(s);
        arma::mat permutedCPM = logCPM.cols(permutationI);
        // arma::mat permutedFC = permutedCPM;
        arma::mat permutedFC(logCPM.n_rows, s-p);
        for (int j=0; j< p; j++){
            int temp = indices[j];
            arma::uvec x1 = arma::linspace<arma::uvec>(j*(sEachp-1),j*(sEachp-1)+sEachp-2, sEachp-1);
            arma::uvec x2 = arma::linspace<arma::uvec>(temp+1, temp+sEachp-1, sEachp-1);
            arma::uvec x3 =  arma::linspace<arma::uvec>(temp, temp, sEachp-1);
            permutedFC.cols(x1) = permutedCPM.cols(x2) - permutedCPM.cols(x3);
            }
         output[i] = permutedFC;


    }

    return(output);
}

// .generate_permutedFC_alt <- function(logCPM, metadata, factor, control, weight, NB, seed){
//
// metadata <- as.data.frame(metadata)
//     sEachp <- length(levels(as.factor(metadata$treatment)))
//
//     nSample <- nrow(metadata)
//     NB <- min(factorial(nSample), NB)
//     set.seed(seed)
//
//     sapply(1:NB, function(x){
//         colIndex <- sample(seq_len(ncol(logCPM)), ncol(logCPM))
//         logCPM <- logCPM[,colIndex]
//         permutedFC <- lapply(seq(1, ncol(logCPM), by = sEachp), function(x){
//             logCPM[,c((x+1):(x+(sEachp-1)))] -logCPM[,x]
//
//         })
//         permutedFC <- do.call(cbind,permutedFC)
//
// # Multiply permuted FCs by gene-wise weights
//             permutedFC * weight
//     }, simplify = FALSE)
//
// }
