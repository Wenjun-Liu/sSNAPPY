#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

List permutedFC_RCPP(const NumericMatrix& logCPM, int NB){

    List output(NB);
    int s = logCPM.cols();
    IntegerVector index = seq(1, s);
    IntegerVector randomI = sample(index, s);
    Eigen::MatrixXd X(as<Eigen::MatrixXd>(logCPM));

    // NumericMatrix randomM = logCPM(_,randomI);
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(s);
    perm.setIdentity();
    std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
    Eigen::MatrixXd test = X * perm;
    return(test);
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
