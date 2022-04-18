#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;

//[[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List ssPertScore_RCPP(const List& BminsI, arma::mat weightedFC, const CharacterVector& expressedG, const CharacterVector& sample) {

    // s number of samples
    int s = weightedFC.n_cols;

    // number of pathways
    int p = BminsI.length();

    //create empty vector as output
    List output(p);

    for (int i=0; i<p; i++){

        if (i % 100 == 0)
            Rcpp::checkUserInterrupt();

        CharacterVector pathwayG = rownames(BminsI[i]);
        IntegerVector offset( pathwayG.length(),1);
        IntegerVector index = match(pathwayG, expressedG) - offset;
        arma::uvec indexA = as<arma::uvec>(index);
        arma::mat subsetFC = weightedFC.rows(indexA);
        arma::mat X(as<arma::mat>(BminsI[i]));

        NumericVector tA(s);

        // loop through each sample
        for (int j=0; j<s; j++){
            arma::vec subset = subsetFC.col(j);

            arma::vec pf = solve(X,-subset,arma::solve_opts::fast);
            arma::vec diff = pf - subset;
            tA[j] = sum(diff);

        }
        tA.names() = sample;
        output[i] = tA;
    }

    output.names() = BminsI.names();
    return(output);
}


