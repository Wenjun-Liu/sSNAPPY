#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List permutedPertScore_RCPP(const arma::mat& X, const CharacterVector& pathwayG, const CharacterVector& expressedG, const List& permutedFC, int newS) {

    // number of permutation
    int NB = permutedFC.length();

    IntegerVector offset(pathwayG.length(),1);
    IntegerVector index = match(pathwayG, expressedG) - offset;
    arma::uvec indexA = as<arma::uvec>(index);

    List output(NB);

    for (int i=0; i<NB; i++){

        if (i % 5 == 0)
            Rcpp::checkUserInterrupt();

        arma::mat thisPermu = permutedFC[i];
        arma::mat subsetFC = thisPermu.rows(indexA);
        arma::vec tA(newS);

        // loop through each sample
        for (int j=0; j<newS; j++){
            arma::vec subset = subsetFC.col(j);

            arma::vec pf = arma::solve(X,-subset,arma::solve_opts::fast);
            arma::vec diff = pf - subset;
            tA(j) = sum(diff);

        }

        // arma::uvec vI = arma::linspace<arma::uvec>(m*newS, ((m+1)*newS)-1, newS);
        //
        // allPerResult.elem(vI) = tA;

        output[i] = tA;

    }

    return(output);
}


