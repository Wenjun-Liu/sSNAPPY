#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;

//[[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List GenePertScore_RCPP(const List& BminsI, arma::mat weightedFC, const CharacterVector& expressedG) {

    // s number of samples
    int s = weightedFC.n_cols;

    // // ng number of genes
    // int ng = weightedFC.n_rows;

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

        arma::mat gene_pert(pathwayG.size(),s);

        // loop through each sample
        for (int j=0; j<s; j++){
            arma::vec subset = subsetFC.col(j);

            arma::vec pf = solve(X,-subset,arma::solve_opts::fast);
            arma::vec diff = pf - subset;
            gene_pert.col(j) = diff;
        }
        output[i] = gene_pert;
    }

    output.names() = BminsI.names();
    return(output);
}


