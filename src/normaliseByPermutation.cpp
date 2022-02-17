#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List permutedFC_RCPP(arma::mat logCPM, int NB, int sEachp, arma::vec weight){

    List output(NB);
    int s = logCPM.n_cols;
    int p = s/sEachp;
    arma::uvec indices = arma::regspace<arma::uvec>(0, sEachp, s-1);

    for (int i=0; i<NB; i++){
        // for each permutation, randomly permutate column numbers
        arma::uvec permutationI = randperm(s);
        // use permutated column numbers as indices to get the permuted logCPM matrix
        arma::mat permutedCPM = logCPM.cols(permutationI);
        // create empty permutedFC matrix with same number of rows as logCPM and number of treated samples columns
        arma::mat permutedFC(logCPM.n_rows, s-p);

        // looping through each patients
        for (int j=0; j< p; j++){
            //get the indices of control sample for this atient
            int temp = indices[j];
            // indices for permuted FC
            arma::uvec x1 = arma::linspace<arma::uvec>(j*(sEachp-1),j*(sEachp-1)+sEachp-2, sEachp-1);
            // indices used to subset permuted logCPM (treated samples)
            arma::uvec x2 = arma::linspace<arma::uvec>(temp+1, temp+sEachp-1, sEachp-1);
            // indices used to subset permuted logCPM (control samples). There's just one control but
            // for the next step the same column will be extracted multiple times
            arma::uvec x3 =  arma::linspace<arma::uvec>(temp, temp, sEachp-1);
            // logFC(treated) -logFC(control) to compute ssFCs for each treated samples
            permutedFC.cols(x1) = permutedCPM.cols(x2) - permutedCPM.cols(x3);
            }

        for (int i=0; i<s-p; i++){
            permutedFC.col(i) = permutedFC.col(i) % weight;
            }

         output[i] = permutedFC;


    }

    return(output);
}


// [[Rcpp::export]]
List permutedPertScore_RCPP(const List& BminsI, const CharacterVector& expressedG, arma::mat LogCPM, int NB, int sEachp, arma::vec weight) {

    // number of pathways
    int pathway = BminsI.length();
    // number of samples
    int s = LogCPM.n_cols;
    // number of matching pairs
    int p = s/sEachp;

    // generate empty list the same length as pahtway topologies as output
    List output(pathway);

    // indices for sslogFC computation. seq(0, number of sample -1, by = samples each patient)
    arma::uvec indices = arma::regspace<arma::uvec>(0, sEachp, s-1);

    List allPermutedFC(NB);

    for (int i=0; i<NB; i++){

        Rcpp::checkUserInterrupt();

        // for each permutation, randomly permutate column numbers
        arma::uvec permutationI = randperm(s);
        // use permutated column numbers as indices to get the permuted logCPM matrix
        arma::mat permutedCPM = LogCPM.cols(permutationI);
        // create empty permutedFC matrix with same number of rows as logCPM and number of treated samples columns
        arma::mat permutedFC(LogCPM.n_rows, s-p);

        // looping through each patients
        for (int j=0; j< p; j++){

            //get the indices of control sample for this atient
            int temp = indices[j];
            // indices for permuted FC
            arma::uvec x1 = arma::linspace<arma::uvec>(j*(sEachp-1),j*(sEachp-1)+sEachp-2, sEachp-1);
            // indices used to subset permuted logCPM (treated samples)
            arma::uvec x2 = arma::linspace<arma::uvec>(temp+1, temp+sEachp-1, sEachp-1);
            // indices used to subset permuted logCPM (control samples). There's just one control but
            // for the next step the same column will be extracted multiple times
            arma::uvec x3 =  arma::linspace<arma::uvec>(temp, temp, sEachp-1);
            // logFC(treated) -logFC(control) to compute ssFCs for each treated samples
            permutedFC.cols(x1) = permutedCPM.cols(x2) - permutedCPM.cols(x3);
        }

        for (int n=0; n<s-p; n++){
            permutedFC.col(n) = permutedFC.col(n) % weight;
        }

        allPermutedFC[i] = permutedFC;

        }

    int newS = s-p;

    int totalPerS = NB*newS;

    for (int i=0; i<pathway; i++){

        Rcpp::checkUserInterrupt();

        CharacterVector pathwayG = rownames(BminsI[i]);
        IntegerVector offset( pathwayG.length(),1);
        IntegerVector index = match(pathwayG, expressedG) - offset;
        arma::uvec indexA = as<arma::uvec>(index);
        arma::mat X(as<arma::mat>(BminsI[i]));

        arma::vec allPerResult(totalPerS);

        for (int m=0; m<NB; m++){

            if (m % 10 == 0)
                Rcpp::checkUserInterrupt();

            arma::mat thisPermu = allPermutedFC[m];
            arma::mat subsetFC = thisPermu.rows(indexA);
            arma::vec tA(newS);

            // loop through each sample
            for (int j=0; j<newS; j++){
                arma::vec subset = subsetFC.col(j);

                arma::vec pf = solve(X,-subset,arma::solve_opts::fast);
                arma::vec diff = pf - subset;
                tA(j) = sum(diff);

            }

            arma::uvec vI = arma::linspace<arma::uvec>(m*newS, ((m+1)*newS)-1, newS);

            allPerResult.elem(vI) = tA;

            }

       output[i] = allPerResult;
    }

    output.names() = BminsI.names();
    return(output);
}
