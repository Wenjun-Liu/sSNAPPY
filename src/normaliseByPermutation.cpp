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

         output[i] = permutedFC;


    }

    return(output);
}

List permutedPertScore_RCPP(const List& BminsI, arma::mat logCPM, int NB, int sEachp) {

    // number of pathways
    int pathway = BminsI.length();
    // number of samples
    int s = logCPM.n_cols;
    // number of matching pairs
    int p = s/sEachp;

    // generate empty list the same length as pahtway topologies as output
    List output(NB);

    // indices for sslogFC computation. seq(0, number of sample -1, by = samples each patient)
    arma::uvec indices = arma::regspace<arma::uvec>(0, sEachp, s-1);

    List permutedFC(NB);

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

        permutedFC[i] = permutedFC;

        }

    for (int i=0; i<pathway; i++){

        }
    List FC(s);

    for (int j=0; j<s; j++){
        // extract the ssFC of the given sample
        NumericVector temp = weightedFC(_,j);
        // set the names of the fc vector to be expressed genes
        temp.names() = rownames(weightedFC);
        FC[j] = temp;
    }

    for (int i=0; i<p; i++){

        Eigen::MatrixXd X(as<Eigen::MatrixXd>(BminsI[i]));

        // pathwayG is the genes contained in this given pathway
        CharacterVector pathwayG = rownames(BminsI[i]);

        NumericVector tA(s);

        // loop through each sample
        for (int j=0; j<s; j++){
            NumericVector ssFC = FC[j];
            Eigen::VectorXd subset(as<Eigen::VectorXd>(ssFC[pathwayG]));

            Eigen::VectorXd pf = X.colPivHouseholderQr().solve(-subset);
            Eigen::VectorXd diff = pf - subset;
            tA[j] = diff.sum();

        }
        tA.names() = colnames(weightedFC);
        output[i] = tA;
    }

    output.names() = BminsI.names();
    return(output);
}
