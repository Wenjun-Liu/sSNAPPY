#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

List ssPertScore_RCPP(const List& BminsI, const NumericMatrix& weightedFC) {
    // double thres = 1e-7;
    int n = BminsI.size();
    int s = weightedFC.ncol();
    CharacterVector genes = rownames(weightedFC);
    std::sort(genes.begin(), genes.end());
    List output(n);
    for(int i=0; i<n; i++){
        Eigen::MatrixXd X(as<Eigen::MatrixXd>(BminsI[i]));
        NumericVector PF(s);
        CharacterVector pathwayG = rownames(BminsI[i]);
        int g = pathwayG.length();
        for (int j=0; j<s; j++){
            NumericVector temp = weightedFC(_,j);
            temp.names() = genes;
            NumericVector v(g);
            for (int m=0; m<g; m++){
                if(std::binary_search(genes.begin(), genes.end(), pathwayG[m])){
                    String index = pathwayG[m];
                    v[m] = temp[index];
                } else {
                    v[m] = 0.0;
                }
            }
            Eigen::VectorXd v1(as<Eigen::VectorXd>(v));
            Eigen::VectorXd pf =X.colPivHouseholderQr().solve(-v1);
            Eigen::VectorXd diff = pf - v1;
            PF[j] = diff.sum();
            PF.names() = colnames(weightedFC);

        }
        output(i) = PF;
        output.names() = BminsI.names();

    }

    return(output);
}



