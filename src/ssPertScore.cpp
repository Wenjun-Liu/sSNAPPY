#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

List ssPertScore_RCPP(List BminsI, NumericMatrix weightedFC) {
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



// .ssPertScore <- function(BminsI, weightedFC){
// sapply(names(BminsI), function(x){
//
//     if (abs(det(BminsI[[x]])) > 1e-7){
//         sapply(colnames(weightedFC), function(y){
//             delE  <- weightedFC[,y]
//             delE  <- delE[rownames(BminsI[[x]])]
// # If any of the pathway gene was not expressed, set the ssLogFC to 0
//             delE  <- replace(delE, is.na(delE), 0)
//             PF <- solve(BminsI[[x]], -delE)
//             x <- sum(PF - delE)
//         })
//     } else {
// # if determinant of the pathway topology is not positive, the equation does not have a unique solution
//         x <- NULL
//     }
//
// }, simplify = FALSE)
// }
