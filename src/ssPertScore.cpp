#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

List ssPertScore(List BminsI, NumericMatrix weightedFC) {
    double thres = 1e-7;
    int n = BminsI.size();
    int s = weightedFC.ncol();
    List output(n);
    for(int i=0; i<n; i++){
        Eigen::MatrixXd X(as<Eigen::MatrixXd>(BminsI[i]));
        double det = X.determinant();
        int g = X.rows();
        // NumericVector PF(s);
        Eigen::MatrixXd PF(g,s);
        if (fabs(det) > thres){
            CharacterVector pathwayG = rownames(BminsI[i]);
            CharacterVector genes = rownames(weightedFC);
            std::sort(genes.begin(), genes.end());
            for (int j=0; j<s; j++){
                NumericVector temp = weightedFC(_,j);
                temp.names() = genes;
                Eigen::VectorXd v(g);
                for (int m=0; m<g; m++){
                    if(std::binary_search(genes.begin(), genes.end(), pathwayG[m])){
                        String index = pathwayG[m];
                        v[m] = temp[index];
                    } else {
                        v[m] = 0.0;
                    }
                }
                // Eigen::VectorXd pf =X.colPivHouseholderQr().solve(-v);
                // Eigen::VectorXd diff = pf - v;
                // PF[s] = diff.sum();
                PF(_,j) =v;

            }
            output(i) = PF;

        }
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
