#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".convert_geno")]]
IntegerMatrix convert_geno(IntegerMatrix g, IntegerMatrix codes)
{
    int n_mar = g.rows();
    int n_ind = g.cols();
    int n_code = codes.cols();
    IntegerMatrix result(n_mar,n_ind);

    if(codes.rows() != n_mar)
        throw std::invalid_argument("nrow(g) != nrow(codes)");

    for(int j=0; j<n_ind; j++) {
        for(int i=0; i<n_mar; i++) {
            if(g(i,j)==NA_INTEGER)
                result(i,j)=NA_INTEGER;
            else {
                bool unassigned=true;
                for(int k=0; k<n_code; k++) {
                    if(g(i,j)==codes(i,k)) {
                        unassigned=false;
                        result(i,j) = k+1;
                    }
                }
                if(unassigned) result(i,j) = 0;
            }
        }
    }

    return result;
}

/*** R
convert_geno <-
function(geno, codes)
{
    u <- unique(c(geno[!is.na(geno)], "H"))
    d <- dimnames(geno)
    for(i in seq(along=u)) {
        geno <- gsub(u[i], as.character(i), geno)
        codes <- gsub(u[i], as.character(i), codes)
    }
    geno <- matrix(as.numeric(geno), nrow=nrow(geno))
    codes <- matrix(as.numeric(codes), nrow=nrow(codes))

    geno <- .convert_geno(geno, codes)
    dimnames(geno) <- d

    geno
}
***/
