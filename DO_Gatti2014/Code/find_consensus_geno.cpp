#include <Rcpp.h>
using namespace Rcpp;

int consensus_geno(IntegerVector g)
{
    int n = g.size();

    // find maximum value
    int maxval = 0;
    int n_missing = 0;
    for(int i=0; i<n; i++) {
        if(g[i] == NA_INTEGER) n_missing++;
        else if(g[i] > maxval) maxval=g[i];
    }
    // return NA if 1/2 or more missing
    if(n_missing*2 >= n)
        return NA_INTEGER;

    // count unique values
    IntegerVector tab(maxval+1);
    for(int i=0; i<n; i++)
        if(g[i] != NA_INTEGER)
            tab[g[i]]++;

    // find maximum
    int maxtab=0;
    int themax=NA_INTEGER;
    for(int i=0; i<=maxval; i++) {
        if(tab[i] > maxtab) {
            maxtab=tab[i];
            themax=i;
        }
        else if(tab[i] == maxtab) {
            themax=NA_INTEGER;
        }
    }
    return themax;
}

// [[Rcpp::export(".find_consensus_geno")]]
IntegerVector find_consensus_geno(IntegerMatrix g)
{
    int n_mar = g.rows();
    IntegerVector result(n_mar);

    for(int i=0; i<n_mar; i++)
        result[i] = consensus_geno(g(i,_));

    return result;
}

/*** R
find_consensus_geno <-
function(mat)
{
    u <- unique(mat[!is.na(mat)])
    for(i in seq(along=u))
        mat <- gsub(u[i], as.character(i), mat)
    mat <- matrix(as.numeric(mat), nrow=nrow(mat))

    result <- .find_consensus_geno(mat)
    result <- as.character(result)
    for(i in seq(along=u))
        result <- gsub(as.character(i), u[i], result)

    result
}
***/
