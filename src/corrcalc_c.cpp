
#include <Rcpp.h>
using namespace Rcpp;

//' @describeIn corrcalc_r c++ application of corrcalc_r
// [[Rcpp::export]]
NumericMatrix corrcalc_c(NumericMatrix matr, int m,
                         NumericVector order_vecti,
                         NumericVector order_vectj)
{
  NumericMatrix out(m, m);

  for (int i1 = 0; i1 < m; i1++) {
    for (int j1 = i1; j1 < m; j1++) {
      int i = order_vecti[i1];
      int j = order_vectj[i1];
      int k = order_vecti[j1];
      int l = order_vectj[j1];

      double matr_ij = matr(i,j);
      double matr_kl = matr(k,l);
      double matr_ik = matr(i,k);
      double matr_il = matr(i,l);
      double matr_jk = matr(j,k);
      double matr_jl = matr(j,l);

      out(i1,j1) =
        (matr_ij*matr_kl/2) * (pow(matr_ik, 2) + pow(matr_il, 2) + pow(matr_jk, 2) + pow(matr_jl, 2)) -
        matr_ij*(matr_ik*matr_il + matr_jk*matr_jl) -
        matr_kl*(matr_ik*matr_jk + matr_il*matr_jl) +
        (matr_ik*matr_jl + matr_il*matr_jk);
    }
  }
  return(out);
}
