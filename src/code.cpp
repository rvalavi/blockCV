#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix nndm_cpp(NumericMatrix X,
                       NumericVector Gjstar,
                       NumericVector Gij,
                       double rmin,
                       double phi,
                       int jmin,
                       int kmin,
                       double min_train){

  double Gst_len = Gjstar.size();
  double Gij_len = Gij.size();

  while( rmin <= phi ) {

    LogicalVector Gst_lte_rmin = Gjstar <= rmin;
    LogicalVector Gij_lte_rmin = Gij <= rmin;
    LogicalVector Gst_ht_rmin = Gjstar > rmin;

    double tcdf = (sum(Gst_lte_rmin) - 1) / Gst_len;
    double pcdf = sum(Gij_lte_rmin) / Gij_len;
    bool tcdf_hte_pcdf = tcdf >= pcdf;

    NumericVector one_row = X(jmin, _);
    LogicalVector index = is_na(one_row);
    double n_col = X.ncol();
    bool in_prop = sum(!index) / n_col > min_train;

    if( tcdf_hte_pcdf && in_prop ) {

      // replace the element with NA
      X(jmin, kmin) = NA_REAL;
      for(int i=0; i<X.nrow() ;i++){
        NumericVector row_i = X( i , _ );
        LogicalVector index = is_na(row_i);
        NumericVector Xi = row_i[!index];
        Gjstar[i] = min(Xi);
      }

      NumericVector Gst_hte_rmin = Gjstar[Gjstar >= rmin];
      rmin = min(Gst_hte_rmin);

      for(int i=0; i<Gjstar.size(); i++){
        if(Gjstar[i] == rmin){
          jmin = i;
          break;
        }
      }
      // this can be coverted to for loop for j and then for for i?
      NumericVector row_jmin = X( jmin, _ );
      for(int i=0; i<row_jmin.size(); i++){
        if(row_jmin[i] == rmin){
          kmin = i;
          break;
        }
      }

    } else if( sum(Gst_ht_rmin) == 0 ){
      break;
    } else {
      NumericVector Gst_h_rmin = Gjstar[Gst_ht_rmin];
      rmin = min(Gst_h_rmin);

      for(int i=0; i<Gjstar.size(); i++){
        if(Gjstar[i] == rmin){
          jmin = i;
          break;
        }
      }
      NumericVector row_jmin = X( jmin, _ );
      for(int i=0; i<row_jmin.size(); i++){
        if(row_jmin[i] == rmin){
          kmin = i;
          break;
        }
      }
    }

  }

  return X;
}
