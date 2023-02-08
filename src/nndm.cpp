#include <Rcpp.h>

#include <limits>
#include <vector>

template<typename T>
class Lightweight_matrix {
private:
  using MatrixType = std::vector<T>;
public:

  explicit Lightweight_matrix(int rows, int columns):
    rows_(rows), columns_(columns), matrix_(rows * columns) {}

  template<typename Matrix>
  explicit Lightweight_matrix(Matrix matrix):
    Lightweight_matrix(matrix.nrow(), matrix.ncol()) {
    for(int i(0); i < rows_; ++i) {
      for(int j(0); j < columns_; ++j) {
        operator()(i, j) = static_cast<T>(matrix(i, j));
      }
    }
  }

  T operator()(int i, int j) const {
    return matrix_[i * columns_ + j];
  }

  T& operator()(int i, int j) {
    return matrix_[i * columns_ + j];
  }

  int ncol() const {
    return columns_;
  }

  int nrow() const {
    return rows_;
  }

private:
  int rows_;
  int columns_;
  MatrixType matrix_;
};

// [[Rcpp::export]]
Rcpp::NumericMatrix nndm_cpp(Rcpp::NumericMatrix X,
                             Rcpp::NumericVector Gij,
                             double phi,
                             double min_train){

  // Value representing NAs (NA_REAL does not work well with custom structures)
  const double na_value(-1.);

  // C++ representation of Gij
  const std::vector<double> Gij_vec(Rcpp::as<std::vector<double>>(Gij));

  // C++ representation of X
  Lightweight_matrix<double> X_cpp(X);
  for(int i=0; i<X_cpp.nrow() ;i++){
    X_cpp(i, i) = na_value;
  }

  // Entry i of this vector contains the index at which Gjstar attains the minimum
  // along row i of X
  std::vector<int> Gjstar_indices(X.ncol());
  double rmin = std::numeric_limits<double>::infinity();
  int imin;
  for(int i=0; i<X_cpp.nrow() ;i++){
    double minimum_i = std::numeric_limits<double>::infinity();
    for(int j=0; j<X_cpp.ncol() ;j++){
      if(i != j && X_cpp(i, j) < minimum_i) {
        minimum_i = X_cpp(i, j);
        Gjstar_indices[i] = j;
      }
    }
    if(minimum_i < rmin) {
      rmin = minimum_i;
      imin = i;
    }
  }

  while( rmin <= phi ) {

    // Number of entries of Gjstar that are less than rmin
    int Gjstar_less_than_rmin(0);
    for(int i=0; i < Gjstar_indices.size(); i++){
      if(X_cpp(i, Gjstar_indices[i]) <= rmin){
        ++Gjstar_less_than_rmin;
      }
    }

    // Number of entries of Gij that are less than rmin
    int Gij_less_than_rmin(0);
    for(int i=0; i < Gij_vec.size(); i++){
      if(Gij_vec[i] <= rmin){
        ++Gij_less_than_rmin;
      }
    }
    const double pcdf = static_cast<double>(Gij_less_than_rmin) / static_cast<double>(Gij_vec.size());
    const double tcdf = static_cast<double>(Gjstar_less_than_rmin - 1) / static_cast<double>(Gjstar_indices.size());

    // Number of non-NA values in X[imin, ]
    int X_non_na(0);
    for(int j=0; j<X_cpp.ncol(); j++){
      if(X_cpp(imin, j) != na_value) {
        ++X_non_na;
      }
    }

    const bool in_prop = static_cast<double>(X_non_na) / static_cast<double>(X_cpp.ncol()) > min_train;

    // Will end up containing the updated rmin value
    double new_rmin = std::numeric_limits<double>::infinity();

    if( (tcdf >= pcdf) && in_prop ) {
      // Signal that the current element is NA
      X_cpp(imin, Gjstar_indices[imin]) = na_value;

      double running_minimum = std::numeric_limits<double>::infinity();
      for(int j=0; j<X_cpp.ncol() ;j++){
        if((X_cpp(imin, j) != na_value) && X_cpp(imin, j) < running_minimum) {
          running_minimum = X_cpp(imin, j);
          Gjstar_indices[imin] = j;
        }
      }

      for(int i=0; i<X_cpp.nrow() ;i++){
        if(X_cpp(i, Gjstar_indices[i]) >= rmin && X_cpp(i, Gjstar_indices[i]) < new_rmin) {
          new_rmin = X_cpp(i, Gjstar_indices[i]);
          imin = i;
        }
      }
    } else if( Gjstar_indices.size() == Gjstar_less_than_rmin ){
      break;
    } else {
      for(int i=0; i<Gjstar_indices.size(); i++){
        if(X_cpp(i, Gjstar_indices[i]) > rmin && X_cpp(i, Gjstar_indices[i]) < new_rmin) {
          new_rmin = X_cpp(i, Gjstar_indices[i]);
          imin = i;
        }
      }
    }

    // Update rmin
    rmin = new_rmin;
  }

  // Update X before returning, since we were working on C++ object
  for(int i=0; i<X.nrow() ;i++){
    for(int j=0; j<X.ncol() ;j++){
      if(X_cpp(i, j) == na_value) {
        X(i, j) = NA_REAL;
      }
    }
  }

  return X;
}
