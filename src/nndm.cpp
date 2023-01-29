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

  Rcpp::NumericMatrix as_rcpp_matrix() const {
    Rcpp::NumericMatrix ret(rows_, columns_);
    for(int i(0); i < rows_; ++i) {
      for(int j(0); j < columns_; ++j) {
        ret(i, j) = operator()(i, j);
      }
    }
    return ret;
  }

private:
  int rows_;
  int columns_;
  MatrixType matrix_;
};

int count_less_than_value(const std::vector<double>& vector,
                          double value){
  int less_than_value(0);
  for(int i=0; i < vector.size(); i++){
    if(vector[i] <= value){
      ++less_than_value;
    }
  }
  return less_than_value;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix nndm_cpp(Rcpp::NumericMatrix X,
                             Rcpp::NumericVector Gjstar,
                             Rcpp::NumericVector Gij,
                             double rmin,
                             double phi,
                             int imin,
                             int jmin,
                             double min_train){

  const std::vector<double> Gij_vec(Rcpp::as<std::vector<double>>(Gij));
  std::vector<double> Gjstar_vec(Rcpp::as<std::vector<double>>(Gjstar));
  const Lightweight_matrix<double> X_cpp(X);
  Lightweight_matrix<int> is_na(X.nrow(), X.ncol());

  for(int i=0; i<is_na.nrow() ;i++){
    for(int j=0; j<is_na.ncol() ;j++){
      is_na(i, j) = static_cast<int>(i == j);
    }
  }

  while( rmin <= phi ) {

    const int Gjstar_less_than_rmin(count_less_than_value(Gjstar_vec, rmin));
    const int Gij_less_than_rmin(count_less_than_value(Gij_vec, rmin));
    const double pcdf = static_cast<double>(Gij_less_than_rmin) / static_cast<double>(Gij_vec.size());
    const double tcdf = static_cast<double>(Gjstar_less_than_rmin - 1) / static_cast<double>(Gjstar_vec.size());

    int X_non_na(0);
    for(int j=0; j<X_cpp.ncol(); j++){
      if(is_na(imin, j) != 1) {
        ++X_non_na;
      }
    }

    const bool in_prop = static_cast<double>(X_non_na) / static_cast<double>(X_cpp.ncol()) > min_train;

    double new_rmin = std::numeric_limits<double>::infinity();

    if( (tcdf >= pcdf) && in_prop ) {
      // Signal that the current element is NA
      is_na(imin, jmin) = 1;

      for(int i=0; i<X_cpp.nrow() ;i++){
        Gjstar_vec[i] = std::numeric_limits<double>::infinity();
        for(int j=0; j<X_cpp.ncol() ;j++){
          if((is_na(i, j) != 1) && X_cpp(i, j) < Gjstar_vec[i]) {
            Gjstar_vec[i] = X_cpp(i, j);
          }
        }
        if(Gjstar_vec[i] >= rmin && Gjstar_vec[i] < new_rmin) {
          new_rmin = Gjstar_vec[i];
          imin = i;
        }
      }
    } else if( Gjstar_vec.size() == Gjstar_less_than_rmin ){
      break;
    } else {
      for(int i=0; i<Gjstar_vec.size(); i++){
        if(Gjstar_vec[i] > rmin && Gjstar_vec[i] < new_rmin) {
          new_rmin = Gjstar_vec[i];
          imin = i;
        }
      }
    }

    rmin = new_rmin;

    for(int j=0; j<X_cpp.ncol(); j++){
      if((is_na(imin, j) != 1) && X_cpp(imin, j) == rmin){
        jmin = j;
        break;
      }
    }
  }

  for(int i=0; i<X.nrow() ;i++){
    for(int j=0; j<X.ncol() ;j++){
      if(is_na(i, j) == 1) {
        X(i, j) = NA_REAL;
      }
    }
  }

  return X;
}
