#ifndef IGHTWEIGHT_MATRIX_H
#define IGHTWEIGHT_MATRIX_H

#include <vector>
#include <algorithm>

// A light weight matrix class for faster calculation with row-major opertaions
template<typename T>
class Lightweight_matrix {
private:
    using MatrixType = std::vector<T>;
public:

    explicit Lightweight_matrix(int rows, int columns):
        rows_(rows), columns_(columns), matrix_(rows * columns) {}

    template<typename Matrix>
    explicit Lightweight_matrix(Matrix matrix):
        rows_(matrix.nrow()), columns_(matrix.ncol()), matrix_(matrix.nrow() * matrix.ncol()) {
        for(int i = 0; i < rows_; ++i) {
            for(int j = 0; j < columns_; ++j) {
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

    T* data() {
        return matrix_.data();
    }

    const T* data() const {
        return matrix_.data();
    }

private:
    int rows_;
    int columns_;
    MatrixType matrix_;
};


// rbind implementation for Lightweight_matrix
template<typename T>
Lightweight_matrix<T> rbind(const Lightweight_matrix<T>& A, const Lightweight_matrix<T>& B)
{
    const int ncol = A.ncol();
    const int nrow_A = A.nrow();
    const int nrow_B = B.nrow();

    if (ncol != B.ncol())
    {
        throw std::invalid_argument("Matrices must have the same number of columns");
    }

    const int total_rows = nrow_A + nrow_B;
    Lightweight_matrix<T> result(total_rows, ncol);

    std::copy(A.data(), A.data() + nrow_A * ncol, result.data());
    std::copy(B.data(), B.data() + nrow_B * ncol, result.data() + nrow_A * ncol);

    return result;
}

#endif /* IGHTWEIGHT_MATRIX_H */
