#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include "Lightweight_matrix.h"


// a function to calculate L1 or L2 distances
inline double distance(Lightweight_matrix<double> &x, Lightweight_matrix<double> &y, int i, int ii, bool L1)
{
    int nc = x.ncol();
    double dist = 0.0;

    if (L1)
    {
        // could be faster by avoiding () operator and accessing data directly, and
        // vectorised operations but potentially won't be necessary here..
        for (int j = 0; j < nc; j++)
        {
            dist += std::abs(x(i, j) - y(ii, j)); // Manhattan distance
        }
    }
    else
    {
        for (int j = 0; j < nc; j++)
        {
            double diff = x(i, j) - y(ii, j);
            dist += diff * diff;
        }
        dist = std::sqrt(dist); // Euclidean distance
    }

    return dist;
}


// [[Rcpp::export]]
Rcpp::NumericVector similarity_cpp(
    const Rcpp::NumericMatrix &train_mat,
    const Rcpp::NumericMatrix &test_mat,
    const Rcpp::NumericMatrix &rand_mat,
    bool L1 = true
) {

    Lightweight_matrix<double> train(train_mat);
    Lightweight_matrix<double> test(test_mat);
    Lightweight_matrix<double> rand(rand_mat);
    // rbind train and test matrices for baseline calc
    Lightweight_matrix<double> full_mat = rbind(train, test);

    // check for possible errors.
    if (train.ncol() != test.ncol()) Rcpp::stop("Number of columns in train and test are different!");
    if (train.ncol() != rand.ncol()) Rcpp::stop("Number of columns in train and raster layers are different!");

    int nr_test = test.nrow();
    int nr_train = train.nrow();
    int nr_rand = rand.nrow();

    // output vector
    std::vector<double> output(nr_test, 0.0);

    // loop over the random samples
    double sum_dist = 0.0;
    for (int i = 0; i < nr_rand; i++)
    {
        double min_dist(std::numeric_limits<double>::infinity());
        // loop over the training samples
        for (int row = 0; row < full_mat.nrow(); row++)
        {
            // calcualte L1/L2 distance
            double dist_rand = distance(full_mat, rand, row, i, L1);
            // get the min distance
            if (dist_rand < min_dist)
            {
                min_dist = dist_rand;
            }
        }
        sum_dist += min_dist;
    }
    // average min distance of a random sample to training samples
    const double base_dist = sum_dist / static_cast<double>(nr_rand);

    // loop over the test samples
    for (int i = 0; i < nr_test; i++)
    {
        double min_test(std::numeric_limits<double>::infinity());
        // loop over the training samples to get the min distance
        for (int row = 0; row < nr_train; row++)
        {
            // calcualte L1/L2 distance
            double dist = distance(train, test, row, i, L1);
            // get the min distance
            if (dist < min_test)
            {
                min_test = dist;
            }
        }
        // expected distance - test distance
        double test_dist = base_dist - min_test;

        output[i] = test_dist;
    }

    return Rcpp::wrap(output);
}

