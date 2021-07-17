//
// Created by Xun Li on 6/3/21.
//
#include <boost/algorithm/string.hpp>
#include "cluster.h"
#include "pam_wrapper.h"
#include "pam.h"

pam_wrapper::pam_wrapper(int k,
            const std::vector<std::vector<double> > &data,
            const std::string &distance_method,
            int maxiter,
            const std::string& initializer,
            double fasttol,
            int seed
)
{
    // get distance matrix
    int n_cols = data.size();
    int num_obs = data[0].size();

    double** matrix = new double*[num_obs];
    int** mask = new int*[num_obs];
    for (int i=0; i<num_obs; ++i) {
        matrix[i] = new double[n_cols];
        mask[i] = new int[n_cols];
        for (int j=0; j<n_cols; ++j) mask[i][j] = 1;
    }
    for (int i=0; i<n_cols; ++i) {
        // the data will be standardized in the caller
        //std::vector<double>& vals = data[i];
        //GenUtils::StandardizeData(vals);
        for (int r=0; r<num_obs; ++r) {
            matrix[r][i] = data[i][r];
        }
    }
    char dist = 'e';
    if (boost::iequals(distance_method, "manhattan")) dist = 'b';
    int transpose = 0; // row wise
    double* weight = new double[n_cols];
    for (int i=0; i<n_cols; ++i) weight[i] = 1.0;

    double** distances = distancematrix(num_obs, n_cols, matrix,  mask, weight, dist, transpose);

    RawDistMatrix dm(distances);

    PAMInitializer* init;
    if (initializer.compare("BUILD")) {
        init = new BUILD(&dm);
    } else {
        init = new LAB(&dm, seed);
    }

    FastPAM pam(num_obs, &dm, init, k, maxiter, fasttol);

    cost = pam.run();

    medoids = pam.getMedoids();
    results = pam.getResults();

    delete init;

    if (weight) delete[] weight;
    if (distances) {
        for (int i = 1; i < num_obs; i++) free(distances[i]);
        free(distances);
    }
    if (matrix) {
        for (int i = 0; i < num_obs; ++i) delete[] matrix[i];
        delete[] matrix;
    }
}