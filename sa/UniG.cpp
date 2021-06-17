//
// Created by Xun Li on 10/1/19.
//
#include <math.h>
#include <vector>

#include "../GenUtils.h"
#include "../weights/GeodaWeight.h"

#include "UniG.h"

UniG::~UniG() {

}

UniG::UniG(int num_obs,
           GeoDaWeight *w,
           const std::vector<double> &_data,
           const std::vector<bool> &_undefs,
           double significance_cutoff,
           int _nCPUs, int _perm, const std::string& _permutation_method, uint64_t _last_seed)
: LISA(num_obs, w, _undefs, significance_cutoff, _nCPUs, _perm, _permutation_method, _last_seed),
          CLUSTER_NOT_SIG(0),
          CLUSTER_HIGHHIGH(1),
          CLUSTER_LOWLOW(2),
          CLUSTER_UNDEFINED(3),
          CLUSTER_NEIGHBORLESS(4),
          data(_data),
          sum_x(0)

{
    labels.push_back("Not significant");
    labels.push_back("High-High");
    labels.push_back("Low-Low");
    labels.push_back("Undefined");
    labels.push_back("Isolated");

    colors.push_back("#eeeeee");
    colors.push_back("#FF0000");
    colors.push_back("#0000FF");
    colors.push_back("#464646");
    colors.push_back("#999999");

    G_defined.resize(num_obs, true);

    for (int i=0; i<num_obs; i++) {
        if (!undefs[i])
            sum_x += data[i];
    }

    Run();
}

void UniG::ComputeLoalSA() {
    for (int i=0; i<num_obs; i++) {
        if (undefs[i]) {
            lag_vec[i] = 0;
            lisa_vec[i] = 0;
            cluster_vec[i] = CLUSTER_UNDEFINED;

        } else {
            if (weights->GetNbrSize(i) == 0) {
                cluster_vec[i] = CLUSTER_NEIGHBORLESS;
            } else {
                double lag = 0;
                const std::vector<long>& nbrs = weights->GetNeighbors(i);
                unsigned int nn = 0;
                for (size_t j=0; j<nbrs.size(); ++j) {
                    if (nbrs[j] != i && !undefs[nbrs[j]]) { // not including the value at the location
                        lag += data[ nbrs[j] ];
                        nn += 1;
                    }
                }
                double xd_i = sum_x - data[i];
                if (xd_i == 0) {
                    G_defined[i] = false;
                    cluster_vec[i] = CLUSTER_UNDEFINED;
                    lisa_vec[i] = 0;
                    continue;
                }
                // row-standardize
                lag = lag / nn;
                lisa_vec[i] = lag / xd_i;
            }
        }
    }
    // mean G value
    unsigned int ng = 0;
    double mean_g = 0;
    for (int i=0; i<num_obs; ++i) {
        if (weights->GetNbrSize(i) == 0 || undefs[i] || G_defined[i] == false) {
            continue;
        }
        mean_g += lisa_vec[i];
        ng += 1;
    }
    mean_g = mean_g / ng;

    // assign cluster
    for (int i=0; i<num_obs; ++i) {
        if (weights->GetNbrSize(i) == 0 || undefs[i] || G_defined[i] == false) {
            continue;
        }
        if (lisa_vec[i] >= mean_g) {
            cluster_vec[i] = CLUSTER_HIGHHIGH;
        } else {
            cluster_vec[i] = CLUSTER_LOWLOW;
        }
    }
}

// Used by permutation_method=="lookup"
void UniG::PermLocalSA(int cnt, int perm, int numNeighbors, const int* permNeighbors,
                                std::vector<double>& permutedSA) {
    int validNeighbors = 0;
    double permutedLag = 0;
    // use permutation to compute the lag
    // compute the lag for binary weights
    for (int cp=0; cp<numNeighbors; cp++) {
        int nb = permNeighbors[cp];
        if (nb >= cnt) nb = nb + 1;
        if (!undefs[nb]) {
            permutedLag += data[nb];
            validNeighbors ++;
        }
    }
    double permutedG = permutedLag;

    if (validNeighbors > 0 && row_standardize) {
        permutedLag /= validNeighbors; // row_standardize
        if ((sum_x - data[cnt]) == 0)
            permutedG = 0;
        else
            permutedG = permutedLag / (sum_x - data[cnt]);
    }
    permutedSA[perm] = permutedG;
}

void UniG::PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors,
                          std::vector<double>& permutedSA) {
    int validNeighbors = 0;
    double permutedLag = 0;
    int numNeighbors = permNeighbors.size();

    // use permutation to compute the lag
    // compute the lag for binary weights
    for (int cp=0; cp<numNeighbors; cp++) {
        int nb = permNeighbors[cp];
        if (!undefs[nb] && nb != cnt) {
            permutedLag += data[nb];
            validNeighbors ++;
        }
    }
    double permutedG = permutedLag;

    if (validNeighbors > 0 && row_standardize) {
        permutedLag /= validNeighbors; // row_standardize
        if ((sum_x - data[cnt]) == 0)
            permutedG = 0;
        else
            permutedG = permutedLag / (sum_x - data[cnt]);
    }
    permutedSA[perm] = permutedG;
}

uint64_t UniG::CountLargerSA(int cnt, const std::vector<double>& permutedSA)
{
    uint64_t countLarger = 0;
    for (int i=0; i<permutations; ++i) {
        if (permutedSA[i] >= lisa_vec[cnt]) {
            countLarger += 1;
        }
    }

    // pick the smallest counts
    if (permutations-countLarger < countLarger) {
        countLarger = permutations-countLarger;
    }
    return countLarger;
}

std::vector<int> UniG::GetClusterIndicators() {
    std::vector<int> clusters(num_obs);
    double cutoff = GetSignificanceCutoff();

    for (int i=0; i<num_obs; i++) {
        if ((const unsigned long)cluster_vec[i] == CLUSTER_UNDEFINED &&
            (const unsigned long)cluster_vec[i] == CLUSTER_NEIGHBORLESS)
            continue;

        if (sig_local_vec[i] > cutoff) {
            clusters[i] = CLUSTER_NOT_SIG;
        } else {
            clusters[i] = cluster_vec[i];
        }
    }
    return clusters;
}

