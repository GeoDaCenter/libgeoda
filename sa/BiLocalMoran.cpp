//
// Created by Xun Li on 2021-09-05.
//

#include "../GenUtils.h"
#include "../weights/GeodaWeight.h"
#include "BiLocalMoran.h"


BiLocalMoran::~BiLocalMoran() {

}

BiLocalMoran::BiLocalMoran(int num_obs,
                 GeoDaWeight *w,
                 const std::vector<double> &_data1,
                 const std::vector<double> &_data2,
                 const std::vector<bool> &_undefs,
                 double significance_cutoff,
                 int _nCPUs, int _perm,
                 const std::string& _permutation_method, uint64_t _last_seed)
: LISA(num_obs, w, _undefs, significance_cutoff, _nCPUs, _perm, _permutation_method, _last_seed),
  CLUSTER_NOT_SIG(0),
  CLUSTER_HIGHHIGH(1),
  CLUSTER_LOWLOW(2),
  CLUSTER_LOWHIGH(3),
  CLUSTER_HIGHLOW(4),
  CLUSTER_UNDEFINED(5),
  CLUSTER_NEIGHBORLESS(6),
  data1(_data1),
  data2(_data2)
{
    labels.push_back("Not significant");
    labels.push_back("High-High");
    labels.push_back("Low-Low");
    labels.push_back("Low-High");
    labels.push_back("High-Low");
    labels.push_back("Undefined");
    labels.push_back("Isolated");

    colors.push_back("#eeeeee");
    colors.push_back("#FF0000");
    colors.push_back("#0000FF");
    colors.push_back("#a7adf9");
    colors.push_back("#f4ada8");
    colors.push_back("#464646");
    colors.push_back("#999999");

    GenUtils::StandardizeData(data1, undefs);
    GenUtils::StandardizeData(data2, undefs);
    Run();
}

void BiLocalMoran::ComputeLoalSA() {
    for (int i=0; i<num_obs; i++) {
        if (undefs[i] || weights->IsMasked(i) == false) {
            lag_vec[i] = 0;
            lisa_vec[i] = 0;
            cluster_vec[i] = CLUSTER_UNDEFINED;

        } else {
            if (weights->GetNbrSize(i) == 0) {
                cluster_vec[i] = CLUSTER_NEIGHBORLESS;
            } else {
                double sp_lag = 0;
                const std::vector<long>& nbrs = weights->GetNeighbors(i);
                unsigned int nn = 0;
                for (size_t j=0; j<nbrs.size(); ++j) {
                    if (nbrs[j] != i &&  !undefs[nbrs[j]]) { // not including the value at the location
                        sp_lag += data2[ nbrs[j] ];
                        nn += 1;
                    }
                }
                sp_lag = sp_lag / nn;
                lag_vec[i] = sp_lag;
                lisa_vec[i] = data1[i] * sp_lag;
                // assign the cluster
                if (data1[i] > 0 && sp_lag < 0) cluster_vec[i] = CLUSTER_HIGHLOW;
                else if (data1[i] < 0 && sp_lag > 0) cluster_vec[i] = CLUSTER_LOWHIGH;
                else if (data1[i] < 0 && sp_lag < 0) cluster_vec[i] = CLUSTER_LOWLOW;
                else cluster_vec[i] = CLUSTER_HIGHHIGH; //data1[i] > 0 && Wdata > 0

            }
        }
    }
}

void BiLocalMoran::PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors,
        std::vector<double>& permutedSA) {
    int validNeighbors = 0;
    double permutedLag = 0;
    int numNeighbors = permNeighbors.size();
    // use permutation to compute the lag
    // compute the lag for binary weights
    for (int cp=0; cp<numNeighbors; cp++) {
        int nb = permNeighbors[cp];
        if (!undefs[nb]) {
            permutedLag += data2[nb];
            validNeighbors ++;
        }
    }
    //NOTE: we shouldn't have to row-standardize or
    // multiply by data1[cnt]
    if (validNeighbors > 0 && row_standardize) {
        permutedLag /= validNeighbors;
    }
    const double localMoranPermuted = permutedLag * data1[cnt];
    permutedSA[perm] = localMoranPermuted;
}

void BiLocalMoran::PermLocalSA(int cnt, int perm, int numNeighbors, const int* permNeighbors,
                                std::vector<double>& permutedSA) {
    int validNeighbors = 0;
    double permutedLag = 0;
    // use permutation to compute the lag
    // compute the lag for binary weights
    for (int cp=0; cp<numNeighbors; cp++) {
        int nb = permNeighbors[cp];
        if (nb >= cnt) nb = nb + 1; // self "cnt" should be excluded, index should be adjusted
        if (!undefs[nb]) {
            permutedLag += data2[nb];
            validNeighbors ++;
        }
    }
    //NOTE: we shouldn't have to row-standardize or
    // multiply by data1[cnt]
    if (validNeighbors > 0 && row_standardize) {
        permutedLag /= validNeighbors;
    }
    const double localMoranPermuted = permutedLag * data1[cnt];
    permutedSA[perm] = localMoranPermuted;
}

uint64_t BiLocalMoran::CountLargerSA(int cnt, const std::vector<double>& permutedSA)
{
    uint64_t countLarger = 0;
    for (int i=0; i<permutations; ++i) {
        if (permutedSA[i] >= lisa_vec[cnt]) {
            countLarger += 1;
        }
    }

    // pick the smallest counts
    if (permutations-countLarger <= countLarger) {
        countLarger = permutations-countLarger;
    }
    return countLarger;
}

std::vector<int> BiLocalMoran::GetClusterIndicators() {
    std::vector<int> clusters(num_obs);
    double cuttoff = GetSignificanceCutoff();
    for (int i=0; i<num_obs; i++) {
        if (sig_local_vec[i] > cuttoff &&
                (const unsigned long)cluster_vec[i] != CLUSTER_UNDEFINED &&
                (const unsigned long)cluster_vec[i] != CLUSTER_NEIGHBORLESS)
        {
            clusters[i] = CLUSTER_NOT_SIG;
        } else {
            clusters[i] = cluster_vec[i];
        }
    }
    return clusters;
}

