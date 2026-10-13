//
// Created by Agent on 2026-04-19.
//

#include <cmath>
#include "../GeoDaSet.h"
#include "../GenUtils.h"
#include "../weights/GeodaWeight.h"
#include "UniLOSH.h"

UniLOSH::~UniLOSH() {

}

UniLOSH::UniLOSH(int num_obs,
                 GeoDaWeight *w,
                 const std::vector<double> &_data,
                 const std::vector<bool> &_undefs,
                 double significance_cutoff,
                 int _nCPUs, int _perm,
                 const std::string& _permutation_method, uint64_t _last_seed,
                 double _a)
: LISA(num_obs, w, _undefs, significance_cutoff, _nCPUs, _perm, _permutation_method, _last_seed),
  CLUSTER_NOT_SIG(0),
  CLUSTER_HETEROGENEOUS(1),
  CLUSTER_HOMOGENEOUS(2),
  CLUSTER_UNDEFINED(5),
  CLUSTER_NEIGHBORLESS(6),
  data(_data),
  a(_a),
  mean_e(0)
{
    labels.push_back("Not significant"); // 0
    labels.push_back("Heterogeneous");   // 1
    labels.push_back("Homogeneous");     // 2
    labels.push_back("");                // 3
    labels.push_back("");                // 4
    labels.push_back("Undefined");       // 5
    labels.push_back("Isolated");        // 6

    colors.push_back("#eeeeee"); // 0
    colors.push_back("#FF00FF"); // 1
    colors.push_back("#00FFFF"); // 2
    colors.push_back("");        // 3
    colors.push_back("");        // 4
    colors.push_back("#464646"); // 5
    colors.push_back("#999999"); // 6

    local_mean.resize(num_obs, 0.0);
    local_residuals.resize(num_obs, 0.0);
    
    if (_perm == 0) {
        calc_significances = false;
    }

    Run();
}

void UniLOSH::ComputeLoalSA() {
    mean_e = 0.0;
    int valid_obs = 0;
    
    // Pass 1: Compute spatially lagged means and local residuals
    for (int i=0; i<num_obs; i++) {
        if (undefs[i] || weights->IsMasked(i) == false) {
            lag_vec[i] = 0;
            lisa_vec[i] = 0;
            cluster_vec[i] = CLUSTER_UNDEFINED;
        } else {
            if (weights->GetNbrSize(i) == 0 || (weights->GetNbrSize(i) == 1 && weights->CheckNeighbor(i, i))) {
                cluster_vec[i] = CLUSTER_NEIGHBORLESS;
            } else {
                const std::vector<long>& nbrs = weights->GetNeighbors(i);
                const std::vector<double>& nbr_w = weights->GetNeighborWeights(i);
                int num_nbrs = nbrs.size();
                
                double sp_lag = 0.0;
                double w_sum = 0.0;
                int valid_nbrs = 0;
                
                for (int j=0; j<num_nbrs; ++j) {
                    if (nbrs[j] != i && nbrs[j] < num_obs && !undefs[nbrs[j]]) {
                        sp_lag += data[nbrs[j]] * nbr_w[j];
                        w_sum += nbr_w[j];
                        valid_nbrs++;
                    }
                }
                
                if (valid_nbrs == 0) {
                    cluster_vec[i] = CLUSTER_NEIGHBORLESS;
                    continue;
                }
                
                if (w_sum > 0 && row_standardize) {
                    sp_lag /= w_sum;
                }
                
                local_mean[i] = sp_lag;
                double diff = std::abs(data[i] - sp_lag);
                if (a == 2.0) {
                    local_residuals[i] = diff * diff;
                } else if (a == 1.0) {
                    local_residuals[i] = diff;
                } else {
                    local_residuals[i] = std::pow(diff, a);
                }
                
                mean_e += local_residuals[i];
                valid_obs++;
            }
        }
    }
    
    if (valid_obs > 0) {
        mean_e /= valid_obs;
    }
    if (mean_e == 0.0) {
        mean_e = 1.0;
    }
    
    // Pass 2: Compute spatially lagged residuals (Hi)
    for (int i=0; i<num_obs; i++) {
        if (undefs[i] || weights->IsMasked(i) == false || weights->GetNbrSize(i) == 0 || (weights->GetNbrSize(i) == 1 && weights->CheckNeighbor(i, i)) || cluster_vec[i] == CLUSTER_NEIGHBORLESS) continue;
        
        const std::vector<long>& nbrs = weights->GetNeighbors(i);
        const std::vector<double>& nbr_w = weights->GetNeighborWeights(i);
        int num_nbrs = nbrs.size();
        
        double lag_e = 0.0;
        double w_sum = 0.0;
        
        for (int j=0; j<num_nbrs; ++j) {
            if (nbrs[j] != i && nbrs[j] < num_obs && !undefs[nbrs[j]]) {
                lag_e += local_residuals[nbrs[j]] * nbr_w[j];
                w_sum += nbr_w[j];
            }
        }
        
        if (w_sum > 0 && row_standardize) {
            lag_e /= w_sum;
        }
        
        lag_vec[i] = lag_e;
        lisa_vec[i] = lag_e / mean_e;
        
        // Initialize clusters (will be filtered by significance later in GetClusterIndicators)
        if (lisa_vec[i] > 1.0) {
            cluster_vec[i] = CLUSTER_HETEROGENEOUS;
        } else {
            cluster_vec[i] = CLUSTER_HOMOGENEOUS;
        }
    }
}

void UniLOSH::CalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start)
{
    GeoDaSet workPermutation(num_obs);
    int max_rand = num_obs - 1;

    for (int cnt = obs_start; cnt <= obs_end; cnt++) {
        if (undefs[cnt] || weights->IsMasked(cnt) == false) {
            sig_cat_vec[cnt] = 6; // undefined
            sig_local_vec[cnt] = -1.0;
            continue;
        }

        int numNeighbors = weights->GetNbrSize(cnt);
        if (weights->CheckNeighbor(cnt, cnt)) {
            numNeighbors -= 1;
        }
        if (numNeighbors <= 0) {
            sig_cat_vec[cnt] = 5; // neighborless cat
            sig_local_vec[cnt] = -1.0;
            continue;
        }

        std::vector<double> nbr_w;
        const std::vector<long>& nbrs = weights->GetNeighbors(cnt);
        const std::vector<double>& raw_nbr_w = weights->GetNeighborWeights(cnt);
        for (size_t j = 0; j < nbrs.size(); ++j) {
            if (nbrs[j] != cnt) {
                nbr_w.push_back(raw_nbr_w[j]);
            }
        }
        std::vector<double> permutedSA(permutations, 0.0);

        for (int perm = 0; perm < permutations; perm++) {
            int rand = 0, newRandom;
            double rng_val;
            while (rand < numNeighbors) {
                rng_val = Gda::ThomasWangHashDouble(seed_start++) * max_rand;
                newRandom = (int)(rng_val < 0.0 ? std::ceil(rng_val - 0.5) : std::floor(rng_val + 0.5));

                if (newRandom != cnt && !workPermutation.Belongs(newRandom) && weights->GetNbrSize(newRandom) > 0) {
                    workPermutation.Push(newRandom);
                    rand++;
                }
            }
            
            double permutedLag = 0.0;
            double w_sum = 0.0;
            for (int cp = 0; cp < numNeighbors; cp++) {
                int nb = workPermutation.Pop();
                if (!undefs[nb]) {
                    permutedLag += local_residuals[nb] * nbr_w[cp];
                    w_sum += nbr_w[cp];
                }
            }
            
            if (w_sum > 0 && row_standardize) {
                permutedLag /= w_sum;
            }
            permutedSA[perm] = permutedLag / mean_e;
        }

        uint64_t countLarger = CountLargerSA(cnt, permutedSA);
        double _sigLocal = (countLarger + 1.0) / (permutations + 1);

        if (_sigLocal <= 0.0001) sig_cat_vec[cnt] = 4;
        else if (_sigLocal <= 0.001) sig_cat_vec[cnt] = 3;
        else if (_sigLocal <= 0.01) sig_cat_vec[cnt] = 2;
        else if (_sigLocal <= 0.05) sig_cat_vec[cnt] = 1;
        else sig_cat_vec[cnt] = 0;

        sig_local_vec[cnt] = _sigLocal;
    }
}

void UniLOSH::PermCalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start)
{
    for (int cnt = obs_start; cnt <= obs_end; cnt++) {
        if (undefs[cnt] || weights->IsMasked(cnt) == false) {
            sig_cat_vec[cnt] = 6; // undefined
            sig_local_vec[cnt] = -1.0;
            continue;
        }

        int numNeighbors = weights->GetNbrSize(cnt);
        if (weights->CheckNeighbor(cnt, cnt)) {
            numNeighbors -= 1;
        }
        if (numNeighbors <= 0) {
            sig_cat_vec[cnt] = 5; // neighborless cat
            sig_local_vec[cnt] = -1.0;
            continue;
        }

        std::vector<double> nbr_w;
        const std::vector<long>& nbrs = weights->GetNeighbors(cnt);
        const std::vector<double>& raw_nbr_w = weights->GetNeighborWeights(cnt);
        for (size_t j = 0; j < nbrs.size(); ++j) {
            if (nbrs[j] != cnt) {
                nbr_w.push_back(raw_nbr_w[j]);
            }
        }
        std::vector<double> permutedSA(permutations, 0.0);

        for (size_t perm = 0; perm < (size_t)permutations; perm++) {
            double permutedLag = 0.0;
            double w_sum = 0.0;
            int* permNeighbors = perm_table[perm];
            for (int cp = 0; cp < numNeighbors; cp++) {
                int nb = permNeighbors[cp];
                if (nb >= cnt) nb = nb + 1;
                if (!undefs[nb]) {
                    permutedLag += local_residuals[nb] * nbr_w[cp];
                    w_sum += nbr_w[cp];
                }
            }
            
            if (w_sum > 0 && row_standardize) {
                permutedLag /= w_sum;
            }
            permutedSA[perm] = permutedLag / mean_e;
        }

        uint64_t countLarger = CountLargerSA(cnt, permutedSA);
        double _sigLocal = (countLarger + 1.0) / (permutations + 1);

        if (_sigLocal <= 0.0001) sig_cat_vec[cnt] = 4;
        else if (_sigLocal <= 0.001) sig_cat_vec[cnt] = 3;
        else if (_sigLocal <= 0.01) sig_cat_vec[cnt] = 2;
        else if (_sigLocal <= 0.05) sig_cat_vec[cnt] = 1;
        else sig_cat_vec[cnt] = 0;

        sig_local_vec[cnt] = _sigLocal;
    }
}

void UniLOSH::PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors,
        std::vector<double>& permutedSA) {
    // Should not be called because we override CalcPseudoP_range
}

void UniLOSH::PermLocalSA(int cnt, int perm, int numNeighbors, const int* permNeighbors,
                                std::vector<double>& permutedSA) {
    // Should not be called because we override PermCalcPseudoP_range
}

uint64_t UniLOSH::CountLargerSA(int cnt, const std::vector<double>& permutedSA)
{
    uint64_t countLarger = 0;
    for (int i=0; i<permutations; ++i) {
        if (permutedSA[i] >= lisa_vec[cnt]) {
            countLarger += 1;
        }
    }

    if (permutations - countLarger <= countLarger) {
        countLarger = permutations - countLarger;
    }
    return countLarger;
}

std::vector<int> UniLOSH::GetClusterIndicators() {
    std::vector<int> clusters(num_obs);
    double cutoff = GetSignificanceCutoff();
    for (int i=0; i<num_obs; i++) {
        if ((sig_local_vec[i] > cutoff || sig_local_vec[i] < 0) &&
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

std::vector<double> UniLOSH::GetLocalMean() {
    return local_mean;
}

std::vector<double> UniLOSH::GetLocalResiduals() {
    return local_residuals;
}
