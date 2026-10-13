//
// Created by Agent on 2026-04-19.
//

#ifndef GEODA_UNILOSH_H
#define GEODA_UNILOSH_H

#include <vector>

#include "LISA.h"

class GeoDaWeight;

class UniLOSH : public LISA {

    const unsigned long CLUSTER_NOT_SIG;
    const unsigned long CLUSTER_HETEROGENEOUS;
    const unsigned long CLUSTER_HOMOGENEOUS;
    const unsigned long CLUSTER_UNDEFINED;
    const unsigned long CLUSTER_NEIGHBORLESS;

public:
    UniLOSH(int num_obs,
            GeoDaWeight* w,
            const std::vector<double>& data,
            const std::vector<bool>& undefs,
            double significance_cutoff,
            int nCPUs, int permutations,
            const std::string& _permutation_method,
            uint64_t last_seed_used,
            double a = 2.0);

    virtual ~UniLOSH();

    virtual void ComputeLoalSA() ;

    virtual void CalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start);

    virtual void PermCalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start);

    virtual void PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors, std::vector<double>& permutedSA);

    virtual void PermLocalSA(int cnt, int perm, int numNeighbors, const int* permNeighbors, std::vector<double>& permutedSA);

    virtual uint64_t CountLargerSA(int cnt, const std::vector<double>& permutedSA);

    virtual std::vector<int> GetClusterIndicators();

    std::vector<double> GetLocalMean();
    std::vector<double> GetLocalResiduals();

protected:
    std::vector<double> data;
    double a;
    double mean_e;
    std::vector<double> local_mean;
    std::vector<double> local_residuals;
};


#endif //GEODA_UNILOSH_H
