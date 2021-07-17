//
// Created by Xun Li on 6/3/21.
//

#ifndef JSGEODA_PAM_WRAPPER_H
#define JSGEODA_PAM_WRAPPER_H

#include <vector>
#include <string>

class pam_wrapper {
public:
    pam_wrapper(int k,
                const std::vector<std::vector<double> > &data,
                const std::string &distance_method,
                int maxiter,
                const std::string& initializer,
                double fasttol,
                int seed
                );

    virtual ~pam_wrapper();

    double GetCost() { return cost;}
    std::vector<int> GetMedoids() { return medoids; }
    std::vector<int> GetClusters() { return results; }

private:
    double cost;
    std::vector<int> medoids;
    std::vector<int> results;
};
#endif //JSGEODA_PAM_WRAPPER_H
