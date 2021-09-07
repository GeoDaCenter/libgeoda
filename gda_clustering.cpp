#include <iostream>
#include <boost/algorithm/string.hpp>
#include "weights/GeodaWeight.h"
#include "clustering/maxp_wrapper.h"
#include "clustering/redcap_wrapper.h"
#include "clustering/azp_wrapper.h"
#include "clustering/schc_wrapper.h"
#include "clustering/spatial_validation.h"
#include "clustering/joincount_ratio.h"
#include "clustering/make_spatial.h"
#include "GenUtils.h"
#include "gda_data.h"
#include "gda_clustering.h"
#include "gda_interface.h"

const std::vector<std::vector<int> > gda_azp_greedy(int p, GeoDaWeight *w,
                                                     const std::vector<std::vector<double> > &_data,
                                                     const std::string& scale_method,
                                                     int inits,
                                                     const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                     const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                     const std::vector<int>& init_regions,
                                                     const std::string &distance_method,
                                                     int rnd_seed)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    // transform data
    int columns = (int)_data.size();
    std::vector<std::vector<double> > data = _data;
    if (!boost::iequals(scale_method, "raw")) {
        for (int i=0; i<columns; i++) {
            gda_transform_inplace(data[i], scale_method);
        }
    }

    azp_greedy_wrapper azp(p, w, data, inits, min_bounds, max_bounds, init_regions, distance_method,
            rnd_seed);

    return azp.GetClusters();
}

const std::vector<std::vector<int> > gda_azp_sa(int p, GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &_data,
                                                const std::string& scale_method,
                                                int inits,
                                                double cooling_rate,
                                                int sa_maxit,
                                                const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                const std::vector<int>& init_regions,
                                                const std::string &distance_method,
                                                int rnd_seed)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    // transform data
    int columns = (int)_data.size();
    std::vector<std::vector<double> > data = _data;
    if (!boost::iequals(scale_method, "raw")) {
        for (int i=0; i<columns; i++) {
            gda_transform_inplace(data[i], scale_method);
        }
    }

    azp_sa_wrapper azp(p, w, data, inits, cooling_rate, sa_maxit, min_bounds, max_bounds, init_regions, distance_method,
                           rnd_seed);

    return azp.GetClusters();
}

const std::vector<std::vector<int> > gda_azp_tabu(int p, GeoDaWeight *w,
                                                  const std::vector<std::vector<double> > &_data,
                                                  const std::string& scale_method,
                                                  int inits,
                                                  int tabu_length,
                                                  int conv_tabu,
                                                  const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                  const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                  const std::vector<int>& init_regions,
                                                  const std::string &distance_method,
                                                  int rnd_seed)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    // transform data
    int columns = (int)_data.size();
    std::vector<std::vector<double> > data = _data;
    if (!boost::iequals(scale_method, "raw")) {
        for (int i=0; i<columns; i++) {
            gda_transform_inplace(data[i], scale_method);
        }
    }

    azp_tabu_wrapper azp(p, w, data, inits, tabu_length, conv_tabu, min_bounds, max_bounds, init_regions,
            distance_method, rnd_seed);

    return azp.GetClusters();
}

const std::vector<std::vector<int> > gda_maxp_greedy(GeoDaWeight *w,
                                                     const std::vector<std::vector<double> > &_data,
                                                     const std::string& scale_method,
                                                     int iterations,
                                                     const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                     const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                     const std::vector<int>& init_regions,
                                                     const std::string &distance_method,
                                                     int rnd_seed,
                                                     int cpu_threads)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    // transform data
    int columns = (int)_data.size();
    std::vector<std::vector<double> > data = _data;
    if (!boost::iequals(scale_method, "raw")) {
        for (int i=0; i<columns; i++) {
            gda_transform_inplace(data[i], scale_method);
        }
    }

    maxp_greedy_wrapper maxp(w, data, iterations, min_bounds, max_bounds, init_regions, distance_method, rnd_seed, cpu_threads);

    return maxp.GetClusters();
}

const std::vector<std::vector<int> > gda_maxp_sa(GeoDaWeight *w,
                                                 const std::vector<std::vector<double> > &_data,
                                                 const std::string& scale_method,
                                              int iterations,
                                              double cooling_rate,
                                              int sa_maxit,
                                              const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                              const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                              const std::vector<int>& init_regions,
                                              const std::string &distance_method,
                                              int rnd_seed,
                                              int cpu_threads)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    // transform data
    int columns = (int)_data.size();
    std::vector<std::vector<double> > data = _data;
    if (!boost::iequals(scale_method, "raw")) {
        for (int i=0; i<columns; i++) {
            gda_transform_inplace(data[i], scale_method);
        }
    }

    maxp_sa_wrapper maxp(w, data, iterations, cooling_rate, sa_maxit, min_bounds, max_bounds, init_regions,
            distance_method, rnd_seed, cpu_threads);

    return maxp.GetClusters();
}

const std::vector<std::vector<int> > gda_maxp_tabu(GeoDaWeight *w,
                                                   const std::vector<std::vector<double> > &_data,
                                                   const std::string& scale_method,
                                                   int iterations,
                                                   int tabu_length,
                                                   int conv_tabu,
                                                   const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                   const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                   const std::vector<int>& init_regions,
                                                   const std::string &distance_method,
                                                   int rnd_seed,
                                                   int cpu_threads)
{
    std::vector<std::vector<int> > result;

    if (w == 0) return result;

    // transform data
    int columns = (int)_data.size();
    std::vector<std::vector<double> > data = _data;
    if (!boost::iequals(scale_method, "raw")) {
        for (int i=0; i<columns; i++) {
            gda_transform_inplace(data[i], scale_method);
        }
    }

    maxp_tabu_wrapper maxp(w, data, iterations, tabu_length, conv_tabu, min_bounds, max_bounds, init_regions,
            distance_method, rnd_seed, cpu_threads);

    return maxp.GetClusters();
}

const std::vector<std::vector<int> > gda_redcap(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &_data,
                                                const std::string& scale_method,
                                                const std::string &redcap_method,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound,
                                                int rand_seed,
                                                int cpu_threads)
{
    std::vector<std::vector<int> > result;
    unsigned int method = 0;
    if (boost::iequals(redcap_method, "firstorder-singlelinkage")) {
        method = 0;
    } else if  (boost::iequals(redcap_method, "fullorder-completelinkage")){
        method = 1;
    } else if  (boost::iequals(redcap_method, "fullorder-averagelinkage")) {
        method = 2;
    } else if  (boost::iequals(redcap_method, "fullorder-singlelinkage")) {
        method = 3;
    } else if  (boost::iequals(redcap_method, "fullorder-wardlinkage")) {
        method = 4;
    }

    if (w == 0 ||  method > 4) return result;

    //if ((int)k > w->num_obs) return result;

    // transform data
    int columns = (int)_data.size();
    std::vector<std::vector<double> > data = _data;
    if (!boost::iequals(scale_method, "raw")) {
        for (int i=0; i<columns; i++) {
            gda_transform_inplace(data[i], scale_method);
        }
    }
    redcap_wrapper redcap(k, w, data, method, distance_method, bound_vals, min_bound, rand_seed, cpu_threads);
    return redcap.GetClusters();
}

const std::vector<std::vector<int> > gda_skater(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &_data,
                                                const std::string& scale_method,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound,
                                                int rand_seed,
                                                int cpu_threads)
{
    return gda_redcap(k, w, _data, scale_method, "firstorder-singlelinkage", distance_method, bound_vals, min_bound, rand_seed, cpu_threads);
}

const std::vector<std::vector<int> > gda_schc(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &_data,
                                                const std::string& scale_method,
                                                const std::string &linkage_method,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound)
{
    std::vector<std::vector<int> > result;
    unsigned int method = 0;
    if (boost::iequals(linkage_method, "single")) {
        method = 0;
    } else if  (boost::iequals(linkage_method, "complete")){
        method = 1;
    } else if  (boost::iequals(linkage_method, "average")) {
        method = 2;
    } else if  (boost::iequals(linkage_method, "ward")) {
        method = 3;
    }

    if (w == 0 ||  method > 4) return result;

    if ((int)k > w->num_obs) return result;

    // transform data
    int columns = (int)_data.size();
    std::vector<std::vector<double> > data = _data;
    if (!boost::iequals(scale_method, "raw")) {
        for (int i=0; i<columns; i++) {
            gda_transform_inplace(data[i], scale_method);
        }
    }

    schc_wrapper schc(k, w, data, method, distance_method, bound_vals, min_bound);
    return schc.GetClusters();
}

double gda_sumofsquares(const std::vector<double>& vals)
{
    std::vector<double> data = vals;
    return  GenUtils::SumOfSquares(data);
}

double gda_totalsumofsquare(const std::vector<std::vector<double> >& vals)
{
    double ssq = 0;
    for (size_t i=0; i<vals.size(); ++i) {
        std::vector<double> data = vals[i];
        GenUtils::StandardizeData(data);
        double ss = gda_sumofsquares(data);
        ssq += ss;
    }
    return ssq;
}

std::vector<double> gda_withinsumofsquare(const std::vector<std::vector<int> >& solution,
                             const std::vector<std::vector<double> >& _data)
{
    size_t cols = _data.size();

    // standardize data
    std::vector<std::vector<double> > data(cols);
    for (size_t c=0; c<cols; ++c) {
        data[c] = _data[c];
        GenUtils::StandardizeData(data[c]);
    }

    std::vector<double> within_ss;
    for (size_t c=0; c<cols; ++c) {
        double ss = 0;
        for (size_t i=0; i<solution.size(); ++i) {
            std::vector<double> vals;
            for (size_t j = 0; j < solution[i].size(); ++j) {
                size_t r = solution[i][j];
                vals.push_back(data[c][r]);
            }
            ss += gda_sumofsquares(vals);
        }
        within_ss.push_back(ss);
    }

    return within_ss;
}

double gda_totalwithinsumofsquare(const std::vector<std::vector<int> >& solution,
                             const std::vector<std::vector<double> >& _data)
{
    double ssq = 0;
    size_t cols = _data.size();

    // standardize data
    std::vector<std::vector<double> > data(cols);
    for (size_t c=0; c<cols; ++c) {
        data[c] = _data[c];
        GenUtils::StandardizeData(data[c]);
    }

    for (size_t c=0; c<cols; ++c) {
        for (size_t i=0; i<solution.size(); ++i) {
            std::vector<double> vals;
            for (size_t j = 0; j < solution[i].size(); ++j) {
                size_t r = solution[i][j];
                vals.push_back(data[c][r]);
            }
            double ss = gda_sumofsquares(vals);
            ssq += ss;
        }
    }
    return ssq;
}


double gda_betweensumofsquare(const std::vector<std::vector<int> >& solution,
                              const std::vector<std::vector<double> >& data)
{
    double totss = gda_totalsumofsquare(data);
    double totwithiness = gda_totalwithinsumofsquare(solution, data);
    double betweenss = totss - totwithiness;
    return betweenss;
}

ValidationResult gda_spatialvalidation(AbstractGeoDa* geoda, const std::vector<int>& clusters,
        GeoDaWeight *w)
{
    int num_obs = geoda->GetNumObs();

    std::vector<std::vector<int> > groups;
    std::map<int, std::vector<int> > cluster_dict;
    for (int i = 0; i < num_obs; ++i) {
        cluster_dict[clusters[i]].push_back(i);
    }
    std::map<int, std::vector<int> >::iterator it;
    for (it = cluster_dict.begin(); it != cluster_dict.end(); ++it) {
        groups.push_back(it->second);
    }
    std::sort(groups.begin(), groups.end(), GenUtils::less_vectors);

    SpatialValidation sv(num_obs, groups, w, geoda->GetMainMap().records, geoda->GetMainMap().shape_type);

    ValidationResult result;
    result.spatially_constrained = sv.IsSpatiallyConstrained();
    result.fragmentation = sv.GetFragmentation();
    result.cluster_fragmentation = sv.GetFragmentationFromClusters();
    result.cluster_diameter = sv.GetDiameterFromClusters();
    result.cluster_compactness = sv.GetCompactnessFromClusters();


    result.joincount_ratio = joincount_ratio(clusters, w);

    return result;
}

std::vector<JoinCountRatio> gda_joincount_ratio(const std::vector<int>& clusters, GeoDaWeight *w)
{
    return joincount_ratio(clusters, w);
}

JoinCountRatio gda_all_joincount_ratio(const std::vector<JoinCountRatio>& items)
{
    return all_joincount_ratio(items);
}

std::vector<int> gda_makespatial(const std::vector<int>& clusters, GeoDaWeight* w)
{
    int num_obs = w->GetNumObs();

    std::vector<std::vector<int> > groups;
    std::map<int, std::vector<int> > cluster_dict;
    for (int i = 0; i < num_obs; ++i) {
        cluster_dict[clusters[i]].push_back(i);
    }
    std::map<int, std::vector<int> >::iterator it;
    for (it = cluster_dict.begin(); it != cluster_dict.end(); ++it) {
        groups.push_back(it->second);
    }
    std::sort(groups.begin(), groups.end(), GenUtils::less_vectors);

    MakeSpatial ms(num_obs, groups, w);
    ms.Run();

    std::vector<int> result(num_obs);
    std::vector<std::vector<int> > new_groups = ms.GetClusters();
    int n_cluster = new_groups.size();

    for (int i=0; i < n_cluster; i++) {
        int c = i + 1;
        for (int j=0; j<new_groups[i].size(); j++) {
            int idx = new_groups[i][j];
            result[idx] = c;
        }
    }

    return result;
}
