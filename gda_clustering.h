#ifndef __JSGEODSA_GDA_CLUSTERING__
#define __JSGEODSA_GDA_CLUSTERING__

#include <vector>
#include <string>

#include "./weights/GeodaWeight.h"

class AbstractGeoDa;

// APIs of clustering
/**
 *
 * @param p
 * @param w
 * @param data
 * @param inits
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @return
 */
const std::vector<std::vector<int> > gda_azp_greedy(int p, GeoDaWeight *w,
                                                    const std::vector<std::vector<double> > &_data,
                                                    const std::string& scale_method,
                                                    int inits,
                                                    const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                    const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                    const std::vector<int>& init_regions,
                                                    const std::string &distance_method,
                                                    int rnd_seed);

/**
 *
 * @param p
 * @param w
 * @param data
 * @param inits
 * @param cooling_rate
 * @param sa_maxit
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @return
 */
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
                                                int rnd_seed);

/**
 *
 * @param p
 * @param w
 * @param data
 * @param inits
 * @param tabu_length
 * @param conv_tabu
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @return
 */
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
                                                  int rnd_seed);

/**
 *
 * @param w
 * @param data
 * @param iterations
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @param cpu_threads
 * @return
 */
const std::vector<std::vector<int> > gda_maxp_greedy(GeoDaWeight *w,
                                                     const std::vector<std::vector<double> > &_data,
                                                     const std::string& scale_method,
                                                     int iterations,
                                                     const std::vector<std::pair<double, std::vector<double> > >& min_bounds,
                                                     const std::vector<std::pair<double, std::vector<double> > >& max_bounds,
                                                     const std::vector<int>& init_regions,
                                                     const std::string &distance_method,
                                                     int rnd_seed,
                                                     int cpu_threads);


/**
 *
 * @param w
 * @param data
 * @param iterations
 * @param cooling_rate
 * @param sa_maxit
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @param cpu_threads
 * @return
 */
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
                                                 int cpu_threads);

/**
 *
 * @param w
 * @param data
 * @param iterations
 * @param tabu_length
 * @param conv_tabu
 * @param min_bounds
 * @param max_bounds
 * @param init_regions
 * @param distance_method
 * @param rnd_seed
 * @param cpu_threads
 * @return
 */
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
                                                   int cpu_threads);

/**
 *
 * @param k
 * @param w
 * @param data
 * @param redcap_method
 * @param distance_method
 * @param bound_vals
 * @param min_bound
 * @param rand_seed
 * @param cpu_threads
 * @return
 */
const std::vector<std::vector<int> > gda_redcap(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &_data,
                                                const std::string& scale_method,
                                                const std::string &redcap_method,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound,
                                                int rand_seed,
                                                int cpu_threads);

/**
 *
 * @param k
 * @param w
 * @param data
 * @param distance_method
 * @param bound_vals
 * @param min_bound
 * @param rand_seed
 * @param cpu_threads
 * @return
 */
const std::vector<std::vector<int> > gda_skater(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &_data,
                                                const std::string& scale_method,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound,
                                                int rand_seed,
                                                int cpu_threads);


/**
 *
 * @param k
 * @param w
 * @param data
 * @param linkage_method
 * @param distance_method
 * @param bound_vals
 * @param min_bound
 * @param rand_seed
 * @param cpu_threads
 * @return
 */
const std::vector<std::vector<int> > gda_schc(unsigned int k,
                                                GeoDaWeight *w,
                                                const std::vector<std::vector<double> > &_data,
                                                const std::string& scale_method,
                                                const std::string &linkage_method,
                                                const std::string &distance_method,
                                                const std::vector<double>& bound_vals,
                                                double min_bound);

/**
 *
 * @param vals
 * @return
 */
double gda_sumofsquares(const std::vector<double>& vals);

/**
 *
 * @param vals
 * @return
 */
double gda_totalsumofsquare(const std::vector<std::vector<double> >& vals);

/**
 *
 * @param solution
 * @param vals
 * @return
 */
std::vector<double> gda_withinsumofsquare(const std::vector<std::vector<int> >& solution,
                             const std::vector<std::vector<double> >& vals);

/**
 *
 * @param solution
 * @param data
 * @return
 */
double gda_betweensumofsquare(const std::vector<std::vector<int> >& solution,
                              const std::vector<std::vector<double> >& data);



struct Fragmentation {
    int n;
    double entropy;
    double std_entropy;
    double simpson;
    double std_simpson;
    int min_cluster_size;
    int max_cluster_size;
    double mean_cluster_size;
    bool is_spatially_contiguous;
    double fraction;

    Fragmentation() : n(0), entropy(0), std_entropy(0), simpson(0), std_simpson(0),
                      min_cluster_size(0), max_cluster_size(0), mean_cluster_size(0),
                      is_spatially_contiguous(true), fraction(0) {}

    Fragmentation& operator = (const Fragmentation& item) {
        n = item.n;
        entropy = item.entropy;
        std_entropy = item.std_entropy;
        simpson = item.simpson;
        std_simpson = item.std_simpson;
        min_cluster_size = item.min_cluster_size;
        max_cluster_size = item.max_cluster_size;
        mean_cluster_size = item.mean_cluster_size;
        is_spatially_contiguous = item.is_spatially_contiguous;
        return *this;
    }
};

struct Compactness {
    double isoperimeter_quotient;
    double area;
    double perimeter;
    Compactness() : isoperimeter_quotient(0), area(0), perimeter(0) {}

    Compactness& operator = (const Compactness& item) {
        isoperimeter_quotient = item.isoperimeter_quotient;
        area = item.area;
        perimeter = item.perimeter;
        return *this;
    }
};

struct Diameter {
    int steps;
    double ratio;
    Diameter() : steps(0), ratio(0) {}

    Diameter& operator = (const Diameter& item) {
        steps = item.steps;
        ratio = item.ratio;
        return *this;
    }
};

struct JoinCountRatio {
    int cluster;
    int n;
    int totalNeighbors;
    int totalJoinCount;
    double ratio;
    JoinCountRatio(): cluster(0), n(0), totalNeighbors(0), totalJoinCount(0),ratio(0) {}
};

struct ValidationResult {
    bool spatially_constrained;
    Fragmentation fragmentation;
    std::vector<Fragmentation> cluster_fragmentation;
    std::vector<Diameter> cluster_diameter;
    std::vector<Compactness> cluster_compactness;
    std::vector<JoinCountRatio> joincount_ratio;
};

/**
 *
 * @param geoda
 * @param clusters
 * @param w
 * @return
 */
ValidationResult gda_spatialvalidation(AbstractGeoDa* geoda, const std::vector<int>& clusters, GeoDaWeight *w);

/**
 * Make spatially constrained clusters from non-spatially constrained clusters
 *
 * @param clusters
 * @param w
 * @return
 */
std::vector<int> gda_makespatial(const std::vector<int>& clusters, GeoDaWeight* w);

/**
 * 
 * @param items
 * @return
 */
JoinCountRatio gda_all_joincount_ratio(const std::vector<JoinCountRatio>& items);

/**
 * 
 * @param clusters
 * @param w
 * @return
 */
std::vector<JoinCountRatio> gda_joincount_ratio(const std::vector<int>& clusters, GeoDaWeight *w);

#endif

