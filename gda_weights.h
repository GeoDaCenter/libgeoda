#ifndef __JSGEODSA_GDA_WEIGHTS__
#define __JSGEODSA_GDA_WEIGHTS__

#include <string>

class AbstractGeoDa;

class GeoDaWeight;

// APIs of weights creation
/**
 *
 * @param geoda
 * @param polyid
 * @param order
 * @param include_lower_order
 * @param precision_threshold
 * @return
 */
GeoDaWeight* gda_queen_weights(AbstractGeoDa* geoda,
                               unsigned int order,
                               bool include_lower_order,
                               double precision_threshold);

GeoDaWeight* gda_rook_weights(AbstractGeoDa* geoda,
                               unsigned int order,
                               bool include_lower_order,
                               double precision_threshold);

GeoDaWeight* gda_knn_weights(AbstractGeoDa* geoda, unsigned int k,
                             double power,
                             bool is_inverse,
                             bool is_arc,
                             bool is_mile,
                             const std::string& kernel,
                             double bandwidth,
                             bool adaptive_bandwidth,
                             bool use_kernel_diagonals,
                             const std::string& polyid);

GeoDaWeight* gda_knn_weights_sub(AbstractGeoDa* geoda, unsigned int k, int start, int end,
                             double power,
                             bool is_inverse,
                             bool is_arc,
                             bool is_mile,
                             const std::string& kernel,
                             double bandwidth,
                             bool adaptive_bandwidth,
                             bool use_kernel_diagonals,
                             const std::string& polyid);

double gda_min_distthreshold(AbstractGeoDa* geoda, bool is_arc, bool is_mile);

GeoDaWeight* gda_distance_weights(AbstractGeoDa* geoda, double dist_thres,
                                  const std::string& polyid,
                                  double power,
                                  bool is_inverse,
                                  bool is_arc,
                                  bool is_mile,
                                  const std::string& kernel,
                                  bool use_kernel_diagonals);

GeoDaWeight* gda_load_gal(const char* weights_path, const std::vector<std::string>& id_vec = std::vector<std::string>());

GeoDaWeight* gda_load_gwt(const char* weights_path, const std::vector<std::string>& id_vec = std::vector<std::string>());

GeoDaWeight* gda_load_swm(const char* weights_path, const std::vector<int>& id_vec = std::vector<int>());
#endif
