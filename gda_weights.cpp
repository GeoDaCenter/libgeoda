#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <regex>
#include <filesystem>
#include <map>

#ifdef _MSC_VER
#include <io.h>
#include <fcntl.h>
#else
#include <locale>
#include <clocale>
#endif

#include "weights/VoronoiUtils.h"
#include "weights/PolysToContigWeights.h"
#include "weights/GalWeight.h"
#include "weights/GeodaWeight.h"
#include "SpatialIndAlgs.h"
#include "gda_interface.h"
#include "gda_weights.h"

GeoDaWeight* contiguity_weights(bool is_queen,
                                AbstractGeoDa* geoda,
                                unsigned int order,
                                bool include_lower_order,
                                double precision_threshold)
{
    if (geoda == 0) return 0;

    int num_obs = geoda->GetNumObs();
    GalWeight* poW = new GalWeight;
    poW->num_obs = num_obs;
    poW->is_symmetric = true;
    poW->symmetry_checked = true;

    if (geoda->GetMapType() == gda::POINT_TYP) {
        std::vector<std::set<int> > nbr_map;
        const std::vector<gda::PointContents*>& centroids = geoda->GetCentroids();
        std::vector<double> x(num_obs), y(num_obs);
        for (int i=0; i<num_obs; ++i) {
            x[i] = centroids[i]->x;
            y[i] = centroids[i]->y;
        }
        Gda::VoronoiUtils::PointsToContiguity(x, y, is_queen, nbr_map);
        poW->gal = Gda::VoronoiUtils::NeighborMapToGal(nbr_map);
        //gda::PointsToContiguity(x, y, is_queen, nbr_map);
        //poW->gal = Gda::NeighborMapToGal(nbr_map);
        if (order > 1) {
            Gda::MakeHigherOrdContiguity(order, num_obs, poW->gal, include_lower_order);
        }

    } else if (geoda->GetMapType() == gda::POLYGON) {
        poW->gal = PolysToContigWeights(geoda->GetMainMap(), is_queen, precision_threshold);
        if (order > 1) {
            Gda::MakeHigherOrdContiguity(order, num_obs, poW->gal, include_lower_order);
        }

    } else {
        // line_type not supported yet, should be detected at script side
        delete poW;
        return 0;
    }

    poW->GetNbrStats();
    return (GeoDaWeight*)poW;
}

GeoDaWeight* gda_queen_weights(AbstractGeoDa* geoda,
                               unsigned int order,
                               bool include_lower_order,
                               double precision_threshold)
{
    bool is_queen = true;
    return contiguity_weights(is_queen, geoda, order, include_lower_order, precision_threshold);
}

GeoDaWeight* gda_rook_weights(AbstractGeoDa* geoda,
                               unsigned int order,
                               bool include_lower_order,
                               double precision_threshold)
{
    bool is_queen = false;
    return contiguity_weights(is_queen, geoda, order, include_lower_order, precision_threshold);
}

GeoDaWeight* gda_knn_weights(AbstractGeoDa* geoda, unsigned int k,
                             double power,
                             bool is_inverse,
                             bool is_arc,
                             bool is_mile,
                             const std::string& kernel,
                             double bandwidth,
                             bool adaptive_bandwidth,
                             bool use_kernel_diagonal,
                             const std::string& polyid)
{
    if (geoda == 0) return 0;

    const std::vector<gda::PointContents*>& centroids = geoda->GetCentroids();

    GwtWeight* poW = SpatialIndAlgs::knn_build(centroids, k, is_arc,
                                               is_mile, is_inverse, power,
                                               kernel, bandwidth, adaptive_bandwidth, use_kernel_diagonal);
    poW->GetNbrStats();
    poW->is_symmetric = false;
    poW->symmetry_checked = true;

    return (GeoDaWeight*)poW;
}

GeoDaWeight* gda_knn_weights_sub(AbstractGeoDa* geoda, unsigned int k, int start, int end,
                             double power,
                             bool is_inverse,
                             bool is_arc,
                             bool is_mile,
                             const std::string& kernel,
                             double bandwidth,
                             bool adaptive_bandwidth,
                             bool use_kernel_diagonal,
                             const std::string& polyid)
{
    if (geoda == 0) return 0;

    const std::vector<gda::PointContents*>& centroids = geoda->GetCentroids();

    GwtWeight* poW = SpatialIndAlgs::knn_build_sub(centroids, k, start, end, is_arc,
                                               is_mile, is_inverse, power,
                                               kernel, bandwidth, adaptive_bandwidth, use_kernel_diagonal);
    poW->GetNbrStats();
    poW->is_symmetric = false;
    poW->symmetry_checked = true;

    return (GeoDaWeight*)poW;
}

double gda_min_distthreshold(AbstractGeoDa* geoda, bool is_arc, bool is_mile)
{
    if (geoda == 0) return 0;

    int num_obs = geoda->GetNumObs();

    const std::vector<gda::PointContents*>& centroids = geoda->GetCentroids();

    std::vector<double> x(num_obs), y(num_obs);
    for (int i=0; i<num_obs; ++i) {
        x[i] = centroids[i]->x;
        y[i] = centroids[i]->y;
    }

    double max_1nn_dist = SpatialIndAlgs::find_max_1nn_dist(x, y, is_arc, is_mile);
    return max_1nn_dist;
}

GeoDaWeight* gda_distance_weights(AbstractGeoDa* geoda, double dist_thres,
                                  const std::string& polyid,
                                  double power,
                                  bool is_inverse,
                                  bool is_arc,
                                  bool is_mile,
                                  const std::string& kernel,
                                  bool use_kernel_diagonal)
{
    if (geoda == 0) return 0;

    int num_obs = geoda->GetNumObs();

    const std::vector<gda::PointContents*>& centroids = geoda->GetCentroids();

    std::vector<double> x(num_obs), y(num_obs);
    for (int i=0; i<num_obs; ++i) {
        x[i] = centroids[i]->x;
        y[i] = centroids[i]->y;
    }
    dist_thres = dist_thres * 1.00000; //m_thres_delta_factor
    GwtWeight* poW = SpatialIndAlgs::thresh_build(x, y, dist_thres, power, is_arc, is_mile,
                                                  kernel, use_kernel_diagonal);

    poW->GetNbrStats();
    poW->is_symmetric = kernel.empty() ? true : false;
    poW->symmetry_checked = true;

    return (GeoDaWeight*)poW;
}

GeoDaWeight* gda_load_weights(const char* weights_path)
{
    // Create and install global locale
#ifdef __WIN32__
    w_char_t wstr[1024];
    std::mbstowcs(wstr, file_path, 1024);
    std::ifstream file(wstr);
#else
    std::ifstream file;
    file.open(weights_path, std::ios::in);
#endif

    if (!(file.is_open() && file.good())) {
        return 0;
    }

    // First determine if header line is correct
    // Can be either: int int string string  (type n_obs filename field)
    // or : int (n_obs)

    int line_cnt = 0;
    bool use_rec_order = false;

    int num1 = 0, num2 = 0, num_obs = 0;
    std::string str, dbf_name, key_field, line;

    std::getline(file, str);
    line_cnt += 1;
    std::string header = str;

    // header format: num1 num2 "dbf_name" key_field
    const std::regex r("^([0-9]+)\\s([0-9]+)\\s(.*?|'.*?'|\".*?\")\\s(.*?|'.*?'|\".*?\")$");
    std::smatch what;
    if(std::regex_match(header, what, r)) {
        // The first sub_match is the whole string; the next
        // sub_match is the first parenthesized expression.
        int n_items = what.size();
        if (n_items != 5) {
            return 0;
        }
        std::string num1_str = what.str(1);
        std::string num2_str = what.str(2);
        dbf_name = what.str(3);
        key_field = what.str(4);

        num1 = std::stoi(num1_str);
        num2 = std::stoi(num2_str);
    } else {
        // header format is illegal
        return 0;
    }

    if (num2 == 0) {
        use_rec_order = true;
        num_obs = num1;
    } else {
        num_obs = num2;
        if (key_field.empty() || key_field == "ogc_fid") {
            use_rec_order = true;
        }
    }

    std::map<std::string, int> id_map;

    if (use_rec_order) {
        // we need to traverse through every second line of the file and
        // record the max and min values.  So long as the max and min
        // values are such that num_obs = (max - min) + 1, we will assume
        // record order is valid.
        int min_val = std::numeric_limits<int>::max();
        int max_val = 0;

        while (!file.eof()) {
            int  obs=0, num_neigh=0;
            // get next non-blank line
            str = "";
            while (str.empty() && !file.eof()) {
                getline(file, str);
                line_cnt++;
            }
            if (file.eof()) {
                continue;
            }
            std::stringstream ss (str, std::stringstream::in | std::stringstream::out);
            ss >> obs >> num_neigh;
            if (obs < min_val) {
                min_val = obs;
            } else if (obs > max_val) {
                max_val = obs;
            }
            if (num_neigh > 0) { // ignore the list of neighbors
                // get next non-blank line
                str = "";
                while (str.empty() && !file.eof()) {
                    getline(file, str);
                    line_cnt++;
                }
                if (file.eof()) continue;
            }
        }
        if (max_val - min_val != num_obs - 1) {
            // num_obs doesn't match
            return 0;
        }
        for (int i=0; i<num_obs; i++) {
            std::string iid;
            iid = std::to_string(i+min_val);
            id_map[ iid ] = i;
        }
    } else {
        // using sequential ids 0,1,2,...
        for (int i=0; i<num_obs; i++) {
            std::string str_id;
            str_id = std::to_string(i);
            id_map[ str_id ] = i;
        }
    }

    // create weights
    GalElement* gal = new GalElement[num_obs];
    file.clear();
    file.seekg(0, std::ios::beg); // reset to beginning
    line_cnt = 0;

    getline(file, str); // skip header line
    line_cnt++;

    std::map<std::string, int>::iterator it;
    while (!file.eof()) {
        int gal_obs = 0, num_neigh=0;
        std::string obs;
        // get next non-blank line
        str = "";
        while (str.empty() && !file.eof()) {
            getline(file, str);
            line_cnt++;
        }
        if (file.eof()) continue;
        std::stringstream ss (str, std::stringstream::in | std::stringstream::out);
        ss >> obs >> num_neigh;
        it = id_map.find(obs);
        if (it == id_map.end()) {
            delete [] gal;
            return 0; // observation doesn't match
        }
        gal_obs = (*it).second; // value
        gal[gal_obs].SetSizeNbrs(num_neigh);

        if (num_neigh > 0) {
            // skip next if no neighbors
            // get next non-blank line
            str = "";
            while (str.empty() && !file.eof()) {
                getline(file, str);
                line_cnt++;
            }
            if (file.eof()) continue;
            std::stringstream ss (str, std::stringstream::in | std::stringstream::out);
            for (int j=0; j<num_neigh; j++) {
                std::string neigh;
                ss >> neigh;
                it = id_map.find(neigh);
                if (it == id_map.end()) {
                    delete [] gal;
                    return 0; // observation doesn't match
                }
                long n_id = (*it).second; // value of id_map[neigh];
                gal[gal_obs].SetNbr(j, n_id);
            }
        }
    }

    file.clear();
    if (file.is_open()) {
        file.close();
    }

    GalWeight *rw = new GalWeight();
    rw->gal = gal;
    return rw;
}
