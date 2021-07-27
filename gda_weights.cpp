#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <regex>
#include <map>

#ifdef _MSC_VER
#include <io.h>
#include <fcntl.h>
#else
#include <locale>
#include <clocale>
#endif

#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>

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

GeoDaWeight* gda_load_gal(const char* weights_path, const std::vector<std::string>& id_vec)
{
    std::ifstream file;
    file.open(weights_path, std::ios::in);

    if (!(file.is_open() && file.good())) {
        //std::cout << "not open" << std::endl;
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
            //std::cout << "n_items != 5" << std::endl;
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
        //std::cout << "Error: header format is illegal" << std::endl;
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

    if (!id_vec.empty() && num_obs != id_vec.size()) {
		//std::cout << "Error: the size of id_vec is not equal to the number of observations." << std::endl;
        return 0;
	}

    std::map<std::string, int> id_map;

    if (use_rec_order) {
        // using sequential ids 0,1,2,...
        for (int i=0; i<num_obs; i++) {
            std::string str_id;
            str_id = std::to_string(i);
            id_map[ str_id ] = i;
        }
    } else {
        for (int i=0; i<num_obs; i++) {
            id_map[ id_vec[i] ] = i;
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
            //std::cout << "Error: observation id (" << obs << ") doesn't match input vector of ids." << std::endl;
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
                    //std::cout << "observation doesn't match2" << std::endl;
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
    rw->num_obs = num_obs;
    rw->is_symmetric = true;
    rw->id_field = key_field;
    rw->GetNbrStats();
    return (GeoDaWeight*)rw;
}

GeoDaWeight* gda_load_gwt(const char* weights_path, const std::vector<std::string>& id_vec)
{
    std::ifstream file;
    file.open(weights_path, std::ios::in);

    if (!(file.is_open() && file.good())) {
        //std::cout << "Error: can't open file." << std::endl; 
        return 0;
    }

    // First determine if header line is correct
	// Can be either: int int string string  (type n_obs filename field)
	// or : int (n_obs)
	
	bool use_rec_order = false;
	std::string str;
	getline(file, str);
	std::stringstream ss(str, std::stringstream::in | std::stringstream::out);
	
	int num1 = 0;
	int num2 = 0;
	int num_obs = 0;	
	std::string dbf_name, key_field;
    
    std::string line;
    getline(ss, line);
    std::string header(line);
    
    // header format: num1 num2 "dbf_name" key_field
    const std::regex r("^([0-9]+)\\s([0-9]+)\\s(.*?|'.*?'|\".*?\")\\s(.*?|'.*?'|\".*?\")$");
    std::smatch what;
    if(std::regex_match(header, what, r)) {
        // The first sub_match is the whole string; the next
        // sub_match is the first parenthesized expression.
        int n_items = what.size();
        if (n_items != 5) {
            //std::cout << "Error: header format is illegal" << std::endl;
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
        //std::cout << "Error: header format is illegal" << std::endl;
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
	
	if (!id_vec.empty() && num_obs != id_vec.size()) {
		//std::cout << "Error: the size of id_vec is not equal to the number of observations." << std::endl;
		return 0;
	}

	std::map<std::string, int> id_map;
	if (use_rec_order) {
        // using sequential ids 0,1,2,...
        for (int i=0; i<num_obs; i++) {
            std::string str_id;
            str_id = std::to_string(i);
            id_map[ str_id ] = i;
        }
    } else {
        for (int i=0; i<num_obs; i++) {
            id_map[ id_vec[i] ] = i;
        }
    }

	file.clear();
	file.seekg(0, std::ios::beg); // reset to beginning
	getline(file, str); // skip header line
	// we need to traverse through every line of the file and
	// record the number of neighbors for each observation.
	std::map<std::string, int>::iterator it;
	std::map<std::string, int> nbr_histogram;
	while (!file.eof()) {
		std::string obs1;
		getline(file, str);
		if (!str.empty()) {
			std::stringstream ss (str, std::stringstream::in | std::stringstream::out);
			ss >> obs1;
			
			it = nbr_histogram.find(obs1);
			if (it == nbr_histogram.end()) {
				nbr_histogram[obs1] = 1;
			} else {
				nbr_histogram[obs1] = (*it).second + 1;
			}
		}
	}
	
	GwtElement* gwt = new GwtElement[num_obs];
	file.clear();
	file.seekg(0, std::ios::beg); // reset to beginning
	getline(file, str); // skip header line
	std::map<std::string, int>::iterator it1;
	std::map<std::string, int>::iterator it2;
	int line_num = 1;
	while (!file.eof()) {
		int gwt_obs1, gwt_obs2;
		std::string obs1, obs2;
        double w_val;
		getline(file, str);
		if (!str.empty()) {
			std::stringstream ss(str, std::stringstream::in | std::stringstream::out);
			ss >> obs1 >> obs2 >> w_val;
			it1 = id_map.find(obs1);
			it2 = id_map.find(obs2);
			if (it1 == id_map.end() || it2 == id_map.end()) {
                // id not found
				delete [] gwt;
				return 0;
			}
			gwt_obs1 = (*it1).second; // value
			gwt_obs2 = (*it2).second; // value
            
			if (gwt[gwt_obs1].empty()) {
                gwt[gwt_obs1].alloc(nbr_histogram[obs1]);
            }
            
			gwt[gwt_obs1].Push(GwtNeighbor(gwt_obs2, w_val));
		}
		line_num++;
	}	
	
	if (file.is_open()) file.close();
	
    GwtWeight *rw = new GwtWeight();
    rw->gwt = gwt;
    rw->num_obs = num_obs;
    rw->is_symmetric = false;
    rw->id_field = key_field;
    rw->GetNbrStats();
    return rw;
}

GeoDaWeight* gda_load_swm(const char* weights_path, const std::vector<int>& id_vec)
{
    std::ifstream istream;
    istream.open(weights_path, std::ios::binary|std::ios::in);  // a text file
    
    if (!(istream.is_open() && istream.good())) {
        //std::cout << "Error: can't open file." << std::endl; 
        return 0;
    }
    // first line
    // ID_VAR_NAME;ESRI_SRS\n
    std::string line;
    getline(istream, line, '\n');
    std::string id_name = line.substr(0, line.find(';'));
    
    int swmType = 0; // old
    bool fixed = false;
    
    if (id_name.find("VERSION")==0) {
        swmType = 1;
        // new format: VERSION@10.1;UNIQUEID@FIELD_ID;
        id_name = line.substr(line.find(';')+1, line.size()-1);
        id_name = id_name.substr(0, id_name.find(';'));
        if (id_name.find("UNIQUEID") ==0 ){
            id_name = id_name.substr(id_name.find('@')+1, id_name.size()-1);
        }
        int pos = line.find("FIXEDWEIGHTS@");
        if (pos > 0) {
            std::string fixed_w = line.substr(pos+13, 4);
            if (boost::iequals(fixed_w, "True")) {
                fixed = true;
            }
        }
    }
    
    // NO_OBS length=4
    uint32_t no_obs = 0;
    istream.read((char*)&no_obs, 4); // reads 4 bytes into

    if (!id_vec.empty() && id_vec.size() != no_obs) {
        // input id_vec doesn't match with NumRecords in swm
		//std::cout << "Error: the size of id_vec is not equal to the number of observations." << std::endl;
        return 0;
    }

    boost::unordered_map<int, uint32_t> id_map;
    
    if (id_name != "Unknown" && !id_vec.empty()) {
        for (size_t i=0; i<id_vec.size(); i++) {
            id_map[id_vec[i]] = i;
        }
    } else {
        for (size_t i=0; i<no_obs; i++) {
            id_map[i] = i;
        }
    }
    
    // ROW_STD length = 4
    uint32_t row_std = 0;
    istream.read((char*)&row_std, 4);
    
    std::vector<std::vector<int> > nbr_ids(no_obs);
    std::vector<std::vector<double> > nbr_ws(no_obs);

    for (size_t i=0; i<no_obs; i++) {
        // origin length = 4
        uint32_t origin = 0;
        istream.read((char*)&origin, 4);
        int o_idx = id_map[origin];
        
        if ( id_map.find(o_idx) == id_map.end() ) {
            // WeightsIntegerKeyNotFoundException(o_idx)
            //std::cout << "Error: observation id (" << o_idx << ") not found." << std::endl;
            return 0;
        }
        
        // no_nghs length = 4
        uint32_t no_nghs = 0;
        istream.read((char*)&no_nghs, 4);
        
        if (no_nghs > 0) {
            if (fixed) {
                uint32_t* n_ids = new uint32_t[no_nghs];
                istream.read ((char*)n_ids, sizeof (uint32_t) * no_nghs);
                
                double _w = 0;
                istream.read((char*)&_w, sizeof(double));
                
                double sum_w;
                istream.read((char*)&sum_w, sizeof(double)); // 8
                
                nbr_ids[o_idx].resize(no_nghs);
                nbr_ws[o_idx].resize(no_nghs);
                for (size_t j=0; j<no_nghs; j++) {
                    if ( id_map.find(n_ids[j]) == id_map.end() ) {
                        // WeightsIntegerKeyNotFoundException(o_idx);
                        //std::cout << "Error: observation id (" << n_ids[j] << ") not found." << std::endl;
                        return 0;
                    }
                    nbr_ids[o_idx][j] = id_map[ n_ids[j] ];
                    nbr_ws[o_idx][j] = _w;
                }
                delete[] n_ids;
                
            } else {
                uint32_t* n_ids = new uint32_t[no_nghs];
                istream.read ((char*)n_ids, sizeof (uint32_t) * no_nghs);
                
                double* n_w = new double[no_nghs];
                istream.read ((char*)n_w, sizeof (double) * no_nghs);
                
                double sum_w;
                istream.read((char*)&sum_w, 8);
                
                nbr_ids[o_idx].resize(no_nghs);
                nbr_ws[o_idx].resize(no_nghs);
                for (size_t j=0; j<no_nghs; j++) {
                    if ( id_map.find(n_ids[j]) == id_map.end() ) {
                        // WeightsIntegerKeyNotFoundException(o_idx);
                        //std::cout << "Error: observation id (" << n_ids[j] << ") not found." << std::endl;
                        return 0;
                    }
                    nbr_ids[ o_idx ][j] = id_map[ n_ids[j] ];
                    nbr_ws[ o_idx ][j] = n_w[j];
                }
                
                delete[] n_w;
                delete[] n_ids;
            }
        }
    }
    
    GalElement* gal = new GalElement[no_obs];
    for (size_t i=0; i<no_obs; i++) {
        int no_nghs = nbr_ids[i].size();
        gal[i].SetSizeNbrs(no_nghs);
        std::vector<int>& n_ids = nbr_ids[i];
        std::vector<double>& n_w = nbr_ws[i];
        for (int j=0; j<no_nghs; j++) {
            int nid = n_ids[j];
            gal[ i ].SetNbr(j, nid, n_w[j]);
        }
    }
    
    istream.close();

    GalWeight *rw = new GalWeight();
    rw->gal = gal;
    rw->num_obs = no_obs;
    rw->is_symmetric = false;
    rw->id_field = id_name;
    rw->GetNbrStats();
    return rw;
}