#include <algorithm>
#include <queue>
#include <iterator>
#include <stack>

#include "../DataUtils.h"
#include "azp.h"

#ifndef __NO_THREAD__
#ifndef __USE_PTHREAD__
#include <boost/system/config.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#else
#include <pthread.h>

        struct maxp_thread_args {
            MaxpRegion *maxp;
            int start;
            int end;
        };

        void* maxp_thread_helper_construction(void* voidArgs)
        {
            maxp_thread_args *args = (maxp_thread_args*)voidArgs;
            args->maxp->RunConstructionRange(args->start, args->end);
            return 0;
        }

        void* maxp_thread_helper_localimprovement(void* voidArgs)
        {
            maxp_thread_args *args = (maxp_thread_args*)voidArgs;
            args->maxp->RunLocalImprovementRange(args->start, args->end);
            return 0;
        }
#endif
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ZoneControl
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ZoneControl::ZoneControl(const std::vector<double>& in_data)
: data(in_data)
{
}

ZoneControl::~ZoneControl()
{
}

void ZoneControl::AddControl(Operation op, Comparator cmp, const double& val)
{
    operations.push_back(op);
    comparators.push_back(cmp);
    comp_values.push_back(val);
}

bool ZoneControl::CheckRemove(int area, boost::unordered_map<int, bool>& candidates)
{
    bool is_valid = true; // default true since no check will yield good cands
    boost::unordered_map<int, bool>::iterator it;
    for (size_t i=0;  i< comparators.size(); ++i) {
        if (comparators[i] != MORE_THAN) {
            continue;
        }

        // get zone value for comparison
        double zone_val = 0;
        if (operations[i] == SUM) {
            double sum = 0;
            for (it=candidates.begin(); it!=candidates.end(); ++it) {
                sum += data[ it->first ];
            }
            sum -= data[area];
            zone_val = sum;
        } else if (operations[i] == MEAN) {
            double sum = 0;
            for (it=candidates.begin(); it!=candidates.end(); ++it) {
                sum += data[it->first];
            }
            sum -= data[area];
            double mean = sum / (double) (candidates.size() - 1);
            zone_val = mean;
        } else if (operations[i] == MAX) {
            double max = data[candidates[0]];
            for (it=candidates.begin(); it!=candidates.end(); ++it) {
                if (max < data[it->first] && it->first != area) {
                    max = data[it->first];
                }
            }
            zone_val = max;
        } else if (operations[i] == MIN) {
            double min = data[candidates[0]];
            for (it=candidates.begin(); it!=candidates.end(); ++it) {
                if (min > data[it->first] && it->first != area) {
                    min = data[it->first];
                }
            }
            zone_val = min;
        }

        // compare zone value
        if (comparators[i] == MORE_THAN) {
            if (zone_val <= comp_values[i]) {
                return false;
            }
        }
    }
    return is_valid;
}

bool ZoneControl::CheckAdd(int area, boost::unordered_map<int, bool>& candidates)
{
    bool is_valid = true; // default true since no check will yield good cands
    boost::unordered_map<int, bool>::iterator it;
    for (size_t i=0;  i< comparators.size(); ++i) {
        if (comparators[i] != LESS_THAN) {
            continue;
        }

        // get zone value for comparison
        double zone_val = 0;
        if (operations[i] == SUM) {
            double sum = 0;
            for (it=candidates.begin(); it!=candidates.end(); ++it) {
                sum += data[ it->first ];
            }
            sum += data[area];
            zone_val = sum;
        } else if (operations[i] == MEAN) {
            double sum = 0;
            for (it=candidates.begin(); it!=candidates.end(); ++it) {
                sum += data[it->first];
            }
            sum += data[area];
            double mean = sum / (double) (candidates.size() + 1);
            zone_val = mean;
        } else if (operations[i] == MAX) {
            double max = data[candidates[0]];
            for (it=candidates.begin(); it!=candidates.end(); ++it) {
                if (max < data[it->first]) {
                    max = data[it->first];
                }
            }
            if (max < data[area]) {
                max = data[area];
            }
            zone_val = max;
        } else if (operations[i] == MIN) {
            double min = data[candidates[0]];
            for (it=candidates.begin(); it!=candidates.end(); ++it) {
                if (min > data[it->first]) {
                    min = data[it->first];
                }
            }
            if (min > data[area]) {
                min = data[area];
            }
            zone_val = min;
        }

        // compare zone value
        if (comparators[i] == LESS_THAN) {
            if (zone_val >= comp_values[i]) {
                return false;
            }
        }
    }
    return is_valid;
}

double ZoneControl::getZoneValue(int i, boost::unordered_map<int, bool>& candidates)
{
    // get zone value for comparison
    double zone_val = 0;
    boost::unordered_map<int, bool>::iterator it;
    if (operations[i] == SUM) {
        double sum = 0;
        for (it=candidates.begin(); it!=candidates.end(); ++it) {
            sum += data[ it->first ];
        }
        zone_val = sum;
    } else if (operations[i] == MEAN) {
        double sum = 0;
        for (it=candidates.begin(); it!=candidates.end(); ++it) {
            sum += data[it->first];
        }
        double mean = sum / (double) candidates.size();
        zone_val = mean;
    } else if (operations[i] == MAX) {
        double max = data[candidates[0]];
        for (it=candidates.begin(); it!=candidates.end(); ++it) {
            if (max < data[it->first]) {
                max = data[it->first];
            }
        }
        zone_val = max;
    } else if (operations[i] == MIN) {
        double min = data[candidates[0]];
        for (it=candidates.begin(); it!=candidates.end(); ++it) {
            if (min > data[it->first]) {
                min = data[it->first];
            }
        }
        zone_val = min;
    }
    return zone_val;
}

bool ZoneControl::SatisfyLowerBound(boost::unordered_map<int, bool>& candidates)
{
    bool is_valid = true; // default true since no check will yield good cands
    boost::unordered_map<int, bool>::iterator it;

    for (size_t i=0;  i< comparators.size(); ++i) {
        if (comparators[i] != MORE_THAN) {
            continue;
        }

        // get zone value for comparison
        double zone_val = getZoneValue(i, candidates);
        // compare zone value
        if (comparators[i] == MORE_THAN) {
            if (zone_val < comp_values[i]) {
                return false; // not yet satisfy lower bound
            }
        }
    }
    return is_valid;
}

bool ZoneControl::CheckBound(boost::unordered_map<int, bool>& candidates)
{
    bool is_valid = true; // default true since no check will yield good cands
    boost::unordered_map<int, bool>::iterator it;

    for (size_t i=0;  i< comparators.size(); ++i) {

        // get zone value for comparison
        double zone_val = getZoneValue(i, candidates);

        // compare zone value
        if (comparators[i] == MORE_THAN) {
            if (zone_val < comp_values[i]) {
                return false; // not yet satisfy lower bound
            }
        } else if (comparators[i] == LESS_THAN) {
            if (zone_val > comp_values[i]) {
                return false; // not yet satisfy lower bound
            }
        }
    }
    return is_valid;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// AreaManager
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

AreaManager::AreaManager(int _n, int _m, GalElement* const _w, double** _data, DistMatrix* const _dist_matrix)
: n(_n), m(_m), w(_w), dist_matrix(_dist_matrix), data(_data)
{
}

double AreaManager::returnDistance2Area(int i, int j)
{
    return dist_matrix->getDistance(i, j);
}

std::vector<double> AreaManager::getDataAverage(const std::set<int>& areaList)
{
    return std::vector<double>();
}

double AreaManager::getDistance2Region(int area, int region, REGION_AREAS& regions)
{
    std::vector<double> d(m,0);
    for (int i=0; i<m; ++i) d[i] = data[area][i];

    // get centroid of region
    if (region_centroids.find(region) == region_centroids.end()) {
        updateRegionCentroids(region, regions);
    }
    std::vector<double>& centroidRegion = region_centroids[region];

    // get distance between area and region
    double dist = DataUtils::EuclideanDistance(d, centroidRegion);

    return dist;
}

void AreaManager::updateRegionCentroids(int region, REGION_AREAS& regions)
{
    boost::unordered_map<int, bool>& areaList = regions[region];
    std::vector<double> centroid(m, 0);
    boost::unordered_map<int, bool>::iterator it;
    for (it = areaList.begin(); it != areaList.end(); ++it) {
        int area_id = it->first;
        for (int j=0; j<m; ++j) {
            centroid[j] += data[area_id][j];
        }
    }
    for (int j=0; j<m; ++j) {
        centroid[j] /= (double)areaList.size();
    }
    region_centroids[region] = centroid;
}

////////////////////////////////////////////////////////////////////////////////
////// RegionMaker
////////////////////////////////////////////////////////////////////////////////
RegionMaker::RegionMaker(int _p, GalElement* const _w,
                         double** _data, // row-wise
                         RawDistMatrix* _dist_matrix,
                         int _n, int _m, const std::vector<ZoneControl>& c,
                         const std::vector<int>& _init_regions,
                         long long seed)
: p(_p), w(_w), data(_data), dist_matrix(_dist_matrix), n(_n), m(_m), controls(c),
am(_n, _m, _w, _data, _dist_matrix),
rng(seed), is_control_satisfied(true), init_regions(_init_regions), objInfo(-1)
{
    if (p < 0) {
        is_control_satisfied = false;
        return;
    }

    // init unassigned areas
    for (int i=0; i<n; ++i) {
        unassignedAreas[i] = true;
    }

    // mark neighborless areas
    AssignAreasNoNeighs();

    // get p start areas using k-means
    if (init_regions.empty()) {
        std::vector<int> seeds = this->kmeansInit();

        // for replicate clusterpy: CA_polygons
        // calif.cluster('azp', ['PCR2002'], 9, wType='rook', dissolve=1)
        // [31, 38, 37, 28, 20, 29, 17, 50, 55]
        //seeds[0] = 31; seeds[1] = 38; seeds[2] = 37; seeds[3]=28; seeds[4]=20;
        //seeds[5] = 29; seeds[6] = 17; seeds[7] = 50; seeds[8] = 55;

        // process and assign the neighbors of p starting areas
        this->setSeeds(seeds);

        // for any other unassigned areas, assign to a region
        while (unassignedAreas.size() != 0) {
            this->constructRegions();
        }

        //  create objectiveFunction object for local improvement
        objective_function = new ObjectiveFunction(n, m, data, w, region2Area);

        // get objective function value
        this->objInfo = objective_function->GetValue();

    } else {
        // // get p start areas using init_regions from input
        InitFromRegion(init_regions);
    }
}


RegionMaker::~RegionMaker()
{
    if (objective_function) {
        delete objective_function;
    }
}

void RegionMaker::Copy(RegionMaker& rm)
{
    this->init_regions = rm.init_regions;
    this->unassignedAreas = rm.unassignedAreas;
    this->assignedAreas = rm.assignedAreas;
    this->areaNoNeighbor = rm.areaNoNeighbor;
    this->area2Region = rm.area2Region;
    this->region2Area = rm.region2Area;
    this->potentialRegions4Area = rm.potentialRegions4Area;
    this->candidateInfo = rm.candidateInfo;
    this->objInfo = rm.objInfo;
    if (objective_function) {
        delete objective_function;
    }
    this->objective_function = new ObjectiveFunction(n, m, data, w, region2Area);
}

void RegionMaker::InitFromRegion(std::vector<int>& init_regions)
{
    // check init regions, index start from 1
    for (size_t i=0; i<init_regions.size(); ++i) {
        if (init_regions[i] > 0 && init_regions[i] <= p) {
            int areaID = i;
            int regionID = init_regions[i] - 1;
            assignAreaStep1(areaID, regionID);
        }
    }
    // for any unassigned areas, create potentialRegions4Area
    if (unassignedAreas.size() != 0) {
        for (int r=0; r<p; ++r) {
            std::set<int> buffer_areas = getBufferingAreas(region2Area[r]);
            std::set<int>::iterator it;
            for (it = buffer_areas.begin(); it != buffer_areas.end(); ++it) {
                int areaID = *it;
                if (unassignedAreas.find(areaID) != unassignedAreas.end()) {
                    potentialRegions4Area[areaID].insert(r);
                }
            }
        }
    }

    // for any other unassigned areas, assign to a region
    while (unassignedAreas.size() != 0) {
        this->constructRegions();
    }

    //  create objectiveFunction object for local improvement
    objective_function = new ObjectiveFunction(n, m, data, w, region2Area);

    // get objective function value
    this->objInfo = objective_function->GetValue();
}

void RegionMaker::AssignAreasNoNeighs()
{
    // w should not be NULL
    for (int i=0; i<n; ++i) {
        if (w[i].Size() == 0) {
            areaNoNeighbor[i] = true;
            assignedAreas[i] = true;
        } else {
            areaNoNeighbor[i] = false;
            unassignedAreas[i] = true;
        }
    }
}

std::vector<int> RegionMaker::kmeansInit()
{
    //Initial p starting observations using K-Means
    std::vector<double> probabilities(n, 0);
    std::vector<double> distances(n, 1);
    double total = n;
    for (int j=0; j<n; ++j) {
        probabilities[j] = distances[j] / total;
    }

    std::vector<int> seeds;
    for (int i=0; i < p; ++i) {
        double random = rng.nextDouble();
        //std::cout << random << std::endl;
        bool find = false;
        double acum = 0;
        int cont = 0;
        while (find == false) {
            double inf = acum;
            double sup = acum + probabilities[cont];
            if (inf <= random && random <= sup) {
                find = true;
                seeds.push_back(cont);
                // for each observation/area, find which seed is the closest
                total = 0;
                for (int j=0; j<n; ++j) {
                    double distancei = 0;
                    for (size_t k=0; k<seeds.size(); ++k) {
                        double d = dist_matrix->getDistance(j, seeds[k]);
                        if (k==0) distancei = d;
                        else if (d <  distancei) {
                            distancei = d;
                        }
                    }
                    total += distancei; // get a sum of new distances
                    distances[j] = distancei;
                }
                for (int j=0; j<n; ++j) {
                    probabilities[j] = distances[j] / total;
                }
            } else {
                cont += 1;
                acum = sup;
            }
        }
    }
    return seeds;
}

bool RegionMaker::IsSatisfyControls()
{
    REGION_AREAS::iterator it;
    for (it = region2Area.begin(); it != region2Area.end(); ++it) {
        for (size_t i=0;  i<controls.size(); ++i) {
            if (controls[i].CheckBound(it->second)  == false) {
                return false;
            }
        }
    }
    return true;
}

void RegionMaker::setSeeds(std::vector<int> seeds)
{
    // Sets the initial seeds for clustering
    if ((int)seeds.size() < p) {
        // pick more to fill up seeds
        std::vector<bool> didx(n, true);
        // remove seeds, and neighborless observations
        for (size_t i=0; i<seeds.size(); ++i) {
            didx[ seeds[i] ] = false;
        }
        for (int i=0; i<n; ++i) {
            if (w[i].Size() == 0) {
                didx[i] = false;
            }
        }
        std::vector<int> cands;
        for (size_t i=0; i<didx.size(); ++i) {
            if (didx[i] == true) {
                cands.push_back(i);
            }
        }
        DataUtils::Shuffle(cands, rng);

        for (size_t i=0; i< (p -seeds.size()); ++i) {
            seeds.push_back(cands[i]);
        }
    }

    // assign each seed with a new region
    int c = 0;
    for (size_t i=0; i< seeds.size(); ++i) {
        assignAreaStep1(seeds[i], c);
        c += 1;
    }

    // grow the region if needed
    if (controls.size() > 0) {
       bool grow_flag = growRegion();

       if (grow_flag == false) {
           // raise exception
           is_control_satisfied = false;
       }
    }

    // find potential
    for (int i=0; i<p; ++i) {
        // check neighbors of areaID that are not been assigned yet
        // and assign neighbor to potential regions
        std::set<int> buffer_areas = getBufferingAreas(region2Area[i]);
        std::set<int>::iterator it;

        for (it = buffer_areas.begin(); it != buffer_areas.end(); ++it) {
            int neigh = *it;
            if (assignedAreas.find(neigh) == assignedAreas.end()) {
                potentialRegions4Area[neigh].insert(i);
            }
        }
    }
}

bool RegionMaker::growRegion()
{
    bool is_valid = false;

    std::map<int, bool> grow_flags;
    for (int i=0; i< p; ++i)  grow_flags[i] = true;

    std::set<int>::iterator it;
    while (is_valid == false) {
        is_valid = true;
        for (int i=0; i< p; ++i) {
            // check if grow is needed
            if (grow_flags[i] == false) {
                continue;
            }
            // still need to grow
            is_valid = false;

            // each time, grow just one area, to avoid dominant grow
            // if two seeds are next to each other
            bool has_assign = false;
            std::set<int> buffer_areas = getBufferingAreas(region2Area[i]);
            for (it = buffer_areas.begin(); !has_assign && it != buffer_areas.end(); ++it){
                int nn = *it;
                if (assignedAreas.find(nn) == assignedAreas.end()) {
                    assignAreaStep1(nn, i);
                    has_assign = true;
                }
            }
            // check if this region satisfy the low bounds
            bool satisfy = true;
            for (size_t j=0; satisfy && j<controls.size(); ++j) {
                if (controls[j].SatisfyLowerBound(region2Area[i])  == false) {
                    satisfy = false;
                }
            }
            // check infinite look
            if (satisfy == false && has_assign == false) {
                // can't grow region to satisfy controls
                return false;
            }
            // update grow flag
            grow_flags[i] = !satisfy;
        }
    }


    return true;
}

void RegionMaker::assignSeeds(int areaID, int regionID)
{
    // Assign an area to a region and updates potential regions for the neighs
    // parameters
    assignAreaStep1(areaID, regionID);

    // check neighbors of areaID that are not been assigned yet
    // and assign neighbor to potential regions
    const std::vector<long>& nbrs = w[areaID].GetNbrs();
    for (size_t i=0; i<nbrs.size(); ++i) {
        int neigh = (int)nbrs[i];
        if (assignedAreas.find(neigh) == assignedAreas.end()) {
            potentialRegions4Area[neigh].insert(regionID);
        }
    }
}

void RegionMaker::assignAreaStep1(int areaID, int regionID)
{
    //  Assgin an area to a region
    region2Area[regionID][areaID] = false;

    // attach region with area
    area2Region[areaID] = regionID;

    // remove areaID from unassignedAreas
    unassignedAreas.erase(areaID);

    // add areaID to assignedAreas
    assignedAreas[areaID] = true;
}

void RegionMaker::constructRegions()
{
    // Construct potential regions per area (not assigned)
    std::map<int, std::set<int> >::iterator it;
    std::set<int>::iterator rit, ait;

    // Process the most feasible area that has shortest distance to a region

    for (it=potentialRegions4Area.begin(); it != potentialRegions4Area.end(); ++it) {
        int areaID = it->first;
        const std::set<int>& regionIDs = it->second;

        for (rit = regionIDs.begin(); rit != regionIDs.end(); ++rit) {
            int region = *rit;
            std::pair<int, int> a_r = std::make_pair(areaID, region);
            double regionDistance = am.getDistance2Region(areaID, region, region2Area);
            // (areaID, region): distance
            candidateInfo[a_r] = regionDistance;
        }
    }

    if (candidateInfo.empty() == false)  {

        // Select and assign the nearest area to a region
        // minimumSelection(RegionMaker)
        // if there are more than one (area, region) pair
        // random select from cands

        std::vector<std::pair<int, int> > cands;
        std::map<std::pair<int, int>, double>::iterator cit;

        // get min dist
        double min_region_distance = std::numeric_limits<double>::max();
        for (cit = candidateInfo.begin(); cit != candidateInfo.end(); ++cit) {
            if (cit->second < min_region_distance)  {
                min_region_distance = cit->second;
            }
        }
        // get all pairs with min dist
        for (cit = candidateInfo.begin(); cit != candidateInfo.end(); ++cit) {
            if (cit->second == min_region_distance) {
                cands.push_back(cit->first);
            }
        }

        // if random select candidate from pairs is needed
        int rnd_sel = cands.size() == 1? 0 : rng.nextInt((int)cands.size());
        std::pair<int, int>& sel_area_region = cands[rnd_sel];
        int aid = sel_area_region.first;
        int rid = sel_area_region.second;

        // assign select areaID to regionID, and process the neighbors of areaID
        // update the centroid of the regionID
        if (this->assignArea(aid, rid)) {

            // remove from candidateInfo with areaID=aid
            std::vector<std::pair<int, int> > removed_cands;
            for (cit = candidateInfo.begin(); cit != candidateInfo.end(); ++cit) {
                if (cit->first.first == aid) {
                    removed_cands.push_back(cit->first);
                }
            }
            for (size_t i=0; i<removed_cands.size(); ++i) {
                candidateInfo.erase(removed_cands[i]);
            }
        } else {
            // can not assign, violate controls
            // just remove current pair
            candidateInfo.erase(sel_area_region);
        }
    }
}



bool RegionMaker::assignArea(int areaID, int regionID)
{
    // Check upper bounds of controls only, since we are assigning
    for (size_t i=0; i<controls.size(); ++i) {
        if (controls[i].CheckAdd(areaID, region2Area[regionID])  == false) {
            return false;
        }
    }

    //Assign an area to a region and updates potential regions for its neighs
    this->assignAreaStep1(areaID, regionID);

    //for neigh in self.neighsMinusAssigned:
    // assign regionID as a potential region for neigh
    const std::vector<long>& neighs = this->w[areaID].GetNbrs();
    for (size_t i=0; i< neighs.size(); ++i) {
        int nn = (int)neighs[i];
        if (assignedAreas.find(nn) == assignedAreas.end()) {
            // for not yet assigned neighbor
            potentialRegions4Area[nn].insert(regionID);
        }
    }
    // remove areaID since its processed
    potentialRegions4Area.erase(areaID);

    // update centroid of the region
    am.updateRegionCentroids(regionID, region2Area);

    return true;
}

std::vector<int> RegionMaker::returnRegions()
{
    std::vector<int> areasId, results;
    boost::unordered_map<int, int>::iterator it;
    for (it = area2Region.begin(); it != area2Region.end(); ++it) {
        areasId.push_back( it->first );
    }
    std::sort(areasId.begin(), areasId.end());
    for (size_t i=0; i<areasId.size(); ++i) {
        results.push_back( area2Region[areasId[i] ]);
    }
    return results;
}

std::set<int> RegionMaker::getBufferingAreas(boost::unordered_map<int, bool>& areas)
{
    std::set<int> buffering_areas;

    //boost::unordered_map<int, bool>& areas  = region2Area[regionID];
    boost::unordered_map<int, bool>::iterator it;
    for (it = areas.begin(); it != areas.end(); ++it) {
        int area = it->first;
        const std::vector<long>& nn  = w[area].GetNbrs();
        for (size_t i=0; i<nn.size(); ++i) {
            if (areas.find((int)nn[i])  == areas.end()) {
                // neighbor not in this region, then the area is at border
                buffering_areas.insert((int)nn[i]);
            }
        }
    }

    return buffering_areas;
}

void RegionMaker::getBorderingAreas(int regionID)
{
    // check every area, and return who has neighbors out of the region
    boost::unordered_map<int, bool>& areas  = region2Area[regionID];
    boost::unordered_map<int, bool>::iterator it;
    for (it = areas.begin(); it != areas.end(); ++it) {
        int area = it->first;
        const std::vector<long>& nn  = w[area].GetNbrs();
        areas[area] = false;  // not a border area by default
        for (size_t i=0; i<nn.size(); ++i) {
            if (areas.find((int)nn[i])  == areas.end()) {
                // neighbor not in this region, the the area is at border
                areas[area] = true;
                break;
            }
        }
    }
}

std::set<int> RegionMaker::getPossibleMove(int area)
{
    std::set<int> moves;
    int myRegion = area2Region[area];

    // Check if moving area out of myRegion satisfy controls
    for (size_t i=0; i<controls.size(); ++i) {
        if (controls[i].CheckRemove(area, region2Area[myRegion]) == false) {
            return moves;
        }
    }

    const std::vector<long>& nn  = w[area].GetNbrs();
    for (size_t i=0; i<nn.size(); ++i) {
        // check neighbors, which region it belongs to
        int nbrRegion = area2Region[(int)nn[i]];
        if (nbrRegion != myRegion) {
            // check if moving area to nbrRegion satisfy controls
            for (size_t i=0; i<controls.size(); ++i) {
                if (controls[i].CheckAdd(area, region2Area[nbrRegion]) == false) {
                    return std::set<int>();
                }
            }
            moves.insert(nbrRegion);
        }
    }
    return moves;
}
////////////////////////////////////////////////////////////////////////////////
////// AZP
////////////////////////////////////////////////////////////////////////////////
void AZP::LocalImproving()
{
    // NOTE: this could be parallelized: when one thread update region A
    // and A exchanges area with another region B, only the threads interact
    // with region B should wait and other threads can run parallely to update
    // regions C,D,E... However, the mutex lock should be carefully designed...

    int improve = 1;
    //std::vector<int> rand_test = {7,4,3,5,3,3,1,1,0,3,3,2,3,4,2,0,0,0};
    //std::vector<int> rand_test1 = {57, 56, 52, 24, 16, 10, 3, 24, 51, 57, 22, 57, 46, 11, 46, 52, 50, 46, 10, 52, 24, 57, 3, 51, 7, 5, 20, 30, 28, 8, 26, 39, 43, 18, 55, 41, 36, 29, 17, 0, 56, 33, 35, 1, 23, 9, 32, 22, 2, 49, 15, 11, 48, 14, 16, 50, 34, 12, 42, 40, 31, 45, 44, 31, 30, 28, 8, 20, 40, 42, 17, 41, 18, 26, 55, 43, 39, 29, 36, 44, 31, 14, 11, 16, 2, 48, 0, 1, 15, 35, 50, 12, 23, 9, 49, 33, 32, 34, 56, 22, 24, 7, 45, 57, 10, 51, 5, 3, 46, 52};
    //int rr = 0, rr1 = 0;
    while (improve == 1) {
        std::vector<int> regions;
        for (int i=0; i<p; ++i) {
            regions.push_back(i);
        }

        while (regions.size() > 0) {
            // step 3
            int randomRegion = 0;
            if (regions.size() > 1) {
                randomRegion = rng.nextInt((int)regions.size());
            }
            //randomRegion = rand_test[rr++];
            int region = regions[randomRegion];
            regions.erase(std::find(regions.begin(), regions.end(), region));

            // step 4
            // get bordering areas in region
            getBorderingAreas(region);

            boost::unordered_map<int, bool>& areas = region2Area[region];
            boost::unordered_map<int, bool>::iterator area_it;
            std::set<int>::iterator move_it;
            improve = 0;

            while (areas.size() > 1) {
                // step 5
                bool nothing_can_do = true;
                bool moved = false;
                for (area_it = areas.begin(); !moved && area_it != areas.end(); ++area_it) {
                    if (area_it->second == true) { // only procee the bordering area
                        nothing_can_do = false;
                        // pick a random bordering area, mark it as processed
                        int randomArea = area_it->first;
                        areas[randomArea] = false;
                        // get possible move of this randomArea
                        std::set<int> possibleMove = getPossibleMove(randomArea);
                        // check obj change before contiguity check
                        for (move_it = possibleMove.begin();
                             !moved && move_it != possibleMove.end();
                             ++move_it)
                        {
                            int move = *move_it;
                            // request lock for region[move], pause threads[move]
                            //std::pair<double, bool> swap = objective_function->TrySwap(randomArea, region, move);
                            //double delta = swap.first;
                            //bool swapped = swap.second;
                            //if (swapped) {
                            std::pair<double, bool> swap = objective_function->TrySwapSA(randomArea, region, move, this->objInfo);
                            double obj = swap.first;
                            bool contiguous = swap.second;
                            if (obj <= this->objInfo && contiguous) {
                                improve = 1;
                                moved = true;

                                this->area2Region[randomArea] = move;
                                this->objInfo = obj;

                                // move happens, update bordering area
                                getBorderingAreas(region);
                                //getBorderingAreas(move);

                                //std::cout << "--- Local improvement (area, region)" << randomArea << "," << move << std::endl;
                                //std::cout << "--- New Objective Function value: " << this->objInfo << std::endl;

                            }
                            // release lock for region[move], awake threads[move]
                        }
                    }
                    if (moved) {
                        break;
                    }
                }
                if (nothing_can_do) {
                    break;
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
////// MaxpRegionMaker
////////////////////////////////////////////////////////////////////////////////
MaxpRegionMaker::MaxpRegionMaker(GalElement* const _w,
                                 double** _data, // row-wise
                                 RawDistMatrix* _dist_matrix,
                                 int _n, int _m, const std::vector<ZoneControl>& c,
                                 const std::vector<int>& _init_areas,
                                 long long seed)
: RegionMaker(-1, _w, _data, _dist_matrix, _n, _m, c, std::vector<int>(), seed),
init_areas(_init_areas)
{
    objective_function = 0;
    p = 0;
    InitSolution();
}

void MaxpRegionMaker::InitSolution()
{
    // init unassigned areas
    for (int i=0; i<n; ++i) {
        unassignedAreas[i] = true;
    }

    // mark neighborless areas
    AssignAreasNoNeighs();

    // get max-p regions
    // get starting position: either from init_areas or randomly selected
    std::set<int> start_areas;
    for (size_t i=0; i<init_areas.size(); ++i) {
        start_areas.insert(init_areas[i]);
    }
    std::vector<int> _candidates;
    for (int i=0; i < n; i++) {
        if (start_areas.empty() || start_areas.find(i) == start_areas.end()) {
            _candidates.push_back(i);
        }
    }
    DataUtils::Shuffle(_candidates, rng);
    for (size_t i=0; i<init_areas.size(); ++i) {
        _candidates.insert(_candidates.begin(), init_areas[i]);
    }
    std::list<int> candidates;
    boost::unordered_map<int, bool> candidates_dict;
    for (int i=0; i<n;i++) {
        candidates.push_back( _candidates[i] );
        candidates_dict[ _candidates[i] ] = false;
    }

    // grow p regions using the starting positions above
    int r = 0;
    bool satisfy_lower_bound = false;

    while (!candidates.empty()) {
        int seed = candidates.front();
        candidates.pop_front();

        // try to grow it till threshold constraint is satisfied
        bool is_growing = true;
        bool reach_ub = false;
        bool reach_lb = false;
        std::set<int>::iterator it;

        // assign this seed with a new region
        boost::unordered_map<int, bool> region;
        region[seed] = true;
        candidates_dict[seed] = true;

        while (is_growing && !reach_ub && !reach_lb) {
            // each time, grow just one area, to avoid dominant grow
            // if two seeds are next to each other
            bool has_assign = false;
            std::set<int> buffer_areas = getBufferingAreas(region);
            for (it = buffer_areas.begin(); !has_assign && it != buffer_areas.end(); ++it) {
                int nn = *it;
                if (candidates_dict[nn] == false) { // not processed
                    // check upper bound if adding this area to region
                    for (size_t j=0; !reach_ub && j<controls.size(); ++j) {
                        if (controls[j].CheckAdd(nn, region)  == false) {
                            // reach upper bound, no need to grow anymore
                            reach_ub = true;
                            is_growing = false;
                        }
                    }
                    if (!reach_ub) {
                        region[nn] = true;
                        candidates.remove(nn);
                        candidates_dict[nn] = true; // processed
                        has_assign = true;
                    }
                }
            }
            // check if this region satisfy the low bounds
            bool satisfy = true;
            for (size_t j=0; satisfy && j<controls.size(); ++j) {
                if (controls[j].SatisfyLowerBound(region)  == false) {
                    satisfy = false;
                }
            }
            // check infinite loop
            if (satisfy == false && has_assign == false) {
                // can't grow region to satisfy controls
                is_growing = false;
            }
            reach_lb = satisfy;
        }

        // the region will be valid only if reaches the lower bound
        if (reach_lb) {
            boost::unordered_map<int, bool>::iterator rit;
            for (rit = region.begin(); rit != region.end(); ++rit) {
                assignAreaStep1(rit->first, r);
            }
            r = r + 1; // create another region
            satisfy_lower_bound = true;
        }
    }
    if (satisfy_lower_bound) {
        // enclave assignment
        // find potential
        if (unassignedAreas.size() > 0) {
            for (size_t i=0; i<region2Area.size(); ++i) {
                // check neighbors of areaID that are not been assigned yet
                // and assign neighbor to potential regions
                std::set<int> buffer_areas = getBufferingAreas(region2Area[i]);
                std::set<int>::iterator it;

                for (it = buffer_areas.begin(); it != buffer_areas.end(); ++it) {
                    int neigh = *it;
                    if (assignedAreas.find(neigh) == assignedAreas.end()) {
                        potentialRegions4Area[neigh].insert(i);
                    }
                }
            }
        }

        // for any other unassigned areas (enclaves), assign to a region
        while (unassignedAreas.size() != 0) {
            this->constructRegions();
        }

        // now we have p regions
        p = (int)region2Area.size();

        //  create objectiveFunction object for local improvement
        objective_function = new ObjectiveFunction(n, m, data, w, region2Area);

        // get objective function value
        this->objInfo = objective_function->GetValue();
    }
}

std::vector<int> MaxpRegionMaker::returnRegions()
{
    // area2Region is not a complete solution, others that are not assigned will
    // have region 0
    std::vector<int> areasId, results(n, 0);
    boost::unordered_map<int, int>::iterator it;
    for (it = area2Region.begin(); it != area2Region.end(); ++it) {
        areasId.push_back( it->first );
        results[ it->first ] = it->second + 1;
    }
    return results;
}

////////////////////////////////////////////////////////////////////////////////
////// MaxpRegion
////////////////////////////////////////////////////////////////////////////////
MaxpRegion::MaxpRegion(int _max_iter, GalElement* const _w,
                       double** _data, // row-wise
                       RawDistMatrix* _dist_matrix,
                       int _n, int _m, const std::vector<ZoneControl>& c,
                       const std::vector<int>& init_regions,
                       long long _seed, int cpu_threads)
: RegionMaker(-1, _w, _data, _dist_matrix, _n, _m, c, std::vector<int>(), _seed),
seed(_seed), init_areas(init_regions), max_iter(_max_iter), cpu_threads(cpu_threads) {
    objective_function = 0;
    largest_p = 0;
    best_of = DBL_MAX;
#ifndef __NO_THREAD__
#ifdef __USE_PTHREAD__
    pthread_mutex_init(&lock, 0);
    pthread_cond_init(&wcond, 0);
#endif
#endif
}

MaxpRegion::~MaxpRegion()
{
#ifndef __NO_THREAD__
#ifdef __USE_PTHREAD__
    pthread_mutex_destroy(&lock);
    pthread_cond_destroy(&wcond);
#endif
#endif
}

void MaxpRegion::Run() {
#ifdef __NO_THREAD__
    // construction phase: find a collection of feasible solution with largest p
    for (int iter=0; iter< max_iter; ++iter) {
        RunConstruction(seed+iter);
    }
    // local improvement phase:
    int i=0;
    std::map<double, std::vector<int> >::iterator it;
    for (it = candidates.begin(); it != candidates.end(); ++it) {
        initial_objectivefunction = it->first;
        std::vector<int> solution = it->second;
        RunAZP(solution, seed+initial_objectivefunction, i++) ;
    }
#else
    PhaseConstructionThreaded();

    std::map<double, std::vector<int> >::iterator it;
    candidate_keys.clear();
    for (it = candidates.begin(); it != candidates.end(); ++it) {
        candidate_keys.push_back(it->first);
    }

    PhaseLocalImprovementThreaded();
#endif

    final_objectivefunction = best_of;
    final_solution = best_result;
}

void MaxpRegion::PhaseConstructionThreaded()
{
    int nCPUs = cpu_threads;
    if (nCPUs <=0) nCPUs = 1;
    int work_chunk = max_iter / nCPUs;
    if (work_chunk == 0) work_chunk = 1;
    int quotient = max_iter / nCPUs;
    int remainder = max_iter % nCPUs;
    int tot_threads = (quotient > 0) ? nCPUs : remainder;

#ifdef __NO_THREAD__
    RunConstructionRange(0, max_iter-1);
#else
#ifdef __USE_PTHREAD__
    pthread_t *threadPool = new pthread_t[nCPUs];
    struct maxp_thread_args *args = new maxp_thread_args[nCPUs];
    for (int i=0; i<tot_threads; i++) {
        int a=0, b=0;
        if (i < remainder) {
            a = i*(quotient+1);
            b = a+quotient;
        } else {
            a = remainder*(quotient+1) + (i-remainder)*quotient;
            b = a+quotient-1;
        }
        args[i].maxp = this;
        args[i].start = a;
        args[i].end = b;
        if (pthread_create(&threadPool[i], NULL, &maxp_thread_helper_construction, &args[i])) {
            perror("Thread create failed.");
        }
    }
    for (int j = 0; j < nCPUs; j++) {
        pthread_join(threadPool[j], NULL);
    }
    delete[] args;
    delete[] threadPool;
#else
    boost::thread_group threadPool;
    for (int i=0; i<tot_threads; i++) {
        int a=0, b=0;
        if (i < remainder) {
            a = i*(quotient+1);
            b = a+quotient;
        } else {
            a = remainder*(quotient+1) + (i-remainder)*quotient;
            b = a+quotient-1;
        }
        boost::thread* worker = new boost::thread(boost::bind(&MaxpRegion::RunConstructionRange, this, a,b));
        threadPool.add_thread(worker);
    }
    threadPool.join_all();
#endif
#endif
}

void MaxpRegion::PhaseLocalImprovementThreaded()
{
    int nCPUs = cpu_threads;
    if (nCPUs <=0) nCPUs = 1;
    int nCandidates = candidates.size();
    int work_chunk = nCandidates / nCPUs;
    if (work_chunk == 0) work_chunk = 1;
    int quotient = nCandidates / nCPUs;
    int remainder = nCandidates % nCPUs;
    int tot_threads = (quotient > 0) ? nCPUs : remainder;
#ifdef __NO_THREAD__
    RunLocalImprovementRange(0, nCandidates-1);
#else
#ifdef __USE_PTHREAD__
    pthread_t *threadPool = new pthread_t[nCPUs];
    struct maxp_thread_args *args = new maxp_thread_args[nCPUs];
    for (int i=0; i<tot_threads; i++) {
        int a=0, b=0;
        if (i < remainder) {
            a = i*(quotient+1);
            b = a+quotient;
        } else {
            a = remainder*(quotient+1) + (i-remainder)*quotient;
            b = a+quotient-1;
        }
        args[i].maxp = this;
        args[i].start = a;
        args[i].end = b;
        if (pthread_create(&threadPool[i], NULL, &maxp_thread_helper_localimprovement, &args[i])) {
            perror("Thread create failed.");
        }
    }
    for (int j = 0; j < nCPUs; j++) {
        pthread_join(threadPool[j], NULL);
    }
    delete[] args;
    delete[] threadPool;
#else
    boost::thread_group threadPool;
    for (int i=0; i<tot_threads; i++) {
        int a=0, b=0;
        if (i < remainder) {
            a = i*(quotient+1);
            b = a+quotient;
        } else {
            a = remainder*(quotient+1) + (i-remainder)*quotient;
            b = a+quotient-1;
        }
        boost::thread* worker = new boost::thread(boost::bind(&MaxpRegion::RunLocalImprovementRange, this, a, b));
        threadPool.add_thread(worker);
    }
    threadPool.join_all();
#endif
#endif
}

void MaxpRegion::RunConstructionRange(int start, int end)
{
    for (int i=start; i <= end; ++i) {
        RunConstruction(seed + i);
    }
}

void MaxpRegion::RunLocalImprovementRange(int start, int end)
{
    for (int i=start; i<=end; ++i) {
        double initial_objectivefunction = candidate_keys[i];
        std::vector<int> solution = candidates[initial_objectivefunction];
        RunAZP(solution, seed+initial_objectivefunction, i);
    }
}

void MaxpRegion::RunConstruction(long long seed)
{
    MaxpRegionMaker rm_local(w, data, dist_matrix, n, m, controls, init_areas, seed);
    int tmp_p = rm_local.GetPRegions();
    double of = rm_local.GetInitObjectiveFunction();
#ifndef __NO_THREAD__
#ifndef __USE_PTHREAD__
    mutex.lock();
#else
    pthread_mutex_lock(&lock);
#endif
#endif
    if (largest_p < tmp_p) {
        candidates.clear(); // new collection for largest p
        largest_p = tmp_p;
    }
    if (largest_p == tmp_p) {
        candidates[of] = rm_local.GetResults(); // could have duplicates
    }
#ifndef __NO_THREAD__
#ifndef __USE_PTHREAD__
    mutex.unlock();
#else
    //pthread_cond_signal(&wcond);
    pthread_mutex_unlock(&lock);
#endif
#endif
}

MaxpGreedy::MaxpGreedy(int _max_iter, GalElement* const _w,
               double** _data, // row-wise
               RawDistMatrix* _dist_matrix,
               int _n, int _m, const std::vector<ZoneControl>& c,
               const std::vector<int>& init_regions,
               long long seed, int cpu_threads)
        : MaxpRegion(_max_iter, _w, _data, _dist_matrix, _n, _m, c, init_regions, seed, cpu_threads)
{
    Run();
}

void MaxpGreedy::RunAZP(std::vector<int>& solution, long long seed, int i)
{
    AZP azp(largest_p, w, data, dist_matrix, n, m, controls, 0, solution, seed);

    std::vector<int> result = azp.GetResults();
    double of = azp.GetFinalObjectiveFunction();
#ifndef __NO_THREAD__
#ifndef __USE_PTHREAD__
    mutex.lock();
#else
    pthread_mutex_lock(&lock);
#endif
#endif

    if (of < best_of) {
        best_result = result;
        best_of = of;
    }

#ifndef __NO_THREAD__
#ifndef __USE_PTHREAD__
    mutex.unlock();
#else
    pthread_mutex_unlock(&lock);
#endif
#endif
}

MaxpSA::MaxpSA(int _max_iter, GalElement* const _w,
               double** _data, // row-wise
               RawDistMatrix* _dist_matrix,
               int _n, int _m, const std::vector<ZoneControl>& c,
               double _alpha, int _sa_iter,
               const std::vector<int>& init_regions,
               long long seed, int cpu_threads)
: MaxpRegion(_max_iter, _w, _data, _dist_matrix, _n, _m, c, init_regions, seed, cpu_threads),
temperature(1.0), alpha(_alpha), sa_iter(_sa_iter)
{
    Run();
}

void MaxpSA::RunAZP(std::vector<int>& solution, long long seed, int i)
{
    AZPSA azp(largest_p, w, data, dist_matrix, n, m, controls, alpha, sa_iter, 0, solution, seed);

    std::vector<int> result = azp.GetResults();
    double of = azp.GetFinalObjectiveFunction();

#ifndef __NO_THREAD__
#ifndef __USE_PTHREAD__
    mutex.lock();
#else
    pthread_mutex_lock(&lock);
#endif
#endif

    if (of < best_of) {
        best_result = result;
        best_of = of;
    }

#ifndef __NO_THREAD__
#ifndef __USE_PTHREAD__
    mutex.unlock();
#else
    pthread_mutex_unlock(&lock);
#endif
#endif
}

MaxpTabu::MaxpTabu(int _max_iter, GalElement* const _w,
               double** _data, // row-wise
               RawDistMatrix* _dist_matrix,
               int _n, int _m, const std::vector<ZoneControl>& c,
               int _tabu_length, int _conv_tabu,
               const std::vector<int>& init_regions,
               long long seed, int cpu_threads)
: MaxpRegion(_max_iter, _w, _data, _dist_matrix, _n, _m, c, init_regions, seed, cpu_threads),
tabuLength(_tabu_length), convTabu(_conv_tabu)
{
    Run();
}

void MaxpTabu::RunAZP(std::vector<int>& solution, long long seed, int i)
{
    if (convTabu == 0) {
        convTabu = std::max(10, n/largest_p); //230 * sqrt(largest_p);
    }

    AZPTabu azp(largest_p, w, data, dist_matrix, n, m, controls, tabuLength, convTabu, 0, solution, seed);

    std::vector<int> result = azp.GetResults();
    double of = azp.GetFinalObjectiveFunction();

#ifndef __NO_THREAD__
#ifndef __USE_PTHREAD__
    mutex.lock();
#else
    pthread_mutex_lock(&lock);
#endif
#endif

    if (of < best_of) {
        best_result = result;
        best_of = of;
    }


#ifndef __NO_THREAD__
#ifndef __USE_PTHREAD__
    mutex.unlock();
#else
    pthread_mutex_unlock(&lock);
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
////// AZP Simulated Annealing
////////////////////////////////////////////////////////////////////////////////
void AZPSA::LocalImproving()
{
    // Openshaw's Simulated Annealing for AZP algorithm
    int totalMoves = 0;
    int acceptedMoves = 0;
    double bestOBJ = this->objInfo;
    double currentOBJ = this->objInfo;
    std::vector<int> bestRegions = this->returnRegions();
    std::vector<int> currentRegions = this->returnRegions();
    REGION_AREAS region2AreaBest = this->region2Area;
    boost::unordered_map<int, int> area2RegionBest = this->area2Region;
    std::set<int>::iterator it;

    int improve = 1;
    while (improve == 1) {
        std::vector<int> regions(p);
        for (int i=0; i<p; ++i) regions[i] = i;

        while (regions.size() > 0) {
            // step 3
            int randomRegion = 0;
            if (regions.size() > 1) {
                randomRegion = rng.nextInt((int)regions.size());
            }
            //randomRegion = rand_test[rr++];
            int region = regions[randomRegion];
            regions.erase(std::find(regions.begin(), regions.end(), region));

            // step 4
            // get bordering areas in region
            getBorderingAreas(region);

            boost::unordered_map<int, bool>& areas = region2Area[region];
            boost::unordered_map<int, bool>::iterator area_it;
            std::set<int>::iterator move_it;
            improve = 0;

            while (areas.size() > 1) {
                // step 5
                bool nothing_can_do = true;
                bool moved = false;
                for (area_it = areas.begin(); !moved && area_it != areas.end(); ++area_it) {
                    if (area_it->second == true) { // only procee the bordering area
                        nothing_can_do = false;
                        // pick a random bordering area, mark it as processed
                        int randomArea = area_it->first;
                        areas[randomArea] = false;
                        // get possible move of this randomArea
                        std::set<int> possibleMove = getPossibleMove(randomArea);
                        // check obj change before contiguity check
                        for (move_it = possibleMove.begin();
                             !moved && move_it != possibleMove.end();
                             ++move_it)
                        {
                            int move = *move_it;
                            // request lock for region[move], pause threads[move]
                            std::pair<double, bool> swap = objective_function->TrySwapSA(randomArea, region, move, bestOBJ);
                            double obj = swap.first;
                            bool contiguous = swap.second;
                            if (obj <= bestOBJ && contiguous) { // means swapped
                                improve = 1;
                                moved = true;
                                this->area2Region[randomArea] = move;
                                this->objInfo = obj;

                                // update SA variables
                                bestOBJ= obj;
                                currentOBJ = obj;
                                bestRegions = this->returnRegions();
                                region2AreaBest = this->region2Area;
                                area2RegionBest = this->area2Region;

                                // move happens, update bordering area
                                getBorderingAreas(region);
                                //getBorderingAreas(move);

                                //std::cout << "--- Local improvement (area, region)" << randomArea << "," << move << std::endl;
                                //std::cout << "--- New Objective Function value: " << this->objInfo << std::endl;
                                // step 4
                            }
                            else if (!(obj <= bestOBJ && contiguous==false)) { // not consider NON-contiguous
                                contiguous = objective_function->checkFeasibility(region, randomArea);
                                if (contiguous) {
                                    double random = rng.nextDouble();
                                    totalMoves += 1;
                                    double sa = std::exp(-(obj - currentOBJ) * n / (currentOBJ * temperature));
                                    if (sa > random && region2Area[region].size() > 1) { // make move if SA satisfies
                                        double new_obj = objective_function->MakeMove(randomArea, region, move);
                                        if (new_obj > 0) { // prevent empty region
                                            moved = true;
                                            this->area2Region[randomArea] = move;
                                            this->objInfo = new_obj;

                                            // update SA variables
                                            acceptedMoves += 1;
                                            currentOBJ = new_obj;
                                            if (new_obj < bestOBJ) {
                                                bestOBJ = new_obj;
                                                currentOBJ = new_obj;
                                                bestRegions = this->returnRegions();
                                                region2AreaBest = this->region2Area;
                                                area2RegionBest = this->area2Region;
                                            }

                                            // move happens, update bordering area
                                            getBorderingAreas(region);
                                            //getBorderingAreas(move);

                                            //std::cout << "--- NON-Local improvement (area, region)" << randomArea << "," << move << std::endl;
                                            //std::cout << "--- New Objective Function value: " << currentOBJ << std::endl;
                                            // step 4
                                        }
                                    }
                                }
                            }
                            // release lock for region[move], awake threads[move]
                        }
                    }
                    if (moved) {
                        break;
                    }
                }
                if (nothing_can_do) {
                    break;
                }
            }
        }
    }
    this->objInfo = bestOBJ;
    this->region2Area = region2AreaBest;
    this->area2Region = area2RegionBest;
    objective_function->UpdateRegions();
}

////////////////////////////////////////////////////////////////////////////////
////// AZP Tabu
////////////////////////////////////////////////////////////////////////////////

// Compares two intervals according to staring times.
bool CompareTabu(std::pair<std::pair<int, int>, double> i1,
                 std::pair<std::pair<int, int>, double> i2)
{
    return (i1.second < i2.second);
}


void AZPTabu::LocalImproving()
{
    // Tabu search algorithm for Openshaws AZP-tabu (1995)
    double aspireOBJ = this->objInfo;  // local best obj
    double currentOBJ = this->objInfo; // global best obj

    boost::unordered_map<std::pair<int, int>, bool> tabuDict;
    std::vector<std::pair<int, int> > tabuList;

    boost::unordered_map<std::pair<int, int>, double>::iterator it;

    // always remember the best result even it's just an interim one
    BasicMemory basicMemory;
    basicMemory.updateBasicMemory(this->objInfo, this->returnRegions());

    int c = 1;
    double epsilon = 1e-10;

    while (c < convTabu) {
        this->allCandidates();

        if (this->neighSolutions.size() == 0) {
            c += convTabu;
        } else {
            int minFound = 0;
            c += 1;

            // find  best feasible move not prohibited
            bool find_global = false;
            std::pair<int, int> move;
            double obj4Move  = 0; // current best objective from feasible moves

            while (!find_global && !neighSolObjs.empty()) {
                double minObj = neighSolObjs.top();
                neighSolObjs.pop();

                for (it = neighSolutions.begin(); it != neighSolutions.end(); ++it) {
                    if (it->second == minObj) {
                        // check if in tabuDict, [a, r] should not be in the tabu list
                        if (tabuDict.find(it->first) != tabuDict.end()) {
                            continue;
                        }
                        // check connectivity
                        int a = it->first.first;
                        int from = area2Region[a];
                        // move "a" to region "r", check if "a" can be removed from 'from'
                        if (objective_function->checkFeasibility(from, a)) {
                            find_global = true;
                            move = it->first;
                            obj4Move = minObj;
                            break;
                        }
                    }
                }
            }

            if (obj4Move == 0 && find_global == false) {
                // no more feasible move
                c = convTabu;
                break;
            }
            // best move in set
            if (currentOBJ - obj4Move >= epsilon) {
                minFound = 1; // use the global best move if improvement can be made
                //std::cout << move.first << "," << area2Region[move.first] << "," << move.second << "," << obj4Move << ",,,,,";
            } else if (tabuList.size() > 0){
                // if no improving move can be made,
                // then check tabu list for aspirational move
                // all valid tabu move should be in neighSolutions
                double best_tabuobj = aspireOBJ;
                std::pair<int, int> best_tabumove;

                for (size_t j=0; j<tabuList.size(); ++j) {
                    std::pair<int, int> m = tabuList[j];
                    if (neighSolutions.find(m) != neighSolutions.end()) {
                        // valid tabu move only exists in neighSolutions
                        double obj = neighSolutions[m];
                        if (obj < best_tabuobj && aspireOBJ - obj >= epsilon) {
                            // also need to check if this move breaks contiguity
                            if (objective_function->checkFeasibility(area2Region[m.first], m.first)) {
                                // tabu move improves local beset objectives
                                best_tabuobj = obj;
                                best_tabumove = m;
                                //break;
                            }
                        }
                    }
                }

                if (aspireOBJ - best_tabuobj >= epsilon) {
                    //std::cout << move.first << "," << area2Region[move.first] << "," << move.second << "," << obj4Move << "," << best_tabumove.first << "," << area2Region[best_tabumove.first] << "," << best_tabumove.second << "," << best_tabuobj <<",";
                    obj4Move = best_tabuobj;
                    move = best_tabumove;
                    minFound = 1;
                    // don't remove it from tabuList, since it's reversed move
                    // will be adde to tabuList
                }
            }

            // if none improvement can be made, then go with the best global move
            // even there is no improvement
            int area = move.first;
            int oldRegion = area2Region[area];
            int region = move.second;

            // Add the reverse of current move to the tabu list.
            // always added to the front and remove the last one
            std::pair<int, int> add_tabu(area, oldRegion);
            if (tabuDict.find(add_tabu) != tabuDict.end()) {
                // if it already existed in tabuList,
                // remove it first then add to the front
                tabuList.erase(std::remove(tabuList.begin(), tabuList.end(), add_tabu), tabuList.end());
            }
            tabuList.insert(tabuList.begin(), add_tabu);
            tabuDict[add_tabu] = true;
            if ((int)tabuList.size() > tabuLength) {
                std::pair<int, int> pop_tabu = tabuList.back();
                tabuDict.erase(pop_tabu);
                tabuList.pop_back();
            }


            // implement move
            region2Area[oldRegion].erase(area);
            region2Area[region].insert(std::make_pair(area, true));
            area2Region[area] = region;

            // update objective
            objective_function->UpdateRegion(region);
            objective_function->UpdateRegion(oldRegion);
            //double ssd = objective_function->GetValue();
            //double raw_ssd = objective_function->GetRawValue();
            //std::cout << area << "," << oldRegion << "," << region << "," << obj4Move << "," << currentOBJ << "," << aspireOBJ << std::endl;

            // update feasible neighboring set

            // save interim best result
            if (minFound == 1) {
                if (currentOBJ - obj4Move > epsilon) {
                    currentOBJ = obj4Move; // update global best, don't get worse
                    basicMemory.updateBasicMemory(obj4Move, this->returnRegions());
                    c = 1; // only reset conv when global objective improved
                }
            }
            this->objInfo = obj4Move;
            aspireOBJ = obj4Move; // update local best
        }
    }
    this->objInfo = basicMemory.objInfo;
    this->regions = basicMemory.regions;
}


void AZPTabu::allCandidates()
{
    // Select neighboring solutions.
    boost::unordered_map<int, bool>::iterator it;
    std::set<int>::iterator moves_it;

    neighSolutions.clear();
    neighSolObjs.clear();

    for (int r=0; r<p; ++r) {
        // get bordering area in each region "r"
        getBorderingAreas(r);
        boost::unordered_map<int, bool>& areas = region2Area[r];
        for (it = areas.begin(); it != areas.end(); ++it) {
            int a = it->first;
            if (it->second == true) {
                // processing boarding area, find possible moves
                std::set<int> moves = getPossibleMove(a);
                for (moves_it = moves.begin(); moves_it != moves.end(); ++moves_it) {
                    int move = *moves_it;
                    double obj = objective_function->TabuSwap(a, r, move);
                    neighSolutions[std::make_pair(a, move)] = obj;
                    neighSolObjs.push(obj);
                }
            }
        }
    }

}

void AZPTabu::updateNeighSolution(int area, int from, int to)
{
    // rebuild  neighSolutions that related to region "from" and "to"
    std::vector<std::pair<int, int> > removed_keys;

    boost::unordered_map<std::pair<int, int>, double>::iterator it;
    for (it = neighSolutions.begin(); it != neighSolutions.end(); ++it) {
        if (it->first.second == from || it ->first.second == to) {
            removed_keys.push_back(it->first);
        }
    }

    for (size_t i=0; i<removed_keys.size(); ++i) {
        neighSolutions.erase(removed_keys[i]);
    }

    boost::unordered_map<int, bool>::iterator ait;
    std::set<int>::iterator moves_it;

    std::vector<int> rr;
    rr.push_back(from);
    rr.push_back(to);

    for (size_t i=0; i<rr.size(); ++i) {
        int r = rr[i];
        getBorderingAreas(r);
        boost::unordered_map<int, bool>& areas = region2Area[r];
        for (ait = areas.begin(); ait != areas.end(); ++ait) {
            int a = ait->first;
            if (ait->second == true) { // boarding area
                std::set<int> moves = getPossibleMove(a);
                for (moves_it = moves.begin(); moves_it != moves.end(); ++moves_it) {
                    int move = *moves_it;
                    if (a != area) { // ignore already just moved area
                        double obj = objective_function->TabuSwap(a, r, move);
                        neighSolutions[std::make_pair(a, move)] = obj;
                    }
                }
            }
        }
    }

    neighSolObjs.clear();
    for (it = neighSolutions.begin(); it != neighSolutions.end(); ++it) {
        neighSolObjs.push(it->second);
    }
}
