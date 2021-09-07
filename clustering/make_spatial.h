/**
 * GeoDa TM, Copyright (C) 2011-2015 by Luc Anselin - all rights reserved
 *
 * This file is part of GeoDa.
 *
 * GeoDa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GeoDa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Created: 8/9/2021 lixun910@gmail.com
 */

#ifndef __GEODA_CENTER_SPATIAL_KMEANS_H__
#define __GEODA_CENTER_SPATIAL_KMEANS_H__

#include <vector>
#include <map>

class GeoDaWeight;
class MakeSpatialCluster;
class MakeSpatial;

class MakeSpatialComponent
{
public:
    MakeSpatialComponent(int cid, const std::vector<int>& elements,
                           GeoDaWeight* weights, std::map<int, int>& cluster_dict);
    virtual ~MakeSpatialComponent();
    
    int GetClusterId() { return cid; }
    void SetClusterId(int cid) { this->cid = cid; }
    
    int GetSize() { return (int)elements.size(); }
    std::vector<int> GetElements() { return elements; }
    
    void Merge(MakeSpatialComponent* comp);
    
    bool Has(int eid);
    
    bool isIsland;
    bool isSingleton;
    bool isSurroundedSingleton;
    
protected:
    int cid;
    std::vector<int> elements;
    GeoDaWeight* weights;
    std::map<int, int>& cluster_dict;
    std::map<int, bool> elements_dict;
};

class MakeSpatialCluster
{
public:
    MakeSpatialCluster(int cid, const std::vector<int>& elements, GeoDaWeight* weights,
                         std::map<int, int>& cluster_dict);
    virtual ~MakeSpatialCluster();
    
    // find singletons that only connect to one cluster
    std::vector<MakeSpatialComponent*> GetSurroundedSingletons();
    
    std::vector<MakeSpatialComponent*> GetComponentsBySize(int component_size);
    
    void MergeComponent(MakeSpatialComponent* from, MakeSpatialComponent* to);
    
    void RemoveComponent(MakeSpatialComponent* comp);
    
    std::vector<int> GetCoreElements();
    
    int GetCoreSize();
    
    int GetComponentSize(int eid);
    
    int GetSmallestComponentSize();
    
    MakeSpatialComponent* GetComponent(int eid);
    
    bool BelongsToCore(int eid);

    std::vector<int> GetComponentSize();
    
protected:
    int cid;
    std::vector<int> elements;
    std::map<int, int>& cluster_dict;
    GeoDaWeight* weights;
    MakeSpatialComponent* core;
    std::vector<MakeSpatialComponent*> components;
    std::map<int, MakeSpatialComponent*> component_dict;
};

class MakeSpatial
{
public:
    MakeSpatial(int num_obs, const std::vector<std::vector<int> >& clusters,
                  GeoDaWeight* weights);
	virtual ~MakeSpatial();
    
    void Run();
    
    std::vector<std::vector<int> > GetClusters();
    
    bool IsValid() { return valid; }
    
protected:
    void UpdateComponent(MakeSpatialComponent* moved_comp,
                         MakeSpatialComponent* target);
    
    void MoveComponent(MakeSpatialComponent* comp);
    
    int GetSmallestComponentSize();
    
    std::vector<MakeSpatialCluster*> GetClustersByComponentSize(int sz);
    
protected:
    int num_obs;
    std::vector<std::vector<int> > clusters;
    GeoDaWeight* weights;
    bool valid;

    int num_clusters;
    
    std::map<int, int> cluster_dict;
    
    std::vector<MakeSpatialCluster*> sk_clusters;
};

#endif
