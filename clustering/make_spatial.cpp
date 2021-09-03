#include <assert.h>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <iostream>
#include <boost/unordered_map.hpp>

#include "../weights/GeodaWeight.h"
#include "make_spatial.h"

MakeSpatialComponent::MakeSpatialComponent(int cid, const std::vector<int>& elements,
                                               GeoDaWeight* weights, std::map<int, int>& cluster_dict)
: cid(cid), elements(elements), weights(weights), cluster_dict(cluster_dict)
{
    int num_elements = (int)elements.size();
    isSingleton = num_elements == 1;
    isIsland = isSingleton && weights->GetNeighbors(elements[0]).empty();
    
    isSurroundedSingleton = false;
    if (isSingleton) {
        std::vector<long> nbrs = weights->GetNeighbors(elements[0]);
        boost::unordered_map<long ,bool> nbr_dict;
        for (int i = 0; i < (int)nbrs.size(); ++i) {
            if (elements[0] != (int)nbrs[i]) {
                nbr_dict[cluster_dict[(int)nbrs[i]]] = true;
            }
        }
        isSurroundedSingleton = nbr_dict.size() == 1;
    }
    
    for (int i = 0; i < num_elements; ++i) {
        elements_dict[elements[i]] = true;
    }
}

MakeSpatialComponent::~MakeSpatialComponent()
{
}

bool MakeSpatialComponent::Has(int eid)
{
    return elements_dict[eid];
}

void MakeSpatialComponent::Merge(MakeSpatialComponent *comp)
{
    std::vector<int> new_elements = comp->GetElements();
    for (int i = 0; i < (int)new_elements.size(); ++i) {
        this->elements.push_back(new_elements[i]);
        this->elements_dict[new_elements[i]] = true;
    }
}

MakeSpatialCluster::MakeSpatialCluster(int cid, const std::vector<int>& elements,
                                           GeoDaWeight* weights,
                                           std::map<int, int>& cluster_dict)
: cid(cid), elements(elements), cluster_dict(cluster_dict), weights(weights), core(0)
{
    int num_elements = (int)elements.size();
    
    // create components from elements
    std::map<int, bool> visited;
    for (int i = 0; i < num_elements; ++i) {
        int eid = elements[i];
        if (visited[eid]) {
            continue;
        }
        std::vector<int> component;
        visited[eid] = true;
        
        // find contiguous neighbors
        std::stack<int> stack;
        stack.push(eid);
        
        while (!stack.empty()) {
            int tmp_id = stack.top();
            stack.pop();
            component.push_back(tmp_id);
            std::vector<long> nbrs = weights->GetNeighbors(tmp_id);
            for (int j = 0; j < (int)nbrs.size(); ++j) {
                int neighbor = (int)nbrs[j];
                if (cluster_dict[neighbor] == this->cid && !visited[neighbor]) {
                    visited[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        
        
        MakeSpatialComponent* c = new MakeSpatialComponent(this->cid, component, weights, cluster_dict);
        components.push_back(c);
        
        // create lookup table
        for (int j = 0; j < (int)component.size(); ++j) {
            int eid = component[j];
            component_dict[eid] = c;
        }
        
        if (core == 0 || core->GetSize() < c->GetSize()) {
            core = c;
        }
    }
}

MakeSpatialCluster::~MakeSpatialCluster()
{
    // clean up components
    for (int i = 0; i < (int)components.size(); ++i) {
        delete components[i];
    }
}

MakeSpatialComponent* MakeSpatialCluster::GetComponent(int eid)
{
    // get component by giving an element id
    return component_dict[eid];
}

std::vector<int> MakeSpatialCluster::GetCoreElements()
{
    return core->GetElements();
}

bool MakeSpatialCluster::BelongsToCore(int eid)
{
    return core->Has(eid);
}

int MakeSpatialCluster::GetCoreSize()
{
    return core->GetSize();
}

std::vector<MakeSpatialComponent*> MakeSpatialCluster::GetSurroundedSingletons()
{
    std::vector<MakeSpatialComponent*> result;
    
    // find singleton components that connect to only one cluster (island?)
    for (int i = 0; i < (int)components.size(); ++i) {
        MakeSpatialComponent* comp = components[i];
        if (comp != core && comp->isSurroundedSingleton) {
            result.push_back(comp);
        }
    }
    
    return result;
}

std::vector<int> MakeSpatialCluster::GetComponentSize()
{
    // return the size of all components (as a set) except the core
    std::set<int> sz_set;
    
    for (int i = 0; i < (int)components.size(); ++i) {
        MakeSpatialComponent* comp = components[i];
        if (comp != this->core) {
            sz_set.insert(comp->GetSize());
        }
    }
    
    std::set<int>::iterator it;
    std::vector<int> result;
    for (it = sz_set.begin(); it != sz_set.end(); ++it) {
        result.push_back(*it);
    }
    return result;
}

int MakeSpatialCluster::GetSmallestComponentSize()
{
    int result = -1;
    for (int i = 0; i < (int)this->components.size(); ++i) {
        if (components[i] != core) {
            int sz = components[i]->GetSize();
            if (result < 0 || sz < result) {
                result = sz;
            }
        }
    }
    return result;
}

std::vector<MakeSpatialComponent*> MakeSpatialCluster::GetComponentsBySize(int component_size)
{
    std::vector<MakeSpatialComponent*> result;
    
    for (int i = 0; i < (int)this->components.size(); ++i) {
        MakeSpatialComponent* comp = this->components[i];
        if (comp != core && comp->GetSize() == component_size) {
            result.push_back(comp);
        }
    }
    
    return result;
}

void MakeSpatialCluster::MergeComponent(MakeSpatialComponent* from, MakeSpatialComponent* to)
{
    for (int i = 0; i < (int)components.size(); ++i) {
        if (components[i] == to) {
            // merge "from" -> "to"
            to->Merge(from);
            
            // update lookup table
            std::vector<int> new_elements = from->GetElements();
            for (int j = 0; j < (int)new_elements.size(); ++j) {
                int eid = new_elements[j];
                component_dict[eid] = to;
            }
            
            // update core component
            if (to->GetSize() > core->GetSize()) {
                core = to;
            }
            break;
        }
    }
}

void MakeSpatialCluster::RemoveComponent(MakeSpatialComponent *comp)
{
    // update lookup table first
    if (comp->GetClusterId() != this->cid) {
        // There is no need to erase elements from lookup dict if the componenet is removed
        // from the same cluster that it has been merged into
        std::vector<int> removed_elements = comp->GetElements();
        for (int i = 0; i < (int)removed_elements.size(); ++i) {
            int eid = removed_elements[i];
            component_dict.erase(eid);
        }
    }
    
    // component has been merged to another component and will be removed from memory!
    for (int i = 0; i < (int)components.size(); ++i) {
        if (components[i] != core && components[i] == comp) {
            delete comp;
            components.erase(components.begin() + i);
            break;
        }
    }
}

MakeSpatial::MakeSpatial(int num_obs, const std::vector<std::vector<int> >& clusters, GeoDaWeight* weights)
: num_obs(num_obs), clusters(clusters), weights(weights), valid(true)
{
    num_clusters = (int)clusters.size();
    
    std::vector<int>::iterator it;
    for (int i=0; i < num_clusters; ++i) {
        std::vector<int> cluster = clusters[i];
        for (it = cluster.begin(); it != cluster.end(); ++it) {
            int j = *it;
            cluster_dict[j] = i;
        }
    }
    
    // check if clusters are valid
    if (cluster_dict.size() != num_obs) {
        valid = false;
    }
    
    // create MakeSpatialCluster for each cluster
    for (int i=0; i < num_clusters; ++i) {
        sk_clusters.push_back(new MakeSpatialCluster(i, clusters[i], weights, cluster_dict));
    }
}

MakeSpatial::~MakeSpatial()
{
    for (int i = 0; i < (int)sk_clusters.size(); ++i) {
        delete sk_clusters[i];
    }
}

struct ClusterSmall {
    bool operator()(MakeSpatialCluster* left, MakeSpatialCluster* right) const {
        return left->GetCoreSize() > right->GetCoreSize();
    }
};

int MakeSpatial::GetSmallestComponentSize()
{
    int result = -1;
    for (int i = 0; i < num_clusters; ++i) {
        int sz = sk_clusters[i]->GetSmallestComponentSize();
        if (result < 0 || (sz > 0 && sz < result) ) {
            result = sz;
        }
    }
    return result;
}

std::vector<MakeSpatialCluster*> MakeSpatial::GetClustersByComponentSize(int sz)
{
    std::vector<MakeSpatialCluster*> result;
    
    for (int i = 0; i < num_clusters; ++i) {
        if (sk_clusters[i]->GetSmallestComponentSize() == sz) {
            result.push_back(sk_clusters[i]);
        }
    }
    
    return result;
}

void MakeSpatial::Run()
{
    if (!valid) return;
    
    // find all surrounded singleton cases, assign to nearby cluster
    for (int i = 0; i < num_clusters; ++i) {
        std::vector<MakeSpatialComponent*> moved_comps = sk_clusters[i]->GetSurroundedSingletons();
        for (int j=0; j < (int)moved_comps.size(); ++j) {
            MakeSpatialComponent* moved_comp = moved_comps[j];
            this->MoveComponent(moved_comp);
        }
    }
    
    // process components (small->large) in clusters
    int N;
    while ((N = GetSmallestComponentSize()) > 0) {
        // get clusters which has components with size = N
        std::vector<MakeSpatialCluster*> cands = GetClustersByComponentSize(N);

        // process cluster by its size (smallest to largets)
        while (!cands.empty()) {
            std::make_heap(cands.begin(), cands.end(), ClusterSmall());
            std::pop_heap(cands.begin(), cands.end());
            
            MakeSpatialCluster* c = cands.back();
            cands.pop_back();
            
            std::vector<MakeSpatialComponent*> moved_comps = c->GetComponentsBySize(N);
            
            for (int j=0; j < (int)moved_comps.size(); ++j) {
                MakeSpatialComponent* moved_comp = moved_comps[j];
                this->MoveComponent(moved_comp);
            }
        }
    }
}

void MakeSpatial::MoveComponent(MakeSpatialComponent *comp)
{
    // move a component by merging it to the largest nearby cluster
    int largest_size = 0;
    MakeSpatialComponent* best_to = 0;
    
    // get target cluster ids (component) from neighbors
    std::set<int> target_clusters;
    std::vector<int> elements = comp->GetElements();
    for (int i = 0; i < (int)elements.size(); ++i) {
        int eid = elements[i];
        std::vector<long> nbrs = weights->GetNeighbors(eid);
        for (int j = 0; j < (int)nbrs.size(); ++j) {
            int nbr = (int)nbrs[j];
            if (!comp->Has(nbr)) {
                int target = this->cluster_dict[nbr];
                // no need "to != from" since the component could be merged to the core
                MakeSpatialComponent* to = sk_clusters[target]->GetComponent(nbr);
                if (to!= NULL && to != comp && to->GetSize() > largest_size) {
                    best_to = to;
                    largest_size = to->GetSize();
                }
            }
        }
    }
    // move component to the largest cluster
    if (best_to != 0) {
        UpdateComponent(comp, best_to);
    } else {
        //std::cout << "Can't move component" << std::endl;
        valid = false;
    }
}


void MakeSpatial::UpdateComponent(MakeSpatialComponent* moved_comp, MakeSpatialComponent* target)
{
    // update cluster_dict first!!
    std::vector<int> elements = moved_comp->GetElements();
    
    for (int i = 0; i < (int)elements.size(); ++i) {
        int eid = elements[i];
        cluster_dict[eid] = target->GetClusterId();
    }
    
    int to = target->GetClusterId();
    int from = moved_comp->GetClusterId();
    
    sk_clusters[to]->MergeComponent(moved_comp, target);
    sk_clusters[from]->RemoveComponent(moved_comp);
}

std::vector<std::vector<int> > MakeSpatial::GetClusters()
{
    // verify clustering result
    int total_core_obs = 0;
    for (int i = 0; i < num_clusters; ++i) {
        MakeSpatialCluster* skc = sk_clusters[i];
        total_core_obs += skc->GetCoreSize();
    }
    
    if (total_core_obs != num_obs) {
        valid = false;
        return this->clusters;
    }
    
    std::vector<std::vector<int> > result;
    for (int i = 0; i < num_clusters; ++i) {
        MakeSpatialCluster* skc = sk_clusters[i];
        result.push_back(skc->GetCoreElements());
    }
    return result;
}
