
#include <cstddef> // for std::ptrdiff_t
#include <limits> // for std::numeric_limits<...>::infinity()
#include <algorithm> // for std::fill_n
#include <stdexcept> // for std::runtime_error
#include <string> // for std::string

#include "fastcluster.h"

using namespace fastcluster;

void fastcluster::MST_linkage_core(const t_index N, const t_float * const D,
                                   cluster_result & Z2) {
/*
    N: integer, number of data points
    D: condensed distance matrix N*(N-1)/2
    Z2: output data structure

    The basis of this algorithm is an algorithm by Rohlf:

    F. James Rohlf, Hierarchical clustering using the minimum spanning tree,
    The Computer Journal, vol. 16, 1973, p. 93–95.
*/
  t_index i;
  t_index idx2;
  doubly_linked_list active_nodes(N);
  auto_array_ptr<t_float> d(N);

  t_index prev_node;
  t_float min;

  // first iteration
  idx2 = 1;
  min = std::numeric_limits<t_float>::infinity();
  for (i=1; i<N; ++i) {
    d[i] = D[i-1];
    if (d[i] < min) {
      min = d[i];
      idx2 = i;
    }
    else if (fc_isnan(d[i]))
      throw (nan_error());
  }
  Z2.append(0, idx2, min);

  for (t_index j=1; j<N-1; ++j) {
    prev_node = idx2;
    active_nodes.remove(prev_node);

    idx2 = active_nodes.succ[0];
    min = d[idx2];
    for (i=idx2; i<prev_node; i=active_nodes.succ[i]) {
      t_float tmp = D_(i, prev_node);
      if (tmp < d[i])
        d[i] = tmp;
      else if (fc_isnan(tmp))
        throw (nan_error());
      if (d[i] < min) {
        min = d[i];
        idx2 = i;
      }
    }
    for (; i<N; i=active_nodes.succ[i]) {
      t_float tmp = D_(prev_node, i);
      if (d[i] > tmp)
        d[i] = tmp;
      else if (fc_isnan(tmp))
        throw (nan_error());
      if (d[i] < min) {
        min = d[i];
        idx2 = i;
      }
    }
    Z2.append(prev_node, idx2, min);
  }
}
