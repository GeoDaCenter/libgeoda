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
 */


// Note: Boost 1.52 and later is required for use of the Boost.Polygon
//   Voronoi Library.  Many thanks to Andrii Sydorchuk for contributing
//   this high-quality Voronoi Diagram library to Boost.
#include <algorithm>
#include <map>
#include <utility>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>
#include <boost/polygon/voronoi.hpp>
#include <boost/polygon/voronoi_builder.hpp>
#include <boost/polygon/voronoi_diagram.hpp>

#include "GalWeight.h"
#include "../GenUtils.h"
#include "../GenGeomAlgs.h"
#include "VoronoiUtils.h"

using namespace boost::polygon;

struct Point {
	int a;
	int b;
	Point (int x, int y) : a(x), b(y) {}
};

struct Segment {
	Point p0;
	Point p1;
	Segment (int x1, int y1, int x2, int y2) : p0(x1, y1), p1(x2, y2) {}
};

namespace boost {
	namespace polygon {

		template <>
		struct geometry_concept<Point> { typedef point_concept type; };

		template <>
		struct point_traits<Point> {
			typedef int coordinate_type;

			static inline coordinate_type get(const Point& point,
											  orientation_2d orient) {
				return (orient == HORIZONTAL) ? point.a : point.b;
			}
		};

		template <>
		struct geometry_concept<Segment> { typedef segment_concept type; };

		template <>
		struct segment_traits<Segment> {
			typedef int coordinate_type;
			typedef Point point_type;

			static inline point_type get(const Segment& segment,
										 direction_1d dir) {
				return dir.to_int() ? segment.p1 : segment.p0;
			}
		};
	}  // polygon
}  // boost

/**
 Note: Input is double centroids, but we then scale up to large integers
 to create the Voronoi diagram.  The output of the Voronoi diagram
 is large doubles.
 */

namespace Gda {
	namespace VoronoiUtils {
		typedef voronoi_builder<int> VB;
		typedef voronoi_diagram<double> VD;

		std::list<int>* getCellList(
							const VD::cell_type& cell,
							std::map<std::pair<int,int>,
							std::list<int>* >& pt_to_id_list,
							std::vector<std::pair<int,int> >& int_pts);
		bool isVertexOutsideBB(const VD::vertex_type& vertex,
							   const double& xmin, const double& ymin,
							   const double& xmax, const double& ymax);
		bool clipEdge(const VD::edge_type& edge,
					  std::vector<std::pair<int,int> >& int_pts,
					  const double& xmin, const double& ymin,
					  const double& xmax, const double& ymax,
					  double& x0, double& y0, double& x1, double& y1);
		bool clipInfiniteEdge(const VD::edge_type& edge,
							  std::vector<std::pair<int,int> >& int_pts,
							  const double& xmin, const double& ymin,
							  const double& xmax, const double& ymax,
							  double& x0, double& y0, double& x1, double& y1);
		bool clipFiniteEdge(const VD::edge_type& edge,
							std::vector<std::pair<int,int> >& int_pts,
							const double& xmin, const double& ymin,
							const double& xmax, const double& ymax,
							double& x0, double& y0, double& x1, double& y1);
	}
}

/** Input: double precision x/y coordinates, indexed by observation record id
 Output: list of list of duplicates
 */
void Gda::VoronoiUtils::FindPointDuplicates(const std::vector<double>& x,
											  const std::vector<double>& y,
										std::list<std::list<int> >& duplicates)
{
	typedef std::pair<int,int> int_pair;
	int num_obs = x.size();
	double x_orig_min=0, x_orig_max=0;
	double y_orig_min=0, y_orig_max=0;
	SampleStatistics::CalcMinMax(x, x_orig_min, x_orig_max);
	SampleStatistics::CalcMinMax(y, y_orig_min, y_orig_max);
	double orig_scale = std::max(x_orig_max-x_orig_min,
											  y_orig_max-y_orig_min);
	if (orig_scale == 0) orig_scale = 1;
	double big_dbl = 1073741824; // 2^30
	double p = (big_dbl/orig_scale);

	std::map<int, std::list<int> > dups;
	std::map<int, std::list<int> >::iterator dups_iter;
	std::map<int_pair, int> pt_map;
	std::map<int_pair, int>::iterator map_iter;
	std::vector<int> x_int(num_obs);
	std::vector<int> y_int(num_obs);
	bool duplicates_exist = false;
	for (int i=0; i<num_obs; i++) {
		x_int[i] = (int) ((x[i]-x_orig_min)*p);
		y_int[i] = (int) ((y[i]-y_orig_min)*p);
		int_pair key(std::make_pair(x_int[i], y_int[i]));
		map_iter = pt_map.find(key);
		if (map_iter == pt_map.end()) {
			pt_map[key] = i;
		} else {
			duplicates_exist = true;
			int ind1 = map_iter->second;
			int ind2 = i;
			dups_iter = dups.find(ind1);
			if (dups_iter == dups.end()) {
				std::list<int> l;
				l.push_back(ind1);
				l.push_back(ind2);
				dups[ind1] = l;
			} else {
				dups_iter->second.push_back(ind2);
			}
		}
	}
	duplicates.clear();
	if (duplicates_exist) {
		for (dups_iter = dups.begin(); dups_iter != dups.end(); dups_iter++) {
			duplicates.push_back(dups_iter->second);
		}
	}
}

std::list<int>* Gda::VoronoiUtils::getCellList(
				const VD::cell_type& cell,
				std::map<std::pair<int,int>, std::list<int>* >& pt_to_id_list,
				std::vector<std::pair<int,int> >& int_pts)
{
	std::map<std::pair<int,int>, std::list<int>* >::iterator iter;
	//iter = pt_to_id_list.find(std::make_pair(cell.point0().x(),
	//										 cell.point0().y()));
	iter = pt_to_id_list.find(int_pts[cell.source_index()]);
	if (iter == pt_to_id_list.end()) {
		return 0;
	}
	return iter->second;
}

bool Gda::VoronoiUtils::isVertexOutsideBB(const VD::vertex_type& vertex,
											const double& xmin,
											const double& ymin,
											const double& xmax,
											const double& ymax)
{
	double x = vertex.x();
	double y = vertex.y();
	return (x < xmin || x > xmax || y < ymin || y > ymax);
}

/** Clip both infinite and finite edges to bounding rectangle.
 return true if intersection or if edge is contained within bounding box,
 otherwise return false */
bool Gda::VoronoiUtils::clipEdge(const VD::edge_type& edge,
								   std::vector<std::pair<int,int> >& int_pts,
								   const double& xmin, const double& ymin,
								   const double& xmax, const double& ymax,
								   double& x0, double& y0,
								   double& x1, double& y1)
{
	if (edge.is_finite()) {
		return clipFiniteEdge(edge, int_pts, xmin, ymin, xmax, ymax,
							  x0, y0, x1, y1);
	} else {
		return clipInfiniteEdge(edge, int_pts, xmin, ymin, xmax, ymax,
								x0, y0, x1, y1);
	}
}

/** Clip infinite edge to bounding rectangle */
bool Gda::VoronoiUtils::clipInfiniteEdge(const VD::edge_type& edge,
									std::vector<std::pair<int,int> >& int_pts,
									const double& xmin, const double& ymin,
									const double& xmax, const double& ymax,
									double& x0, double& y0,
									double& x1, double& y1)
{
    const VD::cell_type& cell1 = *edge.cell();
    const VD::cell_type& cell2 = *edge.twin()->cell();
    double origin_x, origin_y, direction_x, direction_y;
    // Infinite edges could not be created by two segment sites.
    if (cell1.contains_point() && cell2.contains_point()) {
		double p1_x = (double) int_pts[cell1.source_index()].first;
		double p1_y = (double) int_pts[cell1.source_index()].second;
		double p2_x = (double) int_pts[cell2.source_index()].first;
		double p2_y = (double) int_pts[cell2.source_index()].second;
		origin_x = ((p1_x + p2_x) * 0.5);
		origin_y = ((p1_y + p2_y) * 0.5);
		direction_x = (p1_y - p2_y);
		direction_y = (p2_x - p1_x);
    } else {
		// This case should never happen for point maps.
		//LOG_MSG("Warning! one clipInfiniteEdge cells contains a segment!");
		return false;
    }
    double side = xmax - xmin;
    double koef =
		side / (std::max)(fabs(direction_x), fabs(direction_y));
    if (edge.vertex0() == NULL) {
		x0 = origin_x - direction_x * koef;
		y0 = origin_y - direction_y * koef;
    } else {
		x0 = edge.vertex0()->x();
		y0 = edge.vertex0()->y();
    }
    if (edge.vertex1() == NULL) {
		x1 = origin_x + direction_x * koef;
		y1 = origin_y + direction_y * koef;
    } else {
		x1 = edge.vertex1()->x();
		y1 = edge.vertex1()->y();
    }
	return GenGeomAlgs::ClipToBB(x0, y0, x1, y1, xmin, ymin, xmax, ymax);
}

/** Clip finite edge to bounding rectangle */
bool Gda::VoronoiUtils::clipFiniteEdge(const VD::edge_type& edge,
									std::vector<std::pair<int,int> >& int_pts,
									const double& xmin, const double& ymin,
									const double& xmax, const double& ymax,
									double& x0, double& y0,
									double& x1, double& y1)
{
	// we know that edge is finite, so both vertex0 and vertex1 are defined
	x0 = edge.vertex0()->x();
	y0 = edge.vertex0()->y();
	x1 = edge.vertex1()->x();
	y1 = edge.vertex1()->y();
	return GenGeomAlgs::ClipToBB(x0, y0, x1, y1, xmin, ymin, xmax, ymax);
}

/** If false returned, then an unexpected error.  Otherwise, neighbor map
 created successfully.  The presence of duplicates is indicated in
 duplicates_exists and the list of duplicates is filled in.
 */
bool Gda::VoronoiUtils::PointsToContiguity(const std::vector<double>& x,
									const std::vector<double>& y,
									bool queen,
									std::vector<std::set<int> >& nbr_map)
{
	//LOG_MSG("Entering Gda::VoronoiUtils::PointsToContiguity");
	typedef std::pair<int,int> int_pair;
	typedef std::list<int> id_list;

	int num_obs = x.size();
	double x_orig_min=0, x_orig_max=0;
	double y_orig_min=0, y_orig_max=0;
	SampleStatistics::CalcMinMax(x, x_orig_min, x_orig_max);
	SampleStatistics::CalcMinMax(y, y_orig_min, y_orig_max);
	double orig_scale = std::max(x_orig_max-x_orig_min,
											  y_orig_max-y_orig_min);
	if (orig_scale == 0) orig_scale = 1;
	double big_dbl = 1073741824; // 2^30
	double p = (big_dbl/orig_scale);

	// Add 2% offset to the bounding rectangle
	const double bb_pad = 0.02;
	// note data has been translated to origin and scaled
	double bb_xmin = -bb_pad*big_dbl;
	double bb_xmax = (x_orig_max-x_orig_min)*p + bb_pad*big_dbl;
	double bb_ymin = -bb_pad*big_dbl;
	double bb_ymax = (y_orig_max-y_orig_min)*p + bb_pad*big_dbl;

	std::map<int_pair, id_list* > pt_to_id_list;
	// for each unique point, the list of cells at that point
	std::map<int_pair, id_list* >::iterator pt_to_id_list_iter;

	std::map<int_pair, std::set<id_list* > > pt_to_nbr_sets;
	// for each unique point, the set of lists of points that are neighbors
	std::map<int_pair, std::set<id_list* > >::iterator pt_to_nbr_sets_iter;

	std::vector<int_pair> int_pts(num_obs);
	for (int i=0; i<num_obs; i++) {
		int_pts[i].first = (int) ((x[i]-x_orig_min)*p);
		int_pts[i].second = (int) ((y[i]-y_orig_min)*p);
	}

	for (int i=0; i<num_obs; i++) {
		pt_to_id_list_iter = pt_to_id_list.find(int_pts[i]);
		if (pt_to_id_list_iter == pt_to_id_list.end()) {
			pt_to_id_list[int_pts[i]] = new id_list;
		}
		pt_to_id_list[int_pts[i]]->push_back(i);
	}

	nbr_map.clear();
	nbr_map.resize(num_obs);

	VD vd;
	VB vb;
	for (int i=0; i<num_obs; i++) {
		//int index = vb.insert_point(int_pts[i].first, int_pts[i].second);
		vb.insert_point(int_pts[i].first, int_pts[i].second);
	}
	vb.construct(&vd);

	for (VD::const_cell_iterator it = vd.cells().begin();
		 it != vd.cells().end(); ++it) {
		const VD::cell_type &cell = *it;

		int_pair key = int_pts[cell.source_index()];

		std::set<id_list* >& nbr_set = pt_to_nbr_sets[key];

		const VD::edge_type* edge = cell.incident_edge();
		typedef std::list<const VD::vertex_type*> v_list;
		v_list verts;
		do {
			id_list* nbr_list = getCellList(*(edge->twin()->cell()),
											pt_to_id_list, int_pts);
			if (!nbr_list) {
				return false;
			}

			double x0=0, y0=0, x1=0, y1=0;
			if (clipEdge(*edge, int_pts,
						 bb_xmin, bb_ymin, bb_xmax, bb_ymax,
						 x0, y0, x1, y1)) {
				nbr_set.insert(nbr_list);
			}

			if (queen) { // add all cells that share each edge vertex
				if (edge->vertex0() &&
					!isVertexOutsideBB(*edge->vertex0(), bb_xmin, bb_ymin,
									   bb_xmax, bb_ymax)) {
						verts.push_back(edge->vertex0());
					}
				if (edge->vertex1() &&
					!isVertexOutsideBB(*edge->vertex1(), bb_xmin, bb_ymin,
									   bb_xmax, bb_ymax)) {
					verts.push_back(edge->vertex1());
				}
			}
			edge = edge->next();
		} while (edge != cell.incident_edge());

		// add all neighboring cells.  List will be empty if !queen
		for (v_list::iterator it = verts.begin(); it != verts.end(); it++) {
			const VD::edge_type *edge =
				(*it)->incident_edge();
			do {
				id_list* nbr_list = getCellList(*(edge->cell()),
												 pt_to_id_list, int_pts);
				if (!nbr_list) {
					return false;
				}
				nbr_set.insert(nbr_list);

				edge = edge->rot_next();
			} while (edge != (*it)->incident_edge());
		}
	}

	for (int i=0; i<num_obs; i++) {
		std::set<id_list* >& nbr_set = pt_to_nbr_sets[int_pts[i]];
		//LOG_MSG(wxString::Format("obs %d has %d neighbors",
		//        i, nbr_set.size()));

		for (std::set<id_list* >::iterator set_iter = nbr_set.begin();
			 set_iter != nbr_set.end(); set_iter++) {
			for (id_list::iterator list_iter = (*set_iter)->begin();
				 list_iter != (*set_iter)->end(); list_iter++) {
				if (*list_iter != i) nbr_map[i].insert(*list_iter);
			}
		}
		// for each set of points, make every point a neighbor of the
		// other points in the set, excluding self
		id_list* pl = pt_to_id_list[int_pts[i]];
		for (id_list::iterator it=pl->begin(); it!=pl->end(); it++) {
			if (*it != i) nbr_map[i].insert(*it);
		}
	}

	for (pt_to_id_list_iter = pt_to_id_list.begin();
		 pt_to_id_list_iter != pt_to_id_list.end();
		 pt_to_id_list_iter++) {
		delete pt_to_id_list_iter->second;
	}

	return true;
}

GalElement* Gda::VoronoiUtils::NeighborMapToGal(
										std::vector<std::set<int> >& nbr_map)
{
	if (nbr_map.size() == 0) return 0;
	GalElement* gal = new GalElement[nbr_map.size()];
	if (!gal) return 0;
	for (int i=0, iend=nbr_map.size(); i<iend; i++) {
		gal[i].SetSizeNbrs(nbr_map[i].size());
		long cnt = 0;
		for (std::set<int>::iterator it=nbr_map[i].begin();
			 it != nbr_map[i].end(); it++) {
			gal[i].SetNbr(cnt++, *it);
		}
	}
	return gal;
}
