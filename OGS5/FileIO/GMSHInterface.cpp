/*
 * GMSHInterface.cpp
 *
 *  Created on: Apr 29, 2010
 *      Author: TF
 */

#include <fstream>
#include <vector>
#include <list>

// Base
#include "swap.h"

// FileIO
#include <GMSHInterface.h>

// GEOLIB
#include "Point.h"
#include "Polyline.h"
#include "Polygon.h"
#include "SimplePolygonTree.h"
#include "QuadTree.h"

namespace FileIO {

GMSHInterface::GMSHInterface (const std::string &fname) :
	_n_pnt_offset (0),
	_n_lines (0), // _n_line_loops (0),
	_n_plane_sfc(0)
{
	// open file
	_out.open (fname.c_str());
	// check file stream
	if (!_out) {
		std::cerr << "could not open file " << fname << std::endl;
		return;
	}

	_out << "// GMSH input file created by OpenGeoSys" << std::endl;
	_out << std::endl;
	_out.precision (20);
}

GMSHInterface::~GMSHInterface ()
{
	_out.close ();
}

void GMSHInterface::writeGMSHPoints(const std::vector<GEOLIB::Point*>& pnt_vec)
{
	// write points
//	float characteristic_len (10.0);
	size_t n (pnt_vec.size());
	for (size_t k(0); k<n; k++) {
		_out << "Point(" << _n_pnt_offset + k << ") = {" << (*(pnt_vec[k]))[0] << ","
			<< (*(pnt_vec[k]))[1] << "," << (*(pnt_vec[k]))[2]
//			<< "," << characteristic_len
			<< "};" << std::endl;
	}
	_n_pnt_offset += n;
}

void GMSHInterface::writeGMSHPolyline (const GEOLIB::Polyline* ply)
{
	size_t s (ply->getNumberOfPoints());
	// write line segments (= Line) of the polyline
	for (size_t j(0); j<s-1; j++) {
		_out << "Line(" << _n_lines+j << ") = {" << _n_pnt_offset + ply->getPointID(j) << ","
				<< _n_pnt_offset + ply->getPointID(j+1) << "};" << std::endl;
	}
	// write the line segments contained in the polyline (=Line Loop)
	_out << "Line Loop (" << _n_lines + s - 1 << ") = {";
	for (size_t j(0); j<s-2; j++) {
		_out << _n_lines + j << ",";
	}
	_out << _n_lines + s-2 << "};" << std::endl;
	_n_lines += s;
}

void GMSHInterface::writeGMSHPolylines(const std::vector<GEOLIB::Polyline*>& ply_vec)
{
	size_t n (ply_vec.size());
	for (size_t k(0); k<n; k++) {
		// write k-th polyline
		writeGMSHPolyline (ply_vec[k]);
	}
	_out << std::endl;
}

void GMSHInterface::writeGMSHPolygon(const GEOLIB::Polygon& polygon)
{
	writeGMSHPolyline (&polygon);
	_polygon_list.push_back (_n_lines-1);
}

bool GMSHInterface::writeGMSHInputFile(const std::string &proj_name, const GEOLIB::GEOObjects& geo)
{
	std::cerr << "GMSHInterface::writeGMSHInputFile " << std::endl;
	std::cerr << "get data from geo ... " << std::flush;
	// get data from geo
	const std::vector<GEOLIB::Point*> *pnts (geo.getPointVec (proj_name));
	const std::vector<GEOLIB::Polyline*> *plys (geo.getPolylineVec (proj_name));
	std::cerr << "ok" << std::endl;

	// check file stream
	if (!_out) return false;

	// write points
	writeGMSHPoints (*pnts);

	// write Polylines
	writeGMSHPolylines (*plys);
	std::cerr << "ok" << std::endl;

	return true;
}

void GMSHInterface::writePlaneSurface ()
{
	_out << "Plane Surface (" << _n_plane_sfc << ") = {" << std::flush;
	std::list<size_t>::const_iterator it (_polygon_list.begin());
	_out << *it << std::flush;
	for (it++; it != _polygon_list.end(); it++)
		_out << ", " << *it << std::flush;
	_out << "};" << std::endl;
	_n_plane_sfc++;
	_polygon_list.clear();
}

void GMSHInterface::writeAllDataToGMSHInputFile (GEOLIB::GEOObjects& geo,
		std::vector<std::string> const & selected_geometries,
		size_t number_of_point_per_quadtree_node,
		double mesh_density_scaling, double mesh_density_scaling_station_pnts)
{
	// check file stream
	if (! _out) return;

	std::cout << "GMSHInterface::writeGMSHInputFile " << std::endl;

	std::vector<GEOLIB::Point*> all_points;
	std::vector<GEOLIB::Polyline*> all_polylines;
	std::vector<GEOLIB::Point*> all_stations;
	fetchGeometries (geo, selected_geometries, all_points, all_polylines, all_stations);

	// search bounding polygon
	size_t bp_idx (0); // bounding polygon index
	GEOLIB::Polygon* bounding_polygon (getBoundingPolygon(all_polylines, bp_idx));
	if (!bounding_polygon) return;

	// *** QuadTree - determining bounding box
#ifndef NDEBUG
	std::cout << "computing axis aligned bounding box for quadtree ... " << std::flush;
#endif
	// determine axis aligned bounding box
	GEOLIB::Point ll(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0);
	GEOLIB::Point ur(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), 0);
	for (size_t k(0); k<all_points.size(); k++) {
		if ((*(all_points[k]))[0] < ll[0]) ll[0] = (*(all_points[k]))[0];
		if ((*(all_points[k]))[1] < ll[1]) ll[1] = (*(all_points[k]))[1];
		if ((*(all_points[k]))[0] > ur[0]) ur[0] = (*(all_points[k]))[0];
		if ((*(all_points[k]))[1] > ur[1]) ur[1] = (*(all_points[k]))[1];
	}
#ifndef NDEBUG
	std::cout << "ok" << std::endl;
#endif
	// *** QuadTree - create object
#ifndef NDEBUG
	std::cout << "creating quadtree ... " << std::flush;
#endif
	GEOLIB::QuadTree<GEOLIB::Point> quad_tree (ll, ur, number_of_point_per_quadtree_node);
	std::cout << "ok" << std::endl;

	// *** QuadTree - insert points
#ifndef NDEBUG
	std::cout << "inserting " << all_points.size() << " points into quadtree ... " << std::flush;
#endif
	for (size_t k(0); k < all_points.size(); k++) {
		quad_tree.addPoint (all_points[k]);
	}
#ifndef NDEBUG
	std::cout << "ok" << std::endl;
#endif

	// *** QuadTree - insert stations
#ifndef NDEBUG
	std::cout << "inserting " << all_stations.size() << " stations into quadtree ... " << std::flush;
#endif
	for (size_t k(0); k<all_stations.size(); k++) {
		quad_tree.addPoint (all_stations[k]);
	}
#ifndef NDEBUG
	std::cout << "ok" << std::endl;
#endif

	// *** QuadTree - balance
#ifndef NDEBUG
	std::cout << "balancing quadtree ... " << std::flush;
#endif
	quad_tree.balance ();
#ifndef NDEBUG
	std::cout << "ok" << std::endl;
#endif

	// *** GMSH - write all non-station points
	const size_t n (all_points.size());
	for (size_t k(0); k<n; k++) {
		if (bounding_polygon->isPntInPolygon (*(all_points[k]))) {
			GEOLIB::Point ll, ur;
			quad_tree.getLeaf (*(all_points[k]), ll, ur);
			double mesh_density (mesh_density_scaling * (ur[0]-ll[0]));
			_out << "Point(" << _n_pnt_offset + k << ") = {" << (*(all_points[k]))[0] << ","
				<< (*(all_points[k]))[1] << "," << (*(all_points[k]))[2]
				<< "," << mesh_density
				<< "};" << std::endl;
		}
	}
	_n_pnt_offset += n;

	// write the bounding polygon
	writeBoundingPolygon ( bounding_polygon );

	// write all other polylines as constraints
	const size_t n_polylines (all_polylines.size());
	for (size_t k(0); k<n_polylines; k++) {
		if (k != bp_idx) {
			bool begin_line_pnt_inside_polygon (true);
			bool end_line_pnt_inside_polygon (true);

			size_t s (all_polylines[k]->getNumberOfPoints());

			// write line segments (= Line) of the polyline
			for (size_t j(0); j<s-1; j++) {
				// check if line segment is contained in bounding polygon
				bool line_seg_is_already_used (GEOLIB::containsEdge (*(dynamic_cast<GEOLIB::Polyline*>(bounding_polygon)), (all_polylines[k])->getPointID(j), (all_polylines[k])->getPointID(j+1)));
				// check if line segment is contained in a previous polyline
				for (size_t i(0); i<k && !line_seg_is_already_used; i++) {
					line_seg_is_already_used = GEOLIB::containsEdge (*(all_polylines[i]), (all_polylines[k])->getPointID(j), (all_polylines[k])->getPointID(j+1));
				}

				if (!line_seg_is_already_used) {
					// check if first point of polyline is inside bounding polygon
					if (j==0) {
						begin_line_pnt_inside_polygon = bounding_polygon->isPntInPolygon (*(all_polylines[k])->getPoint(j));
					}
					// check if end point of the line is inside bounding polygon
					end_line_pnt_inside_polygon = bounding_polygon->isPntInPolygon (*(all_polylines[k])->getPoint(j+1));

					if (begin_line_pnt_inside_polygon && end_line_pnt_inside_polygon) {
						_out << "Line(" << _n_lines+j << ") = {" << (all_polylines[k])->getPointID(j) << ","
								<< (all_polylines[k])->getPointID(j+1) << "};" << std::endl;
						// write line as constraint
						_out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
					} else {
						if (begin_line_pnt_inside_polygon && !end_line_pnt_inside_polygon) {
							// create new point
							GEOLIB::Point *s (bounding_polygon->getIntersectionPointPolygonLine	(*(all_polylines[k])->getPoint(j), *(all_polylines[k])->getPoint(j+1)));
							if (s != NULL) {
								// write new point as gmsh geo point with mesh density from existing point
								GEOLIB::Point ll, ur;
								quad_tree.getLeaf (*(all_polylines[k])->getPoint(j), ll, ur);
								double mesh_density (0.3*(ur[0]-ll[0])); // scaling with 0.3 - do not know if this is a good value
								_out << "Point(" << _n_pnt_offset << ") = {" << (*s)[0] << ","
									<< (*s)[1] << "," << (*s)[2] << "," << mesh_density
									<< "}; // new end point of polyline " << k << std::endl;
								// write line
								_out << "Line(" << _n_lines+j << ") = {" << (all_polylines[k])->getPointID(j) << ","
									<< _n_pnt_offset << "};" << std::endl;
								// write line as constraint
								_out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
								_n_pnt_offset++;
								delete s;
							}
						}
						if (!begin_line_pnt_inside_polygon && end_line_pnt_inside_polygon) {
							// create new point
							GEOLIB::Point *s (bounding_polygon->getIntersectionPointPolygonLine(*(all_polylines[k])->getPoint(j), *(all_polylines[k])->getPoint(j+1)));
							if (s != NULL) {
								// write new point as gmsh geo point with mesh density from existing point
								GEOLIB::Point ll, ur;
								quad_tree.getLeaf (*(all_polylines[k])->getPoint(j+1), ll, ur);
								double mesh_density (0.3*(ur[0]-ll[0])); // scaling with 0.3 - do not know if this is a good value
								_out << "Point(" << _n_pnt_offset << ") = {" << (*s)[0] << ","
									<< (*s)[1] << "," << (*s)[2] << "," << mesh_density
									<< "}; // new end point of polyline " << k << std::endl;
								// write line
								_out << "Line(" << _n_lines+j << ") = {" << _n_pnt_offset << "," << (all_polylines[k])->getPointID(j+1)
									<< "};" << std::endl;
								// write line as constraint
								_out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
								_n_pnt_offset++;
								delete s;
							}
						}
					}
					begin_line_pnt_inside_polygon = end_line_pnt_inside_polygon;
				}
			}
			// update line counter
			_n_lines += s;
		}
	}

	// write stations as constraints
	_out << "// Stations" << std::endl;
	const size_t n_stations (all_stations.size());
	for (size_t k(0); k<n_stations; k++) {
		GEOLIB::Point ll, ur;
		quad_tree.getLeaf (*(all_stations[k]), ll, ur);
		double mesh_density (mesh_density_scaling_station_pnts * (ur[0]-ll[0]));
		_out << "Point(" << _n_pnt_offset+k << ") = {" << (*(all_stations[k]))[0] << ","
			<< (*(all_stations[k]))[1] << "," << (*(all_stations[k]))[2] << "," << mesh_density
			<< "};" << std::endl;

		_out << "Point {" << _n_pnt_offset+k << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
	}
	_n_pnt_offset += n_stations;

	// write Steiner points
	std::list<GEOLIB::QuadTree<GEOLIB::Point>*> leaf_list;
	quad_tree.getLeafs (leaf_list);
	_out << "// write Steiner points" << std::endl;
	for (std::list<GEOLIB::QuadTree<GEOLIB::Point>*>::const_iterator it (leaf_list.begin());
		it != leaf_list.end(); it++) {
		if ((*it)->getPoints().empty()) {
			// compute point from square
			GEOLIB::Point ll, rr;
			(*it)->getSquarePoints (ll, rr);
			GEOLIB::Point mid_point (0.5*(rr[0]+ll[0]), 0.5*(rr[1]+ll[1]), 0.5*(rr[2]+ll[2]));
			if (bounding_polygon->isPntInPolygon (mid_point)) {
				_out << "Point(" << _n_pnt_offset << ") = {" << mid_point[0] << ","
							<< mid_point[1] << "," << mid_point[2]
							<< "," << 0.5*(rr[0]-ll[0])
							<< "};" << std::endl;
				_out << "Point {" << _n_pnt_offset << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
				_n_pnt_offset++;
			}

		}
	}
	std::cout << "ok" << std::endl;

//	// go through all geometric data sets and get polygons
//	std::vector<GEOLIB::Polygon*> polygon_vec;
//	for (std::vector<std::string>::const_iterator it (geo_names.begin());
//		it != geo_names.end(); it++) {
//		const std::vector<GEOLIB::Polyline*> *polylines (geo.getPolylineVec (*it));
//		for (std::vector<GEOLIB::Polyline*>::const_iterator it (polylines->begin()); it != polylines->end(); it++) {
//			if ((*it)->isClosed ()) {
//				polygon_vec.push_back (new GEOLIB::Polygon (*(*it)));
//			}
//		}
//	}
//
//	// we assume that all polygons are simple polygons
//	// forest consist of (hierarchy) trees
//	std::list<GEOLIB::SimplePolygonTree*> polygon_forest;
//	// create polygon forest
//	for (std::vector<GEOLIB::Polygon*>::iterator
//			polygon_it(polygon_vec.begin()); polygon_it != polygon_vec.end(); polygon_it++) {
//		// get the list and insert the elements as SimplePolygonTree items into the forest
//		const std::list<GEOLIB::Polygon*> simple_polygon_list(
//				(*polygon_it)->getListOfSimplePolygons());
//		for (std::list<GEOLIB::Polygon*>::const_iterator simple_polygon_it(
//				simple_polygon_list.begin()); simple_polygon_it
//				!= simple_polygon_list.end(); simple_polygon_it++) {
//			GEOLIB::SimplePolygonTree *spt(new GEOLIB::SimplePolygonTree(*simple_polygon_it));
//			polygon_forest.push_back(spt);
//		}
//	}
//
//	// create the hierarchy
//	GEOLIB::createPolygonTree(polygon_forest);
//	std::cout << "\"Polygon forest\" consists of "
//			<< polygon_forest.size() << " trees" << std::endl;
//
//	// *** insert additional Points (for instance Stations) and Polylines (Wadis, rivers, ...)
//	// Stations
//	for (std::vector<std::string>::const_iterator it (geo_names.begin());
//		it != geo_names.end(); it++) {
//		const std::vector<GEOLIB::Point*> *stations (geo.getStationVec (*it));
//		if (stations) {
//			// go through all Stations
//			for (std::vector<GEOLIB::Point*>::const_iterator it (stations->begin()); it != stations->end(); it++) {
//				bool nfound (true);
//				// go through all top level polygons / SimplePolygonTrees
//				for (std::list<GEOLIB::SimplePolygonTree*>::iterator polygon_it (polygon_forest.begin());
//					polygon_it != polygon_forest.end() && nfound; polygon_it++) {
//					if ((*polygon_it)->isGeoObjInside (*it)) {
//						(*polygon_it)->insertGeoObj (*it);
//						nfound = false;
//					}
//				}
//			}
//		}
//	}
//	// Polylines
//	for (std::vector<std::string>::const_iterator it (geo_names.begin());
//		it != geo_names.end(); it++) {
//		const std::vector<GEOLIB::Polyline*> *polylines (geo.getPolylineVec (*it));
//		for (std::vector<GEOLIB::Polyline*>::const_iterator it (polylines->begin()); it != polylines->end(); it++) {
//			if (! (*it)->isClosed ()) {
//				bool nfound (true);
//				// go through all top level polygons / SimplePolygonTrees
//				for (std::list<GEOLIB::SimplePolygonTree*>::iterator polygon_it (polygon_forest.begin());
//					polygon_it != polygon_forest.end() && nfound; polygon_it++) {
//					if ((*polygon_it)->isGeoObjInside (*it)) {
//						(*polygon_it)->insertGeoObj (*it);
//						nfound = false;
//					}
//				}
//			}
//		}
//	}

}


void GMSHInterface::writeAllDataToGMSHInputFile (GEOLIB::GEOObjects& geo,
		std::vector<std::string> const & selected_geometries,
		double mesh_density)
{
	// check file stream
	if (! _out) return;

	std::cout << "GMSHInterface::writeGMSHInputFile non adaptive" << std::endl;

	std::vector<GEOLIB::Point*> all_points;
	std::vector<GEOLIB::Polyline*> all_polylines;
	std::vector<GEOLIB::Point*> all_stations;
	fetchGeometries (geo, selected_geometries, all_points, all_polylines, all_stations);

	// search bounding polygon
	size_t bp_idx (0); // bounding polygon index
	GEOLIB::Polygon* bounding_polygon (getBoundingPolygon(all_polylines, bp_idx));
	if (!bounding_polygon) return;

	// *** GMSH - write all non-station points
	const size_t n (all_points.size());
	for (size_t k(0); k<n; k++) {
		if (bounding_polygon->isPntInPolygon (*(all_points[k]))) {
			_out << "Point(" << _n_pnt_offset + k << ") = {" << (*(all_points[k]))[0] << ","
				<< (*(all_points[k]))[1] << "," << (*(all_points[k]))[2]
				<< "," << mesh_density
				<< "};" << std::endl;
		}
	}
	_n_pnt_offset += n;

	std::cout << "write bounding polygon ... " << std::flush;
	// write bounding polygon
	writeBoundingPolygon (bounding_polygon);

	// write all other polylines as constraints
	const size_t n_polylines (all_polylines.size());
	for (size_t k(0); k<n_polylines; k++) {
		if (k != bp_idx) {
			bool begin_line_pnt_inside_polygon (true);
			bool end_line_pnt_inside_polygon (true);

			size_t s (all_polylines[k]->getNumberOfPoints());

			// write line segments (= Line) of the polyline
			for (size_t j(0); j<s-1; j++) {
				// check if line segment is contained in bounding polygon
				bool line_seg_is_already_used (GEOLIB::containsEdge (*(dynamic_cast<GEOLIB::Polyline*>(bounding_polygon)), (all_polylines[k])->getPointID(j), (all_polylines[k])->getPointID(j+1)));
				// check if line segment is contained in a previous polyline
				for (size_t i(0); i<k && !line_seg_is_already_used; i++) {
					line_seg_is_already_used = GEOLIB::containsEdge (*(all_polylines[i]), (all_polylines[k])->getPointID(j), (all_polylines[k])->getPointID(j+1));
				}

				if (!line_seg_is_already_used) {
					// check if first point of polyline is inside bounding polygon
					if (j==0) {
						begin_line_pnt_inside_polygon = bounding_polygon->isPntInPolygon (*(all_polylines[k])->getPoint(j));
					}
					// check if end point of the line is inside bounding polygon
					end_line_pnt_inside_polygon = bounding_polygon->isPntInPolygon (*(all_polylines[k])->getPoint(j+1));

					if (begin_line_pnt_inside_polygon && end_line_pnt_inside_polygon) {
						_out << "Line(" << _n_lines+j << ") = {" << (all_polylines[k])->getPointID(j) << ","
								<< (all_polylines[k])->getPointID(j+1) << "};" << std::endl;
						// write line as constraint
						_out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
					} else {
						if (begin_line_pnt_inside_polygon && !end_line_pnt_inside_polygon) {
							// create new point
							GEOLIB::Point *s (bounding_polygon->getIntersectionPointPolygonLine	(*(all_polylines[k])->getPoint(j), *(all_polylines[k])->getPoint(j+1)));
							if (s != NULL) {
								_out << "Point(" << _n_pnt_offset << ") = {" << (*s)[0] << ","
									<< (*s)[1] << "," << (*s)[2] << "," << mesh_density
									<< "}; // new end point of polyline " << k << std::endl;
								// write line
								_out << "Line(" << _n_lines+j << ") = {" << (all_polylines[k])->getPointID(j) << ","
									<< _n_pnt_offset << "};" << std::endl;
								// write line as constraint
								_out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
								_n_pnt_offset++;
								delete s;
							}
						}
						if (!begin_line_pnt_inside_polygon && end_line_pnt_inside_polygon) {
							// create new point
							GEOLIB::Point *s (bounding_polygon->getIntersectionPointPolygonLine(*(all_polylines[k])->getPoint(j), *(all_polylines[k])->getPoint(j+1)));
							if (s != NULL) {
								_out << "Point(" << _n_pnt_offset << ") = {" << (*s)[0] << ","
									<< (*s)[1] << "," << (*s)[2] << "," << mesh_density
									<< "}; // new end point of polyline " << k << std::endl;
								// write line
								_out << "Line(" << _n_lines+j << ") = {" << _n_pnt_offset << "," << (all_polylines[k])->getPointID(j+1)
									<< "};" << std::endl;
								// write line as constraint
								_out << "Line {" << _n_lines + j << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
								_n_pnt_offset++;
								delete s;
							}
						}
					}
					begin_line_pnt_inside_polygon = end_line_pnt_inside_polygon;
				}
			}
			// update line counter
			_n_lines += s;
		}
	}

	// write stations as constraints
	_out << "// Stations" << std::endl;
	const size_t n_stations (all_stations.size());
	for (size_t k(0); k<n_stations; k++) {
		_out << "Point(" << _n_pnt_offset+k << ") = {" << (*(all_stations[k]))[0] << ","
			<< (*(all_stations[k]))[1] << "," << (*(all_stations[k]))[2] << "," << mesh_density
			<< "};" << std::endl;
		_out << "Point {" << _n_pnt_offset+k << "} In Surface {" << _n_plane_sfc-1 << "};" << std::endl;
	}
	_n_pnt_offset += n_stations;

	std::cout << "ok" << std::endl;
}


bool GMSHInterface::isGMSHMeshFile (const std::string& fname)
{
	std::ifstream input (fname.c_str());

	if (!input) {
		std::cerr << "GMSHInterface::isGMSHMeshFile could not open file " << fname << std::endl;
		return false;
	}

	std::string header_first_line;
	input >> header_first_line;
	if (header_first_line.find ("$MeshFormat") != std::string::npos) {
		// read version
		std::string version;
		getline (input, version);
		getline (input, version);
		std::cerr << "found GMSH mesh file version: " << version << std::endl;
		input.close ();
		return true;
	}

	return false;
}

void GMSHInterface::fetchGeometries (GEOLIB::GEOObjects const & geo,
		std::vector<std::string> const & selected_geometries,
		std::vector<GEOLIB::Point*>& all_points,
		std::vector<GEOLIB::Polyline*>& all_polylines,
		std::vector<GEOLIB::Point*>& all_stations) const
{
	// get names of all available data sources except stations
	std::vector<std::string> geo_names;
	// get station names
	std::vector<std::string> geo_station_names;

	for (std::vector<std::string>::const_iterator it (selected_geometries.begin());
		it != selected_geometries.end(); it++) {
		if ((geo.getPointVecObj (*it))->getType() == GEOLIB::PointVec::POINT) {
			geo_names.push_back (*it);
		} else if ((geo.getPointVecObj (*it))->getType() == GEOLIB::PointVec::STATION) {
			geo_station_names.push_back (*it);
		}
	}

	size_t pnt_offset (0);
	// fetch points and polylines and add them to the vectors, add points to the QuadTree
	for (std::vector<std::string>::const_iterator it (geo_names.begin());
			it != geo_names.end(); it++) {
		// get data from geo
#ifndef NDEBUG
		std::cout << "fetch geometrical data for " << *it << " " << std::flush;
#endif
		const std::vector<GEOLIB::Point*> *pnts (geo.getPointVec (*it));
		const std::vector<GEOLIB::Polyline*> *plys (geo.getPolylineVec (*it));
#ifndef NDEBUG
		std::cout << "ok" << std::endl;
#endif

		// insert points into vector all_points
		all_points.insert (all_points.end(), pnts->begin(), pnts->end());

		for (size_t k(0); k<plys->size(); k++) {
			size_t pos (all_polylines.size());
			// insert new polyline
			all_polylines.push_back (new GEOLIB::Polyline (all_points));
			// copy points
			for (size_t j(0); j<(*plys)[k]->getNumberOfPoints(); j++) {
				// set points of polyline
				(all_polylines[pos])->addPoint (pnt_offset + ((*plys)[k])->getPointID(j));
			}
		}
		pnt_offset += pnts->size();
	}

	for (std::vector<std::string>::const_iterator it (geo_station_names.begin());
		it != geo_station_names.end(); it++) {
		// get data from geo
#ifndef NDEBUG
		std::cout << "fetch station data for " << *it << " " << std::flush;
#endif
		const std::vector<GEOLIB::Point*> *pnts (geo.getPointVec (*it));
#ifndef NDEBUG
		std::cout << "ok" << std::endl;
#endif
		// insert points into vector all_stations
		all_stations.insert (all_stations.end(), pnts->begin(), pnts->end());
	}
}

GEOLIB::Polygon* GMSHInterface::getBoundingPolygon (std::vector<GEOLIB::Polyline*> const & all_polylines, size_t& bp_idx) const
{
	GEOLIB::Polygon* bounding_polygon(NULL);
	const size_t n_polylines (all_polylines.size());
	for (size_t k(0); k<n_polylines; k++) {
		if (all_polylines[k]->isClosed ()) { // == Polygon
			if (bounding_polygon) { // we have already a bounding polygon
				if (! bounding_polygon->isPolylineInPolygon (*(all_polylines[k]))) {
					GEOLIB::Polygon* tmp_polygon (new GEOLIB::Polygon (*(all_polylines[k])));
					if (tmp_polygon->isPolylineInPolygon (*bounding_polygon)) {
						// found new bounding polygon
						delete bounding_polygon;
						bounding_polygon = tmp_polygon;
						bp_idx = k;
					} else {
						std::cerr << "INFO: there is no inclusion relation between the polygons " << k << " and " << bp_idx << std::endl;
					}
				}
			} else {
				bounding_polygon = new GEOLIB::Polygon (*(all_polylines[k]));
				bp_idx = k;
			}
		}
	}

	if (! bounding_polygon) {
		std::cerr << "WARNING: GMSHInterface::writeAllDataToGMSHInputFile: did not found bounding polygon - abort writing" << std::endl;
		return NULL;
	}

	return bounding_polygon;
}

void GMSHInterface::writeBoundingPolygon (GEOLIB::Polygon const * const bounding_polygon )
{
	std::cout << "write bounding polygon ... " << std::flush;
	// write bounding polygon
	size_t s (bounding_polygon->getNumberOfPoints());
	// write line segments (= Line) of the polyline
	for (size_t j(0); j<s-1; j++) {
		_out << "Line(" << _n_lines+j << ") = {" <<  bounding_polygon->getPointID(j) << ","
				<< bounding_polygon->getPointID(j+1) << "};" << std::endl;
	}
	// write the line segments contained in the polyline (=Line Loop)
	_out << "Line Loop (" << _n_lines + s - 1 << ") = {";
	for (size_t j(0); j<s-2; j++) {
		_out << _n_lines + j << ",";
	}
	_out << _n_lines + s-2 << "};" << std::endl;
	_n_lines += s;
	// write plane surface
	_out << "Plane Surface (" << _n_plane_sfc << ") = {" << _n_lines-1 << "};" << std::endl;
	_n_plane_sfc++;
	std::cout << "ok" << std::endl;
}

} // end namespace FileIO
