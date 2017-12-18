/*
 * Polygon.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: TF
 */

#include <cstdlib> // for exit

#include "Polygon.h"

// MathLib
#include "AnalyticalGeometry.h"
#include "MathTools.h"

// Base
#include "quicksort.h"
#include "swap.h"

namespace GEOLIB {

Polygon::Polygon(const Polyline &ply, bool init) :
	Polyline (ply)
{
	if (init) {
		initialise ();
	}
}

Polygon::~Polygon()
{
	// remove polygons from list
	for (std::list<Polygon*>::iterator it (_simple_polygon_list.begin()); it != _simple_polygon_list.end(); it++) {
		delete *it;
	}
}

bool Polygon::initialise ()
{
	if (this->isClosed()) {
		calculateAxisAlignedBoundingBox();
		ensureCWOrientation();
		return true;
	} else {
		std::cerr << "ERROR in Polygon::initialise() - base polyline is not closed" << std::endl;
		return false;
	}
}

bool Polygon::isPntInPolygon (GEOLIB::Point const & pnt) const
{
	GEOLIB::Point min_aabb_pnt (_aabb.getMinPoint());
	GEOLIB::Point max_aabb_pnt (_aabb.getMaxPoint());

	if (pnt[0] < min_aabb_pnt[0] || max_aabb_pnt[0] < pnt[0] || pnt[1] < min_aabb_pnt[1] || max_aabb_pnt[1] < pnt[1])
		return false;

	size_t n_intersections (0);
	GEOLIB::Point s;

	if (_simple_polygon_list.empty ()) {
		const size_t n_nodes (getNumberOfPoints()-1);
		for (size_t k(0); k<n_nodes; k++) {
			if (((*(getPoint(k)))[1] <= pnt[1] && pnt[1] <= (*(getPoint(k+1)))[1]) ||
					((*(getPoint(k+1)))[1] <= pnt[1] && pnt[1] <= (*(getPoint(k)))[1])) {
				switch (getEdgeType(k, pnt)) {
				case EdgeType::TOUCHING:
					return true;
					break;
				case EdgeType::CROSSING:
					n_intersections++;
					break;
				case EdgeType::INESSENTIAL:
					break;
				default:
					// do nothing
					;
				}
			}
		}
		if (n_intersections%2 == 1) return true;
	} else {
		for (std::list<Polygon*>::const_iterator it (_simple_polygon_list.begin());
			it != _simple_polygon_list.end(); ++it) {
			if ((*it)->isPntInPolygon (pnt)) return true;
		}
		return false;
	}
	return false;
}

bool Polygon::isPntInPolygon (double x, double y, double z) const
{
	const GEOLIB::Point pnt(x,y,z);
	return isPntInPolygon (pnt);
}

bool Polygon::isPolylineInPolygon (const Polyline& ply) const
{
	size_t ply_size (ply.getNumberOfPoints()), cnt (0);
	for (size_t k(0); k<ply_size; k++) {
		if (isPntInPolygon (*(ply[k]))) {
			cnt++;
		}
	}
	if (cnt == ply_size)
		return true;
	return false;
}

GEOLIB::Point* Polygon::getIntersectionPointPolygonLine (GEOLIB::Point const & a, GEOLIB::Point const & b) const
{
	GEOLIB::Point* s (new GEOLIB::Point (0,0,0));

	if (_simple_polygon_list.empty ()) {
		const size_t n_nodes (getNumberOfPoints()-1);
		for (size_t k(0); k<n_nodes; k++) {
			if (MATHLIB::lineSegmentIntersect (*(getPoint(k)), *(getPoint(k+1)), a, b, *s)) {
				return s;
			}
		}
	} else {
		for (std::list<Polygon*>::const_iterator it (_simple_polygon_list.begin());
			it != _simple_polygon_list.end(); ++it) {
			const Polygon* polygon (*it);
			const size_t n_nodes_simple_polygon (polygon->getNumberOfPoints()-1);
			for (size_t k(0); k<n_nodes_simple_polygon; k++) {
				if (MATHLIB::lineSegmentIntersect (*(polygon->getPoint(k)), *(polygon->getPoint(k+1)), a, b, *s)) {
					return s;
				}
			}
		}
	}
	delete s;
	return NULL;
}

const std::list<Polygon*>& Polygon::getListOfSimplePolygons()
{
	if (_simple_polygon_list.empty())
		_simple_polygon_list.push_back (this);
	return _simple_polygon_list;
}

void Polygon::computeListOfSimplePolygons ()
{
	_simple_polygon_list.push_back (this);
	splitPolygonAtPoint (_simple_polygon_list.begin());
	splitPolygonAtIntersection (_simple_polygon_list.begin());

	for (std::list<Polygon*>::iterator it (_simple_polygon_list.begin());
		it != _simple_polygon_list.end(); it++) {
		(*it)->initialise ();
	}
}

EdgeType::value Polygon::getEdgeType (size_t k, GEOLIB::Point const & pnt) const
{
	switch (getLocationOfPoint(k, pnt)) {
	case Location::LEFT: {
		const GEOLIB::Point & v (*(getPoint(k)));
		const GEOLIB::Point & w (*(getPoint(k+1)));
		if (v[1] < pnt[1] && pnt[1] <= w[1]) return EdgeType::CROSSING;
		else return EdgeType::INESSENTIAL;
		break;
	}
	case Location::RIGHT: {
		const GEOLIB::Point & v (*(getPoint(k)));
		const GEOLIB::Point & w (*(getPoint(k+1)));
		if (w[1] < pnt[1] && pnt[1] <= v[1]) return EdgeType::CROSSING;
		else return EdgeType::INESSENTIAL;
		break;
	}
	case Location::BETWEEN:
	case Location::SOURCE:
	case Location::DESTINATION:
		return EdgeType::TOUCHING;
		break;
	default:
		return EdgeType::INESSENTIAL;
	}
}

void Polygon::calculateAxisAlignedBoundingBox ()
{
	size_t n_nodes (getNumberOfPoints());
	for (size_t k(0); k<n_nodes; k++) {
		_aabb.update ((*(getPoint(k))));
	}
}

void Polygon::ensureCWOrientation ()
{
	size_t n_nodes (getNumberOfPoints());
	// get the left most upper point
	size_t min_x_max_y_idx (0);	// for orientation check
	for (size_t k(0); k<n_nodes; k++) {
		if ((*(getPoint(k)))[0] <= (*(getPoint(min_x_max_y_idx)))[0]) {
			if ((*(getPoint(k)))[0] < (*(getPoint(min_x_max_y_idx)))[0]) {
				min_x_max_y_idx = k;
			} else {
				if ((*(getPoint(k)))[1] > (*(getPoint(min_x_max_y_idx)))[1]) {
					min_x_max_y_idx = k;
				}
			}
		}
	}
	// determine orientation
	MATHLIB::Orientation orient;
	if (0 < min_x_max_y_idx && min_x_max_y_idx < n_nodes-1) {
		orient = MATHLIB::getOrientation (
			(*(getPoint(min_x_max_y_idx-1)))[0], (*(getPoint(min_x_max_y_idx-1)))[1],
			(*(getPoint(min_x_max_y_idx)))[0], (*(getPoint(min_x_max_y_idx)))[1],
			(*(getPoint(min_x_max_y_idx+1)))[0], (*(getPoint(min_x_max_y_idx+1)))[1]);
	} else {
		if (0 == min_x_max_y_idx) {
			orient = MATHLIB::getOrientation (getPoint(n_nodes-2), getPoint(0), getPoint(1));
		} else {
			orient = MATHLIB::getOrientation (
				(*(getPoint(n_nodes-2)))[0], (*(getPoint(n_nodes-2)))[1],
				(*(getPoint(n_nodes-1)))[0], (*(getPoint(n_nodes-1)))[1],
				(*(getPoint(0)))[0], (*(getPoint(0)))[1]);
		}
	}
	if (orient == MATHLIB::CCW) {
		// switch orientation
		for (size_t k(0); k<n_nodes/2; k++) {
			BASELIB::swap (_ply_pnt_ids[k], _ply_pnt_ids[n_nodes-1-k]);
		}
	}
}

void Polygon::splitPolygonAtIntersection (std::list<Polygon*>::iterator polygon_it)
{
	size_t idx0 (0), idx1 (0);
	while (polygon_it != _simple_polygon_list.end()) {
		GEOLIB::Point *intersection_pnt (new GEOLIB::Point);
		bool is_simple (!MATHLIB::lineSegmentsIntersect (*polygon_it, idx0, idx1, *intersection_pnt));
		if (!is_simple) {
			// adding intersection point to pnt_vec
			size_t intersection_pnt_id (_ply_pnts.size());
			const_cast<std::vector<Point*>& >(_ply_pnts).push_back (intersection_pnt);

			// split Polygon
			if (idx0 > idx1) BASELIB::swap (idx0, idx1);

			GEOLIB::Polygon* polygon0 (new GEOLIB::Polygon((*polygon_it)->getPointsVec(), false));
			for (size_t k(0); k<=idx0; k++) polygon0->addPoint ((*polygon_it)->getPointID (k));
			polygon0->addPoint (intersection_pnt_id);
			for (size_t k(idx1+1); k<(*polygon_it)->getNumberOfPoints(); k++)
				polygon0->addPoint ((*polygon_it)->getPointID (k));
			if (! polygon0->initialise()) {
				std::cerr << "ERROR in Polygon::splitPolygonAtIntersection polygon0" << std::endl;
				exit (1);
			}

			GEOLIB::Polygon* polygon1 (new GEOLIB::Polygon((*polygon_it)->getPointsVec(), false));
			polygon1->addPoint (intersection_pnt_id);
			for (size_t k(idx0+1); k<=idx1; k++) polygon1->addPoint ((*polygon_it)->getPointID (k));
			polygon1->addPoint (intersection_pnt_id);
			if (! polygon1->initialise()) {
				std::cerr << "ERROR in Polygon::splitPolygonAtIntersection polygon1" << std::endl;
				exit (1);
			}

			// remove original polyline and add two new polylines
			std::list<GEOLIB::Polygon*>::iterator polygon0_it, polygon1_it;
			polygon_it = _simple_polygon_list.erase (polygon_it);
			polygon1_it = _simple_polygon_list.insert (polygon_it, polygon1);
			polygon0_it = _simple_polygon_list.insert (polygon1_it, polygon0);

			splitPolygonAtIntersection (polygon0_it);
			splitPolygonAtIntersection (polygon1_it);
		} else {
			delete intersection_pnt;
		}
		++polygon_it;
	}
}

void Polygon::splitPolygonAtPoint (std::list<GEOLIB::Polygon*>::iterator polygon_it)
{
	size_t n ((*polygon_it)->getNumberOfPoints()-1), idx0 (0), idx1(0);
	size_t *id_vec (new size_t[n]), *perm (new size_t[n]);
	for (size_t k(0); k<n; k++) {
		id_vec[k] = (*polygon_it)->getPointID (k);
		perm[k] = k;
	}

	quicksort (id_vec, 0, n, perm);

	for (size_t k(0); k<n-1; k++) {
		if (id_vec[k] == id_vec[k+1]) {
			idx0 = perm[k];
			idx1 = perm[k+1];
			delete [] perm;
			delete [] id_vec;

			if (idx0 > idx1) BASELIB::swap (idx0, idx1);

			// create two closed polylines
			GEOLIB::Polygon* polygon0 (new GEOLIB::Polygon((*polygon_it)->getPointsVec(), false));
			for (size_t k(0); k<=idx0; k++)
				polygon0->addPoint ((*polygon_it)->getPointID (k));
			for (size_t k(idx1+1); k<(*polygon_it)->getNumberOfPoints(); k++)
				polygon0->addPoint ((*polygon_it)->getPointID (k));
			polygon0->initialise();

			GEOLIB::Polygon* polygon1 (new GEOLIB::Polygon((*polygon_it)->getPointsVec(), false));
			for (size_t k(idx0); k<=idx1; k++)
				polygon1->addPoint ((*polygon_it)->getPointID (k));
			polygon1->initialise();

			// remove original polygon and add two new polygons
			std::list<GEOLIB::Polygon*>::iterator polygon0_it, polygon1_it;
			polygon1_it = _simple_polygon_list.insert (_simple_polygon_list.erase (polygon_it), polygon1);
			polygon0_it = _simple_polygon_list.insert (polygon1_it, polygon0);

			splitPolygonAtPoint (polygon0_it);
			splitPolygonAtPoint (polygon1_it);

			return;
		}
	}
	delete [] perm;
	delete [] id_vec;
}


} // end namespace GEOLIB
