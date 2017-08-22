/*
 * EarClippingTriangulation.cpp
 *
 *  Created on: Feb 23, 2011
 *      Author: TF
 */

// STL
#include <vector>

// BASELIB
#include "swap.h"
#include "printList.h"
#include "uniqueListInsert.h"

// MATHLIB
#include "EarClippingTriangulation.h"

namespace MATHLIB {

EarClippingTriangulation::EarClippingTriangulation(const GEOLIB::Polygon* polygon,
		std::list<GEOLIB::Triangle> &triangles, bool rot)
{
	copyPolygonPoints (polygon);

	if (rot) {
		rotate ();
		ensureCWOrientation ();
	}

	initVertexList ();
	initLists ();
	clipEars ();

	std::vector<GEOLIB::Point*> const& ref_pnts_vec (polygon->getPointsVec());
	std::list<GEOLIB::Triangle>::const_iterator it (_triangles.begin());
	while (it != _triangles.end()) {
		const size_t i0 (polygon->getPointID ((*it)[0]));
		const size_t i1 (polygon->getPointID ((*it)[1]));
		const size_t i2 (polygon->getPointID ((*it)[2]));
		triangles.push_back (GEOLIB::Triangle (ref_pnts_vec, i0, i1, i2));
		it++;
	}
}

EarClippingTriangulation::~EarClippingTriangulation()
{
	const size_t n_pnts (_pnts.size());
	for (size_t k(0); k<n_pnts; k++) {
		delete _pnts[k];
	}
}

void EarClippingTriangulation::copyPolygonPoints (const GEOLIB::Polygon* polygon)
{
	// copy points - last point is identical to the first
	size_t n_pnts (polygon->getNumberOfPoints()-1);
	for (size_t k(0); k < n_pnts; k++) {
		_pnts.push_back (new GEOLIB::Point (*(polygon->getPoint(k))));
	}
}

void EarClippingTriangulation::rotate ()
{
	// calculate supporting plane
	Vector plane_normal;
	double d;
	// compute the plane normal
	getNewellPlane(_pnts, plane_normal, d);

	double tol (sqrt(std::numeric_limits<double>::min()));
	if (fabs(plane_normal[0]) > tol || fabs(plane_normal[1]) > tol) {
		// rotate copied points into x-y-plane
		rotatePointsToXY(plane_normal, _pnts);
	}
}

void EarClippingTriangulation::ensureCWOrientation ()
{
	size_t n_pnts (_pnts.size());
	// get the left most upper point
	size_t min_x_max_y_idx (0);	// for orientation check
	for (size_t k(0); k<n_pnts; k++) {
		if ((*(_pnts[k]))[0] <= (*(_pnts[min_x_max_y_idx]))[0]) {
			if ((*(_pnts[k]))[0] < (*(_pnts[min_x_max_y_idx]))[0]) {
				min_x_max_y_idx = k;
			} else {
				if ((*(_pnts[k]))[1] > (*(_pnts[min_x_max_y_idx]))[1]) {
					min_x_max_y_idx = k;
				}
			}
		}
	}
	// determine orientation
	MATHLIB::Orientation orient;
	if (0 < min_x_max_y_idx && min_x_max_y_idx < n_pnts-1) {
		orient = MATHLIB::getOrientation (
			_pnts[min_x_max_y_idx-1], _pnts[min_x_max_y_idx], _pnts[min_x_max_y_idx+1]);
	} else {
		if (0 == min_x_max_y_idx) {
			orient = MATHLIB::getOrientation (_pnts[n_pnts-1], _pnts[0], _pnts[1]);
		} else {
			orient = MATHLIB::getOrientation (_pnts[n_pnts-2], _pnts[n_pnts-1], _pnts[0]);
		}
	}
	if (orient == MATHLIB::CCW) {
		// switch orientation
		for (size_t k(0); k<n_pnts/2; k++) {
			BASELIB::swap (_pnts[k], _pnts[n_pnts-1-k]);
		}
	}
}

bool EarClippingTriangulation::isEar(size_t v0, size_t v1, size_t v2) const
{
	for (std::list<size_t>::const_iterator it (_vertex_list.begin ());
		it != _vertex_list.end(); ++it) {
		if (*it != v0 && *it != v1 && *it != v2) {
			if (isPointInTriangle (_pnts[*it], _pnts[v0], _pnts[v1], _pnts[v2])) {
				return false;
			}
		}
	}
	return true;
}

void EarClippingTriangulation::initVertexList ()
{
	size_t n_pnts (_pnts.size());
	for (size_t k(0); k<n_pnts; k++) _vertex_list.push_back (k);
}

void EarClippingTriangulation::initLists ()
{
	// go through points checking ccw, cw or collinear order and identifying ears
	std::list<size_t>::iterator it (_vertex_list.begin()), prev(_vertex_list.end()), next;
	prev--;
	next = it;
	next++;
	MATHLIB::Orientation orientation;
	while (next != _vertex_list.end()) {
		orientation  = getOrientation (_pnts[*prev], _pnts[*it], _pnts[*next]);
		if (orientation == COLLINEAR) {
			it = _vertex_list.erase (it);
			next++;
		} else {
			if (orientation == CW) {
				_convex_vertex_list.push_back (*it);
				if (isEar (*prev, *it, *next))
					_ear_list.push_back (*it);
			}
			prev = it;
			it = next;
			next++;
		}
	}

	next = _vertex_list.begin();
	orientation = getOrientation (_pnts[*prev], _pnts[*it], _pnts[*next]);
	if (orientation == COLLINEAR) {
		it = _vertex_list.erase (it);
	}
	if (orientation == CW) {
		_convex_vertex_list.push_back (*it);
		if (isEar (*prev, *it, *next))
			_ear_list.push_back (*it);
	}
}

void EarClippingTriangulation::clipEars()
{
	std::list<size_t>::iterator it, prev, next;
	// *** clip an ear
	while (_vertex_list.size() > 3) {
		// pop ear from list
		size_t ear = _ear_list.front();
		_ear_list.pop_front();
		// remove ear tip from _convex_vertex_list
		_convex_vertex_list.remove(ear);

		// remove ear from vertex_list, apply changes to _ear_list, _convex_vertex_list
		bool nfound(true);
		it = _vertex_list.begin();
		prev = _vertex_list.end();
		prev--;
		while (nfound && it != _vertex_list.end()) {
			if (*it == ear) {
				nfound = false;
				it = _vertex_list.erase(it); // remove ear tip
				next = it;
				if (next == _vertex_list.end()) {
					next = _vertex_list.begin();
					prev = _vertex_list.end();
					prev--;
				}
				// add triangle
				_triangles.push_back(GEOLIB::Triangle(_pnts, *prev, ear, *next));

				// check the orientation of prevprev, prev, next
				std::list<size_t>::iterator prevprev;
				if (prev == _vertex_list.begin()) {
					prevprev = _vertex_list.end();
				} else {
					prevprev = prev;
				}
				prevprev--;

				// apply changes to _convex_vertex_list and _ear_list looking "backward"
				MATHLIB::Orientation orientation = getOrientation(_pnts[*prevprev], _pnts[*prev],
						_pnts[*next]);
				if (orientation == CW) {
					BASELIB::uniqueListInsert(_convex_vertex_list, *prev);
					// prev is convex
					if (isEar(*prevprev, *prev, *next)) {
						// prev is an ear tip
						BASELIB::uniqueListInsert(_ear_list, *prev);
					} else {
						// if necessary remove prev
						_ear_list.remove(*prev);
					}
				} else {
					// prev is not convex => reflex or collinear
					_convex_vertex_list.remove(*prev);
					_ear_list.remove(*prev);
					if (orientation == COLLINEAR) {
						prev = _vertex_list.erase(prev);
						if (prev == _vertex_list.begin()) {
							prev = _vertex_list.end();
							prev--;
						} else {
							prev--;
						}
					}
				}

				// check the orientation of prev, next, nextnext
				std::list<size_t>::iterator nextnext,
						help_it(_vertex_list.end());
				help_it--;
				if (next == help_it) {
					nextnext = _vertex_list.begin();
				} else {
					nextnext = next;
					nextnext++;
				}

				// apply changes to _convex_vertex_list and _ear_list looking "forward"
				orientation = getOrientation(_pnts[*prev], _pnts[*next],
						_pnts[*nextnext]);
				if (orientation == CW) {
					BASELIB::uniqueListInsert(_convex_vertex_list, *next);
					// next is convex
					if (isEar(*prev, *next, *nextnext)) {
						// next is an ear tip
						BASELIB::uniqueListInsert(_ear_list, *next);
					} else {
						// if necessary remove *next
						_ear_list.remove(*next);
					}
				} else {
					// next is not convex => reflex or collinear
					_convex_vertex_list.remove(*next);
					_ear_list.remove(*next);
					if (orientation == COLLINEAR) {
						next = _vertex_list.erase(next);
						if (next == _vertex_list.end()) next
								= _vertex_list.begin();
					}
				}

			} else {
				prev = it;
				it++;
			}
		}

	}
	// add last triangle
	next = _vertex_list.begin();
	prev = next;
	next++;
	it = next;
	next++;
	if (getOrientation(_pnts[*prev], _pnts[*it], _pnts[*next]) == CCW)
		_triangles.push_back(GEOLIB::Triangle(_pnts, *prev, *it, *next));
	else
		_triangles.push_back(GEOLIB::Triangle(_pnts, *prev, *next, *it));
}

} // end namespace MATHLIB
