/*
 * EarClippingTriangulation.h
 *
 *  Created on: Feb 23, 2011
 *      Author: TF
 */

#ifndef EARCLIPPINGTRIANGULATION_H_
#define EARCLIPPINGTRIANGULATION_H_

// STL
#include <list>

// GEOLIB
#include "Polygon.h"
#include "Triangle.h"

// MATHLIB
#include "AnalyticalGeometry.h"

namespace MATHLIB {

class EarClippingTriangulation {
public:
	EarClippingTriangulation(const GEOLIB::Polygon* ply, std::list<GEOLIB::Triangle> &triangles, bool rot = true);
	virtual ~EarClippingTriangulation();
private:
	/**
	 * copies the points of the polygon to the vector _pnts
	 */
	inline void copyPolygonPoints (const GEOLIB::Polygon* polygon);
	inline void rotate ();
	inline void ensureCWOrientation ();

	inline bool isEar(size_t v0, size_t v1, size_t v2) const;

	inline void initVertexList ();
	inline void initLists ();
	inline void clipEars ();

	/**
	 * a copy of the polygon points
	 */
	std::vector<GEOLIB::Point*> _pnts;
	std::list<size_t> _vertex_list;
	std::list<size_t> _convex_vertex_list;
	std::list<size_t> _ear_list;

	/**
	 * triangles of the triangulation (maybe in the wrong orientation)
	 */
	std::list<GEOLIB::Triangle> _triangles;
};

} // end namespace MATHLIB

#endif /* EARCLIPPINGTRIANGULATION_H_ */
