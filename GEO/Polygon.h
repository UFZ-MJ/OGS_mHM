/*
 * Polygon.h
 *
 *  Created on: Jun 21, 2010
 *      Author: TF
 */

#ifndef POLYGON_H_
#define POLYGON_H_

// STL
#include <list>

// GEOLIB
#include "AxisAlignedBoundingBox.h"
#include "Polyline.h"

namespace GEOLIB {

/**
 * \ingroup GEOLIB
 */

/**
 * edge classification
 */
class EdgeType {
	public:
		enum value {
			TOUCHING,  //!< TOUCHING
			CROSSING,  //!< CROSSING
			INESSENTIAL//!< INESSENTIAL
		};
};

/**
 *
 */
class Polygon : public Polyline
{
public:
	/**
	 * constructor checks if the given polyline is closed,
	 * and assures that the orientation is clock wise.
	 * @param ply closed Polyline
	 * @param init if true, check if polyline is closed, calculate bounding box
	 * @return
	 */
	Polygon(const Polyline &ply, bool init = true);
	virtual ~Polygon();

	/**
	 *
	 * @return
	 */
	bool initialise ();

	/**
	 * Method checks if the given point is inside the polygon.
	 * The method requires that the polygon has clock wise orientation.
	 * @param pnt the Point
	 * @return if point is inside the polygon true, else false
	 */
	bool isPntInPolygon (const GEOLIB::Point& pnt) const;
	/**
	 * wrapper for method isPntInPolygon (const GEOLIB::Point&)
	 * @param x x coordinate of point
	 * @param y y coordinate of point
	 * @param z z coordinate of point
	 * @return if point is inside the polygon true, else false
	 */
	bool isPntInPolygon (double x, double y, double z) const;
	bool isPolylineInPolygon (const Polyline& ply) const;
	GEOLIB::Point* getIntersectionPointPolygonLine (GEOLIB::Point const & a, GEOLIB::Point const & b) const;
	void computeListOfSimplePolygons ();
	const std::list<Polygon*>& getListOfSimplePolygons ();

private:
	/**
	 * get the type of edge with respect to the given point (2d method!)
	 * @param k number of line segment
	 * @param pnt point that is edge type computed for
	 * @return a value of enum EdgeType
	 */
	EdgeType::value getEdgeType (size_t k, GEOLIB::Point const & pnt) const;

	void calculateAxisAlignedBoundingBox ();
	void ensureCWOrientation ();

	void splitPolygonAtIntersection (std::list<Polygon*>::iterator polygon_it);
	void splitPolygonAtPoint (std::list<Polygon*>::iterator polygon_it);
	std::list<Polygon*> _simple_polygon_list;
	AABB _aabb;
};

}

#endif /* POLYGON_H_ */
