/*
 * \file AxisAlignedBoundingBox.h
 *
 *  Created on: April 22, 2010
 *      Author: TF
 */

#ifndef AXISALIGNEDBOUNDINGBOX_H_
#define AXISALIGNEDBOUNDINGBOX_H_

#include "Point.h"

namespace GEOLIB {

/**
 *
 * \ingroup GEOLIB
 *
 * \brief Class AABB is a bounding box around a given geometric entity
 * */
class AABB
{
public:
	/**
	 * construction of object, initialization the axis aligned bounding box
	 * */
	AABB ();

	void update (GEOLIB::Point const & pnt);
	/**
	 * update axis aligned bounding box
	 */
	void update (double x, double y, double z);

	/**
	 * update axis aligned bounding box
	 */
	void update (const double *pnt)
	{
		update (pnt[0], pnt[1], pnt[2]);
	}

	/**
	 * check if point is in the axis aligned bounding box
	 * (employing containsPoint (double x, double y, double z))
	 */
	bool containsPoint (GEOLIB::Point const & pnt) const;

	/**
	 * wrapper for GEOLIB::Point
	 */
	bool containsPoint (const double *pnt) const;

	/**
	 * check if point described by its coordinates x, y, z is in
	 * the axis aligned bounding box
	 */
	bool containsPoint (double x, double y, double z) const;

	GEOLIB::Point getMinPoint () const { return _min_pnt; }
	GEOLIB::Point getMaxPoint () const { return _max_pnt; }

private:
	GEOLIB::Point _min_pnt;
	GEOLIB::Point _max_pnt;
};

} // end namespace

#endif /* AXISALIGNEDBOUNDINGBOX_H_ */
