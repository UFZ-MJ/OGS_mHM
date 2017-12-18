/*
 * \file AxisAlignedBoundingBox.cpp
 *
 *  Created on: April 22, 2010
 *      Author: TF
 */

#include <limits>
#include <cstddef>
#include "AxisAlignedBoundingBox.h"

namespace GEOLIB {

AABB::AABB ()
{
	for (std::size_t k(0); k<3; k++) {
		_min_pnt[k] = std::numeric_limits<double>::max();
		_max_pnt[k] = std::numeric_limits<double>::min();
	}
}

void AABB::update (GEOLIB::Point const & pnt)
{
	update (pnt[0], pnt[1], pnt[2]);
}

void AABB::update (double x, double y, double z)
{
	if (x < _min_pnt[0]) _min_pnt[0] = x;
	if (_max_pnt[0] < x) _max_pnt[0] = x;
	if (y < _min_pnt[1]) _min_pnt[1] = y;
	if (_max_pnt[1] < y) _max_pnt[1] = y;
	if (z < _min_pnt[2]) _min_pnt[2] = z;
	if (_max_pnt[2] < z) _max_pnt[2] = z;
}

bool AABB::containsPoint (GEOLIB::Point const & pnt) const
{
	return containsPoint (pnt[0], pnt[1], pnt[2]);
}

bool AABB::containsPoint (const double *pnt) const
{
	return containsPoint (pnt[0], pnt[1], pnt[2]);
}

bool AABB::containsPoint (double x, double y, double z) const
{
	if (_min_pnt[0] <= x && x <= _max_pnt[0]) {
		if (_min_pnt[1] <= y && y <= _max_pnt[1]) {
			if (_min_pnt[2] <= z && z <= _max_pnt[2]) {
				return true;
			} else return false;
		} else return false;
	} else return false;
}

}
