/*
 * DivideAndConquerClosestPair.cpp
 *
 *  Created on: Jan 25, 2011
 *      Author: TF
 */

#include <algorithm>

// MATHLIB
#include "MathTools.h"

#include "Point.h"
#include "DivideAndConquerClosestPair.h"

namespace GEOLIB {

DivideAndConquerClosestPair::DivideAndConquerClosestPair(
		std::vector<GEOLIB::Point*> const & pnts, size_t& id0, size_t& id1) :
	ClosestPair (pnts, id0, id1)
{
	for (size_t k(0); k<_pnts.size(); k++) {
		_pnts_with_ids_sorted_by_x.push_back (GEOLIB::PointWithID ((*(_pnts[k]))[0], (*(_pnts[k]))[1], (*(_pnts[k]))[2], k));
		_pnts_with_ids_sorted_by_y.push_back (GEOLIB::PointWithID ((*(_pnts[k]))[0], (*(_pnts[k]))[1], (*(_pnts[k]))[2], k));
	}

	std::sort (_pnts_with_ids_sorted_by_x.begin(), _pnts_with_ids_sorted_by_x.end(), GEOLIB::lessX);
	std::sort (_pnts_with_ids_sorted_by_y.begin(), _pnts_with_ids_sorted_by_y.end(), GEOLIB::lessY);
}

double DivideAndConquerClosestPair::closestPair (size_t s, size_t e, size_t id0, size_t id1)
{
	if (e == s || e == s+1) {
		id0 = s;
		id1 = e;
		return std::numeric_limits<double>::max();
	}

	// simple case - brute force approach
	if (e-s < 3) {
		double sqr_dist (std::numeric_limits<double>::max());
		for (size_t k(s); k<e; k++) {
			for (size_t j(k+1); j<e; j++) {
				double act_dist (MATHLIB::sqrDist (&(_pnts_with_ids_sorted_by_x[k]), &(_pnts_with_ids_sorted_by_x[j])));
				if (act_dist < sqr_dist) {
					id0 = k;
					id1 = j;
					sqr_dist = act_dist;
				}
			}
		}
		return sqr_dist;
	}

	// divide and conquer
	size_t m ((s+e)/2), lid0, lid1, rid0, rid1;
	// recursion
	double d0 (closestPair (s, m, lid0, lid1));
	double d1 (closestPair (m+1, e, rid0, rid1));

	// merge
	double d (std::min (d0,d1));

}

} // end namespace GEOLIB
