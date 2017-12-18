/*
 * DivideAndConquerClosestPair.h
 *
 *  Created on: Jan 25, 2011
 *      Author: TF
 */

#ifndef DIVIDEANDCONQUERCLOSESTPAIR_H_
#define DIVIDEANDCONQUERCLOSESTPAIR_H_

#include "ClosestPair.h"
#include "PointWithID.h"

namespace GEOLIB {

class DivideAndConquerClosestPair: public GEOLIB::ClosestPair {
public:
	DivideAndConquerClosestPair(std::vector<GEOLIB::Point*> const & pnts, size_t& id0, size_t& id1);

private:
	double closestPair (size_t s, size_t e, size_t id0, size_t id1);
	std::vector<GEOLIB::PointWithID> _pnts_with_ids_sorted_by_x;
	std::vector<GEOLIB::PointWithID> _pnts_with_ids_sorted_by_y;
};

}

#endif /* DIVIDEANDCONQUERCLOSESTPAIR_H_ */
