/*
 * PolylineSmoother.h
 *
 *  Created on: Aug 16, 2010
 *      Author: TF
 */

#ifndef POLYLINESMOOTHER_H_
#define POLYLINESMOOTHER_H_

#include "Polyline.h"

namespace GEOLIB {

/**
 * \ingroup GEOLIB
 */

class PolylineSmoother {
public:
	PolylineSmoother();
	virtual ~PolylineSmoother();

	Polyline* execute (const Polyline& in_ply);

protected:
	void computePolylineCharacteristics (const Polyline& ply);
	/**
	 * minimal segment length
	 */
	double _min_segment_length;
	/**
	 * maximal segment length
	 */
	double _max_segment_length;
	/**
	 * average segment length
	 */
	double _average_segment_length;
	double _min_threshold_angle;
	double _max_threshold_angle;
};

}

#endif /* POLYLINESMOOTHER_H_ */
