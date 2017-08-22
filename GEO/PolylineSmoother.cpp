/*
 * PolylineSmoother.cpp
 *
 *  Created on: Aug 16, 2010
 *      Author: TF
 */

#include "PolylineSmoother.h"

namespace GEOLIB {

PolylineSmoother::PolylineSmoother()
{
}

PolylineSmoother::~PolylineSmoother() {
	// TODO Auto-generated destructor stub
}

Polyline* PolylineSmoother::execute (const Polyline& in_ply)
{
	computePolylineCharacteristics (in_ply);
	Polyline *out_ply (new Polyline (in_ply.getPointsVec()));
	// add points to polyline
	out_ply->addPoint (in_ply.getPointID (0));
	for (size_t k(1); k<in_ply.getNumberOfPoints()-1; k++) {
		if (in_ply.getLength (k) - in_ply.getLength(k-1) >= 0.8 * _average_segment_length)
			out_ply->addPoint (in_ply.getPointID (k));
	}
	return out_ply;
}

void PolylineSmoother::computePolylineCharacteristics (const Polyline& ply)
{
	double sqr_max_dist (std::numeric_limits<double>::min());
	double sqr_min_dist (std::numeric_limits<double>::max());

//	double min_coords[3] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
//	double max_coords[3] = {std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min()};

	size_t n_pnts (ply.getNumberOfPoints());
	size_t n_segs (n_pnts - 1);

	_average_segment_length = ply.getLength (n_segs-1) / n_segs;

	double sigma (0.0);
	for (size_t k(0); k<n_segs-1; k++) {
		double dist (ply.getLength(k+1) - ply.getLength(k));
		sigma += (dist - _average_segment_length) * (dist - _average_segment_length);
		if (dist < sqr_min_dist) sqr_min_dist = dist;
		if (dist > sqr_max_dist) sqr_max_dist = dist;
	}
	sigma /= n_segs;

//	for (size_t k(0); k<n_pnts; k++) {
//		for (size_t j(0); j<3; j++) {
//			if ((*(ply[k]))[j] < min_coords[j]) min_coords[j] = (*(ply[k]))[j];
//			if ((*(ply[k]))[j] > max_coords[j]) max_coords[j] = (*(ply[k]))[j];
//		}
//	}

	_min_segment_length = (sqr_min_dist);
	_max_segment_length = (sqr_max_dist);

	std::cerr << "*** charakteristics for polyline with " << n_pnts << " points ********" << std::endl;
	std::cerr << "minimal segment length:" << _min_segment_length << std::endl;
	std::cerr << "maximal segment length:" << _max_segment_length << std::endl;
	std::cerr << "average segment length:" << _average_segment_length << std::endl;

//	std::cerr << "maximal x coord: " << max_coords[0] << ", minimal x coord: " << min_coords[0] << ", max x dist: " << max_coords[0] - min_coords[0] << std::endl;
//	std::cerr << "maximal y coord: " << max_coords[1] << ", minimal y coord: " << min_coords[1] << ", max y dist: " << max_coords[1] - min_coords[1] << std::endl;
//	std::cerr << "maximal z coord: " << max_coords[2] << ", minimal z coord: " << min_coords[2] << ", max z dist: " << max_coords[2] - min_coords[2] << std::endl;
}

}
