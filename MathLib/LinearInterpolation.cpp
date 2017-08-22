/*
 * LinearInterpolation.cpp
 *
 *  Created on: Sep 7, 2010
 *      Author: fischeth
 */

#include "LinearInterpolation.h"
#include "binarySearch.h"

namespace MATHLIB {

LinearInterpolation::LinearInterpolation(const std::vector<double>& supporting_points, const std::vector<double>& values_at_supp_pnts)
	: _supporting_points (supporting_points), _values_at_supp_pnts (values_at_supp_pnts)
{}

LinearInterpolation::LinearInterpolation(const std::vector<double>& supporting_points, const std::vector<double>& values_at_supp_pnts, const std::vector<double>& points_to_interpolate, std::vector<double>& values_at_interpol_pnts)
	: _supporting_points (supporting_points), _values_at_supp_pnts (values_at_supp_pnts)
{
	values_at_interpol_pnts.clear();
	for (size_t k(0); k<points_to_interpolate.size(); k++)
		values_at_interpol_pnts.push_back (this->getValue (points_to_interpolate[k]));
}

LinearInterpolation::~LinearInterpolation()
{}

double LinearInterpolation::getValue ( double pnt_to_interpolate )
{
	if (pnt_to_interpolate < _supporting_points[0] || _supporting_points[_supporting_points.size()-1] < pnt_to_interpolate)
		return std::numeric_limits<double>::min ();

	// search interval
	size_t interval_idx (searchElement (pnt_to_interpolate, 0, _supporting_points.size(), _supporting_points));

	// compute linear interpolation polynom: y = m * x + n
	double m ((_values_at_supp_pnts[interval_idx+1] - _values_at_supp_pnts[interval_idx]) / (_supporting_points[interval_idx+1] - _supporting_points[interval_idx]));
	double n (m * _supporting_points[interval_idx+1] - _values_at_supp_pnts[interval_idx+1]);

	return m * pnt_to_interpolate + n;
}

} // end MATHLIB
