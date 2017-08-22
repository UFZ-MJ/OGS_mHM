/*
 * testPointInTriangle.cpp
 *
 *  Created on: Feb 28, 2011
 *      Author: TF
 */

#include <iostream>

// GEOLIB
#include "Point.h"

// MATHLIB
#include "AnalyticalGeometry.h"

// FileIO
//#include "XMLInterface.h"

int main (int argc, char *argv[])
{
	GEOLIB::Point pnt_a (0.0, 0.0, 0.0);
	GEOLIB::Point pnt_b (1.0, 0.0, 0.0);
	GEOLIB::Point pnt_c (0.0, 1.0, 0.0);

	size_t n_x (10);
	size_t n_y (10);

	for (size_t i(0); i<n_x; i++) {
		for (size_t j(0); j<n_y; j++) {
			GEOLIB::Point pnt_p ((static_cast<double>(i) - n_x/2.0)/n_x, (static_cast<double>(j) - n_y/2.0)/n_y, 0.0);
			std::cout << "(" << pnt_p[0] << ", " << pnt_p[1] << ") is " << std::flush;
			if (MATHLIB::isPointInTriangle (&pnt_p, &pnt_a, &pnt_b, &pnt_c)) {
				std::cout << " in [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]" << std::endl;
			} else {
				std::cout << " *not* in [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]" << std::endl;
			}
		}
		std::cout << std::endl;
	}
}
