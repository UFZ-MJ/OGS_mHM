/**
 * \file testGeo.cpp
 * 2011-01-31 TF Initial implementation
 *
 * Tests for the GEO directory
 */

// ** INCLUDES **
#include <gtest/gtest.h>

#include <string>
#include <vector>

// GEOLIB
#include "GEOObjects.h"
#include "Point.h"
#include "Polyline.h"
#include "Polygon.h"

// FileIO
#include "OGSIOVer4.h"

#include "Configure.h"

TEST(GEO, PointInPolygon)
{
	GEOLIB::GEOObjects* _geo (new GEOLIB::GEOObjects);
	std::string fname(SOURCEPATH);
	fname += "/tests/data/GEO/TestDataPointInPolygon.gli";
	FileIO::readGLIFileV4(fname, _geo);
	const std::vector<GEOLIB::Polyline*>* plys(_geo->getPolylineVec(fname));
	std::vector<GEOLIB::Point*>* pnts_in_polygon(new std::vector<GEOLIB::Point*>);
	std::vector<GEOLIB::Point*>* pnts_outside_of_polygon(new std::vector<GEOLIB::Point*>);
	GEOLIB::Polygon polygon(*((*plys)[0]));
	std::cout << "creating test points ..." << std::flush;
	std::vector<GEOLIB::Point*> pnts;
	for (size_t j(0); j < 320; j++) {
		for (size_t k(0); k < 920; k++) {
			pnts.push_back(new GEOLIB::Point(0.9 + k / 100.0, -1.1 + j / 100.0, 0.0));
		}
	}
	std::cout << pnts.size() << " created" << std::endl;

	ASSERT_EQ (static_cast<int>(pnts.size()), 294400);

	std::cout << "testing points ..." << std::endl;
	const size_t size(pnts.size());
	for (size_t k(0); k < size; k++) {
		if (polygon.isPntInPolygon(*(pnts[k]))) {
			pnts_in_polygon->push_back(pnts[k]);
		} else {
			pnts_outside_of_polygon->push_back(pnts[k]);
		}
	}

	ASSERT_EQ (static_cast<int>(pnts.size()), static_cast<int>(pnts_in_polygon->size() + pnts_outside_of_polygon->size()));
	ASSERT_EQ (static_cast<int>(pnts_in_polygon->size()), 95346);
	ASSERT_EQ (static_cast<int>(pnts_outside_of_polygon->size()), 199054);
}
