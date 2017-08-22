/*
 * GeoIO.cpp
 *
 *  Created on: Sep 29, 2010
 *      Author: TF
 */

#include <sstream>

#include "FEMIO/GeoIO.h"

// FEM
#include "files0.h" // for GetLineFromFile1

namespace FileIO {

void GeoIO::readGeoInfo (GeoInfo* geo_info, std::ifstream& in_str, std::string& geo_name,
			const GEOLIB::GEOObjects& geo_obj, const std::string& unique_geo_name) const
{
	std::stringstream strstream;
	strstream.str(GetLineFromFile1(&in_str));
	std::string geo_type_name;
	strstream >> geo_type_name;

	if (geo_type_name.find("POINT") != std::string::npos) {
		strstream >> geo_name;
		const GEOLIB::Point *pnt(
				(geo_obj.getPointVecObj(unique_geo_name))->getElementByName(geo_name));
		if (pnt == NULL) {
			std::cerr << "ERROR in GeoIO::GeoIO : point name \"" << geo_name
					<< "\" not found!" << std::endl;
			exit(1);
		}
		geo_info->setGeoType(GEOLIB::POINT);
		geo_info->setGeoObj(pnt);
		strstream.clear();
	}

	if (geo_type_name.find("POLYLINE") != std::string::npos) {
		geo_info->setGeoType(GEOLIB::POLYLINE);
		strstream >> geo_name;
		const GEOLIB::Polyline *ply(
				(geo_obj.getPolylineVecObj(unique_geo_name))->getElementByName(geo_name));
		if (ply == NULL) {
			std::cerr << "error in COutput::Read: polyline name \"" << geo_name
					<< "\" not found!" << std::endl;
			exit(1);
		}
		geo_info->setGeoObj(ply);
		strstream.clear();
	}

	if (geo_type_name.find("SURFACE") != std::string::npos) {
		geo_info->setGeoType(GEOLIB::SURFACE);
		strstream >> geo_name;
		strstream.clear();
	}

	if (geo_type_name.find("VOLUME") != std::string::npos) {
		geo_info->setGeoType(GEOLIB::VOLUME);
		strstream >> geo_name;
		strstream.clear();
	}

	if (geo_type_name.find("DOMAIN") != std::string::npos) {
		geo_info->setGeoType(GEOLIB::GEODOMAIN);
		strstream >> geo_name;
		strstream.clear();
	}
}

}
