/*
 * GeoType.h
 *
 *  Created on: Jun 17, 2010
 *      Author: TF
 */

#ifndef GEOTYPE_H_
#define GEOTYPE_H_

#include <string>

namespace GEOLIB {

/**
 * \ingroup GEOLIB
 */

enum GEOTYPE {
	INVALID = 0,
	POINT,   //!< POINT
	POLYLINE,//!< POLYLINE
	SURFACE, //!< SURFACE
	VOLUME,  //!< VOLUME
	GEODOMAIN//!< GEODOMAIN
};

GEOTYPE convertGeoType (const std::string& geo_type_str);

std::string convertGeoTypeToString (GEOTYPE geo_type);

} // end namespace GEOLIB

#endif /* GEOTYPE_H_ */
