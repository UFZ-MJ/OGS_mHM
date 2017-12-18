/*
 * OGSIOVer4.h
 *
 *  Created on: Jan 14, 2010
 *      Author: TF / KR
 */

#ifndef OGSIOVER4_H_
#define OGSIOVER4_H_

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "Point.h"
#include "Polyline.h"
#include "Surface.h"
#include "StringTools.h"

// forward declaration
namespace GEOLIB {
class GEOObjects;
}

namespace FileIO {

/** I/O - routines for the OGS-4 gli file format */

/** method reads geometric objects from file in gli format */
void readGLIFileV4 (const std::string& fname, GEOLIB::GEOObjects* geo);

void writeGLIFileV4 (const std::string& fname, const std::string& proj_name, const GEOLIB::GEOObjects& geo);

void writeAllDataToGLIFileV4 (const std::string& fname, const GEOLIB::GEOObjects& geo);

} // end namespace

#endif /* OGSIOVER4_H_ */
