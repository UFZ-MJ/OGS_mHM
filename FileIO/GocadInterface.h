/*
 * IOGocad.h
 *
 *  Created on: Feb 10, 2010
 *      Author: fischeth
 */

#ifndef GocadInterface_H_
#define GocadInterface_H_

#include <string>
#include <vector>
#include <iostream>

#include "Point.h"
#include "Polyline.h"

namespace GEOLIB {
class GEOObjects;
}

namespace FileIO {

class GocadInterface {
public:
	GocadInterface(const std::string &fname, GEOLIB::GEOObjects *geo_obj);
	~GocadInterface();

private:
	void readObjects ( std::istream &in );
	void readHeader ( std::istream &in );
	void readCoordinateSystem ( std::istream &in );
	bool isProperty ( const std::string &line, std::istream &in );
	bool isGeologicInformation ( const std::string &line );
	void readTSurfData ( std::istream &in );
//	void readPLineData ( std::istream &in );

	std::string _fname;
	std::vector<GEOLIB::Point*>* _pnt_vec;
	std::vector<GEOLIB::Polyline*>* _ply_vec;
	static const size_t MAX_COLS_PER_ROW = 256;
};

} // end namespace

#endif /* GocadInterface_H_ */
