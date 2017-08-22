/**
 * \file SHPInterface.cpp
 * 25/01/2010 KR Initial implementation
 */
#ifndef SHPINTERFACE_H
#define SHPINTERFACE_H

#include <string>

//ShapeLib includes
#include "shapefil.h"

#include "GEOObjects.h"

/**
 * \brief Manages the import of ESRI shape files into GEOLIB.
 */
class SHPInterface
{
public:
	/// Connection between ESRI type system for shape files and OGS GEOLIB.
	enum OGSType
	{
		UNDEFINED	= 0,
		POINT		= 1,
		STATION		= 2,
		POLYLINE	= 3,
		POLYGON		= 4
	};

	/// Constructor
	SHPInterface(GEOLIB::GEOObjects* geoObjects) : _geoObjects(geoObjects) {};

	/// Reads the header of the shape file.
	bool readSHPInfo(const std::string &filename, int &shapeType, int &numberOfEntities);

	/// Reads data from the shape file.
	void readSHPFile(const std::string &filename, OGSType choice, std::string listName);

private:
	/// Reads points into a vector of Point objects.
	void readPoints    (const SHPHandle &hSHP, int numberOfElements, std::string listName);

	/// Reads points into a vector of Point objects and marks them as Station.
	void readStations  (const SHPHandle &hSHP, int numberOfElements, std::string listName);

	/// Reads lines into a vector of Polyline objects.
	void readPolylines (const SHPHandle &hSHP, int numberOfElements, std::string listName);

	/// Reads lines into a vector of Polyline and Surface objects.
	void readPolygons  (const SHPHandle &hSHP, int numberOfElements, std::string listName);

	GEOLIB::GEOObjects* _geoObjects;

};

#endif //SHPINTERFACE_H
