/*
 * GeoInfo.h
 *
 *  Created on: Jun 18, 2010
 *      Author: TF
 */

#ifndef GEOINFO_H_
#define GEOINFO_H_

// STL
//#include <cstddef>
//#include <string>
//#include <limits>

// GEO
#include "GeoType.h"
#include "GeoObject.h"

/**
 * \brief GeoInfo stores the type of the geometric entity and
 * the index within the vector the geometric entity is
 * managed. Possible geometric entities are documented in
 * GeoType.h.
 */
class GeoInfo
{
   public:
      /**
       * standard constructor. You need to set the attributes via
       * setGeoType() and setGeoObj()!
       */
      GeoInfo ();
      /**
       * The constructor of a GeoInfo object initializes the
       * attributes of the object.
       * @param geo_type the type of the geometric entity.
       * Possible geometric entities are documented in GeoType.h.
       * @param geo_obj the pointer to an object of class GeoObject
       */
      GeoInfo(GEOLIB::GEOTYPE geo_type, const GEOLIB::GeoObject* geo_obj = NULL);
      /**
       * virtual destructor - destroys the object
       */
      virtual ~GeoInfo();

      /**
       * getter method for the geo type
       * @sa enum GeoType
       * @return the geo type
       */
      GEOLIB::GEOTYPE getGeoType () const;

      /**
       * get the type as a string for log output
       * @return
       */
      std::string getGeoTypeAsString () const;

      /**
       * getter for the pointer to the object
       * @return
       */
      const GEOLIB::GeoObject* getGeoObj () const;

      /**
       * setter for the geo type
       * @sa enum GeoType
       * @param geo_type type of the geometric entity
       */
      void setGeoType (GEOLIB::GEOTYPE geo_type);
      /**
       * setter for the pointer to the GeoObject object
       * @param geo_obj an instance of class GeoObject
       */
      void setGeoObj (const GEOLIB::GeoObject* geo_obj);

   protected:
      /**
       * type of the geometric entity. @sa enum GeoType
       */
      GEOLIB::GEOTYPE _geo_type;
      /**
       * pointer to geometric object (GEOLIB::Point, GEOLIB::Polyline, GEOLIB::Surface, ...)
       */
      const GEOLIB::GeoObject* _geo_obj;
};
#endif                                            /* GEOINFO_H_ */
