/**
 * \file XMLInterface.h
 * 18/02/2010 KR Initial implementation
 */

#ifndef XMLINTERFACE_H
#define XMLINTERFACE_H

#include "Configure.h"
#include "GEOObjects.h"
#include "FEMCondition.h"

#include <QXmlStreamReader>

class QFile;
class QDomDocument;
class QDomNode;
class QDomElement;

/**
 * \brief Reads and writes GeoObjects to and from XML files.
 * Note: Currently it is possible to store names for objects in the xml-file (the schema DOES support that) but
 * these names will neither be read or written by this interface as this functionality is not yet fully implemented
 * within the new GEOLIB.
 */
class XMLInterface
{
public:
	/**
	 * Constructor
	 * \param geoObjects An GEOObject that into which data will be read or from which data will be written.
	 * \param schemaFile An XML schema file (*.xsd) that defines the structure of a valid data file.
	 */
	XMLInterface(GEOLIB::GEOObjects* geoObjects, const std::string &schemaFile);

	/// As QXMLStreamWriter seems currently unable to include style-file links into xml-files, this method will workaround this issue and include the stylefile link.
	int insertStyleFileDefinition(const QString &fileName) const;

	/// Check if the given xml-file is valid considering the schema-file used in the constructor
	int isValid(const QString &fileName) const;

	/// Sets the schema filename used to check if xml files are valid.
	void setSchema(const std::string &schemaName);

	/// Reads an xml-file containing a GeoSys project.
	/// Project files currently cover only geo- and station data. This will be expanded in the future.
	int readProjectFile(const QString &fileName);

	/// Reads an xml-file containing geometric object definitions into the GEOObjects used in the contructor
	int readGLIFile(const QString &fileName);

	/// Reads an xml-file containing station object definitions into the GEOObjects used in the contructor
	int readSTNFile(const QString &fileName);

	/// Reads an xml-file containing containing FEM Conditions such as Boundary- or Initial Conditions
	int readFEMCondFile(std::vector<FEMCondition*> &conditions, const QString &fileName, const QString &geoName);

	/// Writes a GeoSys project file containing all data that is currently loaded.
	/// Project files currently cover only geo- and station data. This will be expanded in the future.
	int writeProjectFile(const QString &fileName) const;

	/**
	 * Writes geometric data from GEOObjects to an xml-file (using the QString version)
	 * \param file The file into which the data will be written.
	 * \param gliName The name of the GEOOBjects that will be written into the file.
	 */
	void writeGLIFile(const std::string &filename, const std::string &gliName) const
	{
		writeGLIFile (QString::fromStdString(filename), QString::fromStdString(gliName));
	}

	/**
	 * Writes geometric data from GEOObjects to an xml-file
	 * \param file The file into which the data will be written.
	 * \param gliName The name of the GEOOBjects that will be written into the file.
	 */
	void writeGLIFile(const QString &filename, const QString &gliName) const;
	/**
	 * Writes geometric data from GEOObjects to an xml-file
	 * \param filename The filename for the file into which the data will be written.
	 * \param stnName The name of the station vector that will be written into the file.
	 */
	int writeSTNFile(const QString &filename, const QString &stnName) const;
	/**
	 * Writes geometric data from GEOObjects to an xml-file (using QString version)
	 * \param fname The filename for the file into which the data will be written.
	 * \param stn_name The name of the station vector that will be written into the file.
	 */
	int writeSTNFile(std::string const& fname, std::string const &stn_name) const
	{
		return writeSTNFile (QString::fromStdString(fname), QString::fromStdString(stn_name));
	}

	/// Writes borehole-specific data to a station-xml-file.
	void writeBoreholeData(QDomDocument &doc, QDomElement &boreholeTag, GEOLIB::StationBorehole* borehole) const;

private:
	/// Reads GEOLIB::Point-objects from an xml-file
	void readPoints    ( const QDomNode &pointsRoot, std::vector<GEOLIB::Point*> *points, std::map<std::string, size_t> *pnt_names );

	/// Reads GEOLIB::Polyline-objects from an xml-file
	void readPolylines ( const QDomNode &polylinesRoot, std::vector<GEOLIB::Polyline*> *polylines, std::vector<GEOLIB::Point*> *points, const std::vector<size_t> &pnt_id_map, std::map<std::string, size_t> *ply_names );

	/// Reads GEOLIB::Surface-objects from an xml-file
	void readSurfaces  ( const QDomNode &surfacesRoot, std::vector<GEOLIB::Surface*> *surfaces, std::vector<GEOLIB::Point*> *points, const std::vector<size_t> &pnt_id_map, std::map<std::string, size_t> *sfc_names );

	/// Reads GEOLIB::Station- or StationBorehole-objects from an xml-file
	void readStations  ( const QDomNode &stationsRoot, std::vector<GEOLIB::Point*> *stations );

	/// Reads the stratigraphy of a borehole from an xml-file
	void readStratigraphy( const QDomNode &stratRoot, GEOLIB::StationBorehole* borehole );

	/// Read the details of various FEM Conditions from an xml-file
	void readConditions( const QDomNode &condRoot, std::vector<FEMCondition*> &conditions, const QString &geoName, FEMCondition::CondType type);

	/// Checks if a hash for the given data file exists to skip the time-consuming validation part.
	/// If a hash file exists _and_ the hash of the data file is the same as the content of the hash file the validation is skipped
	/// If no hash file exists, the xml-file is validated and a hash file is written if the xml-file was valid.
	bool checkHash(const QString &fileName) const;

	/// Calculates an MD5 hash of the given file.
	QByteArray calcHash(const QString &fileName) const;

	/// Checks if the given file is conform to the given hash.
	bool hashIsGood(const QString &fileName, const QByteArray &hash) const;

	GEOLIB::GEOObjects* _geoObjects;

	std::string _schemaName;
	std::map<size_t, size_t> _idx_map;
};

#endif // XMLINTERFACE_H
