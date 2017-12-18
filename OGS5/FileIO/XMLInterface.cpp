/**
 * \file XMLInterface.cpp
 * 18/02/2010 KR Initial implementation
 */

#include "XMLInterface.h"
#include "DateTools.h"

#include <iostream>
#include <QFileInfo>

#include <QFile>
#include <QTextCodec>
#include <QCryptographicHash>
#include <QtXml/QDomDocument>
#if OGS_QT_VERSION > 45
	#include <QtXmlPatterns/QXmlSchema>
	#include <QtXmlPatterns/QXmlSchemaValidator>
#endif // QT_VERSION > 45

#include <QTime>

XMLInterface::XMLInterface(GEOLIB::GEOObjects* geoObjects, const std::string &schemaFile) 
: _geoObjects(geoObjects), _schemaName(schemaFile)
{
}

int XMLInterface::isValid(const QString &fileName) const
{
#if OGS_QT_VERSION > 45
	QXmlSchema schema;
	schema.load( QUrl::fromLocalFile((QString::fromStdString(_schemaName))) );

	if ( schema.isValid() )
	{
		QXmlSchemaValidator validator( schema );
		if ( validator.validate( QUrl::fromLocalFile((fileName))) )
			return 1;
		else
		{
			std::cout << "XMLInterface::isValid() - XML File is invalid (in reference to schema " << _schemaName << ")." << std::endl;
			return 0;
		}
	}
	else
	{
		std::cout << "XMLInterface::isValid() - Schema " << _schemaName << " is invalid." << std::endl;
		return 0;
	}
#else
    Q_UNUSED (fileName);
	std::cout << "XMLInterface: XML schema validation skipped. Qt 4.6 is required for validation." << std::endl;
	return 1;
#endif // QT_VERSION > 45
}

void XMLInterface::setSchema(const std::string &schemaName)
{
	_schemaName = schemaName;
}

int XMLInterface::readProjectFile(const QString &fileName)
{
	QFile* file = new QFile(fileName);
	QFileInfo fi(fileName);
	QString path = (fi.path().length()>3) ? QString(fi.path() + "/") : fi.path();

	QFileInfo si(QString::fromStdString(_schemaName));
	QString schemaPath(si.absolutePath() + "/");

	if (!file->open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "XMLInterface::readProjectFile() - Can't open xml-file." << std::endl;
		delete file;
		return 0;
	}
	if (!checkHash(fileName)) { delete file; return 0; }

	QDomDocument doc("OGS-PROJECT-DOM");
	doc.setContent(file);
	QDomElement docElement = doc.documentElement(); //OpenGeoSysProject
	if (docElement.nodeName().compare("OpenGeoSysProject"))
	{
		std::cout << "XMLInterface::readProjectFile() - Unexpected XML root." << std::endl;
		delete file;
		return 0;
	}

	QDomNodeList fileList = docElement.childNodes();

	for(int i=0; i<fileList.count(); i++)
    {
		if (fileList.at(i).nodeName().compare("geo") == 0)
		{
			this->setSchema(std::string(schemaPath.toStdString() + "OpenGeoSysGLI.xsd"));
			this->readGLIFile(QString(path + fileList.at(i).toElement().text()));
		}
		else if (fileList.at(i).nodeName().compare("stn") == 0)
		{
			this->setSchema(std::string(schemaPath.toStdString() + "OpenGeoSysSTN.xsd"));
			QDomNodeList childList = fileList.at(i).childNodes();
			for(int j=0; j<childList.count(); j++)
				if (childList.at(j).nodeName().compare("file") == 0)
					this->readSTNFile(QString(path + childList.at(j).toElement().text()));
		}
		/*
		else (fileList.at(i).nodeName().compare("msh") == 0)
		{
			GridAdapter msh(fileList.at(i).toElement().text().toStdString());
			// TODO gridadapter to mesh-models
		}
		*/
	}

	return 1;
}


int XMLInterface::readGLIFile(const QString &fileName)
{
	std::string gliName("[NN]");

	QFile* file = new QFile(fileName);
	if (!file->open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "XMLInterface::readGLIFile() - Can't open xml-file " << fileName.toStdString() << "." << std::endl;
		delete file;
		return 0;
	}
	if (!checkHash(fileName)) { delete file; return 0; }

	std::vector<GEOLIB::Point*>    *points    = new std::vector<GEOLIB::Point*>;
	std::vector<GEOLIB::Polyline*> *polylines = new std::vector<GEOLIB::Polyline*>;
	std::vector<GEOLIB::Surface*>  *surfaces  = new std::vector<GEOLIB::Surface*>;

	std::map<std::string, size_t> *pnt_names  = new std::map<std::string, size_t>;
	std::map<std::string, size_t> *ply_names  = new std::map<std::string, size_t>;
	std::map<std::string, size_t> *sfc_names  = new std::map<std::string, size_t>;

	QDomDocument doc("OGS-GLI-DOM");
	doc.setContent(file);
	QDomElement docElement = doc.documentElement(); //OpenGeoSysGLI
	if (docElement.nodeName().compare("OpenGeoSysGLI"))
	{
		std::cout << "XMLInterface::readGLIFile() - Unexpected XML root." << std::endl;
		delete file;
		return 0;
	}

	QDomNodeList geoTypes = docElement.childNodes();

	for (int i = 0; i < geoTypes.count(); i++) {
		if (geoTypes.at(i).nodeName().compare("name") == 0)
			gliName = geoTypes.at(i).toElement().text().toStdString();
		else if (geoTypes.at(i).nodeName().compare("points") == 0)
		{
			readPoints(geoTypes.at(i), points, pnt_names);
			_geoObjects->addPointVec(points, gliName, pnt_names);
		}
		else if (geoTypes.at(i).nodeName().compare("polylines") == 0)
			readPolylines(geoTypes.at(i), polylines, points, _geoObjects->getPointVecObj(gliName)->getIDMap(), ply_names);
		else if (geoTypes.at(i).nodeName().compare("surfaces") == 0)
			readSurfaces(geoTypes.at(i), surfaces, points, _geoObjects->getPointVecObj(gliName)->getIDMap(), sfc_names);
		else
			std::cout << "Unknown XML-Node found in file." << std::endl;
	}
	delete file;

	if (!polylines->empty()) _geoObjects->addPolylineVec(polylines, gliName, ply_names);
	if (!surfaces->empty())  _geoObjects->addSurfaceVec(surfaces, gliName, sfc_names);
	return 1;
}

void XMLInterface::readPoints( const QDomNode &pointsRoot, std::vector<GEOLIB::Point*> *points, std::map<std::string, size_t> *pnt_names )
{
	char* pEnd;
	QDomElement point = pointsRoot.firstChildElement();
	while (!point.isNull())
	{
		if (point.hasAttribute("id") && point.hasAttribute("x") && point.hasAttribute("y"))
		{
			_idx_map.insert (std::pair<size_t,size_t>(strtol((point.attribute("id")).toStdString().c_str(), &pEnd, 10), points->size()));
			double zVal = (point.hasAttribute("z")) ? strtod((point.attribute("z")).toStdString().c_str(), 0) : 0.0;
			GEOLIB::Point* p = new GEOLIB::Point(strtod((point.attribute("x")).toStdString().c_str(), 0),
												 strtod((point.attribute("y")).toStdString().c_str(), 0),
												 zVal);
			if (point.hasAttribute("name")) pnt_names->insert( std::pair<std::string, size_t>(point.attribute("name").toStdString(), points->size()) );
			points->push_back(p);
		}
		else std::cout << "XMLInterface::readPoints() - Attribute missing in <point> tag ..." << std::endl;

		point = point.nextSiblingElement();
	}
	if (pnt_names->empty()) pnt_names = NULL; // if names-map is empty, set it to NULL because it is not needed
}

void XMLInterface::readPolylines( const QDomNode &polylinesRoot, std::vector<GEOLIB::Polyline*> *polylines, std::vector<GEOLIB::Point*> *points, const std::vector<size_t> &pnt_id_map, std::map<std::string, size_t> *ply_names )
{
	size_t idx(0);
	QDomElement polyline = polylinesRoot.firstChildElement();
	while (!polyline.isNull())
	{
		if (polyline.hasAttribute("id"))
		{
			idx = polylines->size();
			polylines->push_back(new GEOLIB::Polyline(*points));

			if (polyline.hasAttribute("name")) ply_names->insert( std::pair<std::string, size_t>(polyline.attribute("name").toStdString(), idx) );

			QDomElement point = polyline.firstChildElement();
			while (!point.isNull())
			{
				(*polylines)[idx]->addPoint(pnt_id_map[_idx_map[atoi(point.text().toStdString().c_str())]]);
				point = point.nextSiblingElement();
			}
		}
		else std::cout << "XMLInterface::readPolylines() - Attribute missing in <polyline> tag ..." << std::endl;

		polyline = polyline.nextSiblingElement();
	}
	if (ply_names->empty()) ply_names = NULL; // if names-map is empty, set it to NULL because it is not needed
}

void XMLInterface::readSurfaces( const QDomNode &surfacesRoot, std::vector<GEOLIB::Surface*> *surfaces, std::vector<GEOLIB::Point*> *points, const std::vector<size_t> &pnt_id_map, std::map<std::string, size_t> *sfc_names )
{
	QDomElement surface = surfacesRoot.firstChildElement();
	while (!surface.isNull())
	{
		if (surface.hasAttribute("id"))
		{
			surfaces->push_back(new GEOLIB::Surface(*points));
			if (surface.hasAttribute("name")) sfc_names->insert( std::pair<std::string, size_t>(surface.attribute("name").toStdString(), surfaces->size()-1) );

			QDomElement element = surface.firstChildElement();
			while (!element.isNull())
			{
				if (element.hasAttribute("p1") && element.hasAttribute("p2") && element.hasAttribute("p3"))
				{
					size_t p1 = pnt_id_map[_idx_map[atoi((element.attribute("p1")).toStdString().c_str())]];
					size_t p2 = pnt_id_map[_idx_map[atoi((element.attribute("p2")).toStdString().c_str())]];
					size_t p3 = pnt_id_map[_idx_map[atoi((element.attribute("p3")).toStdString().c_str())]];
					surfaces->back()->addTriangle(p1,p2,p3);
				}
				else std::cout << "XMLInterface::readSurfaces() - Attribute missing in <element> tag ..." << std::endl;
				element = element.nextSiblingElement();
			}
		}
		else std::cout << "XMLInterface::readSurfaces() - Attribute missing in <surface> tag ..." << std::endl;

		surface = surface.nextSiblingElement();
	}
	if (sfc_names->empty()) sfc_names = NULL; // if names-map is empty, set it to NULL because it is not needed
}

int XMLInterface::readSTNFile(const QString &fileName)
{
	QFile* file = new QFile(fileName);
	if (!file->open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "XMLInterface::readSTNFile() - Can't open xml-file." << std::endl;
		delete file;
		return 0;
	}
	if (!checkHash(fileName)) { delete file; return 0; }

	QDomDocument doc("OGS-STN-DOM");
	doc.setContent(file);
	QDomElement docElement = doc.documentElement(); //root element, used for identifying file-type
	if (docElement.nodeName().compare("OpenGeoSysSTN"))
	{
		std::cout << "XMLInterface::readSTNFile() - Unexpected XML root." << std::endl;
		delete file;
		return 0;
	}

	QDomNodeList lists = docElement.childNodes();
	for (int i=0; i<lists.count(); i++)
    {
		// read all the station lists
		QDomNodeList stationList = lists.at(i).childNodes();
		std::vector<GEOLIB::Point*> *stations = new std::vector<GEOLIB::Point*>;
		std::string stnName("[NN]");

		for (int j=0; j<stationList.count(); j++)
		{
			if (stationList.at(j).nodeName().compare("name") == 0) stnName = stationList.at(j).toElement().text().toStdString();
			else if (stationList.at(j).nodeName().compare("stations") == 0)  readStations(stationList.at(j), stations);
			else if (stationList.at(j).nodeName().compare("boreholes") == 0) readStations(stationList.at(j), stations);
		}

		GEOLIB::Color* color = GEOLIB::getRandomColor();
		if (!stations->empty()) _geoObjects->addStationVec(stations, stnName, color);
		else delete stations;
	}

	delete file;

	return 1;
}

void XMLInterface::readStations( const QDomNode &stationsRoot, std::vector<GEOLIB::Point*> *stations )
{
	QDomElement station = stationsRoot.firstChildElement();
	while (!station.isNull())
	{
		if (station.hasAttribute("id") && station.hasAttribute("x") && station.hasAttribute("y"))
		{
			std::string stationName("[NN]");
			std::string boreholeDate("0000-00-00");
			double boreholeDepth(0.0);

			QDomNodeList stationFeatures = station.childNodes();
			for(int i=0; i<stationFeatures.count(); i++)
		    {
				if (stationFeatures.at(i).nodeName().compare("name") == 0) stationName = stationFeatures.at(i).toElement().text().toStdString();
				/* add other station features here */

				else if (stationFeatures.at(i).nodeName().compare("bdepth") == 0) boreholeDepth = strtod(stationFeatures.at(i).toElement().text().toStdString().c_str(), 0);
				else if (stationFeatures.at(i).nodeName().compare("bdate") == 0)  boreholeDate  = stationFeatures.at(i).toElement().text().toStdString();
				/* add other borehole features here */
			}

			double zVal = (station.hasAttribute("z")) ? strtod((station.attribute("z")).toStdString().c_str(), 0) : 0.0;

			if (station.nodeName().compare("station") == 0)
			{
				GEOLIB::Station* s = new GEOLIB::Station(strtod((station.attribute("x")).toStdString().c_str(), 0),
														 strtod((station.attribute("y")).toStdString().c_str(), 0),
														 zVal, stationName);
				stations->push_back(s);
			}
			else if (station.nodeName().compare("borehole") == 0)
			{
				GEOLIB::StationBorehole* s = GEOLIB::StationBorehole::createStation(stationName,
																	strtod((station.attribute("x")).toStdString().c_str(), 0),
																	strtod((station.attribute("y")).toStdString().c_str(), 0),
																	zVal, boreholeDepth, boreholeDate);
				/* add stratigraphy to the borehole */
				for(int i=0; i<stationFeatures.count(); i++)
					if (stationFeatures.at(i).nodeName().compare("strat") == 0)
						this->readStratigraphy(stationFeatures.at(i), s);

				stations->push_back(s);
			}

		}
		else std::cout << "XMLInterface::readStations() - Attribute missing in <station> tag ..." << std::endl;
		station = station.nextSiblingElement();
	}
}

void XMLInterface::readStratigraphy( const QDomNode &stratRoot, GEOLIB::StationBorehole* borehole )
{
	//borehole->addSoilLayer((*borehole)[0], (*borehole)[1], (*borehole)[2], "");
	QDomElement horizon = stratRoot.firstChildElement();
	while (!horizon.isNull())
	{
		if (horizon.hasAttribute("id") && horizon.hasAttribute("x") && horizon.hasAttribute("y") && horizon.hasAttribute("z"))
		{
			std::string horizonName("[NN]");

			QDomNodeList horizonFeatures = horizon.childNodes();
			for(int i=0; i<horizonFeatures.count(); i++)
		    {
				if (horizonFeatures.at(i).nodeName().compare("name") == 0) horizonName = horizonFeatures.at(i).toElement().text().toStdString();
				/* add other horizon features here */
			}
			borehole->addSoilLayer(strtod((horizon.attribute("x")).toStdString().c_str(), 0),
				                   strtod((horizon.attribute("y")).toStdString().c_str(), 0),
							       strtod((horizon.attribute("z")).toStdString().c_str(), 0),
							       horizonName);
		}
		else std::cout << "XMLInterface::readStratigraphy() - Attribute missing in <horizon> tag ..." << std::endl;
		horizon = horizon.nextSiblingElement();
	}
}

int XMLInterface::readFEMCondFile(std::vector<FEMCondition*> &conditions, const QString &fileName, const QString &geoName)
{
	QFile* file = new QFile(fileName);
	if (!file->open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "XMLInterface::readFEMCondFile() - Can't open xml-file." << std::endl;
		delete file;
		return 0;
	}
	if (!checkHash(fileName)) { delete file; return 0; }

	QDomDocument doc("OGS-Cond-DOM");
	doc.setContent(file);
	QDomElement docElement = doc.documentElement(); //root element, used for identifying file-type
	if (docElement.nodeName().compare("OpenGeoSysCond"))
	{
		std::cout << "XMLInterface::readFEMCondFile() - Unexpected XML root." << std::endl;
		delete file;
		return 0;
	}

	//std::vector<FEMCondition*> conditions;
	QDomNodeList lists = docElement.childNodes();
	for (int i=0; i<lists.count(); i++)
    {
		if      (lists.at(i).nodeName().compare("BoundaryConditions") == 0) readConditions(lists.at(i), conditions, geoName, FEMCondition::BOUNDARY_CONDITION);
		else if (lists.at(i).nodeName().compare("InitialConditions") == 0)  readConditions(lists.at(i), conditions, geoName, FEMCondition::INITIAL_CONDITION);
		else if (lists.at(i).nodeName().compare("SourceTerms") == 0)        readConditions(lists.at(i), conditions, geoName, FEMCondition::SOURCE_TERM);
	}
	if (!conditions.empty()) return 1;//do something like _geoObjects->addStationVec(stations, stnName, color);
	else
	{
		std::cout << "XMLInterface::readFEMCondFile() - No FEM Conditions found..." << std::endl;
		return 0;
	}

	delete file;

	return 1;
}

void XMLInterface::readConditions( const QDomNode &listRoot, std::vector<FEMCondition*> &conditions, const QString &geoName, FEMCondition::CondType type)
{
	QDomElement cond = listRoot.firstChildElement();
	while (!cond.isNull())
	{
		FEMCondition* c = new FEMCondition(geoName.toStdString(), type);
		QDomNodeList condProperties = cond.childNodes();
		for (int i=0; i<condProperties.count(); i++)
		{
			if (condProperties.at(i).nodeName().compare("ProcessType") == 0) c->setProcessType(convertProcessType(condProperties.at(i).toElement().text().toStdString()));
			if (condProperties.at(i).nodeName().compare("PrimaryVariable") == 0) c->setProcessPrimaryVariable(convertPrimaryVariable(condProperties.at(i).toElement().text().toStdString()));
			if (condProperties.at(i).nodeName().compare("GeoType") == 0)
			{
				QDomNodeList geoProps = condProperties.at(i).childNodes();
				for (int j=0; j<geoProps.count(); j++)
				{
					if (geoProps.at(j).nodeName().compare("geoObject") == 0) c->setGeoType(GEOLIB::convertGeoType(geoProps.at(j).toElement().text().toStdString()));
					if (geoProps.at(j).nodeName().compare("geoName") == 0) c->setGeoName(geoProps.at(j).toElement().text().toStdString());
				}
			}
			if (condProperties.at(i).nodeName().compare("DisType") == 0)
			{
				QDomNodeList distProps = condProperties.at(i).childNodes();
				for (int j=0; j<distProps.count(); j++)
				{
					if (distProps.at(j).nodeName().compare("disName") == 0) c->setProcessDistributionType(FiniteElement::convertDisType(distProps.at(j).toElement().text().toStdString()));
					if (distProps.at(j).nodeName().compare("disValue") == 0) c->setDisValue(strtod(distProps.at(j).toElement().text().toStdString().c_str(),0));
				}
			}
		}
		conditions.push_back(c);
		cond = cond.nextSiblingElement();
	}
}

int XMLInterface::writeProjectFile(const QString &fileName) const
{
	std::fstream stream(fileName.toStdString().c_str(), std::ios::out);
	QFileInfo fi(fileName);
	QString path(fi.absolutePath() + "/");
	if (!stream.is_open())
    {
		std::cout << "XMLInterface::writeProjectFile() - Could not open file...\n";
		return 0;
	}

	stream << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";	// xml definition
	stream << "<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysProject.xsl\"?>\n\n"; // stylefile definition

	QDomDocument doc("OGS-PROJECT-DOM");
	QDomElement root = doc.createElement("OpenGeoSysProject");
	root.setAttribute( "xmlns:ogs", "http://www.opengeosys.net" );
    root.setAttribute( "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance" );
    root.setAttribute( "xsi:noNamespaceSchemaLocation", "http://141.65.34.25/OpenGeoSysProject.xsd" );

	doc.appendChild(root);

	// GLI
	std::vector<std::string> geoNames;
	_geoObjects->getGeometryNames(geoNames);
	for (std::vector<std::string>::const_iterator it(geoNames.begin());	it != geoNames.end(); it++)
	{
		// write GLI file
		QString name(QString::fromStdString(*it));
		this->writeGLIFile(QString(path + name + ".gml"), name);

		// write entry in project file
		QDomElement geoTag = doc.createElement("geo");
		root.appendChild(geoTag);
		QDomElement fileNameTag = doc.createElement("file");
		geoTag.appendChild(fileNameTag);
		QDomText fileNameText = doc.createTextNode(QString(name + ".gml"));
		fileNameTag.appendChild(fileNameText);
	}

	// STN
	std::vector<std::string> stnNames;
	_geoObjects->getStationNames(stnNames);
	for (std::vector<std::string>::const_iterator it(stnNames.begin());	it != stnNames.end(); it++)
	{
		// write STN file
		QString name(QString::fromStdString(*it));

		if (this->writeSTNFile(QString(path + name + ".stn"), name))
		{
			// write entry in project file
			QDomElement geoTag = doc.createElement("stn");
			root.appendChild(geoTag);
			QDomElement fileNameTag = doc.createElement("file");
			geoTag.appendChild(fileNameTag);
			QDomText fileNameText = doc.createTextNode(QString(name + ".stn"));
			fileNameTag.appendChild(fileNameText);
		}
		else std::cout << "XMLInterface::writeProjectFile() -  Error writing file: " << name.toStdString() << std::endl;
	}

	std::string xml = doc.toString().toStdString();
	stream << xml;
	stream.close();
	return 1;
}

void XMLInterface::writeGLIFile(const QString &filename, const QString &gliName) const
{
	QFile file(filename);
	file.open( QIODevice::WriteOnly );
	std::cout << "Writing " << filename.toStdString() << " ... ";

	size_t nPoints=0, nPolylines=0, nSurfaces=0;

	QXmlStreamWriter xml(&file);
	xml.setAutoFormatting(true);
	xml.setCodec(QTextCodec::codecForName("ISO-8859-1"));

	xml.writeStartDocument();

	// to-do: insert stylesheet tag
	// <?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>
	// reserves space at the location where the style-file entry will be places later
	xml.writeCharacters("                                                                                \n\n\n");
	xml.writeStartElement("OpenGeoSysGLI");
	xml.writeNamespace("http://www.opengeosys.net", "ogs");
	xml.writeAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
	xml.writeAttribute("xsi:noNamespaceSchemaLocation", "http://141.65.34.25/OpenGeoSysGLI.xsd");

	xml.writeTextElement("name", gliName);

	// POINTS
	xml.writeStartElement("points");

	const GEOLIB::PointVec *pnt_vec (_geoObjects->getPointVecObj(gliName.toStdString()));
	if (pnt_vec)
	{
		const std::vector<GEOLIB::Point*> *points (pnt_vec->getVector());

		nPoints = points->size();
		for (size_t i=0; i<nPoints; i++)
		{
			xml.writeStartElement("point");
			xml.writeAttribute("id", QString::number(i));
			xml.writeAttribute("x", QString::number((*(*points)[i])[0], 'f'));
			xml.writeAttribute("y", QString::number((*(*points)[i])[1], 'f'));
			xml.writeAttribute("z", QString::number((*(*points)[i])[2], 'f'));

			std::string point_name;
			if (pnt_vec->getNameOfElementByID(i, point_name))
				xml.writeAttribute("name", QString::fromStdString(point_name));
			xml.writeEndElement(); //point
		}
		xml.writeEndElement(); //points
	}
	else std::cout << "Point vector empty, no points written to file." << std::endl;

	// POLYLINES
	const GEOLIB::PolylineVec *ply_vec (_geoObjects->getPolylineVecObj(gliName.toStdString()));
	if (ply_vec)
	{
		const std::vector<GEOLIB::Polyline*> *polylines (ply_vec->getVector());

		if (polylines) {
			xml.writeStartElement("polylines");
			nPolylines = polylines->size();
			for (size_t i=0; i<nPolylines; i++)
			{
				xml.writeStartElement("polyline");
				xml.writeAttribute("id", QString::number(i));

				std::string ply_name("");
				if (ply_vec->getNameOfElementByID(i, ply_name))
					xml.writeAttribute("name", QString::fromStdString(ply_name));

				nPoints = (*polylines)[i]->getNumberOfPoints();
				for (size_t j=0; j<nPoints; j++)
				{
					xml.writeTextElement("pnt", QString::number(((*polylines)[i])->getPointID(j)));
				}
				xml.writeEndElement(); //polyline
			}
			xml.writeEndElement(); //polylines
		}
	}
	else std::cout << "Polyline vector empty, no polylines written to file." << std::endl;

	// SURFACES
	const GEOLIB::SurfaceVec *sfc_vec (_geoObjects->getSurfaceVecObj(gliName.toStdString()));
	if (sfc_vec)
	{
		const std::vector<GEOLIB::Surface*> *surfaces (sfc_vec->getVector());

		if (surfaces) {
			xml.writeStartElement("surfaces");
			nSurfaces = surfaces->size();
			for (size_t i=0; i<nSurfaces; i++)
			{
				xml.writeStartElement("surface");
				xml.writeAttribute("id", QString::number(i));
				std::string sfc_name("");
				if (sfc_vec->getNameOfElementByID(i, sfc_name))
					xml.writeAttribute("name", QString::fromStdString(sfc_name));

				// writing the elements compromising the surface
				size_t nElements = ((*surfaces)[i])->getNTriangles();
				for (size_t j=0; j<nElements; j++)
				{
					xml.writeStartElement("element"); //triangle-element
					xml.writeAttribute("p1", QString::number((*(*(*surfaces)[i])[j])[0]));
					xml.writeAttribute("p2", QString::number((*(*(*surfaces)[i])[j])[1]));
					xml.writeAttribute("p3", QString::number((*(*(*surfaces)[i])[j])[2]));
					xml.writeEndElement(); //triangle-element
				}
				xml.writeEndElement(); //surface
			}
			xml.writeEndElement(); //surfaces
		}
	}
	else std::cout << "Surface vector empty, no surfaces written to file." << std::endl;

	xml.writeEndElement(); // OpenGeoSysGLI

	xml.writeEndDocument();

	file.close();

	insertStyleFileDefinition(filename);
	std::cout << "done." << std::endl;
}


int XMLInterface::writeSTNFile(const QString &filename, const QString &stnName) const
{
	std::fstream stream(filename.toStdString().c_str(), std::ios::out);
	if (!stream.is_open())
    {
		std::cout << "XMLInterface::writeSTNFile() - Could not open file...\n";
		return 0;
	}
	stream << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";	// xml definition
	stream << "<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysSTN.xsl\"?>\n\n"; // stylefile definition

	QDomDocument doc("OGS-STN-DOM");
	QDomElement root = doc.createElement("OpenGeoSysSTN");
	root.setAttribute( "xmlns:ogs", "http://www.opengeosys.net" );
    root.setAttribute( "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance" );
    root.setAttribute( "xsi:noNamespaceSchemaLocation", "http://141.65.34.25/OpenGeoSysSTN.xsd" );

	const std::vector<GEOLIB::Point*> *stations (_geoObjects->getStationVec(stnName.toStdString()));
	bool isBorehole = (static_cast<GEOLIB::Station*>((*stations)[0])->type() == GEOLIB::Station::BOREHOLE) ? true : false;

	doc.appendChild(root);
	QDomElement stationListTag = doc.createElement("stationlist");
	root.appendChild(stationListTag);

	QDomElement listNameTag = doc.createElement("name");
	stationListTag.appendChild(listNameTag);
	QDomText stationListNameText = doc.createTextNode(stnName);
	listNameTag.appendChild(stationListNameText);
	QString listType = (isBorehole) ? "boreholes" : "stations";
	QDomElement stationsTag = doc.createElement(listType);
	stationListTag.appendChild(stationsTag);

	size_t nStations(stations->size());
	for (size_t i=0; i<nStations; i++)
	{
		QString stationType =  (isBorehole) ? "borehole" : "station";
		QDomElement stationTag = doc.createElement(stationType);
		stationTag.setAttribute( "id", QString::number(i) );
		stationTag.setAttribute( "x",  QString::number((*(*stations)[i])[0], 'f') );
		stationTag.setAttribute( "y",  QString::number((*(*stations)[i])[1], 'f') );
		stationTag.setAttribute( "z",  QString::number((*(*stations)[i])[2], 'f') );
		stationsTag.appendChild(stationTag);

		QDomElement stationNameTag = doc.createElement("name");
		stationTag.appendChild(stationNameTag);
		QDomText stationNameText = doc.createTextNode(QString::fromStdString(static_cast<GEOLIB::Station*>((*stations)[i])->getName()));
		stationNameTag.appendChild(stationNameText);

		if (isBorehole) writeBoreholeData(doc, stationTag, static_cast<GEOLIB::StationBorehole*>((*stations)[i]));
	}

	std::string xml = doc.toString().toStdString();
	stream << xml;
	stream.close();
	return 1;
}

void XMLInterface::writeBoreholeData(QDomDocument &doc, QDomElement &boreholeTag, GEOLIB::StationBorehole* borehole) const
{
	QDomElement stationDepthTag = doc.createElement("bdepth");
	boreholeTag.appendChild(stationDepthTag);
	QDomText stationDepthText = doc.createTextNode(QString::number(borehole->getDepth(), 'f'));
	stationDepthTag.appendChild(stationDepthText);
	QDomElement stationDateTag = doc.createElement("bdate");
	boreholeTag.appendChild(stationDateTag);
	QDomText stationDateText = doc.createTextNode(QString::fromStdString(date2string(borehole->getDate())));
	stationDateTag.appendChild(stationDateText);

	std::vector<GEOLIB::Point*> profile = borehole->getProfile();
	std::vector<std::string> soilNames = borehole->getSoilNames();
	size_t nHorizons(profile.size());

	if (nHorizons>1)
	{
		QDomElement stratTag = doc.createElement("strat");
		boreholeTag.appendChild(stratTag);

		for (size_t j=1; j<nHorizons; j++)	/// the first entry in the profile vector is just the position of the borehole
		{
			QDomElement horizonTag = doc.createElement("horizon");
			horizonTag.setAttribute( "id", QString::number(j) );
			horizonTag.setAttribute( "x",  QString::number((*profile[j])[0], 'f') );
			horizonTag.setAttribute( "y",  QString::number((*profile[j])[1], 'f') );
			horizonTag.setAttribute( "z",  QString::number((*profile[j])[2], 'f') );
			stratTag.appendChild(horizonTag);
			QDomElement horizonNameTag = doc.createElement("name");
			horizonTag.appendChild(horizonNameTag);
			QDomText horizonNameText = doc.createTextNode(QString::fromStdString(soilNames[j]));
			horizonNameTag.appendChild(horizonNameText);
		}
	}
}

int XMLInterface::insertStyleFileDefinition(const QString &fileName) const
{
	std::string path = fileName.toStdString();
	std::fstream stream(path.c_str());
	std::string line;
	std::string styleDef("\n<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysGLI.xsl\"?>");

	if (!stream.is_open())
    {
		std::cout << "XMLInterface::insertStyleFileDefinition() - Could not open file...\n";
		return 0;
	}

	stream.seekp(43*sizeof(char),std::ios_base::beg);	// go to the correct position in the stream
	stream.write(styleDef.c_str(), 60*sizeof(char));	// write new line with xml-stylesheet definition
	stream.close();
	return 1;
}

bool XMLInterface::checkHash(const QString &fileName) const
{
	QFileInfo fi(fileName);
	QString md5FileName(fileName + ".md5");

	std::ifstream md5( md5FileName.toStdString().c_str() );
	if (md5.is_open())
	{
		char* md5HashStr = new char[16];
		md5.read(md5HashStr, 16);
		QByteArray md5Hash(md5HashStr, 16);
		delete md5HashStr;
		if (hashIsGood(fileName, md5Hash)) return true;
	}

	if (!this->isValid(fileName)) return false;

	std::cout << "File is valid, writing hashfile..." << std::endl;
	QByteArray hash = calcHash(fileName);
	std::ofstream out( md5FileName.toStdString().c_str(), std::ios::out );
	out.write(hash.data(), 16);
	out.close();
	return true;
}

bool XMLInterface::hashIsGood(const QString &fileName, const QByteArray &hash) const
{
	int hashLength = hash.length();
	QByteArray fileHash = calcHash(fileName);
	if (fileHash.length() != hashLength) return false;
	for (int i=0; i<hashLength; i++)
	{
		if (fileHash[i] != hash[i])
		{
			std::cout << "Hashfile does not match data ... checking file ..." << std::endl;
			return false;
		}
	}
	return true;
}

QByteArray XMLInterface::calcHash(const QString &fileName) const
{
	std::ifstream is(fileName.toStdString().c_str(), std::ios::binary );
	is.seekg (0, std::ios::end);
	int length = is.tellg();
	is.seekg (0, std::ios::beg);
	char * buffer = new char [length];
	is.read (buffer,length);
	is.close();

	QByteArray hash = QCryptographicHash::hash(buffer, QCryptographicHash::Md5);
	delete [] buffer;
	return hash;
}
