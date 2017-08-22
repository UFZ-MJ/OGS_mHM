/*
 * IOGocad.cpp
 *
 *  Created on: Feb 10, 2010
 *      Author: fischeth
 */

#include "GocadInterface.h"
#include "GEOObjects.h"
#include "StringTools.h"

#include <fstream>

using namespace GEOLIB;

namespace FileIO {

GocadInterface::GocadInterface(const std::string &fname, GEOObjects* obj) :
	_fname (fname), _pnt_vec (new std::vector<GEOLIB::Point*>),
	_ply_vec (new std::vector<GEOLIB::Polyline*>)
{
	std::cout << "GocadInterface::GocadInterface open stream from file " << fname << " ... " << std::flush;
	std::ifstream in(fname.c_str());
	if (in) {
		std::cout << "done" << std::endl;
		readObjects (in);
		in.close();
	} else {
		std::cerr << "error opening stream " << std::endl;
	}

	obj->addPointVec (_pnt_vec, _fname);
	obj->addPolylineVec (_ply_vec, _fname);
}

GocadInterface::~GocadInterface()
{}

void GocadInterface::readObjects (std::istream &in)
{
	std::string line;
	while  (in) {
		// read object type
		in >> line;
		if (line.find("GOCAD") != std::string::npos) {
			in >> line;
			if (line.find("TSurf") != std::string::npos) {
				std::cout << "reading GOCAD object TSurf" << std::endl;
				readHeader(in);
				readCoordinateSystem(in);

				char buffer[MAX_COLS_PER_ROW];
				in.getline(buffer, MAX_COLS_PER_ROW);
				line = buffer;
				if (line.find('\r') || line.find('\n')) {
					in.getline(buffer, MAX_COLS_PER_ROW);
					line = buffer;
				}

//				std::cout << "reading optional properties ... " << std::flush;
				// read optional properties
				while (isProperty(line, in) && in) {
					in.getline(buffer, MAX_COLS_PER_ROW);
					line = buffer;
				}
//				std::cout << "done" << std::endl;

//				std::cout << "reading optional geologic information ... "
//						<< std::flush;
				// read optional geologic information
				// gocad file specification sec A.4.1
				while (isGeologicInformation(line) && in) {
					in.getline(buffer, MAX_COLS_PER_ROW);
					line = buffer;
				}
//				std::cout << "done" << std::endl;

				if (line.find("TFACE") != std::string::npos) {
					readTSurfData(in);
				} else
					std::cout << "no points found" << std::endl;

			} else if (line.find("PLine") != std::string::npos) {
				std::cout << "reading GOCAD object PLine" << std::endl;
				readHeader(in);
				readCoordinateSystem(in);
//				readPLineData(in);
			} else
				std::cout << "*** reading GOCAD object " << line
						<< " not implemented" << std::endl;
		}
	}
}

void GocadInterface::readHeader (std::istream &in)
{
	std::string line;
	in >> line; // version
	std::cout << "version " << line << std::endl;
	in >> line;

	// read header of TSurf object
	if (line.find ("HEADER") != std::string::npos) {
//		std::cout << "found header" << std::endl;
//		std::cout << line << std::endl;
		while (in && line.find ("}") == std::string::npos) {
			if (line.find("name:") != std::string::npos) {
				// extract name
//				std::cout << "name of object: " << line.substr (5) << std::endl;
			}
			if (line.find("mesh:") != std::string::npos) {
//				std::cout << "mesh of object: " << line.substr (5) << std::endl;
			}
			if (line.find("cn:") != std::string::npos) {
//				std::cout << "cn of object: " << line.substr (3) << std::endl;
			}
			if (line.find("*solid*color:") != std::string::npos) {
				float col[4];
				col[0] = str2number<float> ((line.substr(13)).c_str());
				in >> col[1] >> col[2] >> col[3];
//				std::cout << "color of object: " << col[0] << "," << col[1] << "," << col[2] << "," << col[3] << std::endl;
			}
			if (line.find("ivolmap:") != std::string::npos) {
//				std::cout << "ivolmap of object: " << line.substr (8) << std::endl;
			}
			if (line.find("imap:") != std::string::npos) {
//				std::cout << "imap of object: " << line.substr (5) << std::endl;
			}
			in >> line;
		}
	}
//	std::cout << "read header done" << std::endl;
}

void GocadInterface::readCoordinateSystem (std::istream &in)
{
	std::string line;
	// read coordinate system
	in >> line;
	if (line.find("GOCAD_ORIGINAL_COORDINATE_SYSTEM") != std::string::npos) {
		std::string name_of_coordinate_system, axis_name[3], axis_unit[3],
				zpositive;
//		std::cout << "reading coordinate system" << std::endl;
		in >> line;
		if (line.find("NAME") != std::string::npos)
			in >> name_of_coordinate_system;
		in >> line;
		if (line.find("AXIS_NAME") != std::string::npos)
			in >> axis_name[0] >> axis_name[1] >> axis_name[2];
		in >> line;
		if (line.find("AXIS_UNIT") != std::string::npos)
			in >> axis_unit[0] >> axis_unit[1] >> axis_unit[2];
		in >> line;
		if (line.find("ZPOSITIVE") != std::string::npos)
			in >> zpositive;
		in >> line;
		if (line.find("END_ORIGINAL_COORDINATE_SYSTEM") != std::string::npos)
//			std::cout << "ok" << std::endl;
			return;

//		std::cout << "name of coordinate system: " << name_of_coordinate_system
//				<< std::endl;
//		std::cout << "axis names: " << axis_name[0] << " " << axis_name[1]
//				<< " " << axis_name[2] << std::endl;
//		std::cout << "axis units: " << axis_unit[0] << " " << axis_unit[1]
//				<< " " << axis_unit[2] << std::endl;
//		std::cout << "z axis is: " << zpositive << std::endl;
	}
}

bool GocadInterface::isProperty (const std::string &line, std::istream &in)
{
	if (line.find("PROPERTIES") != std::string::npos) {
		std::cout << "found following PROPERTIES: " << std::flush;
		std::list<std::string> str_list(splitString(line, ' '));
		std::list<std::string>::iterator it(str_list.begin());
//		it++;
		while (it != str_list.end()) std::cout << *it++ << " " << std::flush;
		std::cout << std::endl;
		return true;
	}
	if (line.find("PROPERTY_CLASS_HEADER") != std::string::npos) {
		std::cout << "found tag PROPERTY_CLASS_HEADER" << std::endl;
		std::list<std::string> str_list(splitString(line, ' '));
		std::list<std::string>::iterator it(str_list.begin());
//		while (it != str_list.end())
//			std::cout << *it++ << " " << std::flush;
//		std::cout << std::endl;

		char buffer[MAX_COLS_PER_ROW];
		in.getline(buffer, MAX_COLS_PER_ROW);
		std::string str (buffer);
		while (in && str.find("}") == std::string::npos) {
			std::list<std::string> property_str_list(splitString(str, ':'));
			it = property_str_list.begin();
//			while (it != property_str_list.end())
//				std::cout << *it++ << " " << std::flush;
//			std::cout << std::endl;
			in.getline(buffer, MAX_COLS_PER_ROW);
			str = buffer;
		}
		return true;
	}
	if (line.find("PROPERTY_CLASSES") != std::string::npos) {
		std::cout << "found tag PROPERTY_CLASSES" << std::endl;
		return true;
	}
	return false;
}

bool GocadInterface::isGeologicInformation (const std::string &line)
{
	std::string geologic_type, geologic_feature, stratigraphic_age;
	unsigned stratigraphic_time;

	// GEOLOGICAL_TYPE, GEOLOGICAL_FEATURE, STRATIGRAPHIC_POSITION
	if (line.find("GEOLOGICAL_TYPE") != std::string::npos) {
		std::cout << "found following GEOLOGICAL_TYPE: " << std::flush;
		std::list<std::string> str_list(splitString(line, ' '));
		std::list<std::string>::iterator it(str_list.begin()++);
		geologic_type = *it;
		std::cout << geologic_type << std::endl;
		return true;
	}
	if (line.find("GEOLOGICAL_FEATURE") != std::string::npos) {
		std::cout << "found following GEOLOGICAL_FEATURE: " << std::flush;
		std::list<std::string> str_list(splitString(line, ' '));
		std::list<std::string>::iterator it(str_list.begin()++);
		geologic_feature = *it;
		std::cout << geologic_feature << std::endl;
		return true;
	}
	if (line.find("STRATIGRAPHIC_POSITION") != std::string::npos) {
		std::cout << "found following STRATIGRAPHIC_POSITION: " << std::flush;
		std::list<std::string> str_list(splitString(line, ' '));
		std::list<std::string>::iterator it(str_list.begin()++);
		stratigraphic_age = *it++;
		char *res;
		stratigraphic_time = strtol (it++->c_str(), &res, 10);
		std::cout << stratigraphic_age << " " << stratigraphic_time << std::endl;
		return true;
	}
	return false;
}

void GocadInterface::readTSurfData (std::istream &in)
{
	char buffer[MAX_COLS_PER_ROW];
	bool end_not_read (true);

	// for index mapping: first arg is file_idx, second arg is idx in vec
	std::map<size_t, size_t> idx_map;

	// tmp array for triangles
	std::vector<size_t*> trgl;

	while (in && end_not_read) {
		in.getline(buffer, MAX_COLS_PER_ROW);
		std::string line(buffer);

		std::list<std::string> str_list(splitString(line, ' '));
		std::list<std::string>::iterator it(str_list.begin());
		if (it->find("VRTX") != std::string::npos) {
//			std::cout << "reading vertex ... " << std::flush;
			it++;
			char *res;
			size_t id (strtol (it++->c_str(), &res, 10));
			idx_map.insert (std::pair<size_t,size_t>(id, _pnt_vec->size()));
			_pnt_vec->push_back (new Point (strtod(it++->c_str(), &res), strtod(it++->c_str(), &res), strtod(it++->c_str(), &res)));
			id = _pnt_vec->size() - 1;
//			std::cout << id << " " << *((*_pnt_vec)[id]) << " " << std::flush;
//			if (it != str_list.end())
//				std::cout << " - additional info: " << std::flush;
//			while (it != str_list.end())
//				std::cout << *it++ << " " << std::flush;
//			std::cout << std::endl;
		} else if (it->find("TRGL") != std::string::npos) {
//			std::cout << "reading triangle ... " << std::flush;
			it++;
			char *res;
			size_t trgl_pos (trgl.size());
			trgl.push_back (new size_t[3]);
			trgl[trgl_pos][0] = strtol(it++->c_str(), &res, 10);
			trgl[trgl_pos][1] = strtol(it++->c_str(), &res, 10);
			trgl[trgl_pos][2] = strtol(it++->c_str(), &res, 10);

//			std::cout << trgl[trgl_pos][0] << " " << trgl[trgl_pos][1] << " " << trgl[trgl_pos][2] << std::flush;
//			if (it != str_list.end())
//				std::cout << " - additional information: " << std::flush;
//			while (it != str_list.end())
//				std::cout << *it++ << " " << std::flush;
//			std::cout << std::endl;
		} else if (it->find("BSTONE") != std::string::npos) {
			// do something
		} else if (it->find("BORDER") != std::string::npos) {
			// do something with border
		} else if (it->find("END") != std::string::npos) {
			end_not_read = false;
		} else end_not_read = false; // read unknown information tag
	}

	// go trough triangles and change idx
	for (size_t k(0); k<trgl.size(); k++) {
		for (size_t j(0); j<3; j++) {
			// idx read from file
			size_t file_idx(trgl[k][j]);
			// set new point idx to the triangle
			trgl[k][j] = idx_map [file_idx];
		}
	}

	for (size_t k(0); k<trgl.size(); k++) {
		// create Polyline
		Polyline *ply (new Polyline (*_pnt_vec));
		for (size_t j(0); j<3; j++) ply->addPoint (trgl[k][j]);
		// add Polyline to vec
		_ply_vec->push_back (ply);
		// delete data structure
		delete [] trgl[k];
	}
//	for (size_t k(0); k<_ply_vec->size(); k++) {
//		std::cout << k << ": ";
//		((*_ply_vec)[k])->write (std::cout);
//		std::cout << std::endl;
//	}
}

//void GocadInterface::readPLineData (std::istream &in)
//{
//	std::string line;
//	in >> line;
//}

} // end namespace GEOLIB
