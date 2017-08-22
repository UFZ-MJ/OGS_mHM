/**
 * \file GMSInterface.cpp
 * 08/06/2010 KR Initial implementation
 *
 */

#include <fstream>
#include "GMSInterface.h"
#include "MSHEnums.h"
#include "StringTools.h"

int GMSInterface::readBoreholesFromGMS(std::vector<GEOLIB::Point*> *boreholes, const std::string &filename)
{
	double depth(0.0);
	std::string line(""), cName(""), sName("");
	std::list<std::string>::const_iterator it;
	GEOLIB::Point* pnt = new GEOLIB::Point();
	GEOLIB::StationBorehole* newBorehole = NULL;
	std::ifstream in( filename.c_str() );

	if (!in.is_open())
    {
		std::cout << "GMSInterface::readBoreholeFromGMS() - Could not open file...\n";
		return 0;
	}

	/* skipping first line because it contains field names */
	getline(in, line);

	/* read all stations */
	while ( getline(in, line) )
	{
		std::list<std::string> fields = splitString(line, '\t');

		if (fields.size() >= 5)
		{
			if (fields.begin()->compare(cName) == 0) // add new layer
			{
				it = fields.begin();
				(*pnt)[0] = strtod((++it)->c_str(), 0);
				(*pnt)[1] = strtod((++it)->c_str(), 0);
				(*pnt)[2] = strtod((++it)->c_str(), 0);
				newBorehole->addSoilLayer((*pnt)[0], (*pnt)[1], (*pnt)[2], sName);
				sName = (*(++it));
				depth=(*pnt)[2];
			}
			else // add new borehole
			{
				if (newBorehole != NULL)
				{
					newBorehole->setDepth((*newBorehole)[2]-depth);
					boreholes->push_back(newBorehole);
				}
				cName = *fields.begin();
				it = fields.begin();
				(*pnt)[0] = strtod((++it)->c_str(), 0);
				(*pnt)[1] = strtod((++it)->c_str(), 0);
				(*pnt)[2] = strtod((++it)->c_str(), 0);
				sName = (*(++it));
				newBorehole = GEOLIB::StationBorehole::createStation(cName, (*pnt)[0], (*pnt)[1], (*pnt)[2], 0);
			}
		}
		else
			std::cout << "GMSInterface::readBoreholeFromGMS() - Error reading format..." << std::endl;
	}
	// write the last borehole from the file
	if (newBorehole != NULL)
	{
		newBorehole->setDepth((*newBorehole)[2]-depth);
		boreholes->push_back(newBorehole);
	}

	in.close();
	return 1;
}

/*
// all boreholes to GMS which each borehole in a single file
void StationIO::writeBoreholesToGMS(const std::vector<GEOLIB::Point*> *stations)
{
	//std::vector<std::string> soilID(1);
	std::vector<std::string> soilID = readSoilIDfromFile("d:/BodeTimeline.txt");
	for (size_t i=0; i<stations->size(); i++)
		StationIO::writeBoreholeToGMS(static_cast<GEOLIB::StationBorehole*>((*stations)[i]), std::string("Borehole-" + static_cast<GEOLIB::StationBorehole*>((*stations)[i])->getName() + ".txt"), soilID);
	StationIO::writeSoilIDTable(soilID, "SoilIDReference.txt");
}
*/
void GMSInterface::writeBoreholesToGMS(const std::vector<GEOLIB::Point*> *stations, const std::string &filename)
{
	std::ofstream out( filename.c_str(), std::ios::out );
	size_t idx = 0;
	std::vector<std::string> soilID = readSoilIDfromFile("d:/BodeTimeline.txt");

	// write header
	out	<< "name" << "\t" << std::fixed << "X" << "\t" << "Y"  << "\t" << "Z" <<  "\t" << "soilID" << std::endl;

	for (size_t j=0; j<stations->size(); j++)
	{

		GEOLIB::StationBorehole* station = static_cast<GEOLIB::StationBorehole*>((*stations)[j]);
		std::vector<GEOLIB::Point*> profile = station->getProfile();
		std::vector<std::string> soilNames  = station->getSoilNames();

		size_t nLayers = profile.size();
		for (size_t i=1; i<nLayers; i++) {

			if ( (i>1) && (soilNames[i].compare(soilNames[i-1]) == 0) ) continue;
			idx = getSoilID(soilID, soilNames[i]);

			out	<< station->getName() << "\t" << std::fixed << (*(profile[i-1]))[0] << "\t"
				<< (*(profile[i-1]))[1]  << "\t" << (*(profile[i-1]))[2] <<  "\t"
				<< idx << std::endl;
		}
		out	<< station->getName() << "\t" << std::fixed << (*(profile[nLayers-1]))[0] << "\t"
				<< (*(profile[nLayers-1]))[1]  << "\t" << (*(profile[nLayers-1]))[2] <<  "\t"
				<< idx << std::endl;	// this line marks the end of the borehole
	}

	out.close();
	GMSInterface::writeSoilIDTable(soilID, "d:/SoilIDReference.txt");
}


int GMSInterface::writeBoreholeToGMS(const GEOLIB::StationBorehole* station, const std::string &filename, std::vector<std::string> &soilID)
{
	std::ofstream out( filename.c_str(), std::ios::out );
	size_t idx = 0;

	// write header
	out	<< "name" << "\t" << std::fixed << "X" << "\t" << "Y"  << "\t" << "Z" <<  "\t" << "soilID" << std::endl;

	std::vector<GEOLIB::Point*> profile = station->getProfile();
	std::vector<std::string> soilNames  = station->getSoilNames();

	// write table
	size_t nLayers = profile.size();
	for (size_t i=1; i<nLayers; i++) {

		if ( (i>1) && (soilNames[i].compare(soilNames[i-1]) == 0) ) continue;
		idx = getSoilID(soilID, soilNames[i]);

		out	<< station->getName() << "\t" << std::fixed << (*(profile[i-1]))[0] << "\t"
			<< (*(profile[i-1]))[1]  << "\t" << (*(profile[i-1]))[2] <<  "\t"
			<< idx << std::endl;
	}
	out	<< station->getName() << "\t" << std::fixed << (*(profile[nLayers-1]))[0] << "\t"
				<< (*(profile[nLayers-1]))[1]  << "\t" << (*(profile[nLayers-1]))[2] <<  "\t"
				<< idx << std::endl;	// this line marks the end of the borehole
	out.close();

    return 1;
}


size_t GMSInterface::getSoilID(std::vector<std::string> &soilID, std::string &soilName)
{
	for (size_t j=0; j<soilID.size(); j++)
	{
		if (soilID[j].compare(soilName) == 0) return j;
	}
	soilID.push_back(soilName);
	return (soilID.size() - 1);
}


int GMSInterface::writeSoilIDTable(const std::vector<std::string> &soilID, const std::string &filename)
{
	std::ofstream out( filename.c_str(), std::ios::out );

	// write header
	out	<< "ID" << "\t" << std::fixed << "Soil name"<< std::endl;

	// write table
	size_t nIDs = soilID.size();
	for (size_t i=0; i<nIDs; i++)
		out	<< i << "\t" << std::fixed << soilID[i] << "\t" << std::endl;
	out.close();

    return 1;
}

std::vector<std::string> GMSInterface::readSoilIDfromFile(const std::string &filename)
{
	std::vector<std::string> soilID;
	std::string line;

	std::ifstream in( filename.c_str() );

	if (in.is_open())
	{
		while ( getline(in, line) )
		{
			trim(line);
			soilID.push_back(line);
		}
	}
	in.close();

	return soilID;
}


Mesh_Group::CFEMesh* GMSInterface::readGMS3DMMesh(std::string filename)
{
	std::string buffer("");

	std::ifstream in(filename.c_str());
	if (!in.is_open())
    {
		std::cout << "GMSInterface::readGMS3DMMesh() - Could not open file..." << std::endl;
		return NULL;
	}

	// Read data from file
	getline(in, buffer); // "MESH3D"
	if (buffer.compare("MESH3D") != 0)
	{
		std::cout << "GMSInterface::readGMS3DMMesh() - Could not read expected file header..." << std::endl;
		return NULL;
	}

	std::cout << "Read GMS-3DM data...";
	Mesh_Group::CFEMesh* mesh (new Mesh_Group::CFEMesh());

	while (!in.eof())
	{
		Mesh_Group::CElem* elem = new Mesh_Group::CElem();
		std::string element_id("");
		in >> element_id;

		if (element_id.compare("E6W") == 0)
		{
			elem->SetElementType(MshElemType::PRISM);
			elem->Read(in, 8);
			mesh->ele_vector.push_back(elem);
		}
		else if (element_id.compare("E4T") == 0)
		{
			elem->SetElementType(MshElemType::TETRAHEDRON);
			elem->Read(in, 8);
			mesh->ele_vector.push_back(elem);

		}
		else if (element_id.compare("E4P") == 0)
		{
			int i(0);
			long node_index[5];
			elem->SetElementType(MshElemType::TETRAHEDRON);
			elem->Read(in, 8);
			mesh->ele_vector.push_back(elem);
			for (size_t j=0; j<4; j++)
				node_index[j] = elem->GetNodeIndex(j);
			in >> node_index[4];
			i = elem->GetPatchIndex();
			elem->SetPatchIndex(node_index[4]-1);
			node_index[4] = i;

			elem->SetNodeIndex(0, node_index[0]);
			elem->SetNodeIndex(1, node_index[1]);
			elem->SetNodeIndex(2, node_index[3]);
			elem->SetNodeIndex(3, node_index[4]);

			Mesh_Group::CElem* elem2 = new Mesh_Group::CElem(mesh->ele_vector.size()-1, elem);
			elem2->SetNodeIndex(0, node_index[1]);
			elem2->SetNodeIndex(1, node_index[2]);
			elem2->SetNodeIndex(2, node_index[3]);
			elem2->SetNodeIndex(3, node_index[4]);
			mesh->ele_vector.push_back(elem2);
		}
		else if (element_id.compare("ND") == 0)
		{
			int i;
			double xyz[3];
			in >> i;
			Mesh_Group::CNode* node = new Mesh_Group::CNode(i-1);
			in >> xyz[0] >> xyz[1] >> xyz[2] >> std::ws;
			node->SetCoordinates(xyz);
			mesh->nod_vector.push_back(node);
		}
		else //default
		{
			std::cout << std::endl << "GMSInterface::readGMS3DMMesh() - Unknown identifier ..." << std::endl;
			return NULL;
		}
	}

	in.close();

	mesh->ConstructGrid();
	mesh->FillTransformMatrix();

	std::cout << "finished" << std::endl;

	return mesh;
}
