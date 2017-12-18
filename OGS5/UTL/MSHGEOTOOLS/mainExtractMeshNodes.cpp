/*
 * extractMeshNodes.cpp
 *
 *  Created on: Dec 6, 2010
 *      Author: TF
 *
 */

#include "ExtractMeshNodes.h"

// Base
#include "StringTools.h"

// FEM
#include "problem.h"

// MSH
#include "msh_mesh.h"
#include "msh_lib.h" // for FEMRead

// GEO
#include "GEOObjects.h"
#include "PolylineVec.h"
#include "Polygon.h"

// FileIO
#include "OGSIOVer4.h"

Problem *aproblem = NULL;

// MATHLIB
#include "Matrix.h"
#include "Vector3.h"

int main (int argc, char *argv[])
{
	if (argc < 5) {
		std::cout << "Usage: " << argv[0] << " --mesh ogs_meshfile --geometry ogs_geometry" << std::endl;
		return -1;
	}

	// *** read mesh
	std::string tmp (argv[1]);
	if (tmp.find ("--mesh") == std::string::npos) {
		std::cout << "could not extract mesh file name" << std::endl;
		return -1;
	}

	tmp = argv[2];
	std::string file_base_name (tmp);
	if (tmp.find (".msh") != std::string::npos)
		file_base_name = tmp.substr (0, tmp.size()-4);

	Mesh_Group::CFEMesh* mesh (FEMRead(file_base_name));
	if (!mesh) {
		std::cerr << "could not read mesh from file " << std::endl;
		return -1;
	}

	// *** read geometry
	tmp = argv[3];
	if (tmp.find ("--geometry") == std::string::npos) {
		std::cout << "could not extract geometry file name" << std::endl;
		return -1;
	}

	GEOLIB::GEOObjects *geo (new GEOLIB::GEOObjects);
	tmp = argv[4];
	FileIO::readGLIFileV4(tmp, geo);

	// *** get Polygon
	const std::vector<GEOLIB::Polyline*>* plys (geo->getPolylineVec (tmp));
	if (!plys) {
		std::cout << "could not get vector of polylines" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	Mesh_Group::ExtractMeshNodes extract_mesh_nodes (mesh);

	std::string fname ("MeshIds.txt");
	std::ofstream out (fname.c_str());

	std::string fname_gli ("MeshNodesAsPnts.gli");
	std::ofstream pnt_out (fname_gli.c_str());
	pnt_out << "#POINTS" << std::endl;

	const size_t n_plys (plys->size());
	for (size_t k(0); k<n_plys; k++) {
		bool closed ((*plys)[k]->isClosed());
		if (!closed) {
			std::cout << "polyline " << k << " is not closed" << std::endl;
		} else {
			GEOLIB::Polygon polygon (*((*plys)[k]));
			extract_mesh_nodes.writeTopSurfaceMeshNodeIDs (out, pnt_out, polygon);
			// write all nodes - not only the surface nodes
//			extract_mesh_nodes.writeMeshNodeIDs (out, pnt_out, polygon);
		}
	}
	pnt_out << "#STOP" << std::endl;
	pnt_out.close ();

	delete mesh;
	delete geo;

	return 0;
}
