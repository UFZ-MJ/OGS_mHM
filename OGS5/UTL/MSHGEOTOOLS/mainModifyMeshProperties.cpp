/*
 * modifyMeshProperties.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 *
 */

#include "ModifyMeshProperties.h"

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

int main (int argc, char *argv[])
{
	if (argc < 7) {
		std::cout << "Usage: " << argv[0] << " --mesh ogs_meshfile --geometry ogs_geometry --material_id id" << std::endl;
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

	bool closed ((*plys)[0]->isClosed());
	if (!closed) {
		std::cout << "polyline is not closed" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	// *** get material id
	tmp = argv[5];
	if (tmp.find ("--material") == std::string::npos) {
		std::cout << "could not read material id" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}
	size_t material_id (atoi (argv[6]));

	GEOLIB::Polygon polygon (*((*plys)[0]));
	Mesh_Group::ModifyMeshProperties modify_mesh_nodes (mesh);

	modify_mesh_nodes.setMaterial (polygon, material_id);

	std::fstream mesh_out;
	mesh_out.open ("MeshWithMaterial.msh", std::fstream::out);

	if (mesh_out.is_open()) {
		std::cout << "write TestMesh.msh" << std::endl;
		mesh->Write (&mesh_out);
	}
	mesh_out.close ();

	return 0;
}
