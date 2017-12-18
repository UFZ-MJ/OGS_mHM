/*!
   MappingMatBetweenTwoMeshes.h

   A class for
   Mapping material data from one OGS mesh to another OGS mesh

   WW 07.09.2011
 */

#ifndef MAPPING_MAT_BETWEEN_TWO_MESHES_H_
#define MAPPING_MAT_BETWEEN_TWO_MESHES_H_

#include <string>
#include <vector>

namespace Mesh_Group
{
class CFEMesh;

using std::string;
using std::vector;

class MappingMatBetweenTwoMeshes
{
public:
	MappingMatBetweenTwoMeshes(const string f_path, const string f_name);
	~MappingMatBetweenTwoMeshes();

	void ReadMeshes_and_MaterialData();   // Init
	void ReadMaterialData(string f_name,  vector<double>& mat_vec);

	void MappingMaterialData();

private:
	/// Data for geometry and grid
	long n_dx;
	long n_dy;
	long n_dz;
	/// Coodinates of the low left corner of the grid
	double x0;
	double y0;
	double z0;
	/// Cell size
	double cell_size;

	/// Material value of elements
	vector<double> perm;
	vector<double> poro;

	CFEMesh* mesh_a;
	CFEMesh* mesh_b;

	string file_name;
	string file_path;

	vector<string> names;
};
} //End of Mesh_Group
#endif
