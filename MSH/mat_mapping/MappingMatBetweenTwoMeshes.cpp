/*!
   MappingMatBetweenTwoMeshes.h

   Defenitions of member functions for
   class MappingMatBetweenTwoMeshes, which is for
   mapping material data from one OGS mesh to another OGS mesh

   WW 07.09.2011
 */
#include "MappingMatBetweenTwoMeshes.h"

#include "msh_lib.h"
#include "msh_mesh.h"

#include <fstream>
#include <sstream>
#include <vector>

namespace Mesh_Group
{
using namespace std;
//----------------------------------------------------------------------
/*!
    Constructor
 */
MappingMatBetweenTwoMeshes
::MappingMatBetweenTwoMeshes(const string f_path, const string f_name)
	: mesh_a(NULL), mesh_b(NULL)
{
	n_dx = 0;
	n_dy = 0;
	n_dz = 0;
	x0 = 0.;
	y0 = 0.;
	z0 = 0.;
	cell_size = 0;

	file_name = f_name;
	file_path = f_path;
}

//----------------------------------------------------------------------
/*!
    Destructor
 */
MappingMatBetweenTwoMeshes::~MappingMatBetweenTwoMeshes()
{
	delete mesh_a;
	delete mesh_b;

	mesh_a = mesh_a = NULL;
}

//----------------------------------------------------------------------
/*!
    Read material data
 */
void MappingMatBetweenTwoMeshes:: ReadMaterialData(string f_name, vector<double>& mat_vec)
{
	int id;
	double val;
	string new_file_name = file_path + f_name + "_new.txt";

	f_name = file_path + f_name;
	/// Open a file that contains the file names of the two meshes
	ifstream ins(f_name.c_str());

	if(!ins.good())
	{
		cout << "File " << f_name << " is not found" << endl;
		exit(1);
	}

	ofstream os(new_file_name.c_str(), ios::trunc);

	string aline;
	std::stringstream ss;
	while(!ins.eof())
	{
		getline(ins, aline);
		if(aline.length() == 0)
			continue;

		if(os.good())
			os << aline << endl;

		if(aline.find("$DATA") != string::npos)
		{
			os.close();

			for(;; )
			{
				getline(ins, aline);
				if(aline.find("#STOP") != string::npos)
					break;
				ss.str(aline);
				ss >> id >> val;
				ss.clear();
				mat_vec.push_back(val);
			}
		}
	}
}

//----------------------------------------------------------------------
/*!
    Read meshes and material data
 */
void MappingMatBetweenTwoMeshes::ReadMeshes_and_MaterialData()
{
	int count;
	string in_file;
	string file_a;
	string file_b;

	in_file = file_name + ".2mesh";
	/// Open a file that contains the file names of the two meshes
	ifstream ins(in_file.c_str());

	if(!ins.good())
	{
		cout << "File " << in_file << " is not found" << endl;
		exit(1);
	}

	string aline;
	std::stringstream ss;

	count = 0;
	while(!ins.eof())
	{
		getline(ins, aline);
		if(aline.length() == 0)
			continue;

		ss.str(aline);
		ss >> aline;
		names.push_back(aline);
		ss.clear();
		count++;

		if(count == 4)
			break;
	}
	ins.close();

	file_a = names[0];
	file_a = file_path + file_a;
	//NULL, NULL: non geometry object, no additional file name
	mesh_a = FEMRead(file_a, NULL, NULL);
	mesh_a->ConstructGrid();

	file_b = names[1];
	file_b = file_path + file_b;
	mesh_b = FEMRead(file_b, NULL, NULL);
	mesh_b->ConstructGrid();

	ReadMaterialData(names[2], perm);
	ReadMaterialData(names[3], poro);
}

//----------------------------------------------------------------------
/*!
    Read meshes
 */
void MappingMatBetweenTwoMeshes::MappingMaterialData()
{
	size_t ii;
	long i, j, k, l, m, n;
	long ix, iy, iz;
	CElem* elem;
	CElem* elem_nb;
	CNode* node;
	double* grav_c;
	double min_xyz[3], max_xyz[3];
	double cell_size_z = 0.0;

	mesh_b->getMinMax_XYZ(min_xyz, max_xyz);

	double vol;
	double cell_size = DBL_MAX;
	int dim = mesh_a->GetMaxElementDim();
	long size_a, size_b;
	size_a = static_cast<long>(mesh_a->ele_vector.size());
	size_b = static_cast<long>(mesh_b->ele_vector.size());

	for(i = 0; i < size_a; i++)
	{
		elem = mesh_a->ele_vector[i];
		vol = elem->GetVolume();
		if(cell_size > vol)
			cell_size = vol;
	}
	cell_size = 10. * pow(cell_size, 1.0 / static_cast<float>(dim)); //0.5*
	cell_size_z = 0.01 * cell_size;

	x0 = min_xyz[0];
	y0 = min_xyz[1];
	z0 = min_xyz[2];

	n_dx = static_cast<long>((max_xyz[0] - x0) / cell_size);
	n_dy = static_cast<long>((max_xyz[1] - y0) / cell_size);
	n_dz = static_cast<long>((max_xyz[2] - z0) / cell_size_z);
	if(n_dx < 1) // if mesh b is  2D
		n_dx = 1;
	if(n_dy < 1) // if mesh b is  2D
		n_dy = 1;
	if(n_dz < 1) // if mesh b is  2D
		n_dz = 1;

	double def_v = -9999.0;
	double* perm_value, * poro_value;
	long size;
	long lx_min, lx_max;
	long ly_min, ly_max;
	long lz_min, lz_max;

	size = n_dx * n_dy * n_dz;

	while(size > LONG_MAX || size < 0) // Integer overflow
	{
		cell_size *= 8;
		cell_size_z *= 2;
		n_dx /= 8;
		n_dy /= 8;
		n_dz /= 2;
		size = n_dx * n_dy * n_dz;
	}

	perm_value = new double[size];
	poro_value = new double[size];
	for(i = 0; i < size; i++)
	{
		perm_value[i] = def_v;
		poro_value[i] = def_v;
	}

	for(i = 0; i < size_a; i++)
	{
		elem = mesh_a->ele_vector[i];
		grav_c = elem->GetGravityCenter();

		ix = static_cast<long>((grav_c[0] - x0) / cell_size);
		iy = static_cast<long>((grav_c[1] - y0) / cell_size);
		iz = static_cast<long>((grav_c[2] - z0) / cell_size_z);

		if(ix < 0 || ix > n_dx - 1)
			continue;
		if(iy < 0 || iy > n_dy - 1)
			continue;
		if(iz < 0 || iz > n_dz - 1)
			continue;

		k = n_dx * n_dy * iz + n_dx * iy + ix;
		perm_value[k] = perm[i];
		poro_value[k] = poro[i];

		for(ii = 0; ii < 3; ii++)
		{
			max_xyz[ii] = DBL_MIN;
			min_xyz[ii] = DBL_MAX;
		}
		for(ii = 0; ii < elem->GetNodesNumber(false); ii++)
		{
			node = elem->GetNode(ii);
			if(max_xyz[0] < node->X())
				max_xyz[0] = node->X();
			if(max_xyz[1] < node->Y())
				max_xyz[1] = node->Y();
			if(max_xyz[2] < node->Z())
				max_xyz[2] = node->Z();

			if(min_xyz[0] > node->X())
				min_xyz[0] = node->X();
			if(min_xyz[1] > node->Y())
				min_xyz[1] = node->Y();
			if(min_xyz[2] > node->Z())
				min_xyz[2] = node->Z();
		}

		// Spur surround occupied by this element
		//
		lz_min = static_cast<long>((min_xyz[2] - z0) / cell_size_z);
		lz_max = static_cast<long>((max_xyz[2] - z0) / cell_size_z);
		ly_min = static_cast<long>((min_xyz[1] - y0) / cell_size);
		ly_max = static_cast<long>((max_xyz[1] - y0) / cell_size);
		lx_min = static_cast<long>((min_xyz[0] - x0) / cell_size);
		lx_max = static_cast<long>((max_xyz[0] - x0) / cell_size);

		if(lx_min < 0)
			lx_min = 0;
		if(lx_min > n_dx - 1)
			lx_min = n_dx - 1;
		if(lx_max < 0)
			lx_max = 0;
		if(lx_max > n_dx - 1)
			lx_max = n_dx - 1;

		if(ly_min < 0)
			ly_min = 0;
		if(ly_min > n_dy - 1)
			ly_min = n_dy - 1;
		if(ly_max < 0)
			ly_max = 0;
		if(ly_max > n_dy - 1)
			ly_max = n_dy - 1;

		if(lz_min < 0)
			lz_min = 0;
		if(lz_min > n_dz - 1)
			lz_min = n_dz - 1;
		if(lz_max < 0)
			lz_max = 0;
		if(lz_max > n_dz - 1)
			lz_max = n_dz - 1;

		for(j = lz_min; j < lz_max; j++)
		{
			for(n = ly_min; n < ly_max; n++)
				for(l = lx_min; l < lx_max; l++)
				{
					m = n_dx * n_dy * j + n_dx * n + l;
					if(fabs(perm_value[m] - def_v) < DBL_MIN)
					{
						perm_value[m] =  perm_value[k];
						poro_value[m] =  poro_value[k];
					}
				}
		}
	}

	vector<double> perm_b(size_b);
	vector<double> poro_b(size_b);

	for(i = 0; i < size_b; i++)
	{
		elem = mesh_b->ele_vector[i];
		grav_c = elem->GetGravityCenter();

		ix = static_cast<long>((grav_c[0] - x0) / cell_size);
		iy = static_cast<long>((grav_c[1] - y0) / cell_size);
		iz = static_cast<long>((grav_c[2] - z0) / cell_size_z);

		if(ix > n_dx - 1)
			continue;
		if(iy > n_dy - 1)
			continue;
		if(iz > n_dz - 1)
			continue;

		if(ix < 0 || ix > n_dx - 1)
		{
			perm_b[i] = def_v;
			poro_b[i] = def_v;
		}
		if(iy < 0 || iy > n_dy - 1)
		{
			perm_b[i] = def_v;
			poro_b[i] = def_v;
		}
		if(iz < 0 || iz > n_dz - 1)
		{
			perm_b[i] = def_v;
			poro_b[i] = def_v;
		}

		k = n_dx * n_dy * iz + n_dx * iy + ix;
		perm_b[i] = perm_value[k];
		poro_b[i] = poro_value[k];
	}

	// Check again
	long counter;
	for(;; )
	{
		counter = 0;
		for(i = 0; i < size_b; i++)
		{
			if(fabs(perm_b[i] - def_v) > DBL_MIN)
			{
				for(j = 0; j < elem->GetFacesNumber(); j++)
				{
					elem_nb = elem->GetNeighbor(j);
					if(!elem_nb)
						continue;
					if(elem_nb->GetDimension() != elem->GetDimension())
						continue;
					l = elem_nb->GetIndex();
					if(fabs(perm_b[l] - def_v) < DBL_MIN)
					{
						perm_b[l] = perm_b[i];
						poro_b[l] = poro_b[i];
					}
				}
				counter++;
				continue;
			}
			elem = mesh_b->ele_vector[i];
			for(j = 0; j < elem->GetFacesNumber(); j++)
			{
				elem_nb = elem->GetNeighbor(j);
				if(!elem_nb)
					continue;
				if(elem_nb->GetDimension() != elem->GetDimension())
					continue;
				l = elem_nb->GetIndex();
				if(fabs(perm_b[l] - def_v) > DBL_MIN)
				{
					perm_b[i] = perm_b[l];
					poro_b[i] = poro_b[l];
					break;
				}
			}
		}
		if(counter == size_b)
			break;
	}

	names[2] = file_path + names[2] + "_new.txt";
	names[3] = file_path + names[3] + "_new.txt";
	ofstream os_k(names[2].c_str(), ios::out | ios::app);
	ofstream os_n(names[3].c_str(), ios::out | ios::app);
	for(i = 0; i < size_b; i++)
	{
		os_k << i << "  " << perm_b[i] * 1.e-4 << endl;
		os_n << i << "  " << poro_b[i] << endl;
	}
	os_k << "#STOP" << endl;
	os_n << "#STOP" << endl;

	delete [] perm_value;
	delete [] poro_value;
}
} // End of Mesh_Group
