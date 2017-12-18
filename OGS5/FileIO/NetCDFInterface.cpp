/**
 * \file NetCDFInterface.cpp
 * 29/07/2010 YW Initial implementations
 */
#include <stdio.h>
#include <vtknetcdf/netcdf.h>
#include "NetCDFInterface.h"

using namespace GEOLIB;
using namespace Mesh_Group;

/* Names of variables. */
#define DIM_RLAT "rlat"
#define DIM_RLON "rlon"
#define LAT_NAME "lat"
#define LON_NAME "lon"
/* Length of the original output NetCDF file. */
#define LEN_ORIGINAL 18

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define HANDLE_ERROR(e) {printf("Error: %s\n", nc_strerror(e)); exit(-1);}

namespace FileIO {

void NetCDFInterface::readNetCDFData(std::string &ifname, std::vector<GEOLIB::Point*>* points_vec, GEOLIB::GEOObjects* obj, size_t &NRLAT, size_t &NRLON)
{
   int ncid;
   int rlat_dimid, rlon_dimid;
   int lat_varid, lon_varid;
   int choose_varid;

   /* The start and count arrays will tell the netCDF library where to read our data. */
   static size_t start2[2], count2[2];
   static size_t start3[3], count3[3];
   static size_t start4[4], count4[4];
   
   /* Error handling. */
   int status;

  /* Open the file. */
   std::cout << "Open NetCDF file " << ifname << "\n" << std::flush;
   status = nc_open(ifname.c_str(), NC_NOWRITE, &ncid);
   if (status != NC_NOERR)   HANDLE_ERROR(status);

   /* Get the lengths of dimensions, such as "rlat" and "rlat". */
   status = nc_inq_dimid(ncid, DIM_RLAT, &rlat_dimid); 
   if (status != NC_NOERR)   HANDLE_ERROR(status);
   status = nc_inq_dimlen(ncid, rlat_dimid, &NRLAT); 
   if (status != NC_NOERR)   HANDLE_ERROR(status);

   status = nc_inq_dimid(ncid, DIM_RLON, &rlon_dimid); 
   if (status != NC_NOERR)   HANDLE_ERROR(status);
   status = nc_inq_dimlen(ncid, rlon_dimid, &NRLON); 
   if (status != NC_NOERR)   HANDLE_ERROR(status);

   /* Get the choosen variable name from the opened NetCDF file. */
   std::string nc_fname; 
   size_t indexChar; 
   indexChar = ifname.find_last_of('/');
   if(indexChar==std::string::npos)
     nc_fname = ifname;
   else
     nc_fname = ifname.substr(indexChar+1);
   size_t len_varname = nc_fname.size() - LEN_ORIGINAL;   
   std::string var_name;
   var_name = nc_fname.substr(0,len_varname);

   /* Get the varids of the netCDF variables, such as longitude, latitude,
    * temperature or precipitation, and so on. */
   status = nc_inq_varid(ncid, LON_NAME, &lon_varid);
   if (status != NC_NOERR)   HANDLE_ERROR(status);
   status = nc_inq_varid(ncid, LAT_NAME, &lat_varid);
   if (status != NC_NOERR)   HANDLE_ERROR(status);
   status = nc_inq_varid(ncid, var_name.c_str(), &choose_varid);
   if (status != NC_NOERR)   HANDLE_ERROR(status);

   /* Program variables to hold the data we will read. We will only
      need enough space to hold one timestep of data; one record.*/
   size_t len_vals = NRLAT * NRLON;
   float *lat_in = new float[len_vals];
   float *lon_in = new float[len_vals];
   float *var_in = new float[len_vals];

   /* Read the data. Since we know the contents of the file we know
    * that the data arrays in this program are the correct size to
    * hold one timestep. */
   count2[0] = NRLAT;
   count2[1] = NRLON;
   count3[0] = 1;
   count3[1] = NRLAT;
   count3[2] = NRLON;
   count4[0] = 1;
   count4[1] = 1;
   count4[2] = NRLAT;
   count4[3] = NRLON;

   for (size_t i = 0; i < 2; i++)
     start2[i] = 0;
   for (size_t i = 0; i < 3; i++)
     start3[i] = 0;
   for (size_t i = 0; i < 4; i++)
     start4[i] = 0;

   /* Read the above netCDF variables with their corresponding varids. */
   status = nc_get_vara_float(ncid, lon_varid, start2, count2, lon_in);
   if (status != NC_NOERR)   HANDLE_ERROR(status);

   status = nc_get_vara_float(ncid, lat_varid, start2, count2, lat_in);
   if (status != NC_NOERR)   HANDLE_ERROR(status);

   if (var_name.compare("T_2M") == 0)
   {
       status = nc_get_vara_float(ncid, choose_varid, start4, count4, var_in);
       if (status != NC_NOERR)   HANDLE_ERROR(status);
   }
   else if (var_name.compare("TOT_PREC") == 0)
   {
       status = nc_get_vara_float(ncid, choose_varid, start3, count3, var_in);
       if (status != NC_NOERR)   HANDLE_ERROR(status);
   }
    else if (var_name.compare("HSURF") == 0) 
   { 
       status = nc_get_vara_float(ncid, choose_varid, start3, count3, var_in); 
       if (status != NC_NOERR)   HANDLE_ERROR(status); 
   }

   /* Close the file. */
   status = nc_close(ncid);
   if (status != NC_NOERR)   HANDLE_ERROR(status);

   printf("Reading netCDF file successfully!\n");

   for (size_t i = 0; i < len_vals; i++)  
	   points_vec->push_back(new Point( lon_in[i], lat_in[i], var_in[i] ));

   delete [] lat_in;
   delete [] lon_in;
   delete [] var_in;

   obj->addPointVec(points_vec, ifname);
}

CFEMesh* NetCDFInterface::createMeshFromPoints(std::vector<GEOLIB::Point*>* points_vec, size_t &NRLAT, size_t &NRLON)
{	
	std::cout << "Converting points data to mesh with quad elements\n";
	CFEMesh* mesh = new CFEMesh();

	// setmesh nodes
	size_t nNodes = points_vec->size();
	for (size_t i=0; i<nNodes; i++)
	{
		Mesh_Group::CNode* node = new Mesh_Group::CNode(i);
	    double coords[3] = {(*(*points_vec)[i])[0], (*(*points_vec)[i])[1], (*(*points_vec)[i])[2]};;
	    node->SetCoordinates(coords);     
	    mesh->nod_vector.push_back(node);
	}

    // set mesh elements
    for (size_t i = 0; i < NRLAT-1; i++)  
    {
      for (size_t j = 0; j < NRLON-1; j++)
      {
        size_t n_Elems = i * NRLON + j;
		Mesh_Group::CElem* elem = new Mesh_Group::CElem();
		elem->setElementProperties(MshElemType::QUAD);
  	    // Assignment for the element nodes
		elem->nodes_index[0] = (long) (n_Elems);
        elem->nodes_index[1] = (long) (n_Elems + 1);
		elem->nodes_index[2] = (long) (n_Elems + NRLON + 1);
		elem->nodes_index[3] = (long) (n_Elems + NRLON);
        // Initialize topological properties
		elem->InitializeMembers();
		mesh->ele_vector.push_back(elem);
	  }
	}
	// Establish topology of a grid
	mesh->ConstructGrid();
	std::cout << "Mesh is generated successfully!\n" << std::endl;
	return mesh;
}

} // end namespace FileIO
