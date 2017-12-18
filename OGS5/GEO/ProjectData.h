/**
 * \file ProjectData.h
 * 25/08/2010 KR Initial implementation
 */

#ifndef PROJECTDATA_H_
#define PROJECTDATA_H_

#include "GEOObjects.h"
#include "msh_mesh.h"
#include "FEMCondition.h"


/**
 * The ProjectData Object contains all the data needed for a certain project, i.e. all
 * geometric data (stored in a GEOObjects-object), all the meshes, FEM Conditions (i.e.
 * Boundary Conditions, Source Terms and Initial Conditions), etc.
 * ProjectData does not administrate any of the objects, it is just a "container class" 
 * to store them all in one place.
 * For each class of object stored in this container exists an add-, get- and remove-method.
 *
 * \sa GEOModels, FEMCondition
 */
class ProjectData
{
public:
	ProjectData();
	virtual ~ProjectData();

	// Returns the GEOObjects containing all points, polylines and surfaces
	GEOLIB::GEOObjects& getGEOObjects() const;

	/// Adds a new mesh
	virtual void addMesh(Mesh_Group::CFEMesh* mesh, std::string &name);

	/// Returns the mesh with the given name.
	const Mesh_Group::CFEMesh* getMesh(const std::string &name) const;

	/// Removes the mesh with the given name.
	virtual bool removeMesh(const std::string &name);

	/// Adds a new FEM Condition
	virtual void addCondition(FEMCondition* cond);

	/// Returns the FEM Condition set on a GeoObject with the given name.
	const FEMCondition* getCondition(const std::string &name) const;

	/// Removes the FEM Condition set on a GeoObject with the given name.
	virtual bool removeCondition(const std::string &name);

	/// Checks if the name of the mesh is already exists, if so it generates a unique name.
	bool isUniqueMeshName(std::string &name);

private:
	GEOLIB::GEOObjects _geoObjects;
	std::map<std::string, Mesh_Group::CFEMesh*> _msh_vec;
	std::vector<FEMCondition*> _cond_vec;
};

#endif //PROJECTDATA_H_
