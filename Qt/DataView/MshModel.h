/**
 * \file MshModel.h
 * 19/10/2009 LB Initial implementation
 * 12/05/2010 KR re-implementation
 *
 */


#ifndef MSHMODEL_H
#define MSHMODEL_H

// ** INCLUDES **
#include "TreeModel.h"
#include "ProjectData.h"
#include "GridAdapter.h"

class VtkMeshSource;

/**
 * The MshModel is a Qt model which represents Mesh objects.
 */
class MshModel : public TreeModel
{
	Q_OBJECT

public:
	MshModel(ProjectData &project, QObject* parent = 0);

	/// Returns the number of columns used for the data list
	int columnCount(const QModelIndex& parent = QModelIndex()) const;

	static Mesh_Group::CFEMesh* loadMeshFromFile(std::string fileName);

public slots:
	/// Adds a new mesh
	void addMesh(Mesh_Group::CFEMesh* mesh, std::string &name);
	/// Returns the mesh with the given index.
	const GridAdapter* getMesh(const QModelIndex &idx) const;
	/// Returns the mesh with the given name.
	const GridAdapter* getMesh(const std::string &name) const;
	/// Removes the mesh with the given index.
	bool removeMesh(const QModelIndex &idx);
	/// Removes the mesh with the given name.
	bool removeMesh(const std::string &name);
	/// Returns the VTK source item for the mesh with the given index.
	VtkMeshSource* vtkSource(const QModelIndex &idx) const;
	/// Returns the VTK source item for the mesh with the given name.
	VtkMeshSource* vtkSource(const std::string &name) const;

private:
	/// Adds the mesh to the GUI-Mesh-Model und -View
	void addMesh(GridAdapter* mesh, std::string &name);

	/// Checks if the name of the mesh is already exists, if so it generates a unique name.
	//bool isUniqueMeshName(std::string &name);
	ProjectData& _project;

signals:
	void meshAdded(MshModel*, const QModelIndex&);
	void meshRemoved(MshModel*, const QModelIndex&);

};

#endif // MSHMODEL_H
