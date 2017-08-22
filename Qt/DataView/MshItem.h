/**
 * \file MshItem.h
 * 17/05/2010 KR Initial implementation
 */

#ifndef MSHITEM_H
#define MSHITEM_H

#include "TreeItem.h"
#include "VtkMeshSource.h"

class GridAdapter;
class VtkMeshSource;

/**
 * \brief A TreeItem containing a mesh and the associated vtk object used in the Mesh Model.
 * \sa TreeItem
 */
class MshItem : public TreeItem
{

public:
	/// Constructor, automatically generates VTK object of the given mesh.
	MshItem(const QList<QVariant> &data, TreeItem *parent, GridAdapter* grid);
	~MshItem();

	/// Returns the mesh as a GridAdapter.
	const GridAdapter* getGrid() const { return this->_meshSource->GetGrid(); };
	/// Returns the VTK object.
	VtkMeshSource* vtkSource() const { return _meshSource; };	


private:
	VtkMeshSource* _meshSource;
};

#endif //MSHITEM_H
