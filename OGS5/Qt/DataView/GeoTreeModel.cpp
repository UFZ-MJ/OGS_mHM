/**
 * \file GeoTreeModel.cpp
 * 2011/02/07 KR Initial implementation
 */

#include "GeoTreeModel.h"
#include "OGSError.h"
#include "GeoTreeItem.h"
#include "GeoObjectListItem.h"

/**
 * Constructor.
 */
GeoTreeModel::GeoTreeModel( QObject *parent )
: TreeModel(parent)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Id" << "x" << "y" << "z";
	_rootItem = new GeoTreeItem(rootData, NULL, NULL);
}

GeoTreeModel::~GeoTreeModel()
{
}

const GEOLIB::GeoObject* GeoTreeModel::objectFromIndex( const QModelIndex& index, QString &geoName ) const
{
	if (index.isValid())
	{
		GeoTreeItem* treeItem = static_cast<GeoTreeItem*>(index.internalPointer());
		//TreeItem* parentItem = treeItem->parentItem();
		//geoName = parentItem->data(0).toString();
		if (treeItem) return treeItem->getGeoObject();
	}
	return NULL;
}

void GeoTreeModel::addPointList(QString geoName, const GEOLIB::PointVec* pointVec)
{
	const std::vector<GEOLIB::Point*> *points = pointVec->getVector();

	QList<QVariant> geoData;
	geoData << QVariant(geoName) << "" << "" << "";
	//BaseItem* grpItem = new BaseItem(listName, stations);
	GeoTreeItem* geo = new GeoTreeItem(geoData, _rootItem);
	_lists.push_back(geo);
	_rootItem->appendChild(geo);

	QList<QVariant> pointData;
	pointData << "Points" << "" << "" << "";
	//BaseItem* grpItem = new BaseItem(listName, stations);
	GeoObjectListItem* pointList = new GeoObjectListItem(pointData, geo, points, GEOLIB::POINT);
	geo->appendChild(pointList);


	size_t nPoints = points->size();

	for (size_t j=0; j<nPoints; j++)
	{
		QList<QVariant> pnt;
		pnt << static_cast<unsigned>(j) << QString::number((*(*points)[j])[0],'f') << QString::number((*(*points)[j])[1],'f') << QString::number((*(*points)[j])[2],'f');
		GeoTreeItem* point = new GeoTreeItem(pnt, pointList, static_cast<GEOLIB::Point*>((*points)[j]));
		pointList->appendChild(point);
	}

	std::cout << "Geometry \"" << geoName.toStdString() << "\" built." << std::endl;
	std::cout << nPoints << " points added." << std::endl;

	reset();
}

void GeoTreeModel::addPolylineList(QString geoName, const GEOLIB::PolylineVec* polylineVec)
{
	int nLists = _rootItem->childCount();
	TreeItem* geo(NULL);
	for (int i=0; i<nLists; i++)
	{
		if (_rootItem->child(i)->data(0).toString().compare(geoName) == 0)
			geo = _rootItem->child(i);
	}

	if (geo == NULL)
	{
		std::cout << "GeoTreeModel::addPolylineList() - Error: No corresponding geometry found..." << std::endl;
		return;
	}

	const std::vector<GEOLIB::Polyline*> *lines = polylineVec->getVector();

	QList<QVariant> plyData;
	plyData << "Polylines" << "" << "" << "";
	//BaseItem* grpItem = new BaseItem(listName, stations);
	GeoObjectListItem* plyList = new GeoObjectListItem(plyData, geo, lines, GEOLIB::POLYLINE);
	geo->appendChild(plyList);
	this->addChildren(plyList, polylineVec, 0, lines->size());
	reset();
}

void GeoTreeModel::appendPolylines(const std::string &name, const GEOLIB::PolylineVec* polylineVec)
{
	for (size_t i=0; i<_lists.size(); i++)
	{
		if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
		{
			for (int j=0; j<_lists[i]->childCount(); j++)
			{
				GeoObjectListItem* parent = static_cast<GeoObjectListItem*>(_lists[i]->child(j));
				if (GEOLIB::POLYLINE == parent->getType())
				{
					this->addChildren(parent, polylineVec, parent->childCount(), polylineVec->getVector()->size());
					reset();
					parent->vtkSource()->Modified();
					return;
				}
			}
		}
	}
	OGSError::box("Error adding polyline to geometry.");
}

void GeoTreeModel::addChildren(GeoObjectListItem* plyList, const GEOLIB::PolylineVec* polyline_vec, size_t start_index, size_t end_index)
{
	const std::vector<GEOLIB::Polyline*> *lines = polyline_vec->getVector();

	for (size_t i=start_index; i<end_index; i++)
	{
		QList<QVariant> line;
		std::string ply_name("");
		line << "Line " + QString::number(i);
		if (polyline_vec->getNameOfElementByID(i, ply_name)) line << QString::fromStdString(ply_name) << "" << "";
		else line << "" << "" << "";


		GeoTreeItem* lineItem = new GeoTreeItem(line, plyList, (*lines)[i]);
		plyList->appendChild(lineItem);

		int nPoints = static_cast<int>((*lines)[i]->getNumberOfPoints());
		for (int j=0; j<nPoints; j++)
		{
			QList<QVariant> pnt;
			pnt << static_cast<int>((*lines)[i]->getPointID(j)) << QString::number((*(*(*lines)[i])[j])[0],'f') << QString::number((*(*(*lines)[i])[j])[1],'f') << QString::number((*(*(*lines)[i])[j])[2],'f');
			TreeItem* child = new TreeItem(pnt, lineItem);
			lineItem->appendChild(child);
		}
	}
	std::cout << end_index-start_index << " polylines added." << std::endl;
}

void GeoTreeModel::addSurfaceList(QString geoName, const GEOLIB::SurfaceVec* surfaceVec)
{
	int nLists = _rootItem->childCount();
	TreeItem* geo(NULL);
	for (int i=0; i<nLists; i++)
	{
		if (_rootItem->child(i)->data(0).toString().compare(geoName) == 0)
			geo = _rootItem->child(i);
	}

	if (geo == NULL)
	{
		std::cout << "GeoTreeModel::addSurfaceList() - Error: No corresponding geometry found..." << std::endl;
		return;
	}

	const std::vector<GEOLIB::Surface*> *surfaces = surfaceVec->getVector();

	QList<QVariant> sfcData;
	sfcData << "Surfaces" << "" << "" << "";
	GeoObjectListItem* sfcList = new GeoObjectListItem(sfcData, geo, surfaces, GEOLIB::SURFACE);
	geo->appendChild(sfcList);
	this->addChildren(sfcList, surfaceVec, 0, surfaces->size());

	reset();
}

void GeoTreeModel::appendSurfaces(const std::string &name, GEOLIB::SurfaceVec* surfaceVec)
{
	for (size_t i=0; i<_lists.size(); i++)
	{
		if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
		{
			for (int j=0; j<_lists[i]->childCount(); j++)
			{
				GeoObjectListItem* parent = static_cast<GeoObjectListItem*>(_lists[i]->child(j));
				if (GEOLIB::SURFACE == parent->getType())
				{
					this->addChildren(parent, surfaceVec, parent->childCount(), surfaceVec->getVector()->size());
					reset();
					parent->vtkSource()->Modified();
					return;
				}
			}
		}
	}
	OGSError::box("Error adding surface to geometry.");
}

void GeoTreeModel::addChildren(GeoObjectListItem* sfcList, const GEOLIB::SurfaceVec* surface_vec, size_t start_index, size_t end_index)
{
	const std::vector<GEOLIB::Surface*> *surfaces = surface_vec->getVector();

	for (size_t i=start_index; i<end_index; i++)
	{
		QList<QVariant> surface;
		std::string sfc_name("");
		surface << "Surface " + QString::number(i);
		if (surface_vec->getNameOfElementByID(i, sfc_name)) surface << QString::fromStdString(sfc_name)  << "" << "";
		else  surface << "" << "" << "";

		GeoTreeItem* surfaceItem = new GeoTreeItem(surface, sfcList, (*surfaces)[i]);
		sfcList->appendChild(surfaceItem);

		const std::vector<GEOLIB::Point*> *nodesVec = (*surfaces)[i]->getPointVec();

		int nElems = static_cast<int>((*surfaces)[i]->getNTriangles());
		for (int j=0; j<nElems; j++)
		{
			QList<QVariant> elem;
			elem << j << static_cast<int>((*(*(*surfaces)[i])[j])[0]) << static_cast<int>((*(*(*surfaces)[i])[j])[1]) << static_cast<int>((*(*(*surfaces)[i])[j])[2]);
			TreeItem* child = new TreeItem(elem, surfaceItem);
			surfaceItem->appendChild(child);

			for (int k=0; k<3; k++)
			{
				QList<QVariant> node;
				node << static_cast<int>((*(*(*surfaces)[i])[j])[k]) << QString::number((*(*nodesVec)[(*(*(*surfaces)[i])[j])[k]])[0],'f') << QString::number((*(*nodesVec)[(*(*(*surfaces)[i])[j])[k]])[1],'f') << QString::number((*(*nodesVec)[(*(*(*surfaces)[i])[j])[k]])[2],'f');
				TreeItem* nchild = new TreeItem(node, child);
				child->appendChild(nchild);
			}
		}
	}
	std::cout << end_index-start_index << " surfaces added." << std::endl;
}

/**
 * Removes the TreeItem with the given name including all its children
 */
void GeoTreeModel::removeGeoList(const std::string &name, GEOLIB::GEOTYPE type)
{
	for (size_t i=0; i<_lists.size(); i++)
	{
		if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
		{
			for (int j=0; j<_lists[i]->childCount(); j++)
			{
				if (type == static_cast<GeoObjectListItem*>(_lists[i]->child(j))->getType())
				{
					QModelIndex index = createIndex(j, 0, _lists[i]->child(j));
					removeRows(0, _lists[i]->child(j)->childCount(), index);
					removeRows(j, 1, parent(index));
					break;
				}
			}
			if (_lists[i]->childCount() == 0)
			{
				_lists.erase(_lists.begin() + i);
				removeRows(i, 1, QModelIndex());
			}
		}
	}
}

vtkPolyDataAlgorithm* GeoTreeModel::vtkSource(const std::string &name, GEOLIB::GEOTYPE type) const
{
	size_t nLists = _lists.size();
	for (size_t i=0; i<nLists; i++)
	{
		if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
		{
			for (int j=0; j<_lists[i]->childCount(); j++)
			{
				GeoObjectListItem* item = dynamic_cast<GeoObjectListItem*>(_lists[i]->child(j));
				if (item->getType() == type)
					return item->vtkSource();
			}
		}
	}
	return NULL;
}
