/*
 * ModifyMeshProperties.h
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 */

#ifndef MODIFYMESHPROPERTIES_H_
#define MODIFYMESHPROPERTIES_H_

#include "Polygon.h"

namespace Mesh_Group {

class CFEMesh;

class ModifyMeshProperties {
public:
	ModifyMeshProperties(CFEMesh* msh);
	virtual ~ModifyMeshProperties();

	void setMaterial (const GEOLIB::Polygon& polygon, size_t mat_id);
private:
	CFEMesh* _mesh;
};

}

#endif /* MODIFYMESHPROPERTIES_H_ */
