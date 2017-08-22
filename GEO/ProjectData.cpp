/**
 * \file ProjectData.cpp
 * 25/08/2010 KR Initial implementation
 */

#include "ProjectData.h"
#include "StringTools.h"


ProjectData::ProjectData()
	: _geoObjects ()
{}

ProjectData::~ProjectData()
{
	for (std::map<std::string, Mesh_Group::CFEMesh*>::iterator it = _msh_vec.begin(); it != _msh_vec.end(); ++it)
	{
		delete it->second;
		//_msh_vec.erase(it);
	}
}

void ProjectData::addMesh(Mesh_Group::CFEMesh* mesh, std::string &name)
{
	isUniqueMeshName(name);
	_msh_vec[name] = mesh;
};

const Mesh_Group::CFEMesh* ProjectData::getMesh(const std::string &name) const
{
	return _msh_vec.find(name)->second;
}

bool ProjectData::removeMesh(const std::string &name)
{
	delete _msh_vec[name];
	size_t result = _msh_vec.erase(name);
	return (result>0);
}

void ProjectData::addCondition(FEMCondition* cond)
{
	_cond_vec.push_back(cond);
};

const FEMCondition* ProjectData::getCondition(const std::string &name) const
{
	for (std::vector<FEMCondition*>::const_iterator it = _cond_vec.begin(); it != _cond_vec.end(); ++it)
	{
		if ((*it)->getGeoName().compare(name) == 0) return *it;
	}
	std::cout << "Error in ProjectData::getCondition() - No condition found with name \"" << name << "\"..." << std::endl;
	return NULL;
}

bool ProjectData::removeCondition(const std::string &name)
{
	for (std::vector<FEMCondition*>::iterator it = _cond_vec.begin(); it != _cond_vec.end(); ++it)
	{
		if ((*it)->getGeoName().compare(name) == 0)
		{
			_cond_vec.erase(it);
			return true;
		}
	}
	std::cout << "Error in ProjectData::getCondition() - No condition found with name \"" << name << "\"..." << std::endl;
	return false;
}

bool ProjectData::isUniqueMeshName(std::string &name)
{
	int count(0);
	bool isUnique(false);
	std::string cpName;

	while (!isUnique)
	{
		isUnique = true;
		cpName = name;

		count++;
		// If the original name already exists we start to add numbers to name for
		// as long as it takes to make the name unique.
		if (count>1) cpName = cpName + "-" + number2str(count);

		for (std::map<std::string, Mesh_Group::CFEMesh*>::iterator it = _msh_vec.begin(); it != _msh_vec.end(); ++it)
		{
			if ( cpName.compare(it->first) == 0 ) isUnique = false;
		}
	}

	// At this point cpName is a unique name and isUnique is true.
	// If cpName is not the original name, "name" is changed and isUnique is set to false,
	// indicating that a vector with the original name already exists.
	if (count>1)
	{
		isUnique = false;
		name = cpName;
	}
	return isUnique;
}
