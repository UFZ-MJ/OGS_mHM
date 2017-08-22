/**
 * \file ConvertVtkToOsg.cpp
 * 25/10/2010 LB Initial implementation
 * 
 * Implementation of ConvertVtkToOsg utility
 */

// ** INCLUDES **
#include "vtkOsgActor.h"

#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkImageMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkGeometryFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkImageDataGeometryFilter.h>

#include <OpenSG/OSGSceneFileHandler.h>

int main (int argc, char const *argv[])
{
	if (argc == 1)
	{
		std::cout << "Usage: ConvertVtkToOsg filename.vt*" << std::endl;
		return 1;
	}

	OSG::osgInit(0, NULL);

	std::string filename (argv[1]);
	std::cout << "Opening file " << filename << " ... " << std::endl << std::flush;
	size_t stringLength = filename.size();
	std::string osgFilename = filename.substr(0, stringLength - 3);
	osgFilename.append("osb");
	
	vtkXMLDataReader* reader = NULL;
	vtkGenericDataObjectReader* oldStyleReader = vtkGenericDataObjectReader::New();
	vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
	if (filename.find("vti", stringLength - 4))
	{
		reader = vtkXMLImageDataReader::New();
		vtkSmartPointer<vtkImageDataGeometryFilter> geoFilter = vtkSmartPointer<vtkImageDataGeometryFilter>::New();
		geoFilter->SetInputConnection(reader->GetOutputPort());
		mapper->SetInputConnection(geoFilter->GetOutputPort());
	}
	if (filename.find("vtr", stringLength - 4) != std::string::npos)
	{
		reader = vtkXMLRectilinearGridReader::New();
		vtkSmartPointer<vtkGeometryFilter> geoFilter = vtkSmartPointer<vtkGeometryFilter>::New();
		geoFilter->SetInputConnection(reader->GetOutputPort());
		mapper->SetInputConnection(geoFilter->GetOutputPort());
	}
	else if (filename.find("vts", stringLength - 4) != std::string::npos)
	{
		reader = vtkXMLStructuredGridReader::New();
		vtkSmartPointer<vtkGeometryFilter> geoFilter = vtkSmartPointer<vtkGeometryFilter>::New();
		geoFilter->SetInputConnection(reader->GetOutputPort());
		mapper->SetInputConnection(geoFilter->GetOutputPort());
	}
	else if (filename.find("vtp", stringLength - 4) != std::string::npos)
	{
		reader = vtkXMLPolyDataReader::New();
		mapper->SetInputConnection(reader->GetOutputPort());
	}
	else if (filename.find("vtu", stringLength - 4) != std::string::npos)
	{
		reader = vtkXMLUnstructuredGridReader::New();
		vtkSmartPointer<vtkGeometryFilter> geoFilter = vtkSmartPointer<vtkGeometryFilter>::New();
		geoFilter->SetInputConnection(reader->GetOutputPort());
		mapper->SetInputConnection(geoFilter->GetOutputPort());
	}
	else if (filename.find("vtk", stringLength - 4) != std::string::npos)
	{
		oldStyleReader->SetFileName(filename.c_str());
		oldStyleReader->Update();
		if(oldStyleReader->IsFilePolyData())
			mapper->SetInputConnection(oldStyleReader->GetOutputPort());
		else
		{
			vtkSmartPointer<vtkGeometryFilter> geoFilter = vtkSmartPointer<vtkGeometryFilter>::New();
			geoFilter->SetInputConnection(oldStyleReader->GetOutputPort());
			mapper->SetInputConnection(geoFilter->GetOutputPort());
		}
	}
	else
	{
		std::cout << "Not a valid vtk file ending (vti, vtr, vts, vtp, vtu, vtk)" << std::endl;
		return 1;
	}

    
	if (filename.find("vtk", stringLength - 4) == std::string::npos)
	{
		reader->SetFileName(filename.c_str());
		reader->Update();
	}

	vtkOsgActor* actor = vtkOsgActor::New();
	actor->SetVerbose(true);
	actor->SetMapper(mapper);
	actor->UpdateOsg();
	OSG::NodePtr root = actor->GetOsgRoot();
	OSG::SceneFileHandler::the().write(root, osgFilename.c_str());
	actor->ClearOsg();
	actor->Delete();
	
	//mapper->Delete(); // TODO crashes
	if (reader) reader->Delete();
	oldStyleReader->Delete();

	OSG::osgExit();
	
	std::cout << "File conversion finished" << std::endl;
	
	return 0;
}
