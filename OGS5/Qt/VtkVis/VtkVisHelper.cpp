/**
 * \file VtkVisHelper.cpp
 * 22/09/2010 LB Initial implementation
 *
 * Implementation of VtkVisHelper class
 */

// ** INCLUDES **
#include "VtkVisHelper.h"

#include <vtkImageData.h>
#include <vtkTexture.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>

#include <QImage>

vtkImageData* VtkVisHelper::QImageToVtkImageData(QImage &img)
{
	size_t imgWidth = img.width(), imgHeight = img.height();
	vtkSmartPointer<vtkUnsignedCharArray> data = vtkSmartPointer<vtkUnsignedCharArray>::New();
	data->SetNumberOfComponents(3);
	data->SetNumberOfTuples( imgWidth*imgHeight );

	for (size_t j=0; j<imgHeight; j++) {
		for (size_t i=0; i<imgWidth; i++) {
			QRgb pix = img.pixel(i,j);
			const float color[3] = { qRed(pix), qGreen(pix), qBlue(pix) };
			data->SetTuple(j*imgWidth+i, color);
		}
	}

	vtkImageData* imgData = vtkImageData::New();
	imgData->SetExtent(0, imgWidth-1, 0, imgHeight-1, 0, 0);
	imgData->SetOrigin(0, 0, 0);
	imgData->SetNumberOfScalarComponents(3);
	imgData->GetPointData()->SetScalars(data);

	return imgData;
}

vtkTexture* VtkVisHelper::QImageToVtkTexture(QImage &img)
{
	
	vtkSmartPointer<vtkImageData> imgData = QImageToVtkImageData(img);

	vtkTexture* texture = vtkTexture::New();
	texture->InterpolateOff();
	texture->RepeatOff();
	//texture->EdgeClampOff();
	//texture->SetBlendingMode(0);
	texture->SetInput(imgData);

	return texture;
}
