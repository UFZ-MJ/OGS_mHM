/*
 * DistributionInfo.cpp
 *
 *  Created on: Sep 28, 2010
 *      Author: fischeth
 */

#include "DistributionInfo.h"

DistributionInfo::DistributionInfo() :
_dis_type (FiniteElement::INVALID_DIS_TYPE)
{}

DistributionInfo::~DistributionInfo()
{}

void DistributionInfo::setProcessDistributionType (FiniteElement::DistributionType dis_type)

{
   _dis_type = dis_type;
}


FiniteElement::DistributionType DistributionInfo::getProcessDistributionType () const
{
   return _dis_type;
}
