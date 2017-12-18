/*
 * \file ProcessInfo.cpp
 *
 *  Created on: Sep 2, 2010
 *      Author: TF
 */

#include <ProcessInfo.h>
#include "rf_pcs.h"

ProcessInfo::ProcessInfo() :
_pcs_type (INVALID_PROCESS), _pcs_pv (INVALID_PV), _pcs (NULL)
{}

ProcessInfo::ProcessInfo (ProcessType pcs_type, PrimaryVariable pcs_pv, CRFProcess* pcs) :
_pcs_type (pcs_type), _pcs_pv (pcs_pv), _pcs (pcs)
{}

void ProcessInfo::setProcessType (ProcessType pcs_type)
{
   _pcs_type = pcs_type;
}


void ProcessInfo::setProcessPrimaryVariable (PrimaryVariable pcs_pv)
{
   _pcs_pv = pcs_pv;
}


void ProcessInfo::setProcess (CRFProcess* pcs)
{
   _pcs = pcs;
}


ProcessType ProcessInfo::getProcessType () const
{
   return _pcs_type;
}


PrimaryVariable ProcessInfo::getProcessPrimaryVariable () const
{
   return _pcs_pv;
}


CRFProcess* ProcessInfo::getProcess ()
{
   return _pcs;
}


ProcessInfo::~ProcessInfo()
{}
