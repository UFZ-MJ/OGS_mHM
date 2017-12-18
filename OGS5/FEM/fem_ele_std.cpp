/*
   The members of class Element definitions.
*/

#include "Configure.h"

// C++ STL
#include <cfloat>
//#include <iostream>
//#include <limits>	// PCH to better use system max and min
// Method
#include "fem_ele_std.h"
#include "mathlib.h"
// Problems
//#include "rf_mfp_new.h"
#include "rf_mmp_new.h"

#include "pcs_dm.h"                               // displacement coupled
#include "rfmat_cp.h"
// Steps
//#include "rf_pcs.h"
//#include "rf_tim_new.h"
#ifndef NEW_EQS                                   //WW. 06.11.2008
// Sytem matrix
#include "matrix.h"
#endif
// Parallel computing
//#include "par_ddc.h"
// MSHLib
//#include "msh_elem.h"
// Solver
#ifdef NEW_EQS
#include "equation_class.h"
using Math_Group::CSparseMatrix;
#endif

extern double gravity_constant;                   // TEST, must be put in input file
#define COMP_MOL_MASS_AIR   28.96                 // kg/kmol WW  28.96
#define COMP_MOL_MASS_WATER 18.016                //WW 18.016
#define GAS_CONSTANT    8314.41                   // J/(kmol*K) WW
#define GAS_CONSTANT_V  461.5                     //WW
#define T_KILVIN_ZERO  273.15                     //WW

using namespace std;

namespace FiniteElement
{

   //========================================================================
   // Element calculation
   //========================================================================
   /**************************************************************************
      GeoSys - Function: Constructor
      Programmaenderungen:
      01/2005   WW    Erste Version
   **************************************************************************/
   CFiniteElementStd:: CFiniteElementStd(CRFProcess *Pcs, const int C_Sys_Flad, const int order)
      : CElement(C_Sys_Flad, order), phase(0), comp(0), SolidProp(NULL),
      FluidProp(NULL), MediaProp(NULL),
      pcs(Pcs), dm_pcs(NULL), HEAD_Flag(false)
   {
      int i;
      int size_m = 20;                            //25.2.2007
      string name2;
      char name1[MAX_ZEILE];
      cpl_pcs = NULL;
      //27.2.2007 WW
      newton_raphson = false;
                                                  //WW
      if (pcs->m_num->nls_method_name.compare("NEWTON_RAPHSON") == 0)
         newton_raphson = true;
      Mass = NULL;
      Mass2 = NULL;
      Laplace = NULL;
      Advection = NULL;
      Storage = NULL;
      Content = NULL;
      StrainCoupling = NULL;
      RHS = NULL;
      //
      NodalVal1 = new double[size_m];
      NodalVal2 = new double[size_m];
      NodalVal3 = new double[size_m];
      NodalVal4 = new double[size_m];
      NodalValC = new double[size_m];
      NodalValC1 = new double[size_m];
      NodalVal_Sat = new double[size_m];
      NodalVal_SatNW = new double[size_m];
      NodalVal_p2 = new double[size_m];
      NodalVal_p20 = new double [size_m];         //AKS
      //NW
      switch (C_Sys_Flad / 10)
      {
         case 1:
            weight_func = new double[2];
            break;
         case 2:
            weight_func = new double[4];
            break;
         case 3:
            weight_func = new double[8];
            break;
      }
      //
      // 27.2.2007. GravityMatrix = NULL;
      m_dom = NULL;
      eqs_rhs = NULL;                             //08.2006 WW
      //
      // 12.12.2007 WW
      for (i = 0; i < 4; i++)
         NodeShift[i] = 0;
      //
      dynamic = false;
      if (pcs->pcs_type_name_vector.size() && pcs->pcs_type_name_vector[0].find(
         "DYNAMIC") != string::npos)
         dynamic = true;

      dm_pcs = NULL;
      heat_phase_change = false;

      idx_vel_disp[0] = idx_vel_disp[1] = idx_vel_disp[2] = -1;
      idx_pres = -1;

      idxS = idx3 = -1;
      if (pcs->primary_variable_name.compare("HEAD") == 0)
         HEAD_Flag = true;
      //SB4218 added
      string pcs_primary = pcs->pcs_primary_function_name[0];
      if (pcs_primary.compare("HEAD") == 0)
         HEAD_Flag = true;
      for (i = 0; i < 9; i++)
         mat[i] = 0.0;

      idx0 = idx1 = 0;                            // column index in the node value data
      LocalShift = 0;
      //	pcsT = pcs->pcs_type_name[0]; // TF
      pcsT = (convertProcessTypeToString (pcs->getProcessType()))[0];
      //	if (pcs->pcs_type_name.find("AIR") != string::npos) //OK // TF commented out
      if (pcs->getProcessType () == AIR_FLOW)     //OK
         pcsT = 'A';
      //	if (pcs->pcs_type_name.find("MULTI") != string::npos) // 24.02.2007 WW
                                                  // 24.02.2007 WW
      if (pcs->getProcessType () == MULTI_PHASE_FLOW)
         pcsT = 'V';                              // Non-isothermal multiphase flow
      switch (pcsT)
      {
         default:
            PcsType = L;
            //WW GravityMatrix = new  SymMatrix(size_m);
            if (dynamic)
            {
               idx0 = pcs->GetNodeValueIndex("PRESSURE_RATE1");
               idx1 = idx0 + 1;
               idx_pres = pcs->GetNodeValueIndex("PRESSURE1");
               idx_vel_disp[0] = pcs->GetNodeValueIndex("VELOCITY_DM_X");
               idx_vel_disp[1] = pcs->GetNodeValueIndex("VELOCITY_DM_Y");
               if (dim == 3)
                  idx_vel_disp[2] = pcs->GetNodeValueIndex("VELOCITY_DM_Z");
            }
            else
            {
               idx0 = pcs->GetNodeValueIndex("PRESSURE1");
               idx1 = idx0 + 1;
            }
            break;
         case 'L':                                // Liquid flow
            PcsType = L;
            // 02.2.2007 GravityMatrix = new  SymMatrix(size_m);
            if (dynamic)
            {
               idx0 = pcs->GetNodeValueIndex("PRESSURE_RATE1");
               idx1 = idx0 + 1;
               idx_pres = pcs->GetNodeValueIndex("PRESSURE1");
               idx_vel_disp[0] = pcs->GetNodeValueIndex("VELOCITY_DM_X");
               idx_vel_disp[1] = pcs->GetNodeValueIndex("VELOCITY_DM_Y");
               if (dim == 3)
                  idx_vel_disp[2] = pcs->GetNodeValueIndex("VELOCITY_DM_Z");
            }
            else
            {
               idx0 = pcs->GetNodeValueIndex("PRESSURE1");
               idx1 = idx0 + 1;
            }
            idx_vel[0] = pcs->GetNodeValueIndex("VELOCITY_X1");
            idx_vel[1] = pcs->GetNodeValueIndex("VELOCITY_Y1");
            idx_vel[2] = pcs->GetNodeValueIndex("VELOCITY_Z1");
            break;
         case 'U':                                // Unconfined flow
            PcsType = U;
            break;
         case 'G':                                // Groundwater flow
            PcsType = G;
            idx0 = pcs->GetNodeValueIndex("HEAD");
            idx1 = idx0 + 1;
                                                  //WW
            idx_vel[0] =  pcs->GetNodeValueIndex("VELOCITY_X1");
                                                  //WW
            idx_vel[1] =  pcs->GetNodeValueIndex("VELOCITY_Y1");
                                                  //WW
            idx_vel[2] =  pcs->GetNodeValueIndex("VELOCITY_Z1");
            break;
         case 'T':                                // Two-phase flow
            PcsType = T;
            break;
         case 'C':                                // Componental flow
            PcsType = C;
            break;
         case 'H':                                // heat transport
            PcsType = H;
            idx0 = pcs->GetNodeValueIndex("TEMPERATURE1");
            idx1 = idx0 + 1;
            break;
         case 'M':                                // Mass transport
            PcsType = M;
            sprintf(name1, "%s", pcs->pcs_primary_function_name[0]);
            name2 = name1;
            idx0 = pcs->GetNodeValueIndex(name2);
            idx1 = idx0 + 1;
            break;
         case 'O':                                // Liquid flow
            PcsType = O;
            break;
         case 'R':                                //OK4104 Richards flow
            // 02.2.2007 GravityMatrix = new  SymMatrix(size_m);
            idx0 = pcs->GetNodeValueIndex("PRESSURE1");
            idx1 = idx0 + 1;
            idxS = pcs->GetNodeValueIndex("SATURATION1") + 1;
            idx_vel[0] = pcs->GetNodeValueIndex("VELOCITY_X1");
            idx_vel[1] = pcs->GetNodeValueIndex("VELOCITY_Y1");
            idx_vel[2] = pcs->GetNodeValueIndex("VELOCITY_Z1");
            if ((int) pcs->dof > 1)               //Dual porosity model. WW
            {
               idxp20 = pcs->GetNodeValueIndex("PRESSURE2");
               idxp21 = idxp20 + 1;
                                                  //WW
               Advection = new Matrix(size_m, size_m);
                                                  //12.12.2007 WW
               for (i = 0; i < pcs->pcs_number_of_primary_nvals; i++)
                  NodeShift[i] = i * pcs->m_msh->GetNodesNumber(false);
            }
            PcsType = R;
            break;
         case 'A':                                // Air (gas) flow
            PcsType = A;
                                                  //OK
            idx0 = pcs->GetNodeValueIndex("PRESSURE1");
            idx1 = idx0 + 1;                      //OK
            break;
         case 'F':                                // Fluid Momentum Process
            PcsType = R;                          // R should include L if the eqn of R is written right.
            break;
         case 'V':                                // 24.02.2007 WW
            // // 02.2.2007 GravityMatrix = new  SymMatrix(size_m);
                                                  //12.12.2007 WW
            for (i = 0; i < pcs->pcs_number_of_primary_nvals; i++)
               NodeShift[i] = i * pcs->m_msh->GetNodesNumber(false);
            //
            idx0 = pcs->GetNodeValueIndex("PRESSURE1");
            idx1 = idx0 + 1;
            idxp20 = pcs->GetNodeValueIndex("PRESSURE2");
            idxp21 = idxp20 + 1;
            idxS = pcs->GetNodeValueIndex("SATURATION1") + 1;
            idx_vel[0] = pcs->GetNodeValueIndex("VELOCITY_X1");
            idx_vel[1] = pcs->GetNodeValueIndex("VELOCITY_Y1");
            idx_vel[2] = pcs->GetNodeValueIndex("VELOCITY_Z1");
            PcsType = V;
            size_m = 40;
            break;
         case 'P':                                // 04.03.2009 PCH
            for (i = 0; i < pcs->pcs_number_of_primary_nvals; i++)
               NodeShift[i] = i * pcs->m_msh->GetNodesNumber(false);
            //
            idx0 = pcs->GetNodeValueIndex("PRESSURE1");
            idx1 = idx0 + 1;
            idxSn0 = pcs->GetNodeValueIndex("SATURATION2");
            idxSn1 = idxSn0 + 1;
            idxS = pcs->GetNodeValueIndex("SATURATION1") + 1;
            idx_vel[0] = pcs->GetNodeValueIndex("VELOCITY_X1");
            idx_vel[1] = pcs->GetNodeValueIndex("VELOCITY_Y1");
            idx_vel[2] = pcs->GetNodeValueIndex("VELOCITY_Z1");
            PcsType = P;
            size_m = 40;
            break;
      }
      if (pcs->Memory_Type == 0)                  // Do not store local matrices
      {
         if (PcsType == V || PcsType == P)        // 04.03.2009 PCH
            Mass2 = new Matrix(size_m, size_m);
         else
            Mass = new Matrix(size_m, size_m);
         Laplace = new Matrix(size_m, size_m);
         if (pcsT == 'H' || pcsT == 'M' || pcsT == 'A')
         {
            Advection = new Matrix(size_m, size_m);
            Storage = new Matrix(size_m, size_m);
            Content = new Matrix(size_m, size_m);
         }
         if (D_Flag)
            StrainCoupling = new Matrix(size_m, 60);
         RHS = new Vec(size_m);
      }
      //
      StiffMatrix = new Matrix(size_m, size_m);
      AuxMatrix = new Matrix(size_m, size_m);
      AuxMatrix1 = new Matrix(size_m, size_m);

      time_unit_factor = pcs->time_unit_factor;

      check_matrices = true;
      //
      SolidProp1 = NULL;
      MediaProp1 = NULL;
      flag_cpl_pcs = false;                       //OK
      // size_m changed
      NodalVal = new double[size_m];
      NodalVal0 = new double[size_m];
   }

   /**************************************************************************
      GeoSys - Function: Destructor
      Programmaenderungen:
      01/2005   WW    Erste Version
   **************************************************************************/
   // Destructor
   CFiniteElementStd::~CFiniteElementStd()
   {
      //  02.2.2007 if(GravityMatrix) delete GravityMatrix;
      // 02.2.2007  GravityMatrix = NULL;

      if(pcs->Memory_Type==0)                     // Do not store local matrices
      {
         if(Mass) delete Mass;
         if(Mass2) delete Mass2;
         if(Laplace) delete Laplace;
         if(Advection) delete Advection;
         if(Storage) delete Storage;
         if(Content) delete Content;
         if(StrainCoupling) delete StrainCoupling;
         if(RHS) delete RHS;
         Mass = NULL;
         Laplace = NULL;
         Advection = NULL;
         Storage = NULL;
         Content = NULL;
         StrainCoupling = NULL;
         RHS = NULL;
      }

      delete StiffMatrix;
      delete AuxMatrix;
      delete AuxMatrix1;

      StiffMatrix = NULL;
      AuxMatrix = NULL;
      AuxMatrix1 = NULL;
      // 27.2.2007 WW
      delete [] NodalVal;
      delete [] NodalVal0;
      delete [] NodalVal1;
      delete [] NodalVal2;
      delete [] NodalVal3;
      delete [] NodalVal4;
      delete [] NodalValC;
      delete [] NodalValC1;
      delete [] NodalVal_Sat;
      delete [] NodalVal_SatNW;
      delete [] NodalVal_p2;
      delete [] NodalVal_p20;                     //AKS
      //NW
      delete [] weight_func;
      weight_func = NULL;
   }
   /**************************************************************************
      GeoSys - Function: SetMemory

      Aufgabe:
            Set memory for local matrices
      Programmaenderungen:
      01/2005   WW    Erste Version

   **************************************************************************/
   void CFiniteElementStd::SetMemory()
   {
      int Size=nnodes;
      if(PcsType==V || PcsType==P)                //4.3.2009 PCH
         Size *= 2;
      ElementMatrix * EleMat = NULL;
      // Prepare local matrices
      // If local matrices are not stored, resize the matrix
      if(pcs->Memory_Type==0)
      {
         if(PcsType==V || PcsType==P)             //04.3.2009 PCH
            Mass2->LimitSize(Size, Size);
         else
            Mass->LimitSize(nnodes, nnodes);      // Mass->LimitSize(nnodes); // unsymmetric in case of Upwinding
         Laplace->LimitSize(Size, Size);
         if(PcsType==H||PcsType==M ||PcsType==A)
         {
            Advection->LimitSize(Size, Size);     //SB4200
            Storage->LimitSize(Size, Size);       //SB4200
            Content->LimitSize(Size, Size);       //SB4209
         }
         if(PcsType==R&&pcs->type==22)            //dual-porosity. WW
            Advection->LimitSize(Size, Size);
         if(D_Flag>0)
            StrainCoupling->LimitSize(Size, dim*nnodesHQ);
         RHS->LimitSize(Size);
      }
      else
      {
         EleMat = pcs->Ele_Matrices[Index];
         // if(PcsType==V) //24.2.2007 WW
         // Mass2 = EleMat->GetMass2();
         Mass = EleMat->GetMass();
         Laplace = EleMat->GetLaplace();
         // Advection, Storage, Content SB4200
         if(PcsType==M)
         {
            Advection = EleMat->GetAdvection();
            Storage = EleMat->GetStorage();
            Content = EleMat->GetContent();
         }
         RHS = EleMat->GetRHS();
         if(D_Flag>0)
            StrainCoupling = EleMat->GetCouplingMatrixB();
         if(D_Flag==41) LocalShift = dim*nnodesHQ;
      }

      //25.2.2007.WW if(GravityMatrix) GravityMatrix->LimitSize(nnodes);

      StiffMatrix->LimitSize(Size, Size);
      AuxMatrix->LimitSize(Size, Size);
      AuxMatrix1->LimitSize(Size, Size);
   }

   /**************************************************************************
      GeoSys - Function: ConfigureCoupling

      Aufgabe:
            Set coupling information for local fem calculation
      Programmaenderungen:
      01/2005   WW    Erste Version
      02/2007   WW    Multi phase flow
       03/2009   PCH	 PS_GLOBAL

   **************************************************************************/
   void  CFiniteElementStd::ConfigureCoupling(CRFProcess* pcs, const int *Shift, bool dyn)
   {
      char pcsT;
      //  pcsT = pcs->pcs_type_name[0]; // TF
                                                  // TF
      ProcessType pcs_type (pcs->getProcessType());
      pcsT = convertProcessTypeToString(pcs_type)[0];

      //  if(pcs->pcs_type_name.find("AIR")!=string::npos) //OK
      if(pcs_type == AIR_FLOW)                    //OK
         pcsT = 'A';
      if(pcs_type == MULTI_PHASE_FLOW)            //24.2.2007 WW
         pcsT = 'V';

      if (D_Flag > 0)
      {
         for (size_t i = 0; i < pcs_vector.size(); i++)
         {
            //			if(pcs_vector[i]->pcs_type_name.find("DEFORMATION")!=string::npos){
            if (isDeformationProcess(pcs_vector[i]->getProcessType()))
            {
               dm_pcs = (CRFProcessDeformation*) pcs_vector[i];
               break;
            }
         }
         if (dyn)
         {
            Idx_dm0[0] = dm_pcs->GetNodeValueIndex("ACCELERATION_X1");
            Idx_dm0[1] = dm_pcs->GetNodeValueIndex("ACCELERATION_Y1");
         }
         else
         {
            Idx_dm0[0] = dm_pcs->GetNodeValueIndex("DISPLACEMENT_X1");
            Idx_dm0[1] = dm_pcs->GetNodeValueIndex("DISPLACEMENT_Y1");
         }
         Idx_dm1[0] = Idx_dm0[0] + 1;
         Idx_dm1[1] = Idx_dm0[1] + 1;
         //     if(problem_dimension_dm==3)
         if (dim == 3)
         {
            if (dyn)
               Idx_dm0[2] = dm_pcs->GetNodeValueIndex("ACCELERATION_Z1");
            else
               Idx_dm0[2] = dm_pcs->GetNodeValueIndex("DISPLACEMENT_Z1");
            Idx_dm1[2] = Idx_dm0[2] + 1;
         }
         if (dm_pcs->type == 41)
         {
            for (size_t i = 0; i < 4; i++)
               NodeShift[i] = Shift[i];
         }
      }

      switch(pcsT)
      {
         default:
            if(T_Flag)
            {
               cpl_pcs = PCSGet("HEAT_TRANSPORT");
               idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
               idx_c1 = idx_c0+1;
            }
            break;
         case 'L':                                // Liquid flow
            if(T_Flag)
            {
               cpl_pcs = PCSGet("HEAT_TRANSPORT");
               idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
               idx_c1 = idx_c0+1;
            }
            break;
         case 'U':                                // Unconfined flow
            break;
         case 'G':                                // Groundwater flow
            if(T_Flag)
            {
               cpl_pcs = PCSGet("HEAT_TRANSPORT");
               idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
               idx_c1 = idx_c0+1;
            }
            break;
         case 'T':                                // Two-phase flow
            if(pcs->pcs_type_number==0)
            {
               cpl_pcs = pcs_vector[pcs->pcs_number+1];
               idx_c0 = cpl_pcs->GetNodeValueIndex("SATURATION2");
               idx_c1 = idx_c0+1;
            }
            else if(pcs->pcs_type_number==1)
            {
               cpl_pcs = pcs_vector[pcs->pcs_number-1];
               idx_c0 = cpl_pcs->GetNodeValueIndex("PRESSURE1");
               idx_c1 = idx_c0+1;
            }
            break;
         case 'C':                                // Componental flow
            break;
         case 'H':                                // heat transport
            //SB CMCD this needs to be fixed
            cpl_pcs = PCSGet("GROUNDWATER_FLOW");
            if(cpl_pcs)                           //WW
            {
               idx_c0 = cpl_pcs->GetNodeValueIndex("HEAD");
               idx_c1 = idx_c0+1;
            }
            else
            {
               cpl_pcs = PCSGet("LIQUID_FLOW");
               if(cpl_pcs == NULL)
               {
                                                  //OK
                  cpl_pcs = PCSGet("RICHARDS_FLOW");
                  if(cpl_pcs)
                                                  //WW
                     idxS = cpl_pcs->GetNodeValueIndex("SATURATION1")+1;
               }
               if(cpl_pcs == NULL)
               {
                                                  //24.042.2004 WW
                  cpl_pcs = PCSGet("MULTI_PHASE_FLOW");
                  if(cpl_pcs)
                                                  //WW
                     idxS = cpl_pcs->GetNodeValueIndex("SATURATION1")+1;
               }
               if(cpl_pcs == NULL)                //23.02.2009 NB 4.9.05
               {
                  cpl_pcs = PCSGet("TWO_PHASE_FLOW");
                  if(cpl_pcs)
                     idxS = cpl_pcs->GetNodeValueIndex("SATURATION1")+1;
               }
               if(cpl_pcs == NULL)                //23.02.2009 NB 4.9.05
               {
                  cpl_pcs = PCSGet("AIR_FLOW");   //23.01.2009 NB
               }

               if (cpl_pcs)                       //MX
               {
                  idx_c0 = cpl_pcs->GetNodeValueIndex("PRESSURE1");
                  idx_c1 = idx_c0+1;
               }
            }
            break;
         case 'M':                                // Mass transport
            if(T_Flag)
            {
               cpl_pcs = PCSGet("HEAT_TRANSPORT");
               idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
               idx_c1 = idx_c0+1;
            }
            break;
         case 'O':                                // Liquid flow
            break;
         case 'R':                                // Richards flow
            if(T_Flag)                            //if(PCSGet("HEAT_TRANSPORT"))
            {
               cpl_pcs = PCSGet("HEAT_TRANSPORT");
               idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
               idx_c1 = idx_c0+1;
            }
            break;
         case 'V':                                // Multi-phase flow. 24.2.2007 WW
            if(T_Flag)                            //if(PCSGet("HEAT_TRANSPORT"))
            {
               cpl_pcs = PCSGet("HEAT_TRANSPORT");
               idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
               idx_c1 = idx_c0+1;
            }
            break;
         case 'A':                                //Gas flow
            if(T_Flag)                            //NB 23.01.2009 4.9.05
            {
               cpl_pcs = PCSGet("HEAT_TRANSPORT");
               idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
               idx_c1 = idx_c0+1;
            }
            break;
         case 'P':                                // Multi-phase flow. 03.03.2009 PCH
            if(T_Flag)
            {
               cpl_pcs = PCSGet("HEAT_TRANSPORT");
               idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
               idx_c1 = idx_c0+1;
            }
            break;

      }
   }

   /*************************************************************************
   FEMLib-Function:
   Task: Set material pointers to the current element
   01/2005 WW Implementation
   03/2005 OK MultiMSH
   11/2005 YD Set cursor of gas
   06/2009 OK MMP test not here (time consuming)
   01/2010 NW Set geo_area here
   **************************************************************************/
   void CFiniteElementStd::SetMaterial(int phase)
   {
      phase = 0;
      //----------------------------------------------------------------------
      // MMP
      int mmp_index=0;
      long group = MeshElement->GetPatchIndex();
      mmp_index = group;
      // Single continua thermal:
      if(msp_vector.size()>0)
      {
         SolidProp = msp_vector[mmp_index];
         SolidProp->m_pcs = pcs;                  //NW
         SolidProp->Fem_Ele_Std = this;           //CMCD for Decovalex
      }

      if(pcs->type==22)                           //WW/YD
      {
         if(pcs->GetContinnumType()== 0)          // Matrix //WW
            mmp_index = 2*group;
         else                                     // fracture //WW
            mmp_index = 2*group+1;
      }
      MediaProp = mmp_vector[mmp_index];
      MediaProp->m_pcs = pcs;
      MediaProp->Fem_Ele_Std = this;
      MeshElement->area = MediaProp->geo_area;    // NW
      //----------------------------------------------------------------------
      // MSP
      // If dual thermal:
      /*
      if(msp_vector.size()>0)
      {
        SolidProp = msp_vector[mmp_index];
        SolidProp->Fem_Ele_Std = this;//CMCD for Decovalex
      }
      */
      if(pcs->type==22)                           //WW
      {
         if(pcs->GetContinnumType()== 0)          // Matrix //WW
            mmp_index = 2*group+1;
         else                                     // fracture //WW
            mmp_index = 2*group;
         MediaProp1 = mmp_vector[mmp_index];
         MediaProp1->m_pcs = pcs;
         MediaProp1->Fem_Ele_Std = this;
         //----------------------------------------------------------------------
         // MSP
         // If dual thermal:
         /*
         if(msp_vector.size()>0)
         {
           SolidProp1 = msp_vector[mmp_index];
           SolidProp1->Fem_Ele_Std = this;//CMCD for Decovalex
         }
         */
      }
      //----------------------------------------------------------------------
      // MFP
      /* Comment out - NW
        if(PCSGet("LIQUID_FLOW")){
          FluidProp = MFPGet("LIQUID");
          if(!FluidProp)
            cout << "Warning: LIQUID was not found in fluid properties." << endl;
        }
      */
      if(mfp_vector.size()>0)
      {
         FluidProp = mfp_vector[0];
         FluidProp->Fem_Ele_Std = this;
      }
                                                  // 03.2009 PCH
      if((PCSGet("RICHARDS_FLOW")&&PCSGet("HEAT_TRANSPORT"))||pcs->type==1212||pcs->type==1313)
      {
         FluidProp = MFPGet("LIQUID");
         FluidProp->Fem_Ele_Std = this;
         //FluidProp = mfp_vector[0];
         GasProp = MFPGet("GAS");
         if (GasProp) GasProp->Fem_Ele_Std = this;
      }
      //----------------------------------------------------------------------
      // MCP
      //----------------------------------------------------------------------
   }

   /*************************************************************************
   FEMLib-Function:
   Task: Line element integration data for CVFEM overland flow
         to move
   Programming:
        6/2007 : JOD
   **************************************************************************/
   void CFiniteElementStd::GetOverlandBasisFunctionMatrix_Line()
   {

      edlluse[0] = 1.0;
      edlluse[1] = -1.0;
      edlluse[2] = -1.0;
      edlluse[3] = 1.0;

      edttuse[0] = 0.0;
      edttuse[1] = 0.0;
      edttuse[2] = 0.0;
      edttuse[3] = 0.0;
      ////MB nur Zeitweise hier
   }
   /*************************************************************************
   FEMLib-Function:
   Task: Quad element integration data for CVFEM overland flow
         to move
   Programming:
            ?    MB
   **************************************************************************/
   void CFiniteElementStd::GetOverlandBasisFunctionMatrix_Quad()
   {

      edlluse[0] = 0.5;
      edlluse[1] = -0.5;
      edlluse[2] = 0.0;
      edlluse[3] = 0.0;
      edlluse[4] = -0.5;
      edlluse[5] = 0.5;
      edlluse[6] = 0.;
      edlluse[7] = 0.;
      edlluse[8] = 0.;
      edlluse[9] = 0.;
      edlluse[10] = 0.5;
      edlluse[11] = -0.5;
      edlluse[12] = 0.;
      edlluse[13] = 0.;
      edlluse[14] = -0.5;
      edlluse[15] = 0.5;

      edttuse[0] = 0.5;
      edttuse[1] = 0.;
      edttuse[2] = 0.;
      edttuse[3] = -0.5;
      edttuse[4] = 0.;
      edttuse[5] = 0.5;
      edttuse[6] = -0.5;
      edttuse[7] = 0.;
      edttuse[8] = 0.;
      edttuse[9] = -0.5;
      edttuse[10] = 0.5;
      edttuse[11] = 0.;
      edttuse[12] = -0.5;
      edttuse[13] = 0.;
      edttuse[14] = 0.;
      edttuse[15] = 0.5;

   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculates consitutive relationships for CVFEM Overland Flow -> swval, swold
         for surface structure
   Programing:
   06/2005 MB Implementation
   04/2007 JOD modifications
   **************************************************************************/
   void CFiniteElementStd::CalcOverlandNLTERMS(double* haa, double* haaOld, double* swval, double* swold)
   {

      if(MediaProp->channel == 1)
         CalcOverlandNLTERMSChannel(haa, haaOld, swval, swold);
      else
         CalcOverlandNLTERMSRills(haa, haaOld, swval, swold);

   }
   /**************************************************************************
   FEMLib-Method:
   Task: Calculates consitutive relationships for CVFEM Overland Flow -> swval, swold
         for surface structure
   Programing:
   06/2007 JOD implementation
   **************************************************************************/
   void CFiniteElementStd::CalcOverlandNLTERMSRills(double* haa, double* haaOld, double* swval, double* swold)
   {
      double WDepth[4], WDepthOld[4];
      double rill_height = MediaProp->rill_height;
      double eps = MediaProp->rill_epsilon;

      for(int i=0; i<nnodes; i++)
      {
         WDepth[i] = haa[i] - Z[i];
         WDepthOld[i] = haaOld[i] - Z[i];
         if (MediaProp->rill_epsilon > 0)
         {
            if(WDepth[i] > 0)
               swval[i] = (WDepth[i] + eps)*(WDepth[i] + eps) / ( WDepth[i] + rill_height + eps) - pow(eps, 2.) / (rill_height + eps);
            else
               swval[i] = 0;

            if(WDepthOld[i] > 0)
                                                  // JOD
               swold[i] = (WDepthOld[i] + eps)*(WDepthOld[i] + eps) / ( WDepthOld[i] + rill_height + eps) - pow(eps, 2.) / (rill_height + eps);
            else
               swold[i] = 0;
         }                                        // end epsilon > 0
         else
         {
            swval[i] = WDepth[i];
            swold[i] = WDepthOld[i];
         }
      }

   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculates consitutive relationships for CVFEM Overland Flow -> swval, swold
         for channel
   Programing:
   06/2007 JOD Implementation
   **************************************************************************/
   void CFiniteElementStd::CalcOverlandNLTERMSChannel(double* haa, double* haaOld, double* swval, double* swold)
   {
      double WDepth[4], WDepthOld[4];
      double eps = MediaProp->rill_epsilon;
      double ratio;
      double xxx;

      for(int i=0; i<2; i++)
      {
         WDepth[i] = haa[i] - Z[i];
         WDepthOld[i] = haaOld[i] - Z[i];
         if (eps > 0)
         {
            ratio = WDepth[i] / eps;
            if (ratio > 1.0)
               swval[i] = WDepth[i];
            else if (ratio > 0.0)
            {
               xxx = 2.0 * (1.0 - ratio);
               swval[i] = WDepth[i] * pow(ratio,xxx);
            }
            else
               swval[i] = 0.0;
            ////////////////////////

            ratio = WDepthOld[i] / eps;
            if (ratio > 1.0)
               swold[i] = WDepthOld[i];
            else if (ratio > 0.0)
            {
               xxx = 2.0 * (1.0 - ratio);
               swold[i] =  WDepthOld[i] * pow(ratio,xxx);
            }
            else
               swold[i] = 0.0;
         }                                        // end epsilon > 0
         else
         {
            swval[i] = WDepth[i];
            swold[i] = WDepthOld[i];
         }
      }                                           //end for

   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculates upstream weighting for CVFEM Overland Flow -> ckwr and iups
   Programing:
   06/2005 MB Implementation
   04/2007 JOD modifications
   **************************************************************************/
   void CFiniteElementStd::CalcOverlandCKWR(double* head, double* ckwr, int* iups)
   {

      double width = MediaProp->overland_width;
      double depth_exp = MediaProp->friction_exp_depth;
      double rill_depth = MediaProp->rill_height;
      int i, j;
      double maxZ;
      double flow_depth;

      for (i = 0; i < nnodes; i++)
      {
         for (j = 0; j < nnodes; j++)
         {
            maxZ = MMax(Z[i],Z[j]);
            if( head[i] > head[j] )
            {
               iups[i*nnodes + j] = i;
               flow_depth = head[i] - maxZ - rill_depth;
            }
            else
            {
               iups[i*nnodes + j]= j;
               flow_depth = head[j] - maxZ - rill_depth;
            }
            ////////////////////////////////////////
            if(flow_depth<0.0)
               ckwr[i*nnodes + j] = 0.0;
            else
            {
               if(MediaProp->channel == 1)
                  ckwr[i*nnodes + j] = flow_depth * pow(flow_depth * width / (2 * flow_depth + width),depth_exp);
               else
                  ckwr[i*nnodes + j] = pow(flow_depth,depth_exp + 1);
            }
         }                                        //end for j
      }                                           //end for i

   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculates upstream weighting for CVFEM Overland Flow -> ckwr and iups
         at node (i,j)
        used in AssemleParabolicEquationNewtonJacobi()
   Programing:
   06/2005 MB Implementation
   04/2007 JOD modifications
   **************************************************************************/
   void CFiniteElementStd::CalcOverlandCKWRatNodes(int i, int j, double* head, double* ckwr, int* iups)
   {
      double width = MediaProp->overland_width;
      double depth_exp = MediaProp->friction_exp_depth;
      double rill_depth = MediaProp->rill_height;
      double flow_depth;
      double maxZ;

      maxZ = MMax(Z[i],Z[j]);
      if(iups[i*nnodes+j] == i)
         flow_depth = head[i] - maxZ - rill_depth;
      else
         flow_depth = head[j] - maxZ - rill_depth;
      ///////////////////////////////////////
      if(flow_depth < 0.0)
         *ckwr = 0;
      else
      {
         if(MediaProp->channel == 1)
            *ckwr = flow_depth * pow(flow_depth * width / ( 2 * flow_depth + width), depth_exp);
         else
            *ckwr = pow(flow_depth, depth_exp + 1);
      }

   }
   /**************************************************************************
   FEMLib-Method:
   Task: calculate upwinded diffusion matric coefficient for CVFEM
         used in AssemleParabolicEquationNewton()
                AssemleParabolicEquationNewtonJacobi()
   Programing:
   06/2007 JOD Implementation
   **************************************************************************/
   void CFiniteElementStd::CalcOverlandUpwindedCoefficients(double** amat, double* ckwr, double axx, double ayy)
   {
      //double** amat;
      double gammaij;

      //amat = (double**) Malloc(nnodes * sizeof(double));
      //for (int i = 0; i < nnodes; i++)
      //  amat[i] = (double*) Malloc(nnodes*sizeof(double));

      //for (int i = 0; i < nnodes; i++)
      //  for (int j = 0; j < nnodes; j++)
      //    amat[i][j]= 0.0;

      for (int i = 0; i < nnodes; i++)
      {
         for (int j = (i+1); j < nnodes; j++)
         {
            gammaij = ckwr[i*nnodes+j] * ((edlluse[i*nnodes+j] * axx) + (edttuse[i*nnodes+j]* ayy));
            amat[i][j]= gammaij;
            amat[j][i]= gammaij;
            amat[i][i]= amat[i][i] - gammaij;
            amat[j][j]= amat[j][j] - gammaij;
         }
      }

      //return amat;

   }
   /**************************************************************************
   FEMLib-Method:
   Task: residual vector for overland CVFEM
         used in AssemleParabolicEquationNewton()
   Programing:
   06/2007 JOD Implementation
   **************************************************************************/
   void CFiniteElementStd::CalcOverlandResidual(double* head, double* swval, double* swold, double ast, double* residual, double** amat)
   {
      double sum;
      double storinit[4], astor[4], rhs[4];

      MNulleVec(astor,4);
      MNulleVec(rhs,nnodes);

      for (int i = 0; i < nnodes; i++)            // storage term
         rhs[i] = -ast * (swval[i] - swold[i]);

      /* if(MediaProp->channel ==1){ // channel, JOD removed, don't know what it was for
          astor[0] = swval[0] * ast;
          astor[1] = swval[1] * ast;
         rhs[0] = swold[0] * ast * HaaOld[0]; // swval ?????
          rhs[1] = swold[1] * ast * HaaOld[1]; // swval ?????
        }
      */
      //Form the residual excluding the right hand side vector

      for(int i = 0; i < nnodes; i++)
      {
         sum = 0.0;
         for(int j = 0; j < nnodes; j++)
         {
            sum = sum + ( amat[i][j]* head[j] );
         }
                                                  // astor = 0, rillDepth??
         storinit[i] = -rhs[i] + astor[i]* (head[i] - Z[i]);
         residual[i] = sum + storinit[i];
      }

   }
   /**************************************************************************
   FEMLib-Method:
   Task: calcukate jacobi overland CVFEM
         used in  AssemleParabolicEquationNewtonJacobi()
   Programing:
   06/2007 JOD Implementation
   **************************************************************************/
   double CFiniteElementStd::CalcOverlandJacobiNodes(int i, int j, double *head, double *headKeep, double akrw, double axx, double ayy, double** amat, double* sumjac )
   {
      double jacobi, gammaij, amatEps, amatKeep;

      gammaij= akrw *( axx*edlluse[i*nnodes+j] + ayy*edttuse[i*nnodes+j] );
      amatEps = gammaij * (head[j] - head[i]);
      amatKeep = amat[i][j] * (headKeep[j] - headKeep[i]);
      jacobi = -(amatEps-amatKeep);

      *sumjac = *sumjac + amatEps;

      return jacobi;

   }

   /**************************************************************************
   FEMLib-Method:
   Task: calculate topology coefficients for overland CVFEM
   Programing:
   08/2006 JOD Implementation
   **************************************************************************/
   void CFiniteElementStd::CalcOverlandCoefficients(double* head, double* axx, double* ayy, double* ast )
   {
      if(MeshElement->geo_type==1)
      {
         CalcOverlandCoefficientsLine(head, axx, ast );
         ayy = 0;
      }
      else if(MeshElement->geo_type==2)
         CalcOverlandCoefficientsQuad(head, axx, ayy, ast );
      else if(MeshElement->geo_type==4)
         CalcOverlandCoefficientsTri(head, axx, ayy, ast );
      else
         std::cout << "Error in CFiniteElementStd::CalcOverlandCoefficients !!!";

   }
   /**************************************************************************
   FEMLib-Method:
   Task:  calculate topology coefficientsfor overland CVFEM, line elements
   Programing:
   08/2006 JOD Implementation
   **************************************************************************/
   void CFiniteElementStd::CalcOverlandCoefficientsLine(double* head, double* axx, double* ast )
   {

      double dx, dy, dzx;
      double delt, dhds;
      double fric, width, eslope, slope_exp;

      fric =  MediaProp->friction_coefficient;
      slope_exp = MediaProp->friction_exp_slope;
      width = MediaProp->overland_width;

      dx = X[1] - X[0];
      dy = Y[1] - Y[0];
      dzx = Z[1] - Z[0];
      delt = sqrt(dx*dx + dy*dy);
      dhds = fabs( (head[0] - head[1]) / delt );

      GetOverlandBasisFunctionMatrix_Line();

      dhds = MMax(1.0e-10,dhds);
      eslope = 1.0 / dhds;
      eslope = pow(eslope, 1 - slope_exp);

      *axx = eslope * fric * width / delt;
      *ast = delt * width /(double) (nnodes * dt);

   }
   /**************************************************************************
   FEMLib-Method:
   Task:  calculate topology coefficientsfor overland CVFEM, rectangles
   Programing:
   08/2006 JOD Implementation
   **************************************************************************/
   void CFiniteElementStd::CalcOverlandCoefficientsQuad(double* head, double* axx, double* ayy, double* ast )
   {

      double dx, dy, dzx, dzy;
      double delt;
      double dhds, GradH[2];
      double fric, eslope, slope_exp;

      fric =  MediaProp->friction_coefficient;
      slope_exp = MediaProp->friction_exp_slope;

      /////////////////////////////
      dx = X[1] - X[0];                           //ell
      dy = Y[3] - Y[0];                           //ett
      dzx = Z[1] - Z[0];
      dzy = Z[3] - Z[0];
      dx = sqrt(dx*dx + dzx*dzx);
      dy = sqrt(dy*dy + dzy*dzy);
      delt = dx * dy;

      GetOverlandBasisFunctionMatrix_Quad();

      GradH[0] = (head[0] - head[1] - head[2] + head[3]) / (2.0*dx);
      GradH[1] = (head[0] + head[1] - head[2] - head[3]) / (2.0*dy) ;
                                                  // dh/ds (dh in the direction of maximum slope)
      dhds = sqrt((GradH[0] * GradH[0]) + (GradH[1] * GradH[1]));
      dhds = MMax(1.0e-10,dhds);
      eslope = 1.0 / dhds;
      eslope = pow(eslope, 1 - slope_exp);

      *axx = eslope * fric * dy/dx;               //ett/ell
      *ayy = eslope * fric * dx/dy;
      *ast = delt /(double) (nnodes * dt );

   }
   /**************************************************************************
   FEMLib-Method:
   Task:  calculate topology coefficientsfor overland CVFEM, triangles
   Programing:
   08/2006 JOD Implementation
   **************************************************************************/
   void CFiniteElementStd::CalcOverlandCoefficientsTri(double* head, double* axx, double* ayy, double* ast )
   {

      double x2, x3, y2, y3;
      double delt, delt2, delt2inv, b[3], g[3];
      double dhds, GradH[2];
      double fric, eslope, slope_exp;

      fric =  MediaProp->friction_coefficient;
      slope_exp = MediaProp->friction_exp_slope;

      x2 = X[1] - X[0];
      x3 = X[2] - X[0];
      y2 = Y[1] - Y[0];
      y3 = Y[2] - Y[0];
      delt = (x2*y3 - x3*y2) * 0.5;
      delt2 = 2.0 * delt;
      delt2inv = 1.0 / delt2;

      /////////////////////  GetOverlandBasisFunctionMatrix_Tri()
      b[0] = (y2-y3) * delt2inv;
      b[1] = y3 * delt2inv;
      b[2] = -y2 * delt2inv;
      g[0] = (x3-x2) * delt2inv;
      g[1] = -x3 * delt2inv;
      g[2] = x2 * delt2inv;

      for(int i=0; i<nnodes; i++)
         for(int j=0; j<nnodes; j++)
      {
         edlluse[i*nnodes + j] = b[i] * b[j];
         edttuse[i*nnodes + j] = g[i] * g[j];
      }
      //////////////////////////

      GradH[0] = ( b[0]*head[0] + b[1]*head[1] +  b[2]*head[2] );
      GradH[1] =( g[0]*head[0] + g[1]*head[1] +  g[2]*head[2] );
                                                  // dh/ds (dh in the direction of maximum slope)
      dhds = sqrt((GradH[0] * GradH[0]) + (GradH[1] * GradH[1]));
      dhds = MMax(1.0e-10,dhds);
      eslope = 1.0 / dhds;

      eslope = pow(eslope, 1 - slope_exp);
      *axx = eslope * fric * delt;
      *ayy = eslope * fric * delt;
      *ast = delt /(double) (nnodes * dt );
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate nodal enthalpy
   Programming: WW 09/2005
   **************************************************************************/
   inline void CFiniteElementStd::CalNodalEnthalpy()
   {
      int i;
      double temp, dT;
      for(i=0;i<nnodes;i++)
      {
         heat_phase_change =
            SolidProp->CheckTemperature_in_PhaseChange(NodalVal0[i], NodalVal1[i]);
         if(heat_phase_change) break;
      }
      if(!heat_phase_change) return;
      // Calculate node enthalpy
      for(i=0;i<nnodes;i++)
      {
         NodalVal_Sat[i] = pcs->GetNodeValue(nodes[i], idxS);
         SetCenterGP();
         temp =  FluidProp->Density()
            *MediaProp->Porosity(Index,pcs->m_num->ls_theta)
            *NodalVal_Sat[i] ;
         // Enthalpy
         dT=0.0;
         NodalVal2[i] = SolidProp->Enthalpy(NodalVal0[i], temp);
         if(fabs(NodalVal1[i]-NodalVal0[i])<1.0e-8)
            dT=1.0e-8;
         NodalVal3[i] = SolidProp->Enthalpy(NodalVal1[i]+dT, temp);
      }
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for mass matrix
   Programing:
   01/2005 WW/OK Implementation
   03/2005 WW Heat transport
   07/2005 WW Change for geometry element object
   08/2005 OK Gas flow
   10/2005 YD/OK: general concept for heat capacity
   11/2005 CMCD Heat capacity function included in mmp
   01/2007 OK Two-phase flow
   **************************************************************************/
   inline double CFiniteElementStd::CalCoefMass()
   {
      int Index = MeshElement->GetIndex();
      double val = 0.0;
      double humi = 1.0;
      double rhov = 0.0;
      double biot_val, poro_val, rho_val, Se;
      CompProperties *m_cp = NULL;

      if(pcs->m_num->ele_mass_lumping)
         ComputeShapefct(1);
      switch(PcsType)
      {
         default:
            std::cout << "Fatal error in CalCoefMass: No valid PCS type" << std::endl;
            break;
         case L:                                  // Liquid flow
                                                  // Is this really needed?
            val = MediaProp->StorageFunction(Index,unit,pcs->m_num->ls_theta);
            // JTARON 2010, needed storage term and fluid compressibility...
            // We derive here the storage at constant strain, or the inverse of Biot's "M" coefficient
            // Assumptions are the most general possible::  Invarience under "pi" (Detournay & Cheng) loading.
            // Se = 1/M = poro/Kf + (alpha-poro)/Ks    ::    Cf = 1/Kf = 1/rho * drho/dp    ::    alpha = 1 - K/Ks
            // Second term (of Se) below vanishes for incompressible grains
            rho_val = FluidProp->Density();
            if(D_Flag>0  && rho_val>MKleinsteZahl)
            {
               biot_val = SolidProp->biot_const;
               poro_val = MediaProp->Porosity(Index,pcs->m_num->ls_theta);

               val += poro_val * (FluidProp->drho_dp / rho_val) \
                  + (biot_val - poro_val) * (1.0 - biot_val) / SolidProp->K;
               // Will handle the dual porosity version later...
            }
            val /= time_unit_factor;
            break;
         case U:                                  // Unconfined flow
            break;
         case G:                                  // MB now Groundwater flow
            if(MediaProp->unconfined_flow_group>0)//OK
               val = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
            else
               val = MediaProp->StorageFunction(Index,unit,pcs->m_num->ls_theta);
            break;
         case T:                                  // Two-phase flow
            // val = (1/rho*n*d_rho/d_p*S + Se*S )
            if(pcs->pcs_type_number==0)
            {
               // PCH cpl_pcs gives a funny process number.
               // It is just the opposite of the phase. So, I get the value the other way around.
               idxS = cpl_pcs->GetNodeValueIndex("SATURATION2");
               for(int i=0;i<nnodes;i++)
                  NodalVal_Sat[i] = cpl_pcs->GetNodeValue(nodes[i],idxS+1);
               Sw = 1.0 - interpolate(NodalVal_Sat);
                                                  // Is this really needed?
               val = MediaProp->StorageFunction(Index,unit,pcs->m_num->ls_theta) * MMax(0.,Sw);

               // JTARON 2010, generalized poroelastic storage. See single phase version in case "L".
               // Se = 1/M = poro/Kf + (alpha-poro)/Ks
               rho_val = FluidProp->Density();
               if(D_Flag>0 && rho_val>MKleinsteZahl)
               {
                  biot_val = SolidProp->biot_const;
                  poro_val = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
                  Se = poro_val * (FluidProp->drho_dp / rho_val) \
                     + (biot_val - poro_val) * (1.0 - biot_val) / SolidProp->K;
                  // The poroelastic portion
                  val += Se*MMax(0.,Sw);
               }

               // If Partial-Pressure-Based model
               if(pcs->PartialPS == 1)
               {
                  // Let's get dPcdSw and apply to the mat_fac
                  CMediumProperties *m_mmp = NULL;
                  m_mmp = mmp_vector[0];

                  // dSedPc always return positive numbers for default case
                  // However, the value should be negative analytically.
                  double dSwdPc = m_mmp->SaturationPressureDependency(Sw, rho_val, 1);
                  val -= poro_val*dSwdPc;
               }
            }
            if(pcs->pcs_type_number==1)
            {
               val = MediaProp->Porosity(Index,pcs->m_num->ls_theta) * MediaProp->geo_area;
               // PCH: geo_area is only used for 1 and 2 dimensions.
               // This is originated from Old RockFlow that handles 1, 2, and 3 dimensional
               // elements in seperate functions that show inconsistency in handling geo_area.
               // For example, geo_area is never used for 3D even in Old RockFlow.
            }
            break;
         case C:                                  // Componental flow
            //OK comp = m_pcs->pcs_type_number;
            //OK coefficient = MPCCalcStorativityNumber(ele,phase,comp,gp);
            break;
            //....................................................................
         case H:                                  // Heat transport
            TG = interpolate(NodalVal1);
            val = MediaProp->HeatCapacity(Index,pcs->m_num->ls_theta,this);
            val /=time_unit_factor;
            break;
            //....................................................................
         case M:                                  // Mass transport //SB4200
                                                  // Porosity
            val = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
            val *= PCSGetEleMeanNodeSecondary_2(Index, pcs->flow_pcs_type, "SATURATION1", 1);
            //   val *= PCSGetEleMeanNodeSecondary(Index, "RICHARDS_FLOW", "SATURATION1", 1);
            m_cp = cp_vec[pcs->pcs_component_number];
                                                  //Retardation Factor
            val *= m_cp->CalcElementRetardationFactorNew(Index, unit, pcs);
            break;
         case O:                                  // Liquid flow
            val = 1.0;
            break;
         case R:                                  // Richards
            PG = interpolate(NodalVal1);          //12.02.2007.  Important! WW
                                                  //WW
            Sw = MediaProp->SaturationCapillaryPressureFunction(-PG,0);
            //     Sw = interpolate(NodalVal_Sat);
            rhow = FluidProp->Density();
            dSdp = MediaProp->SaturationPressureDependency(Sw, rhow, pcs->m_num->ls_theta);
            poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
            // Storativity
            val = MediaProp->StorageFunction(Index,unit,pcs->m_num->ls_theta) *Sw;

            // Fluid compressibility
            if(rhow>0.0)
               val += poro  *Sw* FluidProp->drho_dp / rhow;
            // Capillarity
            val += poro * dSdp;
                                                  //WW
            if(MediaProp->heat_diffusion_model==273)
            {
               //  	     PG = fabs(interpolate(NodalVal1));
               TG = interpolate(NodalValC)+T_KILVIN_ZERO;
               //Rv = GAS_CONSTANT;
               humi = exp(PG/(GAS_CONSTANT_V*TG*rhow));
               rhov = humi*FluidProp->vaporDensity(TG);
               //
               val -= poro * rhov*dSdp/rhow;
               val += (1.0-Sw)*poro*rhov/(rhow*rhow*GAS_CONSTANT_V*TG);
            }
            break;
         case F:                                  // Fluid Momentum
            val = 1.0;
            break;
         case A:                                  // Air (gas) flow
            val = MediaProp->Porosity(Index,pcs->m_num->ls_theta)/interpolate(NodalVal1);
            break;
      }
      return val;
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for mass matrix
   Programing:
   02/2007 WW Multi-phase flow
   05/2008 WW Generalization
   **************************************************************************/
   inline double CFiniteElementStd::CalCoefMass2(int dof_index)
   {
      int Index = MeshElement->GetIndex();
      double val = 0.0;
      double expfactor = 0.0;
      double dens_arg[3];                         //08.05.2008 WW
      bool diffusion = false;                     //08.05.2008 WW
      if(MediaProp->heat_diffusion_model==273&&cpl_pcs)
         diffusion = true;
      dens_arg[1] = 293.15;
      //
      if(pcs->m_num->ele_mass_lumping)
         ComputeShapefct(1);
      switch(dof_index)
      {
         case 0:
            PG = interpolate(NodalVal1);          // Capillary pressure
            dens_arg[0] = PG;                     // Should be P_w in some cases
            if(diffusion)
            {
               TG = interpolate(NodalValC1)+T_KILVIN_ZERO;
               dens_arg[1] = TG;
            }
            Sw = MediaProp->SaturationCapillaryPressureFunction(PG,0);
            rhow = FluidProp->Density(dens_arg);
            dSdp = -MediaProp->SaturationPressureDependency(Sw, rhow, pcs->m_num->ls_theta);
            poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
            // Storativity   28.05.2008
            // val = MediaProp->StorageFunction(Index,unit,pcs->m_num->ls_theta) *Sw;
            // Fluid compressibility
            // val += poro  *Sw* FluidProp->drho_dp / rhow;
            val = poro * dSdp;
            // Coupled (T)
            if(diffusion)
            {
               // Water vapour pressure
               expfactor = COMP_MOL_MASS_WATER/(rhow*GAS_CONSTANT*TG);
               rho_gw = FluidProp->vaporDensity(TG)*exp(-PG*expfactor);
               //
               val -= poro * dSdp*rho_gw/rhow;
               //
               val -= (1.0-Sw)*poro*COMP_MOL_MASS_WATER*rho_gw
                  /(rhow*GAS_CONSTANT*TG*rhow);
            }
            break;
         case 1:                                  //01
            val = 0.0;
            break;
         case 2:                                  //
            // (1-S)n(d rhop_c/d p_c)
            PG2 = interpolate(NodalVal_p2);
            dens_arg[0] = PG2;                    //28.05.2008. WW
            val = 0.;                             //28.05.2008. WW
            if(diffusion)                         //28.05.2008. WW
            {
               val = (1.0-Sw)*COMP_MOL_MASS_WATER*rho_gw
                  /(rhow*GAS_CONSTANT*TG*rhow);
               p_gw = rho_gw*GAS_CONSTANT*TG/COMP_MOL_MASS_WATER;
               dens_arg[0] -= p_gw;
               dens_arg[1] = TG;
            }
            rho_ga = GasProp->Density(dens_arg);  //28.05.2008. WW
            val -= rho_ga*dSdp/rhow;
            val *= poro;
            break;
         case 3:                                  //
            // Water vapour pressure
            if(diffusion)                         //28.05.2008. WW
               val = (1.0-Sw)*poro*GasProp->molar_mass/(GAS_CONSTANT*TG*rhow);
            else
               val = 0.;
            break;
      }
      return val;
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for mass matrix
   Programing:
   03/2009 PCH Multi-phase flow
   **************************************************************************/
   inline double CFiniteElementStd::CalCoefMassPSGLOBAL(int dof_index)
   {
      int Index = MeshElement->GetIndex();
      double val = 0.0;
      double P,T;
      //OK411 double expfactor = 0.0;
      double dens_arg[3];                         //08.05.2008 WW
      bool diffusion = false;                     //08.05.2008 WW
      if(MediaProp->heat_diffusion_model==273&&cpl_pcs)
         diffusion = true;
      dens_arg[1] = 293.15;
      //
      if(pcs->m_num->ele_mass_lumping)
         ComputeShapefct(1);
      switch(dof_index)
      {
         case 0:

            // compressibility also for the wetting phase NB
            poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
            Sw = 1.0-interpolate(NodalVal_SatNW); //Sw = 1-Snw
            P  = interpolate(NodalVal1);          //Pw
            T = interpolate(NodalValC1);

            val = poro*(Sw)*FluidProp->drhodP(P,T)/FluidProp->Density();
            //		cout << FluidProp->fluid_name << " Pressure: " << P << " Temp: " << ": drhodP: " << FluidProp->drhodP(P,T) << " density: " << FluidProp->Density() << endl;
            break;
         case 1:                                  // Snw in the wetting equation
            poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
            val = -poro;
            break;
         case 2:                                  // Pw in the non-wetting equation
            Sw = 1.0-interpolate(NodalVal_SatNW); //Sw = 1 - Snw
                                                  // Pnw = Pw + Pc(Sw)
            P = interpolate(NodalVal1) + MediaProp->CapillaryPressureFunction(0,NULL,0.0,1,Sw);
            // 	P = interpolate(NodalVal1);  // Pw
            T = interpolate(NodalValC1);

            val = poro*(1.-Sw)*GasProp->drhodP(P,T)/GasProp->Density();

            break;
         case 3:                                  // Snw in the non-wetting equation
            poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
            val = poro;
            break;
      }

      return val;
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for mass matrix
   Programing:
   01/2005 WW/OK Implementation
   03/2005 WW Heat transport
   07/2005 WW Change for geometry element object
   last modification:
   **************************************************************************/
   inline double CFiniteElementStd::CalCoefStorage()
   {
      int Index = MeshElement->GetIndex();
      double val = 0.0;
      CompProperties *m_cp =NULL;                 //CMCD
      //CompProperties *m_cp = cp_vec[pcs->pcs_component_number]; //SB4200
      switch(PcsType)
      {
         default:
            std::cout << "Fatal error in CalCoefStorage: No valid PCS type" << std::endl;
            break;
         case L:                                  // Liquid flow
            break;
         case U:                                  // Unconfined flow
            break;
         case G:                                  // MB now Groundwater flow
            break;
         case T:                                  // Two-phase flow
            break;
         case C:                                  // Componental flow
            break;
         case H:                                  // heat transport
            val = 0.0;
            break;
         case M:                                  // Mass transport //SB4200
                                                  //CMCD
            m_cp = cp_vec[pcs->pcs_component_number];
                                                  //Porosity
            val = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
            val *= PCSGetEleMeanNodeSecondary(Index, "RICHARDS_FLOW", "SATURATION1", 1);
                                                  // Decay rate
            val *= m_cp->CalcElementDecayRateNew(Index, pcs);
                                                  //Retardation Factor
            val *= m_cp->CalcElementRetardationFactorNew(Index, unit, pcs);
            break;
         case O:                                  // Liquid flow
            break;
         case R:                                  // Richards
            break;
         case F:                                  // Fluid Momentum
            break;
         case A:                                  // Air (gas) flow
            break;
      }
      return val;
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for Content matrix
   Programing:
   01/2005 WW/OK Implementation
   03/2005 WW Heat transport
   07/2005 WW Change for geometry element object
   last modification:
   **************************************************************************/
   inline double CFiniteElementStd::CalCoefContent()
   {
      int Index = MeshElement->GetIndex();
      double val = 0.0;
      double dS = 0.0;
      double nodeval0, nodeval1;
      //CompProperties *m_cp = NULL; //SB4200
      string name;

      switch(PcsType)
      {
         default:
            std::cout << "Fatal error in CalCoefContent: No valid PCS type" << std::endl;
            break;
         case L:                                  // Liquid flow
            break;
         case U:                                  // Unconfined flow
            break;
         case G:                                  // MB now Groundwater flow
            break;
         case T:                                  // Two-phase flow
            break;
         case C:                                  // Componental flow
            break;
         case H:                                  // heat transport
            break;
         case M:                                  // Mass transport //SB4200
         {
                                                  // Porosity
            val = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
            // Get saturation change:
            nodeval0 = PCSGetEleMeanNodeSecondary_2(Index, pcs->flow_pcs_type, "SATURATION1", 0);
            nodeval1 = PCSGetEleMeanNodeSecondary_2(Index, pcs->flow_pcs_type, "SATURATION1", 1);
            // 	nodeval0 = PCSGetEleMeanNodeSecondary(Index, "RICHARDS_FLOW", "SATURATION1", 0);
            //  nodeval1 = PCSGetEleMeanNodeSecondary(Index, "RICHARDS_FLOW", "SATURATION1", 1);
            dS = nodeval1 - nodeval0;             // 1/dt accounted for in assemble function
            //		if(Index == 195) cout << val << "Sat_old = " << nodeval0 << ", Sa_new: "<< nodeval1<< ", dS: " << dS << endl;
            val*= dS;
            break;
         }
         case O:                                  // Liquid flow
            break;
         case R:                                  // Richards
            break;
         case F:                                  // Fluid Momentum
            break;
         case A:                                  // Air (gas) flow
            break;
      }
      return val;
   }
   /**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for Laplacian matrix
   Programing:
   01/2005 WW/OK Implementation
   02/2005 OK Richards flow
   03/2005 WW Heat transport
   06/2005 OK Overland flow based on CalcEle2DQuad_OF by MB
   07/2005 WW Change for geometry element object
   08/2005 OK Air (gas) flow
   01/2007 OK Two-phase flow
   10/2008 PCH Two-phase flow modified
   **************************************************************************/
   inline void CFiniteElementStd::CalCoefLaplace(bool Gravity, int ip)
   {
      int i=0;
      double dens_arg[3];                         //AKS
      double mat_fac = 1.0;
      double Dpv = 0.0;
      double poro = 0.0;
      double tort = 0.0;
      double humi = 1.0;
      double rhow = 0.0;
      double *tensor = NULL;
      double Hav,manning,chezy,expp,chezy4,Ss,arg;
      static double Hn[9],z[9];
      double GradH[3],Gradz[3],w[3],v1[3],v2[3];
      int nidx1;
      int Index = MeshElement->GetIndex();
      double k_rel;
      ComputeShapefct(1);                         //  12.3.2007 WW
      double variables[3];                        //OK4709

      // For nodal value interpolation
      //======================================================================
      switch(PcsType)
      {
         default:
            break;
         case L:                                  // Liquid flow
            tensor = MediaProp->PermeabilityTensor(Index);
            if (ele_dim != dim)
            {
               Matrix local_tensor(dim,dim);
               Matrix temp_tensor(dim,dim);
               Matrix t_transform_tensor(*MeshElement->tranform_tensor);
               MeshElement->tranform_tensor->GetTranspose(t_transform_tensor);
               Matrix global_tensor(dim,dim);
               for (i=0; i<ele_dim; i++)
                  for (int j=0; j<ele_dim; j++)
                     local_tensor(i,j) = tensor[j+i*ele_dim];
               //cout << "K':" << endl; local_tensor.Write();
               local_tensor.multi(t_transform_tensor, temp_tensor);
               for (i=0; i<dim; i++)
               {
                  for (int j=0; j<dim; j++)
                     for (int k=0; k<dim; k++)
                        global_tensor(i,j)+=(*MeshElement->tranform_tensor)(i,k)*temp_tensor(k,j);
               }
               //cout << "K:" << endl; global_tensor.Write();
               for(i=0; i<dim; i++)
               {
                  for(int j=0; j<dim; j++)
                  {
                     tensor[dim*i+j] = global_tensor(i,j);
                  }
               }
            }
            variables[0] = interpolate(NodalVal1);//OK4709 pressure
            variables[1] = interpolate(NodalValC);//OK4709 temperature
                                                  //OK4709
            mat_fac = FluidProp->Viscosity(variables);
            //OK4709 mat_fac = FluidProp->Viscosity();
            if(gravity_constant<MKleinsteZahl)    // HEAD version
               mat_fac = 1.0;
            if(HEAD_Flag) mat_fac=1.0;
                                                  // Modified LBNL model WW
            if(MediaProp->permeability_stress_mode>1)
            {
               if(cpl_pcs)
                  TG = interpolate(NodalValC1)+T_KILVIN_ZERO;
               else
                  TG = 296.0;
               MediaProp->CalStressPermeabilityFactor(w, TG);
               for(i=0; i<dim; i++)
                  tensor[i*dim+i] *= w[i];
            }
            for(i=0; i<dim*dim; i++)
               mat[i] = tensor[i]/mat_fac;

            break;
         case G:                                  // Groundwater flow
            /* SB4218 - moved to ->PermeabilityTensor(Index);
                    if(MediaProp->permeability_model==2){ //?efficiency
                      for(i=0;i<(int)pcs->m_msh->mat_names_vector.size();i++){
                        if(pcs->m_msh->mat_names_vector[i].compare("PERMEABILITY")==0)
                          break;
                      }

                      mat_fac = MeshElement->mat_vector(i);
                      mat_fac /= FluidProp->Viscosity();
                     for(i=0; i<dim; i++) //WW
                        mat[i*dim+i] = mat_fac;
            }
            else{
            */
            tensor = MediaProp->PermeabilityTensor(Index);
            for(i=0;i<dim*dim;i++)
                                                  //16.10.2009 .WW
                  mat[i] = tensor[i]*time_unit_factor;
            break;
            //..................................................................
         case T:                                  // Two-phase flow
            // PCH Rewriting...
            // PCH Laplace mat_fac is accounted for two phases here.
            // thought to be related to the reference pressure.
            tensor = MediaProp->PermeabilityTensor(Index);
            if(pcs->pcs_type_number==0)
            {
               if(!Gravity)
               {
                  // PCH Laplace mat_fac is accounted for two phases here.
                  // thought to be related to the reference pressure.
                  int numOfPhases = 2;
                  double mat_fac = 0.0;
                  for(int p=0; p<numOfPhases; ++p)
                  {
                     // PCH Check if Capillary term is on
                     if(pcs->ML_Cap == 1)
                        p=1;

                     idxS = cpl_pcs->GetNodeValueIndex("SATURATION2");

                     for(i=0;i<nnodes;i++)
                        NodalVal_Sat[i] = cpl_pcs->GetNodeValue(nodes[i],idxS+1);

                     // Whatever phase, we use Sw for getting relative permeabilities
                     Sw = 1.0-interpolate(NodalVal_Sat);
                     k_rel = MediaProp->PermeabilitySaturationFunction(Sw,p);

                     // Note here mat_fac is += meaning adding two phases
                     mat_fac += time_unit_factor * k_rel / mfp_vector[p]->Viscosity();
                     // If Partial-Pressure-Based model
                     if(pcs->PartialPS == 1)
                        p=1;
                  }
                  for(i=0;i<dim*dim;i++)
                     mat[i] = tensor[i] * mat_fac;
               }
               else                               // Here is only active for Cal_Velocity
               {
                  // This is to calculate velocity
                  //WW					int numOfPhases = 2;
                  double mat_fac = 0.0;

                  idxS = cpl_pcs->GetNodeValueIndex("SATURATION2");

                  for(i=0;i<nnodes;i++)
                     NodalVal_Sat[i] = cpl_pcs->GetNodeValue(nodes[i],idxS+1);

                  Sw = 1.0-interpolate(NodalVal_Sat);
                  k_rel = MediaProp->PermeabilitySaturationFunction(Sw,pcs->pcs_type_number);

                  // Note here mat_fac is += meaning adding two phases
                  mat_fac = time_unit_factor * k_rel / mfp_vector[pcs->pcs_type_number]->Viscosity();

                  for(i=0;i<dim*dim;i++)
                     mat[i] = tensor[i] * mat_fac;
               }
            }
            else if(pcs->pcs_type_number==1)
            {
               // PCH Check if Capillary term is on
               if(pcs->ML_Cap == 1)
               {
                  int phase = pcs->pcs_type_number;

                  idxS = pcs->GetNodeValueIndex("SATURATION2");
                  for(i=0;i<nnodes;i++)
                     NodalVal_Sat[i] = pcs->GetNodeValue(nodes[i],idxS+1);
                  Sw = 1.0-interpolate(NodalVal_Sat);
                  k_rel = MediaProp->PermeabilitySaturationFunction(Sw,phase);

                  // Here only the second phase accounted.
                  // Let's get dPcdSw and apply to the mat_fac
                  CMediumProperties *m_mmp = NULL;
                  m_mmp = mmp_vector[0];
                  double dPcdSw=0.0;
                  if (m_mmp->capillary_pressure_model_values[0] < MKleinsteZahl)
                     dPcdSw = 0.0;
                  else
                     dPcdSw = 1.0/m_mmp->SaturationPressureDependency(Sw, 1000.0, 1);
                  mat_fac = time_unit_factor * k_rel / mfp_vector[phase]->Viscosity()*(-dPcdSw);
                  for(i=0;i<dim*dim;i++)
                     mat[i] = tensor[i] * mat_fac;
               }
               else
               {
                  int phase = pcs->pcs_type_number;

                  idxS = pcs->GetNodeValueIndex("SATURATION2");
                  for(i=0;i<nnodes;i++)
                     NodalVal_Sat[i] = pcs->GetNodeValue(nodes[i],idxS+1);
                  Sw = 1.0-interpolate(NodalVal_Sat);
                  k_rel = MediaProp->PermeabilitySaturationFunction(Sw,phase);

                  // Here only the second phase accounted.
                  mat_fac = time_unit_factor * k_rel / mfp_vector[phase]->Viscosity();
                  for(i=0;i<dim*dim;i++)
                     mat[i] = tensor[i] * mat_fac;
               }
            }

            break;
            //..................................................................
         case C:                                  // Componental flow
            break;
         case H:                                  // heat transport
            if(SolidProp->GetConductModel()==2)   // Boiling model. DECOVALEX THM2
            {
               TG = interpolate(NodalVal1);
               for(i=0; i<dim*dim; i++) mat[i] = 0.0;
               for(i=0; i<dim; i++)
                  mat[i*dim+i] = SolidProp->Heat_Conductivity(TG);
            }
                                                  // DECOVALEX THM1 or Curce 12.09. WW
            else if(SolidProp->GetConductModel()==3||SolidProp->GetConductModel()==4)
            {
               // WW
               PG = interpolate(NodalValC1);
               if(cpl_pcs->type!=1212)
                  PG *= -1.0;
               Sw = MediaProp->SaturationCapillaryPressureFunction(PG,0);
               for(i=0; i<dim*dim; i++) mat[i] = 0.0;
               mat_fac = SolidProp->Heat_Conductivity(Sw);
               for(i=0; i<dim; i++)
                  mat[i*dim+i] = mat_fac;
            }
            //WW        else if(SolidProp->GetCapacityModel()==1 && MediaProp->heat_diffusion_model == 273){
            else if(SolidProp->GetConductModel()==1)
            {
               TG = interpolate(NodalVal1);
               tensor = MediaProp->HeatDispersionTensorNew(ip);
               for(i=0;i<dim*dim;i++)
                  mat[i] = tensor[i];
            }
            else
            {
               tensor = MediaProp->HeatConductivityTensor(Index);
               for(i=0; i<dim*dim; i++)
                  mat[i] = tensor[i];             //mat[i*dim+i] = tensor[i];
            }
            break;
         case M:                                  // Mass transport
            mat_fac = 1.0;                        //MediaProp->Porosity(Index,pcs->m_num->ls_theta); // porosity now included in MassDispersionTensorNew()
            tensor = MediaProp->MassDispersionTensorNew(ip);
            //CB
            //SB->CB I think this does not belong here
            // mat_fac *= PCSGetEleMeanNodeSecondary_2(Index, pcs->flow_pcs_type, "SATURATION1", 1);
            //if(PCSGet("RICHARDS_FLOW"))
            //	    mat_fac *= PCSGetEleMeanNodeSecondary(Index, "RICHARDS_FLOW", "SATURATION1", 1);
            for(i=0;i<dim*dim;i++)
            {
               mat[i] = tensor[i]*mat_fac*time_unit_factor;
            }
            break;
            //------------------------------------------------------------------
         case O:                                  // Overland flow
            //................................................................
            // H - water level
            nidx1 = pcs->GetNodeValueIndex("HEAD")+1;
            Hav = 0.0;
            for(i=0;i<nnodes;i++)
            {
               z[i] = MeshElement->nodes[i]->Z();
               Hn[i] = pcs->GetNodeValue(MeshElement->nodes_index[i],nidx1) - z[i];
               if (Hn[i] < 0.0) {Hn[i] = 0.0;}
               Hav += Hn[i]/(double)nnodes;
            }
            //................................................................
            // Friction coefficient
            tensor = MediaProp->PermeabilityTensor(Index);
                                                  // Manning-coefficient: n
            manning = MediaProp->permeability_tensor[0];
            // ToDo MB MMP function: m_mmp->FrictionCoefficientChezy(gp)
            if (MediaProp->conductivity_model==3) // Chezy-coefficient C
            {
               expp = 1.0/6.0;
               chezy = pow(Hav,expp) / manning;   // f? b >> h gilt: C = H**1/6 n**-1
               // Grad H: grad_N H J^-1
               MMultMatVec(dshapefct,dim,nnodes,Hn,nnodes,v1,dim);
               MMultVecMat(v1,dim,invJacobian,dim,dim,GradH,dim);
               // Grad z: ? s.Z.380ff
               MMultMatVec(dshapefct,dim,nnodes,z,nnodes,v2,dim);
               MMultVecMat(v2,dim,invJacobian,dim,dim,Gradz,dim);
               w[0] = GradH[0] + Gradz[0];
               w[1] = GradH[1] + Gradz[1];
               chezy4 = pow(chezy,4);
               Ss = ((w[0] * w[0]) / chezy4) +  ((w[1] * w[1]) / chezy4);
               Ss = pow(Ss,0.25);
               if (fabs(Ss) < 1.0e-7)
               {
                  Ss = 1.0e-7;
               }
               expp =  5.0 / 3.0;
               arg = (pow(Hav,expp))/(chezy*chezy);
               mat_fac = arg / Ss;
            }
            //................................................................
            // Tensor
            for(i=0;i<dim*dim;i++)
                                                  //ToDo
                  mat[i] = tensor[i]/manning * mat_fac;
            break;
            //------------------------------------------------------------------
         case R:                                  // Richards flow
            // The following line only applies when Fluid Momentum is on
            PG = interpolate(NodalVal1);          //05.01.07 WW
                                                  //05.01.07 WW
            Sw = MediaProp->SaturationCapillaryPressureFunction(-PG,0);

            tensor = MediaProp->PermeabilityTensor(Index);
            mat_fac = time_unit_factor* MediaProp->PermeabilitySaturationFunction(Sw,0) \
               / FluidProp->Viscosity();
                                                  // Modified LBNL model WW
            if(MediaProp->permeability_stress_mode>1)
            {
               if(cpl_pcs)
                  TG = interpolate(NodalValC1)+T_KILVIN_ZERO;
               else
                  TG = 296.0;
               MediaProp->CalStressPermeabilityFactor(w, TG);
               for(i=0; i<dim; i++)
                  tensor[i*dim+i] *= w[i];
            }
            //
            for(i=0; i<dim*dim; i++)
               mat[i] = tensor[i] * mat_fac;
            if(MediaProp->heat_diffusion_model==273&&!Gravity)
            {
               rhow = FluidProp->Density();
               //PG = fabs(interpolate(NodalVal1));
               TG = interpolate(NodalValC)+T_KILVIN_ZERO;
               poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
               tort = MediaProp->TortuosityFunction(Index,unit,pcs->m_num->ls_theta);
               //Rv = GAS_CONSTANT;
               humi = exp(PG/(GAS_CONSTANT_V*TG*rhow));
               //
               Dpv = 2.16e-5*tort*(1-Sw)*poro*pow(TG/T_KILVIN_ZERO, 1.8);
               Dpv *= time_unit_factor*FluidProp->vaporDensity(TG)*humi/(GAS_CONSTANT_V*rhow*TG);
               for(i=0; i<dim; i++)
                  mat[i*dim+i] += Dpv/rhow;
            }
            break;
            //------------------------------------------------------------------
         case A:                                  // Air flow
            dens_arg[0] = interpolate(NodalVal1);
            dens_arg[1] = interpolate(NodalValC1)+T_KILVIN_ZERO;
            dens_arg[2] = Index;
            mat_fac = FluidProp->Viscosity(dens_arg);
            tensor = MediaProp->PermeabilityTensor(Index);
            for(i=0;i<dim*dim;i++)
               mat[i] = tensor[i]/mat_fac;
            break;
            //------------------------------------------------------------------
      }
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for Laplacian matrix
   Programing:
   10/2008 PCH Implementation
   **************************************************************************/
   inline void CFiniteElementStd::CalCoefLaplaceMultiphase(int phase, int ip)
   {
      ip = ip;                                    //OK411

      int i=0;
      double mat_fac = 1.0;
      double *tensor = NULL;
      // static double Hn[9],z[9];
      int Index = MeshElement->GetIndex();
      double k_rel;
      ComputeShapefct(1);                         //  12.3.2007 WW

      // For nodal value interpolation
      //======================================================================
      switch(PcsType)
      {
         default:
            break;

         case T:                                  // Two-phase flow
            // PCH Rewriting...
            // PCH Laplace mat_fac is accounted for two phases here.
            // thought to be related to the reference pressure.
            tensor = MediaProp->PermeabilityTensor(Index);
            if(pcs->pcs_type_number==0)
            {
               // PCH Laplace mat_fac is accounted for two phases here.
               // thought to be related to the reference pressure.
               double mat_fac = 0.0;

               idxS = cpl_pcs->GetNodeValueIndex("SATURATION2");

               for(i=0;i<nnodes;i++)
                  NodalVal_Sat[i] = cpl_pcs->GetNodeValue(nodes[i],idxS+1);
               Sw = 1.0-interpolate(NodalVal_Sat);

               k_rel = MediaProp->PermeabilitySaturationFunction(Sw,phase);

               // Note here mat_fac is += meaning adding two phases
               mat_fac = time_unit_factor * k_rel / mfp_vector[phase]->Viscosity();

               for(i=0;i<dim*dim;i++)
                  mat[i] = tensor[i] * mat_fac;
            }
            else if(pcs->pcs_type_number==1)
            {
               int phase = pcs->pcs_type_number;

               idxS = pcs->GetNodeValueIndex("SATURATION2");
               for(i=0;i<nnodes;i++)
                  NodalVal_Sat[i] = pcs->GetNodeValue(nodes[i],idxS+1);
               Sw = 1.0 - interpolate(NodalVal_Sat);
               k_rel = MediaProp->PermeabilitySaturationFunction(Sw,phase);

               // Here only the second phase accounted.
               mat_fac = time_unit_factor * k_rel / mfp_vector[phase]->Viscosity();
               for(i=0;i<dim*dim;i++)
                  mat[i] = tensor[i] * mat_fac;
            }
            break;
      }
   }

   ///////
   /**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for Laplacian matrix of multi-phase
         flow
   Programing:
   02/2007 WW Implementation
   last modification:
   **************************************************************************/
   inline void CFiniteElementStd::CalCoefLaplace2(bool Gravity,  int dof_index)
   {
      int i=0;
      double *tensor = NULL;
      double mat_fac = 1.0, m_fac=0.;
      double expfactor, D_gw, D_ga;
      expfactor = D_gw = D_ga =0.0;
      double dens_arg[3];                         //08.05.2008 WW
      bool diffusion = false;                     //08.05.2008 WW
      if(MediaProp->heat_diffusion_model==273&&cpl_pcs)
         diffusion = true;
      //
      dens_arg[1] = 293.15;
      //
      int Index = MeshElement->GetIndex();
      //
      ComputeShapefct(1);                         //  12.3.2007 WW
      //======================================================================
      for(i=0; i<dim*dim; i++)
         mat[i] = 0.0;
      switch(dof_index)
      {
         case 0:
            PG = interpolate(NodalVal1);
            Sw = MediaProp->SaturationCapillaryPressureFunction(PG,0);
            //
            tensor = MediaProp->PermeabilityTensor(Index);
            mat_fac = MediaProp->PermeabilitySaturationFunction(Sw,0) \
               / FluidProp->Viscosity();
            for(i=0; i<dim*dim; i++)
               mat[i] = -tensor[i] * mat_fac*time_unit_factor;
            // For velocity caculation
            if(!Gravity)
            {
               dens_arg[0] = PG;                  // Shdould be Pw in some cases
               if(diffusion)
               {
                  TG = interpolate(NodalValC1)+T_KILVIN_ZERO;
                  dens_arg[1] = TG;
               }
               //
               rhow = FluidProp->Density(dens_arg);
               poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
               PG2 = interpolate(NodalVal_p2);
               dens_arg[0] = PG2;
               //
               if(diffusion)
               {
                  tort = MediaProp->TortuosityFunction(Index,unit,pcs->m_num->ls_theta);
                  tort *=(1.0-Sw)*poro*2.16e-5*pow(TG/T_KILVIN_ZERO, 1.8);
                  expfactor = COMP_MOL_MASS_WATER/(rhow*GAS_CONSTANT*TG);
                  rho_gw = FluidProp->vaporDensity(TG)*exp(-PG*expfactor);
                  p_gw = rho_gw*GAS_CONSTANT*TG/COMP_MOL_MASS_WATER;
                  dens_arg[0] -= p_gw;
               }
               //
               rho_ga = GasProp->Density(dens_arg);
               //
               if(diffusion)
               {
                  rho_g = rho_ga+rho_gw;
                  // 1/Mg
                  M_g = (rho_gw/COMP_MOL_MASS_WATER+rho_ga/GasProp->molar_mass)/rho_g;
                  D_gw = tort*rho_g*COMP_MOL_MASS_WATER*GasProp->molar_mass*M_g*M_g/rhow;
                  D_gw *= rho_gw/(rhow*PG2);
                  for(i=0; i<dim; i++)
                     mat[i*dim+i] -= D_gw*time_unit_factor;
               }
            }
            break;
         case 1:
            if(Gravity)
            {
               PG = interpolate(NodalVal1);
               PG2 = interpolate(NodalVal_p2);
               Sw = MediaProp->SaturationCapillaryPressureFunction(PG,0);
               dens_arg[0] = PG;                  // Shdould be Pw in some cases
               if(diffusion)
               {
                  TG = interpolate(NodalValC1)+T_KILVIN_ZERO;
                  dens_arg[1] = TG;
               }
               // Liquid density
               rhow = FluidProp->Density(dens_arg);
               dens_arg[0] = PG2;
               rho_gw = 0.0;
               if(diffusion)
               {
                  expfactor = COMP_MOL_MASS_WATER/(rhow*GAS_CONSTANT*TG);
                  rho_gw = FluidProp->vaporDensity(TG)*exp(-PG*expfactor);
                  p_gw = rho_gw*GAS_CONSTANT*TG/COMP_MOL_MASS_WATER;
                  dens_arg[0] -= p_gw;
               }
               rho_ga = GasProp->Density(dens_arg);
               rho_g = rho_ga+rho_gw;
            }
            tensor = MediaProp->PermeabilityTensor(Index);
            mat_fac = MediaProp->PermeabilitySaturationFunction(Sw,0) \
               / FluidProp->Viscosity();
            m_fac = 0.;
            if(diffusion)
               m_fac = rho_gw*MediaProp->PermeabilitySaturationFunction(Sw,1) \
                  / (GasProp->Viscosity()*rhow);
            if(Gravity)
               mat_fac = mat_fac+m_fac*rho_g/rhow;
            else
               mat_fac += m_fac;
            //
            for(i=0; i<dim*dim; i++)
               mat[i] = tensor[i] * mat_fac*time_unit_factor;
            //
            if((!Gravity)&&diffusion)
            {
               D_gw = tort*COMP_MOL_MASS_WATER*GasProp->molar_mass*M_g*M_g*rho_g/rhow;
               D_gw *= time_unit_factor*p_gw/(PG2*PG2);
               for(i=0; i<dim; i++)
                  mat[i*dim+i] -= D_gw;
            }
            break;
         case 2:
            if(diffusion)
            {
               D_ga = tort*COMP_MOL_MASS_WATER*GasProp->molar_mass*M_g*M_g*rho_g/rhow;
               D_ga *= time_unit_factor*rho_gw/(PG2*rhow);
            }
            else
               D_ga = 0.;
            for(i=0; i<dim; i++)
               mat[i*dim+i] = D_ga;
            break;
         case 3:
            //
            tensor = MediaProp->PermeabilityTensor(Index);
            mat_fac = rho_ga*MediaProp->PermeabilitySaturationFunction(Sw,1) \
               / (GasProp->Viscosity()*rhow);
            //
            if(Gravity)
               //        mat_fac *= rhow/rho_ga;
               mat_fac *= rho_ga/rhow;            //29.04.2009 WW
            //
            for(i=0; i<dim*dim; i++)
               mat[i] = tensor[i] * mat_fac*time_unit_factor;
            if((!Gravity)&&diffusion)
            {
               D_ga = tort*rho_g*COMP_MOL_MASS_WATER*GasProp->molar_mass*M_g*M_g/rhow;
               D_ga *= p_gw/(PG2*PG2);
               for(i=0; i<dim; i++)
                  mat[i*dim+i] += D_ga*time_unit_factor;
            }
            break;
            //------------------------------------------------------------------
      }
   }
   /**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for Laplacian matrix of PS multi-phase
         flow
   Programing:
   03/2009 PCH Implementation
   last modification:
   **************************************************************************/
   inline void CFiniteElementStd::CalCoefLaplacePSGLOBAL(bool Gravity,  int dof_index)
   {
      int i=0;
      double *tensor = NULL;
      double mat_fac = 1.0;                       //OK411 m_fac=0.;
      double k_rel=0.0;
      double mfp_arg[2];

      int Index = MeshElement->GetIndex();
      //
      ComputeShapefct(1);                         //  12.3.2007 WW
      //======================================================================
      for(i=0; i<dim*dim; i++)
         mat[i] = 0.0;
      switch(dof_index)
      {
         case 0:
            tensor = MediaProp->PermeabilityTensor(Index);
            if(pcs->m_num->ele_upwinding == 1)
            {
               // Doing Upwind elements for saturation by divergent of pressure.
               // Pw upwind
               int WhichNode = UpwindElement((int)(pcs->m_num->ele_upwind_method), 0);
               Sw = 1.0 - NodalVal_SatNW[WhichNode];
            }
            else
            {
               Sw = 1.0-interpolate(NodalVal_SatNW);
            }
            k_rel = MediaProp->PermeabilitySaturationFunction(Sw,0);

            mat_fac = k_rel / FluidProp->Viscosity();

            // Since gravity for water phase is handled directly in Assemble_Gravity,
            // no need of any code for water phase here.
            for(i=0;i<dim*dim;i++)
               mat[i] = tensor[i] * mat_fac*time_unit_factor;
            break;
         case 1:
            tensor = MediaProp->PermeabilityTensor(Index);
            mat_fac = 0.0;                        // Snw has no laplace term
            for(i=0; i<dim*dim; i++)
               mat[i] = tensor[i] * mat_fac*time_unit_factor;
            break;
         case 2:
            tensor = MediaProp->PermeabilityTensor(Index);
            if(pcs->m_num->ele_upwinding == 1)
            {
               // Doing Upwind elements for saturation by divergent of pressure.
               // Pnw upwind
               int WhichNode = UpwindElement((int)(pcs->m_num->ele_upwind_method), 0);
               Sw = 1.0 - NodalVal_SatNW[WhichNode];
            }
            else
            {
               Sw = 1.0-interpolate(NodalVal_SatNW);
            }

            k_rel = MediaProp->PermeabilitySaturationFunction(Sw,1);
                                                  // Pnw = Pw + Pc(Sw) //TODO: could cause errors in some cases
            mfp_arg[0] = interpolate(NodalVal1) + MediaProp->CapillaryPressureFunction(0,NULL,0.0,1,Sw);
            mfp_arg[1] = interpolate(NodalValC1); // TEMPERATURE1 in most cases

            mat_fac = k_rel / GasProp->Viscosity(mfp_arg);

            // The density of the non-wetting phase fluid should be considered here.
            // However, the default water phase density should be canceled out simultaneously.
            if(Gravity)
               mat_fac *= GasProp->Density()/FluidProp->Density();

            for(i=0;i<dim*dim;i++)
               mat[i] = tensor[i] * mat_fac*time_unit_factor;
            break;
         case 3:
            // Snw Laplace from Pc in Eqn 2
            tensor = MediaProp->PermeabilityTensor(Index);

            if(pcs->num_type_name.find("dPcdSwGradSnw")!=string::npos)
            {
               double Snw = -1.0;
               if(pcs->m_num->ele_upwinding == 1)
               {
                  // Doing Upwind elements for saturation by divergent of pressure.
                  // Pnw upwind
                  int WhichNode = UpwindElement((int)(pcs->m_num->ele_upwind_method), 1);
                  Snw = NodalVal_SatNW[WhichNode];
               }
               else
               {
                  Snw = interpolate(NodalVal_SatNW);
               }

               CMediumProperties *m_mmp = NULL;
               CElem* thisEle = pcs->m_msh->ele_vector[index];
               int matgrp = thisEle->GetPatchIndex();
               m_mmp = mmp_vector[matgrp];
               //		Sw = 1.0 - Snw;
               Sw = 1.0 - MRange(m_mmp->saturation_res[1],Snw,1.0);
               //		Sw = 1.0 - MRange(m_mmp->saturation_res[1],Snw,1.0-m_mmp->saturation_res[0]);
               k_rel = MediaProp->PermeabilitySaturationFunction(Sw,1);

               double dPcdSw=0.0;
               if (m_mmp->capillary_pressure_model_values[0] < MKleinsteZahl)
                  dPcdSw = 0.0;
               else
                  dPcdSw=m_mmp->PressureSaturationDependency(Sw, 1000.0, 1);

                                                  // Pnw = Pw + Pc(Sw) // TODO: could cause errors in some cases
               mfp_arg[0] = interpolate(NodalVal1) + MediaProp->CapillaryPressureFunction(0,NULL,0.0,1,Sw);
               mfp_arg[1] = interpolate(NodalValC1);
               mat_fac = k_rel / GasProp->Viscosity(mfp_arg)*(-dPcdSw);
            }
            else
               mat_fac = 0.0;

            for(i=0; i<dim*dim; i++)
               mat[i] = tensor[i] * mat_fac*time_unit_factor;
            break;
         case 4:
            // For Vnw
            tensor = MediaProp->PermeabilityTensor(Index);

            double Snw = -1.0;
            if(pcs->m_num->ele_upwinding == 1)
            {
               // Doing Upwind elements for saturation by divergent of pressure.
               // Pnw upwind
               int WhichNode = UpwindElement((int)(pcs->m_num->ele_upwind_method), 1);
               Snw = NodalVal_SatNW[WhichNode];
            }
            else
            {
               Snw = interpolate(NodalVal_SatNW);
            }

            CMediumProperties *m_mmp = NULL;
            CElem* thisEle = pcs->m_msh->ele_vector[index];
            int matgrp = thisEle->GetPatchIndex();
            m_mmp = mmp_vector[matgrp];
            //		Sw = 1.0 - Snw;
            Sw = 1.0 - MRange(m_mmp->saturation_res[1],Snw,1.0);
            //	Sw = 1.0 - MRange(m_mmp->saturation_res[1],Snw,1.0-m_mmp->saturation_res[0]);
            k_rel = MediaProp->PermeabilitySaturationFunction(Sw,1);
            mat_fac = k_rel / GasProp->Viscosity();

            for(i=0; i<dim*dim; i++)
               mat[i] = tensor[i] * mat_fac*time_unit_factor;
            break;
            //------------------------------------------------------------------
      }
   }
   /**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2009 PCH Upwind Material Scheme
   Background:
     Now material property at the upstream node is taken to exclude future
     predition at the current due to abrupt change of material properties.
     This is conservative perticularly the nodes very close to the interface
     that divides two highly different materials. Thus,this is vector
     characterisitics. In our case, pressure gradient can be a determinant
   for permeability and saturation and so other material properties.

   Description:
   0: none yet
   1: Upwind element determined by div of pressure
   **************************************************************************/
   int CFiniteElementStd::UpwindElement(int option, int phase)
   {
      int WhichNodeInTheElement = -1;             // Initialized to be none of elements
      double Pmin = 1.0/DBL_MIN;                  // Just set to be outrageously big.
      int GravityOn = 1;                          // Initialized to be on

      // If no gravity, then set GravityOn to be zero.
      if((coordinate_system)%10!=2&&(!axisymmetry))
         GravityOn = 0;

      if(option==1)                               // If upwind by divergent of pressure
      {
         double PdivAtnode = -1.0/DBL_MIN;        // Meaningless pressure.
         double Pc = 0.0;                         // Set to be no capillary at all.
         int idx_p1 = pcs->GetNodeValueIndex("PRESSURE1");
         int idx_pc = pcs->GetNodeValueIndex("PRESSURE_CAP");

         for(int i=0; i<nnodes; ++i)
         {
            double Pw = pcs->GetNodeValue(nodes[i],idx_p1+1);
            if(phase==0)
            {
               PdivAtnode = Pw;
               if(GravityOn)
                  PdivAtnode -= FluidProp->Density()*gravity_constant;
            }
            else if(phase==1)
            {
               Pc = pcs->GetNodeValue(nodes[i],idx_pc);
               PdivAtnode = Pw+Pc;
               if(GravityOn)
                  PdivAtnode -= GasProp->Density()*gravity_constant;
            }
            else
            {
               std::cout << "Phase number is wrong in UpwindElement."<< std::endl;
               abort();
            }

            if(Pmin > PdivAtnode)
            {
               Pmin = PdivAtnode;
               WhichNodeInTheElement = i;
            }
         }
      }

      if(WhichNodeInTheElement == -1)
      {
         std::cout<<"UpwindElement is failed. Impossible node index!!!"<< std::endl;
         std::cout<<"Pmin = "<<Pmin<< std::endl;
         abort();
      }
      return WhichNodeInTheElement;
   }
   //CB 090507
   /**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2007 CB
   last modification:
   **************************************************************************/
   //inline void CFiniteElementStd::UpwindUnitCoord(int p, int point, int ind, double *rupw, double *supw, double *tupw)
   inline void CFiniteElementStd::UpwindUnitCoord(int p, int point, int ind)
   {
      p = p;                                      //OK411
      //Laufvariablen
      static long i, j;                           //, k, l;
      // static long *element_nodes;
      double scale;
      double alpha[3];
      int gp_r, gp_s, gp_t;
      bool mmp_element_integration_method_maximum = false;

      int upwind_meth;
      double upwind_para ;
      double v[3], v_rst[3];

      //
      ElementValue* gp_ele = ele_gp_value[ind];
      if(pcs->pcs_type_number==1)                 //WW/CB
         gp_ele = ele_gp_value[ind+(long)pcs->m_msh->ele_vector.size()];

      // TF unused:  MshElemType::type eletyp = MeshElement->GetElementType();

      //
      upwind_para = pcs->m_num->ele_upwinding;
      upwind_meth = pcs->m_num->ele_upwind_method;

      // Numerik
      // *rupw = *supw = *tupw = 0.0;
      gp_r = gp_s = gp_t = 0;
      // alpha initialisieren
      MNulleVec(alpha, 3);
      // v initialisieren
      MNulleVec(v,3);

      // get the velocities
      //gp_ele->GetEleVelocity(v);

      // CB: not sure if this is correct, as velocity
      // at each Gauss point is regarded here (within GP loop)
      // while in cel_mmp.cpp velocity is evaluated before the GP loop:
      // CalcVelo3Drst(phase, index, GetTimeCollocationupwind_MMP(), 0., 0., 0., v);
      //v[0] = gp_ele->Velocity(0, point);
      //v[1] = gp_ele->Velocity(1, point);
      //v[2] = gp_ele->Velocity(2, point);
      v[0] = gp_ele->Velocity(0, 0);
      v[1] = gp_ele->Velocity(1, 0);
      v[2] = gp_ele->Velocity(2, 0);
      // this would give v at GP 0
      // but: if(PcsType==T) SetCenterGP(); // CB 11/2007
      // was set in Cal_Velo_2(); which is calculated,
      // when mass matrix is assembled
      // hence V is at element center of gravity
      // otherwise use element averaged v?:
      //for(i=0; i<nGaussPoints; i++)
      //{
      //  v[0] += gp_ele->Velocity(0, i)/(double)nGaussPoints;
      //  v[1] += gp_ele->Velocity(1, i)/(double)nGaussPoints;
      //  v[2] += gp_ele->Velocity(2, i)/(double)nGaussPoints;
      //}

      for (i=0; i<3; i++)
         v[i] *= time_unit_factor;

      //instead of r_upw, etc. we use unit[i], next function sets Gauss Integrals
      //unit[0] = MXPGaussPkt(nGauss, gp_r); -> r_upw ; etc. for unit[i]
      SetGaussPoint(point, gp_r, gp_s, gp_t);     // this sets unit[] to standard coordinates

      // v transformation: Jacobi*v-->v_rst
      // computing the Jacobian at this point results in pressure difference
      // in comparison to the original model of CT
      // However, it seems not necessary, as Jacobian has already
      // been computed in function UpwindAlphaMass
      // computeJacobian(1); // order 1

      // multiply velocity vector with Jacobian matrix
      // Jacobi*v-->v_rst
      // This may need attention to tell different types of elements in the same dimension	// PCH
      for (i=0; i<ele_dim; i++)
      {
         v_rst[i] = 0.0;
         for (j=0; j<ele_dim; j++)
            v_rst[i] += Jacobian[i*dim+j]*v[j];
      }
      //

      // These need to be rewritten according to different types of elements.  // PCH
      if(MBtrgVec(v_rst, ele_dim) > MKleinsteZahl)
      {
         // Upwind-Faktoren
         for(i=0;i<ele_dim;i++)
            alpha[i] = -upwind_para * v_rst[i] / (MBtrgVec(v_rst, ele_dim) + MKleinsteZahl);

         // moving the Gauss points
         if (upwind_meth == 1)                    // limit Gauss point moving on Element domain
         {
            scale = 1.;
            for(i=0;i<ele_dim;i++)
            {
                                                  // Integral over GaussPoints, not used
               if (mmp_element_integration_method_maximum)
               {
                  if (fabs(unit[i] + alpha[i]) > 1.)
                     scale = MMin(scale, (1. - fabs(unit[i])) / fabs(alpha[i]));
               }
               else                               // regard all quantities in the center of element
               {
                  if (fabs(alpha[i]) > 1.)
                     scale = MMin(scale, (1./fabs(alpha[i])) );
               }
            }
            for(i=0;i<ele_dim;i++)
            {
                                                  // Integral over GaussPoints, not used
               if (mmp_element_integration_method_maximum)
                  unit[i] += scale * alpha[i];    // scale is added to unit[i] (=Gaussintegral)
               else                               // regard all quantities in the center of element
                  unit[i] = scale * alpha[i];     // unit[i] (=Gaussintegral)
            }
         }
         else if (upwind_meth == 2)               // limit moving on -1<x<1
         {
            // PCH this has never been used, but the code is only for line, quad, and hex.
            for(i=0;i<ele_dim;i++)
            {
                                                  // Integral ?er GaussPunkte
               if (mmp_element_integration_method_maximum)
                  unit[i] = MRange(-1., unit[i] + alpha[i], 1.);
               else                               // regard all quantities in the center of element
                  unit[i] = MRange(-1., alpha[i], 1.);
            }
         }
         //here only Methods 1 + 2; M3 = MaxMobilUW is done in CalcCoefLaplace
      }

#ifdef OLD_UPWINDING
      //test
      for(i=0;i<ele_dim;i++)
         cout << unit[i] << " ";
      cout << endl;

      double ur, us, ut;

      switch(eletyp)
      {
         case 1:                                  // Line
         {
            // Elementgeometriedaten
            static double detjac, *invjac, jacobi[4];
            double l[3];
            //invjac = GetElementJacobiMatrix(index, &detjac);
            //Calc1DElementJacobiMatrix(ind, invjac, &detjac);  //ind = element id number  // wird das irgendwo gebraucht?
            detjac = computeJacobian(1);          // order
            invjac = invJacobian;

            MNulleVec(l,3);
            l[0] = X[1] - X[0];
            l[1] = Y[1] - Y[0];
            l[2] = Z[1] - Z[0];

            //hier nur Methoden 1 + 2; 3 wird in CalcCoefLaplace erledigt
            if ((upwind_meth == 1) || (upwind_meth == 2))
            {
               if (MBtrgVec(v, 3) > MKleinsteZahl)
               {
                  if (MSkalarprodukt(v, l, 3) > 0.)
                                                  // CB VZ ge?dert!!!
                     *rupw = MRange(-1., -upwind_para , 1.);
                  //*rupw = MRange(-1., upwind_para , 1.);
                  else
                                                  // CB VZ ge?dert!!!
                     *rupw = MRange(-1., upwind_para , 1.);
                  //*rupw = MRange(-1., -upwind_para , 1.); //
               }
               //else { // test
               //    cout << "-vau";
               //    if (MSkalarprodukt(v, l, 3) > 0.)
               //     *rupw = MRange(-1., upwind_para , 1.);
               //    else
               //     *rupw = MRange(-1., -upwind_para , 1.);
               //}
               if(aktueller_zeitschritt==1)       // test
                  *rupw = MRange(-1., upwind_para , 1.);
            }
            // Upwind-Faktor Fully upwinding
         }
         break;
         case 2:                                  // Quadrilateral
         {
            // Elementgeometriedaten
            static double detjac, *invjac, jacobi[4];
            // Elementdaten
            static double v_rs[2];
            // Initialisieren
            MNulleVec(v_rs,2);

            //if (mmp_element_integration_method_maximum) ?? CB: was ist das
            if(1>0)
            {
               gp_r = (int)(point/nGauss);
               gp_s = point%nGauss;
               ur = MXPGaussPkt(nGauss, gp_r);
               us = MXPGaussPkt(nGauss, gp_s);
            }
            else
            {
               ur = 0.0;                          // Alle Groessen nur in Elementmitte betrachten
               us = 0.0;                          // Alle Groessen nur in Elementmitte betrachten
            }

            // Geschwindigkeitstransformation: a,b -> r,s
            //Calc2DElementJacobiMatrix(ind, 0., 0., invjac, &detjac);
            detjac = computeJacobian(1);          // order
            invjac = invJacobian;
            MKopierVec(invjac, jacobi, 4);
            M2Invertiere(jacobi);                 // Jacobi-Matrix
            MMultMatVec(jacobi, 2, 2, v, 2, v_rs, 2);

            if(MBtrgVec(v_rs, 2) > MKleinsteZahl)
            {
               // Upwind-Faktoren
               for(k=0;k<2;k++)
                  alpha[k] = -upwind_para * v_rs[k] / (MBtrgVec(v_rs, 2) + MKleinsteZahl);
            }

            //hier nur Methoden 1 + 2; 3 wird in CalcCoefLaplace erledigt
            if (upwind_meth == 1)
            {
               // Verschiebungen der Gausspunkte auf Element begrenzen
               scale = 1.;
               if (fabs(ur + alpha[0]) > 1.)
                  scale = MMin(scale, (1. - fabs(ur)) / fabs(alpha[0]));
               if (fabs(us + alpha[1]) > 1.)
                  scale = MMin(scale, (1. - fabs(us)) / fabs(alpha[1]));
               *rupw = ur + scale * alpha[0];
               *supw = us + scale * alpha[1];
            }
            else if (upwind_meth == 2)
            {
               // Verschiebungen auf -1<x<1 begrenzen
               *rupw = MRange(-1., ur + alpha[0], 1.);
               *supw = MRange(-1., us + alpha[1], 1.);
            }
         }
         break;
         case 3:                                  // Hexahedra
         {
            /* Elementgeometriedaten */
            static double detjac, *invjac, jacobi[9];
            /* Elementdaten */
            static double v_rst[3];
            // Initialisieren
            MNulleVec(v_rst,3);

            //if (mmp_element_integration_method_maximum) ?? CB: was ist das
            if(1>0)                               // CB: to do ??
            {
               gp_r = (int)(point/(nGauss*nGauss));
               gp_s = (point%(nGauss*nGauss));
               gp_t = gp_s%nGauss;
               gp_s /= nGauss;
               ur = MXPGaussPkt(nGauss, gp_r);
               us = MXPGaussPkt(nGauss, gp_s);
               ut = MXPGaussPkt(nGauss, gp_t);
            }
            else
            {
               ur = 0.0;                          // Alle Groessen nur in Elementmitte betrachten
               us = 0.0;                          // Alle Groessen nur in Elementmitte betrachten
               ut = 0.0;                          // Alle Groessen nur in Elementmitte betrachten
            }

            //Calc3DElementJacobiMatrix(ind, 0., 0., 0., invjac, &detjac);
            detjac = computeJacobian(1);          // order
            invjac = invJacobian;
            MKopierVec(invjac, jacobi, 9);
            M3Invertiere(jacobi);                 /* zurueck zur Jacobi-Matrix */
            MMultMatVec(jacobi, 3, 3, v, 3, v_rst, 3);

            if (MBtrgVec(v_rst, 3) > MKleinsteZahl)
            {
               /* Upwind-Faktoren */
               for (l = 0; l < 3; l++)
                  alpha[l] = -upwind_para * v_rst[l] / MBtrgVec(v_rst, 3) + MKleinsteZahl;
            }
            //hier nur Methoden 1 + 2; 3 wird in CalcCoefLaplace erledigt
            if (upwind_meth == 1)
            {
               // Verschiebungen der Gausspunkte auf Element begrenzen
               scale = 1.;
               if (fabs(ur + alpha[0]) > 1.)
                  scale = MMin(scale, (1. - fabs(ur)) / fabs(alpha[0]));
               if (fabs(us + alpha[1]) > 1.)
                  scale = MMin(scale, (1. - fabs(us)) / fabs(alpha[1]));
               if (fabs(ut + alpha[2]) > 1.)
                  scale = MMin(scale, (1. - fabs(ut)) / fabs(alpha[2]));
               *rupw = ur + scale * alpha[0];     // ist die reihenfolge hier richtig?
               *supw = us + scale * alpha[1];     // scale h?gt hier ja nur von dem letzten if ab..
               *tupw = ut + scale * alpha[2];
            }
            else if (upwind_meth == 2)
            {
               // Verschiebungen auf -1<x<1 begrenzen
               *rupw = MRange(-1., ur + alpha[0], 1.);
               *supw = MRange(-1., us + alpha[1], 1.);
               *tupw = MRange(-1., ut + alpha[2], 1.);
            }
         }
         break;
         case 4:                                  // Triangle
            break;
         case 5:                                  // Tedrahedra
            break;
         case 6:                                  // Prism
            break;
      }

      //test
      for(i=0;i<ele_dim;i++)
         cout << unit[i] << " ";
      cout << endl;
#endif
   }

   //SB4200
   /**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for advection matrix
   Programing:
   01/2005 WW/OK Implementation
   03/2005 WW Heat transport
   07/2005 WW Change for geometry element object
   09/2005 SB
   last modification:
   **************************************************************************/
   inline double CFiniteElementStd::CalCoefAdvection()
   {
      double val = 0.0;
      double dens_arg[3];                         //AKS
      //OK long Index = MeshElement->GetIndex();
      //----------------------------------------------------------------------
      switch(PcsType)
      {
         default:
            cout << "Fatal error in CalCoefAdvection: No valid PCS type" << std::endl;
            break;
         case L:                                  // Liquid flow
            break;
         case U:                                  // Unconfined flow
            break;
         case G:                                  // MB now Groundwater flow
            break;
         case T:                                  // Two-phase flow
            break;
         case C:                                  // Componental flow
            break;
         case H:                                  // heat transport
            if(FluidProp->density_model==14 && MediaProp->heat_diffusion_model==273 && cpl_pcs )
            {
               dens_arg[0]=interpolate(NodalValC1);
               dens_arg[1]=interpolate(NodalVal1)+T_KILVIN_ZERO;
               dens_arg[2]=Index;
               val = FluidProp->SpecificHeatCapacity(dens_arg)*FluidProp->Density(dens_arg);
            }
            else
            {
               val = FluidProp->SpecificHeatCapacity()*FluidProp->Density();
            }
            break;
         case M:                                  // Mass transport //SB4200
            val = 1.0*time_unit_factor;           //*MediaProp->Porosity(Index,pcs->m_num->ls_theta); // Porosity;
            break;
         case O:                                  // Liquid flow
            val = 1.0;
            break;
         case R:                                  // Richards
            break;
         case F:                                  // Fluid Momentum
            break;
         case A:                                  // Air (gas) flow
            val = 1.0/interpolate(NodalVal1);     // 1/p
            break;
      }
      return val;
   }
   /**************************************************************************
      GeoSys - Function: CalCoefStrainCouping

      Aufgabe:
            Calculate coefficient for StrainCouping matrix
      Programmaenderungen:
      01/2005   WW/OK    Erste Version
      07/2005 WW Change for geometry element object
   **************************************************************************/
   inline double CFiniteElementStd::CalCoefStrainCouping()
   {
      double val = 0.0;
      /*
      double r = unit[0];
      double s = unit[1];
      double t = unit[2];
      */

      switch(PcsType)
      {
         default:
            break;
         case L:                                  // Liquid flow
            //
            val = 1.0;
            break;
         case U:                                  // Unconfined flow
            break;
         case G:                                  // Groundwater
            break;
         case T:                                  // Two-phase flow
            break;
         case C:                                  // Componental flow
            break;
         case H:                                  // heat transport
            break;
         case M:                                  // Mass transport
            break;
         case O:                                  // Overland flow
            break;
         case R:                                  // Richard flow
            val = interpolate(NodalVal_Sat);      // Water saturation
         case A:
            break;
            break;
      }
      return val;
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcMass
      Aufgabe:
              Compute mass matrix, i.e. int (N.mat.N). Linear interpolation

      Programming:
      01/2005   WW
   02/2005 OK GEO factor
   01/2010 NW SUPG
   **************************************************************************/
   void CFiniteElementStd::CalcMass()
   {
      int i, j;                                   //OK411 k;
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;
      double fkt,mat_fac;
      // Material
      mat_fac = 1.0;
      double alpha[3], summand[8];
      double vel[3];                              //NW
      //  int indice = MeshElement->GetIndex();
      //  int phase = pcs->pcs_type_number;
      int upwind_method = pcs->m_num->ele_upwind_method;
      MNulleVec(alpha,3);
      MNulleVec(summand,8);

      if(PcsType==T)
      {
         if(upwind_method > 0)
         {
            // CB 11/07 this is to provide the velocity at the element center of gravity
            // call to his function here is also required for upwinding in CalcCoefLaplace
            Cal_Velocity_2();
            UpwindAlphaMass(alpha);               // CB 160507
         }
      }

      ElementValue* gp_ele = ele_gp_value[Index]; //NW

      //----------------------------------------------------------------------
      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         //if(PcsType==T)
         //{
         //  if((upwind_method == 1) || (upwind_method == 2))
         //      UpwindUnitCoord(phase, gp, indice); // phase 0
         //}
         ComputeShapefct(1);                      // Linear interpolation function
         // Material
         mat_fac = CalCoefMass();
         // if(Index < 0) cout << "mat_fac in CalCoeffMass: " << mat_fac << endl;
         // GEO factor
         mat_fac *= fkt;
         // Calculate mass matrix
         if(PcsType==T)
         {
            // upwinding: addiere SUPG-Summanden auf shapefct entsprechend Fkt. Mphi2D_SPG
            if(pcs->m_num->ele_upwind_method > 0)
               UpwindSummandMass(gp, gp_r, gp_s, gp_t, alpha, summand);
            for(i = 0; i < nnodes; i++)
            {
               for(j = 0; j < nnodes; j++)
                                                  // bei CT: phi * omega; phi beinh. uw-fakt.
                  (*Mass)(i,j) += mat_fac * (shapefct[i] + summand[i]) * shapefct[j];
            }
            //TEST OUTPUT
         }
         else
         {
            for (i = 0; i < nnodes; i++)
            {
               for (j = 0; j < nnodes; j++)
               {
                                                  //NW
                  if (pcs->m_num->ele_supg_method == 0)
                     if(j>i) continue;
                  (*Mass)(i,j) += mat_fac *shapefct[i]*shapefct[j];
               }
            }
            if (pcs->m_num->ele_supg_method > 0)  //NW
            {
               vel[0] = gp_ele->Velocity(0, gp);
               vel[1] = gp_ele->Velocity(1, gp);
               vel[2] = gp_ele->Velocity(2, gp);

               double tau = 0;
               CalcSUPGWeightingFunction(vel, gp, tau, weight_func);

               // tau*({v}[dN])^T*[N]
               for (i=0; i<nnodes; i++)
               {
                  for (j=0; j<nnodes; j++)
                     (*Mass)(i,j) += mat_fac*tau*weight_func[i]*shapefct[j];
               }
            }
         }                                        //end else
      }                                           // loop gauss points

                                                  //WW/CB //NW
      if(PcsType!=T && pcs->m_num->ele_supg_method==0)
      {
         for(i = 0; i < nnodes; i++)
         {
            for(j = 0; j < nnodes; j++)
            {
               if(j>i)
                  (*Mass)(i,j) = (*Mass)(j,i);
            }
         }
      }
      // Test Output
      //Mass->Write();
   }
   /**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2009 NW
   last modification:
   **************************************************************************/
   inline double CFiniteElementStd::CalcSUPGCoefficient(double*vel,int ip)
   {
      //--------------------------------------------------------------------
      // Collect following information to determine SUPG coefficient
      // + flow velocity
      // + diffusivity coefficient (scalar)
      // + characteristic element length
      // + (Peclet number)

      double v_mag = sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);
      // Characteristic element length
      double ele_len = CalcSUPGEffectiveElemenetLength(vel);
      // Diffusivity = (effective heat conductivity) / (fluid heat capacity)
      double diff = 0;
      if (PcsType==H)                             //heat
      {
         double *heat_conductivity_tensor = MediaProp->HeatConductivityTensor(MeshElement->GetIndex());

         if((FluidProp->density_model==14))
         {
            double dens_arg[3];                   //AKS
            int Index = MeshElement->GetIndex();
            dens_arg[0]=interpolate(NodalValC1);
            dens_arg[1]=interpolate(NodalVal1)+T_KILVIN_ZERO;
            dens_arg[2] =Index;
            diff = heat_conductivity_tensor[0] / (FluidProp->SpecificHeatCapacity(dens_arg)*FluidProp->Density(dens_arg));
         }
         else
         {
            diff = heat_conductivity_tensor[0] / (FluidProp->SpecificHeatCapacity()*FluidProp->Density());
         }
      }                                           //mass
      else if (PcsType==M)
      {
         double *advection_dispersion_tensor = MediaProp->MassDispersionTensorNew(ip);
         switch (pcs->m_num->ele_supg_method_diffusivity)
         {
            case 1:                               // min
            {
               double min_diff = advection_dispersion_tensor[0];
               for (int i=1; i<dim; i++)
               {
                  if (advection_dispersion_tensor[i]<min_diff) min_diff = advection_dispersion_tensor[i];
               }
               diff = min_diff;
            }
            break;
            case 2:                               // magnitude
            {
               double tmp_diff = 0.0;
               for (int i=0; i<dim; i++)
               {
                  tmp_diff = pow(advection_dispersion_tensor[i], 2.0);
               }
               diff = sqrt(tmp_diff);
            }
            break;
            default:                              //0 or any invalid number: max. in dispersion coefficient
            {
               double max_diff = advection_dispersion_tensor[0];
               for (int i=1; i<dim; i++)
               {
                  if (advection_dispersion_tensor[i]>max_diff) max_diff = advection_dispersion_tensor[i];
               }
               diff = max_diff;
            }
         }
      }

      //--------------------------------------------------------------------
      // Here calculates SUPG coefficient (tau)
      double tau = 0.0;
      switch (pcs->m_num->ele_supg_method)
      {
         case 1:
         {
            // this coefficient matches with the analytical solution in 1D stady state case
            double alpha = 0.5*v_mag*ele_len/diff;// 0.5*Pe
            double func = MLangevin(alpha);
            tau = 0.5*ele_len/v_mag*func;
         }
         break;
         case 2:
         {
            // taking into account time step
            tau = 1.0 / sqrt(pow(2.0/dt ,2.0)+pow(2.0*v_mag/ele_len,2.0));
         }
         break;
      }

      return tau;
   }
   /**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2009 NW
   last modification:
   **************************************************************************/
   inline void CFiniteElementStd::CalcSUPGWeightingFunction(double *vel, int ip, double &tau, double *v_dN)
   {
      if (pcs->m_num->ele_supg_method==0)
      {
         cout << "***Warning in CalcSUPGWeightingFunction(): SPUG is not selected" << endl;
         return;
      }

      // tau
      tau = CalcSUPGCoefficient(vel,ip);

      // {v}[dN]
      for (int i=0; i<nnodes; i++)
         v_dN[i] = 0.0;
      for (int i=0; i<nnodes; i++)
         for (int k=0; k<dim; k++)
            v_dN[i] += dshapefct[k*nnodes+i] * vel[k];
   }
   /**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2009 NW
   last modification:
   **************************************************************************/
   inline double CFiniteElementStd::CalcSUPGEffectiveElemenetLength(double *vel)
   {
      vel = vel;                                  //OK411
      double L = 0.0;
      switch (this->ele_dim)
      {
         case 1:
         {
            L = this->MeshElement->GetVolume();
         }
         break;
         case 2:
         case 3:
         {
            switch (pcs->m_num->ele_supg_method_length)
            {
               case 1:                            //min
               {
                  double min = MeshElement->GetEdge(0)->Length();
                  for (int i=1; i<MeshElement->GetEdgesNumber(); i++)
                  {
                     L = MeshElement->GetEdge(i)->Length();
                     if (L<min) min = L;
                  }
                  L = min;
               }
               break;
               case 2:                            //average
               {
                  double tmp_L=0.0;
                  for (int i=1; i<MeshElement->GetEdgesNumber(); i++)
                  {
                     tmp_L += MeshElement->GetEdge(i)->Length();
                  }
                  L = tmp_L/MeshElement->GetEdgesNumber();
               }
               break;
               case 3:                            // stream line length
               {
                  cout << "***Error: ele_supg_method_length <3> has not been supported yet." << endl;
               }
               break;
               default:                           //0 or any invalid number: max edge length
               {
                  double max = MeshElement->GetEdge(0)->Length();
                  for (int i=1; i<MeshElement->GetEdgesNumber(); i++)
                  {
                     L = MeshElement->GetEdge(i)->Length();
                     if (L>max) max = L;
                  }
                  L = max;
               }
               break;
            }
         }
         break;
      }
      return L;
   }

   //CB 090507
   /**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2007 CB
   last modification:
   **************************************************************************/
   inline void CFiniteElementStd::UpwindAlphaMass(double *alpha)
   {
      //Laufvariablen
      static long i, j;                           //, k, l;
      int no_phases;

      //static long *element_nodes;
      double gp[3], v_rst[3], v_tot[3];
      static double zeta;
      //static double *velovec, vg, v[2], vt[2], v_rs[2];
      //static double alpha_adv[3];
      double upwind_para;
      double upwind_meth;
      // Numerik
      zeta = 0.0;
      gp[0]=0.0;   gp[1]=0.0;   gp[2]=0.0;

      int ind = MeshElement->GetIndex();
      ElementValue* gp_ele = ele_gp_value[ind];
      //  TF ununsed: MshElemType::type eletyp = MeshElement->GetElementType();

      // Elementdaten und globale Modellparameter
      no_phases =(int)mfp_vector.size();
      //
      upwind_para = pcs->m_num->ele_upwinding;
      upwind_meth = pcs->m_num->ele_upwind_method;

      // alpha initialisieren
      MNulleVec(alpha, 3);
      // v initialisieren
      MNulleVec(v_tot,3);

      // get the velocities for phase 0
      v_tot[0] = gp_ele->Velocity(0, 0);
      v_tot[1] = gp_ele->Velocity(1, 0);
      v_tot[2] = gp_ele->Velocity(2, 0);
      // this would only give v at GP 0
      // but: if(PcsType==T) SetCenterGP(); // CB 11/2007
      // was set in Cal_Velo_2(); which is calculated,
      // when mass matrix is assembled
      // hence v is at element center of gravity
      // otherwise use following approximation:
      //for(i=0; i<nGaussPoints; i++) //element averaged v?
      //{
      //  v_tot[0] += gp_ele->Velocity(0, i)/nGaussPoints;
      //  v_tot[1] += gp_ele->Velocity(1, i)/nGaussPoints;
      //  v_tot[2] += gp_ele->Velocity(2, i)/nGaussPoints;
      //}

      // switch to next phase
      gp_ele = ele_gp_value[ind+(long)pcs->m_msh->ele_vector.size()];
      // get the velocities for phases 1 and add
      v_tot[0] += gp_ele->Velocity(0, 0);
      v_tot[1] += gp_ele->Velocity(1, 0);
      v_tot[2] += gp_ele->Velocity(2, 0);
      //for(i=0; i<nGaussPoints; i++) //element averaged v?
      //{
      //  v_tot[0] += gp_ele->Velocity(0, i)/nGaussPoints;
      //  v_tot[1] += gp_ele->Velocity(1, i)/nGaussPoints;
      //  v_tot[2] += gp_ele->Velocity(2, i)/nGaussPoints;
      //}
      for (i=0; i<3; i++)
         v_tot[i] *= time_unit_factor;

      //SetGaussPoint(point, gp_r, gp_s, gp_t);
      // velocity transformation a,b,c -> r,s,t
      computeJacobian(1);                         // order 1
      // multiply velocity vector with Jacobian matrix
      // Jacobi*v-->v_rst
      for (i=0; i<ele_dim; i++)
      {
         v_rst[i] = 0.0;
         for (j=0; j<ele_dim; j++)
            v_rst[i] += Jacobian[i*dim+j]*v_tot[j];
      }

      // Upwind-Factors
      if(MBtrgVec(v_rst, ele_dim) > MKleinsteZahl)// if(lengthOftheVector > tolerance)
      {
         for(i=0;i<ele_dim;i++)
            alpha[i] = -upwind_para * v_rst[i] / (MBtrgVec(v_rst, ele_dim) + MKleinsteZahl);
      }

#ifdef OLD_UPWINDING

      //test
      for(i=0;i<ele_dim;i++)
         cout << alpha[i] << " " ;
      cout << endl;

      switch(eletyp)
      {
         case 1:                                  // Line
         {
            // Elementgeometriedaten
            static double detjac, *invjac, jacobi[4];
            static double l[3];
            //invjac = GetElementJacobiMatrix(index, &detjac);
            //Calc1DElementJacobiMatrix(ind, invjac, &detjac);  //index = element id number
            detjac = computeJacobian(1);          // order
            invjac = invJacobian;
            //element_nodes = ElGetElementNodes(ind);

            MNulleVec(l,3);
            l[0] = X[1] - X[0];
            l[1] = Y[1] - Y[0];
            l[2] = Z[1] - Z[0];

            if (MBtrgVec(v_tot, 3) > MKleinsteZahl)
            {
               if (MSkalarprodukt(v_tot, l, 3) > 0.)
                  zeta = 1.;                      // upwind_para
               else zeta = -1.;                   //-upwind_para
            }

            //aus RF 3.5.06 CT 1D elements: {
            //// detjac = A*L/2
            //vorfk = porosity * detjac * Mdrittel;
            //// Massenmatrix mit SUPG ohne Zeitanteile
            //mass[0] = (2.0 + 1.5 * mms_upwind_parameter * zeta) * vorfk;
            //mass[1] = (1.0 + 1.5 * mms_upwind_parameter * zeta) * vorfk;
            //mass[2] = (1.0 - 1.5 * mms_upwind_parameter * zeta) * vorfk;
            //mass[3] = (2.0 - 1.5 * mms_upwind_parameter * zeta) * vorfk; // }

            // Upwind-Faktor Fully upwinding
            //alpha[0]     = m_pcs->m_num->ele_upwinding * zeta;
            //alpha_adv[0] = m_pcs->m_num->ele_upwinding * zeta;
            alpha[0]     = 1.0 * zeta;            //??
            //alpha_adv[0] = 1.0 * zeta;
            // Advection upwinding
            //if (MTM2_upwind_method == 2) alpha_adv[0] = ele_upwinding * zeta;   /
         }
         break;
         case 2:                                  // Quadrilateral
         {
            // Elementgeometriedaten
            static double detjac, *invjac, jacobi[4];
            // Elementdaten
            static double v_rs[3];

            // Geschwindigkeitstransformation: a,b -> r,s
            //Calc2DElementJacobiMatrix(ind, 0., 0., invjac, &detjac);
            detjac = computeJacobian(1);          // order
            invjac = invJacobian;
            MKopierVec(invjac, jacobi, 4);
            M2Invertiere(jacobi);                 /* Jacobi-Matrix */
            MMultMatVec(jacobi, 2, 2, v_tot, 2, v_rs, 2);

            if(MBtrgVec(v_rs, 2) > MKleinsteZahl)
            {
               // Upwind-Faktoren
               for(k=0;k<2;k++)
               {
                  alpha[k] = -upwind_para * v_rs[k] / (MBtrgVec(v_rs, 2) + MKleinsteZahl);
               }
            }
         }
         break;
         case 3:                                  // Hexahedra
         {
            /* Elementgeometriedaten */
            static double *invjac, jacobi[9], detjac;
            /* Elementdaten */
            //static double v_rst[3];

            if (MBtrgVec(v_tot, 3) > MKleinsteZahl)
            {
               /* Geschwindigkeitstransformation: x,y,z -> r,s,t */
               //Calc3DElementJacobiMatrix(ind, 0., 0., 0., invjac, &detjac);
               detjac = computeJacobian(1);       // order
               invjac = invJacobian;
               MKopierVec(invjac, jacobi, 9);
               M3Invertiere(jacobi);              /* Jacobi-Matrix */
               MMultMatVec(jacobi, 3, 3, v_tot, 3, v_rst, 3);

               /* Upwind-Faktoren */
               for (l = 0; l < 3; l++)
               {
                  alpha[l] = -upwind_para * v_rst[l] / (MBtrgVec(v_rst, 3) + MKleinsteZahl);
               }
            }
         }
         break;
         case 4:                                  // Triangle
            break;
         case 5:                                  // Tedrahedra
            break;
         case 6:                                  // Prism
            break;
      }
      //test
      for(i=0;i<ele_dim;i++)
         cout << alpha[i] << " " ;
      cout << endl;
#endif

   }

   //CB 160507
   /**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2007 CB
   last modification:
   **************************************************************************/
   inline void CFiniteElementStd::UpwindSummandMass(const int gp, int& gp_r, int& gp_s, int& gp_t, double *alpha, double *summand)

   {
      int i, k;
      //
      GetGaussData(gp, gp_r, gp_s, gp_t);         // this sets unit[] to standard values
      GradShapeFunction(dshapefct, unit);
      for(i=0;i<nnodes;i++)
      {
         summand[i] = 0.0;
         for(k=0;k<dim;k++)
            summand[i] += dshapefct[nnodes*k+i]*alpha[k];
         //summand[i] /= (double)nnodes;
      }

#ifdef OLD_UPWINDING

      double u1, u2, u3;
      u1 = u2 = u3 = 0;

      MshElemType::type eletyp = MeshElement->GetElementType();
      switch(eletyp)
      {
         case MshElemType::LINE:
         {
            // Line
            gp_r = gp;
            u1 = MXPGaussPkt(nGauss, gp_r);
            summand[0] = + alpha[0]*(1+u1);       //CB: ?? hab ich mir so gedacht
            summand[1] = - alpha[0]*(1-u1);       //CB: ?? hab ich mir so gedacht
            for(i=0;i<2;i++)
               summand[i] *= 0.5;
         }
         break;
         case MshElemType::QUAD:                  // Quadrilateral
         {
            gp_r = (int)(gp/nGauss);
            gp_s = gp%nGauss;
            u1 = MXPGaussPkt(nGauss, gp_r);
            u2 = MXPGaussPkt(nGauss, gp_s);
            // derived from MPhi2D_SUPG
            summand[0] = + alpha[0]*(1+u2) + alpha[1]*(1+u1);
            summand[1] = - alpha[0]*(1+u2) + alpha[1]*(1-u1);
            summand[2] = - alpha[0]*(1-u2) - alpha[1]*(1-u1);
            summand[3] = + alpha[0]*(1-u2) - alpha[1]*(1+u1);
            for(i=0;i<4;i++)
               summand[i] *= 0.25;
         }
         break;
         case MshElemType::HEXAHEDRON:            // Hexahedra
         {
            gp_r = (int)(gp/(nGauss*nGauss));
            gp_s = (gp%(nGauss*nGauss));
            gp_t = gp_s%nGauss;
            gp_s /= nGauss;
            u1 = MXPGaussPkt(nGauss, gp_r);
            u2 = MXPGaussPkt(nGauss, gp_s);
            u3 = MXPGaussPkt(nGauss, gp_t);
            // derived from MPhi3D_SUPG
            summand[0] = + alpha[0]*(1+u2)*(1+u3) + alpha[1]*(1+u1)*(1+u3) + alpha[2]*(1+u1)*(1+u2);
            summand[1] = - alpha[0]*(1+u2)*(1+u3) + alpha[1]*(1-u1)*(1+u3) + alpha[2]*(1-u1)*(1+u2);
            summand[2] = - alpha[0]*(1-u2)*(1+u3) - alpha[1]*(1-u1)*(1+u3) + alpha[2]*(1-u1)*(1-u2);
            summand[3] = + alpha[0]*(1-u2)*(1+u3) - alpha[1]*(1+u1)*(1+u3) + alpha[2]*(1+u1)*(1-u2);
            summand[4] = + alpha[0]*(1+u2)*(1-u3) + alpha[1]*(1+u1)*(1-u3) - alpha[2]*(1+u1)*(1+u2);
            summand[5] = - alpha[0]*(1+u2)*(1-u3) + alpha[1]*(1-u1)*(1-u3) - alpha[2]*(1-u1)*(1+u2);
            summand[6] = - alpha[0]*(1-u2)*(1-u3) - alpha[1]*(1-u1)*(1-u3) - alpha[2]*(1-u1)*(1-u2);
            summand[7] = + alpha[0]*(1-u2)*(1-u3) - alpha[1]*(1+u1)*(1-u3) - alpha[2]*(1+u1)*(1-u2);
            for(i=0;i<8;i++)
               summand[i] *= 0.125;
         }
         break;
         case MshElemType::TRIANGLE:              // Triangle
         {
            //SamplePointTriHQ(gp, unit);
         }
         break;
         case MshElemType::TETRAHEDRON:           // Tedrahedra
         {
            //SamplePointTet5(gp, unit);
         }
         break;
         case MshElemType::PRISM:                 // Prism
         {
            gp_r = gp%nGauss;
            gp_s = (int)(gp/nGauss);
            gp_t = (int)(nGaussPoints/nGauss);
            //u1 = MXPGaussPktTri(nGauss,gp_r,0); //femlib.cpp statt mathlib.cpp, nicht verfügbar?
            //u2 = MXPGaussPktTri(nGauss,gp_r,1);
            //u3 = MXPGaussPkt(gp_t,gp_s);
         }
         break;
      }
#endif

      return;
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcMass2
      Programming:
      02/2007   WW
   **************************************************************************/
   void CFiniteElementStd::CalcMass2()
   {
      int i, j,in,jn;
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;
      double fkt,mat_fac;
      // Material
      int dof_n = 2;
      mat_fac = 1.0;
      //----------------------------------------------------------------------
      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         ComputeShapefct(1);                      // Linear interpolation function
         for(in=0; in<dof_n; in++)
         {
            for(jn=0; jn<dof_n; jn++)
            {
               // Material
               mat_fac = CalCoefMass2(in*dof_n+jn);
               mat_fac *= fkt;
               // Calculate mass matrix
               for (i = 0; i < nnodes; i++)
               {
                  for (j = 0; j < nnodes; j++)
                     (*Mass2)(i+in*nnodes,j+jn*nnodes) += mat_fac *shapefct[i]*shapefct[j];
               }
            }
         }
      }
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcMassPSGLOBAL
      Programming:
      03/2009   PCH
   **************************************************************************/
   void CFiniteElementStd::CalcMassPSGLOBAL()
   {
      int i, j,in,jn;
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;
      double fkt,mat_fac;
      // Material
      int dof_n = 2;
      mat_fac = 1.0;
      //----------------------------------------------------------------------
      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         ComputeShapefct(1);                      // Linear interpolation function
         for(in=0; in<dof_n; in++)
         {
            for(jn=0; jn<dof_n; jn++)
            {
               // Material
               mat_fac = CalCoefMassPSGLOBAL(in*dof_n+jn);
               mat_fac *= fkt;
               // Calculate mass matrix
               for (i = 0; i < nnodes; i++)
               {
                  for (j = 0; j < nnodes; j++)
                     (*Mass2)(i+in*nnodes,j+jn*nnodes) += mat_fac *shapefct[i]*shapefct[j];
               }
            }
         }
      }
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcLumpedMass
      Aufgabe:
              Compute lumped mass matrix, i.e. int (N.mat.N). Linear interpolation

      Programming:
      01/2005   WW
      02/2005 OK GEO factor
      02/2007   WW Multi-phase flow
      05/2007   WW Axismmetry volume
   01/2010   NW geometrical area
   **************************************************************************/
   void CFiniteElementStd::CalcLumpedMass()
   {
      int i, gp_r, gp_s, gp_t;
      double factor, vol=0.0;
      gp=0;
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      // Initialize
      (*Mass) = 0.0;
      // Volume
      if(axisymmetry)
      {                                           // This calculation should be done in CompleteMesh.
         // However, in order not to destroy the concise of the code,
         // it is put here. Anyway it is computational cheap. WW
         vol = 0.0;
         for (gp = 0; gp < nGaussPoints; gp++)
         {
            //---------------------------------------------------------
            //  Get local coordinates and weights
            //  Compute Jacobian matrix and its determinate
            //---------------------------------------------------------
            vol += GetGaussData(gp, gp_r, gp_s, gp_t);
         }
      }
      else
                                                  //NW multiply geo_area
         vol = MeshElement->GetVolume()*MeshElement->GetFluxArea();
      // Center of the reference element
      SetCenterGP();
      factor = CalCoefMass();
      pcs->timebuffer = factor;                   // Tim Control "Neumann"
      factor *= vol/(double)nnodes;
      for (i=0; i<nnodes; i++)
         (*Mass)(i,i) =  factor;
      //
#ifdef otherLumpedMass
      int i, j;
      int gp_r=0, gp_s=0, gp_t=0;
      double fkt;
      //----------------------------------------------------------------------
      for (i=0; i<nnodes; i++)
      {
         for(j=0; j<ele_dim; j++)
            x2buff[j] = nodes_xyz[j*nnodes+i];
         UnitCoordinates(x2buff);
         fkt = GetGaussData(i, gp_r, gp_s, gp_t)*CalCoefMass();
         (*Mass)(i+in*nnodes,i+jn*nnodes) += fkt ;
      }
      //----------------------------------------------------------------------
#endif
      //TEST OUT
      // Mass->Write();
   }

   ///
   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcLumpedMass2
      Programming:
      02/2007   WW Multi-phase flow
   **************************************************************************/
   void CFiniteElementStd::CalcLumpedMass2()
   {
      int i, in, jn, gp_r, gp_s, gp_t;
      double factor, vol=0.0;
      int dof_n = 2;
      //----------------------------------------------------------------------
      // Volume
      if(axisymmetry)
      {                                           // This calculation should be done in CompleteMesh.
         // However, in order not to destroy the concise of the code,
         // it is put here. Anyway it is computational cheap. WW
         vol = 0.0;
         for (gp = 0; gp < nGaussPoints; gp++)
         {
            //---------------------------------------------------------
            //  Get local coordinates and weights
            //  Compute Jacobian matrix and its determinate
            //---------------------------------------------------------
            vol += GetGaussData(gp, gp_r, gp_s, gp_t);
         }
      }
      else
         vol = MeshElement->GetVolume();
      //----------------------------------------------------------------------
      // Initialize
      (*Mass2) = 0.0;
      // Center of the reference element
      SetCenterGP();
      for(in=0; in<dof_n; in++)
      {
         for(jn=0; jn<dof_n; jn++)
         {
            // Factor
            factor = CalCoefMass2(in*dof_n+jn);
            pcs->timebuffer = factor;             // Tim Control "Neumann"
            // Volume
            factor *= vol/(double)nnodes;
            for (i=0; i<nnodes; i++)
               (*Mass2)(i+in*nnodes,i+jn*nnodes) = factor;
         }
      }
      //TEST OUT
      // Mass2->Write();
   }
   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcLumpedMassPSGLOBAL
      Programming:
      03/2009   PCH PS_GLOBAL for Multi-phase flow
   **************************************************************************/
   void CFiniteElementStd::CalcLumpedMassPSGLOBAL()
   {
      int i, in, jn, gp_r, gp_s, gp_t;
      double factor, vol=0.0;
      int dof_n = 2;
      //----------------------------------------------------------------------
      // Volume
      if(axisymmetry)
      {                                           // This calculation should be done in CompleteMesh.
         // However, in order not to destroy the concise of the code,
         // it is put here. Anyway it is computational cheap. WW
         vol = 0.0;
         for (gp = 0; gp < nGaussPoints; gp++)
         {
            //---------------------------------------------------------
            //  Get local coordinates and weights
            //  Compute Jacobian matrix and its determinate
            //---------------------------------------------------------
            vol += GetGaussData(gp, gp_r, gp_s, gp_t);
         }
      }
      else
         vol = MeshElement->GetVolume();
      //----------------------------------------------------------------------
      // Initialize
      (*Mass2) = 0.0;
      // Center of the reference element
      SetCenterGP();
      for(in=0; in<dof_n; in++)
      {
         for(jn=0; jn<dof_n; jn++)
         {
            // Factor
            factor = CalCoefMassPSGLOBAL(in*dof_n+jn);
            pcs->timebuffer = factor;             // Tim Control "Neumann"
            // Volume
            factor *= vol/(double)nnodes;
            for (i=0; i<nnodes; i++)
               (*Mass2)(i+in*nnodes,i+jn*nnodes) = factor;
         }
      }
      //TEST OUT
      //  Mass2->Write();
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcStorage
      Aufgabe:
              Compute mass matrix, i.e. int (N.mat.N). Linear interpolation

      Programming:
      01/2005   WW
   02/2005 OK GEO factor
   **************************************************************************/
   void CFiniteElementStd::CalcStorage()
   {
      int i, j;
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;
      double fkt,mat_fac;
      // Material
      mat_fac = 1.0;
      //----------------------------------------------------------------------
      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         ComputeShapefct(1);                      // Linear interpolation function
         // Material
         mat_fac = CalCoefStorage();
         // GEO factor
         fkt *= mat_fac;
         // Calculate mass matrix
         for (i = 0; i < nnodes; i++)
            for (j = 0; j < nnodes; j++)
               (*Storage)(i,j) += fkt *shapefct[i]*shapefct[j];

      }
      //TEST OUTPUT
      //  if(Index == 195){cout << "Storage Matrix: " << endl; Storage->Write(); }
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcContent
      Aufgabe:
              Compute Content matrix, i.e. int (N.mat.N). Linear interpolation

      Programming:
      01/2005   WW
   02/2005 OK GEO factor
   **************************************************************************/
   void CFiniteElementStd::CalcContent()
   {
      int i, j;
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;
      double fkt,mat_fac;
      // Material
      mat_fac = 1.0;
      //----------------------------------------------------------------------
      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         ComputeShapefct(1);                      // Linear interpolation function
         // Material
         mat_fac = CalCoefContent();
         // GEO factor
         fkt *= mat_fac;
         // Calculate mass matrix
         for (i = 0; i < nnodes; i++)
            for (j = 0; j < nnodes; j++)
               (*Content)(i,j) += fkt *shapefct[i]*shapefct[j];

      }
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcLaplace
      Aufgabe:
              Compute mass matrix, i.e. int (gradN.mat.gradN). Linear interpolation

      Programming:
      01/2005   WW
      02/2005 OK GEO factor
      02/2007 WW Multi-phase
       03/2009 PCH PS_GLOBAL for Multiphase flow
   **************************************************************************/
   void CFiniteElementStd::CalcLaplace()
   {
      int i, j, k, l, in, jn;

      // ---- Gauss integral
      int gp_r=0, gp_s=0, gp_t;
      gp_t = 0;
      double fkt, water_depth;
      int dof_n = 1;

      int in_times_nnodes, jn_times_nnodes, i_plus_in_times_nnodes, j_plus_jn_times_nnodes, dim_times_k, dim_times_k_plus_l, nnodes_plus_j, l_times_nnodes_plus_j, k_times_nnodes_plus_i;
      double fkt_times_dshapefct__k_times_nnodes_plus_i__;

      if(PcsType==V || PcsType==P) dof_n = 2;     // 03.03 2009 PCH

      //----------------------------------------------------------------------
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         //---------------------------------------------------------
         // Compute geometry
         ComputeGradShapefct(1);                  // Linear interpolation function
         // Calculate mass matrix
         water_depth = 1.0;
         // The following "if" is done by WW
         if(PcsType == G && MediaProp->unconfined_flow_group == 1 && MeshElement->ele_dim == 2 && !pcs->m_msh->hasCrossSection())
         {
            water_depth = 0.0;
            for(i=0; i< nnodes; i++)
               water_depth += (pcs->GetNodeValue(nodes[i],idx1) -Z[i])*shapefct[i];
         }
         fkt *= water_depth;
         //---------------------------------------------------------

         for (in_times_nnodes = 0,in = 0; in < dof_n; in++, in_times_nnodes+=nnodes)
         {

            for (jn_times_nnodes = 0,jn = 0; jn < dof_n; jn++, jn_times_nnodes+=nnodes)
            {
               // Material
               if(dof_n==1)
                  CalCoefLaplace(false,gp);
               else if (dof_n==2)
               {
                  if (PcsType==V)
                     CalCoefLaplace2(false,in*dof_n+jn);
                  else if (PcsType==P)
                     CalCoefLaplacePSGLOBAL(false,in*dof_n+jn);
               }
               //---------------------------------------------------------

               for (i = 0, i_plus_in_times_nnodes = in_times_nnodes; i < nnodes; i++, i_plus_in_times_nnodes++)
               {
                  //i_plus_in_times_nnodes = i + in_times_nnodes;

                  for (j = 0,j_plus_jn_times_nnodes = jn_times_nnodes,
                     nnodes_plus_j = nnodes; j < nnodes; j++, j_plus_jn_times_nnodes++, nnodes_plus_j++)
                  {
                     for (k = 0, dim_times_k=0,k_times_nnodes_plus_i=i; k < dim; k++, dim_times_k+=dim,k_times_nnodes_plus_i+=nnodes)
                     {
                        fkt_times_dshapefct__k_times_nnodes_plus_i__ = fkt * dshapefct[k_times_nnodes_plus_i];
                        for(l=0, dim_times_k_plus_l = dim_times_k,l_times_nnodes_plus_j = j; l< dim; l++, dim_times_k_plus_l++,l_times_nnodes_plus_j+=nnodes)
                        {

                           (*Laplace)(i_plus_in_times_nnodes,j_plus_jn_times_nnodes) += fkt_times_dshapefct__k_times_nnodes_plus_i__ \
                              * mat[dim_times_k_plus_l] * dshapefct[l_times_nnodes_plus_j];
                           /*
                                             (*Laplace)(i+in*nnodes,j+jn*nnodes) += fkt * dshapefct[k*nnodes+i] \
                                                * mat[dim*k+l] * dshapefct[l*nnodes+j];
                                             if(Index < 10) {cout << " i, j, k, l, nnodes, dim: " << i << ", " << j << ", " << k << ", " << l << ", " << nnodes << ", " << dim << ". fkt, dshapefct[k*nnodes+i], mat[dim*k+l], dshapefct[l*nnodes+j]: ";
                                             cout << fkt << ", " << dshapefct[k*nnodes+i] << ", " << mat[dim*k+l] << ", " << dshapefct[l*nnodes+j] << endl;}
                           */
                        }
                     }
                  }                               // j: nodes
               }                                  // i: nodes
            }                                     // dof j
         }                                        // dof i
      }
      //TEST OUTPUT
      // Laplace->Write();
   }
   /**************************************************************************
   FEMLib-Method:
   10/2006 YD Implementation
   01/2007 WW Fundamental changes
   **************************************************************************/
   void CFiniteElementStd:: Assemble_DualTransfer()
   {
      int i,j;
      int gp_r=0, gp_s=0, gp_t=0;
      double W, fkt,mat_fac = 0.;
#if defined(NEW_EQS)
      CSparseMatrix *A = NULL;                    //WW
      if(m_dom)
         A = m_dom->eqs->A;
      else
         A = pcs->eqs_new->A;
#endif

      //Inintialize
      //-------------------------- WW
      W = pcs->continuum_vector[pcs->GetContinnumType()];
      //
      for(i=0;i<nnodes;i++)
      {
                                                  // Pressure 1
         NodalVal3[i] = pcs->GetNodeValue(nodes[i], idx1);
                                                  // Pressure 2
         NodalVal4[i] = pcs->GetNodeValue(nodes[i], idxp21);
      }
      (*Advection) = 0.0;
      //---------------------------------------------------------
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determination
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         mat_fac = CalcCoefDualTransfer();
         mat_fac *= fkt;
         // Material
         ComputeShapefct(1);                      // Linear interpolation function
         // Calculate mass matrix
         for (i = 0; i < nnodes; i++)
         {
            for (j = 0; j < nnodes; j++)
               (*Advection)(i,j) += mat_fac*shapefct[i]*shapefct[j];
         }
      }
      // Add local matrix to global matrix
      // 15.02.2007 WW
      long cshift = pcs->m_msh->GetNodesNumber(false);
      double fm = 1.0/W;
      //
      if(pcs->continuum == 0)
      {
         double ff = 1.0/(1.0-W);
         if(MediaProp->transfer_coefficient<0.0)  // for LBNL
            ff = 1.0;
         for(int i=0;i<nnodes;i++)
         {
            for(int j=0;j<nnodes;j++)
            {
#ifdef NEW_EQS
               (*A)(eqs_number[i], eqs_number[j]+cshift) += -fm*(*Advection)(i,j);
               (*A)(eqs_number[i]+cshift, eqs_number[j]) += -ff*(*Advection)(i,j);
#else
               MXInc(eqs_number[i], eqs_number[j]+cshift, -fm*(*Advection)(i,j));
               MXInc(eqs_number[i]+cshift, eqs_number[j], -ff*(*Advection)(i,j));
#endif
            }
         }
      }
      else
      {
         if(MediaProp->transfer_coefficient<0.0)  // for LBNL
            fm = 1.0;
      }
      //
      (*Advection) *= fm;
      (*Laplace) += (*Advection);
      //
      //-------------------------- WW
   }
   /**************************************************************************
   FEMLib-Method:
   10/2006 YD Implementation
   01/2007 WW Fundamental changes
   **************************************************************************/
   inline double  CFiniteElementStd::CalcCoefDualTransfer()
   {
      double Sm=0.0, Sf=0.0, ExFac=0.0;
      double pm=0.0, pf=0.0, matrix_conductivity, val=0;
      //double* permeability;
      double *permeability = NULL;
      //-------------------------------------------WW
      CMediumProperties *m_matrix = NULL;
      CMediumProperties *f_matrix = NULL;
      if(pcs->GetContinnumType()== 0)
      {
         m_matrix = MediaProp;
         f_matrix = MediaProp1;
      }
      else                                        // fracture //WW
      {
         m_matrix = MediaProp1;
         f_matrix = MediaProp;
      }
      //-------------------------------------------WW
      switch(PcsType)
      {
         default:
            break;
         case R:
            pm = interpolate(NodalVal3);
            pf = interpolate(NodalVal4);
                                                  // Matrix
            Sm = m_matrix->SaturationCapillaryPressureFunction(-pm,0);
                                                  // Fracture
            Sf = f_matrix->SaturationCapillaryPressureFunction(-pf,0);
            permeability = m_matrix->PermeabilityTensor(Index);
            ExFac = m_matrix->transfer_coefficient;
            // Dual by van Genuchten
            if(ExFac>0.0)
               matrix_conductivity = 0.5*(m_matrix->PermeabilitySaturationFunction(Sm,0)
                  +m_matrix->PermeabilitySaturationFunction(Sf,0))
                  /FluidProp->Viscosity();

            else                                  // by LBNL. WW
            {
               double Sf_e = (Sf - f_matrix->saturation_res[phase])/(f_matrix->saturation_max[phase]
                  -f_matrix->saturation_res[phase]);
               matrix_conductivity = Sf_e*m_matrix->PermeabilitySaturationFunction(Sm,0)\
                  /FluidProp->Viscosity();
               ExFac *= -1.0;
            }
            //
            val = time_unit_factor*permeability[0]*matrix_conductivity*ExFac;
            break;
            //---------------------------------------------------------
         case H:

            break;
      }
      return val;
   }

   //SB4200
   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcAdvection
      Aufgabe:  Calculate the advection matrix

      Programming:
      01/2005   WW
      02/2005   OK GEO factor
      09/2005   SB - adapted to advection
      03/2007   WW - Fluid advection with multiphase flow
      05/2008   WW - General densty for multiphase flow
   01/2010   NW - SUPG
   **************************************************************************/
   void CFiniteElementStd::CalcAdvection()
   {
      int i, j, k;
      int gp_r=0, gp_s=0, gp_t;
      double fkt,mat_factor = 0.0;
      double vel[3], dens_aug[3];
      CFluidProperties *m_mfp_g = NULL;
      bool multiphase = false;
                                                  //18.02.2008, 04.09.2008 WW
      if(!cpl_pcs&&(pcs->type!=2)&&(pcs->type!=5)) return;
      if(cpl_pcs&&cpl_pcs->type==1212)
      {
         multiphase = true;
         m_mfp_g =  mfp_vector[1];
         GasProp = MFPGet("GAS");
      }
      ElementValue* gp_ele = ele_gp_value[Index];
      CRFProcess* pcs_fluid_momentum = PCSGet("FLUID_MOMENTUM");

      //Initial values
      gp_t = 0;
      (*Advection)=0.0;

      //----------------------------------------------------------------------
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         mat_factor = CalCoefAdvection();         //T
         //---------------------------------------------------------
         // Compute geometry
         ComputeGradShapefct(1);                  // Linear interpolation function....dNJ-1....var dshapefct
         ComputeShapefct(1);                      // Linear interpolation N....var shapefct
         //---------------------------------------------------------
         //Velocity
         vel[0] = mat_factor*gp_ele->Velocity(0, gp);
         vel[1] = mat_factor*gp_ele->Velocity(1, gp);
         vel[2] = mat_factor*gp_ele->Velocity(2, gp);
         if(multiphase)                           //02/2007 WW
         {
            PG2=interpolate(NodalVal_p2);
            PG=  interpolate(NodalValC1);
            TG=interpolate(NodalVal1)+T_KILVIN_ZERO;
            rhow=FluidProp->Density();
            rho_gw = FluidProp->vaporDensity(TG)*exp(-PG*COMP_MOL_MASS_WATER/(rhow*GAS_CONSTANT*TG));
            p_gw = rho_gw*GAS_CONSTANT*TG/COMP_MOL_MASS_WATER;
            dens_aug[0] = PG2-p_gw;
            dens_aug[1] = TG;
                                                  // 29.05.2008. WW/ 2 Dec 2010 AKS
            rho_g = rho_gw + GasProp->Density(dens_aug);
            mat_factor = rho_g*m_mfp_g->SpecificHeatCapacity();
            vel[0] += mat_factor*gp_ele->Velocity_g(0, gp);
            vel[1] += mat_factor*gp_ele->Velocity_g(1, gp);
            vel[2] += mat_factor*gp_ele->Velocity_g(2, gp);
         }
         // Velocity by Fluid_Momentum - 13.11.2009  PCH
         if(pcs_fluid_momentum)
         {
            CRFProcess* m_pcs = pcs_fluid_momentum;

            vel[0] = mat_factor*m_pcs->GetElementValue(index, m_pcs->GetElementValueIndex("VELOCITY1_X")+1);
            vel[1] = mat_factor*m_pcs->GetElementValue(index, m_pcs->GetElementValueIndex("VELOCITY1_Y")+1);
            vel[2] = mat_factor*m_pcs->GetElementValue(index, m_pcs->GetElementValueIndex("VELOCITY1_Z")+1);
         }

         for (i = 0; i< nnodes; i++)
         {
            for (j = 0; j < nnodes; j++)
               for (k = 0; k < dim; k++)
                  (*Advection)(i,j) += fkt*shapefct[i]*vel[k]
                     *dshapefct[k*nnodes+j];
         }

         if (pcs->m_num->ele_supg_method > 0)     //NW
         {
            vel[0] = gp_ele->Velocity(0, gp);
            vel[1] = gp_ele->Velocity(1, gp);
            vel[2] = gp_ele->Velocity(2, gp);
            if(pcs_fluid_momentum)
            {
               CRFProcess* m_pcs = pcs_fluid_momentum;

               vel[0] = m_pcs->GetElementValue(index, m_pcs->GetElementValueIndex("VELOCITY1_X")+1);
               vel[1] = m_pcs->GetElementValue(index, m_pcs->GetElementValueIndex("VELOCITY1_Y")+1);
               vel[2] = m_pcs->GetElementValue(index, m_pcs->GetElementValueIndex("VELOCITY1_Z")+1);
            }

            double tau = 0;
            CalcSUPGWeightingFunction(vel, gp, tau, weight_func);

            //Calculate mat_factor*tau*({v}[dN])^T*({v}[dN])
            for (i=0; i<nnodes; i++)
            {
               for (j=0; j<nnodes; j++)
                  (*Advection)(i,j) += fkt*mat_factor*tau*weight_func[i]*weight_func[j];
            }
         }
      }
      //TEST OUTPUT
      //  cout << "Advection Matrix: " << endl; Advection->Write();
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcAdvection
      Aufgabe:  Calculate the advection matrix

      Programming:
      12/2005   WW
   ***************************************************************************/
   void CFiniteElementStd::CalcRHS_by_ThermalDiffusion()
   {
      int i, j, k;
      // ---- Gauss integral
      int gp_r=0, gp_s=0, gp_t;
      gp = 0;
      double fkt;
      double Dv = 0.0;
      double Dtv = 0.0;
      double poro = 0.0;
      double tort = 0.0;
      double humi = 1.0;
      double rhov = 0.0;
      double drdT = 0.0;
      double beta = 0.0;
      // 12.12.2007 WW
      long cshift = 0;
      if(pcs->dof>1)
         cshift = NodeShift[pcs->continuum];

      (*Laplace) = 0.0;
      (*Mass) = 0.0;
      //----------------------------------------------------------------------
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

         //---------------------------------------------------------
         // Compute geometry
         ComputeGradShapefct(1);                  // Linear interpolation function
         //---------------------------------------------------------
         // Material
         //	  if(FluidProp->diffusion_model=273)
         //    {
         ComputeShapefct(1);
         double rhow = FluidProp->Density();
         PG = interpolate(NodalVal1);
         TG = interpolate(NodalValC)+T_KILVIN_ZERO;
                                                  //WW
         Sw = MediaProp->SaturationCapillaryPressureFunction(-PG,0);
         poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
         tort = MediaProp->TortuosityFunction(Index,unit,pcs->m_num->ls_theta);
         beta = poro*MediaProp->StorageFunction(Index,unit,pcs->m_num->ls_theta) *Sw;
         //Rv = GAS_CONSTANT;
         humi = exp(PG/(GAS_CONSTANT_V*TG*rhow));
         Dv = tort*(1.0-Sw)*poro*2.16e-5*pow(TG/T_KILVIN_ZERO, 1.8);
         rhov = humi*FluidProp->vaporDensity(TG);
         drdT= (FluidProp->vaporDensity_derivative(TG)*humi-rhov*PG/(GAS_CONSTANT_V*rhow*pow(TG, 2.0)))/rhow;
         Dtv = time_unit_factor*Dv*drdT;

         //    }
         //---------------------------------------------------------
         // Calculate a Laplace
         for (i = 0; i < nnodes; i++)
         {
            for (j = 0; j < nnodes; j++)
            {
               if(j>i) continue;
               for (k = 0; k < dim; k++)
                  (*Laplace)(i,j) +=   fkt*Dtv*dshapefct[k*nnodes+i]
                     * dshapefct[k*nnodes+j];
               (*Mass)(i,j) += fkt*poro*(beta+(1.0-Sw)*drdT)*shapefct[i]*shapefct[j];
            }
         }
      }
      // Symmetry
      for (i = 0; i < nnodes; i++)
      {
         for (j = 0; j < nnodes; j++)
         {
            if(j<=i) continue;
            (*Laplace)(i,j) = (*Laplace)(j,i);
         }
      }
      cshift += NodeShift[problem_dimension_dm];
      for (i = 0; i < nnodes; i++)
      {
         for (j = 0; j < nnodes; j++)
         {
            (*RHS)(i) -= (*Laplace)(i,j)*(NodalValC[j]+T_KILVIN_ZERO);
            (*RHS)(i) += (*Mass)(i,j)*(NodalValC1[j]-NodalValC[j])/dt;
         }
         eqs_rhs[cshift + eqs_number[i]]
            += (*RHS)(i);
      }

      //TEST OUTPUT
      // Laplace->Write();
      // Mass->Write();
      // RHS->Write();
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: Coordinates for high order nodes
      Aufgabe:
              Compute the strain couping matrix

      Programming:
      02/2007   WW
   **************************************************************************/
   void CFiniteElementStd::SetHighOrderNodes()
   {
      int i=0;
      setOrder(2);
      // Swap cordinates in case of (x, 0.0, z) only for 2D problem
      if(coordinate_system%10==2)                 // Z has number
      {
         switch(dim)
         {
            case 1:
               for(i=0; i<nNodes; i++)
               {
                  X[i] = MeshElement->nodes[i]->Z();
                  Y[i] = MeshElement->nodes[i]->Y();
                  Z[i] = MeshElement->nodes[i]->X();
               }
               break;
            case 2:
               for(i=0; i<nNodes; i++)
               {
                  X[i] = MeshElement->nodes[i]->X();
                  Y[i] = MeshElement->nodes[i]->Z();
                  Z[i] = MeshElement->nodes[i]->Y();
               }
               break;
            case 3:
               for(i=nnodes; i<nnodesHQ; i++)
               {
                  X[i] = MeshElement->nodes[i]->X();
                  Y[i] = MeshElement->nodes[i]->Y();
                  Z[i] = MeshElement->nodes[i]->Z();
               }

         }
      }
      else
      {
         if(dim==1||dim==2)
         {
            for(i=nnodes; i<nnodesHQ; i++)
            {
               X[i] = MeshElement->nodes[i]->X();
               Y[i] = MeshElement->nodes[i]->Y();
               Z[i] = MeshElement->nodes[i]->Z();
            }
         }
      }
   }
   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: CalcStrainCoupling
      Aufgabe:
              Compute the strain couping matrix

      Programming:
      01/2005   WW
   **************************************************************************/
   void CFiniteElementStd::CalcStrainCoupling()
   {
      int i,k,l,kl, gp, gp_r, gp_s, gp_t;
      double fkt, du=0.0;
      SetHighOrderNodes();
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         for (gp = 0; gp < nGaussPoints; gp++)
         {
            fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

            ComputeGradShapefct(2);
            ComputeShapefct(1);
            ComputeShapefct(2);
            //
            fkt *= CalCoefStrainCouping();
            for(i=0; i<dim; i++ )
            {
               for (k=0;k<nnodes;k++)
               {
                  for (l=0;l<nnodesHQ;l++)
                  {
                     kl = nnodesHQ*i+l;
                     du  = dshapefctHQ[kl];
                     if(i==0&&axisymmetry) du += shapefctHQ[l]/Radius;
                     (*StrainCoupling)(k, kl) += shapefct[k] * du * fkt;
                  }
               }
            }
         }
      }
      setOrder(1);
      // StrainCoupling->Write();
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: Assemby_Gravity
      Aufgabe:
              Assemble the contribution of gravity to RHS in Darcy flow
              to the global system

      Programming:
      01/2005   WW/OK
      08/2006   WW Re-implement
      02/2007   WW Multi-phase flow
   **************************************************************************/
   // Local assembly
   void  CFiniteElementStd::Assemble_Gravity()
   {
      if((coordinate_system)%10!=2)               //NW: exclude (!axisymmetry)
      {
         // 27.2.2007 WW (*GravityMatrix) = 0.0;
         return;
      }
      int i, ii, k;
      // ---- Gauss integral
      int gp_r=0, gp_s=0, gp_t;
      gp_t = 0;
      double fkt, rho;                            //, rich_f;
      double k_rel_iteration;
      // GEO
      //NW  double geo_fac = MediaProp->geo_area;
      if(!FluidProp->CheckGravityCalculation()) return;
      long cshift = 0;                            //WW
      //
      //
      int dof_n = 1;                              // 27.2.2007 WW
      if(PcsType==V || PcsType==P) dof_n = 2;

      //WW 05.01.07
      cshift = 0;
      if(pcs->dof>1)
         cshift = NodeShift[pcs->continuum];

      //rich_f = 1.0;
      //if(PcsType==R) rich_f = -1.0; //WW

      k_rel_iteration = 1.0;

      for (i = 0; i < dof_n*nnodes; i++)
         NodalVal[i] = 0.0;

      // (*GravityMatrix) = 0.0;
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determination
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         //---------------------------------------------------------
         // Compute geometry
         //---------------------------------------------------------
         ComputeGradShapefct(1);                  // Linear interpolation function
         ComputeShapefct(1);                      // Moved from CalCoefLaplace(). 12.3.2007 WW
         // Material
         rho = FluidProp->Density();              //Index,unit,pcs->m_num->ls_theta
         if(gravity_constant<MKleinsteZahl)       // HEAD version
            rho = 1.0;
         else if(HEAD_Flag) rho = 1.0;
         else
            rho *= gravity_constant;
         fkt *= rho;                              //*rich_f;
         //
         for(ii=0; ii<dof_n; ii++)
         {
            if(dof_n==1)
            {
               if(PcsType==T)
                  CalCoefLaplace(false);
               else
                  CalCoefLaplace(true);
            }

            if(dof_n==2)
            {
               if(PcsType==V)
                  CalCoefLaplace2(true, ii*dof_n+1);
               else if(PcsType==P)
                  CalCoefLaplacePSGLOBAL(true, ii*dof_n);
            }
            // Calculate mass matrix
            for (i = 0; i < nnodes; i++)
            {
               for (k = 0; k < dim; k++)
                  NodalVal[i+ii*nnodes] -= fkt*dshapefct[k*nnodes+i]*mat[dim*k+dim-1];
            }
         }
      }
      //
      cshift += NodeShift[problem_dimension_dm];  // 05.01.07 WW
      int ii_sh = 0;
      for(ii=0; ii<dof_n; ii++)                   // 07.02.07 WW
      {
         cshift += NodeShift[ii];
         ii_sh = ii*nnodes;
         for (i=0;i<nnodes;i++)
         {
            eqs_rhs[cshift + eqs_number[i]]
               += k_rel_iteration * NodalVal[i+ii_sh];
            //NW not necessary to multiply geo_area(geo_fac) here. It's already multiplied in ComputeJacobian() through fkt.
            //          eqs_rhs[cshift + eqs_number[i]]
            //                  += k_rel_iteration* geo_fac*NodalVal[i+ii_sh];
            (*RHS)(i+LocalShift+ii_sh) += NodalVal[i+ii_sh];
         }
      }
      //TEST OUTPUT
      //RHS->Write();
   }
   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: Assemby_Gravity
      Aufgabe:
              Assemble the contribution of gravity to RHS in Darcy flow
              to the global system for the pressure equation of multiphase flow

      Programming:
      10/2008   PCH
   **************************************************************************/
   // Local assembly
   void  CFiniteElementStd::Assemble_Gravity_Multiphase()
   {
      if((coordinate_system)%10!=2&&(!axisymmetry))
      {
         // 27.2.2007 WW (*GravityMatrix) = 0.0;
         return;
      }
      int i, ii, k;
      // ---- Gauss integral
      int gp_r=0, gp_s=0, gp_t;
      gp_t = 0;
      double fkt, rho;                            //, rich_f;
      double k_rel_iteration;
      // GEO
      double geo_fac = MediaProp->geo_area;
      if(!FluidProp->CheckGravityCalculation()) return;
      long cshift = 0;                            //WW
      //
      //
      int dof_n = 1;                              // 27.2.2007 WW
      if(PcsType==V) dof_n = 2;

      //WW 05.01.07
      cshift = 0;
      if(pcs->dof>1)
         cshift = NodeShift[pcs->continuum];

      //rich_f = 1.0;
      //if(PcsType==R) rich_f = -1.0; //WW

      k_rel_iteration = 1.0;

      for (i = 0; i < dof_n*nnodes; i++)
         NodalVal[i] = 0.0;

      // (*GravityMatrix) = 0.0;
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determination
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         //---------------------------------------------------------
         // Compute geometry
         //---------------------------------------------------------
         ComputeGradShapefct(1);                  // Linear interpolation function
         ComputeShapefct(1);                      // Moved from CalCoefLaplace(). 12.3.2007 WW
         // Material
         // PCH
         // Pressure equation is the sum of all pressures for all the phases
         // so is gravity. Thus, this sumation should be considered depending on
         // solving the pressure or saturation equation.
         // Thus, the gravity term for the presure equations should cover
         // two phases.
         if(pcs->pcs_type_number==0)              // if the pressure equation,
         {
            int numOfPhases = 2;
            for(int p=0; p< numOfPhases; ++p)
            {
               rho = mfp_vector[p]->Density();
               if(gravity_constant<MKleinsteZahl) // HEAD version
                  rho = 1.0;
               else if(HEAD_Flag) rho = 1.0;
               else
                  rho *= gravity_constant;

                                                  // Initialization
               fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
               fkt *= rho;                        //*rich_f;

               for(ii=0; ii<dof_n; ii++)
               {
                  // CalCoefLaplace does the sumation mat_fac for twophase flow.
                  // This is not right for two phase gravity terms, because the equation
                  // is not summable this way. It should be seperated.
                  if(dof_n==1)
                     CalCoefLaplaceMultiphase(p);
                  if(dof_n==2)
                     CalCoefLaplace2(false, ii*dof_n+1);

                  // Calculate mass matrix
                  for (i = 0; i < nnodes; i++)
                     for (k = 0; k < dim; k++)
                        NodalVal[i+ii*nnodes] -= fkt*dshapefct[k*nnodes+i]*mat[dim*k+dim-1];
               }
            }
         }
         else
         {
            rho = mfp_vector[1]->Density();
            if(gravity_constant<MKleinsteZahl)    // HEAD version
               rho = 1.0;
            else if(HEAD_Flag) rho = 1.0;
            else
               rho *= gravity_constant;
            fkt *= rho;                           //*rich_f;

            for(ii=0; ii<dof_n; ii++)
            {
               // CalCoefLaplace does the sumation mat_fac for twophase flow.
               // This is not right for two phase gravity terms, because the equation
               // is not summable this way. It should be seperated.
               if(dof_n==1)
                  CalCoefLaplace(false);
               if(dof_n==2)
                  CalCoefLaplace2(false, ii*dof_n+1);

               // Calculate mass matrix
               for (i = 0; i < nnodes; i++)
                  for (k = 0; k < dim; k++)
                     NodalVal[i+ii*nnodes] -= fkt*dshapefct[k*nnodes+i]*mat[dim*k+dim-1];
            }
         }
      }

      //
      cshift += NodeShift[problem_dimension_dm];  // 05.01.07 WW
      int ii_sh = 0;
      for(ii=0; ii<dof_n; ii++)                   // 07.02.07 WW
      {
         cshift += NodeShift[ii];
         ii_sh = ii*nnodes;
         for (i=0;i<nnodes;i++)
         {
            eqs_rhs[cshift + eqs_number[i]]
               += k_rel_iteration* geo_fac*NodalVal[i+ii_sh];
            (*RHS)(i+LocalShift+ii_sh) += NodalVal[i+ii_sh];
         }
      }
      //TEST OUTPUT
      //RHS->Write();
   }
   ////////////////////////////////////////////////////////////////
   /*
   void  CFiniteElementStd::Assemble_Gravity()
   {
     if((coordinate_system)%10!=2)
        return;
     int i, j, k, l;
     // ---- Gauss integral
     int gp, gp_r=0, gp_s=0, gp_t;
     gp_t = 0;
     double fkt, rho;
     double k_rel_iteration;
   // GEO
   double geo_fac = MediaProp->geo_area;

   k_rel_iteration = 1.0;

   (*GravityMatrix) = 0.0;
   // Loop over Gauss points
   for (gp = 0; gp < nGaussPoints; gp++)
   {
   //---------------------------------------------------------
   //  Get local coordinates and weights
   //  Compute Jacobian matrix and its determination
   //---------------------------------------------------------
   fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

   //---------------------------------------------------------
   // Compute geometry
   //---------------------------------------------------------
   ComputeGradShapefct(1); // Linear interpolation function

   // Material
   CalCoefLaplace(true);
   rho = FluidProp->Density(Index,unit,pcs->m_num->ls_theta);
   if(gravity_constant<MKleinsteZahl) // HEAD version
   rho = 1.0;
   else if(HEAD_Flag) rho = 1.0;
   else
   rho *= gravity_constant;

   fkt *= rho;
   // Calculate mass matrix
   for (i = 0; i < nnodes; i++)
   for (j = 0; j < nnodes; j++)
   {
   if(j>i) continue;
   for (k = 0; k < dim; k++)
   {
   for(l=0; l<dim; l++)
   (*GravityMatrix)(i,j) += fkt*dshapefct[k*nnodes+i]
   *mat[dim*k+l]* dshapefct[l*nnodes+j];
   }
   }
   }

   //TEST OUTPUT
   //GravityMatrix->Write();

   double* G_coord = NULL;
   if((coordinate_system)/10==1)
   G_coord = X;
   else if((coordinate_system)/10==2)
   G_coord = Y;
   else if((coordinate_system)/10==3)
   G_coord = Z;

   for (i = 0; i < nnodes; i++)
   {
   NodalVal[i] = 0.0;
   for (j = 0; j < nnodes; j++)
   NodalVal[i] -= (*GravityMatrix)(i,j)* G_coord[j];
   }

   for (i=0;i<nnodes;i++)
   {
   pcs->eqs->b[NodeShift[problem_dimension_dm] + eqs_number[i]]
   += k_rel_iteration* geo_fac*NodalVal[i];
   (*RHS)(i+LocalShift) += NodalVal[i];
   }
   //TEST OUTPUT
   //RHS->Write();
   }
   */

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: Velocity calulation

      Programming:  WW
      08/2005
      03/2007   WW  Multi-phase flow
      01/2010   NW  Fix multi-dimensional case
      02/2010   WW  Fix a bug in velocity of the first phase
   **************************************************************************/
   // Local assembly
   void  CFiniteElementStd::Cal_Velocity()
   {
      int i, j, k;
      static double vel[3], vel_g[3];
      double dens_arg[3];                         //AKS
      // ---- Gauss integral
      int gp_r=0, gp_s=0, gp_t;
      double coef = 0.0;
      int dof_n = 1;
      if(PcsType==V || PcsType==P) dof_n = 2;
      //
      gp_t = 0;

      // Get room in the memory for local matrices
      SetMemory();
      // Set material
      SetMaterial();

      ElementValue* gp_ele = ele_gp_value[Index];

      //gp_ele->Velocity = 0.0; // CB commented and inserted below due to conflict with transport calculation, needs velocities
      // Loop over Gauss points
      k = (coordinate_system)%10;
      if(PcsType==T)                              //WW/CB
      {
         if(pcs->pcs_type_number==0)
         {
                                                  // gas pressure
            idx1 = pcs->GetNodeValueIndex("PRESSURE1")+1;
            for(i=0; i<nnodes; i++)
               NodalVal[i] = pcs->GetNodeValue(nodes[i], idx1);
         }
         else if (pcs->pcs_type_number==1)
         {
            idxp21 = pcs->GetNodeValueIndex("PRESSURE_CAP");
                                                  // gas pressure
            idx1 = cpl_pcs->GetNodeValueIndex("PRESSURE1")+1;
            gp_ele = ele_gp_value[Index+(long)pcs->m_msh->ele_vector.size()];
            for(i=0; i<nnodes; i++)
               // P_l = P_g - P_cap
               NodalVal[i] = cpl_pcs->GetNodeValue(nodes[i], idx1) - pcs->GetNodeValue(nodes[i], idxp21);
         }
      }
      else
      {
         // This should be enough for Vw in PS_GLOBAL as well,
         // since the first primary variable is Pw.
         for(i=0; i<nnodes; i++)
         {
            NodalVal[i] = pcs->GetNodeValue(nodes[i], idx1);
            NodalVal1[i] = NodalVal[i];
         }
      }
      //
      if(PcsType==V)
      {
         gp_ele->Velocity_g = 0.0;                //WW
         for(i=0; i<nnodes; i++)
         {
            // 02.2010. WW
            NodalVal2[i] = pcs->GetNodeValue(nodes[i], idxp21);
            NodalVal[i] = NodalVal2[i]  - NodalVal[i];
         }
      }
      if(PcsType==P)
      {
         gp_ele->Velocity_g = 0.0;                //PCH
         // Just get Pnw, which is the secondary variable in PS_GLOBAL
         int idx_pn = pcs->GetNodeValueIndex("PRESSURE2");
         for(i=0; i<nnodes; i++)
         {
            NodalVal1[i] = pcs->GetNodeValue(nodes[i], idx_pn);
         }
      }
      gp_ele->Velocity = 0.0;                     // CB inserted here and commented above due to conflict with transport calculation, needs
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determination
         //---------------------------------------------------------
         GetGaussData(gp, gp_r, gp_s, gp_t);

         //---------------------------------------------------------
         // Compute geometry
         //---------------------------------------------------------
         ComputeGradShapefct(1);                  // Linear interpolation function
         ComputeShapefct(1);                      // Moved from CalCoefLaplace(). 12.3.2007 WW
                                                  //WW/CB
         if((PcsType==T)&&(pcs->pcs_type_number==1))
            flag_cpl_pcs = true;
         // Material
         if(dof_n==1)
            CalCoefLaplace(true);
         else if (dof_n==2 && PcsType==V)         // PCH 05.2009
            CalCoefLaplace2(true,0);
         else if (dof_n==2 && PcsType==P)         // PCH 05.2009
            CalCoefLaplacePSGLOBAL(true,0);
                                                  //WW/CB
         if((PcsType==T)&&(pcs->pcs_type_number==1))
            flag_cpl_pcs = false;
         // Velocity
         for (i = 0; i < dim; i++)
         {
            vel[i] = 0.0;
            for(j=0; j<nnodes; j++)
               vel[i] += NodalVal[j]*dshapefct[i*nnodes+j];
            //			 vel[i] += fabs(NodalVal[j])*dshapefct[i*nnodes+j];
         }
         if(PcsType==V || PcsType==P)             // PCH 05.2009
         {
            for (i = 0; i < dim; i++)
            {
               vel_g[i] = 0.0;
               for(j=0; j<nnodes; j++)
                                                  // Change   NodalVal2 to NodalVal1. 02.2010. WW
                  vel_g[i] += NodalVal2[j]*dshapefct[i*nnodes+j];
            }
         }
         // Gravity term
                                                  //NW
         if(k==2&&(!HEAD_Flag)&&FluidProp->CheckGravityCalculation())
         {
            if((FluidProp->density_model==14))
            {
               dens_arg[0] = interpolate(NodalVal1);
               dens_arg[1] = interpolate(NodalValC)+T_KILVIN_ZERO;
               dens_arg[2] = Index;
               coef  =  gravity_constant*FluidProp->Density(dens_arg);
            }
            else
            {
               coef  =  gravity_constant*FluidProp->Density();
            }
            if(dim==3&&ele_dim==2)
            {
               vel[dim-1] += coef;                //NW local permeability tensor is already transformed to global one in CalCoefLaplace()
               if (PcsType==V || PcsType==P)
               {
                  for(i=0; i<dim; i++)
                  {
                     for(j=0; j<ele_dim; j++)
                     {
                        if(PcsType==V)
                           vel_g[i] += rho_g*gravity_constant*(*MeshElement->tranform_tensor)(i, k)
                              *(*MeshElement->tranform_tensor)(2, k);
                        if(PcsType==P)            // PCH 05.2009
                           vel_g[i] += coef*GasProp->Density()/FluidProp->Density()*(*MeshElement->tranform_tensor)(i, k)
                              *(*MeshElement->tranform_tensor)(2, k);
                     }
                  }
               }
            }                                     // To be correctted
            else
            {
               if(PcsType==V)
               {
                  vel[dim-1] += coef;
                  vel_g[dim-1] += gravity_constant*rho_ga;
               }
               else if(PcsType==P)                // PCH 05.2009
               {
                  vel[dim-1] -= coef;
                  vel_g[dim-1] += gravity_constant*GasProp->Density();
               }
               else
                  vel[dim-1] += coef;
            }
         }
         //
         if(PcsType==V)
         {
            for (i = 0; i < dim; i++)             // 02.2010. WW
            {
               for(j=0; j<dim; j++)
                  gp_ele->Velocity(i, gp) += mat[dim*i+j]*vel[j]/time_unit_factor;
            }
            CalCoefLaplace2(true,3);
            for (i = 0; i < dim; i++)
            {
               for(j=0; j<dim; j++)
                  gp_ele->Velocity_g(i, gp) -= mat[dim*i+j]*vel_g[j]/time_unit_factor;
            }
         }
         else                                     // 02.2010. WW
         {
            for (i = 0; i < dim; i++)
            {
               for(j=0; j<dim; j++)
                  //              gp_ele->Velocity(i, gp) -= mat[dim*i+j]*vel[j];  // unit as that given in input file
                                                  //SI Unit
                  gp_ele->Velocity(i, gp) -= mat[dim*i+j]*vel[j]/time_unit_factor;
            }

         }
         if(PcsType==P)                           // PCH 05.2009
         {
            // Juse use the coefficient of PSGLOBAL Pressure-based velocity (4)
            CalCoefLaplacePSGLOBAL(true,4);
            for (i = 0; i < dim; i++)
            {
               for(j=0; j<dim; j++)
                  gp_ele->Velocity_g(i, gp) -= mat[dim*i+j]*vel_g[j]/time_unit_factor;
            }
         }
         //
      }
      //
      if(pcs->Write_Matrix)
      {
         (*pcs->matrix_file) << "### Element: " << Index << endl;
         (*pcs->matrix_file) << "---Velocity of water " << endl;
         gp_ele->Velocity.Write(*pcs->matrix_file);
         if(gp_ele->Velocity_g.Size()>0)
         {
            (*pcs->matrix_file) << "---Velocity of gas " << endl;
            gp_ele->Velocity_g.Write(*pcs->matrix_file);
         }
      }
      // gp_ele->Velocity.Write();
   }

   /***************************************************************************
      GeoSys - Funktion: Cal_GP_Velocity_FM
      CFiniteElementStd:: Velocity calulation in gauss points from
      node velocities obtained by fluid momentum for one element

      Programming:  SB
      09/2009	first version
   **************************************************************************/
   void  CFiniteElementStd::Cal_GP_Velocity_FM(int *i_ind)
   {
      int i, i_dim;
      static double vel_g_old[3]=
      {
         0.0,0.0,0.0
      }
      , vel_g[3]=
      {
         0.0,0.0,0.0
      };
      // ---- Gauss integral
      int gp_r=0, gp_s=0, gp_t=0;
      double fkt=0.0;                             //OK411 coef = 0.0
      int i_idx;
      // Get fluid_momentum process
      CRFProcess *m_pcs_fm =  PCSGet("FLUID_MOMENTUM");

      ElementValue* gp_ele = ele_gp_value[Index];

      // Gauss point loop
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         // Get gauss point data
         // GetGaussData(gp, gp_r, gp_s, gp_t);
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute the shape function for interpolation within element
         ComputeShapefct(1);

         // Save former gp velocity
         for(i_dim=0;i_dim<dim;i_dim++) vel_g_old[i_dim] = gp_ele->Velocity(i_dim,gp);

         // Interpolate velocity from nodes to gauss point for all three velocity components
         for(i_dim=0;i_dim< dim;i_dim++)
         {
            // Get  velocities from FLUID_MOMENTUM process in element nodes:
            i_idx = i_ind[i_dim];
            for(i=0; i<nnodes; i++)
            {
               NodalVal[i] = m_pcs_fm->GetNodeValue(nodes[i], i_idx);
                                                  //dirty fix for permebility to conductivity
               NodalVal[i] = NodalVal[i] /gravity_constant/1000.0*0.001;
            }
            vel_g[i_dim] = interpolate(NodalVal);
         }                                        // end for dim

         // Set gauss point velocity
         for(i_dim=0; i_dim<dim; i_dim++)
            gp_ele->Velocity(i_dim, gp) = vel_g[i_dim];

         /* 	  // Write out differences:
              if((Index < 100)&&(Index > 0)&&(gp < 3)){
              cout << " Element: " << Index << ", GP: " << gp << ": ";
              cout << "vel_fem: " ;
              for(i_dim=0;i_dim<dim;i_dim++) cout << vel_g_old[i_dim] << "  ";
         //  	  cout << "vel_FM: " ;
         //	  for(i_dim=0;i_dim<dim;i_dim++) cout << vel_g[i_dim] << "  ";
         //	  cout << "vel_diff: " ;
         //	  for(i_dim=0;i_dim<dim;i_dim++) cout << vel_g_old[i_dim]-vel_g[i_dim] << "  ";
              cout << endl;
              }
         */
      }                                           // end gauss point loop

      // Output
      // gp_ele->Velocity.Write();
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: Velocity calulation

      Programming:  WW
      08/2005
      03/2007   WW  Multi-phase flow
      11/2007   CB  this function was only introduced to allow the calculation of
                    the element center of gravity velocity for upwinding

   **************************************************************************/
   // Local assembly
   void  CFiniteElementStd::Cal_Velocity_2()
   {
      int i, j, k;
      static double vel[3], vel_g[3];
      // ---- Gauss integral
      int gp_r=0, gp_s=0, gp_t;
      double coef = 0.0;
      int dof_n = 1;
      if(PcsType==V) dof_n = 2;
      //
      gp_t = 0;

      // Get room in the memory for local matrices
      SetMemory();
      // Set material
      SetMaterial();

      ElementValue* gp_ele = ele_gp_value[Index];

      // Loop over Gauss points
      k = (coordinate_system)%10;
      if(PcsType==T)                              //WW/CB
      {
         if(pcs->pcs_type_number==0)
         {
                                                  // gas pressure
            idx1 = pcs->GetNodeValueIndex("PRESSURE1")+1;
            for(i=0; i<nnodes; i++)
               NodalVal[i] = pcs->GetNodeValue(nodes[i], idx1);
         }
         else if (pcs->pcs_type_number==1)
         {
            idxp21 = pcs->GetNodeValueIndex("PRESSURE_CAP");
                                                  // gas pressure
            idx1 = cpl_pcs->GetNodeValueIndex("PRESSURE1")+1;
            gp_ele = ele_gp_value[Index+(long)pcs->m_msh->ele_vector.size()];
            for(i=0; i<nnodes; i++)
               // P_l = P_g - P_cap
               NodalVal[i] = cpl_pcs->GetNodeValue(nodes[i], idx1) - pcs->GetNodeValue(nodes[i], idxp21);
         }
      }
      else
      {
         for(i=0; i<nnodes; i++)
            NodalVal[i] = pcs->GetNodeValue(nodes[i], idx1);
      }
      //
      if(PcsType==V)
      {
         for(i=0; i<nnodes; i++)
         {
            NodalVal[i] -= pcs->GetNodeValue(nodes[i], idxp21);
            NodalVal1[i] = pcs->GetNodeValue(nodes[i], idxp21);
         }
      }
      //
      gp_ele->Velocity = 0.0;
      //

      gp = 0;

      //for (gp = 0; gp < nGaussPoints; gp++)
      //{
      //---------------------------------------------------------
      //  Get local coordinates and weights
      //  Compute Jacobian matrix and its determination
      //---------------------------------------------------------

      GetGaussData(gp, gp_r, gp_s, gp_t);
      // calculate the velocity at the element center of gravity
      if(PcsType==T) SetCenterGP();               // CB 11/2007

      //---------------------------------------------------------
      // Compute geometry
      //---------------------------------------------------------
      ComputeGradShapefct(1);                     // Linear interpolation function
      ComputeShapefct(1);                         // Moved from CalCoefLaplace(). 12.3.2007 WW
      if((PcsType==T)&&(pcs->pcs_type_number==1)) //WW/CB
         flag_cpl_pcs = true;
      // Material
      if(dof_n==1)
         CalCoefLaplace(true);
      else if (dof_n==2)
         CalCoefLaplace2(true,0);
      if((PcsType==T)&&(pcs->pcs_type_number==1)) //WW/CB
         flag_cpl_pcs = false;

      // Velocity
      for (i = 0; i < dim; i++)
      {
         vel[i] = 0.0;
         for(j=0; j<nnodes; j++)
            vel[i] += NodalVal[j]*dshapefct[i*nnodes+j];
         //			 vel[i] += fabs(NodalVal[j])*dshapefct[i*nnodes+j];
      }
      if(PcsType==V)
      {
         for (i = 0; i < dim; i++)
         {
            vel_g[i] = 0.0;
            for(j=0; j<nnodes; j++)
               vel_g[i] += NodalVal1[j]*dshapefct[i*nnodes+j];
         }
      }
      // Gravity term
      if(k==2&&(!HEAD_Flag))
      {
         coef  =  gravity_constant*FluidProp->Density();
         if(dim==3&&ele_dim==2)
         {
            for(i=0; i<dim; i++)
            {
               for(j=0; j<ele_dim; j++)
               {
                  vel[i] += coef*(*MeshElement->tranform_tensor)(i, k)
                     *(*MeshElement->tranform_tensor)(2, k);
                  if(PcsType==V)
                     vel_g[i] += rho_g*gravity_constant*(*MeshElement->tranform_tensor)(i, k)
                        *(*MeshElement->tranform_tensor)(2, k);
               }
            }
         }                                        // To be correctted
         else
         {
            if(PcsType==V)
            {
               vel[dim-1] -= coef;
               vel_g[dim-1] += gravity_constant*rho_g;
            }
            else
               vel[dim-1] += coef;
         }
      }
      for (i = 0; i < dim; i++)
      {
         for(j=0; j<dim; j++)
            //            gp_ele->Velocity(i, gp) -= mat[dim*i+j]*vel[j];  // unit as that given in input file
            gp_ele->Velocity(i, gp) -= mat[dim*i+j]*vel[j]/time_unit_factor;
      }
      //
      if(PcsType==V)
      {
         CalCoefLaplace2(true,3);
         coef = rhow/rho_ga;
         for (i = 0; i < dim; i++)
         {
            for(j=0; j<dim; j++)
               gp_ele->Velocity_g(i, gp) -= coef*mat[dim*i+j]*vel_g[j]/time_unit_factor;
         }
      }
      //
      //   cout << gp << " " << vel[0] << " " << vel[1] << " " << vel[2] << endl; //Test
      //} // for (gp = 0;...
      //
      if(pcs->Write_Matrix)
      {
         (*pcs->matrix_file) << "### Element: " << Index << endl;
         (*pcs->matrix_file) << "---Velocity of water " << endl;
         gp_ele->Velocity.Write(*pcs->matrix_file);
         if(gp_ele->Velocity_g.Size()>0)
         {
            (*pcs->matrix_file) << "---Velocity of gas " << endl;
            gp_ele->Velocity_g.Write(*pcs->matrix_file);
         }
      }
      // gp_ele->Velocity.Write();
   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: Assemby_Gravity
      Aufgabe:
              Assemble the contribution of known gradient of hydraulic head or
            pressure and gravity to RHS in Darcy flow
              to the global system

      Programming:
      05/2005   PCH
      09/2005   PCH
   **************************************************************************/
   // Local assembly
   void  CFiniteElementStd::AssembleRHS(int dimension)
   {
      // ---- Gauss integral
      int gp_r=0, gp_s=0, gp_t;
      gp_t = 0;
      double fkt, fktG, rho;

      // Declare two known properties on node
      // Since I declare these variables locally, the object of Vec should handle destruction nicely
      // when this local function is done so that I don't bother with memory leak.

      // Initialize Pressure from the value already computed previously.
      CRFProcess* m_pcs = NULL;
      for (size_t i = 0; i < pcs_vector.size(); ++i)
      {
         m_pcs = pcs_vector[i];
         //		if(m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos) // TF
         if (m_pcs->getProcessType() == LIQUID_FLOW)
         {
            PcsType = L;
            break;
            //		} else if (m_pcs->pcs_type_name.find("RICHARDS_FLOW") != string::npos) { // TF
         }
         else if (m_pcs->getProcessType() == RICHARDS_FLOW)
         {
            PcsType = R;
            break;
            //		} else if (m_pcs->pcs_type_name.find("GROUNDWATER_FLOW") // TF
         }
         else if (m_pcs->getProcessType() == GROUNDWATER_FLOW)
         {
            PcsType = G;
            break;
         }
      }
      // Update the process for proper coefficient calculation.
      pcs = m_pcs;
      int nidx1;
      //	if (!(m_pcs->pcs_type_name.find("GROUNDWATER_FLOW") != string::npos)) // TF
      if (!(m_pcs->getProcessType () == GROUNDWATER_FLOW))
         nidx1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1;
      else                                        // then, this is GROUNDWATER_FLOW
      {
         nidx1 = m_pcs->GetNodeValueIndex("HEAD") + 1;
         HEAD_Flag = 1;
         PcsType = G;
      }

      for (int i = 0; i < nnodes; ++i)
      {
         NodalVal[i] = 0.0;
         NodalVal1[i] = m_pcs->GetNodeValue(nodes[i], nidx1);
         NodalVal2[i] = 0.0;
      }

      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determination
         //---------------------------------------------------------
         fktG = fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

         //---------------------------------------------------------
         // Compute geometry
         //---------------------------------------------------------
         ComputeGradShapefct(1);                  // Linear interpolation function
         ComputeShapefct(1);                      // Linear interpolation function

         // Material
         CalCoefLaplace(true);

         // Calculate vector that computes dNj/dx*Ni*Pressure(j)
         // These index are very important.
         rho = FluidProp->Density();
         // Let's get the viscosity too.
         // 09/2010 TF compiler complains about unused value
         //		CFluidProperties *FluidProp = mfp_vector[0];
         if(gravity_constant<MKleinsteZahl)       // HEAD version
            rho = 1.0;
         else if(HEAD_Flag)
         {
            // FS/WW 21.05.2010
             //fkt = fkt*rho * gravity_constant/FluidProp->Viscosity();
            rho = 1.0;
         }
         else
            rho *= gravity_constant;
         //			rho *= gravity_constant/FluidProp->Viscosity();		// This seems to divide viscosity two times. Thus, wrong.

         fktG *= rho;
         for (int i = 0; i < nnodes; i++)
            for (int j = 0; j < nnodes; j++)
               for (int k = 0; k < dim; k++)
               {
                  NodalVal[i]  -= fkt*dshapefct[k*nnodes+j]
                     *mat[dim*dimension+k]* shapefct[i] * NodalVal1[j]; //NW  dshapefct[dimension*nnodes+j] -> dshapefct[k*nnodes+j]
            //	*************************************
            //FS/WW 21.05.2010
                  if (HEAD_Flag)
                      continue;
            //***************************************
                  NodalVal2[i] += fktG*dshapefct[k*nnodes+j]
                     *mat[dim*dimension+k]* shapefct[i] * MeshElement->nodes[j]->Z();  //NW  dshapefct[dimension*nnodes+j] -> dshapefct[k*nnodes+j]
         }
      }

      // Just influence when it's the gravitational direction in the case of Liquid_Flow
      // Thus, it needs one more switch to tell Liquid_Flow and Groundwater_Flow.
      int IsGroundwaterIntheProcesses = 0;
      for (size_t i = 0; i < pcs_vector.size(); ++i)
      {
         m_pcs = pcs_vector[i];
         //		if (m_pcs->pcs_type_name.find("GROUNDWATER_FLOW") != string::npos) // TF
         if (m_pcs->getProcessType () == GROUNDWATER_FLOW)
            IsGroundwaterIntheProcesses = 1;
      }

      // Checking the coordinateflag for proper solution.
      int checkZaxis = 0;
      int coordinateflag = pcs->m_msh->GetCoordinateFlag();
      if( (coordinateflag == 12) || (coordinateflag == 22 && dimension == 1) ||
         (coordinateflag == 32 && dimension == 2) )
         checkZaxis = 1;                          // Then, this gotta be z axis.

      // Compansate the gravity term along Z direction
      if(checkZaxis && IsGroundwaterIntheProcesses == 0 )
         for (int i = 0; i < nnodes; i++)
            NodalVal[i] -= NodalVal2[i];

      // Store the influence into the global vectors.
      m_pcs = PCSGet("FLUID_MOMENTUM");
      for (int i=0;i<nnodes;i++)
#if defined(NEW_EQS)                        //WW
         m_pcs->eqs_new->b[eqs_number[i]] += NodalVal[i];
#else
      m_pcs->eqs->b[eqs_number[i]] += NodalVal[i];
#endif
      // OK. Let's add gravity term that incorporates the density coupling term.
      // This is convenient. The function is already written in RF.
      //Assemble_Gravity();
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Assemble local matrices of parabolic equation to the global system
   Programing:
   01/2005 WW/OK Implementation
   05/2005 WW Dynamics and others -> new equation type
   02/2007 WW Mono-scheme for dual porosity flow
   02/2007 WW Mult-phase flow
   **************************************************************************/
   void CFiniteElementStd::AssembleParabolicEquation()
   {
      int i,j, ii, jj;
      // NUM
      double relax0, relax1;
      //----------------------------------------------------------------------
      long cshift = 0;                            //WW 05.01.07
#if defined(NEW_EQS)
      CSparseMatrix *A = NULL;                    //WW
      if(m_dom)
         A = m_dom->eqs->A;
      else
         A = pcs->eqs_new->A;
#endif

      //WW 05.01.07
      relax0 = pcs->m_num->nls_relaxation;        //WW
      relax1 = 1.0;
      if(relax0<DBL_MIN)
         relax0 = 1.0;
      relax1 = 1.0-relax0;
      //
      cshift = 0;
      if(pcs->dof>1)
         cshift = NodeShift[pcs->continuum];
      //----------------------------------------------------------------------
      // Dynamic
      dynamic = false;
      double *p_n = NULL;
      double fac1, fac2;
      double beta1 = 0.0;
      if(pcs->pcs_type_name_vector.size()&&pcs->pcs_type_name_vector[0].find("DYNAMIC")==0)
      {
         dynamic = true;
         if(pcs->m_num->CheckDynamic())           // why NUM, it is PCS
            beta1  = pcs->m_num->GetDynamicDamping_beta1();
      }
      //----------------------------------------------------------------------
      // Initialize.
      // if (pcs->Memory_Type==2) skip the these initialization
      if(PcsType==V || PcsType==P)                //PCH
         (*Mass2) = 0.0;
      else
         (*Mass) = 0.0;
      (*Laplace) = 0.0;
      //----------------------------------------------------------------------
      // GEO
      // double geo_fac = MediaProp->geo_area;
      //----------------------------------------------------------------------
      // Calculate matrices
      // Mass matrix..........................................................
      if(PcsType==V)                              //WW
      {
         if(pcs->m_num->ele_mass_lumping)
            CalcLumpedMass2();
         else
            CalcMass2();
      }
      else if(PcsType==P)                         //PCH
      {
         if(pcs->m_num->ele_mass_lumping)
            CalcLumpedMassPSGLOBAL();
         else
            CalcMassPSGLOBAL();
      }
      else
      {
         if(pcs->m_num->ele_mass_lumping)
            CalcLumpedMass();
         else
            CalcMass();
      }
      // Laplace matrix.......................................................
      CalcLaplace();
      if(RD_Flag)                                 //YD /WW
         Assemble_DualTransfer();
      if(pcs->Tim->time_control_name.find("NEUMANN")!=string::npos)
         pcs->timebuffer /= mat[0];               //YD
      //======================================================================
      // Assemble global matrix
      //----------------------------------------------------------------------
      // Time discretization
      // ToDo PCS time step
      double dt_inverse = 0.0;
      if(dt<MKleinsteZahl)
      {
         cout<<"\n Zeitschritt ist Null ! Abbruch !"<<endl;
         //abort(); //WW. return;
         return;                                  //exit(NULL); //TK
      }
      else
         dt_inverse = 1.0 / dt;
      //----------------------------------------------------------------------
      // Assemble local left matrix:
      // [C]/dt + theta [K] non_linear_function for static problems
      // [C] + beta1*dt [K] for dynamic problems: ? different equation type
      if(dynamic)
      {
         fac1 = 1.0;
         fac2 = beta1*dt;
      }
      else
      {
         fac1 = dt_inverse;
         fac2 = relax0;                           //unterrelaxation WW theta* non_linear_function_iter; //*geo_fac;
      }

      //Mass matrix
      if(pcs->PartialPS != 1)                     // PCH if not partial-pressure-based
      {
         if(PcsType==V || PcsType==P)             //PCH
            *StiffMatrix    = *Mass2;
         else
            *StiffMatrix    = *Mass;
         (*StiffMatrix) *= fac1;
      }
      // Laplace matrix
      // PCH to reduce PDE to ODE in Saturation model
                                                  // PCH: If equation 2 in Two-phase flow.
      if(pcs->pcs_type_number==1 && pcs->ML_Cap != 0)
      {                                           // then, Laplace equation is no need. Only solve for ODE
      }
      else
      {
         *AuxMatrix      = *Laplace;
      }
      (*AuxMatrix)   *= fac2;
      *StiffMatrix   += *AuxMatrix;
      //----------------------------------------------------------------------
      // Add local matrix to global matrix
      if(PcsType==V || PcsType==P)                // For DOF>1: 03.03.2009 PCH
      {
         int ii_sh, jj_sh;
         long i_sh, j_sh=0;
         for(ii=0;ii<pcs->dof;ii++)
         {
            i_sh = NodeShift[ii];
            ii_sh = ii*nnodes;
            for(jj=0;jj<pcs->dof;jj++)
            {
               j_sh = NodeShift[jj];
               jj_sh = jj*nnodes;
               for(i=0;i<nnodes;i++)
               {
                  for(j=0;j<nnodes;j++)
                  {
#ifdef NEW_EQS
                     (*A)(i_sh+eqs_number[i], j_sh+eqs_number[j]) += \
                        (*StiffMatrix)(i+ii_sh,j+jj_sh);
#else
                     MXInc(i_sh+eqs_number[i], j_sh+eqs_number[j],\
                        (*StiffMatrix)(i+ii_sh,j+jj_sh));
#endif
                  }
               }
            }
         }
      }
      else
      {
                                                  //WW 05.01.07
         cshift += NodeShift[problem_dimension_dm];
         for(i=0;i<nnodes;i++)
         {
            for(j=0;j<nnodes;j++)
            {
#ifdef NEW_EQS
               (*A)(cshift+eqs_number[i], cshift+eqs_number[j]) += \
                  (*StiffMatrix)(i,j);
#else
               MXInc(cshift+eqs_number[i], cshift+eqs_number[j],\
                  (*StiffMatrix)(i,j));
#endif
            }
         }
      }
      //======================================================================
      // Assemble local RHS vector:
      // ( [C]/dt - (1.0-theta) [K] non_linear_function ) u0  for static problems
      // ( [C] + beta1*dt [K] ) dp  for dynamic problems
      if(dynamic)
      {
         fac1 = -1.0;
         fac2 = beta1*dt;
      }
      else
      {
         fac1 = dt_inverse;
         fac2 = relax1;                           // Unerrelaxation. WW  (1.0-theta) * non_linear_function_t0; //*geo_fac;
      }

      // Mass - Storage
      if(pcs->PartialPS != 1)                     // PCH if not partial-pressure-based
      {
         if(PcsType==V || PcsType==P)             //PCH
            *AuxMatrix1 = *Mass2;
         else
            *AuxMatrix1 = *Mass;
         (*AuxMatrix1) *= fac1;
      }
      //Laplace - Diffusion
      // Laplace matrix
                                                  // PCH: If equation 2 in Two-phase flow.
      if(pcs->pcs_type_number==1 && pcs->ML_Cap != 0)
      {                                           // then, Laplace equation is no need. Only solve for ODE
      }
      else
      {
         *AuxMatrix      = *Laplace;
      }
      (*AuxMatrix)  *= fac2;
      *AuxMatrix1   -= *AuxMatrix;
      // 07.01.07 WW
      int idx = idx0;
      if(pcs->continuum==1)
         idx = idxp20;
      for (i=0;i<nnodes; i++)
      {
         NodalVal0[i] = pcs->GetNodeValue(nodes[i],idx);
         NodalVal[i] = 0.0;
      }
      if(PcsType==V)                              // For DOF>1: 27.2.2007 WW
      {
         for (i=0;i<nnodes; i++)
         {
            NodalVal0[i+nnodes] = pcs->GetNodeValue(nodes[i],idxp20);
            NodalVal[i+nnodes] = 0.0;
         }
      }
      else if(PcsType==P)                         // For DOF>1:
      {
         for (i=0;i<nnodes; i++)
         {
            NodalVal0[i+nnodes] = pcs->GetNodeValue(nodes[i],idxSn0);
            NodalVal[i+nnodes] = 0.0;
         }
      }
      AuxMatrix1->multi(NodalVal0, NodalVal);
      // PCH: Type III (Cauchy boundary conditions) in need, it should be added here.

      //
      if(dynamic)
      {
         // Velocity of pressure of the previous step
         p_n = dm_pcs->GetAuxArray();
         for (i=0;i<nnodes; i++)
            NodalVal0[i] = p_n[nodes[i]+NodeShift[problem_dimension_dm]];
         Mass->multi(NodalVal0, NodalVal, -1.0);
         //p_n+vp*dt
         for (i=0;i<nnodes; i++)
         {
            NodalVal0[i] *= dt;
            NodalVal0[i] += pcs->GetNodeValue(nodes[i],idx_pres);
         }
         Laplace->multi(NodalVal0, NodalVal, -1.0);
      }
      //
      if(PcsType==V || PcsType==P)                // For DOF>1: 03.03.2009 PCH
      {
         int ii_sh;
         long i_sh;
         for(ii=0;ii<pcs->dof;ii++)
         {
            i_sh = NodeShift[ii];
            ii_sh = ii*nnodes;
            for (i=0;i<nnodes;i++)
            {
               eqs_rhs[i_sh + eqs_number[i]] += NodalVal[i+ii_sh];
               (*RHS)(i+LocalShift+ii_sh) +=  NodalVal[i+ii_sh];
            }
         }
      }
      else
      {
         for (i=0;i<nnodes;i++)
         {
            eqs_rhs[cshift + eqs_number[i]] += NodalVal[i];
            (*RHS)(i+LocalShift) +=  NodalVal[i];
         }
      }
      //
   }

   //SB4200
   /**************************************************************************
   FEMLib-Method:
   Task: Assemble local matrices of parabolic equation to the global system
   Programing:
   01/2005 WW/OK Implementation
   05/2005 WW Dynamics and others
   09/2005 SB Adapted from AssembleParabolicEquation to assemble transport equation
   01/2010 NW Steady state
   **************************************************************************/
   void CFiniteElementStd::AssembleMixedHyperbolicParabolicEquation()
   {
      int i,j;
      ElementMatrix * EleMat = NULL;              //SB-3
      // NUM
      double theta = pcs->m_num->ls_theta;        //OK
#if defined(NEW_EQS)
      CSparseMatrix *A = NULL;                    //WW
      if(m_dom)
         A = m_dom->eqs->A;
      else
         A = pcs->eqs_new->A;
#endif

      //----------------------------------------------------------------------
      unit[0] = unit[1] = unit[2] = 0.0;
      // Non-linearities
      //  double non_linear_function_iter = 1.0; //OK MediaProp->NonlinearFlowFunction(Index,unit,theta);
      //  double non_linear_function_t0   = 1.0; //OK MediaProp->NonlinearFlowFunction(Index,unit,0.0);
      double fac_mass, fac_laplace, fac_advection, fac_storage, fac_content;
      //if(((aktueller_zeitschritt==1)||(pcs->tim_type_name.compare("TRANSIENT")==0))){   //SB-3
                                                  //SB-3
      if(((aktueller_zeitschritt==1)||(pcs->Memory_Type == 0)))
      {
         // Initialize.
         (*Mass) = 0.0;
         (*Laplace) = 0.0;
         (*Advection) = 0.0;
         (*Storage) = 0.0;
         (*Content) = 0.0;
         //----------------------------------------------------------------------
         // GEO
         // double geo_fac = MediaProp->geo_area;
         //----------------------------------------------------------------------
         // Calculate matrices
         // Mass matrix..........................................................
                                                  //NW
         if(this->pcs->tim_type_name.compare("STEADY")!=0)
         {
            if(pcs->m_num->ele_mass_lumping)
               CalcLumpedMass();
            else
               CalcMass();
         }
         // Laplace matrix.......................................................
         CalcLaplace();
         // Advection matrix.....................................................
         CalcAdvection();
         // Calc Storage Matrix for decay
         CalcStorage();
         // Calc Content Matrix for  saturation changes
         CalcContent();

         // Store matrices to memory for steady state element matrices     //SB-3
         if(pcs->Memory_Type > 0)
         {
            EleMat = pcs->Ele_Matrices[Index];
            EleMat->SetMass_notsym(Mass);
            EleMat->SetLaplace(Laplace);
            EleMat->SetAdvection(Advection);
            EleMat->SetStorage(Storage);
            EleMat->SetContent(Content);
         }

      }                                           //SB-3
      else
      {
         if(Index < 1) cout << "        Skipping calculation of element matrices " << endl;
         // Get Element Matrices
         EleMat = pcs->Ele_Matrices[Index];
         Mass = EleMat->GetMass_notsym();
         Laplace = EleMat->GetLaplace();
         Advection = EleMat->GetAdvection();
         Storage = EleMat->GetStorage();
         Content = EleMat->GetContent();

      }                                           //pcs->tim_type    //SB-3
      //======================================================================
      // Assemble global matrix
      //----------------------------------------------------------------------
      // Time discretization
      // ToDo PCS time step
      double dt_inverse = 0.0;
      if(dt<MKleinsteZahl)
      {
         cout<<"\n Zeitschritt ist Null ! Abbruch !"<<endl;
         return;
      }
      else
         dt_inverse = 1.0 / dt;
      //----------------------------------------------------------------------
      // Assemble local left matrix:
      // [C]/dt + theta [K] non_linear_function for static problems

      fac_mass = dt_inverse;                      //*geo_fac;
      fac_laplace = theta ;                       //* non_linear_function_iter; //*geo_fac;
      fac_advection = theta;
      fac_storage = theta;
      fac_content = theta*dt_inverse;

      //Mass matrix
      *StiffMatrix    = *Mass;
      (*StiffMatrix) *= fac_mass;
      // Laplace matrix
      *AuxMatrix      = *Laplace;
      (*AuxMatrix)   *= fac_laplace;
      *StiffMatrix   += *AuxMatrix;
      // Advection matrix
      *AuxMatrix      = *Advection;
      (*AuxMatrix)   *= fac_advection;
      *StiffMatrix   += *AuxMatrix;
      // Storage matrix
      *AuxMatrix      = *Storage;
      (*AuxMatrix)   *= fac_storage;
      *StiffMatrix   += *AuxMatrix;
      // Content matrix
      *AuxMatrix      = *Content;
      (*AuxMatrix)   *= fac_content;
      *StiffMatrix   += *AuxMatrix;

      //----------------------------------------------------------------------
      // Add local matrix to global matrix
      for(i=0;i<nnodes;i++)
      {
         for(j=0;j<nnodes;j++)
         {
#ifdef NEW_EQS
                                                  //WW
            (*A)(NodeShift[problem_dimension_dm]+eqs_number[i],
               NodeShift[problem_dimension_dm]+eqs_number[j]) += (*StiffMatrix)(i,j);
#else
            MXInc(NodeShift[problem_dimension_dm]+eqs_number[i],\
               NodeShift[problem_dimension_dm]+eqs_number[j],\
               (*StiffMatrix)(i,j));
#endif
         }
      }
      //======================================================================
      // Assemble local RHS vector:
      // ( [C]/dt - (1.0-theta) [K] non_linear_function ) u0  for static problems
      // ( [C] + beta1*dt [K] ) dp  for dynamic problems

      fac_mass = dt_inverse;                      //*geo_fac;
      fac_laplace = -(1.0-theta);                 // * non_linear_function_t0; //*geo_fac;
      fac_advection = -(1.0-theta);
      fac_storage = -(1.0-theta);                 //*lambda
      fac_content = -(1.0-theta)*dt_inverse;

      // Mass - Storage
      *AuxMatrix1    = *Mass;
      (*AuxMatrix1) *= fac_mass;
      //Laplace - Diffusion
      *AuxMatrix     = *Laplace;
      (*AuxMatrix)  *= fac_laplace;
      *AuxMatrix1   += *AuxMatrix;
      // Advection
      *AuxMatrix     = *Advection;
      (*AuxMatrix)  *= fac_advection;
      *AuxMatrix1   += *AuxMatrix;
      // Storage
      *AuxMatrix     = *Storage;
      (*AuxMatrix)  *= fac_storage;
      *AuxMatrix1   += *AuxMatrix;
      // Content
      *AuxMatrix     = *Content;
      (*AuxMatrix)  *= fac_content;
      *AuxMatrix1   += *AuxMatrix;

      for (i=0;i<nnodes; i++)
      {
         NodalVal1[i] = pcs->GetNodeValue(nodes[i],idx0);
         NodalVal[i] = 0.0;
      }
      AuxMatrix1->multi(NodalVal1, NodalVal);     //AuxMatrix1 times vector NodalVal1 = NodalVal
      //----------------------------------------------------------------------
      for (i=0;i<nnodes;i++)
      {
         eqs_rhs[NodeShift[problem_dimension_dm] + eqs_number[i]] += NodalVal[i];
         (*RHS)(i+LocalShift) +=  NodalVal[i];
      }
      //----------------------------------------------------------------------
      //Debug output
      /*
       if(Index < 10){
       cout << " Element Number " << Index << endl;
       cout << " Mass matrix" << endl;
       Mass->Write();
       cout << " Advection matrix" << endl;
       Advection->Write();
       cout << " Dispersion matrix" << endl;
       Laplace->Write();
       cout << " Storage matrix" << endl;
       Storage->Write();
      cout << " Content matrix" << endl;
      Content->Write();
      cout << " Left matrix" << endl;
      StiffMatrix->Write();
      cout << " Right matrix" << endl;
      AuxMatrix1->Write();
      cout << "RHS: " << endl ;
      for (i=0;i<nnodes; i++) cout << "| " << NodalVal[i] << " |" << endl;
      cout << " initial concentrations" << endl;
      for (i=0;i<nnodes; i++) cout << "| " << NodalVal1[i] << " |" << endl;
      //	cout << " RHS vector: " << endl;
      //	for (i=0;i<nnodes; i++) cout << "| " <<  (double)(*RHS)(i+LocalShift) << " |" << endl;
      }
      */

   }
   /**************************************************************************
   FEMLib-Method:
   Task: Assemble local matrices of parabolic equation to the global system
   Comment: Based on hydrosphere, CVFE Method, noch lange nicht allgemein,
   Programing:
   06/2005 MB Implementation
   06/2007 JOD Separation of 1D channel and overland flow
               Introduction of rill depth
            Surface structure with parameter rill_epsilon in st-file
   **************************************************************************/
   void CFiniteElementStd::AssembleParabolicEquationNewton()
   {
      double haaOld[4], haa[4];
      int nidx;
      double axx = 0, ayy = 0, ast=0.0, ckwr[16];
      double swval[4], swold[4];
      double residual[4];
      double **jacobian ;
      double **amat;
      int iups[16];

      jacobian = (double**) Malloc(nnodes * sizeof(double));
      amat = (double**) Malloc(nnodes * sizeof(double));
      for (int i = 0; i < nnodes; i++)
      {
         jacobian[i] = (double*) Malloc(nnodes*sizeof(double));
         amat[i] = (double*) Malloc(nnodes*sizeof(double));
      }

      //////////////////////////// initialize with 0
      MNulleMat(ckwr, nnodes, nnodes);
      MNulleMat(edlluse, nnodes, nnodes);
      MNulleMat(edttuse, nnodes, nnodes);
      for (int i = 0; i < nnodes; i++)
         for (int j = 0; j < nnodes; j++)
      {
         jacobian[i][j] = 0;
         amat[i][j] = 0;
      }

      /////////////////////////// fetch head (depth)
      nidx = pcs->GetNodeValueIndex("HEAD");

      for(int i=0;i<nnodes;i++)
      {
         haa[i] = pcs->GetNodeValue(nodes[i],nidx + 1);
         haaOld[i] = pcs->GetNodeValue(nodes[i],nidx);
      }
      ///////////////////////////// assemble upwinded coefficients
      CalcOverlandCoefficients(haa, &axx, &ayy, &ast);
      // compute axx, ayy, ast  basis functions edlluse, edttuse (element topology (with friction coef and inv. headdiff))
      CalcOverlandNLTERMS(haa, haaOld, swval, swold);
      // compute swval, swold, introduces surface structure in storage term
      CalcOverlandCKWR(haa, ckwr, iups);
      //compute ckwr, iups,  upstream weighting, hydraulic radius for channel
      CalcOverlandUpwindedCoefficients(amat, ckwr, axx, ayy);
      //Form elemental matrix
      /////////////////////////// form residual vector and jacobi matrix
      CalcOverlandResidual(haa, swval, swold, ast, residual, amat);
      AssembleParabolicEquationNewtonJacobian(jacobian, haa, haaOld, axx, ayy, amat, ast, swold, residual, iups);
      /////////////////////////// store
      for(int i = 0; i < nnodes; i++)
      {
#if defined(NEW_EQS)                     //WW
         pcs->eqs_new->b[NodeShift[problem_dimension_dm] + eqs_number[i]] -= residual[i];
#else
         pcs->eqs->b[NodeShift[problem_dimension_dm] + eqs_number[i]] -= residual[i];
#endif
         for(int j=0;j<nnodes;j++)
#if defined(NEW_EQS)                     //WW
            (*pcs->eqs_new->A)(NodeShift[problem_dimension_dm]+eqs_number[i], NodeShift[problem_dimension_dm]+eqs_number[j])
            += jacobian[i][j] ;                   //WW
#else
         MXInc( NodeShift[problem_dimension_dm]+eqs_number[i], NodeShift[problem_dimension_dm]+eqs_number[j], jacobian[i][j] );
#endif
      }

      for(int i = 0; i < nnodes; i++)
      {
         free(jacobian[i]);
         free(amat[i]);
      }
      free(jacobian);
      free(amat);

   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculates jacobi matrix for AssembleParabolicEquationNewton()
         be carefull with epsilon
   Programing:
   06/2007 JOD Implementation
   **************************************************************************/
   void CFiniteElementStd::AssembleParabolicEquationNewtonJacobian(double** jacob, double* haa, double* hOld, double axx, double ayy, double** amat, double ast, double* swold, double* residual, int* iups)
   {

      // double** jacob;
      double hEps[4], hKeep[4], swval_eps[4];
      double sumjac, stor_eps, akrw, remember;
      double epsilon = 1.e-7;                     // be carefull, like in primary variable dependent source terms (critical depth, normal depth)

      /* jacob = (double**) Malloc(nnodes * sizeof(double));
       for (int i = 0; i < nnodes; i++)
        jacob[i] = (double*) Malloc(nnodes*sizeof(double));
       for (int i = 0; i < nnodes; i++)
         for (int j = 0; j < nnodes; j++)
           jacob[i][j]= 0.0;
      */
      for(int i = 0; i < nnodes; i++)
      {
         hEps[i] = haa[i] + epsilon;
         hKeep[i] = haa[i];
      }

      CalcOverlandNLTERMS(hEps, hOld, swval_eps, swold);
      // compute swval_eps, swold, introduces surface structure in storage term

      for(int i = 0; i < nnodes; i++)             // Form jacobian !
      {
         remember = haa[i];
         haa[i] = hEps[i];
         sumjac = 0.0;

         for(int j = 0; j < nnodes; j++)
         {
            if(i != j)                            // nondiagonal
            {
               CalcOverlandCKWRatNodes(i, j, haa, &akrw, iups);
               //compute ckwr, iups,  upstream weighting, hydraulic radius for channel
               jacob[j][i] = CalcOverlandJacobiNodes(i, j, haa, hKeep, akrw, axx, ayy, amat, &sumjac) / epsilon;
               // if(MediaProp->channel ==1)
               // sumjac +=  swval_eps[i] * ast * (Haa[i] - Hold[i]);
            }                                     //end if (i!=j)
         }                                        //end j

         //Compute diagonal for row i, Lump the storage term
         stor_eps = ast * (swval_eps[i] - swold[i]);
         sumjac = sumjac + stor_eps;
         jacob[i][i] = (sumjac - residual[i]) / epsilon ;
         haa[i] = remember;
      }                                           // end i

      // return jacob;

   }

   /***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: Assemby_strainCPL
      Aufgabe:
              Assemble local metrices of strain coupling
              to the global system

      Programming:
      01/2005   WW/OK
      05/2005   WW dyn
      07/2005   WW Change due to geometry element object
   **************************************************************************/
   void CFiniteElementStd::Assemble_strainCPL()
   {
      int i, j;
      double *u_n = NULL;                         // Dynamic
      double fac;
      int Residual = -1;
#if defined(NEW_EQS)
      CSparseMatrix *A = NULL;
      if(m_dom)
         A = m_dom->eqsH->A;
      else
         A = pcs->eqs_new->A;
#endif

      fac = 1.0 / dt;
      if(D_Flag != 41)
         Residual = 0;
      else                                        // Mono
      {
         if(pcs_deformation>100)                  // Pls
            Residual = 1;
      }
      if(dynamic)
      {
         Residual = 2;
         fac = pcs->m_num->GetDynamicDamping_beta1()*dt;
         u_n = dm_pcs->GetAuxArray();
      }
      if(MediaProp->storage_model==7)             //RW/WW
         fac *= MediaProp->storage_model_values[0];
      //
      for (i=nnodes;i<nnodesHQ;i++)
         nodes[i] = MeshElement->nodes_index[i];
      (*StrainCoupling) = 0.0;
      CalcStrainCoupling();
      //	if(D_Flag != 41&&aktueller_zeitschritt>1)
      if(Residual>=0)
      {                                           // Incorparate this after the first time step
         if(Residual==0)                          // Partitioned
         {
            for (i=0;i<nnodesHQ;i++)
            {
               NodalVal2[i] = -fac*(dm_pcs->GetNodeValue(nodes[i],Idx_dm1[0])-dm_pcs->GetNodeValue(nodes[i],Idx_dm0[0]));
               NodalVal3[i] = -fac*(dm_pcs->GetNodeValue(nodes[i],Idx_dm1[1])-dm_pcs->GetNodeValue(nodes[i],Idx_dm0[1]));
               if(dim==3)                         // 3D.
                  NodalVal4[i] = -fac*(dm_pcs->GetNodeValue(nodes[i],Idx_dm1[2])
                     -dm_pcs->GetNodeValue(nodes[i],Idx_dm0[2]));
            }
         }
         else if(Residual==1)                     //Mono
         {
            // du is stored in u_0
            for (i=0;i<nnodesHQ;i++)
            {
               NodalVal2[i] = -fac*pcs->GetNodeValue(nodes[i],Idx_dm0[0]);
               NodalVal3[i] = -fac*pcs->GetNodeValue(nodes[i],Idx_dm0[1]);
               if(dim==3)                         // 3D.
                  NodalVal4[i] = -fac*pcs->GetNodeValue(nodes[i],Idx_dm0[2]);
            }
         }
         else if(Residual==2)                     //Mono dynamic
         {
            // da is stored in a_0
            // v_{n+1} = v_{n}+a_n*dt+beta1*dt*da
            // a_n is in dm_pcs->ARRAY
            for (i=0;i<nnodesHQ;i++)
            {
               NodalVal2[i] = -(pcs->GetNodeValue(nodes[i],idx_vel_disp[0])
                  +fac*pcs->GetNodeValue(nodes[i],Idx_dm0[0])
                  +u_n[nodes[i]]*dt);
               NodalVal3[i] = -(pcs->GetNodeValue(nodes[i],idx_vel_disp[1])
                  +fac*pcs->GetNodeValue(nodes[i],Idx_dm0[1])
                  +u_n[nodes[i]+NodeShift[1]]*dt);
               if(dim==3)                         // 3D.
                  NodalVal4[i] = -(pcs->GetNodeValue(nodes[i],idx_vel_disp[2])
                     +fac*pcs->GetNodeValue(nodes[i],Idx_dm0[2])
                     +u_n[nodes[i]+NodeShift[2]]*dt);

            }
         }

         for (i=0;i<nnodes; i++)
         {
            NodalVal[i] = 0.0;
            for (j=0;j<nnodesHQ; j++)
            {
               NodalVal[i] += (*StrainCoupling)(i,j)*NodalVal2[j];
               NodalVal[i] += (*StrainCoupling)(i,j+nnodesHQ)*NodalVal3[j];
               if(dim==3)                         // 3D.
                  NodalVal[i] += (*StrainCoupling)(i,j+2*nnodesHQ)*NodalVal4[j];
            }
         }
         // Add RHS
         for (i=0;i<nnodes;i++)
         {
            eqs_rhs[NodeShift[problem_dimension_dm] + eqs_number[i]]
               += NodalVal[i];
            (*RHS)(i+LocalShift) +=  NodalVal[i];
         }
      }
      // Monolithic scheme.
      if(D_Flag == 41)
      {
         // if Richard, StrainCoupling should be multiplied with -1.
         for(i=0;i<nnodes;i++)
         {
            for(j=0;j<nnodesHQ;j++)
            {
#ifdef NEW_EQS
               (*A)(NodeShift[problem_dimension_dm] + eqs_number[i],
                  eqs_number[j]+NodeShift[0]) += (*StrainCoupling)(i,j)*fac;
               (*A)(NodeShift[problem_dimension_dm] + eqs_number[i],
                  eqs_number[j]+NodeShift[1]) += (*StrainCoupling)(i,j+nnodesHQ)*fac;
               if(problem_dimension_dm==3)
                  (*A)(NodeShift[problem_dimension_dm] + eqs_number[i],
                  eqs_number[j]+NodeShift[2]) += (*StrainCoupling)(i,j+2*nnodesHQ)*fac;
#else
               MXInc(NodeShift[problem_dimension_dm] + eqs_number[i],
                  eqs_number[j]+NodeShift[0],(*StrainCoupling)(i,j)*fac);
               MXInc(NodeShift[problem_dimension_dm] + eqs_number[i],
                  eqs_number[j]+NodeShift[1], (*StrainCoupling)(i,j+nnodesHQ)*fac);
               if(problem_dimension_dm==3)
                  MXInc(NodeShift[problem_dimension_dm] + eqs_number[i],
                     eqs_number[j]+NodeShift[2], (*StrainCoupling)(i,j+2*nnodesHQ)*fac);
#endif
            }
         }
      }
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Assemble local mass matrices to the global system
   Programing:
   05/2005 PCH Implementation
   **************************************************************************/
   void CFiniteElementStd::AssembleMassMatrix(int option)
   {
      // Calculate matrices
      // Mass matrix..........................................................
      // ---- Gauss integral
      int gp;
      int gp_r=0,gp_s=0,gp_t=0;
      double fkt,mat_fac;
      // Material
      mat_fac = 1.0;

#if defined(NEW_EQS)
      CSparseMatrix *A = NULL;                    //PCH
      if(m_dom)
         A = m_dom->eqs->A;
      else
         A = pcs->eqs_new->A;
#endif
      //----------------------------------------------------------------------
      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

         // Compute geometry
         ComputeShapefct(1);                      // Linear interpolation function

         if(option == 0)                          // The consistent method
         {
            // Calculate mass matrix
            for (int i = 0; i < nnodes; i++)
               for (int j = 0; j < nnodes; j++)
            {
               //			if(j>i) continue;
               (*Mass)(i,j) += fkt *shapefct[i]*shapefct[j];
            }
         }
         else if(option == 1)                     // The lumped method
         {
            // Calculate mass matrix
            for (int i = 0; i < nnodes; i++)
               for (int j = 0; j < nnodes; j++)
            {
               (*Mass)(i,i) += fkt *shapefct[i]*shapefct[j];
            }
         }
      }

      //----------------------------------------------------------------------
      // Add local matrix to global matrix
      if(PcsType==V || PcsType==P)                // For DOF>1: 03.03.2009 PCH
      {
         int ii_sh, jj_sh;
         long i_sh, j_sh=0;
         for(int ii=0;ii<pcs->dof;ii++)
         {
            i_sh = NodeShift[ii];
            ii_sh = ii*nnodes;
            for(int jj=0;jj<pcs->dof;jj++)
            {
               j_sh = NodeShift[jj];
               jj_sh = jj*nnodes;
               for(int i=0;i<nnodes;i++)
               {
                  for(int j=0;j<nnodes;j++)
                  {
#ifdef NEW_EQS
                     (*A)(i_sh+eqs_number[i], j_sh+eqs_number[j]) += \
                        (*Mass)(i+ii_sh,j+jj_sh);
#else
                     MXInc(i_sh+eqs_number[i], j_sh+eqs_number[j],\
                        (*Mass)(i+ii_sh,j+jj_sh));
#endif
                  }
               }
            }
         }
      }
      else
      {
         int cshift = 0;
                                                  //WW 05.01.07
         cshift += NodeShift[problem_dimension_dm];
         for(int i=0;i<nnodes;i++)
         {
            for(int j=0;j<nnodes;j++)
            {
#ifdef NEW_EQS
               (*A)(cshift+eqs_number[i], cshift+eqs_number[j]) += \
                  (*Mass)(i,j);
#else
               MXInc(cshift+eqs_number[i], cshift+eqs_number[j],\
                  (*Mass)(i,j));
#endif
            }
         }
      }
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Config material and knowns for local assembly
   Programing:
   08/2008 WW Implementation
   **************************************************************************/
   void CFiniteElementStd::Config()
   {
      int i, nn;
      //----------------------------------------------------------------------
      //OK index = m_dom->elements[e]->global_number;
      index = Index;
      //----------------------------------------------------------------------
      nn = nnodes;
                                                  // ?2WW
      if(pcs->type==41||pcs->type==4) nn = nnodesHQ;
      //----------------------------------------------------------------------
      // For DDC WW
#ifdef NEW_EQS
      eqs_rhs = pcs->eqs_new->b;
#else
      eqs_rhs = pcs->eqs->b;
#endif
      // EQS indices
      if(m_dom)                                   //WW
      {
         eqs_rhs = m_dom->eqs->b;
         for(i=0;i<nn;i++)
            eqs_number[i] = element_nodes_dom[i]; //WW
         if(pcs->dof>1)                           //12.12.2007 WW
         {
            for(i=0; i<pcs->dof; i++)
               NodeShift[i]=i*m_dom->nnodes_dom;
         }
      }
      else
      {                                           //OK4111
         for(i=0;i<nn;i++)
            eqs_number[i] = MeshElement->nodes[i]->GetEquationIndex();
      }
      //----------------------------------------------------------------------
      // Get room in the memory for local matrices
      SetMemory();
      //----------------------------------------------------------------------
      // Set material
      SetMaterial();
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
                                                  // ?2WW
      if((D_Flag==41&&pcs_deformation>100)||dynamic)
         dm_pcs = (process::CRFProcessDeformation*)pcs;
      //----------------------------------------------------------------------
      // Initialize RHS
      if(pcs->Memory_Type>0)
      {
         for(i=LocalShift;(size_t)i<RHS->Size();i++)
            (*RHS)(i) = 0.0;
      }
      else
         (*RHS) = 0.0;
      //----------------------------------------------------------------------
      // Node value of the previous time step
      int idx00 = idx0;                           //----------WW 05.01.07
      int idx11 = idx1;
      if(pcs->GetContinnumType()==1)
      {
         idx00 = idxp20;
         idx11 = idxp21;
      }
      for(i=0;i<nnodes;i++)
      {
         NodalVal0[i] = pcs->GetNodeValue(nodes[i],idx00);
         NodalVal1[i] = pcs->GetNodeValue(nodes[i],idx11);
      }                                           //----------WW 05.01.07
      if(PcsType==V)                              // 25.2.2007
      {
         for(i=0;i<nnodes;i++)
            NodalVal_p2[i] = pcs->GetNodeValue(nodes[i],idxp21);
      }
      if(PcsType==P)
      {
         for(i=0;i<nnodes;i++)
            NodalVal_SatNW[i] = pcs->GetNodeValue(nodes[i],idxSn1);
      }                                           //----------WW 05.01.07
      if(cpl_pcs)                                 // ?2WW: flags are necessary
      {
         for(i=0;i<nnodes;i++)
         {
            NodalValC[i] = cpl_pcs->GetNodeValue(nodes[i],idx_c0);
            NodalValC1[i] = cpl_pcs->GetNodeValue(nodes[i],idx_c1);
            if(cpl_pcs->type==1212)
               NodalVal_p2[i] = cpl_pcs->GetNodeValue(nodes[i],idx_c1+2);
            NodalVal_p20[i] = cpl_pcs->GetNodeValue(nodes[i],idx_c0+2);
         }
      }
   }
   /**************************************************************************
   FEMLib-Method:
   Task: Assemble local matrices to the global system
   Programing:
   01/2005 WW Implementation
   02/2005 OK Richards flow
   02/2005 WW Matrix output
   03/2005 WW Heat transport
   04/2005 OK MSH
   05/2005 OK regional PCS
   08/2005 OK Air (gas) flow
   10/2005 OK DDC
   06/2005 WW Adjustment in DDC
   07/2007 WW Nonisothermal multi-phase flow
   10/2007 OK Two-phase flow
   08/2008 WW Extract the configuration of material properties and knowns as
   a single inline function
   **************************************************************************/
   void CFiniteElementStd::Assembly()
   {
      int i;
      Config();                                   //26.08.2008
      //======================================================================
      switch(PcsType)
      {
         //....................................................................
         case L:                                  // Liquid flow
            AssembleParabolicEquation();
            Assemble_Gravity();
            if(dm_pcs)
               Assemble_strainCPL();
            break;
            //....................................................................
            //case U: // Unconfined flow  //  part of Groundwater flow mmp keyword ($UNCONFINED)
            //....................................................................
         case G:                                  // Groundwater flow
            AssembleParabolicEquation();
            //RHS->Write();
            if(dm_pcs)
               Assemble_strainCPL();
            break;
            //....................................................................
         case T:                                  // Two-phase flow
            if(pcs->pcs_type_number==0)
            {
               // Start partial-pressure-based model
               pcs->PartialPS = 0;

               AssembleParabolicEquation();

               if(pcs->PartialPS == 1)            // If it is partial-pressure-based
                  AssembleRHSVector();
               //			PrintTheSetOfElementMatrices("Pressure1");
               AssembleCapillaryEffect();
               Assemble_Gravity_Multiphase();
            }
            else if(pcs->pcs_type_number==1)
            {
               // Turn off the partial-pressure-based model for Snw equation
               pcs->PartialPS = 0;

               pcs->ML_Cap = 0;
               AssembleParabolicEquation();
               pcs->ML_Cap = 0;

               AssembleRHSVector();
               Assemble_Gravity_Multiphase();
            }
            break;
            //....................................................................
         case C:                                  // Componental flow
            for(i=0;i<nnodes;i++)
               NodalVal_Sat[i] = pcs->GetNodeValue(nodes[i], idxS);
            break;
            //....................................................................
         case H:                                  // Heat transport
            heat_phase_change = false;            // ?2WW
            //  if(SolidProp->GetCapacityModel()==2) // Boiling model
            //    CalNodalEnthalpy();
                                                  //CMCD4213
            AssembleMixedHyperbolicParabolicEquation();
            if(FluidProp->density_model==14 && MediaProp->heat_diffusion_model==273 && cpl_pcs )
               Assemble_RHS_HEAT_TRANSPORT();     // This include when need pressure terms n dp/dt + nv.Nabla p//AKS
            if(MediaProp->evaporation==647)
               Assemble_RHS_HEAT_TRANSPORT2();    //AKS
            break;
            //....................................................................
         case M:                                  // Mass transport
                                                  //SB4200
            AssembleMixedHyperbolicParabolicEquation();
            break;
            //....................................................................
         case O:                                  // Overland flow
            if(pcs->m_num->nls_method == 0)       //PICARD
               AssembleParabolicEquation();       //OK
            else
               AssembleParabolicEquationNewton(); //NEWTON
            break;
            //....................................................................
         case R:                                  // Richards flow
            if(MediaProp->heat_diffusion_model==273)
               CalcRHS_by_ThermalDiffusion();
            AssembleParabolicEquation();          //OK
            Assemble_Gravity();
            if(dm_pcs)
               Assemble_strainCPL();
            break;
            //....................................................................
         case F:                                  // Fluid Momentum - Assembly handled in Assembly in Fluid_Momentum file
            break;
            //....................................................................
         case A:                                  // Air (gas) flow
                                                  //To account advection like term nv.Nabla p
            AssembleMixedHyperbolicParabolicEquation();
                                                  //AKS
            if(MediaProp->heat_diffusion_model==273 && cpl_pcs )
               Assemble_RHS_AIR_FLOW();           // n*drho/dt + Nabla.[rho*k/mu rho g]//AKS
            break;
         case V:                                  // Multi-phase flow 24.02.2007 WW
            AssembleParabolicEquation();
            Assemble_Gravity();
            if(cpl_pcs&&MediaProp->heat_diffusion_model==273)
               Assemble_RHS_T_MPhaseFlow();
            if(dm_pcs)
               Assemble_RHS_M();
            break;

         case P:                                  // PS_GLOBAL for Multi-phase flow 03.03 2009 PCH
            AssembleParabolicEquation();
            PrintTheSetOfElementMatrices("Laplace");
            if(pcs->num_type_name.find("DirectPc")!=string::npos)
            {
               Assemble_RHS_Pc();
               PrintTheSetOfElementMatrices("RHS_Pc");
            }
            Assemble_Gravity();

            if(dm_pcs)
               Assemble_RHS_M();
            Assemble_RHS_T_PSGlobal();
            break;

            //....................................................................
         default:
            cout << "Fatal error: No valid PCS type" << endl;
            break;
      }

      //----------------------------------------------------------------------
      // Irregulaere Knoten eliminieren
      //----------------------------------------------------------------------
      // Output matrices
      if(pcs->Write_Matrix)
      {
         (*pcs->matrix_file) << "### Element: " << Index << endl;
         (*pcs->matrix_file) << "---Mass matrix: " << endl;
         if(Mass)
            Mass->Write(*pcs->matrix_file);
         else if(Mass2)
            Mass2->Write(*pcs->matrix_file);
         (*pcs->matrix_file) << "---Laplacian matrix: " << endl;
         Laplace->Write(*pcs->matrix_file);
         if(Advection)
         {
                                                  //CMCD
            (*pcs->matrix_file) << "---Advective matrix: " << endl;
            Advection->Write(*pcs->matrix_file);
         }
         if(StrainCoupling)
         {
            (*pcs->matrix_file) << "---Strain couping matrix: " << endl;
            StrainCoupling->Write(*pcs->matrix_file);
         }
         (*pcs->matrix_file) << "---RHS: " <<endl;
         RHS->Write(*pcs->matrix_file);
         (*pcs->matrix_file) <<endl;
         (*pcs->matrix_file) << "Stiffness: " <<endl;
         StiffMatrix->Write(*pcs->matrix_file);
         (*pcs->matrix_file) <<endl;
      }
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Assemble local matrices to the global system
   Programing:
   01/2005 WW Implementation
   02/2005 OK Richards flow
   02/2005 WW Matrix output
   03/2005 WW Heat transport
   08/2005 PCH for Fluid_Momentum
   last modification:
   **************************************************************************/
   void  CFiniteElementStd::Assembly(int option, int dimension)
   {
      int i,nn;
      //----------------------------------------------------------------------
#ifdef PARALLEL
      index = m_dom->elements[e]->global_number;
#else
      index = Index;
#endif
      //----------------------------------------------------------------------

      nn = nnodes;
      // PCH should check the following line carefully.
      if(pcs->type==41||pcs->type==4) nn = nnodesHQ;

#ifdef NEW_EQS                              //PCH
      eqs_rhs = pcs->eqs_new->b;
#else
      eqs_rhs = pcs->eqs->b;
#endif

      for(i=0;i<nn;i++)
      {
#ifdef PARALLEL
         eqs_number[i] = MeshElement->domain_nodes[i];
#else
         eqs_number[i] = MeshElement->nodes[i]->GetEquationIndex();
#endif
      }

      // Get room in the memory for local matrices
      SetMemory();

      // Set material
      SetMaterial();

      // Initialize.
      // if (pcs->Memory_Type==2) skip the these initialization
      (*Mass) = 0.0;
      (*Laplace) = 0.0;
      if(pcs->Memory_Type>0)
      {
         for(i=LocalShift;(size_t)i<RHS->Size();i++)
            (*RHS)(i) = 0.0;
      }
      else
         (*RHS) = 0.0;

      // Fluid Momentum
      AssembleMassMatrix(option);                 // This is exactly same with CalcMass().
      AssembleRHS(dimension);
      //Output matrices
      if(pcs->Write_Matrix)
      {
         for (i = 0; i < nnodes; i++)
            (*RHS)(i) = NodalVal[i];
         (*pcs->matrix_file) << "### Element: " << Index << endl;
         (*pcs->matrix_file) << "---Mass matrix: " << endl;
         Mass->Write(*pcs->matrix_file);
         (*pcs->matrix_file) << "---Laplacian matrix: " << endl;
         Laplace->Write(*pcs->matrix_file);
         (*pcs->matrix_file) << "---RHS: " <<endl;
         RHS->Write(*pcs->matrix_file);
         (*pcs->matrix_file) <<endl;
         (*pcs->matrix_file) << "Stiffness: " <<endl;
         StiffMatrix->Write(*pcs->matrix_file);
         (*pcs->matrix_file) <<endl;
      }

   }
   /**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   18/02/2006 WW Implementation
   **************************************************************************/
   void CFiniteElementStd::ExtropolateGauss(CRFProcess *m_pcs, const int idof)
   {
      int i, j, gp, gp_r, gp_s, gp_t, idx_v2=0;
      int i_s, i_e, ish;
      double EV, EV1=0.0, varx=0.0;
      double r=0.0;
      //
      MshElemType::type ElementType = MeshElement->GetElementType();

      if(m_pcs->type==1212 || m_pcs->type==1313)  // Multi-phase flow 03.2009 PCH
      {
         switch(idof)
         {
            case 0:
               idx_v2 = m_pcs->GetNodeValueIndex("VELOCITY_X2");
               break;
            case 1:
               idx_v2 = m_pcs->GetNodeValueIndex("VELOCITY_Y2");
               break;
            case 2:
               idx_v2 = m_pcs->GetNodeValueIndex("VELOCITY_Z2");
               break;
         }

      }
      // For strain and stress extropolation all element types
      // Number of elements associated to nodes
      for(i=0; i<nnodes; i++)
         dbuff[i] = (double)MeshElement->nodes[i]->connected_elements.size();
      //
      gp_r=gp_s=gp_t=gp=0;
      ElementValue* gp_ele = ele_gp_value[Index];
      //
      for(gp=0; gp<nGaussPoints; gp++)
      {
         if(ElementType==MshElemType::QUAD||ElementType==MshElemType::HEXAHEDRON)
         {
            if(ElementType==MshElemType::QUAD)
            {
               gp_r = (int)(gp/nGauss);
               gp_s = gp%nGauss;
               gp_t = 0;
            }
            else if(ElementType==MshElemType::HEXAHEDRON)
            {
               gp_r = (int)(gp/(nGauss*nGauss));
               gp_s = (gp%(nGauss*nGauss));
               gp_t = gp_s%nGauss;
               gp_s /= nGauss;
            }
            i = GetLocalIndex(gp_r, gp_s, gp_t);
            if(i==-1) continue;
         }
         else
            i = gp;
         NodalVal1[i] = gp_ele->Velocity(idof,gp)*time_unit_factor;
         //
         //
                                                  // PCH 05.2009
         if(m_pcs->type==1212 || m_pcs->type==1313)
            NodalVal2[i] =gp_ele->Velocity_g(idof,gp)*time_unit_factor;
      }

      if(ElementType==MshElemType::QUAD||ElementType==MshElemType::HEXAHEDRON)
      {
         Xi_p = 0.0;
         for (gp = 0; gp < nGauss; gp++)
         {
            r = MXPGaussPkt(nGauss, gp);
            if(fabs(r)>Xi_p) Xi_p = fabs(r);
         }
         r = 1.0/Xi_p;
         Xi_p = r;
      }
      //
      i_s=0;
      i_e=nnodes;
      ish=0;
      if(ElementType==MshElemType::TETRAHEDRON)   // tet
      {
         i_s=1;
         i_e=nnodes+1;
         ish=1;
      }
      //---------------------------------------------------------
      // Mapping Gauss point strains to nodes and update nodes
      // strains:
      //---------------------------------------------------------
      for(i=0; i<nnodes; i++)
      {
         EV = EV1 = varx = 0.0;
         SetExtropoGaussPoints(i);
         //
         ComputeShapefct(1);                      // Linear interpolation function
         for(j=i_s; j<i_e; j++)
            EV += NodalVal1[j]*shapefct[j-ish];
         // Average value of the contribution of ell neighbor elements
         EV /= dbuff[i];
         EV += m_pcs->GetNodeValue(nodes[i],idx_vel[idof]);
         m_pcs->SetNodeValue (nodes[i], idx_vel[idof], EV);
         //
                                                  // Multi-phase flow PCH 05.2009
         if(m_pcs->type==1212 || m_pcs->type==1313)
         {
            for(j=i_s; j<i_e; j++)
               EV1 += NodalVal2[j]*shapefct[j-ish];
            //
            EV1 /= dbuff[i];
            EV1 += m_pcs->GetNodeValue(nodes[i],idx_v2);
            m_pcs->SetNodeValue (nodes[i], idx_v2, EV1);
         }
         //
      }
   }

   /***********************************************************************
    27.03.2007 WW
   ***********************************************************************/
   void CFiniteElementStd::CalcSatution()
   {
      int i, j, gp, gp_r, gp_s, gp_t, idx_cp, idx_S;
      int i_s, i_e, ish;
      //  int l1,l2,l3,l4; //, counter;
      double sign, eS=0.0;
      double r=0.0;
      //
      MshElemType::type ElementType = MeshElement->GetElementType();
      //----------------------------------------------------------------------
      // Media
      int mmp_index=0;
      long group = MeshElement->GetPatchIndex();
      mmp_index = group;
      //
      if(pcs->type==22)
      {
         if(pcs->GetContinnumType()== 0)          // Matrix //WW
            mmp_index = 2*group;
         else                                     // fracture //WW
            mmp_index = 2*group+1;
      }
      MediaProp = mmp_vector[mmp_index];
      MediaProp->m_pcs = pcs;
      MediaProp->Fem_Ele_Std = this;
      //
      sign = -1.0;
      idx_cp = pcs->GetNodeValueIndex("PRESSURE1")+1;
      idx_S =  pcs->GetNodeValueIndex("SATURATION1")+1;
                                                  // Dual Richards
      if(pcs->type==22&&pcs->GetContinnumType()==1)
      {
         idx_cp = pcs->GetNodeValueIndex("PRESSURE2")+1;
         idx_S =  pcs->GetNodeValueIndex("SATURATION2")+1;
      }
      if(pcs->type==1212)
         sign = 1.0;
      //
      for(i=0; i<nnodes; i++)
      {
         // Number of elements associated to nodes
         dbuff[i] = (double)MeshElement->nodes[i]->connected_elements.size();
                                                  // pressure
         NodalVal0[i] = sign*pcs->GetNodeValue(nodes[i],idx_cp);
      }
      //
      gp_r=gp_s=gp_t=gp=0;
      //
      for(gp=0; gp<nGaussPoints; gp++)
      {
         SetGaussPoint(gp, gp_r, gp_s, gp_t);
         if(ElementType==2||ElementType==3)
         {
            i = GetLocalIndex(gp_r, gp_s, gp_t);
            if(i==-1) continue;
         }
         else
            i = gp;
         //
         if(i>nnodes) continue;
         ComputeShapefct(1);
         //
         PG = interpolate(NodalVal0);
         NodalVal_Sat[i] = MediaProp->SaturationCapillaryPressureFunction(PG,0);
      }

      if(ElementType==MshElemType::QUAD || ElementType==MshElemType::HEXAHEDRON)
      {
         Xi_p = 0.0;
         for (gp = 0; gp < nGauss; gp++)
         {
            r = MXPGaussPkt(nGauss, gp);
            if(fabs(r)>Xi_p) Xi_p = fabs(r);
         }
         r = 1.0/Xi_p;
         Xi_p = r;
      }
      //
      i_s=0;
      i_e=nnodes;
      ish=0;
      if(ElementType==MshElemType::TETRAHEDRON)   // tet
      {
         i_s=1;
         i_e=nnodes+1;
         ish=1;
      }
      //---------------------------------------------------------
      // Mapping Gauss point strains to nodes and update nodes
      // strains:
      //---------------------------------------------------------
      for(i=0; i<nnodes; i++)
      {
         eS = 0.0;
         SetExtropoGaussPoints(i);
         //
         ComputeShapefct(1);                      // Linear interpolation function
         for(j=i_s; j<i_e; j++)
            eS += NodalVal_Sat[j]*shapefct[j-ish];
         // Average value of the contribution of ell neighbor elements
         eS /= dbuff[i];
         eS += pcs->GetNodeValue(nodes[i],idx_S);
         // In case the node is on the material interface
         if(eS>1.0)
            eS = 1.0;
         //
         pcs->SetNodeValue (nodes[i], idx_S, eS);
      }
   }
   /**************************************************************************
   FEMLib-Method:
   Task: Caculate material parameter at element nodes for output
   Programing:
   04/2007 WW Implementation
   **************************************************************************/
   void CFiniteElementStd::CalcNodeMatParatemer()
   {
      int i, j, k, gp_r, gp_s, gp_t, idx_perm[3], idxp=0;
      int i_s, i_e, ish;
      double w[3], r=0.0, nval=0.0;
      //
      MshElemType::type ElementType = MeshElement->GetElementType();
      //----------------------------------------------------------------------
      gp = 0;
      index = Index;
      w[0] = w[1] = w[2] = 1.0;
      //----------------------------------------------------------------------
      setOrder(1);
      // Set material
      SetMaterial();
      //----------------------------------------------------------------------
      // Node value of the previous time step
      int idx11 = idx1;
      if(pcs->GetContinnumType()==1)
         idx11 = idxp21;
      for(i=0;i<nnodes;i++)
         NodalVal1[i] = pcs->GetNodeValue(nodes[i],idx11);
      if(PcsType==V)
      {
         for(i=0;i<nnodes;i++)
            NodalVal_p2[i] = pcs->GetNodeValue(nodes[i],idxp21);
      }
      if(PcsType==P)                              // 4.3.2009 PCH
      {
         for(i=0;i<nnodes;i++)
            NodalVal_SatNW[i] = pcs->GetNodeValue(nodes[i],idxSn1);
      }
      if(cpl_pcs)
      {
         for(i=0;i<nnodes;i++)
         {
            NodalValC[i] = cpl_pcs->GetNodeValue(nodes[i],idx_c0);
            NodalValC1[i] = cpl_pcs->GetNodeValue(nodes[i],idx_c1);
            if(cpl_pcs->type==1212)
               NodalVal_p2[i] = cpl_pcs->GetNodeValue(nodes[i],idx_c1+2);
                                                  //AKS
            NodalVal_p20[i] = pcs->GetNodeValue(nodes[i],idx_c0+2);
         }
      }
      //
      if((pcs->additioanl2ndvar_print>0)&&(pcs->additioanl2ndvar_print<3))
      {
         idx_perm[0] = pcs->GetNodeValueIndex("PERMEABILITY_X1");
         idx_perm[1] = pcs->GetNodeValueIndex("PERMEABILITY_Y1");
         if(dim==3)                               // 3D
            idx_perm[2] = pcs->GetNodeValueIndex("PERMEABILITY_Z1");
      }
      if(pcs->additioanl2ndvar_print>1)
         idxp = pcs->GetNodeValueIndex("POROSITY");
      // Number of elements associated to nodes
      for(i=0; i<nnodes; i++)
         dbuff[i] = (double)MeshElement->nodes[i]->connected_elements.size();
      //
      gp_r=gp_s=gp_t=gp=0;
      //
      for(gp=0; gp<nGaussPoints; gp++)
      {
         SetGaussPoint(gp, gp_r, gp_s, gp_t);
         if (ElementType==MshElemType::QUAD || ElementType==MshElemType::HEXAHEDRON)
         {
            i = GetLocalIndex(gp_r, gp_s, gp_t);
            if(i==-1) continue;
         }
         else
            i = gp;
         //
         if(i>nnodes) continue;
         ComputeShapefct(1);
         PG = interpolate(NodalVal1);
         //
         if((pcs->additioanl2ndvar_print>0)&&(pcs->additioanl2ndvar_print<3))
         {
            double* tensor = MediaProp->PermeabilityTensor(Index);
                                                  // Modified LBNL model
            if( MediaProp->permeability_stress_mode==2||MediaProp->permeability_stress_mode==3)
            {
               if(cpl_pcs)
                  TG = interpolate(NodalValC1)+T_KILVIN_ZERO;
               else
                  TG = 293.15;
               MediaProp->CalStressPermeabilityFactor(w, TG);
               for(j=0; j<dim; j++)
                  tensor[j*dim+j] *= w[j];
            }
            NodalVal2[i] = tensor[0];             // w[0];
            NodalVal3[i] = tensor[dim+1];         // w[1]; //
            if(dim==3)
               NodalVal4[i] = tensor[2*dim+2];    // w[2]; //
         }
         // Porosity
         if(pcs->additioanl2ndvar_print>1)
                                                  //MediaProp->Porosity(this);
            NodalVal0[i] = MediaProp->Porosity(MeshElement->index, 1.0);
      }
      //
      if (ElementType==MshElemType::QUAD || ElementType==MshElemType::HEXAHEDRON)
      {
         Xi_p = 0.0;
         for (gp = 0; gp < nGauss; gp++)
         {
            r = MXPGaussPkt(nGauss, gp);
            if(fabs(r)>Xi_p) Xi_p = fabs(r);
         }
         r = 1.0/Xi_p;
         Xi_p = r;
      }
      //
      i_s=0;
      i_e=nnodes;
      ish=0;
      if(ElementType==MshElemType::TETRAHEDRON)   // tet
      {
         i_s=1;
         i_e=nnodes+1;
         ish=1;
      }
      //---------------------------------------------------------
      // Mapping Gauss point strains to nodes and update nodes
      // strains:
      //---------------------------------------------------------
      for(i=0; i<nnodes; i++)
      {
         SetExtropoGaussPoints(i);
         //
         ComputeShapefct(1);                      // Linear interpolation function
         if((pcs->additioanl2ndvar_print>0)&&(pcs->additioanl2ndvar_print<3))
         {
            w[0] = w[1] = w[2] = 0.0;
            for(j=i_s; j<i_e; j++)
            {
               w[0] += NodalVal2[j]*shapefct[j-ish];
               w[1] += NodalVal3[j]*shapefct[j-ish];
               if(dim==3)
                  w[2] += NodalVal4[j]*shapefct[j-ish];
            }
            // Average value of the contribution of ell neighbor elements
            for(k=0; k<dim; k++)
            {
               w[k] /= dbuff[i];
               w[k] += pcs->GetNodeValue(nodes[i],idx_perm[k]);
               //
               pcs->SetNodeValue (nodes[i], idx_perm[k], w[k]);
            }
         }
         if(pcs->additioanl2ndvar_print>1)
         {
            nval = 0.0;
            for(j=i_s; j<i_e; j++)
               nval += NodalVal0[j]*shapefct[j-ish];
            nval /=  dbuff[i];
            nval += pcs->GetNodeValue(nodes[i],idxp);
            //
            pcs->SetNodeValue (nodes[i], idxp, nval);
         }
      }
   }

   //WW 08/2007
   ElementValue::ElementValue(CRFProcess* m_pcs, CElem* ele):pcs(m_pcs)
   {
      int NGPoints=0, NGP = 0;
      int ele_dim;

      MshElemType::type ele_type = ele->GetElementType();
      ele_dim = ele->GetDimension();

      NGP = GetNumericsGaussPoints(ele_type);
      if(ele_type==MshElemType::LINE)
                                                  //OKWW
            NGPoints = m_pcs->m_num->ele_gauss_points;
      else if(ele_type==MshElemType::TRIANGLE)
         NGPoints=3;
      else if(ele_type==MshElemType::TETRAHEDRON)
         NGPoints=15;
      else NGPoints = (int)pow((double)NGP, (double)ele_dim);

      //WW Velocity.resize(m_pcs->m_msh->GetCoordinateFlag()/10, NGPoints);
      Velocity.resize(3, NGPoints);
      Velocity = 0.0;
      if(pcs->type ==1212 || pcs->type == 1313)   // 15.3.2007 Multi-phase flow WW
      {
         Velocity_g.resize(3, NGPoints);
         Velocity_g = 0.0;
      }
   }
   //WW 08/2007
   void ElementValue::getIPvalue_vec(const int IP, double * vec)
   {
      for(int i=0; (size_t)i<Velocity.Rows(); i++) vec[i] = Velocity(i, IP);
   }

   /**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2006 YD Implementation
   last modification:
   **************************************************************************/
   void ElementValue::GetEleVelocity(double * vec)
   {
      for(int i=0; (size_t)i<Velocity.Rows(); i++)
      {
         vec[i] = 0.0;
         for(int j=0; (size_t)j<Velocity.Cols(); j++)
            vec[i] += Velocity(i, j);
         vec[i] /= Velocity.Cols();
      }
   }
   //WW
   ElementValue::~ElementValue()
   {
      Velocity.resize(0,0);
      Velocity_g.resize(0,0);
   }

   /**************************************************************************
   FEMLib-Method:
   01/2006 OK Implementation
   **************************************************************************/
   //void CFiniteElementStd::AssembleLHSMatrix()
   void CFiniteElementStd::AssembleParabolicEquationRHSVector()
   {
      int i;
      //----------------------------------------------------------------------
      // TIM
      double dt_inverse = 0.0;
      dt_inverse = 1.0 / dt;
      //----------------------------------------------------------------------
      // Initialize
      // if (pcs->Memory_Type==2) skip the these initialization
      (*Mass) = 0.0;
      (*Laplace) = 0.0;
      //----------------------------------------------------------------------
      // Calculate matrices
      // Mass matrix..........................................................
      if(pcs->m_num->ele_mass_lumping)
         CalcLumpedMass();
      else
         CalcMass();
      // Laplace matrix.......................................................
      CalcLaplace();
      //----------------------------------------------------------------------
      // Assemble local LHS matrix:
      // [C]/dt + theta [K]
      //Mass matrix
      *StiffMatrix    = *Mass;
      (*StiffMatrix) *= dt_inverse;
      // Laplace matrix
      *AuxMatrix      = *Laplace;
      *StiffMatrix   += *AuxMatrix;
      //----------------------------------------------------------------------
      for (i=0;i<nnodes; i++)
      {
         NodalVal1[i] = pcs->GetNodeValue(nodes[i],idx1);
         NodalVal[i] = 0.0;
      }
      //----------------------------------------------------------------------
      StiffMatrix->multi(NodalVal1, NodalVal);
      //----------------------------------------------------------------------
#ifdef NEW_EQS
      eqs_rhs = pcs->eqs_new->b;                  //WW
      if(m_dom)
         eqs_rhs = m_dom->eqs->b;                 //WW
#else
      eqs_rhs = pcs->eqs->b;                      //WW
#endif
      for (i=0;i<nnodes;i++)
      {
         eqs_number[i] = MeshElement->nodes[i]->GetEquationIndex();
         eqs_rhs[eqs_number[i]] +=  NodalVal[i];
      }
      //----------------------------------------------------------------------
   }

   ///////
   /**************************************************************************
   FEMLib-Method:
   Task: Calculate  coefficient of temperature induced RHS of multi-phase
         flow
   Programing:
   02/2007 WW Implementation
   last modification:
   **************************************************************************/
   inline double CFiniteElementStd::CalCoef_RHS_T_MPhase(int dof_index)
   {
      double val = 0.0, D_gw=0.0, D_ga=0.0;
      double expfactor=0.0,dens_arg[3];
      int Index = MeshElement->GetIndex();
      ComputeShapefct(1);
      //======================================================================
      switch(dof_index)
      {
         case 0:
            PG = interpolate(NodalVal1);
            Sw = MediaProp->SaturationCapillaryPressureFunction(PG,0);
            TG = interpolate(NodalValC1)+T_KILVIN_ZERO;
            TG0 = interpolate(NodalValC)+T_KILVIN_ZERO;
            PG2 = interpolate(NodalVal_p2);
            rhow = FluidProp->Density();
            poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
            expfactor = COMP_MOL_MASS_WATER/(rhow*GAS_CONSTANT*TG);
            rho_gw = FluidProp->vaporDensity(TG)*exp(-PG*expfactor);
            //
            drho_gw_dT = (FluidProp->vaporDensity_derivative(TG)
               +PG*expfactor*FluidProp->vaporDensity(TG)/TG)*exp(-PG*expfactor);
            val = (1.-Sw)*poro*drho_gw_dT/rhow;
            //
            if(SolidProp)
               val -= (1.0-poro)*((1-Sw)*rho_gw/rhow+Sw)*SolidProp->Thermal_Expansion();
            //
            // val += n*(1.0-rho_gw/rhow)*(dSw/dT)
            val *= (TG-TG0);
            break;
         case 1:
            //
            val = -(1.-Sw)*poro*drho_gw_dT/rhow;
            //
            if(SolidProp)
               val -= (1.0-poro)*(1-Sw)*rho_ga*SolidProp->Thermal_Expansion()/rhow;
            //
            // val -= n*rho_ga/rhow)*(dSw/dT)
            //---------------------------------------------------------------
            val *= (TG-TG0);
            break;
         case 2:
            //------------------------------------------------------------------------
            // From grad (p_gw/p_g)
            tort = MediaProp->TortuosityFunction(Index,unit,pcs->m_num->ls_theta);
            tort *=(1.0-Sw)*poro*2.16e-5*pow(TG/T_KILVIN_ZERO, 1.8);
            //
            p_gw = rho_gw*GAS_CONSTANT*TG/COMP_MOL_MASS_WATER;
            dens_arg[0] = PG2-p_gw;
            dens_arg[1] = TG;
            rho_ga = GasProp->Density(dens_arg);  //AKS SEP 2010  //(PG2-p_gw)*GasProp->molar_mass/(GAS_CONSTANT*TG);
            rho_g = rho_ga+rho_gw;
            // 1/Mg
            M_g = (rho_gw/COMP_MOL_MASS_WATER+rho_ga/GasProp->molar_mass)/rho_g;
            D_gw = tort*rho_g*COMP_MOL_MASS_WATER*GasProp->molar_mass*M_g*M_g/rhow;
            val = D_gw*drho_gw_dT*GAS_CONSTANT*TG/(COMP_MOL_MASS_WATER*PG2)*time_unit_factor;
            break;
         case 3:
            //---------------------------------------------------------------
            //
            D_ga = tort*rho_g*COMP_MOL_MASS_WATER*GasProp->molar_mass*M_g*M_g/rhow;
            // From grad (p_gw/p_g)
            val = -D_ga*drho_gw_dT*GAS_CONSTANT*TG/(COMP_MOL_MASS_WATER*PG2)*time_unit_factor;

            break;
            //------------------------------------------------------------------
      }
      return val;
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate coefficient of temperature induced RHS of PSGlobal scheme

   Programing:

   last modification:
   **************************************************************************/
   inline double CFiniteElementStd::CalCoef_RHS_T_PSGlobal(int dof_index)
   {
      double val = 0.0;                           //OK411 D_gw=0.0, D_ga=0.0;
      //OK411 double expfactor=0.0;
      double P,T;
      int Index = MeshElement->GetIndex();
      ComputeShapefct(1);
      //======================================================================
      switch(dof_index)
      {
         case 0:
            val = 0.;
            break;
         case 1:
            poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
            Sw = 1.0-interpolate(NodalVal_SatNW);
                                                  // Pnw = Pw + Pc(Sw)
            P  = interpolate(NodalVal1) + MediaProp->CapillaryPressureFunction(0,NULL,0.0,1,Sw);
            // 	P  = interpolate(NodalVal1);  // Pw
            T  = interpolate(NodalValC1);

            val = -(1.-Sw)*poro*GasProp->drhodT(P,T)/GasProp->Density();
            break;
         case 2:
            val = 0.;
            break;
         case 3:
            val = 0.;
            break;
            //------------------------------------------------------------------
      }
      return val;
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate  Capillary pressure for RHS in the global scheme
   Programing:
   03/2009 PCH Implementation
   last modification:
   **************************************************************************/
   inline void CFiniteElementStd::CalCoef_RHS_Pc(int dof_index)
   {
      int i=0;
      double *tensor = NULL;
      double mat_fac = 1.0;                       //OK411 m_fac=0.;
      double k_rel=0.0;

      int Index = MeshElement->GetIndex();
      //
      ComputeShapefct(1);                         //  12.3.2007 WW
      //======================================================================
      for(i=0; i<dim*dim; i++)
         mat[i] = 0.0;

      switch(dof_index)
      {
         case 0:
            break;
         case 1:
            break;
         case 2:
            tensor = MediaProp->PermeabilityTensor(Index);
            Sw = 1.0-interpolate(NodalVal_SatNW);
            k_rel = MediaProp->PermeabilitySaturationFunction(Sw,1);
            mat_fac = k_rel / GasProp->Viscosity()*time_unit_factor;

            for(i=0;i<dim*dim;i++)
               mat[i] = tensor[i] * mat_fac;
            break;
         case 3:
            break;
            //------------------------------------------------------------------
      }
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate  coefficient of Pc induced RHS of multi-phase
         flow
   Programing:
   03/2009 PCH Implementation
   last modification:
   **************************************************************************/
   inline double CFiniteElementStd::CalCoef_RHS_PSGLOBAL(int dof_index)
   {
      double val = 0.0;
      double k_rel=0.0;                           //OK411 mat_fac=0.0;
      int i=0;
      //OK411 int Index = MeshElement->GetIndex();
      //OK411 double *tensor = NULL;
      ComputeShapefct(1);
      for(i=0; i<dim*dim; i++)
         mat[i] = 0.0;

      switch(dof_index)
      {
         case 0:
            val = 0.0;
            break;
         case 1:
            val = 0.0;
            break;
         case 2:
            Sw = 1.0-interpolate(NodalVal_SatNW);
            k_rel = MediaProp->PermeabilitySaturationFunction(Sw,1);
            val = k_rel / GasProp->Viscosity()*time_unit_factor;
            break;
         case 3:
            val = 0.0;
            break;
            //------------------------------------------------------------------
      }
      return val;
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate right hand terms temperature coupled term and body force
   Programing:
   05/2010 AKS Implementation
   last modification:
   **************************************************************************/

   inline double CFiniteElementStd::CalCoef_RHS_AIR_FLOW(int dof_index)
   {
      double val=0.0;
      int Index = MeshElement->GetIndex();
      PG=interpolate(NodalVal1);
      TG=interpolate(NodalValC1)+T_KILVIN_ZERO;
      TG0=interpolate(NodalValC)+T_KILVIN_ZERO;
      switch(dof_index)
      {
         case 0:
            val = -MediaProp->Porosity(Index,pcs->m_num->ls_theta)/TG;
            val *= (TG-TG0);
            break;

         case 1:
            val = -1.0/TG;
            break;

      }
      return val;
   }
   /**************************************************************************
   FEMLib-Method:
   Task: Calculate RHS of pressure coupled term
   Programing:
   05/2010 AKS Implementation
   last modification:
   **************************************************************************/

   inline double CFiniteElementStd::CalCoef_RHS_HEAT_TRANSPORT(int dof_index)
   {
      double val=0.0, rho_g=0.0, rho_0=0.0;
      int Index = MeshElement->GetIndex();
      double dens_arg[3];
      dens_arg[0] = interpolate(NodalValC1);
      dens_arg[1] = interpolate(NodalVal1)+T_KILVIN_ZERO;
      dens_arg[2] = Index;
      rho_g=FluidProp->Density(dens_arg);
      dens_arg[0] = 4.0e6;
      dens_arg[1] = 120+T_KILVIN_ZERO;
      rho_0=FluidProp->Density(dens_arg);

      switch(dof_index)
      {
         case 0:
            val = (interpolate(NodalValC1)-interpolate(NodalValC))*MediaProp->Porosity(Index,pcs->m_num->ls_theta)*rho_g/rho_0;
            break;

         case 1:
            val = rho_g/rho_0;
            val -=1.0;                            // term coresponding to the 'Viscour dissipation'
            break;

      }
      return val;
   }
   /***************************************************************************
      GeoSys - Funktion:
             Assemble_RHS_T_MPhaseFlow
      Programming:
      02/2007   WW
   **************************************************************************/
   void CFiniteElementStd::Assemble_RHS_T_MPhaseFlow()
   {
      int i, j, k, ii;
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;
      double fkt, fac;
      // Material
      int dof_n = 2;
      //----------------------------------------------------------------------
      for (i = 0; i < dof_n*nnodes; i++) NodalVal[i] = 0.0;
      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         ComputeGradShapefct(1);                  // Linear interpolation function
         ComputeShapefct(1);                      // Linear interpolation function
         for(ii=0; ii<dof_n; ii++)
         {
            // Material
            fac = fkt*CalCoef_RHS_T_MPhase(ii)/dt;
            // Calculate THS
            for (i = 0; i < nnodes; i++)
               NodalVal[i+ii*nnodes] += fac *shapefct[i];
         }
         // grad T
         for(ii=0; ii<dof_n; ii++)
         {
            // Material
            fac = fkt*CalCoef_RHS_T_MPhase(ii+dof_n);
            // Calculate THS
            for (i = 0; i < nnodes; i++)
            {
               for (j = 0; j < nnodes; j++)
               {
                  for (k = 0; k < dim; k++)
                     NodalVal[i+ii*nnodes] +=
                        fac*dshapefct[k*nnodes+i]*dshapefct[k*nnodes+j]
                        *(NodalValC1[j]+T_KILVIN_ZERO);
               }
            }
         }
      }
      int ii_sh;
      long i_sh;
      for(ii=0;ii<pcs->dof;ii++)
      {
         i_sh = NodeShift[ii];
         ii_sh = ii*nnodes;
         for (i=0;i<nnodes;i++)
         {
            eqs_rhs[i_sh + eqs_number[i]] -= NodalVal[i+ii_sh];
            (*RHS)(i+LocalShift+ii_sh) -=  NodalVal[i+ii_sh];
         }
      }
      //
   }

   /***************************************************************************
   GeoSys - Function: Assemble_RHS_T_PSGlobal
   Programming:
    09/2009
   **************************************************************************/
   void CFiniteElementStd::Assemble_RHS_T_PSGlobal()
   {
      int i, ii;                                  //OK411 j, k,
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;
      double fkt, fac;
      // Material
      int dof_n = 2;
      double temp[10];

      for (i=0;i<10;i++) temp[i]=0;               // remove

      //----------------------------------------------------------------------
      for (i = 0; i < dof_n*nnodes; i++) NodalVal[i] = 0.0;
      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         ComputeGradShapefct(1);                  // Linear interpolation function
         ComputeShapefct(1);                      // Linear interpolation function

         for(ii=0; ii<dof_n; ii++)
         {
            // Material
            fac = fkt*CalCoef_RHS_T_PSGlobal(ii)/dt;
            // Calculate THS
            for (i = 0; i < nnodes; i++)
            {
               NodalVal[i+ii*nnodes] += fac *shapefct[i];
                                                  //remove
               temp[i+ii*nnodes] += fac *shapefct[i];
            }
         }

      }
      int ii_sh;
      long i_sh;
      for(ii=0;ii<pcs->dof;ii++)
      {
         i_sh = NodeShift[ii];
         ii_sh = ii*nnodes;
         for (i=0;i<nnodes;i++)
         {
            eqs_rhs[i_sh + eqs_number[i]] -= NodalVal[i+ii_sh];
            (*RHS)(i+LocalShift+ii_sh) -=  NodalVal[i+ii_sh];
         }
      }
      //
   }

   /***************************************************************************
      GeoSys - Funktion:
             Assemble_RHS_Pc
      Programming:
      03/2009   PCH
   **************************************************************************/
   void CFiniteElementStd::Assemble_RHS_Pc()
   {
      int i, j, k, l, ii;
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;
      double fkt;
      // Material
      int dof_n = 2;
      int ndx_p_cap = pcs->GetNodeValueIndex("PRESSURE_CAP");
      //----------------------------------------------------------------------

      double temp[20];

      for (i = 0; i < dof_n*nnodes; i++)
      {
         temp[i] = NodalVal[i] = 0.0;
         NodalVal1[i] = 0.0;
      }
      for (i = 0; i < nnodes; i++)
         temp[i+dof_n] = NodalVal1[i+dof_n] = -pcs->GetNodeValue(nodes[i],ndx_p_cap);

      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         ComputeGradShapefct(1);                  // Linear interpolation function
         ComputeShapefct(1);                      // Linear interpolation function
         // grad Pc
         for(ii=0; ii<dof_n; ii++)
         {
            // Material
            CalCoef_RHS_Pc(ii+dof_n);
            // Calculate Pc
            for (i = 0; i < nnodes; i++)
            {
               for (j = 0; j < nnodes; j++)
               {
                  for (k = 0; k < dim; k++)
                     for (l = 0; l < dim; l++)
                        NodalVal[dof_n+i] += fkt*mat[dim*k+l]*dshapefct[k*nnodes+i]*dshapefct[l*nnodes+j]
                           *NodalVal1[j+dof_n];
               }
            }
         }
      }

      for(i=0; i<2*nnodes; ++i)
         temp[i]=NodalVal[i];

      int ii_sh;
      long i_sh;
      for(ii=0;ii<pcs->dof;ii++)
      {
         i_sh = NodeShift[ii];
         ii_sh = ii*nnodes;
         for (i=0;i<nnodes;i++)
         {
            eqs_rhs[i_sh + eqs_number[i]] += NodalVal[i+ii_sh];
            (*RHS)(i+LocalShift+ii_sh) +=  NodalVal[i+ii_sh];
         }
      }
      //
   }

   /***************************************************************************
      GeoSys - Funktion:
             Assemble_RHS_M
      Programming:
      02/2007   WW
   **************************************************************************/
   void CFiniteElementStd::Assemble_RHS_M()
   {
      int i, ii;
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;
      double fkt, fac, grad_du=0.0;
      // Material
      int dof_n = 2;
      //----------------------------------------------------------------------
      for (i = 0; i < dof_n*nnodes; i++) NodalVal[i] = 0.0;
      for (i=nnodes;i<nnodesHQ;i++)
         nodes[i] = MeshElement->nodes_index[i];
      for (i=0;i<nnodesHQ;i++)
      {
         NodalVal2[i] = ( dm_pcs->GetNodeValue(nodes[i],Idx_dm1[0])
            -dm_pcs->GetNodeValue(nodes[i],Idx_dm0[0]));
         NodalVal3[i] = ( dm_pcs->GetNodeValue(nodes[i],Idx_dm1[1])
            -dm_pcs->GetNodeValue(nodes[i],Idx_dm0[1]));
         if(dim==3)                               // 3D.
            NodalVal4[i] = ( dm_pcs->GetNodeValue(nodes[i],Idx_dm1[2])
               -dm_pcs->GetNodeValue(nodes[i],Idx_dm0[2]));
      }
      //======================================================================
      SetHighOrderNodes();
      //
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         ComputeShapefct(1);                      // Linear interpolation function
         //ComputeShapefct(2);
         ComputeGradShapefct(2);
         grad_du = 0.0;
         for (i=0;i<nnodesHQ;i++)
         {
            grad_du += dshapefctHQ[i]*NodalVal2[i]+dshapefctHQ[i+nnodesHQ]*NodalVal3[i];
            if(dim==3)                            // 3D.
               grad_du += dshapefctHQ[i+nnodesHQ*2]*NodalVal4[i];
         }
         grad_du /= dt;
         for(ii=0; ii<dof_n; ii++)
         {
            // Material
            fac = fkt*grad_du*CalCoef_RHS_M_MPhase(ii);
            // Calculate MHS
            for (i = 0; i < nnodes; i++)
               NodalVal[i+ii*nnodes] += fac *shapefct[i];
         }
      }
      //
      int ii_sh;
      long i_sh;
      for(ii=0;ii<pcs->dof;ii++)
      {
         i_sh = NodeShift[ii];
         ii_sh = ii*nnodes;
         for (i=0;i<nnodes;i++)
         {
            eqs_rhs[i_sh + eqs_number[i]] -= NodalVal[i+ii_sh];
            (*RHS)(i+LocalShift+ii_sh) -=  NodalVal[i+ii_sh];
         }
      }
      setOrder(1);
      //
   }

   /***************************************************************************
      GeoSys - Funktion:
      Assemble_RHS_AIR_FLOW
      Programming:
       05/2010 AKS
   **************************************************************************/

   void CFiniteElementStd::Assemble_RHS_AIR_FLOW()
   {
      int i, j, k, ii;                            //KR,idxd;
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;                   //KR ,z_sum;
      double vel[3];                              //KR,rhoz[3];
      double fkt, fac,mat_fac,fluid_density;
      double dens_arg[3];                         //08.05.2008 WW
      double *tensor = NULL;
      //KR CFEMesh* m_msh;
      int GravityOn = 1;                          // Initialized to be on
      // If no gravity, then set GravityOn to be zero.
      if((coordinate_system)%10!=2&&(!axisymmetry))
         GravityOn = 0;
      // Material
      int dof_n = 1;
      //----------------------------------------------------------------------
      for (i = 0; i < dof_n*nnodes; i++) NodalVal[i] = 0.0;
      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         ComputeGradShapefct(1);                  // Linear interpolation function
         ComputeShapefct(1);                      // Linear interpolation function
         ElementValue* gp_ele = ele_gp_value[Index];

         for(ii=0; ii<dof_n; ii++)
         {
            // Material
            fac = CalCoef_RHS_AIR_FLOW(ii)/dt;

            for (i = 0; i < nnodes; i++)
               NodalVal[i+ii*nnodes] += fac*fkt*shapefct[i];
         }

         // grad T
         for(ii=0; ii<dof_n; ii++)
         {

            fac = CalCoef_RHS_AIR_FLOW(ii+1);
            //Velocity
            vel[0] = fac*gp_ele->Velocity(0, gp);
            vel[1] = fac*gp_ele->Velocity(1, gp);
            vel[2] = fac*gp_ele->Velocity(2, gp);

            for (i = 0; i< nnodes; i++)
            {
               for (j = 0; j < nnodes; j++)
               {
                  for (k = 0; k < dim; k++)
                     NodalVal[i+ii*nnodes] += fkt*shapefct[i]*vel[k]*dshapefct[k*nnodes+j]*(NodalValC1[j]+T_KILVIN_ZERO);
               }
            }
         }

         //Body force term
         if(GravityOn)
         {
            dens_arg[0] = interpolate(NodalVal1);
            dens_arg[1] = interpolate(NodalValC1)+T_KILVIN_ZERO;
            dens_arg[2] = Index;
            fluid_density=FluidProp->Density(dens_arg);
            mat_fac = FluidProp->Viscosity(dens_arg);
            tensor = MediaProp->PermeabilityTensor(Index);
            for(i=0;i<dim*dim;i++)
               mat[i] = tensor[i]/mat_fac;
            for(ii=0; ii<dof_n; ii++)
            {

               for (i = 0; i < nnodes; i++)
               {
                  for (k = 0; k < dim; k++)
                     NodalVal[i+ii*nnodes] -= fkt*fluid_density*gravity_constant*mat[dim*k+dim-1]*dshapefct[k*nnodes+i];
               }
            }
         }

      }
      int ii_sh;
      long i_sh;
      for(ii=0;ii<pcs->dof;ii++)
      {
         i_sh = NodeShift[ii];
         ii_sh = ii*nnodes;
         for (i=0;i<nnodes;i++)
         {
            eqs_rhs[i_sh + eqs_number[i]] -= NodalVal[i+ii_sh];
            (*RHS)(i+LocalShift+ii_sh) -=  NodalVal[i+ii_sh];
         }
      }
   }

   /***************************************************************************
      GeoSys - Funktion:
      Assemble_RHS_HEAT_TRANSPORT: This include when need pressure terms n dp/dt + nv.Nabla p
      Programming:
      05/2010   AKS
   **************************************************************************/

   void CFiniteElementStd::Assemble_RHS_HEAT_TRANSPORT()
   {
      int i, j, k, ii;
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;
      double vel[3];
      double fkt=0.0, fac=0.0;
      // Material
      int dof_n = 1;
      //----------------------------------------------------------------------
      for (i = 0; i < dof_n*nnodes; i++) NodalVal[i] = 0.0;
      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         ComputeGradShapefct(1);                  // Linear interpolation function
         ComputeShapefct(1);                      // Linear interpolation function
         ElementValue* gp_ele = ele_gp_value[Index];

         for(ii=0; ii<dof_n; ii++)
         {
            // Material
            fac = CalCoef_RHS_HEAT_TRANSPORT(ii)/dt;

            for (i = 0; i < nnodes; i++)
               NodalVal[i+ii*nnodes] += fac*fkt*shapefct[i];

         }

         // grad P

         for(ii=0; ii<dof_n; ii++)
         {
            fac = CalCoef_RHS_HEAT_TRANSPORT(ii+1);
            //Velocity
            vel[0] = fac*gp_ele->Velocity(0, gp);
            vel[1] = fac*gp_ele->Velocity(1, gp);
            vel[2] = fac*gp_ele->Velocity(2, gp);

            for (i = 0; i< nnodes; i++)
            {
               for (j = 0; j < nnodes; j++)
               {
                  for (k = 0; k < dim; k++)
                     NodalVal[i+ii*nnodes] += fkt*vel[k]*shapefct[i]*dshapefct[k*nnodes+j]*NodalValC1[j];
               }
            }
         }

      }
      int ii_sh;
      long i_sh;
      for(ii=0;ii<pcs->dof;ii++)
      {
         i_sh = NodeShift[ii];
         ii_sh = ii*nnodes;
         for (i=0;i<nnodes;i++)
         {
            eqs_rhs[i_sh + eqs_number[i]] -= NodalVal[i+ii_sh];
            (*RHS)(i+LocalShift+ii_sh) -=  NodalVal[i+ii_sh];
         }
      }
   }
   /**************************************************************************
   FEMLib-Method:
   Task: Calculate RHS of pressure coupled term
   Programing:
   05/2010 AKS Implementation
   last modification:
   **************************************************************************/

   inline double CFiniteElementStd::CalCoef_RHS_HEAT_TRANSPORT2(int dof_index)
   {
      int i;
      // TF unused variable - comment fix compile warning
//      ElementValue* gp_ele = ele_gp_value[Index];
      double *tensor = NULL;
      double val = 0.0,mat_fac;
      // TF unused variable - comment fix compile warning
//      double Tc=647.096;
      double H_vap=0.0,dens_arg[3];
      ComputeShapefct(1);
      PG = interpolate(NodalValC1);
      PG2 = interpolate(NodalVal_p2);
      TG = interpolate(NodalVal1)+T_KILVIN_ZERO;
      PG0 = interpolate(NodalValC);
      PG20 = interpolate(NodalVal_p20);
      dens_arg[1] = TG;
      dens_arg[0] = PG2-PG;
      rhow = FluidProp->Density(dens_arg);
      Sw =MediaProp->SaturationCapillaryPressureFunction(PG,0);
      dSdp = -MediaProp->SaturationPressureDependency(Sw, rhow, pcs->m_num->ls_theta);
      poro = MediaProp->Porosity(Index,pcs->m_num->ls_theta);
      if(MediaProp->evaporation==647)
         H_vap = -2257000;                        //pow((Tc - TG),0.38)*2.5397E+5;//It is specific you can change thi value as you chaning fluid from water
      for (i=0;i<dim*dim;i++) mat[i]=0.0;
      switch(dof_index)
      {
         case 0:
            val = H_vap*rhow*poro*dSdp;
            val *= (PG-PG0);
            return val;
            break;

         case 1:
            val=0.0;
            return val;
            break;

         case 2:
            tensor = MediaProp->PermeabilityTensor(Index);
            mat_fac = H_vap*rhow*MediaProp->PermeabilitySaturationFunction(Sw,0)/FluidProp->Viscosity();
            for(i=0;i<dim*dim;i++)
               mat[i] = tensor[i]*mat_fac*time_unit_factor;
            break;

         case 3:
            tensor = MediaProp->PermeabilityTensor(Index);
            mat_fac = -H_vap*rhow*MediaProp->PermeabilitySaturationFunction(Sw,0)/FluidProp->Viscosity();
            for(i=0;i<dim*dim;i++)
               mat[i] = tensor[i]*mat_fac*time_unit_factor;
            break;

      }
      return 0.;                                  //WW
   }
   /***************************************************************************
      GeoSys - Funktion:
      Assemble_RHS_HEAT_TRANSPORT2: This include when need pressure terms n dp/dt + nv.Nabla p
      Programming:
      05/2010   AKS
   **************************************************************************/

   void CFiniteElementStd::Assemble_RHS_HEAT_TRANSPORT2()
   {
      int i, j, k, ii;
      // ---- Gauss integral
      int gp_r=0,gp_s=0,gp_t=0;
      // TF unused variable - comment fix compile warning
//      double *tensor = NULL,
      double dens_arg[3];
      // TF unused variable - comment fix compile warning
//      double H_vap=0;
      // TF unused variable - comment fix compile warning
//      double Tc=647.096;
      double fkt=0.0, fac=0.0;                    //WW,mat_fac;
      // Material
      int dof_n = 1;
      //----------------------------------------------------------------------
      for (i = 0; i < dof_n*nnodes; i++) NodalVal[i] = 0.0;
      //======================================================================
      // Loop over Gauss points
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determinate
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         // Compute geometry
         ComputeGradShapefct(1);                  // Linear interpolation function
         ComputeShapefct(1);                      // Linear interpolation function
         // TF unused variable - comment fix compile warning
//         ElementValue* gp_ele = ele_gp_value[Index];
         int dof_n = 2;
         int GravityOn = 1;                       // Initialized to be on
         // If no gravity, then set GravityOn to be zero.
         if((coordinate_system)%10!=2&&(!axisymmetry))
            GravityOn = 0;
         TG = interpolate(NodalVal1)+T_KILVIN_ZERO;
         PG=interpolate(NodalValC1);
         PG2=interpolate(NodalVal_p2);
         dens_arg[1] = TG;
         dens_arg[0] = PG2-PG;

         for(ii=0; ii<dof_n; ii++)
         {
            // Material
            fac = fkt*CalCoef_RHS_HEAT_TRANSPORT2(ii)/dt;
            // Calculate THS
            for (i = 0; i < nnodes; i++)
               NodalVal[i] += fac *shapefct[i];
         }

         // grad pc
         for(ii=0; ii<dof_n-1; ii++)
         {
            // Material
            CalCoef_RHS_HEAT_TRANSPORT2(ii+dof_n);
            for (i = 0; i< nnodes; i++)
            {
               for (j = 0; j < nnodes; j++)
               {
                  for (k = 0; k < dim; k++)
                     NodalVal[i]  += fkt*mat[dim*k+k]*dshapefct[k*nnodes+i]*dshapefct[k*nnodes+j]*NodalValC1[j];
               }
            }
         }

         // grad pc
         for(ii=0; ii<dof_n-1; ii++)
         {
            // Material
            CalCoef_RHS_HEAT_TRANSPORT2(ii+dof_n+1);
            for (i = 0; i< nnodes; i++)
            {
               for (j = 0; j < nnodes; j++)
               {
                  for (k = 0; k < dim; k++)
                     NodalVal[i] += fkt*mat[dim*k+k]*dshapefct[k*nnodes+i]*dshapefct[k*nnodes+j]*NodalVal_p2[j];
               }
            }
         }

         //gravity
         if(GravityOn)
         {
            CalCoef_RHS_HEAT_TRANSPORT2(2);
            for (i = 0; i < nnodes; i++)
            {
               for (k = 0; k < dim; k++)
                  NodalVal[i] -= fkt*mat[dim*k+dim-1]*FluidProp->Density(dens_arg)*gravity_constant*dshapefct[k*nnodes+i];
            }
         }

      }
      int ii_sh;
      long i_sh;
      for(ii=0;ii<pcs->dof;ii++)
      {
         i_sh = NodeShift[ii];
         ii_sh = ii*nnodes;
         for (i=0;i<nnodes;i++)
         {
            eqs_rhs[i_sh + eqs_number[i]] -= NodalVal[i+ii_sh];
            (*RHS)(i+LocalShift+ii_sh) -=  NodalVal[i+ii_sh];
         }
      }
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate  coefficient of displacement induced RHS of multi-phase
         flow
   Programing:
   02/2007 WW Implementation
   05/2008 WW Generalization
   last modification:
   **************************************************************************/
   inline double CFiniteElementStd::CalCoef_RHS_M_MPhase(int dof_index)
   {
      double val = 0.0;
      double expfactor=0.0;
      double dens_aug[3];
      dens_aug[1] = 293.15;
      bool diffusion = false;                     //08.05.2008 WW
      if(MediaProp->heat_diffusion_model==273&&cpl_pcs)
         diffusion = true;
      //======================================================================
      switch(dof_index)
      {
         case 0:
            PG = interpolate(NodalVal1);
            Sw = MediaProp->SaturationCapillaryPressureFunction(PG,0);
            if(diffusion)
            {
               TG = interpolate(NodalValC1)+T_KILVIN_ZERO;
               dens_aug[1] = TG;
            }
            //
            dens_aug[0] = PG;
            rhow = FluidProp->Density(dens_aug);
            val = Sw;
            //
            PG2 = interpolate(NodalVal_p2);
            dens_aug[0] = PG2;
            //
            if(diffusion)
            {
               expfactor = COMP_MOL_MASS_WATER/(rhow*GAS_CONSTANT*TG);
               rho_gw = FluidProp->vaporDensity(TG)*exp(-PG*expfactor);
               p_gw = rho_gw*GAS_CONSTANT*TG/COMP_MOL_MASS_WATER;
               dens_aug[0] -= p_gw;
            }
            rho_ga = GasProp->Density(dens_aug);
            if(diffusion)
               val += (1.0-Sw)*rho_gw/rhow;
            break;
         case 1:
            val =  (1.0-Sw)*rho_ga/rhow;
            break;
            //------------------------------------------------------------------
      }
      return val;
   }

   /**************************************************************************
   PCSLib-Method:
   01/2007 OK Implementation
   **************************************************************************/

   /**************************************************************************
   PCSLib-Method:
   01/2007 OK Implementation
   02/2009 PCH modified to handle known Snw in the mass matrix
   **************************************************************************/
   void CFiniteElementStd::AssembleRHSVector()
   {
      int i;
      int idx_fv=0, idx_pw = 0, idx_pc=0;
      double NodalVal_FV[20];
      //OK411 double FV;
      CRFProcess* pcs_p = NULL;
      CRFProcess* pcs_s = NULL;
      //----------------------------------------------------------------------
      // Initializations
      for(i=0;i<nnodes;i++)
         NodalVal[i] = 0.0;

      double temp[8];

      switch(PcsType)
      {
         //....................................................................
         case T:                                  // Two-phase flow
            if(pcs->PartialPS == 0)               // If not partial-pressure-based
               (*Laplace) = 0.0;
            else
               (*Mass) = 0.0;
            break;
         default:
            break;
            //....................................................................
      }
      //----------------------------------------------------------------------
      // Field variables
      switch(PcsType)
      {
         //....................................................................
         case T:                                  // Two-phase flow
            if(pcs->PartialPS == 0)
            {
               pcs_p = pcs_vector[0];
               pcs_s = pcs_vector[1];

               idx_pw = pcs_p->GetNodeValueIndex("PRESSURE1");
               idx_pc = pcs_p->GetNodeValueIndex("PRESSURE_CAP");
               idx_fv = pcs_s->GetNodeValueIndex("SATURATION2");
               CMediumProperties *m_mmp = NULL;
               m_mmp = mmp_vector[0];
               for(i=0;i<nnodes;i++)
               {
                  Sw = 1.0 - pcs_s->GetNodeValue(nodes[i],idx_fv+1);
                  double Pw = pcs_p->GetNodeValue(nodes[i],idx_pw+1);
                  double Pc = pcs_p->GetNodeValue(nodes[i],idx_pc+1);
                  if(pcs->ML_Cap == 0)            // If ODE method for Snw,
                     NodalVal_FV[i] = -(Pw + Pc);
                  else                            // If PDE method,
                     NodalVal_FV[i] = -Pw;
               }
            }
            else
            {
            }
            break;
         default:
            break;
            //....................................................................
      }

      //----------------------------------------------------------------------
      // Element matrices
      switch(PcsType)
      {
         //....................................................................
         case T:                                  // Two-phase flow
            if(pcs->PartialPS == 0)
               CalcLaplace();
            else
               CalcMass();
            break;
         default:
            break;
            //....................................................................
      }
      //----------------------------------------------------------------------
      // Calc RHS contribution
      switch(PcsType)
      {
         //....................................................................
         case T:                                  // Two-phase flow
            if(pcs->PartialPS == 0)
               Laplace->multi(NodalVal_FV,NodalVal);
            else
               Mass->multi(NodalVal_FV,NodalVal);
            break;
         default:
            break;
            //....................................................................
      }

      for(i=0;i<nnodes;i++)
         temp[i]=NodalVal[i];

      //----------------------------------------------------------------------
      // Store RHS contribution
      for(i=0;i<nnodes;i++)
      {
         //CB 04008
#ifdef NEW_EQS
         pcs->eqs_new->b[NodeShift[problem_dimension_dm]+eqs_number[i]] += NodalVal[i];
#else
         pcs->eqs->b[NodeShift[problem_dimension_dm]+eqs_number[i]] += NodalVal[i];
#endif
         (*RHS)(i+LocalShift) += NodalVal[i];
      }
      //----------------------------------------------------------------------
      //RHS->Write();
   }

   /**************************************************************************
   PCSLib-Method:
   09/2008 PCH Implementation
   **************************************************************************/
   void CFiniteElementStd::AssembleCapillaryEffect()
   {
      int i;
      int idx_pc=0;                               //OK411 idx_w=0, idx_pw=0, idx_fv=0, idx_nw=0,
      double NodalVal_FV[20];
      //OK411 double FV;
      //OK411 CRFProcess* m_pcs_cpl = NULL;
      CRFProcess* pcs_p = NULL;
      //OK411 CRFProcess* pcs_s = NULL;
      //OK411 CMediumProperties *m_mmp = NULL;

      //----------------------------------------------------------------------
      // Initializations
      for(i=0;i<nnodes;i++)
         NodalVal[i] = 0.0;

      // Setting the laplace matrix initialized
      (*Laplace) = 0.0;

      if(pcs->pcs_type_number==0)
      {
         pcs_p = pcs_vector[0];

         idx_pc = pcs_p->GetNodeValueIndex("PRESSURE_CAP");

         for(i=0;i<nnodes;i++)
         {
            double Pc = pcs_p->GetNodeValue(nodes[i],idx_pc+1);
            NodalVal_FV[i] = -Pc;
         }
      }
      else if(pcs->pcs_type_number==1)
      {
      }
      else
         printf("Something's wrong!!!\n");

      //----------------------------------------------------------------------
      // Element matrices and RHS calculation
      if(pcs->pcs_type_number==0)
      {
         // Laplace should be calculated according to the nonwetting phase
         // Then I need one switch to tell CalcCoefLaplace to use the nonwetting parameter only
         pcs->ML_Cap = 1;
         CalcLaplace();
         Laplace->multi(NodalVal_FV,NodalVal);
         pcs->ML_Cap = 0;
      }
      else if(pcs->pcs_type_number==1)
      {
         CalcLaplace();
         Laplace->multi(NodalVal_FV,NodalVal);
      }
      else
         printf("Something's wrong!!!\n");

      //----------------------------------------------------------------------
      // Store RHS contribution
      for(i=0;i<nnodes;i++)
      {
         //CB 04008
#ifdef NEW_EQS
         pcs->eqs_new->b[NodeShift[problem_dimension_dm]+eqs_number[i]] += NodalVal[i];
#else
         pcs->eqs->b[NodeShift[problem_dimension_dm]+eqs_number[i]] += NodalVal[i];
#endif
         (*RHS)(i+LocalShift) += NodalVal[i];
      }
      //----------------------------------------------------------------------
      //RHS->Write();
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate the energy norm error
   Programing:
   25.08.2008 WW Implementation
   last modification:
   **************************************************************************/
   void CFiniteElementStd::CalcEnergyNorm(const double *x_n1, double &err_norm0,
      double &err_normn)
   {
      int i, dof_n = 1;
      // NUM
      double rtol, atol;
      //----------------------------------------------------------------------
      //
      Config();
      //

      //
      rtol = pcs->Tim->GetRTol();
      atol = pcs->Tim->GetATol();
      //_new
      for(i=0; i<pcs->dof; i++)
         NodeShift[i]=i*pcs->m_msh->GetNodesNumber(false);

      //----------------------------------------------------------------------
      //  double beta1 = 0.0;
      //----------------------------------------------------------------------
      // Initialize.
      // if (pcs->Memory_Type==2) skip the these initialization
      if(PcsType==V || PcsType==P)                // 03.2009 PCH
         (*Mass2) = 0.0;
      else
         (*Mass) = 0.0;
      (*Laplace) = 0.0;
      //----------------------------------------------------------------------
      // GEO
      // double geo_fac = MediaProp->geo_area;
      //----------------------------------------------------------------------
      // Calculate matrices
      // Mass matrix..........................................................
      if(PcsType==V)
      {
         if(pcs->m_num->ele_mass_lumping)
            CalcLumpedMass2();
         else
            CalcMass2();
      }
      else if(PcsType==P)                         // 03.2009 PCH
      {
         if(pcs->m_num->ele_mass_lumping)
            CalcLumpedMassPSGLOBAL();
         else
            CalcMassPSGLOBAL();
      }
      else
      {
         if(pcs->m_num->ele_mass_lumping)
            CalcLumpedMass();
         else
            CalcMass();
      }
      // Laplace matrix.......................................................
      CalcLaplace();
      if(PcsType==V || PcsType==P)                // 03.2009 PCH
         *AuxMatrix1 = *Mass2;
      else
         *AuxMatrix1 = *Mass;
      (*AuxMatrix1) *= 1.0/dt;
      //Laplace - Diffusion
      *AuxMatrix1   += *Laplace;
      //
      int idx = idx1;
      if(pcs->continuum==1)
         idx = idxp21;

      //--------------------------------------------------------------
      //1. Error epsilon
      for (i=0;i<nnodes; i++)
      {
         NodalVal0[i] = fabs(pcs->GetNodeValue(nodes[i],idx)-x_n1[nodes[i]+NodeShift[pcs->continuum]]);
         NodalVal[i] = 0.0;
      }
      if(PcsType==V)                              //
      {
         dof_n = 2;
         //
         // _new for(i=0; i<pcs->pcs_number_of_primary_nvals; i++)
         // _new NodeShift[i] = i*pcs->m_msh->GetNodesNumber(false);
         //
         for (i=0;i<nnodes; i++)
         {
            NodalVal0[i+nnodes] = fabs(pcs->GetNodeValue(nodes[i],idxp21)-x_n1[nodes[i]+NodeShift[1]]);
            NodalVal[i+nnodes] = 0.0;
         }
      }
      else if(PcsType==P)
      {
         dof_n = 2;
         //
         // _new for(i=0; i<pcs->pcs_number_of_primary_nvals; i++)
         // _new NodeShift[i] = i*pcs->m_msh->GetNodesNumber(false);
         //
         for (i=0;i<nnodes; i++)
         {
            NodalVal0[i+nnodes] = fabs(pcs->GetNodeValue(nodes[i],idxSn1)-x_n1[nodes[i]+NodeShift[1]]);
            NodalVal[i+nnodes] = 0.0;
         }
      }
      //
      AuxMatrix1->multi(NodalVal0, NodalVal);

      // Error epsilon
      for (i=0;i<nnodes*dof_n; i++)
         err_norm0 += NodalVal0[i]*NodalVal[i];
      //
      //--------------------------------------------------------------
      //2. Error e_n
      for (i=0;i<nnodes; i++)
      {
         NodalVal0[i] = atol+rtol*max(fabs(pcs->GetNodeValue(nodes[i],idx)),fabs(x_n1[nodes[i]+NodeShift[pcs->continuum]]));
         NodalVal[i] = 0.0;
      }
      if(PcsType==V)                              //
      {
         for (i=0;i<nnodes; i++)
         {
            NodalVal0[i+nnodes] = atol+rtol*max(fabs(pcs->GetNodeValue(nodes[i],idxp21)), fabs(x_n1[nodes[i]+NodeShift[1]]));
            NodalVal[i+nnodes] = 0.0;
         }
      }
      else if(PcsType==P)                         // 03.2009 PCH
      {
         for (i=0;i<nnodes; i++)
         {
            NodalVal0[i+nnodes] = atol+rtol*max(fabs(pcs->GetNodeValue(nodes[i],idxSn1)), fabs(x_n1[nodes[i]+NodeShift[1]]));
            NodalVal[i+nnodes] = 0.0;
         }
      }
      //
      //
      AuxMatrix1->multi(NodalVal0, NodalVal);

      // Error epsilon
      for (i=0;i<nnodes*dof_n; i++)
         err_normn += NodalVal0[i]*NodalVal[i];
      //
      //
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Calculate the energy norm error
   Programing:
   25.09.2008 WW Implementation
   last modification:
   **************************************************************************/
   void CFiniteElementStd::CalcEnergyNorm_Dual(const double *x_n1, double &err_norm0,
      double &err_normn)
   {
      double rtol, atol;
      //----------------------------------------------------------------------
      //
      //
      rtol = pcs->Tim->GetRTol();
      atol = pcs->Tim->GetATol();
      //
      int i,j;
      int gp_r=0, gp_s=0, gp_t=0;
      double W, fkt,mat_fac = 0.;

      //Inintialize
      //-------------------------- WW
      W = pcs->continuum_vector[pcs->GetContinnumType()];
      //
      for(i=0;i<nnodes;i++)
      {
                                                  // Pressure 1
         NodalVal3[i] = pcs->GetNodeValue(nodes[i], idx1);
                                                  // Pressure 2
         NodalVal4[i] = pcs->GetNodeValue(nodes[i], idxp21);
      }
      (*Advection) = 0.0;
      //---------------------------------------------------------
      for (gp = 0; gp < nGaussPoints; gp++)
      {
         //---------------------------------------------------------
         //  Get local coordinates and weights
         //  Compute Jacobian matrix and its determination
         //---------------------------------------------------------
         fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
         mat_fac = CalcCoefDualTransfer();
         mat_fac *= fkt;
         // Material
         ComputeShapefct(1);                      // Linear interpolation function
         // Calculate mass matrix
         for (i = 0; i < nnodes; i++)
         {
            for (j = 0; j < nnodes; j++)
               (*Advection)(i,j) += mat_fac*shapefct[i]*shapefct[j];
         }
      }
      // Add local matrix to global matrix
      long cshift = pcs->m_msh->GetNodesNumber(false);
      //
      double fm = 1.0/W;
      double ff = 1.0/(1.0-W);
      if(MediaProp->transfer_coefficient<0.0)     // for LBNL
         ff = 1.0;
      //
      //--------------------------------------------------------------
      //1. Error epsilon
      for (i=0;i<nnodes; i++)
      {
         NodalVal0[i] = fabs(NodalVal3[i]-x_n1[nodes[i]])
            -fabs(NodalVal4[i]-x_n1[nodes[i]+cshift]);
         NodalVal[i] = 0.0;
      }
      //
      //
      AuxMatrix1->multi(NodalVal0, NodalVal);

      // Error epsilon
      for (i=0;i<nnodes; i++)
         err_norm0 += (fm*(NodalVal3[i]-x_n1[nodes[i]])-
            ff*(NodalVal4[i]-x_n1[nodes[i]+cshift]))
            *NodalVal[i];
      //
      //--------------------------------------------------------------
      //2. Error e_n
      for (i=0;i<nnodes; i++)
      {
         NodalVal0[i] = max(NodalVal3[i],x_n1[nodes[i]])
            -max(NodalVal4[i],x_n1[nodes[i]+cshift]);
         NodalVal[i] = 0.0;
      }
      //
      AuxMatrix1->multi(NodalVal0, NodalVal);
      for (i=0;i<nnodes; i++)
         err_normn += (fm*(atol+rtol*max(NodalVal3[i],x_n1[nodes[i]]))-
            ff*(atol+rtol*max(NodalVal4[i],x_n1[nodes[i]+cshift])))
            *NodalVal[i];
      //
      //
   }
   /**************************************************************************
   PCSLib-Method:
   02/2009 PCH Implementation
   **************************************************************************/
   void CFiniteElementStd::PrintTheSetOfElementMatrices(std::string mark)
   {
      // Output matrices
      if(pcs->Write_Matrix)
      {
         (*pcs->matrix_file) << "### Mark: " << mark << endl;

         (*pcs->matrix_file) << "### Element: " << Index << endl;
         (*pcs->matrix_file) << "---Mass matrix: " << endl;
         if(Mass)
            Mass->Write(*pcs->matrix_file);
         else if(Mass2)
            Mass2->Write(*pcs->matrix_file);
         (*pcs->matrix_file) << "---Laplacian matrix: " << endl;
         Laplace->Write(*pcs->matrix_file);

         (*pcs->matrix_file) << "---AuxMatrix1 matrix: " << endl;
         AuxMatrix1->Write(*pcs->matrix_file);    // PCH for debug
         if(Advection)
         {
                                                  //CMCD
            (*pcs->matrix_file) << "---Advective matrix: " << endl;
            Advection->Write(*pcs->matrix_file);
         }
         if(StrainCoupling)
         {
            (*pcs->matrix_file) << "---Strain couping matrix: " << endl;
            StrainCoupling->Write(*pcs->matrix_file);
         }
         (*pcs->matrix_file) << "---RHS: " <<endl;
         RHS->Write(*pcs->matrix_file);
         (*pcs->matrix_file) <<endl;
         (*pcs->matrix_file) << "Stiffness: " <<endl;
         StiffMatrix->Write(*pcs->matrix_file);
         (*pcs->matrix_file) <<endl;
      }
   }

}                                                 // end namespace


//////////////////////////////////////////////////////////////////////////

using FiniteElement::ElementValue;
vector<ElementValue*> ele_gp_value;
