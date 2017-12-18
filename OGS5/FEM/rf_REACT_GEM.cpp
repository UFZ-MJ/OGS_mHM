/**
 * rf_REACT_GEM.cpp
 * Haibing Shao 25.03.08
 * haibing.shao@ufz.de
 * GEM Reaction Package
 * based on the PSI node-GEM source code
 * using the node-GEM code from Paul Sherrer Institute (PSI)
 * current maintainer: Georg Kosakowski
 * georg.kosakowski@psi.ch
 * last changes: 01 April 2010 start of doxygen documentation
 *	       03 April 2010 start of cleaning code and reducing memory consumption/communication for MPI!
 *              05.10.2010  extend doxygen documentation, kintetics with solid solutions and richards flow coupling
 *
 * Code description
 *
 * The files rf_REACT_GEM.cpp and rf_REACT_GEM.h contain the main core modules of the OpenGeosys - GEMIPM2K coupling.
 * GEMIPM2K is the calculation kernel of GEMS-PSI (http://gems.web.psi.ch). GEMS-PSI executables for various platforms are freely
 * availabe for download.
 * The the kernel GEMIPM2K source code is available on request.
 * GEMIPM2K is currently coupled in a non-iterave sequential way to groundwater flow and transport. The coupling to the Richards flow
 * module is under development.
 *
*/

// There is a name conflict between stdio.h and the MPI C++ binding
// with respect to the names SEEK_SET, SEEK_CUR, and SEEK_END.  MPI
// wants these in the MPI namespace, but stdio.h will #define these
// to integer values.  #undef'ing these can cause obscure problems
// with other include files (such as iostream), so we instead use
// #error to indicate a fatal error.  Users can either #undef
// the names before including mpi.h or include mpi.h *before* stdio.h
// or iostream.
#ifdef USE_MPI_GEMS
#include "mpi.h"//Parallel Computing Support
#include "par_ddc.h"
// HS 07.01.2008: Comment the following 2 lines on LiClus.
// int size;
// int myrank;
#endif
#include "Configure.h"

#include "rf_REACT_GEM.h"
#include "rf_pcs.h"
#include "rfmat_cp.h"
#include "msh_node.h"
#include "msh_elem.h"
// -----------------------
#include "files0.h"
//Headers for shuffling
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>

#include <cstring>
//
#ifdef _WIN32
#include "direct.h" // on win32 and win64 platform
#else
#include "unistd.h" // on unix/linux platform
#include "stdlib.h"
#endif

#ifdef GEM_REACT
/**
REACT_GEM() is the main constructor. It initializes several variable with default values.
*/
REACT_GEM::REACT_GEM ( void )
{
   m_Node = new TNode();
   REACT_GEM::dch_input_file_path = "calcite-dch.dat";
   REACT_GEM::ipm_input_file_path = "calcite-ipm.dat";
   REACT_GEM::dbr_input_file_path = "calcite-dbr-0-0000.dat";
   REACT_GEM::dbr_bc_input_file_path = "calcite-dbr-0-0001.dat";
   REACT_GEM::init_input_file_path= "calcite-init.dat";
   REACT_GEM::init_input_file_sig = "-init.dat";

   nIC = 0;
   nDC = 0;
   nPH = 0;
   nPS = 0;
   nNodes = 0;
   nElems = 0;
   idx_water = -1;
   idx_oxygen = -1;
   idx_hydrogen = -1;

   initialized_flag = 0;
   heatflag = 0;
   flowflag = 0;
   flag_node_element_based = 0;                   //0-node based; 1-elem based;
   flag_porosity_change = 1   ;                   //0-not coupled;1=coupled;
   min_possible_porosity=1.e-4;                   // minimum porostiy in case of changing porosity: avoid zero porosity
   max_possible_porosity=0.9999;                  // max porosity
   flag_coupling_hydrology = 1;                   //0-not coupled;1=coupled;
   flag_permeability_porosity=0;                  //0-no coupling; 1-Kozeny-Carman; 2-Kozeny-Carman normalized;
   m_gem_temperature=298.15;                      //default gem temperature is 25°C, in Kelvin!
   m_gem_pressure=1.0e+5;                         // default pressure 1 bar
   flag_iterative_scheme = 0;                     //0-not iteration;1=iteration;
   // flag for different iterative scheme
   // 0 - sequential non-iterative scheme
   // 1 - standard iterative scheme
   // 2 - symetric iterative scheme
   // 3 - strang splitting scheme
   flag_transport_b=0;                            // default is transport full speciation; 1: transport only dissolved components of b vector
   gem_mass_scale=1.0e-0;                         // GEMS default mass scaling parameter
   flag_gem_smart=0;                              // should be zero if working with kinetics
   m_IterDoneIndexSort = 100;                     // every 100 timesteps we sort NodeIndex for MPI
   m_ShuffleGems = 0;                             // do not randomize GEMS execution
   m_max_failed_nodes = 5; //default number of max allowed nodes to fail
   m_diff_gems=0.0;
   mp_nodeTypes = new long;
   * ( mp_nodeTypes ) = 0;

   m_FluidProp = NULL;
}


REACT_GEM::~REACT_GEM ( void )
{
   if ( initialized_flag > 0 )
   {

      delete [] m_xDC;
      delete [] m_gam;
      delete  [] m_xPH;
      delete []m_aPH;
      delete []m_vPS;
      delete []m_mPS;
      delete []m_bPS;
      delete []m_xPA;
      delete []m_dul;
      delete []m_dll;
      delete []m_uIC;
      delete []m_bIC;
      delete []m_bIC_dummy;
      delete []m_rMB;
      delete []m_xDC_pts;
		delete []m_soluteB_pts;
		delete []m_bIC_pts;
      delete []m_xDC_MT_delta;
      delete []m_xDC_Chem_delta;
      delete []m_NodeHandle;

      delete [] m_NodeStatusCH;
      delete [] m_IterDone;
      delete [] m_IterDoneIndex;
      delete [] m_IterDoneCumulative;
      delete [] m_T;
      delete [] m_P;
      delete [] m_Vs;
      delete [] m_Ms;
      delete [] m_Gs;
      delete [] m_Hs;
      delete [] m_IC;
      delete [] m_pH;
      delete [] m_pe;
      delete [] m_Eh;
      delete [] m_porosity;
		delete [] m_porosity_initial;
      delete [] m_excess_water;
      delete [] m_excess_gas;
      delete [] m_Node_Volume;
      delete [] m_gas_volume;
      delete [] m_fluid_volume;
      delete [] m_soluteB;
      delete [] m_soluteB_buff;
      // delete MPI buffer--------
      delete [] m_NodeHandle_buff;
      delete [] m_NodeStatusCH_buff;
      delete [] m_IterDone_buff;

      delete [] m_Vs_buff;
      delete [] m_Ms_buff;
      delete []  m_Gs_buff;
      delete [] m_Hs_buff;
      delete [] m_IC_buff;
      delete [] m_pH_buff;
      delete [] m_pe_buff;
      delete [] m_Eh_buff;

      delete [] m_xDC_buff;
      delete [] m_xPH_buff;
      delete [] m_xPA_buff,delete [] m_excess_water_buff;
      delete [] m_excess_gas_buff;
      delete [] m_porosity_buff;
      delete [] m_boundary;
      delete [] m_gas_volume_buff;
      delete [] m_fluid_volume_buff;
      delete [] m_dul_buff;
      delete [] m_dll_buff;
      delete [] m_xDC_pts_buff ;
      delete [] m_xDC_MT_delta_buff ;
      delete [] m_xDC_Chem_delta_buff;
      delete [] m_aPH_buff;
      delete [] m_bIC_buff;
      delete [] m_bIC_dummy_buff;
      delete [] m_porosity_Elem_buff;
      delete [] m_porosity_Elem;
      // -------------------------
      delete mp_nodeTypes;

      delete [] mol_phase;
      delete [] omega_phase;
      delete [] omega_components;

      delete [] dmdt;
      delete [] omega_phase_buff;                 // this we need for kinetics
      delete [] mol_phase_buff;                   // this we need for kinetics
      delete [] omega_components_buff;            // this we need for kinetics

      delete [] dmdt_buff;

      m_FluidProp = NULL;
      m_pcs = NULL;
      m_flow_pcs = NULL;
      m_kin.clear();
   }
   delete  m_Node;
}


/**
short REACT_GEM::Init_Nodes ( string Project_path )

Initialization of the GEM TNode Class

Here we read the files needed as input for initializing GEMIPM2K
The easiest way to prepare them is to use GEMS-PSI code (GEM2MT module)

**/
short REACT_GEM::Init_Nodes ( string Project_path )
{
   long ii=0, i=0, in=0;

   // Creating TNode structure accessible trough node pointer
   // Here we read the files needed as input for initializing GEMIPM2K
   // The easiest way to prepare them is to use GEMS-PSI code (GEM2MT module)
   if ( Load_Init_File ( Project_path ) )
   {
      // The init file is successfully loaded
      // Getting direct access to DataCH structure in GEMIPM2K memory
      dCH = m_Node->pCSD();
      if ( !dCH )
         return 3;

      // Getting direct access to work node DATABR structure which
      // exchanges data between GEMIPM and FMT parts
      dBR = m_Node->pCNode();
      if ( !dBR )
         return 4;

      // Extracting data bridge array sizes
      nIC = dCH->nICb;                            //Num of Independent components
      nDC = dCH->nDCb;                            //Num of Chemical species in the reactive part
      nPH = dCH->nPHb;                            //Num of Phases
      nPS = dCH->nPSb;                            //Num of multicomponent phases; ASSERT(nPS < nPH)

      // get the index of water
      idx_water = FindWater_xDC();
      // imediately check
      if ( idx_water == -1 ) return 1;
      // get index for H and OH
      idx_oxygen = Findoxygen_bIC();
      // imediately check
      if ( idx_oxygen == -1 ) return 1;
      idx_hydrogen = Findhydrogen_bIC();
      // imediately check
      if ( idx_hydrogen == -1 ) return 1;

      heatflag = GetHeatFlag_MT();                // Get heatflag
      flowflag = GetFlowType_MT();                // Get flow flag

      // get m_flow_pcs already, then check the flag:
      if ( flag_coupling_hydrology == 1 )
      {
         // need to couple to flow process;
         // mark the flag
         m_flow_pcs->flag_couple_GEMS = 1;
      }
      else
      {
         m_flow_pcs->flag_couple_GEMS = 0;        //default is set to not-coupled!
      }

      gem_pressure_flag=0;                        // use not a user defined pressure
      // Get number of Nodes
      nNodes = GetNodeNumber_MT();
      // Get number of Elems
      nElems = GetElemNumber_MT();
      // kg44 03 april 2010 first we take all variables that we need for all nodes!

      // Allocating work memory for FMT part (here only chemical variables)
      m_NodeHandle = new long [nNodes];
      m_NodeHandle_buff = new long [nNodes];

      m_NodeStatusCH = new long [nNodes];
      m_NodeStatusCH_buff = new long [nNodes];

      m_IterDone = new long [nNodes];
      m_IterDone_buff = new long[nNodes];

      m_IterDoneCumulative = new long [nNodes];
      m_IterDoneIndex = new long[nNodes];

      // MPI Buffer Variable---------------

      m_boundary = new int [nNodes];              // this marks boundary nodes with fixed concentrations!?
      m_T  = new double [nNodes];
      m_P  = new double [nNodes];

      m_Vs = new double [nNodes];
      m_Vs_buff = new double[nNodes];

      m_Ms = new double [nNodes];
      m_Ms_buff = new double[nNodes];

      m_Gs = new double [nNodes];
      m_Gs_buff = new double[nNodes];

      m_Hs = new double [nNodes];
      m_Hs_buff = new double[nNodes];

      m_IC = new double [nNodes];
      m_IC_buff = new double[nNodes];

      m_pH = new double [nNodes];
      m_pH_buff = new double[nNodes];

      m_pe = new double [nNodes];
      m_pe_buff = new double[nNodes];

      m_Eh = new double [nNodes];
      m_Eh_buff = new double[nNodes];

      m_porosity     = new double [nNodes];
      m_porosity_buff = new double[nNodes];
		m_porosity_initial     = new double [nNodes];

      m_excess_water = new double [nNodes];
      m_excess_water_buff = new double [nNodes];

      m_excess_gas = new double [nNodes];
      m_excess_gas_buff = new double [nNodes];

      m_Node_Volume  = new double [nNodes];

      m_fluid_volume  = new double [nNodes];
      m_fluid_volume_buff  = new double [nNodes];

      m_gas_volume  = new double [nNodes];
      m_gas_volume_buff  = new double [nNodes];

      m_porosity_Elem = new double [nElems];
      m_porosity_Elem_buff = new double [nElems];

      m_soluteB = new double [nNodes*nIC];
      m_soluteB_buff = new double [nNodes*nIC];

      m_bIC = new double [nNodes*nIC];
      m_bIC_buff = new double [nNodes*nIC];

      m_bIC_dummy = new double [nNodes*nIC];
      m_bIC_dummy_buff = new double [nNodes*nIC];

      m_dul = new double [nNodes*nDC];
      m_dll = new double [nNodes*nDC];
      m_dul_buff = new double [nNodes*nDC];
      m_dll_buff = new double [nNodes*nDC];

      m_xDC = new double [nNodes*nDC];
      m_xDC_buff = new double [nNodes*nDC];

      m_aPH = new double [nNodes*nPH];            // surface area for surface species ..input to GEMS!
      m_aPH_buff = new double [nNodes*nPH];
      m_xPH = new double [nNodes*nPH];            // amount of carrier...used for smart initial aproximation
      m_xPH_buff = new double [nNodes*nPH];
      m_xPA = new double [nNodes*nPS];
      m_xPA_buff = new double [nNodes*nPS];

      // ----------------------------------
      // this is for kinetics
      Kinetic_GEMS m_kin;                         // new kinetic vector

      omega_phase = new double [nNodes*nPH];
      omega_phase_buff = new double [nNodes*nPH];
      mol_phase = new double [nNodes*nPH];
      mol_phase_buff = new double [nNodes*nPH];
      dmdt       = new double [nNodes*nPH];
      dmdt_buff       = new double [nNodes*nPH];
      omega_components = new double [nNodes*nDC];
      omega_components_buff = new double [nNodes*nDC];

      m_xDC_pts = new double [nNodes*nDC];
		m_soluteB_pts = new double [nNodes*nIC];
		m_bIC_pts = new double [nNodes*nIC];
      m_xDC_MT_delta = new double [nNodes*nDC];
      m_xDC_Chem_delta =  new double [nNodes*nDC];

      m_xDC_pts_buff = new double [nNodes*nDC];
      m_xDC_MT_delta_buff = new double [nNodes*nDC];
      m_xDC_Chem_delta_buff =  new double [nNodes*nDC];

      // kg44 03 april 2010 ...from here on, most data is only necessary once (check if this is needed for all nodes!)
      m_rMB = new double [nIC];                   // buffer removed...we do not need the mass balance residuals globally
      m_uIC = new double [nIC];                   // chemical potentials...in current code not needed, but maybe later?

      m_gam = new double [nDC];

      m_vPS = new double [nPS];
      m_mPS = new double [nPS];
      m_bPS = new double [nIC*nPS];

      // ------------------------
      for ( in = 0; in < nNodes ; in++ )
      {
         m_boundary[in]=0;                        // cout << "boundary init ok"<<endl;
         m_NodeHandle[in] = 0;
         m_NodeStatusCH[in] = 0;
         m_IterDone[in] = 0;
         m_IterDoneCumulative[in]=0;
         m_IterDoneIndex[in]=in;                  //initial values order of Nodes
         m_NodeHandle_buff[in] = 0;
         m_NodeStatusCH_buff[in] = 0;
         m_IterDone_buff[in] = 0;

         m_T[in] = 298.15;                        //equivalent to 25°C
         m_P[in] = 1.0e+5;
         m_Vs[in] = 0.0;
         m_Ms[in] = 0.0;
         m_Gs[in] = 0.0;
         m_Hs[in] = 0.0;
         m_IC[in] = 0.0;
         m_pH[in] = 0.0;
         m_pe[in] = 0.0;
         m_Eh[in] = 0.0;
         m_porosity[in]=0.0;
			m_porosity_initial[in]=0.0;
         m_fluid_volume[in]=0.0;
         m_gas_volume[in]=0.0;

         m_Vs_buff[in] = 0.0;
         m_Ms_buff[in] = 0.0;
         m_Gs_buff[in] = 0.0;
         m_Hs_buff[in] = 0.0;
         m_IC_buff[in] = 0.0;
         m_pH_buff[in] = 0.0;
         m_pe_buff[in] = 0.0;
         m_Eh_buff[in] = 0.0;
         m_porosity_buff[in]= 0.0;
         m_fluid_volume_buff[in]=0.0;
         m_gas_volume_buff[in]=0.0;

         m_excess_water[in] = 0.0;
         m_excess_gas[in] = 0.0;

         m_Node_Volume[in]  = REACT_GEM::GetNodeAdjacentVolume ( in );
         m_excess_water_buff[in] = 0.0;
         m_excess_gas_buff[in] = 0.0;

         for ( ii = 0; ii < nIC ; ii++ )
         {
            m_bIC[in*nIC + ii ] = 0.0;
            m_bIC_dummy [ in*nIC + ii ] = 0.0;

            m_soluteB[in*nIC + ii ] = 0.0;
            m_soluteB_buff[ in*nIC + ii ] = 0.0;
				m_soluteB_pts[in*nIC+ii ] = 0.0;
				m_bIC_pts[in*nIC+ii ] = 0.0;
            m_rMB[ii ] = 0.0;
            m_uIC[ii ] = 0.0;

            m_bIC_buff[in*nIC + ii ] = 0.0;
            m_bIC_dummy_buff [ in*nIC + ii ] = 0.0;
         }

         for ( ii = 0; ii < nDC ; ii++ )
         {
            m_xDC[in*nDC+ii ] = 0.0;
            m_gam[ii ] = 0.0;
            m_dul[in*nDC+ii ] = 1.0e+10;          // this should be a large number, because after scaling to 1kg in Gems it should be 1.e+6
            m_dll[in*nDC+ii ] = 0.0;              // zero is ok
            m_xDC_pts[in*nDC+ii ] = 0.0;
            m_xDC_MT_delta[in*nDC+ii ] = 0.0;
            m_xDC_Chem_delta[in*nDC+ii ] = 0.0;

            m_xDC_buff[in*nDC+ii ] = 0.0;
            m_dul_buff[in*nDC+ii ] = 0.0;
            m_dll_buff[in*nDC+ii ] = 0.0;
            m_xDC_pts_buff[in*nDC+ii ] = 0.0;
            m_xDC_MT_delta_buff[in*nDC+ii ] = 0.0;
            m_xDC_Chem_delta_buff[in*nDC+ii ] = 0.0;
            omega_components[in*nDC+ii ] = 0.0;
            omega_components_buff[in*nDC+ii ] = 0.0;
         }

         for ( ii = 0; ii < nPH ; ii++ )
         {
            m_aPH[in*nPH+ii ] = 0.0;
            m_xPH[in*nPH+ii ] = 0.0;

            m_xPH_buff[in*nPH+ii ] = 0.0;
            m_aPH_buff[in*nPH+ii ] = 0.0;

            omega_phase[in*nPH+ii ] = 0.0;
            omega_phase_buff[in*nPH+ii ] = 0.0;
            mol_phase[in*nPH+ii ] = 0.0;
            mol_phase_buff[in*nPH+ii ] = 0.0;

            dmdt[in*nPH+ii ] = 0.0;
            dmdt_buff[in*nPH+ii ] = 0.0;

         }

         for ( ii = 0; ii < nPS ; ii++ )
         {
            m_vPS[ii ] = 0.0;
            m_mPS[ii ] = 0.0;
            m_xPA[in*nPS+ii ] = 0.0;
            m_xPA_buff[in*nPS+ii ] = 0.0;
         }

         for ( ii = 0; ii < nIC ; ii++ )
         {
            for ( int jj = 0; jj < nPS ; jj++ )
            {
               m_bPS[ii*nPS+jj ] = 0.0;
            }
         }
      }

      for ( in = 0; in < nElems ; in++ )
      {
         m_porosity_Elem[in] = 0.0;
         m_porosity_Elem_buff[in] = 0.0;
      }
      nNodes = GetNodeNumber_MT();
      nElems = GetElemNumber_MT();
      // run Gems for all nodes exactly once! ---will be done in the Init-Run again with IC and BC values
#ifdef USE_MPI_GEMS
      // MPI initialization.
      // So here is going to distribute the task.
      MPI_Bcast ( &nNodes, 1, MPI_LONG, 0, MPI_COMM_WORLD );
      // here "myrank" is the index of the CPU Processes, and "mysize" is the number of CPU Processes
      for ( ii = myrank; ii < nNodes ; ii+= mysize )
      {
         i=m_IterDoneIndex[ii];                   // m_IterDoneIndex can be sortet to account for nodes that spend a lot of time for gems
         //		cout << " rank, calc node, ii " << myrank << " " << in <<" " <<ii<< endl;
#else
         for ( ii = 0; ii < nNodes; ii++ )
         {
            i=m_IterDoneIndex[ii];                // m_IterDoneIndex can be sortet to account for nodes that spend a lot of time for gems
            // cout <<"iteration "<< ii << " node " << in <<endl;
#endif
            dBR->NodeStatusCH = NEED_GEM_AIA;
            m_NodeStatusCH[i] = m_Node->GEM_run ( false );

            if ( ! ( m_NodeStatusCH[i] == OK_GEM_AIA || m_NodeStatusCH[i] == OK_GEM_SIA ) )
            {
               dBR->NodeStatusCH = NEED_GEM_AIA;
               m_NodeStatusCH[i] = m_Node->GEM_run ( false );
               if ( m_NodeStatusCH[i] == ERR_GEM_AIA )
               {
                  cout << "critical error: Initial GEMs run after Read GEMS failed at node " << i << endl;
					m_Node->GEM_write_dbr ( "dbr_for_crash_node_init.txt" );
					m_Node->GEM_print_ipm ( "ipm_for_crash_node_init.txt" );
#ifdef USE_MPI_GEMS
                  MPI_Finalize();                 //make sure MPI exits
#endif

                  exit ( 1 );
               }
				else if ( m_NodeStatusCH[i] == BAD_GEM_AIA )
               {
                  cout << "error: Initial GEMs run after Read GEMS gives bad result..proceed in any case " << i << endl;
               }

                                                  // this we need also for restart runs
               REACT_GEM::GetReactInfoFromGEM ( i );
            }
            else                                  // init arrays with value from databridge file ..all nodes should get the same
            {
                                                  // this we need also for restart runs
               REACT_GEM::GetReactInfoFromGEM ( i );

            }
#ifdef USE_MPI_GEMS
            REACT_GEM::CopyToMPIBuffer ( i );     //copy data to MPI buffer in any case
#endif
         }                                        // end for loop for all nodes
#ifdef USE_MPI_GEMS
         // For MPI scheme, gather the data here.
         REACT_GEM::GetGEMResult_MPI();
         REACT_GEM::CleanMPIBuffer();
#endif
         // m_Node->na->GEM_write_dbr ( "dbr_for_crash_node1.txt" );
         // m_Node->na->GEM_print_ipm ( "ipm_for_crash_node1.txt" );

         return 0;                                //successed
      }
      else
      {
#ifdef USE_MPI_GEMS
         MPI_Finalize();                          //make sure MPI exits
#endif

         cout << "Error loading initial files to GEMS" <<endl;
         exit ( 1 );
         return 1;
      }
   }

   short REACT_GEM::Init_RUN()
   {
      nNodes = GetNodeNumber_MT();
      nElems = GetElemNumber_MT();
      long ii=0,in = 0,i=0,j;
      //        CompProperties *m_cp = NULL;
      CRFProcess* this_pcs;

      double BCValue = 0.0;

      string cstr;

      for ( ii=0;ii<(int) m_kin.size();ii++ )
      {
         m_kin[ii].phase_number=-1;
         m_kin[ii].dc_counter=0;                  //this is the starting dc
         // phase numbers are not yet set...do it!
         for ( j=0;j<nPH;j++ )
         {
            cstr.assign ( dCH->PHNL[j] );
            cout << m_kin[ii].phase_name  << " gems phase: " << cstr << " "<< " component start number " <<  m_kin[ii].dc_counter << endl;
            if ( m_kin[ii].phase_name == cstr )
            {
               m_kin[ii].phase_number=j;
               break;
            }
                                                  //add the number of components to the counter..after test for phasename...as this counter gives the starting position of the dependent components
            m_kin[ii].dc_counter+=dCH->nDCinPH[j] ;
         }

         if ( m_kin[ii].phase_number < 0 || m_kin[ii].phase_number >= nPH )
         {
            cout << " GEMS: Error in Phase kinetics..check input " << endl;
#ifdef USE_MPI_GEMS
            MPI_Finalize();                       //make sure MPI exits
#endif

            exit ( 1 );
         }
         else
         {
            cout << "GEM Kinetics phase number:  " <<m_kin[ii].phase_name<< " in phase number " << m_kin[ii].phase_number<<endl;
         }
      }
      // Marking species with fixed concentrations (in boundary nodes)
      // first we check the boundary nodes with fixed concentrations
      // this is adopted from rf_REACT_BRNS...we only look for the first species, as with GEMS we should define boundary conditions for ALL species
      for ( i=0 ; i < nNodes ; i++ )
      {
         this_pcs = NULL;
         //			m_cp = cp_vec[0];

         // Get the pointer to the proper PCS.
         this_pcs = PCSGet ( "MASS_TRANSPORT" );
         if ( this_pcs )
         {
            // BC printing
            if ( IsThisPointBCIfYesStoreValue ( i, this_pcs, BCValue ) )
            {
               // If this node is on the fixed boudnary for this component
               cout << "Node " << i <<", Comp " << this_pcs->pcs_primary_function_name[0] << " ,Value " << BCValue << " is boundary node" <<endl;
               m_boundary[i] = 1;
            }
            else
            {
               // If this node is NOT on the fixed boudnary for this component
               m_boundary[i] = 0;
            }
		}
		else   // not getting the pointer to the proper PCS.
         {
            cout << this_pcs->pcs_primary_function_name[0]<<"!!! In InitGEMS, can not find corresponding PCS for checking boundary conditions! "<<endl;
            abort();
         }

      }

      // we have to make sure saturation is properly defined..otherwise we can not calculate concentrations properly
      if ( flowflag == 3 )                        //this is only for Richards flow
      {
         cout << "GEM-INIT: " ;
         if ( m_flow_pcs->saturation_switch == true )
         {
            cout << "CalcSaturationRichards " << endl;
                                                  //is true here correct?
            m_flow_pcs->CalcSaturationRichards ( 1, true );
         }                                        // JOD
         else
         {
                                                  //WW
            m_flow_pcs->CalcSecondaryVariablesUnsaturatedFlow ( true );
            cout << " CalcSecondaryVariablesUnsaturatedFlow" << endl;
         }
      }

      //if (!m_flow_pcs->GetRestartFlag() >=2) {
      // now we do the calculation!
      // for GEM_REACT we also need internal information on porosity!....

      if ( ( m_flow_pcs->GetRestartFlag() >=2 ) )
      {
         if ( !ReadReloadGem() ) abort();
      }
	else
	{
		cout << "Attentione GEMS users: Initial kinetics calculated without restart! This probably kills kinetics, as phases in m_xDC are not yet properly initialized!" << endl;
		cout << "No upper or lower constrains set during equilibration!...If your setup requires constrains, please contact georg.kosakowski@psi.ch" << endl;
	}

      GetInitialReactInfoFromMassTransport ( 1 ); //get the initial values from MD - last time step!

#ifdef USE_MPI_GEMS
      // MPI initialization.
      // So here is going to distribute the task.
      MPI_Bcast ( &nNodes, 1, MPI_LONG, 0, MPI_COMM_WORLD );
      // here "myrank" is the index of the CPU Processes, and "mysize" is the number of CPU Processes
      for ( ii = myrank; ii < nNodes ; ii+= mysize )
      {
         in=m_IterDoneIndex[ii];                  // m_IterDoneIndex can be sortet to account for nodes that spend a lot of time for gems
         //		cout << " rank, calc node, ii " << myrank << " " << in <<" " <<ii<< endl;
#else
         for ( ii = 0; ii < nNodes; ii++ )
         {
            in=m_IterDoneIndex[ii];               // m_IterDoneIndex can be sortet to account for nodes that spend a lot of time for gems
            // cout <<"iteration "<< ii << " node " << in <<endl;
#endif

            // only when porosity is coupled; kg44...no..should be done always! as porosity might be spatially varying
            // if ( flag_porosity_change > 0 )
            //{
            // if ( m_porosity[in]<min_possible_porosity ) m_porosity[in]= min_possible_porosity;
            // if ( m_porosity[in]>max_possible_porosity ) m_porosity[in]= max_possible_porosity;
            // here get mass of
            // Convert to concentration NOT NECESSARY FOR FIRST RUN AS ICs should be given in MOLS...but necessary for reload case!!
            // cout << m_porosity[in] << endl;
            // cout << " restart flag " <<  m_flow_pcs->GetRestartFlag() << endl;
            if ( ( m_flow_pcs->GetRestartFlag() >=2 ) ) REACT_GEM::ConcentrationToMass ( in ,1 );
            //}
            // now we calculate kinetic constraints for GEMS!
            REACT_GEM::CalcLimitsInitial ( in );  //switched of for debug...
            //Get data
            REACT_GEM::SetReactInfoBackGEM ( in );// this is necessary, otherwise the correct data is not available

		// m_Node->na->GEM_write_dbr ( "dbr_for_crash_node0.txt" );
            // m_Node->na->GEM_print_ipm ( "ipm_for_crash_node0.txt" );
            // Order GEM to run
            dBR->NodeStatusCH = NEED_GEM_AIA;

            m_NodeStatusCH[in] = m_Node->GEM_run ( false );

            if ( ! ( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) )
            {
               cout << "Initial GEMs run first pass failed at node " << in ;
               dBR->NodeStatusCH = NEED_GEM_AIA;
               m_NodeStatusCH[in] = m_Node->GEM_run ( false );
               if ( ! ( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) )
               {
                  cout << " Error: Init Loop failed when running GEM on Node #" << in << "." << endl;
                  cout << "Returned Error Code: " << m_NodeStatusCH[in] << endl;
				m_Node->GEM_write_dbr ( "dbr_for_crash_node.txt" );
#ifdef USE_MPI_GEMS
                  MPI_Finalize();                 //make sure MPI exits
#endif

                  exit ( 1 );
               }
               else
               {
                  cout << " sucess with second try.... "<<  endl;

               }

            }                                     // end loop if initial gems run fails

            //Get data
            REACT_GEM::GetReactInfoFromGEM ( in );//

            // scale data so that second pass gives the normalized volume of 1m^3
            for ( j=0 ; j < nDC; j++ )
            {
               m_xDC[in*nDC+j] /= m_Vs[in] ;
               if ( m_xDC[in*nDC+j] < 0.0 )
                  m_xDC[in*nDC+j]=0.0;

            }
		for ( j=0;j<nIC;j++ )
		{
			m_bIC[in*nIC+j]/= m_Vs[in]; //This is then for b vector
		}

            //Get data
            REACT_GEM::SetReactInfoBackGEM ( in );// this is necessary, otherwise the correct data is not available
            // m_Node->na->GEM_write_dbr ( "dbr_for_crash_node0.txt" );
            // m_Node->na->GEM_print_ipm ( "ipm_for_crash_node0.txt" );
            // Order GEM to run
            dBR->NodeStatusCH = NEED_GEM_AIA;

            m_NodeStatusCH[in] = m_Node->GEM_run ( false );

            if ( ! ( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) )
            {
               cout << "Initial GEMs run second pass failed at node " << in ;

               dBR->NodeStatusCH = NEED_GEM_AIA;
               m_NodeStatusCH[in] = m_Node->GEM_run (false );
               if ( ! ( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) )
               {
                  cout << " Error: Init Loop failed when running GEM on Node #" << in << "." << endl;
                  cout << "Returned Error Code: " << m_NodeStatusCH[in] << endl;
				m_Node->GEM_write_dbr ( "dbr_for_crash_node.txt" );
#ifdef USE_MPI_GEMS
                  MPI_Finalize();                 //make sure MPI exits
#endif

                  exit ( 1 );
               }
               else
               {
				cout << " sucess with second try.... "<<  endl;

               }
            }                                     // end loop if initial gems run fails

            REACT_GEM::GetReactInfoFromGEM ( in );// this we need also for restart runs



		// calculate the chemical porosity
		if ( m_flow_pcs->GetRestartFlag() <2 ) REACT_GEM::CalcPorosity ( in ); //during init it should be always done, except for restart !!!

		REACT_GEM::CalcReactionRate ( in, m_T[in] ); //moved it after porosity calculation because of chrunchflow kinetics (model 4)!
            // Convert to concentration: should be done always...
		REACT_GEM::MassToConcentration ( in,0 /*no failed node*/ ); // concentrations are now based on the fluid-gas volumes...

#ifdef USE_MPI_GEMS
            REACT_GEM::CopyToMPIBuffer ( in );    //copy data to MPI buffer in any case
#endif
         }                                        // end for loop for all nodes
#ifdef USE_MPI_GEMS
         // For MPI scheme, gather the data here.
         REACT_GEM::GetGEMResult_MPI();
         REACT_GEM::CleanMPIBuffer();
#endif



	for ( in = 0; in < nNodes; in++ )
	{
		if ( flag_transport_b==0 )
		{
			// we should also push back the initial system including boundary nodes to avoid inconsistencies

			REACT_GEM::SetDCValue_MT ( in , 0, & ( m_xDC[in*nDC] ) ); // old timestep
			REACT_GEM::SetDCValue_MT ( in , 1, & ( m_xDC[in*nDC] ) ); // new timestep
		}
		else
		{

		}
		REACT_GEM::SetBValue_MT ( in , 0, & ( m_soluteB[in*nIC] ) ); // old timestep
		REACT_GEM::SetBValue_MT ( in , 1, & ( m_soluteB[in*nIC] ) ); // new timestep

	}

	// when switch is on;
	ConvPorosityNodeValue2Elem ( 0 ); // old timestep: update element porosity and push back values
	ConvPorosityNodeValue2Elem ( 1 ); // new timestep: update element porosity and push back values
	CopyCurBPre();

         cout << "Initial Running GEM  to get the correct porosities successful. "  << endl;

         return 0;
      }

      string REACT_GEM::Get_Init_File_Path ( void )
      {
         return init_input_file_path;
      }

      string REACT_GEM::Get_IPM_File_Path ( void )
      {
         return ipm_input_file_path;
      }

      string REACT_GEM::Get_DBR_File_Path ( void )
      {
         return dbr_input_file_path;
      }

      string REACT_GEM::Get_DCH_File_Path ( void )
      {
         return dch_input_file_path;
      }

      int REACT_GEM::Set_IPM_FILE_PATH ( string m_path )
      {
         REACT_GEM::ipm_input_file_path = m_path;
         return 0;
      }

      int REACT_GEM::Set_DBR_FILE_PATH ( string m_path )
      {
         REACT_GEM::dbr_input_file_path = m_path;
         return 0;
      }

      int REACT_GEM::Set_DCH_FILE_PATH ( string m_path )
      {
         REACT_GEM::dch_input_file_path = m_path;
         return 0;
      }

      int REACT_GEM::Set_Init_File_Path ( string m_path )
      {
         REACT_GEM::init_input_file_path = m_path;
         return 0;
      }
      bool REACT_GEM::Load_Init_File ( string m_Project_path )
      {
         string init_path;
         char *buffer=NULL;

         init_path = m_Project_path.append ( REACT_GEM::init_input_file_path );

#ifdef _WIN32
                                                  // keep this on windows
         if ( init_path.rfind ( "\\" ) == string::npos )
#else
                                                  // keep this on linux
            if ( init_path.rfind ( "/" ) == string::npos )
#endif
         {
#ifdef _WIN32
            if ( ( buffer = _getcwd ( NULL, 0 ) ) == NULL )
#else
               if ( ( buffer = getcwd ( NULL, 0 ) ) == NULL )
#endif
                  perror ( "_getcwd error" );
            else
            {
#ifdef _WIN32
               init_path.insert ( 0, "\\" );      // keep this on window
#else
               init_path.insert ( 0, "/" );       // keep this on linux
#endif
               init_path.insert ( 0, buffer );
            }
         }

         if ( buffer ) free ( buffer );

         if ( m_Node->GEM_init ( init_path.c_str() , mp_nodeTypes , false ) )
         {
            return 0;                             // error occured during reading the files
         }
         else
         {
            return 1;                             // read init file successed
         }
      }

      short REACT_GEM::GetInitialReactInfoFromMassTransport ( int timelevel )
      {
         heatflag = GetHeatFlag_MT();
         flowflag = GetFlowType_MT();
         REACT_GEM::nNodes = GetNodeNumber_MT();

         for ( long node_i=0; node_i < nNodes ; node_i++ )
         {
            //get temperature from MT
            m_T[node_i] = REACT_GEM::GetTempValue_MT ( node_i, timelevel );

            //get pressure from MT
            m_P[node_i] = REACT_GEM::GetPressureValue_MT ( node_i, timelevel );
            //get Independent and dependent Component value from MT
            if ( flag_transport_b == 0 ) REACT_GEM::GetDCValue_MT ( node_i, timelevel, m_xDC+node_i*nDC, m_xDC_pts+node_i*nDC, m_xDC_MT_delta+node_i*nDC );
		else if ( ( flag_transport_b ==1 ) && ( m_flow_pcs->GetRestartFlag() <2 ) ) REACT_GEM::GetBValue_MT ( node_i, timelevel, m_bIC+node_i*nIC ); //do this not for restart...overwrites values!!!
            // else {cout << "RESTART. no initial values for GEMS initialization read from mass transport!" << endl;  }
            // Convert to mole values
            // if(conv_concentration == 1)	REACT_GEM::ConcentrationToMass( node_i );

            // Setting Solid Phase Component // HS: Solid does not move.
            // REACT_GEM::GetSoComponentValue_MT(node_i, timelevel, m_xPH+node_i*nPH );

         }

         return 0;
      }

      short REACT_GEM::GetReactInfoFromMassTransport ( int timelevel )
      {
         heatflag = GetHeatFlag_MT();
         flowflag = GetFlowType_MT();
         REACT_GEM::nNodes = GetNodeNumber_MT();

         for ( long node_i=0; node_i < nNodes ; node_i++ )
         {
            //get temperature from MT
            m_T[node_i] = REACT_GEM::GetTempValue_MT ( node_i, timelevel );

            //get pressure from MT
            m_P[node_i] = REACT_GEM::GetPressureValue_MT ( node_i, timelevel );
            //get Independent and dependent Component value from MT
            if ( flag_transport_b == 0 ) REACT_GEM::GetDCValue_MT ( node_i, timelevel, m_xDC+node_i*nDC, m_xDC_pts+node_i*nDC, m_xDC_MT_delta+node_i*nDC );
            else REACT_GEM::GetBValue_MT ( node_i, timelevel, m_soluteB+node_i*nIC );
            // Convert to mole values
            // if(conv_concentration == 1)	REACT_GEM::ConcentrationToMass( node_i );

            // Setting Solid Phase Component // HS: Solid does not move.
            // REACT_GEM::GetSoComponentValue_MT(node_i, timelevel, m_xPH+node_i*nPH );

         }

         return 0;
      }

      short REACT_GEM::SetReactInfoBackMassTransport ( int timelevel )
      {

         for ( long in=0; in < nNodes ; in++ )
         {

            // Setting Temperature // disabled by HS. temperature is NOT the output from chemistry.
            // REACT_GEM::SetTempValue_MT(in,timelevel,m_T[in]);

            // Setting Pressure // disabled by HS. pressure is NOT the output from chemistry.
            // REACT_GEM::SetPressureValue_MT(in,timelevel,m_P[in]);

            // if (m_pcs->m_msh->nod_vector[in]->onBoundary() == false) {
            // Setting Independent Component
            if ( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA )
            {
			if ( flag_transport_b == 0 ) REACT_GEM::SetDCValue_MT ( in , timelevel , & ( m_xDC[in*nDC] ) );
			else REACT_GEM::SetBValue_MT ( in , timelevel , & ( m_soluteB[in*nIC] ) );
            }

		// Set the extra water as source/sink term; not for boundary nodes

		if ( flag_coupling_hydrology>0 && !m_boundary[in] ) REACT_GEM::SetSourceSink_MT ( in, dt /*in sec*/ );
         }
#ifdef USE_MPI_GEMS
	if ( flag_coupling_hydrology>0 )
	{
		m_flow_pcs->SetSTWaterGemSubDomain ( myrank ); // necessary for domain decomposition
	}
#endif
	if ( flag_porosity_change>0 ) ConvPorosityNodeValue2Elem ( timelevel ); // new timestep :update element porosity and push back values
	return 0;
}

      void REACT_GEM::GetReactInfoFromGEM ( long in )
      {
         m_Node->GEM_to_MT ( m_NodeHandle[in], m_NodeStatusCH[in], m_IterDone[in],
            m_Vs[in], m_Ms[in], m_Gs[in], m_Hs[in], m_IC[in], m_pH[in], m_pe[in], m_Eh[in],
            m_rMB, m_uIC, m_xDC+in*nDC, m_gam, m_xPH+in*nPH, m_vPS, m_mPS,
            m_bPS, m_xPA+in*nPS, m_aPH+in*nPH );

         return;
      }

      void REACT_GEM::SetReactInfoBackGEM ( long in )
      {
         long i;
         // Setting input data for GEMIPM

         // set charge to zero

         //	cout << *m_xDC+in*nDC << "   " << *m_bIC_dummy+in*nIC << endl;

         if ( flag_transport_b == 0 )
         {
            // Using the overloaded version of GEM_from_MT() to load the data	// HS 10.07.2007
            for ( i=0;i<nIC;i++ )
            {
               m_bIC_dummy[in*nIC+i]=0.0;         //old B vector: this should be always zero
            }
            m_Node->GEM_from_MT ( m_NodeHandle[in], m_NodeStatusCH[in],
               m_T[in], m_P[in], m_Vs[in], m_Ms[in],
               m_bIC_dummy+in*nIC, m_dul+in*nDC, m_dll+in*nDC, m_aPH+in*nPH  ,m_xDC+in*nDC );
	}
	else   //here we insert the actual B vector
         {
            m_Node->GEM_from_MT ( m_NodeHandle[in], m_NodeStatusCH[in],
               m_T[in], m_P[in], m_Vs[in], m_Ms[in],
               m_bIC+in*nIC, m_dul+in*nDC, m_dll+in*nDC, m_aPH+in*nPH );
         }
         //	cout << m_xDC+in*nDC << endl;
         // set charge to zero
         m_Node->pCNode()->bIC[nIC-1]=0.0;

         return;
      }

      short REACT_GEM::Run_MainLoop ( )
      {
         nNodes = GetNodeNumber_MT();
         nElems = GetElemNumber_MT();
         long /*i,j,ii,*/in,ii,idummy,node_fail=0, repeated_fail=0;
         double oldvolume;

#ifdef USE_MPI_GEMS
         // MPI initialization.
         // So here is going to distribute the task.
         MPI_Bcast ( &nNodes, 1, MPI_LONG, 0, MPI_COMM_WORLD );
         // here "myrank" is the index of the CPU Processes, and "mysize" is the number of CPU Processes
         repeated_fail=0;                         // make sure this is zero at the beginning ....if this is greater 3 program stops because number of failed nodes per CPU is exeeded!
         for ( ii = myrank; ii < nNodes ; ii+= mysize )
         {
            in=m_IterDoneIndex[ii];               // m_IterDoneIndex can be sortet to account for nodes that spend a lot of time for gems
            //		cout << " rank, calc node, ii " << myrank << " " << in <<" " <<ii<< endl;
#else
                                                  // randomize order of Gems executions !
            if ( m_ShuffleGems ) ShuffleIterations ( m_IterDoneIndex, nNodes );
            for ( ii = 0; ii < nNodes; ii++ )
            {
               in=m_IterDoneIndex[ii];            // m_IterDoneIndex can be sortet to account for nodes that spend a lot of time for gems
               // cout <<"iteration "<< ii << " node " << in <<endl;
#endif
		// if ( (CalcSoluteBDelta(in) > m_diff_gems)) { cout << "DEBUG GEMS: node " <<in << " " <<   CalcSoluteBDelta(in) << endl;}
		if ( m_boundary[in] || ( CalcSoluteBDelta ( in ) < m_diff_gems ) )   // do this only without calculation if  on a boundary or differences very small!!!
		{
			// cout << "DEBUG GEMS: node " <<in << " " <<   CalcSoluteBDelta(in) << " " << m_diff_gems << endl;
#ifdef USE_MPI_GEMS
			REACT_GEM::CopyToMPIBuffer ( in ); // copy old values to buffer..otherwise we loose them
#endif
		}
		else   // calculate if not on a boundary or transport causes differences in solutes
		{
			//if(in == 100) { m_Node->na->GEM_write_dbr("dbr_for_crash_node_1.txt");  m_Node->na->GEM_print_ipm("ipm_for_crash_node_1.txt");}
			// Convert from concentration
			REACT_GEM::ConcentrationToMass ( in,1 ); // I believe this is save for MPI
			// now we calculate kinetic constraints for GEMS!
			REACT_GEM::CalcLimits ( in );
			//Get data
			REACT_GEM::SetReactInfoBackGEM ( in ); // this should be also save for MPI
			// take values from old B volume for comparison
			oldvolume=m_Vs[in];

			// Order GEM to run
			//if(in == 100) { m_Node->na->GEM_write_dbr("dbr_for_crash_node_2.txt"); m_Node->na->GEM_print_ipm("ipm_for_crash_node_2.txt");}
			if ( flag_gem_smart )
			{
				dBR->NodeStatusCH = NEED_GEM_SIA;
				m_NodeStatusCH[in] = m_Node->GEM_run ( true );
			}
			else
			{
				dBR->NodeStatusCH = NEED_GEM_AIA; // first try without simplex using old solution
				m_NodeStatusCH[in] = m_Node->GEM_run ( false );
			}


			if ( ! ( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) )   // ups...failed..try again without SIMPLEX
			{
				dBR->NodeStatusCH = NEED_GEM_AIA;
				m_NodeStatusCH[in] = m_Node->GEM_run ( false );
			}

// test for bad GEMS and for volume changes bigger than 10% ...maximum 5 failed nodes per process.....
			if (
			    ! ( m_NodeStatusCH[in] == OK_GEM_AIA || m_NodeStatusCH[in] == OK_GEM_SIA ) ||
			    ( ( ( abs ( oldvolume-dBR->Vs ) /oldvolume ) >0.1 ) && ( flowflag != 3 ) ) // not for Richards flow
			)
			{
				cout << "Error: Main Loop failed when running GEM on Node #" << in << "." << endl << "Returned Error Code: " << m_NodeStatusCH[in] ;
				cout << " or GEM weird result at node " << in << " volume " <<  dBR->Vs << " old volume " <<oldvolume << endl;
				cout  << " continue with last good solution for this node" << endl;
				m_Node->GEM_write_dbr ( "dbr_for_crash_node_fail.txt" );
				m_Node->GEM_print_ipm ( "ipm_for_crash_node_fail.txt" );
				// exit ( 1 );
				node_fail=1;
				repeated_fail +=1;
				if ( repeated_fail >m_max_failed_nodes )
				{
					cout << "GEFMS: "<<repeated_fail << "nodes failed this timestep, check chemical system!" << endl;
#ifdef USE_MPI_GEMS
                        MPI_Finalize();           //make sure MPI exits
#endif
                        exit ( 1 );
                     }                            // we do not tolerate more than three failed nodes -> something wrong with chemistry/time step?

                     idummy=1;                    // end the while loop because we have a result
                     // cout << " GEM sucess with scaling of : " << scalfac << endl;

                  }                               // end loop if initial gems run fails
                  else
                  {

				// this is if gem run is ok
				REACT_GEM::GetReactInfoFromGEM ( in ); //from here on we have buffer values if GEMS_MPI is defined
				// check if result is ok..we do it via volume

//		if (abs((m_xDC[in*nDC+70] + m_xDC[in*nDC+69]+ m_xDC[in*nDC+122]+m_xDC[in*nDC+125]+ m_xDC[in*nDC+26]+ m_xDC[in*nDC+34]+ (2.0*m_xDC[in*nDC+35])+ (3.0*m_xDC[in*nDC+36]))- m_Node->pCNode()->bIC[4]) > 1.0e-10) {
//	if (abs((m_xDC_pts[in*nDC+70]*m_porosity[in] - m_xDC[in*nDC+70])/ m_xDC[in*nDC+70] )>1.0e-4 ) {
//		m_Node->na->GEM_write_dbr ( "dbr_for_crash_node127.txt" );
//		cout << "DEBUG: node "<< in << " Cl- " << m_xDC[in*nDC+70] <<" ClO4-  "<< m_xDC[in*nDC+69]<< " surf w1 "<< m_xDC[in*nDC+122]<< " surf w2"<< m_xDC[in*nDC+125]<<" FeCl+ "<< m_xDC[in*nDC+26]<< " FeCl+2 "<< m_xDC[in*nDC+34]<<" FeCl2+ "<< m_xDC[in*nDC+35]<<" FeCl3@ "<< m_xDC[in*nDC+36]<<endl;
//		cout << "DEBUG: residual " << (m_xDC[in*nDC+70] + m_xDC[in*nDC+69]+ m_xDC[in*nDC+122]+m_xDC[in*nDC+125]+ m_xDC[in*nDC+26]+ m_xDC[in*nDC+34]+ (2.0*m_xDC[in*nDC+35])+ (3.0*m_xDC[in*nDC+36]))- m_Node->pCNode()->bIC[4] << endl;
//                cout << " DEBUG node " << in << " rel error " << m_Node->pCNode()->rMB[4]/m_Node->pCNode()->bIC[4] << " previous amount " << m_xDC_pts[in*nDC+70]*m_porosity[in] << " now " << m_xDC[in*nDC+70]<< " rel change " << (m_xDC_pts[in*nDC+70]*m_porosity[in] - m_xDC[in*nDC+70])/ m_xDC[in*nDC+70] <<endl;
//         }

				/*   moved the volume check in order to be able to proceed at some nodes....
				                               if ( ( m_Vs[in]/oldvolume ) >2.0 ) {
				                                        cout << "GEM weird result at node " << in << " volume " << m_Vs[in] << " old volume " <<oldvolume << endl;
				                                        m_Node->na->GEM_write_dbr ( "dbr_for_crash_node_4.txt" );
				                                        m_Node->na->GEM_print_ipm ( "ipm_for_crash_node_4.txt" );
				#ifdef USE_MPI_GEMS
				                                        MPI_Finalize();  //make sure MPI exits
				#endif

				                                        exit ( 1 );
				                                }
				*/
			}
			if ( node_fail <1 )   // this needs to be done with buffer variables for MPI
			{
				// CALC kintetic constrains
				REACT_GEM::CalcReactionRate ( in, m_T[in] ); // check for kinetics is done in the subroutine for each species separately

				// test kinetics!
                     //Get data
                     //		REACT_GEM::SetReactInfoBackGEM ( in ); // this should be also save for MPI
                     //		m_NodeStatusCH[in] = m_Node->GEM_run ( scalfac, false );
                     // CALC kintetic constrains
                     //		cout << " second pass kinetics m_NodeStatusCH[in] " << m_NodeStatusCH[in]<< endl;
                     //		REACT_GEM::CalcReactionRate ( in, m_T[in], m_P[in] ); // check for kinetics is done in the subroutine for each species separately
                     // test kinetics finished
                     // calculate the chemical porosity
                     if ( flag_porosity_change ==1 )
                     {
                        REACT_GEM::CalcPorosity ( in );
                     }

                     // Convert to concentration  ..
				REACT_GEM::MassToConcentration ( in, 0 /*not a failed node */ );

			}
			else
			{
// seems not to work...better copy old values!!!    REACT_GEM::MassToConcentration ( in, 1 /* node failed */); //This is for failed nodes!!!!!
				RestoreOldSolution ( in );
				node_fail=0;

			}
#ifdef USE_MPI_GEMS
			REACT_GEM::CopyToMPIBuffer ( in ); //copy data to MPI buffer in any case
#endif
		} //end if check for boundary node
	} // end for loop for all nodes
#ifdef USE_MPI_GEMS
            // For MPI scheme, gather the data here.
            REACT_GEM::GetGEMResult_MPI();
            REACT_GEM::CleanMPIBuffer();
            // we sort it every timestep...this can be improved by a parallel sorting and/or more effective sorting algorithm
            if ( aktueller_zeitschritt%m_IterDoneIndexSort == 0 )
            {
               if ( myrank == 0 )  cout << "GEMS: Run_Main resorting. timestep " << aktueller_zeitschritt << endl;
               SortIterations ( m_IterDoneCumulative,m_IterDoneIndex,nNodes );
               for ( ii=0;ii<nNodes;ii++ )
               {
                  //			cout << "Index: " << m_IterDoneIndex[ii] << " Iterations: " << m_IterDoneCumulative[ii] << endl;
                  m_IterDoneCumulative[ii]=0;
               }
            }
#endif
            // this is done in main loop?        REACT_GEM::SetReactInfoBackMassTransport(1);
            cout << " GEM  run successful. "  << endl;
            return 0;
         }

         int REACT_GEM::GetHeatFlag_MT ( void )
         {
            //heat transport
            for ( size_t i=0; i < pcs_vector.size() ; i++ )
            {
               m_pcs = pcs_vector[i];
               //                if ( m_pcs->pcs_type_name.compare ( "HEAT_TRANSPORT" ) == 0 ) {
                                                  // TF
               if ( m_pcs->getProcessType() == HEAT_TRANSPORT )
               {
                  return 1;
               }
            }
            return 0;
         }

         int REACT_GEM::GetFlowType_MT ( void )
         {
            //flow type
            for ( size_t i=0; i < pcs_vector.size(); i++ )
            {
               m_pcs = pcs_vector[i];
               //                if ( m_pcs->pcs_type_name.compare ( "GROUNDWATER_FLOW" ) ==0 ) {
               if ( m_pcs->getProcessType() == GROUNDWATER_FLOW)
               {
                  m_flow_pcs = m_pcs;
                  return 1;
                  //                } else if ( m_pcs->pcs_type_name.compare ( "LIQUID_FLOW" ) ==0 ) {
               }
               else if ( m_pcs->getProcessType() == LIQUID_FLOW)
               {
                  m_flow_pcs = m_pcs;
                  return 2;
                  //                } else if ( m_pcs->pcs_type_name.compare ( "RICHARDS_FLOW" ) ==0 ) {
               }
               else if ( m_pcs->getProcessType() == RICHARDS_FLOW)
               {
                  m_flow_pcs = m_pcs;
                  return 3;
                  //                } else if ( m_pcs->pcs_type_name.compare ( "MULTI_PHASE_FLOW" ) ==0 ) {
               }
               else if ( m_pcs->getProcessType() == MULTI_PHASE_FLOW )
               {
               }
               m_flow_pcs = m_pcs;
               return 4;
            }
            return 0;
         }

         void REACT_GEM::GetFluidProperty_MT ( void )
         {
            m_FluidProp = MFPGet ( "LIQUID" );
         }

         long REACT_GEM::GetNodeNumber_MT ( void )
         {
            long number;
            //------------read number of nodes--------------
            for ( size_t i=0; i < pcs_vector.size(); i++ )
            {
               m_pcs = pcs_vector[i];
               //		if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) ==0 ) {
               if ( m_pcs->getProcessType() == MASS_TRANSPORT )
               {
                  number = ( long ) m_pcs->m_msh->GetNodesNumber ( false );
                  return number;
               }
            }
            //------------end of reading number of nodes----
            return 0;
         }

         long REACT_GEM::GetElemNumber_MT ( void )
         {
            long number;
            //------------read number of elems--------------
            for (size_t i=0; i < pcs_vector.size(); i++ )
            {
               m_pcs = pcs_vector[i];
               //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) ==0 ) {
               if ( m_pcs->getProcessType() == MASS_TRANSPORT)
               {
                  number = ( long ) m_pcs->m_msh->ele_vector.size();
                  return number;
               }
            }
            return 0;
         }

         double REACT_GEM::GetTempValue_MT ( long node_Index, int timelevel )
         {
            int indx;
            double temp;

            if ( heatflag == 1 )
            {
               m_pcs = PCSGet ( "HEAT_TRANSPORT" );

               indx = m_pcs->GetNodeValueIndex ( "TEMPERATURE1" ) +timelevel;
               temp = m_pcs->GetNodeValue ( node_Index, indx );

               //sysT[i] = m_pcs->GetNodeValue(i, indx1);
               //if (sysT0[i] <273.15) sysT0[i] += 273.15;  //ToDo �C->K
               //if (sysT[i] <273.15) sysT[i] += 273.15;  //ToDo �C->K
            }
            else
            {
               temp = m_gem_temperature;
            }
            return temp;
         }
         short REACT_GEM::SetTempValue_MT ( long node_Index, int timelevel, double temp )
         {
            int indx;

            if ( heatflag == 1 )
            {
               m_pcs = PCSGet ( "HEAT_TRANSPORT" );

               indx = m_pcs->GetNodeValueIndex ( "TEMPERATURE1" ) +timelevel;
               m_pcs->SetNodeValue ( node_Index, indx, temp );

               //sysT[i] = m_pcs->GetNodeValue(i, indx1);
               //if (sysT0[i] <273.15) sysT0[i] += 273.15;  //ToDo �C->K
               //if (sysT[i] <273.15) sysT[i] += 273.15;  //ToDo �C->K
               return 1;
            } else
            return 0;
         }

         double REACT_GEM::GetPressureValue_MT ( long node_Index, int timelevel )
         {
            //Get pressure value
            double pressure;
            int indx;
            pressure = 0.0;

            if ( flowflag > 0 )
            {
               GetFluidProperty_MT();
               switch ( flowflag )
               {
                  case 1:                         // for "GROUNDWATER_FLOW";
				if ( aktueller_zeitschritt > 0 )   // not for first time step
                     {

                        if ( gem_pressure_flag<1 )
                        {
                           // just set to 1.0 bar.
                           pressure = m_gem_pressure ;
                           break;
                        }
                        else
                        {
                           indx = m_flow_pcs->GetNodeValueIndex ( "HEAD" );
                                                  // The unit of HEAD is in meters
                           pressure = m_flow_pcs->GetNodeValue ( node_Index, indx + timelevel );
                           // cout << pressure << " " << " timelevel "<<timelevel;
                           // change the pressure unit from hydraulic head to Pa.
                           pressure = Pressure_M_2_Pa ( pressure , m_FluidProp->Density() );
                           //y cout << pressure << " density" << m_FluidProp->Density()<<endl;
                           // add atmospheric pressure
                        }

                        if ( pressure <= 0.0      /*valcumm suction in groundwater is not so realistic*/
                           || pressure > 1.0e+15  /*some very high pressure*/
                           )
                        {
                                                  // then set it to 1.0 bar = 1.0e5 Pa;
                           cout << " high pressure " << pressure << endl;
                           pressure = 1.0e+05;
                        }
                        break;
                     }
                     else
                     {
                        pressure=m_gem_pressure;  // set to m_gem_pressure default is 1e5P
                        break;
                     }

                  case 2:                         // for "LIQUID_FLOW", not tested!!!
                                                  // not for first time step
                     if ( aktueller_zeitschritt > 0 )
                     {
                        cout << "Liquid flow not supported - quitting program" << endl;
#ifdef USE_MPI_GEMS
                        MPI_Finalize();           //make sure MPI exits
#endif
                        exit ( 1 );

                        indx = m_flow_pcs->GetNodeValueIndex ( "PRESSURE1" ) +timelevel;
                                                  // The unit of HEAD is in meters
                        pressure = m_flow_pcs->GetNodeValue ( node_Index, indx );

                        // change the pressure unit from meters of water to bar.
                        // pressure = Pressure_M_2_Bar ( pressure , m_FluidProp->Density() );
                        // add atmospheric pressure
                        if ( pressure < 0.0       /*valcumm suction in groundwater is not so realistic*/
                           || pressure > 1.0e+15  /*some very high pressure*/
                           ) pressure = 1.0e+5;   // then set it to 1.0 bar;
                        break;
                     }
                     else
                     {
                        // just set to 1.0 bar.
                        pressure = 1.0e+05 ;
                        break;
                     }
                  case 3:                         // for "RICHARDS_FLOW", not tested!!!
                                                  // not for first time step
                     if ( aktueller_zeitschritt > 0 )
                     {
                        pressure = m_gem_pressure;
                        //	indx = m_flow_pcs->GetNodeValueIndex ( "PRESSURE1" ) +timelevel;
                        //	pressure = m_flow_pcs->GetNodeValue ( node_Index, indx ); // The unit of HEAD is in meters

                        // change the pressure unit from meters of water to bar.
                        //	pressure = Pressure_M_2_Bar ( pressure , m_FluidProp->Density() );
                        // add atmospheric pressure
                        //	pressure +=1.0;
                        if ( pressure < 0.0       /*valcumm suction in groundwater is not so realistic*/
                           || pressure > 1.0e+15  /*some very high pressure*/
                           ) pressure = 1.0e+5;   // then set it to 1.0 bar;
                        break;
                     }
                     else
                     {
                        // just set to 1.0 bar.
                        pressure = m_gem_pressure;
                        if ( pressure < 0.0       /*valcumm suction in groundwater is not so realistic*/
                           || pressure > 1.0e+15  /*some very high pressure*/
                           ) pressure = 1.0e+5;   // then set it to 1.0 bar;
                        break;
                     }
                  case 4:                         // MULTIPHASE ....not tested
                                                  // not for first time step
                     if ( aktueller_zeitschritt > 0 )
                     {
                        indx = m_flow_pcs->GetNodeValueIndex ( "PRESSURE1" );
                                                  // The unit of HEAD is in meters
                        pressure = m_flow_pcs->GetNodeValue ( node_Index, indx +timelevel );

                        // change the pressure unit from meters of water to bar.
                        // pressure = Pressure_M_2_Bar ( pressure , m_FluidProp->Density() );
                        // add atmospheric pressure

                        if ( pressure < 0.0       /*valcumm suction in groundwater is not so realistic*/
                           || pressure > 1.0e+15  /*some very high pressure*/
                           ) pressure = 1.0e+5;   // then set it to 1.0 bar;
                        break;
                     }
                     else
                     {
                        // just set to 1.0 bar.
                        pressure = 1.0e+05 ;
                        break;
                     }
                  default:
#ifdef USE_MPI_GEMS
                     if ( myrank == 0 /*should be set to root*/ )
#endif
                        cout<<"Error: Not implemented for the flow in GEM case!!!"<<endl;
                     pressure = 1.0e+05;
                     break;
               }                                  // end of switch case;
            }                                     // end of if (flow_flag);
            else
            {
               // if no valid flow pcs existing;
               cout<< "Warning: No valid flow process!!" <<endl;
            }
            return pressure;
         }

         short REACT_GEM::SetPressureValue_MT ( long node_Index, int timelevel, double pressure )
         {
            //Set pressure value
            int indx;
            indx = 0;
            if ( flowflag > 0 )
            {
               switch ( flowflag )
               {
                  case 1:
                     m_pcs = PCSGet ( "GROUNDWATER_FLOW" );
                     pressure = Pressure_Pa_2_M ( pressure ,  m_FluidProp->Density() );
                     indx = m_pcs->GetNodeValueIndex ( "HEAD" ) +timelevel;
                     m_pcs->SetNodeValue ( node_Index, indx, pressure );
                     break;
                  case 2:
                     m_pcs = PCSGet ( "LIQUID_FLOW" );
                     indx = m_pcs->GetNodeValueIndex ( "PRESSURE1" ) +timelevel;

                     m_pcs->SetNodeValue ( node_Index, indx, pressure );
                     break;
                  case 3:
                     // do nothing for Richards flow
                     break;
                     //				m_pcs = PCSGet ( "RICHARDS_FLOW" );
                     //				indx = m_pcs->GetNodeValueIndex ( "PRESSURE1" ) +timelevel;
                     //				pressure = Pressure_Bar_2_Pa ( pressure );
                     //				m_pcs->SetNodeValue ( node_Index, indx, pressure );
                     //				break;
                  case 4:
                     m_pcs = PCSGet ( "MULTI_PHASE_FLOW" );
                     indx = m_pcs->GetNodeValueIndex ( "PRESSURE1" ) +timelevel;

                     m_pcs->SetNodeValue ( node_Index, indx, pressure );
                  default:
#ifdef USE_MPI_GEMS
                     if ( myrank == 0 /*should be set to root*/ )
#endif
                        cout<< "Error: Not implemented for the flow in GEM case!!!" <<endl;
                     break;
               }
            }
            else
            {
#ifdef USE_MPI_GEMS
               if ( myrank == 0 /*should be set to root*/ )
#endif
                  cout<< "Warning: No valid flow process!!" <<endl;
               return 0;
            }
            return 1;
         }

         double REACT_GEM::GetComponentValue_MT ( long node_Index, string m_component, int timelevel )
         {
            double m_comp_value;
            m_comp_value = -1.0;
            for ( size_t i=0; i < pcs_vector.size(); i++ )
            {
               m_pcs = pcs_vector[i];
               //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) == 0 ) {
               if ( m_pcs->getProcessType() == MASS_TRANSPORT)
               {
                  if ( strcmp ( m_pcs->pcs_primary_function_name[0],m_component.c_str() ) == 0 )
                  {
                     m_comp_value = m_pcs->GetNodeValue ( node_Index,m_pcs->GetNodeValueIndex ( m_pcs->pcs_primary_function_name[0] ) +timelevel );
                  }
               }
            }
            if ( m_comp_value != -1.0 )
            {
               return m_comp_value;
            }
            else
            {
#ifdef USE_MPI_GEMS
               if ( myrank == 0 /*should be set to root*/ )
#endif
                  cout<< "Error: Corresponding Component NOT FOUND!!!" <<endl;
               return m_comp_value;
            }
         }

         short REACT_GEM::GetDCValue_MT ( long node_Index, int timelevel, double* m_DC, double* m_DC_pts ,double* m_DC_MT_delta )
         {
            string str;
            double /*DC_MT_pre,*/ DC_MT_cur;

            for ( int i=0; i < nDC ; i++ )
            {
               m_pcs = pcs_vector[i+1];           // dangerous!!
               //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) == 0 ) {
               if ( m_pcs->getProcessType() == MASS_TRANSPORT )
               {
                  //if ( m_pcs->m_msh->nod_vector[node_Index]->onBoundary() == false ) // do not update values for boundary node?

                  str = m_pcs->pcs_primary_function_name[0];
                  if ( str.compare ( "pH" ) != 0 && str.compare ( "pe" ) != 0 && str.compare ( "Eh" ) != 0 && str.compare ( "NodePorosity" ) != 0 )
                  {
                     // Get previous iteration mass transport concentration value
                     // DC_MT_pre = m_pcs->GetNodeValue ( node_Index,m_pcs->GetNodeValueIndex ( str ) +0 );
                     // Get current iteration mass transport concentration value
                     DC_MT_cur = m_pcs->GetNodeValue ( node_Index,m_pcs->GetNodeValueIndex ( str ) +timelevel );

                     * ( m_DC+i ) = DC_MT_cur;
                  }
               }
            }

            return 1;
         }

         short REACT_GEM::GetBValue_MT ( long node_Index, int timelevel, double* m_soluteB )
         {
            string str;
            double /*DC_MT_pre,*/ DC_B_cur;

            for ( int i=0; i < nIC ; i++ )
            {
               m_pcs = pcs_vector[i+1];           // dangerous!!
               //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) == 0 ) {
               if ( m_pcs->getProcessType() == MASS_TRANSPORT)
               {
                  //if ( m_pcs->m_msh->nod_vector[node_Index]->onBoundary() == false ) // do not update values for boundary node?

                  str = m_pcs->pcs_primary_function_name[0];

                  // Get previous iteration mass transport concentration value
                  // DC_MT_pre = m_pcs->GetNodeValue ( node_Index,m_pcs->GetNodeValueIndex ( str ) +0 );
                  // Get current iteration mass transport concentration value
                  DC_B_cur = m_pcs->GetNodeValue ( node_Index,m_pcs->GetNodeValueIndex ( str ) +timelevel );

                  * ( m_soluteB+i ) = DC_B_cur;

               }
            }

            return 1;
         }

                                                  // same functionality as GetComponentValue_MT
         double REACT_GEM::GetDCValueSpecies_MT ( long node_Index, int timelevel, int iDc )
         {
            string str;
            double /*DC_MT_pre,*/ DC_MT_cur=0.0;

            m_pcs = pcs_vector[iDc+1];            // dangerous!!
            //        if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) == 0 ) {
            if ( m_pcs->getProcessType() == MASS_TRANSPORT)
            {

               str = m_pcs->pcs_primary_function_name[0];
               if ( str.compare ( "pH" ) != 0 && str.compare ( "pe" ) != 0 && str.compare ( "Eh" ) != 0 && str.compare ( "NodePorosity" ) != 0 )
               {
                  // Get previous iteration mass transport concentration value
                  // DC_MT_pre = m_pcs->GetNodeValue ( node_Index,m_pcs->GetNodeValueIndex ( str ) +0 );
                  // Get current iteration mass transport concentration value
                  DC_MT_cur = m_pcs->GetNodeValue ( node_Index,m_pcs->GetNodeValueIndex ( str ) +timelevel );

               }
               else
               {
                  cout << "Error in GetDCValueSpecies_MT ... return zero value" << endl;
                  DC_MT_cur=0.0;
               }
            }
            else
            {
               cout << "Error in GetDCValueSpecies_MT ... return zero value" << endl;
               DC_MT_cur=0.0;
            }

            return DC_MT_cur;
         }

         short REACT_GEM::GetSoComponentValue_MT ( long node_Index, int timelevel, double* m_Phase )
         {
            string str;
            int x_Component = 0;
            for (size_t i=0; i < pcs_vector.size() ; i++ )
            {
               m_pcs = pcs_vector[i];
               //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) == 0 ) {
               if ( m_pcs->getProcessType() == MASS_TRANSPORT)
               {

                  x_Component = -1;

                                                  //get the name of compound from MT;
                  str = m_pcs->pcs_primary_function_name[0];
                                                  //get the index of certain compound, -1: no match
                  x_Component = m_Node->Ph_name_to_xDB ( str.c_str() );
                  if ( x_Component > -1 )
                  {
                     * ( m_Phase+x_Component ) = m_pcs->GetNodeValue ( node_Index,m_pcs->GetNodeValueIndex ( str ) +timelevel );
                  }
                  else
                  {
                     //DisplayErrorMsg("Error: Corresponding Component NOT FOUND in GEM part!!");
                     //return 0;
                  }
               }
            }
            //DisplayErrorMsg("Error: MASS TRANSPORT NOT FOUND!!");
            return 1;
         }
         short REACT_GEM::SetDCValue_MT ( long node_Index, int timelevel, double* m_DC )
         {
            string str;
            for ( int i=0; i < nDC ; i++ )
            {

               m_pcs = pcs_vector[i+1];

               //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) == 0 ) {
               if ( m_pcs->getProcessType() == MASS_TRANSPORT)
               {
                  str = m_pcs->pcs_primary_function_name[0];
                  if ( str.compare ( "pH" ) != 0 && str.compare ( "pe" ) != 0 && str.compare ( "Eh" ) != 0 && str.compare ( "NodePorosity" ) != 0 )
                  {
                     if ( flag_iterative_scheme > 0 )
                     {
                        if ( CPGetMobil ( m_pcs->GetProcessComponentNumber() ) > 0 )
                        {
                           // m_pcs->eqs->b[node_Index] += m_xDC_Chem_delta[node_Index*nDC+i] / dt ;
                           m_pcs->SetNodeValue ( node_Index , m_pcs->GetNodeValueIndex ( str ) +timelevel , * ( m_DC+i ) );
                        }
                        else
                        {
                           m_pcs->SetNodeValue ( node_Index , m_pcs->GetNodeValueIndex ( str ) +timelevel , * ( m_DC+i ) );
                        }
                     }
                     else
                     {
                        m_pcs->SetNodeValue ( node_Index , m_pcs->GetNodeValueIndex ( str ) +timelevel , * ( m_DC+i ) );
                     }
                  }
               }
            }

            return 1;
         }

         short REACT_GEM::SetBValue_MT ( long node_Index, int timelevel, double* m_soluteB )
         {
            string str;
            for ( int i=0; i < nIC ; i++ )
            {

               m_pcs = pcs_vector[i+1];

               //                if ( m_pcs->pcs_type_name.compare ( "MASS_TRANSPORT" ) == 0 ) {
               if ( m_pcs->getProcessType() == MASS_TRANSPORT)
               {
                  str = m_pcs->pcs_primary_function_name[0];
                  if ( flag_iterative_scheme > 0 )
                  {
                     if ( CPGetMobil ( m_pcs->GetProcessComponentNumber() ) > 0 )
                     {
                        // m_pcs->eqs->b[node_Index] += m_xDC_Chem_delta[node_Index*nDC+i] / dt ;
                        m_pcs->SetNodeValue ( node_Index , m_pcs->GetNodeValueIndex ( str ) +timelevel , * ( m_soluteB+i ) );
                     }
                     else
                     {
                        m_pcs->SetNodeValue ( node_Index , m_pcs->GetNodeValueIndex ( str ) +timelevel , * ( m_soluteB+i ) );
                     }
                  }
                  else
                  {
                     m_pcs->SetNodeValue ( node_Index , m_pcs->GetNodeValueIndex ( str ) +timelevel , * ( m_soluteB+i ) );
                  }
               }
            }

            return 1;
         }

         // i_timestep 0: old timestep 1: new timestep
         int REACT_GEM::SetPorosityValue_MT ( long ele_Index,  double m_porosity_Elem , int i_timestep )
         {

            int idx;
            double old_porosity ;

            idx = -1;
            old_porosity = 0.0;

            CRFProcess* m_pcs = NULL;

            if ( flowflag > 0 )
            {
               GetFluidProperty_MT();
               switch ( flowflag )
               {
                  case 1:
                     m_pcs     = PCSGet ( "GROUNDWATER_FLOW" );
                     idx       = m_pcs->GetElementValueIndex ( "POROSITY" );

                     // set new porosity;
                     m_pcs->SetElementValue ( ele_Index, idx+i_timestep , m_porosity_Elem );
                     break;
                  case 2:
                     m_pcs = PCSGet ( "LIQUID_FLOW" );
                     idx=m_pcs->GetElementValueIndex ( "POROSITY" );
                     // always write into the new step
                     m_pcs->SetElementValue ( ele_Index,idx+i_timestep,m_porosity_Elem );
                     break;
                  case 3:
                     m_pcs = PCSGet ( "RICHARDS_FLOW" );
                     idx=m_pcs->GetElementValueIndex ( "POROSITY" );
                     // always write into the new step
                     m_pcs->SetElementValue ( ele_Index,idx+i_timestep,m_porosity_Elem );
                     break;
                  case 4:                         // kg44: do we have to update POROSITY_IL and POROSITY_SW?
                     m_pcs = PCSGet ( "MULTI_PHASE_FLOW" );
                     idx=m_pcs->GetElementValueIndex ( "POROSITY1" );
                     // always write into the new step
                     m_pcs->SetElementValue ( ele_Index,idx+i_timestep,m_porosity_Elem );
                     break;
                  default:
#ifdef USE_MPI_GEMS
                     if ( myrank == 0 /*should be set to root*/ )
#endif
                        cout << "Error: Not implemented for the flow in GEM case!!!" <<endl;
                     break;
               }
            }
            return 1;
         }

         int REACT_GEM::SetSourceSink_MT ( long in, double time_step_size /*in sec*/ ) {

         Water_ST_GEMS m_st;

         GetFluidProperty_MT();
         switch ( flowflag )
         {
            case 1:                               // groundwater flow
               m_st.index_node = in ;
               m_st.water_st_value  = m_excess_water[in]  / time_step_size;
                                                  // normalize with node volume
               m_st.water_st_value *= m_Node_Volume[in];
               m_flow_pcs->Water_ST_vec.push_back ( m_st );
               return 1;
               break;
            case 2:                               // liquid flow
               m_st.index_node = in ;
               m_st.water_st_value  = m_excess_water[in]  / time_step_size;
                                                  // normalize with node volume
               m_st.water_st_value *= m_Node_Volume[in];
               m_flow_pcs->Water_ST_vec.push_back ( m_st );
               return 1;
               break;
            case 3:                               // Richards flow
               m_st.index_node = in ;
               m_st.water_st_value  = m_excess_water[in]  / time_step_size;
                                                  // normalize with node volume
               m_st.water_st_value *= m_Node_Volume[in];
               m_flow_pcs->Water_ST_vec.push_back ( m_st );
               return 1;
               break;
            case 4:                               // multiphase flow...works with case 1 ...pressure saturation scheme
               m_st.index_node = in ;
               m_st.water_st_value  = m_excess_water[in]  / time_step_size;
                                                  // normalize with node volume
               m_st.water_st_value *= m_Node_Volume[in];
               m_flow_pcs->Water_ST_vec.push_back ( m_st );
               return 1;
               break;
            default:
#ifdef USE_MPI_GEMS
               if ( myrank == 0 /*should be set to root*/ )
#endif
                  cout<< "Error: Not implemented for the flow in GEM case!!!"<<endl;
               break;
         }
         return 0;
      }

      int REACT_GEM::FindWater_xDC ( void )
      {
         // initialization
         int rt = -1;
         int i;

         // loop over all the xDC names, and find the one that is water
         for ( i = 0 ; i < nDC ; i++ )
         {
            if ( dCH->ccDC[i] == 'W' )
            {
               rt = i ;
               return rt;
            }
         }
         return rt;
      }

      int REACT_GEM::Findhydrogen_bIC ( void )
      {
         // initialization
         int rt = -1;
         int i;

         // loop over all the bIC names, and find the one that is hydrogen
         for ( i = 0 ; i < nIC ; i++ )
         {
            if ( dCH->ccIC[i] == 'h' )
            {
               rt = i ;
               return rt;
            }
         }
         return rt;
      }

      int REACT_GEM::Findoxygen_bIC ( void )
      {
         // initialization
         int rt = -1;
         int i;

         // loop over all the bIC names, and find the one that is oxygen
         for ( i = 0 ; i < nIC ; i++ )
         {
            if ( dCH->ccIC[i] == 'o' )
            {
               rt = i ;
               return rt;
            }
         }
         return rt;
      }



      double REACT_GEM::Pressure_Pa_2_Bar ( double Pre_in_Pa )
      {
         return Pre_in_Pa / 1.0e+5;
      }

      double REACT_GEM::Pressure_Bar_2_Pa ( double Pre_in_Bar )
      {
         return Pre_in_Bar * 1.0e+5;
      }

      double REACT_GEM::Pressure_M_2_Bar ( double Pre_in_M, double flu_density )
      {
         return Pre_in_M * 9.81 * flu_density/1.0e+5  ;
      }

      double REACT_GEM::Pressure_Bar_2_M ( double Pre_in_Bar, double flu_density )
      {
         return Pre_in_Bar * 1.0e5 / 9.81 / flu_density ;
      }

      double REACT_GEM::Pressure_M_2_Pa ( double Pre_in_M, double flu_density )
      {
         return Pre_in_M * 9.81 * flu_density  ;
      }

      double REACT_GEM::Pressure_Pa_2_M ( double Pre_in_Pa, double flu_density )
      {
         return Pre_in_Pa / 9.81 / flu_density ;
      }

      double REACT_GEM::GetNodeAdjacentVolume ( long Idx_Node )
      {
         double volume;
         long Idx_Ele;
         int number_of_nodes;
         volume = 0.0;
         number_of_nodes = 0;

         CNode* m_Node;
         CElem* m_Elem;

         // get the pointer to current node;
         m_Node =  m_pcs->m_msh->nod_vector[Idx_Node];

         // loop over all the elements that adjacent to this node;
         for ( int i=0 ; i < ( long ) m_Node->connected_elements.size() ; i++ )
         {
            // get the index of current element;
            Idx_Ele = m_Node->connected_elements[i];

            // get the pointer of this element;
            m_Elem = m_pcs->m_msh->ele_vector[Idx_Ele];

            // get the number of nodes in this element;
            // given argument "false" means giving node number instead of Gauss points;
            number_of_nodes = m_Elem->GetNodesNumber ( false );

            // taking part of volume from this element;
            volume += m_Elem->GetVolume() / number_of_nodes ;
         }

         return volume;
      }

      // i_timestep: 0: old timestep 1: new timestep
      void REACT_GEM::ConvPorosityNodeValue2Elem ( int i_timestep )
      {
         long i,idx_Node;
         int j, number_of_nodes;
         double pormin=2.0,pormax=0.0;
         CNode* m_Node;
         CElem* m_Elem;

         for ( i=0 ; i < nElems ; i++ )
         {
            m_Elem =  m_pcs->m_msh->ele_vector[i];

            // first set the parameters to zero;
            m_porosity_Elem[i] = 0.0;
            number_of_nodes = ( int ) m_Elem->GetNodesNumber ( false );
            // then get the values from nodes
            for ( j = 0 ; j < number_of_nodes ; j++ )
            {
                                                  // get the connected nodes;
               idx_Node = m_Elem->GetNodeIndex ( j );
               m_Node = m_pcs->m_msh->nod_vector[idx_Node];
                                                  // this is arithmetric mean
               m_porosity_Elem[i] += m_porosity[idx_Node] / number_of_nodes;
               // here we use harmonic mean, as porosity is used for permeability/diffusivity changes....flux in the element is strongly influenced by the minimum values
               // cout << " porosity " << idx_Node <<" "<<m_porosity[idx_Node] << endl;
               // m_porosity_Elem[i] += 1.0/m_porosity[idx_Node] ; // this is for harmonic mean
            }
            //	cout << " Porosity: Element "<< i << " " <<m_porosity_Elem[i] << endl;

            //m_porosity_Elem[i] = (double) number_of_nodes / m_porosity_Elem[i];

                                                  // upper limit of porosity
            if ( m_porosity_Elem[i] >= max_possible_porosity ) m_porosity_Elem[i]=max_possible_porosity;
                                                  //lower limit of porosity..
            if ( m_porosity_Elem[i] <= min_possible_porosity ) m_porosity_Elem[i]=min_possible_porosity;

            pormin=min ( pormin,m_porosity_Elem[i] );
            pormax=max ( pormax,m_porosity_Elem[i] );

            // push back porosities
            SetPorosityValue_MT ( i , m_porosity_Elem[i] , i_timestep );
         }
         cout <<"min, max porosity: "<< pormin << " " << pormax <<endl;
      }

      void REACT_GEM::CalcPorosity ( long in )
      {
         int k;                                   //must be int for Ph_Volume !

         m_porosity[in]=0.0;
         for ( k=0;k<dCH->nPHb;k++ )
         {
            if ( ( dCH->ccPH[k] == 's' ) || ( dCH->ccPH[k] == 'x' ) || ( dCH->ccPH[k] == 'd' ) ) m_porosity[in] += m_Node->Ph_Volume ( k ) ;
            // * 1.0e-6;    // transform cm3 -> m3 !!
            // second not here in the loop. move to the place where divided by m_Vs.
         }

         // normalized by m_Vs
         //	m_porosity[in] = 1.0 - m_porosity[in] / ( m_Vs[in] * 1.0e-6 /*convert to cm3 here*/) ;
                                                  //kg44 this is the correct way to do it
         m_porosity[in] = 1.0 - ( m_porosity[in] ) ;

         //	cout <<" porosity:" << m_porosity[in] << " node: "<< in <<endl;
         // checking whether out of bounary. only needed if calc porosity is done

         if ( m_porosity[in] >= max_possible_porosity )
         {
            m_porosity[in]=max_possible_porosity; // upper limit of porosity
         }
         if ( m_porosity[in] <= min_possible_porosity )
         {
            m_porosity[in]=min_possible_porosity; //lower limit of porosity..
         }

         //	cout <<" porosity:" << m_porosity[in] << " node: "<< in <<" Vs "<<m_Vs[in]<<endl;
         return;
      }

void REACT_GEM::MassToConcentration ( long in /*idx of node*/ ,int i_failed )   //attention second argument is not timestep. I is a flag that indicates if we deal with failed nodes!!!!...do not get data from GEMS for this nodes
{
	// converting the value from moles to the value in mol/m^3 water.
	long i,j,k;
	int idx;
	double gas_volume,fluid_volume;
	double skal_faktor;

	// get the fluid volume
	if ( i_failed )
	{
		fluid_volume=m_fluid_volume[in];
		gas_volume=m_gas_volume[in];
	}
	else
	{
		fluid_volume=0.0;
		gas_volume=0.0;
		for ( k=0;k<dCH->nPHb;k++ )
		{
			if ( dCH->ccPH[k] == 'a' ) fluid_volume += m_Node->Ph_Volume ( k );
			if ( dCH->ccPH[k] == 'g' ) gas_volume += m_Node->Ph_Volume ( k );

			// * 1.0e-6;    // transform cm3 -> m3 !!
		}
	}

	if ( ( fluid_volume <= 0.0 ) )
	{
		cout <<"fluid volume negative or zero" << fluid_volume  << " node " << in << " " << m_fluid_volume[in] << endl;
#ifdef USE_MPI_GEMS
         MPI_Finalize();                          //make sure MPI exits
#endif
         exit ( 1 );
      }
      if ( ( gas_volume <0.0 ) )
      {
         cout <<"gas volume negative" << gas_volume  << " node " << in <<endl;
#ifdef USE_MPI_GEMS
         MPI_Finalize();                          //make sure MPI exits
#endif

         exit ( 1 );
      }

      // store volumes for reuse during concentration_to_mass
      m_fluid_volume[in]=fluid_volume;
      m_gas_volume[in]=gas_volume ;

      // calculate excess water: volume of fluid phase bigger/smaller than porosity;
      // requires updated porosity (calc porosity before calculation of execc water
      CRFProcess* m_pcs = NULL;
      GetFluidProperty_MT();
      switch ( flowflag )
      {
         case 1:                                  // groundwater flow
            m_excess_water[in] = m_fluid_volume[in]- m_porosity[in];

                                                  //mol amount of water in first phase
			skal_faktor= ( m_porosity[in] ) /m_fluid_volume[in]; //mol amount of water in first phase

            m_fluid_volume[in] = m_porosity[in];
            break;
         case 2:                                  // liquid flow
            m_excess_water[in] = m_fluid_volume[in]- m_porosity[in];
                                                  //mol amount of water in first phase
            m_xDC[in*nDC + idx_water]*= m_porosity[in]/m_fluid_volume[in];
                                                  // this should give the same value as the above statement ....is used in concentration2mass and mass2concentration
            m_xPA[0] *= m_porosity[in]/m_fluid_volume[in];
            // change also stored phase volume!!!!
            m_fluid_volume[in] = m_porosity[in];
            break;
         case 3:                                  // Richards flow
            m_pcs = PCSGet ( "RICHARDS_FLOW" );
            idx=m_pcs->GetNodeValueIndex ( "SATURATION1" );
            m_excess_water[in] = m_fluid_volume[in]- m_porosity[in]*m_pcs->GetNodeValue ( in,idx+1 );
            m_excess_gas[in] = m_gas_volume[in]- m_porosity[in]* ( 1.0-m_pcs->GetNodeValue ( in,idx+1 ) );
			skal_faktor= ( m_porosity[in]*m_pcs->GetNodeValue ( in,idx+1 ) ) /m_fluid_volume[in]; //mol amount of water in first phase
			m_fluid_volume[in] = m_porosity[in]*m_pcs->GetNodeValue ( in,idx+1 );

            break;
         case 4:                                  // multiphase flow...works with case 1 ...pressure saturation scheme
            m_pcs = PCSGet ( "MULTI_PHASE_FLOW" );
            idx=m_pcs->GetNodeValueIndex ( "SATURATION0" );
            m_excess_water[in] = m_fluid_volume[in]- m_porosity[in]*m_pcs->GetNodeValue ( in,idx+1 );
            m_excess_gas[in] = m_gas_volume[in]- m_porosity[in]* ( 1.0-m_pcs->GetNodeValue ( in,idx+1 ) );
            if ( m_fluid_volume[in] > 0.0 )
            {
               m_xDC[in*nDC + idx_water] *= m_pcs->GetNodeValue ( in,idx+1 ) * m_porosity[in]/m_fluid_volume[in];
            }
            else
            {
               m_xDC[in*nDC + idx_water] = 0.0;
            }

            break;
         default:
#ifdef USE_MPI_GEMS
            if ( myrank == 0 /*should be set to root*/ )
#endif
               cout<< "Error: Not implemented for the flow in GEM case!!!" <<endl;
            break;
      }
#ifdef USE_MPI_GEMS
      //	if ( fabs ( m_excess_water_buff[in] ) >= 0.01 ) cout << "node "<< in <<" m_excess_water" << m_excess_water_buff[in] <<endl;
      //	if ( fabs ( m_excess_gas_buff[in] ) >= 0.01 ) cout << "node "<< in <<" m_excess_gas" << m_excess_water_buff[in] <<endl;
#else
      //	if ( fabs ( m_excess_water[in] ) >= 0.01 ) cout << "node "<< in <<" m_excess_water " << m_excess_water[in] <<endl;
      //	if ( fabs ( m_excess_gas[in] ) >= 0.01 ) cout << "node "<< in <<" m_excess_gas " << m_excess_water[in] <<endl;
#endif
      // if positive, then source, otherwise sink.
      // change amount of fluid .... at the moment only water!!!!!!

      // if ( fabs ( m_excess_water[in] ) >= 1.e-10 ) cout << "node "<< in <<"m_excess_water" << m_excess_water[in] <<endl;
      // do we need to scale this accoring to the node volume?..I think so: it is done while transering this to GEMS

      for ( j=0 ; j < nIC; j++ )
      {
         i = in*nIC + j;
         m_bIC[i] -= m_bPS[j];                    // B vector without solute
         m_bPS[j] *=skal_faktor;                  // newly scaled first phase

      }
      for ( j=0 ; j <= idx_water; j++ )
      {
         i = in*nDC + j;
         m_xDC[i]*=skal_faktor;                   //newly scaled xDC including water /excluding rest
      }

      // here we do not need to account for different flow processes .... as long as we work only with one wetting phase (water!!!) ...for liqid CO2 we have to change something ....
      //****************************************************************************************
      if ( flag_transport_b ==0 )
      {
         for ( j=0 ; j < nDC; j++ )
         {
            i = in*nDC + j;
            // a better way to do would be for all mobile species...or?
            if ( ( dCH->ccDC[j] == 'S' ) || ( dCH->ccDC[j] == 'E' ) || ( dCH->ccDC[j] == 'T' ) || ( dCH->ccDC[j] == 'L' ) )
            {
               // only for solutions species in the water. and H+ as well as electron
               m_xDC[i] /= m_fluid_volume[in] ;
            }
            // else
            //  { cout << "Error calculating node volume " <<endl; }

            // also check if xDC values is negative
            if ( m_xDC[i] < 0.0 )
               m_xDC[i]=0.0;
            // cout << "mass2conc: l " << l << " j "<< j << " i " << i << " por " << m_porosity[l] << " m_xDC " <<  m_xDC[i] << endl;
         }
      }                                           // here we adjust the b vector
      else
      {
         /*		for ( j=0 ; j < nDC; j++ )
              {
                 i = l*nDC + j;
                 // a better way to do would be for all mobile species...or?
                 if ( ( dCH->ccDC[j] == 'S' ) || ( dCH->ccDC[j] == 'E' ) || ( dCH->ccDC[j] == 'T' ) || ( dCH->ccDC[j] == 'L' ) )
                 {
                    // only for solutions species in the water. and H+ as well as electron
                    m_xDC[i] /= water_volume ;
                 }
              } */
         for ( j=0 ; j < nIC; j++ )
         {
            i = in*nIC + j;

            // correct for water and transform bPS into concentrations
                                                  //carrier for zero(first phase)  is normally water!
            if ( idx_hydrogen == j )   m_soluteB[i] = m_bPS[j] - ( 2.0*m_xDC[in*nDC + idx_water] ) ;
            else if ( idx_oxygen == j )  m_soluteB[i] = m_bPS[j]- m_xDC[in*nDC + idx_water]  ;
            else  m_soluteB[i] = m_bPS[j] ;

            m_soluteB[i] /=m_fluid_volume[in];    // now these are the concentrations
         }

      }
      return;
   }

   void REACT_GEM::ConcentrationToMass ( long l /*idx of node*/, int i_timestep ) {
   // converting the value from mol/m^3 water to moles.
   long i,j;
   double water_volume=0.0;
   int idx;
   string ErrorOut;
   // get the water volume

   //	cout <<"water volume " << water_volume << endl;
   // I think here it has to be different than in MassToConcentration module, as after hydraulics, the volume is limited to porosity
   // As we work with an unit volume of "one", it should be save to directly take porosity as volume....
   // of course we have to differenciate between the different flow models
   CRFProcess* m_pcs = NULL;
   GetFluidProperty_MT();
   switch ( flowflag )
   {
      case 1:                                     // groundwater flow
         water_volume=m_fluid_volume[l];          // kg44: reuse stored fluid volume from last timestep!
         break;
      case 2:                                     // liquid flow
         water_volume=m_fluid_volume[l];
         break;
      case 3:                                     // Richards flow
			m_pcs = PCSGet ( "RICHARDS_FLOW" );
			idx=m_pcs->GetNodeValueIndex ( "SATURATION1" );
			//	water_volume=m_porosity[l]*m_pcs->GetNodeValue ( l,idx+i_timestep );
			water_volume=m_fluid_volume[l]*m_pcs->GetNodeValue ( l,idx+i_timestep );
         break;
      case 4:                                     // multiphase flow...works with case 1 ...pressure saturation scheme
         m_pcs = PCSGet ( "MULTI_PHASE_FLOW" );
         idx=m_pcs->GetNodeValueIndex ( "SATURATION0" );
         water_volume=m_porosity[l]*m_pcs->GetNodeValue ( l,idx+i_timestep );
         break;
      default:
#ifdef USE_MPI_GEMS
         if ( myrank == 0 /*should be set to root*/ )
#endif
            cout<<"Error: Not implemented for the flow in GEM case!!!"<<endl;
         break;
   }

   if ( ( water_volume < 1.0e-6 ) || ( water_volume > 1.0 ) )
   {
      cout <<"conctomass water volume " << water_volume << " at node: "<< l << endl;
#ifdef USE_MPI_GEMS
      MPI_Finalize();                             //make sure MPI exits
#endif

      exit ( 1 );
   }

   if ( flag_transport_b == 0 )
   {
      for ( j=0 ; j < nDC; j++ )
      {
         i = l*nDC + j;

         if ( ( dCH->ccDC[j] == 'S' ) || ( dCH->ccDC[j] == 'E' ) || ( dCH->ccDC[j] == 'T' ) || ( dCH->ccDC[j] == 'L' ) )
         {
            // only for solutions species in the water, as well as H+ and electrons
            m_xDC[i] *= water_volume ;
         }

         // also check if xDC values is negative
         if ( m_xDC[i] < 0.0 )
         {
            //			m_xDC[i] = abs ( m_xDC[i] );  //kg44 test if taking absolute value is more robust than setting component to zero!!!
            //m_xDC[i] = 0.0; //this is very bad if this is a base species....
            //			cout << "GEM: conc negativ for species " << j << " at node " << l << "value: " << m_xDC[i] ;
                                                  // make sure it is greater/eq zero and take value from last timestep ---> second argument is zero!
            m_xDC[i] = fabs ( GetDCValueSpecies_MT ( l, 0, ( int ) j ) );

            //			cout << " corrected to: " << m_xDC[i] << endl;
         }
         // cout << "conc2mass: l " << l << " j "<< j << " i " << i << " por " << m_porosity[l] << " m_xDC " <<  m_xDC[i] << endl;
      }
   }                                              // transport of b species...we have to construct the new b vector!
   else
   {
      for ( j=0 ; j < nIC; j++ )
      {
         i = l*nIC + j;
         // also check if xDC values is negative
         if ( m_soluteB[i] < 0.0 )
         {
                                                  // make sure it is greater/eq zero and take value from last timestep ---> second argument is zero!
            m_soluteB[i] = fabs ( GetDCValueSpecies_MT ( l, 0, ( int ) j ) );
         }
         m_soluteB[i] *= water_volume ;
                                                  //carrier for zero(first phase)  is normally water!
         if ( idx_hydrogen == j )   m_soluteB[i] += ( 2.0*m_xDC[l*nDC + idx_water] ) ;
         else if ( idx_oxygen == j )  m_soluteB[i] +=  m_xDC[l*nDC + idx_water]  ;

         m_bIC[i]+=m_soluteB[i];                  //updated B vector for GEMS
			// here we check again if the vector if negative or smaller than a minimum amount...if yes we add some stuff....
			// adding 10-6 Mol/m^3 should be save, as this corresponds aprox to 10-9 Mol/kg ..which is well above the accuracy for most calculations
			if ( m_bIC[i] <= 1.0e-8 )  m_bIC[i]=1e-6;
         //                        cout <<  " i " << i << " " << m_bIC[i] << endl;
      }

   }
   return;
}


void REACT_GEM::CopyCurXDCPre ( void )
{
   long i;
   for ( i=0 ; i < nNodes*nDC ; i++ )
   {
      m_xDC_pts[i] = m_xDC[i];
   }
	return;
}


void REACT_GEM::UpdateXDCChemDelta ( void )
{
   long i;
   for ( i=0 ; i < nNodes*nDC ; i++ )
   {
      m_xDC_Chem_delta[i] = m_xDC[i] - m_xDC_pts[i];
   }
	return;
}

void REACT_GEM::CopyCurBPre ( void )
{
	long i;
	for ( i=0 ; i < nNodes*nIC ; i++ )
	{
		m_soluteB_pts[i] = m_soluteB[i];
		m_bIC_pts[i] = m_bIC[i];
	}
	return;
}

double REACT_GEM::CalcSoluteBDelta ( long in )
{
	long i;
	double dummy=0;
	for ( i=0 ; i < nIC-1 ; i++ )
	{
		dummy = max ( dummy,abs ( m_soluteB[in*nIC+i] - m_soluteB_pts[in*nIC+i] ) /m_soluteB[in*nIC+i] );
	}
	return dummy;
}

void REACT_GEM::RestoreOldSolution ( long in )
{
	long i;
	for ( i=0 ; i < nIC-1 ; i++ )
	{
		m_soluteB[in*nIC+i] = m_soluteB_pts[in*nIC+i];
		m_bIC[in*nIC+i] = m_bIC_pts[in*nIC+i];
	}
	for ( i=0 ; i < nDC ; i++ )
	{
		m_xDC[in*nDC+i] = m_xDC_pts[in*nDC+i];
	}
	return;
}

// retrun the new permeability based on original permeability and old/new porosity
double REACT_GEM::KozenyCarman ( double k0, double n0, double n )
{
   double rt = 0.0;

   if ( k0 < 1.0e-20 || k0 > 1.0 || n0 <=0 || n0 >= 1 || n <=0 || n >= 1 )
      return 0.0;
   else
   {
      rt = k0 ;

      rt *=pow ( n / n0 , 3 );
   }

   return rt;
}


// retrun the new permeability based on original permeability and old/new porosity

double REACT_GEM::KozenyCarman_normalized ( double k0, double n0, double n )
{
   double rt = 0.0;

   if ( k0 < 1.0e-20 || k0 > 1.0 || n0 <=0 || n0 >= 1 || n <=0 || n >= 1 )
      return 0.0;
   else
   {
      rt = k0 ;

      rt *=pow ( n / n0 , 3 );

      rt *=pow ( ( 1 - n0 ) / ( 1 - n ) , 2 );
   }

   return rt;
}


bool GEMRead ( string base_file_name, REACT_GEM *m_GEM_p )
{
   cout << "GEMRead" << endl;
   char line[MAX_ZEILE];
   string sub_line;
   string line_string;
   ios::pos_type position;
   //========================================================================
   // file handling
   string gem_file_name;
   gem_file_name = base_file_name + GEM_FILE_EXTENSION;
   ifstream gem_file ( gem_file_name.data(),ios::in );
   if ( !gem_file.good() )
   {
      cout << "! Error in GEMRead: No GEM data !" << endl;
      return false;
   }
   gem_file.seekg ( 0L,ios::beg );
   //========================================================================
   // keyword loop
   while ( !gem_file.eof() )
   {
      gem_file.getline ( line,MAX_ZEILE );
      line_string = line;
      if ( line_string.find ( "#STOP" ) !=string::npos )
         return true;
      //----------------------------------------------------------------------
                                                  // keyword found
      if ( line_string.find ( "#GEM_PROPERTIES" ) !=string::npos )
      {
         position = m_GEM_p->Read ( &gem_file );
         gem_file.seekg ( position,ios::beg );
      }                                           // keyword found
   }                                              // eof
   return true;
}


ios::pos_type REACT_GEM::Read ( std::ifstream *gem_file )
{
   Kinetic_GEMS d_kin;                            // dummy kinetic vector
   int j;
   // Initialization----------------
   string sub_line;
   string line_string;
   string delimiter ( " " );
   bool new_keyword = false;
   string hash ( "#" );
   ios::pos_type position;
   string sub_string;
   bool new_subkeyword = false;
   string dollar ( "$" );
   string delimiter_type ( ":" );
   std::stringstream in;
   // ------------------------------

   // Loop over all the key words----------------------
   while ( !new_keyword )
   {
      new_subkeyword = false;
      position = gem_file->tellg();
      line_string = GetLineFromFile1 ( gem_file );
      if ( line_string.size() < 1 ) break;
      if ( line_string.find ( hash ) !=string::npos )
      {
         new_keyword = true;
         break;
      }
      // Key word "$GEM_INIT_FILE" .......................
      if ( line_string.find ( "$GEM_INIT_FILE" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> init_input_file_path;
         in.clear();
         continue;
      }
      // Key word "$GEM_INIT_FILE" .......................
      if ( line_string.find ( "$TRANSPORT_B" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> flag_transport_b;
         in.clear();
         continue;
      }

      // ......................................................
      // Key word "$FLAG_POROSITY_CHANGE" .......................
      if ( line_string.find ( "$FLAG_POROSITY_CHANGE" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> flag_porosity_change;
         in.clear();
         continue;
      }
      // ......................................................
      // Key word "$MIN_POROSITY" .............................
      if ( line_string.find ( "$MIN_POROSITY" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> min_possible_porosity;
         in.clear();
         continue;
      }
      // ......................................................
      // Key word "$MAX_POROSITY" .............................
      if ( line_string.find ( "$MAX_POROSITY" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> max_possible_porosity;
         in.clear();
         continue;
      }
      // ......................................................
      // Key word "$FLAG_COUPLING_HYDROLOGY" ..................
      if ( line_string.find ( "$FLAG_COUPLING_HYDROLOGY" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> flag_coupling_hydrology;
         in.clear();
         continue;
      }
      // ......................................................
      // Key word "$PERMEABILITY_POROSITY_MODEL" ..............
      if ( line_string.find ( "$PERMEABILITY_POROSITY_MODEL" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> flag_permeability_porosity;
         in.clear();
         continue;
      }
      // ......................................................
      // Key word "$ITERATIVE_SCHEME" .........................
      if ( line_string.find ( "$ITERATIVE_SCHEME" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> flag_iterative_scheme;
         in.clear();
         continue;
      }
      // ......................................................
      // Key word "$GEM_SIMPLEX" .........................
      if ( line_string.find ( "$GEM_SMART" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> flag_gem_smart;
         in.clear();
         continue;
      }
      // ......................................................
      // ......................................................
      // Key word "$ITERATIVE_SCHEME" .........................
      if ( line_string.find ( "$TEMPERATURE_GEM" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> m_gem_temperature;
         in.clear();
         continue;
      }
      // Key word "$ITERATIVE_SCHEME" .........................
      if ( line_string.find ( "$PRESSURE_GEM" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> m_gem_pressure;
         in.clear();
         gem_pressure_flag=1;
         continue;
      }
      // ......................................................
      // Key word "$INDEX_SORT_MPI"   only for MPI: every M_IterDoneIndexSort the NodeIndex is sorted according to gems iterations .........................
      if ( line_string.find ( "$INDEX_SORT_MPI" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> m_IterDoneIndexSort;
         in.clear();
         continue;
      }
      // ......................................................
      // ......................................................
      // Key word "$SHUFFLE_GEMS" randomize the order of gems executions to exclude effects due to small systematic changes .........................
      if ( line_string.find ( "$SHUFFLE_GEMS" ) !=string::npos )
      {
         // subkeyword found
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> m_ShuffleGems;
         in.clear();
         continue;
      }
      // ......................................................
		// Key word "$MAX_FAILED_NODES" limits the number of failed nodes .........................
		if ( line_string.find ( "$MAX_FAILED_NODES" ) !=string::npos )
		{
			// subkeyword found
			in.str ( GetLineFromFile1 ( gem_file ) );
			in >> m_max_failed_nodes;
			in.clear();
			continue;
		}
		// Key word "$MAX_FAILED_NODES" limits the number of failed nodes .........................
		if ( line_string.find ( "$MY_SMART_GEMS" ) !=string::npos )
		{
			// subkeyword found
			in.str ( GetLineFromFile1 ( gem_file ) );
			in >> m_diff_gems;
			in.clear();
			continue;
		}                // ......................................................

      // kg44 26.11.2008 read in parameters for kinetics and GEM
                                                  // subkeyword found
      if ( line_string.find ( "$KINETIC_GEM" ) !=string::npos )
      {

         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> d_kin.phase_name >> d_kin.kinetic_model;
			if ( d_kin.kinetic_model >=1 && d_kin.kinetic_model<=4 )
         {
            cout << " found kinetics " << d_kin.kinetic_model << endl;
            in >> d_kin.n_activities;
            if ( d_kin.n_activities >10 )
            {
               cout << "To many dependent species for GEM kinetic model "<< d_kin.n_activities<< endl;
#ifdef USE_MPI_GEMS
               MPI_Finalize();                    //make sure MPI exits
#endif
               exit ( 1 );
            }
            //                cout <<" activities " << n_activities << endl;
            // first general kinetic parameters
            //	0,1,2  double E_acid,E_neutral,E_base; // activation energies
            in >> d_kin.kinetic_parameters[0] >> d_kin.kinetic_parameters[1] >> d_kin.kinetic_parameters[2];
            //			cout << kinetic_parameters[0] << kinetic_parameters[1] << kinetic_parameters[1]<<endl;

            //      3-5  double k_acid, k_neutral,k_base; // dissolution/precipitation rate constants
            in >> d_kin.kinetic_parameters[3] >> d_kin.kinetic_parameters[4] >> d_kin.kinetic_parameters[5];
            //			cout << kinetic_parameters[3] << kinetic_parameters[4] << kinetic_parameters[5]<<endl;

            //      6-11  double q1,p1,q2,q3,p2,p3; // exponents for omega
            in >> d_kin.kinetic_parameters[6] >> d_kin.kinetic_parameters[7] >> d_kin.kinetic_parameters[8] >>d_kin.kinetic_parameters[9] >> d_kin.kinetic_parameters[10] >> d_kin.kinetic_parameters[11];
            for ( j=0;j<d_kin.n_activities;j++ )
            {

               in >> d_kin.active_species[j];
               //				cout << active_species[j] ;
               //      12,13,14  double n_1, n_2,n_3; // exponents for acidic, neutral and base cases for species one

               in >> d_kin.kinetic_parameters[j+12] >> d_kin.kinetic_parameters[j+13] >> d_kin.kinetic_parameters[j+14];
            }

         }
         in.clear();
         // next line is surface area
         in.str ( GetLineFromFile1 ( gem_file ) );
         in >> d_kin.surface_model;
         cout << d_kin.surface_model << endl;
			if ( d_kin.surface_model >= 1 || d_kin.surface_model <= 3 )  // surface model 1, 2 and 3....only one parameter...
         {
            in >> d_kin.surface_area[0];          // surface: m*m / mol
            cout << "surface area " << d_kin.surface_area[0] << endl;
         }
         in.clear();
         // push back vector
         m_kin.push_back ( d_kin );
      }                                           // subkeyword found

   }
   // End of looping over all the key words-----------
   return position;
}


/** calculate Reaction rates at a node in for all solids (for which kinetics is defined)!
* kg44 30.07.2009
* in: node temp: temperature
*/
void REACT_GEM::CalcReactionRate ( long in, double temp )
{

   int idx=0,i,ii;
   long j,k;
   double rrn=0.0, rrb=0.0,rra=0.0, sa=0.0;
   double R=8.31451070;                           // molar gas konstant [J K-1 mol-1]
   double aa=1.0,ab=1.0,ac=1.0;                   // activity products ...species are input from material file

   double sactivity;                              // dummy variable for extracting activities
   const char *species;

   /**
   int kinetic_model;  // only 1 = GEMS implemented right now
        int n_activities;  // number of species for activities
        string active_species[10];  // name for species ...maximum 10 names
   double kinetic_parameters[32];
   0,1,2  double E_acid,E_neutral,E_base; // activation energies
      3-5  double k_acid, k_neutral,k_base; // dissolution/precipitation rate constants at standart konditions
      6-11  double p1,q1,p2,q2,p3,q3; // exponents for omega
      12,13,14  double n_1, n_2,n_3; // exponents for acidic neutral and base cases for species one
      append for each species another set of n_1, n_2  n_3 (up to 10 sets -> up to ten species)

   */

   // loop over all kinetic vectors and do for the defined phases and get rate for each phase ...
   // this algorithm assumes that we have correct ordering of phases and components that belong into the phase

	for ( ii=0;ii < ( int ) m_kin.size();ii++ )
   {
      k=m_kin[ii].phase_number;

      if ( m_kin[ii].kinetic_model > 0 )          // do it only if kinetic model is defined take model
         // kinetic_model==1 dissolution+precipitation kinetics
         // kinetic_model==2 only dissolution (no precipitation)
         // kinetic_mocel==3 only precipitation (no dissolution)
      {
         mol_phase[in*nPH+k]=0.0;
         omega_phase[in*nPH+k]=0.0;
         dmdt[in*nPH+k]=0.0;

         for ( j=m_kin[ii].dc_counter;j<m_kin[ii].dc_counter+dCH->nDCinPH[k];j++ )
         {
                                                  // do not include surface complexation species!
            if ( ! ( dCH -> ccDC[j]  == '0' ) && ! ( dCH->ccDC[j]  == 'X' ) && ! ( dCH->ccDC[j]  == 'Y' ) && ! ( dCH->ccDC[j]  == 'Z' ) )
            {
               //				omega_phase[k] += CalcSaturationIndex ( j, in,tempC,press ); // loop over all components of the phase
                                                  // we need this later for solid solutions....
               omega_components[in*nDC+j] = m_Node->DC_a ( j );
                                                  // loop over all components of the phase
               omega_phase[in*nPH+k] += m_Node->DC_a ( j );

               mol_phase[in*nPH+k] += m_xDC[in*nDC+j];
            }
         }

         //		 	cout << omega_phase[k] << " " <<  mol_phase[k] << endl; // debug

			sa=REACT_GEM::SurfaceAreaPh ( ii,in ); // value for surface area in m^2...(specific surface area multiplied with volume of the phase)
			if ( m_kin[ii].kinetic_model == 4 ) // in the next part we try to mimic Crunchflow (at least partially)..could be also done in CalcLimits if this makes the code easier to read
			{
				if ( m_porosity_initial[in] < min_possible_porosity ) m_porosity_initial[in]=m_porosity[in]; //make sure m_porosity_initial is not zero! ...does not work properly with RESTART!!!!!
				if ( m_porosity_initial[in] < min_possible_porosity ) m_porosity_initial[in]=min_possible_porosity;
				if ( omega_phase[in*nPH+k] < 1.0 ) // this is the dissolution case
				{
					// in a first try we let it like it is!
				}
				else if ( omega_phase[in*nPH+k] > 1.0 ) // this is the precipitation case
				{

					sa *= pow ( m_porosity[in]/m_porosity_initial[in],0.66666666667 );
				}

			}

         aa=1.0;
         ab=1.0;
         ac=1.0;                                  // reset values for each phase!
         for ( i=0;i<m_kin[ii].n_activities;i++ )
         {
            species=m_kin[ii].active_species[i].c_str();
                                                  // loop over all the names in the list
            idx= m_Node-> DC_name_to_xCH ( species );
            if ( idx < 0 )
            {
               cout << "CalcReactionRate: no DC-name "<< m_kin[ii]. active_species[i] <<" found" << endl;
#ifdef USE_MPI_GEMS
               MPI_Finalize();                    //make sure MPI exits
#endif

               exit ( 1 );
            }
            else
            {
               // 	cout << "activities " <<aa << " " << ab << " " << ac <<" " << temp << endl;
               sactivity=m_Node->DC_a ( idx );    // extract activities (not activity coefficients!)
               // cout << "Activity " << sactivity << pow(10.0,sactivity)<< endl;
               aa *= pow ( sactivity,m_kin[ii].kinetic_parameters[12+i] );
               ab *= pow ( sactivity,m_kin[ii].kinetic_parameters[13+i] );
               ac *= pow ( sactivity,m_kin[ii].kinetic_parameters[14+i] );
               //    *Y_m,     // Molalities of aqueous species and sorbates [0:Ls-1]
               //    *Y_la,    // log activity of DC in multi-component phases (mju-mji0) [0:Ls-1]
               //    *Y_w,     // Mass concentrations of DC in multi-component phases,%(ppm)[Ls]
               //    *Gamma,   // DC activity coefficients in molal or other phase-specific scale [0:L-1]
            }
         }
         // 	cout << "activities " <<aa << " " << ab << " " << ac <<" " << temp << endl;
         // terms for each case
         rra= exp ( -1.0* m_kin[ii].kinetic_parameters[0]/R * ( 1.0/temp + 1.0/298.15 ) ) * aa * pow ( ( 1.0-pow ( omega_phase[in*nPH+k],m_kin[ii].kinetic_parameters[6] ) ),m_kin[ii].kinetic_parameters[7] );

         rrn= exp ( -1.0*m_kin[ii].kinetic_parameters[1]/R * ( 1.0/temp + 1.0/298.15 ) ) * ab * pow ( ( 1.0-pow ( omega_phase[in*nPH+k],m_kin[ii].kinetic_parameters[8] ) ),m_kin[ii].kinetic_parameters[9] );

         rrb= exp ( -1.0*m_kin[ii].kinetic_parameters[2]/R * ( 1.0/temp + 1.0/298.15 ) ) * ac * pow ( ( 1.0-pow ( omega_phase[in*nPH+k],m_kin[ii].kinetic_parameters[10] ) ),m_kin[ii].kinetic_parameters[11] );

                                                  // rate is scaled to the total amount available via the surface area
         dmdt[in*nPH+k]= -1.0 * sa * ( pow ( 10.0,m_kin[ii].kinetic_parameters[3] ) * rra + pow ( 10.0,m_kin[ii].kinetic_parameters[4] ) * rrn + pow ( 10.0,m_kin[ii].kinetic_parameters[5] ) * rrb );
         // cout << "dmdt " <<dmdt[in*nPH+k] << " sa "<< sa <<" rra "<< rra << " rrn "<< rrn << " rrb "<< rrb << " m_gam " << m_gam[idx]<< " dmdt " << dmdt[in*nPH+k]<< endl;

         // test for NaN!! ---seems necessary as sometimes rra, rrn, rrb get Inf! ---seems enough to test the upper limit---this test does not resolve the real problem ;-)...probably pow(0.0,0.0) for rra,rrn,rrb ?
         if ( ! ( dmdt[in*nPH+k]<=1.0 ) && ! ( dmdt[in*nPH+k]>1.0 ) )
         {
            cout << "dmdt " <<dmdt << " is NaN " << " sa "<< sa <<" rra "<< rra << " rrn "<< rrn << " rrb "<< rrb << endl;
            dmdt[in*nPH+k]= 0.0;                  // no change!
         }

      }                                           //end if for kinetic model >1
   }

   return ;
}

/**
 * REACT_GEM::CalcLimitsInitial ( long in )
 * This is part of the OGS-GEMS kinetic implementation. It should be called during initialization phase, when
 * no information from the previous timestep is available. All kinetically controlled phases are set to their
 * initial values by assigning dll (lower limit) and dul (upper limit) to the xDC values from a restart file or from IC files (in
 * case transport is done with full speciation).
 * In case the simulation starts and no xDC restart values are available, one should create such a restart file e.g. by conducting a
 * a equilibrium simulation with one time-step and - if necessary - adjucst the values for the kinetically controlled phases e.g.
 * with an editor or by scripts.

 */
void REACT_GEM::CalcLimitsInitial ( long in )
{

   long ii,k,j;

   for ( j=0;j<nDC;j++ )
   {
      m_dll[in*nDC+j]=0.0;                        // set to zero
      m_dul[in*nDC+j]=1.0e+10;                    // very high number
   }

	if ( ! ( m_flow_pcs->GetRestartFlag() >=2 ) )
	{
		return;
	}

	for ( ii=0;ii< ( int ) m_kin.size();ii++ )
   {
      k=m_kin[ii].phase_number;

      if ( m_kin[ii].kinetic_model > 0 )          // do it only if kinetic model is defined take model
      {
         // kinetic_model==1 dissolution+precipitation kinetics
         // kinetic_model==2 only dissolution (no precipitation)
         // kinetic_mocel==3 only precipitation (no dissolution)
// we test if restart flag is set....if no....this will not work, as x_dc might be not correct
         for ( j=m_kin[ii].dc_counter;j<m_kin[ii].dc_counter+dCH->nDCinPH[k];j++ )
         {
            if ( ( dCH -> ccDC[j]  == '0' ) || ( dCH->ccDC[j]  == 'X' ) || ( dCH->ccDC[j]  == 'Y' ) || ( dCH->ccDC[j]  == 'Z' ) )
            {
               m_dll[in*nDC+j]=0.0;               // set to zero
               m_dul[in*nDC+j]=1.0e+10;           // very high number
            }
            else
            {
               m_dul[in*nDC+j]= m_xDC[in*nDC+j];
               m_dll[in*nDC+j]= m_xDC[in*nDC+j];
					if ( m_dll[in*nDC+j] < 0.0 ) m_dll[in*nDC+j]=0.0; // no negative masses allowed
               if ( m_dul[in*nDC+j] < 0.0 ) m_dul[in*nDC+j]=0.0;
               if ( m_dll[in*nDC+j] > m_dul[in*nDC+j] ) m_dll[in*nDC+j]=m_dul[in*nDC+j];
            }
         }
      }                                           //end kinetic model

   }                                              // end loop over phases

   return;
}

/** In this function we calculate the actual upper and lower metastability constraints for the GEMS solution
*  from the phase reaction rates (calculated at the previous time step)
*/
void REACT_GEM::CalcLimits ( long in )
{

   long ii,k,j;

   for ( j=0;j<nDC;j++ )
   {
      m_dll[in*nDC+j]=0.0;                        // set to zero
      m_dul[in*nDC+j]=1.0e+6;                     // very high number
   }

	for ( ii=0;ii< ( int ) m_kin.size();ii++ )
   {
      k=m_kin[ii].phase_number;

      if ( m_kin[ii].kinetic_model > 0 )          // do it only if kinetic model is defined take model
      {
         // kinetic_model==1 dissolution+precipitation kinetics
         // kinetic_model==2 only dissolution (no precipitation)
         // kinetic_mocel==3 only precipitation (no dissolution)
         for ( j=m_kin[ii].dc_counter;j<m_kin[ii].dc_counter+dCH->nDCinPH[k];j++ )
         {
// cout << "Kin debug " << in << " " << m_xDC[in*nDC+j] << " " << omega_phase[in*nPH+k] << " " << mol_phase[in*nPH+k]<< endl;
            // surface complexation species are not kinetically controlled -- 0 is old way...X is new way in DCH files
            if ( ( dCH -> ccDC[j]  == '0' ) || ( dCH->ccDC[j]  == 'X' ) || ( dCH->ccDC[j]  == 'Y' ) || ( dCH->ccDC[j]  == 'Z' ) )
            {
               m_dll[in*nDC+j]=0.0;               // set to zero
               m_dul[in*nDC+j]=1.0e+6;            // very high number
            }
            else
            {
               if ( omega_phase[in*nPH+k] > 1.0001 )
               {
                  // cout << mol_phase[in*nPH+k] << " " << dmdt[in*nPH+k]*dt << " " <<omega_components[in*nDC+j] <<" " << omega_phase[in*nPH+k]<< endl;
                  m_dul[in*nDC+j]= ( mol_phase[in*nPH+k] + dmdt[in*nPH+k]*dt ) * omega_components[in*nDC+j] /omega_phase[in*nPH+k] ;
                  // cout << m_dul[in*nDC+j] << endl;
                  if ( ! ( m_dul[in*nDC+j]<=1.0 ) && ! ( m_dul[in*nDC+j]>1.0 ) )
                  {

                                                  // no change!
                     m_dul[in*nDC+j]= m_xDC[in*nDC+j];
                  }
                                                  // give some freedom
                  m_dll[in*nDC+j]= m_xDC[in*nDC+j];
                  //					m_dll[in*nDC+j]= m_dul[in*nDC+j];

               }
					else if ( omega_phase[in*nPH+k] < 0.9999 )
               {
                  m_dll[in*nDC+j]= ( mol_phase[in*nPH+k] + dmdt[in*nPH+k]*dt ) *omega_components[in*nDC+j]  /omega_phase[in*nPH+k] ;
                  if ( ! ( m_dll[in*nDC+j]<=1.0 ) && ! ( m_dll[in*nDC+j]>1.0 ) )
                  {

                                                  // no change!
                     m_dll[in*nDC+j]= m_xDC[in*nDC+j];
                  }
                                                  // give some freedom
                  m_dul[in*nDC+j]= m_xDC[in*nDC+j];
                  //					m_dul[in*nDC+j]= m_dll[in*nDC+j];

               }
               else
               {
                  m_dul[in*nDC+j]= m_xDC[in*nDC+j];
                  m_dll[in*nDC+j]= m_xDC[in*nDC+j];
               }
               // do some corrections
               // kinetic_model==2 only dissolution (no precipitation)
               // kinetic_mocel==3 only precipitation (no dissolution)
					if ( ( m_kin[ii].kinetic_model==2 ) && ( m_dul[in*nDC+j] > m_xDC[in*nDC+j] ) ) {m_dul[in*nDC+j]= m_xDC[in*nDC+j];}
               if ( ( m_kin[ii].kinetic_model==3 ) && ( m_dll[in*nDC+j] < m_xDC[in*nDC+j] ) ) m_dll[in*nDC+j]= m_xDC[in*nDC+j];

               if ( ( m_xDC[in*nDC+j] < 1.0e-6 ) && ( omega_phase[in*nPH+k] >=1.0001 ) && ( m_dul[in*nDC+j]<1.0e-6 ) )
               {
                  m_dul[in*nDC+j]=1.0e-6;         // allow some kind of precipitation...based on saturation index for component value...here we set 10-6 mol per m^3 ..which is maybe 10-10 per litre ...?
                  m_dll[in*nDC+j]=m_dul[in*nDC+j];
               }
                                                  // no negative masses allowed
               if ( m_dll[in*nDC+j] < 0.0 ) m_dll[in*nDC+j]=0.0;
                                                  // no negative masses allowed..give some freedom
               if ( m_dul[in*nDC+j] <= 0.0 ) m_dul[in*nDC+j]=1.0e-6;
               if ( m_dll[in*nDC+j] > m_dul[in*nDC+j] ) m_dll[in*nDC+j]=m_dul[in*nDC+j];

            }

            //	if ((fabs((m_dul[in*nDC+j]- m_dll[in*nDC+j]))>0.0)) cout << "Kinetics for component no. " << j << " at node " << in << " m_xDC "  <<  m_xDC[in*nDC+j] << " m_dll, mdul " << m_dll[in*nDC+j] << " " << m_dul[in*nDC+j] << " diff " << m_dul[in*nDC+j]- m_dll[in*nDC+j] << endl; // give some debug output for kinetics
         }
      }                                           //end kinetic model

   }                                              // end loop over phases

   return;
}


// simplest case....scaling with a specific surface area per volume mineral
// phasenr: index for phase compnr: index of component which belongs to the phase and for which a specific surface area is defined
double REACT_GEM::SurfaceAreaPh ( long kin_phasenr,long in )
{
   double surf_area=0.0;

                                                  // now it is volume of the phase ..with respect to the unit volume in m^3
   surf_area=m_Node->Ph_Volume ( m_kin[kin_phasenr].phase_number) ;

   if ( m_kin[kin_phasenr].surface_model == 1 )
   {
                                                  // multiplication with specific surface area gives area
      surf_area *= m_kin[kin_phasenr].surface_area[0];
   }
   else if ( m_kin[kin_phasenr].surface_model == 2 )
   {
                                                  // constant surface area
      surf_area= m_kin[kin_phasenr].surface_area[0];
   }
	else if ( m_kin[kin_phasenr].surface_model == 3 )
	{
		surf_area *= m_kin[kin_phasenr].surface_area[0]/m_porosity[in];  // multiplication with specific surface area and division by porosity
	}
   else
   {
      surf_area=0.0;                              // no kinetics...safe solution
   }

   // cout << "phase " << kin_phasenr << " area default " << m_kin[kin_phasenr].surface_area[0] << " model " << m_kin[kin_phasenr].surface_model <<  " surface area: " << surf_area << endl;
   return surf_area;
}


#ifdef USE_MPI_GEMS
void REACT_GEM::GetGEMResult_MPI ( void )
{
   // Now gather the calculated values------------------------------------------------------------------------------
   MPI_Allreduce ( m_NodeHandle_buff, m_NodeHandle, nNodes, MPI_LONG, MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce ( m_NodeStatusCH_buff, m_NodeStatusCH, nNodes, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_IterDone_buff, m_IterDone, nNodes, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

   MPI_Allreduce ( m_Vs_buff, m_Vs, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_Ms_buff, m_Ms, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_Gs_buff, m_Gs, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_Hs_buff, m_Hs, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_IC_buff, m_IC, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_pH_buff, m_pH, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_pe_buff, m_pe, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_Eh_buff, m_Eh, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_porosity_buff, m_porosity, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_fluid_volume_buff, m_fluid_volume, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_gas_volume_buff, m_gas_volume, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_excess_water_buff, m_excess_water, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_excess_gas_buff, m_excess_gas, nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   MPI_Allreduce ( m_bIC_buff, m_bIC, nNodes*nIC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_bIC_dummy_buff, m_bIC_dummy, nNodes*nIC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_soluteB_buff, m_soluteB, nNodes*nIC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   MPI_Allreduce ( m_xDC_buff, m_xDC, nNodes*nDC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   MPI_Allreduce ( m_dll_buff, m_dll, nNodes*nDC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_dul_buff, m_dul, nNodes*nDC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_aPH_buff, m_aPH, nNodes*nPH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_xPH_buff, m_xPH, nNodes*nPH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_xPA_buff, m_xPA, nNodes*nPS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   MPI_Allreduce ( m_xDC_pts_buff, m_xDC_pts, nNodes*nDC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_xDC_MT_delta_buff, m_xDC_MT_delta, nNodes*nDC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( m_xDC_Chem_delta_buff, m_xDC_Chem_delta, nNodes*nDC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( omega_phase_buff, omega_phase, nNodes*nPH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( mol_phase_buff, mol_phase, nNodes*nPH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   MPI_Allreduce ( dmdt_buff, dmdt, nNodes*nPH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce ( omega_components_buff, omega_components, nNodes*nDC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   // --------------------------------------------------------------------------------------------------------------

}


void REACT_GEM::CleanMPIBuffer ( void )
{
	long in, ii;
   for ( in = 0; in < nNodes ; in++ )
   {

      m_IterDoneCumulative[in] += m_IterDone[in];
      m_NodeHandle_buff[in] = 0;
      m_NodeStatusCH_buff[in] = 0;
      m_IterDone_buff[in] = 0;

      m_Vs_buff[in] = 0.0;
      m_Ms_buff[in] = 0.0;
      m_Gs_buff[in] = 0.0;
      m_Hs_buff[in] = 0.0;
      m_IC_buff[in] = 0.0;
      m_pH_buff[in] = 0.0;
      m_pe_buff[in] = 0.0;
      m_Eh_buff[in] = 0.0;
      m_porosity_buff[in]= 0.0;
      m_fluid_volume_buff[in]= 0.0;
      m_gas_volume_buff[in]= 0.0;

      m_excess_water_buff[in] = 0.0;
      m_excess_gas_buff[in] = 0.0;

      for ( ii = 0; ii < nIC ; ii++ )
      {

         m_soluteB_buff[in*nIC + ii ] = 0.0;
         m_bIC_buff[in*nIC + ii ] = 0.0;
         m_bIC_dummy_buff [ in*nIC + ii ] = 0.0;
      }

      for ( ii = 0; ii < nDC ; ii++ )
      {

         m_xDC_buff[in*nDC+ii ] = 0.0;
         m_dul_buff[in*nDC+ii ] = 0.0;
         m_dll_buff[in*nDC+ii ] = 0.0;
         m_xDC_pts_buff[in*nDC+ii ] = 0.0;
         m_xDC_MT_delta_buff[in*nDC+ii ] = 0.0;
         m_xDC_Chem_delta_buff[in*nDC+ii ] = 0.0;
         omega_components_buff[in*nDC+ii ] = 0.0;

      }

      for ( ii = 0; ii < nPH ; ii++ )
      {

         m_xPH_buff[in*nPH+ii ] = 0.0;
         m_aPH_buff[in*nPH+ii ] = 0.0;
         omega_phase_buff[in*nPH+ii ] = 0.0;
         mol_phase_buff[in*nPH+ii ] = 0.0;
         dmdt_buff[in*nPH+ii ] = 0.0;
      }

      for ( ii = 0; ii < nPS ; ii++ )
      {
         m_xPA_buff[in*nPS+ii ] = 0.0;
      }

   }

   for ( in = 0; in < nElems ; in++ )
   {
      m_porosity_Elem_buff[in] = m_porosity_Elem[in];
   }
   return;
}


void REACT_GEM::CopyToMPIBuffer ( long in )
{
	long ii;
   m_NodeHandle_buff[in] = m_NodeHandle[in];
   m_NodeStatusCH_buff[in] = m_NodeStatusCH[in];
   m_IterDone_buff[in] = m_IterDone[in];

   m_Vs_buff[in] = m_Vs[in];
   m_Ms_buff[in] = m_Ms[in];
   m_Gs_buff[in] = m_Gs[in];
   m_Hs_buff[in] = m_Hs[in];
   m_IC_buff[in] = m_IC[in];
   m_pH_buff[in] = m_pH[in];
   m_pe_buff[in] = m_pe[in];
   m_Eh_buff[in] = m_Eh[in];
   m_porosity_buff[in]= m_porosity[in];
   m_fluid_volume_buff[in]= m_fluid_volume[in];
   m_gas_volume_buff[in]= m_gas_volume[in];

   m_excess_water_buff[in] = m_excess_water[in];
   m_excess_gas_buff[in] = m_excess_gas[in];

   for ( ii = 0; ii < nIC ; ii++ )
   {

      m_bIC_buff[in*nIC + ii ] = m_bIC[in*nIC + ii ];
      m_soluteB_buff[in*nIC + ii ] = m_soluteB[in*nIC + ii ];

      m_bIC_dummy_buff [ in*nIC + ii ] = m_bIC_dummy[ in*nIC + ii ];
   }

   for ( ii = 0; ii < nDC ; ii++ )
   {

      m_xDC_buff[in*nDC+ii ] = m_xDC[in*nDC+ii ];
      m_dul_buff[in*nDC+ii ] = m_dul[in*nDC+ii ];
      m_dll_buff[in*nDC+ii ] = m_dll[in*nDC+ii ];
      m_xDC_pts_buff[in*nDC+ii ] = m_xDC_pts[in*nDC+ii ];
      m_xDC_MT_delta_buff[in*nDC+ii ] = m_xDC_MT_delta[in*nDC+ii ];
      m_xDC_Chem_delta_buff[in*nDC+ii ] = m_xDC_Chem_delta[in*nDC+ii ];
      omega_components_buff[in*nDC+ii ] = omega_components[in*nDC+ii ];

   }

   for ( ii = 0; ii < nPH ; ii++ )
   {

      m_xPH_buff[in*nPH+ii ] = m_xPH[in*nPH+ii ];
      m_aPH_buff[in*nPH+ii ] = m_aPH[in*nPH+ii ];
      omega_phase_buff[in*nPH+ii ] = omega_phase[in*nPH+ii ];
      mol_phase_buff[in*nPH+ii ] = mol_phase[in*nPH+ii ];
      dmdt_buff[in*nPH+ii ] = dmdt[in*nPH+ii ];
   }

   for ( ii = 0; ii < nPS ; ii++ )
   {
      m_xPA_buff[in*nPS+ii ] = m_xPA[in*nPS+ii ];
   }

}


// this is a "slow" bubblesort, but it gives back the indexes of our sort
void REACT_GEM::SortIterations ( long *iterations, long *indexes, long len )
{
   long i,j,iTmp;
   // Set the index to match the current sort
   // order of the names
   for ( i = 0; i < len; i++ )
      indexes[i] = i;

   // Bubblesort. Sort the indexes
   for ( i = 0; i < len; i++ )
   {
      for ( j = 0; j < len; j++ )
      {
         if ( iterations[indexes[i]] < iterations[indexes[j]] )
         {
            // swap the indexes
            iTmp = indexes[i];
            indexes[i] = indexes[j];
            indexes[j] = iTmp;
         }
      }
   }
}
#endif                                            // end MPI

// function coming from GEMS coupling...extract node based porosities..does only
// work with GEMS coupling and if a species "POROSITY" or "NodePorosity" is defined
// georg.kosakowski@psi.ch 02.11.2009

double  REACT_GEM::GetNodePorosityValue ( long node_Index )
{

   double node_poros=-1.0;                        //default value is negative...that should not happen

   node_poros=m_porosity[node_Index];

   return node_poros;

}


// this is a "shuffle" algorithm that randomizes the order of the indexes
void REACT_GEM::ShuffleIterations ( long* indexes, long len )
{
   long i;
   vector<long> mydata ( indexes,indexes + len ); //constructor

   //Seeding random number generator
   std::srand ( std::time ( 0 ) );

   //assigning the numbers
   //        for(long i=0;i<len;++i)
   //        mydata[i]=indexes[i];

   //Shuffling the numbers
                                                  // using build in random generator
   std::random_shuffle ( mydata.begin(), mydata.end() );

   //Displaying the shuffled numbers
   for ( i=0;i<len;++i )
   {
      indexes[i]=mydata[i];
      //cout << i << " " <<indexes[i];
   }

}


//taken from rf_REACT_BRNS
int REACT_GEM::IsThisPointBCIfYesStoreValue ( long index, CRFProcess* m_pcs, double& value )
{
	for ( long p=0; p< ( int ) m_pcs->bc_node_value.size(); ++p )
   {
      if ( index == m_pcs->bc_node_value[p]->msh_node_number )
      {
         value = m_pcs->bc_node_value[p]->node_value;
         return 1;                                // Yes, found it.
      }
   }

   return 0;
}


int REACT_GEM::WriteReloadGem()
{
   long i,j;

   // node porosity ...necessary for restart
   // first test if m_poorosity exists (possible if no gem process exists)
   if ( !m_porosity ) return 2;
   string m_file_name = FileName +"_"+"m_porosity"+"_gem.asc";
   ofstream os ( m_file_name.c_str(), ios::trunc|ios::out );
   if ( !os.good() )
   {
      cout << "Failure to open file: "<<m_file_name << endl;
      abort();
   }
   os.precision ( 15 );                           // 15 digits accuracy seems enough? more fields are filled up with random numbers!
   os.setf ( ios_base::scientific,ios_base::floatfield );

   for ( i=0; i<nNodes; i++ )
   {
      os<<m_porosity[i] <<"  ";
      os<<endl;
   }
   os.close();
   // now write the fluid_volumes
   if ( !m_fluid_volume ) return 2;
   string m_file_name2 = FileName +"_"+"m_fluid_volume"+"_gem.asc";
   ofstream os2 ( m_file_name2.c_str(), ios::trunc|ios::out );
   if ( !os2.good() )
   {
      cout << "Failure to open file: "<<m_file_name2 << endl;
      abort();
   }
   os2.precision ( 15 );                          // 15 digits accuracy seems enough? more fields are filled up with random numbers!
   os2.setf ( ios_base::scientific,ios_base::floatfield );

   for ( i=0; i<nNodes; i++ )
   {
      os2<<m_fluid_volume[i] <<"  ";
      os2<<endl;
   }
   os2.close();

   // write b vector
   if ( !m_bIC ) return 2;
   string m_file_name3 = FileName +"_"+"m_bIC"+"_gem.asc";
   ofstream os3 ( m_file_name3.c_str(), ios::trunc|ios::out );
   if ( !os3.good() )
   {
      cout << "Failure to open file: "<<m_file_name3 << endl;
      abort();
   }
   os3.precision ( 15 );                          // 15 digits accuracy seems enough? more fields are filled up with random numbers!
   os3.setf ( ios_base::scientific,ios_base::floatfield );

   for ( i=0; i<nNodes; i++ )
   {
      for ( j=0;j<nIC;j++ )
      {
         os3<<m_bIC[i*nIC + j ]<< "  ";
         os3<<m_soluteB[i*nIC + j ]<< "  ";
         os3<<endl;
      }
   }
   os3.close();
   // write xDC vector
   if ( !m_xDC ) return 2;

   string m_file_name4 = FileName +"_"+"m_xDC"+"_gem.asc";
   ofstream os4 ( m_file_name4.c_str(), ios::trunc|ios::out );
   if ( !os4.good() )
   {
      cout << "Failure to open file: "<<m_file_name4 << endl;
      abort();
   }
   os4.precision ( 15 );                          // 15 digits accuracy seems enough? more fields are filled up with random numbers!
   os4.setf ( ios_base::scientific,ios_base::floatfield );

   for ( i=0; i<nNodes; i++ )
   {
      for ( j=0;j<nDC;j++ )
      {
         os4<<m_xDC[i*nDC + j ]<< "  ";
         os4<<endl;
      }
   }
   os4.close();

   return 1;
}


int REACT_GEM::ReadReloadGem()
{
   long i,j;
   string m_file_name = FileName +"_"+"m_porosity"+"_gem.asc";

   ifstream is ( m_file_name.c_str(), ios::in );
   if ( !is.good() )
   {
      cout << "Failure to open file: "<<m_file_name << endl;
      abort();
   }
   //
   // node porosity ...necessary for restart

   for ( i=0; i<nNodes; i++ )
   {

      is>>m_porosity[i];
      is>>ws;
   }
   is.close();

   string m_file_name2 = FileName +"_"+"m_fluid_volume"+"_gem.asc";

   ifstream is2 ( m_file_name2.c_str(), ios::in );
   if ( !is2.good() )
   {
      cout << "Failure to open file: "<<m_file_name2 << endl;
      abort();
   }

   for ( i=0; i<nNodes; i++ )
   {

      is2>>m_fluid_volume[i];
      is2>>ws;
   }
   is2.close();

   string m_file_name3 = FileName +"_"+"m_bIC"+"_gem.asc";

   ifstream is3 ( m_file_name3.c_str(), ios::in );
   if ( !is3.good() )
   {
      cout << "Failure to open file: "<<m_file_name3 << endl;
      abort();
   }

   for ( i=0; i<nNodes; i++ )
   {
      for ( j=0;j<nIC;j++ )
      {
         is3>>m_bIC[i*nIC + j ];
         is3 >> m_soluteB[i*nIC + j ];
         is3>>ws;
      }
   }
   is3.close();

   string m_file_name4 = FileName +"_"+"m_xDC"+"_gem.asc";

   ifstream is4 ( m_file_name4.c_str(), ios::in );
   if ( !is4.good() )
   {
      cout << "Failure to open file: "<<m_file_name4 << endl;
      abort();
   }

   for ( i=0; i<nNodes; i++ )
   {
      for ( j=0;j<nDC;j++ )
      {
         is4>>m_xDC[i*nDC + j ];
         is4>>ws;
      }
   }
   is4.close();

   return 1;
}


void REACT_GEM::WriteVTKGEMValues ( fstream &vtk_file )
{

  long i,j,k;
  double bdummy=0.0;
   // this is point data

    for ( j=0 ; j < nIC; j++ )
                {         
                vtk_file << "SCALARS " << dCH->ICNL[j] << " double 1" << endl;
                vtk_file << "LOOKUP_TABLE default" <<endl;                    
                //....................................................................
                for ( k=0;k<nNodes;k++ )                                              
                {                                                                     
                        bdummy=m_soluteB[k*nIC + j] * m_fluid_volume[k] ; //soluteB contains volume based concentrations
                        // now we have to add water                                                                     
                        if ( idx_hydrogen == j )        bdummy += ( 2.0*m_xDC[k*nDC + idx_water] ) ;   //carrier for zero(first phase)  is normally water!                                                                                                                                                
                        else if ( idx_oxygen == j )  bdummy +=  m_xDC[k*nDC + idx_water]  ;                                                          

                        bdummy+=m_bIC[k*nIC + j]; //add the solids
                        vtk_file <<" "<<  bdummy << endl; // and output
                }
    }

   // loop over speciation vector!
   for ( i=0;i<nDC;i++ )
   {
      vtk_file << "SCALARS " << dCH->DCNL[i] << " double 1" << endl;
      vtk_file << "LOOKUP_TABLE default" <<endl;
      //....................................................................
      for ( j=0;j<nNodes;j++ )
      {
         vtk_file <<" "<<  m_xDC[j*nDC + i ] << endl;
      }
   }
   // eh, pe, pH, Nodeporosity
   vtk_file << "SCALARS " << " pH " << " double 1" << endl;
   vtk_file << "LOOKUP_TABLE default" <<endl;
   //....................................................................
   for ( j=0;j<nNodes;j++ )
   {
      vtk_file <<" "<<  m_pH[j] << endl;
   }
   vtk_file << "SCALARS " << " pe " << " double 1" << endl;
   vtk_file << "LOOKUP_TABLE default" <<endl;
   //....................................................................
   for ( j=0;j<nNodes;j++ )
   {
      vtk_file <<" "<<  m_pe[j] << endl;
   }
   vtk_file << "SCALARS " << " Eh " << " double 1" << endl;
   vtk_file << "LOOKUP_TABLE default" <<endl;
   //....................................................................
   for ( j=0;j<nNodes;j++ )
   {
      vtk_file <<" "<<  m_Eh[j] << endl;
   }
   vtk_file << "SCALARS " << " NodePorosity " << " double 1" << endl;
   vtk_file << "LOOKUP_TABLE default" <<endl;
   //....................................................................
   for ( j=0;j<nNodes;j++ )
   {
      vtk_file <<" "<<  m_porosity[j] << endl;
   }
   vtk_file << "SCALARS " << " Fluidvolume " << " double 1" << endl;
   vtk_file << "LOOKUP_TABLE default" <<endl;
   //....................................................................
   for ( j=0;j<nNodes;j++ )
   {
      vtk_file <<" "<<  m_fluid_volume[j] << endl;
   }
   vtk_file << "SCALARS " << " ExcessVolume " << " double 1" << endl;
   vtk_file << "LOOKUP_TABLE default" <<endl;
   //....................................................................
   for ( j=0;j<nNodes;j++ )
   {
      vtk_file <<" "<<  m_excess_water[j] << endl;
   }
	vtk_file << "SCALARS " << " NodeVolume " << " double 1" << endl;
	vtk_file << "LOOKUP_TABLE default" <<endl;
	//....................................................................
	for ( j=0;j<nNodes;j++ )
	{
		vtk_file <<" "<<  m_Vs[j] << endl;
	}

	return;
}
#endif
