//-------------------------------------
// rf_REACT_GEM.h
// Haibing Shao 23.05.07
// haibing.shao@ufz.de
// GEM Reaction Package
// based on the PSI node-GEM source code
// using the node-GEM code from Paul Sherrer Institute (PSI)
//-------------------------------------
#ifndef RF_REACT_GEM_H
#define RF_REACT_GEM_H

#include "Configure.h"

#ifdef GEM_REACT

#include <time.h>
#include <cmath>
#include <string>
#include "rf_pcs.h"
// #include "rfmat_cp.h"
#include "node.h"
#include "rf_mfp_new.h"

/**
 * class REACT_GEM for coupling OGS with GEMS
 */
class REACT_GEM
{
   private:
      string ipm_input_file_path;
      string dbr_input_file_path;
      string dbr_bc_input_file_path;
      string dch_input_file_path;
      string init_input_file_path;
      string init_input_file_sig;

   public:
      REACT_GEM(void);
      ~REACT_GEM(void);

		/// Instance of TNode class
		/// HS: 06.2007 Only set one TNode here, repeatedly use its resources for GEM calculation
		TNode* m_Node; 

      // DATABR structure for exchange with GEMIPM
      DATACH* dCH;                                //pointer to DATACH
      DATABR* dBR;                                //pointer to DATABR

      // Read function for gem input file
      ios::pos_type Read(ifstream *gem_file);

      // Number of ICs, DCs, Phases and Phases-solutions kept in the node or elem;
      long nIC, nDC, nPH, nPS;

      // Data structure for each node to carry the chemical information (real FMT problems consider many nodes)
      // Names are consistent with the DataBridge structure (also see "\GEM\databr.h")
      long *m_NodeHandle, *m_NodeStatusCH, *m_IterDone, *m_IterDoneCumulative, *m_IterDoneIndex;
      long m_IterDoneIndexSort, m_ShuffleGems;
      // this is for porosity calculated on volume of solids
      double *m_porosity;
		 /// this is used for kinetic law number 4
		double  *m_porosity_initial;

		/// this we need for porosity coupling to groundwater flow & multiphase flow
		double *m_fluid_volume, *m_gas_volume;

		/// indexes, which one in the xDC vector is water, oxygen or hydrogen
		int idx_water, idx_hydrogen, idx_oxygen;

      double *m_T, *m_P, *m_Vs, *m_Ms,
         *m_Gs, *m_Hs, *m_IC, *m_pH, *m_pe,
         *m_Eh;

      double *m_xDC, *m_gam, *m_xPH, *m_aPH, *m_vPS, *m_mPS, *m_bPS,
         *m_xPA, *m_dul, *m_dll, *m_bIC, *m_bIC_dummy, *m_rMB, *m_uIC;

      double  *m_porosity_Elem, *m_porosity_Elem_buff;

		/// data for transport of IC
		double *m_soluteB, *m_soluteB_buff, *m_soluteB_pts, *m_bIC_pts;

      // previous time step DC values
      double *m_xDC_pts;                          // previous time step Concentration;
      double *m_xDC_MT_delta;                     // delta C from Mass Transport;
      double *m_xDC_Chem_delta;                   // delta C from Chemistry;

      double *m_excess_water;                     // excess water in m3/s for each node;
      double *m_excess_gas;                       // excess gas in m3/s for each node;
      double *m_saturation;
      double *m_Node_Volume;                      // Volume around the node;

                                                  // this we need for kinetics
      double *mol_phase, *omega_phase,*omega_components;

      double *dmdt;                               // kinetically controlled rates
      void CalcLimits ( long in);
      void CalcLimitsInitial ( long in);
      int *m_boundary;                            //holds marker for boundary nodes

      CFluidProperties *m_FluidProp;
      CRFProcess *m_pcs;                          // pointer to the PCS Class.
      CRFProcess *m_flow_pcs;                     // pointer to the flow PCS.


		/** Initialization of the GEM TNode Class
		*  return: 0-ok;
		*          1-loading init file failure;
		*          3-dch problem;
		*           4-dbr problem;
                */
      short Init_Nodes(string Project_path);      // Initialization of the GEM TNode Class
      //  return: 0-ok;
      //          1-loading init file failure;
      //          3-dch problem;
      //          4-dbr problem;

      short Init_RUN();                           // Run the node-GEM
      //  return: 0-ok;5-GEM does not converge
      short Run_MainLoop();
      //  return: 0-ok;5-GEM does not converge

      string Get_Init_File_Path(void);
      string Get_DCH_File_Path(void);
      string Get_IPM_File_Path(void);
      string Get_DBR_File_Path(void);

      int Set_Init_File_Path(string m_path);
      int Set_DCH_FILE_PATH(string m_path);
      int Set_IPM_FILE_PATH(string m_path);
      int Set_DBR_FILE_PATH(string m_path);

      bool Load_Init_File(string m_Project_path);
      long* mp_nodeTypes;

      //---flags------
      int initialized_flag;                       //0 - not initialized 1 - initialized
      int flag_node_element_based;                // 0 - node based; 1 - element based;
      int flag_iterative_scheme;                  // 0 - sequential non-iterative scheme;
      // 1 - standard iterative scheme;
      // 2 - symetric iterative scheme;
      // 3 - strang splitting scheme;
      int heatflag;                               //0-initialized and not heat transport;1-heat_transport;
      int flowflag;                               //0-initialized;1-GROUNDWATER_FLOW;2-LIQUID_FLOW;3-RICHARDS_FLOW;4-FLOW;
      int flag_porosity_change;                   //0-porosity change not coupled into transport; 1=coupled;
      int flag_coupling_hydrology;                //0-without coupling; 1=with coupling;
      int flag_permeability_porosity;             //0-no coupling; 1-Kozeny-Carman; 2-Kozeny-Carman normalized;
      int flag_gem_smart;                         // shall we work with faster simplex for GEM?
      int gem_pressure_flag;                      //shall we give a constant user defined pressure to gems?
      int flag_transport_b;                       //1: transport only dissolved components of b vector; 0: transport full speciation
		long m_max_failed_nodes; ///maximum number of failed nodes
      //--------------

      long nNodes;                                // number of all nodes;
      long nElems;                                // number of all elements;
      int GetHeatFlag_MT(void);
      int GetFlowType_MT(void);
      long GetNodeNumber_MT(void);
      long GetElemNumber_MT(void);
      void GetFluidProperty_MT(void);

      short GetInitialReactInfoFromMassTransport(int timelevel);
      short GetReactInfoFromMassTransport(int timelevel);
      short SetReactInfoBackMassTransport(int timelevel);
      void GetReactInfoFromGEM(long in);
      void SetReactInfoBackGEM(long in);
      // necessary for reload with gems
      int WriteReloadGem();
      int ReadReloadGem();

      double GetTempValue_MT(long node_Index, int timelevel);
      double GetPressureValue_MT(long node_Index, int timelevel);
      double GetComponentValue_MT(long node_Index, string m_component, int timelevel);
      short GetDCValue_MT(long node_Index, int timelevel, double* m_DC, double* m_DC_pts, double* m_DC_MT_delta);
      short GetBValue_MT ( long node_i, int timelevel, double *m_soluteB);
      short GetSoComponentValue_MT(long node_Index, int timelevel, double* m_Phase);
      double GetDCValueSpecies_MT ( long node_Index, int timelevel, int iDc );
      short SetTempValue_MT(long node_Index, int timelevel, double temp);
      short SetPressureValue_MT(long node_Index, int timelevel, double pressure);
      short SetDCValue_MT(long node_Index, int timelevel, double* m_DC);
      short SetBValue_MT(long node_Index, int timelevel, double* m_soluteB);

		int IsThisPointBCIfYesStoreValue ( long index, CRFProcess* m_pcs, double& value );/// taken from rf_REACT_BRNS

		/// Copy current values into previous time step values
		void CopyCurXDCPre ( void );
		void UpdateXDCChemDelta ( void );
		void CopyCurBPre ( void );
		double CalcSoluteBDelta ( long in );
		double m_diff_gems;
		void RestoreOldSolution ( long in );
		/// this is only for porosity interpolation to elemens
      void ConvPorosityNodeValue2Elem(int i_timestep);
      void CalcPorosity(long in);

      double min_possible_porosity, max_possible_porosity;
      void ScaleVolume_Water(long in);

      // Set porosity in Mass Transport
      int SetPorosityValue_MT(long ele_Index, double m_porosity_Elem, int i_timestep);
      int SetSourceSink_MT(long in, double time_step_size /*in sec*/);

      // find which one in xDC vector is water
      int FindWater_xDC(void);
      int Findhydrogen_bIC ( void );
      int Findoxygen_bIC ( void );
      //kg44 11/2008 for kinetics
		void CalcReactionRate ( long node, double temp );
		double SurfaceAreaPh ( long kin_phasenr,long in );

      // concentration related
      void ConcentrationToMass (long l /*idx of node*/,int i_timestep);
      void MassToConcentration (long l /*idx of node*/,int i_timestep);

      // Unit conversion for pressures
      double Pressure_Pa_2_Bar(double Pre_in_Pa);
      double Pressure_Bar_2_Pa(double Pre_in_Bar);
      double Pressure_M_2_Bar(double Pre_in_M, double flu_density );
      double Pressure_Bar_2_M(double Pre_in_Bar, double flu_density );
      double Pressure_M_2_Pa ( double Pre_in_M, double flu_density );
      double Pressure_Pa_2_M ( double Pre_in_Pa, double flu_density );
      // Calculate the volume of the nodes;
      // given argument is the index of one particular node;
      double GetNodeAdjacentVolume(long Idx_Node);

      // Permeability-Porosity relationship--------------------------------
      // they return the new permeability value
      double KozenyCarman( double k0 /*original permeability*/,
         double n0 /*original porosity*/,
         double n  /*new porosity*/);
      double KozenyCarman_normalized( double k0 /*original permeability*/,
         double n0 /*original porosity*/,
         double n  /*new porosity*/);
      // ------------------------------------------------------------------

      // GEMS mass scaling parameter
      double gem_mass_scale;
      // GEM temperature (without coupling to temperature)
      double m_gem_temperature;
      // GEM pressure (needed for Richards flow)
      double m_gem_pressure;

		/// Definition of buffer variables for MPI
      long *m_NodeHandle_buff, *m_NodeStatusCH_buff, *m_IterDone_buff;
                                                  // porosity buffer
      double *m_porosity_buff, *m_fluid_volume_buff, *m_gas_volume_buff;
      double *m_Vs_buff, *m_Ms_buff, *m_Gs_buff, *m_Hs_buff, *m_IC_buff, *m_pH_buff, *m_pe_buff, *m_Eh_buff;
      double *m_xDC_buff, *m_xPH_buff,*m_aPH_buff,*m_xPA_buff,*m_excess_water_buff,*m_excess_gas_buff,*m_dul_buff, *m_dll_buff, *m_Node_Volume_buff, *m_saturation_buff,*m_bIC_buff,*m_bIC_dummy_buff, *m_xDC_pts_buff, *m_xDC_MT_delta_buff, *m_xDC_Chem_delta_buff;

                                                  // this we need for kinetics
      double *omega_phase_buff, *mol_phase_buff, *dmdt_buff, *omega_components_buff;

#ifdef USE_MPI_GEMS
      // MPI implementation

      void CleanMPIBuffer(void);
      void CopyToMPIBuffer(long in);

      void GetGEMResult_MPI(void);

      void SortIterations(long *iterations, long *indexes, long len);
#endif

      int Random(long n);
      void ShuffleIterations ( long *indexes, long len );

      double GetNodePorosityValue( long node_Index);

      void WriteVTKGEMValues(fstream &vtk_file);

      typedef struct
      {
         //kg44 25.11.2008 kinetics...for coupling with GEMS
         //
         string phase_name;
         int phase_number;
         int dc_counter;
         int kinetic_model;                       // only 1 = GEMS implemented right now
         int n_activities;                        // number of species for activities
         string active_species[10];               // name for species ...maximum 10 names

/**	this vector holds the kinetic material parameters
*      0,1,2  double E_acid,E_neutral,E_base; // activation energies
*      3-5  double k_acid, k_neutral,k_base; // dissolution/precipitation rate constants
*      6-11  double p1,q1,p2,q2,p2,q2; // exponents for omega
*      12,13, 14  double n_1, n_2, n_3; // exponents for acidic neutral and base cases for species one
*      append for each species another set of n_1, n_2, n_3 (up to 10 sets -> up to ten species)
*/
			double kinetic_parameters[41];
         int surface_model;                       // currently only 1 implemented
         double surface_area[10];
      } Kinetic_GEMS;

      vector<Kinetic_GEMS> m_kin;

};

#define GEM_FILE_EXTENSION ".gem"
/** This is the function for reading the OGS-GEM specific parameters.
*  Todo: if .gems not exist, decouple gems
*/
extern bool GEMRead(string base_file_name, REACT_GEM *m_GEM_p);
#endif
#endif
