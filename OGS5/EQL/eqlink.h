//*******************************************************
//Class CEqlink: 
//First implementation:  MX 04/2006
//******************************************************

// #include "stdafx.h" /* MFC */

#ifdef CHEMAPP

#include <vector>
using namespace std;

//General gas constant; the following value was verified to be the same
//than is used within chemapp.
  const static double R_GAS = 8.31441;
//  const static double lambda=0.001;
  const static int verbosity = 0;


  class CEqlink 
  {private:
   friend class REACT;
//*******************************************************
//Part 1: members and functions for prepare()
//Based on parameterfile-s.h in EQLINK ( Helge C. Moog)
//******************************************************
// Short Description:
/*   ------------------
    Variables        :
   ------------------
       RR_GAS       : General gas constant
     *     PCSTOI_GAS   : Matrix(i,j) for composition of consituents of
     *     PCSTOI_AQ      gas- and aqueous phase
     *     PCSTOI_SOLIDS  Matrix(i,j,k) for all constituents of all
                          solid phases with respect to system components
     *     PCMMASS_GAS  : Vector(i) for mol. weights of constituents of gas-
     *     PCMMASS_AQ     ... aqueous phases
     *     PCMMASS_SOLIDS Matrix(i,j) for mol. weights of constituents of
                          all solid phases
     *     NPHAS        : total number of phases contained in data file
     *     NMIXPHAS     : Number of mixed phases in parameterfile
     *     NSCOM        : ... of system components in data file
     *     NPCON        : ... of constituents OF ALL PHASES
     *     PMOD         : Array of names of models for all phases
    *     SCNAME       : Array of NSCOM system component names
     *     PCNAME_GAS   : Vector(i) of phase constituents names for all
     *     PCNAME_AQ      constituents in gas- and aqueous phase
     *     PCNAME_SOLIDS  Matrix(i,j) of phase consituent names for all
                          solid phases
                          All these names are limited to a max. length
                          of 24 by chemapp
     *     PNAME        : Array for phase names contained in data file:
     *                    "gas", "aqueous", and all the mineral phases
     *     PMOLVOL_AQ   : Array for molar volumes of all constituents of
                          aqueous phase
     *     PMOLVOL_SOLIDS Matrix(i,j) ... all solid phases.
     *     PSTATNAME    : Name for phase stati, necessary for CHEMAPP
     *     FILE_DB      : File name of parameterfile. if not in same
     *                    directory like executables, directory must be
     *                    included.
     *     KIND_OF_PFILE: Two kinds of parameterfile are to be distin-
                          guished:
                          1)
                            SCNAME(1)(1:2) = 'EA'
                            SCNAME(2)(1:1) = 'H'
                            SCNAME(3)(1)   = 'O'
                          2)
                            SCNAME(1)(1:3) = 'H2O'
                            SCNAME(2)(1:4) = 'H<+>'
                          No other kinds are allowed in the present
                          version.
  ******************************************************/
// History of modifications
//   1.0.0. first implementaion

//Further variables *********************
      double ***PCSTOI_SOLIDS,**PCSTOI_GAS, **PCSTOI_AQ;   //stoi coefficient
      double **PCMMASS_SOLIDS,*PCMMASS_GAS,*PCMMASS_AQ;    //mol mass (g/mol)
      double **PMOLVOL_SOLIDS, *PMOLVOL_GAS, *PMOLVOL_AQ;  //molar volume (cc/mol)
      long   NPHAS, NMIXPHAS, NSCOM, KIND_OF_PFILE;
  public:
	  long *NPCON;
  private:
      char   *FILE_DB;//120
      char   **PMOD, **SCNAME, **PCNAME_GAS, **PCNAME_AQ, ***PCNAME_SOLIDS;  
      char   **PNAME, *PSTATNAME[4]; //24
      char   TUnit[3],PUnit[3],VUnit[3],MUnit[4],EUnit[3];
//      char   *TUnit, *PUnit,*VUnit,*MUnit,*EUnit;
//	  int	 verbosity;
   	  long   nodenumber;   // number of nodes, on which reactions are calculated 
	  int	 heatflag, flowflag;
	  int	 action;

/*     *******************************************************
     VERSION
 ### memEQLINK-s.h    : 1.0.0  (2006-04-03
     EDITOR           : Helge C. Moog
     *******************************************************
     * Short Description:
     * ------------------
     *     Variables for common block 'memEQLINK'. This common block
     *     serves to control or define calculations of chemapp.
     *
     * Calling Program(s):
     * -------------------
     *     initEQLINK, EQLINK, EQCALC, EQCHECK
     *
     * Variables        :
     * ------------------
 DOUBLE PRECISION

     memEQ_Time
      vector(i) for the point in time for which calculational
C       results are stored in memory. The unit of time is up to the
C       calling routine. It usually applies:
C       Unit = y (LOPOS) or = s (MARNIE)
C         i = number of node (iNode)
C
C     memEQ_pH
C       vector(i) to store calculated pH
C         i = number of node (iNode)
C
C     memEQ_molH2O
C       Vector(i) to store calculated moles of free water
C         i = number of node (iNode)
C
C     memEQ_kgH2O
C       vector(i) for kg H2O (free water)
C         i = number of node (iNode)
C
C     memEQ_mol_pcon_gas, memEQ_mol_pcon_aq
C       matrix(i,j) for moles of constituents in gas- and aqueous phase
C         i = number of node (iNode)
C         j = index number of phaseconstituent in the list
C
C     memEQ_mol_pcon_solids
C       matrix(i,j,k) for moles of constituents in all solid phases
C       phases.
C         i = number of node (iNode)
C         j = number of phase
C         k = index number of phaseconstituent in the list
C
C     memEQ_mol_scom_gas,
C     memEQ_mol_scom_aq,
C     memEQ_mol_scom_solids
C     memEQ_mol_scom_total
C       matrix(i,j) for moles of system components in the gas-, aqueous
C       and all solid phases and summed up for all phases.
C       i = number of node (iNode)
C       j = number of system component
C       Note that saving these vectors is redundant as this information
C       could be retrieved by adding up the system components contained
C       memEQ_mol_gas, memEQ_mol_aq, and memEQ_mol_solids. However,
C       preprocessing done by eqcheck is entirely in terms of system
C       components. Therefore, both kind of information is saved and
C       will reveal which shall be mainntained and which not.
C
C INTEGER
C
C     maxnSys
C       maximum number of allowed systems (=segments or nodes)
C       This variable may be removed with C++
C
C     maxn_iSysCalc_failed
C       maximum number of failed equilibrium calculations in a system
C       before the entire simulation is STOPped. This value is
C       determined by the calling program in subroutine initEQLINK.
C
C     maxNPHAS_exist
C       realistic total number of phases which might evolve within one
C       system during a complete run of the calling program.
C
C     maxNPCON_gas_exist
C       realistic number of phase constituents for the gaseous phase,
C       which might evolve within a given system cumulatively at all
C       times.
C
C     maxNPCON_aq_exist
C       ... same for the aq. phase.
C
C     maxNPCON_solids_exist
C       ... same for all solid phases. For stoichiometric solid phases
C       the number of phase constituents is 1. Only with mixed solid
C       phases this number is >1. At the moment it is assumed that
C       future mixed solid phases are not more complex as 5
C       constituents...
C
C     maxn_solids_exist
C       maximum number of solid phases which are likely to occure within
C       a run of the calling programm in a given system. It does not
C       mean the number of solid phases which may occure
C       ***simultaneously***, but also includes those with only
C       intermediate stability.
C
C     iTimeStep
C       vector(i) to count number of EQLINK-calls for any system
C         i = number of node (iNode)
C
C     iSysCalc
C       vector(i) which counts (successfull) calculations done in each
C       system
C         i = number of node (iNode)
C
C     iSysCalc_failed
C       vector(i) which counts when chemapp fails to calculate the
C       equilibrium for the present system
C         i = number of node (iNode)
C       When this counter attains a pre-determined value
C       (maxn_iSysCalc_failed) then the entire simulation is STOPped.
C
C     iSysCalc_repeated
C       vector(i) which counts how often a non-target calculation had
C       to be repeated as target calculation because in the former case
C       calculated system volume deviated too much from to its input
C       value.
C       If this value is getting too large one should consider to restart
C       the whole calculation with enforced target calculation
C       (option(3) = 1). In the present version EQCALC independently
C       switches to the target-mode (option(3) = 2) for the given system
C       as soon as this counter attains a value of 100 or more and
C       represents a proportion of >= 10% of all calculations.
C         i = number of node (iNode)
C
C     iSysPrevResults
C       vector(i) which counts return of previous results in each system
C         i = number of node (iNode)
C
C     iSeqCalc
C       vector(i) which counts the number of calculations within a
C       sequence.
C         i = number of node (iNode)
C
C     nSeqCalc
C       vector(i) for the number of calculations which are performed
C       before the next extrapolation. CHEMAPP is called at least
C       nSeqCalc times, before checking whether it is appropriate to
C       a new calculation. It must not <1. This variable may be
C       modified at runtime from within EQLINK-related subroutines.
C       In the present version no credit as taken from this possiblity.
C       In future versions it may be decided that nSeqCalc is increased
C       under certain conditions.
C         i = number of node (iNode)
C
C     iNPCON_gas_exist,iNPCON_aq_exist
C       vector(i) for the total number of gas- or aq. phase constituents
C       and solid phases which for a given system has been at any time
C       or is at the present time existing
C         i = number of node (iNode)
C
C     iNPCON_solids_exist
C       matrix(i,j) for the total number of phase constituents in all
C       solid phases which for a given system has been at any time
C       or is at the present time existing
C         i = number of node (iNode)
C         j = number of (solid) phase
C
C     gas_pcon_exist,aq_pcon_exist
C       matrix(i,j) to store a list of gas- or aqeous phase constituents
C       which for a given system has been at any time or is at the
C       present time existing
C         i = number of node (iNode)
C         j = index number in the list
C
C     solids_pcon_exist
C       matrix(i,j,k) to store a list of solid phases which for a given
C       system has been at any time or is at the present time existing
C         i = number of node (iNode)
C         j = index number of solid phase in the list
C         k = index number of solid phase consituent in the list
C
C     channel
C       channel to which all outputs of EQLINK-related subroutines
C       (incl. CHEMAPP itself!) shall be directed.
C       Fortran-speciality, necessary in C++?
C
C     verbosity
C       defines to which extent EQLINK-related subroutines deliver
C       messages with write-commands. It is set by the calling program.
C         = 0: no messages at all!
C         = 1: only error messages
C         = 2: errors and warnings
C         = 3: errors, warnings and notes
C         = 4: errors, warnings, notes, and debug infos
C             (= ALL write-commands are executed)
C
C     pstat
C       Array for status of all phases.
C       Index number = number of phase in parameterfile
C         1 = (ALWAYS) GAS
C         2 = (ALWAYS) aqueous solution
C         >3= all other phases
C       Three stati are possible:
C       pstat(i) = 1: entered (= active)
C                  2: dormant (= activity of phase is
C                                calculated, but it doesn't form
C                  3: eliminated (= phase is excluded from calculation
*/
//*******************************************************
//History of modifications
// 1.0.0.
//   The precursor for this file was 1.1.6
//C Define common block
//      COMMON /memEQLINK/
     double *memEQ_Time,
     *memEQ_pH,
     *memEQ_molH2O,
     *memEQ_kgH2O,
     **memEQ_mol_pcon_gas,
     **memEQ_mol_pcon_aq,
     ***memEQ_mol_pcon_solids,
     **memEQ_mol_scom_gas,
     **memEQ_mol_scom_aq,
     **memEQ_mol_scom_solids,
     **memEQ_mol_scom_total,
	 *mol_scom_total_save,
	 *sysT0,*sysT,
	 *sysP,
	 *sysV;
     int *iTimeStep, 
		 *iSysCalc,  //??
           *iSysCalc_failed,
           *iSysCalc_repeated,
      *iSysPrevResults,
      *iSeqCalc,  //for future extrapolation
      *nSeqCalc,
      *iNPCON_gas_exist,
      *iNPCON_aq_exist,
      **iNPCON_solids_exist,
      **gas_pcon_exist,
      **aq_pcon_exist,
      ***solids_pcon_exist,
      *pstat;
	  int *option; //option[5];
	  int *scom_elim_index;   //channel,verbosity,

    //Parameters for EQCALC
	 double *mol_scom_total,
            *molH2O_in, 
            *molH2O_out, *kgH2O_out,
            *pH_out, *Eh_out,
            *sysP_out, *sysV_out,
            **PVOL_out, 
            **PMASS_out,
            **mol_pcon_gas_out,
            **mol_pcon_aq_out,
            ***mol_pcon_solids_out,
            **mol_scom_gas_out,
            **mol_scom_aq_out,
            **mol_scom_solids_out;
	 int iPrevActions,
		 PrevActions[4];
	  

      
  public:
	bool flag_chemapp;     // flag if *.chm file exists 
	int maxn_iSysCalc_failed;

  public:
    CEqlink(void);
    ~CEqlink(void);
  
  void TestCHEMAPPParameterFile(string); //Test file exist
  void prepare(string filename);   //reading the database of CHEMAPP 
  CEqlink * CEqlink::GetREACTION(void);
  void GetPconNameAq(void);
  double GetPconAq_mol_amount(const long, const long);
  void initEQLINK(void);
  void DestroyMemory();
  void callCHEMAPP(string filename);

  //++++++++++ioEQLINK++++++++++++++
  void AllocateMemory_IO(void);
  void DeleteMemory_IO(void);
  void ExecuteEQLINK();
  void callEQCHECK(long node);
  void callEQCALC(long node);
  void callEQGETPRE(long node);
  void SetResultsBackMassTransport(void);

  //++++++++EQCalc++++++++++++++++++
  void AllocateMemory_EQCalc(void);
  void DeleteMemory_EQCalc(void);

};
extern vector <CEqlink*> Eqlink_vec;
extern vector <string> PconName_Aq;
extern void CHMRead(string filename);

#endif