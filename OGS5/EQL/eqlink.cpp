#ifdef CHEMAPP

// #include "stdafx.h" /* MFC */


#include "eqlink.h"
#include "files0.h"
#include <iostream>
#include <fstream>
//#include "geo_strings.h"
#include "cacint.h"
#include <cstring>
#include "rf_pcs.h"
#include "rfmat_cp.h"

vector <CEqlink*> Eqlink_vec;
vector <string> PconName_Aq;
vector <double>Pcon_Aq_mol_vec;

#define no_DEBUG_CHEMAPP  //for debug
 static double lambda=0.001;

/**************************************************************************
FEMLib - Object: 
Task: 
Programing:
01/2006 MX Implementation
last modified
**************************************************************************/
CEqlink::CEqlink(void){
//	CRFProcess * m_pcs = NULL;
	flag_chemapp = false;
	//Save level of verbosity, 
//	verbosity = 5;
	maxn_iSysCalc_failed = 100; //const static
	nodenumber = 1;
	heatflag = 0;
	flowflag = 0;
	action = 1;
    option = new int [5];
	option[1] = 2;  //gas, 1=ENTERED, 2=DORMANT, 3=ELIMINATED
	option[2] = 1;  //enforced new calculation
	option[3] = 0;  //system volume is not considered a target
	option[4] = 0;  //No estimate allowed

}
/* Destructor */
CEqlink::~CEqlink(void){
    delete option;
    option=NULL;
}

/**************************************************************************
FEMLib - Object: 
Task: 
Programing:
04/2006 MX Implementation
last modified
**************************************************************************/
void CEqlink::DestroyMemory(){

	string str;
	int i,j;

//    delete [] TUnit;
//    delete [] PUnit;
//    delete [] VUnit;
//    delete [] MUnit;
//    delete [] EUnit;


	for (i=0; i<NPHAS; i++){
	  for (j=0; j<NPCON[i]; j++){
	    str = PNAME[i];
  	    if (str.find("GAS")!=string::npos){
  		  if (PCNAME_GAS) delete []PCNAME_GAS[j];
	    } else if (str.find("AQUEOUS")!=string::npos){
		  if (PCNAME_AQ) delete [] PCNAME_AQ[j];
	    } else {
		  if (PCNAME_SOLIDS[i]) 
			  delete []PCNAME_SOLIDS[i][j];
	    }


  	    if (str.find("GAS")!=string::npos){
		  if (PCSTOI_GAS[j]) delete [] PCSTOI_GAS[j];
	    } 
		else if (str.find("AQUEOUS")!=string::npos){
		  if (PCSTOI_AQ[j]) delete [] PCSTOI_AQ[j];   //[j] = new double[NSCOM];
	    } 
		else {
		  if (PCSTOI_SOLIDS[i][j]) delete []PCSTOI_SOLIDS[i][j]; //[j]=new double[NSCOM];
	    }
		

	  }  //j =NPCON[i]

	  //Get names of models for each phase
	  if (str.find("GAS")!=string::npos){
		delete []PCNAME_GAS; // = new char *[NPCON[i]];
		delete []PCMMASS_GAS; // = new double [NPCON[i]];
		delete []PCSTOI_GAS; // = new double *[NPCON[i]];
		delete []PMOLVOL_GAS; // = new double [NPCON[i]];
	  } else if (str.find("AQUEOUS")!=string::npos){
		delete []PCNAME_AQ; // = new char *[NPCON[i]];
		delete []PCMMASS_AQ; // = new double [NPCON[i]];
		delete []PCSTOI_AQ; // = new double *[NPCON[i]];
		delete []PMOLVOL_AQ; // = new double [NPCON[i]];
	  } else {
		delete []PCNAME_SOLIDS[i]; // = new char *[NPCON[i]];
		delete []PCMMASS_SOLIDS[i]; // = new double [NPCON[i]];
		delete []PCSTOI_SOLIDS[i]; // = new double* [NPCON[i]];
		delete []PMOLVOL_SOLIDS[i]; // = new double [NPCON[i]];
	  }
    }
	delete []NPCON; //  = new long [NPHAS];
//	NPCON =NULL;
	delete []PCNAME_SOLIDS; //  = new char **[NPHAS];
	PCNAME_SOLIDS =NULL;
	delete []PCMMASS_SOLIDS; //  = new double *[NPHAS];
	PCMMASS_SOLIDS =NULL;
	delete []PCSTOI_SOLIDS; //  = new double**[NPHAS];
	PCSTOI_SOLIDS =NULL;
	delete []PMOLVOL_SOLIDS; // = new double *[NPHAS];
	PMOLVOL_SOLIDS =NULL;

	for (int i=0; i<NSCOM; i++){
		delete []SCNAME[i];
	}
	delete[] SCNAME;
	SCNAME=NULL;
	for (int i=0; i<NPHAS; i++){
		delete []PMOD[i];
		delete []PNAME[i];
	}
	delete []PMOD;
	PMOD=NULL;
	delete []PNAME;
	PNAME=NULL;

	DeleteMemory_IO();
    DeleteMemory_EQCalc();  //test MX

  for(int i=0;i<(int)Eqlink_vec.size();i++){
    delete Eqlink_vec[i];
  }
  Eqlink_vec.clear();

 
}


/**************************************************************************
FEMLib - Object: 
Task: 
Programing:
01/2006 MX Implementation
last modified
**************************************************************************/
void CEqlink::TestCHEMAPPParameterFile(string file_base_name){
 
  string chm_file_name = file_base_name + CHEMAPP_REACTION_EXTENSION;
  ifstream chm_file (chm_file_name.data(),ios::in);
  if (chm_file.good()){ 
	flag_chemapp = true;
	chm_file.close();
  }
  else {
	DisplayErrorMsg("!!! Warning: no  CHEMAPP datafile *.chm !!!");
    flag_chemapp = false;
	chm_file.close();
  }
}

/**************************************************************************
FEMLib - Object: 
Task: 
Programing:
01/2006 MX Implementation
last modified
**************************************************************************/
void CHMRead(string filename){

  CEqlink *eq = new CEqlink(); 
  eq->TestCHEMAPPParameterFile(filename); // Test if *.chm file is present
  if (eq->flag_chemapp){
    eq->flag_chemapp = true;
    Eqlink_vec.push_back(eq);
  }
  else { 
    delete eq;
    eq=NULL;
  }
}

/**************************************************************************
  GeoSys - Funktion: CEqlink::GetREACTION(void)
 
   Aufgabe:
     Call EQLINK	
  
  Ergebnis:
   
                                                                         
  Programmaenderungen:
     064/2006   MX  Erste Version
**************************************************************************/
CEqlink * CEqlink::GetREACTION(void)
{
  CEqlink *eq = NULL;
  /* Tests */
  if(Eqlink_vec.capacity() > 0)
	eq = Eqlink_vec[0];
  return eq;
}

/**************************************************************************
  GeoSys - Funktion: CEqlink::GetPconNameAq(void)
 
   Aufgabe:
     Get PCON_aq name vector	
  
  Ergebnis:
   
                                                                         
  Programmaenderungen:
     08/2008   MX  Erste Version
**************************************************************************/
void CEqlink::GetPconNameAq(void)
{
  int i;
	CEqlink *eq = NULL;
  string pconName_aq;
  eq = GetREACTION();
  for(i=0;i<NPCON[1];i++){
    pconName_aq = eq->PCNAME_AQ[i];
    PconName_Aq.push_back(pconName_aq);
  }
}

/**************************************************************************
  GeoSys - Funktion: CEqlink::GetPconNameAq(void)
 
   Aufgabe:
     Get PCON_aq mol vector	
  
  Ergebnis:
   
                                                                         
  Programmaenderungen:
     08/2008   MX  Erste Version
**************************************************************************/
double CEqlink::GetPconAq_mol_amount(const long node, const long pcon_indx)
{
  CEqlink *eq = NULL;
  eq = GetREACTION();

  return eq->mol_pcon_aq_out[node][pcon_indx];
}

/**************************************************************************
  GeoSys - Funktion: CEqlink::callCHEMAPP(void)
 
   Aufgabe:
     Call EQLINK	
  
  Ergebnis:
   
                                                                         
  Programmaenderungen:
     064/2006   MX  Erste Version
**************************************************************************/
void CEqlink::callCHEMAPP(string filename){

	prepare(filename);
    //AllocateMemory_IO();
	//initEQLINK();
    //AllocateMemory_EQCalc();
	ExecuteEQLINK();
}
/**************************************************************************
FEMLib - Object: 
Task: 
Programing:
01/2006 MX Implementation based on Helge Moog prepare.f
last modified
**************************************************************************/
void CEqlink::prepare(string file_base_name){

	long int ivers, noerr;
    long int  unitno=0;
	long i, j, k;
	char *FILE_DB;
	char dstr[25];
	double temp[30],dtemp, Tdtemp;
	string str;
	char strpcon[25];
	FILE_DB = new char[20];
//	TUnit =new char[3];
//    PUnit =new char[3];
//    VUnit =new char[3];
//    MUnit =new char[4];
//    EUnit =new char[3];

    string chm_file_name = file_base_name + CHEMAPP_REACTION_EXTENSION;
    strcpy(FILE_DB, chm_file_name.c_str());

	//Initialise ChemApp
	tqini(&noerr);

  /* After all following invocations of ChemApp subroutines, the error
     variable will be checked for a non-zero value, in which case an 
     error has occurred and the program is aborted */
	if (noerr) {
		DisplayErrorMsg("!!! Fehler: Initializing CHEMAPP !!!");
		abort();
	} 

  /* Get version number of ChemApp*/
    tqvers(&ivers,&noerr);
	if (noerr) {
		DisplayErrorMsg("!!! Fehler: Initializing CHEMAPP !!!");
		abort();
	} 

	tqgio("FILE", &unitno, &noerr);
	if (noerr) {
		DisplayErrorMsg("!!! Fehler: tqgio by CHEMAPP !!!");
		abort();
	} 

	tqopna(FILE_DB,unitno,&noerr);
	if (noerr) {
		DisplayErrorMsg("!!! Fehler: Initializing CHEMAPP by tqopna !!!");
		abort();
	} 
  /* Read data-file */
    tqrfil(&noerr);
	if (noerr) {
		DisplayErrorMsg("!!! Fehler: Initializing CHEMAPP in tqrfil !!!");
		abort();
	} 
	tqclos(unitno,&noerr);

	  /* Get system units */
	tqgsu("Pressure",dstr,&noerr);
	tqgsu("Volume",dstr,&noerr);
	tqgsu("Temperature",dstr,&noerr);
	tqgsu("Amount",dstr,&noerr);

	/*Get the number of system components*/
	tqnosc(&NSCOM, &noerr);

	//Get the names of all system components
	SCNAME = new char *[NSCOM];
	for (i=0; i<NSCOM; i++){
		SCNAME[i]=new char [25];
		tqgnsc(i+1, SCNAME[i], &noerr);
	}

    //Get total number of phases contained in the parameterfile
	tqnop(&NPHAS, &noerr);
	#ifdef _DEBUG_CHEMAPP
		printf("  Number of system components =  %ld \n", NSCOM);
		printf("  Number of phases            =  %ld \n", NPHAS);
	#endif

	//Set amount unit to gram, so we get the molecular mass in g/mol
	tqcsu("Amount", "gram", &noerr);
	if (noerr) {
		printf("  !!! Fehler: ChemApp reports error no. %lf \n", noerr);
		DisplayErrorMsg("Possible cause: only zeros in stoichiometric");
		DisplayErrorMsg(" matrix of a phase consituent");
		abort();
	} 

	NMIXPHAS = 0;

	//Collect information for all phases in parameterfile
	PNAME = new char *[NPHAS];  //phase name
	PMOD  = new char *[NPHAS];  //names of models for each phase
	NPCON = new long [NPHAS];
	PCNAME_SOLIDS = new char **[NPHAS];
	PCMMASS_SOLIDS = new double *[NPHAS];
	PCSTOI_SOLIDS = new double**[NPHAS];
	PMOLVOL_SOLIDS= new double *[NPHAS];
	for (i=0; i<NPHAS; i++){
	  PNAME[i] = new char [25];
	  tqgnp(i+1, PNAME[i], &noerr);


	  //Get names of models for each phase
	  PMOD[i] = new char [25];
	  tqmodl(i+1, PMOD[i],&noerr);

	  //Get number of phase constituents for each phase
	  tqnopc(i+1,&NPCON[i], &noerr);
	  if (NPCON[i]>1) NMIXPHAS=NMIXPHAS+1;  // count mixed phases

	  str = PNAME[i];
	  if (str.find("GAS")!=string::npos){
		PCNAME_GAS = new char *[NPCON[i]];
		PCMMASS_GAS = new double [NPCON[i]];
		PCSTOI_GAS = new double *[NPCON[i]];
		PMOLVOL_GAS = new double [NPCON[i]];
	  } else if (str.find("AQUEOUS")!=string::npos){
		PCNAME_AQ = new char *[NPCON[i]];
		PCMMASS_AQ = new double [NPCON[i]];
		PCSTOI_AQ = new double *[NPCON[i]];
		PMOLVOL_AQ = new double [NPCON[i]];
	  } else {
		PCNAME_SOLIDS[i] = new char *[NPCON[i]];
		PCMMASS_SOLIDS[i] = new double [NPCON[i]];
		PCSTOI_SOLIDS[i] = new double* [NPCON[i]];
		PMOLVOL_SOLIDS[i] = new double [NPCON[i]];
	  }

	  //Get names of phase constituents, their composition with
	  //respect to system components and their molecular mass.
	  for (j=0; j<NPCON[i]; j++){
		tqgnpc(i+1, j+1, strpcon, &noerr);
  	    if (str.find("GAS")!=string::npos){
  		  PCNAME_GAS[j] = new char [25];
		  strcpy(PCNAME_GAS[j], strpcon);
	    } else if (str.find("AQUEOUS")!=string::npos){
		  PCNAME_AQ[j] = new char [25];
		  strcpy(PCNAME_AQ[j], strpcon);
	    } else {
		  PCNAME_SOLIDS[i][j] = new char [35];
		  strcpy(PCNAME_SOLIDS[i][j], strpcon);
	    }


		tqstpc(i+1, j+1, temp, &dtemp, &noerr);
  	    if (str.find("GAS")!=string::npos){
  		  PCMMASS_GAS[j] = dtemp;
		  PCSTOI_GAS[j] = new double[NSCOM];
		  for (k=0; k<NSCOM;k++)
			PCSTOI_GAS[j][k] = temp[k];
	    } 
		else if (str.find("AQUEOUS")!=string::npos){
		  PCMMASS_AQ[j] = dtemp;
		  PCSTOI_AQ[j] = new double[NSCOM];
		  for (k=0; k<NSCOM;k++)
			PCSTOI_AQ[j][k] = temp[k];
	    } 
		else {
		  PCMMASS_SOLIDS[i][j] = dtemp;
		  PCSTOI_SOLIDS[i][j]=new double[NSCOM];
		  for (k=0; k<NSCOM;k++)
			PCSTOI_SOLIDS[i][j][k]=temp[k];
	    }
		

		//Get molar volume of phase constituents
		tqgdpc("V",i+1,j+1,&dtemp,&noerr);
		tqgdpc("T",i+1,j+1,&Tdtemp,&noerr);
//		PMOLVOL_SOLIDS[i][j] =dtemp*R*1.0E+06; 
		if (str.find("GAS")!=string::npos){
  		  PMOLVOL_GAS[j] = dtemp*R_GAS*Tdtemp;  //TODO, test, l/mol
	    } 
		else if (str.find("AQUEOUS")!=string::npos){
		  PMOLVOL_AQ[j] = dtemp*R_GAS*1.0E+06;
	    } 
		else {
		  PMOLVOL_SOLIDS[i][j] =dtemp*R_GAS*1.0E+06;  //cc/mol
	    }  //if
	  }  //j =NPCON[i]
	}    //i=NPHAS
	  //Set amount unit back to moles for calculations
	tqcsu("Amount", "mol", &noerr);
	tqcsu("Volume", "m3", &noerr);
	tqcsu("Pressure", "Pa", &noerr);
	tqsetc("P", 0, 0, 1.0e+05, &i, &noerr);

	//Get units for calculation results of chemapp
	tqgsu("Temperature", TUnit, &noerr);
	tqgsu("Pressure", PUnit, &noerr);
	tqgsu("Volume", VUnit, &noerr);
	tqgsu("Energy", EUnit, &noerr);
	tqgsu("Amount", MUnit, &noerr);

/* Initialize array for allowed names of phase stati.
  Necessary to be able to change phase stati in initEQLINK and NATHAN
  because the corresponding CHEMAPP-command expects these strings as
  parameters.*/
	for (i=1;i<4;i++)
      PSTATNAME[i] = new char [15];
	strcpy(PSTATNAME[1], "ENTERED");
	strcpy(PSTATNAME[2],"DORMANT");
    strcpy(PSTATNAME[3],"ELIMINATED");

	//Now check parameterfile
	if (strcmp(PCNAME_AQ[0],"H2O")!=0 || strcmp(PCNAME_AQ[1],"H<+>")!=0 || strcmp(PCNAME_AQ[2],"OH<->")!=0) {
		DisplayErrorMsg("Fehler: Illegal names in parameterfile !!!");
		DisplayErrorMsg("The first 3 phase constituents of aqueous");
		DisplayErrorMsg("phase must be: H2O, H<+>, OH<->");
		abort();
	}
	KIND_OF_PFILE = 0;
	if (strcmp(SCNAME[0],"EA")==0 && strcmp(SCNAME[1],"H")==0 && strcmp(SCNAME[2],"O")==0 )
	  KIND_OF_PFILE = 1;
	if (strcmp(SCNAME[0],"H2O")==0 && strcmp(SCNAME[1],"H<+>")==0)
	  KIND_OF_PFILE = 2;
	if (KIND_OF_PFILE == 0){
		DisplayErrorMsg("PREPARE (Fatal Error), Illegal names in parameterfile.");
		DisplayErrorMsg("Two possible designators for system");
		DisplayErrorMsg("1. EA, H, O or 2. H2O, H<+>");
		printf("Current: %s,%s,%s",SCNAME[0],SCNAME[1],SCNAME[2]);
		abort();
	}
}

/**************************************************************************
  GeoSys - Funktion: CEqlink::initEQLINK(void)
 
   Aufgabe:
     Call EQLINK	
  
  Ergebnis:
   
                                                                         
  Programmaenderungen:
     064/2006   MX  Erste Version
**************************************************************************/
void CEqlink::initEQLINK(void){
  
  int i,j,k,m,n;
  long noerr=0;
  CRFProcess *m_pcs = NULL;
  int indx0, indx1;
  double temperature=0.0;
  int timelevel = 1;
  int comp =0;
  double temp=0.0;
  string str,str1;

  if (aktueller_zeitschritt==0){
  // heat transport
  for(i=0;i<(int)pcs_vector.size();i++){
	m_pcs = pcs_vector[i];
	if (m_pcs->pcs_type_name.compare("HEAT_TRANSPORT")==0){ 
		heatflag = 1;
		break;
	}
  }
  // heat transport
  for(i=0;i<(int)pcs_vector.size();i++){
	m_pcs = pcs_vector[i];
	if (m_pcs->pcs_type_name.compare("GROUNDWATER_FLOW")==0){ 
		flowflag = 1;
		break;
	}
	else if (m_pcs->pcs_type_name.compare("LIQUID_FLOW")==0){ 
		flowflag = 2;
		break;
	}
	else if (m_pcs->pcs_type_name.compare("RICHARDS_FLOW")==0){ 
		flowflag = 3;
		break;
	}
	else if (m_pcs->pcs_type_name.find("FLOW")==0){ 
		flowflag = 4;
		break;
	}
  }

  
  // read number of nodes
  for(i=0;i<(int)pcs_vector.size();i++){
	m_pcs = pcs_vector[i];
	if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT")==0){ 
		nodenumber = (long) m_pcs->m_msh->GetNodesNumber(false);
		break;
	}
  }
  if (aktueller_zeitschritt==0) AllocateMemory_IO();   //Moved to CallEQLINK, only once
  for(i=0;i<nodenumber;i++){
	nSeqCalc[i] = 1;
  }
  for (i=0;i<NPHAS;i++){
	pstat[i] = 1;  //!By default, set all phases 'ENTERED' (=1)
  }

  #ifdef _DEBUG_CHEMAPP
    DisplayErrorMsg("Physic unit defined after initialization");
	
//	tqgsu("Temperature", TUnit, &noerr);
	printf("\n  Temperature   %s \n", TUnit);
	printf("  Pressure      %s \n", PUnit);
	printf("  Volume        %s \n", VUnit);
	printf("  Mass          %s \n", MUnit);
    printf("  Energy        %s \n", EUnit);
  #endif
 }
//Get temperature value 
 if(heatflag > 0){
     m_pcs = PCSGet("HEAT_TRANSPORT");
     timelevel =0;
     indx0 = m_pcs->GetNodeValueIndex("TEMPERATURE1")+timelevel;
     timelevel =1;
     indx1 = indx0 +1;  //m_pcs->GetNodeValueIndex("TEMPERATURE1")+timelevel;

	 for(i=0;i<this->nodenumber;i++){
		sysT0[i] = m_pcs->GetNodeValue(i, indx0); 
        sysT[i] = m_pcs->GetNodeValue(i, indx1); 
		if (sysT0[i] <273.15) sysT0[i] += 273.15;  //ToDo °C->K
        if (sysT[i] <273.15) sysT[i] += 273.15;  //ToDo °C->K
	 }
  }
  else {
    for(i=0;i<this->nodenumber;i++)
      sysT0[i] = sysT[i] = 298.15;
  }

  //Get pressure value 
  if(flowflag > 0){
	  switch(flowflag){
        case 1:
		  m_pcs = PCSGet("GROUNDWATER_FLOW");
		  indx0 = m_pcs->GetNodeValueIndex("HEAD")+timelevel;
		  for(i=0;i<this->nodenumber;i++){
//MX test			sysP[i] = m_pcs->GetNodeValue(i, indx0);
			sysP[i] = 1.0e5;  //m_pcs->GetNodeValue(i, indx0); 
			if (sysP[i]<1.0e5) sysP[i] = 1.0e5;
		  }
          break;
        case 2:
		  m_pcs = PCSGet("LIQUID_FLOW");
		  indx0 = m_pcs->GetNodeValueIndex("PRESSURE1")+timelevel;
		  for(i=0;i<this->nodenumber;i++){
			sysP[i] = m_pcs->GetNodeValue(i, indx0); 
		  	if (sysP[i]<1.0e5) sysP[i] = 1.0e5;
		  }
          break;
        case 3:
		  m_pcs = PCSGet("RICHARDS_FLOW");
		  indx0 = m_pcs->GetNodeValueIndex("PRESSURE1")+timelevel;
		  for(i=0;i<this->nodenumber;i++){
			sysP[i] = m_pcs->GetNodeValue(i, indx0);
			if (sysP[i]<1.0e5) sysP[i] = 1.0e5;
		  }
          break;
        case 4:
		  DisplayErrorMsg("Error: Not implemented for the flow in CHEMAPP case!!!");
		  sysP[i] = 1.0e5;   //TODO  MX
          break;
      }
	}
	else {
	  DisplayErrorMsg("Warning: No valid flow process!!");
	  exit(1);
	}

	//Get value from mass transport
	if (KIND_OF_PFILE == 1) {

    for(i=0;i<nodenumber;i++){
		for (m =1; m<NSCOM;m++){
		  memEQ_mol_scom_solids[i][m]=0.0;
		}
	}

	  for(j=0;j<(int)pcs_vector.size();j++){ // for all processes
		m_pcs = pcs_vector[j];
		if(m_pcs->pcs_type_name.compare("MASS_TRANSPORT")==0){ // if it is a mass transport process
			comp = m_pcs->pcs_component_number; // get component number
			str=m_pcs->pcs_primary_function_name[0];
			if (cp_vec[comp]->mobil){
				for (k=0; k<NSCOM;k++){   //pH  -> H+
//				    if (strcmp(m_pcs->pcs_primary_function_name[0],"pH") ==0 && strcmp(SCNAME[k],"H") ==0){
				    if (strcmp(m_pcs->pcs_primary_function_name[0],"pH") ==0){
					  for(i=0;i<nodenumber;i++){
						memEQ_mol_scom_aq[i][k] =  \
							m_pcs->GetNodeValue(i,m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) \
							+timelevel);
						memEQ_mol_scom_aq[i][k] += 0.0;  //pow(10.0, -memEQ_mol_scom_aq[i][k]);
					  } //for i
					  break;
					}  //if
					else if (strcmp(m_pcs->pcs_primary_function_name[0],SCNAME[k]) ==0) {
					  for(i=0;i<nodenumber;i++){
						memEQ_mol_scom_aq[i][k] =  \
							m_pcs->GetNodeValue(i,m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) \
							+timelevel);
					  }
					  break;
					}  //else if
				}  //for k
				for(i=0;i<nodenumber;i++){
				// get concentration values
				  if((memEQ_mol_scom_aq[i][j] < 0.0)  ) {  //MX, 01.2005 val_in[comp][i]&& (strcmp(name[comp],"pe")!= 0)
					if(abs(memEQ_mol_scom_aq[i][j]) > MKleinsteZahl) {
					  DisplayMsg(" Neg. conc for component "); DisplayLong((long) comp);
					  DisplayMsg(" at node "); DisplayLong((long)i);
					  DisplayMsg("; conc = "); DisplayDouble(memEQ_mol_scom_aq[i][j],0,0);
					  DisplayMsgLn(" ");
				    }
				    memEQ_mol_scom_aq[i][j] = 0.0 ;
				  }
				}  //for i
			}  //if mobile
			else if (cp_vec[comp]->mobil==0){ //immobile
				for (k=2; k<NPHAS;k++){  //0 -gas, 1-aq
					for (n =0;n<NPCON[k]; n++){
					  str1 = PCNAME_SOLIDS[k][n];
					  if (str1.find(str)!=string::npos){
						for(i=0;i<nodenumber;i++){
						    temp =  \
							  m_pcs->GetNodeValue(i,m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) \
							  +timelevel);
						    for (m =1; m<NSCOM;m++){
							  memEQ_mol_scom_solids[i][m]+= temp *PCSTOI_SOLIDS[k][n][m];
							}
						}
					    break;
					  } //if
					}  //for n
				}	//for k		
			}//if immobile
		}//if pcs = mass transport
	  }	//processes
	}  //KIND_OF_PFILE ==1
	else if (KIND_OF_PFILE == 2) {

      for(i=0;i<nodenumber;i++){
		for (m =0; m<NSCOM;m++){
		  memEQ_mol_scom_gas[i][m]= 0.0;
          memEQ_mol_scom_aq[i][m]=0.0;
          memEQ_mol_scom_solids[i][m]=0.0;
		}
	  }

	  for(j=0;j<(int)pcs_vector.size();j++){ // for all processes
		m_pcs = pcs_vector[j];
		if(m_pcs->pcs_type_name.compare("MASS_TRANSPORT")==0){ // if it is a mass transport process
			comp = m_pcs->pcs_component_number; // get component number
			str=m_pcs->pcs_primary_function_name[0];
			if (cp_vec[comp]->mobil){
				for (k=0; k<NPCON[1];k++){   //pH  -> H+
				    if (strcmp(m_pcs->pcs_primary_function_name[0],"pH") ==0 && strcmp(PCNAME_AQ[k],"H") ==0){
					  for(i=0;i<nodenumber;i++){
						memEQ_mol_pcon_aq[i][k] =  \
							m_pcs->GetNodeValue(i,m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) \
							+timelevel);
						memEQ_mol_pcon_aq[i][k] = 0.0;  //pow(10.0, -memEQ_mol_scom_aq[i][k]);
					  } //for i
					  break;
					}  //if
					else if (strcmp(m_pcs->pcs_primary_function_name[0],PCNAME_AQ[k]) ==0) {
					  for(i=0;i<nodenumber;i++){
						memEQ_mol_pcon_aq[i][k] =  \
							m_pcs->GetNodeValue(i,m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) \
							+timelevel);
                        for (m =0; m<NSCOM;m++){
							  memEQ_mol_scom_aq[i][m]+= memEQ_mol_pcon_aq[i][k]*PCSTOI_AQ[k][m];
					    }
					  }
					  break;
					}  //else if
				}  //for k
				for(i=0;i<nodenumber;i++){
				// get concentration values
				  if((memEQ_mol_pcon_aq[i][j] < 0.0)  ) {  //MX, 01.2005 val_in[comp][i]&& (strcmp(name[comp],"pe")!= 0)
					if(abs(memEQ_mol_scom_aq[i][j]) > MKleinsteZahl) {
					  DisplayMsg(" Neg. conc for component "); DisplayLong((long) comp);
					  DisplayMsg(" at node "); DisplayLong((long)i);
					  DisplayMsg("; conc = "); DisplayDouble(memEQ_mol_scom_aq[i][j],0,0);
					  DisplayMsgLn(" ");
				    }
				    memEQ_mol_pcon_aq[i][j] = 0.0 ;
				  }
				}  //for i
			}  //if mobile
			else if (cp_vec[comp]->mobil==0){ //immobile
				for (k=2; k<NPHAS;k++){  //0 -gas, 1-aq
					for (n =0;n<NPCON[k]; n++){
					  str1 = PCNAME_SOLIDS[k][n];
//For debug only                        if (str.find("MgSO4_OW")!=string::npos)
//                          k=k;
					  if (str1.find(str)!=string::npos){
						for(i=0;i<nodenumber;i++){
						    temp =  \
							  m_pcs->GetNodeValue(i,m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) \
							  +timelevel);
                            memEQ_mol_pcon_solids[i][k-2][n]= temp;  //unit [mol/l soil]
						    for (m =0; m<NSCOM;m++){
							  memEQ_mol_scom_solids[i][m]+= temp *PCSTOI_SOLIDS[k][n][m];
							}
						}
					    break;
					  } //if
					}  //for n
				}	//for k		
			}//if immobile
		}//if pcs = mass transport
	  }	//processes
	  //DisplayErrorMsg("Warning: The function for current chm file is not implemented now!!");
	 // exit(1);
	}
	else {
	  DisplayErrorMsg("Warning: No valid database file (.chm)!!");
	  exit(1);
	}

	// water 
    if (KIND_OF_PFILE == 1)
	  for(i=0;i<nodenumber;i++) molH2O_in[i]=55.5083491;  //TODO MX, 1kg water, gas?

	//total system component amount in each node
	for(i=0;i<nodenumber;i++){
	  for (j=0; j<NSCOM;j++) 
		memEQ_mol_scom_total[i][j]=memEQ_mol_scom_gas[i][j]  \
							  + memEQ_mol_scom_aq[i][j]   \
							  + memEQ_mol_scom_solids[i][j];
	}
}


/**************************************************************************
  GeoSys - Funktion: CEqlink::AllocateMemory_IO(void)
 
   Aufgabe:
     Call EQLINK	
  
  Ergebnis:
   
                                                                         
  Programmaenderungen:
     064/2006   MX  Erste Version
**************************************************************************/
void CEqlink::AllocateMemory_IO(void){
    int i,j;
	iTimeStep = new int[nodenumber]();
	iSysCalc  = new int[nodenumber]();
	iSysCalc_repeated = new int[nodenumber]();
	iSysPrevResults = new int[nodenumber]();
	iSeqCalc = new int[nodenumber]();
	nSeqCalc = new int[nodenumber]();
	iNPCON_gas_exist = new int[nodenumber]();

	memEQ_mol_pcon_gas = new double *[nodenumber];
	gas_pcon_exist = new int *[nodenumber];
	for (i=0;i<nodenumber;i++){
		memEQ_mol_pcon_gas[i] = new double [NPCON[0]]();
		gas_pcon_exist[i] = new int [NPCON[0]]();
	}

	iNPCON_aq_exist = new int[nodenumber]();
	memEQ_mol_pcon_aq = new double *[nodenumber];
	aq_pcon_exist = new int *[nodenumber];
	for (i=0;i<nodenumber;i++){
		memEQ_mol_pcon_aq[i] = new double [NPCON[1]]();
		aq_pcon_exist[i] = new int [NPCON[1]]();
	}

	iNPCON_solids_exist = new int *[nodenumber];
	solids_pcon_exist = new int **[nodenumber];
	memEQ_mol_pcon_solids = new double **[nodenumber];
	for (i=0;i<nodenumber;i++){
		iNPCON_solids_exist[i] = new int [NPHAS-2]();   //be careful with the NPHAS
		solids_pcon_exist[i] = new int *[NPHAS-2];
		memEQ_mol_pcon_solids[i] = new double *[NPHAS-2];
		for (j=0;j<NPHAS-2;j++){
		  solids_pcon_exist[i][j] =new int [NPCON[j+2]]();
		  memEQ_mol_pcon_solids[i][j] =new double [NPCON[j+2]]();
		}
	}

	memEQ_Time = new double [nodenumber]();
	memEQ_pH = new double [nodenumber]();
	memEQ_molH2O = new double [nodenumber]();
	memEQ_kgH2O = new double [nodenumber]();

	memEQ_mol_scom_gas = new double *[nodenumber];
	memEQ_mol_scom_aq = new double *[nodenumber];
	memEQ_mol_scom_solids = new double *[nodenumber];
	memEQ_mol_scom_total = new double *[nodenumber];
	for (i=0;i<nodenumber;i++){
		memEQ_mol_scom_gas[i] = new double [NSCOM]();
		memEQ_mol_scom_aq[i] = new double [NSCOM]();
		memEQ_mol_scom_solids[i] = new double [NSCOM]();
		memEQ_mol_scom_total[i] = new double [NSCOM]();
	}

	pstat = new int [NPHAS]();
	
	scom_elim_index = new int [NSCOM]();
	mol_scom_total_save = new double [NSCOM]();
	mol_scom_total = new double [NSCOM]();  //??? **
	molH2O_in = new double [nodenumber]();
	sysT0 = new double [nodenumber]();
    sysT = new double [nodenumber]();
	sysP = new double [nodenumber]();
	sysV = new double [nodenumber]();
}

/**************************************************************************
  GeoSys - Funktion: CEqlink::DeleteMemory_IO(void)
 
   Aufgabe:
     Call EQLINK	
  
  Ergebnis:
   
                                                                         
  Programmaenderungen:
     064/2006   MX  Erste Version
**************************************************************************/
void CEqlink::DeleteMemory_IO(void){
    int i,j;
	delete [] iTimeStep;  // = new int[nodenumber]();
	iTimeStep =NULL;
	delete [] iSysCalc;  //  = new int[nodenumber]();
	delete [] iSysCalc_repeated;  // = new int[nodenumber]();
	delete [] iSysPrevResults;  // = new int[nodenumber]();
	delete [] iSeqCalc;  // = new int[nodenumber]();
	delete [] nSeqCalc;  // = new int[nodenumber]();
	delete [] iNPCON_gas_exist;  // = new int[nodenumber]();

	for (i=0;i<nodenumber;i++){
		delete [] memEQ_mol_pcon_gas[i]; // = new double [NPCON[0]]();
		delete [] gas_pcon_exist[i]; // = new int [NPCON[0]]();
	}
	delete [] memEQ_mol_pcon_gas;
	delete [] gas_pcon_exist;

	delete [] iNPCON_aq_exist;  // = new int[nodenumber]();
//	memEQ_mol_pcon_aq = new double *[nodenumber];
//	aq_pcon_exist = new int *[nodenumber];
	for (i=0;i<nodenumber;i++){
		delete [] memEQ_mol_pcon_aq[i]; //= new double [NPCON[1]]();
		delete [] aq_pcon_exist[i]; // = new int [NPCON[1]]();
	}
	delete [] memEQ_mol_pcon_aq;
	delete [] aq_pcon_exist;

//	iNPCON_solids_exist = new int *[nodenumber];
//	solids_pcon_exist = new int **[nodenumber];
//	memEQ_mol_pcon_solids = new double **[nodenumber];
	for (i=0;i<nodenumber;i++){
		delete [] iNPCON_solids_exist[i]; // = new int [NPHAS-2]();   //be careful with the NPHAS
//		solids_pcon_exist[i] = new int *[NPHAS-2];
//		memEQ_mol_pcon_solids[i] = new double *[NPHAS-2];
		for (j=0;j<NPHAS-2;j++){
		  delete [] solids_pcon_exist[i][j]; // =new int [NPCON[j+2]]();
		  delete [] memEQ_mol_pcon_solids[i][j]; // =new double [NPCON[j+2]]();
		}
		delete [] solids_pcon_exist[i];
		delete [] memEQ_mol_pcon_solids[i];
	}
	delete [] solids_pcon_exist;
	delete [] memEQ_mol_pcon_solids;
	delete [] iNPCON_solids_exist;

	delete [] memEQ_Time;  // = new double [nodenumber]();
	delete [] memEQ_pH;  // = new double [nodenumber]();
	delete [] memEQ_molH2O;  // = new double [nodenumber]();
	delete [] memEQ_kgH2O;  // = new double [nodenumber]();

	for (i=0;i<nodenumber;i++){
		delete [] memEQ_mol_scom_gas[i];  // = new double [NSCOM]();
		delete [] memEQ_mol_scom_aq[i];  // = new double [NSCOM]();
		delete [] memEQ_mol_scom_solids[i];  // = new double [NSCOM]();
		delete [] memEQ_mol_scom_total[i];  // = new double [NSCOM]();
	}
	delete [] memEQ_mol_scom_gas;
	delete [] memEQ_mol_scom_aq;
	delete [] memEQ_mol_scom_solids;
	delete [] memEQ_mol_scom_total;

	delete [] pstat;  // = new int [NPHAS]();

	delete [] scom_elim_index;  // = new int [NSCOM]();
	delete [] mol_scom_total_save;  // = new double [NSCOM]();
	delete [] mol_scom_total;
	delete [] molH2O_in;
	delete [] sysT0;
    delete [] sysT;
	delete [] sysP;
	delete [] sysV;

}

/**************************************************************************
  GeoSys - Funktion: CEqlink::ExecuteEQLINK(void)
 
   Aufgabe:
     Call EQLINK	
  
  Ergebnis:
   
                                                                         
  Programmaenderungen:
     064/2006   MX  Erste Version
**************************************************************************/
void CEqlink::ExecuteEQLINK(void){
  long i;

  initEQLINK();
  
  if (aktueller_zeitschritt==0) AllocateMemory_EQCalc(); 

  cout<<endl << "      Execute ChemApp reactions  "<<endl;
  for (i=0; i<nodenumber;i++){
	action =1;
	while (action !=0){
//		if (action==1) callEQCHECK(i);
//		if (action==2) callEQCALC(i);  //action =0;
        callEQCALC(i);  
		if (action==3) callEQGETPRE(i);
	}
  }
  SetResultsBackMassTransport();
  

}
void CEqlink::callEQCHECK(long node){
	int i,j;
	double sumd_mol_scom;
	double dtemp=0.0;
	iPrevActions =0;

	//+++++++++++++++++++++++++++++
	//check input parameters
	if ((option[1] !=1) && option[1] !=3) {
	  DisplayErrorMsg("Error: option[1] must be 1 or 3!!");
	  DisplayErrorMsg("option[1] = 1 means gas entered");
	  DisplayErrorMsg("option[1] = 3 means gas eliminated!!");
	  exit(1);
	}
	if (sysP[node]<1.0e5){
	  DisplayErrorMsg("Error: Input system pressure should be no less than 1atm!!");
	  exit(1);
	}
	if (sysV[node] < 0.0){
	  DisplayErrorMsg("Error: Input system volume < 0.0E+00!!!");
	  exit(1);
	}
	if (option[3] ==1 && sysV[node] < 0.0){
	  DisplayErrorMsg("Error: Target volume calculation cannot be");
	  DisplayErrorMsg("enforced when input system volume is .le. 0!!!");
	  exit(1);
	}
    //+++++++++++++++++++++++++++++
	//enforced calculation
	if (option[2] ==1){
		action =2;
		if (verbosity >= 2){
			DisplayErrorMsg("Note(EQCHECK): A new calculation of thermodynamic equilibrium is enforced.");
		}
		return;
	}

	if (iSeqCalc[node] < nSeqCalc[node]){
		action =2;
		if (verbosity >= 2){
			DisplayErrorMsg("Note(EQCHECK): A new calculation of thermodynamic equilibrium is enforced.");
			DisplayErrorMsg("Less than nSeqCalc calculations have been done within the present sequence");
		}
		return;
	}

	if (sysV[node] > 0.0){
		action =2;
		if (verbosity >= 2){
			DisplayErrorMsg("Note(EQCHECK: A new calculation of thermodynamic equilibrium is enforced.");
			DisplayErrorMsg("Input system volume is > zero.");
		}
		return;
	}
	//+++++++++++++++++++++++++++++
    //  Check whether system components disappeared or appeared in which case
    //  calculation is enforced.
	if (KIND_OF_PFILE ==1) j=2;
	else j=1;

	sumd_mol_scom =0.0;
	i=1;
	j=j-1;

	while (i==1 && j<NSCOM){
		j +=j;
		if ((mol_scom_total[j] !=0.0 && memEQ_mol_scom_total[node][j] ==0.0) \
			||(mol_scom_total[j] ==0.0 && memEQ_mol_scom_total[node][j] !=0.0))
			i=0;
		else {
			sumd_mol_scom += abs(mol_scom_total[j]-memEQ_mol_scom_total[node][j]);
		}
	}
	if (i==0) {
		iSeqCalc[node] =0;
		//	errvec[1] =3;
		option[4] =0;
		
		if (verbosity >= 2){
			DisplayErrorMsg("Note(EQCHECK): A new calculation of thermodynamic equilibrium is enforced.");
			printf("  System component %d %s \n",j, SCNAME[j]);
			DisplayErrorMsg("appeared or disappeared!");
		}
		action =2;
		return;
	}
  //+++++++++++++++++++++++++++++
	if (sumd_mol_scom > lambda){
		if (verbosity >= 2){
			DisplayErrorMsg("Note(EQCHECK): A new calculation of thermodynamic equilibrium is enforced.");
			DisplayErrorMsg("Integral change of input system ");
			printf(" component inventory exceeds lambda: %d \n",sumd_mol_scom);
		}
		action =2;
		return;
	}
  //+++++++++++++++++++++++++++++
	if (memEQ_molH2O[node] > 0.0){
		dtemp = abs(molH2O_in[node]/memEQ_molH2O[node]-1.0);
		if (dtemp >= lambda){
		  if (verbosity >= 2){
			DisplayErrorMsg("Note(EQCHECK): A new calculation of thermodynamic equilibrium is enforced.");
			DisplayErrorMsg("Change of input system free water exceeds lambda.");
		  }
		  action =2;
		  return;
		}
	}

	  //+++++++++++++++++++++++++++++
	if (verbosity >= 2){
		DisplayErrorMsg("Note(EQCHECK): A new calculation of thermodynamic equilibrium is enforced.");
		DisplayErrorMsg("Change of input system free water exceeds lambda.");
	}
	action =3;
	if (iPrevActions >4	) {
		DisplayErrorMsg("Fatal error(EQCHECK): STOPPED.");
		exit(1);
	} 
	iPrevActions +=1;
	PrevActions[iPrevActions]=1;
	
}


void CEqlink::callEQCALC(long node){

  if (node ==40 && aktueller_zeitschritt==4 ){
  node =node;
  cout<<"    Timestep: "<< aktueller_zeitschritt<<endl;
}

    int i,j,k,m,n;
    long numerr=0,numcon=0, numvers=0;
    double dtemp=0.0,*vals;  //,act_H2O, act_H, act_OH, kw;
    double esum =0.0, sumsolutescharge=0.0;
    bool OnceMore=true;
    string str;
    vals = new double[3]();

//    cout<<"   callEQCALC--Start  "<<node<<endl;
//    if (aktueller_zeitschritt==0) AllocateMemory_EQCalc(); //Moved to CallEQLINK, only once
    #ifdef _DEBUG_CHEMAPP 
		DisplayErrorMsg("******Report from callEQCALC--Start******");
        tqvers(&numvers, &numerr); 
        cout<<"ChemApp version is: "<<numvers<<endl;
        cout<<"Node : "<<node<<endl;
    #endif
//        printf("\nCurrent time and node was: %d %ld. \n",aktuelle_zeit, node);

    if (aktueller_zeitschritt>=1){
        for (i=0;i<NSCOM;i++){
            mol_scom_gas_out[node][i]=0.0;
            mol_scom_aq_out[node][i]=0.0;
            mol_scom_solids_out[node][i]=0.0;
        }
        for (i=0;i<NPCON[0];i++){
            mol_pcon_gas_out[node][i]=0.0;
        }
        for (i=0;i<NPCON[1];i++){
            mol_pcon_aq_out[node][i]=0.0;
        }
        for (i=2;i<NPHAS;i++){
          for (j=0;i<NPCON[i];i++){
            mol_pcon_solids_out[node][i][j]=0.0;
          }
        }
        molH2O_out[node] = 0.0;  
        kgH2O_out[node] = 0.0;
        pH_out[node] = 0.0;
        Eh_out[node] = 0.0;

	    sysP_out[node] = 0.0;
	    sysV_out[node] = 0.0;
    } //if (aktueller_zeitschritt)

   //Eliminate system components (ToDo)

  //Remove all earlier conditions but keep system units as they are
    tqremc(-2, &numerr);
  //Set condition for incoming mass of water
    tqsetc("ia",2,1,molH2O_in[node],&numcon, &numerr);
  //MX  tqshow(&numerr);

  //set incoming amount of scom in mol, except EA 
    if (KIND_OF_PFILE ==1) j=1;
    else j=0;
    for (i=j;i<NSCOM;i++){
//MX      if (memEQ_mol_scom_total[node][i] !=0.0)
    if (memEQ_mol_scom_total[node][i] <0.0)
		cout<<"!!!Warning: Negative total system component of  "<<SCNAME[i]<<i+1<<" |"<<scientific<<memEQ_mol_scom_total[node][i]<<endl;
	if (memEQ_mol_scom_total[node][i] > 0.0)
        tqsetc("ia",0,i+1,memEQ_mol_scom_total[node][i],&numcon, &numerr);
    }
    #ifdef _DEBUG_CHEMAPP
      tqshow(&numerr); //MX
      cout << " " <<endl;
      for (i=0;i<NSCOM;i++){
        tqgetr("A ",0,i+1,&dtemp, &numerr);
        if (dtemp != 0.0)
       //   cout<<"400 |  "<<i+1<<"|IA"<<"|"<<scientific<<dtemp<<"          <<"<<SCNAME[i]<<endl;
          cout<<"400 |  "<<i+1<<" |"<<scientific<<dtemp<<"           <<"<<SCNAME[i]<<endl;
      }

    #endif 

  //Set status of gas phase
    m=0;
    if (option[1] != pstat[1]){ 
        tqcsp(1,PSTATNAME[option[1]],&numerr);
        pstat[1] = option[1];
        m=1;
    }

    //Set T,P
    tqsetc("T",0,0,sysT[node],&numcon, &numerr);
    tqsetc("P",0,0,sysP[node],&numcon, &numerr);

    //Check whether too many calculations had to be repeated
    if (option[3] != 2)
      if (iSysCalc_repeated[node]> 100  && (iSysCalc_repeated[node]/iSysCalc[node] > 0.1)) {
        option[3] = 2;  //volume as target calculation, must with gas!
      }
    if (option[3] == 2){
        tqsetc("V",0,0,sysV[node],&numcon, &numerr);
        dtemp = sysV[node];
    }
    else dtemp = 0.0;

    //++++++++++Calculate the equilibrium++++++++++
    while (OnceMore){
        if (verbosity > 3) 
          tqshow(&numerr);
        if ((iSysCalc[node] > 0) && m==0 && option[4]==1){
          if (dtemp == 0.0) 
            tqcen(" ",0,0, vals,&numerr);
          else {
            vals[1]=sysP[node];
            tqcen("P",0,0, vals,&numerr);
          }
        }
        else {
          if (dtemp == 0.0){
            #ifdef _DEBUG_CHEMAPP
              tqcel(" ",0,0, vals,&numerr);
            #else
              tqce(" ",0,0, vals,&numerr);
            #endif
          }
          else {
            vals[1]=sysP[node];
            tqcen("P",0,0, vals,&numerr);
          }
        }

        //Check whether calculated system volume is sufficiently
        // close to SysVol_in! If not, new calculation is needed
        if ((numerr == 0) && (sysV[node] > 0.0)){
            tqgetr("VT",0,0,&sysV_out[node],&numerr);
            tqgetr("P ",0,0,&sysP_out[node],&numerr);
            dtemp = abs(sysV_out[node]/sysV[node]-1.0);
            if (verbosity>3) {
              printf("\nNode    : %d \n",node);
              printf("sysV_in : %d \n",sysV[node]);
              printf("sysV_out: %d \n",sysP_out[node]);
              printf("Change rate : %d \n",dtemp);
            }
            if (dtemp > (0.01*lambda)) {  //! Repeat calculation with volume as target variable
              tqsetc("V",0,0,sysV[node],&numcon, &numerr);
              dtemp = sysV[node];
              OnceMore = true;
              if (verbosity>0) {
		        DisplayErrorMsg("******EQCALC (Warning):******");
		        DisplayErrorMsg("Unrestricted calculated system volume");
		        DisplayErrorMsg("deviates too much from input system volume.");
                printf("\nThreshold value was: %d \n",0.01*lambda);
		        DisplayErrorMsg("Calculation will be repeated with");
		        DisplayErrorMsg("input system volume as target variable.");
                iSysCalc_repeated[node] = iSysCalc_repeated[node] + 1;
              }
            }
            else 
              OnceMore = false;
        }
        else 
          OnceMore = false;

    }//while

    //Get CHEMAPP results
    if (numerr !=0 || ( aktuelle_zeit == 400 && node == 144) ) {
		        DisplayErrorMsg("******EQCALC fails ******");
//		        DisplayErrorMsg("Unrestricted calculated system volume");
//		        DisplayErrorMsg("deviates too much from input system volume.");
				cout << " ChemApp error No. is " << numerr <<endl;
                cout << " Aktuelle_zeit is " << aktuelle_zeit <<endl;
				cout << " Node is " << node <<endl;
				
				tqshow(&numerr); //MX
	if (node == 145) 	exit(1);  //MX for debug only!!

//		        DisplayErrorMsg("Calculation will be repeated with");
//		        DisplayErrorMsg("input system volume as target variable.");
}
    if (!numerr){
        iSysCalc[node] = iSysCalc[node]+1;
        iSeqCalc[node] = iSeqCalc[node]+1;

        tqgetr("P ",0,0,&sysP_out[node],&numerr);
        tqgetr("VT ",0,0,&sysV_out[node],&numerr);
        tqgetr("T ",0,0,&sysT0[node],&numerr);

        if (KIND_OF_PFILE == 1) k=1;
        else k=0;
        for (i=0;i<NPHAS;i++){
            tqgetr("V ",i+1,0,&PVOL_out[node][i],&numerr);
            tqcsu("Amount ", "gram ", &numerr);
            tqgetr("A  ",i+1,0,&PMASS_out[node][i],&numerr);
            tqcsu("Amount ", "mol ", &numerr);
            if (PMASS_out[node][i]>0.0 ||i>=2) {  //solid phase value should be set with current data!!
                for (j=0;j<NPCON[i];j++){
                    tqgetr("A ",i+1,j+1,&dtemp, &numerr);
                    if (dtemp !=0.0 ||i>=2){
                      str=PNAME[i];
                      if (str.find("GAS") !=string::npos){
                        mol_pcon_gas_out[node][j]=dtemp;
                        for (n=k;n<NSCOM;n++)
                          mol_scom_gas_out[node][n] += dtemp*PCSTOI_GAS[j][n];
                      }
                      else if (str.find("AQUEOUS") !=string::npos){
                        mol_pcon_aq_out[node][j]=dtemp;

                 /*  Calculate charge balance and sum of charges in
                   aqueous solution after calculation. It is meaningful
                   only with scoms = elements. The new version of chemapp
                   allegedly is able to return charge balance also for
                    parameterfiles with scom=species.*/
                        if (KIND_OF_PFILE == 1){
                          esum += dtemp*PCSTOI_AQ[j][0];  //charge balance
                          sumsolutescharge += dtemp*abs(PCSTOI_AQ[j][0]);
                        }
                        if (j==0){
                          kgH2O_out[node] = dtemp*PCMMASS_AQ[j]/1000.;
                          molH2O_out[node] = dtemp;
                        }
                        for (n=k;n<NSCOM;n++)
                          mol_scom_aq_out[node][n] += dtemp*PCSTOI_AQ[j][n];
                      } //else if
                      else {  //solids
                        mol_pcon_solids_out[node][i-2][j] = dtemp;
						if (dtemp > 0.0)
                          for (n=k;n<NSCOM;n++)
                            mol_scom_solids_out[node][n] += dtemp*PCSTOI_SOLIDS[i][j][n];

                      }
                    } //if (dtemp!=0)
                } //for j
            }  //if 
            else {
                if (NPHAS==2) molH2O_out[node] = 0.0;  //We are in the aqueous phase!??
            }
            
        } //for NPHAS

        if (PMASS_out[node][1]>0.0) {
            tqgetr("AC",2,2,&dtemp, &numerr);
            pH_out[node] = - log10(dtemp);
            if (pstat[1]==2){
              dtemp =0.0;
              // tqgetr("Eh",2,0,&dtemp, &numerr);   //for Eh output // 12.02.2009 HS, MX, disabled this. 
              Eh_out[node] = dtemp;
            }
        }
        else pH_out[node] = - 999.00;  //no water
      action=0;
    } //if (!numerr) CHEMAPP calculation OK
    delete [] vals;
    vals =NULL;
}

/**************************************************************************
  GeoSys - Funktion: CEqlink::AllocateMemory_IO(void)
 
   Aufgabe:
     Call EQLINK	
  
  Ergebnis:
   
                                                                         
  Programmaenderungen:
     064/2006   MX  Erste Version
**************************************************************************/
void CEqlink::AllocateMemory_EQCalc(void){

    int i,j;

    molH2O_out = new double [nodenumber]();
    kgH2O_out = new double [nodenumber]();
    pH_out = new double [nodenumber]();
    Eh_out = new double [nodenumber]();

	sysP_out = new double [nodenumber]();
	sysV_out = new double [nodenumber]();

	PVOL_out = new double *[nodenumber];
	PMASS_out = new double *[nodenumber];
	for (i=0;i<nodenumber;i++){
		PVOL_out[i] = new double [NPHAS]();
		PMASS_out[i] = new double [NPHAS]();
	}

	mol_scom_gas_out = new double *[nodenumber];
	mol_scom_aq_out = new double *[nodenumber];
	mol_scom_solids_out = new double *[nodenumber];
	for (i=0;i<nodenumber;i++){
		mol_scom_gas_out[i] = new double [NSCOM]();
		mol_scom_aq_out[i] = new double [NSCOM]();
		mol_scom_solids_out[i] = new double [NSCOM]();
	}

	mol_pcon_gas_out = new double *[nodenumber];
	for (i=0;i<nodenumber;i++){
		mol_pcon_gas_out[i] = new double [NPCON[0]]();
	}

	mol_pcon_aq_out = new double *[nodenumber];
	for (i=0;i<nodenumber;i++){
		mol_pcon_aq_out[i] = new double [NPCON[1]]();
	}

	mol_pcon_solids_out = new double **[nodenumber];
	for (i=0;i<nodenumber;i++){
		mol_pcon_solids_out[i] = new double *[NPHAS-2];
		for (j=0;j<NPHAS-2;j++){
		  mol_pcon_solids_out[i][j] =new double [NPCON[j+2]]();
		}
	}

}

void CEqlink::DeleteMemory_EQCalc(void){

    int i,j;

    delete [] molH2O_out;
    delete [] kgH2O_out;
    delete [] pH_out;
    delete [] Eh_out;

	delete [] sysP_out;
	delete [] sysV_out;

	for (i=0;i<nodenumber;i++){
		delete [] PVOL_out[i];  
		delete [] PMASS_out[i];  
	}
	delete [] PVOL_out;
	delete [] PMASS_out;

	for (i=0;i<nodenumber;i++){
		delete [] mol_scom_gas_out[i];  
		delete [] mol_scom_aq_out[i];  
		delete [] mol_scom_solids_out[i]; 
	}
	delete [] mol_scom_gas_out;
	delete [] mol_scom_aq_out;
	delete [] mol_scom_solids_out;

	for (i=0;i<nodenumber;i++){
		delete [] mol_pcon_gas_out[i]; // = new double [NPCON[0]]();
	}
	delete [] mol_pcon_gas_out;

	for (i=0;i<nodenumber;i++){
		delete [] mol_pcon_aq_out[i]; //= new double [NPCON[1]]();
	}
	delete [] mol_pcon_aq_out;

	for (i=0;i<nodenumber;i++){
		for (j=0;j<NPHAS-2;j++){
		  delete [] mol_pcon_solids_out[i][j]; // =new double [NPCON[j+2]]();
		}
		delete [] mol_pcon_solids_out[i];
	}
	delete [] mol_pcon_solids_out;
}

void CEqlink::callEQGETPRE(long node){
	action =0;
    DisplayErrorMsg("******callEQGETPRE (Warning):******");
    DisplayErrorMsg("******Not implemented function!!!******");

}

void CEqlink::SetResultsBackMassTransport(void){
    int i,j,k,n;
    CRFProcess *m_pcs = NULL;
    int idx0, idx1;
    double temperature=0.0;
    int timelevel = 1;
    int comp =0;
    double temp=0.0;
    string str,str1;

    if (KIND_OF_PFILE == 1) {

	  for(j=0;j<(int)pcs_vector.size();j++){ // for all processes
		m_pcs = pcs_vector[j];
		if(m_pcs->pcs_type_name.compare("MASS_TRANSPORT")==0){ // if it is a mass transport process
			comp = m_pcs->pcs_component_number; // get component number
			str=m_pcs->pcs_primary_function_name[0];
            idx0=m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]);
            idx1 = idx0+1;
			if (cp_vec[comp]->mobil){
				for (k=0; k<NSCOM;k++){   //pH  -> H+
//MX				    if (strcmp(m_pcs->pcs_primary_function_name[0],"pH") ==0 && strcmp(SCNAME[k],"H") ==0){				    
					if (strcmp(m_pcs->pcs_primary_function_name[0],"pH") ==0 && strcmp(SCNAME[k],"EA") ==0){
					  for(i=0;i<nodenumber;i++){
                        m_pcs->SetNodeValue(i, idx1,pH_out[i]);
					  } //for i
					  break;
					}  //if
					else if (strcmp(m_pcs->pcs_primary_function_name[0],SCNAME[k]) ==0) {
					  for(i=0;i<nodenumber;i++){
                        if (strcmp(SCNAME[k],"O")==0) 
                            mol_scom_aq_out[i][k] -= molH2O_out[i]; //H2O extract
					    if (strcmp(SCNAME[k],"H")==0) 
                            mol_scom_aq_out[i][k] -= 2*molH2O_out[i]; //H2O extract

                        m_pcs->SetNodeValue(i, idx1,mol_scom_aq_out[i][k]);
						temp =  \
							m_pcs->GetNodeValue(i,m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) \
							+timelevel);
					  }
					  break;
					}  //else if
				}  //for k
				for(i=0;i<nodenumber;i++){
				// get concentration values
				  if((m_pcs->GetNodeValue(i,idx1) < 0.0)  ) {  //MX, 01.2005 val_in[comp][i]&& (strcmp(name[comp],"pe")!= 0)
					if(abs(m_pcs->GetNodeValue(i,idx1)) > MKleinsteZahl) {
					  DisplayMsg(" Neg. conc for component "); DisplayLong((long) comp);
					  DisplayMsg(" at node "); DisplayLong((long)i);
					  DisplayMsg("; conc = "); DisplayDouble(m_pcs->GetNodeValue(i,idx1),0,0);
					  DisplayMsgLn(" ");
				    }
				    memEQ_mol_scom_aq[i][j] = 0.0 ;
				  }
				}  //for i
			}  //if mobile
			else if (cp_vec[comp]->mobil==0){ //immobile
                for (k=0; k<NSCOM;k++){  
//				    if (strcmp(m_pcs->pcs_primary_function_name[0],"pH") ==0 && strcmp(SCNAME[k],"H") ==0){
				    if (strcmp(m_pcs->pcs_primary_function_name[0],"pH") ==0 ){
					  for(i=0;i<nodenumber;i++){
                        m_pcs->SetNodeValue(i, idx1,pH_out[i]);
					  } //for i
					  continue;
					}  //if
				    if (strcmp(m_pcs->pcs_primary_function_name[0],"Eh") ==0){
					  for(i=0;i<nodenumber;i++){
                        m_pcs->SetNodeValue(i, idx1,Eh_out[i]);
					  } //for i
					  break;
					}  //if
                } //for k
				for (k=2; k<NPHAS;k++){  //0 -gas, 1-aq
					for (n =0;n<NPCON[k]; n++){
					  str1 = PCNAME_SOLIDS[k][n];
					  if (str1.find(str)!=string::npos){
						for(i=0;i<nodenumber;i++){
                            m_pcs->SetNodeValue(i, idx1,mol_pcon_solids_out[i][k-2][n]);
						    temp =  \
							  m_pcs->GetNodeValue(i,m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) \
							  +1);
						}
					    break;
					  } //if
					}  //for n
				}	//for k		
			}//if immobile
		}//if pcs = mass transport
	  }	//processes
    } //if (KIND_OF_PFILE == 1)
    else if (KIND_OF_PFILE == 2) {
	  for(j=0;j<(int)pcs_vector.size();j++){ // for all processes
		m_pcs = pcs_vector[j];
		if(m_pcs->pcs_type_name.compare("MASS_TRANSPORT")==0){ // if it is a mass transport process
			comp = m_pcs->pcs_component_number; // get component number
			str=m_pcs->pcs_primary_function_name[0];
            idx0=m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]);
            idx1 = idx0+1;
			if (cp_vec[comp]->mobil){
				for (k=0; k<NPCON[1];k++){   //pH  -> H+
				    if (strcmp(m_pcs->pcs_primary_function_name[0],"pH") ==0 && strcmp(PCNAME_AQ[k],"H<+>") ==0){
					  for(i=0;i<nodenumber;i++){
                        m_pcs->SetNodeValue(i, idx1,pH_out[i]);
					  } //for i
					  break;
					}  //if
					else if (strcmp(m_pcs->pcs_primary_function_name[0],PCNAME_AQ[k]) ==0) {
					  for(i=0;i<nodenumber;i++){
                        m_pcs->SetNodeValue(i, idx1, mol_pcon_aq_out[i][k]); //  mol_scom_aq_out[i][k]);
						temp =  \
							m_pcs->GetNodeValue(i,m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) \
							+timelevel);
					  }
					  break;
					}  //else if
				}  //for k
				for(i=0;i<nodenumber;i++){
				// get concentration values
				  if((m_pcs->GetNodeValue(i,idx1) < 0.0)  ) {  //MX, 01.2005 val_in[comp][i]&& (strcmp(name[comp],"pe")!= 0)
					if(abs(m_pcs->GetNodeValue(i,idx1)) > MKleinsteZahl) {
					  DisplayMsg(" Neg. conc for component "); DisplayLong((long) comp);
					  DisplayMsg(" at node "); DisplayLong((long)i);
					  DisplayMsg("; conc = "); DisplayDouble(m_pcs->GetNodeValue(i,idx1),0,0);
					  DisplayMsgLn(" ");
				    }
				    memEQ_mol_scom_aq[i][j] = 0.0 ;
				  }
				}  //for i
			}  //if mobile
			else if (cp_vec[comp]->mobil==0){ //immobile
                for (k=0; k<NSCOM;k++){  
				    if (strcmp(m_pcs->pcs_primary_function_name[0],"pH") ==0 && strcmp(SCNAME[k],"H") ==0){
					  for(i=0;i<nodenumber;i++){
                        m_pcs->SetNodeValue(i, idx1,pH_out[i]);
					  } //for i
					  continue;
					}  //if
				    if (strcmp(m_pcs->pcs_primary_function_name[0],"Eh") ==0){
					  for(i=0;i<nodenumber;i++){
                        m_pcs->SetNodeValue(i, idx1,Eh_out[i]);
					  } //for i
					  break;
					}  //if
                } //for k
				for (k=2; k<NPHAS;k++){  //0 -gas, 1-aq
					for (n =0;n<NPCON[k]; n++){
					  str1 = PCNAME_SOLIDS[k][n];
					  if (str1.find(str)!=string::npos){
						for(i=0;i<nodenumber;i++){
                            m_pcs->SetNodeValue(i, idx1,mol_pcon_solids_out[i][k-2][n]);
						    temp =  \
							  m_pcs->GetNodeValue(i,m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) \
							  +1);
						}
					    break;
					  } //if
					}  //for n
				}	//for k		
			}//if immobile
		}//if pcs = mass transport
	  }	//processes

    }

/*for(comp=0; comp<this->number_of_comp;comp++){
	name = this->name[comp];
    m_pcs = PCSGet("MASS_TRANSPORT",name); //???
	idx = m_pcs->GetNodeValueIndex(name)+1;
	for(i=0;i<this->nodenumber;i++)
			m_pcs->SetNodeValue(i, idx,this->val_out[comp][i]);
}
*/


}

#endif