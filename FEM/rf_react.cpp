/*
rf_react.cpp
Reaction package to go with MTM2
*/

#include <cfloat>
#include <signal.h>
#include "display.h"
#include "makros.h"
#include "rf_pcs.h"
#include "rf_pcs.h"
#include "mathlib.h"
#include "rfmat_cp.h"
#include "rf_react.h"
#include "rf_ic_new.h"
#include "stdio.h"
#include "files0.h"
#include "rf_pcs.h"                               //OK_MOD"
#include <iostream>
#include <fstream>
#include "rf_tim_new.h"
#include "rf_mmp_new.h"
#include "rf_kinreact.h"
// Elem object
#include "fem_ele_std.h"
#ifdef CHEMAPP
#include "eqlink.h"
#endif

#include <vector>
using namespace std;

void destroy_react(void *item);

#define PHREEQC
int NR_INORG_COMP_EQU;
int NR_MIN_EQU;
extern char *file_name;
extern char *crdat;                               /* MX */
//REACTION_MODEL *rcml=NULL;
extern double gravity_constant;

// MDL: new coupling
#ifdef LIBPHREEQC
extern "C"
{
#include "phreeqc.h"
}


string libphreeqc_print;
#endif

vector <REACT*> REACT_vec;

/**************************************************************************
   ROCKFLOW - Funktion: ExecuteReactionsPHREEQC

   Aufgabe:
   Berechnet chemische Reaktionen zwischen den einzelnen Komponenten
   allererste VErsion

   Programmaenderungen:
   06/2003     SB         Erste Version
   06/2003     MX         Read and Reaction functions
   09/2003     SB         Included temperature
11/2003     SB         included time step fro PHEEQC kinetics
bugfix for many species (large files)
adapted to faster output file setup

**************************************************************************/
void REACT::ExecuteReactionsPHREEQC(void)
{

   long i,  ok = 0;
   // CRFProcess *m_pcs = NULL;
   // static REACTION_MODEL *rcml;
   FILE *indatei, *fphinp, *fsel_out=NULL;
   char fsout[80];

   DisplayMsgLn("ExecuteReactionsPHREEQC:");

   /* Initialize arrays of concentrations and array for reaction rates (list, in pretimeloop)*/
   GetTransportResults();

   /* Perform reaction step */
   /* --------------------------------------------------------------------------*/
   if(flag_pqc)
   {
      indatei=fopen(crdat,"r");
      /* do only at first timestep */
      if (indatei != NULL)
      {
         //		if (aktueller_zeitschritt<2){
         //		rcml = CreateReactionModel();
         ok= this->ReadReactionModel(indatei);
         //		}
      }
      else
      {
         printf("ERROR: Can not open file *.pqc");
         exit(1);
      }

      /* Read the input file (*.pqc) and set up the input file for PHREEQC ("phinp.dat")*/
      /* set input file name */
      strcpy(fsout, this->outfile);

      fphinp = fopen("phinp.dat", "w");
      if ((fphinp) && ok)
      {
         for(i=0;i<this->nodenumber;i++)
         {
            if(this->rateflag[i] > 0)
            {
               rewind(indatei);
               ok = ReadInputPhreeqc(i,indatei, fphinp);
            }
         }
         fclose(indatei);
         fclose(fphinp);
      }
      else
      {
         DisplayMsgLn("The file phinput.dat could not be opened !");
         exit(1);
      }

      /* Extern Program call to PHREEQC */
      if(ok) ok = Call_Phreeqc();
      if(ok ==0) exit(1);

      /* Set up the output values for rockflow after Phreeqc reaction*/
      fsel_out=fopen(fsout,"r");
      if ((ok) && !fsel_out)
      {
         DisplayMsgLn("The selected output file doesn't exist!!!");
         exit(1);
      }
      else if(ok)
      {
         ok = ReadOutputPhreeqc(fsout);
         if(!ok) DisplayMsgLn(" Error in call to PHREEQC !!!");
         fclose(fsel_out);
      }
   }                                              /* if flag */

   /* Calculate Rates */
   CalculateReactionRates();

   SetConcentrationResults();

   /* determine where to calculate the chemistry */
   // CalculateReactionRateFlag();

   /* pH and pe constant or variable */
   // ResetpHpe(rc, rcml);

   /* test output*/
   /*
    for(comp=0; comp<this->number_of_comp;comp++){
       DisplayMsg("component : "); DisplayLong(comp);DisplayMsg(", name = "); DisplayMsg(this->name[comp]);DisplayMsgLn(". ");
       DisplayMsgLn("val_in: ");
       for(i=0;i<this->nodenumber;i++){	DisplayDouble(this->val_in[comp][i],0,0); DisplayMsg(", ");}
       DisplayMsgLn(" ");
       DisplayMsgLn("val_out: ");
       for(i=0;i<this->nodenumber;i++){	DisplayDouble(this->val_out[comp][i],0,0); DisplayMsg(", ");}
       DisplayMsgLn(" ");

       DisplayMsgLn("rate : ");
   for(i=0;i<this->nodenumber;i++){	DisplayDouble(this->rate[comp][i],0,0); DisplayMsg(", ");}
   DisplayMsgLn(" ");
   }
   DisplayMsgLn("rateflag : ");
   for(i=0;i<this->nodenumber;i++){	DisplayDouble((double) this->rateflag[i],0,0); DisplayMsg(", ");}
   DisplayMsgLn(" ");
   */

}                                                 /* End of ExecuteReactionsPHREEQC */


/**************************************************************************
   ROCKFLOW - Funktion: ExecuteReactionsPHREEQCNew

   Aufgabe:
   Berechnet chemische Reaktionen zwischen den einzelnen Komponenten
   allererste Version

   Programmaenderungen:
   06/2003     SB         Erste Version
   06/2003     MX         Read and Reaction functions
   09/2003     SB         Included temperature
11/2003     SB         included time step fro PHEEQC kinetics
bugfix for many species (large files)
adapted to faster output file setup
01/2006     SB         ReImplementation as Class, IO streaming, bugfixes

**************************************************************************/
void REACT::ExecuteReactionsPHREEQCNew(void)
{

   long i, ii,  ok = 0;

   std::cout << "   ExecuteReactionsPHREEQCNew:" << std::endl;

   /* File handling - GeoSys input file */
   std::ifstream pqc_file (this->file_name_pqc.data(),ios::in);
   if (!pqc_file.good())
   {
      std::cout << "! Error in ExecuteReactionsPHREEQCNew: no Input File (*.pqc) found !" << std::endl;
      //         exit(1);
   }
   //	File handling - data exchange file to phreeqc, input to PHREEQC
   std::ofstream outfile (this->outfile_name.data(),ios::out);
   if(!outfile.is_open())
   {
      std::cout << "Error: Outfile phinp.dat could not be opened for writing " << std::endl;
      //        exit(1);
   }

   // Set up reaction model
   if((int)this->pqc_names.size()==0)
   {
      ok = this->ReadReactionModelNew(&pqc_file);
      if(!ok)
         std::cout << "Error setting up reaction model" << std::endl;
   }

   // Check for nodes without reactions
   if((int)this->check_no_reaction_nodes == false)
   {
      ok = this->CheckNoReactionNodes();
      if(!ok)
         std::cout << "Error when checking for nodes without reactions" << std::endl;
   }
   /* Read the input file (*.pqc) and set up the input file for PHREEQC ("phinp.dat")*/
   // Write input data block to PHREEQC for each node
   ii = 0;
   for(i=0;i<this->nodenumber;i++)
      if(this->rateflag[i] > 0)
   {
      pqc_file.seekg(0L,ios_base::beg);
      ok = WriteInputPhreeqc(i, /*&pqc_file,*/ &outfile);
      ii++;
   }

   //  Close *.pqc input file
   pqc_file.close();
   //  Close phinp.dat file with input for phreeqc.exe
   outfile.close();

   /* Extern Program call to PHREEQC */
   if(ok) ok = Call_Phreeqc();
   if(ok ==0)
   {
      std::cout << " Error executing PHREEQC.exe - Stopping "<< std::endl;
      std::cout.flush();
      //        exit(1);
   }

   if(ok)
   {
      ok = ReadOutputPhreeqcNew();
      if(!ok)
         std::cout << " Error in call to PHREEQC !!!" << std::endl;
   }

   std::cout << " Calculated equilibrium geochemistry at " << ii << " nodes." << std::endl;

   /* Calculate Rates */
   // CalculateReactionRates();

   /* determine where to calculate the chemistry */
   // CalculateReactionRateFlag();

   /* pH and pe constant or variable */
   // ResetpHpe(rc, rcml);

}                                                 /* End of ExecuteReactionsPHREEQCNew */


/**************************************************************************/
/* Constructor */
REACT::REACT(void)
{
   name = NULL;
   val_in=NULL;
   val_out=NULL;
   rate=NULL;
   rateflag = NULL;
   heatflag = 0;
   temperature = -1.0;
   nodenumber = 0;
   flag_pqc = false;
   this->check_no_reaction_nodes = false;

   rcml_number_of_master_species=0;
   rcml_number_of_equi_phases=0;
   rcml_number_of_ion_exchanges=0;
   rcml_number_of_gas_species=0;
   rcml_number_of_kinetics=0;
   rcml_pH_flag=1;
   rcml_pe_flag=1;
   rcml_heat_flag = 0;
   rcml_pH_charge = 0;
   rcml_number_of_pqcsteps = 1;                   //standard = 1;
   outfile = NULL;
   outfile_name = "phinp.dat";
   results_file_name = "phout_sel.dat";
   gamma_Hplus=-1.0;
}


/* Destructor */
REACT::~REACT(void)
{
}


/**************************************************************************
   ROCKFLOW - Funktion: CreateREACT

   Aufgabe:
   Stellt Datenstruktur für Ratenkommunikation zwischen Raecations and MTM2
   zur verfügung

   Programmaenderungen:
   06/2003     SB         Erste Version
   10/2010 TF restructured method a little bit, changed access to process type
**************************************************************************/
void REACT::CreateREACT(void)
{
   int np = 0;

   //rc = new REACT;
   size_t vector_size (pcs_vector.size());
   for (size_t i = 0; i < vector_size; i++)
   {
      //		CRFProcess *m_pcs = pcs_vector[i];
      ProcessType pcs_type (pcs_vector[i]->getProcessType());
      //		if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0)
      if (pcs_type == MASS_TRANSPORT)
         np++;
      //		if (m_pcs->pcs_type_name.compare("HEAT_TRANSPORT") == 0)
      if (pcs_type == HEAT_TRANSPORT)
         heatflag = 1;
      nodenumber = pcs_vector[i]->m_msh->GetNodesNumber(false);
      elenumber = (long) pcs_vector[i]->m_msh->ele_vector.size();

#ifdef REACTION_ELEMENT
                                                  //MX 06/2006
      if (nodenumber<elenumber) nodenumber = elenumber;
#endif
   }
   number_of_comp = np;

   np = np + heatflag;                            // allocate one for temperature
   // allocate memory for the concentration arrays
   name = (const char **) Malloc(sizeof(char *) * np);
   val_in = (double **) Malloc(sizeof(double *) * np);
   val_out = (double **) Malloc(sizeof(double *) * np);
   rate = (double **) Malloc(sizeof(double *) * np);
   rateflag = (int *) Malloc(sizeof(int) * nodenumber);

   for (int comp = 0; comp < np; comp++)
   {
      val_in[comp] = (double *) Malloc(sizeof(double) * nodenumber);
      val_out[comp] = (double *) Malloc(sizeof(double) * nodenumber);
      rate[comp] = (double *) Malloc(sizeof(double) * nodenumber);
      name[comp] = (char *) Malloc(sizeof(char) * 80);
   }

   // REACT_vec.push_back(rc);
}


/**************************************************************************
   ROCKFLOW - Funktion: DestroyREACT

   Aufgabe:
   Zerstört Datenstructur REACT

   Programmaenderungen:
   06/2003     SB         Erste Version
**************************************************************************/
void DestroyREACT(void)
{
   REACT* rc = NULL;
   rc=rc->GetREACT();
   /* Free memeory for concentration arrays */
   if(rc != NULL)
   {
      rc->val_in   = (double **)Free(rc->val_in);
      rc->val_out  = (double **)Free(rc->val_in);
      rc->rate     = (double **)Free(rc->val_in);
      rc->name     = (const char **)  Free(rc->name);
      rc->rateflag = (int *) Free(rc->rateflag);
      rc = (REACT *) Free(rc);
   }
}


/**************************************************************************
   ROCKFLOW - Funktion: REACT::InitREACT

   Aufgabe:
   Inserts reaction rate of component at index in field RECT->rate

   Programmaenderungen:
   06/2003     SB         Erste Version
   10/2010 TF restructured method a little bit, changed access to process type
**************************************************************************/
void REACT::InitREACT(void)
{
   int comp;
   long i;
   CRFProcess *m_pcs = NULL;
   int timelevel = 1;                             // concentrations are in new timelevel
   //	int phase = 0; //single phase so far // TF not used

   /* Initialize arrays of concentrations and array for reaction rates */
   size_t np (pcs_vector.size());
   for (size_t j = 0; j < np; j++)                // for all processes
   {
      m_pcs = pcs_vector[j];
      //		if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0) { // if it is a mass transport process
                                                  // if it is a mass transport process
      if (m_pcs->getProcessType () == MASS_TRANSPORT)
      {
         comp = m_pcs->pcs_component_number;      // get component number
                                                  // get component name
         this->name[comp] = (char *) cp_vec[comp]->compname.data();
         for (i = 0; i < this->nodenumber; i++)
         {
            // get concentration values
            this->val_in[comp][i] = m_pcs->GetNodeValue(i,
               m_pcs->GetNodeValueIndex(
               m_pcs->pcs_primary_function_name[0])
               + timelevel);
            if ((this->val_in[comp][i] < 0.0) && (strcmp(this->name[comp],
               "pe") != 0))                       //MX, 01.2005
            {
               if (abs(this->val_in[comp][i]) > MKleinsteZahl)
               {
                  DisplayMsg(" Neg. conc for component ");
                  DisplayLong((long) comp);
                  DisplayMsg(" at node ");
                  DisplayLong((long) i);
                  DisplayMsg("; conc = ");
                  DisplayDouble(this->val_in[comp][i], 0, 0);
                  DisplayMsgLn(" ");
               }
               this->val_in[comp][i] = 0.0 * this->val_in[comp][i];
            }
            this->val_out[comp][i] = this->val_in[comp][i];
            this->rate[comp][i] = 0.0;
         }
      }
   }

   for (i = 0; i < this->nodenumber; i++)
      this->rateflag[i] = 1;
   this->countsteps = 50;

   // Get Temperature values
   if (this->heatflag > 0)
   {
      this->name[np - 2] = "temp";                //MX CMCD
      m_pcs = PCSGet("HEAT_TRANSPORT");
      int index = m_pcs->GetNodeValueIndex("TEMPERATURE1");
      for (i = 0; i < this->nodenumber; i++)
                                                  //OK PCSGetNODTemperature1L(i)//MX CMCD -2, Liquid Flow, Heat Transport
         this->val_in[np - 2][i] = m_pcs->GetNodeValue(i, index);
   }
}


/**************************************************************************
   ROCKFLOW - Funktion: REACT::InitREACT0

   Aufgabe:
   Inserts reaction rate of component at index in field RECT->rate

   Programmaenderungen:
   06/2006     MX         Erste Version
   10/2010 TF restructured method a little bit, changed access to process type
**************************************************************************/
void REACT::InitREACT0()
{
   long i;
   CRFProcess *pcs = NULL;
   CInitialCondition* ic = NULL;

   int timelevel = 1;                             // concentrations are in new timelevel
   //	int phase = 0; //single phase so far

   /* Initialize arrays of concentrations and array for reaction rates */
   size_t np (pcs_vector.size());
   for (size_t j = 0; j < np; j++)                // for all processes
   {
      pcs = pcs_vector[j];
      //		if (pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0) { // if it is a mass transport process
      if (pcs->getProcessType () == MASS_TRANSPORT)
      {
         int comp = pcs->pcs_component_number;    // get component number
                                                  // get component name
         this->name[comp] = (char *) cp_vec[comp]->compname.data();

         if (comp >= 0 && (cp_vec[comp]->mobil > 0))
            pcs->InterpolateTempGP(pcs, this->name[comp]);

         if (comp >= 0 && (cp_vec[comp]->mobil == 0))
         {
            ic = ic_vector[j];
            ic->setProcess(pcs);
            if (ic->getProcessPrimaryVariable()
               == convertPrimaryVariable(this->name[comp]))
            {
               int nidxe = pcs->GetElementValueIndex(this->name[comp]);
               ic->SetEle(nidxe);                 //set IC to elements for reactions //MX 06/2006
               ic->SetEle(nidxe+1);
            }                                     //for
         }

         for (i = 0; i < this->elenumber; i++)
         {
            // get concentration values
            this->val_in[comp][i] = pcs->GetElementValue(i,
               pcs->GetElementValueIndex(
               pcs->pcs_primary_function_name[0])
               + timelevel);
            if ((this->val_in[comp][i] < 0.0) && (strcmp(this->name[comp],
               "pe") != 0))                       //MX, 01.2005
            {
               if (abs(this->val_in[comp][i]) > MKleinsteZahl)
               {
                  DisplayMsg(" Neg. conc for component ");
                  DisplayLong((long) comp);
                  DisplayMsg(" at node ");
                  DisplayLong((long) i);
                  DisplayMsg("; conc = ");
                  DisplayDouble(this->val_in[comp][i], 0, 0);
                  DisplayMsgLn(" ");
               }
               this->val_in[comp][i] = 0.0 * this->val_in[comp][i];
            }
            this->val_out[comp][i] = this->val_in[comp][i];
            this->rate[comp][i] = 0.0;
         }
      }
   }

   for (i = 0; i < this->elenumber; i++)
      this->rateflag[i] = 1;
   this->countsteps = 50;

   // Get Temperature values for elements
   if (this->heatflag > 0)
   {
      this->name[np - 2] = "temp";
      pcs = PCSGet("HEAT_TRANSPORT");
      int idxT = pcs->GetElementValueIndex("TEMPERATURE1") + 1;
      pcs->InterpolateTempGP(pcs, "TEMPERATURE1");
      for (i = 0; i < this->elenumber; i++)
         this->val_in[np - 2][i] = pcs->GetElementValue(i, idxT);
   }
}

                                                  //MX
void CRFProcess::InterpolateTempGP(CRFProcess *m_pcs, std::string name)
{
   MshElemType::type EleType;
   int j;
   long i, enode;
   long group;
   double T_ele;
   double GP[3];
   static double Node_T[8];
   int index1;                                    //idxp,idxcp,idxS;
   CMediumProperties* m_mmp = NULL;
   CElem* elem =NULL;
   index1 = m_pcs->GetElementValueIndex(name)+1;  //->fem->interpolate(

   for (i=0;i<(long) m_pcs->m_msh->ele_vector.size();i++)
   {
      elem = m_pcs->m_msh->ele_vector[i];
      m_pcs->GetAssembler();

      // Activated Element
      group = elem->GetPatchIndex();
      m_mmp = mmp_vector[group];
      m_mmp->m_pcs=m_pcs;                         //m_pcs_mmp
      EleType = elem->GetElementType();
      if(EleType==MshElemType::TRIANGLE)          // Traingle
      {
         GP[0] = GP[1] = 0.1/0.3;
         GP[2] = 0.0;
      }
      else if(EleType==MshElemType::TETRAHEDRON)
         GP[0] = GP[1] = GP[2] = 0.25;
      else
         GP[0] = GP[1] = GP[2] = 0.0;

      m_pcs->fem->ConfigElement(elem);
      m_pcs->fem->setUnitCoordinates(GP);
      m_pcs->fem->ComputeShapefct(1);             // Linear
      for(j=0; j<elem->GetVertexNumber(); j++)
      {
         enode = elem->GetNodeIndex(j);
                                                  //m_pcs_mmp
         Node_T[j] =  m_pcs->GetNodeValue(enode,m_pcs->GetNodeValueIndex(name)+1);
      }
      T_ele = fem->interpolate(Node_T);
      m_pcs->SetElementValue(i,index1,T_ele);
      m_pcs->SetElementValue(i,index1-1,T_ele);
   }
}


                                                  //MX
void CRFProcess::ExtropolateTempGP(CRFProcess *m_pcs, std::string name)
{
   MshElemType::type EleType;
   int j;
   long i, enode, nn;
   long group;
   double GP[3];
   static double Node_T[8];
   double T_sum=0.0;
   int index1, index_nod;                         //idxp,idxcp,idxS;
   CMediumProperties* m_mmp = NULL;
   CElem* elem =NULL;

   index1 = m_pcs->GetElementValueIndex(name)+1;  //->fem->interpolate(
   index_nod = m_pcs->GetNodeValueIndex(name)+1;

   for (i = 0; i < m_msh->GetNodesNumber(false); i++)
      SetNodeValue(i,index_nod, 0.0);

   for (i=0;i<(long) m_pcs->m_msh->ele_vector.size();i++)
   {
      elem = m_pcs->m_msh->ele_vector[i];
      m_pcs->GetAssembler();

      // Activated Element
      group = elem->GetPatchIndex();
      m_mmp = mmp_vector[group];
      m_mmp->m_pcs=m_pcs;                         //m_pcs_mmp
      EleType = elem->GetElementType();
      if(EleType==MshElemType::TRIANGLE)          // Traingle
      {
         GP[0] = GP[1] = 0.1/0.3;
         GP[2] = 0.0;
      }
      else if(EleType==MshElemType::TETRAHEDRON)
         GP[0] = GP[1] = GP[2] = 0.25;
      else
         GP[0] = GP[1] = GP[2] = 0.0;

      m_pcs->fem->ConfigElement(elem);
      m_pcs->fem->setUnitCoordinates(GP);
      m_pcs->fem->ComputeShapefct(1);             // Linear
      for(j=0; j<elem->GetVertexNumber(); j++)
      {
         enode = elem->GetNodeIndex(j);
         Node_T[j] =  m_pcs->GetElementValue(i,index1);
         T_sum = m_pcs->GetNodeValue(enode, index_nod);
         m_pcs->SetNodeValue(enode,index_nod, T_sum+Node_T[j]);
      }
   }                                              //for

   // Average
   for (i = 0; i <(long)m_msh->GetNodesNumber(false); i++)
   {
      T_sum = m_pcs->GetNodeValue(i, index_nod);
      nn = (int) m_msh->nod_vector[i]->connected_elements.size();
      if(nn==0) nn =1;
      T_sum /= (double)nn;
      m_pcs->SetNodeValue(i,index_nod, T_sum);
   }
}


/**************************************************************************
   ROCKFLOW - Funktion: CalculateReactionRates

   Aufgabe:
   Calculates the reaction rates by: rate =  (val_out-val_in)/dt

   Programmaenderungen:
   06/2003     SB         Erste Version
**************************************************************************/
void REACT::CalculateReactionRates(void)
{
   /* Calculate Rates */

   int comp, i;
   double teta, help=0.0, maxi, dc, delta_t;

   CTimeDiscretization *m_tim = NULL;
   if(time_vector.size()>0)
      m_tim = time_vector[0];
   else
      std::cout << "Error in MPCCalcCharacteristicNumbers: no time discretization data !" << std::endl;
   delta_t = m_tim->CalcTimeStep();
   //OK_TIM delta_t = GetDt(aktueller_zeitschritt-1l);
   teta = 0.0;                                    /* teta = 0; concentration before reactions */
   for(comp=0; comp<this->number_of_comp; comp++)
   {
      maxi = 0.0;
      for(i=0;i<this->nodenumber;i++)
      {
         maxi = max(maxi, this->val_out[comp][i]);
         this->rate[comp][i] = (this->val_out[comp][i] - this->val_in[comp][i])/delta_t;

         help = (teta*this->val_out[comp][i] + (1.0-teta)*this->val_in[comp][i]);
         if(help < MKleinsteZahl)                 // umgekehrt wichten
         {
            help = (1.0-teta)*this->val_out[comp][i] + teta*this->val_in[comp][i];
            if(help < MKleinsteZahl)
               this->rate[comp][i] = 0.0;
            else
               this->rate[comp][i] = this->rate[comp][i] /help;
         }
         else
            this->rate[comp][i] = this->rate[comp][i] /help;

      }

      /* Reaction - Number dC */
      for(i=0;i<this->nodenumber;i++)
      {
         if(fabs(maxi)< MKleinsteZahl)
            dc = 0;
         else
            dc = this->rate[comp][i]*delta_t/maxi;
         //SB:todo		MTM2SetElementDamkohlerNumNew(i,0,comp,dc);
      }
   }                                              /* end for(comp=... */
}


/**************************************************************************
   ROCKFLOW - Funktion: SetConcentrationResults

   Aufgabe:
   Save concentrations after reaction in concentration array

   Programmaenderungen:
   06/2003     SB         Erste Version
**************************************************************************/
void REACT::SetConcentrationResults(void)
{
   /*  */
   int comp, timelevel, np, idx;
   long i;
   string name;
   CRFProcess *m_pcs = NULL;

   //Fix Process here ??
   //Rücksrache mit Sebastian

   timelevel = 1;                                 // concentrations are in new timelevel
   np = (int)pcs_vector.size();
   /*  SB4218
    for(comp=0; comp<this->number_of_comp;comp++){
       for(i=0;i<this->nodenumber;i++){
         name = this->name[comp];
         idx = m_pcs->GetNodeValueIndex(name)+1;
         m_pcs->SetNodeValue(i, idx,this->val_out[comp][i]);
       }
    }
    */
   for(comp=0; comp<this->number_of_comp;comp++)
   {
      name = this->name[comp];
      m_pcs = PCSGet("MASS_TRANSPORT",name);      //???
      idx = m_pcs->GetNodeValueIndex(name)+1;
      for(i=0;i<this->nodenumber;i++)
         m_pcs->SetNodeValue(i, idx,this->val_out[comp][i]);
   }

}


/**************************************************************************
   ROCKFLOW - Funktion: SetConcentrationResultsEle

   Aufgabe:
   Save concentrations after reaction in concentration array

   Programmaenderungen:
   06/2006     MX         Erste Version
**************************************************************************/
void REACT::SetConcentrationResultsEle(void)
{
   int comp, timelevel, np, idx;
   long i;
   string name;
   CRFProcess *m_pcs = NULL;

   timelevel = 1;                                 // concentrations are in new timelevel
   np = (int)pcs_vector.size();

   for(comp=0; comp<this->number_of_comp;comp++)
   {
      name = this->name[comp];
      m_pcs = PCSGet("MASS_TRANSPORT",name);      //???
      idx = m_pcs->GetElementValueIndex(name)+1;

      for(i=0;i<this->elenumber;i++)
         m_pcs->SetElementValue(i, idx,this->val_out[comp][i]);

      //Extropolate element data (center) to nodes
      if (CPGetMobil(m_pcs->GetProcessComponentNumber())>= 0)
         m_pcs->ExtropolateTempGP(m_pcs, name);

   }
}


REACT* REACT::GetREACT(void)
{
   REACT *rc = NULL;
   /* Tests */
   if(REACT_vec.capacity() > 0)
      rc = REACT_vec[0];
   return rc;
}


/**************************************************************************
   ROCKFLOW - Funktion: ReadReactionModel

   Aufgabe:
   Liest die Eingabedatei *.pqc und geben die Model-values zurück

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R REACTION_MODEL *rcml:Zeiger des REACTION_MODELs
   E FILE *File :Eigabedatei

   Ergebnis:
0 bei Fehler oder Ende aufgrund Dateitest, sonst 1

Programmaenderungen:
06/2003     MX         Erste Version
10/2010 TF restructured method a little bit, changed access to process type
**************************************************************************/
int REACT::ReadReactionModel(FILE *File)
{
   int n_master_species = 0;
   int n_equi_phases = 0;
   int n_ion_exchange = 0;
   int beginn = 0, p = 0;
   int pos;
   int nj;
   char str[256], *sub, *sub1;

   FILE *indatei = NULL;
   indatei = File;

   CRFProcess* m_pcs = NULL;

   /* Open input file  and read the reaction model values*/
   if (indatei == NULL)                           /*input dateien does not exist*/
   {
      DisplayMsgLn("");
      DisplayMsgLn(" The input file *.pqc does not exist!!!");
      exit(1);
   }                                              /*end if*/

   /* zeilenweise lesen */
   while (fgets(str, 256, indatei))
   {
      /* Schleife ueber Keyword Solution */
      pos = 0;
      beginn = 1;
      p = 0;
      while ((!strstr(str, "END")) && (!strstr(str, "#ende"))
         && StringReadStr(&sub, str, &p))
      {
         if (strcmp(sub, "SOLUTION") == 0)
         {
            while (fgets(str, 256, indatei) && (!strstr(str, "#ende")))
            {
               if (strstr(str, "# comp"))
               {
                  StringReadStr(&sub, str, &p);
                  if (!strcmp(sub, "pH") == 0 && !strcmp(sub, "pe") == 0)
                  {
                     n_master_species += 1;
                  }
               }
               if (strstr(str, "# temp"))
               {
                  size_t no_processes (pcs_vector.size());
                  for (size_t i = 0; i < no_processes; i++)
                  {
                     m_pcs = pcs_vector[i];
                     //							if (m_pcs->pcs_type_name.find("HEAT_TRANSPORT") != string::npos) {
                     if (m_pcs->getProcessType () == HEAT_TRANSPORT)
                     {
                        //			  if(GetRFProcessProcessingAndActivation("HT"))
                        this->rcml_heat_flag = 1;
                        break;
                     }
                  }
               }
            }
            this->rcml_number_of_master_species = n_master_species;
         }
         else if (strcmp(sub, "EQUILIBRIUM_PHASES") == 0)
         {
            while (fgets(str, 256, indatei) && (!strstr(str, "#ende")))
            {
               if (strstr(str, "# comp"))
               {
                  StringReadStr(&sub, str, &p);
                  n_equi_phases += 1;
               }
            }
            this->rcml_number_of_equi_phases = n_equi_phases;
         }                                        /*end if*/
         else if (strcmp(sub, "EXCHANGE") == 0)
         {
            while (fgets(str, 256, indatei) && (!strstr(str, "#ende")))
            {
               if (strstr(str, "# comp"))
               {
                  StringReadStr(&sub, str, &p);
                  n_ion_exchange += 1;
               }
            }
            this->rcml_number_of_ion_exchanges = n_ion_exchange;
         }                                        /*end if*/
         //SB:neu
         else if (strcmp(sub, "SELECTED_OUTPUT") == 0)
         {
            while (fgets(str, 256, indatei) && (!strstr(str, "#ende")))
            {
               if (strstr(str, "-file"))
               {
                  DisplayMsgLn("-file in *.pqc found");
                  p = 0;
                  StringReadStr(&sub, str, &p);
                  StringReadStr(&sub1, &str[p], &p);
                  this->outfile = sub1;
               }
            }
         }

         else
         {
            while (fgets(str, 256, indatei) && (!strstr(str, "#ende")))
            {
            }
            break;
            //        return 1;
         }
      }                                           /*end while 2*/
   }                                              /*end while 1*/
   nj = rcml_number_of_master_species + rcml_number_of_equi_phases
      + rcml_number_of_ion_exchanges;
   if (nj + 2 != this->number_of_comp)
   {
      DisplayMsgLn(
         "!!!Error:Number of components in file *.pqc is not equal to that in file *.rfd!");
      //      fclose(indatei);
      //      return 0;
   }
   //    fclose(indatei);
   return 1;
}


/**************************************************************************
   ROCKFLOW - Funktion: ReadReactionModelNew

   Aufgabe:
   Liest die Eingabedatei *.pqc und gibt die Parameter der Schnittstelle zurück

    indices of vectors pqc_name, pqc_index and pqc_process
    n1 = number_of_master_species
    | master species 1   | <- 0
    |      ...           |
    |master species n1   | <- n1-1
|       pH           | <- n1
|       H+           | <- n1+1
|       pe           | <- n1+2
|equilibrium phase 1 | <- n1+3
|      ...           |
|equilibrium phase n2| <- n1+3+n2-1
|exchange species 1  | <- n1+3+n2
|      ...           |
|exchange species n3 | <- n1+3+n2+n3-1

Programmaenderungen:
01/2006     SB         Erste Version
10/2010 TF restructured method a little bit, changed access to process type
**************************************************************************/
int REACT::ReadReactionModelNew( ifstream *pqc_infile)
{

   int n_master_species=0;
   int n_equi_phases=0;
   int n_ion_exchange=0;
   int n_gas_species=0;
   int idx;
   int nj;
   char line[MAX_ZEILE];
   string line_string, dummy, speciesname;
   CRFProcess* m_pcs = NULL;
   std::stringstream in;
   int pH_found = 0;                              // Hplus_found = 0, count = -1;

   /* File handling */
   pqc_infile->seekg(0L,ios::beg);

   /* zeilenweise lesen */
   while(!pqc_infile->eof())
   {
      pqc_infile->getline(line,MAX_ZEILE);
      line_string = line;
      if(line_string.find("#STOP")!=string::npos)
         break;

      /* Schleife ueber Keyword Solution */
                                                  // keyword found
      if(line_string.find("SOLUTION")!=string::npos)
      {

         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile->getline(line,MAX_ZEILE);
            line_string = line;
            if (line_string.find("# comp")!=string::npos)
            {
               if (line_string.find("pH") == string::npos && line_string.find("pe") == string::npos)
               {
                  n_master_species +=1;           // this excludes now pH and pe
                  in.str(line_string);
                  in >> speciesname ;
                  m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
                  if(m_pcs == NULL)
                  {
                     cout << " PHREEQC SOLUTION SPECIES not in GeoSys components - Stopping" << endl;
                     cout.flush(); exit(1);
                  }
                                                  // old timelevel
                  idx = m_pcs->GetNodeValueIndex(speciesname)+1;
                  in.clear();
                                                  //store species name in vector names
                  pqc_names.push_back(speciesname);
                                                  // store process name in vector processes
                  pqc_process.push_back(m_pcs->pcs_number);
                  pqc_index.push_back(idx);
                  /*                    if(speciesname.compare("H+") == 0){
                                          Hplus_found = 1;
                                          cout << "H+ found in GeoSys species" << endl;
                                          this->pqc_Hplus_index = n_master_species-1; // Store index for H+ in vectors
                                      } */
               }
               else if (line_string.find("pH") != string::npos)
               {
                  pH_found = 1;
                  //                    cout << " pH found in GeoSys species" << endl;
                  if (line_string.find("charge") != string::npos) this->rcml_pH_charge=1;
               }
            }
            if (line_string.find("# temp")!=string::npos)
            {
               // check if heat transport process is calculated in GeoSys
               size_t no_processes (pcs_vector.size());
               for (size_t i = 0; i < no_processes; i++)
               {
                  m_pcs = pcs_vector[i];
                  //					if (m_pcs->pcs_type_name.find("HEAT_TRANSPORT") != string::npos) {
                  if (m_pcs->getProcessType() == HEAT_TRANSPORT)
                  {
                     this->rcml_heat_flag = 1;
                     break;
                  }
               }
            }
            if(line_string.find("temp")!=string::npos)
            {
               if(this->rcml_heat_flag < 1)       // if no heat transport is calculated
               {
                  in.str(line_string);
                                                  // save temperature
                  in >> speciesname >> this->temperature;
                  in.clear();
               }
            }
         }
         this->rcml_number_of_master_species = n_master_species;

         // Handle pH, H+ and pe
         speciesname = "pH";
         m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
                                                  // old timelevel
         idx = m_pcs->GetNodeValueIndex(speciesname)+1;
         pqc_names.push_back(speciesname);        //store species name in vector names
         pqc_process.push_back(m_pcs->pcs_number);// store process name in vector processes
         pqc_index.push_back(idx);
         // check, if H+ is a GEoSys transport process
         speciesname = "H+";
         m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
         if(m_pcs != NULL)
         {
                                                  // old timelevel
            idx = m_pcs->GetNodeValueIndex(speciesname)+1;
            pqc_names.push_back(speciesname);     //store species name in vector names
                                                  // store process name in vector processes
            pqc_process.push_back(m_pcs->pcs_number);
            pqc_index.push_back(idx);
            this->gamma_Hplus =1.0;
         }
         else
         {
            // Store dummy values anyway, as otherwise indexing in Write becomes too complicated
            pqc_names.push_back(speciesname);
            pqc_process.push_back(-1);
            pqc_index.push_back(-1);
         }
         // Treat pe
         speciesname = "pe";
         m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
                                                  // old timelevel
         idx = m_pcs->GetNodeValueIndex(speciesname)+1;
         pqc_names.push_back(speciesname);        //store species name in vector names
         pqc_process.push_back(m_pcs->pcs_number);// store process name in vector processes
         pqc_index.push_back(idx);
      }
      /* Schleife ueber Keyword EQUILIBRIUM PHASES */
                                                  // keyword found
      if(line_string.find("EQUILIBRIUM_PHASES")!=string::npos)
      {
         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile->getline(line,MAX_ZEILE);
            line_string = line;
            if (line_string.find("# comp")!=string::npos)
            {
               n_equi_phases +=1;
               in.str(line_string);
               in >> speciesname ;
               in.clear();
               m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
                                                  // old timelevel
               idx = m_pcs->GetNodeValueIndex(speciesname)+1;
               pqc_names.push_back(speciesname);  //store species name in vector names
                                                  // store process name in vector processes
               pqc_process.push_back(m_pcs->pcs_number);
               pqc_index.push_back(idx);

            }
         }
         this->rcml_number_of_equi_phases = n_equi_phases;
      }

      /* Schleife ueber Keyword EXCHANGE */
                                                  // keyword found
      if(line_string.find("EXCHANGE")!=string::npos)
      {

         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile->getline(line,MAX_ZEILE);
            line_string = line;
            if (line_string.find("# comp") !=string::npos)
            {
               n_ion_exchange +=1;
               in.str(line_string);
               in >> speciesname ;
               in.clear();
               m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
                                                  // old timelevel
               idx = m_pcs->GetNodeValueIndex(speciesname)+1;
               pqc_names.push_back(speciesname);  //store species name in vector names
                                                  // store process name in vector processes
               pqc_process.push_back(m_pcs->pcs_number);
               pqc_index.push_back(idx);
            }
         }
         this->rcml_number_of_ion_exchanges = n_ion_exchange;
      }
      /* Scleife über Keyword GAS_PHASE */
                                                  // keyword found
      if(line_string.find("GAS_PHASE")!=string::npos)
      {

         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile->getline(line,MAX_ZEILE);
            line_string = line;
            if (line_string.find("# comp") !=string::npos)
            {
               n_gas_species +=1;
               in.str(line_string);
               in >> speciesname ;
               in.clear();
               m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
                                                  // old timelevel
               idx = m_pcs->GetNodeValueIndex(speciesname)+1;
               pqc_names.push_back(speciesname);  //store species name in vector names
                                                  // store process name in vector processes
               pqc_process.push_back(m_pcs->pcs_number);
               pqc_index.push_back(idx);
            }
         }
         this->rcml_number_of_gas_species = n_gas_species;
      }
      /* Schleife ueber Keyword SELECTED_OUTPUT */
                                                  // keyword found
      if(line_string.find("SELECTED_OUTPUT")!=string::npos)
      {
         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile->getline(line,MAX_ZEILE);
            line_string = line;
            if (line_string.find("-file") !=string::npos)
            {
               in.str(line_string);
               in >> dummy >> this->results_file_name;
            }
         }
      }

   }                                              /*end while */

   if((this->gamma_Hplus > 0)) cout << " pH found and H+ found, pH for PHREEQC input is calculated using gamma_H+ and [H+] " << gamma_Hplus << endl;

   nj=rcml_number_of_master_species + rcml_number_of_equi_phases + rcml_number_of_ion_exchanges + rcml_number_of_gas_species;
   cout << " Found in *.pqc file: " << rcml_number_of_master_species << " master species (excluding pH, H+ and pe), " ;
   cout << rcml_number_of_equi_phases << " equilibrium phases,  " << rcml_number_of_ion_exchanges << " ion exchangers and " ;
   cout << rcml_number_of_gas_species << " gas species " << endl;
   //    for(i=0; i< (int) pqc_names.size();i++)
   //        cout << pqc_names[i] << ", " << pqc_index[i] << ", " << pqc_process[i] << endl;

   return 1;
}


/**************************************************************************
   ROCKFLOW - Funktion: ReadInputPhreeqc

   Aufgabe:
   Liest die Eingabedatei *.pqc und erstellt aktuelles PHREEQC-Eingabefile
                                                                          */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *dateiname: Dateiname ohne Extension
                                                                          */
/* Ergebnis:
   0 bei Fehler oder Ende aufgrund Dateitest, sonst 1

   Programmaenderungen:
   06/2003     SB         Erste Version
   06/2003     MX         Read function realisation
**************************************************************************/
int REACT::ReadInputPhreeqc(long index, FILE *fpqc, FILE *Fphinp)
{

   int i=0, np;
   int nn=0,nep=0, nj=0, nk=0;
   int beginn=0, p=0, found;
   int pH_flag=-1, pe_flag=-1, iheat;
   int pos,j;
   double dvalue, d_help;
   char str[256],sub[256], s[256];
   char *s1, *s2;

   s1 = s2 = NULL;
   FILE *indatei=NULL, *f;
   np = this->number_of_comp;
   // crdat = (char *) Malloc((int)(int)strlen(file_name));
   // crdat = strcat(strcpy(crdat,file_name),CHEM_REACTION_EXTENSION);  /*MX:0603*/

   /* Open input and output file */
   // indatei = fopen(crdat, "r");   /*input dateien*/
   indatei = fpqc;                                /*input dateien*/
   f=Fphinp;                                      /*output dateien*/

   if(indatei == NULL)
   {
      DisplayMsgLn("Erro:The input file *.pqc doesn't exist!!!");
      return 0;
   }

   /* zeilenweise lesen */
   while(fgets(str,256,indatei))
   {
      //    DisplayMsgLn("");
      //    DisplayMsgLn(str);

      pos=0;
      beginn=1;
      p=0;
      while ((!strstr(str, "END")) && (!strstr(str, "#ende")) && StrOnlyReadStr(sub, str, f, TFString, &p))
      {
         LineFeed(f);
         found=0;

         /* Schleife ueber Keyword Solution */
         /*-------------------------------------------------------------------------------------*/
         if (!strcmp(sub, "SOLUTION"))
         {
            found=1;
            FilePrintString(f, "SOLUTION  ");
            FilePrintInt(f, index+1);
            LineFeed(f);

            FilePrintString(f, "#GRID  ");
            FilePrintInt(f, index+1);
            LineFeed(f);

            //           while (fgets(str,256,indatei) && (!StrTestDollar(str, &pos))){
            while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
            {
               if ((!strstr(str, "# comp")) && (!strstr(str, "# temp")))
               {
                  //               DisplayMsgLn(" #ende and # not found");
                  FilePrintString(f, str);

               }
               else
               {
                  //                DisplayMsgLn(" # comp found ");
                  //sscanf(str, "%s");
                  StrReadStr(s, str, f, TFString, &pos);

                  /* SB: temperature introduced */
                  if(strcmp(s,"temp")==0)
                  {
                     if(rcml_heat_flag == 1)
                     {
                        iheat = np;
                        //                        printf("Temperature is %lf", val_in[iheat][index]);
                        /* convert temperature from Kelvin to degree celsius for PHREEQC */
                                                  //MX 03.05
                        if (val_in[iheat][index]<273.0) val_in[iheat][index] += 273.15;
                        FilePrintDouble(f, val_in[iheat][index] - 273.15);
                        FilePrintString(f, " # temp ");
                        LineFeed(f);
                     }
                     /* else write line to file as read in */
                     else
                     {
                        FilePrintString(f, str);
                     };

                  };
                  /* SB: end temperature */

                  /* find the concentration of each component from val_in*/
                  for (j=0;j<np;j++)
                  {
                     if (strcmp(s,name[j])==0)
                     {
                        if(strcmp("pH",name[j])!=0)
                        {
                           nn += 1;
                           d_help = val_in[j][index];
                           //                        printf("Component %s: concentration %lf", s, d_help);
                           if(fabs(d_help) < 1.0e-019) d_help = 0.0;
                           FilePrintDouble(f, d_help);
                           FilePrintString(f, " # comp ");
                           FilePrintInt(f, j+1);
                           LineFeed(f);
                           if (strcmp("pe",name[j])==0)
                              pe_flag=0;
                           break;
                        }
                        /*		                        else{  // pH in RF input file
                                                pH_flag=0;
                                                printf("Component %s: concentration %lf", s, val_in[j][index]);
                                                FilePrintDouble(f, val_in[j][index]);
                                                FilePrintString(f, " charge ");
                                                FilePrintString(f, " # comp ");
                                                FilePrintInt(f, j+1);
                                                LineFeed(f);
                                                break;
                                             }
                        */
                     }

                     else if (strcmp(s,"pH")==0)  /* if pH will change with the reaction! */
                     {
                        for (i=0;i<np;i++)
                        {
                           if (strcmp(s,name[i])==0)
                           {
                                                  /* pH in RF input file */
                              if(strcmp("pH",name[i])==0)
                              {

                                 pH_flag=0;
                                 //                           printf("Component %s: concentration %lf", s, val_in[i][index]);
                                 FilePrintDouble(f, val_in[i][index]);
                                 //MX						   FilePrintString(f, " charge ");

                                 if (strstr(str, "charge")) FilePrintString(f, " charge ");
                                 FilePrintString(f, " # comp ");
                                 FilePrintInt(f, rcml_number_of_master_species+1);
                                 LineFeed(f);
                                 break;
                              }
                           }
                        }                         /*end for*/

                        if (pH_flag<0)
                        {

                           DisplayMsgLn("pH is not included in the transport but will be calculated in the reaction");
                           pH_flag=1;
                           p=0;
                           StrReadDouble(&dvalue, &str[p+=pos], f, TFDouble, &pos);
                           StrReadStr(s, &str[p+=pos], f, TFString, &pos);
                           FilePrintString(f, " # comp ");
                           FilePrintInt(f, rcml_number_of_master_species+1);
                           LineFeed(f);
                           name[rcml_number_of_master_species+1]="pH";
                           break;
                        }
                        break;
                     }

                     else if (strcmp(s,"pe")==0)  /* if pe will change with the reaction! */
                     {
                        for (i=0;i<np;i++)
                        {
                           if (strcmp(s,name[i])==0)
                           {
                                                  /* pe in RF input file */
                              if(strcmp("pe",name[i])==0)
                              {
                                 pe_flag=0;
                                 //                           printf("Component %s: concentration %lf", s, val_in[i][index]);
                                 FilePrintDouble(f, val_in[i][index]);
                                 FilePrintString(f, " # comp ");
                                 FilePrintInt(f, rcml_number_of_master_species+2);
                                 LineFeed(f);
                                 break;
                              }
                           }
                        }                         /*end for*/

                        if (pe_flag<0)
                        {
                           DisplayMsgLn("pe is not included in the transport but will be calculated in the reaction");
                           pe_flag=1;
                           StrReadDouble(&dvalue, &str[pos], f, TFDouble, &pos);
                           FilePrintString(f, " # comp ");
                           FilePrintInt(f, rcml_number_of_master_species+2);
                           LineFeed(f);
                           break;
                        }
                        break;
                     }

                  }                               /*end-for*/
               }                                  /*end-else*/
            }                                     /*end while -3*/
            if (nn != rcml_number_of_master_species)
            {
               FilePrintString(f, "Warnung:Concentration of master species not found in *.rfd file!!!");
            }
            pos=0;
            beginn=1;
            p=0;
         }                                        /*if_SOLUTION*/

         /* Schleife ueber Keyword SOLUTION_SPECIES */
         /*-------------------------------------------------------------------------------------*/
         else if (!strcmp(sub, "SOLUTION_SPECIES"))
         {
            found=1;
            FilePrintString(f, "SOLUTION_SPECIES  ");
            FilePrintInt(f, index+1);
            LineFeed(f);
            while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
            {
               FilePrintString(f, str);
            }
         }

         /* Schleife ueber Keyword SOLUTION_MASTER_SPECIES */
         /*-------------------------------------------------------------------------------------*/
         else if (!strcmp(sub, "SOLUTION_MASTER_SPECIES"))
         {
            found=1;
            FilePrintString(f, "SOLUTION_MASTER_SPECIES  ");
            FilePrintInt(f, index+1);
            LineFeed(f);
            while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
            {
               FilePrintString(f, str);
            }
         }

         /* Schleife ueber Keyword EQUILIBRIUM_PHASES */
         /*-------------------------------------------------------------------------------------*/
         else if (!strcmp(sub, "EQUILIBRIUM_PHASES"))
         {
            found=1;
            FilePrintString(f, "EQUILIBRIUM_PHASES  ");
            FilePrintInt(f, index+1);
            LineFeed(f);
            while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
            {
               pos=0;
               p=0;
               StrReadStr(s, &str[p+=pos], f, TFString, &pos);
               if (strcmp(s,"")==0) break;
               StrReadDouble(&dvalue, &str[pos], f, TFDouble, &pos);

               /* find the mass of each phase from val_in*/
               for (i=0;i<np;i++)
               {
                  if (strcmp(s,name[i])==0)
                  {
                     nep += 1;
                     //                        printf("  Phase %s:  %lf \n", s, val_in[i][index]);
                     FilePrintDouble(f, val_in[i][index]);
                     FilePrintString(f, " # comp ");
                     FilePrintInt(f, i+1);
                     LineFeed(f);
                     break;
                  }
               }                                  /*end for*/
            }                                     /*end while*/
            if (nep != rcml_number_of_equi_phases)
            {
               FilePrintString(f, "Warning: One or more solid phase(s) in EQUILIBRIUM_PHASES not found in *.rfd file!!!");
               //                    FilePrintString(f, "Value(s) in *.pqc file used!!!");
               //                    FilePrintString(f, str);
            }
         }                                        /*end if*/

         /* Schleife ueber Keyword KINETICS */
         /*-------------------------------------------------------------------------------------*/
         else if (!strcmp(sub, "KINETICS"))
         {
            found=1;
            FilePrintString(f, "KINETICS  ");
            FilePrintInt(f, index+1);
            LineFeed(f);
            while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
            {
               if(strstr(str,"-steps"))
               {
                  /* zerlegen und RockFlow timestep einfügen */
                                                  ///OK
                  sscanf(str," %s %lf %s %i %s", s1,&dvalue,s,&j,s2);
                  /* steps in reaction model eintragen */
                  this->rcml_number_of_pqcsteps = j;

                  if(index > -1)
                  {
                     FilePrintString(f, " -steps ");
                     CTimeDiscretization *m_tim = NULL;
                     if(time_vector.size()>0)
                        m_tim = time_vector[0];
                     else
                        cout << "Error in MPCCalcCharacteristicNumbers: no time discretization data !" << endl;
                     FilePrintDouble(f,m_tim->CalcTimeStep());
                     //OK_TIM	              FilePrintDouble(f, GetDt(aktueller_zeitschritt-1l));
                     FilePrintString(f, " in ");
                     FilePrintInt(f, j);
                     FilePrintString(f, " steps ");
                     LineFeed(f);
                  }

               }
               else
                  FilePrintString(f, str);
            }
         }

         /* Schleife ueber Keyword RATES */
         /*-------------------------------------------------------------------------------------*/
         else if (!strcmp(sub, "RATES"))
         {
            found=1;
            FilePrintString(f, "RATES  ");
            FilePrintInt(f, index+1);
            LineFeed(f);
            while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
            {
               FilePrintString(f, str);
            }
         }

         /*  Keyword PHASE */
         /*-------------------------------------------------------------------------------------*/
         else if (!strcmp(sub, "PHASE"))
         {
            found=1;
            FilePrintString(f, "PHASE  ");
            FilePrintInt(f, index+1);
            LineFeed(f);
            while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
            {
               FilePrintString(f, str);
            }
         }

         /*  Keyword EXCHANGE */
         /*-------------------------------------------------------------------------------------*/
         else if (!strcmp(sub, "EXCHANGE"))
         {
            found=1;
            FilePrintString(f, "EXCHANGE  ");
            FilePrintInt(f, index+1);
            LineFeed(f);
            while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
            {
               pos=0;
               p=0;
               StrReadStr(s, &str[p+=pos], f, TFString, &pos);

               /* find the mass of each phase from val_in*/
               for (i=0;i<np;i++)
               {
                  if (strcmp(s,name[i])==0)
                  {
                     nj += 1;
                     //                        printf("  Phase %s:  %lf \n", s, val_in[i][index]);
                     FilePrintDouble(f, val_in[i][index]);
                     FilePrintString(f, " # comp ");
                     FilePrintInt(f, i+1);
                     LineFeed(f);
                     break;
                  }
               }                                  /*end for*/
            }                                     /*end while*/
            if (nj != rcml_number_of_ion_exchanges)
            {
               FilePrintString(f, "Warning: Problem in reading EXCHANGE  in *.rfd or *.pqc file!!!");
               //                    FilePrintString(f, "Value(s) in *.pqc file used!!!");
               //                    FilePrintString(f, str);
            }
         }                                        /*end if*/

         /*  Keyword PRINT */
         /*-------------------------------------------------------------------------------------*/
         else if (!strcmp(sub, "PRINT"))
         {
            found=1;
            if (index==0)
            {
               FilePrintString(f, "PRINT  ");
               LineFeed(f);
               while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
               {
                  FilePrintString(f, str);
               }
            }
            else
               while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
            {
            }
         }
         /*  Keyword SELECTED_OUTPUT */
         /*-------------------------------------------------------------------------------------*/
         else if (!strcmp(sub, "SELECTED_OUTPUT"))
         {
            found=1;
            if (index==0)
            {
               FilePrintString(f, "SELECTED_OUTPUT  ");
               LineFeed(f);
               while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
               {
                  FilePrintString(f, str);
                  p=0;
                  StrOnlyReadStr(sub, str, f, TFString, &p);
                  if (!strcmp(sub, "-file"))
                  {
                     StrOnlyReadStr(sub, &str[p], f, TFString, &p);
                     /* SB: moved to structure rcml
                                         strcpy(fsout, sub);
                     */
                  }
               }
            }
            else
               while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
            {
            }
         }

         /*  Keyword USER_PUNCH */
         /*-------------------------------------------------------------------------------------*/
         else if (!strcmp(sub, "USER_PUNCH"))
         {
            found=1;
            if (index==0)
            {
               FilePrintString(f, "USER_PUNCH  ");
               LineFeed(f);

               /*------------head------------------------------*/
               fprintf(f, " -head");
               for (j=0; j<rcml_number_of_master_species; j++)
               {
                  if((strcmp("pH",name[j])!=0) &&(strcmp("pe",name[j])!=0))
                     fprintf(f, " %s", name[j]);
               }
               fprintf(f, " pH pe");

               nj=rcml_number_of_master_species;
               nk=rcml_number_of_equi_phases;

               for (j=nj+2; j<nj+rcml_number_of_equi_phases+2; j++)
               {
                  fprintf(f, " %s", name[j]);
               }

               for (j=nj+nk+2; j<nj+nk+rcml_number_of_ion_exchanges+2; j++)
               {
                  fprintf(f, " %s", name[j]);
               }
               LineFeed(f);

               /*------------master speciese---------------------*/
               fprintf(f, " 10 PUNCH");
               for (j=0; j<nj-1; j++)
               {
                  if((strcmp("pH",name[j])!=0) && (strcmp("pe",name[j])!=0))
                     fprintf(f, " TOT(%c%s%c),", 34, name[j], 34);
               }
               fprintf(f, " TOT(%c%s%c)", 34, name[nj-1], 34);

               /*------------pH pe-------------------------------*/
               LineFeed(f);
               fprintf(f, " 20 PUNCH");
               fprintf(f, " -LA(%c%s%c),", 34, "H+", 34);
               fprintf(f, " -LA(%c%s%c)", 34, "e-", 34);

               /*------------equilibrium phases-------------------*/
               if (rcml_number_of_equi_phases>0)
               {
                  LineFeed(f);
                  fprintf(f, " 40 PUNCH");
                  for (j=nj+2; j<nj+rcml_number_of_equi_phases+1; j++)
                  {
                     fprintf(f, " EQUI(%c%s%c),", 34, name[j], 34);
                  }
                  fprintf(f, " EQUI(%c%s%c),", 34, name[nj+rcml_number_of_equi_phases+1], 34);
               }

               /*------------exchange-------------------*/
               if (rcml_number_of_ion_exchanges>0)
               {
                  LineFeed(f);
                  fprintf(f, " 60 PUNCH");
                  for (j=nj+nk+2; j<nj+nk+rcml_number_of_ion_exchanges+1; j++)
                  {
                     fprintf(f, " MOL(%c%s%c),", 34, name[j], 34);
                  }
                  fprintf(f, " MOL(%c%s%c),", 34, name[nj+nk+rcml_number_of_ion_exchanges+1], 34);
               }
               while (fgets(str,256,indatei) && (!strstr(str, "#ende")))
               {
                  //                     FilePrintString(f, str);
               }
            }
            else
               while (fgets(str,256,indatei) && ((!strstr(str, "#ende"))||(!strstr(str, "END"))))
            {
            }
         }                                        /* end if_USER_PUNCH */
         /* SB: added: time discretisation: Keyword is "-steps" */
         else if (!strcmp(sub, "-steps"))
         {
            found = 1;
                                                  ///OK
            sscanf(str," %s %lf %s %i %s", s1, &dvalue,s,&j,s2);
            if(index > -1)
            {
               FilePrintString(f, " -steps ");
               CTimeDiscretization *m_tim = NULL;
               if(time_vector.size()>0)
                  m_tim = time_vector[0];
               else
                  cout << "Error in MPCCalcCharacteristicNumbers: no time discretization data !" << endl;
               FilePrintDouble(f,m_tim->CalcTimeStep());
               //OK_TIM              FilePrintDouble(f, GetDt(aktueller_zeitschritt-1l));
               FilePrintString(f, " in ");
               FilePrintInt(f, j);
               FilePrintString(f, " steps ");
               LineFeed(f);
            }
         }
         /* Not found */
         else
         {
            if (found !=1)
            {
               fprintf(f, " %s %s ", str, " in the *.pqc file is an unknown keyword!!!");
               exit(1);
            }
         }
      }                                           /* end while - 2 */
   }                                              /* end while - 1 */

   LineFeed(f);
   fprintf(f, "END");
   LineFeed(f);
   LineFeed(f);
   //    if (indatei !=NULL) fclose(indatei);

   return 1;
}                                                 /* end-if */


/**************************************************************************
   ROCKFLOW - Funktion: WriteInputPhreeqc

   Aufgabe:
   Liest die Eingabedatei *.pqc und erstellt aktuelles PHREEQC-Eingabefile
                                                                          */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *dateiname: Dateiname ohne Extension
                                                                          */
/* Ergebnis:
   0 bei Fehler oder Ende aufgrund Dateitest, sonst 1

   Programmaenderungen:
   06/2003     SB         Erste Version
   06/2003     MX         Read function realisation
   01/2006     SB         Reimplementation, C++ class, IO, bugfixes
**************************************************************************/
int REACT::WriteInputPhreeqc(long index, /*ifstream *pqc_iinfile,*/ ofstream *out_file){

char line[MAX_ZEILE];
std::stringstream in;
string name, line_string, speciesname, dummy;
CRFProcess* m_pcs = NULL;
int i, ii, idx, n1, n2, n3, n4, count=-1;
double dval, dval1, sat_index=0.0;
double z, h, dens, press, partial_press, volume, temp=-1.0, mm;

//  cout << " WriteInputPhreeqc for node " << index << endl;
cout.flush();
/* File handling - rewind file */
// pqc_infile->seekg(0L,ios_base::beg);
//  pqc_infile->close();
ifstream pqc_infile (this->file_name_pqc.data(),ios::in);
pqc_infile.seekg(0L,ios::beg);

// precision output file
out_file->setf(ios::scientific,ios::floatfield);
out_file->precision(12);

/* zeilenweise lesen */
while(!pqc_infile.eof())
{
   pqc_infile.getline(line,MAX_ZEILE);
   line_string = line;
   if(line_string.find("#STOP")!=string::npos)
      break;
   //-------------------------------------------------------------------------------------------------------------
   /* Schleife ueber Keyword Solution */
   if(line_string.find("SOLUTION")!=string::npos) // keyword found
   {
      *out_file << "SOLUTION " << index+1 << " #New Version " << endl;
      *out_file << "#GRID " << index+1 << endl;
      while(line_string.find("#ende")==string::npos)
      {
         pqc_infile.getline(line,MAX_ZEILE);
         line_string = line;
         if (line_string.find("# comp")!=string::npos)
         {
            if (line_string.find("pH") == string::npos && line_string.find("pe") == string::npos)
            {
               // Component found; write name and concentration of component
               count++;
               /*                in.str(line_string);
                               in >> speciesname ;
                               m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
                               idx = m_pcs->GetNodeValueIndex(speciesname)+1; // old timelevel
                               dval = m_pcs->GetNodeValue(index,idx);
               */
               speciesname = pqc_names[count];
               dval = pcs_vector[pqc_process[count]]->GetNodeValue(index,pqc_index[count]);
               /*                cout << "Testing index vectors: " << speciesname << ", " << m_pcs->pcs_number << ", "<< idx <<",  With vectors: " << pqc_names[count] << ", " << pcs_vector[pqc_process[count]]->pcs_number << ", " << pqc_index[count];
                               cout << " Values: " << dval << ", " << dval1 << endl;
               */
               if(speciesname.compare("pe"))      // if this is not pe
                  if(dval < 1.0e-19) dval = 0.0;
               //                if(speciesname.compare("pH"))
               *out_file << speciesname << "       " << dval << "     # comp " <<endl;
               //					if(index <2) cout << speciesname << " " << dval << endl;
               //                else
               //                    *out_file << speciesname << "       " << dval << " charge " << "       # comp " <<endl;
               //                in.clear();
            }
         }
         else if (line_string.find("# temp")!=string::npos)
         {
            // check if heat transport process is calculated in GeoSys
            if(this->rcml_heat_flag > 0)
            {
               m_pcs = PCSGet("HEAT_TRANSPORT");
               idx = m_pcs->GetNodeValueIndex("TEMPERATURE1");
               dval = m_pcs->GetNodeValue(index,idx);
               if (dval<273.0) dval += 273.15;    //change from °C to Kelvin if necessary
               dval -=273.15;                     // Input to PHREEQC is in °C
               *out_file << "temp " << dval << "  # temp " << endl;
               temp = dval;                       // save for gas phase input
            }
         }
         else                                     // Write units and temperature in the standard case
         {
            if (line_string.find("pH") == string::npos && line_string.find("pe") == string::npos && line_string.find("#ende") == string::npos)
               *out_file << line_string << endl;
         }
      }                                           // end while

      // special treat pH, and pe
      n1 = this->rcml_number_of_master_species;
      count++;
      if(count != n1) cout << "Error in index of pqc_vectors !" << endl;
      dval = pcs_vector[pqc_process[count]]->GetNodeValue(index,pqc_index[count]);
      //		if(index <2)         cout << " pH: " << dval;
      // H+
      count++;
      if(this->gamma_Hplus > 0)                   // pH and H+ in GeoSys species, calculate pH from H+
      {
         dval1 = pcs_vector[pqc_process[n1+1]]->GetNodeValue(index,pqc_index[n1+1]);
         dval = -log10(dval1*gamma_Hplus);
         //            if(index<2) cout << " .  Resetting pH to: " << dval << "; MOL[H+]= " << dval1 << ", gamma_H+ = " << gamma_Hplus;
      }
      if(this->rcml_pH_charge > 0)
         *out_file << "pH" << "       " << dval << " charge " << "       # comp " <<endl;
      else
         *out_file << "pH" << "       " << dval << "       # comp " <<endl;
      //SB to do screen output				if(index <2) cout << "  pH: " << dval << ", " << pcs_vector[pqc_process[count]]->pcs_number << endl;
      // write pe
      count++;
      dval = pcs_vector[pqc_process[n1+2]]->GetNodeValue(index,pqc_index[n1+2]);
      *out_file << "pe" << "       " << dval << "       # comp " <<endl;
      //SB to do screen output		if(index <2)  cout << "  pe: " << dval << ", " << pcs_vector[pqc_process[count]]->pcs_number << endl;

      *out_file << line_string << endl;
   }                                              // end SOLUTION
   //-------------------------------------------------------------------------------------------------------------
   /* Schleife ueber Keyword EQUILIBRIUM PHASES */
                                                  // keyword found
   if(line_string.find("EQUILIBRIUM_PHASES")!=string::npos)
   {
      *out_file << endl << "EQUILIBRIUM_PHASES   " << index+1 << endl;
      while(line_string.find("#ende")==string::npos)
      {
         pqc_infile.getline(line,MAX_ZEILE);
         line_string = line;
         if (line_string.find("# comp")!=string::npos)
         {
            count++;
            in.str(line_string);
            in >> speciesname >> sat_index ;

            /*                m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
                            idx = m_pcs->GetNodeValueIndex(speciesname)+1; // old timelevel
                            dval = m_pcs->GetNodeValue(index,idx);
            */
            //                if(count != (n1+3)) cout << " Error in index pqc " << endl;
            speciesname = pqc_names[count];
            dval = pcs_vector[pqc_process[count]]->GetNodeValue(index,pqc_index[count]);
            if(dval < 1.0e-19) dval = 0.0;

            //                cout << "Testing index vectors: " << speciesname << ", " << m_pcs->pcs_number << ", "<< idx <<",  With vectors: " << pqc_names[count] << ", " << pcs_vector[pqc_process[count]]->pcs_number << ", " << pqc_index[count];
            //                if(index <2) cout << " EQ-Species " << speciesname << " " << dval << endl;
            *out_file << speciesname << " " << sat_index << "  " << dval << "       # comp " <<endl;
            in.clear();
         }
         else
            *out_file << line_string << endl;
      }
   }                                              // end EQUILIBRIUM PHASES
   //-------------------------------------------------------------------------------------------------------------
   /* Schleife ueber Keyword EXCHANGE */
   if(line_string.find("EXCHANGE")!=string::npos) // keyword found
   {
      *out_file << endl << "EXCHANGE   " <<  index+1 << endl;
      while(line_string.find("#ende")==string::npos)
      {
         pqc_infile.getline(line,MAX_ZEILE);
         line_string = line;
         if (line_string.find("# comp") !=string::npos)
         {
            count++;
            /*
                            in.str(line_string);
                            in >> speciesname ;
                            m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
                            idx = m_pcs->GetNodeValueIndex(speciesname)+1; // old timelevel
                            dval = m_pcs->GetNodeValue(index,idx);
            */
            speciesname = pqc_names[count];
            dval = pcs_vector[pqc_process[count]]->GetNodeValue(index,pqc_index[count]);
            if(dval < 1.0e-19) dval = 0.0;
            /*                cout << "Testing index vectors: " << speciesname << ", " << m_pcs->pcs_number << ", "<< idx <<",  With vectors: " << pqc_names[count] << ", " << pcs_vector[pqc_process[count]]->pcs_number << ", " << pqc_index[count];
                            cout << " Values: " << dval << ", " << dval1 << endl;
            */
            *out_file << speciesname << "       " << dval << "       # comp " <<endl;
            //                in.clear();
         }
         else
            *out_file << line_string << endl;
      }
   }                                              // end EXCHANGE
   //-------------------------------------------------------------------------------------------------------------
   /* Schleife ueber Keyword GAS_PHASE */
   if(line_string.find("GAS_PHASE")!=string::npos)// keyword found
   {
      *out_file << endl << "GAS_PHASE   " <<  index+1 << endl;

      // get necessary values for conversion of molar concentrations to partial pressures, and to calculate total pressure and total volume

      // get height of node z
      CFEMesh* m_msh = fem_msh_vector[0];         //SB: ToDo hart gesetzt
      Mesh_Group::CNode* m_nod = NULL;
      m_nod = m_msh->nod_vector[index];
      z = m_msh->nod_vector[index]->Z();

      // get piezometric hight h
      m_pcs = PCSGet("GROUNDWATER_FLOW");
      if(m_pcs == NULL) cout << "   Error - no flow process found!" << endl;
      idx = m_pcs->GetNodeValueIndex("HEAD")+1;
      h = m_pcs->GetNodeValue(index,idx);

      // get fluid density
      dens = mfp_vector[0]->Density();

      // calculate pressure in [Pa]
      press = dens * gravity_constant * (h-z);
      // cout << " Pressure: " << press << " = " << dens << " * " << gravity_constant << " * ( " << h << " - " << z << " ) " << endl;

      // get temperature in [°C]
      if(rcml_heat_flag < 1) temp = this->temperature;

      // get molar masses of gas phase
      mm=0.0;                                     // mm is total molar mass of gas phase in [mol]
      ii= rcml_number_of_master_species + 3 + rcml_number_of_ion_exchanges + rcml_number_of_equi_phases;
      for(i=ii;i<ii+rcml_number_of_gas_species;i++)
      {
         speciesname = this->pqc_names[i];
         //			cout << "Testing index vectors: " << speciesname << ",   With vectors: " << pqc_names[i] << ", " << pcs_vector[pqc_process[i]]->pcs_number << ", " << pqc_index[i];
         dval = pcs_vector[pqc_process[i]]->GetNodeValue(index,pqc_index[i]);
         //			cout << dval << endl;
         mm += dval;
      }

      //  calculate Volume of gas phase in [mol * Pa * m^3 / K / mol * K / Pa = m^3 ]
      volume = mm * 8.314472 * (273.15 + temp) /press;

      while(line_string.find("#ende")==string::npos)
      {
         pqc_infile.getline(line,MAX_ZEILE);
         line_string = line;
         if (line_string.find("-pressure") !=string::npos)
                                                  // pressure in atmospheres
            *out_file << "        -pressure       " << press/101325.0 << endl;
         else if (line_string.find("-volume") !=string::npos)
                                                  // volume in Liters
            *out_file << "        -volume       " << volume*1000.0 << endl;
         else if (line_string.find("-temperature") !=string::npos)
                                                  // temperature in °Celsius
            *out_file << "        -temperature       " << temp << endl;
         else if (line_string.find("# comp") !=string::npos)
         {
            count++;
            speciesname = pqc_names[count];
            dval = pcs_vector[pqc_process[count]]->GetNodeValue(index,pqc_index[count]);
            if(dval < 1.0e-19) dval = 0.0;
            // cout << "Testing index vectors: " << speciesname << " , with vectors: " << pqc_names[count] << ", " << pcs_vector[pqc_process[count]]->pcs_number << ", " << pqc_index[count];
            // cout << " Molar mass: " << dval  << ", component / total mass: " << dval/mm << endl;
            // dval is molar mass of gas component, convert to partial pressure [Pa]
            if (mm > 0.0)
               partial_press = press * dval/mm;
            else
               partial_press = 0.0;
            // cout << " partial pressure: " << partial_press  << ", partial/total press: " << (partial_press/press) << endl;

            *out_file << "        " << speciesname << "       " << partial_press/101325.0 << "       # comp " <<endl;
         }
         else
            *out_file << line_string << endl;     //write line unchanged
      }
   }                                              // end GAS_PHASE
   //-------------------------------------------------------------------------------------------------------------
   /* Schleife ueber Keyword SELECTED_OUTPUT */
                                                  // keyword found
   if(line_string.find("SELECTED_OUTPUT")!=string::npos)
   {
      if(index < 1)
      {
         *out_file << endl << "SELECTED_OUTPUT" << endl;
         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile.getline(line,MAX_ZEILE);
            line_string = line;
            if (line_string.find("-file") !=string::npos)
               *out_file << "-file " << this->results_file_name << endl;
            else
               *out_file << line_string << endl;
         }
      }
   }                                              // end SELECTED OUTPUT
   //-------------------------------------------------------------------------------------------------------------
   /* Schleife ueber Keyword PRINT */
   if(line_string.find("PRINT")!=string::npos)    // keyword found
   {
      if(index < 1)
      {
         *out_file << endl << "PRINT" << endl;
         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile.getline(line,MAX_ZEILE);
            line_string = line;
            *out_file << line_string << endl;
         }
      }
   }                                              // end PRINT
   //-------------------------------------------------------------------------------------------------------------
   /* Schleife ueber Keyword USER_PUNCH */
                                                  // keyword found
   if(line_string.find("USER_PUNCH")!=string::npos)
   {
      if(index < 1)
      {
         *out_file << endl << "USER_PUNCH" << endl;
         // Write Header
         n1 = this->rcml_number_of_master_species;
         n2 = this->rcml_number_of_equi_phases;
         n3 = this->rcml_number_of_ion_exchanges;
         n4 = this->rcml_number_of_gas_species;
         *out_file << "-head ";
         for(i=0;i<n1; i++) *out_file << " " << pqc_names[i];
         *out_file << " pH ";
         *out_file << " H+ ";
         *out_file << " pe ";
         for(i=n1+3; i<n1+3+n2; i++) *out_file << " " << pqc_names[i];
         for(i=n1+3+n2;i<n1+3+n2+n3; i++) *out_file << " " << pqc_names[i];
         for(i=n1+3+n2+n3; i<n1+3+n2+n3+n4; i++) *out_file << " " << pqc_names[i];
         *out_file << endl;
         // Write master species
         *out_file << " 10 PUNCH ";
         for(i=0;i<n1;i++)
         {
            if(pqc_names[i].compare("H+")==0)
                                                  // extra treat H+
               *out_file << " MOL(\"" << pqc_names[i] << "\"),";
            else
                                                  // without pH and pe here
               *out_file << " TOT(\"" << pqc_names[i] << "\"),";
         }
         *out_file << endl;
         // Write pH and pe
         *out_file << " 20 PUNCH " << " -LA(\"H+\"), ";
         *out_file << " MOL(\"H+\"), ";
         *out_file << "  -LA(\"e-\")" << endl;
         // Write equilibrium phases
         if(n2 > 0)
         {
            *out_file << " 40 PUNCH ";
            for(i=n1+3;i<n1+3+n2;i++) *out_file << " EQUI(\"" << pqc_names[i] << "\"),";
            *out_file << endl;
         }
         // Write ion exchangers
         if(n3 > 0)
         {
            *out_file << " 60 PUNCH ";
            for(i=n1+3+n2;i<n1+3+n2+n3;i++) *out_file << " MOL(\"" << pqc_names[i] << "\"),";
            *out_file << endl;
         }
         // Write gas phase species
         if(n4 > 0)
         {
            *out_file << " 70 PUNCH ";
            for(i=n1+3+n2+n3; i<n1+3+n2+n3+n4; i++) *out_file << " GAS(\"" << pqc_names[i] << "\"),";
            *out_file << endl;
         }
      }                                           // end if index < 1

      // search for end of USER_PUNCH data block in *.pqc input file
      while(!pqc_infile.eof())
      {
         pqc_infile.getline(line,MAX_ZEILE);
         line_string = line;
         if((line_string.find("#ende")!=string::npos) || (line_string.find("END")!=string::npos))
            break;
      }
   }                                              // end USER_PUNCH
   //-------------------------------------------------------------------------------------------------------------
   if(line_string.find("-steps")!=string::npos)   // keyword found
   {
      in.str(line_string);
      in >> dummy >> dval >> this->rcml_number_of_pqcsteps >> dummy;
      CTimeDiscretization *m_tim = NULL;
      if(time_vector.size()>0)
         m_tim = time_vector[0];
      else
         cout << "Error in WriteInputPhreeqc: no time discretization data !" << endl;
      dval = m_tim->CalcTimeStep();
      *out_file << "-steps " << dval << " in " << this->rcml_number_of_pqcsteps << " steps" << endl;
   }                                              // end -steps
   //-------------------------------------------------------------------------------------------------------------
   if(line_string.find("KNOBS")!=string::npos)
   {
      if(index < 1)
      {
         *out_file << endl << "KNOBS" << endl;
         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile.getline(line,MAX_ZEILE);
            line_string = line;
            *out_file << line_string << endl;
         }
      }
   }

}                                                 /*end while zeilenweises lesen */


*out_file << "END" << endl << endl;

pqc_infile.close();
//    out_file.close();

return 1;
}


/**************************************************************************
   ROCKFLOW - Funktion: Call_Phreeqc

   Aufgabe:
   Ruft PHREEQC auf

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)

   Ergebnis:
   0 bei Fehler oder Ende aufgrund Dateitest, sonst 1

Programmaenderungen:
06/2003     SB         Erste Version
**************************************************************************/

int REACT::Call_Phreeqc(void)
{

   const char *m_phreeqc;
   //  m_phreeqc="phrqc phinp.dat  phinp.out  phreeqc.dat";
   m_phreeqc="phreeqc phinp.dat  phinp.out  phreeqc.dat";

#ifdef PHREEQC
   if (!system(m_phreeqc))
   {
      //    DisplayMsgLn("Phreeqc runs succesfully! ");
      return 1;
   }
   else
   {
      DisplayMsgLn("Warnung: Phreeqc doesn't run properly!!! ");
      exit(1);
   }
#endif

#ifndef PHREEQC
   return 1;
#endif
}


/**************************************************************************
   ROCKFLOW - Funktion: ReadOutputPhreeqc

   Aufgabe:
   Liest Ergebnisse der PHREEQC-Berechnungen aus PHREEQC-Ausdgabedatei

  Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *File: Dateiname der PHREEQC-Ausdgabedatei
   R char **val_out: Zeiger für chemische Conzentration
   E char **name: Componentnamen

Ergebnis:
0 bei Fehler oder Ende aufgrund Dateitest, sonst 1

Programmaenderungen:
06/2003     MX/SB         Erste Version
11/2003     SB            Bugfix for large files, read long lines
switched to iostreams
************************************************************************************************/
int REACT::ReadOutputPhreeqc(char* fout)
{

   int ok = 0;
   int nj,nk, ntot;
   int i, j, ii, zeilenlaenge=10000, anz;
   char str[4000];
   double dval;

   ifstream ein;
   ein.open(fout);
   if (! ein)
   {
      DisplayMsgLn("The selected output file doesn't exist!!!");
      return 0;
   }
   /* get total number of species in PHREEQC output file */
   ntot = rcml_number_of_master_species+2 + rcml_number_of_equi_phases + rcml_number_of_ion_exchanges;
   /* get lines to skip */
   anz = this->rcml_number_of_pqcsteps;

   ein.getline(str,zeilenlaenge);                 /* lies header-Zeile */

   for (i=0; i<this->nodenumber; i++)
   {
      if(this->rateflag[i] > 0)
      {
         /* skip one line, if keyword steps larger than 1 even more lines */
         for(j=0;j<anz;j++)
         {
            for(ii=0;ii<ntot;ii++)
            {
               ein >> dval;
            }
         }
         //		if(1 == 1){
         /*-----------Read the concentration of all master species and pH pe values-------*/
         for (j=0; j<rcml_number_of_master_species+2; j++)
         {
            if(ein >> dval)
               this->val_out[j][i] = dval;
         }

         /*--------------------Read the concentration of all equilibrium phases -------*/
         nj=rcml_number_of_master_species;
         for (j=nj+2; j<nj+2+rcml_number_of_equi_phases; j++)
         {
            if(ein >> dval)
               this->val_out[j][i] = dval;
         }

         /*--------------------Read the concentration of all ion exchangers -------*/
         nk=rcml_number_of_equi_phases;
         for (j=nj+nk+2; j<nj+nk+2+rcml_number_of_ion_exchanges; j++)
         {
            if(ein >> dval)
               this->val_out[j][i] = dval;
         }
         //         }

      }                                           //if rateflag
      else
      {
         for(j=0;j<ntot;j++)
            this->val_out[j][i] =  this->val_in[j][i];
      }
   }

   ok=1;
   return ok;
}


/**************************************************************************
   ROCKFLOW - Funktion: ReadOutputPhreeqcNew

   Aufgabe:
   Liest Ergebnisse der PHREEQC-Berechnungen aus PHREEQC-Ausdgabedatei

  Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *File: Dateiname der PHREEQC-Ausdgabedatei
   R char **val_out: Zeiger für chemische Conzentration
   E char **name: Componentnamen

Ergebnis:
0 bei Fehler oder Ende aufgrund Dateitest, sonst 1

Programmaenderungen:
06/2003     MX/SB         Erste Version
11/2003     SB            Bugfix for large files, read long lines
switched to iostreams
01/2006     SB            ReImplementation, C++ classes, IO, bugfixes
************************************************************************************************/
int REACT::ReadOutputPhreeqcNew(void)
{

   int ok = 0;
   int ntot;
   int index, j, ii, zeilenlaenge=10000, anz;
   char str[4000];
   double dval, dval1;
   string speciesname;
   CRFProcess* m_pcs = NULL;
   int n1, n2, n3, n4, dix=0;
   CTimeDiscretization *m_tim = NULL;

   // Get time step number
   if(time_vector.size()>0)
   {
      m_tim = time_vector[0];
      if(m_tim->step_current == 0) dix = -1;
   }

   ifstream ein (this->results_file_name.data(),ios::in);
   if (!ein)
   {
      cout << "The selected output file doesn't exist!!!" << endl;
      return 0;
   }
   n1 = this->rcml_number_of_master_species;
   n2 = this->rcml_number_of_equi_phases;
   n3 = this->rcml_number_of_ion_exchanges;
   n4 = this->rcml_number_of_gas_species;

   /* get total number of species in PHREEQC output file */
   ntot = rcml_number_of_master_species + 3 + rcml_number_of_equi_phases + rcml_number_of_ion_exchanges + rcml_number_of_gas_species;
   /* get lines to skip */
   anz = this->rcml_number_of_pqcsteps;

   ein.getline(str,zeilenlaenge);                 /* lies header-Zeile */

   for (index=0; index<this->nodenumber; index++)
   {
      if(this->rateflag[index] > 0)
      {
         /* skip one line, if keyword steps larger than 1 even more lines */
         for(j=0;j<anz;j++)
         {
            for(ii=0;ii<ntot;ii++)
            {
               ein >> dval;
            }
         }
         //		if(1 == 1){
         /*-----------Read the concentration of all master species and pH pe values-------*/
         for (j=0; j<n1; j++)
         {
            if(ein >> dval)
            {
               //					this->val_out[j][i] = dval;
               //                    speciesname = pqc_names[j];
               //                    m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
               //                	idx = m_pcs->GetNodeValueIndex(speciesname)+1;
               //        			m_pcs->SetNodeValue(index,idx,dval);
               pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,dval);
               //					if(index <2) cout << " Read aqu. for " << pqc_names[j] << " " << dval << endl;
            }
         }
         /* Read pH and pe */
         if(ein >> dval)                          // read pH
         {
            j = n1;
            pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,dval);
            //				if(index <2) cout << " Read for pH: " << dval << ", ";
         }
         if(ein >> dval)                          // read H+
         {
            j++;
            if(this->gamma_Hplus > 0)
            {
               m_pcs = pcs_vector[pqc_process[j]];
               //				if(index<2) cout << " H+: " <<  dval << ", ";
               pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,dval);
            }
         }
         if(ein >> dval)                          // read pe
         {
            j++;
            m_pcs = pcs_vector[pqc_process[j]];
            //				if(index <2) cout << " pe: " <<  dval << endl;
            pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,dval);
         }
         /*--------------------Read the concentration of all equilibrium phases -------*/
         for (j=n1+3; j<n1+3+n2; j++)
         {
            if(ein >> dval)
            {
               /*
                               speciesname = pqc_names[j];
                               m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
                                 idx = m_pcs->GetNodeValueIndex(speciesname)+1;
                              m_pcs->SetNodeValue(index,idx,dval);
               */
               pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,dval);
               //				if(index <2)  cout << " Read equi. for " << pqc_names[j] << " " << dval << endl;
            }
         }
         /*--------------------Read the concentration of all ion exchangers -------*/
         for (j=n1+3+n2; j<n1+3+n2+n3; j++)
         {
            if(ein >> dval)
            {
               /*                speciesname = pqc_names[j];
                               m_pcs = PCSGet("MASS_TRANSPORT",speciesname);
                                 idx = m_pcs->GetNodeValueIndex(speciesname)+1;
                              m_pcs->SetNodeValue(index,idx,dval);
               */
               pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,dval);
               //                cout << " Read ex. for " << pqc_names[j] << " " << dval << endl;
            }

         }

         /*--------------------Read the concentration of all gas phase species -------*/
         for (j=n1+3+n2+n3; j<n1+3+n2+n3+n4; j++)
         {
            if(ein >> dval)
            {
               pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,dval);
               if(index <2) cout << " Read gas phase for " << pqc_names[j] << " " << dval << endl;
            }
         }

      }                                           //if rateflag
      // Determine new gamma_Hplus
      if(this->gamma_Hplus > 0)
      {
         // Calculate new gamma_Hplus
                                                  // molality H+
         dval =  pcs_vector[pqc_process[n1+1]]->GetNodeValue(index,pqc_index[n1+1]+dix);
                                                  //pH
         dval1 = pcs_vector[pqc_process[n1]]->GetNodeValue(index,pqc_index[n1]+dix);
         dval1 = pow(10.0,-dval1);                // activity H+ from pH
         this->gamma_Hplus = dval1/dval;
         //        cout << " New gamma_Hplus: " << gamma_Hplus << endl;
      }

   }                                              // end for(index...

   ok=1;
   ein.close();
   return ok;
}


/**************************************************************************
   ROCKFLOW - Funktion: REACT::TestPHREEQC

   Aufgabe:
   Testet, ob ein externes PHREEQC- Eingabefile da ist
   allererste Version

   Programmaenderungen:
   09/2003     SB         Erste Version
   01/2006     SB         New File handling
**************************************************************************/
void REACT::TestPHREEQC(string file_base_name)
{

   file_name_pqc = file_base_name + CHEM_REACTION_EXTENSION;
   ifstream pqc_file (file_name_pqc.data(),ios::in);
   if (pqc_file.good())
      flag_pqc = true;
   pqc_file.close();
}


/**************************************************************************
   ROCKFLOW - Funktion: CalculateREactionRateFlag

   Aufgabe:
   Sets all reaction rates to Zero

   Programmaenderungen:
   03/2004     SB         Erste Version
**************************************************************************/
void REACT::CalculateReactionRateFlag(void)
{

   long i, ni, comp, np, j;
   double rate, schwellwert;
   int level=1;
   static int counti;
   int *help;
   double *helprates;
   // Knoten * node;
   // long *neighbor_nodes=NULL;
   // int anz_neighbor_nodes;

   /* Wann soll gerechnet werden: rel. Genauigkeit der Konzentration = deltaC/Cin /Zeitschrittlänge */
   schwellwert = 1.0e-12;
   if(dt>0.0)  schwellwert = 1.0e-4/dt;

   ni = this->nodenumber;
   np = this->number_of_comp;

   /* Calculate chemmistry all rc->countsteps timesteps as well as always in the first and second timestep */
   counti++;
   // test = (aktueller_zeitschritt+1) % rc->count7steps;
   if(((aktueller_zeitschritt+1) % this->countsteps) == 0) counti = 0;
   if(aktueller_zeitschritt < 2) counti = 0;
   if(counti == 0)
   {
      for(i=0;i<ni;i++)
         this->rateflag[i] = 1;
   }
   else
   {
      //determine for each node separately
      DisplayMsgLn(" sum of rates ");

      helprates = (double*) Malloc(np*sizeof(double));
      for(j=0; j<np; j++) helprates[j]=0.0;

      /* Go through all nodes */
      for(i=0;i<ni;i++)
      {
         rate = 0.0;
         for(comp=0; comp<np; comp++)
            rate = rate + fabs(this->rate[comp][i]*this->val_in[comp][i]);
         DisplayLong(i); DisplayMsg(": "); DisplayDouble(rate,0,0); DisplayMsgLn(" ");

         if(rate > schwellwert)
            this->rateflag[i] = 1;
         else
            this->rateflag[i] = 0;
      }                                           //end for(i<ni
      DisplayMsgLn("  ");

      help = (int*) Malloc(ni*sizeof(int));
      for(i=0;i<ni;i++)
         help[i] = this->rateflag[i];

      for(i=0;i<ni;i++)
      {
         if(this->rateflag[i] > 0)
         {
            //get sourrounding elements and nodes
            /* Aussuchen der Nachbarknoten */
            SetNeighborNodesActive (i,level,help);
            //		  DisplayMsgLn(" ");
         }
      }

      //set nodes active finally
      j=0;
      for(i=0;i<ni;i++)
      {
         DisplayMsg(" i: "); DisplayLong(i);
         DisplayMsg(",   this->rf[i]: "); DisplayLong(this->rateflag[i]);
         DisplayMsg(",   this->rate[i]: "); DisplayDouble(this->rate[0][i],0,0);
         DisplayMsg(",   help[i]: ");DisplayLong(help[i]); DisplayMsgLn(" ");
         if(help[i] > 0)
         {
            this->rateflag[i] = 1;
            j++;
         }
      }
      DisplayMsg("Total number of nodes, at which chemistry is calculated: ");DisplayLong(j); DisplayMsgLn(" ");

      //give back storage
      help = (int*) Free(help);
      helprates = (double*) Free(helprates);
   }

}


/**************************************************************************
   ROCKFLOW - Funktion: REACT::GetTransportResults

   Aufgabe:
   Inserts reaction rate of component at index in field RECT->rate

   Programmaenderungen:
   03/2004     SB         Erste Version
   10/2010 TF restructured method a little bit, changed access to process type
**************************************************************************/
void REACT::GetTransportResults(void)
{
   CRFProcess *pcs = NULL;
   int timelevel = 1;                             // concentrations are in new timelevel
   //	int phase = 0; //single phase so far
   size_t np(pcs_vector.size());

   for (size_t j = 0; j < np; j++)                // for all processes
   {
      pcs = pcs_vector[j];
      //		if (pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0) { // if it is a mass transport process
      if (pcs->getProcessType() == MASS_TRANSPORT)
      {
         int comp = pcs->pcs_component_number;    // get component number
         for (long i = 0; i < nodenumber; i++)
         {
            // get concentration values
            val_in[comp][i] = pcs->GetNodeValue(i,
               pcs->GetNodeValueIndex(
               pcs->pcs_primary_function_name[0])
               + timelevel);
                                                  //MX, 01.2005
            if ((val_in[comp][i] < 0.0) && (strcmp(name[comp], "pe") != 0))
            {
               if (abs(val_in[comp][i]) > MKleinsteZahl)
               {
                  DisplayMsg(" Neg. conc for component ");
                  DisplayLong((long) comp);
                  DisplayMsg(" at node ");
                  DisplayLong((long) i);
                  DisplayMsg("; conc = ");
                  DisplayDouble(val_in[comp][i], 0, 0);
                  DisplayMsgLn(" ");
               }
               val_in[comp][i] = 0.0 * val_in[comp][i];
            }
         }
      }
   }

   if (heatflag > 0)
   {
      name[np - 2] = "temp";                      //MX CMCD
      pcs = PCSGet("HEAT_TRANSPORT");
      long index = pcs->GetNodeValueIndex("TEMPERATURE1");
      for (long i = 0; i < this->nodenumber; i++)
                                                  //OK PCSGetNODTemperature1L(i)//MX CMCD -2, Liquid Flow, Heat Transport
         this->val_in[np - 2][i] = pcs->GetNodeValue(i, index);
   }
}


                                                  //ToDo - reactivate with new structure
void REACT::SetNeighborNodesActive(long startnode, long level, int* help)
{
   startnode = startnode;
   level = level;
   help = help;
   /*OK411
   //  long *kanten=NULL;
     int anz_n;
   //  long nd1[2];
     int j=0,k;
     long *knoten;
   //  long* neighbor_nodes =NULL;
     long *elems1d,*elems2d,*elems3d;
     int num_elems1d,num_elems2d,num_elems3d,num_elems;

    // DisplayMsg(" RSNNA: startnode: "); DisplayLong(startnode); DisplayMsg(", level: "); DisplayLong(level); DisplayMsgLn("");
   // Ende rekursiv
   if(level == 0) {
   help[startnode] = 1;
   //	  DisplayMsg(" Knoten, Level 0 "); DisplayLong(startnode); DisplayMsgLn(" ");
   }
   else {

   elems1d = GetNode1DElems(NodeNumber[startnode],&num_elems1d); //SB: liefert zu viele Elemente zurück
   if(num_elems1d > 0) num_elems1d =2;
   elems2d = GetNode2DElems(NodeNumber[startnode],&num_elems2d);
   if(num_elems2d > 4) num_elems2d =4;
   elems3d = GetNode3DElems(NodeNumber[startnode],&num_elems3d);
   if(num_elems3d > 8) num_elems3d =8;
   num_elems=num_elems1d+num_elems2d+num_elems3d;

   // Alle benachbarten Knoten holen
   for (j=0; j<num_elems1d; j++) {
   knoten = ElGetElementNodes(elems1d[j]);
   anz_n = ElGetElementNodesNumber(elems1d[j]);
   for(k=0;k<anz_n;k++){
   if(knoten[k]<0) {
   DisplayMsgLn(" Error");
   }
   SetNeighborNodesActive(knoten[k],level-1,help);
   }
   }
   for (j=0; j<num_elems2d; j++) {
   knoten = ElGetElementNodes(elems2d[j]);
   anz_n = ElGetElementNodesNumber(elems2d[j]);
   for(k=0;k<anz_n;k++){
   SetNeighborNodesActive(knoten[k],level-1,help);
   }
   }
   for (j=0; j<num_elems3d; j++) {
   knoten = ElGetElementNodes(elems3d[j]);
   anz_n = ElGetElementNodesNumber(elems3d[j]);
   for(k=0;k<anz_n;k++){
   SetNeighborNodesActive(knoten[k],level-1,help);
   }
   }

   } //endifs
   */
}


void REACT::ResetpHpe(void)
{
   int i, j;

   for (j=0; j<rcml_number_of_master_species+2; j++)
   {
      if(strcmp(this->name[j], "pe")==0)
         for(i=0;i<this->nodenumber;i++) this->val_out[j][i] = this->val_in[j][i];
      if(strcmp(this->name[j], "pH")==0)
         for(i=0;i<this->nodenumber;i++) this->val_out[j][i] = this->val_in[j][i];
   }
}


double MATCalcIonicStrengthNew(long index)
{
   double ion_strength=0.0;
   index = index;
   /*OK411
     int component, j, n=0, nn, timelevel, z_i;
     static long *element_nodes;
   //  long component;
     double conc[100];
   //  if (rcml)
     REACT *rc = NULL; //SB
     rc->GetREACT();
      n=rc->rcml_number_of_master_species;

   //  n = get_rcml_number_of_master_species(rcml);
   for ( component=0; component<n; component++) {
   conc[component] =0.0;
   if (ElGetElementActiveState(index)){
   nn = ElGetElementNodesNumber(index);
   element_nodes = ElGetElementNodes(index);
   for (j=0;j<nn;j++) {
   timelevel=1;
   conc[component] += PCSGetNODConcentration(element_nodes[j],component,timelevel);
   }
   conc[component] /= (double)nn;
   element_nodes = NULL;
   }

   // calculate the solid phase: volume =v_mi*Ci
   //      conc[i]=CalcElementMeanConcentration (index, i, timelevel, ergebnis);
   CompProperties *m_cp = cp_vec[component];
   z_i = m_cp->valence;

   ion_strength += 0.5 * conc[component]*z_i*z_i;

   }
   */
   return ion_strength;
}



void RCRead(std::string filename)
{
   REACT *rc = new REACT();                       //SB
   rc->TestPHREEQC(filename);                     // Test if *.pqc file is present
   if (rc->flag_pqc)                              //MX
      REACT_vec.push_back(rc);
   else
      delete rc;
   rc =NULL;
}


/**************************************************************************
   ROCKFLOW - Funktion: ExecuteReactionsPHREEQC

   Aufgabe:
   Berechnet chemische Reaktionen zwischen den einzelnen Komponenten
   allererste VErsion

   Programmaenderungen:
   06/2006     MX         Erste Version

**************************************************************************/
void REACT::ExecuteReactionsPHREEQC0(void)
{

   long i,  ok = 0;
   FILE *indatei, *fphinp, *fsel_out=NULL;
   char fsout[80];

   //MDL:
   DisplayMsgLn("ExecuteReactionsPHREEQC0:");
   // DisplayMsgLn("ExecuteReactionsPHREEQC:");

   if (aktuelle_zeit>0)
      GetTransportResults2Element();

   /* Perform reaction step */
   /* --------------------------------------------------------------------------*/
   if(flag_pqc)
   {
      indatei=fopen(crdat,"r");
      if (indatei != NULL)
      {
         ok= this->ReadReactionModel(indatei);
      }
      else
      {
         printf("ERROR: Can not open file *.pqc");
         exit(1);
      }

      /* Read the input file (*.pqc) and set up the input file for PHREEQC ("phinp.dat")*/
      /* set input file name */
      strcpy(fsout, this->outfile);

      fphinp = fopen("phinp.dat", "w");
      if ((fphinp) && ok)
      {
         for(i=0;i<this->elenumber;i++)
         {
            if(this->rateflag[i] > 0)
            {
               rewind(indatei);
               ok = ReadInputPhreeqc(i,indatei, fphinp);
            }
         }
         fclose(indatei);
         fclose(fphinp);
      }
      else
      {
         DisplayMsgLn("The file phinput.dat could not be opened !");
         exit(1);
      }

      /* Extern Program call to PHREEQC */
      if(ok) ok = Call_Phreeqc();
      if(ok ==0) exit(1);

      /* Set up the output values for rockflow after Phreeqc reaction*/
      fsel_out=fopen(fsout,"r");
      if ((ok) && !fsel_out)
      {
         DisplayMsgLn("The selected output file doesn't exist!!!");
         exit(1);
      }
      else if(ok)
      {
         ok = ReadOutputPhreeqc(fsout);
         if(!ok) DisplayMsgLn(" Error in call to PHREEQC !!!");
         fclose(fsel_out);
      }
   }                                              /* if flag */

   /* Calculate Rates */
   CalculateReactionRates();

   // SetConcentrationResults();
   SetConcentrationResultsEle();
}                                                 /* End of ExecuteReactionsPHREEQC */


/**************************************************************************
   ROCKFLOW - Funktion: REACT::GetTransportResults2Element

   Aufgabe:
   Get node results after transport and interpolate to elements (center)

   Programmaenderungen:
   06/2006     MX         Erste Version
   10/2010 TF restructured method a little bit, changed access to process type
**************************************************************************/
void REACT::GetTransportResults2Element()
{
   CRFProcess *pcs = NULL;
   int timelevel = 1;                             // concentrations are in new timelevel
   //	int phase = 0; //single phase so far

   size_t np (pcs_vector.size());
   for (size_t j = 0; j < np; j++)                // for all processes
   {
      pcs = pcs_vector[j];
      //		if (pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0) { // if it is a mass transport process
                                                  // if it is a mass transport process
      if (pcs->getProcessType () == MASS_TRANSPORT)
      {
         int comp = pcs->pcs_component_number;    // get component number
         std::string cname (pcs->pcs_primary_function_name[0]);
         if (CPGetMobil(pcs->GetProcessComponentNumber()) > 0)
            pcs->InterpolateTempGP(pcs, cname);   //line Interpolate to elements for ions

         for (long i = 0; i < elenumber; i++)
         {
            // get concentration values
            val_in[comp][i] = pcs->GetElementValue(i,
               pcs->GetElementValueIndex(cname) + timelevel);
                                                  //MX, 01.2005
            if ((val_in[comp][i] < 0.0) && (strcmp(name[comp], "pe") != 0))
            {
               if (abs(val_in[comp][i]) > MKleinsteZahl)
               {
                  DisplayMsg(" Neg. conc for component ");
                  DisplayLong((long) comp);
                  DisplayMsg(" at node ");
                  DisplayLong((long) i);
                  DisplayMsg("; conc = ");
                  DisplayDouble(val_in[comp][i], 0, 0);
                  DisplayMsgLn(" ");
               }
               val_in[comp][i] = 0.0 * val_in[comp][i];
            }
         }
      }
   }

   // Get Temperature values for elements
   if (this->heatflag > 0)
   {
      this->name[np - 2] = "temp";                //MX CMCD
      pcs = PCSGet("HEAT_TRANSPORT");
      int idxT = pcs->GetElementValueIndex("TEMPERATURE1") + 1;
      pcs->InterpolateTempGP(pcs, "TEMPERATURE1");
      for (long i = 0; i < this->elenumber; i++)
                                                  //OK PCSGetNODTemperature1L(i)//MX CMCD -2, Liquid Flow, Heat Transport
         this->val_in[np - 2][i] = pcs->GetElementValue(i, idxT);

   }
}


/**************************************************************************
PCSLib-Method:
        SB Initialization of REACT structure for rate exchange between MTM2 and Reactions
07/2007 OK Encapsulation
**************************************************************************/
void REACTInit()
{
   REACT *rc = NULL;                              //SB
   //  rc->TestPHREEQC(); // Test if *.pqc file is present
   rc = rc->GetREACT();
   if(rc)                                         //OK
   {
      if(rc->flag_pqc)
      {
         if(cp_vec.size()>0)
         {                                        //OK
#ifdef REACTION_ELEMENT
            rc->CreateREACT();                    //SB
            rc->InitREACT0();
            rc->ExecuteReactionsPHREEQC0();
            REACT_vec.clear();
            REACT_vec.push_back(rc);
#else
            //--------------------------------------------------
            // HB, for the GEM chemical reaction engine 05.2007
#ifdef REAC_GEM
            REACT_GEM *p_REACT_GEM = NULL;
            p_REACT_GEM->REACT_GEM();
#else
            //--------------------------------------------------
            rc->CreateREACT();                    //SB
            rc->InitREACT();
#ifdef LIBPHREEQC
            cout << "MDL Calling ExecuteReactionsPHREEQCNewLib" << endl;
            rc->ExecuteReactionsPHREEQCNewLib();
#else
            rc->ExecuteReactionsPHREEQCNew();
#endif                                //LIBPHREEQC
            REACT_vec.clear();
            REACT_vec.push_back(rc);
#endif
#endif
         }
      }
      //  delete rc;
   }
#ifdef CHEMAPP
   CEqlink *eq=NULL;
   eq = eq->GetREACTION();
   if(cp_vec.size()>0  && eq)                     //MX
   {
      eq->TestCHEMAPPParameterFile(pcs_vector[0]->file_name_base);
      if (eq->flag_chemapp)
      {
         eq->callCHEMAPP(pcs_vector[0]->file_name_base);
      }
   }
#endif
   //  delete rc;
}


/***********************************************************************************************
   ROCKFLOW - Funktion: CheckNoReactionNodes

   Aufgabe:
   Checks, if in KinReactData (from kinetic reactions input file) geometries
   are defined, where no reactions are to be calculated. These are then
   transferred node wise to vector rateflag

   Ergebnis:
   0 bei Fehler oder Ende aufgrund Dateitest, sonst 1

Programmaenderungen:
12/2007     SB	         Erste Version
************************************************************************************************/
int REACT::CheckNoReactionNodes(void)
{

   int ok = 0;
   long l;

   cout << " CheckNoReactionNodes " << endl;

   if(aktueller_zeitschritt < 2)                  //do not in the very first calculation before first time step and in the first time step
   {
      this->check_no_reaction_nodes = false;
   }
   else
   {
      CKinReactData *m_krd = NULL;
      if(KinReactData_vector.size()>0)
         m_krd = KinReactData_vector[0];
      if(m_krd == NULL)
      {
         // no KinReactData specified in *.krc file
         cout << "No CKinReactData available in CheckNoReactionNodes" << endl;
      }
      else
      {
         if(m_krd->is_a_CCBC.size() > 0)          // reaction nodes specified in krc input file
         {
            CFEMesh* m_msh = fem_msh_vector[0];   //SB: ToDo hart gesetzt
            if(m_msh == NULL) {cout << "No mesh in CheckNoReactionNodes" << endl; exit(1);}
            // Initialize vector is_a_CCBC
                                                  // node 1 needed for phreeqc-input
            for(l=1; l< (long)m_msh->nod_vector.size();l++)
            {
                                                  // rateflag == 0 means no reactions calculated
               if(m_krd->is_a_CCBC[l] == true) this->rateflag[l] = 0;
               // cout << " Node " << l << " is " << this->rateflag[l] << endl;
            }
         }
      }

      this->check_no_reaction_nodes = true;
   }

   ok=1;
   return ok;
}


// MDL: here the new functions begin
#ifdef LIBPHREEQC

/**************************************************************************
   ROCKFLOW - Funktion: WriteInputPhreeqcLib

   Aufgabe:
   Read input data *.pqc and form the buffer (stringstream) for libphreeqc

   Formal parameters:
   In:   index = number of simulation (==node)
   Out:  out_buff = stringstream with the simulation's phreeqc input
         nl = number of line of the current input

Ergebnis:
0 bei Fehler oder Ende aufgrund Dateitest, sonst 1

Programmaenderungen:
04/2009     MDL        First Version
**************************************************************************/
int REACT::WriteInputPhreeqcLib(long index, stringstream* out_buff, int* nl)
{
   char line[MAX_ZEILE];
   std::stringstream in;
   string name, line_string, speciesname, dummy;
   CRFProcess* m_pcs = NULL;
   int i, ii, idx, n1, n2, n3, n4, count=-1;
   int nline;
   double dval, dval1;
   double z, h, dens, press, partial_press, volume, temp=-1.0, mm;

   cout.flush();
   //  cout << "WriteInputPhreeqcLib starting ..." << endl;
   ifstream pqc_infile (this->file_name_pqc.data(),ios::in);
   pqc_infile.seekg(0L,ios::beg);

   // precision output file
   out_buff->setf(ios::scientific,ios::floatfield);
   out_buff->precision(12);

   // reset the (partial) counter of lines
   nline=0;
   //  cout << "Write: *nline =" << *nline << "; nline++ ="<< nline++ << endl;
   // linewise read
   while(!pqc_infile.eof())
   {
      pqc_infile.getline(line,MAX_ZEILE);
      line_string = line;
      if(line_string.find("#STOP")!=string::npos)
         break;
      //-------------------------------------------------------------------------------------------------------------
      /* loop for Keyword SOLUTION */
                                                  // keyword found
      if(line_string.find("SOLUTION")!=string::npos)
      {
         *out_buff << "SOLUTION " << index+1 << endl;
         nline++;
         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile.getline(line,MAX_ZEILE);
            line_string = line;
            if (line_string.find("# comp")!=string::npos)
            {
               if (line_string.find("pH") == string::npos && line_string.find("pe") == string::npos)
               {
                  // Component found; write name and concentration of component
                  count++;
                  speciesname = pqc_names[count];
                  dval = pcs_vector[pqc_process[count]]->GetNodeValue(index,pqc_index[count]);

                  if(speciesname.compare("pe"))   // if this is not pe
                     if(dval < 1.0e-19) dval = 0.0;
                  *out_buff << speciesname << " " << dval << endl;
                  nline++;
               }
            }
            else if (line_string.find("# temp")!=string::npos)
            {
               // check if heat transport process is calculated in GeoSys
               if(this->rcml_heat_flag > 0)
               {
                  m_pcs = PCSGet("HEAT_TRANSPORT");
                  idx = m_pcs->GetNodeValueIndex("TEMPERATURE1");
                  dval = m_pcs->GetNodeValue(index,idx);
                  if (dval<273.0) dval += 273.15; //change from °C to Kelvin if necessary
                  dval -=273.15;                  // Input to PHREEQC is in °C
                  *out_buff << "temp " << dval << endl ;
                  nline++;
                  temp = dval;                    // save for gas phase input
               }
            }
            else
            {                                     // Write units and temperature in the standard case
               if (line_string.find("pH") == string::npos && line_string.find("pe") == string::npos && line_string.find("#ende") == string::npos)
               {
                  *out_buff << line_string << endl;
                  nline++;
               }
            }
         }                                        // end while

         // special treat pH, and pe
         n1 = this->rcml_number_of_master_species;
         count++;
         if(count != n1)
            cout << "Error in index of pqc_vectors !" << endl;
         dval = pcs_vector[pqc_process[count]]->GetNodeValue(index,pqc_index[count]);
         count++;
         if(this->gamma_Hplus > 0)                // pH and H+ in GeoSys species, calculate pH from H+
         {
            dval1 = pcs_vector[pqc_process[n1+1]]->GetNodeValue(index,pqc_index[n1+1]);
            dval = -log10(dval1*gamma_Hplus);
         }
         if(this->rcml_pH_charge > 0)
         {
            *out_buff << "pH " << dval << " charge " << endl;
            nline++;
         }
         else
         {
            *out_buff << "pH " << dval << endl;
            nline++;
         }
         // pe
         count++;
         dval = pcs_vector[pqc_process[n1+2]]->GetNodeValue(index,pqc_index[n1+2]);
         *out_buff << "pe " << dval << endl;
         nline++;
         if (line_string.find("#ende")!=0)
         {
            *out_buff << line_string << endl;
            nline++;
         }
      }                                           // end SOLUTION

      // Keyword EQUILIBRIUM PHASES
      if(line_string.find("EQUILIBRIUM_PHASES")!=string::npos)
      {                                           // keyword found
         *out_buff << "PURE " << index+1 << endl;
         nline++;
         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile.getline(line,MAX_ZEILE);
            line_string = line;
            if (line_string.find("# comp")!=string::npos)
            {
               count++;
               speciesname = pqc_names[count];
               dval = pcs_vector[pqc_process[count]]->GetNodeValue(index,pqc_index[count]);
               if(dval < 1.0e-19)
                  dval = 0.0;
               *out_buff << speciesname << " 0.0 " << dval << endl;
               nline++;
            }
            else
            if (line_string.find("#ende")!=0)
               *out_buff << line_string << endl;
         }
      }                                           // end EQUILIBRIUM PHASES

      // Keyword EXCHANGE
      if(line_string.find("EXCHANGE")!=string::npos)
      {                                           // keyword found
         *out_buff << "EXCHANGE " <<  index+1 << endl;
         nline++;
         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile.getline(line,MAX_ZEILE);
            line_string = line;
            if (line_string.find("# comp") !=string::npos)
            {
               count++;
               speciesname = pqc_names[count];
               dval = pcs_vector[pqc_process[count]]->GetNodeValue(index,pqc_index[count]);
               if(dval < 1.0e-19)
                  dval = 0.0;
               *out_buff << speciesname << "       " << dval << endl;
               nline++;
            }
            else
            if (line_string.find("#ende")!=0)
            {
               *out_buff << line_string << endl;
               nline++;
            }
         }
      }                                           // end EXCHANGE

      // Keyword GAS_PHASE
      if(line_string.find("GAS_PHASE")!=string::npos)
      {                                           // keyword found
         *out_buff << "GAS_PHASE " <<  index+1 << endl;
         nline++;

         // get necessary values for conversion of molar concentrations to partial pressures, and to calculate total pressure and total volume

         // get height of node z
         CFEMesh* m_msh = fem_msh_vector[0];      //SB: ToDo hart gesetzt
         Mesh_Group::CNode* m_nod = NULL;
         m_nod = m_msh->nod_vector[index];
         z = m_msh->nod_vector[index]->Z();

         // get piezometric hight h
         m_pcs = PCSGet("GROUNDWATER_FLOW");
         if(m_pcs == NULL)
            cout << "   Error - no flow process found!" << endl;
         idx = m_pcs->GetNodeValueIndex("HEAD")+1;
         h = m_pcs->GetNodeValue(index,idx);

         // get fluid density
         dens = mfp_vector[0]->Density();

         // calculate pressure in [Pa]
         press = dens * gravity_constant * (h-z);

         // get temperature in [°C]
         if(rcml_heat_flag < 1)
            temp = this->temperature;

         // get molar masses of gas phase
         mm=0.0;                                  // mm is total molar mass of gas phase in [mol]
         ii= rcml_number_of_master_species + 3 + rcml_number_of_ion_exchanges + rcml_number_of_equi_phases;
         for(i=ii;i<ii+rcml_number_of_gas_species;i++)
         {
            speciesname = this->pqc_names[i];
            dval = pcs_vector[pqc_process[i]]->GetNodeValue(index,pqc_index[i]);
            mm += dval;
         }

         //  calculate Volume of gas phase in [mol * Pa * m^3 / K / mol * K / Pa = m^3 ]
         volume = mm * 8.314472 * (273.15 + temp) /press;

         while(line_string.find("#ende")==string::npos)
         {
            pqc_infile.getline(line,MAX_ZEILE);
            line_string = line;
            if (line_string.find("-pressure") !=string::npos)
            {
                                                  // pressure in atmospheres
               *out_buff << "-pressure " << press/101325.0 << endl;
               nline++;
            }
            else if (line_string.find("-volume") !=string::npos)
            {
                                                  // volume in Liters
               *out_buff << "-volume " << volume*1000.0 << endl;
               nline++;
            }
            else if (line_string.find("-temperature") !=string::npos)
            {
                                                  // temperature in °Celsius
               *out_buff << "-temperature " << temp << endl;
               nline++;
            }
            else if (line_string.find("# comp") !=string::npos)
            {
               count++;
               speciesname = pqc_names[count];
               dval = pcs_vector[pqc_process[count]]->GetNodeValue(index,pqc_index[count]);
               if(dval < 1.0e-19)
                  dval = 0.0;
               if (mm > 0.0)
               {
                  partial_press = press * dval/mm;
               }
               else
               {
                  partial_press = 0.0;
               }

               *out_buff << " " << speciesname << " " << partial_press/101325.0 << endl;
               nline++;
            }
            else
            if (line_string.find("#ende")!=0)
            {
               *out_buff << line_string << endl;  //write line unchanged
               nline++;
            }
         }
      }                                           // end GAS_PHASE

      // Keyword SELECTED_OUTPUT
      if(line_string.find("SELECTED_OUTPUT")!=string::npos)
      {                                           // keyword found
         if(index < 1)                            // this block has to appear just in the first solution
         {
            *out_buff << "SELECTED_OUTPUT" << endl;
            nline++;
            while(line_string.find("#ende")==string::npos)
            {
               pqc_infile.getline(line,MAX_ZEILE);
               line_string = line;
               if (line_string.find("-file") !=string::npos)
               {
                  *out_buff << "-file " << this->results_file_name << endl;
                  nline++;
               }
               else
               if (line_string.find("#ende")!=0)
               {
                  *out_buff << line_string << endl;
                  nline++;
               }
            }
         }
      }                                           // end SELECTED OUTPUT

      // Keyword PRINT
      if(line_string.find("PRINT")!=string::npos)
      {                                           // keyword found
         if(index < 1)                            // this block has to appear just in the first solution
         {
            *out_buff << "PRINT" << endl;
            nline++;
            while(line_string.find("#ende")==string::npos)
            {
               pqc_infile.getline(line,MAX_ZEILE);
               line_string = line;
               if (line_string.find("#libprint")!=string::npos)
               {
                  libphreeqc_print="T";
               }
               else
               {
                  if (line_string.find("#ende")!=0)
                  {
                     *out_buff << line_string << endl;
                     nline++;
                  }
               }
            }
         }
      }                                           // end PRINT

      // Keyword USER_PUNCH
      if(line_string.find("USER_PUNCH")!=string::npos)
      {                                           // keyword found
         if(index < 1)
         {
            *out_buff << "USER_PUNCH" << endl;
            nline++;
            // Write Header
            n1 = this->rcml_number_of_master_species;
            n2 = this->rcml_number_of_equi_phases;
            n3 = this->rcml_number_of_ion_exchanges;
            n4 = this->rcml_number_of_gas_species;
            *out_buff << "-head ";
            for(i=0;i<n1; i++)
               *out_buff << " " << pqc_names[i];
            *out_buff << " pH ";
            *out_buff << " H+ ";
            *out_buff << " pe ";
            for(i=n1+3; i<n1+3+n2; i++)
               *out_buff << " " << pqc_names[i];
            for(i=n1+3+n2;i<n1+3+n2+n3; i++)
               *out_buff << " " << pqc_names[i];
            for(i=n1+3+n2+n3; i<n1+3+n2+n3+n4; i++)
               *out_buff << " " << pqc_names[i];
            *out_buff << endl;
            nline++;
            // Write everything in just 1 line
            *out_buff << "10 PUNCH ";
            for(i=0;i<n1;i++)
            {
               if(pqc_names[i].compare("H+")==0)
                                                  // extra treat H+
                  *out_buff << " MOL(\"" << pqc_names[i] << "\"),";
               else
                                                  // without pH and pe here
                  *out_buff << " TOT(\"" << pqc_names[i] << "\"),";
            }

            *out_buff << " -LA(\"H+\"),  MOL(\"H+\"),  -LA(\"e-\")";

            // Write equilibrium phases
            if(n2 > 0)
            {
               for(i=n1+3;i<n1+3+n2;i++)
                  *out_buff << ", EQUI(\"" << pqc_names[i] << "\")";
            }

            // Write ion exchangers
            if(n3 > 0)
            {
               for(i=n1+3+n2;i<n1+3+n2+n3;i++)
                  *out_buff << ", MOL(\"" << pqc_names[i] << "\")";
            }

            // Write gas phase species
            if(n4 > 0)
            {
               for(i=n1+3+n2+n3; i<n1+3+n2+n3+n4; i++)
                  *out_buff << ", GAS(\"" << pqc_names[i] << "\")";
            }

            // MDL: now the endl
            *out_buff << endl;
            nline++;
         }                                        // end if index < 1

         // search for end of USER_PUNCH data block in *.pqc input file
         while(!pqc_infile.eof())
         {
            pqc_infile.getline(line,MAX_ZEILE);
            line_string = line;
            if((line_string.find("#ende")!=string::npos) || (line_string.find("END")!=string::npos))
               break;
         }
      }                                           // end USER_PUNCH

      if(line_string.find("-steps")!=string::npos)
      {                                           // keyword found
         in.str(line_string);
         in >> dummy >> dval >> this->rcml_number_of_pqcsteps >> dummy;
         CTimeDiscretization *m_tim = NULL;
         if(time_vector.size()>0)
            m_tim = time_vector[0];
         else
            cout << "Error in WriteInputPhreeqcLib: no time discretization data !" << endl;
         dval = m_tim->CalcTimeStep();
         *out_buff << "-steps " << dval << " in " << this->rcml_number_of_pqcsteps << " steps" << endl;
         nline++;
      }                                           // end -steps

      if(line_string.find("KNOBS")!=string::npos)
      {
         if(index < 1)
         {
            *out_buff << "KNOBS" << endl;
            nline++;
            while(line_string.find("#ende")==string::npos)
            {
               pqc_infile.getline(line,MAX_ZEILE);
               line_string = line;
               *out_buff << line_string << endl;
               nline++;
            }
         }
      }

   }                                              // end of "while" linewise read
   *out_buff << "END" << endl;
   nline++;
   *nl = nline;
   pqc_infile.close();

   return 1;
}


/**************************************************************************
   ROCKFLOW - Funktion: ExecuteReactionsPHREEQCNewLib

   04/2009     MDL         First Version

**************************************************************************/
void REACT::ExecuteReactionsPHREEQCNewLib(void)
{
   long i, ii,  ok = 0;
   int nl, nline, npunch;

   // File handling - GeoSys input file
   ifstream pqc_file (this->file_name_pqc.data(),ios::in);
   if (!pqc_file.good())
   {
      cout << "! Error in ExecuteReactionsPHREEQCNew: no Input File (*.pqc) found !" << endl;
   }

   // data exchange buffer to libphreeqc, as stringstream
   stringstream out_buff;

   // Set up reaction model
   if((int)this->pqc_names.size()==0)
   {
      ok = this->ReadReactionModelNew(&pqc_file);
      if(!ok)
         cout << "Error setting up reaction model" << endl;
   }

   // Check for nodes without reactions
   if((int)this->check_no_reaction_nodes == false)
   {
      ok = this->CheckNoReactionNodes();
      if(!ok)
         cout << "Error when checking for nodes without reactions" << endl;
   }

   // Read the input file (*.pqc) and set up the input for PHREEQC;
   // Write input data block to PHREEQC for each node;
   // sum number of lines of input.
   cout << endl << "Preparing phreeqc's input...";
   ii = 0;
   nline = 0;
   // Should libphreeqc print to file? Defaults to no (=="F")
   libphreeqc_print="F";

   for(i=0;i<this->nodenumber;i++)
      if(this->rateflag[i] > 0)
   {
      pqc_file.seekg(0L,ios_base::beg);
      ok = WriteInputPhreeqcLib(i, &out_buff, &nl);
      ii++;
      nline = nline+nl;
   }

   if(ok)
      cout << " [OK]" << endl;
   else
      cout << " [ERROR]" << endl;

   // TODO: don't try to run phreeqc if output is not ok

#ifdef MDL_DEBUG
   cout << "MDL_DEBUG: libphreeqc_print = " << libphreeqc_print << endl;
   cout << "MDL_DEBUG: final input is " << nline << " lines, nodes " << ii << "/" << this->nodenumber <<endl;
   cout << endl << "MDL_DEBUG: ****** Current input:" << endl << out_buff.str() << endl << "MDL_DEBUG: ****** end of input" << endl;
#endif

   //  Close *.pqc input file
   pqc_file.close();

   npunch = this->rcml_number_of_master_species + this->rcml_number_of_equi_phases +
      this->rcml_number_of_ion_exchanges+ this->rcml_number_of_gas_species + 3;

   double* out_vec = new double [ npunch * ii ];

   // call to libphreeqc
   std::cout << endl << std::endl;
   ok = Call_PhreeqcLib(ii,npunch, nline, &out_buff, out_vec);

   if(ok)
   {
      cout << "Reading phreeqc's output...";
      ok = ReadOutputPhreeqcNewLib(out_vec);
      if(!ok)
         std::cout << " [ERROR] ExecuteReactionsPHREEQCNewLib: Error reading phreeqc's output" << std::endl;
      else
         std::cout << " [OK]" << std::endl;
   }

   // deallocation
   delete [] out_vec;

   /* Calculate Rates */
   // CalculateReactionRates();

   /* determine where to calculate the chemistry */
   // CalculateReactionRateFlag();

   /* pH and pe constant or variable */
   // ResetpHpe(rc, rcml);

}                                                 /* End of ExecuteReactionsPHREEQCNewLib */


/**************************************************************************
   ROCKFLOW - Funktion: ReadOutputPhreeqcNewLib

   Aufgabe:
   Liest Ergebnisse der PHREEQC-Berechnungen aus PHREEQC-Ausdgabedatei

   04/2009     MDL
************************************************************************************************/
int REACT::ReadOutputPhreeqcNewLib(double* pqc_vec)
{
   int ok = 0;
   int ntot;
   int index, j;
   double dval, dval1;
   CRFProcess* m_pcs = NULL;
   int n1, n2, n3, n4, dix=0;
   CTimeDiscretization *m_tim = NULL;

   // Get time step number
   if(time_vector.size()>0)
   {
      m_tim = time_vector[0];
      if(m_tim->step_current == 0)
         dix = -1;
   }

   // get total number of species in the vector
   n1 = this->rcml_number_of_master_species;
   n2 = this->rcml_number_of_equi_phases;
   n3 = this->rcml_number_of_ion_exchanges;
   n4 = this->rcml_number_of_gas_species;

   ntot = n1 + 3 + n2 + n3 + n4;

   for (index=0; index<this->nodenumber; index++)
   {
      if(this->rateflag[index] > 0)
      {
         // concentrations of master species and pH pe values
         for (j=0; j<n1; j++)
         {
            pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,pqc_vec[index*ntot+j]);
         }

         // pH
         j = n1;
         pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,pqc_vec[index*ntot+j]);

         // H+
         j++;
         if(this->gamma_Hplus > 0)
         {
            m_pcs = pcs_vector[pqc_process[j]];
            pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,pqc_vec[index*ntot+j]);
         }
         // pe
         j++;
         m_pcs = pcs_vector[pqc_process[j]];
         pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,pqc_vec[index*ntot+j]);

         // concentrations of all equilibrium phases
         if (n2>0)
         {
            for (j=n1+3; j<n1+3+n2; j++)
            {
               pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,pqc_vec[index*ntot+j]);
            }
         }

         // concentrations of all ion exchangers
         if (n3 > 0)
         {
            for (j=n1+3+n2; j<n1+3+n2+n3; j++)
            {
               pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,pqc_vec[index*ntot+j]);
            }
         }

         // concentrations of all gas phase species
         if (n4 > 0)
         {
            for (j=n1+3+n2+n3; j<n1+3+n2+n3+n4; j++)
            {
               pcs_vector[pqc_process[j]]->SetNodeValue(index,pqc_index[j]+dix,pqc_vec[index*ntot+j]);
            }
         }
      }                                           //endif rateflag

      // Determine new gamma_Hplus
      if(this->gamma_Hplus > 0)
      {
         // Calculate new gamma_Hplus
                                                  // molality H+
         dval =  pcs_vector[pqc_process[n1+1]]->GetNodeValue(index,pqc_index[n1+1]+dix);
                                                  //pH
         dval1 = pcs_vector[pqc_process[n1]]->GetNodeValue(index,pqc_index[n1]+dix);
         dval1 = pow(10.0,-dval1);                // activity H+ from pH
         this->gamma_Hplus = dval1/dval;
      }

   }                                              // end for (index...

   ok=1;

   return ok;
}


/**************************************************************************
   ROCKFLOW - Funktion: Call_PhreeqcLib

   Ergebnis:
   0 bei Fehler oder Ende aufgrund Dateitest, sonst 1

   Programmaenderungen:
   04/2009     MDL         First Version
 **************************************************************************/

int REACT::Call_PhreeqcLib(int nsim, int npunch, int nline, stringstream* pqc_buffer, double* pqc_out)
{

   int i, out;
   string line;
   size_t length;

#ifdef MDL_DEBUG
   cout <<"MDL_DEB:: Call_PhreeqcLib ->  nsim = "<< nsim << "; npunch = "  << npunch << "; nline = " << nline << endl;
#endif

   // copy the stringstream to a char ** for the c library
   char** buffer = new char* [nline];

   for (i=0; i<nline; i++)
   {
      getline(*pqc_buffer, line);
      buffer[i] = new char[line.size()+1];

      length = line.copy(buffer[i],line.size()+1,0);
      buffer[i][length]='\0';
   }

   // prepare args for Phreeqcmain (must be char **)
   char** libargs = new char*[5];
   for (i=0;i<5;i++)
      libargs[i] = new char[30];

   strcpy(libargs[0],"libphreeqc");
   strcpy(libargs[1],"FIXME");
   strcpy(libargs[2],"phinp.out");
   strcpy(libargs[3],"phreeqc.dat");
   strcpy(libargs[4],libphreeqc_print.c_str());

#ifdef MDL_DEBUG
   cout << "libargs assigned\n";
#endif

   // call to libphreeqc
   // NB: Phreeqcmain returns 1 if everything is fine, 0 if errors
   out = Phreeqcmain( 5, libargs, nsim, npunch,
      nline, buffer, pqc_out);

   // deallocations
   for (i=0;i<nline;i++)
      delete [] buffer[i];

   delete [] buffer;

   for (i=0; i<5;i++)
      delete [] libargs[i];
   delete [] libargs;

   if (out == 1)
      return(1);                                  // ok == 1 for Geosys conventions
   else
   {
      DisplayMsgLn("Warning: libphreeqc doesn't run properly!!! ");
      exit(1);
   }

}
#endif                                            // LIBPHREEQC
