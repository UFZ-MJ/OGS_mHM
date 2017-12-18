/**************************************************************************
FEMLib - Object: TIM
Task:
Programing:
08/2004 OK Implementation
last modified:
**************************************************************************/
// C++ STL
//#include <math.h>
//#include <iostream>
#include <cfloat>
// FEM-Makros
#include "makros.h"
// GeoSys-GeoLib
#include "files0.h"
// GeoSys-FEMLib
#include "rf_tim_new.h"
//#include "rf_pcs.h"
#include "rf_mmp_new.h"
#include "fem_ele_std.h"
#include "mathlib.h"
// kg44 not found #include "elements.h"
#include "rfmat_cp.h"
//WW #include "elements.h" //set functions for stability criteria
// ToDo
double aktuelle_zeit;
size_t aktueller_zeitschritt = 0;
double dt = 0.0;
int rwpt_numsplits = -1;                          //JTARON 2010
//==========================================================================
std::vector<CTimeDiscretization*>time_vector;
/**************************************************************************
FEMLib-Method:
Task: OBJ constructor
Programing:
08/2004 OK Implementation
**************************************************************************/
CTimeDiscretization::CTimeDiscretization(void)
:Write_tim_discrete(false),tim_discrete(NULL)     //YD
{
   step_current = 0;
   time_start = 0.0;
   time_end = 1.0;
   time_type_name = "CONSTANT";                   //OK
   time_control_name = "";                        //kg44
   time_unit = "SECOND";
   max_time_step = 1.e10;                         //YD
   min_time_step = 0;                             //YD
   courant_desired = 0.5;                         //JTARON
   courant_initial = 1.e-6;                       //JTARON
   courant_static = 0;                            //JTARON
   repeat = false;                                //OK/YD
   step_current = 0;                              //WW
   this_stepsize = 0.;                            //WW
   dt_sum = .0;                                   //WW
   relative_error = 1.e-4;                        //26.08.2008. WW
   absolute_error = 1.e-10;                       //26.08.2008. WW
   h_max = 6;                                     //27.08.2008. WW
   h_min = 0.2;                                   //27.08.2008. WW
   hacc = 0.;                                     //27.08.2008. WW
   erracc = 0.;                                   //27.08.2008. WW
   tsize_ctrl_type = -1;                          //27.08.2008. WW
}


/**************************************************************************
FEMLib-Method:
Task: OBJ destructor
Programing:
08/2004 OK Implementation
**************************************************************************/
CTimeDiscretization::~CTimeDiscretization(void)
{
   if(tim_discrete)                               //YD
   {
      tim_discrete->close();
      if(tim_discrete) delete tim_discrete;
      time_step_vector.clear();
      time_adapt_tim_vector.clear();
      time_adapt_coe_vector.clear();
   }
}

std::ios::pos_type GetNextSubKeyword(std::ifstream* file,std::string* line, bool* keyword)
{
   char buffer[MAX_ZEILE];
   std::ios::pos_type position;
   position = file->tellg();
   *keyword = false;
   std::string line_complete;
   int i,j;
   // Look for next subkeyword
   while(!(line_complete.find("$")!=std::string::npos)&&(!file->eof()))
   {
      file->getline(buffer,MAX_ZEILE);
      line_complete = buffer;
      if(line_complete.find("#")!=std::string::npos)
      {
         *keyword = true;
         return position;
      }
                                                  //Anf�ngliche Leerzeichen �berlesen, i=Position des ersten Nichtleerzeichens im string
      i = (int) line_complete.find_first_not_of(" ",0);
      j = (int) line_complete.find(";",i);        //Nach Kommentarzeichen ; suchen. j = Position des Kommentarzeichens, j=-1 wenn es keines gibt.
      if(j<0)
         j = (int)line_complete.length();
      //if(j!=i) break;						 //Wenn das erste nicht-leerzeichen ein Kommentarzeichen ist, zeile �berlesen. Sonst ist das eine Datenzeile
      if(i!=-1)
         *line = line_complete.substr(i,j-i);     //Ab erstem nicht-Leerzeichen bis Kommentarzeichen rauskopieren in neuen substring, falls Zeile nicht leer ist
   }
   return position;
}


/**************************************************************************
FEMLib-Method:
Task: OBJ read function
Programing:
08/2004 OK Implementation
11/2004 OK string streaming by SB for lines
10/2005 YD Time Controls
08/2008 WW General classic time step size control (PI control)
**************************************************************************/
std::ios::pos_type CTimeDiscretization::Read(std::ifstream *tim_file)
{
   std::string sub_line;
   std::string line_string;
   std::string delimiter(" ");
   bool new_keyword = false;
   std::string hash("#");
   std::ios::pos_type position;
   std::string sub_string;
   bool new_subkeyword = false;
   std::string dollar("$");
   int no_time_steps = 0;
   double time_step_length;
   std::ios::pos_type position_subkeyword;
   std::stringstream line;
   std::string line_complete;
   int iter_times;                                //YD
   double multiply_coef;                          //YD
   int i;
   CRFProcess* m_pcs = NULL;
   //    m_pcs = PCSGet("RICHARDS_FLOW");
   m_pcs = PCSGet("GROUNDWATER_FLOW");            //kg44 changed default

   //========================================================================
   // Schleife ueber alle Phasen bzw. Komponenten
   while(!new_keyword)
   {
      if(new_subkeyword)
         tim_file->seekg(position,std::ios::beg);
      new_subkeyword = false;
      position = GetNextSubKeyword(tim_file,&line_string,&new_keyword);
      if(new_keyword)
         return position;
      /*
          position = tim_file->tellg();
          if(new_subkeyword)
            tim_file->seekg(position_subkeyword,ios::beg);
          new_subkeyword = false;
          tim_file->getline(buffer,MAX_ZEILE);
          line_string = buffer;
         if(line_string.size()<1) // empty line
            continue;
          if(Keyword(line_string))
            return position;
      */
      //....................................................................

                                                  // subkeyword found
      if(line_string.find("$PCS_TYPE")!=std::string::npos)
      {
         line.str(GetLineFromFile1(tim_file));
         line >> pcs_type_name;
         line.clear();
         m_pcs = PCSGet(pcs_type_name);           // kg44 inserted to overwrite default Richards_flow
         // this works only of pcs_type is read before adaption
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$TIME_START")!=std::string::npos)
      {
         line.str(GetLineFromFile1(tim_file));
         line >> time_start;
         line.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$TIME_END")!=std::string::npos)
      {
         line.str(GetLineFromFile1(tim_file));
         line >> time_end;
         line.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$TIME_UNIT")!=std::string::npos)
      {
         *tim_file>>time_unit>>std::ws;           //WW unit of time
         continue;
      }
      //....................................................................
      /* //WW
      if(line_string.find("$TIME_FIXED_POINTS")!=string::npos) { // subkeyword found
       int no_fixed_points;
        double fixed_point;
       line.str(GetLineFromFile1(tim_file));
        line >> no_fixed_points;
       line.clear();
       for(i=0;i<no_fixed_points;i++) {
          line.str(GetLineFromFile1(tim_file));
          line >> fixed_point;
        fixed_point_vector.push_back(fixed_point);
      line.clear();
      }
      continue;
      }
      */
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$TIME_STEPS")!=std::string::npos)
      {
         while((!new_keyword)||(!new_subkeyword)||(!tim_file->eof()))
         {
            position = tim_file->tellg();
            line_string = GetLineFromFile1(tim_file);
            if(line_string.find("#")!=std::string::npos)
            {
               return position;
            }
            if(line_string.find("$")!=std::string::npos)
            {
               new_subkeyword = true;
               break;
            }
            line.str(line_string);
            line >> no_time_steps;
            line >> time_step_length;
            for(i=0;i<no_time_steps;i++)
               time_step_vector.push_back(time_step_length);
            line.clear();
         }
      }
                                                  // subkeyword found
      if(line_string.find("$TIME_SPLITS")!=std::string::npos)
      {
         line.str(GetLineFromFile1(tim_file));
         line >> rwpt_numsplits;
         line.clear();
         continue;
      }
                                                  // 25.08.2008. WW
      if(line_string.find("$CRITICAL_TIME")!=std::string::npos)
      {
         while((!new_keyword)||(!new_subkeyword)||(!tim_file->eof()))
         {
            position = tim_file->tellg();
            line_string = GetLineFromFile1(tim_file);
            if(line_string.find("#")!=std::string::npos)
            {
               return position;
            }
            if(line_string.find("$")!=std::string::npos)
            {
               new_subkeyword = true;
               break;
            }
            line.str(line_string);
            double  crtime;
            line >> crtime;
            critical_time.push_back(crtime);
            line.clear();
         }
      }
                                                  // subkeyword found
      if(line_string.find("$TIME_CONTROL")!=std::string::npos)
      {
         while((!new_keyword)||(!new_subkeyword)||(!tim_file->eof()))
         {
            position = tim_file->tellg();
            line_string = GetLineFromFile1(tim_file);

            if(line_string.find("#")!=std::string::npos)
            {
               return position;
            }
            if(line_string.find("$")!=std::string::npos)
            {
               new_subkeyword = true;
               break;
            }
            line.str(line_string);
            line >> time_control_name;
            line.clear();

            // 26.08.2008. WW
            if(time_control_name=="PI_AUTO_STEP_SIZE")
            {
               line.str(GetLineFromFile1(tim_file));
               line >> tsize_ctrl_type>>relative_error>>absolute_error>>this_stepsize;
                                                  //13.03.2008. WW
               int real_type = (int)(tsize_ctrl_type/10);
               if(real_type<10&&real_type>0)      //
               {
                  tsize_ctrl_type = real_type;
                  line >> h_min>> h_max >> max_time_step;
               }
               else
                  max_time_step = 0.0;
               line.clear();
            }
            // 12.03.2010 JTARON
            if(time_control_name=="COURANT")
            {
               line_string = GetLineFromFile1(tim_file);
               line.str(line_string);
               //	      courant_desired = desired Courant number
               //		  courant_initial = first time step size
               //        courant_static = 0  --> variable velocity (adapt in time)
               //        courant_static > 0  --> steady velocity (# of timesteps to recalculate...
               //                                                 i.e. =2, calculate first 2 time steps only, then use last value as constant in remaining simulation)
               line >> courant_desired >> courant_initial >> courant_static;
               line.clear();
            }
            // 26.08.2008. WW
            if(time_control_name=="STEP_SIZE_RESTRICTION")
            {
               line.str(GetLineFromFile1(tim_file));
               line >> h_min >> h_max;
               line.clear();
            }
            if(time_control_name=="NEUMANN")
            {
               line.clear();
            }
            if(time_control_name=="ERROR_CONTROL_ADAPTIVE")
            {
               m_pcs->adaption = true;
               line.clear();
            }
            if(time_control_name=="SELF_ADAPTIVE")
            {
               //m_pcs->adaption = true; JOD removed
               //WW minish = 10;
               while((!new_keyword)||(!new_subkeyword)||(!tim_file->eof()))
               {
                  position = tim_file->tellg();
                  line_string = GetLineFromFile1(tim_file);
                  if(line_string.find("#")!=std::string::npos)
                  {
                     return position;
                  }
                  if(line_string.find("$")!=std::string::npos)
                  {
                     new_subkeyword = true;
                     break;
                  }
                  if(line_string.find("MAX_TIME_STEP")!=std::string::npos)
                  {
                     *tim_file >> line_string;
                     max_time_step = strtod(line_string.data(),NULL);
                     line.clear();
                     // kg44 should not break break;
                  }
                  if(line_string.find("MIN_TIME_STEP")!=std::string::npos)
                  {
                     *tim_file >> line_string;
                     min_time_step = strtod(line_string.data(),NULL);
                     line.clear();
                     // kg44 should not break break;
                  }
                  /*  //WW
                  if(line_string.find("MINISH")!=string::npos){
                  *tim_file >> line_string;
                  minish = strtod(line_string.data(),NULL);
                  line.clear();
                  }
                  */
                  if(line_string.find("M")==std::string::npos)
                  {
                     line.str(line_string);
                     line >> iter_times;
                     line >> multiply_coef;
                     time_adapt_tim_vector.push_back(iter_times);
                     time_adapt_coe_vector.push_back(multiply_coef);
                     line.clear();
                  }
               }                                  // end of while loop adaptive
            }                                     // end of if "SELF_ADAPTIVE"
         }                                        // end of while
      }                                           // end of "TIME_CONTROL"
      //....................................................................
      /* //WW
      if(line_string.find("$SUBSTEPS")!=string::npos) { // subkeyword found JOD 4.7.10
        *tim_file>>sub_steps>>ws;
        continue;
      }
      */
      //....................................................................
   }                                              // end of while(!new_keyword)

   return position;
}


/**************************************************************************
FEMLib-Method:
Task: OBJ read function
Programing:
08/2004 OK Implementation
01/2005 OK Boolean type
01/2005 OK Destruct before read
**************************************************************************/
bool TIMRead(std::string file_base_name)
{
   //----------------------------------------------------------------------
   //OK  TIMDelete();
   //----------------------------------------------------------------------
   CTimeDiscretization *m_tim = NULL;
   char line[MAX_ZEILE];
   std::string sub_line;
   std::string line_string;
   std::ios::pos_type position;
   //========================================================================
   // File handling
   std::string tim_file_name = file_base_name + TIM_FILE_EXTENSION;
   std::ifstream tim_file (tim_file_name.data(),std::ios::in);
   if (!tim_file.good())
      return false;
   tim_file.seekg(0L,std::ios::beg);
   //========================================================================
   // Keyword loop
   std::cout << "TIMRead" << std::endl;
   while (!tim_file.eof())
   {
      tim_file.getline(line,MAX_ZEILE);
      line_string = line;
      if(line_string.find("#STOP")!=std::string::npos)
         return true;
      //----------------------------------------------------------------------
                                                  // keyword found
      if(line_string.find("#TIME_STEPPING")!=std::string::npos)
      {
         m_tim = new CTimeDiscretization();
         position = m_tim->Read(&tim_file);
         m_tim->time_current = m_tim->time_start;
         //----------------------------------------------------------------------
         if(m_tim->Write_tim_discrete)            //YD Write out Time Steps & Iterations
         {
            std::string m_file_name = file_base_name + "_TimeDiscrete.txt";
            m_tim->tim_discrete = new std::fstream(m_file_name.c_str(),std::ios::trunc|std::ios::out);
            std::fstream *tim_dis =  m_tim->tim_discrete;
            *tim_dis << " No. Time  Tim-Disc  Iter"<< std::endl;
            if (!m_tim->tim_discrete->good())
               std::cout << "Warning : Time-Discrete files are not found" << std::endl;
         }
         //----------------------------------------------------------------------
         time_vector.push_back(m_tim);
         tim_file.seekg(position,std::ios::beg);
      }                                           // keyword found
   }                                              // eof
   return true;
}


/**************************************************************************
FEMLib-Method:
Task: master write function
01/2005 OK Implementation
06/2009 OK Write only existing data
**************************************************************************/
void TIMWrite(std::string base_file_name)
{
   //----------------------------------------------------------------------
   if((int)time_vector.size()<1)
      return;
   //----------------------------------------------------------------------
   CTimeDiscretization *m_tim = NULL;
   std::string sub_line;
   std::string line_string;
   //========================================================================
   // File handling
   std::string tim_file_name = base_file_name + TIM_FILE_EXTENSION;
   std::fstream tim_file (tim_file_name.data(),std::ios::trunc|std::ios::out);
   tim_file.setf(std::ios::scientific,std::ios::floatfield);
   tim_file.precision(12);
   if (!tim_file.good()) return;
   tim_file.seekg(0L,std::ios::beg);
   //========================================================================
   // OUT vector
   tim_file << "GeoSys-TIM: ------------------------------------------------\n";
   for(int i=0;i<(int)time_vector.size();i++)
   {
      m_tim = time_vector[i];
      m_tim->Write(&tim_file);
   }
   tim_file << "#STOP";
   tim_file.close();
}


/**************************************************************************
FEMLib-Method:
01/2004 OK Implementation
05/2009 OK $TIME_CONTROL
**************************************************************************/
void CTimeDiscretization::Write(std::fstream*tim_file)
{
   int i;
   //--------------------------------------------------------------------
   // KEYWORD
   *tim_file  << "#TIME_STEPPING" << std::endl;
   //--------------------------------------------------------------------
   // PCS_TYPE
   *tim_file  << " $PCS_TYPE" << std::endl;
   *tim_file  << "  " << pcs_type_name << std::endl;
   //--------------------------------------------------------------------
   *tim_file  << " $TIME_START" << std::endl;
   *tim_file  << "  " << time_start << std::endl;
   //--------------------------------------------------------------------
   *tim_file  << " $TIME_END" << std::endl;
   *tim_file  << "  " << time_end << std::endl;
   //--------------------------------------------------------------------
   if(time_control_name.size()==0)
   {
      *tim_file  << " $TIME_STEPS" << std::endl;
      for(i=0;i<(int)time_step_vector.size();i++)
      {
         *tim_file  << "  " << 1 << " " << time_step_vector[i] << std::endl;
      }
   }
   //--------------------------------------------------------------------
   if(time_control_name.size()>0)
   {
      *tim_file  << " $TIME_CONTROL" << std::endl;
      if(time_control_name=="COURANT_MANIPULATE")
      {
         *tim_file  << "  " << time_control_name << std::endl;
         *tim_file  << "   " << time_control_manipulate << std::endl;
      }
      if(time_control_name=="PI_AUTO_STEP_SIZE")
      {
         *tim_file  << "  " << time_control_name << std::endl;
         *tim_file  << "   " << tsize_ctrl_type << " " << relative_error << " " << absolute_error << " " << this_stepsize << std::endl;
      }
      if(time_control_name=="STEP_SIZE_RESTRICTION")
      {
         *tim_file  << "  " << time_control_name << std::endl;
         *tim_file  << "   " << h_min << " " << h_max << std::endl;
      }
      if(time_control_name=="NEUMANN")
      {
         *tim_file  << "  " << time_control_name << std::endl;
      }
      if(time_control_name=="ERROR_CONTROL_ADAPTIVE")
      {
         *tim_file  << "  " << time_control_name << std::endl;
      }
      if(time_control_name=="SELF_ADAPTIVE")
      {
         *tim_file  << "  MAX_TIME_STEP " << max_time_step << std::endl;
         *tim_file  << "  MIM_TIME_STEP " << min_time_step << std::endl;
      }
   }
   //--------------------------------------------------------------------
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
08/2004 OK Implementation
08/2008 WW Force t+dt be indentical to the time for output or other special time
        WW Auto time step size control
**************************************************************************/
double CTimeDiscretization::CalcTimeStep(double crt_time)
{
   // time_step_length = 0.0;
   int no_time_steps = (int)time_step_vector.size();
   if(no_time_steps>0)
      time_step_length = time_step_vector[0];
   // Standard time stepping
   if(step_current<no_time_steps)
   {
      time_step_length = time_step_vector[step_current];
   }
   // Time step controls
   if( (time_control_name == "NEUMANN" ) || (time_control_name == "SELF_ADAPTIVE" ))
   {
      if(aktuelle_zeit < MKleinsteZahl && repeat == false)
         time_step_length = FirstTimeStepEstimate();
      else if( time_control_name == "NEUMANN" )
         time_step_length = NeumannTimeControl();
      else if(time_control_name == "SELF_ADAPTIVE")
         time_step_length = SelfAdaptiveTimeControl();
   }
   if(time_control_name == "ERROR_CONTROL_ADAPTIVE")
   {
      if(aktuelle_zeit < MKleinsteZahl)
         time_step_length = AdaptiveFirstTimeStepEstimate();
      else
         time_step_length = ErrorControlAdaptiveTimeControl();
   }
   if(time_control_name == "PI_AUTO_STEP_SIZE")   //WW
      time_step_length = this_stepsize;

   if(time_control_name == "COURANT")             // JTARON 2010, a more straight-forward Courant control
   {
      if(step_current==0)                         // first time step, set an initial
         time_step_length = courant_initial;
      else if(courant_static==0)
         time_step_length = CourantTimeControl();
      else if(step_current < courant_static+1)
      {
         time_step_length = CourantTimeControl();
         courant_initial = time_step_length;
      }
      else
         time_step_length = courant_initial;      // take constant value from step_current = 0
   }

   //--------
   //WW. Force t+dt be indentical to the time for output or other special time
   // 25.08.2008
   for(int i=0; i<(int)critical_time.size(); i++)
   {
      if((crt_time<critical_time[i])&&(crt_time+time_step_length>critical_time[i]))
      {
                                                  // _new, 23.09.2008.
         if(fabs( critical_time[i]-crt_time)>DBL_EPSILON)
            time_step_length = critical_time[i]-crt_time;
         std::cout << "Time step set to " << time_step_length << " in order to match critical times!"<< std::endl;
         break;
      }
   }
   //
   return time_step_length;
}


/**************************************************************************
FEMLib-Method: Operator
Task:
Programing:
08/2008 WW Implementation
**************************************************************************/
CTimeDiscretization::CTimeDiscretization(const CTimeDiscretization& a_tim, std::string pcsname)
{
   int i;
   safty_coe = a_tim.safty_coe;
   dt_sum = a_tim.dt_sum;
   this_stepsize = a_tim.this_stepsize;
   file_base_name = a_tim.file_base_name;
   time_start = a_tim.time_start;
   time_end = a_tim.time_end;
   time_current = a_tim.time_current;
   time_control_manipulate = a_tim.time_control_manipulate;
   step_current = a_tim.step_current;
   repeat = a_tim.repeat;
   pcs_type_name = pcsname;                       // by argument
   time_type_name = a_tim.time_type_name;
   time_control_name = a_tim.time_control_name;
   time_unit = a_tim.time_unit;
   iter_times = a_tim.iter_times;
   multiply_coef = a_tim.multiply_coef;
   max_time_step = a_tim.max_time_step;
   min_time_step = a_tim.min_time_step;
   Write_tim_discrete = a_tim.Write_tim_discrete;
   tim_discrete = a_tim.tim_discrete;
   nonlinear_iteration_error = a_tim.nonlinear_iteration_error;
   //
   time_step_vector.clear();
   time_adapt_tim_vector.clear();
   time_adapt_coe_vector.clear();
   for(i=0; i<(int)a_tim.time_step_vector.size(); i++)
      time_step_vector.push_back(a_tim.time_step_vector[i]);
   for(i=0; i<(int)a_tim.time_adapt_tim_vector.size(); i++)
      time_adapt_tim_vector.push_back(a_tim.time_adapt_tim_vector[i]);
   for(i=0; i<(int)a_tim.time_adapt_coe_vector.size(); i++)
      time_adapt_coe_vector.push_back(a_tim.time_adapt_coe_vector[i]);
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
12/2004 OK Implementation
last modified:
**************************************************************************/
                                                  //kg44 const string gave trouble for me
CTimeDiscretization* TIMGet(const std::string &pcs_type_name)
{
   CTimeDiscretization *m_tim = NULL;
   int i;
   int no_times = (int)time_vector.size();
   for(i=0;i<no_times;i++)
   {
      m_tim = time_vector[i];
      if(m_tim->pcs_type_name.compare(pcs_type_name)==0)
         return time_vector[i];
   }
   return NULL;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2005 OK Implementation
last modified:
**************************************************************************/
void TIMDelete()
{
   long i;
   int no_tim =(int)time_vector.size();
   for(i=0;i<no_tim;i++)
   {
      delete time_vector[i];
   }
   time_vector.clear();
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2005 OK Implementation
last modified:
**************************************************************************/
void TIMDelete(std::string pcs_type_name)
{
   long i;
   CTimeDiscretization *m_tim = NULL;
   int no_tim =(int)time_vector.size();
   for(i=0;i<no_tim;i++)
   {
      m_tim = TIMGet(pcs_type_name);
      if(!m_tim)                                  //OK
         continue;
      if(m_tim->pcs_type_name.compare(pcs_type_name)==0)
      {
         delete time_vector[i];
         time_vector.erase(time_vector.begin()+i);
      }
   }
}


/**************************************************************************
FEMLib-Method:
Task: Neumann estimation
Programing:
10/2005 YD Implementation
**************************************************************************/
double CTimeDiscretization::FirstTimeStepEstimate(void)
{
   CMediumProperties* m_mmp = NULL;
   CRFProcess* m_pcs = NULL;
   CElem* elem = NULL;
   int idxS;
   long group;
   double GP[3];
   static double Node_Sat[8];
   double buffer;
   int no_time_steps;
   //WW  int no_processes =(int)pcs_vector.size();
   CFluidProperties *m_mfp = NULL;
   m_mfp = MFPGet("LIQUID");                      //WW
   double density_fluid = m_mfp->Density();       //WW

   for (size_t n_p = 0; n_p < pcs_vector.size(); n_p++)
   {
      m_pcs = pcs_vector[n_p];
      CFiniteElementStd* fem = m_pcs->GetAssembler();

      time_step_length = min_time_step;           // take min time step as conservative best guess for testing
      //		switch (m_pcs->pcs_type_name[0]) {
      switch (m_pcs->getProcessType())            // TF
      {
         //		case 'G': // kg44 groudnwater flow ---if steady state, time step should be greater zero...transient flow does not work with adaptive stepping
         case GROUNDWATER_FLOW:                   // TF, if steady state, time step should be greater zero...transient flow does not work with adaptive stepping
            time_step_length = min_time_step;     // take min time step as conservative best guess for testing
            break;
            //		case 'M': // kg44 Mass transport ---if steady state, time step should be greater zero..
         case MASS_TRANSPORT:                     // TF, if steady state, time step should be greater zero..
            time_step_length = min_time_step;     // take min time step as conservative best guess for testing
            break;
            //		case 'R': // Richards
         case RICHARDS_FLOW:                      // TF
         {
            idxS = m_pcs->GetNodeValueIndex("SATURATION1");
            no_time_steps = 1000000000;           //OK (int)(1.0e10);
            time_step_length = 1.e10;
            size_t mmp_vector_size = mmp_vector.size();
            for (size_t m = 0; m < mmp_vector_size; m++)
               m_mmp = mmp_vector[m];
            const size_t size (m_pcs->m_msh->ele_vector.size());
            for (size_t i = 0; i < size; i++)
            {
               elem = m_pcs->m_msh->ele_vector[i];
               if (elem->GetMark())               // Element selected
               {
                  // Activated Element
                  group = elem->GetPatchIndex();
                  m_mmp = mmp_vector[group];
                  m_mmp->m_pcs = m_pcs;
                  MshElemType::type EleType = elem->GetElementType();
                                                  // Triangle
                  if (EleType == MshElemType::TRIANGLE)
                  {
                     GP[0] = GP[1] = 0.1 / 0.3;
                     GP[2] = 0.0;
                  } else if (EleType == MshElemType::TETRAHEDRON)
                  GP[0] = GP[1] = GP[2] = 0.25;
                  else
                     GP[0] = GP[1] = GP[2] = 0.0;
               }
               const int vertex_number (elem->GetVertexNumber());
               for (int j = 0; j < vertex_number; j++)
                  Node_Sat[j] = m_pcs->GetNodeValue(elem->GetNodeIndex(j),
                     idxS);
               buffer = m_mmp->SaturationPressureDependency(fem->interpolate(
                  Node_Sat), density_fluid, m_pcs->m_num->ls_theta);
               buffer *= 0.5 * elem->GetVolume() * elem->GetVolume();
               buffer *= m_mmp->porosity_model_values[0]
                  * mfp_vector[0]->Viscosity()
                  / m_mmp->permeability_tensor[0];
               buffer /= m_pcs->time_unit_factor;
               time_step_length = MMin(time_step_length, buffer);
            }                                     // ele_vector

            if (time_step_length < MKleinsteZahl)
            {
               std::cout << "Warning : Time Control Step Wrong, dt = 0.0 " << std::endl;
               time_step_length = 1.e-6;
            }
            std::cout << "Neumann Time Step: " << time_step_length << std::endl;
            time_step_length_neumann = 1.e10;
            time_step_length = MMin(time_step_length, max_time_step);
            if (Write_tim_discrete)
               *tim_discrete << aktueller_zeitschritt << "  " << aktuelle_zeit
                  << "   " << time_step_length << "  " << m_pcs->iter_lin
                  << std::endl;
            break;
         }
         default:
            std::cout << "CTimeDiscretization::FirstTimeStepEstimate default case" << std::endl;
      }
   }
   return time_step_length;
}


/**************************************************************************
FEMLib-Method:
Task: Module controls adaptive timestep based on a simple Courant condition
Programing:
12.03.2010 JTARON implementation
**************************************************************************/
double CTimeDiscretization::CourantTimeControl(void)
{
   CRFProcess *m_pcs = NULL;
   CFEMesh* m_msh = NULL;
   ElementValue* gp_ele =NULL;
   CElem* m_ele =NULL;
   CNode* m_nod1 = NULL;
   CNode* m_nod2 = NULL;

   long iel;
   int i, j, k, nds;
   double vel[3], dx[3], dx_temp[3];
   double vel_mag, new_dt, char_len;
   double final_dt = 1.e-6;

   m_pcs = PCSGet(pcs_type_name);
   m_msh = FEMGet(pcs_type_name);

   // Loop over elements
   for (iel=0; iel<(long)m_msh->ele_vector.size();iel++)
   {
      // velocities first
      gp_ele = ele_gp_value[iel];
      gp_ele->GetEleVelocity(vel);                // the ELEMENTAL velocity vector

      // positions next
      m_ele = m_msh->ele_vector[iel];
      nds = m_ele->GetNodesNumber(false);
      for(i=0; i<nds; i++)                        // max x,y,z lengths of element
      {
         for(j=i+1; j<nds; j++)
         {
            m_nod1 = m_msh->nod_vector[m_ele->nodes_index[i]];
            m_nod2 = m_msh->nod_vector[m_ele->nodes_index[j]];

            dx_temp[0] = (m_nod1->X()-m_nod2->X())*(m_nod1->X()-m_nod2->X());
            dx_temp[1] = (m_nod1->Y()-m_nod2->Y())*(m_nod1->Y()-m_nod2->Y());
            dx_temp[2] = (m_nod1->Z()-m_nod2->Z())*(m_nod1->Z()-m_nod2->Z());

            if(i==0&&j==1)
               VCopy(dx,dx_temp,3);               //					copy dx_temp onto dx: first iteration
            else
            {
               for(k=0; k<3; k++)
                  if(dx_temp[k]>dx[k])
                     dx[k]=dx_temp[k];            //					take max dx values (x,y,z) subsequent iterations
            }
         }
      }
      for(i=0; i<3; i++)                          //							Hold square root until after internal loop, for efficiency only
         dx[i]=sqrt(dx[i]);


      //    vel_mag = MVectorlength(vel[0],vel[1],vel[2]); //	velocity magnitude
      vel_mag = sqrt (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
      for(i=0; i<3; i++)
         vel[i] = fabs(vel[i])/vel_mag;           //							normalized velocity
      char_len = PointProduction(vel,dx);         //			dot product, gives projection of element size onto velocity vector
      new_dt = courant_desired*char_len/vel_mag;  //		Courant timestep, this element

      if(iel==0)
         final_dt = new_dt;
      else if(new_dt<final_dt)                    //  take minimum value in domain
         final_dt = new_dt;
   }
   time_step_length = final_dt;

   if(time_step_length < MKleinsteZahl)
   {
      std::cout << "Fatal error: Courant time step is less than machine tolerance, setting dt = 1.e-6" << std::endl;
      time_step_length = 1.e-6;
   }

   if(Write_tim_discrete)
      *tim_discrete<<aktueller_zeitschritt<<"  "<<aktuelle_zeit<<"   "<<time_step_length<< "  "<<m_pcs->iter_lin <<std::endl;
   return time_step_length;
}


/**************************************************************************
FEMLib-Method:
Task: Nuemann Control
Programing:
10/2005 YD Implementation
**************************************************************************/
double CTimeDiscretization::NeumannTimeControl(void)
{
   CRFProcess* m_pcs = NULL;

   for (size_t n_p = 0; n_p < pcs_vector.size(); n_p++)
   {
      m_pcs = pcs_vector[n_p];
      //		switch (m_pcs->pcs_type_name[0]) {
      //				case 'R': // Richards
      switch (m_pcs->getProcessType())            // TF
      {
         case RICHARDS_FLOW:
            time_step_length = time_step_length_neumann;
            break;
         default:
            std::cout << "Fatal error: No valid PCS type" << std::endl;
            break;
      }
   }

   std::cout << "Neumann Time Step: " << time_step_length << std::endl;
   time_step_length_neumann = 1.e10;
   if (Write_tim_discrete)
      *tim_discrete << aktueller_zeitschritt << "  " << aktuelle_zeit
         << "   " << time_step_length << "  " << m_pcs->iter_lin << std::endl;
   return time_step_length;
}


/**************************************************************************
FEMLib-Method:
Task: Self adaptive method
Programing:
10/2005 YD Implementation
03/2008 HS KG Implementation for Groundwater flow and mass transport
10/2010 KG updates
**************************************************************************/
double CTimeDiscretization::SelfAdaptiveTimeControl ( void )
{
   int imflag=1, iprocs=0;
   int iterdum=1;
	double my_max_time_step=0.0;
   CRFProcess* m_pcs = NULL;                      //YDToDo: m_pcs should be member



  // First calculate maximum time step according to Neumann criteria
#ifdef GEM_REACT
  my_max_time_step= MMin(max_time_step,MaxTimeStep());
  std::cout<<"Self_Adaptive Time Step: max time step "<<my_max_time_step<< std::endl;
#else
// kg44 This does not work in this way with multiple mass tranport processes!
   if ( repeat )
   {
      std::cout << "   TIM step is repeated" << std::endl;
      m_pcs = PCSGet ( pcs_type_name );
      m_pcs->PrimaryVariableReload();
   }
#endif
                                                  // TF
   ProcessType pcs_type (convertProcessType (pcs_type_name));
   for (size_t n_p = 0; n_p < pcs_vector.size(); n_p++)
   {
       m_pcs = pcs_vector[n_p];
       
       if (m_pcs->getProcessType() == pcs_type)  //compare process type and type name from Tim object
       {
           iprocs++;
           //			switch (m_pcs->pcs_type_name[0]) {
           switch (m_pcs->getProcessType())         // TF
           {
               default:
                   std::cout << "Fatal error: No valid PCS type" << std::endl;
                   break;
                   //			case 'R': // Richards
               case RICHARDS_FLOW:                   // TF
                   if ( (imflag>0) && ( m_pcs->iter_lin  >= time_adapt_tim_vector[time_adapt_tim_vector.size()-1] ) )
                   {
                       imflag=0;
                       time_step_length = time_step_length * time_adapt_coe_vector[time_adapt_tim_vector.size()-1];
                   }
                   if ((imflag == 1) && (m_pcs->iter_lin <= time_adapt_tim_vector[0]))
                   {
                       imflag = 2;
                       time_step_length = time_step_length * time_adapt_coe_vector[0];
                   }
                   break;
                   //			case 'G': //Groundwater flow
               case GROUNDWATER_FLOW:                // TF
                   // iterdum=MMax(iterdum,m_pcs->iter);
                   imflag = 1;
                   if ( (imflag>0) && ( m_pcs->iter_lin  >= time_adapt_tim_vector[1] ) )
                   {
                       imflag=0; std::cout << "Self adaptive time step: to many iterations for Groundwater flow" << std::endl;
                   }
                   if (((imflag == 1) && (m_pcs->iter_lin <= time_adapt_tim_vector[0])))
                   {
                       imflag = 2;
                   }
                   break;
                   //			case 'M': // Mass transport
               case MASS_TRANSPORT:                  // TF
                   iterdum = std::max(iterdum, m_pcs->iter_lin);
                   if ( (imflag>0) && ( m_pcs->iter_lin  >= time_adapt_tim_vector[1] ) )
                   {
                       imflag=0; std::cout << "Self adaptive time step: to many iterations for Transport " << m_pcs->iter_lin << " "<< time_adapt_tim_vector[1] << std::endl;
                   }
                   if ( ((imflag == 1) && ( m_pcs->iter_lin  <= time_adapt_tim_vector[0] ) ))
                   {
                       imflag=2; 
                   }
                   break;
           }
       }
   }


	if (imflag==0 && (time_adapt_coe_vector[time_adapt_tim_vector.size()-1]*time_step_length )>min_time_step) {time_step_length = time_step_length * time_adapt_coe_vector[time_adapt_tim_vector.size()-1];} //timestep smaller
        else if (imflag==2){ time_step_length = time_step_length * time_adapt_coe_vector[0];} //timestep bigger

   // BUG my_max_time_step is not necessarily initialised
   time_step_length = MMin ( time_step_length,my_max_time_step );
   time_step_length = MMax ( time_step_length,min_time_step );

   std::cout<<"Self_Adaptive Time Step: "<<" imflag " << imflag << " dr " << time_step_length<<" max iterations: " << iterdum << " number of evaluated processes: " << iprocs << std::endl;
   if ( Write_tim_discrete )
      *tim_discrete<<aktueller_zeitschritt<<"  "<<aktuelle_zeit<<"   "<<time_step_length<< "  "<<m_pcs->iter_lin<<std::endl;
   //}
   return time_step_length;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:CMCD 03/2006
**************************************************************************/
double CTimeDiscretization::CheckCourant(void)
{
   long index;
   long group;
   double velocity[3]={0.,0.,0.};
   double porosity, vg, advective_velocity, length, courant;
   CRFProcess* m_pcs = NULL;
   m_pcs = PCSGetFluxProcess();
   if (!m_pcs) return 0.0;
   int pcs_no = m_pcs->pcs_number;
   CMediumProperties* m_mmp = NULL;
   CElem* elem =NULL;
   ElementValue* gp_ele;
   long critical_element_no = -1;
   double recommended_time_step = 0.0;
   double stable_time_step = 0.0;
   //  int edx;

   for (index=0;index< (long)m_pcs->m_msh->ele_vector.size();index++)
   {
      elem = m_pcs->m_msh->ele_vector[index];
      length = elem->GetRepLength();
      group = elem->GetPatchIndex();
      m_mmp = mmp_vector[group];
      m_mmp->m_pcs = m_pcs;
      porosity = m_mmp->Porosity(m_mmp->Fem_Ele_Std);
      gp_ele = ele_gp_value[index];               //to get gp velocities
      gp_ele->getIPvalue_vec(pcs_no, velocity);
      vg = MBtrgVec(velocity,3);
      advective_velocity = vg/porosity;
                                                  //kg44 avoid zero velocity..otherwise stable_time_step is a problem
      if (advective_velocity<DBL_EPSILON) advective_velocity=DBL_EPSILON;
      courant = dt * advective_velocity/length;
      elem->SetCourant(courant);
      //    edx = m_pcs->GetElementValueIndex("COURANT"); //kg44 does this work?
      //    m_pcs->SetElementValue(index,edx,courant);    // kg44 seems not to work
      stable_time_step = (1./courant)*dt;
      if (index == 0) recommended_time_step = stable_time_step;
      if (stable_time_step < recommended_time_step)
      {
         recommended_time_step = stable_time_step;
         critical_element_no = index;
      }
   }
   std::cout<<"Courant time step control, critical element = "<<critical_element_no<<" Recomended time step "<<recommended_time_step<<std::endl;
   return recommended_time_step;
}


/**************************************************************************
FEMLib-Method:
Task: Neumann estimation
Programing:
04/2006 YD Implementation
**************************************************************************/
double CTimeDiscretization::AdaptiveFirstTimeStepEstimate(void)
{
   CNumerics *m_num (num_vector[0]);
   CElem* elem = NULL;
   static double Node_p[8];
   double p_ini, buff = 0.0;
   int no_time_steps;
   safty_coe = 5.0;
   p_ini = 1.0e-10;

   for (size_t n_p = 0; n_p < pcs_vector.size(); n_p++)
   {
      CRFProcess* m_pcs = pcs_vector[n_p];
      CFiniteElementStd* fem = m_pcs->GetAssembler();

      //  switch(m_pcs->pcs_type_name[0]){
      switch (m_pcs->getProcessType())            // TF
      {
         //  case 'R': // Richards
         case RICHARDS_FLOW:                      // TF
         {
            int idxp = m_pcs->GetNodeValueIndex("PRESSURE1") + 1;
            no_time_steps = 1000000000;           //OK (int)(1e10);
            time_step_length = 1.e10;
            for (size_t i = 0; i < m_pcs->m_msh->ele_vector.size(); i++)
            {
               elem = m_pcs->m_msh->ele_vector[i];
               for (int j = 0; j < elem->GetVertexNumber(); j++)
                  Node_p[j]
                     = m_pcs->GetNodeValue(elem->GetNodeIndex(j), idxp);
               p_ini = MMax(fabs(fem->interpolate(Node_p)), p_ini);
            }
            buff = safty_coe * sqrt(m_num->nls_error_tolerance / p_ini);
            buff /= m_pcs->time_unit_factor;
            time_step_length = MMin(time_step_length, buff);
            if (time_step_length < MKleinsteZahl)
            {
               std::cout << "Warning : Time Control Step Wrong, dt = 0.0 " << std::endl;
               time_step_length = 1.0e-8;
            }
            std::cout << "Error Control Time Step: " << time_step_length << std::endl;
            if (Write_tim_discrete)
               *tim_discrete << aktueller_zeitschritt << "  " << aktuelle_zeit
                  << "   " << time_step_length << "  " << m_pcs->iter_lin
                  << std::endl;
            break;
         }
         //		case 'M': // kg44 mass transport
         case MASS_TRANSPORT:                     // TF
            time_step_length = min_time_step;     // take min time step as conservative best guess for testing
            break;
            //		case 'G': // kg44 groudnwater flow ---if steady state, time step should be greater zeor...transient flow does not work with adaptive stepping
         case GROUNDWATER_FLOW:                   // if steady state, time step should be greater zero ... transient flow does not work with adaptive stepping
            time_step_length = min_time_step;     // take min time step as conservative best guess for testing
            break;
         default:
            break;
      }
   }
   return time_step_length;
}


/**************************************************************************
FEMLib-Method:
Task: Error control adaptive method
Programing:
04/2006 YD Implementation
**************************************************************************/
double CTimeDiscretization::ErrorControlAdaptiveTimeControl(void)
{
   CRFProcess* m_pcs = NULL;
   double rmax = 5.0;
   double rmin = 0.5;
   double safty_coe = 0.8;

   for (size_t n_p = 0; n_p < pcs_vector.size(); n_p++)
   {
      m_pcs = pcs_vector[n_p];
      //		switch (m_pcs->pcs_type_name[0]) {
      switch (m_pcs->getProcessType())            // TF
      {
         default:
            std::cout << "Fatal error: No valid PCS type" << std::endl;
            break;
            //		case 'R': // Richards, accepted and refused time step
         case RICHARDS_FLOW:                      // accepted and refused time step
            //nonlinear_iteration_error = m_pcs->nonlinear_iteration_error;
            if (repeat)
            {
               time_step_length *= MMax(safty_coe * sqrt(
                  m_pcs->m_num->nls_error_tolerance
                  / nonlinear_iteration_error), rmin);
            }
            else
            {
               time_step_length *= MMin(safty_coe * sqrt(
                  m_pcs->m_num->nls_error_tolerance
                  / nonlinear_iteration_error), rmax);
            }
            std::cout << "Error_Self_Adaptive Time Step: " << time_step_length
               << std::endl;
            if (Write_tim_discrete)
               *tim_discrete << aktueller_zeitschritt << "  " << aktuelle_zeit
                  << "   " << time_step_length << "  " << m_pcs->iter_lin
                  << std::endl;
      }
   }
   return time_step_length;
}


/**************************************************************************
FEMLib-Method:
Task:  Check the time of the process in the case: different process has
       different time step size
Return boolean value: skip or execute the process
Programing:
06/2007 WW Implementation
09/2007 WW The varable of the time step accumulation as a  member
**************************************************************************/
double CTimeDiscretization::CheckTime(double const c_time, const double dt0)
{
   double pcs_step;
   double time_forward;
   bool ontime = false;
   if((int)time_vector.size()==1)
      return dt0;
   //
   //WW please check +1
   //OK   double pcs_step = time_step_vector[step_current+1];
   if(time_step_vector.size()>0)                  // 16.09.2008. WW
   {
                                                  //OK
      if(step_current>=(int)time_step_vector.size())
                                                  //OK
         pcs_step = time_step_vector[(int)time_step_vector.size()-1];
      else
                                                  //OK
            pcs_step = time_step_vector[step_current];
   }
   else
      pcs_step = this_stepsize;                   // 16.09.2008. WW
   time_forward = c_time - time_current-pcs_step;
   if(time_forward>0.0||fabs(time_forward)<MKleinsteZahl)
   {
      time_current += pcs_step;
      //WW. 02.02.2009    step_current++;
      this_stepsize = dt_sum+dt0;
      ontime = true;
      dt_sum = 0.0;
   }
   /*
   // HS-WW: 04.01.2010, bugfix, if not ontime, set this_stepsize to zero
   if ( time_forward < 0.0 && fabs(time_forward)>MKleinsteZahl)
   {
     //check if current time step is critical time
     bool isCriticalTime = false;
     for (int i=0; i<(int)critical_time.size(); i++) {
       if (critical_time[i]-c_time>MKleinsteZahl)
         break;
       if (fabs(c_time-critical_time[i])<MKleinsteZahl) {
         isCriticalTime = true;
   break;
   }
   }

   if (!isCriticalTime)
   this_stepsize = 0.0;
   }
   */
   if((fabs(pcs_step-time_end)<DBL_MIN)&&fabs(c_time-time_end)<DBL_MIN)
   {
      this_stepsize = dt_sum+dt0;
      ontime = true;
      dt_sum = 0.0;
   }
   if(!ontime)
   {
      dt_sum += dt0;
      //this_stepsize = 0.0;    //20.03.2009. WW
   }
   if(pcs_step>time_end)                          // make output for other processes
   {
      dt_sum = 0.0;
      this_stepsize = 0.0;
   }
   return this_stepsize;
}


/**************************************************************************
FEMLib-Method:
Task:  Used to force time steps matching the times requried by output or
       boundary
Programing:
08/2008 WW Implementation
**************************************************************************/
void CTimeDiscretization::FillCriticalTime()
{
   for (size_t i = 0; i < out_vector.size(); i++)
   {
      COutput *a_out = out_vector[i];
      for (size_t j = 0; j < a_out->getTimeVector().size(); j++)
      {
         bool done = false;
         for (size_t k = 0; k < critical_time.size(); k++)
         {
            if (fabs(critical_time[k] - a_out->getTimeVector()[j]) < DBL_MIN)
            {
               done = true;
               break;
            }
         }
         if (!done)
            critical_time.push_back(a_out->getTimeVector()[j]);
      }
   }
   // Sort
   for (size_t i = 0; i < critical_time.size(); i++)
   {
      for (size_t j = i; j < critical_time.size(); j++)
      {
         if (critical_time[i] > critical_time[j])
         {
            //				double val = critical_time[i];
            //				critical_time[i] = critical_time[j];
            //				critical_time[j] = val;
            std::swap (critical_time[i], critical_time[j]);
         }
      }
   }
}


/**************************************************************************
FEMLib-Method:
Programing:
09/2007 WW Implementation
**************************************************************************/
bool IsSynCron()
{
   int i, count = 0;
   for(i=0; i<(int)time_vector.size(); i++)
   {
      if(time_vector[i]->dt_sum<DBL_MIN)
         count++;
   }
   if(count==(int)time_vector.size())
      return true;
   else
      return false;
}


/**************************************************************************
FEMLib-Method:
Task:  construct time_step_target_vector from ic-/bc-curves (time curves)
Return boolean value:
Programing:
12/2007 KG44 Implementation
**************************************************************************/
/* bool CTimeDiscretization::GetTimeStepTargetVector() {

  bool have_vector=false;
  int no_times, i,j, anz;
  StuetzStellen *s = NULL;

   if (anz_kurven<=0) return have_vector;
// first get the time curves
    for (i;i<anz_kurven;i++) {
       anz = kurven[i].anz_stuetzstellen;
       s = kurven[i].stuetzstellen;
for (j;j<anz;j++){
time_step_target_vector.push_back(s[j].punkt);
}
}

return have_vector;
} */
#ifdef GEM_REACT
double CTimeDiscretization::MaxTimeStep()
{
   long i;
   double max_diff_time_step=1.0e+100,Dm,dummy,max_adv_time_step=1.0e+100;
   double theta=0.0;                              // direction zero...no anisotropy
   CRFProcess* this_pcs=NULL;
   CElem* melem=NULL;

   CMediumProperties *m_mat_mp = NULL;
   // Get the pointer to a proper PCS. ..we assume that all transport processes use the same diffusion coefficient

   this_pcs = PCSGet ( "MASS_TRANSPORT" );
   long nElems = ( long ) this_pcs->m_msh->ele_vector.size();
   int component = this_pcs->pcs_component_number;
   int group;

   CompProperties *m_cp = cp_vec[component];

   dummy=CheckCourant();                          // courant number
   std::cout << " Advective Time Step " << dummy << " " ;
                                                  //only do if Courant number bigger than zero
   if (dummy>DBL_EPSILON) max_adv_time_step=std::min(max_diff_time_step, dummy);

   // find Neumann for first mass transport process
   for(i=0;i<nElems;i++)
   {
      group = this_pcs->m_msh->ele_vector[i]->GetPatchIndex();
      m_mat_mp = mmp_vector[group];

      melem =  this_pcs->m_msh->ele_vector[i];
      //		cout << m_mat_mp->Porosity(i,theta) << " " << melem->representative_length << endl;
                                                  // KG44 attention DM needs to be multiplied with porosity!
      Dm = m_mat_mp->Porosity(i,theta) * m_cp->CalcDiffusionCoefficientCP(i,theta,this_pcs);
      // calculation of typical length

      dummy=( 0.5 * (melem->representative_length*melem->representative_length))/Dm;
      max_diff_time_step=std::min(max_diff_time_step, dummy);
      //	std::cout << "Neumann criteria: " << max_diff_time_step << " i " << i << std::endl;
   }

   std::cout << " Diffusive Time Step " << max_diff_time_step  << std::endl;

   return std::min(max_diff_time_step, max_adv_time_step);
}
#endif                                            // end of GEM_REACT
