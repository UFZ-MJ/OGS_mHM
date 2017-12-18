/**************************************************************************
FEMLib - Object: TIM
Task: class implementation
Programing:
06/2004 OK Implementation
12/2008 WW PI time step control
last modified:
**************************************************************************/
#ifndef rf_tim_new_INC
#define rf_tim_new_INC
// C++ STL
//#include <fstream>
#include <iostream>
//#include <string>
//#include <vector>

//----------------------------------------------------------------
class CRFProcess;                                 //21.08.2008. WW
class CTimeDiscretization
{
   private:
      double safty_coe;
      double dt_sum;                              // 17.09.2007 WW
      // For PI time step control. Aug-Nov.2008. by WW
      //Begin of data section for PI Time control ------------------------
   public:                                        //OK
      double this_stepsize;
      double relative_error;
      double absolute_error;
   private:
      double h_min;
      double h_max;
      double hacc;
      double erracc;
   public:                                        //OK
      int tsize_ctrl_type;
   private:
      std::vector<double> critical_time;
      //End of data section for PI Time control ------------------------
      friend bool IsSynCron();                    //WW
   public:
      std::string file_base_name;
      // TIM
      std::vector<double>time_step_vector;
      std::vector<int> time_adapt_tim_vector;
      std::vector<double>time_adapt_coe_vector;
      //WW vector<double>fixed_point_vector;

      //WW vector<double> time_step_target_vector; // kg44 for adaptive steps..intermediate time target that need to be reached
      double time_start;
      double time_end;
      double time_current;
      double time_control_manipulate;             //CMCD
      int step_current;
      bool repeat;                                //OK/YD
      // PCS
      std::string pcs_type_name;                  //OK
      // NUM
      std::string time_type_name;                 //OK
      std::string time_control_name;
      std::string time_unit;                      //WW
      double iter_times;                          //YD
      double multiply_coef;                       //YD
      double max_time_step;                       //YD
      double min_time_step;
      double courant_desired;                     //JTARON
      double courant_initial;                     //JTARON
      int courant_static;                         //JTARON
      //WW double minish; // JOD
      //WW int sub_steps; // JOD 4.7.10
      bool Write_tim_discrete;                    //YD
      std::fstream *tim_discrete;                 //YD
      double nonlinear_iteration_error;           //OK/YD
      //WW double max_adaptive_factor; // kg44
      //WW double max_adaptive_concentration_change; // kg44
   public:
      CTimeDiscretization(void);
      //21.08.2008. WW
      CTimeDiscretization(const CTimeDiscretization& a_tim, std::string pcsname);
      ~CTimeDiscretization(void);
      std::ios::pos_type Read(std::ifstream*);
      void Write(std::fstream*);
      double time_step_length_neumann;            //YD
      double time_step_length;                    //YD
      double CalcTimeStep(double crt_time=0.0);   // Add argument double crt_time. 25.08.2008. WW
      double FirstTimeStepEstimate();
      double AdaptiveFirstTimeStepEstimate();
      // For PI time step control. Aug-Nov.2008. by WW
      //Begin of function section for PI Time control ------------------------
      int GetTimeStepCrtlType() const {return tsize_ctrl_type;}
      double GetTimeStep() const {return this_stepsize;}
      void SetTimeStep( double hnew)  {this_stepsize = hnew;}
      double GetRTol() const { return relative_error;}
      double GetATol() const { return absolute_error;}
      double GetMaximumTSizeRestric() const { return h_max;}
      double GetMinimumTSizeRestric() const { return h_min;}
      double GetHacc() const { return hacc;}
      double GetErracc() const { return erracc;}
      void SetHacc(const double hacc_val)  { hacc = hacc_val;}
      void setErracc(const double erracc_val) { erracc = erracc_val;}
      double CheckTime(double const c_time, const double dt0);
      void FillCriticalTime();                    //21.08.2008.
      //Begin of function section for PI Time control ------------------------
      double ErrorControlAdaptiveTimeControl();
      double NeumannTimeControl();
      double CourantTimeControl();                // JTARON 2010
      double SelfAdaptiveTimeControl();
#ifdef GEM_REACT
      double MaxTimeStep();
#endif
      //
      //WW bool GetTimeStepTargetVector(); // kg44
      double CheckCourant();                      //CMCD
};

extern std::vector<CTimeDiscretization*> time_vector;
extern bool TIMRead(std::string);
// extern CTimeDiscretization* TIMGet(const string&);
                                                  //kg44 const string made trouble for me
extern CTimeDiscretization* TIMGet(const std::string& pcs_type_name);

extern void TIMWrite(std::string);
extern bool IsSynCron();                          //WW
extern void TIMDelete();
extern void TIMDelete(std::string);
#define TIM_FILE_EXTENSION ".tim"
//ToDo
extern double aktuelle_zeit;
extern size_t aktueller_zeitschritt;
extern double dt;
extern int rwpt_numsplits;                        // JTARON 2010, for specifying sub time step for random walker in .tim input file
#endif
