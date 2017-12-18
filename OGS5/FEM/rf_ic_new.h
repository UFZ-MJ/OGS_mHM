/**************************************************************************
FEMLib - Object: Initial Conditions IC
Task: class implementation
Programing:
08/2004 OK Implementation
last modified
**************************************************************************/
#ifndef rf_ic_new_INC
#define rf_ic_new_INC

#define IC_FILE_EXTENSION ".ic"

// C++ STL
//#include <fstream>
//#include <string>
//#include <vector>

// FEM
#include "GeoInfo.h"                              // TF
#include "ProcessInfo.h"                          // KR
#include "DistributionInfo.h"                     // TF

//#include "rf_pcs.h"

class CNodeValue;

/**
 * class for handling initial conditions
 */
class CInitialCondition : public ProcessInfo, public GeoInfo, public DistributionInfo
{
   private:
      size_t SubNumber;                           //WW
      std::vector<int> subdom_index;              //WW
      std::vector<double> subdom_ic;              //WW
      // Coefficents for linear distribution function
      // f = a0+b0*x+c0*y+d0*z
      double *a0, *b0, *c0, *d0;                  //WW
      std::string fname;                          //17.11.2009. PCH
   private:
      // REMOVE CANDIDATE
      std::string geo_name;                       // TF 05/2010
   public:
      int GetNumDom() const
      {
         return (int) subdom_index.size();
      }                                           //WW
      int GetDomain(const int dom_index) const
      {
         return subdom_index[dom_index];
      }                                           //WW
      //int mat_type; //MX
      // DIS
      std::vector<CNodeValue*> node_value_vector;
      void SetDomain(int);
      void SetByNodeIndex(int);                   // 19.11.2009 PCH
      void SetPolyline(int);
      void SetSurface(int);
      double DistributionFuntion(const int dom_i, const double x, const double y,
         const double z)                          //WW
      {
         return a0[dom_i] + b0[dom_i] * x + c0[dom_i] * y + d0[dom_i] * z;
      }
      void SetPoint(int);                         //MX
      //void SetMaterialDomain(int); //MX
      double gradient_ref_depth;
      double gradient_ref_depth_value;
      double gradient_ref_depth_gradient;
      std::string rfr_file_name;                  //OK
      CInitialCondition();
      ~CInitialCondition();
      /**
       * read initial condition from stream
       * @param in input stream from file
       * @param geo_obj object of class GEOObjects that manages the geometric entities
       * @param unique_name the name of the project to access the right geometric entities
       * @return the new position in the stream after reading
       */
      std::ios::pos_type Read(std::ifstream* in, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);
      void Write(std::fstream*) const;
      void Set(int);
      void SetEle(int);                           //MX
      void SetDomainEle(int);                     //MX
	  Mesh_Group::CFEMesh* m_msh;
      CNodeValue* m_node;
};

class CInitialConditionGroup
{
   public:
      std::string pcs_type_name;                  //OK
      std::string pcs_pv_name;                    //OK
      std::vector<CNodeValue*>group_vector;
};

extern std::vector<CInitialConditionGroup*>ic_group_vector;
extern std::vector<CInitialCondition*>ic_vector;
/**
 * read file that stores initial conditions
 * @param file_base_name base file name (without extension) containing the initial conditions
 * @param geo_obj object of class GEOObjects managing the geometric entities
 * @param unique_name unique name to access the geometric entities in geo_obj
 * @return true if initial conditions found in file, else false
 */
bool ICRead(const std::string& file_base_name, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);
extern void ICWrite(std::string);
extern void ICDelete();
extern CInitialCondition* ICGet(std::string);     //OK
#endif
