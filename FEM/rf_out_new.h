/**************************************************************************
FEMLib - Object: OUT
Task: class implementation
Programing:
06/2004 OK Implementation
last modified:
**************************************************************************/
#ifndef rf_out_new_INC
#define rf_out_new_INC
// C++ STL
//#include <list>
//#include <fstream>
//#include <string>
//#include <vector>

// FEMLib
//#include "rf_pcs.h"
//#include <sstream>        // for istringstream (ME)
#include "vtk.h"

// FEM
#include "GeoInfo.h"
#include "ProcessInfo.h"
#include "DistributionInfo.h"

namespace GEOLIB
{
   class GEOObjects;
}


namespace Mesh_Group{class CFEMesh;}
//using Mesh_Group::CFEMesh;

class COutput: public GeoInfo, public ProcessInfo, public DistributionInfo
{
   public:
      COutput();
      COutput (size_t id);
      /**
       * method initializes process and mesh attributes
       */
      void init ();
      ~COutput(void);

      /**
       * scaling factor for values
       * @param amplifier - a double value for scaling data
       */
      void setAmplifier(double amplifier)
      {
         out_amplifier = amplifier;
      }

      CRFProcess* GetPCS(const std::string&);     //OK
      CRFProcess* GetPCS();                       // 09/2010 TF
      CRFProcess* GetPCS_ELE(const std::string&); //OK

      /**
       * checking the consistency of the output data as specified in the input file
       * This means up to now, that data for missing processes is not written.
       */
      void checkConsistency();                    // CB (refactored by TF)

      void GetNodeIndexVector(std::vector<int>&); //OK
      void SetNODFluxAtPLY();                     //OK

      // ELE values
      const std::vector<std::string>& getElementValueVector () const { return _ele_value_vector; }
                                                  //OK
      void GetELEValuesIndexVector(std::vector<int>&);

      /**
       *
       * @return
       */
      const std::vector<std::string>& getRandomWalkParticleTracingValueVector() const { return _rwpt_value_vector; }

      /**
       * ToDo remove after transition to new GEOLIB - REMOVE CANDIDATE
       * getGeoName returns a string used as id for geometric entity
       * @return the value of attribute geo_name in case of
       * geo_type_name == POLYLINE or geo_type_name = SURFACE
       * If geo_type_name == POINT the id of the point is returned.
       */
      const std::string& getGeoName() const;      // TF 05/2010

      CFEMesh* getMesh ()                         // TF
      {
         return m_msh;
      }

      /**
       * read from file stream
       * @param in input file stream
       * @param geo_obj object of class GEOObjects that manages the geometric entities
       * @param unique_name the name of the project to access the right geometric entities
       * @return the new position in the stream after reading
       */
      std::ios::pos_type Read(std::ifstream& in, const GEOLIB::GEOObjects& geo_obj,
         const std::string& unique_name);

      void Write(std::fstream*);
      // TF not used (at the moment?) REMOVE CANDIDATE
      //    int GetPointClose(CGLPoint);
      void WriteTimeCurveData(std::fstream &);
      void WriteTimeCurveHeader(std::fstream &);
      void NODWriteDOMDataTEC();
      void WriteTECHeader(std::fstream&, int, std::string);
      void WriteTECNodeData(std::fstream&);
      void WriteTECElementData(std::fstream&, int);
      double NODWritePLYDataTEC(int);
      void NODWritePNTDataTEC(double, int);
      void ELEWriteDOMDataTEC();
      void WriteELEValuesTECHeader(std::fstream&);
      void WriteELEValuesTECData(std::fstream&);
      void NODWriteSFCDataTEC(int);
      void NODWriteSFCAverageDataTEC(double, int);//OK
      void WriteDataVTK(int);                     //GK
      void WriteVTKHeader(std::fstream&, int);
      void WriteVTKNodeData(std::fstream&);
      void WriteVTKElementData(std::fstream&);
      void WriteVTKValues(std::fstream&);
      void WriteRFO();                            //OK
      void WriteRFOHeader(std::fstream&);         //OK
      void WriteRFONodes(std::fstream&);          //OK
      void WriteRFOElements(std::fstream&);       //OK
      void WriteRFOValues(std::fstream&);         //OK
      void NODWriteLAYDataTEC(int);               //OK
      void ELEWriteSFC_TEC();                     //OK
      void ELEWriteSFC_TECHeader(std::fstream&);  //OK
      void ELEWriteSFC_TECData(std::fstream&);    //OK
      void CalcELEFluxes();                       //OK
      void ELEWritePLY_TEC();                     //OK
      void ELEWritePLY_TECHeader(std::fstream&);  //OK
      void ELEWritePLY_TECData(std::fstream&);    //OK
      void TIMValue_TEC(double);                  //OK
      double NODFlux(long);                       //OK
      void PCONWriteDOMDataTEC();                 //MX
      void WriteTECNodePCONData(std::fstream &);  //MX

      void setTime (double time) { _time = time; }
      /**
       * get time returns the value of attribute time
       * @return
       */
      double getTime () const { return _time; }

      const std::vector<double>& getTimeVector () const { return time_vector; }
      std::string& getFileBaseName () { return file_base_name; }

      size_t getNSteps () const { return nSteps; }
      /**
       * constructs/adds the output file name using geo_name,
       * process type, mesh type
       * @param fname a reference to the constructed file name
       * @param geo switch on/off geo info in file name (default = on)
       * @param process switch on/off process info in file name (default = on)
       * @param mesh switch on/off mesh info in file name (default = on)
       */
      void addInfoToFileName(std::string& fname, bool geo = true, bool process =
         true, bool mesh = true) const;           // 09/2010 TF

      std::vector<std::string> _nod_value_vector;
      // MAT values
      std::vector<std::string> mmp_value_vector;  //OK
      std::vector<std::string> mfp_value_vector;  //OK

      CRFProcess* m_pcs;                          //OK

      //	std::vector<double>& getRWPTTimeVector () { return rwpt_time_vector; }
      std::vector<double>& getRWPTTimeVector () { return time_vector; }

   private:
      friend void OUTData(double, int step);

      //	std::vector<double> rwpt_time_vector; //JTARON, needed because outputs are treated differently in RWPT

      // MSH
      std::string msh_type_name;                  //OK

      // TIM
      std::string tim_type_name;                  // STEPS or TIMES ?
      std::vector<double> time_vector;
      double _time;

      /**
       * the position in the global vector out_vector, used only in NODWritePLYDataTEC
       */
      size_t _id;

      std::string file_base_name;
      double out_amplifier;                       //WW to amplify output
                                                  //WW/OK
      inline void WriteELEVelocity(std::iostream &vtk_file);

      CFEMesh* m_msh;
      int nSteps;                                 // After each nSteps, make output

      CVTK* vtk;
      // GEO
      /**
       * the id of the geometric object as string REMOVE CANDIDATE
       */
      std::string geo_name;                       // TF 05/2010

      // File status
      bool _new_file_opened;                      //WW

      // DAT
      /**
       * this attribute stores the output format
       */
      std::string dat_type_name;

      // ELE value
      std::vector<std::string> _ele_value_vector;

      // RWPT values
      std::vector<std::string> _rwpt_value_vector;

      // PCON values
      std::vector<std::string> _pcon_value_vector;
};

extern std::vector<COutput*> out_vector;

/**
 * read file that stores information about output
 * @param file_base_name base file name (without extension)
 * @param geo_obj object of class GEOObjects managing the geometric entities
 * @param unique_name unique name to access the geometric entities in geo_obj
 * @return true if file reading was successful, else false
 */
bool OUTRead(const std::string& file_base_name, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);

extern void OUTWrite(std::string);
#define OUT_FILE_EXTENSION ".out"
extern void OUTData(double, const int step);
extern void OUTDelete();
extern COutput* OUTGet(const std::string &);
extern void OUTCheck(void);                       // new SB
extern COutput* OUTGetRWPT(const std::string &);  //JTARON
#endif
