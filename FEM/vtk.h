#ifndef VTK_INC
#define VTK_INC

#include <string>
#include <vector>
//#include "rf_out_new.h"
#include "MSHEnums.h"

class COutput;

namespace Mesh_Group
{
   class CFEMesh;
}

using namespace Mesh_Group;

typedef struct
{
   double timestep;
   std::string vtk_file;
} VTK_Info;

class CVTK
{
   public:
      std::vector<VTK_Info> vec_dataset;
      std::string pvd_file_name;
      std::string pvd_vtk_file_name_base;
      double useBinary;

      enum VTK_XML_DATA_TYPE { Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32, Float64 };

   protected:
      //for binary output
      bool isInitialized;
      VTK_XML_DATA_TYPE type_UChar;
      VTK_XML_DATA_TYPE type_Int;
      VTK_XML_DATA_TYPE type_UInt;
      VTK_XML_DATA_TYPE type_Long;
      VTK_XML_DATA_TYPE type_Double;
      int SIZE_OF_BLOCK_LENGTH_TAG;
      bool isLittleEndian;                        //Endian(byte order)

   public:
      CVTK(void){isInitialized=false;};
      virtual ~CVTK(void){};

   protected:
      //PVD
      bool WriteHeaderOfPVD(std::fstream &fin);
      bool WriteEndOfPVD(std::fstream &fin);
      bool WriteDatasetOfPVD(std::fstream &fin, double timestep, const std::string &vtkfile);

      //VTU
      void InitializeVTU();
      unsigned char GetVTKCellType(const MshElemType::type ele_type);
      bool WriteDataArrayHeader(std::fstream &fin, VTK_XML_DATA_TYPE data_type, const std::string &str_name, int nr_components, const std::string &str_format, long offset=-1);
      bool WriteDataArrayFooter(std::fstream &fin);
      inline bool WriteMeshNodes(std::fstream &fin, bool output_data, CFEMesh *m_msh, long &offset);
      inline bool WriteMeshElementConnectivity(std::fstream &fin, bool output_data, CFEMesh *m_msh, long &offset, long &sum_ele_components);
      inline bool WriteMeshElementOffset(std::fstream &fin, bool output_data, CFEMesh *m_msh, long &offset);
      inline bool WriteMeshElementType(std::fstream &fin, bool output_data, CFEMesh *m_msh, long &offset);
      inline bool WriteNodalValue(std::fstream &fin, bool output_data, COutput *out, CFEMesh *m_msh, long &offset);
      inline bool WriteElementValue(std::fstream &fin, bool output_data, COutput *out, CFEMesh *m_msh, long &offset);

      //util
      template <typename T> void write_value_binary(std::fstream &fin, T val);
      bool IsLittleEndian();

   public:
      //PVD
      bool InitializePVD(const std::string &file_base_name, const std::string &pcs_type_name, bool binary=false);
      bool UpdatePVD(const std::string &pvdfile, const std::vector<VTK_Info> &vec_vtk);
      bool CreateDirOfPVD(const std::string &pvdfile);

      //VTU
      bool WriteXMLUnstructuredGrid(const std::string &vtkfile, COutput *out, const int time_step_number);

};
#endif
