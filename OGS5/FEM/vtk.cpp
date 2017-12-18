#include "vtk.h"
#include <fstream>
#if defined(WIN32)
#include <direct.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "makros.h"
#include "rf_mmp_new.h"
#include "fem_ele_std.h"                          // for element velocity

using namespace std;

const std::string INDEX_STR = "  ";
const std::string velocity_name[3][4] =
{
   {
      "VELOCITY_X1", "VELOCITY_Y1", "VELOCITY_Z1", "NODAL_VELOCITY1"
   }
   ,
   {
      "VELOCITY_X2", "VELOCITY_Y2", "VELOCITY_Z2", "NODAL_VELOCITY2"
   }
   ,
   {
      "VELOCITY1_X", "VELOCITY1_Y", "VELOCITY1_Z", "GL_NODAL_VELOCITY1"
   }
};

//#################################################################################################
// Functions for Paraview Data File (PVD)

bool CVTK::InitializePVD(const string &file_base_name, const string &pcs_type_name, bool binary)
{
   //PVD
   this->vec_dataset.clear();
   this->pvd_file_name = file_base_name;
   if(pcs_type_name.size()>0)                     // PCS
      this->pvd_file_name += "_" + pcs_type_name;
   this->pvd_file_name += ".pvd";
   //VTK
   int ibs = (int)file_base_name.find_last_of("\\");
   int is = (int)file_base_name.find_last_of("/");
                                                  //OK411
   if (ibs != (int)string::npos  || is != (int)string::npos)
   {
      int ibegin = ibs; if (is > ibs) ibegin = is;
      ibegin+=1;
      this->pvd_vtk_file_name_base = file_base_name.substr(ibegin);
   }
   else
   {
      this->pvd_vtk_file_name_base = file_base_name;
   }
   if (pcs_type_name.size()>0)                    // PCS
      this->pvd_vtk_file_name_base += "_" + pcs_type_name;

   this->useBinary = binary;

   return true;
}

bool CVTK::WriteHeaderOfPVD(std::fstream &fin)
{
   fin << "<?xml version=\"1.0\"?>" << std::endl;
   fin << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"  << std::endl;
   fin << INDEX_STR << "<Collection>" << std::endl;
   return true;
}

bool CVTK::WriteEndOfPVD(std::fstream &fin)
{
   fin << INDEX_STR << "</Collection>" << std::endl;
   fin << "</VTKFile>" << std::endl;
   return true;
}

bool CVTK::WriteDatasetOfPVD(std::fstream &fin, double timestep, const std::string &vtkfile)
{
   fin.setf(ios::scientific,std::ios::floatfield);
   fin.precision(12);
   fin << INDEX_STR << INDEX_STR << "<DataSet timestep=\"" << timestep << "\" group=\"\" part=\"0\" file=\"" << vtkfile << "\"/>" << std::endl;
   return true;
}


bool CVTK::UpdatePVD(const string &pvdfile, const vector<VTK_Info> &vec_vtk)
{

   fstream fin(pvdfile.data(), ios::out);
   if (!fin.good()) return false;

   CVTK::WriteHeaderOfPVD(fin);
   for (int i=0; i<(int)vec_vtk.size(); i++)
      CVTK::WriteDatasetOfPVD(fin, vec_vtk[i].timestep, vec_vtk[i].vtk_file);
   CVTK::WriteEndOfPVD(fin);

   fin.close();

   return true;
}


bool CVTK::CreateDirOfPVD(const string &pvdfile)
{
   string pvd_dir_path = pvdfile + ".d";
#if defined(WIN32)
   if (mkdir(pvd_dir_path.c_str()) == -1)
   {
#else
      if (mkdir(pvd_dir_path.c_str(), 0777) == -1)
      {
#endif
         //error
         cout << "***ERROR: Fail to create a PVD directory: " << pvd_dir_path << endl;
         return false;
      }
      return true;
   }

   //#################################################################################################
   // Functions for VTU files (XML UnstructuredGrid file)

   unsigned char CVTK::GetVTKCellType(const MshElemType::type ele_type)
   {
      unsigned char cell_type = 0;

      switch(ele_type)
      {
         case MshElemType::LINE:                  // vtk_line=3
            cell_type = 3;
            break;
         case MshElemType::QUAD:                  // quadrilateral=9
            cell_type = 9;
            break;
         case MshElemType::HEXAHEDRON:            // hexahedron=12
            cell_type = 12;
            break;
         case MshElemType::TRIANGLE:              // triangle=5
            cell_type = 5;
            break;
         case MshElemType::TETRAHEDRON:           // tetrahedron=10
            cell_type = 10;
            break;
         case MshElemType::PRISM:                 // wedge=13
            cell_type = 13;
            break;
         default:
            std::cerr << "***ERROR: NO CORRESPONDING VTK CELL TYPE FOUND. (ELEMENT TYPE=" << ele_type << ")" << std::endl;
      }
      return cell_type;
   }

   void CVTK::InitializeVTU()
   {
      //if (this->useBinary) {
      //======================================================================
      //# Set machine dependent stuff
      //Data type
      if (sizeof(unsigned char)==1)
         type_UChar = CVTK::UInt8;
      else if (sizeof(unsigned char)==2)
         type_UChar = CVTK::UInt16;
      if (sizeof(int)==4)
         type_Int = CVTK::Int32;
      else if (sizeof(int)==8)
         type_Int = CVTK::Int64;
      if (sizeof(unsigned int)==4)
         type_UInt = CVTK::UInt32;
      else if (sizeof(unsigned int)==8)
         type_UInt = CVTK::UInt64;
      if (sizeof(long)==4)
         type_Long = CVTK::Int32;
      else if (sizeof(long)==8)
         type_Long = CVTK::Int64;
      if (sizeof(double)==4)
         type_Double = CVTK::Float32;
      else if (sizeof(double)==8)
         type_Double = CVTK::Float64;
      //
      SIZE_OF_BLOCK_LENGTH_TAG = sizeof(unsigned int);
      //Endian(byte order)
      isLittleEndian = IsLittleEndian();
      //}

      this->isInitialized = true;
   }

   bool CVTK::WriteDataArrayHeader(std::fstream &fin, VTK_XML_DATA_TYPE data_type, const std::string &str_name, int nr_components, const std::string &str_format, long offset)
   {
      std::string str_data_type;
      switch (data_type)
      {
         case CVTK::Int8: str_data_type = "Int8"; break;
         case CVTK::UInt8: str_data_type = "UInt8"; break;
         case CVTK::Int16: str_data_type = "Int16"; break;
         case CVTK::UInt16: str_data_type = "UInt16"; break;
         case CVTK::Int32: str_data_type = "Int32"; break;
         case CVTK::UInt32: str_data_type = "UInt32"; break;
         case CVTK::Int64: str_data_type = "Int64"; break;
         case CVTK::UInt64: str_data_type = "UInt64"; break;
         case CVTK::Float32: str_data_type = "Float32"; break;
         case CVTK::Float64: str_data_type = "Float64"; break;
      }
      fin << "        <DataArray type=\"" << str_data_type << "\"";
      if (str_name!="")
         fin << " Name=\"" << str_name << "\"";
      if (nr_components>1)
         fin << " NumberOfComponents=\"" << nr_components << "\"";
      fin << " format=\"" << str_format << "\"";
      if (useBinary)
         fin << " offset=\"" << offset << "\" /";
      fin << ">" << std::endl;

      return true;
   }

   bool CVTK::WriteDataArrayFooter(std::fstream &fin)
   {
      if (!this->useBinary)
         fin << "        </DataArray>" << std::endl;

      return true;
   }

   bool CVTK::WriteXMLUnstructuredGrid(const std::string &vtkfile, COutput *out, const int time_step_number)
   {
      if (!this->isInitialized)
         this->InitializeVTU();

      //-------------------------------------------------------------------------
      //# Setup file stream
      //-------------------------------------------------------------------------
      std::fstream fin;
      if (this->useBinary)
         fin.open(vtkfile.data(), std::ios::out|std::ios::binary);
      else
         fin.open(vtkfile.data(), std::ios::out);

      if (!fin.good())
      {
         std::cout << "***Warning: Cannot open the output file, " << vtkfile << std::endl;
         return false;
      }

      if (!this->useBinary)
      {
         fin.setf(std::ios::scientific,std::ios::floatfield);
         fin.precision(12);
      }

      //-------------------------------------------------------------------------
      //# Output
      //-------------------------------------------------------------------------
      //  CFEMesh *msh = out->GetMSH();
      CFEMesh *msh = out->getMesh();
      long offset = 0;

      string str_format;
      if (!this->useBinary)
         str_format = "ascii";
      else
         str_format = "appended";
      bool data_out = !useBinary;

      //# Header
      fin << "<?xml version=\"1.0\"?>" << std::endl;
      fin << "<!-- Time step: " << time_step_number << " | Time: " << out->getTime() << " -->" << std::endl;
      fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
      if (!this->useBinary || isLittleEndian)
      {
         fin << " byte_order=\"LittleEndian\"";
      }
      else
      {
         fin << " byte_order=\"BigEndian\"";
      }
      fin << ">" << std::endl;
      //  fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"  << endl;

      //# Unstructured Grid information
      fin << "  <UnstructuredGrid>" << endl;
      fin << "    <Piece NumberOfPoints=\"" << msh->GetNodesNumber(false) << "\" NumberOfCells=\"" << msh->ele_vector.size() << "\">" << std::endl;
      //....................................................................
      // Nodes
      //OK411 CNode *nod = NULL;
      fin << "      <Points>" << std::endl;
      WriteDataArrayHeader(fin, type_Double, "", 3, str_format, offset);
      WriteMeshNodes(fin, data_out, msh, offset);
      WriteDataArrayFooter(fin);
      fin << "      </Points>" << endl;
      //....................................................................
      // Elements
      //OK411 CElem * ele = NULL;
      fin << "      <Cells>" << endl;
      //connectivity
      WriteDataArrayHeader(fin, type_Long, "connectivity", 0, str_format, offset);
      long sum_ele_components = 0;
      WriteMeshElementConnectivity(fin, data_out, msh, offset, sum_ele_components);
      WriteDataArrayFooter(fin);
      //offset
      WriteDataArrayHeader(fin, type_Long, "offsets", 0, str_format, offset);
      WriteMeshElementOffset(fin, data_out, msh, offset);
      WriteDataArrayFooter(fin);
      //type
      WriteDataArrayHeader(fin, type_UChar, "types", 0, str_format, offset);
      WriteMeshElementType(fin, data_out, msh, offset);
      WriteDataArrayFooter(fin);
      fin << "      </Cells>" << endl;
      //....................................................................
      // Nodal values
      if (out->_nod_value_vector.size() > 0)
      {
         fin << "      <PointData Scalars=\"" << out->_nod_value_vector[0] << "\">" << endl;
      }
      else
      {
         fin << "      <PointData Scalars=\"scalars\">" << endl;
      }
      WriteNodalValue(fin, data_out, out, msh, offset);
      fin << "      </PointData>" << endl;

      //======================================================================
      //....................................................................
      // Element values
      fin << "      <CellData>" << endl;
      WriteElementValue(fin, data_out, out, msh, offset);
      fin << "      </CellData>" << endl;
      fin << "    </Piece>" << endl;
      fin << "  </UnstructuredGrid>" << endl;

      //======================================================================
      // Raw data (for binary mode)
      if (useBinary)
      {
         fin << "  <AppendedData encoding=\"raw\">" << endl;
         fin << "    _";

         //Node
         this->WriteMeshNodes(fin, true, msh, offset);
         //Element
         //conncectivity
         this->WriteMeshElementConnectivity(fin, true, msh, offset, sum_ele_components);
         //offset
         this->WriteMeshElementOffset(fin, true, msh, offset);
         //type
         this->WriteMeshElementType(fin, true, msh, offset);
         // Nodal values
         this->WriteNodalValue(fin, true, out, msh, offset);
         // Elemental values
         this->WriteElementValue(fin, true, out, msh, offset);

         fin << endl;
         fin << "  </AppendedData>" << endl;
      }

      fin << "</VTKFile>" << endl;
      fin.close();

      return true;
   }

   bool CVTK::IsLittleEndian()
   {
      int x = 0x00000001;
      if (*(char*)&x) return true;                //am little
      else return false;                          //am big
   }

   template <typename T> void CVTK::write_value_binary(std::fstream &fin, T val)
   {
      fin.write((const char *)&val, sizeof(T));
   }

   bool CVTK::WriteMeshNodes(std::fstream &fin, bool output_data, CFEMesh *msh, long &offset)
   {
      if (output_data)
      {
         CNode *nod = NULL;
         if (!useBinary)
         {
            for (long i=0; i<(long)msh->nod_vector.size(); i++)
            {
               nod = msh->nod_vector[i];
               fin << "          " << nod->X() << " " << nod->Y() << " " << nod->Z() << endl;
            }
         }
         else
         {
                                                  //OK411
            write_value_binary<unsigned int>(fin, sizeof(double)*3*(long)msh->nod_vector.size());
            for (long i=0; i<(long)msh->nod_vector.size(); i++)
            {
               nod = msh->nod_vector[i];
               write_value_binary(fin,  nod->X());
               write_value_binary(fin,  nod->Y());
               write_value_binary(fin,  nod->Z());
            }
         }
      }
      else if (useBinary)
      {
                                                  //OK411
         offset += (long)msh->nod_vector.size()*sizeof(double)*3 + SIZE_OF_BLOCK_LENGTH_TAG;
      }

      return true;
   }

   bool CVTK::WriteMeshElementConnectivity(std::fstream &fin, bool output_data, CFEMesh *msh, long &offset, long &sum_ele_components)
   {
      if (output_data)
      {
         CElem * ele = NULL;
         if (!useBinary)
         {
            for (long i=0; i<(long)msh->ele_vector.size(); i++)
            {
               ele = msh->ele_vector[i];
               fin << "          ";
               for (long j=0; j<ele->GetNodesNumber(false); j++)
                  fin << ele->GetNodeIndex(j) << " ";
               fin << endl;
            }
         }
         else
         {
            write_value_binary<unsigned int>(fin, sizeof(long)*sum_ele_components);
            for (long i=0; i<(long)msh->ele_vector.size(); i++)
            {
               ele = msh->ele_vector[i];
               for (long j=0; j<msh->ele_vector[i]->GetNodesNumber(false); j++)
                  write_value_binary<long>(fin, ele->GetNodeIndex(j));
            }
         }
      }
      else
      {
         if (useBinary)
         {
            sum_ele_components = 0;
            for (long i=0; i<(long)msh->ele_vector.size(); i++)
            {
               sum_ele_components += msh->ele_vector[i]->GetNodesNumber(false);
            }
            offset += sum_ele_components*sizeof(long) + SIZE_OF_BLOCK_LENGTH_TAG;
         }
      }

      return true;
   }

   bool CVTK::WriteMeshElementOffset(std::fstream &fin, bool output_data, CFEMesh *msh, long &offset)
   {
      if (output_data)
      {
         CElem * ele = NULL;

         if (!useBinary)
         {
            fin << "          ";
            long ele_offset = 0;
            for (long i=0; i<(long)msh->ele_vector.size(); i++)
            {
               ele = msh->ele_vector[i];
               ele_offset += ele->GetNodesNumber(false);
               fin << ele_offset << " ";
            }
            fin << endl;
         }
         else
         {
                                                  //OK411
            write_value_binary<unsigned int>(fin, sizeof(long)*(long)msh->ele_vector.size());
            long ele_offset = 0;
            for (long i=0; i<(long)msh->ele_vector.size(); i++)
            {
               ele = msh->ele_vector[i];
               ele_offset += ele->GetNodesNumber(false);
               write_value_binary(fin, ele_offset);
            }
         }
      }
      else
      {
         if (useBinary)
         {
                                                  //OK411
            offset += (long)msh->ele_vector.size()*sizeof(long) + SIZE_OF_BLOCK_LENGTH_TAG;
         }
      }

      return true;
   }

   bool CVTK::WriteMeshElementType(std::fstream &fin, bool output_data, CFEMesh *msh, long &offset)
   {
      if (output_data)
      {
         CElem * ele = NULL;
         if (!useBinary)
         {
            fin << "          ";
            for(long i=0;i<(long)msh->ele_vector.size();i++)
            {
               ele = msh->ele_vector[i];
               fin << (int)this->GetVTKCellType(ele->GetElementType()) << " ";
            }
            fin << endl;
         }
         else
         {
                                                  //OK411
            write_value_binary<unsigned int>(fin, sizeof(unsigned char)*(long)msh->ele_vector.size());
            for(long i=0;i<(long)msh->ele_vector.size();i++)
            {
               ele = msh->ele_vector[i];
               write_value_binary(fin, this->GetVTKCellType(ele->GetElementType()));
            }
         }
      }
      else
      {
         if (useBinary)
         {
                                                  //OK411
            offset += (long)msh->ele_vector.size()*sizeof(unsigned char) + SIZE_OF_BLOCK_LENGTH_TAG;
         }
      }

      return true;
   }


   bool CVTK::WriteNodalValue(std::fstream &fin, bool output_data, COutput *out, CFEMesh *msh, long &offset)
   {
      CRFProcess* m_pcs = NULL;
      std::vector<int> NodeIndex(out->_nod_value_vector.size());
      //  if (out->m_pcs == NULL && out->pcs_type_name.compare("NO_PCS")!=0)
      //  	out->m_pcs = PCSGet(out->pcs_type_name);
      if (out->m_pcs == NULL && out->getProcessType() == NO_PCS)
         out->m_pcs = PCSGet(out->getProcessType());
      if (out->m_pcs != NULL)
         m_pcs = out->m_pcs;

      string str_format;
      if (!this->useBinary)
         str_format = "ascii";
      else
         str_format = "appended";

      bool outNodeVelocity = false;

      //Nodal values
      for (int i = 0; i < (int) out->_nod_value_vector.size(); i++)
      {
         //is velocity
         if (out->_nod_value_vector[i].find("VELOCITY") != string::npos)
         {
            outNodeVelocity = true;
            continue;
         }
         //    if (out->m_pcs == NULL || out->pcs_type_name.compare("NO_PCS")==0)
         if (out->m_pcs == NULL || out->getProcessType() == NO_PCS)
            m_pcs = PCSGet(out->_nod_value_vector[i], true);
         if (!m_pcs)
            continue;
         NodeIndex[i] = m_pcs->GetNodeValueIndex(out->_nod_value_vector[i]);
         if (NodeIndex[i] < 0)
            continue;

         if (!useBinary || !output_data)
         {
            WriteDataArrayHeader(fin, type_Double, out->_nod_value_vector[i],
               0, str_format, offset);
         }

         if (output_data)
         {
            for (size_t j = 0; j < m_pcs->GetPrimaryVNumber(); j++)
            {
               if (out->_nod_value_vector[i].compare(
                  m_pcs->pcs_primary_function_name[j]) == 0)
               {
                  NodeIndex[i]++;
                  break;
               }
            }
            if (!useBinary)
            {
               fin << "          ";
               for (int j = 0; j < msh->GetNodesNumber(false); j++)
               {
                  fin << m_pcs->GetNodeValue(msh->nod_vector[j]->GetIndex(),
                     NodeIndex[i]) << " ";
               }
               fin << endl;
            }
            else
            {
               write_value_binary<unsigned int> (fin, sizeof(double)
                  * msh->GetNodesNumber(false));
               for (int j = 0; j < msh->GetNodesNumber(false); j++)
               {
                  write_value_binary(fin, m_pcs->GetNodeValue(
                     msh->nod_vector[j]->GetIndex(), NodeIndex[i]));
               }
            }
         }
         else
         {
            offset += msh->GetNodesNumber(false) * sizeof(double)
               + SIZE_OF_BLOCK_LENGTH_TAG;
         }

         if (!useBinary || !output_data)
         {
            WriteDataArrayFooter(fin);
         }
      }

      // Nodal velocities
      if (outNodeVelocity)
      {
         unsigned int velocity_id = 0;
         for (int i = 0; i < (int) out->_nod_value_vector.size(); i++)
         {
            if (out->_nod_value_vector[i].find("VELOCITY_X1") != string::npos)
            {
               if (out->m_pcs == NULL)
                  m_pcs = PCSGet(out->_nod_value_vector[i], true);
               velocity_id = 0;
            } else if (out->_nod_value_vector[i].find("VELOCITY_X2")
               != string::npos)
            {
               if (out->m_pcs == NULL)
                  m_pcs = PCSGet(out->_nod_value_vector[i], true);
               velocity_id = 1;
            } else if (out->_nod_value_vector[i].find("VELOCITY1_X")
               != string::npos)
            {
               if (out->m_pcs == NULL)
                  m_pcs = PCSGet(out->_nod_value_vector[i], true);
               velocity_id = 2;
            }
            else
            {
               continue;
            }
            if (!m_pcs)
               continue;

            if (!useBinary || !output_data)
            {
               WriteDataArrayHeader(fin, this->type_Double,
                  velocity_name[velocity_id][3], 3, str_format, offset);
            }
            if (output_data)
            {
               int ix, iy, iz;
               ix = m_pcs->GetNodeValueIndex(velocity_name[velocity_id][0]);
               iy = m_pcs->GetNodeValueIndex(velocity_name[velocity_id][1]);
               iz = m_pcs->GetNodeValueIndex(velocity_name[velocity_id][2]);
               if (!useBinary)
               {
                  fin << "          ";
                  for (int j = 0l; j < msh->GetNodesNumber(false); j++)
                  {
                     fin << m_pcs->GetNodeValue(
                        msh->nod_vector[j]->GetIndex(), ix) << " ";
                     fin << m_pcs->GetNodeValue(
                        msh->nod_vector[j]->GetIndex(), iy) << " ";
                     fin << m_pcs->GetNodeValue(
                        msh->nod_vector[j]->GetIndex(), iz) << " ";
                  }
                  fin << endl;
               }
               else
               {
                  write_value_binary<unsigned int> (fin, sizeof(double)
                     * msh->GetNodesNumber(false) * 3);
                  for (int j = 0l; j < msh->GetNodesNumber(false); j++)
                  {
                     write_value_binary(fin, m_pcs->GetNodeValue(
                        msh->nod_vector[j]->GetIndex(), ix));
                     write_value_binary(fin, m_pcs->GetNodeValue(
                        msh->nod_vector[j]->GetIndex(), iy));
                     write_value_binary(fin, m_pcs->GetNodeValue(
                        msh->nod_vector[j]->GetIndex(), iz));
                  }
               }
            }
            else
            {
               offset += msh->GetNodesNumber(false) * 3* sizeof( double) +SIZE_OF_BLOCK_LENGTH_TAG;
            }

            if (!useBinary || !output_data)
            {
               WriteDataArrayFooter(fin);
            }
         }
      }
      return true;
   }

   bool CVTK::WriteElementValue(std::fstream &fin, bool output_data, COutput *out, CFEMesh *msh, long &offset)
   {
      std::vector<int> ele_value_index_vector(out->getElementValueVector().size());
      if (ele_value_index_vector.size() > 0)      // GetELEValuesIndexVector() should check this!
         out->GetELEValuesIndexVector(ele_value_index_vector);
      CRFProcess* m_pcs = NULL;
      CElem* ele = NULL;

      string str_format;
      if (!this->useBinary)
         str_format = "ascii";
      else
         str_format = "appended";

      bool outEleVelocity = false;

      //Element values
      for (int i = 0; i < (int) ele_value_index_vector.size(); i++)
      {
         if (out->getElementValueVector()[i].find("VELOCITY") != string::npos)
         {
            outEleVelocity = true;
            continue;
         }
         m_pcs = out->GetPCS_ELE(out->getElementValueVector()[i]);

         if (!useBinary || !output_data)
         {
            WriteDataArrayHeader(fin, this->type_Double,
               out->getElementValueVector()[i], 0, str_format, offset);
         }
         if (output_data)
         {
            if (!useBinary)
            {
               fin << "          ";
               for (long j = 0; j < (long) msh->ele_vector.size(); j++)
               {
                  fin << m_pcs->GetElementValue(j, ele_value_index_vector[i])
                     << " ";
               }
               fin << endl;
            }
            else
            {
               write_value_binary<unsigned int> (fin, sizeof(double)
                  * (long) msh->ele_vector.size());//OK411
               for (long j = 0; j < (long) msh->ele_vector.size(); j++)
               {
                  write_value_binary(fin, m_pcs->GetElementValue(j,
                     ele_value_index_vector[i]));
               }
            }
         }
         else
         {
            offset += (long) msh->ele_vector.size() * sizeof(double)
               + SIZE_OF_BLOCK_LENGTH_TAG;        //OK411
         }
         if (!useBinary || !output_data)
         {
            WriteDataArrayFooter(fin);
         }
      }

      //Element veolocity
      if (outEleVelocity)
      {
         if (!useBinary || !output_data)
         {
            WriteDataArrayHeader(fin, this->type_Double, "ELEMENT_VELOCITY", 3,
               str_format, offset);
         }
         if (output_data)
         {
            if (!useBinary)
            {
               fin << "          ";
               static double ele_vel[3] = { 0.0, 0.0, 0.0 };
               for (long i = 0; i < (long) msh->ele_vector.size(); i++)
               {
                  ele_gp_value[i]->getIPvalue_vec(0, ele_vel);
                  fin << ele_vel[0] << " ";
                  fin << ele_vel[1] << " ";
                  fin << ele_vel[2] << " ";
               }
               fin << endl;
            }
            else
            {
               static double ele_vel[3] = { 0.0, 0.0, 0.0 };
               write_value_binary<unsigned int> (fin, sizeof(double) * 3*
                  (long )msh->ele_vector.size()); //OK411
               for(long i=0;i<(long)msh->ele_vector.size();i++)
               {
                  ele_gp_value[i]->getIPvalue_vec(0, ele_vel);
                  write_value_binary(fin, ele_vel[0]);
                  write_value_binary(fin, ele_vel[1]);
                  write_value_binary(fin, ele_vel[2]);
               }
            }
         }
         else
         {
                                                  //OK411
            offset += (long)msh->ele_vector.size()*sizeof(double)*3 +SIZE_OF_BLOCK_LENGTH_TAG;
         }
         if (!useBinary || !output_data)
         {
            WriteDataArrayFooter(fin);
         }

         //    if(out->pcs_type_name.compare("FLUID_MOMENTUM")==0)
         if(out->getProcessType () == FLUID_MOMENTUM)
         {
            if (!useBinary || !output_data)
            {
               WriteDataArrayHeader(fin, this->type_Double,"GLOBAL_VELOCITY", 3, str_format, offset);
            }
            if (output_data)
            {
               CRFProcess* pch_pcs = PCSGet("FLUID_MOMENTUM");
               if (!this->useBinary)
               {
                  fin << "          ";
                  for(long i=0;i<(long)msh->ele_vector.size();i++)
                  {
                     fin << pch_pcs->GetElementValue(i, pch_pcs->GetElementValueIndex("VELOCITY1_X")+1) << " ";
                     fin << pch_pcs->GetElementValue(i, pch_pcs->GetElementValueIndex("VELOCITY1_Y")+1) << " ";
                     fin << pch_pcs->GetElementValue(i, pch_pcs->GetElementValueIndex("VELOCITY1_Z")+1) << " ";
                  }
                  fin << endl;
               }
               else
               {
                                                  //OK411
                  write_value_binary<unsigned int>(fin, sizeof(double)*3*(long)msh->ele_vector.size());
                  for(long i=0;i<(long)msh->ele_vector.size();i++)
                  {
                     write_value_binary(fin, pch_pcs->GetElementValue(i, pch_pcs->GetElementValueIndex("VELOCITY1_X")+1));
                     write_value_binary(fin, pch_pcs->GetElementValue(i, pch_pcs->GetElementValueIndex("VELOCITY1_Y")+1));
                     write_value_binary(fin, pch_pcs->GetElementValue(i, pch_pcs->GetElementValueIndex("VELOCITY1_Z")+1));
                  }
               }
            }
            else
            {
                                                  //OK411
               offset += (long)msh->ele_vector.size()*sizeof(double)*3 + SIZE_OF_BLOCK_LENGTH_TAG;
            }
            if (!useBinary || !output_data)
            {
               WriteDataArrayFooter(fin);
            }
         }
      }
      //Material information
      if(mmp_vector.size() > 1)
      {
         if (!useBinary || !output_data)
         {
            WriteDataArrayHeader(fin, this->type_Int, "MatGroup", 0, str_format, offset);
         }
         if (output_data)
         {
            if (!this->useBinary)
            {
               fin << "          ";
               for(long i=0;i<(long)msh->ele_vector.size();i++)
               {
                  ele = msh->ele_vector[i];
                  fin << ele->GetPatchIndex() << " ";
               }
               fin << endl;
            }
            else
            {
                                                  //OK411
               write_value_binary<unsigned int>(fin, sizeof(int)*(long)msh->ele_vector.size());
               for (long i=0;i<(long)msh->ele_vector.size();i++)
                  write_value_binary(fin, msh->ele_vector[i]->GetPatchIndex());
            }
         }
         else
         {
                                                  //OK411
            offset += (long)msh->ele_vector.size()*sizeof(int) + SIZE_OF_BLOCK_LENGTH_TAG;
         }
         if (!useBinary || !output_data)
         {
            WriteDataArrayFooter(fin);
         }
      }
      return true;
   }
