/**************************************************************************
GeoSys - Object:
Task:
Programing:
01/2005 OK/TK Implementation
last modified:
**************************************************************************/

#include "gs_project.h"
// C++ STL
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

// FEM-Makros
#include "makros.h"
extern std::string GetLineFromFile1(std::ifstream*);
#include "files0.h"
// GeoSys-GEOLib
#include "geo_lib.h"
// MSHLib
#include "msh_lib.h"
// GeoSys-FEMLib
#include "rf_out_new.h"
#include "rf_bc_new.h"
#include "rf_st_new.h"
#include "rf_ic_new.h"
#include "rf_mmp_new.h"
#include "rf_pcs.h"
#include "rf_tim_new.h"
#include "rf_fct.h"
#include "rf_out_new.h"
#include "rf_mfp_new.h"
#include "rf_msp_new.h"
#include "rfmat_cp.h"
extern bool RFDOpen(std:: string file_name_base);
#include "problem.h"

// FileIO
#include "OGSIOVer4.h"

//==========================================================================
std::vector<CGSProject*>gsp_vector;
std::string g_gsp_path;
std::string g_gsp_base;
std::string gsp_path_base;

/**************************************************************************
GeoSys-Method:
Task:
Programing:
01/2005 OK Implementation based on ExecuteProjectFile (TK)
last modification:
**************************************************************************/

// 05/2010 TF not needed for compilation
//void GSPRead(string gsp_path_base_this)
//{
//  CGSProject* m_gs_project = NULL;
//  //----------------------------------------------------------------------
//  // File handling
//  string gsp_this;
//  gsp_this = gsp_path_base_this + ".gsp";
//  FILE* gsp_file_this = NULL;
//  gsp_file_this = fopen(gsp_this.c_str(),"r");
//  if(gsp_file_this)
//    fclose(gsp_file_this);
//  gsp_file_this = fopen(gsp_this.c_str(),"r");
//  char char_line[MAX_ZEILE];
//  //
//  ifstream gsp_file(gsp_this.data(),ios::in);
//  if (!gsp_file.good()){
//    cout << "Warning: no project data *.gsp file is missing" << endl;
//    return;
//  }
//  gsp_file.clear();
//  gsp_file.seekg(0,ios::beg); // rewind the file
//  //----------------------------------------------------------------------
//  // File analysis
//  int pos;
//  pos = (int)gsp_path_base_this.find_last_of('\\');
//  g_gsp_path = gsp_path_base_this.substr(0,pos);
//  //----------------------------------------------------------------------
//
//  GEOLIB::GEOObjects *geo_obj (NULL);
//  std::string unique_fname;
//  unique_fname += ".gli";
//
////OK  while(!gsp_file.eof()){
//  bool its_true = true;
//  while(its_true){
//    //gsp_file >> gsp_member_file;
//    fgets(char_line,MAX_ZEILE,gsp_file_this);
//    string gsp_member_file = char_line;
//    if(gsp_member_file.find("STOP")!=string::npos)
//      return;
//    pos = (int)gsp_member_file.find_last_of('.');
//    string gsp_member_file_base = gsp_member_file.substr(0,pos);
//    string gsp_member_file_extension = gsp_member_file.substr(pos+1,string::npos);
//    string gsp_member_path_base = g_gsp_path + "\\" + gsp_member_file_base;
//    //....................................................................
//    if(gsp_member_file_extension.find("gli")!=string::npos){
//
//    	// old data structures
//    	GEOLIB_Read_GeoLib(gsp_member_path_base);
//
//    	// new data structures
//    	if (!geo_obj) geo_obj = new GEOLIB::GEOObjects ();
//    	unique_fname = gsp_member_file_extension;
//    	FileIO::readGLIFileV4 (unique_fname, geo_obj);
//
//    	m_gs_project = new CGSProject;
//	  m_gs_project->path = g_gsp_path; //?
//	  m_gs_project->base = gsp_member_file_base;
//	  m_gs_project->type = gsp_member_file_extension; //?
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("rfi")==0){
//      MSHOpen(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("pcs")==0){
//      PCSRead(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//      pcs_created = true;
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("tim")==0){
//      TIMRead(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("out")==0){
//      OUTRead(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("num")==0){
//      NUMRead(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("ic")==0){
//      ICRead(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("bc")==0){
//
//    	// 05/2010 TF
////    	BCRead(gsp_member_path_base);
//        BCRead(gsp_member_path_base, geo_obj, unique_fname);
//
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("st")==0){
//      STRead(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("mfp")==0){
//      MFPRead(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("msp")==0){
//      MSPRead(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("mmp")==0){
//      MMPRead(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("cp")==0){
//      CPRead(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("rfd")==0){
//      RFDOpen(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    else if(gsp_member_file_extension.compare("fct")==0){
//      FCTRead(gsp_member_path_base);
//	  m_gs_project = new CGSProject;
//	  m_gs_project->base = gsp_member_file_base;
//	  gsp_vector.push_back(m_gs_project);
//    }
//    //....................................................................
//    gsp_file.ignore(MAX_ZEILE,'\n');
//  }
//  //----------------------------------------------------------------------
//}

/**************************************************************************
GeoSys-Method:
Task: Liest Projectmember-Liste und schreibt die Komponenten in das
      GSP-File
Programing:
01/2004 TK project related file handling
01/2005 OK V4 implementation
last modification:
**************************************************************************/
void GSPWrite()
{
  int i;
  std::string sub_line;
  std::string line_string;
  CGSProject* m_gsp = NULL;
  //========================================================================
  // File handling
  std::string gsp_this;
  gsp_this = g_gsp_path + g_gsp_base + GSP_FILE_EXTENSION;
  FILE* gsp_file_this = NULL;
  gsp_file_this = fopen(gsp_this.c_str(),"r");
  if(gsp_file_this)
    fclose(gsp_file_this);
  gsp_file_this = fopen(gsp_this.c_str(),"w");
  //
  std::fstream gsp_file (gsp_this.data(),std::ios::trunc|std::ios::out);
  //gsp_file.close();
  //gsp_file.open(gsp_this.data(),ios::trunc|ios::out);
  if (!gsp_file.good())
    return;
  gsp_file.seekg(0L,std::ios::beg);
  //========================================================================
  gsp_file << "#PROJECT_MEMBERS" << std::endl;
  //========================================================================
  // GSP vector
  std::string line;
  int gsp_vector_size =(int)gsp_vector.size();
  for(i=0;i<gsp_vector_size;i++){
    m_gsp = gsp_vector[i];
    gsp_file << m_gsp->base << "." << m_gsp->type << std::endl;
  }
  gsp_file << "#STOP";
  gsp_file.close();
}

/**************************************************************************
GeoSys-Method:
Task:
Programing:
01/2005 OK GSP function
last modification:
**************************************************************************/
void GSPRemoveMember(std::string file_type)
{
  int i;
  CGSProject* m_gsp = NULL;
  int gsp_vector_size = (int)gsp_vector.size();
  for(i=0;i<gsp_vector_size;i++){
    m_gsp = gsp_vector[i];
    if(m_gsp->type.compare(file_type)==0){
      delete m_gsp;
      gsp_vector.erase((gsp_vector.begin()+i));
      return;
    }
  }
}

/**************************************************************************
GeoSys-Method:
Task: Liest ?ergbenen File-pfad und erg?zt das Projektfile
Programing:
11/2003 TK project related file handling
01/2005 OK GSP function
last modification:
**************************************************************************/
void GSPAddMember(std::string base_plus_type)
{
  CGSProject* m_gsp = NULL;
  // String analysis
  size_t pos = base_plus_type.find_last_of('.');
  std::string base = base_plus_type.substr(0,pos);
  std::string type = base_plus_type.substr(pos+1,std::string::npos);
  // Check GSP member and remove existing GSP type
  GSPRemoveMember(type);
  // Add GSP member
  m_gsp = new CGSProject;
  m_gsp->path = g_gsp_path;
  m_gsp->base = base;
  m_gsp->type = type;
  gsp_vector.push_back(m_gsp);
}

/**************************************************************************
GeoSys-Method:
Task:
Programing:
01/2005 OK Implementation
last modification:
**************************************************************************/
void GSPAddMemberNew(std::string path_base_orig,std::string path_base_copy,std::string type)
{
  // Copy file to GSP directory
  std::string path_base_type_orig = path_base_orig + "." + type;
  std::string path_base_type_copy = path_base_copy + "." + type;
  char input_text[MAX_ZEILE];
  FILE *gsp_member_file_orig = NULL;
  FILE *gsp_member_file_copy = NULL;
  gsp_member_file_orig = fopen(path_base_type_orig.c_str(),"rt");
  gsp_member_file_copy = fopen(path_base_type_copy.c_str(),"w+t");
  if(!gsp_member_file_orig)//OK4104
    return;
  while(!feof(gsp_member_file_orig)){
    fgets(input_text,MAX_ZEILE,gsp_member_file_orig);
    FilePrintString(gsp_member_file_copy,input_text);
  }
  fclose(gsp_member_file_orig);
  fclose(gsp_member_file_copy);
  // Check GSP member and remove existing GSP type
  GSPRemoveMember(type);
  // Add GSP member
  int pos;
  char a = '\\';
  pos = (int)path_base_copy.find_last_of(a);
  std::string path = path_base_copy.substr(0,pos);
  std::string base = path_base_copy.substr(pos+1,std::string::npos);
  CGSProject* m_gsp = NULL;
  m_gsp = new CGSProject;
  m_gsp->path = path;
  m_gsp->base = base;
  m_gsp->type = type;
  gsp_vector.push_back(m_gsp);
}

/**************************************************************************
GeoSys-Method:
Task:
Programing:
01/2005 OK Implementation
05/2005 OK MSHWrite
last modification:
**************************************************************************/
void GSPWriteData()
{
  int i;
  std::string path_base;
  CGSProject* m_gsp = NULL;
  int gsp_vector_size =(int)gsp_vector.size();
  for(i=0;i<gsp_vector_size;i++){
    m_gsp = gsp_vector[i];
    if(m_gsp->type.compare("gli")==0){
      path_base = g_gsp_path + m_gsp->base;
      GEOWrite(path_base);
    }
    if(m_gsp->type.compare("msh")==0){
      path_base = g_gsp_path + m_gsp->base;
      MSHWrite(path_base);
    }
    else if(m_gsp->type.compare("pct")==0){     // PCH
      path_base = g_gsp_path + m_gsp->base;
      DATWriteParticleFile(aktueller_zeitschritt);
    }
    else if(m_gsp->type.compare("pcs")==0){
      path_base = g_gsp_path + m_gsp->base;
      PCSWrite(path_base.data());
    }
    else if(m_gsp->type.compare("num")==0){
      path_base = g_gsp_path + m_gsp->base;
      NUMWrite(path_base.data());
    }
    else if(m_gsp->type.compare("tim")==0){
      path_base = g_gsp_path + m_gsp->base;
      TIMWrite(path_base.data());
    }
    else if(m_gsp->type.compare("out")==0){
      path_base = g_gsp_path + m_gsp->base;
      OUTWrite(path_base.data());
    }
    else if(m_gsp->type.compare("ic")==0){
      path_base = g_gsp_path + m_gsp->base;
      ICWrite(path_base.data());
    }
    else if(m_gsp->type.compare("bc")==0){
      path_base = g_gsp_path + m_gsp->base;
      BCWrite(path_base.data());
    }
    else if(m_gsp->type.compare("st")==0){
      path_base = g_gsp_path + m_gsp->base;
      STWrite(path_base.data());
    }
    else if(m_gsp->type.compare("mfp")==0){
      path_base = g_gsp_path + m_gsp->base;
      MFPWrite(path_base.data());
    }
    else if(m_gsp->type.compare("msp")==0){
      path_base = g_gsp_path + m_gsp->base;
      MSPWrite(path_base.data());
    }
    else if(m_gsp->type.compare("mmp")==0){
      path_base = g_gsp_path + m_gsp->base;
      MMPWrite(path_base.data());
    }
    else if(m_gsp->type.compare("mcp")==0){
      path_base = g_gsp_path + m_gsp->base;
      CPWrite(path_base.data(),0);
    }
    else if(m_gsp->type.compare("fct")==0){
      path_base = g_gsp_path + m_gsp->base;
      FCTWrite(path_base.data());
    }
  }
}

/**************************************************************************
GeoSys-Method:
Task:
Programing:
01/2005 OK GSP function
last modification:
**************************************************************************/
CGSProject* GSPGetMember(std::string file_type)
{
  int i;
  CGSProject* m_gsp = NULL;
  int gsp_vector_size = (int)gsp_vector.size();
  for(i=0;i<gsp_vector_size;i++){
    m_gsp = gsp_vector[i];
    if(m_gsp->type.compare(file_type)==0){
      return m_gsp;
    }
  }
  return NULL;
}

/**************************************************************************
GeoSys-Method:
Task:
Programing:
01/2005 OK GSP function
last modification:
**************************************************************************/
bool GSPSimulatorReady()
{
  std::string problem_type = PCSProblemType();
  if(!GSPGetMember("gli"))
    return false;
  if(!GSPGetMember("rfi")&&(!GSPGetMember("msh")))
    return false;
  if(!GSPGetMember("pcs"))
    return false;
  if(!GSPGetMember("num"))
    return false;
  if(!GSPGetMember("tim"))
    return false;
  if(!GSPGetMember("ic"))
    return false;
  if(!GSPGetMember("bc"))
    return false;
  if(!GSPGetMember("mmp"))
    return false;
  //......................................................................
  if(problem_type.find("LIQUID_FLOW")!=std::string::npos){
    if(!GSPGetMember("mfp"))
      return false;
  }
  //......................................................................
  if(problem_type.find("DEFORMATION")!=std::string::npos){
    if(!GSPGetMember("msp"))
      return false;
  }
  //......................................................................
  if(problem_type.find("HEAT")!=std::string::npos){
    if(!GSPGetMember("msp"))
      return false;
  }
  //......................................................................
  if(problem_type.find("MASS")!=std::string::npos){
    if(!GSPGetMember("mcp"))
      return false;
  }
  //......................................................................
  return true;
}


/**************************************************************************
GeoSys-Method:
06/2009 OK Implementation
**************************************************************************/
// 05/2010 TF not needed for compilation
//bool GSPReadData(string gsp_member_path_base)
//{
//  //----------------------------------------------------------------------
//  GEOLIB_Read_GeoLib(gsp_member_path_base);
//  //......................................................................
//  FEMRead(gsp_member_path_base);
//  //......................................................................
//  PCSRead(gsp_member_path_base);
//  TIMRead(gsp_member_path_base);
//  OUTRead(gsp_member_path_base);
//  NUMRead(gsp_member_path_base);
//  //......................................................................
//  ICRead(gsp_member_path_base);
//  BCRead(gsp_member_path_base);
//  STRead(gsp_member_path_base);
//  //......................................................................
//  MFPRead(gsp_member_path_base);
//  MSPRead(gsp_member_path_base);
//  MMPRead(gsp_member_path_base);
//  CPRead(gsp_member_path_base);
//  //......................................................................
//  RFDOpen(gsp_member_path_base);
//  FCTRead(gsp_member_path_base);
//  //----------------------------------------------------------------------
//  return true;
//}
