/**************************************************************************
GeoSys-Project
Task: class implementation
Programing:
01/2005 OK Implementation
last modified:
**************************************************************************/
#ifndef gs_project_INC
#define gs_project_INC

#include "Configure.h"
// C++ STL
#include <string>
#include <vector>

// GEOLib
// FEMLib
//----------------------------------------------------------------
// ToDo CGSProject
class CGSProject
{
  private:
  public:
    std::string path;
    std::string base;
    std::string type;
};

extern std::vector<CGSProject*>gsp_vector;
extern void GSPRead(std::string);
extern void GSPWrite(); 
extern void GSPAddMember(std::string);
//extern void GSPRemoveMemberstd::(string);
extern void GSPWriteData();
extern CGSProject* GSPGetMember(std::string);
extern void GSPAddMemberNew(std::string path_base_orig,std::string path_base_copy,std::string type);
extern bool GSPSimulatorReady();
extern bool GSPReadData(std::string); //OK

extern std::string g_gsp_path;
extern std::string g_gsp_base;
#define GSP_FILE_EXTENSION ".gsp"

#endif
