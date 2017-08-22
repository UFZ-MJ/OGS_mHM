/**************************************************************************
MSHLib - Object:
Task:
Programing:
08/2005 WW/OK Encapsulation from rf_ele_msh
last modified
**************************************************************************/
#ifndef msh_lib_INC
#define msh_lib_INC

#include "msh_mesh.h"

#ifdef USE_TOKENBUF
#include "tokenbuf.h"
#endif

using Mesh_Group::CFEMesh;
using Mesh_Group::CElem;
using Mesh_Group::CNode;

extern std::vector<Mesh_Group::CFEMesh*> fem_msh_vector;
extern CFEMesh* FEMGet(const std::string &msh_name);
                                                  //OK
extern void MSHCreateNOD2ELERelations(Mesh_Group::CFEMesh*);

extern CFEMesh* FEMRead(const std::string& , GEOLIB::GEOObjects* geo_obj = NULL, std::string* unique_name = NULL);

extern void MSHWrite(std::string);
extern void CompleteMesh();                       //WW
extern bool CompleteMesh(std::string);            //OK
extern void FEMDeleteAll();
extern void MSHCalcMinMaxMidCoordinates();        //OK
extern double msh_x_min,msh_x_max;                //OK
extern double msh_y_min,msh_y_max;                //OK
extern double msh_z_min,msh_z_max;                //OK
extern double msh_x_mid,msh_y_mid,msh_z_mid;      //OK
// Might be removed
void Read_RFI(std::istream& msh_file, CFEMesh* m_msh);
extern void MSHAssignMATGroup2Elements(std::string);
extern void MSHCreateQuadsFromPLY(CGLPolyline*,int);
//OK411 extern void MSHCreatePrismsFromTriangles();
extern void MSHCreateNodes();
extern void MSHDeleteDoubleElements(int);
extern long MSHMarkDoubleElementsType(int);
extern void RFIWriteTecplot();
extern void MSHWriteTecplot();
extern void MSHAssignMATGroup2LineElements();
extern void MSHAssignMATGroup2TrisElements(std::string);
extern void MSHAssignMATGroup2QuadElements();
extern void MSHAssignMATGroup2TetsElements();
extern void MSHAssignMATGroup2PrisElements();
extern void MSHAssignMATGroup2PrisElementsNew();
extern void MSH2MATPris();
extern void MSHAssignMATGroup2HexsElements();
extern void MSHDestroy();
extern void MSHDelete(std::string);
extern void DATWriteRFIFile(const char *file_name);
extern void DATWriteParticleFile(int);            // PCH
extern void MSHWriteVOL2TEC(std::string);         //OK
//KR extern bool msh_file_binary; //OK
extern void GMSH2MSH(const char*,CFEMesh*);
                                                  //TK
extern void Mesh_Single_Surface(std::string surface_name, const char *file_name_const_char);
//extern void Select_Nodes_Elements_by_TINFile(const char *file_name_const_char);
extern void Clear_Selected_Nodes_Elements();
extern void GMSH2TIN(const char *file_name_const_char);
extern void MSHLayerWriteTecplot();               //OK

                                                  //OK
extern CFEMesh* MSHGet(const std::string &mat_type_name);
                                                  //OK
extern CFEMesh* MSHGet(const std::string &pcs_type_name,const std::string &mat_type_name);
extern CFEMesh* MSHGetGEO(std::string);           //OK
extern int MSHSetMaxMMPGroups();                  //OK
extern bool MSHTestMATGroups();                   //OK
#ifdef RFW_FRACTURE
extern bool MSHGetCommonNodes(CElem*, CElem*, vector<CNode*>&);
extern void MSHSetFractureElements(void);
extern void MSHResetFractureElements(void);
                                                  //RFW
extern long MSHWhatElemIsPointIn(double x, double y, long index);
#endif
extern void MSHDefineMobile(CRFProcess*);         //OK411
extern void MSHMoveNODUcFlow (CRFProcess*);       //OK411
extern long* MSHGetNodesClose(long*,CGLPolyline*);//OK411

//extern bool IsPointInSurface(Surface*,CGLPoint*); //OK411

extern long* GetPointsIn(Surface*,long*);         //OK411
                                                  //OK411
extern void GEOGetNodesInMaterialDomain(CFEMesh*, const int, std::vector<long>&, bool);
extern void SetRFIPointsClose(CGLLine*);          //OK411
                                                  //OK411
extern void MSHGetNodesClose(std::vector<long>&,CGLPoint*);

#endif
