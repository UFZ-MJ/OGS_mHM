/**************************************************************************
 MSHLib - Object:
 Task:
 Programing:
 08/2005 OK Encapsulated from mshlib
 **************************************************************************/

#include <string>
#include <vector>
#include <list>

// GEOLib
#include "geo_sfc.h"

// MSHLib
#include "msh_lib.h"
// PCSLib
#include "rf_mmp_new.h"

/**************************************************************************
 GeoSysGUI-Function: MSHAssignMATGroup2Elements
 Programing:
 05/2003 OK Implementation
 last modified:
 **************************************************************************/
void MSHAssignMATGroup2LineElements(void)
{
   /*OK411
    long j;
    int k,l;
    long *element_nodes;
    CGLPolyline *m_polyline= NULL;
    CGLPoint m_point;
    long *ply_nodes;
    long no_ply_nodes;
    //-----------------------------------------------------------------------
    CMediumProperties *m_mat_mp = NULL;
    int m;
   int no_mat_mp = (int)mmp_vector.size();
   for(m=0;m<no_mat_mp;m++){
   m_mat_mp = mmp_vector[m];
   //.....................................................................
   if(m_mat_mp->geo_type_name.compare("POLYLINE")==0){
   m_polyline = GEOGetPLYByName(m_mat_mp->geo_name);//CC 10/05
   ply_nodes = MSHGetNodesClose(&no_ply_nodes,m_polyline);//CC 10/05
   for(j=0;j<ElListSize();j++){
   if((ElGetElementType(j)==1)){ // Lines
   element_nodes = ElGetElementNodes(j);
   for(k=0;k<no_ply_nodes;k++) {
   if(element_nodes[0]==ply_nodes[k]){
   for(l=0;l<no_ply_nodes;l++) {
   if(element_nodes[1]==ply_nodes[l]){
   ElSetElementGroupNumber(j,m_mat_mp->number);
   }
   }
   }
   }
   }
   }
   }
   }
   */
}


/**************************************************************************
 GeoSysGUI-Function:
 Programing:
 05/2003 OK Implementation
 08/2005 OK MSH
 10/2005 OK MAT-GEO
 last modified:
 **************************************************************************/
void MSHAssignMATGroup2TrisElements(std::string msh_name)
{
   // Tests
   CFEMesh* m_msh (FEMGet(msh_name));
   if (!m_msh)
   {
      std::cout << "Warning: MSHAssignMATGroup2TrisElements: no MSH data: "
         << msh_name << std::endl;
      return;
   }
   //-----------------------------------------------------------------------
   // Initialize MAT groups
   for (size_t j = 0; j < m_msh->ele_vector.size(); j++)
   {
      if (m_msh->ele_vector[j]->GetElementType() == MshElemType::TRIANGLE)
         m_msh->ele_vector[j]->SetPatchIndex(-1);
   }
   //-----------------------------------------------------------------------
   // TF m_sfc->PointInSurface(&m_point) returns always false
   //	int m;
   //	int k;
   //	CGLPoint m_pnt;
   //	CMediumProperties *m_mmp = NULL;
   //	Surface *m_sfc = NULL;
   //	double* xyz;
   //	for (m = 0; m < (int) mmp_vector.size(); m++) {
   //		m_mmp = mmp_vector[m];
   //		if (m_mmp->geo_type_name.compare("SURFACE") == 0) {
   //			for (j = 0; j < (long) m_msh->ele_vector.size(); j++) {
   //				m_ele = m_msh->ele_vector[j];
   //				if (m_ele->GetElementType() == 4) {
   //					xyz = m_ele->GetGravityCenter();
   //					m_pnt.x = xyz[0];
   //					m_pnt.y = xyz[1];
   //					m_pnt.z = xyz[2];
   //					for (k = 0; k < (int) m_mmp->geo_name_vector.size(); k++) {
   //						m_sfc = GEOGetSFCByName(m_mmp->geo_name_vector[k]);
   //						if (m_sfc) {
   //							if (m_sfc->PointInSurface(&m_pnt)) {
   //								m_ele->SetPatchIndex(m_mmp->number);
   //							}
   //						}
   //					}
   //				}
   //			}
   //		}
   //	}
}


void MSHAssignMATGroup2QuadElements(void)
{
   /*OK411
    long j;
    int nn,k;
    double x,y,z;
    long *element_nodes;
    Surface *m_surface = NULL;
    CGLPoint m_point;
    //-----------------------------------------------------------------------
    CMediumProperties *m_mat_mp = NULL;
    int m;
    int no_mat_mp = (int)mmp_vector.size();
   for(m=0;m<no_mat_mp;m++){
   m_mat_mp = mmp_vector[m];
   if(m_mat_mp->geo_type_name.compare("SURFACE")==0){
   for(j=0;j<ElListSize();j++){
   if(ElGetElementType(j)==2){
   element_nodes = ElGetElementNodes(j);
   nn = ElNumberOfNodes[ElGetElementType(j)-1];
   x=0.0; y=0.0; z=0.0;
   for(k=0;k<nn;k++) {
   x += GetNodeX(element_nodes[k]);
   y += GetNodeY(element_nodes[k]);
   z += GetNodeZ(element_nodes[k]);
   }
   x /= double(nn);
   y /= double(nn);
   z /= double(nn);
   m_point.x = x;
   m_point.y = y;
   m_point.z = z;
   m_surface = GEOGetSFCByName(m_mat_mp->geo_name);//CC 10/05
   if(m_surface){
   if(IsPointInSurface(m_surface,&m_point)){//CC 10/05
   ElSetElementGroupNumber(j,m_mat_mp->number);
   }
   }
   }
   }
   }
   }
   */
}


void MSHAssignMATGroup2TetsElements(void)
{
   /*OK411
    long j;
    int nn,k;
    double x,y,z;
    long *element_nodes;
    CGLVolume *m_volume = NULL;
    CGLPoint m_point;
    //-----------------------------------------------------------------------
    CMediumProperties *m_mat_mp = NULL;
    int m;
    int no_mat_mp = (int)mmp_vector.size();
   for(m=0;m<no_mat_mp;m++){
   m_mat_mp = mmp_vector[m];
   if(m_mat_mp->geo_type_name.compare("VOLUME")==0){
   for(j=0;j<ElListSize();j++){
   if(ElGetElementType(j)==5){
   element_nodes = ElGetElementNodes(j);
   nn = ElNumberOfNodes[ElGetElementType(j)-1];
   x=0.0; y=0.0; z=0.0;
   for(k=0;k<nn;k++) {
   x += GetNodeX(element_nodes[k]);
   y += GetNodeY(element_nodes[k]);
   z += GetNodeZ(element_nodes[k]);
   }
   x /= double(nn);
   y /= double(nn);
   z /= double(nn);
   m_point.x = x;
   m_point.y = y;
   m_point.z = z;
   vector<CGLVolume*>::iterator p_vol;//CC
   p_vol = volume_vector.begin();
   while(p_vol!=volume_vector.end()) {
   m_volume = *p_vol;
   if(m_volume->PointInVolume(&m_point,0)){
   ElSetElementGroupNumber(j,m_volume->mat_group);
   }
   ++p_vol;
   }
   }
   }
   }
   }
   */
}


void MSHAssignMATGroup2PrisElements(void)
{
   /*OK411
    long j;
    int nn,k;
    double x,y,z;
    long *element_nodes;
    CGLVolume *m_volume = NULL;
    CGLPoint m_point;
    //-----------------------------------------------------------------------
    CMediumProperties *m_mat_mp = NULL;
    int m;
    int no_mat_mp = (int)mmp_vector.size();
   for(m=0;m<no_mat_mp;m++){
   m_mat_mp = mmp_vector[m];
   if(m_mat_mp->geo_type_name.compare("VOLUME")==0){
   for(j=0;j<ElListSize();j++){
   if(ElGetElementType(j)==6){
   element_nodes = ElGetElementNodes(j);
   nn = ElNumberOfNodes[ElGetElementType(j)-1];
   x=0.0; y=0.0; z=0.0;
   for(k=0;k<nn;k++) {
   x += GetNodeX(element_nodes[k]);
   y += GetNodeY(element_nodes[k]);
   z += GetNodeZ(element_nodes[k]);
   }
   x /= double(nn);
   y /= double(nn);
   z /= double(nn);
   m_point.x = x;
   m_point.y = y;
   m_point.z = z;
   vector<CGLVolume*>::iterator p_vol;//CC
   p_vol = volume_vector.begin();
   while(p_vol!=volume_vector.end()) {
   m_volume = *p_vol;
   if(m_volume->PointInVolume(&m_point,0)){
   ElSetElementGroupNumber(j,m_volume->mat_group);
   }
   ++p_vol;
   }
   }
   }
   }
   }
   */
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void MSHAssignMATGroup2PrisElementsNew(void)
{
   /*OK411
    long i,j;
    CGLVolume *m_vol = NULL;
    vector<CGLVolume*>::iterator p_vol;
    p_vol = volume_vector.begin();
    long *nodes = NULL;
    double x,y,z;
    CGLPoint m_point;
    int mmp_group = 0;
    CMediumProperties* m_mmp = NULL;
    long counter = 0;
   //----------------------------------------------------------------------
   // Initializations
   for(i=0;i<ElementListLength;i++){
   if(ElGetElementType(i)==6){
   ElSetElementGroupNumber(i,-1);
   }
   }
   //----------------------------------------------------------------------
   while(p_vol!=volume_vector.end()){
   m_vol = *p_vol;
   m_vol->mat_group = mmp_group;
   //--------------------------------------------------------------------
   for(i=0;i<ElementListLength;i++){
   if(ElGetElementType(i)==6){
   //.....................................................................
   // Element center point
   nodes = ElGetElementNodes(i);
   x=0.0; y=0.0; z=0.0;
   for(j=0;j<6;j++) {
   x += GetNodeX(nodes[j]);
   y += GetNodeY(nodes[j]);
   z += GetNodeZ(nodes[j]);
   }
   x /= double(6);
   y /= double(6);
   z /= double(6);
   m_point.x = x;
   m_point.y = y;
   m_point.z = z;
   if(m_vol->PointInVolume(&m_point,0)){
   ElSetElementGroupNumber(i,m_vol->mat_group);
   counter++;
   }
   }
   else{
   ElSetElementGroupNumber(i,-1);
   }
   }
   //--------------------------------------------------------------------
   m_mmp = mmp_vector[m_vol->mat_group];
   m_mmp->number = m_vol->mat_group;
   if(m_mmp)
   m_mmp->WriteTecplot("ToDo");
   //--------------------------------------------------------------------
   mmp_group++;
   ++p_vol;
   }
   */
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void MSH2MATPris(void)
{
   /*OK411
    long j;
    int nn,k,m;
    double x,y,z;
    long *element_nodes;
    CGLVolume *m_vol = NULL;
    CGLPoint m_point;
    CMediumProperties* m_mmp = NULL;
    //-----------------------------------------------------------------------
    vector<CGLVolume*>::iterator p_vol;//CC
    int mmp_vector_size = (int)mmp_vector.size();
   for(m=0;m<mmp_vector_size;m++){
   m_mmp = mmp_vector[m];
   p_vol = volume_vector.begin();//CC
   while(p_vol!=volume_vector.end()) {//CC
   m_vol = *p_vol;
   if(m_mmp->geo_name.compare(m_vol->name)==0){
   m_vol->mat_group_name = m_mmp->geo_name;
   m_vol->mat_group = m;
   break;
   }
   ++p_vol;
   }
   }
   //======================================================================
   // test file
   string mat_test_file_name = "MAT_test_file.txt";
   fstream mat_test_file (mat_test_file_name.data(),ios::trunc|ios::out);
   mat_test_file.setf(ios::scientific,ios::floatfield);
   mat_test_file.precision(12);
   mat_test_file.seekg(0L,ios::beg);
   //-----------------------------------------------------------------------
   for(j=0;j<ElListSize();j++){
   if(fmod((double)j,100.)<1e-3)
   mat_test_file << j << endl;
   if(ElGetElementType(j)==6){
   ElSetElementGroupNumber(j,-1);
   element_nodes = ElGetElementNodes(j);
   nn = ElNumberOfNodes[ElGetElementType(j)-1];
   x=0.0; y=0.0; z=0.0;
   for(k=0;k<nn;k++) {
   x += GetNodeX(element_nodes[k]);
   y += GetNodeY(element_nodes[k]);
   z += GetNodeZ(element_nodes[k]);
   }
   x /= double(nn);
   y /= double(nn);
   z /= double(nn);
   m_point.x = x;
   m_point.y = y;
   m_point.z = z;
   vector<CGLVolume*>::iterator p_vol;//CC
   p_vol = volume_vector.begin();//CC
   while(p_vol!=volume_vector.end()) {//CC
   m_vol = *p_vol;
   if(m_vol->PointInVolume(&m_point,0)){
   ElSetElementGroupNumber(j,m_vol->mat_group);
   }
   ++p_vol;
   }
   }
   }
   //======================================================================
   // Write failed elements
   string elements_failed = "elements_failed.txt";
   fstream elements_failed_file (elements_failed.data(),ios::trunc|ios::out);
   elements_failed_file.seekg(0L,ios::beg);
   long element_mat_group;
   for(j=0;j<ElListSize();j++){
   element_mat_group = ElGetElementGroupNumber(j);
   if(element_mat_group<0){
   elements_failed_file << "Element number: " << j \
   << ", Element type: " << ElGetElementType(j) << endl;
   }
   }
   elements_failed_file.close();
   //======================================================================
   // Write Tecplot files
   for(m=0;m<mmp_vector_size;m++){
   m_mmp = mmp_vector[m];
   m_mmp->WriteTecplot("ToDo");
   }
   */
}


void MSHAssignMATGroup2HexsElements(void)
{
   /*OK411
    long j;
    int nn,k;
    double x,y,z;
    long *element_nodes;
    CGLVolume *m_volume = NULL;
    CGLPoint m_point;
    //-----------------------------------------------------------------------
    CMediumProperties *m_mat_mp = NULL;
    int m;
    int no_mat_mp = (int)mmp_vector.size();
   for(m=0;m<no_mat_mp;m++){
   m_mat_mp = mmp_vector[m];
   if(m_mat_mp->geo_type_name.compare("VOLUME")==0){
   for(j=0;j<ElListSize();j++){
   if(ElGetElementType(j)==3){
   element_nodes = ElGetElementNodes(j);
   nn = ElNumberOfNodes[ElGetElementType(j)-1];
   x=0.0; y=0.0; z=0.0;
   for(k=0;k<nn;k++) {
   x += GetNodeX(element_nodes[k]);
   y += GetNodeY(element_nodes[k]);
   z += GetNodeZ(element_nodes[k]);
   }
   x /= double(nn);
   y /= double(nn);
   z /= double(nn);
   m_point.x = x;
   m_point.y = y;
   m_point.z = z;
   vector<CGLVolume*>::iterator p_vol;//CC
   p_vol = volume_vector.begin();
   while(p_vol!=volume_vector.end()) {
   m_volume = *p_vol;
   if(m_volume->PointInVolume(&m_point,0)){
   ElSetElementGroupNumber(j,m_volume->mat_group);
   }
   ++p_vol;
   }
   }
   }
   }
   }
   */
}

void MSHAssignMATGroup2Elements(std::string msh_name)
{
   MSHAssignMATGroup2LineElements();
   MSHAssignMATGroup2TrisElements(msh_name);
   MSHAssignMATGroup2QuadElements();
   MSHAssignMATGroup2TetsElements();
   MSHAssignMATGroup2PrisElements();
   MSHAssignMATGroup2HexsElements();
   MSHSetMaxMMPGroups();
}


/**************************************************************************
 MSHLib-Method:
 01/2006 OK Implementation
 06/2009 OK Bug fix
 **************************************************************************/
int MSHSetMaxMMPGroups()
{
   int i;
   long j;
   CFEMesh* m_msh = NULL;
   //----------------------------------------------------------------------
   int msh_max_mmp_groups;
   for (i = 0; i < (int) fem_msh_vector.size(); i++)
   {
      m_msh = fem_msh_vector[i];
      m_msh->max_mmp_groups = 0;
      msh_max_mmp_groups = 0;
      for (j = 0; j < (long) m_msh->ele_vector.size(); j++)
      {
         if ((m_msh->ele_vector[j]->GetPatchIndex() + 1)
            > msh_max_mmp_groups)                 //OK
            msh_max_mmp_groups++;
      }
      m_msh->max_mmp_groups = msh_max_mmp_groups;
   }
   //----------------------------------------------------------------------
   int g_msh_max_mmp_groups = 0;
   for (i = 0; i < (int) fem_msh_vector.size(); i++)
   {
      if (m_msh->max_mmp_groups > g_msh_max_mmp_groups)
         g_msh_max_mmp_groups++;
   }
   //----------------------------------------------------------------------
   return g_msh_max_mmp_groups;
}


/**************************************************************************
 MSHLib-Method:
 07/2007 OK Implementation
 **************************************************************************/
bool MSHTestMATGroups()
{
   int g_max_mmp_groups = MSHSetMaxMMPGroups();
   if (g_max_mmp_groups > (int) mmp_vector.size())
   {
      std::cout << "Error: not enough MMP data";
      return false;                               //abort();
   }
   return true;
}
