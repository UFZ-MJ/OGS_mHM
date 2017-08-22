/**************************************************************************
MSHLib - Object:
Task:
Programing:
08/2005 OK Encapsulated from mshlib
**************************************************************************/

#include "math.h"
// C++
#include <string>
#include <vector>

// MSHLib
#include "msh_lib.h"
#include "msh_elem.h"
#include "mathlib.h"
// PCSLib
#include "rf_pcs.h"
#include "gs_project.h"
#include "rf_mmp_new.h"                           //OK

/**************************************************************************
MshLib-Method:
Task: Create QUAD elements from line elements
Programing:
02/2004 OK Implementation
**************************************************************************/
void MSHCreateQuadsFromLine(CGLLine *m_line)
{
   m_line = m_line;
   /*OK411
     int i,k;
     int col=0;
     long j;
     long **matrix = NULL;
     double x,y,x1,y1;
     long msh_node_number;
     double dist;
     double eps = 1e-3;
     int no_msh_nodes = m_line->no_msh_nodes;
     //---------------------------------------------------------------------
   // Node columns
   ///matrix = new long*[number_of_msh_nodes];
   matrix = (long**) Malloc(no_msh_nodes*sizeof(long));
   for(i=0;i<no_msh_nodes;i++) {
   msh_node_number = m_line->msh_nodes[i];
   x1 = GetNodeX(msh_node_number);
   y1 = GetNodeY(msh_node_number);
   col=0;
   matrix[i] = NULL;
   for(j=0;j<NodeListLength;j++) {
   x = GetNodeX(j); // ToDo NodeNumber[j]
   y = GetNodeY(j);
   dist = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
   if(dist<eps) {
   col++;
   matrix[i] = (long*) Realloc(matrix[i],col*sizeof(long));
   matrix[i][col-1] = j; // ToDoNodeNumber[j];
   }
   }
   }
   //-----------------------------------------------------------------------
   // Create Quads
   long *nodes = NULL;
   for(i=0;i<no_msh_nodes-1;i++) {
   for(k=0;k<col-1;k++) {
   nodes = (long*) Malloc(sizeof(long)*ElNumberOfNodes[1]);
   nodes[0] = matrix[i][k];
   nodes[1] = matrix[i+1][k];
   nodes[2] = matrix[i+1][k+1];
   nodes[3] = matrix[i][k+1];
   //OK      ELECreateTopology(2,-1,0,ElementListLength);
   ElSetElementNodes(ElementListLength-1,nodes);
   ElSetElementGroupNumber(ElementListLength-1,m_line->mat_group);
   anz_2D++;
   }
   }
   //-----------------------------------------------------------------------
   // Memory
   for(i=0;i<no_msh_nodes-1;i++) {
   matrix[i] = (long*) Free(matrix[i]);
   }
   matrix = (long**) Free(matrix);
   */
}


/**************************************************************************
MshLib-Method:
Task: Create QUAD elements from polylines
Programing:
02/2004 OK Implementation
**************************************************************************/
void MSHCreateQuadsFromPLY(CGLPolyline *m_polyline, int msh_type)
{
   switch (msh_type)
   {
      case 1:
         CGLLine *m_line = NULL;
         std::vector<CGLLine*>::iterator pl;      //CC

         // create polyline lines, if necessary
         if (m_polyline->getLineVector().size() == 0)
            m_polyline->ComputeLines();
         pl = m_polyline->getLineVector().begin();//CC
                                                  //CC
         while (pl != m_polyline->getLineVector().end())
         {
            m_line = *pl;
            m_line->mat_group = m_polyline->getMatGroup();
            if (m_line->no_msh_nodes == 0)
               SetRFIPointsClose(m_line);         //CC 10/05
            MSHCreateQuadsFromLine(m_line);
            pl++;
         }
         break;
   }
}


/**************************************************************************
MSHLib-Method:
Task:
Programing:
07/2005 MB Implementation
11/2005 WW/MB element class
03/2006 OK mat_group
**************************************************************************/
void CFEMesh::CreatePriELEFromTri(int no_layer,double layer_thickness)
{
   no_layer = no_layer;
   layer_thickness = layer_thickness;
   int j;
   int hempel, hampel;
   long i;
   int k;
   long size;
   long no_tri_elements = (long)ele_vector.size();
   long no_tri_nodes = (long)nod_vector.size();
   Mesh_Group::CElem* m_tri_ele = NULL;
   Mesh_Group::CElem* m_ele = NULL;
   //----------------------------------------------------------------------
   // Create MSH
   MSHDelete("PRIS_from_TRI");
   CFEMesh* m_msh_pri (new CFEMesh(_geo_obj, _geo_name));
   m_msh_pri->pcs_name = "PRIS_from_TRI";
   m_msh_pri->setElementType (MshElemType::PRISM);
   //m_msh_pri->no_msh_layer = no_layer;
   m_msh_pri->setNumberOfMeshLayers (no_layer);
   //----------------------------------------------------------------------
   // Create Prism elements
   size = (no_layer + 1) * no_tri_nodes;
   m_msh_pri->nod_vector.resize(size);
   for(j=0;j<size;j++)
      m_msh_pri->nod_vector[j] = NULL;
   for(j=0;j<no_layer;j++)
   {
      for(i=0;i<no_tri_elements;i++)
      {
         //Elements
         m_tri_ele = ele_vector[i];
         m_ele = new  Mesh_Group::CElem;
                                                  //OK
         m_ele->SetPatchIndex((int)mmp_vector.size()+j);
         m_ele->SetElementType(MshElemType::PRISM);
         m_ele->nnodes = 6;
         m_ele->nodes_index.resize(m_ele->nnodes);
         //Set indices
         m_ele->nodes_index[0] = m_tri_ele->GetNodeIndex(0)+ j*no_tri_nodes;
         m_ele->nodes_index[1] = m_tri_ele->GetNodeIndex(1) + j*no_tri_nodes;
         m_ele->nodes_index[2] = m_tri_ele->GetNodeIndex(2) + j*no_tri_nodes;
         m_ele->nodes_index[3] = m_ele->GetNodeIndex(0)+ no_tri_nodes;
         m_ele->nodes_index[4] = m_ele->GetNodeIndex(1)+ no_tri_nodes;
         m_ele->nodes_index[5] = m_ele->GetNodeIndex(2)+ no_tri_nodes;
         //Nodes
         hempel = 0;
         hampel = 0;
         for(k=0;k<m_ele->nnodes;k++)
         {
            if(m_msh_pri->nod_vector[m_ele->GetNodeIndex(k)]==NULL)
            {
               m_ele->nodes[k] = new CNode(m_ele->GetNodeIndex(k));
               m_msh_pri->nod_vector[m_ele->GetNodeIndex(k)] = m_ele->nodes[k];
               if(k>2)
               {
                  hempel = 3;
                  hampel = 1;
               }
               m_ele->nodes[k]->SetX(nod_vector[m_tri_ele->GetNodeIndex(k-hempel)]->X());
               m_ele->nodes[k]->SetY(nod_vector[m_tri_ele->GetNodeIndex(k-hempel)]->Y());
               m_ele->nodes[k]->SetZ(nod_vector[m_tri_ele->GetNodeIndex(k-hempel)]->Z() - (j+hampel)*layer_thickness);
            }
            else
            {
               m_ele->nodes[k] = m_msh_pri->nod_vector[m_ele->GetNodeIndex(k)];
            }
         }                                        //end for m_ele->nnodes
         m_msh_pri->ele_vector.push_back(m_ele);
      }                                           //end for no_tri_elements
   }                                              //end for no_layers
   //----------------------------------------------------------------------
   if(m_msh_pri->ele_vector.size()>0)
   {
      m_msh_pri->ConstructGrid();
      fem_msh_vector.push_back(m_msh_pri);
   }
   else
      delete m_msh_pri;
}


/**************************************************************************
MshLib-Method: CreateLines
Task: Create linear 1-D elements from polylines
Programing:
05/2005 OK Implementation
**************************************************************************/
void CFEMesh::CreateLineELEFromPLY(CGLPolyline *m_ply)
{
   m_ply = m_ply;
   //WWW
#ifdef TODO
   long i;
   double dx,dy,dz;
   CGLPoint* m_pnt = NULL;
   //----------------------------------------------------------------------
   // Create LINE elements
   FiniteElement::CElement* m_quad_ele = NULL;
   for(i=0;i<nr;i++)
   {
      m_ele = new FiniteElement::CElement;
      m_ele->GetPatchIndex() = (int)mmp_vector.size();
      m_ele->type_name = "line";
      m_ele->ElementType = 1;
      m_ele->nnodes = 2;
      m_ele->nodes = new long[2];
      m_ele->nodes_index[0] = i;
      m_ele->nodes_index[1] = i+1;
      ele_vector.push_back(m_ele);
   }
   //----------------------------------------------------------------------
   // Create LINE nodes
   CMSHNodes* m_nod = NULL;
   m_pnt = m_pnt->Get(m_ply->point_vector[0]->gli_point_id);
   m_ply->point_vector[0]->x = m_pnt->x;
   m_ply->point_vector[0]->y = m_pnt->y;
   m_ply->point_vector[0]->z = m_pnt->z;
   m_pnt = m_pnt->Get(m_ply->point_vector[1]->gli_point_id);
   m_ply->point_vector[1]->x = m_pnt->x;
   m_ply->point_vector[1]->y = m_pnt->y;
   m_ply->point_vector[1]->z = m_pnt->z;
   dx = (m_ply->point_vector[1]->x - m_ply->point_vector[0]->x) / (double)nr;
   dy = (m_ply->point_vector[1]->y - m_ply->point_vector[0]->y) / (double)nr;
   dz = (m_ply->point_vector[1]->z - m_ply->point_vector[0]->z) / (double)nr;
   for(i=0;i<nr+1;i++)
   {
      m_nod = new CMSHNodes();
      //m_nod->nodenumber = i + j*no_quad_elements;
      m_nod->x = m_ply->point_vector[0]->x + i*dx;
      m_nod->y = m_ply->point_vector[0]->y + i*dy;
      m_nod->z = m_ply->point_vector[0]->z + i*dz;
      nod_vector.push_back(m_nod);
   }
   //----------------------------------------------------------------------
#endif
}


/**************************************************************************
MSHLib-Method:
Task:
Programing:
04/2005 OK Implementation
11/2005 MB ELE
**************************************************************************/
void CFEMesh::CreateLineELEFromQuad(int m_numberofprismlayers,double m_thicknessofprismlayer,int m_iMATGroup)
{
   m_iMATGroup = m_iMATGroup;
   int j;
   long i, k;
   double x, y, z;
   long size;
   long no_quad_elements = (long)ele_vector.size();
   //long no_quad_nodes = (long)nod_vector.size();
   double* center = NULL;
   Mesh_Group::CElem* m_quad_ele = NULL;
   Mesh_Group::CElem* m_ele = NULL;

   //----------------------------------------------------------------------
   // Create MSH
   MSHDelete("LINE_from_QUAD");
   CFEMesh* m_msh_line(new CFEMesh(_geo_obj, _geo_name));
   m_msh_line->pcs_name = "LINE_from_QUAD";
   m_msh_line->setElementType (MshElemType::LINE);
   m_msh_line->setNumberOfMeshLayers (m_numberofprismlayers);
   //----------------------------------------------------------------------
   // Create LINE elements
   size = (m_numberofprismlayers+1) * no_quad_elements;
   m_msh_line->nod_vector.resize(size);

   for(j=0;j<size;j++)
   {
      m_msh_line->nod_vector[j]=NULL;
   }

   //for(j=0;j<m_numberofprismlayers;j++){
   //  for(i=0;i<no_quad_elements;i++){
   for(i=0;i<no_quad_elements;i++)                //Elements
   {
      for(j=0;j<m_numberofprismlayers;j++)
      {
         //Elements
         m_quad_ele = ele_vector[i];
         m_ele = new Mesh_Group::CElem;
         //m_ele->GetPatchIndex() = m_msh_line->mat_group;
         m_ele->SetIndex((i*m_numberofprismlayers) + j);
         m_ele->SetElementType(MshElemType::LINE);
         m_ele->nnodes = 2;
         m_ele->nodes_index.resize(m_ele->nnodes);
         //Set indices
         m_ele->nodes_index[0] = j + i * (m_numberofprismlayers+1);
         m_ele->nodes_index[1] = m_ele->nodes_index[0]+ 1;
         m_msh_line->ele_vector.push_back(m_ele);
         //
         center = m_quad_ele->GetGravityCenter();
         x = center[0];
         y = center[1];
         z = center[2];
         //Nodes
         for(k=0;k<m_ele->nnodes;k++)
         {
                                                  //if new node
            if(m_msh_line->nod_vector[m_ele->GetNodeIndex(k)]==NULL)
            {
               m_ele->nodes[k] = new CNode(m_ele->GetNodeIndex(k));
               m_msh_line->nod_vector[m_ele->GetNodeIndex(k)] = m_ele->nodes[k];
               //Set Coordinates
               m_ele->nodes[k]->SetX(x);
               m_ele->nodes[k]->SetY(y);
               z = z + ((k+j) * m_thicknessofprismlayer);
               m_ele->nodes[k]->SetZ(z);
            }                                     // end if new node
            else
            {
               m_ele->nodes[k] = m_msh_line->nod_vector[m_ele->GetNodeIndex(k)];
            }
         }                                        // end for nnodes
      }                                           // end for no_layers
   }                                              // end for quad elements
   if(m_msh_line->ele_vector.size()>0)
      fem_msh_vector.push_back(m_msh_line);
   else
      delete m_msh_line;
}


/**************************************************************************
MSHLib-Method:
Task:
Programing:
05/2005 OK Implementation based on 2D_2 by KM
**************************************************************************/
void CFEMesh::CreateQuadELEFromSFC(Surface*m_sfc)
{
   m_sfc = m_sfc;
   ///WWW
#ifdef TODO

   long ir,is;
   double x[4],y[4],z[4];
   CGLPolyline* m_ply = NULL;
   //----------------------------------------------------------------------
   //
                                                  //CC
   vector<CGLPolyline*>::iterator p_sfc_ply = m_sfc->polyline_of_surface_vector.begin();
                                                  //CC
   while(p_sfc_ply!=m_sfc->polyline_of_surface_vector.end())
   {
      m_ply = *p_sfc_ply;
      if(m_ply)
      {
         x[0] = m_ply->point_vector[0]->x;
         x[1] = m_ply->point_vector[1]->x;
         x[2] = m_ply->point_vector[2]->x;
         x[3] = m_ply->point_vector[3]->x;
         y[0] = m_ply->point_vector[0]->y;
         y[1] = m_ply->point_vector[1]->y;
         y[2] = m_ply->point_vector[2]->y;
         y[3] = m_ply->point_vector[3]->y;
         z[0] = m_ply->point_vector[0]->z;
         z[1] = m_ply->point_vector[1]->z;
         z[2] = m_ply->point_vector[2]->z;
         z[3] = m_ply->point_vector[3]->z;
      }
      ++p_sfc_ply;
   }
   //----------------------------------------------------------------------
   // Nodes
   int k;
   double N4[4];
   double u[2];
   CMSHNodes* m_nod = NULL;
   for(is=0;is<ns+1;is++)
   {
      u[1] = 1. - 2.*(double)is/(double)ns;
      for(ir=0;ir<(nr+1);ir++)
      {
         m_nod = new CMSHNodes();
         u[0] = 1. - 2.*(double)ir/(double)nr;
         ShapeFunctionQuad(N4,u);
         for(k=0;k<4;k++)
         {
            m_nod->x += N4[k]*x[k];
            m_nod->y += N4[k]*y[k];
            m_nod->z += N4[k]*z[k];
         }
         nod_vector.push_back(m_nod);
      }
   }
   //----------------------------------------------------------------------
   // Elements
   FiniteElement::CElement* m_ele = NULL;
   for(is=0;is<ns;is++)
   {
      for(ir=0;ir<nr;ir++)
      {
         m_ele = new FiniteElement::CElement();
         m_ele->type_name = "quad";
         m_ele->nnodes = 4;
                                                  //m_sfc->mat_group;
         m_ele->GetPatchIndex() = (int)mmp_vector.size();
         m_ele->nodes = new long[4];
         m_ele->nodes_index[0] = is*(nr+1) + ir;
         m_ele->nodes_index[1] = is*(nr+1) + ir+1;
         m_ele->nodes_index[2] = (is+1)*(nr+1) + ir+1;
         m_ele->nodes_index[3] = (is+1)*(nr+1) + ir;
         m_ele->nodes_xyz = new double[12];
         for(k=0;k<m_ele->nnodes;k++)
         {
            m_ele->nodes_xyz[k]                 = nod_vector[m_ele->nodes_index[k]]->x;
            m_ele->nodes_xyz[k+m_ele->nnodes]   = nod_vector[m_ele->nodes_index[k]]->y;
            m_ele->nodes_xyz[k+2*m_ele->nnodes] = nod_vector[m_ele->nodes_index[k]]->z;
         }
         ele_vector.push_back(m_ele);
      }
   }
   //----------------------------------------------------------------------
#endif
}


/**************************************************************************
MSHLib-Method:
Task:
Programing:
04/2005 OK Implementation
**************************************************************************/
void CFEMesh::AppendLineELE()
{
   //WWW
#ifdef TODO
   int j;
   long i;
   //----------------------------------------------------------------------
   // Create LINE elements
   long no_line_elements = (long)ele_vector.size();
   CFEMesh* m_msh_quad = NULL;
   m_msh_quad = FEMGet("QUADfromSFC");
   long no_line_elements_base = (long)m_msh_quad->ele_vector.size();
   for(j=0;j<no_layer;j++)
   {
      for(i=0;i<no_line_elements_base;i++)
      {
         m_ele = new FiniteElement::CElement;
         m_ele->GetPatchIndex() = mat_group;
         m_ele->type_name = "line";
         m_ele->ElementType = 1;
         m_ele->nnodes = 2;
         m_ele->nodes = new long[2];
         m_ele->nodes_index[0] = no_line_elements + i + j*no_line_elements_base;
         m_ele->nodes_index[1] = m_ele->nodes_index[0]+ no_line_elements_base;
         ele_vector.push_back(m_ele);
      }
   }
   //----------------------------------------------------------------------
   // Create LINE nodes
   CMSHNodes* m_nod = NULL;
   double dz = z_min/(double)no_layer;
   long no_line_nodes = (long)nod_vector.size();
   long no_line_nodes_base = (long)m_msh_quad->ele_vector.size();
   long i_start = no_line_nodes - no_line_nodes_base;
   for(j=1;j<no_layer+1;j++)
   {
      for(i=i_start;i<no_line_nodes;i++)
      {
         m_nod = new CMSHNodes();
         m_nod->x = nod_vector[i]->x;
         m_nod->y = nod_vector[i]->y;
         m_nod->z = nod_vector[i]->z + j*dz;
         nod_vector.push_back(m_nod);
      }
   }
   //----------------------------------------------------------------------
#endif
}


/**************************************************************************
MshLib-Method: CreateLines
Task: Create linear 1-D elements from triangulated meshes
Programing:
01/2004 OK Implementation
04/2005 OK MSH method for case 1
04/2007 OK OO-FEM
07/2007 NW MSH method for case 3,4
**************************************************************************/
void CFEMesh::CreateLineELEFromPLY(CGLPolyline *m_polyline,int type,CFEMesh*m_msh_ply)
{
   CGLLine *m_line=NULL;
   std::list<CGLLine*>::const_iterator pl;
   //  int hits;
   long i,j,k;
   //WW  long *nodes_unsorted = NULL;
   //WW  double *node_distances = NULL;
   std::list<CGLLine*>msh_line_list;
   std::list<CGLLine*>::iterator pl1;
   std::list<CGLLine*>::iterator pl2;
   std::list<CGLLine*>::iterator pl3;
   long m_point11,m_point12,m_point21,m_point22;
   CGLLine *m_line1,*m_line2;
   //  bool hitP0,hitP1,hitP2;
   //  double v1[3],v2[3];
   //  double angle;
   //WW  double eps_angle = 1.;
   //OK
   std::vector<long>nodes_vector;
   CNode* m_nod = NULL;
   CElem* m_ele = NULL;
   CEdge* m_edg = NULL;
   vec<long>node_indeces(3);
   long no_elements;
   long no_nodes;
   double m_nod_x,m_nod_y,m_nod_z;
   std::vector<long>elements_vector;
   vec<CNode*>edge_nodes(3);
   vec<CEdge*>ele_edges_vector(15);
   std::vector<long>ele_vector_at_ply;
   std::vector<long>nod_vector_at_ply;

   //======================================================================
   switch(type)
   {
      //--------------------------------------------------------------------
      case 0:                                     // simply sort by distance
         GetNODOnPLY(m_polyline,nod_vector_at_ply);
         //------------------------------------------------------------------
         // Create nodes
         if(m_msh_ply)
         {
            for(i=0;i<(long)nod_vector_at_ply.size();i++)
            {
               m_nod = nod_vector[i];
               no_nodes = (long)m_msh_ply->nod_vector.size();
               m_nod_x = nod_vector[nod_vector_at_ply[i]]->X();
               m_nod_y = nod_vector[nod_vector_at_ply[i]]->Y();
               m_nod_z = nod_vector[nod_vector_at_ply[i]]->Z();
               m_nod = new CNode(no_nodes,m_nod_x,m_nod_y,m_nod_z);
               m_msh_ply->nod_vector.push_back(m_nod);
            }
         }
         //------------------------------------------------------------------
         // Create elements
         for(i=0;i<(long)nod_vector_at_ply.size();i++)
         {
            no_elements = (long)m_msh_ply->ele_vector.size();
            m_ele = new CElem(no_elements);
            m_ele->setElementProperties(MshElemType::LINE);
            m_ele->nodes_index[0] = i;
            m_ele->nodes_index[1] = i+1;
            m_ele->SetPatchIndex((int)mmp_vector.size());
            m_msh_ply->ele_vector.push_back(m_ele);
         }
         /*OK
               // 1 - sort by distance
               nodes_unsorted = m_polyline->MSHGetNodesClose(&no_nodes);
               pt1[0] = m_polyline->point_vector[0]->x;
               pt1[1] = m_polyline->point_vector[0]->y;
               pt1[2] = m_polyline->point_vector[0]->z;
               node_distances = new double[no_nodes];
               for(i=0;i<no_nodes;i++) {
                 pt2[0] = GetNodeX(nodes_unsorted[i]);
                 pt2[1] = GetNodeY(nodes_unsorted[i]);
                 pt2[2] = GetNodeZ(nodes_unsorted[i]);
         node_distances[i] = MCalcDistancePointToPoint(pt1,pt2);
         }
         nodes = TOLSortNodes1(nodes_unsorted,node_distances,no_nodes);
         delete [] node_distances;
         // 2 - create polyline MSH nodes
         m_polyline->msh_nodes_vector.clear();
         for(i=0;i<no_nodes;i++){
         m_polyline->msh_nodes_vector.push_back(nodes[i]);
         }
         // 3 - create line elements
         m_polyline->MSHCreateLines();
         */
         break;
         //--------------------------------------------------------------------
      case 1:                                     // based on lines
         /*OK
               // 0 - create polyline lines, if necessary
               if(m_polyline->line_list.size()==0)
                 m_polyline->ComputeLines(m_polyline);
               // 1 - get MSH nodes close to polyline->lines
               pl = m_polyline->line_list.begin();
               while(pl!=m_polyline->line_list.end()) {
                 m_line = *pl;
                 m_line->epsilon = m_polyline->epsilon;
                 if(ElListSize()>0)
                   m_line->SetRFIPointsClose();
         if((long)ele_vector.size()>0)
         SetLINPointsClose(m_line);
         pl++;
         }
         // 2 - create line elements for polyline->lines
         pl = m_polyline->line_list.begin();
         while(pl!=m_polyline->line_list.end()) {
         m_line = *pl;
         m_line->mat_group = m_polyline->mat_group;
         if(ElListSize()>0)
         m_line->CreateMSHLines();
         if((long)ele_vector.size()>0)
         m_msh_ply->CreateLineELEFromLIN(m_line);
         pl++;
         }
         // 3 - create new MSH nodes for polyline->lines
         pl = m_polyline->line_list.begin();
         while(pl!=m_polyline->line_list.end()) {
         m_line = *pl;
         for(i=0;i<m_line->no_msh_nodes-1;i++){
         m_nod = new CMSHNodes();
         m_nod->origin_rfi_node_number = m_line->msh_nodes[i];
         m_nod->x = nod_vector[m_line->msh_nodes[i]]->x;
         m_nod->y = nod_vector[m_line->msh_nodes[i]]->y;
         m_nod->z = nod_vector[m_line->msh_nodes[i]]->z;
         m_msh_ply->nod_vector.push_back(m_nod);
         }
         pl++;
         }
         pl--;
         m_line = *pl; // last node
         if(m_line){
         m_nod = new CMSHNodes();
         m_nod->origin_rfi_node_number = m_line->msh_nodes[m_line->no_msh_nodes-1];
         m_nod->x = nod_vector[m_line->msh_nodes[m_line->no_msh_nodes-1]]->x;
         m_nod->y = nod_vector[m_line->msh_nodes[m_line->no_msh_nodes-1]]->y;
         m_nod->z = nod_vector[m_line->msh_nodes[m_line->no_msh_nodes-1]]->z;
         m_msh_ply->nod_vector.push_back(m_nod);
         }
         // 4 - renumber elements
         for(i=0;i<(long)m_msh_ply->ele_vector.size();i++){
         m_ele = m_msh_ply->ele_vector[i];
         m_ele->nodes_index[0] = i;
         m_ele->nodes_index[1] = i+1;
         }
         */
         break;
         //--------------------------------------------------------------------
      case 2:                                     // based on triangles
         msh_line_list.clear();
         m_polyline->getLineVector().clear();
         ele_vector_at_ply.clear();
         GetELEOnPLY(m_polyline,ele_vector_at_ply);
         for(i=0;i<(long)ele_vector_at_ply.size();i++)
         {
            std::cout << ele_vector_at_ply[i] << std::endl;
            m_ele = ele_vector[ele_vector_at_ply[i]];
            m_line = new CGLLine();
            m_line->m_point1 = new CGLPoint();
            m_line->m_point2 = new CGLPoint();
            m_ele->GetEdges(ele_edges_vector);
            for(j=0;j<(int)m_ele->GetEdgesNumber();j++)
            {
               m_edg = ele_edges_vector[j];
               if(m_edg->GetMark())
               {
                  m_edg->GetNodes(edge_nodes);
                  m_line->m_point1->id = edge_nodes[0]->GetIndex();
                  m_line->m_point2->id = edge_nodes[1]->GetIndex();
               }
            }
            m_line->gli_line_id = i;
            msh_line_list.push_back(m_line);
         }
         //------------------------------------------------------------------
         //OK m_polyline->line_list.clear();
         pl1 = msh_line_list.begin();
         while(pl1!=msh_line_list.end())
         {
            m_line = *pl1;
            //OK m_line->m_point1->x = GetNodeX(m_line->m_point1->id);
            m_line->m_point1->x = nod_vector[m_line->m_point1->id]->X();
            //OK m_line->m_point1->y = GetNodeY(m_line->m_point1->id);
            m_line->m_point1->y = nod_vector[m_line->m_point1->id]->Y();
            //OK m_line->m_point2->x = GetNodeX(m_line->m_point2->id);
            m_line->m_point2->x = nod_vector[m_line->m_point2->id]->X();
            //OK m_line->m_point2->y = GetNodeY(m_line->m_point2->id);
            m_line->m_point2->y = nod_vector[m_line->m_point2->id]->Y();
            m_polyline->getLineVector().push_back(m_line);
            ++pl1;
         }
         //------------------------------------------------------------------
         // Remove double elements
         pl1 = msh_line_list.begin();
         while(pl1!=msh_line_list.end())
         {
            m_line1 = *pl1;
            m_point11 = m_line1->m_point1->id;
            m_point12 = m_line1->m_point2->id;
            pl2=msh_line_list.begin();
            while(pl2!=msh_line_list.end())
            {
               m_line2 = *pl2;
               if (m_line1->gli_line_id==m_line2->gli_line_id)
               {
                  pl2++;
                  continue;
               }
               m_point21 = m_line2->m_point1->id;
               m_point22 = m_line2->m_point2->id;
               if( ( (m_point11==m_point21)&&(m_point12==m_point22) ) \
                  ||( (m_point11==m_point22)&&(m_point12==m_point21) ) ) \
               { \
                  pl3 = msh_line_list.erase(pl2);
                  pl2 = pl3;
               }
               else
                  pl2++;
            }
            pl1++;
         }
         //------------------------------------------------------------------
         // Create elements
         pl1=msh_line_list.begin();
         while(pl1!=msh_line_list.end())
         {
            m_line = *pl1;
            m_line->mat_group = m_polyline->getMatGroup();
            m_line->no_msh_nodes = 2;
            m_line->msh_nodes = new long[2];
            m_line->msh_nodes[0] = m_line->m_point1->id;
            m_line->msh_nodes[1] = m_line->m_point2->id;
            //m_line->CreateMSHLines();
            //OK    ELECreateTopology(1,-1,0,ElementListLength);
            if(m_msh_ply)
               no_elements = (long)m_msh_ply->ele_vector.size();
            else
               no_elements = (long)ele_vector.size();
            m_ele = new CElem(no_elements);
            m_ele->setElementProperties(MshElemType::LINE);
            //OK        ElSetElementNodes(ElementListLength-1,m_line->msh_nodes);
            m_ele->nodes_index[0] = m_line->m_point1->id;
            m_ele->nodes_index[1] = m_line->m_point2->id;
            //OK        ElSetElementGroupNumber(ElementListLength-1,m_line->mat_group);
            m_ele->SetPatchIndex(m_line->mat_group);
            if(m_msh_ply)
               m_msh_ply->ele_vector.push_back(m_ele);
            else
               ele_vector.push_back(m_ele);
            //OK        anz_1D++;
            pl1++;
         }
         //------------------------------------------------------------------
         // Create nodes
         if(m_msh_ply)
         {
            for(i=0;i<(long)m_msh_ply->ele_vector.size();i++)
            {
               m_ele = m_msh_ply->ele_vector[i];
               for(k=0;k<m_ele->GetNodesNumber(false);k++)
               {
                  if(m_msh_ply->NodeExists(m_ele->nodes_index[k]))
                     continue;
                  no_nodes = (long)m_msh_ply->nod_vector.size();
                  m_nod_x = nod_vector[m_ele->nodes_index[k]]->X();
                  m_nod_y = nod_vector[m_ele->nodes_index[k]]->Y();
                  m_nod_z = nod_vector[m_ele->nodes_index[k]]->Z();
                  m_nod = new CNode(no_nodes,m_nod_x,m_nod_y,m_nod_z);
                  //m_msh_local->nod_vector[j]->SetEquationIndex(j);
                  //m_msh_local->Eqs2Global_NodeIndex[j] = m_msh_local->nod_vector[j]->GetIndex();
                  m_msh_ply->nod_vector.push_back(m_nod);
                                                  // renumber
                  m_ele->nodes_index[k] = no_nodes;
               }
            }
         }
         break;
         //====================================================================
      case 3:                                     // based on triangles / create new mesh (NW)
      {
         m_polyline->getLineVector().clear();
         ele_vector_at_ply.clear();

         //------------------------------------------------------------------
         // Search elements that locate on the polyline
         GetELEOnPLY(m_polyline, ele_vector_at_ply);

         //// DEBUG CODE: highlight collected elements
         //for(i=0;i<(long)ele_vector_at_ply.size();i++) {
         //  ele_vector[ele_vector_at_ply[i]]->selected=1;
         //}

         //------------------------------------------------------------------
         // Check marked edges: really necessary or not
         this->CheckMarkedEdgesOnPolyLine(m_polyline, ele_vector_at_ply);

         //------------------------------------------------------------------
         // Create line elements
         this->CreateLineElementsFromMarkedEdges(m_msh_ply, ele_vector_at_ply);

         //------------------------------------------------------------------
         // Create nodes into new mesh
         if(m_msh_ply)
         {
            for (i=0; i<(long)m_msh_ply->ele_vector.size(); i++)
            {
               m_ele = m_msh_ply->ele_vector[i];
               m_ele->nodes.resize(m_ele->GetNodesNumber(false));

               for (j=0; j<m_ele->GetNodesNumber(false); j++)
               {
                  if (m_msh_ply->NodeExists(m_ele->nodes_index[j]))
                  {
                     for (k=0; k<(long)m_msh_ply->nod_vector.size(); k++)
                     {
                        if(m_msh_ply->nod_vector[k]->GetIndex() == (size_t)m_ele->nodes_index[j])
                        {
                           m_ele->nodes[j] = m_msh_ply->nod_vector[k];
                           m_ele->nodes_index[j] = k;
                           break;
                        }
                     }
                  }
                  else
                  {
                     no_nodes = (long)m_msh_ply->nod_vector.size();
                     m_nod_x = nod_vector[m_ele->nodes_index[j]]->X();
                     m_nod_y = nod_vector[m_ele->nodes_index[j]]->Y();
                     m_nod_z = nod_vector[m_ele->nodes_index[j]]->Z();
                     m_nod = new CNode(no_nodes,m_nod_x,m_nod_y,m_nod_z);
                     m_nod->SetIndex(m_ele->nodes_index[j]);
                     m_msh_ply->nod_vector.push_back(m_nod);

                                                  // renumber
                     m_ele->nodes_index[j] = no_nodes;
                     m_ele->nodes[j] = m_nod;
                  }
               }
            }
            for (i=0; i<(long)m_msh_ply->nod_vector.size(); i++)
            {
               m_msh_ply->nod_vector[i]->SetIndex(i);
            }
         }
      }

      break;
      //====================================================================
      case 4:                                     // based on triangles / existing mesh (NW)
      {
         m_polyline->getLineVector().clear();
         ele_vector_at_ply.clear();

         //------------------------------------------------------------------
         // Search elements on the polyline
         this->GetELEOnPLY(m_polyline, ele_vector_at_ply);

         //// DEBUG CODE: highlight the collected elements
         //for(i=0;i<(long)ele_vector_at_ply.size();i++) {
         //  this->ele_vector[ele_vector_at_ply[i]]->selected=1;
         //}

         //------------------------------------------------------------------
         // Check edge if it's really necessary or not
         this->CheckMarkedEdgesOnPolyLine(m_polyline, ele_vector_at_ply);

         //------------------------------------------------------------------
         // Create line elements
         long old_element_size = (long) m_msh_ply->ele_vector.size();
         this->CreateLineElementsFromMarkedEdges(m_msh_ply, ele_vector_at_ply);

         //------------------------------------------------------------------
         // Create nodes into existing mesh
         if(m_msh_ply)
         {
            for(i=old_element_size; i<(long)m_msh_ply->ele_vector.size(); i++)
            {
               m_ele = m_msh_ply->ele_vector[i];
               m_ele->nodes.resize(m_ele->GetNodesNumber(false));

               for(j=0; j<m_ele->GetNodesNumber(false); j++)
               {
                  m_nod_x = nod_vector[m_ele->nodes_index[j]]->X();
                  m_nod_y = nod_vector[m_ele->nodes_index[j]]->Y();
                  m_nod_z = nod_vector[m_ele->nodes_index[j]]->Z();

                  long exist_node_no = 0;
                  if (m_msh_ply->HasSameCoordinatesNode(nod_vector[m_ele->nodes_index[j]], exist_node_no))
                  {
                     CNode* exist_node = m_msh_ply->nod_vector[exist_node_no];
                     m_ele->nodes[j] = exist_node;
                     m_ele->nodes_index[j] = exist_node->GetIndex();
                  }
                  else
                  {
                     no_nodes = (long)m_msh_ply->nod_vector.size();
                     m_nod = new CNode(no_nodes,m_nod_x,m_nod_y,m_nod_z);
                     m_nod->SetIndex(no_nodes);
                     m_msh_ply->nod_vector.push_back(m_nod);

                                                  // renumber
                     m_ele->nodes_index[j] = no_nodes;
                     m_ele->nodes[j] = m_nod;
                  }
               }
            }
         }
      }

      break;
      //====================================================================
   }
}


/**************************************************************************
MSHLib-Method:
Task:
Programing:
04/2005 OK Implementation
**************************************************************************/
void CFEMesh::CreateHexELEFromQuad(int no_layer,double layer_thickness)
{
   int j;
   int hempel, hampel;
   long i;
   int k;
   long size;
   long no_quad_elements = (long)ele_vector.size();
   long no_quad_nodes = (long)nod_vector.size();
   Mesh_Group::CElem* m_quad_ele = NULL;
   Mesh_Group::CElem* m_ele = NULL;

   //----------------------------------------------------------------------
   // Create MSH
   MSHDelete("HEX_from_QUAD");
   CFEMesh* m_msh_hex (new CFEMesh(_geo_obj, _geo_name));
   m_msh_hex->pcs_name = "HEX_from_QUAD";
   m_msh_hex->setElementType (MshElemType::HEXAHEDRON);
   m_msh_hex->setNumberOfMeshLayers (no_layer);
   //----------------------------------------------------------------------
   // Create HEX elements
   size = (no_layer + 1) * no_quad_nodes;
   m_msh_hex->nod_vector.resize(size);

   for(j=0;j<size;j++)
      m_msh_hex->nod_vector[j] = NULL;

   for(j=0;j<no_layer;j++)
   {
      for(i=0;i<no_quad_elements;i++)
      {
         //Elements
         m_quad_ele = ele_vector[i];
         m_ele = new Mesh_Group::CElem;
         //m_ele->SetPatchIndex(j);
         m_ele->SetIndex((j*no_quad_elements) + i);
         m_ele->SetElementType(MshElemType::HEXAHEDRON);
         m_ele->nnodes = 8;
         m_ele->nodes_index.resize(m_ele->nnodes);
         //Set indices
         m_ele->nodes_index[0] = m_quad_ele->GetNodeIndex(0) + j*no_quad_nodes;
         m_ele->nodes_index[1] = m_quad_ele->GetNodeIndex(1) + j*no_quad_nodes;
         m_ele->nodes_index[2] = m_quad_ele->GetNodeIndex(2) + j*no_quad_nodes;
         m_ele->nodes_index[3] = m_quad_ele->GetNodeIndex(3) + j*no_quad_nodes;
         m_ele->nodes_index[4] = m_ele->GetNodeIndex(0)+ no_quad_nodes;
         m_ele->nodes_index[5] = m_ele->GetNodeIndex(1)+ no_quad_nodes;
         m_ele->nodes_index[6] = m_ele->GetNodeIndex(2)+ no_quad_nodes;
         m_ele->nodes_index[7] = m_ele->GetNodeIndex(3)+ no_quad_nodes;
         //Nodes
         hempel = 0;
         hampel = 0;
         for(k=0;k<m_ele->nnodes;k++)
         {
            if(m_msh_hex->nod_vector[m_ele->GetNodeIndex(k)]==NULL)
            {
               m_ele->nodes[k] = new CNode(m_ele->GetNodeIndex(k));
               m_msh_hex->nod_vector[m_ele->GetNodeIndex(k)] = m_ele->nodes[k];
               if(k>3)
               {
                  hempel = 4;
                  hampel = 1;
               }
               m_ele->nodes[k]->SetX(nod_vector[m_quad_ele->GetNodeIndex(k-hempel)]->X());
               m_ele->nodes[k]->SetY(nod_vector[m_quad_ele->GetNodeIndex(k-hempel)]->Y());
               m_ele->nodes[k]->SetZ(nod_vector[m_quad_ele->GetNodeIndex(k-hempel)]->Z() - (j+hampel)*layer_thickness);
            }
            else
            {
               m_ele->nodes[k] = m_msh_hex->nod_vector[m_ele->GetNodeIndex(k)];
            }
         }                                        //end for m_ele->nnodes
         //msh_no_hexs++;
         m_msh_hex->ele_vector.push_back(m_ele);
      }                                           //end for no_quad_elements
   }                                              //end for no_layers
   //----------------------------------------------------------------------
   if(m_msh_hex->ele_vector.size()>0)
      fem_msh_vector.push_back(m_msh_hex);
   else
      delete m_msh_hex;

}


/**************************************************************************
MSHLib-Method:
Task:
Programing:
04/2005 OK Implementation based on MSH2RFI by WW
08/2005 WW Re-implememtation
10/2005 TK proper ordering and closing of gaps
**************************************************************************/
void GMSH2MSH(const char* filename,CFEMesh* m_msh)
{
   long id;
   long i = 0;
   int NumNodes = 0;
   int NumElements = 0;
   double x, y, z;
   std::string strbuffer;

   //WW  bool quad=false;
   //WW  CRFProcess* m_pcs = NULL;
   CNode* node = NULL;
   CElem* elem = NULL;
   std::ifstream msh_file(filename, std::ios::in);
   getline(msh_file, strbuffer);                  // Node keyword

   // OLD GMSH  FORMAT----------------------------------------------------------------------
   if (strbuffer.compare("$NOD") == 0)
   {
      while (strbuffer.compare("$ENDELM") != 0)
      {
         msh_file >> NumNodes >> std::ws;
         //....................................................................
         // Node data
         for (i = 0; i < NumNodes; i++)
         {
            msh_file >> id >> x >> y >> z >> std::ws;

            node = new CNode(id, x, y, z);
            m_msh->nod_vector.push_back(node);
         }

         getline(msh_file, strbuffer);            // End Node keyword
         //....................................................................
         // Element data
         getline(msh_file, strbuffer);            // Element keyword
         msh_file >> NumElements >> std::ws;
         for (i = 0; i < NumElements; i++)
         {
            elem = new CElem(i);
            elem->Read(msh_file, 2);
            m_msh->ele_vector.push_back(elem);
         }
         getline(msh_file, strbuffer);            // END keyword

         // ordering nodes and closing gaps TK
         std::vector<int> gmsh_id;
         long new_node_id;
         int counter = 0;
         int diff = 0;
         int j = 0;
         for (i = 0; i < (int) m_msh->nod_vector.size(); i++)
         {
            diff = m_msh->nod_vector[i]->GetIndex() - counter;
            if (diff == 0)
            {
               gmsh_id.push_back(i);
               counter++;
            }
            else
            {
               for (j = 0; j < diff; j++)
               {
                  gmsh_id.push_back(i);
                  counter++;
               }
               i--;
            }
         }

         for (i = 0; i < (int) m_msh->ele_vector.size(); i++)
         {
            for (j = 0; j < (int) m_msh->ele_vector[i]->GetVertexNumber(); j++)
            {
               new_node_id = gmsh_id[m_msh->ele_vector[i]->GetNodeIndex(j)
                  + 1];
               //m_msh->ele_vector[i]->nodes[j]->SetIndex(new_node_id);/*global*/
                                                  /*local*/
               m_msh->ele_vector[i]->nodes_index[j] = new_node_id;
            }
         }
         for (i = 0; i < (int) m_msh->nod_vector.size(); i++)
         {
            m_msh->nod_vector[i]->SetIndex(i);
         }
         // END OF: ordering nodes and closing gaps TK

      }                                           /*End while*/
   }
   // END old GMSH Format----------------------------------------------------------------------

   // NEW 2008 GMSH  FORMAT----------------------------------------------------------------------
   if (strbuffer.compare("$MeshFormat") == 0)
   {
      getline(msh_file, strbuffer);               // version-number file-type data-size
      getline(msh_file, strbuffer);               //$EndMeshFormat
      getline(msh_file, strbuffer);               //$Nodes Keywords

      while (strbuffer.compare("$EndElements") != 0)
      {
         // Node data
         msh_file >> NumNodes >> std::ws;
         for (i = 0; i < NumNodes; i++)
         {
            msh_file >> id >> x >> y >> z >> std::ws;
            node = new CNode(id, x, y, z);
            m_msh->nod_vector.push_back(node);
         }
         getline(msh_file, strbuffer);            // End Node keyword $EndNodes

         // Element data
         getline(msh_file, strbuffer);            // Element keyword $Elements
         msh_file >> NumElements >> std::ws;      // number-of-elements
         for (i = 0; i < NumElements; i++)
         {
            elem = new CElem(i);
            elem->Read(msh_file, 7);
            if (elem->GetElementType() != MshElemType::INVALID)
               m_msh->ele_vector.push_back(elem);
         }
         getline(msh_file, strbuffer);            // END keyword

         // correct indices TF
         const size_t n_elements(m_msh->ele_vector.size());
         for (size_t k(0); k < n_elements; k++)
         {
            m_msh->ele_vector[k]->SetIndex(k);
         }

         // ordering nodes and closing gaps TK
         std::vector<int> gmsh_id;
         long new_node_id;
         int counter = 0;
         int diff = 0;
         int j = 0;
         for (i = 0; i < (int) m_msh->nod_vector.size(); i++)
         {
            diff = m_msh->nod_vector[i]->GetIndex() - counter;
            if (diff == 0)
            {
               gmsh_id.push_back(i);
               counter++;
            }
            else
            {
               for (j = 0; j < diff; j++)
               {
                  gmsh_id.push_back(i);
                  counter++;
               }
               i--;
            }
         }

         for (i = 0; i < (int) m_msh->ele_vector.size(); i++)
         {
            for (j = 0; j < (int) m_msh->ele_vector[i]->GetVertexNumber(); j++)
            {
               new_node_id = gmsh_id[m_msh->ele_vector[i]->GetNodeIndex(j)
                  + 1];
               //m_msh->ele_vector[i]->nodes[j]->SetIndex(new_node_id);/*global*/
                                                  /*local*/
               m_msh->ele_vector[i]->nodes_index[j] = new_node_id;
            }
         }
         for (i = 0; i < (int) m_msh->nod_vector.size(); i++)
         {
            m_msh->nod_vector[i]->SetIndex(i);
         }
         // END OF: ordering nodes and closing gaps TK

      }                                           /*End while*/
   }
   // END New 2008 GMSH Format----------------------------------------------------------------------

   //  m_msh->ConstructGrid(); // TF

   msh_file.close();
}


/**************************************************************************
MSHLib-Method:
Task:   Takes Surface, meshs it and saves it with defined file name

        surface_name         = name of surface for meshing
        file_name_const_char = path + name without extension

        Attention: You have to remove/delete files extern.

Programing:
12/2005 TK implementation
**************************************************************************/
void Mesh_Single_Surface(std::string surface_name, const char *file_name_const_char)
{

   //Searching Methods
   int i=0, j=0, k=0;
   std::string Name;
   std::vector<CGLPoint*> surface_points_searchvector;
   CGLPoint *m_point = NULL;

   for (i=0; i<(int)surface_points_searchvector.size();i++)
   {
      delete surface_points_searchvector[i];
   }
   surface_points_searchvector.clear();

   //Get SFC
   for (i=0;i<(int)surface_vector.size();i++)
   {
      if (surface_vector[i]->name.find(surface_name)==0)
      {
         for (j=0;j<(int)surface_vector[i]->polyline_of_surface_vector.size();j++)
         {
            if (j==0)
            {
               for (k=0;k<(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size();k++)
               {
                  m_point = new CGLPoint;
                  m_point->nb_of_ply = i;
                  m_point->id = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->id;
                  m_point->x  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->x;
                  m_point->y  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->y;
                  m_point->z  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->z;
                  surface_points_searchvector.push_back(m_point);
                  if (k==(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1 &&
                     surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->id == surface_vector[i]->polyline_of_surface_vector[j]->point_vector[0]->id)
                  {
                     surface_points_searchvector.erase(surface_points_searchvector.begin()+k);
                  }
               }
            }

            else
            {
               if (m_point->id != surface_vector[i]->polyline_of_surface_vector[j]->point_vector[0]->id &&
                  m_point->id == surface_vector[i]->polyline_of_surface_vector[j]->point_vector[(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1]->id)
               {
                  for (k=1;k<(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size();k++)
                  {
                     m_point = new CGLPoint;
                     m_point->nb_of_ply = i;
                     m_point->id = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1-k]->id;
                     m_point->x  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1-k]->x;
                     m_point->y  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1-k]->y;
                     m_point->z  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size()-1-k]->z;
                     surface_points_searchvector.push_back(m_point);
                  }
               }
               else
               {
                  for (k=1;k<(int)surface_vector[i]->polyline_of_surface_vector[j]->point_vector.size();k++)
                  {
                     m_point = new CGLPoint;
                     m_point->nb_of_ply = i;
                     m_point->id = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->id;
                     m_point->x  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->x;
                     m_point->y  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->y;
                     m_point->z  = surface_vector[i]->polyline_of_surface_vector[j]->point_vector[k]->z;
                     surface_points_searchvector.push_back(m_point);
                  }
               }
            }
         }

         if (surface_points_searchvector[0]->id == surface_points_searchvector[surface_points_searchvector.size()-1]->id)
            surface_points_searchvector.erase(surface_points_searchvector.begin());

         //Write GMSH_GEO_FILE of marked Surface and mesh it
         std::string m_strFileNameGEO = file_name_const_char;
         m_strFileNameGEO = m_strFileNameGEO + ".geo";
         file_name_const_char = m_strFileNameGEO.data();
         FILE *geo_file=NULL;
         geo_file = fopen(file_name_const_char, "w+t");
         long id=0;
         double density = 1e100;
         double topologic_distance;

         for (k=0;k<(int)surface_points_searchvector.size();k++)
         {
            if (k==0)
            {
               density =   EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[k+1]->x,
                  surface_points_searchvector[k+1]->y,
                  surface_points_searchvector[k+1]->z);
               topologic_distance = EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[surface_points_searchvector.size()-1]->x,
                  surface_points_searchvector[surface_points_searchvector.size()-1]->y,
                  surface_points_searchvector[surface_points_searchvector.size()-1]->z);
               if (topologic_distance < density) density = topologic_distance;
            }

            if (k>0 && k<(int)surface_points_searchvector.size()-1)
            {
               density =   EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[k-1]->x,
                  surface_points_searchvector[k-1]->y,
                  surface_points_searchvector[k-1]->z);
               topologic_distance = EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[k+1]->x,
                  surface_points_searchvector[k+1]->y,
                  surface_points_searchvector[k+1]->z);
               if (topologic_distance < density) density = topologic_distance;
            }

            if (k == (int)surface_points_searchvector.size()-1)
            {
               density =   EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[k-1]->x,
                  surface_points_searchvector[k-1]->y,
                  surface_points_searchvector[k-1]->z);
               topologic_distance = EuklVek3dDistCoor ( surface_points_searchvector[k]->x,
                  surface_points_searchvector[k]->y,
                  surface_points_searchvector[k]->z,
                  surface_points_searchvector[0]->x,
                  surface_points_searchvector[0]->y,
                  surface_points_searchvector[0]->z);
               if (topologic_distance < density) density = topologic_distance;
            }

            id++;
            fprintf(geo_file,"%s","Point(");
            fprintf(geo_file,"%li",id);
            fprintf(geo_file,"%s",") = {");
            fprintf(geo_file,"%20.14f",surface_points_searchvector[k]->x);
            fprintf(geo_file,"%s",", ");
            fprintf(geo_file,"%20.14f",surface_points_searchvector[k]->y);
            fprintf(geo_file,"%s",", ");
            fprintf(geo_file,"%20.14f",surface_points_searchvector[k]->z);
            fprintf(geo_file,"%s",", ");
            fprintf(geo_file,"%g",density/4);
            fprintf(geo_file,"%s\n","};");

         }
         id=0;
         for (k=0;k<(int)surface_points_searchvector.size()-1;k++)
         {
            id++;
            fprintf(geo_file,"%s","Line(");
            fprintf(geo_file,"%li",id);
            fprintf(geo_file,"%s",") = {");
            fprintf(geo_file,"%d",k+1);
            fprintf(geo_file,"%s",", ");
            fprintf(geo_file,"%d",k+2);
            fprintf(geo_file,"%s\n","};");

         }
         id++;
         fprintf(geo_file,"%s","Line(");
         fprintf(geo_file,"%li",id);
         fprintf(geo_file,"%s",") = {");
         fprintf(geo_file,"%ld",id);
         fprintf(geo_file,"%s",", 1");
         fprintf(geo_file,"%s\n","};");

         fprintf(geo_file,"%s","Line Loop(");
         fprintf(geo_file,"%li",id+1);
         fprintf(geo_file,"%s",") = {");
         for (k=0;k<(int)surface_points_searchvector.size();k++)
         {
            fprintf(geo_file,"%i",k+1);
            if (k<(int)surface_points_searchvector.size()-1)
               fprintf(geo_file,"%s",", ");
         }
         fprintf(geo_file,"%s\n","};");

         fprintf(geo_file,"%s","Plane Surface(");
         fprintf(geo_file,"%li",id+2);
         fprintf(geo_file,"%s",") = {");
         fprintf(geo_file,"%ld",id+1);
         fprintf(geo_file,"%s\n","};");

         fprintf(geo_file,"%s","Physical Surface(");
         fprintf(geo_file,"%li",id+3);
         fprintf(geo_file,"%s",") = {");
         fprintf(geo_file,"%ld",id+2);
         fprintf(geo_file,"%s\n","};");

         fclose(geo_file);

         std::string m_strExecuteGEO = "gmsh " + m_strFileNameGEO +" -2";

         // PCH & TK: Workaround for the old problem.

         //system(m_strExecute);
         //remove(file_name_const_char);
         //END Meshing
         break;
      }
   }
}


///**************************************************************************
//MSHLib-Method:
//Task:    file_name_const_char = Path + Name  without extension of GMSH *msh-File
//
//Programing:
//12/2005 TK implementation
//**************************************************************************/
//void Select_Nodes_Elements_by_TINFile(const char *file_name_const_char)
//{
//   int i=0, j=0, k=0;
////READ GMSH-File and fill local Element Vector
//  vector<Mesh_Group::CFEMesh*>check_msh_vector;
////  Mesh_Group::CFEMesh* m_check_elements;
////WW  char text[1024];
//  long id_elem;
//
//  string m_strFileNameTIN = file_name_const_char;
//  m_strFileNameTIN = m_strFileNameTIN + ".tin";
//  file_name_const_char = m_strFileNameTIN.data();
//  ifstream tin_file (file_name_const_char,ios::in);
//  ifstream tin2check (file_name_const_char,ios::in);
////  m_check_elements = new Mesh_Group::CFEMesh;
////Loop over all generated triangles of surface
//  CGLPoint point;
//  double angle_sum, dist;
//  double tolerance = 0.001;
//  double tri_point1[3],tri_point2[3],tri_point3[3],checkpoint[3];
//  double tri_x[3],tri_y[3],tri_z[3];
//  double min_mesh_dist=0.0;
//  double sfc_min[3],sfc_max[3];
//
//  while (!tin2check.eof())
//  {
//    i=tin2check.tellg();
//    tin2check>>id_elem>>tri_point1[0]>>tri_point1[1]>>tri_point1[2]>>tri_point2[0]>>tri_point2[1]>>tri_point2[2]>>tri_point3[0]>>tri_point3[1]>>tri_point3[2];
//
//    if (i==0)
//    {
//     sfc_min[0]= tri_point1[0];
//     sfc_min[1]= tri_point1[1];
//     sfc_min[2]= tri_point1[2];
//     sfc_max[0]= tri_point1[0];
//     sfc_max[1]= tri_point1[1];
//     sfc_max[2]= tri_point1[2];
//     if (tri_point1[0] < sfc_min[0]) sfc_min[0] = tri_point1[0];
//     if (tri_point2[0] < sfc_min[0]) sfc_min[0] = tri_point2[0];
//     if (tri_point3[0] < sfc_min[0]) sfc_min[0] = tri_point3[0];
//     if (tri_point1[0] > sfc_max[0]) sfc_max[0] = tri_point1[0];
//     if (tri_point2[0] > sfc_max[0]) sfc_max[0] = tri_point2[0];
//     if (tri_point3[0] > sfc_max[0]) sfc_max[0] = tri_point3[0];
//     if (tri_point1[1] < sfc_min[1]) sfc_min[1] = tri_point1[1];
//     if (tri_point2[1] < sfc_min[1]) sfc_min[1] = tri_point2[1];
//     if (tri_point3[1] < sfc_min[1]) sfc_min[1] = tri_point3[1];
//     if (tri_point1[1] > sfc_max[1]) sfc_max[1] = tri_point1[1];
//     if (tri_point2[1] > sfc_max[1]) sfc_max[1] = tri_point2[1];
//     if (tri_point3[1] > sfc_max[1]) sfc_max[1] = tri_point3[1];
//     if (tri_point1[2] < sfc_min[2]) sfc_min[2] = tri_point1[2];
//     if (tri_point2[2] < sfc_min[2]) sfc_min[2] = tri_point2[2];
//     if (tri_point3[2] < sfc_min[2]) sfc_min[2] = tri_point3[2];
//     if (tri_point1[2] > sfc_max[2]) sfc_max[2] = tri_point1[2];
//     if (tri_point2[2] > sfc_max[2]) sfc_max[2] = tri_point2[2];
//     if (tri_point3[2] > sfc_max[2]) sfc_max[2] = tri_point3[2];
//    }
//    else
//    {
//     if (tri_point1[0] < sfc_min[0]) sfc_min[0] = tri_point1[0];
//     if (tri_point2[0] < sfc_min[0]) sfc_min[0] = tri_point2[0];
//     if (tri_point3[0] < sfc_min[0]) sfc_min[0] = tri_point3[0];
//     if (tri_point1[0] > sfc_max[0]) sfc_max[0] = tri_point1[0];
//     if (tri_point2[0] > sfc_max[0]) sfc_max[0] = tri_point2[0];
//     if (tri_point3[0] > sfc_max[0]) sfc_max[0] = tri_point3[0];
//     if (tri_point1[1] < sfc_min[1]) sfc_min[1] = tri_point1[1];
//     if (tri_point2[1] < sfc_min[1]) sfc_min[1] = tri_point2[1];
//     if (tri_point3[1] < sfc_min[1]) sfc_min[1] = tri_point3[1];
//     if (tri_point1[1] > sfc_max[1]) sfc_max[1] = tri_point1[1];
//     if (tri_point2[1] > sfc_max[1]) sfc_max[1] = tri_point2[1];
//     if (tri_point3[1] > sfc_max[1]) sfc_max[1] = tri_point3[1];
//     if (tri_point1[2] < sfc_min[2]) sfc_min[2] = tri_point1[2];
//     if (tri_point2[2] < sfc_min[2]) sfc_min[2] = tri_point2[2];
//     if (tri_point3[2] < sfc_min[2]) sfc_min[2] = tri_point3[2];
//     if (tri_point1[2] > sfc_max[2]) sfc_max[2] = tri_point1[2];
//     if (tri_point2[2] > sfc_max[2]) sfc_max[2] = tri_point2[2];
//     if (tri_point3[2] > sfc_max[2]) sfc_max[2] = tri_point3[2];
//    }
//  }
//  tin2check.close();
//
//  CFEMesh* m_msh = NULL;
//  m_msh = new CFEMesh();
//  CNode* node = NULL;
//  fem_msh_vector.push_back(m_msh);
//  int temp_mesh = (long)fem_msh_vector.size();
//  //Loop over all meshes
//    for(j=0;j<(long)fem_msh_vector.size()-1;j++)
//    {
//    //Loop over all edges
//        for(i=0;i<(long)fem_msh_vector[j]->edge_vector.size();i++)
//        {
//            if (j==0 && i==0){
//              min_mesh_dist = fem_msh_vector[j]->edge_vector[i]->Length();
//            }
//            else{
//              if (min_mesh_dist  > fem_msh_vector[j]->edge_vector[i]->Length())
//                  min_mesh_dist = fem_msh_vector[j]->edge_vector[i]->Length();
//            }
//        }
//        tolerance = min_mesh_dist;
//    //Loop over all mesh nodes
//        for(i=0;i<(long)fem_msh_vector[j]->nod_vector.size();i++)
//        {
//            checkpoint[0] = fem_msh_vector[j]->nod_vector[i]->X();
//            checkpoint[1] = fem_msh_vector[j]->nod_vector[i]->Y();
//            checkpoint[2] = fem_msh_vector[j]->nod_vector[i]->Z();
//            node = new CNode(i,checkpoint[0],checkpoint[1],checkpoint[2]);
//            if((checkpoint[0]>=sfc_min[0]-tolerance && checkpoint[0]<=sfc_max[0]+tolerance )&&
//               (checkpoint[1]>=sfc_min[1]-tolerance && checkpoint[1]<=sfc_max[1]+tolerance )&&
//               (checkpoint[2]>=sfc_min[2]-tolerance && checkpoint[2]<=sfc_max[2]+tolerance ) )
//            {
//                fem_msh_vector[temp_mesh-1]->nod_vector.push_back(node);
//            }
//        }
//    }
//
//
//  while (!tin_file.eof())
//  {
//    //tin_file.getline(text, 1024);
//    tin_file>>id_elem>>tri_point1[0]>>tri_point1[1]>>tri_point1[2]>>tri_point2[0]>>tri_point2[1]>>tri_point2[2]>>tri_point3[0]>>tri_point3[1]>>tri_point3[2];
//         tri_x[0]=tri_point1[0];
//         tri_x[1]=tri_point2[0];
//         tri_x[2]=tri_point3[0];
//         tri_y[0]=tri_point1[1];
//         tri_y[1]=tri_point2[1];
//         tri_y[2]=tri_point3[1];
//         tri_z[0]=tri_point1[2];
//         tri_z[1]=tri_point2[2];
//         tri_z[2]=tri_point3[2];
//    //Loop over all preselected mesh nodes
//        for(i=0;i<(long)fem_msh_vector[temp_mesh-1]->nod_vector.size();i++)
//        {
//            point.x = fem_msh_vector[temp_mesh-1]->nod_vector[i]->X();
//            point.y = fem_msh_vector[temp_mesh-1]->nod_vector[i]->Y();
//            point.z = fem_msh_vector[temp_mesh-1]->nod_vector[i]->Z();
//            checkpoint[0] = fem_msh_vector[temp_mesh-1]->nod_vector[i]->X();
//            checkpoint[1] = fem_msh_vector[temp_mesh-1]->nod_vector[i]->Y();
//            checkpoint[2] = fem_msh_vector[temp_mesh-1]->nod_vector[i]->Z();
//            dist = MCalcDistancePointToPlane(checkpoint,tri_point1,tri_point2,tri_point3);
//            if (k==0) fem_msh_vector[temp_mesh-1]->nod_vector[i]->epsilon = dist;
//            else
//            {
//                if (fem_msh_vector[temp_mesh-1]->nod_vector[i]->epsilon > dist)
//                    fem_msh_vector[temp_mesh-1]->nod_vector[i]->epsilon = dist;
//            }
//                if (dist<=tolerance && dist>=-tolerance)
//                {
//                  angle_sum = AngleSumPointInsideTriangle(checkpoint,tri_point1,tri_point2,tri_point3, min_mesh_dist);
//                  if(angle_sum>359)
//                  fem_msh_vector[temp_mesh-1]->nod_vector[i]->selected = 1;
//                }
//        }
//        // Loop over all mesh elements
//        for(i=0;i<(long)fem_msh_vector[temp_mesh-2]->ele_vector.size();i++)
//        {
//            if (fem_msh_vector[temp_mesh-2]->ele_vector[i]->GetElementType() == 2 ||
//                fem_msh_vector[temp_mesh-2]->ele_vector[i]->GetElementType() == 4) /*Quad or Tri*/
//            {
//            double* xyz;
//            xyz = fem_msh_vector[temp_mesh-2]->ele_vector[i]->GetGravityCenter();
//            checkpoint[0] = xyz[0];
//            checkpoint[1] = xyz[1];
//            checkpoint[2] = xyz[2];
//            angle_sum = AngleSumPointInsideTriangle(checkpoint,tri_point1,tri_point2,tri_point3, min_mesh_dist/10);
//
//            if(angle_sum>359.8)
//            fem_msh_vector[temp_mesh-2]->ele_vector[i]->selected = 1;
//            }
//        }
//  }
//
////WW  int a = (int)fem_msh_vector[temp_mesh-1]->nod_vector.size();
//  int index;
//  //Loop over all meshes
//    for(j=0;j<(long)fem_msh_vector.size()-1;j++)
//    {
//    //Loop over selected nodes
//        for(i=0;i<(long)fem_msh_vector[temp_mesh-1]->nod_vector.size();i++)
//        {
//            index = fem_msh_vector[temp_mesh-1]->nod_vector[i]->GetIndex();
//            if(index < (int)fem_msh_vector[j]->nod_vector.size())
//            {
//            if ((fem_msh_vector[temp_mesh-1]->nod_vector[i]->GetIndex() == fem_msh_vector[j]->nod_vector[index]->GetIndex())
//                &&
//                fem_msh_vector[temp_mesh-1]->nod_vector[i]->selected==1
//                &&
//                (fem_msh_vector[temp_mesh-1]->nod_vector[i]->X() == fem_msh_vector[j]->nod_vector[index]->X())
//                &&
//                (fem_msh_vector[temp_mesh-1]->nod_vector[i]->Y() == fem_msh_vector[j]->nod_vector[index]->Y())
//                &&
//                (fem_msh_vector[temp_mesh-1]->nod_vector[i]->Z() == fem_msh_vector[j]->nod_vector[index]->Z()))
//            {
//                fem_msh_vector[j]->nod_vector[index]->selected = 1;
//            }
//            }
//        }
//
//        // Loop over all mesh elements
//        Math_Group::vec<long> node_index(20);
//        for(i=0;i<(long)fem_msh_vector[j]->ele_vector.size();i++)
//        {
//            fem_msh_vector[j]->ele_vector[i]->GetNodeIndeces(node_index);
//
//            if (fem_msh_vector[j]->ele_vector[i]->GetElementType() == 0) /*Quad*/
//            {
//                if (fem_msh_vector[j]->nod_vector[node_index[0]]->selected == 1 ||
//                    fem_msh_vector[j]->nod_vector[node_index[1]]->selected == 1 ||
//                    fem_msh_vector[j]->nod_vector[node_index[2]]->selected == 1 ||
//                    fem_msh_vector[j]->nod_vector[node_index[3]]->selected == 1 )
//                    fem_msh_vector[j]->ele_vector[i]->selected = 1;
//            }
//
//            if (fem_msh_vector[j]->ele_vector[i]->GetElementType() == 0) /*TRI*/
//            {
//                if (fem_msh_vector[j]->nod_vector[node_index[0]]->selected == 1 &&
//                    fem_msh_vector[j]->nod_vector[node_index[1]]->selected == 1 &&
//                    fem_msh_vector[j]->nod_vector[node_index[2]]->selected == 1 )
//                    fem_msh_vector[j]->ele_vector[i]->selected = 1;
//            }
//        }
//    }
//
//
//    fem_msh_vector.erase(fem_msh_vector.begin()+ temp_mesh-1);
//
//}

/**************************************************************************
MSHLib-Method:
Task:    Clears the selection and set flag selected = 0;
Programing:
12/2005 TK implementation
**************************************************************************/
void Clear_Selected_Nodes_Elements()
{
   int i=0, j=0;
   for(j=0;j<(long)fem_msh_vector.size();j++)
   {
      for(i=0;i<(long)fem_msh_vector[j]->nod_vector.size();i++)
      {
         fem_msh_vector[j]->nod_vector[i]->epsilon = 0.0;
         fem_msh_vector[j]->nod_vector[i]->selected= 0;
      }

      for(i=0;i<(long)fem_msh_vector[j]->ele_vector.size();i++)
      {
         fem_msh_vector[j]->ele_vector[i]->selected= 0;
      }
   }
}


/**************************************************************************
GeoSys-Method:
Task:
Programing:
11/2005 MB
**************************************************************************/
void CFEMesh::SetMSHPart(std::vector<long>&elements_active, long StrangNumber)
{
   int j;
   int k;
   long i;
   long size;
   CElem* m_ele = NULL;
   CFEMesh* m_msh_strang = NULL;
   CNode* m_nod = NULL;
   bool found = false;
   StrangNumber= StrangNumber;
   size = (long)elements_active.size();

   m_msh_strang = FEMGet("MSH_Strang");

   // Create MSH
   if (!m_msh_strang)
   {
      m_msh_strang = new CFEMesh(_geo_obj, _geo_name);
      m_msh_strang->pcs_name = "MSH_Strang";
      m_msh_strang->setElementType (MshElemType::LINE);
      m_msh_strang->setNumberOfMeshLayers (getNumberOfMeshLayers());
      //Resize
      m_msh_strang->ele_vector.resize(size);
      m_msh_strang->Eqs2Global_NodeIndex.resize(size+1);

      fem_msh_vector.push_back(m_msh_strang);
   }

   m_msh_strang->nod_vector.resize(0);

   //Elements
   for(i=0;i<(long)elements_active.size();i++)
   {
      //Elements
      m_ele = ele_vector[elements_active[i]];
      m_msh_strang->ele_vector[i] = m_ele;
      //Nodes
      for(k=0;k<m_ele->nnodes;k++)
      {
         m_nod = m_ele->GetNode(k);
         found = false;
         for(j=0; j<(int)m_msh_strang->nod_vector.size(); j++)
         {
            if(*m_msh_strang->nod_vector[j]==*m_nod)
               found = true;
         }
         if(!found)
            m_msh_strang->nod_vector.push_back(m_nod);
      }
      // cout<< "Element"  << i << "  " << m_ele->GetNodeIndex(0)<< "  "  << m_ele->GetNodeIndex(1) <<  endl;
   }

   for(i=0; i<(long)m_msh_strang->nod_vector.size(); i++)
   {
      m_msh_strang->nod_vector[i]->SetEquationIndex(i);
                                                  //+ test;
      m_msh_strang->Eqs2Global_NodeIndex[i] = m_msh_strang->nod_vector[i]->GetIndex();
      //cout<< " " << i << "  " << m_msh_strang->Eqs2Global_NodeIndex[i] <<  endl;
   }

}


/**************************************************************************
MSHLib-Method:
Task:   GMSH 2 TIN
        const char *file_name_const_char = Location and Name of GMSH FILE
Programing:
02/2006 TK implementation
**************************************************************************/
//void GMSH2TIN(const char *file_name_const_char)
//{
//   int i=0, k=0;
////READ GMSH-File and fill local Element Vector
//  vector<Mesh_Group::CFEMesh*>check_msh_vector;
//  Mesh_Group::CFEMesh* m_check_elements;
//  char text[1024];
//  long nbnod, nbelm;
//  long node_id;
//  long pnt;
//  double x,y,z;
//
//
//  string m_strFileNameGMSH = file_name_const_char;
//  string m_strFileNameTIN = file_name_const_char;
//  m_strFileNameGMSH = m_strFileNameGMSH + ".msh";
//  m_strFileNameTIN = m_strFileNameTIN + ".tin";
//  file_name_const_char = m_strFileNameGMSH.data();
//  ifstream gmsh_file (file_name_const_char,ios::in);
//
//  m_check_elements = new Mesh_Group::CFEMesh;
//  while (!gmsh_file.eof())
//  {
//    gmsh_file.getline(text, 1024);
//    if (!strncmp(text,"$NOD",4)){
//        gmsh_file>>nbnod>>ws;
//        m_check_elements->nod_vector.resize(nbnod);
//        for(i=0;i<nbnod;i++)
//        {
//            gmsh_file>>node_id>>x>>y>>z>>ws;
//            m_check_elements->nod_vector[i] = new CNode(i,x,y,z);
//        }
//
//    }
//    if (!strncmp(text,"$ELM",4)){
//        gmsh_file>>nbelm>>ws;
//        m_check_elements->ele_vector.resize(nbelm);
//        for(i=0;i<nbelm;i++)
//        {
//            m_check_elements->ele_vector[i] = new CElem(i);
//	        m_check_elements->ele_vector[i]->Read(gmsh_file, 2);
//        }
//
//    }
//   if (!strncmp(text,"$ENDELM",7)) break;
//  }
//
////Loop over all generated triangles of surface
//  double tri_point1[3],tri_point2[3],tri_point3[3];
//  FILE *tin_file=NULL;
//  file_name_const_char = m_strFileNameTIN.data();
//  tin_file = fopen(file_name_const_char, "w+t");
//  for(k=0;k<(int)m_check_elements->ele_vector.size();k++)
//  {
//        pnt = m_check_elements->ele_vector[k]->GetNodeIndex(0);
//         tri_point1[0] = m_check_elements->nod_vector[pnt]->X();
//         tri_point1[1] = m_check_elements->nod_vector[pnt]->Y();
//         tri_point1[2] = m_check_elements->nod_vector[pnt]->Z();
//        pnt = m_check_elements->ele_vector[k]->GetNodeIndex(1);
//         tri_point2[0] = m_check_elements->nod_vector[pnt]->X();
//         tri_point2[1] = m_check_elements->nod_vector[pnt]->Y();
//         tri_point2[2] = m_check_elements->nod_vector[pnt]->Z();
//        pnt = m_check_elements->ele_vector[k]->GetNodeIndex(2);
//         tri_point3[0] = m_check_elements->nod_vector[pnt]->X();
//         tri_point3[1] = m_check_elements->nod_vector[pnt]->Y();
//         tri_point3[2] = m_check_elements->nod_vector[pnt]->Z();
//
//		    fprintf(tin_file,"%i ",k);
//            fprintf(tin_file,"%lf %lf %lf ",tri_point1[0], tri_point1[1], tri_point1[2]);
//            fprintf(tin_file,"%lf %lf %lf ",tri_point2[0], tri_point2[1], tri_point2[2]);
//            fprintf(tin_file,"%lf %lf %lf\n",tri_point3[0], tri_point3[1], tri_point3[2]);
//
//  }
//  fclose(tin_file);
//}

/**************************************************************************
MSHLib-Method:
Task:   CheckMarkedEdgesOnPolyLine
Programing:
05/2007 NW implementation
**************************************************************************/
void CFEMesh::CheckMarkedEdgesOnPolyLine(CGLPolyline*m_polyline, std::vector<long> &ele_vector_at_ply)
{
   CElem* m_ele = NULL;
   CEdge* m_edg = NULL;
   vec<CEdge*>ele_edges_vector(15);
   vec<CNode*>edge_nodes(3);

   //------------------------------------------------------------------
   // Make a list of elements to check
   std::vector<int> chk_ele_list;
   for (int i=0;i<(long)ele_vector_at_ply.size();i++)
   {
      m_ele = ele_vector[ele_vector_at_ply[i]];
      m_ele->GetEdges(ele_edges_vector);

      int marked_count = 0;
      for(int j=0;j<(int)m_ele->GetEdgesNumber();j++)
      {
         m_edg = ele_edges_vector[j];
         if(m_edg->GetMark()) marked_count++;
      }

      if (marked_count == m_ele->GetEdgesNumber())
      {
         //if all edges of one element are marked, those edges need to be checked.
         chk_ele_list.push_back(i);
      }
   }

   //------------------------------------------------------------------
   // Check continuity of node indexes that compose one edge
   const double search_tolerance = 1.0e-10;
   for (int i=0; i<(int)chk_ele_list.size(); i++)
   {
      m_ele = ele_vector[ele_vector_at_ply[chk_ele_list[i]]];
      m_ele->GetEdges(ele_edges_vector);
      for(int j=0; j<(int)m_ele->GetEdgesNumber(); j++)
      {
         m_edg = ele_edges_vector[j];
         if(!m_edg->GetMark()) continue;

         int index[2] = {-2,-2};
         for (int k=0;k<(int)m_polyline->point_vector.size();k++)
         {
            double distance1 = EuklVek3dDistCoor(m_edg->GetNode(0)->X(), m_edg->GetNode(0)->Y(), m_edg->GetNode(0)->Z(),
               m_polyline->point_vector[k]->x, m_polyline->point_vector[k]->y, m_polyline->point_vector[k]->z);
            double distance2 = EuklVek3dDistCoor(m_edg->GetNode(1)->X(), m_edg->GetNode(1)->Y(), m_edg->GetNode(1)->Z(),
               m_polyline->point_vector[k]->x, m_polyline->point_vector[k]->y, m_polyline->point_vector[k]->z);

            if (search_tolerance > distance1)
            {
               index[0] = k;
            }
            else if (search_tolerance > distance2)
            {
               index[1] = k;
            }
         }
         if (index[0] >= 0 && index[1] >= 0)
         {
            if (abs(index[0]-index[1]) != 1)      // not continuity
            {
               m_edg->SetMark(false);             //remove
            }
         }
      }
   }

   //------------------------------------------------------------------
   //Check the connectivity of marked edge with others
   for (int i=0; i<(int)chk_ele_list.size(); i++)
   {
      m_ele = ele_vector[ele_vector_at_ply[chk_ele_list[i]]];
      m_ele->GetEdges(ele_edges_vector);

      int marked_count = 0;
      for(int j=0; j<(int)m_ele->GetEdgesNumber(); j++)
      {
         m_edg = ele_edges_vector[j];
         if(m_edg->GetMark()) marked_count++;
      }
      if (marked_count < m_ele->GetEdgesNumber())
      {
         continue;
      }

      vec<CNode*>n_edge_nodes(3);
      vec<CEdge*>n_ele_edges_vector(15);
      int node_use[3] = {0,0,0};
      for (int j=0;j<(long)ele_vector_at_ply.size();j++)
      {
         if (i==j) continue;
         CElem* n_ele = ele_vector[ele_vector_at_ply[j]];
         n_ele->GetEdges(n_ele_edges_vector);
         CEdge* n_edg;

         for(int k=0;k<(int)n_ele->GetEdgesNumber();k++)
         {
            n_edg = n_ele_edges_vector[k];
            if(!n_edg->GetMark()) continue;

            //exclude the elements sharing the marked edge
            bool ret=true;
            for (int l=0;l<(int)m_ele->GetEdgesNumber();l++)
            {
               if (!ele_edges_vector[l]->GetMark()) continue;
               if (n_edg == ele_edges_vector[l])
               {
                  ret = false;
                  break;
               }
            }
            if (!ret) continue;

            n_edg->GetNodes(n_edge_nodes);
            for (int l=0;l<3;l++)
            {
               if (n_edge_nodes[0]->GetIndex() == (size_t)m_ele->nodes_index[l]
                  || n_edge_nodes[1]->GetIndex() == (size_t)m_ele->nodes_index[l])
               {
                  node_use[l]+=1;
               }
            }
         }

         long node1 = -1;
         long node2 = -1;
         for (int j=0;j<3;j++)
         {
            if (node_use[j] == 0) continue;
            if (node1 < 0)
            {
               node1 = m_ele->nodes_index[j];
            }
            else
            {
               node2 = m_ele->nodes_index[j];
            }
         }

         if (node1 >= 0 && node2 >= 0)
         {
            for (int j=0;j<m_ele->GetEdgesNumber();j++)
            {
               m_edg = ele_edges_vector[j];
               m_edg->GetNodes(edge_nodes);
               if ((edge_nodes[0]->GetIndex() == (size_t)node1 &&  edge_nodes[1]->GetIndex() == (size_t)node2)
                  || (edge_nodes[1]->GetIndex() == (size_t)node1 &&  edge_nodes[0]->GetIndex() == (size_t)node2))
               {
                  m_edg->SetMark(false);
               }
            }
         }
      }
   }
}


/**************************************************************************
MSHLib-Method:
Task:   CreateLineElementsFromMarkedEdges
Programing:
05/2007 NW implementation
**************************************************************************/
void CFEMesh::CreateLineElementsFromMarkedEdges(CFEMesh*m_msh_ply, std::vector<long> &ele_vector_at_ply)
{
   CElem* m_ele = NULL;
   CEdge* m_edg = NULL;
   vec<CEdge*>ele_edges_vector(15);
   vec<CNode*>edge_nodes(3);
   long no_elements;

   // Create line elements
   std::vector<CEdge*> vct_used_edge;
   for (int i=0;i<(long)ele_vector_at_ply.size();i++)
   {
      m_ele = ele_vector[ele_vector_at_ply[i]];
      m_ele->GetEdges(ele_edges_vector);

      for(int j=0;j<(int)m_ele->GetEdgesNumber();j++)
      {
         m_edg = ele_edges_vector[j];
         if(!m_edg->GetMark())
         {
            continue;
         }
         bool done = false;
         for (int k=0;k<(int)vct_used_edge.size();k++)
         {
            if (vct_used_edge[k] == m_edg)
            {
               done = true;
               break;
            }
         }
         if (done)
         {
            continue;
         }
         else
         {
            vct_used_edge.push_back(m_edg);
         }

         if(m_msh_ply)
         {
            no_elements = (long)m_msh_ply->ele_vector.size();
         }
         else
         {
            no_elements = (long)ele_vector.size();
         }
         CElem* new_ele = new CElem(no_elements);
         new_ele->setElementProperties(MshElemType::LINE);
         m_edg->GetNodes(edge_nodes);
         new_ele->nodes_index[0] = edge_nodes[0]->GetIndex();
         new_ele->nodes_index[1] = edge_nodes[1]->GetIndex();
         new_ele->SetPatchIndex(m_ele->GetPatchIndex()+1);

         if(m_msh_ply)
         {
            m_msh_ply->ele_vector.push_back(new_ele);
         }
         else
         {
            ele_vector.push_back(new_ele);
         }
      }
   }
}


/**************************************************************************
MSHLib-Method:
Task:   HasSameCoordinateNode
Programing:
05/2007 NW implementation
**************************************************************************/
bool CFEMesh::HasSameCoordinatesNode(CNode* src_nod, long &nod_no)
{
   double dist,distmin;
   double pnt_src[3];
   double pnt_nod[3];
   long number;
   const double search_tolerance = 1.e-10;

   pnt_src[0] = src_nod->X();
   pnt_src[1] = src_nod->Y();
   pnt_src[2] = src_nod->Z();
   number = -1;
   distmin = 1.e+300;

   for(long i=0; i<(long)this->nod_vector.size(); i++)
   {
      pnt_nod[0] = nod_vector[i]->X();
      pnt_nod[1] = nod_vector[i]->Y();
      pnt_nod[2] = nod_vector[i]->Z();
      dist = EuklVek3dDist(pnt_src,pnt_nod);
      if(dist < distmin)
      {
         distmin = dist;
         number = i;
      }
   }

   if (distmin < search_tolerance)
   {
      nod_no = number;
      return true;
   }
   return false;
}
