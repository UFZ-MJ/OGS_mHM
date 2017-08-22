/**************************************************************************
MSHLib - Object:
Task:
Programing:
08/2005 WW/OK Encapsulation from rf_ele_msh
last modified
**************************************************************************/

// C++
#include <cmath>
// MSHLib
#include "msh_edge.h"
#include "mathlib.h"                              //OK

//========================================================================
namespace Mesh_Group
{
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   03/2010 TF moved initialization of attributes to initialization list, added docu
   **************************************************************************/
   CEdge::CEdge(size_t Index, bool quadr)
      : CCore(Index), nodes_of_edges(3), joint (0), velocity (NULL)
   {
      quadratic = quadr;
      // Assume that each edge has three nodes, third node is middle point
      for(size_t i=0; i<3; i++) nodes_of_edges[i] = NULL;
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   CEdge::~CEdge()
   {
      nodes_of_edges.resize(0);
   }
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CEdge::operator = (CEdge& ed)
   {
      boundary_type = ed.boundary_type;
      index = ed.index;
      mark = ed.mark;
      for(size_t i=0; i<nodes_of_edges.Size(); i++)
         nodes_of_edges[i] = ed.nodes_of_edges[i];
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   double CEdge::Length()
   {
      double dx, dy, dz;
      dx = nodes_of_edges[1]->X()-nodes_of_edges[0]->X();
      dy = nodes_of_edges[1]->Y()-nodes_of_edges[0]->Y();
      dz = nodes_of_edges[1]->Z()-nodes_of_edges[0]->Z();
      return sqrt(dx*dx+dy*dy+dz*dz);
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   bool CEdge::operator == (CEdge& ed)
   {
      int identical;

      // Compare two ends
      identical=0;
      for(int i=0; i<2; i++)
      {
         if(nodes_of_edges[i] == ed.nodes_of_edges[i])
            identical++;
      }
      if(identical==2)
         return true;

      identical=0;
      for(int i=0; i<2; i++)
      {
         if(nodes_of_edges[1-i] == ed.nodes_of_edges[i])
            identical++;
      }
      if(identical==2)
         return true;

      return false;
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CEdge::Write(std::ostream& osm) const
   {
      osm<<"Edge: "<< index<<std::endl;
      for(size_t i=0; i<nodes_of_edges.Size(); i++)
      {
         osm<<"Node: "<< i<<std::endl;
         nodes_of_edges[i]->Write(osm);
      }
      osm<<std::endl;
   }

   /**************************************************************************
   MSHLib-Method:
   08/2006 OK Implementation
   **************************************************************************/
   void CEdge::SetNormalVector(double*ele_normal_vector,double*normal_vector)
   {
      // Element normal vector
      // Edge vector
      double edge_vector[3];
      GetEdgeVector(edge_vector);
      // Edge normal vector
      CrossProduction(edge_vector,ele_normal_vector,normal_vector);
      NormalizeVector(normal_vector,3);
   }

   /**************************************************************************
   MSHLib-Method:
   08/2006 OK Implementation
   **************************************************************************/
   void CEdge::GetEdgeVector(double*edge_vector)
   {
      // Edge vector
      //double edge_vector[3];
      edge_vector[0] = nodes_of_edges[1]->X()-nodes_of_edges[0]->X();
      edge_vector[1] = nodes_of_edges[1]->Y()-nodes_of_edges[0]->Y();
      edge_vector[2] = nodes_of_edges[1]->Z()-nodes_of_edges[0]->Z();
      //return edge_vector;
   }

   /**************************************************************************
   MSHLib-Method:
   08/2006 OK Implementation
   **************************************************************************/
   void CEdge::GetEdgeMidPoint(double*edge_vector)
   {
      // Edge vector
      //double edge_vector[3];
      edge_vector[0] = 0.5*(nodes_of_edges[1]->X()+nodes_of_edges[0]->X());
      edge_vector[1] = 0.5*(nodes_of_edges[1]->Y()+nodes_of_edges[0]->Y());
      edge_vector[2] = 0.5*(nodes_of_edges[1]->Z()+nodes_of_edges[0]->Z());
      //return edge_vector;
   }

}                                                 // namespace Mesh_Group
