/**************************************************************************
MSHLib - Object:
Task:
Programing:
08/2005 WW/OK Encapsulated from mshlib
**************************************************************************/

#include <iomanip>
#include <string>
#include <vector>

// MSHLib
#include "msh_node.h"

//========================================================================
namespace Mesh_Group
{
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   CNode::CNode(size_t Index, double x, double y, double z) :
   CCore(Index), epsilon (0.0), free_surface (-1), selected (0),
      patch_area (-1.0), crossroad (0), eqs_index (-1)
   {
      coordinate[0] = x;
      coordinate[1] = y;
      coordinate[2] = z;
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   10/2009 NW Implementation
   **************************************************************************/
   CNode::CNode(size_t Index, const CNode* parent) :
   CCore(Index), epsilon (0.0),
      free_surface (-1), selected (0), patch_area (-1.0), crossroad (0), eqs_index (-1)
   {
      coordinate[0] = parent->coordinate[0];
      coordinate[1] = parent->coordinate[1];
      coordinate[2] = parent->coordinate[2];
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CNode::operator = (const CNode& n)
   {
      boundary_type = n.boundary_type;
      index = n.index;
      mark = n.mark;
      eqs_index = n.eqs_index;
      coordinate[0] = n.coordinate[0];
      coordinate[1] = n.coordinate[1];
      coordinate[2] = n.coordinate[2];
   }
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   bool CNode::operator == (const CNode& n)
   {
      if(index == n.index)
         return true;
      else
         return false;
   }
   /**************************************************************************
   MSHLib-Method:
   06/2005 WW Implementation
   03/2006 OK patch_area
   **************************************************************************/
   void CNode::Write(std::ostream& osm) const
   {
      osm.setf(std::ios::scientific, std::ios::floatfield);
      std::string deli(" ");
      std::setw(14);
      osm.precision(14);
      osm << index << deli
         << coordinate[0] << deli
         << coordinate[1] << deli
         << coordinate[2] << deli;

      if(patch_area>0.0)
      {
         osm << "$AREA" << deli << patch_area;
      }
      osm << std::endl;
   }
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CNode::SetCoordinates(const double* argCoord)
   {
      coordinate[0] = argCoord[0];
      coordinate[1] = argCoord[1];
      coordinate[2] = argCoord[2];
   }

}                                                 // namespace Mesh_Group
