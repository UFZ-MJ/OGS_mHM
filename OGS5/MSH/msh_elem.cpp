/**************************************************************************
MSHLib - Object:
Task:
Programing:
08/2005 WW/OK Encapsulation from rf_ele_msh
last modified
**************************************************************************/

#include <float.h>                                //WW
#include <cstdlib>                                //WW
#include <cmath>
#include "mathlib.h"
// MSHLib
//WW#include "MSHEnums.h" // KR 2010/11/15
#include "msh_elem.h"

#ifdef RFW_FRACTURE
#include "rf_pcs.h"                               //RFW 06/2006
#endif

#ifdef USE_TOKENBUF
#include "tokenbuf.h"

typedef struct
{
   char *name;
   int geo_type;
   int nnodes;
   int nnodesHQ;
   int ele_dim;
   int nfaces;
   int nedges;
} elem_descr_t;

elem_descr_t elem_descr[] =
{
   { "line",  1,  2,  3,  1,  2,  0 },
   { "quad",  2,  4,  9,  2,  4,  4 },
   { "hex",   3,  8, 20,  3,  6, 12 },
   { "tri",   4,  3,  6,  2,  3,  3 },
   { "tet",   5,  4, 10,  3,  4,  6 },
   { "pri",   6,  6, 15,  3,  5,  9 },
   { NULL,   -1, -1, -1, -1, -1, -1 }
};
#endif

//========================================================================
namespace Mesh_Group
{
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   CElem::CElem(size_t Index) : CCore(Index), normal_vector(NULL)
   {
      grid_adaptation = -1;
      nnodes = 0;
      nnodesHQ = 0;
      ele_dim = 1;                                // Dimension of element
      patch_index = 0;
      //
      volume = 0.0;
      face_index=-1;
      no_faces_on_surface = 0;
      owner = NULL;
      tranform_tensor = NULL;
      angle = NULL;
      gravity_center[0] = gravity_center[1] = gravity_center[2] = 0.0;
      normal_vector = NULL;
      area = 1.0;                                 //WW
   }
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   CElem::CElem():CCore(0), normal_vector(NULL)
   {
      selected = 0;
      matgroup_view = 0;
      grid_adaptation = -1;
      nnodes = 0;
      nnodesHQ = 0;
      ele_dim = 1;                                // Dimension of element
      patch_index = 0;
      //
      volume = 0.0;
      face_index=-1;
      no_faces_on_surface = 0;
      owner = NULL;
      nodes.resize(8);                            // Nodes of face
      tranform_tensor = NULL;
      angle = NULL;
      gravity_center[0] = gravity_center[1] = gravity_center[2] = 0.0;
      area = 1.0;                                 //WW area = 1.0
      normal_vector = NULL;
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   CElem::CElem(size_t Index, CElem* onwer, int Face) :
   CCore(Index), normal_vector(NULL), owner(onwer)
   {
      int i, j, k, n, ne;
      static int faceIndex_loc[10];
      static int edgeIndex_loc[10];
      no_faces_on_surface=0;
      n = owner->GetElementFaceNodes(Face, faceIndex_loc);
      face_index = Face;
      gravity_center[0] = gravity_center[1] = gravity_center[2] = 0.0;
      tranform_tensor = NULL;
      angle = NULL;
      normal_vector = NULL;
      //
      switch(owner->geo_type)
      {
         //case MshElemType::LINE:  // 1-D bar element //KR need not be processed
         case MshElemType::QUAD:                  // 2-D quadrilateral element
            this->setElementProperties(MshElemType::LINE, true);
            break;
         case MshElemType::HEXAHEDRON:            // 3-D hexahedral element
            this->setElementProperties(MshElemType::QUAD, true);
            break;
         case MshElemType::TRIANGLE:              // 2-D triagular element
            this->setElementProperties(MshElemType::LINE, true);
            break;
         case MshElemType::TETRAHEDRON:           // 3-D tetrahedral element
            this->setElementProperties(MshElemType::TRIANGLE, true);
            break;
         case MshElemType::PRISM:                 // 3-D prismatic element
            if (Face<2)                           // top or bottom face of the prism
               this->setElementProperties(MshElemType::TRIANGLE, true);
            else                                  // side of the prism
               this->setElementProperties(MshElemType::QUAD, true);
            break;
         default:
            std::cerr << "CElem::CElem MshElemType not handled" << std::endl;
      }

      patch_index =  owner->patch_index;
      quadratic = owner->quadratic;
      nodes_index.resize(n);
      nodes.resize(n);

      boundary_type='B';
      for(i=0; i<n; i++)
      {
         nodes_index[i] =
            owner->nodes_index[faceIndex_loc[i]];
         nodes[i] = owner->nodes[faceIndex_loc[i]];
                                                  //18.02.2009. cf. changes in mapping & generation. WW
         if((nodes[i]->boundary_type != '0')&&(nodes[i]->boundary_type != '1'))
            nodes[i]->boundary_type = 'B';
      }
      // Face edges
      ne = owner->GetEdgesNumber();
      edges.resize(nnodes);
      edges_orientation.resize(nnodes);
      edges_orientation = 1;
      for(i=0; i<nnodes; i++)
      {
         k = (i+1)%nnodes;
         for(j=0; j<ne; j++)
         {
            owner->GetLocalIndicesOfEdgeNodes(j, edgeIndex_loc);
            if( (faceIndex_loc[i]==edgeIndex_loc[0]&&
               faceIndex_loc[k]==edgeIndex_loc[1])||
               (faceIndex_loc[i]==edgeIndex_loc[1]&&
               faceIndex_loc[k]==edgeIndex_loc[0]) )
            {
               edges[i] = owner->edges[j];
               if(faceIndex_loc[i]==edgeIndex_loc[1]&&
                  faceIndex_loc[k]==edgeIndex_loc[0] )
                  edges_orientation[i] = -1;
               edges[i]->boundary_type = 'B';
               break;
            }
         }
      }
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW/OK Implementation
   **************************************************************************/
   CElem::CElem(size_t Index, CElem* m_ele_parent) :
   CCore(Index), normal_vector(NULL)
   {
      int i;
      //  static int faceIndex_loc[10];
      //  static int edgeIndex_loc[10];
      gravity_center[0] = gravity_center[1] = gravity_center[2] = 0.0;
      tranform_tensor = NULL;
      angle = NULL;
      this->setElementProperties(m_ele_parent->geo_type);

      patch_index =  m_ele_parent->patch_index;
      quadratic = m_ele_parent->quadratic;

      boundary_type='I';
      //----------------------------------------------------------------------
      // Initialize topological properties
      neighbors.resize(nfaces);
      for(i=0; i<nfaces; i++)
         neighbors[i] = NULL;
      edges.resize(nedges);
      edges_orientation.resize(nedges);
      for(i=0; i<nedges; i++)
      {
         edges[i] = NULL;
         edges_orientation[i] = 1;
      }
      area = 1.0;                                 //WW
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   CElem::~CElem()
   {
      nodes_index.resize(0);
      nodes.resize(0);
      edges.resize(0);
      neighbors.resize(0);
      mat_vector.resize(0);
      edges_orientation.resize(0);
      owner=NULL;
      if(tranform_tensor) delete tranform_tensor;
      tranform_tensor = NULL;
      if(angle) delete [] angle;
      angle = NULL;
      tranform_tensor = NULL;
      if (normal_vector)
         delete[] normal_vector;
      normal_vector = NULL;
   }
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CElem::FillTransformMatrix()
   {
      if(tranform_tensor) return;

      double xx[3];
      double yy[3];
      double zz[3];
      tranform_tensor = new Matrix(3,3);
      if(geo_type==MshElemType::LINE)
      {
         // x"_vec
         xx[0] = nodes[1]->X()-nodes[0]->X();
         xx[1] = nodes[1]->Y()-nodes[0]->Y();
         xx[2] = nodes[1]->Z()-nodes[0]->Z();
         NormalizeVector(xx,3);
         // an arbitrary vector
         for(size_t i=0; i<3; i++)
            yy[i] = 0.0;
                                                  //WW. 06.11.2007
         if(fabs(xx[0])>0.0 && fabs(xx[1])+fabs(xx[2]) < DBL_MIN)
            yy[2] = 1.0;
         else if (fabs(xx[1])>0.0 && fabs(xx[0])+fabs(xx[2])<DBL_MIN)
            yy[0] = 1.0;
         else if (fabs(xx[2])>0.0 && fabs(xx[0])+fabs(xx[1])<DBL_MIN)
            yy[1] = 1.0;
         else
         {
            for(size_t i=0; i<3; i++)
            {
               if(fabs(xx[i])>0.0)
               {
                  yy[i] = - xx[i];
                  break;
               }
            }
         }
         // z"_vec
         CrossProduction(xx,yy,zz);
         NormalizeVector(zz,3);
         // y"_vec
         CrossProduction(zz,xx,yy);
         NormalizeVector(yy,3);
      }
      else if (geo_type==MshElemType::QUAD || geo_type==MshElemType::TRIANGLE)
      {
         // x"_vec
         xx[0] = nodes[1]->X()-nodes[0]->X();
         xx[1] = nodes[1]->Y()-nodes[0]->Y();
         xx[2] = nodes[1]->Z()-nodes[0]->Z();
         NormalizeVector(xx,3);
         // a vector on the plane
         yy[0] = nodes[2]->X()-nodes[1]->X();
         yy[1] = nodes[2]->Y()-nodes[1]->Y();
         yy[2] = nodes[2]->Z()-nodes[1]->Z();
         // z"_vec. off plane
         CrossProduction(xx,yy,zz);
         NormalizeVector(zz,3);
         // y"_vec
         CrossProduction(zz,xx,yy);
         NormalizeVector(yy,3);
      }
      for(size_t i=0; i<3; i++)
      {
         (*tranform_tensor)(i,0) = xx[i];
         (*tranform_tensor)(i,1) = yy[i];
         (*tranform_tensor)(i,2) = zz[i];
      }
   }
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   02/2006 PCH Implementation
   **************************************************************************/
   double CElem::getTransformTensor(int idx)
   {
      //WW
      int i = idx%3;
      int j = idx/3;
      return (*tranform_tensor)(i,j);
      //return MatT[idx];
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CElem::SetFace()
   {
      nodes.resize(8);
      nodes_index.resize(8);
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CElem:: SetFace(CElem* onwer, const int Face)
   {
      int i, n;
      static int nodeIndex_loc[8];
      no_faces_on_surface=0;
      owner = onwer;
      n = owner->GetElementFaceNodes(Face, nodeIndex_loc);
      face_index = Face;
      switch(owner->geo_type)
      {
         //case MshElemType::LINE:  // 1-D bar element
         //case MshElemType::QUAD: // 2-D quadrilateral element
         case MshElemType::HEXAHEDRON:            // 3-D hexahedral element
            this->setElementProperties(MshElemType::QUAD, true);
            break;
            //case MshElemType::TRIANGLE:  // 2-D triagular element
         case MshElemType::TETRAHEDRON:           // 3-D tetrahedral element
            this->setElementProperties(MshElemType::TRIANGLE, true);
            break;
         case MshElemType::PRISM:
            if(Face<2)
               this->setElementProperties(MshElemType::TRIANGLE, true);
            else
               this->setElementProperties(MshElemType::QUAD, true);
            break;                                // 3-D prismatic element
         default:
            std::cerr << "CElem::SetFace MshElemType not handled" << std::endl;
      }

      for(i=0; i<n; i++)
      {
         nodes[i] =
            owner->nodes[nodeIndex_loc[i]];
      }
   }
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   11/2010 KR moved to MSHEnums.h
   **************************************************************************/
   std::string CElem::GetName() const
   {
      return MshElemType2String(geo_type);
   }

/**************************************************************************
 MSHLib-Method:
 Programing:
 06/2005 WW Implementation
 08/2005 WW/OK Extension for GMS/SOL files
 10/2008 OK FEFLOW
 **************************************************************************/
void CElem::Read(std::istream& is, int fileType)
{
	//fileType=0: msh
	//fileType=1: rfi
	//fileType=2: gmsh
	//fileType=3: GMS
	//fileType=4: SOL
	//fileType=5: FLAC 3D. WW
	//fileType=5x: FLAC3D //MR
	//fileType=6: FEFLOW //OK
	//fileType=7: GMSH 2008 Version//TK
	int idummy, et;
	std::string buffer, name;
	idummy = et = -1;
	int j = 0;
	int gmsh_patch_index; //TKOK
	int nb_tags; //TK

	//   is.ignore(numeric_limits<int>::max(), '\n');
	//----------------------------------------------------------------------
	// 1 Reading element type data
	switch (fileType) {
	//....................................................................
	case 0: // msh
		is >> index >> patch_index;
		is >> buffer;

		if (buffer.find("-1") != std::string::npos) {
			grid_adaptation = strtol(buffer.data(), NULL, 0);
			is >> name;
		} else
			name = buffer;

		geo_type = String2MshElemType(name);
		break;
		//....................................................................
	case 1: // rfi
		is >> index >> patch_index >> name;
		geo_type = String2MshElemType(name);
		break;
		//....................................................................
	case 2: // gmsh
		is >> index >> et >> gmsh_patch_index >> idummy >> nnodes;
		patch_index = gmsh_patch_index - 1; //OK
		switch (et) {
		case 1:
			geo_type = MshElemType::LINE;
			break;
		case 2:
			geo_type = MshElemType::TRIANGLE;
			break;
		case 3:
			geo_type = MshElemType::QUAD;
			break;
		case 4:
			geo_type = MshElemType::TETRAHEDRON;
			break;
		case 5:
			geo_type = MshElemType::HEXAHEDRON;
			break;
		case 6:
			geo_type = MshElemType::PRISM;
			break;
		}
		index--;
		break;
	case 7: // GMSH 2008
		is >> index >> et >> nb_tags >> gmsh_patch_index;
		patch_index = gmsh_patch_index -1;
		for (j = 0; j < nb_tags - 1; j++) {
			is >> idummy;
		}
		switch (et) {
		case 1:
			geo_type = MshElemType::LINE;
			nnodes = 2;
			break;
		case 2:
			geo_type = MshElemType::TRIANGLE;
			nnodes = 3;
			break;
		case 3:
			geo_type = MshElemType::QUAD;
			nnodes = 4;
			break;
		case 4:
			geo_type = MshElemType::TETRAHEDRON;
			nnodes = 4;
			break;
		case 5:
			geo_type = MshElemType::HEXAHEDRON;
			nnodes = 8;
			break;
		case 6:
			geo_type = MshElemType::PRISM;
			nnodes = 6;
			break;
		case 15:
			geo_type = MshElemType::INVALID;
			nnodes = 1;
			break;
		default:
			geo_type = MshElemType::INVALID;
		}
		index--;
		break;
		//....................................................................
	case 3: // GMS
		geo_type = MshElemType::TRIANGLE;
		break;
		//....................................................................
	case 4: // gmsh
		geo_type = MshElemType::TRIANGLE;
		break;
		//....................................................................
	case 5: // FLAC 3D. 14.01.2008 WW
		geo_type = MshElemType::HEXAHEDRON;
		break;
	case 56: // FLAC3D - pri (Wedge)        //MR
		geo_type = MshElemType::PRISM; //MR
		fileType = 5; //MR
		break; //MR
	}

	// TF
	if (geo_type == MshElemType::INVALID) {
		// read rest of line
		std::string tmp;
		getline(is, tmp);
		return;
	}

	//----------------------------------------------------------------------
	// 2 Element configuration
	this->setElementProperties(geo_type);

	//----------------------------------------------------------------------
	// 3 Reading element node data
	switch (fileType) {
	case 0: // msh
		for (int i = 0; i < nnodes; i++)
			is >> nodes_index[i];
		break;
	case 1: // rfi
		for (int i = 0; i < nnodes; i++)
			is >> nodes_index[i];
		break;
	case 2: // gmsh
		for (int i = 0; i < nnodes; i++) {
			is >> nodes_index[i];
			nodes_index[i] -= 1;
		}
		break;
	case 7: // GMSH 2008
		if (et != 15) {
			for (int i = 0; i < nnodes; i++) {
				is >> nodes_index[i];
				nodes_index[i] -= 1;
			}
		} else {
			// eat rest of line
			std::string dummy;
			is >> dummy;
		}
		break;
	case 3: // GMS
		for (int i = 0; i < nnodes; i++) {
			is >> nodes_index[i];
			nodes_index[i] -= 1;
		}
		break;
	case 4: // SOL
		for (int i = 0; i < nnodes; i++) {
			is >> nodes_index[i];
			nodes_index[i] -= 1;
		}
		is >> patch_index;
		break;
	case 5: // FLAC 3D. 14.01.2008. WW
		for (int i = 0; i < nnodes; i++) {
			is >> nodes_index[i];
			nodes_index[i] -= 1;
		}
		break;
	case 6: // FEFLOLW
		for (int i = 0; i < nnodes; i++) {
			is >> nodes_index[i];
			nodes_index[i] -= 1;
		}
		break;
	case 8: // GMS_3DM
		is >> idummy;
		for (int i = 0; i < nnodes; i++) {
			is >> nodes_index[i];
			nodes_index[i] -= 1;
		}
		is >> patch_index;
		patch_index -= 1;
		break;
	}
	is >> std::ws;
	//----------------------------------------------------------------------
	// Initialize topological properties
	neighbors.resize(nfaces);
	for (int i = 0; i < nfaces; i++)
		neighbors[i] = NULL;
	edges.resize(nedges);
	edges_orientation.resize(nedges);
	for (int i = 0; i < nedges; i++) {
		edges[i] = NULL;
		edges_orientation[i] = 1;
	}
}

   void CElem::InitializeMembers()
   {
      // Initialize topological properties
      neighbors.resize(nfaces);
      for(int i=0; i<nfaces; i++)
         neighbors[i] = NULL;
      edges.resize(nedges);
      edges_orientation.resize(nedges);
      for(int i=0; i<nedges; i++)
      {
         edges[i] = NULL;
         edges_orientation[i] = 1;
      }
   }

#ifdef USE_TOKENBUF
   void CElem::Read(TokenBuf* tokenbuf, int fileType)
   {
      //fileType=0: msh
      //fileType=1: rfi
      //fileType=2: gmsh
      //fileType=3: GMS
      //fileType=4: SOL
      int idummy, et;
      string buffer, name;
      idummy=et=-1;
      char line_buf[LINE_MAX];
      char str1[LINE_MAX], remains[LINE_MAX];

      //   is.ignore(numeric_limits<int>::max(), '\n');
      //----------------------------------------------------------------------
      // 1 Reading element type data
      switch(fileType)
      {
         //....................................................................
         case 0:                                  // msh
            tokenbuf->get_non_empty_line(line_buf, LINE_MAX);
            sscanf(line_buf, "%ld %ld %s %[0-9a-zA-Z ]", &index, &patch_index, str1, remains);
            if(!strcmp(str1, "-1"))
            {
               grid_adaptation = strtol(str1, NULL, 0);
               sscanf(remains, "%s %[0-9a-zA-Z ]", str1, remains);
               name = std::string(str1);
            }
            else
               name = std::string(str1);

            for(int i=0; elem_descr[i].name != NULL; i++)
            {
               if(!strcmp(str1, elem_descr[i].name))
               {
                  geo_type = elem_descr[i].geo_type;
                  nnodes   = elem_descr[i].nnodes;
                  nnodesHQ = elem_descr[i].nnodesHQ;
                  ele_dim  = elem_descr[i].ele_dim;
                  nfaces   = elem_descr[i].nfaces;
                  nedges   = elem_descr[i].nedges;
                  nodes_index.resize(nnodes);
                  for(int j=0; j<nnodes; j++)
                  {
                     sscanf(remains, "%ld %[0-9] ", &nodes_index[j], remains);
                  }
                  break;
               }
            }

            break;
            //....................................................................
      }
      //----------------------------------------------------------------------
      // Initialize topological properties
      neighbors.resize(nfaces);
      for(int i=0; i<nfaces; i++)
         neighbors[i] = NULL;
      edges.resize(nedges);
      edges_orientation.resize(nedges);
      for(int i=0; i<nedges; i++)
      {
         edges[i] = NULL;
         edges_orientation[i] = 1;
      }
   }
#endif

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CElem::WriteIndex(std::ostream &os) const
   {
      //Comment for GUI WW if(quadratic) nn = nnodesHQ;
      os << index << "  " << patch_index << "  " << GetName() << "  ";
      for(int i=0; i<nnodes; i++)
         os << nodes_index[i] << "  ";
      os << std::endl;
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CElem::WriteIndex_TEC(std::ostream &os) const
   {
      std::string deli = "  ";
      if(geo_type==MshElemType::LINE)
         os<<nodes_index[0]+1<<deli<<nodes_index[1]+1<<deli<<nodes_index[1]+1<<deli<<nodes_index[0]+1;

      else if(geo_type==MshElemType::TRIANGLE)
      {
         os<<nodes_index[0]+1<<deli<<nodes_index[1]+1<<deli<<nodes_index[2]+1<<deli<<nodes_index[0]+1;
      }
      else if(geo_type==MshElemType::PRISM)
      {
         os<<nodes_index[0]+1<<deli<<nodes_index[0]+1<<deli<<nodes_index[1]+1<<deli<<nodes_index[2]+1<<deli
            <<nodes_index[3]+1<<deli<<nodes_index[3]+1<<deli<<nodes_index[4]+1<<deli<<nodes_index[5]+1<<deli;
      }
      else
      {

         for(int i=0; i<nnodes; i++)
            os<<nodes_index[i]+1<<deli;
      }
      os<<'\n';                                   // GK44 for io performance do not flush buffer with endl...
   }
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CElem::WriteAll(std::ostream &os) const
   {
      std::string deli = "  ";
      os<<index<<deli<<patch_index<<deli<<GetName()<<deli;
      //if(index==0)
      os<<"Index X Y Z: "<<std::endl;
      for(size_t i=0; i<nodes.Size(); i++)
      {
         os<<nodes_index[i]
            <<deli<<nodes[i]->X()
            <<deli<<nodes[i]->Y()
            <<deli<<nodes[i]->Z()<<std::endl;
      }

   }
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CElem::WriteNeighbors(std::ostream &os) const
   {
      os<<"Neighbors of "<<index<<std::endl;
      for(int i=0; i<nfaces; i++)
         neighbors[i]->WriteAll(os);
      os<<"End neighbors of "<<index<<std::endl<<std::endl;;
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CElem:: MarkingAll(bool makop)
   {
      int i;
      mark=makop;
      int SizeV = nnodes;
      if(quadratic) SizeV = nnodesHQ;
      for (i=0; i< SizeV;i++)
         nodes[i]->mark = mark;
      for (i=0; i< nedges;i++)
         edges[i]->mark = mark;
   }
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CElem::SetNodes(vec<CNode*>&  ele_nodes, bool ReSize)
   {
      int SizeV = nnodes;
      if(quadratic) SizeV = nnodesHQ;
      if(ReSize)
      {
         nodes.resize(SizeV);
         nodes_index.resize(SizeV);
      }
      for (int i=0; i< SizeV;i++)
      {
         nodes[i] = ele_nodes[i];
         nodes_index[i] = nodes[i]->GetIndex();
      }
   }
   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void  CElem::GetLocalIndicesOfEdgeNodes(const int Edge, int *EdgeNodes)
   {
      switch(geo_type)
      {
         case MshElemType::LINE:
            break;                                // 1-D bar element
         case MshElemType::QUAD:                  // 2-D quadrilateral element
            EdgeNodes[0] = Edge;
            EdgeNodes[1] = (Edge+1)%4;
            break;
         case MshElemType::HEXAHEDRON:            // 3-D hexahedral element
            if(Edge<8)
            {
               EdgeNodes[0] = Edge;
               EdgeNodes[1] = (Edge+1)%4+4*(int)(Edge/4);
            }
            else
            {
               EdgeNodes[0] = Edge%4;
               EdgeNodes[1] = Edge%4+4;
            }
            break;
         case MshElemType::TRIANGLE:              // 2-D triagular element
            EdgeNodes[0] = Edge;
            EdgeNodes[1] = (Edge+1)%3;
            break;
         case MshElemType::TETRAHEDRON:           // 3-D tetrahedra
            if(Edge<3)
            {
               EdgeNodes[0] = Edge;
               EdgeNodes[1] = (Edge+1)%3;
            }
            else
            {
               EdgeNodes[0] = 3;
               EdgeNodes[1] = (Edge+1)%3;
            }

            break;
         case MshElemType::PRISM:                 // 3-D prismatic element
            if(Edge<6)
            {
               EdgeNodes[0] = Edge;
               EdgeNodes[1] = (Edge+1)%3+3*(int)(Edge/3);
            }
            else
            {
               EdgeNodes[0] = Edge%3;
               EdgeNodes[1] = Edge%3+3;
            }
            break;
         default:
            std::cerr << "CElem::GetLocalIndicesOfEdgeNodes MshElemType not handled" << std::endl;
      }
   }
   /**************************************************************************
   GetElementFaceNodes
   Task: Get local indeces of an element face nodes
   Return: number of nodes of a face
   Programing:
   06/2004 WW
   **************************************************************************/
   int CElem::GetElementFaces1D(int *FaceNode)
   {
      FaceNode[0] = 0;
      FaceNode[1] = 1;
      return 2;
   }
   /**************************************************************************
   GetElementFaceNodesTri
   Task: Get local indeces of a traingle element face nodes
   Augs.:
           const int Face :  Local index of element face
           const int order:  1 Linear. 2, quadratic
           int *FaceNode  :  Local index of face nodes
   Return: number of nodes of a face
   Programing:
   09/2004 WW
   **************************************************************************/
   int CElem::GetElementFacesTri(int Face, int *FaceNode)
   {
      if(!quadratic)
      {
         FaceNode[0] = Face;
         FaceNode[1] = (Face+1)%3;
         return 2;
      }
      else
      {
         FaceNode[0] = Face;
         FaceNode[1] = (Face+1)%3;
         FaceNode[2] = Face+3;
         return 3;
      }
   }
   /**************************************************************************
   GetElementFaceNodesQuad
   Task: Get local indeces of a quadralateral element face nodes
   Augs.:
           const int Face :  Local index of element face
           const int order:  1 Linear. 2, quadratic
           int *FaceNode  :  Local index of face nodes
   Return: number of nodes of a face
   Programing:
   09/2004 WW
   **************************************************************************/
   int CElem::GetElementFacesQuad(int Face, int *FaceNode)
   {
      if(!quadratic)
      {
         FaceNode[0] = Face;
         FaceNode[1] = (Face+1)%4;
         return 2;
      }
      else
      {
         FaceNode[0] = Face;
         FaceNode[1] = (Face+1)%4;
         FaceNode[2] = Face+4;
         return 3;
      }
   }
   /**************************************************************************
   GetElementFaceNodesHex
   Task: Get local indeces of a hexahedra element face nodes
   Augs.:
           const int Face :  Local index of element face
           const int order:  1 Linear. 2, quadratic
           int *FaceNode  :  Local index of face nodes
   Return: number of nodes of a face
   Programing:
   09/2004 WW
   **************************************************************************/
   int CElem::GetElementFacesHex(int Face, int *FaceNode)
   {
      int nn=4, k = 0;
      if(quadratic) nn = 8;
      switch(Face)
      {
         case 0:
            for(k=0; k<4; k++)
               FaceNode[k] = k;
            if(quadratic)
            {
               for(k=0; k<4; k++)
                  FaceNode[k+4] = k+8;
            }
            break;
         case 1:
            for(k=0; k<4; k++)
               FaceNode[k] = k+4;
            if(quadratic)
            {
               for(k=0; k<4; k++)
                  FaceNode[k+4] = k+12;
            }
            break;
         case 2:
            FaceNode[0] = 0;
            FaceNode[1] = 4;
            FaceNode[2] = 5;
            FaceNode[3] = 1;
            if(quadratic)
            {
               FaceNode[4] = 16;
               FaceNode[5] = 12;
               FaceNode[6] = 17;
               FaceNode[7] = 8;
            }
            break;
         case 3:
            FaceNode[0] = 1;
            FaceNode[1] = 5;
            FaceNode[2] = 6;
            FaceNode[3] = 2;
            if(quadratic)
            {
               FaceNode[4] = 17;
               FaceNode[5] = 13;
               FaceNode[6] = 18;
               FaceNode[7] = 9;
            }

            break;
         case 4:
            FaceNode[0] = 2;
            FaceNode[1] = 6;
            FaceNode[2] = 7;
            FaceNode[3] = 3;
            if(quadratic)
            {
               FaceNode[4] = 18;
               FaceNode[5] = 14;
               FaceNode[6] = 19;
               FaceNode[7] = 10;
            }
            break;
         case 5:
            FaceNode[0] = 0;
            FaceNode[1] = 3;
            FaceNode[2] = 7;
            FaceNode[3] = 4;
            if(quadratic)
            {
               FaceNode[4] = 11;
               FaceNode[5] = 19;
               FaceNode[6] = 15;
               FaceNode[7] = 16;
            }
            break;
      }
      return nn;
   }
   /**************************************************************************
   GetElementFaceNodesTet
   Task: Get local indeces of a Tedrahedra element face nodes
   Augs.:
           const int Face :  Local index of element face
           const int order:  1 Linear. 2, quadratic
           int *FaceNode  :  Local index of face nodes
   Return: number of nodes of a face
   Programing:
   09/2004 WW
   **************************************************************************/
   int CElem::GetElementFacesTet(int Face, int *FaceNode)
   {
      int nn=3;
      if(quadratic) nn =6;
      switch(Face)
      {
         case 0:
            FaceNode[0] = 1;
            FaceNode[1] = 2;
            FaceNode[2] = 3;
            if(quadratic)
            {
               FaceNode[3] = 5 ;
               FaceNode[4] = 8;
               FaceNode[5] = 7;
            }
            break;
         case 1:
            FaceNode[0] = 3;
            FaceNode[1] = 2;
            FaceNode[2] = 0;
            if(quadratic)
            {
               FaceNode[3] = 8 ;
               FaceNode[4] = 6;
               FaceNode[5] = 9;
            }
            break;
         case 2:
            FaceNode[0] = 1;
            FaceNode[1] = 3;
            FaceNode[2] = 0;
            if(quadratic)
            {
               FaceNode[3] = 7 ;
               FaceNode[4] = 9;
               FaceNode[5] = 4;
            }
            break;
         case 3:
            FaceNode[0] = 0;
            FaceNode[1] = 2;
            FaceNode[2] = 1;
            if(quadratic)
            {
               FaceNode[3] = 6 ;
               FaceNode[4] = 5;
               FaceNode[5] = 4;
            }
            break;

      }
      return nn;
   }

   /**************************************************************************
   GetElementFaceNodesPri
   Task: Get local indeces of a prismal element face nodes
   Augs.:
           int Face :  Local index of element face
           int *FaceNode  :  Local index of face nodes
   Return: number of nodes of a face
   Programing:
   09/2004 WW
   **************************************************************************/
   int CElem::GetElementFacesPri(int Face, int *FaceNode)
   {
      int nn=3, k = 0;
      switch(Face)
      {
         case 0:
            nn = 3;
            for(k=0; k<3; k++)
               FaceNode[k] = k;
            if(quadratic)
            {
               for(k=0; k<3; k++)
                  FaceNode[k+3] = k+6;
               nn = 6;
            }
            break;
         case 1:
            for(k=0; k<3; k++)
               FaceNode[k] = k+3;
            nn = 3;
            if(quadratic)
            {
               for(k=0; k<3; k++)
                  FaceNode[k+3] = k+9;
               nn = 6;
            }
            break;
         case 2:
            FaceNode[0] = 1;
            FaceNode[1] = 2;
            FaceNode[2] = 5;
            FaceNode[3] = 4;
            nn = 4;
            if(quadratic)
            {
               FaceNode[4] = 7 ;
               FaceNode[5] = 14;
               FaceNode[6] = 10;
               FaceNode[7] = 13;
               nn = 8;
            }
            break;
         case 3:
            FaceNode[0] = 5;
            FaceNode[1] = 2;
            FaceNode[2] = 0;
            FaceNode[3] = 3;
            nn = 4;
            if(quadratic)
            {
               FaceNode[4] = 14 ;
               FaceNode[5] =  8;
               FaceNode[6] = 12;
               FaceNode[7] = 10;
               nn = 8;
            }
            break;
         case 4:
            FaceNode[0] = 0;
            FaceNode[1] = 1;
            FaceNode[2] = 4;
            FaceNode[3] = 3;
            nn = 4;
            if(quadratic)
            {
               FaceNode[4] = 6 ;
               FaceNode[5] = 13;
               FaceNode[6] = 9;
               FaceNode[7] = 12;
               nn = 8;
            }
            break;

      }
      return nn;
   }
   /**************************************************************************
   GetElementFaces
   Task: set element faces (Geometry)
   Augs.:
           const int Face :  Local index of element face
           const int order:  1 Linear. 2, quadratic
           int *FaceNode  :  Local index of face nodes

   Programing:
   09/2004 WW
   **************************************************************************/
   int CElem::GetElementFaceNodes(int Face, int *FacesNode)
   {
      switch(geo_type)
      {
         case MshElemType::LINE:                  // 1-D bar element
            return GetElementFaces1D(FacesNode);
         case MshElemType::QUAD:                  // 2-D quadrilateral element
            return GetElementFacesQuad(Face, FacesNode);
         case MshElemType::HEXAHEDRON:            // 3-D hexahedral element
            return GetElementFacesHex(Face, FacesNode);
         case MshElemType::TRIANGLE:              // 2-D triagular element
            return GetElementFacesTri(Face, FacesNode);
         case MshElemType::TETRAHEDRON:           // 3-D tetrahedral element
            return GetElementFacesTet(Face, FacesNode);
         case MshElemType::PRISM:
            return GetElementFacesPri(Face, FacesNode);
            // 3-D prismatic element
         default:
            std::cerr << "CElem::GetElementFaceNodes MshElemType not handled" << std::endl;
      }
      return 0;
   }

   /**************************************************************************
   FindFaceEdges(const int LocalFaceIndex, vec<CEdge*>&  face_edges)
   Task: set element faces (Geometry)
   Augs.:
           const int LocalFaceIndex :  Local index of element face
           vec<CEdge*>&  face_edges :  Found face edges
   retrun: number of face edges

   Programing:
   07/2005 WW
   **************************************************************************/
   /* KR not used
   int CElem::FindFaceEdges(const int LocalFaceIndex, vec<CEdge*>&  face_edges)
   {
       int i;
      i=0;
      switch(geo_type)
      {
         case MshElemType::LINE:
              break; // 1-D bar element
         case MshElemType::QUAD: // 2-D quadrilateral element
           face_edges[0] = edges[LocalFaceIndex];
   return 1;
   break;
   case MshElemType::HEXAHEDRON: // 3-D hexahedral element
   if(LocalFaceIndex<2)
   {
   for(i=0; i<4; i++)
   face_edges[i] = edges[LocalFaceIndex*4+i];
   }
   else if(LocalFaceIndex==2)
   {
   face_edges[0] = edges[0];
   face_edges[1] = edges[9];
   face_edges[2] = edges[4];
   face_edges[3] = edges[8];
   }
   else if(LocalFaceIndex==3)
   {
   face_edges[0] = edges[1];
   face_edges[1] = edges[9];
   face_edges[2] = edges[5];
   face_edges[3] = edges[10];
   }
   else if(LocalFaceIndex==4)
   {
   face_edges[0] = edges[2];
   face_edges[1] = edges[10];
   face_edges[2] = edges[6];
   face_edges[3] = edges[11];
   }
   else if(LocalFaceIndex==5)
   {
   face_edges[0] = edges[7];
   face_edges[1] = edges[8];
   face_edges[2] = edges[3];
   face_edges[3] = edges[11];
   }

   return 4;
   break;
   case MshElemType::TRIANGLE:  // 2-D triagular element
   face_edges[0] = edges[LocalFaceIndex];
   return 1;
   break;
   case MshElemType::TETRAHEDRON:  // 3-D tetrahedra
   if(LocalFaceIndex==0)
   {
   face_edges[0] = edges[1];
   face_edges[1] = edges[4];
   face_edges[2] = edges[3];
   }
   else if(LocalFaceIndex==1)
   {
   face_edges[0] = edges[2];
   face_edges[1] = edges[4];
   face_edges[2] = edges[5];
   }
   else if(LocalFaceIndex==2)
   {
   face_edges[0] = edges[0];
   face_edges[1] = edges[3];
   face_edges[2] = edges[5];
   }
   else if(LocalFaceIndex==3)
   {
   face_edges[0] = edges[0];
   face_edges[1] = edges[1];
   face_edges[2] = edges[2];
   }
   return 3;
   break;
   case MshElemType::PRISM: // 3-D prismatic element
   if(LocalFaceIndex<2)
   {
   for(i=0; i<3; i++)
   face_edges[i] = edges[LocalFaceIndex*3+i];
   return 3;
   }
   else if(LocalFaceIndex==2)
   {
   face_edges[0] = edges[1];
   face_edges[1] = edges[8];
   face_edges[2] = edges[4];
   face_edges[3] = edges[9];
   return 4;
   }
   else if(LocalFaceIndex==3)
   {
   face_edges[0] = edges[2];
   face_edges[1] = edges[9];
   face_edges[2] = edges[5];
   face_edges[3] = edges[7];
   return 4;
   }
   else if(LocalFaceIndex==4)
   {
   face_edges[0] = edges[0];
   face_edges[1] = edges[8];
   face_edges[2] = edges[3];
   face_edges[3] = edges[7];
   return 4;
   }

   break;
   }
   return 0;
   }
   */

   /**************************************************************************
   MSHLib-Method:
   Task:´For elements with straight edges and surfaces
   Programing:
   06/2005 WW Implementation
   **************************************************************************/
   void CElem::ComputeVolume()
   {
      double x1buff[3];
      double x2buff[3];
      double x3buff[3];
      double x4buff[3];
      volume = 0.0;

      if(geo_type!=MshElemType::LINE)
      {
         x1buff[0] = nodes[0]->X();
         x1buff[1] = nodes[0]->Y();
         x1buff[2] = nodes[0]->Z();

         x2buff[0] = nodes[1]->X();
         x2buff[1] = nodes[1]->Y();
         x2buff[2] = nodes[1]->Z();

         x3buff[0] = nodes[2]->X();
         x3buff[1] = nodes[2]->Y();
         x3buff[2] = nodes[2]->Z();
      }

      switch(geo_type)
      {
         case MshElemType::LINE:                  // Line
            x2buff[0] = nodes[nnodes-1]->X()-nodes[0]->X();
            x2buff[1] = nodes[nnodes-1]->Y()-nodes[0]->Y();
            x2buff[2] = nodes[nnodes-1]->Z()-nodes[0]->Z();
            volume = sqrt(x2buff[0]*x2buff[0]+x2buff[1]*x2buff[1]+x2buff[2]*x2buff[2]);
                                                  //CMCD kg44 reactivated
            representative_length = sqrt(x2buff[0]*x2buff[0]+x2buff[1]*x2buff[1]+x2buff[2]*x2buff[2]) ;
            break;
         case MshElemType::TRIANGLE:              // Triangle
            volume = ComputeDetTri(x1buff, x2buff, x3buff);
                                                  //kg44 reactivated
            representative_length = sqrt(volume)*4.0;
            break;
         case MshElemType::QUAD:                  // Quadrilateral
            x4buff[0] = nodes[3]->X();
            x4buff[1] = nodes[3]->Y();
            x4buff[2] = nodes[3]->Z();

            volume =  ComputeDetTri(x1buff, x2buff, x3buff)
               +ComputeDetTri(x3buff, x4buff, x1buff);
            representative_length = sqrt(volume); //kg44 reactivated
            break;
         case MshElemType::TETRAHEDRON:           // Tedrahedra
            x4buff[0] = nodes[3]->X();
            x4buff[1] = nodes[3]->Y();
            x4buff[2] = nodes[3]->Z();

            volume =  ComputeDetTex(x1buff, x2buff, x3buff, x4buff);
                                                  //kg44 reactivated
            representative_length = sqrt(volume)*6.0;
            break;
         case MshElemType::HEXAHEDRON:            // Hexehadra
            x1buff[0] = nodes[4]->X();
            x1buff[1] = nodes[4]->Y();
            x1buff[2] = nodes[4]->Z();

            x2buff[0] = nodes[7]->X();
            x2buff[1] = nodes[7]->Y();
            x2buff[2] = nodes[7]->Z();

            x3buff[0] = nodes[5]->X();
            x3buff[1] = nodes[5]->Y();
            x3buff[2] = nodes[5]->Z();

            x4buff[0] = nodes[0]->X();
            x4buff[1] = nodes[0]->Y();
            x4buff[2] = nodes[0]->Z();
            volume  = ComputeDetTex(x1buff, x2buff, x3buff, x4buff);

            x1buff[0] = nodes[5]->X();
            x1buff[1] = nodes[5]->Y();
            x1buff[2] = nodes[5]->Z();

            x2buff[0] = nodes[3]->X();
            x2buff[1] = nodes[3]->Y();
            x2buff[2] = nodes[3]->Z();

            x3buff[0] = nodes[1]->X();
            x3buff[1] = nodes[1]->Y();
            x3buff[2] = nodes[1]->Z();

            x4buff[0] = nodes[0]->X();
            x4buff[1] = nodes[0]->Y();
            x4buff[2] = nodes[0]->Z();
            volume  += ComputeDetTex(x1buff, x2buff, x3buff, x4buff);

            x1buff[0] = nodes[5]->X();
            x1buff[1] = nodes[5]->Y();
            x1buff[2] = nodes[5]->Z();

            x2buff[0] = nodes[7]->X();
            x2buff[1] = nodes[7]->Y();
            x2buff[2] = nodes[7]->Z();

            x3buff[0] = nodes[3]->X();
            x3buff[1] = nodes[3]->Y();
            x3buff[2] = nodes[3]->Z();

            x4buff[0] = nodes[0]->X();
            x4buff[1] = nodes[0]->Y();
            x4buff[2] = nodes[0]->Z();
            volume  += ComputeDetTex(x1buff, x2buff, x3buff, x4buff);

            x1buff[0] = nodes[5]->X();
            x1buff[1] = nodes[5]->Y();
            x1buff[2] = nodes[5]->Z();

            x2buff[0] = nodes[7]->X();
            x2buff[1] = nodes[7]->Y();
            x2buff[2] = nodes[7]->Z();

            x3buff[0] = nodes[6]->X();
            x3buff[1] = nodes[6]->Y();
            x3buff[2] = nodes[6]->Z();

            x4buff[0] = nodes[2]->X();
            x4buff[1] = nodes[2]->Y();
            x4buff[2] = nodes[2]->Z();
            volume  += ComputeDetTex(x1buff, x2buff, x3buff, x4buff);

            x1buff[0] = nodes[1]->X();
            x1buff[1] = nodes[1]->Y();
            x1buff[2] = nodes[1]->Z();

            x2buff[0] = nodes[3]->X();
            x2buff[1] = nodes[3]->Y();
            x2buff[2] = nodes[3]->Z();

            x3buff[0] = nodes[5]->X();
            x3buff[1] = nodes[5]->Y();
            x3buff[2] = nodes[5]->Z();

            x4buff[0] = nodes[2]->X();
            x4buff[1] = nodes[2]->Y();
            x4buff[2] = nodes[2]->Z();
            volume  += ComputeDetTex(x1buff, x2buff, x3buff, x4buff);

            x1buff[0] = nodes[3]->X();
            x1buff[1] = nodes[3]->Y();
            x1buff[2] = nodes[3]->Z();

            x2buff[0] = nodes[7]->X();
            x2buff[1] = nodes[7]->Y();
            x2buff[2] = nodes[7]->Z();

            x3buff[0] = nodes[5]->X();
            x3buff[1] = nodes[5]->Y();
            x3buff[2] = nodes[5]->Z();

            x4buff[0] = nodes[2]->X();
            x4buff[1] = nodes[2]->Y();
            x4buff[2] = nodes[2]->Z();

            volume  += ComputeDetTex(x1buff, x2buff, x3buff, x4buff);
                                                  //kg44 reactivated
            representative_length = pow(volume,1./3.);
            break;
         case MshElemType::PRISM:                 // Prism
            x4buff[0] = nodes[3]->X();
            x4buff[1] = nodes[3]->Y();
            x4buff[2] = nodes[3]->Z();
            volume =  ComputeDetTex(x1buff, x2buff, x3buff, x4buff);

            x1buff[0] = nodes[1]->X();
            x1buff[1] = nodes[1]->Y();
            x1buff[2] = nodes[1]->Z();

            x2buff[0] = nodes[4]->X();
            x2buff[1] = nodes[4]->Y();
            x2buff[2] = nodes[4]->Z();

            x3buff[0] = nodes[2]->X();
            x3buff[1] = nodes[2]->Y();
            x3buff[2] = nodes[2]->Z();

            x4buff[0] = nodes[3]->X();
            x4buff[1] = nodes[3]->Y();
            x4buff[2] = nodes[3]->Z();
            volume  += ComputeDetTex(x1buff, x2buff, x3buff, x4buff);

            x1buff[0] = nodes[2]->X();
            x1buff[1] = nodes[2]->Y();
            x1buff[2] = nodes[2]->Z();

            x2buff[0] = nodes[4]->X();
            x2buff[1] = nodes[4]->Y();
            x2buff[2] = nodes[4]->Z();

            x3buff[0] = nodes[5]->X();
            x3buff[1] = nodes[5]->Y();
            x3buff[2] = nodes[5]->Z();

            x4buff[0] = nodes[3]->X();
            x4buff[1] = nodes[3]->Y();
            x4buff[2] = nodes[3]->Z();
            volume  += ComputeDetTex(x1buff, x2buff, x3buff, x4buff);
                                                  // kg44 reactivated ---------Here the direction of flow needs to be taken into account, we need rep length in x,y,z direction
            representative_length = pow(volume,1./3.);
            break;
         default:
            std::cerr << "CElem::ComputeVolume MshElemType not handled" << std::endl;
      }
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   01/2006 YD Implementation
   **************************************************************************/
   void CElem::FaceNormal(int index0, int index1, double* face)
   {
      int i;
      double xx[3];
      double yy[3];
      double zz[3];
      if (GetElementType() == MshElemType::TRIANGLE || GetElementType() == MshElemType::QUAD)
      {
         //----plane normal----------------------------
         // tranform_tensor = new Matrix(3,3);
         // face"_vec
         face[0] = nodes[index1]->X()-nodes[index0]->X();
         face[1] = nodes[index1]->Y()-nodes[index0]->Y();
         face[2] = nodes[index1]->Z()-nodes[index0]->Z();
         // x"_vec
         xx[0] = nodes[1]->X()-nodes[0]->X();
         xx[1] = nodes[1]->Y()-nodes[0]->Y();
         xx[2] = nodes[1]->Z()-nodes[0]->Z();
         //NormalizeVector(xx,3);
         // a vector on the plane
         yy[0] = nodes[2]->X()-nodes[1]->X();
         yy[1] = nodes[2]->Y()-nodes[1]->Y();
         yy[2] = nodes[2]->Z()-nodes[1]->Z();
         // z"_vec
         CrossProduction(xx,yy,zz);
         NormalizeVector(zz,3);
         // y"_vec
         CrossProduction(face,zz,yy);
         NormalizeVector(yy,3);
         for(i=0; i<3; i++)
         {
            face[i] = yy[i];
         }
      }
   }

   /**************************************************************************
   MSHLib-Method:
   06/2006 OK Implementation
   **************************************************************************/
   void CElem::SetNormalVector()
   {
      double v1[3],v2[3];
      if(!normal_vector)
         normal_vector = new double[3];           //WW

      if (this->GetElementType() == MshElemType::TRIANGLE)
      {
         v1[0] = nodes[1]->X() - nodes[0]->X();
         v1[1] = nodes[1]->Y() - nodes[0]->Y();
         v1[2] = nodes[1]->Z() - nodes[0]->Z();
         v2[0] = nodes[2]->X() - nodes[0]->X();
         v2[1] = nodes[2]->Y() - nodes[0]->Y();
         v2[2] = nodes[2]->Z() - nodes[0]->Z();
         CrossProduction(v1,v2,normal_vector);
         NormalizeVector(normal_vector,3);
      }
   }

   // KR 2010/11/16
   void CElem::setElementProperties(MshElemType::type t, bool isFace)
   {
      switch(t)
      {
         case MshElemType::LINE:
            nnodes = 2;
            nnodesHQ = 3;
            ele_dim = 1;
            geo_type = MshElemType::LINE;
            nfaces = 2;
            nedges = 0;
            break;
         case MshElemType::QUAD:
            nnodes = 4;
            nnodesHQ = (isFace) ? 8 : 9;          // if a QUAD is the face of a hex it has 8 nodes, otherwise it has 9
            ele_dim = 2;
            geo_type = MshElemType::QUAD;
            nfaces = 4;
            nedges = 4;
            break;
         case MshElemType::HEXAHEDRON:
            nnodes = 8;
            nnodesHQ = 20;
            ele_dim = 3;
            nfaces = 6;
            nedges = 12;
            geo_type = MshElemType::HEXAHEDRON;
            break;
         case MshElemType::TRIANGLE:
            nnodes = 3;
            nnodesHQ = 6;
            ele_dim = 2;
            geo_type = MshElemType::TRIANGLE;
            nfaces = 3;
            nedges = 3;
            break;
         case MshElemType::TETRAHEDRON:
            nnodes = 4;
            nnodesHQ = 10;
            ele_dim = 3;
            geo_type = MshElemType::TETRAHEDRON;
            nfaces = 4;
            nedges = 6;
            break;
         case MshElemType::PRISM:
            nnodes = 6;
            nnodesHQ = 15;
            ele_dim = 3;
            geo_type = MshElemType::PRISM;
            nfaces = 5;
            nedges = 9;
            break;
         default:
            std::cerr << "CElem::setElementProperties MshElemType not handled" << std::endl;
      }
      this->nodes_index.resize(nnodes);
   }


//  WW. 03.2009
void CElem::WriteGmsh(std::ostream& os,  const int sdom_idx) const
{

   int et=1;
   int ntags=3;
   std::string deli = " ";

   int nn = GetNodesNumber(false);
   if(quadratic)
   {
	   nn = this->GetNodesNumber_H();
      switch(geo_type)
      {
         case MshElemType::LINE:
            et = 8;
            break;    //Line
		 case MshElemType::QUAD:
            et = 10;
            break;    //Quad
		 case MshElemType::HEXAHEDRON:
            et = 12;
            break;    //Hex
		 case MshElemType::TRIANGLE:
            et = 9;
            break;    //Tri
		 case MshElemType::TETRAHEDRON:
            et = 11;
            break;    //Tet
		 case MshElemType::PRISM:
            et = 18;
            break;    //Pris
		 default:
            break;    
      }
   }
   else
   {
      switch(geo_type)
      {
	  case MshElemType::LINE:
            et = 1;
            break;    //Line
	  case MshElemType::QUAD:
            et = 3;
            break;    //Quad
	  case MshElemType::HEXAHEDRON:
            et = 5;
            break;    //Hex
	  case MshElemType::TRIANGLE:
            et = 2;
            break;    //Tri
	  case MshElemType::TETRAHEDRON:
            et = 4;
            break;    //Tet
	  case MshElemType::PRISM:
            et = 6;
            break;    //Pris
      default:
            break;    //Pris
      }
   }
   os<<index+1<<deli<<et<<deli<<ntags<<deli<<GetPatchIndex()+1<<deli<<GetPatchIndex()+1<<deli<<sdom_idx<<deli;
   for(int i=0; i<nn; i++)
      os<<nodes[i]->index + 1<<deli;
   os<<std::endl;
}


}                                                 // namespace Mesh_Group
