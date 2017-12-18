/**************************************************************************
 MSHLib - Object:
 Task:
 Programing:
 08/2005 WW/OK Encapsulation from rf_ele_msh
 last modified
 **************************************************************************/
#include "Configure.h"

#include <cmath>
#include <cfloat>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>                                //WW
#ifdef NON_PROCESS
#include <algorithm>    // std::replace
#endif

#ifndef NON_GEO                                   //WW
#include "gs_project.h"
#endif                                            //#ifndef NON_GEO

// GEOLib
//#include "geo_pnt.h"
//#include "geo_ply.h"

// MathLib
#include "Vector3.h"
#include "MathTools.h"
//WW #include "msh_lib.h"
#include "rf_mmp_new.h"

// MSHLib
#include "msh_mesh.h"
#ifdef BENCHMARKING
#include "benchtimer.h"
#endif
#include "rf_random_walk.h"
// For surface integration. WW. 29.02.2009

#include "mathlib.h"
// FEM
#include "files0.h"
#include "fem_ele.h"

using FiniteElement::CElement;

// PCSLib
extern std::string GetLineFromFile1(std::ifstream*);

#define noMSH_CHECK

size_t max_dim = 0;                               //OK411

//========================================================================
namespace Mesh_Group
{
    using namespace std;

   /**************************************************************************
    FEMLib-Method:
    Task:
    Programing:
    03/2005 OK Implementation
    **************************************************************************/
   CFEMesh::CFEMesh(GEOLIB::GEOObjects* geo_obj, std::string* geo_name) :
   max_mmp_groups (0), _geo_obj (geo_obj), _geo_name (geo_name),
      _n_msh_layer (0), _cross_section (false),
      _msh_n_lines (0), _msh_n_quads (0), _msh_n_hexs (0),
      _msh_n_tris (0), _msh_n_tets (0), _msh_n_prisms (0),
      _min_edge_length (1e-3), _axisymmetry (false)
   {
      useQuadratic = false;
      coordinate_system = 1;

      max_ele_dim = 0;                            //NW
      pcs_name = "NotSpecified";                  //WW
      PT=NULL;                                    // WW+TK
      fm_pcs=NULL;                                //WW
      // 1.11.2007 WW
#ifdef NEW_EQS
      sparse_graph = NULL;
      sparse_graph_H = NULL;
#endif
      map_counter = 0;                            //21.01.2009 WW
      mapping_check = false;                      //23.01.2009 WW
      has_multi_dim_ele = false;                  //NW
      top_surface_checked = false;                // 07.06.2010.  WW
   }

   // Copy-Constructor for CFEMeshes.
   // Programming: 2010/11/10 KR
   CFEMesh::CFEMesh(const CFEMesh* old_mesh)
   {
      if (old_mesh)
      {
         std::cout << "Copying mesh object ... ";

         // mesh nodes
         size_t nNodes = old_mesh->nod_vector.size();
         for (size_t i=0; i<nNodes; i++)
         {
            Mesh_Group::CNode* node = new Mesh_Group::CNode(i);
            double coords[3] = { old_mesh->nod_vector[i]->X(), old_mesh->nod_vector[i]->Y(), old_mesh->nod_vector[i]->Z() };
            node->SetCoordinates(coords);
            this->nod_vector.push_back(node);
         }

         //  mesh elements
         size_t nElems = old_mesh->ele_vector.size();
         for (size_t i=0; i<nElems; i++)
         {
            Mesh_Group::CElem* elem = new Mesh_Group::CElem();
            elem->SetElementType(old_mesh->ele_vector[i]->GetElementType());
            elem->SetPatchIndex(old_mesh->ele_vector[i]->GetPatchIndex());

            size_t nElemNodes = old_mesh->ele_vector[i]->nodes_index.Size();
            elem->SetNodesNumber(nElemNodes);
            elem->nodes_index.resize(nElemNodes);
            for (size_t j=0; j<nElemNodes; j++)
            {
               elem->SetNodeIndex(j, old_mesh->ele_vector[i]->GetNodeIndex(j));
            }

            this->ele_vector.push_back(elem);
         }

         pcs_name = "NotSpecified";
         this->setNumberOfMeshLayers(old_mesh->getNumberOfMeshLayers());

         std::cout << "done." << std::endl;
      }
      else
      {
         std::cout << "Error in CFEMEsh::CFEMesh() - Given parameter is not a valid pointer... " << std::endl;
      }
   }

   /**************************************************************************
    FEMLib-Method:
    Task:
    Programing:
    03/2005 OK Implementation
    07/2005 WW Changes due to the geometry objects
    01/2006 YD Changes for face normal
    03/2010 changed long to size_t to avoid casts
    **************************************************************************/
   CFEMesh::~CFEMesh(void)
   {
      // delete nodes
      for (size_t i = 0; i < nod_vector.size(); i++)
         delete nod_vector[i];
      nod_vector.clear();
      // Edges
      for (size_t i = 0; i < edge_vector.size(); i++)
         delete edge_vector[i];
      edge_vector.clear();
      // Surface faces
      for (size_t i = 0; i < face_vector.size(); i++)
         delete face_vector[i];
      face_vector.clear();
      // Element
      for (size_t i = 0; i < ele_vector.size(); i++)
         delete ele_vector[i];
      ele_vector.clear();
      // normal  YD
      for (size_t i = 0; i < face_normal.size(); i++)
         delete face_normal[i];
      face_normal.clear();

#ifndef NON_GEO                             //WW
      delete PT;                                  // PCH
#endif
      // 1.11.2007 WW
#ifdef NEW_EQS
      if(sparse_graph) delete sparse_graph;
      if(sparse_graph_H) delete sparse_graph_H;
      sparse_graph = NULL;
      sparse_graph_H = NULL;
#endif

   }

   void CFEMesh::setElementType (MshElemType::type type)
   {
      _ele_type = type;
   }

   void CFEMesh::setNumberOfMeshLayers(size_t n_msh_layer)
   {
      _n_msh_layer = n_msh_layer;
   }

   size_t CFEMesh::getNumberOfMeshLayers () const
   {
      return _n_msh_layer;
   }

   bool CFEMesh::hasCrossSection () const
   {
      return _cross_section;
   }

   size_t CFEMesh::getNumberOfLines () const
   {
      return _msh_n_lines;
   }

   size_t CFEMesh::getNumberOfQuads () const
   {
      return _msh_n_quads;
   }

   size_t CFEMesh::getNumberOfHexs () const
   {
      return _msh_n_hexs;
   }

   size_t CFEMesh::getNumberOfTris () const
   {
      return _msh_n_tris;
   }

   size_t CFEMesh::getNumberOfTets () const
   {
      return _msh_n_tets;
   }

   size_t CFEMesh::getNumberOfPrisms () const
   {
      return _msh_n_prisms;
   }

   double CFEMesh::getMinEdgeLength () const
   {
      return _min_edge_length;
   }
   /**************************************************************************
    FEMLib-Method:
    Task:
    Programing:
    03/2005 OK Implementation
    07/2005 WW Changes due to the geometry objects
    08/2005 WW/MB Keyword CrossSection
    09/2005 WW 2D-3D flag
    12/2005 OK MAT_TYPE
    **************************************************************************/
   std::ios::pos_type CFEMesh::Read(std::ifstream *fem_file)
   {
      std::string line_string;
      bool new_keyword = false;
      std::ios::pos_type position;
      double x, y, z;
      this->max_ele_dim = 0;                      //NW

      

      // Keyword loop
      while (!new_keyword)
      {
         position = fem_file->tellg();
         getline(*fem_file, line_string);

         //If Gmsh. 05.2013. WW
		 if(line_string.find("$MeshFormat") != string::npos) //
		 {
            getline(*fem_file, line_string);
            getline(*fem_file, line_string);
            getline(*fem_file, line_string);
            size_t no_nodes, ibuff;
            *fem_file >> no_nodes >> std::ws;
            std::string s;
            for (size_t i = 0; i < no_nodes; i++)
            {
               *fem_file >> ibuff >> x >> y >> z;
               CNode* newNode (new CNode(ibuff, x, y, z));
               nod_vector.push_back(newNode);
            }
                  
            getline(*fem_file, line_string);
            size_t no_elements;
            *fem_file >> no_elements >> std::ws;
            for (size_t i = 0; i < no_elements; i++)
            {
               CElem* newElem (new CElem(i));
               newElem->Read(*fem_file, 2);
               setElementType (newElem->geo_type);
               if (newElem->GetPatchIndex() > max_mmp_groups)
                  max_mmp_groups = newElem->GetPatchIndex();
               if (newElem->GetDimension() > this->max_ele_dim)
                  this->max_ele_dim = newElem->GetDimension();
               ele_vector.push_back(newElem);
            }

            break;
		 }



         if (fem_file->fail())
            break;
         if (line_string.find("#") != std::string::npos)
         {
            new_keyword = true;
            break;
         }
                                                  // subkeyword found
         if (line_string.find("$PCS_TYPE") != std::string::npos)
         {
            *fem_file >> pcs_name >> std::ws;     //WW
            continue;
         }
                                                  // subkeyword found
         if (line_string.find("$GEO_NAME") != std::string::npos)
         {
            *fem_file >> geo_name >> std::ws;     //WW
            continue;
         }
                                                  //OK9_4310
         if (line_string.find("$GEO_TYPE") != std::string::npos)
         {
            *fem_file >> geo_type_name >> geo_name >> std::ws;
            continue;
         }
                                                  // subkeyword found
         if (line_string.find("$AXISYMMETRY") != std::string::npos)
         {
            _axisymmetry = true;
            continue;
         }
                                                  // subkeyword found
         if (line_string.find("$CROSS_SECTION") != std::string::npos)
         {
            _cross_section = true;
            continue;
         }
                                                  // subkeyword found
         if (line_string.find("$NODES") != std::string::npos)
         {
            size_t no_nodes, ibuff;
            *fem_file >> no_nodes >> std::ws;
            std::string s;
            for (size_t i = 0; i < no_nodes; i++)
            {
               *fem_file >> ibuff >> x >> y >> z;
               CNode* newNode (new CNode(ibuff, x, y, z));
               nod_vector.push_back(newNode);
               position = fem_file->tellg();
               *fem_file >> s;
               if (s.find("$AREA") != std::string::npos)
               {
                  *fem_file >> newNode->patch_area;
               } else
               fem_file->seekg(position, std::ios::beg);
               *fem_file >> std::ws;
            }
            continue;
         }
                                                  // subkeyword found
         if (line_string.find("$ELEMENTS") != std::string::npos)
         {
            size_t no_elements;
            *fem_file >> no_elements >> std::ws;
            for (size_t i = 0; i < no_elements; i++)
            {
               CElem* newElem (new CElem(i));
               newElem->Read(*fem_file);
               setElementType (newElem->geo_type);//CC02/2006
               if (newElem->GetPatchIndex() > max_mmp_groups)
                  max_mmp_groups = newElem->GetPatchIndex();
                                                  //NW
               if (newElem->GetDimension() > this->max_ele_dim)
                  this->max_ele_dim = newElem->GetDimension();
               ele_vector.push_back(newElem);
            }
            continue;
         }
                                                  // subkeyword found
         if (line_string.find("$LAYER") != std::string::npos)
         {
            *fem_file >> _n_msh_layer >> std::ws;
            continue;
         }
      }
      return position;
   }

#ifdef NON_PROCESS
   // Read Gmsh. 05.2013. WW 
   void CFEMesh::ReadGmsh(std::ifstream *fem_file)
   {
      std::string line_string;
      double x, y, z;
      max_ele_dim = 0;                     

            getline(*fem_file, line_string);
            getline(*fem_file, line_string);
            getline(*fem_file, line_string);
            size_t no_nodes, ibuff;
            *fem_file >> no_nodes >> std::ws;
            std::string s;

		    size_t max_id = 0; 
            for (size_t i = 0; i < no_nodes; i++)
            {
               *fem_file >> ibuff >> x >> y >> z >> ws;
               CNode* newNode (new CNode(i, x, y, z));
              
			   newNode->SetEquationIndex(ibuff-1); 
               nod_vector.push_back(newNode);

			   if(ibuff-1 > max_id)
				   max_id = ibuff-1;
            }

			vector<size_t> id_map(max_id+1);
            for (size_t i = 0; i < no_nodes; i++)
            {
				id_map[nod_vector[i]->GetEquationIndex()] = nod_vector[i]->GetIndex();
			}

                  
            getline(*fem_file, line_string);
            getline(*fem_file, line_string);
            size_t no_elements;
            *fem_file >> no_elements >> std::ws;
            for (size_t i = 0; i < no_elements; i++)
            {
               CElem* newElem (new CElem(i));
               newElem->Read(*fem_file, 7);

			   for(int k=0; k<newElem->GetNodesNumber(false); k++)
			   {
				   size_t o_node_id = newElem->GetNodeIndex(k);
                   size_t new_node_id = id_map[o_node_id];
				   newElem->SetNodeIndex(k, new_node_id);
			   }
               setElementType (newElem->geo_type);
               if (newElem->GetPatchIndex() > max_mmp_groups)
                  max_mmp_groups = newElem->GetPatchIndex();
               if (newElem->GetDimension() > this->max_ele_dim)
                  this->max_ele_dim = newElem->GetDimension();
               ele_vector.push_back(newElem);
            }
			id_map.clear();
   }
#endif

#ifdef USE_TOKENBUF
   int
      CFEMesh::Read(TokenBuf *tokenbuf)
   {
      std::string sub_line;
      std::string line_string;
      char line_buf[LINE_MAX];
      char str1[LINE_MAX], str2[LINE_MAX];
      bool new_keyword = false;
      std::string hash("#");
      std::ios::pos_type position;
      std::string sub_string,sub_string1;
      long i, ibuff;
      long no_elements;
      long no_nodes;
      double x,y,z;
      CNode* newNode = NULL;
      CElem* newElem = NULL;
#ifdef BENCHMARKING
      BenchTimer read_timer;

      read_timer.start();
#endif

      //========================================================================
      // Keyword loop
      while (!new_keyword && !tokenbuf->done())
      {
         tokenbuf->get_non_empty_line(line_buf, LINE_MAX);
         if(tokenbuf->done()) break;
         line_string = std::string(line_buf);
         if(line_string.find(hash)!=std::string::npos)
         {
            new_keyword = true;
            break;
         }
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("$PCS_TYPE")!=std::string::npos)
         {
            tokenbuf->get_non_empty_line(line_buf, LINE_MAX);
            pcs_name = std::string(line_buf);
            continue;
         }
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("$GEO_NAME")!=std::string::npos)
         {
            tokenbuf->get_non_empty_line(line_buf, LINE_MAX);
            geo_name = std::string(line_buf);
            continue;
         }
         //....................................................................
                                                  //OK9_4310
         if(line_string.find("$GEO_TYPE")!=std::string::npos)
         {
            tokenbuf->get_non_empty_line(line_buf, LINE_MAX);
            sscanf(line_buf, "%s %s", str1, str2);
            geo_type_name = std::string(str1);
            geo_name = std::string(str2);
            continue;
         }
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("$AXISYMMETRY")!=std::string::npos)
         {
            _axisymmetry=true;
            continue;
         }
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("$CROSS_SECTION")!=string::npos)
         {
            cross_section=true;
            continue;
         }
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("$NODES")!=std::string::npos)
         {
            tokenbuf->get_non_empty_line(line_buf, LINE_MAX);
            sscanf(line_buf, "%ld", &no_nodes);
            for(i=0;i<no_nodes;i++)
            {
               tokenbuf->get_non_empty_line(line_buf, LINE_MAX);
               sscanf(line_buf, "%ld %lf %lf %lf", &ibuff, &x, &y, &z);
               newNode = new CNode(ibuff,x,y,z);
               nod_vector.push_back(newNode);
            }
            continue;
         }
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("$ELEMENTS")!=std::string::npos)
         {
            tokenbuf->get_non_empty_line(line_buf, LINE_MAX);
            sscanf(line_buf, "%ld", &no_elements);
            for(i=0;i<no_elements;i++)
            {
               newElem = new CElem(i);
               newElem->Read(tokenbuf);
               setElementType (newElem->geo_type);//CC02/2006
               if(newElem->GetPatchIndex()>max_mmp_groups)
                  max_mmp_groups = newElem->GetPatchIndex();
               ele_vector.push_back(newElem);
            }
            continue;
         }
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("$LAYER")!=std::string::npos)
         {
            tokenbuf->get_non_empty_line(line_buf, LINE_MAX);
            sscanf(line_buf, "%ld", &_n_msh_layer);
            continue;
         }
         //....................................................................
      }
      //========================================================================

#ifdef BENCHMARKING
      read_timer.stop();
      std::cout << "Reading mesh took " << read_timer.time_s() << " s." << std::endl;
#endif
      return 0;
   }
#endif

   /**************************************************************************
    FEMLib-Method: ConnectedElements2Node
    Task:
    Programing:
    04/2007 WW Cut from Construct grid
    **************************************************************************/
   void CFEMesh::ConnectedElements2Node(bool quadratic)
   {
      // Set neighbors of node
      for (size_t e = 0; e < nod_vector.size(); e++)
         nod_vector[e]->connected_elements.clear();

      bool done = false;

      size_t ele_vector_size (ele_vector.size());
      for (size_t e = 0; e < ele_vector_size; e++)
      {
         CElem* celem_e = ele_vector[e];
         if (! celem_e->GetMark())
            continue;
         size_t n_nodes(   static_cast<size_t> (celem_e->GetNodesNumber(quadratic)) );
         for (size_t i = 0; i < n_nodes; i++)
         {
            done = false;
            Mesh_Group::CNode* cnode_i (nod_vector[celem_e->GetNodeIndex(i)]);
            size_t n_connected_elements (cnode_i->connected_elements.size());
            for (size_t j = 0; j < n_connected_elements; j++)
            {
               if (e == static_cast<size_t> (cnode_i->connected_elements[j]))
               {
                  done = true;
                  break;
               }
            }
            if (!done)
               cnode_i->connected_elements.push_back(e);
         }
      }

      //	CElem* thisElem0 = NULL;
      //	for (size_t e = 0; e < (long) ele_vector.size(); e++) {
      //		thisElem0 = ele_vector[e];
      //		if (!thisElem0->GetMark())
      //			continue;
      //		size_t n_nodes(static_cast<size_t>(thisElem0->GetNodesNumber(quadratic)));
      //		for (size_t i = 0; i < n_nodes; i++) {
      //			done = false;
      //			ni = thisElem0->GetNodeIndex(i);
      //			for (j = 0; j < (int) nod_vector[ni]->connected_elements.size(); j++) {
      //				if (e == nod_vector[ni]->connected_elements[j]) {
      //					done = true;
      //					break;
      //				}
      //			}
      //			if (!done)
      //				nod_vector[ni]->connected_elements.push_back(e);
      //		}
      //	}
   }

   /**************************************************************************
    FEMLib-Method: Construct grid
    Task: Establish topology of a grid
    Programing:
    05/2005 WW Implementation
    02/2006 YD Add 1D line neighbor element set
    01/2010 NW Changed to determine the coordinate system by domain dimensions
    **************************************************************************/
   void CFEMesh::ConstructGrid()
   {
      int counter;
      int i, j, k, ii, jj, m0, m, n0, n;
      int nnodes0, nedges0, nedges;
      long e, ei, ee, e_size, e_size_l;
      bool done;
      double x_sum, y_sum, z_sum;

      int edgeIndex_loc0[2];
      int edgeIndex_loc[2];
      int faceIndex_loc0[10];
      int faceIndex_loc[10];
      vec<CNode*> e_nodes0(20);
      //	vec<long> node_index_glb(20);
      //	vec<long> node_index_glb0(20);
      vec<int> Edge_Orientation(15);
      vec<CEdge*> Edges(15);
      vec<CEdge*> Edges0(15);
      vec<CElem*> Neighbors(15);
      vec<CElem*> Neighbors0(15);

      vec<CNode*> e_edgeNodes0(3);
      vec<CNode*> e_edgeNodes(3);
      CElem* thisElem0 = NULL;
      CElem* thisElem = NULL;

      //Elem->nodes not initialized

      e_size = (long) ele_vector.size();
      NodesNumber_Linear = nod_vector.size();

      Edge_Orientation = 1;
      //----------------------------------------------------------------------
      // Set neighbors of node
      ConnectedElements2Node();
      //----------------------------------------------------------------------

      //----------------------------------------------------------------------
      // Compute neighbors and edges
      for (e = 0; e < e_size; e++)
      {
         thisElem0 = ele_vector[e];
         nnodes0 = thisElem0->nnodes;             // Number of nodes for linear element
         //		thisElem0->GetNodeIndeces(node_index_glb0);
         const vec<long>& node_index_glb0 (thisElem0->GetNodeIndeces());
         thisElem0->GetNeighbors(Neighbors0);
         for (i = 0; i < nnodes0; i++)            // Nodes
            e_nodes0[i] = nod_vector[node_index_glb0[i]];
         m0 = thisElem0->GetFacesNumber();
         // neighbors
         for (i = 0; i < m0; i++)                 // Faces
         {
            if (Neighbors0[i])
               continue;
            n0 = thisElem0->GetElementFaceNodes(i, faceIndex_loc0);
            done = false;
            for (k = 0; k < n0; k++)
            {
               e_size_l
                  = (long) e_nodes0[faceIndex_loc0[k]]->connected_elements.size();
               for (ei = 0; ei < e_size_l; ei++)
               {
                  ee = e_nodes0[faceIndex_loc0[k]]->connected_elements[ei];
                  if (ee == e)
                     continue;
                  thisElem = ele_vector[ee];
                  const vec<long>& node_index_glb (thisElem->GetNodeIndeces());
                  thisElem->GetNeighbors(Neighbors);
                  m = thisElem->GetFacesNumber();

                  for (ii = 0; ii < m; ii++)      // Faces
                  {
                     n = thisElem->GetElementFaceNodes(ii, faceIndex_loc);
                     if (n0 != n)
                        continue;
                     counter = 0;
                     for (j = 0; j < n0; j++)
                     {
                        for (jj = 0; jj < n; jj++)
                        {
                           if (node_index_glb0[faceIndex_loc0[j]]
                              == node_index_glb[faceIndex_loc[jj]])
                           {
                              counter++;
                              break;
                           }
                        }
                     }
                     if (counter == n)
                     {
                        Neighbors0[i] = thisElem;
                        Neighbors[ii] = thisElem0;
                        thisElem->SetNeighbor(ii, thisElem0);
                        done = true;
                        break;
                     }
                  }
                  if (done)
                     break;
               }
               if (done)
                  break;
            }
         }
         thisElem0->SetNeighbors(Neighbors0);
         //------------neighbor of 1D line
         if (thisElem0->geo_type == 1)            //YD
         {
            ii = 0;
            for (i = 0; i < m0; i++)
            {
               n0 = thisElem0->GetElementFaceNodes(i, faceIndex_loc0);
               for (k = 0; k < n0; k++)
               {
                  e_size_l
                     = (long) e_nodes0[faceIndex_loc0[k]]->connected_elements.size();
                  for (ei = 0; ei < e_size_l; ei++)
                  {
                     ee
                        = e_nodes0[faceIndex_loc0[k]]->connected_elements[ei];
                     thisElem = ele_vector[ee];
                     if (e_size_l == 2 && thisElem->GetIndex()
                        != thisElem0->GetIndex())
                     {
                        Neighbors0[i] = thisElem;
                        Neighbors[ii] = thisElem;
                        // thisElem->SetNeighbor(ii, thisElem0);   //?? Todo YD
                        ii++;
                     }
                  }
               }
            }
            thisElem0->SetNeighbors(Neighbors0);
         }
         // --------------------------------
         // Edges
         nedges0 = thisElem0->GetEdgesNumber();
         thisElem0->GetEdges(Edges0);
         for (i = 0; i < nedges0; i++)
         {
            thisElem0->GetLocalIndicesOfEdgeNodes(i, edgeIndex_loc0);
            // Check neighbors
            done = false;
            for (k = 0; k < 2; k++)
            {
               e_size_l
                  = (long) e_nodes0[edgeIndex_loc0[k]]->connected_elements.size();
               for (ei = 0; ei < e_size_l; ei++)
               {
                  ee = e_nodes0[edgeIndex_loc0[k]]->connected_elements[ei];
                  if (ee == e)
                     continue;
                  thisElem = ele_vector[ee];
                  const vec<long>& node_index_glb (thisElem->GetNodeIndeces());
                  nedges = thisElem->GetEdgesNumber();
                  thisElem->GetEdges(Edges);
                  // Edges of neighbors
                  for (ii = 0; ii < nedges; ii++)
                  {
                     thisElem->GetLocalIndicesOfEdgeNodes(ii, edgeIndex_loc);
                     if ((node_index_glb0[edgeIndex_loc0[0]]
                        == node_index_glb[edgeIndex_loc[0]]
                        && node_index_glb0[edgeIndex_loc0[1]]
                        == node_index_glb[edgeIndex_loc[1]])
                        || (node_index_glb0[edgeIndex_loc0[0]]
                        == node_index_glb[edgeIndex_loc[1]]
                        && node_index_glb0[edgeIndex_loc0[1]]
                        == node_index_glb[edgeIndex_loc[0]]))
                     {
                        if (Edges[ii])
                        {
                           Edges0[i] = Edges[ii];
                           Edges[ii]->GetNodes(e_edgeNodes);
                           if ((size_t)node_index_glb0[edgeIndex_loc0[0]]
                              == e_edgeNodes[1]->GetIndex()
                              && (size_t)node_index_glb0[edgeIndex_loc0[1]]
                              == e_edgeNodes[0]->GetIndex())
                              Edge_Orientation[i] = -1;
                           done = true;
                           break;
                        }
                     }
                  }                               //  for(ii=0; ii<nedges; ii++)
                  if (done)
                     break;
               }                                  // for(ei=0; ei<e_size_l; ei++)
               if (done)
                  break;
            }                                     //for(k=0;k<2;k++)
            if (!done)                            // new edges and new node
            {
               Edges0[i] = new CEdge((long) edge_vector.size());
               Edges0[i]->SetOrder(false);
               e_edgeNodes0[0] = e_nodes0[edgeIndex_loc0[0]];
               e_edgeNodes0[1] = e_nodes0[edgeIndex_loc0[1]];
               e_edgeNodes0[2] = NULL;
               Edges0[i]->SetNodes(e_edgeNodes0);
               edge_vector.push_back(Edges0[i]);
            }                                     // new edges
         }                                        //  for(i=0; i<nedges0; i++)
         //
         // Set edges and nodes
         thisElem0->SetOrder(false);
         thisElem0->SetEdgesOrientation(Edge_Orientation);
         thisElem0->SetEdges(Edges0);
         // Resize is true
         thisElem0->SetNodes(e_nodes0, true);
      }                                           // Over elements
      // Set faces on surfaces and others
      _msh_n_lines = 0;                           // Should be members of mesh
      _msh_n_quads = 0;
      _msh_n_hexs = 0;
      _msh_n_tris = 0;
      _msh_n_tets = 0;
      _msh_n_prisms = 0;
      for (e = 0; e < e_size; e++)
      {
         thisElem0 = ele_vector[e];
         switch (thisElem0->GetElementType())
         {
            case MshElemType::LINE:
               _msh_n_lines++;
               break;
            case MshElemType::QUAD:
               _msh_n_quads++;
               break;
            case MshElemType::HEXAHEDRON:
               _msh_n_hexs++;
               break;
            case MshElemType::TRIANGLE:
               _msh_n_tris++;
               break;
            case MshElemType::TETRAHEDRON:
               _msh_n_tets++;
               break;
            case MshElemType::PRISM:
               _msh_n_prisms++;
               break;
            default:
               std::cerr << "CFEMesh::ConstructGrid MshElemType not handled" << std::endl;
         }
         // Compute volume meanwhile
         thisElem0->ComputeVolume();

         if (thisElem0->GetElementType() == MshElemType::LINE)
            continue;                             // line element
         //		thisElem0->GetNodeIndeces(node_index_glb0);
         //		const vec<long>& node_index_glb0 (thisElem0->GetNodeIndeces()); // compiler said: unused variable // TF
         thisElem0->GetNeighbors(Neighbors0);
         m0 = thisElem0->GetFacesNumber();

         // Check face on surface
         for (i = 0; i < m0; i++)                 // Faces
         {
            if (Neighbors0[i])
               continue;
            CElem* newFace = new CElem((long) face_vector.size(), thisElem0, i);
            //          thisElem0->boundary_type='B';
            thisElem0->no_faces_on_surface++;
            face_vector.push_back(newFace);
            Neighbors0[i] = newFace;
         }
         thisElem0->SetNeighbors(Neighbors0);

      }
      NodesNumber_Quadratic = (long) nod_vector.size();
      if ((_msh_n_hexs + _msh_n_tets + _msh_n_prisms) > 0)
         max_ele_dim = 3;
      else if ((_msh_n_quads + _msh_n_tris) > 0)
         max_ele_dim = 2;
      else
         max_ele_dim = 1;

      // check if this mesh includes multi-dimensional elements
      if (max_ele_dim == 2 && _msh_n_lines > 0)   //NW
      {
         this->has_multi_dim_ele = true;
      } else if (max_ele_dim == 3 && (_msh_n_quads + _msh_n_tris + _msh_n_lines)
         > 0)
      {
         this->has_multi_dim_ele = true;
      }

      //----------------------------------------------------------------------
      // Node information
      // 1. Default node index <---> eqs index relationship
      // 2. Coordiate system flag
      x_sum = 0.0;
      y_sum = 0.0;
      z_sum = 0.0;
      Eqs2Global_NodeIndex.clear();
      double xyz_max[3] =                         //NW
      {
         -DBL_MAX, -DBL_MAX, -DBL_MAX
      };
      double xyz_min[3] =                         //NW
      {
         DBL_MAX, DBL_MAX, DBL_MAX
      };

      for (e = 0; e < (long) nod_vector.size(); e++)
      {
         nod_vector[e]->SetEquationIndex(e);
         Eqs2Global_NodeIndex.push_back(nod_vector[e]->GetIndex());
         x_sum += fabs(nod_vector[e]->X());
         y_sum += fabs(nod_vector[e]->Y());
         z_sum += fabs(nod_vector[e]->Z());
         if (nod_vector[e]->X() > xyz_max[0])
            xyz_max[0] = nod_vector[e]->X();
         if (nod_vector[e]->Y() > xyz_max[1])
            xyz_max[1] = nod_vector[e]->Y();
         if (nod_vector[e]->Z() > xyz_max[2])
            xyz_max[2] = nod_vector[e]->Z();
         if (nod_vector[e]->X() < xyz_min[0])
            xyz_min[0] = nod_vector[e]->X();
         if (nod_vector[e]->Y() < xyz_min[1])
            xyz_min[1] = nod_vector[e]->Y();
         if (nod_vector[e]->Z() < xyz_min[2])
            xyz_min[2] = nod_vector[e]->Z();
      }
      double xyz_dim[3];                          //NW
      xyz_dim[0] = xyz_max[0] - xyz_min[0];
      xyz_dim[1] = xyz_max[1] - xyz_min[1];
      xyz_dim[2] = xyz_max[2] - xyz_min[2];

      //check dimension of the domain to select appropriate coordinate system
      if (xyz_dim[0] > 0.0 && xyz_dim[1] < MKleinsteZahl && xyz_dim[2]
         < MKleinsteZahl)
         coordinate_system = 10;
      else if (xyz_dim[1] > 0.0 && xyz_dim[0] < MKleinsteZahl && xyz_dim[2]
         < MKleinsteZahl)
         coordinate_system = 11;
      else if (xyz_dim[2] > 0.0 && xyz_dim[0] < MKleinsteZahl && xyz_dim[1]
         < MKleinsteZahl)
         coordinate_system = 12;
      else if (xyz_dim[0] > 0.0 && xyz_dim[1] > 0.0 && xyz_dim[2] < MKleinsteZahl)
         coordinate_system = 21;
      else if (xyz_dim[0] > 0.0 && xyz_dim[2] > 0.0 && xyz_dim[1] < MKleinsteZahl)
         coordinate_system = 22;
      else if (xyz_dim[0] > 0.0 && xyz_dim[1] > 0.0 && xyz_dim[2] > 0.0)
         coordinate_system = 32;

      // 1D in 2D
      if (_msh_n_lines > 0)
      {
         if (xyz_dim[0] > 0.0 && xyz_dim[1] > 0.0 && xyz_dim[2] < MKleinsteZahl)
            coordinate_system = 32;
         if (xyz_dim[0] > 0.0 && xyz_dim[2] > 0.0 && xyz_dim[1] < MKleinsteZahl)
            coordinate_system = 32;
      }

      max_dim = coordinate_system / 10 - 1;
      //----------------------------------------------------------------------
      // Gravity center
      for (e = 0; e < e_size; e++)
      {
         thisElem0 = ele_vector[e];
         nnodes0 = thisElem0->nnodes;
         for (i = 0; i < nnodes0; i++)            // Nodes
         {
            thisElem0->gravity_center[0] += thisElem0->nodes[i]->X();
            thisElem0->gravity_center[1] += thisElem0->nodes[i]->Y();
            thisElem0->gravity_center[2] += thisElem0->nodes[i]->Z();
         }
         thisElem0->gravity_center[0] /= (double) nnodes0;
         thisElem0->gravity_center[1] /= (double) nnodes0;
         thisElem0->gravity_center[2] /= (double) nnodes0;
      }
      //----------------------------------------------------------------------

      //TEST WW
      // For sparse matrix
      ConnectedNodes(false);
      //
      e_nodes0.resize(0);
      //	node_index_glb.resize(0);
      //	node_index_glb0.resize(0);
      Edge_Orientation.resize(0);
      Edges.resize(0);
      Edges0.resize(0);
      Neighbors.resize(0);
      Neighbors0.resize(0);
      e_edgeNodes0.resize(0);
      e_edgeNodes.resize(0);
   }

   /**************************************************************************
    FEMLib-Method: GenerateHighOrderNodes()
    Task:
    Programing:
    07/2007 WW Implementation
    01/2010 NW Case: a mesh with line elements
    **************************************************************************/
   void CFEMesh::GenerateHighOrderNodes()
   {
      int i, j, k, ii;
      int nnodes0, nedges0, nedges;
      long e, ei, ee, e_size, e_size_l;
      int edgeIndex_loc0[2];
      bool done;
      double x0 = 0.0, y0 = 0.0, z0 = 0.0;        //OK411

      // Set neighbors of node. All elements, even in deactivated subdomains, are taken into account here.
      for (e = 0; e < (long) nod_vector.size(); e++)
         nod_vector[e]->connected_elements.clear();
      done = false;
      for (e = 0; e < (long) ele_vector.size(); e++)
      {
         CElem *thisElem0 = ele_vector[e];
         for (i = 0; i < thisElem0->GetNodesNumber(false); i++)
         {
            done = false;
            long ni = thisElem0->GetNodeIndex(i);
            for (j = 0; j < (int) nod_vector[ni]->connected_elements.size(); j++)
            {
               if (e == nod_vector[ni]->connected_elements[j])
               {
                  done = true;
                  break;
               }
            }
            if (!done)
               nod_vector[ni]->connected_elements.push_back(e);
         }
      }
      //
      CNode *aNode = NULL;
      vec<CNode*> e_nodes0(20);
      vec<CNode*> e_nodes(20);
      CElem *thisElem0 = NULL;
      CElem *thisElem = NULL;
      CEdge *thisEdge0 = NULL;
      CEdge *thisEdge = NULL;
      //----------------------------------------------------------------------
      // Loop over elements (except for line elements)
      e_size = (long) ele_vector.size();
      for (e = 0; e < e_size; e++)
      {
         thisElem0 = ele_vector[e];
         if (thisElem0->GetElementType() == MshElemType::LINE)
            continue;                             //NW

         nnodes0 = thisElem0->nnodes;             // Number of nodes for linear element
         //      thisElem0->GetNodeIndeces(node_index_glb0);
         for (i = 0; i < nnodes0; i++)            // Nodes
            e_nodes0[i] = thisElem0->GetNode(i);
         // --------------------------------
         // Edges
         nedges0 = thisElem0->GetEdgesNumber();
         // Check if there is any neighbor that has new middle points
         for (i = 0; i < nedges0; i++)
         {
            thisEdge0 = thisElem0->GetEdge(i);
            thisElem0->GetLocalIndicesOfEdgeNodes(i, edgeIndex_loc0);
            // Check neighbors
            done = false;
            for (k = 0; k < 2; k++)
            {
               e_size_l
                  = (long) e_nodes0[edgeIndex_loc0[k]]->connected_elements.size();
               for (ei = 0; ei < e_size_l; ei++)
               {
                  ee = e_nodes0[edgeIndex_loc0[k]]->connected_elements[ei];
                  if (ee == e)
                     continue;
                  thisElem = ele_vector[ee];
                  nedges = thisElem->GetEdgesNumber();
                  // Edges of neighbors
                  for (ii = 0; ii < nedges; ii++)
                  {
                     thisEdge = thisElem->GetEdge(ii);
                     if (*thisEdge0 == *thisEdge)
                     {
                        aNode = thisEdge->GetNode(2);
                        if (aNode)                // The middle point exist
                        {
                           e_nodes0[nnodes0] = aNode;
                           nnodes0++;
                           done = true;
                           break;
                        }
                     }
                  }                               //  for(ii=0; ii<nedges; ii++)
                  if (done)
                     break;
               }                                  // for(ei=0; ei<e_size_l; ei++)
               if (done)
                  break;
            }                                     //for(k=0;k<2;k++)
            if (!done)
            {
               aNode = new CNode((long) nod_vector.size());
               aNode->SetX(0.5 * (thisEdge0->GetNode(0)->X()
                  + thisEdge0->GetNode(1)->X()));
               aNode->SetY(0.5 * (thisEdge0->GetNode(0)->Y()
                  + thisEdge0->GetNode(1)->Y()));
               aNode->SetZ(0.5 * (thisEdge0->GetNode(0)->Z()
                  + thisEdge0->GetNode(1)->Z()));
               e_nodes0[nnodes0] = aNode;
               thisEdge0->SetNode(2, aNode);
               nnodes0++;
               nod_vector.push_back(aNode);
            }
         }                                        //  for(i=0; i<nedges0; i++)

         // No neighors or no neighbor has new middle point
         //
                                                  // Quadrilateral
         if (thisElem0->GetElementType() == MshElemType::QUAD)
         {
            x0 = y0 = z0 = 0.0;
            aNode = new CNode((long) nod_vector.size());
            e_nodes0[nnodes0] = aNode;
            nnodes0 = thisElem0->nnodes;
            for (i = 0; i < nnodes0; i++)         // Nodes
            {
               x0 += e_nodes0[i]->X();
               y0 += e_nodes0[i]->Y();
               z0 += e_nodes0[i]->Z();
            }
            x0 /= (double) nnodes0;
            y0 /= (double) nnodes0;
            z0 /= (double) nnodes0;
            aNode->SetX(x0);
            aNode->SetY(y0);
            aNode->SetZ(z0);
            nod_vector.push_back(aNode);
         }
         // Set edges and nodes
         thisElem0->SetOrder(true);
         // Resize is true
         thisElem0->SetNodes(e_nodes0, true);
      }                                           // Over elements

      // Setup 1d line elements at the end
      if (_msh_n_lines > 0)
      {
         for (e = 0; e < e_size; e++)
         {
            thisElem0 = ele_vector[e];
            if (thisElem0->GetElementType() != MshElemType::LINE)
               continue;

            nnodes0 = thisElem0->nnodes;
            for (i = 0; i < nnodes0; i++)
               e_nodes0[i] = thisElem0->GetNode(i);

            done = false;

            for (i = 0; i < thisElem0->GetFacesNumber(); i++)
            {
               thisElem = thisElem0->GetNeighbor(i);
               // look for adjacent solid elements
               if (thisElem->GetElementType() == MshElemType::LINE)
                  continue;

               for (j = 0; j < thisElem->nnodes; j++)
                  e_nodes[j] = thisElem->GetNode(j);
               nedges = thisElem->GetEdgesNumber();
               // search a edge connecting to this line element
               for (j = 0; j < nedges; j++)
               {
                  thisEdge = thisElem->GetEdge(j);
                  thisElem->GetLocalIndicesOfEdgeNodes(j, edgeIndex_loc0);
                  // Check neighbors
                  for (k = 0; k < 2; k++)
                  {
                     //OK411 CNode *tmp_nod = e_nodes[edgeIndex_loc0[k]];
                     e_size_l
                        = (long) e_nodes[edgeIndex_loc0[k]]->connected_elements.size();
                     for (ei = 0; ei < e_size_l; ei++)
                     {
                        ee
                           = e_nodes[edgeIndex_loc0[k]]->connected_elements[ei];
                        if (ele_vector[ee] != thisElem0)
                           continue;
                        //the edge is found now
                        aNode = thisEdge->GetNode(2);
                        if (aNode)                // The middle point exist
                        {
                           e_nodes0[nnodes0] = aNode;
                           nnodes0++;
                           done = true;
                           break;
                        }
                        if (done)
                           break;
                     }                            // for(ei=0; ei<e_size_l; ei++)
                     if (done)
                        break;
                  }                               //for(k=0;k<2;k++)
                  if (done)
                     break;
               }                                  //  for(i=0; i<nedges0; i++)
               if (done)
                  break;
            }
            if (!done)
            {
               aNode = new CNode((long) nod_vector.size());
               for (i = 0; i < nnodes0; i++)      // Nodes
               {
                  x0 += e_nodes0[i]->X();
                  y0 += e_nodes0[i]->Y();
                  z0 += e_nodes0[i]->Z();
               }
               x0 /= (double) nnodes0;
               y0 /= (double) nnodes0;
               z0 /= (double) nnodes0;
               aNode->SetX(x0);
               aNode->SetY(y0);
               aNode->SetZ(z0);
               e_nodes0[nnodes0] = aNode;
               nnodes0++;
               nod_vector.push_back(aNode);
            }
            thisElem0->SetOrder(true);
            thisElem0->SetNodes(e_nodes0, true);
         }
      }
      //
      NodesNumber_Quadratic = (long) nod_vector.size();
      for (e = NodesNumber_Linear; (size_t)e < NodesNumber_Quadratic; e++)
      {
         nod_vector[e]->SetEquationIndex(e);
         Eqs2Global_NodeIndex.push_back(nod_vector[e]->GetIndex());
      }
      for (e = 0; e < e_size; e++)
      {
         thisElem0 = ele_vector[e];
         for (i = thisElem0->nnodes; i < thisElem0->nnodesHQ; i++)
         {
            done = false;
            aNode = thisElem0->GetNode(i);
            for (k = 0; k < (int) aNode->connected_elements.size(); k++)
            {
               if (e == aNode->connected_elements[k])
               {
                  done = true;
                  break;
               }
            }
            if (!done)
               aNode->connected_elements.push_back(e);
         }
      }

      // For sparse matrix
      ConnectedNodes(true);
      ConnectedElements2Node(true);
      //
      e_nodes0.resize(0);

      // Test	WW
      /*
       fstream n_out;
       n_out.open("node.txt", ios::out );
       for(e=0; e<NodesNumber_Quadratic; e++)
       nod_vector[e]->Write(n_out);
       n_out.close();
       */
   }

   /**************************************************************************
    FEMLib-Method:
    Task:  Renumbering nodes corresponding to the activate of elements
    Programing:
    09/2005 WW Implementation
    05/2007 WW 1D in 2D
    **************************************************************************/
   void CFEMesh::FillTransformMatrix()
   {
      CElem* elem = NULL;
#ifndef NON_PROCESS                         //05.01.2011. WW
                                                  // PCH
      CRFProcess* m_pcs = PCSGet("FLUID_MOMENTUM");
#endif                                      //#ifndef NON_PROCESS

      //
      if ((_msh_n_hexs + _msh_n_tets + _msh_n_prisms) == ele_vector.size())
         return;
      else if (coordinate_system != 32 && !this->has_multi_dim_ele)
      {
#ifndef NON_PROCESS                      //05.01.2011. WW
         if (m_pcs)
            ;                                     // Need to do FillTransformMatrix	// PCH
         else
            return;
#endif                                   //#ifndef NON_PROCESS
      }
      bool tilted = false;
      if (coordinate_system == 32 || coordinate_system == 21 || coordinate_system
         == 22)
         tilted = true;
      if (!tilted)
         return;
      for (size_t i = 0; i < ele_vector.size(); i++)
      {
         elem = ele_vector[i];
         if (elem->GetMark())                     // Marked for use
         {
            if (coordinate_system == 21 || coordinate_system == 22)
            {
               if (elem->GetElementType() == MshElemType::LINE)
                  elem->FillTransformMatrix();
            }
            else
            {
               if (elem->GetElementType() == MshElemType::LINE || elem->GetElementType() == MshElemType::QUAD
                  || elem->GetElementType() == MshElemType::TRIANGLE)
                  elem->FillTransformMatrix();
            }
         }
      }
   }

   /**************************************************************************
    FEMLib-Method:
    Task:  Renumbering nodes corresponding to the activiate of elements
    Programing:
    05/2005 WW Implementation
    **************************************************************************/
   void CFEMesh::RenumberNodesForGlobalAssembly()
   {
      int i;
      long l, el, el_0;

      CElem* elem = NULL;
      Eqs2Global_NodeIndex.clear();
      for (l = 0; l < (long) nod_vector.size(); l++)
         nod_vector[l]->SetEquationIndex(-1);

      el_0 = 0;
      // Lower order
      for (l = 0; l < (long) ele_vector.size(); l++)
      {
         elem = ele_vector[l];
         if (elem->GetMark())                     // Marked for use
         {
            for (i = 0; i < elem->GetVertexNumber(); i++)
            {
               if (elem->nodes[i]->GetEquationIndex() < 0)
               {
                  elem->nodes[i]->SetEquationIndex(el_0);
                  Eqs2Global_NodeIndex.push_back(elem->nodes[i]->GetIndex());
                  el_0++;
               } else
               continue;
            }
         }
      }
      el = el_0;
      if (!getOrder())
      {
         NodesNumber_Linear = el_0;
         NodesNumber_Quadratic = el;
         return;
      }
      // High order
      for (l = 0; l < (long) ele_vector.size(); l++)
      {
         elem = ele_vector[l];
         if (elem->GetMark())                     // Marked for use
         {

            for (i = elem->GetVertexNumber(); i < elem->GetNodesNumber(true); i++)
            {
               if (elem->nodes[i]->GetEquationIndex() < 0)
               {
                  elem->nodes[i]->SetEquationIndex(el);
                  Eqs2Global_NodeIndex.push_back(elem->nodes[i]->GetIndex());
                  el++;
               } else
               continue;
            }
         }
      }
      NodesNumber_Linear = el_0;
      NodesNumber_Quadratic = el;
   }
   /**************************************************************************
    FEMLib-Method:
    Task:
    Programing:
    03/2005 OK Implementation
    07/2005 WW Write by member methods of geometry objects.
    12/2005 OK MAT_TYPE
    **************************************************************************/
   void CFEMesh::Write(std::fstream*fem_msh_file) const
   {
      long i;
      //--------------------------------------------------------------------
      //KEYWORD
      *fem_msh_file << "#FEM_MSH" << std::endl;
      //--------------------------------------------------------------------
      // PCS
      *fem_msh_file << " $PCS_TYPE" << std::endl;
      *fem_msh_file << "  ";
      *fem_msh_file << pcs_name << std::endl;
      //--------------------------------------------------------------------
      // MAT
      if (geo_name.size() > 0)
      {
         *fem_msh_file << " $GEO_TYPE" << std::endl;
         *fem_msh_file << "  ";
                                                  //OK10_4310
         *fem_msh_file << geo_type_name << " " << geo_name << std::endl;
      }
      //--------------------------------------------------------------------
      // NODES
      *fem_msh_file << " $NODES" << std::endl;
      *fem_msh_file << "  ";
                                                  //WW
      *fem_msh_file << GetNodesNumber(false) << std::endl;
      for (i = 0; i < (long) nod_vector.size(); i++)
         nod_vector[i]->Write(*fem_msh_file);     //WW
      //--------------------------------------------------------------------
      // ELEMENTS
      *fem_msh_file << " $ELEMENTS" << std::endl;
      *fem_msh_file << "  ";
      *fem_msh_file << (long) ele_vector.size() << std::endl;
      for (i = 0; i < (long) ele_vector.size(); i++)
      {
         ele_vector[i]->SetIndex(i);              //20.01.06 WW/TK
         ele_vector[i]->WriteIndex(*fem_msh_file);//WW
      }
      //--------------------------------------------------------------------
      *fem_msh_file << " $LAYER" << std::endl;
      *fem_msh_file << "  ";
      *fem_msh_file << _n_msh_layer << std::endl;
      //--------------------------------------------------------------------
      *fem_msh_file << "#STOP"<<endl;
   }

#ifndef NON_GEO
   /**************************************************************************
    FEMLib-Method:
    Task: Ermittelt den nahliegenden existierenden Knoten
    Programing:
    03/2010 TF implementation based on long CFEMesh::GetNODOnPNT(CGLPoint*m_pnt)
    by OK, WW
    **************************************************************************/
   long CFEMesh::GetNODOnPNT(const GEOLIB::Point* const pnt) const
   {
      double sqr_dist(0.0), distmin(std::numeric_limits<double>::max());
      long number(-1);

      for (size_t i=0; i<(size_t)NodesInUsage(); i++)
      {
         sqr_dist = 0.0;
         sqr_dist += (nod_vector[i]->X() - (*pnt)[0]) * (nod_vector[i]->X()
            - (*pnt)[0]);
         sqr_dist += (nod_vector[i]->Y() - (*pnt)[1]) * (nod_vector[i]->Y()
            - (*pnt)[1]);
         sqr_dist += (nod_vector[i]->Z() - (*pnt)[2]) * (nod_vector[i]->Z()
            - (*pnt)[2]);
         if (sqr_dist < distmin)
         {
            distmin = sqr_dist;
            number = i;
         }
      }
      return number;
   }

   /**************************************************************************
    FEMLib-Method:
    Task: Ermittelt das nahliegende Element
    Programing:
    03/2010 TF implementation based on long CFEMesh::GetNODOnPNT(CGLPoint*m_pnt)
    by MB
    **************************************************************************/
   long CFEMesh::GetNearestELEOnPNT(const GEOLIB::Point* const pnt) const
   {
      long nextele(-1);
      double dist(std::numeric_limits<double>::max()), dist1;
      double* center(NULL);

      for (size_t i = 0; i < ele_vector.size(); i++)
      {
         center = ele_vector[i]->GetGravityCenter();
         dist1 = 0.0;
         for (size_t k(0); k<3; k++) dist1 += (center[k]-(*pnt)[k]) * (center[k]-(*pnt)[k]);
         if (dist1 < dist)
         {
            dist = dist1;
            nextele = i;
         }
      }
      return nextele;
   }

   //WW. (x1-x0).(x2-x0)
   inline double dotProduction(const double *x1, const double *x2,
      const double *x0)
   {
      return (x1[0] - x0[0]) * (x2[0] - x0[0]) + (x1[1] - x0[1])
         * (x2[1] - x0[1]) + (x1[2] - x0[2]) * (x2[2] - x0[2]);
   }

   /**************************************************************************
    FEMLib-Method: GetNodesOnArc
    Task:  To get fe nodes on a Arc defined by
    start point, center point, end point
    Programing:
    01/2005 WW Implementation
    05/2005 WW Transplant to this object from GEOLIB
    **************************************************************************/
   void CFEMesh::GetNodesOnArc(CGLPolyline*m_ply, std::vector<long>&msh_nod_vector)
   {
      long i;

      //Obtain fem node for groupvector
      const long nNodes = NodesInUsage();
      const int SizeCGLPoint = (int) m_ply->point_vector.size();
      double r1, r2, a0, a1, a2;
      const double Tol = 1.0e-5;
      double xa[3], xb[3], xc[3], xn[3];
      msh_nod_vector.clear();
      if (SizeCGLPoint > 3)
      {
         std::cout << "More than 3 points that an arc needed" << std::endl;
         abort();
      }
      else if (SizeCGLPoint == 0)
      {
         std::cout << "No point are given for this arc" << std::endl;
         abort();
      }

      xa[0] = m_ply->point_vector[0]->x;
      xa[1] = m_ply->point_vector[0]->y;
      xa[2] = m_ply->point_vector[0]->z;
      xb[0] = m_ply->point_vector[2]->x;
      xb[1] = m_ply->point_vector[2]->y;
      xb[2] = m_ply->point_vector[2]->z;
      xc[0] = m_ply->point_vector[1]->x;
      xc[1] = m_ply->point_vector[1]->y;
      xc[2] = m_ply->point_vector[1]->z;

      r1 = sqrt(MATHLIB::sqrDist(xa, xc));
      r2 = sqrt(MATHLIB::sqrDist(xb, xc));

      if (fabs(r1 - r2) > m_ply->epsilon)
      {
         std::cout << "Start point and end point do not have identical "
            << "distance to the center of the arc" << std::endl;
         abort();
      }
      a0 = acos(dotProduction(xa, xb, xc) / (r1 * r1));

      // Check nodes by comparing distance
      for (i = 0; i < nNodes; i++)
      {
         xn[0] = nod_vector[i]->X();
         xn[1] = nod_vector[i]->Y();
         xn[2] = nod_vector[i]->Z();
         r2 = sqrt(MATHLIB::sqrDist(xn, xc));
         if (fabs(r2 - r1) < m_ply->epsilon)
         {
            if (a0 < Tol)                         // Closed arc
               msh_nod_vector.push_back(nod_vector[i]->GetIndex());
            else
            {
               a1 = acos(dotProduction(xa, xn, xc) / (r1 * r2));
               a2 = acos(dotProduction(xb, xn, xc) / (r1 * r2));
               if (fabs(a1 + a2 - a0) < Tol)
                  msh_nod_vector.push_back(nod_vector[i]->GetIndex());
            }
         }
      }

      m_ply->GetPointOrderByDistance();
   }

   /**************************************************************************
    FEMLib-Method: GetNodesOnArc
    Task:  To get fe nodes on a Arc defined by start point, center point, end point
    Programing:
    01/2005 WW Implementation
    05/2005 WW Transplant to this object from GEOLIB
    03/2010 TF change to new datastructure GEOLIB::Polyline - some improvements
    **************************************************************************/
   //void CFEMesh::GetNodesOnArc(const GEOLIB::Point* a, const GEOLIB::Point* m,
   //		const GEOLIB::Point* b, std::vector<size_t>& msh_nod_vector) const
   //{
   //	msh_nod_vector.clear();
   //	MATHLIB::Vector v_am(*a, *m);
   //	MATHLIB::Vector v_bm(*b, *m);
   //	double sqr_radius1 = MATHLIB::scpr(v_am.getData(), v_am.getData(), 3);
   //	double sqr_radius2 = MATHLIB::scpr(v_bm.getData(), v_bm.getData(), 3);
   //
   //	// since we describe a circle arc, the radii should match very good
   //	// relative error criterion
   //	if (fabs(sqr_radius1 - sqr_radius2)
   //			> std::numeric_limits<double>::epsilon() * sqr_radius1) {
   //		std::cout << "Start point and end point do not have identical "
   //				<< "distance to the center of the arc" << std::endl;
   //		abort();
   //	}
   //
   //	// compute the angle at centre
   //	double a0 = acos(MATHLIB::scpr(v_am.getData(), v_bm.getData(), 3) / sqr_radius1);
   //
   //	// obtain fem node for groupvector
   //	const size_t nNodes = NodesInUsage();
   //	double tol(1e-5);
   //	// small arc - put all mesh nodes within the annulus into msh_nod_vector
   //	if (a0 < tol) {
   //		// Check nodes by comparing distance
   //		for (size_t i = 0; i < nNodes; i++) {
   //			double xn[3];
   //			xn[0] = nod_vector[i]->X();
   //			xn[1] = nod_vector[i]->Y();
   //			xn[2] = nod_vector[i]->Z();
   //			// compute vector middle point of circle arc and mesh node
   //			for (size_t k = 0; k < 3; k++)
   //				xn[k] -= (*m)[k];
   //
   //			sqr_radius2 = MATHLIB::scpr(xn, xn, 3);
   //			if (fabs(sqr_radius2 - sqr_radius1) < _min_edge_length) {
   //				msh_nod_vector.push_back(nod_vector[i]->GetIndex());
   //			}
   //		}
   //	} else { // arc a0 > tol
   //		for (size_t i = 0; i < nNodes; i++) {
   //			double xn[3];
   //			xn[0] = nod_vector[i]->X();
   //			xn[1] = nod_vector[i]->Y();
   //			xn[2] = nod_vector[i]->Z();
   //			// compute vector middle point of circle arc and mesh node
   //			for (size_t k = 0; k < 3; k++)
   //				xn[k] -= (*m)[k];
   //
   //			sqr_radius2 = MATHLIB::scpr(xn, xn, 3);
   //			if (fabs(sqr_radius2 - sqr_radius1) < _min_edge_length) {
   //				double a1 = acos(MATHLIB::scpr(v_am.getData(), xn, 3) / sqr_radius1);
   //				double a2 = acos(MATHLIB::scpr(xn, v_bm.getData(), 3) / sqr_radius1);
   //				if (fabs(a1 + a2 - a0) < tol)
   //					msh_nod_vector.push_back(nod_vector[i]->GetIndex());
   //			}
   //		}
   //	}
   //}

   /**************************************************************************
    FEMLib-Method:
    Task: Ermittelt den nahliegenden existierenden Knoten
    Programing:
    03/2005 OK Implementation (based on ExecuteSourceSinkMethod11 by CT)
    07/2005 WW Node object is replaced
    10/2005 OK test
    **************************************************************************/
   void CFEMesh::GetNODOnPLY(CGLPolyline* m_ply, std::vector<long>&msh_nod_vector) const
   {
      if (m_ply->point_vector.size() == 0)
         return;

      //	_min_edge_length = m_ply->epsilon;

      long j, k, l;
      double pt1[3], line1[3], line2[3];
      double mult_eps = 1.0;
      double dist1p, dist2p;
      // 07/2010 TF declare and initialize
      double *length (new double[m_ply->point_vector.size()]);
      double laenge;
      long anz_relevant = 0;
      typedef struct
      {
         long knoten;
         long abschnitt;
         double laenge;
      } INFO;
      INFO *relevant = NULL;

      // 07/2010 TF - variables only for swapping
      //	long knoten_help;
      //	double laenge_help;
      m_ply->getSBuffer().clear();
      m_ply->getIBuffer().clear();
      msh_nod_vector.clear();
      //
      /* */
      for (k = 0; k < (long) m_ply->point_vector.size() - 1; k++)
      {
         line1[0] = m_ply->point_vector[k]->x;
         line1[1] = m_ply->point_vector[k]->y;
         line1[2] = m_ply->point_vector[k]->z;
         line2[0] = m_ply->point_vector[k + 1]->x;
         line2[1] = m_ply->point_vector[k + 1]->y;
         line2[2] = m_ply->point_vector[k + 1]->z;
         length[k] = MCalcDistancePointToPoint(line2, line1);
      }

      /* Wiederholen bis zumindest ein Knoten gefunden wurde */
      while (anz_relevant == 0)
      {
         /* Schleife ueber alle Knoten */
         for (j = 0; j < NodesInUsage(); j++)
         {
            /* Schleife ueber alle Punkte des Polygonzuges */
            for (k = 0; k < (long) m_ply->point_vector.size() - 1; k++)
            {
               /* ??? */
               pt1[0] = nod_vector[j]->X();
               pt1[1] = nod_vector[j]->Y();
               pt1[2] = nod_vector[j]->Z();
               line1[0] = m_ply->point_vector[k]->x;
               line1[1] = m_ply->point_vector[k]->y;
               line1[2] = m_ply->point_vector[k]->z;
               line2[0] = m_ply->point_vector[k + 1]->x;
               line2[1] = m_ply->point_vector[k + 1]->y;
               line2[2] = m_ply->point_vector[k + 1]->z;
               /* Ist der Knoten nah am Polygonabschnitt? */
               if (MCalcDistancePointToLine(pt1, line1, line2) <= mult_eps
                  * m_ply->epsilon)
               {
                  /* Im folgenden wird mit der Projektion weitergearbeitet */
                  MCalcProjectionOfPointOnLine(pt1, line1, line2, pt1);
                  /* Abstand des Punktes zum ersten Punkt des Polygonabschnitts */
                  dist1p = MCalcDistancePointToPoint(line1, pt1);
                  /* Abstand des Punktes zum zweiten Punkt des Polygonabschnitts */
                  dist2p = MCalcDistancePointToPoint(line2, pt1);
                  /* Ist der Knoten innerhalb des Intervalls? */
                  /* bis rf3807: if ((length[k] - dist1p - dist2p + MKleinsteZahl)/(length[k] + dist1p + dist2p + MKleinsteZahl) > -MKleinsteZahl){ */
                  if ((dist1p + dist2p - length[k]) <= mult_eps
                     * m_ply->epsilon)
                  {
                     // For boundara conditions. WW
                     m_ply->getSBuffer().push_back(dist1p);
                                                  //Section index
                     m_ply->getIBuffer().push_back(k);
                     anz_relevant++;
                     /* Feld anpassen */
                     //nodes_all = (long *) Realloc(nodes_all,sizeof(long)*anz_relevant);
                     relevant = (INFO *) Realloc(relevant, sizeof(INFO)
                        * anz_relevant);
                     /* Ablegen von Knotennummer und Position */
                     //nodes_all[anz_relevant-1] = j;
                     msh_nod_vector.push_back(nod_vector[j]->GetIndex());
                     /* Position ermitteln */
                     laenge = 0.;
                     for (l = 0; l < k; l++)
                        laenge += length[l];
                     /* Ablegen von Knotennummer und Position */
                     relevant[anz_relevant - 1].knoten = j;
                     relevant[anz_relevant - 1].laenge = laenge + dist1p;
                     /* Suche am Polygon abbrechen, naechster Knoten */
                     k = (long) m_ply->point_vector.size();
                  }
               }                                  /* endif */
            }                                     /* Ende Schleife ueber Polygonabschnitte */
         }                                        /* Ende Schleife ueber Knoten */
         if (anz_relevant == 0)
            mult_eps *= 2.;
      }                                           /* Ende Schleife Wiederholungen */
      if (mult_eps > 1.)
         std::cout << "!!! Epsilon increased in sources!" << std::endl;

      delete [] length;

      //	/* Schleife ueber alle Knoten; sortieren nach Reihenfolge auf dem Abschnitt (zyklisches Vertauschen, sehr lahm)*/
      //	int weiter;
      //	double w1, w2;
      //	do {
      //		weiter = 0;
      //		for (k = 0; k < anz_relevant - 1; k++) {
      //			w1 = relevant[k].laenge;
      //			w2 = relevant[k + 1].laenge;
      //			if (w1 > w2) { /* Die Eintraege vertauschen */
      //				// 07/2010 TF
      //				std::swap (relevant[k].knoten, relevant[k+1].knoten);
      //				std::swap (relevant[k].laenge, relevant[k+1].laenge);
      ////				knoten_help = relevant[k].knoten;
      ////				laenge_help = relevant[k].laenge;
      ////				relevant[k].knoten = relevant[k + 1].knoten;
      ////				relevant[k].laenge = relevant[k + 1].laenge;
      ////				relevant[k + 1].knoten = knoten_help;
      ////				relevant[k + 1].laenge = laenge_help;
      //				weiter = 1;
      //			}
      //		}
      //	} while (weiter);
      relevant = (INFO*) Free(relevant);

      //WW
      m_ply->GetPointOrderByDistance();
      //----------------------------------------------------------------------
#ifdef MSH_CHECK
      cout << "MSH nodes at polyline:" << endl;
      for(j=0;j<(int)msh_nod_vector.size();j++)
      {
         cout << msh_nod_vector [j] << endl;
      }
#endif
      //----------------------------------------------------------------------
   }

   /**************************************************************************
    FEMLib-Method:
    Task: Ermittelt den nahliegenden existierenden Knoten
    Programing:
    03/2005 OK Implementation (based on ExecuteSourceSinkMethod11 by CT)
    07/2005 WW Node object is replaced
    10/2005 OK test
    03/2010 TF adaption to new data GEO-structures, changed the algorithm
    **************************************************************************/
   void CFEMesh::GetNODOnPLY(const GEOLIB::Polyline* const ply,
      std::vector<size_t>& msh_nod_vector)
   {
      msh_nod_vector.clear();

      // search for nodes along polyline in previous computed polylines
      std::vector<MeshNodesAlongPolyline>::const_iterator it (_mesh_nodes_along_polylines.begin());
      for (; it != _mesh_nodes_along_polylines.end(); it++)
      {
         if (it->getPolyline() == ply)
         {
            const std::vector<size_t> node_ids (it->getNodeIDs ());

            size_t n_valid_nodes (0);
            if (useQuadratic)
               n_valid_nodes = node_ids.size();
            else
               n_valid_nodes = it->getNumberOfLinearNodes();

            for (size_t k(0); k<n_valid_nodes; k++)
               msh_nod_vector.push_back (node_ids[k]);
            std::cout << "****** access " << msh_nod_vector.size() << " buffered nodes for polyline " << ply << std::endl;
            return;
         }
      }

      // compute nodes (and supporting points) along polyline
      _mesh_nodes_along_polylines.push_back (MeshNodesAlongPolyline (ply, this));
      const std::vector<size_t> node_ids (_mesh_nodes_along_polylines[_mesh_nodes_along_polylines.size()-1].getNodeIDs ());

      size_t n_valid_nodes (0);
      if (useQuadratic)
         n_valid_nodes = node_ids.size();
      else
         n_valid_nodes = _mesh_nodes_along_polylines[_mesh_nodes_along_polylines.size()-1].getNumberOfLinearNodes();

      for (size_t k(0); k<n_valid_nodes; k++)
         msh_nod_vector.push_back (node_ids[k]);
      std::cout << "****** computed " << n_valid_nodes << " nodes for polyline " << ply << " - " << NodesInUsage() << std::endl;
   }

   const MeshNodesAlongPolyline& CFEMesh::GetMeshNodesAlongPolyline(const GEOLIB::Polyline* const ply)
   {
      // search for nodes along polyline in previous computed polylines
      std::vector<MeshNodesAlongPolyline>::const_iterator it (_mesh_nodes_along_polylines.begin());
      for (; it != _mesh_nodes_along_polylines.end(); it++)
      {
         if (it->getPolyline() == ply)
         {
            return *it;
         }
      }
      // compute nodes (and supporting points for interpolation) along polyline
      _mesh_nodes_along_polylines.push_back (MeshNodesAlongPolyline (ply, this));
      return _mesh_nodes_along_polylines[_mesh_nodes_along_polylines.size()-1];
   }

   void CFEMesh::getPointsForInterpolationAlongPolyline (const GEOLIB::Polyline* const ply, std::vector<double>& points)
   {
      // search for nodes along polyline in previous computed polylines
      std::vector<MeshNodesAlongPolyline>::const_iterator it (_mesh_nodes_along_polylines.begin());
      for (; it != _mesh_nodes_along_polylines.end(); it++)
      {
         if (it->getPolyline() == ply)
         {
            // copy points from object into vector
            for (size_t k(0); k<it->getDistOfProjNodeFromPlyStart().size(); k++)
               points.push_back (it->getDistOfProjNodeFromPlyStart()[k]);
            return;
         }
      }

      // compute nodes (and points according the nodes) along polyline
      _mesh_nodes_along_polylines.push_back (MeshNodesAlongPolyline (ply, this));
      // copy supporting points from object into vector
      for (size_t k(0); k<(_mesh_nodes_along_polylines.back()).getDistOfProjNodeFromPlyStart().size(); k++)
         points.push_back (it->getDistOfProjNodeFromPlyStart()[k]);
   }

   void CFEMesh::GetNODOnPLY(const GEOLIB::Polyline* const ply, std::vector<long>& msh_nod_vector)
   {
      std::vector<size_t> tmp_msh_node_vector;
      GetNODOnPLY (ply, tmp_msh_node_vector);
      for (size_t k(0); k<tmp_msh_node_vector.size(); k++)
         msh_nod_vector.push_back (tmp_msh_node_vector[k]);
   }

   /**************************************************************************
    MSHLib-Method:
    Task:
    Programing:
    04/2005 OK
    last modification:
    **************************************************************************/
   void CFEMesh::GetNODOnSFC(Surface*m_sfc, std::vector<long>&msh_nod_vector)
   {
      msh_nod_vector.clear();
      //----------------------------------------------------------------------
      switch (m_sfc->type)
      {
         //....................................................................
         case 0:                                  // Surface polygon
            GetNODOnSFC_PLY(m_sfc, msh_nod_vector);
            break;
         case 1:                                  // TIN
            if (!m_sfc->TIN)
            {
               return;
            }
            GetNODOnSFC_TIN(m_sfc, msh_nod_vector);
            break;
            //....................................................................
         case 2:                                  // 2 vertical polylines
            GetNODOnSFC_Vertical(m_sfc,msh_nod_vector);
            break;
         case 3:                                  // projection on xy plane (all mesh points above and below the surface) //MB
            GetNODOnSFC_PLY_XY(m_sfc, msh_nod_vector);
            break;
            //....................................................................
         case 100:
            GetNodesOnCylindricalSurface(m_sfc, msh_nod_vector);
            break;
         case 4:                                  // layer polyline, all z
            GetNODOnSFC_PLY_Z(m_sfc,msh_nod_vector);
            break;
      }
   }

   /**************************************************************************
    MSHLib-Method:
    Task: Get nodes on plane surface
    Programing:
    03/2010 TF
    last modification:
    **************************************************************************/
   void CFEMesh::GetNODOnSFC(const GEOLIB::Surface* sfc, std::vector<size_t>& msh_nod_vector) const
   {
      msh_nod_vector.clear();

      for (size_t j(0); j<(size_t)NodesInUsage(); j++)
      {
         if (sfc->isPntInBV ( (nod_vector[j])->getData() ))
         {
            if (! sfc->isPntInSfc((nod_vector[j])->getData()))
               msh_nod_vector.push_back (nod_vector[j]->GetIndex());
         }
      }
   }

   /**************************************************************************
    MSHLib-Method:
    Task: Get nodes on plane surface by comparing the area of polygon computed
    by triangles, which are formed by node and the gravity center
    with edges of polygon, respectively
    Programing:
    09/2004 WW Implementation
    04/2005 OK MSH
    07/2005 WW Node object is replaced
    last modification:
    **************************************************************************/
   void CFEMesh::GetNODOnSFC_PLY(Surface*m_sfc, std::vector<long>&msh_nod_vector)
   {
      long i, j, k;
      int nPointsPly = 0;
      double gC[3], p1[3], p2[3];
      double Area1, Area2;
      double Tol = m_sfc->epsilon;
      CGLPolyline* m_ply = NULL;
      std::vector<CGLPolyline*>::iterator p_ply;  //CC
      // Init
      msh_nod_vector.clear();
      //----------------------------------------------------------------------
      // nodes close to first polyline
                                                  //CC
      p_ply = m_sfc->polyline_of_surface_vector.begin();
                                                  //CC
      while (p_ply != m_sfc->polyline_of_surface_vector.end())
      {
         m_ply = *p_ply;
         nPointsPly = (int) m_ply->point_vector.size();
         //....................................................................
         // Gravity center of this polygon
         for (i = 0; i < 3; i++)
            gC[i] = 0.0;
         for (i = 0; i < nPointsPly; i++)
         {
            gC[0] += m_ply->point_vector[i]->x;
            gC[1] += m_ply->point_vector[i]->y;
            gC[2] += m_ply->point_vector[i]->z;
         }
         for (i = 0; i < 3; i++)
            gC[i] /= (double) nPointsPly;
         //....................................................................
         // Area of this polygon by the grativity center
         Area1 = 0.0;
         for (i = 0; i < nPointsPly; i++)
         {
            p1[0] = m_ply->point_vector[i]->x;
            p1[1] = m_ply->point_vector[i]->y;
            p1[2] = m_ply->point_vector[i]->z;
            if (i < nPointsPly - 1)
            {
               p2[0] = m_ply->point_vector[i + 1]->x;
               p2[1] = m_ply->point_vector[i + 1]->y;
               p2[2] = m_ply->point_vector[i + 1]->z;
            }
            else
            {
               p2[0] = m_ply->point_vector[0]->x;
               p2[1] = m_ply->point_vector[0]->y;
               p2[2] = m_ply->point_vector[0]->z;
            }
            Area1 += fabs(ComputeDetTri(p1, gC, p2));
         }
         //....................................................................
         // Check nodes by comparing area
         for (j = 0; j < NodesInUsage(); j++)
         {
            Area2 = 0.0;
            gC[0] = nod_vector[j]->X();
            gC[1] = nod_vector[j]->Y();
            gC[2] = nod_vector[j]->Z();
            for (i = 0; i < nPointsPly; i++)
            {
               p1[0] = m_ply->point_vector[i]->x;
               p1[1] = m_ply->point_vector[i]->y;
               p1[2] = m_ply->point_vector[i]->z;
               k = i + 1;
               if (i == nPointsPly - 1)
                  k = 0;
               p2[0] = m_ply->point_vector[k]->x;
               p2[1] = m_ply->point_vector[k]->y;
               p2[2] = m_ply->point_vector[k]->z;
               Area2 += fabs(ComputeDetTri(p1, gC, p2));
            }
            if (fabs(Area1 - Area2) < Tol)
               msh_nod_vector.push_back(nod_vector[j]->GetIndex());
         }
         p_ply++;
      }
   }

   /**************************************************************************
    MSHLib-Method:
    Task: Get nodes on plane surface by comparing the area of polygon computed
    by triangles, which are formed by node and the gravity center
    with edges of polygon, respectively
    Programing:
    08/2005 MB based on GetNODOnSFC_PLY
    03/2009 WW Case only search the specified nodes given in msh_nod_vector
    03/2009 WW Efficiency improvement
    last modification:
    **************************************************************************/
   void CFEMesh::GetNODOnSFC_PLY_XY(Surface*m_sfc, std::vector<long>&msh_nod_vector,
      bool givenNodesOnSurface)
   {
      long i, j, k;
      int nPointsPly = 0;
      double gC[3], p1[3], p2[3];
      double Area1, Area2;
      double Tol = m_sfc->epsilon;
      CGLPolyline* m_ply = NULL;
      std::vector<CGLPolyline*>::iterator p_ply;  //CC
      // Init
      //----19.03.2009. WW
      CNode *a_node = NULL;
      std:: vector<long> temp_v;
      //
      if (givenNodesOnSurface)
      {
         temp_v.resize((long) msh_nod_vector.size());
         temp_v = msh_nod_vector;
      }
      p1[2] = p2[2] = 0.;
      //----19.03.2009. WW
      //
      msh_nod_vector.clear();
      //----------------------------------------------------------------------
      // nodes close to first polyline
                                                  //CC
      p_ply = m_sfc->polyline_of_surface_vector.begin();
                                                  //CC
      while (p_ply != m_sfc->polyline_of_surface_vector.end())
      {
         m_ply = *p_ply;
         nPointsPly = (int) m_ply->point_vector.size();
         //....................................................................
         // Grativity center of this polygon
         for (i = 0; i < 3; i++)
            gC[i] = 0.0;
         for (i = 0; i < nPointsPly; i++)
         {
            gC[0] += m_ply->point_vector[i]->x;
            gC[1] += m_ply->point_vector[i]->y;
         }
         for (i = 0; i < 3; i++)
            gC[i] /= (double) nPointsPly;
         //....................................................................
         // Area of this polygon by the grativity center
         Area1 = 0.0;
         for (i = 0; i < nPointsPly; i++)
         {
            p1[0] = m_ply->point_vector[i]->x;
            p1[1] = m_ply->point_vector[i]->y;
            k = i + 1;
            if (i == nPointsPly - 1)
               k = 0;
            p2[0] = m_ply->point_vector[k]->x;
            p2[1] = m_ply->point_vector[k]->y;
            Area1 += fabs(ComputeDetTri(p1, gC, p2));
         }
         //....................................................................
         // Check nodes by comparing area
         //------- 19.03.2009. WW ----------------
         if (givenNodesOnSurface)
         {
            for (j = 0; j < (long) temp_v.size(); j++)
            {
               Area2 = 0.0;
               a_node = nod_vector[temp_v[j]];
               gC[0] = a_node->X();
               gC[1] = a_node->Y();
               for (i = 0; i < nPointsPly; i++)
               {
                  p1[0] = m_ply->point_vector[i]->x;
                  p1[1] = m_ply->point_vector[i]->y;
                  k = i + 1;
                  if (i == nPointsPly - 1)
                     k = 0;
                  p2[0] = m_ply->point_vector[k]->x;
                  p2[1] = m_ply->point_vector[k]->y;
                  Area2 += fabs(ComputeDetTri(p1, gC, p2));
               }
               if (fabs(Area1 - Area2) < Tol)
                  msh_nod_vector.push_back(a_node->GetIndex());
            }
         }
         //-----------------------------------------------
         else
         {
            for (j = 0; j < NodesInUsage(); j++)
            {
               Area2 = 0.0;
               a_node = nod_vector[j];            //19.03.2009. WW

               gC[0] = a_node->X();
               gC[1] = a_node->Y();
               for (i = 0; i < nPointsPly; i++)
               {
                  p1[0] = m_ply->point_vector[i]->x;
                  p1[1] = m_ply->point_vector[i]->y;
                  k = i + 1;
                  if (i == nPointsPly - 1)
                     k = 0;
                  p2[0] = m_ply->point_vector[k]->x;
                  p2[1] = m_ply->point_vector[k]->y;
                  Area2 += fabs(ComputeDetTri(p1, gC, p2));
               }
               if (fabs(Area1 - Area2) < Tol)
                  msh_nod_vector.push_back(a_node->GetIndex());
            }
         }
         p_ply++;
      }
      //
      if (givenNodesOnSurface)                    //19.03.2009. WW
         temp_v.clear();
   }

   /**************************************************************************
    MSHLib-Method:
    Task:
    Programing:
    04/2005 OK
    07/2005 WW Node object is replaced
    04/2006 TK new method
    01/2010 NW use epsilon specified in GEO as tolerance
    last modification:
    **************************************************************************/
   void CFEMesh::GetNODOnSFC_TIN(Surface*m_sfc, std::vector<long>&msh_nod_vector)
   {
      long i = 0, k = 0, m = 0;
      double angle_sum, dist;
      double tolerance = 0.001;
      double min_mesh_dist = 0.0;
      double tri_point1[3], tri_point2[3], tri_point3[3], checkpoint[3];
      double sfc_min[3] = {m_sfc->TIN->Triangles[0]->x[0], m_sfc->TIN->Triangles[0]->y[0], m_sfc->TIN->Triangles[0]->z[0] };
      double sfc_max[3] = {m_sfc->TIN->Triangles[0]->x[0], m_sfc->TIN->Triangles[0]->y[0], m_sfc->TIN->Triangles[0]->z[0] };

      CTriangle *m_triangle = NULL;
      //----------------------------------------------------------------------
      // Create Bounding BOX = MIN/MAX of X/Y/Z
      //----------------------------------------------------------------------
      //Loop over all generated triangles of surface
      for (m = 0; m < (long) m_sfc->TIN->Triangles.size(); m++)
      {
         m_triangle = m_sfc->TIN->Triangles[m];
         tri_point1[0] = m_triangle->x[0];
         tri_point1[1] = m_triangle->y[0];
         tri_point1[2] = m_triangle->z[0];
         tri_point2[0] = m_triangle->x[1];
         tri_point2[1] = m_triangle->y[1];
         tri_point2[2] = m_triangle->z[1];
         tri_point3[0] = m_triangle->x[2];
         tri_point3[1] = m_triangle->y[2];
         tri_point3[2] = m_triangle->z[2];
         if (m == 0)
         {
            sfc_min[0] = tri_point1[0];
            sfc_min[1] = tri_point1[1];
            sfc_min[2] = tri_point1[2];
            sfc_max[0] = tri_point1[0];
            sfc_max[1] = tri_point1[1];
            sfc_max[2] = tri_point1[2];
            if (tri_point1[0] < sfc_min[0])
               sfc_min[0] = tri_point1[0];
            if (tri_point2[0] < sfc_min[0])
               sfc_min[0] = tri_point2[0];
            if (tri_point3[0] < sfc_min[0])
               sfc_min[0] = tri_point3[0];
            if (tri_point1[0] > sfc_max[0])
               sfc_max[0] = tri_point1[0];
            if (tri_point2[0] > sfc_max[0])
               sfc_max[0] = tri_point2[0];
            if (tri_point3[0] > sfc_max[0])
               sfc_max[0] = tri_point3[0];
            if (tri_point1[1] < sfc_min[1])
               sfc_min[1] = tri_point1[1];
            if (tri_point2[1] < sfc_min[1])
               sfc_min[1] = tri_point2[1];
            if (tri_point3[1] < sfc_min[1])
               sfc_min[1] = tri_point3[1];
            if (tri_point1[1] > sfc_max[1])
               sfc_max[1] = tri_point1[1];
            if (tri_point2[1] > sfc_max[1])
               sfc_max[1] = tri_point2[1];
            if (tri_point3[1] > sfc_max[1])
               sfc_max[1] = tri_point3[1];
            if (tri_point1[2] < sfc_min[2])
               sfc_min[2] = tri_point1[2];
            if (tri_point2[2] < sfc_min[2])
               sfc_min[2] = tri_point2[2];
            if (tri_point3[2] < sfc_min[2])
               sfc_min[2] = tri_point3[2];
            if (tri_point1[2] > sfc_max[2])
               sfc_max[2] = tri_point1[2];
            if (tri_point2[2] > sfc_max[2])
               sfc_max[2] = tri_point2[2];
            if (tri_point3[2] > sfc_max[2])
               sfc_max[2] = tri_point3[2];
         }
         else
         {
            if (tri_point1[0] < sfc_min[0])
               sfc_min[0] = tri_point1[0];
            if (tri_point2[0] < sfc_min[0])
               sfc_min[0] = tri_point2[0];
            if (tri_point3[0] < sfc_min[0])
               sfc_min[0] = tri_point3[0];
            if (tri_point1[0] > sfc_max[0])
               sfc_max[0] = tri_point1[0];
            if (tri_point2[0] > sfc_max[0])
               sfc_max[0] = tri_point2[0];
            if (tri_point3[0] > sfc_max[0])
               sfc_max[0] = tri_point3[0];
            if (tri_point1[1] < sfc_min[1])
               sfc_min[1] = tri_point1[1];
            if (tri_point2[1] < sfc_min[1])
               sfc_min[1] = tri_point2[1];
            if (tri_point3[1] < sfc_min[1])
               sfc_min[1] = tri_point3[1];
            if (tri_point1[1] > sfc_max[1])
               sfc_max[1] = tri_point1[1];
            if (tri_point2[1] > sfc_max[1])
               sfc_max[1] = tri_point2[1];
            if (tri_point3[1] > sfc_max[1])
               sfc_max[1] = tri_point3[1];
            if (tri_point1[2] < sfc_min[2])
               sfc_min[2] = tri_point1[2];
            if (tri_point2[2] < sfc_min[2])
               sfc_min[2] = tri_point2[2];
            if (tri_point3[2] < sfc_min[2])
               sfc_min[2] = tri_point3[2];
            if (tri_point1[2] > sfc_max[2])
               sfc_max[2] = tri_point1[2];
            if (tri_point2[2] > sfc_max[2])
               sfc_max[2] = tri_point2[2];
            if (tri_point3[2] > sfc_max[2])
               sfc_max[2] = tri_point3[2];
         }
      }
      //----------------------------------------------------------------------
      // Create Local Search Vector
      // Only nodes inside searching box
      //----------------------------------------------------------------------

      CFEMesh* m_msh_aux (new CFEMesh(_geo_obj,_geo_name));
      CNode* node = NULL;

      tolerance = m_sfc->epsilon;                 //NW
      // NW commented out below. Minimum edge length doesn't work for some cases
      ////Loop over all edges
      //     for(i=0;i<(long)edge_vector.size();i++)
      //     {
      //         if (j==0 && i==0){
      //           min_mesh_dist = edge_vector[i]->Length();
      //         }
      //         else{
      //           if (min_mesh_dist  > edge_vector[i]->Length())
      //               min_mesh_dist =  edge_vector[i]->Length();
      //         }
      //     }
      //     tolerance = min_mesh_dist;

      // TF expand bounding box with epsilon environment
      for (size_t k(0); k<3; k++) {
    	  sfc_min[k] -= tolerance;
    	  sfc_max[k] += tolerance;
      }

      //Loop over all mesh nodes
      for (i = 0; i < NodesInUsage(); i++)        //NW cannot use nod_vector.size() because of higher order elements
      {
         checkpoint[0] = nod_vector[i]->X();
         checkpoint[1] = nod_vector[i]->Y();
         checkpoint[2] = nod_vector[i]->Z();
         node = new CNode(i, checkpoint[0], checkpoint[1], checkpoint[2]);
         if ((checkpoint[0] >= sfc_min[0] && checkpoint[0] <= sfc_max[0])
        		 && (checkpoint[1] >= sfc_min[1] && checkpoint[1] <= sfc_max[1])
        		 && (checkpoint[2] >= sfc_min[2] && checkpoint[2] <= sfc_max[2]))
         {
            m_msh_aux->nod_vector.push_back(node);
         }
      }

      //----------------------------------------------------------------------
      // Search preselected Nodes within TIN Triangles
      //----------------------------------------------------------------------
      for (m = 0; m < (long) m_sfc->TIN->Triangles.size(); m++)
      {
         m_triangle = m_sfc->TIN->Triangles[m];
         tri_point1[0] = m_triangle->x[0];
         tri_point1[1] = m_triangle->y[0];
         tri_point1[2] = m_triangle->z[0];
         tri_point2[0] = m_triangle->x[1];
         tri_point2[1] = m_triangle->y[1];
         tri_point2[2] = m_triangle->z[1];
         tri_point3[0] = m_triangle->x[2];
         tri_point3[1] = m_triangle->y[2];
         tri_point3[2] = m_triangle->z[2];
         //Loop over all preselected mesh nodes
         for (i = 0; i < (long) m_msh_aux->nod_vector.size(); i++)
         {
            checkpoint[0] = m_msh_aux->nod_vector[i]->X();
            checkpoint[1] = m_msh_aux->nod_vector[i]->Y();
            checkpoint[2] = m_msh_aux->nod_vector[i]->Z();
            dist = MCalcDistancePointToPlane(checkpoint, tri_point1,
               tri_point2, tri_point3);
            if (k == 0)
               m_msh_aux->nod_vector[i]->epsilon = dist;
            else
            {
               if (m_msh_aux->nod_vector[i]->epsilon > dist)
                  m_msh_aux->nod_vector[i]->epsilon = dist;
            }
            if (dist <= tolerance && dist >= -tolerance)
            {
               angle_sum = AngleSumPointInsideTriangle(checkpoint, tri_point1,
                  tri_point2, tri_point3, min_mesh_dist);
               if (angle_sum > 359)
                  m_msh_aux->nod_vector[i]->selected = 1;
            }
         }
      }

      //----------------------------------------------------------------------
      // Identify the preselected nodes of the search vector and copy to msh_nod_vector
      // TODO: Works only for one mesh!!!
      //----------------------------------------------------------------------
      int index;
      //Loop over selected nodes
      for (i = 0; i < (long) m_msh_aux->nod_vector.size(); i++)
      {
         index = m_msh_aux->nod_vector[i]->GetIndex();
         if (index < (int) nod_vector.size())
         {
            if ((m_msh_aux->nod_vector[i]->GetIndex()
               == nod_vector[index]->GetIndex())
               && m_msh_aux->nod_vector[i]->selected == 1
               && (m_msh_aux->nod_vector[i]->X() == nod_vector[index]->X())
               && (m_msh_aux->nod_vector[i]->Y() == nod_vector[index]->Y())
               && (m_msh_aux->nod_vector[i]->Z() == nod_vector[index]->Z()))
            {
               msh_nod_vector.push_back(nod_vector[index]->GetIndex());
            }
         }
      }
      //----------------------------------------------------------------------
      // Delete Search Vector at the end of fem_msh_vector
      // TODO: Proper delete by MSHDelete!!!
      //----------------------------------------------------------------------
      for (i = 0; i < (long) m_msh_aux->nod_vector.size(); i++)
      {
         delete m_msh_aux->nod_vector[i];
      }
   }

   /**************************************************************************
    GeoLib-Method:
    Task: Get nodes on cylindrical surface by comparing the area of
    triangles form by nodes and two axis points
    Programing:
    10/2004 WW Implementation
    05/2005 WW Transplant to this object from GEOLIB
    04/2006 WW Case of quadratic elements
    last modification:
    **************************************************************************/
   void CFEMesh::GetNodesOnCylindricalSurface(Surface*m_sfc, std::vector<long>& NodesS)
   {
      int k, l, nf;
      long i, j, m, fnode;
      const int nNodes = NodesInUsage();
      int faceIndex_loc[10];
      double gC[3], p1[3], p2[3];
      double dist, R, dc1, dc2;
      CElem* elem = NULL;
      CNode* cnode = NULL;
      NodesS.clear();
      //m_sfc->epsilon = 1.0e-6;
      p1[0] = m_sfc->polygon_point_vector[0]->x;
      p1[1] = m_sfc->polygon_point_vector[0]->y;
      p1[2] = m_sfc->polygon_point_vector[0]->z;

      p2[0] = m_sfc->polygon_point_vector[1]->x;
      p2[1] = m_sfc->polygon_point_vector[1]->y;
      p2[2] = m_sfc->polygon_point_vector[1]->z;

      dist = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1]
         - p2[1]) + (p1[2] - p2[2]) * (p1[2] - p2[2]));

      // Check nodes by comparing area
      for (j = 0; j < nNodes; j++)
      {
         cnode = nod_vector[j];
         gC[0] = cnode->X();
         gC[1] = cnode->Y();
         gC[2] = cnode->Z();

         dc1 = (p2[0] - p1[0]) * (gC[0] - p1[0]) + (p2[1] - p1[1]) * (gC[1]
            - p1[1]) + (p2[2] - p1[2]) * (gC[2] - p1[2]);
         dc2 = (p2[0] - p1[0]) * (gC[0] - p2[0]) + (p2[1] - p1[1]) * (gC[1]
            - p2[1]) + (p2[2] - p1[2]) * (gC[2] - p2[2]);
         if (dc1 < 0.0)
            continue;
         if (dc2 > 0.0)
            continue;

         R = 2.0 * fabs(ComputeDetTri(p1, gC, p2)) / dist;

         if (fabs(R - m_sfc->Radius) < m_sfc->epsilon)
            NodesS.push_back(cnode->GetIndex());
      }
      bool done = false;
      int counter = 0;
      int hs;
      long NodesS_size = (long) NodesS.size();
      //
      if (useQuadratic)
      {
         // Face elements are only in quadrilaterals or triangles
         for (i = 0; i < NodesS_size; i++)
         {
            cnode = nod_vector[NodesS[i]];
            for (j = 0; j < (long) cnode->connected_elements.size(); j++)
            {
               elem = ele_vector[cnode->connected_elements[j]];
               for (k = 0; k < elem->GetFacesNumber(); k++)
               {
                  nf = elem->GetElementFaceNodes(k, faceIndex_loc);
                  counter = 0;
                  hs = (int) (nf / 2);
                  for (l = 0; l < hs; l++)        // loop over face vertices
                  {
                     fnode = elem->GetNodeIndex(faceIndex_loc[l]);
                     for (m = 0; m < NodesS_size; m++)
                     {
                        if (fnode == NodesS[m])
                           counter++;
                     }
                  }
                  //
                  if (counter == hs)              // face is found on surface
                  {
                     for (l = hs; l < nf; l++)    // loop over face vertices
                     {
                        fnode = elem->GetNodeIndex(faceIndex_loc[l]);
                        done = false;
                        for (m = 0; m < (long) NodesS.size(); m++)
                        {
                           if (fnode == NodesS[m])
                           {
                              done = true;
                              break;
                           }
                        }
                        if (!done)
                           NodesS.push_back(fnode);
                     }
                  }
               }
            }
         }
      }
   }

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   04/2005 OK
   last modification:
   07/2010 TF small modifications concerning coding style
   **************************************************************************/
   void CFEMesh::GetNODOnSFC_Vertical(Surface*m_sfc, std::vector<long>&msh_nod_vector)
   {
      long *nodes_array = NULL;
      long no_nodes = 0;

      // nodes close to first polyline
      std::vector<CGLPolyline*>::const_iterator p_ply (m_sfc->polyline_of_surface_vector.begin());
      while (p_ply != m_sfc->polyline_of_surface_vector.end())
      {
         //OK41 nodes_array = m_polyline->MSHGetNodesCloseXY(&no_nodes);
                                                  //CC 10/05
         nodes_array = MSHGetNodesClose(&no_nodes, *p_ply);
         break;
      }

      // using triangles
      CGLPolyline* m_polyline1 = NULL;
      CGLPolyline* m_polyline2 = NULL;
      p_ply = m_sfc->polyline_of_surface_vector.begin();
      while (p_ply != m_sfc->polyline_of_surface_vector.end())
      {
         m_polyline1 = *p_ply;
         ++p_ply;
         m_polyline2 = *p_ply;
         break;
      }
      long no_points = (long) m_polyline1->point_vector.size();

      CGLPoint m_node;
      double xp[3], yp[3], zp[3];

      long i, j;
      for (j = 0; j < no_nodes; j++)
      {
         //OK m_node.x = GetNodeX(nodes_array[j]);
         //OK m_node.y = GetNodeY(nodes_array[j]);
         //OK m_node.z = GetNodeZ(nodes_array[j]);
         for (i = 0; i < no_points - 1; i++)
         {
            // first triangle of quad
            xp[0] = m_polyline1->point_vector[i]->x;
            yp[0] = m_polyline1->point_vector[i]->y;
            zp[0] = m_polyline1->point_vector[i]->z;
            xp[1] = m_polyline1->point_vector[i + 1]->x;
            yp[1] = m_polyline1->point_vector[i + 1]->y;
            zp[1] = m_polyline1->point_vector[i + 1]->z;
            xp[2] = m_polyline2->point_vector[i]->x;
            yp[2] = m_polyline2->point_vector[i]->y;
            zp[2] = m_polyline2->point_vector[i]->z;
                                                  //CC 10/05
            if (m_node.IsInsideTriangle(xp, yp, zp))
            {
               msh_nod_vector.push_back(nodes_array[j]);
            }
            // second triangle of quad
            xp[0] = m_polyline2->point_vector[i]->x;
            yp[0] = m_polyline2->point_vector[i]->y;
            zp[0] = m_polyline2->point_vector[i]->z;
            xp[1] = m_polyline2->point_vector[i + 1]->x;
            yp[1] = m_polyline2->point_vector[i + 1]->y;
            zp[1] = m_polyline2->point_vector[i + 1]->z;
            xp[2] = m_polyline1->point_vector[i + 1]->x;
            yp[2] = m_polyline1->point_vector[i + 1]->y;
            zp[2] = m_polyline1->point_vector[i + 1]->z;
                                                  //CC 10/05
            if (m_node.IsInsideTriangle(xp, yp, zp))
            {
               msh_nod_vector.push_back(nodes_array[j]);
            }
         }                                        // no_points
      }                                           // no_nodes
   }
   /**************************************************************************
    FEMLib-Method:
    Task: Ermittelt den nahliegenden existierenden Knoten
    Programing:
    03/2005 OK Implementation (based on ExecuteSourceSinkMethod11 by CT)
    last modification:
    **************************************************************************/
   void CFEMesh::GetNODOnPLY_XY(CGLPolyline*m_ply, std::vector<long>&msh_nod_vector)
   {
      long j, k, l;
      double pt1[3], line1[3], line2[3], pt0[3];
      double mult_eps = 1.0;
      double dist1p, dist2p, *length, laenge;
      long anz_relevant = 0;
      typedef struct
      {
         long knoten;
         long abschnitt;
         double laenge;
      } INFO;
      INFO *relevant = NULL;
      int weiter;
      double w1, w2;
      long knoten_help;
      double laenge_help;
      m_ply->getSBuffer().clear();
      m_ply->getIBuffer().clear();
      msh_nod_vector.clear();
      //
      length = (double*) Malloc(sizeof(double)
         * (long) m_ply->point_vector.size());
      pt0[0] = m_ply->point_vector[0]->x;
      pt0[1] = m_ply->point_vector[0]->y;
      pt0[2] = 0.0;
      /* */
      for (k = 0; k < (long) m_ply->point_vector.size() - 1; k++)
      {
         line1[0] = m_ply->point_vector[k]->x;
         line1[1] = m_ply->point_vector[k]->y;
         line1[2] = 0.0;
         line2[0] = m_ply->point_vector[k + 1]->x;
         line2[1] = m_ply->point_vector[k + 1]->y;
         line2[2] = 0.0;
         length[k] = MCalcDistancePointToPoint(line2, line1);
      }
      /* Wiederholen bis zumindest ein Knoten gefunden wurde */
      while (anz_relevant == 0)
      {
         /* Schleife ueber alle Knoten */
         for (j = 0; j < (long) nod_vector.size(); j++)
         {
            /* Schleife ueber alle Punkte des Polygonzuges */
            for (k = 0; k < (long) m_ply->point_vector.size() - 1; k++)
            {
               /* ??? */
               pt1[0] = nod_vector[j]->X();
               pt1[1] = nod_vector[j]->Y();
               pt1[2] = 0.0;
               line1[0] = m_ply->point_vector[k]->x;
               line1[1] = m_ply->point_vector[k]->y;
               line1[2] = 0.0;
               line2[0] = m_ply->point_vector[k + 1]->x;
               line2[1] = m_ply->point_vector[k + 1]->y;
               line2[2] = 0.0;
               /* Ist der Knoten nah am Polygonabschnitt? */
               if (MCalcDistancePointToLine(pt1, line1, line2) <= mult_eps
                  * m_ply->epsilon)
               {
                  /* Im folgenden wird mit der Projektion weitergearbeitet */
                  MCalcProjectionOfPointOnLine(pt1, line1, line2, pt1);
                  /* Abstand des Punktes zum ersten Punkt des Polygonabschnitts */
                  dist1p = MCalcDistancePointToPoint(line1, pt1);
                  /* Abstand des Punktes zum zweiten Punkt des Polygonabschnitts */
                  dist2p = MCalcDistancePointToPoint(line2, pt1);
                  /* Ist der Knoten innerhalb des Intervalls? */
                  /* bis rf3807: if ((length[k] - dist1p - dist2p + MKleinsteZahl)/(length[k] + dist1p + dist2p + MKleinsteZahl) > -MKleinsteZahl){ */
                  if ((dist1p + dist2p - length[k]) <= mult_eps
                     * m_ply->epsilon)
                  {
                     // For boundara conditions. WW
                     m_ply->getSBuffer().push_back(dist1p);
                     m_ply->getIBuffer().push_back(k);
                     anz_relevant++;
                     /* Feld anpassen */
                     //nodes_all = (long *) Realloc(nodes_all,sizeof(long)*anz_relevant);
                     relevant = (INFO *) Realloc(relevant, sizeof(INFO)
                        * anz_relevant);
                     /* Ablegen von Knotennummer und Position */
                     //nodes_all[anz_relevant-1] = j;
                     msh_nod_vector.push_back(j);
                     /* Position ermitteln */
                     laenge = 0.;
                     for (l = 0; l < k; l++)
                        laenge += length[l];
                     /* Ablegen von Knotennummer und Position */
                     relevant[anz_relevant - 1].knoten = j;
                     relevant[anz_relevant - 1].laenge = laenge + dist1p;
                     /* Suche am Polygon abbrechen, naechster Knoten */
                     k = (long) m_ply->point_vector.size();
                  }
               }                                  /* endif */
            }                                     /* Ende Schleife ueber Polygonabschnitte */
         }                                        /* Ende Schleife ueber Knoten */
         if (anz_relevant == 0)
            mult_eps *= 2.;
      }                                           /* Ende Schleife Wiederholungen */
      if (mult_eps > 1.)
         std::cout << "!!! Epsilon increased in sources!" << std::endl;
      /* Schleife ueber alle Knoten; sortieren nach Reihenfolge auf dem Abschnitt (zyklisches Vertauschen, sehr lahm)*/
      do
      {
         weiter = 0;
         for (k = 0; k < anz_relevant - 1; k++)
         {
            w1 = relevant[k].laenge;
            w2 = relevant[k + 1].laenge;
            if (w1 > w2)                          /* Die Eintraege vertauschen */
            {
               knoten_help = relevant[k].knoten;
               laenge_help = relevant[k].laenge;
               relevant[k].knoten = relevant[k + 1].knoten;
               relevant[k].laenge = relevant[k + 1].laenge;
               relevant[k + 1].knoten = knoten_help;
               relevant[k + 1].laenge = laenge_help;
               weiter = 1;
            }
         }
      } while (weiter);
      relevant = (INFO*) Free(relevant);
   }

   /**************************************************************************
    MSHLib-Method:
    Task:
    Programing:
    09/2005 OK Implementation
    09/2005 OK Epsilon
    10/2005 OK Delete existing layer polylines
    02/2006 CC polyline id
    **************************************************************************/
   void CFEMesh::CreateLayerPolylines(CGLPolyline* m_ply)
   {
      long i;
      CGLPolyline* m_polyline = NULL;
      char layer_number[3];
      //---------------------------------------------------------------------
      // Delete existing layer polylines
      std::string ply_lay_name = m_ply->getName() + "_L";
      for (int k = 0; k < (int) polyline_vector.size(); k++)
      {
         m_polyline = polyline_vector[k];
         if (m_polyline->getName().find(ply_lay_name) != std::string::npos)
         {
            GEORemovePLY(m_polyline);
            //GEORemovePolyline(polyline_vector.begin()+(k-l));
            k--;
         }
      }
      //---------------------------------------------------------------------
      //
      std::vector<long> ply_nod_vector;
      std::vector<long> ply_nod_vector_dummy;
      GetNODOnPLY_XY(m_ply, ply_nod_vector);
      //nodes = MSHGetNodesCloseXY(&no_nodes); //OK41
      long nodes_per_layer = (long) nod_vector.size() / (_n_msh_layer + 1);
      int ply_nod_vector_layer = (int) ply_nod_vector.size() / (_n_msh_layer + 1);
      //---------------------------------------------------------------------
      // Create layer polylines
      //polyline id CC8888---------------------------------------------------
      long size = 0;
      CGLPolyline *ms_polyline = NULL;
      long number_of_polylines = (long) polyline_vector.size();
      if (number_of_polylines == 0)
         size = 0;
      else
      {
    	  std::vector<CGLPolyline*>::iterator ps = polyline_vector.begin();
         while (ps != polyline_vector.end())
         {
            ms_polyline = *ps;
            ++ps;
         }
         size = ms_polyline->getID() + 1;
      }

      if (ply_nod_vector_layer < 1)
         return;

      //---------------------------------------------------------------------
      // Create layer polylines
      //......................................................................
      m_polyline = new CGLPolyline;
      sprintf(layer_number, "%ld", 0L);
      //CString names =  m_ply->name + "_L" + layer_number;
      // m_polyline->name = names;

      // 10/2010 TF
      //	m_polyline->ply_name = m_ply->getName(); //CC/TK8888
      //	m_polyline->ply_name.append("_L"); //CC/TK8888
      //	m_polyline->ply_name.append(layer_number);//CC/TK8888
      std::string tmp_name (m_ply->getName() + "_L");
      tmp_name.append (layer_number);

      //	m_polyline->name = m_polyline->ply_name.data();//CC/TK8888 // TF
      m_polyline->setName (tmp_name);             // TF

      m_polyline->setDataType (1);
      m_polyline->setID(size);                    //CC8888 / TF
      m_polyline->epsilon = m_ply->epsilon;       //OK
      //	m_polyline->ply_data = m_polyline->getName () + ".ply";//CC

      for (i = 0; i < ply_nod_vector_layer; i++)
      {
         CGLPoint* point (new CGLPoint(nod_vector[ply_nod_vector[i]]->getData()));
         m_polyline->point_vector.push_back(point);
      }
      m_polyline->SetPointOrderByDistance(m_ply->point_vector[0]);
      polyline_vector.push_back(m_polyline);
      m_polyline->WritePointVector(m_polyline->getName());
      m_polyline->WriteTecplot(" ");              //OK41
      //......................................................................
      for (size_t j = 1; j < (_n_msh_layer + 1); j++)
      {
         m_polyline = new CGLPolyline;
         sprintf(layer_number, "%ld", j);
         //		m_polyline->ply_name = m_ply->name.data();//CC/TK8888
         //		m_polyline->ply_name.append("_L");//CC/TK8888
         //		m_polyline->ply_name.append(layer_number);//CC/TK8888
         //		m_polyline->name = m_polyline->ply_name.data();//CC/TK8888

         std::string tmp_name (m_ply->getName() + "_L");
         tmp_name.append (layer_number);
         m_polyline->setName (tmp_name);

         m_polyline->setDataType (1);
                                                  //OK
         m_polyline->epsilon = _min_edge_length / 2.;
         for (i = 0; i < ply_nod_vector_layer; i++)
         {
            CGLPoint* point (new CGLPoint (nod_vector[ply_nod_vector[i] + j * nodes_per_layer]->getData()));
            m_polyline->point_vector.push_back(point);
         }
         //OK    m_polyline->SortPointVectorByDistance();
         m_polyline->SetPointOrderByDistance(m_ply->point_vector[0]);
         polyline_vector.push_back(m_polyline);
         m_polyline->WritePointVector(m_polyline->getName ());
         m_polyline->WriteTecplot(" ");           //OK41
      }
   }

   /**************************************************************************
    MSHLib-Method:
    Task: Copies the selected nodes to a msh_node_vector
    Programing:
    12/2005 TK Implementation
    **************************************************************************/
   void CFEMesh::CopySelectedNodes(std::vector<long>&msh_nod_vector)
   {
      int i = 0, j = 0;
      // Init
      msh_nod_vector.clear();

      //Loop over all meshes
      for (j = 0; j < (long) fem_msh_vector.size(); j++)
      {
         //Loop over all mesh nodes
         for (i = 0; i < (long) fem_msh_vector[j]->nod_vector.size(); i++)
         {
            if (fem_msh_vector[j]->nod_vector[i]->selected == 1)
            {
               msh_nod_vector.push_back(
                  fem_msh_vector[j]->nod_vector[i]->GetIndex());
            }
         }
      }
   }

   /**************************************************************************
    MSHLib-Method:
    08/2006 OK Implementation
    03/2010 TF change to new data structures, changed algorithm
    **************************************************************************/
   void CFEMesh::GetELEOnPLY(const GEOLIB::Polyline* ply, std::vector<size_t>& ele_vector_ply)
   {
      vec<CEdge*> ele_edges_vector(15);
      vec<CNode*> edge_nodes(3);

      std::vector<size_t> nodes_near_ply;

      // get mesh nodes near the polyline
      GetNODOnPLY(ply, nodes_near_ply);

      // clear the given vector
      ele_vector_ply.clear();

      // loop over all elements
      for (size_t i=0; i<ele_vector.size(); i++)
      {
         ele_vector[i]->GetEdges (ele_edges_vector);
         size_t n_edges (ele_vector[i]->GetEdgesNumber());
         // loop over all edges of the i-th element
         for (size_t j=0; j<n_edges; j++)
         {
            ele_edges_vector[j]->GetNodes(edge_nodes);
            size_t selected (0);
            // get all elements having an edge in common with ply
            for (size_t k=0; k<nodes_near_ply.size(); k++)
            {
               if (edge_nodes[0]->GetIndex() == nodes_near_ply[k])
                  selected++;
               if (edge_nodes[1]->GetIndex() == nodes_near_ply[k])
                  selected++;
            }
            if (selected == 2)
            {
               ele_vector_ply.push_back(ele_vector[i]->GetIndex());
            }
         }
      }
   }
   /**************************************************************************
   MSHLib-Method:
   Task: All nodes vertical to a polyline
   02/2009 OK
   **************************************************************************/
   void CFEMesh::GetNODOnSFC_PLY_Z(Surface*m_sfc,std::vector<long>&msh_nod_vector)
   {
	   std::vector<CGLPolyline*>::iterator p_ply;
      CGLPolyline* m_ply = NULL;
      // .................................................................
      // nodes close to first polyline
      p_ply = m_sfc->polyline_of_surface_vector.begin();
      while(p_ply!=m_sfc->polyline_of_surface_vector.end())
      {
         m_ply = *p_ply;
         GetNODOnPLY_XY(m_ply,msh_nod_vector);
         break;
      }
   }
#endif                                         //WW#ifndef NON_GEO
   //-------------------------------------------------------------------------

   /**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   05/2005 OK Implementation: hex elements to line nodes
    last modification:
    **************************************************************************/
   // REMOVE CANDIDATE
   void CFEMesh::SetELE2NODTopology()
   {
      /* //WW TODO
       int k;
       long j;
       double xr[4],yr[4],zr[4];
       CGLPoint m_pnt;
       CGLPoint m_pnt1,m_pnt2,m_pnt3,m_pnt4;
       FiniteElement::CElement* m_ele = NULL;
       CFEMesh* m_msh_cond = NULL;
       CMSHNodes* m_mod = NULL;
       //----------------------------------------------------------------------
       m_msh_cond = FEMGet("RICHARDS_FLOW");
      for(long i=0;i<(long)ele_vector.size();i++){
      m_ele = ele_vector[i];
      for(k=0;k<4;k++){
      xr[k] = nod_vector[m_ele->nodes_index[k]]->x;
      yr[k] = nod_vector[m_ele->nodes_index[k]]->y;
      zr[k] = nod_vector[m_ele->nodes_index[k]]->z;
      }
      for(j=0;j<(long)m_msh_cond->nod_vector.size();j++){
      m_mod = m_msh_cond->nod_vector[j];
      m_mod->nodenumber = j;
      m_pnt.x = m_msh_cond->nod_vector[j]->x;
      m_pnt.y = m_msh_cond->nod_vector[j]->y;
      m_pnt.z = m_msh_cond->nod_vector[j]->z;
      if(m_pnt.PointInRectangle(xr,yr,zr)){
      ele_vector[i]->nod_vector.push_back(m_mod);
      //cout << i << " " << j << endl;
      }
      }
      }
      */
   }

   /**************************************************************************
    GeoSys-Method:
    Task:
    Programing:
    07/2005 OK Implementation
    **************************************************************************/
   std::ios::pos_type CFEMesh::GMSReadTIN(std::ifstream *tin_file)
   {
	   std::string line_string, s_buff;
      std::stringstream in;
      std::ios::pos_type position;
      int i, ibuf;
      double d_buf;
      i = 0;
      ibuf = 0;
      d_buf = 0.0;
      std::string line;
      std::string sub_line;
      long no_vertexes = 0;
      CNode* m_nod = NULL;
      long no_triangles = 0;
      CElem* m_ele = NULL;
      double xyz[3];
      //========================================================================
      while (sub_line.compare("ENDT"))
      {
         getline(*tin_file, s_buff);              //WW
         in.str(s_buff);
         in >> sub_line;
         //................................................................
                                                  // TNAM "PriTIN_1gr"
         if (sub_line.find("TNAM") != std::string::npos)
         {
            in >> sub_line;
         }
         //................................................................
                                                  // VERT 3173
         if (sub_line.find("VERT") != std::string::npos)
         {
            in >> no_vertexes;
            in.clear();
            for (i = 0; i < no_vertexes; i++)
            {
               m_nod = new CNode(i);
               getline(*tin_file, s_buff);        //WW
               in.str(s_buff);
               in >> xyz[0] >> xyz[1] >> xyz[2];
               m_nod->SetCoordinates(xyz);
               nod_vector.push_back(m_nod);
               in.clear();
            }
         }
         //................................................................
         if (sub_line.find("TRI") != std::string::npos)// TRI 6117
         {
            // Evaluate ele_type
            in >> no_triangles;
            in.clear();
            /* OKWW
             for(i=0;i<no_triangles;i++){
             m_ele = new FiniteElement::CElement();
             m_ele->type_name = "tri";
             m_ele->ElementType = ele_type;
             m_ele->nnodes = 3;
             m_ele->nodes = new long[3];
             in.str(GetLineFromFile1(tin_file)); // 3169	3168	3173
             in >> ele_nod_number;
             m_ele->nodes_index[0] = ele_nod_number-1;
             in >> ele_nod_number;
            m_ele->nodes_index[1] = ele_nod_number-1;
            in >> ele_nod_number;
            m_ele->nodes_index[2] = ele_nod_number-1;
            in.clear();
            m_ele->nodes_xyz = new double[9];
            for(k=0;k<m_ele->nnodes;k++){
            m_ele->nodes_xyz[k]                 = nod_vector[m_ele->nodes_index[k]]->x;
            m_ele->nodes_xyz[k+m_ele->nnodes]   = nod_vector[m_ele->nodes_index[k]]->y;
            m_ele->nodes_xyz[k+2*m_ele->nnodes] = nod_vector[m_ele->nodes_index[k]]->z;
            }
            ele_vector.push_back(m_ele);
            }
            */
            for (i = 0; i < no_triangles; i++)
            {
               m_ele = new CElem(i);
               m_ele->geo_type = MshElemType::TRIANGLE;
               m_ele->Read(*tin_file, 3);
               ele_vector.push_back(m_ele);
            }
         }
         //................................................................
         in.clear();
      }
      return position;
   }

   /**************************************************************************
    GeoSys-Method:
    Task:
    Programing:
    02/2005 OK Implementation (MMP groups)
    02/2005 OK Activate from vector
    08/2005 WW Changes due to geometry objects applied
    see also: BuildActiveElementsArray
    **************************************************************************/
   void CFEMesh::SetActiveElements(std::vector<long>&elements_active)
   {
      long i;
      //-----------------------------------------------------------------------
      for (i = 0; i < (long) this->ele_vector.size(); i++)
      {
         ele_vector[i]->MarkingAll(false);
      }
      //-----------------------------------------------------------------------
      for (i = 0; i < (long) elements_active.size(); i++)
      {
         ele_vector[elements_active[i]]->MarkingAll(true);
      }
      //-----------------------------------------------------------------------
      // Inactivate element with -1 MMP group
   }

   /**************************************************************************
    MSHLib-Method:
    Task:
    const int NLayers        :      Number of layers (start mesh)
    const int Layer          :      Layer number of the layer to be refined
    const int SUBLayers      :      Number of sublayers
    Programing:
    10/2003 WW/MB Erste Version (PrismRefine)
    08/2004 MB NElementsPerLayer = msh_no_pris / NLayers;
    09/2005 OK MSH
    ***************************************************************************/
   void CFEMesh::PrismRefine(const int Layer, const int subdivision)
   {
      const int nn = 6;
      int i, j, nes;
      int iii;
      long *element_nodes = NULL;
      //WW  static long ele_nodes[6];
      static double nx[6], ny[6], nz[6];
      //WW  static int newNode0[3], newNode1[3];
      double dx[3], dy[3], dz[3];
      double newz;
      //WW  static double  newx0[3],newy0[3],newz0[3];
      //WW  static double  newx1[3],newy1[3],newz1[3];
      int NumNodesNew, NumElementNew;
      // KR not needed long *knoten = NULL;
      int NNodesPerRow = 0;
      int NElementsPerLayer = 0;
      int row;
      int NRowsToShift;
      int NRows;
      int count;
      int CountNLayers;
      CNode* m_nod = NULL;
      double xyz[3];
      //----------------------------------------------------------------------
      const int NSubLayers = subdivision + 1;     //OK
      const long NumElement0 = (long) ele_vector.size();
                                                  //NodeListSize() / (NLayers+1);
      NNodesPerRow = (long) nod_vector.size() / (_n_msh_layer + 1);
                                                  //msh_no_pris / NLayers;
      NElementsPerLayer = (long) ele_vector.size() / _n_msh_layer;
      row = Layer;
      NRows = _n_msh_layer + 1;
      NumNodesNew = (long) nod_vector.size() - 1; //NodeListSize()-1;
                                                  //ElListSize()-1;
      NumElementNew = (long) ele_vector.size() - 1;
      //----------------------------------------------------------------------
      long nod_vector_size_add = NNodesPerRow * subdivision;
      for (i = 0; i < nod_vector_size_add; i++)
      {
         m_nod = new CNode(i);
         nod_vector.push_back(m_nod);
      }
      //nod_vector.resize(nod_vector_size_new);
      //----------------------------------------------------------------------
      // Initialisierung der Knoten flags
      for (i = 0; i < (long) nod_vector.size(); i++)
      {
         nod_vector[i]->SetMark(true);
      }
      //======================================================================
      CElem* m_ele = NULL;
      CElem* m_ele_new = NULL;
      for (int ne = 0; ne < NumElement0; ne++)
      {
         m_ele = ele_vector[ne];
         if (m_ele->GetElementType() == MshElemType::PRISM)
         {
            //element_nodes = m_ele->nodes;
            CountNLayers = _n_msh_layer;
            for (i = 0; i < nn; i++)
            {
               nx[i] = nod_vector[m_ele->nodes_index[i]]->X();
               ny[i] = nod_vector[m_ele->nodes_index[i]]->Y();
               nz[i] = nod_vector[m_ele->nodes_index[i]]->Z();
            }
            nes = 0;
            for (i = 0; i < 3; i++)
            {
               if (element_nodes[i] >= (row - 1) * NNodesPerRow
                  && element_nodes[i] <= (row * NNodesPerRow) - 1)
               {
                  nes++;
               }
            }
            if (nes == 3)
            {
               for (i = 0; i < 3; i++)
               {
                  dx[i] = (nx[i + 3] - nx[i]) / (float) NSubLayers;
                  dy[i] = (ny[i + 3] - ny[i]) / (float) NSubLayers;
                  dz[i] = (nz[i + 3] - nz[i]) / (float) NSubLayers;
               }
               //----------------------------------------------------------------
               // Create new nodes
                                                  // Loop over SubLayers
               for (iii = 0; iii < NSubLayers - 1; iii++)
               {
                  //..............................................................
                  // neue Knoten ganz unten
                  for (i = 0; i < 3; i++)
                  {
                     //if(NODGetFreeSurfaceFlag(element_nodes[i])==0){
                     if (nod_vector[element_nodes[i]]->GetMark())
                     {
                        //m_nod = new CMSHNodes(); //kno = (Knoten *)CreateNodeGeometry();
                        m_nod
                           = nod_vector[(m_ele->nodes_index[i]
                           + ((CountNLayers + 2) - row)
                           * NNodesPerRow)];
                        xyz[0]
                           = nod_vector[m_ele->nodes_index[i]
                           + ((CountNLayers + 1) - row)
                           * NNodesPerRow]->X();
                        xyz[1]
                           = nod_vector[m_ele->nodes_index[i]
                           + ((CountNLayers + 1) - row)
                           * NNodesPerRow]->Y();
                        xyz[2]
                           = nod_vector[m_ele->nodes_index[i]
                           + ((CountNLayers + 1) - row)
                           * NNodesPerRow]->Z();
                        //PlaceNode(kno,(element_nodes[i] + ((CountNLayers+2) - row) * NNodesPerRow));
                        m_nod->SetCoordinates(xyz);
                        nod_vector[(m_ele->nodes_index[i] + ((CountNLayers
                           + 2) - row) * NNodesPerRow)] = m_nod;
                     }
                  }
                  //..............................................................
                  // neues Element ganz unten
                  m_ele_new = new CElem();
                  //m_ele_new = m_ele;
                  m_ele_new->nnodes = m_ele->nnodes;
                  m_ele_new->patch_index = m_ele->GetPatchIndex();
                  m_ele_new->SetElementType(m_ele->GetElementType());
                  // KR memory leak! knoten = new long[6];
                  for (j = 0; j < 3; j++)
                  {
                     m_ele_new->nodes_index[j] = element_nodes[j] + ((CountNLayers + 1) - row) * NNodesPerRow;
                     m_ele_new->nodes_index[j+3] = element_nodes[j] + ((CountNLayers + 2) - row) * NNodesPerRow;
                  }
                  // KR fixed memory leak, knoten not needed
                  //for (j = 0; j < 6; j++) {
                  //	m_ele_new->nodes_index[j] = knoten[j];
                  //}
                  ele_vector.push_back(m_ele_new);
                  //..............................................................
                  /* "rowx hochziehen"   */
                  /* loop ?er die betroffenen rows   */
                  NRowsToShift = NRows - Layer;
                  count = 0;
                  for (i = NRowsToShift; i > 0; i--)
                  {
                     if (i != 1)
                     {
                        count++;
                        for (j = 0; j < 3; j++)
                        {
                           //if(NODGetFreeSurfaceFlag(element_nodes[j])==0){
                           if (nod_vector[element_nodes[j]]->GetMark())
                           {
                              m_nod = nod_vector[element_nodes[j + 3]
                                 + NNodesPerRow * (iii + i - 1)];
                              m_nod->SetZ(
                                 nod_vector[element_nodes[j]
                                 + NNodesPerRow * (iii + i
                                 - 1)]->Z());
                              //SetNodeZ(element_nodes[j+3] + NNodesPerRow*(iii+i-1), GetNodeZ(element_nodes[j] + NNodesPerRow*(iii+i-1)));
                           }
                        }
                     }
                     else
                     {
                        for (j = 0; j < 3; j++)
                        {
                           //if(NODGetFreeSurfaceFlag(element_nodes[j])==0)   {
                           if (nod_vector[element_nodes[j]]->GetMark())
                           {
                              newz = nod_vector[element_nodes[j]]->Z()
                                 + (dz[j] * (iii + 1));
                              //newz = GetNodeZ(element_nodes[j]) + (dz[j]*(iii+1));
                              nod_vector[element_nodes[j] + (i)
                                 * NNodesPerRow * (iii + 1)]->SetZ(
                                 newz);
                              //SetNodeZ(element_nodes[j] + (i) * NNodesPerRow *(iii+1), newz);
                           }
                        }
                     }
                  }                               /* end for Rows to shift */
                  //..............................................................
                  if (iii == NSubLayers - 2)
                  {
                     for (j = 0; j < 3; j++)
                     {
                        //NODSetFreeSurfaceFlag(element_nodes[j], 33);
                        nod_vector[element_nodes[j]]->SetMark(false);
                     }
                  }
                  CountNLayers++;
               }                                  /* End Loop over SubLayers  */
            }                                     /* End if nes==3 */
         }                                        /* Elementtyp ==6 */
      }                                           /* Element loop */
      //======================================================================
      _n_msh_layer += subdivision;
      //----------------------------------------------------------------------
   }

   /**************************************************************************
    MSHLib-Method:
    Programing:
    09/2005 TK/OK Implementation ab 4.2
    ***************************************************************************/
   //void CFEMesh::EdgeLengthMinMax()
   //{
   //  int j;
   //  long i;
   //  double edge_length;
   //  CNode* m_nod1 = NULL;
   //  CNode* m_nod2 = NULL;
   //  CElem* m_ele = NULL;
   //  for(i=0;i<(long)ele_vector.size();i++){
   //    m_ele = ele_vector[i];
   //    for(j=0;j<m_ele->nnodes-1;j++){
   //      m_nod1 = nod_vector[m_ele->nodes_index[j]];
   //      m_nod2 = nod_vector[m_ele->nodes_index[j+1]];
   //      edge_length = sqrt(((m_nod1->X()-m_nod2->X())*(m_nod1->X()-m_nod2->X()))+
   //                         ((m_nod1->Y()-m_nod2->Y())*(m_nod1->Y()-m_nod2->Y()))+
   //                         ((m_nod1->Z()-m_nod2->Z())*(m_nod1->Z()-m_nod2->Z())));
   //      if(i==0&&j==0)
   //      {
   //        _min_edge_length = edge_length;
   //        max_edge_length = edge_length;
   //      }
   //      else
   //      {
   //        if(_min_edge_length > edge_length)
   //          _min_edge_length = edge_length;
   //        if (max_edge_length < edge_length)
   //          max_edge_length = edge_length;
   //      }
   //    }
   //    m_nod1 = nod_vector[m_ele->nodes_index[m_ele->nnodes-1]];
   //    m_nod2 = nod_vector[m_ele->nodes_index[0]];
   //    edge_length = sqrt(((m_nod1->X()-m_nod2->X())*(m_nod1->X()-m_nod2->X()))+
   //                       ((m_nod1->Y()-m_nod2->Y())*(m_nod1->Y()-m_nod2->Y()))+
   //                       ((m_nod1->Z()-m_nod2->Z())*(m_nod1->Z()-m_nod2->Z())));
   //    if(_min_edge_length > edge_length)
   //      _min_edge_length = edge_length;
   //    if (max_edge_length < edge_length)
   //        max_edge_length = edge_length;
   //  }
   //}

   // TF the following two methods are not used, at least in the standard config
   /**************************************************************************
    MSHLib-Method:
    Task:
    Programing:
    09/2005 OK Implementation
    **************************************************************************/
   // TF m_sfc->PointInSurface(&m_point) returns always false
   // REMOVE CANDIDATE
   //void CFEMesh::SetMATGroupFromVOLLayer(CGLVolume*m_vol) {
   //	CElem* m_ele = NULL;
   //	//......................................................................
   ////	CGLPoint m_pnt;
   //	Surface* m_sfc = NULL;
   //	vector<Surface*>::const_iterator p_sfc;
   //	p_sfc = m_vol->surface_vector.begin();
   //	m_sfc = *p_sfc;
   //	//......................................................................
   //
   ////	long ep_layer = (long) ele_vector.size() / _n_msh_layer;
   ////	long jb = (m_vol->layer - 1) * ep_layer;
   ////	long je = jb + ep_layer;
   ////	double* xy;
   ////	for (long j = jb; j < je; j++) {
   ////		// Point in surface
   ////		m_ele = ele_vector[j];
   ////		xy = m_ele->GetGravityCenter();
   ////		m_pnt.x = xy[0];
   ////		m_pnt.y = xy[1];
   ////		if (m_sfc->PointInSurface(&m_pnt))
   ////			m_ele->SetPatchIndex(m_vol->mat_group);
   ////	}
   //}

   /**************************************************************************
    MSHLib-Method:
    Task:
    Programing:
    09/2005 OK Implementation
    **************************************************************************/
   // REMOVE CANDIDATE
   // TF uses method SetMATGroupFromVOLLayer, which can be removed!?
   //void CFEMesh::SetMATGroupsFromVOLLayer() {
   //	CGLVolume* m_vol;
   //	vector<CGLVolume*>::const_iterator p_vol = volume_vector.begin();
   //	while (p_vol != volume_vector.end()) {
   //		m_vol = *p_vol;
   //		//....................................................................
   //		// LAYER
   //		SetMATGroupFromVOLLayer(m_vol);
   //		//....................................................................
   //		++p_vol;
   //	}
   //}

   /**************************************************************************
    FEMLib-Method:
    Task:
    Programing:
    11/2005 MB Implementation based on NodeExists.
    03/2010 TF changed data type of loop variables
    **************************************************************************/
   bool CFEMesh::NodeExists(size_t node)
   {
      size_t no_nodes (nod_vector.size());

      for (size_t i=0; i<no_nodes; i++)
      {
         if (node == (size_t)Eqs2Global_NodeIndex[i])
            return true;
      }
      return false;
   }

   /**************************************************************************
    MSHLib-Method:
    Task:
    Programing:
    10/2005 OK Implementation
    02/2006 WW Ordering and remove bugs
    **************************************************************************/
   void CFEMesh::ConnectedNodes(bool quadratic) const
   {
      bool exist = false;
#define noTestConnectedNodes
      //----------------------------------------------------------------------
      for (size_t i = 0; i < nod_vector.size(); i++)
      {
         CNode * nod = nod_vector[i];
         size_t n_connected_elements (nod->connected_elements.size());
         for (size_t j = 0; j < n_connected_elements; j++)
         {
            CElem *ele = ele_vector[nod->connected_elements[j]];
            size_t n_quadratic_node (static_cast<size_t>(ele->GetNodesNumber(quadratic)));
            for (size_t l = 0; l < n_quadratic_node; l++)
            {
               exist = false;
               size_t n_connected_nodes (nod->connected_nodes.size());
                                                  //WW
               for (size_t k = 0; k < n_connected_nodes; k++)
               {
                  if (nod->connected_nodes[k] == ele->nodes_index[l])
                  {
                     exist = true;
                     break;
                  }
               }
               if (!exist)                        //WW
                  nod->connected_nodes.push_back(ele->nodes_index[l]);
            }
         }
      }

      // Sorting. WW
                                                  //WW
      for (size_t i = 0; i < nod_vector.size(); i++)
      {
         CNode *nod = nod_vector[i];
         size_t n_connected_nodes = nod->connected_nodes.size();
         for (size_t k = 0; k < n_connected_nodes; k++)
         {
            for (size_t l = k; l < n_connected_nodes; l++)
            {
               if (nod->connected_nodes[l] < nod->connected_nodes[k])
               {
                  long n = nod->connected_nodes[k];
                  nod->connected_nodes[k] = nod->connected_nodes[l];
                  nod->connected_nodes[l] = n;
               }
            }
         }

      }
      //----------------------------------------------------------------------
#ifdef TestConnectedNodes
      for(i=0;i<(long)nod_vector.size();i++)
      {
         nod = nod_vector[i];
         cout << (int)nod->connected_nodes.size() << ": ";
         for(m=0;m<(int)nod->connected_nodes.size();m++)
         {
            cout << nod->connected_nodes[m] << " ";
         }
         cout << endl;
      }
#endif
   }

   /**************************************************************************
    GeoLib-Method:
    Task:
    Programing:
    03/2004 OK Implementation
    11/2005 OK MSH
    03/2006 CC
    **************************************************************************/
   // 07/2010 TF commented out since method Surface::IsPointInSurface returns always false
   //void CFEMesh::CreateSurfaceTINfromTri(Surface*m_sfc) {
   //	int j;
   //	CGLPoint m_point;
   //	CTriangle *m_triangle;
   //	CElem* m_ele = NULL;
   //	double* xyz;
   //	vec<long> node_indeces(3);
   //	m_sfc->TIN = new CTIN;//CC
   //	//----------------------------------------------------------------------
   //	for (long i = 0; i < (long) ele_vector.size(); i++) {
   //		m_ele = ele_vector[i];
   //		if (m_ele->GetElementType() == 4) { // use only triangles
   //			xyz = m_ele->GetGravityCenter();
   //			m_point.x = xyz[0];
   //			m_point.y = xyz[1];
   //			m_point.z = xyz[2];
   //			if (IsPointInSurface(m_sfc, &m_point)) {
   //				m_triangle = new CTriangle;
   //				m_triangle->number = (long) m_sfc->TIN->Triangles.size();
   //				m_ele->GetNodeIndeces(node_indeces);
   //				for (j = 0; j < 3; j++) {
   //					m_triangle->x[j] = nod_vector[node_indeces[j]]->X();
   //					m_triangle->y[j] = nod_vector[node_indeces[j]]->Y();
   //					m_triangle->z[j] = nod_vector[node_indeces[j]]->Z();
   //				}
   //				m_sfc->TIN->Triangles.push_back(m_triangle);
   //				m_sfc->TIN->name = m_sfc->name;//CC
   //			} // element found
   //		} // triangle
   //	} // ele_vector
   //}

   /**************************************************************************
    GeoLib-Method:
    Task: Geometric / topological method
    Programing:
    01/2005 OK Implementation
    01/2005 CC add variable :long i for the function
    08/2005 CC move from Geolib to Mshlib
    **************************************************************************/
   //void CFEMesh::CreateLayerSurfaceTINsfromPris(Surface*m_sfc) {
   //	int j;
   //	CGLPoint m_point;
   //	CTriangle *m_tri;
   //	CTriangle *m_tri_new;
   //	CElem* m_ele = NULL;
   //	double* xyz;
   //	vec<long> node_indeces(6);
   //	//---------------------------------------------------------------------
   //	if (!m_sfc) {
   //		return;
   //	}
   //	//---------------------------------------------------------------------
   //	// Create layer surfaces
   //	char layer_number[3];
   //	Surface *m_sfc_layer = NULL;
   //	CTIN *m_TIN = NULL;
   //	string sfc_layer_name = m_sfc->name + "_layer_";
   //	for (int l = 0; l < _n_msh_layer + 1; l++) {
   //		m_sfc_layer = new Surface;
   //		m_sfc_layer->type_name = "TIN";
   //		sprintf(layer_number, "%i", l);
   //		m_sfc_layer->name = sfc_layer_name + layer_number;
   //		m_sfc_layer->data_name = m_sfc_layer->name + ".tin";
   //		m_TIN = new CTIN;
   //		m_TIN->name = m_sfc_layer->name;
   //		m_sfc_layer->TIN = m_TIN;
   //		surface_vector.push_back(m_sfc_layer);
   //	}
   //	//----------------------------------------------------------------------
   //	Surface *m_sfc_layer_0 = NULL;
   //	sfc_layer_name = m_sfc->name + "_layer_0";
   //	m_sfc_layer_0 = GEOGetSFCByName(sfc_layer_name);
   //	if (!m_sfc_layer_0) {
   //		return;
   //	}
   ////	long no_nodes_per_layer = (long) nod_vector.size() / (_n_msh_layer + 1);
   ////	long no_elements_per_layer = (long) ele_vector.size() / (_n_msh_layer);
   ////	for (long i = 0; i < no_elements_per_layer; i++) {
   ////		m_ele = ele_vector[i];
   ////		if (m_ele->GetElementType() == 6) { // prism
   ////			xyz = m_ele->GetGravityCenter();
   ////			m_point.x = xyz[0];
   ////			m_point.y = xyz[1];
   ////			m_point.z = xyz[2];
   ////			if (IsPointInSurface(m_sfc, &m_point)) {
   ////				m_tri = new CTriangle;
   ////				m_tri->number = (long) m_sfc_layer_0->TIN->Triangles.size();
   ////				m_ele->GetNodeIndeces(node_indeces);
   ////				for (j = 0; j < 3; j++) {
   ////					m_tri->x[j] = nod_vector[node_indeces[j]]->X();
   ////					m_tri->y[j] = nod_vector[node_indeces[j]]->Y();
   ////					m_tri->z[j] = nod_vector[node_indeces[j]]->Z();
   ////				}
   ////				m_sfc_layer_0->TIN->Triangles.push_back(m_tri);
   ////				m_sfc_layer_0->TIN->name = m_sfc_layer_0->name;//CC
   ////				for (int l = 0; l < _n_msh_layer; l++) {
   ////					sprintf(layer_number, "%i", l + 1);
   ////					sfc_layer_name = m_sfc->name + "_layer_" + layer_number;
   ////					m_sfc_layer = GEOGetSFCByName(sfc_layer_name);
   ////					if (!m_sfc_layer) {
   ////						return;
   ////					}
   ////					m_tri_new = new CTriangle;
   ////					for (j = 0; j < 3; j++) {
   ////						m_tri_new->number = m_tri->number;
   ////						m_tri_new->x[j] = m_tri->x[j];
   ////						m_tri_new->y[j] = m_tri->y[j];
   ////						m_tri_new->z[j] = nod_vector[node_indeces[j + 3] + l
   ////								* no_nodes_per_layer]->Z();
   ////					}
   ////					m_sfc_layer->TIN->Triangles.push_back(m_tri_new);
   ////					m_sfc_layer->TIN->name = m_sfc_layer->name;//CC
   ////				}
   ////			} // element found
   ////		} // triangle
   ////	} // ele_vector
   //	//---------------------------------------------------------------------
   //}

   /**************************************************************************
    FEMLib-Method:
    Task:  Renumbering nodes corresponding to the activiate of elements
    Programing:
    01/2006 YD Implementation
    **************************************************************************/
   void CFEMesh::FaceNormal()
   {
      int i, j;
      int idx0_face, idx1_face, idx_owner, index0 = 0, index1 = 0;

      CElem* elem = NULL;
      CElem* elem_face = NULL;

      // if(coordinate_system!=32)
      //   return;
      if ((long) face_normal.size() > 0)
         return;                                  //WW
      //------------------------
      for (i = 0; i < (long) face_vector.size(); i++)
      {
         elem_face = face_vector[i];
         elem = face_vector[i]->GetOwner();

         if (elem->GetElementType() == MshElemType::LINE)
            return;

         //WW      int no_face_vertex = face_vector[i]->GetVertexNumber();
         int no_owner_vertex = face_vector[i]->GetOwner()->GetVertexNumber();
         idx0_face = face_vector[i]->GetNodeIndex(0);
         idx1_face = face_vector[i]->GetNodeIndex(1);
		 double* normal = new double[3];
         for (j = 0; j < no_owner_vertex; j++)
         {
            idx_owner = face_vector[i]->GetOwner()->GetNodeIndex(j);
            if (idx0_face == idx_owner)
               index0 = j;
         }
         for (j = 0; j < no_owner_vertex; j++)
         {
            idx_owner = face_vector[i]->GetOwner()->GetNodeIndex(j);
            if (idx1_face == idx_owner)
               index1 = j;
         }
         if (elem->GetMark())
         {
            if ((index1 - index0) >= 1)
               elem->FaceNormal(index0, index1, normal);
            else
               elem->FaceNormal(index1, index0, normal);

         }
         face_normal.push_back(normal);
         elem_face->ComputeVolume();
      }
   }

#ifndef NON_GEO
   /**************************************************************************
    MSHLib-Method:
    Programing:
    02/2006 OK Implementation
    **************************************************************************/
   void CFEMesh::CreateLineELEFromTri()
   {
      int k;
      long i;                                     //,e;
      //WW  double v1[3],v2[3],v3[3];
      //  double patch_area;
      double x, y, z;
      double x0, y0, z0;
      double x1, y1, z1;
      double dl;
      CNode* m_nod = NULL;
      //WW  CNode* m_nod1 = NULL;
      //WW  CNode* m_nod2 = NULL;
      CNode* m_nod_line = NULL;
      CElem* m_tri_ele = NULL;
      CElem* m_ele = NULL;
      //  CElem* m_ele1 = NULL;
      //----------------------------------------------------------------------
      // 1 - Element normal vector (for 2D elements only)
      FillTransformMatrix();
      //----------------------------------------------------------------------
      // 2 - Node patch area
      SetNODPatchAreas();
      //----------------------------------------------------------------------
      // 3 - Intersection nodes
      //OKSetNetworkIntersectionNodes();
      //----------------------------------------------------------------------
      // 4 - Create MSH
      MSHDelete("LINE_from_TRI");
      CFEMesh* m_msh_line (new CFEMesh(_geo_obj,_geo_name));
      m_msh_line->pcs_name = "LINE_from_TRI";
      m_msh_line->setElementType (MshElemType::LINE);
      m_msh_line->_n_msh_layer = 10;              // User-defined
      double element_length = -0.1;               // User-defined
      dl = element_length * m_msh_line->_n_msh_layer;
      //----------------------------------------------------------------------
      // 4.1 - Line nodes
      for (i = 0; i < (long) nod_vector.size(); i++)
      {
         m_nod = nod_vector[i];
         //OKif(m_nod->selected)
         //OKcontinue;
         if ((int) m_nod->connected_elements.size() == 0)
            continue;
         m_tri_ele = ele_vector[m_nod->connected_elements[0]];
         //....................................................................
         // Node normal vector
         x0 = m_nod->X();
         y0 = m_nod->Y();
         z0 = m_nod->Z();
         //    x1 = x0 + m_tri_ele->normal_vector[0]*dl;
         //    y1 = y0 + m_tri_ele->normal_vector[1]*dl;
         //    z1 = z0 + m_tri_ele->normal_vector[2]*dl;
                                                  //WW
         x1 = x0 + (*m_tri_ele->tranform_tensor)(2, 0) * dl;
                                                  //WW
         y1 = y0 + (*m_tri_ele->tranform_tensor)(2, 1) * dl;
                                                  //WW
         z1 = z0 + (*m_tri_ele->tranform_tensor)(2, 2) * dl;
         //....................................................................
         for (size_t j = 0; j < m_msh_line->_n_msh_layer + 1; j++)
         {
            x = x0 + (x1 - x0) * (j) / m_msh_line->_n_msh_layer;
            y = y0 + (y1 - y0) * (j) / m_msh_line->_n_msh_layer;
            z = z0 + (z1 - z0) * (j) / m_msh_line->_n_msh_layer;
            m_nod_line = new CNode((long) m_msh_line->nod_vector.size(), x, y,
               z);
            m_msh_line->nod_vector.push_back(m_nod_line);
         }
      }
      //----------------------------------------------------------------------
      // 4.2 - Line elements
      long i_count = 0;
      for (i = 0; i < (long) nod_vector.size(); i++)
      {
         m_nod = nod_vector[i];
         //....................................................................
         // Intersection node
         if (m_nod->selected)
            continue;
         if ((int) m_nod->connected_elements.size() == 0)
            continue;
         m_tri_ele = ele_vector[m_nod->connected_elements[0]];
         //....................................................................
         // Line elements
         for (size_t j = 0; j < m_msh_line->_n_msh_layer; j++)
         {
            m_ele = new Mesh_Group::CElem;
            m_ele->SetIndex((long) m_msh_line->ele_vector.size());
            m_ele->SetElementType(MshElemType::LINE);
            m_ele->nnodes = 2;
            m_ele->nodes_index.resize(m_ele->nnodes);
            //....................................................................
            // Line element nodes
            for (k = 0; k < m_ele->nnodes; k++)
            {
               m_ele->nodes_index[k] = i_count * m_msh_line->_n_msh_layer + j
                  + k + i_count;
               m_ele->nodes[k] = m_msh_line->nod_vector[m_ele->nodes_index[k]];
            }
            //....................................................................
            m_msh_line->ele_vector.push_back(m_ele);
         }
         i_count++;
      }
      //----------------------------------------------------------------------
      if (m_msh_line->ele_vector.size() > 0)
         fem_msh_vector.push_back(m_msh_line);
      else
         delete m_msh_line;
      //----------------------------------------------------------------------
      CGSProject* m_gsp = NULL;
      m_gsp = GSPGetMember("gli");
      if (m_gsp)
         MSHWrite(m_gsp->path + "test");
      else
         MSHWrite("test");
   }
   /**************************************************************************
    MSHLib-Method:
    Programing:
    04/2006 OK Implementation
    **************************************************************************/
   void CFEMesh::CreateLineELEFromTriELE()
   {
      int k;
      long i;
      double x, y, z;
      double x0, y0, z0;
      double x1, y1, z1;
      double dl;
      double* gravity_center;
      CNode* m_nod_line = NULL;
      CElem* m_tri_ele = NULL;
      CElem* m_ele = NULL;
      //----------------------------------------------------------------------
      // 1 - Element normal vector (for 2D elements only)
      FillTransformMatrix();
      //----------------------------------------------------------------------
      // 2 - Create MSH
      MSHDelete("LINE_from_TRI");
      CFEMesh* m_msh_line (new CFEMesh(_geo_obj,_geo_name));
      m_msh_line->pcs_name = "LINE_from_TRI";
      m_msh_line->setElementType  (MshElemType::LINE);
      m_msh_line->_n_msh_layer = 20;              // User-defined
      double element_length = -0.05;              // User-defined
      dl = element_length * m_msh_line->_n_msh_layer;
      //----------------------------------------------------------------------
      // 3.1 - Line nodes
      for (i = 0; i < (long) ele_vector.size(); i++)
      {
         m_tri_ele = ele_vector[i];
         //....................................................................
         // Element normal vector
         gravity_center = m_tri_ele->GetGravityCenter();
         x0 = gravity_center[0];
         y0 = gravity_center[1];
         z0 = gravity_center[2];
         //    x1 = x0 + m_tri_ele->normal_vector[0]*dl;
         //    y1 = y0 + m_tri_ele->normal_vector[1]*dl;
         //    z1 = z0 + m_tri_ele->normal_vector[2]*dl;

                                                  //WW
         x1 = x0 + (*m_tri_ele->tranform_tensor)(2, 0) * dl;
                                                  //WW
         y1 = y0 + (*m_tri_ele->tranform_tensor)(2, 1) * dl;
                                                  //WW
         z1 = z0 + (*m_tri_ele->tranform_tensor)(2, 2) * dl;

         //....................................................................
         for (size_t j = 0; j < m_msh_line->_n_msh_layer + 1; j++)
         {
            x = x0 + (x1 - x0) * (j) / m_msh_line->_n_msh_layer;
            y = y0 + (y1 - y0) * (j) / m_msh_line->_n_msh_layer;
            z = z0 + (z1 - z0) * (j) / m_msh_line->_n_msh_layer;
            m_nod_line = new CNode((long) m_msh_line->nod_vector.size(), x, y,
               z);
            m_msh_line->nod_vector.push_back(m_nod_line);
         }
      }
      //----------------------------------------------------------------------
      // 3.2 - Line elements
      long i_count = 0;
      for (i = 0; i < (long) ele_vector.size(); i++)
      {
         m_tri_ele = ele_vector[i];
         //....................................................................
         // Line elements
         for (size_t j = 0; j < m_msh_line->_n_msh_layer; j++)
         {
            m_ele = new Mesh_Group::CElem;
            m_ele->SetIndex((long) m_msh_line->ele_vector.size());
            m_ele->SetElementType(MshElemType::LINE);
            m_ele->nnodes = 2;
                                                  //OK4310
            m_ele->SetPatchIndex((int) mmp_vector.size());
            m_ele->nodes_index.resize(m_ele->nnodes);
            //....................................................................
            // Line element nodes
            for (k = 0; k < m_ele->nnodes; k++)
            {
               m_ele->nodes_index[k] = i_count * m_msh_line->_n_msh_layer + j
                  + k + i_count;
               m_ele->nodes[k] = m_msh_line->nod_vector[m_ele->nodes_index[k]];
            }
            //....................................................................
            m_msh_line->ele_vector.push_back(m_ele);
         }
         i_count++;
      }
      //----------------------------------------------------------------------
      if (m_msh_line->ele_vector.size() > 0)
         fem_msh_vector.push_back(m_msh_line);
      else
         delete m_msh_line;
      //----------------------------------------------------------------------
      CGSProject* m_gsp = NULL;
      m_gsp = GSPGetMember("gli");
      if (m_gsp)
         MSHWrite(m_gsp->path + "test");
      else
         MSHWrite("test");
   }

   /**************************************************************************
    MSHLib-Method:
    08/2006 OK Implementation
    **************************************************************************/
   void CFEMesh::GetELEOnPLY(CGLPolyline* m_ply, std::vector<long>&ele_vector_ply) const
   {
#ifdef MSH_CHECK
      cout << "CFEMesh::GetELEOnPLY" << endl;
#endif
      long i;
      int j;
      int k;
      CElem* m_ele = NULL;
      CEdge* m_edg = NULL;
      //WW  CNode* m_nod = NULL;
      vec<CEdge*> ele_edges_vector(15);
      std::vector<long> nodes_vector_ply;
      vec<long> ele_nodes(8);
      //int edge_node_numbers[2];
      vec<CNode*> edge_nodes(3);
      long edge_node_0, edge_node_1;
      long nn;
      //----------------------------------------------------------------------
      GetNODOnPLY(m_ply, nodes_vector_ply);
      //----------------------------------------------------------------------
#ifdef MSH_CHECK
      cout << "Elements at polyline: " << endl;
#endif
      switch (m_ply->getType())
      {
         //....................................................................
         case 0:                                  // PNT-TOP
            //....................................................................
         case 2:                                  // PNT-TOP CC!!!
            //..................................................................
            // All elements having 2 points in common with m_ply
            /*
             for(i=0;i<(long)ele_vector.size();i++)
             {
             m_ele = ele_vector[i];
             m_ele->selected = 0;
             m_ele->GetNodeIndeces(ele_nodes);
             for(j=0;j<(int)m_ele->GetNodesNumber(false);j++)
             {
             for(k=0;k<(long)nodes_vector_ply.size();k++)
             {
             if(ele_nodes[j]==nodes_vector_ply[k])
            {
            m_ele->selected++;
            }
            }
            }
            if(m_ele->selected==2)
            m_ele->SetMark(true);
            }
            */
            //..................................................................
            // All elements having an edge in common with m_ply
            for (i = 0; i < (long) ele_vector.size(); i++)
            {
               m_ele = ele_vector[i];
               m_ele->SetMark(false);
               m_ele->selected = 0;
               m_ele->GetEdges(ele_edges_vector);
               for (j = 0; j < (int) m_ele->GetEdgesNumber(); j++)
               {
                  m_edg = ele_edges_vector[j];
                  m_edg->SetMark(false);
               }
            }
            for (i = 0; i < (long) ele_vector.size(); i++)
            {
               m_ele = ele_vector[i];
               m_ele->SetMark(false);
               m_ele->GetEdges(ele_edges_vector);
               for (j = 0; j < (int) m_ele->GetEdgesNumber(); j++)
               {
                  m_edg = ele_edges_vector[j];
                  //m_ele->GetLocalIndicesOfEdgeNodes(j,edge_node_numbers);
                  m_edg->GetNodes(edge_nodes);
                  m_ele->selected = 0;
                  for (k = 0; k < (long) nodes_vector_ply.size(); k++)
                  {
                     nn = nodes_vector_ply[k];
                     //if(edge_node_numbers[0]==nodes_vector_ply[k])
                     edge_node_0 = edge_nodes[0]->GetIndex();
                     if (edge_nodes[0]->GetIndex() == (size_t)nodes_vector_ply[k])
                        m_ele->selected++;
                     //if(edge_node_numbers[1]==nodes_vector_ply[k])
                     edge_node_1 = edge_nodes[1]->GetIndex();
                     if (edge_nodes[1]->GetIndex() == (size_t)nodes_vector_ply[k])
                        m_ele->selected++;
                  }
                  if (m_ele->selected == 2)
                  {
                     m_ele->SetMark(true);
                     m_edg->SetMark(true);
                  }
               }
            }
            break;
            //....................................................................
         case 1:                                  // PLY-RAS
            break;
         default:
            std::cout << "Warning in CFEMesh::GetELEOnPLY: case not implemented" << std::endl;
      }
      //----------------------------------------------------------------------
      ele_vector_ply.clear();
      vec<long> node_indeces(20);
      for (i = 0; i < (long) ele_vector.size(); i++)
      {
         m_ele = ele_vector[i];
         m_ele->GetEdges(ele_edges_vector);
         if (m_ele->GetMark())
         {
            ele_vector_ply.push_back(m_ele->GetIndex());
#ifdef MSH_CHECK
            cout << "Element: " << m_ele->GetIndex() << ", Nodes:";
#endif
            m_ele->GetNodeIndeces(node_indeces);
#ifdef MSH_CHECK
            for(j=0;j<(int)m_ele->GetNodesNumber(false);j++)
            {
               cout << " " << node_indeces[j];
            }
#endif
         }
         for (j = 0; j < (int) m_ele->GetEdgesNumber(); j++)
         {
            m_edg = ele_edges_vector[j];
            if (m_edg->GetMark())
            {
               m_edg->GetNodes(edge_nodes);
#ifdef MSH_CHECK
               cout << ", Edge nodes: " << edge_nodes[0]->GetIndex() << "," << edge_nodes[1]->GetIndex() << endl;
#endif
            }
         }
      }
      nodes_vector_ply.clear();
   }

   /**************************************************************************
    MSHLib-Method:
    Programing:
    05/2006 OK Implementation
    08/2006 YD
    **************************************************************************/
   void CFEMesh::CreateLineELEFromSFC()
   {
      int k, i_k;
      double x0, y0, z0;
      double z, dz;
      Surface* m_sfc = NULL;
      CNode* m_nod = NULL;
      CElem* m_ele = NULL;
      CColumn* m_col = NULL;
      //  CGLLine* m_lin = NULL;
      CSoilProfile* m_prf = NULL;                 //YD
      //======================================================================
      i_k = 0;
      dz = -0.05;
      long i_count = 0;
      for (long i = 0; i < (long) surface_vector.size(); i++)
      {
         m_sfc = surface_vector[i];
         m_sfc->CalcCenterPoint();
         m_col = COLGet(m_sfc->name);
         if (!m_col)
            return;
         m_prf = profile_vector[m_sfc->profile_code - 1];
         //--------------------------------------------------------------------
         // NOD
         x0 = m_sfc->center_point[0];
         y0 = m_sfc->center_point[1];
         z0 = m_sfc->center_point[2];
         for (size_t j = 0; j < _n_msh_layer + 1; j++)
         {
            z = z0 + dz * (j);
            m_nod = new CNode((long) nod_vector.size(), x0, y0, z);
            nod_vector.push_back(m_nod);
         }
         //--------------------------------------------------------------------
         // ELE
         for (size_t j = 0; j < _n_msh_layer; j++)
         {
            //..................................................................
            m_ele = new Mesh_Group::CElem;
            m_ele->SetIndex((long) ele_vector.size());
            m_ele->SetElementType(MshElemType::LINE);
            m_ele->nnodes = 2;
            m_ele->nodes_index.resize(m_ele->nnodes);
            m_ele->SetPatchIndex(-1);
            m_ele->gravity_center[0] = 0.0;
            m_ele->gravity_center[1] = 0.0;
            m_ele->gravity_center[2] = 0.0;
            //..................................................................
            // Line element nodes
            for (k = 0; k < m_ele->nnodes; k++)
            {
               // if(k == 0) i_k=1;           //YD: Right habd rule
               // if(k == 1) i_k=0;
               m_ele->nodes_index[k] = i_count * _n_msh_layer + j + i_k
                  + i_count;
               m_ele->nodes[k] = nod_vector[m_ele->nodes_index[k]];
               m_ele->gravity_center[0]
                  += nod_vector[m_ele->nodes_index[k]]->X()
                  / m_ele->nnodes;
               m_ele->gravity_center[1]
                  += nod_vector[m_ele->nodes_index[k]]->Y()
                  / m_ele->nnodes;
               m_ele->gravity_center[2]
                  += nod_vector[m_ele->nodes_index[k]]->Z()
                  / m_ele->nnodes;
            }
            //..................................................................
            // MAT
            /*
             for(k=0;k<(int)m_col->line_vector.size();k++)
             {
             m_lin = m_col->line_vector[k];
             if((abs(m_ele->gravity_center[2])>m_lin->m_point1->z)&&(abs(m_ele->gravity_center[2])<m_lin->m_point2->z))
             {
             m_ele->SetPatchIndex(m_lin->mat_group);
             }
             }
             */

            for (k = 0; k < (int) m_prf->soil_layer_thickness.size() - 1; k++)
            {
               if ((fabs(m_ele->gravity_center[2])
                  > m_prf->soil_layer_thickness[k]) && (fabs(
                  m_ele->gravity_center[2])
                  < m_prf->soil_layer_thickness[k + 1]))
                  m_ele->SetPatchIndex(m_prf->soil_type[k]);
            }
            //..................................................................
            ele_vector.push_back(m_ele);
            //..................................................................
         }
         i_count++;
         //--------------------------------------------------------------------
      }
   }
#endif                                         //#ifndef NON_GEO

   /**************************************************************************
    MSHLib-Method:
    Programing:
    03/2006 OK Implementation
    **************************************************************************/
   void CFEMesh::SetNODPatchAreas()
   {
      long i, e;
      int j, k;
      int n1 = 0, n2 = 0;
      double v1[3], v2[3], v3[3];
      double patch_area;
      double x0, y0, z0;
      double* gravity_center;
      CNode* m_nod = NULL;
      CNode* m_nod1 = NULL;
      CNode* m_nod2 = NULL;
      CElem* m_ele = NULL;
      //----------------------------------------------------------------------
      for (i = 0; i < (long) nod_vector.size(); i++)
      {
         m_nod = nod_vector[i];                   // this node
         patch_area = 0.0;
         //....................................................................
         // triangle neighbor nodes
         for (j = 0; j < (int) m_nod->connected_elements.size(); j++)
         {
            e = m_nod->connected_elements[j];
            m_ele = ele_vector[e];
            for (k = 0; k < 3; k++)
            {
               if (m_ele->GetNodeIndex(k) == i)
               {
                  switch (k)
                  {
                     case 0:
                        n1 = 2;
                        n2 = 1;
                        break;
                     case 1:
                        n1 = 0;
                        n2 = 2;
                        break;
                     case 2:
                        n1 = 1;
                        n2 = 0;
                        break;
                  }
               }
            }
            //..................................................................
            gravity_center = m_ele->GetGravityCenter();
            v2[0] = gravity_center[0] - m_nod->X();
            v2[1] = gravity_center[1] - m_nod->Y();
            v2[2] = gravity_center[2] - m_nod->Z();
            //..................................................................
            m_nod1 = nod_vector[m_ele->GetNodeIndex(n1)];
            x0 = 0.5 * (m_nod1->X() - m_nod->X());
            y0 = 0.5 * (m_nod1->Y() - m_nod->Y());
            z0 = 0.5 * (m_nod1->Z() - m_nod->Z());
            v1[0] = x0 - m_nod->X();
            v1[1] = y0 - m_nod->Y();
            v1[2] = z0 - m_nod->Z();
            CrossProduction(v1, v2, v3);
            patch_area += 0.5 * MBtrgVec(v3, 3);
            //..................................................................
            m_nod2 = nod_vector[m_ele->GetNodeIndex(n2)];
            x0 = 0.5 * (m_nod2->X() - m_nod->X());
            y0 = 0.5 * (m_nod2->Y() - m_nod->Y());
            z0 = 0.5 * (m_nod2->Z() - m_nod->Z());
            v1[0] = x0 - m_nod->X();
            v1[1] = y0 - m_nod->Y();
            v1[2] = z0 - m_nod->Z();
            CrossProduction(v1, v2, v3);
            patch_area += 0.5 * MBtrgVec(v3, 3);
            //..................................................................
         }
         m_nod->patch_area = patch_area;
      }
   }

   /**************************************************************************
    MSHLib-Method:
    Programing:
    03/2006 OK Implementation
    **************************************************************************/
   void CFEMesh::SetNetworkIntersectionNodes()
   {
      long i, e;
      int j, k;
      double v3[3], nr1[3], nr2[3];
      //WW  double* gravity_center;
      CNode* m_nod = NULL;
      CElem* m_ele = NULL;
      CElem* m_ele1 = NULL;
      //----------------------------------------------------------------------
      // Is node intersection node
      for (i = 0; i < (long) nod_vector.size(); i++)
      {
         m_nod = nod_vector[i];
         m_nod->selected = false;
      }
      double eps = 1e-3;
      for (i = 0; i < (long) nod_vector.size(); i++)
      {
         m_nod = nod_vector[i];                   // this node
         if ((int) m_nod->connected_elements.size() == 0)
            continue;
         m_ele = ele_vector[m_nod->connected_elements[0]];
         for (k = 0; k < 3; k++)
            nr1[k] = (*m_ele->tranform_tensor)(2, k);
         //....................................................................
         // Compare element normal vectors
         for (j = 1; j < (int) m_nod->connected_elements.size(); j++)
         {
            e = m_nod->connected_elements[j];
            m_ele1 = ele_vector[e];
            for (k = 0; k < 3; k++)
               nr2[k] = (*m_ele1->tranform_tensor)(2, k);
            CrossProduction(nr1, nr2, v3);
            if (MBtrgVec(v3, 3) > eps)
               m_nod->selected = true;
         }
      }
      // Count non-intersection nodes
      long no_non_intersection_nodes = 0;
      for (i = 0; i < (long) nod_vector.size(); i++)
      {
         m_nod = nod_vector[i];
         if (m_nod->selected)
            continue;
         no_non_intersection_nodes++;
      }
   }

#ifndef NON_PROCESS                            // 05.03.2010 WW
#ifdef NEW_EQS                                 // 1.11.2007 WW
   /**************************************************************************
    MSHLib-Method:
    Programing:
    11/2007 WW Implementation
    **************************************************************************/
   void CFEMesh::CreateSparseTable()
   {
      // Symmetry case is skipped.
      // 1. Sparse_graph_H for high order interpolation. Up to now, deformation
      if(NodesNumber_Linear!=NodesNumber_Quadratic)
         sparse_graph_H = new SparseTable(this, true);
      // 2. M coupled with other processes with linear element
      if(sparse_graph_H)
      {
         if((int)pcs_vector.size()>1)
            sparse_graph = new SparseTable(this, false);
      }
      // 3. For process with linear elements
      else
         sparse_graph = new SparseTable(this, false);

      //sparse_graph->Write();
      //  sparse_graph_H->Write();
      //
      //ofstream Dum("sparse.txt", ios::out);
      //sparse_graph_H->Write(Dum);
   }
#endif                                         //#ifndef NON_PROCESS  // 05.03.2010 WW
#endif

	//---------------------------------------------------------------------------
	/*!
	\brief Import the MODFlow grid into OGS
	The finite difference grid for MODFlow is read, converted into
	hexahedra, and the the coverted elements are written in an OGS syntax mesh file.

	\param fname The file name.

	10/2009 WW
	//    By WW
	* 01/2010 TF changed signature to save a copy
	* 				improved implementation a little bit (still a lot to do)
	*/
	//---------------------------------------------------------------------------
void CFEMesh::ImportMODFlowGrid(std::string const & fname)
{
	size_t nrows;
	size_t nlayers;
	size_t l;
	std::string aline;
	std::stringstream ss;

	std::vector<double> xx, yy, zz, eval;
	std::vector<bool> layer_act_flag;

	std::ifstream ins(fname.c_str());

	size_t ncols;
	ins >> ncols;
	xx.resize(ncols + 1);
	for (l = 0; l < ncols + 1; l++)
		ins >> xx[l] >> std::ws;
	ins >> nrows;
	yy.resize(nrows + 1);
	for (l = 0; l < nrows + 1; l++)
		ins >> yy[nrows - l] >> std::ws;
	ins >> nlayers;
	zz.resize(nlayers + 1);
	for (size_t l = 0; l < nlayers + 1; l++)
		ins >> zz[l];
	ins >> std::ws;

	//getline(ins, aline);
	bool flag;
	for (size_t i = 0; i < nlayers; i++) {
		for (size_t l = 0; l < nrows; l++) {

			getline(ins, aline);
			ss.str(aline);
			for (size_t k = 0; k < ncols; k++) {
				ss >> flag;
				layer_act_flag.push_back(flag);
			}
			ss.clear();
		}
		getline(ins, aline);
	}

	getline(ins, aline);

	ins >> ncols;
	ins >> nrows;
	ins >> nlayers;

	getline(ins, aline);
	getline(ins, aline);

	double el;
	for (size_t i = 0; i < nlayers + 1; i++) {
		for (l = 0; l < ncols * nrows; l++) {
			ins >> el;
			eval.push_back(el);
		}
	}

	int ii, jj;
	long counter = 0, lnn = (ncols + 1) * (nrows + 1);
	long nn[8];

	std::vector<long> marked_nodes((nlayers + 1) * (ncols + 1) * (nrows + 1));
	marked_nodes[l] = -1;

	for (size_t i = 0; i < nlayers; i++) {
		for (l = 0; l < nrows; l++) {
			for (size_t k = 0; k < ncols; k++) {
				if (layer_act_flag[counter]) {
					nn[0] = i * lnn + l * ncols + k;
					nn[1] = i * lnn + (l + 1) * ncols + k;
					nn[2] = i * lnn + (l + 1) * ncols + k + 1;
					nn[3] = i * lnn + l * ncols + k + 1;

					nn[4] = (i + 1) * lnn + l * ncols + k;
					nn[5] = (i + 1) * lnn + (l + 1) * ncols + k;
					nn[6] = (i + 1) * lnn + (l + 1) * ncols + k + 1;
					nn[7] = (i + 1) * lnn + l * ncols + k + 1;

					for (ii = 0; ii < 8; ii++)
						marked_nodes[nn[ii]] = 1;

				}
				counter++;
			}
		}
	}

	// Create nodes:
	CNode *node;
	counter = 0;
	long ri, ci, li, nnly;
	for (size_t l = 0; l < marked_nodes.size(); l++) {
		li = l / lnn;
		nnly = l % lnn;

		ri = nnly / ncols;
		ci = nnly % ncols;

		if (marked_nodes[l] > 0) {
			marked_nodes[l] = counter;
			node = new CNode(counter, xx[ci], yy[ri], 0.);
			nod_vector.push_back(node);
			counter++;
		}
	}
	xx.clear();
	yy.clear();

	// Element
	CElem *elem;
	double z;
	counter = 0;
	for (size_t i = 0; i < nlayers; i++) {
		for (l = 0; l < nrows; l++) {
			for (size_t k = 0; k < ncols; k++) {
				if (layer_act_flag[counter]) {

					elem = new CElem((long) ele_vector.size());
					ele_vector.push_back(elem);

					elem->setElementProperties(MshElemType::HEXAHEDRON);

					elem->patch_index = i;
					elem->nodes_index.resize(elem->nnodes);
					elem->nodes.resize(elem->nnodes);

					elem->boundary_type = 'I';
					// Initialize topological properties
					elem->neighbors.resize(elem->nfaces);
					for (ii = 0; ii < elem->nfaces; ii++)
						elem->neighbors[ii] = NULL;
					elem->edges.resize(elem->nedges);
					elem->edges_orientation.resize(elem->nedges);
					for (ii = 0; ii < elem->nedges; ii++) {
						elem->edges[ii] = NULL;
						elem->edges_orientation[ii] = 1;
					}
					elem->area = 1.0;

					nn[0] = i * lnn + l * ncols + k;
					nn[1] = i * lnn + (l + 1) * ncols + k;
					nn[2] = i * lnn + (l + 1) * ncols + k + 1;
					nn[3] = i * lnn + l * ncols + k + 1;

					nn[4] = (i + 1) * lnn + l * ncols + k;
					nn[5] = (i + 1) * lnn + (l + 1) * ncols + k;
					nn[6] = (i + 1) * lnn + (l + 1) * ncols + k + 1;
					nn[7] = (i + 1) * lnn + l * ncols + k + 1;
					for (ii = 0; ii < 8; ii++) {
						elem->SetNodeIndex(ii, marked_nodes[nn[ii]]);
						elem->nodes[ii] = nod_vector[elem->GetNodeIndex(ii)];
					}

					ri = nrows - l;
					for (jj = 0; jj < 2; jj++) {
						li = nrows * ncols * (i + jj) + k * nrows + ri;
						for (ii = 4* jj ; ii < 4* (jj +1); ii++) {
							z = nod_vector[marked_nodes[nn[ii]]]->Z()+ eval[li];
							nod_vector[marked_nodes[nn[ii]]]->SetZ(z);
						}
					}

				}
				counter++;
			}
		}
	}

	eval.clear();
	layer_act_flag.clear();
	marked_nodes.clear();

	ConstructGrid();
	for(l=0; l<nod_vector.size(); l++) {
			node = nod_vector[l];
			z = node->Z()/(double)node->connected_elements.size();
			node->SetZ(z);
	}
	ins.close();
}

   //---------------------------------------------------------------------------
   /*!
      \brief Covert GIS shapfile defined cells into triangle/quadrilateral elements

      The GIS shapfile defined cells is read, converted into
       triangle/quadrilateral elements, and the the coverted mesh are written in an OGS syntax file.

      \param fname The file name.

       12/2009 WW
       01/2010 TF changed signature to const & in order to save a string copy
   added serveral std:: in order to avoid name space polution
   */
   //---------------------------------------------------------------------------
   void CFEMesh::ConvertShapeCells(std::string const & fname)
   {
      int i;
      long l, k, counter;
      double  x, y;

      ReadShapeFile(fname);

      long ll, nsize;
      long neighbor[4];

      nsize = (nrows+1) * (ncols+1);
      std::vector<bool> mark(nsize);
      std::vector<long> node_index(nsize);

      CElem *elem;
      CNode *point = NULL;
      bool onboundary;

#define no_need_quad

      for(l=0; l<nsize; l++)
      {
         mark[l] = false;
         node_index[l] = -1;
      }
      for(l=0; l<nrows; l++)
      {
         ll = nrows-1-l;
         for(k=0; k<ncols; k++)
         {
            counter = l*ncols+k;
            if(fabs(zz[counter]-ndata_v)>DBL_MIN)
            {
               elem = new CElem((long)ele_vector.size());
               ele_vector.push_back(elem);

               // Boundary
               onboundary = false;
               if(k==0||(k==ncols-1)||l==0||(l==nrows-1))
                  onboundary = true;
               else
               {
                  neighbor[0] =  l*ncols+k-1;
                  neighbor[1] =  l*ncols+k+1;
                  neighbor[2] =  (l-1)*ncols+k;
                  neighbor[3] =  (l+1)*ncols+k;
                  for(i=0; i<4; i++)
                  {
                     if(fabs(zz[ neighbor[i]]-ndata_v)<DBL_MIN)
                     {
                        onboundary = true;
                        break;
                     }

                  }
               }

               elem->nnodes = 4;
               elem->nodes_index.resize(elem->nnodes);
               elem->nodes.resize(elem->nnodes);

#ifdef need_quad
               elem->nnodesHQ = 9;
               elem->ele_dim = 2;
               elem->nfaces = 4;
               elem->nedges = 4;
               elem->geo_type = MshElemType::QUAD;

               elem->patch_index = 0;

               elem->boundary_type='I';
               // Initialize topological properties
               elem->neighbors.resize(elem->nfaces);
               for(i=0; i<elem->nfaces; i++)
                  elem->neighbors[i] = NULL;
               elem->edges.resize(elem->nedges);
               elem->edges_orientation.resize(elem->nedges);
               for(i=0; i<elem->nedges; i++)
               {
                  elem->edges[i] = NULL;
                  elem->edges_orientation[i] = 1;
               }
               elem->area = 1.0;
#endif

               // nn-->   neighbor
               neighbor[0] =  l*(ncols+1)+k;
               neighbor[1] =  (l+1)*(ncols+1)+k;
               neighbor[2] =  (l+1)*(ncols+1)+k+1;
               neighbor[3] =  l*(ncols+1)+k+1;
               for(i=0; i<4; i++)
               {
                  elem->SetNodeIndex(i, neighbor[i]);
                  mark[neighbor[i]] = true;
               }

            }
         }
      }

      // Nodes
      for(l=0; l<=nrows; l++)
      {
         ll = nrows-l;
         for(k=0; k<=ncols; k++)
         {
            counter = l*(ncols+1)+k;
            if(mark[counter])
            {
               x = x0+csize*k;
               y = y0+csize*ll;
               point =  new CNode((long)nod_vector.size(), x, y);
               nod_vector.push_back(point);
               node_index[counter] = point->GetIndex();
            }

         }
      }

      for(l=0; l<(long)ele_vector.size(); l++)
      {
         elem = ele_vector[l];
         for(i=0; i<4; i++)
         {
            elem->nodes[i] = nod_vector[node_index[elem->GetNodeIndex(i)]];
            elem->SetNodeIndex(i, elem->nodes[i]->GetIndex());

         }
      }

      /// To Karsten: If tri elements, activate the following
#ifndef need_quad
      //-------------------------------------------------------------
      // Up here is the quadrilateral elements

      // Elems
      //
      std::vector<CElem*> ele_tri_v;
      std::vector<CNode*> nodes_e(3);
      CElem *ele_tri;
      for(l=0; l<(long)ele_vector.size(); l++)
      {
         elem = ele_vector[l];

         ele_tri = new CElem((long)ele_tri_v.size());
         ele_tri_v.push_back(ele_tri);

         ele_tri->nnodes = 3;
         ele_tri->nnodesHQ = 6;
         ele_tri->ele_dim = 2;
         ele_tri->nfaces = 3;
         ele_tri->nedges = 3;
         ele_tri->geo_type = MshElemType::TRIANGLE;

         ele_tri->patch_index = 0;
         ele_tri->nodes_index.resize(ele_tri->nnodes);
         ele_tri->nodes.resize(ele_tri->nnodes);

         ele_tri->boundary_type='I';
         // Initialize topological properties
         ele_tri->neighbors.resize(ele_tri->nfaces);
         for(i=0; i<ele_tri->nfaces; i++)
            ele_tri->neighbors[i] = NULL;
         ele_tri->edges.resize(ele_tri->nedges);
         ele_tri->edges_orientation.resize(ele_tri->nedges);
         for(i=0; i<ele_tri->nedges; i++)
         {
            ele_tri->edges[i] = NULL;
            ele_tri->edges_orientation[i] = 1;
         }
         ele_tri->area = 1.0;

         for(i=0; i<3; i++)
         {
            ele_tri->nodes[i] = elem->nodes[i];
            ele_tri->nodes_index[i] = ele_tri->nodes[i]->GetIndex();

            nodes_e[i] = elem->nodes[(i+2)%4];

         }

         // Modify the existing quad.
         elem->nnodes = 3;
         elem->nnodesHQ = 6;
         elem->ele_dim = 2;
         elem->nfaces = 3;
         elem->nedges = 3;
         elem->geo_type = MshElemType::TRIANGLE;

         elem->patch_index = 0;
         elem->nodes.resize(elem->nnodes);
         elem->nodes_index.resize(elem->nnodes);

         elem->boundary_type='I';
         // Initialize topological properties
         elem->neighbors.resize(elem->nfaces);
         for(i=0; i<elem->nfaces; i++)
            elem->neighbors[i] = NULL;
         elem->edges.resize(elem->nedges);
         elem->edges_orientation.resize(elem->nedges);
         for(i=0; i<elem->nedges; i++)
         {
            elem->edges[i] = NULL;
            elem->edges_orientation[i] = 1;
         }
         elem->area = 1.0;

         for(i=0; i<3; i++)
         {
            elem->nodes[i] = nodes_e[i] ;
            elem->nodes_index[i] = elem->nodes[i]->GetIndex();

         }

      }

      for(l=0; l<(long)ele_tri_v.size(); l++)
         ele_vector.push_back(ele_tri_v[l]);
      ele_tri_v.clear();
      //-------------------------------------------------------------
      ///  To Karsten: If quad elements, comment ends here
#endif                                      //ifndef need_quad

      ConstructGrid();
      node_index.clear();
      mark.clear();

   }
   //---------------------------------------------------------------------------
   /*!
      \brief Read GIS shapfile

      \param fname The file name.

       03/2010 WW
   */
   inline void CFEMesh::ReadShapeFile(std::string const & fname)
   {
      long l;

      std::string aline;
      std::stringstream ss;

      zz.clear();

      std::ifstream ins(fname.c_str());
      if(!ins.good())
      {
         std::cout << "Can not find file " << std::endl;
         return ;
      }

      getline(ins, aline);
      ss.str(aline);
      ss>> aline >> ncols;
      ss.clear();

      getline(ins, aline);
      ss.str(aline);
      ss>> aline >> nrows;
      ss.clear();

      getline(ins, aline);
      ss.str(aline);
      ss>> aline >> x0;
      ss.clear();

      getline(ins, aline);
      ss.str(aline);
      ss>> aline >> y0;
      ss.clear();

      getline(ins, aline);
      ss.str(aline);
      ss>> aline >> csize;
      ss.clear();

      getline(ins, aline);
      ss.str(aline);
      ss>> aline >> ndata_v;
      ss.clear();

      zz.resize(nrows * ncols);
      for(l=0; l<(long)zz.size(); l++)
         ins>>zz[l];
      ins.close();
   }

   /*!
      \brief Read GIS shapfile that stores the precipitation data

       Assume the precipitation data is stored in a GIS shapfile. This funtion read the data
       and then performs the numerical integration in order to take the precipitation
       as the Neumman boundary conditions and transform them into the finite element node values.

      \param fname The input file name.
      \param ofname The output file name.
      \param ratio The ration of precipitation to the infiltration.

   03/2010  WW

   */
   void CFEMesh::Precipitation2NeumannBC(std::string const & fname, std::string const & ofname, bool is2D, double ratio)
   {
      int k;
      long i, nx, ny;
      double  x, y;
      double node_val[8];

      CNode *node;
      CElem* elem = NULL;
      CElement* fem = NULL;
      fem = new CElement(GetCoordinateFlag());

      std::vector<double> val;
      val.resize(NodesNumber_Linear);
      for(i=0; i<(long)nod_vector.size(); i++)
      {
         nod_vector[i]->SetMark(false);
         val[i] = 0.0;
      }

      //
      ReadShapeFile(fname);

      //		CElem *elem (face_vector[i]);
      //		CElem *own_elem (elem->owner);
      //
     
      long ele_num = (long)face_vector.size();
      if(is2D)
         ele_num = (long)ele_vector.size();

      for(i=0; i<ele_num; i++)
      {
         if(is2D)
           elem = ele_vector[i]; 
         else
           elem = face_vector[i];
         if(!elem->GetMark())
            continue;

         for(k=0; k<elem->nnodes; k++)
            node_val[k] = 0.0;

         for(k=0; k<elem->nnodes; k++)
         {
            node = elem->nodes[k];
            x = node->X();
            y = node->Y();

            nx = (long)((x-x0)/csize);
            ny = (long)((y-y0)/csize);
            ny = nrows-ny;
            if(ny<0) ny = 0;
            if(ny>nrows) ny = nrows;

            if(nx*csize+x0>=x)  nx -= 1;
            if(ny*csize+y0>=y)  ny -= 1;
            if(nx>=ncols-1) nx = ncols-2;
            if(ny>=nrows-1) ny = nrows-2;
            if(nx<0) nx = 0;
            if(ny<0) ny = 0;

            node_val[k] = zz[ncols*ny+nx];
            if(fabs(node_val[k]-ndata_v)<DBL_MIN)
               node_val[k] = 0.;
         }

         elem->ComputeVolume();
         fem->setOrder(getOrder()+1);
         fem->ConfigElement(elem);
         fem->FaceIntegration(node_val);
         for(k=0; k<elem->nnodes; k++)
         {
            node = elem->nodes[k];
            node->SetMark(true);
            val[node->GetIndex()] += node_val[k];
         }
      }
      
      if(is2D)
      {
         std::ofstream ofile(ofname.c_str(), std::ios::trunc);
         ofile.setf(std::ios::scientific,std::ios::floatfield);
         ofile.precision(14);
         for(i=0; i<(long)val.size(); i++)
         {
            if(fabs(val[i])>DBL_MIN)
              ofile<<i<<"  "<<val[i]<<std::endl;
         }
         ofile<<"#STOP"<<std::endl;
         ofile.close();   

         return;
      }



      //
      std::ofstream ofile_bin(ofname.c_str(), std::ios::trunc|std::ios::binary);
      ofile_bin.setf(std::ios::scientific,std::ios::floatfield);
      ofile_bin.precision(14);

      long counter = 0;
      for(i=0; i<(long)nod_vector.size(); i++)
      {
         if(!nod_vector[i]->GetMark())
            continue;
         counter++;
         val[i] *= ratio;                    // Assuming the unit of precipitation is mm/day
      }
      ofile_bin.write((char*)(&counter), sizeof(counter));

      for(i=0; i<(long)nod_vector.size(); i++)
      {
         node = nod_vector[i];
         if(!node->GetMark())
            continue;
         nx = node->GetIndex();
         ofile_bin.write((char*)(&nx), sizeof(nx));
         ofile_bin.write((char*)(&val[i]), sizeof(val[i]));
      }

      ofile_bin.close();
      delete fem;
      fem = NULL;
      val.clear();
   }

   /*!
       Find element nodes on the top surface of a mesh domain
       07.06.2010
       By WW
   */
   void CFEMesh::MarkInterface_mHM_Hydro_3D(bool quad, bool top_bottom_nodes_outp)
   {
      int k;
      long i;
      CElem *elem;
      CElem *own_elem;
      CNode *node;
      double cent[3];
      double fac;
      double tol = 1.e-9;
      enum top_bot {top, bottom, none};
      vector<top_bot> node_mark;

      if(top_bottom_nodes_outp)  
      {  
         /// For output z coordinate of all nodes on the top surface
         /// 13.08.2010. WW
         node_mark.resize(GetNodesNumber(quad));
         for(i=0; i<GetNodesNumber(quad); i++)
            node_mark[i] = none;
      }

      for(i=0; i<(long)face_vector.size(); i++)
      {
         elem = face_vector[i];

         own_elem = elem->owner;

         //// In element
         //		// compute center of mass
         cent[0] = cent[1] = cent[2] = 0.;
         //		const int nodes_number_of_element (own_elem->GetNodesNumber(false));
         for(k=0; k<own_elem->GetNodesNumber(false); k++)
         {
            node = own_elem->nodes[k];
            cent[0] += node->X();
            cent[1] += node->Y();
            cent[2] += node->Z();
         }
         for(k=0; k<3; k++)
            cent[k] /= (double)own_elem->GetNodesNumber(false);

         node = elem->nodes[0];
         cent[0] -= node->X();
         cent[1] -= node->Y();
         cent[2] -= node->Z();
         NormalizeVector(cent,3);
         elem->ComputeVolume();
         elem->FillTransformMatrix();
         /// Compute the normal to this surface element
         fac =  cent[0]*(*elem->tranform_tensor)(0,2)
            +cent[1]*(*elem->tranform_tensor)(1,2)
            +cent[2]*(*elem->tranform_tensor)(2,2);
         if(fac>0.0) fac = -1.0;
         else fac = 1.0;
         //////

         /// If n.z>0
         if((*elem->tranform_tensor)(2,2)*fac>tol)
         {
            elem->SetMark(true);      /// Only mark nodes on the top
            for(k=0; k<3; k++)
               (*elem->tranform_tensor)(k,2) *= fac;

            if(top_bottom_nodes_outp) 
            { 
               for(k=0; k<elem->GetNodesNumber(quad); k++)
                  node_mark[elem->nodes[k]->GetIndex()] = top;
            }
         }
         else if ((*elem->tranform_tensor)(2,2)*fac<-tol)
         {
            elem->SetMark(false);
            for(k=0; k<3; k++)
               (*elem->tranform_tensor)(k,2) *= fac;

            if(top_bottom_nodes_outp) 
            { 
               for(k=0; k<elem->GetNodesNumber(quad); k++)
                  node_mark[elem->nodes[k]->GetIndex()] = bottom;
            }
         }
         else
            elem->SetMark(false);
      }

      if(top_bottom_nodes_outp)
      {
         string ccc = FileName + "_nodes_on_top.asc";
         ofstream ofile_asci(ccc.c_str(), ios::trunc);
         ccc = FileName + "_nodes_on_bottom.asc";
         ofstream ofile_asci1(ccc.c_str(), ios::trunc);

         for(i=0; i<(long)nod_vector.size(); i++)
         {
            node = nod_vector[i];
            if(node_mark[i] == top)
              ofile_asci<< node->GetIndex()<<" "<< "0. " <<endl;
            else if(node_mark[i] == bottom)
              ofile_asci1<< node->GetIndex()<<" "<< "0. " <<endl;
         }
         ofile_asci<<"#STOP"<<endl;
         ofile_asci1<<"#STOP"<<endl;
         ofile_asci.clear();
         ofile_asci1.clear();
         ofile_asci.close();
         ofile_asci1.close();
      }

   }

   /*!
      \brief Transform GIS shapfile stored precitation data into finite element node values

       Assume the precipitation data is stored in GIS shapfiles. This funtion read the data
       and then performs the numerical integration for each  GIS shapfile.

       06/2010  WW

   */
   void CFEMesh::mHM2NeumannBC()
   {
      double ratio, step;

      std::string aline;
      std::stringstream ss;

      std::string fname = FileName+".pcp";

      std::ifstream ins(fname.c_str());
      if(!ins.good())
      {
         std::cout<<"Can not open file "<<fname<<std::endl;
         return;
      }

      ConstructGrid();

      MarkInterface_mHM_Hydro_3D();

      std::string key, uname, ofname;
      //char stro[1024];

      getline(ins, aline);
      ss.str(aline);
      ss>> key >> uname;
      ss.clear();

      getline(ins, aline);
      ss.str(aline);
      ss>> key >> ratio;
      ss.clear();

      step = 0.;

      std::string infiltration_files;
      infiltration_files = FileName+".ifl";
      std::ofstream infil(infiltration_files.c_str(), std::ios::trunc);
	  int counter = 0;
      while(!ins.eof())
      {
         getline(ins, aline);
         ss.str(aline);
         ss>> key;
         ss.clear();

         if(key.size()==0)                        // An empty line
            continue;

         if(key.find("#STOP")!=std::string::npos)
            break;

         //sprintf(stro, "%f",step);
         // ofname = stro;

		 // Write the headers of shape file
		 if (counter == 0)
		 {
			 const std::string asc_fname = FilePath + key;
             std::ifstream ins_asc(asc_fname.c_str());
			 std::string keyw;
			 double value;
			 for (int k=0; k< 6; k++)
			 {
                ins_asc >> keyw >> value >> std::ws;
				infil << keyw << " " << value << std::endl;
			 }
			 ins_asc.close();
		 }
		 counter++;

         ofname = FilePath+key+".bin";
         infil<<step<<" "<<key+".bin"<<std::endl;

         key = FilePath+key;
         Precipitation2NeumannBC(key, ofname, false, ratio);

         step += 1.0;

      }
      infil<<"#STOP"<<std::endl;
      infil.close();

   }

   /*!
    Comptute \int {f} a dA on top surface.

      WW. 29.11.2010
   */
   void CFEMesh::TopSurfaceIntegration()
   {
      int k;
      long i, nx;
      double node_val[8];

      CNode *node;
      CElem* elem = NULL;
      CElement* fem = NULL;

      ConstructGrid();
      MarkInterface_mHM_Hydro_3D();

      fem = new CElement(GetCoordinateFlag());
      std::vector<double> val;
      val.resize(NodesNumber_Linear);
      for(i=0; i<(long)nod_vector.size(); i++)
      {
         nod_vector[i]->SetMark(false);
         val[i] = 0.0;
      }

      //
      //	//
      std::string ofname = FileName+"_top_surface_Neumann_BC.txt";
      std::ofstream ofile_asci(ofname.c_str(), std::ios::trunc);
      ofile_asci.setf(std::ios::scientific,std::ios::floatfield);
      ofile_asci.precision(14);
      //	//
      //

      for(i=0; i<(long)face_vector.size(); i++)
      {
         elem = face_vector[i];
         if(!elem->GetMark())
            continue;

         for(k=0; k<elem->nnodes; k++)
            node_val[k] = 1.0;

         elem->ComputeVolume();
         fem->setOrder(getOrder()+1);
         fem->ConfigElement(elem);
         fem->FaceIntegration(node_val);
         for(k=0; k<elem->nnodes; k++)
         {
            node = elem->nodes[k];
            node->SetMark(true);
            val[node->GetIndex()] += node_val[k];
         }
      }

      for(i=0; i<(long)nod_vector.size(); i++)
      {
         node = nod_vector[i];
         if(!node->GetMark())
            continue;
         nx = node->GetIndex();

         ofile_asci<< nx <<" "<< val[i] << std::endl;

      }
      ofile_asci<<"#STOP "<< std::endl;

      ofile_asci.close();
      delete fem;
      fem = NULL;
      val.clear();
   }


   void CFEMesh::Output_Z_TopSurface(bool select_points)
   {
      int k;
      long i, nx;
      vector<CNode *> points;


      CNode *node;
      CElem* elem = NULL;
      CElement* fem = NULL;

      ConstructGrid();
      MarkInterface_mHM_Hydro_3D();

      fem = new CElement(GetCoordinateFlag());
      std::vector<double> val;
      val.resize(NodesNumber_Linear);
      for(i=0; i<(long)nod_vector.size(); i++)
      {
         nod_vector[i]->SetMark(false);
         val[i] = 0.0;
      }

      //
      //	//
      std::string ofname = FileName+"_z_coordinate_of_top_surface.txt";
      std::ofstream ofile_asci(ofname.c_str(), std::ios::trunc);
      ofile_asci.setf(std::ios::scientific,std::ios::floatfield);
      ofile_asci.precision(14);
      //	//
      //

      for(i=0; i<(long)face_vector.size(); i++)
      {
         elem = face_vector[i];
         if(!elem->GetMark())
            continue;

         for(k=0; k<elem->nnodes; k++)
         {
            node = elem->nodes[k];
            node->SetMark(true);
         }
      }

      if(select_points)
	  {
          std::string aline;
          std::stringstream ss;
          double xx, yy;
		  cout<<"Input file of points: "<<endl;
          string filename;
          cin>>filename;

		  filename = FilePath + filename;
		  ifstream is_pnt(filename.c_str());

		  while(!is_pnt.eof())
		  {
             getline(is_pnt, aline);
			 if(aline.find("#STOP")!=string::npos)
                break;

             ss.str(aline);
             ss>> xx>>yy;
             ss.clear();
             CNode *new_pnt = new CNode(points.size()-1, xx, yy, 0.);
			 points.push_back(new_pnt);
		  }

		  size_t j = 0;
          for(i=0; i<(long)nod_vector.size(); i++)
          {
            node = nod_vector[i];
            if(!node->GetMark())
               continue;
            nx = node->GetIndex();

            for(j=0; j<points.size(); j++)
			{
                xx = points[j]->X()-node->X();
                yy = points[j]->Y()-node->Y();
//                if(sqrt(xx*xx+yy*yy)<sqrt(DBL_EPSILON))                   
                if(sqrt(xx*xx+yy*yy)<0.1)                   
                   ofile_asci<< nx <<" "<< node->Z() << std::endl;
			}
         }

/*
          for(j=0; j<points.size(); j++)
          {
              xx = points[j]->X();
              yy = points[j]->X();

             for(i=0; i<(long)nod_vector.size(); i++)
             {
               node = nod_vector[i];
               if(!node->GetMark())
                  continue;
               nx = node->GetIndex();

               dx = xx - node->X();
               dy = yy - node->Y();
               if(sqrt(dx*dx+dy*dy)<sqrt(DBL_EPSILON))                   
//               if(sqrt(xx*xx+yy*yy)<0.1)
			   {
				   ofile_asci<< nx <<" "<< node->Z() << std::endl;
				   break;

			   }
			}
         }

*/

	  }
	  else
	  {
         for(i=0; i<(long)nod_vector.size(); i++)
         {
            node = nod_vector[i];
            if(!node->GetMark())
               continue;
            nx = node->GetIndex();

            ofile_asci<< nx <<" "<< node->Z() << std::endl;

         }
	  }
      ofile_asci<<"#STOP "<< std::endl;



      ofile_asci.close();
      delete fem;
      fem = NULL;
      val.clear();
   }


#ifdef USE_HydSysMshGen
   /**************************************************************************
   MSHLib-Method:
   Task:
   05/2009 WW
   **************************************************************************/
   void  CFEMesh::HydroSysMeshGenerator(string fname, const int nlayers, const double thickness, int mapping)
   {
      fstream gs_out;
      int k, mat_num;
      long i;
      string name = fname+"_hrdosys.msh";
      string deli = " ";
      CNode *a_node = NULL;
      CElem *an_ele = NULL;

      gs_out.open(name.c_str(), ios::out);
      gs_out.setf(ios::scientific, ios::floatfield);
      setw(14);
      gs_out.precision(14);

      if(mapping)
      {
         // To do
      }

      gs_out<<"#FEM_MSH\n$PCS_TYPE\nOVERLAND_FLOW\n$NODES"<<endl;
      gs_out<<(long)nod_vector.size()<<endl;
      for(i=0; i<(long)nod_vector.size(); i++)
         nod_vector[i]->Write(gs_out);
      gs_out<<"$ELEMENTS"<<endl;
      gs_out<<(long)ele_vector.size()<<endl;
      mat_num = 0;
      for(i=0; i<(long)ele_vector.size(); i++)
      {
         an_ele = ele_vector[i];
         an_ele->WriteGSmsh(gs_out);
         k =  an_ele->GetPatchIndex();
         if(k>mat_num) mat_num = k;
      }
      mat_num++;

      gs_out<<"#FEM_MSH\n$PCS_TYPE\n RICHARDS_FLOW\n$GEO_TYPE\nPOLYLINE REGIONAL\n$GEO_NAME\nREGIONAL"<<endl;
      gs_out<<"$NODES\n"<<(nlayers+1)*(long)nod_vector.size()<<endl;
      double  seg = thickness/(double)nlayers;
      double depth;
      long l, size_nodes_msh_t;
      size_nodes_msh_t = (long)nod_vector.size();
      for(i=0; i<size_nodes_msh_t; i++)
      {
         for(k=0; k<=nlayers; k++)
         {
            depth = seg*(double)k;
            a_node = nod_vector[i];
            gs_out<<k+i*(nlayers+1)<<deli
               <<a_node->X()<<deli
               <<a_node->Y()<<deli
               <<a_node->Z()-depth<<deli<<endl;
         }
      }
      gs_out<<"$ELEMENTS"<<endl;
      gs_out<<size_nodes_msh_t*nlayers<<endl;
      l = 0;
      int mat_index;
      for(i=0; i<size_nodes_msh_t; i++)
      {
         // mat_index = ele_vector[nod_vector[i]->connected_elements[0]]->GetPatchIndex();
         mat_index = mat_num;
         for(k=0; k<nlayers; k++)
         {
            l = k+(nlayers+1)*i;
            gs_out<<k+nlayers*i<<deli<<mat_index<<deli<<"line"<<deli;
            gs_out<< l<<deli<<l+1<<endl;
            mat_index++;
         }
      }
      gs_out<<"$LAYER\n"<<nlayers<<endl;
      //
      gs_out<<"$BORDERS"<<endl;
      gs_out<<"SECTOR_GROUND\n"<<size_nodes_msh_t<<endl;
      for(i=0; i<size_nodes_msh_t; i++)
      {
         k=nlayers;
         gs_out<<k+i*(nlayers+1)<<deli<<endl;
      }

      mat_num += nlayers;
      gs_out<<"#FEM_MSH\n$PCS_TYPE\nGROUNDWATER_FLOW\n$NODES"<<endl;
      gs_out<<(long)nod_vector.size()<<endl;
      for(i=0; i<(long)nod_vector.size(); i++)
      {
         a_node = nod_vector[i];
         a_node->SetZ(a_node->Z()-thickness);
         a_node->Write(gs_out);
         a_node->SetZ(a_node->Z()+thickness);
      }
      gs_out<<"$ELEMENTS"<<endl;
      gs_out<<(long)ele_vector.size()<<endl;
      for(i=0; i<(long)ele_vector.size(); i++)
      {
         an_ele = ele_vector[i];
         an_ele->SetPatchIndex(an_ele->GetPatchIndex()+mat_num);
         an_ele->WriteGSmsh(gs_out);
      }
      gs_out<<"$BORDERS"<<endl;
      gs_out<<"SECTOR_SOIL\n"<<(long)nod_vector.size()<<endl;
      for(i=0; i<(long)nod_vector.size(); i++)
         gs_out<<i<<deli<<endl;

      gs_out<<"#STOP"<<endl;
      gs_out.close();
   }
#endif



#ifdef NON_PROCESS
void  CFEMesh::TwoDRechargeDatato3DSurface(std::string data_file, std::string file_path,
	                                       const double ratio, const int unit_type)
{
  std::string aline;
  std::stringstream ss;
  std::string key;
  std::ifstream ins;
  size_t i, ncols0, ncols, nrows;
  CNode *node;

  vector<string> names;

  if(max_ele_dim > 2)
  {
      cerr << "This is not a 2D mesh. "; 
	  exit(1);
  }

  double x_min, y_min, x_max, y_max;
  x_min = DBL_MAX;
  y_min = DBL_MAX;
  x_max = DBL_MIN;
  y_max = DBL_MIN;

  size_t nnodes = nod_vector.size();
  CElem * elem;

  //ins.open(data_file.c_str());
  //getline(ins, aline);  
  //std::replace( aline.begin(), aline.end(), '/', '_'); 
  std::string sub_key;

  double *x;


  x_min = DBL_MAX;
  y_min = DBL_MAX;
  x_max = DBL_MIN;
  y_max = DBL_MIN;


  ins.open(data_file.c_str());
  getline(ins, aline);  

  ss.str(aline);

  ss >> key >> key ;
  vector<string> year_names;
  vector<size_t> year_mark;

  int counter = 0;
  int counter_y = 0;
  while(!ss.eof())
  {
     ss >> key >> sub_key;
	 if (!key.empty() && !sub_key.empty()  )
	 {
         key = key + sub_key + ".asc";
		 names.push_back(key);
         if(unit_type == 2) // Year
		 {
         
             bool exsit = false; 
             for(size_t k=0; k<year_names.size(); k++)
     		 {
	     		 if(year_names[k].find(sub_key) != string::npos)
				 {
                       exsit = true;
                       break;
				 }
			 }
             if(!exsit)
			 {
                 if(counter > 0) 
		           counter_y++;
                 year_mark.push_back(counter_y);
                 year_names.push_back(sub_key);
			 }
			 else
                year_mark.push_back(counter_y);

 		 }

		 counter++;
	 }		  
  }
  ss.clear();
  ncols0 = counter;

  if(unit_type == 2) // Year
  {
	  names.clear();
	  names.resize(year_names.size());
      for(i=0; i<year_names.size(); i++)
	  {
          names[i] = "year"+year_names[i]+".asc";
	  }
  }



  ncols = names.size();



  getline(ins, aline);  
  ss.str(aline);
  ss >> nrows;
  ss.clear();
 
  x = new double[2*nrows]; 

  //double *rchg_data =  new double[nrows * ncols];   
  vector<double*> rchg_data(ncols);
  vector<long> ele_ID(nrows); 
  for(i=0; i<nrows; i++)
  {
	  ele_ID[i] = -1;
  }
  for(i=0; i<ncols; i++)
  {
     double *new_row = new double[nrows];
     rchg_data[i] = new_row; 
  }


  if(unit_type == 2) // Year
  {
     for(i=0; i<ncols; i++)
     {
        double *new_row = new double[nrows];
        for(size_t k=0; k < nrows; k++)
           rchg_data[i][k] = 0.0; 
     }
  }

  cout <<"*** Reading data ... "<<endl;
  double xx, yy;
  for(i=0; i<nrows; i++)
  {
     getline(ins, aline);  
     ss.str(aline);
     ss >> xx;
	 ss >> yy;
	  /*
	  if(xx < x_min)
        x_min = xx; 
	  if(xx > x_max)
        x_max = xx; 
	  if(yy < y_min)
        y_min = yy; 
	  if(yy > x_max)
        y_max = yy; 
     */
     x[2*i] = xx;
     x[2*i+1] = yy;

     for(size_t k=0; k<ncols0; k++)
	 {
        double val;
        ss >> val; 
        if(unit_type == 2)
		{
            double *col_dat = rchg_data[year_mark[k]];
			col_dat[i] += val;
		}
		else
		{
            double *col_dat = rchg_data[k];
		    col_dat[i] = val; //rchg_data[i*ncols+k];
		}
	 }
	 ss.clear();
  }


  // Search element IDs.
  double volume, volume0;
  double x0buff[3];
  double x1buff[3];
  double x2buff[3];
  double x3buff[3];
  double x4buff[3];


  x0buff[2] = 0.;
  x1buff[2] = 0.;
  x2buff[2] = 0.;
  x3buff[2] = 0.;
  x4buff[2] = 0.;

  cout <<"*** Finding element ID (takes long time)... "<<endl;
  for(size_t j=0; j<nrows; j++)
  {
     x0buff[0] = x[2*j];
     x0buff[1] = x[2*j+1];

     for(i=0; i<ele_vector.size(); i++)
     {
        elem = ele_vector[i];


		CNode *n1 = nod_vector[elem->GetNodeIndex(0)];
		CNode *n2 = nod_vector[elem->GetNodeIndex(1)];
		CNode *n3 = nod_vector[elem->GetNodeIndex(2)];
		double const *x1 = n1->getData();
		double const *x2 = n2->getData();
		double const *x3 = n3->getData();

        for(int k=0; k<2; k++)
        {
            x1buff[k] = x1[k];
            x2buff[k] = x2[k];
			x3buff[k] = x3[k];
     	}

		if(elem->GetElementType() == MshElemType::TRIANGLE)
		{
            volume =   ComputeDetTri(x1buff, x2buff, x0buff)
                     + ComputeDetTri(x2buff, x3buff, x0buff)
                     + ComputeDetTri(x3buff, x1buff, x0buff);

			volume0 = ComputeDetTri(x1buff, x2buff, x3buff);
		} 
		else if(elem->GetElementType() == MshElemType::QUAD)
		{
       		CNode *n4 = nod_vector[elem->GetNodeIndex(3)];
		    double const *x4 = n4->getData();
            for(int k=0; k<2; k++)
			{
				x4buff[k] = x4[k];
			}
            volume =   ComputeDetTri(x1buff, x2buff, x0buff)
                     + ComputeDetTri(x2buff, x3buff, x0buff)
                     + ComputeDetTri(x3buff, x4buff, x0buff)
                     + ComputeDetTri(x4buff, x1buff, x0buff);
            volume0 =   ComputeDetTri(x1buff, x2buff, x3buff)
                     + ComputeDetTri(x3buff, x4buff, x1buff);
		}

		if(fabs(volume - volume0) < DBL_EPSILON)
		{
			ele_ID[j] =  static_cast<long> (i);

			//cout <<j<< "ele ID  "<<ele_ID[j] << endl;
            elem->SetVolume(volume0);
            break;
		}
     }
  }

  

  std::string infiltration_files;
  infiltration_files = data_file+"_monthly.ifl";
  if(unit_type == 2) // Year
     infiltration_files = data_file+"_annual.ifl";

  std::ofstream infil(infiltration_files.c_str(), std::ios::trunc);

  infil <<"DIRECT"<<endl;


  cout <<"*** Performing face integration ... "<<endl;

  ofstream os_r;
  setw(14);
  os_r.precision(14);
  
  ofstream os; // VTK
  setw(14);
  os.precision(14);

  double node_val[8];

  CElement* fem = NULL;
  fem = new CElement(23);


  const size_t nn = nod_vector.size(); 
  std::vector<double> val(nn);
  std::vector<double> e_recharge(ele_vector.size());
  //std::vector<size_t> node_id(nn);
  //for(i=0; i<nn; i++)
 // {
  //    node_id[i] = nod_vector[i]->GetEquationIndex();
 // }

  for(i=0; i<ncols; i++)
  {

     //string fname = names[i] + ".bin"; 
	 //os.open(fname.c_str(), ios::trunc|ios::binary);


     double *col_dat = rchg_data[i];
	 for(size_t k=0; k<nn; k++)
     {
        val[k] = 0.0;
     }



////////////////////

	 for (size_t k=0; k<nrows; k++)
	 {
        if(ele_ID[k] == -1 )
           continue;

        elem =  ele_vector[ele_ID[k]]; 
		e_recharge[ele_ID[k]] = col_dat[k];

         if(fabs(col_dat[k]) < DBL_EPSILON)
            continue;

         for(int m=0; m<elem->nnodes; m++)
         {
            node_val[m] = col_dat[k];
         }
		 
         fem->setOrder(1);
		 fem->SetMeshElement(elem);
         fem->SetNNodes();
		 fem->ConfigNumerics(elem->GetElementType());
         fem->FaceIntegration(node_val);
		 
         for(int m=0; m<elem->nnodes; m++)
         {
            node = nod_vector[elem->GetNodeIndex(m)];
			val[node->GetIndex()] += node_val[m];

         }
		 
	 }


     double sum =  0.;
	 size_t counter = 0;
	 for(size_t k=0; k<nn; k++)
     {
		sum += val[k];

        if(fabs(val[k]) > DBL_MIN )
           counter++; 
     }
	 if(fabs(sum) < DBL_MIN)
		 continue;

     infil<<i<<" "<<names[i] <<std::endl;


     string fname = file_path + names[i];
	 os_r.open(fname.c_str(), ios::trunc);

	 os_r << counter << " 1\n";
	 for(size_t k=0; k<nn; k++)
     {
        if(fabs(val[k]) > DBL_MIN )
		{ 
           val[k] *= ratio;
	       os_r << nod_vector[k]->GetEquationIndex() <<" "<< val[k] <<"\n";
		}
     }
	 os_r << endl;

	 os_r.clear();
	 os_r.close();




	 fname = file_path + names[i] + ".vtk"; 
	 os.open(fname.c_str(), ios::trunc);
     os<<"# vtk DataFile Version 4.0\nRecharge2Source\nASCII\n"<<endl;
     os<<"DATASET UNSTRUCTURED_GRID"<<endl;
     os<<"POINTS "<<nn<<" double"<<endl;
	 for(size_t k=0; k<nn; k++)
     {
        node = nod_vector[k];
        os<<node->X()<<" "<<node->Y()<<" "<<node->Z()<<endl;
     }
     long size = (long)ele_vector.size();
     for(size_t k=0; k<(long)ele_vector.size(); k++)
        size += ele_vector[k]->GetNodesNumber(false);
     os<<"\nCELLS "<<(long)ele_vector.size()<<" "<<size<<endl;
     // CELLs
     for(size_t k=0; k<(long)ele_vector.size(); k++)
     {
        elem = ele_vector[k];
        os<<elem->GetNodesNumber(false)<<"  ";
        for(int m=0; m<elem->GetNodesNumber(false); m++)
			os << elem->GetNodeIndex(m)<< " ";
        os << endl;
      }
      os << endl;

      // CELL types
      os << "CELL_TYPES " << (int) ele_vector.size() << endl;
      for(size_t k=0; k<(long)ele_vector.size(); k++)
      {
         elem = ele_vector[k];
		 if(elem->GetElementType() == MshElemType::TRIANGLE)
			os<<5<< " "<<endl;
		 else if(elem->GetElementType() == MshElemType::QUAD)
			os<<9<< " "<<endl;
      }
      os << endl;
      // Coordinates
      os<<"POINT_DATA "<<nn<<endl;
      os<<"SCALARS Integrated_Value double 1\nLOOKUP_TABLE default"<<endl;
	 for(size_t k=0; k<nn; k++)
         os<<val[k]<<endl;
       os<<endl;
     size_t size_ne = ele_vector.size();
     os<<"CELL_DATA "<<size_ne<<endl;
     os<<"SCALARS Recharge double 1\nLOOKUP_TABLE default"<<endl;
     for(size_t k=0; k<size_ne; k++)
        os<<e_recharge[k]<<endl;




	  os.clear();
	  os.close();


	 /*
     ofile_bin.write((char*)(&counter), sizeof(counter));

      for(i=0; i<(long)nod_vector.size(); i++)
      {
         node = nod_vector[i];
         if(!node->GetMark())
            continue;
         nx = node->GetIndex();
         ofile_bin.write((char*)(&nx), sizeof(nx));
         ofile_bin.write((char*)(&val[i]), sizeof(val[i]));
      }

      ofile_bin.close();

	 */





  } 
  

  infil<<"#STOP"<<std::endl;
  infil.close();

  for(i=0; i<ncols; i++)
  {
    delete [] rchg_data[i];   
  }

  delete [] x;
  delete fem;
  fem = NULL;
  val.clear();


}


/*! \brief As the function name

*/
void CFEMesh::Write_Surface_GMSH(std::string ofile)
{
   CNode *node;
   CElem *elem;
   long i;
   string fname;

   // For output
   long nn_sfc, ne_sfc;
   nn_sfc = 0;
   ne_sfc = 0;


   // For output surface
   fname =  ofile+"_sfc.msh";
   ofstream os(fname.c_str(),ios::trunc);
   setw(14);
   os.precision(14);

   for(i=0; i<(long)nod_vector.size(); i++)
	   nod_vector[i]->SetMark(false);

   for(i=0; i<(long)face_vector.size(); i++)
   {
      elem = face_vector[i];
	  if(!elem->GetMark())
         continue;

	  for(int k=0; k< elem->GetNodesNumber(false); k++)
	  {
		  elem->GetNode(k)->SetMark(true);
	  }
      ne_sfc++;
   }
   for(i=0; i<(long)nod_vector.size(); i++)
   {
	   if(nod_vector[i]->GetMark())
         nn_sfc++;
   }
   os<<"$MeshFormat\n2 0 8\n$EndMeshFormat\n$Nodes"<<endl;
   os<<nn_sfc<<endl;
   for(i=0; i<(long)nod_vector.size(); i++)
   {
      node = nod_vector[i];
	  if(!node->GetMark())
         continue;
      os<<node->GetIndex()+1<<" "<<node->X()<<" "<<node->Y()<<" "<<node->Z()<<endl;
   }
   os<<"$EndNodes\n$Elements"<<endl;
   os<<ne_sfc<<endl;
   for(i=0; i<(long)face_vector.size(); i++)
   {
      elem = face_vector[i];
	  if(!elem->GetMark())
         continue;
      elem->WriteGmsh(os);
   }
   os<<"$EndElements"<<endl;

   os.close();
   ///
}


#endif



}


// namespace Mesh_Group
