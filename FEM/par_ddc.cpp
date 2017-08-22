/**************************************************************************
PARLib - Object:
Task:
Programing:
07/2004 OK Implementation
07/2004 OK Version 1 untill 3.9.17OK6
07/2004 OK Version 2 from   3.9.17OK7
07/2006 WW Local topology, High order nodes
last modified:
**************************************************************************/
//---- MPI Parallel --------------
// There is a name conflict between stdio.h and the MPI C++ binding
// with respect to the names SEEK_SET, SEEK_CUR, and SEEK_END.  MPI
// wants these in the MPI namespace, but stdio.h will #define these
// to integer values.  #undef'ing these can cause obscure problems
// with other include files (such as iostream), so we instead use
// #error to indicate a fatal error.  Users can either #undef
// the names before including mpi.h or include mpi.h *before* stdio.h
// or iostream.
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS)
//#undef SEEK_SET  //WW
//#undef SEEK_END  //WW
//#undef SEEK_CUR  //WW
#include<mpi.h>
int size;
int myrank;
int mysize;
char t_fname[3];
double time_ele_paral;
#endif
//---- MPI Parallel --------------

#include <math.h>
// C++ STL
#include <iostream>
using namespace std;
#include "par_ddc.h"
// FEM-Makros
#include "makros.h"
#ifndef NEW_EQS                                   //WW. 11.2008
#include "matrix.h"
#endif
#include "files0.h"
#include "rf_num_new.h"
#include "gs_project.h"
#ifdef NEW_EQS
// Solver WW
#include "matrix_class.h"
#include "equation_class.h"
#endif
vector<CPARDomain*>dom_vector;
vector<int> node_connected_doms;                  //This will be removed after sparse class is finished WW

/**************************************************************************
STRLib-Method:
Task:
Programing:
07/2004 OK Implementation
**************************************************************************/
bool KeywordFound(string line)
{
   string hash("#");
   if(line.find(hash)!=string::npos)
      return true;
   else
      return false;
}


/**************************************************************************
STRLib-Method:
Task:
Programing:
07/2004 OK Implementation
**************************************************************************/
bool SubKeywordFound(string line)
{
   string dollar("$");
   if(line.find(dollar)!=string::npos)
      return true;
   else
      return false;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
07/2004 OK Implementation
07/2004 OK Version 1 untill 3.9.17OK6
07/2004 OK Version 2 from   3.9.17OK7
10/2005 OK cout
**************************************************************************/
void DOMRead(string file_base_name)
{
   //----------------------------------------------------------------------
   cout << "DOMRead: ";
   //----------------------------------------------------------------------
   CPARDomain *m_dom = NULL;
   char line[MAX_ZEILE];
   string sub_line;
   string line_string;
   string ddc_file_name;
   ios::pos_type position;
   //----------------------------------------------------------------------
   // File handling
   ddc_file_name = file_base_name + DDC_FILE_EXTENSION;
   ifstream ddc_file (ddc_file_name.data(),ios::in);
   if(!ddc_file.good())
   {
      cout << "no DDC file" << endl;
      return;
   }
   ddc_file.seekg(0L,ios::beg);
   //----------------------------------------------------------------------
   // Keyword loop
   while (!ddc_file.eof())
   {
      ddc_file.getline(line,MAX_ZEILE);
      line_string = line;
      //----------------------------------------------------------------------
                                                  // keyword found
      if(line_string.find("#DOMAIN")!=string::npos)
      {
         m_dom = new CPARDomain();
         position = m_dom->Read(&ddc_file);
         dom_vector.push_back(m_dom);
         ddc_file.seekg(position,ios::beg);
      }                                           // keyword found
   }                                              // eof
   //----------------------------------------------------------------------
   cout << dom_vector.size() << " domains" << endl;
   //----------------------------------------------------------------------
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
09/2007 WW Implementation
**************************************************************************/
void CountDoms2Nodes(CRFProcess *m_pcs)
{
   int i,k;
   long j, n_index=0;
   CPARDomain *m_dom = NULL;
   CNode *anode = NULL;
   CElem *elem = NULL;
   CFEMesh *a_msh = m_pcs->m_msh;

   //Average of nodal Neumann BCs contributed by nodes from different domains
   bool quad = false;
   if(m_pcs->type==4||m_pcs->type==41) quad = true;
   long nsize = m_pcs->m_msh->GetNodesNumber(quad);
   node_connected_doms.resize(nsize);

   for(j=0; j<nsize; j++)
      node_connected_doms[j] = 0;
   for(i=0;i<(int)dom_vector.size();i++)
   {
      m_dom = dom_vector[i];
      for(j=0; j<nsize; j++)
         a_msh->nod_vector[j]->SetMark(false);
      for(j=0;j<(long)m_dom->elements.size();j++)
      {
         elem = a_msh->ele_vector[m_dom->elements[j]];
         if(elem->GetMark())
         {
            for(k=0; k<elem->GetNodesNumber(quad); k++)
            {
               n_index = elem->GetNodeIndex(k);
               anode = elem->GetNode(k);
               if(!anode->GetMark())
               {
                  node_connected_doms[n_index] += 1;
                  anode->SetMark(true);
               }
            }
         }
      }
   }
   //
   for(j=0; j<(long)a_msh->ele_vector.size(); j++)
   {
      elem = a_msh->ele_vector[j];
      if(elem->GetMark())
      {
         for(k=0; k<elem->GetNodesNumber(quad); k++)
         {
            anode = elem->GetNode(k);
            anode->SetMark(true);
         }
      }
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
07/2004 OK Implementation
05/2006 WW Fix bugs and add the case of DOF>1
07/2006 WW Find nodes of all neighbors of each node
09/2007 WW Check the nodes whether it is belong to the deactivated elements
10/2010 TF changed access to process type
**************************************************************************/
void DOMCreate()
{
   size_t no_domains = dom_vector.size();
   if(no_domains==0)
      return;
   CPARDomain *m_dom = NULL;
   CRFProcess* m_pcs = NULL;
   bool quadr = false;                            //WW
   size_t i;
   //WW ----- D
   //----------------------------------------------------------------------
   for(i=0;i<pcs_vector.size();i++)
   {
      m_pcs = pcs_vector[i];
      //    if(m_pcs->pcs_type_name.find("DEFORMATION")!=string::npos) { // TF 10/2010
      if(isDeformationProcess (m_pcs->getProcessType ()))
      {
         quadr = true;
         break;
      }
   }
   if(!quadr)
      m_pcs = pcs_vector[0];

   //----------------------------------------------------------------------
   // Create domain nodes
   cout << "->Create DOM" << endl;
   /* // Comment by WW
   for(i=0;i<no_domains;i++){
     m_dom = dom_vector[i];
     m_dom->m_msh = fem_msh_vector[0]; //OK:ToDo
   }
   */

   //----------------------------------------------------------------------
   // Create domain nodes
   cout << "  Create domain nodes" << endl;
   for(i=0;i<no_domains;i++)
   {
      m_dom = dom_vector[i];
      m_dom->m_msh = m_pcs->m_msh;
      cout << "    Domain:" << m_dom->ID << endl;
      m_dom->CreateNodes();
   }
   //----------------------------------------------------------------------
   // Create domain elements
   cout << "  Create domain elements" << endl;
   for(i=0;i<no_domains;i++)
   {
      m_dom = dom_vector[i];
      cout << "    Domain:" << m_dom->ID << endl;
      m_dom->CreateElements(quadr);
   }
   // For find nodes connected to node WW
   long j;
   long nsize = m_pcs->m_msh->GetNodesNumber(true);
   //node_connected_doms.resize(nsize);
   for(j=0; j<nsize; j++)
      node_connected_doms[j] = 0;
   for(i=0;i<dom_vector.size();i++)
   {
      m_dom = dom_vector[i];
      for(j=0;j<(long)m_dom->nodes.size();j++)
         node_connected_doms[m_dom->nodes[j]] += 1;
   }
   //
   // Find nodes of all neighbors of each node. // WW
   // Local topology. WW
   cout << "  Find nodes on borders" << endl;
   FindNodesOnInterface(m_pcs->m_msh, quadr);
   cout << "  Find the connected nodes for each node" << endl;
#ifndef USE_MPI                                //WW
   for(i=0;i<no_domains;i++)
   {
#else
      i = myrank;                                 //WW
#endif
      m_dom = dom_vector[i];
      m_dom->NodeConnectedNodes();
#ifndef USE_MPI                             //WW
   }
#endif
   //----------------------------------------------------------------------
   // Create domain EQS
   cout << "  Create domain EQS" << endl;
#ifdef USE_MPI
   i = myrank;
#else
   for(i=0;i<no_domains;i++)
   {
#endif
      m_dom = dom_vector[i];
      cout << "    Domain:" << m_dom->ID << endl;
#ifdef NEW_EQS
      m_dom->CreateEQS();
#else
      m_dom->CreateEQS(m_pcs);
#endif
      //
#ifndef USE_MPI
   }
#endif
   //----------------------------------------------------------------------
}


//////////////////////////////////////////////////////////////////////////
// CPARDomain

/**************************************************************************
FEMLib-Method:
Task:
Programing:
07/2004 OK Implementation
**************************************************************************/
CPARDomain::CPARDomain(void)
{
   ID =(int)dom_vector.size();
   for(int i=0; i<4; i++)                         //WW
      shift[i] = 0;
   quadratic = false;                             //WW
#ifdef NEW_EQS
   sparse_graph = NULL;
   sparse_graph_H = NULL;
   eqs = NULL;
   eqsH = NULL;
#endif
#if defined(USE_MPI)                           // 13.12.2007 WW
   t_border_nodes = NULL;
   t_border_nodes_size = t_border_nodes_sizeH = 0;
   //
#if defined(NEW_BREDUCE)
   receive_cnt_b = new int[mysize];
   receive_disp_b = new int[mysize];
#endif
   receive_cnt_i = new int[mysize];
   receive_disp_i = new int[mysize];
   receive_cnt = new int[mysize];
   receive_disp = new int[mysize];
#endif
}


CPARDomain::~CPARDomain(void)
{
   elements.clear();
   nodes.clear();
   nodes_inner.clear();
   nodes_halo.clear();
   long i;

   // WW
   for(i=0; i< (long)element_nodes_dom.size(); i++)
   {
      delete [] element_nodes_dom[i];
      element_nodes_dom[i] = NULL;
   }
   for(long i=0; i< (long)node_conneted_nodes.size(); i++)
   {
      delete [] node_conneted_nodes[i];
      node_conneted_nodes[i] = NULL;
   }

   //
#ifdef NEW_EQS
   if(eqs) delete eqs;
   if(eqsH) delete eqsH;
   if(sparse_graph) delete sparse_graph;
   if(sparse_graph_H) delete sparse_graph_H;
#endif

#if defined(USE_MPI)                           // 13.12.2007 WW
   //
   if(t_border_nodes) delete [] t_border_nodes;
   t_border_nodes = NULL;
#if defined(NEW_BREDUCE)
   delete [] receive_cnt_b;
   delete [] receive_disp_b;
   receive_cnt_b = NULL;
   receive_disp_b = NULL;
#endif
   delete [] receive_cnt_i;
   delete [] receive_disp_i;
   delete [] receive_cnt;
   delete [] receive_disp;
   receive_cnt_i = NULL;
   receive_disp_i = NULL;
   receive_cnt = NULL;
   receive_disp = NULL;
#endif
}


/**************************************************************************
FEMLib-Method:
Task: Release partial memory
Programing:
012/2007 WW Implementation
**************************************************************************/
#if defined(USE_MPI)
void CPARDomain::ReleaseMemory()
{
   // WW
   for(long i=0; i< (long)element_nodes_dom.size(); i++)
   {
      delete [] element_nodes_dom[i];
      element_nodes_dom[i] = NULL;
   }
   for(long i=0; i< (long)node_conneted_nodes.size(); i++)
   {
      delete [] node_conneted_nodes[i];
      node_conneted_nodes[i] = NULL;
   }
   if(eqs) delete eqs;
   if(eqsH) delete eqsH;
   if(sparse_graph) delete sparse_graph;
   if(sparse_graph_H) delete sparse_graph_H;
   if(t_border_nodes) delete [] t_border_nodes;
   t_border_nodes = NULL;
#if defined(NEW_BREDUCE)
   delete [] receive_cnt_b;
   delete [] receive_disp_b;
   receive_cnt_b = NULL;
   receive_disp_b = NULL;
#endif
   delete [] receive_cnt_i;
   delete [] receive_disp_i;
   delete [] receive_cnt;
   delete [] receive_disp;
   receive_cnt_i = NULL;
   receive_disp_i = NULL;
   receive_cnt = NULL;
   receive_disp = NULL;
}
#endif
/**************************************************************************
FEMLib-Method:
Task: ST read function
Programing:
07/2004 OK Implementation
07/2004 OK Version 1 untill 3.9.17OK6
07/2004 OK Version 2 from   3.9.17OK7
**************************************************************************/
ios::pos_type CPARDomain::Read(ifstream *ddc_file)
{
   char line[MAX_ZEILE];
   string sub_line;
   string sub_string;
   string cut_string;
   string line_string;
   string delimiter(" ");
   string delimiter_type(";");
   bool new_subkeyword = false;
   string dollar("$");
   bool new_keyword = false;
   string hash("#");
   ios::pos_type position;
   long i;
   //  CElem* m_ele = NULL;
   //======================================================================
   while (!new_keyword)
   {
      position = ddc_file->tellg();
      ddc_file->getline(line,MAX_ZEILE);
      line_string = line;
      if(line_string.find(hash)!=string::npos)
      {
         new_keyword = true;
         break;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$ELEMENTS")!=string::npos)
      {
         while ((!new_keyword)&&(!new_subkeyword))
         {
            position = ddc_file->tellg();
            ddc_file->getline(line,MAX_ZEILE);
            line_string = line;
            if(line_string.find(hash)!=string::npos)
            {
               new_keyword = true;
               break;
            }
            if(line_string.find(dollar)!=string::npos)
            {
               new_subkeyword = true;
               break;
            }
            i = strtol(line,NULL,0);
            elements.push_back(i);
         }
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$NODES_INNER")!=string::npos)
      {
         while (!new_keyword)
         {
            position = ddc_file->tellg();
            ddc_file->getline(line,MAX_ZEILE);
            line_string = line;
            if(line_string.find(hash)!=string::npos)
            {
               new_keyword = true;
               break;
            }
            if(line_string.find(dollar)!=string::npos)
            {
               new_subkeyword = true;
               break;
            }
            i = strtol(line,NULL,0);
            nodes_inner.push_back(i);
         }
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$NODES_BORDER")!=string::npos)
      {
         while (!new_keyword)
         {
            position = ddc_file->tellg();
            ddc_file->getline(line,MAX_ZEILE);
            line_string = line;
            if(line_string.find(hash)!=string::npos)
            {
               new_keyword = true;
               break;
            }
            if(line_string.find(dollar)!=string::npos)
            {
               new_subkeyword = true;
               break;
            }
            i = strtol(line,NULL,0);
            nodes_halo.push_back(i);
         }
      }
      //....................................................................
   }
   //======================================================================
   return position;
}


/**************************************************************************
PARLib-Method:
Task:
Programing:
07/2004 OK Implementation
**************************************************************************/
void CPARDomain::CreateNodes()
{
   long j,k;
   long no_nodes_halo, no_nodes_inner;
   //----------------------------------------------------------------------
   no_nodes_inner = (long)nodes_inner.size();
   for(j=0;j<no_nodes_inner;j++)
   {
      k = nodes_inner[j];
      nodes.push_back(k);
      //cout << nodes[j] << endl;
   }
   //cout << "---" << endl;
   no_nodes_halo = (long)nodes_halo.size();
   for(j=0;j<no_nodes_halo;j++)
   {
      k = nodes_halo[j];
      nodes.push_back(k);
      //cout << nodes[no_nodes_inner+j] << endl;
   }
   nnodes_dom = no_nodes_halo+no_nodes_inner;     //WW
   nnodesHQ_dom = nnodes_dom;

   //----------------------------------------------------------------------
}


/**************************************************************************
PARLib-Method:
Task:
Programing:
07/2004 OK Implementation
05/2006 WW Fix bugs and add the case of DOF>1
09/2007 WW Improve the extremly slow performance of the old code
**************************************************************************/
void CPARDomain::CreateElements(const bool quadr)
{
   //----------------------------------------------------------------------
   if(!m_msh)
      return;
   //----------------------------------------------------------------------
   long i,k;
   int j, nNodes, nNodesHQ;
   long *elem_nodes=NULL;
   Mesh_Group::CElem* m_ele = NULL;
   Mesh_Group::CNode* m_nod = NULL;
   //*** Buffer for acceleration. 14.09.2007 WW:
                                                  // As long buffer
   node_connected_doms.resize((long)m_msh->nod_vector.size());
   for(k=0; k<(long)m_msh->nod_vector.size(); k++)
      node_connected_doms[k] = -1;
   for(k=0; k<(long)nodes.size(); k++)
   {
      i = nodes[k];
      m_nod = m_msh->nod_vector[i];
      node_connected_doms[m_nod->GetIndex()] = k;
   }
   //***
   //----------------------------------------------------------------------
   for(i=0;i<(long)elements.size();i++)
   {
      if(elements[i]>(long)m_msh->ele_vector.size())
      {
         cout << "Warning: no ELE data" << '\n';
         continue;
      }
      m_ele = m_msh->ele_vector[elements[i]];
      nNodes = m_ele->GetNodesNumber(false);      //WW
      nNodesHQ = m_ele->GetNodesNumber(quadr);
      // cout << i << " " << elements[i] << ": ";
      elem_nodes = new long[nNodesHQ];            //WW
      element_nodes_dom.push_back(elem_nodes);    //WW
      for(j=0;j<nNodes;j++)
      {
         m_nod = m_ele->GetNode(j);
                                                  //14.09.2007 WW
         elem_nodes[j] = node_connected_doms[m_nod->GetIndex()];
      }
      //------------------WW
      if(!quadr) continue;
      for(j=nNodes;j<nNodesHQ;j++)
      {
         m_nod = m_ele->GetNode(j);
         //14.09.2007 WW
         k = m_nod->GetIndex();
         if(node_connected_doms[k]>-1)
            elem_nodes[j] = node_connected_doms[k];
         else
         {
            elem_nodes[j] = (long)nodes.size();
            node_connected_doms[k] = elem_nodes[j];
            nodes.push_back(m_nod->GetIndex());
         }
      }
      nnodesHQ_dom = (long) nodes.size();
      //------------------WW
      // cout << endl;
   }
   //
   //----------------------------------------------------------------------
}


/**************************************************************************
MSHLib-Method:
Task:
Programing:
06/2006 WW Implementation
09/2007 WW Improve computational efficiency
**************************************************************************/
void CPARDomain::NodeConnectedNodes()
{
   int k, i_buff;
   long i, j, long_buff;
   CNode* m_nod = NULL;
   vector<long> nodes2node;
   // node_connected_doms as buffer to accelarate the computation
   // 14.09.2007 WW
   for(i=0; i<(long)m_msh->nod_vector.size(); i++)
      node_connected_doms[i] = -1;
   for(j=0; j<(long)nodes.size(); j++)
   {
      i = nodes[j];
      m_nod = m_msh->nod_vector[i];
      node_connected_doms[m_nod->GetIndex()] = j;
   }

   //----------------------------------------------------------------------
   for(i=0;i<(long)nodes.size();i++)
   {
      m_nod = m_msh->nod_vector[nodes[i]];
      nodes2node.clear();
      for(k=0; k<(int)m_nod->connected_nodes.size(); k++)
      {
         long_buff = m_nod->connected_nodes[k];
         j = node_connected_doms[long_buff];
         if(j>-1)
            nodes2node.push_back(j);
      }
      i_buff = (int)nodes2node.size();
      long  *nodes_to_node = new long[i_buff];
      for(k=0; k<i_buff; k++)
         nodes_to_node[k] = nodes2node[k];
      node_conneted_nodes.push_back(nodes_to_node);
      num_nodes2_node.push_back(i_buff);
   }
   //----------------------------------------------------------------------
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
07/2004 OK Implementation
**************************************************************************/
long CPARDomain::GetDOMNode(long global_node)
{
   long i;
   long no_nodes = (long)nodes.size();
   for(i=0;i<no_nodes;i++)
   {
      if(nodes[i]==global_node)
         return i;
   }
   return -1;
}


#ifdef NEW_EQS
/**************************************************************************
PARLib-Method:
Task:      Create local subdomain sparse matrix
Programing:
12/2007 WW Implementation
**************************************************************************/
void CPARDomain::CreateSparseTable()
{
   // Symmetry case is skipped.
   // 1. Sparse_graph_H for high order interpolation. Up to now, deformation
   if(nnodesHQ_dom!=nnodes_dom)
      sparse_graph_H = new SparseTable(*this, true);
   // 2. M coupled with other processes with linear element
   if(sparse_graph_H)
   {
      if((int)pcs_vector.size()>1)
         sparse_graph = new SparseTable(*this, false);
   }
   // 3. For process with linear elements
   else
      sparse_graph = new SparseTable(*this, false);

   //sparse_graph->Write();
   //sparse_graph_H->Write();
   //
   //ofstream Dum("sparse.txt", ios::out);
   //sparse_graph_H->Write(Dum);}

   /*

     //TEST
   string test = "rank";
   char stro[1028];
   sprintf(stro, "%d",myrank);
   string test1 = test+(string)stro+"sparse.txt";
      ofstream Dum(test1.c_str(), ios::out); // WW
      sparse_graph->Write(Dum);   Dum.close();

   */

}


/**************************************************************************
PARLib-Method:
Task:  Create local equation system
Programing:
12/2007 WW Implementation
**************************************************************************/
void CPARDomain::CreateEQS()
{
   size_t dof_nonDM=1, dof_DM=1;

   const size_t pcs_vector_size (pcs_vector.size());
   for(size_t i=0;i<pcs_vector_size;i++)
   {
      CRFProcess *pcs (pcs_vector[i]);
      if(pcs->type==22)                           // Monolithic TH2
         dof_nonDM = pcs->GetPrimaryVNumber();
      if(pcs->type==4 || pcs->type==41)           // Deformation
         dof_DM = pcs->GetPrimaryVNumber();
      else                                        // Monolithic scheme for the process with linear elements
      {
         if(dof_nonDM < pcs->GetPrimaryVNumber())
            dof_nonDM = pcs->GetPrimaryVNumber();
      }
   }

   CreateSparseTable();

   if(sparse_graph)
      eqs = new Linear_EQS(*sparse_graph, dof_nonDM);
   if(sparse_graph_H)
      eqsH = new Linear_EQS(*sparse_graph_H, dof_DM);
}


/**************************************************************************
PARLib-Method:
Task:
Programing:
12/2007 WW Implementation
**************************************************************************/
void CPARDomain::InitialEQS(CRFProcess *m_pcs)
{
   long size = nnodes_dom;
   Linear_EQS *this_eqs = NULL;
   if(m_pcs->type==4||m_pcs->type==41)
   {
      size = nnodesHQ_dom;
      this_eqs = eqsH;
   }
   else
      this_eqs = eqs;
   for(size_t i=0; i<m_pcs->GetPrimaryVNumber(); i++)
      shift[i] = i*size;
   //
   this_eqs->Initialize();
}


//
#else
/**************************************************************************
PARLib-Method:
Task:
Programing:
07/2004 OK Implementation
05/2006 WW Fix bugs and add the case of DOF>1
**************************************************************************/
void CPARDomain::CreateEQS(CRFProcess *m_pcs)
{

   long no_nodes = (long)nodes.size();
   int dof = m_pcs->GetPrimaryVNumber();

   if(m_pcs->type==4||m_pcs->type==41)
   {
      for(size_t i=0; i<m_pcs->GetPrimaryVNumber(); i++)
         shift[i] = i*nnodesHQ_dom;

   }
   if(m_pcs->type==4)
   {
      eqs = CreateLinearSolverDim(m_pcs->m_num->ls_storage_method,dof,dof*nnodesHQ_dom);
      //    InitializeLinearSolver(eqs,m_num);
   }
   else if(m_pcs->type==41)
   {
      eqs = CreateLinearSolverDim(m_pcs->m_num->ls_storage_method,dof,dof*nnodesHQ_dom+nnodes_dom);

   }
   else
                                                  //WW.
      eqs = CreateLinearSolver(m_pcs->m_num->ls_storage_method, no_nodes);
   //  eqs = CreateLinearSolver(0,no_nodes);
   //InitializeLinearSolver(m_dom->eqs,NULL,NULL,NULL,m_dom->lsp_name);
   InitLinearSolver(eqs);
}
#endif
/**************************************************************************
FEMLib-Method:
Task:
Programing:
07/2004 OK Implementation
**************************************************************************/
void CPARDomain::CalcElementMatrices(CRFProcess* m_pcs)
{
   m_pcs = m_pcs;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
07/2004 OK Implementation
**************************************************************************/
/* //WW
void CPARDomain::AssembleMatrix(CRFProcess* m_pcs)
{
  long i;
  //----------------------------------------------------------------------
  SetZeroLinearSolver(eqs);
//MXDumpGLS("AssembleMatrix1.txt",1,eqs->b,eqs->x);
  //----------------------------------------------------------------------
  long no_elements = (long)elements.size();
  for(i=0;i<no_elements;i++){
    // virtual function PCSAssembleMatrix(i)
//MakeElementEntryEQS_ASM(elements[i]->global_number,eqs->b,NULL,this);
//WW    MakeElementEntryEQS_ASM(i,eqs->b,NULL,this,m_pcs);
//MXDumpGLS("AssembleMatrix1.txt",1,eqs->b,eqs->x);
}
}
*/
/**************************************************************************
FEMLib-Method:
Task:
Programing:
07/2004 OK Implementation
**************************************************************************/
bool NodeExists(long node,vector<long>node_vector)
{
   long i;
   long no_nodes = (long)node_vector.size();
   for(i=0;i<no_nodes;i++)
   {
      if(node==node_vector[i])
         return true;
   }
   return false;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2005 OK Implementation
09/2005 OK MSH
last modification:
**************************************************************************/
void CPARDomain::WriteTecplot(string msh_name)
{
   long i;
   string element_type;
   //----------------------------------------------------------------------
   // GSP
   CGSProject* m_gsp = NULL;
   m_gsp = GSPGetMember("msh");
   if(!m_gsp)
   {
      return;
   }
   //--------------------------------------------------------------------
   // file handling
   char name[10];
   sprintf(name,"%i",ID);
   string dom_name = "DOMAIN";
   dom_name += name;
   string dom_file_name = m_gsp->path + dom_name + TEC_FILE_EXTENSION;
   fstream dom_file (dom_file_name.data(),ios::trunc|ios::out);
   dom_file.setf(ios::scientific,ios::floatfield);
   dom_file.precision(12);
   //--------------------------------------------------------------------
   // MSH
   CFEMesh* m_msh = NULL;
   CNode* m_nod = NULL;
   CElem* m_ele = NULL;
   m_msh = FEMGet(msh_name);
   if(!m_msh)
      return;
   //--------------------------------------------------------------------
   if (!dom_file.good()) return;
   dom_file.seekg(0L,ios::beg);
   //--------------------------------------------------------------------
   for(i=0;i<(long)elements.size();i++)
   {
      m_ele = m_msh->ele_vector[elements[i]];
      if(!m_ele)
         continue;
      switch(m_ele->GetElementType())
      {
         case MshElemType::LINE:
            element_type = "ET = QUADRILATERAL";
            break;
         case MshElemType::QUAD:
            element_type = "ET = QUADRILATERAL";
            break;
         case MshElemType::HEXAHEDRON:
            element_type = "ET = BRICK";
            break;
         case MshElemType::TRIANGLE:
            element_type = "ET = TRIANGLE";
            break;
         case MshElemType::TETRAHEDRON:
            element_type = "ET = TETRAHEDRON";
            break;
         case MshElemType::PRISM:
            element_type = "ET = BRICK";
            break;
         default:
            std::cerr << "CPARDomain::WriteTecplot MshElemType not handled" << std::endl;
      }
   }
   //--------------------------------------------------------------------
   dom_file << "VARIABLES = X,Y,Z,DOM" << endl;
   long no_nodes = (long)m_msh->nod_vector.size();
   dom_file << "ZONE T = " << dom_name << ", " \
      << "N = " << no_nodes << ", " \
      << "E = " << (long)elements.size() << ", " \
      << "F = FEPOINT" << ", " << element_type << endl;
   //......................................................................
   for(i=0;i<no_nodes;i++)
   {
      m_nod = m_msh->nod_vector[i];
      dom_file \
         << m_nod->X() << " " << m_nod->Y() << " " << m_nod->Z() << " " \
         << ID << endl;
   }
   //......................................................................
   for(i=0;i<(long)elements.size();i++)
   {
      m_ele = m_msh->ele_vector[elements[i]];
      if(!m_ele)
         continue;
      switch(m_ele->GetElementType())
      {
         case MshElemType::LINE:
            dom_file \
               << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[0]+1 << endl;
            element_type = "ET = QUADRILATERAL";
            break;
         case MshElemType::QUAD:
            dom_file \
               << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[2]+1 << " " << m_ele->nodes_index[3]+1 << endl;
            element_type = "ET = QUADRILATERAL";
            break;
         case MshElemType::HEXAHEDRON:
            dom_file \
               << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[2]+1 << " " << m_ele->nodes_index[3]+1 << " " \
               << m_ele->nodes_index[4]+1 << " " << m_ele->nodes_index[5]+1 << " " << m_ele->nodes_index[6]+1 << " " << m_ele->nodes_index[7]+1 << endl;
            element_type = "ET = BRICK";
            break;
         case MshElemType::TRIANGLE:
            dom_file \
               << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[2]+1 << endl;
            element_type = "ET = TRIANGLE";
            break;
         case MshElemType::TETRAHEDRON:
            dom_file \
               << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[2]+1 << " " << m_ele->nodes_index[3]+1 << endl;
            element_type = "ET = TETRAHEDRON";
            break;
         case MshElemType::PRISM:
            dom_file \
               << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[2]+1 << " " \
               << m_ele->nodes_index[3]+1 << " " << m_ele->nodes_index[3]+1 << " " << m_ele->nodes_index[4]+1 << " " << m_ele->nodes_index[5]+1 << endl;
            element_type = "ET = BRICK";
            break;
         default:
            std::cerr << "CPARDomain::WriteTecplot MshElemType not handled" << std::endl;
      }
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
10/2005 OK Implementation
last modification:
**************************************************************************/
void DOMWriteTecplot(string msh_name)
{
   CPARDomain *m_dom = NULL;
   for(int i=0;i<(int)dom_vector.size();i++)
   {
      m_dom = dom_vector[i];
      m_dom->WriteTecplot(msh_name);
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
07/2006 WW Implementation
last modification:
**************************************************************************/
void FindNodesOnInterface(CFEMesh *m_msh, bool quadr)
{
   // int k;
   long i, j, nnodes_gl, g_index, nnodes_l;
   long l_buff = 0, l_buff1 = 0;

   long *elem_nodes=NULL;
   //
   Mesh_Group::CElem* m_ele = NULL;
   //
   CPARDomain *m_dom = NULL;
   vector<long> boundary_nodes;
   vector<long> boundary_nodes_HQ;
   vector<long> inner_nodes_HQ;
   vector<long> inner_nodes;
   vector<long> long_buffer;
   vector<long> bc_buffer;
   vector<long> dom_bc_buffer;
   vector<long> dom_bc_bufferHQ;
   //
   nnodes_gl = m_msh->GetNodesNumber(true);
   nnodes_l = m_msh->GetNodesNumber(false);
   //
   bc_buffer.resize(nnodes_gl);
   //
#if defined(USE_MPI)                           //13.12.2007
   long overlapped_entry_size = 0;
#endif
   for(i=0;i<nnodes_gl;i++)
   {
      if(node_connected_doms[i]>1.0)
      {
         // Mapp BC entry-->global node array
         bc_buffer[i] = (long)long_buffer.size();
         long_buffer.push_back(m_msh->nod_vector[i]->GetIndex());
#ifdef USE_MPI
         if(m_msh->nod_vector[i]->GetIndex()<m_msh-> GetNodesNumber(false))
            overlapped_entry_size = (long)long_buffer.size();
#endif
      }
      else
         bc_buffer[i] = -1;
   }
   //
#if defined(USE_MPI)                           // 13.12.2007
   // Total border nodes
   m_dom = dom_vector[myrank];
   m_dom->t_border_nodes_size = overlapped_entry_size;
   m_dom->t_border_nodes_sizeH = (long)long_buffer.size();
   m_dom->t_border_nodes = new long[m_dom->t_border_nodes_sizeH];
   for(i=0;i<m_dom->t_border_nodes_sizeH;i++)
      m_dom->t_border_nodes [i] = long_buffer[i];
#endif

   // Sort
#ifndef USE_MPI
   for(int k=0;k<(int)dom_vector.size();k++)
   {
      int myrank = k;
#endif
      m_dom = dom_vector[myrank];
      //
      boundary_nodes.clear();
      inner_nodes.clear();
      boundary_nodes_HQ.clear();
      inner_nodes_HQ.clear();
      long_buffer.clear();
      long_buffer.resize((long)m_dom->nodes.size());
      dom_bc_buffer.clear();
      dom_bc_bufferHQ.clear();

      for(i=0;i<(long)m_dom->nodes.size();i++)
      {
         g_index = m_dom->nodes[i];
         if(node_connected_doms[ g_index]>1.0)
         {
            if(g_index>=nnodes_l)
            {
               boundary_nodes_HQ.push_back(i);
               dom_bc_bufferHQ.push_back(bc_buffer[g_index]);
               long_buffer[i] = -(long)boundary_nodes_HQ.size()-nnodes_gl;
            }
            else
            {
               boundary_nodes.push_back(i);
               dom_bc_buffer.push_back(bc_buffer[g_index]);
               long_buffer[i] = -(long)boundary_nodes.size();
            }
         }
         else
         {
            if(g_index>=nnodes_l)
            {
               long_buffer[i] = (long)inner_nodes_HQ.size()+nnodes_gl;
               inner_nodes_HQ.push_back(i);
            }
            else
            {
               long_buffer[i] = (long)inner_nodes.size();
               inner_nodes.push_back(i);
            }
         }
      }
      //
      m_dom->num_inner_nodes = (long)inner_nodes.size();
      m_dom->num_inner_nodesHQ = (long)inner_nodes_HQ.size();
      m_dom->num_boundary_nodes = (long)boundary_nodes.size();
      m_dom->num_boundary_nodesHQ = (long)boundary_nodes_HQ.size();
      // Sort for high order nodes

      m_dom->nodes_inner.clear();
      m_dom->nodes_halo.clear();
      for(i=0;i<m_dom->num_inner_nodes;i++)
         m_dom->nodes_inner.push_back(m_dom->nodes[inner_nodes[i]]);
      for(i=0;i<(long)inner_nodes_HQ.size();i++)
         m_dom->nodes_inner.push_back(m_dom->nodes[inner_nodes_HQ[i]]);
      //
      for(i=0;i<m_dom->num_boundary_nodes;i++)
         m_dom->nodes_halo.push_back(m_dom->nodes[boundary_nodes[i]]);
      for(i=0;i<(long)boundary_nodes_HQ.size();i++)
         m_dom->nodes_halo.push_back(m_dom->nodes[boundary_nodes_HQ[i]]);
      //
      m_dom->nodes.clear();
      m_dom->nodes.resize(m_dom->nnodesHQ_dom);
      // First interior nodes, then interface nodes
      j=0;
      for(i=0;i<m_dom->num_inner_nodes;i++)
      {
         m_dom->nodes[i] = m_dom->nodes_inner[i];
         //       m_dom->nodes_inner[i] = i;
      }
      j += m_dom->num_inner_nodes;
      for(i=0;i<m_dom->num_boundary_nodes;i++)
      {
         m_dom->nodes[i+j] = m_dom->nodes_halo[i];
         //      m_dom->nodes_halo[i] = i+j;
      }
      j += m_dom->num_boundary_nodes;
      for(i=0;i<(long)inner_nodes_HQ.size();i++)
      {
         m_dom->nodes[i+j] = m_dom->nodes_inner[i+m_dom->num_inner_nodes];
         //     m_dom->nodes_inner[i+m_dom->num_inner_nodes] = i+j;
      }
      j += (long)inner_nodes_HQ.size();
      for(i=0;i<(long)boundary_nodes_HQ.size();i++)
      {
         m_dom->nodes[i+j] = m_dom->nodes_halo[i+m_dom->num_boundary_nodes];
         //      m_dom->nodes_halo[i+m_dom->num_boundary_nodes] = i+j;
      }

      for(i=0;i<(long)m_dom->nodes.size();i++)
      {
         l_buff = long_buffer[i];
         if(l_buff<0)                             //interface nodes
         {
            if(-l_buff-nnodes_gl>0)               //HQ nodes
               l_buff1 = m_dom->nnodes_dom+(long)inner_nodes_HQ.size()-l_buff-nnodes_gl-1;
            else
               l_buff1 = (long)inner_nodes.size()-l_buff-1;
            //             l_buff1 = m_dom->num_inner_nodesHQ-l_buff-1;
         }
         else
         {
            if(l_buff-nnodes_gl>=0)               //HQ nodes
               l_buff1 = m_dom->nnodes_dom +l_buff-nnodes_gl;
            else
               l_buff1 = l_buff;

         }
         long_buffer[i] = l_buff1;
      }
      //
#ifdef USE_MPI                              //WW
      m_dom->FillBorderNodeConnectDom(node_connected_doms);
#endif
      //
      m_dom->nodes_inner.clear();
      m_dom->nodes_halo.clear();
      // Mapping the local index to global BC array, overlapped_entry.
      //
      for(i=0;i<m_dom->num_boundary_nodes;i++)
         m_dom->nodes_halo.push_back(dom_bc_buffer[i]);
      for(i=0;i<(long)boundary_nodes_HQ.size();i++)
         m_dom->nodes_halo.push_back(dom_bc_bufferHQ[i]);
      //----------------------------------------------------------------------
      for(i=0;i<(long)m_dom->elements.size();i++)
      {
         m_ele = m_msh->ele_vector[m_dom->elements[i]];
         elem_nodes = m_dom->element_nodes_dom[i];
         for(j=0;j<m_ele->GetNodesNumber(quadr);j++)
         {
            l_buff = elem_nodes[j];
            elem_nodes[j]=long_buffer[l_buff];
         }
      }

#ifndef USE_MPI
   }
#endif
}


/**************************************************************************
DDCLib-Function
07/2007 OK Encapsulation
10/2010 TF changed access to process type
***************************************************************************/
void DDCCreate()
{
   //----------------------------------------------------------------------
   // DDC
   if(dom_vector.size()>0)
   {
      //WW ----- Domain decomposition ------------------
      int i;
      int no_processes =(int)pcs_vector.size();
      CRFProcess* m_pcs = NULL;
      bool DOF_gt_one = false;
      //----------------------------------------------------------------------
      for(i=0;i<no_processes;i++)
      {
         m_pcs = pcs_vector[i];
         //if(m_pcs->pcs_type_name.find("DEFORMATION")!=string::npos) { // TF 10/2010
         if(m_pcs->getProcessType () == DEFORMATION || m_pcs->getProcessType() == DEFORMATION_FLOW)
         {
            DOF_gt_one = true;
            break;
         }
      }
      if(!DOF_gt_one)
         m_pcs = pcs_vector[0];
      // -----------------------
      DOMCreate();
      //
      for(i=0;i<no_processes;i++)
      {
         m_pcs = pcs_vector[i];
         // Config boundary conditions for domain decomposition
         m_pcs->SetBoundaryConditionSubDomain();  //WW
      }
      //
      node_connected_doms.clear();
   }
   // PA PCSProcessDependencies();
}


#if defined(USE_MPI)                              //WW
//------------------------For parallel solvers------------------------------
/*************************************************************************
GeoSys-Function:
Task:
Programming:
02/2008 WW Implementation
**************************************************************************/
void CPARDomain::FillBorderNodeConnectDom(vector<int> allnodes_doms)
{
   long i, ig;
   int k;
   b_start[0] = num_inner_nodes;
   b_end[0] = num_inner_nodes+num_boundary_nodes;
   nq = 1;
   //
   if(nnodesHQ_dom>nnodes_dom)
   {
      //
      nq = 2;
      b_start[1] = b_end[0] + num_inner_nodesHQ;
      b_end[1] = nnodesHQ_dom;
   }
   //
   //
   for(k=0; k<nq; k++)
   {
      for(i=b_start[k]; i<b_end[k]; i++)
      {
         ig = nodes[i];
         bnode_connected_dom.push_back(allnodes_doms[ig]);
      }
   }

   /*
   string test = "rank";
   static char stro[102];

   sprintf(stro, "%d",myrank);
   string test1 = test+(string)stro+"dom.txt";
   ofstream Dum(test1.c_str(), ios::out); // WW
   Dum<<b_start[0]<<endl;
   Dum<<b_end[0]<<endl;
   Dum<<b_start[1]<<endl;
   Dum<<b_end[1]<<endl;
   for(i=0;i<bnode_connected_dom.size();i++)
   Dum<<bnode_connected_dom[i]<<endl;
   exit(1);
   */
}


/*************************************************************************
GeoSys-Function:
Task:
Programming:
12/2007 WW Implementation
**************************************************************************/
void CPARDomain::ConfigEQS(CNumerics *m_num, const long n, bool quad)
{
   int i;
   long dim = 0;
   i_start[0] = 0;
   i_end[0] = num_inner_nodes;                    //Number of interior nodes
   b_start[0] =0;
   b_end[0] = num_boundary_nodes;
   n_shift[0] = num_inner_nodes;
   long inner_size = num_inner_nodes;
   long border_size = num_boundary_nodes;
   quadratic = quad;
   //
   double cpu_time_local = -MPI_Wtime();

   /*

   //TEST_MPI
   string test = "rank";
   static char stro[102];

   sprintf(stro, "%d",myrank);
   string test1 = test+(string)stro+"Assemble.txt";
   ofstream Dum(test1.c_str(), ios::out); // WW
   Dum<<"Mysize" <<mysize<<"  "<<quadratic<<endl;

   Dum.close();
   MPI_Finalize();
   exit(1);
   */

   if(quadratic)
   {
      n_loc = nnodesHQ_dom;
      nq = 2;
      n_bc = t_border_nodes_sizeH;
      i_start[1] = i_end[0]+num_boundary_nodes;
      i_end[1] = i_start[1]+num_inner_nodesHQ;    //Number of interior nodes
      //
      b_start[1] = b_end[0];
      b_end[1] = b_start[1]+num_boundary_nodesHQ;
      n_shift[1] =  n_shift[0]+ num_inner_nodesHQ;
      //
      dof = eqsH->DOF();
      eqsH->SetDomain(this);
      eqsH->ConfigNumerics(m_num, n);
      inner_size += num_inner_nodesHQ;
      border_size += num_boundary_nodesHQ;
   }
   else
   {
      n_loc = nnodes_dom;
      nq = 1;
      n_bc = t_border_nodes_size;
      dof = eqs->DOF();
      eqs->SetDomain(this);
      eqs->ConfigNumerics(m_num, n);
   }
   //  Concatenate index
   inner_size *= dof;
   border_size *= dof;
   for(i=0; i<mysize; i++)
   {
      receive_cnt[i] = 1;
      receive_disp[i] = i;
   }
   //
   // receive_cnt_i[]: number of subdomain inner nodes in the concatenated array
   MPI_Allgatherv ( &inner_size, 1, MPI_INT, receive_cnt_i, receive_cnt, receive_disp,
      MPI_INT, MPI_COMM_WORLD );
   inner_size = 0;
   for(i=0; i<mysize; i++)
   {
      receive_disp_i[i] = inner_size;
      inner_size += receive_cnt_i[i];
   }
#if defined(NEW_BREDUCE)
   // receive_cnt_b[]: number of subdomain border nodes in the concatenated array
   MPI_Allgatherv ( &border_size, 1, MPI_INT, receive_cnt_b, receive_cnt, receive_disp,
      MPI_INT, MPI_COMM_WORLD );
   border_size = 0;
   for(i=0; i<mysize; i++)
   {
      receive_disp_b[i] = border_size;
      border_size += receive_cnt_b[i];
   }
   //
   cpu_time_local += MPI_Wtime();
   if(quadratic)
   {
      eqsH->f_buffer[(int)eqsH->f_buffer.size()-3] = new double[border_size];
      eqsH->cpu_time += cpu_time_local;
   }
   else
   {
      eqs->f_buffer[(int)eqs->f_buffer.size()-3] = new double[border_size];
      eqs->cpu_time += cpu_time_local;
   }
#endif
   //
   // dim = n_loc*dof;
   // MPI_Allreduce(&dim,  &max_dimen, 1, MPI_INT,  MPI_MAX, MPI_COMM_WORLD);
}


/*************************************************************************
GeoSys-Function:
Task: for parallel solver
  HM monolithic case is to be considered.
Programming:
02/2008 WW Implementation
**************************************************************************/
double CPARDomain::Dot_Border_Vec(const double *vec_x, const double *vec_y)
{
   long i, l_buff;
   int ii, k;
   long b_shift[2];
   double val, fac;
   val = 0.;
   //
   for(k=0; k<nq; k++)
   {
      for(i=b_start[k];i<b_end[k];i++)
      {
         fac =  1.0/(double)bnode_connected_dom[i];
         for(ii=0; ii<dof; ii++)
         {
            l_buff = i+n_loc*ii+n_shift[k];
            val += fac*vec_x[l_buff]*vec_y[l_buff];
         }
      }
   }
   return val;
}


/*************************************************************************
GeoSys-Function:
Task: Parallel solver

  HM monolithic case is to be considered.
Programming:
07/2006 WW Implementation
12/2007 WW Revise
**************************************************************************/
double CPARDomain::Dot_Interior(const double *localr0,  const double *localr1)
{
   long i;
   int ii, k;
   double val;
   //

   /*
   if(dof>1)
   {
   //TEST

   string test = "rank";
   char stro[64];
   sprintf(stro, "%d",myrank);
   string test1 = test+(string)stro+"dom.txt";

   ofstream Dum(test1.c_str(), ios::out);
   Dum<<" nnodesHQ_dom  "<< nnodesHQ_dom<<endl;

   Dum<<" nq "<<nq <<endl;

   for(k=0; k<nq; k++)
   {
   Dum<<" i_start[k]  "<<i_start[k] <<endl;

   Dum<<"  i_end[k] "<<i_end[k] <<endl;

   for(i=i_start[k];i<i_end[k];i++)
   {
   for(ii=0; ii<dof; ii++)
   {
   //
   val += localr0[i+n_loc*ii]*localr0[i+n_loc*ii];

   Dum<<"[i+n_loc*ii] "<< i+n_loc*ii <<" localr0[i+n_loc*ii] "<< localr0[i+n_loc*ii]<<endl;

   }
   }
   }
   exit(1);

   }
   */

   val = 0.0;
   if(!localr1)
   {
      for(k=0; k<nq; k++)
      {
         for(i=i_start[k];i<i_end[k];i++)
         {
            for(ii=0; ii<dof; ii++)
               //
               val += localr0[i+n_loc*ii]*localr0[i+n_loc*ii];
         }
      }
   }
   else
   {
      for(k=0; k<nq; k++)
      {
         for(i=i_start[k];i<i_end[k];i++)
         {
            for(ii=0; ii<dof; ii++)
               val += localr0[i+n_loc*ii]*localr1[i+n_loc*ii];
         }
      }
   }

   /*
   //TEST
   Dum.close();
   if(nq>1)
    exit(1);
   */

   return val;
}


/*************************************************************************
GeoSys-Function:
Task: Parallel solver
    n: Dimension of the global EQS
  HM monolithic case is to be considered.
Programming:
06/2006 WW Implementation
12/2007 WW Revise
**************************************************************************/
void CPARDomain::Global2Local(const double *global_x, double *local_x, const long n )
{
   long i, ig;
   int ii;
   //
   //
   long n_global = (long)n/dof;
   for(i=0;i<n_loc*dof;i++)
      local_x[i] = 0.;
   for(i=0;i<n_loc;i++)
   {
      ig = nodes[i];
      for(ii=0; ii<dof; ii++)
         local_x[i+n_loc*ii] = global_x[ig+n_global*ii];
   }
}


/*************************************************************************
GeoSys-Function:
Task: Parallel solver
      n: Dimension of the global EQS
  HM monolithic case is to be considered.
Programming:
06/2006 WW Implementation
12/2007 WW Revise
02/2008 WW Revise
**************************************************************************/
void CPARDomain:: Local2Global(const double *local_x, double *global_x, const long n )
{
   long i, ig, b_index;
   int ii, k;
   double fac = 0.;
   //
   //
   long n_global = (long)n/dof;
   //
   for(i=0; i<n; i++)
      global_x[i] = 0.;
   //
   for(k=0; k<nq; k++)
   {
      for(i=i_start[k];i<i_end[k];i++)
      {
         ig = nodes[i];
         for(ii=0; ii<dof; ii++)
            global_x[ig+n_global*ii] = local_x[i+n_loc*ii];
      }
   }
   //
   for(k=0; k<nq; k++)
   {
      for(i=b_start[k];i<b_end[k];i++)
      {
         b_index = i+n_shift[k];
         fac =  1.0/(double)bnode_connected_dom[i];
         //
         ig = nodes[b_index];
         for(ii=0; ii<dof; ii++)
            global_x[ig+n_global*ii] += fac*local_x[b_index+n_loc*ii];
      }
   }
   //
}


/*************************************************************************
GeoSys-Function:
Task: Parallel solver

  HM monolithic case is to be considered.
Programming:
12/2007 WW
**************************************************************************/
void CPARDomain::Global2Border(const double *x, double *local_x, const long n )
{
   long i, ig;
   int ii, k;
   //
   long  nnodes_g = (long)n/dof;
   // BC
   for(i=0;i<dof*n_bc;i++)
      local_x[i] = 0.0;
   //
   for(i=0;i<n_bc;i++)
   {
      k = t_border_nodes[i];
      for(ii=0; ii<dof; ii++)
         local_x[i+n_bc*ii]= x[k+nnodes_g*ii];
   }
}


/*************************************************************************
GeoSys-Function:
Task: Parallel solver

  HM monolithic case is to be considered.
Programming:
12/2007 WW
**************************************************************************/
void CPARDomain::Border2Global(const double *local_x, double *x, const long n)
{
   long i, ig;
   int ii, k;
   //
   long nnodes_g = (long)n/dof;
   //
   for(i=0;i<n_bc;i++)
   {
      k = t_border_nodes[i];
      for(ii=0; ii<dof; ii++)
         x[k+nnodes_g*ii] = local_x[i+n_bc*ii];
   }
}


/*************************************************************************
GeoSys-Function:
Task: Parallel solver

  HM monolithic case is to be considered.
Programming:
07/2006 WW Implementation
12/2007 WW Revise
**************************************************************************/
void CPARDomain::Local2Border(const double *local_x, double *border_x)
{
   long i, ig;
   int ii, k;
   //
   // BC
   for(i=0;i<dof*n_bc;i++)
      border_x[i] = 0.0;
   //
   for(k=0; k<nq; k++)
   {
      for(i=b_start[k];i<b_end[k];i++)
      {
         ig = nodes_halo[i];
         for(ii=0; ii<dof; ii++)
            border_x[ig+n_bc*ii] = local_x[i+n_shift[k]+n_loc*ii];
      }
   }
}


/*************************************************************************
GeoSys-Function:
Task: Parallel solver

  HM monolithic case is to be considered.
Programming:
07/2006 WW Implementation
12/2007 WW Revise
**************************************************************************/
void CPARDomain::Border2Local(const double *border_x, double *local_x)
{
   long i, ig;
   int ii, k;
   //
   //
   for(k=0; k<nq; k++)
   {
      for(i=b_start[k];i<b_end[k];i++)
      {
         ig = nodes_halo[i];
         for(ii=0; ii<dof; ii++)
            local_x[i+n_shift[k]+n_loc*ii] = border_x[ig+n_bc*ii];
      }
   }
}


/*\!
 ********************************************************************
   Concatenate the inertanal entries of local subdomain solution
   Programm:
   12/2007 WW
********************************************************************/
//#define NEW_BREDUCE2
void CPARDomain::CatInnerX(double *global_x, const double *local_x, const long n)
{
   long i, j, ig;
   int ii, k;
   long counter = 0;
   double *x_i, *x_g;
   Linear_EQS *eq = NULL;
   if(quadratic) eq = eqsH;
   else          eq = eqs;
   //
   x_i = eq->f_buffer[0];
   x_g = eq->f_buffer[(long)eq->f_buffer.size()-1];
   //
   long n_global = (long)n/dof;
   //
#if defined(NEW_BREDUCE2)
   // Not finished
   // Due to the parallel computing of dom topology, not all num_inner_node of
   //   dom_vector[j] is caculated.
   // for(i=0; i<eq->A->Dim();i++)
   //   x_i[i] = 0.;
   //
   for(k=0; k<nq; k++)
   {
      for(i=i_start[k];i<i_end[k];i++)
      {
         for(ii=0; ii<dof; ii++)
         {
            x_i[counter] = local_x[i+n_loc*ii];
            counter++;                            //
         }
      }
   }
   // Concatentate
   MPI_Allgatherv (x_i, counter, MPI_DOUBLE, x_g, receive_cnt_i, receive_disp_i,
      MPI_DOUBLE, MPI_COMM_WORLD );
   //
   // Mapping to the golbal x
   CPARDomain * a_dom;
   for(j=0; j<mysize; j++)
   {
      counter = receive_disp_i[j];
      // Problem from here
      if(j==myrank)
         a_dom = this;
      else;
      a_dom = dom_vector[j];
      a_dom->nq = 1;
      a_dom->i_start[0] = 0;
      a_dom->i_end[0] =  a_dom->num_inner_nodes;  //Number of interior nodes
      if(quadratic)
      {
         a_dom->nq = 2;
         a_dom->i_start[1] = a_dom->i_end[0]+a_dom->num_boundary_nodes;
         a_dom->i_end[1] = a_dom->i_start[1]+a_dom->num_inner_nodesHQ;
      }

      for(k=0; k<a_dom->nq; k++)                  // This should come from different processors
      {
         for(i=a_dom->i_start[k];i<a_dom->i_end[k];i++)
         {
            ig = a_dom->nodes[i];
            for(ii=0; ii<dof; ii++)
            {
               global_x[ig+n_global*ii] = x_g[counter];
               counter++;
            }
         }
      }
   }

#else
   Local2Global(local_x, x_g, n);
   MPI_Allreduce( x_g, global_x, n, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
}


#if defined(NEW_BREDUCE)
/*\!
 ********************************************************************
   Reduce border entries by concatenating
   Programm:
   12/2007 WW
********************************************************************/
void CPARDomain::ReduceBorderV(double *local_x)
{
   long i, j, ig;
   int ii, k;
   long counter = 0;
   double *x_b, *x_cat, *x_g;
   Linear_EQS *eq = NULL;
   if(quadratic) eq = eqsH;
   else          eq = eqs;
   //
   x_g = eq->f_buffer[(long)eq->f_buffer.size()-1];
   x_cat = eq->f_buffer[(int)eq->f_buffer.size()-2];
   x_b = &local_x[eq->Dim()];
   //
   for(k=0; k<nq; k++)
   {
      for(i=b_start[k];i<b_end[k];i++)
      {
         for(ii=0; ii<dof; ii++)
         {
                                                  //;
            x_g[counter] = local_x[i+n_shift[k]+n_loc*ii];
            counter++;
         }
      }
   }
   // Concatentate
   MPI_Allgatherv (x_g, counter, MPI_DOUBLE, x_cat, receive_cnt_b, receive_disp_b,
      MPI_DOUBLE, MPI_COMM_WORLD );
   for(i=0;i<dof*n_bc;i++)
      x_b[i] = 0.0;
   //
   CPARDomain * a_dom;
   for(j=0; j<mysize; j++)
   {
      a_dom = dom_vector[j];
      counter = receive_disp_b[j];
      for(k=0; k<a_dom->nq; k++)
      {
         for(i=a_dom->b_start[k];i<a_dom->b_end[k];i++)
         {
            ig = a_dom->nodes_halo[i];
            for(ii=0; ii<dof; ii++)
            {
               x_b[ig+n_bc*ii] += x_cat[counter];
               counter++;
            }
         }
      }
   }
}
#endif                                            //if defined(NEW_BREDUCE)
/********************************************************************
   As the title
   Programm:
   12/2007 WW
********************************************************************/
void CPARDomain::PrintEQS_CPUtime(ostream &os)
{
   if(eqs)
   {
      os<<"CPU time elapsed in linear solver for linear elements: "
         <<eqs->GetCPUtime()<<endl;
   }
   if(eqsH)
   {
      os<<"CPU time elapsed in linear solver for quadratic elements: "
         <<eqsH->GetCPUtime()<<endl;
   }
}
#endif                                            //// if defined(USE_MPI)
