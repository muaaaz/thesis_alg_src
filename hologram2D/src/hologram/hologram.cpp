/***************************************************************************
 *   Copyright (C) 2005 by Siegfried Nijssen                               *
 *   snijssen@informatik.uni-freiburg.de                                   *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2010,2011 by Adriana Prado and Baptiste Jeudy           *
 *   baptiste.jeudy at univ-st-etienne fr                                  *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2015 by Romain Deville                                  *
 *   romain.deville[at]insa-lyon.fr                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 ***************************************************************************/

//================================ INCLUDE ===================================//
// Include library class
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stack>
// Include project class
#include "hologram.hpp"
#include "global.hpp"
#include "searchtreenode.hpp"

//=============================== NAMESPACE ==================================//
using namespace std;

//============================ STRUCT & TYPEDEF ==============================//


//=============================== VARIABLES ==================================//
// Will be used in cmpGToken when extending spatFirst or tempFirst
GToken tmpT1 = GToken();
GToken tmpT2 = GToken();


int de_arr[10000];
int sel_arr[10000];

int de = 0;
int roll_depth = 0;
//================================= CLASS ====================================//
//---- PUBLIC ----------------------------------------------------------------//
// Public CONSTANTS __________________________________________________________//
// Public Constructor/Desctructor ____________________________________________//
Hologram::Hologram():
  nbPatternAllClass(0),
  nbClosedPat(0),
  nbTotalClosedPat(0),
  freqPatternId(0),
  firstTick(0),
  totalTick(0),
  canonicalTick(0),
  extensionTick(0),
  subgraphisoTick(0),
  current_class_id(0),
  nodes_counter(0)
{
  /*
   * TODO : RD
   * Copy Desc
   */
  roll_depth = 0;
  vocabPattern      = new GVocab();
  v_NbPatternByClass.clear();
  v_ClassName.clear();
  v_nbClosedPatClass.clear();
}

Hologram::Hologram(float _minf, GClassDB* _pClassDB,vector<string> classes_names):
  nbPatternAllClass(0),
  nbClosedPat(0),
  nbTotalClosedPat(0),
  freqPatternId(0),
  firstTick(0),
  totalTick(0),
  canonicalTick(0),
  extensionTick(0),
  subgraphisoTick(0),
  current_class_id(1),
  v_ClassName(classes_names)
{
  /*
   * TODO : RD
   * Copy Desc
   */
  roll_depth = 0;
  minF = _minf;
  pClassDB = _pClassDB;
  vocabPattern      = new GVocab();
  v_NbPatternByClass.clear();
  v_ClassName.clear();
  v_nbClosedPatClass.clear();
}
Hologram::~Hologram()
{
  /*
   * TODO : RD
   * Copy Desc
   */
  delete vocabPattern;
  v_NbPatternByClass.clear();
  v_ClassName.clear();
  v_nbClosedPatClass.clear();
}

// Accessor __________________________________________________________________//
// Mutator ___________________________________________________________________//
// Public Methods ____________________________________________________________//
void Hologram::set_class_id(int id)
{
  current_class_id = id;
}

void Hologram::initialize( )
{
  number_of_classes =  pClassDB->number_of_classes;
  number_of_graphs =  pClassDB->number_of_graphs;
  class_count = pClassDB->class_count;

  v_NbPatternByClass.resize( number_of_classes + 1, 0 );
  v_nbClosedPatClass.resize( number_of_classes + 1, 0 );
  //for ( uint iClass = 0; iClass < v_GClassDB.size(); iClass++ )
  //  v_ClassName.push_back( v_GClassDB.at(iClass)->className );
}
// End of Grima::initialize( int nbClass )

//===========================================================================//
//                                   MAIN                                    //
//===========================================================================//
void print(vector<int> vec)
{
  for(int i=0;i< int(vec.size());i++)
    cout<<vec[i]<<' ';
  cout<<endl;
}

int Hologram::processMining( )
{
  /*
   * TODO : RD
   * Copy desc
   */

  if(pattern_buffer_set.size())
  {
    cerr<<"est not clear\n";
    exit(1);
  }
  
  memset(de_arr,0,sizeof de_arr);
  memset(sel_arr,0,sizeof sel_arr);
  
  cerr<<"start processMining\n";
  // Initialize pattern id to zero
  int minFreq      = 0;
  int returnStatus = 0;
  nbClosedPat      = 0;

  v_PatternCurrClass.clear();
  
  

  // cerr<<"number_of_classes: "<<number_of_classes<<endl;
  // cerr<<"number_of_graphs: "<<number_of_graphs<<endl;
  // cerr<<"class count: ";
  // print(class_count);
  
  /*************************************************************************/
  /*   REAL STUFF BEGIN HERE                                               */
  /*************************************************************************/
  // ============ ALL MINING STUFF BEGIN HERE !!!!!

  if ( minF <= 1 )
    minFreq = (int) round( minF * class_count[current_class_id] );
  else
    minFreq = minF;

  firstTick    = clock();
  returnStatus = search( pClassDB->v_ClassGraphs, pClassDB->m_TokenData, minFreq);
  totalTick    += clock() - firstTick;
 
  unbuffer_class_patterns();
  
  //  sort( vocabPattern->v_AllPatterns.begin(), vocabPattern->v_AllPatterns.end(), GpatComptLt() );

  for ( uint iPat = 0; iPat < uint(vocabPattern->v_AllPatterns.size()) ; iPat++ )
    vocabPattern->v_AllPatterns.at(iPat)->pGraph->graphID = iPat;

  return returnStatus;
}
// End of Grima::processMining( double minF , GClassDB *pClassDB )

void Hologram::saveData( )
{
  vocabPattern->saveVocab( PARAM.OUTDIR + PARAM.PATFILE, v_ClassName, v_ReturnStatus );
  ofstream patFile;
  patFile.open( PARAM.OUTDIR + PARAM.PAT_STAT_FILE );
  // for ( uint iClass = 0; iClass < v_ReturnStatus.size(); iClass++ )
  // {
  //   if ( v_ReturnStatus.at(iClass) == -1 || timeOutOverride )
  //     patFile << "# INCOMPLETE MINING !! TIMEOUT REACHED (" << PARAM.TIMEOUT
  //             << "H) for class " << v_ClassName.at(iClass) << endl;
  //   else if ( v_ReturnStatus.at(iClass) == -2  )
  //     patFile << "# INCOMPLETE MINING !! NB PAT LIMIT REACHED (" << PARAM.NBPATLIMIT
  //             << ") for class " << v_ClassName.at(iClass) << endl;
  // }
  patFile << "nb_total_pat," << vocabPattern->v_AllPatterns.size() << endl;
  patFile << "nb_total_closed_pat_class," << vocabPattern->v_AllPatterns.size() - nbTotalClosedPat << endl;
  patFile << "nb_pat_in_all_class," << nbPatternAllClass << endl;
  return ;
  for ( int iClass = 1; iClass <= number_of_classes; iClass++ )
  {
    patFile << "nb_pat_class_" << v_ClassName.at(iClass)
            << "," << v_NbPatternByClass.at(iClass) << endl;
    patFile << "nb_closed_pat_class_" << v_ClassName.at(iClass)
            << "," << v_NbPatternByClass.at(iClass) - v_nbClosedPatClass.at(iClass) << endl;
    patFile << "prop_closed_pat_class_" << v_ClassName.at(iClass)
            << "," << (double) v_nbClosedPatClass.at(iClass) / max(1,v_NbPatternByClass.at(iClass))  << endl;
  }
  patFile << "canonical_time_sec,"   << (double)canonicalTick / CLOCKS_PER_SEC << endl;
  patFile << "subgraphiso_time_sec," << (double)subgraphisoTick / CLOCKS_PER_SEC << endl;
  patFile << "extension_time_sec,"   << (double)extensionTick / CLOCKS_PER_SEC << endl;
  patFile << "total_time_sec,"       << (double)totalTick / CLOCKS_PER_SEC << endl;
  patFile.close();

  // for ( uint iClass = 0; iClass < v_ReturnStatus.size(); iClass++ )
  // {
  //   if ( v_ReturnStatus.at(iClass) == -1 || timeOutOverride )
  //     patFile << "# INCOMPLETE MINING !! TIMEOUT REACHED (" << PARAM.TIMEOUT
  //             << "H) for class " << v_ClassName.at(iClass) << endl;
  //   else if ( v_ReturnStatus.at(iClass) == -2  )
  //     patFile << "# INCOMPLETE MINING !! NB PAT LIMIT REACHED (" << PARAM.NBPATLIMIT
  //             << ") for class " << v_ClassName.at(iClass) << endl;
  // }
  cout  << "nb_total_pat," << vocabPattern->v_AllPatterns.size() << endl;
  cout  << "nb_total_closed_pat_class," << vocabPattern->v_AllPatterns.size() - nbTotalClosedPat << endl;
  cout  << "nb_pat_in_all_class," << nbPatternAllClass << endl;
  for ( int iClass = 1; iClass <= number_of_classes; iClass++ )
  {
    cout  << "nb_pat_class_" << v_ClassName.at(iClass)
          << "," << v_NbPatternByClass.at(iClass) << endl;
    cout  << "nb_closed_pat_class_" << v_ClassName.at(iClass)
          << "," << v_NbPatternByClass.at(iClass) - v_nbClosedPatClass.at(iClass) << endl;
    cout  << "prop_closed_pat_class_" << v_ClassName.at(iClass)
          << "," << (double) v_nbClosedPatClass.at(iClass) / max(1,v_NbPatternByClass.at(iClass)) << endl;
  }
  cout  << "canonical_time_sec,"   << (double)canonicalTick / CLOCKS_PER_SEC << endl;
  cout  << "subgraphiso_time_sec," << (double)subgraphisoTick / CLOCKS_PER_SEC << endl;
  cout  << "extension_time_sec,"   << (double)extensionTick / CLOCKS_PER_SEC << endl;
  cout  << "mapping_extension_time_sec,"   << (double)mappingExtTick / CLOCKS_PER_SEC << endl;
  cout  << "total_time_sec,"       << (double)totalTick / CLOCKS_PER_SEC << endl;
  cout << endl;
  cerr<<"end saving\n";
}

void Hologram::unbuffer_class_patterns()
{
  set<Pattern_buffer>::iterator it = pattern_buffer_set.begin();
  for( ; it != pattern_buffer_set.end() ; it++)
  {
    v_NbPatternByClass.at( current_class_id )++;
    if ( it->is_closed )
    {
      nbClosedPat++;
      nbTotalClosedPat++;
    }
    it->pat->pGraph->graphID = freqPatternId;
    freqPatternId++;
    v_PatternCurrClass.push_back(it->pat);
    string tmp_string = it->pat->getCanonincalString();
    unordered_map< string , bool >::iterator it1 =  all_patterns.find(tmp_string);
    if( it1 == all_patterns.end() )
    {
      all_patterns[tmp_string] = 1;
      vocabPattern->v_AllPatterns.push_back( it->pat );    
    }
    tmp_string.clear();
    tmp_string.shrink_to_fit();
  }
  pattern_buffer_set.clear();

  vocabPattern->v_PatternByClass.push_back( v_PatternCurrClass );
}

void print(GTokenData ob)
{
  for ( int i=0; i < int(ob.v_SparseOcc.size()) ; ++i)
  {
    
    if(!ob.v_SparseOcc.at(i).size())
      continue;
    cerr<<"G "<<ob.v_SparseOcc.at(i).graphID<<" C "<<ob.v_SparseOcc.at(i).pGraph->classID<<":\n";
    GSparseSet SparseOccTmp = ob.v_SparseOcc.at(i);
    uint domainSize = SparseOccTmp.size();

    for ( uint iDom = 0; iDom < domainSize  ; iDom++ )
    {
      GSparseSet::mapEdge mEdge = SparseOccTmp.at(iDom);
      cerr<<mEdge.nodeFrom<<' '<<mEdge.nodeDest<<endl;
    }
    
  }
}

void print( map<GToken, GTokenData, GTokenGt> &m_TokenData)
{

  map<GToken, GTokenData, GTokenGt>::iterator it = m_TokenData.begin();
  
  while(it != m_TokenData.end())
  {
    cerr<<"token: "<<it->first.nodeLabelFrom<<' '<<it->first.nodeLabelDest<<'\n';//<<it->second.v_SparseOcc.size()<<" :\n";
    print(it->second);
    it++;
  }
}

void print_report()
{
  cerr<<"de_arr; \n"; 
  for(int i=0;i<10000;++i)
    if(de_arr[i])
      cerr<<i<<' '<<de_arr[i]<<endl;
  cerr<<endl;
  
  cerr<<"sel_arr; \n";
  for(int i=0;i<10000;++i)
    if(sel_arr[i])
      cerr<<i<<' '<<sel_arr[i]<<endl;
  cerr<<endl;
}

int Hologram::search( vector<GGraph*>                   &v_Graphs,
                   map<GToken, GTokenData, GTokenGt> &m_TokenData,
                   GGlobFreq                         minFreq )
{
  /*
   * TODO : RD
   * Copy Desc
   */

  root = new searchtreenode();
  root->N_node = 1;
  root->nodeID = nodes_counter++;
  
  map<GToken, GTokenData, GTokenGt>::iterator it = m_TokenData.begin();
  cerr<<"reach the initial loop, minfreq = "<<minFreq<<", edges mapsize = "<<m_TokenData.size()<<" \n";

  while ( it != m_TokenData.end() )
  {
    // For all edges store in database
    if ( it->second.get_class_freq(current_class_id) >= int(minFreq) )
    {
      GExtensionData tmp;
      tmp.nbOcc     = 0;
      tmp.frequency = it->second.freq;
      
      for ( uint iGraph = 0 ; iGraph < it->second.v_SparseOcc.size() ; iGraph++ )
        tmp.nbOcc += it->second.v_SparseOcc.at(iGraph).size();
      
      root->valid_extenstions.push_back(make_pair(it->first,tmp));
    }
    // Go to next edge possible edge
    it++;
  }

  cerr<<"reach the main loop, root map size "<<root->valid_extenstions.size()<<endl;
  int budget = PARAM.TIMEOUT;
  cerr<<"budget: "<<budget<<endl;
  cerr<<"mining class number: "<<current_class_id<<endl;
  // the main loop
  while(true)
  {
    if( (clock() - firstTick ) / CLOCKS_PER_SEC >= (round( PARAM.TIMEOUT*60 / number_of_classes ) ) )
    {
      cerr<<"Time budget has finished\n";
      break;
    }

    // if(PARAM.NBPATLIMIT != -1 && number_of_mined_patterns > PARAM.NBPATLIMIT)
    // {
    //   cerr<<"Reached the limit in number of patterns\n";
    //   break;
    // }

    roll_depth = de = 0;
    
    delta = 0;
    N_delta = 0;

    // the father of the selected node:
    last_father = NULL;
    searchtreenode* selcted_node = select(root);
    
    if(selcted_node == NULL)
    {
     // cerr<<"bad sel\n";
      // if the selection lead to a node which iss fully expanded or 
      //has to possible expansions, we should delete it
      // if this node is the root, then we a comnplete search has been done
      if(last_father->parents.size() == 0 || last_father == root)
      {
        cerr<<"root is dead, root children numer: "<< root->children_nodes->size() <<"\n";
        root = NULL;
        //print_report();
        return 0;
      }

      // update the ancestors with dleta = 0, and delete 
      // the node from the search tree
      update_ancestors(last_father,delta);
 
      // generate the pattern of the deleted node so can remove it from the hash table
      GPattern* currentPattern = new GPattern();
      currentPattern->pGraph->graphID   = freqPatternId;
      currentPattern->pGraph->className = "FrequentPattern";
      build_pattern(last_father,currentPattern);
      string canonical_code_string = currentPattern->getCanonincalString();
      
      auto it = nodes_pointers.find(canonical_code_string);
      if(it != nodes_pointers.end())
        it->second = NULL;

      canonical_code_string.clear();
      canonical_code_string.shrink_to_fit();
      if(last_father->children_nodes->size())
      { 
        cerr<<"pad delete, children number: "<<last_father->children_nodes->size()<<"\n";
        exit(1);
      }
      delete_tree_node(last_father);
      delete currentPattern;
      budget++;
      continue;
    }

    //cerr<<"good sel\n";
    GToken ext;
    GExtensionData tmp_GExtensionData;
    // expand the node and save the edge that has been added
    searchtreenode* exp_node = expand(selcted_node,ext,tmp_GExtensionData);
    
    // if the node has no expansions, mark it as fully expanded, it could be deleted next iteration
    if(exp_node == NULL)
    {
      //cerr<<"bad expand\n";
      budget++;
      continue;
    }
    //cerr<<"good expand\n";
    // build the pattern
    GPattern* currentPattern = new GPattern();
    currentPattern->pGraph->graphID   = freqPatternId;
    currentPattern->pGraph->className = "FrequentPattern";
    build_pattern(selcted_node,currentPattern);

    // add the new node to the search tree
    selcted_node->children_nodes->insert(make_pair(ext,exp_node));
    
    // don't add the parent here, it has benn added in the constructor of the node
    // exp_node->parents.push_back(selcted_node);

    GExtensionData firstExtensionData;
    GTokenData tmp_GTokenData;

    // copy the GTokenData cuz we pass it as referance
    if(selcted_node == root)
      tmp_GTokenData = GTokenData (m_TokenData[ext]);
    else
      tmp_GTokenData = GTokenData (selcted_node->node_tokenData);
    
    roll_depth = 0;
    do_update = true;
    delta = 0;
    //cerr<<"start roll\n";
    // in the first level of the rollout, we compute the possible expansions of a pattern
    roll_out( exp_node,
              selcted_node,
              true,
              v_Graphs,
              minFreq,
              currentPattern,
              ext,
              tmp_GTokenData,
              tmp_GExtensionData,
              firstExtensionData);
    
    
    // average reward
    delta = delta / max(1,N_delta);

    // debugging info
    //if(roll_depth > 100)
    //  cerr<<"L: "<<roll_depth<<endl;

    if(do_update)
    {
      //cerr<<"L: "<<roll_depth<<endl;
      update_ancestors(exp_node,delta);
    }
    else
    {
      budget++;
    }
    // memory cleaning
    for(int i=0; i < int(tmp_GTokenData.v_SparseOcc.size()) ; ++i)
    {
      tmp_GTokenData.v_SparseOcc[i].clear();
    }
    tmp_GTokenData.v_SparseOcc.clear();
    tmp_GTokenData.v_SparseOcc.shrink_to_fit();
    delete currentPattern;
  }

  cerr<<"Budget has finished\n";
  // budget is finished
  //print_report();
  return 0;
}

// End of Grima::search( vector<GGraph*>                   &v_Graphs,
//                       map<GToken, GTokenData, GTokenGt> &m_TokenData,
//                       GGlobFreq                         minFreq )


inline double Hologram::UCB(searchtreenode* cur, searchtreenode* child)
{
  //cerr<<"child->N_node: "<<child->N_node<<endl;
  //return double(rand()%100) / double(100);
  double ret = child->Q + 2 * PARAM.C_p * sqrt( 2 * log(double(cur->N_node)) / double(child->N_node)  );
  return ret;// logaritem base??
}

searchtreenode* Hologram::best_child(searchtreenode* cur,GToken& ext)
{

    //this function will return the a pointer to the best child of the current MCTS node
    map< GToken ,searchtreenode*, GTokenGt>::iterator it;
    it = cur->children_nodes->begin();

    double cur_UCB = UCB(cur,it->second);
    double mx_UCB = cur_UCB;
    searchtreenode* ret = it->second;
    ext = it->first;
    
    it++;
    
    //cerr<<"children number: "<<cur->children_nodes->size()<<endl;;
    for(;it != cur->children_nodes->end() ;it++)
    {
      cur_UCB = UCB(cur,it->second);
      //cerr<<cur_UCB<<endl;
      if( cur_UCB > mx_UCB)
      {
        mx_UCB = cur_UCB;
        ret = it->second;
        ext = it->first;
      }

    }
   
    // ret could be NULL 
    return ret;
}

searchtreenode* Hologram::select(searchtreenode* cur)
{
  int mx_depth = int(1e6);
  int cc = 1;
  //cerr<<"start selcetion\n";
  while(mx_depth--){ // some stopping condition should we change this?
    //cerr<<"go deeper "<<cc<<" node address: "<<cur<<endl;
    if(!cur->is_fully_expanded)
    {
      //cerr<<"sel return "<<cur<<endl;
      return cur;
    }
    
    // if there is no children
    //cerr<<cur->children_nodes->size()<<endl;
    if(cur->children_nodes->size() == 0){
      last_father = cur;
      //cerr<<last_father<<" will be deleted\n";
      return NULL;
    }
    
    // debigging use only
    de++;
    ++sel_arr[cc];
    cc++;
    
    // save the edge added to get the best child
    GToken ext;
    last_father = cur;
    cur = best_child(cur,ext);
    if(cur == NULL)
    {
      //cerr<<last_father<<" will be deleted\n";
      return NULL;
    }
    //cerr<<"go deeper "<<cc<<" node address: "<<cur<<endl;
    // build the pattern
    //pPattern->push_back(ext,false);
  }

  // this sould never be excuted
  cerr<<"select depth not enough\n";
  exit(1);

  return cur; //should change this???
}

// we we want to expand a node we expect that teh node has it valid_extenstions ready but not all the children_nodes

searchtreenode* Hologram::expand(searchtreenode* cur,GToken& ext,GExtensionData& tmp)
{
  
  if(cur->valid_extenstions.size() == 0){
    cur->is_fully_expanded = true;
    cur->valid_extenstions.shrink_to_fit();
    return NULL;
  }

  // chose the last expansion, the vector is shuffled so it's a random choice
  ext = GToken(cur->valid_extenstions.back().first);
  tmp = GExtensionData(cur->valid_extenstions.back().second);
  
  //remove the selected choise from the list of available expansions
  cur->valid_extenstions.pop_back();
  
  if(cur->valid_extenstions.size() == 0)
  {
    cur->is_fully_expanded = true;
    cur->valid_extenstions.shrink_to_fit();
  }

  searchtreenode* new_child = new searchtreenode(cur,ext,nodes_counter++);
  //cerr<<"cereate node: "<<new_child<<endl;
  return new_child;
}

double Hologram::WRAcc(GTokenData& tokenData,int classID,int support)
{
  //return double(rand()%100) / double(100);

  double p_d = 0;

  for( int i = 0 ; i < int(tokenData.v_SparseOcc.size()) ; ++i)
  {
    if(tokenData.v_SparseOcc.at(i).pGraph->classID == classID 
        && tokenData.v_SparseOcc.at(i).size() > 0)
    {  
      p_d = p_d + 1.0;
    }
  }

  p_d /= double(max(1,support));
  
  double p = class_count[classID] / number_of_graphs;

  double ret = ( support * (p_d - p) ) / number_of_graphs;
  //double ret = p_d - p;

 // double ret = log( 1 + double(support) / double(number_of_graphs) ) * (p_d - p)   ;
  
  return ret;
  
}

void Hologram::add_parent(searchtreenode* old, searchtreenode* parent,const GToken& lastExt)
{

  if( parent->children_nodes->find(lastExt) != parent->children_nodes->end() )
  {
    parent->children_nodes->erase(parent->children_nodes->find(lastExt));
  }
  parent->children_nodes->insert(make_pair(lastExt,old));
  
  // wot?
  if(old->parents.back() != parent)
    old->parents.push_back(parent);
  
}

void Hologram::build_pattern(searchtreenode* selcted_node, GPattern* pPattern)
{
  //cerr<<"start build\n";
  stack<GToken> token_stack;
  while (selcted_node->parents.size() > 0)
  {
    token_stack.push(selcted_node->lastExt);
    selcted_node = selcted_node->parents[0];
  }
  //cerr<<token_stack.size()<<endl;
  while(!token_stack.empty())
  {
    pPattern->push_back(token_stack.top(),false);
    token_stack.pop();
  }
  //cerr<<"end build\n";
}

int Hologram::roll_out(searchtreenode* cur,
                        searchtreenode* parent,
                        bool rollout_first_level,
                        const vector<GGraph*>& v_Graphs,
                        const GGlobFreq       minFreq,    //Mininmum global frequency
                        GPattern*        pPattern,
                        const GToken&    lastExt, 
                        GTokenData      &tokenData,
                        GExtensionData  &suppData,   // Tmp variable, supposed frequency
                        GExtensionData  &prevData )
{
  //cerr<<"here0\n";
  N_delta++;
  GGlobFreq nbOcc        = 0;
  // current frequency of P
  GGlobFreq currentFreq  = 0;
  // Store last graph tID to compute frequency
  bool firstOcc          = true;
  uint lastOccGraphID    = 0;
  GTid lastOccGraphMemId = NULL;
  clock_t  firstTickTracker  = 0;

  // Add the extension to the pattern
  pPattern->push_back( lastExt, false );
  //cerr<<"here1\n";
  if(rollout_first_level)
  {
    firstTickTracker = clock();
    string canonical_code_string = pPattern->getCanonincalString();
    auto it = nodes_pointers.find(canonical_code_string);
    //cerr<<"here2\n";
    if(it == nodes_pointers.end())
    {
      //cerr<<"here3\n";
      //cerr<<"create new node "<<cur<<"\n";
      nodes_pointers[canonical_code_string] = cur;
      canonical_code_string.clear();
      canonical_code_string.shrink_to_fit();
    }
    else if (it->second == NULL)
    {
      //cerr<<"here4\n";
      if(cur->children_nodes->size())
      { 
        cerr<<"cur pad delete, children number: "<<cur->children_nodes->size()<<"\n";
      }
      delete_tree_node (cur);
      do_update = false;
      
      canonical_code_string.clear();
      canonical_code_string.shrink_to_fit();
      return 0;
    }
    else
    {
      //cerr<<"here5\n";
      //fine the connection to current node
      auto tmp_it = parent->children_nodes->find(lastExt);
      // delete it
      parent->children_nodes->erase(tmp_it);

      add_parent(it->second, parent, lastExt);
      
      if(cur->children_nodes->size())
      { 
        cerr<<"cur pad delete, children number: "<<cur->children_nodes->size()<<"\n";
      }

      delete_tree_node (cur);
      do_update = false;
      
      canonical_code_string.clear();
      canonical_code_string.shrink_to_fit();
      return 0;
    }
    canonicalTick += clock() - firstTickTracker;
  }
  //cerr<<"here6\n";
  //cerr<<"de: "<<de<<endl;
  //if(de > 20)
  //  return 0;
  // if this is the first time we go into this node "of it's not in the first level of the rool out"
  //if(!cur->occ_list_is_computed)
  //{
  // Object with map of possible extensions
  GExtensionCollect extCollect( minFreq, current_class_id );
  // Object that store pattern and occurences
  GSubgraphIso subGraphIso( pPattern, pPattern->v_OccList );
  // Clear occurence list before searching new occurency as pPattern should
  // already be stored in GVocab.
  subGraphIso.clearOccList();
  /* Number of sparse set (i.e. number of graphs in which first edge of pattern
  * is frequent */
  int my_last = -1;
  /// this is the most expensive part of the code!
  //cerr<<"-------\n";

  //bool del = 0;
  int nbGraphSparse = tokenData.v_SparseOcc.size();
  for ( int iGraph = 0 ; iGraph < nbGraphSparse; iGraph++ )
  {
     
    // For each sparseset
    // As the sparse set is modified, just create copy for the FOR LOOP
    
    GSparseSet SparseOccTmp = tokenData.v_SparseOcc.at(iGraph);
    
    // Get domain size
    for ( uint idx = 0; idx < uint(SparseOccTmp.size()) ; idx++ )
    {
      
      GSparseSet::mapEdge mEdge = SparseOccTmp.at(idx);
      // Find if this edge wear an occurence of the pattern
      firstTickTracker = clock();

      
      bool findOcc = subGraphIso.run( SparseOccTmp.pGraph,
                                      mEdge.nodeFrom,
                                      mEdge.edgeId   );
      subgraphisoTick += clock() - firstTickTracker;

      if ( findOcc )
      {
        
        // If edge wears pattern
        // Increment occurences counter
        nbOcc++;
        if ( firstOcc )
        {
          firstOcc = false;
          // If first occurenc of pattern in graph
          // Increment frequency counter to 1
          currentFreq++;
          // Get graph ID of this occurence
          lastOccGraphID    = SparseOccTmp.graphID;
          lastOccGraphMemId = SparseOccTmp.pGraph;

        }
        else if ( lastOccGraphID != SparseOccTmp.graphID
                  && lastOccGraphMemId == SparseOccTmp.pGraph )
        {
          cerr << "ERROR - Error While managing graphID in sparseset" << endl;
        }
        else if ( lastOccGraphID != SparseOccTmp.graphID
                  && lastOccGraphMemId != SparseOccTmp.pGraph )
        {
          // If last occurence was in another graph
          // Increment frequency counter
          currentFreq++;
          // Get graph ID of this occurence
          lastOccGraphID    = SparseOccTmp.graphID;
          lastOccGraphMemId = SparseOccTmp.pGraph;
        }
        else if(iGraph != my_last)
          currentFreq++;
        // Compute all extension of this occurence
        firstTickTracker = clock();
        subGraphIso.current_occurrence_class_id = SparseOccTmp.pGraph->classID;
        extCollect.process( subGraphIso );
        extensionTick += clock() - firstTickTracker;
        my_last = iGraph;
      }
      else
      {
        // If edge does not wear pattern "Remove" edge from domain in the sparse set
        tokenData.v_SparseOcc.at(iGraph).remove(mEdge) ;
        //if(!(de%100))
       // cerr<<"*";
        //del = 1;
      }
    }

    tokenData.v_SparseOcc.at(iGraph).data.shrink_to_fit();

  }
 

  for ( int iGraph = 0 ; iGraph < int(tokenData.v_SparseOcc.size()) ; iGraph++ )
  {
    if (tokenData.v_SparseOcc[iGraph].size() == 0)
    {
      tokenData.v_SparseOcc[iGraph] = tokenData.v_SparseOcc.back();
      tokenData.v_SparseOcc.pop_back();
      --iGraph;
    }
  }

  tokenData.v_SparseOcc.shrink_to_fit();
  //cerr<<"here7\n";
  //if(del)
  //  cerr<<"\n";
  if(rollout_first_level)
  {
  // save the occ list in the node
    cur->node_tokenData = GTokenData(tokenData);
    cur->occ_list_is_computed = true;
  }

  //here tokenData contains the occlist of the current patarn
  // and extCollect contains all the possible expansions
  /////////////////////////////////////////////////////

  mappingExtTick += extCollect.mapExtTick;

  //cerr<<"L: "<<de<<endl; 
  

  if ( currentFreq != suppData.frequency )
  {
   
    cerr<<"L: "<<de<<" is first level: "<<rollout_first_level<<endl; 
    cerr << "Computed Frequency not the same as supposed one" << endl;
    cerr << "# Suppose freq : " << suppData.frequency << " VS "<<currentFreq <<" VS "<<tokenData.size()<< endl;
    cerr << "# Suppose occ  : " << suppData.nbOcc     << endl;
    cerr << pPattern;
    cerr << endl;
    print(tokenData);
    exit( EXIT_FAILURE ); 
  }
  if ( nbOcc != suppData.nbOcc && pPattern->v_Tokens.at(0).angle > 0  )
  {
    //return 0;
    cerr<<"L: "<<de<<endl; 
    cerr << "Computed occurency not the same as supposed one" << endl;
    cerr << "# Suppose freq : " << suppData.frequency << endl;
    cerr << "# Suppose occ  : " << suppData.nbOcc     << endl;
    uint i = 0;
    while ( i < suppData.v_OccList.size() )
    {
      cerr << "o "
          << suppData.v_OccList.at(i) << " "
          << suppData.v_OccList.at(i+1) << " "
          << suppData.v_OccList.at(i+2) << " " << endl;
      i+=3;
    }
    cerr << pPattern;
    cerr << endl;
    //return 0;
    exit( EXIT_FAILURE );
  }
  
  /*************************************************************************/
  
  double cur_Q = WRAcc(tokenData,current_class_id,currentFreq);
  delta += cur_Q;
  //cerr<<"here8\n";
  // save hthe pattern
  //cerr<<"start add, nplimit: "<<PARAM.NBPATLIMIT<<"\n";
  bool add_pattern = false;
  if( int(pattern_buffer_set.size()) < PARAM.NBPATLIMIT  )
  {
    add_pattern = true;
  }
  else if(pattern_buffer_set.begin()->evaluation < cur_Q)
  {
    add_pattern = true;
  }

  if(add_pattern)
  {
    string tmp_string = pPattern->getCanonincalString();
    unordered_map< string , bool >::iterator it = class_patterns.find(tmp_string);
    if(it == class_patterns.end())
    {
      class_patterns[tmp_string] = 1;

      //not found
      GPattern * tmpPat = new GPattern( pPattern );
      bool is_closed = prevData.frequency == suppData.frequency 
                        && prevData.nbOcc == suppData.nbOcc ;

      pattern_buffer_set.insert( Pattern_buffer(tmpPat,is_closed,cur_Q) );
      
      if( int(pattern_buffer_set.size()) > PARAM.NBPATLIMIT )
      {
        tmpPat = pattern_buffer_set.begin()->pat;
        pattern_buffer_set.erase(pattern_buffer_set.begin());
        delete tmpPat;
      }

    }

    tmp_string.clear();
    tmp_string.shrink_to_fit();
  }
  //cerr<<"end add\n";
  //cerr<<"here13\n";
  /********************* store the possible expansions ********************/
  vector<pair<GToken, GExtensionData> > tmp_vec;

  for ( map<GToken, GExtensionData, GTokenGt>::iterator it= extCollect.m_Extensions.begin();
        it != extCollect.m_Extensions.end(); it++ )
  {
    if ( it->second.current_class_frequency >= int(minFreq) )
    {
      // save the possible moves on the mape
      if ( it->first.angle != GNOANGLE )
      {
        tmp_vec.push_back(make_pair(it->first,it->second));
      }
    }
    
  }

  //no need for this
  extCollect.m_Extensions.clear();

  // we shuffle to make the chose of the action in the simulation random. Good?
  random_shuffle ( tmp_vec.begin(), tmp_vec.end() );
  
  if(rollout_first_level)
  {
    cur->valid_extenstions = tmp_vec;
  }

  //cerr<<lastExt;
  //if(!(de%100) || de < 10)
  if( 0 && de > 100 )
    cerr<<"L: "<<de<<" B: "<< tmp_vec.size()
    <<" Memory: "<<tokenData.size()<<" Frequancy: "<<currentFreq
    <<" Occ: "<<nbOcc<<endl;
  //cerr.flush();
  roll_depth = max(roll_depth,de);
  if(0)
  {
    cerr<<pPattern;
    cerr<<"possible expansions:\n";
    for(int i=0; i < int(tmp_vec.size()) ; ++i )
      cerr<<tmp_vec[i].first;
    cerr<<"--------------------\n";
    print(tokenData);
  }
  //cerr.flush();
  
  //}
  // if it's not the first time we visite this node: restore the occ list from the node
  //else{
  //  cerr<<"get the saved data\n";
  //  tokenData = GTokenData(cur->node_tokenData);
  //  cerr<<"valid expansions size: "<<cur->valid_extenstions.size()<<endl;
 // }
  

  if(tmp_vec.size() == 0 )
  {
    //cerr<<"No valid expantions!\n";
    pPattern->pop_back( false );
    ++de_arr[roll_depth];
    roll_depth--;
    de--;
    return 0;
    //wot?
  }
  ++roll_depth;
  ++de;

  // chose the last element
  pair<GToken,GExtensionData> sel_move = tmp_vec.back();

  int ret;
 
  ret = roll_out(NULL,
                NULL,
                false,
                v_Graphs, 
                minFreq, 
                pPattern,
                sel_move.first,
                tokenData, 
                sel_move.second,
                suppData);
  
  de--;
  roll_depth--;
  return ret;
}

void Hologram::update_ancestors(searchtreenode* cur, double _delta)
{ 
  if(cur == NULL )
    return;
  
  cur->Q = (cur->N_node*cur->Q+_delta)/max(1,cur->N_node+1);
  if(cur->Q < -1e9)
    cerr<<"very low quality2\n";
  cur->N_node = cur->N_node+1;
  
  if (cur->parents.size() == 0)
  {
    // the root
    return;
  }
  for(auto it : cur->parents)
  {
    update_ancestors(it,_delta);
  }
}

void Hologram::delete_tree_node(searchtreenode* cur)
{
  // delete connections from paretns
  //cerr<<"start deleting "<<cur<<"\n";

  for(auto dad : cur->parents)
  {
    vector< map<GToken, searchtreenode*, GTokenGt>::iterator > it_vec;
    for(auto it = dad->children_nodes->begin() ; 
        it != dad->children_nodes->end() ;
        it++)
    {
      if(it->second == cur)
      {
        //dad->children_nodes->erase(it);
        it_vec.push_back(it);
        //break;
      }
    }
    for(auto it : it_vec)
    {
       dad->children_nodes->erase(it);
    }
    it_vec.clear();
    it_vec.shrink_to_fit();
  }

  // delete connections from sons "no need"

  // we will never delete a node with a son

  for(auto it = cur->children_nodes->begin(); 
    it != cur->children_nodes->end() ; 
    it++ )
  {
    for(int i=0; i < int(it->second->parents.size()) ; i++)
    {
      if(it->second->parents[i] == cur)
      {
        //cerr<<"DELDEL\n";
        it->second->parents[i] = it->second->parents.back();
        it->second->parents.pop_back();
      }
    }
  }
  //cerr<<"deleting\n";
  
  // delete node
  delete cur;
  //cerr<<"end deleting\n";  
}

void Hologram::clean()
{
  //cerr<<"here00\n";
  unordered_map< string , searchtreenode* >::iterator it = nodes_pointers.begin();
  //cerr<<"here0\n";
  for(it = nodes_pointers.begin();
      it != nodes_pointers.end(); it++)
  {
    if(it->second != NULL)
    {
      delete_tree_node (it->second);
    }
  }
  
  class_patterns.clear();
  nodes_pointers.clear();
  pattern_buffer_set.clear();

  return;
 
}

void Hologram::delete_search_subtree(searchtreenode* cur)
{

  map< GToken ,searchtreenode*, GTokenGt>::iterator it;

  for(it = cur->children_nodes->begin();
      it != cur->children_nodes->end() ;it++)
  {
    delete_search_subtree(it->second);
  }

  search_tree_nodes[cur->nodeID] = cur;

}
