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
// Include project class
#include "MCTSgrima.hpp"
#include "global.hpp"
#include "MCTS.hpp"

//=============================== NAMESPACE ==================================//
using namespace std;

//============================ STRUCT & TYPEDEF ==============================//

//=============================== VARIABLES ==================================//
// Will be used in cmpGToken when extending spatFirst or tempFirst
GToken tmpT1 = GToken();
GToken tmpT2 = GToken();

//================================= CLASS ====================================//
//---- PUBLIC ----------------------------------------------------------------//
// Public CONSTANTS __________________________________________________________//
// Public Constructor/Desctructor ____________________________________________//
MCTSGrima::MCTSGrima():
  nbPatternAllClass(0),
  nbClosedPat(0),
  nbTotalClosedPat(0),
  freqPatternId(0),
  firstTick(0),
  totalTick(0),
  canonicalTick(0),
  extensionTick(0),
  subgraphisoTick(0)
{
  /*
   * TODO : RD
   * Copy Desc
   */
  vocabPattern      = new GVocab();
  v_NbPatternByClass.clear();
  v_ClassName.clear();
  v_nbClosedPatClass.clear();
}

MCTSGrima::MCTSGrima(float _minf,GClassDB* _pClassDB):
  nbPatternAllClass(0),
  nbClosedPat(0),
  nbTotalClosedPat(0),
  freqPatternId(0),
  firstTick(0),
  totalTick(0),
  canonicalTick(0),
  extensionTick(0),
  subgraphisoTick(0)
{
  /*
   * TODO : RD
   * Copy Desc
   */
  minF = _minf;
  pClassDB = _pClassDB;
  vocabPattern      = new GVocab();
  v_NbPatternByClass.clear();
  v_ClassName.clear();
  v_nbClosedPatClass.clear();
}
MCTSGrima::~MCTSGrima()
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
void MCTSGrima::initNbPatternByClass( vector<GClassDB*> v_GClassDB )
{
  v_NbPatternByClass.resize( v_GClassDB.size(), 0 );
  v_nbClosedPatClass.resize( v_GClassDB.size(), 0 );
  for ( uint iClass = 0; iClass < v_GClassDB.size(); iClass++ )
    v_ClassName.push_back( v_GClassDB.at(iClass)->className );
}
// End of Grima::initNbPatternByClass( int nbClass )

//===========================================================================//
//                                   MAIN                                    //
//===========================================================================//
int MCTSGrima::processMining( )
{
  /*
   * TODO : RD
   * Copy desc
   */
  
  // Initialize pattern id to zero
  int minFreq      = 0;
  int returnStatus = 0;
  nbClosedPat  = 0;

  v_PatternCurrClass.clear();

  /*************************************************************************/
  /*   REAL STUFF BEGIN HERE                                               */
  /*************************************************************************/
  // ============ ALL MINING STUFF BEGIN HERE !!!!!
  // Begin of recursive cal
  //  int nbPattern = 0;
  //  uint lastSizePattern = vocabPattern->v_Patterns.size() ;
  if ( minF <= 1 )
    minFreq = (int) round( minF * pClassDB->v_ClassGraphs.size() );
  else
    minFreq = minF;

  firstTick    = clock();
  returnStatus = search( pClassDB->v_ClassGraphs, pClassDB->m_TokenData, minFreq);
  totalTick   += clock() - firstTick;

  vocabPattern->v_PatternByClass.push_back( v_PatternCurrClass );

  //  sort( vocabPattern->v_AllPatterns.begin(), vocabPattern->v_AllPatterns.end(), GpatComptLt() );

  for ( uint iPat = 0; iPat < vocabPattern->v_AllPatterns.size(); iPat++ )
    vocabPattern->v_AllPatterns.at(iPat)->pGraph->graphID = iPat;

  return returnStatus;
}
// End of Grima::processMining( double minF , GClassDB *pClassDB )

void MCTSGrima::saveData( bool timeOutOverride )
{
  vocabPattern->saveVocab( PARAM.OUTDIR + PARAM.PATFILE, v_ClassName, v_ReturnStatus );
  ofstream patFile;
  patFile.open( PARAM.OUTDIR + PARAM.PAT_STAT_FILE );
  for ( uint iClass = 0; iClass < v_ReturnStatus.size(); iClass++ )
  {
    if ( v_ReturnStatus.at(iClass) == -1 || timeOutOverride )
      patFile << "# INCOMPLETE MINING !! TIMEOUT REACHED (" << PARAM.TIMEOUT
              << "H) for class " << v_ClassName.at(iClass) << endl;
    else if ( v_ReturnStatus.at(iClass) == -2  )
      patFile << "# INCOMPLETE MINING !! NB PAT LIMIT REACHED (" << PARAM.NBPATLIMIT
              << ") for class " << v_ClassName.at(iClass) << endl;
  }
  patFile << "nb_total_pat," << vocabPattern->v_AllPatterns.size() << endl;
  patFile << "nb_total_closed_pat_class," << vocabPattern->v_AllPatterns.size() - nbTotalClosedPat << endl;
  patFile << "nb_pat_in_all_class," << nbPatternAllClass << endl;
  for ( uint iClass = 0; iClass < v_NbPatternByClass.size(); iClass++ )
  {
    patFile << "nb_pat_class_" << v_ClassName.at(iClass)
            << "," << v_NbPatternByClass.at(iClass) << endl;
    patFile << "nb_closed_pat_class_" << v_ClassName.at(iClass)
            << "," << v_NbPatternByClass.at(iClass) - v_nbClosedPatClass.at(iClass) << endl;
    patFile << "prop_closed_pat_class_" << v_ClassName.at(iClass)
            << "," << (double) v_nbClosedPatClass.at(iClass) / v_NbPatternByClass.at(iClass)  << endl;
  }
  patFile << "canonical_time_sec,"   << (double)canonicalTick / CLOCKS_PER_SEC << endl;
  patFile << "subgraphiso_time_sec," << (double)subgraphisoTick / CLOCKS_PER_SEC << endl;
  patFile << "extension_time_sec,"   << (double)extensionTick / CLOCKS_PER_SEC << endl;
  patFile << "total_time_sec,"       << (double)totalTick / CLOCKS_PER_SEC << endl;
  patFile.close();

  for ( uint iClass = 0; iClass < v_ReturnStatus.size(); iClass++ )
  {
    if ( v_ReturnStatus.at(iClass) == -1 || timeOutOverride )
      patFile << "# INCOMPLETE MINING !! TIMEOUT REACHED (" << PARAM.TIMEOUT
              << "H) for class " << v_ClassName.at(iClass) << endl;
    else if ( v_ReturnStatus.at(iClass) == -2  )
      patFile << "# INCOMPLETE MINING !! NB PAT LIMIT REACHED (" << PARAM.NBPATLIMIT
              << ") for class " << v_ClassName.at(iClass) << endl;
  }
  cout  << "nb_total_pat," << vocabPattern->v_AllPatterns.size() << endl;
  cout  << "nb_total_closed_pat_class," << vocabPattern->v_AllPatterns.size() - nbTotalClosedPat << endl;
  cout  << "nb_pat_in_all_class," << nbPatternAllClass << endl;
  for ( uint iClass = 0; iClass < v_NbPatternByClass.size(); iClass++ )
  {
    cout  << "nb_pat_class_" << v_ClassName.at(iClass)
          << "," << v_NbPatternByClass.at(iClass) << endl;
    cout  << "nb_closed_pat_class_" << v_ClassName.at(iClass)
          << "," << v_NbPatternByClass.at(iClass) - v_nbClosedPatClass.at(iClass) << endl;
    cout  << "prop_closed_pat_class_" << v_ClassName.at(iClass)
          << "," << (double) v_nbClosedPatClass.at(iClass) / v_NbPatternByClass.at(iClass)  << endl;
  }
  cout  << "canonical_time_sec,"   << (double)canonicalTick / CLOCKS_PER_SEC << endl;
  cout  << "subgraphiso_time_sec," << (double)subgraphisoTick / CLOCKS_PER_SEC << endl;
  cout  << "extension_time_sec,"   << (double)extensionTick / CLOCKS_PER_SEC << endl;
  cout  << "mapping_extension_time_sec,"   << (double)mappingExtTick / CLOCKS_PER_SEC << endl;
  cout  << "total_time_sec,"       << (double)totalTick / CLOCKS_PER_SEC << endl;
  cout << endl;
}
// End of Grima::saveData( string outDir, int returnStatus )

//---- PROTECTED  ------------------------------------------------------------//
// Protected CONSTANTS _______________________________________________________//
// Protected Methods _________________________________________________________//

//---- PRIVATE ---------------------------------------------------------------//
// Private CONSTANTS _________________________________________________________//
// Private Methods ___________________________________________________________//
int MCTSGrima::search( vector<GGraph*>                   &v_Graphs,
                   map<GToken, GTokenData, GTokenGt> &m_TokenData,
                   GGlobFreq                         minFreq )
{
  /*
   * TODO : RD
   * Copy Desc
   */

  MCTS_node* root = new MCTS_node();

  int returnStatus = 0;
  // Instanciate iterator
  map<GToken, GTokenData, GTokenGt>::iterator it = m_TokenData.begin();
  // Instanciate current pattern object
  GPattern currentPattern;
  currentPattern.pGraph->graphID   = freqPatternId;
  currentPattern.pGraph->className = "FrequentPattern";
  while ( it != m_TokenData.end() )
  {
    // For all canonical code of 1-edge store in database
    if ( it->second.freq >= minFreq )
    {
      GExtensionData tmp; 
      GExtensionData first;
      tmp.nbOcc     = 0;
      tmp.frequency = it->second.freq;
      for ( uint iGraph = 0 ; iGraph < it->second.v_SparseOcc.size() ; iGraph++ )
        tmp.nbOcc += it->second.v_SparseOcc.at(iGraph).size;
      root->valid_extenstions->insert(make_pair(it->first,tmp));
      root->children_nodes->insert(make_pair(it->first,(MCTS_node*)(NULL)));

      // Apply recursive call
      //returnStatus = search( v_Graphs, minFreq, &currentPattern,
       //                      it->first, it->second, tmp, first, iClassDB );

      if ( returnStatus != 0 )
        return returnStatus;
    }
    // Go to next edge possible edge
    it++;
  }

  int budget = 1994;
  while(budget--){
    GPattern currentPattern;
    currentPattern.pGraph->graphID   = freqPatternId;
    currentPattern.pGraph->className = "FrequentPattern";
    MCTS_node* sel_node = select(root);
    MCTS_node* exp_node = expand(sel_node);

    ///we have to make the expand function save the action that cuase the expansion. wot??
    currentPattern.push_back(/*the move we did in expand*/) // wot?
    

    double delta = rool_out(exp_node,true,v_Graphs,minF,&currentPattern,
                            );
    update_ancestors(exp_node,delta);
  }

  return returnStatus;
}

// End of Grima::search( vector<GGraph*>                   &v_Graphs,
//                       map<GToken, GTokenData, GTokenGt> &m_TokenData,
//                       GGlobFreq                         minFreq )


inline double MCTSGrima::UCB(MCTS_node* cur, MCTS_node* child){
  return child->Q + 2*PARAM.C_p*sqrt(2*(log(cur->N_node))/child->N_node); // logaritem base??
}

MCTS_node* MCTSGrima::best_child(MCTS_node* cur){
    //this functio will return the best child of the current MCTS node
    MCTS_node* ret = NULL;
    double cur_UCB = -1, mx_UCB = -1 ;
    map< GToken ,MCTS_node*, GTokenGt>::iterator it;
    for(it = cur->children_nodes->begin();it != cur->children_nodes->end() ;it++)
    {
      cur_UCB = UCB(cur,it->second);
      if( cur_UCB > mx_UCB)
        mx_UCB = cur_UCB, ret = it->second;
    }
    return ret;
}

MCTS_node* MCTSGrima::select(MCTS_node* cur){
  
  int mx_depth = 200;
  while(mx_depth--){ // some stopping condition
    if(!cur->is_fully_expanded)
      return cur;
    cur = best_child(cur);
  }

  return cur; //should change this???
}

MCTS_node* MCTSGrima::expand(MCTS_node* cur){
  if(cur->is_fully_expanded){
    //do soething that make the algorithm dont come here again
  }
  MCTS_node* new_child;
  // we should have the list of the possible extentions here
  // choose as the algorithm says

  map< GToken ,MCTS_node*, GTokenGt>::iterator it = cur->children_nodes->begin();
  while(true){
    if(it == cur->children_nodes->end()){
      cur->is_fully_expanded = true;
      break;
    }
    if(it->second == NULL)
      

    it++;
  }
  if(cur->is_fully_expanded){
    //do soething that make the algorithm dont come here again
  }

  return new_child;
}

double MCTSGrima::rool_out(MCTS_node* cur, 
                    bool first_level,
                    vector<GGraph*>* v_Graphs,
                    GGlobFreq       &minFreq,    //Mininmum global frequency
                    GPattern*        pPattern,
                    GTokenData      &tokenData,
                    GExtensionData  &suppData,   // Tmp variable, supposed frequency
                    GExtensionData  &prevData )
{

  // pPattern the expanded pattern
  // tokenData is the occlist of the parent pattern

  // calculate the possible valid extention and but them in a map

  // Object with map of possible extensions
  GExtensionCollect extCollect( minFreq );
  // Object that store pattern and occurences
  GSubgraphIso subGraphIso( pPattern, pPattern->v_OccList );
  // Clear occurence list before searching new occurency as pPattern should
  // already be stored in GVocab.
  subGraphIso.clearOccList();
  /* Number of sparse set (i.e. number of graphs in which first edge of pattern
   * is frequent */
  
  //calculate the feruancy of a pattern and the possible expansion at the same time
  // Nb of occurence of Pattern ( To check if P is close)
  GGlobFreq nbOcc        = 0;
  // current frequency of P
  GGlobFreq currentFreq  = 0;
  // Store last graph tID to compute frequency
  bool firstOcc          = true;
  uint lastOccGraphID    = 0;
  GTid lastOccGraphMemId = NULL;
  clock_t  firstTickTracker  = 0;
  int nbGraphSparse = tokenData.v_SparseOcc.size();
  for ( int iGraph = 0 ; iGraph < nbGraphSparse; iGraph++ )
  {
    // For each sparseset
    // As the sparse set is modified, just create copy for the FOR LOOP
    GSparseSet SparseOccTmp = tokenData.v_SparseOcc.at(iGraph);
    // Get domain size
    uint domainSize = SparseOccTmp.size;
    for ( uint iDom = 0; iDom < domainSize  ; iDom++ )
    {
      // For each edge in the domain of the sparse set
      // Get map value of edge wearing pattern
      uint map                  = SparseOccTmp.atDom(iDom);
      GSparseSet::mapEdge mEdge = SparseOccTmp.atMap(map);
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
        // Compute all extension of this occurence
        firstTickTracker = clock();
        extCollect.process( subGraphIso );
        extensionTick += clock() - firstTickTracker;
      }
      else
      {
        // If edge does not wear pattern "Remove" edge from domain in the sparse set
        tokenData.v_SparseOcc.at(iGraph).remove(map) ;
      }
    }
  }


  //here tokenData contains the occlist of the current patarn
  // and extCollect contains all the possible expansions
  /////////////////////////////////////////////////////

  mappingExtTick += extCollect.mapExtTick;

  if ( currentFreq != suppData.frequency )
  {
    cerr << "Computed Frequency not the same as supposed one" << endl;
    cerr << "# Suppose freq : " << suppData.frequency << endl;
    cerr << "# Suppose occ  : " << suppData.nbOcc     << endl;
    cerr << pPattern;
    cerr << endl;
    exit( EXIT_FAILURE );
  }
  if ( nbOcc != suppData.nbOcc && pPattern->v_Tokens.at(0).angle > 0  )
  {
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
    exit( EXIT_FAILURE );
  }

  /*************************************************************************/
  {
    // === Fourth step : Output pattern and occurences in output file
    // Increment pattern Id
    uint countPat = 0;
    bool find     = false;
    while ( countPat < vocabPattern->v_AllPatterns.size() && !find )
    {
      find = vocabPattern->comp( pPattern, vocabPattern->v_AllPatterns[countPat] );
      if ( !find )
        countPat++;
    }

    GPattern * tmpPat = new GPattern( *pPattern );
    tmpPat->pGraph = new GGraph( *pPattern->pGraph );

    if ( !find )
    {
      //v_NbPatternByClass.at( currClassIdx )++;

      if ( prevData.frequency == suppData.frequency
          && prevData.nbOcc == suppData.nbOcc )
      {
        nbClosedPat++;
        nbTotalClosedPat++;
      }
      pPattern->pGraph->graphID = freqPatternId;
      freqPatternId++;
      vocabPattern->v_AllPatterns.push_back( tmpPat );
    }
    v_PatternCurrClass.push_back(tmpPat);

    if ( PARAM.DEBUG_MODE )
    {
      // Write pattern in output file
      cout << "# Suppose freq : " << suppData.frequency << endl;
      cout << "# Suppose occ  : " << suppData.nbOcc     << endl;
      cout << "p " << freqPatternId << endl;
      cout << pPattern->v_Tokens << endl;
    }
  }

  /********************* Make recurive call to extend P ********************/

  cur->valid_extenstions->clear();

  for ( map<GToken, GExtensionData, GTokenGt>::iterator it= extCollect.m_Extensions.begin();it != extCollect.m_Extensions.end(); it++ )
  {
    if ( it->second.frequency >= minFreq )
    {
      
      // Recursive call
      if ( it->first.angle != GNOANGLE )
      {
        //isLeaf = search( v_Graphs, minFreq, pPattern,
        //                 it->first, tokenData, it->second, suppData,iClassDB );
        
        GToken lastExt = it->first;
         //possible expantion loop
      
        // Add the extension to the pattern
        pPattern->push_back( lastExt, false );

        // === Second step : Check if P is canonical
        firstTickTracker = clock();
        //  cout << pPattern->v_Tokens << endl;

        bool isCanonical = pPattern->isCanonincal();

        canonicalTick += clock() - firstTickTracker;

        //check some stuff
        {
          //----------------------------------------------------------------------------
          // CANONICAL TEST STOP
          if ( !isCanonical )
          {
            // If extension is not canonical, we should not go deeper
            pPattern->pop_back( false );
            // Exit search. We will find this pattern later
            continue;
          }
          //----------------------------------------------------------------------------
          // PATTERN NB TIME NODE LIMIT
          if ( PARAM.TIMELIMIT != -1
              && ( pPattern->maxTCoord - pPattern->minTCoord > PARAM.TIMELIMIT - 1 ) )
          {
            // If extension is not canonical, we should not go deeper
            pPattern->pop_back( false );
            // Exit search. We will find this pattern later
            continue;
          }
          //----------------------------------------------------------------------------
          // PATTERN NB NODE LIMIT STOP
          if ( PARAM.NODELIMIT != -1
              && (int)pPattern->pGraph->v_Nodes.size() > PARAM.NODELIMIT )
          {
            // If extension is not canonical, we should not go deeper
            pPattern->pop_back( false );
            // Exit search. We will find this pattern later
            continue;
          }
          if ( PARAM.DIRECT_HOLEMINING
                    && lastExt.nodeLabelFrom == PARAM.MAXLBL && lastExt.nodeLabelDest == PARAM.MAXLBL
                    && lastExt.direction == gForward )
          {
            // If extension is not canonical, we should not go deeper
            pPattern->pop_back( false );
            // Exit search. We will find this pattern later
            continue;
          }

          //----------------------------------------------------------------------------
          // TIMEOUT STOP
          if ( PARAM.TIMEOUT != -1
              && ( clock() - firstTick ) / CLOCKS_PER_SEC >= round(PARAM.TIMEOUT*3600) )
          {
            cout << " Mining timeout reached during class : " << v_Graphs->at(0)->className << endl;
            cout << " Time of mining : " << ( clock() - firstTick ) / CLOCKS_PER_SEC << endl;
            pPattern->pop_back( false );
            return -1;
          }
          //----------------------------------------------------------------------------
          // NB PATTERN STOP
          if ( PARAM.NBPATLIMIT != -1
              && (int)v_PatternCurrClass.size() == PARAM.NBPATLIMIT )
          {
            cout << " Limit of number of pattern reached during class : " << v_Graphs->at(0)->className << endl;
            pPattern->pop_back( false );
            return -2;
          }
        }


        cur->valid_extenstions->insert(make_pair(lastExt,it->second));

        // Restore the value of pPattern
        pPattern->pop_back( false );
      }
      
    }
  }
  
  // chosse one child and the clear children list ?

  map<GToken,GExtensionData,GTokenGt> ::iterator it1;// = get_random_element();// yet to be writen
    //go for a random child
  pPattern->push_back( it1->first, false );
  MCTS_node* new_child = new MCTS_node();
  cur->children_nodes->operator[](it1->first) = new_child;
  
  double ret;
  if(true) // a stopping condition needed
    ret = rool_out(new_child,
                  false, 
                  v_Graphs, minFreq, pPattern,
                  tokenData, it1->second,suppData);


  // because we only save the first node occList other wise no memory will be enough
  // which means that here is where the actual expansion is happening
  if(!first_level){
    //clear everything that you allocate here
    cur->valid_extenstions->clear();
    cur->children_nodes->clear();
    delete new_child;
  }

  return ret;
}

void MCTSGrima::update_ancestors(MCTS_node* cur, double delta){ //wot??
  while(cur!=NULL){
    cur->Q = (cur->N_node*cur->Q+delta)/(cur->N_node+1);
    cur->N_node = cur->N_node+1;
    cur = cur->parent;
  }
}

