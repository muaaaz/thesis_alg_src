/***************************************************************************
 *   Copyright (C) 2005 by Siegfried Nijssen                               *
 *   snijssen[at]informatik.uni-freiburg.de                                   *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2010,2011 by Adriana Prado and Baptiste Jeudy           *
 *   baptiste.jeudy[at]univ-st-etienne fr                                  *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2016 by Romain Deville                                  *
 *   romain.deville[at]insa-lyon.fr                                        *
 * ----------------------------------------------------------------------- *
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
// -- C++ STL Lib
#include <fstream>
#include <algorithm>
// Include project class
#include "../global.hpp"
#include "gdatabase.hpp"

//=============================== NAMESPACE ==================================//
//================================ METHODS ===================================//
//================================= CLASS ====================================//
// ----- PUBLIC --------------------------------------------------------------//
// Public CONSTANTS __________________________________________________________//
// Public Constructor/Desctructor ____________________________________________//
GDatabase::GDatabase()
{
  // Default constructor
  v_GClassDB.resize(0);
}
// End of GDatabase::GDatabase():

GDatabase::~GDatabase()
{
  // Default desctructor
  for ( uint i = 0; i < v_GClassDB.size(); i++ )
    delete v_GClassDB.at(i);
  v_GClassDB.clear();
}
// End of GDatabase::~GDatabase()

// Accessor __________________________________________________________________//
// Mutator ___________________________________________________________________//
// Public Methods ____________________________________________________________//
void GDatabase::createGrapheDB( string &graphFile )
{
  /*
   * @brief createGrapheDB
   * Load graphFile into v_GClassDB.
   * If hole = true, make hole into graph, i.e. if all neighbors of a node have
   * the same label, then this node is deleted (set to GNONODEID).
   * @param graphFile : path to file that store graph
   * @param hole      : boolean to specify if it's needed to make hole in graph
   */
  if ( PARAM.DEBUG_MODE )
    cout << "#- Create Graph Database." << endl;

  if ( ! ifstream( graphFile ) )
    error::fileNotExistError( graphFile );

  read( graphFile );

  PARAM.MAXLBL += 1;
  for ( uint i = 0; i < v_GClassDB.size(); i++ )
  {
    if ( PARAM.INDIRECT_HOLEMINING || PARAM.DIRECT_HOLEMINING || PARAM.FRAMEHOLE )
      processHoleInGrid( v_GClassDB.at(i) );
    processSparseSet( v_GClassDB.at(i) );
  }
}
// End of GDatabase::createGrapheDB( string &graphFile, bool hole )

// ----- PROTECTED  ----------------------------------------------------------//
// Protected CONSTANTS _______________________________________________________//
// Protected Methods _________________________________________________________//
void GDatabase::read( string filename )
{
  /*
   * @brief read
   * Read the file that store graph described in the filename
   * @param filename : File path where graph are stored
   * @return A vector of structure that store graph by class with class infos.
   */
  if ( PARAM.DEBUG_MODE )
    cout << "#- Parse graphs files." << endl;

  // Variable initialization
  uint              suppGraphID; // Supposed graph ID
  GGraph            *graph;
  vector<string>    v_Tokens;
  GClassDB          *tmpGClassDB = new GClassDB();
  GReader           *data = new GReader( filename );
  int tmpClassID = 0;

  // Read first line of type "t <tId> <tClassName>
  data->getLine( &v_Tokens );

  suppGraphID = 0;

  while ( !data->finished() )
  {
    // Instanciate new object graph
    graph = new GGraph();

    // Read one graph in file to store it
    graph->read( data, v_Tokens, suppGraphID, tmpGClassDB->className  );
   

    if (graph->className != tmpGClassDB->className )
    {
      ++tmpClassID;
      // If graph belong to a new classname
      //if ( tmpGClassDB->v_ClassGraphs.size() != 0 )
      //  v_GClassDB.push_back( tmpGClassDB );
      //else
      //  delete tmpGClassDB;

      //tmpGClassDB = new GClassDB();
      tmpGClassDB->className = graph->className;
      suppGraphID = 1;
    }
    graph->classID = tmpClassID;
    // Add graph to class database
    tmpGClassDB->v_ClassGraphs.push_back( graph );
  }
  if ( tmpGClassDB->v_ClassGraphs.size() != 0 )
    v_GClassDB.push_back( tmpGClassDB );
  // Free memory by closing the file.
  delete data;
}
// End of GDatabase::read( string filename )

void GDatabase::processSparseSet( GClassDB *pClassDB )
{
  /*
   * TODO : RD
   * Copy desc from header
   */
  for ( uint iGraph = 0; iGraph < pClassDB->v_ClassGraphs.size(); iGraph++ )
  {
    vector<GGraphNode> v_Nodes = pClassDB->v_ClassGraphs.at(iGraph)->v_Nodes;
    for ( uint iNode = 0; iNode < v_Nodes.size(); iNode++)
      if ( v_Nodes.at(iNode).label != GNOLABEL )
        processEdgeSparseSet( pClassDB, iGraph, iNode, v_Nodes );
  }
}
// End of GDatabase::processSparseSet( GClassDB *pClassDB )

void GDatabase::processEdgeSparseSet( GClassDB           *pClassDB,
                                      uint               iGraph,
                                      uint               iNode,
                                      vector<GGraphNode> &v_Nodes )
{
  /*
   * TODO : RD
   * Copy Desc from header
   */
  vector<GGraphEdge> v_Edges = v_Nodes.at(iNode).v_Edges;
  for ( uint iEdge = 0; iEdge < v_Edges.size(); iEdge++ )
  {
    if ( v_Edges.at(iEdge).destNodeId != GNONODEID )
    {
      GNodeLabel nodeLabelFrom = v_Nodes.at(iNode).label;
      GNodeLabel nodeLabelDest = v_Nodes.at( v_Edges.at(iEdge).destNodeId ).label ;
      if ( nodeLabelFrom >= nodeLabelDest )
      {
        GEdgeLabel edgeLabel = v_Edges.at(iEdge).label;
        uint       iNodeDest = v_Edges.at(iEdge).destNodeId;
        if ( iEdge < 4 )
          addEdgeSpatialSparseSet( nodeLabelFrom, nodeLabelDest, edgeLabel,
                                   pClassDB, iGraph, iNode, iNodeDest, iEdge );
        else
          addEdgeTemporalSparseSet( nodeLabelFrom, nodeLabelDest, edgeLabel,
                                    pClassDB, iGraph, iNode, iNodeDest, iEdge );
      }
    }
  }
}
// End of GDatabase::processEdgeSparseSet( GClassDB           *pClassDB,
//                                         uint               iGraph,
//                                         uint               iNode,
//                                         vector<GGraphNode> &v_Nodes )

void GDatabase::addEdgeSpatialSparseSet( GNodeLabel nodeLabelFrom,
                                         GNodeLabel nodeLabelDest,
                                         GEdgeLabel edgeLabel,
                                         GClassDB   *pClassDB,
                                         uint iGraph,
                                         uint iNodeFrom,
                                         uint iNodeDest,
                                         uint iEdge )
{
  // If first edge is canonical (i.e. label from > label dest)
  pair<GToken,GTokenData> p_Data;
  p_Data.second.freq    = 0;
  p_Data.first.nodeFrom = 0;
  p_Data.first.nodeDest = 1;
  p_Data.first.angle    = 0;
  p_Data.first.nodeLabelFrom = nodeLabelFrom;
  p_Data.first.nodeLabelDest = nodeLabelDest;
  p_Data.first.direction     = gForward;
  p_Data.first.edgeLabel     = edgeLabel;

  // Try to insert edge code in m_TokenData
  pair< map< GToken, GTokenData, GTokenGt>::iterator, bool>
      p = pClassDB->m_TokenData.insert( p_Data );
  // Store this edge into sparse set
  GSparseSet::mapEdge sparseEdge;
  sparseEdge.nodeFrom = iNodeFrom;
  sparseEdge.nodeDest = iNodeDest;
  sparseEdge.edgeId   = iEdge;
  if ( p.first->second.v_SparseOcc.size() == 0
       || p.first->second.v_SparseOcc.back().graphID
       != pClassDB->v_ClassGraphs.at(iGraph)->graphID)
  {
    GGraph *pGraph = pClassDB->v_ClassGraphs.at(iGraph);
    GSparseSet sparseSet( pGraph->graphID, pGraph );
    p.first->second.v_SparseOcc.push_back( sparseSet );
    p.first->second.freq++;
  }
  p.first->second.v_SparseOcc.back().add( sparseEdge );
}
// End of GDatabase::addEdgeSpatialSparseSet( GNodeLabel nodeLabelFrom,
//                                            GNodeLabel nodeLabelDest,
//                                            GEdgeLabel edgeLabel,
//                                            GClassDB   *pClassDB,
//                                            uint iGraph,
//                                            uint iNodeFrom,
//                                            uint iNodeDest,
//                                            uint iEdge )

void GDatabase::addEdgeTemporalSparseSet( GNodeLabel nodeLabelFrom,
                                          GNodeLabel nodeLabelDest,
                                          GEdgeLabel edgeLabel,
                                          GClassDB   *pClassDB,
                                          uint iGraph,
                                          uint iNodeFrom,
                                          uint iNodeDest,
                                          uint iEdge )
{
  // If first edge is canonical (i.e. label from > label dest)
  pair<GToken,GTokenData> p_Data;
  p_Data.second.freq    = 0;
  p_Data.first.nodeFrom = 0;
  p_Data.first.nodeDest = 1;

  if ( iEdge == 4 )
    p_Data.first.angle    = -1;
  else
    p_Data.first.angle    = -2;

  p_Data.first.nodeLabelFrom = nodeLabelFrom;
  p_Data.first.nodeLabelDest = nodeLabelDest;
  p_Data.first.direction     = gForward;
  p_Data.first.edgeLabel     = edgeLabel;

  // Try to insert edge code in m_TokenData
  pair< map< GToken, GTokenData, GTokenGt>::iterator, bool>
      p = pClassDB->m_TokenData.insert( p_Data );
  // Store this edge into sparse set
  GSparseSet::mapEdge sparseEdge;
  sparseEdge.nodeFrom = iNodeFrom;
  sparseEdge.nodeDest = iNodeDest;
  sparseEdge.edgeId   = iEdge;
  if ( p.first->second.v_SparseOcc.size() == 0
       || p.first->second.v_SparseOcc.back().graphID
       != pClassDB->v_ClassGraphs.at(iGraph)->graphID)
  {
    GGraph *pGraph = pClassDB->v_ClassGraphs.at(iGraph);
    GSparseSet sparseSet( pGraph->graphID, pGraph );
    p.first->second.v_SparseOcc.push_back( sparseSet );
    p.first->second.freq++;
  }
  p.first->second.v_SparseOcc.back().add( sparseEdge );
}
// End of GDatabase::addEdgeTemporalSparseSet( GNodeLabel nodeLabelFrom,
//                                             GNodeLabel nodeLabelDest,
//                                             GEdgeLabel edgeLabel,
//                                             GClassDB   *pClassDB,
//                                             uint iGraph,
//                                             uint iNodeFrom,
//                                             uint iNodeDest,
//                                             uint iEdge )

void GDatabase::processHoleInGrid( GClassDB* pClassDB )
{
  /*
   * @brief processHoleInGrid
   * Process hole in grid by calling isNodeToDeleted2D (if grid as tdim=1) or
   * isNodeToDelete2D_t (if grid as tdim > 1) for each node for each graph.
   * This will construct a vector of node "to delete" (i.e. set their label
   * to GNOLABEL and remove their edges).
   * Then, for each node "to delete", call setNodeToNoLabel.
   * @param pClassDB : Class database of grids
   */
  vector<GNodeID> v_idNodeToDelete;
  vector<GNodeID> v_TmpIdNodeToDelete;
  if ( PARAM.DEBUG_MODE )
    cout << "#- Process hole in graphs of class " << pClassDB->className << endl;

  // For all graph of the database
  for ( uint idxGraph = 0;
        idxGraph < pClassDB->v_ClassGraphs.size();
        idxGraph++ )
  {
    if ( PARAM.INDIRECT_HOLEMINING || PARAM.DIRECT_HOLEMINING )
    {
      v_TmpIdNodeToDelete = pClassDB->v_ClassGraphs.at(idxGraph)->processHole();
      v_idNodeToDelete.insert( v_idNodeToDelete.end(),
                               v_TmpIdNodeToDelete.begin(),
                               v_TmpIdNodeToDelete.end() );
    }
    if ( PARAM.FRAMEHOLE )
    {
      v_TmpIdNodeToDelete = pClassDB->v_ClassGraphs.at(idxGraph)->processFrameHole();
      v_idNodeToDelete.insert( v_idNodeToDelete.end(),
                               v_TmpIdNodeToDelete.begin(),
                               v_TmpIdNodeToDelete.end() );
    }
    sort(v_idNodeToDelete.begin(), v_idNodeToDelete.end() );
    pClassDB->v_ClassGraphs.at(idxGraph)->setNodeToNoLabel( v_idNodeToDelete );
    v_idNodeToDelete.clear();
    v_TmpIdNodeToDelete.clear();
  }
}
// End of GDatabase::processHoleInGrid( GClassDB* pClassDB )

// ----- PRIVATE -------------------------------------------------------------//
// Private CONSTANTS _________________________________________________________//
// Private Methods ___________________________________________________________//

//============================== OPERATOR OVERLOAD  ==========================//
