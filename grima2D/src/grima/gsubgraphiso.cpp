/***************************************************************************
 *   Copyright (C) 2005 by Siegfried Nijssen                               *
 *   snijssen[at]informatik.uni-freiburg.de                                *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2010,2011 by Adriana Prado and Baptiste Jeudy           *
 *   baptiste.jeudy[at]univ-st-etienne fr                                  *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2014 by Romain Deville                                  *
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
#include <iostream>
// Include project class
#include "gsubgraphiso.hpp"
//#include "../holeG/hgrimaangle.hpp"

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//
//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
//---- PUBLIC --------------------------------------------------------------//
// Public CONSTANTS ________________________________________________________//
// Public Constructor/Desctructor __________________________________________//
GSubgraphIso::GSubgraphIso( GPattern *pPat, vector<GNodeID> *v_OccList ):
  pPattern( pPat ),
  largeGraph( NULL ),
  offset(1)
{
  /*
   * @brief GSubgraphIso
   * Default constructor
   * Define the pattern that must be found and the vector that store occurence
   * @param pPat, the pattern that must be found
   * @param v_OccList, A vector of GNodeId that store the occurence
   */
  v_LargeGraphIdx = v_OccList;
  v_LastEdgeUsed.resize( pPattern->pGraph->v_Nodes.size() );
  v_PatternIdx.resize( 0, GNONODEID );
}
// End of GSubgraphIso::GSubgraphIso( GPattern *pPat,
//                                            vector<GNodeID> *v_OccList ):
GSubgraphIso::~GSubgraphIso()
{
  // Default Destructor
  // Nothing to do here
}
// End of GSubgraphIso::~GSubgraphIso()

// Accessor ________________________________________________________________//
// Mutator _________________________________________________________________//
// Public Methods __________________________________________________________//
void GSubgraphIso::clearOccList()
{
  /*
   * @brief clearOccList
   * Method that fill v_LargeGraphIdx with GNONODEID
   */
  v_LargeGraphIdx->clear();
//  for ( uint i = 0; i < v_LargeGraphIdx->size(); i++ )
//    v_LargeGraphIdx->at(i) = GNONODEID;
}
// End of GSubgraphIso::clearOccList()

bool GSubgraphIso::run( GGraph* pNewLargeGraph,
                        GNodeID nodeLargeIdFrom,
                        GEdgeID edgeLargeId )
{
  /*
   * @brief run
   * Function that returns true if an occurence of the pattern is found in
   * pNewLargeGraph. largeId and edgeId indicate where the occurence must be
   * look for in largeGraph. The first edge of the pattern must already match
   * @param pNewLargeGraph :pointer to the large graph
   * @param nodeLargeIdFrom: node ID from of the first edge in the large graph
   * @param edgeLargeId: edge position in the list of edge of the node from
   * @return TRUE if edge wear occurence, FALSE else.
   */
  // Get the first token of the pattern to start isomorphism test
  GToken &token = pPattern->v_Tokens.at(0);

  // If this is the first call of the function, ie there is no largeGraph.
  if ( largeGraph != NULL )
  {
    // Fill v_PatternIndex (and v_LargeGraphIdx if needed) with GNONODEID
    for ( uint i = offset; i < offset +  pPattern->pGraph->v_Nodes.size(); i++)
    {
      // For each nodes of the largeGraph
      if ( v_LargeGraphIdx->at(i) != GNONODEID )
      {
        if ( v_PatternIdx.at( v_LargeGraphIdx->at(i) ) != ( i - offset) )
          cerr << "ERROR - INTERNAL ERROR While managing isomorphism node index." << endl;

        v_PatternIdx.at( v_LargeGraphIdx->at(i) ) = GNONODEID;

        // Delete last occ (wich was incomplete)
        if ( !lastOccValid)
          v_LargeGraphIdx->at(i) = GNONODEID;
      }
    }
    // If last occurence "wear" the pattern, then update offset
    if ( lastOccValid )
      offset += patternNbNode + 1;
  }
  patternNbNode    = pPattern->pGraph->v_Nodes.size();

  // Resize v_LargeGraphIdx if needed and fill new space with GNONODEID
  if ( v_LargeGraphIdx->size() < offset + patternNbNode )
    v_LargeGraphIdx->resize( offset + patternNbNode, GNONODEID );

  // If working on new largeGraph
  if ( largeGraph != pNewLargeGraph )
  {
    // Assign current the newLargeGraph as the current largeGraph and update
    // Number of node
    largeGraph = pNewLargeGraph;

    // Resize v_PatternIdx and full with GNONODEID if needed
    if ( v_PatternIdx.size() != largeGraph->v_Nodes.size())
      v_PatternIdx.resize( largeGraph->v_Nodes.size(), GNONODEID );
  }
  
  // Set the id of the largeGraph in the v_LargeGraphIdx
  v_LargeGraphIdx->at(offset-1) = largeGraph->graphID;
  
  // Get the corresponding edge in largeGraph
  GGraphEdge &edge =
      largeGraph->v_Nodes.at(nodeLargeIdFrom).v_Edges.at(edgeLargeId);
  
  // Map node of P into largeGraph
  // Token.nodeFrom in P maps nodeLargeId in the graph largeGraph
  v_PatternIdx[nodeLargeIdFrom] = token.nodeFrom;
  v_LargeGraphIdx->at(offset + token.nodeFrom) = nodeLargeIdFrom;
  // And token.nodeDest in P mapts to edge.nodeDest in largeGraph
  v_PatternIdx[edge.destNodeId] = token.nodeDest;
  v_LargeGraphIdx->at(offset + token.nodeDest) = edge.destNodeId;

  // Keep track of the last edge used in the match
  v_LastEdgeUsed.at(token.nodeFrom) = edgeLargeId;
  v_LastEdgeUsed.at(token.nodeDest) = edge.eqEdgeId;

  if ( pPattern->v_Tokens.size() == 1 )
    lastOccValid = true;
  else if ( subRun() )
    lastOccValid = true;
  else
  {
    lastOccValid = false;
    // Erase last occ which was incomplete
    for ( uint i = offset; i < offset +  pPattern->pGraph->v_Nodes.size(); i++)
    {
      // For each nodes of the largeGraph
      if ( v_LargeGraphIdx->at(i) != GNONODEID )
      {
        v_PatternIdx.at( v_LargeGraphIdx->at(i) ) = GNONODEID;
        v_LargeGraphIdx->at(i) = GNONODEID;
      }
    }
    v_LargeGraphIdx->at( offset -1 ) = GNONODEID;
  }

  return lastOccValid;
}
// End of GSubgraphIso::run( GGraph* pNewLargeGraph,
//                           GNodeID nodeLargeIdFrom,
//                           GEdgeID edgeLargeId )

GEdgeAngle GSubgraphIso::getTokenAngleValue( uint    prevEdgeFrom,
                                             GNodeID nodeFrom,
                                             GNodeID nodeDest )
{
  /*
   * TODO : RD
   * Write Desc
   */
  if ( prevEdgeFrom < 4 )
  {
    // If previous edge is spatial
    uint currEdgePos = largeGraph->getEdgePositionFrom( nodeFrom, nodeDest );
    if ( currEdgePos >= 4 )
    {
      if ( currEdgePos ==  4 )
        return -1;
      else
        return -2;
      // return currEdgePos;
    }
    else if ( prevEdgeFrom == 0 )
    {
      if ( currEdgePos == 0 )
        return 0;
      else if ( currEdgePos == 1 )
        return 3;
      else if ( currEdgePos == 2 )
        return 2;
      else if ( currEdgePos == 3 )
        return 1;
    }
    else if ( prevEdgeFrom == 1 )
    {
      if ( currEdgePos == 0 )
        return 1;
      else if ( currEdgePos == 1 )
        return 0;
      else if ( currEdgePos == 2 )
        return 3;
      else if ( currEdgePos == 3 )
        return 2;
    }
    else if ( prevEdgeFrom == 2 )
    {
      if ( currEdgePos == 0 )
        return 2;
      else if ( currEdgePos == 1 )
        return 1;
      else if ( currEdgePos == 2 )
        return 0;
      else if ( currEdgePos == 3 )
        return 3;
    }
    else if ( prevEdgeFrom == 3 )
    {
      if ( currEdgePos == 0 )
        return 3;
      else if ( currEdgePos == 1 )
        return 2;
      else if ( currEdgePos == 2 )
        return 1;
      else if ( currEdgePos == 3 )
        return 0;
    }
    //    return ( currEdgePos - prevEdgeFrom ) % 4 ;

  }
  else if ( pPattern->v_Tokens.at(0).angle < 0 )
  {
    uint currEdgePos = largeGraph->getEdgePositionFrom( nodeFrom, nodeDest );
    if ( currEdgePos ==  4 )
      return -1;
    else if ( currEdgePos ==  5 )
      return -2;
    else
      return GNOANGLE;
  }
  else
  {
    // Previous edge is temporal, need to find first spatial edge that arrive on
    // nodeFrom
    uint iToken;
    GNodeID nodeLargeFrom ;
    GNodeID nodeLargeDest ;

    iToken = findNotTemporalToken( v_PatternIdx.at( nodeFrom ) );

    if ( pPattern->v_Tokens.at(iToken).nodeFrom == 0 && iToken != 0 )
    {
      // Previous edge was temporal starting from 0, then take edge (0,1)
      nodeLargeFrom = v_LargeGraphIdx->at( offset + 1 );
      nodeLargeDest = v_LargeGraphIdx->at( offset + 0 );
    }
    else
    {
      // Previous edge was temporal not starting from 0
      nodeLargeFrom = v_LargeGraphIdx->at( offset + pPattern->v_Tokens.at(iToken).nodeFrom );
      nodeLargeDest = v_LargeGraphIdx->at( offset + pPattern->v_Tokens.at(iToken).nodeDest );
    }
    uint currSpatEdgePos = largeGraph->getEdgePositionFrom( nodeFrom, nodeDest );
    uint prevSpatEdgePos = largeGraph->getEdgePositionFrom( nodeLargeFrom, nodeLargeDest );

    if ( currSpatEdgePos >= 4 )
    {
      if ( currSpatEdgePos ==  4 )
        return -1;
      else
        return -2;
    }
    else if ( prevSpatEdgePos == 0 )
    {
      if ( currSpatEdgePos == 0 )
        return 2;
      else if ( currSpatEdgePos == 1 )
        return 1;
      else if ( currSpatEdgePos == 2 )
        return 0;
      else if ( currSpatEdgePos == 3 )
        return 3;
    }
    else if ( prevSpatEdgePos == 1 )
    {
      if ( currSpatEdgePos == 0 )
        return 3;
      else if ( currSpatEdgePos == 1 )
        return 2;
      else if ( currSpatEdgePos == 2 )
        return 1;
      else if ( currSpatEdgePos == 3 )
        return 0;
    }
    else if ( prevSpatEdgePos == 2 )
    {
      if ( currSpatEdgePos == 0 )
        return 0;
      else if ( currSpatEdgePos == 1 )
        return 3;
      else if ( currSpatEdgePos == 2 )
        return 2;
      else if ( currSpatEdgePos == 3 )
        return 1;
    }
    else if ( prevSpatEdgePos == 3 )
    {
      if ( currSpatEdgePos == 0 )
        return 1;
      else if ( currSpatEdgePos == 1 )
        return 0;
      else if ( currSpatEdgePos == 2 )
        return 3;
      else if ( currSpatEdgePos == 3 )
        return 2;
    }
  }
  return GNOANGLE;
}
// End of GSubgraphIso::getAngleValue( uint    prevEdgeFrom,
//                                     GNodeID nodeFrom,
//                                     GNodeID nodeDest )

//---- PROTECTED  ----------------------------------------------------------//
// Protected CONSTANTS _____________________________________________________//
// Protected Methods _______________________________________________________//

//---- PRIVATE -------------------------------------------------------------//
// Private CONSTANTS _______________________________________________________//
// Private Methods _________________________________________________________//
bool GSubgraphIso::subRun()
{
  /*
   * TODO : RD
   * Copy desc
   */

  GEdgeAngle prevEdgeFrom;
  uint iToken;
  uint smallGraphPos = 1;
  GNodeID largeNodeIdFrom;
  GNodeID largeNodeIdDest;
  GGraphEdge edge;

  while ( true )
  {
    if ( smallGraphPos == pPattern->v_Tokens.size() )
    {
      /* If position of the subgraph is the same of the size of the pattern
       * that means matching succeed
       */
      return true;
    }
    // Else try tso match the the next edge of pattern
    GToken &token = pPattern->v_Tokens.at(smallGraphPos);
    if ( token.direction == gForward )
    {
      if ( token.nodeFrom == 0 )
      {
        // If first extension starting from node 0 in pattern
        largeNodeIdFrom = v_LargeGraphIdx->at(offset+1);
        largeNodeIdDest = v_LargeGraphIdx->at(offset+0);
        prevEdgeFrom = largeGraph->getEdgePositionFrom(
              largeNodeIdFrom, largeNodeIdDest );
      }
      else if ( token.direction == gForward )
      {
        // Find previous edge in v_Tokens
        iToken = 0;
        while ( pPattern->v_Tokens.at(iToken).nodeDest != token.nodeFrom )
          iToken++;

        largeNodeIdFrom = v_LargeGraphIdx->at(offset+pPattern->v_Tokens.at(iToken).nodeFrom);
        largeNodeIdDest = v_LargeGraphIdx->at(offset+pPattern->v_Tokens.at(iToken).nodeDest);
        prevEdgeFrom = largeGraph->getEdgePositionFrom(
              largeNodeIdFrom, largeNodeIdDest );
      }

      if ( prevEdgeFrom == GNOANGLE )
      {
        // If there was an error in finding previous edge orientation
        cerr << "ERROR - Error, subgraphiso find an edge with GNOANGLE"  << endl;
        exit ( EXIT_FAILURE );
      }

      edge = getEdge( token, prevEdgeFrom );
      if ( edge.destNodeId == GNONODEID )
      {
        // If edge is null, not an valid extension
        return false;
      }
      if ( edge.label == token.edgeLabel
           && largeGraph->v_Nodes.at(edge.destNodeId).label == token.nodeLabelDest
           && v_PatternIdx.at(edge.destNodeId) == GNONODEID)
      {
        // If edge in large graph is the same as the pattern ie extension is valid
        // Store mapping node
        v_PatternIdx.at(edge.destNodeId) = token.nodeDest;
        v_LargeGraphIdx->at(offset+token.nodeDest) = edge.destNodeId;
        // Increase size of subgraph counter
        smallGraphPos++;
      }
      else
      {
        // Edge in the large graph not the same of the pattern
        // ie extension not valid
        return false;
      }
    }
    else
    {
      // token.direction == grimaBackward
      // Just to be sure, verify that edge exist in largeGraph.
      if ( largeGraph->edgeExist(
             v_LargeGraphIdx->at(offset+token.nodeFrom),
             v_LargeGraphIdx->at(offset+token.nodeDest) ) )
        smallGraphPos++;
      else
        return false;
    }
  }
}
// End of GSubgraphIso::subRun()

GGraphEdge GSubgraphIso::getEdge( GToken& token, int prevPosition )
{
  /*
   * @brief getEdge
   * Function that return edge in function of its angle in pattern and the
   * previous edge angle in the large graph
   * @param token: token with angle, ie angle of this edge from the previous one
   * @param prevPosition : orientation of the previous edge in the large graph
   * @return return the edge that extent pattern in the large graph
   */

  GGraphEdge edge;
  if ( token.angle == -1 )
  {
    // Edge is forward in time
    edge = largeGraph->v_Nodes.at( v_LargeGraphIdx->at( offset + token.nodeFrom ) ).v_Edges.at(4);
  }
  else if ( token.angle == -2 )
  {
    // Edge is backward in time
    edge = largeGraph->v_Nodes.at( v_LargeGraphIdx->at( offset + token.nodeFrom ) ).v_Edges.at(5);
  }
  else
  {
    // Previous edge is temporal
    if ( prevPosition >= 4 )
    {
      uint   iToken = findNotTemporalToken( token.nodeFrom );
      GToken tmpToken = pPattern->v_Tokens.at(iToken);
      if ( tmpToken.nodeFrom == 0 && iToken != 0 )
      {
        prevPosition = largeGraph->getEdgePositionFrom(
              v_LargeGraphIdx->at( offset + 1 ),
              v_LargeGraphIdx->at( offset + 0 ) );
      }
      else
      {
        prevPosition = largeGraph->getEdgePositionFrom(
              v_LargeGraphIdx->at( offset + tmpToken.nodeFrom ),
              v_LargeGraphIdx->at( offset + tmpToken.nodeDest ) );
      }
    }
    switch ( token.angle )
    {
    // CCW = Counter ClockWise | CW = ClockWise
    case 0 :
      // Next edge at 90CCW
      edge = largeGraph->v_Nodes.at(
            v_LargeGraphIdx->at(offset+token.nodeFrom) ).v_Edges.at(
            ( prevPosition + 2 ) % 4 );
      break;
    case 1 :
      // Next edge at 90CCW
      edge = largeGraph->v_Nodes.at(
            v_LargeGraphIdx->at(offset+token.nodeFrom) ).v_Edges.at(
            ( prevPosition + 1 ) % 4 );
      break;
    case 2 :
      // Next edge at 180CCW
      edge = largeGraph->v_Nodes.at(
            v_LargeGraphIdx->at(offset+token.nodeFrom) ).v_Edges.at(
            ( prevPosition ) % 4 );
      break;
    case 3 :
      // Next edge at 90CW
      edge = largeGraph->v_Nodes.at(
            v_LargeGraphIdx->at(offset+token.nodeFrom) ).v_Edges.at(
            ( prevPosition + 3 ) % 4 );
      break;
    default :
      cerr << "ERROR - getEdge in GSubGraphIso unable to find edge. "
           << "Wrong angle value !" << endl;
      exit( EXIT_FAILURE );
      break;
    }
  }
  return edge;
}
// End of GSubgraphIso::getEdge( GToken& token, int prevPosition )

int GSubgraphIso::findNotTemporalToken( GNodeID patNodeFrom )
{
  uint iToken = 0;
  while ( pPattern->v_Tokens.at(iToken).nodeDest != patNodeFrom
         || pPattern->v_Tokens.at(iToken).direction != gForward )
    iToken++;

  if ( pPattern->v_Tokens.at(iToken).angle < 0 )
  {
    if ( pPattern->v_Tokens.at(iToken).nodeFrom == 0 )
      return iToken;
    else
      return findNotTemporalToken( pPattern->v_Tokens.at(iToken).nodeFrom );
  }
  else
    return iToken;
}
