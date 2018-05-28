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
#include <algorithm>
// Include project class
#include "gpattern.hpp"
#include "gcanonicaldfscomputer.hpp"
#include "../global.hpp"

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//
//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
//---- PUBLIC --------------------------------------------------------------//
// Public CONSTANTS ________________________________________________________//
// Public Structure & Typedef ______________________________________________//
// Public Constructor/Desctructor __________________________________________//
GPattern::GPattern():
  minTCoord(0),
  maxTCoord(0),
  nodeMinTCoord(-1),
  nodeMaxTCoord(-1)
{
  // Default constructor
  pGraph    = new GGraph();
  v_OccList = new vector<GNodeID>();
  v_OccList->clear();
  // Set default size of OccList to 100 in order to avoid null pointer
  //v_OccList->resize( 100, GNONODEID );
}

GPattern::GPattern(GPattern* pat):
  minTCoord(0),
  maxTCoord(0),
  nodeMinTCoord(-1),
  nodeMaxTCoord(-1)
{
  v_Tokens = vector<GToken>(pat->v_Tokens);
  pGraph    = new GGraph(pat->pGraph);
  v_OccList = new vector<GNodeID>(*pat->v_OccList);
  minTCoord = pat->minTCoord;
  maxTCoord = pat->maxTCoord;
  nodeMinTCoord = pat->nodeMinTCoord;
  nodeMaxTCoord = pat->nodeMaxTCoord;
  
}
// End of GPattern::GPattern()

GPattern::~GPattern()
{
  // Default destructor
  delete pGraph;
  delete v_OccList;
  v_Tokens.clear();
  v_Tokens.shrink_to_fit();
}
// End of GPattern::~GPattern()

// Accessor ________________________________________________________________//
// Mutator _________________________________________________________________//
// Public Methods __________________________________________________________//
void GPattern::read( GReader *data, vector<string> &v_Token )
{
  // === Initialisation

  // While reading "V <nodeId> <nodeLabel>" line
  while ( v_Token[0] == "c" )
  {
    // Extract only data information parsed in v_Tokens
    char tmp[v_Token[1].size()];
    char* token = NULL;
    string substr =  v_Token[1].substr(1,v_Token[1].size()-2).c_str();
    v_Token.clear();
    strcpy( tmp, substr.c_str());
    token = strtok(tmp, ",");
    while (token != NULL)
    {
      v_Token.push_back(token);
      token = strtok(NULL, ",");
    }

    // Get the edge token
    GToken tokenDFS;
    tokenDFS.nodeFrom = atoi(v_Token[1].c_str());
    tokenDFS.nodeDest = atoi(v_Token[2].c_str());
    tokenDFS.angle    = atoi(v_Token[3].c_str());
    tokenDFS.nodeLabelFrom = atoi(v_Token[4].c_str());
    tokenDFS.edgeLabel     = atoi(v_Token[5].c_str());
    tokenDFS.nodeLabelDest = atoi(v_Token[6].c_str());
    if ( atoi( v_Token[0].c_str() ) == 0 )
      tokenDFS.direction = gForward;
    else
      tokenDFS.direction = gBackward;

    push_back( tokenDFS, false  );
    data->getLine( &v_Token );
  }
  // As we reach the end of the graph, we follow on the next one.
  while ( !data->finished() && v_Token[0] != "p" )
    data->getLine( &v_Token );
}
// End of GPattern::read( GReader *data, vector<string> &v_Token )


void GPattern::push_back( const  GToken &v_Ext, bool canonTest )
{
  /*
   * TODO : RD
   * Copy Desc
   */
  GEdgeAngle prevEdgeFrom;
  uint nbToken;
  uint iToken;
  int nodeDestXCoord;
  int nodeDestYCoord;
  int nodeDestTCoord;

  if ( v_Ext.direction != gBackward )
  {
    if ( v_Tokens.size() == 0 )
      setup( v_Ext );

    // Add extension to pattern code
    v_Tokens.push_back( v_Ext ) ;
    nbToken = v_Tokens.size();

    nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord;
    nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord;
    nodeDestTCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).tCoord;


    if ( nbToken == 1 && v_Ext.angle == 0 )
    {
      // If pattern code has only one edge, this mean for the moment, only
      // one node at (0,0,0) as been added. Thus, if this code is not a temporal
      // one, we need to add second node at (1,0,0).
      nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord + 1;
    }
    else if ( v_Ext.angle == -1 )
    {
      // Extension is forward in time
      nodeDestTCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).tCoord + 1;
    }
    else if ( v_Ext.angle == -2 )
    {
      // Extension is backward in time
      nodeDestTCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).tCoord - 1;
    }
    else if ( nbToken > 1 )
    {
      if ( v_Tokens.at(0).angle >= 0 )
      {
        if ( v_Ext.nodeFrom == 0 )
        {
          // If pattern does not continue from continue from last edge wich is not a backward edge
          prevEdgeFrom = pGraph->getEdgePositionFrom( 1, 0 );
        }
        else
        {
          if ( v_Tokens.at(nbToken-2).nodeDest == v_Ext.nodeFrom
               &&  v_Tokens.at(nbToken-2).direction != gBackward  )
          {
            // If pattern does not continue from continue from last edge wich is not a backward edge
            prevEdgeFrom = pGraph->getEdgePositionFrom(
                  v_Tokens.at(nbToken-2).nodeFrom, v_Ext.nodeFrom );
          }
          else
          {
            // Else, pattern continue from previously seen edge
            uint iToken = 0;
            while ( v_Tokens.at(iToken).nodeDest != v_Ext.nodeFrom )
              iToken++;

            prevEdgeFrom = pGraph->getEdgePositionFrom(
                  v_Tokens.at(iToken).nodeFrom, v_Ext.nodeFrom );
          }
        }
        if ( prevEdgeFrom >= 4 )
        {
          uint iTokenPrev = findNotTemporalToken( v_Ext.nodeFrom ) ;
          if ( v_Tokens.at(iTokenPrev).nodeFrom == 0
               && v_Tokens.at(iTokenPrev).nodeDest != 1 )
            prevEdgeFrom = pGraph->getEdgePositionFrom( 1, 0 );
          else
            prevEdgeFrom = pGraph->getEdgePositionFrom( v_Tokens.at(iTokenPrev).nodeFrom,
                                                        v_Tokens.at(iTokenPrev).nodeDest );
        }
        switch ( prevEdgeFrom )
        {
        case 0:
          /*
           * B--A
           */
          if ( v_Ext.angle == 0 )
            nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord + 1 ;
          else if ( v_Ext.angle == 1 )
            nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord + 1;
          else if ( v_Ext.angle == 2 )
            nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord - 1 ;
          else if ( v_Ext.angle == 3 )
            nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord - 1;
          break;
        case 1:
          /*
           * B
           * |
           * A
           */
          if ( v_Ext.angle == 0 )
            nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord - 1;
          else if ( v_Ext.angle == 1 )
            nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord + 1;
          else if ( v_Ext.angle == 2 )
            nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord + 1 ;
          else if ( v_Ext.angle == 3 )
            nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord - 1;
          break;
        case 2:
          /*
           * A--B
           */
          if ( v_Ext.angle == 0 )
            nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord - 1 ;
          else if ( v_Ext.angle == 1 )
            nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord - 1;
          else if ( v_Ext.angle == 2 )
            nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord + 1 ;
          else if ( v_Ext.angle == 3 )
            nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord + 1;
          break;
        case 3:
          /*
           * A
           * |
           * B
           */
          if ( v_Ext.angle == 0 )
            nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord + 1;
          else if ( v_Ext.angle == 1 )
            nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord - 1;
          else if ( v_Ext.angle == 2 )
            nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord - 1 ;
          else if ( v_Ext.angle == 3 )
            nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord + 1;
          break;
        default:
          cerr << "ERROR - Error in gpattern, trying to push an extension with "
               << "angle of GNOANGLE" << endl;
          exit( EXIT_FAILURE );
          break;
        }
      }
      else
      {
        // First check if there is already a spatial edge
        iToken = 0;
        while ( iToken < nbToken && v_Tokens.at(iToken).angle != 0 )
          iToken++;
        if ( iToken == nbToken - 1 )
        {
          // If there is no already existing spatial edge except the last one added
          nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord + 1;
        }
        else
        {
          // As it already exist spatial edge, this one will be take as reference
          // Find previous spatial edge
          uint iPrevToken  = 1;
          uint prevSpatPos = -1;
          while ( iPrevToken < nbToken
                  && v_Tokens.at(iPrevToken).nodeDest != v_Ext.nodeFrom )
            iPrevToken++;
          if ( iPrevToken == v_Tokens.size() )
            // If edge start from node that have no previous spatial edge then
            // take first spatial edge as reference
            prevSpatPos = pGraph->getEdgePositionFrom( v_Tokens.at(iToken).nodeFrom,
                                                       v_Tokens.at(iToken).nodeDest );
          else
            prevSpatPos = pGraph->getEdgePositionFrom( v_Tokens.at(iPrevToken).nodeFrom,
                                                       v_Tokens.at(iPrevToken).nodeDest );

          switch ( prevSpatPos )
          {
          case 0:
            /*
               * B--A
               */
            if ( v_Ext.angle == 0 )
              nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord + 1 ;
            else if ( v_Ext.angle == 1 )
              nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord + 1;
            else if ( v_Ext.angle == 2 )
              nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord - 1 ;
            else if ( v_Ext.angle == 3 )
              nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord - 1;
            break;
          case 1:
            /*
               * B
               * |
               * A
               */
            if ( v_Ext.angle == 0 )
              nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord - 1;
            else if ( v_Ext.angle == 1 )
              nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord + 1;
            else if ( v_Ext.angle == 2 )
              nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord + 1 ;
            else if ( v_Ext.angle == 3 )
              nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord - 1;
            break;
          case 2:
            /*
               * A--B
               */
            if ( v_Ext.angle == 0 )
              nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord - 1 ;
            else if ( v_Ext.angle == 1 )
              nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord - 1;
            else if ( v_Ext.angle == 2 )
              nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord + 1 ;
            else if ( v_Ext.angle == 3 )
              nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord + 1;
            break;
          case 3:
            /*
               * A
               * |
               * B
               */
            if ( v_Ext.angle == 0 )
              nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord + 1;
            else if ( v_Ext.angle == 1 )
              nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord - 1;
            else if ( v_Ext.angle == 2 )
              nodeDestYCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).yCoord - 1 ;
            else if ( v_Ext.angle == 3 )
              nodeDestXCoord = pGraph->v_Nodes.at( v_Ext.nodeFrom ).xCoord + 1;
            break;
          default:
            cerr << "ERROR - Error in gpattern, trying to push an extension with "
                 << "angle of GNOANGLE" << endl;
            exit( EXIT_FAILURE );
            break;
          }
        }
      }
    }
    pGraph->addNode( v_Ext.nodeLabelDest,
                     nodeDestXCoord, nodeDestYCoord, nodeDestTCoord );
    pGraph->addEdge( v_Ext.nodeFrom, v_Ext.nodeDest, v_Ext.edgeLabel );

    if ( (int)nodeMaxTCoord == -1 && nodeDestTCoord >= maxTCoord )
    {
      maxTCoord     = nodeDestTCoord;
      nodeMaxTCoord = v_Ext.nodeDest;
    }
    else if ( nodeDestTCoord > maxTCoord )
    {
      maxTCoord     = nodeDestTCoord;
      nodeMaxTCoord = v_Ext.nodeDest;
    }
    if ( (int)nodeMinTCoord == -1 && nodeDestTCoord <= minTCoord )
    {
      minTCoord = nodeDestTCoord;
      nodeMinTCoord = v_Ext.nodeDest;
    }
    else if ( nodeDestTCoord < minTCoord )
    {
      minTCoord = nodeDestTCoord;
      nodeMinTCoord = v_Ext.nodeDest;
    }

    if ( PARAM.INDUCED &&  v_Tokens.size() > 2 && !canonTest )
    {
      vector<GToken> v_BackwardToken;
      for ( uint iNodeFrom = 0; iNodeFrom < pGraph->v_Nodes.size(); iNodeFrom++ )
      {
        int prevPosFrom;
        int deltaX = pGraph->v_Nodes.at(iNodeFrom).xCoord - pGraph->v_Nodes.at(v_Ext.nodeDest).xCoord;
        int deltaY = pGraph->v_Nodes.at(iNodeFrom).yCoord - pGraph->v_Nodes.at(v_Ext.nodeDest).yCoord;
        int deltaT = pGraph->v_Nodes.at(iNodeFrom).tCoord - pGraph->v_Nodes.at(v_Ext.nodeDest).tCoord;
        if ( abs( deltaX ) + abs ( deltaY ) + abs( deltaT ) == 1 )
        {
          GToken tmpExt;
          tmpExt.direction = gBackward;
          tmpExt.nodeFrom  = v_Ext.nodeDest;
          tmpExt.nodeDest  = iNodeFrom;
          tmpExt.nodeLabelFrom = pGraph->v_Nodes.at(v_Ext.nodeDest).label;
          tmpExt.nodeLabelDest = pGraph->v_Nodes.at(iNodeFrom).label;
          tmpExt.edgeLabel = 0;

          if ( deltaT == 1 )
            tmpExt.angle     = -1;
          else if ( deltaT == -1 )
            tmpExt.angle     = -2;
          else
          {
            if ( v_Ext.nodeFrom == 0 )
              prevPosFrom = pGraph->getEdgePositionFrom( 1, 0 );
            else
              prevPosFrom = pGraph->getEdgePositionFrom( v_Ext.nodeFrom,
                                                         v_Ext.nodeDest );
            if ( prevPosFrom >= 4  )
            {
              prevPosFrom = findNotTemporalToken( v_Ext.nodeFrom );
              prevPosFrom = pGraph->getEdgePositionFrom( v_Tokens.at(prevPosFrom).nodeFrom,
                                                         v_Tokens.at(prevPosFrom).nodeDest );
            }
            if ( deltaX == -1 )
            {
              if ( prevPosFrom == 0 )
                tmpExt.angle = 2;
              else if ( prevPosFrom == 1 )
                tmpExt.angle = 3;
              else if ( prevPosFrom == 2 )
                tmpExt.angle = 0;
              else if ( prevPosFrom == 3 )
                tmpExt.angle = 1;
            }
            else if ( deltaX == 1 )
            {
              if ( prevPosFrom == 0 )
                tmpExt.angle = 0;
              else if ( prevPosFrom == 1 )
                tmpExt.angle = 1;
              else if ( prevPosFrom == 2 )
                tmpExt.angle = 2;
              else if ( prevPosFrom == 3 )
                tmpExt.angle = 3;
            }
            else if ( deltaY == 1 )
            {
              if ( prevPosFrom == 0 )
                tmpExt.angle = 1;
              else if ( prevPosFrom == 1 )
                tmpExt.angle = 2;
              else if ( prevPosFrom == 2 )
                tmpExt.angle = 3;
              else if ( prevPosFrom == 3 )
                tmpExt.angle = 0;
            }
            else if ( deltaY == -1 )
            {
              if ( prevPosFrom == 0 )
                tmpExt.angle = 3;
              else if ( prevPosFrom == 1 )
                tmpExt.angle = 0;
              else if ( prevPosFrom == 2 )
                tmpExt.angle = 1;
              else if ( prevPosFrom == 3 )
                tmpExt.angle = 2;
            }
          }
          v_BackwardToken.push_back( tmpExt );
        }
      }

      sort(v_BackwardToken.begin(), v_BackwardToken.end(), cmpGToken );

      for ( uint iExt = 0; iExt < v_BackwardToken.size(); iExt ++ )
      {
        GToken tmpExt = v_BackwardToken.at(iExt);
        if ( !pGraph->edgeExist( tmpExt.nodeFrom, tmpExt.nodeDest) )
        {
            push_back( tmpExt, canonTest );
        }
      }
    }
  }
  else
  {
    if ( !pGraph->edgeExist( v_Ext.nodeFrom, v_Ext.nodeDest) )
    {
      v_Tokens.push_back( v_Ext );
      if ( PARAM.TORGRAPH )
      {
        int deltaX = pGraph->v_Nodes.at(v_Ext.nodeDest).xCoord - pGraph->v_Nodes.at(v_Ext.nodeFrom).xCoord;
        int deltaY = pGraph->v_Nodes.at(v_Ext.nodeDest).yCoord - pGraph->v_Nodes.at(v_Ext.nodeFrom).yCoord;
        int deltaT = pGraph->v_Nodes.at(v_Ext.nodeDest).tCoord - pGraph->v_Nodes.at(v_Ext.nodeFrom).tCoord;
        // Check neighborhood condition
        if ( abs( deltaX ) + abs ( deltaY ) + abs( deltaT ) == 1 )
          pGraph->addEdge( v_Ext.nodeFrom, v_Ext.nodeDest, v_Ext.edgeLabel );
        else
          pGraph->addEdgeTor( v_Ext.nodeFrom, v_Ext.nodeDest, v_Ext.edgeLabel );
      }
      else
        pGraph->addEdge( v_Ext.nodeFrom, v_Ext.nodeDest, v_Ext.edgeLabel );
    }
    else
    {
      cerr << "ERROR - In gPattern, trying to add already existing backward edge "
           << "between " << v_Ext.nodeFrom << " and " << v_Ext.nodeDest << endl;
    }
  }
}
// End of void GPattern::push_back( const  GToken &v_Ext )

void GPattern::pop_back( bool canonTest )
{
  /*
   * @brief pop_back
   * Remove last extension(s) from the pattern. Undo what was done in last
   * push_back. If backward edge was added, then it's removed.
   */
  GToken &token = v_Tokens.back();
  if ( PARAM.INDUCED && !canonTest )
  {
    while ( token.direction == gBackward )
    {
      pGraph->deleteEdge( token.nodeFrom, token.nodeDest );
      v_Tokens.pop_back();
      token = v_Tokens.back();
    }
  }
  pGraph->deleteEdge( token.nodeFrom, token.nodeDest );
  v_Tokens.pop_back();

  if ( token.direction == gForward )
    pGraph->deleteNode();

  if ( token.nodeDest == nodeMinTCoord || token.nodeDest == nodeMaxTCoord )
  {
    maxTCoord = 0;
    minTCoord = 0;
    nodeMaxTCoord = -1;
    nodeMinTCoord = -1;
    for ( uint iNode = 0; iNode < pGraph->v_Nodes.size(); iNode++ )
    {
      if ( (int)nodeMaxTCoord == -1 && pGraph->v_Nodes.at(iNode).tCoord >= maxTCoord )
      {
        maxTCoord     = pGraph->v_Nodes.at(iNode).tCoord;
        nodeMaxTCoord = iNode;
      }
      else if ( pGraph->v_Nodes.at(iNode).tCoord > maxTCoord )
      {
        maxTCoord     = pGraph->v_Nodes.at(iNode).tCoord;
        nodeMaxTCoord = iNode;
      }
      if ( (int)nodeMinTCoord == -1 && pGraph->v_Nodes.at(iNode).tCoord <= minTCoord )
      {
        minTCoord = pGraph->v_Nodes.at(iNode).tCoord;
        nodeMinTCoord = iNode;
      }
      else if ( pGraph->v_Nodes.at(iNode).tCoord < minTCoord )
      {
        minTCoord = pGraph->v_Nodes.at(iNode).tCoord;
        nodeMinTCoord = iNode;
      }
    }
  }
}
// End of GPattern::pop_back()

bool GPattern::isCanonincal()
{
  /*
   * @brief isCanonincal
   * Test if pattern code is canonical
   * @return TRUE if Pattern is canonical, FALSE if it's not.
   */
  GCanonicalDFSComputer computer( this );
  bool res = computer.isCanonincal();
  return res;
}
// End of GPattern::isCanonincal()

GPattern* GPattern::getCanonincal()
{
  /*
   * @brief isCanonincal
   * Test if pattern code is canonical
   * @return TRUE if Pattern is canonical, FALSE if it's not.
   */
  GCanonicalDFSComputer computer( this );
  return computer.getCanonincal();
}
// End of GPattern::isCanonincal()

void GPattern::printOcc()
{
  uint i     = 0;
  uint occId = 0;
  while ( i < v_OccList->size()
          && v_OccList->at(i) != GNONODEID )
  {

    cout << "o " << pGraph->graphID
         << " "  << v_OccList->at(i++)
         << " " << occId;
    occId++;
    for ( uint k = 0 ; k < pGraph->v_Nodes.size(); k++ )
    {
      cout << " " << v_OccList->at(i);
      cout.flush();
      i++;
    }
    cout << endl;
  }
}

//---- PROTECTED  ------------------------------------------------------------//
// Protected CONSTANTS _______________________________________________________//
// Protected Methods _________________________________________________________//

//---- PRIVATE ---------------------------------------------------------------//
// Private CONSTANTS _________________________________________________________//
// Private Methods ___________________________________________________________//
void GPattern::setup( const GToken& v_Ext )
{
  /*
   * @brief setup
   * Method call by push_back for the first token. Clean node in pGraph, add the
   * two first nodes and edge between this of the pattern based on token.
   * @param token : DFS Code of the first edge.
   */

  if ( !pGraph->v_Nodes.empty() )
    pGraph->deleteNode();
  pGraph->addNode( v_Ext.nodeLabelFrom, 0, 0, 0 );
}
// End of GPattern::setup( const GToken& v_Ext )

int GPattern::findNotTemporalToken( GNodeID patNodeFrom )
{
  uint iToken = 0;
  while ( v_Tokens.at(iToken).nodeDest != patNodeFrom
          || v_Tokens.at(iToken).direction != gForward )
    iToken++;

  if ( v_Tokens.at(iToken).angle < 0 )
  {
    if ( v_Tokens.at(iToken).nodeFrom == 0 )
      return iToken;
    else
      return findNotTemporalToken( v_Tokens.at(iToken).nodeFrom );
  }
  else
    return iToken;
}

//============================== OPERATOR OVERLOAD  ==========================//
ostream& operator<<(ostream& stream, GPattern *pPattern)
{
  GGraph graph            = pPattern->pGraph;
  vector<GToken> v_Tokens = pPattern->v_Tokens;

  int nbNodes = graph.v_Nodes.size();

  stream << "p " << graph.graphID << " " << graph.className << endl;
  stream << "n " << v_Tokens.size() << endl
         << v_Tokens;

  // === Output nodes
  for ( int i = 0; i < nbNodes; ++i )
    stream << "v " << i << " "
           << graph.v_Nodes.at(i).label << " "
           << graph.v_Nodes.at(i).xCoord << " "
           << graph.v_Nodes.at(i).yCoord << " "
           << graph.v_Nodes.at(i).tCoord << " "
           << endl;

  // === Output edges
  for ( int i = 0; i < nbNodes ; ++i )
  {
    vector<GGraphEdge> v_Edges = graph.v_Nodes[i].v_Edges;
    for ( int j = 0, m = v_Edges.size(); j < m; ++j )
    {
      if ( v_Edges[j].destNodeId != GNONODEID )
        stream << "e " << i
               << " " << v_Edges.at(j).destNodeId
               << " " << v_Edges.at(j).label <<  endl;
    }
  }

  if ( PARAM.DEBUG_MODE )
    pPattern->printOcc();
  return stream;
}
// End of operator<<(ostream& stream, GPattern *pPattern)

ostream& operator<<(ostream& stream, vector<GToken> v_Tokens)
{
  // === Output Tokens
  for ( uint i = 0 ; i < v_Tokens.size() ; i++ )
    stream << v_Tokens.at(i);
  return stream;
}
// End of ostream& operator<<(ostream& stream, vector<GToken> v_Tokens )

ostream& operator<<(ostream& stream, GToken token )
{
  stream << "c ("
         << token.direction     << ","
         << token.nodeFrom      << ","
         << token.nodeDest      << ","
         << token.angle         << ","
         << token.nodeLabelFrom << ","
         << token.edgeLabel     << ","
         << token.nodeLabelDest << ")" << endl;
  return stream;
}



// End of ostream& operator<<(ostream& stream, GToken tokens )
