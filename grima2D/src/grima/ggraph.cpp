/***************************************************************************
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
// Include project class#include <vector>
#include <algorithm>
// Include project class
#include "../global.hpp"
#include "ggraph.hpp"

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//
//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//
//================================= CLASS ====================================//
//---- PUBLIC --------------------------------------------------------------//
// Public CONSTANTS ________________________________________________________//
// Public Constructor/Desctructor __________________________________________//
GGraph::GGraph()
{
  // Default constructor
  // Nothing to do here.
}
// End of GGraph::GGraph()

GGraph::GGraph(GGraph *src )
{
  // Copy constructor
  this->className  = src->className;
  this->graphID    = src->graphID;
  this->dim.xdim   = src->dim.xdim;
  this->dim.ydim   = src->dim.ydim;
  this->dim.tdim   = src->dim.tdim;
  this->v_Nodes.resize( src->v_Nodes.size() );
  copy( src->v_Nodes.begin(), src->v_Nodes.end(), this->v_Nodes.begin() );
}
// End of GGraph::GGraph()

GGraph::~GGraph()
{
  // Default destructor
  for ( uint i = 0; i < v_Nodes.size(); i++ )
  {
    v_Nodes[i].v_Edges.clear();
    v_Nodes[i].v_Edges.shrink_to_fit();
  }

  v_Nodes.clear();
  v_Nodes.shrink_to_fit();

  v_NbOccPattern.clear();
  v_NbOccPattern.shrink_to_fit();
}
// End of GGraph::~GGraph()

// Accessor ________________________________________________________________//
// Mutator _________________________________________________________________//
// Public Methods __________________________________________________________//
void GGraph::read( GReader        *data,
                   vector<string> &v_GraphIdToken,
                   uint           &suppGraphID,
                   string         suppClassName )
{
  /*
   * TODO : RD
   * Copy desc
   */
  // === Initialisation
  // List of code
  vector<string> v_Tokens;
  // Line where number of edge is specified
  uint lineNbEdge;
  // Counter to check dimension are respected
  int suppXCoord = 0;
  int suppYCoord = 0;
  int suppTCoord = 0;
  // Current edge id if not set in the file. Also use to check if nbEdge given
  // at line "n" is correct
  GEdgeID nbEdges       = 0;
  GEdgeID currentEdgeId = 0;
  // Vector used when managing torgraph
  vector<GNodeID> v_NodeTop;
  vector<GNodeID> v_NodeBottom;
  vector<GNodeID> v_NodeLeft;
  vector<GNodeID> v_NodeRight;

  // === Graph information
  // Get graph id
  readGraphID( v_GraphIdToken, suppGraphID, suppClassName, data->get_line_nb() );

  // === Edge Informations
  // Read line of type "n"
  data->getLine( &v_Tokens );
  nbEdges    = readNbEdges( v_Tokens, data->get_line_nb() );
  lineNbEdge = data->get_line_nb();

  // === Dim information
  data->getLine( &v_Tokens );
  readDimension( v_Tokens, data->get_line_nb() );

  // === Node management
  // Read line of type "v"
  data->getLine( &v_Tokens );
  if ( v_Tokens.at(0) != "v" )
    error::structError( "v", "n", data->get_line_nb() );
  while ( !data->finished() && v_Tokens.at(0) == "v" )
  {
    if ( ( suppTCoord >= PARAM.INITIALTIME )
         && ( ( PARAM.NBTIMETOREAD != -1 && suppTCoord < PARAM.NBTIMETOREAD )
              || PARAM.NBTIMETOREAD == -1 ) )
    {
      // Parse vertex line
      readVertex( v_Tokens, suppXCoord, suppYCoord, suppTCoord, data->get_line_nb() );
      // If torgraph is set, manage edge node
      if ( PARAM.TORGRAPH )
      {
        if ( v_Nodes.back().yCoord == 0 )
          v_NodeBottom.push_back( v_Nodes.size() - 1  );
        else if ( ( v_Nodes.back().yCoord % ( dim.ydim - 1 ) ) == 0 )
          v_NodeTop.push_back( v_Nodes.size() - 1  );
        if ( v_Nodes.back().xCoord == 0 )
          v_NodeLeft.push_back( v_Nodes.size() - 1  );
        else if ( ( v_Nodes.back().xCoord % ( dim.xdim -1 ) ) == 0 )
          v_NodeRight.push_back( v_Nodes.size() - 1  );
      }
    }
    // Update supposed coord of following node
    suppXCoord++;
    if ( suppXCoord == dim.xdim )
    {
      suppXCoord = 0;
      suppYCoord++;
    }
    if ( suppYCoord == dim.ydim )
    {
      suppXCoord = 0;
      suppYCoord = 0;
      suppTCoord++;
    }
    // Read next line
    data->getLine( &v_Tokens );
  }
  if ( PARAM.NBTIMETOREAD != -1 || PARAM.INITIALTIME != 0 )
  {
    if ( PARAM.INITIALTIME != 0  && PARAM.NBTIMETOREAD == -1 )
      dim.tdim = dim.tdim - PARAM.INITIALTIME;
    else
      dim.tdim = PARAM.NBTIMETOREAD - PARAM.INITIALTIME;
  }

  // === Edge management
  // Read line of type "e"
  if ( PARAM.NOEDGELABEL )
  {
    currentEdgeId = computeEdge();
    while ( !data->finished() && v_Tokens.at(0) == "e" )
      data->getLine( &v_Tokens );
  }
  else
  {
    if ( v_Tokens.at(0) != "e" )
      error::structError( "e", "v", data->get_line_nb() );

    while ( !data->finished() && v_Tokens.at(0) == "e" )
    {
      readEdge( v_Tokens, currentEdgeId, data->get_line_nb() );
      // Read next line
      data->getLine( &v_Tokens );
    }
  }

  if ( (!PARAM.NOEDGELABEL) && (( PARAM.NBTIMETOREAD != -1 && suppTCoord < PARAM.NBTIMETOREAD ) || PARAM.NBTIMETOREAD == -1 ))
    if ( currentEdgeId != nbEdges )
      error::nbEdgeError( currentEdgeId, nbEdges, lineNbEdge );

  if ( PARAM.TORGRAPH )
  {
    for ( uint iNode = 0; iNode < v_NodeTop.size(); iNode++ )
      addEdgeTor( v_NodeBottom.at(iNode), v_NodeTop.at(iNode), 0 );
    for ( uint iNode = 0; iNode < v_NodeLeft.size(); iNode++ )
      addEdgeTor( v_NodeLeft.at(iNode), v_NodeRight.at(iNode), 0 );
  }

  // As we reach the end of the graph, we follow on the next one.
  if ( !data->finished() && v_Tokens[0] != "t" )
    error::structError( "t", "e", data->get_line_nb() );

  v_GraphIdToken = v_Tokens;
}
// End of GGraph::read( GReader *data, vector<string> &v_GraphIdToken )

void GGraph::addNode( GNodeLabel nodeLabel ,
                      int xcoord,
                      int ycoord,
                      int tcoord )
{
  /*
   * @brief addNode
   * Check if no other node share same coordinate, if so, print an error and
   * exit, else add node at the end of the list of node v_Nodes
   * @param nodeLabel: label of node added
   */
  // Create graph node with label
  vector<GGraphNode>::iterator it;
  GGraphNode graphNode = GGraphNode();
  graphNode.label  = nodeLabel;
  graphNode.xCoord = xcoord;
  graphNode.yCoord = ycoord;
  graphNode.tCoord = tcoord;
  // Check if there is no other node with same coord
  it = find( v_Nodes.begin(), v_Nodes.end(), graphNode );
  if ( it != v_Nodes.end() )
    error::addNodeError( this, graphNode, it - v_Nodes.begin() );

  // Add the node to the list v_Nodes
  v_Nodes.push_back( graphNode );
}
// End of GGraph::addNode( GNodeLabel nodeLabel ,
//                         uint xcoord,
//                         uint ycoord,
//                         uint tcoord )

void GGraph::deleteNode()
{
  /*
   * @brief deleteNode
   * Remove the last node in v_Nodes
   */
  v_Nodes.pop_back();
}
// End of GGraph::deleteNode()

void GGraph::setNodeToNoLabel( uint i )
{
  /*
   * @brief setNodeToNoLabel
   * Set node with ID 'pos' label to GNOLABEL and remove edge that link to
   * that node.
   * NB : It aim to 'remove' this node from the graph to not be considered
   * during mining.
   * @param pos : position of node to 'remove' from v_Nodes
   */
  if ( !PARAM.KEEPHOLE )
  {
    vector<GGraphEdge> v_tmpEdges;
    v_tmpEdges = v_Nodes.at(i).v_Edges;
    for ( uint idxEdge = 0; idxEdge < v_tmpEdges.size(); idxEdge++ )
    {
      if ( v_tmpEdges.at( idxEdge ).destNodeId != GNONODEID
           &&  v_Nodes.at( v_tmpEdges.at( idxEdge ).destNodeId ).label != GNOLABEL )
        deleteEdge( i, v_tmpEdges.at( idxEdge ).destNodeId );
    }
  }
  else
  {
    v_Nodes.at( i ).label = PARAM.MAXLBL;
  }
}
// End of GGraph::setNodeToNoLabel( uint i )

void GGraph::setNodeToNoLabel( vector<GNodeID> v_NodeId )
{
  /*
   * @brief setNodeToNoLabel
   * Set node with ID 'pos' label to GNOLABEL and remove edge that link to
   * that node.
   * NB : It aim to 'remove' this node from the graph to not be considered
   * during mining.
   * @param pos : position of node to 'remove' from v_Nodes
   */
  if ( !PARAM.KEEPHOLE )
  {
    vector<GGraphEdge> v_tmpEdges;
    for ( uint i = 0; i < v_NodeId.size(); i++ )
    {
      v_tmpEdges = v_Nodes.at(v_NodeId.at(i)).v_Edges;

      for ( uint idxEdge = 0; idxEdge < v_tmpEdges.size(); idxEdge++ )
      {
        if ( v_tmpEdges.at( idxEdge ).destNodeId != GNONODEID
             &&  v_Nodes.at( v_tmpEdges.at( idxEdge ).destNodeId ).label != GNOLABEL )
          deleteEdge( v_NodeId.at(i), v_tmpEdges.at( idxEdge ).destNodeId );
      }
    }
  }
  else
  {
    for ( uint i =0; i < v_NodeId.size(); i++ )
    {
      v_Nodes.at( v_NodeId.at(i) ).label = PARAM.MAXLBL;
    }
  }
}
// End of GGraph::setNodeToNoLabel( uint i )

int GGraph::computeEdge()
{
  int nbEdge = 0;
  vector<GNodeID> v_Neighbor;
  for ( uint iNode = 0; iNode < v_Nodes.size(); iNode++ )
  {
    v_Neighbor = computeNeighbor( iNode );
    for ( uint iNodeDest = 0; iNodeDest < v_Neighbor.size(); iNodeDest++ )
    {
      if ( !edgeExist( iNode, v_Neighbor.at(iNodeDest) ) )
      {
        addEdge( iNode, v_Neighbor.at(iNodeDest), 0 );
        nbEdge++;
      }
    }
  }
  return nbEdge;
}
// End of int GGraph::computeEdge()

void GGraph::addEdge( GNodeID    nodeIdFrom,
                      GNodeID    nodeIdDest,
                      GEdgeLabel edgeLabel,
                      int        line )
{
  /*
   * @brief addEdge
   * Check if nodeFrom and nodeDest are neighbord, if not, print an error
   * and exit, else compute position of edge in the list of edge in both
   * nodeFrom and nodeDest and add them at previously computed position
   * @param nodeIdFrom : ID of the node from which start the edges
   * @param nodeIdDest : ID of the node from which end the edges
   * @param edgeLabel  : Label of the edge
   * @param line       : line where edge is declared false (OPTIONNAL, use while
   *    reading graph file in read function )
   */
  int positionFrom = 0;
  int positionDest = 0;
  int deltaX = v_Nodes.at(nodeIdDest).xCoord - v_Nodes.at(nodeIdFrom).xCoord;
  int deltaY = v_Nodes.at(nodeIdDest).yCoord - v_Nodes.at(nodeIdFrom).yCoord;
  int deltaT = v_Nodes.at(nodeIdDest).tCoord - v_Nodes.at(nodeIdFrom).tCoord;

  // Check neighborhood condition
  if ( abs( deltaX ) + abs ( deltaY ) + abs( deltaT ) > 1 )
    error::addEdgeNotNeighborError( nodeIdFrom, nodeIdDest, line );

  // Set position of edge in list v_Edges
  if ( deltaX == 1 )
  {
    positionFrom = 2;
    positionDest = 0;
  }
  else if ( deltaX == -1 )
  {
    positionFrom = 0;
    positionDest = 2;
  }
  else if ( deltaY == 1 )
  {
    positionFrom = 1;
    positionDest = 3;
  }
  else if ( deltaY == -1 )
  {
    positionFrom = 3;
    positionDest = 1;
  }
  else if ( deltaT == 1 )
  {
    positionFrom = 4;
    positionDest = 5;
  }
  else if ( deltaT == -1 )
  {
    positionFrom = 5;
    positionDest = 4;
  }

  // Add edge to node From and node Dest
  v_Nodes.at(nodeIdFrom).v_Edges.at(positionFrom).label      = edgeLabel;
  v_Nodes.at(nodeIdFrom).v_Edges.at(positionFrom).destNodeId = nodeIdDest;
  v_Nodes.at(nodeIdFrom).v_Edges.at(positionFrom).eqEdgeId   = positionDest;

  v_Nodes.at(nodeIdDest).v_Edges.at(positionDest).label      = edgeLabel;
  v_Nodes.at(nodeIdDest).v_Edges.at(positionDest).destNodeId = nodeIdFrom;
  v_Nodes.at(nodeIdDest).v_Edges.at(positionDest).eqEdgeId   = positionFrom;
}
// End of GGraph::addEdge( GNodeID    nodeIdFrom,
//                         GNodeID    nodeIdDest,
//                         GEdgeLabel edgeLabel )

void GGraph::addEdgeTor( GNodeID    nodeIdFrom,
                         GNodeID    nodeIdDest,
                         GEdgeLabel edgeLabel )
{
  /*
   * @brief addEdge
   * Check if nodeFrom and nodeDest are neighbord, if not, print an error
   * and exit, else compute position of edge in the list of edge in both
   * nodeFrom and nodeDest and add them at previously computed position
   * @param nodeIdFrom : ID of the node from which start the edges
   * @param nodeIdDest : ID of the node from which end the edges
   * @param edgeLabel  : Label of the edge
   * @param line       : line where edge is declared false (OPTIONNAL, use while
   *    reading graph file in read function )
   */
  int positionFrom = 0;
  int positionDest = 0;
  int deltaX = v_Nodes.at(nodeIdDest).xCoord - v_Nodes.at(nodeIdFrom).xCoord;
  int deltaY = v_Nodes.at(nodeIdDest).yCoord - v_Nodes.at(nodeIdFrom).yCoord;
  int deltaT = v_Nodes.at(nodeIdDest).tCoord - v_Nodes.at(nodeIdFrom).tCoord;

  if ( deltaT == 1  )
  {
    positionFrom = 4;
    positionDest = 5;
  }
  else if ( deltaT == -1 )
  {
    positionFrom = 5;
    positionDest = 4;
  }
  else if ( deltaX > 0 )
  {
    positionFrom = 0;
    positionDest = 2;
  }
  else if ( deltaX < 0 )
  {
    positionFrom = 2;
    positionDest = 0;
  }
  else if ( deltaY > 0  )
  {
    positionFrom = 3;
    positionDest = 1;
  }
  else if ( deltaY < 0  )
  {
    positionFrom = 1;
    positionDest = 3;
  }
  else
  {
    cerr << "ERROR - Impossible to add torEdge with GNOANGLE" << endl;
    exit( EXIT_FAILURE );
  }


  // Add edge to node From and node Dest
  v_Nodes.at(nodeIdFrom).v_Edges.at(positionFrom).label      = edgeLabel;
  v_Nodes.at(nodeIdFrom).v_Edges.at(positionFrom).destNodeId = nodeIdDest;
  v_Nodes.at(nodeIdFrom).v_Edges.at(positionFrom).eqEdgeId   = positionDest;

  v_Nodes.at(nodeIdDest).v_Edges.at(positionDest).label      = edgeLabel;
  v_Nodes.at(nodeIdDest).v_Edges.at(positionDest).destNodeId = nodeIdFrom;
  v_Nodes.at(nodeIdDest).v_Edges.at(positionDest).eqEdgeId   = positionFrom;
}
// End of GGraph::addEdge( GNodeID    nodeIdFrom,
//                         GNodeID    nodeIdDest,
//                         GEdgeLabel edgeLabel )

void GGraph::deleteEdge( GNodeID nodeIdFrom, GNodeID nodeIdDest )
{
  /*
   * @brief deleteEdge
   * Find if edge between nodeIdFrom and nodeIdDest exist. If not, print an
   * error and exit, else, compute position of this edge in the list v_Edges
   * of both nodeIdFrom and nodeDest and set it de GNONODEID.
   * @param nodeIdFrom : Node from which edge to delete start
   * @param nodeIdDest : Node from which edge to delete end
   */

  int  pos      = 0 ;
  int  equivPos = 0 ;

  if ( !edgeExist( nodeIdFrom, nodeIdDest, pos ) )
    error::removeEdgeNotExistError( nodeIdFrom, nodeIdDest );

  // Set equivalent position depending on if edge is in x,y plane or not
  if ( pos < 4 )
    // If edge in x,y plane
    equivPos = ( pos + 2 ) % 4;
  else if ( pos == 4 )
    equivPos = pos + 1;
  else
    equivPos= pos - 1 ;

  // Update edge for node from
  v_Nodes.at(nodeIdFrom).v_Edges.at(pos).destNodeId = GNONODEID;
  v_Nodes.at(nodeIdFrom).v_Edges.at(pos).eqEdgeId   = GNONODEID;
  v_Nodes.at(nodeIdFrom).v_Edges.at(pos).label      = GNONODEID;

  // Update equivalent edge on node destination
  v_Nodes.at(nodeIdDest).v_Edges.at(equivPos).destNodeId = GNONODEID;
  v_Nodes.at(nodeIdDest).v_Edges.at(equivPos).eqEdgeId   = GNONODEID;
  v_Nodes.at(nodeIdDest).v_Edges.at(equivPos).label      = GNONODEID;
}
// End of GGraph::deleteEdge( GNodeID nodeIdFrom, GNodeID nodeIdDest )


bool GGraph::edgeExist( GNodeID nodeIdFrom, GNodeID nodeIdDest )
{
  /*
   * @brief edgeExist
   * Search if there is an edge between nodeIdFrom and nodeIdDest, return true
   * if such edge found, else return false.
   * @param nodeIdFrom : Node from which edge to find start
   * @param nodeIdDest : Node from which edge to find end
   * @return FALSE if edge is not equal to GNOEDGEID, else return TRUE
   */
  int nbEdges = v_Nodes.at(nodeIdFrom).v_Edges.size();
  int pos = 0;
  while ( pos != nbEdges )
  {
    // While edge not found or end of edge list not reached
    if ( v_Nodes.at(nodeIdFrom).v_Edges.at(pos).destNodeId == nodeIdDest )
      return true;
    else
      pos++;
  }
  return false;
}
// End of bool GGraph::edgeExist( GNodeID nodeIdFrom, GNodeID nodeIdDest )

bool GGraph::edgeExist( GNodeID nodeIdFrom, GNodeID nodeIdDest, int &pos )
{
  /*
   * @brief edgeExist
   * Search if there is an edge between nodeIdFrom and nodeIdDest, return true
   * if such edge found, else return false.
   * @param nodeIdFrom : Node from which edge to find start
   * @param nodeIdDest : Node from which edge to find end
   * @param pos        : Variable to know position of edge if found
   * @return FALSE if edge is not equal to GNOEDGEID, else return TRUE
   */
  int nbEdges = v_Nodes.at(nodeIdFrom).v_Edges.size();
  pos = 0;
  while ( pos != nbEdges )
  {
    // While edge not found or end of edge list not reached
    if ( v_Nodes.at(nodeIdFrom).v_Edges.at(pos).destNodeId == nodeIdDest )
      return true;
    else
      pos++;
  }
  return false;
}
// End of bool GGraph::edgeExist( GNodeID nodeIdFrom, GNodeID nodeIdDest, int &pos )

vector<GNodeID> GGraph::processHole()
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
  if ( PARAM.DEBUG_MODE )
    cout << "- Process hole in graphs." << endl;

  uint nbNeigh;
  vector<GNodeID> v_IdNodeToDelete;
  vector<GNodeID> v_Neighbor;

  for ( uint iNode = 0; iNode < v_Nodes.size(); iNode++ )
  {
    // No need to start from node 0 and to finish at node.size()
    nbNeigh  = 0;
    v_Neighbor = computeNeighborHole( iNode );
    for ( uint iNeigh = 0; iNeigh < v_Neighbor.size(); iNeigh++ )
      if ( v_Neighbor.at(iNeigh) < v_Nodes.size()
           && v_Nodes.at(iNode).label == v_Nodes.at( v_Neighbor.at(iNeigh) ).label )
        nbNeigh++;
    if ( PARAM.DIRECT_HOLEMINING )
    {
      if ( dim.tdim == 1 && nbNeigh == 5 )
        v_IdNodeToDelete.push_back( iNode ); // Only 1 frame (i.e. 2D grid)
      else if ( dim.tdim != 1 )
      {
        if ( ( v_Nodes.at(iNode).tCoord == PARAM.INITIALTIME ) && nbNeigh == 6 )
          v_IdNodeToDelete.push_back( iNode ); // First Frame
        if ( ( v_Nodes.at(iNode).tCoord  == PARAM.NBTIMETOREAD - 1 ) && nbNeigh == 6 )
          v_IdNodeToDelete.push_back( iNode ); // Last frame
        else if ( ( ( v_Nodes.at(iNode).tCoord  > PARAM.INITIALTIME )
                    || ( v_Nodes.at(iNode).tCoord  < PARAM.NBTIMETOREAD ) )
                  && nbNeigh == 7 )
          v_IdNodeToDelete.push_back( iNode );
        else if ( ( v_Nodes.at(iNode).tCoord  < PARAM.INITIALTIME )
                  || ( v_Nodes.at(iNode).tCoord  > PARAM.NBTIMETOREAD ) )
        {
          cout << "ERROR - Node " << iNode << " does not have valid temporal coordinate !" << endl;
          exit( EXIT_FAILURE );
        }
      }
    }
    else if ( PARAM.INDIRECT_HOLEMINING )
    {
      if ( dim.tdim == 1 && nbNeigh == 9 )
        v_IdNodeToDelete.push_back( iNode );
      else if ( dim.tdim != 1 )
      {
        if  ( ( ( v_Nodes.at(iNode).tCoord  == PARAM.INITIALTIME )
                || ( v_Nodes.at(iNode).tCoord  == PARAM.NBTIMETOREAD ) )
              && nbNeigh == v_Neighbor.size() )
          v_IdNodeToDelete.push_back( iNode );
        else if ( ( ( v_Nodes.at(iNode).tCoord  > PARAM.INITIALTIME )
                    || ( v_Nodes.at(iNode).tCoord  < PARAM.NBTIMETOREAD ) )
                  && nbNeigh == v_Neighbor.size() )
          v_IdNodeToDelete.push_back( iNode );
        else if ( ( v_Nodes.at(iNode).tCoord  < PARAM.INITIALTIME )
                  || ( v_Nodes.at(iNode).tCoord  > PARAM.NBTIMETOREAD ) )
        {
          cout << "ERROR - Node " << iNode << " does not have valid temporal coordinate !" << endl;
          exit( EXIT_FAILURE );
        }
      }
    }
  }

  return v_IdNodeToDelete;

  //  for ( uint idxNode = 0; idxNode < v_IdNodeToDelete.size(); idxNode++ )
  //this->setNodeToNoLabel( v_IdNodeToDelete.at( idxNode ) );
}
// End of GDatabase::processHoleInGrid( GClassDB* pClassDB )

vector<GNodeID> GGraph::processFrameHole()
{

  if ( PARAM.DEBUG_MODE )
    cout << "- Process frame hole in graphs." << endl;

  vector<GNodeID> v_IdNodeToDelete;
  for ( uint i = 0; i < uint(dim.xdim*dim.ydim); i++ )
  {
    for ( uint j = 1 ; j < uint(dim.tdim); j++ )
    {
      if ( v_Nodes.at( i + ( (j-1)*dim.xdim*dim.ydim ) ).label
           == v_Nodes.at( i+(j*dim.xdim*dim.ydim) ).label
           && v_Nodes.at( i + ((j-1)*dim.xdim*dim.ydim) ).label
           == v_Nodes.at( i+j*dim.xdim*dim.ydim ).label )
        v_IdNodeToDelete.push_back( i+j*dim.xdim*dim.ydim );
    }
  }
  return v_IdNodeToDelete;
}

uint GGraph::getEdgePositionFrom( GNodeID nodeFrom, GNodeID nodeDest )
{
  /*
   * TODO : RD
   * Copy Desc
   */
  uint iEdge = 0;
  while ( iEdge < v_Nodes.at(nodeFrom).v_Edges.size() )
  {
    if ( v_Nodes.at(nodeFrom).v_Edges.at(iEdge).destNodeId == nodeDest )
      return iEdge;
    iEdge++;
  }

  return GNOANGLE;
}
// End of GGraph::getEdgePositionFrom( GNodeID nodeFrom, GNodeID nodeDest )

// Private CONSTANTS _________________________________________________________//
// Private Methods ___________________________________________________________//
void GGraph::readGraphID( vector<string> &v_Tokens,
                          uint           &suppGraphID,
                          string         suppClassName,
                          int            line )
{
  /*
   * TODO : RD
   * Copy Desc
   */
  if ( v_Tokens[0] != "t" )
    error::structError( "t", line );
  else if ( !isUInt( v_Tokens.at(1) ) )
    error::structValue_UINT_Error( "graph ID", line );

  graphID = stoul( v_Tokens.at(1) );

  if ( v_Tokens.size() == 3 )
    className = v_Tokens.at(2);
  else
    className = "";

  if ( PARAM.DEBUG_MODE )
    cout << "Reading graph " << suppGraphID << " of class " << className << endl;

  if ( className != suppClassName )
    suppGraphID = 0;

  if ( graphID != suppGraphID )
    error::graphIDError( graphID, suppGraphID, line );
  else
    suppGraphID++;


}
// End of GGraph::readGraphID( vector<string> &v_Tokens,
//                             uint           &suppGraphID,
//                             string         suppClassName,
//                             int            line )

void GGraph::readDimension( vector<string> &v_Tokens, int line )
{
  /*
   * @brief readDimension
   * Check if line start with right char and if value of dimensiont are unsigned
   * int, if not write error and exit, else set them
   * @param v_Tokens : String token to parse
   * @param line     : Line actually read by GReader
   */
  if ( v_Tokens.at(0) != "dim" )
    error::structError( "dim", "t", line );
  else
  {
    if ( v_Tokens.size() == 3 )
    {
      if ( !isUInt( v_Tokens.at(1) )
           || !isUInt( v_Tokens.at(2) ) )
        error::structValue_UINT_Error( "dimensions", line );
      else
      {
        dim.xdim = stoul( v_Tokens.at(1) );
        dim.ydim = stoul( v_Tokens.at(2) );
        dim.tdim = 1;
      }
    }
    else
    {

      if ( !isUInt( v_Tokens.at(1) )
           || !isUInt( v_Tokens.at(2) )
           || !isUInt( v_Tokens.at(3) ) )
        error::structValue_UINT_Error( "dimensions", line );
      else
      {
        dim.xdim = stoul( v_Tokens.at(1) );
        dim.ydim = stoul( v_Tokens.at(2) );
        dim.tdim = stoul( v_Tokens.at(3) );
      }
    }
  }

  if ( dim.xdim == 0 || dim.ydim == 0 || dim.tdim == 0 )
    error::structValue_PosUINT_Error( "dimensions", line );
}
// End of GGraph::readDimension( vector<string> &v_Tokens, int line )

GEdgeID GGraph::readNbEdges( vector<string> &v_Tokens, int line )
{
  /*
   * @brief readNbEdges
   * Check if line start with right car and if value of number of edge is
   * unsigned int, if not write error and exit, else return the number of edge
   * set in the file
   * @param v_Tokens : String token to parse
   * @param line     : Line actually read by GReader
   * @return The number of edge
   */
  if ( v_Tokens.at(0) != "n" )
    error::structError( "n", "dim", line );
  else if ( !isUInt( v_Tokens.at(1) ) )
    error::structValue_UINT_Error( "edge number", line );

  return stoul( v_Tokens.at(1) );
}
// End of GGraph::readNbEdges( vector<string> &v_Tokens, int line )

void GGraph::readVertex( vector<string> &v_Tokens,
                         int &suppXCoord,
                         int &suppYCoord,
                         int &suppTCoord,
                         int line )
{
  /*
   * @brief readVertex
   * Check if line that start with 'v', print error if line does not respect
   * specified format (like all coord must be listed, value must be unsigned
   * integer or coord must be in the write order).
   * @param v_Tokens   : Tokens of the string to be parsed
   * @param suppXCoord : Supposed coordinate along X of the node parsed
   * @param suppYCoord : Supposed coordinate along Y of the node parsed
   * @param suppTCoord : Supposed coordinate along T of the node parsed
   * @param line       : Line actually read by GReader
   */
  if ( v_Tokens[0] != "v" )
    error::structError( "v", line );

  if ( stoi( v_Tokens.at(2) ) > PARAM.MAXLBL )
    PARAM.MAXLBL = stoi( v_Tokens.at(2) );

  addNode( stoul( v_Tokens.at(2) ), suppXCoord, suppYCoord, suppTCoord );
}
// End of GGraph::readVertex( vector<string> &v_Tokens,
//                            uint           &suppXCoord,
//                            uint           &suppYCoord,
//                            uint           &suppTCoord,
//                            int            line )

void GGraph::readEdge( vector<string> &v_Tokens,
                       uint           &currentEdgeId,
                       int            line )
{
  /*
   * @brief readEdge
   * Check if line that start with 'e', print error if line does not respect
   * specified format
   * @param v_Tokens      : Tokens of the string to be parsed
   * @param currentEdgeId : Currently computed/read edgeID
   * @param line          : Line actually read by GReader
   */
  uint       idx;
  // Create edge needed variables
  GEdgeID    edgeId;
  GNodeID    nodeIdFrom;
  GNodeID    nodeIdDest;
  GEdgeLabel edgeLabel;

  if ( v_Tokens[0] != "e" )
    error::structError( "e", line  ) ;
  else if ( v_Tokens.size() < 4 )
    error::edgeMissingValError( line );
  else if ( !isUInt( v_Tokens.at(1) )
            || !isUInt( v_Tokens.at(2) )
            || !isUInt( v_Tokens.at(3) ) )
    error::edgeValue_UINT_Error( line );

  nodeIdFrom = stoul( v_Tokens.at(1) );
  nodeIdDest = stoul( v_Tokens.at(2) );
  edgeLabel  = stoul( v_Tokens.at(3) );
  // Part of edge management that handle if input file has edgeID or not
  if ( v_Tokens.size() == 5 )
  {
    if ( !isUInt( v_Tokens.at(4) ) )
      error::edgeValue_UINT_Error( line );

    // If edge ID given : Set edgeID
    edgeId = stoul( v_Tokens.at(4) );
    if ( edgeId > currentEdgeId )
      // If new edge : Update currentEdgeId counter
      currentEdgeId = edgeId++;
  }
  else
  {
    // Find if edge already exists
    idx = 0;
    while ( idx < v_Nodes.at(nodeIdDest).v_Edges.size()
            && ( v_Nodes.at(nodeIdDest).v_Edges.at(idx).destNodeId != nodeIdFrom
                 || v_Nodes.at(nodeIdDest).v_Edges.at(idx).label != edgeLabel ) )
      idx++;
    if ( idx == v_Nodes.at(nodeIdDest).v_Edges.size() )
      // If edge does not exist yet : increment current edge id
      currentEdgeId++;
  }
  addEdge( nodeIdFrom, nodeIdDest, edgeLabel, line );
}
// End of  GGraph::readEdge( vector<string> &v_Tokens,
//                           uint           &currentEdgeId,
//                           int            line )

vector<GNodeID> GGraph::computeNeighbor( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor;
  GNodeID nodeLeft      = iNode - 1;
  GNodeID nodeRight     = iNode + 1;
  GNodeID nodeTop       = iNode + dim.xdim;
  GNodeID nodeBottom    = iNode - dim.xdim;
  GNodeID nodeForward   = iNode + ( dim.xdim * dim.ydim );
  GNodeID nodeBackward  = iNode - ( dim.xdim * dim.ydim );;

  if ( nodeLeft > iNode && nodeLeft < v_Nodes.size()
       && areNeighbor( v_Nodes.at(iNode), v_Nodes.at(nodeLeft) ) )
    v_Neighbor.push_back(nodeLeft);
  if ( nodeRight > iNode && nodeRight< v_Nodes.size()
       && areNeighbor( v_Nodes.at(iNode), v_Nodes.at(nodeRight) ) )
    v_Neighbor.push_back(nodeRight);
  if ( nodeTop > iNode && nodeTop < v_Nodes.size()
       && areNeighbor( v_Nodes.at(iNode), v_Nodes.at(nodeTop) ) )
    v_Neighbor.push_back(nodeTop);
  if ( nodeBottom > iNode && nodeBottom < v_Nodes.size()
       && areNeighbor( v_Nodes.at(iNode), v_Nodes.at(nodeBottom) ) )
    v_Neighbor.push_back(nodeBottom);
  if ( nodeForward > iNode && nodeForward < v_Nodes.size()
       && areNeighbor( v_Nodes.at(iNode), v_Nodes.at(nodeForward) ) )
    v_Neighbor.push_back(nodeForward);
  if ( nodeBackward > iNode && nodeBackward < v_Nodes.size()
       && areNeighbor( v_Nodes.at(iNode), v_Nodes.at(nodeBackward) ) )
    v_Neighbor.push_back(nodeBackward);

  return v_Neighbor;
}

vector<GNodeID> GGraph::computeNeighborHole( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor;
  if ( PARAM.DIRECT_HOLEMINING )
    v_Neighbor = computeDirectNeighbor( iNode );
  else if ( PARAM.INDIRECT_HOLEMINING )
    v_Neighbor = computeIndirectNeighbor( iNode );
  return v_Neighbor;
}

vector<GNodeID> GGraph::computeDirectNeighbor( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor;
  if ( dim.tdim != 1 )
    v_Neighbor = getDirectNeighbor2Dt( iNode );
  else
    v_Neighbor = getDirectNeighbor2D( iNode );
  return v_Neighbor;
}

vector<GNodeID> GGraph::computeIndirectNeighbor( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor;
  if ( dim.tdim != 1 )
    v_Neighbor = getIndirectNeighbor2Dt( iNode );
  else
    v_Neighbor = getIndirectNeighbor2D( iNode );
  return v_Neighbor;
}

vector<GNodeID> GGraph::getDirectNeighbor2D( GNodeID iNode )
{
  if ( v_Nodes.at(iNode).xCoord == 0
       && v_Nodes.at(iNode).yCoord == 0 )
    return getDirectNeighbor2D_BL( iNode ); // Bottom Left
  else if ( v_Nodes.at(iNode).xCoord == dim.xdim-1
            && v_Nodes.at(iNode).yCoord == 0 )
    return getDirectNeighbor2D_BR( iNode ); // Bottom right
  else if ( v_Nodes.at(iNode).xCoord == 0
            && v_Nodes.at(iNode).yCoord == dim.ydim-1 )
    return getDirectNeighbor2D_TL( iNode ); // Top left
  else if ( v_Nodes.at(iNode).xCoord == dim.xdim-1
            && v_Nodes.at(iNode).yCoord == dim.ydim-1 )
    return getDirectNeighbor2D_TR( iNode ); // Top right
  else if ( v_Nodes.at(iNode).yCoord == 0 )
    return getDirectNeighbor2D_B( iNode );  // Bottom
  else if ( v_Nodes.at(iNode).xCoord == dim.xdim-1 )
    return getDirectNeighbor2D_R( iNode );  // Right
  else if ( v_Nodes.at(iNode).yCoord == dim.ydim-1 )
    return getDirectNeighbor2D_T( iNode );  // Top
  else if ( v_Nodes.at(iNode).xCoord == 0 )
    return getDirectNeighbor2D_L( iNode );  // Left
  else
    return getDirectNeighbor2D_C( iNode );  // Center
}

vector<GNodeID> GGraph::getDirectNeighbor2Dt( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor;
  GNodeID currNode;
  currNode = iNode;

  v_Neighbor = getDirectNeighbor2D( currNode );
  if ( v_Nodes.at(iNode).tCoord == PARAM.INITIALTIME )
    v_Neighbor.push_back( iNode + ( dim.xdim * dim.ydim ) ); // Node Forward
  else if ( v_Nodes.at(iNode).tCoord == PARAM.NBTIMETOREAD-1 )
    v_Neighbor.push_back( iNode - ( dim.xdim * dim.ydim ) ); // Node Backward
  else
  {
    v_Neighbor.push_back( iNode + ( dim.xdim * dim.ydim ) ); // Node Forward
    v_Neighbor.push_back( iNode - ( dim.xdim * dim.ydim ) ); // Node Backward
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getIndirectNeighbor2D( GNodeID iNode )
{
  if ( v_Nodes.at(iNode).xCoord == 0
       && v_Nodes.at(iNode).yCoord == 0 )
    return getIndirectNeighbor2D_BL( iNode ); // Bottom Left
  else if ( v_Nodes.at(iNode).xCoord == dim.xdim-1
            && v_Nodes.at(iNode).yCoord == 0 )
    return getIndirectNeighbor2D_BR( iNode ); // Bottom right
  else if ( v_Nodes.at(iNode).xCoord == 0
            && v_Nodes.at(iNode).yCoord == dim.ydim-1 )
    return getIndirectNeighbor2D_TL( iNode ); // Top left
  else if ( v_Nodes.at(iNode).xCoord == dim.xdim-1
            && v_Nodes.at(iNode).yCoord == dim.ydim-1 )
    return getIndirectNeighbor2D_TR( iNode ); // Top right
  else if ( v_Nodes.at(iNode).yCoord == 0 )
    return getIndirectNeighbor2D_B( iNode );  // Bottom
  else if ( v_Nodes.at(iNode).xCoord == dim.xdim-1 )
    return getIndirectNeighbor2D_R( iNode );  // Right
  else if ( v_Nodes.at(iNode).yCoord == dim.ydim-1 )
    return getIndirectNeighbor2D_T( iNode );  // Top
  else if ( v_Nodes.at(iNode).xCoord == 0 )
    return getIndirectNeighbor2D_L( iNode );  // Left
  else
    return getIndirectNeighbor2D_C( iNode );  // Center
}

vector<GNodeID> GGraph::getIndirectNeighbor2Dt( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor;
  vector<GNodeID> v_NeighborTime;
  GNodeID currNode;

  currNode = iNode;
  v_Neighbor = getIndirectNeighbor2D( currNode );

  if ( v_Nodes.at(iNode).tCoord == PARAM.INITIALTIME )
  {
    currNode = iNode + ( dim.xdim * dim.ydim ); // Node Forward
    v_NeighborTime = getIndirectNeighbor2D( currNode );
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTime.begin(), v_NeighborTime.end() );
  }
  else if ( v_Nodes.at(iNode).tCoord == PARAM.NBTIMETOREAD-1 )
  {
    currNode = iNode - ( dim.xdim * dim.ydim ); // Node Backward
    v_NeighborTime = getIndirectNeighbor2D( currNode );
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTime.begin(), v_NeighborTime.end() );
  }
  else
  {
    currNode = iNode + ( dim.xdim * dim.ydim ); // Node Forward
    v_NeighborTime = getIndirectNeighbor2D( currNode );
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTime.begin(), v_NeighborTime.end() );

    currNode = iNode - ( dim.xdim * dim.ydim ); // Node Backward
    v_NeighborTime = getIndirectNeighbor2D( currNode );
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTime.begin(), v_NeighborTime.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getIndirectNeighbor2D_C( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                           // This Node
    iNode - 1,                       // Left
    iNode + 1,                       // Right
    iNode + dim.xdim,                // Top
    iNode - dim.xdim,                // Bottom
    iNode + dim.xdim - 1,            // Top Left
    iNode + dim.xdim + 1,            // Top Right
    iNode - dim.xdim - 1,            // Bottom Left
    iNode - dim.xdim + 1,            // Bottom Right
  };
  return v_Neighbor;
}

vector<GNodeID> GGraph::getIndirectNeighbor2D_BL( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode + 1,            // Right
    iNode + dim.xdim,     // Top
    iNode + dim.xdim + 1, // Top Right
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode + dim.xdim -1    ,                 // Left
      iNode + 2*dim.xdim - 1 ,                 // Top Left
      iNode + ( dim.xdim * (dim.ydim-1) ),     // Bottom
      iNode + ( dim.xdim * (dim.ydim-1) ) + 1, // Bottom Right
      (GNodeID)dim.ydim * dim.xdim -1 ,        // Bottom Left
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getIndirectNeighbor2D_BR( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode - 1,            // Left
    iNode + dim.xdim,     // Top
    iNode + dim.xdim - 1, // Top Left
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      0,                                       // Right
      iNode + 1 ,                              // Top Right
      iNode + ( dim.xdim * (dim.ydim-1) ),     // Bottom
      (GNodeID)(dim.xdim * (dim.ydim -1)),     // Bottom Right
      iNode + ( dim.xdim * (dim.ydim-1) ) - 1, // Bottom Left
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getIndirectNeighbor2D_TL( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode + 1,            // Right
    iNode - dim.xdim,     // Bottom
    iNode - dim.xdim + 1, // Bottom Right
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode + dim.xdim - 1,                    // Left
      iNode - 1,                               // Bottom Left
      iNode - ( dim.xdim * (dim.ydim-1) ),     // Top
      iNode - ( dim.xdim * (dim.ydim-1) ) + 1, // Top Right
      (GNodeID)(dim.xdim - 1),                 // Top Left
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
  return v_Neighbor;
}

vector<GNodeID> GGraph::getIndirectNeighbor2D_TR( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode - 1,            // Left
    iNode - dim.xdim,     // Bottom
    iNode - dim.xdim - 1, // Bottom Left
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode - dim.xdim + 1,                    // Right
      iNode - 2*dim.xdim + 1,                  // Bottom Right
      iNode - ( dim.xdim * (dim.ydim-1) ),     // Top
      0,                                       // Top Right
      iNode - ( dim.xdim * (dim.ydim-1) ),     // Top Left
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getIndirectNeighbor2D_B( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode - 1,            // Left
    iNode + 1,            // Right
    iNode + dim.xdim,     // Top
    iNode + dim.xdim - 1, // Top Left
    iNode + dim.xdim + 1, // Top Right
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode + ( dim.xdim * (dim.ydim-1) ),     // Bottom
      iNode + ( dim.xdim * (dim.ydim-1) ) + 1, // Bottom Right
      iNode + ( dim.xdim * (dim.ydim-1) ) - 1, // Bottom Left
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getIndirectNeighbor2D_R( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode - 1,            // Left
    iNode + dim.xdim,     // Top
    iNode - dim.xdim,     // Bottom
    iNode + dim.xdim - 1, // Top Left
    iNode - dim.xdim - 1, // Bottom Left
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode - dim.xdim + 1,   // Right
      iNode - 2*dim.xdim + 1, // Bottom Right
      iNode + 1,              // Top Right
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }

  return v_Neighbor;
}

vector<GNodeID> GGraph::getIndirectNeighbor2D_T( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode - 1,            // Left
    iNode + 1,            // Right
    iNode - dim.xdim,     // Bottom
    iNode - dim.xdim - 1, // Bottom Left
    iNode - dim.xdim + 1, // Bottom Right
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode - ( dim.xdim * (dim.ydim-1) ),     // Top
      iNode - ( dim.xdim * (dim.ydim-1) ) + 1, // Top Right
      iNode - ( dim.xdim * (dim.ydim-1) ) - 1, // Top Left
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getIndirectNeighbor2D_L( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode + 1,            // Right
    iNode + dim.xdim,     // Top
    iNode - dim.xdim,     // Bottom
    iNode + dim.xdim + 1, // Top Right
    iNode - dim.xdim + 1, // Bottom Right
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode + dim.xdim - 1,   // Left
      iNode + 2*dim.xdim - 1, // Top Left
      iNode - 1,              // Bottom Left
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getDirectNeighbor2D_C( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                           // This Node
    iNode - 1,                       // Left
    iNode + 1,                       // Right
    iNode + dim.xdim,                // Top
    iNode - dim.xdim,                // Bottom
  };
  return v_Neighbor;
}

vector<GNodeID> GGraph::getDirectNeighbor2D_BL( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode + 1,            // Right
    iNode + dim.xdim,     // Top
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode + dim.xdim -1    ,                 // Left
      iNode + ( dim.xdim * (dim.ydim-1) ),     // Bottom
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getDirectNeighbor2D_BR( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode - 1,            // Left
    iNode + dim.xdim,     // Top
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      0,                                       // Right
      iNode + ( dim.xdim * (dim.ydim-1) ),     // Bottom
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getDirectNeighbor2D_TL( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode + 1,            // Right
    iNode - dim.xdim,     // Bottom
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode + dim.xdim - 1,                    // Left
      iNode - ( dim.xdim * (dim.ydim-1) ),     // Top
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
  return v_Neighbor;
}

vector<GNodeID> GGraph::getDirectNeighbor2D_TR( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode - 1,            // Left
    iNode - dim.xdim,     // Bottom
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode - dim.xdim + 1,                    // Right
      iNode - ( dim.xdim * (dim.ydim-1) ),     // Top
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getDirectNeighbor2D_B( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode - 1,            // Left
    iNode + 1,            // Right
    iNode + dim.xdim,     // Top
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode + ( dim.xdim * (dim.ydim-1) ),     // Bottom
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getDirectNeighbor2D_R( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode - 1,            // Left
    iNode + dim.xdim,     // Top
    iNode - dim.xdim,     // Bottom
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode - dim.xdim + 1,   // Right
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }

  return v_Neighbor;
}

vector<GNodeID> GGraph::getDirectNeighbor2D_T( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode - 1,            // Left
    iNode + 1,            // Right
    iNode - dim.xdim,     // Bottom
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode - ( dim.xdim * (dim.ydim-1) ),     // Top
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

vector<GNodeID> GGraph::getDirectNeighbor2D_L( GNodeID iNode )
{
  vector<GNodeID> v_Neighbor = {
    iNode,                // This Node
    iNode + 1,            // Right
    iNode + dim.xdim,     // Top
    iNode - dim.xdim,     // Bottom
  };
  if ( PARAM.TORGRAPH )
  {
    vector<GNodeID> v_NeighborTor = {
      iNode + dim.xdim - 1,   // Left
    };
    v_Neighbor.insert( v_Neighbor.end(), v_NeighborTor.begin(), v_NeighborTor.end() );
  }
  return v_Neighbor;
}

bool GGraph::areNeighbor( GGraphNode nodeFrom, GGraphNode nodeDest )
{
  /*
   * TODO : RD =================================================================
   * Copy description from header
   * ===========================================================================
   */
  int deltaX = nodeFrom.xCoord - nodeDest.xCoord;
  int deltaY = nodeFrom.yCoord - nodeDest.yCoord;
  int deltaT = nodeFrom.tCoord - nodeDest.tCoord;
  // Check neighborhood condition
  return ( abs( deltaX ) + abs ( deltaY ) + abs( deltaT ) == 1 );
}
// End of GridGenerator::areNeighbor( GGraphNode nodeFrom, GGraphNode nodeDest )

bool GGraph::areNeighborHole( GGraphNode nodeFrom, GGraphNode nodeDest )
{
  /*
   * TODO : RD =================================================================
   * Copy description from header
   * ===========================================================================
   */
  int deltaX = nodeFrom.xCoord - nodeDest.xCoord;
  int deltaY = nodeFrom.yCoord - nodeDest.yCoord;
  int deltaT = nodeFrom.tCoord - nodeDest.tCoord;

  if ( abs( deltaT ) <= 1 )
  {
    if ( nodeFrom.yCoord == 0 && nodeFrom.xCoord == 0 ) // Top left corner
    {
      if ( nodeDest.yCoord == -1 && nodeDest.xCoord == dim.xdim-1 ) // Bottom right corner
        return true;
      if ( nodeDest.yCoord == dim.ydim-1 && nodeDest.xCoord == 0 ) // Bottom left corner
        return true;
      if ( nodeDest.yCoord == dim.ydim-1 && nodeDest.xCoord == 1 ) // Right of bottom left corner
        return true;
      if ( nodeDest.yCoord == 0 && nodeDest.xCoord == dim.xdim-1 ) // Top right corner
        return true;
      if ( nodeDest.yCoord == 1 && nodeDest.xCoord == dim.xdim-1 )  // Below of Top right corner
        return true;
    }
    else if ( nodeFrom.yCoord == 0 && nodeFrom.xCoord == dim.xdim-1  ) // Top right corner
    {
      if ( nodeDest.yCoord == dim.ydim-1 && nodeDest.xCoord == dim.xdim-1 ) // Bottom right corner
        return true;
      if ( nodeDest.yCoord == dim.ydim-1 && nodeDest.xCoord == dim.xdim-2 ) // Left of Bottom right corner
        return true;
      if ( nodeDest.yCoord == dim.ydim-1 && nodeDest.xCoord == 0 ) // Bottom left corner
        return true;
      if ( nodeDest.yCoord == 0 && nodeDest.xCoord == 0 ) // Top left corner
        return true;
      if ( nodeDest.yCoord == 1 && nodeDest.xCoord == 0 )  // Belowof of Top left corner
        return true;
    }
    else if ( nodeFrom.yCoord == dim.ydim-1 && nodeFrom.xCoord == 0 ) // Bottom left corner
    {
      if ( nodeDest.yCoord == 0 && nodeDest.xCoord == dim.xdim-1 ) // Top right corner
        return true;
      if ( nodeDest.yCoord == dim.ydim-1 && nodeDest.xCoord == dim.xdim-1 ) // Bottom right corner
        return true;
      if ( nodeDest.yCoord == dim.ydim-2 && nodeDest.xCoord == dim.xdim-1 ) // Above Bottom right corner
        return true;
      if ( nodeDest.yCoord == 0 && nodeDest.xCoord == 0 ) // Top left corner
        return true;
      if ( nodeDest.yCoord == 0 && nodeDest.xCoord == 1 )  // Right of of Top left corner
        return true;
    }
    else if ( nodeFrom.yCoord == dim.ydim-1 && nodeFrom.xCoord == dim.xdim-1 ) // Bottom right cornder
    {
      if ( nodeDest.yCoord == 0 && nodeDest.xCoord == 0 ) // Top left corner
        return true;
      if ( nodeDest.yCoord == 0 && nodeDest.xCoord == dim.xdim-1 ) // Top right corner
        return true;
      if ( nodeDest.yCoord == 0 && nodeDest.xCoord == dim.xdim-2 ) // Left of Top right corner
        return true;
      if ( nodeDest.yCoord == dim.ydim-1 && nodeDest.xCoord == 0 ) // Bottom left corner
        return true;
      if ( nodeDest.yCoord == dim.ydim-2 && nodeDest.xCoord == 0 ) // Above of Bottom left corner
        return true;
    }
    else if ( nodeFrom.yCoord == 0  ) // Left line
    {
      if ( nodeDest.yCoord == dim.ydim-1 && nodeDest.xCoord == nodeFrom.xCoord-1 ) // Above of Bottom left corner
        return true;
      if ( nodeDest.yCoord == dim.ydim-1 && nodeDest.xCoord == nodeFrom.xCoord ) // Above of Bottom left corner
        return true;
      if ( nodeDest.yCoord == dim.ydim-1 && nodeDest.xCoord == nodeFrom.xCoord+1 ) // Above of Bottom left corner
        return true;
    }
    else if ( nodeFrom.yCoord == dim.ydim-1  ) // Right line
    {
      if ( nodeDest.yCoord == 0 && nodeDest.xCoord == nodeFrom.xCoord-1 ) // Above of Bottom left corner
        return true;
      if ( nodeDest.yCoord == 0 && nodeDest.xCoord == nodeFrom.xCoord ) // Above of Bottom left corner
        return true;
      if ( nodeDest.yCoord == 0 && nodeDest.xCoord == nodeFrom.xCoord+1 ) // Above of Bottom left corner
        return true;
    }
    else if ( nodeFrom.xCoord == 0 ) // Top line
    {
      if ( nodeDest.yCoord == nodeFrom.yCoord-1 && nodeDest.xCoord == dim.xdim-1 ) // Above of Bottom left corner
        return true;
      if ( nodeDest.yCoord == nodeFrom.yCoord && nodeDest.xCoord == dim.xdim-1 ) // Above of Bottom left corner
        return true;
      if ( nodeDest.yCoord == nodeFrom.yCoord+1 && nodeDest.xCoord == dim.xdim-1 ) // Above of Bottom left corner
        return true;
    }
    else if ( nodeFrom.xCoord == dim.xdim-1 ) // Bottom line
    {
      if ( nodeDest.yCoord == nodeFrom.yCoord-1 && nodeDest.xCoord == 0 ) // Above of Bottom left corner
        return true;
      if ( nodeDest.yCoord == nodeFrom.yCoord && nodeDest.xCoord == 0 ) // Above of Bottom left corner
        return true;
      if ( nodeDest.yCoord == nodeFrom.yCoord+1 && nodeDest.xCoord == 0 ) // Above of Bottom left corner
        return true;
    }
  }
  return ( abs( deltaX ) <= 1 &&  abs ( deltaY ) <= 1 && abs( deltaT ) <= 1 );
}
// End of GridGenerator::areNeighbor( GGraphNode nodeFrom, GGraphNode nodeDest )

bool GGraph::isUInt( const string &s )
{
  /*
   * @brief isInt
   * Test if string is an unsigned int
   * @param s : Input string
   * @return True if s is an uint false elsewhere
   */
  if( s.empty() || !isdigit( s[0] ) )
    return false ;

  char * p ;
  strtol( s.c_str(), &p, 10) ;

  return (*p == 0) ;
}
// End of GGraph::isUInt( const string &s )

//============================== OPERATOR OVERLOAD  ==========================//
ostream& operator<<( ostream& stream, GGraph &graph )
{
  /*
   * @brief operator <<
   * Overload operator << to pring graph (ie list of node and edges.
   * @param stream : output stream
   * @param graph : graph to print
   * @return
   */
  int nbEdgePlane =
      ( graph.dim.xdim - 1 ) * ( graph.dim.ydim )
      + ( graph.dim.xdim ) * ( graph.dim.ydim - 1 );

  int nbEdges = nbEdgePlane * graph.dim.tdim
                + ( graph.dim.xdim * graph.dim.ydim * ( graph.dim.tdim - 1 ) );

  int nbNodes = graph.v_Nodes.size();

  stream << "t " << graph.graphID << " " << graph.className << endl;
  stream << "dim " << graph.dim.xdim << " " << graph.dim.ydim << " " << graph.dim.tdim << endl;
  stream << "n " << nbEdges << endl;

  // === Output nodes
  // Write node on stream
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
    // For each node
    vector<GGraphEdge> v_Edges = graph.v_Nodes[i].v_Edges;
    for ( int j = 0, m = v_Edges.size(); j < m; ++j )
    {
      if ( v_Edges[j].destNodeId != GNONODEID )
        // Write edge on the output
        stream << "e " << i
               << " " << v_Edges.at(j).destNodeId
               << " " << v_Edges.at(j).label  <<  endl;
    }
  }
  return stream;
}
// End of operator << overload
