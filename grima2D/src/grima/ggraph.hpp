/***************************************************************************
 *   Copyright (C) 2005 by Siegfried Nijssen                               *
 *   snijssen[at]informatik.uni-freiburg.de                                *
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

//================================= IFDEF ====================================//
#ifndef GGRAPH_HPP
#define GGRAPH_HPP

//================================ INCLUDE ===================================//
// Include library class
// Include project class
#include "gglobal.hpp"
#include "greader.hpp"

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//
/**
 * @brief GGraphEdge
 * Describe an edge with its label, its ID and the destination node
 */
struct GGraphEdge {
  /// Edge label
  GEdgeLabel label;
  /// Node destination ID
  GNodeID   destNodeId;
  /**
   * The equivalent edge is :
   * - ( i + 2 ) % 4 in nodes[from].edges[i] for edge in (x,y) plane
   * - 6 if edge is forward in time and 5 if edge is backward in timer
   */
  GEdgeID   eqEdgeId;

  /// Default Constructor
  GGraphEdge():
    label     ( GNOEDGEID ),
    destNodeId( GNONODEID ),
    eqEdgeId  ( GNOEDGEID )
  {
    // Nothing to to here
  }

  /// Operator equal to compare two edges
  bool operator==( const GGraphEdge& m ) const
  {
    return ( ( m.label==label ) && ( m.destNodeId == destNodeId ) );
  }
};
// End of struct GGraphEdge

/**
 * @brief GGraphNode
 * Describe a graph node with its label and a list of edge that go from this
 * node and its coordinate.
 */
struct GGraphNode {
  /// Node label
  GNodeLabel label;
  /// Node coord along the x-axis
  int xCoord;
  /// Node coord along the y-axis
  int yCoord;
  /// Node coord along the t-axis
  int tCoord;
  /**
   * Vector that list edges with this node as origin in clockwise order for the
   * fourth first edge. Forward in time for the fifth edge and backward in time
   * for sixth edge.
   */
  vector<GGraphEdge>  v_Edges;

  /// Default constructor
  GGraphNode():
    label( GNONODEID ),
    xCoord( 0 ),
    yCoord( 0 ),
    tCoord( 0 )
  {
    // Create null edge object
    GGraphEdge nullEdge;
    // Resize list of edges for this node
    v_Edges.resize( 6, nullEdge );
  }

  /// == Operator overload to compare node when searching if node exist
  bool operator==( const GGraphNode& node )
  {
    return ( label == node.label
             && xCoord == node.xCoord
             && yCoord == node.yCoord
             && tCoord == node.tCoord );
  }
};
// End of struct GGraphNode

/**
 * @brief GGraphDim
 * Store the graph dimension along all axis
 */
struct GGraphDim {
  /// Dimension along x-axis, i.e. nb node in x
  int xdim;
  /// Dimension along x-axis, i.e. nb node in y
  int ydim;
  /// Dimension along x-axis, i.e. nb node in t
  int tdim;

  /// Default constructor
  GGraphDim():
    xdim( 0 ),
    ydim( 0 ),
    tdim( 0 )
  {
    // Nothing to put here
  }
};
// End of struct GGraphDim

//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//
//================================= CLASS ====================================//
/**
 * @brief GGraph
 * Manage graph object with a list of nodes and edges and some function to
 * manipulate graph and get informations
 */
class GGraph
{
  //---- PUBLIC --------------------------------------------------------------//
public:
  // Public CONSTANTS ________________________________________________________//
  // Public Structure & Typedef ______________________________________________//
  // Public Variables ________________________________________________________//
  /// List of graph nodes that compose the graph
  vector<GGraphNode> v_Nodes;
  //  /// Vector that store the number of occurence of a pattern
  vector<int> v_NbOccPattern;
  /// Class of the current graph
  string className;
  int classID;
  /// Id of the current graph
  uint graphID;
  /// Graph Dimension
  GGraphDim dim;

  // Public Constructor/Desctructor __________________________________________//
  /// Default constructor
  GGraph();
  /// Copy constructor
  GGraph( GGraph* src );
  /// Default Destructor
  ~GGraph();

  // Accessor ________________________________________________________________//
  // Mutator _________________________________________________________________//
  // Public Methods __________________________________________________________//
  /**
   * TODO: RD
   * Write Desc
   */
  void read( GReader        *data,
             vector<string> &v_GraphIdToken ,
             uint           &suppGraphID,
             string         suppClassName );

  /**
   * @brief addNode
   * Check if no other node share same coordinate, if so, print an error and
   * exit, else add node at the end of the list of node v_Nodes
   * @param nodeLabel: label of node added
   */
  void addNode( GNodeLabel nodeLabel, int xcoord, int ycoord, int tcoord );

  /**
   * @brief deleteNode
   * Remove the last node in v_Nodes
   */
  void deleteNode();

  /**
   * @brief setNodeToNoLabel
   * Set node with ID 'pos' label to GNOLABEL and remove edge that link to
   * that node.
   * NB : It aim to 'remove' this node from the graph to not be considered
   * during mining.
   * @param pos : position of node to 'remove' from v_Nodes
   */
  void setNodeToNoLabel( uint pos );
  void setNodeToNoLabel( vector<GNodeID> v_NodeId );

  /**
   * TODO
   */
  int computeEdge();

  /**
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
  void addEdge( GNodeID    nodeIdFrom,
                GNodeID    nodeIdDest,
                GEdgeLabel edgeLabel ,
                int        line = -1 );

  /**
   * TODO
   */
  void addEdgeTor( GNodeID    nodeIdFrom,
                   GNodeID    nodeIdDest,
                   GEdgeLabel edgeLabel );

  /**
   * @brief deleteEdge
   * Find if edge between nodeIdFrom and nodeIdDest exist. If not, print an
   * error and exit, else, compute position of this edge in the list v_Edges
   * of both nodeIdFrom and nodeDest and set it de GNONODEID.
   * @param nodeIdFrom : Node from which edge to delete start
   * @param nodeIdDest : Node from which edge to delete end
   */
  void deleteEdge( GNodeID nodeIdFrom, GNodeID nodeIdDest );

  /**
   * @brief edgeExist
   * Search if there is an edge between nodeIdFrom and nodeIdDest, return true
   * if such edge found, else return false.
   * @param nodeIdFrom : Node from which edge to find start
   * @param nodeIdDest : Node from which edge to find end
   * @return FALSE if edge is not equal to GNOEDGEID, else return TRUE
   */
  bool edgeExist( GNodeID nodeIdFrom, GNodeID nodeIdDest );
  /**
   * @brief edgeExist
   * Search if there is an edge between nodeIdFrom and nodeIdDest, return true
   * if such edge found, else return false.
   * @param nodeIdFrom : Node from which edge to find start
   * @param nodeIdDest : Node from which edge to find end
   * @param pos        : Variable to know position of edge if found
   * @return FALSE if edge is not equal to GNOEDGEID, else return TRUE
   */
  bool edgeExist( GNodeID nodeIdFrom, GNodeID nodeIdDest, int &pos );

  /**
   * TODO
   */
  vector<GNodeID> processHole();

  vector<GNodeID> processFrameHole();

  /**
   * TODO
   */
  uint getEdgePositionFrom( GNodeID nodeFrom, GNodeID nodeDest );

  //---- PROTECTED  ----------------------------------------------------------//
protected:
  // Protected CONSTANTS _____________________________________________________//
  // Protected Variables _____________________________________________________//
  // Protected Methods _______________________________________________________//

  //---- PRIVATE -------------------------------------------------------------//
private:
  // Private CONSTANTS _______________________________________________________//
  // Private Structure _______________________________________________________//
  /**
   * @brief The error struct
   * Structure with multiple parameterize method to pring error message
   */
  struct error {
    /**
     * @brief structError
     * Print the message when a line does not start with the right char
     * @param err  : Character wanted
     * @param line : Line wrong in the graph file
     */
    static void structError( string err, int line )
    {
      cerr << "ERROR - There's no line beginning with a '"
           << err << "' at line" << line  << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief structError
     * Print the message when a line does not start with the right char after
     * a specific line
     * @param err  : Character wanted
     * @param prev : Character of the previous line
     * @param line : Line wrong in the graph file
     */
    static void structError( string err, string prev, int line )
    {
      cerr << "ERROR - There's no line beginning with a '"
           << err << "' after '" << prev << "' line at line" << line << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief structValue_UINT_Error
     * Print the message when a variable is not an unsigned int in the graph file
     * @param err  : Variable that raise this message
     * @param line : Line wrong in the graph file
     */
    static void structValue_UINT_Error( string err, int line )
    {
      cerr << "ERROR - " << err << " MUST BE UNSIGNED INT at line " << line << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief structValue_UINT_Error
     * Print the message when a variable is not an unsigned int > 0 in the graph file
     * @param err  : Variable that raise this message
     * @param line : Line wrong in the graph file
     */
    static void structValue_PosUINT_Error( string err, int line )
    {
      cerr << "ERROR - " << err << " MUST BE UNSIGNED INT > 0 at line " << line << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief graphIDError
     * Print a message when two graph of the same class share the same ID
     * @param graphID : Graphe ID shared between two graph
     * @param line    : Line of the file where error occurs
     */
    static void graphIDError( uint graphID, uint suppGraphID, int line )
    {
      cerr << "ERROR - Two or more graph of the same classe share the same ID : "
           << graphID << ". The second one occurs at line : " << line
           << ". Its GraphID should be : " << suppGraphID << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief nodeMissingValError
     * Print the message when node line does not have enough variable set
     * @param nodeID : ID of the node that raise this message
     * @param line : Line wrong in the graph file
     */
    static void nodeMissingValError( string  nodeID, int line )
    {
      cerr << "ERROR - Some informations are missing missing for the node "
           << nodeID << " at line " << line << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief nodeValueError
     * Print the message when node variable are not UINT
     * @param nodeID : ID of the node that raise this message
     * @param line : Line wrong in the graph file
     */
    static void nodeValue_UINT_Error( string nodeID, int line )
    {
      cerr << "ERROR - Values of node " << nodeID
           << " MUST BE UNSIGNED INT at line " << line << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief nodeCoordError
     * Print the message when coordinate of a node are not those that are
     * supposed to be
     * @param nodeID : ID of the node
     * @param x      : x coordinate supposed
     * @param y      : y coordinate supposed
     * @param t      : t coordinate supposed
     * @param line : Line wrong in the graph file
     */
    static void nodeCoordError( string nodeID, int x, int y, int t, int line )
    {
      cerr << "ERROR - Coordinate of node " << nodeID
           << " SHOULD BE (" << x << "," << y << "," << t << ") at line " << line << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief edgeMissingValError
     * * Print the message when edge line does not have enough variable set
     * @param line : Line wrong in the graph file
     */
    static void edgeMissingValError( int line )
    {
      cerr << "ERROR - Some information are missing for edge at line " << line << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief edgeValueError
     * Print the message when node variable are not UINT
     * @param line : Line wrong in the graph file
     */
    static void edgeValue_UINT_Error( int line )
    {
      cerr << "ERROR - Values of edge MUST BE UNSIGNED INT at line " << line << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief nbEdgeError
     * Print the message when number of edge declared in the file are not
     * the same as computed while reading line that start with 'e'
     * @param real     : Real value computed
     * @param declared : Declared value in the file
     * @param line : Line wrong in the graph file
     */
    static void nbEdgeError( uint real, uint declared, int line )
    {
      cerr << "ERROR - The number of edges specified at line " << line
           << " which is " << declared << " is not the same as computed when "
           << "adding edge (i.e. while reading line starting with 'e'." << endl
           << "Number of edge computed is : " << real << endl
           << "Please review :" << endl
           << "- Either your edges declarations (line starting with 'e')" << endl
           << "- Either the number of delcared edge at line : " << line << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief addNodeError
     * Print the message when trying to add node with coordinate that are the
     * same that another node
     * @param graph     : Graph on which there is a probleme
     * @param errorNode : Node wich have wrong coordinate
     * @param nodeID    : ID of the first node that have those coordinate
     */
    static void addNodeError( GGraph     *graph,
                              GGraphNode errorNode,
                              uint       nodeID )
    {
      cerr << "ERROR - Trying to add node with coordinate that are already "
           << "given to another node. Here are information about those node : " << endl
           << "Graph className : " << graph->className << endl
           << "Graph ID : " << graph->graphID << endl
           << "Existing node ID : " << nodeID << " with coord ("
           << errorNode.xCoord << ","
           << errorNode.yCoord << ","
           << errorNode.tCoord << ")" << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief addEdgeNotNeighborError
     * Print the message when trying to add an edge between to node that are
     * not neighbor
     * @param nodeIdFrom : Id of the node from which the edge start
     * @param nodeIdDest : Id of the node to which the edge end
     * @param line       : Line wrong in the graph file (OPTIONNAL)
     */
    static void addEdgeNotNeighborError( GNodeID  nodeIdFrom,
                                         GNodeID  nodeIdDest,
                                         int      line = -1 )
    {
      cerr << "ERROR - Trying to add edge between nodes " << nodeIdFrom
           << " and node " << nodeIdDest << " is invalid because those node are "
           << "not neighbor !" << endl;
      if ( line != -1 )
        cerr << "This error occur while reading graph file at line " << line << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief removeEdgeNotExistError
     * Print the message when trying to remove an edge wich is already set
     * to GNOEDGEID.
     * @param nodeIdFrom : Id of the node from which the edge start
     * @param nodeIdDest : Id of the node to which the edge end
     */
    static void removeEdgeNotExistError( GNodeID nodeIdFrom,
                                         GNodeID nodeIdDest )
    {
      cerr << "ERROR - Trying to remove non-existant edge between node "
           << nodeIdFrom << " and " << nodeIdDest << endl;
      exit( EXIT_FAILURE );
    }
  };
  // End of struct error

  // Private Variables _______________________________________________________//
  // Private Methods _________________________________________________________//
  /**
   * TODO : RD
   * Write Desc
   */
  void readGraphID( vector<string> &v_Tokens,
                    uint           &suppGraphID,
                    string         suppClassName,
                    int            line );

  /**
   * @brief readDimension
   * Check if line start with right char and if value of dimensiont are unsigned
   * int, if not write error and exit, else set them
   * @param v_Tokens : String token to parse
   * @param line     : Line actually read by GReader
   */
  void readDimension( vector<string> &v_Tokens, int line );

  /**
   * @brief readNbEdges
   * Check if line start with right car and if value of number of edge is
   * unsigned int, if not write error and exit, else return the number of edge
   * set in the file
   * @param v_Tokens : String token to parse
   * @param line     : Line actually read by GReader
   * @return The number of edge
   */
  GEdgeID readNbEdges( vector<string> &v_Tokens, int line );

  /**
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
  void readVertex( vector<string> &v_Tokens,
                   int &suppXCoord,
                   int &suppYCoord,
                   int &suppTCoord ,
                   int line );

  /**
   * @brief readEdge
   * Check if line that start with 'e', print error if line does not respect
   * specified format
   * @param v_Tokens      : Tokens of the string to be parsed
   * @param currentEdgeId : Currently computed/read edgeID
   * @param line          : Line actually read by GReader
   */
  void readEdge( vector<string> &v_Tokens,
                 uint           &currentEdgeId,
                 int            line );

  /**
   * TODO
   */
  vector<GNodeID> computeNeighbor( GNodeID iNode );
  /**
   * TODO
   */
  vector<GNodeID> computeDirectNeighbor( GNodeID iNode );
  /**
   * TODO
   */
  vector<GNodeID> computeIndirectNeighbor( GNodeID iNode );

  /**
   * TODO
   */
  vector<GNodeID> computeNeighborHole( GNodeID iNode );

  /**
   * TODO
   */
  vector<GNodeID> getIndirectNeighbor2D(GNodeID iNode );

  /**
   * TODO
   */
  vector<GNodeID> getIndirectNeighbor2Dt(GNodeID iNode);

  /**
   * TODO
   */
  vector<GNodeID> getDirectNeighbor2D(GNodeID iNode );

  /**
   * TODO
   */
  vector<GNodeID> getDirectNeighbor2Dt(GNodeID iNode);

  /**
   * TODO
   */
  vector<GNodeID> getIndirectNeighbor2D_C ( GNodeID iNode ); // Bottom Left Node in 2D grids
  vector<GNodeID> getIndirectNeighbor2D_BL( GNodeID iNode ); // Bottom Left Node in 2D grids
  vector<GNodeID> getIndirectNeighbor2D_BR( GNodeID iNode ); // Bottom Right Node in 2D grids
  vector<GNodeID> getIndirectNeighbor2D_TL( GNodeID iNode ); // Top Left Node in 2D grids
  vector<GNodeID> getIndirectNeighbor2D_TR( GNodeID iNode ); // Top Left Node in 2D grids
  vector<GNodeID> getIndirectNeighbor2D_B ( GNodeID iNode ); // Bottom Node in 2D grids
  vector<GNodeID> getIndirectNeighbor2D_R ( GNodeID iNode ); // Left Node in 2D grids
  vector<GNodeID> getIndirectNeighbor2D_T ( GNodeID iNode ); // Top Node in 2D grids
  vector<GNodeID> getIndirectNeighbor2D_L ( GNodeID iNode ); // Righ Node in 2D grids

  vector<GNodeID> getDirectNeighbor2D_C ( GNodeID iNode ); // Bottom Left Node in 2D grids
  vector<GNodeID> getDirectNeighbor2D_BL( GNodeID iNode ); // Bottom Left Node in 2D grids
  vector<GNodeID> getDirectNeighbor2D_BR( GNodeID iNode ); // Bottom Right Node in 2D grids
  vector<GNodeID> getDirectNeighbor2D_TL( GNodeID iNode ); // Top Left Node in 2D grids
  vector<GNodeID> getDirectNeighbor2D_TR( GNodeID iNode ); // Top Left Node in 2D grids
  vector<GNodeID> getDirectNeighbor2D_B ( GNodeID iNode ); // Bottom Node in 2D grids
  vector<GNodeID> getDirectNeighbor2D_R ( GNodeID iNode ); // Left Node in 2D grids
  vector<GNodeID> getDirectNeighbor2D_T ( GNodeID iNode ); // Top Node in 2D grids
  vector<GNodeID> getDirectNeighbor2D_L ( GNodeID iNode ); // Righ Node in 2D grids

  /**
   * TODO
   */
  bool areNeighbor( GGraphNode nodeFrom, GGraphNode nodeDest );

  /**
   * TODO
   */
  bool areNeighborHole( GGraphNode nodeFrom, GGraphNode nodeDest );

  /**
   * @brief isInt
   * Test if string is an unsigned int
   * @param s : Input string
   * @return True if s is an uint false elsewhere
   */
  bool isUInt( const string &s );
};

//============================== OPERATOR OVERLOAD  ==========================//
/**
 * @brief operator <<
 * Overload operator << to pring graph (ie list of node and edges.
 * @param stream : output stream
 * @param graph : graph to print
 * @return
 */
ostream& operator<<(ostream& stream, GGraph &graph );

//============================== OPERATOR OVERLOAD  ==========================//
//================================= END IFDEF ================================//
#endif // GGRAPH_HPP
