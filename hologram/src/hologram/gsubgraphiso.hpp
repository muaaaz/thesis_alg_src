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

//================================= IFDEF ====================================//
#ifndef GRIMASUBGRAPHISO_HPP
#define GRIMASUBGRAPHISO_HPP

//================================ INCLUDE ===================================//
// Include library class
// Include project class
#include "gpattern.hpp"

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//
//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//


//================================= CLASS ====================================//
/**
 * @brief The GSubgraphIso class
 * Class that will allow to do isomorphisms test of a pattern in a largegraph
 * from first canonical edges.
 */
class GSubgraphIso
{
  //---- PUBLIC --------------------------------------------------------------//
public:
  // Public CONSTANTS ________________________________________________________//
  // Public Structure & Typedef ______________________________________________//
  // Public Variables ________________________________________________________//
  /// The pattern to find
  GPattern *pPattern;
  /// A large graph from the database
  GGraph *largeGraph;
  /// Number of node of pattern
  uint patternNbNode;


  int current_occurrence_class_id;
  /**
   * Integer to store offset in list v_OccList to know from where to put
   * nodes id in the large graph that are isomorph to pattern
   */
  uint offset;
  /**
   * Boolean that informe if the edge is the origin of the pattern.
   * If FALSE : Edge occurence does not "wear" pattern
   * If TRUE  : Edge occurence does "wear" pattern (i.e. we can construct
   *      pattern from this edge
   */
  bool lastOccValid;
  /**
   * Vectors that make a mapping between nodes of pattern and nodes in the
   * large graph from database.
   * If i is a node of the database graph, pattern_index[i] is the
   * corresponding node in the pattern (or GNONODEID there is no
   * correspondance). If there is a correspondance, then:
   * -  (*largegraph_index)[offset +pattern_index[i]] = i
   * -  pattern_index[(*largegraph_index)[offset + i]] = i */
  vector<GNodeID> v_PatternIdx;
  vector<GNodeID> *v_LargeGraphIdx;
  /// Vector that will let to know which was the last edge used in the pattern
  vector<GEdgeID> v_LastEdgeUsed;

  // Public Constructor/Desctructor __________________________________________//
  /**
   * @brief GSubgraphIso
   * Default constructor
   * Define the pattern that must be found and the vector that store occurence
   * @param pPat, the pattern that must be found
   * @param v_OccList, A vector of GNodeId that store the occurence
   */
  GSubgraphIso( GPattern *pPat, vector<GNodeID> *v_OccList );
  /// Desctructor
  ~GSubgraphIso();

  // Accessor ________________________________________________________________//
  // Mutator _________________________________________________________________//
  // Public Methods __________________________________________________________//
  /**
   * @brief clearOccList
   * Method that fill v_LargeGraphIdx with GNONODEID
   */
  void clearOccList();

  /**
   * @brief run
   * Function that returns true if an occurence of the pattern is found in
   * pNewLargeGraph. largeId and edgeId indicate where the occurence must be
   * look for in largeGraph. The first edge of the pattern must already match
   * @param pNewLargeGraph :pointer to the large graph
   * @param nodeLargeIdFrom: node ID from of the first edge in the large graph
   * @param edgeLargeId: edge position in the list of edge of the node from
   * @return TRUE if edge wear occurence, FALSE else.
   */
  bool run( GGraph* pNewLargeGraph,  GNodeID nodeLargeIdFrom,
            GEdgeID edgeLargeId );

  /**
   * TODO : RD
   * Write Desc
   */
  GEdgeAngle getTokenAngleValue( uint    prevEdgeFrom,
                                 GNodeID nodeFrom,
                                 GNodeID nodeDest );


  //---- PROTECTED  ----------------------------------------------------------//
protected:
  // Protected CONSTANTS _____________________________________________________//
  // Protected Structure _____________________________________________________//
  // Protected Variables _____________________________________________________//
  // Protected Methods _______________________________________________________//

  //---- PRIVATE -------------------------------------------------------------//
private:
  // Private CONSTANTS _______________________________________________________//
  // Private Structure _______________________________________________________//
  // Private Variables _______________________________________________________//
  // Private Methods _________________________________________________________//
  /**
   * TODO : RD
   * Write desc
   */
  bool subRun();

  /**
   * TODO : RD
   * Write desc
   */
  int findNotTemporalToken( GNodeID patNodeFrom );

  /**
   * @brief getEdge
   * Function that return edge in function of its angle in pattern and the
   * previous edge angle in the large graph
   * @param token: token with angle, ie angle of this edge from the previous one
   * @param prevPosition : orientation of the previous edge in the large graph
   * @return return the edge that extent pattern in the large graph
   */
  GGraphEdge getEdge( GToken &token, int prevPosition );

};

//============================== OPERATOR OVERLOAD  ==========================//
//================================= END IFDEF ================================//
#endif // GSUBGRAPHISO_HPP
