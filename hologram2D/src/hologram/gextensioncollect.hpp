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
#ifndef GRIMAEXTENSIONCOLLECT_HPP
#define GRIMAEXTENSIONCOLLECT_HPP

//================================ INCLUDE ===================================//
// Include library class
#include <map>
#include <list>
#include <cmath>
// Include project class
#include "gpattern.hpp"
#include "gsubgraphiso.hpp"
//#include "../holeG/hgrimagraphecode.hpp"

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//
/**
 * @brief The GExtensionData struct
 * Structure that store information about possible extension
 */
struct GExtensionData {
  /// Frequency of extension
  GGlobFreq frequency;
  int current_class_frequency;
  /// Vector that store graph were occurency
  vector<GTid>    v_Graphs;
  vector<GNodeID> v_OccList;
  /// Number of occurences of the extension
  GGlobFreq nbOcc;
  /// Graph ID to know when we change graph for frequency computing
  GTid  tId; // used to count the frequency correctly
};

//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//
//================================= CLASS ====================================//
/**
 * @brief The GExtensionCollect class
 * Class that compute all possible extension from a subgraph, check if these
 * extension are valids
 */
class GExtensionCollect
{
  //---- PUBLIC --------------------------------------------------------------//
public:
  // Public CONSTANTS ________________________________________________________//
  // Public Structure & Typedef ______________________________________________//
  // Public Variables ________________________________________________________//
  /// Map that store all possible extension and their code
  map<GToken,GExtensionData,GTokenGt> m_Extensions;
  /// Minimum global threshold
  GGlobFreq minFreqG;
  clock_t mapExtTick;
  int current_class_id;
  // Public Constructor/Desctructor __________________________________________//
  /// Default constructor
  GExtensionCollect( GGlobFreq minFG , int ID);

  /// Desctructor
  ~GExtensionCollect();

  // Accessor ________________________________________________________________//
  // Mutator _________________________________________________________________//
  // Public Methods __________________________________________________________//
  /**
   * TODO : RD
   * Write desc
   */
  void process( GSubgraphIso &subGraphIso );

  //---- PROTECTED  ----------------------------------------------------------//
protected:
  // Protected CONSTANTS _____________________________________________________//
  // Protected Variables _____________________________________________________//
  // Protected Methods _______________________________________________________//

  //---- PRIVATE -------------------------------------------------------------//
private:
  // Private CONSTANTS _______________________________________________________//
  // Private Structure _______________________________________________________//
  // Private Variables _______________________________________________________//
  // Private Methods _________________________________________________________//
  void getListExtension(
      vector<GGraphEdge> &v_Edges,
      vector<GGraphNode> &v_LargeNodes,
      pair<GToken,GExtensionData> &p_NewExt,
      vector<GGraphNode> &v_PatternNodes,
      GSubgraphIso &subGraphIso,
      GNodeID nodeLargeI,
      GNodeID nodeLargeJ,
      uint prevEdgeFrom,
      uint iNode,
      GPattern *pPattern,
      list< pair<GToken, GExtensionData> > &lp_ExtList );

  void getListExtensionTempFirst(vector<GGraphEdge> v_Edges,
                                 vector<GGraphNode> &v_LargeNodes ,
                                 pair<GToken, GExtensionData> p_NewExt,
                                 vector<GGraphNode> &v_PatternNodes,
                                 GSubgraphIso &subGraphIso,
                                 GNodeID nodeLargeI,
                                 GNodeID nodeLargeJ,
                                 uint prevEdgeFrom,
                                 uint iNode,
                                 GPattern *pPattern,
                                 list< pair<GToken, GExtensionData> > &lp_ExtList
                                 );

  void getListExtensionSpatFirst( vector<GGraphEdge> v_Edges,
                                  vector<GGraphNode> &v_LargeNodes ,
                                  pair<GToken, GExtensionData> p_NewExt,
                                  vector<GGraphNode> &v_PatternNodes,
                                  GSubgraphIso &subGraphIso,
                                  GNodeID nodeLargeI,
                                  GNodeID nodeLargeJ,
                                  uint prevEdgeFrom,
                                  uint iNode,
                                  GPattern *pPattern,
                                  list< pair<GToken, GExtensionData> > &lp_ExtList
                                  );


};

//============================== OPERATOR OVERLOAD  ==========================//
//================================= END IFDEF ================================//
#endif // GEXTENSIONCOLLECT_HPP
