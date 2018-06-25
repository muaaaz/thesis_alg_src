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
#ifndef GGLOBAL_HPP
#define GGLOBAL_HPP

//================================ INCLUDE ===================================//
// Include library class
#include <string>
#include <string.h>
// Include project class
#include "global.hpp"

//=============================== NAMESPACE ==================================//
using namespace std;

//================================= CLASS ====================================//
class GGraph;

//============================ STRUCT & TYPEDEF ==============================//
typedef unsigned int GNodeLabel;
typedef unsigned int GEdgeLabel;
typedef unsigned int GNodeID;
typedef unsigned int GEdgeID;
typedef unsigned int GGlobFreq;
typedef int GEdgeAngle;
typedef GGraph *GTid;
typedef string *GTclass;

#define GNONODEID ((GNodeID) -1)
#define GNOEDGEID ((GEdgeLabel) -1)
#define GNOLABEL  ((GNodeLabel) -1)
#define GNOANGLE  ((GEdgeAngle) -99)

/**
 * @brief The GTokenDirection enum
 * Enum that list the two possible direction of an edge.
 * - gForward : Edge that go to a new node
 * - gBackward : Edge that go to an already know node.
 */
enum GTokenDirection { gForward = 0 , gBackward = 1 };

/**
 * @brief The GToken struct
 * Structure of token as describe in the class description
 */
struct GToken {
  /** Information to know edge direction
   * -> backward means that this edge will close the pattern
   */
  GTokenDirection direction;
  /// Id of the node from (N_i)
  GNodeID    nodeFrom;
  /// Id of the node destination (N_j)
  GNodeID    nodeDest;
  /// Angle between (N_i-1, N_i) & (N_i, N_j)
  GEdgeAngle angle;
  /// Label of the node from
  GNodeLabel nodeLabelFrom;
  /// Label of the edge
  GEdgeLabel edgeLabel;
  /// Label of the node destination
  GNodeLabel nodeLabelDest;
  /// Default constructor
  GToken():
    direction( gForward ),
    nodeFrom ( GNONODEID ),
    nodeDest ( GNONODEID ),
    angle    ( GNOANGLE  ),
    nodeLabelFrom( GNOLABEL ),
    edgeLabel    ( GNOLABEL ),
    nodeLabelDest( GNOLABEL )
  {
    // Nothing to do here
  }

  bool operator == (const GToken ob){
    return direction == ob.direction &&
            nodeFrom == ob.nodeFrom &&
            nodeDest == ob.nodeDest &&
            angle == ob.angle &&
            nodeLabelFrom == ob.nodeLabelFrom &&
            edgeLabel == ob.edgeLabel &&
            nodeLabelDest == ob.nodeLabelDest;
  }
};

extern struct GToken tmpT1;
extern struct GToken tmpT2;


//=== END IFDEF ================================================================
#endif // GRIMAGLOBAL_HPP
