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
#ifndef GRIMACANONICALDFSCOMPUTER_HPP
#define GRIMACANONICALDFSCOMPUTER_HPP

//================================ INCLUDE ===================================//
// Include library class
#include <vector>
#include <map>
// Include project class
#include "gglobal.hpp"
#include "ggraph.hpp"
#include "gpattern.hpp"

//=============================== NAMESPACE ==================================//
using namespace std;

//============================ STRUCT & TYPEDEF ==============================//
/**
 * @brief The GCanonicalToken struct
 * Structure that store possible canonical extension and their correspondance
 * nodes in the large graph
 */
struct GCanonicalToken {
  /// DFSToken code
  GToken codeToken;
  /// Equiv node from in large graph
  GNodeID nodeLargeI;
  /// Equiv node dest in large graph
  GNodeID nodeLargeJ;
  /// Default Constructor
  GCanonicalToken():
    nodeLargeI( GNONODEID ),
    nodeLargeJ( GNONODEID )
  {
    // Nothing do to here
  }
};

/**
 * @brief The GEdge struct
 * Structure that describe an edge by the label of it node from and node dest.
 * Use for canonical comparison
 */
struct GEdge {
  /// Node from label
  GNodeLabel nodeI;
  /// Node dest label
  GNodeLabel nodeJ;
};

/**
 * @brief The GCanonicalDFSTokenLt struct
 * Structure that overload operator () to compare two if
 * GCanonicalToken1 < GCanonicalToken2
 */
struct HGCanonicalDFSTokenLt {
  /**
   * @brief operator ()
   * Operator () overload that compare t1.codeToken & t2.codeToken
   * @param t1 : DFS "canonical" Code 1
   * @param t2 : DFS "canonical" Code 2
   * @return TRUE  if t1.codeToken < t2.codeToken
   *         FALSE if t1codeToken >= t2.codeToken
   */
  bool operator()( const GCanonicalToken &t1,
                   const GCanonicalToken &t2 ) const
  {
    return cmpGTokenCanonTest( t1.codeToken, t2.codeToken ) < 0;
    // The order of the nodeId does not matter.
  }
};

/**
 * @brief The GCanonicalDFSTokenGt struct
 * Structure that overload operator () to compare two if
 * GCanonicalToken1 < GCanonicalToken2
 */
//struct HGCanonicalDFSTokenGt {
//  /**
//   * @brief operator ()
//   * Operator () overload that compare t1.codeToken & t2.codeToken
//   * @param t1 : DFS "canonical" Code 1
//   * @param t2 : DFS "canonical" Code 2
//   * @return TRUE  if t1.codeToken > t2.codeToken
//   *         FALSE if t1codeToken <= t2.codeToken
//   */
//  bool operator()( const GCanonicalToken &t1, const GCanonicalToken &t2 ) const
//  {
//    return cmpGTokenCa ( t1.codeToken, t2.codeToken ) > 0;
//    // The order of the nodeId does not matter.
//  }
//};

//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
/**
 * @brief The GCanonicalDFSComputer class
 * Class that compute possible code for a pattern to know if it's canonical.
 */
class GCanonicalDFSComputer
{
  //---- PUBLIC --------------------------------------------------------------//
public:
  // Public CONSTANTS ________________________________________________________//
  // Public Structure & Typedef ______________________________________________//
  // Public Variables ________________________________________________________//
  // Public Constructor/Desctructor __________________________________________//
  /// Default constructor
  GCanonicalDFSComputer(GPattern *pPat, bool _compute_code );
  /// Desctructor
  ~GCanonicalDFSComputer();

  // Accessor ________________________________________________________________//
  // Mutator _________________________________________________________________//
  // Public Methods __________________________________________________________//
  /**
   * @brief isCanonincal
   * Function that check if code of Pattern is canonical in three step.
   * Step One   : If the pattern is only one edge, then it's canonical
   * Step Two   : Check in the full DFS code of the pattern if there's no clue
   *              to directly know that pattern is not canonical.
   * Step Three : If step one and step two passed, then begin a recursion to
   *              construct all possible DFS Code of the pattern to find if
   *              there's a DFS Code greater than the current one.
   * @return TRUE if pattern is canonical, FALSE else.
   */
  bool isCanonincal();
  GPattern* getCanonincal();

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
  /// A pointer to the pattern to test
  GPattern *pPattern;
  /// A version of the new pattern
  GPattern *pTestNewPat;
  
  GPattern *largestNewPat;
  /// A vector of DFSToken with all possible extenstion for the pattern
  vector<GCanonicalToken> v_StackCode;
  /// A vector of DFSToken buffer, used to keep equivalent backward edge added
  vector<GCanonicalToken> v_Buffer;
  /// A vector that is the possible new code that can represent the pattern
  vector<GCanonicalToken> v_NewCode;
  /**
   * Vectors of Node for the mapping between :
   * - The new pattern nodes (v_CodeIndex)
   * - The old pattern nodes (v_GraphIndex)
   */
  vector<GNodeID> v_CodeIndex;
  vector<GNodeID> v_GraphIndex;
  /// Counter of nodes in the new pattern
  int nbNodes;
  /// Counter of the position in the pattern
  int position;
  /// Boolean to know if we stop the test because the pattern is not canonical
  bool stop;

  bool largerPatternFound;

  bool cocoboolean ;

  int debug_c;

  bool compute_code;
  // Private Methods _________________________________________________________//
  /**
   * @brief recurse
   * Methods that construct all possible extension of the current new pattern
   * and store them in v_StackCode.
   * If prevAngle==GNOANGLE,
   *  construct all extenstion from nodeLargeI and nodeLargeJ
   * Else
   *  construct only all exteion from nodeLargeJ
   * Once done, it will sort these extenstion and call subRecurse() to test it.
   * @param prevAngle  : The orientation of the previous edge
   * @param nodeLargeI  : The node from of the last construct edge
   * @param nodeLargeJ  : The node dest of the last construct edge
   */
  void recurse( GEdgeAngle prevAngle, GNodeID nodeLargeI,
                GNodeID nodeLargeJ );

  void get_recurse( GEdgeAngle prevAngle, GNodeID nodeLargeI,
                    GNodeID nodeLargeJ );
  /**
   * @brief subRecurse
   * Methods that will add the greatest extension previously construct from
   * recurse() to the new DFS Code and compare it to the equivalent line code of
   * the pattern to test.
   * If the pattern DFS code is less than the new one
   *  the pattern DFS Code is not canonical, exit the test
   * Else if the pattern DFS code is greater than the new one
   *  Exit subRecurse() to test another extension
   * Call recurse to extends the current new pattern from the last node or all
   * the previous one
   */
  void subRecurse();


  void get_subRecurse();


  GEdgeAngle getTokenAngleValue(uint    prevEdgeFrom,
                                 GNodeID nodeFrom,
                                 GNodeID nodeDest);

  int findNotTemporalToken(GNodeID patNodeFrom , GPattern *pPattern);

};

//============================== OPERATOR OVERLOAD  ==========================//
//================================= END IFDEF ================================//
#endif // GCANONICALDFSCOMPUTER_HPP
