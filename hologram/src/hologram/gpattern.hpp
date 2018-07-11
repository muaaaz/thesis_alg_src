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
 *   Copyright (C) 2018 by Muaz Twaty                                      *
 *   muaz.sy123[at]gmail.com                                               *
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
#ifndef GRIMAPATTERN_HPP
#define GRIMAPATTERN_HPP

//================================ INCLUDE ===================================//
// Include library class
// Include project class
#include "ggraph.hpp"
#include "gdatabase.hpp"

//#include "../holeG/hgrimagraph.hpp"
//#include "../holeG/hgrimagraphecode.hpp"

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//
//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
/**
 * @brief The GPattern class
 * Class that inherite from GGraphCode the list of token (i.e. DFS code)
 * and allow to manage patterns.
 */
class GPattern
{
  //---- PUBLIC --------------------------------------------------------------//
public:
  // Public CONSTANTS ________________________________________________________//
  // Public Structure & Typedef ______________________________________________//
  // Public Variables ________________________________________________________//
  vector<GToken> v_Tokens;
  /// The graph code, maintained by push_back and pop_back
  GGraph * pGraph;
  /** List of occurrences of current pattern
    * Such that :
    * First  : Graph tid
    * Second : List of GNodeID that wear the pattern */
  vector<GNodeID> *v_OccList;
  /// Min time coord
  int minTCoord;
  /// Max time coord
  int maxTCoord;
  /// Node that handle minTCoord
  GNodeID nodeMinTCoord;
  /// Node that handle maxTCoord
  GNodeID nodeMaxTCoord;

  // Public Constructor/Desctructor __________________________________________//
  /// Default constructor
  GPattern();

  GPattern(GPattern* pat);
  /// Destructor
  ~GPattern();

  // Accessor ________________________________________________________________//
  // Mutator _________________________________________________________________//
  // Public Methods __________________________________________________________//
  void read( GReader *data, vector<string> &v_Token);

  /**
   * @brief push_back
   * Add an extension to the pattern, update graph and code
   * If extension allow to complete pattern, then add backward edge, s.t :
   *
   *    B          A--B             A--B
   *    | + A--B =    |   => Hence  |  |  where A->A is the first edges
   * A->A          A->A             A->A
   *
   * @param v_Ext : Code extension to add
   */
  void push_back(const GToken &v_Ext , bool canonTest );


  /**
   * @brief pop_back
   * Remove last extension(s) from the pattern. Undo what was done in last
   * push_back. If backward edge was added, then it's removed.
   */
  void pop_back( bool canonTest );

  /**
   * @brief isCanonincal
   * Test if pattern code is canonical
   * @return TRUE if Pattern is canonical, FALSE if it's not.
   */
  bool isCanonincal();

  GPattern* getCanonincalPattern();

  string token_to_string(GToken token);

  string getCanonincalString();
  
  void printOcc();


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
  /**
   * @brief setup
   * Method call by push_back for the first token. Clean node in pGraph, add the
   * two first nodes and edge between this of the pattern based on token.
   * @param token : DFS Code of the first edge.
   */
  void setup( const GToken &token );

  /**
   * TODO : RD
   */
  int findNotTemporalToken( GNodeID patNodeFrom );
};

//============================== OPERATOR OVERLOAD  ==========================//
struct GpatComptLt
{
  bool operator()( const GPattern *pat1, const GPattern *pat2 )
  {
    uint minSize = min( pat1->v_Tokens.size(), pat2->v_Tokens.size() );
    uint i = 0;
    int cmp = 0;

    while ( i < minSize )
    {
      cmp = cmpGTokenCanonTest( pat1->v_Tokens.at(i), pat2->v_Tokens.at(i) );
      if ( cmp > 0 )
        return true;
      else if ( cmp < 0 )
        return false;
      else
        i++;
    }
    if ( pat1->v_Tokens.size() > pat2->v_Tokens.size() )
      return true;
    else
      return false;
  }
};

ostream& operator<<(ostream& stream, GToken token );

ostream& operator<<(ostream& stream, vector<GToken> v_Tokens );

ostream& operator<<(ostream& stream, GPattern *pPattern );

//================================= END IFDEF ================================//
#endif // GPATTERN_HPP
