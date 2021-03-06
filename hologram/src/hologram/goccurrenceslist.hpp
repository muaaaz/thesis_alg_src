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
#ifndef HOCCURRENCESLIST_HPP
#define HOCCURRENCESLIST_HPP

//================================ INCLUDE ===================================//
// Include library class
#include <vector>
//#include <array>
//#include <map>
//#include <unistd.h>
// Include project class
#include "gglobal.hpp"
//#include "../holeG/hgrimagraphecode.hpp"

//=============================== NAMESPACE ==================================//
using namespace std;
//============================ STRUCT & TYPEDEF ==============================//
//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
/**
 * @brief The OCCURRENCES LIST class
 * Class that will managed occurence of pattern without having to duplicate it
 */
class GOccurrencesList
{
  //---- PUBLIC --------------------------------------------------------------//
public:
  // Public CONSTANTS ________________________________________________________//
  // Public Structure & Typedef ______________________________________________//
  /// Structure that define an edge store in map
  struct mapEdge {
    GNodeID nodeFrom;
    GNodeID nodeDest;
    GEdgeID edgeId;
    uint    element;
    bool operator == (const mapEdge ob) const {
      return nodeFrom == ob.nodeFrom &&
             nodeDest == ob.nodeDest &&
             edgeId   == ob.edgeId &&
             element  == ob.element ;
             
    }
  };

  /// Domain that store element
  //vector<uint>    v_Domain;
  /// Map that store value of the element
  vector<mapEdge> data;
  /// Size of unremoved domain
  /// Ggraph ID
  uint graphID;
  /// Graph Memory ID
  GGraph *pGraph;

  /// Default constructor
  GOccurrencesList();

  GOccurrencesList( uint graphId, GGraph *p_Graph );

  /// Default destructor
  ~GOccurrencesList();

  int size()
  {
    return data.size();
  }
  // Accessor ________________________________________________________________//
  /**
   * TODO : RD
   * Write desc
   */
  mapEdge at( uint i );


  // Public Methods _________________________________________________________//
  /**
   * TODO : RD
   * Write desc
   */
  void add( mapEdge edge );

  /**
   * TODO : RD
   * Write Desc
   */
  void add( GNodeID from, GNodeID dest, GEdgeID edge );

  /**
   * TODO : RD
   * Write Desc
   */
  void remove( mapEdge e );

  /**
   * TODO : RD
   * Write Desc
   */
  void remove( uint i );


  void clear();

  //  bool contains( mapEdge e )
  //  {
  //    /*
  //       * Function that return search if mapEdge object is in the domain,
  //       * if not return -1
  //       * else
  //       *   check if element in domain is before size, and return result.
  //       */

  //    int idx = find (e);
  //    if ( idx != -1 )
  //      return v_Map[idx].element < size;
  //    else
  //      return false;
  //  }
  //  // End of contains( mapEdge e )

  //  bool contains(uint element)
  //  {
  //    /*
  //       * Function that check if element is in active domain, ie if element
  //       * position is less than size.
  //       */

  //    return v_Map[v_Domain[element]].element < size ;

  //  }
  //  // End of contains(sparseElement element)

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
   * TODO : RD
   * Write desc
   */
  int find( mapEdge e );

  /**
   * TODO : RD
   * Write Desc
   */
  void swap( uint i, uint j );

  /**
   * TODO : RD
   * Write desc
   */
  void swap( mapEdge ei, mapEdge ej );

};

//============================== OPERATOR OVERLOAD  ==========================//
//================================= END IFDEF ================================//
#endif // HOCCURRENCESLIST_HPP
