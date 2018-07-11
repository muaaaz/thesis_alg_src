/***************************************************************************
 *   Copyright (C) 2005 by Siegfried Nijssen                               *
 *   snijssen[at]informatik.uni-freiburg.de                                *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2010,2011 by Adriana Prado and Baptiste Jeudy           *
 *   baptiste.jeudy[at]univ-st-etienne fr                                  *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2014 by Romain Deville                                  *
 *   romain.deville[a]insa-lyon.fr                                         *
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
 ***************************************************************************/

//=== INCLUDE ==================================================================
// Include library class
#include <unistd.h>
#include <iostream>
#include <ostream>
// Include project class
#include "goccurrenceslist.hpp"

//=============================== NAMESPACE ==================================//
using namespace std;
//============================ STRUCT & TYPEDEF ==============================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
//---- PUBLIC ----------------------------------------------------------------//
// Public CONSTANTS __________________________________________________________//
// Public Constructor/Desctructor ____________________________________________//
GOccurrencesList::GOccurrencesList():
  graphID(0)
{
  // Default Constructor
  data.resize(0);
}
// End of GOccurrencesList::GOccurrencesList()

GOccurrencesList::GOccurrencesList(uint graphId , GGraph *p_Graph)
{
  pGraph  = p_Graph;
  graphID = graphId;
  data.resize(0);
}
// End of GOccurrencesList::GOccurrencesList( uint graph ):

GOccurrencesList::~GOccurrencesList()
{
  /// Default destructor
  data.clear();
  data.shrink_to_fit();
}
// End of GOccurrencesList::~GOccurrencesList()

// Accessor __________________________________________________________________//
GOccurrencesList::mapEdge GOccurrencesList::at( uint i )
{
  /*
   * @brief atMat
   * TODO : RD
   * Coy desc from header
   * @param i
   * @return
   */
  return data.at(i);
}
// End of GOccurrencesList::atMat( uint i )

// End of GOccurrencesList::atDom( uint i )

// Mutator ___________________________________________________________________//

// End of GOccurrencesList::setSize( uint newSize )

// Public Methods ____________________________________________________________//
void GOccurrencesList::add( mapEdge edge )
{
  /*
   * TODO : RD
   * Copy desc from head
   */
  // Insert edge
  edge.element = size();
  data.push_back( edge );
  // update size
}
// End of GOccurrencesList::add( mapEdge edge )

void GOccurrencesList::add( GNodeID from, GNodeID dest, GEdgeID edge )
{
  /*
   * TODO : RD
   * Copy Desc
   */
  // Instanciate object mapEdge
  mapEdge e;
  e.nodeDest = dest;
  e.nodeFrom = from;
  e.edgeId   = edge;
  e.element  = size();
  
  data.push_back(e);
  // update size
}
// End of add(GNodeID from, GNodeID dest, GEdgeID edge )


void GOccurrencesList::remove( mapEdge e )
{
  /*
   * TODO : RD
   * Copy desc
   */
  int idx = find( e );
  if ( idx != -1 )
  {
    swap( idx, size()-1 );
    data.pop_back();
  }
}
// End of GOccurrencesList::remove( mapEdge e )

void GOccurrencesList::remove( uint i )
{
  /*
   * TODO : RD
   * Copy Desc
   */
  if(i < uint(size()) )
  {
    swap( i, size()-1 );
    data.pop_back();
  }
}
// End of GOccurrencesList::remove( int i )

//---- PROTECTED  ------------------------------------------------------------//
// Protected CONSTANTS _______________________________________________________//
// Protected Methods _________________________________________________________//

//---- PRIVATE ---------------------------------------------------------------//
// Private CONSTANTS _________________________________________________________//
// Private Methods ___________________________________________________________//
int GOccurrencesList::find( mapEdge e )
{
  /*
   * TODO : RD
   * Copy Desc
   */
  if ( size() == 0 )
    return -1;
  else
    for ( uint i=0 ; i < uint(size()); ++i )
      if ( data[i].nodeDest == e.nodeDest && data[i].nodeFrom == e.nodeFrom )
        return i;
  // If not find, return -1
  return -1;
}
// End of GOccurrencesList::find( mapEdge e )

void GOccurrencesList::swap( uint i, uint j )
{
  /*
   * TODO : RD
   * Copy Desc
   */
  mapEdge tmp = mapEdge(data[i]);
  data[i] = mapEdge(data[j]);
  data[j] =  mapEdge(tmp);

}
// End of GOccurrencesList::swap( uint i, uint j )

void GOccurrencesList::swap( mapEdge ei, mapEdge ej )
{
  /*
   * TODO : RD
   * Copy desc
   */
  uint i = find( ei );
  uint j = find( ej );
  swap( i, j );
}

void GOccurrencesList::clear()
{
  /*
   * TODO : RD
   * Copy desc
   */
  data.clear();
  data.shrink_to_fit();
}

//============================== OPERATOR OVERLOAD  ==========================//
