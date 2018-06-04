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
#include "gsparseset.hpp"

//=============================== NAMESPACE ==================================//
using namespace std;
//============================ STRUCT & TYPEDEF ==============================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
//---- PUBLIC ----------------------------------------------------------------//
// Public CONSTANTS __________________________________________________________//
// Public Constructor/Desctructor ____________________________________________//
GSparseSet::GSparseSet():
  graphID(0)
{
  // Default Constructor
  data.resize(0);
}
// End of HSparseSet::SparseSet()

GSparseSet::GSparseSet(uint graphId , GGraph *p_Graph)
{
  /**
   * @brief GSparseSet
   * Overloaded constructor
   * @param graphId : Id of the graph to which is associated the sparseset
   */
  pGraph  = p_Graph;
  graphID = graphId;
  data.resize(0);
}
// End of GSparseSet::GSparseSet( uint graph ):

GSparseSet::~GSparseSet()
{
  /// Default destructor
  data.clear();
  data.shrink_to_fit();
}
// End of GSparseSet::~GSparseSet()

// Accessor __________________________________________________________________//
GSparseSet::mapEdge GSparseSet::at( uint i )
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
// End of GSparseSet::atMat( uint i )

// End of GSparseSet::atDom( uint i )

// Mutator ___________________________________________________________________//

// End of GSparseSet::setSize( uint newSize )

// Public Methods ____________________________________________________________//
void GSparseSet::add( mapEdge edge )
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
// End of GSparseSet::add( mapEdge edge )

void GSparseSet::add( GNodeID from, GNodeID dest, GEdgeID edge )
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


void GSparseSet::remove( mapEdge e )
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
// End of GSparseSet::remove( mapEdge e )

void GSparseSet::remove( uint i )
{
  /*
   * TODO : RD
   * Copy Desc
   */
  if(i < size())
  {
    swap( i, size()-1 );
    data.pop_back();
  }
}
// End of GSparseSet::remove( int i )

//---- PROTECTED  ------------------------------------------------------------//
// Protected CONSTANTS _______________________________________________________//
// Protected Methods _________________________________________________________//

//---- PRIVATE ---------------------------------------------------------------//
// Private CONSTANTS _________________________________________________________//
// Private Methods ___________________________________________________________//
int GSparseSet::find( mapEdge e )
{
  /*
   * TODO : RD
   * Copy Desc
   */
  if ( size() == 0 )
    return -1;
  else
    for ( uint i=0 ; i < size(); ++i )
      if ( data[i].nodeDest == e.nodeDest && data[i].nodeFrom == e.nodeFrom )
        return i;
  // If not find, return -1
  return -1;
}
// End of GSparseSet::find( mapEdge e )

void GSparseSet::swap( uint i, uint j )
{
  /*
   * TODO : RD
   * Copy Desc
   */
  mapEdge tmp = mapEdge(data[i]);
  data[i] = mapEdge(data[j]);
  data[j] =  mapEdge(tmp);

}
// End of GSparseSet::swap( uint i, uint j )

void GSparseSet::swap( mapEdge ei, mapEdge ej )
{
  /*
   * TODO : RD
   * Copy desc
   */
  uint i = find( ei );
  uint j = find( ej );
  swap( i, j );
}

void GSparseSet::clear()
{
  /*
   * TODO : RD
   * Copy desc
   */
  data.clear();
  data.shrink_to_fit();
}

//============================== OPERATOR OVERLOAD  ==========================//
