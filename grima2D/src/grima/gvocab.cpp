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

//================================ INCLUDE ===================================//
// Include library class
#include <vector>
#include <algorithm>
// Include prohect class
#include "global.hpp"
#include "gvocab.hpp"
#include "gsubgraphiso.hpp"

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//
//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
//---- PUBLIC --------------------------------------------------------------//
// Public CONSTANTS ________________________________________________________//
// Public Constructor/Desctructor __________________________________________//
GVocab::GVocab()
{
  // Default constructor
  // Nothing to do here
}
// End of GVocab::GVocab()

GVocab::~GVocab()
{
  // Default destructor
  for ( uint i = 0; i < v_AllPatterns.size(); i++ )
    delete v_AllPatterns.at(i);
  v_AllPatterns.clear();
  v_AllPatterns.shrink_to_fit();
  for ( uint i = 0; i < v_PatternByClass.size(); i++ )
    v_PatternByClass.at(i).clear();
  v_PatternByClass.clear();
  v_PatternByClass.shrink_to_fit();
  

}
// End of GVocab::~GVocab()

// Accessor ________________________________________________________________//
// Mutator _________________________________________________________________//
// Public Methods __________________________________________________________//
bool GVocab::comp( GPattern* pPat1, GPattern* pPat2 )
{
  /*
   * @brief comp
   * Method that compare 2 pattern (pPat1 and pPat2).
   * @param pPat1 ; First pattern to compare
   * @param pPat2 : Second pattern to compare
   * @return TRUE if pPat1 < pPat2, FALSE if pPat2 < pPat1
   */
  for ( uint i = 0 ; i < pPat1->v_Tokens.size(); i++)
  {
    if ( cmpGToken( pPat1->v_Tokens[i], pPat2->v_Tokens[i] ) != 0 )
      return false;
  }
  return true;
}
// End of GVocab::comp( GPattern* pPat1, GPattern* pPat2 )

void GVocab::saveVocab( string filename, vector<string> v_ClassName, vector<int> v_ReturnStatus )
{
  ofstream patFile;
  patFile.open( filename );
  for ( uint iClass = 0; iClass < v_ReturnStatus.size(); iClass++ )
  {
    if ( v_ReturnStatus.at(iClass) == -1  )
      patFile << "# INCOMPLETE MINING !! TIMEOUT REACHED (" << PARAM.TIMEOUT
              << "H) for class " << v_ClassName.at(iClass) << endl;
    else if ( v_ReturnStatus.at(iClass) == -2  )
      patFile << "# INCOMPLETE MINING !! NB PAT LIMIT REACHED (" << PARAM.NBPATLIMIT
              << ") for class " << v_ClassName.at(iClass) << endl;
  }
  patFile << "NbFreqPattern : " << v_AllPatterns.size() << endl;
  for ( uint iPat = 0; iPat < v_AllPatterns.size(); iPat++ )
    patFile << v_AllPatterns.at(iPat);
  patFile.close();

  for ( uint iClass = 0; iClass < v_PatternByClass.size(); iClass++ )
  {
    for ( uint iPat = 0; iPat < v_PatternByClass.at(iClass).size(); iPat++ )
      v_PatternByClass.at(iClass).at(iPat)->pGraph->graphID = iPat;

    patFile.open( filename + "." + v_ClassName.at(iClass) );
    if ( v_ReturnStatus.at(iClass) == -1  )
      patFile << "# INCOMPLETE MINING !! TIMEOUT REACHED (" << PARAM.TIMEOUT
              << "H) for class " << v_ClassName.at(iClass) << endl;
    else if ( v_ReturnStatus.at(iClass) == -2  )
      patFile << "# INCOMPLETE MINING !! NB PAT LIMIT REACHED (" << PARAM.NBPATLIMIT
              << ") for class " << v_ClassName.at(iClass) << endl;
    patFile << "NbFreqPattern : " <<  v_PatternByClass.at(iClass).size() << endl;
    for ( uint iPat = 0; iPat < v_PatternByClass.at(iClass).size(); iPat++ )
      patFile << v_PatternByClass.at(iClass).at(iPat);
    patFile.close();
  }
}
//============================== OPERATOR OVERLOAD  ==========================//
