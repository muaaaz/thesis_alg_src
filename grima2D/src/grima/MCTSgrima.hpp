/***************************************************************************
 *   Copyright (C) 2005 by Siegfried Nijssen                               *
 *   snijssen@informatik.uni-freiburg.de                                   *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2010,2011 by Adriana Prado and Baptiste Jeudy           *
 *   baptiste.jeudy at univ-st-etienne fr                                  *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2015 by Romain Deville                                  *
 *   romain.deville[at]insa-lyon.fr                                        *
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
#ifndef MCTSGRIMA_HPP
#define MCTSGRIMA_HPP

//================================ INCLUDE ===================================//
// Include library class
#include "time.h"
// Include project class
#include "gglobal.hpp"
#include "ggraph.hpp"
#include "gdatabase.hpp"
#include "gvocab.hpp"
#include "gpattern.hpp"
#include "gextensioncollect.hpp"
#include "gsubgraphiso.hpp"
#include "gvocab.hpp"

#include "MCTS.hpp"



//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//
//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
/**
 * @brief The G class
 * Class that store vocabulary of pattern and allow to process grima mining
 * algorithm for each class.
 */
class MCTSGrima
{
  //---- PUBLIC --------------------------------------------------------------//
public:
  // Public CONSTANTS ________________________________________________________//
  // Public variables ________________________________________________________//
  /// Vocabulary of pattern
  GVocab* vocabPattern;
  /// Vector to store number of pattern for each class
  vector<int> v_NbPatternByClass;
  vector<GPattern *> v_PatternCurrClass;
  vector<int> v_ReturnStatus;

  /// Vector to class name
  vector<string> v_ClassName;
  /// Numb of pattern that belong to all class
  int nbPatternAllClass;
  /// Number of closed pattern, i.e. pattern that have over-pattern
  /// (i.e. pattern with at least one more edge) which have exacly same freq
  /// and same occurency
  int nbClosedPat;
  int nbTotalClosedPat;
  vector<int> v_nbClosedPatClass;
  /// Current idx of the class that is mined
  /// Id of the last frequent pattern that will be incremented for each new pattern
  int freqPatternId;
  /// First clock tick in case of TIMEOUT Specified
  clock_t firstTick;
  /// Total time to mine
  clock_t totalTick;
  /// Canonical time
  clock_t canonicalTick;
  /// Extension time
  clock_t mappingExtTick;
  clock_t extensionTick;
  /// Subgraphiso time
  clock_t subgraphisoTick;


  float minF;
  GClassDB* pClassDB;
  // Public Structure & Typedef ______________________________________________//
  // Public Constructor/Desctructor __________________________________________//
  /**
   * @brief Grima
   * Default constructor
   * Initialize variables and create new GVocab object
   */
  MCTSGrima();

  MCTSGrima(float _minf,GClassDB* _pClassDB);
  /**
   * @brief ~G
   * Default Destructor
   * Free memory by deleting GVocab object and clearing vectors
   */
  ~MCTSGrima();

  // Accessor ________________________________________________________________//
  // Mutator _________________________________________________________________//
  // Public Methods __________________________________________________________//
  /**
   * TODO
   */
  void initNbPatternByClass(vector<GClassDB *> v_GClassDB );

  /**
   * @brief processMining
   * Method that apply mining process algorithm by calling search trough all
   * images withing the same class in pClassDB for a user set min support.
   * @param minF
   * @param pClassDB
   * @param classIdx
   */
  int processMining( );

  /**
   * @brief saveData
   * TODO
   * @param filename
   * @param returnStatus
   */
  void saveData(bool timeOutOverride );


  inline double UCB(MCTS_node* cur, MCTS_node* child);

  MCTS_node* best_child(MCTS_node* cur);

  MCTS_node* select(MCTS_node* cur); // static??

  MCTS_node* expand(MCTS_node* cur);

  double rool_out(MCTS_node* cur, 
                    bool first_level,
                    vector<GGraph*>* v_Graphs,
                    GGlobFreq       &minFreq,    //Mininmum global frequency
                    GPattern*        pPattern,
                    GTokenData      &tokenData,
                    GExtensionData  &suppData,   // Tmp variable, supposed frequency
                    GExtensionData  &prevData);

  void update_ancestors(MCTS_node* cur, double delta); 

  //---- PROTECTED  ----------------------------------------------------------//
protected:
  // Protected CONSTANTS _____________________________________________________//
  // Protected Structure & Typedef ___________________________________________//
  // Protected Variables _____________________________________________________//
  // Protected Methods _______________________________________________________//

  //---- PRIVATE -------------------------------------------------------------//
private:
  // Private CONSTANTS _______________________________________________________//
  // Private Structure & Typedef _____________________________________________//
  // Private Variables _______________________________________________________//
  // Private Methods _________________________________________________________//
  /**
   * TODO : RD
   * Write Desc
   */
  int search( vector<GGraph*>                   &v_Graphs,
              map<GToken, GTokenData, GTokenGt> &m_TokenData,
              GGlobFreq                         minFreq);

  
};

//============================== OPERATOR OVERLOAD  ==========================//
//================================= END IFDEF ================================//
#endif // MCTSGRIMA_HPP