/***************************************************************************
 *   Copyright (C) 2018 by Muaz Twaty                                      *
 *   muaz.sy123[at]gmail.com                                              
 * 
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
#ifndef MCTS_HPP
#define MCTS_HPP

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

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//
//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
/**
 * @brief The MCTS class
 * Class that store the information stored in a MCTS tree node
 */

class Searchtreenode
{
  //---- PUBLIC --------------------------------------------------------------//
public:

  // Public CONSTANTS ________________________________________________________//
  // Public variables ________________________________________________________//
  /// pointers to parents nodes
  long long nodeID;
  vector<Searchtreenode*> parents;
  // the valise of the N parameter of the MCTS algorithm
  int N_node;
  // the value of the Q function "the evaluation function"
  double Q;
  //
  bool is_fully_expanded;
  //
  bool occ_list_is_computed;
  
  GToken lastExt;
  
  // possible extentions
  map<GToken, Searchtreenode*, GTokenGt>* children_nodes;
  
  // occurance lists
  //GTokenData tokenData;

  // remaining coninical extenstions
  vector<pair<GToken, GExtensionData> > valid_extenstions;
  GTokenData      node_tokenData;

  

  // Public Structure & Typedef ______________________________________________//
  // Public Constructor/Desctructor __________________________________________//
  /**
   * @brief MCTS
   * Default constructor
   * Initialize variables and create new GVocab object
   */
  Searchtreenode();

  Searchtreenode(Searchtreenode* _parent_,GToken ext,long long ID);

  /**
   * @brief ~MCTS
   * Default Destructor
   * Free memory by deleting GVocab object and clearing vectors
   */
  ~Searchtreenode();

  // Accessor ________________________________________________________________//
  // Mutator _________________________________________________________________//
  // Public Methods __________________________________________________________//
  /**
   * TODO
   */
  
  
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

  
  
};

//============================== OPERATOR OVERLOAD  ==========================//
//================================= END IFDEF ================================//


#endif // MCTS_HPP
