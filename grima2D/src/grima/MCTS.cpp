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

//================================ INCLUDE ===================================//
// Include library class
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <algorithm>
// Include project class
#include "MCTS.hpp"
#include "global.hpp"

//=============================== NAMESPACE ==================================//
using namespace std;

//============================ STRUCT & TYPEDEF ==============================//

//=============================== VARIABLES ==================================//
// Will be used in cmpGToken when extending spatFirst or tempFirst


//================================= CLASS ====================================//
//---- PUBLIC ----------------------------------------------------------------//
// Public CONSTANTS __________________________________________________________//
// Public Constructor/Desctructor ____________________________________________//


MCTS_node::MCTS_node():
  parent(NULL),
  N_node(0),
  Q(-1),
  is_fully_expanded(false),
  occ_list_is_computed(false)
{
  /*
   * TODO : RD
   * Copy Desc
   */
  children_nodes    = new map<GToken, MCTS_node*,     GTokenGt>();
  valid_extenstions = vector<pair<GToken, GExtensionData> > ();
}

MCTS_node::MCTS_node(MCTS_node* _parent_):
  parent(_parent_),
  N_node(0),
  Q(-1),
  is_fully_expanded(false),
  occ_list_is_computed(false)
{
  /*
   * TODO : RD
   * Copy Desc
   */
  children_nodes    = new map<GToken, MCTS_node*,     GTokenGt>();
  valid_extenstions = vector<pair<GToken, GExtensionData> > ();
}



MCTS_node::~MCTS_node()
{
  /*
   * TODO : RD
   * Copy Desc
   */

  node_tokenData.v_SparseOcc.clear();
  valid_extenstions.clear();
  delete children_nodes;
}
