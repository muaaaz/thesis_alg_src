/***************************************************************************
 * Copyright (C) 2016 Deville Romain                                       *
 * Contact : romain.deville[at]insa-lyon.fr                                   *
 *-------------------------------------------------------------------------*
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, see <http://www.gnu.org/licenses/>   *
 *   or you can write to the                                               *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 ***************************************************************************/

//=== IFDEF ====================================================================
#ifndef GLOBAL_HPP
#define GLOBAL_HPP

//================================ INCLUDE ===================================//
// Include library class
// -- C++ STD LIB
#include <string>
#include <vector>
#include <math.h>
// Include project class
//=============================== NAMESPACE ==================================//
using namespace std;

//========================== STRUCT & TYPEDEF ================================//
/// Shortcut for unsigned int
typedef unsigned int uint;

/**
 * @brief Store with option was called by user
 */
struct param {
  /// Store input file
  string INFILE; // REQUIRED
  /// Store output directory were files will be saved
  string OUTDIR; // REQUIRED
  /// Store output pattern file were pattern will be saved
  string PATFILE;
  /// Store output pattern statistic file were pattern will be saved
  string PAT_STAT_FILE;
  /// Int to specify time limit for pattern
  int TIMEOUT;
  /// Store test file if specified
  bool TORGRAPH;
  /// Store boolean to know if grid will have hole
  bool INDIRECT_HOLEMINING;
  /// Store boolean to know if grid will have hole
  bool DIRECT_HOLEMINING;
  /// Store boolean to know if we mine with node-induced
  bool INDUCED;
  /// Specify that edge does not have label
  bool NOEDGELABEL;
  bool KEEPHOLE;
  /// Int to specify time limit to read for the graph
  int NBTIMETOREAD;
  /// Int to specify time limit for pattern
  int TIMELIMIT;
  /// Int to specify node number limit for pattern
  int NODELIMIT;
  /// Int to specify maximum number of pattern to min
  int NBPATLIMIT;
  /// Int to specify temporal dimension to start
  int INITIALTIME;
  ///
  int FRAMEHOLE;
  /// Store minimum frequency threshold
  float MINFREQ;
  /// Show debug informations
  bool DEBUG_MODE;
  /// Show warning informations
  bool WARNING_MODE;
  ///
  int MAXLBL;
  bool TEMP_FIRST;
  bool SPAT_FIRST;

  //MCTS
  //
  double C_p;


  param():
    INFILE(""),
    OUTDIR(""),
    PATFILE(""),
    PAT_STAT_FILE(""),
    TIMEOUT(-1),
    TORGRAPH(false),
    INDIRECT_HOLEMINING(false),
    DIRECT_HOLEMINING(false),
    INDUCED(false),
    NOEDGELABEL(false),
    KEEPHOLE(false),
    NBTIMETOREAD(-1),
    TIMELIMIT(-1),
    NODELIMIT(-1),
    NBPATLIMIT(-1),
    INITIALTIME(0),
    FRAMEHOLE(false),
    MINFREQ(0.5),
    DEBUG_MODE(false),
    WARNING_MODE(false),
    MAXLBL(-1),
    TEMP_FIRST(false),
    SPAT_FIRST(true),
    C_p(0.5)

  {
    // Nothing to put here
  }

  ~param()
  {
    // Nothing to put here as variable deletion are already managed by C++ or std
  }
};
// End of struct param

//=============================== VARIABLES ==================================//
/// Boolean to know if output debug
extern struct param PARAM;

//================================ METHODS ===================================//
//=============================== CONSTANTS ==================================//
//=========================== OPERATOR OVERLOAD  =============================//
//================================= END IFDEF ================================//
#endif
