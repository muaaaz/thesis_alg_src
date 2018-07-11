/***************************************************************************
 *   Copyright (C) 2005 by Siegfried Nijssen                               *
 *   snijssen[at]informatik.uni-freiburg.de                                *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2010,2011 by Adriana Prado and Baptiste Jeudy           *
 *   baptiste.jeudy[at]univ-st-etienne fr                                  *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2016 by Romain Deville                                  *
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
#ifndef DATA_HPP
#define DATA_HPP

//================================ INCLUDE ===================================//
// Include library class
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
// Include project class

//=============================== NAMESPACE ==================================//
using namespace std;

//========================== STRUCT & TYPEDEF ================================//
//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//
//================================= CLASS ====================================//
/**
 * @brief GReader class
 * Handle reading line for the file that store graph, split line into token
 * to be parsed later
 */
class GReader
{
  //---- PUBLIC --------------------------------------------------------------//
public:
  // Public CONSTANTS ________________________________________________________//
  // Public Structure & Typedef ______________________________________________//
  // Public Variables ________________________________________________________//
  // Public Constructor/Desctructor __________________________________________//
  /**
   * @brief GReader Default constructor
   * Open file fn, test if file is open, if not exit with return code 1
   * @param fn : File name where graph are describe
   */
  GReader( const string& fn );
  /**
   * @brief GReader Overloaded constructor
   * Open file fn, test if file is open, if not exit with return code 1
   * @param fn : File name where graph are describe
   * @param stream
   */
  GReader( const string& fn, ostream *stream );

  /**
   * @brief ~GReader Default Destructor
   * Close file fn
   */
  ~GReader();

  // Accessor ________________________________________________________________//
  /**
   * @brief get_line_nb
   * Return line_nb
   * @return line_nb
   */
  int get_line_nb(); // Added by Baptiste

  /**
   * @brief getLine
   * Get next line of file and update the vector of string tokens (ie. a vector
   * with each chars separeted by space in input file.
   * @param tokens  Vector of string that will handle char.
   */
  void getLine( vector<string>* tokens );

  // Mutator _________________________________________________________________//
  // Public Methods __________________________________________________________//
  /**
   * @brief finished
   * Check if the file is completely read, hence, read is finished
   * @return TRUE if EOF, FALSE else.
   */
  bool finished();

  //---- PROTECTED  ----------------------------------------------------------//
protected:
  // Protected CONSTANTS _____________________________________________________//
  // Protected Structure _____________________________________________________//
  // Protected Variables _____________________________________________________//
  // Protected Methods _______________________________________________________//

  //---- PRIVATE -------------------------------------------------------------//
private:
  // Private Structure _______________________________________________________//
  struct error {
    /**
     * @brief openError
     * Print the message when unable to open file fn
     * @param fn : File unable to be open
     */
    static void openError( string fn )
    {
      cerr << "ERROR - Unable to open input file " << fn << endl;
      exit( EXIT_FAILURE );
    }
    /**
     * @brief nbLineError
     * Print an error when the line read is too long
     * @param line_nb : Line that raise issue
     * @param s       : Maximum size of that line
     */
    static void nbLineError( int line_nb, int s )
    {
      cerr << "Error input line "<< line_nb << " (line too long ?)" << endl;
      cerr << "Max line length: " << s << endl;
      exit ( EXIT_FAILURE );
    }
  };

  // Private CONSTANTS _______________________________________________________//
  // Private Variables _______________________________________________________//
  /// Nb line in the file
  int      line_nb; // Added by Baptiste
  /// Input File Stream
  ifstream in;
  /// Output File stream
  ostream *out;

  // Private Methods _________________________________________________________//
};

//============================== OPERATOR OVERLOAD  ==========================//
//================================= END IFDEF ================================//
#endif
