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

//================================ INCLUDE ===================================//
// Include library class
// Include project class
#include "greader.hpp"

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//
//=============================== VARIABLES ==================================//
//================================= CLASS ====================================//
//---- PUBLIC ----------------------------------------------------------------//
// Public CONSTANTS __________________________________________________________//

// Public Constructor/Desctructor ____________________________________________//
GReader::GReader( const string& fn )
{
  /*
   * @brief GReader Default constructor
   * Open file fn, test if file is open, if not exit with return code 1
   * @param fn : File name where graph are describe
   */
  in.open( fn.c_str(), ifstream::in );
  if ( !in.is_open() )
    error::openError( fn );
  out     = 0;
  line_nb = 0;
}
// End of GReader::GReader( const string& fn )

GReader::GReader(const string& fn, ostream *stream):
  out(stream)
{
  /*
   * @brief GReader Overloaded constructor
   * Open file fn, test if file is open, if not exit with return code 1
   * @param fn : File name where graph are describe
   * @param stream
   */
  in.open( fn.c_str(), ifstream::in );
  if ( !in.is_open() )
    error::openError( fn );
  line_nb = 0;
}
// End of GReader::GReader(const string& fn, ostream *stream)

GReader::~GReader()
{
  /*
   * @brief ~GReader Default Destructor
   * Close file fn
   */
  in.close() ;
}
// End of GReader::~GReader()

// Accessor __________________________________________________________________//
int GReader::get_line_nb()
{
  /*
   * Function that return the current number of line
   */
  return line_nb;
}
// End of GReader::get_line_nb()

bool GReader::finished()
{
  /*
   * @brief finished
   * Check if the file is completely read, hence, read is finished
   * @return TRUE if EOF, FALSE else.
   */
  return in.eof();
}
// End of GReader::finished()

// Mutator ___________________________________________________________________//
// Public Methods ____________________________________________________________//
void GReader::getLine( vector<string>* tokens )
{
  /*
   * @brief getLine
   * Get next line of file and update the vector of string tokens (ie. a vector
   * with each chars separeted by space in input file.
   * @param tokens  Vector of string that will handle char.
   */

  // Max size of line
  const int s = 2048;
  // Create char of size s
  char line[s];

  tokens->clear();
  do {
    do {
      in.getline(line, s);
      // If end of file or line greater than s
      if ( in.fail() && !in.eof() )
        error::nbLineError( line_nb, s );
      line_nb++;
      // If comment line
      if (line[0] == '#')
      {
        if (out)
          *out << line << endl;
      }
    } while (line[0] == '#');

    char* token = NULL;
    char* ptr = line;

    // Parse line
    token = strtok(ptr, " ");
    while (token != NULL)
    {
      tokens->push_back(token);
      token = strtok(NULL, " ");
    }
  } while (tokens->size() == 0 && !in.eof());
}

//---- PROTECTED  ------------------------------------------------------------//
// Protected CONSTANTS _______________________________________________________//
// Protected Structure _______________________________________________________//
// Protected Variables _______________________________________________________//
// Protected Methods _________________________________________________________//

//---- PRIVATE ---------------------------------------------------------------//
// Private CONSTANTS _________________________________________________________//
// Private Structure _________________________________________________________//
// Private Variables _________________________________________________________//
// Private Methods ___________________________________________________________//

//============================== OPERATOR OVERLOAD  ==========================//
