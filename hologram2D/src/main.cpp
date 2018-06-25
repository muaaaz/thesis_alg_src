/***************************************************************************
 *   Copyright (C) 2016 by Romain Deville                                  *
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

//================================ INCLUDE ===================================//
// Include library class
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
// Include project class
#include "global.hpp"

#include "hologram/gdatabase.hpp"
#include "hologram/hologram.hpp"


// for testing mod
// #define TESTCAN

#ifdef TESTCAN

#include "hologram/gpattern.hpp"
#include "hologram/gcanonicaldfscomputer.hpp"
#include "hologram/gglobal.hpp"

#endif


//=============================== NAMESPACE ==================================//
using namespace std;

//============================ STRUCT & TYPEDEF ==============================//
//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//
/**
 * @brief help
 * Print usage of the program if user input is invalid or directly called by
 * user
 */
void help()
{
  /*
   * TODO : RD
   * Write usage
   */
  cout << "TODO : RD " << endl
       << "Write HELP !!!" << endl;
  cout << "{           \"debug\",       no_argument, 0, \'D\'},"<< endl <<
          "{         \"warning\",       no_argument, 0, \'W\'},"<< endl <<
          "{(in Hour)\"timeout\",       no_argument, 0, \'T\'},"<< endl <<
          "{             \"tor\",       no_argument, 0, \'t\'},"<< endl <<
          "{            \"hole\",       no_argument, 0, \'h\'},"<< endl <<
          "{         \"induced\",       no_argument, 0, \'i\'},"<< endl <<
          "{     \"noedgelabel\",       no_argument, 0, \'n\'},"<< endl <<
          "{      \"gameOfLife\",       no_argument, 0, \'g\'},"<< endl <<
          "{ \"gameOfLifeAlive\",       no_argument, 0, \'a\'},"<< endl <<
          "{  \"gameOfLifeDead\",       no_argument, 0, \'d\'},"<< endl <<
          "{   \"sizeTimeLimit\", required_argument, 0, \'S\'},"<< endl <<
          "{     \"nbNodeLimit\", required_argument, 0, \'N\'},"<< endl <<
          "{      \"nbPatLimit\", required_argument, 0, \'P\'},"<< endl <<
          "{\"initialTimeStart\",       no_argument, 0, \'I\'},"<< endl <<
          "{  \"sizeTimeToRead\", required_argument, 0, \'s\'},"<< endl <<
          "{       \"frequency\", required_argument, 0, \'f\'},"<< endl <<
          "{         \"pattern\", required_argument, 0, \'p\'},"<< endl;
}
// End of help()

/**
 * @brief parse argument give by user
 * @param argc: # of argument
 * @param argv: vector of argument
 * @return 0 if no error, 1 otherwise
 */
void parseArg( int argc, char **argv, param &PARAM )
{
  string   infile;
  int      c;

  while ( true )
  {
    static struct option long_options[] =
    {
      // These options are distinguish by their indices.
    {           "debug",       no_argument, 0, 'D'},
    {         "warning",       no_argument, 0, 'W'},
    {         "timeout",       no_argument, 0, 'T'},
    {             "tor",       no_argument, 0, 't'},
    {    "indirecthole",       no_argument, 0, 'h'},
    {      "directhole",       no_argument, 0, 'H'},
    {         "induced",       no_argument, 0, 'i'},
    {     "noedgelabel",       no_argument, 0, 'n'},
    {   "sizeTimeLimit", required_argument, 0, 'S'},
    {     "nbNodeLimit", required_argument, 0, 'N'},
    {      "nbPatLimit", required_argument, 0, 'P'},
    {    "keepNodeHole",       no_argument, 0, 'K'},
    {"initialTimeStart", required_argument, 0, 'I'},
    {  "sizeTimeToRead", required_argument, 0, 's'},
    {       "frequency", required_argument, 0, 'f'},
    {       "frameHole",       no_argument, 0, 'F'},
    {       "outputDir", required_argument, 0, 'o'},
    {         "pattern", required_argument, 0, 'p'},
    {     "patternStat", required_argument, 0, 'x'},
    {  "exclude100Freq",       no_argument, 0, 'X'},
    {         "classID", required_argument, 0, 'c'},
    {            "help",       no_argument, 0,  0 },
    {                 0,                 0, 0,  0 }
  };
    // getopt_long stores the option index here.
    int option_index = 0;
    c = getopt_long( argc, argv, "DWthXHinKT:S:N:P:I:s:f:Fo:p:x:", long_options, &option_index );
    // Detect the end of the options.
    if ( c == -1 )
      break;
    switch ( c )
    {
    case 0:
      if ( (string)long_options[option_index].name == "help" )
      {
        help();
        exit( EXIT_SUCCESS );
      }
      break;
    case 'D':
      PARAM.DEBUG_MODE = true;
      break;
    case 'W':
      PARAM.WARNING_MODE = true;
      break;
    case 't':
      PARAM.TORGRAPH = true;
      break;
    case 'h':
      PARAM.INDIRECT_HOLEMINING = true;
      break;
      case 'H':
        PARAM.DIRECT_HOLEMINING = true;
        break;
    case 'i':
      PARAM.INDUCED  = true;
      break;
    case 'n':
      PARAM.NOEDGELABEL = true;
      break;
    case 'S':
      if ( atoi(optarg) > 0 )
        PARAM.TIMELIMIT = atoi(optarg);
      break;
    case 's':
      PARAM.NBTIMETOREAD = atoi(optarg);
      break;
    case 'N':
      PARAM.NODELIMIT = atoi(optarg);
      break;
    case 'P':
      PARAM.NBPATLIMIT = atoi(optarg);
      break;
    case 'T':
      PARAM.TIMEOUT = atoi(optarg);
      break;
    case 'I':
      PARAM.INITIALTIME = atoi(optarg);
      break;
    case 'f':
      PARAM.MINFREQ  = atof(optarg);
      break;
    case 'F':
      PARAM.FRAMEHOLE = true;
      break;
    case 'K':
      PARAM.KEEPHOLE = true;
      break;
    case 'o':
      PARAM.OUTDIR  = optarg;
      break;
    case 'p':
      PARAM.PATFILE = optarg;
      break;
    case 'x':
      PARAM.PAT_STAT_FILE = optarg;
      break;
    default:
      help();
      exit( EXIT_FAILURE );
    }
  }
#ifdef _DEBUG
  PARAM.DEBUG_MODE = true;
  PARAM.WARNING_MODE = true;
#endif
  // Print help if to wrong number of arg
  if ( argc - optind != 1 )
  {
    help();
    exit( EXIT_FAILURE );
  }

  // Check if input file exist and is openable, if not exist, else set
  // INPUTFILE
  infile = argv[optind];
  if ( ! ifstream( infile ) )
  {
    cerr << "ERROR - file " << infile << " does not exist ." << endl;
    help();
    exit ( EXIT_FAILURE );
  }
  else
    PARAM.INFILE=infile;

  // Check if output dir is define
  if ( PARAM.OUTDIR == "" )
  {
    cerr << "ERROR - Output dir MUST be specified with option -o <outDir>" << endl;
    cerr << "Pattern will then be print in stdout" << endl;
    help();
    exit( EXIT_FAILURE );
  }
  else if ( PARAM.OUTDIR.back() != '/' )
    PARAM.OUTDIR += '/';

  // Check if pattern file is define, if not print a warning with default path
  if ( PARAM.PATFILE == "" )
  {
    PARAM.PATFILE = "pattern.gri";
    cerr << "WARNING - Output pattern file is not be specified with option -p <fileName>" << endl;
    cerr << "Pattern will then be print to " << PARAM.OUTDIR + PARAM.PATFILE << endl;
  }

  // Check if pattern file is define, if not print a warning with default path
  if ( PARAM.PAT_STAT_FILE == "" )
  {
    PARAM.PAT_STAT_FILE = "pattern_stat.csv";
    cerr << "WARNING - Output pattern stat file is not be specified with option -x <fileName>" << endl;
    cerr << "Pattern will then be print to " << PARAM.OUTDIR + PARAM.PAT_STAT_FILE << endl;
  }

  if ( PARAM.DIRECT_HOLEMINING == true && PARAM.INDIRECT_HOLEMINING == true )
  {
    cerr << "ERROR - Indirect and direct hole are set. Please choose one with option : " << endl
         << "-H : For Direct Hole" << endl
         << "-h : For indirect Hole" << endl;
    exit( EXIT_FAILURE ) ;
  }

}
// End of parseArg( int argc, char **argv )

// ===== CONSTANT DEFINITION =================================================//
param PARAM = param();

//============================================================================//
//                                  MAIN                                      //
//============================================================================//
/**
  * @brief Main function
  * Read graph file
  *   - Determine number of graph
  *   - Determine max number of label for color
  * IF one outfile (i.e. no option '-s')
  *   - Write latex header to outfile.tex
  *   - Write color definition to outfile.tex
  *   - Write one graph per page in tikz format resize to fit page
  * ELSE
  *   - Write latex header to outfile-header.tex
  *   - Write color definition to outfile-header.tex
  *   FOR each graph
  *     - Write one graph in tikz format resize to fit page in outfile-i.tex
  *     - Write an input to this file to outfile-header.tex
  * @param argc # of arg
  * @param argv vector of arg
  * @return 0 if no issu, 1 otherwise
*/

// ./hologram2D -o /udd/mtwaty/Downloads/thesis_alg_src/hologram2D/output -T 1 -f 0.5 in.gri

int incremental_counter = 0;

int main( int argc, char **argv )
{
  srand(time(0));
  cout << "#==== START OF Hologram !" << endl;
  parseArg( argc, argv, PARAM );

  cerr <<"input file "<< PARAM.INFILE << endl;
  cerr <<"output file "<< PARAM.OUTDIR << endl;
  cerr <<"time "<< PARAM.TIMEOUT << endl;
  cerr <<"freq "<< PARAM.MINFREQ << endl;

  GDatabase graphDB;
  graphDB.createGrapheDB( PARAM.INFILE );
  
  Hologram hologram(PARAM.MINFREQ,graphDB.v_GClassDB.at(0),graphDB.v_ClassName);

  hologram.initialize();
  for ( int iClassDB = 1; iClassDB <= hologram.pClassDB->number_of_classes; iClassDB++ )
  {
    hologram.set_class_id(iClassDB);
    hologram.processMining();
    hologram.clean();
  }
  
  hologram.saveData();
  cout << "#==== END OF Hologram !" << endl;
  exit(0);
}
// End of main( int argc, char **argv )
