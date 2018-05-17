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
#ifndef GDATABASE_HPP
#define GDATABASE_HPP

//================================ INCLUDE ===================================//
// Include library class
#include <vector>
#include <map>
// Include project class
#include "../global.hpp"
#include "gglobal.hpp"
#include "greader.hpp"
#include "ggraph.hpp"
#include "gsparseset.hpp"

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYPEDEF ==============================//

/**
 * @brief cmpGToken
 * TODO : RD
 * Write Desc
 * @param t1
 * @param t2
 * @return
 */
inline int cmpGTokenCanonTest ( const GToken &t1, const GToken &t2 )
{
  if ( t1.direction != t2.direction )
    return t1.direction - t2.direction;
  else if ( t1.nodeFrom != t2.nodeFrom )
    return t1.nodeFrom - t2.nodeFrom;
  else if ( t1.nodeDest != t2.nodeDest )
    return t1.nodeDest - t2.nodeDest;
  else if ( t1.angle != t2.angle )
    return t1.angle - t2.angle;
  else if ( t1.nodeLabelFrom != t2.nodeLabelFrom )
    return t1.nodeLabelFrom - t2.nodeLabelFrom;
  else if ( t1.edgeLabel != t2.edgeLabel )
    return t1.edgeLabel - t2.edgeLabel;
  else if ( t1.nodeLabelDest != t2.nodeLabelDest )
    return t1.nodeLabelDest - t2.nodeLabelDest;
  else
    return 0;

//  return memcmp ( (void*) &t1, (void*) &t2,
//                  sizeof( GTokenDirection )
//                  + sizeof( GNodeID ) * 2
//                  + sizeof( GEdgeAngle )
//                  + sizeof( GNodeLabel ) * 2
//                  + sizeof( GEdgeLabel ) );
}

inline int cmpGToken ( const GToken &t1, const GToken &t2 )
{
  tmpT1.direction = t1.direction;
  tmpT1.nodeFrom = t1.nodeFrom ;
  tmpT1.nodeDest = t1.nodeDest;
  tmpT1.nodeLabelFrom = t1.nodeLabelFrom;
  tmpT1.nodeLabelDest = t1.nodeLabelDest;
  tmpT1.edgeLabel = t1.edgeLabel;


  tmpT2.direction = t2.direction;
  tmpT2.nodeFrom = t2.nodeFrom ;
  tmpT2.nodeDest = t2.nodeDest;
  tmpT2.nodeLabelFrom = t2.nodeLabelFrom;
  tmpT2.nodeLabelDest = t2.nodeLabelDest;
  tmpT2.edgeLabel = t2.edgeLabel;

  if ( PARAM.SPAT_FIRST )
  {
    tmpT1.angle = t1.angle + 2;
    tmpT2.angle = t2.angle + 2;
  }
  else if ( PARAM.TEMP_FIRST )
  {
    tmpT1.angle = t1.angle + 0;
    tmpT2.angle = t2.angle + 0;
  }

  return memcmp ( (void*) &tmpT1, (void*) &tmpT2,
                  sizeof( GTokenDirection )
                  + sizeof( GNodeID ) * 2
                  + sizeof( GEdgeAngle )
                  + sizeof( GNodeLabel ) * 2
                  + sizeof( GEdgeLabel ) );
}

/**
 * TODO : RD
 * Write Desc
 */
struct GTokenGt {
  /**
   * @brief operator ()
   * TODO : RD
   * Write desc
   */
  bool operator()( const GToken &t1, const GToken &t2 ) const
  {
    return  cmpGToken ( t1, t2 ) > 0 ;
  }
};

/**
 * @brief The GTokenData struct
 * TODO : RD
 * Write Desc
 */
struct GTokenData {
  /// List of all occurence through a vector of sparseset, one sparset per grid
  vector<GSparseSet> v_SparseOcc;
  /// Absolute frequency of edge ( <= that the number of graph )
  int size()
  {
    int ret = 0;
    for ( int i=0;i<v_SparseOcc.size();++i)
      ret += v_SparseOcc[i].data.size();
    
    return ret;
  }
  uint               freq;
};

struct D_GTokenData {
  
  GTokenData token;
  vector<bool> is_difference_v;
};

/**
 * @brief GClassDB
 * Store graph by classname.
 */
struct GClassDB {
  /// Name of the class
  string className;
  /// Vector that store all graph of the class
  vector<GGraph*> v_ClassGraphs;
  /// List of all 1-edge canonical DFS code occure
  map<GToken,GTokenData,GTokenGt> m_TokenData;

  /// Default constructor
  GClassDB():
    className("")
  {
    v_ClassGraphs.resize(0);
  }

  /// Default desctructor
  ~GClassDB()
  {
    for ( uint i = 0; i < v_ClassGraphs.size(); i++ )
      delete v_ClassGraphs[i];
    v_ClassGraphs.clear();
    m_TokenData.clear();
  }
};
// End of struct GClassDB

//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
/**
 * @brief The GDatabase class
 * Class that handle the graph file, read it and process it to create database
 */
class GDatabase
{
  //---- PUBLIC --------------------------------------------------------------//
public:
  // Public CONSTANTS ________________________________________________________//
  // Public Structure & Typedef ______________________________________________//
  // Public Variables_________________________________________________________//
  /// Vector of GClassDB object, ie each class with their graph
  vector<GClassDB*> v_GClassDB;

  // Public Constructor/Desctructor __________________________________________//
  /// Default constructor
  GDatabase();

  /// Default Desctructor
  ~GDatabase();

  // Accessor ________________________________________________________________//
  // Mutator _________________________________________________________________//
  // Public Methods __________________________________________________________//
  /**
   * @brief createGrapheDB
   * Load graphFile into v_GClassDB.
   * If hole = true, make hole into graph, i.e. if all neighbors of a node have
   * the same label, then this node is deleted (set to GNONODEID).
   * @param graphFile : path to file that store graph
   * @param hole      : boolean to specify if it's needed to make hole in graph
   */
  void createGrapheDB(string &graphFile);

  // ----- PROTECTED  ----------------------------------------------------------//
protected:
  // Protected CONSTANTS _____________________________________________________//
  // Protected Variables _____________________________________________________//
  // Protected Methods _______________________________________________________//


  // ----- PRIVATE -----------------------------------------------------------//
private:
  // Private CONSTANTS _______________________________________________________//
  // Private Structure _______________________________________________________//
  /**
   * @brief The error struct
   * Structure with multiple parameterize method to pring error message
   */
  struct error {
    /**
     * @brief fileNotExistError
     * Print a message when file does not exist, thus unable to be open
     * @param file : File that raise issue
     */
    static void fileNotExistError( string file )
    {
      cerr << "ERROR - File " << file
           << " does not exist. Can't be read to create graph database" << endl;
      exit( EXIT_FAILURE );
    }
  };
  // Private Variables _______________________________________________________//
  // Private Methods _________________________________________________________//
  /**
   * @brief read
   * Read the file that store graph described in the filename
   * @param filename : File path where graph are stored
   * @return A vector of structure that store graph by class with class infos.
   */
  void read( string filename );

  /**
   * TODO : RD
   * Write Desc
   */
  void processSparseSet( GClassDB *pClassDB );

  /**
   * TODO : RD
   * Write Desc
   */
  void processEdgeSparseSet( GClassDB           *pClassDB,
                             uint               iGraph,
                             uint               iNode ,
                             vector<GGraphNode> &v_Nodes );

  /**
   * TODO : RD
   * Write Desc
   */
  void addEdgeSpatialSparseSet( GNodeLabel nodeLabelFrom,
                                GNodeLabel nodeLabelDest,
                                GEdgeLabel edgeLabel,
                                GClassDB   *pClassDB,
                                uint iGraph,
                                uint iNodeFrom,
                                uint iNodeDest ,
                                uint iEdge);

  void addEdgeTemporalSparseSet( GNodeLabel nodeLabelFrom,
                                 GNodeLabel nodeLabelDest,
                                 GEdgeLabel edgeLabel,
                                 GClassDB   *pClassDB,
                                 uint iGraph,
                                 uint iNodeFrom,
                                 uint iNodeDest ,
                                 uint iEdge);

  /**
   * @brief processHoleInGrid
   * Process hole in grid by calling isNodeToDeleted2D (if grid as tdim=1) or
   * isNodeToDelete2D_t (if grid as tdim > 1) for each node for each graph.
   * This will construct a vector of node "to delete" (i.e. set their label
   * to GNOLABEL and remove their edges).
   * Then, for each node "to delete", call setNodeToNoLabel.
   * @param pClassDB : Class database of grids
   */
  void processHoleInGrid( GClassDB* pClassDB );
};

//============================== OPERATOR OVERLOAD  ==========================//
//================================= END IFDEF ================================//
#endif // GDATABASE_H
