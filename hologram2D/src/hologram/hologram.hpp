/***************************************************************************
 *   Copyright (C) 2005 by Siegfried Nijssen                               *
 *   snijssen@informatik.uni-freiburg.de                                   *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2010,2011 by Adriana Prado and Baptiste Jeudy           *
 *   baptiste.jeudy at univ-st-etienne fr                                  *
 * ----------------------------------------------------------------------- *
 *   Copyright (C) 2015 by Romain Deville                                  *
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

//================================= IFDEF ====================================//
#ifndef HOLOGRAM_HPP
#define HOLOGRAM_HPP

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
#include "searchtreenode.hpp"
#include <set>
#include <unordered_map>

//=============================== NAMESPACE ==================================//
//============================ STRUCT & TYP-EDEF ==============================//



struct Pattern_buffer
{
  GPattern* pat;
  bool is_closed;
  double evaluation;

  Pattern_buffer()
  {
    pat = NULL, is_closed = 1, evaluation = 0.0;
  }

  Pattern_buffer(GPattern* p, bool closed, double ev)
  {
    pat = p, is_closed = closed, evaluation = ev;
  }

  bool operator < ( const Pattern_buffer& ob) const
  {
    if( evaluation != ob.evaluation )
      return evaluation < ob.evaluation;
    else if(is_closed != ob.is_closed)
      return is_closed < ob.is_closed;
    
    uint minSize = min( pat->v_Tokens.size(), ob.pat->v_Tokens.size() );
    uint i = 0;
    int cmp = 0;

    while ( i < minSize )
    {
      cmp = cmpGTokenCanonTest( pat->v_Tokens.at(i), ob.pat->v_Tokens.at(i) );
      if ( cmp > 0 )
        return true;
      else if ( cmp < 0 )
        return false;
      else
        i++;
    }
    if ( pat->v_Tokens.size() > ob.pat->v_Tokens.size() )
      return true;
    else
      return false;

  }

};

//=============================== VARIABLES ==================================//
//================================ METHODS ===================================//

//================================= CLASS ====================================//
/**
 * @brief The G class
 * Class that store vocabulary of pattern and allow to process grima mining
 * algorithm for each class.
 */
class Hologram
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

  unordered_map< string , searchtreenode* > nodes_pointers;

  unordered_map< string , bool > class_patterns;

  unordered_map< string , bool > all_patterns;
  
  set<Pattern_buffer> pattern_buffer_set;

  map<long long, searchtreenode*> search_tree_nodes;

  double delta;
  int N_delta;

  vector<int> class_count;

  int number_of_classes;
  int number_of_graphs;

  searchtreenode* last_father;

  int roll_depth;

  bool do_update;

  int current_class_id;

  searchtreenode* root;

  long long nodes_counter;
  // Public Structure & Typedef ______________________________________________//
  // Public Constructor/Desctructor __________________________________________//
  /**
   * @brief Grima
   * Default constructor
   * Initialize variables and create new GVocab object
   */
  Hologram();

  Hologram(float _minf,GClassDB* _pClassDB,vector<string> classes_names);
  /**
   * @brief ~G
   * Default Destructor
   * Free memory by deleting GVocab object and clearing vectors
   */
  ~Hologram();

  // Accessor ________________________________________________________________//
  // Mutator _________________________________________________________________//
  // Public Methods __________________________________________________________//
  /**
   * TODO
   */
  void set_class_id(int id);

  void initialize( );

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
  void saveData( );

  void unbuffer_class_patterns();

  inline double UCB(searchtreenode* cur, searchtreenode* child);

  searchtreenode* best_child(searchtreenode* cur,GToken& ext);

  searchtreenode* select(searchtreenode* cur);

  searchtreenode* expand(searchtreenode* cur,GToken& ext,GExtensionData& tmp);

  void add_parent(searchtreenode* cur,searchtreenode* parent,const GToken& lastExt);

  void build_pattern(searchtreenode* selcted_node, GPattern* pPattern);

  int roll_out( searchtreenode* cur, 
                searchtreenode* parent,
                bool rollout_first_level,
                const vector<GGraph*>& v_Graphs,
                const GGlobFreq       minFreq,    //Mininmum global frequency
                GPattern*        pPattern,
                const GToken&    lastExt, 
                GTokenData      &tokenData,
                GExtensionData  &suppData,   // Tmp variable, supposed frequency
                GExtensionData  &prevData);

  void update_ancestors(searchtreenode* cur, double delta); 

  double WRAcc(GTokenData& tokenData,int classID,int support); 

  void delete_tree_node(searchtreenode* cur);

  void clean();

  void delete_search_subtree(searchtreenode* cur);
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
#endif // HOLOGRAM_HPP
