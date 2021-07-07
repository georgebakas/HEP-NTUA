#ifndef FIDUCIALCUT_H
#define FIDUCIALCUT_H

#include <map>
#include <string>
#include <vector>
#include "logger.h"
#include "user.defined.h"

using namespace std;

class event_set;
class observable_set;

typedef enum {
    UPPER = 0,
    LOWER,
    WINDOW,
    GAP,
} cut_type;

class fiducialcut{
 private:

 public:
////////////////////
//  constructors  //
////////////////////

  fiducialcut();
  fiducialcut(string & _name, event_set & _esi, observable_set & _osi);

////////////////
//  elements  //
////////////////

  observable_set *osi;
  event_set *esi;

  string name;

  string observable;
  
  cut_type type;

  vector<string> type_particle;
  vector<vector<int> > no_particles;

  vector<int> order_min;
  vector<int> order_max;

  double cut_value_lower;
  double cut_value_upper;

  int n_combination;
  vector<int> cut_combination;
  int cut;
  
///////////////
//  methods  //
///////////////

  void apply_fiducialcut(int i_a);


  /*
  void determine_n_partonlevel(vector<vector<int> > & type_parton);

  map<string, int> observed_object;
  vector<string> object_list;
  vector<int> object_category;

  map<string, int> observed_object_selection;
  vector<string> object_list_selection;
  vector<int> object_category_selection;

  vector<define_particle_set> pda;  // particle definition all
  vector<define_particle_set> pds;  // particle definition selection

  int n_parton_nu;

  void define_basic_object_list();
  void define_specific_object_list(user_defined & user);
  */

  friend std::ostream & operator << ( std::ostream &, const fiducialcut &);
};
#endif
