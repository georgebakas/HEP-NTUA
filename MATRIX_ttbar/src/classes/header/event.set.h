#ifndef EVENTSET_H
#define EVENTSET_H

#include <map>
#include <string>
#include <vector>
#include "logger.h"
#include "user.defined.h"

using namespace std;

class event_set{
 private:

 public:
////////////////////
//  constructors  //
////////////////////
  event_set();

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

  vector<string> name_fiducial_cut;
  vector<string> fiducial_cut_list;
  vector<fiducialcut> fiducial_cut;

  void define_basic_object_list();
  void define_specific_object_list(user_defined & user);
};
#endif
