#ifndef USERDEFINED_H
#define USERDEFINED_H

#include <map>
#include <string>
//#include <vector>

using namespace std;

class user_defined{
 private:

 public:
////////////////////
//  constructors  //
////////////////////
  user_defined();
  //  user_defined(string _process_class, string _subprocess, string _type_perturbative_order, string _type_contribution, string _type_correction, int _order_alpha_s, int _order_alpha_e, int _order_interference);

  vector<string> switch_name;
  vector<int> switch_value;
  map<string, int> switch_map;

  vector<string> cut_name;
  vector<double> cut_value;
  map<string, int> cut_map;

  vector<string> int_name;
  vector<int> int_value;
  map<string, int> int_map;

  vector<string> double_name;
  vector<double> double_value;
  map<string, int> double_map;

  vector<string> string_name;
  vector<string> string_value;
  map<string, int> string_map;

  // probably only particle name needed !!!

  vector<string> particle_name;
  vector<string> particle_value;
  map<string, int> particle_map;
};
#endif
