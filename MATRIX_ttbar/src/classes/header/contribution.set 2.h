#ifndef CONTRIBUTIONSET_H
#define CONTRIBUTIONSET_H

#include <map>
#include <string>
#include <vector>

using namespace std;

class contribution_set{
 private:

 public:
////////////////////
//  constructors  //
////////////////////
  contribution_set();
  contribution_set(string _process_class, vector<string> _decay, string _type_perturbative_order, string _type_contribution, string _type_correction, int _order_alpha_s, int _order_alpha_e, int _order_interference, string _subprocess);

  void determination_order_alpha_s_born();
  void determination_class_contribution();
  
  void readin_hadronic_process();
  void readin_subprocess();
  //  void determine_basic_subprocess();
  //  void determine_subprocess();

  void readin_subprocess_from_name(string & this_processname, vector<int> & this_type_parton, int temp_switch);
  void fill_code_particle();
  void fill_name_particle();

  void newscheme_readin_subprocess_from_name(string & this_processname, vector<int> & this_type_parton, int temp_switch);
  void newscheme_fill_code_particle();
  void newscheme_fill_name_particle();

  void oldscheme_readin_subprocess_from_name(string & this_processname, vector<int> & this_type_parton, int temp_switch);
  void oldscheme_fill_code_particle();
  void oldscheme_fill_name_particle();

  string path_MUNICH;
  
  int switch_to_newscheme;

  string process_class;
  string basic_process_class;

  vector<string> decay;

  string type_perturbative_order;
  string type_contribution;
  string type_correction;
  int contribution_order_alpha_s;
  int contribution_order_alpha_e;
  int contribution_order_interference;

 
  string subprocess;
  string basic_subprocess; // probably not needed...

  vector<vector<int> > type_parton;
  vector<vector<int> > basic_type_parton;

  // new: to replace no_prc/no_map and o_prc/o_map in observable_set/phasespace_set
  vector<int> no_process_parton; // no_prc
  vector<vector<int> > swap_parton; // o_prc

  vector<int> order_contribution_alpha_s;
  vector<int> order_contribution_alpha_e;
  vector<int> order_contribution_interference;
  vector<int> order_phasespace_alpha_s;
  vector<int> order_phasespace_alpha_e;
  vector<int> order_phasespace_interference;

  int n_particle;
  int process_type;
  int n_particle_born;
  int n_jet_born;
  int n_photon_born;
  int order_alpha_s_born;
  
  map<string,int> code_particle;
  map<int,string> name_particle;

  map<int, int> lepton_exchange;

  int class_contribution_born;
  int class_contribution_CS_collinear;
  int class_contribution_CS_virtual;
  int class_contribution_CS_real;
  int class_contribution_QT_collinear;
  int class_contribution_QT_virtual;
  int class_contribution_QT_real;
  int class_contribution_NJ_collinear;
  int class_contribution_NJ_virtual;
  int class_contribution_NJ_real;
  int class_contribution_collinear;
  int class_contribution_virtual;
  int class_contribution_real;
  int class_contribution_loopinduced;
  int class_contribution_IRcut_implicit;
  int class_contribution_IRcut_explicit;
  int class_contribution_CS;
  int class_contribution_QT;
  int class_contribution_NJ;
  int class_contribution_QT_CS;
  int class_contribution_NJ_CS;
  
  //  string run_mode;
};
#endif
