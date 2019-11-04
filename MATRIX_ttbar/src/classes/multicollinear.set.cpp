#include "../include/classes.cxx"
////////////////////
//  constructors  //
////////////////////
multicollinear_set::multicollinear_set(){}

//multicollinear_set::multicollinear_set(multicollinear_set ms){}

multicollinear_set::multicollinear_set(int _no_prc, vector<int> _type_parton, vector<vector<int> > _pdf){

  no_prc = _no_prc;
  type_parton = _type_parton;
  pdf = _pdf;

  name = "C";
  all_name = vector<string> (pdf.size(), name);
  type_splitting.resize(1, 0);
  charge_factor.resize(1, 1.);
  in_collinear.resize(1, vector<int> (3, 0));
  no_emitter.resize(1, 0);
  type_correction.resize(1, 0);
  massive.resize(1, 0);
  pair.resize(1, vector<int> (2, 0));
  no_spectator.resize(1, 0);
  no_OL_entry.resize(1, 0);
  no_BLHA_entry.resize(1, 0);

  n_emission.resize(3, 0);
  emission.resize(3, 0);

  previous_splitting.resize(1, 0);
  no_endpoint.resize(4, -1);
  type_splitting_full.resize(3, 0);
  type_splitting_leg.resize(3);
  n_emission_leg_type.resize(3, vector<int> (3, 0));
  //  order_emitter.resize(2, 0);
}

multicollinear_set::multicollinear_set(string _name, vector<string> _all_name, vector<int> _type_splitting, vector<vector<int> > _in_collinear, int _no_prc, vector<int> _type_parton, vector<vector<int> > _pdf, vector<double> _charge_factor, vector<int> _no_emitter, vector<int> _no_spectator, vector<vector<int> > _pair, vector<int> _type_correction, vector<int> _massive){
  //vector<vector<int> > pdf, 
  no_prc = _no_prc;
  type_parton = _type_parton;
  pdf = _pdf;

  name = _name;
  all_name = _all_name;
  type_splitting = _type_splitting;
  charge_factor = _charge_factor;
  in_collinear = _in_collinear;
  no_emitter = _no_emitter;
  type_correction = _type_correction;
  massive = _massive;

  // maybe modify this part: spectators are method-specific !!!
  pair = _pair;
  // needs to be modified - e.g. loop over spectators
  no_spectator = _no_spectator;

  for (int i_c = 0; i_c < no_spectator.size(); i_c++){
    if (no_spectator[i_c] == 0){no_OL_entry[i_c] = 0;}
    else if (no_emitter[i_c] < no_spectator[i_c]){no_OL_entry[i_c] = (no_spectator[i_c] * (no_spectator[i_c] - 1)) / 2 + no_emitter[i_c];}
    else if (no_emitter[i_c] > no_spectator[i_c]){no_OL_entry[i_c] = (no_emitter[i_c] * (no_emitter[i_c] - 1)) / 2 + no_spectator[i_c];}
    else {cout << "Should not happen !" << endl;}
    
    
    // added on 20141111
    if (no_spectator[i_c] == 0){no_BLHA_entry[i_c] = 0;}
    else if (no_emitter[i_c] < no_spectator[i_c]){no_BLHA_entry[i_c] = ((no_spectator[i_c] - 1) * (no_spectator[i_c] - 2)) / 2 + (no_emitter[i_c] - 1);}
    else if (no_emitter[i_c] > no_spectator[i_c]){no_BLHA_entry[i_c] = ((no_emitter[i_c] - 1) * (no_emitter[i_c] - 2)) / 2 + (no_spectator[i_c] - 1);}
  }
}

///////////////////////
//  access elements  //
///////////////////////
/*
string multicollinear_set::name() const {return _name;}
const vector<string> multicollinear_set::all_name() const {return _all_name;}
int multicollinear_set::type_splitting() const {return _type_splitting;}
const vector<int> multicollinear_set::in_collinear() const {return _in_collinear;}
int multicollinear_set::no_prc() const {return _no_prc;}
const vector<int> multicollinear_set::type_parton() const {return _type_parton;}
//const vector<int> multicollinear_set::pdf() const {return _pdf;}
const vector<vector<int> > multicollinear_set::pdf() const {return _pdf;}
double multicollinear_set::charge_factor() const {return _charge_factor;}
int multicollinear_set::no_emitter() const {return _no_emitter;}
int multicollinear_set::no_spectator() const {return _no_spectator;}
int multicollinear_set::no_OL_entry() const {return _no_OL_entry;}
const vector<int> multicollinear_set::pair() const {return _pair;}
int multicollinear_set::type_correction() const {return _type_correction;}
int multicollinear_set::massive() const {return _massive;}
*/
///////////////
//  methods  //
///////////////
