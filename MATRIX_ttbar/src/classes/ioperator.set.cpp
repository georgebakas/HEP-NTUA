#include "../include/classes.cxx"
////////////////////
//  constructors  //
////////////////////
ioperator_set::ioperator_set(){}
ioperator_set::ioperator_set(string name, int type, vector<int> pair, int no_prc, vector<int> type_parton, double charge_factor, int no_emitter, int no_spectator, int type_correction, int massive){
  _name = name;
  _type = type;
  _pair = pair;
  _no_prc = no_prc;
  _type_parton = type_parton;
  _charge_factor = charge_factor;
  _no_emitter = no_emitter;
  _no_spectator = no_spectator;
  if (no_spectator == 0){_no_OL_entry = 0;}
  else if (no_emitter < no_spectator){_no_OL_entry = (no_spectator * (no_spectator - 1)) / 2 + no_emitter;}
  else if (no_emitter > no_spectator){_no_OL_entry = (no_emitter * (no_emitter - 1)) / 2 + no_spectator;}
  else {cout << "Should not happen !" << endl;}

  _type_correction = type_correction;
  _massive = massive;

  if (_no_emitter < _no_spectator){no_BLHA_entry = ((_no_spectator - 1) * (_no_spectator - 2)) / 2 + (_no_emitter - 1);}
  else if (_no_emitter > _no_spectator){no_BLHA_entry = ((_no_emitter - 1) * (_no_emitter - 2)) / 2 + (_no_spectator - 1);}

  // added on 20181022
  if (type_parton[no_emitter] == 0){_colour_factor = C_A;}
  else if (abs(type_parton[no_emitter]) < 7){_colour_factor = C_F;}
  else {_colour_factor = 0.;}

}

///////////////////////
//  access elements  //
///////////////////////
string ioperator_set::name() const {return _name;}
int ioperator_set::type() const {return _type;}
const vector<int> ioperator_set::pair() const {return _pair;}
int ioperator_set::no_prc() const {return _no_prc;}
const vector<int> ioperator_set::type_parton() const {return _type_parton;}
double ioperator_set::charge_factor() const {return _charge_factor;}
double ioperator_set::colour_factor() const {return _colour_factor;}
int ioperator_set::no_emitter() const {return _no_emitter;}
int ioperator_set::no_spectator() const {return _no_spectator;}
int ioperator_set::no_OL_entry() const {return _no_OL_entry;}
int ioperator_set::type_correction() const {return _type_correction;}
int ioperator_set::massive() const {return _massive;}

//const vector<int> ioperator_set::type_parton() const {return _type_parton;}
/*
double ioperator_set::symmetry_factor() const {return _symmetry_factor;}
int ioperator_set::no_map() const {return _no_map;}
const vector<int> ioperator_set::o_map() const {return _o_map;}
int ioperator_set::no_prc() const {return _no_prc;}
const vector<int> ioperator_set::o_prc() const {return _o_prc;}
int ioperator_set::n_channel() const {return _n_channel;}
int ioperator_set::sum_channel() const {return _sum_channel;}
*/
/*
int ioperator_set::type_splitting() const {return _type_splitting;}
*/
//int ioperator_set::binary_A_emitter() const {return _binary_A_emitter;}
//int ioperator_set::binary_A_spectator() const {return _binary_A_spectator;}
/*
const vector<vector<int> > ioperator_set::dx_pa() const {return _dx_pa;}
const vector<vector<double> > ioperator_set::colourmatrix() const {return _colourmatrix;}
const vector<vector<int> > ioperator_set::spinorder() const {return _spinorder;}
const vector<int> ioperator_set::fckm() const {return _fckm;}
const vector<int> ioperator_set::data() const {return _data;}
double  ioperator_set::exp_pdf() const {return _exp_pdf;}
const vector<int> ioperator_set::phasespace() const {return _phasespace;}
*/
///////////////
//  methods  //
///////////////
