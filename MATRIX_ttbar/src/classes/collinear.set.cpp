#include "../include/classes.cxx"
////////////////////
//  constructors  //
////////////////////
collinear_set::collinear_set(){}
collinear_set::collinear_set(string name, vector<string> all_name, int type, vector<int> in_collinear, int no_prc, vector<int> type_parton, vector<vector<int> > all_pdf, double charge_factor, int no_emitter, int no_spectator, vector<int> pair, int type_correction, int massive){
  //vector<vector<int> > all_pdf, 
  _name = name;
  _all_name = all_name;
  _type = type;
  _pair = pair;
  _in_collinear = in_collinear;
  _no_prc = no_prc;
  _type_parton = type_parton;
  _all_pdf = all_pdf;
  _charge_factor = charge_factor;
  _no_emitter = no_emitter;
  _no_spectator = no_spectator;
  if (no_spectator == 0){_no_OL_entry = 0;}
  else if (no_emitter < no_spectator){_no_OL_entry = (no_spectator * (no_spectator - 1)) / 2 + no_emitter;}
  else if (no_emitter > no_spectator){_no_OL_entry = (no_emitter * (no_emitter - 1)) / 2 + no_spectator;}
  else {cout << "Should not happen !" << endl;}

  _type_correction = type_correction;
  _massive = massive;

  // added on 20141111
  if (_no_spectator == 0){no_BLHA_entry = 0;}
  else if (_no_emitter < _no_spectator){no_BLHA_entry = ((_no_spectator - 1) * (_no_spectator - 2)) / 2 + (_no_emitter - 1);}
  else if (_no_emitter > _no_spectator){no_BLHA_entry = ((_no_emitter - 1) * (_no_emitter - 2)) / 2 + (_no_spectator - 1);}

  // added on 20181017
  if (type_parton[no_emitter] == 0){_colour_factor = C_A;}
  else if (abs(type_parton[no_emitter]) < 7){_colour_factor = C_F;}
  else {_colour_factor = 0.;}
  
}


///////////////////////
//  access elements  //
///////////////////////
string collinear_set::name() const {return _name;}
const vector<string> collinear_set::all_name() const {return _all_name;}
int collinear_set::type() const {return _type;}
const vector<int> collinear_set::in_collinear() const {return _in_collinear;}
int collinear_set::no_prc() const {return _no_prc;}
const vector<int> collinear_set::type_parton() const {return _type_parton;}
//const vector<int> collinear_set::all_pdf() const {return _all_pdf;}
const vector<vector<int> > collinear_set::all_pdf() const {return _all_pdf;}
double collinear_set::charge_factor() const {return _charge_factor;}
double collinear_set::colour_factor() const {return _colour_factor;}
int collinear_set::no_emitter() const {return _no_emitter;}
int collinear_set::no_spectator() const {return _no_spectator;}
int collinear_set::no_OL_entry() const {return _no_OL_entry;}
const vector<int> collinear_set::pair() const {return _pair;}
int collinear_set::type_correction() const {return _type_correction;}
int collinear_set::massive() const {return _massive;}
//const vector<int> collinear_set::type_parton() const {return _type_parton;}
/*
double collinear_set::symmetry_factor() const {return _symmetry_factor;}
int collinear_set::no_map() const {return _no_map;}
const vector<int> collinear_set::o_map() const {return _o_map;}
int collinear_set::no_prc() const {return _no_prc;}
const vector<int> collinear_set::o_prc() const {return _o_prc;}
int collinear_set::n_channel() const {return _n_channel;}
int collinear_set::sum_channel() const {return _sum_channel;}
*/
/*
int collinear_set::type_splitting() const {return _type_splitting;}
*/
//int collinear_set::binary_A_emitter() const {return _binary_A_emitter;}
//int collinear_set::binary_A_spectator() const {return _binary_A_spectator;}
/*
const vector<vector<int> > collinear_set::dx_pa() const {return _dx_pa;}
const vector<vector<double> > collinear_set::colourmatrix() const {return _colourmatrix;}
const vector<vector<int> > collinear_set::spinorder() const {return _spinorder;}
const vector<int> collinear_set::fckm() const {return _fckm;}
const vector<int> collinear_set::data() const {return _data;}
double  collinear_set::exp_pdf() const {return _exp_pdf;}
const vector<int> collinear_set::phasespace() const {return _phasespace;}
*/
///////////////
//  methods  //
///////////////
