#include "../include/classes.cxx"
////////////////////
//  constructors  //
////////////////////
dipole_set::dipole_set(){}
dipole_set::dipole_set(string name, vector<int> type_parton, vector<int> basic_type_parton, int type_dipole, int type_splitting, int no_R_emitter_1, int no_R_emitter_2, int no_R_spectator, int no_A_emitter, int no_A_spectator, int type_correction){
  static Logger logger("dipole_set::dipole_set");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  logger << LOG_DEBUG_VERBOSE << "initialization dipole_candidate" << endl;
  // initialization dipole_candidate

  _name = name;
  _type_parton = type_parton;
  _basic_type_parton = basic_type_parton;

  _type_dipole = type_dipole;
  _type_splitting = type_splitting;
  _no_R_emitter_1 = no_R_emitter_1;
  _no_R_emitter_2 = no_R_emitter_2;
  _no_R_spectator = no_R_spectator;
  _no_A_emitter = no_A_emitter;
  _no_A_spectator = no_A_spectator;
  _binary_R_emitter_1 = intpow(2, no_R_emitter_1 - 1);
  _binary_R_emitter_2 = intpow(2, no_R_emitter_2 - 1);
  _binary_R_spectator = intpow(2, no_R_spectator - 1);
  _binary_A_emitter = intpow(2, no_A_emitter - 1);
  _binary_A_spectator = intpow(2, no_A_spectator - 1);
  _type_correction = type_correction;

  _dx_pa.resize(type_parton.size(), vector<int> (1));
  for (int i_p = 0; i_p < type_parton.size(); i_p++){_dx_pa[i_p][0] = type_parton[i_p];}
  /*
  _phasespace = phasespace;
  _colourmatrix = colourmatrix;
  _spinorder = spinorder;
  _fckm = fckm;
  _data = data;
  _exp_pdf = exp_pdf;
  */

  if (_no_A_emitter < _no_A_spectator){no_BLHA_entry = ((_no_A_spectator - 1) * (_no_A_spectator - 2)) / 2 + (_no_A_emitter - 1);}
  else if (_no_A_emitter > _no_A_spectator){no_BLHA_entry = ((_no_A_emitter - 1) * (_no_A_emitter - 2)) / 2 + (_no_A_spectator - 1);}

  // added on 20181024
  if (_type_parton[_no_A_emitter] == 0){_colour_factor = C_A;}
  else if (abs(_type_parton[_no_A_emitter]) < 7){_colour_factor = C_F;}
  else {_colour_factor = 0.;}
  
  //  cout << "_no_R_emitter_2 = " << _no_R_emitter_2 << endl;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
dipole_set::dipole_set(const dipole_set & QEW_dipole_candidate, int no_map, vector<int> o_map, int no_prc, vector<int> o_prc, int n_channel, int sum_channel, double charge_factor, double symmetry_factor, int massive){
  static Logger logger("dipole_set::dipole_set");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  logger << LOG_DEBUG_VERBOSE << "initialization selected dipole" << endl;
   // initialization selected dipole

  _name = QEW_dipole_candidate._name;
  _type_parton = QEW_dipole_candidate._type_parton;
  _basic_type_parton = QEW_dipole_candidate._basic_type_parton;

  _charge_factor = charge_factor;
  _symmetry_factor = symmetry_factor;
  _no_map = no_map;
  _o_map = o_map;
  _no_prc = no_prc;
  _o_prc = o_prc;
  _n_channel = n_channel;
  _sum_channel = sum_channel;

  _type_dipole = QEW_dipole_candidate._type_dipole;
  _type_splitting = QEW_dipole_candidate._type_splitting;
  _no_R_emitter_1 = QEW_dipole_candidate._no_R_emitter_1;
  _no_R_emitter_2 = QEW_dipole_candidate._no_R_emitter_2;
  _no_R_spectator = QEW_dipole_candidate._no_R_spectator;
  _no_A_emitter = QEW_dipole_candidate._no_A_emitter;
  _no_A_spectator = QEW_dipole_candidate._no_A_spectator;
  _binary_R_emitter_1 = QEW_dipole_candidate._binary_R_emitter_1;
  _binary_R_emitter_2 = QEW_dipole_candidate._binary_R_emitter_2;
  _binary_R_spectator = QEW_dipole_candidate._binary_R_spectator;
  _binary_A_emitter = QEW_dipole_candidate._binary_A_emitter;
  _binary_A_spectator = QEW_dipole_candidate._binary_A_spectator;
  _type_correction = QEW_dipole_candidate._type_correction;
  _massive = massive;

  _dx_pa = QEW_dipole_candidate._dx_pa;
  /*
  _phasespace = phasespace;
  _colourmatrix = colourmatrix;
  _spinorder = spinorder;
  _fckm = fckm;
  _data = data;
  _exp_pdf = exp_pdf;
  */
  no_BLHA_entry = QEW_dipole_candidate.no_BLHA_entry;
  /*
  if (_no_A_emitter < _no_A_spectator){no_BLHA_entry = ((_no_A_spectator - 1) * (_no_A_spectator - 2)) / 2 + (_no_A_emitter - 1);}
  else if (_no_A_emitter > _no_A_spectator){no_BLHA_entry = ((_no_A_emitter - 1) * (_no_A_emitter - 2)) / 2 + (_no_A_spectator - 1);}
  */

  // added on 20181024
  if (_type_parton[_no_A_emitter] == 0){_colour_factor = C_A;}
  else if (abs(_type_parton[_no_A_emitter]) < 7){_colour_factor = C_F;}
  else {_colour_factor = 0.;}
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
 }
dipole_set::dipole_set(string name, vector<int> type_parton, vector<int> basic_type_parton, double symmetry_factor, int no_map, vector<int> o_map, int no_prc, vector<int> o_prc, int n_channel){
   // initialization real contribution
  _name = name;
  _type_parton = type_parton;
  _basic_type_parton = basic_type_parton;

  _charge_factor = 1.;
  _symmetry_factor = symmetry_factor;
  _no_map = no_map;
  _o_map = o_map;
  _no_prc = no_prc;
  _o_prc = o_prc;
  _n_channel = n_channel;
  _sum_channel = n_channel;

  _type_dipole = 0;
  _type_splitting = 0;
  _no_R_emitter_1 = 0;
  _no_R_emitter_2 = 0;
  _no_R_spectator = 0;
  _no_A_emitter = 0;
  _no_A_spectator = 0;
  _binary_R_emitter_1 = 0;
  _binary_R_emitter_2 = 0;
  _binary_R_spectator = 0;
  _binary_A_emitter = 0;
  _binary_A_spectator = 0;

  _type_correction = 0;
  _massive = 0;

  _dx_pa.resize(type_parton.size(), vector<int> (1));
  for (int i_p = 0; i_p < type_parton.size(); i_p++){_dx_pa[i_p][0] = type_parton[i_p];}
  /*
  _phasespace = phasespace;
  _colourmatrix = colourmatrix;
  _spinorder = spinorder;
  _fckm = fckm;
  _data = data;
  _exp_pdf = exp_pdf;
  */

  no_BLHA_entry = 0;

  // added on 20181024
  if (_type_parton[_no_A_emitter] == 0){_colour_factor = C_A;}
  else if (abs(_type_parton[_no_A_emitter]) < 7){_colour_factor = C_F;}
  else {_colour_factor = 0.;}
    
}
dipole_set::dipole_set(string name, double charge_factor, double symmetry_factor, vector<vector<int> > dx_pa, vector<int> phasespace, vector<int> o_prc, vector<int> o_map, vector<vector<double> > colourmatrix, vector<vector<int> > spinorder, vector<int> fckm, vector<int> data, int no_prc, int no_map, int n_channel, double exp_pdf, vector<int> type_parton, vector<int> basic_type_parton, int type_dipole, int type_splitting, int no_R_emitter_1, int no_R_emitter_2, int no_R_spectator, int no_A_emitter, int no_A_spectator, int type_correction, int massive){
  _name = name;
  _charge_factor = charge_factor;
  _symmetry_factor = symmetry_factor;
  _dx_pa = dx_pa;
  _phasespace = phasespace;
  _o_prc = o_prc;
  _o_map = o_map;
  _colourmatrix = colourmatrix;
  _spinorder = spinorder;
  _fckm = fckm;
  _data = data;
  _no_prc = no_prc;
  _no_map = no_map;
  _n_channel = n_channel;
  _exp_pdf = exp_pdf;
  _type_parton = type_parton;
  _basic_type_parton = basic_type_parton;
  _type_dipole = type_dipole;
  _type_splitting = type_splitting;
  _no_R_emitter_1 = no_R_emitter_1;
  _no_R_emitter_2 = no_R_emitter_2;
  _no_R_spectator = no_R_spectator;
  _no_A_emitter = no_A_emitter;
  _no_A_spectator = no_A_spectator;

  _binary_R_emitter_1 = intpow(2, no_R_emitter_1 - 1);
  _binary_R_emitter_2 = intpow(2, no_R_emitter_2 - 1);
  _binary_R_spectator = intpow(2, no_R_spectator - 1);
  _binary_A_emitter = intpow(2, no_A_emitter - 1);
  _binary_A_spectator = intpow(2, no_A_spectator - 1);

  _type_correction = type_correction;
  _massive = massive;

  if (_no_A_emitter < _no_A_spectator){no_BLHA_entry = ((_no_A_spectator - 1) * (_no_A_spectator - 2)) / 2 + (_no_A_emitter - 1);}
  else if (_no_A_emitter > _no_A_spectator){no_BLHA_entry = ((_no_A_emitter - 1) * (_no_A_emitter - 2)) / 2 + (_no_A_spectator - 1);}

  // added on 20181024
  if (_type_parton[_no_A_emitter] == 0){_colour_factor = C_A;}
  else if (abs(_type_parton[_no_A_emitter]) < 7){_colour_factor = C_F;}
  else {_colour_factor = 0.;}
  
}

///////////////////////
//  access elements  //
///////////////////////
string dipole_set::name() const {return _name;}
const vector<int> dipole_set::type_parton() const {return _type_parton;}
const vector<int> dipole_set::basic_type_parton() const {return _basic_type_parton;}
double dipole_set::charge_factor() const {return _charge_factor;}
double dipole_set::colour_factor() const {return _colour_factor;}
double dipole_set::symmetry_factor() const {return _symmetry_factor;}
int dipole_set::no_map() const {return _no_map;}
const vector<int> dipole_set::o_map() const {return _o_map;}
int dipole_set::no_prc() const {return _no_prc;}
const vector<int> dipole_set::o_prc() const {return _o_prc;}
int dipole_set::n_channel() const {return _n_channel;}
int dipole_set::sum_channel() const {return _sum_channel;}

int dipole_set::type_dipole() const {return _type_dipole;}
int dipole_set::type_splitting() const {return _type_splitting;}
int dipole_set::no_R_emitter_1() const {return _no_R_emitter_1;}
int dipole_set::no_R_emitter_2() const {return _no_R_emitter_2;}
int dipole_set::no_R_spectator() const {return _no_R_spectator;}
int dipole_set::no_A_emitter() const {return _no_A_emitter;}
int dipole_set::no_A_spectator() const {return _no_A_spectator;}
int dipole_set::binary_R_emitter_1() const {return _binary_R_emitter_1;}
int dipole_set::binary_R_emitter_2() const {return _binary_R_emitter_2;}
int dipole_set::binary_R_spectator() const {return _binary_R_spectator;}
int dipole_set::binary_A_emitter() const {return _binary_A_emitter;}
int dipole_set::binary_A_spectator() const {return _binary_A_spectator;}

int dipole_set::type_correction() const {return _type_correction;}
int dipole_set::massive() const {return _massive;}


const vector<vector<int> > dipole_set::dx_pa() const {return _dx_pa;}
const vector<vector<double> > dipole_set::colourmatrix() const {return _colourmatrix;}
const vector<vector<int> > dipole_set::spinorder() const {return _spinorder;}
const vector<int> dipole_set::fckm() const {return _fckm;}
const vector<int> dipole_set::data() const {return _data;}
double  dipole_set::exp_pdf() const {return _exp_pdf;}
const vector<int> dipole_set::phasespace() const {return _phasespace;}

///////////////
//  methods  //
///////////////
