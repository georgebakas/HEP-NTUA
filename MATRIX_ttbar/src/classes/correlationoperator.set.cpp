#include "../include/classes.cxx"
////////////////////
//  constructors  //
////////////////////
correlationoperator::correlationoperator(){}
correlationoperator::correlationoperator(string _name, int _type, vector<int> _pair, int _no_prc, vector<int> _type_parton, double _charge_factor, int _no_emitter, int _no_spectator, int _type_correction, int _type_combination, int _massive){
  name = _name;
  type = _type;
  pair = _pair;
  no_prc = _no_prc;
  type_parton = _type_parton;
  charge_factor = _charge_factor;
  no_emitter = _no_emitter;
  no_spectator = _no_spectator;
  type_combination = _type_combination;
  type_correction = _type_correction;
  massive = _massive;

  if (no_spectator == 0){no_OL_entry = 0;}
  else if (no_emitter < no_spectator){no_OL_entry = (no_spectator * (no_spectator - 1)) / 2 + no_emitter;}
  else if (no_emitter > no_spectator){no_OL_entry = (no_emitter * (no_emitter - 1)) / 2 + no_spectator;}
  else {no_OL_entry = 0;}

  if (no_emitter < no_spectator){no_BLHA_entry = ((no_spectator - 1) * (no_spectator - 2)) / 2 + (no_emitter - 1);}
  else if (no_emitter > no_spectator){no_BLHA_entry = ((no_emitter - 1) * (no_emitter - 2)) / 2 + (no_spectator - 1);}
  else {no_BLHA_entry = 0;}

  // added on 20181128
  if (_type_parton[_no_emitter] == 0){colour_factor = C_A;}
  else if (abs(_type_parton[_no_emitter]) < 7){colour_factor = C_F;}
  else {colour_factor = 0.;}
  
}

///////////////////////
//  access elements  //
///////////////////////

///////////////
//  methods  //
///////////////
