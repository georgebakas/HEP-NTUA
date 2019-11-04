#include "../include/classes.cxx"

////////////////////
//  constructors  //
////////////////////
user_defined::user_defined(){
  switch_name.resize(1);
  switch_value.resize(1);
  //  switch_map;

  cut_name.resize(1);
  cut_value.resize(1);
  //  cut_map;

  int_name.resize(1);
  int_value.resize(1);
  //  int_map;

  double_name.resize(1);
  double_value.resize(1);
  //  double_map;

  string_name.resize(1);
  string_value.resize(1);
  //  string_map;

  //  particle_name.resize(1);
  //  particle_value.resize(1);
}
  /*
user_defined::user_defined(string _process_class, string _subprocess, string _type_perturbative_order, string _type_contribution, string _type_correction, int _order_alpha_s, int _order_alpha_e, int _order_interference){
  Logger logger("user_defined::user_defined");
  logger << LOG_DEBUG << "started" << endl;
  logger << LOG_DEBUG << "The same information is also contained in observable_set !!!" << endl;

  process_class = _process_class;
  subprocess = _subprocess;

  type_perturbative_order = _type_perturbative_order;
  type_contribution = _type_contribution;
  type_correction = _type_correction;

  contribution_order_alpha_s = _order_alpha_s;
  contribution_order_alpha_e = _order_alpha_e;
  contribution_order_interference = _order_interference;

  logger << LOG_DEBUG << "before contribution output" << endl;

  logger << LOG_DEBUG << "process_class = " << process_class << endl;
  logger << LOG_DEBUG << "subprocess = " << subprocess << endl;
  logger << LOG_DEBUG << "type_perturbative_order = " << type_perturbative_order << endl;
  logger << LOG_DEBUG << "type_contribution = " << type_contribution << endl;
  logger << LOG_DEBUG << "type_correction = " << type_correction << endl;
  
  logger << LOG_DEBUG << "contribution_order_alpha_s = " << contribution_order_alpha_s << endl;
  logger << LOG_DEBUG << "contribution_order_alpha_e = " << contribution_order_alpha_e << endl;
  logger << LOG_DEBUG << "contribution_order_interference = " << contribution_order_interference << endl;
  
  logger << LOG_DEBUG << "after contribution output" << endl;
  logger << LOG_DEBUG << "finished" << endl;
}
*/

