#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

////////////////////
//  constructors  //
////////////////////
inputparameter_set::inputparameter_set(){
  Logger logger("inputparameter_set::inputparameter_set");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  
  path_MUNICH = get_path();
  csi.path_MUNICH = path_MUNICH;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
inputparameter_set::inputparameter_set(string basic_process_class, string & subprocess){//, vector<string> & readin){
  Logger logger("inputparameter_set::inputparameter_set");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  
  path_MUNICH = get_path();
  csi.path_MUNICH = path_MUNICH;
  
  vector<string> readin;
  parameter_readin(subprocess, readin);

  /*
  int isystem = 0; 
  string xorder = "";
  */
  csi.basic_process_class = basic_process_class;
  csi.subprocess = subprocess;

  //  observable_set oset;
  //  user_defined user;

  //////////////////////////
  //  unknown parameters  //
  //////////////////////////

  path_to_main = "../../../../";
  // ???

  process_number = 0;
  // ???

  n_particle = 0;
  // ???

  test_output = 0;
  // ???

  tau_0 = 0.;
  // ???

  ckm_choice = 0;
  // not relevant by now - but maybe later...

  select_contribution = 0;
  // ???

  n_moments = 0;
  // ??? should be replaced by input file !!!

  dir_directory = "";
  // ???

  //  vector<steering_optimization> steering;
  // !!! removed

  //  n_madgraph = 0;
  // removed !!!

  no_contribution = 0;
  // ??? meaning ???

  //  filename_process = "";
  // !!! removed

  impulsname.resize(0);
  // ??? remove



  //////////////////////////////////////////////////////////////
  //  selection of process and contribution to be calculated  //
  //////////////////////////////////////////////////////////////

  /*
  csi.process_class = "";
  csi.decay.resize(0);
  csi.type_perturbative_order = "";
  csi.type_contribution = "";
  csi.type_correction = "";
  csi.contribution_order_alpha_s = 0;
  csi.contribution_order_alpha_e = 0;
  csi.contribution_order_interference = 0;
  */
  // default for empty csi initilization

  ////////////////////////////////////////////////
  //  switches to steer calculation of results  //
  ////////////////////////////////////////////////

  switch_result = 1;
  switch_distribution = 1;
  switch_moment = 0;

  ///////////////////////////////////////////////
  //  switches to steer output of calculation  //
  ///////////////////////////////////////////////

  switch_output_execution = 1;
  switch_output_integration = 1;
  switch_output_maxevent = 1;
  switch_output_comparison = 1;
  switch_output_gnuplot = 1;
  switch_output_proceeding = 1;
  switch_output_weights = 1;
  switch_output_result = 1;
  switch_output_moment = 0;
  switch_output_distribution = 1;

  switch_output_testpoint = 0;
  switch_output_cancellation_check = 0;
  switch_output_cutinfo = 0;

  switch_testcut = 0; // ???

  switch_console_output_runtime = 2;
  switch_console_output_tau_0 = 1;
  switch_console_output_techcut_RA = 1;
  switch_console_output_phasespace_issue = 1;
  switch_console_output_ME2_issue = 1;

 

  ////////////////////////////////////////////////////////////////////////
  //  unit of calculation output / result output / distribution output  //
  ////////////////////////////////////////////////////////////////////////

  unit_calculation = "fb";
  unit_result = "fb";
  unit_distribution = "fb";

  ///////////////////////
  //  beam parameters  //
  ///////////////////////

  E = 0.;
  coll_choice = 0;
  //  pdf_set = 0;
  // ??? only via LHAPDF ???

  //  pdf_content_modify = -1;
  // ???

  pdf_selection.resize(0);
  pdf_disable.resize(0);
  LHAPDFname = "";
  LHAPDFsubset = 0;

  ////////////////////////////////////////
  //  jet and jet-algorithm parameters  //
  ////////////////////////////////////////

  jet_algorithm = 0;
  jet_R_definition = 0;
  jet_R = 0.;
  parton_y_max = 1.e99;
  parton_eta_max = 1.e99;
  jet_algorithm_selection.resize(0);
  jet_algorithm_disable.resize(0);
  N_f = 0;
  N_f_active = 0;
  N_quarks = 6;
  N_nondecoupled = 0; 
  // !!! should be calculated from chosen PDF set


  ///////////////////////////////////////
  //  photon-recombination parameters  //
  ///////////////////////////////////////

  photon_recombination = 0;
  photon_R_definition = 0;
  photon_R = 0.;
  photon_E_threshold_ratio = 0.;
  photon_jet_algorithm = 0;
  photon_recombination_selection.resize(0);
  photon_recombination_disable.resize(0);

  photon_photon_recombination = 0;
  photon_photon_recombination_R = 0.;

  ///////////////////////////////////
  //  photon-isolation parameters  //
  ///////////////////////////////////

  frixione_isolation = 0;
  frixione_n = 0.;
  frixione_epsilon = 0.;
  frixione_fixed_ET_max = 0.;
  frixione_delta_0 = 0.;
  frixione_jet_removal = 0;

  /////////////////////////////////
  //  qT-subtraction parameters  //
  /////////////////////////////////

  switch_qTcut = 0;
  n_qTcut = 0;
  min_qTcut = 0.;
  step_qTcut = 0.;
  max_qTcut = 0.;
  binning_qTcut = "linear";
  selection_qTcut = "";

  selection_qTcut_distribution = "";
  selection_no_qTcut_distribution = "";
  selection_qTcut_result = "";
  selection_no_qTcut_result = "";
  selection_qTcut_integration = "";
  selection_no_qTcut_integration = "";
  //  value_qTcut;
  // !!! no input parameter !!!
  //  no_qTcut_distribution;
  // !!! no input parameter !!!
  //  value_qTcut_distribution;
  // !!! no input parameter !!!

  ////////////////////////////////////////////////////
  //  technical switches for selected contributions //
  ////////////////////////////////////////////////////

  switch_KP = 0;
  switch_VI = 0;
  switch_VI_bosonic_fermionic = 0;
  switch_polenorm = 0;
  // ??? default value ???

  switch_H1gg = 0;
  switch_H2 = 0;
 
  switch_old_qT_version = 1;

  switch_RS = 0;
  switch_RS_mapping = 0;
  switch_CM = 0;
  // ???

  switch_OL = 1;
  // ???


  switch_yuk = 0;
  order_y = 0;  // yukawa coupling squared

  switch_NJcut = 0;
  switch_NJcut_axes = 0;
  switch_NJcut_axes_energy = 0;
  switch_NJcut_measure = 0;



  ///////////////////////////////
  //  scale-choice parameters  //
  ///////////////////////////////

  scale_fact = 0.;
  scale_ren = 0.;
  dynamic_scale = 0;
  prefactor_reference = 1.;

  ///////////////////////////////////////
  //  scale variation parameters - CV  //
  ///////////////////////////////////////

  switch_CV = 0;
  n_scales_CV = 0;
  dynamic_scale_CV = 0;
  variation_mu_ren_CV = 0;
  // ???

  variation_mu_fact_CV = 0;
  // ???

  variation_factor_CV = 0;
  central_scale_CV = 0;
  // ???

  prefactor_CV = 1.0;


  //////////////////////////////////////
  //  weight-optimization parameters  //
  //////////////////////////////////////

  switch_n_events_opt = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - multi-channel Monte Carlo for final-state variables  //
  ////////////////////////////////////////////////////////////////////////////////////////////

  switch_MC = 0;
  n_alpha_steps = 0;
  n_alpha_events = 0;
  n_alpha_epc = 0;
  MCweight_min = 0;
  MCweight_limit_min = 0.;
  MCweight_limit_max = 0.;

  switch_use_alpha_after_IS = 0;
  switch_step_mode_grid = 0;

  MCweight_in_contribution = "";
  // ??? needed ??? maybe location of grid to be read in...

  MCweight_in_directory = "";
  // ??? needed ???
  
  /////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - IS for final-state variables  //
  /////////////////////////////////////////////////////////////////////

  switch_IS_MC = 0;
  switch_IS_mode_phasespace = 0;
  n_IS_events = 100000;
  n_IS_events_factor = 10;
  n_IS_steps = 20;
  n_IS_gridsize = 64;
  n_IS_gridsize_p = n_IS_gridsize;
  n_IS_gridsize_f = n_IS_gridsize;
  n_IS_gridsize_t_t = n_IS_gridsize;
  n_IS_gridsize_t_phi = n_IS_gridsize;
  n_IS_gridsize_d_cth = n_IS_gridsize;
  n_IS_gridsize_d_phi = n_IS_gridsize;

  /////////////////////////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - multi-channel Monte Carlo for CMS-energy mapping  //
  /////////////////////////////////////////////////////////////////////////////////////////

  switch_MC_tau = 0;

  //////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - IS for CMS-energy mapping  //
  //////////////////////////////////////////////////////////////////

  switch_IS_tau = 0;
  n_tau_steps = 0;
  n_tau_events = 0;
  n_tau_bins = 0;

  ////////////////////////////////////////////////////////////
  //  weight-optimization parameters - IS for x1x2 mapping  //
  ////////////////////////////////////////////////////////////

  switch_IS_x1x2 = 0;
  n_x1x2_steps = 0;
  n_x1x2_events = 0;
  n_x1x2_bins = 0;

  ///////////////////////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - IS for z1 and z2 mappings (collinear emission)  //
  ///////////////////////////////////////////////////////////////////////////////////////

  switch_IS_z1z2 = 0;
  n_z1z2_steps = 0;
  n_z1z2_events = 0;
  n_z1z2_bins = 0;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - multi-channel Monte Carlo for x-values in dipole mappings  //
  //////////////////////////////////////////////////////////////////////////////////////////////////

  switch_MC_x_dipole = 0;

  ////////////////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - IS for xy/zuv values in dipole mappings  //
  ////////////////////////////////////////////////////////////////////////////////

  n_IS_gridsize_xy = 1;
  n_IS_gridsize_zuv = 1;





  

  /////////////////////////////////////////
  //  phase-space generation parameters  //
  /////////////////////////////////////////

  nu = 0.;
  // removed ???

  nuxs = 0.;
  nuxt = 0.;
  exp_pdf = 0.;
  exp_pT = 0.;
  // removed ???

  exp_y = 0.;
  // removed ???

  exp_ij_k_y = 0.;
  exp_ij_k_z = 0.;
  exp_ij_a_x = 0.;
  exp_ij_a_z = 0.;
  exp_ai_k_x = 0.;
  exp_ai_k_u = 0.;
  exp_ai_b_x = 0.;
  exp_ai_b_v = 0.;

  ////////////////////////////////////////
  //  technical integration parameters  //
  ////////////////////////////////////////

  mass0 = 0.;
  map_technical_s = 0.;
  map_technical_t = 0.;
  map_technical_x = 0.;
  cut_technical = 0.;

  //////////////////////////////
  //  integration parameters  //
  //////////////////////////////

  n_events_max = 0;
  n_events_min = 0;
  sigma_normalization = 0.;
  sigma_normalization_deviation = 0.;
  n_step = 0;
  zwahl = 0;

  ////////////////////////////////////////////////////////////////////////////////////////
  //  parameters directly forwarded to OpenLoops - after the default settings are done  //
  ////////////////////////////////////////////////////////////////////////////////////////

  OL_parameter.resize(0);
  OL_value.resize(0);



  vector<string> user_variable;
  vector<string> user_variable_additional;
  vector<string> user_value;

  get_userinput_from_readin(user_variable, user_variable_additional, user_value, readin);

  ///////////////////////////////////////////////////////////////////////////
  //  some input that needs to be known before interpreting the remainder  //
  ///////////////////////////////////////////////////////////////////////////

  // determination of  'max_perturbative_QCD_order' from the input files
  // determination of  the chosen threshold for 'logger' output from the input files
  //  // determination of  switch_TSV


  max_perturbative_QCD_order = 0;
  for (int i = 0; i < user_variable.size(); i++){
    if (user_variable[i] == "type_perturbative_order"){
      //    if (user_variable[i] == "perturbative_QCD_order"){
      if      (user_value[i] == "LO"){if (max_perturbative_QCD_order < 0){max_perturbative_QCD_order = 0;}}
      else if (user_value[i] == "NLO"){if (max_perturbative_QCD_order < 1){max_perturbative_QCD_order = 1;}}
      else if (user_value[i] == "NNLO"){if (max_perturbative_QCD_order < 2){max_perturbative_QCD_order = 2;}}
    }
    else if (user_variable[i] == "output_level"){
      output_level = user_value[i].c_str();
      if (output_level == "DEBUG_VERBOSE"){Log::setLogThreshold(LOG_DEBUG_VERBOSE);}
      if (output_level == "DEBUG_POINT"){Log::setLogThreshold(LOG_DEBUG_POINT);}
      if (output_level == "DEBUG"){Log::setLogThreshold(LOG_DEBUG);}
      if (output_level == "INFO"){Log::setLogThreshold(LOG_INFO);}
      if (output_level == "WARN"){Log::setLogThreshold(LOG_WARN);}
      if (output_level == "ERROR"){Log::setLogThreshold(LOG_ERROR);}
      if (output_level == "FATAL"){Log::setLogThreshold(LOG_FATAL);}
      logger << LOG_DEBUG_VERBOSE << "output_level = " << output_level << endl;
    }
    //    else if (user_variable[i] == "switch_TSV"){switch_TSV = atoi(user_value[i].c_str());}
  }


  ////////////////////////////////////////////////////////////////////////////
  //  determination of particular type_perturbative_order to be calculated  //
  ////////////////////////////////////////////////////////////////////////////

  present_type_perturbative_order = -1;
  type_perturbative_order_counter = 0;
  collection_type_perturbative_order.resize(1, "all");
  for (int i = 0; i < user_variable.size(); i++){
    if (user_variable[i] == "type_perturbative_order"){
      int flag = -1;
      for (int is = 0; is < type_perturbative_order_counter + 1; is++){if (user_value[i] == collection_type_perturbative_order[is]){flag = is; break;}}
      logger << LOG_DEBUG_VERBOSE << "flag = " << flag << endl;
      if (flag == -1){
	type_perturbative_order_counter++;
	collection_type_perturbative_order.push_back(user_value[i]);
	present_type_perturbative_order = type_perturbative_order_counter;
      }
      else {present_type_perturbative_order = flag;}
    }
  }
  if (present_type_perturbative_order > -1){csi.type_perturbative_order = collection_type_perturbative_order[present_type_perturbative_order];}
  logger << LOG_DEBUG << "number of type_perturbative_order's: " << collection_type_perturbative_order.size() << endl;
  for (int is = 0; is < type_perturbative_order_counter + 1; is++){logger << LOG_DEBUG << "collection_type_perturbative_order[" << setw(2) << is << "] = " << collection_type_perturbative_order[is] << endl;}
  logger << LOG_DEBUG << "present_type_perturbative_order:     " << present_type_perturbative_order << endl;
  logger << LOG_DEBUG << "type_perturbative_order:             " << csi.type_perturbative_order << endl;
  logger.newLine(LOG_DEBUG);

  ///////////////////////
  //  beam parameters  //
  ///////////////////////

  contribution_LHAPDFname.resize(collection_type_perturbative_order.size());
  contribution_LHAPDFsubset.resize(collection_type_perturbative_order.size());


  //////////////////////////////////////////////////////////////////////
  //  determination of particular type_contribution to be calculated  //
  //////////////////////////////////////////////////////////////////////

  present_type_contribution = -1;
  type_contribution_counter = 0;
  collection_type_contribution.resize(1, "all");
  for (int i = 0; i < user_variable.size(); i++){
    if (user_variable[i] == "type_contribution"){
      int flag = -1;
      for (int is = 0; is < type_contribution_counter + 1; is++){if (user_value[i] == collection_type_contribution[is]){flag = is; break;}}
      logger << LOG_DEBUG_VERBOSE << "flag = " << flag << endl;
      if (flag == -1){
	type_contribution_counter++;
	collection_type_contribution.push_back(user_value[i]);
	present_type_contribution = type_contribution_counter;
      }
      else {present_type_contribution = flag;}
    }
  }
  if (present_type_contribution > -1){csi.type_contribution = collection_type_contribution[present_type_contribution];}
  logger << LOG_DEBUG << "number of type_contributions: " << collection_type_contribution.size() << endl;
  for (int is = 0; is < type_contribution_counter + 1; is++){logger << LOG_DEBUG << "collection_type_contribution[" << setw(2) << is << "] = " << collection_type_contribution[is] << endl;}
  logger << LOG_DEBUG << "present_type_contribution:    " << present_type_contribution << endl;
  logger << LOG_DEBUG << "type_contribution:            " << csi.type_contribution << endl;
  logger.newLine(LOG_DEBUG);

  ////////////////////////////////////////////////////////////////////
  //  determination of particular type_correction to be calculated  //
  ////////////////////////////////////////////////////////////////////

  present_type_correction = -1;
  type_correction_counter = 0;
  collection_type_correction.resize(1, "all");
  for (int i = 0; i < user_variable.size(); i++){
    if (user_variable[i] == "type_correction"){
      int flag = -1;
      for (int is = 0; is < type_correction_counter + 1; is++){if (user_value[i] == collection_type_correction[is]){flag = is; break;}}
      logger << LOG_DEBUG_VERBOSE << "flag = " << flag << endl;
      if (flag == -1){
	type_correction_counter++;
	collection_type_correction.push_back(user_value[i]);
	present_type_correction = type_correction_counter;
      }
      else {present_type_correction = flag;}
    }
  }
  if (present_type_correction > -1){csi.type_correction = collection_type_correction[present_type_correction];}
  logger << LOG_DEBUG << "number of type_correction's: " << collection_type_correction.size() << endl;
  for (int is = 0; is < type_correction_counter + 1; is++){logger << LOG_DEBUG << "collection_type_correction[" << setw(2) << is << "] = " << collection_type_correction[is] << endl;}
  logger << LOG_DEBUG << "present_type_correction:     " << present_type_correction << endl;
  logger << LOG_DEBUG << "type_correction:             " << csi.type_correction << endl;
  logger.newLine(LOG_DEBUG);



  ////////////////////////////////////////
  //  scale variation parameters - TSV  //
  ////////////////////////////////////////

  switch_TSV = 0;
  for (int i = 0; i < user_variable.size(); i++){
    if (user_variable[i] == "switch_TSV"){switch_TSV = atoi(user_value[i].c_str());}
  }

  int present_scaleset = -1;
  int counter_scaleset = -1;
  name_set_TSV.resize(0);
  for (int i = 0; i < user_variable.size(); i++){
    if (user_variable[i] == "name_set_TSV" || user_variable[i] == "scaleset"){
      int flag = -1;
      for (int i_s = 0; i_s < counter_scaleset + 1; i_s++){
	if (user_value[i] == name_set_TSV[i_s]){flag = i_s; break;}
      }
      logger << LOG_DEBUG_VERBOSE << "flag = " << flag << endl;
      if (flag == -1){
	counter_scaleset++;
	name_set_TSV.push_back(user_value[i]);
	present_scaleset = counter_scaleset;
      }
      else {present_scaleset = flag;}
    }
  }
  
  logger << LOG_INFO << "number of scalesets: " << name_set_TSV.size() << endl;
  for (int i_s = 0; i_s < name_set_TSV.size(); i_s++){
    logger << LOG_INFO << "name_set_TSV[" << setw(2) << i_s << "] = " << name_set_TSV[i_s] << endl;
  }
  logger.newLine(LOG_DEBUG);

  n_set_TSV = name_set_TSV.size();
  logger << LOG_DEBUG << "n_set_TSV = " << n_set_TSV << endl;

  central_scale_TSV.resize(n_set_TSV);
  central_scale_ren_TSV.resize(n_set_TSV);
  central_scale_fact_TSV.resize(n_set_TSV);
  relative_central_scale_TSV.resize(n_set_TSV);
  relative_central_scale_ren_TSV.resize(n_set_TSV);
  relative_central_scale_fact_TSV.resize(n_set_TSV);

  factor_scale_TSV.resize(n_set_TSV);
  factor_scale_ren_TSV.resize(n_set_TSV);
  factor_scale_fact_TSV.resize(n_set_TSV);

  dynamic_scale_TSV.resize(n_set_TSV);
  dynamic_scale_ren_TSV.resize(n_set_TSV);
  dynamic_scale_fact_TSV.resize(n_set_TSV);
  min_qTcut_TSV.resize(n_set_TSV);
  max_qTcut_TSV.resize(n_set_TSV);
  min_qTcut_distribution_TSV.resize(n_set_TSV);
  max_qTcut_distribution_TSV.resize(n_set_TSV);

  no_central_scale_ren_TSV.resize(n_set_TSV);
  // !!! no input parameter !!!

  no_central_scale_fact_TSV.resize(n_set_TSV);
  // !!! no input parameter !!!

  double empty_double = -1.;
  double empty_int = -1;
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    central_scale_TSV[i_s] = empty_double;
    central_scale_ren_TSV[i_s] = empty_double;
    central_scale_fact_TSV[i_s] = empty_double;
    relative_central_scale_TSV[i_s] = empty_double;
    relative_central_scale_ren_TSV[i_s] = empty_double;
    relative_central_scale_fact_TSV[i_s] = empty_double;
    dynamic_scale_TSV[i_s] = empty_int;
    dynamic_scale_ren_TSV[i_s] = empty_int;
    dynamic_scale_fact_TSV[i_s] = empty_int;
    factor_scale_TSV[i_s] = empty_int;
    factor_scale_ren_TSV[i_s] = empty_int;
    factor_scale_fact_TSV[i_s] = empty_int;
    min_qTcut_TSV[i_s] = empty_double;
    max_qTcut_TSV[i_s] = empty_double;
    min_qTcut_distribution_TSV[i_s] = empty_double;
    max_qTcut_distribution_TSV[i_s] = empty_double;
  }



  // to deal with differences between two selected scale choices, acting as a local subtraction for each other:

  int present_scalediffset = -1;
  int counter_scalediffset = -1;
  name_diff_set_TSV.resize(0);
  for (int i = 0; i < user_variable.size(); i++){
    if (user_variable[i] == "name_diff_set_TSV"){
      int flag = -1;
      for (int i_s = 0; i_s < counter_scalediffset + 1; i_s++){
	if (user_value[i] == name_diff_set_TSV[i_s]){flag = i_s; break;}
      }
      logger << LOG_DEBUG_VERBOSE << "flag = " << flag << endl;
      if (flag == -1){
	counter_scalediffset++;
	name_diff_set_TSV.push_back(user_value[i]);
	present_scalediffset = counter_scalediffset;
      }
      else {present_scalediffset = flag;}
    }
  }
  
  n_diff_set_TSV = name_diff_set_TSV.size();

  /*
  logger << LOG_INFO << "number of scalediffsets: " << n_diff_set_TSV << endl;
  for (int i_s = 0; i_s < name_diff_set_TSV.size(); i_s++){
    logger << LOG_INFO << "name_diff_set_TSV[" << setw(2) << i_s << "] = " << name_diff_set_TSV[i_s] << endl;
  }
  logger.newLine(LOG_INFO);
  */

  name_diff_set_plus_TSV.resize(n_diff_set_TSV, "");
  name_diff_set_minus_TSV.resize(n_diff_set_TSV, "");
  no_diff_set_plus_TSV.resize(n_diff_set_TSV, -1);
  no_diff_set_minus_TSV.resize(n_diff_set_TSV, -1);

  n_extended_set_TSV = n_set_TSV + n_diff_set_TSV;

  n_scale_TSV.resize(n_extended_set_TSV);
  n_scale_ren_TSV.resize(n_extended_set_TSV);
  n_scale_fact_TSV.resize(n_extended_set_TSV);

  switch_distribution_TSV.resize(n_extended_set_TSV);
  switch_moment_TSV.resize(n_extended_set_TSV);

  max_n_integrand_TSV.resize(n_extended_set_TSV);
  // !!! no input parameter !!!

    logger << LOG_DEBUG_VERBOSE << "SSSSS   n_extended_set_TSV = " << n_extended_set_TSV << endl;
  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
    n_scale_TSV[i_s] = empty_int;
    n_scale_ren_TSV[i_s] = empty_int;
    n_scale_fact_TSV[i_s] = empty_int;

    switch_distribution_TSV[i_s] = empty_int;
    switch_moment_TSV[i_s] = empty_int;
    logger << LOG_DEBUG_VERBOSE << "SSSSS   switch_distribution_TSV[" << i_s << "] = " << switch_distribution_TSV[i_s] << endl;
  }

  name_extended_set_TSV.resize(n_extended_set_TSV, "");
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    name_extended_set_TSV[i_s] = name_set_TSV[i_s];
  }
  for (int i_s = 0; i_s < n_diff_set_TSV; i_s++){
    name_extended_set_TSV[n_set_TSV + i_s] = name_diff_set_TSV[i_s];
  }

  switch_reference = "";
  
  name_reference_TSV = "" ;
  no_reference_TSV = -1;
  no_scale_ren_reference_TSV = -1;
  no_scale_fact_reference_TSV = -1;
  no_qTcut_reference_TSV = 0;



  // end



  logger << LOG_DEBUG << "esi.object_list.size() = " << esi.object_list.size() << endl;
  for (int i_e = 0; i_e < esi.object_list.size(); i_e++){
    logger << LOG_DEBUG << "esi.object_list[" << i_e << "] = " << esi.object_list[i_e] << endl;
  }
  //  oset.initialization_object_generic();
  // to be proceeded !!!
 
  /////////////////////////////////////////////////////////////////////////////////////
  //  check for user-defined particles -> needed before further proceeding of input  //
  /////////////////////////////////////////////////////////////////////////////////////

  for (int i = 0; i < user_variable.size(); i++){
    if (user_variable[i] == "user_particle"){
      if (user.particle_map[user_variable_additional[i]] == 0){
	user.particle_name.push_back(user_variable_additional[i]);
	user.particle_map[user_variable_additional[i]] = user.particle_name.size() - 1;
	user.particle_value.push_back(user_value[i]);
      }
      else {
	user.particle_value[user.particle_map[user_variable_additional[i]]] = user_value[i];
      }
    }
  }

  for (int i_p = 0; i_p < user.particle_name.size(); i_p++){
    logger << LOG_DEBUG << "user.particle_name[" << i_p << "] = " << user.particle_name[i_p] << endl;
  }

  esi.define_specific_object_list(user);

  logger << LOG_DEBUG << "After adding user-defined particles:" << endl;

  logger << LOG_DEBUG << "esi.object_list.size() = " << esi.object_list.size() << endl;
  for (int i_e = 0; i_e < esi.object_list.size(); i_e++){
    logger << LOG_DEBUG << "esi.object_list[" << i_e << "] = " << esi.object_list[i_e] << endl;
  }




  ////////////////////////////////////////////////////////////
  //  proceed the input from all 'file_paramter.dat' files  //
  ////////////////////////////////////////////////////////////

  int temp_type_perturbative_order = -1;
  int temp_type_contribution = -1;
  int temp_type_correction = -1;

  for (int i = 0; i < user_variable.size(); i++){
    logger << LOG_DEBUG_VERBOSE << "user_variable[" << i << "] = " << user_variable[i] << endl;

    //////////////////////////
    //  unknown parameters  //
    //////////////////////////

    if (user_variable[i] == "path_to_main"){path_to_main = user_value[i];}
    // ???

    else if (user_variable[i] == "ckm_choice"){ckm_choice = atoi(user_value[i].c_str());}
    // not relevant by now - but maybe later...

    else if (user_variable[i] == "n_moments"){n_moments = atoi(user_value[i].c_str());}
    // ??? should be replaced by input file !!!

    //      else if (user_variable[i] == "n_madgraph"){n_madgraph = atoi(user_value[i].c_str());} // not needed anymore ???

    else if (user_variable[i] == "test_output"){test_output = atoi(user_value[i].c_str());}

    //    else if (user_variable[i] == "switch_alpha_CMS"){switch_alpha_CMS = atoi(user_value[i].c_str());}


    //////////////////////////////////////////////////////////////
    //  selection of process and contribution to be calculated  //
    //////////////////////////////////////////////////////////////

    else if (user_variable[i] == "process_class"){csi.process_class = user_value[i];}
    else if (user_variable[i] == "decay"){csi.decay.push_back(user_value[i].c_str());}
    else if (user_variable[i] == "type_perturbative_order"){
      temp_type_perturbative_order = -1;
      for (int i_to = 0; i_to < collection_type_perturbative_order.size(); i_to++){
	if (user_value[i] == "all"){temp_type_perturbative_order = present_type_perturbative_order; break;}
	else if (user_value[i] == collection_type_perturbative_order[i_to]){temp_type_perturbative_order = i_to; break;}
      }
      logger << LOG_DEBUG << "user_variable[" << i << "] = " << setw(32) << left << user_variable[i] << " = " << setw(32) << user_value[i] << setw(32) << "temp_type_perturbative_order = " << temp_type_perturbative_order << endl;
    }
    else if (user_variable[i] == "type_contribution"){
      temp_type_contribution = -1;
      for (int i_tc = 0; i_tc < collection_type_contribution.size(); i_tc++){
	if (user_value[i] == "all"){temp_type_contribution = present_type_contribution; break;}
	else if (user_value[i] == collection_type_contribution[i_tc]){temp_type_contribution = i_tc; break;}
      }
      logger << LOG_DEBUG << "user_variable[" << i << "] = " << setw(32) << left << user_variable[i] << " = " << setw(32) << user_value[i] << setw(32) << "temp_type_contribution = " << temp_type_contribution << endl;
    }
    else if (user_variable[i] == "type_correction"){
      temp_type_correction = -1;
      for (int i_tc = 0; i_tc < collection_type_correction.size(); i_tc++){
	if (user_value[i] == "all"){temp_type_correction = present_type_correction; break;}
	else if (user_value[i] == collection_type_correction[i_tc]){temp_type_correction = i_tc; break;}
      }
      logger << LOG_DEBUG << "user_variable[" << i << "] = " << setw(32) << left << user_variable[i] << " = " << setw(32) << user_value[i] << setw(32) << "temp_type_correction = " << temp_type_correction << endl;
    }

    ////////////////////////////////////////////////
    //  switches to steer calculation of results  //
    ////////////////////////////////////////////////

    else if (user_variable[i] == "switch_distribution"){switch_distribution = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "distributions"){switch_distribution = atoi(user_value[i].c_str());} // to be removed...
    else if (user_variable[i] == "switch_moment"){switch_moment = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_result"){switch_result = atoi(user_value[i].c_str());}
    
    ///////////////////////////////////////////////
    //  switches to steer output of calculation  //
    ///////////////////////////////////////////////

    else if (user_variable[i] == "switch_output_execution"){switch_output_execution = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_integration"){switch_output_integration = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_maxevent"){switch_output_maxevent = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_cancellation_check"){switch_output_cancellation_check = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_comparison"){switch_output_comparison = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_gnuplot"){switch_output_gnuplot = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_proceeding"){switch_output_proceeding = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_weights"){switch_output_weights = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_result"){switch_output_result = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_moment"){switch_output_moment = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_distribution"){switch_output_distribution = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_testpoint"){switch_output_testpoint = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_output_cutinfo"){switch_output_cutinfo = atoi(user_value[i].c_str());}

    else if (user_variable[i] == "switch_console_output_runtime"){switch_console_output_runtime = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_console_output_tau_0"){switch_console_output_tau_0 = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_console_output_techcut_RA"){switch_console_output_techcut_RA = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_console_output_phasespace_issue"){switch_console_output_phasespace_issue = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_console_output_ME2_issue"){switch_console_output_ME2_issue = atoi(user_value[i].c_str());}

    else if (user_variable[i] == "switch_testcut"){switch_testcut = atoi(user_value[i].c_str());}
    
    ////////////////////////////////////////////////////////////////////////
    //  unit of calculation output / result output / distribution output  //
    ////////////////////////////////////////////////////////////////////////

    else if (user_variable[i] == "unit_calculation"){unit_calculation = user_value[i];}
    else if (user_variable[i] == "unit_result"){unit_result = user_value[i];}
    else if (user_variable[i] == "unit_distribution"){unit_distribution = user_value[i];}

    ///////////////////////
    //  beam parameters  //
    ///////////////////////

    else if (user_variable[i] == "E"){E = atof(user_value[i].c_str());}
    else if (user_variable[i] == "coll_choice"){coll_choice = atoi(user_value[i].c_str());}
    //    else if (user_variable[i] == "pdf_set"){pdf_set = atoi(user_value[i].c_str());}
    // ???

    //    else if (user_variable[i] == "pdf_content_modify"){pdf_content_modify = atoi(user_value[i].c_str());}
    // ???

    else if (user_variable[i] == "pdf_selection"){
      //  Reasonable non-overlapping combinations (flavour-blind):
      //  qq,qxqx,gg,aa,qqx,qxq,gq,qg,gqx,qxg,aq,qa,aqx,qxa,ga,ag
      //  Q := q + qx  ->  e.g.  QQ = qqx + qxq + qq + qxqx
      //  Q/q/qx can be replaced by defined flavour.
      //  D := d + dx, U := u + ux, S := s + sx, C := c + cx, B := b + bx, T := t + tx
      if (user_value[i] == "clear"){pdf_selection.clear();}
      else {pdf_selection.push_back(user_value[i]);}
    }
    else if (user_variable[i] == "pdf_disable"){
      //  Partons whose pdf's can be switched off:
      //  g,a,Q,q,qx,D,d,dx,U,u,ux,S,s,sx,C,c,cx,B,b,bx,T,t,tx
      if (user_value[i] == "clear"){pdf_disable.clear();}
      else {pdf_disable.push_back(user_value[i]);}
    }
    else if (user_variable[i] == "LHAPDFname"){
      contribution_LHAPDFname[temp_type_perturbative_order] = user_value[i];
      logger << LOG_DEBUG << "user_variable[" << i << "] = " << setw(32) << left << user_variable[i] << " = " << setw(32) << user_value[i] << setw(32) << "contribution_LHAPDFname[temp_type_perturbative_order = " << temp_type_perturbative_order << "] = " << contribution_LHAPDFname[temp_type_perturbative_order] << endl;
    }
    else if (user_variable[i] == "LHAPDFsubset"){
      contribution_LHAPDFsubset[temp_type_perturbative_order] = atoi(user_value[i].c_str());
      logger << LOG_DEBUG << "user_variable[" << i << "] = " << setw(32) << left << user_variable[i] << " = " << setw(32) << user_value[i] << setw(32) << "contribution_LHAPDFsubset[temp_type_perturbative_order = " << temp_type_perturbative_order << "] = " << contribution_LHAPDFsubset[temp_type_perturbative_order] << endl;
    }

    ////////////////////////////////////////
    //  jet and jet-algorithm parameters  //
    ////////////////////////////////////////

    else if (user_variable[i] == "jet_algorithm"){jet_algorithm = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "jet_R_definition"){jet_R_definition = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "jet_R"){jet_R = atof(user_value[i].c_str());}
    else if (user_variable[i] == "parton_y_max"){parton_y_max = atof(user_value[i].c_str());}
    else if (user_variable[i] == "parton_eta_max"){parton_eta_max = atof(user_value[i].c_str());}
    else if (user_variable[i] == "jet_algorithm_selection"){
      //  Partons that can be selected for undergoing the jet algorithm:
      //  g,a,Q,q,qx,D,d,dx,U,u,ux,S,s,sx,C,c,cx,B,b,bx,T,t,tx
      if (user_value[i] == "clear"){jet_algorithm_selection.clear();}
      else {jet_algorithm_selection.push_back(user_value[i]);}
    }
    else if (user_variable[i] == "jet_algorithm_disable"){
      //  Partons that can be excluded from the jet algorithm:
      //  g,a,Q,q,qx,D,d,dx,U,u,ux,S,s,sx,C,c,cx,B,b,bx,T,t,tx
      if (user_value[i] == "clear"){jet_algorithm_disable.clear();}
      else {jet_algorithm_disable.push_back(user_value[i]);}
    }
    else if (user_variable[i] == "N_f"){N_f = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "N_f_active"){N_f_active = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "N_quarks"){N_quarks = atoi(user_value[i].c_str());}
    // !! actually set by the Standard Model

    else if (user_variable[i] == "N_nondecoupled"){N_nondecoupled = atoi(user_value[i].c_str());}
    // !! actually internal parameter

    ///////////////////////////////////////
    //  photon-recombination parameters  //
    ///////////////////////////////////////

    else if (user_variable[i] == "photon_recombination"){photon_recombination = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "photon_R_definition"){photon_R_definition = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "photon_R"){photon_R = atof(user_value[i].c_str());}
    else if (user_variable[i] == "photon_E_threshold_ratio"){photon_E_threshold_ratio = atof(user_value[i].c_str());}
    else if (user_variable[i] == "photon_jet_algorithm"){photon_jet_algorithm = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "photon_recombination_selection"){
      //  Partons that can be selected for undergoing the photon recombination:
      if (user_value[i] == "clear"){photon_recombination_selection.clear();}
      else {photon_recombination_selection.push_back(user_value[i]);}
    }
    else if (user_variable[i] == "photon_recombination_disable"){
      //  Partons that can be excluded from the photon recombination:
      if (user_value[i] == "clear"){photon_recombination_disable.clear();}
      else {photon_recombination_disable.push_back(user_value[i]);}
    }
    else if (user_variable[i] == "photon_photon_recombination"){photon_photon_recombination = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "photon_photon_recombination_R"){photon_photon_recombination_R = atof(user_value[i].c_str());}


    
    ///////////////////////////////////
    //  photon-isolation parameters  //
    ///////////////////////////////////
    
    else if (user_variable[i] == "frixione_isolation"){frixione_isolation = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "frixione_n"){frixione_n = atof(user_value[i].c_str());}
    else if (user_variable[i] == "frixione_epsilon"){frixione_epsilon = atof(user_value[i].c_str());}
    else if (user_variable[i] == "frixione_fixed_ET_max"){frixione_fixed_ET_max = atof(user_value[i].c_str());}
    else if (user_variable[i] == "frixione_delta_0"){frixione_delta_0 = atof(user_value[i].c_str());}
    else if (user_variable[i] == "frixione_jet_removal"){frixione_jet_removal = atoi(user_value[i].c_str());}

    /////////////////////////////////
    //  qT-subtraction parameters  //
    /////////////////////////////////

    else if (user_variable[i] == "switch_qTcut"){switch_qTcut = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "n_qTcut"){n_qTcut = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "min_qTcut"){min_qTcut = atof(user_value[i].c_str());}
    else if (user_variable[i] == "step_qTcut"){step_qTcut = atof(user_value[i].c_str());}
    else if (user_variable[i] == "max_qTcut"){max_qTcut = atof(user_value[i].c_str());}
    else if (user_variable[i] == "binning_qTcut"){binning_qTcut = user_value[i];}
    else if (user_variable[i] == "selection_qTcut"){
      if (user_value[i] == "clear"){selection_qTcut = "";}
      else {selection_qTcut = user_value[i];}
    }
    else if (user_variable[i] == "selection_qTcut_distribution"){
      if (user_value[i] == "clear"){selection_qTcut_distribution = "";}
      else {selection_qTcut_distribution = user_value[i];}
    }
    else if (user_variable[i] == "selection_no_qTcut_distribution"){
      if (user_value[i] == "clear"){selection_no_qTcut_distribution = "";}
      else {selection_no_qTcut_distribution = user_value[i];}
    }
    else if (user_variable[i] == "selection_qTcut_result"){
      if (user_value[i] == "clear"){selection_qTcut_result = "";}
      else {selection_qTcut_result = user_value[i];}
    }
    else if (user_variable[i] == "selection_no_qTcut_result"){
      if (user_value[i] == "clear"){selection_no_qTcut_result = "";}
      else {selection_no_qTcut_result = user_value[i];}
    }
    else if (user_variable[i] == "selection_qTcut_integration"){
      if (user_value[i] == "clear"){selection_qTcut_integration = "";}
      else {selection_qTcut_integration = user_value[i];}
    }
    else if (user_variable[i] == "selection_no_qTcut_integration"){
      if (user_value[i] == "clear"){selection_no_qTcut_integration = "";}
      else {selection_no_qTcut_integration = user_value[i];}
    }

    
    ///////////////////////////////////////////
    //  N-jettiniess-subtraction parameters  //
    ///////////////////////////////////////////

    else if (user_variable[i] == "switch_NJcut"){switch_NJcut = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_NJcut_axes"){switch_NJcut_axes = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_NJcut_axes_energy"){switch_NJcut_axes_energy = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "switch_NJcut_measure"){switch_NJcut_measure = atoi(user_value[i].c_str());}

    

    ////////////////////////////////////////////////////
    //  technical switches for selected contributions //
    ////////////////////////////////////////////////////

    else if (user_variable[i] == "switch_polenorm"){switch_polenorm = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "polenorm_switch"){switch_polenorm = atoi(user_value[i].c_str());} // to be removed...

   else if (user_variable[i] == "switch_old_qT_version"){switch_old_qT_version = atoi(user_value[i].c_str());}

    ///////////////////////////////
    //  user-defined parameters  //
    ///////////////////////////////

    else if (user_variable[i] == "user_switch"){
      if (user.switch_map[user_variable_additional[i]] == 0){
	user.switch_name.push_back(user_variable_additional[i]);
	user.switch_map[user_variable_additional[i]] = user.switch_name.size() - 1;
	user.switch_value.push_back(atoi(user_value[i].c_str()));
      }
      else {
	user.switch_value[user.switch_map[user_variable_additional[i]]] = atoi(user_value[i].c_str());
      }
    }

    else if (user_variable[i] == "user_cut"){
      if (user.cut_map[user_variable_additional[i]] == 0){
	user.cut_name.push_back(user_variable_additional[i]);
	user.cut_map[user_variable_additional[i]] = user.cut_name.size() - 1;
	user.cut_value.push_back(atof(user_value[i].c_str()));
      }
      else {
	user.cut_value[user.cut_map[user_variable_additional[i]]] = atof(user_value[i].c_str());
      }
    }

    else if (user_variable[i] == "user_int"){
      if (user.int_map[user_variable_additional[i]] == 0){
	user.int_name.push_back(user_variable_additional[i]);
	user.int_map[user_variable_additional[i]] = user.int_name.size() - 1;
	user.int_value.push_back(atoi(user_value[i].c_str()));
      }
      else {
	user.int_value[user.int_map[user_variable_additional[i]]] = atoi(user_value[i].c_str());
      }
    }

    else if (user_variable[i] == "user_double"){
      if (user.double_map[user_variable_additional[i]] == 0){
	user.double_name.push_back(user_variable_additional[i]);
	user.double_map[user_variable_additional[i]] = user.double_name.size() - 1;
	user.double_value.push_back(atof(user_value[i].c_str()));
      }
      else {
	user.double_value[user.double_map[user_variable_additional[i]]] = atof(user_value[i].c_str());
      }
    }

    else if (user_variable[i] == "user_string"){
      if (user.string_map[user_variable_additional[i]] == 0){
	user.string_name.push_back(user_variable_additional[i]);
	user.string_map[user_variable_additional[i]] = user.string_name.size() - 1;
	user.string_value.push_back(user_value[i]);
      }
      else {
	user.string_value[user.string_map[user_variable_additional[i]]] = atof(user_value[i].c_str());
      }
    }

    ////////////////////////////////////////////////////
    //  basic event-selection criteria for particles  //
    ////////////////////////////////////////////////////

    else if (user_variable[i] == "n_observed_min"){
      if (esi.observed_object[user_variable_additional[i]] > 0){
	esi.pda[esi.observed_object[user_variable_additional[i]]].n_observed_min = atoi(user_value[i].c_str());
      }
      else {
	logger << LOG_FATAL << "Object " << user_variable_additional[i] << " has not been defined (n_observed_min) !" << endl;
	//	exit(1);
      }
    }
    else if (user_variable[i] == "n_observed_max"){
      if (esi.observed_object[user_variable_additional[i]] > 0){
	esi.pda[esi.observed_object[user_variable_additional[i]]].n_observed_max = atoi(user_value[i].c_str());
      }
      else {
	logger << LOG_FATAL << "Object " << user_variable_additional[i] << " has not been defined (n_observed_max) !" << endl;
	//	exit(1);
      }
    }
    else if (user_variable[i] == "define_pT"){
      if (esi.observed_object[user_variable_additional[i]] > 0){
	esi.pda[esi.observed_object[user_variable_additional[i]]].define_pT = atof(user_value[i].c_str());
      }
      else {
	logger << LOG_FATAL << "Object " << user_variable_additional[i] << " has not been defined (define_pT) !" << endl;
	//	exit(1);
      }
    }
    else if (user_variable[i] == "define_ET"){
      if (esi.observed_object[user_variable_additional[i]] > 0){
	esi.pda[esi.observed_object[user_variable_additional[i]]].define_ET = atof(user_value[i].c_str());
      }
      else {
	logger << LOG_FATAL << "Object " << user_variable_additional[i] << " has not been defined (define_ET) !" << endl;
	//	exit(1);
      }
    }
    else if (user_variable[i] == "define_eta"){
      if (esi.observed_object[user_variable_additional[i]] > 0){
	esi.pda[esi.observed_object[user_variable_additional[i]]].define_eta = atof(user_value[i].c_str());
      }
      else {
	logger << LOG_FATAL << "Object " << user_variable_additional[i] << " has not been defined (define_eta) !" << endl;
	//	exit(1);
      }
    }
    else if (user_variable[i] == "define_y"){
      if (esi.observed_object[user_variable_additional[i]] > 0){
	esi.pda[esi.observed_object[user_variable_additional[i]]].define_y = atof(user_value[i].c_str());
      }
      else {
	logger << LOG_FATAL << "Object " << user_variable_additional[i] << " has not been defined (define_y) !" << endl;
	//	exit(1);
      }
    }

    else if (user_variable[i] == "fiducial_cut"){
      logger << LOG_INFO << "fiducial_cut[" << i << "] = " << user_value[i] << endl;
      name_fiducial_cut.push_back(user_value[i]);
    }



    ///////////////////////////////
    //  scale-choice parameters  //
    ///////////////////////////////

    else if (user_variable[i] == "scale_fact"){scale_fact = atof(user_value[i].c_str());}
    else if (user_variable[i] == "scale_ren"){scale_ren = atof(user_value[i].c_str());}
    else if (user_variable[i] == "dynamic_scale"){dynamic_scale = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "prefactor_reference"){prefactor_reference = atof(user_value[i].c_str());}

    ///////////////////////////////////////
    //  scale variation parameters - CV  //
    ///////////////////////////////////////

    else if (user_variable[i] == "switch_CV"){switch_CV = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "variation_CV"){switch_CV = atoi(user_value[i].c_str());} // !!! to be removed
    else if (user_variable[i] == "n_scales_CV"){n_scales_CV = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "dynamic_scale_CV"){dynamic_scale_CV = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "variation_mu_ren_CV"){variation_mu_ren_CV = atoi(user_value[i].c_str());}
  // ???

    else if (user_variable[i] == "variation_mu_fact_CV"){variation_mu_fact_CV = atoi(user_value[i].c_str());}
  // ???

    else if (user_variable[i] == "variation_factor_CV"){variation_factor_CV = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "central_scale_CV"){central_scale_CV = atoi(user_value[i].c_str());}
  // ???

    else if (user_variable[i] == "prefactor_CV"){prefactor_CV = atof(user_value[i].c_str());}

    //////////////////////////////
    //  integration parameters  //
    //////////////////////////////

    else if (user_variable[i] == "sigma_normalization"){sigma_normalization = atof(user_value[i].c_str());}
    else if (user_variable[i] == "sigma_LO"){sigma_normalization = atof(user_value[i].c_str());} // to be removed...
    else if (user_variable[i] == "zwahl"){zwahl = atoi(user_value[i].c_str());}

    ////////////////////////////////////////////////////////////////////////////////////////
    //  parameters directly forwarded to OpenLoops - after the default settings are done  //
    ////////////////////////////////////////////////////////////////////////////////////////

    else if (user_variable[i] == "OL"){
      if (user_variable_additional[i] != ""){
	OL_parameter.push_back(user_variable_additional[i]);
	OL_value.push_back(user_value[i]);
      }
    }

    ////////////////////////////////////////
    //  scale variation parameters - TSV  //
    ////////////////////////////////////////

    else if (user_variable[i] == "switch_TSV"){switch_TSV = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "name_set_TSV" || user_variable[i] == "scaleset"){ // to be removed...
      for (int i_s = 0; i_s < n_set_TSV; i_s++){
	if (user_value[i] == name_set_TSV[i_s]){present_scaleset = i_s; break;}
      }
    }
    else if (user_variable[i] == "central_scale_TSV"){central_scale_TSV[present_scaleset] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "central_scale_ren_TSV"){central_scale_ren_TSV[present_scaleset] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "central_scale_fact_TSV"){central_scale_fact_TSV[present_scaleset] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "relative_central_scale_TSV"){relative_central_scale_TSV[present_scaleset] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "relative_central_scale_ren_TSV"){relative_central_scale_ren_TSV[present_scaleset] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "relative_central_scale_fact_TSV"){relative_central_scale_fact_TSV[present_scaleset] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "n_scale_TSV"){n_scale_TSV[present_scaleset] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "n_scale_ren_TSV"){n_scale_ren_TSV[present_scaleset] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "n_scale_fact_TSV"){n_scale_fact_TSV[present_scaleset] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "factor_scale_TSV"){factor_scale_TSV[present_scaleset] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "factor_scale_ren_TSV"){factor_scale_ren_TSV[present_scaleset] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "factor_scale_fact_TSV"){factor_scale_fact_TSV[present_scaleset] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "dynamic_scale_TSV"){dynamic_scale_TSV[present_scaleset] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "dynamic_scale_ren_TSV"){dynamic_scale_ren_TSV[present_scaleset] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "dynamic_scale_fact_TSV"){dynamic_scale_fact_TSV[present_scaleset] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "min_qTcut_TSV"){min_qTcut_TSV[present_scaleset] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "max_qTcut_TSV"){max_qTcut_TSV[present_scaleset] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "switch_distribution_TSV"){switch_distribution_TSV[present_scaleset] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "min_qTcut_distribution_TSV"){min_qTcut_TSV[present_scaleset] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "max_qTcut_distribution_TSV"){max_qTcut_TSV[present_scaleset] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "switch_moment_TSV"){switch_moment_TSV[present_scaleset] = atoi(user_value[i].c_str());}

    else if (user_variable[i] == "name_diff_set_TSV"){ // to be removed...
      for (int i_s = 0; i_s < n_diff_set_TSV; i_s++){
	if (user_value[i] == name_diff_set_TSV[i_s]){present_scalediffset = i_s; break;}
      }
    }
    else if (user_variable[i] == "name_diff_set_plus_TSV"){name_diff_set_plus_TSV[present_scalediffset] = user_value[i];}
    else if (user_variable[i] == "name_diff_set_minus_TSV"){name_diff_set_minus_TSV[present_scalediffset] = user_value[i];}

    else if (user_variable[i] == "switch_reference"){switch_reference = user_value[i];}
    else if (user_variable[i] == "name_reference_TSV"){name_reference_TSV = user_value[i];}
    else if (user_variable[i] == "no_reference_TSV"){no_reference_TSV = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "no_scale_ren_reference_TSV"){no_scale_ren_reference_TSV = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "no_scale_fact_reference_TSV"){no_scale_fact_reference_TSV = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "no_qTcut_reference_TSV"){no_qTcut_reference_TSV = atoi(user_value[i].c_str());}

    ///////////////////////////////////////////////
    //  selection of run_mode (grid, time, run)  //
    ///////////////////////////////////////////////
    else if (user_variable[i] == "run_mode"){run_mode = user_value[i];}






    ///////////////////////////////////////////////////////////////////////////////////////
    //  the following input is only used if the 'present' contribution                   //
    //  matches the 'temp' contribution selected at the respective point in input file:  //
    // - type_perturbative_order                                                         //
    // - type_contribution                                                               //
    // - type_correction                                                                 //
    ///////////////////////////////////////////////////////////////////////////////////////

    else if ((temp_type_perturbative_order == present_type_perturbative_order && 
	      temp_type_contribution == present_type_contribution && 
	      temp_type_correction == present_type_correction)){
      logger << LOG_DEBUG << "user_variable[" << i << "] = " << left << setw(32) << user_variable[i] << " = " << setw(32) << user_value[i] << endl;

      //////////////////////////////////////////////////
      //  selection of contribution to be calculated  //
      //////////////////////////////////////////////////

      if      (user_variable[i] == "contribution_order_alpha_s"){csi.contribution_order_alpha_s = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "contribution_order_alpha_e"){csi.contribution_order_alpha_e = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "contribution_order_interference"){csi.contribution_order_interference = atoi(user_value[i].c_str());}

      ////////////////////////////////////////////////////
      //  technical switches for selected contributions //
      ////////////////////////////////////////////////////

      else if (user_variable[i] == "switch_KP"){switch_KP = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "switch_VI"){switch_VI = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "switch_VI_bosonic_fermionic"){switch_VI = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "switch_H1gg"){switch_H1gg = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "switch_H2"){switch_H2 = atoi(user_value[i].c_str());}
      
      else if (user_variable[i] == "switch_CM"){switch_CM = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "switch_OL"){switch_OL = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "switch_RS"){switch_RS = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "switch_RS_mapping"){switch_RS_mapping = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "switch_yuk"){switch_yuk = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "order_y"){order_y = atoi(user_value[i].c_str());}

      //////////////////////////////////////
      //  weight-optimization parameters  //
      //////////////////////////////////////

      else if (user_variable[i] == "switch_n_events_opt"){switch_n_events_opt = atoi(user_value[i].c_str());}

      ////////////////////////////////////////////////////////////////////////////////////////////
      //  weight-optimization parameters - multi-channel Monte Carlo for final-state variables  //
      ////////////////////////////////////////////////////////////////////////////////////////////

      else if (user_variable[i] == "switch_MC"){switch_MC = atoi(user_value[i].c_str());}

      else if (user_variable[i] == "n_alpha_steps"){n_alpha_steps = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_alpha_events"){n_alpha_events = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_alpha_epc"){n_alpha_epc = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "MCweight_min"){MCweight_min = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "MCweight_limit_min"){MCweight_limit_min = atof(user_value[i].c_str());}
      else if (user_variable[i] == "MCweight_limit_max"){MCweight_limit_max = atof(user_value[i].c_str());}
      else if (user_variable[i] == "MCweight_in_contribution"){MCweight_in_contribution = user_value[i];}
      else if (user_variable[i] == "switch_use_alpha_after_IS"){switch_use_alpha_after_IS = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "switch_step_mode_grid"){switch_step_mode_grid = atoi(user_value[i].c_str());}

      // ??? needed ???

      else if (user_variable[i] == "MCweight_in_directory"){MCweight_in_directory = user_value[i];}
      // ??? needed ??? maybe location of grid to be read in...


      /////////////////////////////////////////////////////////////////////
      //  weight-optimization parameters - IS for final-state variables  //
      /////////////////////////////////////////////////////////////////////

      else if (user_variable[i] == "switch_IS_MC"){switch_IS_MC = atoi(user_value[i].c_str());}

      else if (user_variable[i] == "switch_IS_mode_phasespace"){switch_IS_mode_phasespace = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "MCweight_IS"){switch_IS_mode_phasespace = atoi(user_value[i].c_str());}

      else if (user_variable[i] == "n_IS_events"){n_IS_events = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_IS_events_factor"){n_IS_events_factor = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_IS_steps"){n_IS_steps = atoi(user_value[i].c_str());}
      
      else if (user_variable[i] == "n_IS_gridsize"){n_IS_gridsize = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_IS_gridsize_p"){n_IS_gridsize_p = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_IS_gridsize_f"){n_IS_gridsize_f = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_IS_gridsize_t_t"){n_IS_gridsize_t_t = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_IS_gridsize_t_phi"){n_IS_gridsize_t_phi = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_IS_gridsize_d_cth"){n_IS_gridsize_d_cth = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_IS_gridsize_d_phi"){n_IS_gridsize_d_phi = atoi(user_value[i].c_str());}

      /////////////////////////////////////////////////////////////////////////////////////////
      //  weight-optimization parameters - multi-channel Monte Carlo for CMS-energy mapping  //
      /////////////////////////////////////////////////////////////////////////////////////////

      else if (user_variable[i] == "switch_MC_tau"){switch_MC_tau = atoi(user_value[i].c_str());}

      //////////////////////////////////////////////////////////////////
      //  weight-optimization parameters - IS for CMS-energy mapping  //
      //////////////////////////////////////////////////////////////////

      else if (user_variable[i] == "switch_IS_tau"){switch_IS_tau = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_tau_steps"){n_tau_steps = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_tau_events"){n_tau_events = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_tau_bins"){n_tau_bins = atoi(user_value[i].c_str());}

      ////////////////////////////////////////////////////////////
      //  weight-optimization parameters - IS for x1x2 mapping  //
      ////////////////////////////////////////////////////////////

      else if (user_variable[i] == "switch_IS_x1x2"){switch_IS_x1x2 = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_x1x2_steps"){n_x1x2_steps = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_x1x2_events"){n_x1x2_events = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_x1x2_bins"){n_x1x2_bins = atoi(user_value[i].c_str());}

      ///////////////////////////////////////////////////////////////////////////////////////
      //  weight-optimization parameters - IS for z1 and z2 mappings (collinear emission)  //
      ///////////////////////////////////////////////////////////////////////////////////////

      else if (user_variable[i] == "switch_IS_z1z2"){switch_IS_z1z2 = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_z1z2_steps"){n_z1z2_steps = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_z1z2_events"){n_z1z2_events = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_z1z2_bins"){n_z1z2_bins = atoi(user_value[i].c_str());}

      //////////////////////////////////////////////////////////////////////////////////////////////////
      //  weight-optimization parameters - multi-channel Monte Carlo for x-values in dipole mappings  //
      //////////////////////////////////////////////////////////////////////////////////////////////////

      else if (user_variable[i] == "switch_MC_x_dipole"){switch_MC_x_dipole = atoi(user_value[i].c_str());}

      ////////////////////////////////////////////////////////////////////////////////
      //  weight-optimization parameters - IS for xy/zuv values in dipole mappings  //
      ////////////////////////////////////////////////////////////////////////////////

      else if (user_variable[i] == "n_IS_gridsize_xy"){n_IS_gridsize_xy = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_IS_gridsize_zuv"){n_IS_gridsize_zuv = atoi(user_value[i].c_str());}

      /////////////////////////////////////////
      //  phase-space generation parameters  //
      /////////////////////////////////////////

      else if (user_variable[i] == "nu"){nu = atof(user_value[i].c_str());}
      // removed ???

      else if (user_variable[i] == "nuxs"){nuxs = atof(user_value[i].c_str());}
      else if (user_variable[i] == "nuxt"){nuxt = atof(user_value[i].c_str());}
      else if (user_variable[i] == "exp_pdf"){exp_pdf = atof(user_value[i].c_str());}
      else if (user_variable[i] == "exp_pT"){exp_pT = atof(user_value[i].c_str());}
      // removed ???

      else if (user_variable[i] == "exp_y"){exp_y = atof(user_value[i].c_str());}
      // removed ???

      else if (user_variable[i] == "exp_ij_k_y"){exp_ij_k_y = atof(user_value[i].c_str());}
      else if (user_variable[i] == "exp_ij_k_z"){exp_ij_k_z = atof(user_value[i].c_str());}
      else if (user_variable[i] == "exp_ij_a_x"){exp_ij_a_x = atof(user_value[i].c_str());}
      else if (user_variable[i] == "exp_ij_a_z"){exp_ij_a_z = atof(user_value[i].c_str());}
      else if (user_variable[i] == "exp_ai_k_x"){exp_ai_k_x = atof(user_value[i].c_str());}
      else if (user_variable[i] == "exp_ai_k_u"){exp_ai_k_u = atof(user_value[i].c_str());}
      else if (user_variable[i] == "exp_ai_b_x"){exp_ai_b_x = atof(user_value[i].c_str());}
      else if (user_variable[i] == "exp_ai_b_v"){exp_ai_b_v = atof(user_value[i].c_str());}

      ////////////////////////////////////////
      //  technical integration parameters  //
      ////////////////////////////////////////

      else if (user_variable[i] == "mass0"){mass0 = atof(user_value[i].c_str());}
      else if (user_variable[i] == "map_technical_s"){map_technical_s = atof(user_value[i].c_str());}
      else if (user_variable[i] == "map_technical_t"){map_technical_t = atof(user_value[i].c_str());}
      else if (user_variable[i] == "map_technical_x"){map_technical_x = atof(user_value[i].c_str());}
      else if (user_variable[i] == "cut_technical"){cut_technical = atof(user_value[i].c_str());}

      //////////////////////////////
      //  integration parameters  //
      //////////////////////////////

      else if (user_variable[i] == "sigma_normalization_deviation"){sigma_normalization_deviation = atof(user_value[i].c_str());}
      else if (user_variable[i] == "sigma_LO_deviation"){sigma_normalization_deviation = atof(user_value[i].c_str());} // to be removed...
      else if (user_variable[i] == "n_events_max"){n_events_max = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_events_min"){n_events_min = atoi(user_value[i].c_str());}
      else if (user_variable[i] == "n_step"){n_step = atoi(user_value[i].c_str());}
    }

  }

  LHAPDFname = contribution_LHAPDFname[present_type_perturbative_order];
  LHAPDFsubset = contribution_LHAPDFsubset[present_type_perturbative_order];

  csi.readin_hadronic_process();
  csi.readin_subprocess();

  esi.name_fiducial_cut = name_fiducial_cut;
  /* // ???
  for (int i_f = 0; i_f < name_fiducial_cut.size(); i_f++){
    (esi.fiducial_cut).push_back(fiducialcut(name_fiducial_cut[i_f], esi));
  }
  */
  
  // Could be partially or completely shifted to TSV initialization !!!

  for (int i_s = 0; i_s < n_diff_set_TSV; i_s++){
    for (int j_s = 0; j_s < n_set_TSV; j_s++){
      if (name_diff_set_plus_TSV[i_s] == name_set_TSV[j_s]){no_diff_set_plus_TSV[i_s] = j_s;}
      if (name_diff_set_minus_TSV[i_s] == name_set_TSV[j_s]){no_diff_set_minus_TSV[i_s] = j_s;}
    }
    if (no_diff_set_plus_TSV[i_s] == -1){
      logger << LOG_FATAL << "name_diff_set_TSV[" << i_s << "] = " << name_diff_set_TSV[i_s] << "   has an invalid plus scale: name_diff_set_plus_TSV[" << i_s << "] = " << name_diff_set_plus_TSV[i_s] << endl;
      exit(1);
    }
    if (no_diff_set_minus_TSV[i_s] == -1){
      logger << LOG_FATAL << "name_diff_set_TSV[" << i_s << "] = " << name_diff_set_TSV[i_s] << "   has an invalid minus scale: name_diff_set_minus_TSV[" << i_s << "] = " << name_diff_set_minus_TSV[i_s] << endl;
      exit(1);
    }
  }

  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    logger << LOG_INFO << "name_set_TSV[" << i_s << "] = " << name_set_TSV[i_s] << endl;
  }
  logger.newLine(LOG_INFO);
  for (int i_s = 0; i_s < n_diff_set_TSV; i_s++){
    logger << LOG_INFO << "name_diff_set_TSV[" << i_s << "] = " << name_diff_set_TSV[i_s] << endl;
  }
  logger.newLine(LOG_INFO);
  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
    logger << LOG_INFO << "name_extended_set_TSV[" << i_s << "] = " << name_extended_set_TSV[i_s] << endl;
  }
  logger.newLine(LOG_INFO);

  // set switch_reference if not selected: could be shifted to observable_set !!!
  if (switch_reference == ""){
    if (switch_TSV && switch_CV){
      if (name_reference_TSV != ""){switch_reference = "TSV";}
      else {switch_reference = "CV";}
    }
    else if (switch_TSV){switch_reference = "TSV";}
    else if (switch_CV){switch_reference = "CV";}
  }
  // stop if switch_reference if not reasonably set (no results would be produced):
  if (switch_reference == ""){
    logger << LOG_FATAL << "No reasonable reference scale is set." << endl;
    exit(1);
  }

  // A number of variables need to be filled !!!
  if (switch_TSV){
    //    logger << LOG_INFO << "name_reference_TSV          = " << name_reference_TSV << endl;
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      //      logger << LOG_INFO << "name_extended_set_TSV[" << i_s << "] = " << name_extended_set_TSV[i_s] << endl;
      if (name_reference_TSV == name_extended_set_TSV[i_s]){
	no_reference_TSV = i_s;
	//	logger << LOG_INFO << "no_reference_TSV          = " << no_reference_TSV << endl;
	break;
      }
    }
    
    if (no_reference_TSV){
      logger << LOG_INFO << "default output from TSV:" << endl;
      logger << LOG_INFO << "name_reference_TSV          = " << name_reference_TSV << endl;
      logger << LOG_INFO << "no_reference_TSV            = " << no_reference_TSV << endl;
      logger << LOG_INFO << "no_scale_ren_reference_TSV  = " << no_scale_ren_reference_TSV << endl;
      logger << LOG_INFO << "no_scale_fact_reference_TSV = " << no_scale_fact_reference_TSV << endl;
      logger << LOG_INFO << "no_qTcut_reference_TSV      = " << no_qTcut_reference_TSV << endl;
    }
    /*
    else {
      logger << LOG_INFO << "No reference scale defined. Standard used for default output." << endl;
      logger << LOG_INFO << "No reference scale defined. Standard used for default output." << endl;
      // central value of first TSV scale is used as reference in summary:
      int x_s = 0;
      string temp_name_reference_TSV = osi_name_extended_set_TSV[x_s];
      oset.no_reference_TSV = x_s;
      oset.no_scale_ren_reference_TSV = oset.no_central_scale_ren_TSV[x_s];
      oset.no_scale_fact_reference_TSV = oset.no_central_scale_fact_TSV[x_s];

    }
*/
    logger.newLine(LOG_INFO);
  }
  //  logger << LOG_FATAL << "name_extended_set_TSV[" << i_s << "] = " << name_extended_set_TSV[i_s] << endl;

  /*
  csi.type_parton.resize(1);
  //  subprocess_readin(csi.process_class, csi.subprocess, csi.decay, csi.process_type, csi.n_particle, csi.type_parton);
  csi.readin_basic_subprocess();
  csi.determine_subprocess();
  */
  // !!! simplify via csi constructor !!!



  /*
  oset.initialization_object_generic();
  // shifted here

  oset.initialization_path(path_to_main);

  oset.initialization_unit(unit_calculation, unit_result, unit_distribution);
  // all oset references should be shifted elsewhere !!!


  map<string,int> code_particle;
  fill_code_particle(code_particle);
  map<int,string> name_particle;
  fill_name_particle(name_particle);
  */

  logger << LOG_INFO << setw(35) << "switch_TSV" << " = " << switch_TSV << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void inputparameter_set::parameter_readin_file(string filename, vector<string> & readin, bool essential, int retry=1) {
  Logger logger("inputparameter_set::parameter_readin_file");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  char LineBuffer[256];
  errno=0;

  for (int attempts=0; attempts<retry; attempts++) {
    ifstream in_new_parameter(filename.c_str());
    if (!in_new_parameter) {
      logger << LOG_WARN << filename << " could not be opened" << endl;
      logger << LOG_WARN << "Error code: " << errno << " (" << strerror(errno) << ")" << endl;
    } else {
      while (in_new_parameter.getline(LineBuffer, 256)){readin.push_back(LineBuffer);}
      in_new_parameter.close();
      break;
    }
    if (attempts<retry-1) {
      // workaround for NFS issues..
      logger << LOG_WARN << "retrying.." << endl;
      usleep(1000000);
    }
  }
  // assert out if the file_parameter.dat to be read in was essential
  assert(!essential || errno==0);
  logger << LOG_DEBUG << "filename = " << filename << "   " << readin.size() << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void inputparameter_set::parameter_readin(string & subprocess, vector<string> & readin){
  Logger logger("inputparameter_set::parameter_readin");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  
  logger << LOG_INFO << "subprocess = " << subprocess << endl;

  if (subprocess == ""){
    char LineBuffer[256];
    // default parameter file - should not be changed: the following ones should be used to overwrite any parameters set here!
    string filename;
    filename = "../setup/file_parameter.dat";
    //  cout << "filename = " << filename << endl;
    ifstream in_new_parameter_default(filename.c_str());
    while (in_new_parameter_default.getline(LineBuffer, 256)){readin.push_back(LineBuffer);}
    in_new_parameter_default.close();
    //  for (int i = 0; i < readin.size(); i++){cout << "default:    readin[" << setw(3) << "] = " << readin[i] << endl;}
    // general parameter file
    filename = "../file_parameter.dat";
    //  cout << "filename = " << filename << endl;
    ifstream in_new_parameter_gen(filename.c_str());
    while (in_new_parameter_gen.getline(LineBuffer, 256)){readin.push_back(LineBuffer);}
    in_new_parameter_gen.close();
    //  for (int i = 0; i < readin.size(); i++){cout << "gen:        readin[" << setw(3) << "] = " << readin[i] << endl;}
    // individual parameter file
    filename = "file_parameter.dat";
    //  cout << "filename = " << filename << endl;
    ifstream in_new_parameter_con(filename.c_str());
    while (in_new_parameter_con.getline(LineBuffer, 256)){readin.push_back(LineBuffer);}
    //  for (int i = 0; i < readin.size(); i++){cout << "individual: readin[" << setw(3) << "] = " << readin[i] << endl;}
    in_new_parameter_con.close();
  }
  else {
    string filename;
    // default parameter file - should not be changed: the following ones should be used overwrite any parameters set here!
    filename = "../setup/file_parameter.dat";
    parameter_readin_file(filename, readin, false, 1);
    
    filename = "../../../../setup/file_parameter.dat";
    parameter_readin_file(filename, readin, true, 3);
    
    // general parameter file
    filename = "../../../../file_parameter.dat";
    parameter_readin_file(filename, readin, true, 3);
    
    // perturbative_order.subtraction_method parameter file
    filename = "../../../file_parameter.dat";
    parameter_readin_file(filename, readin, false, 1);
    
    // coupling_order parameter file
    filename = "../../file_parameter.dat";
    parameter_readin_file(filename, readin, false, 1);
    
    // contribution.correction parameter file
    filename = "../file_parameter.dat";
    parameter_readin_file(filename, readin, false, 1);
    
    // individual rundirectory parameter file
    filename = "file_parameter.dat";
    parameter_readin_file(filename, readin, true, 3);
    
    // individual rundirectory subprocess parameter file
    filename = "log/file_parameter." + subprocess + ".dat";
    parameter_readin_file(filename, readin, false, 1);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void inputparameter_set::get_userinput_from_readin(vector<string> & user_variable, vector<string> & user_variable_additional, vector<string> & user_value, vector<string> & readin){
  Logger logger("get_userinput_from_readin");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string readindata;
  for (int i = 0; i < readin.size(); i++){
    logger << LOG_DEBUG << "readin[" << setw(4) << i << "].size() = " << readin[i].size() << endl;
    readindata = readin[i][0];
    if (readindata != "/" && readindata != "#" && readindata != "%"){
      int start = 0;
      user_variable.push_back("");
      user_variable_additional.push_back("");
      user_value.push_back("");
      for (int j = 0; j < readin[i].size(); j++){
	if (start == 0 || start == 1){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 0){}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    user_variable[user_variable.size() - 1].push_back(readin[i][j]);
	    if (start != 1){start = 1;}
	  }
	  else {start++;}
	}
	else if (start == 2){
	  if (readin[i][j] == '='){start = 5;}
	  else if ((readin[i][j] == ' ') || (readin[i][j] == char(9))){}
	  else if (start == 2){
	    start++;
	    j--;
	  }
	  else {
	    logger << LOG_ERROR << "Incorrect input in line " << i << endl;
	    user_variable.erase(user_variable.end(), user_variable.end());
	    user_variable_additional.erase(user_variable_additional.end(), user_variable_additional.end());
	    user_value.erase(user_value.end(), user_value.end());
	    break;
	  }
	}
	else if (start == 3){
	  if ((readin[i][j] != '=')){
	    user_variable_additional[user_variable_additional.size() - 1].push_back(readin[i][j]);
	  }
	  else {start++; j--;} // should be the same as shifting (start == 4) here !!!
	}
	else if (start == 4){
	  // additional: should be allowed to contain ' ' !!!
	  if (readin[i][j] == '='){
	    start = 5;
	    logger << LOG_DEBUG << "before: ---" << user_variable_additional[user_variable_additional.size() - 1] << "---" << endl;
	    for (int i_s = user_variable_additional[user_variable_additional.size() - 1].size() - 1; i_s > 0; i_s--){
	      if (user_variable_additional[user_variable_additional.size() - 1][i_s] == ' ' ||
		  user_variable_additional[user_variable_additional.size() - 1][i_s] == char(9)){
		user_variable_additional[user_variable_additional.size() - 1].erase(user_variable_additional[user_variable_additional.size() - 1].end() - 1, user_variable_additional[user_variable_additional.size() - 1].end());
	      }
	      else {break;}
	    }
	    logger << LOG_DEBUG << "after:  ---" << user_variable_additional[user_variable_additional.size() - 1] << "---" << endl;
	  }
	  //	  else if ((readin[i][j] == ' ') || (readin[i][j] == char(9))){}
	  else {
	    logger << LOG_ERROR << "Incorrect input in line " << i << endl;
	    user_variable.erase(user_variable.end(), user_variable.end());
	    user_variable_additional.erase(user_variable_additional.end(), user_variable_additional.end());
	    user_value.erase(user_value.end(), user_value.end());
	    break;
	  }
	}
	else if (start == 5 || start == 6){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 5){}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    user_value[user_value.size() - 1].push_back(readin[i][j]);
	    if (start != 6){start = 6;}
	  }
	  else {start++;}
	}

	else if (start == 5 || start == 6 || start == 7){
	  // start == 5: before beginning of user_value
	  // start == 6: user_value has started
	  // start == 7: 
	  logger << LOG_DEBUG_VERBOSE << "readin[" << i << "][" << j << "] = " << readin[i][j] << endl;
	  if (readin[i][j] == '%'){break;}
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 5){}
	  else if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start > 5){start = 7;}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    if (start == 7){user_value[user_value.size() - 1].push_back(' ');}
	    user_value[user_value.size() - 1].push_back(readin[i][j]);
	    start = 6;
	  }
	  else {start++;}
	}

	/*
	else if (start == 5 || start == 6){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 5){}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    user_value[user_value.size() - 1].push_back(readin[i][j]);
	    if (start != 6){start = 6;}
	  }
	  else {start++;}
	}
	*/

	else {break;}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



ostream & operator << (ostream & s_isi, const inputparameter_set & isi){
  int width = 40;

  s_isi << left << "Input from 'file_parameter.dat'" << endl;
  s_isi << left << "-------------------------------" << endl;
  s_isi << endl;

  s_isi << left << "unknown parameters" << endl;
  s_isi << left << "------------------" << endl;
  s_isi << left << setw(width) << "path_to_main" << " = " << isi.path_to_main << endl;
  s_isi << left << setw(width) << "process_number" << " = " << isi.process_number << endl;
  s_isi << left << setw(width) << "n_particle" << " = " << isi.n_particle << endl;
  s_isi << left << setw(width) << "test_output" << " = " << isi.test_output << endl;
  s_isi << left << setw(width) << "tau_0" << " = " << isi.tau_0 << endl;
  s_isi << left << setw(width) << "ckm_choice" << " = " << isi.ckm_choice << endl;
  s_isi << left << setw(width) << "select_contribution" << " = " << isi.select_contribution << endl;
  s_isi << left << setw(width) << "n_moments" << " = " << isi.n_moments << endl;
  s_isi << left << setw(width) << "dir_directory" << " = " << isi.dir_directory << endl;
  s_isi << left << setw(width) << "no_contribution" << " = " << isi.no_contribution << endl;
  s_isi << left << setw(width) << "impulsname.size()" << " = " << isi.impulsname.size() << endl;
  for (int i_i = 0; i_i < isi.impulsname.size(); i_i++){
    s_isi << left << setw(width - 5) << "impulsname" << "[" << setw(3) << i_i << "]" << " = " << isi.impulsname[i_i] << endl;
  }
  s_isi << endl;


  s_isi << left << "selection of process and contribution to be calculated (csi)" << endl;
  s_isi << left << "------------------------------------------------------------" << endl;

  s_isi << left << setw(width) << "csi.process_class" << " = " << isi.csi.process_class << endl;
  s_isi << left << setw(width) << "csi.decay.size()" << " = " << isi.csi.decay.size() << endl;
  for (int i_i = 0; i_i < isi.csi.decay.size(); i_i++){
    s_isi << left << setw(width - 5) << "csi.decay" << "[" << setw(3) << i_i << "]" << " = " << isi.csi.decay[i_i] << endl;
  }
  s_isi << left << setw(width) << "csi.type_perturbative_order" << " = " << isi.csi.type_perturbative_order << endl;
  s_isi << left << setw(width) << "csi.type_contribution" << " = " << isi.csi.type_contribution << endl;
  s_isi << left << setw(width) << "csi.type_correction" << " = " << isi.csi.type_correction << endl;
  s_isi << left << setw(width) << "csi.contribution_order_alpha_s" << " = " << isi.csi.contribution_order_alpha_s << endl;
  s_isi << left << setw(width) << "csi.contribution_order_alpha_e" << " = " << isi.csi.contribution_order_alpha_e << endl;
  s_isi << left << setw(width) << "csi.contribution_order_interference" << " = " << isi.csi.contribution_order_interference << endl;
  s_isi << left << setw(width) << "csi.subprocess" << " = " << isi.csi.subprocess << endl;
  s_isi << endl;


  s_isi << left << "switches to steer calculation of results" << endl;
  s_isi << left << "----------------------------------------" << endl;
  s_isi << left << setw(width) << "switch_result" << " = " << isi.switch_result << endl;
  s_isi << left << setw(width) << "switch_distribution" << " = " << isi.switch_distribution << endl;
  s_isi << left << setw(width) << "switch_moment" << " = " << isi.switch_moment << endl;
  s_isi << endl;


  s_isi << left << "switches to steer output of calculation" << endl;
  s_isi << left << "---------------------------------------" << endl;
  s_isi << left << setw(width) << "switch_output_execution" << " = " << isi.switch_output_execution << endl;
  s_isi << left << setw(width) << "switch_output_integration" << " = " << isi.switch_output_integration << endl;
  s_isi << left << setw(width) << "switch_output_maxevent" << " = " << isi.switch_output_maxevent << endl;
  s_isi << left << setw(width) << "switch_output_comparison" << " = " << isi.switch_output_comparison << endl;
  s_isi << left << setw(width) << "switch_output_gnuplot" << " = " << isi.switch_output_gnuplot << endl;
  s_isi << left << setw(width) << "switch_output_proceeding" << " = " << isi.switch_output_proceeding << endl;
  s_isi << left << setw(width) << "switch_output_weights" << " = " << isi.switch_output_weights << endl;
  s_isi << left << setw(width) << "switch_output_result" << " = " << isi.switch_output_result << endl;
  s_isi << left << setw(width) << "switch_output_moment" << " = " << isi.switch_output_moment << endl;
  s_isi << left << setw(width) << "switch_output_distribution" << " = " << isi.switch_output_distribution << endl;
  s_isi << left << setw(width) << "switch_output_cancellation_check" << " = " << isi.switch_output_cancellation_check << endl;
  s_isi << left << setw(width) << "switch_output_testpoint" << " = " << isi.switch_output_testpoint << endl;
  s_isi << left << setw(width) << "switch_output_cutinfo" << " = " << isi.switch_output_cutinfo << endl;
  s_isi << left << setw(width) << "switch_console_output_runtime" << " = " << isi.switch_console_output_runtime << endl;
  s_isi << left << setw(width) << "switch_console_output_tau_0" << " = " << isi.switch_console_output_tau_0 << endl;
  s_isi << left << setw(width) << "switch_console_output_techcut_RA" << " = " << isi.switch_console_output_techcut_RA << endl;
  s_isi << left << setw(width) << "switch_console_output_phasespace_issue" << " = " << isi.switch_console_output_phasespace_issue << endl;
  s_isi << left << setw(width) << "switch_console_output_ME2_issue" << " = " << isi.switch_console_output_ME2_issue << endl;
  s_isi << left << setw(width) << "switch_testcut" << " = " << isi.switch_testcut << endl;
  s_isi << endl;


  s_isi << left << "unit of calculation output / result output / distribution output" << endl;
  s_isi << left << "----------------------------------------------------------------" << endl;
  s_isi << left << setw(width) << "unit_calculation" << " = " << isi.unit_calculation << endl;
  s_isi << left << setw(width) << "unit_result" << " = " << isi.unit_result << endl;
  s_isi << left << setw(width) << "unit_distribution" << " = " << isi.unit_distribution << endl;
  s_isi << endl;


  s_isi << left << "beam parameters" << endl;
  s_isi << left << "---------------" << endl;
  s_isi << left << setw(width) << "E" << " = " << isi.E << endl;
  s_isi << left << setw(width) << "coll_choice" << " = " << isi.coll_choice << endl;
  //  s_isi << left << setw(width) << "pdf_set" << " = " << isi.pdf_set << endl;
  //  s_isi << left << setw(width) << "pdf_content_modify" << " = " << isi.pdf_content_modify << endl;
  s_isi << left << setw(width) << "pdf_selection.size()" << " = " << isi.pdf_selection.size() << endl;
  for (int i_i = 0; i_i < isi.pdf_selection.size(); i_i++){
    s_isi << left << setw(width - 5) << "pdf_selection" << "[" << setw(3) << i_i << "]" << " = " << isi.pdf_selection[i_i] << endl;
  }
  s_isi << left << setw(width) << "pdf_disable.size()" << " = " << isi.pdf_disable.size() << endl;
  for (int i_i = 0; i_i < isi.pdf_disable.size(); i_i++){
    s_isi << left << setw(width - 5) << "pdf_disable" << "[" << setw(3) << i_i << "]" << " = " << isi.pdf_disable[i_i] << endl;
  }
  s_isi << left << setw(width) << "LHAPDFname" << " = " << isi.LHAPDFname << endl;
  s_isi << left << setw(width) << "LHAPDFsubset" << " = " << isi.LHAPDFsubset << endl;
  s_isi << left << setw(width) << "contribution_LHAPDFname.size()" << " = " << isi.contribution_LHAPDFname.size() << endl;
  for (int i_i = 0; i_i < isi.contribution_LHAPDFname.size(); i_i++){
    s_isi << left << setw(width - 5) << "contribution_LHAPDFname" << "[" << setw(3) << i_i << "]" << " = " << isi.contribution_LHAPDFname[i_i] << endl;
  }
  s_isi << left << setw(width) << "contribution_LHAPDFsubset.size()" << " = " << isi.contribution_LHAPDFsubset.size() << endl;
  for (int i_i = 0; i_i < isi.contribution_LHAPDFsubset.size(); i_i++){
    s_isi << left << setw(width - 5) << "contribution_LHAPDFsubset" << "[" << setw(3) << i_i << "]" << " = " << isi.contribution_LHAPDFsubset[i_i] << endl;
  }
  s_isi << endl;


  s_isi << left << "jet and jet-algorithm parameters" << endl;
  s_isi << left << "--------------------------------" << endl;
  s_isi << left << setw(width) << "jet_algorithm" << " = " << isi.jet_algorithm << endl;
  s_isi << left << setw(width) << "jet_R_definition" << " = " << isi.jet_R_definition << endl;
  s_isi << left << setw(width) << "jet_R" << " = " << isi.jet_R << endl;
  s_isi << left << setw(width) << "parton_y_max" << " = " << isi.parton_y_max << endl;
  s_isi << left << setw(width) << "parton_eta_max" << " = " << isi.parton_eta_max << endl;
  s_isi << left << setw(width) << "jet_algorithm_selection.size()" << " = " << isi.jet_algorithm_selection.size() << endl;
  for (int i_i = 0; i_i < isi.jet_algorithm_selection.size(); i_i++){
    s_isi << left << setw(width - 5) << "jet_algorithm_selection" << "[" << setw(3) << i_i << "]" << " = " << isi.jet_algorithm_selection[i_i] << endl;
  }
  s_isi << left << setw(width) << "jet_algorithm_disable.size()" << " = " << isi.jet_algorithm_disable.size() << endl;
  for (int i_i = 0; i_i < isi.jet_algorithm_disable.size(); i_i++){
    s_isi << left << setw(width - 5) << "jet_algorithm_disable" << "[" << setw(3) << i_i << "]" << " = " << isi.jet_algorithm_disable[i_i] << endl;
  }
  s_isi << left << setw(width) << "N_f" << " = " << isi.N_f << endl;
  s_isi << left << setw(width) << "N_f_active" << " = " << isi.N_f_active << endl;
  s_isi << left << setw(width) << "N_quarks" << " = " << isi.N_quarks << endl;
  s_isi << left << setw(width) << "N_nondecoupled" << " = " << isi.N_nondecoupled << endl;
  s_isi << endl;


  s_isi << left << "photon-recombination parameters" << endl;
  s_isi << left << "-------------------------------" << endl;
  s_isi << left << setw(width) << "photon_recombination" << " = " << isi.photon_recombination << endl;
  s_isi << left << setw(width) << "photon_R_definition" << " = " << isi.photon_R_definition << endl;
  s_isi << left << setw(width) << "photon_R" << " = " << isi.photon_R << endl;
  s_isi << left << setw(width) << "photon_E_threshold_ratio" << " = " << isi.photon_E_threshold_ratio << endl;
  s_isi << left << setw(width) << "photon_jet_algorithm" << " = " << isi.photon_jet_algorithm << endl;
  s_isi << left << setw(width) << "photon_recombination_selection.size()" << " = " << isi.photon_recombination_selection.size() << endl;
  for (int i_i = 0; i_i < isi.photon_recombination_selection.size(); i_i++){
    s_isi << left << setw(width - 5) << "photon_recombination_selection" << "[" << setw(3) << i_i << "]" << " = " << isi.photon_recombination_selection[i_i] << endl;
  }
  s_isi << left << setw(width) << "photon_recombination_disable.size()" << " = " << isi.photon_recombination_disable.size() << endl;
  for (int i_i = 0; i_i < isi.photon_recombination_disable.size(); i_i++){
    s_isi << left << setw(width - 5) << "photon_recombination_disable" << "[" << setw(3) << i_i << "]" << " = " << isi.photon_recombination_disable[i_i] << endl;
  }
  s_isi << left << setw(width) << "photon_photon_recombination" << " = " << isi.photon_photon_recombination << endl;
  s_isi << left << setw(width) << "photon_photon_recombination_R" << " = " << isi.photon_photon_recombination_R << endl;
  s_isi << endl;

  s_isi << left << "photon-isolation parameters" << endl;
  s_isi << left << "---------------------------" << endl;
  s_isi << left << setw(width) << "frixione_isolation" << " = " << isi.frixione_isolation << endl;
  s_isi << left << setw(width) << "frixione_n" << " = " << isi.frixione_n << endl;
  if (isi.frixione_isolation == 2){
    s_isi << left << setw(width) << "frixione_fixed_ET_max" << " = " << isi.frixione_fixed_ET_max << endl;
  }
  else {
    s_isi << left << setw(width) << "frixione_epsilon" << " = " << isi.frixione_epsilon << endl;
  }
  s_isi << left << setw(width) << "frixione_delta_0" << " = " << isi.frixione_delta_0 << endl;
  s_isi << left << setw(width) << "frixione_jet_removal" << " = " << isi.frixione_jet_removal << endl;
  s_isi << endl;

  s_isi << left << "qT-subtraction parameters" << endl;
  s_isi << left << "-------------------------" << endl;
  s_isi << left << setw(width) << "switch_qTcut" << " = " << isi.switch_qTcut << endl;
  s_isi << left << setw(width) << "n_qTcut" << " = " << isi.n_qTcut << endl;
  s_isi << left << setw(width) << "min_qTcut" << " = " << isi.min_qTcut << endl;
  s_isi << left << setw(width) << "step_qTcut" << " = " << isi.step_qTcut << endl;
  s_isi << left << setw(width) << "binning_qTcut" << " = " << isi.binning_qTcut << endl;
  s_isi << left << setw(width) << "selection_qTcut_integration" << " = " << isi.selection_qTcut_integration << endl;
  s_isi << left << setw(width) << "selection_no_qTcut_integration" << " = " << isi.selection_no_qTcut_integration << endl;
  s_isi << left << setw(width) << "selection_qTcut_result" << " = " << isi.selection_qTcut_result << endl;
  s_isi << left << setw(width) << "selection_no_qTcut_result" << " = " << isi.selection_no_qTcut_result << endl;
  s_isi << left << setw(width) << "selection_qTcut_distribution" << " = " << isi.selection_qTcut_distribution << endl;
  s_isi << left << setw(width) << "selection_no_qTcut_distribution" << " = " << isi.selection_no_qTcut_distribution << endl;
  s_isi << endl;

  s_isi << left << "N-jettiness-subtraction parameters" << endl;
  s_isi << left << "----------------------------------" << endl;
  s_isi << left << setw(width) << "switch_NJcut" << " = " << isi.switch_NJcut << endl;
  s_isi << left << setw(width) << "switch_NJcut_axes" << " = " << isi.switch_NJcut_axes << endl;
  s_isi << left << setw(width) << "switch_NJcut_axes_energy" << " = " << isi.switch_NJcut_axes_energy << endl;
  s_isi << left << setw(width) << "switch_NJcut_measure" << " = " << isi.switch_NJcut_measure << endl;
  s_isi << endl;


  s_isi << left << "technical switches for selected contributions" << endl;
  s_isi << left << "---------------------------------------------" << endl;
  s_isi << left << setw(width) << "switch_KP" << " = " << isi.switch_KP << endl;
  s_isi << left << setw(width) << "switch_VI" << " = " << isi.switch_VI << endl;
  s_isi << left << setw(width) << "switch_VI_bosonic_fermionic" << " = " << isi.switch_VI_bosonic_fermionic << endl;
  s_isi << left << setw(width) << "switch_polenorm" << " = " << isi.switch_polenorm << endl;
  s_isi << left << setw(width) << "switch_H1gg" << " = " << isi.switch_H1gg << endl;
  s_isi << left << setw(width) << "switch_H2" << " = " << isi.switch_H2 << endl;
  s_isi << left << setw(width) << "switch_old_qT_version" << " = " << isi.switch_old_qT_version << endl;
  s_isi << left << setw(width) << "switch_RS" << " = " << isi.switch_RS << endl;
  s_isi << left << setw(width) << "switch_RS_mapping" << " = " << isi.switch_RS_mapping << endl;
  s_isi << left << setw(width) << "switch_CM" << " = " << isi.switch_CM << endl;
  s_isi << left << setw(width) << "switch_OL" << " = " << isi.switch_OL << endl;
  s_isi << endl;

  s_isi << left << "scale-choice parameters" << endl;
  s_isi << left << "-----------------------" << endl;
  s_isi << left << setw(width) << "scale_fact" << " = " << isi.scale_fact << endl;
  s_isi << left << setw(width) << "scale_ren" << " = " << isi.scale_ren << endl;
  s_isi << left << setw(width) << "dynamic_scale" << " = " << isi.dynamic_scale << endl;
  s_isi << left << setw(width) << "prefactor_reference" << " = " << isi.prefactor_reference << endl;
  s_isi << endl;

  s_isi << left << "scale variation parameters - CV" << endl;
  s_isi << left << "-------------------------------" << endl;
  s_isi << left << setw(width) << "switch_CV" << " = " << isi.switch_CV << endl;
  s_isi << left << setw(width) << "n_scales_CV" << " = " << isi.n_scales_CV << endl;
  s_isi << left << setw(width) << "dynamic_scale_CV" << " = " << isi.dynamic_scale_CV << endl;
  s_isi << left << setw(width) << "variation_mu_ren_CV" << " = " << isi.variation_mu_ren_CV << endl;
  s_isi << left << setw(width) << "variation_mu_fact_CV" << " = " << isi.variation_mu_fact_CV << endl;
  s_isi << left << setw(width) << "variation_factor_CV" << " = " << isi.variation_factor_CV << endl;
  s_isi << left << setw(width) << "central_scale_CV" << " = " << isi.central_scale_CV << endl;
  s_isi << left << setw(width) << "prefactor_CV" << " = " << isi.prefactor_CV << endl;
  s_isi << endl;

  s_isi << left << "weight-optimization parameters" << endl;
  s_isi << left << "------------------------------" << endl;
  s_isi << left << setw(width) << "switch_n_events_opt" << " = " << isi.switch_n_events_opt << endl;
  s_isi << endl;

  s_isi << left << "weight-optimization parameters - multi-channel Monte Carlo for final-state variables" << endl;
  s_isi << left << "------------------------------------------------------------------------------------" << endl;
  s_isi << left << setw(width) << "switch_MC" << " = " << isi.switch_MC << endl;
  s_isi << left << setw(width) << "n_alpha_steps" << " = " << isi.n_alpha_steps << endl;
  s_isi << left << setw(width) << "n_alpha_events" << " = " << isi.n_alpha_events << endl;
  s_isi << left << setw(width) << "n_alpha_epc" << " = " << isi.n_alpha_epc << endl;
  s_isi << left << setw(width) << "MCweight_min" << " = " << isi.MCweight_min << endl;
  s_isi << left << setw(width) << "MCweight_limit_min" << " = " << isi.MCweight_limit_min << endl;
  s_isi << left << setw(width) << "MCweight_limit_max" << " = " << isi.MCweight_limit_max << endl;
  s_isi << left << setw(width) << "MCweight_in_contribution" << " = " << isi.MCweight_in_contribution << endl;
  s_isi << left << setw(width) << "MCweight_in_directory" << " = " << isi.MCweight_in_directory << endl;
  s_isi << endl;

  s_isi << left << "weight-optimization parameters - IS for final-state variables" << endl;
  s_isi << left << "-------------------------------------------------------------" << endl;
  s_isi << left << setw(width) << "switch_IS_MC" << " = " << isi.switch_IS_MC << endl;
  s_isi << left << setw(width) << "switch_IS_mode_phasespace" << " = " << isi.switch_IS_mode_phasespace << endl;
  s_isi << left << setw(width) << "n_IS_events" << " = " << isi.n_IS_events << endl;
  s_isi << left << setw(width) << "n_IS_events_factor" << " = " << isi.n_IS_events_factor << endl;
  s_isi << left << setw(width) << "n_IS_steps" << " = " << isi.n_IS_steps << endl;
  s_isi << left << setw(width) << "n_IS_gridsize" << " = " << isi.n_IS_gridsize << endl;
  s_isi << left << setw(width) << "n_IS_gridsize_p" << " = " << isi.n_IS_gridsize_p << endl;
  s_isi << left << setw(width) << "n_IS_gridsize_f" << " = " << isi.n_IS_gridsize_f << endl;
  s_isi << left << setw(width) << "n_IS_gridsize_t_t" << " = " << isi.n_IS_gridsize_t_t << endl;
  s_isi << left << setw(width) << "n_IS_gridsize_t_phi" << " = " << isi.n_IS_gridsize_t_phi << endl;
  s_isi << left << setw(width) << "n_IS_gridsize_d_cth" << " = " << isi.n_IS_gridsize_d_cth << endl;
  s_isi << left << setw(width) << "n_IS_gridsize_d_phi" << " = " << isi.n_IS_gridsize_d_phi << endl;
  s_isi << endl;

  s_isi << left << "weight-optimization parameters - multi-channel Monte Carlo for CMS-energy mapping" << endl;
  s_isi << left << "---------------------------------------------------------------------------------" << endl;
  s_isi << left << setw(width) << "switch_MC_tau" << " = " << isi.switch_MC_tau << endl;
  s_isi << endl;

  s_isi << left << "weight-optimization parameters - IS for CMS-energy mapping" << endl;
  s_isi << left << "----------------------------------------------------------" << endl;
  s_isi << left << setw(width) << "switch_IS_tau" << " = " << isi.switch_IS_tau << endl;
  s_isi << left << setw(width) << "n_tau_steps" << " = " << isi.n_tau_steps << endl;
  s_isi << left << setw(width) << "n_tau_events" << " = " << isi.n_tau_events << endl;
  s_isi << left << setw(width) << "n_tau_bins" << " = " << isi.n_tau_bins << endl;
  s_isi << endl;

  s_isi << left << "weight-optimization parameters - IS for x1x2 mapping" << endl;
  s_isi << left << "----------------------------------------------------" << endl;
  s_isi << left << setw(width) << "switch_IS_x1x2" << " = " << isi.switch_IS_x1x2 << endl;
  s_isi << left << setw(width) << "n_x1x2_steps" << " = " << isi.n_x1x2_steps << endl;
  s_isi << left << setw(width) << "n_x1x2_events" << " = " << isi.n_x1x2_events << endl;
  s_isi << left << setw(width) << "n_x1x2_bins" << " = " << isi.n_x1x2_bins << endl;
  s_isi << endl;

  s_isi << left << "weight-optimization parameters - IS for z1 and z2 mappings (collinear emission)" << endl;
  s_isi << left << "-------------------------------------------------------------------------------" << endl;
  s_isi << left << setw(width) << "switch_IS_z1z2" << " = " << isi.switch_IS_z1z2 << endl;
  s_isi << left << setw(width) << "n_z1z2_steps" << " = " << isi.n_z1z2_steps << endl;
  s_isi << left << setw(width) << "n_z1z2_events" << " = " << isi.n_z1z2_events << endl;
  s_isi << left << setw(width) << "n_z1z2_bins" << " = " << isi.n_z1z2_bins << endl;
  s_isi << endl;

  s_isi << left << "weight-optimization parameters - multi-channel Monte Carlo for x-values in dipole mappings" << endl;
  s_isi << left << "------------------------------------------------------------------------------------------" << endl;
  s_isi << left << setw(width) << "switch_MC_x_dipole" << " = " << isi.switch_MC_x_dipole << endl;
  s_isi << endl;

  s_isi << left << "" << endl;
  s_isi << left << "------------------------------------------------------------" << endl;
  s_isi << left << setw(width) << "n_IS_gridsize_xy" << " = " << isi.n_IS_gridsize_xy << endl;
  s_isi << left << setw(width) << "n_IS_gridsize_zuv" << " = " << isi.n_IS_gridsize_zuv << endl;
  s_isi << endl;

  s_isi << left << "phase-space generation parameters" << endl;
  s_isi << left << "---------------------------------" << endl;
  s_isi << left << setw(width) << "nu" << " = " << isi.nu << endl;
  s_isi << left << setw(width) << "nuxs" << " = " << isi.nuxs << endl;
  s_isi << left << setw(width) << "nuxt" << " = " << isi.nuxt << endl;
  s_isi << left << setw(width) << "exp_pdf" << " = " << isi.exp_pdf << endl;
  s_isi << left << setw(width) << "exp_pT" << " = " << isi.exp_pT << endl;
  s_isi << left << setw(width) << "exp_y" << " = " << isi.exp_y << endl;
  s_isi << left << setw(width) << "exp_ij_k_y" << " = " << isi.exp_ij_k_y << endl;
  s_isi << left << setw(width) << "exp_ij_k_z" << " = " << isi.exp_ij_k_z << endl;
  s_isi << left << setw(width) << "exp_ij_a_x" << " = " << isi.exp_ij_a_x << endl;
  s_isi << left << setw(width) << "exp_ij_a_z" << " = " << isi.exp_ij_a_z << endl;
  s_isi << left << setw(width) << "exp_ai_k_x" << " = " << isi.exp_ai_k_x << endl;
  s_isi << left << setw(width) << "exp_ai_k_u" << " = " << isi.exp_ai_k_u << endl;
  s_isi << left << setw(width) << "exp_ai_b_x" << " = " << isi.exp_ai_b_x << endl;
  s_isi << left << setw(width) << "exp_ai_b_v" << " = " << isi.exp_ai_b_v << endl;
  s_isi << endl;

  s_isi << left << "technical integration parameters" << endl;
  s_isi << left << "--------------------------------" << endl;
  s_isi << left << setw(width) << "mass0" << " = " << isi.mass0 << endl;
  s_isi << left << setw(width) << "map_technical_s" << " = " << isi.map_technical_s << endl;
  s_isi << left << setw(width) << "map_technical_t" << " = " << isi.map_technical_t << endl;
  s_isi << left << setw(width) << "map_technical_x" << " = " << isi.map_technical_x << endl;
  s_isi << left << setw(width) << "cut_technical" << " = " << isi.cut_technical << endl;
  s_isi << endl;

  s_isi << left << "------------------------------------------------------------" << endl;
  s_isi << left << setw(width) << "sigma_normalization" << " = " << isi.sigma_normalization << endl;
  s_isi << left << setw(width) << "zwahl" << " = " << isi.zwahl << endl;
  s_isi << left << setw(width) << "sigma_normalization_deviation" << " = " << isi.sigma_normalization_deviation << endl;
  s_isi << left << setw(width) << "n_events_max" << " = " << isi.n_events_max << endl;
  s_isi << left << setw(width) << "n_events_min" << " = " << isi.n_events_min << endl;
  s_isi << left << setw(width) << "n_step" << " = " << isi.n_step << endl;
  s_isi << endl;

  // more reasonable organization of that output: !!!
  s_isi << left << "parameters directly forwarded to OpenLoops - after the default settings are done" << endl;
  s_isi << left << "--------------------------------------------------------------------------------" << endl;
  s_isi << left << setw(width) << "OL_parameter.size()" << " = " << isi.OL_parameter.size() << endl;
  for (int i_i = 0; i_i < isi.OL_parameter.size(); i_i++){
    s_isi << left << setw(width - 5) << "OL_parameter" << "[" << setw(3) << i_i << "]" << " = " << isi.OL_parameter[i_i] << endl;
  }
  s_isi << left << setw(width) << "OL_value.size()" << " = " << isi.OL_value.size() << endl;
  for (int i_i = 0; i_i < isi.OL_value.size(); i_i++){
    s_isi << left << setw(width - 5) << "OL_value" << "[" << setw(3) << i_i << "]" << " = " << isi.OL_value[i_i] << endl;
  }
  s_isi << endl;

  s_isi << left << "scale variation parameters - TSV" << endl;
  s_isi << left << "--------------------------------" << endl;
  s_isi << left << setw(width) << "switch_TSV" << " = " << isi.switch_TSV << endl;
  s_isi << left << setw(width) << "name_set_TSV.size()" << " = " << isi.name_set_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.name_set_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "name_set_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.name_set_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "n_set_TSV" << " = " << isi.n_set_TSV << endl;
  s_isi << left << setw(width) << "central_scale_TSV.size()" << " = " << isi.central_scale_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.central_scale_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "central_scale_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.central_scale_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "central_scale_ren_TSV.size()" << " = " << isi.central_scale_ren_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.central_scale_ren_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "central_scale_ren_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.central_scale_ren_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "central_scale_fact_TSV.size()" << " = " << isi.central_scale_fact_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.central_scale_fact_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "central_scale_fact_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.central_scale_fact_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "relative_central_scale_TSV.size()" << " = " << isi.relative_central_scale_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.relative_central_scale_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "relative_central_scale_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.relative_central_scale_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "relative_central_scale_ren_TSV.size()" << " = " << isi.relative_central_scale_ren_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.relative_central_scale_ren_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "relative_central_scale_ren_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.relative_central_scale_ren_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "relative_central_scale_fact_TSV.size()" << " = " << isi.relative_central_scale_fact_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.relative_central_scale_fact_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "relative_central_scale_fact_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.relative_central_scale_fact_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "n_scale_TSV.size()" << " = " << isi.n_scale_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.n_scale_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "n_scale_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.n_scale_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "n_scale_ren_TSV.size()" << " = " << isi.n_scale_ren_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.n_scale_ren_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "n_scale_ren_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.n_scale_ren_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "n_scale_fact_TSV.size()" << " = " << isi.n_scale_fact_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.n_scale_fact_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "n_scale_fact_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.n_scale_fact_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "factor_scale_TSV.size()" << " = " << isi.factor_scale_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.factor_scale_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "factor_scale_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.factor_scale_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "factor_scale_ren_TSV.size()" << " = " << isi.factor_scale_ren_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.factor_scale_ren_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "factor_scale_ren_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.factor_scale_ren_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "factor_scale_fact_TSV.size()" << " = " << isi.factor_scale_fact_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.factor_scale_fact_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "factor_scale_fact_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.factor_scale_fact_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "dynamic_scale_TSV.size()" << " = " << isi.dynamic_scale_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.dynamic_scale_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "dynamic_scale_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.dynamic_scale_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "dynamic_scale_ren_TSV.size()" << " = " << isi.dynamic_scale_ren_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.dynamic_scale_ren_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "dynamic_scale_ren_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.dynamic_scale_ren_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "dynamic_scale_fact_TSV.size()" << " = " << isi.dynamic_scale_fact_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.dynamic_scale_fact_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "dynamic_scale_fact_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.dynamic_scale_fact_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "min_qTcut_TSV.size()" << " = " << isi.min_qTcut_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.min_qTcut_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "min_qTcut_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.min_qTcut_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "max_qTcut_TSV.size()" << " = " << isi.max_qTcut_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.max_qTcut_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "max_qTcut_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.max_qTcut_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "switch_distribution_TSV.size()" << " = " << isi.switch_distribution_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.switch_distribution_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "switch_distribution_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.switch_distribution_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "min_qTcut_distribution_TSV.size()" << " = " << isi.min_qTcut_distribution_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.min_qTcut_distribution_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "min_qTcut_distribution_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.min_qTcut_distribution_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "max_qTcut_distribution_TSV.size()" << " = " << isi.max_qTcut_distribution_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.max_qTcut_distribution_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "max_qTcut_distribution_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.max_qTcut_distribution_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "switch_moment_TSV.size()" << " = " << isi.switch_moment_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.switch_moment_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "switch_moment_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.switch_moment_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "max_n_integrand_TSV.size()" << " = " << isi.max_n_integrand_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.max_n_integrand_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "max_n_integrand_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.max_n_integrand_TSV[i_i] << endl;
  }
  s_isi << left << setw(width) << "no_central_scale_ren_TSV.size()" << " = " << isi.no_central_scale_ren_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.no_central_scale_ren_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "no_central_scale_ren_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.no_central_scale_ren_TSV[i_i] << endl;
  }
  s_isi << endl;
  s_isi << left << setw(width) << "no_central_scale_fact_TSV.size()" << " = " << isi.no_central_scale_fact_TSV.size() << endl;
  for (int i_i = 0; i_i < isi.no_central_scale_fact_TSV.size(); i_i++){
    s_isi << left << setw(width - 5) << "no_central_scale_fact_TSV" << "[" << setw(3) << i_i << "]" << " = " << isi.no_central_scale_fact_TSV[i_i] << endl;
  }
  s_isi << endl;

  s_isi << left << "user-defined parameters - user" << endl;
  s_isi << left << "------------------------------" << endl;
  for (int i_i = 0; i_i < isi.user.switch_name.size(); i_i++){
    s_isi << left << setw(width - 5) << "user.switch_name" << "[" << setw(3) << i_i << "]" << " = " << isi.user.switch_name[i_i] << " -> " << isi.user.switch_value[i_i] << "   " << endl;//isi.user.switch_map[(isi.user.switch_name)[i_i]] << endl;
    /*
    string temp_s = isi.user.switch_name[i_i];
    map<string, int> switch_map = isi.user.switch_map;
    int temp = switch_map[temp_s];
    int temp2 = (isi.user.switch_map)[temp_s];
    s_isi << temp << endl;
    s_isi << temp2 << endl;
    */    
  }
  s_isi << endl;
  for (int i_i = 0; i_i < isi.user.cut_name.size(); i_i++){
    s_isi << left << setw(width - 5) << "user.cut_name" << "[" << setw(3) << i_i << "]" << " = " << isi.user.cut_name[i_i] << " -> " << isi.user.cut_value[i_i] << "   " << endl;//isi.user.cut_map[isi.user.cut_name[i_i]] << endl;
  }
  s_isi << endl;
  for (int i_i = 0; i_i < isi.user.int_name.size(); i_i++){
    s_isi << left << setw(width - 5) << "user.int_name" << "[" << setw(3) << i_i << "]" << " = " << isi.user.int_name[i_i] << " -> " << isi.user.int_value[i_i] << "   " << endl;//isi.user.int_map[isi.user.int_name[i_i]] << endl;
  }
  s_isi << endl;
  for (int i_i = 0; i_i < isi.user.double_name.size(); i_i++){
    s_isi << left << setw(width - 5) << "user.double_name" << "[" << setw(3) << i_i << "]" << " = " << isi.user.double_name[i_i] << " -> " << isi.user.double_value[i_i] << endl;//"   " << isi.user.double_map[isi.user.double_name[i_i]] << endl;
  }
  s_isi << endl;
  for (int i_i = 0; i_i < isi.user.string_name.size(); i_i++){
    s_isi << left << setw(width - 5) << "user.string_name" << "[" << setw(3) << i_i << "]" << " = " << isi.user.string_name[i_i] << " -> " << isi.user.string_value[i_i] << endl;//"   " << isi.user.string_map[isi.user.string_name[i_i]] << endl;
  }
  s_isi << endl;


  return s_isi;

  /*
  s_isi << left << "" << endl;
  s_isi << left << "------------------------------------------------------------" << endl;
  s_isi << endl;

  s_isi << left << setw(width) << "" << " = " << isi. << endl;

  s_isi << left << setw(width) << "xxx.size()" << " = " << isi.xxx.size() << endl;
  for (int i_i = 0; i_i < isi.xxx.size(); i_i++){
    s_isi << left << setw(width - 5) << "xxx" << "[" << setw(3) << i_i << "]" << " = " << isi.xxx[i_i] << endl;
  }
  */
}
