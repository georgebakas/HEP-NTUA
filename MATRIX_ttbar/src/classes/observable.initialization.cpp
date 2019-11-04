#include "classes.cxx"
#include "definitions.phasespace.set.cxx"
#include "more.h"

// Arguments should be redundant !!!
void observable_set::initialization_complete(inputparameter_set & isi, phasespace_set & psi){
  Logger logger("observable_set::initialization_complete");
  logger << LOG_DEBUG << "started" << endl;
  logger << LOG_DEBUG << "csi->type_contribution "<<csi->type_contribution << endl;

  if (csi->type_contribution == "born" ||
      csi->type_contribution == "RT" ||
      csi->type_contribution == "L2I" ||
      csi->type_contribution == "loop" ||
      csi->type_contribution == "L2RT"){
    initialization_tree();
    initialization_runtime_partonlevel();
    perform_selection_content_pdf();
  }

  if (csi->type_contribution == "RA" ||
      csi->type_contribution == "RRA" ||
      csi->type_contribution == "L2RA"){
    initialization_RA(xmunich->dipole);
    initialization_runtime_partonlevel();
    perform_selection_content_pdf();
  }
  if (csi->type_contribution == "CA" ||
      csi->type_contribution == "RCA" ||
      csi->type_contribution == "L2CA"){
    CA_collinear = &xmunich->collinear;
    if (csi->type_correction == "QCD"){determine_collinear_QCD(psi);}
    if (csi->type_correction == "QEW"){determine_collinear_QEW(psi);}
    initialization_CA();
    initialization_runtime_partonlevel();
    perform_selection_content_pdf_collinear();
  }
  if (csi->type_contribution == "VA" ||
      csi->type_contribution == "RVA" ||
      csi->type_contribution == "L2VA"){
    VA_ioperator = &xmunich->ioperator;
    if (csi->type_correction == "QCD" || csi->type_correction == "MIX"){determine_ioperator_QCD(psi);}
    if (csi->type_correction == "QEW" || csi->type_correction == "MIX"){determine_ioperator_QEW(psi);}
    initialization_VA();
    initialization_runtime_partonlevel();
    perform_selection_content_pdf();
  }
  if (csi->type_contribution == "CT" ||
      csi->type_contribution == "VT"){
    determine_CX_QCD(psi, 1);
    initialization_QT();
    initialization_runtime_partonlevel();
  }
  if (csi->type_contribution == "VT"){
    //  workaround for HH checks: actually not needed in VT contribution
    //  determine_ioperator_QCD(psi);
  }
  if (csi->type_contribution == "CT2" ||
      csi->type_contribution == "VT2"){
    determine_CX_QCD(psi, 2);
    initialization_QT();
    initialization_runtime_partonlevel();
    // Some version of perform_selection_content_pdf(); ??? !!!
  }
  
  initialization_TSV();
  initialization_CV();
  
  if (csi->type_contribution == "CT" ||
      csi->type_contribution == "VT" ||
      csi->type_contribution == "CT2" ||
      csi->type_contribution == "VT2"){
    initialization_CX_ncollinear();
  }
  
  initialization_distribution();
  
  if (csi->type_contribution == "VT"){
    //  workaround for HH checks: actually not needed in VT contribution
    //  osi.initialization_specific_VA();
  }
  if (csi->type_contribution == "CT" ||
      csi->type_contribution == "VT" ||
      csi->type_contribution == "CT2" ||
      csi->type_contribution == "VT2"){
    initialization_QT_coefficients();
  }
 
  initialization_LHAPDF();
  if (switch_OL){initialization_OpenLoops_process(psi);}

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_path(inputparameter_set & isi){
  static Logger logger("observable_set::initialization_path");
  logger << LOG_DEBUG << "called" << endl;

  path_to_main = isi.path_to_main;

  logger << LOG_DEBUG << "finished" << endl;
}



void observable_set::initialization_unit(inputparameter_set & isi){
  static Logger logger("observable_set::initialization_unit (isi)");
  logger << LOG_DEBUG << "called" << endl;

  unit_calculation = isi.unit_calculation;
  unit_result = isi.unit_result;
  unit_distribution = isi.unit_distribution;

  logger << LOG_DEBUG << "unit_calculation  = " << unit_calculation << endl;
  logger << LOG_DEBUG << "unit_result       = " << unit_result << endl;
  logger << LOG_DEBUG << "unit_distribution = " << unit_distribution << endl;

  // used in the output files in 'result' and 'distribution':
  unit_factor_calculation = determine_unit_factor(unit_calculation);
  unit_factor2_calculation = pow(unit_factor_calculation, 2);
  // only needed for final collection of 'result':
  unit_factor_result = determine_unit_factor(unit_result) / unit_factor_calculation;
  // only needed for final collection of 'distribution':
  unit_factor_distribution = determine_unit_factor(unit_distribution) / unit_factor_calculation;
  
  logger << LOG_DEBUG << "unit_calculation  = " << unit_calculation << "   ->   " << unit_calculation << " / fb = " << unit_factor_calculation << endl;
  logger << LOG_DEBUG << "unit_result       = " << unit_result << "   ->   " << unit_result << " / " << unit_calculation << " = " << unit_factor_result << endl;
  logger << LOG_DEBUG << "unit_distribution = " << unit_distribution << "   ->   " << unit_distribution << " / " << unit_calculation << " = " << unit_factor_distribution << endl;

  logger << LOG_DEBUG << "finished" << endl;
}



double observable_set::determine_unit_factor(string & unit){
  static Logger logger("determine_unit_factor");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG << "unit = " << unit << endl;

  double unit_factor;
  // determines the relative factor wrt. fb:
  if      (unit == "fb"){unit_factor = 1.;}
  else if (unit == "pb"){unit_factor = 1.e-3;}
  else if (unit == "nb"){unit_factor = 1.e-6;}
  else if (unit == "Âµb"){unit_factor = 1.e-9;}
  else if (unit == "mb"){unit_factor = 1.e-12;}
  else if (unit == "ab"){unit_factor = 1.e3;}
  else if (unit == "zb"){unit_factor = 1.e6;}
  else if (unit == "yb"){unit_factor = 1.e9;}
  else {
    logger << LOG_FATAL << "No valid unit chosen: " << unit << " !" << endl;
    exit(1);
  }

  logger << LOG_DEBUG << "unit_factor = " << unit_factor << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return unit_factor;
}



void observable_set::initialization_object_event_selection(inputparameter_set & isi){
  static Logger logger("observable_set::initialization_object_event_selection (isi)");
  logger << LOG_DEBUG << "called" << endl;

  photon_recombination = isi.photon_recombination;
  photon_R_definition = isi.photon_R_definition;
  photon_R = isi.photon_R;
  photon_R2 = pow(photon_R, 2);
  photon_E_threshold_ratio = isi.photon_E_threshold_ratio;
  photon_jet_algorithm = isi.photon_jet_algorithm;
  // check if treated correctly !!!
  photon_recombination_selection = isi.photon_recombination_selection;
  photon_recombination_disable = isi.photon_recombination_disable;
  determine_selection_content_particle(photon_recombination_list, photon_recombination_selection, photon_recombination_disable);
  for (int i_j = 0; i_j < photon_recombination_list.size(); i_j++){
    logger << LOG_DEBUG << "photon_recombination_list[" << i_j << "] = " << photon_recombination_list[i_j] << endl;
  }
  photon_photon_recombination = isi.photon_photon_recombination;
  photon_photon_recombination_R = isi.photon_photon_recombination_R;

  // check if there is a photon to be recombined:
  int flag_photon = 0;
  for (int i_a = 0; i_a < csi->type_parton.size(); i_a++){
    logger << LOG_DEBUG << "csi->type_parton[" << i_a << "].size() = " << csi->type_parton[i_a].size() << endl;
    for (int i_p = 1; i_p < csi->type_parton[i_a].size(); i_p++){
      //    for (int i_p = 1; i_p < csi->type_parton.size(); i_p++){
      logger << LOG_DEBUG << "csi->type_parton[" << i_a << "][" << i_p << "] = " << csi->type_parton[i_a][i_p] << endl;
      //      if (csi->type_parton[i_a][i_p] == 21){flag_photon++;}// break;} // photon should be 22 !!!
      if (csi->type_parton[i_a][i_p] == 22){flag_photon++;}// break;}
    }
  }

  logger << LOG_DEBUG << "flag_photon = " << flag_photon << endl;

  if (!flag_photon){
    photon_recombination = 0;
    photon_recombination_list.clear();
  }

  frixione_isolation = isi.frixione_isolation;
  frixione_n = isi.frixione_n;
  frixione_epsilon = isi.frixione_epsilon;
  frixione_fixed_ET_max = isi.frixione_fixed_ET_max;
  frixione_delta_0 = isi.frixione_delta_0;
  frixione_jet_removal = isi.frixione_jet_removal;

  jet_algorithm = isi.jet_algorithm;
  jet_R_definition = isi.jet_R_definition;
  jet_R = isi.jet_R;
  jet_R2 = pow(jet_R, 2);
  parton_y_max = isi.parton_y_max;
  parton_eta_max = isi.parton_eta_max;
  jet_algorithm_selection = isi.jet_algorithm_selection;
  jet_algorithm_disable = isi.jet_algorithm_disable;
  logger << LOG_DEBUG << "jet_algorithm_list" << endl;
  if (jet_algorithm_selection.size() == 0){
    jet_algorithm_selection.push_back("g+Q");
    logger << LOG_INFO << "jet_algorithm_selection is set to its default value 'g+Q'" << endl;
  }
  determine_selection_content_particle(jet_algorithm_list, jet_algorithm_selection, jet_algorithm_disable);

  for (int i_j = 0; i_j < jet_algorithm_list.size(); i_j++){
    logger << LOG_DEBUG << "jet_algorithm_list[" << i_j << "] = " << jet_algorithm_list[i_j] << endl;
  }

  logger << LOG_DEBUG << "finished" << endl;
}


/*
void observable_set::initialization_object_process(){
  static Logger logger("observable_set::initialization_object_process (---)");
  logger << LOG_DEBUG << "called" << endl;

  determine_object_definition();
  //    determine_object_definition(n_observed_min, n_observed_max, define_pT, define_ET, define_eta, define_y, esi.observed_object, esi.object_list);

  esi.determine_n_partonlevel(csi->type_parton);
  
  for (int i_e = 0; i_e < esi.pda.size(); i_e++){n_partonlevel[i_e] = esi.pda[i_e].n_partonlevel;}

  logger << LOG_DEBUG << "finished" << endl;
}
*/


void observable_set::initialization_filename(){
  static Logger logger("observable_set::initialization_filename (---)");
  logger << LOG_DEBUG << "called" << endl;

  int isystem = 0; 
  string xorder = "";
  //  string temp_directory;

  //  string prozessdatname; -> name_process
  //  map<int, string> datname;
  //  fill_datname(datname);

  string prozessgnuplot = "";
  map<int, string> gnuplotname;
  fill_gnuplotname(gnuplotname);

  if (process_type == 2){
    //    name_process = datname[csi->type_parton[0][1]] + datname[csi->type_parton[0][2]] + "_";
    //    prozessgnuplot = gnuplotname[csi->type_parton[0][1]] + gnuplotname[csi->type_parton[0][2]] + "_";
    prozessgnuplot = gnuplotname[csi->type_parton[0][1]] + gnuplotname[csi->type_parton[0][2]] + "&rightarrow ";
    for (int i = 3; i < csi->type_parton[0].size(); i++){
      //      name_process = name_process + datname[csi->type_parton[0][i]];
      prozessgnuplot = prozessgnuplot + gnuplotname[csi->type_parton[0][i]];
    }
  }
  else if (process_type == 1){
    //    name_process = datname[csi->type_parton[0][0]] + "_";
    prozessgnuplot = gnuplotname[csi->type_parton[0][0]] + "&rightarrow ";
    for (int i = 1; i < csi->type_parton[0].size(); i++){
      //      name_process = name_process + datname[csi->type_parton[0][i]];
      prozessgnuplot = prozessgnuplot + gnuplotname[csi->type_parton[0][i]];
    }
  }

  name_process = csi->subprocess;

  logger << LOG_DEBUG << "process_type = " << process_type << endl;
  logger << LOG_DEBUG << "name_process = " << name_process << endl;
  logger << LOG_DEBUG << "csi->subprocess = " << csi->subprocess << endl;



  //  string dir_madgraph = "madgraph";
  // subcontribution + "_"+ in essentially all filename_...s



  string dir_integration = "integration";
  if (switch_output_integration){system_execute(logger, "mkdir " + dir_integration, isystem);}
  filename_integration = dir_integration + "/integration_" + name_process + ".txt";
  logger << LOG_DEBUG << "filename_integration = " << filename_integration << endl;
  //switch_output_integration && 
  if (switch_CV){
    //    temp_directory = dir_integration + "/CV";
    if (switch_output_integration){system_execute(logger, "mkdir " + dir_integration + "/CV", isystem);}
    filename_integration_CV = dir_integration + "/CV" + "/integration_" + name_process + ".txt";
  }
  if (switch_TSV){
    filename_integration_TSV.resize(n_extended_set_TSV);
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      logger << LOG_DEBUG << "i_s = " << i_s << endl;
      logger << LOG_DEBUG << "n_extended_set_TSV = " << n_extended_set_TSV << endl;
      logger << LOG_DEBUG << "name_extended_set_TSV.size() = " << name_extended_set_TSV.size() << endl;
      logger << LOG_DEBUG << "name_extended_set_TSV[" << i_s << "] = " << name_extended_set_TSV[i_s] << endl;
      //      string temp_directory = dir_integration + "/" + name_extended_set_TSV[i_s];
      //      logger << LOG_DEBUG << "temp_directory = " << temp_directory << endl;
      if (switch_output_integration){system_execute(logger, "mkdir " + dir_integration + "/" + name_extended_set_TSV[i_s], isystem);}
      filename_integration_TSV[i_s] = dir_integration + "/" + name_extended_set_TSV[i_s] + "/integration_" + name_process + ".txt";
    }
  }

  string dir_result = "result";
  if (switch_output_result){system_execute(logger, "mkdir " + dir_result, isystem);}
  filename_result = dir_result + "/result_" + name_process + ".dat";
  if (switch_TSV){
    filename_result_TSV.resize(n_extended_set_TSV);
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      //      string temp_directory = dir_result + "/" + name_extended_set_TSV[i_s];
      if (switch_output_result){system_execute(logger, "mkdir " + dir_result + "/" + name_extended_set_TSV[i_s], isystem);}
      filename_result_TSV[i_s] = dir_result + "/" + name_extended_set_TSV[i_s] + "/result_" + name_process + ".txt";
    }
  }

  string dir_maxevent = "maxevent";
  if (switch_output_maxevent){system_execute(logger, "mkdir " + dir_maxevent, isystem);}
  filename_maxevent = dir_maxevent + "/maxevent_" + name_process + ".txt";


  string dir_comparison = "comparison";
  if (switch_output_comparison){system_execute(logger, "mkdir " + dir_comparison, isystem);}
  filename_comparison = dir_comparison + "/comparison_" + name_process + ".txt";


  string dir_execution = "execution";
  if (switch_output_execution){system_execute(logger, "mkdir " + dir_execution, isystem);}
  filename_execution = dir_execution + "/execution_" + name_process + ".dat";


  string dir_proceeding = "proceeding";
  if (switch_output_proceeding){system_execute(logger, "mkdir " + dir_proceeding, isystem);}
  filename_proceeding = dir_proceeding + "/proceeding_" + name_process + ".dat";
  filename_proceeding_2 = dir_proceeding + "/proceeding_2_" + name_process + ".dat";


  string dir_gnuplot = "gnuplot";
  if (switch_output_gnuplot){system_execute(logger, "mkdir " + dir_gnuplot, isystem);}
  filename_gnuplot = dir_gnuplot + "/gnuplot_" + name_process + ".dat";

  filename_makegnuplot = dir_gnuplot + "/gnuplot_" + name_process + ".gp";
  replace(filename_makegnuplot.begin(), filename_makegnuplot.end(), '~', 'x');

  if (name_process != ""){
    // **************************************************************************
    // *                                                                        *
    // *  creation of gnuplot deviation check file                              *
    // *                                                                        *
    // **************************************************************************

    ifstream in_gnuplot("../../../../visualize.error/mask.visualize.error.dat");
    //  ifstream in_gnuplot("gnuplot/sigmaplot.maske");
    vector<string> readin;
    //    readin = vs0;
    //    int counter = 0;
    char LineBuffer[256];
    while (in_gnuplot.getline(LineBuffer, 256)){readin.push_back(LineBuffer);}
    ofstream out_makegnuplot;
    out_makegnuplot.open(filename_makegnuplot.c_str(), ofstream::out | ofstream::trunc);
    for (int i = 0; i < readin.size(); i++){
      if (i == 10){
	out_makegnuplot << readin[i] + readin[i + 1] << endl;
	i++;
      }
      else if (i == 18){
	out_makegnuplot << readin[i] << " " << char(34) << char(92) << "$" << prozessgnuplot  << char(92) << "$" << char(34) << " " << readin[i + 1] << endl;
	i++;
      }
      else if (i == 120){
	out_makegnuplot << readin[i] << " '../" << filename_gnuplot << "' " << readin[i + 1] << endl;
	i++;
      }
      else{out_makegnuplot << readin[i] << endl;}
    }
    out_makegnuplot.close();
    
    xorder = "chmod 700 " + filename_makegnuplot;
    logger << LOG_DEBUG << xorder << endl;
    isystem = system(xorder.c_str());
    logger << LOG_DEBUG << xorder << "   execution status: " << isystem << "." << endl;
  }






  string dir_moment = "moment";
  if (switch_output_moment){system_execute(logger, "mkdir " + dir_moment, isystem);}
  filename_moment.resize(n_moments);
  for (int nm = 0; nm < n_moments; nm++){
    char nmstring[255];
    sprintf(nmstring, "%d", nm);
    filename_moment[nm] = dir_moment + "/moment_" + nmstring + "_" + name_process + ".dat";
  }


  if (switch_TSV){
    filename_moment_TSV.resize(n_extended_set_TSV);
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      string temp_directory = dir_moment + "/" + name_extended_set_TSV[i_s];
      if (switch_output_moment){system_execute(logger, "mkdir " + temp_directory, isystem);}
      filename_moment_TSV[i_s] = temp_directory + "/moment_" + name_process + ".txt";
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void observable_set::initialization_switches(inputparameter_set & isi){
  static Logger logger("observable_set::initialization_switches (isi)");
  logger << LOG_DEBUG << "called" << endl;

  switch_result = isi.switch_result;
  switch_distribution = isi.switch_distribution;
  switch_moment = isi.switch_moment;

  switch_output_execution = isi.switch_output_execution;
  switch_output_integration = isi.switch_output_integration;
  switch_output_maxevent = isi.switch_output_maxevent;
  switch_output_comparison = isi.switch_output_comparison;
  switch_output_gnuplot = isi.switch_output_gnuplot;
  switch_output_proceeding = isi.switch_output_proceeding;
  switch_output_weights = isi.switch_output_weights;
  switch_output_result = isi.switch_output_result;
  switch_output_moment = isi.switch_output_moment;
  switch_output_distribution = isi.switch_output_distribution;

  switch_output_testpoint = isi.switch_output_testpoint;
  switch_output_cancellation_check = isi.switch_output_cancellation_check;
  switch_output_cutinfo = isi.switch_output_cutinfo;

  switch_console_output_runtime = isi.switch_console_output_runtime;
  switch_console_output_tau_0 = isi.switch_console_output_tau_0;
  switch_console_output_techcut_RA = isi.switch_console_output_techcut_RA;
  switch_console_output_phasespace_issue = isi.switch_console_output_phasespace_issue;
  switch_console_output_ME2_issue = isi.switch_console_output_ME2_issue;

  switch_testcut = isi.switch_testcut;

  switch_polenorm = isi.switch_polenorm;

  switch_old_qT_version = isi.switch_old_qT_version;

  switch_RS = isi.switch_RS;
  switch_RS_mapping = isi.switch_RS_mapping;

  switch_VI = isi.switch_VI;
  switch_VI_bosonic_fermionic = isi.switch_VI_bosonic_fermionic;
  switch_KP = isi.switch_KP;
  switch_CM = isi.switch_CM;
  switch_OL = isi.switch_OL;
  switch_yuk = isi.switch_yuk;
  order_y = isi.order_y;
    
  switch_H1gg = isi.switch_H1gg;
  switch_H2 = isi.switch_H2;

  if (switch_KP != 0 && switch_KP != 1 && switch_KP != 2){
    logger << LOG_FATAL << "Invalid value of switch_KP!" << endl;
    exit(1);
  }

  if (switch_VI != 0 && switch_VI != 1 && switch_VI != 2){
    logger << LOG_FATAL << "Invalid value of switch_VI!" << endl;
    exit(1);
  }

  if (switch_H1gg != 0 && switch_H1gg != 1){
    logger << LOG_FATAL << "Invalid value of switch_H1gg!" << endl;
    exit(1);
  }

  if (switch_H2 != 0 && switch_H2 != 1){
    logger << LOG_FATAL << "Invalid value of switch_H2!" << endl;
    exit(1);
  }


  logger << LOG_DEBUG << "finished" << endl;
}



void observable_set::initialization_qTcut(inputparameter_set & isi){
  static Logger logger("observable_set::initialization_qTcut (isi)");
  logger << LOG_DEBUG << "called" << endl;

  switch_qTcut = isi.switch_qTcut;
  n_qTcut = isi.n_qTcut;
  min_qTcut = isi.min_qTcut;
  max_qTcut = isi.max_qTcut;
  step_qTcut = isi.step_qTcut;
  binning_qTcut = isi.binning_qTcut;
  selection_qTcut = isi.selection_qTcut;
  selection_qTcut_distribution = isi.selection_qTcut_distribution;
  selection_no_qTcut_distribution = isi.selection_no_qTcut_distribution;
  selection_qTcut_result = isi.selection_qTcut_result;
  selection_no_qTcut_result = isi.selection_no_qTcut_result;
  selection_qTcut_integration = isi.selection_qTcut_integration;
  selection_no_qTcut_integration = isi.selection_no_qTcut_integration;

  switch_NJcut = isi.switch_NJcut;
  switch_NJcut_axes = isi.switch_NJcut_axes;
  switch_NJcut_axes_energy = isi.switch_NJcut_axes_energy;
  switch_NJcut_measure = isi.switch_NJcut_measure;

  
  if (!switch_qTcut){
    n_qTcut = 1;
    min_qTcut = 0;
    step_qTcut = 0;
    max_qTcut = 0;
    value_qTcut = vector<double> (1, 0.);
    logger << LOG_INFO << "No qTcut variation applied!" << endl;
    logger << LOG_INFO << setw(40) << "switch_qTcut" << "     =     " << setw(20) << switch_qTcut << endl;
    logger << LOG_INFO << setw(40) << "n_qTcut" << "     =     " << setw(20) << n_qTcut << endl;
    logger << LOG_INFO << setw(40) << "min_qTcut" << "     =     " << setw(20) << min_qTcut << endl;
    logger << LOG_INFO << setw(40) << "max_qTcut" << "     =     " << setw(20) << max_qTcut << endl;
    logger << LOG_INFO << setw(40) << "step_qTcut" << "     =     " << setw(20) << step_qTcut << endl;
    logger << LOG_INFO << setw(40) << "binning_qTcut" << "     =     " << setw(20) << binning_qTcut << endl;
    logger << LOG_INFO << setw(40) << "selection_qTcut" << "     =     " << setw(20) << selection_qTcut << endl;
    
  }
  else {
    if (switch_qTcut && !n_qTcut){exit(1);}
    //  value_qTcut = isi.value_qTcut;

    if (binning_qTcut == "linear"){
      // possible ways to set input parameters:
      // min_qTcut, n_qTcut, step_qTcut, max_qTcut
      if (min_qTcut != 0. && n_qTcut != 0 && step_qTcut != 0. && max_qTcut != 0.){
	double temp_max_qTcut = min_qTcut + (n_qTcut - 1) * step_qTcut;
	if (max_qTcut != temp_max_qTcut){logger << LOG_FATAL << "Inconsistently over-defined qTcut input!" << endl; exit(1);}
      }
      // min_qTcut, n_qTcut, step_qTcut
      else if (min_qTcut != 0. && n_qTcut != 0 && step_qTcut != 0.){
	max_qTcut = min_qTcut + (n_qTcut - 1) * step_qTcut;
      }
      // min_qTcut, n_qTcut, max_qTcut
      else if (min_qTcut != 0. && n_qTcut != 0 && max_qTcut != 0.){
	step_qTcut = (max_qTcut - min_qTcut) / (n_qTcut - 1);
      }
      // min_qTcut, step_qTcut, max_qTcut
      else if (min_qTcut != 0. && step_qTcut != 0. && max_qTcut != 0.){
	n_qTcut = (max_qTcut - min_qTcut) / step_qTcut + 1;
      }
      // n_qTcut, step_qTcut, max_qTcut
      else if (min_qTcut != 0. && n_qTcut != 0 && step_qTcut != 0. && max_qTcut != 0.){
	min_qTcut = max_qTcut - (n_qTcut - 1) * step_qTcut;
      }
      else {
	logger << LOG_FATAL << "Inconsistent qTcut input for linear binning!" << endl;
	exit(1);
      }
      for (int i_q = 0; i_q < n_qTcut; i_q++){
	value_qTcut.push_back(min_qTcut + i_q * step_qTcut);
      }
    }
    
    else if (binning_qTcut == "logarithmic"){
      // possible ways to set input parameters:
      // min_qTcut, n_qTcut, max_qTcut
      if (min_qTcut != 0. && n_qTcut != 0 && max_qTcut != 0.){
	step_qTcut = 0.; // no constant step width in logarithmic binning!
      }
      else {
	logger << LOG_FATAL << "Inconsistent qTcut input for logarithmic binning!" << endl;
	exit(1);
      }
      for (int i_q = 0; i_q < n_qTcut; i_q++){
	value_qTcut.push_back(min_qTcut * exp10(log10(max_qTcut / min_qTcut) * double(i_q) / (n_qTcut - 1)));
      }
      value_qTcut[n_qTcut - 1] = max_qTcut;
    }

    else if (binning_qTcut == "irregular"){
      // possible ways to set input parameters:

      // selection_qTcut
      if (selection_qTcut != ""){
	vector<string> vs_selection(1);
	for (int i_b = 0; i_b < selection_qTcut.size(); i_b++){
	  logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << "   vs_selection.size() = " << vs_selection.size() << endl;
	  if (selection_qTcut[i_b] == ':'){vs_selection.push_back("");}
	  else if (selection_qTcut[i_b] != ':'){vs_selection[vs_selection.size() - 1].push_back(selection_qTcut[i_b]);}
	}
	n_qTcut = vs_selection.size();
	value_qTcut.resize(n_qTcut);
	for (int i_b = 0; i_b < value_qTcut.size(); i_b++){
	  logger << LOG_DEBUG_VERBOSE << "vs_selection[" << i_b << "] = " << setprecision(8) << setw(15) << vs_selection[i_b] << endl;
	  value_qTcut[i_b] = atof(vs_selection[i_b].c_str());
	}
	min_qTcut = value_qTcut[0];
	max_qTcut = value_qTcut[n_qTcut - 1];
      }
      else {
	logger << LOG_FATAL << "Inconsistent qTcut input for logarithmic binning!" << endl;
	exit(1);
      }
    }
  }

  logger << LOG_INFO << setw(40) << "n_qTcut" << "     =     " << setw(20) << n_qTcut << endl;
  logger << LOG_INFO << setw(40) << "min_qTcut" << "     =     " << setw(20) << min_qTcut << endl;
  logger << LOG_INFO << setw(40) << "max_qTcut" << "     =     " << setw(20) << max_qTcut << endl;
  logger << LOG_INFO << setw(40) << "step_qTcut" << "     =     " << setw(20) << step_qTcut << endl;
  logger << LOG_INFO << setw(40) << "binning_qTcut" << "     =     " << setw(20) << binning_qTcut << endl;
  logger << LOG_INFO << setw(40) << "selection_qTcut" << "     =     " << setw(20) << selection_qTcut << endl;
  for (int i_b = 0; i_b < value_qTcut.size(); i_b++){
    logger << LOG_INFO << "value_qTcut[" << setw(3) << i_b << "] = " << value_qTcut[i_b] << endl;
  }


      

  if (csi->type_contribution == "" ||
      csi->type_contribution == "CT" ||
      csi->type_contribution == "CJ" ||
      csi->type_contribution == "RT" ||
      csi->type_contribution == "CT2" ||
      csi->type_contribution == "CJ2" ||
      csi->type_contribution == "RVA" ||
      csi->type_contribution == "RCA" ||
      csi->type_contribution == "RRA" ||
      csi->type_contribution == "L2RT" ||
      csi->type_contribution == "L2CT"){
    active_qTcut = 1;
    output_n_qTcut = n_qTcut;
  }
  else {
    active_qTcut = 0;
    output_n_qTcut = 1;
  }

  logger << LOG_INFO << "selection_qTcut_distribution = " << selection_qTcut_distribution << endl;
  logger << LOG_INFO << "selection_no_qTcut_distribution = " << selection_no_qTcut_distribution << endl;
  logger << LOG_INFO << "selection_qTcut_result = " << selection_qTcut_result << endl;
  logger << LOG_INFO << "selection_no_qTcut_result = " << selection_no_qTcut_result << endl;
  logger << LOG_INFO << "selection_qTcut_integration = " << selection_qTcut_integration << endl;
  logger << LOG_INFO << "selection_no_qTcut_integration = " << selection_no_qTcut_integration << endl;
  
  for (int i_m = 0; i_m < 3; i_m++){
    string temp_name;
    string temp_selection_qTcut;
    string temp_selection_no_qTcut;
    vector<int> temp_no_qTcut;
    vector<double> temp_value_qTcut;
    if (i_m == 0){temp_selection_qTcut = selection_qTcut_distribution; temp_selection_no_qTcut = selection_no_qTcut_distribution; temp_name = "distribution";}
    else if (i_m == 1){temp_selection_qTcut = selection_qTcut_result; temp_selection_no_qTcut = selection_no_qTcut_result; temp_name = "result";}
    else if (i_m == 2){temp_selection_qTcut = selection_qTcut_integration; temp_selection_no_qTcut = selection_no_qTcut_integration; temp_name = "integration";}
    else {logger << LOG_ERROR << "Wrong entry in selection_qTcut!" << endl;}
    
    if (active_qTcut){
      if (temp_selection_qTcut == "" && temp_selection_no_qTcut == ""){
	temp_value_qTcut.resize(1, 0.);
	temp_no_qTcut.resize(1, 0);
	logger << LOG_INFO << temp_name << " output generated for lowest qTcut value." << endl;
      }
      else if (temp_selection_qTcut == "all" || temp_selection_no_qTcut == "all"){
	temp_value_qTcut = value_qTcut;
	temp_no_qTcut.resize(n_qTcut);
	for (int i_q = 0; i_q < n_qTcut; i_q++){temp_no_qTcut[i_q] = i_q;}
	logger << LOG_INFO << temp_name << " output generated for all qTcut values." << endl;
      }
      else if (temp_selection_qTcut != "" && temp_selection_no_qTcut == ""){
	logger << LOG_INFO << temp_name << " output generated for selected qTcut values." << endl;
	
	vector<string> vs_selection(1);
	for (int i_b = 0; i_b < temp_selection_qTcut.size(); i_b++){
	  logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << "   vs_selection.size() = " << vs_selection.size() << endl;
	  if (temp_selection_qTcut[i_b] == ':'){vs_selection.push_back("");}
	  else if (temp_selection_qTcut[i_b] != ':'){vs_selection[vs_selection.size() - 1].push_back(temp_selection_qTcut[i_b]);}
	}
	//    int n_selection = vs_selection.size();
	temp_no_qTcut.resize(vs_selection.size());
	temp_value_qTcut.resize(vs_selection.size());
	
	for (int i_b = 0; i_b < temp_value_qTcut.size(); i_b++){
	  logger << LOG_DEBUG_VERBOSE << "vs_selection[" << i_b << "] = " << setprecision(8) << setw(15) << vs_selection[i_b] << endl;
	  temp_value_qTcut[i_b] = atof(vs_selection[i_b].c_str());
	  int flag = n_qTcut;
	  for (int i_q = 0; i_q < n_qTcut; i_q++){
	    if (abs(temp_value_qTcut[i_b] - value_qTcut[i_q]) < 1.e-12 * value_qTcut[i_q]){flag = i_q; break;}
	  }
	  if (flag == n_qTcut){temp_value_qTcut.erase(temp_value_qTcut.begin() + i_b); i_b--;}
	  else {temp_no_qTcut[i_b] = flag;}
	}
      }

      else if (temp_selection_qTcut == "" && temp_selection_no_qTcut != ""){
	logger << LOG_INFO << temp_name << " output generated for selected no_qTcut values." << endl;
	vector<string> vs_selection(1);
	for (int i_b = 0; i_b < temp_selection_no_qTcut.size(); i_b++){
	  logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << "   vs_selection.size() = " << vs_selection.size() << endl;
	  if (temp_selection_no_qTcut[i_b] == ':'){vs_selection.push_back("");}
	  else if (temp_selection_no_qTcut[i_b] != ':'){vs_selection[vs_selection.size() - 1].push_back(temp_selection_no_qTcut[i_b]);}
	}
	//    int n_selection = vs_selection.size();
	temp_no_qTcut.resize(vs_selection.size());
	temp_value_qTcut.resize(vs_selection.size());

	for (int i_b = 0; i_b < temp_no_qTcut.size(); i_b++){
	  logger << LOG_DEBUG << "vs_selection[" << i_b << "] = " << setprecision(8) << setw(15) << vs_selection[i_b] << endl;
	  temp_no_qTcut[i_b] = atoi(vs_selection[i_b].c_str());
	  temp_value_qTcut[i_b] = value_qTcut[temp_no_qTcut[i_b]];
	}
      }

      else if (temp_selection_qTcut != "" && temp_selection_no_qTcut != ""){
	logger << LOG_FATAL << "Inconsistent qTcut input for " << temp_name << "!" << endl;
	exit(1);
      }
      
      logger << LOG_DEBUG << temp_name << "   temp_value_qTcut.size() = " << temp_value_qTcut.size() << endl;
      for (int i_b = 0; i_b < temp_value_qTcut.size(); i_b++){
	logger << LOG_INFO << "value_qTcut_" << temp_name << "[" << i_b << "] = " << setw(15) << setprecision(8) << temp_value_qTcut[i_b] << " -> " << setw(3) << temp_no_qTcut[i_b] << endl;
      }
    }

    else {
      temp_value_qTcut.resize(1, 0.);
      temp_no_qTcut.resize(1, 0);
      for (int i_b = 0; i_b < temp_value_qTcut.size(); i_b++){
	logger << LOG_INFO << "value_qTcut_" << temp_name << "[" << i_b << "] = " << setw(15) << setprecision(8) << temp_value_qTcut[i_b] << " -> " << setw(3) << temp_no_qTcut[i_b] << endl;
      }
    }
    if (i_m == 0){no_qTcut_distribution = temp_no_qTcut; value_qTcut_distribution = temp_value_qTcut;}
    else if (i_m == 1){no_qTcut_result = temp_no_qTcut; value_qTcut_result = temp_value_qTcut;}
    else if (i_m == 2){no_qTcut_integration = temp_no_qTcut; value_qTcut_integration = temp_value_qTcut;}

    //    logger << LOG_INFO << "temp_name = " << temp_name << endl;
    //    logger << LOG_INFO << "temp_selection_qTcut = " << temp_selection_qTcut << endl;
  }

  if (csi->type_contribution == "RVA" ||
      csi->type_contribution == "L2RT" ||
      csi->type_contribution == "L2RJ"){
    counter_killed_qTcut.resize(n_qTcut, 0);
    counter_acc_qTcut.resize(n_qTcut, 0);
  }


  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_LHAPDF_parameters(inputparameter_set & isi){
  static Logger logger("observable_set::initialization_LHAPDF (isi)");
  logger << LOG_DEBUG << "called" << endl;

  N_f_active = isi.N_f_active;
  pdf_selection = isi.pdf_selection;
  pdf_disable = isi.pdf_disable;
  //  pdf_content_modify = isi.pdf_content_modify;

  LHAPDFname = isi.LHAPDFname;
  LHAPDFsubset = isi.LHAPDFsubset;

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_integration_parameters(inputparameter_set & isi){
  static Logger logger("observable_set::initialization_integration_parameters (isi)");
  logger << LOG_DEBUG << "called" << endl;

  sigma_normalization = isi.sigma_normalization;
  sigma_normalization_deviation = isi.sigma_normalization_deviation;

  coll_choice = isi.coll_choice;

  n_moments = isi.n_moments;

  N_f = isi.N_f;
  N_quarks = isi.N_quarks;

  check_vanishing_ME2_end = 0;
  flag_vanishing_ME2 = 0;
  n_event_vanishing_ME2 = 100; // temporary !!!

  temp_n_step = 0;

  logger << LOG_DEBUG << "finished" << endl;
}
/*
void observable_set::initialization_generic(inputparameter_set & isi){
  static Logger logger("observable_set::initialization_generic (isi)");
  logger << LOG_DEBUG << "called" << endl;


  logger << LOG_DEBUG << "finished" << endl;
}
// to be shifted !!!
*/


/*
void observable_set::initialization_xxx(inputparameter_set & isi){
  static Logger logger("observable_set::initialization_xxx (isi)");
  logger << LOG_DEBUG << "called" << endl;

  logger << LOG_DEBUG << "finished" << endl;
}
*/

//#include "OV.observable.initialization.cpp"














  




void observable_set::initialization_distribution(vector<xdistribution> _dat, vector<dddistribution> _dddat){
  static Logger logger("observable_set::initialization_distribution");
  logger << LOG_DEBUG << "started" << endl;

  dat = _dat;
  dddat = _dddat;

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_masses(vector<double> _M, vector<double> _M2){
  static Logger logger("observable_set::initialization_masses");
  logger << LOG_DEBUG << "started" << endl;

  M = _M;
  M2 = _M2;

  for (int i = 0; i < 26; i++){logger << LOG_DEBUG << "M[" << setw(2) << i << "]  = " << right << setprecision(15) << setw(23) << M[i] << endl;}

  logger << LOG_DEBUG << "finished" << endl;

}

void observable_set::initialization_integration(phasespace_set & psi){//int _n_particle, vector<int> _csi->type_parton){
  static Logger logger("observable_set::initialization_integration");
  logger << LOG_DEBUG << "started" << endl;

  //  n_particle = _n_particle;
  //  int_end = 0;
  logger << LOG_DEBUG << "filename_proceeding = " << filename_proceeding << endl;


  relation_pc_ps.resize(n_pc);
  for (int i_c = 0; i_c < n_pc; i_c++){
    if (n_ps == 1){relation_pc_ps[i_c] = 0;}
    else {relation_pc_ps[i_c] = i_c;}
  }

  cut_ps.resize(n_ps, 0);
  first_non_cut_ps = 0;

  p_parton.resize(n_ps, vector<fourvector> (n_particle + 3, fourvector ()));
  boost = 0.;
  
  logger << LOG_DEBUG << "csi->type_parton.size() = " << csi->type_parton.size() << endl;
  for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){logger << LOG_DEBUG << "csi->type_parton[0][" << i_p << "] = " << csi->type_parton[0][i_p] << endl;}
  //  csi->type_parton.resize(n_ps, vector<int> (n_particle + 3, 0));
  /*
  cout << "csi->type_parton.size() = " << csi->type_parton.size() << endl;
  for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){cout << "csi->type_parton[0][" << i_p << "] = " << csi->type_parton[0][i_p] << endl;}
  exit(1);
  */
  //  csi->type_parton[0] = _csi->type_parton;
  mass_parton.resize(n_ps, vector<double> (n_particle + 3, 0.));
  mass2_parton.resize(n_ps, vector<double> (n_particle + 3, 0.));

  for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){
    mass_parton[0][i_p] = M[abs(csi->type_parton[0][i_p])];
    mass2_parton[0][i_p] = M2[abs(csi->type_parton[0][i_p])];
  }

  // determination if massless or massive subtraction is needed:
  massive_QCD = 0; 
  for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){
    if (mass_parton[0][i_p] > 0. && abs(csi->type_parton[0][i_p]) < 7){massive_QCD = 1; break;}
  }
  logger << LOG_DEBUG << "massive_QCD = " << massive_QCD << endl;

  massive_QEW = 0; 
  for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){
    if (mass_parton[0][i_p] > 0. && (abs(csi->type_parton[0][i_p]) < 7 || abs(csi->type_parton[0][i_p]) == 11 || abs(csi->type_parton[0][i_p]) == 13 || abs(csi->type_parton[0][i_p]) == 15 || abs(csi->type_parton[0][i_p]) == 24)){massive_QEW = 1; break;}
  }
  logger << LOG_DEBUG << "massive_QEW = " << massive_QEW << endl;

  logger << LOG_DEBUG << "n_ps = " << n_ps << endl;
  logger << LOG_DEBUG << "n_particle = " << n_particle << endl;
  


  moment_symm.resize(n_moments, 0);
  moment.resize(n_moments, vector<double> (n_ps, 0.));
  directed_moment.resize(n_moments, vector<vector<double> > (3, vector<double> (n_ps, 0.)));

  max_integrand = 0.;

  integrand = 0.;
  integrand_CV.resize(n_scales_CV, 0.);
  integrand_qTcut_CV.resize(n_qTcut, integrand_CV);

  integrand_D.resize(3, vector<double> (n_ps, 0.));
  integrand_D_CV.resize(n_scales_CV, integrand_D);
  integrand_D_qTcut_CV.resize(n_qTcut, integrand_D_CV);
  

  change_cut.resize(n_qTcut + 1, 0);
  //  change_cut.resize(n_qTcut, 0);

  this_psp_weight = 0.;
  this_psp_weight2 = 0.;
  this_psp_weight_CV.resize(n_qTcut, vector<double> (n_scales_CV, 0.));
  this_psp_weight2_CV.resize(n_qTcut, vector<double> (n_scales_CV, 0.));
  
  step_sum_weight = 0.;
  step_sum_weight2 = 0.;
  step_sum_weight_CV.resize(n_qTcut, vector<double> (n_scales_CV, 0.));
  step_sum_weight2_CV.resize(n_qTcut, vector<double> (n_scales_CV, 0.));
  
  full_sum_weight = 0.;
  full_sum_weight2 = 0.;
  full_sum_weight_CV.resize(n_qTcut, vector<double> (n_scales_CV, 0.));
  full_sum_weight2_CV.resize(n_qTcut, vector<double> (n_scales_CV, 0.));
  
  this_psp_moment.resize(n_moments, 0.);
  this_psp_moment_CV.resize(n_moments, vector<vector<double> > (n_qTcut, vector<double> (n_scales_CV, 0.)));
  this_psp_moment2.resize(n_moments);
  this_psp_moment2_CV.resize(n_moments, vector<vector<double> > (n_qTcut, vector<double> (n_scales_CV, 0.)));

  step_sum_moment.resize(n_moments, 0.);
  step_sum_moment_CV.resize(n_moments, vector<vector<double> > (n_qTcut, vector<double> (n_scales_CV, 0.)));
  step_sum_moment2.resize(n_moments);
  step_sum_moment2_CV.resize(n_moments, vector<vector<double> > (n_qTcut, vector<double> (n_scales_CV, 0.)));

  full_sum_moment.resize(n_moments, 0.); 
  full_sum_moment2.resize(n_moments, 0.);
  full_sum_moment_CV.resize(n_moments, vector<vector<double> > (n_qTcut, vector<double> (n_scales_CV, 0.)));
  full_sum_moment2_CV.resize(n_moments, vector<vector<double> > (n_qTcut, vector<double> (n_scales_CV, 0.)));


  logger << LOG_DEBUG << "before switch_distribution" << endl;

  if (switch_distribution){
    bin.resize(dat.size() + dddat.size(), vector<int> (n_ps));
    bin_max.resize(dat.size() + dddat.size(), vector<int> (n_ps));
    
    bin_weight.resize(dat.size() + dddat.size());
    bin_weight2.resize(dat.size() + dddat.size());
    bin_counts.resize(dat.size() + dddat.size());
    bin_weight_CV.resize(n_scales_CV);
    bin_weight2_CV.resize(n_scales_CV);
    for (int i_d = 0; i_d < dat.size(); i_d++){
      bin_weight[i_d].resize(dat[i_d].n_bins);
      bin_weight2[i_d].resize(dat[i_d].n_bins);
      bin_counts[i_d].resize(dat[i_d].n_bins);
    }
    for (int i_d = 0; i_d < dddat.size(); i_d++){
      bin_weight[dat.size() + i_d].resize(dddat[i_d].n_bins);
      bin_weight2[dat.size() + i_d].resize(dddat[i_d].n_bins);
      bin_counts[dat.size() + i_d].resize(dddat[i_d].n_bins);
    }
    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      bin_weight_CV[i_s] = bin_weight;
      bin_weight2_CV[i_s] = bin_weight2;
    }
  }

  logger << LOG_DEBUG << "before particle_event.resize" << endl;

  particle_event.resize(esi.object_list.size(), vector<vector<particle> > (n_ps));
  
  particle_event[0][0].resize(n_particle + 3);
  for (int i_a = 1; i_a < n_ps; i_a++){particle_event[0][i_a].resize(n_particle + 3 - 1);}
  n_object.resize(esi.object_list.size(), vector<int> (n_ps));

  // !!! needed ???
  for (int i_o = 0; i_o < particle_event.size(); i_o++){
    for (int i_a = 0; i_a < particle_event[i_o].size(); i_a++){
      particle_event[i_o][i_a].reserve(n_particle + 3);
    }
  }


  logger << LOG_DEBUG_VERBOSE << "particle_event.size() = " << particle_event.size() << endl;
  for (int i_o = 0; i_o < particle_event.size(); i_o++){
    logger << LOG_DEBUG_VERBOSE << "particle_event[" << i_o << "].size() = " << particle_event[i_o].size() << endl;
    for (int i_a = 0; i_a < particle_event[i_o].size(); i_a++){
      logger << LOG_DEBUG_VERBOSE << "particle_event[" << i_o << "][" << i_a << "].size() = " << particle_event[i_o][i_a].size() << endl;
      for (int i_p = 0; i_p < particle_event[i_o][i_a].size(); i_p++){
	//	logger << LOG_DEBUG_VERBOSE << "particle_event[" << i_o << "][" << i_a << "][" << i_p << "] = " << particle_event[i_o][i_a][i_p] << endl;
	//	particle_event[i_o][i_a].resize(n_particle + 3 - 1);}
      }
    }
  }
  logger << LOG_DEBUG << "after particle_event.resize" << endl;

  
  recombination_history.resize(n_ps);
  
  
  start_p_parton = p_parton;
  //  start_p_event = p_event;
  //  start_p_event_object = p_event_object;
  start_particle_event = particle_event;
  start_n_object = n_object;
  
  h = 0;
  min = 0;
  sec = 0;
  time_counter = 0;

  sec_import = 0;

  /*
  if (csi->type_contribution == "RA" || 
      csi->type_contribution == "RRA"){
  */
  if (csi->class_contribution_CS_real){
    initialization_specific_RA();
  }
  /*
  if (csi->type_contribution == "CA" || 
      csi->type_contribution == "RCA"){
  */
    if (csi->class_contribution_CS_collinear){
    initialization_specific_CA();
  }
    /*
  if (csi->type_contribution == "VA" || 
      csi->type_contribution == "RVA"){
    */
    if (csi->class_contribution_CS_virtual){
    initialization_specific_VA();
  }
  if (csi->type_contribution == "VT" || 
      csi->type_contribution == "VT2" || 
      csi->type_contribution == "VJ" || 
      csi->type_contribution == "VJ2" || 
      csi->type_contribution == "L2VT" || 
      csi->type_contribution == "L2VJ"){
    initialization_specific_VT();
    if (switch_old_qT_version){initialization_specific_QT();}
  }
  if (csi->type_contribution == "CT2" ||
      csi->type_contribution == "CJ2"){
    initialization_specific_VT();
    initialization_specific_QT();
    initialization_specific_CT(psi);
  }
  if (csi->type_contribution == "CT" ||
      csi->type_contribution == "CJ" || 
      csi->type_contribution == "L2CT" ||
      csi->type_contribution == "L2CJ"){
    initialization_specific_CT(psi);
    if (switch_old_qT_version){initialization_specific_QT();}
  }

  // from "initialization.integration.cxx":
  
  ofstream out_maxevent;
  ofstream out_comparison;
  ofstream out_integration;
  ofstream out_result;
  ofstream out_gnuplot;
  ofstream out_proceeding;
  ofstream out_distribution;
  ofstream out_execution;

  logger << LOG_DEBUG << "before in_proceeding" << endl;

  int newstart = 0;
  ifstream in_proceeding(filename_proceeding.c_str());  
  if (in_proceeding){
    perform_proceeding_in(psi);
    logger << LOG_INFO << "int_end = " << int_end << endl;

    /*
    psi_i_gen = psi_i_acc + psi_i_rej;
    out_gnuplot.open(filename_gnuplot.c_str(), ofstream::out | ofstream::app);  
    out_gnuplot.close();
    perform_proceeding_check(psi);
    out_integration.open(filename_integration.c_str(), ofstream::out | ofstream::app);  
    */
    if (int_end == 2){
      int_end = 0;
      newstart = 1;
      system_execute(logger, "rm " + filename_proceeding);
    }
    else {
      //      psi_i_gen = psi_i_acc + psi_i_rej;
      out_gnuplot.open(filename_gnuplot.c_str(), ofstream::out | ofstream::app);  
      out_gnuplot.close();
      perform_proceeding_check(psi);
      out_integration.open(filename_integration.c_str(), ofstream::out | ofstream::app);  
      
      if (int_end == 1){
	/*	
	if (psi_switch_n_events_opt == 2){
	  int_end = 0;
	  out_integration << "******************************************************* resumption ******************************************************" << endl;
	}
	else {
*/
	  logger << LOG_INFO << "integration already finished" << endl; 
	  
	  if (switch_output_distribution){
	    output_distribution_complete(psi);
	    logger << LOG_INFO << "distribution output has been written anew." << endl; 
	  }
	  /*
	if (switch_output_integration){
	  output_step_integration(psi);
	  output_step_integration_TSV(psi);
	}
	*/
	  if (switch_output_result){
	    calculate_intermediate_result(psi);
	    output_step_result(psi);
	    output_step_result_TSV(psi);
	    logger << LOG_INFO << "result output has been written anew." << endl; 
	  }
	  
	  // needs to be updated !!! MC_....end_optimization need to be adapted !!!
	  // grid output should also be repeated !!!
	  /*
	  if (psi_switch_n_events_opt == 2){
	    logger << LOG_INFO << "psi_MC_opt_end = " << psi_MC_opt_end << endl;
	    if (psi_MC_opt_end > 0){psi.weight_output_MCweight_optimization();}
	    psi_random_manager.writeout_weights();
	    logger << LOG_INFO << "psi_tau_opt_end = " << psi_tau_opt_end << endl;
	    if (psi_tau_opt_end > 0){output_weight_vegas(psi_filename_tauweight, psi_tau_alpha);}
	    logger << LOG_INFO << "psi_x1x2_opt_end = " << psi_x1x2_opt_end << endl;
	    if (psi_x1x2_opt_end > 0){output_weight_vegas(psi_filename_x1x2weight, psi_x1x2_alpha);}
	    logger << LOG_INFO << "psi_z1z2_opt_end = " << psi_z1z2_opt_end << endl;
	    if (psi_z1z2_opt_end > 0){
	      for (int i_z = 1; i_z < 3; i_z++){
		output_weight_vegas(psi_filename_z1z2weight[i_z], psi_z1z2_alpha[i_z]);
	      }
	    }
	    logger << LOG_INFO << "grid output has been written anew." << endl; 
	  }
	  */
	  exit(0);
	  //	}
      }
      else {

	out_integration << "******************************************************* resumption ******************************************************" << endl;
	out_execution.open(filename_execution.c_str(), ofstream::out | ofstream::trunc);
	//	out_execution << "started" << endl;
	out_execution << "restarted" << endl;
	out_execution.close();
      }
    }
    out_integration.close();
  }
  else {
    //new
    newstart = 1;
  }
  if (newstart){
    //new
    out_integration.open(filename_integration.c_str(), ofstream::out | ofstream::trunc);  
    header_integration(out_integration, psi);
    out_integration.close();
    if (switch_TSV){
      for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	out_integration.open((filename_integration_TSV[i_s]).c_str(), ofstream::out | ofstream::trunc);
	out_result.open((filename_result_TSV[i_s]).c_str(), ofstream::out | ofstream::trunc);
	header_integration(out_integration, psi);
	out_integration.close();
	out_result.close();  
      }
    }
    out_execution.open(filename_execution.c_str(), ofstream::out | ofstream::trunc);
    out_execution.close();
    out_gnuplot.open(filename_gnuplot.c_str(), ofstream::out | ofstream::trunc);
    out_gnuplot.close();
    if (switch_CV){
      out_integration.open(filename_integration_CV.c_str(), ofstream::out | ofstream::trunc);
      header_integration(out_integration, psi);
      out_integration.close();      
    }

    out_result.open(filename_result.c_str(), ofstream::out | ofstream::trunc);  
    out_result.close();  
    out_maxevent.open(filename_maxevent.c_str(), ofstream::out | ofstream::trunc);  
    out_maxevent.close();  
 
  }
  in_proceeding.close();

  logger << LOG_DEBUG << "after in_proceeding" << endl;

  /*
  out_result.open(filename_result.c_str(), ofstream::out | ofstream::trunc);  
  out_result.close();  
  out_maxevent.open(filename_maxevent.c_str(), ofstream::out | ofstream::trunc);  
  out_maxevent.close();  
  */
  out_comparison.open(filename_comparison.c_str(), ofstream::out | ofstream::trunc);  
  out_comparison.close();
  
  logger << LOG_DEBUG << "finished" << endl;
}




void observable_set::initialization_runtime_partonlevel(){
  static Logger logger("observable_set::initialization_runtime_partonlevel");
  logger << LOG_DEBUG << "started" << endl;

  logger << LOG_DEBUG << "esi.object_list.size() = " << esi.object_list.size() << endl;

  /*
  ps_n_partonlevel.resize(n_ps, vector<int> (esi.object_list.size(), 0));

  determine_phasespace_object_partonlevel();
  for (int i_p = 1; i_p < esi.object_list.size(); i_p++){
    for (int i_a = 0; i_a < n_ps; i_a++){
      if (ps_n_partonlevel[i_a][i_p] > esi.pda[i_p].n_partonlevel){esi.pda[i_p].n_partonlevel = ps_n_partonlevel[i_a][i_p];}
    }
  }
  //  ps_relevant_n_partonlevel;
  */




  ps_runtime_jet_algorithm.resize(n_ps);
  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i_p = 3; i_p < csi->type_parton[i_a].size(); i_p++){
      for (int i_l = 0; i_l < jet_algorithm_list.size(); i_l++){
	if (csi->type_parton[i_a][i_p] == jet_algorithm_list[i_l]){
	  if (csi->type_parton[i_a][i_p] == 22){
	    photon_jet_algorithm = 1;
	    //	    if (!frixione_isolation){
	    // added to have photons in jet_algorithm (needs clarification, under which circumstances photons can become jets...)
	    //	    ps_runtime_jet_algorithm[i_a].push_back(i_p);
	    // removed again -> would lead to doubled photons in jet algorithm...
	    // would, however, be useful -> maybe adapt at the point where the protojets are selected (photon are added later by now).
	  }
	  else {ps_runtime_jet_algorithm[i_a].push_back(i_p);}
	}
      }
    }
    stringstream temp;
    for (int i_l = 0; i_l < ps_runtime_jet_algorithm[i_a].size(); i_l++){temp << setw(5) << ps_runtime_jet_algorithm[i_a][i_l];}
    logger << LOG_INFO << "ps_runtime_jet_algorithm[" << i_a << "].size() = " << ps_runtime_jet_algorithm[i_a].size() << " --- " << temp.str() << endl;
  }

  ps_runtime_photon_recombination.resize(n_ps);
  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i_p = 3; i_p < csi->type_parton[i_a].size(); i_p++){
      logger << LOG_DEBUG << "photon_recombination_list.size() = " << photon_recombination_list.size() << endl;
      for (int i_l = 0; i_l < photon_recombination_list.size(); i_l++){
	logger << LOG_DEBUG << "csi->type_parton[" << i_a << "][" << i_p << "] = " << csi->type_parton[i_a][i_p] << "   photon_recombination_list[" << i_l << "] = " << photon_recombination_list[i_l] << endl;
	if (csi->type_parton[i_a][i_p] == photon_recombination_list[i_l]){
	  ps_runtime_photon_recombination[i_a].push_back(i_p);
	}
      }
    }
    stringstream temp;
    for (int i_l = 0; i_l < ps_runtime_photon_recombination[i_a].size(); i_l++){temp << setw(5) << ps_runtime_photon_recombination[i_a][i_l];}
    logger << LOG_INFO << "ps_runtime_photon_recombination[" << i_a << "].size() = " << ps_runtime_photon_recombination[i_a].size() << " --- " << temp.str() << endl;
  }


  ps_runtime_photon.resize(n_ps);
  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i_p = 3; i_p < csi->type_parton[i_a].size(); i_p++){
      if (csi->type_parton[i_a][i_p] == 22){
	ps_runtime_photon[i_a].push_back(i_p);
      }
    }
    stringstream temp;
    for (int i_l = 0; i_l < ps_runtime_photon[i_a].size(); i_l++){temp << setw(5) << ps_runtime_photon[i_a][i_l];}
    logger << LOG_INFO << "ps_runtime_photon[" << i_a << "].size() = " << ps_runtime_photon[i_a].size() << " --- " << temp.str() << endl;
  }


  //  ps_runtime_original containts partons (type_parton numbers) that enter event selection with their (potentially dressed) parton-level momenta:
  ps_runtime_original.resize(n_ps);
  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i_p = 3; i_p < csi->type_parton[i_a].size(); i_p++){
      int flag = 0;
      for (int i_l = 0; i_l < ps_runtime_jet_algorithm[i_a].size(); i_l++){
	if (i_p == ps_runtime_jet_algorithm[i_a][i_l]){flag = 1; break;}
      }
      for (int i_l = 0; i_l < ps_runtime_photon[i_a].size(); i_l++){
	if (i_p == ps_runtime_photon[i_a][i_l]){flag = 1; break;}
      }

      if (!flag){
	ps_runtime_original[i_a].push_back(i_p);
      }
    }
    stringstream temp;
    for (int i_l = 0; i_l < ps_runtime_original[i_a].size(); i_l++){temp << setw(5) << ps_runtime_original[i_a][i_l];}
    logger << LOG_INFO << "ps_runtime_original[" << i_a << "].size() = " << ps_runtime_original[i_a].size() << " --- " << temp.str() << endl;
  }



  ps_runtime_missing.resize(n_ps);
  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i_p = 3; i_p < csi->type_parton[i_a].size(); i_p++){
      if (csi->type_parton[i_a][i_p] == 12 ||
	  csi->type_parton[i_a][i_p] == 14 ||
	  csi->type_parton[i_a][i_p] == 16 ||
	  csi->type_parton[i_a][i_p] == -12 ||
	  csi->type_parton[i_a][i_p] == -14 ||
	  csi->type_parton[i_a][i_p] == -16){
	ps_runtime_missing[i_a].push_back(i_p);
      }
    }
    stringstream temp;
    for (int i_l = 0; i_l < ps_runtime_missing[i_a].size(); i_l++){temp << setw(5) << ps_runtime_missing[i_a][i_l];}
    logger << LOG_INFO << "ps_runtime_missing[" << i_a << "].size() = " << ps_runtime_missing[i_a].size() << " --- " << temp.str() << endl;
  }


  logger.newLine(LOG_INFO);
  logger << LOG_INFO << setw(5) << "ps" << "   " << setw(15) << "photon" << "   " << setw(15) << "photon" << "   " << setw(15) << "jet" << "   " << setw(15) << "missing" << "   " << setw(15) << "original" << endl;
  logger << LOG_INFO << setw(5) << "" << "   " << setw(15) << "" << "   " << setw(15) << "recombination" << "   " << setw(15) << "algorithm" << "   " << setw(15) << "" << setw(15) << "" << endl;
  //  logger << LOG_INFO << setw(15) << "phasespace" << "   " << setw(25) << "photon" << "   " << setw(25) << "photon_recombination" << "   " << setw(25) << "jet_algorithm" << endl;
  for (int i_a = 0; i_a < n_ps; i_a++){
    stringstream temp_ph;
    for (int i_i = 0; i_i < ps_runtime_photon[i_a].size(); i_i++){
      temp_ph << ps_runtime_photon[i_a][i_i];
      if (i_i < ps_runtime_photon[i_a].size() - 1){temp_ph << "  ";}
    }
    stringstream temp_pr;
    for (int i_i = 0; i_i < ps_runtime_photon_recombination[i_a].size(); i_i++){
      temp_pr << ps_runtime_photon_recombination[i_a][i_i];
      if (i_i < ps_runtime_photon_recombination[i_a].size() - 1){temp_pr << "  ";}
    }
    stringstream temp_ja;
    for (int i_i = 0; i_i < ps_runtime_jet_algorithm[i_a].size(); i_i++){
      temp_ja << ps_runtime_jet_algorithm[i_a][i_i];
      if (i_i < ps_runtime_jet_algorithm[i_a].size() - 1){temp_ja << "  ";}
    }
    stringstream temp_mi;
    for (int i_i = 0; i_i < ps_runtime_missing[i_a].size(); i_i++){
      temp_mi << ps_runtime_missing[i_a][i_i];
      if (i_i < ps_runtime_missing[i_a].size() - 1){temp_mi << "  ";}
    }
    stringstream temp_og;
    for (int i_i = 0; i_i < ps_runtime_original[i_a].size(); i_i++){
      temp_og << ps_runtime_original[i_a][i_i];
      if (i_i < ps_runtime_original[i_a].size() - 1){temp_og << "  ";}
    }
    logger << LOG_INFO << setw(5) << i_a << "   " << setw(15) << temp_ph.str() << "   " << setw(15) << temp_pr.str() << "   " << setw(15) << temp_ja.str() << "   " << setw(15) << temp_mi.str() << "   " << setw(15) << temp_og.str() << endl;
  }
  logger.newLine(LOG_INFO);


  ps_n_partonlevel.resize(n_ps, vector<int> (esi.object_list.size(), 0));

  determine_phasespace_object_partonlevel();

  for (int i_p = 1; i_p < esi.object_list.size(); i_p++){
    for (int i_a = 0; i_a < n_ps; i_a++){
      if (ps_n_partonlevel[i_a][i_p] > esi.pda[i_p].n_partonlevel){esi.pda[i_p].n_partonlevel = ps_n_partonlevel[i_a][i_p];}
    }
  }


  logger << LOG_INFO << "Before removal of non-contributiong subprocesses:" << endl;
  logger.newLine(LOG_INFO);

  logger << LOG_INFO << "All objects:" << endl;
  logger.newLine(LOG_INFO);
  logger << LOG_INFO << setw(5) << right << "no" << setw(10) << "object" << setw(10) << "access" << setw(10) << "n_min" << setw(10) << "n_max" << setw(10) << "pT" << setw(10) << "ET" << setw(10) << "|eta|" << setw(10) << "|y|" << setw(5) << "" << setw(10) << "n_partonic"<< endl;
  for (int i = 0; i < esi.object_list.size(); i++){
    logger << LOG_INFO << setw(5) << right << i << setw(10) << esi.object_list[i] << setw(10) << access_object[esi.object_list[i]] << setw(10) << esi.pda[i].n_observed_min << setw(10) << esi.pda[i].n_observed_max << setw(10) << esi.pda[i].define_pT << setw(10) << esi.pda[i].define_ET << setw(10) << esi.pda[i].define_eta << setw(10) << esi.pda[i].define_y << setw(5) << "" << setw(10) << esi.pda[i].n_partonlevel << endl;
  }
  /*
  logger.newLine(LOG_INFO);
  logger << LOG_INFO << setw(5) << right << "no" << setw(10) << "object" << setw(10) << "n_part" << setw(5) << "" << "partons that can enter to respective object (phase-space dependent)" << endl;
  for (int i = 0; i < esi.object_list.size(); i++){
    stringstream temp;
    for (int i_a = 0; i_a < n_ps; i_a++){ 
      for (int i_j = 0; i_j < ps_runtime_order_inverse[i_a][i].size(); i_j++){
	temp << setw(1) << ps_runtime_order_inverse[i_a][i][i_j] << " ";
      }
      if (i_a < n_ps - 1){temp << "- ";}
    }
    logger << LOG_INFO << setw(5) << right << i << setw(10) << esi.object_list[i] << setw(10) << esi.pda[i].n_partonlevel << setw(5) << "" << temp.str() << endl;
  }
  logger.newLine(LOG_INFO);
  */


  // perform these checks before adapting n_observed_min etc. to n_partonlevel (if requirements cannot be fulfilled, runs won't start and write 0 as XS).

  int flag_remove = 0;
  vector<int> remove_phasespace(n_ps, 0);
  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i = 1; i < esi.object_list.size(); i++){
      if ((esi.pda[i].n_observed_min > ps_n_partonlevel[i_a][i]) && esi.object_category[i] != -1){
	logger << LOG_FATAL << esi.object_list[i] << ":   esi.pda.n_observed_min[" << i << "] = " << esi.pda[i].n_observed_min << " > " << ps_n_partonlevel[i_a][i] << " =   n_partonlevel[esi.observed_object[esi.object_list[" << i << "] = " << esi.object_list[i] << "] = " << i << "]" << endl;
	logger << LOG_FATAL << "Dipole " << i_a << " does not have a sufficient parton content !!!" << endl;

	//... check things like photons mimicking jets etc., which could let the dipole contribute (mainly photons treated as jets, collinear e+e- pairs as photons, etc.

	remove_phasespace[i_a] = 1;
	flag_remove = 1;
	break;
      }
    }
  }
  for (int i_a = 0; i_a < n_ps; i_a++){
    if (csi->type_contribution == "RA"){logger << LOG_INFO << "Phasespace " << setw(2) << i_a << ": " << setw(15) << (*RA_dipole)[i_a].name() << "   remove_phasespace = " << remove_phasespace[i_a] << endl;}
    else {logger << LOG_INFO << "Phasespace " << setw(2) << i_a << ":    remove_phasespace = " << remove_phasespace[i_a] << endl;}
  }
  if (accumulate(remove_phasespace.begin(), remove_phasespace.end(), 0) == remove_phasespace.size()){
    logger.newLine(LOG_FATAL);
    logger << LOG_FATAL << "No contribution can result from this subprocess due to parton content." << endl;
    logger.newLine(LOG_FATAL);
    int_end = 1;
  }

  if (flag_remove){
    // treatment of removed dipoles has to be added. By now, they are simply taken into account, which might cause problems !!!
    // Ideally, they should be removed completely (including phase-space generation)

    /*
    logger << LOG_FATAL << "(*RA_dipole).size() = " << (*RA_dipole).size() << endl;
    for (int i_a = n_ps - 1; i_a >= 0; i_a--){
      if (remove_dipole[i_a] == 1){(*RA_dipole).erase((*RA_dipole).begin() + i_a);}
    }
    logger << LOG_FATAL << "(*RA_dipole).size() = " << (*RA_dipole).size() << endl;
    //    initialization_RA((*RA_dipole));
    
    n_ps = (*RA_dipole).size();
    n_pc = (*RA_dipole).size();
    n_pz = 1;
    
    csi->type_parton.resize(n_ps);
    for (int i_a = 0; i_a < (*RA_dipole).size(); i_a++){
      csi->type_parton[i_a] = (*RA_dipole)[i_a].type_parton();
    }
    */
  }



  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i = 1; i < esi.object_list.size(); i++){
      int flag = 0;
      for (int i_j = 0; i_j < ps_runtime_order_inverse[i_a][i].size(); i_j++){
	for (int i_p = 0; i_p < ps_runtime_original[i_a].size(); i_p++){
	  if (ps_runtime_order_inverse[i_a][i][i_j] == ps_runtime_original[i_a][i_p]){
	    flag++;
	    break;
	  }
	}
	if (flag){break;}
      }
    
      if (flag){      
	if (esi.object_category[i] != -1 && esi.pda[i].n_observed_max < ps_n_partonlevel[i_a][i]){
	  if (esi.pda[i].define_pT == 0. && 
	      esi.pda[i].define_ET == 0. && 
	      esi.pda[i].define_y > 1.e10 &&
	      esi.pda[i].define_eta > 1.e10){
	    logger << LOG_INFO << "No event will survive event selection because of esi.object_list[" << i << "] = " << esi.object_list[i] << " -> Don't start integration." << endl;
	    int_end = 1;
	  }
	}
      }
    }
  }

  determine_object_partonlevel(); // needs modification !!!  User input should never be overwritten, even if contribution gives zero !!!

  determine_equivalent_object(); // Might need adaptation accordingly !!! Could be avoided by checking full object list for n_observed_min - n_partonlevel correspondence !!!
  determine_relevant_object();
  determine_runtime_object();


  logger << LOG_INFO << "After determine_runtime_object in initialization_runtime_partonlevel" << endl;

  logger.newLine(LOG_INFO);
  logger << LOG_INFO << setw(5) << right << "no" << setw(10) << "object" << setw(10) << "n_part" << setw(5) << "" << "partons that can enter to respective object (phase-space dependent)" << endl;
  for (int i = 0; i < esi.object_list.size(); i++){
    stringstream temp;
    for (int i_a = 0; i_a < n_ps; i_a++){ 
      for (int i_j = 0; i_j < ps_runtime_order_inverse[i_a][i].size(); i_j++){
	temp << setw(1) << ps_runtime_order_inverse[i_a][i][i_j] << " ";
      }
      if (i_a < n_ps - 1){temp << "- ";}
    }
    logger << LOG_INFO << setw(5) << right << i << setw(10) << esi.object_list[i] << setw(10) << esi.pda[i].n_partonlevel << setw(5) << "" << temp.str() << endl;
  }
  logger.newLine(LOG_INFO);


  /*
  // replaced by the procedure above !!!
  int flag_remove = 0;
  vector<int> remove_dipole(n_ps, 0);
  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i = 1; i < esi.object_list_selection.size(); i++){
      if ((esi.pds[i].n_observed_min > ps_n_partonlevel[i_a][esi.observed_object[esi.object_list_selection[i]]]) && esi.object_category[equivalent_no_object[i]] != -1){
	logger << LOG_FATAL << esi.object_list_selection[i] << ":   esi.pds.n_observed_min[" << i << "] = " << esi.pds[i].n_observed_min << " > " << ps_n_partonlevel[i_a][esi.observed_object[esi.object_list_selection[i]]] << " =   n_partonlevel[esi.observed_object[esi.object_list_selection[" << i << "] = " << esi.object_list_selection[i] << "] = " << esi.observed_object[esi.object_list_selection[i]] << "]" << endl;
	logger << LOG_FATAL << "Dipole " << i_a << " does not have a sufficient parton content !!!" << endl;
	//... check things like photons mimicking jets etc., which could let the dipole contribute (mainly photons treated as jets, collinear e+e- pairs as photons, etc.

	remove_dipole[i_a] = 1;
	flag_remove = 1;
	break;
      }
    }
  }
  */

 


  /*
  // replaced by above procedure !!!
  for (int i = 1; i < esi.object_list_selection.size(); i++){
    if ((esi.pds[i].n_observed_min > esi.pda[esi.observed_object[esi.object_list_selection[i]]].n_partonlevel) && esi.object_category[equivalent_no_object[i]] != -1){
      logger << LOG_FATAL << esi.object_list_selection[i] << ":   esi.pds.n_observed_min[" << i << "] = " << esi.pds[i].n_observed_min << " > " << esi.pds[esi.observed_object[esi.object_list_selection[i]]].n_partonlevel << " =   n_partonlevel[esi.observed_object[esi.object_list_selection[" << i << "] = " << esi.object_list_selection[i] << "] = " << esi.observed_object[esi.object_list_selection[i]] << "]" << endl;
      logger << LOG_FATAL << "No contribution can result from this subprocess - modification needed due to type conversion (e.g. jet/photon recombination alogorithm)" << endl;
      logger << LOG_FATAL << "Check if correct! Dipoles might have a different parton content, which is not covered by now !!!" << endl;
      int_end = 1;
      //      exit(0);
    }
  }
  */

  // Check if an upper limit on n_partonlevel can be realized (useful to technically veto e.g. *all* bjets in the massive case):
  for (int i_r = 1; i_r < esi.object_list_selection.size(); i_r++){
    int flag = 0;
    for (int i_j = 0; i_j < runtime_original.size(); i_j++){
      if (i_r == runtime_original[i_j]){flag++; break;}
    }
    logger << LOG_DEBUG << "esi.object_list_selection[" << i_r << "] = " << esi.object_list_selection[i_r] << " enters event selection with its original momentum." << endl;
    /*
    logger << LOG_DEBUG << "esi.object_list_selection[" << i_r << "] = " << esi.object_list_selection[i_r] << ": esi.pds.define_pT[i_r] = " << esi.pds.define_pT[i_r] << endl;
    logger << LOG_DEBUG << "esi.object_list_selection[" << i_r << "] = " << esi.object_list_selection[i_r] << ": esi.pds.define_ET[i_r] = " << esi.pds.define_ET[i_r] << endl;
    logger << LOG_DEBUG << "esi.object_list_selection[" << i_r << "] = " << esi.object_list_selection[i_r] << ": esi.pds.define_y[i_r] = " << esi.pds.define_y[i_r] << endl;
    logger << LOG_DEBUG << "esi.object_list_selection[" << i_r << "] = " << esi.object_list_selection[i_r] << ": esi.pds.define_eta[i_r] = " << esi.pds.define_eta[i_r] << endl;
    */
    if (esi.pds[i_r].n_observed_max < esi.pda[esi.observed_object[esi.object_list_selection[i_r]]].n_partonlevel && esi.object_category[equivalent_no_object[i_r]] != -1 &&
	//	if (esi.pds.n_observed_max[i_r] < n_partonlevel[esi.observed_object[esi.object_list_selection[i_r]]] &&
	esi.pds[i_r].define_pT == 0. &&
	esi.pds[i_r].define_ET == 0. &&
	esi.pds[i_r].define_y > 1.e10 &&
	esi.pds[i_r].define_eta > 1.e10){
      logger << LOG_DEBUG << "No event will survive event selection because of esi.object_list_selection[" << i_r << "] = " << esi.object_list_selection[i_r] << " -> Don't start integration." << endl;
      int_end = 1;
    }
  }



  for (int i_p = 0; i_p < esi.object_list.size(); i_p++){
    access_object[esi.object_list[i_p]] = equivalent_no_object[no_relevant_object[equivalent_object[esi.object_list[i_p]]]];
    logger << LOG_INFO << setw(5) << i_p << "   " << setw(10) << esi.object_list[i_p] << "   equivalent_object -> " << setw(5) << equivalent_object[esi.object_list[i_p]] << "   no_relevant_object - > " << no_relevant_object[equivalent_object[esi.object_list[i_p]]] << "   equivalent_no_object -> " << equivalent_no_object[no_relevant_object[equivalent_object[esi.object_list[i_p]]]] << endl;
  }

  output_initialization_runtime_partonlevel();


  for (int i_f = 0; i_f < esi.name_fiducial_cut.size(); i_f++){
    (esi.fiducial_cut).push_back(fiducialcut(esi.name_fiducial_cut[i_f], esi, *this));
  }
  for (int i_f = 0; i_f < esi.fiducial_cut.size(); i_f++){
    logger << LOG_INFO << "fiducial_cut[" << i_f << "] = " << esi.fiducial_cut[i_f].name << endl << esi.fiducial_cut[i_f] << endl;
  }



  
  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::output_initialization_runtime_partonlevel(){
  static Logger logger("observable_set::output_initialization_runtime_partonlevel");
  logger << LOG_DEBUG << "started" << endl;


  logger << LOG_INFO << "All objects:" << endl;
  logger.newLine(LOG_INFO);
  logger << LOG_INFO << setw(5) << right << "no" << setw(10) << "object" << setw(10) << "access" << setw(10) << "n_min" << setw(10) << "n_max" << setw(10) << "pT" << setw(10) << "ET" << setw(10) << "|eta|" << setw(10) << "|y|" << setw(5) << "" << setw(10) << "n_partonic"<< endl;
  for (int i = 0; i < esi.object_list.size(); i++){
    logger << LOG_INFO << setw(5) << right << i << setw(10) << esi.object_list[i] << setw(10) << access_object[esi.object_list[i]] << setw(10) << esi.pda[i].n_observed_min << setw(10) << esi.pda[i].n_observed_max << setw(10) << esi.pda[i].define_pT << setw(10) << esi.pda[i].define_ET << setw(10) << esi.pda[i].define_eta << setw(10) << esi.pda[i].define_y << setw(5) << "" << setw(10) << esi.pda[i].n_partonlevel << endl;
  }

  logger.newLine(LOG_INFO);

  logger << LOG_DEBUG << endl << "Access generic objects by name: " << endl;
  logger << LOG_DEBUG << setw(5) << right << "no" << setw(20) << "object_list[no]" << setw(30) << "access_object[.<--.]" << setw(20) << "object_list[no]" << setw(30) << "equivalent_object[.<--.]" << setw(30) << "no_relevant_object[.<--.]" << setw(30) << "equivalent_no_object[.<--.]" << endl;
  //"n_observed_min" << setw(20) << "n_observed_max" << setw(20) << "define_pT" << setw(20) << "define_ET" << setw(20) << "define_eta"  << setw(20) << "define_y"<< endl;
  for (int i = 1; i < esi.object_list.size(); i++){
    logger << LOG_DEBUG << setw(5) << right << i << setw(20) << esi.object_list[i] << setw(30) << access_object[esi.object_list[i]] << setw(20) << esi.object_list[i] << setw(30) << equivalent_object[esi.object_list[i]] << setw(30) << no_relevant_object[equivalent_object[esi.object_list[i]]] << setw(30) << equivalent_no_object[no_relevant_object[equivalent_object[esi.object_list[i]]]] << endl;
    //setw(20) << esi.pds[i].n_observed_min << setw(20) << esi.pds[i].n_observed_max << setw(20) << esi.pds[i].define_pT << setw(20) << esi.pds[i].define_ET << setw(20) << esi.pds[i].define_eta << setw(20) << esi.pds[i].define_y << endl;
  }
  logger.newLine(LOG_DEBUG);



  logger << LOG_INFO << "Relevant objects:" << endl;
  logger.newLine(LOG_INFO);
  logger << LOG_INFO << setw(5) << right << "no" << setw(10) << "object" << setw(10) << "access" << setw(10) << "n_min" << setw(10) << "n_max" << setw(10) << "pT" << setw(10) << "ET" << setw(10) << "|eta|"  << setw(10) << "|y|"<< endl;
  //  logger << LOG_INFO << setw(10) << "no" << setw(20) << "object"  << setw(20) << "access_object" << setw(20) << "n_observed_min" << setw(20) << "n_observed_max" << setw(20) << "define_pT" << setw(20) << "define_ET" << setw(20) << "define_eta"  << setw(20) << "define_y"<< endl;
  for (int i = 0; i < esi.object_list_selection.size(); i++){
    //    logger << LOG_INFO << setw(20) << esi.object_list_selection[i] << setw(20) << n_observed_min[esi.pds.object_category[i]] << setw(20) << n_observed_max[esi.pds.object_category[i]] << setw(20) << define_pT[esi.pds.object_category[i]] << setw(20) << define_eta[esi.pds.object_category[i]] << setw(20) << define_y[esi.pds.object_category[i]] << endl;
    logger << LOG_INFO << setw(5) << right << i << setw(10) << esi.object_list_selection[i] << setw(10) << access_object[esi.object_list_selection[i]] << setw(10) << esi.pds[i].n_observed_min << setw(10) << esi.pds[i].n_observed_max << setw(10) << esi.pds[i].define_pT << setw(10) << esi.pds[i].define_ET << setw(10) << esi.pds[i].define_eta << setw(10) << esi.pds[i].define_y << endl;
    //    logger << LOG_INFO << setw(10) << i << setw(20) << esi.object_list_selection[i] << setw(20) << access_object[esi.object_list_selection[i]] << setw(20) << esi.pds[i].n_observed_min << setw(20) << esi.pds[i].n_observed_max << setw(20) << esi.pds[i].define_pT << setw(20) << esi.pds[i].define_ET << setw(20) << esi.pds[i].define_eta << setw(20) << esi.pds[i].define_y << endl;
  }
  logger.newLine(LOG_INFO);

  logger << LOG_DEBUG << "Relevant objects:   runtime_jet_recombination   collects partons that enter jet recombination." << endl;
  logger.newLine(LOG_DEBUG);

  for (int i = 0; i < runtime_jet_recombination.size(); i++){
    logger << LOG_INFO << "runtime_jet_recombination[" << setw(3) << right << i << "] = " << setw(3) << runtime_jet_recombination[i] << "[" << setw(10) << esi.object_list_selection[runtime_jet_recombination[i]] << "]" << endl;
  }
  logger.newLine(LOG_INFO);
  
  logger << LOG_DEBUG << "Relevant objects:   runtime_photon_isolation   collects partons that could become isolated photons a la Frixione." << endl;
  logger.newLine(LOG_DEBUG);

  for (int i = 0; i < runtime_photon_isolation.size(); i++){
    logger << LOG_INFO << "runtime_photon_isolation[" << setw(3) << right << i << "] = " << setw(3) << runtime_photon_isolation[i] << "[" << setw(10) << esi.object_list_selection[runtime_photon_isolation[i]] << "]" << endl;
  }
  logger.newLine(LOG_INFO);
  
  /*
  //  runtime_photon_recombination is never used !!!
  logger << LOG_DEBUG << "Relevant objects:   runtime_photon_recombination   collects partons that enter photon recombination." << endl;
  logger.newLine(LOG_DEBUG);

  for (int i = 0; i < runtime_photon_recombination.size(); i++){
    logger << LOG_INFO << "runtime_photon_recombination[" << setw(3) << right << i << "] = " << setw(3) << runtime_photon_recombination[i] << "[" << setw(10) << esi.object_list_selection[runtime_photon_recombination[i]] << "]" << endl;
  }
  logger.newLine(LOG_INFO);
  */
  /*
  logger << LOG_DEBUG << "Relevant objects:   runtime_missing   collects partons that enter event selection with their parton-level momenta." << endl;
  logger.newLine(LOG_DEBUG);

  for (int i = 0; i < runtime_missing.size(); i++){
    logger << LOG_INFO << "runtime_missing[" << setw(3) << right << i << "] = " << setw(3) << runtime_missing[i] << "[" << setw(10) << esi.object_list_selection[runtime_missing[i]] << "]" << endl;
  }
  logger.newLine(LOG_INFO);
  */

  logger << LOG_DEBUG << "Relevant objects:   runtime_original   collects partons that enter event selection with their parton-level momenta." << endl;
  logger.newLine(LOG_DEBUG);

  for (int i = 0; i < runtime_original.size(); i++){
    logger << LOG_INFO << "runtime_original[" << setw(3) << right << i << "] = " << setw(3) << runtime_original[i] << "[" << setw(10) << esi.object_list_selection[runtime_original[i]] << "]" << endl;
  }
  logger.newLine(LOG_INFO);
  


  logger << LOG_DEBUG << "Relevant objects:   runtime_order   collects partons that belong to the respective object class." << endl;
  logger.newLine(LOG_DEBUG);

  for (int i = 0; i < esi.object_list_selection.size(); i++){
    stringstream sstemp;
    for (int j = 0; j < runtime_order[i].size(); j++){sstemp << setw(3) << runtime_order[i][j] << "[" << setw(6) << csi->type_parton[0][runtime_order[i][j]] << "]";}// << setw(3) << runtime_order[i][csi->type_parton[0][i]] << " ["
    logger << LOG_INFO << "runtime_order[" << setw(3) << right << i << " -> " << setw(10) << esi.object_list_selection[i] << "] -> " << sstemp.str() << endl;
  }
  logger.newLine(LOG_INFO);

  logger << LOG_DEBUG << "Relevant objects:   runtime_order_inverse   collects object classes a parton belongs to." << endl;
  logger.newLine(LOG_DEBUG);

  for (int i = 3; i < csi->type_parton[0].size(); i++){
    stringstream sstemp;
    for (int j = 0; j < runtime_order_inverse[i].size(); j++){sstemp << setw(10) << esi.object_list_selection[runtime_order_inverse[i][j]] << " [" << setw(3) << runtime_order_inverse[i][j] << "]";}
    logger << LOG_INFO << "runtime_order_inverse[csi->type_parton[" << setw(3) << right << i << "] = " << setw(5) << csi->type_parton[0][i] << "] -> " << sstemp.str() << endl;
  }
  logger.newLine(LOG_INFO);



  logger << LOG_DEBUG << endl << "Object definition after definition: " << endl;
  logger << LOG_DEBUG << setw(5) << "no" << setw(20) << "object" << setw(20) << "n_observed_min" << setw(20) << "n_observed_max" << setw(20) << "define_pT" << setw(20) << "define_ET" << setw(20) << "define_eta"  << setw(20) << "define_y"<< endl;
  for (int i = 1; i < esi.object_list.size(); i++){
    logger << LOG_DEBUG << setw(5) << i << setw(20) << esi.object_list[i] << setw(20) << esi.pda[i].n_observed_min << setw(20) << esi.pda[i].n_observed_max << setw(20) << esi.pda[i].define_pT << setw(20) << esi.pda[i].define_ET << setw(20) << esi.pda[i].define_eta << setw(20) << esi.pda[i].define_y << endl;
  }
  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << endl << "Relevant object definition after definition: " << endl;
  logger << LOG_DEBUG << setw(5) << "no" << setw(20) << "object" << setw(20) << "n_observed_min" << setw(20) << "n_observed_max" << setw(20) << "define_pT" << setw(20) << "define_ET" << setw(20) << "define_eta"  << setw(20) << "define_y"<< endl;
  for (int i = 1; i < esi.object_list_selection.size(); i++){
    logger << LOG_DEBUG << setw(5) << i << setw(20) << esi.object_list_selection[i] << setw(20) << esi.pds[i].n_observed_min << setw(20) << esi.pds[i].n_observed_max << setw(20) << esi.pds[i].define_pT << setw(20) << esi.pds[i].define_ET << setw(20) << esi.pds[i].define_eta << setw(20) << esi.pds[i].define_y << endl;
  }
  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG << endl << "Relation between relevant objects and generic objects: " << endl;
  logger << LOG_DEBUG << setw(5) << "no" << setw(25) << "esi.object_list_selection[no]" << " " << setw(47) << "no_relevant_object[esi.object_list_selection[no]]" << setw(25) << "equivalent_no_object[no]" << endl;
  //"n_observed_min" << setw(20) << "n_observed_max" << setw(20) << "define_pT" << setw(20) << "define_ET" << setw(20) << "define_eta"  << setw(20) << "define_y"<< endl;
  for (int i = 1; i < esi.object_list_selection.size(); i++){
    logger << LOG_DEBUG << setw(5) << i << setw(25) << esi.object_list_selection[i] << " " << setw(47) << no_relevant_object[esi.object_list_selection[i]] << setw(25) << equivalent_no_object[i] << endl;
    //setw(20) << esi.pds[i].n_observed_min << setw(20) << esi.pds[i].n_observed_max << setw(20) << esi.pds[i].define_pT << setw(20) << esi.pds[i].define_ET << setw(20) << esi.pds[i].define_eta << setw(20) << esi.pds[i].define_y << endl;
  }
  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG << "finished" << endl;
}
