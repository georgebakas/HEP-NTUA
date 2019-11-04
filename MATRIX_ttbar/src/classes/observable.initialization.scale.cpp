#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

// Relevant for TSV variation:

void observable_set::initialization_basic_TSV(inputparameter_set & isi){
  Logger logger("observable_set::initialization_basic_TSV (isi)");
  logger << LOG_DEBUG << "called" << endl;

  double empty_double = -1.;
  double empty_int = -1;

  switch_TSV = isi.switch_TSV;
  n_set_TSV = isi.n_set_TSV;
  switch_distribution_at_all_TSV = 0;

  if (!switch_TSV){n_set_TSV = 0;}
  if (!switch_TSV){logger << LOG_DEBUG << "switch_TSV = " << switch_TSV << " -> finished" << endl; return;}

  name_set_TSV = isi.name_set_TSV;



  // needed for scale-difference contributions:
  name_diff_set_TSV = isi.name_diff_set_TSV;
  n_diff_set_TSV = isi.n_diff_set_TSV;
  name_diff_set_plus_TSV = isi.name_diff_set_plus_TSV;
  name_diff_set_minus_TSV = isi.name_diff_set_minus_TSV;
  no_diff_set_plus_TSV = isi.no_diff_set_plus_TSV;
  no_diff_set_minus_TSV = isi.no_diff_set_minus_TSV;

  name_extended_set_TSV = isi.name_extended_set_TSV;
  n_extended_set_TSV = isi.n_extended_set_TSV;

  // needed for reference definition via TSV:
  name_reference_TSV = isi.name_reference_TSV;
  no_reference_TSV = isi.no_reference_TSV;
  no_scale_ren_reference_TSV = isi.no_scale_ren_reference_TSV;
  no_scale_fact_reference_TSV = isi.no_scale_fact_reference_TSV;
  no_qTcut_reference_TSV = isi.no_qTcut_reference_TSV;


  switch_distribution_TSV = isi.switch_distribution_TSV;
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    if (switch_distribution_TSV[i_s] > 0){switch_distribution_at_all_TSV = 1; break;}
  }
  
  switch_moment_TSV = isi.switch_moment_TSV;

  max_n_integrand_TSV = isi.max_n_integrand_TSV;

  // probably not used:
  min_qTcut_TSV = isi.min_qTcut_TSV;
  max_qTcut_TSV = isi.max_qTcut_TSV;
  min_qTcut_distribution_TSV = isi.min_qTcut_distribution_TSV;
  max_qTcut_distribution_TSV = isi.max_qTcut_distribution_TSV;

  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    if (min_qTcut_TSV[i_s] == empty_double){min_qTcut_TSV[i_s] = min_qTcut;}
    if (max_qTcut_TSV[i_s] == empty_double){max_qTcut_TSV[i_s] = min_qTcut + n_qTcut * step_qTcut;}
    if (min_qTcut_distribution_TSV[i_s] == empty_double){min_qTcut_distribution_TSV[i_s] = min_qTcut;}
    if (max_qTcut_distribution_TSV[i_s] == empty_double){max_qTcut_distribution_TSV[i_s] = min_qTcut + n_qTcut * step_qTcut;}
  }


    logger << LOG_DEBUG_VERBOSE << "SSSSS   n_extended_set_TSV = " << n_extended_set_TSV << endl;
  // also needed for scale-difference contributions:
    // switch_distribution may not be set to 0 in summary routine !!!
  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
    logger << LOG_DEBUG_VERBOSE << "SSSSS   switch_distribution_TSV[" << i_s << "] = " << switch_distribution_TSV[i_s] << endl;
    if (switch_distribution_TSV[i_s] == empty_int){switch_distribution_TSV[i_s] = switch_distribution;}
    if (switch_moment_TSV[i_s] == empty_int){switch_moment_TSV[i_s] = n_moments;}

    if (switch_distribution_TSV[i_s] == 0){max_n_integrand_TSV[i_s] = 1;}
    else {max_n_integrand_TSV[i_s] = 3;}

    logger << LOG_DEBUG_VERBOSE << "SSSSS   switch_distribution_TSV[" << i_s << "] = " << switch_distribution_TSV[i_s] << endl;
  }

  // initialization of scale-related quantities

  // central_scale_TSV
  // central_scale_ren_TSV
  // central_scale_fact_TSV
  // relative_central_scale_TSV
  // relative_central_scale_ren_TSV
  // relative_central_scale_fact_TSV
  // n_scale_TSV
  // n_scale_ren_TSV
  // n_scale_fact_TSV
  // factor_scale_TSV
  // factor_scale_ren_TSV
  // factor_scale_fact_TSV
  // dynamic_scale_TSV
  // dynamic_scale_ren_TSV
  // dynamic_scale_fact_TSV
  // no_central_scale_ren_TSV
  // no_central_scale_fact_TSV

  // max_n_scale_ren_TSV
  // max_n_scale_fact_TSV

  // relative_scale_ren_TSV
  // relative_scale_fact_TSV

  // n_dynamic_scale_TSV
  // switch_dynamic_scale_TSV

  central_scale_TSV = isi.central_scale_TSV;
  central_scale_ren_TSV = isi.central_scale_ren_TSV;
  central_scale_fact_TSV = isi.central_scale_fact_TSV;

  relative_central_scale_TSV = isi.relative_central_scale_TSV;
  relative_central_scale_ren_TSV = isi.relative_central_scale_ren_TSV;
  relative_central_scale_fact_TSV = isi.relative_central_scale_fact_TSV;

  n_scale_TSV = isi.n_scale_TSV;
  n_scale_ren_TSV = isi.n_scale_ren_TSV;
  n_scale_fact_TSV = isi.n_scale_fact_TSV;

  factor_scale_TSV = isi.factor_scale_TSV;
  factor_scale_ren_TSV = isi.factor_scale_ren_TSV;
  factor_scale_fact_TSV = isi.factor_scale_fact_TSV;
  
  dynamic_scale_TSV = isi.dynamic_scale_TSV;
  dynamic_scale_ren_TSV = isi.dynamic_scale_ren_TSV;
  dynamic_scale_fact_TSV = isi.dynamic_scale_fact_TSV;
  
  no_central_scale_ren_TSV = isi.no_central_scale_ren_TSV;
  no_central_scale_fact_TSV = isi.no_central_scale_fact_TSV;

  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    if (n_scale_ren_TSV[i_s] == empty_int){n_scale_ren_TSV[i_s] = n_scale_TSV[i_s];}
    if (n_scale_ren_TSV[i_s] == empty_int){logger << LOG_FATAL << "n_scale_ren_TSV[" << i_s << "] has not been defined!" << endl; exit(1);}
    if (n_scale_ren_TSV[i_s] % 2 == 0){logger << LOG_FATAL << "n_scale_ren_TSV[" << i_s << "] must be an odd number!" << endl; exit(1);}
    no_central_scale_ren_TSV[i_s] = (n_scale_ren_TSV[i_s] - 1) / 2;

    if (n_scale_fact_TSV[i_s] == empty_int){n_scale_fact_TSV[i_s] = n_scale_TSV[i_s];}
    if (n_scale_fact_TSV[i_s] == empty_int){logger << LOG_FATAL << "n_scale_fact_TSV[" << i_s << "] has not been defined!" << endl; exit(1);}
    if (n_scale_fact_TSV[i_s] % 2 == 0){logger << LOG_FATAL << "n_scale_fact_TSV[" << i_s << "] must be an odd number!" << endl; exit(1);}
    no_central_scale_fact_TSV[i_s] = (n_scale_fact_TSV[i_s] - 1) / 2;
  }



  if (name_reference_TSV == ""){
    logger << LOG_INFO << "No reference scale defined. Standard value used." << endl;
    logger << LOG_INFO << "name_extended_set_TSV.size() = " << name_extended_set_TSV.size() << endl;
    int x_s = 0;
    name_reference_TSV = name_extended_set_TSV[x_s];
    no_reference_TSV = x_s;
  }
  if (no_scale_ren_reference_TSV == -1){no_scale_ren_reference_TSV = no_central_scale_ren_TSV[no_reference_TSV];}
  if (no_scale_fact_reference_TSV == -1){no_scale_fact_reference_TSV = no_central_scale_fact_TSV[no_reference_TSV];}
  if (no_qTcut_reference_TSV == -1){no_qTcut_reference_TSV = 0;}

  logger << LOG_INFO << "default output from TSV:" << endl;
  logger << LOG_INFO << "name_reference_TSV          = " << name_reference_TSV << endl;
  logger << LOG_INFO << "no_reference_TSV            = " << no_reference_TSV << endl;
  logger << LOG_INFO << "no_scale_ren_reference_TSV  = " << no_scale_ren_reference_TSV << endl;
  logger << LOG_INFO << "no_scale_fact_reference_TSV = " << no_scale_fact_reference_TSV << endl;
  logger << LOG_INFO << "no_qTcut_reference_TSV      = " << no_qTcut_reference_TSV << endl;
    logger.newLine(LOG_INFO);





  // set  n_scale_TSV / n_scale_ren_TSV / n_scale_fact_TSV  for scale-difference contributions (and exit if not consistent):
  for (int i_s = 0; i_s < n_diff_set_TSV; i_s++){
    int i_es = n_set_TSV + i_s;
    if (n_scale_TSV[no_diff_set_plus_TSV[i_s]] == n_scale_TSV[no_diff_set_minus_TSV[i_s]]){n_scale_TSV[i_es] = n_scale_TSV[no_diff_set_plus_TSV[i_s]];}
    else {logger << LOG_FATAL << "n_scale_TSV of " << name_diff_set_TSV[i_s] << " is not well-defined: n_scale_TSV[" << no_diff_set_plus_TSV[i_s] << "] = " << n_scale_TSV[no_diff_set_plus_TSV[i_s]] << " =/= " << n_scale_TSV[no_diff_set_minus_TSV[i_s]] << " = n_scale_TSV[" << no_diff_set_minus_TSV[i_s] << "]" << endl; exit(1);}

    if (n_scale_ren_TSV[no_diff_set_plus_TSV[i_s]] == n_scale_ren_TSV[no_diff_set_minus_TSV[i_s]]){n_scale_ren_TSV[i_es] = n_scale_ren_TSV[no_diff_set_plus_TSV[i_s]];}
    else {logger << LOG_FATAL << "n_scale_ren_TSV of " << name_diff_set_TSV[i_s] << " is not well-defined: n_scale_ren_TSV[" << no_diff_set_plus_TSV[i_s] << "] = " << n_scale_ren_TSV[no_diff_set_plus_TSV[i_s]] << " =/= " << n_scale_ren_TSV[no_diff_set_minus_TSV[i_s]] << " = n_scale_ren_TSV[" << no_diff_set_minus_TSV[i_s] << "]" << endl; exit(1);}

    if (n_scale_fact_TSV[no_diff_set_plus_TSV[i_s]] == n_scale_fact_TSV[no_diff_set_minus_TSV[i_s]]){n_scale_fact_TSV[i_es] = n_scale_fact_TSV[no_diff_set_plus_TSV[i_s]];}
    else {logger << LOG_FATAL << "n_scale_fact_TSV of " << name_diff_set_TSV[i_s] << " is not well-defined: n_scale_fact_TSV[" << no_diff_set_plus_TSV[i_s] << "] = " << n_scale_fact_TSV[no_diff_set_plus_TSV[i_s]] << " =/= " << n_scale_fact_TSV[no_diff_set_minus_TSV[i_s]] << " = n_scale_fact_TSV[" << no_diff_set_minus_TSV[i_s] << "]" << endl; exit(1);}
  }

  // reason for max_n_scale... not clear !!! probably not used...
  max_n_scale_ren_TSV = 0;
  max_n_scale_fact_TSV = 0;
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    if (n_scale_ren_TSV[i_s] > max_n_scale_ren_TSV){max_n_scale_ren_TSV = n_scale_ren_TSV[i_s];}
    if (n_scale_fact_TSV[i_s] > max_n_scale_fact_TSV){max_n_scale_fact_TSV = n_scale_fact_TSV[i_s];}
  }
  logger << LOG_DEBUG << "max_n_scale_ren_TSV = " << max_n_scale_ren_TSV << endl;
  logger << LOG_DEBUG << "max_n_scale_fact_TSV = " << max_n_scale_fact_TSV << endl;


  // also needed for scale-difference contributions:
  relative_scale_ren_TSV.resize(n_extended_set_TSV, vector<double> (max_n_scale_ren_TSV, 0.));
  relative_scale_fact_TSV.resize(n_extended_set_TSV, vector<double> (max_n_scale_fact_TSV, 0.));


  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    if (dynamic_scale_ren_TSV[i_s] == empty_int){dynamic_scale_ren_TSV[i_s] = dynamic_scale_TSV[i_s];}
    if (dynamic_scale_ren_TSV[i_s] == empty_int){
      dynamic_scale_ren_TSV[i_s] = 0;
      logger << LOG_WARN << "dynamic_scale_ren_TSV[" << i_s << "] has been set to default value (0) !" << endl;
    }

    if (dynamic_scale_fact_TSV[i_s] == empty_int){dynamic_scale_fact_TSV[i_s] = dynamic_scale_TSV[i_s];}
    if (dynamic_scale_fact_TSV[i_s] == empty_int){
      dynamic_scale_fact_TSV[i_s] = 0;
      logger << LOG_WARN << "dynamic_scale_fact_TSV[" << i_s << "] has been set to default value (0) !" << endl;
    }

    if (relative_central_scale_ren_TSV[i_s] == empty_double){relative_central_scale_ren_TSV[i_s] = relative_central_scale_TSV[i_s];}
    if (central_scale_ren_TSV[i_s] == empty_double){central_scale_ren_TSV[i_s] = central_scale_TSV[i_s];}

    if (relative_central_scale_ren_TSV[i_s] == empty_double && central_scale_ren_TSV[i_s] == empty_double){
      logger << LOG_FATAL << "relative_central_scale_ren_TSV[" << i_s << "] and central_scale_ren_TSV[" << i_s << "] have not been defined!" << endl; exit(1);
    }
    else if (relative_central_scale_ren_TSV[i_s] == empty_double && central_scale_ren_TSV[i_s] != empty_double){
      relative_central_scale_ren_TSV[i_s] = central_scale_ren_TSV[i_s] / (scale_ren / prefactor_reference);
      // !!!     relative_central_scale_ren_TSV[i_s] = central_scale_ren_TSV[i_s] / scale_ren;
    }
    else if (relative_central_scale_ren_TSV[i_s] != empty_double && central_scale_ren_TSV[i_s] == empty_double){
      central_scale_ren_TSV[i_s] = relative_central_scale_ren_TSV[i_s] * (scale_ren / prefactor_reference);
      // !!!      central_scale_ren_TSV[i_s] = relative_central_scale_ren_TSV[i_s] * scale_ren;
    }
    else if (relative_central_scale_ren_TSV[i_s] != empty_double && central_scale_ren_TSV[i_s] != empty_double){
      logger << LOG_FATAL << "Both alternative definitions relative_central_scale_ren_TSV[" << i_s << "] and central_scale_ren_TSV[" << i_s << "] have been defined!" << endl; exit(1);
    }


    if (relative_central_scale_fact_TSV[i_s] == empty_double){relative_central_scale_fact_TSV[i_s] = relative_central_scale_TSV[i_s];}
    if (central_scale_fact_TSV[i_s] == empty_double){central_scale_fact_TSV[i_s] = central_scale_TSV[i_s];}

    if (relative_central_scale_fact_TSV[i_s] == empty_double && central_scale_fact_TSV[i_s] == empty_double){logger << LOG_FATAL << "relative_central_scale_fact_TSV[" << i_s << "] and central_scale_fact_TSV[" << i_s << "] have not been defined!" << endl; exit(1);}
    else if (relative_central_scale_fact_TSV[i_s] == empty_double && central_scale_fact_TSV[i_s] != empty_double){
      relative_central_scale_fact_TSV[i_s] = central_scale_fact_TSV[i_s] / (scale_fact / prefactor_reference);
      // !!!      relative_central_scale_fact_TSV[i_s] = central_scale_fact_TSV[i_s] / scale_fact;
    }
    else if (relative_central_scale_fact_TSV[i_s] != empty_double && central_scale_fact_TSV[i_s] == empty_double){
      central_scale_fact_TSV[i_s] = relative_central_scale_fact_TSV[i_s] * (scale_fact / prefactor_reference);
      // !!!      central_scale_fact_TSV[i_s] = relative_central_scale_fact_TSV[i_s] * scale_fact;
    }
    else if (relative_central_scale_fact_TSV[i_s] != empty_double && central_scale_fact_TSV[i_s] != empty_double){
      logger << LOG_FATAL << "Both alternative definitions relative_central_scale_fact_TSV[" << i_s << "] and central_scale_fact_TSV[" << i_s << "] have been defined!" << endl; exit(1);
    }

    if (factor_scale_ren_TSV[i_s] == empty_int){factor_scale_ren_TSV[i_s] = factor_scale_TSV[i_s];}
    if (factor_scale_ren_TSV[i_s] == empty_int){logger << LOG_FATAL << "factor_scale_ren_TSV[" << i_s << "] has not been defined!" << endl; exit(1);}

    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
      if (n_scale_ren_TSV[i_s] > 1){relative_scale_ren_TSV[i_s][i_r] = pow(10., log10(double(factor_scale_ren_TSV[i_s])) * double(2. * i_r / (n_scale_ren_TSV[i_s] - 1) - 1));}
      else {relative_scale_ren_TSV[i_s][i_r] = 1.;}
    }

    if (factor_scale_fact_TSV[i_s] == empty_int){factor_scale_fact_TSV[i_s] = factor_scale_TSV[i_s];}
    if (factor_scale_fact_TSV[i_s] == empty_int){logger << LOG_FATAL << "factor_scale_fact_TSV[" << i_s << "] has not been defined!" << endl; exit(1);}
    for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
      if (n_scale_fact_TSV[i_s] > 1){relative_scale_fact_TSV[i_s][i_f] = pow(10., log10(double(factor_scale_fact_TSV[i_s])) * double(2. * i_f / (n_scale_fact_TSV[i_s] - 1) - 1));}
      else {relative_scale_fact_TSV[i_s][i_f] = 1.;}
    }

    if (dynamic_scale_ren_TSV[i_s] == empty_int){dynamic_scale_ren_TSV[i_s] = dynamic_scale_TSV[i_s];}
    if (dynamic_scale_ren_TSV[i_s] == empty_int){logger << LOG_FATAL << "dynamic_scale_ren_TSV[" << i_s << "] has not been defined!" << endl; exit(1);}

    if (dynamic_scale_fact_TSV[i_s] == empty_int){dynamic_scale_fact_TSV[i_s] = dynamic_scale_TSV[i_s];}
    if (dynamic_scale_fact_TSV[i_s] == empty_int){logger << LOG_FATAL << "dynamic_scale_fact_TSV[" << i_s << "] has not been defined!" << endl; exit(1);}

    if (n_scale_ren_TSV[i_s] > max_n_scale_ren_TSV){max_n_scale_ren_TSV = n_scale_ren_TSV[i_s];}
    if (n_scale_fact_TSV[i_s] > max_n_scale_fact_TSV){max_n_scale_fact_TSV = n_scale_fact_TSV[i_s];}
  }



  // also needed for scale-difference contributions (check if identity for vector<double> comparison works here !!!):
  
  for (int i_s = 0; i_s < n_diff_set_TSV; i_s++){
    int i_es = n_set_TSV + i_s;
    if (relative_scale_ren_TSV[no_diff_set_plus_TSV[i_s]] == relative_scale_ren_TSV[no_diff_set_minus_TSV[i_s]]){relative_scale_ren_TSV[i_es] = relative_scale_ren_TSV[no_diff_set_plus_TSV[i_s]];}
    else {
      for (int i_r = 0; i_r < relative_scale_ren_TSV[i_s].size(); i_r++){
	logger << LOG_FATAL << "relative_scale_ren_TSV of " << name_diff_set_TSV[i_s] << " is not well-defined: " 
	       << "relative_scale_ren_TSV[" << no_diff_set_plus_TSV[i_s] << "][" << i_r << "] = " << relative_scale_ren_TSV[no_diff_set_plus_TSV[i_s]][i_r] << " =/= " 
	       << relative_scale_ren_TSV[no_diff_set_minus_TSV[i_s]][i_r] << " = relative_scale_ren_TSV[" << no_diff_set_minus_TSV[i_s] << "][" << i_r << "]" << endl; 
      }
      exit(1);
    }
  }
  
  for (int i_s = 0; i_s < n_diff_set_TSV; i_s++){
    int i_es = n_set_TSV + i_s;
    if (relative_scale_fact_TSV[no_diff_set_plus_TSV[i_s]] == relative_scale_fact_TSV[no_diff_set_minus_TSV[i_s]]){relative_scale_fact_TSV[i_es] = relative_scale_fact_TSV[no_diff_set_plus_TSV[i_s]];}
    else {
      for (int i_f = 0; i_f < relative_scale_fact_TSV[i_s].size(); i_f++){
	logger << LOG_FATAL << "relative_scale_fact_TSV of " << name_diff_set_TSV[i_s] << " is not well-defined: " 
	       << "relative_scale_fact_TSV[" << no_diff_set_plus_TSV[i_s] << "][" << i_f << "] = " << relative_scale_fact_TSV[no_diff_set_plus_TSV[i_s]][i_f] << " =/= " 
	       << relative_scale_fact_TSV[no_diff_set_minus_TSV[i_s]][i_f] << " = relative_scale_fact_TSV[" << no_diff_set_minus_TSV[i_s] << "][" << i_f << "]" << endl; 
      }
      exit(1);
    }
  }
  

  n_dynamic_scale_TSV = 1;
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    if (dynamic_scale_ren_TSV[i_s] > n_dynamic_scale_TSV + 1){n_dynamic_scale_TSV = dynamic_scale_ren_TSV[i_s] + 1;}
    if (dynamic_scale_fact_TSV[i_s] > n_dynamic_scale_TSV + 1){n_dynamic_scale_TSV = dynamic_scale_fact_TSV[i_s] + 1;}
  }

  switch_dynamic_scale_TSV.resize(n_dynamic_scale_TSV, 0);
  for (int i_d = 0; i_d < n_dynamic_scale_TSV; i_d++){
    for (int i_s = 0; i_s < n_set_TSV; i_s++){
      if (dynamic_scale_ren_TSV[i_s] == i_d || dynamic_scale_fact_TSV[i_s] == i_d){
	switch_dynamic_scale_TSV[i_d] = 1;
	break;
      }
    }
  }

  // initialization of scale-related quantities done



  logger << LOG_INFO << "The following dynamic scales are used:" << endl;
  for (int i_d = 0; i_d < n_dynamic_scale_TSV; i_d++){
    if (switch_dynamic_scale_TSV[i_d] == 1){
      logger << LOG_INFO << "switch_dynamic_scale_TSV[" << i_d << "] = " << switch_dynamic_scale_TSV[i_d] << endl;
    }
  }
  logger.newLine(LOG_INFO); 

  logger << LOG_INFO << "The following scalesets have been defined:" << endl;
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    logger << LOG_INFO << endl;
    logger << LOG_INFO << "scaleset " << i_s << ":" << endl;
    logger << LOG_INFO << setw(35) << "name_set_TSV" << " = " << name_set_TSV[i_s] << endl;
    logger << LOG_INFO << setw(35) << "dynamic_scale_ren_TSV" << " = " << dynamic_scale_ren_TSV[i_s] << endl;
    logger << LOG_INFO << setw(35) << "relative_central_scale_ren_TSV" << " = " << relative_central_scale_ren_TSV[i_s] << endl;
    if (dynamic_scale_ren_TSV[i_s] == 0){
      logger << LOG_INFO << setw(35) << "central_scale_ren_TSV" << " = " << central_scale_ren_TSV[i_s] << endl;
    }
    logger << LOG_INFO << setw(35) << "n_scale_ren_TSV" << " = " << n_scale_ren_TSV[i_s] << endl;
    logger << LOG_INFO << setw(35) << "factor_scale_ren_TSV" << " = " << factor_scale_ren_TSV[i_s] << endl;
    logger << LOG_INFO << setw(35) << "no_central_scale_ren_TSV" << " = " << no_central_scale_ren_TSV[i_s] << endl;
    logger << LOG_INFO << setw(35) << "dynamic_scale_fact_TSV" << " = " << dynamic_scale_fact_TSV[i_s] << endl;
    logger << LOG_INFO << setw(35) << "relative_central_scale_fact_TSV" << " = " << relative_central_scale_fact_TSV[i_s] << endl;
    if (dynamic_scale_fact_TSV[i_s] == 0){
      logger << LOG_INFO << setw(35) << "central_scale_fact_TSV" << " = " << central_scale_fact_TSV[i_s] << endl;
    }
    logger << LOG_INFO << setw(35) << "n_scale_fact_TSV" << " = " << n_scale_fact_TSV[i_s] << endl;
    logger << LOG_INFO << setw(35) << "factor_scale_fact_TSV" << " = " << factor_scale_fact_TSV[i_s] << endl;
    logger << LOG_INFO << setw(35) << "no_central_scale_fact_TSV" << " = " << no_central_scale_fact_TSV[i_s] << endl;
    logger.newLine(LOG_INFO); 
    logger << LOG_INFO << setw(35) << "min_qTcut_TSV" << " = " << min_qTcut_TSV[i_s] << endl;
    logger << LOG_INFO << setw(35) << "max_qTcut_TSV" << " = " << max_qTcut_TSV[i_s] << endl;
    logger << LOG_INFO << setw(35) << "switch_distribution_TSV" << " = " << switch_distribution_TSV[i_s] << endl;
    if (switch_distribution_TSV[i_s] != 0){
      logger << LOG_INFO << setw(35) << "min_qTcut_distribution_TSV" << " = " << min_qTcut_distribution_TSV[i_s] << endl;
      logger << LOG_INFO << setw(35) << "max_qTcut_distribution_TSV" << " = " << max_qTcut_distribution_TSV[i_s] << endl;
    }
    logger << LOG_INFO << setw(35) << "switch_moment_TSV" << " = " << switch_moment_TSV[i_s] << endl;
  }
  

  logger << LOG_DEBUG << "finished" << endl;
}



void observable_set::initialization_TSV(){
  Logger logger("observable_set::initialization_TSV");
  logger << LOG_DEBUG << "started" << endl;

  // Should actually be shifted elsewhere (not only needed in TSV ???):
  QCD_order = contribution_order_alpha_s;
  // !!!    alpha_S = LHAPDF::alphasPDF(scale_ren / prefactor_reference);
  alpha_S = LHAPDF::alphasPDF(scale_ren);
  value_ME2term.resize(n_pc, 0.); // needed also without TSV !!!
  
  if (!switch_TSV){logger << LOG_DEBUG << "switch_TSV = " << switch_TSV << " -> finished" << endl; return;}
  
  logger << LOG_DEBUG << "define_scale_calculation_TSV started" << endl;

  
  needed_scale2_ren = 0;
  needed_scale2_fact = 0;
  if (csi->class_contribution_CS_collinear ||
      csi->class_contribution_QT_collinear ||
      csi->class_contribution_NJ_collinear ||
      csi->class_contribution_QT_virtual ||
      csi->class_contribution_NJ_virtual ||
      csi->type_contribution == "CT2" ||
      csi->type_contribution == "VT2" ||
      csi->type_contribution == "CJ2" ||
      csi->type_contribution == "VJ2"){
    /*
  if (type_contribution == "CA" ||
      type_contribution == "RCA" ||
      type_contribution == "RCJ" ||
      type_contribution == "L2CA" ||
      type_contribution == "CT" ||
      type_contribution == "CT2" ||
      type_contribution == "L2CT" ||
      type_contribution == "VT" ||
      type_contribution == "VT2" ||
      type_contribution == "L2VT" ||
      type_contribution == "CJ"  ||
      type_contribution == "CJ2" ||
      type_contribution == "L2CJ" ||
      type_contribution == "VJ"  ||
      type_contribution == "VJ2" ||
      type_contribution == "L2VJ"){
    */
    needed_scale2_fact = 1;
  }
  

  // initialization of all renormalization-scale connected variables

  // max_dyn_ren
  // dynamic_scale_ren_TSV
  // value_relative_scale_ren
  // value_relative_scale2_ren
  // value_scale_ren
  // value_scale2_ren
  // no_value_ren_TSV

  max_dyn_ren = 0;
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    if (dynamic_scale_ren_TSV[i_s] > max_dyn_ren){
      max_dyn_ren = dynamic_scale_ren_TSV[i_s];
    }
  }
  n_scale_dyn_ren.resize(max_dyn_ren + 1, 0);
  ///  n_scale_dyn_ren.resize(max_dyn_ren + 1);
  
  value_relative_scale_ren.resize(max_dyn_ren + 1);
  value_relative_scale2_ren.resize(max_dyn_ren + 1);
  
  value_alpha_S_TSV.resize(n_ps);
  value_alpha_S_TSV[0].resize(max_dyn_ren + 1);
  
  value_relative_factor_alpha_S.resize(n_ps);
  value_relative_factor_alpha_S[0].resize(max_dyn_ren + 1);

  value_scale_ren.resize(n_ps);
  value_scale_ren[0].resize(max_dyn_ren + 1);
  
  value_scale2_ren.resize(n_ps);
  value_scale2_ren[0].resize(max_dyn_ren + 1);
  
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    int v_sr = dynamic_scale_ren_TSV[i_s];
    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
      int flag = 0;
      for (int v_xr = 0; v_xr < value_relative_scale_ren[v_sr].size(); v_xr++){
	if (abs(value_relative_scale_ren[v_sr][v_xr] - relative_central_scale_ren_TSV[i_s] * relative_scale_ren_TSV[i_s][i_r]) / value_relative_scale_ren[v_sr][v_xr] < 1.e-15){
	  flag = 1; 
	  break;
	}
      }
      if (flag == 0){
	value_relative_scale_ren[v_sr].push_back(relative_central_scale_ren_TSV[i_s] * relative_scale_ren_TSV[i_s][i_r]);
	value_relative_scale2_ren[v_sr].push_back(pow(relative_central_scale_ren_TSV[i_s] * relative_scale_ren_TSV[i_s][i_r], 2));
	if (v_sr == 0){
	  value_alpha_S_TSV[0][v_sr].push_back(LHAPDF::alphasPDF(central_scale_ren_TSV[i_s] * relative_scale_ren_TSV[i_s][i_r]));
	  value_relative_factor_alpha_S[0][v_sr].push_back(pow(LHAPDF::alphasPDF(central_scale_ren_TSV[i_s] * relative_scale_ren_TSV[i_s][i_r]) / LHAPDF::alphasPDF(scale_ren), contribution_order_alpha_s));
	  // !!!     value_relative_factor_alpha_S[0][v_sr].push_back(pow(LHAPDF::alphasPDF(central_scale_ren_TSV[i_s] * relative_scale_ren_TSV[i_s][i_r]) / LHAPDF::alphasPDF(scale_ren / prefactor_reference), contribution_order_alpha_s));
	  value_scale_ren[0][v_sr].push_back(central_scale_ren_TSV[i_s] * relative_scale_ren_TSV[i_s][i_r]);
	  value_scale2_ren[0][v_sr].push_back(pow(central_scale_ren_TSV[i_s] * relative_scale_ren_TSV[i_s][i_r], 2));
	}
	else {
	  value_alpha_S_TSV[0][v_sr].push_back(0.);
	  value_relative_factor_alpha_S[0][v_sr].push_back(0.);
	  value_scale_ren[0][v_sr].push_back(0.);
	  value_scale2_ren[0][v_sr].push_back(0.);
	}
      }
    }
  }

  for (int i_a = 1; i_a < n_ps; i_a++){
    value_alpha_S_TSV[i_a] = value_alpha_S_TSV[0];
    value_relative_factor_alpha_S[i_a] = value_relative_factor_alpha_S[0];
    value_scale_ren[i_a] = value_scale_ren[0];
    value_scale2_ren[i_a] = value_scale2_ren[0];
  }
  
  logger << LOG_DEBUG << "value_relative_scale_ren before 'sort':" << endl;
  for (int v_sr = 0; v_sr < max_dyn_ren + 1; v_sr++){
    for (int v_xr = 0; v_xr < value_relative_scale_ren[v_sr].size(); v_xr++){
      logger << LOG_DEBUG << "value_relative_scale_ren[" << setw(2) << v_sr << "][" << setw(2) << v_xr << "] = " << value_relative_scale_ren[v_sr][v_xr] << endl;
    }
  }
  logger.newLine(LOG_DEBUG);

  for (int v_sr = 0; v_sr < max_dyn_ren + 1; v_sr++){
    sort(value_relative_scale_ren[v_sr].begin(), value_relative_scale_ren[v_sr].end());
    n_scale_dyn_ren[v_sr] = value_relative_scale_ren[v_sr].size();
  }

  logger << LOG_DEBUG << "value_relative_scale_ren after 'sort':" << endl;
  for (int v_sr = 0; v_sr < max_dyn_ren + 1; v_sr++){
    for (int v_xr = 0; v_xr < value_relative_scale_ren[v_sr].size(); v_xr++){
      logger << LOG_DEBUG << "value_relative_scale_ren[" << setw(2) << v_sr << "][" << setw(2) << v_xr << "] = " << value_relative_scale_ren[v_sr][v_xr] << endl;
    }
  }
  logger.newLine(LOG_DEBUG);


  
  no_value_ren_TSV.resize(n_set_TSV);
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    int v_sr = dynamic_scale_ren_TSV[i_s];
    no_value_ren_TSV[i_s].resize(value_relative_scale_ren[v_sr].size());
    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
      for (int v_xr = 0; v_xr < value_relative_scale_ren[v_sr].size(); v_xr++){
	if (abs(value_relative_scale_ren[v_sr][v_xr] - relative_central_scale_ren_TSV[i_s] * relative_scale_ren_TSV[i_s][i_r]) / value_relative_scale_ren[v_sr][v_xr] < 1.e-15){
	  no_value_ren_TSV[i_s][i_r] = v_xr;
	}
      }
    }
  }

  logger << LOG_DEBUG << "value_relative_scale_ren.size() = " << value_relative_scale_ren.size() << endl;
  for (int v_sr = 0; v_sr < max_dyn_ren + 1; v_sr++){
    logger << LOG_DEBUG << "n_scale_dyn_ren[" << v_sr << "] = " << n_scale_dyn_ren[v_sr] << endl;
    logger << LOG_DEBUG << "value_relative_scale_ren[" << v_sr << "].size() = " << value_relative_scale_ren[v_sr].size() << endl;
    for (int v_xr = 0; v_xr < value_relative_scale_ren[v_sr].size(); v_xr++){
      logger << LOG_DEBUG << "value_relative_scale_ren[" << setw(2) << v_sr << "][" << setw(2) << v_xr << "] = " << value_relative_scale_ren[v_sr][v_xr] << endl;
    }
  logger.newLine(LOG_DEBUG);
  }
  logger.newLine(LOG_DEBUG);

  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    logger << LOG_DEBUG << name_set_TSV[i_s] << "   dynamic_scale_ren_TSV[" << setw(2) << i_s << "] = " << dynamic_scale_ren_TSV[i_s] << endl;
    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
      logger << LOG_DEBUG << "no_value_ren_TSV[" << setw(2) << i_s << "][" << setw(2) << i_r << "] = " << no_value_ren_TSV[i_s][i_r] << endl;
    }
  logger.newLine(LOG_DEBUG);
  }

  // initialization of all renormalization-scale connected variables done
  logger << LOG_DEBUG << "initialization of all renormalization-scale connected variables done" << endl;

  // initialization of all factorization-scale connected variables

  // max_dyn_fact
  // n_scale_dyn_fact
  // value_relative_scale_fact
  // value_relative_scale2_fact
  // value_relative_logscale2_fact
  // value_scale_fact
  // value_scale2_fact
  // value_pdf_factor
  // value_no_central_scale_fact   ->   doesn't seem to be used !!!
  // value_central_scale_fact
  // value_central_logscale2_fact
  // no_value_fact_TSV

  // initialization of variables relevant for debugging KP terms against Sherpa

  // value_pdf_factor_1
  // value_pdf_factor_2
  // value_pdf_factor_combination_1
  // value_pdf_factor_combination_2

  max_dyn_fact = 0;
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    if (dynamic_scale_fact_TSV[i_s] > max_dyn_fact){
      max_dyn_fact = dynamic_scale_fact_TSV[i_s];
    }
  }
  n_scale_dyn_fact.resize(max_dyn_fact + 1, 0);
  //!!!  n_scale_dyn_fact.resize(max_dyn_fact + 1);
  
  value_relative_scale_fact.resize(max_dyn_fact + 1);
  
  value_relative_scale2_fact.resize(max_dyn_fact + 1);

  value_relative_logscale2_fact.resize(max_dyn_fact + 1);
  
  value_scale_fact.resize(n_ps);
  value_scale_fact[0].resize(max_dyn_fact + 1);
  
  value_scale2_fact.resize(n_ps);
  value_scale2_fact[0].resize(max_dyn_fact + 1);
  
  value_pdf_factor.resize(n_pc, vector<vector<vector<vector<double> > > > (n_pz));

  // initialization of variables relevant for debugging KP terms against Sherpa
  /*
  value_pdf_factor_1.resize(n_pc, vector<vector<vector<vector<double> > > > (n_pz));
  value_pdf_factor_2.resize(n_pc, vector<vector<vector<vector<double> > > > (n_pz));
  */
  if (csi->class_contribution_CS_collinear){
    //    if (csi->type_contribution == "CA" || csi->type_contribution == "RCA"){
    value_pdf_factor_combination_1.resize(n_pc, vector<vector<vector<vector<vector<double> > > > > (n_pz));
    value_pdf_factor_combination_2.resize(n_pc, vector<vector<vector<vector<vector<double> > > > > (n_pz));
  }
  // initialization of variables relevant for debugging KP terms against Sherpa done
  

  if (n_pc > 0){
    value_pdf_factor[0][0].resize(max_dyn_fact + 1);
    
    // initialization of variables relevant for debugging KP terms against Sherpa
    /*
    value_pdf_factor_1[0][0].resize(max_dyn_fact + 1);
    value_pdf_factor_2[0][0].resize(max_dyn_fact + 1);
    */
    if (csi->class_contribution_CS_collinear){
    //    if (csi->type_contribution == "CA" || csi->type_contribution == "RCA"){
      value_pdf_factor_combination_1[0][0].resize(max_dyn_fact + 1);
      value_pdf_factor_combination_2[0][0].resize(max_dyn_fact + 1);
    }
    // initialization of variables relevant for debugging KP terms against Sherpa done
  }

  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    int v_sf = dynamic_scale_fact_TSV[i_s];
    for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
      int flag = 0;
      for (int v_xf = 0; v_xf < value_relative_scale_fact[v_sf].size(); v_xf++){
	if (abs(value_relative_scale_fact[v_sf][v_xf] - relative_central_scale_fact_TSV[i_s] * relative_scale_fact_TSV[i_s][i_f]) / value_relative_scale_fact[v_sf][v_xf] < 1.e-15){
	  flag = 1; 
	  break;
	}
      }
      if (flag == 0){
	value_relative_scale_fact[v_sf].push_back(relative_central_scale_fact_TSV[i_s] * relative_scale_fact_TSV[i_s][i_f]);
	value_relative_scale2_fact[v_sf].push_back(pow(relative_central_scale_fact_TSV[i_s] * relative_scale_fact_TSV[i_s][i_f], 2));
	value_relative_logscale2_fact[v_sf].push_back(2 * log(relative_central_scale_fact_TSV[i_s] * relative_scale_fact_TSV[i_s][i_f]));
	if (v_sf == 0){
	  value_scale_fact[0][v_sf].push_back(central_scale_fact_TSV[i_s] * relative_scale_fact_TSV[i_s][i_f]);
	  value_scale2_fact[0][v_sf].push_back(pow(central_scale_fact_TSV[i_s] * relative_scale_fact_TSV[i_s][i_f], 2));
	}
	else {
	  value_scale_fact[0][v_sf].push_back(0.);
	  value_scale2_fact[0][v_sf].push_back(0.);
	}
	if (n_pc > 0){
	  value_pdf_factor[0][0][v_sf].push_back(vector<double> (3, 0.));
	  
	  // initialization of variables relevant for debugging KP terms against Sherpa
	  if (csi->class_contribution_CS_collinear){
	    value_pdf_factor_combination_1[0][0][v_sf].resize(value_relative_scale_fact[v_sf].size(), vector<vector<double> > (3));
	    value_pdf_factor_combination_2[0][0][v_sf].resize(value_relative_scale_fact[v_sf].size(), vector<vector<double> > (3));
	  }
	  // initialization of variables relevant for debugging KP terms against Sherpa done
	}
      }
    }
  }

  for (int i_c = 0; i_c < n_pc; i_c++){
    for (int i_z = 0; i_z < n_pz; i_z++){
      value_pdf_factor[i_c][i_z] = value_pdf_factor[0][0];

      // initialization of variables relevant for debugging KP terms against Sherpa
      if (csi->class_contribution_CS_collinear){
	value_pdf_factor_combination_1[i_c][i_z] = value_pdf_factor_combination_1[0][0];
	value_pdf_factor_combination_2[i_c][i_z] = value_pdf_factor_combination_2[0][0];
      }
      // initialization of variables relevant for debugging KP terms against Sherpa done

    }
  }
  logger << LOG_DEBUG << "after Sherpa3" << endl;


  logger << LOG_DEBUG << "value_relative_scale_fact before 'sort':" << endl;
  for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
    for (int v_xf = 0; v_xf < value_relative_scale_fact[v_sf].size(); v_xf++){
      logger << LOG_DEBUG << "value_relative_scale_fact[" << setw(2) << v_sf << "][" << setw(2) << v_xf << "] = " << value_relative_scale_fact[v_sf][v_xf] << endl;
    }
  }
  logger.newLine(LOG_DEBUG);
 
  for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
    sort(value_relative_scale_fact[v_sf].begin(), value_relative_scale_fact[v_sf].end());
    // montone functions: same sorting automatically:
    sort(value_relative_scale2_fact[v_sf].begin(), value_relative_scale2_fact[v_sf].end());
    sort(value_relative_logscale2_fact[v_sf].begin(), value_relative_logscale2_fact[v_sf].end());
    
    n_scale_dyn_fact[v_sf] = value_relative_scale_fact[v_sf].size();
  }
 
  logger << LOG_DEBUG << "value_relative_scale_fact after 'sort':" << endl;
  for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
    for (int v_xf = 0; v_xf < value_relative_scale_fact[v_sf].size(); v_xf++){
      logger << LOG_DEBUG << "value_relative_scale_fact[" << setw(2) << v_sf << "][" << setw(2) << v_xf << "] = " << value_relative_scale_fact[v_sf][v_xf] << endl;
    }
  }
  logger.newLine(LOG_DEBUG);


  
  value_no_central_scale_fact.resize(max_dyn_fact + 1, -1);
  value_central_scale_fact.resize(max_dyn_fact + 1);
  value_central_logscale2_fact.resize(max_dyn_fact + 1);
  
  no_value_fact_TSV.resize(n_set_TSV);
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    int v_sf = dynamic_scale_fact_TSV[i_s];
    no_value_fact_TSV[i_s].resize(n_scale_dyn_fact[v_sf]);
    for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
      for (int v_xf = 0; v_xf < value_relative_scale_fact[v_sf].size(); v_xf++){
	if (abs(value_relative_scale_fact[v_sf][v_xf] - relative_central_scale_fact_TSV[i_s] * relative_scale_fact_TSV[i_s][i_f]) / value_relative_scale_fact[v_sf][v_xf] < 1.e-15){
	  no_value_fact_TSV[i_s][i_f] = v_xf;

	  // ???
	  if (abs(value_relative_scale_fact[v_sf][v_xf] - 1.) < 1.e-15){
	    value_no_central_scale_fact[v_sf] = v_xf;
	    //	    logger << LOG_DEBUG << "   relative_scale_fact_TSV[" << i_s << "][" << i_f << "] = " << relative_scale_fact_TSV[i_s][i_f] << "   value_no_central_scale_fact[" << v_sf << "] = " << value_no_central_scale_fact[v_sf] << endl;
	  }
	  // ???
	}
      }
    }
  }

  logger << LOG_DEBUG << "value_relative_scale_fact.size() = " << value_relative_scale_fact.size() << endl;
  for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
    logger << LOG_DEBUG << "n_scale_dyn_fact[" << v_sf << "] = " << n_scale_dyn_fact[v_sf] << endl;
    logger << LOG_DEBUG << "value_relative_scale_fact[" << v_sf << "].size() = " << value_relative_scale_fact[v_sf].size() << endl;
    for (int v_xf = 0; v_xf < value_relative_scale_fact[v_sf].size(); v_xf++){
      logger << LOG_DEBUG << "value_relative_scale_fact[" << setw(2) << v_sf << "][" << setw(2) << v_xf << "] = " << value_relative_scale_fact[v_sf][v_xf] << endl;
    }
  logger.newLine(LOG_DEBUG);
  }
  logger.newLine(LOG_DEBUG);

  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    logger << LOG_DEBUG << name_set_TSV[i_s] << "   dynamic_scale_fact_TSV[" << setw(2) << i_s << "] = " << dynamic_scale_fact_TSV[i_s] << endl;
    for (int i_r = 0; i_r < n_scale_fact_TSV[i_s]; i_r++){
      logger << LOG_DEBUG << "no_value_fact_TSV[" << setw(2) << i_s << "][" << setw(2) << i_r << "] = " << no_value_fact_TSV[i_s][i_r] << endl;
    }
  logger.newLine(LOG_DEBUG);
  }

  /*
  for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
    logger << LOG_DEBUG << "value_no_central_scale_fact[" << v_sf << "] = " << value_no_central_scale_fact[v_sf] << endl;
  }
  */
  
  for (int i_a = 1; i_a < n_ps; i_a++){
    value_scale_fact[i_a] = value_scale_fact[0];
    value_scale2_fact[i_a] = value_scale2_fact[0];
  }

  // initialization of all factorization-scale connected variables done

  
  // initialization of value_ME2term_...
  logger << LOG_DEBUG << "initialization of value_ME2term_..." << endl;
  
  //  value_ME2term.resize(n_pc, 0.); needed also without TSV !!!
  value_ME2term_ren.resize(n_pc, vector<vector<double> > (max_dyn_ren + 1));
  value_ME2term_fact.resize(n_pc, vector<vector<vector<double> > > (n_pz, vector<vector<double> > (max_dyn_fact + 1)));
  
  for (int i_c = 0; i_c < n_pc; i_c++){
    for (int v_sr = 0; v_sr < max_dyn_ren + 1; v_sr++){  //for (int i_m = 0; i_m < value_relative_scale_ren[v_sr].size(); i_m++){
      value_ME2term_ren[i_c][v_sr].resize(n_scale_dyn_ren[v_sr], 0.);
    }
    for (int i_z = 0; i_z < n_pz; i_z++){
      for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	value_ME2term_fact[i_c][i_z][v_sf].resize(n_scale_dyn_fact[v_sf], 0.);
      }
    }
  }

  // initialization of value_ME2term_... done
  logger << LOG_DEBUG << "initialization of value_ME2term_... done" << endl;


  /*
  // initialization of new 'combination' variables: could be used to reduce the number of needed dyn_ren - dyn_fact combinations
  logger << LOG_DEBUG << "initialization of new 'combination' variables: could be used to reduce the number of needed dyn_ren - dyn_fact combinations" << endl;
  
  vector<vector<int> > combination_TSV(n_set_TSV, vector<int> (2, 0));
  logger << LOG_DEBUG << "combination_TSV.size() = " << combination_TSV.size() << endl;
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    logger << LOG_DEBUG << "combination_TSV[" << i_s << "].size() = " << combination_TSV[i_s].size() << endl;
    combination_TSV[i_s][0] = dynamic_scale_ren_TSV[i_s];
    combination_TSV[i_s][1] = dynamic_scale_fact_TSV[i_s];
  }
  logger << LOG_DEBUG << "combination_TSV done" << endl;
  
  vector<vector<int> > value_combination(0);
  vector<int> relation_combination(n_set_TSV);
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    int flag = 0;
    for (int i_v = 0; i_v < value_combination.size(); i_v++){
      if (combination_TSV[i_s] == value_combination[i_v]){flag = i_v; break;}
    }
    logger << LOG_DEBUG << "i_s = " << i_s << endl;
    if (flag == 0){
      value_combination.push_back(combination_TSV[i_s]);
      relation_combination[i_s] = flag;
    }
  }
  
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    logger << LOG_DEBUG_VERBOSE << "combination_TSV[" << setw(2) << i_s << "][0] = " << combination_TSV[i_s][0] << "   [1] = " << combination_TSV[i_s][1] << "   relation_combination = " << relation_combination[i_s] << endl;
  }
  for (int i_v = 0; i_v < value_combination.size(); i_v++){
    logger << LOG_DEBUG_VERBOSE << "value_combination[" << setw(2) << i_v << "][0] = " << value_combination[i_v][0] << "   [1] = " << value_combination[i_v][1] << endl;
  }

  // initialization of new 'combination' variables: could be used to reduce the number of needed dyn_ren - dyn_fact combinations done
  logger << LOG_DEBUG << "initialization of new 'combination' variables: could be used to reduce the number of needed dyn_ren - dyn_fact combinations done" << endl;
  */
  









  // initialization of 'integrand, weight, moment and variations

  // ps_integrand_qTcut_TSV
  // integrand_qTcut_TSV
  // sum_weight_TSV
  // sum_weight2_TSV
  // fullsum_weight_TSV
  // fullsum_weight2_TSV
  // sum_weight_qTcut_TSV
  // sum_weight2_qTcut_TSV
  // fullsum_weight_qTcut_TSV
  // fullsum_weight2_qTcut_TSV

  // ps_moment_TSV
  // moment_TSV
  // sum_moment_TSV
  // sum_moment2_TSV
  // fullsum_moment_TSV
  // fullsum_moment2_TSV
  // sum_moment_qTcut_TSV
  // sum_moment2_qTcut_TSV
  // fullsum_moment_qTcut_TSV
  // fullsum_moment2_qTcut_TSV

  // Xsection_TSV
  // Xsection_delta_TSV


  ps_integrand_qTcut_TSV.resize(n_qTcut, vector<vector<vector<vector<vector<double> > > > > (n_extended_set_TSV));
  integrand_qTcut_TSV.resize(n_qTcut, vector<vector<vector<double> > > (n_extended_set_TSV));

  ps_integrand_TSV.resize(n_extended_set_TSV);
  integrand_TSV.resize(n_extended_set_TSV);
  sum_weight_TSV.resize(n_extended_set_TSV);
  sum_weight2_TSV.resize(n_extended_set_TSV);
  fullsum_weight_TSV.resize(n_extended_set_TSV);
  fullsum_weight2_TSV.resize(n_extended_set_TSV);
  sum_weight_qTcut_TSV.resize(n_extended_set_TSV);
  sum_weight2_qTcut_TSV.resize(n_extended_set_TSV);
  fullsum_weight_qTcut_TSV.resize(n_extended_set_TSV);
  fullsum_weight2_qTcut_TSV.resize(n_extended_set_TSV);

  // only if moments are switched on !!!
  ps_moment_TSV.resize(n_extended_set_TSV);
  moment_TSV.resize(n_extended_set_TSV);
  sum_moment_TSV.resize(n_extended_set_TSV);
  sum_moment2_TSV.resize(n_extended_set_TSV);
  fullsum_moment_TSV.resize(n_extended_set_TSV);
  fullsum_moment2_TSV.resize(n_extended_set_TSV);
  sum_moment_qTcut_TSV.resize(n_extended_set_TSV);
  sum_moment2_qTcut_TSV.resize(n_extended_set_TSV);
  fullsum_moment_qTcut_TSV.resize(n_extended_set_TSV);
  fullsum_moment2_qTcut_TSV.resize(n_extended_set_TSV);
  
  Xsection_TSV.resize(n_extended_set_TSV);
  Xsection_delta_TSV.resize(n_extended_set_TSV);

  logger << LOG_DEBUG << "fullsum_moment2_qTcut_TSV done" << endl;
  

  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
      ps_integrand_qTcut_TSV[i_q][i_s].resize(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (max_n_integrand_TSV[i_s], vector<double> (n_pc, 0.)))); // check if really ps !!!
      integrand_qTcut_TSV[i_q][i_s].resize(n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.));
    }

    ps_integrand_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (max_n_integrand_TSV[i_s], vector<double> (n_pc, 0.)))); // check if really ps !!!
    integrand_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.));
    sum_weight_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.));
    sum_weight2_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.));
    fullsum_weight_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.));
    fullsum_weight2_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.));
    if (active_qTcut){
      sum_weight_qTcut_TSV[i_s].resize(n_qTcut, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.)));
      sum_weight2_qTcut_TSV[i_s].resize(n_qTcut, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.)));
      fullsum_weight_qTcut_TSV[i_s].resize(n_qTcut, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.)));
      fullsum_weight2_qTcut_TSV[i_s].resize(n_qTcut, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.)));
    }
    
    if (switch_moment_TSV[i_s] > 0){
      ps_moment_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<vector<vector<vector<double> > > > (n_moments, vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (max_n_integrand_TSV[i_s], vector<double> (n_ps, 0.)))));
      sum_moment_TSV[i_s].resize(n_moments, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.)));
      sum_moment2_TSV[i_s].resize(n_moments, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.)));
      fullsum_moment_TSV[i_s].resize(n_moments, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.)));
      fullsum_moment2_TSV[i_s].resize(n_moments, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.)));
      
      if (active_qTcut){
	sum_moment_qTcut_TSV[i_s].resize(n_moments, vector<vector<vector<double> > > (n_qTcut, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.))));
	sum_moment2_qTcut_TSV[i_s].resize(n_moments, vector<vector<vector<double> > > (n_qTcut, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.))));
	fullsum_moment_qTcut_TSV[i_s].resize(n_moments, vector<vector<vector<double> > > (n_qTcut, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.))));
	fullsum_moment2_qTcut_TSV[i_s].resize(n_moments, vector<vector<vector<double> > > (n_qTcut, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.))));
      }
    }

    Xsection_TSV[i_s].resize(output_n_qTcut, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.)));
    Xsection_delta_TSV[i_s].resize(output_n_qTcut, vector<vector<double> > (n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.)));
  }
  
  // initialization of 'integrand, weight, moment and variations done




  // initialization of pointers for renormalization/factorization dependent variables
  logger << LOG_DEBUG << "initialization of pointers for renormalization/factorization dependent variables" << endl;

  // pointer_scale_ren.resize
  // pointer_scale2_ren
  // pointer_relative_factor_alpha_S
  // pointer_scale_fact
  // pointer_scale2_fact
  // pointer_pdf_factor

  pointer_scale_ren.resize(n_ps, vector<vector<double*> > (n_set_TSV));
  pointer_scale2_ren.resize(n_ps, vector<vector<double*> > (n_set_TSV));
  pointer_alpha_S_TSV.resize(n_ps, vector<vector<double*> > (n_set_TSV));
  pointer_relative_factor_alpha_S.resize(n_ps, vector<vector<double*> > (n_set_TSV));
  pointer_scale_fact.resize(n_ps, vector<vector<double*> > (n_set_TSV));
  pointer_scale2_fact.resize(n_ps, vector<vector<double*> > (n_set_TSV));
  
  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i_s = 0; i_s < n_set_TSV; i_s++){
      pointer_alpha_S_TSV[i_a][i_s].resize(n_scale_ren_TSV[i_s]);
      pointer_relative_factor_alpha_S[i_a][i_s].resize(n_scale_ren_TSV[i_s]);
      pointer_scale_ren[i_a][i_s].resize(n_scale_ren_TSV[i_s]);
      pointer_scale2_ren[i_a][i_s].resize(n_scale_ren_TSV[i_s]);
      pointer_scale_fact[i_a][i_s].resize(n_scale_fact_TSV[i_s]);
      pointer_scale2_fact[i_a][i_s].resize(n_scale_fact_TSV[i_s]);
      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	pointer_scale_ren[i_a][i_s][i_r] = &value_scale_ren[i_a][dynamic_scale_ren_TSV[i_s]][no_value_ren_TSV[i_s][i_r]];
	pointer_scale2_ren[i_a][i_s][i_r] = &value_scale2_ren[i_a][dynamic_scale_ren_TSV[i_s]][no_value_ren_TSV[i_s][i_r]];
	pointer_alpha_S_TSV[i_a][i_s][i_r] = &value_alpha_S_TSV[i_a][dynamic_scale_ren_TSV[i_s]][no_value_ren_TSV[i_s][i_r]];
	pointer_relative_factor_alpha_S[i_a][i_s][i_r] = &value_relative_factor_alpha_S[i_a][dynamic_scale_ren_TSV[i_s]][no_value_ren_TSV[i_s][i_r]];
	pointer_scale_fact[i_a][i_s][i_r] = &value_scale_fact[i_a][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_r]];
	pointer_scale2_fact[i_a][i_s][i_r] = &value_scale2_fact[i_a][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_r]];
	logger << LOG_DEBUG << setw(40) << right << "pointer_scale_ren[" << setw(2) << i_a << "][" << setw(2) << i_s << "][" << setw(2) << i_r << "] = " << *pointer_scale_ren[i_a][i_s][i_r] << endl;
      }
    }
  }

  pointer_pdf_factor.resize(n_pc, vector<vector<vector<vector<double*> > > > (n_pz, vector<vector<vector<double*> > > (n_set_TSV)));
  for (int i_c = 0; i_c < n_pc; i_c++){
    for (int i_z = 0; i_z < n_pz; i_z++){
      for (int i_s = 0; i_s < n_set_TSV; i_s++){
	pointer_pdf_factor[i_c][i_z][i_s].resize(n_scale_fact_TSV[i_s], vector<double*> (3));
	for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	  for (int i_h = 0; i_h < 3; i_h++){
	    pointer_pdf_factor[i_c][i_z][i_s][i_f][i_h] = &value_pdf_factor[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h];
	  }
	}
      }
    }
  }


  // initialization of variables relevant for debugging KP terms against Sherpa

  logger << LOG_DEBUG << "before pointer_pdf_factor_1/2 for debugging against Sherpa" << endl;
  /*
  pointer_pdf_factor_1.resize(n_pc, vector<vector<vector<vector<double*> > > > (n_pz, vector<vector<vector<double*> > > (n_set_TSV)));
  pointer_pdf_factor_2.resize(n_pc, vector<vector<vector<vector<double*> > > > (n_pz, vector<vector<vector<double*> > > (n_set_TSV)));
  for (int i_c = 0; i_c < n_pc; i_c++){
    for (int i_z = 0; i_z < n_pz; i_z++){
      for (int i_s = 0; i_s < n_set_TSV; i_s++){
	pointer_pdf_factor_1[i_c][i_z][i_s].resize(n_scale_fact_TSV[i_s], vector<double*> (3));
	pointer_pdf_factor_2[i_c][i_z][i_s].resize(n_scale_fact_TSV[i_s], vector<double*> (3));
	for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	  for (int i_h = 0; i_h < 3; i_h++){
	    pointer_pdf_factor_1[i_c][i_z][i_s][i_f][i_h] = &value_pdf_factor_1[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h];
	    pointer_pdf_factor_2[i_c][i_z][i_s][i_f][i_h] = &value_pdf_factor_2[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h];
	  }
	}
      }
    }
  }
  */
  if (csi->class_contribution_CS_collinear){
    //    if (csi->type_contribution == "CA" || csi->type_contribution == "RCA"){
    logger << LOG_DEBUG_VERBOSE << "Initialization phase of value_pdf_factor_combination_1:" << endl;
    for (int i_c = 0; i_c < n_pc; i_c++){
      for (int i_z = 0; i_z < n_pz; i_z++){
	logger << LOG_DEBUG_VERBOSE << "max_dyn_fact = " << max_dyn_fact << endl;
	for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	  value_pdf_factor_combination_1[i_c][i_z][v_sf].resize(n_scale_dyn_fact[v_sf], vector<vector<double> > (3));
	  value_pdf_factor_combination_2[i_c][i_z][v_sf].resize(n_scale_dyn_fact[v_sf], vector<vector<double> > (3));
	  logger << LOG_DEBUG_VERBOSE << "n_scale_dyn_fact[" << v_sf << "] = " << n_scale_dyn_fact[v_sf] << endl;
	  for (int i_m = 0; i_m < n_scale_dyn_fact[v_sf]; i_m++){
	    for (int i_h = 0; i_h < 3; i_h++){
	      value_pdf_factor_combination_1[i_c][i_z][v_sf][i_m][i_h].resize(CA_combination_pdf[i_c].size(), 0.);
	      value_pdf_factor_combination_2[i_c][i_z][v_sf][i_m][i_h].resize(CA_combination_pdf[i_c].size(), 0.);
	    }
	  }
	}
      }
    }
    
    pointer_pdf_factor_combination_1.resize(n_pc, vector<vector<vector<vector<vector<double*> > > > > (n_pz, vector<vector<vector<vector<double*> > > > (n_set_TSV)));
    pointer_pdf_factor_combination_2.resize(n_pc, vector<vector<vector<vector<vector<double*> > > > > (n_pz, vector<vector<vector<vector<double*> > > > (n_set_TSV)));
    for (int i_c = 0; i_c < n_pc; i_c++){
      for (int i_z = 0; i_z < n_pz; i_z++){
	for (int i_s = 0; i_s < n_set_TSV; i_s++){
	  logger << LOG_DEBUG_VERBOSE << "n_scale_fact_TSV[" << i_s << "] = " << n_scale_fact_TSV[i_s] << endl;
	  pointer_pdf_factor_combination_1[i_c][i_z][i_s].resize(n_scale_fact_TSV[i_s], vector<vector<double*> > (3));
	  pointer_pdf_factor_combination_2[i_c][i_z][i_s].resize(n_scale_fact_TSV[i_s], vector<vector<double*> > (3));
	  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	    for (int i_h = 0; i_h < 3; i_h++){
	      pointer_pdf_factor_combination_1[i_c][i_z][i_s][i_f][i_h].resize(CA_combination_pdf[i_c].size());
	      pointer_pdf_factor_combination_2[i_c][i_z][i_s][i_f][i_h].resize(CA_combination_pdf[i_c].size());
	      for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
		pointer_pdf_factor_combination_1[i_c][i_z][i_s][i_f][i_h][i_i] = &value_pdf_factor_combination_1[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h][i_i];
		pointer_pdf_factor_combination_2[i_c][i_z][i_s][i_f][i_h][i_i] = &value_pdf_factor_combination_2[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h][i_i];
	      }
	    }
	  }
	}
      }
    }

    logger << LOG_DEBUG_VERBOSE << "Initialization phase of value_pdf_factor_combination_1: finished" << endl;

    for (int i_c = 0; i_c < n_pc; i_c++){
      for (int i_z = 0; i_z < n_pz; i_z++){
	for (int i_s = 0; i_s < n_set_TSV; i_s++){
	  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	    for (int i_h = 0; i_h < 3; i_h++){
	      for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
		//		pointer_pdf_factor_combination_1[i_c][i_z][i_s][i_f][i_h][i_i]
		logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_1[" << i_c << "][" << i_z << "][" << dynamic_scale_fact_TSV[i_s] << "][" << no_value_fact_TSV[i_s][i_f] << "][" << i_h << "][" << i_i << "] = " << value_pdf_factor_combination_1[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h][i_i] << endl;
		//		pointer_pdf_factor_combination_2[i_c][i_z][i_s][i_f][i_h][i_i] = &value_pdf_factor_combination_2[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h][i_i];
		logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_2[" << i_c << "][" << i_z << "][" << dynamic_scale_fact_TSV[i_s] << "][" << no_value_fact_TSV[i_s][i_f] << "][" << i_h << "][" << i_i << "] = " << value_pdf_factor_combination_2[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h][i_i] << endl;
	      }
	    }
	  }
	}
      }
    }

  /*
	if (i_v == 0 && relation_pc_ps[i_c] > 0){
	  value_pdf_factor_combination_1[i_c][i_z][i_v] = value_pdf_factor_combination_1[0][i_z][i_v];
	  value_pdf_factor_combination_2[i_c][i_z][i_v] = value_pdf_factor_combination_2[0][i_z][i_v];
	  // check
	  for (int i_m = 0; i_m < n_scale_dyn_fact[i_v]; i_m++){
	    for (int i_h = 0; i_h < 3; i_h++){
	      value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].resize(CA_combination_pdf[i_c].size(), 0.);
	      value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][i_h].resize(CA_combination_pdf[i_c].size(), 0.);
	      logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].size() = " << value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].size() << endl;
	    }
	  }
	  // !!!
	}
	else {
	  value_pdf_factor_combination_1[i_c][i_z][i_v].resize(n_scale_dyn_fact[i_v], vector<vector<double> > (3));
	  value_pdf_factor_combination_2[i_c][i_z][i_v].resize(n_scale_dyn_fact[i_v], vector<vector<double> > (3));
	  for (int i_m = 0; i_m < n_scale_dyn_fact[i_v]; i_m++){
	    for (int i_h = 0; i_h < 3; i_h++){
	      value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].resize(CA_combination_pdf[i_c].size(), 0.);
	      value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][i_h].resize(CA_combination_pdf[i_c].size(), 0.);
	      logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].size() = " << value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].size() << endl;
	    }
	  }
	}
      }
    }
  }
  */




  }

  // initialization of variables relevant for debugging KP terms against Sherpa done

  // initialization of pointers for renormalization/factorization dependent variables done
  logger << LOG_DEBUG << "initialization of pointers for renormalization/factorization dependent variables done" << endl;



  // initialization of pointer_ME2term
  // check if pointer_ME2term is actually used !!!
  logger << LOG_DEBUG << "initialization of pointer_ME2term" << endl;
 
  pointer_ME2term.resize(n_pc, vector<vector<vector<vector<double*> > > > (n_pz, vector<vector<vector<double*> > > (n_set_TSV)));
  for (int i_c = 0; i_c < n_pc; i_c++){
    for (int i_z = 0; i_z < n_pz; i_z++){
      for (int i_s = 0; i_s < n_set_TSV; i_s++){
	pointer_ME2term[i_c][i_z][i_s].resize(n_scale_ren_TSV[i_s]);
	for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	  pointer_ME2term[i_c][i_z][i_s][i_r].resize(n_scale_fact_TSV[i_s]);
	  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	    if (csi->type_contribution == "born" ||
		csi->type_contribution == "loop" ||
		csi->type_contribution == "L2I" ||
		csi->type_contribution == "RT" ||
		csi->type_contribution == "L2RT" ||
		csi->type_contribution == "RJ" ||
		csi->type_contribution == "L2RJ" ||
		csi->type_contribution == "RA" ||
		csi->type_contribution == "RRA" ||
		csi->type_contribution == "RRJ" ||
		csi->type_contribution == "L2RA" ||
		csi->type_contribution == "CT" ||
		csi->type_contribution == "L2CT" ||
		csi->type_contribution == "CJ"   ||
		csi->type_contribution == "L2CJ"){
	      pointer_ME2term[i_c][i_z][i_s][i_r][i_f] = &value_ME2term[i_c];
	    }
	    //	    else if (csi->class_contribution_CS_collinear){
	    else if (csi->type_contribution == "CA" ||
		     csi->type_contribution == "RCA" ||
		     csi->type_contribution == "RCJ" ||
		     csi->type_contribution == "L2CA" ||
		     csi->type_contribution == "CT2" ||
		     csi->type_contribution == "CJ2"){
	      pointer_ME2term[i_c][i_z][i_s][i_r][i_f] = &value_ME2term_fact[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]];
	    }
	    else if (csi->type_contribution == "VA" ||
		     csi->type_contribution == "RVA" ||
		     csi->type_contribution == "RVJ" ||
		     csi->type_contribution == "L2VA" ||
		     csi->type_contribution == "VT" ||
		     csi->type_contribution == "VT2" ||
		     csi->type_contribution == "L2VT" ||
		     csi->type_contribution == "VJ"  ||
		     csi->type_contribution == "VJ2" ||
		     csi->type_contribution == "L2VJ"){ // ???
	      pointer_ME2term[i_c][i_z][i_s][i_r][i_f] = &value_ME2term_ren[i_c][dynamic_scale_ren_TSV[i_s]][no_value_ren_TSV[i_s][i_r]];
	      
	      logger << LOG_DEBUG_VERBOSE << "pointer_ME2term[" << i_c << "][" << i_z << "][" << i_s << "][" << i_r << "][" << i_f << "] = &value_ME2term_ren[" << i_c << "][dynamic_scale_ren_TSV[" << i_s << "] = " << dynamic_scale_ren_TSV[i_s] << "][no_value_ren_TSV[" << i_s << "][" << i_r << "] = " << no_value_ren_TSV[i_s][i_r] << "];" << endl;
	    }
	    else {logger << LOG_DEBUG << "No valid subcontribution!" << endl; exit(1);}
	  }
	}
      }
    }
  }

  // initialization of pointer_ME2term done
  logger << LOG_DEBUG << "initialization of pointer_ME2term done" << endl;
  



  /*
  // initialization of min-max-qTcut values for TSV variation (not used so far !!!)
  logger << LOG_DEBUG << "initialization of min-max-qTcut values for TSV variation" << endl;
  
  no_min_qTcut_TSV.resize(n_set_TSV);
  no_max_qTcut_TSV.resize(n_set_TSV);
  no_min_qTcut_distribution_TSV.resize(n_set_TSV);
  no_max_qTcut_distribution_TSV.resize(n_set_TSV);
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    if (switch_qTcut == 0){
      no_min_qTcut_TSV[i_s] = 0;
      no_max_qTcut_TSV[i_s] = 1;
      no_min_qTcut_distribution_TSV[i_s] = 0;
      no_max_qTcut_distribution_TSV[i_s] = 1;
    }
    else {
      
      for (int i_q = 0; i_q < n_qTcut; i_q++){
	if (min_qTcut_TSV[i_s] > value_qTcut[i_q]){
	  no_min_qTcut_TSV[i_s] = i_q - 1;
	  break;
	}
      }
      for (int i_q = n_qTcut - 1; i_q >= 0; i_q--){
	if (max_qTcut_TSV[i_s] < value_qTcut[i_q]){
	  no_max_qTcut_TSV[i_s] = i_q + 1;
	  break;
	}
      }
      for (int i_q = 0; i_q < n_qTcut; i_q++){
	if (min_qTcut_distribution_TSV[i_s] > value_qTcut[i_q]){
	  no_min_qTcut_distribution_TSV[i_s] = i_q - 1;
	  break;
	}
      }
      for (int i_q = n_qTcut - 1; i_q >= 0; i_q--){
	if (max_qTcut_distribution_TSV[i_s] < value_qTcut[i_q]){
	  no_max_qTcut_distribution_TSV[i_s] = i_q + 1;
	  break;
	}
      }
    }
    logger << LOG_DEBUG << "no_min_qTcut_TSV[" << i_s << "] = " << no_min_qTcut_TSV[i_s] << endl;
    logger << LOG_DEBUG << "no_max_qTcut_TSV[" << i_s << "] = " << no_max_qTcut_TSV[i_s] << endl;
    logger << LOG_DEBUG << "no_min_qTcut_distribution_TSV[" << i_s << "] = " << no_min_qTcut_distribution_TSV[i_s] << endl;
    logger << LOG_DEBUG << "no_max_qTcut_distribution_TSV[" << i_s << "] = " << no_max_qTcut_distribution_TSV[i_s] << endl;
  }
  
  logger << LOG_DEBUG << "after min-max-qTcut determination" << endl;
  
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    logger << LOG_DEBUG << setw(3) << no_min_qTcut_TSV[i_s] << " <= no_qTcut_TSV <= " << setw(3) << no_max_qTcut_TSV[i_s] << "   " << setw(3) << no_min_qTcut_distribution_TSV[i_s] << " <= no_qTcut_distribution_TSV <= " << setw(3) << no_max_qTcut_distribution_TSV[i_s] << endl;
  }

  // initialization of min-max-qTcut values for TSV variation (not used so far !!!) done
  logger << LOG_DEBUG << "initialization of min-max-qTcut values for TSV variation done" << endl;
  */


  logger << LOG_DEBUG << "finished" << endl;
}





// Relevant for CV variation and reference result:




void observable_set::initialization_basic_CV(inputparameter_set & isi){
  Logger logger("observable_set::initialization_basic_CV (isi)");
  logger << LOG_DEBUG << "called" << endl;

  /*
  double empty_double = -1.;
  double empty_int = -1;
  */

  scale_ren = isi.scale_ren;
  scale_fact = isi.scale_fact;
  dynamic_scale = isi.dynamic_scale;
  prefactor_reference = isi.prefactor_reference;


  switch_CV = isi.switch_CV;
  variation_mu_ren_CV = isi.variation_mu_ren_CV;
  variation_mu_fact_CV = isi.variation_mu_fact_CV;
  prefactor_CV = isi.prefactor_CV;

  dynamic_scale_CV = isi.dynamic_scale_CV;
  n_scales_CV = isi.n_scales_CV;
  variation_factor_CV = isi.variation_factor_CV;

  scale_ren *= prefactor_reference;
  scale_fact *= prefactor_reference;


  // new
  if (switch_CV != 0){
    if (switch_CV == 1){variation_mu_ren_CV = 1; variation_mu_fact_CV = 1;}
    // equal variation of mu_ren and mu_fact
    else if (switch_CV == 2){variation_mu_ren_CV = 1; variation_mu_fact_CV = 0;}
    //  variation of mu_ren
    else if (switch_CV == 3){variation_mu_ren_CV = 0; variation_mu_fact_CV = 1;}
    //  variation of mu_fact
    else if (switch_CV == 4){variation_mu_ren_CV = 1; variation_mu_fact_CV = -1;}
    // antipodal variation of mu_ren and mu_fact
    else if (switch_CV == 5){variation_mu_ren_CV = 1; variation_mu_fact_CV = 1;}
    // 7-point variation of mu_ren and mu_fact
    else if (switch_CV == 6){variation_mu_ren_CV = 1; variation_mu_fact_CV = 1;}
    // 9-point variation of mu_ren and mu_fact
    else {logger << LOG_ERROR << "No correct switch_CV selected!" << endl; exit(1);}
  }
  else {
    if (variation_mu_ren_CV == 0 && variation_mu_fact_CV == 0){}
    if (variation_mu_ren_CV == 1 && variation_mu_fact_CV == 1){switch_CV = 1;}
    // equal variation of mu_ren and mu_fact
    else if (variation_mu_ren_CV == 1 && variation_mu_fact_CV == 0){switch_CV = 2;}
    //  variation of mu_ren
    else if (variation_mu_ren_CV == 0 && variation_mu_fact_CV == 1){switch_CV = 3;}
    //  variation of mu_fact
    else if (variation_mu_ren_CV == 1 && variation_mu_fact_CV == -1){switch_CV = 4;}
    // antipodal variation of mu_ren and mu_fact
    else {switch_CV = 0; variation_mu_ren_CV = 0; variation_mu_fact_CV = 0;}
    //    else {logger << LOG_ERROR << "No correct switch_CV selected!" << endl; exit(1);}
  }






  no_central_scale_CV = (n_scales_CV - 1) / 2; // !!! could be made element of observable_set
  directory_name_scale_CV.resize(n_scales_CV);
  for (int i_s = 0; i_s < n_scales_CV; i_s++){
    int j_s = n_scales_CV - i_s - 1;
    stringstream temp_ss;
    temp_ss << "scale";
    if (variation_mu_ren_CV == 1){temp_ss << "." << i_s;}
    if (variation_mu_ren_CV == 0){temp_ss << "." << no_central_scale_CV;}
    if (variation_mu_ren_CV == -1){temp_ss << "." << j_s;}

    if (variation_mu_fact_CV == 1){temp_ss << "." << i_s;}
    if (variation_mu_fact_CV == 0){temp_ss << "." << no_central_scale_CV;}
    if (variation_mu_fact_CV == -1){temp_ss << "." << j_s;}
    directory_name_scale_CV[i_s] = temp_ss.str();
  }

  logger << LOG_DEBUG << "switch_CV = " << switch_CV << endl;
  logger << LOG_DEBUG << "variation_mu_ren_CV = " << variation_mu_ren_CV << endl;
  logger << LOG_DEBUG << "variation_mu_fact_CV = " << variation_mu_fact_CV << endl;

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_CV(){
  Logger logger("observable_set::initialization_CV (---)");
  logger << LOG_DEBUG << "started" << endl;

  ME2 = 0.;
  // why here ???

  if (switch_CV == 1){id_scales = 1;}
  else {id_scales = 0;}

  central_scale_CV = no_central_scale_CV;
  // should be identical !!! check if needed !!!

  //  dynamic_scale = _dynamic_scale;

  mu_ren.resize(3);
  mu_fact.resize(3);
  mu_ren[0] = 1.;
  mu_ren[1] = scale_ren;
  mu_ren[2] = pow(scale_ren, 2);
  mu_fact[0] = log(mu_fact[2]);
  mu_fact[1] = scale_fact;
  mu_fact[2] = pow(scale_fact, 2);

  logger << LOG_DEBUG << "mu_fact[1] = " << mu_fact[1] << endl;
  // should be calculated here:
  /*
  mu_fact_CV = _mu_fact_CV;
  mu_ren_CV = _mu_ren_CV;
  */
  //  alpha_S = _alpha_S;
  alpha_S = LHAPDF::alphasPDF(mu_ren[1]);
  logger << LOG_DEBUG << "alpha_S = " << alpha_S << endl;

  //  no_central_scale_CV = (n_scales_CV - 1) / 2;

  logger << LOG_DEBUG << "n_scales_CV = " << n_scales_CV << endl;

  if (switch_CV){
    scale_ren_CV.resize(n_scales_CV);
    scale_fact_CV.resize(n_scales_CV);
    
    mu_ren_CV.resize(3, vector<double> (n_scales_CV));
    mu_fact_CV.resize(3, vector<double> (n_scales_CV));
    
    rel_scale_factor_CV.resize(n_scales_CV);
    
    rel_scale_factor_ren_CV.resize(n_scales_CV);
    rel_scale_factor_fact_CV.resize(n_scales_CV);
    
    if (switch_CV < 5){
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	if (n_scales_CV > 1){rel_scale_factor_CV[i_s] = pow(10., log10(double(variation_factor_CV)) * double(2. * i_s / (n_scales_CV - 1) - 1));}
	else {rel_scale_factor_CV[i_s] = 1.;}
	// !!!
	rel_scale_factor_CV[i_s] = rel_scale_factor_CV[i_s] * prefactor_CV / prefactor_reference;
	// !!!
      }
      
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	int j_s = n_scales_CV - i_s - 1;
	//    stringstream temp_ss;
	//    temp_ss << "scale." << i_s;
	if (variation_mu_ren_CV == 1){
	  //      scale_ren_CV[i_s] = rel_scale_factor_CV[i_s] * scale_ren;
	  rel_scale_factor_ren_CV[i_s] = rel_scale_factor_CV[i_s];
	  //  temp_ss << "." << i_s;
	}
	if (variation_mu_ren_CV == 0){
	  //      scale_ren_CV[i_s] = scale_ren;
	  rel_scale_factor_ren_CV[i_s] = prefactor_CV / prefactor_reference;
	  //	  rel_scale_factor_ren_CV[i_s] = 1.;
	  //  temp_ss << "." << no_central_scale_CV;
	}
	if (variation_mu_ren_CV == -1){
	  //      scale_ren_CV[i_s] = rel_scale_factor_CV[j_s] * scale_ren;
	  rel_scale_factor_ren_CV[i_s] = rel_scale_factor_CV[j_s];
	  //  temp_ss << "." << j_s;
	}
	
	if (variation_mu_fact_CV == 1){
	  //      scale_fact_CV[i_s] = rel_scale_factor_CV[i_s] * scale_fact;
	  rel_scale_factor_fact_CV[i_s] = rel_scale_factor_CV[i_s];
	  //  temp_ss << "." << i_s;
	}
	if (variation_mu_fact_CV == 0){
	  //      scale_fact_CV[i_s] = scale_fact;
	  rel_scale_factor_fact_CV[i_s] = prefactor_CV / prefactor_reference;
	  //	  rel_scale_factor_fact_CV[i_s] = 1.;
	  //  temp_ss << "." << no_central_scale_CV;
	}
	if (variation_mu_fact_CV == -1){
	  //      scale_fact_CV[i_s] = rel_scale_factor_CV[j_s] * scale_fact;
	  rel_scale_factor_fact_CV[i_s] = rel_scale_factor_CV[j_s];
	  //  temp_ss << "." << j_s;
	}
      
	//    directory_name_scale_CV[i_s] = temp_ss.str();
	/*
	  mu_ren_CV[0][i_s] = 1.;
	  mu_ren_CV[1][i_s] = scale_ren_CV[i_s];
	  mu_ren_CV[2][i_s] = pow(scale_ren_CV[i_s], 2);
	  
	  mu_fact_CV[1][i_s] = scale_fact_CV[i_s];
	  mu_fact_CV[2][i_s] = pow(scale_fact_CV[i_s], 2);
	  mu_fact_CV[0][i_s] = log(mu_fact_CV[2][i_s]);
	*/
      }
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	logger << LOG_DEBUG << "TEST   rel_scale_factor_fact_CV[" << i_s << "] = " << rel_scale_factor_fact_CV[i_s] << endl;
      }
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	logger << LOG_DEBUG << "TEST   rel_scale_factor_ren_CV[" << i_s << "] = " << rel_scale_factor_ren_CV[i_s] << endl;
      }
    }
  
      // 7-pt variation
    else if (switch_CV == 5 && n_scales_CV != 1) {
      assert(n_scales_CV==7);
      
      rel_scale_factor_ren_CV[0] = 1.0/variation_factor_CV;
      rel_scale_factor_ren_CV[1] = 1.0/variation_factor_CV;
      rel_scale_factor_ren_CV[2] = variation_factor_CV;
      rel_scale_factor_ren_CV[3] = 1;
      rel_scale_factor_ren_CV[4] = 1;
      rel_scale_factor_ren_CV[5] = 1;
      rel_scale_factor_ren_CV[6] = variation_factor_CV;
      
      rel_scale_factor_fact_CV[0] = 1.0/variation_factor_CV;
      rel_scale_factor_fact_CV[1] = 1;
      rel_scale_factor_fact_CV[2] = 1;
      rel_scale_factor_fact_CV[3] = 1;
      rel_scale_factor_fact_CV[4] = variation_factor_CV;
      rel_scale_factor_fact_CV[5] = 1.0/variation_factor_CV;
      rel_scale_factor_fact_CV[6] = variation_factor_CV;

      
      //!!!
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	rel_scale_factor_ren_CV[i_s] = rel_scale_factor_ren_CV[i_s] * prefactor_CV / prefactor_reference;
	rel_scale_factor_fact_CV[i_s] = rel_scale_factor_fact_CV[i_s] * prefactor_CV / prefactor_reference;
      }
      
    }
    // 9-pt variation
    else if (switch_CV == 6 && n_scales_CV != 1) {
      assert(n_scales_CV==9);
      
      rel_scale_factor_ren_CV[0] = 1.0/variation_factor_CV;
      rel_scale_factor_ren_CV[1] = 1.0/variation_factor_CV;
      rel_scale_factor_ren_CV[2] = 1.0/variation_factor_CV;
      rel_scale_factor_ren_CV[3] = variation_factor_CV;
      rel_scale_factor_ren_CV[4] = 1;
      rel_scale_factor_ren_CV[5] = 1;
      rel_scale_factor_ren_CV[6] = 1;
      rel_scale_factor_ren_CV[7] = variation_factor_CV;
      rel_scale_factor_ren_CV[8] = variation_factor_CV;
      
      rel_scale_factor_fact_CV[0] = variation_factor_CV;
      rel_scale_factor_fact_CV[1] = 1.0/variation_factor_CV;
      rel_scale_factor_fact_CV[2] = 1;
      rel_scale_factor_fact_CV[3] = 1;
      rel_scale_factor_fact_CV[4] = 1;
      rel_scale_factor_fact_CV[5] = variation_factor_CV;
      rel_scale_factor_fact_CV[6] = 1.0/variation_factor_CV;
      rel_scale_factor_fact_CV[7] = variation_factor_CV;
      rel_scale_factor_fact_CV[8] = 1.0/variation_factor_CV;

      
      //!!!
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	rel_scale_factor_ren_CV[i_s] = rel_scale_factor_ren_CV[i_s] * prefactor_CV / prefactor_reference;
	rel_scale_factor_fact_CV[i_s] = rel_scale_factor_fact_CV[i_s] * prefactor_CV / prefactor_reference;
      }
      
    }
    else {
      logger << LOG_FATAL << "No valid CV settings. Aborted." << endl;
      exit(1);
    }

    //  if (switch_CV > 0) {
    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      
      scale_ren_CV[i_s] = rel_scale_factor_ren_CV[i_s] * scale_ren;
      scale_fact_CV[i_s] = rel_scale_factor_fact_CV[i_s] * scale_fact;
      
      mu_ren_CV[0][i_s] = 1.;
      mu_ren_CV[1][i_s] = scale_ren_CV[i_s];
      mu_ren_CV[2][i_s] = pow(scale_ren_CV[i_s], 2);

      mu_fact_CV[1][i_s] = scale_fact_CV[i_s];
      mu_fact_CV[2][i_s] = pow(scale_fact_CV[i_s], 2);
      mu_fact_CV[0][i_s] = log(mu_fact_CV[2][i_s]);
    }

    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      logger << LOG_DEBUG << "TEST   scale_fact_CV[" << i_s << "] = " << scale_fact_CV[i_s] << endl;
    }
    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      logger << LOG_DEBUG << "TEST   scale_ren_CV[" << i_s << "] = " << scale_ren_CV[i_s] << endl;
    }
    
  }
    
  logger << LOG_DEBUG << "mu_fact[1] = " << mu_fact[1] << endl;

  // !!! QCD_order == contribution_order_alpha_s removed -> csi->contribution_order_alpha_s

  logger << LOG_DEBUG << "csi->contribution_order_alpha_s = " << csi->contribution_order_alpha_s << endl;

  if (switch_CV){
    alpha_S_CV.resize(csi->contribution_order_alpha_s + 1, vector<double> (n_scales_CV));
    rel_alpha_S_CV.resize(csi->contribution_order_alpha_s + 1, vector<double> (n_scales_CV));
    
    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      alpha_S_CV[0][i_s] = 1.;
      rel_alpha_S_CV[0][i_s] = 1.;
      if (csi->contribution_order_alpha_s == 0){continue;}



    /*
    // !!!
    if (variation_mu_ren_CV == 0){
      alpha_S_CV[1][i_s] = alpha_S;
      for (int i_o = 2; i_o < alpha_S_CV.size(); i_o++){alpha_S_CV[i_o][i_s] = pow(alpha_S, i_o);}
      for (int i_o = 1; i_o < rel_alpha_S_CV.size(); i_o++){rel_alpha_S_CV[i_o][i_s] = 1.;}
    }
    else {
// !!!
    */
      if (variation_mu_ren_CV == 0){
	// !!!
	alpha_S_CV[1][i_s] = LHAPDF::alphasPDF(mu_ren_CV[1][i_s]);
	for (int i_o = 2; i_o < alpha_S_CV.size(); i_o++){alpha_S_CV[i_o][i_s] = pow(alpha_S_CV[1][i_s], i_o);}
	for (int i_o = 1; i_o < rel_alpha_S_CV.size(); i_o++){rel_alpha_S_CV[i_o][i_s] = alpha_S_CV[1][i_s] / alpha_S;}
	// !!!
      }
      else {
	if (variation_mu_ren_CV == 1){
	  alpha_S_CV[1][i_s] = LHAPDF::alphasPDF(mu_ren_CV[1][i_s]);
	}
	else if (variation_mu_ren_CV == -1){
	  int j_s = n_scales_CV - i_s - 1;
	  alpha_S_CV[1][i_s] = LHAPDF::alphasPDF(mu_ren_CV[1][j_s]);
	}
	for (int i_o = 2; i_o < alpha_S_CV.size(); i_o++){alpha_S_CV[i_o][i_s] = pow(alpha_S_CV[1][i_s], i_o);}
	rel_alpha_S_CV[1][i_s] = alpha_S_CV[1][i_s] / alpha_S;
	for (int i_o = 2; i_o < rel_alpha_S_CV.size(); i_o++){rel_alpha_S_CV[i_o][i_s] = pow(rel_alpha_S_CV[1][i_s], i_o);}
      }
    }
  


    // new scale calculation not yet tested !!! only output:
    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      if (csi->contribution_order_alpha_s == 0){continue;}
      logger << LOG_INFO << "alpha_S(" << setw(23) << setprecision(15) << mu_ren_CV[1][i_s] << ") = " << setw(23) << setprecision(15) << alpha_S_CV[1][i_s] << "   rel_alpha_S = " << setw(23) << setprecision(15) << rel_alpha_S_CV[1][i_s] << endl;
    }
    logger.newLine(LOG_INFO);

    map_value_scale_fact_CV.resize(n_scales_CV);
  }

  logger << LOG_DEBUG << "mu_fact[1] = " << mu_fact[1] << endl;


  initialization_scale_fact_CV();

  // !!! most likely only to have same vector size:
  value_mu_fact = value_mu_fact_rel;

  if (switch_CV){
    if (value_mu_fact[0].size() != 0){
      for (int i_s = 0; i_s < value_mu_fact[0].size(); i_s++){
	if (variation_mu_fact_CV == 0){value_mu_fact[0][i_s] = value_mu_fact_central[0];}
	else if (variation_mu_fact_CV == 1){value_mu_fact[0][i_s] = value_mu_fact_central[0] * value_mu_fact_rel[0][i_s];}
	else if (variation_mu_fact_CV == -1){value_mu_fact[0][i_s] = value_mu_fact_central[0] * value_mu_fact_rel[0][i_s];}
	//!!!      else if (variation_mu_fact_CV == -1){value_mu_fact[0][i_s] = value_mu_fact_central[0] / value_mu_fact_rel[0][i_s];}
	logger << LOG_DEBUG << "TEST   value_mu_fact[0][" << i_s << "] = " << value_mu_fact[0][i_s] << endl;
      }
    }

    map_value_scale_ren_CV.resize(n_scales_CV);
  }
  
  initialization_scale_ren_CV();
  // !!! most likely only to have same vector size:
  value_mu_ren = value_mu_ren_rel;
  value_alpha_S = value_mu_ren;
  value_factor_alpha_S = value_mu_ren;

  
  if (value_mu_ren[0].size() != 0){
    for (int i_s = 0; i_s < value_mu_ren[0].size(); i_s++){
      if (switch_CV){
	if (variation_mu_ren_CV == 0){value_mu_ren[0][i_s] = value_mu_ren_central[0];}
	else if (variation_mu_ren_CV == 1){value_mu_ren[0][i_s] = value_mu_ren_central[0] * value_mu_ren_rel[0][i_s];}
	else if (variation_mu_ren_CV == -1){value_mu_ren[0][i_s] = value_mu_ren_central[0] * value_mu_ren_rel[0][i_s];}
      }
      //!!!      else if (variation_mu_ren_CV == -1){value_mu_ren[0][i_s] = value_mu_ren_central[0] / value_mu_ren_rel[0][i_s];}
      value_alpha_S[0][i_s] = LHAPDF::alphasPDF(value_mu_ren[0][i_s]);
      value_factor_alpha_S[0][i_s] = pow(value_alpha_S[0][i_s] / alpha_S, contribution_order_alpha_s);

      logger << LOG_DEBUG << "TEST   value_mu_ren[0][" << i_s << "] = " << value_mu_ren[0][i_s] << "   value_alpha_S[0][" << i_s << "] = " << value_alpha_S[0][i_s] << endl;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "contribution_order_alpha_s = " << contribution_order_alpha_s << endl;


  var_mu_fact = mu_fact[1];
  logger << LOG_DEBUG << "mu_fact[1] = " << mu_fact[1] << endl;
  logger << LOG_DEBUG << "var_mu_fact = " << var_mu_fact << endl;
  var_mu_ren = mu_ren[1];
  var_alpha_S_reference = LHAPDF::alphasPDF(var_mu_ren);
  var_rel_alpha_S = 1.;
  
  if (switch_CV){
    // should be removed when ...value.. is completely installed! ???
    // used at all ???
    var_mu_ren_CV.resize(n_scales_CV);
    var_mu_fact_CV.resize(n_scales_CV);
    var_alpha_S_CV.resize(n_scales_CV);
    var_rel_alpha_S_CV.resize(n_scales_CV);

    var_mu_fact_CV = mu_fact_CV[1];
    var_mu_ren_CV = mu_ren_CV[1];
    var_alpha_S_CV = alpha_S_CV[contribution_order_alpha_s];
    var_rel_alpha_S_CV = rel_alpha_S_CV[contribution_order_alpha_s];
  }

  logger << LOG_DEBUG_VERBOSE << "xafter value_mu_ren[0].size() = " << value_mu_ren[0].size() << endl;


  pdf_factor.resize(3);
  
  if (switch_CV){
    pdf_factor_CV.resize(n_scales_CV, vector<double> (3));
  }


  Xsection = 0;
  Xsection_delta = 0;

  if (switch_CV){
    Xsection_CV.resize(n_qTcut, vector<double> (n_scales_CV, 0.));
    Xsection_delta_CV.resize(n_qTcut, vector<double> (n_scales_CV, 0.));
  }


  for (int i = 0; i < value_mu_fact_rel.size(); i++){
    logger << LOG_DEBUG << "value_mu_fact_central[" << setw(2) << i << "] = " << setw(25) << value_mu_fact_central[i] << endl;
    for (int j = 0; j < value_mu_fact_rel[i].size(); j++){
      logger << LOG_DEBUG << "value_mu_fact_rel[" << setw(2) << i << "][" << setw(2) << j << "] = " << setw(25) << value_mu_fact_rel[i][j] << endl;
    }
  }
  if (value_mu_fact[0].size() != 0){
    for (int j = 0; j < value_mu_fact_rel[0].size(); j++){
      logger << LOG_DEBUG << "value_mu_fact[" << setw(2) << 0 << "][" << setw(2) << j << "] = " << setw(25) << value_mu_fact[0][j] << endl;
    }
  }
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "map_value_scale_fact = " << map_value_scale_fact << endl;
  logger << LOG_DEBUG << endl;
  if (switch_CV){
    for (int s = 0; s < n_scales_CV; s++){
      logger << LOG_DEBUG << "map_value_scale_fact_CV[" << s << "] = " << map_value_scale_fact_CV[s] << endl;
    }
    logger << LOG_DEBUG << endl;
    logger.newLine(LOG_DEBUG);
  }

  for (int i = 0; i < value_mu_ren_rel.size(); i++){
    logger << LOG_DEBUG << "value_mu_ren_central[" << setw(2) << i << "] = " << setw(25) << value_mu_ren_central[i] << endl;
    for (int j = 0; j < value_mu_ren_rel[i].size(); j++){
      logger << LOG_DEBUG << "value_mu_ren_rel[" << setw(2) << i << "][" << setw(2) << j << "] = " << setw(25) << value_mu_ren_rel[i][j] << endl;
    }
  }
  if (value_mu_ren[0].size() != 0){
    for (int j = 0; j < value_mu_ren_rel[0].size(); j++){
      logger << LOG_DEBUG << "value_mu_ren[" << setw(2) << 0 << "][" << setw(2) << j << "] = " << setw(25) << value_mu_ren[0][j] << endl;
    }
    for (int j = 0; j < value_mu_ren_rel[0].size(); j++){
      logger << LOG_DEBUG << "value_alpha_S[" << setw(2) << 0 << "][" << setw(2) << j << "] = " << setw(25) << value_alpha_S[0][j] << endl;
    }
    for (int j = 0; j < value_mu_ren_rel[0].size(); j++){
      logger << LOG_DEBUG << "value_factor_alpha_S[" << setw(2) << 0 << "][" << setw(2) << j << "] = " << setw(25) << value_factor_alpha_S[0][j] << endl;
    }
  }
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "map_value_scale_ren = " << map_value_scale_ren << endl;
  logger << LOG_DEBUG << endl;
  if (switch_CV){
    for (int s = 0; s < n_scales_CV; s++){
      logger << LOG_DEBUG << "map_value_scale_ren_CV[" << s << "] = " << map_value_scale_ren_CV[s] << endl;
    }
    logger << LOG_DEBUG << endl;
  }
  logger << LOG_DEBUG << "id_scales = " << id_scales << endl;


  logger << LOG_DEBUG << "finished" << endl;
}



void observable_set::initialization_scale_fact_CV(){
  Logger logger("observable_set::initialization_scale_fact_CV");
  logger << LOG_DEBUG << "started" << endl;

  int max_dyn = 0;
  if (dynamic_scale > max_dyn){max_dyn = dynamic_scale;}

  if (switch_CV){if (dynamic_scale_CV > max_dyn){max_dyn = dynamic_scale_CV;}}

  value_mu_fact_rel.resize(max_dyn + 1);
  value_mu_fact_central.resize(max_dyn + 1);
  value_mu_fact_central[dynamic_scale] = mu_fact[1];

  // !!!  if (switch_CV){value_mu_fact_central[dynamic_scale_CV] = mu_fact_CV[1][no_central_scale_CV];}
  // !!!
  if (switch_CV){
    if (dynamic_scale != dynamic_scale_CV){
      value_mu_fact_central[dynamic_scale_CV] = mu_fact_CV[1][no_central_scale_CV];
    }
  }
  // !!!

  logger << LOG_DEBUG << "TEST   value_mu_fact_central[" << dynamic_scale << "] = " << value_mu_fact_central[dynamic_scale] << endl;
  if (switch_CV){
    logger << LOG_DEBUG << "TEST   CV   value_mu_fact_central[" << dynamic_scale_CV << "] = " << value_mu_fact_central[dynamic_scale_CV] << endl;
  }

  logger << LOG_DEBUG << "mu_fact[1] = " << mu_fact[1] << endl;
  for (int i = 0; i < value_mu_fact_central.size(); i++){
    logger << LOG_DEBUG << "value_mu_fact_central[" << dynamic_scale << "] = " << value_mu_fact_central[dynamic_scale] << endl;
  }
  vector<double> rel_scale_CV(n_scales_CV);
  int flag = 0;
  logger << LOG_DEBUG << "TEST   value_mu_fact_rel[" << dynamic_scale << "].size() = " << value_mu_fact_rel[dynamic_scale].size() << endl;
  for (int i = 0; i < value_mu_fact_rel[dynamic_scale].size(); i++){
    if (abs(value_mu_fact_rel[dynamic_scale][i] - mu_fact[1] / value_mu_fact_central[dynamic_scale]) < 1.e-15){flag = 1; break;}
  }
  if (flag == 0){value_mu_fact_rel[dynamic_scale].push_back(mu_fact[1] / value_mu_fact_central[dynamic_scale]);}

  logger << LOG_DEBUG << "TEST   value_mu_fact_rel[" << dynamic_scale << "].size() = " << value_mu_fact_rel[dynamic_scale].size() << endl;
  for (int i_s = 0; i_s < value_mu_fact_rel[dynamic_scale].size(); i_s++){
    logger << LOG_DEBUG << "TEST   value_mu_fact_rel[" << dynamic_scale << "][" << i_s << "] = " << value_mu_fact_rel[dynamic_scale][i_s] << endl;
  }

  if (switch_CV){
    for (int s = 0; s < n_scales_CV;s++){
      logger << LOG_DEBUG << "CV TEST   s = " << s << "   value_mu_fact_rel[" << dynamic_scale_CV << "].size() = " << value_mu_fact_rel[dynamic_scale_CV].size() << endl;
      flag = 0;
      //!!!      if (n_scales_CV > 1){rel_scale_CV[s] = pow(10., log10(double(variation_factor_CV)) * double(2. * s / (n_scales_CV - 1) - 1));}
      if (switch_CV == 5){
	rel_scale_CV[0] = 1. / variation_factor_CV;
	rel_scale_CV[1] = 1.;
	rel_scale_CV[2] = 1.;
	rel_scale_CV[3] = 1.;
	rel_scale_CV[4] = variation_factor_CV;
	rel_scale_CV[5] = 1. / variation_factor_CV;
	rel_scale_CV[6] = variation_factor_CV;
      }
      else if (switch_CV == 6){
	rel_scale_CV[0] = variation_factor_CV;
	rel_scale_CV[1] = 1. / variation_factor_CV;
	rel_scale_CV[2] = 1.;
	rel_scale_CV[3] = 1.;
	rel_scale_CV[4] = 1.;
	rel_scale_CV[5] = variation_factor_CV;
	rel_scale_CV[6] = 1. / variation_factor_CV;
	rel_scale_CV[7] = variation_factor_CV;
	rel_scale_CV[8] = 1. / variation_factor_CV;
      }
      // !!!
      else if (variation_mu_fact_CV == 0){
	rel_scale_CV[s] = 1.;
      }
      else if (variation_mu_fact_CV == -1){
	if (n_scales_CV > 1){rel_scale_CV[s] = pow(10., log10(double(variation_factor_CV)) * double(2. * (n_scales_CV - 1 - s) / (n_scales_CV - 1) - 1));}
	else {rel_scale_CV[s] = 1.;}
      }
      else {
	if (n_scales_CV > 1){rel_scale_CV[s] = pow(10., log10(double(variation_factor_CV)) * double(2. * s / (n_scales_CV - 1) - 1));}
	else {rel_scale_CV[s] = 1.;}
      }
      /*
      //!!!
      else {
	if (n_scales_CV > 1){rel_scale_CV[s] = pow(10., log10(double(variation_factor_CV)) * double(2. * s / (n_scales_CV - 1) - 1));}
	else {rel_scale_CV[s] = 1.;}
      }
      */

      for (int i = 0; i < value_mu_fact_rel[dynamic_scale_CV].size(); i++){
	if (abs(value_mu_fact_rel[dynamic_scale_CV][i] / rel_scale_CV[s] - mu_fact_CV[1][no_central_scale_CV] / value_mu_fact_central[dynamic_scale_CV]) < 1.e-15){flag = 1; break;}
      }
      if (flag == 0){value_mu_fact_rel[dynamic_scale_CV].push_back(rel_scale_CV[s] * mu_fact_CV[1][no_central_scale_CV] / value_mu_fact_central[dynamic_scale_CV]);}
      logger << LOG_DEBUG << "CV TEST   s = " << s << "   value_mu_fact_rel[" << dynamic_scale_CV << "].size() = " << value_mu_fact_rel[dynamic_scale_CV].size() << endl;
    }

    for (int i_s = 0; i_s < value_mu_fact_rel[dynamic_scale_CV].size(); i_s++){
      logger << LOG_DEBUG << "TEST   value_mu_fact_rel[" << dynamic_scale_CV << "][" << i_s << "] = " << value_mu_fact_rel[dynamic_scale_CV][i_s] << endl;
    }

  }



  for (int j = 0; j < max_dyn + 1; j++){sort (value_mu_fact_rel[j].begin(), value_mu_fact_rel[j].end());}

  if (switch_CV){
    for (int i_s = 0; i_s < value_mu_fact_rel[dynamic_scale_CV].size(); i_s++){
      logger << LOG_DEBUG << "TEST   value_mu_fact_rel[" << dynamic_scale_CV << "][" << i_s << "] = " << value_mu_fact_rel[dynamic_scale_CV][i_s] << endl;
    }
  }

  for (int i = 0; i < value_mu_fact_rel[dynamic_scale].size(); i++){
    if (abs(value_mu_fact_rel[dynamic_scale][i] - mu_fact[1] / value_mu_fact_central[dynamic_scale]) < 1.e-15){map_value_scale_fact = i;}
  }
  logger << LOG_DEBUG << "TEST   dynamic_scale = " << dynamic_scale << "   map_value_scale_fact = " << map_value_scale_fact << endl;



  if (switch_CV){
    for (int s = 0; s < n_scales_CV;s++){
      for (int i = 0; i < value_mu_fact_rel[dynamic_scale_CV].size(); i++){
	if (abs(value_mu_fact_rel[dynamic_scale_CV][i] / rel_scale_CV[s] - mu_fact_CV[1][no_central_scale_CV] / value_mu_fact_central[dynamic_scale_CV]) < 1.e-15){map_value_scale_fact_CV[s] = i;}
      }
      logger << LOG_DEBUG << "TEST   dynamic_scale_CV = " << dynamic_scale_CV << "   map_value_scale_fact_CV[" << s << "] = " << map_value_scale_fact_CV[s] << endl;
    }
  }


  
  // to be checked later today !!!
  /*
  // 7-pt variation
  if (switch_CV == 5) {
//    for (int sd = 1; sd < value_mu_fact_rel.size(); sd++){
      int sd = dynamic_scale_CV;
      // !!!      int sd = dynamic_scale;
      assert(value_mu_fact_rel[sd].size()==7 || value_mu_fact_rel[sd].size()==1);
      if (value_mu_fact_rel[sd].size()==7) {
        value_mu_fact_rel[sd][0] = 1.0/variation_factor_CV;
        value_mu_fact_rel[sd][1] = 1;
        value_mu_fact_rel[sd][2] = 1;
        value_mu_fact_rel[sd][3] = 1;
        value_mu_fact_rel[sd][4] = variation_factor_CV;
        value_mu_fact_rel[sd][5] = 1.0/variation_factor_CV;
        value_mu_fact_rel[sd][6] = variation_factor_CV;
      }
      /*
 else {
        // for warmup runs
        value_mu_fact_rel[sd][0] = 1;
      }
  *//*
  //  }
  }
  // 9-pt variation
  if (switch_CV == 6) {
    //for (int sd = 1; sd < value_mu_fact_rel.size(); sd++){
      int sd = dynamic_scale_CV;
      // !!!      int sd = dynamic_scale;
      assert(value_mu_fact_rel[sd].size()==9 || value_mu_fact_rel[sd].size()==1);
      if (value_mu_fact_rel[sd].size()==9) {
        value_mu_fact_rel[sd][0] = variation_factor_CV;
        value_mu_fact_rel[sd][1] = 1.0/variation_factor_CV;
        value_mu_fact_rel[sd][2] = 1;
        value_mu_fact_rel[sd][3] = 1;
        value_mu_fact_rel[sd][4] = 1;
        value_mu_fact_rel[sd][5] = variation_factor_CV;
        value_mu_fact_rel[sd][6] = 1.0/variation_factor_CV;
        value_mu_fact_rel[sd][7] = variation_factor_CV;
        value_mu_fact_rel[sd][8] = 1.0/variation_factor_CV;
      }/*
 else {
        // for warmup runs
        value_mu_fact_rel[sd][0] = 1;
      }
    *//*
    //}
  }
      */

  /*
  // !!!
  // not needed for static scales, since central value of old (?) implementation is already scaled
  if (dynamic_scale_CV) {
    for (int ss=0; ss < value_mu_fact_rel[dynamic_scale].size(); ss++) {
      value_mu_fact_rel[dynamic_scale][ss] *= prefactor_CV / prefactor_reference;
      // ???
    }
  }
  */
  

  if (dynamic_scale){// && !dynamic_scale_CV) {
    for (int ss = 0; ss < value_mu_fact_rel[dynamic_scale].size(); ss++) {
      value_mu_fact_rel[dynamic_scale][ss] *= prefactor_reference;
      // ???
    }
  }

  logger << LOG_DEBUG << "dynamic_scale_CV = " << dynamic_scale_CV << endl;
  if (!dynamic_scale && (switch_CV && dynamic_scale_CV)){
    for (int ss = 0; ss < value_mu_fact_rel[dynamic_scale_CV].size(); ss++) {
      value_mu_fact_rel[dynamic_scale_CV][ss] *= prefactor_CV;
      // ???
    }
  }
  logger << LOG_DEBUG << "dynamic_scale_CV = " << dynamic_scale_CV << endl;
 

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_scale_ren_CV(){
  Logger logger("observable_set::initialization_scale_ren_CV");
  logger << LOG_DEBUG << "started" << endl;

  int max_dyn = 0;
  if (dynamic_scale > max_dyn){max_dyn = dynamic_scale;}

  if (switch_CV){if (dynamic_scale_CV > max_dyn){max_dyn = dynamic_scale_CV;}}

  value_mu_ren_rel.resize(max_dyn + 1);
  value_mu_ren_central.resize(max_dyn + 1);
  value_mu_ren_central[dynamic_scale] = mu_ren[1];

  // !!!  if (switch_CV){value_mu_ren_central[dynamic_scale_CV] = mu_ren_CV[1][no_central_scale_CV];}
  // !!!
  if (switch_CV){
    if (dynamic_scale != dynamic_scale_CV){
      value_mu_ren_central[dynamic_scale_CV] = mu_ren_CV[1][no_central_scale_CV];
    }
  }
  // !!!
  
  logger << LOG_DEBUG_VERBOSE << "TEST   value_mu_ren_central[" << dynamic_scale << "] = " << value_mu_ren_central[dynamic_scale] << endl;
  if (switch_CV){
    logger << LOG_DEBUG_VERBOSE << "TEST   CV   value_mu_ren_central[" << dynamic_scale_CV << "] = " << value_mu_ren_central[dynamic_scale_CV] << endl;
  }
  


  // Could be a member variable...
  vector<double> rel_scale_CV;
  
  int flag = 0;
  logger << LOG_DEBUG_VERBOSE << "TEST   value_mu_ren_rel[" << dynamic_scale << "].size() = " << value_mu_ren_rel[dynamic_scale].size() << endl;
  for (int i = 0; i < value_mu_ren_rel[dynamic_scale].size(); i++){
    if (abs(value_mu_ren_rel[dynamic_scale][i] - mu_ren[1] / value_mu_ren_central[dynamic_scale]) < 1.e-15){flag = 1; break;}
  }
  if (flag == 0){value_mu_ren_rel[dynamic_scale].push_back(mu_ren[1] / value_mu_ren_central[dynamic_scale]);}
  
  logger << LOG_DEBUG_VERBOSE << "TEST   value_mu_ren_rel[" << dynamic_scale << "].size() = " << value_mu_ren_rel[dynamic_scale].size() << endl;
  for (int i_s = 0; i_s < value_mu_ren_rel[dynamic_scale].size(); i_s++){
    logger << LOG_DEBUG_VERBOSE << "TEST   value_mu_ren_rel[" << dynamic_scale << "][" << i_s << "] = " << value_mu_ren_rel[dynamic_scale][i_s] << endl;
  }
  
  if (switch_CV){
    rel_scale_CV.resize(n_scales_CV);
    for (int s = 0; s < n_scales_CV;s++){
      logger << LOG_DEBUG_VERBOSE << "CV TEST   s = " << s << "   value_mu_ren_rel[" << dynamic_scale_CV << "].size() = " << value_mu_ren_rel[dynamic_scale_CV].size() << endl;
      flag = 0;
      if (switch_CV == 5){
        rel_scale_CV[0] = 1. / variation_factor_CV;
        rel_scale_CV[1] = 1. / variation_factor_CV;
        rel_scale_CV[2] = variation_factor_CV;
        rel_scale_CV[3] = 1.;
        rel_scale_CV[4] = 1.;
        rel_scale_CV[5] = 1.;
        rel_scale_CV[6] = variation_factor_CV;
      }
      else if (switch_CV == 6){
        rel_scale_CV[0] = 1. / variation_factor_CV;
        rel_scale_CV[1] = 1. / variation_factor_CV;
        rel_scale_CV[2] = 1. / variation_factor_CV;
        rel_scale_CV[3] = variation_factor_CV;
        rel_scale_CV[4] = 1.;
        rel_scale_CV[5] = 1.;
        rel_scale_CV[6] = 1.;
        rel_scale_CV[7] = variation_factor_CV;
        rel_scale_CV[8] = variation_factor_CV;
      }
      // !!!
      else if (variation_mu_ren_CV == 0){
	rel_scale_CV[s] = 1.;
      }
      else if (variation_mu_ren_CV == -1){
	if (n_scales_CV > 1){rel_scale_CV[s] = pow(10., log10(double(variation_factor_CV)) * double(2. * (n_scales_CV - 1 - s) / (n_scales_CV - 1) - 1));}
	else {rel_scale_CV[s] = 1.;}
      }
      else {
	if (n_scales_CV > 1){rel_scale_CV[s] = pow(10., log10(double(variation_factor_CV)) * double(2. * s / (n_scales_CV - 1) - 1));}
	else {rel_scale_CV[s] = 1.;}
      }

      for (int i = 0; i < value_mu_ren_rel[dynamic_scale_CV].size(); i++){
	if (abs(value_mu_ren_rel[dynamic_scale_CV][i] / rel_scale_CV[s] - mu_ren_CV[1][no_central_scale_CV] / value_mu_ren_central[dynamic_scale_CV]) < 1.e-15){flag = 1; break;}
      }
      if (flag == 0){value_mu_ren_rel[dynamic_scale_CV].push_back(rel_scale_CV[s] * mu_ren_CV[1][no_central_scale_CV] / value_mu_ren_central[dynamic_scale_CV]);}
      logger << LOG_DEBUG_VERBOSE << "CV TEST   s = " << s << "   value_mu_ren_rel[" << dynamic_scale_CV << "].size() = " << value_mu_ren_rel[dynamic_scale_CV].size() << endl;
    }
    
    for (int i_s = 0; i_s < value_mu_ren_rel[dynamic_scale_CV].size(); i_s++){
      logger << LOG_DEBUG_VERBOSE << "TEST   value_mu_ren_rel[" << dynamic_scale_CV << "][" << i_s << "] = " << value_mu_ren_rel[dynamic_scale_CV][i_s] << endl;
    }
    
  }


  
  for (int j = 0; j < max_dyn + 1; j++){sort (value_mu_ren_rel[j].begin(), value_mu_ren_rel[j].end());}
  for (int i = 0; i < value_mu_ren_rel[dynamic_scale].size(); i++){
    if (abs(value_mu_ren_rel[dynamic_scale][i] - mu_ren[1] / value_mu_ren_central[dynamic_scale]) < 1.e-15){map_value_scale_ren = i;}
  }
  logger << LOG_DEBUG_VERBOSE << "TEST   dynamic_scale = " << dynamic_scale << "   map_value_scale_ren = " << map_value_scale_ren << endl;



  if (switch_CV){
    for (int s = 0; s < n_scales_CV;s++){
      for (int i = 0; i < value_mu_ren_rel[dynamic_scale_CV].size(); i++){
	if (abs(value_mu_ren_rel[dynamic_scale_CV][i] / rel_scale_CV[s] - mu_ren_CV[1][no_central_scale_CV] / value_mu_ren_central[dynamic_scale_CV]) < 1.e-15){map_value_scale_ren_CV[s] = i;}
      }
      logger << LOG_DEBUG_VERBOSE << "TEST   dynamic_scale_CV = " << dynamic_scale_CV << "   map_value_scale_ren_CV[" << s << "] = " << map_value_scale_ren_CV[s] << endl;
    }
  }
  

  // to be checked later today !!!
  /*
  // 7-pt variation
  if (switch_CV == 5) {
    //for (int sd = 1; sd < value_mu_ren_rel.size(); sd++){
      int sd = dynamic_scale_CV;
      //      int sd = dynamic_scale;
      assert(value_mu_ren_rel[sd].size()==7 || value_mu_ren_rel[sd].size()==1);
      if (value_mu_ren_rel[sd].size()==7) {
        value_mu_ren_rel[sd][0] = 1.0/variation_factor_CV;
        value_mu_ren_rel[sd][1] = 1.0/variation_factor_CV;
        value_mu_ren_rel[sd][2] = variation_factor_CV;
        value_mu_ren_rel[sd][3] = 1;
        value_mu_ren_rel[sd][4] = 1;
        value_mu_ren_rel[sd][5] = 1;
        value_mu_ren_rel[sd][6] = variation_factor_CV;
      }
      /*
 else {
        // for warmup runs
        value_mu_ren_rel[sd][0] = 1;
      }
  *//*
    //}
  }
  // 9-pt variation
  if (switch_CV == 6) {
    //for (int sd = 1; sd < value_mu_ren_rel.size(); sd++){
      int sd = dynamic_scale_CV;
      //      int sd = dynamic_scale;
      assert(value_mu_ren_rel[sd].size()==9 || value_mu_ren_rel[sd].size()==1);
      if (value_mu_ren_rel[sd].size()==9) {
        value_mu_ren_rel[sd][0] = 1.0/variation_factor_CV;
        value_mu_ren_rel[sd][1] = 1.0/variation_factor_CV;
        value_mu_ren_rel[sd][2] = 1.0/variation_factor_CV;
        value_mu_ren_rel[sd][3] = variation_factor_CV;
        value_mu_ren_rel[sd][4] = 1;
        value_mu_ren_rel[sd][5] = 1;
        value_mu_ren_rel[sd][6] = 1;
        value_mu_ren_rel[sd][7] = variation_factor_CV;
        value_mu_ren_rel[sd][8] = variation_factor_CV;
      }
      /*
 else {
        // for warmup runs
        value_mu_ren_rel[sd][0] = 1;
      }
    *//*
    //}
  }
      */

  if (dynamic_scale){// && !dynamic_scale_CV) {
    for (int ss = 0; ss < value_mu_ren_rel[dynamic_scale].size(); ss++) {
      value_mu_ren_rel[dynamic_scale][ss] *= prefactor_reference;
      // ???
    }
  }
  if (!dynamic_scale && (switch_CV && dynamic_scale_CV)){
    for (int ss = 0; ss < value_mu_ren_rel[dynamic_scale_CV].size(); ss++) {
      value_mu_ren_rel[dynamic_scale_CV][ss] *= prefactor_CV;
      // ???
    }
  }
 
  /*
  // !!!
  // not needed for static scales, since central value of old (?) implementation is already scaled
  if (dynamic_scale_CV) {
    for (int ss=0; ss < value_mu_ren_rel[dynamic_scale].size(); ss++) {
      value_mu_ren_rel[dynamic_scale][ss] *= prefactor_CV / prefactor_reference;
      // ???
    }
  }
  */
  logger << LOG_DEBUG << "finished" << endl;
}



