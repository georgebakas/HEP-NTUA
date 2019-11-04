#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_subprocess::readin_result_TSV(string & result_moment){
  Logger logger("summary_subprocess::readin_result_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // shifted here in order to reduze temporary memory consumption
  if (ycontribution->output_n_qTcut < osi->n_qTcut){
    result_run_TSV.resize(n_seed, vector<vector<vector<vector<vector<double> > > > > (osi->n_moments + 1));
    deviation_run_TSV.resize(n_seed, vector<vector<vector<vector<vector<double> > > > > (osi->n_moments + 1));
    for (int i_z = 0; i_z < n_seed; i_z++){
      for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
	result_run_TSV[i_z][i_m].resize(ycontribution->output_n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV));
	deviation_run_TSV[i_z][i_m].resize(ycontribution->output_n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV));
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  result_run_TSV[i_z][i_m][ycontribution->output_n_qTcut - 1][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  deviation_run_TSV[i_z][i_m][ycontribution->output_n_qTcut - 1][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	}
      }
    }
  }
  else {
    result_run_TSV.resize(n_seed, result_TSV);
    deviation_run_TSV.resize(n_seed, deviation_TSV);
  }
    
  int x_m = 0;
  for (int i_z = 0; i_z < n_seed; i_z++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      string result_file = result_moment + "_" + name + ".txt";
      logger << LOG_DEBUG << "TSV: result_file = " << result_file << endl;
      string filename = "../" + ycontribution->directory[i_z] + "/result/" + osi->name_extended_set_TSV[i_s] + "/" + result_file;
      logger << LOG_DEBUG << "TSV: filename = " << filename << endl;
      
      vector<string> readin;
      vector<vector<string> > readin_data;
      char LineBuffer[128];
      ifstream in_result(filename.c_str());  
      while (in_result.getline(LineBuffer, 128)){readin.push_back(LineBuffer);}
      in_result.close();
      for (int i = 0; i < readin.size(); i++){
	vector<string> temp_s(1);
	int temp_comment = 0;
	for (int j_s = 0; j_s < readin[i].size(); j_s++){
	  if (readin[i][0] == '#'){temp_comment = 1; break;}
	  else if (readin[i][j_s] == ' ' || readin[i][j_s] == char(9)){
	    if (temp_s[temp_s.size() - 1].size() == 0){}
	    else{temp_s.push_back("");}
	  }
	  else {temp_s[temp_s.size() - 1].push_back(readin[i][j_s]);}
	}
	if (temp_comment == 0){
	  if (temp_s[temp_s.size() - 1] == ""){temp_s.erase(temp_s.end());}
	  readin_data.push_back(temp_s);
	}
      }
      
      int default_size = 1 + ycontribution->output_n_qTcut * (osi->n_scale_ren_TSV[i_s] * osi->n_scale_fact_TSV[i_s]);

      if (default_size == readin_data.size()){
	N_event_run_TSV[i_z] = atoll(readin_data[0][0].c_str());
	for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
	  int proc = i_q * osi->n_scale_ren_TSV[i_s] * osi->n_scale_fact_TSV[i_s];
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] = atof(readin_data[proc + i_r * osi->n_scale_fact_TSV[i_s] + i_f + 1][4].c_str());
	      deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] = atof(readin_data[proc + i_r * osi->n_scale_fact_TSV[i_s] + i_f + 1][5].c_str());
	      logger << LOG_DEBUG_VERBOSE << "result_run_TSV[" << i_z << "][" << x_m << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << setw(23) << setprecision(15) << result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << " +- = " << setw(23) << setprecision(15) << deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << endl;
	    }
	  }
	}
      }
      else {
	logger << LOG_DEBUG << "TSV: result_file = " << filename << " has not the expected content:   default_size = " << default_size << " != " << readin_data.size() << " = readin_data.size()" << endl;
      }
      int x_q = 0;
      int x_r = 0;
      int x_f = 0;
      logger << LOG_DEBUG << "Readin of " << i_z << " done:   " << "N_event_run_TSV[" << i_z << "] = " << setw(15) << N_event_run_TSV[i_z] << "   result_run_TSV[" << i_z << "][" << x_m << "][" << x_q << "][" << i_s << "][" << x_r << "][" << x_f << "] = " << setw(23) << setprecision(15) << result_run_TSV[i_z][x_m][x_q][i_s][x_r][x_f] << " +- = " << setw(23) << setprecision(15) << deviation_run_TSV[i_z][x_m][x_q][i_s][x_r][x_f] << endl;
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combination_result_TSV(string & result_moment){
  Logger logger("summary_subprocess::combination_result_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG << left << name << endl;

  int x_m = 0;

  // 0 replaced by reference value:
  int x_s = osi->no_reference_TSV;
  int x_r = (osi->n_scale_ren_TSV[x_s] - 1) / 2;
  int x_f = (osi->n_scale_fact_TSV[x_s] - 1) / 2;

  for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
    n_parallel_runs_reference_TSV[i_q] = 0;
    for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
      if (N_event_run_TSV[i_z] != 0 && 
	  result_run_TSV[i_z][x_m][i_q][x_s][x_r][x_f] != 0. && 
	  deviation_run_TSV[i_z][x_m][i_q][x_s][x_r][x_f] != 0.){
	n_parallel_runs_reference_TSV[i_q]++;
	if (!remove_run_qTcut_TSV[i_q][i_z]){
	  n_valid_parallel_runs_reference_TSV[i_q]++;
	}
      }
    }
  }

  removal_exceptional_runs_result_TSV();

  if (ygeneric->average_factor == 0){combine_result_alternative_TSV();}  // removal of exceptional runs not included !!!
  else if (ygeneric->average_factor == 1){combine_result_original_TSV();}  // removal of exceptional runs included
  else {combine_result_hybrid_TSV();}  // removal of exceptional runs not included !!!

  // higher qTcut-values are identified with i_q = 0 (for qTcut independent runs); could be improved...
  if (ycontribution->output_n_qTcut < osi->n_qTcut){
    for (int i_q = ycontribution->output_n_qTcut; i_q < osi->n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    result_TSV[x_m][i_q][i_s][i_r][i_f] = result_TSV[x_m][0][i_s][i_r][i_f];
	    deviation_TSV[x_m][i_q][i_s][i_r][i_f] = deviation_TSV[x_m][0][i_s][i_r][i_f];
	  }
	}
      }
    }
  }

  //  for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
  for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
	    if (N_event_run_TSV[i_z] != 0 && 
		result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] != 0. && 
		deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] != 0.){
	      n_parallel_runs_TSV[x_m][i_q][i_s][i_r][i_f]++;
	      if (!remove_run_qTcut_TSV[i_q][i_z]){
		n_valid_parallel_runs_TSV[x_m][i_q][i_s][i_r][i_f]++;
	      }
	    }
	  }

	  chi2_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	  if (n_parallel_runs_TSV[x_m][i_q][i_s][i_r][i_f] > 1){
	    for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
	      if (!(result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0. && deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0.)){
		chi2_TSV[x_m][i_q][i_s][i_r][i_f] += pow((result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] - result_TSV[x_m][i_q][i_s][i_r][i_f]) / deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f], 2.);
	      }
	    }
	    chi2_TSV[x_m][i_q][i_s][i_r][i_f] = chi2_TSV[x_m][i_q][i_s][i_r][i_f] / n_parallel_runs_TSV[x_m][i_q][i_s][i_r][i_f];
	  }
	}
      }
    }
  }

  // new location: to clear result_run_TSV before runtime determination:
  // store reference result in result_run !
  logger << LOG_DEBUG << "osi->name_reference_TSV = " << osi->name_reference_TSV << endl;
  
  /*
  int x_m = 0;
  // 0 replaced by reference value:
  int x_s = osi->no_reference_TSV;
  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    logger << LOG_DEBUG << "osi->name_extended_set_TSV[" << i_s << "] = " << osi->name_extended_set_TSV[i_s] << endl;
    if (osi->name_reference_TSV == osi->name_extended_set_TSV[i_s]){x_s = i_s; break;}
  }
  int x_r = osi->no_scale_ren_reference_TSV;
  int x_f = osi->no_scale_fact_reference_TSV;
  */
  int x_q = 0;
  //osi->no_qTcut_reference_TSV;
  if (ycontribution->active_qTcut && ygeneric->no_qTcut_runtime_estimate){
    x_q = ygeneric->no_qTcut_runtime_estimate;
  }
  logger << LOG_DEBUG << "x_m = " << x_m << endl;
  logger << LOG_DEBUG << "x_q = " << x_q << endl;
  logger << LOG_DEBUG << "x_s = " << x_s << endl;
  logger << LOG_DEBUG << "x_r = " << x_r << endl;
  logger << LOG_DEBUG << "x_f = " << x_f << endl;
  
  for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
    result_run[i_z] = result_run_TSV[i_z][x_m][x_q][x_s][x_r][x_f];
    deviation_run[i_z] = deviation_run_TSV[i_z][x_m][x_q][x_s][x_r][x_f];
  }
  // untile here.

  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::removal_exceptional_runs_result_TSV(){
  Logger logger("summary_subprocess::removal_exceptional_runs_result_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int x_m = 0;
  N_event_TSV = accumulate(N_event_run_TSV.begin(), N_event_run_TSV.end(), 0);
  Nd_event_TSV = 0.;
  for (int i_z = 0; i_z < n_seed; i_z++){
    Nd_event_TSV += double(N_event_run_TSV[i_z]);
  }

  vector<vector<double> > deviation_measure_qTcut(ycontribution->output_n_qTcut, vector<double> (n_seed, 0.));

  //  if (n_parallel_runs_reference_TSV[0] > 0 && N_event_TSV > 0){
  if (n_parallel_runs_reference_TSV[0] > 0 && Nd_event_TSV > 0.){
    for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      if (ygeneric->deviation_tolerance_factor != 0.){
	// remove runs that indicate numerical instabilities from the combination procedure
	// (separately for each qTcut, but not for all different scales)
	double tolerance_factor = ygeneric->deviation_tolerance_factor;
	// 0 replaced by reference value:
	int x_s = osi->no_reference_TSV;
	int x_r = (osi->n_scale_ren_TSV[x_s] - 1) / 2;
	int x_f = (osi->n_scale_fact_TSV[x_s] - 1) / 2;
	vector<double> deviation_measure(n_seed, 0.);
	double min_deviation_measure = 1.e99;
	double max_deviation_measure = 0.;
	for (int i_z = 0; i_z < n_seed; i_z++){
	  logger << LOG_DEBUG << "result_run_TSV[i_z][x_m][i_q][x_s][x_r][x_f] = " << setw(23) << setprecision(15) << result_run_TSV[i_z][x_m][i_q][x_s][x_r][x_f] << " +- " << setw(23) << setprecision(15) << deviation_run_TSV[i_z][x_m][i_q][x_s][x_r][x_f] << "   ---   " << N_event_run_TSV[i_z] << endl;
	  if (N_event_run_TSV[i_z] > 0){
	    deviation_measure[i_z] = deviation_run_TSV[i_z][x_m][i_q][x_s][x_r][x_f] * sqrt(N_event_run_TSV[i_z]);
	  }
	}
	for (int i_z = 0; i_z < n_seed; i_z++){
	  if (N_event_run_TSV[i_z] > 0 && deviation_measure[i_z] < min_deviation_measure){min_deviation_measure = deviation_measure[i_z];}
	  if (N_event_run_TSV[i_z] > 0 && deviation_measure[i_z] > max_deviation_measure){max_deviation_measure = deviation_measure[i_z];}
	}
	N_valid_event_TSV[i_q] = 0;
	for (int i_z = 0; i_z < n_seed; i_z++){
	  deviation_measure_qTcut[i_q][i_z] = deviation_measure[i_z] / min_deviation_measure;
	  if (deviation_measure[i_z] > tolerance_factor * min_deviation_measure){remove_run_qTcut_TSV[i_q][i_z] = 1;}
	  else {N_valid_event_TSV[i_q] += N_event_run_TSV[i_z];}
	}
      }
    }
  }

  /*
  for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
    logger << LOG_INFO << setw(20) << name << "   N_event_TSV = " << setw(12) << N_event_TSV << "   N_valid_event_TSV[i_q] = " << setw(12) << N_valid_event_TSV[i_q] << "   N_event_TSV = " << setw(12) << N_event_TSV << endl;
      }
  */

  if (ygeneric->deviation_tolerance_factor != 0.){
    // information output about removed runs
    for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      if (accumulate(remove_run_qTcut_TSV[i_q].begin(), remove_run_qTcut_TSV[i_q].end(), 0)){
	stringstream temp_ss;
	temp_ss << setw(8) << ycontribution->type_contribution << "   " << setw(15) << name << "   i_q = " << setw(3) << i_q << "   erased: " << setw(3) << accumulate(remove_run_qTcut_TSV[i_q].begin(), remove_run_qTcut_TSV[i_q].end(), 0) << " (" << setw(4) << n_parallel_runs_reference_TSV[i_q] << ")   ---   ";
	///	temp_ss << setw(8) << ycontribution->type_contribution << "   " << setw(15) << name << "   i_q = " << setw(3) << i_q << "   erased: " << setw(3) << accumulate(remove_run_qTcut_TSV[i_q].begin(), remove_run_qTcut_TSV[i_q].end(), 0) << " (" << setw(3) << remove_run_qTcut_TSV[i_q].size() << ")   ---   ";
	for (int i_z = 0; i_z < n_seed; i_z++){
	  if (remove_run_qTcut_TSV[i_q][i_z]){
	    temp_ss << setprecision(3) << setw(5) << deviation_measure_qTcut[i_q][i_z] << " ";
	  }
	}
	logger << LOG_INFO << temp_ss.str() << endl;
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combine_result_original_TSV(){
  Logger logger("summary_subprocess::combine_result_original_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int x_m = 0;
  ///  N_event_TSV = accumulate(N_event_run_TSV.begin(), N_event_run_TSV.end(), 0);

  logger << LOG_DEBUG_VERBOSE << "TSV:   N_event_TSV = " << N_event_TSV << endl;
  logger << LOG_DEBUG_VERBOSE << "ycontribution->output_n_qTcut = " << ycontribution->output_n_qTcut << endl;

  //  vector<vector<double> > deviation_measure_qTcut(ycontribution->output_n_qTcut, vector<double> (n_seed, 0.));

  if (n_parallel_runs_reference_TSV[0] == 0 || Nd_event_TSV == 0.){
    //  if (n_parallel_runs_reference_TSV[0] == 0 || N_event_TSV == 0){
    for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    result_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    logger << LOG_DEBUG_VERBOSE << "TSV: no result   " << setw(20) << left << name << ":   result_TSV[" << setw(2) << right << x_m << "][" << setw(2) << right << i_q << "][" << setw(2) << right << i_s << "][" << setw(2) << right << i_r << "][" << setw(2) << right << i_f << "] = " << setw(23) << setprecision(15) << result_TSV[x_m][i_q][i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << deviation_TSV[x_m][i_q][i_s][i_r][i_f] << endl;
	  }
	}
      }
    }
  }
  else {
    for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      /*///
      if (ygeneric->deviation_tolerance_factor != 0.){
	// remove runs that indicate numerical instabilities from the combination procedure
	// (separately for each qTcut, but not for all different scales)
	double tolerance_factor = ygeneric->deviation_tolerance_factor;
	// 0 replaced by reference value:
	int x_s = osi->no_reference_TSV;
	int x_r = (osi->n_scale_ren_TSV[x_s] - 1) / 2;
	int x_f = (osi->n_scale_fact_TSV[x_s] - 1) / 2;
	vector<double> deviation_measure(n_seed, 0.);
	double min_deviation_measure = 1.e99;
	double max_deviation_measure = 0.;
	for (int i_z = 0; i_z < n_seed; i_z++){
	  if (N_event_run_TSV[i_z] > 0){
	    deviation_measure[i_z] = deviation_run_TSV[i_z][x_m][i_q][x_s][x_r][x_f] * sqrt(N_event_run_TSV[i_z]);
	  }
	}
	for (int i_z = 0; i_z < n_seed; i_z++){
	  if (N_event_run_TSV[i_z] > 0 && deviation_measure[i_z] < min_deviation_measure){min_deviation_measure = deviation_measure[i_z];}
	  if (N_event_run_TSV[i_z] > 0 && deviation_measure[i_z] > max_deviation_measure){max_deviation_measure = deviation_measure[i_z];}
	}
	N_event_TSV = 0;
	for (int i_z = 0; i_z < n_seed; i_z++){
	  deviation_measure_qTcut[i_q][i_z] = deviation_measure[i_z] / min_deviation_measure;
	  if (deviation_measure[i_z] > tolerance_factor * min_deviation_measure){remove_run_qTcut_TSV[i_q][i_z] = 1;}
	  else {N_event_TSV += N_event_run_TSV[i_z];}
	}
      }
      *////

      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    result_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    double dev = 0.;
	    double res = 0.;
	    for (int i_z = 0; i_z < n_seed; i_z++){
	      // simplify !!!
	      if (remove_run_qTcut_TSV[i_q][i_z]){
		res += 0.;
		dev += 0.;
		logger << LOG_DEBUG_VERBOSE << "removed: result_run_TSV[" << i_z << "][" << x_m << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << setw(23) << setprecision(15) << result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << " +- = " << setw(23) << setprecision(15) << deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << endl;
	      }
	      else if (N_event_run_TSV[i_z] != 0 && result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] != 0. && deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] != 0.){
		logger << LOG_DEBUG_VERBOSE << "normal:  result_run_TSV[" << i_z << "][" << x_m << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << setw(23) << setprecision(15) << result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << " +- = " << setw(23) << setprecision(15) << deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << endl;
		res += N_event_run_TSV[i_z] * result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f];
		dev += pow((N_event_run_TSV[i_z] - 1) * deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f], 2) + N_event_run_TSV[i_z] * pow(result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f], 2);
	      }
	      else if (result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0. && deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0.){
		res += 0.;
		dev += 0.;
		logger << LOG_DEBUG_VERBOSE << "empty:   result_run_TSV[" << i_z << "][" << x_m << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << setw(23) << setprecision(15) << result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << " +- = " << setw(23) << setprecision(15) << deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << endl;
	      }
	      else {logger << LOG_FATAL << "Should not happen!" << endl;}
	    }

	    if (N_valid_event_TSV[i_q] != 0 && res != 0. && dev != 0.){
	      deviation_TSV[x_m][i_q][i_s][i_r][i_f] = sqrt((N_valid_event_TSV[i_q] * dev - pow(res , 2)) / N_valid_event_TSV[i_q]) / (N_valid_event_TSV[i_q] - 1);
	      result_TSV[x_m][i_q][i_s][i_r][i_f] = res / N_valid_event_TSV[i_q];
	    }
	    else {
	      result_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	      deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    }

	    logger << LOG_DEBUG_VERBOSE << setw(20) << left << name << ":   result_TSV[" << setw(2) << right << x_m << "][" << setw(2) << right << i_q << "][" << setw(2) << right << i_s << "][" << setw(2) << right << i_r << "][" << setw(2) << right << i_f << "] = " << setw(23) << setprecision(15) << result_TSV[x_m][i_q][i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << deviation_TSV[x_m][i_q][i_s][i_r][i_f] << endl;
	  }
	}
      }
    }
  }

  /*///
  if (ygeneric->deviation_tolerance_factor != 0.){
    // information output about removed runs
    for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      if (accumulate(remove_run_qTcut_TSV[i_q].begin(), remove_run_qTcut_TSV[i_q].end(), 0)){
	stringstream temp_ss;
	temp_ss << setw(8) << ycontribution->type_contribution << "   " << setw(15) << name << "   i_q = " << setw(3) << i_q << "   erased: " << setw(3) << accumulate(remove_run_qTcut_TSV[i_q].begin(), remove_run_qTcut_TSV[i_q].end(), 0) << " (" << setw(3) << remove_run_qTcut_TSV[i_q].size() << ")   ---   ";
	for (int i_z = 0; i_z < n_seed; i_z++){
	  if (remove_run_qTcut_TSV[i_q][i_z]){
	    temp_ss << setprecision(3) << setw(5) << deviation_measure_qTcut[i_q][i_z] << " ";
	  }
	}
	logger << LOG_INFO << temp_ss.str() << endl;
      }
    }
  }
  *////

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combine_result_hybrid_TSV(){
  Logger logger("summary_subprocess::combine_result_hybrid_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int n_conservative = ycontribution->average_factor;
  int temp_n_conservative = 0;
  int temp_n_runs = 0;
  vector<int> temp_list_no_run;
  for (int i_z = 0; i_z < n_seed; i_z++){
    if (N_event_run_TSV[i_z] != 0){
      temp_n_runs++;
      temp_list_no_run.push_back(i_z);
    }
  }
  if (temp_n_runs < n_conservative){temp_n_conservative = temp_n_runs;}
  else {temp_n_conservative = n_conservative;}
  if (ycontribution->average_factor == 0){temp_n_conservative = temp_n_runs;}

  vector<vector<vector<vector<vector<vector<double> > > > > > conservative_result_run_TSV(n_seed, result_TSV);
  vector<vector<vector<vector<vector<vector<double> > > > > > conservative_deviation_run_TSV(n_seed, deviation_TSV);
  vector<long long> conservative_N_all(temp_n_conservative, 0);
  vector<long long> conservative_no_results(temp_n_conservative, 0);

  for (int i_z = 0; i_z < temp_list_no_run.size(); i_z++){
    int conservative_i_z = i_z % temp_n_conservative;
    if (!(result_run_TSV[i_z][0][0][0][0][0] == 0. && deviation_run_TSV[i_z][0][0][0][0][0] == 0.)){
      conservative_N_all[conservative_i_z] += N_event_run_TSV[i_z];
      conservative_no_results[conservative_i_z]++;
    }
  }

  int x_m = 0;
  for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  result_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	  deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
  
	  for (int i_z = 0; i_z < temp_list_no_run.size(); i_z++){
	    int conservative_i_z = i_z % temp_n_conservative;
	    if (!(result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0. && deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0.)){
	      conservative_result_run_TSV[conservative_i_z][x_m][i_q][i_s][i_r][i_f] += N_event_run_TSV[i_z] * result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f];
	      conservative_deviation_run_TSV[conservative_i_z][x_m][i_q][i_s][i_r][i_f] += pow((N_event_run_TSV[i_z] - 1) * deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f], 2) + N_event_run_TSV[i_z] * pow(result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f], 2);
      	    }
	  }
	  for (int conservative_i_z = 0; conservative_i_z < temp_n_conservative; conservative_i_z++){
	    conservative_deviation_run_TSV[conservative_i_z][x_m][i_q][i_s][i_r][i_f] = sqrt((conservative_N_all[conservative_i_z] * conservative_deviation_run_TSV[conservative_i_z][x_m][i_q][i_s][i_r][i_f] - pow(conservative_result_run_TSV[conservative_i_z][x_m][i_q][i_s][i_r][i_f], 2)) / conservative_N_all[conservative_i_z]) / (conservative_N_all[conservative_i_z] - 1); 
	    conservative_result_run_TSV[conservative_i_z][x_m][i_q][i_s][i_r][i_f] = conservative_result_run_TSV[conservative_i_z][x_m][i_q][i_s][i_r][i_f] / conservative_N_all[conservative_i_z];
	  }
	  for (int conservative_i_z = 0; conservative_i_z < temp_n_conservative; conservative_i_z++){
	    if (conservative_no_results[conservative_i_z] == 0){
	      conservative_result_run_TSV[conservative_i_z][x_m][i_q][i_s][i_r][i_f] = 0.;
	      conservative_deviation_run_TSV[conservative_i_z][x_m][i_q][i_s][i_r][i_f] = 0.;
	    }
	    logger << LOG_DEBUG_VERBOSE << "TSV[" << setw(2) << right << i_s << "] " << setw(20) << left << name << ":   cons_res_subp_TSV[" << conservative_i_z << "][" << x_m << "][" << setw(2) << right << i_q << "][" << setw(2) << right << i_s << "][" << setw(2) << right << i_r << "][" << setw(2) << right << i_f << "] = " << setw(13) << right << setprecision(6) << conservative_result_run_TSV[conservative_i_z][x_m][i_q][i_s][i_r][i_f] << " +- " << setw(13) << setprecision(6) << conservative_deviation_run_TSV[conservative_i_z][x_m][i_q][i_s][i_r][i_f] << endl;
	  }
	}
      }
    }
  }

  N_event_TSV = 0;
  for (int conservative_i_z = 0; conservative_i_z < temp_n_conservative; conservative_i_z++){
    N_event_TSV += conservative_N_all[conservative_i_z];
  }

  logger << LOG_DEBUG_VERBOSE << "TSV:   N_event_TSV = " << N_event_TSV << endl;
  if (N_event_TSV != 0){
    for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    result_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    for (int i_z = 0; i_z < temp_n_conservative; i_z++){
	      if (N_event_run_TSV[i_z] != 0 && conservative_result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] != 0. && conservative_deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] != 0.){
		double temp_singleweight = 1. / pow(conservative_deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f], 2);
		logger << LOG_DEBUG_VERBOSE << "z_result_TSV[" << i_z << "][" << x_m << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << setw(23) << setprecision(15) << conservative_result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << " +- = " << setw(23) << setprecision(15) << conservative_deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << endl;
		result_TSV[x_m][i_q][i_s][i_r][i_f] += temp_singleweight * conservative_result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f];
		deviation_TSV[x_m][i_q][i_s][i_r][i_f] += temp_singleweight;
	      }
	      else if (N_event_run_TSV[i_z] != 0 && conservative_result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0. && conservative_deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0.){
		result_TSV[x_m][i_q][i_s][i_r][i_f] += 0.;
		deviation_TSV[x_m][i_q][i_s][i_r][i_f] += 0.;
		logger << LOG_DEBUG_VERBOSE << "empty: z_result_TSV[" << i_z << "][" << x_m << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << setw(23) << setprecision(15) << conservative_result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << " +- = " << setw(23) << setprecision(15) << conservative_deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << endl;
	      }
	    }
	    if (N_event_TSV != 0 && result_TSV[x_m][i_q][i_s][i_r][i_f] != 0. && deviation_TSV[x_m][i_q][i_s][i_r][i_f] != 0.){
	      result_TSV[x_m][i_q][i_s][i_r][i_f] = result_TSV[x_m][i_q][i_s][i_r][i_f] / deviation_TSV[x_m][i_q][i_s][i_r][i_f];
	      deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 1. / sqrt(deviation_TSV[x_m][i_q][i_s][i_r][i_f]);
	    }
	    else {
	      result_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	      deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    }
	    logger << LOG_DEBUG_VERBOSE << "TSV:   " << setw(20) << left << name << ":   x_result_TSV[" << setw(2) << right << x_m << "][" << setw(2) << right << i_q << "][" << setw(2) << right << i_s << "][" << setw(2) << right << i_r << "][" << setw(2) << right << i_f << "] = " << setw(23) << setprecision(15) << result_TSV[x_m][i_q][i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << deviation_TSV[x_m][i_q][i_s][i_r][i_f] << endl;
	  }
	}
      }
    }
  }
  else if (N_event_TSV == 0){
    for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    result_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    logger << LOG_DEBUG_VERBOSE << "TSV: no result   " << setw(20) << left << name << ":   x_result_TSV[" << setw(2) << right << x_m << "][" << setw(2) << right << i_q << "][" << setw(2) << right << i_s << "][" << setw(2) << right << i_r << "][" << setw(2) << right << i_f << "] = " << setw(23) << setprecision(15) << result_TSV[x_m][i_q][i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << deviation_TSV[x_m][i_q][i_s][i_r][i_f] << endl;
	  }
	}
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combine_result_alternative_TSV(){
  Logger logger("summary_subprocess::combine_result_alternative_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int x_m = 0;
  N_event_TSV = accumulate(N_event_run_TSV.begin(), N_event_run_TSV.end(), 0);
  logger << LOG_DEBUG_VERBOSE << "TSV:   N_event_TSV = " << N_event_TSV << endl;
  if (N_event_TSV != 0){
     for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    result_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    for (int i_z = 0; i_z < n_seed; i_z++){
	      if (N_event_run_TSV[i_z] != 0 && result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] != 0. && deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] != 0.){
		double temp_singleweight = 1. / pow(deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f], 2);
		logger << LOG_DEBUG_VERBOSE << "z_result_TSV[" << i_z << "][" << x_m << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << setw(23) << setprecision(15) << result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << " +- = " << setw(23) << setprecision(15) << deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << endl;
		result_TSV[x_m][i_q][i_s][i_r][i_f] += temp_singleweight * result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f];
		deviation_TSV[x_m][i_q][i_s][i_r][i_f] += temp_singleweight;
	      }
	      else if (N_event_run_TSV[i_z] != 0 && result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0. && deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0.){
		result_TSV[x_m][i_q][i_s][i_r][i_f] += 0.;
		deviation_TSV[x_m][i_q][i_s][i_r][i_f] += 0.;
		logger << LOG_DEBUG_VERBOSE << "empty: z_result_TSV[" << i_z << "][" << x_m << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << setw(23) << setprecision(15) << result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << " +- = " << setw(23) << setprecision(15) << deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] << endl;
	      }
	    }
	    if (N_event_TSV != 0 && result_TSV[x_m][i_q][i_s][i_r][i_f] != 0. && deviation_TSV[x_m][i_q][i_s][i_r][i_f] != 0.){
	      result_TSV[x_m][i_q][i_s][i_r][i_f] = result_TSV[x_m][i_q][i_s][i_r][i_f] / deviation_TSV[x_m][i_q][i_s][i_r][i_f];
	      deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 1. / sqrt(deviation_TSV[x_m][i_q][i_s][i_r][i_f]);
	    }
	    else {
	      result_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	      deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    }
	    logger << LOG_DEBUG_VERBOSE << "TSV:   " << setw(20) << left << name << ":   x_result_TSV[" << setw(2) << right << x_m << "][" << setw(2) << right << i_q << "][" << setw(2) << right << i_s << "][" << setw(2) << right << i_r << "][" << setw(2) << right << i_f << "] = " << setw(23) << setprecision(15) << result_TSV[x_m][i_q][i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << deviation_TSV[x_m][i_q][i_s][i_r][i_f] << endl;
	  }
	}
      }
    }
  }
  else if (N_event_TSV == 0){
    for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    result_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    deviation_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    logger << LOG_DEBUG_VERBOSE << "TSV: no result   " << setw(20) << left << name << ":   x_result_TSV[" << setw(2) << right << x_m << "][" << setw(2) << right << i_q << "][" << setw(2) << right << i_s << "][" << setw(2) << right << i_r << "][" << setw(2) << right << i_f << "] = " << setw(23) << setprecision(15) << result_TSV[x_m][i_q][i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << deviation_TSV[x_m][i_q][i_s][i_r][i_f] << endl;
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void summary_subprocess::readin_result_CV(string & result_moment){
  Logger logger("summary_subprocess::readin_result_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_z = 0; i_z < n_seed; i_z++){
    string result_file = result_moment + "_" + name + ".dat";
    logger << LOG_DEBUG_VERBOSE << "result_file = " << result_file << endl;
    string filename = "../" + ycontribution->directory[i_z] + "/result/" + result_file;
    logger << LOG_DEBUG << "Input:   filename = " << filename << endl;
    vector<string> readin;
    char LineBuffer[128];
    ifstream in_result(filename.c_str());  
    while (in_result.getline(LineBuffer, 128)){readin.push_back(LineBuffer);}
    in_result.close();
    logger << LOG_DEBUG_VERBOSE << "default_size_CV  = " << default_size_CV << endl;
    logger << LOG_DEBUG_VERBOSE << "readin.size() = " << readin.size() << endl;
    if (readin.size() > default_size_CV){
      /*
      logger << LOG_DEBUG_VERBOSE << "error2_time_run.size() = " << error2_time_run.size() << endl;
      logger << LOG_DEBUG_VERBOSE << "time_per_event_run.size() = " << time_per_event_run.size() << endl;
      logger << LOG_DEBUG_VERBOSE << "error2_time.size() = " << error2_time.size() << endl;
      logger << LOG_DEBUG_VERBOSE << "time_per_event.size() = " << time_per_event.size() << endl;
      */
      /*
      error2_time_run[i_z] = atof(readin[default_size_CV + 2].c_str());
      time_per_event_run[i_z] = atof(readin[default_size_CV + 1].c_str());
      used_runtime_run[i_z] = atof(readin[default_size_CV].c_str());
      */
      /* 
      logger << LOG_DEBUG_VERBOSE << "error2_time_run[" << i_p << "][" << i_z << "] = " << error2_time_run[i_z] << endl;
      logger << LOG_DEBUG_VERBOSE << "time_per_event_run[" << i_p << "][" << i_z << "] = " << time_per_event_run[i_z] << endl;
      logger << LOG_DEBUG_VERBOSE << "used_runtime_run[" << i_p << "][" << i_z << "] = " << used_runtime_run[i_z] << endl;
      */  
      readin.erase(readin.begin() + default_size_CV, readin.end());
    }
      
    if (readin.size() == default_size_CV){
      used_n_event_run[i_z] = atoll(readin[0].c_str());
      N_event_run_CV[i_z] = atoll(readin[0].c_str());

      result_run[i_z] = atof(readin[1].c_str());
      deviation_run[i_z] = atof(readin[2].c_str());
      int proc = 3;

      for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  result_run_CV[i_q][i_s][i_z] = atof(readin[proc + 2 * (i_s + osi->n_scales_CV * i_q)].c_str());
	  deviation_run_CV[i_q][i_s][i_z] = atof(readin[proc + 1 + 2 * (i_s + osi->n_scales_CV * i_q)].c_str());
	  //	  logger << LOG_INFO << setw(20) << name << setw(15) << ycontribution->type_contribution << "   result_run_CV[" << i_q << "][" << i_s << "][" << i_z << "] = " << result_run_CV[i_q][i_s][i_z] << " +/- " << deviation_run_CV[i_q][i_s][i_z] << endl;
	}
      }
    }
    else {
      if (readin.size() != 0){logger << LOG_ERROR << result_file << "   error" << endl;}
      else {logger << LOG_DEBUG << result_file << "   empty" << endl;}
      used_n_event_run[i_z] = 0;
      N_event_run_CV[i_z] = 0;
      result_run[i_z] = 0.;
      deviation_run[i_z] = 0.;
      for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  result_run_CV[i_q][i_s][i_z] = 0.;
	  deviation_run_CV[i_q][i_s][i_z] = 0.;
	}
      }
    }
  }

  int x_q = 0;
  int x_s = (osi->n_scales_CV - 1) / 2;

  for (int i_z = 0; i_z < n_seed; i_z++){
    if (N_event_run_CV[i_z] > 0){
      logger << LOG_INFO << setw(20) << name << setw(15) << ycontribution->type_contribution << "   result_run_CV[" << x_q << "][" << x_s << "][" << i_z << "] = " << result_run_CV[x_q][x_s][i_z] << " +/- " << deviation_run_CV[x_q][x_s][i_z] << endl;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combination_result_CV(string & result_moment){
  Logger logger("summary_subprocess::combination_result_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG << left << name << endl;

  int x_s = (osi->n_scales_CV - 1) / 2;

  for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
    n_parallel_runs_reference_CV[i_q] = 0;
    for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
      logger << LOG_DEBUG_VERBOSE << "N_event_run_CV[" << setw(4) << i_z << "] = " << N_event_run_CV[i_z] << endl;
      if (N_event_run_CV[i_z] != 0 && 
	  result_run_CV[i_q][x_s][i_z] != 0. && 
	  deviation_run_CV[i_q][x_s][i_z] != 0.){
	n_parallel_runs_reference_CV[i_q]++;
	if (!remove_run_qTcut_CV[i_q][i_z]){
	  n_valid_parallel_runs_reference_CV[i_q]++;
	}
      }
    }
  }


  for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
    ///    logger << LOG_INFO << "n_parallel_runs_reference_CV[" << i_q << "] = " << setw(4) << n_parallel_runs_reference_CV[i_q] << "   n_valid_parallel_runs_reference_CV[" << i_q << "] = " << setw(4) << n_valid_parallel_runs_reference_CV[i_q] << endl;
  }

  ///  logger << LOG_INFO << "ygeneric->average_factor = " << ygeneric->average_factor << endl;

  removal_exceptional_runs_result_CV();

  if (ygeneric->average_factor == 0){combine_result_alternative_CV();}  // bug: results zero !!!
  else if (ygeneric->average_factor == 1){combine_result_original_CV();}  // removal of exceptional runs included
  else {combine_result_hybrid_CV();}  // removal of exceptional runs not included !!!

  int x_q = 0;
  //  int x_s = (osi->n_scales_CV - 1) / 2;
  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
    logger << LOG_INFO << "result_CV[" << x_q << "][" << i_s << "] = " << result_CV[x_q][i_s] << endl;
  }

  // higher qTcut-values are identified with i_q = 0 (for qTcut independent runs); could be improved...
  if (ycontribution->output_n_qTcut < osi->n_qTcut){
    for (int i_q = ycontribution->output_n_qTcut; i_q < osi->n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	result_CV[i_q][i_s] = result_CV[0][i_s];
	deviation_CV[i_q][i_s] = deviation_CV[0][i_s];
      }
    }
  }


  for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
	if (N_event_run_CV[i_z] != 0 && 
	    result_run_CV[i_q][i_s][i_z] != 0. && 
	    deviation_run_CV[i_q][i_s][i_z] != 0.){
	  n_parallel_runs_CV[i_q][i_s]++;
	  if (!remove_run_qTcut_CV[i_q][i_z]){
	    n_valid_parallel_runs_CV[i_q][i_s]++;
	  }
	}
      }
    
      chi2_CV[i_q][i_s] = 0.;
      if (n_parallel_runs_CV[i_q][i_s] > 1){
	for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
	  if (!(result_run_CV[i_q][i_s][i_z] == 0. && deviation_run_CV[i_q][i_s][i_z] == 0.)){
	    chi2_CV[i_q][i_s] += pow((result_run_CV[i_q][i_s][i_z] - result_CV[i_q][i_s]) / deviation_run_CV[i_q][i_s][i_z], 2.);
	  }
	}
	chi2_CV[i_q][i_s] = chi2_CV[i_q][i_s] / n_parallel_runs_CV[i_q][i_s];
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::removal_exceptional_runs_result_CV(){
  Logger logger("summary_subprocess::removal_exceptional_runs_result_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  N_event_CV = accumulate(N_event_run_CV.begin(), N_event_run_CV.end(), 0);
  //  Nd_event_CV = accumulate(N_event_run_CV.begin(), N_event_run_CV.end(), 0);
  Nd_event_CV = 0.;
  for (int i_z = 0; i_z < n_seed; i_z++){
    Nd_event_CV += double(N_event_run_CV[i_z]);
  }

  vector<vector<double> > deviation_measure_qTcut(osi->n_qTcut, vector<double> (n_seed, 0.));

  ///  logger << LOG_INFO << left << setw(10) << ycontribution->type_contribution << setw(20) << name << "N_event_CV = " << N_event_CV << "   Nd_event_CV = " << Nd_event_CV << "   n_parallel_runs_reference_CV[0] = " << setw(4) << n_parallel_runs_reference_CV[0] << "   n_valid_parallel_runs_reference_CV[0] = " << setw(4) << n_valid_parallel_runs_reference_CV[0] << endl;

  //  if (n_parallel_runs_reference_CV[0] > 0 && N_event_CV > 0){
  if (n_parallel_runs_reference_CV[0] > 0 && Nd_event_CV > 0.){
    for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      if (ygeneric->deviation_tolerance_factor != 0.){
	// remove runs that indicate numerical instabilities from the combination procedure
	// (separately for each qTcut, but not for all different scales)
	double tolerance_factor = ygeneric->deviation_tolerance_factor;
	int i_s = 0;
	if (osi->switch_CV == 1 || osi->switch_CV == 2 || osi->switch_CV == 3 || osi->switch_CV == 4){i_s = (osi->n_scales_CV - 1) / 2;}
	else if (osi->switch_CV == 5){i_s = 3;}
	else if (osi->switch_CV == 6){i_s = 4;}
	else {logger << LOG_FATAL << "Illegal value of switch_CV = " << osi->switch_CV << "."; exit(1);}
	vector<double> deviation_measure(n_seed, 0.);
	double min_deviation_measure = 1.e99;
	double max_deviation_measure = 0.;
	for (int i_z = 0; i_z < n_seed; i_z++){
	  if (N_event_run_CV[i_z] > 0){
	    deviation_measure[i_z] = deviation_run_CV[i_q][i_s][i_z] * sqrt(N_event_run_CV[i_z]);
	  }
	}
	for (int i_z = 0; i_z < n_seed; i_z++){
	  if (N_event_run_CV[i_z] > 0 && deviation_measure[i_z] < min_deviation_measure){min_deviation_measure = deviation_measure[i_z];}
	  if (N_event_run_CV[i_z] > 0 && deviation_measure[i_z] > max_deviation_measure){max_deviation_measure = deviation_measure[i_z];}
	}
	N_valid_event_CV[i_q] = 0;
	for (int i_z = 0; i_z < n_seed; i_z++){
	  deviation_measure_qTcut[i_q][i_z] = deviation_measure[i_z] / min_deviation_measure;
	  if (deviation_measure[i_z] > tolerance_factor * min_deviation_measure){remove_run_qTcut_CV[i_q][i_z] = 1;}
	  else {N_valid_event_CV[i_q] += N_event_run_CV[i_z];}
	  ///	  logger << LOG_INFO << left << setw(10) << ycontribution->type_contribution << setw(20) << name << "   N_valid_event_CV[" << setw(3) << i_q << "] = " << N_valid_event_CV[i_q] << "   N_event_run_CV[" << i_z << "] = " << N_event_run_CV[i_z] << endl;
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combine_result_original_CV(){
  Logger logger("summary_subprocess::combine_result_original_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ///  N_event_CV = accumulate(N_event_run_CV.begin(), N_event_run_CV.end(), 0);

  //  vector<vector<double> > deviation_measure_qTcut(osi->n_qTcut, vector<double> (n_seed, 0.));

  ///  logger << LOG_INFO << "n_parallel_runs_reference_CV[0] = " << n_parallel_runs_reference_CV[0] << "   N_event_CV = " << N_event_CV << endl;

  if (n_parallel_runs_reference_CV[0] == 0 || Nd_event_CV == 0.){
    //  if (n_parallel_runs_reference_CV[0] == 0 || N_event_CV == 0){
    for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	result_CV[i_q][i_s] = 0.;
	deviation_CV[i_q][i_s] = 0.;
      }
    }
  }
  else {
    for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
      /****
      if (ygeneric->deviation_tolerance_factor != 0.){
	// remove runs that indicate numerical instabilities from the combination procedure
	// (separately for each qTcut, but not for all different scales)
	double tolerance_factor = ygeneric->deviation_tolerance_factor;
	int i_s = 0;
	if (osi->switch_CV == 1 || osi->switch_CV == 2 || osi->switch_CV == 3 || osi->switch_CV == 4){i_s = (osi->n_scales_CV - 1) / 2;}
	else if (osi->switch_CV == 5){i_s = 3;}
	else if (osi->switch_CV == 6){i_s = 4;}
	else {logger << LOG_FATAL << "Illegal value of switch_CV = " << osi->switch_CV << "."; exit(1);}
	vector<double> deviation_measure(n_seed, 0.);
	double min_deviation_measure = 1.e99;
	double max_deviation_measure = 0.;
	for (int i_z = 0; i_z < n_seed; i_z++){
	  if (N_event_run_CV[i_z] > 0){
	    deviation_measure[i_z] = deviation_run_CV[i_q][i_s][i_z] * sqrt(N_event_run_CV[i_z]);
	  }
	}
	for (int i_z = 0; i_z < n_seed; i_z++){
	  if (N_event_run_CV[i_z] > 0 && deviation_measure[i_z] < min_deviation_measure){min_deviation_measure = deviation_measure[i_z];}
	  if (N_event_run_CV[i_z] > 0 && deviation_measure[i_z] > max_deviation_measure){max_deviation_measure = deviation_measure[i_z];}
	}
	N_event_CV = 0;
	for (int i_z = 0; i_z < n_seed; i_z++){
	  deviation_measure_qTcut[i_q][i_z] = deviation_measure[i_z] / min_deviation_measure;
	  if (deviation_measure[i_z] > tolerance_factor * min_deviation_measure){remove_run_qTcut_CV[i_q][i_z] = 1;}
	  else {N_event_CV += N_event_run_CV[i_z];}
	}
      }
      *////

      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	result_CV[i_q][i_s] = 0.;
	deviation_CV[i_q][i_s] = 0.;

	double dev = 0.;
	double res = 0.;
	for (int i_z = 0; i_z < n_seed; i_z++){
	  /*
	  if (i_q == 0 && i_s == 3){
	    logger << LOG_INFO << "remove_run_qTcut_CV[" << i_q << "][" << i_z << "] = " << remove_run_qTcut_CV[i_q][i_z] 
		   << "   N_event_run_CV[" << i_z << "] = " << N_event_run_CV[i_z]
		   << "   result_run_CV[" << i_q << "][" << i_s << "][" << i_z << "] = " << result_run_CV[i_q][i_s][i_z]
		   << "   deviation_run_CV[" << i_q << "][" << i_s << "][" << i_z << "] = " << deviation_run_CV[i_q][i_s][i_z] << endl;

	  }
	  */

	  if (remove_run_qTcut_CV[i_q][i_z]){
	    res += 0.;
	    dev += 0.;
	    logger << LOG_DEBUG_VERBOSE << "removed: result_run_CV[" << i_s << "][" << i_q << "][" << i_z << "] = " << setw(23) << setprecision(15) << result_run_CV[i_q][i_s][i_z] << " +- = " << setw(23) << setprecision(15) << deviation_run_CV[i_q][i_s][i_z] << endl;
	  }
	  else if (N_event_run_CV[i_z] != 0 && result_run_CV[i_q][i_s][i_z] != 0. && deviation_run_CV[i_q][i_s][i_z] != 0.){
	    logger << LOG_DEBUG_VERBOSE << "normal:  result_run_CV[" << i_s << "][" << i_q << "][" << i_z << "] = " << setw(23) << setprecision(15) << result_run_CV[i_q][i_s][i_z] << " +- = " << setw(23) << setprecision(15) << deviation_run_CV[i_q][i_s][i_z] << endl;

	    res += N_event_run_CV[i_z] * result_run_CV[i_q][i_s][i_z];
	    dev += pow((N_event_run_CV[i_z] - 1) * deviation_run_CV[i_q][i_s][i_z], 2) + N_event_run_CV[i_z] * pow(result_run_CV[i_q][i_s][i_z], 2);
	  }
	  else if (result_run_CV[i_q][i_s][i_z] == 0. && deviation_run_CV[i_q][i_s][i_z] == 0.){
	    res += 0.;
	    dev += 0.;
	  }
	  else {
	    logger << LOG_FATAL << name << endl;
	    logger << LOG_FATAL << "N[" << i_z << "] = " << N_event_run_CV[i_z] << endl;
	    logger << LOG_FATAL << "result_run_CV[" << i_q << "][" << i_s << "][" << i_z << "] = " << result_run_CV[i_q][i_s][i_z] << endl;
	    logger << LOG_FATAL << "deviation_run_CV[" << i_q << "][" << i_s << "][" << i_z << "] = " << deviation_run_CV[i_q][i_s][i_z] << endl;
	    logger << LOG_FATAL << "Should not happen!" << endl;
	  }
	  ///	  if (i_q == 0){logger << LOG_INFO << "res = " << res << "   i_z = " << i_z << endl;}
	}
	
	if (i_q == 0){
	  ///	  logger << LOG_INFO << setw(20) << name << "   N_valid_event_CV[" << i_q << "] = " << N_valid_event_CV[i_q] << "   res = " << res << "   dev = " << dev << endl;
	}

	if (N_valid_event_CV[i_q] != 0 && res != 0. && dev != 0.){
	  deviation_CV[i_q][i_s] = sqrt((N_valid_event_CV[i_q] * dev - pow(res, 2)) / N_valid_event_CV[i_q]) / (N_valid_event_CV[i_q] - 1); 
	  result_CV[i_q][i_s] = res / N_valid_event_CV[i_q];
	}
	else {
	  result_CV[i_q][i_s] = 0.;
	  deviation_CV[i_q][i_s] = 0.;
	}
      }
    }
  }

  if (ycontribution->output_n_qTcut < osi->n_qTcut){
    for (int i_q = ycontribution->output_n_qTcut; i_q < osi->n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	result_CV[i_q][i_s] = result_CV[0][i_s];
	deviation_CV[i_q][i_s] = deviation_CV[0][i_s];
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void summary_subprocess::combine_result_alternative_CV(){
  Logger logger("summary_subprocess::combine_result_alternative_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double temp_singleweight;
  //  for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
  for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      result_CV[i_q][i_s] = 0.;
      deviation_CV[i_q][i_s] = 0.;
      for (int i_z = 0; i_z < n_seed; i_z++){
	if (N_event_run_CV[i_z] != 0 && deviation_run_CV[i_q][i_s][i_z] != 0. && result_run_CV[i_q][i_s][i_z] != 0.){
	  temp_singleweight = 1. / pow(deviation_run_CV[i_q][i_s][i_z], 2);
	  result_CV[i_q][i_s] += temp_singleweight * result_run_CV[i_q][i_s][i_z];
	  deviation_CV[i_q][i_s] += temp_singleweight;
	}
      }
      if (N_event_CV != 0){
	result_CV[i_q][i_s] = result_CV[i_q][i_s] / deviation_CV[i_q][i_s];
	deviation_CV[i_q][i_s] = 1. / sqrt(deviation_CV[i_q][i_s]);
	if (result_CV[i_q][i_s] != result_CV[i_q][i_s]){
	  logger << LOG_DEBUG_VERBOSE << N_event_CV << endl;
	  for (int i_z = 0; i_z < n_seed; i_z++){logger << LOG_DEBUG_VERBOSE << N_event_run_CV[i_z] << endl;}
	}
      }
      else {
	result_CV[i_q][i_s] = 0.;
	deviation_CV[i_q][i_s] = 0.;
      }

      logger << LOG_DEBUG << "CV[" << setw(2) << right << i_s << "] " << setw(20) << left << name << ":   result_subprocess_CV[" << setw(2) << right << i_q << "][" << setw(2) << right << i_s << "] = " << setw(13) << right << setprecision(6) << result_CV[i_q][i_s] << " +- " << setw(13) << setprecision(6) << deviation_CV[i_q][i_s] << endl;
    }
  }
  
  if (ycontribution->output_n_qTcut < osi->n_qTcut){
    for (int i_q = ycontribution->output_n_qTcut; i_q < osi->n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	result_CV[i_q][i_s] = result_CV[0][i_s];
      }
    }
  }

 logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combine_result_hybrid_CV(){
  Logger logger("summary_subprocess::combine_result_hybrid_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int n_conservative = ycontribution->average_factor;
  int temp_n_conservative = 0;

  int temp_n_runs = 0;
  vector<int> temp_list_no_run;
  for (int i_z = 0; i_z < n_seed; i_z++){
    if (N_event_run_CV[i_z] != 0){
      temp_n_runs++;
      temp_list_no_run.push_back(i_z);
    }
  }

  if (temp_n_runs < n_conservative){
    temp_n_conservative = temp_n_runs;
  }
  else {
    temp_n_conservative = n_conservative;
  }
  vector<vector<vector<double> > > conservative_result_run_CV(ycontribution->output_n_qTcut, vector<vector<double> > (osi->n_scales_CV, vector<double> (temp_n_conservative, 0.)));
  vector<vector<vector<double> > > conservative_deviation_run_CV(ycontribution->output_n_qTcut, vector<vector<double> > (osi->n_scales_CV, vector<double> (temp_n_conservative, 0.)));
  vector<vector<vector<double> > > conservative_chi2_CV(ycontribution->output_n_qTcut, vector<vector<double> > (osi->n_scales_CV, vector<double> (temp_n_conservative, 0.)));
  vector<long long> conservative_N_all(temp_n_conservative, 0);
  vector<long long> conservative_no_results(temp_n_conservative, 0);

  for (int i_z = 0; i_z < temp_list_no_run.size(); i_z++){
    int conservative_i_z = i_z % temp_n_conservative;
    if (!(result_run_CV[0][0][i_z] == 0. && deviation_run_CV[0][0][i_z] == 0.)){
      conservative_N_all[conservative_i_z] += N_event_run_CV[i_z];
      conservative_no_results[conservative_i_z]++;
    }
  }

  for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      result_CV[i_q][i_s] = 0.;
      deviation_CV[i_q][i_s] = 0.;
      for (int i_z = 0; i_z < temp_list_no_run.size(); i_z++){
	int conservative_i_z = i_z % temp_n_conservative;
	if (!(result_run_CV[i_q][i_s][i_z] == 0. && deviation_run_CV[i_q][i_s][i_z] == 0.)){
	  conservative_result_run_CV[i_q][i_s][conservative_i_z] += N_event_run_CV[i_z] * result_run_CV[i_q][i_s][i_z];
	  conservative_deviation_run_CV[i_q][i_s][conservative_i_z] += pow((N_event_run_CV[i_z] - 1) * deviation_run_CV[i_q][i_s][i_z], 2) + N_event_run_CV[i_z] * pow(result_run_CV[i_q][i_s][i_z], 2);
	}
      }

      for (int conservative_i_z = 0; conservative_i_z < temp_n_conservative; conservative_i_z++){
	conservative_deviation_run_CV[i_q][i_s][conservative_i_z] = sqrt((conservative_N_all[conservative_i_z] * conservative_deviation_run_CV[i_q][i_s][conservative_i_z] - pow(conservative_result_run_CV[i_q][i_s][conservative_i_z], 2)) / conservative_N_all[conservative_i_z]) / (conservative_N_all[conservative_i_z] - 1); 
	conservative_result_run_CV[i_q][i_s][conservative_i_z] = conservative_result_run_CV[i_q][i_s][conservative_i_z] / conservative_N_all[conservative_i_z];
	conservative_chi2_CV[i_q][i_s][conservative_i_z] = 0.;
      }

      for (int i_z = 0; i_z < temp_list_no_run.size(); i_z++){
	int conservative_i_z = i_z % temp_n_conservative;
	if (!(result_run_CV[i_q][i_s][i_z] == 0. && deviation_run_CV[i_q][i_s][i_z] == 0.)){
	  conservative_chi2_CV[i_q][i_s][conservative_i_z] += pow((result_run_CV[i_q][i_s][i_z] - conservative_result_run_CV[i_q][i_s][conservative_i_z]) / deviation_run_CV[i_q][i_s][i_z], 2.);
	}
      }

      for (int conservative_i_z = 0; conservative_i_z < temp_n_conservative; conservative_i_z++){
	conservative_chi2_CV[i_q][i_s][conservative_i_z] = conservative_chi2_CV[i_q][i_s][conservative_i_z] / conservative_no_results[conservative_i_z];
	if (conservative_no_results[conservative_i_z] == 0){
	  conservative_result_run_CV[i_q][i_s][conservative_i_z] = 0.;
	  conservative_deviation_run_CV[i_q][i_s][conservative_i_z] = 0.;
	  conservative_chi2_CV[i_q][i_s][conservative_i_z] = 0.;
	}
	logger << LOG_DEBUG << "CV[" << setw(2) << right << i_s << "] " << setw(20) << left << name << ":   cons_res_subp_CV[" << setw(2) << right << i_q << "][" << setw(2) << right << i_s << "][" << conservative_i_z << "] = " << setw(13) << right << setprecision(6) << conservative_result_run_CV[i_q][i_s][conservative_i_z] << " +- " << setw(13) << setprecision(6) << conservative_deviation_run_CV[i_q][i_s][conservative_i_z] << "   chi2 = " << setw(6) << setprecision(4) << conservative_chi2_CV[i_q][i_s][conservative_i_z] << endl;

      }
    }
  }
  


 

  N_event_CV = 0;
  for (int conservative_i_z = 0; conservative_i_z < temp_n_conservative; conservative_i_z++){
    N_event_CV += conservative_N_all[conservative_i_z];
  }

  double temp_singleweight = 0.;
  for (int i_q = 0; i_q < ycontribution->output_n_qTcut; i_q++){
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      result_CV[i_q][i_s] = 0.;
      deviation_CV[i_q][i_s] = 0.;
      for (int i_z = 0; i_z < temp_n_conservative; i_z++){
	if (conservative_N_all[i_z] != 0 && conservative_deviation_run_CV[i_q][i_s][i_z] != 0. && conservative_result_run_CV[i_q][i_s][i_z] != 0.){
	  temp_singleweight = 1. / pow(conservative_deviation_run_CV[i_q][i_s][i_z], 2);
	  result_CV[i_q][i_s] += temp_singleweight * conservative_result_run_CV[i_q][i_s][i_z];
	  deviation_CV[i_q][i_s] += temp_singleweight;
	}
      }
      if (N_event_CV != 0){
	result_CV[i_q][i_s] = result_CV[i_q][i_s] / deviation_CV[i_q][i_s];
	deviation_CV[i_q][i_s] = 1. / sqrt(deviation_CV[i_q][i_s]);
	if (result_CV[i_q][i_s] != result_CV[i_q][i_s]){
	  logger << LOG_DEBUG_VERBOSE << N_event_CV << endl;
	  for (int i_z = 0; i_z < temp_n_conservative; i_z++){logger << LOG_DEBUG_VERBOSE << conservative_N_all[i_z] << endl;}
	}
      }
      else {
	result_CV[i_q][i_s] = 0.;
	deviation_CV[i_q][i_s] = 0.;
      }
    }
  }
  
  if (ycontribution->output_n_qTcut < osi->n_qTcut){
    for (int i_q = ycontribution->output_n_qTcut; i_q < osi->n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	result_CV[i_q][i_s] = result_CV[0][i_s];
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


