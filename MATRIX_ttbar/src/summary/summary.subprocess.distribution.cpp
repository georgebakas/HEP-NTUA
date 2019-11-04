#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_subprocess::initialization_distribution_TSV(){
  Logger logger("summary_subprocess::initialization_distribution_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  pid_t pid = getpid();
  stringstream temp_pid;
  temp_pid << "cat /proc/" << pid << "/status | grep VmSize";
  string temp_s_pid = temp_pid.str();
  logger << LOG_INFO << "begin initialization_distribution_TSV:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  int temp_i = system(temp_s_pid.c_str());


  distribution_result_TSV.resize(osi->extended_distribution.size());
  distribution_deviation_TSV.resize(osi->extended_distribution.size());
  distribution_chi2_TSV.resize(osi->extended_distribution.size());

  distribution_N_total_TSV.resize(osi->extended_distribution.size());
  distribution_N_total_binwise_TSV.resize(osi->extended_distribution.size());

  logger << LOG_DEBUG << ycontribution->type_contribution << "." << ycontribution->type_correction << "   " << name << "   ycontribution->selection_n_qTcut = " << ycontribution->selection_n_qTcut << endl;

  //  vector<vector<vector<vector<vector<vector<long long> > > > > > temp_vvvvvvll_dbqsrf(osi->extended_distribution.size());
  vector<vector<long long> > temp_vvll_db(osi->extended_distribution.size());
  vector<vector<vector<long long> > > temp_vvvll_dbq(osi->extended_distribution.size());
  vector<vector<vector<int> > > temp_vvvi_dbq(osi->extended_distribution.size());
  
  logger << LOG_INFO << "resize step: extended_distribution started:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());

  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    temp_vvll_db[i_d].resize(osi->extended_distribution[i_d].n_bins);
    temp_vvvll_dbq[i_d].resize(osi->extended_distribution[i_d].n_bins);
    temp_vvvi_dbq[i_d].resize(osi->extended_distribution[i_d].n_bins);
    //    temp_vvvvvvll_dbqsrf[i_d].resize(osi->extended_distribution[i_d].n_bins);
    
    distribution_result_TSV[i_d].resize(osi->extended_distribution[i_d].n_bins);
    distribution_deviation_TSV[i_d].resize(osi->extended_distribution[i_d].n_bins);
    distribution_chi2_TSV[i_d].resize(osi->extended_distribution[i_d].n_bins);

    // total is also binwise, but after removal of runs...
    distribution_N_total_TSV[i_d].resize(osi->extended_distribution[i_d].n_bins, vector<long long> (ycontribution->selection_n_qTcut, 0));
    distribution_N_total_binwise_TSV[i_d].resize(osi->extended_distribution[i_d].n_bins, vector<long long> (ycontribution->selection_n_qTcut, 0));
    
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
      temp_vvvll_dbq[i_d][i_b].resize(ycontribution->selection_n_qTcut, 0);
      temp_vvvi_dbq[i_d][i_b].resize(ycontribution->selection_n_qTcut, 0);
      //      temp_vvvvvvll_dbqsrf[i_d][i_b].resize(ycontribution->selection_n_qTcut);
	  
      distribution_result_TSV[i_d][i_b].resize(ycontribution->selection_n_qTcut);
      distribution_deviation_TSV[i_d][i_b].resize(ycontribution->selection_n_qTcut);
      distribution_chi2_TSV[i_d][i_b].resize(ycontribution->selection_n_qTcut);
      for (int x_q = 0; x_q < ycontribution->selection_n_qTcut; x_q++){
	//	temp_vvvvvvll_dbqsrf[i_d][i_b][x_q].resize(osi->n_extended_set_TSV);

	distribution_result_TSV[i_d][i_b][x_q].resize(osi->n_extended_set_TSV);
	distribution_deviation_TSV[i_d][i_b][x_q].resize(osi->n_extended_set_TSV);
	distribution_chi2_TSV[i_d][i_b][x_q].resize(osi->n_extended_set_TSV);
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  if (!osi->switch_distribution_TSV[i_s]){continue;}
	  //	  temp_vvvvvvll_dbqsrf[i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<long long> (osi->n_scale_fact_TSV[i_s], 0));

	  distribution_result_TSV[i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  distribution_deviation_TSV[i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  distribution_chi2_TSV[i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	}
      }
    }
  }

  logger << LOG_INFO << "resize step: extended_distribution done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());

  distribution_group_N_TSV.resize(ycontribution->extended_directory.size(), 0);
  distribution_group_run_result_TSV.resize(ycontribution->extended_directory.size());
  distribution_group_run_deviation_TSV.resize(ycontribution->extended_directory.size());
  distribution_group_run_N_TSV.resize(ycontribution->extended_directory.size());
  distribution_group_run_N_binwise_TSV.resize(ycontribution->extended_directory.size());
  distribution_group_run_removal_run_qTcut_TSV.resize(ycontribution->extended_directory.size());
  for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
    distribution_group_run_result_TSV[i_m].resize(ycontribution->extended_directory[i_m].size(), distribution_result_TSV);
    distribution_group_run_deviation_TSV[i_m].resize(ycontribution->extended_directory[i_m].size(), distribution_deviation_TSV);
    distribution_group_run_N_binwise_TSV[i_m].resize(ycontribution->extended_directory[i_m].size(), temp_vvvll_dbq);
    distribution_group_run_N_TSV[i_m].resize(ycontribution->extended_directory[i_m].size(), 0);
    distribution_group_run_removal_run_qTcut_TSV[i_m].resize(ycontribution->extended_directory[i_m].size(), temp_vvvi_dbq);
  }

  logger << LOG_INFO << "resize step: extended_distribution group_run done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());

  distribution_group_result_TSV.resize(ycontribution->extended_directory.size(), distribution_result_TSV);
  distribution_group_deviation_TSV.resize(ycontribution->extended_directory.size(), distribution_deviation_TSV);
  distribution_group_chi2_TSV.resize(ycontribution->extended_directory.size(), distribution_chi2_TSV);
  distribution_group_N_binwise_TSV.resize(ycontribution->extended_directory.size(), temp_vvvll_dbq);
  distribution_group_N_total_TSV.resize(ycontribution->extended_directory.size(), temp_vvvll_dbq);
  distribution_group_N_total_binwise_TSV.resize(ycontribution->extended_directory.size(), temp_vvvll_dbq);
  distribution_group_counter_removal_run_qTcut_TSV.resize(ycontribution->extended_directory.size(), temp_vvvi_dbq);
  distribution_group_counter_nonzero_run_qTcut_TSV.resize(ycontribution->extended_directory.size(), temp_vvvi_dbq);

  distribution_group_counter_run_TSV.resize(ycontribution->extended_directory.size(), 0);
    
  logger << LOG_INFO << "resize step: extended_distribution group done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
  
void summary_subprocess::readin_distribution_TSV(){
  Logger logger("summary_subprocess::readin_distribution_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  
  pid_t pid = getpid();
  stringstream temp_pid;
  temp_pid << "cat /proc/" << pid << "/status | grep VmSize";
  string temp_s_pid = temp_pid.str();
  logger << LOG_INFO << "Before initialization_distribution_TSV:   " << ycontribution->infix_contribution << " -- " << ycontribution->infix_order_contribution <<  " -- " << ycontribution->infix_path_contribution << " --- " << name << "   " << temp_s_pid << endl;
  int temp_i = system(temp_s_pid.c_str());
  						       
  initialization_distribution_TSV();
  
  logger << LOG_INFO << "After initialization_distribution_TSV:    " << ycontribution->infix_contribution << " -- " << ycontribution->infix_order_contribution <<  " -- " << ycontribution->infix_path_contribution << " --- " << name << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());
  						       

  /*
	int x_s = osi->no_reference_TSV;
	int x_r = (osi->n_scale_ren_TSV[x_s] - 1) / 2;
	int x_f = (osi->n_scale_fact_TSV[x_s] - 1) / 2;
  */
  for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
    for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
      string distribution_file;
      distribution_file = "distribution_" + name + ".txt";
      logger << LOG_DEBUG_VERBOSE << "distribution_file = " << distribution_file << endl;
      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	if (!osi->switch_distribution_TSV[i_s]){continue;}
	//	if (osi->switch_distribution_TSV[i_s] == 1){
	  //	  for (int x_q = 0; x_q < ycontribution->selection_n_qTcut; x_q++){
	    // order i_s <-> x_q switched !!!
	    
	distribution_file = "distribution_" + name + ".txt";
	string filepath = "../" + ycontribution->extended_directory[i_m][i_z] + "/distribution/" + osi->name_extended_set_TSV[i_s] + "/" + distribution_file;
	
	int no_attempts = 1;
	
	while (no_attempts < 4){
	  
	  vector<string> readin;
	  vector<vector<string> > readin_data;
	  char LineBuffer[128];
	  ifstream in_result(filepath.c_str());  
	  while (in_result.getline(LineBuffer, 128)){readin.push_back(LineBuffer);}
	  in_result.close();
	  logger << LOG_DEBUG_VERBOSE << "number of input lines in " << filepath << " = " << readin.size() << endl;
	  if (readin.size() > 0){
	    logger << LOG_DEBUG << "filepath = " << filepath << endl;
	    for (int i = 0; i < readin.size(); i++){
	      vector<string> temp_s(1);
	      int temp_comment = 0;
	      for (int i_s = 0; i_s < readin[i].size(); i_s++){
		if (readin[i][0] == '#'){temp_comment = 1; break;}
		else if (readin[i][i_s] == ' ' || readin[i][i_s] == char(9)){
		  if (temp_s[temp_s.size() - 1].size() == 0){}
		  else{temp_s.push_back("");}
		}
		else {temp_s[temp_s.size() - 1].push_back(readin[i][i_s]);}
	      }
	      if (temp_comment == 0){
		if (temp_s[temp_s.size() - 1] == ""){temp_s.erase(temp_s.end());}
		readin_data.push_back(temp_s);
	      }
	    }
	    
	    int default_size = 0;
	    int temp_distfactor = 0;
	    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){temp_distfactor += osi->extended_distribution[i_d].n_bins;}
	    //	    default_size = ycontribution->selection_n_qTcut * (1 + temp_distfactor * osi->n_scale_ren_TSV[i_s] * osi->n_scale_fact_TSV[i_s]);
	    default_size = 1 + ycontribution->selection_n_qTcut * (temp_distfactor * osi->n_scale_ren_TSV[i_s] * osi->n_scale_fact_TSV[i_s]);
	    logger << LOG_DEBUG << "filepath = " << filepath << endl;
	    if (default_size != readin_data.size() && default_size + 1 != readin_data.size()){
	      logger << LOG_INFO << "no_attempts = " << no_attempts << "   default_size = " << default_size << " != " << readin_data.size() << " = readin_data.size()" << endl;
	      no_attempts++;
	      if (no_attempts < 4){
		logger << LOG_INFO << "Wait for 10 seconds and re-try!" << endl;
		sleep(10);
	      }
	      continue;
	    }
	    ///	    default_size = readin_data.size();distribution_group_run_N_TSV[i_m][i_z]
	    if (default_size == readin_data.size() || default_size + 1 == readin_data.size()){
	      //	    if (default_size == readin_data.size()){
	      
	      if (i_s == osi->no_reference_TSV){
		distribution_group_counter_run_TSV[i_m]++;
	      }
	      
	      distribution_group_run_N_TSV[i_m][i_z] = atoll(readin_data[0][0].c_str());
	      if (default_size + 1 == readin_data.size()){n_ps = atoll(readin_data[1][0].c_str());} // n_ps
	      logger << LOG_DEBUG_VERBOSE << "distribution_group_run_N_TSV[" << i_m << "][" << i_z << "] = " << distribution_group_run_N_TSV[i_m][i_z] << endl;
	      for (int x_q = 0; x_q < ycontribution->selection_n_qTcut; x_q++){
		// x_q shifted here !!!
		int proc = 1;
		if (default_size + 1 == readin_data.size()){proc++;}
		
		for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
		  if (!ygeneric->switch_output_distribution[i_d]){proc += osi->extended_distribution[i_d].n_bins * osi->n_scale_ren_TSV[i_s] * osi->n_scale_fact_TSV[i_s]; continue;}
		  int temp_factor_qTcut = temp_distfactor * osi->n_scale_ren_TSV[i_s] * osi->n_scale_fact_TSV[i_s];
		  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
		    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
		      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
			//  double temp = atof(readin_data[proc + x_q * temp_factor_qTcut + i_b + osi->extended_distribution[i_d].n_bins * (i_r * osi->n_scale_fact_TSV[i_s] + i_f)][0].c_str());
			//  distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q] = (long long)temp;
			distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q] = atol(readin_data[proc + x_q * temp_factor_qTcut + i_b + osi->extended_distribution[i_d].n_bins * (i_r * osi->n_scale_fact_TSV[i_s] + i_f)][0].c_str());
			distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f] = atof(readin_data[proc + x_q * temp_factor_qTcut + i_b + osi->extended_distribution[i_d].n_bins * (i_r * osi->n_scale_fact_TSV[i_s] + i_f)][1].c_str());
			distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f] = atof(readin_data[proc + x_q * temp_factor_qTcut + i_b + osi->extended_distribution[i_d].n_bins * (i_r * osi->n_scale_fact_TSV[i_s] + i_f)][2].c_str());
			/*
			if (name == "du_tt~du" && osi->extended_distribution[i_d].xdistribution_name == "m_ttx-pT_t"){
			  if (i_s == x_s && i_r == x_r && i_f == x_f){
			  logger << LOG_INFO
				 << "DGRR_TSV[" << setw(3) << i_m << "][" << setw(3) << i_z << "][" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "][" << setw(3) << i_s << "][" << setw(3) << i_r << "][" << setw(3) << i_f << "] = "
				 << setprecision(15) << setw(23) << distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f] 
				 << " +- "
				 << setprecision(15) << setw(23) << distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f]
				 << " ( "
				 << setprecision(15) << setw(23) << distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q]
				 << " ) "
				 << endl;
			  }
			}
			*/
			if (i_r == 0 && i_f == 0 && i_b == 0 && x_q == 0){
			  logger << LOG_DEBUG
				 << "DGRR_TSV[" << setw(3) << i_m << "][" << setw(3) << i_z << "][" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "][" << setw(3) << i_s << "][" << setw(3) << i_r << "][" << setw(3) << i_f << "] = "
				 << setprecision(15) << setw(23) << distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f] 
				 << " +- "
				 << setprecision(15) << setw(23) << distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f]
				 << endl;
			}
		      }
		    }
		  }
		  proc += osi->extended_distribution[i_d].n_bins * osi->n_scale_ren_TSV[i_s] * osi->n_scale_fact_TSV[i_s];
		  //		  if (name == "du_tt~du" && osi->extended_distribution[i_d].xdistribution_name == "m_ttx-pT_t"){logger.newLine(LOG_INFO);}
		}
	      }
	      break;
	    }
	    else {
	      logger << LOG_DEBUG_VERBOSE << "TSV: result_file= " << filepath << " has not the expected content!" <<  endl;
	      break;
	    }
	  }
	  else {
	    logger << LOG_DEBUG_VERBOSE << "TSV: result_file = " << filepath << " is empty!" <<  endl;
	    break;
	  }
	  //	}
	}
	
      }
    }
  }
  
  /*
  //  first_DS  could be replaced by the reference TSV scale !!!
  int first_DS = 0;
  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
  if (osi->switch_distribution_TSV[i_s] != 0){first_DS = i_s; break;}}
  */
  
  for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
    for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
      distribution_group_N_TSV[i_m] += distribution_group_run_N_TSV[i_m][i_z];
      logger << LOG_DEBUG_VERBOSE << "distribution_group_N_TSV[" << i_m << "] = " << setw(15) << distribution_group_N_TSV[i_m] << "   distribution_group_run_N_TSV[" << i_m << "][" << i_z << "] = " << setw(15) << distribution_group_run_N_TSV[i_m][i_z] << endl;
      // with no runs removed at [i_d][i_b][x_q], one should later get: distribution_group_run_N_TSV[i_m] == distribution_group_N_total_TSV[i_m][i_d][i_b][x_q] !!!
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combination_distribution_TSV(){
  Logger logger("summary_subprocess::combination_distribution_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;


  vector<long long> temp_events_counted_completely(ycontribution->extended_directory.size(), 0);
  //  vector<vector<long long> > temp_events_counted_completely(ycontribution->extended_directory.size(), vector<long long> (osi->extended_distribution.size(), 0));
  int x_d = 0;
  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (osi->extended_distribution[i_d].xdistribution_type == "XS"){x_d = i_d; break;}
  }
  logger << LOG_DEBUG_VERBOSE << osi->extended_distribution[x_d] << endl;
  for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
    for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
      if (distribution_group_run_N_TSV[i_m][i_z] > 0){
	temp_events_counted_completely[i_m] += distribution_group_run_N_binwise_TSV[i_m][i_z][x_d][0][0];
      }
    }
    logger << LOG_DEBUG << name << "   temp_events_counted_completely[" << i_m << "] = " << temp_events_counted_completely[i_m] << endl;
  }

  

  // 0 replaced by reference value:
  int x_s = osi->no_reference_TSV;
  int x_r = (osi->n_scale_ren_TSV[x_s] - 1) / 2;
  int x_f = (osi->n_scale_fact_TSV[x_s] - 1) / 2;

  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    logger << LOG_DEBUG_VERBOSE << left << setw(20) << name << "  osi->extended_distribution[" << i_d << "].xdistribution_name = " << osi->extended_distribution[i_d].xdistribution_name << endl;
    /*
    if (name == "du_tt~du" && osi->extended_distribution[i_d].xdistribution_name == "m_ttx-pT_t"){
      logger << LOG_INFO << name << "   " << osi->extended_distribution[i_d].xdistribution_name << endl;
    }
    */
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
      for (int x_q = 0; x_q < ycontribution->selection_n_qTcut; x_q++){

	double tolerance_factor = ygeneric->deviation_tolerance_factor;

	for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
	  distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q] = 0;
	  vector<double> deviation_measure(ycontribution->extended_directory[i_m].size(), 0.);
	  vector<double> deviation_measure_sorted;
	  double min_deviation_measure = 1.e99;
	  double max_deviation_measure = 0.;
	  for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
	    distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] += distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q];
	    if (distribution_group_run_N_TSV[i_m][i_z] > 0){
	      deviation_measure[i_z] = sqrt(distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][x_s][x_r][x_f] - pow(distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][x_s][x_r][x_f], 2) / distribution_group_run_N_TSV[i_m][i_z]) / (distribution_group_run_N_TSV[i_m][i_z] - 1) * sqrt(distribution_group_run_N_TSV[i_m][i_z]);
	      if (deviation_measure[i_z] > 0.){
		distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q]++;
		if (deviation_measure[i_z] < min_deviation_measure){min_deviation_measure = deviation_measure[i_z];}
		if (deviation_measure[i_z] > max_deviation_measure){max_deviation_measure = deviation_measure[i_z];}
		deviation_measure_sorted.push_back(deviation_measure[i_z]);
	      }
	      /*
	      if (name == "du_tt~du" && osi->extended_distribution[i_d].xdistribution_name == "m_ttx-pT_t"){
		logger << LOG_INFO << name << "   distribution_group_run_result_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "][x_s][x_r][x_f] = " << setw(15) << setprecision(8) << distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][x_s][x_r][x_f] << " +- " << setw(15) << setprecision(8) << distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][x_s][x_r][x_f] << "   ( " << setw(10) << distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q] << " )" << endl;
	      }
	      */
	      //	      else {logger << LOG_WARN << "deviation_measure[" << i_z << "] = " << deviation_measure[i_z] << endl;}
	    }
	  }
	  /*
	  if (name == "du_tt~du" && osi->extended_distribution[i_d].xdistribution_name == "m_ttx-pT_t"){
	    logger << LOG_INFO << name << "   distribution_group_N_binwise_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "] = " << distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] << endl;
	    //	    logger.newLine(LOG_INFO);
	  }
	  */
	  
	  logger << LOG_DEBUG_VERBOSE << "i_m = " << i_m << "   deviation_measure_sorted.size() = " << deviation_measure_sorted.size() << endl;

	  if (deviation_measure_sorted.size() > 0){
	    sort(deviation_measure_sorted.begin(), deviation_measure_sorted.end());
	    // in order to ignore runs with too low error estimates:
	    min_deviation_measure = deviation_measure_sorted[deviation_measure_sorted.size() / 10];	  

	    logger << LOG_DEBUG_VERBOSE << "distribution_group_N_binwise_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "] = " << distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] << endl;
	    logger << LOG_DEBUG_VERBOSE << "distribution_group_N_TSV[" << i_m << "] = " << distribution_group_N_TSV[i_m] << endl;
	    
    	    double temp_fraction = double(distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q]) / (distribution_group_N_TSV[i_m] * n_ps);
	    tolerance_factor = ygeneric->deviation_tolerance_factor * (1. - log10(temp_fraction));
	    //	    tolerance_factor = tolerance_factor * (1. - log10(temp_fraction));
	    logger << LOG_DEBUG_VERBOSE << "temp_fraction = " << temp_fraction << endl;
	    logger << LOG_DEBUG_VERBOSE << "tolerance_factor = " << tolerance_factor << endl;

	    /*
	    long long temp_events_in_bin = 0;
	    for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
	      temp_events_in_bin += distribution_group_run_N_TSV[i_m][i_z];
	    }
	    // Doesn't work well because temp_events_in_bin only counts the number of accepted events, which is not the same as what is counted for the distributions because dipole events are counted separately...
	    */
	    
	    for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
	      // replace absolute threshold 1000 by 100 * average(N_events in bin) of all runs contributing ???
	      // try 25 instead of 100
	      //	      if (distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] > 25 * temp_events_in_bin / distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q] && deviation_measure[i_z] > tolerance_factor * min_deviation_measure){
	      
	      //	      logger << LOG_DEBUG_VERBOSE << "distribution_group_N_binwise_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "] = " << distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] << endl;
	      logger << LOG_DEBUG_VERBOSE << "distribution_group_N_binwise_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "] = " << distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] << endl;
	      logger << LOG_DEBUG_VERBOSE << "  temp_events_counted_completely[" << i_m << "][ ][ ][ ] = " << temp_events_counted_completely[i_m] << endl;
	      logger << LOG_DEBUG_VERBOSE << "distribution_group_counter_nonzero_run_qTcut_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "]  = " << distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q] << endl;

	      logger << LOG_DEBUG_VERBOSE << "deviation_measure[" << setw(4) << i_z << "] = " << deviation_measure[i_z] << endl;
	      logger << LOG_DEBUG_VERBOSE << "tolerance_factor        = " << tolerance_factor << endl;
	      logger << LOG_DEBUG_VERBOSE << "min_deviation_measure   = " << min_deviation_measure << endl;
	      logger << LOG_DEBUG_VERBOSE << "tolerance_factor * min_deviation_measure = " << tolerance_factor * min_deviation_measure << endl;
	      /*
	      if (name == "du_tt~du" && osi->extended_distribution[i_d].xdistribution_name == "m_ttx-pT_t"){
		logger << LOG_INFO << "distribution_group_N_binwise_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "] = " << distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] << endl;
		logger << LOG_INFO << "  temp_events_counted_completely[" << i_m << "][ ][ ][ ] = " << temp_events_counted_completely[i_m] << endl;
		logger << LOG_INFO << "distribution_group_counter_nonzero_run_qTcut_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "]  = " << distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q] << endl;
		
		logger << LOG_INFO << "deviation_measure[" << setw(4) << i_z << "] = " << deviation_measure[i_z] << endl;
		logger << LOG_INFO << "tolerance_factor        = " << tolerance_factor << endl;
		logger << LOG_INFO << "min_deviation_measure   = " << min_deviation_measure << endl;
		logger << LOG_INFO << "tolerance_factor * min_deviation_measure = " << tolerance_factor * min_deviation_measure << endl;
	      }
	      */
	      /*
	      if (name == "du_tt~du" && osi->extended_distribution[i_d].xdistribution_name == "m_ttx-pT_t"){
		logger << LOG_INFO << "distribution_group_N_binwise_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "] = " << distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] << endl;
		logger << LOG_INFO << "distribution_group_counter_nonzero_run_qTcut_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "] = " << distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q] << endl;
		logger << LOG_INFO << "temp_events_counted_completely[" << i_m << "] / distribution_group_counter_nonzero_run_qTcut_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "] = " << temp_events_counted_completely[i_m] / distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q] << endl;
		logger << LOG_INFO << "deviation_measure[i_z] = " << deviation_measure[i_z] << " > " << tolerance_factor * min_deviation_measure << " = tolerance_factor * min_deviation_measure" << endl;    
	      }
	      */
	      
	      if (distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] > 25 * temp_events_counted_completely[i_m] / distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q] && deviation_measure[i_z] > tolerance_factor * min_deviation_measure){

		//temp_events_counted_completely[i_m]
		
		//	      if (distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] > 1000 && deviation_measure[i_z] > tolerance_factor * min_deviation_measure){
		logger << LOG_INFO << "Only run is removed --- should not happen, of course !!!" << endl;
		logger << LOG_DEBUG_VERBOSE << "distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] = " << distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] << endl;
		logger << LOG_DEBUG_VERBOSE << "temp_events_counted_completely[i_m] = " << temp_events_counted_completely[i_m] << endl;
		//		logger << LOG_DEBUG_VERBOSE << "temp_events_in_bin = " << temp_events_in_bin << endl;
		logger << LOG_DEBUG_VERBOSE << "distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q] = " << distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q] << endl;
		logger << LOG_DEBUG_VERBOSE << "deviation_measure[i_z] = " << deviation_measure[i_z] << endl;
		logger << LOG_DEBUG_VERBOSE << "tolerance_factor = " << tolerance_factor << endl;
		logger << LOG_DEBUG_VERBOSE << "min_deviation_measure = " << min_deviation_measure << endl;
		//		logger << LOG_DEBUG_VERBOSE << "(deviation_measure[i_z] > tolerance_factor * min_deviation_measure) = " << (deviation_measure[i_z] > tolerance_factor * min_deviation_measure) << endl;
		//	  logger << LOG_DEBUG_VERBOSE << "(distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] > 25 * temp_events_in_bin / distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q]) = " << (distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] > 25 * temp_events_in_bin / distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q]) << endl;
		
		distribution_group_run_removal_run_qTcut_TSV[i_m][i_z][i_d][i_b][x_q] = 1;
		distribution_group_counter_removal_run_qTcut_TSV[i_m][i_d][i_b][x_q]++;
		/*
		if (name == "du_tt~du" && osi->extended_distribution[i_d].xdistribution_name == "m_ttx-pT_t"){
		  logger << LOG_INFO << "distribution_group_counter_removal_run_qTcut_TSV[" << i_m << "][" << i_d << "][" << i_b << "][" << x_q << "]  = " << distribution_group_counter_removal_run_qTcut_TSV[i_m][i_d][i_b][x_q] << endl;
		}
		*/
	      }
	      else {
		distribution_group_N_total_TSV[i_m][i_d][i_b][x_q] += distribution_group_run_N_TSV[i_m][i_z];
		distribution_group_N_total_binwise_TSV[i_m][i_d][i_b][x_q] += distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q];
	      }
	    }
	  }
	}

	for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
	  distribution_N_total_TSV[i_d][i_b][x_q] += distribution_group_N_total_TSV[i_m][i_d][i_b][x_q];
	  distribution_N_total_binwise_TSV[i_d][i_b][x_q] += distribution_group_N_total_binwise_TSV[i_m][i_d][i_b][x_q];
	}
 
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  if (!osi->switch_distribution_TSV[i_s]){continue;}
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      if (ycontribution->average_factor == 1){
		combine_distribution_conservative_TSV(x_q, i_d, i_b, i_s, i_r, i_f);
	      }
	      else {
		//	      if (0 > 100){
		combine_distribution_hybrid_TSV(x_q, i_d, i_b, i_s, i_r, i_f);
		/*
		// hybrid version should include both extreme cases.
		// check why the standard conservative method does not work !!!
		if (ycontribution->average_factor == 0){combine_distribution_aggressive_TSV(x_z, x_q, i_d, i_b, i_s, i_r, i_f);}
		//		    else if (ycontribution->average_factor == 1){combine_distribution_conservative_TSV(x_z, x_q, i_d, i_b, i_s, i_r, i_f);}
		else {combine_distribution_hybrid_TSV(x_q, i_d, i_b, i_s, i_r, i_f);}
		*/
	      }
	      /*
	      else {
		combine_distribution_conservative_TSV(x_q, i_d, i_b, i_s, i_r, i_f);
	      }
	      */
	      
	      //  weighted combination of results from different extended_directory contributions:
	      //  should receive a mechanism to prevent underestimated errors from low-statistics contributions to get much weight !!!

	      double result_subprocess = 0.;
	      double deviation_subprocess = 0.;
	      vector<double> deltasigma(ycontribution->extended_directory.size(), 0.);
	      vector<double> singleweight(ycontribution->extended_directory.size(), 0.);
	      for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
		/*
		if (name == "du_tt~du" && osi->extended_distribution[i_d].xdistribution_name == "m_ttx-pT_t"){
		  if (i_s == x_s && i_r == x_r && i_f == x_f){
		    logger << LOG_INFO << "DGR_TSV[" << setw(3) << i_m << "][" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "][" << setw(3) << i_s << "][" << setw(3) << i_r << "][" << setw(3) << i_f << "] = " << setw(15) << setprecision(8) << distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] << " +- " << setw(15) << setprecision(9) << distribution_group_deviation_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] << endl;
		  }
		}
		*/
		logger << LOG_DEBUG_VERBOSE << "distribution_group_result_TSV[" << setw(3) << i_m << "][" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "][" << setw(3) << i_s << "][" << setw(3) << i_r << "][" << setw(3) << i_f << "] = " << setw(15) << setprecision(8) << distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] << " +- " << setw(15) << setprecision(9) << distribution_group_deviation_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] << endl;
		// simplified version available ???
		//		if (distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] > 0){
		/*
		if (name == "du_tt~du" && osi->extended_distribution[i_d].xdistribution_name == "m_ttx-pT_t"){
		  if (i_s == x_s && i_r == x_r && i_f == x_f){
		    logger << LOG_INFO << "distribution_group_N_total_binwise_TSV[" << setw(3) << i_m << "][" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "] = " << distribution_N_total_binwise_TSV[i_d][i_b][x_q] << "/ 100 > " << "distribution_N_total_binwise_TSV[" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "] / 100" << endl;
		  }
		}
		*/

	        if (distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] > 0 &&
		    distribution_group_deviation_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] != 0. &&
		    distribution_group_N_total_binwise_TSV[i_m][i_d][i_b][x_q] > distribution_N_total_binwise_TSV[i_d][i_b][x_q] / 100){
		  ///		    distribution_group_N_total_TSV[i_m][i_d][i_b][x_q] > distribution_N_total_TSV[i_d][i_b][x_q] / 100){
		  //		singleweight[i_m] = 1. / distribution_group_deviation2_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f];
		  singleweight[i_m] = 1. / pow(distribution_group_deviation_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f], 2);
		  result_subprocess += singleweight[i_m] * distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f];
		  deviation_subprocess += singleweight[i_m];
		}
	      }
	      if (deviation_subprocess != 0.){
		distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f] = result_subprocess / deviation_subprocess;
		/// simplified version that does not take into account the difference of central values in determination of new standard deviation...
		distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f] = 1. / sqrt(deviation_subprocess);
		/*
		double final_deviation = 0.;
		for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
		  final_deviation += (singleweight[i_m] + 4 * pow(singleweight[i_m] * (distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] - distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f]), 2)) / pow(deviation_subprocess, 2);
		}	   
		distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f] = sqrt(final_deviation);
		*/
	      }
	      else {
		distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f] = 0.;
		distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f] = 0.;
	      }
	      /*
		if (name == "du_tt~du" && osi->extended_distribution[i_d].xdistribution_name == "m_ttx-pT_t"){
		  if (i_s == x_s && i_r == x_r && i_f == x_f){
		    logger << LOG_INFO << "DR_TSV[" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "][" << setw(3) << i_s << "][" << setw(3) << i_r << "][" << setw(3) << i_f << "] = " << setw(15) << setprecision(8) << distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f] << " +- " << setw(15) << setprecision(9) << distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f] << endl;
		  //" ( " << setw(10) << distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] << " ) "
		  }
		}
	      */
	      
	      logger << LOG_DEBUG_VERBOSE << "           distribution_result_TSV[" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "][" << setw(3) << i_s << "][" << setw(3) << i_r << "][" << setw(3) << i_f << "] = " << setw(15) << setprecision(8) << distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f] << " +- " << setw(15) << setprecision(9) << distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f] << endl;
	    
	      //  chi2 calculation from [i_m] combination:
	      distribution_chi2_TSV[i_d][i_b][x_q][i_s][i_r][i_f] = 0.;
	      logger << LOG_DEBUG_VERBOSE << "chi2 evaluation in " << ycontribution->resultdirectory << " -- " << ycontribution->infix_contribution << " -- " << ycontribution->infix_order_contribution <<  " -- " << ycontribution->infix_path_contribution << " --- " << name << endl;
	      logger << LOG_DEBUG_VERBOSE << "combined result = " << setw(15) << setprecision(9) << distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f] << " +- " << setw(15) << setprecision(8) << distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f] << endl;

	      int temp_counter = 0;
	      for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
		if (distribution_group_N_TSV[i_m] > 0){
		  if (!(distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] == 0. && distribution_group_deviation_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] == 0.)){
		    temp_counter++;
		    distribution_chi2_TSV[i_d][i_b][x_q][i_s][i_r][i_f] += pow((distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] - distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f]) / distribution_group_deviation_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f], 2.);
		  }
		}
	      }
	      if (temp_counter > 1){
		distribution_chi2_TSV[i_d][i_b][x_q][i_s][i_r][i_f] = distribution_chi2_TSV[i_d][i_b][x_q][i_s][i_r][i_f] / temp_counter;
	      }
	      else {
		//  chi2 cannot be reasonably defined for only one run.
		distribution_chi2_TSV[i_d][i_b][x_q][i_s][i_r][i_f] = 0.;
	      }
	      logger << LOG_DEBUG_VERBOSE << "distribution_chi2_TSV[" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "][" << setw(3) << i_s << "][" << setw(3) << i_r << "][" << setw(3) << i_f << "] = " << distribution_chi2_TSV[i_d][i_b][x_q][i_s][i_r][i_f] << " (" << temp_counter << ")" << endl;

	      
	      if (x_s == i_s && x_r == i_r && x_f == i_f){
		for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
		  if (osi->extended_distribution[i_d].xdistribution_name == "total_rate" && distribution_group_chi2_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] > 10. && distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] > 1000){
		  logger << LOG_INFO << "WARNING : too large  chi2_dof = " << distribution_group_chi2_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] << "  for  " << name << "  in  " << ycontribution->infix_order_contribution << endl;
		  //		logger << LOG_INFO << "WARNING : too large  chi2_dof = " << distribution_chi2_TSV[i_d][i_b][x_q][i_s][i_r][i_f] << "  for  " << name << "  in  " << ycontribution->infix_order_contribution << endl;
		    logger << LOG_INFO << "chi2_TSV[" << setw(3) << i_m << "][" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "][" << setw(3) << i_s << "][" << setw(3) << i_r << "][" << setw(3) << i_f << "] = " << showpoint << setw(15) << setprecision(8) << distribution_group_chi2_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] << " (" << temp_counter << ")" << endl;
		    
		    logger << LOG_INFO<< "DGR_TSV[" << setw(3) << i_m << "][" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "][" << setw(3) << i_s << "][" << setw(3) << i_r << "][" << setw(3) << i_f << "] = " << showpoint << setw(15) << setprecision(8) << distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] << " +- " << setw(15) << setprecision(9) << distribution_group_deviation_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] << endl;
		    
		    for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
		      double temp_result = distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][x_s][x_r][x_f] / distribution_group_run_N_TSV[i_m][i_z];
		      double temp_deviation = sqrt((distribution_group_run_N_TSV[i_m][i_z] * distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][x_s][x_r][x_f] - pow(temp_result, 2)) / distribution_group_run_N_TSV[i_m][i_z]) / (distribution_group_run_N_TSV[i_m][i_z] - 1);
		      logger << LOG_INFO << "DGRR_TSV[" << setw(3) << i_m << "][" << setw(3) << i_z << "][" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_q << "][" << setw(3) << i_s << "][" << setw(3) << i_r << "][" << setw(3) << i_f << "] = " << showpoint << setw(17) << setprecision(10) << temp_result << " +- " << temp_deviation << "   ( " << showpoint << setw(15) << setprecision(8) << distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q] << " - " << setw(8) << distribution_group_run_N_TSV[i_m][i_z] << " )   [ " << distribution_group_run_removal_run_qTcut_TSV[i_m][i_z][i_d][i_b][x_q] << " ] " << endl;
		      
		      //	    deviation_measure[i_z] = sqrt(distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][x_s][x_r][x_f] - pow(distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][x_s][x_r][x_f], 2) / distribution_group_run_N_TSV[i_m][i_z]) / (distribution_group_run_N_TSV[i_m][i_z] - 1) * sqrt(distribution_group_run_N_TSV[i_m][i_z]);
		    }
		  }
		}
	      }

	    }
	  }
	}
      }
    }
  }

  // new:
  // some more such clean-ups to avoid too much memory consumption !!!
    pid_t pid = getpid();
  stringstream temp_pid;
  temp_pid << "cat /proc/" << pid << "/status | grep VmSize";
  string temp_s_pid = temp_pid.str();
  logger << LOG_INFO << "before clearing via swap:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  int temp_i = system(temp_s_pid.c_str());

  vector<vector<vector<long long> > > ().swap(distribution_N_total_TSV);
  logger << LOG_INFO << "swap(distribution_N_total_TSV) done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());


  vector<vector<vector<long long> > > ().swap(distribution_N_total_binwise_TSV);
  logger << LOG_INFO << "swap(distribution_N_total_binwise_TSV) done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());


  vector<vector<vector<vector<long long> > > > ().swap(distribution_group_N_total_TSV);
  logger << LOG_INFO << "swap(distribution_group_N_total_TSV) done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());


  vector<vector<vector<vector<long long> > > > ().swap(distribution_group_N_total_binwise_TSV);
    logger << LOG_INFO << "swap(distribution_group_N_total_binwise_TSV) done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());


  vector<vector<long long> > ().swap(distribution_group_run_N_TSV);
  logger << LOG_INFO << "swap(distribution_group_run_N_TSV) done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());


  vector<vector<vector<vector<vector<int> > > > > ().swap(distribution_group_run_removal_run_qTcut_TSV);
  logger << LOG_INFO << "swap(distribution_group_run_removal_run_qTcut_TSV) done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());


  vector<vector<vector<vector<vector<long long> > > > > ().swap(distribution_group_run_N_binwise_TSV);
  logger << LOG_INFO << "swap(distribution_group_run_N_binwise_TSV) done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());


  vector<vector<vector<vector<vector<vector<vector<vector<double> > > > > > > > ().swap(distribution_group_run_result_TSV);
  logger << LOG_INFO << "swap(distribution_group_run_result_TSV) done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());


  vector<vector<vector<vector<vector<vector<vector<vector<double> > > > > > > > ().swap(distribution_group_run_deviation_TSV);
  logger << LOG_INFO << "distribution_group_run_deviation_TSV done:   " << ycontribution->infix_contribution << "   " << temp_s_pid << endl;
  temp_i = system(temp_s_pid.c_str());


  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combine_distribution_conservative_TSV(int x_q, int i_d, int i_b, int i_s, int i_r, int i_f){
  Logger logger("summary_subprocess::combine_distribution_conservative_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
    double temp_result_subprocess = 0.;
    double temp_deviation2_subprocess = 0.;
    for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
      logger << LOG_DEBUG_VERBOSE << "distribution_group_run_removal_run_qTcut_TSV[" << i_m << "][" << i_z << "][" << i_d << "][" << i_b << "][" << x_q << "] = " << distribution_group_run_removal_run_qTcut_TSV[i_m][i_z][i_d][i_b][x_q] << "   distribution_group_run_N_binwise_TSV = " << distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q] << endl;
      if (!distribution_group_run_removal_run_qTcut_TSV[i_m][i_z][i_d][i_b][x_q] && distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q] > 0){
	///      if (!distribution_group_run_removal_run_qTcut_TSV[i_m][i_z][i_d][i_b][x_q] && distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f] > 0){
	// better: use only 'reference' version thereof ???
	// repeated for all i_s, i_r, i_f !!! (should be counted only once - if run is not removed ??? ) -> comment out !!!
	//	distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] += distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q];
	///	distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] += distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f];
	temp_result_subprocess += distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f];
	  // should be named sum_result...
	temp_deviation2_subprocess += distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f];
	  // should be named sum_deviation²...
      }
    }
    // difference to CV??? Maybe because of errors in CV counter !!! ???
    if (distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q] != 0){
      ///    if (distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] != 0){
      distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] = temp_result_subprocess / distribution_group_N_total_TSV[i_m][i_d][i_b][x_q];
      // distribution_group_N_total_TSV[i_m][i_d][i_b][x_q]  instead of  distribution_group_N_TSV[i_m]  since the latter contains also results from removed runs
      //  distribution_group_deviation2_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] = (distribution_group_N_TSV[i_m] * temp_deviation2_subprocess - pow(temp_result_subprocess, 2)) / distribution_group_N_TSV[i_m] / pow(distribution_group_N_TSV[i_m] - 1, 2); // possible replacement - check if correct !!!
      distribution_group_deviation_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] = sqrt((distribution_group_N_total_TSV[i_m][i_d][i_b][x_q] * temp_deviation2_subprocess - pow(temp_result_subprocess, 2)) / distribution_group_N_total_TSV[i_m][i_d][i_b][x_q]) / (distribution_group_N_total_TSV[i_m][i_d][i_b][x_q] - 1);
    }
    else {
      distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] = 0.;
      distribution_group_deviation_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] = 0.;
    }




    //  chi2 calculation:
    distribution_group_chi2_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] = 0.;
    logger << LOG_DEBUG_VERBOSE << "chi2 evaluation in " << ycontribution->resultdirectory << " -- " << ycontribution->infix_contribution << " -- " << ycontribution->infix_order_contribution <<  " -- " << ycontribution->infix_path_contribution << " --- " << name << endl;
    logger << LOG_DEBUG_VERBOSE << "combined result = " << setw(15) << setprecision(9) << distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] << " +- " << setw(15) << setprecision(9) << distribution_group_deviation_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] << endl;

    if (distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q] > 1){
      for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
	if (!(distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f] == 0. && distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f] == 0.)){
	  // expectation value -> temp_run_result
	  // sum (result) ->distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f]
	  double temp_run_result = distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f] / distribution_group_run_N_TSV[i_m][i_z]; 
	  // standard deviation -> temp_run_deviation
	  // sum (result²) ->distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f]
	  double temp_run_deviation = sqrt(distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f] - pow(distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f], 2) / distribution_group_run_N_TSV[i_m][i_z]) / (distribution_group_run_N_TSV[i_m][i_z] - 1);
	  distribution_group_chi2_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] += pow((temp_run_result - distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f]) / temp_run_deviation, 2.);
	}
      }
      distribution_group_chi2_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] = distribution_group_chi2_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] / distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q];
    }
    else {
      // chi2 cannot be reasonably defined for only one run.
      distribution_group_chi2_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] = 0.;
    }
  }
    
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combine_distribution_hybrid_TSV(int x_q, int i_d, int i_b, int i_s, int i_r, int i_f){
  Logger logger("summary_subprocess::combine_distribution_hybrid_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  /*////
  int n_conservative = ycontribution->average_factor;
  int temp_n_conservative = 0;
  int temp_n_runs = 0;
  vector<int> temp_list_no_run;
  for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
    if (!remove_run_bin[i_z]){
      if (distribution_run_N_TSV[i_z] != 0){
	temp_n_runs++;
	temp_list_no_run.push_back(i_z);
      }
    }
  }
  if (temp_n_runs < n_conservative){temp_n_conservative = temp_n_runs;}
  else {temp_n_conservative = n_conservative;}
  if (ycontribution->average_factor == 0){temp_n_conservative = temp_n_runs;}

  vector<double> conservative_result_run_TSV(temp_n_conservative, 0.);
  vector<double> conservative_deviation_run_TSV(temp_n_conservative, 0.);
  vector<double> conservative_chi2_TSV(temp_n_conservative, 0.);
  vector<long long> conservative_N_all(temp_n_conservative, 0);
  vector<long long> conservative_no_results(temp_n_conservative, 0);
  
  for (int i_z = 0; i_z < temp_list_no_run.size(); i_z++){
    int conservative_i_z = i_z % temp_n_conservative;

    if (!(distribution_run_result_TSV[temp_list_no_run[i_z]][i_d][i_b][x_q][i_s][i_r][i_f] == 0. && distribution_run_deviation_TSV[temp_list_no_run[i_z]][i_d][i_b][x_q][i_s][i_r][i_f] == 0.)){
      ///      if (!remove_run_bin[temp_list_no_run[i_z]]){
      temp_N_binwise += distribution_run_N_binwise_TSV[temp_list_no_run[i_z]][i_d][i_b][x_q][i_s][i_r][i_f];
      conservative_N_all[conservative_i_z] += distribution_run_N_TSV[temp_list_no_run[i_z]];
      conservative_no_results[conservative_i_z]++;
      // input different for distributions
      conservative_result_run_TSV[conservative_i_z] += distribution_run_result_TSV[temp_list_no_run[i_z]][i_d][i_b][x_q][i_s][i_r][i_f];
      // distribution_run_result_TSV   should be sum over weights
      conservative_deviation_run_TSV[conservative_i_z] += distribution_run_deviation_TSV[temp_list_no_run[i_z]][i_d][i_b][x_q][i_s][i_r][i_f];
      // distribution_run_result_TSV   should be sum over squared weights
      // logger << LOG_DEBUG << "conservative_N_all[" << setw(3) << conservative_i_z << "] = " << setw(10) << conservative_N_all[conservative_i_z] << "   distribution_run_N_TSV[" << setw(3) << i_z << "] = " << setw(10) << distribution_run_N_TSV[i_z] << "   result = " << distribution_run_result_TSV[i_z][i_d][i_b][x_q][i_s][i_r][i_f] << " +- " << distribution_run_deviation_TSV[i_z][i_d][i_b][x_q][i_s][i_r][i_f] << endl;
	///      }
	///      else {logger << LOG_INFO << "Should never happen !!!" << endl;}
    }
  }

  for (int conservative_i_z = 0; conservative_i_z < temp_n_conservative; conservative_i_z++){
    conservative_deviation_run_TSV[conservative_i_z] = sqrt((conservative_N_all[conservative_i_z] * conservative_deviation_run_TSV[conservative_i_z] - pow(conservative_result_run_TSV[conservative_i_z], 2)) / conservative_N_all[conservative_i_z]) / (conservative_N_all[conservative_i_z] - 1); 
    conservative_result_run_TSV[conservative_i_z] = conservative_result_run_TSV[conservative_i_z] / conservative_N_all[conservative_i_z];
  }
  for (int conservative_i_z = 0; conservative_i_z < temp_n_conservative; conservative_i_z++){
    //    conservative_chi2_TSV[conservative_i_z] = conservative_chi2_TSV[conservative_i_z] / conservative_no_results[conservative_i_z];
    if (conservative_no_results[conservative_i_z] == 0){
      conservative_result_run_TSV[conservative_i_z] = 0.;
      conservative_deviation_run_TSV[conservative_i_z] = 0.;
    }
  }
  
  combination_N_event_all = 0;
  for (int conservative_i_z = 0; conservative_i_z < temp_n_conservative; conservative_i_z++){
    combination_N_event_all += conservative_N_all[conservative_i_z];
    if (conservative_N_all[conservative_i_z] > 0){no_results++;}
  }
  if (combination_N_event_all != 0){
    combination_result[x_z] = 0.;
    combination_deviation[x_z] = 0.;
    for (int i_z = 0; i_z < temp_n_conservative; i_z++){
      if (distribution_run_N_TSV[i_z] != 0 && conservative_result_run_TSV[i_z] != 0. && conservative_deviation_run_TSV[i_z] != 0.){
	double temp_singleweight = 1. / pow(conservative_deviation_run_TSV[i_z], 2);
	logger << LOG_DEBUG_VERBOSE << "z_result_TSV[" << i_z << "] = " << setw(23) << setprecision(15) << conservative_result_run_TSV[i_z] << " +- = " << setw(23) << setprecision(15) << conservative_deviation_run_TSV[i_z] << endl;
	combination_result[x_z] += temp_singleweight * conservative_result_run_TSV[i_z];
	combination_deviation[x_z] += temp_singleweight;
      }
      else if (distribution_run_N_TSV[i_z] != 0 && conservative_result_run_TSV[i_z] == 0. && conservative_deviation_run_TSV[i_z] == 0.){
	combination_result[x_z] += 0.;
	combination_deviation[x_z] += 0.;
	logger << LOG_DEBUG_VERBOSE << "empty: z_result_TSV[" << i_z << "] = " << setw(23) << setprecision(15) << conservative_result_run_TSV[i_z] << " +- = " << setw(23) << setprecision(15) << conservative_deviation_run_TSV[i_z] << endl;
      }
    }
    if (combination_N_event_all != 0 && combination_result[x_z] != 0. && combination_deviation[x_z] != 0.){
      combination_result[x_z] = combination_result[x_z] / combination_deviation[x_z];
      combination_deviation[x_z] = 1. / sqrt(combination_deviation[x_z]);
    }
    else {
      combination_result[x_z] = 0.;
      combination_deviation[x_z] = 0.;
    }
  }
  else if (combination_N_event_all == 0){
    combination_result[x_z] = 0.;
    combination_deviation[x_z] = 0.;
  }
  
  for (int conservative_i_z = 0; conservative_i_z < temp_n_conservative; conservative_i_z++){
    if (!(distribution_run_result_TSV[conservative_i_z][i_d][i_b][x_q][i_s][i_r][i_f] == 0. && 
	  distribution_run_deviation_TSV[conservative_i_z][i_d][i_b][x_q][i_s][i_r][i_f] == 0.)){
      conservative_chi2_TSV[conservative_i_z] += pow((conservative_result_run_TSV[conservative_i_z] - combination_result[x_z]) / conservative_deviation_run_TSV[conservative_i_z], 2.);
    }
    conservative_chi2_TSV[conservative_i_z] = conservative_chi2_TSV[conservative_i_z] / temp_n_conservative;
  }
  *////  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}












void summary_subprocess::initialization_distribution_CV(){
  Logger logger("summary_subprocess::initialization_distribution_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  distribution_result_CV.resize(osi->extended_distribution.size());
  distribution_deviation_CV.resize(osi->extended_distribution.size());
  distribution_chi2_CV.resize(osi->extended_distribution.size());

  distribution_N_total_CV.resize(osi->extended_distribution.size());
  distribution_N_total_binwise_CV.resize(osi->extended_distribution.size());

  // old version:
  /////  distribution_N_CV.resize(osi->extended_distribution.size(), 0);
  /////  distribution_N_binwise_CV.resize(osi->extended_distribution.size());

  // in analogy to TSV:

  vector<vector<long long> > temp_vvll_db(osi->extended_distribution.size());
  //  vector<vector<long long> > temp_vvll_db(osi->extended_distribution.size());
  vector<vector<int> > temp_vvi_db(osi->extended_distribution.size());

  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    temp_vvll_db[i_d].resize(osi->extended_distribution[i_d].n_bins, 0);
    //    temp_vvll_db[i_d].resize(osi->extended_distribution[i_d].n_bins);
    temp_vvi_db[i_d].resize(osi->extended_distribution[i_d].n_bins, 0);
    
    distribution_result_CV[i_d].resize(osi->extended_distribution[i_d].n_bins, vector<double> (osi->n_scales_CV, 0.));
    distribution_deviation_CV[i_d].resize(osi->extended_distribution[i_d].n_bins, vector<double> (osi->n_scales_CV, 0.));
    distribution_chi2_CV[i_d].resize(osi->extended_distribution[i_d].n_bins, vector<double> (osi->n_scales_CV, 0.));

    // total is also binwise, but after removal of runs...
    distribution_N_total_CV[i_d].resize(osi->extended_distribution[i_d].n_bins, 0);
    distribution_N_total_binwise_CV[i_d].resize(osi->extended_distribution[i_d].n_bins, 0);

    /////    distribution_N_binwise_CV[i_d].resize(osi->extended_distribution[i_d].n_bins, 0);
  }

  distribution_group_N_CV.resize(ycontribution->extended_directory.size(), 0);
  distribution_group_run_result_CV.resize(ycontribution->extended_directory.size());
  distribution_group_run_deviation_CV.resize(ycontribution->extended_directory.size());
  distribution_group_run_N_CV.resize(ycontribution->extended_directory.size());
  distribution_group_run_N_binwise_CV.resize(ycontribution->extended_directory.size());
  distribution_group_run_removal_run_CV.resize(ycontribution->extended_directory.size());

  for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
    distribution_group_run_result_CV[i_m].resize(ycontribution->extended_directory[i_m].size(), distribution_result_CV);
    distribution_group_run_deviation_CV[i_m].resize(ycontribution->extended_directory[i_m].size(), distribution_deviation_CV);
    distribution_group_run_N_binwise_CV[i_m].resize(ycontribution->extended_directory[i_m].size(), temp_vvll_db);
    distribution_group_run_N_CV[i_m].resize(ycontribution->extended_directory[i_m].size(), 0);
    distribution_group_run_removal_run_CV[i_m].resize(ycontribution->extended_directory[i_m].size(), temp_vvi_db);
  }

  distribution_group_result_CV.resize(ycontribution->extended_directory.size(), distribution_result_CV);
  distribution_group_deviation_CV.resize(ycontribution->extended_directory.size(), distribution_deviation_CV);
  distribution_group_chi2_CV.resize(ycontribution->extended_directory.size(), distribution_chi2_CV);
  distribution_group_N_binwise_CV.resize(ycontribution->extended_directory.size(), temp_vvll_db);
  distribution_group_N_total_CV.resize(ycontribution->extended_directory.size(), temp_vvll_db);
  distribution_group_N_total_binwise_CV.resize(ycontribution->extended_directory.size(), temp_vvll_db);
  distribution_group_counter_removal_run_CV.resize(ycontribution->extended_directory.size(), temp_vvi_db);
  distribution_group_counter_nonzero_run_CV.resize(ycontribution->extended_directory.size(), temp_vvi_db);



  
  /*
  // old version:
  distribution_run_N_CV.resize(ycontribution->directory.size(), distribution_N_CV);
  //  (roughly) equivalent of distribution_run_N_TSV
  distribution_run_N_binwise_CV.resize(ycontribution->directory.size(), distribution_N_binwise_CV);
  distribution_run_result_CV.resize(ycontribution->directory.size(), distribution_result_CV);
  distribution_run_deviation_CV.resize(ycontribution->directory.size(), distribution_deviation_CV);
  //  distribution_N_TSV = 0;
  // !!! missing: equivalent of distribution_N_TSV

  all_distribution_N_binwise_CV.resize(osi->extended_distribution.size());

  
  combination_all_distribution_N_binwise_CV.resize(2, vector<vector<long long> > (osi->extended_distribution.size()));
  // newly introduced: equivalent of  all_distribution_N_binwise_TSV  and  combination_all_distribution_N_binwise_TSV

  //  distribution_run_deltasigma.resize(ycontribution->directory.size(), distribution_result_CV);
  //  distribution_run_weight.resize(ycontribution->directory.size(), distribution_result_CV);
  */


  

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void summary_subprocess::readin_distribution_CV(){
  Logger logger("summary_subprocess::readin_distribution_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

 initialization_distribution_CV();


  for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
    for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  string filename = "../" + ycontribution->extended_directory[i_m][i_z] + "/distribution/CV/" + osi->directory_name_scale_CV[i_s] + "/" + osi->extended_distribution[i_d].xdistribution_name + "_" + name + ".dat";
	  //	  string filename = "../" + ycontribution->directory[i_z] + "/distribution/CV/" + osi->directory_name_scale_CV[i_s] + "/" + osi->extended_distribution[i_d].xdistribution_name + "_" + name + ".dat";
	  vector<string> readin;
	  char LineBuffer[128];
	  ifstream in_result(filename.c_str());  
	  while (in_result.getline(LineBuffer, 128)){readin.push_back(LineBuffer);}
	  in_result.close();
	  logger << LOG_DEBUG_VERBOSE << filename << "   number of lines: " << setw(5) << readin.size() << endl;
	  if (readin.size() != 0){
	    // old version: to be adapted !!!
	    // should be identical for each distribution, thus independent of [i_d] (is overwritten for each i_d):
	    distribution_group_run_N_CV[i_m][i_z] = atoll(readin[0].c_str());
	    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	      distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b] = atoll(readin[1 + i_b * 3].c_str());

	      // temporary (because N is doubled):
	      distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b] = distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b] / 2;
	      
	      //  N[i_z][i_p][i_d][i_b] = N[i_z][i_p][i_d][i_b] / (osi->n_scales_CV + 1); // !!! ???
	      distribution_group_run_result_CV[i_m][i_z][i_d][i_b][i_s] = atof(readin[2 + i_b * 3].c_str());
	      distribution_group_run_deviation_CV[i_m][i_z][i_d][i_b][i_s] = atof(readin[3 + i_b * 3].c_str());
	      
	      // only to deal with incompletely filled distribution files (N->0)
	      if (distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b] == 0 && distribution_group_run_result_CV[i_m][i_z][i_d][i_b][i_s] != 0.){
		logger << LOG_INFO << name << "   distribution_group_run_N_binwise_CV corrected." << endl;
		distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b] = 1;
	      }
	      // by here.
	    }

	  }
	  else {
	    distribution_group_run_N_CV[i_m][i_z] = 0;
	    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	      distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b] = 0;
	      distribution_group_run_result_CV[i_m][i_z][i_d][i_b][i_s] = 0.;
	      distribution_group_run_deviation_CV[i_m][i_z][i_d][i_b][i_s] = 0.;
	    }
	  }
	}
      }

    }
  }

  for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
    for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
      distribution_group_N_CV[i_m] += distribution_group_run_N_CV[i_m][i_z];
      // with no runs removed at [i_d][i_b], one should later get: distribution_group_run_N_CV[i_m] == distribution_group_N_total_CV[i_m][i_d][i_b] !!!
    }
  }


  
  /* old result:  
  for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	string filename = "../" + ycontribution->directory[i_z] + "/distribution/CV/" + osi->directory_name_scale_CV[i_s] + "/" + osi->extended_distribution[i_d].xdistribution_name + "_" + name + ".dat";
	vector<string> readin;
	char LineBuffer[128];
	ifstream in_result(filename.c_str());  
	while (in_result.getline(LineBuffer, 128)){readin.push_back(LineBuffer);}
	in_result.close();
	logger << LOG_INFO << filename << "   number of lines: " << setw(5) << readin.size() << endl;
	if (readin.size() != 0){
	  distribution_run_N_CV[i_z][i_d] = atoll(readin[0].c_str());
	  for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	    distribution_run_N_binwise_CV[i_z][i_d][i_b] = atoll(readin[1 + i_b * 3].c_str());
	    //  N[i_z][i_p][i_d][i_b] = N[i_z][i_p][i_d][i_b] / (osi->n_scales_CV + 1); // !!! ???
	    distribution_run_result_CV[i_z][i_d][i_b][i_s] = atof(readin[2 + i_b * 3].c_str());
	    distribution_run_deviation_CV[i_z][i_d][i_b][i_s] = atof(readin[3 + i_b * 3].c_str());
	    //  distribution_run_deviation_CV[i_z][i_d][i_b][i_s] = distribution_run_deviation_CV[i_z][i_d][i_b][i_s] * 2; // !!! ???

	    if (distribution_run_N_binwise_CV[i_z][i_d][i_b] == 0 && distribution_run_result_CV[i_z][i_d][i_b][i_s] != 0.){
	      logger << LOG_INFO << name << "   distribution_run_N_binwise_CV corrected." << endl;
	      distribution_run_N_binwise_CV[i_z][i_d][i_b] = 1;
	    }
	    
	  }
	}
	else {
	  distribution_run_N_CV[i_z][i_d] = 0;
	  for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	    distribution_run_N_binwise_CV[i_z][i_d][i_b] = 0;
	    distribution_run_result_CV[i_z][i_d][i_b][i_s] = 0.;
	    distribution_run_deviation_CV[i_z][i_d][i_b][i_s] = 0.;
	  }
	}
	if (i_s == 0){
	  logger << LOG_INFO << name << "   distribution_run_CV[" << i_z << "][" << i_d << "][" << 0 << "][" << i_s << "] = " << distribution_run_result_CV[i_z][i_d][0][i_s] << endl;
	}
      }
    }
  }
*/

  // part of old implementation that should not be needed any more...
  /*
  //  distribution_N_TSV-equivalent filled here ???
  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    distribution_N_CV[i_d] = 0;
    for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
      distribution_N_CV[i_d] += distribution_run_N_CV[i_z][i_d];
    }
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
      distribution_N_binwise_CV[i_d][i_b] = 0;
      for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
	distribution_N_binwise_CV[i_d][i_b] += distribution_run_N_binwise_CV[i_z][i_d][i_b];
      }
    }
  }


  //  new: definition of
  //  all_distribution_N_binwise_CV
  //  combination_all_distribution_N_binwise_CV
  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    all_distribution_N_binwise_CV[i_d].resize(osi->extended_distribution[i_d].n_bins, 0);
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
      for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
	all_distribution_N_binwise_CV[i_d][i_b] += distribution_run_N_binwise_CV[i_z][i_d][i_b];
      }
      for (int x_z = 0; x_z < 2; x_z++){
	combination_all_distribution_N_binwise_CV[x_z][i_d].resize(osi->extended_distribution[i_d].n_bins, 0);
	for (int i_z = 0; i_z < ycontribution->directory.size(); i_z++){
	  if ((x_z == 0 && ycontribution->directory_extra[i_z] == 0) || (x_z == 1 && ycontribution->directory_extra[i_z] == 1)){
	    combination_all_distribution_N_binwise_CV[x_z][i_d][i_b] += distribution_run_N_binwise_CV[i_z][i_d][i_b];
	  }
	}
      }
    }
  }
*/
    
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combination_distribution_CV(){
  Logger logger("summary_subprocess::combination_distribution_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    logger << LOG_DEBUG << left << setw(20) << name << "  osi->extended_distribution[" << i_d << "].xdistribution_name = " << osi->extended_distribution[i_d].xdistribution_name << endl;
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){

      double tolerance_factor = ygeneric->deviation_tolerance_factor;

      // replace by reference value ???
      int x_s = (osi->n_scales_CV - 1) / 2;
      
      for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
	distribution_group_counter_nonzero_run_CV[i_m][i_d][i_b] = 0;
	vector<double> deviation_measure(ycontribution->extended_directory[i_m].size(), 0.);
	vector<double> deviation_measure_sorted;
	double min_deviation_measure = 1.e99;
	double max_deviation_measure = 0.;

	if (ycontribution->type_contribution == "RRA" && name == "du~_emmumepvm~gg"){
	  logger << LOG_DEBUG << name << endl;
	}
						      
	for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
	  /*
	  if (ycontribution->type_contribution == "VT2"){
	    if (distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b] > 0){
	      logger << LOG_DEBUG << "distribution_group_run_result_CV[" << setw(3) << i_m << "][" << setw(3) << i_z << "][" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << x_s << "] = "
		     << setw(15) << setprecision(8)
		     << distribution_group_run_result_CV[i_m][i_z][i_d][i_b][x_s] / distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b]
		     << " +- "
		     << setw(15) << setprecision(9) << sqrt(distribution_group_run_deviation_CV[i_m][i_z][i_d][i_b][x_s] - pow(distribution_group_run_result_CV[i_m][i_z][i_d][i_b][x_s], 2) / distribution_group_run_N_CV[i_m][i_z]) / (distribution_group_run_N_CV[i_m][i_z] - 1)
		     << "   "
		     << setw(8)
		     << distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b]
		     << endl;	// simplified version available ???
	    }
	  }
	  */
	  
	  distribution_group_N_binwise_CV[i_m][i_d][i_b] += distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b];
	  if (distribution_group_run_N_CV[i_m][i_z] > 0){
	    deviation_measure[i_z] = sqrt(distribution_group_run_deviation_CV[i_m][i_z][i_d][i_b][x_s] - pow(distribution_group_run_result_CV[i_m][i_z][i_d][i_b][x_s], 2) / distribution_group_run_N_CV[i_m][i_z]) / (distribution_group_run_N_CV[i_m][i_z] - 1) * sqrt(distribution_group_run_N_CV[i_m][i_z]);

	    if (deviation_measure[i_z] > 0.){
	      distribution_group_counter_nonzero_run_CV[i_m][i_d][i_b]++;
	      if (deviation_measure[i_z] < min_deviation_measure){min_deviation_measure = deviation_measure[i_z];}
	      if (deviation_measure[i_z] > max_deviation_measure){max_deviation_measure = deviation_measure[i_z];}
	      deviation_measure_sorted.push_back(deviation_measure[i_z]);
	    }
	    //	      else {logger << LOG_WARN << "deviation_measure[" << i_z << "] = " << deviation_measure[i_z] << endl;}
	  }
	}

	if (deviation_measure_sorted.size() > 0){
	  sort(deviation_measure_sorted.begin(), deviation_measure_sorted.end());
	  // in order to ignore runs with too low error estimates:
	  min_deviation_measure = deviation_measure_sorted[deviation_measure_sorted.size() / 10];	  


	  logger << LOG_DEBUG_VERBOSE << "distribution_group_N_binwise_CV[" << i_m << "][" << i_d << "][" << i_b << "] = " << distribution_group_N_binwise_CV[i_m][i_d][i_b] << endl;
	  
	  double temp_fraction = double(distribution_group_N_binwise_CV[i_m][i_d][i_b]) / distribution_group_N_CV[i_m];
	  tolerance_factor = ygeneric->deviation_tolerance_factor * (1. - log10(temp_fraction));
	  if (tolerance_factor < ygeneric->deviation_tolerance_factor){tolerance_factor = ygeneric->deviation_tolerance_factor;}

	  for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
	    if (deviation_measure[i_z] > tolerance_factor * min_deviation_measure &&
		distribution_group_N_binwise_CV[i_m][i_d][i_b] > 1000){
	      distribution_group_run_removal_run_CV[i_m][i_z][i_d][i_b] = 1;
	      distribution_group_counter_removal_run_CV[i_m][i_d][i_b]++;
	    }
	    else {
	      distribution_group_N_total_CV[i_m][i_d][i_b] += distribution_group_run_N_CV[i_m][i_z];
	      distribution_group_N_total_binwise_CV[i_m][i_d][i_b] += distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b];
	      if (distribution_group_N_binwise_CV[i_m][i_d][i_b] <= 1000){
		logger << LOG_DEBUG << setw(20) << ycontribution->infix_order_contribution << "   " << setw(20) << name << "TFE [" << i_m << "] -> " << setw(6) << ygeneric->phasespace_optimization[i_m] << ": [" << i_d << "] = " << setw(25) << osi->extended_distribution[i_d].xdistribution_name << "   [" << setw(3) << i_b << "] = " << osi->extended_distribution[i_d].bin_edge[i_b] << "   " << endl;

	      }
	    }
	  }
	}
      }
      

      for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
	distribution_N_total_CV[i_d][i_b] += distribution_group_N_total_CV[i_m][i_d][i_b];
	distribution_N_total_binwise_CV[i_d][i_b] += distribution_group_N_total_binwise_CV[i_m][i_d][i_b];
      }
 


      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	if (ycontribution->average_factor == 1){
	  combine_distribution_conservative_CV(i_d, i_b, i_s);
	}
	else {
	  //	  combine_distribution_hybrid_CV(i_d, i_b, i_s);
	}
      }

      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      //  weighted combination of results from different extended_directory contributions:
      //  should receive a mechanism to prevent underestimated errors from low-statistics contributions to get much weight !!!
      double result_subprocess = 0.;
      double deviation_subprocess = 0.;
      vector<double> deltasigma(ycontribution->extended_directory.size(), 0.);
      vector<double> singleweight(ycontribution->extended_directory.size(), 0.);
      for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
	logger << LOG_DEBUG_VERBOSE << "distribution_group_result_CV[" << setw(3) << i_m << "][" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << i_s << "] = " << setw(15) << setprecision(8) << distribution_group_result_CV[i_m][i_d][i_b][i_s] << " +- " << setw(15) << setprecision(9) << distribution_group_deviation_CV[i_m][i_d][i_b][i_s] << endl;

	//		if (distribution_group_N_binwise_CV[i_m][i_d][i_b][i_s] > 0){
	//	if (distribution_group_N_binwise_CV[i_m][i_d][i_b] > 0){
	if (distribution_group_N_binwise_CV[i_m][i_d][i_b] > 0 &&
	    distribution_group_deviation_CV[i_m][i_d][i_b][i_s] != 0. &&
	    distribution_group_N_total_binwise_CV[i_m][i_d][i_b] > distribution_N_total_binwise_CV[i_d][i_b] / 100){
	  ///	    distribution_group_N_total_CV[i_m][i_d][i_b] > distribution_N_total_CV[i_d][i_b] / 100){
	  //		singleweight[i_m] = 1. / distribution_group_deviation2_CV[i_m][i_d][i_b][i_s];
	  singleweight[i_m] = 1. / pow(distribution_group_deviation_CV[i_m][i_d][i_b][i_s], 2);
	  result_subprocess += singleweight[i_m] * distribution_group_result_CV[i_m][i_d][i_b][i_s];
	  deviation_subprocess += singleweight[i_m];
	}
	else if (distribution_group_N_total_binwise_CV[i_m][i_d][i_b] < distribution_N_total_binwise_CV[i_d][i_b] / 100){
	  ///	else if (distribution_group_N_total_CV[i_m][i_d][i_b] <= distribution_N_total_CV[i_d][i_b] / 100){
	  logger << LOG_DEBUG << setw(20) << ycontribution->infix_order_contribution << "   " << setw(20) << name << "OCU [" << i_m << "] -> " << setw(6) << ygeneric->phasespace_optimization[i_m] << ": [" << i_d << "] = " << setw(25) << osi->extended_distribution[i_d].xdistribution_name << "   [" << setw(3) << i_b << "] = " << osi->extended_distribution[i_d].bin_edge[i_b] << "   " << endl;

	}
      }
      if (deviation_subprocess != 0.){
	distribution_result_CV[i_d][i_b][i_s] = result_subprocess / deviation_subprocess;
	/// simplified version that does not take into account the difference of central values in determination of new standard deviation...
	distribution_deviation_CV[i_d][i_b][i_s] = 1. / sqrt(deviation_subprocess);
	/*
	  double final_deviation = 0.;
	  for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
	  final_deviation += (singleweight[i_m] + 4 * pow(singleweight[i_m] * (distribution_group_result_CV[i_m][i_d][i_b][i_s] - distribution_result_CV[i_d][i_b][i_s]), 2)) / pow(deviation_subprocess, 2);
	  }	   
	  distribution_deviation_CV[i_d][i_b][i_s] = sqrt(final_deviation);
	*/
      }
      else {
	distribution_result_CV[i_d][i_b][i_s] = 0.;
	distribution_deviation_CV[i_d][i_b][i_s] = 0.;
      }


      //	if (ycontribution->type_contribution == "RRA" && name == "du~_emmumepvm~gg"){
      if (i_d == 0){
	logger << LOG_DEBUG
	       << setw(5) << ycontribution->type_contribution
	       << "   "
	       << setw(20) << name
	       << "   [i_d][i_b][i_s]"
	       << "[" << setw(3) << i_d << "]"
	       << "[" << setw(3) << i_d << "]"
	       << "[" << setw(3) << i_b << "]"
	       << "[" << setw(3) << i_s << "]"
	       << "   "
	       << setw(23) << setprecision(15) << distribution_result_CV[i_d][i_b][i_s]
	       << " +- "
	       << setw(23) << setprecision(15) << distribution_deviation_CV[i_d][i_b][i_s]
	       << endl;
      }

      
      logger << LOG_DEBUG_VERBOSE << "           distribution_result_CV[" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << i_s << "] = " << setw(15) << setprecision(8) << distribution_result_CV[i_d][i_b][i_s] << " +- " << setw(15) << setprecision(9) << distribution_deviation_CV[i_d][i_b][i_s] << endl;
	    
      //  chi2 calculation from [i_m] combination:
      distribution_chi2_CV[i_d][i_b][i_s] = 0.;
      logger << LOG_DEBUG_VERBOSE << "chi2 evaluation in " << ycontribution->resultdirectory << " -- " << ycontribution->infix_contribution << " -- " << ycontribution->infix_order_contribution <<  " -- " << ycontribution->infix_path_contribution << " --- " << name << endl;
      logger << LOG_DEBUG_VERBOSE << "combined result = " << setw(15) << setprecision(9) << distribution_result_CV[i_d][i_b][i_s] << " +- " << setw(15) << setprecision(8) << distribution_deviation_CV[i_d][i_b][i_s] << endl;
      
      int temp_counter = 0;
      for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
	if (distribution_group_N_CV[i_m] > 0){
	  if (!(distribution_group_result_CV[i_m][i_d][i_b][i_s] == 0. && distribution_group_deviation_CV[i_m][i_d][i_b][i_s] == 0.)){
	    temp_counter++;
	    distribution_chi2_CV[i_d][i_b][i_s] += pow((distribution_group_result_CV[i_m][i_d][i_b][i_s] - distribution_result_CV[i_d][i_b][i_s]) / distribution_group_deviation_CV[i_m][i_d][i_b][i_s], 2.);
	  }
	}
      }
      if (temp_counter > 1){
	distribution_chi2_CV[i_d][i_b][i_s] = distribution_chi2_CV[i_d][i_b][i_s] / temp_counter;
      }
      else {
	//  chi2 cannot be reasonably defined for only one run.
	distribution_chi2_CV[i_d][i_b][i_s] = 0.;
      }
      logger << LOG_DEBUG_VERBOSE << "distribution_chi2_CV[" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(3) << i_s << "] = " << distribution_chi2_CV[i_d][i_b][i_s] << " (" << temp_counter << ")" << endl;
      }
      
    }
  }

  // new:
  // some more such clean-ups to avoid too much memory consumption !!!

  vector<vector<long long> > ().swap(distribution_group_run_N_CV);
  // new: clear ->
  vector<vector<vector<vector<int> > > > ().swap(distribution_group_run_removal_run_CV);
  // <-
  vector<vector<vector<vector<long long> > > > ().swap(distribution_group_run_N_binwise_CV);
  vector<vector<vector<vector<vector<double> > > > > ().swap(distribution_group_run_result_CV);
  vector<vector<vector<vector<vector<double> > > > > ().swap(distribution_group_run_deviation_CV);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_subprocess::combine_distribution_conservative_CV(int i_d, int i_b, int i_s){
  Logger logger("summary_subprocess::combine_distribution_conservative_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // !!! just for now: x_z is not needed as a parameter...

  for (int i_m = 0; i_m < ycontribution->extended_directory.size(); i_m++){
    double temp_result_subprocess = 0.;
    double temp_deviation2_subprocess = 0.;
    for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
      if (!distribution_group_run_removal_run_CV[i_m][i_z][i_d][i_b] && distribution_group_run_N_binwise_CV[i_m][i_z][i_d][i_b] > 0){
	temp_result_subprocess += distribution_group_run_result_CV[i_m][i_z][i_d][i_b][i_s];
	  // should be named sum_result...
	temp_deviation2_subprocess += distribution_group_run_deviation_CV[i_m][i_z][i_d][i_b][i_s];
	  // should be named sum_deviation²...
      }
    }
    if (distribution_group_N_total_CV[i_m][i_d][i_b] != 0){
      //    if (distribution_group_N_binwise_CV[i_m][i_d][i_b] != 0){
      ///    if (distribution_group_N_binwise_CV[i_m][i_d][i_b][i_s] != 0){
      distribution_group_result_CV[i_m][i_d][i_b][i_s] = temp_result_subprocess / distribution_group_N_total_CV[i_m][i_d][i_b];
      // distribution_group_N_total_CV[i_m][i_d][i_b]  instead of  distribution_group_N_CV[i_m]  since the latter contains also results from removed runs
      //  distribution_group_deviation2_CV[i_m][i_d][i_b][i_s] = (distribution_group_N_CV[i_m] * temp_deviation2_subprocess - pow(temp_result_subprocess, 2)) / distribution_group_N_CV[i_m] / pow(distribution_group_N_CV[i_m] - 1, 2); // possible replacement - check if correct !!!
      distribution_group_deviation_CV[i_m][i_d][i_b][i_s] = sqrt((distribution_group_N_total_CV[i_m][i_d][i_b] * temp_deviation2_subprocess - pow(temp_result_subprocess, 2)) / distribution_group_N_total_CV[i_m][i_d][i_b]) / (distribution_group_N_total_CV[i_m][i_d][i_b] - 1);
    }
    else {
      distribution_group_result_CV[i_m][i_d][i_b][i_s] = 0.;
      distribution_group_deviation_CV[i_m][i_d][i_b][i_s] = 0.;
    }

    if (i_d == 0){
      logger << LOG_DEBUG
	     << setw(5) << ycontribution->type_contribution
	     << "   "
	     << setw(20) << name
	     << "   group"
	     << "[" << setw(3) << i_m << "]"
	     << "[" << setw(3) << i_d << "]"
	     << "[" << setw(3) << i_b << "]"
	     << "[" << setw(3) << i_s << "]"
	     << "   "
	     << setw(23) << setprecision(15) << distribution_group_result_CV[i_m][i_d][i_b][i_s]
	     << " +- "
	     << setw(23) << setprecision(15) << distribution_group_deviation_CV[i_m][i_d][i_b][i_s]
	     << endl;
    }


    //  chi2 calculation:
    distribution_group_chi2_CV[i_m][i_d][i_b][i_s] = 0.;
    logger << LOG_DEBUG_VERBOSE << "chi2 evaluation in " << ycontribution->resultdirectory << " -- " << ycontribution->infix_contribution << " -- " << ycontribution->infix_order_contribution <<  " -- " << ycontribution->infix_path_contribution << " --- " << name << endl;
    logger << LOG_DEBUG_VERBOSE << "combined result = " << setw(15) << setprecision(9) << distribution_group_result_CV[i_m][i_d][i_b][i_s] << " +- " << setw(15) << setprecision(9) << distribution_group_deviation_CV[i_m][i_d][i_b][i_s] << endl;

    if (distribution_group_counter_nonzero_run_CV[i_m][i_d][i_b] > 1){
      for (int i_z = 0; i_z < ycontribution->extended_directory[i_m].size(); i_z++){
	if (!(distribution_group_run_result_CV[i_m][i_z][i_d][i_b][i_s] == 0. && distribution_group_run_deviation_CV[i_m][i_z][i_d][i_b][i_s] == 0.)){
	  // expectation value -> temp_run_result
	  // sum (result) ->distribution_group_run_result_CV[i_m][i_z][i_d][i_b][i_s]
	  double temp_run_result = distribution_group_run_result_CV[i_m][i_z][i_d][i_b][i_s] / distribution_group_run_N_CV[i_m][i_z]; 
	  // standard deviation -> temp_run_deviation
	  // sum (result²) ->distribution_group_run_deviation_CV[i_m][i_z][i_d][i_b][i_s]
	  double temp_run_deviation = sqrt(distribution_group_run_deviation_CV[i_m][i_z][i_d][i_b][i_s] - pow(distribution_group_run_result_CV[i_m][i_z][i_d][i_b][i_s], 2) / distribution_group_run_N_CV[i_m][i_z]) / (distribution_group_run_N_CV[i_m][i_z] - 1);
	  distribution_group_chi2_CV[i_m][i_d][i_b][i_s] += pow((temp_run_result - distribution_group_result_CV[i_m][i_d][i_b][i_s]) / temp_run_deviation, 2.);
	}
      }
      distribution_group_chi2_CV[i_m][i_d][i_b][i_s] = distribution_group_chi2_CV[i_m][i_d][i_b][i_s] / distribution_group_counter_nonzero_run_CV[i_m][i_d][i_b];
    }
    else {
      // chi2 cannot be reasonably defined for only one run.
      distribution_group_chi2_CV[i_m][i_d][i_b][i_s] = 0.;
    }
  }
    
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


