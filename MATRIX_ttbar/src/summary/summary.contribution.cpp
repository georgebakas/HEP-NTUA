#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"
////////////////////
//  constructors  //
////////////////////
summary_contribution::summary_contribution(){
  /*
 // each 'summary_contribution' belongs to a specific 'summary_list'
  ylist = &_list;
  // only one 'summary_generic'
  ygeneric = ylist->ygeneric;
  // only one 'observable_set'
  osi = &ygeneric->oset;
  */
  /*
  subprocess.resize(1);
  xsubprocess.resize(1);
  subgroup_no_member.resize(1);

  in_contribution_order_alpha_s = 0;
  in_contribution_order_alpha_e = 0;
  interference = 0;
  photon_induced = 0;
  */
}

summary_contribution::summary_contribution(string _type_contribution, summary_list & _list){
  Logger logger("summary_contribution::summary_contribution");
  logger << LOG_DEBUG << "called" << endl;

  // each 'summary_contribution' belongs to a specific 'summary_list'
  ylist = &_list;
  // only one 'summary_generic'
  ygeneric = ylist->ygeneric;
  // only one 'observable_set'
  osi = &ygeneric->oset;

  logger << LOG_INFO << "ylist->resultdirectory = " << ylist->resultdirectory << endl;

  type_contribution = _type_contribution;
  subprocess.resize(1);
  xsubprocess.resize(1);
  subgroup_no_member.resize(1);

  in_contribution_order_alpha_s = 0;
  in_contribution_order_alpha_e = 0;
  interference = 0;
  photon_induced = 0;

  logger << LOG_DEBUG << "finished" << endl;
}


void summary_contribution::readin_contribution_remove_run(){
  Logger logger("summary_contribution::readin_contribution_remove_run");

  int n_seed = directory.size();
  int default_size_CV = 1 + 2;
  if (active_qTcut){default_size_CV += 2 * osi->n_qTcut * osi->n_scales_CV;}
  else {default_size_CV += 2 * osi->n_scales_CV;}

  xsubprocess.resize(subprocess.size());
  for (int i_p = 0; i_p < xsubprocess.size(); i_p++){
    xsubprocess[i_p] = summary_subprocess(subprocess[i_p], n_seed, default_size_CV, *this);
    /*
    xsubprocess[i_p].average_factor = average_factor;
    xsubprocess[i_p].active_qTcut = active_qTcut;
    xsubprocess[i_p].output_n_qTcut = output_n_qTcut;
    xsubprocess[i_p].selection_n_qTcut = selection_n_qTcut;
    */
  }
}


void summary_contribution::readin_runtime_contribution(){
  Logger logger("summary_contribution::readin_runtime_contribution");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int n_seed = directory.size();
  logger << LOG_DEBUG_VERBOSE << "XXX   n_seed = " << n_seed << endl;

  int default_size_CV =  1 + 2;
  if (active_qTcut){default_size_CV += 2 * osi->n_qTcut * osi->n_scales_CV;}
  else {default_size_CV += 2 * osi->n_scales_CV;}

  logger << LOG_DEBUG_VERBOSE << "default_size_CV = " << default_size_CV << endl;
  /*
  xsubprocess.resize(subprocess.size());
  for (int i_p = 0; i_p < xsubprocess.size(); i_p++){
    xsubprocess[i_p] = summary_subprocess(subprocess[i_p], n_seed, default_size_CV, *this);
    xsubprocess[i_p].average_factor = average_factor;
    xsubprocess[i_p].active_qTcut = active_qTcut;
    xsubprocess[i_p].output_n_qTcut = output_n_qTcut;
    xsubprocess[i_p].selection_n_qTcut = selection_n_qTcut;
  }
  */

  for (int i_p = 1; i_p < subprocess.size(); i_p++){
    string temp = "result";
    xsubprocess[i_p].readin_runtime(temp);
    logger << LOG_DEBUG << "xsubprocess[" << i_p << "].name = " << xsubprocess[i_p].name << "   in1 = " << xsubprocess[i_p].type_parton[0][1] << "   in2 = " << xsubprocess[i_p].type_parton[0][2] << endl;
  }

  // reset to zero for [0] component !!!
  xsubprocess[0].error2_time = 0.;
  xsubprocess[0].used_runtime = 0.;
  xsubprocess[0].used_n_event = 0;

  for (int i_p = 1; i_p < subprocess.size(); i_p++){
    logger << LOG_INFO << "CROSSCHECK   " << xsubprocess[i_p].name << "   xsubprocess[" << i_p << "].error2_time = " << xsubprocess[i_p].error2_time << endl;
    xsubprocess[0].error2_time += xsubprocess[i_p].error2_time;
    xsubprocess[0].used_runtime += xsubprocess[i_p].used_runtime;
    xsubprocess[0].used_n_event += xsubprocess[i_p].used_n_event;
    //    cout << "used_runtime[" << i_p << "] = " << xsubprocess[i_p].used_runtime << endl;
  }
  cout << "xsubprocess_0_used_runtime = " << xsubprocess[0].used_runtime << endl;

  logger << LOG_INFO << "SUMRUNTIME   " << type_contribution << "   " << setw(15) << xsubprocess[0].name << "   xsubprocess[0].error2_time = " << xsubprocess[0].error2_time << "   xsubprocess[0].used_runtime = " << xsubprocess[0].used_runtime << endl;

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_contribution::output_extrapolated_runtime_directories(int i_m, int & list_contribution_counter){
  Logger logger("summary_contribution:output_extrapolated_runtime_directories");
  logger << LOG_DEBUG << "called" << endl;

  //////////////////////////////////////////////////
  //  extrapolated runtime - new run directories  //
  //////////////////////////////////////////////////

  /*  int n_run_phasespace_optimization = ygeneric->phasespace_optimization.size();
  if (ygeneric->phasespace_optimization.size() > 1){n_run_phasespace_optimization = ygeneric->phasespace_optimization.size() - 1;}
  for (int i_m = 0; i_m < n_run_phasespace_optimization; i_m++){
  */
  
  // actually needed only once:
  stringstream list_path;
  list_path << ygeneric->final_resultdirectory + "/newrundir";
  list_path << "/" << ylist->type_perturbative_order;
  if (ylist->type_subtraction_method != "---"){list_path << "." << ylist->type_subtraction_method;}
  //  list_path << ygeneric->phasespace_optimization[i_m];
  system_execute(logger, "mkdir " + list_path.str());
  logger << LOG_INFO << "list_path = " << list_path.str() << endl;
  
  stringstream temp_photon_induced;
  if (ylist->photon_induced){temp_photon_induced << "a";}
  list_path << "/" << temp_photon_induced.str() << ylist->in_contribution_order_alpha_s << ylist->in_contribution_order_alpha_e;
  system_execute(logger, "mkdir " + list_path.str());
  logger << LOG_INFO << "list_path = " << list_path.str() << endl;
    
  /*
  stringstream list_path_list_start;
  list_path_list_start << ygeneric->final_resultdirectory + "/newrundir/start/list.run.subprocesses." + ylist->xcontribution[0].infix_order_contribution + ".txt";
  ofstream out_runlist_new;
  logger << LOG_INFO << "list_path_list_start.str() = " << list_path_list_start.str() << endl;
  out_runlist_new.open((list_path_list_start.str()).c_str(), ofstream::out | ofstream::app);  
  out_runlist_new << "./start/MUNICH." + infix_order_contribution + ".start.sh" << endl;
  out_runlist_new.close();
  */

  ofstream out_runscript_contribution;
  string runscript_name = ygeneric->final_resultdirectory + "/newrundir/start" + ygeneric->phasespace_optimization[i_m] + "/MUNICH." + infix_order_contribution + ".start.sh";
  out_runscript_contribution.open(runscript_name.c_str(), ofstream::out | ofstream::trunc);  
  out_runscript_contribution << "echo '" << infix_path_contribution << ygeneric->phasespace_optimization[i_m] << "'" << endl;
  out_runscript_contribution << "cd " << infix_path_contribution << ygeneric->phasespace_optimization[i_m] << endl;
  
  ofstream out_runlist_contribution;
  string runlist_name = ygeneric->final_resultdirectory + "/newrundir/start" + ygeneric->phasespace_optimization[i_m] + "/list.run.subprocesses." + infix_order_contribution + ".txt";
  out_runlist_contribution.open(runlist_name.c_str(), ofstream::out | ofstream::trunc);  
  
  ofstream out_plainrunlist_contribution;
  string plainrunlist_name = ygeneric->final_resultdirectory + "/newrundir/start" + ygeneric->phasespace_optimization[i_m] + "/runlist." + infix_order_contribution + ".dat";
  out_plainrunlist_contribution.open(plainrunlist_name.c_str(), ofstream::out | ofstream::trunc);  
  
  // could be done only once for the contribution (as a member variable of contribution) !!!
  //  vector<int> xsubprocess[i_p].max_number_of_jobs_for_one_channel(xsubprocess.size(), ygeneric->run_min_number_of_jobs_for_one_channel); 
  
  stringstream temp_path;
  temp_path << list_path.str() << "/" << infix_contribution << ygeneric->phasespace_optimization[i_m];
  system_execute(logger, "mkdir " + temp_path.str());
  //  int 
  max_number_of_jobs_for_one_channel = ygeneric->run_min_number_of_jobs_for_one_channel;
  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    xsubprocess[i_p].max_number_of_jobs_for_one_channel = ygeneric->run_min_number_of_jobs_for_one_channel;
    double temp_njc = xsubprocess[i_p].extrapolated_runtime / ygeneric->run_max_time_per_job;
    if (temp_njc > ygeneric->run_min_number_of_jobs_for_one_channel){xsubprocess[i_p].max_number_of_jobs_for_one_channel = int(temp_njc) + 1;}
    else {xsubprocess[i_p].max_number_of_jobs_for_one_channel = ygeneric->run_min_number_of_jobs_for_one_channel;}
    if (xsubprocess[i_p].max_number_of_jobs_for_one_channel > max_number_of_jobs_for_one_channel){max_number_of_jobs_for_one_channel = xsubprocess[i_p].max_number_of_jobs_for_one_channel;}
  }
  
  //  run.<xxx>/file_parameter.dat
  system_execute(logger, "mkdir " + temp_path.str() + "/log");
  for (int i_r = 0; i_r < max_number_of_jobs_for_one_channel; i_r++){
    stringstream thisrun_path;
    thisrun_path << temp_path.str() << "/run." << i_r;
    system_execute(logger, "mkdir " + thisrun_path.str());
    system_execute(logger, "mkdir " + thisrun_path.str() + "/log");
    ofstream out_file_parameter_run;
    string name_file_parameter_run = thisrun_path.str() + "/file_parameter.dat";
    out_file_parameter_run.open(name_file_parameter_run.c_str(), ofstream::out | ofstream::trunc);  
    write_infile_int(out_file_parameter_run, "zwahl", i_r);
    out_file_parameter_run.close();
  }
  
  //  <runscript> and <processlist> of contribution
  //  int list_contribution_counter = 0;
  stringstream ss_runlist_contribution;
  for (int i_r = 0; i_r < max_number_of_jobs_for_one_channel; i_r++){
    out_runscript_contribution << "echo 'run." << i_r << "'" << endl;
    out_runscript_contribution << "cd run." << i_r << endl;
    out_runscript_contribution << "../../../../MUNICH.start.sh" << endl;
    out_runscript_contribution << "cd .." << endl;
    
    for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
      if (xsubprocess[i_p].extrapolated_sigma_normalization == 0.){continue;}
      if (i_r < xsubprocess[i_p].max_number_of_jobs_for_one_channel){
	ss_runlist_contribution << right << setw(5) << ++list_contribution_counter << "   " << left << setw(25) << infix_path_contribution + ygeneric->phasespace_optimization[i_m] << "   run." << left << setw(4) << i_r << "   " << setw(20) << xsubprocess[i_p].name << endl;
      }
    }
  }
  out_runscript_contribution << "cd ../../../" << endl;
  out_runscript_contribution << "echo" << endl;
  out_runscript_contribution << endl;
  out_runscript_contribution.close();
  
  system_execute(logger, "chmod 700 " + runscript_name);
  
  out_runlist_contribution << right << setw(8) << "" << setw(30) << left << infix_order_contribution + ygeneric->phasespace_optimization[i_m] << "     " << right << setw(5) << list_contribution_counter << " ( " << setw(5) << max_number_of_jobs_for_one_channel << " )" << endl;
  out_runlist_contribution << right << setw(8) << "========" << setw(30) << "==============================" << "=====" << "=====" << "===" << "=====" << "==" << endl;
  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    out_runlist_contribution << right << setw(8) << "" << setw(30) << left << xsubprocess[i_p].name << "     " << right << setw(5) << xsubprocess[i_p].max_number_of_jobs_for_one_channel << " ( " << setw(5) << xsubprocess[i_p].max_number_of_jobs_for_one_channel << " )" << endl;
  }
  out_runlist_contribution << endl;
  out_runlist_contribution << right << setw(5) << "no" << "   " << left << setw(30) << "path to contribution" << "   " << "run " << left << setw(5) << "" << "  " << setw(20) << "subprocess" << endl;
  out_runlist_contribution << right << setw(8) << "========" << setw(30) << "==============================" << "===" << "=====" << "=====" << "==" << "====================" << endl;
  out_runlist_contribution << ss_runlist_contribution.str();
  ss_runlist_contribution.clear();
  out_runlist_contribution.close();
  
  for (int i_r = 0; i_r < max_number_of_jobs_for_one_channel; i_r++){
    out_plainrunlist_contribution << infix_path_contribution << ygeneric->phasespace_optimization[i_m] << "/run." << i_r << endl;
  }
  out_plainrunlist_contribution.close();
  
  //  run.<xxx>/subprocesslist.dat
  ofstream out_subprocesslist;
  for (int i_r = 0; i_r < max_number_of_jobs_for_one_channel; i_r++){
    stringstream thisrun_path;
    thisrun_path << temp_path.str() << "/run." << i_r;
    string subprocesslist_name = thisrun_path.str() + "/subprocesslist.dat";
    out_subprocesslist.open(subprocesslist_name.c_str(), ofstream::out | ofstream::trunc);  
    for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
      if (xsubprocess[i_p].extrapolated_sigma_normalization == 0.){continue;}
      if (i_r < xsubprocess[i_p].max_number_of_jobs_for_one_channel){
	out_subprocesslist << xsubprocess[i_p].name << endl;
      }
    }
    out_subprocesslist.close();
  }
    

  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    long long n_event_per_job = (long long)(xsubprocess[i_p].extrapolated_n_event) / xsubprocess[i_p].max_number_of_jobs_for_one_channel + 1;
    double sigma_normalization_deviation_per_job = xsubprocess[i_p].extrapolated_sigma_normalization_deviation * sqrt(xsubprocess[i_p].max_number_of_jobs_for_one_channel);
    for (int i_r = -1; i_r < xsubprocess[i_p].max_number_of_jobs_for_one_channel; i_r++){
      //    for (int i_r = 0; i_r < xsubprocess[i_p].max_number_of_jobs_for_one_channel; i_r++){
      stringstream thisrun_path;
      if (i_r == -1){ thisrun_path << temp_path.str() << "/log";}
      else {thisrun_path << temp_path.str() << "/run." << i_r << "/log";}
      
      if (xsubprocess[i_p].extrapolated_sigma_normalization == 0.){continue;}
      
      //  run.<xxx>/log/file_parameter.<process>.dat
      ofstream out_file_parameter_process;
      string name_file_parameter_process = thisrun_path.str() + "/file_parameter." + xsubprocess[i_p].name + ".dat";
      out_file_parameter_process.open(name_file_parameter_process.c_str(), ofstream::out | ofstream::trunc);  
      
      stringstream temp_ext_runtime;
      int length_sec = int(log10(xsubprocess[i_p].extrapolated_runtime));
      temp_ext_runtime << setw(4 + length_sec) << showpoint << xsubprocess[i_p].extrapolated_runtime;
      out_file_parameter_process << "# runtime = " << temp_ext_runtime.str() << endl;
      
      long long rounded_n = n_event_per_job + ygeneric->run_min_n_step - n_event_per_job % ygeneric->run_min_n_step;
      long long selected_n = rounded_n;
      if (type_contribution != "VT2" &&
	  type_contribution != "RVA" &&
	  type_contribution != "L2VA" &&
	  type_contribution != "L2VT"){
	if (ygeneric->run_min_n_event > rounded_n){selected_n = ygeneric->run_min_n_event;}
      }
      
      write_infile_long_long(out_file_parameter_process, "n_step", selected_n / ygeneric->run_min_n_step);
      write_infile_long_long(out_file_parameter_process, "n_events_min", selected_n);
      write_infile_long_long(out_file_parameter_process, "n_events_max", ygeneric->run_factor_max * selected_n);
      write_infile_double(out_file_parameter_process, "sigma_normalization", xsubprocess[i_p].extrapolated_sigma_normalization);
      write_infile_double(out_file_parameter_process, "sigma_normalization_deviation", sigma_normalization_deviation_per_job);
      
      out_file_parameter_process.close();
      
      //      counter_subprocess++;
      //      counter_subprocess_complete++;
      //      out_runlist_complete << right << setw(5) << counter_subprocess_complete << "   " << left << setw(25) << infix_path_contribution << ygeneric->phasespace_optimization[i_m] << "   run." << left << setw(3) << i_r << "   " << xsubprocess[i_p].name << endl;
    }
  }
    
    
  //  }

  // introduce universal file_parameter.<subprocess>.dat files to avoid repeated copies of it !!! xxxxx
  /*
  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    long long n_event_per_job = (long long)(xsubprocess[i_p].extrapolated_n_event) / xsubprocess[i_p].max_number_of_jobs_for_one_channel + 1;
    double sigma_normalization_deviation_per_job = xsubprocess[i_p].extrapolated_sigma_normalization_deviation * sqrt(xsubprocess[i_p].max_number_of_jobs_for_one_channel);
    //    for (int i_r = 0; i_r < xsubprocess[i_p].max_number_of_jobs_for_one_channel; i_r++){
      stringstream thisrun_path;
      thisrun_path << temp_path.str() << "/log";
      
      if (xsubprocess[i_p].extrapolated_sigma_normalization == 0.){continue;}

      //  run.<xxx>/log/file_parameter.<process>.dat
      ofstream out_file_parameter_process;
      string name_file_parameter_process = thisrun_path.str() + "/file_parameter." + xsubprocess[i_p].name + ".dat";
      out_file_parameter_process.open(name_file_parameter_process.c_str(), ofstream::out | ofstream::trunc);  

      stringstream temp_ext_runtime;
      int length_sec = int(log10(xsubprocess[i_p].extrapolated_runtime));
      temp_ext_runtime << setw(4 + length_sec) << showpoint << xsubprocess[i_p].extrapolated_runtime;
      out_file_parameter_process << "# runtime = " << temp_ext_runtime.str() << endl;

      long long rounded_n = n_event_per_job + ygeneric->run_min_n_step - n_event_per_job % ygeneric->run_min_n_step;
      long long selected_n = rounded_n;
      if (type_contribution != "VT2" &&
	  type_contribution != "RVA" &&
	  type_contribution != "L2VA" &&
	  type_contribution != "L2VT"){
	if (ygeneric->run_min_n_event > rounded_n){selected_n = ygeneric->run_min_n_event;}
      }
	
      write_infile_long_long(out_file_parameter_process, "n_step", selected_n / ygeneric->run_min_n_step);
      write_infile_long_long(out_file_parameter_process, "n_events_min", selected_n);
      write_infile_long_long(out_file_parameter_process, "n_events_max", ygeneric->run_factor_max * selected_n);
      write_infile_double(out_file_parameter_process, "sigma_normalization", xsubprocess[i_p].extrapolated_sigma_normalization);
      write_infile_double(out_file_parameter_process, "sigma_normalization_deviation", sigma_normalization_deviation_per_job);
      
      out_file_parameter_process.close();

      //      counter_subprocess++;
      //      counter_subprocess_complete++;
      //      out_runlist_complete << right << setw(5) << counter_subprocess_complete << "   " << left << setw(25) << infix_path_contribution << ygeneric->phasespace_optimization[i_m] << "   run." << left << setw(3) << i_r << "   " << xsubprocess[i_p].name << endl;
      //    }
  }
  */
  
  
  //  out_runlist_complete << right << setw(5) << "" << setw(20) << infix_order_contribution << "   " << list_contribution_counter << endl;
  //  out_runscript_list << "./start/MUNICH." + infix_order_contribution + ".start.sh" << endl;
  
  logger << LOG_DEBUG << "finished" << endl;
}




