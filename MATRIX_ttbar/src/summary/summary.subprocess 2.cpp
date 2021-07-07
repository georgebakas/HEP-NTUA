#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

////////////////////
//  constructors  //
////////////////////
summary_subprocess::summary_subprocess(){}
summary_subprocess::summary_subprocess(string _name, summary_contribution & _contribution){

  ycontribution = &_contribution;
  ygeneric = ycontribution->ygeneric;
  osi = &ygeneric->oset;

  name = _name;
}

//#include "old.summary.subprocess.cpp"

summary_subprocess::summary_subprocess(string _name, int _n_seed, int _default_size_CV, summary_contribution & _contribution){
  Logger logger("summary_subprocess::summary_subprocess");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // each 'summary_subprocess' belongs to a specific 'summary_contribution'
  ycontribution = &_contribution;
  // only one 'summary_generic'
  ygeneric = ycontribution->ygeneric;
  // only one 'observable_set'
  osi = &ygeneric->oset;

  name = _name;
  n_seed = _n_seed;
  default_size_CV = _default_size_CV;



  no_results = 0;
  N_all = 0;

  result = 0.;
  deviation = 0.;
  chi2 = 0.;
  error2_time = 0.;
  time_per_event = 0.;
  used_runtime = 0.;
  used_n_event = 0;


  logger << LOG_INFO << "initialization:  " << setw(20) << name << "  in  " << setw(15) << ycontribution->type_contribution << endl;

  //  remove_run.resize(n_seed, 0);
  remove_run_qTcut_CV.resize(osi->n_qTcut, vector<int> (n_seed, 0));

  N.resize(n_seed, 0);

  result_run.resize(n_seed, 0.);
  deviation_run.resize(n_seed, 0.);

  error2_time_run.resize(n_seed, 0.);
  time_per_event_run.resize(n_seed, 0.);
  used_runtime_run.resize(n_seed, 0.);
  used_n_event_run.resize(n_seed, 0);

  N_event_run_CV.resize(n_seed, 0);
  N_event_CV = 0;
  Nd_event_CV = 0.;
  N_valid_event_CV.resize(osi->n_qTcut, 0);

  n_parallel_runs_reference_CV.resize(osi->n_qTcut, 0);
  n_valid_parallel_runs_reference_CV.resize(osi->n_qTcut, 0);

  result_CV.resize(osi->n_qTcut, vector<double> (osi->n_scales_CV, 0.));
  deviation_CV.resize(osi->n_qTcut, vector<double> (osi->n_scales_CV, 0.));
  chi2_CV.resize(osi->n_qTcut, vector<double> (osi->n_scales_CV, 0.));
  n_parallel_runs_CV.resize(osi->n_qTcut, vector<int> (osi->n_scales_CV, 0));
  n_valid_parallel_runs_CV.resize(osi->n_qTcut, vector<int> (osi->n_scales_CV, 0));

  result_run_CV.resize(ycontribution->output_n_qTcut, vector<vector<double> > (osi->n_scales_CV, vector<double> (n_seed, 0.)));
  deviation_run_CV.resize(ycontribution->output_n_qTcut, vector<vector<double> > (osi->n_scales_CV, vector<double> (n_seed, 0.)));
  //  result_run_CV.resize(osi->n_qTcut, vector<vector<double> > (osi->n_scales_CV, vector<double> (n_seed, 0.)));
  //  deviation_run_CV.resize(osi->n_qTcut, vector<vector<double> > (osi->n_scales_CV, vector<double> (n_seed, 0.)));

  if (osi->switch_TSV){
    //    no_results_TSV = 0;
    N_event_run_TSV.resize(n_seed, 0);
    N_event_TSV = 0;
    Nd_event_TSV = 0.;
    N_valid_event_TSV.resize(osi->n_qTcut, 0);

    n_parallel_runs_reference_TSV.resize(osi->n_qTcut, 0);
    n_valid_parallel_runs_reference_TSV.resize(osi->n_qTcut, 0);

    chi2_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<double> > > > (osi->n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV)));
    result_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<double> > > > (osi->n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV)));
    deviation_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<double> > > > (osi->n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV)));
    n_parallel_runs_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<int> > > > (osi->n_qTcut, vector<vector<vector<int> > > (osi->n_extended_set_TSV)));
    n_valid_parallel_runs_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<int> > > > (osi->n_qTcut, vector<vector<vector<int> > > (osi->n_extended_set_TSV)));
    
    for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
      for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  result_TSV[i_m][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  deviation_TSV[i_m][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  chi2_TSV[i_m][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  n_parallel_runs_TSV[i_m][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<int> (osi->n_scale_fact_TSV[i_s], 0));
	  n_valid_parallel_runs_TSV[i_m][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<int> (osi->n_scale_fact_TSV[i_s], 0));
	}
      }
    }
    
    logger << LOG_DEBUG << "ycontribution->output_n_qTcut < osi->n_qTcut = " << ycontribution->output_n_qTcut << " < " << osi->n_qTcut << endl;

    remove_run_qTcut_TSV.resize(ycontribution->output_n_qTcut, vector<int> (n_seed, 0));
    //    remove_run_qTcut_TSV.resize(osi->n_qTcut, vector<int> (n_seed, 0));

    /*
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
    */

  }

  vector<string> vs0;
  int process_type = 2;
  int n_particle = 2;
  type_parton.resize(1);
  subprocess_readin(osi->process_class, name, vs0, process_type, n_particle, type_parton);
  // needed ??? can most likely be removed !!!

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



///////////////////////
//  access elements  //
///////////////////////



///////////////
//  methods  //
///////////////

// Determines for each subprocess (for the selected scales, qTcut values, etc.):
// error2_time
// time_per_event
// used_runtime
// used_n_event
void summary_subprocess::readin_runtime(string & result_moment){
  Logger logger("summary_subprocess::readin_runtime");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  double time_per_event_time = 0.;
  double used_n_event_time = 0.;
  double used_runtime_time = 0.;

  //if (osi->switch_TSV && osi->name_reference_TSV != ""){
  logger << LOG_DEBUG << "osi->switch_TSV = " << osi->switch_TSV << endl;
  logger << LOG_DEBUG << "osi->switch_CV = " << osi->switch_CV << endl;
  
  if (osi->switch_TSV){remove_run_qTcut = remove_run_qTcut_TSV;}
  else if (osi->switch_CV){remove_run_qTcut = remove_run_qTcut_CV;}
  else {logger << LOG_ERROR << "Should not happen! Run-time determination is not correctly selected." << endl; exit(1);}

  logger << LOG_DEBUG << "ygeneric->directory_runtime_estimate = " << ygeneric->directory_runtime_estimate << endl;
  if (ygeneric->directory_runtime_estimate != ""){
    string result_file = result_moment + "_" + name + ".dat";

    logger << LOG_DEBUG << "ycontribution->directory[0] = " << ycontribution->directory[0] << endl;
    string addtopath = "";
    for (int i_s = ycontribution->directory[0].size() - 1; i_s >= 0; i_s--){
      if (ycontribution->directory[0][i_s] == '/'){
	addtopath = ycontribution->directory[0].substr(0, i_s + 1);
	logger << LOG_DEBUG << "addtopath = " << addtopath << endl;
	break;
      }
    }

    string filename = "../" + addtopath + ygeneric->directory_runtime_estimate + "/result/" + result_file;
    logger << LOG_INFO << "filename = " << filename << endl;

    vector<string> readin;
    char LineBuffer[128];
    ifstream in_result(filename.c_str());  
    while (in_result.getline(LineBuffer, 128)){readin.push_back(LineBuffer);}
    in_result.close();
    //    double error2_time_time = 0.;
    //    double time_per_event_time = 0.;
    //    double N_time = 0.;

    if (readin.size() > 5){
      //      error2_time_time = atof(readin[readin.size() - 2].c_str());
      //      time_per_event_time = atof(readin[readin.size() - 3].c_str());
      used_runtime_time = atof(readin[readin.size() - 4].c_str());
      used_n_event_time = atoll(readin[0].c_str());
      //      N_time = atoll(readin[0].c_str());
    }
    logger << LOG_INFO << "used_runtime_time = " << used_runtime_time << endl;


     /*
    // modify runtime estimates from time directory:
    logger << LOG_INFO << name << ": used_runtime from runs:" << used_runtime << endl;
    logger << LOG_INFO << name << ": used_n_event from runs:" << used_n_event << endl;

    logger << LOG_INFO << name << ": used_runtime from time:" << used_runtime_time << endl;
    logger << LOG_INFO << name << ": used_n_event from time:" << used_n_event_time << endl;
    */
  }


  for (int i_z = 0; i_z < n_seed; i_z++){
    string result_file = result_moment + "_" + name + ".dat";
    logger << LOG_DEBUG_VERBOSE << "result_file = " << result_file << endl;
    string filename = "../" + ycontribution->directory[i_z] + "/result/" + result_file;
    logger << LOG_DEBUG_VERBOSE << "Input: xfilename = " << filename << endl;
    vector<string> readin;
    char LineBuffer[128];
    ifstream in_result(filename.c_str());  
    while (in_result.getline(LineBuffer, 128)){readin.push_back(LineBuffer);}
    in_result.close();
    logger << LOG_DEBUG << "default_size_CV  = " << default_size_CV << endl;
    logger << LOG_DEBUG << "readin.size() = " << readin.size() << endl;
    logger << LOG_DEBUG << "osi->n_scales_CV = " << osi->n_scales_CV << endl;
    logger << LOG_DEBUG << "osi->n_qTcut = " << osi->n_qTcut << endl;
    //    logger << LOG_DEBUG << "osi->active_qTcut = " << osi->active_qTcut << endl;
    logger << LOG_DEBUG << "ygeneric->no_qTcut_runtime_estimate = " << ygeneric->no_qTcut_runtime_estimate << endl;
    logger << LOG_DEBUG << "ycontribution->selection_n_qTcut = " << ycontribution->selection_n_qTcut << endl;
    logger << LOG_DEBUG << "ycontribution->output_n_qTcut = " << ycontribution->output_n_qTcut << endl;
    logger << LOG_DEBUG << "ycontribution->active_qTcut = " << ycontribution->active_qTcut << endl;


    if (readin.size() > 5){
      error2_time_run[i_z] = atof(readin[readin.size() - 2].c_str());
      time_per_event_run[i_z] = atof(readin[readin.size() - 3].c_str());
      used_runtime_run[i_z] = atof(readin[readin.size() - 4].c_str());

      /*
      error2_time_run[i_z] = atof(readin[default_size_CV + 2].c_str());
      time_per_event_run[i_z] = atof(readin[default_size_CV + 1].c_str());
      used_runtime_run[i_z] = atof(readin[default_size_CV].c_str());
      */      
      //      readin.erase(readin.end() - 3, readin.end());
      
      used_n_event_run[i_z] = atoll(readin[0].c_str());
      N[i_z] = atoll(readin[0].c_str());
      

      // runtime correction from directory_runtime_estimate:
      if (used_runtime_time != 0.){
	double new_used_runtime = used_n_event_run[i_z] / used_n_event_time * used_runtime_time;
	time_per_event_run[i_z] = time_per_event_run[i_z] * (new_used_runtime / used_runtime_run[i_z]);
	error2_time_run[i_z] = error2_time_run[i_z] * (new_used_runtime / used_runtime_run[i_z]);
	logger << LOG_INFO << name << "(" << setw(4) << i_z << ") : used_runtime from runs:   " << used_runtime_run[i_z] << endl;
	logger << LOG_INFO << name << "(" << setw(4) << i_z << ") : used_runtime from time:   " << new_used_runtime << endl;
	used_runtime_run[i_z] = new_used_runtime;
      }


      //      double default_result = atof(readin[1].c_str());
      double default_deviation = atof(readin[2].c_str());
      /*
	result_run[i_z] = atof(readin[1].c_str());
	deviation_run[i_z] = atof(readin[2].c_str());
      */
      //      logger << LOG_INFO << "TIME RESULT FROM std  [" << setw(3) << i_z << "] = " << setw(23) << setprecision(15) << result_run[i_z] << " +- " << setw(23) << setprecision(15) << deviation_run[i_z] << endl;

      
      // TSV should be used if 
      //      name_reference_TSV              =       FS
      //      no_scale_ren_reference_TSV      =       1
      //      no_scale_fact_reference_TSV     =       1
      //      no_qTcut_reference_TSV          =       0
      // is set !!!

      logger << LOG_DEBUG << "osi->switch_TSV = " << osi->switch_TSV << endl;
      logger << LOG_DEBUG << "osi->switch_CV = " << osi->switch_CV << endl;
      logger << LOG_DEBUG << "osi->name_reference_TSV = " << osi->name_reference_TSV << endl;
      logger << LOG_DEBUG << "ycontribution->active_qTcut = " << ycontribution->active_qTcut << endl;
      logger << LOG_DEBUG << "ygeneric->no_qTcut_runtime_estimate = " << ygeneric->no_qTcut_runtime_estimate << endl;


     if (osi->switch_TSV && osi->name_reference_TSV != ""){
       //	double default_result = atof(readin[1].c_str());
       ///	double default_deviation = atof(readin[2].c_str());
	//	logger << LOG_DEBUG << "default_result = " << default_result << endl;
	logger << LOG_DEBUG << "default_deviation = " << default_deviation << endl;
	/*///
	logger << LOG_DEBUG << "osi->name_reference_TSV = " << osi->name_reference_TSV << endl;

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

	result_run[i_z] = result_run_TSV[i_z][x_m][x_q][x_s][x_r][x_f];
	deviation_run[i_z] = deviation_run_TSV[i_z][x_m][x_q][x_s][x_r][x_f];
*/
	  
	logger << LOG_INFO << "TIME RESULT FROM TSV  [" << setw(3) << i_z << "] = " << setw(23) << setprecision(15) << result_run[i_z] << " +- " << setw(23) << setprecision(15) << deviation_run[i_z] << endl;

	//	logger << LOG_DEBUG << "default_result = " << default_result << endl;
	logger << LOG_DEBUG << "default_deviation = " << default_deviation << endl;
	logger << LOG_DEBUG << "deviation_run[" << i_z << "] = " << deviation_run[i_z] << endl;

	if (error2_time_run[i_z] != 0. && default_deviation != 0.){error2_time_run[i_z] = error2_time_run[i_z] * pow(deviation_run[i_z] / default_deviation, 2);}
     }

      // later: else if ...
      else if (osi->switch_CV && ycontribution->active_qTcut && ygeneric->no_qTcut_runtime_estimate){

	///	if (ycontribution->active_qTcut && ygeneric->no_qTcut_runtime_estimate){
	int temp_position = 3;
	if (osi->switch_CV){
	  temp_position += osi->n_scales_CV * 2 * ygeneric->no_qTcut_runtime_estimate + ((osi->n_scales_CV - 1) / 2) * 2;
	}
	logger << LOG_DEBUG << "temp_position = " << temp_position << endl;
	logger << LOG_DEBUG << "readin.size() = " << readin.size() << endl;
	//	double default_result = atof(readin[1].c_str());
	///	double default_deviation = atof(readin[2].c_str());

	result_run[i_z] = atof(readin[temp_position].c_str());
	deviation_run[i_z] = atof(readin[temp_position + 1].c_str());

	logger << LOG_INFO << "TIME RESULT FROM CV   [" << setw(3) << i_z << "] = " << setw(23) << setprecision(15) << result_run[i_z] << " +- " << setw(23) << setprecision(15) << deviation_run[i_z] << endl;

	//	logger << LOG_DEBUG << "default_result = " << default_result << endl;
	logger << LOG_DEBUG << "default_deviation = " << default_deviation << endl;
	logger << LOG_DEBUG << "deviation_run[" << i_z << "] = " << deviation_run[i_z] << endl;


	if (error2_time_run[i_z] != 0. && default_deviation != 0.){error2_time_run[i_z] = error2_time_run[i_z] * pow(deviation_run[i_z] / default_deviation, 2);}

      }
      else if (osi->switch_CV){
	//	double default_result = atof(readin[1].c_str());
	///	double default_deviation = atof(readin[2].c_str());
	int temp_position = 3 + ((osi->n_scales_CV - 1) / 2) * 2;
	result_run[i_z] = atof(readin[temp_position].c_str());
	deviation_run[i_z] = atof(readin[temp_position + 1].c_str());
	if (error2_time_run[i_z] != 0. && default_deviation != 0.){error2_time_run[i_z] = error2_time_run[i_z] * pow(deviation_run[i_z] / default_deviation, 2);}
      }
      else {
	logger << LOG_ERROR << "Should not happen! Run-time determination is not correctly selected." << endl;}

    }
    else {
      if (readin.size() != 0){logger << LOG_ERROR << result_file << "   error" << endl;}
      else {logger << LOG_DEBUG << result_file << "   empty" << endl;}
      used_n_event_run[i_z] = 0;
      N[i_z] = 0;
      result_run[i_z] = 0.;
      deviation_run[i_z] = 0.;
      error2_time_run[i_z] = 0.;
      time_per_event_run[i_z] = 0.;
      used_runtime_run[i_z] = 0.;
    }
  }
  logger << LOG_DEBUG << "Averaging of convergence estimators" << endl;
  int counter_e2_t = 0;
  int counter_t_n = 0;
  error2_time = 0.;
  time_per_event = 0.;
  used_runtime = 0.;
  used_n_event = 0;
  
  int this_no_qTcut_runtime_estimate = 0;
  if (ycontribution->active_qTcut){this_no_qTcut_runtime_estimate = ygeneric->no_qTcut_runtime_estimate;}

  for (int i_z = 0; i_z < n_seed; i_z++){
    //    logger << LOG_INFO << remove_run_qTcut[this_no_qTcut_runtime_estimate][i_z] << " remove_run_qTcut[" << this_no_qTcut_runtime_estimate << "][" << i_z << "] = " << remove_run_qTcut[this_no_qTcut_runtime_estimate][i_z] << endl;
    if (!remove_run_qTcut[this_no_qTcut_runtime_estimate][i_z]){
      if (error2_time_run[i_z] != 0.){
	error2_time += error2_time_run[i_z];
	counter_e2_t++;
      }
      if (time_per_event_run[i_z] != 0.){
	time_per_event += time_per_event_run[i_z];
	counter_t_n++;
      }
    }
    else {
      logger << LOG_INFO << ycontribution->type_contribution << "   " << setw(15) << name << "   run " << i_z << "   removed in time calculation:   remove_run_qTcut[" << this_no_qTcut_runtime_estimate << "][" << i_z << "] = " << remove_run_qTcut[this_no_qTcut_runtime_estimate][i_z] << endl;
    }
    // still count the runtimes and event if run is removed.
    used_runtime += used_runtime_run[i_z];
    used_n_event += used_n_event_run[i_z];
    //      cout << "used_runtime_run[" << i_p << "][" << i_z << "] = " << used_runtime_run[i_z] << endl;
  }
  if (counter_e2_t > 0){error2_time = sqrt(error2_time / counter_e2_t);}
  if (counter_t_n > 0){time_per_event = time_per_event / counter_t_n;}


  logger << LOG_INFO << "SUMRUNTIME   " << ycontribution->type_contribution << "   " << setw(15) << name << "   error2_time = " << error2_time << endl;

  logger << LOG_DEBUG << "Averaging of convergence estimators done!" << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


// What are these two functions needed for ???
void summary_subprocess::combine_result_original(){
  Logger logger("summary_subprocess::combine_result_original");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  
  N_all = 0;
  no_results = 0;
  result = 0.;
  deviation = 0.;
  for (int i_z = 0; i_z < n_seed; i_z++){
    if (!(result_run[i_z] == 0. && deviation_run[i_z] == 0.)){
      N_all += N[i_z];
      result += N[i_z] * result_run[i_z];
      deviation += pow((N[i_z] - 1) * deviation_run[i_z], 2) + N[i_z] * pow(result_run[i_z], 2);
      no_results++;
    }
  }
  deviation = sqrt((N_all * deviation - pow(result, 2)) / N_all) / (N_all - 1); 
  result = result / N_all;
  for (int i_z = 0; i_z < n_seed; i_z++){
    if (!(result_run[i_z] == 0. && deviation_run[i_z] == 0.)){
      chi2 += pow((result_run[i_z] - result) / deviation_run[i_z], 2.);
    }
  }
  chi2 = chi2 / no_results;
  for (int i_z = 0; i_z < n_seed; i_z++){
    if (N[i_z] != 0 && deviation_run[i_z] != 0. && result_run[i_z] != 0.){
      logger << LOG_DEBUG << "run " << right << setw(2) << i_z << ":   " << showpoint << setw(13) << right << setprecision(6) << result_run[i_z] << " +- " << setw(13) << setprecision(6) << deviation_run[i_z] << "   " << setw(11) << right << setprecision(4) << (result_run[i_z] - result) / deviation_run[i_z] << " sigma   " << setw(15) << N[i_z] << endl;
    }
  }
  logger << LOG_DEBUG << endl << setw(10) << left << "old av:   " << setw(13) << right << setprecision(6) << result << " +- " << setw(13) << setprecision(6) << deviation << endl << endl;
  //    logger << LOG_DEBUG << "old version: " << i_p << "   " << result << " +- " << deviation << endl;
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void summary_subprocess::combine_result_alternative(){
  Logger logger("summary_subprocess::combine_result_alternative");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ///////////////////////////////////////////////////////////////
  //  combination of different-seed runs for all subprocesses  //
  ///////////////////////////////////////////////////////////////
  double temp_singleweight;

  N_all = 0;
  result = 0.;
  deviation = 0.;
  no_results = 0;
  for (int i_z = 0; i_z < n_seed; i_z++){
    if (N[i_z] != 0 && deviation_run[i_z] != 0. && result_run[i_z] != 0.){
      temp_singleweight = 1. / pow(deviation_run[i_z], 2);
      result += temp_singleweight * result_run[i_z];
      deviation += temp_singleweight;
      N_all += N[i_z];
      no_results++;
    }
  }
  if (N_all != 0){
    result = result / deviation;
    deviation = 1. / sqrt(deviation);
    if (result != result){
      logger << LOG_DEBUG_VERBOSE << N_all << endl;
      for (int i_z = 0; i_z < n_seed; i_z++){logger << LOG_DEBUG_VERBOSE << N[i_z] << endl;}
    }
  }
  else{
    result = 0.;
    deviation = 0.;
  }
  for (int i_z = 0; i_z < n_seed; i_z++){
    if (N[i_z] != 0 && deviation_run[i_z] != 0. && result_run[i_z] != 0.){
      logger << LOG_DEBUG << "run " << right << setw(2) << i_z << ":   " << showpoint << setw(13) << right << setprecision(6) << result_run[i_z] << " +- " << setw(13) << setprecision(6) << deviation_run[i_z] << "   " << setw(11) << right << setprecision(4) << (result_run[i_z] - result) / deviation_run[i_z] << " sigma   " << setw(15) << N[i_z] << " EVENTS" << endl;
    }
  }
  logger << LOG_DEBUG << endl << setw(10) << left << "average:  " << setw(13) << right << setprecision(6) << result << " +- " << setw(13) << setprecision(6) << deviation << endl << endl;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}







//ycontribution->n_bin_distribution_modasym -> osi->extended_distribution[i_d].n_bins

