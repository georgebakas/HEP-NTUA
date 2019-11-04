#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

int check_proceeding_in(vector<string> & readin, int temp_check_size, int counter, int & int_end){
  Logger logger("check_proceeding_in");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
  if (temp_check_size > readin.size()){
    int_end = 2;
    logger << LOG_INFO << "int_end = 2   temp_check_size = " << temp_check_size << " > " << readin.size() << " = readin.size()  @ counter = " << counter << endl;
  }
  else {
    logger << LOG_INFO << "temp_check_size = " << temp_check_size << " <= " << readin.size() << " = readin.size()" << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished --- int_end = " << int_end << endl;
  return int_end;
  /*
  if (temp_check_size > readin.size()){int_end = 2; logger << LOG_INFO << "int_end = 2   temp_check_size = " << temp_check_size << " > " << readin.size() << " = readin.size()  @ counter = " << counter << endl; return;}
  else{logger << LOG_INFO << "temp_check_size = " << temp_check_size << " <= " << readin.size() << " = readin.size()" << endl;}
  */
}

  
void observable_set::perform_proceeding_in(phasespace_set & psi){
  Logger logger("observable_set::perform_proceeding_in");
  logger << LOG_DEBUG << "started" << endl;

  /*
  logger << LOG_DEBUG << "begin" << endl;
  for (int j = 0; j < psi.tau_MC_tau_gamma[0].size(); j++){
    logger << LOG_DEBUG << "tau_MC_tau_gamma[0][" << j << "]" << " = " << setprecision(20) << setw(28) << psi.tau_MC_tau_gamma[0][j] << "   " << double2hexastr(psi.tau_MC_tau_gamma[0][j]) << endl;
  }
  */

  ifstream in_proceeding(filename_proceeding.c_str());  
  char LineBuffer[128];
  vector<string> readin;
  while (in_proceeding.getline(LineBuffer, 128)) {
    readin.push_back(LineBuffer);
  }
  in_proceeding.close();
  for (int i = readin.size() - 1; i >=0; i--){
    if (readin[i][0] == '#'){
      cout << readin[i] << endl;
      readin.erase(readin.begin() + i);
    }
  }
  logger << LOG_INFO << "proceeding_in size: " << readin.size() << endl;



  logger << LOG_INFO << "Start checking expected content of proceeding_in" << endl;

  int counter = 0;
  int temp_check_size = 3 + 3 + 5 + 2;
  //  int temp_check_size = 11;
  if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
  if (switch_CV){
    temp_check_size += 2 * n_qTcut * n_scales_CV;
    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
  }
  if (n_moments != 0){
    for (int nm = 0; nm < n_moments; nm++){
      temp_check_size += 2;
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      if (switch_CV != 0){
	temp_check_size += 2 * n_qTcut * n_scales_CV;
	if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      }
    }
  }
  if (switch_distribution){
    int dist_numbers = 3;
    if (switch_CV){dist_numbers += 2 * n_scales_CV;}
    for (int i = 0; i < dat.size(); i++){
      temp_check_size += dist_numbers * dat[i].n_bins;
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
    }
    for (int i = 0; i < (dddat).size(); i++){
      temp_check_size += dist_numbers * (dddat)[i].n_bins;
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
    }
  }
  if (switch_TSV){
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      temp_check_size += 2 * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s];
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      if (active_qTcut == 1){
	temp_check_size += 2 * n_qTcut * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s];
	if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      }
    }
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      temp_check_size += 2 * n_moments * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s];
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
    }
    int next_no_qTcut = 0;
    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
      if (next_no_qTcut == no_qTcut_distribution.size()){continue;}
      if (i_q != no_qTcut_distribution[next_no_qTcut]){continue;}
      for (int i_d = 0; i_d < dat.size(); i_d++){
	for (int i_b = 0; i_b < dat[i_d].n_bins; i_b++){
	  temp_check_size += 1;
	  if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	    temp_check_size += 2 * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s];
	    if (temp_check_size > readin.size()){int_end = 2; logger << LOG_INFO << "int_end = 2   temp_check_size = " << temp_check_size << " > " << readin.size() << " = readin.size()  @ counter = " << counter << endl; return;}
	    else{logger << LOG_INFO << "temp_check_size = " << temp_check_size << " <= " << readin.size() << " = readin.size()" << endl;}
	  }
	}
      }
      for (int i_d = 0; i_d < dddat.size(); i_d++){
	// not needed here !!!	int i_ddd = dat.size() + i_d;
	for (int i_b = 0; i_b < dddat[i_d].n_bins; i_b++){
	  temp_check_size += 1;
	  if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	    temp_check_size += 2 * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s];
	    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	  }
	}
      }
      next_no_qTcut++;
    }
  }
  temp_check_size += 1;
  if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}



  
  int temp_psi_tau_opt_end = atoi(readin[temp_check_size - 1].c_str());
  cout << "# temp_psi_tau_opt_end = " << temp_psi_tau_opt_end << endl;

  if (temp_psi_tau_opt_end != -1){
    temp_check_size += psi.IS_tau.n_gridsize;
    //    temp_check_size += psi_n_tau_bins;
    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
  }
  if (temp_psi_tau_opt_end == 0){
    temp_check_size += 4 * psi.IS_tau.n_gridsize;
    //    temp_check_size += 4 * psi_n_tau_bins;
    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
  }

  temp_check_size += 1;
  if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}

  cout << "# temp_psi_x1x2_opt_end = " << readin[temp_check_size - 1] << endl;
  int temp_psi_x1x2_opt_end = atoi(readin[temp_check_size - 1].c_str());
  cout << "# temp_psi_x1x2_opt_end = " << temp_psi_x1x2_opt_end << endl;

  if (temp_psi_x1x2_opt_end != -1){
    temp_check_size += psi.IS_x1x2.n_gridsize;
    //    temp_check_size += psi_n_x1x2_bins;
    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
  }
  if (temp_psi_x1x2_opt_end == 0){
    temp_check_size += 4 * psi.IS_x1x2.n_gridsize;
    //    temp_check_size += 4 * psi_n_x1x2_bins;
    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
  }
  temp_check_size += 1;
  if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}

  int temp_psi_z1z2_opt_end = atoi(readin[temp_check_size - 1].c_str());
  cout << "# temp_psi_z1z2_opt_end = " << temp_psi_z1z2_opt_end << endl;

  if (csi->class_contribution_CS_collinear){
  if (temp_psi_z1z2_opt_end != -1){
    for (int i_z = 1; i_z < 3; i_z++){
      temp_check_size += psi_n_z1z2_bins;
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      if (temp_psi_z1z2_opt_end == 0){
      temp_check_size += 4 * psi_n_z1z2_bins;
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      }
    }
  }
  }
  
  logger << LOG_INFO << "before psi_random_manager.proceeding_in" << endl;
  // check if this works fine !!!
  //  int temp_temp_check_size = temp_check_size;
  psi_random_manager.check_proceeding_in(int_end, temp_check_size, readin);
  //  temp_check_size += (temp_check_size - temp_temp_check_size);
  logger << LOG_INFO << "after psi_random_manager.proceeding_in" << endl;

  temp_check_size += 1;
  if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}

  int temp_psi_MC_opt_end = atoi(readin[temp_check_size - 1].c_str());
  //  int temp_psi_MC_opt_end = atoi(readin[temp_check_size - 2].c_str());
  cout << "# temp_psi_MC_opt_end = " << temp_psi_MC_opt_end << endl;
  ///  int temp_MC_n_channel = atoi(readin[temp_check_size - 1].c_str());
  ///  cout << "# temp_MC_n_channel = " << temp_MC_n_channel << endl;


  /*
  // -> psi.MC_phasespace.end_optimization
  if (psi.MC_phasespace.active_optimization != -1){
    temp_check_size += 1;
    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
    int temp_psi_MC_phasespace_end_optimization = atoi(readin[temp_check_size - 1].c_str());
    /// ???  int temp_psi_MC_phasespace_end_optimization = atoi(readin[temp_check_size - 2].c_str());
    cout << "# temp_psi_MC_phasespace_end_optimization = " << temp_psi_MC_phasespace_end_optimization << endl;

    temp_check_size += psi.MC_phasespace.n_channel;
    //    temp_check_size += psi_MC_n_channel;
    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}

    //temp_psi_MC_opt_end ???
    //  if (temp_psi_MC_opt_end == 0){
    if (psi.MC_phasespace.end_optimization == 0){
      temp_check_size += 1;
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      int temp_n_opt_steps = atoi(readin[temp_check_size - 1].c_str());
      cout << "# temp_n_opt_steps = " << temp_n_opt_steps << endl;

      temp_check_size += 3 * psi.MC_phasespace.n_channel;
      //      temp_check_size += 3 * psi_MC_n_channel;
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}

      //      if (temp_check_size != readin.size()){
      {
	temp_check_size += 2;
	if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	
	temp_check_size += temp_n_opt_steps * (psi.MC_phasespace.n_channel + 1); // ??? + 1
	//	temp_check_size += temp_n_opt_steps * (psi_MC_n_channel + 1); // ??? + 1
	//      temp_check_size += temp_n_opt_steps * psi_MC_n_channel + 1; // ??? + 1
	if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      }
    }
  }
*/


  // -> psi.MC_phasespace.end_optimization
  if (psi.MC_phasespace.active_optimization != -1){
    temp_check_size += 1;
    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
    int temp_psi_MC_phasespace_end_optimization = atoi(readin[temp_check_size - 1].c_str());
    logger << LOG_INFO << "# temp_psi_MC_phasespace_end_optimization = " << temp_psi_MC_phasespace_end_optimization << endl;

    temp_check_size += psi.MC_phasespace.n_channel;
    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}

    if (temp_psi_MC_phasespace_end_optimization == 0){
      //    if (psi.MC_phasespace.end_optimization == 0){
      temp_check_size += 1;
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      int temp_n_opt_steps = atoi(readin[temp_check_size - 1].c_str());
      logger << LOG_INFO << "# temp_n_opt_steps = " << temp_n_opt_steps << endl;
      temp_check_size += 3 * psi.MC_phasespace.n_channel;
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      //      if (temp_check_size != readin.size()){
      {
	temp_check_size += 2;
	if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	temp_check_size += temp_n_opt_steps * (psi.MC_phasespace.n_channel + 1); // ??? + 1
	logger << LOG_INFO << "temp_n_opt_steps = " << temp_n_opt_steps << endl;
	logger << LOG_INFO << "psi.MC_phasespace.n_channel = " << psi.MC_phasespace.n_channel << endl;
	if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      }
    }
  }


  // -> psi.MC_tau.end_optimization
  if (psi.MC_tau.active_optimization != -1){
    temp_check_size += 1;
    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
    int temp_psi_MC_tau_end_optimization = atoi(readin[temp_check_size - 1].c_str());
    logger << LOG_INFO << "# temp_psi_MC_tau_end_optimization = " << temp_psi_MC_tau_end_optimization << endl;

    temp_check_size += psi.MC_tau.n_channel;
    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
    if (temp_psi_MC_tau_end_optimization == 0){
      //    if (psi.MC_tau.end_optimization == 0){
      temp_check_size += 1;
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      int temp_n_opt_steps = atoi(readin[temp_check_size - 1].c_str());
      logger << LOG_INFO << "# temp_n_opt_steps = " << temp_n_opt_steps << endl;
      temp_check_size += 3 * psi.MC_tau.n_channel;
      if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      //      if (temp_check_size != readin.size()){
      {
	temp_check_size += 2;
	if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	temp_check_size += temp_n_opt_steps * (psi.MC_tau.n_channel + 1); // ??? + 1
	logger << LOG_INFO << "temp_n_opt_steps = " << temp_n_opt_steps << endl;
	logger << LOG_INFO << "psi.MC_tau.n_channel = " << psi.MC_tau.n_channel << endl;
	if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
      }
    }
  }


  if (csi->class_contribution_CS_real){
    for (int i_a = 1; i_a < psi.MC_x_dipole.size(); i_a++){
      // -> psi.MC_x_dipole[i_a].end_optimization
      if (psi.MC_x_dipole[i_a].active_optimization != -1){
	temp_check_size += 1;
	if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	int temp_psi_MC_x_dipole_i_a_end_optimization = atoi(readin[temp_check_size - 1].c_str());
	logger << LOG_INFO << "# temp_psi_MC_x_dipole_i_a_end_optimization = " << temp_psi_MC_x_dipole_i_a_end_optimization << endl;
	
	temp_check_size += psi.MC_x_dipole[i_a].n_channel;
	if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	if (temp_psi_MC_x_dipole_i_a_end_optimization == 0){
	  //	if (psi.MC_x_dipole[i_a].end_optimization == 0){
	  temp_check_size += 1;
	  if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	  int temp_n_opt_steps = atoi(readin[temp_check_size - 1].c_str());
	  logger << LOG_INFO << "# temp_n_opt_steps = " << temp_n_opt_steps << endl;
	  temp_check_size += 3 * psi.MC_x_dipole[i_a].n_channel;
	  if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	  //      if (temp_check_size != readin.size()){
	  {
	    temp_check_size += 2;
	    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	    temp_check_size += temp_n_opt_steps * (psi.MC_x_dipole[i_a].n_channel + 1); // ??? + 1
	    if (check_proceeding_in(readin, temp_check_size, counter, int_end) == 2){return;}
	  }
	}
      }
      
    }
  }

  
  // Improve this last paragraph !!!
  // Add MC_tau etc. !!!

  logger << LOG_INFO << "Checking expected content of proceeding_in succesfully finished." << endl;

  logger << LOG_INFO << "temp_check_size = " << temp_check_size << " ? " << readin.size() << " = readin.size()" << endl;

  











  int proc = 0;
  for (int i = 0; i < 3; i++){
    psi_sran[i] = hexastr2double(readin[i]);
  }
  proc += 3;

  h = atoi(readin[proc].c_str());
  min = atoi(readin[proc + 1].c_str());
  sec = atoi(readin[proc + 2].c_str());
  sec_import = sec + 60 * min + 3600 * h;
  // former typo !!!  sec_import = sec + 60 * min + 600 * h;
  proc += 3;

  logger << LOG_INFO << "Imported runtime: " << h << " h " << min << " min " << sec << " sec   ->   sec_import = " << sec_import << endl;
  /* 
  //  old order:
  psi_i_rej = atol(readin[proc].c_str());
  psi_i_acc = atol(readin[proc + 1].c_str());
  psi_i_nan = atoi(readin[proc + 2].c_str());
  proc += 3;
  */
  psi.i_gen = atol(readin[proc].c_str());
  psi.i_acc = atol(readin[proc + 1].c_str());
  psi.i_rej = atol(readin[proc + 2].c_str());
  psi.i_tec = atoi(readin[proc + 3].c_str());
  psi.i_nan = atoi(readin[proc + 4].c_str());
  // add psi_i_gen and psi_i_tec and change order !!! gen - acc - rej - tec - nan ???
  proc += 5;
  psi.last_step_mode = *psi.i_step_mode;

  full_sum_weight = hexastr2double(readin[proc]);
  full_sum_weight2 = hexastr2double(readin[proc + 1]);
  proc += 2;

  if (switch_CV){
    // !!! Reduce output by putting output_n_qTcut !!!
    for (int c = 0; c < n_qTcut; c++){
      for (int s = 0; s < n_scales_CV; s++){
	full_sum_weight_CV[c][s] = hexastr2double(readin[proc + 2 * (c * n_scales_CV + s)]);
	full_sum_weight2_CV[c][s] = hexastr2double(readin[proc + 1 + 2 * (c * n_scales_CV + s)]);
      }
    }
    proc += 2 * n_qTcut * n_scales_CV;
  }

  if (n_moments != 0){
    for (int nm = 0; nm < n_moments; nm++){
      full_sum_moment[nm] = hexastr2double(readin[proc]);
      full_sum_moment2[nm] = hexastr2double(readin[proc + 1]);
      proc += 2;
      if (switch_CV != 0){
	for (int c = 0; c < n_qTcut; c++){
	  for (int s = 0; s < n_scales_CV; s++){
	    full_sum_moment_CV[nm][c][s] = hexastr2double(readin[proc + 2 * (c * n_scales_CV + s)]);
	    full_sum_moment2_CV[nm][c][s] = hexastr2double(readin[proc + 1 + 2 * (c * n_scales_CV + s)]);
	  }
	}
	proc += 2 * n_qTcut * full_sum_moment_CV[nm][0].size();
      }
    }
  }
  if (switch_distribution){
    int dist_numbers = 3;
    if (switch_CV){dist_numbers += 2 * n_scales_CV;}
    for (int i_d = 0; i_d < dat.size(); i_d++){
      for (int n = 0; n < dat[i_d].n_bins; n++){
	bin_counts[i_d][n] = atoi(readin[proc + dist_numbers * n].c_str());
	bin_weight[i_d][n] = hexastr2double(readin[proc + dist_numbers * n + 1]);
	bin_weight2[i_d][n] = hexastr2double(readin[proc + dist_numbers * n + 2]);
	/*
	bin_counts[i_d][n] = atoi(readin[proc + dist_numbers * n + 2].c_str());
	bin_weight[i_d][n] = hexastr2double(readin[proc + dist_numbers * n]);
	bin_weight2[i_d][n] = hexastr2double(readin[proc + dist_numbers * n + 1]);
	*/
	if (switch_CV){
	  for (int s = 0; s < n_scales_CV; s++){
	    bin_weight_CV[s][i_d][n] = hexastr2double(readin[proc + dist_numbers * n + 2 * s + 3]);
	    bin_weight2_CV[s][i_d][n] = hexastr2double(readin[proc + dist_numbers * n + 2 * s + 4]);
	  }
	}
      }
      proc += dist_numbers * dat[i_d].n_bins;
    }

    for (int i_ddd = 0; i_ddd < (dddat).size(); i_ddd++){
      for (int n = 0; n < (dddat)[i_ddd].n_bins; n++){
	bin_counts[dat.size() + i_ddd][n] = atoi(readin[proc + dist_numbers * n].c_str());
	bin_weight[dat.size() + i_ddd][n] = hexastr2double(readin[proc + dist_numbers * n + 1]);
	bin_weight2[dat.size() + i_ddd][n] = hexastr2double(readin[proc + dist_numbers * n + 2]);
	/*
	bin_counts[dat.size() + i_ddd][n] = atoi(readin[proc + dist_numbers * n + 2].c_str());
	bin_weight[dat.size() + i_ddd][n] = hexastr2double(readin[proc + dist_numbers * n]);
	bin_weight2[dat.size() + i_ddd][n] = hexastr2double(readin[proc + dist_numbers * n + 1]);
	*/
	if (switch_CV){
	  for (int s = 0; s < n_scales_CV; s++){
	    bin_weight_CV[s][dat.size() + i_ddd][n] = hexastr2double(readin[proc + dist_numbers * n + 2 * s + 3]);
	    bin_weight2_CV[s][dat.size() + i_ddd][n] = hexastr2double(readin[proc + dist_numbers * n + 2 * s + 4]);
	  }
	}
      }
      proc += dist_numbers * (dddat)[i_ddd].n_bins;
    }
  }


  if (switch_TSV){
    logger << LOG_INFO << "# TSV input started" << endl;
    logger << LOG_INFO << "# TSV size start: proc = " << proc << endl;
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	  fullsum_weight_TSV[i_s][i_r][i_f] = hexastr2double(readin[proc + 2 * (i_r * n_scale_fact_TSV[i_s] + i_f)]);
	  fullsum_weight2_TSV[i_s][i_r][i_f] = hexastr2double(readin[proc + 2 * (i_r * n_scale_fact_TSV[i_s] + i_f) + 1]);
	}
      }
      proc += 2 * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s];

      if (active_qTcut == 1){
	for (int i_q = 0; i_q < n_qTcut; i_q++){
	  for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	      fullsum_weight_qTcut_TSV[i_s][i_q][i_r][i_f] = hexastr2double(readin[proc + 2 * (i_q * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s] + i_r * n_scale_fact_TSV[i_s] + i_f)]);
	      fullsum_weight2_qTcut_TSV[i_s][i_q][i_r][i_f] = hexastr2double(readin[proc + 2 * (i_q * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s] + i_r * n_scale_fact_TSV[i_s] + i_f) + 1]);
	    }
	  }
	}
	proc += 2 * n_qTcut * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s];
      }
    }
    
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      for (int i_m = 0; i_m < n_moments; i_m++){
	for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	    fullsum_moment_TSV[i_s][i_m][i_r][i_f] = hexastr2double(readin[proc + 2 * (i_m * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s] + i_r * n_scale_fact_TSV[i_s] + i_f)]);
	    fullsum_moment2_TSV[i_s][i_m][i_r][i_f] = hexastr2double(readin[proc + 2 * (i_m * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s] + i_r * n_scale_fact_TSV[i_s] + i_f) + 1]);
	  }
	}
      }
      proc += 2 * n_moments * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s];
    }
    
    logger << LOG_INFO << "no_qTcut_distribution.size() = " << no_qTcut_distribution.size() << endl;

    for (int x_q = 0; x_q < no_qTcut_distribution.size(); x_q++){
      logger << LOG_INFO << "no_qTcut_distribution[" << x_q << "] = " << no_qTcut_distribution[x_q] << "   selection_qTcut_distribution[" << x_q << "] = " << no_qTcut_distribution[x_q] << endl;
    }

    int next_no_qTcut = 0;
    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
      //    for (int i_q = 0; i_q < n_qTcut; i_q++){
      logger << LOG_INFO << "no_qTcut_distribution[" << next_no_qTcut << "] = " << no_qTcut_distribution[next_no_qTcut] << endl;
      if (next_no_qTcut == no_qTcut_distribution.size()){continue;}
      if (i_q != no_qTcut_distribution[next_no_qTcut]){continue;}
      logger << LOG_INFO << "readin at i_q = " << i_q << "   proc = " << proc << endl;
      for (int i_d = 0; i_d < dat.size(); i_d++){
	for (int i_b = 0; i_b < dat[i_d].n_bins; i_b++){
	  bin_count_TSV[i_q][i_d][i_b] = atoi(readin[proc].c_str());
	  proc++;
	  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][i_b] = hexastr2double(readin[proc + 2 * (i_r * n_scale_fact_TSV[i_s] + i_f)]);
		bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][i_b] = hexastr2double(readin[proc + 2 * (i_r * n_scale_fact_TSV[i_s] + i_f) + 1]);
	      }
	    }
	    proc = proc + 2 * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s];
	  }
	}
      }
      for (int i_d = 0; i_d < dddat.size(); i_d++){
	int i_ddd = dat.size() + i_d;
	for (int i_b = 0; i_b < dddat[i_d].n_bins; i_b++){
	  bin_count_TSV[i_q][i_ddd][i_b] = atoi(readin[proc].c_str());
	  proc++;
	  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		bin_weight_TSV[i_s][i_q][i_r][i_f][i_ddd][i_b] = hexastr2double(readin[proc + 2 * (i_r * n_scale_fact_TSV[i_s] + i_f)]);
		bin_weight2_TSV[i_s][i_q][i_r][i_f][i_ddd][i_b] = hexastr2double(readin[proc + 2 * (i_r * n_scale_fact_TSV[i_s] + i_f) + 1]);
	      }
	    }
	    proc = proc + 2 * n_scale_ren_TSV[i_s] * n_scale_fact_TSV[i_s];
	  }
	}
      }
      next_no_qTcut++;
    }
    logger << LOG_INFO << "# TSV size end: proc = " << proc << endl;
    logger << LOG_INFO << "# TSV input finished" << endl;
  }


  

  // One could avoid all the following output for the case of a completely finalized optimization phase (like in all standard non-grid runs) !!!
  // Introduce a general switch like e.g. 'opt_end', and the following is only read if that switch is 0 (or the other way round) !!!

  psi.IS_tau.end_optimization = atoi(readin[proc].c_str());
  logger << LOG_INFO << "# psi.IS_tau.end_optimization = " << psi.IS_tau.end_optimization << endl;
  proc += 1;
  if (psi.IS_tau.end_optimization != -1){
    for (int i = 0; i < psi.IS_tau.n_gridsize; i++){
      psi.IS_tau.alpha[i] = hexastr2double(readin[proc + i]);
    }
    proc += psi.IS_tau.n_gridsize;
    psi.IS_tau.beta[0] = psi.IS_tau.alpha[0];
    for (int i = 1; i < psi.IS_tau.beta.size(); i++){psi.IS_tau.beta[i] = psi.IS_tau.beta[i - 1] + psi.IS_tau.alpha[i];}

  }
  if (psi.IS_tau.end_optimization == 0){
    for (int i = 0; i < psi.IS_tau.n_gridsize; i++){
      psi.IS_tau.n_acc_channel[i] = atoi(readin[proc + 4 * i].c_str());
      psi.IS_tau.n_rej_channel[i] = atoi(readin[proc + 4 * i + 1].c_str());
      psi.IS_tau.sum_channel_weight[i] = hexastr2double(readin[proc + 4 * i + 2]);
      psi.IS_tau.sum_channel_weight2[i] = hexastr2double(readin[proc + 4 * i + 3]);
    }
    proc += 4 * psi.IS_tau.n_gridsize;
  }
  /*
  psi_tau_opt_end = atoi(readin[proc].c_str());
  logger << LOG_INFO << "# psi_tau_opt_end = " << psi_tau_opt_end << endl;
  proc += 1;
  if (psi_tau_opt_end != -1){
    for (int i = 0; i < psi_n_tau_bins; i++){
      psi_tau_alpha[i] = hexastr2double(readin[proc + i]);
    }
    proc += psi_n_tau_bins;
    psi_tau_beta[0] = psi_tau_alpha[0];
    for (int i = 1; i < psi_tau_beta.size(); i++){psi_tau_beta[i] = psi_tau_beta[i - 1] + psi_tau_alpha[i];}

  }
  if (psi_tau_opt_end == 0){
    for (int i = 0; i < psi_n_tau_bins; i++){
      psi_tau_counts_channel[i] = atoi(readin[proc + 4 * i].c_str());
      psi_tau_cuts_channel[i] = atoi(readin[proc + 4 * i + 1].c_str());
      //      psi_tau_counts_channel[i] = atof(readin[proc + 4 * i].c_str());
      //      psi_tau_cuts_channel[i] = atof(readin[proc + 4 * i + 1].c_str());
      psi_sum_tau_channel_weight[i] = hexastr2double(readin[proc + 4 * i + 2]);
      psi_sum_tau_channel_weight2[i] = hexastr2double(readin[proc + 4 * i + 3]);
      //      psi_sum_tau_channel_weight[i] = atof(readin[proc + 4 * i + 2].c_str());
      //      psi_sum_tau_channel_weight2[i] = atof(readin[proc + 4 * i + 3].c_str());
    }
    proc += 4 * psi_n_tau_bins;
  }
  */
  logger << LOG_INFO << "# TSV tau_opt end: proc = " << proc << endl;



  psi.IS_x1x2.end_optimization = atoi(readin[proc].c_str());
  logger << LOG_INFO << "# psi.IS_x1x2.end_optimization = " << psi.IS_x1x2.end_optimization << endl;
  proc += 1;
  if (psi.IS_x1x2.end_optimization != -1){
    for (int i = 0; i < psi.IS_x1x2.n_gridsize; i++){
      psi.IS_x1x2.alpha[i] = hexastr2double(readin[proc + i]);
    }
    proc += psi.IS_x1x2.n_gridsize;
    psi.IS_x1x2.beta[0] = psi.IS_x1x2.alpha[0];
    for (int i = 1; i < psi.IS_x1x2.beta.size(); i++){psi.IS_x1x2.beta[i] = psi.IS_x1x2.beta[i - 1] + psi.IS_x1x2.alpha[i];}

  }
  if (psi.IS_x1x2.end_optimization == 0){
    for (int i = 0; i < psi.IS_x1x2.n_gridsize; i++){
      psi.IS_x1x2.n_acc_channel[i] = atoi(readin[proc + 4 * i].c_str());
      psi.IS_x1x2.n_rej_channel[i] = atoi(readin[proc + 4 * i + 1].c_str());
      psi.IS_x1x2.sum_channel_weight[i] = hexastr2double(readin[proc + 4 * i + 2]);
      psi.IS_x1x2.sum_channel_weight2[i] = hexastr2double(readin[proc + 4 * i + 3]);
    }
    proc += 4 * psi.IS_x1x2.n_gridsize;
  }


  /*
  psi_x1x2_opt_end = atoi(readin[proc].c_str());
  logger << LOG_INFO << "# psi_x1x2_opt_end = " << psi_x1x2_opt_end << endl;
  proc += 1;

  if (psi_x1x2_opt_end != -1){
    for (int i = 0; i < psi_n_x1x2_bins; i++){
      psi_x1x2_alpha[i] = hexastr2double(readin[proc + i]);
      //      psi_x1x2_alpha[i] = atof(readin[proc + 2 * i].c_str());
      //      psi_x1x2_beta[i] = atof(readin[proc + 2 * i + 1].c_str());
    }
    proc += psi_n_x1x2_bins;
    //    proc += 2 * psi_n_x1x2_bins;
    psi_x1x2_beta[0] = psi_x1x2_alpha[0];
    for (int i = 1; i < psi_x1x2_beta.size(); i++){psi_x1x2_beta[i] = psi_x1x2_beta[i - 1] + psi_x1x2_alpha[i];}
  }

  if (psi_x1x2_opt_end == 0){
    for (int i = 0; i < psi_n_x1x2_bins; i++){
      psi_x1x2_counts_channel[i] = atoi(readin[proc + 4 * i].c_str());
      psi_x1x2_cuts_channel[i] = atoi(readin[proc + 4 * i + 1].c_str());
      //      psi_x1x2_counts_channel[i] = atof(readin[proc + 4 * i].c_str());
      //      psi_x1x2_cuts_channel[i] = atof(readin[proc + 4 * i + 1].c_str());
      psi_sum_x1x2_channel_weight[i] = hexastr2double(readin[proc + 4 * i + 2]);
      psi_sum_x1x2_channel_weight2[i] = hexastr2double(readin[proc + 4 * i + 3]);
      //      psi_sum_x1x2_channel_weight[i] = atof(readin[proc + 4 * i + 2].c_str());
      //      psi_sum_x1x2_channel_weight2[i] = atof(readin[proc + 4 * i + 3].c_str());
    }
    proc += 4 * psi_n_x1x2_bins;
  }
  */
  logger << LOG_INFO << "# TSV x1x2_opt end: proc = " << proc << endl;

  // temporary !!!
  if (csi->class_contribution_CS_collinear){
    //  only for backwards-compatibility reasons -> can be removed !!!
    psi.IS_z1z2[1].end_optimization = atoi(readin[proc].c_str());
    psi.IS_z1z2[2].end_optimization = atoi(readin[proc].c_str());
  logger << LOG_INFO << "# psi.IS_z1z2[1].end_optimization = " << psi.IS_z1z2[1].end_optimization << endl;
  logger << LOG_INFO << "# psi.IS_z1z2[2].end_optimization = " << psi.IS_z1z2[2].end_optimization << endl;
  }
  proc += 1;

  
  if (csi->class_contribution_CS_collinear){
    for (int i_z = 1; i_z < 3; i_z++){
      //  later...
      //  add this line !!!
      //  psi.IS_z1z2[i_z].end_optimization = atoi(readin[proc].c_str());
      //  proc += 1;
      
      if (psi.IS_z1z2[i_z].end_optimization != -1){
	for (int i = 0; i < psi_n_z1z2_bins; i++){
	  psi.IS_z1z2[i_z].alpha[i] = hexastr2double(readin[proc + i]);
	}
	proc += psi.IS_z1z2[i_z].n_gridsize;
	
	psi.IS_z1z2[i_z].beta[0] = psi.IS_z1z2[i_z].alpha[0];
	for (int i = 1; i < psi.IS_z1z2[i_z].beta.size(); i++){psi.IS_z1z2[i_z].beta[i] = psi.IS_z1z2[i_z].beta[i - 1] + psi.IS_z1z2[i_z].alpha[i];}
	
	if (psi.IS_z1z2[i_z].end_optimization == 0){
	  for (int i = 0; i < psi_n_z1z2_bins; i++){
	    psi.IS_z1z2[i_z].n_acc_channel[i] = atoi(readin[proc + 4 * i].c_str());
	    psi.IS_z1z2[i_z].n_rej_channel[i] = atoi(readin[proc + 4 * i + 1].c_str());
	    psi.IS_z1z2[i_z].sum_channel_weight[i] = hexastr2double(readin[proc + 4 * i + 2]);
	    psi.IS_z1z2[i_z].sum_channel_weight2[i] = hexastr2double(readin[proc + 4 * i + 3]);
	  }
	  proc += 4 * psi.IS_z1z2[i_z].n_gridsize;
	}
      }
    }
    logger << LOG_INFO << "# TSV z1z2_opt end: proc = " << proc << endl;
  }
 
  /*
  psi_z1z2_opt_end = atoi(readin[proc].c_str());

  logger << LOG_INFO << "# psi_z1z2_opt_end = " << psi_z1z2_opt_end << endl;
  proc += 1;

  if (psi_z1z2_opt_end != -1){
    for (int i_z = 1; i_z < 3; i_z++){
      for (int i = 0; i < psi_n_z1z2_bins; i++){
	psi_z1z2_alpha[i_z][i] = hexastr2double(readin[proc + i]);
	//	psi_z1z2_alpha[i_z][i] = atof(readin[proc + 2 * i].c_str());
	//	psi_z1z2_beta[i_z][i] = atof(readin[proc + 2 * i + 1].c_str());
      }
      proc += psi_n_z1z2_bins;
      //      proc += 2 * psi_n_z1z2_bins;

      psi_z1z2_beta[i_z][0] = psi_z1z2_alpha[i_z][0];
      for (int i = 1; i < psi_z1z2_beta[i_z].size(); i++){psi_z1z2_beta[i_z][i] = psi_z1z2_beta[i_z][i - 1] + psi_z1z2_alpha[i_z][i];}

      if (psi_z1z2_opt_end == 0){
	for (int i = 0; i < psi_n_z1z2_bins; i++){
	  psi_z1z2_counts_channel[i_z][i] = atoi(readin[proc + 4 * i].c_str());
	  psi_z1z2_cuts_channel[i_z][i] = atoi(readin[proc + 4 * i + 1].c_str());
	  //	  psi_z1z2_counts_channel[i_z][i] = atof(readin[proc + 4 * i].c_str());
	  //	  psi_z1z2_cuts_channel[i_z][i] = atof(readin[proc + 4 * i + 1].c_str());
	  psi_sum_z1z2_channel_weight[i_z][i] = hexastr2double(readin[proc + 4 * i + 2]);
	  psi_sum_z1z2_channel_weight2[i_z][i] = hexastr2double(readin[proc + 4 * i + 3]);
	  //	  psi_sum_z1z2_channel_weight[i_z][i] = atof(readin[proc + 4 * i + 2].c_str());
	  //	  psi_sum_z1z2_channel_weight2[i_z][i] = atof(readin[proc + 4 * i + 3].c_str());
	}
	proc += 4 * psi_n_z1z2_bins;
      }
    }
  }
  logger << LOG_INFO << "# TSV z1z2_opt end: proc = " << proc << endl;

   */
  /*
  logger << LOG_DEBUG << "before random_manager" << endl;
  for (int j = 0; j < psi.tau_MC_tau_gamma[0].size(); j++){
    logger << LOG_DEBUG << "tau_MC_tau_gamma[0][" << j << "]" << " = " << setprecision(20) << setw(28) << psi.tau_MC_tau_gamma[0][j] << "   " << double2hexastr(psi.tau_MC_tau_gamma[0][j]) << endl;
  }
  */
  psi_random_manager.proceeding_in(proc, readin);

  logger << LOG_DEBUG << "after random_manager" << endl;
  for (int j = 0; j < psi.tau_MC_tau_gamma[0].size(); j++){
    logger << LOG_DEBUG << "tau_MC_tau_gamma[0][" << j << "]" << " = " << setprecision(20) << setw(28) << psi.tau_MC_tau_gamma[0][j] << "   " << double2hexastr(psi.tau_MC_tau_gamma[0][j]) << endl;
  }

  logger << LOG_INFO << "# TSV random_manager end: proc = " << proc << endl;




  // Solve issues that MC_phasespace does not work correctly after resumption !!!
  // Extend for MC_tau and MC_x_dipole !!!

  // int_end -> MC_opt_end

  psi_MC_opt_end = atoi(readin[proc].c_str());
  proc += 1;


  if (psi.MC_phasespace.active_optimization != -1){
  
    psi.MC_phasespace.end_optimization = atoi(readin[proc].c_str());
    proc += 1;
    for (int j = 0; j < psi.MC_phasespace.n_channel; j++){
      psi.MC_phasespace.alpha[j] = hexastr2double(readin[proc + j]);
    }
    proc += psi.MC_phasespace.n_channel;

    psi.MC_phasespace.beta[0] = psi.MC_phasespace.alpha[0];
    for (int i = 1; i < psi.MC_phasespace.beta.size(); i++){psi.MC_phasespace.beta[i] = psi.MC_phasespace.beta[i - 1] + psi.MC_phasespace.alpha[i];}
    
    //    if (psi_MC_opt_end == 0){
    if (psi.MC_phasespace.end_optimization == 0){
      int n_opt_steps = atoi(readin[proc].c_str());
      proc += 1;
      for (int j = 0; j < psi.MC_phasespace.n_channel; j++){
	psi.MC_phasespace.n_acc_channel[j] = atoi(readin[proc + 3 * j].c_str());
	psi.MC_phasespace.n_rej_channel[j] = atoi(readin[proc + 1 + 3 * j].c_str());
	psi.MC_phasespace.w_channel_step_sum[j] = hexastr2double(readin[proc + 2 + 3 * j]);
      }
      proc += 3 * psi.MC_phasespace.n_channel;
      
      //      if (proc != readin.size()){
      {
	// x_alpha_it_min needs to be added here !!!
	psi.MC_phasespace.i_alpha_it = atoi(readin[proc].c_str());
	psi.MC_phasespace.x_alpha_it_min = atoi(readin[proc + 1].c_str());
	proc += 2;

	psi.MC_phasespace.diff_w.resize(n_opt_steps);
	for (int i = 0; i < n_opt_steps; i++){
	  psi.MC_phasespace.diff_w[i] = hexastr2double(readin[proc + i]);
	}
	proc += n_opt_steps;
	
	if (n_opt_steps != 0){
	  //	proc += n_opt_steps;
	  //	proc = proc + n_opt_steps;
	  // How is this taken into account ??? !!! 
	  psi.MC_phasespace.alpha_it.resize(n_opt_steps);
	  for (int i = 0; i < n_opt_steps; i++){
	    psi.MC_phasespace.alpha_it[i].resize(psi.MC_phasespace.n_channel);
	    for (int j = 0; j < psi.MC_phasespace.n_channel; j++){
	      psi.MC_phasespace.alpha_it[i][j] = hexastr2double(readin[proc + psi.MC_phasespace.n_channel * i + j]);
	    }
	  }
	  //	  proc += n_opt_steps * (psi.MC_phasespace.n_channel + 1);
	  proc += n_opt_steps * psi.MC_phasespace.n_channel;
	}
      }
      logger << LOG_INFO << "end MC_phasespace" << endl;
    }
  }

  logger << LOG_DEBUG << "before MC_tau" << endl;
  for (int j = 0; j < psi.tau_MC_tau_gamma[0].size(); j++){
    logger << LOG_DEBUG << "tau_MC_tau_gamma[0][" << j << "]" << " = " << setprecision(20) << setw(28) << psi.tau_MC_tau_gamma[0][j] << "   " << double2hexastr(psi.tau_MC_tau_gamma[0][j]) << endl;
  }

  if (psi.MC_tau.active_optimization != -1){
  
    psi.MC_tau.end_optimization = atoi(readin[proc].c_str());
    proc += 1;
    logger << LOG_DEBUG << "psi.MC_tau.end_optimization = " << psi.MC_tau.end_optimization << endl;
    logger << LOG_DEBUG << "psi.MC_tau.n_channel = " << psi.MC_tau.n_channel << endl;
    
    for (int j = 0; j < psi.MC_tau.n_channel; j++){
      psi.MC_tau.alpha[j] = hexastr2double(readin[proc + j]);
    }
    proc += psi.MC_tau.n_channel;

    psi.MC_tau.beta[0] = psi.MC_tau.alpha[0];
    for (int i = 1; i < psi.MC_tau.beta.size(); i++){psi.MC_tau.beta[i] = psi.MC_tau.beta[i - 1] + psi.MC_tau.alpha[i];}

    for (int j = 0; j < psi.MC_tau.n_channel; j++){
      logger << LOG_DEBUG << "psi.MC_tau.alpha[" << j << "] = " << setprecision(20) << setw(28) << psi.MC_tau.alpha[j] << "   " << double2hexastr(psi.MC_tau.alpha[j]) << "   " << "psi.MC_tau.beta[" << j << "] = " << setprecision(20) << setw(28) << psi.MC_tau.beta[j] << "   " << double2hexastr(psi.MC_tau.beta[j]) << endl;
    }
    
    //    if (psi_MC_opt_end == 0){
    if (psi.MC_tau.end_optimization == 0){
      int n_opt_steps = atoi(readin[proc].c_str());
      proc += 1;
      for (int j = 0; j < psi.MC_tau.n_channel; j++){
	psi.MC_tau.n_acc_channel[j] = atoi(readin[proc + 3 * j].c_str());
	psi.MC_tau.n_rej_channel[j] = atoi(readin[proc + 1 + 3 * j].c_str());
	psi.MC_tau.w_channel_step_sum[j] = hexastr2double(readin[proc + 2 + 3 * j]);
      }
      proc += 3 * psi.MC_tau.n_channel;
      
      if (proc != readin.size()){
	// x_alpha_it_min needs to be added here !!!
	psi.MC_tau.i_alpha_it = atoi(readin[proc].c_str());
	psi.MC_tau.x_alpha_it_min = atoi(readin[proc + 1].c_str());
	proc += 2;

	psi.MC_tau.diff_w.resize(n_opt_steps);
	for (int i = 0; i < n_opt_steps; i++){
	  psi.MC_tau.diff_w[i] = hexastr2double(readin[proc + i]);
	}
	proc += n_opt_steps;
	
	if (n_opt_steps != 0){
	  //	proc += n_opt_steps;
	  //	proc = proc + n_opt_steps;
	  // How is this taken into account ??? !!! 
	  psi.MC_tau.alpha_it.resize(n_opt_steps);
	  for (int i = 0; i < n_opt_steps; i++){
	    psi.MC_tau.alpha_it[i].resize(psi.MC_tau.n_channel);
	    for (int j = 0; j < psi.MC_tau.n_channel; j++){
	      psi.MC_tau.alpha_it[i][j] = hexastr2double(readin[proc + psi.MC_tau.n_channel * i + j]);
	    }
	  }
	  //	  proc += n_opt_steps * (psi.MC_tau.n_channel + 1);
	  proc += n_opt_steps * psi.MC_tau.n_channel;
	}
      }
      logger << LOG_INFO << "end MC_tau" << endl;
    }
  }


  logger << LOG_DEBUG << "after MC_tau" << endl;
  for (int j = 0; j < psi.tau_MC_tau_gamma[0].size(); j++){
    logger << LOG_DEBUG << "tau_MC_tau_gamma[0][" << j << "]" << " = " << setprecision(20) << setw(28) << psi.tau_MC_tau_gamma[0][j] << "   " << double2hexastr(psi.tau_MC_tau_gamma[0][j]) << endl;
  }


  if (csi->class_contribution_CS_real){
    for (int i_a = 1; i_a < psi.MC_x_dipole.size(); i_a++){
      if (psi.MC_x_dipole[i_a].active_optimization != -1){
	psi.MC_x_dipole[i_a].end_optimization = atoi(readin[proc].c_str());
	proc += 1;
	for (int j = 0; j < psi.MC_x_dipole[i_a].n_channel; j++){
	  psi.MC_x_dipole[i_a].alpha[j] = hexastr2double(readin[proc + j]);
	}
	proc += psi.MC_x_dipole[i_a].n_channel;
	
	psi.MC_x_dipole[i_a].beta[0] = psi.MC_x_dipole[i_a].alpha[0];
	for (int i = 1; i < psi.MC_x_dipole[i_a].beta.size(); i++){psi.MC_x_dipole[i_a].beta[i] = psi.MC_x_dipole[i_a].beta[i - 1] + psi.MC_x_dipole[i_a].alpha[i];}
	
	//    if (psi_MC_opt_end == 0){
	if (psi.MC_x_dipole[i_a].end_optimization == 0){
	  int n_opt_steps = atoi(readin[proc].c_str());
	  proc += 1;
	  for (int j = 0; j < psi.MC_x_dipole[i_a].n_channel; j++){
	    psi.MC_x_dipole[i_a].n_acc_channel[j] = atoi(readin[proc + 3 * j].c_str());
	    psi.MC_x_dipole[i_a].n_rej_channel[j] = atoi(readin[proc + 1 + 3 * j].c_str());
	    psi.MC_x_dipole[i_a].w_channel_step_sum[j] = hexastr2double(readin[proc + 2 + 3 * j]);
	  }
	  proc += 3 * psi.MC_x_dipole[i_a].n_channel;
	  
	  if (proc != readin.size()){
	    // x_alpha_it_min needs to be added here !!!
	    psi.MC_x_dipole[i_a].i_alpha_it = atoi(readin[proc].c_str());
	    psi.MC_x_dipole[i_a].x_alpha_it_min = atoi(readin[proc + 1].c_str());
	    proc += 2;
	    
	    psi.MC_x_dipole[i_a].diff_w.resize(n_opt_steps);
	    for (int i = 0; i < n_opt_steps; i++){
	      psi.MC_x_dipole[i_a].diff_w[i] = hexastr2double(readin[proc + i]);
	    }
	    proc += n_opt_steps;
	    
	    if (n_opt_steps != 0){
	      //	proc += n_opt_steps;
	      //	proc = proc + n_opt_steps;
	      // How is this taken into account ??? !!! 
	      psi.MC_x_dipole[i_a].alpha_it.resize(n_opt_steps);
	      for (int i = 0; i < n_opt_steps; i++){
		psi.MC_x_dipole[i_a].alpha_it[i].resize(psi.MC_x_dipole[i_a].n_channel);
		for (int j = 0; j < psi.MC_x_dipole[i_a].n_channel; j++){
		  psi.MC_x_dipole[i_a].alpha_it[i][j] = hexastr2double(readin[proc + psi.MC_x_dipole[i_a].n_channel * i + j]);
		}
	      }
	      //	      proc += n_opt_steps * (psi.MC_x_dipole[i_a].n_channel + 1);
	      proc += n_opt_steps * psi.MC_x_dipole[i_a].n_channel;
	    }
	  }
	  logger << LOG_INFO << "end MC_x_dipole[" << i_a << "]" << endl;
	}
      }
    }  
  }


  
  logger << LOG_INFO << "# TSV MC end: proc = " << proc << endl;

  logger << LOG_DEBUG << "end" << endl;
  for (int j = 0; j < psi.tau_MC_tau_gamma[0].size(); j++){
    logger << LOG_DEBUG << "tau_MC_tau_gamma[0][" << j << "]" << " = " << setprecision(20) << setw(28) << psi.tau_MC_tau_gamma[0][j] << "   " << double2hexastr(psi.tau_MC_tau_gamma[0][j]) << endl;
  }


  logger << LOG_DEBUG << "finished" << endl;
}


void observable_set::perform_proceeding_out(phasespace_set & psi){
  Logger logger("observable_set::perform_proceeding_out");
  logger << LOG_DEBUG << "started" << endl;

  logger << LOG_DEBUG << "filename_proceeding = " << filename_proceeding << endl;
  ofstream out_proceeding;
  /*
  string test = filename_proceeding + "x";
  cout << "test = " << test << endl;
  out_proceeding.open(test.c_str(), ofstream::out | ofstream::trunc);  
*/
  out_proceeding.open(filename_proceeding.c_str(), ofstream::out | ofstream::trunc);  
  //  out_proceeding.open(filename_proceeding.c_str(), ofstream::out | ofstream::trunc);  
  out_proceeding << "# random number generator (3)" << endl;
  for (int i = 0; i < 3; i++){
    out_proceeding << double2hexastr(psi_sran[i]) << endl;
    //    out_proceeding << setprecision(25) << psi_sran[i] << endl;
  }
  out_proceeding << "# runtime: h (1) - min (1) - sec (1)" << endl;
  out_proceeding << h << endl;
  out_proceeding << min << endl;
  out_proceeding << sec << endl;

  /* 
  //  old order:
  out_proceeding << psi_i_rej << endl;
  out_proceeding << psi_i_acc << endl;
  out_proceeding << psi_i_nan << endl;
  */
  out_proceeding << "# events: i_gen (1) - i_acc (1) - i_rej (1) - i_tec (1) - i_nan (1)" << endl;
  out_proceeding << psi.i_gen << endl;
  out_proceeding << psi.i_acc << endl;
  out_proceeding << psi.i_rej << endl;
  out_proceeding << psi.i_tec << endl;
  out_proceeding << psi.i_nan << endl;
  // add psi_i_gen and psi_i_tec and change order !!! gen - acc - rej - tec - nan ???

  out_proceeding << "# cross section: full_sum_weight - full_sum_weight2 (1 x 2)" << endl;
  out_proceeding << double2hexastr(full_sum_weight) << endl;
  out_proceeding << double2hexastr(full_sum_weight2) << endl;

  if (switch_CV){
    // !!! Reduce output by putting output_n_qTcut !!!
    out_proceeding << "# cross section CV: full_sum_weight - full_sum_weight2 (" << n_qTcut << " x " << n_scales_CV << " x 2)" << endl;
    for (int i_q = 0; i_q < n_qTcut; i_q++){
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	out_proceeding << double2hexastr(full_sum_weight_CV[i_q][i_s]) << endl;
	out_proceeding << double2hexastr(full_sum_weight2_CV[i_q][i_s]) << endl;
      }
    }
  }
  
  if (n_moments){
    out_proceeding << "# moments: full_sum_moment - full_sum_moment2 (" << n_moments << " x 2)" << endl;
    for (int nm = 0; nm < n_moments; nm++){
      out_proceeding << double2hexastr(full_sum_moment[nm]) << endl;
      out_proceeding << double2hexastr(full_sum_moment2[nm]) << endl;
    }
    if (switch_CV){
      out_proceeding << "# cross section CV: full_sum_moment - full_sum_moment2 (" << n_moments << " x "  << n_qTcut << " x " << n_scales_CV << " x 2)" << endl;
      for (int nm = 0; nm < n_moments; nm++){
	for (int i_q = 0; i_q < n_qTcut; i_q++){
	  for (int i_s = 0; i_s < n_scales_CV; i_s++){
	    out_proceeding << double2hexastr(full_sum_moment_CV[nm][i_q][i_s]) << endl;
	    out_proceeding << double2hexastr(full_sum_moment2_CV[nm][i_q][i_s]) << endl;
	  }
	}
      }
    }
  }

  if (switch_distribution){
    // Could be merged by using 'extended_dat' !!!
    for (int i_d = 0; i_d < dat.size(); i_d++){
      out_proceeding << "# distribution CV: bin_counts - bin_weight - bin_weight2 (" << dat[i_d].n_bins << " x 3 + " << n_scales_CV << " x " << dat[i_d].n_bins << " x 2)" << endl;
      for (int i_n = 0; i_n < dat[i_d].n_bins; i_n++){
	out_proceeding << setprecision(18) << bin_counts[i_d][i_n] << endl;
	out_proceeding << double2hexastr(bin_weight[i_d][i_n]) << endl;
	out_proceeding << double2hexastr(bin_weight2[i_d][i_n]) << endl;
	if (switch_CV){
	  for (int i_s = 0; i_s < n_scales_CV; i_s++){
	    out_proceeding << double2hexastr(bin_weight_CV[i_s][i_d][i_n]) << endl;
	    out_proceeding << double2hexastr(bin_weight2_CV[i_s][i_d][i_n]) << endl;
	  }
	}
      }
    }
    for (int i_ddd = 0; i_ddd < (dddat).size(); i_ddd++){
      out_proceeding << "# dddistribution CV: bin_counts - bin_weight - bin_weight2 (" << dddat[i_ddd].n_bins << " x 3 + " << n_scales_CV << " x " << dddat[i_ddd].n_bins << " x 2)" << endl;
      for (int i_n = 0; i_n < (dddat)[i_ddd].n_bins; i_n++){
	out_proceeding << setprecision(18) << bin_counts[dat.size() + i_ddd][i_n] << endl;
	out_proceeding << double2hexastr(bin_weight[dat.size() + i_ddd][i_n]) << endl;
	out_proceeding << double2hexastr(bin_weight2[dat.size() + i_ddd][i_n]) << endl;
	if (switch_CV){
	  for (int i_s = 0; i_s < n_scales_CV; i_s++){
	    out_proceeding << double2hexastr(bin_weight_CV[i_s][dat.size() + i_ddd][i_n]) << endl;
	    out_proceeding << double2hexastr(bin_weight2_CV[i_s][dat.size() + i_ddd][i_n]) << endl;
	  }
	}
      }
    }
  }

  if (switch_TSV){
    out_proceeding << "# TSV output" << endl;
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      out_proceeding << "# cross section TSV (" << name_set_TSV[i_s] << "): fullsum_weight_TSV - fullsum_weight2_TSV (" << n_scale_fact_TSV[i_s] << " x " << n_scale_ren_TSV[i_s] << " x 2)" << endl;
      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	  out_proceeding << double2hexastr(fullsum_weight_TSV[i_s][i_r][i_f]) << endl;
	  out_proceeding << double2hexastr(fullsum_weight2_TSV[i_s][i_r][i_f]) << endl;
	}
      }
      if (active_qTcut == 1){
      out_proceeding << "# cross section qTcut TSV (" << name_set_TSV[i_s] << "): fullsum_weight_TSV - fullsum_weight2_TSV (" << n_qTcut << " x " << n_scale_fact_TSV[i_s] << " x " << n_scale_ren_TSV[i_s] << " x 2)" << endl;
	for (int i_q = 0; i_q < n_qTcut; i_q++){
	  for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	      out_proceeding << double2hexastr(fullsum_weight_qTcut_TSV[i_s][i_q][i_r][i_f]) << endl;
	      out_proceeding << double2hexastr(fullsum_weight2_qTcut_TSV[i_s][i_q][i_r][i_f]) << endl;
	    }
	  }
	}
      }
    }
    
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      out_proceeding << "# moments TSV (" << name_set_TSV[i_s] << "): fullsum_weight_TSV - fullsum_weight2_TSV (" << n_moments << " x " << n_scale_fact_TSV[i_s] << " x " << n_scale_ren_TSV[i_s] << " x 2)" << endl;
      for (int i_m = 0; i_m < n_moments; i_m++){
	for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	    out_proceeding << double2hexastr(fullsum_moment_TSV[i_s][i_m][i_r][i_f]) << endl;
	    out_proceeding << double2hexastr(fullsum_moment2_TSV[i_s][i_m][i_r][i_f]) << endl;
	  }
	}
      }
    }
    //    for (int i_q = 0; i_q < 1; i_q++){
    int next_no_qTcut = 0;
    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
      //    for (int i_q = 0; i_q < n_qTcut; i_q++){
      if (next_no_qTcut == no_qTcut_distribution.size()){continue;}
      if (i_q != no_qTcut_distribution[next_no_qTcut]){continue;}

      if (active_qTcut == 0){out_proceeding << "#   no qTcut" << endl;}
      else {out_proceeding << "#   qTcut[" << setw(3) << i_q << "] = " << setw(15) << setprecision(8) << value_qTcut[i_q] << endl;}

      for (int i_d = 0; i_d < dat.size(); i_d++){
	out_proceeding << "#" << char(9) << dat[i_d].xdistribution_name << endl;
	for (int i_b = 0; i_b < dat[i_d].n_bins; i_b++){
	  out_proceeding << bin_count_TSV[i_q][i_d][i_b] << endl;
	  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		out_proceeding << double2hexastr(bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][i_b]) << endl;
		out_proceeding << double2hexastr(bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][i_b]) << endl;
	      }
	    }
	  }
	}
      }
      for (int i_d = 0; i_d < dddat.size(); i_d++){
	int i_ddd = dat.size() + i_d;
	out_proceeding << "#" << char(9) << dddat[i_d].name << endl;
	for (int i_b = 0; i_b < dddat[i_d].n_bins; i_b++){
	  out_proceeding << bin_count_TSV[i_q][i_ddd][i_b] << endl;
	  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		out_proceeding << double2hexastr(bin_weight_TSV[i_s][i_q][i_r][i_f][i_ddd][i_b]) << endl;
		out_proceeding << double2hexastr(bin_weight2_TSV[i_s][i_q][i_r][i_f][i_ddd][i_b]) << endl;
	      }
	    }
	  }
	}
      }

      next_no_qTcut++;
    }
  }

  // One could avoid all the following output for the case of a completely finalized optimization phase (like in all standard non-grid runs) !!!
  // Introduce a general switch like e.g. 'opt_end', and the following is only written if that switch is 0 (or the other way round) !!!
  
  out_proceeding << "# IS_tau - end_optimization (1) - alpha (" << psi.IS_tau.alpha.size() << ")" << endl;
  out_proceeding << psi.IS_tau.end_optimization << endl;
  if (psi.IS_tau.end_optimization != -1){
    for (int j = 0; j < psi.IS_tau.alpha.size(); j++){
      out_proceeding << double2hexastr(psi.IS_tau.alpha[j]) << endl;
    }
  }
  if (psi.IS_tau.end_optimization == 0){
    out_proceeding << "# IS_tau - n_acc_channel - n_rej_channel - sum_channel_weight - sum_channel_weight2 (" << psi.IS_tau.alpha.size() << " x 4)" << endl;
    for (int j = 0; j < psi.IS_tau.alpha.size(); j++){
      out_proceeding << psi.IS_tau.n_acc_channel[j] << endl;
      out_proceeding << psi.IS_tau.n_rej_channel[j] << endl;
      out_proceeding << double2hexastr(psi.IS_tau.sum_channel_weight[j]) << endl;
      out_proceeding << double2hexastr(psi.IS_tau.sum_channel_weight2[j]) << endl;
    }
  }
  /*
  out_proceeding << "# tau_opt output" << endl;
  out_proceeding << psi_tau_opt_end << endl;
  if (psi_tau_opt_end != -1){
    for (int j = 0; j < psi_tau_alpha.size(); j++){
      out_proceeding << double2hexastr(psi_tau_alpha[j]) << endl;
    }
  }
  if (psi_tau_opt_end == 0){
    for (int j = 0; j < psi_tau_alpha.size(); j++){
      out_proceeding << psi_tau_counts_channel[j] << endl;
      out_proceeding << psi_tau_cuts_channel[j] << endl;
      out_proceeding << double2hexastr(psi_sum_tau_channel_weight[j]) << endl;
      out_proceeding << double2hexastr(psi_sum_tau_channel_weight2[j]) << endl;
    }
  }
  */

  out_proceeding << "# IS_x1x2 - end_optimization (1) - alpha (" << psi.IS_x1x2.alpha.size() << ")" << endl;
  out_proceeding << psi.IS_x1x2.end_optimization << endl;
  if (psi.IS_x1x2.end_optimization != -1){
    for (int j = 0; j < psi.IS_x1x2.alpha.size(); j++){
      out_proceeding << double2hexastr(psi.IS_x1x2.alpha[j]) << endl;
    }
  }
  if (psi.IS_x1x2.end_optimization == 0){
    out_proceeding << "# IS_x1x2 - n_acc_channel - n_rej_channel - sum_channel_weight - sum_channel_weight2 (" << psi.IS_x1x2.alpha.size() << " x 4)" << endl;
    for (int j = 0; j < psi.IS_x1x2.alpha.size(); j++){
      out_proceeding << psi.IS_x1x2.n_acc_channel[j] << endl;
      out_proceeding << psi.IS_x1x2.n_rej_channel[j] << endl;
      out_proceeding << double2hexastr(psi.IS_x1x2.sum_channel_weight[j]) << endl;
      out_proceeding << double2hexastr(psi.IS_x1x2.sum_channel_weight2[j]) << endl;
    }
  }

  /*
  out_proceeding << "# x1x2_opt output" << endl;
  out_proceeding << psi_x1x2_opt_end << endl;
  if (psi_x1x2_opt_end != -1){
    for (int j = 0; j < psi_x1x2_alpha.size(); j++){
      out_proceeding << double2hexastr(psi_x1x2_alpha[j]) << endl;
    }
  }
  if (psi_x1x2_opt_end == 0){
    for (int j = 0; j < psi_x1x2_alpha.size(); j++){
      out_proceeding << psi_x1x2_counts_channel[j] << endl;
      out_proceeding << psi_x1x2_cuts_channel[j] << endl;
      out_proceeding << double2hexastr(psi_sum_x1x2_channel_weight[j]) << endl;
      out_proceeding << double2hexastr(psi_sum_x1x2_channel_weight2[j]) << endl;
    }
  }
  */
  
  out_proceeding << "# IS_z1z2.end_optimization output" << endl;

  logger << LOG_INFO << "csi->class_contribution_CS_collinear = " << csi->class_contribution_CS_collinear << endl;
  
  if (csi->class_contribution_CS_collinear){
    //  only for backwards-compatibility reasons -> can be removed !!!
    out_proceeding << psi.IS_z1z2[1].end_optimization << endl;
    // later...    out_proceeding << psi.IS_z1z2.end_optimization << endl;
    for (int i_z = 1; i_z < 3; i_z++){
      out_proceeding << "# IS_z1z2[" << i_z << "] - end_optimization - alpha (" << psi.IS_z1z2[i_z].alpha.size() << " x 4)" << endl;
      //  add this line !!!
      //      out_proceeding << psi.IS_z1z2[i_z].end_optimization << endl;
      if (psi.IS_z1z2[i_z].end_optimization != -1){
	//      if (psi_z1z2_opt_end != -1){
	for (int j = 0; j < psi.IS_z1z2[i_z].alpha.size(); j++){
	  out_proceeding << double2hexastr(psi.IS_z1z2[i_z].alpha[j]) << endl;
	}
	if (psi.IS_z1z2[i_z].end_optimization == 0){
	  for (int j = 0; j < psi.IS_z1z2[i_z].alpha.size(); j++){
	    out_proceeding << psi.IS_z1z2[i_z].n_acc_channel[j] << endl;
	    out_proceeding << psi.IS_z1z2[i_z].n_rej_channel[j] << endl;
	    out_proceeding << double2hexastr(psi.IS_z1z2[i_z].sum_channel_weight[j]) << endl;
	    out_proceeding << double2hexastr(psi.IS_z1z2[i_z].sum_channel_weight2[j]) << endl;
	  }
	}
      }
    }
  }
  else {
    // only for backwards-compatibility reasons -> can be removed !!!
    out_proceeding << "-1" << endl;
  }
  /*
  out_proceeding << "# z1z2_opt output" << endl;

  out_proceeding << psi_z1z2_opt_end << endl;

  if (psi_z1z2_opt_end != -1){
    for (int i_z = 1; i_z < 3; i_z++){
      //      if (psi_z1z2_opt_end != -1){
      for (int j = 0; j < psi_z1z2_alpha[i_z].size(); j++){
	out_proceeding << double2hexastr(psi_z1z2_alpha[i_z][j]) << endl;
      }
      if (psi_z1z2_opt_end == 0){
	for (int j = 0; j < psi_z1z2_alpha[i_z].size(); j++){
	  out_proceeding << psi_z1z2_counts_channel[i_z][j] << endl;
	  out_proceeding << psi_z1z2_cuts_channel[i_z][j] << endl;
	  out_proceeding << double2hexastr(psi_sum_z1z2_channel_weight[i_z][j]) << endl;
	  out_proceeding << double2hexastr(psi_sum_z1z2_channel_weight2[i_z][j]) << endl;
	}
      }
    }
  }
  */

  
  out_proceeding << "# random_manager" << endl;
  psi_random_manager.proceeding_out(out_proceeding);

  // Solve issues that MC_phasespace does not work correctly after resumption !!!
  // Extend for MC_tau and MC_x_dipole !!!

  out_proceeding << "# MC_opt output" << endl;
  out_proceeding << psi_MC_opt_end << endl;

  if (psi.MC_phasespace.active_optimization != -1){
    out_proceeding << "# MC_phasespace - end_optimization (1) - alpha (" << psi.MC_phasespace.alpha.size() << ")" << endl;
    out_proceeding << psi.MC_phasespace.end_optimization << endl;
    ///  out_proceeding << psi.MC_phasespace.alpha.size() << endl;
    for (int j = 0; j < psi.MC_phasespace.alpha.size(); j++){
      out_proceeding << double2hexastr(psi.MC_phasespace.alpha[j]) << endl;
    }
    //  if (psi_MC_opt_end == 0){
    if (psi.MC_phasespace.end_optimization == 0){
      int n_opt_steps = psi.MC_phasespace.alpha_it.size();
      out_proceeding << "# MC_phasespace - n_opt_steps (1) - n_acc_channel - n_rej_channel - w_sum_channel_weight (3 x " << psi.MC_phasespace.alpha.size() << ")" << endl;
      out_proceeding << n_opt_steps << endl;
      for (int j = 0; j < psi.MC_phasespace.alpha.size(); j++){
	//     out_proceeding << psi_counts_channel[j] << endl;
	//      out_proceeding << psi_cuts_channel[j] << endl;
	out_proceeding << psi.MC_phasespace.n_acc_channel[j] << endl;
	out_proceeding << psi.MC_phasespace.n_rej_channel[j] << endl;
	out_proceeding << double2hexastr(psi.MC_phasespace.w_channel_step_sum[j]) << endl;
      }
      // x_alpha_it_min needs to be added here !!!
      out_proceeding << "# MC_phasespace - i_alpha_it (1) - x_alpha_it_min (1) - diff_w (" << n_opt_steps << ")" << endl;
      
      out_proceeding << setprecision(18) << psi.MC_phasespace.i_alpha_it << endl;
      out_proceeding << setprecision(18) << psi.MC_phasespace.x_alpha_it_min << endl;
      for (int i = 0; i < n_opt_steps; i++){
	out_proceeding << double2hexastr(psi.MC_phasespace.diff_w[i]) << endl;
      }
      out_proceeding << "# MC_phasespace - alpha_it (" << n_opt_steps << " x " << psi.MC_phasespace.alpha.size() << ")" << endl;
      for (int i = 0; i < n_opt_steps; i++){
	for (int j = 0; j < psi.MC_phasespace.alpha.size(); j++){
	  out_proceeding << double2hexastr(psi.MC_phasespace.alpha_it[i][j]) << endl;
	}
      }
      //      out_proceeding << "# MC_phasespace - weight optimization not yet finished!" << endl;
    }
  }
  else {
    out_proceeding << "# MC_phasespace - no_optimization due to active_optimization = " << psi.MC_phasespace.active_optimization << endl;
  }


  
  if (psi.MC_tau.active_optimization != -1){
    out_proceeding << "# MC_tau - end_optimization - alpha (" << psi.MC_tau.alpha.size() << ")" << endl;
    out_proceeding << psi.MC_tau.end_optimization << endl;
    ///  out_proceeding << psi.MC_tau.alpha.size() << endl;
    for (int j = 0; j < psi.MC_tau.alpha.size(); j++){
      out_proceeding << double2hexastr(psi.MC_tau.alpha[j]) << endl;
    }
    //  if (psi_MC_opt_end == 0){
    if (psi.MC_tau.end_optimization == 0){
      int n_opt_steps = psi.MC_tau.alpha_it.size();
      out_proceeding << n_opt_steps << endl;
      for (int j = 0; j < psi.MC_tau.alpha.size(); j++){
	//     out_proceeding << psi_counts_channel[j] << endl;
	//      out_proceeding << psi_cuts_channel[j] << endl;
	out_proceeding << psi.MC_tau.n_acc_channel[j] << endl;
	out_proceeding << psi.MC_tau.n_rej_channel[j] << endl;
	out_proceeding << double2hexastr(psi.MC_tau.w_channel_step_sum[j]) << endl;
      }
      // x_alpha_it_min needs to be added here !!!
      out_proceeding << setprecision(18) << psi.MC_tau.i_alpha_it << endl;
      out_proceeding << setprecision(18) << psi.MC_tau.x_alpha_it_min << endl;
      for (int i = 0; i < n_opt_steps; i++){
	out_proceeding << double2hexastr(psi.MC_tau.diff_w[i]) << endl;
      }
      for (int i = 0; i < n_opt_steps; i++){
	for (int j = 0; j < psi.MC_tau.alpha.size(); j++){
	  out_proceeding << double2hexastr(psi.MC_tau.alpha_it[i][j]) << endl;
	}
      }
      out_proceeding << "# MC_tau - weight optimization not yet finished!" << endl;
    }
  }
  else {
    out_proceeding << "# MC_tau - no_optimization due to active_optimization = " << psi.MC_tau.active_optimization << endl;
  }

  if (csi->class_contribution_CS_real){
    for (int i_a = 1; i_a < psi.MC_x_dipole.size(); i_a++){
      if (psi.MC_x_dipole[i_a].active_optimization != -1){
	out_proceeding << "# MC_x_dipole[" << i_a << "] - end_optimization - alpha (" << psi.MC_x_dipole[i_a].alpha.size() << ")" << endl;
	out_proceeding << psi.MC_x_dipole[i_a].end_optimization << endl;
	///  out_proceeding << psi.MC_x_dipole[i_a].alpha.size() << endl;
	for (int j = 0; j < psi.MC_x_dipole[i_a].alpha.size(); j++){
	  out_proceeding << double2hexastr(psi.MC_x_dipole[i_a].alpha[j]) << endl;
	}
	//  if (psi_MC_opt_end == 0){
	if (psi.MC_x_dipole[i_a].end_optimization == 0){
	  int n_opt_steps = psi.MC_x_dipole[i_a].alpha_it.size();
	  out_proceeding << n_opt_steps << endl;
	  for (int j = 0; j < psi.MC_x_dipole[i_a].alpha.size(); j++){
	    //     out_proceeding << psi_counts_channel[j] << endl;
	    //      out_proceeding << psi_cuts_channel[j] << endl;
	    out_proceeding << psi.MC_x_dipole[i_a].n_acc_channel[j] << endl;
	    out_proceeding << psi.MC_x_dipole[i_a].n_rej_channel[j] << endl;
	    out_proceeding << double2hexastr(psi.MC_x_dipole[i_a].w_channel_step_sum[j]) << endl;
	  }
	  // x_alpha_it_min needs to be added here !!!
	  out_proceeding << setprecision(18) << psi.MC_x_dipole[i_a].i_alpha_it << endl;
	  out_proceeding << setprecision(18) << psi.MC_x_dipole[i_a].x_alpha_it_min << endl;
	  for (int i = 0; i < n_opt_steps; i++){
	    out_proceeding << double2hexastr(psi.MC_x_dipole[i_a].diff_w[i]) << endl;
	  }
	  for (int i = 0; i < n_opt_steps; i++){
	    for (int j = 0; j < psi.MC_x_dipole[i_a].alpha.size(); j++){
	      out_proceeding << double2hexastr(psi.MC_x_dipole[i_a].alpha_it[i][j]) << endl;
	    }
	  }
	  out_proceeding << "# MC_x_dipole[" << i_a << "] - weight optimization not yet finished!" << endl;
	}
      }
      else {
	out_proceeding << "# MC_x_dipole[i_a] - no_optimization due to active_optimization = " << psi.MC_x_dipole[i_a].active_optimization << endl;
      }
      
    }
  }
      
  out_proceeding.close();

  logger << LOG_DEBUG << "finished" << endl;
}


void observable_set::perform_proceeding_check(phasespace_set & psi){
  Logger logger("observable_set::perform_proceeding_check");
  logger << LOG_DEBUG << "started" << endl;

  double Xsection_delta = sqrt((full_sum_weight2 - pow(full_sum_weight, 2) / psi_i_gen)) / (psi_i_gen - 1);
  double Xsection = full_sum_weight / psi_i_gen;
  double weights2_NLO_LO = abs(Xsection_delta / sigma_normalization);

  logger << LOG_INFO << "Xsection_delta = " << Xsection_delta << endl;
  logger << LOG_INFO << "Xsection = " << Xsection << endl;
  logger << LOG_INFO << "weights2_NLO_LO = " << weights2_NLO_LO << endl;
  logger << LOG_INFO << "sigma_LO = " << sigma_normalization << endl;
  logger << LOG_INFO << "n_events_max = " << psi_n_events_max << endl;
  logger << LOG_INFO << "n_step = " << psi_n_step << endl;

  ///  logger << LOG_INFO << "jets[1] = " << jets[1] << endl;
  ///  logger << LOG_INFO << "sigma_LO_deviation = " << sigma_normalization_deviation << endl;
  ///  logger << LOG_INFO << "(jets[1] >= n_events_max - n_step) = " << (jets[1] >= psi_n_events_max - psi_n_step) << endl;
  ///  logger << LOG_INFO << "(jets[1] >= n_events_min && weights2_NLO_LO < sigma_LO_deviation) = " << (jets[1] >= psi_n_events_min && weights2_NLO_LO < sigma_normalization_deviation) << endl;

  ///  if ((jets[1] >= psi_n_events_max - psi_n_step) || 
  ///      (jets[1] >= psi_n_events_min && weights2_NLO_LO < sigma_normalization_deviation)){int_end = 1;}
  logger << LOG_INFO << "psi_i_acc = " << psi_i_acc << endl;
  logger << LOG_INFO << "sigma_LO_deviation = " << sigma_normalization_deviation << endl;
  logger << LOG_INFO << "switch_output_proceeding = " << switch_output_proceeding << endl;
  logger << LOG_INFO << "(psi_i_acc >= n_events_max) = " << (psi_i_acc >= psi_n_events_max) << endl;
  logger << LOG_INFO << "(psi_i_acc >= n_events_max - n_step) = " << (psi_i_acc >= psi_n_events_max - psi_n_step) << endl;
  logger << LOG_INFO << "(psi_i_acc >= n_events_min && weights2_NLO_LO < sigma_LO_deviation) = " << (psi_i_acc >= psi_n_events_min && weights2_NLO_LO < sigma_normalization_deviation) << endl;

  if (((psi_i_acc >= psi_n_events_max) && switch_output_proceeding == 1) || 
      ((psi_i_acc >= psi_n_events_max - psi_n_step) && switch_output_proceeding == 2) || 
      ((psi_i_acc >= psi_n_events_min) && (weights2_NLO_LO < sigma_normalization_deviation) && (psi_switch_n_events_opt != 2))){int_end = 1;}
  else {
	//      out_integration << "******************************************************* resumption ******************************************************" << endl;
  }

  logger << LOG_INFO << "int_end = " << int_end << endl;

  logger << LOG_DEBUG << "finished" << endl;
}
