#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_contribution::readin_result_contribution_CV(string result_moment, int x_m){
  Logger logger("summary_contribution::readin_result_contribution_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  if (yresultdirectory != ylist->resultdirectory){logger << LOG_INFO << "yresultdirectory = " << yresultdirectory << "  =/=  " << ylist->resultdirectory << " =  ylist->resultdirectory" << endl;}

  //  int n_subgroup = ygeneric->subgroup.size();
  //  string result_file;
  //  char LineBuffer[128];

  logger.newLine(LOG_INFO);
  logger << LOG_DEBUG << "readin_result_contribution_CV started:" << endl;
  logger << LOG_DEBUG << "ylist->resultdirectory = " << ylist->resultdirectory << endl;
  logger << LOG_DEBUG << "type_perturbative_order = " << type_perturbative_order << endl;
  logger << LOG_DEBUG << "type_subtraction_method = " << type_subtraction_method << endl;
  logger << LOG_DEBUG << "in_contribution_order_alpha_s = " << in_contribution_order_alpha_s << endl;
  logger << LOG_DEBUG << "in_contribution_order_alpha_e = " << in_contribution_order_alpha_e << endl;
  logger << LOG_DEBUG << "type_contribution = " << type_contribution << endl;
  logger << LOG_DEBUG << "type_correction = " << type_correction << endl;
  logger << LOG_DEBUG << "interference = " << interference << endl;
  logger << LOG_DEBUG << "photon_induced = " << photon_induced << endl;
  logger << LOG_DEBUG << "infix_contribution = " << infix_contribution << endl;

  logger << LOG_INFO << "infix_order_contribution = " << infix_order_contribution << endl;

  for (int i_d = 0; i_d < directory.size(); i_d++){logger << LOG_DEBUG << "directory[" << setw(3) << "] = " << directory[i_d] << endl;}
  logger << LOG_DEBUG << "xsubprocess.size() = " << xsubprocess.size() << endl;
  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){logger << LOG_DEBUG << "xsubprocess[" << setw(3) << "] = " << xsubprocess[i_p].name << endl;}
  logger << LOG_DEBUG << "ygeneric->subgroup.size() = " << ygeneric->subgroup.size() << endl;
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG << "ygeneric->subgroup[" << setw(3) << "] = " << ygeneric->subgroup[i_g] << endl;}
  logger << LOG_DEBUG << "subgroup_no_member = " << subgroup_no_member.size() << endl;
  logger << LOG_DEBUG << "ygeneric->subgroup.size() = " << ygeneric->subgroup.size() << endl;
  logger << LOG_DEBUG << "result_moment = " << result_moment << endl;
  /*
  //  int active_qTcut = 0;
  int output_n_qTcut = 0;
  if (type_contribution == "CT" ||
      type_contribution == "RT" ||
      type_contribution == "CT2" ||
      type_contribution == "RVA" ||
      type_contribution == "RCA" ||
      type_contribution == "RRA"){
    osi->active_qTcut = 1;
    output_n_qTcut = osi->n_qTcut;
  }
  else {
    osi->active_qTcut = 0;
    output_n_qTcut = 1;
  }
  */
  /*
  string name_variation_CV;
  if      (osi->switch_CV == 1){name_variation_CV = "equal";}
  else if (osi->switch_CV == 2){name_variation_CV = "ren";}
  else if (osi->switch_CV == 3){name_variation_CV = "fac";}
  else if (osi->switch_CV == 4){name_variation_CV = "antipodal";}
  else if (osi->switch_CV == 5){name_variation_CV = "7-point";}
  else if (osi->switch_CV == 6){name_variation_CV = "9-point";}
  */
  //  string readindata;

  logger << LOG_DEBUG_VERBOSE << "osi->n_qTcut        = " << osi->n_qTcut << endl;
  logger << LOG_DEBUG_VERBOSE << "osi->min_qTcut     = " << osi->min_qTcut << endl;
  logger << LOG_DEBUG_VERBOSE << "osi->value_qTcut.size()     = " << osi->value_qTcut.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "osi->n_scales_CV    = " << osi->n_scales_CV << endl;
  //osi->value_qTcut



  int default_size_CV =  1 + 2;
  /*
  logger << LOG_DEBUG_VERBOSE << "osi->active_qTcut    = " << osi->active_qTcut << endl;
  if (type_contribution == "CT" ||
      type_contribution == "RT" ||
      type_contribution == "CT2" ||
      type_contribution == "RVA" ||
      type_contribution == "RCA" ||
      type_contribution == "RRA"){
    osi->active_qTcut = 1;
  }
  else {
    osi->active_qTcut = 0;
  }
  */
  //  logger << LOG_DEBUG_VERBOSE << "osi->active_qTcut    = " << osi->active_qTcut << endl;
  //  if (osi->active_qTcut == 0){default_size_CV += 2 * osi->n_scales_CV;}
  //  else {default_size_CV += 2 * osi->n_qTcut * osi->n_scales_CV;}
  if (active_qTcut){default_size_CV += 2 * osi->n_qTcut * osi->n_scales_CV;}
  else {default_size_CV += 2 * osi->n_scales_CV;}

  logger << LOG_DEBUG_VERBOSE << "default_size_CV = " << default_size_CV << endl;



  result_CV.resize(osi->n_moments + 1, vector<vector<vector<double> > > (ygeneric->subgroup.size(), vector<vector<double> > (osi->n_qTcut, vector<double> (osi->n_scales_CV, 0.))));
  deviation_CV.resize(osi->n_moments + 1, vector<vector<vector<double> > > (ygeneric->subgroup.size(), vector<vector<double> > (osi->n_qTcut, vector<double> (osi->n_scales_CV, 0.))));




  //  int n_xsubprocess = xsubprocess.size();

  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    xsubprocess[i_p].readin_result_CV(result_moment);
    logger << LOG_DEBUG << "xsubprocess[" << i_p << "].name = " << xsubprocess[i_p].name << "   in1 = " << xsubprocess[i_p].type_parton[0][1] << "   in2 = " << xsubprocess[i_p].type_parton[0][2] << endl;
  }

  ///////////////////////////////////////////////////////////////
  //  combination of different-seed runs for all subprocesses  //
  ///////////////////////////////////////////////////////////////
  //  double temp_singleweight;
  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    logger << LOG_DEBUG_VERBOSE << left << xsubprocess[i_p].name << ":" << endl;
    xsubprocess[i_p].combine_result_alternative();
    //    xsubprocess[i_p].combine_result_original(oset);
  }

  /////////////////////////////////////////////////////////////////////////////////
  //  combination of different-seed runs for all subprocesses with CV variation  //
  /////////////////////////////////////////////////////////////////////////////////
  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    logger << LOG_DEBUG_VERBOSE << left << xsubprocess[i_p].name << ":" << endl;
    xsubprocess[i_p].combination_result_CV(result_moment);
    /*
    if (average_factor == 0){xsubprocess[i_p].combine_result_alternative_CV();}
    else if (average_factor == 1){xsubprocess[i_p].combine_result_original_CV();}
    else {xsubprocess[i_p].combine_result_hybrid_CV();}
    */
    /*
    vector<vector<vector<double> > > ().swap(xsubprocess[i_p].result_run_CV);
    vector<vector<vector<double> > > ().swap(xsubprocess[i_p].deviation_run_CV);
    */

  }

  for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      xsubprocess[0].result_CV[i_q][i_s] = 0.;
      double dev = 0.;
      for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
	xsubprocess[0].result_CV[i_q][i_s] += xsubprocess[i_p].result_CV[i_q][i_s];
	dev += pow(xsubprocess[i_p].deviation_CV[i_q][i_s], 2);
      }
      xsubprocess[0].deviation_CV[i_q][i_s] = sqrt(dev);
    }
  }


  
  /////////////////////////////////////////////////////////////////////////////////////////
  //  calculation of cross sections for subprocess groups (i_g = 0 -> all subprocesses)  //
  /////////////////////////////////////////////////////////////////////////////////////////
  /*
  vector<string> latex_subgroup = ygeneric->subgroup;
  get_latex_subgroup(latex_subgroup, ygeneric->subgroup);
  */

  /*
  double sub_result[ygeneric->subgroup.size()];
  double sub_deviation[ygeneric->subgroup.size()];
  */

  //  vector<vector<vector<double> > > result_subgroup_CV(ygeneric->subgroup.size(), vector<vector<double> > (osi->n_qTcut, vector<double> (osi->n_scales_CV)));
  //  vector<vector<vector<double> > > deviation_subgroup_CV(ygeneric->subgroup.size(), vector<vector<double> > (osi->n_qTcut, vector<double> (osi->n_scales_CV)));

  ///////////////////////////////////////////////////////
  //  calculation of subgroup results in CV variation  //
  ///////////////////////////////////////////////////////
    for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	  double dev = 0.;
	  for (int i_m = 0; i_m < subgroup_no_member[i_g].size(); i_m++){
	    /*
	    logger << LOG_DEBUG << "i_q[" << setw(2) << right << i_q << "] i_s[" << setw(2) << right << i_s << "] i_g = [" << setw(2) << right << i_g << "]  i_m = [" << setw(2) << right << i_m << "]" << endl;
	    logger << LOG_DEBUG << "x_m[" << setw(2) << right << x_m << "] x_c[" << setw(2) << right << x_c << endl;
	    logger << LOG_DEBUG << "ygeneric->subgroup.size() = " << ygeneric->subgroup.size() << endl;
	    logger << LOG_DEBUG << "subgroup_no_member[" << i_g << "].size() = " << subgroup_no_member[i_g].size() << endl;
	    logger << LOG_DEBUG << "subgroup_no_member[" << i_g << "][" << i_m << "] = " << subgroup_no_member[i_g][i_m] << endl;
	    */
	    result_CV[x_m][i_g][i_q][i_s] += xsubprocess[subgroup_no_member[i_g][i_m]].result_CV[i_q][i_s];
	    dev += pow(xsubprocess[subgroup_no_member[i_g][i_m]].deviation_CV[i_q][i_s], 2);
	  }
	  deviation_CV[x_m][i_g][i_q][i_s] = sqrt(dev);
	}
      }
    }
  logger << LOG_DEBUG_VERBOSE << "x_m = " << x_m << endl;
  logger << LOG_DEBUG_VERBOSE << "result_CV.size()    = " << setw(3) << result_CV.size() << "   deviation_CV.size() = " << setw(3) << deviation_CV.size() << endl;


  ////////////////////
  //  output begin  //
  ////////////////////

  if ((x_m == 0 || x_m == osi->n_moments)){// || (infix_order_contribution == "LO")){
    /*
    logger << LOG_DEBUG << "xsubprocess.size() = " << xsubprocess.size() << endl;
    //    logger << LOG_DEBUG << "result_subprocess.size() = " << result_subprocess.size() << endl;
    //    logger << LOG_DEBUG << "deviation_subprocess.size() = " << deviation_subprocess.size() << endl;

    ofstream outfile_check;
    string filename;
    string signame = "XS_" + infix_order_contribution;
    string dsigname = "dXS_" + infix_order_contribution;
    
    //////////////////////////////////////////
    //  subprocess-wise check output begin  //
    //////////////////////////////////////////
    filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/check." + result_moment + "." + infix_order_contribution + ".dat";
    // + infix_order_contribution + "
    outfile_check.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    outfile_check << showpoint;
    for (int i_p = 1; i_p < subprocess.size(); i_p++){
      outfile_check << subprocess[i_p] << endl << endl;
      outfile_check << setw(15) << setprecision(8) << osi->unit_factor_result * xsubprocess[i_p].result << setw(15) << setprecision(8) << osi->unit_factor_result * xsubprocess[i_p].deviation << endl;

      logger << LOG_DEBUG << "result_subprocess_CV[" << i_p << "][0][" << osi->no_central_scale_CV << "] = " << xsubprocess[i_p].result_CV[0][osi->no_central_scale_CV] << endl;
      if (osi->n_qTcut != 0){outfile_check << setw(15) << setprecision(8) << osi->unit_factor_result * xsubprocess[i_p].result_CV[0][osi->no_central_scale_CV] << setw(15) << setprecision(8) << osi->unit_factor_result * xsubprocess[i_p].deviation_CV[0][osi->no_central_scale_CV] << setw(8) << setprecision(2) << endl;}
      outfile_check << endl;
    }
    stringstream muscale_ss;
    muscale_ss << setprecision(6) << osi->scale_ren;
    string muscale = muscale_ss.str();
    
    vector<double> rel_scale_CV(osi->n_scales_CV);
    vector<double> scale_CV(osi->n_scales_CV);
    logger << LOG_DEBUG << "osi->n_scales_CV = " << osi->n_scales_CV << endl;
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      if (osi->n_scales_CV > 1){rel_scale_CV[i_s] = pow(10,log10(double(osi->variation_factor_CV)) * double(2. * i_s / (osi->n_scales_CV - 1) - 1));}
      else {rel_scale_CV[i_s] = 1;}
      scale_CV[i_s] = rel_scale_CV[i_s] * osi->scale_ren;
    }
    */

    //    if (ygeneric->switch_output_subprocess > 1){output_subprocesses_result_plot_qTcut_CV();}
    /*
    if (ygeneric->switch_output_subprocess > 1){
      ofstream outfile_subprocess;
      ///////////////////////////////////////////////////////////////////
      //  subprocess qTcut output for CV variation (all qTcut values)  //
      ///////////////////////////////////////////////////////////////////
      for (int i_p = 1; i_p < subprocess.size(); i_p++){
	if (x_m == osi->n_moments){filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/subprocesses/plot.qTcut." + subprocess[i_p] + ".dat";}
	else {filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/subprocesses/plot." + result_moment + ".qTcut." + subprocess[i_p] + ".dat";}
	logger << LOG_DEBUG << "Output: xfilename =  " << filename << endl;
	outfile_subprocess.open(filename.c_str(), ofstream::out | ofstream::trunc);  
	//      outfile_subprocess << showpoint;
	for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	  outfile_subprocess << noshowpoint << setw(10) << setprecision(4) << osi->value_qTcut[i_q];
	  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	    outfile_subprocess << showpoint << setw(15) << setprecision(8) << osi->unit_factor_result * xsubprocess[i_p].result_CV[i_q][i_s] << setw(15) << setprecision(8) << osi->unit_factor_result * xsubprocess[i_p].deviation_CV[i_q][i_s];
	  }
	  outfile_subprocess << endl;
	}
	outfile_subprocess.close();
      }
    }
    */

    //    if (ygeneric->switch_output_subprocess > 1){output_subprocesses_result_plot_CV();}
    /*
    if (ygeneric->switch_output_subprocess > 1){
      ofstream outfile_subprocess;
      ////////////////////////////////////////////////////////////////
      //  subprocess output for CV variation (osi->min_qTcut value)  //
      ////////////////////////////////////////////////////////////////
      logger.newLine(LOG_DEBUG);
      logger << LOG_DEBUG << "subprocess output for CV variation (osi->min_qTcut value)" << endl;
      for (int i_p = 1; i_p < subprocess.size(); i_p++){
	if (x_m == osi->n_moments){filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/subprocesses/plot.CV." + subprocess[i_p] + ".dat";}
	else {filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/subprocesses/plot." + result_moment + ".CV." + subprocess[i_p] + ".dat";}
	logger << LOG_DEBUG << "Output: xfilename = " << filename << endl;
	outfile_subprocess.open(filename.c_str(), ofstream::out | ofstream::trunc);  
	outfile_subprocess << showpoint;
	int i_q = 0;
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  outfile_subprocess << setw(15) << setprecision(8) << rel_scale_CV[i_s] << setw(15) << setw(15) << setprecision(8) << osi->unit_factor_result * xsubprocess[i_p].result_CV[i_q][i_s] << setw(15) << setprecision(8) << osi->unit_factor_result * xsubprocess[i_p].deviation_CV[i_q][i_s] << endl;
	}
	outfile_subprocess.close();
      }
      logger << LOG_DEBUG << "output osi->min_qTcut variation - each subprocess finished" << endl;
    }
    */

    /*
    if (ygeneric->switch_output_subprocess){
      ofstream outfile_result;
      ////////////////////////////////////////////////////////////////////////////////////////////
      //  subprocess qTcut output for CV variation (all qTcut values) in human-readable format  //
      ////////////////////////////////////////////////////////////////////////////////////////////
      if (osi->n_qTcut != 0){
	//      if (x_m == osi->n_moments){filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/result.qTcut." + infix_order_contribution + ".txt";}
	//      cout << "YYY infix_contribution = " << infix_contribution << endl;
	//      cout << "YYY ylist->resultdirectory = " << ylist->resultdirectory << endl;
	if (x_m == osi->n_moments){filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/overview.qTcut." + infix_order_contribution + ".txt";}
	else {filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/result." + result_moment + ".qTcut." + infix_order_contribution + ".txt";}
	outfile_result.open(filename.c_str(), ofstream::out | ofstream::trunc);  
	for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	  //	logger << LOG_DEBUG_VERBOSE << "i_q = " << i_q << endl;
	  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	    //	  logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
	    
	    if (osi->dynamic_scale_CV == 0){outfile_result << "mu / " + muscale + " = " << showpoint << setprecision(6) << osi->mu_fact_CV[1][i_s] / osi->scale_ren << "     ";}
	    else {outfile_result << "mu / mu_DS" << osi->dynamic_scale_CV << " = " << showpoint << setprecision(6) << osi->mu_fact_CV[1][i_s] / osi->scale_ren << " [" << i_s << "]    ";}
	    
	    outfile_result << setw(8) << "qTcut = " << noshowpoint << setw(8) << setprecision(4) << osi->value_qTcut[i_q] << endl;
	    outfile_result << endl;
	    outfile_result << left << setw(20) << "subprocess" << right << setw(26) << signame << " +- " << right << setw(22) << dsigname << " [fb]         ";
	    outfile_result << setw(4) << "chi2" << "   " << setw(5) << "n_run";
	    outfile_result << endl;
	    int temp_size_result = int(log10(abs(osi->unit_factor_result * xsubprocess[0].result_CV[i_q][i_s])));
	    int temp_size_deviation = int(log10(abs(osi->unit_factor_result * xsubprocess[0].deviation_CV[i_q][i_s])));
	    outfile_result << noshowpoint << setw(20) << left << subprocess[0] << "" << right << setw(26) << setprecision(10) << osi->unit_factor_result * xsubprocess[0].result_CV[i_q][i_s] << " +- " << right << setw(22) << setprecision(10) << osi->unit_factor_result * xsubprocess[0].deviation_CV[i_q][i_s] << "     ";
	    outfile_result << endl;
	    outfile_result << endl;
	    
	    for (int i_p = 1; i_p < subprocess.size(); i_p++){
	      //	    logger << LOG_DEBUG_VERBOSE << "i_p = " << i_p << endl;
	      //	  for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	      int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_result * xsubprocess[i_p].result_CV[i_q][i_s])));
	      if (temp_size_result_subprocess < -3){temp_size_result_subprocess -= 5;}
	      int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_result * xsubprocess[i_p].deviation_CV[i_q][i_s])));
	      if (temp_size_deviation_subprocess < -3){temp_size_deviation_subprocess -= 5;}
	      
	      // this segfaults..
	      // cout << setprecision(10+2147483643-7) << 10.0/3 << endl;
	      if (xsubprocess[i_p].result_CV[i_q][i_s] == 0) {
		temp_size_result_subprocess = 0;
	      }
	      if (xsubprocess[i_p].deviation_CV[i_q][i_s] == 0) {
		temp_size_deviation_subprocess = 0;
	      }
	      
	      outfile_result << noshowpoint << setw(20) << left << subprocess[i_p] << "" << right << setw(26) << setprecision(10 + temp_size_result_subprocess - temp_size_result) << showpoint << osi->unit_factor_result * xsubprocess[i_p].result_CV[i_q][i_s] << " +- " << right << setw(22) << setprecision(10 + temp_size_deviation_subprocess - temp_size_deviation) << osi->unit_factor_result * xsubprocess[i_p].deviation_CV[i_q][i_s] << "     ";
	      
	      if (xsubprocess[i_p].no_results > 1){outfile_result << "    " << setw(9) << setprecision(4) << xsubprocess[i_p].chi2_CV[i_q][i_s] << "   " << setw(5) << xsubprocess[i_p].no_results;}
	      else {outfile_result << "    " << setw(9) << setprecision(4) << "---" << "   " << setw(5) << xsubprocess[i_p].no_results;}
	      outfile_result << endl;
	    }
	    outfile_result << endl;
	  }
	}
	outfile_result.close();
	logger << LOG_DEBUG << filename << " finished." << endl;
      }
    }    
    */

    if (ygeneric->switch_output_subprocess > 1){
      ofstream outfile_result;
      /*
      if (osi->n_qTcut != 0){
	//////////////////////////////////////////////////////////////////////////////////////////
	//  subgroup qTcut output for CV variation (all qTcut values) in human-readable format  //
	//////////////////////////////////////////////////////////////////////////////////////////
	if (x_m == osi->n_moments){filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/group/result.qTcut." + infix_order_contribution + ".dat";}
	else {filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/group/result." + result_moment + ".qTcut." + infix_order_contribution + ".dat";}
	outfile_result.open(filename.c_str(), ofstream::out | ofstream::trunc);  
	outfile_result << showpoint;
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  outfile_result << "mu / " + muscale + " = " << setprecision(6) << osi->mu_fact_CV[1][i_s] / osi->scale_ren << endl;
	  outfile_result << endl;
	  outfile_result << "-------------------------------------------------------------" << endl;
	  outfile_result << endl;
	  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	    outfile_result << ygeneric->subgroup[i_g] << endl << endl;
	    outfile_result << setw(15) << "qTcut" << setw(16) << signame << " +- " << setw(12) << dsigname;
	    if (osi->unit_factor_result == 1.){outfile_result << " [fb]";}
	    else if (osi->unit_factor_result == 1.e-3){outfile_result << " [pb]";}
	    //	  if (no_results[i_g] > 1){outfile_result << setw(9) << "chi2";}
	    outfile_result << endl;
	    for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	      //	    outfile_result << setw(15) << setprecision(6) << osi->value_qTcut[i_q] << setw(16) << result_CV[x_m][x_c][i_g][i_q][i_s] << setw(16) << setprecision(4) << deviation_CV[x_m][x_c][i_g][i_q][i_s];
	      outfile_result << noshowpoint << setw(15) << setprecision(4) << osi->value_qTcut[i_q] << setw(16) << osi->unit_factor_result * result_CV[x_m][i_g][i_q][i_s] << setw(16) << setprecision(4) << osi->unit_factor_result * deviation_CV[x_m][i_g][i_q][i_s];
	      //	    if (no_results[i_g] > 1){outfile_result << setw(14) << setprecision(2) << chi2_variation[i_g][i_q][i_s];}
	      outfile_result << endl;
	    }
	    outfile_result << endl;
	    outfile_result << "-------------------------------------------------------------" << endl;
	    outfile_result << endl;
	  }
	}
	outfile_result.close();
	logger << LOG_DEBUG << filename << " finished." << endl;
	
	ofstream outfile_moment;
	if (x_m == osi->n_moments){filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/plot.qTcut." + infix_order_contribution + ".dat";}
	else {filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/plot." + result_moment + ".qTcut." + infix_order_contribution + ".dat";}
	logger << LOG_DEBUG << "Output : xfilename = " << filename << endl;
	outfile_moment.open(filename.c_str(), ofstream::out | ofstream::trunc);  
	outfile_moment << showpoint;
	outfile_moment << "# " << setw(15) << left << "qTcut" << "  ";
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  outfile_moment << "mu_R @ " << right << setprecision(2) << setw(4) << showpoint << i_s << " -- " << "mu_F @ " << setw(4) << showpoint << i_s << "    ";
	}
	int i_g = 0;
	//	for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	  outfile_moment << noshowpoint << setw(10) << setprecision(4) << setw(15) << osi->value_qTcut[i_q];
	  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	    //	  outfile_moment << setw(15) << setprecision(8) << result_CV[x_m][x_c][0][i_q][i_s] << setw(15) << setprecision(8) << deviation_CV[x_m][x_c][0][i_q][i_s];
	    outfile_moment << showpoint << setw(15) << setprecision(8) << osi->unit_factor_result * result_CV[x_m][i_g][i_q][i_s] << setw(15) << setprecision(8) << osi->unit_factor_result * deviation_CV[x_m][i_g][i_q][i_s]; 
	    //xxxxx ???? i_g
	  }
	  outfile_moment << endl;
	}
	outfile_moment.close();
	
	logger << LOG_DEBUG << filename << " finished." << endl;
	
	
	if (x_m == osi->n_moments){filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/plot.CV." + infix_order_contribution + ".dat";}
	else {filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/plot." + result_moment + ".CV." + infix_order_contribution + ".dat";}
	logger << LOG_DEBUG << "Output : xfilename = " << filename << endl;
	outfile_moment.open(filename.c_str(), ofstream::out | ofstream::trunc);  
	outfile_moment << showpoint;
	int i_q = 0;
	//    outfile_moment << setw(10) << osi->value_qTcut[i_q];
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  outfile_moment << setw(15) << setprecision(8) << rel_scale_CV[i_s] << setw(15) << setprecision(8) << osi->unit_factor_result * result_CV[x_m][0][i_q][i_s] << setw(15) << setprecision(8) << osi->unit_factor_result * deviation_CV[x_m][0][i_q][i_s] << endl;
	  //	outfile_moment << setw(15) << setprecision(8) << rel_scale_CV[i_s] << setw(15) << setprecision(8) << scale_CV[i_s] << setw(15) << setprecision(8) << result_CV[x_m][0][i_q][i_s] << setw(15) << setprecision(8) << deviation_CV[x_m][0][i_q][i_s] << endl;
	}
	outfile_moment.close();
	logger << LOG_DEBUG << filename << " finished." << endl;
      */	
	
      string filename;
      ofstream outfile_moment;
      
      if (x_m == osi->n_moments){filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_order_contribution + "/group/plot.qTcut." + infix_order_contribution + ".dat";}
      else {filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_order_contribution + "/group/plot." + result_moment + ".qTcut." + infix_order_contribution + ".dat";}
      logger << LOG_DEBUG << "Output : xfilename = " << filename << endl;
      outfile_moment.open(filename.c_str(), ofstream::out | ofstream::trunc);  
      outfile_moment << showpoint;
      for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	outfile_moment << noshowpoint << setw(10) << setprecision(4) << osi->value_qTcut[i_q];
	for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	    //	    outfile_moment << setw(15) << setprecision(8) << result_CV[x_m][i_g][i_q][i_s] << setw(15) << setprecision(8) << deviation_CV[x_m][i_g][i_q][i_s];
	    outfile_moment << setw(15) << setprecision(8) << osi->unit_factor_result * result_CV[x_m][i_g][i_q][i_s] << setw(15) << setprecision(8) << osi->unit_factor_result * deviation_CV[x_m][i_g][i_q][i_s];
	  }
	}
	outfile_moment << endl;
      }
      outfile_moment.close();
      logger << LOG_DEBUG << filename << " finished." << endl;
    }
  }


  logger << LOG_DEBUG << "output osi->min_qTcut variation finished" << endl;

  //  if ((x_m == 0 || x_m == osi->n_moments)){
  if (ygeneric->switch_output_overview > 2){output_result_overview_qTcut_CV();}
  if (ygeneric->switch_output_result > 2){output_result_qTcut_CV();}
  if (ygeneric->switch_output_plot > 2){output_result_plot_qTcut_CV();}
  if (ygeneric->switch_output_plot > 2){output_result_plot_CV();}

  if (ygeneric->switch_output_plot > 3){output_subprocesses_result_plot_qTcut_CV();}
  if (ygeneric->switch_output_plot > 3){output_subprocesses_result_plot_CV();}
  //output_result_check_CV(result_moment, x_m);
  //  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_contribution::readin_result_contribution_TSV(string result_moment, int x_m){ 
  Logger logger("summary_contribution::readin_result_contribution_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  if (yresultdirectory != ylist->resultdirectory){logger << LOG_INFO << "yresultdirectory = " << yresultdirectory << "  =/=  " << ylist->resultdirectory << " =  ylist->resultdirectory" << endl;}

  //  int n_subgroup = ygeneric->subgroup.size();
  //  string result_file;

  logger << LOG_INFO << "infix_order_contribution = " << infix_order_contribution << endl;

  // only debug output -> temporary !!!
  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "readin_result_contribution_TSV started:" << endl;
  logger << LOG_DEBUG << "type_perturbative_order = " << type_perturbative_order << endl;
  logger << LOG_DEBUG << "type_subtraction_method = " << type_subtraction_method << endl;
  logger << LOG_DEBUG << "in_contribution_order_alpha_s = " << in_contribution_order_alpha_s << endl;
  logger << LOG_DEBUG << "in_contribution_order_alpha_e = " << in_contribution_order_alpha_e << endl;
  logger << LOG_DEBUG << "type_contribution = " << type_contribution << endl;
  logger << LOG_DEBUG << "type_correction = " << type_correction << endl;
  logger << LOG_DEBUG << "interference = " << interference << endl;
  logger << LOG_DEBUG << "photon_induced = " << photon_induced << endl;
  for (int i_d = 0; i_d < directory.size(); i_d++){logger << LOG_DEBUG << "directory[" << setw(3) << "] = " << directory[i_d] << endl;}
  logger << LOG_DEBUG << "subprocess.size() = " << subprocess.size() << endl;
  for (int i_p = 1; i_p < subprocess.size(); i_p++){logger << LOG_DEBUG << "subprocess[" << setw(3) << "] = " << subprocess[i_p] << endl;}
  logger << LOG_DEBUG << "ygeneric->subgroup.size() = " << ygeneric->subgroup.size() << endl;
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG << "ygeneric->subgroup[" << setw(3) << "] = " << ygeneric->subgroup[i_g] << endl;}
  logger << LOG_DEBUG << "subgroup_no_member = " << subgroup_no_member.size() << endl;
  logger << LOG_DEBUG << "ygeneric->subgroup.size() = " << ygeneric->subgroup.size() << endl;
  logger << LOG_DEBUG << "result_moment = " << result_moment << endl;

  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    for (int i_v = 0; i_v < ygeneric->name_scale_variation_TSV.size(); i_v++){
      // only if filled ???
      system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[i_v] + "/" + ylist->resultdirectory);
      // only if filled ???
      system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[i_v] + "/" + ylist->resultdirectory + "/" + infix_contribution);
    }
  }

  //  string readindata;
  //  int n_subprocess = subprocess.size();

  result_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<vector<double> > > > > (ygeneric->subgroup.size(), vector<vector<vector<vector<double> > > > (osi->n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV))));
  deviation_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<vector<double> > > > > (ygeneric->subgroup.size(), vector<vector<vector<vector<double> > > > (osi->n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV))));

  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
    logger << LOG_DEBUG_VERBOSE << "result_TSV[" << i_m << "].size() = " << result_TSV[i_m].size() << endl;
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      logger << LOG_DEBUG_VERBOSE << "result_TSV[" << i_m << "][" << i_g << "].size() = " << result_TSV[i_m][i_g].size() << endl;
      /////      for (int i_q = 0; i_q < output_n_qTcut; i_q++){
      for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	logger << LOG_DEBUG_VERBOSE << "result_TSV[" << i_m << "][" << i_g << "][" << i_q << "].size() = " << result_TSV[i_m][i_g][i_q].size() << endl;
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  logger << LOG_DEBUG_VERBOSE << "result_TSV[" << i_m << "][" << i_g << "][" << i_q << "][" << i_s << "].size() = " << result_TSV[i_m][i_g][i_q][i_s].size() << endl;
	  logger << LOG_DEBUG_VERBOSE << "osi->n_scale_ren_TSV[" << i_s << "] = " << osi->n_scale_ren_TSV[i_s] << endl;
	  result_TSV[i_m][i_g][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  deviation_TSV[i_m][i_g][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  logger << LOG_DEBUG_VERBOSE << "result_TSV[" << i_m << "][" << i_g << "][" << i_q << "][" << i_s << "].size() = " << result_TSV[i_m][i_g][i_q][i_s].size() << endl;
	  logger << LOG_DEBUG_VERBOSE << "result_TSV[" << i_m << "][" << i_g << "][" << i_q << "][" << i_s << "][0].size() = " << result_TSV[i_m][i_g][i_q][i_s][0].size() << endl;
	}
      }
    }
  }









  /*
  vector<string> latex_subprocess = subprocess;
  get_latex_subprocess(latex_subprocess, subprocess);
  */
  //  int n_seed = directory.size();

  logger << LOG_DEBUG_VERBOSE << "xsubprocess.size() = " << xsubprocess.size() << endl;


  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    xsubprocess[i_p].readin_result_TSV(result_moment);
    xsubprocess[i_p].combination_result_TSV(result_moment);
    vector<vector<vector<vector<vector<vector<double> > > > > > ().swap(xsubprocess[i_p].result_run_TSV);
    vector<vector<vector<vector<vector<vector<double> > > > > > ().swap(xsubprocess[i_p].deviation_run_TSV);
  }
  /*
  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    vector<vector<vector<vector<vector<vector<double> > > > > > ().swap(xsubprocess[i_p].result_run_TSV);
    vector<vector<vector<vector<vector<vector<double> > > > > > ().swap(xsubprocess[i_p].deviation_run_TSV);
  }
  */
  /////////////////////////////////////
  //  alternative calculation begin  //
  /////////////////////////////////////
  /*
  for (int i_p = 1; i_p < subprocess.size(); i_p++){
    xsubprocess[i_p].combine_result_alternative_TSV();
  }
  */
  /*
  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    logger << LOG_DEBUG_VERBOSE << left << subprocess[i_p] << ":" << endl;
    if (average_factor == 0){xsubprocess[i_p].combine_result_alternative_TSV();}
    // original version had been switched off (should be identical to hybrid version with average_factor == 1) !!!
    else if (average_factor == 1){xsubprocess[i_p].combine_result_original_TSV();}
    else {xsubprocess[i_p].combine_result_hybrid_TSV();}
  }
  */

  /*
  //////////////////////////////////////////////////
  //  calculate chi2 for individual subprocesses  //
  //////////////////////////////////////////////////
  logger << LOG_DEBUG_VERBOSE << "calculate chi2 for individual subprocesses" << endl;
  /////  for (int i_q = 0; i_q < output_n_qTcut; i_q++){
  for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
	    xsubprocess[i_p].chi2_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	    if (xsubprocess[i_p].n_parallel_runs_TSV[x_m][i_q][i_s][i_r][i_f] > 0){
	      for (int i_z = 0; i_z < directory.size(); i_z++){
		if (!(xsubprocess[i_p].result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0. && xsubprocess[i_p].deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] == 0.)){
		  xsubprocess[i_p].chi2_TSV[x_m][i_q][i_s][i_r][i_f] += pow((xsubprocess[i_p].result_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f] - xsubprocess[i_p].result_TSV[x_m][i_q][i_s][i_r][i_f]) / xsubprocess[i_p].deviation_run_TSV[i_z][x_m][i_q][i_s][i_r][i_f], 2.);
		}
	      }
	      xsubprocess[i_p].chi2_TSV[x_m][i_q][i_s][i_r][i_f] = xsubprocess[i_p].chi2_TSV[x_m][i_q][i_s][i_r][i_f] / xsubprocess[i_p].n_parallel_runs_TSV[x_m][i_q][i_s][i_r][i_f];
	    }
	  }
	}
      }
    }
  }
  */

  /*
  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    vector<vector<vector<vector<vector<vector<double> > > > > > ().swap(xsubprocess[i_p].result_run_TSV);
    vector<vector<vector<vector<vector<vector<double> > > > > > ().swap(xsubprocess[i_p].deviation_run_TSV);
  }
  */

  ////////////////////////////////////////////////////////////
  //  calculate result summed over individual subprocesses  //
  ////////////////////////////////////////////////////////////
  logger << LOG_DEBUG_VERBOSE << "calculate result summed over individual subprocesses" << endl;

  for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  xsubprocess[0].result_TSV[x_m][i_q][i_s][i_r][i_f] = 0.;
	  double dev = 0.;
	  for (int i_p = 1; i_p < subprocess.size(); i_p++){
	    xsubprocess[0].result_TSV[x_m][i_q][i_s][i_r][i_f] += xsubprocess[i_p].result_TSV[x_m][i_q][i_s][i_r][i_f];
	    dev += pow(xsubprocess[i_p].deviation_TSV[x_m][i_q][i_s][i_r][i_f] , 2);
	  }
	  xsubprocess[0].deviation_TSV[x_m][i_q][i_s][i_r][i_f] = sqrt(dev);
	}
      }
    }
  }


 
  /////////////////////////////////////////////////////////////////////////////////////////
  //  calculation of cross sections for subprocess groups (i_g = 0 -> all subprocesses)  //
  /////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////
  //  calculation of subgroup results in TSV variation  //
  ////////////////////////////////////////////////////////
  if (osi->switch_TSV != 0){
    for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	      double dev = 0.;
	      logger << LOG_DEBUG_VERBOSE << "subgroup_no_member[" << i_g << "].size() = " << subgroup_no_member[i_g].size() << endl;
	      for (int i_m = 0; i_m < subgroup_no_member[i_g].size(); i_m++){
		logger << LOG_DEBUG_VERBOSE << "   i_q = " << i_q << "   i_s = " << i_s << "   i_r = " << i_r << "   i_f = " << i_f << "   i_g = " << i_g << "   i_m = " << i_m << "   x_m = " << x_m << endl;
		logger << LOG_DEBUG_VERBOSE << "subgroup_no_member[i_g][i_m] = " << subgroup_no_member[i_g][i_m] << endl;
		logger << LOG_DEBUG_VERBOSE << "result_TSV.size() = " << result_TSV.size() << endl;
		logger << LOG_DEBUG_VERBOSE << "result_TSV[x_m].size() = " << result_TSV[x_m].size() << endl;
		logger << LOG_DEBUG_VERBOSE << "result_TSV[x_m][i_g].size() = " << result_TSV[x_m][i_g].size() << endl;
		logger << LOG_DEBUG_VERBOSE << "result_TSV[x_m][i_g][i_q].size() = " << result_TSV[x_m][i_g][i_q].size() << endl;
		logger << LOG_DEBUG_VERBOSE << "result_TSV[x_m][i_g][i_q][i_s].size() = " << result_TSV[x_m][i_g][i_q].size() << endl;
		logger << LOG_DEBUG_VERBOSE << "result_TSV[x_m][i_g][i_q][i_s][i_r].size() = " << result_TSV[x_m][i_g][i_q][i_s].size() << endl;

		result_TSV[x_m][i_g][i_q][i_s][i_r][i_f] += xsubprocess[subgroup_no_member[i_g][i_m]].result_TSV[x_m][i_q][i_s][i_r][i_f];
		dev += pow(xsubprocess[subgroup_no_member[i_g][i_m]].deviation_TSV[x_m][i_q][i_s][i_r][i_f], 2);
	      }

	      deviation_TSV[x_m][i_g][i_q][i_s][i_r][i_f] = sqrt(dev);
	      stringstream temp_res;
	      temp_res << setw(23) << setprecision(15) << result_TSV[x_m][i_g][i_q][i_s][i_r][i_f];
	      stringstream temp_dev;
	      temp_dev << setw(23) << setprecision(15) << deviation_TSV[x_m][i_g][i_q][i_s][i_r][i_f];
	      logger << LOG_DEBUG_VERBOSE << "XXX   result_TSV[" << x_m << "][" << i_g << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
	    }
	  }
	}
      }
    }
  }
  
  if (ygeneric->switch_output_overview > 2){output_result_overview_qTcut_TSV();}
  if (ygeneric->switch_output_result > 2){output_result_qTcut_TSV();}
  if (ygeneric->switch_output_plot > 2){output_result_plot_qTcut_TSV();}

  if (ygeneric->switch_output_plot > 3){output_subprocesses_result_plot_qTcut_TSV();}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



