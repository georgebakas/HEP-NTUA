#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_contribution::output_result_overview_qTcut_TSV(){ 
  Logger logger("summary_contribution::output_result_overview_qTcut_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  if (yresultdirectory != ylist->resultdirectory){logger << LOG_INFO << "yresultdirectory = " << yresultdirectory << "  =/=  " << ylist->resultdirectory << " =  ylist->resultdirectory" << endl;}

  string filename_complete;
  string filename_ren;
  string filename_fact;
  string filename_equal;
  string filename_antipodal;
  
  ofstream outfile_complete;
  ofstream outfile_ren;
  ofstream outfile_fact;
  ofstream outfile_equal;
  ofstream outfile_antipodal;

  stringstream temp_ss_separation_one;
  stringstream temp_ss_separation_two;
  int counter_separation = 15 + 15 + 15 + 25 + 15 + 4 + 15 + 3 + 10 + 3 + 5 + 3 + 9;
  for (int i_x = 0; i_x < counter_separation; i_x++){
    temp_ss_separation_one << "-";
    temp_ss_separation_two << "=";
  }
  
  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[0] + "/" + ylist->resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
      filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[1] + "/" + ylist->resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
      filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[2] + "/" + ylist->resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
      filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[3] + "/" + ylist->resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
      filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[4] + "/" + ylist->resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
      
      logger << LOG_INFO << "FILENAME (complete):   " << filename_complete << endl;

      outfile_complete.open(filename_complete.c_str(), ofstream::out | ofstream::trunc);  
      outfile_ren.open(filename_ren.c_str(), ofstream::out | ofstream::trunc);  
      outfile_fact.open(filename_fact.c_str(), ofstream::out | ofstream::trunc);  
      outfile_equal.open(filename_equal.c_str(), ofstream::out | ofstream::trunc);  
      outfile_antipodal.open(filename_antipodal.c_str(), ofstream::out | ofstream::trunc);  

      stringstream name_ss;
      name_ss << "dXS +- err (subprocess)";

      stringstream header_ss;
      header_ss << left
		<< setw(15) << "qTcut"
		<< setw(5) << "mu_R/" << setw(10) << osi->name_extended_set_TSV[i_s] 
		<< setw(5) << "mu_F/" << setw(10) << osi->name_extended_set_TSV[i_s] 
		<< setw(25) << "subprocess"
		<< right << setw(34) << name_ss.str()
		<< "   "
		<< right << setw(10) << "chi2_dof" 
		<< "   " 
		<< right << setw(5) << "n_run"
		<< "   " 
		<< setw(9) << "(removed)";

      for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	outfile_complete << header_ss.str() << endl;
	outfile_ren << header_ss.str() << endl;
	outfile_fact << header_ss.str() << endl;
	outfile_equal << header_ss.str() << endl;
	outfile_antipodal << header_ss.str() << endl;

	outfile_complete << endl;
	outfile_ren << endl;
	outfile_fact << endl;
	outfile_equal << endl;
	outfile_antipodal << endl;

	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    int temp_size_result = int(log10(abs(osi->unit_factor_result * xsubprocess[0].result_TSV[i_m][i_q][i_s][i_r][i_f])));
	    for (int i_p = 0; i_p < subprocess.size(); i_p++){
	      int temp_size_result_c = int(log10(abs(osi->unit_factor_result * xsubprocess[i_p].result_TSV[i_m][i_q][i_s][i_r][i_f])));
	      if (temp_size_result_c > temp_size_result){temp_size_result = temp_size_result_c;}
	    }
	    int temp_size_deviation = int(log10(abs(osi->unit_factor_result * xsubprocess[0].deviation_TSV[i_m][i_q][i_s][i_r][i_f])));
	    if (osi->unit_factor_result * xsubprocess[0].deviation_TSV[i_m][i_q][i_s][i_r][i_f] >= 1.){temp_size_deviation++;}

	    logger << LOG_DEBUG << "CENTRE " << setw(15) << osi->unit_factor_result * xsubprocess[0].result_TSV[i_m][i_q][i_s][i_r][i_f] 
		   << " --- " << setw(2) << temp_size_result << "     "
		   << setw(15) << osi->unit_factor_result * xsubprocess[0].deviation_TSV[i_m][i_q][i_s][i_r][i_f] 
		   << " --- " << setw(2) << temp_size_deviation << endl;


	    for (int i_p = 0; i_p < subprocess.size(); i_p++){
     	      int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_result * xsubprocess[i_p].result_TSV[i_m][i_q][i_s][i_r][i_f])));
	      int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_result * xsubprocess[i_p].deviation_TSV[i_m][i_q][i_s][i_r][i_f])));
	      if (osi->unit_factor_result * xsubprocess[i_p].deviation_TSV[i_m][i_q][i_s][i_r][i_f] >= 1.){temp_size_deviation_subprocess++;}
	      if (xsubprocess[i_p].result_TSV[i_m][i_q][i_s][i_r][i_f] == 0){temp_size_result_subprocess = 0;}
	      if (xsubprocess[i_p].deviation_TSV[i_m][i_q][i_s][i_r][i_f] == 0){temp_size_deviation_subprocess = 0;}
	      int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
	      int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
	      if (xsubprocess[i_p].result_TSV[i_m][i_q][i_s][i_r][i_f] == 0.){setw_result = 1; setw_deviation = 1;}
	      else if (abs(osi->unit_factor_result * xsubprocess[i_p].result_TSV[i_m][i_q][i_s][i_r][i_f]) < 1.){setw_result--;}

	      logger << LOG_DEBUG_VERBOSE << "TEMP   " 
		     << setw(15) << osi->unit_factor_result * xsubprocess[i_p].result_TSV[i_m][i_q][i_s][i_r][i_f] 
		     << " --- " << setw(2) << setw_result << "     "
		     << setw(15) << osi->unit_factor_result * xsubprocess[i_p].deviation_TSV[i_m][i_q][i_s][i_r][i_f] 
		     << " --- " << setw(2) << setw_deviation << endl;

	      stringstream temp_ss;
	      temp_ss << left;
	      if (active_qTcut){temp_ss << setw(15) << setprecision(8) << osi->value_qTcut[i_q];}
	      else {temp_ss << setw(15) << "independent";}
	      temp_ss << setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r]
		      << setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f];
	      if (i_p == 0){temp_ss << setw(25) << infix_order_contribution;}
	      else {temp_ss << setw(25) << subprocess[i_p];}
	      temp_ss << right << setw(15) << setprecision(setw_result) << showpoint 
		      << osi->unit_factor_result * xsubprocess[i_p].result_TSV[i_m][i_q][i_s][i_r][i_f] 
		      << " +- " 
		      << right << setw(15) << setprecision(setw_deviation) << showpoint
		      << osi->unit_factor_result * xsubprocess[i_p].deviation_TSV[i_m][i_q][i_s][i_r][i_f];
	      if (i_p == 0){temp_ss << setw(3 + 10 + 3 + 5 + 9);}
	      else {
		temp_ss	<< "   "
			<< right << showpoint << setprecision(4) << setw(10) << xsubprocess[i_p].chi2_TSV[i_m][i_q][i_s][i_r][i_f] 
			<< "   "
			<< setw(5) << right << xsubprocess[i_p].n_parallel_runs_TSV[i_m][i_q][i_s][i_r][i_f]
			<< "   "
			<< "   (" << setw(4) << right << accumulate(xsubprocess[i_p].remove_run_qTcut_TSV[i_q].begin(), xsubprocess[i_p].remove_run_qTcut_TSV[i_q].end(), 0) << ")";
	      }

	      outfile_complete << temp_ss.str() << endl;
	      if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str() << endl;}
	      if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str() << endl;}
	      if (i_r == i_f){outfile_equal << temp_ss.str() << endl;}
	      if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str() << endl;}

	      if (i_p == 0){
		outfile_complete << temp_ss_separation_one.str() << endl;
		if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss_separation_one.str() << endl;}
		if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss_separation_one.str() << endl;}
		if (i_r == i_f){outfile_equal << temp_ss_separation_one.str() << endl;}
		if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss_separation_one.str() << endl;}
	      }
	    }
	    outfile_complete << endl;
	    if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << endl;}
	    if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << endl;}
	    if (i_r == i_f){outfile_equal << endl;}
	    if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << endl;}
	  }
	}
	outfile_complete << temp_ss_separation_two.str() << endl << endl;
	outfile_ren << temp_ss_separation_two.str() << endl << endl;
	outfile_fact << temp_ss_separation_two.str() << endl << endl;
	outfile_equal << temp_ss_separation_two.str() << endl << endl;
	outfile_antipodal << temp_ss_separation_two.str() << endl << endl;
      }
      outfile_complete.close();
      outfile_ren.close();
      outfile_fact.close();
      outfile_equal.close();
      outfile_antipodal.close();
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_contribution::output_result_overview_qTcut_CV(){
  Logger logger("summary_contribution::output_result_overview_qTcut_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string filename;
  string signame = "XS_" + infix_order_contribution;
  string dsigname = "dXS_" + infix_order_contribution;

  stringstream muscale_ss;
  muscale_ss << setprecision(6) << osi->scale_ren;
  string muscale = muscale_ss.str();
 
  stringstream temp_ss_separation_one;
  stringstream temp_ss_separation_two;
  int counter_separation = 15 + 15 + 15 + 25 + 15 + 4 + 15 + 3 + 10 + 3 + 5 + 3 + 9;
  for (int i_x = 0; i_x < counter_separation; i_x++){
    temp_ss_separation_one << "-";
    temp_ss_separation_two << "=";
  }

  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    ofstream outfile_result;
    filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
    
    logger << LOG_INFO << "FILENAME:   " << filename << endl;
    
    outfile_result.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    
    stringstream name_ss;
    name_ss << "dXS +- err (subprocess)";
    
    stringstream header_ss;
    header_ss << left
	      << setw(15) << "qTcut"
	      << setw(5) << "mu_R/" << setw(10) << "CV" 
	      << setw(5) << "mu_F/" << setw(10) << "CV"
	      << setw(25) << "subprocess"
	      << right << setw(34) << name_ss.str()
	      << "   "
	      << right << setw(10) << "chi2_dof" 
	      << "   " 
	      << right << setw(5) << "n_run"
	      << "   " 
	      << setw(9) << "(removed)";
    
    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
      outfile_result << header_ss.str() << endl;
      outfile_result << endl;
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	int temp_size_result = int(log10(abs(osi->unit_factor_result * xsubprocess[0].result_CV[i_q][i_s])));
	for (int i_p = 0; i_p < subprocess.size(); i_p++){
	  int temp_size_result_c = int(log10(abs(osi->unit_factor_result * xsubprocess[i_p].result_CV[i_q][i_s])));
	  if (temp_size_result_c > temp_size_result){temp_size_result = temp_size_result_c;}
	}
	int temp_size_deviation = int(log10(abs(osi->unit_factor_result * xsubprocess[0].deviation_CV[i_q][i_s])));
	if (osi->unit_factor_result * xsubprocess[0].deviation_CV[i_q][i_s] >= 1.){temp_size_deviation++;}
	
	for (int i_p = 0; i_p < subprocess.size(); i_p++){
	  int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_result * xsubprocess[i_p].result_CV[i_q][i_s])));
	  int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_result * xsubprocess[i_p].deviation_CV[i_q][i_s])));
	  if (osi->unit_factor_result * xsubprocess[i_p].deviation_CV[i_q][i_s] >= 1.){temp_size_deviation_subprocess++;}
	  if (xsubprocess[i_p].result_CV[i_q][i_s] == 0){temp_size_result_subprocess = 0;}
	  if (xsubprocess[i_p].deviation_CV[i_q][i_s] == 0){temp_size_deviation_subprocess = 0;}
	  int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
	  int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
	  if (xsubprocess[i_p].result_CV[i_q][i_s] == 0.){setw_result = 1; setw_deviation = 1;}
	  else if (abs(osi->unit_factor_result * xsubprocess[i_p].result_CV[i_q][i_s]) < 1.){setw_result--;}
  
	  stringstream temp_ss;
	  temp_ss << left;
	  if (active_qTcut){temp_ss << setw(15) << setprecision(8) << osi->value_qTcut[i_q];}
	  temp_ss << setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
		  << setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s];
	  if (i_p == 0){temp_ss << setw(25) << infix_order_contribution;}
	  else {temp_ss << setw(25) << subprocess[i_p];}
	  temp_ss << showpoint 
		  << right << setw(15) << setprecision(setw_result) << osi->unit_factor_result * xsubprocess[i_p].result_CV[i_q][i_s]
		  << " +- " 
		  << right << setw(15) << setprecision(setw_deviation) << osi->unit_factor_result * xsubprocess[i_p].deviation_CV[i_q][i_s];
	  if (i_p == 0){temp_ss << setw(3 + 10 + 3 + 5 + 9);}
	  else {
	    temp_ss << "   "
		    << right << showpoint << setprecision(4) << setw(10) << xsubprocess[i_p].chi2_CV[i_q][i_s]
		    << "   "
		    << setw(5) << right << xsubprocess[i_p].n_parallel_runs_CV[i_q][i_s]
	      //		    << setw(5) << right << xsubprocess[i_p].no_results
		    << "   "
		    << "   (" << setw(4) << right << accumulate(xsubprocess[i_p].remove_run_qTcut_CV[i_q].begin(), xsubprocess[i_p].remove_run_qTcut_CV[i_q].end(), 0) << ")";
	  }
	  
	  outfile_result << temp_ss.str() << endl;

	  if (i_p == 0){
	    outfile_result << temp_ss_separation_one.str() << endl;
	  }
	}
	outfile_result << endl;
      }
      outfile_result << temp_ss_separation_two.str() << endl << endl;
    }
    outfile_result.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_contribution::output_result_qTcut_TSV(){
  Logger logger("summary_contribution::output_result_qTcut_TSV");
  logger << LOG_DEBUG << "called" << endl;

  string filename_complete;
  string filename_ren;
  string filename_fact;
  string filename_equal;
  string filename_antipodal;
  
  ofstream outfile_complete;
  ofstream outfile_ren;
  ofstream outfile_fact;
  ofstream outfile_equal;
  ofstream outfile_antipodal;

  logger << LOG_INFO << "resultdirectory = " << resultdirectory << endl;
  logger << LOG_INFO << "processname = " << processname << endl;
  logger << LOG_INFO << "type_perturbative_order = " << type_perturbative_order << endl;
  logger << LOG_INFO << "type_subtraction_method = " << type_subtraction_method << endl;
  logger << LOG_INFO << "in_contribution_order_alpha_s = " << in_contribution_order_alpha_s << endl;
  logger << LOG_INFO << "in_contribution_order_alpha_e = " << in_contribution_order_alpha_e << endl;

  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/complete/" + ylist->resultdirectory + "result" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
      filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/ren/" + ylist->resultdirectory + "/result" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
      filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/fact/" + ylist->resultdirectory + "/result" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
      filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/equal/" + ylist->resultdirectory + "/result" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
      filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/antipodal/" + ylist->resultdirectory + "/result" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
            
      logger << LOG_INFO << "FILENAME (complete):   " << filename_complete << endl;
      
      outfile_complete.open(filename_complete.c_str(), ofstream::out | ofstream::trunc);  
      outfile_ren.open(filename_ren.c_str(), ofstream::out | ofstream::trunc);  
      outfile_fact.open(filename_fact.c_str(), ofstream::out | ofstream::trunc);  
      outfile_equal.open(filename_equal.c_str(), ofstream::out | ofstream::trunc);  
      outfile_antipodal.open(filename_antipodal.c_str(), ofstream::out | ofstream::trunc);  

      stringstream name_ss;
      name_ss << "dXS +- err (" << infix_order_contribution << ")";
      
      stringstream header_ss;
      header_ss << left
		<< setw(15) << "qTcut"
		<< setw(5) << "mu_R/" << setw(10) << osi->name_extended_set_TSV[i_s] 
		<< setw(5) << "mu_F/" << setw(10) << osi->name_extended_set_TSV[i_s] 
		<< setw(34) << name_ss.str();
      if (active_qTcut){
	header_ss << setw(15) << "err(stat)" 
		  << setw(15) << "err(extr)";
      }
      
      outfile_complete << header_ss.str() << endl;
      outfile_ren << header_ss.str() << endl;
      outfile_fact << header_ss.str() << endl;
      outfile_equal << header_ss.str() << endl;
      outfile_antipodal << header_ss.str() << endl;

      for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	outfile_complete << endl;
	outfile_ren << endl;
	outfile_fact << endl;
	outfile_equal << endl;
	outfile_antipodal << endl;

	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    stringstream temp_ss;
	    temp_ss << left;
	    if (active_qTcut){temp_ss << setw(15) << osi->value_qTcut[i_q];}
	    else {temp_ss << setw(15) << "independent";}
	    temp_ss << setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r] 
		    << setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f] 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_q][i_s][i_r][i_f] 
		    << " +- "
		    << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_q][i_s][i_r][i_f];
	    outfile_complete << temp_ss.str() << endl;
	    if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str() << endl;}
	    if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str() << endl;}
	    if (i_r == i_f){outfile_equal << temp_ss.str() << endl;}
	    if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str() << endl;}
	  }
	}
      }
      outfile_complete.close();
      outfile_ren.close();
      outfile_fact.close();
      outfile_equal.close();
      outfile_antipodal.close();
    }
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}


void summary_contribution::output_result_qTcut_CV(){ 
  Logger logger("summary_contribution::output_result_qTcut_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string filename;
  ofstream outfile;

  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/result" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";

    logger << LOG_INFO << "FILENAME:   " << filename << endl;
    
    outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    
    stringstream name_ss;
    name_ss << "dXS +- err (" << infix_order_contribution << ")";
    
    stringstream header_ss;
    header_ss << left
	      << setw(15) << "qTcut"
	      << setw(5) << "mu_R/" << setw(10) << "CV" 
	      << setw(5) << "mu_F/" << setw(10) << "CV" 
	      << setw(34) << name_ss.str();
    
    outfile << header_ss.str() << endl;

    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
      outfile << endl;
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	stringstream temp_ss;
	temp_ss << noshowpoint << left;
	if (active_qTcut){temp_ss << setw(15) << osi->value_qTcut[i_q];}
	else {outfile << setw(15) << "independent";}
	temp_ss << setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
		<< setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s]
		<< showpoint 
		<< setprecision(8) << setw(15) << osi->unit_factor_result * result_CV[i_m][0][i_q][i_s]
		<< " +- "
		<< setprecision(8) << setw(15) << osi->unit_factor_result * deviation_CV[i_m][0][i_q][i_s];
	outfile << temp_ss.str() << endl;
      }  
    }
  }
  outfile.close();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}





void summary_contribution::output_result_plot_qTcut_TSV(){ 
  Logger logger("summary_contribution::output_result_plot_qTcut_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string filename_complete;
  string filename_ren;
  string filename_fact;
  string filename_equal;
  string filename_antipodal;
  
  ofstream outfile_complete;
  ofstream outfile_ren;
  ofstream outfile_fact;
  ofstream outfile_equal;
  ofstream outfile_antipodal;

  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[0] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + infix_order_contribution + ".dat";
      filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[1] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + infix_order_contribution + ".dat";
      filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[2] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + infix_order_contribution + ".dat";
      filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[3] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + infix_order_contribution + ".dat";
      filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[4] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + infix_order_contribution + ".dat";

      outfile_complete.open(filename_complete.c_str(), ofstream::out | ofstream::trunc);  
      outfile_ren.open(filename_ren.c_str(), ofstream::out | ofstream::trunc);  
      outfile_fact.open(filename_fact.c_str(), ofstream::out | ofstream::trunc);  
      outfile_equal.open(filename_equal.c_str(), ofstream::out | ofstream::trunc);  
      outfile_antipodal.open(filename_antipodal.c_str(), ofstream::out | ofstream::trunc);  
      
      outfile_complete << left << setw(15) << "# qTcut";
      outfile_ren << left << setw(15) << "# qTcut";
      outfile_fact << left << setw(15) << "# qTcut";
      outfile_equal << left << setw(15) << "# qTcut";
      outfile_antipodal << left << setw(15) << "# qTcut";
      
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  stringstream temp_ss;
	  temp_ss << "mu_R = " 
		  << setw(4) << setprecision(2) << showpoint << osi->relative_scale_ren_TSV[i_s][i_r] 
		  << " -- " 
		  << "mu_F = " 
		  << setw(4) << setprecision(2) << showpoint << osi->relative_scale_fact_TSV[i_s][i_f] << "    ";
	  
	  outfile_complete << setw(30) << temp_ss.str() << left;
	  if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << setw(30) << temp_ss.str() << left;}
	  if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << setw(30) << temp_ss.str() << left;}
	  if (i_r == i_f){outfile_equal << setw(30) << temp_ss.str() << left;}
	  if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << setw(30) << temp_ss.str() << left;}
	}
      }
      outfile_complete << endl;
      outfile_ren << endl;
      outfile_fact << endl;
      outfile_equal << endl;
      outfile_antipodal << endl;
      
      if (!active_qTcut){
	outfile_complete << setw(15) << "0";
	outfile_equal << setw(15) << "0";
	outfile_antipodal << setw(15) << "0";
	outfile_ren << setw(15) << "0";
	outfile_fact << setw(15) << "0";
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    stringstream temp_ss;
	    temp_ss << left 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][0][i_s][i_r][i_f] 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][0][i_s][i_r][i_f];
	    outfile_complete << temp_ss.str();
	    if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str();}
	    if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact <<temp_ss.str();}
	    if (i_r == i_f){outfile_equal << temp_ss.str();}
	    if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str();}
	  }
	}
	outfile_complete << endl;
	outfile_ren << endl;
	outfile_fact << endl;
	outfile_equal << endl;
	outfile_antipodal << endl;
      }

      for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	outfile_complete << setw(15) << osi->value_qTcut[i_q];
	outfile_equal << setw(15) << osi->value_qTcut[i_q];
	outfile_antipodal << setw(15) << osi->value_qTcut[i_q];
	outfile_ren << setw(15) << osi->value_qTcut[i_q];
	outfile_fact << setw(15) << osi->value_qTcut[i_q];
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    stringstream temp_ss;
	    temp_ss << left << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_q][i_s][i_r][i_f] 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_q][i_s][i_r][i_f];
	    outfile_complete << temp_ss.str();
	    if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str();}
	    if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str();}
	    if (i_r == i_f){outfile_equal << temp_ss.str();}
	    if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str();}
	  }
	}
	outfile_complete << endl;
	outfile_ren << endl;
	outfile_fact << endl;
	outfile_equal << endl;
	outfile_antipodal << endl;
      }
      outfile_complete.close();
      outfile_ren.close();
      outfile_fact.close();
      outfile_equal.close();
      outfile_antipodal.close();
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_contribution::output_result_plot_qTcut_CV(){
  Logger logger("summary_contribution::output_result_plot_qTcut_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string filename;

  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    ofstream outfile;
    ///////////////////////////////////////////////////////////////////
    //  subprocess qTcut output for CV variation (all qTcut values)  //
    ///////////////////////////////////////////////////////////////////
    filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + infix_order_contribution + ".dat";
    
    logger << LOG_DEBUG << "FILENAME:   " << filename << endl;
    
    outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    outfile << showpoint;
    outfile << left << setw(15) << "# qTcut";
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      stringstream temp_ss;
      temp_ss << "mu_R = " << setprecision(2) << setw(4) << i_s
	      << " -- " 
	      << "mu_F = " << setw(4) << i_s
	      << "    ";
      outfile << setw(30) << temp_ss.str() << left;
    }
    outfile << endl;
    for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
      outfile << noshowpoint << left << setw(15) << setprecision(4) << osi->value_qTcut[i_q];
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	outfile << showpoint 
		<< setw(15) << setprecision(8) << osi->unit_factor_result * result_CV[i_m][0][i_q][i_s] 
		<< setw(15) << setprecision(8) << osi->unit_factor_result * deviation_CV[i_m][0][i_q][i_s];
      }
      outfile << endl;
    }
    outfile.close();
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



/* doubled !!!
void summary_contribution::output_result_plot_CV(){ 
  Logger logger("summary_contribution::output_result_plot_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string filename;
  
  vector<double> rel_scale_CV(osi->n_scales_CV);
  vector<double> scale_CV(osi->n_scales_CV);
  logger << LOG_DEBUG << "osi->n_scales_CV = " << osi->n_scales_CV << endl;
  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
    if (osi->n_scales_CV > 1){rel_scale_CV[i_s] = pow(10,log10(double(osi->variation_factor_CV)) * double(2. * i_s / (osi->n_scales_CV - 1) - 1));}
    else {rel_scale_CV[i_s] = 1;}
    scale_CV[i_s] = rel_scale_CV[i_s] * osi->scale_ren;
  }
  
  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    ofstream outfile_subprocess;
    ////////////////////////////////////////////////////////////////
    //  subprocess output for CV variation (osi->min_qTcut value)  //
    ////////////////////////////////////////////////////////////////
    logger.newLine(LOG_DEBUG);
    logger << LOG_DEBUG << "subprocess output for CV variation (osi->min_qTcut value)" << endl;
    filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + "/plot" + ygeneric->infix_name_moment[i_m] + ".CV." + infix_order_contribution + ".dat";

    logger << LOG_INFO << "FILENAME:   " << filename << endl;

    outfile_subprocess.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    outfile_subprocess << showpoint;
    int i_q = 0;
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      outfile_subprocess << setw(15) << setprecision(8) << rel_scale_CV[i_s] 
			 << setw(15) << setw(15) << setprecision(8) << osi->unit_factor_result * result_CV[i_m][0][i_q][i_s] 
			 << setw(15) << setprecision(8) << osi->unit_factor_result * deviation_CV[i_m][0][i_q][i_s] << endl;
    }
    outfile_subprocess.close();
  }
 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/





// contains actually subprocess output !!!

void summary_contribution::output_subprocesses_result_plot_qTcut_TSV(){ 
  Logger logger("summary_contribution::output_subprocesses_result_plot_qTcut_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string filename_complete;
  string filename_ren;
  string filename_fact;
  string filename_equal;
  string filename_antipodal;
  
  ofstream outfile_complete;
  ofstream outfile_ren;
  ofstream outfile_fact;
  ofstream outfile_equal;
  ofstream outfile_antipodal;

  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    for (int i_p = 1; i_p < subprocess.size(); i_p++){
      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[0] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + subprocess[i_p] + ".dat";
	filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[1] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + subprocess[i_p] + ".dat";
	filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[2] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + subprocess[i_p] + ".dat";
	filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[3] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + subprocess[i_p] + ".dat";
	filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[4] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + subprocess[i_p] + ".dat";

	outfile_complete.open(filename_complete.c_str(), ofstream::out | ofstream::trunc);  
	outfile_ren.open(filename_ren.c_str(), ofstream::out | ofstream::trunc);  
	outfile_fact.open(filename_fact.c_str(), ofstream::out | ofstream::trunc);  
	outfile_equal.open(filename_equal.c_str(), ofstream::out | ofstream::trunc);  
	outfile_antipodal.open(filename_antipodal.c_str(), ofstream::out | ofstream::trunc);  
	
	outfile_complete << left << setw(15) << "# qTcut";
	outfile_ren << left << setw(15) << "# qTcut";
	outfile_fact << left << setw(15) << "# qTcut";
	outfile_equal << left << setw(15) << "# qTcut";
	outfile_antipodal << left << setw(15) << "# qTcut";

	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    stringstream temp_ss;
	    temp_ss << "mu_R = " 
		    << setw(4) << setprecision(2) << showpoint << osi->relative_scale_ren_TSV[i_s][i_r] 
		    << " -- " 
		    << "mu_F = " 
		    << setw(4) << setprecision(2) << showpoint << osi->relative_scale_fact_TSV[i_s][i_f] << "    ";

	    outfile_complete << setw(30) << temp_ss.str() << left;
	    if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << setw(30) << temp_ss.str() << left;}
	    if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << setw(30) << temp_ss.str() << left;}
	    if (i_r == i_f){outfile_equal << setw(30) << temp_ss.str() << left;}
	    if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << setw(30) << temp_ss.str() << left;}
	  }
	}
	outfile_complete << endl;
	outfile_ren << endl;
	outfile_fact << endl;
	outfile_equal << endl;
	outfile_antipodal << endl;

	if (!active_qTcut){
  	  outfile_complete << setw(15) << "0";
  	  outfile_equal << setw(15) << "0";
  	  outfile_antipodal << setw(15) << "0";
  	  outfile_ren << setw(15) << "0";
  	  outfile_fact << setw(15) << "0";
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      stringstream temp_ss;
	      temp_ss << left 
		      << setprecision(8) << setw(15) << osi->unit_factor_result * xsubprocess[i_p].result_TSV[i_m][0][i_s][i_r][i_f] 
		      << setprecision(8) << setw(15) << osi->unit_factor_result * xsubprocess[i_p].deviation_TSV[i_m][0][i_s][i_r][i_f];
	      outfile_complete << temp_ss.str();
	      if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str();}
	      if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact <<temp_ss.str();}
	      if (i_r == i_f){outfile_equal << temp_ss.str();}
	      if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str();}
	    }
	  }
	  outfile_complete << endl;
	  outfile_ren << endl;
	  outfile_fact << endl;
	  outfile_equal << endl;
	  outfile_antipodal << endl;
	}

	for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
  	  outfile_complete << setw(15) << osi->value_qTcut[i_q];
  	  outfile_equal << setw(15) << osi->value_qTcut[i_q];
  	  outfile_antipodal << setw(15) << osi->value_qTcut[i_q];
  	  outfile_ren << setw(15) << osi->value_qTcut[i_q];
  	  outfile_fact << setw(15) << osi->value_qTcut[i_q];
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      stringstream temp_ss;
	      temp_ss << left << setprecision(8) << setw(15) << osi->unit_factor_result * xsubprocess[i_p].result_TSV[i_m][i_q][i_s][i_r][i_f] 
		      << setprecision(8) << setw(15) << osi->unit_factor_result * xsubprocess[i_p].deviation_TSV[i_m][i_q][i_s][i_r][i_f];
	      outfile_complete << temp_ss.str();
	      if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str();}
	      if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str();}
	      if (i_r == i_f){outfile_equal << temp_ss.str();}
	      if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str();}
	    }
	  }
	  outfile_complete << endl;
	  outfile_ren << endl;
	  outfile_fact << endl;
	  outfile_equal << endl;
	  outfile_antipodal << endl;
	}
	outfile_complete.close();
	outfile_ren.close();
	outfile_fact.close();
	outfile_equal.close();
	outfile_antipodal.close();
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


// contains actually subprocess output !!!

void summary_contribution::output_subprocesses_result_plot_qTcut_CV(){
  Logger logger("summary_contribution::output_subprocesses_result_plot_qTcut_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string filename;

  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    ofstream outfile_subprocess;
    ///////////////////////////////////////////////////////////////////
    //  subprocess qTcut output for CV variation (all qTcut values)  //
    ///////////////////////////////////////////////////////////////////
    for (int i_p = 1; i_p < subprocess.size(); i_p++){
      filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + subprocess[i_p] + ".dat";

      logger << LOG_DEBUG << "FILENAME:   " << filename << endl;

      outfile_subprocess.open(filename.c_str(), ofstream::out | ofstream::trunc);  
      outfile_subprocess << showpoint;
      outfile_subprocess << left << setw(15) << "# qTcut";
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	stringstream temp_ss;
	temp_ss << "mu_R = " << setprecision(2) << setw(4) << i_s
	  //	showpoint << osi->relative_scale_ren_TSV[i_s][i_r]
		<< " -- " 
		<< "mu_F = " << setw(4) << i_s
	  //    showpoint << osi->relative_scale_fact_TSV[i_s][i_f] 
		<< "    ";
	outfile_subprocess << setw(30) << temp_ss.str() << left;
	//	outfile_subprocess << "mu_R @ " << right << setprecision(2) << setw(4) << showpoint << i_s << " -- " << "mu_F @ " << setw(4) << showpoint << i_s << "    ";
      }
      outfile_subprocess << endl;
      //      outfile_subprocess << showpoint;
      for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	outfile_subprocess << noshowpoint << left << setw(15) << setprecision(4) << osi->value_qTcut[i_q];
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  outfile_subprocess 
	    << showpoint << setw(15) << setprecision(8) << osi->unit_factor_result * xsubprocess[i_p].result_CV[i_q][i_s] 
	    << setw(15) << setprecision(8) << osi->unit_factor_result * xsubprocess[i_p].deviation_CV[i_q][i_s];
	}
	outfile_subprocess << endl;
      }
      outfile_subprocess.close();
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}












void summary_contribution::output_result_plot_CV(){ 
  Logger logger("summary_contribution::output_result_plot_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string filename;
  
  vector<double> rel_scale_CV(osi->n_scales_CV);
  vector<double> scale_CV(osi->n_scales_CV);
  logger << LOG_DEBUG << "osi->n_scales_CV = " << osi->n_scales_CV << endl;
  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
    if (osi->n_scales_CV > 1){rel_scale_CV[i_s] = pow(10,log10(double(osi->variation_factor_CV)) * double(2. * i_s / (osi->n_scales_CV - 1) - 1));}
    else {rel_scale_CV[i_s] = 1;}
    scale_CV[i_s] = rel_scale_CV[i_s] * osi->scale_ren;
  }
  
  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    ofstream outfile_subprocess;
    ////////////////////////////////////////////////////////////////
    //  subprocess output for CV variation (osi->min_qTcut value)  //
    ////////////////////////////////////////////////////////////////
    logger.newLine(LOG_DEBUG);
    logger << LOG_DEBUG << "subprocess output for CV variation (osi->min_qTcut value)" << endl;
    filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + "/plot" + ygeneric->infix_name_moment[i_m] + ".CV." + infix_order_contribution + ".dat";

    logger << LOG_INFO << "FILENAME:   " << filename << endl;

    outfile_subprocess.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    outfile_subprocess << showpoint;
    int i_q = 0;
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      outfile_subprocess << setw(15) << setprecision(8) << rel_scale_CV[i_s] 
			 << setw(15) << setw(15) << setprecision(8) << osi->unit_factor_result * result_CV[i_m][0][i_q][i_s] 
			 << setw(15) << setprecision(8) << osi->unit_factor_result * deviation_CV[i_m][0][i_q][i_s] << endl;
    }
    outfile_subprocess.close();
  }
 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}





// contains actually subprocess output !!!

void summary_contribution::output_subprocesses_result_plot_CV(){ 
  //string result_moment, int x_m
  Logger logger("summary_contribution::output_subprocesses_result_plot_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string filename;

  
  vector<double> rel_scale_CV(osi->n_scales_CV);
  vector<double> scale_CV(osi->n_scales_CV);
  logger << LOG_DEBUG << "osi->n_scales_CV = " << osi->n_scales_CV << endl;
  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
    if (osi->n_scales_CV > 1){rel_scale_CV[i_s] = pow(10,log10(double(osi->variation_factor_CV)) * double(2. * i_s / (osi->n_scales_CV - 1) - 1));}
    else {rel_scale_CV[i_s] = 1;}
    scale_CV[i_s] = rel_scale_CV[i_s] * osi->scale_ren;
  }
  

  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    //  if ((x_m == 0 || x_m == osi->n_moments)){
      ofstream outfile_subprocess;
      ////////////////////////////////////////////////////////////////
      //  subprocess output for CV variation (osi->min_qTcut value)  //
      ////////////////////////////////////////////////////////////////
      logger.newLine(LOG_DEBUG);
      logger << LOG_DEBUG << "subprocess output for CV variation (osi->min_qTcut value)" << endl;
      for (int i_p = 1; i_p < subprocess.size(); i_p++){
	/*
	if (x_m == osi->n_moments){filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot.CV." + subprocess[i_p] + ".dat";}
	else {filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".CV." + subprocess[i_p] + ".dat";}
	*/
	filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot" + ygeneric->infix_name_moment[i_m] + ".CV." + subprocess[i_p] + ".dat";
	logger << LOG_INFO << "XXX filename = " << filename << endl;
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
 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}










void summary_contribution::output_result_check_CV(){ 
  //string result_moment, int x_m
  Logger logger("summary_contribution::output_check_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    //  if ((x_m == 0 || x_m == osi->n_moments)){// || (infix_order_contribution == "LO")){
    //    logger << LOG_DEBUG << "n_subprocess = " << n_subprocess << endl;
    //    logger << LOG_DEBUG << "result_subprocess.size() = " << result_subprocess.size() << endl;
    //    logger << LOG_DEBUG << "deviation_subprocess.size() = " << deviation_subprocess.size() << endl;

    ofstream outfile_check;
    string filename;
    string signame = "XS_" + infix_order_contribution;
    string dsigname = "dXS_" + infix_order_contribution;
    
    //////////////////////////////////////////
    //  subprocess-wise check output begin  //
    //////////////////////////////////////////
    filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + ylist->resultdirectory + "/" + infix_contribution + "/check" + ygeneric->infix_name_moment[i_m] + "." + infix_order_contribution + ".dat";
    logger << LOG_INFO << "XXX filename = " << filename << endl;
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
    /*
    vector<double> rel_scale_CV(osi->n_scales_CV);
    vector<double> scale_CV(osi->n_scales_CV);
    logger << LOG_DEBUG << "osi->n_scales_CV = " << osi->n_scales_CV << endl;
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      if (osi->n_scales_CV > 1){rel_scale_CV[i_s] = pow(10,log10(double(osi->variation_factor_CV)) * double(2. * i_s / (osi->n_scales_CV - 1) - 1));}
      else {rel_scale_CV[i_s] = 1;}
      scale_CV[i_s] = rel_scale_CV[i_s] * osi->scale_ren;
    }
    */
  }
 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



