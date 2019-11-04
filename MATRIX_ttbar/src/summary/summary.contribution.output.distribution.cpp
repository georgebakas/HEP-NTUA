#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_contribution::output_distribution_overview_qTcut_TSV(){
  Logger logger("summary_contribution::output_distribution_overview_qTcut_TSV");
  logger << LOG_DEBUG << "called" << endl;

  stringstream header_ss;
  header_ss << left
	    << setw(16) << "bin"
	    << setw(16) << "qTcut"
	    << setw(25) << "subprocess"
	    << right 
	    << setw(36) << "dXS +- err (subprocess)"
	    << "   "
	    << right << setw(10) << "chi2_dof" 
	    << "   " 
	    << right << setw(5) << "n_run"
	    << " " 
	    << right << setw(5) << "[total]"
	    << "   " 
	    << setw(9) << "(removed)"
	    << "   "
	    << setw(12) << right << "events / bin";

  stringstream temp_ss_separation_one;
  stringstream temp_ss_separation_two;
  int counter_separation = 16 + 16 + 25 + 16 + 4 + 16 + 3 + 10 + 3 + 5 + 3 + 9;
  for (int i_x = 0; i_x < counter_separation; i_x++){
    temp_ss_separation_one << "-";
    temp_ss_separation_two << "=";
  }

  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    if (!osi->switch_distribution_TSV[i_s]){continue;}

    if (i_s != osi->no_reference_TSV){
      logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << " != " << osi->no_reference_TSV << " = osi->no_reference_TSV" << endl;
      continue;
    }

    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){

	if (i_r != osi->no_scale_ren_reference_TSV || i_f != osi->no_scale_fact_reference_TSV){
	  logger << LOG_DEBUG_VERBOSE << "i_r = " << i_r << " != " << osi->no_scale_ren_reference_TSV << " = no_scale_ren_reference_TSV" << "  ||  " << endl;
	  logger << LOG_DEBUG_VERBOSE << "i_f = " << i_f << " != " << osi->no_scale_fact_reference_TSV << " = osi->no_scale_fact_reference_TSV" << endl;
	  continue;
	}

	for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	    if (!ygeneric->switch_output_distribution[i_d]){continue;}
	    string filename_distribution_overview = "overview." + osi->extended_distribution[i_d].xdistribution_name + "." + infix_order_contribution + ".dat";
	    string filepath_distribution_overview = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_overview;
	    logger << LOG_DEBUG << "filename_distribution_overview = " << filename_distribution_overview << endl;
	    logger << LOG_DEBUG << "filepath_distribution_overview = " << filepath_distribution_overview << endl;
	    
	    ofstream out_overviewfile;
	    out_overviewfile.open(filepath_distribution_overview.c_str(), ofstream::out | ofstream::trunc);  
	    out_overviewfile << header_ss.str() << endl;
	    out_overviewfile << endl;

	    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	      int counter_empty_bin = 0;
	      for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
		// check if bin is empty (to reduce output):
		if (xsubprocess[0].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f] == 0.){counter_empty_bin++; continue;}

		int temp_size_result = int(log10(abs(osi->unit_factor_distribution * xsubprocess[0].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f])));
		for (int i_p = 0; i_p < subprocess.size(); i_p++){
		  int temp_size_result_c = int(log10(abs(osi->unit_factor_distribution * xsubprocess[i_p].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f])));
		  if (temp_size_result_c > temp_size_result){temp_size_result = temp_size_result_c;}
		}
		int temp_size_deviation = int(log10(abs(osi->unit_factor_distribution * xsubprocess[0].distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f])));
		if (osi->unit_factor_distribution * xsubprocess[0].distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f] >= 1.){temp_size_deviation++;}

		for (int i_p = 0; i_p < xsubprocess.size(); i_p++){
		  int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_distribution * xsubprocess[i_p].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f])));
		  int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_distribution * xsubprocess[i_p].distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f])));
		  if (osi->unit_factor_distribution * xsubprocess[i_p].distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f] >= 1.){temp_size_deviation_subprocess++;}
		  if (xsubprocess[i_p].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f] == 0){temp_size_result_subprocess = 0;}
		  if (xsubprocess[i_p].distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f] == 0){temp_size_deviation_subprocess = 0;}
		  int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
		  int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
		  if (xsubprocess[i_p].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f] == 0.){setw_result = 1; setw_deviation = 1;}
		  else if (abs(osi->unit_factor_distribution * xsubprocess[i_p].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f]) < 1.){setw_result--;}


		  if (extended_directory.size() > 1 || i_p == 0){
		    out_overviewfile << left << noshowpoint 
				     << setw(16) << setprecision(8) << osi->extended_distribution[i_d].bin_edge[i_b];
		    if (active_qTcut){out_overviewfile << setw(16) << setprecision(8) << osi->value_qTcut_distribution[x_q];}
		    else {out_overviewfile << setw(16) << "independent";}
		    if (i_p == 0){out_overviewfile << setw(25) << infix_order_contribution;}
		    else {out_overviewfile << setw(25) << xsubprocess[i_p].name;}
		    out_overviewfile << showpoint << right 
				     << setw(16) << setprecision(setw_result) 
				     << osi->unit_factor_distribution * xsubprocess[i_p].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f]
				     << " +- "
				     << setw(16) << setprecision(setw_deviation)
				     << osi->unit_factor_distribution * xsubprocess[i_p].distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f];
		    if (i_p > 0){
		      out_overviewfile << "   "
				       << right << showpoint << setprecision(4) << setw(10) << xsubprocess[i_p].distribution_chi2_TSV[i_d][i_b][x_q][i_s][i_r][i_f] //xsubprocess[i_p].chi2_TSV[i_m][i_q][i_s][i_r][i_f]
			;
		    }
		    out_overviewfile << endl;
		  }
		  if (i_p > 0){
		    if (xsubprocess[i_p].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f] != 0.){
		    for (int i_m = 0; i_m < extended_directory.size(); i_m++){
		      if (extended_directory.size() > 1){
			out_overviewfile << left << noshowpoint 
					 << setw(16) << ""
					 << setw(16) << ""
					 << setw(23) << right << extended_directory_name[i_m]
					 << setw(2) << "";
		      }
		      else {
			out_overviewfile << left << noshowpoint 
					 << setw(16) << setprecision(8) << osi->extended_distribution[i_d].bin_edge[i_b];
			if (active_qTcut){out_overviewfile << setw(16) << setprecision(8) << osi->value_qTcut_distribution[x_q];}
			else {out_overviewfile << setw(16) << "independent";}
			if (i_p == 0){out_overviewfile << setw(25) << infix_order_contribution;}
			else {out_overviewfile << setw(25) << xsubprocess[i_p].name;}
		      }
			
		      out_overviewfile << showpoint << right 
				       << setw(16) << setprecision(setw_result) 
				       << osi->unit_factor_distribution * xsubprocess[i_p].distribution_group_result_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f]
				       << " +- "
				       << setw(16) << setprecision(setw_deviation)
				       << osi->unit_factor_distribution * xsubprocess[i_p].distribution_group_deviation_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f];
		      out_overviewfile << "   "
				       << right << showpoint << setprecision(4) << setw(10) << xsubprocess[i_p].distribution_group_chi2_TSV[i_m][i_d][i_b][x_q][i_s][i_r][i_f] //xsubprocess[i_p].chi2_TSV[i_m][i_q][i_s][i_r][i_f]
				       << "   "
				       << setw(5) << right << xsubprocess[i_p].distribution_group_counter_nonzero_run_qTcut_TSV[i_m][i_d][i_b][x_q]
			///				     << setw(5) << right << xsubprocess[i_p].counter_nonzero_distribution_run_qTcut_TSV[i_d][i_b][x_q]
				       << " ["
				       << setw(5) << right << xsubprocess[i_p].distribution_group_counter_run_TSV[i_m]
				       << "]   "
				       << "   (" << setw(4) << right << xsubprocess[i_p].distribution_group_counter_removal_run_qTcut_TSV[i_m][i_d][i_b][x_q] << ")"
		      ///				     << "   (" << setw(4) << right << xsubprocess[i_p].counter_remove_distribution_run_qTcut_TSV[i_d][i_b][x_q] << ")";
				       << "   "
				       << setw(12) << right << xsubprocess[i_p].distribution_group_N_binwise_TSV[i_m][i_d][i_b][x_q]
				       << endl;
		    }
		    }
		    if (extended_directory.size() > 1){
		      out_overviewfile << temp_ss_separation_one.str() << endl;
		    }
		    
		  }
		  
		  if (i_p == 0){
		    out_overviewfile << temp_ss_separation_one.str() << endl;
		  }
		  
		}
		out_overviewfile << endl;
	      }
	      if (counter_empty_bin < selection_n_qTcut){out_overviewfile << temp_ss_separation_two.str() << endl << endl;}
	    }
	    out_overviewfile.close();
	    logger << LOG_DEBUG_VERBOSE << "filepath_distribution_overview = " << filepath_distribution_overview << "   finished." << endl;
	    
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}






void summary_contribution::output_distribution_overview_qTcut_CV(){
  Logger logger("summary_contribution::output_distribution_overview_qTcut_CV");
  logger << LOG_DEBUG << "called" << endl;

  /////////////////////////////////////////////////////////////////////////////
  //  output of distributions from respective contribution for CV variation  //
  /////////////////////////////////////////////////////////////////////////////

  stringstream header_ss;
  header_ss << left
	    << setw(16) << "bin"
	    << setw(16) << "qTcut"
	    << setw(25) << "subprocess"
	    << right 
	    << setw(34) << "dXS +- err (subprocess)"
	    << "   "
	    << right << setw(10) << "chi2_dof" 
	    << "   " 
	    << right << setw(5) << "n_run"
	    << "   " 
	    << setw(9) << "(removed)";

  stringstream temp_ss_separation_one;
  stringstream temp_ss_separation_two;
  int counter_separation = 16 + 16 + 25 + 16 + 4 + 16 + 3 + 10 + 3 + 5 + 3 + 9;
  for (int i_x = 0; i_x < counter_separation; i_x++){
    temp_ss_separation_one << "-";
    temp_ss_separation_two << "=";
  }

  int i_q = 0;
  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      string filename_distribution_overview = "overview." + osi->extended_distribution[i_d].xdistribution_name + "." + infix_order_contribution + ".dat";
      string filepath_distribution_overview = ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_overview;
      logger << LOG_DEBUG << "filename_distribution_overview = " << filename_distribution_overview << endl;
      logger << LOG_DEBUG << "filepath_distribution_overview = " << filepath_distribution_overview << endl;
      
      ofstream out_overviewfile;
      out_overviewfile.open(filepath_distribution_overview.c_str(), ofstream::out | ofstream::trunc);  
      out_overviewfile << header_ss.str() << endl;
      out_overviewfile << endl;
      
      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	int counter_empty_bin = 0;
	// check if bin is empty (to reduce output):
	if (xsubprocess[0].distribution_result_CV[i_d][i_b][i_s] == 0.){counter_empty_bin++; continue;}
	
	int temp_size_result = int(log10(abs(osi->unit_factor_distribution * xsubprocess[0].distribution_result_CV[i_d][i_b][i_s])));
	for (int i_p = 0; i_p < subprocess.size(); i_p++){
	  int temp_size_result_c = int(log10(abs(osi->unit_factor_distribution * xsubprocess[i_p].distribution_result_CV[i_d][i_b][i_s])));
	  if (temp_size_result_c > temp_size_result){temp_size_result = temp_size_result_c;}
	}
	int temp_size_deviation = int(log10(abs(osi->unit_factor_distribution * xsubprocess[0].distribution_deviation_CV[i_d][i_b][i_s])));
	if (osi->unit_factor_distribution * xsubprocess[0].distribution_deviation_CV[i_d][i_b][i_s] >= 1.){temp_size_deviation++;}

	for (int i_p = 0; i_p < xsubprocess.size(); i_p++){
	  int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_distribution * xsubprocess[i_p].distribution_result_CV[i_d][i_b][i_s])));
	  int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_distribution * xsubprocess[i_p].distribution_deviation_CV[i_d][i_b][i_s])));
	  if (osi->unit_factor_distribution * xsubprocess[i_p].distribution_deviation_CV[i_d][i_b][i_s] >= 1.){temp_size_deviation_subprocess++;}
	  if (xsubprocess[i_p].distribution_result_CV[i_d][i_b][i_s] == 0){temp_size_result_subprocess = 0;}
	  if (xsubprocess[i_p].distribution_deviation_CV[i_d][i_b][i_s] == 0){temp_size_deviation_subprocess = 0;}
	  int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
	  int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
	  if (xsubprocess[i_p].distribution_result_CV[i_d][i_b][i_s] == 0.){setw_result = 1;setw_deviation = 1;}
	  else if (abs(osi->unit_factor_distribution * xsubprocess[i_p].distribution_result_CV[i_d][i_b][i_s]) < 1.){setw_result--;}



	  
	  if (extended_directory.size() > 1 || i_p == 0){
	    out_overviewfile << left << noshowpoint 
			     << setw(16) << setprecision(8) << osi->extended_distribution[i_d].bin_edge[i_b];
	    if (active_qTcut){out_overviewfile << setw(16) << setprecision(8) << osi->value_qTcut_distribution[i_q];}
	    else {out_overviewfile << setw(16) << "independent";}
	    if (i_p == 0){out_overviewfile << setw(25) << infix_order_contribution;}
	    else {out_overviewfile << setw(25) << xsubprocess[i_p].name;}
	    out_overviewfile << showpoint << right 
			     << setw(16) << setprecision(setw_result)
			     << osi->unit_factor_distribution * xsubprocess[i_p].distribution_result_CV[i_d][i_b][i_s]
			     << " +- "
			     << setw(16) << setprecision(setw_deviation)
			     << osi->unit_factor_distribution * xsubprocess[i_p].distribution_deviation_CV[i_d][i_b][i_s];
	    if (i_p > 0){
	      out_overviewfile << "   "
			       << right << showpoint << setprecision(4) << setw(10) << xsubprocess[i_p].distribution_chi2_CV[i_d][i_b][i_s]
		;
	      /*
			       << "   "
			       << setw(5) << right << xsubprocess[i_p].counter_nonzero_distribution_run_CV[i_d][i_b]
			       << "   "
			       << "   (" << setw(4) << right << xsubprocess[i_p].counter_remove_distribution_run_CV[i_d][i_b]
			       << ")";
	      */
	    }
	    out_overviewfile << endl;
	  }

	  if (i_p > 0){
	    if (xsubprocess[i_p].distribution_result_CV[i_d][i_b][i_s] != 0.){
	      for (int i_m = 0; i_m < extended_directory.size(); i_m++){
		if (extended_directory.size() > 1){
		  out_overviewfile << left << noshowpoint 
				   << setw(16) << ""
				   << setw(16) << ""
				   << setw(23) << right << extended_directory_name[i_m]
				   << setw(2) << "";
		}
		else {
		  out_overviewfile << left << noshowpoint 
				   << setw(16) << setprecision(8) << osi->extended_distribution[i_d].bin_edge[i_b];
		  if (active_qTcut){out_overviewfile << setw(16) << setprecision(8) << osi->value_qTcut_distribution[i_q];}
		  else {out_overviewfile << setw(16) << "independent";}
		  if (i_p == 0){out_overviewfile << setw(25) << infix_order_contribution;}
		  else {out_overviewfile << setw(25) << xsubprocess[i_p].name;}
		}
		
		out_overviewfile << showpoint << right 
				 << setw(16) << setprecision(setw_result) 
				 << osi->unit_factor_distribution * xsubprocess[i_p].distribution_group_result_CV[i_m][i_d][i_b][i_s]
				 << " +- "
				 << setw(16) << setprecision(setw_deviation)
				 << osi->unit_factor_distribution * xsubprocess[i_p].distribution_group_deviation_CV[i_m][i_d][i_b][i_s];
		out_overviewfile << "   "
				 << right << showpoint << setprecision(4) << setw(10) << xsubprocess[i_p].distribution_group_chi2_CV[i_m][i_d][i_b][i_s] //xsubprocess[i_p].chi2_TSV[i_m][i_q][i_s][i_r][i_f]
				 << "   "
				 << setw(5) << right << xsubprocess[i_p].distribution_group_counter_nonzero_run_CV[i_m][i_d][i_b]
				 << "   "
				 << "   (" << setw(4) << right << xsubprocess[i_p].distribution_group_counter_removal_run_CV[i_m][i_d][i_b] << ")"
				 << "   "
				 << setw(12) << right << xsubprocess[i_p].distribution_group_N_binwise_CV[i_m][i_d][i_b]
				 << endl;
	      }
	    }
	    if (extended_directory.size() > 1){
	      out_overviewfile << temp_ss_separation_one.str() << endl;
	    }
	    
	  }



	  
	  if (i_p == 0){
	    out_overviewfile << temp_ss_separation_one.str() << endl;
	  }
	}
	out_overviewfile << endl;
      }
      out_overviewfile.close();
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}






void summary_contribution::output_distribution_norm_qTcut_TSV(){
  Logger logger("summary_contribution::output_distribution_norm_qTcut_TSV");
  logger << LOG_DEBUG << "called" << endl;

  //////////////////////////////////////////////////////////////////////////////
  //  output of distributions from respective contribution for TSV variation  //
  //////////////////////////////////////////////////////////////////////////////

  for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      if (!osi->switch_distribution_TSV[i_s]){continue;}
      //      if (osi->switch_distribution_TSV[i_s] == 0){continue;}
      //      if (osi->switch_distribution_TSV[i_s] != 0){
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  //	  if (osi->switch_distribution_TSV[i_s] != 0){
	  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	      if (!ygeneric->switch_output_distribution[i_d]){continue;}
	      /*
	      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
		if (distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] != distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f]){
		  cout << "distribution_result_TSV[" << i_g << "][" << i_d << "][" << i_b << "][" << x_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] << " +- " << distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] << endl;
		  //		  exit(1);
		}
	      }
	      */
	      if (i_d < osi->dat.size()){
		
		/////////////////////////////////////////
		//  singly-differential distributions  //
		/////////////////////////////////////////
		
		string filename_distribution_plot = "norm." + osi->extended_distribution[i_d].xdistribution_name + ".." + infix_order_contribution  + ".dat";
		string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_plot;
		logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
		ofstream out_plotfile;
		out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
		for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
		  out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[i_b] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] << endl;
		}
		out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[osi->extended_distribution[i_d].n_bins] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][x_q][i_s][i_r][i_f] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][x_q][i_s][i_r][i_f] << endl;// endpoint
		out_plotfile.close();
		logger << LOG_DEBUG_VERBOSE << "singly-differential output finished" << endl;
	      }
	      else {
		/*
		/////////////////////////////////////////
		//  doubly-differential distributions  //
		/////////////////////////////////////////
		
		int i_ddd = i_d - osi->dat.size();
		logger << LOG_DEBUG_VERBOSE << "osi->dddat[" << i_ddd << "].d_1 = " << setw(2) << osi->dddat[i_ddd].d_1 << "   " << setw(20) << osi->dddat[i_ddd].distribution_1.xdistribution_name << "   " << setw(20) << osi->dddat[i_ddd].distribution_1.xdistribution_type << endl;
		logger << LOG_DEBUG_VERBOSE << "osi->dddat[" << i_ddd << "].d_2 = " << setw(2) << osi->dddat[i_ddd].d_2 << "   " << setw(20) << osi->dddat[i_ddd].distribution_2.xdistribution_name << "   " << setw(20) << osi->dddat[i_ddd].distribution_1.xdistribution_type << endl;
		logger << LOG_DEBUG_VERBOSE << "double-differential output to be defined" << endl;
		
		if (osi->dddat[i_ddd].distribution_1.xdistribution_type == "pTratio" ||
		    osi->dddat[i_ddd].distribution_1.xdistribution_type == "multiplicity"){
		  
		  //////////////////////////////////////////////////////////////////////////////////////////////////////
		  //  sum over re-normalized pTratio-distibution to reconstruct the singly-differential distribution  //
		  //////////////////////////////////////////////////////////////////////////////////////////////////////
		  
		  stringstream filename_distribution_plot_ss;
		  filename_distribution_plot_ss << "norm.rec." << osi->dddat[i_ddd].distribution_1.xdistribution_name << "-" << osi->dddat[i_ddd].distribution_2.xdistribution_name << ".." << infix_order_contribution << ".dat";
		  string filename_distribution_plot = filename_distribution_plot_ss.str();
		  string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_plot;
		  logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
		  ofstream out_plotfile;
		  out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
		  for (int i_b2 = 0; i_b2 < osi->dddat[i_ddd].distribution_2.n_bins; i_b2++){
		    double temp_result = 0;
		    double temp_deviation = 0;
		    for (int i_b1 = 0; i_b1 < osi->dddat[i_ddd].distribution_1.n_bins; i_b1++){
		      int i_b = i_b1 * osi->dddat[i_ddd].distribution_2.n_bins + i_b2;
		      temp_result += distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_1.bin_width[i_b1];
		      temp_deviation += pow(distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_1.bin_width[i_b1], 2);
		    }
		    temp_deviation = sqrt(temp_deviation);
		    out_plotfile << setw(10) << osi->dddat[i_ddd].distribution_2.bin_edge[i_b2] << setw(16) << setprecision(8) << osi->unit_factor_distribution * temp_result << setw(16) << setprecision(8) << osi->unit_factor_distribution * temp_deviation << endl;
		  }
		  out_plotfile.close();
		  
		  /////////////////////////////////////////////////////////////
		  //  state re-normalized distibutions for each pTratio-bin  //
		  /////////////////////////////////////////////////////////////
		  
		  for (int i_b1 = 0; i_b1 < osi->dddat[i_ddd].distribution_1.n_bins; i_b1++){
		    stringstream filename_distribution_plot_ss;
		    filename_distribution_plot_ss << "norm.split." << osi->dddat[i_ddd].name << "_" << osi->dddat[i_ddd].distribution_1.bin_edge[i_b1] << "-" << osi->dddat[i_ddd].distribution_1.bin_edge[i_b1 + 1] << ".." << infix_order_contribution << ".dat";
		    string filename_distribution_plot = filename_distribution_plot_ss.str();
		    string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_plot;
		    logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
		    ofstream out_plotfile;
		    out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
		    for (int i_b2 = 0; i_b2 < osi->dddat[i_ddd].distribution_2.n_bins; i_b2++){
		      int i_b = i_b1 * osi->dddat[i_ddd].distribution_2.n_bins + i_b2;
		      out_plotfile << setw(10) << osi->dddat[i_ddd].distribution_2.bin_edge[i_b2] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_1.bin_width[i_b1] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_1.bin_width[i_b1] << endl;
		    }
		    int i_b = i_b1 * osi->dddat[i_ddd].distribution_2.n_bins + (osi->dddat[i_ddd].distribution_2.n_bins - 1);
		    out_plotfile << setw(10) << osi->dddat[i_ddd].distribution_2.bin_edge[osi->dddat[i_ddd].distribution_2.n_bins] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_2.bin_width[osi->dddat[i_ddd].distribution_2.n_bins - 1] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_2.bin_width[osi->dddat[i_ddd].distribution_2.n_bins - 1] << endl;
		  }
		  out_plotfile.close();
		  
		  ////////////////////////////////////////////////////////////////////
		  //  state re-normalized distibutions for all pTratio <>-binnings  //
		  ////////////////////////////////////////////////////////////////////
		  
		  for (int i_b1 = 0; i_b1 < osi->dddat[i_ddd].distribution_1.n_bins; i_b1++){
		    stringstream filename_distribution_plot_lt_ss;
		    stringstream filename_distribution_plot_ge_ss;
		    filename_distribution_plot_lt_ss << "norm.two." << osi->dddat[i_ddd].name << "_lt" << osi->dddat[i_ddd].distribution_1.bin_edge[i_b1] << ".." << infix_order_contribution << ".dat";
		    filename_distribution_plot_ge_ss << "norm.two." << osi->dddat[i_ddd].name << "_ge" << osi->dddat[i_ddd].distribution_1.bin_edge[i_b1] << ".." << infix_order_contribution << ".dat";
		    string filename_distribution_plot_lt = filename_distribution_plot_lt_ss.str();
		    string filename_distribution_plot_ge = filename_distribution_plot_ge_ss.str();
		    string filepath_distribution_plot_lt = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_plot_lt;
		    string filepath_distribution_plot_ge = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_plot_ge;
		    //		    output_dddistribution_lt_ge(i_ddd, i_b1, x_q, i_s, i_r, i_f, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], filepath_distribution_plot_lt, filepath_distribution_plot_ge, osi->unit_factor_distribution, oset);
		  }
		}
		*/
	      }
	    }
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}





void summary_contribution::output_distribution_plot_qTcut_TSV(){
  Logger logger("summary_contribution::output_distribution_plot_qTcut_TSV");
  logger << LOG_DEBUG << "called" << endl;

  //////////////////////////////////////////////////////////////////////////////
  //  output of distributions from respective contribution for TSV variation  //
  //////////////////////////////////////////////////////////////////////////////

  for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      if (!osi->switch_distribution_TSV[i_s]){continue;}
      //      if (osi->switch_distribution_TSV[i_s] == 0){continue;}
      //      if (osi->switch_distribution_TSV[i_s] != 0){
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  //	  if (osi->switch_distribution_TSV[i_s] != 0){
	  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	      if (!ygeneric->switch_output_distribution[i_d]){continue;}
	      /*
	      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
		if (distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] != distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f]){
		  cout << "distribution_result_TSV[" << i_g << "][" << i_d << "][" << i_b << "][" << x_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] << " +- " << distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] << endl;
		  //		  exit(1);
		}
	      }
	      */
	      if (i_d < osi->dat.size()){
		
		/////////////////////////////////////////
		//  singly-differential distributions  //
		/////////////////////////////////////////
		
		string filename_distribution_plot = "plot." + osi->extended_distribution[i_d].xdistribution_name + ".." + infix_order_contribution  + ".dat";
		string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_plot;
		logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
		ofstream out_plotfile;
		out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
		for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
		  out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[i_b]
			       << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[i_b] 
			       << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[i_b] 
			       << endl;
		}
		out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[osi->extended_distribution[i_d].n_bins]
			     << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][x_q][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[osi->extended_distribution[i_d].n_bins - 1] 
			     << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][x_q][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[osi->extended_distribution[i_d].n_bins - 1] 
			     << endl;// endpoint
		out_plotfile.close();
		logger << LOG_DEBUG_VERBOSE << "singly-differential output finished" << endl;
	      }
	      else {
		
		/////////////////////////////////////////
		//  doubly-differential distributions  //
		/////////////////////////////////////////
		
		int i_ddd = i_d - osi->dat.size();
		logger << LOG_DEBUG_VERBOSE << "osi->dddat[" << i_ddd << "].d_1 = " << setw(2) << osi->dddat[i_ddd].d_1 << "   " << setw(20) << osi->dddat[i_ddd].distribution_1.xdistribution_name << "   " << setw(20) << osi->dddat[i_ddd].distribution_1.xdistribution_type << endl;
		logger << LOG_DEBUG_VERBOSE << "osi->dddat[" << i_ddd << "].d_2 = " << setw(2) << osi->dddat[i_ddd].d_2 << "   " << setw(20) << osi->dddat[i_ddd].distribution_2.xdistribution_name << "   " << setw(20) << osi->dddat[i_ddd].distribution_1.xdistribution_type << endl;
		logger << LOG_DEBUG_VERBOSE << "double-differential output to be defined" << endl;
		
		if (osi->dddat[i_ddd].distribution_1.xdistribution_type == "pTratio" ||
		    osi->dddat[i_ddd].distribution_1.xdistribution_type == "multiplicity"){
		  
		  //////////////////////////////////////////////////////////////////////////////////////////////////////
		  //  sum over re-normalized pTratio-distibution to reconstruct the singly-differential distribution  //
		  //////////////////////////////////////////////////////////////////////////////////////////////////////
		  
		  stringstream filename_distribution_plot_ss;
		  filename_distribution_plot_ss << "plot.rec." << osi->dddat[i_ddd].distribution_1.xdistribution_name << "-" << osi->dddat[i_ddd].distribution_2.xdistribution_name << ".." << infix_order_contribution << ".dat";
		  string filename_distribution_plot = filename_distribution_plot_ss.str();
		  string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_plot;
		  logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
		  ofstream out_plotfile;
		  out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
		  for (int i_b2 = 0; i_b2 < osi->dddat[i_ddd].distribution_2.n_bins; i_b2++){
		    double temp_result = 0;
		    double temp_deviation = 0;
		    for (int i_b1 = 0; i_b1 < osi->dddat[i_ddd].distribution_1.n_bins; i_b1++){
		      int i_b = i_b1 * osi->dddat[i_ddd].distribution_2.n_bins + i_b2;
		      temp_result += distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_1.bin_width[i_b1];
		      temp_deviation += pow(distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_1.bin_width[i_b1], 2);
		    }
		    temp_deviation = sqrt(temp_deviation);
		    out_plotfile << setw(10) << osi->dddat[i_ddd].distribution_2.bin_edge[i_b2] << setw(16) << setprecision(8) << osi->unit_factor_distribution * temp_result << setw(16) << setprecision(8) << osi->unit_factor_distribution * temp_deviation << endl;
		  }
		  out_plotfile.close();
		  
		  /////////////////////////////////////////////////////////////
		  //  state re-normalized distibutions for each pTratio-bin  //
		  /////////////////////////////////////////////////////////////
		  
		  for (int i_b1 = 0; i_b1 < osi->dddat[i_ddd].distribution_1.n_bins; i_b1++){
		    stringstream filename_distribution_plot_ss;
		    filename_distribution_plot_ss << "plot.split." << osi->dddat[i_ddd].name << "_" << osi->dddat[i_ddd].distribution_1.bin_edge[i_b1] << "-" << osi->dddat[i_ddd].distribution_1.bin_edge[i_b1 + 1] << ".." << infix_order_contribution << ".dat";
		    string filename_distribution_plot = filename_distribution_plot_ss.str();
		    string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_plot;
		    logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
		    ofstream out_plotfile;
		    out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
		    for (int i_b2 = 0; i_b2 < osi->dddat[i_ddd].distribution_2.n_bins; i_b2++){
		      int i_b = i_b1 * osi->dddat[i_ddd].distribution_2.n_bins + i_b2;
		      out_plotfile << setw(10) << osi->dddat[i_ddd].distribution_2.bin_edge[i_b2] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_1.bin_width[i_b1] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_1.bin_width[i_b1] << endl;
		    }
		    int i_b = i_b1 * osi->dddat[i_ddd].distribution_2.n_bins + (osi->dddat[i_ddd].distribution_2.n_bins - 1);
		    out_plotfile << setw(10) << osi->dddat[i_ddd].distribution_2.bin_edge[osi->dddat[i_ddd].distribution_2.n_bins] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_2.bin_width[osi->dddat[i_ddd].distribution_2.n_bins - 1] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] * osi->dddat[i_ddd].distribution_2.bin_width[osi->dddat[i_ddd].distribution_2.n_bins - 1] << endl;
		  }
		  out_plotfile.close();
		  
		  ////////////////////////////////////////////////////////////////////
		  //  state re-normalized distibutions for all pTratio <>-binnings  //
		  ////////////////////////////////////////////////////////////////////
		  
		  for (int i_b1 = 0; i_b1 < osi->dddat[i_ddd].distribution_1.n_bins; i_b1++){
		    stringstream filename_distribution_plot_lt_ss;
		    stringstream filename_distribution_plot_ge_ss;
		    filename_distribution_plot_lt_ss << "plot.two." << osi->dddat[i_ddd].name << "_lt" << osi->dddat[i_ddd].distribution_1.bin_edge[i_b1] << ".." << infix_order_contribution << ".dat";
		    filename_distribution_plot_ge_ss << "plot.two." << osi->dddat[i_ddd].name << "_ge" << osi->dddat[i_ddd].distribution_1.bin_edge[i_b1] << ".." << infix_order_contribution << ".dat";
		    string filename_distribution_plot_lt = filename_distribution_plot_lt_ss.str();
		    string filename_distribution_plot_ge = filename_distribution_plot_ge_ss.str();
		    string filepath_distribution_plot_lt = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_plot_lt;
		    string filepath_distribution_plot_ge = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/" + filename_distribution_plot_ge;
		    //		    output_dddistribution_lt_ge(i_ddd, i_b1, x_q, i_s, i_r, i_f, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], filepath_distribution_plot_lt, filepath_distribution_plot_ge, osi->unit_factor_distribution, oset);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




void summary_contribution::output_distribution_plot_CV(){
  Logger logger("summary_contribution::output_distribution_plot_CV");
  logger << LOG_DEBUG << "called" << endl;


  vector<double> n_bin_distribution(osi->extended_distribution.size());
  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    n_bin_distribution[i_d] = osi->extended_distribution[i_d].n_bins;
  }
  // can most likely be just replaced by n_bin_distribution_modasym or osi->extended_distribution[i_d].n_bins !!!

  ////////////////////////////
  //  output for subgroups  // (only 0 = all)
  ////////////////////////////
  // too slow // ofstream out_checksum;
  // too slow // string filename_checksum = vs + "" + ygeneric->final_resultdirectory + "/" + scalename + "checksum." + distdat.xdistribution_name + "." + infix_order_contribution + observable + ".dat";
  // too slow // checksum.open(filename_checksum.c_str(), ofstream::out | ofstream::trunc);  
  //  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
  for (int i_g = 0; i_g < 1; i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	//	logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
	if ((osi->extended_distribution[i_d].xdistribution_name).substr(0, 4) == "asym"){
	  double result_fw = 0.;
	  double result_bw = 0.;
	  double deviation_fw = 0.;
	  double deviation_bw = 0.;
	  ofstream out_asymfile;
	  string filename_asym = ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s] + "/" + osi->extended_distribution[i_d].xdistribution_name + "." + infix_order_contribution + ".dat";
	  logger << LOG_DEBUG_VERBOSE << "filename_asym = " << filename_asym << endl;
	  out_asymfile.open(filename_asym.c_str(), ofstream::out | ofstream::trunc);  
	  if (osi->extended_distribution[i_d].xdistribution_type == "asym1" || osi->extended_distribution[i_d].xdistribution_type == "asym2"){
	    ofstream out_plotfile;
	    string filename_plot = ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s] + "/plot." + osi->extended_distribution[i_d].xdistribution_name + "." + infix_order_contribution + ".dat";
	    logger << LOG_DEBUG_VERBOSE << "infix_order_contribution = " << infix_order_contribution << endl;
	    logger << LOG_DEBUG_VERBOSE << "filename_plot = " << filename_plot << endl;
	    out_plotfile.open(filename_plot.c_str(), ofstream::out | ofstream::trunc);  
	    // too slow // checksum << setw(16) << "interval" << " " << setw(16) << "backward cs." << setw(16) << "forward cs." << setw(16) << "cross section" << endl;
	    for (int i_b = 0; i_b < n_bin_distribution[i_d]; i_b++){
	      if (osi->extended_distribution[i_d].xdistribution_type == "asym1"){
		result_fw = distribution_result_CV[0][i_d][n_bin_distribution[i_d] + i_b][i_s] * osi->fakeasymfactor[i_d];
		result_bw = distribution_result_CV[0][i_d][i_b][i_s] * osi->fakeasymfactor[i_d];
		deviation_fw = distribution_deviation_CV[0][i_d][n_bin_distribution[i_d] + i_b][i_s] * osi->fakeasymfactor[i_d];
		deviation_bw = distribution_deviation_CV[0][i_d][i_b][i_s] * osi->fakeasymfactor[i_d];
	      }
	      else if (osi->extended_distribution[i_d].xdistribution_type == "asym2"){
		result_fw = distribution_result_CV[0][i_d][2 * i_b + 1][i_s] * osi->fakeasymfactor[i_d];
		result_bw = distribution_result_CV[0][i_d][2 * i_b + 0][i_s] * osi->fakeasymfactor[i_d];
		deviation_fw = distribution_deviation_CV[0][i_d][2 * i_b + 1][i_s] * osi->fakeasymfactor[i_d];
		deviation_bw = distribution_deviation_CV[0][i_d][2 * i_b + 0][i_s] * osi->fakeasymfactor[i_d];
	      }
	      out_asymfile << setw(25) << "backw. cross section: " << "   [" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << ";" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b + 1] << "]   " << setw(16) << setprecision(8) << result_bw * osi->extended_distribution[i_d].step << setw(16) << setprecision(8) << deviation_bw * osi->extended_distribution[i_d].step << endl;
	      out_asymfile << setw(25) << "forw. cross section: " << "   [" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << ";" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b + 1] << "]   " << setw(16) << setprecision(8) << result_fw * osi->extended_distribution[i_d].step << setw(16) << setprecision(8) << deviation_fw * osi->extended_distribution[i_d].step << endl;
	      out_asymfile << setw(25) << "forw.-backw. asymmetry: " << "   [" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << ";" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b + 1] << "]   " << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	      out_plotfile << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	      // too slow // checksum << "[ " << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << "; " << setw(5) << setprecision(4) << bin_edge[i_d][n_bin_distribution[i_d]] << "] ";
	      double sum_result_fw = 0., sum_result_bw = 0.;
	      for (int c = i_b; c < n_bin_distribution[i_d]; c++){
		if (osi->extended_distribution[i_d].xdistribution_type == "asym1"){
		  sum_result_fw += distribution_result_CV[0][i_d][n_bin_distribution[i_d] + c][i_s] * osi->fakeasymfactor[i_d];
		  sum_result_bw += distribution_result_CV[0][i_d][c][i_s] * osi->fakeasymfactor[i_d];
		}
		else if (osi->extended_distribution[i_d].xdistribution_type == "asym2"){
		  sum_result_fw += distribution_result_CV[0][i_d][2 * c + 1][i_s] * osi->fakeasymfactor[i_d];
		  sum_result_bw += distribution_result_CV[0][i_d][2 * c + 0][i_s] * osi->fakeasymfactor[i_d];
		}
	      }
	      // too slow // checksum << setw(16) << setprecision(8) << sum_result_bw * osi->extended_distribution[i_d].step << setw(16) << setprecision(8) << sum_result_fw * osi->extended_distribution[i_d].step << setw(16) << setprecision(8) << (sum_result_fw + sum_result_bw) * osi->extended_distribution[i_d].step << endl;
	    }
	    out_plotfile << setw(5) << setprecision(4) << osi->extended_distribution[i_d].end << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	    out_asymfile.close();
	    out_plotfile.close();
	  }
	  else {
	    result_fw = distribution_result_CV[0][i_d][1][i_s] * osi->fakeasymfactor[i_d];
	    result_bw = distribution_result_CV[0][i_d][0][i_s] * osi->fakeasymfactor[i_d];
	    deviation_fw = distribution_deviation_CV[0][i_d][1][i_s] * osi->fakeasymfactor[i_d];
	    deviation_bw = distribution_deviation_CV[0][i_d][0][i_s] * osi->fakeasymfactor[i_d];
	    out_asymfile << setw(25) << "backw. cross section: " << setw(16) << setprecision(8) << result_bw << setw(16) << setprecision(8) << deviation_bw << endl;
	    out_asymfile << setw(25) << "forw. cross section: " << setw(16) << setprecision(8) << result_fw << setw(16) << setprecision(8) << deviation_fw << endl;
	    out_asymfile << setw(25) << "forw.-backw. asymmetry: " << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	    // too slow // checksum << setw(16) << "backward cs." << setw(16) << "forward cs." << setw(16) << "cross section" << endl;
	    // too slow // checksum << setw(16) << setprecision(8) << result_bw << setw(16) << setprecision(8) << result_fw << setw(16) << setprecision(8) << (result_fw + result_bw) << endl;
	    out_asymfile.close();
	  }
	}
	else {
	  ofstream out_plotfile;
	  string filename_plot = ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s] + "/" + ylist->resultdirectory + "/" + infix_contribution + "/plot." + osi->extended_distribution[i_d].xdistribution_name + "." + infix_order_contribution + ".dat";
	  logger << LOG_DEBUG_VERBOSE << "filename_plot = " << filename_plot << endl;
	  out_plotfile.open(filename_plot.c_str(), ofstream::out | ofstream::trunc);  
	  // too slow // checksum << setw(16) << "interval" << " " << setw(16) << "cross section" << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << "distribution_result_CV.size() = " << distribution_result_CV.size() << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << "distribution_result_CV[" << i_g << "].size() = " << distribution_result_CV[i_g].size() << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << "distribution_result_CV[" << i_g << "][" << i_d << "].size() = " << distribution_result_CV[i_g][i_d].size() << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << "distribution_deviation_CV.size() = " << distribution_deviation_CV.size() << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << "distribution_deviation_CV[" << i_g << "].size() = " << distribution_deviation_CV[i_g].size() << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << "distribution_deviation_CV[" << i_g << "][" << i_d << "].size() = " << distribution_deviation_CV[i_g][i_d].size() << endl;
	  for (int i_b = 0; i_b < distribution_result_CV[i_g][i_d].size(); i_b++){
	    //	    logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	    //	  for (int i_b = 0; i_b < distribution_result_CV[i_g].size(); i_b++){
	    //	    logger << LOG_DEBUG_VERBOSE << "distribution_result_CV[" << i_g << "][" << i_d << "][" << i_b << "].size() = " << distribution_result_CV[i_g][i_d][i_b].size() << endl;
	    out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[i_b] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_CV[i_g][i_d][i_b][i_s] / osi->extended_distribution[i_d].bin_width[i_b] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_CV[i_g][i_d][i_b][i_s] / osi->extended_distribution[i_d].bin_width[i_b] << endl;
	    // norm ->	    out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[i_b] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_CV[i_g][i_d][i_b][i_s] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_CV[i_g][i_d][i_b][i_s] << endl;



	    // too slow // checksum << "[ " << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << "; " << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[n_bin_distribution[i_d]] << "] ";
	    double sum_result = 0.;
	    for (int j_b = i_b; j_b < n_bin_distribution[i_d]; j_b++){sum_result += distribution_result_CV[i_g][i_d][j_b][i_s];}
	    // norm ->	    for (int j_b = i_b; j_b < n_bin_distribution[i_d]; j_b++){sum_result += distribution_result_CV[i_g][i_d][j_b][i_s] / osi->extended_distribution[i_d].bin_width[i_b];}
	    // too slow // checksum << setw(16) << setprecision(8) << sum_result * osi->extended_distribution[i_d].step << endl;
	  }
	  //	  logger << LOG_DEBUG_VERBOSE << "osi->extended_distribution[i_d].end = " << osi->extended_distribution[i_d].end << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << "distribution_result_CV[" << i_g << "][" << i_d << "][" << distribution_result_CV[i_g][i_d].size() - 1 << "][" << i_s << "] = " << distribution_result_CV[i_g][i_d][distribution_result_CV[i_g][i_d].size() - 1][i_s] << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << "distribution_deviation_CV[" << i_g << "][" << i_d << "][" << distribution_deviation_CV[i_g][i_d].size() - 1 << "][" << i_s << "] = " << distribution_deviation_CV[i_g][i_d][distribution_deviation_CV[i_g][i_d].size() - 1][i_s] << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << " = " <<  << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << " = " <<  << endl;
	  out_plotfile << setw(10) << osi->extended_distribution[i_d].end << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_CV[i_g][i_d][distribution_result_CV[i_g][i_d].size() - 1][i_s] / osi->extended_distribution[i_d].bin_width[osi->extended_distribution[i_d].n_bins - 1] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_CV[i_g][i_d][distribution_deviation_CV[i_g][i_d].size() - 1][i_s] / osi->extended_distribution[i_d].bin_width[osi->extended_distribution[i_d].n_bins - 1] << endl; // endpoint
	  // norm ->	    out_plotfile << setw(10) << osi->extended_distribution[i_d].end << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_CV[i_g][i_d][distribution_result_CV[i_g][i_d].size() - 1][i_s] << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_CV[i_g][i_d][distribution_deviation_CV[i_g][i_d].size() - 1][i_s] << endl; // endpoint
	  //	  logger << LOG_DEBUG_VERBOSE << "end i_s = " << i_s << endl;
	  out_plotfile.close();
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_contribution::output_distribution_asymmetry_CV(){
  Logger logger("summary_contribution::output_distribution_asymmetry_CV");
  logger << LOG_DEBUG << "called" << endl;

  /*
  logger << LOG_DEBUG_VERBOSE << "asym calculation begins" << endl;
 
  ///////////////////////////////
  //  output for subprocesses  //
  ///////////////////////////////
  for (int i_p = 0; i_p < xsubprocess.size(); i_p++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      ///// too slow // ofstream out_checksum;
      ///// too slow // string filename_checksum = ylist->resultdirectory + "/" + infix_order_contribution + "/" + scalename + "checksum." + osi->extended_distribution[i_d].xdistribution_name + "." + subprocess[i_p] + ".dat";
      ///// too slow // checksum.open(filename_checksum.c_str(), ofstream::out | ofstream::trunc);  
      if ((osi->extended_distribution[i_d].xdistribution_name).substr(0, 4) == "asym"){
	double result_fw, result_bw, deviation_fw, deviation_bw;
	ofstream out_asymfile;
	string filename_asym = ylist->resultdirectory + "/" + infix_order_contribution + "/" + scalename + osi->extended_distribution[i_d].xdistribution_name + "." + subprocess[i_p] + ".dat";
	logger << LOG_DEBUG_VERBOSE << "filename_asym = " << filename_asym << endl;
	///      out_asymfile.open(filename_asym.c_str(), ofstream::out | ofstream::trunc);  
	if (osi->extended_distribution[i_d].xdistribution_type == "asym1" || osi->extended_distribution[i_d].xdistribution_type == "asym2"){
	  ///	ofstream out_plotfile;
	  ///	string filename_plot = ylist->resultdirectory + "/" + infix_order_contribution + "/" + scalename + "plot." + osi->extended_distribution[i_d].xdistribution_name + "." + subprocess[i_p] + ".dat";
	  ///	logger << LOG_DEBUG_VERBOSE << "filename_plot = " << filename_plot << endl;
	  ///	out_plotfile.open(filename_plot.c_str(), ofstream::out | ofstream::trunc);  
	  ///// too slow // checksum << setw(16) << "interval" << " " << setw(16) << "backward cs." << setw(16) << "forward cs." << setw(16) << "cross section" << endl;
	  for (int i_b = 0; i_b < n_bin_distribution[i_d]; i_b++){
	    if (osi->extended_distribution[i_d].xdistribution_type == "asym1"){
	      //xsubprocess[i_p].distribution_result_CV[i_d][i_b][i_s]
	    result_fw = xsubprocess[i_p].distribution_result_CV[i_d][n_bin_distribution[i_d] + i_b][i_s] * osi->fakeasymfactor[i_d];
	    result_bw = xsubprocess[i_p].distribution_result_CV[i_d][i_b][i_s] * osi->fakeasymfactor[i_d];
	    deviation_fw = xsubprocess[i_p].distribution_deviation_CV[i_d][n_bin_distribution[i_d] + i_b][i_s] * osi->fakeasymfactor[i_d];
	    deviation_bw = xsubprocess[i_p].distribution_deviation_CV[i_d][i_b][i_s] * osi->fakeasymfactor[i_d];
	  }
	  else if (osi->extended_distribution[i_d].xdistribution_type == "asym2"){
	    result_fw = xsubprocess[i_p].distribution_result_CV[i_d][2 * i_b + 1][i_s] * osi->fakeasymfactor[i_d];
	    result_bw = xsubprocess[i_p].distribution_result_CV[i_d][2 * i_b + 0][i_s] * osi->fakeasymfactor[i_d];
	    deviation_fw = xsubprocess[i_p].distribution_deviation_CV[i_d][2 * i_b + 1][i_s] * osi->fakeasymfactor[i_d];
	    deviation_bw = xsubprocess[i_p].distribution_deviation_CV[i_d][2 * i_b + 0][i_s] * osi->fakeasymfactor[i_d];
	  }
	  ///	  out_asymfile << setw(25) << "backw. cross section: ";
	  ///	  out_asymfile << "   [" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << ";" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b + 1] << "]   " << setw(16) << setprecision(8) << result_bw * osi->extended_distribution[i_d].step << setw(16) << setprecision(8) << deviation_bw * osi->extended_distribution[i_d].step;// * (bin_edge[i_d][n_bin_distribution[i_d]] - bin_edge[i_d][0]);
	  ///	  out_asymfile << endl;
	  ///	  out_asymfile << setw(25) << "forw. cross section: ";
	  ///	  out_asymfile << "   [" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << ";" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b + 1] << "]   " << setw(16) << setprecision(8) << result_fw * osi->extended_distribution[i_d].step << setw(16) << setprecision(8) << deviation_fw * osi->extended_distribution[i_d].step;// * (bin_edge[i_d][n_bin_distribution[i_d]] - bin_edge[i_d][0]);
	  ///	  out_asymfile << endl;
	  ///	  out_asymfile << setw(25) << "forw.-backw. asymmetry: ";
	  ///	  out_asymfile << "   [" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << ";" << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b + 1] << "]   " << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2);
	  ///	  out_asymfile << endl;
	  ///	  out_plotfile << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	  ///// too slow // checksum << "[ " << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << "; " << setw(5) << setprecision(4) << bin_edge[i_d][n_bin_distribution[i_d]] << "] ";
	  double sum_result_fw = 0., sum_result_bw = 0.;
	  for (int j_b = i_b; j_b < n_bin_distribution[i_d]; j_b++){
	    if (osi->extended_distribution[i_d].xdistribution_type == "asym1"){
	      sum_result_fw += xsubprocess[i_p].distribution_result_CV[i_d][n_bin_distribution[i_d] + j_b][i_s] * osi->fakeasymfactor[i_d];
	      sum_result_bw += xsubprocess[i_p].distribution_result_CV[i_d][j_b][i_s] * osi->fakeasymfactor[i_d];
	    }
	    else if (osi->extended_distribution[i_d].xdistribution_type == "asym2"){
	      sum_result_fw += xsubprocess[i_p].distribution_result_CV[i_d][2 * j_b + 1][i_s] * osi->fakeasymfactor[i_d];
	      sum_result_bw += xsubprocess[i_p].distribution_result_CV[i_d][2 * j_b + 0][i_s] * osi->fakeasymfactor[i_d];
	    }
	  }
	  ///// too slow // checksum << setw(16) << setprecision(8) << sum_result_bw * osi->extended_distribution[i_d].step << setw(16) << setprecision(8) << sum_result_fw * osi->extended_distribution[i_d].step << setw(16) << setprecision(8) << (sum_result_fw + sum_result_bw) * osi->extended_distribution[i_d].step << endl;
	}
	///	out_plotfile << setw(5) << setprecision(4) << osi->extended_distribution[i_d].end << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	///	out_plotfile.close();
	///	out_asymfile.close();
      }
      else {
	result_fw = xsubprocess[i_p].distribution_result_CV[i_d][1][i_s] * osi->fakeasymfactor[i_d];
	result_bw = xsubprocess[i_p].distribution_result_CV[i_d][0][i_s] * osi->fakeasymfactor[i_d];
	deviation_fw = xsubprocess[i_p].distribution_deviation_CV[i_d][1][i_s] * osi->fakeasymfactor[i_d];
	deviation_bw = xsubprocess[i_p].distribution_deviation_CV[i_d][0][i_s] * osi->fakeasymfactor[i_d];
	///	out_asymfile << setw(25) << "backw. cross section: " << setw(16) << setprecision(8) << result_bw << setw(16) << setprecision(8) << deviation_bw << endl;
	///	out_asymfile << setw(25) << "forw. cross section: " << setprecision(8) << result_fw << setw(16) << setprecision(8) << deviation_fw << endl;
	///	out_asymfile << setw(25) << "forw.-backw. asymmetry: " << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	///// too slow // checksum << setw(16) << "backward cs." << setw(16) << "forward cs." << setw(16) << "cross section" << endl;
	///// too slow // checksum << setw(16) << setprecision(8) << result_bw << setw(16) << setprecision(8) << result_fw << setw(16) << setprecision(8) << (result_fw + result_bw) << endl;
	///	out_asymfile.close();
      }
    }
    else {
      ofstream out_plotfile;
      string filename_plot = ylist->resultdirectory + "/" + infix_order_contribution + "/" + scalename + "plot." + osi->extended_distribution[i_d].xdistribution_name + "." + subprocess[i_p] + ".dat";
      logger << LOG_DEBUG_VERBOSE << "filename_plot = " << filename_plot << endl;
      out_plotfile.open(filename_plot.c_str(), ofstream::out | ofstream::trunc);  
      ///// too slow // checksum << setw(16) << "interval" << " " << setw(16) << "cross section" << endl;
      for (int i_b = 0; i_b < n_bin_distribution[i_d]; i_b++){
	out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[i_b] << setw(16) << setprecision(8) << xsubprocess[i_p].distribution_result_CV[i_d][i_b][i_s] << setw(16) << setprecision(8) << xsubprocess[i_p].distribution_deviation_CV[i_d][i_b][i_s] << endl;
	///// too slow // checksum << "[ " << setw(5) << setprecision(4) << osi->extended_distribution[i_d].bin_edge[i_b] << "; " << setw(5) << setprecision(4) << bin_edge[i_d][n_bin_distribution[i_d]] << "] ";
	///// too slow // double sum_result = 0.;
	///// too slow // for (int j_b = i_b; j_b < n_bin_distribution[i_d]; j_b++){sum_result += xsubprocess[i_p].distribution_result_CV[i_d][j_b][i_s];}
	///// too slow // checksum << setw(16) << setprecision(8) << sum_result * osi->extended_distribution[i_d].step << endl;
      }
      ///      out_plotfile << setw(10) << osi->extended_distribution[i_d].end << setw(16) << setprecision(8) << xsubprocess[i_p].distribution_result_CV[i_d][n_bin_distribution[i_d] - 1] << setw(16) << setprecision(8) << xsubprocess[i_p].distribution_deviation_CV[i_d][n_bin_distribution[i_d] - 1] << endl; // endpoint
      out_plotfile.close();
    }
    ///// too slow // checksum.close();
  }
  logger << LOG_DEBUG_VERBOSE << "asym calculation ends" << endl;
*/  

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



