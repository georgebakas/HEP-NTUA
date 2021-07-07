#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_list::output_distribution_overview_TSV(){
  Logger logger("summary_list::output_distribution_overview_TSV");
  logger << LOG_DEBUG << "called" << endl;

  stringstream header_ss;
  header_ss << left
	    << setw(16) << "bin"
	    << setw(16) << "qTcut"
	    << setw(25) << "contribution"
	    << right << setw(34) << "dXS +- err (contribution)";

  stringstream temp_ss_separation_one;
  stringstream temp_ss_separation_two;
  int counter_separation = 16 + 16 + 25 + 16 + 4 + 16;
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
	    string filename_distribution_overview = "overview." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[0].infix_order_contribution + ".dat";
	    string filepath_distribution_overview = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + resultdirectory + "/" + filename_distribution_overview;
	    
	    logger << LOG_DEBUG << "filename_distribution_overview = " << filename_distribution_overview << endl;
	    logger << LOG_DEBUG << "filepath_distribution_overview = " << filepath_distribution_overview << endl;
	    
	    ofstream out_overviewfile;
	    out_overviewfile.open(filepath_distribution_overview.c_str(), ofstream::out | ofstream::trunc);  
	    out_overviewfile << header_ss.str() << endl;
	    out_overviewfile << endl;
	    
	    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){

	      if (distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] != 0.){
		out_overviewfile << left
				 << setw(16) << setprecision(8) << noshowpoint 
				 << osi->extended_distribution[i_d].bin_edge[i_b];
		if (active_qTcut){out_overviewfile << setw(16) << "extrapolation";}
		else {out_overviewfile << setw(16) << "independent";}
		//		  out_overviewfile << setw(25) << xcontribution[0].infix_order_contribution;
		out_overviewfile << setw(25) << resultdirectory;
		out_overviewfile << showpoint << right
				 << setw(16) << setprecision(9) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f]
				 << " +- " 
				 << setw(16) << setprecision(9) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] 
				 << endl;
		out_overviewfile << temp_ss_separation_one.str() << endl;
		out_overviewfile << endl;
	      }
	      
	      int counter_empty_bin = 0;
	      for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
		// check if bin is empty (to reduce output):
		if (xcontribution[0].distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] == 0.){counter_empty_bin++; continue;}

		int temp_size_result = int(log10(abs(osi->unit_factor_distribution * xcontribution[0].distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f])));
		//		int temp_size_result_list = temp_size_result; // ???
		for (int i_c = 0; i_c < xcontribution.size(); i_c++){
		  int y_q = 0;
		  if (xcontribution[i_c].active_qTcut){y_q = x_q;}
		  int temp_size_result_c = int(log10(abs(osi->unit_factor_distribution * xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f])));
		  if (temp_size_result_c > temp_size_result){temp_size_result = temp_size_result_c;}
		}
		int temp_size_deviation = int(log10(abs(osi->unit_factor_distribution * xcontribution[0].distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f])));
		if (osi->unit_factor_distribution * xcontribution[0].distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] >= 1.){temp_size_deviation++;}

		for (int i_c = 0; i_c < xcontribution.size(); i_c++){
		  int y_q = 0;
		  if (xcontribution[i_c].active_qTcut){y_q = x_q;}
		  
		  int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_distribution * xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f])));
		  int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f])));
		  if (osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] >= 1.){temp_size_deviation_subprocess++;}
		  if (xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] == 0){temp_size_result_subprocess = 0;}
		  if (xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] == 0){temp_size_deviation_subprocess = 0;}
		  int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
		  int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
		  if (xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] == 0.){setw_result = 1; setw_deviation = 1;}
		  else if (abs(osi->unit_factor_distribution * xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f]) < 1.){setw_result--;}
		  
		  out_overviewfile << left
				   << setw(16) << setprecision(8) << noshowpoint 
				   << osi->extended_distribution[i_d].bin_edge[i_b];
		  if (xcontribution[i_c].active_qTcut){out_overviewfile << noshowpoint << setw(16) << osi->value_qTcut_distribution[x_q];}
		  else {out_overviewfile << setw(16) << "independent";}
		  if (i_c == 0){out_overviewfile << setw(25) << xcontribution[i_c].infix_order_contribution;}
		  else {out_overviewfile << setw(25) << xcontribution[i_c].infix_contribution;}
		  out_overviewfile << showpoint << right
				   << setw(16) << setprecision(setw_result) << osi->unit_factor_distribution * xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f]
				   << " +- " 
				   << setw(16) << setprecision(setw_deviation) << osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] 
		    /*
				   << "    "
				   << " wr " << setw_result 
				   << " wd " << setw_deviation
		    */
				   << endl;
		  if (i_c == 0){
		    out_overviewfile << temp_ss_separation_one.str() << endl;
		  }
		}
		out_overviewfile << endl;
	      }
	      if (counter_empty_bin < selection_n_qTcut){out_overviewfile << temp_ss_separation_two.str() << endl << endl;}
	    }
	    out_overviewfile.close();
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}        



void summary_list::output_distribution_overview_CV(){
  Logger logger("summary_list::output_distribution_overview_CV");
  logger << LOG_DEBUG << "called" << endl;

  stringstream header_ss;
  header_ss << left
	    << setw(16) << "bin"
	    << setw(16) << "qTcut"
	    << setw(25) << "contribution"
	    << right << setw(34) << "dXS +- err (contribution)";

  stringstream temp_ss_separation_one;
  stringstream temp_ss_separation_two;
  int counter_separation = 16 + 16 + 25 + 16 + 4 + 16;
  for (int i_x = 0; i_x < counter_separation; i_x++){
    temp_ss_separation_one << "-";
    temp_ss_separation_two << "=";
  }

  int i_q = 0;
  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	string filename_distribution_overview = "overview." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[0].infix_order_contribution + ".dat";
	string filepath_distribution_overview = ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s] + "/" + resultdirectory + "/" + filename_distribution_overview;
	logger << LOG_DEBUG << "filename_distribution_overview = " << filename_distribution_overview << endl;
	logger << LOG_DEBUG << "filepath_distribution_overview = " << filepath_distribution_overview << endl;
	
	ofstream out_overviewfile;
	out_overviewfile.open(filepath_distribution_overview.c_str(), ofstream::out | ofstream::trunc);  
	out_overviewfile << header_ss.str() << endl;
	out_overviewfile << endl;
	    
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	  /*
	      if (distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] != 0.){
		out_overviewfile << left
				 << setw(16) << setprecision(8) << noshowpoint 
				 << osi->extended_distribution[i_d].bin_edge[i_b];
		if (active_qTcut){out_overviewfile << setw(16) << "extrapolation";}
		else {out_overviewfile << setw(16) << "independent";}
		//		  out_overviewfile << setw(25) << xcontribution[0].infix_order_contribution;
		out_overviewfile << setw(25) << resultdirectory;
		out_overviewfile << showpoint << right
				 << setw(16) << setprecision(9) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f]
				 << " +- " 
				 << setw(16) << setprecision(9) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] 
				 << endl;
		out_overviewfile << temp_ss_separation_one.str() << endl;
		out_overviewfile << endl;
	      }
	  */  
	  int counter_empty_bin = 0;
	  // check if bin is empty (to reduce output):
	  if (xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s] == 0.){counter_empty_bin++; continue;}
	  
	  int temp_size_result = int(log10(abs(osi->unit_factor_distribution * xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s])));
	  //		int temp_size_result_list = temp_size_result; // ???
	  for (int i_c = 0; i_c < xcontribution.size(); i_c++){
	    //	int y_q = 0;
	    //		if (xcontribution[i_c].active_qTcut){y_q = x_q;}
	    int temp_size_result_c = int(log10(abs(osi->unit_factor_distribution * xcontribution[i_c].distribution_result_CV[i_g][i_d][i_b][i_s])));
	    if (temp_size_result_c > temp_size_result){temp_size_result = temp_size_result_c;}
	  }
	  int temp_size_deviation = int(log10(abs(osi->unit_factor_distribution * xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s])));
	  if (osi->unit_factor_distribution * xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s] >= 1.){temp_size_deviation++;}
	  
	  for (int i_c = 0; i_c < xcontribution.size(); i_c++){
	    //		int y_q = 0;
	    //		if (xcontribution[i_c].active_qTcut){y_q = x_q;}
	    
	    int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_distribution * xcontribution[i_c].distribution_result_CV[i_g][i_d][i_b][i_s])));
	    int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_CV[i_g][i_d][i_b][i_s])));
	    if (osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_CV[i_g][i_d][i_b][i_s] >= 1.){temp_size_deviation_subprocess++;}
	    if (xcontribution[i_c].distribution_result_CV[i_g][i_d][i_b][i_s] == 0){temp_size_result_subprocess = 0;}
	    if (xcontribution[i_c].distribution_deviation_CV[i_g][i_d][i_b][i_s] == 0){temp_size_deviation_subprocess = 0;}
	    int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
	    int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
	    if (xcontribution[i_c].distribution_result_CV[i_g][i_d][i_b][i_s] == 0.){setw_result = 1;setw_deviation = 1;}
	    else if (abs(osi->unit_factor_distribution * xcontribution[i_c].distribution_result_CV[i_g][i_d][i_b][i_s]) < 1.){setw_result--;}
	    
	    out_overviewfile << left
			     << setw(16) << setprecision(8) << noshowpoint 
			     << osi->extended_distribution[i_d].bin_edge[i_b];
	    if (xcontribution[i_c].active_qTcut){out_overviewfile << noshowpoint << setw(16) << osi->value_qTcut[i_q];}
	    else {out_overviewfile << setw(16) << "independent";}
	    if (i_c == 0){out_overviewfile << setw(25) << xcontribution[i_c].infix_order_contribution;}
	    else {out_overviewfile << setw(25) << xcontribution[i_c].infix_contribution;}
	    out_overviewfile << showpoint << right
			     << setw(16) << setprecision(setw_result) << osi->unit_factor_distribution * xcontribution[i_c].distribution_result_CV[i_g][i_d][i_b][i_s]
			     << " +- " 
			     << setw(16) << setprecision(setw_deviation) << osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_CV[i_g][i_d][i_b][i_s] 
	      /*
			     << "    "
			     << " wr " << setw_result 
			     << " wd " << setw_deviation
	      */
			     << endl;
	    if (i_c == 0){
	      out_overviewfile << temp_ss_separation_one.str() << endl;
	    }
	  }
	  out_overviewfile << endl;
	  if (counter_empty_bin < selection_n_qTcut){out_overviewfile << temp_ss_separation_two.str() << endl << endl;}
	}
	out_overviewfile.close();
      }
    }
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}        










void summary_list::output_contribution_distribution_norm_TSV(){
  Logger logger("summary_list::output_contribution_distribution_norm_TSV");
  logger << LOG_DEBUG << "called" << endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  output of distributions from respective contribution for TSV variation, at extrapolation qTcut -> 0  //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    if (!osi->switch_distribution_TSV[i_s]){continue;}
    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	    if (!ygeneric->switch_output_distribution[i_d]){continue;}
	    string filename_distribution_plot = "norm." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[0].infix_order_contribution + ".dat";
	    logger << LOG_DEBUG_VERBOSE << "filename_distribution_plot              = " << filename_distribution_plot << endl;
	    string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + resultdirectory + "/" + filename_distribution_plot;
	    logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
	    ofstream out_plotfile;
	    out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
	    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	      // !!! should be only done only for final output, not for the results of subcontribution, which may be negative !!!
	      //		      if (plotmode == 1){if (distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] < 0){distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] = 0.;}}
	      out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[i_b] 
			   << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] 
			   << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] << endl;
	    }
	    out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[osi->extended_distribution[i_d].n_bins] 
			 << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][i_s][i_r][i_f] 
			 << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][i_s][i_r][i_f] << endl;// endpoint
	    out_plotfile.close();
	    logger << LOG_DEBUG_VERBOSE << "singly-differential output finished" << endl;
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}        



void summary_list::output_contribution_distribution_plot_TSV(){
  Logger logger("summary_list::output_contribution_distribution_plot_TSV");
  logger << LOG_DEBUG << "called" << endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  output of distributions from respective contribution for TSV variation, at extrapolation qTcut -> 0  //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    if (!osi->switch_distribution_TSV[i_s]){continue;}
    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	    if (!ygeneric->switch_output_distribution[i_d]){continue;}
	    string filename_distribution_plot = "plot." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[0].infix_order_contribution + ".dat";
	    logger << LOG_DEBUG_VERBOSE << "filename_distribution_plot              = " << filename_distribution_plot << endl;
	    string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + resultdirectory + "/" + filename_distribution_plot;
	    logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
	    ofstream out_plotfile;
	    out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
	    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	      // !!! should be only done only for final output, not for the results of subcontribution, which may be negative !!!
	      //		      if (plotmode == 1){if (distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] < 0){distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] = 0.;}}
	      out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[i_b] 
			   << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[i_b] 
			   << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[i_b]
			   << endl;
	    }
	    out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[osi->extended_distribution[i_d].n_bins] 
			 << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[osi->extended_distribution[i_d].n_bins - 1] 
			 << setw(16) << setprecision(8) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[osi->extended_distribution[i_d].n_bins - 1] << endl;// endpoint
	    out_plotfile.close();
	    logger << LOG_DEBUG_VERBOSE << "singly-differential output finished" << endl;
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}        



void summary_list::output_contribution_distribution_norm_qTcut_TSV(){
  Logger logger("summary_list::output_contribution_distribution_plot_qTcut_TSV");
  logger << LOG_DEBUG << "called" << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  output of distributions from respective contribution for TSV variation, at selected qTcut values  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
    stringstream qTcut_ss;
    qTcut_ss << "qTcut-" << osi->value_qTcut_distribution[x_q];
    //qTcut_ss << "qTcut" << osi->no_qTcut_distribution[x_q];
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      if (!osi->switch_distribution_TSV[i_s]){continue;}
      //      if (osi->switch_distribution_TSV[i_s]){
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  string directory_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + resultdirectory + "/" + qTcut_ss.str();
	  system_execute(logger, "mkdir " + directory_distribution_plot);

	  //	      for (int j_c = 0; j_c < xcontribution.size(); j_c++){
	  for (int j_c = 0; j_c < 1; j_c++){
	    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
		if (!ygeneric->switch_output_distribution[i_d]){continue;}
		//		    string filename_distribution_plot = "plot." + qTcut_ss.str() + "." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[j_c].infix_order_contribution + ".dat";
		//		    string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + resultdirectory + "/" + filename_distribution_plot;
		string filename_distribution_plot = "norm." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[j_c].infix_order_contribution + ".dat";
		string filepath_distribution_plot = directory_distribution_plot + "/" + filename_distribution_plot;
		logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
		ofstream out_plotfile;
		out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
		for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
		  // !!! should be only done only for final output, not for the results of subcontribution, which may be negative !!!
		  //		      if (plotmode == 1){if (xcontribution[j_c].distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] < 0){xcontribution[j_c].distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] = 0.;}}
		  out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[i_b] << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[j_c].distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[j_c].distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] << endl;
		}
		out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[osi->extended_distribution[i_d].n_bins] << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[j_c].distribution_result_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][x_q][i_s][i_r][i_f] << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[j_c].distribution_deviation_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][x_q][i_s][i_r][i_f] << endl;// endpoint
		out_plotfile.close();
		logger << LOG_DEBUG_VERBOSE << "singly-differential output finished" << endl;
		/*
		    ofstream out_overviewfile;
		    string filename_distribution_overview = "overview." + qTcut_ss.str() + "." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[0].infix_order_contribution + ".dat";
		    logger << LOG_DEBUG_VERBOSE << "filename_distribution_overview = " << filename_distribution_overview << endl;
		    string filepath_distribution_overview = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + resultdirectory + "/" + filename_distribution_overview;
		    logger << LOG_DEBUG_VERBOSE << "filepath_distribution_overview = " << filepath_distribution_overview << endl;
		    out_overviewfile.open(filepath_distribution_overview.c_str(), ofstream::out | ofstream::trunc);  
		    ///		    for (int i_b = 0; i_b < xcontribution[0].distribution_result_CV[i_g][i_d].size(); i_b++){
		    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
		      for (int i_c = 0; i_c < xcontribution.size(); i_c++){
			int y_q = 0;
			if (xcontribution[i_c].active_qTcut){y_q = x_q;}
			string sign = "";
			string minussign = " ";
			if (xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] >= 0){sign = " "; minussign = "";}
			out_overviewfile << left << setw(16) << setprecision(8) << noshowpoint << osi->extended_distribution[i_d].bin_edge[i_b] << setw(20) << left << xcontribution[i_c].infix_order_contribution << sign << showpoint << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] << minussign << "+-  " << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] << endl;
			if (i_c == 0){out_overviewfile << endl;}
		      }
		      out_overviewfile << endl;
		    }
		    out_overviewfile.close();
		*/
	      }
	    }
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}        



void summary_list::output_contribution_distribution_plot_qTcut_TSV(){
  Logger logger("summary_list::output_contribution_distribution_plot_qTcut_TSV");
  logger << LOG_DEBUG << "called" << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  output of distributions from respective contribution for TSV variation, at selected qTcut values  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  logger << LOG_INFO << "resultdirectory = " << resultdirectory << "   selection_n_qTcut = " << selection_n_qTcut << endl;
  for (int i_c = 0; i_c < xcontribution.size(); i_c++){
    logger << LOG_INFO << "xcontribution[" << i_c << "].selection_n_qTcut = " << xcontribution[i_c].selection_n_qTcut << endl;
  }
  for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
    stringstream qTcut_ss;
    qTcut_ss << "qTcut-" << osi->value_qTcut_distribution[x_q];
    //qTcut_ss << "qTcut" << osi->no_qTcut_distribution[x_q];
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      if (!osi->switch_distribution_TSV[i_s]){continue;}
      //      if (osi->switch_distribution_TSV[i_s]){
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  string directory_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + resultdirectory + "/" + qTcut_ss.str();
	  logger << LOG_INFO << "directory_distribution_plot = " << directory_distribution_plot << endl;
	  system_execute(logger, "mkdir " + directory_distribution_plot);

	  //	      for (int j_c = 0; j_c < xcontribution.size(); j_c++){
	  for (int j_c = 0; j_c < 1; j_c++){
	    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
		if (!ygeneric->switch_output_distribution[i_d]){continue;}
		//		    string filename_distribution_plot = "plot." + qTcut_ss.str() + "." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[j_c].infix_order_contribution + ".dat";
		//		    string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + resultdirectory + "/" + filename_distribution_plot;
		string filename_distribution_plot = "plot." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[j_c].infix_order_contribution + ".dat";
		string filepath_distribution_plot = directory_distribution_plot + "/" + filename_distribution_plot;
		logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
		ofstream out_plotfile;
		out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
		for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
		  // !!! should be only done only for final output, not for the results of subcontribution, which may be negative !!!
		  //		      if (plotmode == 1){if (xcontribution[j_c].distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] < 0){xcontribution[j_c].distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] = 0.;}}
		  out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[i_b] 
			       << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[j_c].distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[i_b] 
			       << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[j_c].distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[i_b] << endl;
		}
		out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[osi->extended_distribution[i_d].n_bins] 
			     << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[j_c].distribution_result_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][x_q][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[osi->extended_distribution[i_d].n_bins - 1] 
			     << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[j_c].distribution_deviation_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][x_q][i_s][i_r][i_f] / osi->extended_distribution[i_d].bin_width[osi->extended_distribution[i_d].n_bins - 1] << endl;// endpoint
		out_plotfile.close();
		logger << LOG_DEBUG_VERBOSE << "singly-differential output finished" << endl;
		/*
		    ofstream out_overviewfile;
		    string filename_distribution_overview = "overview." + qTcut_ss.str() + "." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[0].infix_order_contribution + ".dat";
		    logger << LOG_DEBUG_VERBOSE << "filename_distribution_overview = " << filename_distribution_overview << endl;
		    string filepath_distribution_overview = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + resultdirectory + "/" + filename_distribution_overview;
		    logger << LOG_DEBUG_VERBOSE << "filepath_distribution_overview = " << filepath_distribution_overview << endl;
		    out_overviewfile.open(filepath_distribution_overview.c_str(), ofstream::out | ofstream::trunc);  
		    ///		    for (int i_b = 0; i_b < xcontribution[0].distribution_result_CV[i_g][i_d].size(); i_b++){
		    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
		      for (int i_c = 0; i_c < xcontribution.size(); i_c++){
			int y_q = 0;
			if (xcontribution[i_c].active_qTcut){y_q = x_q;}
			string sign = "";
			string minussign = " ";
			if (xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] >= 0){sign = " "; minussign = "";}
			out_overviewfile << left << setw(16) << setprecision(8) << noshowpoint << osi->extended_distribution[i_d].bin_edge[i_b] << setw(20) << left << xcontribution[i_c].infix_order_contribution << sign << showpoint << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] << minussign << "+-  " << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] << endl;
			if (i_c == 0){out_overviewfile << endl;}
		      }
		      out_overviewfile << endl;
		    }
		    out_overviewfile.close();
		*/
	      }
	    }
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}        



void summary_list::output_contribution_distribution_check_TSV(){
  Logger logger("summary_list::output_contribution_distribution_check_TSV");
  logger << LOG_DEBUG << "called" << endl;

  //////////////////////////////////////////////////////////////////////////////
  //  output to check if distributions sum to cross section in TSV variation  //
  //////////////////////////////////////////////////////////////////////////////
  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    if (!osi->switch_distribution_TSV[i_s]){continue;}
    for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
      ofstream out_check;
      stringstream qTcut_ss;
      qTcut_ss << "qTcut" << osi->no_qTcut_distribution[x_q];
      string filename_check = ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/scale.1.1/" + resultdirectory + "/" + "check." + qTcut_ss.str() + ".distribution.txt";		     
      //    string filename_check = ygeneric->final_resultdirectory + "/" + resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" +"check.distribution.txt";		     
      logger << LOG_INFO << "filename_check = " << filename_check << endl;
      out_check.open(filename_check.c_str(), ofstream::out | ofstream::trunc);  

      for (int i_c = 0; i_c < xcontribution.size(); i_c++){
	logger << LOG_INFO << "YYY xcontribution[" << i_c << "].type_contribution = " << xcontribution[i_c].type_contribution << endl;
	int y_q = 0;
	if (xcontribution[i_c].active_qTcut){y_q = x_q;}
	for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
  	  // check output: check if sum over all bins gives a result <= XSection 
	  // print output only at central scales!!!
	  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	    if (!ygeneric->switch_output_distribution[i_d]){continue;}
	    //  if (osi->extended_distribution[i_d].xdistribution_number < 100){
	    if (osi->extended_distribution[i_d].typeCumulative == CUMULATIVE_NONE){
	      out_check << setw(5) << "distribution: " << setw(25) << osi->extended_distribution[i_d].xdistribution_name;
	      double temp_weight = 0.;
	      double temp_deviation = 0.;
	      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
		temp_weight += xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][osi->no_central_scale_ren_TSV[i_s]][osi->no_central_scale_fact_TSV[i_s]];
		temp_deviation += pow(xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][osi->no_central_scale_ren_TSV[i_s]][osi->no_central_scale_fact_TSV[i_s]], 2);
	      }
	      temp_weight = temp_weight * osi->extended_distribution[i_d].step;
	      temp_deviation = sqrt(temp_deviation) * osi->extended_distribution[i_d].step;
	      out_check << setw(16) << setprecision(8) << osi->unit_factor_distribution * temp_weight << " +- " << setw(16) << setprecision(8) << osi->unit_factor_distribution * temp_deviation << endl;
	    }
	    //  else if (osi->extended_distribution[i_d].xdistribution_number < 200){
	    else if (osi->extended_distribution[i_d].typeCumulative == CUMULATIVE_LOWER){
	      out_check << setw(5) << "distribution: " << setw(25) << osi->extended_distribution[i_d].xdistribution_name;
	      out_check << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_result_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][y_q][i_s][osi->no_central_scale_ren_TSV[i_s]][osi->no_central_scale_fact_TSV[i_s]] << " +- " << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][y_q][i_s][osi->no_central_scale_ren_TSV[i_s]][osi->no_central_scale_fact_TSV[i_s]] << endl;
	    }
	    //  else if (osi->extended_distribution[i_d].xdistribution_number < 300){
	    else if (osi->extended_distribution[i_d].typeCumulative == CUMULATIVE_UPPER){
	      out_check << setw(5) << "distribution: " << setw(25) << osi->extended_distribution[i_d].xdistribution_name;
	      out_check << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_result_TSV[i_g][i_d][0][y_q][i_s][osi->no_central_scale_ren_TSV[i_s]][osi->no_central_scale_fact_TSV[i_s]] << " +- " << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][0][y_q][i_s][osi->no_central_scale_ren_TSV[i_s]][osi->no_central_scale_fact_TSV[i_s]] << endl;
	    }
	    //  else if (osi->extended_distribution[i_d].xdistribution_number < 400){
	    else if (osi->extended_distribution[i_d].typeCumulative == CUMULATIVE_BOTH){
	      out_check << setw(5) << "distribution: " << setw(25) << osi->extended_distribution[i_d].xdistribution_name;
	      out_check << setw(5) << "0" << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_result_TSV[i_g][i_d][0][y_q][i_s][osi->no_central_scale_ren_TSV[i_s]][osi->no_central_scale_fact_TSV[i_s]] << " +- " << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][0][y_q][i_s][osi->no_central_scale_ren_TSV[i_s]][osi->no_central_scale_fact_TSV[i_s]] << endl;
	      out_check << setw(5) << "max" << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_result_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][y_q][i_s][osi->no_central_scale_ren_TSV[i_s]][osi->no_central_scale_fact_TSV[i_s]] << " +- " << setw(16) << setprecision(8) << osi->unit_factor_distribution * xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][osi->extended_distribution[i_d].n_bins - 1][y_q][i_s][osi->no_central_scale_ren_TSV[i_s]][osi->no_central_scale_fact_TSV[i_s]] << endl;
	    }
	  }
	}
      }
      out_check.close();
    }
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}        
