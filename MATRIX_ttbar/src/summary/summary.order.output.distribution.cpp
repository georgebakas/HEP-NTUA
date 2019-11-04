#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_order::output_distribution_overview_TSV(){
  Logger logger("summary_order::output_distribution_overview_TSV");
  logger << LOG_DEBUG << "called" << endl;

  stringstream header_ss;
  header_ss << left
	    << setw(16) << "bin"
	    << setw(16) << "qTcut"
	    << setw(25) << "list"
	    << right << setw(34) << "dXS +- err (list)";

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
	    string filename_distribution_overview = "overview." + osi->extended_distribution[i_d].xdistribution_name + "." + resultdirectory + ".dat";
	    string filepath_distribution_overview = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + filename_distribution_overview;
	    
	    logger << LOG_DEBUG_VERBOSE << "filename_distribution_overview = " << filename_distribution_overview << endl;
	    logger << LOG_DEBUG << "OUTPUT        " << "filepath_distribution_overview = " << filepath_distribution_overview << endl;
	    
	    ofstream out_overviewfile;
	    out_overviewfile.open(filepath_distribution_overview.c_str(), ofstream::out | ofstream::trunc);  
	    out_overviewfile << header_ss.str() << endl;
	    out_overviewfile << endl;
	    
	    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	      
	      if (distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] != 0.){

		int temp_size_result = int(log10(abs(osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f])));
		int temp_size_result_order = temp_size_result;
		for (int i_l = 0; i_l < xlist.size(); i_l++){
		  int temp_size_result_c = int(log10(abs(osi->unit_factor_distribution * xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f])));
		  if (temp_size_result_c > temp_size_result){temp_size_result = temp_size_result_c;}
		}
		int temp_size_deviation = int(log10(abs(osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f])));
		if (osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] >= 1.){temp_size_deviation++;}
		int setw_result = 9 + temp_size_result_order - temp_size_result;
		int setw_deviation = 9 + temp_size_deviation - temp_size_result - 1;
		if (distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] == 0.){setw_result = 1; setw_deviation = 1;}
		else if (abs(osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f]) < 1.){setw_result--;}


		
		// original from here on
		out_overviewfile << left
				 << setw(16) << setprecision(8) << noshowpoint 
				 << osi->extended_distribution[i_d].bin_edge[i_b];
		if (active_qTcut){out_overviewfile << setw(16) << "extrapolation";}
		else {out_overviewfile << setw(16) << "independent";}
		//		  out_overviewfile << setw(25) << xcontribution[0].infix_order_contribution;
		out_overviewfile << setw(25) << resultdirectory;
		out_overviewfile << showpoint << right
				 << setw(16) << setprecision(setw_result) << osi->unit_factor_distribution * distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f]
				 << " +- " 
				 << setw(16) << setprecision(setw_deviation) << osi->unit_factor_distribution * distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] 
				 << endl;
		out_overviewfile << temp_ss_separation_one.str() << endl;
		// original until here


		
		for (int i_l = 0; i_l < xlist.size(); i_l++){
		  int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_distribution * xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f])));
		  int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_distribution * xlist[i_l]->distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f])));
		  if (osi->unit_factor_distribution * xlist[i_l]->distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] >= 1.){temp_size_deviation_subprocess++;}
		  if (xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] == 0){temp_size_result_subprocess = 0;}
		  if (xlist[i_l]->distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] == 0){temp_size_deviation_subprocess = 0;}
		  int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
		  int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
		  if (xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] == 0.){setw_result = 1; setw_deviation = 1;}
		  else if (abs(osi->unit_factor_distribution * xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f]) < 1.){setw_result--;}
		  
		  out_overviewfile << left
				   << setw(16) << setprecision(8) << noshowpoint 
				   << osi->extended_distribution[i_d].bin_edge[i_b];
		  if (xlist[i_l]->active_qTcut){out_overviewfile << setw(16) << "extrapolation";}
		  else {out_overviewfile << setw(16) << "independent";}
		  out_overviewfile << setw(25) << xlist[i_l]->xcontribution[0].infix_order_contribution;
		  out_overviewfile << showpoint << right
				   << setw(16) << setprecision(setw_result) << osi->unit_factor_distribution * xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f]
				   << " +- " 
				   << setw(16) << setprecision(setw_deviation) << osi->unit_factor_distribution * xlist[i_l]->distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] 
				   << endl;
		}
		out_overviewfile << endl;

	      }

	      if (active_qTcut){
	      int counter_empty_bin = 0;
	      for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
		// check if bin is empty (to reduce output):
		if (distribution_result_qTcut_TSV[x_q][i_g][i_d][i_b][i_s][i_r][i_f] == 0.){counter_empty_bin++; continue;}
		
		int temp_size_result = int(log10(abs(osi->unit_factor_distribution * distribution_result_qTcut_TSV[x_q][i_g][i_d][i_b][i_s][i_r][i_f])));
		int temp_size_result_order = temp_size_result;
		for (int i_l = 0; i_l < xlist.size(); i_l++){
		  int y_q = 0;
		  if (xlist[i_l]->active_qTcut){y_q = x_q;}
		  int temp_size_result_c = int(log10(abs(osi->unit_factor_distribution * xlist[i_l]->distribution_result_qTcut_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f])));
		  if (temp_size_result_c > temp_size_result){temp_size_result = temp_size_result_c;}
		}
		int temp_size_deviation = int(log10(abs(osi->unit_factor_distribution * distribution_deviation_qTcut_TSV[x_q][i_g][i_d][i_b][i_s][i_r][i_f])));
		if (osi->unit_factor_distribution * distribution_deviation_qTcut_TSV[x_q][i_g][i_d][i_b][i_s][i_r][i_f] >= 1.){temp_size_deviation++;}
		int setw_result = 9 + temp_size_result_order - temp_size_result;
		int setw_deviation = 9 + temp_size_deviation - temp_size_result - 1;
		if (distribution_result_qTcut_TSV[x_q][i_g][i_d][i_b][i_s][i_r][i_f] == 0.){setw_result = 1; setw_deviation = 1;}
		else if (abs(osi->unit_factor_distribution * distribution_result_qTcut_TSV[x_q][i_g][i_d][i_b][i_s][i_r][i_f]) < 1.){setw_result--;}
		
		out_overviewfile << left
				 << setw(16) << setprecision(8) << noshowpoint 
				 << osi->extended_distribution[i_d].bin_edge[i_b];
		//			if (active_qTcut){out_overviewfile << noshowpoint << setw(16) << osi->value_qTcut_distribution[x_q];}
		if (active_qTcut){out_overviewfile << noshowpoint << setw(16) << osi->value_qTcut_distribution[x_q];}
		else {out_overviewfile << setw(16) << "independent";}
		out_overviewfile << setw(25) << resultdirectory;
		out_overviewfile << showpoint << right
				 << setw(16) << setprecision(setw_result) << osi->unit_factor_distribution * distribution_result_qTcut_TSV[x_q][i_g][i_d][i_b][i_s][i_r][i_f]
				 << " +- " 
				 << setw(16) << setprecision(setw_deviation) << osi->unit_factor_distribution * distribution_deviation_qTcut_TSV[x_q][i_g][i_d][i_b][i_s][i_r][i_f] 
		  /*
				 << "    "
				 << " wr " << setw_result 
				 << " wd " << setw_deviation
		  */
				 << endl;
		
		out_overviewfile << temp_ss_separation_one.str() << endl;
		
		for (int i_l = 0; i_l < xlist.size(); i_l++){
		  int y_q = 0;
		  if (xlist[i_l]->active_qTcut){y_q = x_q;}

		  int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_distribution * xlist[i_l]->distribution_result_qTcut_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f])));
		  int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_distribution * xlist[i_l]->distribution_deviation_qTcut_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f])));
		  if (osi->unit_factor_distribution * xlist[i_l]->distribution_deviation_qTcut_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] >= 1.){temp_size_deviation_subprocess++;}
		  if (xlist[i_l]->distribution_result_qTcut_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] == 0){temp_size_result_subprocess = 0;}
		  if (xlist[i_l]->distribution_deviation_qTcut_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] == 0){temp_size_deviation_subprocess = 0;}
		  int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
		  int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
		  if (xlist[i_l]->distribution_result_qTcut_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] == 0.){setw_result = 1; setw_deviation = 1;}
		  else if (abs(osi->unit_factor_distribution * xlist[i_l]->distribution_result_qTcut_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f]) < 1.){setw_result--;}
		  
		  out_overviewfile << left
				   << setw(16) << setprecision(8) << noshowpoint 
				   << osi->extended_distribution[i_d].bin_edge[i_b];
		  //			if (active_qTcut){out_overviewfile << noshowpoint << setw(16) << osi->value_qTcut_distribution[x_q];}
		  if (xlist[i_l]->active_qTcut){out_overviewfile << noshowpoint << setw(16) << osi->value_qTcut_distribution[x_q];}
		  else {out_overviewfile << setw(16) << "independent";}
		  out_overviewfile << setw(25) << xlist[i_l]->xcontribution[0].infix_order_contribution;
		  /*
		    if (i_l == 0){out_overviewfile << setw(25) << xlist[i_l]->xcontribution[0].infix_order_contribution;}
		    else {out_overviewfile << setw(25) << xlist[i_l]->xcontribution[0].infix_contribution;}
		  */
		  out_overviewfile << showpoint << right
				   << setw(16) << setprecision(setw_result) << osi->unit_factor_distribution * xlist[i_l]->distribution_result_qTcut_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f]
				   << " +- " 
				   << setw(16) << setprecision(setw_deviation) << osi->unit_factor_distribution * xlist[i_l]->distribution_deviation_qTcut_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] 
		    /*
				   << "    "
				   << " wr " << setw_result 
				   << " wd " << setw_deviation
		    */
				   << endl;
		}
		out_overviewfile << endl;
	      }
	      if (counter_empty_bin < selection_n_qTcut){out_overviewfile << temp_ss_separation_two.str() << endl << endl;}
	      }
	    }
	    out_overviewfile.close();
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_order::output_distribution_overview_CV(){
  Logger logger("summary_order::output_distribution_overview_CV");
  logger << LOG_DEBUG << "called" << endl;

  stringstream header_ss;
  header_ss << left
	    << setw(16) << "bin"
	    << setw(16) << "qTcut"
	    << setw(25) << "list"
	    << right << setw(34) << "dXS +- err (list)";

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
	string filename_distribution_overview = "overview." + osi->extended_distribution[i_d].xdistribution_name + "." + resultdirectory + ".dat";
	string filepath_distribution_overview = ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s] + "/" + filename_distribution_overview;
	logger << LOG_DEBUG_VERBOSE << "filename_distribution_overview = " << filename_distribution_overview << endl;
	logger << LOG_DEBUG << "filepath_distribution_overview = " << filepath_distribution_overview << endl;
	
	ofstream out_overviewfile;
	out_overviewfile.open(filepath_distribution_overview.c_str(), ofstream::out | ofstream::trunc);  
	out_overviewfile << header_ss.str() << endl;
	out_overviewfile << endl;
	    
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	  int counter_empty_bin = 0;
	  // check if bin is empty (to reduce output):
	  if (distribution_result_CV[i_g][i_d][i_b][i_s] == 0.){counter_empty_bin++; continue;}
	  
	  int temp_size_result = int(log10(abs(osi->unit_factor_distribution * distribution_result_CV[i_g][i_d][i_b][i_s])));
	  int temp_size_result_order = temp_size_result;
	  for (int i_l = 0; i_l < xlist.size(); i_l++){
	    int temp_size_result_c = int(log10(abs(osi->unit_factor_distribution * xlist[i_l]->xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s])));
	    if (temp_size_result_c > temp_size_result){temp_size_result = temp_size_result_c;}
	  }
	  int temp_size_deviation = int(log10(abs(osi->unit_factor_distribution * distribution_deviation_CV[i_g][i_d][i_b][i_s])));
	  if (osi->unit_factor_distribution * distribution_deviation_CV[i_g][i_d][i_b][i_s] >= 1.){temp_size_deviation++;}
 	  int setw_result = 9 + temp_size_result_order - temp_size_result;
	  int setw_deviation = 9 + temp_size_deviation - temp_size_result - 1;
	  if (distribution_result_CV[i_g][i_d][i_b][i_s] == 0.){setw_result = 1; setw_deviation = 1;}
	  else if (abs(osi->unit_factor_distribution * distribution_result_CV[i_g][i_d][i_b][i_s]) < 1.){setw_result--;}
	  
	  out_overviewfile << left
			   << setw(16) << setprecision(8) << noshowpoint 
			   << osi->extended_distribution[i_d].bin_edge[i_b];
	  if (active_qTcut){out_overviewfile << noshowpoint << setw(16) << osi->value_qTcut_distribution[i_q];}
	  else {out_overviewfile << setw(16) << "independent";}
	  out_overviewfile << setw(25) << resultdirectory;
	  out_overviewfile << showpoint << right
			   << setw(16) << setprecision(setw_result) << osi->unit_factor_distribution * distribution_result_CV[i_g][i_d][i_b][i_s]
			   << " +- " 
			   << setw(16) << setprecision(setw_deviation) << osi->unit_factor_distribution * distribution_deviation_CV[i_g][i_d][i_b][i_s] 
			   << endl;
	  
	  out_overviewfile << temp_ss_separation_one.str() << endl;
	  //	  out_overviewfile << "xlist.size() = " << xlist.size() << endl;
	  for (int i_l = 0; i_l < xlist.size(); i_l++){
	    //	    out_overviewfile << "xlist[" << i_l << "] = " << xlist[i_l]->resultdirectory << endl;
	    int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_distribution * xlist[i_l]->xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s])));
	    int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_distribution * xlist[i_l]->xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s])));
	    if (osi->unit_factor_distribution * xlist[i_l]->xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s] >= 1.){temp_size_deviation_subprocess++;}
	    if (xlist[i_l]->xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s] == 0){temp_size_result_subprocess = 0;}
	    if (xlist[i_l]->xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s] == 0){temp_size_deviation_subprocess = 0;}
	    
	    int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
	    int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
	    if (xlist[i_l]->xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s] == 0.){
	      setw_result = 1;
	      setw_deviation = 1;
	    }
	    else if (abs(osi->unit_factor_distribution * xlist[i_l]->xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s]) < 1.){
	      setw_result--;
	    }
	    
	    out_overviewfile << left
			     << setw(16) << setprecision(8) << noshowpoint 
			     << osi->extended_distribution[i_d].bin_edge[i_b];
	    if (xlist[i_l]->active_qTcut){out_overviewfile << noshowpoint << setw(16) << osi->value_qTcut[i_q];}
	    else {out_overviewfile << setw(16) << "independent";}
	    out_overviewfile << setw(25) << xlist[i_l]->resultdirectory;
	    out_overviewfile << showpoint << right
			     << setw(16) << setprecision(setw_result) << osi->unit_factor_distribution * xlist[i_l]->xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s]
			     << " +- " 
			     << setw(16) << setprecision(setw_deviation) << osi->unit_factor_distribution * xlist[i_l]->xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s] 
	      /*
			     << "    "
			     << " wr " << setw_result 
			     << " wd " << setw_deviation
	      */
			     << endl;
	    /*
	    if (i_c == 0){
	      out_overviewfile << temp_ss_separation_one.str() << endl;
	    }
	    */
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




















void summary_order::output_sddistribution_TSV(string & identifier, vector<vector<vector<vector<double> > > > & this_distribution_result_TSV, vector<vector<vector<vector<double> > > > & this_distribution_deviation_TSV, int i_d, string & subdirectory){
  Logger logger("summary_order::output_sddistribution_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  /////////////////////////////////////////
  //  singly-differential distributions  //
  /////////////////////////////////////////
  
  int plotmode = 0;

  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    if (!ygeneric->oset.switch_distribution_TSV[i_s]){continue;}
    for (int i_r = 0; i_r < ygeneric->oset.n_scale_ren_TSV[i_s]; i_r++){
      for (int i_f = 0; i_f < ygeneric->oset.n_scale_fact_TSV[i_s]; i_f++){
	string filename_distribution_plot = identifier + "." + ygeneric->oset.extended_distribution[i_d].xdistribution_name + ".." + resultdirectory + ".dat";
	string filename_full = ygeneric->scalename_TSV[i_s][i_r][i_f] + subdirectory + "/" + filename_distribution_plot;
	logger << LOG_DEBUG_VERBOSE << "filename_full              = " << filename_full << endl;
	ofstream out_plotfile;
	out_plotfile.open(filename_full.c_str(), ofstream::out | ofstream::trunc);  
	double temp_result = 0;
	double temp_deviation = 0;
	for (int i_b = 0; i_b < ygeneric->oset.extended_distribution[i_d].n_bins; i_b++){
	  if (identifier == "norm"){
	    temp_result = this_distribution_result_TSV[i_b][i_s][i_r][i_f];
	    temp_deviation = this_distribution_deviation_TSV[i_b][i_s][i_r][i_f];
	  }
	  else if (identifier == "plot"){
	    temp_result = this_distribution_result_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.extended_distribution[i_d].bin_width[i_b];
	    temp_deviation = this_distribution_deviation_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.extended_distribution[i_d].bin_width[i_b];
	  }
	  if (plotmode == 1){if (temp_result < 0.){temp_result = 0.;}}

	  out_plotfile << setw(10) << ygeneric->oset.extended_distribution[i_d].bin_edge[i_b] 
		       << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_result
		       << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_deviation << endl;
	}
	out_plotfile << setw(10) << ygeneric->oset.extended_distribution[i_d].bin_edge[ygeneric->oset.extended_distribution[i_d].n_bins] 
		     << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_result
		     << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_deviation << endl; 
	out_plotfile.close();
	logger << LOG_DEBUG_VERBOSE << "singly-differential output finished" << endl;
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void summary_order::output_dddistribution_reconstruct_sdd_TSV(string & identifier, vector<vector<vector<vector<double> > > > & this_distribution_result_TSV, vector<vector<vector<vector<double> > > > & this_distribution_deviation_TSV, int i_d, string & subdirectory){
  Logger logger("summary_order::output_dddistribution_reconstruct_sdd_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  /////////////////////////////////////////
  //  doubly-differential distributions  //
  /////////////////////////////////////////
  int i_ddd = i_d - ygeneric->oset.dat.size();

  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    if (!ygeneric->oset.switch_distribution_TSV[i_s]){continue;}
    for (int i_r = 0; i_r < ygeneric->oset.n_scale_ren_TSV[i_s]; i_r++){
      for (int i_f = 0; i_f < ygeneric->oset.n_scale_fact_TSV[i_s]; i_f++){
	/*
	  if (ygeneric->oset.dddat[i_ddd].distribution_1.xdistribution_type == "pTratio" ||
	  ygeneric->oset.dddat[i_ddd].distribution_1.xdistribution_type == "multiplicity" ||
	  ygeneric->oset.dddat[i_ddd].distribution_1.xdistribution_type == "abseta"){
	*/
	if (1 > 0){
	  //////////////////////////////////////////////////////////////////////////////////////////////////////
	  //  sum over re-normalized pTratio-distibution to reconstruct the singly-differential distribution  //
	  //////////////////////////////////////////////////////////////////////////////////////////////////////
	  stringstream filename_distribution_plot_ss;
	  filename_distribution_plot_ss << identifier << ".rec." << ygeneric->oset.dddat[i_ddd].distribution_2.xdistribution_name << ".from." << ygeneric->oset.dddat[i_ddd].name << ".." << resultdirectory << ".dat";
	  string filename_distribution_plot = filename_distribution_plot_ss.str();
	  string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + subdirectory + "/" + filename_distribution_plot;
	  logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
	  ofstream out_plotfile;
	  out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
	  for (int i_b2 = 0; i_b2 < ygeneric->oset.dddat[i_ddd].distribution_2.n_bins; i_b2++){
	    double temp_result = 0;
	    double temp_deviation = 0;
	    for (int i_b1 = 0; i_b1 < ygeneric->oset.dddat[i_ddd].distribution_1.n_bins; i_b1++){
	      int i_b = i_b1 * ygeneric->oset.dddat[i_ddd].distribution_2.n_bins + i_b2;
	      if (identifier == "norm"){
		temp_result += this_distribution_result_TSV[i_b][i_s][i_r][i_f];// / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2];
		temp_deviation += pow(this_distribution_deviation_TSV[i_b][i_s][i_r][i_f], 2);// / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2], 2);
	      }
	      else if (identifier == "plot"){
		temp_result += this_distribution_result_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2];
		temp_deviation += pow(this_distribution_deviation_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2], 2);
	      }
	    }
	    temp_deviation = sqrt(temp_deviation);
	    out_plotfile << setw(10) << ygeneric->oset.dddat[i_ddd].distribution_2.bin_edge[i_b2] 
			 << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_result 
			 << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_deviation << endl;
	  }
	  out_plotfile.close();
	}
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void summary_order::output_dddistribution_split_in_first_sdd_TSV(string & identifier, vector<vector<vector<vector<double> > > > & this_distribution_result_TSV, vector<vector<vector<vector<double> > > > & this_distribution_deviation_TSV, int i_d, string & subdirectory){
  Logger logger("summary_order::output_dddistribution_split_in_first_sdd_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int i_ddd = i_d - ygeneric->oset.dat.size();

  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    if (!ygeneric->oset.switch_distribution_TSV[i_s]){continue;}
    for (int i_r = 0; i_r < ygeneric->oset.n_scale_ren_TSV[i_s]; i_r++){
      for (int i_f = 0; i_f < ygeneric->oset.n_scale_fact_TSV[i_s]; i_f++){

	/////////////////////////////////////////////////////////////
	//  state re-normalized distibutions for each pTratio-bin  //
	//  state re-normalized distibutions for each abseta-bin   //
	/////////////////////////////////////////////////////////////
	
	for (int i_b1 = 0; i_b1 < ygeneric->oset.dddat[i_ddd].distribution_1.n_bins; i_b1++){
	  stringstream filename_distribution_plot_ss;
	  filename_distribution_plot_ss << identifier << ".split." << ygeneric->oset.dddat[i_ddd].name << "_" << ygeneric->oset.dddat[i_ddd].distribution_1.bin_edge[i_b1] << "-" << ygeneric->oset.dddat[i_ddd].distribution_1.bin_edge[i_b1 + 1] << ".." << resultdirectory << ".dat";
	  string filename_distribution_plot = filename_distribution_plot_ss.str();
	  string filepath_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + subdirectory + "/" + filename_distribution_plot;
	  logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot              = " << filepath_distribution_plot << endl;
	  ofstream out_plotfile;
	  out_plotfile.open(filepath_distribution_plot.c_str(), ofstream::out | ofstream::trunc);  
	  double temp_result = 0;
	  double temp_deviation = 0;
	  for (int i_b2 = 0; i_b2 < ygeneric->oset.dddat[i_ddd].distribution_2.n_bins; i_b2++){
	    int i_b = i_b1 * ygeneric->oset.dddat[i_ddd].distribution_2.n_bins + i_b2;

	    if (identifier == "norm.norm"){
	      temp_result = this_distribution_result_TSV[i_b][i_s][i_r][i_f];
	      temp_deviation = this_distribution_deviation_TSV[i_b][i_s][i_r][i_f];
	    }
	    else if (identifier == "norm.plot"){
	      temp_result = this_distribution_result_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2];
	      temp_deviation = this_distribution_deviation_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2];
	    }
	    else if (identifier == "plot.norm"){
	      temp_result = this_distribution_result_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	      temp_deviation = this_distribution_deviation_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	    }
	    else if (identifier == "plot.plot"){
	      temp_result = this_distribution_result_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	      temp_deviation = this_distribution_deviation_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	    }
	    out_plotfile << setw(10) << ygeneric->oset.dddat[i_ddd].distribution_2.bin_edge[i_b2]
			 << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_result
			 << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_deviation 
			 << endl;
	  }
	  out_plotfile << setw(10) << ygeneric->oset.dddat[i_ddd].distribution_2.bin_edge[ygeneric->oset.dddat[i_ddd].distribution_2.n_bins]
			 << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_result
			 << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_deviation 
			 << endl; 
	  out_plotfile.close();
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_order::output_dddistribution_ge_lt_in_first_sdd_TSV(string & identifier, vector<vector<vector<vector<double> > > > & this_distribution_result_TSV, vector<vector<vector<vector<double> > > > & this_distribution_deviation_TSV, int i_d, string & subdirectory){
  Logger logger("summary_order::output_dddistribution_ge_lt_in_first_sdd_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int i_ddd = i_d - ygeneric->oset.dat.size();

  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    if (!ygeneric->oset.switch_distribution_TSV[i_s]){continue;}
    for (int i_r = 0; i_r < ygeneric->oset.n_scale_ren_TSV[i_s]; i_r++){
      for (int i_f = 0; i_f < ygeneric->oset.n_scale_fact_TSV[i_s]; i_f++){
	
	////////////////////////////////////////////////////////////////////
	//  state re-normalized distibutions for all pTratio <>-binnings  //
	////////////////////////////////////////////////////////////////////
	
	for (int i_b1 = 0; i_b1 < ygeneric->oset.dddat[i_ddd].distribution_1.n_bins; i_b1++){
	  stringstream filename_distribution_plot_lt_ss;
	  stringstream filename_distribution_plot_ge_ss;
	  filename_distribution_plot_lt_ss << identifier << ".two." << ygeneric->oset.dddat[i_ddd].name << "_lt" << ygeneric->oset.dddat[i_ddd].distribution_1.bin_edge[i_b1] << ".." << resultdirectory << ".dat";
	  filename_distribution_plot_ge_ss << identifier << ".two." << ygeneric->oset.dddat[i_ddd].name << "_ge" << ygeneric->oset.dddat[i_ddd].distribution_1.bin_edge[i_b1] << ".." << resultdirectory << ".dat";
	  string filename_distribution_plot_lt = filename_distribution_plot_lt_ss.str();
	  string filename_distribution_plot_ge = filename_distribution_plot_ge_ss.str();
	  string filepath_distribution_plot_lt = ygeneric->scalename_TSV[i_s][i_r][i_f] + subdirectory + "/" + filename_distribution_plot_lt;
	  string filepath_distribution_plot_ge = ygeneric->scalename_TSV[i_s][i_r][i_f] + subdirectory + "/" + filename_distribution_plot_ge;

	  logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot_lt              = " << filepath_distribution_plot_lt << endl;
	  logger << LOG_DEBUG_VERBOSE << "filepath_distribution_plot_ge              = " << filepath_distribution_plot_ge << endl;
	  
	  ////////////////////////////////////////////////////////////////////
	  //  state re-normalized distibutions for all pTratio <>-binnings  //
	  ////////////////////////////////////////////////////////////////////
	  
	  //  for (int i_b1 = 1; i_b1 < ygeneric->oset.dddat[i_ddd].distribution_1.n_bins; i_b1++){
	  ofstream out_plotfile_lt;
	  ofstream out_plotfile_ge;
	  out_plotfile_lt.open(filepath_distribution_plot_lt.c_str(), ofstream::out | ofstream::trunc);  
	  out_plotfile_ge.open(filepath_distribution_plot_ge.c_str(), ofstream::out | ofstream::trunc);  
	  for (int i_b2 = 0; i_b2 < ygeneric->oset.dddat[i_ddd].distribution_2.n_bins; i_b2++){
	    double temp_result = 0;
	    double temp_deviation = 0;
	    for (int j_b1 = 0; j_b1 < i_b1; j_b1++){
	      //      cout << "lt: ygeneric->oset.dddat[" << i_ddd << "].distribution_1.bin_edge[" << j_b1 << "] = " << ygeneric->oset.dddat[i_ddd].distribution_1.bin_edge[j_b1] << "   ygeneric->oset.dddat[" << i_ddd << "].distribution_1.bin_width[" << j_b1 << "] = " << ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[j_b1] << endl;
	      int i_b = j_b1 * ygeneric->oset.dddat[i_ddd].distribution_2.n_bins + i_b2;
	      if (identifier == "norm" || identifier == "norm.norm"){
		temp_result += this_distribution_result_TSV[i_b][i_s][i_r][i_f];
		temp_deviation += pow(this_distribution_deviation_TSV[i_b][i_s][i_r][i_f], 2);
	      }
	      else if (identifier == "plot" || identifier == "norm.plot"){
		temp_result += this_distribution_result_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2];
		temp_deviation += pow(this_distribution_deviation_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2], 2);
	      }
	      else if (identifier == "plot.norm"){
		temp_result += this_distribution_result_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
		temp_deviation += this_distribution_deviation_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	      }
	      else if (identifier == "plot.plot"){
		temp_result = +this_distribution_result_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
		temp_deviation = +this_distribution_deviation_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	      }
	    }
	    temp_deviation = sqrt(temp_deviation);
	    out_plotfile_lt << setw(10) << ygeneric->oset.dddat[i_ddd].distribution_2.bin_edge[i_b2] 
			    << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_result 
			    << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_deviation << endl;
	    if (i_b2 == ygeneric->oset.dddat[i_ddd].distribution_2.n_bins - 1){
	      out_plotfile_lt << setw(10) << ygeneric->oset.dddat[i_ddd].distribution_2.bin_edge[i_b2 + 1] 
			      << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_result 
			      << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_deviation << endl;}
	    
	    temp_result = 0;
	    temp_deviation = 0;
	    for (int j_b1 = i_b1; j_b1 < ygeneric->oset.dddat[i_ddd].distribution_1.n_bins; j_b1++){
	      //      cout << "ge: ygeneric->oset.dddat[" << i_ddd << "].distribution_1.bin_edge[" << j_b1 << "] = " << ygeneric->oset.dddat[i_ddd].distribution_1.bin_edge[j_b1] << "   ygeneric->oset.dddat[" << i_ddd << "].distribution_1.bin_width[" << j_b1 << "] = " << ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[j_b1] << endl;
	      int i_b = j_b1 * ygeneric->oset.dddat[i_ddd].distribution_2.n_bins + i_b2;
	      if (identifier == "norm" || identifier == "norm.norm"){
		temp_result += this_distribution_result_TSV[i_b][i_s][i_r][i_f];// / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2];
		temp_deviation += pow(this_distribution_deviation_TSV[i_b][i_s][i_r][i_f], 2);// / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2], 2);
	      }
	      else if (identifier == "plot" || identifier == "norm.plot"){
		temp_result += this_distribution_result_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2];
		temp_deviation += pow(this_distribution_deviation_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2], 2);
	      }
	      else if (identifier == "plot.norm"){
		temp_result += this_distribution_result_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
		temp_deviation += this_distribution_deviation_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	      }
	      else if (identifier == "plot.plot"){
		temp_result = +this_distribution_result_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
		temp_deviation = +this_distribution_deviation_TSV[i_b][i_s][i_r][i_f] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	      }
	    }
	    temp_deviation = sqrt(temp_deviation);
	    out_plotfile_ge << setw(10) << ygeneric->oset.dddat[i_ddd].distribution_2.bin_edge[i_b2] 
			    << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_result 
			    << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_deviation << endl;
	    if (i_b2 == ygeneric->oset.dddat[i_ddd].distribution_2.n_bins - 1){
	      out_plotfile_ge << setw(10) << ygeneric->oset.dddat[i_ddd].distribution_2.bin_edge[i_b2 + 1] 
			      << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_result 
			      << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_deviation << endl;}
	  }
	  out_plotfile_lt.close();
	  out_plotfile_ge.close();
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_order::output_distribution_TSV(){
  Logger logger("output_distribution_TSV");
  logger << LOG_DEBUG << "called" << endl;

  string subdirectory = "";

  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < ygeneric->oset.extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      string identifier = "";
      identifier = "norm";
      output_sddistribution_TSV(identifier, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], i_d, subdirectory);
      identifier = "plot";
      output_sddistribution_TSV(identifier, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], i_d, subdirectory);
    }
  }

  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < ygeneric->oset.extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      if (i_d < ygeneric->oset.dat.size()){continue;}
      string identifier = "";
      identifier = "norm";
      output_dddistribution_reconstruct_sdd_TSV(identifier, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], i_d, subdirectory);
      identifier = "plot";
      output_dddistribution_reconstruct_sdd_TSV(identifier, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], i_d, subdirectory);
    }
  }

  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < ygeneric->oset.extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      if (i_d < ygeneric->oset.dat.size()){continue;}
      string identifier = "";
      identifier = "norm.norm";
      output_dddistribution_split_in_first_sdd_TSV(identifier, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], i_d, subdirectory);
      identifier = "norm.plot";
      output_dddistribution_split_in_first_sdd_TSV(identifier, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], i_d, subdirectory);
      identifier = "plot.norm";
      output_dddistribution_split_in_first_sdd_TSV(identifier, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], i_d, subdirectory);
      identifier = "plot.plot";
      output_dddistribution_split_in_first_sdd_TSV(identifier, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], i_d, subdirectory);
    }
  }

  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < ygeneric->oset.extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      if (i_d < ygeneric->oset.dat.size()){continue;}
      string identifier = "";
      identifier = "norm";
      output_dddistribution_ge_lt_in_first_sdd_TSV(identifier, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], i_d, subdirectory);
      identifier = "plot";
      output_dddistribution_ge_lt_in_first_sdd_TSV(identifier, distribution_result_TSV[i_g][i_d], distribution_deviation_TSV[i_g][i_d], i_d, subdirectory);
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_order::output_distribution_qTcut_TSV(){
  Logger logger("summary_order::output_distribution_qTcut_TSV");
  logger << LOG_DEBUG << "called" << endl;

  if (!active_qTcut){return;}

  for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
    stringstream qTcut_ss;
    qTcut_ss << "qTcut-" << ygeneric->oset.value_qTcut_distribution[x_q];
    string subdirectory = "/" + qTcut_ss.str();
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      if (!ygeneric->oset.switch_distribution_TSV[i_s]){continue;}
      for (int i_r = 0; i_r < ygeneric->oset.n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < ygeneric->oset.n_scale_fact_TSV[i_s]; i_f++){
	  string directory_distribution_plot = ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + qTcut_ss.str();
	  system_execute(logger, "mkdir " + directory_distribution_plot);
	}
      }
    }

    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      for (int i_d = 0; i_d < ygeneric->oset.extended_distribution.size(); i_d++){
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	string identifier = "";
	identifier = "norm";
	output_sddistribution_TSV(identifier, distribution_result_qTcut_TSV[x_q][i_g][i_d], distribution_deviation_qTcut_TSV[x_q][i_g][i_d], i_d, subdirectory);
	identifier = "plot";
	output_sddistribution_TSV(identifier, distribution_result_qTcut_TSV[x_q][i_g][i_d], distribution_deviation_qTcut_TSV[x_q][i_g][i_d], i_d, subdirectory);
      }
    }
    
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      for (int i_d = 0; i_d < ygeneric->oset.extended_distribution.size(); i_d++){
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	if (i_d < ygeneric->oset.dat.size()){continue;}
	string identifier = "";
	identifier = "norm";
	output_dddistribution_reconstruct_sdd_TSV(identifier, distribution_result_qTcut_TSV[x_q][i_g][i_d], distribution_deviation_qTcut_TSV[x_q][i_g][i_d], i_d, subdirectory);
	identifier = "plot";
	output_dddistribution_reconstruct_sdd_TSV(identifier, distribution_result_qTcut_TSV[x_q][i_g][i_d], distribution_deviation_qTcut_TSV[x_q][i_g][i_d], i_d, subdirectory);
      }
    }
    
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      for (int i_d = 0; i_d < ygeneric->oset.extended_distribution.size(); i_d++){
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	if (i_d < ygeneric->oset.dat.size()){continue;}
	string identifier = "";
	identifier = "norm.norm";
	output_dddistribution_split_in_first_sdd_TSV(identifier, distribution_result_qTcut_TSV[x_q][i_g][i_d], distribution_deviation_qTcut_TSV[x_q][i_g][i_d], i_d, subdirectory);
	identifier = "norm.plot";
	output_dddistribution_split_in_first_sdd_TSV(identifier, distribution_result_qTcut_TSV[x_q][i_g][i_d], distribution_deviation_qTcut_TSV[x_q][i_g][i_d], i_d, subdirectory);
	identifier = "plot.norm";
	output_dddistribution_split_in_first_sdd_TSV(identifier, distribution_result_qTcut_TSV[x_q][i_g][i_d], distribution_deviation_qTcut_TSV[x_q][i_g][i_d], i_d, subdirectory);
	identifier = "plot.plot";
	output_dddistribution_split_in_first_sdd_TSV(identifier, distribution_result_qTcut_TSV[x_q][i_g][i_d], distribution_deviation_qTcut_TSV[x_q][i_g][i_d], i_d, subdirectory);
      }
    }
    
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      for (int i_d = 0; i_d < ygeneric->oset.extended_distribution.size(); i_d++){
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	if (i_d < ygeneric->oset.dat.size()){continue;}
	string identifier = "";
	identifier = "norm";
	output_dddistribution_ge_lt_in_first_sdd_TSV(identifier, distribution_result_qTcut_TSV[x_q][i_g][i_d], distribution_deviation_qTcut_TSV[x_q][i_g][i_d], i_d, subdirectory);
	identifier = "plot";
	output_dddistribution_ge_lt_in_first_sdd_TSV(identifier, distribution_result_qTcut_TSV[x_q][i_g][i_d], distribution_deviation_qTcut_TSV[x_q][i_g][i_d], i_d, subdirectory);
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}





void summary_order::output_sddistribution_CV(string & identifier, vector<vector<double> > & this_distribution_result_CV, vector<vector<double> > & this_distribution_deviation_CV, int i_d, string & subdirectory){
  Logger logger("summary_order::output_sddistribution_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  /////////////////////////////////////////
  //  singly-differential distributions  //
  /////////////////////////////////////////
  
  int plotmode = 0;

  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
    string filename_distribution_plot = identifier + "." + osi->extended_distribution[i_d].xdistribution_name + ".." + resultdirectory + ".dat";
    string filename_full = ygeneric->final_resultdirectory + "/CV/"  + osi->directory_name_scale_CV[i_s] + subdirectory + "/" + filename_distribution_plot;
    logger << LOG_DEBUG_VERBOSE << "filename_full              = " << filename_full << endl;
    ofstream out_plotfile;
    out_plotfile.open(filename_full.c_str(), ofstream::out | ofstream::trunc);  
    double temp_result = 0;
    double temp_deviation = 0;
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
      if (identifier == "norm"){
	temp_result = this_distribution_result_CV[i_b][i_s];
	temp_deviation = this_distribution_deviation_CV[i_b][i_s];
      }
      else if (identifier == "plot"){
	temp_result = this_distribution_result_CV[i_b][i_s] / osi->extended_distribution[i_d].bin_width[i_b];
	temp_deviation = this_distribution_deviation_CV[i_b][i_s] / osi->extended_distribution[i_d].bin_width[i_b];
      }
      if (plotmode == 1){if (temp_result < 0.){temp_result = 0.;}}

      out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[i_b] 
		   << setw(16) << setprecision(8) << osi->unit_factor_distribution * temp_result
		   << setw(16) << setprecision(8) << osi->unit_factor_distribution * temp_deviation << endl;
    }
    out_plotfile << setw(10) << osi->extended_distribution[i_d].bin_edge[osi->extended_distribution[i_d].n_bins] 
		 << setw(16) << setprecision(8) << osi->unit_factor_distribution * temp_result
		 << setw(16) << setprecision(8) << osi->unit_factor_distribution * temp_deviation << endl; 
    out_plotfile.close();
    logger << LOG_DEBUG_VERBOSE << "singly-differential output finished" << endl;
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_order::output_distribution_CV(){
  Logger logger("output_distribution_CV");
  logger << LOG_DEBUG << "called" << endl;

  string subdirectory = "";

  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      string identifier = "";
      identifier = "norm";
      output_sddistribution_CV(identifier, distribution_result_CV[i_g][i_d], distribution_deviation_CV[i_g][i_d], i_d, subdirectory);
      identifier = "plot";
      output_sddistribution_CV(identifier, distribution_result_CV[i_g][i_d], distribution_deviation_CV[i_g][i_d], i_d, subdirectory);
    }
  }

  /*
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      if (i_d < osi->dat.size()){continue;}
      string identifier = "";
      identifier = "norm";
      output_dddistribution_reconstruct_sdd_CV(identifier, distribution_result_CV[i_g][i_d], distribution_deviation_CV[i_g][i_d], i_d, subdirectory);
      identifier = "plot";
      output_dddistribution_reconstruct_sdd_CV(identifier, distribution_result_CV[i_g][i_d], distribution_deviation_CV[i_g][i_d], i_d, subdirectory);
    }
  }
  */
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      if (i_d < osi->dat.size()){continue;}
      string identifier = "";
      identifier = "norm.norm";
      output_dddistribution_split_in_first_sdd_CV(identifier, distribution_result_CV[i_g][i_d], distribution_deviation_CV[i_g][i_d], i_d, subdirectory);
      identifier = "norm.plot";
      output_dddistribution_split_in_first_sdd_CV(identifier, distribution_result_CV[i_g][i_d], distribution_deviation_CV[i_g][i_d], i_d, subdirectory);
      identifier = "plot.norm";
      output_dddistribution_split_in_first_sdd_CV(identifier, distribution_result_CV[i_g][i_d], distribution_deviation_CV[i_g][i_d], i_d, subdirectory);
      identifier = "plot.plot";
      output_dddistribution_split_in_first_sdd_CV(identifier, distribution_result_CV[i_g][i_d], distribution_deviation_CV[i_g][i_d], i_d, subdirectory);
    }
  }
  /*
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      if (i_d < osi->dat.size()){continue;}
      string identifier = "";
      identifier = "norm";
      output_dddistribution_ge_lt_in_first_sdd_CV(identifier, distribution_result_CV[i_g][i_d], distribution_deviation_CV[i_g][i_d], i_d, subdirectory);
      identifier = "plot";
      output_dddistribution_ge_lt_in_first_sdd_CV(identifier, distribution_result_CV[i_g][i_d], distribution_deviation_CV[i_g][i_d], i_d, subdirectory);
    }
  }
  */
    
  logger << LOG_DEBUG << "finished" << endl;
}


void summary_order::output_dddistribution_split_in_first_sdd_CV(string & identifier, vector<vector<double> > & this_distribution_result_CV, vector<vector<double> > & this_distribution_deviation_CV, int i_d, string & subdirectory){
  Logger logger("summary_order::output_dddistribution_split_in_first_sdd_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int i_ddd = i_d - ygeneric->oset.dat.size();

  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){

    /////////////////////////////////////////////////////////////
    //  state re-normalized distibutions for each pTratio-bin  //
    //  state re-normalized distibutions for each abseta-bin   //
    /////////////////////////////////////////////////////////////
    
    for (int i_b1 = 0; i_b1 < ygeneric->oset.dddat[i_ddd].distribution_1.n_bins; i_b1++){
      //    string filename_distribution_plot = identifier + "." + osi->extended_distribution[i_d].xdistribution_name + ".." + resultdirectory + ".dat";
      stringstream filename_distribution_plot_ss;
      filename_distribution_plot_ss << identifier << ".split." << ygeneric->oset.dddat[i_ddd].name << "_" << ygeneric->oset.dddat[i_ddd].distribution_1.bin_edge[i_b1] << "-" << ygeneric->oset.dddat[i_ddd].distribution_1.bin_edge[i_b1 + 1] << ".." << resultdirectory << ".dat";
      string filename_distribution_plot = filename_distribution_plot_ss.str();
      string filename_full = ygeneric->final_resultdirectory + "/CV/"  + osi->directory_name_scale_CV[i_s] + subdirectory + "/" + filename_distribution_plot;
      ofstream out_plotfile;
      out_plotfile.open(filename_full.c_str(), ofstream::out | ofstream::trunc);  
      
      
      double temp_result = 0;
      double temp_deviation = 0;
      for (int i_b2 = 0; i_b2 < ygeneric->oset.dddat[i_ddd].distribution_2.n_bins; i_b2++){
	int i_b = i_b1 * ygeneric->oset.dddat[i_ddd].distribution_2.n_bins + i_b2;
	/*
	if (identifier == "norm.norm"){
	  temp_result = this_distribution_result_CV[i_b][i_s] * ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2] * ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	  temp_deviation = this_distribution_deviation_CV[i_b][i_s] * ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2] * ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	}
	else if (identifier == "norm.plot"){
	  temp_result = this_distribution_result_CV[i_b][i_s] * ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	  temp_deviation = this_distribution_deviation_CV[i_b][i_s] * ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	}
	else if (identifier == "plot.norm"){
	  temp_result = this_distribution_result_CV[i_b][i_s] * ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2];
	  temp_deviation = this_distribution_deviation_CV[i_b][i_s] * ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2];
	}
	else if (identifier == "plot.plot"){
	  temp_result = this_distribution_result_CV[i_b][i_s];
	  temp_deviation = this_distribution_deviation_CV[i_b][i_s];
	}
	*/
	
	if (identifier == "norm.norm"){
	  temp_result = this_distribution_result_CV[i_b][i_s];
	  temp_deviation = this_distribution_deviation_CV[i_b][i_s];
	}
	else if (identifier == "norm.plot"){
	  temp_result = this_distribution_result_CV[i_b][i_s] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2];
	  temp_deviation = this_distribution_deviation_CV[i_b][i_s] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2];
	}
	else if (identifier == "plot.norm"){
	  temp_result = this_distribution_result_CV[i_b][i_s] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	  temp_deviation = this_distribution_deviation_CV[i_b][i_s] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	}
	else if (identifier == "plot.plot"){
	  temp_result = this_distribution_result_CV[i_b][i_s] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	  temp_deviation = this_distribution_deviation_CV[i_b][i_s] / ygeneric->oset.dddat[i_ddd].distribution_2.bin_width[i_b2] / ygeneric->oset.dddat[i_ddd].distribution_1.bin_width[i_b1];
	}
	
	out_plotfile << setw(10) << ygeneric->oset.dddat[i_ddd].distribution_2.bin_edge[i_b2]
		     << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_result
		     << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_deviation 
		     << endl;
      }
      out_plotfile << setw(10) << ygeneric->oset.dddat[i_ddd].distribution_2.bin_edge[ygeneric->oset.dddat[i_ddd].distribution_2.n_bins]
		   << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_result
		   << setw(16) << setprecision(8) << ygeneric->oset.unit_factor_distribution * temp_deviation 
		   << endl; 
      out_plotfile.close();
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



