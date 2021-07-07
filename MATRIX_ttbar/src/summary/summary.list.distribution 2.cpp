#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_list::collect_contribution_distribution_CV(){
  Logger logger("summary_list::collect_contribution_distribution_CV");
  logger << LOG_DEBUG << "called" << endl;

  logger << LOG_DEBUG << "resultdirectory = " << resultdirectory << endl;

  string filename;
  vector<string> readin;

  logger << LOG_DEBUG << "osi->extended_distribution.size() = " << osi->extended_distribution.size() << endl;

  xcontribution[0].distribution_result_CV.resize(ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->extended_distribution.size()));
  xcontribution[0].distribution_deviation_CV.resize(ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->extended_distribution.size()));
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      xcontribution[0].distribution_result_CV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins, vector<double> (osi->n_scales_CV, 0.));
      xcontribution[0].distribution_deviation_CV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins, vector<double> (osi->n_scales_CV, 0.));
    }
  }

  //////////////////////////////////////////
  //  creation of output CV directories  //
  //////////////////////////////////////////
  // system_execute(logger, "mkdir " + ygeneric->final_resultdirectory);
  //  system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV");
  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
    system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s]);
    system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s] + "/" + resultdirectory);
    for (int i_c = 1; i_c < xcontribution.size(); i_c++){
      system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s] + "/" + resultdirectory + "/" + xcontribution[i_c].infix_contribution);
      //      system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s] + "/" + resultdirectory + "/" + xcontribution[i_c].infix_order_contribution);
    }
  }

  //  double xosi->fakeasymfactor = 0;
  /*
  vector<vector<double> > bin_edge(osi->extended_distribution.size());
  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    bin_edge[i_d] = osi->extended_distribution[i_d].bin_edge;
    /*
    bin_edge[i_d].resize(osi->extended_distribution[i_d].n_bins + 1);
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins + 1; i_b++){bin_edge[i_d][i_b] = osi->extended_distribution[i_d].start + i_b * osi->extended_distribution[i_d].step;}
  *//*
  }
 */
  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    logger << LOG_DEBUG_VERBOSE << "xcontribution[" << i_c << "].subprocess.size() = " << xcontribution[i_c].subprocess.size() << endl;
    logger << LOG_DEBUG_VERBOSE << "xcontribution[" << i_c << "].xsubprocess.size() = " << xcontribution[i_c].xsubprocess.size() << endl;
    //    xcontribution[i_c].readin_distribution_contribution_CV(bin_edge);
    xcontribution[i_c].readin_distribution_contribution_CV();
  }
    
  // n_contribution -> infix... ???
  //  logger << LOG_DEBUG_VERBOSE << "ygeneric->subgroup.size() = " << ygeneric->subgroup.size() << endl;
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    //    logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      //      logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	//	logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	// ???	double dev = 0.;
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  double dev = 0.;
	  //	  logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
	  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
	    //	    logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << endl;
	    if (xcontribution[i_c].distribution_result_CV[i_g][i_d][i_b][i_s] != 0){
	      logger << LOG_DEBUG_VERBOSE << "distribution: result_CV[" << setw(2) << i_c << "][" << setw(2) << i_g << "][" << setw(2) << i_d << "][" << setw(2) << i_b << "][" << setw(2) << i_s << "] = " << setw(23) << setprecision(15) << xcontribution[i_c].distribution_result_CV[i_g][i_d][i_b][i_s]<< " +- " << setw(23) << setprecision(15) << xcontribution[i_c].distribution_deviation_CV[i_g][i_d][i_b][i_s] << endl;
	    }
	    xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s] += xcontribution[i_c].distribution_result_CV[i_g][i_d][i_b][i_s];
	    dev += pow(xcontribution[i_c].distribution_deviation_CV[i_g][i_d][i_b][i_s], 2);
	  }
	  xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s] = sqrt(dev);
	  //	  logger << LOG_DEBUG_VERBOSE << "result_CV[0][" << setw(2) << i_g << "][" << setw(2) << i_d << "][" << setw(2) << i_b << "][" << setw(2) << i_s << "] = " << setw(23) << setprecision(15) << result_CV[0][i_g][i_d][i_b][i_s]<< " +- " << setw(23) << setprecision(15) << deviation_CV[0][i_g][i_d][i_b][i_s]<< endl;

	  //	  list_result_CV[i_g][i_d][i_b][i_s] = result_CV[0][i_g][i_d][i_b][i_s];
	  //	  list_deviation_CV[i_g][i_d][i_b][i_s] = deviation_CV[0][i_g][i_d][i_b][i_s];
	  //	  logger << LOG_DEBUG_VERBOSE << "list_result_CV[" << setw(2) << i_g << "][" << setw(2) << i_d << "][" << setw(2) << i_b << "][" << setw(2) << i_s << "] = " << setw(23) << setprecision(15) << list_result_CV[i_g][i_d][i_b][i_s]<< " +- " << setw(23) << setprecision(15) << list_deviation_CV[i_g][i_d][i_b][i_s]<< endl;
	  /*
	    list_result_CV[i_m][i_g][i_q][i_s] = result_CV[i_m][0][i_g][i_q][i_s];
	    list_deviation_CV[i_m][i_g][i_q][i_s] = deviation_CV[i_m][0][i_g][i_q][i_s];
	  */
	}
	//	logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
      }
      //      logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
    }
    //    logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
  }

  for (int i_g = 0; i_g < 1; i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	if ((osi->extended_distribution[i_d].xdistribution_name).substr(0, 4) == "asym"){
	  /*
	  double result_fw, result_bw, deviation_fw, deviation_bw;
	  ofstream out_asymfile;
	  string filename_asym = ygeneric->final_resultdirectory + "/" + resultdirectory + "/" + scalename_CV[i_s] + "/" + osi->extended_distribution[i_d].xdistribution_name + "." + infix_order_contribution[0] + ".dat";
	  logger << LOG_DEBUG_VERBOSE << "filename_asym = " << filename_asym << endl;
	  out_asymfile.open(filename_asym.c_str(), ofstream::out | ofstream::trunc);  
	  if (osi->extended_distribution[i_d].xdistribution_type == "asym1" || osi->extended_distribution[i_d].xdistribution_type == "asym2"){
	    ofstream out_plotfile;
	    string filename_plot = ygeneric->final_resultdirectory + "/" + resultdirectory + "/" + scalename_CV[i_s] + "/plot." + osi->extended_distribution[i_d].xdistribution_name + "." + infix_order_contribution[0] + ".dat";
	    logger << LOG_DEBUG_VERBOSE << "filename_plot = " << filename_plot << endl;
	    out_plotfile.open(filename_plot.c_str(), ofstream::out | ofstream::trunc);  
	    // too slow // checksum << setw(15) << "interval" << " " << setw(15) << "backward cs." << setw(15) << "forward cs." << setw(15) << "cross section" << endl;
	    for (int i_b = 0; i_b < n_bin_distribution[i_d]; i_b++){
	      if (osi->extended_distribution[i_d].xdistribution_type == "asym1"){
		result_fw = result_CV[0][i_g][i_d][n_bin_distribution[i_d] + i_b][i_s] * osi->fakeasymfactor[i_d];
		result_bw = result_CV[0][i_g][i_d][i_b][i_s] * osi->fakeasymfactor[i_d];
		deviation_fw = deviation_CV[0][i_g][i_d][n_bin_distribution[i_d] + i_b][i_s] * osi->fakeasymfactor[i_d];
		deviation_bw = deviation_CV[0][i_g][i_d][i_b][i_s] * osi->fakeasymfactor[i_d];
	      }
	      else if (osi->extended_distribution[i_d].xdistribution_type == "asym2"){
		result_fw = result_CV[0][i_g][i_d][2 * i_b + 1][i_s] * osi->fakeasymfactor[i_d];
		result_bw = result_CV[0][i_g][i_d][2 * i_b + 0][i_s] * osi->fakeasymfactor[i_d];
		deviation_fw = deviation_CV[0][i_g][i_d][2 * i_b + 1][i_s] * osi->fakeasymfactor[i_d];
		deviation_bw = deviation_CV[0][i_g][i_d][2 * i_b + 0][i_s] * osi->fakeasymfactor[i_d];
	      }
	      out_asymfile << setw(25) << "backw. cross section: " << "   [" << setw(5) << setprecision(4) << bin_edge[i_d][i_b] << ";" << setw(5) << setprecision(4) << bin_edge[i_d][i_b + 1] << "]   " << setw(15) << setprecision(8) << result_bw * osi->extended_distribution[i_d].step << setw(15) << setprecision(8) << deviation_bw * osi->extended_distribution[i_d].step << endl;
	      out_asymfile << setw(25) << "forw. cross section: " << "   [" << setw(5) << setprecision(4) << bin_edge[i_d][i_b] << ";" << setw(5) << setprecision(4) << bin_edge[i_d][i_b + 1] << "]   " << setw(15) << setprecision(8) << result_fw * osi->extended_distribution[i_d].step << setw(15) << setprecision(8) << deviation_fw * osi->extended_distribution[i_d].step << endl;
	      out_asymfile << setw(25) << "forw.-backw. asymmetry: " << "   [" << setw(5) << setprecision(4) << bin_edge[i_d][i_b] << ";" << setw(5) << setprecision(4) << bin_edge[i_d][i_b + 1] << "]   " << setw(15) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(15) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	      out_plotfile << setw(5) << setprecision(4) << bin_edge[i_d][i_b] << setw(15) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(15) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	      double sum_result_fw = 0., sum_result_bw = 0.;
	      for (int c = i_b; c < n_bin_distribution[i_d]; c++){
		if (osi->extended_distribution[i_d].xdistribution_type == "asym1"){
		  sum_result_fw += result_CV[0][i_g][i_d][n_bin_distribution[i_d] + c][i_s] * osi->fakeasymfactor[i_d];
		  sum_result_bw += result_CV[0][i_g][i_d][c][i_s] * osi->fakeasymfactor[i_d];
		}
		else if (osi->extended_distribution[i_d].xdistribution_type == "asym2"){
		  sum_result_fw += result_CV[0][i_g][i_d][2 * c + 1][i_s] * osi->fakeasymfactor[i_d];
		  sum_result_bw += result_CV[0][i_g][i_d][2 * c + 0][i_s] * osi->fakeasymfactor[i_d];
		}
	      }
	    }
	    out_plotfile << setw(5) << setprecision(4) << osi->extended_distribution[i_d].end << setw(15) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(15) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	    out_asymfile.close();
	    out_plotfile.close();
	  }
	  else {
	    result_fw = result_CV[0][i_g][i_d][1][i_s] * osi->fakeasymfactor[i_d];
	    result_bw = result_CV[0][i_g][i_d][0][i_s] * osi->fakeasymfactor[i_d];
	    deviation_fw = deviation_CV[0][i_g][i_d][1][i_s] * osi->fakeasymfactor[i_d];
	    deviation_bw = deviation_CV[0][i_g][i_d][0][i_s] * osi->fakeasymfactor[i_d];
	    out_asymfile << setw(25) << "backw. cross section: " << setw(15) << setprecision(8) << result_bw << setw(15) << setprecision(8) << deviation_bw << endl;
	    out_asymfile << setw(25) << "forw. cross section: " << setw(15) << setprecision(8) << result_fw << setw(15) << setprecision(8) << deviation_fw << endl;
	    out_asymfile << setw(25) << "forw.-backw. asymmetry: " << setw(15) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(15) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	    out_asymfile.close();
	  }
	  */
	}
	else {
	  ofstream out_plotfile;
	  string filename_plot = ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s] + "/" + resultdirectory + "/plot." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[0].infix_order_contribution + ".dat";
	  logger << LOG_DEBUG_VERBOSE << "filename_plot = " << filename_plot << endl;
	  out_plotfile.open(filename_plot.c_str(), ofstream::out | ofstream::trunc);  
	  for (int i_b = 0; i_b < xcontribution[0].distribution_result_CV[i_g][i_d].size(); i_b++){
	    out_plotfile << left << setw(10) << noshowpoint << osi->extended_distribution[i_d].bin_edge[i_b] << setw(15) << setprecision(8) << osi->unit_factor_result * xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s] / osi->extended_distribution[i_d].bin_width[i_b] << setw(15) << setprecision(8) << osi->unit_factor_result * xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s] / osi->extended_distribution[i_d].bin_width[i_b] << endl;
	    // norm ->	    out_plotfile << left << setw(10) << noshowpoint << osi->extended_distribution[i_d].bin_edge[i_b] << setw(15) << setprecision(8) << osi->unit_factor_result * xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s] << setw(15) << setprecision(8) << osi->unit_factor_result * xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s] << endl;
	    double sum_result = 0.;
	    for (int j_b = i_b; j_b < osi->extended_distribution[i_d].n_bins; j_b++){sum_result += xcontribution[0].distribution_result_CV[i_g][i_d][j_b][i_s];}
	    // norm ->	    for (int j_b = i_b; j_b < osi->extended_distribution[i_d].n_bins; j_b++){sum_result += xcontribution[0].distribution_result_CV[i_g][i_d][j_b][i_s] / osi->extended_distribution[i_d].bin_width[i_b];}
	  }
	  out_plotfile << left << setw(10) << noshowpoint << osi->extended_distribution[i_d].end << setw(15) << setprecision(8) << showpoint << osi->unit_factor_result * xcontribution[0].distribution_result_CV[i_g][i_d][xcontribution[0].distribution_result_CV[i_g][i_d].size() - 1][i_s] / osi->extended_distribution[i_d].bin_width[osi->extended_distribution[i_d].n_bins - 1] << setw(15) << setprecision(8) << osi->unit_factor_result * xcontribution[0].distribution_deviation_CV[i_g][i_d][xcontribution[0].distribution_deviation_CV[i_g][i_d].size() - 1][i_s] / osi->extended_distribution[i_d].bin_width[osi->extended_distribution[i_d].n_bins - 1] << endl; // endpoint
	  // norm ->	    out_plotfile << left << setw(10) << noshowpoint << osi->extended_distribution[i_d].end << setw(15) << setprecision(8) << showpoint << osi->unit_factor_result * xcontribution[0].distribution_result_CV[i_g][i_d][xcontribution[0].distribution_result_CV[i_g][i_d].size() - 1][i_s] << setw(15) << setprecision(8) << osi->unit_factor_result * xcontribution[0].distribution_deviation_CV[i_g][i_d][xcontribution[0].distribution_deviation_CV[i_g][i_d].size() - 1][i_s] << endl; // endpoint
	  out_plotfile.close();


	  /*
	  if (xcontribution.size() > 1){
	    ofstream out_overviewfile;
	    string filename_overview = ygeneric->final_resultdirectory + "/CV/" + osi->directory_name_scale_CV[i_s] + "/" + resultdirectory + "/overview." + osi->extended_distribution[i_d].xdistribution_name + "." + xcontribution[0].infix_order_contribution + ".dat";
	    logger << LOG_DEBUG_VERBOSE << "filename_overview = " << filename_overview << endl;
	    out_overviewfile.open(filename_overview.c_str(), ofstream::out | ofstream::trunc);  
	    for (int i_b = 0; i_b < xcontribution[0].distribution_result_CV[i_g][i_d].size(); i_b++){
	      for (int i_c = 0; i_c < xcontribution.size(); i_c++){
		string sign = "";
		string minussign = " ";
		if (xcontribution[i_c].distribution_result_CV[i_g][i_d][i_b][i_s] >= 0){sign = " "; minussign = "";}
		out_overviewfile << left << setw(15) << setprecision(8) << noshowpoint << osi->extended_distribution[i_d].bin_edge[i_b] << setw(20) << left << xcontribution[i_c].infix_order_contribution << sign << showpoint << setw(15) << setprecision(8) << osi->unit_factor_result * xcontribution[i_c].distribution_result_CV[i_g][i_d][i_b][i_s] << minussign << "+-  " << setw(15) << setprecision(8) << osi->unit_factor_result * xcontribution[i_c].distribution_deviation_CV[i_g][i_d][i_b][i_s] << endl;
		if (i_c == 0){out_overviewfile << endl;}
	      }
	      out_overviewfile << endl;
	    }
	    out_overviewfile.close();
	  }
	  */
	}
      }
    }
  }

  if (ygeneric->switch_output_overview > 1){output_distribution_overview_CV();}

  logger << LOG_DEBUG << "finished" << endl;
}


void summary_list::collect_contribution_distribution_TSV(){
  Logger logger("summary_list::collect_contribution_distribution_TSV");
  logger << LOG_DEBUG << "called" << endl;

  for (int i_c = 0; i_c < xcontribution.size(); i_c++){
    logger << LOG_INFO << "xcontribution[" << i_c << "].type_contribution = " << xcontribution[i_c].type_contribution << endl;


    
    ///////////////////////////////////////////////////////////////////
    //  resize xcontribution[i_c].distribution_result/deviation_TSV  //
    ///////////////////////////////////////////////////////////////////

    xcontribution[i_c].distribution_result_TSV.resize(ygeneric->subgroup.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->extended_distribution.size()));
    xcontribution[i_c].distribution_deviation_TSV.resize(ygeneric->subgroup.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->extended_distribution.size()));
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	xcontribution[i_c].distribution_result_TSV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
	xcontribution[i_c].distribution_deviation_TSV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	  xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b].resize(selection_n_qTcut);
	  xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b].resize(selection_n_qTcut);
	  for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
	    xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][x_q].resize(osi->n_extended_set_TSV);
	    xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][x_q].resize(osi->n_extended_set_TSV);
	    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	      if (!osi->switch_distribution_TSV[i_s]){continue;}
	      //	      if (osi->switch_distribution_TSV[i_s] != 0){
	      xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	      xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	      //logger << LOG_DEBUG_VERBOSE << "xcontribution[" << i_c << "].distribution_result_TSV[" << i_g << "][" << i_d << "][" << i_b << "][" << x_q << "].size() = " << xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][x_q].size() << endl;
	      //logger << LOG_DEBUG_VERBOSE << "xcontribution[" << i_c << "].distribution_deviation_TSV[" << i_g << "][" << i_d << "][" << i_b << "][" << x_q << "].size() = " << xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][x_q].size() << endl;
	      //    }
	    }
	  }
	}
      }
    }
  }
  logger << LOG_DEBUG_VERBOSE << "xcontribution[xxx].distribution_result_TSV   resized!" << endl;



  //////////////////////////////////////////
  //  declaration of subgroup_result_TSV  //
  //////////////////////////////////////////

  distribution_result_qTcut_TSV.resize(ygeneric->subgroup.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->extended_distribution.size()));
  distribution_deviation_qTcut_TSV.resize(ygeneric->subgroup.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->extended_distribution.size()));
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      distribution_result_qTcut_TSV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
      distribution_deviation_qTcut_TSV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	distribution_result_qTcut_TSV[i_g][i_d][i_b].resize(selection_n_qTcut);
	distribution_deviation_qTcut_TSV[i_g][i_d][i_b].resize(selection_n_qTcut);
	for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
	  distribution_result_qTcut_TSV[i_g][i_d][i_b][x_q].resize(osi->n_extended_set_TSV);
	  distribution_deviation_qTcut_TSV[i_g][i_d][i_b][x_q].resize(osi->n_extended_set_TSV);
	  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    if (!osi->switch_distribution_TSV[i_s]){continue;}
	    //	    if (osi->switch_distribution_TSV[i_s] != 0){
	    distribution_result_qTcut_TSV[i_g][i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    distribution_deviation_qTcut_TSV[i_g][i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    //	    }
	  }
	}
      }
    }
  }



  //  system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/" + resultdirectory);
  
  //////////////////////////////////////////
  //  creation of output TSV directories  //
  //////////////////////////////////////////

  vector<vector<vector<string> > > scalename_TSV(osi->n_extended_set_TSV); 
  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    if (!osi->switch_distribution_TSV[i_s]){continue;}
    scalename_TSV[i_s].resize(osi->n_scale_ren_TSV[i_s], vector<string> (osi->n_scale_fact_TSV[i_s]));
    //    if (osi->switch_distribution_TSV[i_s] != 0){
    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	stringstream scalepart;
	scalepart << "scale." << i_r << "." << i_f;
	scalename_TSV[i_s][i_r][i_f] = ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + scalepart.str();
	//	  scalename_TSV[i_s][i_r][i_f] = ygeneric->final_resultdirectory + "/" + resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + scalepart.str();
	system_execute(logger, "mkdir " + ygeneric->scalename_TSV[i_s][i_r][i_f]);
	system_execute(logger, "mkdir " + ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + resultdirectory);
	for (int i_c = 1; i_c < xcontribution.size(); i_c++){
	  system_execute(logger, "mkdir " + ygeneric->scalename_TSV[i_s][i_r][i_f] + "/" + resultdirectory + "/" + xcontribution[i_c].infix_contribution);
	  //	    logger << LOG_INFO << "mkdir " << ygeneric->scalename_TSV[i_s][i_r][i_f] << "/" << resultdirectory << "/" << xcontribution[i_c].infix_contribution << endl;
	}
      }
    }
  }
  //}
  
    //  double xosi->fakeasymfactor = 0;
    /*
  vector<vector<double> > bin_edge(osi->extended_distribution.size());
  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    bin_edge[i_d].resize(osi->extended_distribution[i_d].n_bins + 1);
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins + 1; i_b++){bin_edge[i_d][i_b] = osi->extended_distribution[i_d].start + i_b * osi->extended_distribution[i_d].step;}
  }
    */

  ////////////////////////////////////////////////////////////////////
  //  read-in xcontribution[i_c].distribution_result/deviation_TSV  //
  ////////////////////////////////////////////////////////////////////
 
  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    //    xcontribution[i_c].readin_distribution_contribution_TSV(bin_edge);
    xcontribution[i_c].readin_distribution_contribution_TSV();
  }

  ////////////////////////////////////////////////////////////////////
  //  calculate xcontribution[0].distribution_result/deviation_TSV  //
  ////////////////////////////////////////////////////////////////////

  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
	  logger << LOG_DEBUG_VERBOSE << "start   x_q = " << x_q << endl;
	  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    if (!osi->switch_distribution_TSV[i_s]){continue;}
	    //	    if (osi->switch_distribution_TSV[i_s] != 0){
	    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
		double dev = 0.;
		for (int i_c = 1; i_c < xcontribution.size(); i_c++){
		  int y_q = 0;
		  if (xcontribution[i_c].active_qTcut){y_q = x_q;}
		  if (xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] != 0){
		    logger << LOG_DEBUG_VERBOSE << "distribution: result_TSV[" << setw(2) << i_c << "][" << setw(2) << i_g << "][" << setw(2) << i_d << "][" << setw(2) << i_b << "][" << setw(2) << y_q << "][" << setw(2) << i_s << "][" << setw(2) << i_r << "][" << setw(2) << i_f << "] = " << setw(23) << setprecision(15) << xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f] << endl;
		  }
		  xcontribution[0].distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] += xcontribution[i_c].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f];
		  dev += pow(xcontribution[i_c].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f], 2);
		  logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << "   x_q = " << x_q << "   y_q = " << y_q << endl;
		}
		xcontribution[0].distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] = sqrt(dev);

		distribution_result_qTcut_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] = xcontribution[0].distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f];
		distribution_deviation_qTcut_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] = xcontribution[0].distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f];
	      }
	      //   }
	    }
	  }
	  logger << LOG_DEBUG_VERBOSE << "end     x_q = " << x_q << endl;
	}
      }
    }
  }


  ////////////////////////////////////////////////
  //  resize distribution_result/deviation_TSV  //
  ////////////////////////////////////////////////

  distribution_result_TSV.resize(ygeneric->subgroup.size(), vector<vector<vector<vector<vector<double> > > > > (osi->extended_distribution.size()));
  distribution_deviation_TSV.resize(ygeneric->subgroup.size(), vector<vector<vector<vector<vector<double> > > > > (osi->extended_distribution.size()));
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      distribution_result_TSV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
      distribution_deviation_TSV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	distribution_result_TSV[i_g][i_d][i_b].resize(osi->n_extended_set_TSV);
	distribution_deviation_TSV[i_g][i_d][i_b].resize(osi->n_extended_set_TSV);
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  if (!osi->switch_distribution_TSV[i_s]){continue;}
	  //	  if (osi->switch_distribution_TSV[i_s] != 0){
	  distribution_result_TSV[i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  distribution_deviation_TSV[i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  //	  }
	}
      }
    }
  }
  

  ////////////////////////////////
  //  extrapolation qTcut -> 0  //
  ////////////////////////////////

  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	if (active_qTcut && selection_n_qTcut > 1){
	  extrapolation_TSV(osi->value_qTcut_distribution, xcontribution[0].distribution_result_TSV[i_g][i_d][i_b], xcontribution[0].distribution_deviation_TSV[i_g][i_d][i_b], distribution_result_TSV[i_g][i_d][i_b], distribution_deviation_TSV[i_g][i_d][i_b], ygeneric->switch_extrapolation_distribution);
	}
	else {
	  distribution_result_TSV[i_g][i_d][i_b] = xcontribution[0].distribution_result_TSV[i_g][i_d][i_b][0];
	  distribution_deviation_TSV[i_g][i_d][i_b] = xcontribution[0].distribution_deviation_TSV[i_g][i_d][i_b][0];
	}
      }
    }
  }


  ////////////////////////////////////////
  //  output for TSV variation started  //
  ////////////////////////////////////////

  logger << LOG_INFO << "output for TSV variation started" << endl;
  output_info();
  logger << LOG_INFO << "active_qTcut = " << active_qTcut << endl;
  if (ygeneric->switch_output_overview > 1){output_distribution_overview_TSV();}
  logger << LOG_INFO << "output_distribution_overview_TSV done  " << resultdirectory << " !" << endl;
  if (ygeneric->switch_output_plot > 1){output_contribution_distribution_norm_TSV();}
  logger << LOG_INFO << "output_contribution_distribution_norm_TSV done  " << resultdirectory << " !" << endl;
  if (ygeneric->switch_output_plot > 1){output_contribution_distribution_plot_TSV();}
  logger << LOG_INFO << "output_contribution_distribution_plot_TSV done  " << resultdirectory << " !" << endl;
  if (active_qTcut){if (ygeneric->switch_output_plot > 1){output_contribution_distribution_norm_qTcut_TSV();}}
  logger << LOG_INFO << "output_contribution_distribution_norm_qTcut_TSV done  " << resultdirectory << " !" << endl;
  if (active_qTcut){if (ygeneric->switch_output_plot > 1){output_contribution_distribution_plot_qTcut_TSV();}}
  logger << LOG_INFO << "output_contribution_distribution_plot_qTcut_TSV done for  " << resultdirectory << " !" << endl;
  //  output_contribution_distribution_check_TSV();
  
  logger << LOG_DEBUG << "finished" << endl;
}



