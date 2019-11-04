#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

//vector<vector<double> > & bin_edge
void summary_contribution::readin_distribution_contribution_CV(){
  Logger logger("summary_contribution::readin_distribution_contribution_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  if (yresultdirectory != ylist->resultdirectory){logger << LOG_INFO << "yresultdirectory = " << yresultdirectory << "  =/=  " << ylist->resultdirectory << " =  ylist->resultdirectory" << endl;}

  //  string result_file;
  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "readin_distribution_contribution_CV started:" << endl;
  logger << LOG_DEBUG << "ylist->resultdirectory = " << ylist->resultdirectory << endl;
  logger << LOG_DEBUG << "type_perturbative_order = " << type_perturbative_order << endl;
  logger << LOG_DEBUG << "type_subtraction_method = " << type_subtraction_method << endl;
  logger << LOG_DEBUG << "in_contribution_order_alpha_s = " << in_contribution_order_alpha_s << endl;
  logger << LOG_DEBUG << "in_contribution_order_alpha_e = " << in_contribution_order_alpha_e << endl;
  logger << LOG_DEBUG << "type_contribution = " << type_contribution << endl;
  logger << LOG_DEBUG << "type_correction = " << type_correction << endl;
  logger << LOG_DEBUG << "interference = " << interference << endl;
  logger << LOG_DEBUG << "photon_induced = " << photon_induced << endl;

  logger << LOG_DEBUG << "infix_order_contribution = " << infix_order_contribution << endl;

  logger << LOG_DEBUG << "osi->extended_distribution.size() = " << osi->extended_distribution.size() << endl;
  logger << LOG_DEBUG << "ygeneric->final_resultdirectory = " << ygeneric->final_resultdirectory << endl;
  logger << LOG_DEBUG << "ylist->resultdirectory = " << ylist->resultdirectory << endl;

  vector<double> n_bin_distribution(osi->extended_distribution.size());

  n_bin_distribution_modasym.resize(osi->extended_distribution.size());

  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    n_bin_distribution[i_d] = osi->extended_distribution[i_d].n_bins;
    n_bin_distribution_modasym[i_d] = osi->extended_distribution[i_d].n_bins;
    if (osi->extended_distribution[i_d].xdistribution_type == "asym1" || osi->extended_distribution[i_d].xdistribution_type == "asym2"){n_bin_distribution_modasym[i_d] = 2 * n_bin_distribution_modasym[i_d];}
  }
  logger << LOG_DEBUG_VERBOSE << "n_bin_distribution_modasym.size() = " << n_bin_distribution_modasym.size() << endl;



  ///////////////////////////////////////////
  //  declaration of subprocess_result_CV  //
  ///////////////////////////////////////////
  xsubprocess[0].distribution_result_CV.resize(osi->extended_distribution.size());
  xsubprocess[0].distribution_deviation_CV.resize(osi->extended_distribution.size());
  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    xsubprocess[0].distribution_result_CV[i_d].resize(n_bin_distribution_modasym[i_d], vector<double> (osi->n_scales_CV, 0.));
    xsubprocess[0].distribution_deviation_CV[i_d].resize(n_bin_distribution_modasym[i_d], vector<double> (osi->n_scales_CV, 0.));
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //  read-in of results, squared results and event numbers from CV distribution files  //
  ////////////////////////////////////////////////////////////////////////////////////////

  logger << LOG_DEBUG << "xsubprocess.size() = " << xsubprocess.size() << endl;
  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    logger << LOG_INFO << "READIN   " << infix_order_contribution << "   xsubprocess[" << setw(3) << i_p << "].name = " << xsubprocess[i_p].name << endl;
    xsubprocess[i_p].readin_distribution_CV();
  }

  logger << LOG_DEBUG_VERBOSE << "calculation begins" << endl;

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  calculation of results, squared results and event numbers by just adding up everything  //
  //////////////////////////////////////////////////////////////////////////////////////////////
  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    logger << LOG_INFO << "COMBINATION   " << infix_order_contribution << "   xsubprocess[" << setw(3) << i_p << "].name = " << xsubprocess[i_p].name << endl;
    xsubprocess[i_p].combination_distribution_CV();
  }

  logger << LOG_DEBUG_VERBOSE << "original/alternative calculation finished" << endl;



  distribution_result_CV.resize(ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->extended_distribution.size()));
  distribution_deviation_CV.resize(ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->extended_distribution.size()));
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      distribution_result_CV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins, vector<double> (osi->n_scales_CV, 0.));
      distribution_deviation_CV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins, vector<double> (osi->n_scales_CV, 0.));
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //  calculation of results and deviations for each subgroup of processes -> distribution_result_CV  //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      logger << LOG_DEBUG_VERBOSE << "osi->extended_distribution[" << setw(3) << i_d << "].xdistribution_name = " << osi->extended_distribution[i_d].xdistribution_name << endl;
      distribution_result_CV[i_g][i_d].resize(n_bin_distribution_modasym[i_d], vector<double> (osi->n_scales_CV, 0.));
      distribution_deviation_CV[i_g][i_d].resize(n_bin_distribution_modasym[i_d], vector<double> (osi->n_scales_CV, 0.));
      for (int i_b = 0; i_b < n_bin_distribution_modasym[i_d]; i_b++){
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  for (int i_p = 0; i_p < subgroup_no_member[i_g].size(); i_p++){
	    distribution_result_CV[i_g][i_d][i_b][i_s] += xsubprocess[subgroup_no_member[i_g][i_p]].distribution_result_CV[i_d][i_b][i_s];
	    distribution_deviation_CV[i_g][i_d][i_b][i_s] += pow(xsubprocess[subgroup_no_member[i_g][i_p]].distribution_deviation_CV[i_d][i_b][i_s], 2);
	  }
	  distribution_deviation_CV[i_g][i_d][i_b][i_s] = sqrt(distribution_deviation_CV[i_g][i_d][i_b][i_s]);
	  if (i_s == osi->n_scales_CV / 2){
	    logger << LOG_DEBUG_VERBOSE << "distribution_result_CV[" << i_g << "][" << setw(3) << i_d << "][" << setw(3) << i_b << "][" << setw(2) << i_s << "] = " << setw(25) << distribution_result_CV[i_g][i_d][i_b][i_s] << " +- " << distribution_deviation_CV[i_g][i_d][i_b][i_s] << endl;
	  }
	}
      }
    }
  }

  /////////////////////////////////////////////////////////////////////
  //  sum of subprocesses stored as entry 0 in subprocess_result_CV  //
  /////////////////////////////////////////////////////////////////////
  for (int i_q = 0; i_q < 1; i_q++){
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	  double dev = 0.;
	  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
	    xsubprocess[0].distribution_result_CV[i_d][i_b][i_s] += xsubprocess[i_p].distribution_result_CV[i_d][i_b][i_s];
	    dev += pow(xsubprocess[i_p].distribution_deviation_CV[i_d][i_b][i_s], 2);
	    //	    logger << LOG_DEBUG_VERBOSE << "devcheck: subprocess_deviation_CV[" << i_p << "][i_d][i_b][i_s] = " << osi->unit_factor_distribution * xsubprocess[i_p].distribution_deviation_CV[i_d][i_b][i_s] << "   dev = " << osi->unit_factor_distribution * dev << endl;
	  }
	  xsubprocess[0].distribution_deviation_CV[i_d][i_b][i_s] = sqrt(dev);
	  //	  logger << LOG_DEBUG_VERBOSE << "xdevcheck: subprocess_deviation_CV[0][i_d][i_b][i_s] = " << osi->unit_factor_distribution * subprocess_deviation_CV[0][i_d][i_b][i_s] << "   dev = " << osi->unit_factor_distribution * dev << endl;
	}
      }
    }
  }



  if (ygeneric->switch_output_overview > 2){output_distribution_overview_qTcut_CV();}
  if (ygeneric->switch_output_plot > 2){output_distribution_plot_CV();}

  //  deactivated: might need to be updated if needed again !!!
  //  output_distribution_asymmetry_CV();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
















//vector<vector<double> > & bin_edge
void summary_contribution::readin_distribution_contribution_TSV(){
  Logger logger("summary_contribution::readin_distribution_contribution_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  if (yresultdirectory != ylist->resultdirectory){logger << LOG_INFO << "yresultdirectory = " << yresultdirectory << "  =/=  " << ylist->resultdirectory << " =  ylist->resultdirectory" << endl;}

  logger << LOG_DEBUG << "type_contribution = " << type_contribution << "   osi->output_n_qTcut = " << osi->output_n_qTcut << endl;
  logger << LOG_DEBUG << "type_contribution = " << type_contribution << "   output_n_qTcut = " << output_n_qTcut << endl;
  logger << LOG_DEBUG << "type_contribution = " << type_contribution << "   selection_n_qTcut = " << selection_n_qTcut << endl;

  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "readin_distribution_contribution_TSV started:" << endl;
  logger << LOG_DEBUG << "ylist->resultdirectory = " << ylist->resultdirectory << endl;
  logger << LOG_DEBUG << "type_perturbative_order = " << type_perturbative_order << endl;
  logger << LOG_DEBUG << "type_subtraction_method = " << type_subtraction_method << endl;
  logger << LOG_DEBUG << "in_contribution_order_alpha_s = " << in_contribution_order_alpha_s << endl;
  logger << LOG_DEBUG << "in_contribution_order_alpha_e = " << in_contribution_order_alpha_e << endl;
  logger << LOG_DEBUG << "type_contribution = " << type_contribution << endl;
  logger << LOG_DEBUG << "type_correction = " << type_correction << endl;
  logger << LOG_DEBUG << "interference = " << interference << endl;
  logger << LOG_DEBUG << "photon_induced = " << photon_induced << endl;

  //  string result_file;
  //  string readindata;

  //  int n_seed = directory.size();
  //  int n_subprocess = subprocess.size();
  /*

  logger << LOG_DEBUG_VERBOSE << "directory.size() = " << directory.size() << endl;
  for (int ic = 0; ic < directory.size(); ic++){
    logger << LOG_DEBUG << "directory[" << ic << "] = " << directory[ic] << endl;
  }
  logger << LOG_DEBUG << endl;

  logger << LOG_DEBUG_VERBOSE << "osi->extended_distribution.size() = " << osi->extended_distribution.size() << endl;
  */

  vector<int> n_bin_distribution(osi->extended_distribution.size(), 0);

  n_bin_distribution_modasym.resize(osi->extended_distribution.size(), 0);
  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    n_bin_distribution[i_d] = osi->extended_distribution[i_d].n_bins;
    n_bin_distribution_modasym[i_d] = osi->extended_distribution[i_d].n_bins;
    if (osi->extended_distribution[i_d].xdistribution_type == "asym1" || osi->extended_distribution[i_d].xdistribution_type == "asym2"){n_bin_distribution_modasym[i_d] = 2 * n_bin_distribution_modasym[i_d];}
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  calculation of results and deviations for each subprocess -> subprocess_result_TSV from singlerun_result_TSV  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i_p = 1; i_p < subprocess.size(); i_p++){
    logger << LOG_INFO << left << "READIN        " << infix_contribution << "   xsubprocess[" << setw(3) << i_p << "].name = " << xsubprocess[i_p].name << right << "   (" << active_qTcut << " - " << setw(3) << selection_n_qTcut << " / " << output_n_qTcut << ")" << endl;
    //    logger << LOG_INFO << "              " << type_contribution << "." << type_correction << "   " << "   active_qTcut = " << active_qTcut << "   output_n_qTcut = " << output_n_qTcut << "   selection_n_qTcut = " << selection_n_qTcut << endl;

    xsubprocess[i_p].readin_distribution_TSV();

    logger << LOG_INFO << left << "COMBINATION   " << infix_contribution << "   xsubprocess[" << setw(3) << i_p << "].name = " << xsubprocess[i_p].name << right << "   (" << active_qTcut << " - " << setw(3) << selection_n_qTcut << " / " << output_n_qTcut << ")" << endl;

    for (int i_m = 0; i_m < extended_directory.size(); i_m++){
      for (int i_z = 0; i_z < extended_directory[i_m].size(); i_z++){
	for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	  if (!ygeneric->switch_output_distribution[i_d]){continue;}
	  for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
	    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
		for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
		  if (i_r == 1 && i_f == 1){
		    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
		      logger << LOG_DEBUG
			     << "distribution_group_run_result_TSV[" << i_m << "][" << setw(4) << i_z << "][" << setw(2) << i_d << "][" << setw(3) << i_b << "][" << setw(2) << x_q << "][" << setw(2) << i_s << "][" << setw(1) << i_r << "][" << setw(1) << i_f << "] = "
			     << xsubprocess[i_p].distribution_group_run_result_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f]
			     << " +- "
			     << xsubprocess[i_p].distribution_group_run_deviation_TSV[i_m][i_z][i_d][i_b][x_q][i_s][i_r][i_f]
			     << "   ("
			     << xsubprocess[i_p].distribution_group_run_N_binwise_TSV[i_m][i_z][i_d][i_b][x_q]
			     << ")"
			     << endl;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    xsubprocess[i_p].combination_distribution_TSV();
    
    int x_s = osi->no_reference_TSV;
    int x_r = osi->no_scale_ren_reference_TSV;
    int x_f = osi->no_scale_fact_reference_TSV;

    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      if (osi->extended_distribution[i_d].xdistribution_name != "total_rate"){continue;}
      for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	  logger << LOG_INFO
		 << "distribution_result_TSV"
		 << "[" << osi->extended_distribution[i_d].xdistribution_name << "][" << setw(3) << i_b << "][" << setw(2) << x_q << "][" << setw(2) << x_s << "][" << setw(1) << x_r << "][" << setw(1) << x_f << "] = "
	    //		 << "[" << setw(2) << i_d << "][" << setw(3) << i_b << "][" << setw(2) << x_q << "][" << setw(2) << x_s << "][" << setw(1) << x_r << "][" << setw(1) << x_f << "] = "
		 << xsubprocess[i_p].distribution_result_TSV[i_d][i_b][x_q][x_s][x_r][x_f]
		 << " +- "
		 << xsubprocess[i_p].distribution_deviation_TSV[i_d][i_b][x_q][x_s][x_r][x_f]
		 << endl;
	}
      }
    }

    /*
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      if (osi->extended_distribution[i_d].xdistribution_name != "total_rate"){continue;}
      for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
		logger << LOG_INFO
		       << "distribution_result_TSV"
		       << "[" << setw(2) << i_d << "][" << setw(3) << i_b << "][" << setw(2) << x_q << "][" << setw(2) << i_s << "][" << setw(1) << i_r << "][" << setw(1) << i_f << "] = "
		       << xsubprocess[i_p].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f]
		       << " +- "
		       << xsubprocess[i_p].distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f]
		       << endl;
		
	      }
	    }
	  }
	}
      }
    }
     */
  }

  // subprocess_result_TSV[0] has to be resized !!!
  xsubprocess[0].distribution_result_TSV.resize(osi->extended_distribution.size());
  xsubprocess[0].distribution_deviation_TSV.resize(osi->extended_distribution.size());
  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    xsubprocess[0].distribution_result_TSV[i_d].resize(osi->extended_distribution[i_d].n_bins);
    xsubprocess[0].distribution_deviation_TSV[i_d].resize(osi->extended_distribution[i_d].n_bins);
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
      xsubprocess[0].distribution_result_TSV[i_d][i_b].resize(selection_n_qTcut);
      xsubprocess[0].distribution_deviation_TSV[i_d][i_b].resize(selection_n_qTcut);
      for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
	xsubprocess[0].distribution_result_TSV[i_d][i_b][x_q].resize(osi->n_extended_set_TSV);
	xsubprocess[0].distribution_deviation_TSV[i_d][i_b][x_q].resize(osi->n_extended_set_TSV);
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  if (!osi->switch_distribution_TSV[i_s]){continue;}
	  xsubprocess[0].distribution_result_TSV[i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  xsubprocess[0].distribution_deviation_TSV[i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	}
      }
    }
  }

  //////////////////////////////////////////////////////////////////////
  //  sum of subprocesses stored as entry 0 in subprocess_result_TSV  //
  //////////////////////////////////////////////////////////////////////

  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
      for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  if (!osi->switch_distribution_TSV[i_s]){continue;}
	  //if (osi->switch_distribution_TSV[i_s] == 0){continue;}
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      double dev = 0.;
	      for (int i_p = 0; i_p < subprocess.size(); i_p++){
		xsubprocess[0].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f] += xsubprocess[i_p].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f];
		dev += pow(xsubprocess[i_p].distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f], 2);
	      }
	      xsubprocess[0].distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f] = sqrt(dev);
	    }
	  }
	}
      }
    }
  }

  //////////////////////////////////////////
  //  declaration of subgroup_result_TSV  //
  //////////////////////////////////////////

  distribution_result_TSV.resize(ygeneric->subgroup.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->extended_distribution.size()));
  distribution_deviation_TSV.resize(ygeneric->subgroup.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->extended_distribution.size()));
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      distribution_result_TSV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
      distribution_deviation_TSV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
      for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	distribution_result_TSV[i_g][i_d][i_b].resize(selection_n_qTcut);
	distribution_deviation_TSV[i_g][i_d][i_b].resize(selection_n_qTcut);
	for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
	  distribution_result_TSV[i_g][i_d][i_b][x_q].resize(osi->n_extended_set_TSV);
	  distribution_deviation_TSV[i_g][i_d][i_b][x_q].resize(osi->n_extended_set_TSV);
	  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    if (!osi->switch_distribution_TSV[i_s]){continue;}
	    //	    if (osi->switch_distribution_TSV[i_s] != 0){
	    distribution_result_TSV[i_g][i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    //	    }
	  }
	}
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //  calculation of results and deviations for each subgroup of processes -> distribution_result_TSV  //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
      for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  if (!osi->switch_distribution_TSV[i_s]){continue;}
	  //	  if (osi->switch_distribution_TSV[i_s] == 0){continue;}
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
		for (int i_p = 0; i_p < subgroup_no_member[i_g].size(); i_p++){
		  distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] += xsubprocess[subgroup_no_member[i_g][i_p]].distribution_result_TSV[i_d][i_b][x_q][i_s][i_r][i_f];
		  distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] += pow(xsubprocess[subgroup_no_member[i_g][i_p]].distribution_deviation_TSV[i_d][i_b][x_q][i_s][i_r][i_f], 2);
		}
		distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] = sqrt(distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f]);
	      }
	    }
	  }
	}
      }
    }
  }
  
  //////////////////////////////////////////////
  //  test output of distribution_result_TSV  //
  //////////////////////////////////////////////

  /*
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
      if (i_d == 0){
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	  for (int x_q = 0; x_q < selection_n_qTcut; x_q++){
	    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    if (!osi->switch_distribution_TSV[i_s]){continue;}
//	      if (osi->switch_distribution_TSV[i_s] != 0){
		for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
		  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
		    if (distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] != 0.){
		      logger << LOG_DEBUG_VERBOSE << "distribution_result_TSV[" << setw(2) << i_g << "][" << setw(2) << i_d << "][" << setw(2) << i_b << "][" << setw(2) << x_q << "][" << setw(2) << i_s << "][" << setw(2) << i_r << "][" << setw(2) << i_f << "] = " << setw(23) << setprecision(15) << distribution_result_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << distribution_deviation_TSV[i_g][i_d][i_b][x_q][i_s][i_r][i_f] << endl;
		    }
		  }
		}
//	      }
	    }
	  }
	}
      }
    }
  }
  */

  if (ygeneric->switch_output_overview > 2){output_distribution_overview_qTcut_TSV();}
  if (ygeneric->switch_output_plot > 2){output_distribution_plot_qTcut_TSV();}

  pid_t pid = getpid();
  stringstream temp_pid;
  temp_pid << "cat /proc/" << pid << "/status | grep VmSize";
  string temp_s_pid = temp_pid.str();
  logger << LOG_INFO << "begin swapping xsubprocess content:   " << temp_s_pid << endl;
  int temp_i = system(temp_s_pid.c_str());

  for (int i_p = 1; i_p < xsubprocess.size(); i_p++){
    vector<vector<vector<vector<vector<vector<double> > > > > > ().swap(xsubprocess[i_p].distribution_result_TSV);
    vector<vector<vector<vector<vector<vector<double> > > > > > ().swap(xsubprocess[i_p].distribution_deviation_TSV);
    vector<vector<vector<vector<vector<vector<double> > > > > > ().swap(xsubprocess[i_p].distribution_chi2_TSV);
    
    vector<vector<vector<vector<int> > > > ().swap(xsubprocess[i_p].distribution_group_counter_removal_run_qTcut_TSV);
    vector<vector<vector<vector<int> > > > ().swap(xsubprocess[i_p].distribution_group_counter_nonzero_run_qTcut_TSV);

    vector<vector<vector<vector<vector<long long> > > > > ().swap(xsubprocess[i_p].distribution_group_run_N_binwise_TSV);

    vector<vector<vector<vector<long long> > > > ().swap(xsubprocess[i_p].distribution_group_N_binwise_TSV);
    vector<vector<vector<vector<vector<vector<vector<double> > > > > > > ().swap(xsubprocess[i_p].distribution_group_result_TSV);
    vector<vector<vector<vector<vector<vector<vector<double> > > > > > > ().swap(xsubprocess[i_p].distribution_group_deviation_TSV);
    vector<vector<vector<vector<vector<vector<vector<double> > > > > > > ().swap(xsubprocess[i_p].distribution_group_chi2_TSV);

    logger << LOG_INFO << "swap(distribution_N_total_binwise_TSV) done:   " << infix_contribution << "   " << xsubprocess[i_p].name << endl;
    temp_i = system(temp_s_pid.c_str());
  }
  
  logger << LOG_DEBUG_VERBOSE << "contribution_TSV finished" << endl;
}



