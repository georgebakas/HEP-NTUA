#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_list::collect_contribution_result_CV(){
  Logger logger("summary_list::collect_contribution_result_CV");
  logger << LOG_DEBUG << "called" << endl;

  string filename;

  system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV");
  system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV);

  logger << LOG_DEBUG<< "XXX ygeneric->final_resultdirectory/CV/ygeneric->name_variation_CV = " << ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV << endl;
  logger << LOG_DEBUG<< "resultdirectory = " << resultdirectory << endl;

  system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory);
  //  system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/group");
  //  system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/subprocesses");

  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/" + xcontribution[i_c].infix_contribution);
    //    system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/" + xcontribution[i_c].infix_contribution + "/subprocesses");
    //    system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/" + xcontribution[i_c].infix_contribution + "/group");
  }

  string muscale;
  if (osi->scale_ren == osi->msi.M_W){muscale = "M_W";}
  else if (osi->scale_ren == osi->msi.M_Z){muscale = "M_Z";}
  //  else if (osi->scale_ren == osi->msi.M_H){muscale = "M_H";}
  //  else if (osi->scale_ren == osi->msi.M_t){muscale = "M_t";}

  string result_moment;

  xcontribution[0].result_CV.resize(osi->n_moments + 1, vector<vector<vector<double> > > (ygeneric->subgroup.size(), vector<vector<double> > (osi->n_qTcut, vector<double> (osi->n_scales_CV, 0.))));
  xcontribution[0].deviation_CV.resize(osi->n_moments + 1, vector<vector<vector<double> > > (ygeneric->subgroup.size(), vector<vector<double> > (osi->n_qTcut, vector<double> (osi->n_scales_CV, 0.))));

  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
    stringstream temp_ss;
    if (i_m == osi->n_moments){temp_ss << "result";}
    else {temp_ss << "moment_" << i_m;}
    result_moment = temp_ss.str();
    for (int i_c = 1; i_c < xcontribution.size(); i_c++){
      xcontribution[i_c].readin_result_contribution_CV(result_moment, i_m);
    }
  }

  logger << LOG_DEBUG_VERBOSE << "calculate full contributions" << endl;

  //////////////////////////////////////
  //   calculate full contributions   //
  //////////////////////////////////////

  //  if (osi->n_scales_CV != 0 && osi->n_qTcut != 0){
  if (osi->switch_CV){
    for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	    double dev = 0.;
	    for (int i_c = 1; i_c < xcontribution.size(); i_c++){
	      xcontribution[0].result_CV[i_m][i_g][i_q][i_s] += xcontribution[i_c].result_CV[i_m][i_g][i_q][i_s];
	      dev += pow(xcontribution[i_c].deviation_CV[i_m][i_g][i_q][i_s], 2);
	      //	      logger << LOG_DEBUG_VERBOSE << "   i_m = " << i_m << "   i_s = " << i_s << "   i_q = " << i_q << "   i_g = " << i_g << "   i_c = " << i_c << endl;
	    }
	    xcontribution[0].deviation_CV[i_m][i_g][i_q][i_s] = sqrt(dev);
	    if (!(i_m == 0 || i_m == osi->n_moments)){
	      xcontribution[0].deviation_CV[i_m][i_g][i_q][i_s] = (xcontribution[0].result_CV[i_m][i_g][i_q][i_s] * xcontribution[0].deviation_CV[0][0][i_q][i_s] + xcontribution[0].result_CV[0][0][i_q][i_s] * xcontribution[0].deviation_CV[i_m][i_g][i_q][i_s]) / pow(xcontribution[0].result_CV[0][0][i_q][i_s], 2);
	      xcontribution[0].result_CV[i_m][i_g][i_q][i_s] = xcontribution[0].result_CV[i_m][i_g][i_q][i_s] / xcontribution[0].result_CV[0][0][i_q][i_s];
	    }
	  }
	}
      }
    }
  }
  
  //////////////////////////////
  //   calculation finished   //
  //////////////////////////////
  logger << LOG_DEBUG_VERBOSE << "calculation finished" << endl;

  /////////////////////////////////////////////////////////////////////
  //  resize result_CV (to be filled with extrapolation qTcut -> 0)  //
  /////////////////////////////////////////////////////////////////////

  result_CV.resize(osi->n_moments + 1, vector<vector<double> > (ygeneric->subgroup.size(), vector<double> (osi->n_scales_CV, 0.)));
  deviation_CV.resize(osi->n_moments + 1, vector<vector<double> > (ygeneric->subgroup.size(), vector<double> (osi->n_scales_CV, 0.)));
  deviation_extrapolation_CV.resize(osi->n_moments + 1, vector<vector<double> > (ygeneric->subgroup.size(), vector<double> (osi->n_scales_CV, 0.)));
  deviation_statistics_CV.resize(osi->n_moments + 1, vector<vector<double> > (ygeneric->subgroup.size(), vector<double> (osi->n_scales_CV, 0.)));

  ////////////////////////////////
  //  extrapolation qTcut -> 0  //
  ////////////////////////////////

  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
    if (i_m > 0){continue;} // no moments implemented so far !!!
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      if (active_qTcut){
	extrapolation_CV(ygeneric->switch_extrapolation_result);
	//	extrapolation_CV(osi->value_qTcut, xcontribution[0].result_CV[i_m][i_g], xcontribution[0].deviation_CV[i_m][i_g], result_CV[i_m][i_g], deviation_CV[i_m][i_g], deviation_statistics_CV[i_m][i_g], deviation_extrapolation_CV[i_m][i_g], 2);
      }
      else {
	result_CV[i_m][i_g] = xcontribution[0].result_CV[i_m][i_g][0];
	deviation_CV[i_m][i_g] = xcontribution[0].deviation_CV[i_m][i_g][0];
	deviation_statistics_CV[i_m][i_g] = xcontribution[0].deviation_CV[i_m][i_g][0];
	//	deviation_extrapolation_CV[i_m][i_g] = vector<double> (osi->n_scales_CV, 0.);
      }
    }
  }


  /*
  vector<double> rel_scale_CV(osi->n_scales_CV);
  vector<double> scale_CV(osi->n_scales_CV);
  if (osi->n_scales_CV != 0 && osi->n_qTcut != 0){
    for (int s = 0; s < osi->n_scales_CV; s++){
      if (osi->n_scales_CV > 1){rel_scale_CV[s] = pow(10., log10(double(osi->variation_factor_CV)) * double(2. * s / (osi->n_scales_CV - 1) - 1));}
      else {rel_scale_CV[s] = 1;}
      scale_CV[s] = rel_scale_CV[s] * osi->scale_ren;
    }
  }
  */

  
  if (ygeneric->switch_output_overview > 1){output_contribution_result_overview_qTcut_CV();}
  //  if (ygeneric->switch_output_result > 1){output_contribution_result_CV();}
  if (ygeneric->switch_output_result > 1){output_contribution_result_qTcut_CV();}
  if (ygeneric->switch_output_plot > 1){output_contribution_result_plot_CV();}
  if (ygeneric->switch_output_plot > 1){output_contribution_result_plot_qTcut_CV();}

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::collect_contribution_result_TSV(){
  Logger logger("summary_list::collect_contribution_result_TSV");
  logger << LOG_DEBUG << "called" << endl;

  string filename;

  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s]);
  }
  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    //      system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/" + resultdirectory + "/" + osi->name_extended_set_TSV[i_s]);
    for (int i_v = 0; i_v < ygeneric->name_scale_variation_TSV.size(); i_v++){
      system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[i_v]);
      system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[i_v] + "/" + resultdirectory);
    }
  }

  string muscale;
  if (osi->scale_ren == osi->msi.M_W){muscale = "M_W";}
  else if (osi->scale_ren == osi->msi.M_Z){muscale = "M_Z";}
  //  else if (osi->scale_ren == M_H){muscale = "M_H";}
  //  else if (osi->scale_ren == M_t){muscale = "M_t";}

  string result_moment;

  xcontribution[0].result_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<vector<double> > > > > (ygeneric->subgroup.size(), vector<vector<vector<vector<double> > > > (osi->n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV))));
  xcontribution[0].deviation_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<vector<double> > > > > (ygeneric->subgroup.size(), vector<vector<vector<vector<double> > > > (osi->n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV))));
  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  for (int i_c = 0; i_c < xcontribution.size(); i_c++){
	    xcontribution[0].result_TSV[i_m][i_g][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    xcontribution[0].deviation_TSV[i_m][i_g][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  }
	}
      }
    }
  }

  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
    stringstream temp_ss;
    if (i_m == osi->n_moments){temp_ss << "result";}
    else {temp_ss << "moment_" << i_m;}
    result_moment = temp_ss.str();
    for (int i_c = 1; i_c < xcontribution.size(); i_c++){
      xcontribution[i_c].readin_result_contribution_TSV(result_moment, i_m);
    }
  }

  logger << LOG_DEBUG_VERBOSE << "calculate full contributions" << endl;

  //////////////////////////////////////
  //   calculate full contributions   //
  //////////////////////////////////////

  if (osi->switch_TSV){
    for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
      logger << LOG_DEBUG_VERBOSE << " TSV: calculate full contributions start" << endl;
      if (i_m > 0){continue;} // no moments implemented so far !!!
      for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
		double dev = 0.;
		for (int i_c = 1; i_c < xcontribution.size(); i_c++){
		  xcontribution[0].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f] += xcontribution[i_c].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f];
		  dev += pow(xcontribution[i_c].deviation_TSV[i_m][i_g][i_q][i_s][i_r][i_f], 2);
		}
		//		logger << LOG_DEBUG_VERBOSE << "dev = " << dev << endl;
		xcontribution[0].deviation_TSV[i_m][i_g][i_q][i_s][i_r][i_f] = sqrt(dev);
		if (!(i_m == 0 || i_m == osi->n_moments)){
		  logger << LOG_DEBUG_VERBOSE << "Should not happen!" << endl;
		  xcontribution[0].deviation_TSV[i_m][i_g][i_q][i_s][i_r][i_f] = (xcontribution[0].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f] * xcontribution[0].deviation_TSV[0][0][i_q][i_s][i_r][i_f] + xcontribution[0].result_TSV[0][0][i_q][i_s][i_r][i_f] * xcontribution[0].deviation_TSV[i_m][i_g][i_q][i_s][i_r][i_f]) / pow(xcontribution[0].result_TSV[0][0][i_q][i_s][i_r][i_f], 2);
		  xcontribution[0].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f] = xcontribution[0].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f] / xcontribution[0].result_TSV[0][0][i_q][i_s][i_r][i_f];
		}
	      }
	    }
	  }
	}
      }
    }
  }

  //////////////////////////////
  //   calculation finished   //
  //////////////////////////////
  logger << LOG_DEBUG_VERBOSE << "calculation finished" << endl;



  //////////////////////////////////////////////////////////////////////
  //  resize result_TSV (to be filled with extrapolation qTcut -> 0)  //
  //////////////////////////////////////////////////////////////////////

  result_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<double> > > > (ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->n_extended_set_TSV)));
  deviation_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<double> > > > (ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->n_extended_set_TSV)));
  deviation_extrapolation_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<double> > > > (ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->n_extended_set_TSV)));
  deviation_statistics_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<double> > > > (ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->n_extended_set_TSV)));
  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	result_TSV[i_m][i_g][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	deviation_TSV[i_m][i_g][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	deviation_extrapolation_TSV[i_m][i_g][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	deviation_statistics_TSV[i_m][i_g][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
      }
    }
  }

  ////////////////////////////////
  //  extrapolation qTcut -> 0  //
  ////////////////////////////////

  //  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
  //    if (i_m > 0){continue;} // no moments implemented so far !!!
  //    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
  if (active_qTcut && osi->n_qTcut > 1){
    extrapolation_TSV(ygeneric->switch_extrapolation_result);
    //    extrapolation_TSV(osi->value_qTcut, xcontribution[0].result_TSV[i_m][i_g], xcontribution[0].deviation_TSV[i_m][i_g], result_TSV[i_m][i_g], deviation_statistics_TSV[i_m][i_g], deviation_TSV[i_m][i_g], deviation_extrapolation_TSV[i_m][i_g], 2);
  }
  else {
    for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
      if (i_m > 0){continue;} // no moments implemented so far !!!
      for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	result_TSV[i_m][i_g] = xcontribution[0].result_TSV[i_m][i_g][0];
	deviation_TSV[i_m][i_g] = xcontribution[0].deviation_TSV[i_m][i_g][0];
	deviation_statistics_TSV[i_m][i_g] = xcontribution[0].deviation_TSV[i_m][i_g][0];
	//	deviation_extrapolation_TSV[i_m][i_g] = 0.;
      }
    }
  }



  /*
  //  only relevant if moment scans are performed (in that case, no output files are created for the moments)
  int start_i_m = 0;
  if (osi->user.switch_map["n_scan_variable"] != 0){
    logger << LOG_DEBUG<< "user_switch_map[n_scan_variable] = " << osi->user.switch_map["n_scan_variable"] << "   user_switch_value[n_scan_variable] = " << osi->user.switch_value[osi->user.switch_map["n_scan_variable"]] << endl;
  }
  if (osi->user.switch_map["n_scan_gridsize"] != 0){logger << LOG_DEBUG<<"user_switch_map[n_scan_gridsize] = " << osi->user.switch_map["n_scan_gridsize"] << "   user_switch_value[n_scan_gridsize] = " << osi->user.switch_value[osi->user.switch_map["n_scan_gridsize"]] << endl;}

  if (osi->n_moments > 50){
    logger << LOG_DEBUG<< "Moments are considered as a parameter scan!" << endl;
    start_i_m = osi->n_moments;
  }
  logger << LOG_DEBUG_VERBOSE << "start_i_m = " << start_i_m << endl;
  */



  if (ygeneric->switch_output_overview > 1){output_contribution_result_overview_qTcut_TSV();}
  //  if (ygeneric->switch_output_result > 1){output_contribution_result_TSV();}
  if (ygeneric->switch_output_result > 1){output_contribution_result_qTcut_TSV();}
  if (ygeneric->switch_output_plot > 1){output_contribution_result_plot_TSV();}
  if (ygeneric->switch_output_plot > 1){output_contribution_result_plot_qTcut_TSV();}

  logger << LOG_DEBUG << "finished" << endl;
}



