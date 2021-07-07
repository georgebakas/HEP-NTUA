#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_order::collect_distribution_CV(){
  Logger logger("summary_order::collect_distribution_CV");
  logger << LOG_DEBUG << "called" << endl;

  distribution_result_CV.resize(ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->extended_distribution.size()));
  distribution_deviation_CV.resize(ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->extended_distribution.size()));
  for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
    for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
      if (!ygeneric->switch_output_distribution[i_d]){continue;}
      distribution_result_CV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins, vector<double> (osi->n_scales_CV, 0.));
      distribution_deviation_CV[i_g][i_d].resize(osi->extended_distribution[i_d].n_bins, vector<double> (osi->n_scales_CV, 0.));
    }
  }
  
  if (combination.size() == 1){
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	    double dev = 0.;
	    for (int i_l = 0; i_l < xlist.size(); i_l++){logger << LOG_DEBUG_VERBOSE << "i_l = " << i_l << endl; 
	      ///	    for (int i_l = 0; i_l < contribution_file.size(); i_l++){
	      ///	      int x_l = mapping_contribution_file[contribution_file[i_l]];
	      distribution_result_CV[i_g][i_d][i_b][i_s] += xlist[i_l]->xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s];
	      dev += pow(xlist[i_l]->xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s], 2);
	    }
	    distribution_deviation_CV[i_g][i_d][i_b][i_s] = sqrt(dev);

	    if (osi->extended_distribution[i_d].xdistribution_name == "total_rate"){
	      stringstream temp_res;
	      temp_res << setw(23) << setprecision(15) << distribution_result_CV[i_g][i_d][i_b][i_s];
	      stringstream temp_dev;
	      temp_dev << setw(23) << setprecision(15) << distribution_deviation_CV[i_g][i_d][i_b][i_s];
	      //	      if (i_b == 0 && i_s == 0){
		logger << LOG_INFO << resultdirectory << "   result_CV[" << i_g << "][" << i_d << "][" << i_b << "][" << i_s << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
		//	      }
	      logger << LOG_DEBUG_VERBOSE << "XXX   result_CV[" << i_g << "][" << i_d << "][" << i_b << "][" << i_s << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
	    }
	  }
	}
      }
    }
  }
  else if (combination.size() > 1){
    vector<vector<vector<vector<vector<double> > > > > factor_order_result_CV(combination.size(), vector<vector<vector<vector<double> > > > (ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->extended_distribution.size())));
    vector<vector<vector<vector<vector<double> > > > > factor_order_deviation_CV(combination.size(), vector<vector<vector<vector<double> > > > (ygeneric->subgroup.size(), vector<vector<vector<double> > > (osi->extended_distribution.size())));
    for (int i_c = 0; i_c < combination.size(); i_c++){
      for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	  if (!ygeneric->switch_output_distribution[i_d]){continue;}
	  factor_order_result_CV[i_c][i_g][i_d].resize(osi->extended_distribution[i_d].n_bins, vector<double> (osi->n_scales_CV, 0.));
	  factor_order_deviation_CV[i_c][i_g][i_d].resize(osi->extended_distribution[i_d].n_bins, vector<double> (osi->n_scales_CV, 0.));
	}
      }
    }


    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	    for (int i_c = 0; i_c < combination.size(); i_c++){logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << endl;
	      double result = 0.;
	      double dev2 = 0.;
	      for (int j_c = 0; j_c < combination[i_c].size(); j_c++){
		int i_l = combination[i_c][j_c];
		result += xlist[i_l]->xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s];
		dev2 += pow(xlist[i_l]->xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s], 2);
	      }
	      factor_order_result_CV[i_c][i_g][i_d][i_b][i_s] = result;
	      factor_order_deviation_CV[i_c][i_g][i_d][i_b][i_s] = sqrt(dev2);
	    }
	    /*	      
	    stringstream temp_res;
	    temp_res << setw(23) << setprecision(15) << distribution_result_CV[i_g][i_d][i_b][i_s];
	    stringstream temp_dev;
	    temp_dev << setw(23) << setprecision(15) << distribution_deviation_CV[i_g][i_d][i_b][i_s];
	    logger << LOG_DEBUG_VERBOSE << "XXX   result_CV[" << i_o << "][" << i_g << "][" << i_d << "][" << i_b << "][" << i_s << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
	    */
	  }
	}
      }
    }
    
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	    double finalized_result = 0.;
	    double finalized_deviation2 = 0.;
	    double index_LO = 0;
	    double temp_result = 0.;
	    double temp_deviation2 = 0.;
	    for (int i_c = 0; i_c < combination.size(); i_c++){
	      if (combination_type[i_c] == 0){
		if (i_c != 0){
		  finalized_result += temp_result;
		  temp_result = 0.;
		  finalized_deviation2 += temp_deviation2;
		  temp_deviation2 = 0.;
		}
		index_LO = i_c;
		temp_result = factor_order_result_CV[i_c][i_g][i_d][i_b][i_s];
		temp_deviation2 += pow(factor_order_deviation_CV[i_c][i_g][i_d][i_b][i_s], 2);
	      }
	      else if (combination_type[i_c] == 1){
		if (factor_order_result_CV[index_LO][i_g][i_d][i_b][i_s] == 0.){
		  temp_result += factor_order_result_CV[i_c][i_g][i_d][i_b][i_s];
		  // simplification in case of multiplicative combination:
		  temp_deviation2 += pow(factor_order_deviation_CV[i_c][i_g][i_d][i_b][i_s], 2);
		}
		else {
		  temp_result *= (1. + factor_order_result_CV[i_c][i_g][i_d][i_b][i_s] / factor_order_result_CV[index_LO][i_g][i_d][i_b][i_s]);
		  // simplification in case of multiplicative combination:
		  temp_deviation2 += pow(factor_order_deviation_CV[i_c][i_g][i_d][i_b][i_s], 2);
		}
	      }
	    }
	    distribution_result_CV[i_g][i_d][i_b][i_s] = finalized_result + temp_result;
	    distribution_deviation_CV[i_g][i_d][i_b][i_s] = sqrt(finalized_deviation2 + temp_deviation2);
	  }
	}
      }
    }

  }
  
  if (ygeneric->switch_output_overview > 0){output_distribution_overview_CV();}
  if (ygeneric->switch_output_plot > 0){output_distribution_CV();}
  /*
  if (ygeneric->switch_output_overview > 0){output_distribution_overview_TSV();}
  if (ygeneric->switch_output_plot > 0){output_distribution_qTcut_TSV();}
  if (ygeneric->switch_output_plot > 0){output_distribution_TSV();}
  */
  
  logger << LOG_DEBUG << "finished" << endl;
}

//void summary_order::collect_distribution_TSV(string & final_resultdirectory, vector<string> & subgroup, vector<xdistribution> & extended_distribution, vector<double> & fakeasymfactor, observable_set & oset, vector<vector<vector<string> > > & scalename_TSV){
void summary_order::collect_distribution_TSV(){
  //  removed: vector<double> & fakeasymfactor
  Logger logger("summary_order::collect_distribution_TSV");
  logger << LOG_DEBUG << "called" << endl;

  /*
  vector<vector<double> > bin_edge(osi->extended_distribution.size());
  for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
    if (!ygeneric->switch_output_distribution[i_d]){continue;}
    bin_edge[i_d].resize(osi->extended_distribution[i_d].n_bins + 1);
    for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins + 1; i_b++){bin_edge[i_d][i_b] = osi->extended_distribution[i_d].start + i_b * osi->extended_distribution[i_d].step;}
  }
  */
    
  ////////////////////////////////////////////////
  //  resize distribution_result/deviation_TSV  //
  ////////////////////////////////////////////////

  logger << LOG_DEBUG_VERBOSE << "before   distribution_result_TSV" << endl;

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
	  distribution_result_TSV[i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  distribution_deviation_TSV[i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	}
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "after   distribution_result_TSV" << endl;

  //////////////////////////////////////////////////////
  //  resize distribution_result/deviation_qTcut_TSV  //
  //////////////////////////////////////////////////////

  logger << LOG_DEBUG_VERBOSE << "before   distribution_result_qTcut_TSV" << endl;
  
  distribution_result_qTcut_TSV.resize(selection_n_qTcut, vector<vector<vector<vector<vector<vector<double> > > > > > (ygeneric->subgroup.size(), vector<vector<vector<vector<vector<double> > > > > (osi->extended_distribution.size())));
  distribution_deviation_qTcut_TSV.resize(selection_n_qTcut, vector<vector<vector<vector<vector<vector<double> > > > > > (ygeneric->subgroup.size(), vector<vector<vector<vector<vector<double> > > > > (osi->extended_distribution.size())));
  for (int i_q = 0; i_q < selection_n_qTcut; i_q++){
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	distribution_result_qTcut_TSV[i_q][i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
	distribution_deviation_qTcut_TSV[i_q][i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	  distribution_result_qTcut_TSV[i_q][i_g][i_d][i_b].resize(osi->n_extended_set_TSV);
	  distribution_deviation_qTcut_TSV[i_q][i_g][i_d][i_b].resize(osi->n_extended_set_TSV);
	  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    if (!osi->switch_distribution_TSV[i_s]){continue;}
	    distribution_result_qTcut_TSV[i_q][i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    distribution_deviation_qTcut_TSV[i_q][i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  }
	}
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "after   distribution_result_qTcut_TSV" << endl;

  ///////////////////////////////////////////////////////////////
  //  calculate distribution_result/deviation_TSV  //
  ///////////////////////////////////////////////////////////////

  //  for (int i_o = 0; i_o < yorder.size(); i_o++){logger << LOG_DEBUG_VERBOSE << "i_o = " << i_o << endl;
  if (combination.size() == 1){
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    if (!osi->switch_distribution_TSV[i_s]){continue;}

	    /*
	    logger << LOG_DEBUG_VERBOSE << "distribution_result_TSV[i_g][i_d][i_b].size() = " << distribution_result_TSV[i_g][i_d][i_b].size() << endl;
	    logger << LOG_DEBUG_VERBOSE << "distribution_deviation_TSV[i_g][i_d][i_b].size() = " << distribution_deviation_TSV[i_g][i_d][i_b].size() << endl;
	    logger << LOG_DEBUG_VERBOSE << "osi->n_scale_fact_TSV.size() = " << osi->n_scale_fact_TSV.size() << endl;
	    logger << LOG_DEBUG_VERBOSE << "osi->n_scale_ren_TSV.size() = " << osi->n_scale_ren_TSV.size() << endl;
	    for (int i_l = 0; i_l < xlist.size(); i_l++){
	      logger << LOG_DEBUG_VERBOSE << "xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b].size() = " << xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "xlist[i_l]->distribution_deviation_TSV[i_g][i_d][i_b].size() = " << xlist[i_l]->distribution_deviation_TSV[i_g][i_d][i_b].size() << endl;
	    }
	    logger << LOG_DEBUG_VERBOSE << "ygeneric->switch_output_scaleset.size() = " << ygeneric->switch_output_scaleset.size() << endl;
	    logger << LOG_DEBUG_VERBOSE << "ygeneric->switch_output_scaleset[i_s] = " << ygeneric->switch_output_scaleset[i_s] << endl;
	    logger << LOG_DEBUG_VERBOSE << "osi->switch_distribution_TSV.size() = " << osi->switch_distribution_TSV.size() << endl;
	    logger << LOG_DEBUG_VERBOSE << "osi->switch_distribution_TSV[i_s] = " << osi->switch_distribution_TSV[i_s] << endl;
	    */
	    
	    distribution_result_TSV[i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    distribution_deviation_TSV[i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));

	    /*
	    logger << LOG_DEBUG_VERBOSE << "distribution_result_TSV[i_g][i_d][i_b][i_s].size() = " << distribution_result_TSV[i_g][i_d][i_b][i_s].size() << endl;
	    logger << LOG_DEBUG_VERBOSE << "distribution_deviation_TSV[i_g][i_d][i_b][i_s].size() = " << distribution_deviation_TSV[i_g][i_d][i_b][i_s].size() << endl;
	    */
	    
	    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){logger << LOG_DEBUG_VERBOSE << "i_r = " << i_r << endl;
	      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){logger << LOG_DEBUG_VERBOSE << "i_f = " << i_f << endl;
		double dev2 = 0.;
		for (int i_l = 0; i_l < xlist.size(); i_l++){logger << LOG_DEBUG_VERBOSE << "i_l = " << i_l << endl; 
		  //		  int x_l = mapping_contribution_file[contribution_file[i_l]];  logger << LOG_DEBUG_VERBOSE << "x_l = " << x_l << endl; 
		  /*
		  logger << LOG_DEBUG_VERBOSE << "xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b][i_s].size() = " << xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b][i_s].size() << endl;
		  logger << LOG_DEBUG_VERBOSE << "xlist[i_l]->distribution_deviation_TSV[i_g][i_d][i_b][i_s].size() = " << xlist[i_l]->distribution_deviation_TSV[i_g][i_d][i_b][i_s].size() << endl;
		  */

		  distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] += xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f];
		  dev2 += pow(xlist[i_l]->distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f], 2);
		}
		//		logger << LOG_DEBUG_VERBOSE << "before distribution_deviation_TSV" << endl;
 		distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] = sqrt(dev2);
		//		logger << LOG_DEBUG_VERBOSE << "after distribution_deviation_TSV" << endl;
		/*
		stringstream temp_res;
		temp_res << setw(23) << setprecision(15) << distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f];
		stringstream temp_dev;
		temp_dev << setw(23) << setprecision(15) << distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f];
		logger << LOG_DEBUG_VERBOSE << "result_distribution_TSV[" << i_g << "][" << i_d << "][" << i_b << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
		*/
	      }
	    }
	  }
	}
      }
    }
  }
  else if (combination.size() > 1){
    logger << LOG_DEBUG_VERBOSE << "multiplicative combination started" << endl;
    vector<vector<vector<vector<vector<vector<vector<double> > > > > > > factor_order_result_TSV(combination.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (ygeneric->subgroup.size(), vector<vector<vector<vector<vector<double> > > > > (osi->extended_distribution.size())));
    vector<vector<vector<vector<vector<vector<vector<double> > > > > > > factor_order_deviation_TSV(combination.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (ygeneric->subgroup.size(), vector<vector<vector<vector<vector<double> > > > > (osi->extended_distribution.size())));
    for (int i_c = 0; i_c < combination.size(); i_c++){
      for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	  if (!ygeneric->switch_output_distribution[i_d]){continue;}
	  factor_order_result_TSV[i_c][i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
	  factor_order_deviation_TSV[i_c][i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
	  for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	    factor_order_result_TSV[i_c][i_g][i_d][i_b].resize(osi->n_extended_set_TSV);
	    factor_order_deviation_TSV[i_c][i_g][i_d][i_b].resize(osi->n_extended_set_TSV);
	    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	      if (!osi->switch_distribution_TSV[i_s]){continue;}
	      factor_order_result_TSV[i_c][i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	      factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    }
	  }
	}
      }
    }
    
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    if (!osi->switch_distribution_TSV[i_s]){continue;}
	    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){logger << LOG_DEBUG_VERBOSE << "i_r = " << i_r << endl; 
	      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){logger << LOG_DEBUG_VERBOSE << "i_f = " << i_f << endl; 
		for (int i_c = 0; i_c < combination.size(); i_c++){logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << endl;
		  double result = 0.;
		  double dev2 = 0.;
		  for (int j_c = 0; j_c < combination[i_c].size(); j_c++){
		    int i_l = combination[i_c][j_c];
		    result += xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f];
		    dev2 += pow(xlist[i_l]->distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f], 2);
		  }
		  factor_order_result_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f] = result;
		  factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f] = sqrt(dev2);
		}
	      }
	    }
	  }
	}
      }
    }
    
    /*
    // original multiplicative implementation:
    for (int i_c = 0; i_c < combination.size(); i_c++){logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << endl;
      for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
	for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
	  if (!ygeneric->switch_output_distribution[i_d]){continue;}
	  for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl; 
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    if (!osi->switch_distribution_TSV[i_s]){continue;}
	      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){logger << LOG_DEBUG_VERBOSE << "i_r = " << i_r << endl; 
		for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){logger << LOG_DEBUG_VERBOSE << "i_f = " << i_f << endl; 
		  double result = 0.;
		  double dev2 = 0.;
		  for (int i_l = 0; i_l < xlist.size(); i_l++){logger << LOG_DEBUG_VERBOSE << "i_l = " << i_l << endl; 
		    //		    int x_l = mapping_contribution_file[contribution_file[combination[i_c][i_l]]];
		    result += xlist[i_l]->distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f];
		    dev2 += pow(xlist[i_l]->distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f], 2);
		  }
		  if (i_c == 0){
		    factor_order_result_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f] = result;
		    factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f] = sqrt(dev2);
		  }
		  else {
		    if (factor_order_result_TSV[0][i_g][i_d][i_b][i_s][i_r][i_f] != 0){
		      factor_order_result_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f] = 1. + result / factor_order_result_TSV[0][i_g][i_d][i_b][i_s][i_r][i_f];
		      factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f] = sqrt(dev2 * pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_s][i_r][i_f], 2) + pow(factor_order_deviation_TSV[0][i_g][i_d][i_b][i_s][i_r][i_f] * result, 2)) / pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_s][i_r][i_f], 2);
		      // check formula !!!
		    }
		    else {
		      factor_order_result_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f] = 0.;
		      factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f] = 0.;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    */


    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    if (!osi->switch_distribution_TSV[i_s]){continue;}
	    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){logger << LOG_DEBUG_VERBOSE << "i_r = " << i_r << endl; 
	      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){logger << LOG_DEBUG_VERBOSE << "i_f = " << i_f << endl; 
		double finalized_result = 0.;
		double finalized_deviation2 = 0.;
		double index_LO = 0;
		double temp_result = 0.;
		double temp_deviation2 = 0.;
		for (int i_c = 0; i_c < combination.size(); i_c++){
		  if (combination_type[i_c] == 0){
		    if (i_c != 0){
		      finalized_result += temp_result;
		      temp_result = 0.;
		      finalized_deviation2 += temp_deviation2;
		      temp_deviation2 = 0.;
		    }
		    index_LO = i_c;
		    temp_result = factor_order_result_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f];
		    temp_deviation2 += pow(factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f], 2);
		  }
		  else if (combination_type[i_c] == 1){
		    if (factor_order_result_TSV[index_LO][i_g][i_d][i_b][i_s][i_r][i_f] == 0.){
		      temp_result += factor_order_result_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f];
		      // simplification in case of multiplicative combination:
		      temp_deviation2 += pow(factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f], 2);
		    }
		    else {
		      temp_result *= (1. + factor_order_result_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f] / factor_order_result_TSV[index_LO][i_g][i_d][i_b][i_s][i_r][i_f]);
		      // simplification in case of multiplicative combination:
		      temp_deviation2 += pow(factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f], 2);
		    }
		  }
		}
		distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] = finalized_result + temp_result;
		distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] = sqrt(finalized_deviation2 + temp_deviation2);
	      }
	    }
	  }
	}
      }
    }
      //      }
      /*
    // original multiplicative implementation:
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl; 
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  if (!osi->switch_distribution_TSV[i_s]){continue;}
	    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){logger << LOG_DEBUG_VERBOSE << "i_r = " << i_r << endl; 
	      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){logger << LOG_DEBUG_VERBOSE << "i_f = " << i_f << endl; 
		double dev2 = 0.;
		for (int i_c = 0; i_c < combination.size(); i_c++){
		  if (i_c == 0){ // could be replaced by something based on order_combination_type ... so far, one cannot simply add a contribution to the product !!!
		    distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] += factor_order_result_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f];
		  }
		  else {
		    distribution_result_TSV[i_g][i_d][i_b][i_s][i_r][i_f] *= factor_order_result_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f];
		  }
		  double temp_dev = factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_s][i_r][i_f];
		  logger << LOG_DEBUG_VERBOSE << "temp_dev[" << i_c << "][" << i_c << "] = " << temp_dev << endl;
		  for (int j_c = 0; j_c < combination.size(); j_c++){
		    if (i_c == j_c){continue;}
		    temp_dev *= factor_order_result_TSV[j_c][i_g][i_d][i_b][i_s][i_r][i_f];
		    logger << LOG_DEBUG_VERBOSE << "temp_dev[" << i_c << "][" << j_c << "] = " << temp_dev << endl;
		    // check formula !!!
		  }
		  logger << LOG_DEBUG_VERBOSE << "temp_dev[" << i_c << "] = " << temp_dev << endl;
		  dev2 += pow(temp_dev, 2);
		}
		distribution_deviation_TSV[i_g][i_d][i_b][i_s][i_r][i_f] = sqrt(dev2);
	      }
	    }
	  }
	}
      }
      */
  
    logger << LOG_DEBUG_VERBOSE << "multiplicative combination finished" << endl;
  }



  /////////////////////////////////////////////////////////////////////
  //  calculate distribution_result/deviation_qTcut_TSV  //
  /////////////////////////////////////////////////////////////////////

  //  for (int i_o = 0; i_o < yorder.size(); i_o++){logger << LOG_DEBUG_VERBOSE << "i_o = " << i_o << endl;
  if (combination.size() == 1){
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	  for (int i_q = 0; i_q < selection_n_qTcut; i_q++){logger << LOG_DEBUG_VERBOSE << "i_q = " << i_q << endl;
	    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
	      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	      if (!osi->switch_distribution_TSV[i_s]){continue;}
	      distribution_result_qTcut_TSV[i_q][i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	      distribution_deviation_qTcut_TSV[i_q][i_g][i_d][i_b][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){logger << LOG_DEBUG_VERBOSE << "i_r = " << i_r << endl;
		for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){logger << LOG_DEBUG_VERBOSE << "i_f = " << i_f << endl;
		  double dev = 0.;
		  for (int i_l = 0; i_l < xlist.size(); i_l++){logger << LOG_DEBUG_VERBOSE << "i_l = " << i_l << endl; 
		    //		    int x_l = mapping_contribution_file[contribution_file[i_l]];  logger << LOG_DEBUG_VERBOSE << "x_l = " << x_l << endl; 
		    int y_q = 0;
		    if (xlist[i_l]->xcontribution[0].active_qTcut){y_q = i_q;}
		    logger << LOG_DEBUG_VERBOSE << "y_q = " << y_q << endl; 
		    //		    logger << LOG_DEBUG_VERBOSE << "distribution_result_TSV[" << i_g << "][" << i_d << "][" << i_b << "].size() = " << distribution_result_qTcut_TSV[i_g][i_d][i_b].size() << endl;
		    logger << LOG_DEBUG_VERBOSE << "*xlist[" << i_l << "].xcontribution[0].distribution_result_TSV[" << i_g << "][" << i_d << "][" << i_b << "].size() = " << xlist[i_l]->xcontribution[0].distribution_result_TSV[i_g][i_d][i_b].size() << endl;
		    logger << LOG_DEBUG_VERBOSE << "*xlist[" << i_l << "].xcontribution[0].distribution_deviation_TSV[" << i_g << "][" << i_d << "][" << i_b << "].size() = " << xlist[i_l]->xcontribution[0].distribution_deviation_TSV[i_g][i_d][i_b].size() << endl;
		    distribution_result_qTcut_TSV[i_q][i_g][i_d][i_b][i_s][i_r][i_f] += xlist[i_l]->xcontribution[0].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f];
		    dev += pow(xlist[i_l]->xcontribution[0].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f], 2);
		  }
		  distribution_deviation_qTcut_TSV[i_q][i_g][i_d][i_b][i_s][i_r][i_f] = sqrt(dev);
		  stringstream temp_res;
		  temp_res << setw(23) << setprecision(15) << distribution_result_qTcut_TSV[i_q][i_g][i_d][i_b][i_s][i_r][i_f];
		  stringstream temp_dev;
		  temp_dev << setw(23) << setprecision(15) << distribution_deviation_qTcut_TSV[i_q][i_g][i_d][i_b][i_s][i_r][i_f];
		  logger << LOG_DEBUG_VERBOSE << "XXX   result_TSV[" << i_g << "][" << i_d << "][" << i_b << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  else if (combination.size() > 1){
    logger << LOG_DEBUG << "multiplicative combination started" << endl;

    logger << LOG_DEBUG << "combination_type.size() = " << combination_type.size() << endl;
    for (int i_c = 0; i_c < combination_type.size(); i_c++){
      logger << LOG_DEBUG << "combination_type[" << i_c << "] = " << combination_type[i_c] << endl;
    }

    logger << LOG_DEBUG << "combination.size() = " << combination.size() << endl;
    for (int i_c = 0; i_c < combination.size(); i_c++){
      logger << LOG_DEBUG << "combination[" << i_c << "].size() = " << combination[i_c].size() << endl;
      for (int j_c = 0; j_c < combination[i_c].size(); j_c++){
	logger << LOG_DEBUG << "combination[" << i_c << "][" << j_c << "] = " << combination[i_c][j_c] << endl;
      }
    }

  

    
    vector<vector<vector<vector<vector<vector<vector<vector<double> > > > > > > > factor_order_result_TSV(combination.size(), vector<vector<vector<vector<vector<vector<vector<double> > > > > > > (ygeneric->subgroup.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->extended_distribution.size())));
    vector<vector<vector<vector<vector<vector<vector<vector<double> > > > > > > > factor_order_deviation_TSV(combination.size(), vector<vector<vector<vector<vector<vector<vector<double> > > > > > > (ygeneric->subgroup.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->extended_distribution.size())));
    for (int i_c = 0; i_c < combination.size(); i_c++){
      factor_order_result_TSV[i_c].resize(ygeneric->subgroup.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->extended_distribution.size()));
      factor_order_deviation_TSV[i_c].resize(ygeneric->subgroup.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->extended_distribution.size()));
      for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){
	for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){
	  if (!ygeneric->switch_output_distribution[i_d]){continue;}
	  factor_order_result_TSV[i_c][i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
	  factor_order_deviation_TSV[i_c][i_g][i_d].resize(osi->extended_distribution[i_d].n_bins);
	  for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){
	    factor_order_result_TSV[i_c][i_g][i_d][i_b].resize(selection_n_qTcut);
	    factor_order_deviation_TSV[i_c][i_g][i_d][i_b].resize(selection_n_qTcut);
	    for (int i_q = 0; i_q < selection_n_qTcut; i_q++){
	      factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q].resize(osi->n_extended_set_TSV);
	      factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q].resize(osi->n_extended_set_TSV);
	      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
		if (!ygeneric->switch_output_scaleset[i_s]){continue;}
		if (!osi->switch_distribution_TSV[i_s]){continue;}
		factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
		factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	      }
	    }
	  }
	}
      }
    }
    
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	  for (int i_q = 0; i_q < selection_n_qTcut; i_q++){logger << LOG_DEBUG_VERBOSE << "i_q = " << i_q << endl;
	    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
	      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	      if (!osi->switch_distribution_TSV[i_s]){continue;}
	      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
		for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
		  
		  for (int i_c = 0; i_c < combination.size(); i_c++){

		    double result = 0.;
		    double dev2 = 0.;
		    
		    for (int j_c = 0; j_c < combination[i_c].size(); j_c++){
		      int i_l = combination[i_c][j_c];
			//	for (int i_l = 0; i_l < xlist.size(); i_l++){
		      //		    for (int i_l = 0; i_l < order_contribution_file[i_c].size(); i_l++){
		      //		      int x_l = mapping_contribution_file[contribution_file[combination[i_c][i_l]]];
		      //		      logger << LOG_DEBUG_VERBOSE << "XYZ2: mapping_contribution_file[order_contribution_file[" << i_o << "][order_combination[" << i_o << "][" << i_c << "][" << i_l << "] = " << combination[i_c][i_l] << "] = " << contribution_file[combination[i_c][i_l]] << "] = " << x_l << endl;;
		      //		      logger << LOG_INFO << "order_contribution_file[" << i_c << "][" << i_l << "] = " << order_contribution_file[i_c][i_l] << endl;
		      //		      logger << LOG_DEBUG << "list_result_TSV[" << x_l << "][" << i_m << "][" << i_g << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << *xlist[x_l].xcontribution[0].distribution_result_TSV[i_g][i_d][i_b][i_q][i_s][i_r][i_f] << endl;
		      
		      int y_q = 0;
		      if (xlist[i_l]->xcontribution[0].active_qTcut){y_q = i_q;}
		      
		      result += xlist[i_l]->xcontribution[0].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f];
		      dev2 += pow(xlist[i_l]->xcontribution[0].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f], 2);
		    }
		    factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = result;
		    factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = sqrt(dev2);


		    /*
		    if (i_c == 0){
		      factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = result;
		      factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = sqrt(dev2);
		      
		    }
		    else {
		      if (factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] != 0){
			factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = 1. + result / factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f];
			factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = sqrt(dev2 * pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f], 2) + pow(factor_order_deviation_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] * result, 2)) / pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f], 2);
		      }
		      else {
			factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = 0.;
			factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = 0.;
		      }			//		      factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = (sqrt(dev2) * factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] + factor_order_deviation_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] * abs(result)) / pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f], 2);
		    }

		    */

		    
			/*
			  logger << LOG_DEBUG << "sqrt(dev) = " << sqrt(dev2) << endl;
			  logger << LOG_DEBUG << "factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = " << factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] << endl;
			  logger << LOG_DEBUG << "factor_order_deviation_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = " << factor_order_deviation_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] << endl;
			  logger << LOG_DEBUG << "abs(result) = " << abs(result) << endl;
			  logger << LOG_DEBUG << "1. / pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f], 2) = " << 1. / pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f], 2) << endl;
			*/
		    /*
		      stringstream temp_res;
		      temp_res << setw(23) << setprecision(15) << factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f];
		      stringstream temp_dev;
		      temp_dev << setw(23) << setprecision(15) << factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f];
		      logger << LOG_DEBUG_VERBOSE << "XXX   factor_order_result_TSV[" << i_c << "][" << i_m << "][" << i_g << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
		    */
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    /*
// original multiplicative implementation:
for (int i_c = 0; i_c < combination.size(); i_c++){
      for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
	for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
	  if (!ygeneric->switch_output_distribution[i_d]){continue;}
	  for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	    for (int i_q = 0; i_q < selection_n_qTcut; i_q++){logger << LOG_DEBUG_VERBOSE << "i_q = " << i_q << endl;
	      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl; 
	      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	      if (!osi->switch_distribution_TSV[i_s]){continue;}
		for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
		  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
		    double result = 0.;
		    double dev2 = 0.;
		    for (int i_l = 0; i_l < xlist.size(); i_l++){
		      //		    for (int i_l = 0; i_l < order_contribution_file[i_c].size(); i_l++){
		      //		      int x_l = mapping_contribution_file[contribution_file[combination[i_c][i_l]]];
		      //		      logger << LOG_DEBUG_VERBOSE << "XYZ2: mapping_contribution_file[order_contribution_file[" << i_o << "][order_combination[" << i_o << "][" << i_c << "][" << i_l << "] = " << combination[i_c][i_l] << "] = " << contribution_file[combination[i_c][i_l]] << "] = " << x_l << endl;;
		      //		      logger << LOG_INFO << "order_contribution_file[" << i_c << "][" << i_l << "] = " << order_contribution_file[i_c][i_l] << endl;
		      //		      logger << LOG_DEBUG << "list_result_TSV[" << x_l << "][" << i_m << "][" << i_g << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << *xlist[x_l].xcontribution[0].distribution_result_TSV[i_g][i_d][i_b][i_q][i_s][i_r][i_f] << endl;
		      
		      int y_q = 0;
		      if (xlist[i_l]->xcontribution[0].active_qTcut){y_q = i_q;}
		      
		      result += xlist[i_l]->xcontribution[0].distribution_result_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f];
		      dev2 += pow(xlist[i_l]->xcontribution[0].distribution_deviation_TSV[i_g][i_d][i_b][y_q][i_s][i_r][i_f], 2);
		    }
		    if (i_c == 0){
		      factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = result;
		      factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = sqrt(dev2);
		      
		    }
		    else {
		      if (factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] != 0){
			factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = 1. + result / factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f];
			/*
			  logger << LOG_DEBUG << "sqrt(dev) = " << sqrt(dev2) << endl;
			  logger << LOG_DEBUG << "factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = " << factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] << endl;
			  logger << LOG_DEBUG << "factor_order_deviation_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = " << factor_order_deviation_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] << endl;
			  logger << LOG_DEBUG << "abs(result) = " << abs(result) << endl;
			  logger << LOG_DEBUG << "1. / pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f], 2) = " << 1. / pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f], 2) << endl;
    *//*
			factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = sqrt(dev2 * pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f], 2) + pow(factor_order_deviation_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] * result, 2)) / pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f], 2);
		      }
		      else {
			factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = 0.;
			factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = 0.;
		      }			//		      factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f] = (sqrt(dev2) * factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] + factor_order_deviation_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f] * abs(result)) / pow(factor_order_result_TSV[0][i_g][i_d][i_b][i_q][i_s][i_r][i_f], 2);
		    }
		    /*
		      stringstream temp_res;
		      temp_res << setw(23) << setprecision(15) << factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f];
		      stringstream temp_dev;
		      temp_dev << setw(23) << setprecision(15) << factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f];
		      logger << LOG_DEBUG_VERBOSE << "XXX   factor_order_result_TSV[" << i_c << "][" << i_m << "][" << i_g << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
      *//*
		  }
		}
	      }
	    }
	  }
	}
      }
    }
	*/
    
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
      for (int i_d = 0; i_d < osi->extended_distribution.size(); i_d++){logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << endl;
	if (!ygeneric->switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < osi->extended_distribution[i_d].n_bins; i_b++){logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << endl;
	  for (int i_q = 0; i_q < selection_n_qTcut; i_q++){logger << LOG_DEBUG_VERBOSE << "i_q = " << i_q << endl;
	    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
	      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	      if (!osi->switch_distribution_TSV[i_s]){continue;}
	      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
		for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
		  double dev2 = 0.;

		  
		  for (int i_c = 0; i_c < combination.size(); i_c++){
		    if (i_c == 0){ // could be replaced by something based on order_combination_type ... so far, one cannot simply add a contribution to the product !!!
		      distribution_result_qTcut_TSV[i_q][i_g][i_d][i_b][i_s][i_r][i_f] += factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f];
		    }
		    else {
		      distribution_result_qTcut_TSV[i_q][i_g][i_d][i_b][i_s][i_r][i_f] *= factor_order_result_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f];
		    }
		    double temp_dev = factor_order_deviation_TSV[i_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f];
		    logger << LOG_DEBUG_VERBOSE << "temp_dev[" << i_c << "][" << i_c << "] = " << temp_dev << endl;
		    for (int j_c = 0; j_c < combination.size(); j_c++){
		      if (i_c == j_c){continue;}
		      temp_dev *= factor_order_result_TSV[j_c][i_g][i_d][i_b][i_q][i_s][i_r][i_f];
		      logger << LOG_DEBUG_VERBOSE << "temp_dev[" << i_c << "][" << j_c << "] = " << temp_dev << endl;
		    }
		    logger << LOG_DEBUG_VERBOSE << "temp_dev[" << i_c << "] = " << temp_dev << endl;
		    //		    distribution_deviation_qTcut_TSV[i_q][i_g][i_d][i_b][i_s][i_r][i_f] += temp_dev;
		    dev2 += pow(temp_dev, 2);
		  }
		  distribution_deviation_qTcut_TSV[i_q][i_g][i_d][i_b][i_s][i_r][i_f] = sqrt(dev2);
		}
	      }
	    }
	  }
	}
      }
    }
      
    logger << LOG_DEBUG_VERBOSE << "multiplicative combination finished" << endl;
  }




  if (ygeneric->switch_output_overview > 0){output_distribution_overview_TSV();}
  if (ygeneric->switch_output_plot > 0){output_distribution_qTcut_TSV();}
  if (ygeneric->switch_output_plot > 0){output_distribution_TSV();}

  
  logger << LOG_DEBUG << "finished" << endl;
}



