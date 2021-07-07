#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_order::collect_result_TSV(vector<string> & subgroup){
  Logger logger("summary_order::collect_result_TSV");
  logger << LOG_DEBUG << "called" << endl;

  result_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<double> > > > (subgroup.size(), vector<vector<vector<double> > > (osi->n_extended_set_TSV)));
  deviation_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<double> > > > (subgroup.size(), vector<vector<vector<double> > > (osi->n_extended_set_TSV)));
  
  result_qTcut_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<vector<double> > > > > (subgroup.size(), vector<vector<vector<vector<double> > > > (output_n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV))));
  deviation_qTcut_TSV.resize(osi->n_moments + 1, vector<vector<vector<vector<vector<double> > > > > (subgroup.size(), vector<vector<vector<vector<double> > > > (output_n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV))));

  logger << LOG_DEBUG << "multiplicative decision: combination.size() = " << combination.size() << endl;
  if (combination.size() == 1){
    logger << LOG_DEBUG << "result_qTcut_TSV.size() = " << result_qTcut_TSV.size() << endl;
    for (int i_g = 0; i_g < subgroup.size(); i_g++){
      logger << LOG_DEBUG << "i_g = " << i_g << endl;
      for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
	logger << LOG_DEBUG << "i_m = " << i_m << endl;
	if (i_m > 0){continue;} // no moments implemented so far !!!
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  result_TSV[i_m][i_g][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  deviation_TSV[i_m][i_g][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      double dev = 0.;
	      for (int i_l = 0; i_l < xlist.size(); i_l++){
		//		  int x_l = mapping_contribution_file[contribution_file[i_l]];
		// later via pointer !!!
		result_TSV[i_m][i_g][i_s][i_r][i_f] += xlist[i_l]->result_TSV[i_m][i_g][i_s][i_r][i_f];
		dev += pow(xlist[i_l]->deviation_TSV[i_m][i_g][i_s][i_r][i_f], 2);
	      }
	      deviation_TSV[i_m][i_g][i_s][i_r][i_f] = sqrt(dev);
	      stringstream temp_res;
	      temp_res << setw(23) << setprecision(15) << result_TSV[i_m][i_g][i_s][i_r][i_f];
	      stringstream temp_dev;
	      temp_dev << setw(23) << setprecision(15) << deviation_TSV[i_m][i_g][i_s][i_r][i_f];
	      logger << LOG_DEBUG << "result_TSV[" << i_m << "][" << i_g << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
	    }
	  }
	}
	
	for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    result_qTcut_TSV[i_m][i_g][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    deviation_qTcut_TSV[i_m][i_g][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
		double dev2 = 0.;
		for (int i_l = 0; i_l < xlist.size(); i_l++){
		  //		    int x_l = mapping_contribution_file[contribution_file[i_l]];
		  //		    logger << LOG_DEBUG_VERBOSE << "mapping_contribution_file[yorder[" << i_o << "].contribution_file[" << i_l << "] = " << contribution_file[i_l] << "] = " << x_l << endl;
		  int y_q = 0;
		  if (xlist[i_l]->xcontribution[0].active_qTcut){y_q = i_q;}
		  
		  result_qTcut_TSV[i_m][i_g][i_q][i_s][i_r][i_f] += xlist[i_l]->xcontribution[0].result_TSV[i_m][i_g][y_q][i_s][i_r][i_f];
		  dev2 += pow(xlist[i_l]->xcontribution[0].deviation_TSV[i_m][i_g][y_q][i_s][i_r][i_f], 2);
		}
		deviation_qTcut_TSV[i_m][i_g][i_q][i_s][i_r][i_f] = sqrt(dev2);
		stringstream temp_res;
		temp_res << setw(23) << setprecision(15) << result_qTcut_TSV[i_m][i_g][i_q][i_s][i_r][i_f];
		stringstream temp_dev;
		temp_dev << setw(23) << setprecision(15) << deviation_qTcut_TSV[i_m][i_g][i_q][i_s][i_r][i_f];
		logger << LOG_DEBUG_VERBOSE << "result_qTcut_TSV[" << i_m << "][" << i_g << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
	      }
	    }
	  }
	}
      }
    }
  }
  else if (combination.size() > 1){
    logger << LOG_DEBUG << "multiplicative combination started" << endl;
    vector<vector<vector<vector<vector<vector<double> > > > > > factor_order_result_TSV(combination.size(), vector<vector<vector<vector<vector<double> > > > > (osi->n_moments + 1, vector<vector<vector<vector<double> > > > (subgroup.size(), vector<vector<vector<double> > > (osi->n_extended_set_TSV))));
    vector<vector<vector<vector<vector<vector<double> > > > > > factor_order_deviation_TSV(combination.size(), vector<vector<vector<vector<vector<double> > > > > (osi->n_moments + 1, vector<vector<vector<vector<double> > > > (subgroup.size(), vector<vector<vector<double> > > (osi->n_extended_set_TSV))));
    
    vector<vector<vector<vector<vector<vector<vector<double> > > > > > > factor_order_result_qTcut_TSV(combination.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->n_moments + 1, vector<vector<vector<vector<vector<double> > > > > (subgroup.size(), vector<vector<vector<vector<double> > > > (output_n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV)))));
    vector<vector<vector<vector<vector<vector<vector<double> > > > > > > factor_order_deviation_qTcut_TSV(combination.size(), vector<vector<vector<vector<vector<vector<double> > > > > > (osi->n_moments + 1, vector<vector<vector<vector<vector<double> > > > > (subgroup.size(), vector<vector<vector<vector<double> > > > (output_n_qTcut, vector<vector<vector<double> > > (osi->n_extended_set_TSV)))));
    for (int i_c = 0; i_c < combination.size(); i_c++){
      for (int i_g = 0; i_g < subgroup.size(); i_g++){
	for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
	  if (i_m > 0){continue;} // no moments implemented so far !!!
	  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    factor_order_result_TSV[i_c][i_m][i_g][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    factor_order_deviation_TSV[i_c][i_m][i_g][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
		double result = 0.;
		double dev2 = 0.;
		for (int i_l = 0; i_l < xlist.size(); i_l++){
		  //		    int x_l = mapping_contribution_file[contribution_file[combination[i_c][i_l]]];
		  result += xlist[i_l]->result_TSV[i_m][i_g][i_s][i_r][i_f];
		  dev2 += pow(xlist[i_l]->deviation_TSV[i_m][i_g][i_s][i_r][i_f], 2);
		}
		if (i_c == 0){
		  factor_order_result_TSV[i_c][i_m][i_g][i_s][i_r][i_f] = result;
		  factor_order_deviation_TSV[i_c][i_m][i_g][i_s][i_r][i_f] = sqrt(dev2);
		  
		}
		else {
		  factor_order_result_TSV[i_c][i_m][i_g][i_s][i_r][i_f] = 1. + result / factor_order_result_TSV[0][i_m][i_g][i_s][i_r][i_f];
		  factor_order_deviation_TSV[i_c][i_m][i_g][i_s][i_r][i_f] = sqrt(dev2 * pow(factor_order_result_TSV[0][i_m][i_g][i_s][i_r][i_f], 2) + pow(factor_order_deviation_TSV[0][i_m][i_g][i_s][i_r][i_f] * result, 2)) / pow(factor_order_result_TSV[0][i_m][i_g][i_s][i_r][i_f], 2);
		}
		stringstream temp_res;
		temp_res << setw(23) << setprecision(15) << factor_order_result_TSV[i_c][i_m][i_g][i_s][i_r][i_f];
		stringstream temp_dev;
		temp_dev << setw(23) << setprecision(15) << factor_order_deviation_TSV[i_c][i_m][i_g][i_s][i_r][i_f];
		logger << LOG_INFO << "XXXM   factor_order_result_TSV[" << i_c << "][" << i_m << "][" << i_g << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
	      }
	    }
	  }
	  
	  for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	      factor_order_result_qTcut_TSV[i_c][i_m][i_g][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	      factor_order_deviation_qTcut_TSV[i_c][i_m][i_g][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
	      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
		for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
		  double result = 0.;
		  double dev2 = 0.;
		  for (int i_l = 0; i_l < xlist.size(); i_l++){
		    //		    for (int i_l = 0; i_l < order_contribution_file[i_c].size(); i_l++){
		    //		      int x_l = mapping_contribution_file[contribution_file[combination[i_c][i_l]]];
		    //		      logger << LOG_DEBUG << "mapping_contribution_file[yorder[" << i_o << "].contribution_file[yorder[" << i_o << "].combination[" << i_c << "][" << i_l << "] = " << combination[i_c][i_l] << "] = " << contribution_file[combination[i_c][i_l]] << "] = " << x_l << endl;;
		    //		      logger << LOG_INFO << "order_contribution_file[" << i_c << "][" << i_l << "] = " << order_contribution_file[i_c][i_l] << endl;
		    logger << LOG_DEBUG << "list_result_TSV[" << i_l << "][" << i_m << "][" << i_g << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << xlist[i_l]->xcontribution[0].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f] << endl;
		    result += xlist[i_l]->xcontribution[0].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f];
		    dev2 += pow(xlist[i_l]->xcontribution[0].deviation_TSV[i_m][i_g][i_q][i_s][i_r][i_f], 2);
		  }
		  if (i_c == 0){
		    factor_order_result_qTcut_TSV[i_c][i_m][i_g][i_q][i_s][i_r][i_f] = result;
		    factor_order_deviation_qTcut_TSV[i_c][i_m][i_g][i_q][i_s][i_r][i_f] = sqrt(dev2);
		    
		  }
		  else {
		    factor_order_result_qTcut_TSV[i_c][i_m][i_g][i_q][i_s][i_r][i_f] = 1. + result / factor_order_result_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f];
		    /*
		      logger << LOG_DEBUG << "sqrt(dev) = " << sqrt(dev2) << endl;
		      logger << LOG_DEBUG << "factor_order_result_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f] = " << factor_order_result_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f] << endl;
		      logger << LOG_DEBUG << "factor_order_deviation_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f] = " << factor_order_deviation_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f] << endl;
		      logger << LOG_DEBUG << "abs(result) = " << abs(result) << endl;
		      logger << LOG_DEBUG << "1. / pow(factor_order_result_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f], 2) = " << 1. / pow(factor_order_result_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f], 2) << endl;
		    */
		    factor_order_deviation_qTcut_TSV[i_c][i_m][i_g][i_q][i_s][i_r][i_f] = sqrt(dev2 * pow(factor_order_result_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f], 2) + pow(factor_order_deviation_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f] * result, 2)) / pow(factor_order_result_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f], 2);
		    //		      factor_order_deviation_qTcut_TSV[i_c][i_m][i_g][i_q][i_s][i_r][i_f] = (sqrt(dev2) * factor_order_result_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f] + factor_order_deviation_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f] * abs(result)) / pow(factor_order_result_qTcut_TSV[0][i_m][i_g][i_q][i_s][i_r][i_f], 2);
		  }
		  
		  stringstream temp_res;
		  temp_res << setw(23) << setprecision(15) << factor_order_result_qTcut_TSV[i_c][i_m][i_g][i_q][i_s][i_r][i_f];
		  stringstream temp_dev;
		  temp_dev << setw(23) << setprecision(15) << factor_order_deviation_qTcut_TSV[i_c][i_m][i_g][i_q][i_s][i_r][i_f];
		  logger << LOG_DEBUG_VERBOSE << "XXXM   factor_order_result_qTcut_TSV[" << i_c << "][" << i_m << "][" << i_g << "][" << i_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
		}
	      }
	    }
	  }
	}
      }
    }
    
    for (int i_g = 0; i_g < subgroup.size(); i_g++){
      logger << LOG_DEBUG << "i_g = " << i_g << endl;
      for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
	logger << LOG_DEBUG << "i_m = " << i_m << endl;
	if (i_m > 0){continue;} // no moments implemented so far !!!
	
	for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	  if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	  logger << LOG_DEBUG << "i_s = " << i_s << endl; 
	  result_TSV[i_m][i_g][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s]));
	  deviation_TSV[i_m][i_g][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s]));
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      double dev2 = 0.;
	      for (int i_c = 0; i_c < combination.size(); i_c++){
		if (i_c == 0){ // could be replaced by something based on order_combination_type ... so far, one cannot simply add a contribution to the product !!!
		  result_TSV[i_m][i_g][i_s][i_r][i_f] += factor_order_result_TSV[i_c][i_m][i_g][i_s][i_r][i_f];
		}
		else {
		  result_TSV[i_m][i_g][i_s][i_r][i_f] *= factor_order_result_TSV[i_c][i_m][i_g][i_s][i_r][i_f];
		}
		double temp_dev = factor_order_deviation_TSV[i_c][i_m][i_g][i_s][i_r][i_f];
		for (int j_c = 0; j_c < combination.size(); j_c++){
		  if (i_c == j_c){continue;}
		  temp_dev *= factor_order_result_TSV[j_c][i_m][i_g][i_s][i_r][i_f];
		}
		dev2 += pow(temp_dev, 2);
	      }
	      deviation_TSV[i_m][i_g][i_s][i_r][i_f] = sqrt(dev2);
	    }
	  }
	}
	
	for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	  logger << LOG_DEBUG << "i_q = " << i_q << endl;
	  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	    logger << LOG_DEBUG << "i_s = " << i_s << endl; 
	    result_qTcut_TSV[i_m][i_g][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s]));
	    deviation_qTcut_TSV[i_m][i_g][i_q][i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s]));
	    for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
		double dev2 = 0.;
		for (int i_c = 0; i_c < combination.size(); i_c++){
		  if (i_c == 0){ // could be replaced by something based on order_combination_type ... so far, one cannot simply add a contribution to the product !!!
		    result_qTcut_TSV[i_m][i_g][i_q][i_s][i_r][i_f] += factor_order_result_qTcut_TSV[i_c][i_m][i_g][i_q][i_s][i_r][i_f];
		  }
		  else {
		    result_qTcut_TSV[i_m][i_g][i_q][i_s][i_r][i_f] *= factor_order_result_qTcut_TSV[i_c][i_m][i_g][i_q][i_s][i_r][i_f];
		  }
		  double temp_dev = factor_order_deviation_qTcut_TSV[i_c][i_m][i_g][i_q][i_s][i_r][i_f];
		  cout << "temp_dev[" << i_c << "][" << i_c << "] = " << temp_dev << endl;
		  for (int j_c = 0; j_c < combination.size(); j_c++){
		    if (i_c == j_c){continue;}
		    temp_dev *= factor_order_result_qTcut_TSV[j_c][i_m][i_g][i_q][i_s][i_r][i_f];
		    cout << "temp_dev[" << i_c << "][" << j_c << "] = " << temp_dev << endl;
		  }
		  cout << "temp_dev[" << i_c << "] = " << temp_dev << endl;
		  //		    deviation_qTcut_TSV[i_m][i_g][i_q][i_s][i_r][i_f] += temp_dev;
		  dev2 += pow(temp_dev, 2);
		}
		deviation_qTcut_TSV[i_m][i_g][i_q][i_s][i_r][i_f] = sqrt(dev2);
		
	      }
	    }
	  }
	}
      }
    }
    logger << LOG_DEBUG << "multiplicative combination finished" << endl;
  }
  logger << LOG_DEBUG << "additive/multiplicative combination finished" << endl;
  
  if (ygeneric->switch_output_overview > 0){output_result_overview_TSV();}
  if (ygeneric->switch_output_result > 0){output_result_TSV();}
  if (ygeneric->switch_output_plot > 0){output_result_plot_TSV();}
  if (ygeneric->switch_output_plot > 0){output_result_plot_qTcut_TSV();}
  
  logger << LOG_DEBUG << "finished" << endl;
}



void summary_order::collect_result_CV(vector<string> & subgroup){
  Logger logger("summary_order::collect_result_CV");
  logger << LOG_DEBUG << "called" << endl;

  // extrapolated result
  result_CV.resize(osi->n_moments + 1, vector<vector<double> > (subgroup.size(), vector<double> (osi->n_scales_CV, 0.)));
  deviation_CV.resize(osi->n_moments + 1, vector<vector<double> > (subgroup.size(), vector<double> (osi->n_scales_CV, 0.)));

  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      for (int i_g = 0; i_g < subgroup.size(); i_g++){
	double dev = 0.;
	/*
	for (int i_l = 0; i_l < contribution_file.size(); i_l++){
	  int x_l = mapping_contribution_file[contribution_file[i_l]];
	  result_CV[i_m][i_g][i_s] += xlist[x_l].result_CV[i_m][i_g][i_s];
	  dev += pow(xlist[x_l].deviation_CV[i_m][i_g][i_s], 2);
	}
	*/
	for (int i_l = 0; i_l < xlist.size(); i_l++){
	  result_CV[i_m][i_g][i_s] += xlist[i_l]->result_CV[i_m][i_g][i_s];
	  dev += pow(xlist[i_l]->deviation_CV[i_m][i_g][i_s], 2);
	}
	deviation_CV[i_m][i_g][i_s] = sqrt(dev);
	stringstream temp_res;
	temp_res << setw(23) << setprecision(15) << result_CV[i_m][i_g][i_s];
	stringstream temp_dev;
	temp_dev << setw(23) << setprecision(15) << deviation_CV[i_m][i_g][i_s];
	logger << LOG_DEBUG_VERBOSE << "result_CV[" << i_m << "][" << i_g << "][" << i_s << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
      }
    }
  }

  // qTcut result
  result_qTcut_CV.resize(osi->n_moments + 1, vector<vector<vector<double> > > (subgroup.size(), vector<vector<double> > (osi->n_qTcut, vector<double> (osi->n_scales_CV, 0.))));
  deviation_qTcut_CV.resize(osi->n_moments + 1, vector<vector<vector<double> > > (subgroup.size(), vector<vector<double> > (osi->n_qTcut, vector<double> (osi->n_scales_CV, 0.))));

  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	for (int i_g = 0; i_g < subgroup.size(); i_g++){
	  double dev = 0.;
	  for (int i_l = 0; i_l < xlist.size(); i_l++){
	    //	    result_qTcut_CV[i_m][i_g][i_q][i_s] += xlist[i_l]->result_qTcut_CV[i_m][i_g][i_q][i_s];
	    //	    dev += pow(xlist[i_l]->deviation_qTcut_CV[i_m][i_g][i_q][i_s], 2);
	    result_qTcut_CV[i_m][i_g][i_q][i_s] += xlist[i_l]->xcontribution[0].result_CV[i_m][i_g][i_q][i_s];
	    dev += pow(xlist[i_l]->xcontribution[0].deviation_CV[i_m][i_g][i_q][i_s], 2);
	  }
	  /*
	  for (int i_l = 0; i_l < contribution_file.size(); i_l++){
	    int x_l = mapping_contribution_file[contribution_file[i_l]];
	    result_qTcut_CV[i_m][i_g][i_q][i_s] += xlist[x_l].xcontribution[0].result_qTcut_CV[i_m][i_g][i_q][i_s];
	    dev += pow(xlist[x_l].xcontribution[0].deviation_qTcut_CV[i_m][i_g][i_q][i_s], 2);
	  }
	  */
	  deviation_qTcut_CV[i_m][i_g][i_q][i_s] = sqrt(dev);
	  stringstream temp_res;
	  temp_res << setw(23) << setprecision(15) << result_qTcut_CV[i_m][i_g][i_q][i_s];
	  stringstream temp_dev;
	  temp_dev << setw(23) << setprecision(15) << deviation_qTcut_CV[i_m][i_g][i_q][i_s];
	  logger << LOG_DEBUG_VERBOSE << "result_qTcut_CV[" << i_m << "][" << i_g << "][" << i_q << "][" << i_s << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
	}
      }
    }
  }
  
  if (ygeneric->switch_output_overview > 0){output_result_overview_CV();}
  if (ygeneric->switch_output_result > 0){output_result_CV();}
  if (ygeneric->switch_output_plot > 0){output_result_plot_CV();}
  if (ygeneric->switch_output_plot > 0){output_result_plot_qTcut_CV();}
  
  logger << LOG_DEBUG << "finished" << endl;
}

