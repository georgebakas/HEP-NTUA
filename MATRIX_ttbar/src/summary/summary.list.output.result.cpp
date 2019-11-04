#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_list::output_contribution_result_overview_qTcut_TSV(){
  Logger logger("summary_list::output_contribution_result_overview_qTcut_TSV");
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

  stringstream temp_ss_separation_one;
  stringstream temp_ss_separation_two;
  int counter_separation = 15 + 15 + 15 + 25 + 15 + 4 + 15;
  for (int i_x = 0; i_x < counter_separation; i_x++){
    temp_ss_separation_one << "-";
    temp_ss_separation_two << "=";
  }
  
  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/complete/" + resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";
      filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/ren/" + resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";
      filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/fact/" + resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";
      filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/equal/" + resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";
      filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/antipodal/" + resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";

      logger << LOG_INFO << "FILENAME (complete):   " << filename_complete << endl;

      outfile_complete.open(filename_complete.c_str(), ofstream::out | ofstream::trunc);  
      outfile_ren.open(filename_ren.c_str(), ofstream::out | ofstream::trunc);  
      outfile_fact.open(filename_fact.c_str(), ofstream::out | ofstream::trunc);  
      outfile_equal.open(filename_equal.c_str(), ofstream::out | ofstream::trunc);  
      outfile_antipodal.open(filename_antipodal.c_str(), ofstream::out | ofstream::trunc);  

      stringstream name_ss;
      name_ss << "dXS +- err (contribution)";

      stringstream header_ss;
      header_ss << left
		<< setw(15) << "qTcut"
		<< setw(5) << "mu_R/" << setw(10) << osi->name_extended_set_TSV[i_s] 
		<< setw(5) << "mu_F/" << setw(10) << osi->name_extended_set_TSV[i_s] 
		<< setw(25) << "contribution"
		<< setw(34) << right << name_ss.str();
      
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
	  stringstream temp_ss;
	  temp_ss << left;
	  if (active_qTcut){temp_ss << setw(15) << "extrapolation";}
	  else {temp_ss << setw(15) << "independent";}
	  temp_ss << setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r] 
		  << setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f];
	  temp_ss << setw(25) << xcontribution[0].infix_order_contribution
		  << right << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_s][i_r][i_f] 
		  << setw(4) << " +- "
		  << right << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][i_r][i_f];
	  
	  outfile_complete << temp_ss.str() << endl;
	  if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str() << endl;}
	  if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str() << endl;}
	  if (i_r == i_f){outfile_equal << temp_ss.str() << endl;}
	  if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str() << endl;}
	  
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

      for (int i_q = 0; i_q < xcontribution[0].output_n_qTcut; i_q++){
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
	    int temp_size_result = int(log10(abs(osi->unit_factor_result * xcontribution[0].result_TSV[i_m][0][i_q][i_s][i_r][i_f])));
	    for (int i_c = 0; i_c < xcontribution.size(); i_c++){
	      int temp_size_result_c = int(log10(abs(osi->unit_factor_result * xcontribution[i_c].result_TSV[i_m][0][i_q][i_s][i_r][i_f])));
	      if (temp_size_result_c > temp_size_result){temp_size_result = temp_size_result_c;}
	    }
	    int temp_size_deviation = int(log10(abs(osi->unit_factor_result * xcontribution[0].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f])));
	    if (osi->unit_factor_result * xcontribution[0].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f] >= 1.){temp_size_deviation++;}

	    for (int i_c = 0; i_c < xcontribution.size(); i_c++){
     	      int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_result * xcontribution[i_c].result_TSV[i_m][0][i_q][i_s][i_r][i_f])));
	      int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_result * xcontribution[i_c].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f])));
	      if (osi->unit_factor_result * xcontribution[i_c].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f] >= 1.){temp_size_deviation_subprocess++;}
	      if (xcontribution[i_c].result_TSV[i_m][0][i_q][i_s][i_r][i_f] == 0){temp_size_result_subprocess = 0;}
	      if (xcontribution[i_c].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f] == 0){temp_size_deviation_subprocess = 0;}
	      int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
	      int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
	      if (xcontribution[i_c].result_TSV[i_m][0][i_q][i_s][i_r][i_f] == 0.){setw_result = 1; setw_deviation = 1;}
	      else if (abs(osi->unit_factor_result * xcontribution[i_c].result_TSV[i_m][0][i_q][i_s][i_r][i_f]) < 1.){setw_result--;}

	      stringstream temp_ss;
	      temp_ss << left;
	      if (active_qTcut){temp_ss << setw(15) << setprecision(8) << osi->value_qTcut[i_q];}
	      else {temp_ss << setw(15) << "independent";}
	      temp_ss << setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r] 
		      << setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f];
	      if (i_c == 0){temp_ss << setw(25) << xcontribution[i_c].infix_order_contribution;}
	      else {temp_ss << setw(25) << xcontribution[i_c].infix_contribution;}
	      temp_ss << right << setw(15) << setprecision(setw_result) << showpoint 
		      << osi->unit_factor_result * xcontribution[i_c].result_TSV[i_m][0][i_q][i_s][i_r][i_f] 
		      << setw(4) << " +- "
		      << right << setw(15) << setprecision(setw_deviation) 
		      << osi->unit_factor_result * xcontribution[i_c].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f];
	      
	      outfile_complete << temp_ss.str() << endl;
	      if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str() << endl;}
	      if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str() << endl;}
	      if (i_r == i_f){outfile_equal << temp_ss.str() << endl;}
	      if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str() << endl;}

	      if (i_c == 0){
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

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::output_contribution_result_overview_qTcut_CV(){
  Logger logger("summary_list::output_contribution_result_overview_qTcut_CV");
  logger << LOG_DEBUG << "called" << endl;

  string filename;
  ofstream outfile;

  stringstream temp_ss_separation_one;
  stringstream temp_ss_separation_two;
  int counter_separation = 15 + 15 + 15 + 25 + 15 + 4 + 15;
  for (int i_x = 0; i_x < counter_separation; i_x++){
    temp_ss_separation_one << "-";
    temp_ss_separation_two << "=";
  }

  vector<double> rel_scale_CV(osi->n_scales_CV);
  vector<double> scale_CV(osi->n_scales_CV);
  if (osi->n_scales_CV != 0 && osi->n_qTcut != 0){
    for (int s = 0; s < osi->n_scales_CV; s++){
      if (osi->n_scales_CV > 1){rel_scale_CV[s] = pow(10., log10(double(osi->variation_factor_CV)) * double(2. * s / (osi->n_scales_CV - 1) - 1));}
      else {rel_scale_CV[s] = 1;}
      scale_CV[s] = rel_scale_CV[s] * osi->scale_ren;
    }
  }

  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/overview" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";
    logger << LOG_INFO << "FILENAME:   " << filename << endl;
    
    outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    
    stringstream name_ss;
    name_ss << "dXS +- err (contribution)";

    stringstream header_ss;
    header_ss << left
	      << setw(15) << "qTcut"
	      << setw(5) << "mu_R/" << setw(10) << "CV"
	      << setw(5) << "mu_F/" << setw(10) << "CV" 
	      << setw(25) << "contribution"
	      << setw(34) << right << name_ss.str();

    outfile << header_ss.str() << endl;
    
    outfile << endl;
    
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      stringstream temp_ss;
      temp_ss << noshowpoint << left;
      if (active_qTcut){temp_ss << setw(15) << "extrapolation";}
      else {temp_ss << setw(15) << "independent";}
      temp_ss << setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
	      << setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s]
	      << setw(25) << xcontribution[0].infix_order_contribution
	      << showpoint 
	      << right << setprecision(8) << setw(15) << osi->unit_factor_result * result_CV[i_m][0][i_s]
	      << setw(4) << " +- "
	      << right << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_CV[i_m][0][i_s];
      
      outfile << temp_ss.str() << endl;
      
      outfile << endl;
    }
    
    outfile << temp_ss_separation_two.str() << endl << endl;

    for (int i_q = 0; i_q < xcontribution[0].output_n_qTcut; i_q++){
      outfile << header_ss.str() << endl;
      
      outfile << endl;
      
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	int temp_size_result = int(log10(abs(osi->unit_factor_result * xcontribution[0].result_CV[i_m][0][i_q][i_s])));
	for (int i_c = 0; i_c < xcontribution.size(); i_c++){
	  int temp_size_result_c = int(log10(abs(osi->unit_factor_result * xcontribution[i_c].result_CV[i_m][0][i_q][i_s])));
	  if (temp_size_result_c > temp_size_result){temp_size_result = temp_size_result_c;}
	}
	int temp_size_deviation = int(log10(abs(osi->unit_factor_result * xcontribution[0].deviation_CV[i_m][0][i_q][i_s])));
	if (osi->unit_factor_result * xcontribution[0].deviation_CV[i_m][0][i_q][i_s] >= 1.){temp_size_deviation++;}
	
	for (int i_c = 0; i_c < xcontribution.size(); i_c++){
	  int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_result * xcontribution[i_c].result_CV[i_m][0][i_q][i_s])));
	  int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_result * xcontribution[i_c].deviation_CV[i_m][0][i_q][i_s])));
	  if (osi->unit_factor_result * xcontribution[i_c].deviation_CV[i_m][0][i_q][i_s] >= 1.){temp_size_deviation_subprocess++;}
	  if (xcontribution[i_c].result_CV[i_m][0][i_q][i_s] == 0){temp_size_result_subprocess = 0;}
	  if (xcontribution[i_c].deviation_CV[i_m][0][i_q][i_s] == 0){temp_size_deviation_subprocess = 0;}
	  int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
	  int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
	  if (xcontribution[i_c].result_CV[i_m][0][i_q][i_s] == 0.){setw_result = 1;setw_deviation = 1;}
	  else if (abs(osi->unit_factor_result * xcontribution[i_c].result_CV[i_m][0][i_q][i_s]) < 1.){setw_result--;}

	  stringstream temp_ss;
	  temp_ss << noshowpoint << left;
	  if (active_qTcut){temp_ss << setw(15) << setprecision(8) << osi->value_qTcut[i_q];}
	  else {temp_ss << setw(15) << "independent";}
	  temp_ss << setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
		  << setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s];
	  if (i_c == 0){temp_ss << setw(25) << xcontribution[i_c].infix_order_contribution;}
	  else {temp_ss << setw(25) << xcontribution[i_c].infix_contribution;}
	  temp_ss << showpoint 
		  << right << setw(15) << setprecision(setw_result)
		  << osi->unit_factor_result * xcontribution[i_c].result_CV[i_m][0][i_q][i_s]
		  << setw(4) << " +- "
		  << right << setw(15) << setprecision(setw_deviation) 
		  << osi->unit_factor_result * xcontribution[i_c].deviation_CV[i_m][0][i_q][i_s];

	  outfile << temp_ss.str() << endl;

	  if (i_c == 0){
	    outfile << temp_ss_separation_one.str() << endl;
	  }
	}
	outfile << endl;
      }
      outfile << temp_ss_separation_two.str() << endl << endl;
    }
  }
  outfile.close();

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::output_contribution_result_qTcut_TSV(){
  Logger logger("summary_list::output_contribution_result_qTcut_TSV");
  logger << LOG_DEBUG << "called" << endl;

  int start_contribution = 0;
  int end_contribution = 1;

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

  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      for (int i_c = start_contribution; i_c < end_contribution; i_c++){
	filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/complete/" + resultdirectory + "/result" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[i_c].infix_order_contribution + ".dat";
	filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/ren/" + resultdirectory + "/result" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[i_c].infix_order_contribution + ".dat";
	filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/fact/" + resultdirectory + "/result" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[i_c].infix_order_contribution + ".dat";
	filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/equal/" + resultdirectory + "/result" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[i_c].infix_order_contribution + ".dat";
	filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/antipodal/" + resultdirectory + "/result" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[i_c].infix_order_contribution + ".dat";

	logger << LOG_INFO << "FILENAME (complete):   " << filename_complete << endl;

	outfile_complete.open(filename_complete.c_str(), ofstream::out | ofstream::trunc);  
	outfile_ren.open(filename_ren.c_str(), ofstream::out | ofstream::trunc);  
	outfile_fact.open(filename_fact.c_str(), ofstream::out | ofstream::trunc);  
	outfile_equal.open(filename_equal.c_str(), ofstream::out | ofstream::trunc);  
	outfile_antipodal.open(filename_antipodal.c_str(), ofstream::out | ofstream::trunc);  

	stringstream name_ss;
	name_ss << "dXS +- err (" << xcontribution[i_c].infix_order_contribution << ")";

	stringstream header_ss;
	header_ss << left;
	if (xcontribution[i_c].active_qTcut){header_ss << setw(15) << "qTcut";}
	header_ss << setw(5) << "mu_R/" << setw(10) << osi->name_extended_set_TSV[i_s] 
		  << setw(5) << "mu_F/" << setw(10) << osi->name_extended_set_TSV[i_s] 
		  << setw(34) << name_ss.str();
	if (i_c == 0 && xcontribution[i_c].active_qTcut){
	  header_ss << setw(15) << "err(stat)" 
		    << setw(15) << "err(extr)";
	}

	outfile_complete << header_ss.str() << endl;
	outfile_ren << header_ss.str() << endl;
	outfile_fact << header_ss.str() << endl;
	outfile_equal << header_ss.str() << endl;
	outfile_antipodal << header_ss.str() << endl;

	if (i_c == 0 || !(xcontribution[i_c].active_qTcut)){
	  outfile_complete << endl;
	  outfile_ren << endl;
	  outfile_fact << endl;
	  outfile_equal << endl;
	  outfile_antipodal << endl;
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      stringstream temp_ss;
	      temp_ss << left;
	      if (active_qTcut){temp_ss << setw(15) << "extrapolation";}
	      else {temp_ss << setw(15) << "independent";}
	      temp_ss	<< noshowpoint
			<< setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r]
			<< setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f]
			<< showpoint 
			<< setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_s][i_r][i_f]
			<< " +- "
			<< setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][i_r][i_f];
	      if (xcontribution[i_c].active_qTcut){
		temp_ss << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_statistics_TSV[i_m][0][i_s][i_r][i_f]
			<< setprecision(8) << setw(15) << osi->unit_factor_result * deviation_extrapolation_TSV[i_m][0][i_s][i_r][i_f];
	      }

	      outfile_complete << temp_ss.str() << endl;
	      if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str() << endl;}
	      if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str() << endl;}
	      if (i_r == i_f){outfile_equal << temp_ss.str() << endl;}
	      if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str() << endl;}
	    }
	  }
	}
	
	for (int i_q = 0; i_q < xcontribution[i_c].output_n_qTcut; i_q++){
	  outfile_complete << endl;
	  outfile_ren << endl;
	  outfile_fact << endl;
	  outfile_equal << endl;
	  outfile_antipodal << endl;

	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      stringstream temp_ss;
	      temp_ss << left;
	      if (xcontribution[i_c].active_qTcut){temp_ss << setw(15) << osi->value_qTcut[i_q];}
	      temp_ss << setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r] 
		      << setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f] 
		      << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].result_TSV[i_m][0][i_q][i_s][i_r][i_f] 
			<< " +- "
		      << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f];

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
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::output_contribution_result_qTcut_CV(){
  Logger logger("summary_list::output_contribution_result_qTcut_CV");
  logger << LOG_DEBUG << "called" << endl;

  int start_contribution = 0;
  int end_contribution = 1;
  //  if (ygeneric->switch_output_result > 2){end_contribution = xcontribution.size();}

  string filename;
  ofstream outfile;

  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    for (int i_c = start_contribution; i_c < end_contribution; i_c++){
      filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/result" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[i_c].infix_order_contribution + ".dat";

      logger << LOG_INFO << "FILENAME:   " << filename << endl;
      
      outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  

      stringstream name_ss;
      name_ss << "XS +- err (" << xcontribution[i_c].infix_order_contribution << ")";

      stringstream header_ss;
      header_ss << left
		<< setw(15) << "qTcut"
		<< setw(5) << "mu_R/" << setw(10) << "CV" 
		<< setw(5) << "mu_F/" << setw(10) << "CV" 
		<< setw(34) << name_ss.str();
      if (i_c == 0 && xcontribution[i_c].active_qTcut){
	header_ss << setw(15) << "err(stat)" 
		  << setw(15) << "err(extr)";
      }
      outfile << header_ss.str() << endl;
      if (i_c == 0 || !(xcontribution[i_c].active_qTcut)){
	outfile << endl;
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  outfile << left;
	  if (active_qTcut){outfile << setw(15) << "extrapolation";}
	  else {outfile << setw(15) << "independent";}
	  outfile << noshowpoint
		  << setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
		  << setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s]
		  << showpoint 
		  << setprecision(8) << setw(15) << osi->unit_factor_result * result_CV[i_m][0][i_s] 
		  << " +- "
		  << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_CV[i_m][0][i_s];
	  if (xcontribution[i_c].active_qTcut){
	    outfile << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_statistics_CV[i_m][0][i_s]
		    << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_extrapolation_CV[i_m][0][i_s];
	  }
	  outfile << endl;
	}
      }

      if (xcontribution[i_c].active_qTcut){
	for (int i_q = 0; i_q < xcontribution[i_c].output_n_qTcut; i_q++){
	  outfile << endl;
	  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	    stringstream temp_ss;
	    temp_ss << noshowpoint
		    << left;
	    temp_ss << setw(15) << osi->value_qTcut[i_q];
	    temp_ss << setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
		    << setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s]
		    << showpoint 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].result_CV[i_m][0][i_q][i_s]
		    << " +- "
		    << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].deviation_CV[i_m][0][i_q][i_s];
	    outfile << temp_ss.str() << endl;
	  }  
	}
      }
      outfile.close();
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}







// Should be expressed from <list> results !!!

void summary_list::output_contribution_result_plot_qTcut_TSV(){
  Logger logger("summary_list::output_contribution_result_plot_qTcut_TSV");
  logger << LOG_DEBUG << "called" << endl;

  int start_contribution = 0;
  int end_contribution = 1;
  //  if (ygeneric->switch_output_plot > 2){end_contribution = xcontribution.size();}

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

  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      for (int i_c = start_contribution; i_c < end_contribution; i_c++){
	filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/complete/" + resultdirectory + "/plot.qTcut" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[i_c].infix_order_contribution + ".dat";
	filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/ren/" + resultdirectory + "/plot.qTcut" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[i_c].infix_order_contribution + ".dat";
	filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/fact/" + resultdirectory + "/plot.qTcut" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[i_c].infix_order_contribution + ".dat";
	filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/equal/" + resultdirectory + "/plot.qTcut" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[i_c].infix_order_contribution + ".dat";
	filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/antipodal/" + resultdirectory + "/plot.qTcut" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[i_c].infix_order_contribution + ".dat";

	logger << LOG_INFO << "FILENAME (complete):   " << filename_complete << endl;

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
	    temp_ss << "mu_R = " << setprecision(2) << setw(4) << showpoint << osi->relative_scale_ren_TSV[i_s][i_r] << " -- " 
		    << "mu_F = " << setw(4) << showpoint << osi->relative_scale_fact_TSV[i_s][i_f] << "    ";
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

	// write out qTcut = 0 value for qTcut-independent results:

	if (!(xcontribution[i_c].active_qTcut)){
  	  outfile_complete << setw(15) << "0";
  	  outfile_equal << setw(15) << "0";
  	  outfile_antipodal << setw(15) << "0";
  	  outfile_ren << setw(15) << "0";
  	  outfile_fact << setw(15) << "0";
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      stringstream temp_ss;
	      temp_ss << left 
		      << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].result_TSV[i_m][0][0][i_s][i_r][i_f] 
		      << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].deviation_TSV[i_m][0][0][i_s][i_r][i_f];

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
	      temp_ss << left << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].result_TSV[i_m][0][i_q][i_s][i_r][i_f] 
		      << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f];
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

  logger << LOG_DEBUG << "finished" << endl;
}



// Should be expressed from <list> results !!!

void summary_list::output_contribution_result_plot_qTcut_CV(){
  Logger logger("summary_list::output_contribution_result_plot_qTcut_CV");
  logger << LOG_DEBUG << "called" << endl;

  int start_contribution = 0;
  int end_contribution = 1;

  string filename;

  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    for (int i_c = start_contribution; i_c < end_contribution; i_c++){
      ofstream outfile;
      
      if (osi->n_qTcut > 1){
	///////////////////////////////////////
	//  output of usual qTcut variation  //
	///////////////////////////////////////
	filename = ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + xcontribution[i_c].infix_order_contribution + ".dat";

	logger << LOG_INFO << "FILENAME:   " << filename << endl;

	outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  

	outfile << noshowpoint;
	outfile << left << setw(15) << "# qTcut";
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  stringstream temp_ss;
	  temp_ss << "mu_R = " 
		  << setprecision(2) << setw(4) << osi->rel_scale_factor_ren_CV[i_s] 
		  << " -- " 
		  << "mu_F = " 
		  << setprecision(2) << setw(4) << osi->rel_scale_factor_fact_CV[i_s] << "    ";
	  outfile << setw(30) << temp_ss.str() << left;
	}
	outfile << endl;
	
	if (!(xcontribution[i_c].active_qTcut)){
	  outfile << setw(15) << "0";
	  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	    stringstream temp_ss;
	    temp_ss << left 
		    << showpoint 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].result_CV[i_m][0][0][i_s]
		    << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].deviation_CV[i_m][0][0][i_s];
	    outfile << temp_ss.str();
	  }
	  outfile << endl;
	}
	
	for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	  outfile << noshowpoint 
		  << setprecision(4) << setw(15) << osi->value_qTcut[i_q];
	  for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	    outfile << showpoint
		    << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].result_CV[i_m][0][i_q][i_s] 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[i_c].deviation_CV[i_m][0][i_q][i_s];
	  }
	  outfile << endl;
	}
	outfile.close();
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::output_contribution_result_plot_TSV(){
  Logger logger("summary_list::output_contribution_result_plot_TSV");
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

  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/complete/" + resultdirectory + "/plot" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";
      filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/ren/" + resultdirectory + "/plot" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";
      filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/fact/" + resultdirectory + "/plot" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";
      filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/equal/" + resultdirectory + "/plot" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";
      filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/antipodal/" + resultdirectory + "/plot" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";

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
	  temp_ss << showpoint 
		  << "mu_R = " 
		  << setw(4) << setprecision(2) << osi->relative_scale_ren_TSV[i_s][i_r] 
		  << " -- " 
		  << "mu_F = " 
		  << setw(4) << setprecision(2) << osi->relative_scale_fact_TSV[i_s][i_f] << "    ";

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

      outfile_complete << setw(15) << "0";
      outfile_equal << setw(15) << "0";
      outfile_antipodal << setw(15) << "0";
      outfile_ren << setw(15) << "0";
      outfile_fact << setw(15) << "0";

      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  stringstream temp_ss;
	  temp_ss << left 
		  << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_s][i_r][i_f] 
		  << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][i_r][i_f];

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
      
      outfile_complete.close();
      outfile_ren.close();
      outfile_fact.close();
      outfile_equal.close();
      outfile_antipodal.close();
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::output_contribution_result_plot_CV(){
  Logger logger("summary_list::output_contribution_result_plot_CV");
  logger << LOG_DEBUG << "called" << endl;

  string filename;
  
  vector<double> rel_scale_CV(osi->n_scales_CV);
  vector<double> scale_CV(osi->n_scales_CV);
  if (osi->n_scales_CV != 0 && osi->n_qTcut != 0){
    for (int s = 0; s < osi->n_scales_CV; s++){
      if (osi->n_scales_CV > 1){rel_scale_CV[s] = pow(10., log10(double(osi->variation_factor_CV)) * double(2. * s / (osi->n_scales_CV - 1) - 1));}
      else {rel_scale_CV[s] = 1;}
      scale_CV[s] = rel_scale_CV[s] * osi->scale_ren;
    }
  }
  
  int i_q = 0;

  int start_i_m = 0;
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    ofstream outfile;

    ////////////////////////////////////
    //  output of usual CV variation  //
    ////////////////////////////////////
    filename = ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/plot" + ygeneric->infix_name_moment[i_m] + ".CV." + xcontribution[0].infix_order_contribution + ".dat";

    logger << LOG_INFO << "FILENAME:   " << filename << endl;

    outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  

    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      outfile << setw(15) << setprecision(8) << rel_scale_CV[i_s] 
	      << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[0].result_CV[i_m][0][i_q][i_s] 
	      << setprecision(8) << setw(15) << osi->unit_factor_result * xcontribution[0].deviation_CV[i_m][0][i_q][i_s] 
	      << endl;
    }
    outfile.close();
    
    ////////////////////////////////////////////////////////////////////////
    //   output of summed NLO contributions on Born phase-space (BVIKP)   //
    ////////////////////////////////////////////////////////////////////////
    vector<int> temp_exclude_real;
    for (int i_c = 0; i_c < xcontribution.size(); i_c++){
      if (xcontribution[i_c].type_contribution == "RA"){temp_exclude_real.push_back(i_c);}
    }
    if (temp_exclude_real.size() != 0){
      filename = ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/plot" + ygeneric->infix_name_moment[i_m] + ".CV." + xcontribution[0].infix_order_contribution + ".VIKP.dat";

      logger << LOG_INFO << "FILENAME:   " << filename << endl;

      outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  
      
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	logger << LOG_DEBUG_VERBOSE << "i_q = " << i_q << "   i_s = " << i_s << endl;
	double temp_result = 0.;
	double temp_deviation2 = 0.;
	for (int i_c = 1; i_c < xcontribution.size(); i_c++){
	  //	for (int i_c = 1; i_c < temp_result_CV.size(); i_c++){
	  int flag = temp_exclude_real.size();
	  for (int j_c = 0; j_c < temp_exclude_real.size(); j_c++){if (temp_exclude_real[j_c] == i_c){flag = i_c;}}
	  if (flag == temp_exclude_real.size()){
	    temp_result += xcontribution[i_c].result_CV[i_m][0][i_q][i_s];
	    temp_deviation2 += pow(xcontribution[i_c].deviation_CV[i_m][0][i_q][i_s], 2);
	  }
	}
	outfile << setw(15) << setprecision(8) << rel_scale_CV[i_s] << setprecision(8) << setw(15) << osi->unit_factor_result * temp_result << setprecision(8) << setw(15) << osi->unit_factor_result * sqrt(temp_deviation2) << endl;
      }
      outfile.close();
      logger << LOG_DEBUG_VERBOSE << filename << " closed." << endl;
    }

    //////////////////////////////////////////////////////////////////
    //  output of extrapolated (or qTcut-independent) CV variation  //
    //////////////////////////////////////////////////////////////////
    filename = ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/plot" + ygeneric->infix_name_moment[i_m] + ".CV.extrapolated." + xcontribution[0].infix_order_contribution + ".dat";

    logger << LOG_INFO << "FILENAME:   " << filename << endl;

    outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  

    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      logger << LOG_DEBUG_VERBOSE << "i_q = " << i_q << "   i_s = " << i_s << endl;
      outfile << noshowpoint
	      << setw(15) << setprecision(8) << rel_scale_CV[i_s] 
	      << showpoint 
	      << setprecision(8) << setw(15) << osi->unit_factor_result * result_CV[i_m][0][i_s] 
	      << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_CV[i_m][0][i_s] 
	      << endl;
    }
    outfile.close();
    
    ///////////////////////////////////////////////////////////////////////////
    //  standard output of extrapolated (or qTcut-independent) CV variation  //
    ///////////////////////////////////////////////////////////////////////////
    filename = ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/" + resultdirectory + "/plot" + ygeneric->infix_name_moment[i_m] + "." + xcontribution[0].infix_order_contribution + ".dat";

    logger << LOG_INFO << "FILENAME:   " << filename << endl;

    outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  

    outfile << left << setw(15) << "# qTcut";
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      stringstream temp_ss;
      temp_ss << noshowpoint
	      << "mu_R = " 
	      << setprecision(2) << setw(4) << showpoint << osi->rel_scale_factor_ren_CV[i_s] 
	      << " -- " 
	      << "mu_F = " 
	      << setprecision(2) << setw(4) << showpoint << osi->rel_scale_factor_fact_CV[i_s] << "    ";
      outfile << setw(30) << temp_ss.str() << left;
    }
    outfile << endl;

    outfile << noshowpoint
	    << setw(15) << osi->value_qTcut[i_q];
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      outfile << showpoint 
	      << setprecision(8) << setw(15) << osi->unit_factor_result * result_CV[i_m][0][i_s] 
	      << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_CV[i_m][0][i_s];
    }
    outfile << endl;
    outfile.close();
  }

  logger << LOG_DEBUG << "finished" << endl;
}





  /*
  if (osi->n_moments > 50){
    logger << LOG_DEBUG<< "Moments are considered as a parameter scan!" << endl;
    start_i_m = osi->n_moments;
  }
  */
