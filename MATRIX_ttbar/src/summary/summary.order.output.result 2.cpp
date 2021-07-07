#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_order::output_result_overview_TSV(){
  Logger logger("summary_order::output_result_overview_TSV");
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
      system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s]);

      filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/complete/overview" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/ren/overview" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/fact/overview" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/equal/overview" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/antipodal/overview" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";

      outfile_complete.open(filename_complete.c_str(), ofstream::out | ofstream::trunc);  
      outfile_ren.open(filename_ren.c_str(), ofstream::out | ofstream::trunc);  
      outfile_fact.open(filename_fact.c_str(), ofstream::out | ofstream::trunc);  
      outfile_equal.open(filename_equal.c_str(), ofstream::out | ofstream::trunc);  
      outfile_antipodal.open(filename_antipodal.c_str(), ofstream::out | ofstream::trunc);  
      
      stringstream name_ss;
      name_ss << "dXS +- err (list)";

      stringstream header_ss;
      header_ss << left
		<< setw(15) << "qTcut"
		<< setw(5) << "mu_R/" << setw(10) << osi->name_extended_set_TSV[i_s] 
		<< setw(5) << "mu_F/" << setw(10) << osi->name_extended_set_TSV[i_s]
		<< setw(25) << "list"
		<< right 
		<< setw(34) << name_ss.str();
      
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
	  int temp_size_result_0 = int(log10(abs(osi->unit_factor_result * result_TSV[i_m][0][i_s][i_r][i_f])));
	  int temp_size_result = int(log10(abs(osi->unit_factor_result * result_TSV[i_m][0][i_s][i_r][i_f])));
	  for (int i_l = 0; i_l < xlist.size(); i_l++){
	    int temp_size_result_l = int(log10(abs(osi->unit_factor_result * xlist[i_l]->result_TSV[i_m][0][i_s][i_r][i_f])));
	    if (temp_size_result_l > temp_size_result){temp_size_result = temp_size_result_l;}
	  }
	  int temp_size_deviation = int(log10(abs(osi->unit_factor_result * deviation_TSV[i_m][0][i_s][i_r][i_f])));
	  if (osi->unit_factor_result * deviation_TSV[i_m][0][i_s][i_r][i_f] >= 1.){temp_size_deviation++;}
	  
	  stringstream temp_ss;
	  temp_ss << left;
	  if (active_qTcut){temp_ss << setw(15) << "extrapolation";}
	  else {temp_ss << setw(15) << "independent";}
	  temp_ss << noshowpoint
		  << setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r] 
		  << setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f] 
		  << setw(25) << resultdirectory
		  << showpoint 
		  << right << setw(15) << setprecision(10 + temp_size_result_0 - temp_size_result) << showpoint 
		  << osi->unit_factor_result * result_TSV[i_m][0][i_s][i_r][i_f] 
		  << " +- "
		  << right << setw(15) << setprecision(10 + temp_size_deviation - temp_size_deviation - (temp_size_result - temp_size_deviation + 1)) 
		  << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][i_r][i_f];
	  
	  outfile_complete << temp_ss.str() << endl;
	  if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str() << endl;}
	  if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str() << endl;}
	  if (i_r == i_f){outfile_equal << temp_ss.str() << endl;}
	  if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str() << endl;}
	  
	  outfile_complete << temp_ss_separation_one.str() << endl;
	  if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss_separation_one.str() << endl;}
	  if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss_separation_one.str() << endl;}
	  if (i_r == i_f){outfile_equal << temp_ss_separation_one.str() << endl;}
	  if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss_separation_one.str() << endl;}
	  
	  for (int i_l = 0; i_l < xlist.size(); i_l++){
	    int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_result * xlist[i_l]->result_TSV[i_m][0][i_s][i_r][i_f])));
	    int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_result * xlist[i_l]->deviation_TSV[i_m][0][i_s][i_r][i_f])));
	    if (osi->unit_factor_result * xlist[i_l]->deviation_TSV[i_m][0][i_s][i_r][i_f] >= 1.){temp_size_deviation_subprocess++;}
	    if (xlist[i_l]->result_TSV[i_m][0][i_s][i_r][i_f] == 0){temp_size_result_subprocess = 0;}
	    if (xlist[i_l]->deviation_TSV[i_m][0][i_s][i_r][i_f] == 0){temp_size_deviation_subprocess = 0;}
	    int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
	    int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
	    if (xlist[i_l]->result_TSV[i_m][0][i_s][i_r][i_f] == 0.){setw_result = 1;setw_deviation = 1;}
	    else if (abs(osi->unit_factor_result * xlist[i_l]->result_TSV[i_m][0][i_s][i_r][i_f]) < 1.){setw_result--;}
	    
	    stringstream temp_ss;
	    temp_ss << left;
	    if (xlist[i_l]->active_qTcut){temp_ss << setw(15) << "extrapolation";}
	    else {temp_ss << setw(15) << "independent";}
	    
	    temp_ss << noshowpoint
		    << setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r] 
		    << setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f] 
		    << setw(25) << xlist[i_l]->xcontribution[0].infix_order_contribution
		    << showpoint 
		    << right << setw(15) << setprecision(setw_result) << showpoint 
		    << osi->unit_factor_result * xlist[i_l]->result_TSV[i_m][0][i_s][i_r][i_f] 
		    << " +- "
		    << right << setw(15) << setprecision(setw_deviation) 
		    << osi->unit_factor_result * xlist[i_l]->deviation_TSV[i_m][0][i_s][i_r][i_f];

	    outfile_complete << temp_ss.str() << endl;
	    if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str() << endl;}
	    if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str() << endl;}
	    if (i_r == i_f){outfile_equal << temp_ss.str() << endl;}
	    if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str() << endl;}
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
      
      if (active_qTcut){
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
	      int temp_size_result_0 = int(log10(abs(osi->unit_factor_result * result_qTcut_TSV[i_m][0][i_q][i_s][i_r][i_f])));
	      int temp_size_result = int(log10(abs(osi->unit_factor_result * result_qTcut_TSV[i_m][0][i_q][i_s][i_r][i_f])));
	      for (int i_l = 0; i_l < xlist.size(); i_l++){
		int temp_size_result_l = int(log10(abs(osi->unit_factor_result * xlist[i_l]->xcontribution[0].result_TSV[i_m][0][i_q][i_s][i_r][i_f])));
		if (temp_size_result_l > temp_size_result){temp_size_result = temp_size_result_l;}
	      }
	      int temp_size_deviation = int(log10(abs(osi->unit_factor_result * deviation_qTcut_TSV[i_m][0][i_q][i_s][i_r][i_f])));
	      if (osi->unit_factor_result * deviation_qTcut_TSV[i_m][0][i_q][i_s][i_r][i_f] >= 1.){temp_size_deviation++;}
	      
	      stringstream temp_ss;
	      temp_ss << left 
		      << noshowpoint
		      << setw(15) << osi->value_qTcut[i_q]
		      << setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r] 
		      << setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f] 
		      << setw(25) << resultdirectory
		      << showpoint 
		      << right << setw(15) << setprecision(10 + temp_size_result_0 - temp_size_result) << showpoint 
		      << osi->unit_factor_result * result_qTcut_TSV[i_m][0][i_q][i_s][i_r][i_f] 
		      << " +- "
		      << right << setw(15) << setprecision(10 + temp_size_deviation - temp_size_deviation - (temp_size_result - temp_size_deviation + 1)) 
		      << osi->unit_factor_result * deviation_qTcut_TSV[i_m][0][i_q][i_s][i_r][i_f];
	      
	      outfile_complete << temp_ss.str() << endl;
	      if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str() << endl;}
	      if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str() << endl;}
	      if (i_r == i_f){outfile_equal << temp_ss.str() << endl;}
	      if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str() << endl;}
	      
	      outfile_complete << temp_ss_separation_one.str() << endl;
	      if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss_separation_one.str() << endl;}
	      if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss_separation_one.str() << endl;}
	      if (i_r == i_f){outfile_equal << temp_ss_separation_one.str() << endl;}
	      if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss_separation_one.str() << endl;}
	      
	      for (int i_l = 0; i_l < xlist.size(); i_l++){
		int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_result * xlist[i_l]->xcontribution[0].result_TSV[i_m][0][i_q][i_s][i_r][i_f])));
		int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_result * xlist[i_l]->xcontribution[0].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f])));
		if (osi->unit_factor_result * xlist[i_l]->xcontribution[0].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f] >= 1.){temp_size_deviation_subprocess++;}
		if (xlist[i_l]->xcontribution[0].result_TSV[i_m][0][i_q][i_s][i_r][i_f] == 0){temp_size_result_subprocess = 0;}
		if (xlist[i_l]->xcontribution[0].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f] == 0){temp_size_deviation_subprocess = 0;}
		int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
		int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
		if (xlist[i_l]->xcontribution[0].result_TSV[i_m][0][i_q][i_s][i_r][i_f] == 0.){setw_result = 1; setw_deviation = 1;}
		else if (abs(osi->unit_factor_result * xlist[i_l]->xcontribution[0].result_TSV[i_m][0][i_q][i_s][i_r][i_f]) < 1.){setw_result--;}
		
		stringstream temp_l_ss;
		temp_l_ss << noshowpoint
			  << left << setw(15) << osi->value_qTcut[i_q]
			  << setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r] 
			  << setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f] 
			  << setw(25) << xlist[i_l]->xcontribution[0].infix_order_contribution
			  << showpoint 
			  << right << setw(15) << setprecision(setw_result) << showpoint 
			  << osi->unit_factor_result * xlist[i_l]->xcontribution[0].result_TSV[i_m][0][i_q][i_s][i_r][i_f] 
			  << " +- "
			  << right << setw(15) << setprecision(setw_deviation) 
			  << osi->unit_factor_result * xlist[i_l]->xcontribution[0].deviation_TSV[i_m][0][i_q][i_s][i_r][i_f];

		outfile_complete << temp_l_ss.str() << endl;
		if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_l_ss.str() << endl;}
		if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_l_ss.str() << endl;}
		if (i_r == i_f){outfile_equal << temp_l_ss.str() << endl;}
		if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_l_ss.str() << endl;}
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



void summary_order::output_result_overview_CV(){
  Logger logger("summary_order::output_result_overview_CV");
  logger << LOG_DEBUG << "called" << endl;

  string filename;
  ofstream outfile;
  stringstream temp_ss_separation_one;
  stringstream temp_ss_separation_two;
  int counter_separation = 0 + 15 + 15 + 15 + 25 + 15 + 4 + 15;
  for (int i_x = 0; i_x < counter_separation; i_x++){
    temp_ss_separation_one << "-";
    temp_ss_separation_two << "=";
  }
  
  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV");

    filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/overview" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
    outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    
    stringstream name_ss;
    name_ss << "dXS +- err (list)";
    
    stringstream header_ss;
    header_ss << left
	      << setw(15) << "qTcut"
	      << setw(5) << "mu_R/" << setw(10) << "CV"
	      << setw(5) << "mu_F/" << setw(10) << "CV"
	      << setw(25) << "list"
	      << right 
	      << setw(34) << name_ss.str();
    outfile << header_ss.str() << endl;
    outfile << endl;
    
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      int temp_size_result_0 = int(log10(abs(osi->unit_factor_result * result_CV[i_m][0][i_s])));
      int temp_size_result = int(log10(abs(osi->unit_factor_result * result_CV[i_m][0][i_s])));
      for (int i_l = 0; i_l < xlist.size(); i_l++){
	int temp_size_result_l = int(log10(abs(osi->unit_factor_result * xlist[i_l]->result_CV[i_m][0][i_s])));
	if (temp_size_result_l > temp_size_result){temp_size_result = temp_size_result_l;}
      }
      int temp_size_deviation = int(log10(abs(osi->unit_factor_result * deviation_CV[i_m][0][i_s])));
      if (osi->unit_factor_result * deviation_CV[i_m][0][i_s] >= 1.){temp_size_deviation++;}
      int setw_result = 9 + temp_size_result_0 - temp_size_result;
      int setw_deviation = 9 + temp_size_deviation - temp_size_result - 1;
      if (result_CV[i_m][0][i_s] == 0.){setw_result = 1; setw_deviation = 1;}
      else if (abs(osi->unit_factor_result * result_CV[i_m][0][i_s]) < 1.){setw_result--;}
      
      stringstream temp_ss;
      temp_ss << left << noshowpoint;
      if (active_qTcut){temp_ss << setw(15) << "extrapolation";}
      else {temp_ss << setw(15) << "independent";}
      temp_ss << setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
	      << setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s]
 	      << setw(25) << resultdirectory
	      << showpoint 
	      << right << setw(15) << setprecision(setw_result) << showpoint 
	      << osi->unit_factor_result * result_CV[i_m][0][i_s] 
	      << " +- "
	      << right << setw(15) << setprecision(setw_deviation) 
	      << osi->unit_factor_result * deviation_CV[i_m][0][i_s];
      outfile << temp_ss.str() << endl;
      
      outfile << temp_ss_separation_one.str() << endl;
      
      for (int i_l = 0; i_l < xlist.size(); i_l++){
	int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_result * xlist[i_l]->result_CV[i_m][0][i_s])));
	int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_result * xlist[i_l]->deviation_CV[i_m][0][i_s])));
	if (osi->unit_factor_result * xlist[i_l]->deviation_CV[i_m][0][i_s] >= 1.){temp_size_deviation_subprocess++;}
	if (xlist[i_l]->result_CV[i_m][0][i_s] == 0){temp_size_result_subprocess = 0;}
	if (xlist[i_l]->deviation_CV[i_m][0][i_s] == 0){temp_size_deviation_subprocess = 0;}
	int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
	int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);
	if (xlist[i_l]->result_CV[i_m][0][i_s] == 0.){setw_result = 1; setw_deviation = 1;}
	else if (abs(osi->unit_factor_result * xlist[i_l]->result_CV[i_m][0][i_s]) < 1.){setw_result--;}

	stringstream temp_ss;
	temp_ss << left;
	if (xlist[i_l]->active_qTcut){temp_ss << setw(15) << "extrapolation";}
	else {temp_ss << setw(15) << "independent";}
	temp_ss << noshowpoint
		<< setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
		<< setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s]
		<< setw(25) << xlist[i_l]->xcontribution[0].infix_order_contribution
		<< showpoint 
		<< right << setw(15) << setprecision(setw_result) 
		<< osi->unit_factor_result * xlist[i_l]->result_CV[i_m][0][i_s] 
		<< " +- "
		<< right << setw(15) << setprecision(setw_deviation) 
		<< osi->unit_factor_result * xlist[i_l]->deviation_CV[i_m][0][i_s];
	outfile << temp_ss.str() << endl;
      }
      outfile << endl;
    }
  
    outfile << temp_ss_separation_two.str() << endl << endl;
    
    if (active_qTcut){
      for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	outfile << header_ss.str() << endl;
	outfile << endl;
	
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  int temp_size_result_0 = int(log10(abs(osi->unit_factor_result * result_qTcut_CV[i_m][0][i_q][i_s])));
	  int temp_size_result = int(log10(abs(osi->unit_factor_result * result_qTcut_CV[i_m][0][i_q][i_s])));
	  for (int i_l = 0; i_l < xlist.size(); i_l++){
	    int temp_size_result_l = int(log10(abs(osi->unit_factor_result * xlist[i_l]->xcontribution[0].result_CV[i_m][0][i_q][i_s])));
	    if (temp_size_result_l > temp_size_result){temp_size_result = temp_size_result_l;}
	  }
	  int temp_size_deviation = int(log10(abs(osi->unit_factor_result * deviation_qTcut_CV[i_m][0][i_q][i_s])));
	  if (osi->unit_factor_result * deviation_qTcut_CV[i_m][0][i_q][i_s] >= 1.){temp_size_deviation++;}
	  
	  
	  stringstream temp_ss;
	  temp_ss << left 
		  << noshowpoint
		  << setw(15) << osi->value_qTcut[i_q]
		  << setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
		  << setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s]
 		  << setw(25) << resultdirectory
		  << showpoint 
		  << right << setw(15) << setprecision(10 + temp_size_result_0 - temp_size_result) << showpoint 
   		  << osi->unit_factor_result * result_qTcut_CV[i_m][0][i_q][i_s] 
		  << " +- "
		  << right << setw(15) << setprecision(10 + temp_size_deviation - temp_size_deviation - (temp_size_result - temp_size_deviation + 1)) 
		  << osi->unit_factor_result * deviation_qTcut_CV[i_m][0][i_q][i_s];

	  outfile << temp_ss.str() << endl;
	  
	  outfile << temp_ss_separation_one.str() << endl;
	  
	  for (int i_l = 0; i_l < xlist.size(); i_l++){
	    int temp_size_result_subprocess = int(log10(abs(osi->unit_factor_result * xlist[i_l]->xcontribution[0].result_CV[i_m][0][i_q][i_s])));
	    //	    if (temp_size_result_subprocess < -3){temp_size_result_subprocess -= 5;}
	    int temp_size_deviation_subprocess = int(log10(abs(osi->unit_factor_result * xlist[i_l]->xcontribution[0].deviation_CV[i_m][0][i_q][i_s])));
	    if (osi->unit_factor_result * xlist[i_l]->xcontribution[0].deviation_CV[i_m][0][i_q][i_s] >= 1.){temp_size_deviation_subprocess++;}
	    //	    if (temp_size_deviation_subprocess < -3){temp_size_deviation_subprocess -= 5;}
	    if (xlist[i_l]->xcontribution[0].result_CV[i_m][0][i_q][i_s] == 0){temp_size_result_subprocess = 0;}
	    if (xlist[i_l]->xcontribution[0].deviation_CV[i_m][0][i_q][i_s] == 0){temp_size_deviation_subprocess = 0;}
	    
	    int setw_result = 9 + temp_size_result_subprocess - temp_size_result;
	    int setw_deviation = 9 + temp_size_deviation_subprocess - temp_size_deviation - (temp_size_result - temp_size_deviation + 1);

	    stringstream temp_l_ss;
	    temp_l_ss << left 
		      << noshowpoint
		      << setw(15) << osi->value_qTcut[i_q]
		      << setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
		      << setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s]
		      << setw(25) << xlist[i_l]->xcontribution[0].infix_order_contribution
		      << showpoint 
		      << right << setw(15) << setprecision(setw_result) << showpoint 
		      << osi->unit_factor_result * xlist[i_l]->xcontribution[0].result_CV[i_m][0][i_q][i_s] 
		      << " +- "
		      << right << setw(15) << setprecision(setw_deviation) 
		      << osi->unit_factor_result * xlist[i_l]->xcontribution[0].deviation_CV[i_m][0][i_q][i_s];

	    outfile << temp_l_ss.str() << endl;
	  }
	  outfile << endl;
	}
    	outfile << temp_ss_separation_two.str() << endl << endl;
      }
    }
    outfile.close();
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_order::output_result_TSV(){
  Logger logger("summary_order::output_result_TSV");
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
      system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s]);

      filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/complete/result" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/ren/result" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/fact/result" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/equal/result" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/antipodal/result" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";

      outfile_complete.open(filename_complete.c_str(), ofstream::out | ofstream::trunc);  
      outfile_ren.open(filename_ren.c_str(), ofstream::out | ofstream::trunc);  
      outfile_fact.open(filename_fact.c_str(), ofstream::out | ofstream::trunc);  
      outfile_equal.open(filename_equal.c_str(), ofstream::out | ofstream::trunc);  
      outfile_antipodal.open(filename_antipodal.c_str(), ofstream::out | ofstream::trunc);  
      
      stringstream header_ss;
      header_ss << noshowpoint << left
		<< setw(15) << "qTcut"
		<< setw(5) << "mu_R/" << setw(10) << osi->name_extended_set_TSV[i_s] 
		<< setw(5) << "mu_F/" << setw(10) << osi->name_extended_set_TSV[i_s] 
		<< "XS +- err (" << resultdirectory << ")";
      
      outfile_complete << header_ss.str() << endl << endl;
      outfile_ren << header_ss.str() << endl << endl;
      outfile_fact << header_ss.str() << endl << endl;
      outfile_equal << header_ss.str() << endl << endl;
      outfile_antipodal << header_ss.str() << endl << endl;
      
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  stringstream temp_ss;
	  temp_ss << left ;
	  if (active_qTcut){temp_ss << setw(15) << "extrapolation";}
	  else {temp_ss << setw(15) << "independent";}
	  //		 temp_ss << setw(15) << "0"
	  temp_ss << setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r] 
		  << setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f] 
		  << showpoint 
		  << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_s][i_r][i_f] 
		  << " +- "
		  << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][i_r][i_f];
	  
	  outfile_complete << temp_ss.str() << endl;
	  if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str() << endl;}
	  if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str() << endl;}
	  if (i_r == i_f){outfile_equal << temp_ss.str() << endl;}
	  if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str() << endl;}
	}
      }
      
      if (active_qTcut){
	for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	  outfile_complete << endl;
	  outfile_ren << endl;
	  outfile_fact << endl;
	  outfile_equal << endl;
	  outfile_antipodal << endl;
	  
	  for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	      stringstream temp_ss;
	      temp_ss << noshowpoint
		      << left << setw(15) << osi->value_qTcut[i_q]
		      << setw(15) << setprecision(8) << osi->relative_scale_ren_TSV[i_s][i_r] 
		      << setw(15) << setprecision(8) << osi->relative_scale_fact_TSV[i_s][i_f] 
		      << showpoint 
		      << setprecision(8) << setw(15) << osi->unit_factor_result * result_qTcut_TSV[i_m][0][i_q][i_s][i_r][i_f] 
		      << " +- "
		      << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_qTcut_TSV[i_m][0][i_q][i_s][i_r][i_f];
	      
	      outfile_complete << temp_ss.str() << endl;
	      if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str() << endl;}
	      if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str() << endl;}
	      if (i_r == i_f){outfile_equal << temp_ss.str() << endl;}
	      if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str() << endl;}
	    }
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



void summary_order::output_result_CV(){
  Logger logger("summary_order::output_result_CV");
  logger << LOG_DEBUG << "called" << endl;

  string filename;
  ofstream outfile;

  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/CV");
    
    filename = "" + ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/result" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
    outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    
    stringstream header_ss;
    header_ss << noshowpoint << left
	      << setw(15) << "qTcut"
	      << setw(5) << "mu_R/" << setw(10) << "CV" 
	      << setw(5) << "mu_F/" << setw(10) << "CV"
	      << "XS +- err (" << resultdirectory << ")";
    
    outfile << header_ss.str() << endl << endl;
    
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      stringstream temp_ss;
      temp_ss << left;
      if (active_qTcut){temp_ss << setw(15) << "extrapolation";}
      else {temp_ss << setw(15) << "independent";}
      //	      << setw(15) << "0"
      temp_ss << setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
	      << setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s]
	      << showpoint 
	      << setprecision(8) << setw(15) << osi->unit_factor_result * result_CV[i_m][0][i_s]
	      << " +- "
	      << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_CV[i_m][0][i_s];
      
      outfile << temp_ss.str() << endl;
    }
    
    if (active_qTcut){
      for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	outfile << endl;
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  stringstream temp_ss;
	  temp_ss << noshowpoint
		  << left 
		  << setw(15) << osi->value_qTcut[i_q]
		  << setw(15) << setprecision(8) << osi->rel_scale_factor_ren_CV[i_s]
		  << setw(15) << setprecision(8) << osi->rel_scale_factor_fact_CV[i_s]
		  << showpoint 
		  << setprecision(8) << setw(15) << osi->unit_factor_result * result_qTcut_CV[i_m][0][i_q][i_s]
		  << " +- "
		  << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_qTcut_CV[i_m][0][i_q][i_s];
	  
	  outfile << temp_ss.str() << endl;
	}
      }  
      outfile.close();
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_order::output_result_plot_TSV(){
  Logger logger("summary_order::output_result_plot_TSV");
  logger << LOG_DEBUG << "called" << endl;

  string filename_complete;
  string filename_ren;
  string filename_fact;
  string filename_equal;
  string filename_antipodal;
  string filename_seven_point;
  string filename_nine_point;
  
  ofstream outfile_complete;
  ofstream outfile_ren;
  ofstream outfile_fact;
  ofstream outfile_equal;
  ofstream outfile_antipodal;
  ofstream outfile_seven_point;
  ofstream outfile_nine_point;

  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < osi->n_moments + 1; i_m++){
    for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
      if (!ygeneric->switch_output_scaleset[i_s]){continue;}
      system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s]);

      filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/complete/plot" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/ren/plot" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/fact/plot" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/equal/plot" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/antipodal/plot" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_seven_point = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[5] + "/plot" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_nine_point = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/" + ygeneric->name_scale_variation_TSV[6] + "/plot" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";

      logger << LOG_DEBUG_VERBOSE << "filename_complete  = " << filename_complete  << endl;
      logger << LOG_DEBUG_VERBOSE << "filename_ren = " << filename_ren << endl;
      logger << LOG_DEBUG_VERBOSE << "filename_fact = " << filename_fact << endl;
      logger << LOG_DEBUG_VERBOSE << "filename_antipodal = " << filename_antipodal << endl;
      logger << LOG_DEBUG_VERBOSE << "filename_seven_point = " << filename_seven_point << endl;
      logger << LOG_DEBUG_VERBOSE << "filename_nine_point = " << filename_nine_point << endl;

      
      outfile_complete.open(filename_complete.c_str(), ofstream::out | ofstream::trunc);  
      outfile_ren.open(filename_ren.c_str(), ofstream::out | ofstream::trunc);  
      outfile_fact.open(filename_fact.c_str(), ofstream::out | ofstream::trunc);  
      outfile_equal.open(filename_equal.c_str(), ofstream::out | ofstream::trunc);  
      outfile_antipodal.open(filename_antipodal.c_str(), ofstream::out | ofstream::trunc);  
      outfile_seven_point.open(filename_seven_point.c_str(), ofstream::out | ofstream::trunc);  
      outfile_nine_point.open(filename_nine_point.c_str(), ofstream::out | ofstream::trunc);  
     
      outfile_complete << left << setw(15) << "#";
      outfile_ren << left << setw(15) << "#";
      outfile_fact << left << setw(15) << "#";
      outfile_equal << left << setw(15) << "#";
      outfile_antipodal << left << setw(15) << "#";
      
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  stringstream temp_ss;
	  temp_ss << noshowpoint
		  << left 
		  << "mu_R = " 
		  << setprecision(2) << setw(4) << showpoint << osi->relative_scale_ren_TSV[i_s][i_r] 
		  << " -- " 
		  << "mu_F = " 
		  << setw(4) << showpoint << osi->relative_scale_fact_TSV[i_s][i_f] << "    ";

	  outfile_complete << setw(30) << temp_ss.str() << left;
	  if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << setw(30) << temp_ss.str() << left;}
	  if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << setw(30) << temp_ss.str() << left;}
	  if (i_r == i_f){outfile_equal << setw(30) << temp_ss.str() << left;}
	  if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << setw(30) << temp_ss.str() << left;}
	}
      }

      double result_min_nine_point = +1.e99;
      double result_max_nine_point = -1.e99;
      double x_r_min_nine_point = 0;
      double x_r_max_nine_point = 0;
      double x_f_min_nine_point = 0;
      double x_f_max_nine_point = 0;
      
      double result_min_seven_point = +1.e99;
      double result_max_seven_point = -1.e99;
      double x_r_min_seven_point = 0;
      double x_r_max_seven_point = 0;
      double x_f_min_seven_point = 0;
      double x_f_max_seven_point = 0;
      
      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  if (result_TSV[i_m][0][i_s][i_r][i_f] < result_min_nine_point){
	    result_min_nine_point = result_TSV[i_m][0][i_s][i_r][i_f];
	    x_r_min_nine_point = i_r;
	    x_f_min_nine_point = i_f;
	  }
	  if (result_TSV[i_m][0][i_s][i_r][i_f] > result_max_nine_point){
	    result_max_nine_point = result_TSV[i_m][0][i_s][i_r][i_f];
	    x_r_max_nine_point = i_r;
	    x_f_max_nine_point = i_f;
	  }

	  if ((i_r == 0 && i_f == 2) || (i_r == 2 && i_f == 0)){continue;}

	  if (result_TSV[i_m][0][i_s][i_r][i_f] < result_min_seven_point){
	    result_min_seven_point = result_TSV[i_m][0][i_s][i_r][i_f];
	    x_r_min_seven_point = i_r;
	    x_f_min_seven_point = i_f;
	  }
	  if (result_TSV[i_m][0][i_s][i_r][i_f] > result_max_seven_point){
	    result_max_seven_point = result_TSV[i_m][0][i_s][i_r][i_f];
	    x_r_max_seven_point = i_r;
	    x_f_max_seven_point = i_f;
	  }
  	}
      }

      outfile_seven_point << left << setw(15) << "#" << noshowpoint
			  << setw(30) << "central scale"
			  << setw(30) << "minimum (7-point)"
			  << setw(30) << "maximum (7-point)"
			  << endl;
      
      outfile_nine_point << left << setw(15) << "#" << noshowpoint
			 << setw(30) << "central scale"
			 << setw(30) << "minimum (9-point)"
			 << setw(30) << "maximum (9-point)"
			 << endl;
	
      outfile_seven_point << left << setw(15) << "#" << noshowpoint
			  << "mu_R = " 
			  << setprecision(2) << setw(4) << showpoint << osi->relative_scale_ren_TSV[i_s][1] 
			  << " -- " 
			  << "mu_F = " 
			  << setw(4) << showpoint << osi->relative_scale_fact_TSV[i_s][1] << "    "
			  << "mu_R = " 
			  << setprecision(2) << setw(4) << showpoint << osi->relative_scale_ren_TSV[i_s][x_r_min_seven_point] 
			  << " -- " 
			  << "mu_F = " 
			  << setw(4) << showpoint << osi->relative_scale_fact_TSV[i_s][x_f_min_seven_point] << "    "
			  << "mu_R = " 
			  << setprecision(2) << setw(4) << showpoint << osi->relative_scale_ren_TSV[i_s][x_r_max_seven_point] 
			  << " -- " 
			  << "mu_F = " 
			  << setw(4) << showpoint << osi->relative_scale_fact_TSV[i_s][x_f_max_seven_point] << "    "
			  << endl;
      
      outfile_nine_point << left << setw(15) << "#" << noshowpoint
			 << "mu_R = " 
			 << setprecision(2) << setw(4) << showpoint << osi->relative_scale_ren_TSV[i_s][1] 
			 << " -- " 
			 << "mu_F = " 
			 << setw(4) << showpoint << osi->relative_scale_fact_TSV[i_s][1] << "    "
			 << "mu_R = " 
			 << setprecision(2) << setw(4) << showpoint << osi->relative_scale_ren_TSV[i_s][x_r_min_nine_point] 
			 << " -- " 
			 << "mu_F = " 
			 << setw(4) << showpoint << osi->relative_scale_fact_TSV[i_s][x_f_min_nine_point] << "    "
			 << "mu_R = " 
			 << setprecision(2) << setw(4) << showpoint << osi->relative_scale_ren_TSV[i_s][x_r_max_nine_point] 
			 << " -- " 
			 << "mu_F = " 
			 << setw(4) << showpoint << osi->relative_scale_fact_TSV[i_s][x_f_max_nine_point] << "    "
			 << endl;
      
      outfile_complete << endl;
      outfile_ren << endl;
      outfile_fact << endl;
      outfile_equal << endl;
      outfile_antipodal << endl;
      
      //  write out only qTcut = 0 value for qTcut-independent results and extrapolated value for qTcut-dependent result:

      logger << LOG_DEBUG_VERBOSE << "osi->n_qTcut - 1 = " << osi->n_qTcut - 1 << endl;
      for (int i = 0; i < 2; i++){
	logger << LOG_DEBUG_VERBOSE << "i = " << i << endl;

	if (i == 0){
      	  outfile_complete << setw(15) << "0";
	  outfile_equal << setw(15) << "0";
	  outfile_antipodal << setw(15) << "0";
	  outfile_ren << setw(15) << "0";
	  outfile_fact << setw(15) << "0";
	  outfile_seven_point << setw(15) << "0";
	  outfile_nine_point << setw(15) << "0";
	}
	else {
	  logger << LOG_DEBUG_VERBOSE << "osi->value_qTcut.size() = " << osi->value_qTcut.size() << endl;

	  outfile_complete << noshowpoint << setw(15) << osi->value_qTcut[osi->n_qTcut - 1];
	  outfile_equal << noshowpoint << setw(15) << osi->value_qTcut[osi->n_qTcut - 1];
	  outfile_antipodal << noshowpoint << setw(15) << osi->value_qTcut[osi->n_qTcut - 1];
	  outfile_ren << noshowpoint << setw(15) << osi->value_qTcut[osi->n_qTcut - 1];
	  outfile_fact << noshowpoint << setw(15) << osi->value_qTcut[osi->n_qTcut - 1];
	  outfile_seven_point << noshowpoint << setw(15) << osi->value_qTcut[osi->n_qTcut - 1];
	  outfile_nine_point << noshowpoint << setw(15) << osi->value_qTcut[osi->n_qTcut - 1];
	}
	logger << LOG_DEBUG_VERBOSE << "i = " << i << endl;
	  	
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    stringstream temp_ss;
	    temp_ss << left << showpoint 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_s][i_r][i_f] 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][i_r][i_f];
	    
	    outfile_complete << temp_ss.str();
	    if (i_f == osi->no_central_scale_fact_TSV[i_s]){outfile_ren << temp_ss.str();}
	    if (i_r == osi->no_central_scale_ren_TSV[i_s]){outfile_fact << temp_ss.str();}
	    if (i_r == i_f){outfile_equal << temp_ss.str();}
	    if (i_r == osi->n_scale_fact_TSV[i_s] - i_f - 1){outfile_antipodal << temp_ss.str();}
	  }
	}
	logger << LOG_DEBUG_VERBOSE << "i = " << i << endl;

	outfile_seven_point << left << showpoint 
			    << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_s][1][1] 
			    << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][1][1]
			    << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_s][x_r_min_seven_point][x_f_min_seven_point] 
			    << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][x_r_min_seven_point][x_f_min_seven_point] 
			    << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_s][x_r_max_seven_point][x_f_max_seven_point] 
			    << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][x_r_max_seven_point][x_f_max_seven_point];
	
	outfile_nine_point << left << showpoint 
			   << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_s][1][1] 
			   << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][1][1]
			   << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_s][x_r_min_nine_point][x_f_min_nine_point] 
			   << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][x_r_min_nine_point][x_f_min_nine_point] 
			   << setprecision(8) << setw(15) << osi->unit_factor_result * result_TSV[i_m][0][i_s][x_r_max_nine_point][x_f_max_nine_point] 
			   << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_TSV[i_m][0][i_s][x_r_max_nine_point][x_f_max_nine_point];
	logger << LOG_DEBUG_VERBOSE << "i = " << i << endl;
	
	outfile_complete << endl;
	outfile_ren << endl;
	outfile_fact << endl;
	outfile_equal << endl;
	outfile_antipodal << endl;
	outfile_seven_point << endl;
	outfile_nine_point << endl;
      }
      
      outfile_complete.close();
      outfile_ren.close();
      outfile_fact.close();
      outfile_equal.close();
      outfile_antipodal.close();
      outfile_seven_point.close();
      outfile_nine_point.close();
    }
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}



void summary_order::output_result_plot_qTcut_TSV(){
  Logger logger("summary_order::output_result_plot_qTcut_TSV");
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
      system_execute(logger, "mkdir " + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s]);

      filename_complete = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/complete/plot.qTcut" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_ren = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/ren/plot.qTcut" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_fact = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/fact/plot.qTcut" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_equal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/equal/plot.qTcut" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
      filename_antipodal = "" + ygeneric->final_resultdirectory + "/" + osi->name_extended_set_TSV[i_s] + "/antipodal/plot.qTcut" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";

      outfile_complete.open(filename_complete.c_str(), ofstream::out | ofstream::trunc);  
      outfile_ren.open(filename_ren.c_str(), ofstream::out | ofstream::trunc);  
      outfile_fact.open(filename_fact.c_str(), ofstream::out | ofstream::trunc);  
      outfile_equal.open(filename_equal.c_str(), ofstream::out | ofstream::trunc);  
      outfile_antipodal.open(filename_antipodal.c_str(), ofstream::out | ofstream::trunc);  
      
      outfile_complete << left << setw(15) << "#";
      outfile_ren << left << setw(15) << "#";
      outfile_fact << left << setw(15) << "#";
      outfile_equal << left << setw(15) << "#";
      outfile_antipodal << left << setw(15) << "#";

      for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	  stringstream temp_ss;
	  temp_ss << "mu_R = " 
		  << setprecision(2) << setw(4) << showpoint << osi->relative_scale_ren_TSV[i_s][i_r] 
		  << " -- " 
		  << "mu_F = " 
		  << setw(4) << showpoint << osi->relative_scale_fact_TSV[i_s][i_f] << "    ";

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
      
      if (!(active_qTcut)){
	outfile_complete << setw(15) << "0";
	outfile_equal << setw(15) << "0";
	outfile_antipodal << setw(15) << "0";
	outfile_ren << setw(15) << "0";
	outfile_fact << setw(15) << "0";
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    stringstream temp_ss;
	    temp_ss << left << showpoint 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * result_qTcut_TSV[i_m][0][0][i_s][i_r][i_f] 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_qTcut_TSV[i_m][0][0][i_s][i_r][i_f];

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

      for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	int x_q = 0;
	if (active_qTcut){x_q = i_q;}
	outfile_complete << noshowpoint << setw(15) << osi->value_qTcut[i_q];
	outfile_equal << noshowpoint << setw(15) << osi->value_qTcut[i_q];
	outfile_antipodal << noshowpoint << setw(15) << osi->value_qTcut[i_q];
	outfile_ren << noshowpoint << setw(15) << osi->value_qTcut[i_q];
	outfile_fact << noshowpoint << setw(15) << osi->value_qTcut[i_q];
	for (int i_r = 0; i_r < osi->n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < osi->n_scale_fact_TSV[i_s]; i_f++){
	    stringstream temp_ss;
	    temp_ss << left << showpoint 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * result_qTcut_TSV[i_m][0][x_q][i_s][i_r][i_f] 
		    << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_qTcut_TSV[i_m][0][x_q][i_s][i_r][i_f];

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
      outfile_complete.close();
      outfile_ren.close();
      outfile_fact.close();
      outfile_equal.close();
      outfile_antipodal.close();
    }
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}



void summary_order::output_result_plot_CV(){
  Logger logger("summary_order::output_result_plot_CV");
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

  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
    ////////////////////////////////////
    //  output of usual CV variation  //
    ////////////////////////////////////
    ofstream outfile;
    filename = ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/plot" + ygeneric->infix_name_moment[i_m] + ".CV." + resultdirectory + ".dat";
    logger << LOG_DEBUG_VERBOSE << "FILENAME:   " << filename << endl;
    outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    int i_q = 0;
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      // rel_scale_CV does not containt reasonable scale information !!!
      outfile << setw(15) << setprecision(8) << rel_scale_CV[i_s] 
	      << setprecision(8) << setw(15) << osi->unit_factor_result * result_qTcut_CV[i_m][0][i_q][i_s] 
	      << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_qTcut_CV[i_m][0][i_q][i_s] << endl;
    }
    outfile.close();

    //////////////////////////////////////////////////////////////////
    //  output of extrapolated (or qTcut-independent) CV variation  //
    //////////////////////////////////////////////////////////////////
    filename = ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/plot" + ygeneric->infix_name_moment[i_m] + ".CV.extrapolated." + resultdirectory + ".dat";
    logger << LOG_DEBUG_VERBOSE << "FILENAME:   " << filename << endl;
    outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      // rel_scale_CV does not containt reasonable scale information !!!
      outfile << setw(15) << setprecision(8) << rel_scale_CV[i_s] 
	      << setprecision(8) << setw(15) << osi->unit_factor_result * result_CV[i_m][0][i_s] 
	      << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_CV[i_m][0][i_s] 
	      << endl;
      if (i_m == osi->n_moments){
	logger << LOG_DEBUG_VERBOSE << setw(15) << setprecision(8) << rel_scale_CV[i_s] << setprecision(8) << setw(15) << "result_CV[" << i_m << "][0][" << i_s << "] = " << result_CV[i_m][0][i_s] << setprecision(8) << setw(15) << deviation_CV[i_m][0][i_s] << endl;
      }
    }
    outfile.close();

    ///////////////////////////////////////////////////////////////////////////
    //  output of extrapolated (or qTcut-independent) result in plot format  //
    ///////////////////////////////////////////////////////////////////////////
    filename = ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/plot" + ygeneric->infix_name_moment[i_m] + "." + resultdirectory + ".dat";
    logger << LOG_DEBUG_VERBOSE << "xfilename i_m = " << i_m << "   " << filename << endl;
    outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  
    outfile << showpoint << left 
	    << setw(15) << "# qTcut";
    for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
      stringstream temp_ss;
      temp_ss << "mu_R = " 
	      << setprecision(2) << setw(4) << showpoint << osi->rel_scale_factor_ren_CV[i_s] 
	      << " -- " 
	      << "mu_F = " 
	      << setprecision(2) << setw(4) << showpoint << osi->rel_scale_factor_fact_CV[i_s] << "    ";
      outfile << setw(30) << temp_ss.str() << left;
    }
    outfile << endl;

    for (int i = 0; i < 2; i++){
      if (i == 0){
	outfile << setw(15) << "0";
      }
      else {
	outfile << noshowpoint << setw(15) << osi->value_qTcut[osi->n_qTcut - 1];
      }
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	outfile << showpoint 
		<< setprecision(8) << setw(15) << osi->unit_factor_result * result_CV[i_m][0][i_s] 
		<< setprecision(8) << setw(15) << osi->unit_factor_result * deviation_CV[i_m][0][i_s];
      }
      outfile << endl;
    }
    outfile.close();
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_order::output_result_plot_qTcut_CV(){
  Logger logger("summary_order::output_result_plot_qTcut_CV");
  logger << LOG_DEBUG << "called" << endl;

  ofstream outfile;
  string filename;

  if (osi->n_qTcut > 1){
    for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
      filename = ygeneric->final_resultdirectory + "/CV/" + ygeneric->name_variation_CV + "/plot" + ygeneric->infix_name_moment[i_m] + ".qTcut." + resultdirectory + ".dat";      outfile.open(filename.c_str(), ofstream::out | ofstream::trunc);  
      
      outfile << showpoint << left 
	      << setw(15) << "# qTcut";
      for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	stringstream temp_ss;
	temp_ss << "mu_R = " 
		<< setprecision(2) << setw(4) << showpoint << osi->rel_scale_factor_ren_CV[i_s] 
		<< " -- " 
		<< "mu_F = " 
		<< setprecision(2) << setw(4) << showpoint << osi->rel_scale_factor_fact_CV[i_s] << "    ";

	outfile << setw(30) << temp_ss.str() << left;
      }
      outfile << endl;

      if (result_qTcut_CV[i_m][0][0][0] == result_qTcut_CV[i_m][0][osi->n_qTcut - 1][0]){
	outfile << setw(15) << "0.";
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  outfile << showpoint 
		  << setprecision(8) << setw(15) << osi->unit_factor_result * result_qTcut_CV[i_m][0][0][i_s] 
		  << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_qTcut_CV[i_m][0][0][i_s];
	}
	outfile << endl;
      }
      
      for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	outfile << noshowpoint << setw(15) << osi->value_qTcut[i_q];
	for (int i_s = 0; i_s < osi->n_scales_CV; i_s++){
	  outfile << showpoint 
		  << setprecision(8) << setw(15) << osi->unit_factor_result * result_qTcut_CV[i_m][0][i_q][i_s] 
		  << setprecision(8) << setw(15) << osi->unit_factor_result * deviation_qTcut_CV[i_m][0][i_q][i_s];
	}
	outfile << endl;
      }
      outfile.close();
      logger << LOG_DEBUG_VERBOSE << filename << " closed." << endl;
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}
