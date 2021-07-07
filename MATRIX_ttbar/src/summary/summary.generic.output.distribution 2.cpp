#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void summary_generic::output_distribution_table_order(){
  Logger logger("summary_generic::get_summary");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

    // create table output here (without K-factors etc.) !!!

  for (int i_sb = 0; i_sb < outpath_scaleband.size(); i_sb++){
    for (int x_q = 0; x_q < osi_value_qTcut_distribution.size() + 1; x_q++){
      // only copy&paste: more elegant soultion for directory_qTcut !!!
      string directory_qTcut;
      string in_directory_qTcut;
      string out_directory_qTcut;
      if (x_q == osi_value_qTcut_distribution.size()){
	directory_qTcut = "";
	in_directory_qTcut = "";
	out_directory_qTcut = "";
      }
      else {
	stringstream qTcut_ss;
	qTcut_ss << "qTcut-" << osi_value_qTcut_distribution[x_q];
	directory_qTcut = "/" + qTcut_ss.str();
	in_directory_qTcut = "";
	out_directory_qTcut = directory_qTcut;
      }
      
      logger << LOG_DEBUG << "x_q = " << x_q << "   directory_qTcut = " << directory_qTcut << endl;
      for (int i_d = 0; i_d < oset.extended_distribution.size(); i_d++){
	if (!switch_output_distribution_table[i_d]){continue;}

	vector<string> name_plot;
	string temp_sdd = "." + oset.extended_distribution[i_d].xdistribution_name;
	
	name_plot.push_back("table.plot" + temp_sdd + ".dat");
	name_plot.push_back("table.norm" + temp_sdd + ".dat");
	
	if (i_d >= osi_dat.size()){
	  int i_ddd = i_d - osi_dat.size();
	  
	  // recombined distributions (essentially for validation)
	  stringstream name_rec;
	  name_rec << ".rec." << osi_dddat[i_ddd].distribution_2.xdistribution_name << ".from." << osi_dddat[i_ddd].name;
	  //	  string s0;
	  name_plot.push_back("table.plot" + name_rec.str() + ".tex");
	  name_plot.push_back("table.norm" + name_rec.str()  + ".tex");
	  
	  for (int i_b1 = 0; i_b1 < osi_dddat[i_ddd].distribution_1.n_bins; i_b1++){
	    stringstream name_split_bin;
	    name_split_bin << "_" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1] << "-" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1 + 1];
	    stringstream name_split_lt;
	    name_split_lt << "_lt" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1];
	    stringstream name_split_ge;
	    name_split_ge << "_ge" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1];
	    
	    name_plot.push_back("table.norm.norm.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("table.norm.plot.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("table.plot.norm.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("table.plot.plot.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    
	    name_plot.push_back("table.norm.two." + osi_dddat[i_ddd].name + name_split_lt.str() + ".tex");
	    name_plot.push_back("table.norm.two." + osi_dddat[i_ddd].name + name_split_ge.str() + ".tex");
	    name_plot.push_back("table.plot.two." + osi_dddat[i_ddd].name + name_split_lt.str() + ".tex");
	    name_plot.push_back("table.plot.two." + osi_dddat[i_ddd].name + name_split_ge.str() + ".tex");
	  }
	}
	


	
	vector<string> output_order(4);
	output_order[0] = "LO";
	output_order[1] = "NLO.QCD";
	output_order[2] = "nNLO.QCD+gg";
	output_order[3] = "NNLO.QCD";

	vector<int> no_output_order(output_order.size(), -1);
	for (int j_o = 0; j_o < output_order.size(); j_o++){ 
	  for (int i_o = 0; i_o < yorder.size(); i_o++){
	    if (output_order[j_o] == yorder[i_o].resultdirectory){no_output_order[j_o] = i_o; break;}
	  }
	}

	for (int j_o = 0; j_o < output_order.size(); j_o++){ 
	  logger << LOG_INFO << left
		 << setw(4) << j_o
		 << setw(20) << output_order[j_o]
		 << " -> "
		 << setw(4) << no_output_order[j_o] << endl;
	}

	int n_output_version = name_plot.size();
	for (int i_x = 0; i_x < n_output_version; i_x++){
	  logger << LOG_INFO << "i_x = " << i_x << endl;

	  string outfilename_latex_table = final_resultdirectory + "/" + outpath_scaleband[i_sb] + directory_qTcut + "/" + name_plot[i_x];

	  logger << LOG_INFO << "outfilename_latex_table = " <<  outfilename_latex_table << endl;

	  ofstream outf;
	  outf.open(outfilename_latex_table.c_str(), ofstream::out | ofstream::trunc);

	  outf << char(92) << "renewcommand" << char(92) << "arraystretch{1.5}" << endl;
	  outf << char(92) << "begin{table}" << endl;
	  outf << char(92) << "begin{center}" << endl;
	  outf << char(92) << "begin{tabular}{|c|";
	  for (int j_o = 0; j_o < output_order.size(); j_o++){
	    int i_o = no_output_order[j_o];
	    if (i_o == -1){continue;}
	    outf << "c|";
	  }
	  outf << "}" << endl;
	  outf << char(92) << "hline" << endl;
	  outf << "$" << char(92) << "sqrt{s}" << char(92) << ", (" << char(92) << "mathrm{TeV})$ &" << endl;
	  for (int j_o = 0; j_o < output_order.size(); j_o++){
	    int i_o = no_output_order[j_o];
	    if (i_o == -1){continue;}

	    if (output_order[j_o] == "LO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{LO}}";}
	    else if (output_order[j_o] == "nLO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nLO}}";}
	    else if (output_order[j_o] == "nnLO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nnLO}}";}
	    else if (output_order[j_o] == "NLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{NLO}}";}
	    else if (output_order[j_o] == "nNLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nNLO}}";}
	    else if (output_order[j_o] == "nNLO.QCD+gg"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nNLO+gg}}";}
	    else if (output_order[j_o] == "NNLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{NNLO}}";}
	    else {outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_order[j_o] << "}}$";}
	    outf << char(92) << ", (" << char(92) << "mathrm{" << oset.unit_distribution << "})$";
	    
	    if (j_o == output_order.size() - 1){outf << " " << char(92) << char(92) << endl;}
	    else {outf << " &" << endl;}
	  }
	  outf << char(92) << "hline" << endl;

	  // better: variable... -> doubled lines for plotting reasons etc.
	  /*
	  for (int i_b = 0; i_b < oset.extended_distribution[i_d].n_bins; i_b++){
	    outf << oset.extended_distribution[i_d].bin_edge[i_b] << " &" << endl;
	  */
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "].size() = " << scaleband_variable[i_sb].size() << endl;
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "][" << i_d << "].size() = " << scaleband_variable[i_sb][i_d].size() << endl;
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "][" << i_d << "][" << i_x << "].size() = " << scaleband_variable[i_sb][i_d][i_x].size() << endl;

	  for (int i_b = 0; i_b < scaleband_variable[i_sb][i_d][i_x].size() - 1; i_b++){
	    //	  for (int i_b = 0; i_b < scaleband_variable[i_sb][i_d][i_x].size(); i_b++){
	    outf << oset.E_CMS / 1000 << " " << char(92) << "nameobservable{" << scaleband_variable[i_sb][i_d][i_x][i_b] << "} &" << endl;
	    
	    for (int j_o = 0; j_o < output_order.size(); j_o++){
	      int i_o = no_output_order[j_o];
	      if (i_o == -1){continue;}

	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "].size() = " << scaleband_central_deviation[i_sb].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "].size() = " << scaleband_central_deviation[i_sb][i_d].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "][" << i_o << "].size() = " << scaleband_central_deviation[i_sb][i_d][i_o].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "][" << i_o << "][" << i_x << "].size() = " << scaleband_central_deviation[i_sb][i_d][i_o][i_x].size() << endl;
	      
	      outf << "$" << output_result_deviation(scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], scaleband_central_deviation[i_sb][i_d][i_o][i_x][i_b], 1)
		   << "^{" << output_latex_percent(scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], scaleband_maximum_result[i_sb][i_d][i_o][i_x][i_b], 1) << "}"
		   << "_{" << output_latex_percent(scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], scaleband_minimum_result[i_sb][i_d][i_o][i_x][i_b], 1) << "}$";
	      /*
	      outf << setw(16) << setprecision(8) << scaleband_central_result[i_sb][i_d][i_o][i_x][i_b] << "("
		   << setw(16) << setprecision(8) << scaleband_central_deviation[i_sb][i_d][i_o][i_x][i_b] << ")"
		   << "^{" << (scaleband_maximum_result[i_sb][i_d][i_o][i_x][i_b] / scaleband_central_result[i_sb][i_d][i_o][i_x][i_b] - 1.) * 100. << char(92) << "%}"
		   << "_{" << (scaleband_minimum_result[i_sb][i_d][i_o][i_x][i_b] / scaleband_central_result[i_sb][i_d][i_o][i_x][i_b] - 1.) * 100. << char(92) << "%}";
	      */	      
	      if (j_o == output_order.size() - 1){
		outf << " " << char(92) << char(92) << endl;
		outf << char(92) << "hline" << endl;
	      }
	      else {outf << " &" << endl;}
	      //  if (i < minresult.size() - 1){outf << setw(16) << setprecision(8) << variable[i + 1] << setw(16) << setprecision(8) << minresult[i_b] << "   " << setw(16) << setprecision(8) << mindeviation[i_b] << endl;}
	    }
	  }
	  outf << char(92) << "end{tabular}" << endl;
	  outf << char(92) << "end{center}" << endl;
	  outf << char(92) << "caption{" << char(92) << "captiontext}" << endl;
	  outf << char(92) << "end{table}" << endl;
	  outf.close();
	}
      }
    }
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}



void summary_generic::output_distribution_table_Kfactor(){
  Logger logger("summary_generic::get_summary");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_sb = 0; i_sb < outpath_scaleband.size(); i_sb++){
    // create K-factor table output here !!!


    for (int i_y = 0; i_y < 2; i_y++){
      string output_reference;
      string output_reference_latex;
      vector<string> output_order;

      if (i_y == 0){
	output_reference = "NLO.QCD";
	output_reference_latex = "NLO";
	
	output_order.resize(5);
	output_order[0] = "LO";
	output_order[1] = "NLO.QCD";
	output_order[2] = "nNLO.QCD";
	output_order[3] = "nNLO.QCD+gg";
	output_order[4] = "NNLO.QCD";
      }
      
      else if (i_y == 1){
	output_reference = "nNLO.QCD";
	output_reference_latex = "nNLO";
	
	output_order.resize(4);
	output_order[0] = "nnLO";
	output_order[1] = "nNLO.QCD";
	output_order[2] = "nNLO.QCD+gg";
	output_order[3] = "NNLO.QCD";
      }
      


    for (int x_q = 0; x_q < osi_value_qTcut_distribution.size() + 1; x_q++){
      // only copy&paste: more elegant soultion for directory_qTcut !!!
      string directory_qTcut;
      string in_directory_qTcut;
      string out_directory_qTcut;
      if (x_q == osi_value_qTcut_distribution.size()){
	directory_qTcut = "";
	in_directory_qTcut = "";
	out_directory_qTcut = "";
      }
      else {
	stringstream qTcut_ss;
	qTcut_ss << "qTcut-" << osi_value_qTcut_distribution[x_q];
	directory_qTcut = "/" + qTcut_ss.str();
	in_directory_qTcut = "";
	out_directory_qTcut = directory_qTcut;
      }
      
      logger << LOG_DEBUG << "x_q = " << x_q << "   directory_qTcut = " << directory_qTcut << endl;
      for (int i_d = 0; i_d < oset.extended_distribution.size(); i_d++){
	if (!switch_output_distribution_table[i_d]){continue;}
	
	vector<string> name_plot;
	string temp_sdd = "." + oset.extended_distribution[i_d].xdistribution_name;
	
	name_plot.push_back("Kfactor" + output_reference_latex + ".plot" + temp_sdd + ".dat");
	name_plot.push_back("Kfactor" + output_reference_latex + ".norm" + temp_sdd + ".dat");
	
	if (i_d >= osi_dat.size()){
	  int i_ddd = i_d - osi_dat.size();
	  
	// recombined distributions (essentially for validation)
	  stringstream name_rec;
	  name_rec << ".rec." << osi_dddat[i_ddd].distribution_2.xdistribution_name << ".from." << osi_dddat[i_ddd].name;
	  //	  string s0;
	  name_plot.push_back("Kfactor" + output_reference_latex + ".plot" + name_rec.str() + ".tex");
	  name_plot.push_back("Kfactor" + output_reference_latex + ".norm" + name_rec.str()  + ".tex");
	  
	  for (int i_b1 = 0; i_b1 < osi_dddat[i_ddd].distribution_1.n_bins; i_b1++){
	    stringstream name_split_bin;
	    name_split_bin << "_" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1] << "-" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1 + 1];
	    stringstream name_split_lt;
	    name_split_lt << "_lt" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1];
	    stringstream name_split_ge;
	    name_split_ge << "_ge" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1];
	    
	    name_plot.push_back("Kfactor" + output_reference_latex + ".norm.norm.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("Kfactor" + output_reference_latex + ".norm.plot.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("Kfactor" + output_reference_latex + ".plot.norm.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("Kfactor" + output_reference_latex + ".plot.plot.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    
	    name_plot.push_back("Kfactor" + output_reference_latex + ".norm.two." + osi_dddat[i_ddd].name + name_split_lt.str() + ".tex");
	    name_plot.push_back("Kfactor" + output_reference_latex + ".norm.two." + osi_dddat[i_ddd].name + name_split_ge.str() + ".tex");
	    name_plot.push_back("Kfactor" + output_reference_latex + ".plot.two." + osi_dddat[i_ddd].name + name_split_lt.str() + ".tex");
	    name_plot.push_back("Kfactor" + output_reference_latex + ".plot.two." + osi_dddat[i_ddd].name + name_split_ge.str() + ".tex");
	  }
	}
	


	int no_ref = -1;
	vector<int> no_output_order(output_order.size(), -1);
	for (int j_o = 0; j_o < output_order.size(); j_o++){ 
	  for (int i_o = 0; i_o < yorder.size(); i_o++){
	    if (output_order[j_o] == yorder[i_o].resultdirectory){no_output_order[j_o] = i_o; break;}
	  }
	  if (output_reference == output_order[j_o]){no_ref = no_output_order[j_o];}
	}

	for (int j_o = 0; j_o < output_order.size(); j_o++){ 
	  logger << LOG_INFO << left
		 << setw(4) << j_o
		 << setw(20) << output_order[j_o]
		 << " -> "
		 << setw(4) << no_output_order[j_o] << endl;
	}

	int n_output_version = name_plot.size();
	for (int i_x = 0; i_x < n_output_version; i_x++){
	  logger << LOG_INFO << "i_x = " << i_x << endl;

	  string outfilename_latex_Kfactor = final_resultdirectory + "/" + outpath_scaleband[i_sb] + directory_qTcut + "/" + name_plot[i_x];

	  logger << LOG_INFO << "outfilename_latex_Kfactor = " <<  outfilename_latex_Kfactor << endl;

	  ofstream outf;
	  outf.open(outfilename_latex_Kfactor.c_str(), ofstream::out | ofstream::trunc);

	  outf << char(92) << "renewcommand" << char(92) << "arraystretch{1.5}" << endl;
	  outf << char(92) << "begin{table}" << endl;
	  outf << char(92) << "begin{center}" << endl;
	  outf << char(92) << "begin{tabular}{|c|";
	  for (int j_o = 0; j_o < output_order.size(); j_o++){
	    int i_o = no_output_order[j_o];
	    if (i_o == -1){continue;}
	    if (output_order[j_o] != output_reference){outf << "c|";}
	  }
	  outf << "}" << endl;
	  outf << char(92) << "hline" << endl;
	  outf << "$" << char(92) << "sqrt{s}" << char(92) << ", (" << char(92) << "mathrm{TeV})$ &" << endl;
	  for (int j_o = 0; j_o < output_order.size(); j_o++){
	    int i_o = no_output_order[j_o];
	    if (i_o == -1){continue;}

	    if (output_order[j_o] == output_reference){
	      if (j_o == output_order.size() - 1){outf << " " << char(92) << char(92) << endl;}

	    }
	    else {
	      if (output_order[j_o] == "LO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{LO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "nLO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nLO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "nnLO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nnLO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "NLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{NLO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "nNLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nNLO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "nNLO.QCD+gg"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nNLO+gg}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "NNLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{NNLO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else {outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_order[j_o] << "}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      //	      outf << char(92) << ", (" << char(92) << "mathrm{" << oset.unit_distribution << "})$";
	      outf << "$";
	      if (j_o == output_order.size() - 1){outf << " " << char(92) << char(92) << endl;}
	      else {outf << " &" << endl;}
	    }
	  
	  }
	  outf << char(92) << "hline" << endl;

	  // better: variable... -> doubled lines for plotting reasons etc.
	  /*
	  for (int i_b = 0; i_b < oset.extended_distribution[i_d].n_bins; i_b++){
	    outf << oset.extended_distribution[i_d].bin_edge[i_b] << " &" << endl;
	  */
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "].size() = " << scaleband_variable[i_sb].size() << endl;
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "][" << i_d << "].size() = " << scaleband_variable[i_sb][i_d].size() << endl;
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "][" << i_d << "][" << i_x << "].size() = " << scaleband_variable[i_sb][i_d][i_x].size() << endl;

	  for (int i_b = 0; i_b < scaleband_variable[i_sb][i_d][i_x].size() - 1; i_b++){
	    //	  for (int i_b = 0; i_b < scaleband_variable[i_sb][i_d][i_x].size(); i_b++){
	    outf << oset.E_CMS / 1000 << " " << char(92) << "nameobservable{" << scaleband_variable[i_sb][i_d][i_x][i_b] << "} &" << endl;
	    
	    for (int j_o = 0; j_o < output_order.size(); j_o++){
	      int i_o = no_output_order[j_o];
	      if (i_o == -1){continue;}

	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "].size() = " << scaleband_central_deviation[i_sb].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "].size() = " << scaleband_central_deviation[i_sb][i_d].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "][" << i_o << "].size() = " << scaleband_central_deviation[i_sb][i_d][i_o].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "][" << i_o << "][" << i_x << "].size() = " << scaleband_central_deviation[i_sb][i_d][i_o][i_x].size() << endl;
	      
	      if (output_order[j_o] == output_reference){
		//		outf << "$0$";
		if (j_o == output_order.size() - 1){
		  outf << " " << char(92) << char(92) << endl;
		  outf << char(92) << "hline" << endl;
		}
	      }
	      else {
		//		double temp_double = scaleband_central_result[i_sb][i_d][no_ref][i_x][i_b] + scaleband_central_result[i_sb][i_d][i_o][i_x][i_b];
		outf << "$" << output_latex_percent(scaleband_central_result[i_sb][i_d][no_ref][i_x][i_b], scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], 1) << "$";
	      
		if (j_o == output_order.size() - 1){
		  outf << " " << char(92) << char(92) << endl;
		  outf << char(92) << "hline" << endl;
		}
		else {outf << " &" << endl;}
		//  if (i < minresult.size() - 1){outf << setw(16) << setprecision(8) << variable[i + 1] << setw(16) << setprecision(8) << minresult[i_b] << "   " << setw(16) << setprecision(8) << mindeviation[i_b] << endl;}
	      }
	    }
	  }
	  outf << char(92) << "end{tabular}" << endl;
	  outf << char(92) << "end{center}" << endl;
	  outf << char(92) << "caption{" << char(92) << "captiontext}" << endl;
	  outf << char(92) << "end{table}" << endl;
	  outf.close();
	}
      }
    }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_generic::output_distribution_table_crosssection_Kfactor(){
  Logger logger("summary_generic::get_summary");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_sb = 0; i_sb < outpath_scaleband.size(); i_sb++){
    // create K-factor table output here !!!


    for (int i_y = 0; i_y < 2; i_y++){
      string output_reference;
      string output_reference2;
      string output_reference_latex;
      vector<string> output_order;

      if (i_y == 0){
	output_reference = "NLO.QCD";
	output_reference2 = "LO";
	output_reference_latex = "NLO";
	
	output_order.resize(4);
	//	output_order.resize(5);
	output_order[0] = "LO";
	output_order[1] = "NLO.QCD";
	output_order[2] = "nNLO.QCD+gg";
	output_order[3] = "NNLO.QCD";
	//	output_order[2] = "nNLO.QCD";
	//	output_order[3] = "nNLO.QCD+gg";
	//	output_order[4] = "NNLO.QCD";
      }
      
      else if (i_y == 1){
	output_reference = "nNLO.QCD";
	output_reference2 = "nnLO";
	output_reference_latex = "nNLO";
	
	output_order.resize(4);
	output_order[0] = "nnLO";
	output_order[1] = "nNLO.QCD";
	output_order[2] = "nNLO.QCD+gg";
	output_order[3] = "NNLO.QCD";
      }
      


    for (int x_q = 0; x_q < osi_value_qTcut_distribution.size() + 1; x_q++){
      // only copy&paste: more elegant soultion for directory_qTcut !!!
      string directory_qTcut;
      string in_directory_qTcut;
      string out_directory_qTcut;
      if (x_q == osi_value_qTcut_distribution.size()){
	directory_qTcut = "";
	in_directory_qTcut = "";
	out_directory_qTcut = "";
      }
      else {
	stringstream qTcut_ss;
	qTcut_ss << "qTcut-" << osi_value_qTcut_distribution[x_q];
	directory_qTcut = "/" + qTcut_ss.str();
	in_directory_qTcut = "";
	out_directory_qTcut = directory_qTcut;
      }
      
      logger << LOG_DEBUG << "x_q = " << x_q << "   directory_qTcut = " << directory_qTcut << endl;
      for (int i_d = 0; i_d < oset.extended_distribution.size(); i_d++){
	if (!switch_output_distribution_table[i_d]){continue;}
	
	vector<string> name_plot;
	string temp_sdd = "." + oset.extended_distribution[i_d].xdistribution_name;
	
	name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".plot" + temp_sdd + ".dat");
	name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".norm" + temp_sdd + ".dat");
	
	if (i_d >= osi_dat.size()){
	  int i_ddd = i_d - osi_dat.size();
	  
	// recombined distributions (essentially for validation)
	  stringstream name_rec;
	  name_rec << ".rec." << osi_dddat[i_ddd].distribution_2.xdistribution_name << ".from." << osi_dddat[i_ddd].name;
	  //	  string s0;
	  name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".plot" + name_rec.str() + ".tex");
	  name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".norm" + name_rec.str()  + ".tex");
	  
	  for (int i_b1 = 0; i_b1 < osi_dddat[i_ddd].distribution_1.n_bins; i_b1++){
	    stringstream name_split_bin;
	    name_split_bin << "_" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1] << "-" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1 + 1];
	    stringstream name_split_lt;
	    name_split_lt << "_lt" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1];
	    stringstream name_split_ge;
	    name_split_ge << "_ge" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1];
	    
	    name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".norm.norm.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".norm.plot.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".plot.norm.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".plot.plot.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    
	    name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".norm.two." + osi_dddat[i_ddd].name + name_split_lt.str() + ".tex");
	    name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".norm.two." + osi_dddat[i_ddd].name + name_split_ge.str() + ".tex");
	    name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".plot.two." + osi_dddat[i_ddd].name + name_split_lt.str() + ".tex");
	    name_plot.push_back("Xsection.Kfactor" + output_reference_latex + ".plot.two." + osi_dddat[i_ddd].name + name_split_ge.str() + ".tex");
	  }
	}
	


	int no_ref = -1;
	int no_ref2 = -1;
	vector<int> no_output_order(output_order.size(), -1);
	for (int j_o = 0; j_o < output_order.size(); j_o++){ 
	  for (int i_o = 0; i_o < yorder.size(); i_o++){
	    if (output_order[j_o] == yorder[i_o].resultdirectory){no_output_order[j_o] = i_o; break;}
	  }
	  if (output_reference == output_order[j_o]){no_ref = no_output_order[j_o];}
	  if (output_reference2 == output_order[j_o]){no_ref2 = no_output_order[j_o];}
	}

	for (int j_o = 0; j_o < output_order.size(); j_o++){ 
	  logger << LOG_INFO << left
		 << setw(4) << j_o
		 << setw(20) << output_order[j_o]
		 << " -> "
		 << setw(4) << no_output_order[j_o] << endl;
	}

	int n_output_version = name_plot.size();
	for (int i_x = 0; i_x < n_output_version; i_x++){
	  logger << LOG_INFO << "i_x = " << i_x << endl;

	  string outfilename_latex_Kfactor = final_resultdirectory + "/" + outpath_scaleband[i_sb] + directory_qTcut + "/" + name_plot[i_x];

	  logger << LOG_INFO << "outfilename_latex_Kfactor = " <<  outfilename_latex_Kfactor << endl;

	  ofstream outf;
	  outf.open(outfilename_latex_Kfactor.c_str(), ofstream::out | ofstream::trunc);

	  outf << char(92) << "renewcommand" << char(92) << "arraystretch{1.5}" << endl;
	  outf << char(92) << "begin{table}" << endl;
	  outf << char(92) << "begin{center}" << endl;
	  outf << char(92) << "begin{tabular}{|c|";
	  for (int j_o = 0; j_o < output_order.size(); j_o++){
	    int i_o = no_output_order[j_o];
	    if (i_o == -1){continue;}
	    //	    if (output_order[j_o] != output_reference)
	    else {outf << "c|";}
	  }
	  outf << "}" << endl;
	  outf << char(92) << "hline" << endl;
	  outf << "$" << char(92) << "sqrt{s}" << char(92) << ", (" << char(92) << "mathrm{TeV})$ &" << endl;

	  for (int j_o = 0; j_o < output_order.size(); j_o++){
	    int i_o = no_output_order[j_o];
	    if (i_o == -1){continue;}

	    if (output_order[j_o] == "LO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{LO}}";}
	    else if (output_order[j_o] == "nLO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nLO}}";}
	    else if (output_order[j_o] == "nnLO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nnLO}}";}
	    else if (output_order[j_o] == "NLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{NLO}}";}
	    else if (output_order[j_o] == "nNLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nNLO}}";}
	    else if (output_order[j_o] == "nNLO.QCD+gg"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nNLO+gg}}";}
	    else if (output_order[j_o] == "NNLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{NNLO}}";}
	    else {outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_order[j_o] << "}}$";}
	    outf << char(92) << ", (" << char(92) << "mathrm{" << oset.unit_distribution << "})$";
	    
	    if (j_o == output_order.size() - 1){outf << " " << char(92) << char(92) << endl;}
	    else {outf << " &" << endl;}
	  }

	  outf << " &" << endl;
	  for (int j_o = 0; j_o < output_order.size(); j_o++){
	    int i_o = no_output_order[j_o];
	    if (i_o == -1){continue;}

    	    if (output_order[j_o] == output_reference){
	      outf << "$" << char(92) << "frac{" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}}{" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference2 << "}}}-1";
	      outf << "$";
	      if (j_o == output_order.size() - 1){outf << " " << char(92) << char(92) << "[1ex]" << endl;}
	      else {outf << " &" << endl;}
	    }
	    else {
	      if (output_order[j_o] == "LO"){outf << "$" << char(92) << "frac{" << char(92) << "sigma_{" << char(92) << "mathrm{LO}}}{" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}}-1";}
	      else if (output_order[j_o] == "nLO"){outf << "$" << char(92) << "frac{" << char(92) << "sigma_{" << char(92) << "mathrm{nLO}}}{" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}}-1";}
	      else if (output_order[j_o] == "nnLO"){outf << "$" << char(92) << "frac{" << char(92) << "sigma_{" << char(92) << "mathrm{nnLO}}}{" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}}-1";}
	      else if (output_order[j_o] == "NLO.QCD"){outf << "$" << char(92) << "frac{" << char(92) << "sigma_{" << char(92) << "mathrm{NLO}}}{" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}}-1";}
	      else if (output_order[j_o] == "nNLO.QCD"){outf << "$" << char(92) << "frac{" << char(92) << "sigma_{" << char(92) << "mathrm{nNLO}}}{" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}}-1";}
	      else if (output_order[j_o] == "nNLO.QCD+gg"){outf << "$" << char(92) << "frac{" << char(92) << "sigma_{" << char(92) << "mathrm{nNLO+gg}}}{" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}}-1";}
	      else if (output_order[j_o] == "NNLO.QCD"){outf << "$" << char(92) << "frac{" << char(92) << "sigma_{" << char(92) << "mathrm{NNLO}}}{" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}}-1";}
	      else {outf << "$" << char(92) << "frac{" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_order[j_o] << "}}}{" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}}-1";}
	      //	      outf << char(92) << ", (" << char(92) << "mathrm{" << oset.unit_distribution << "})$";
	      /*
	      if (output_order[j_o] == "LO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{LO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "nLO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nLO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "nnLO"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nnLO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "NLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{NLO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "nNLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nNLO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "nNLO.QCD+gg"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{nNLO+gg}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else if (output_order[j_o] == "NNLO.QCD"){outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{NNLO}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      else {outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_order[j_o] << "}}/" << char(92) << "sigma_{" << char(92) << "mathrm{" << output_reference_latex << "}}-1";}
	      //	      outf << char(92) << ", (" << char(92) << "mathrm{" << oset.unit_distribution << "})$";
	      */
	      outf << "$";
	      if (j_o == output_order.size() - 1){outf << " " << char(92) << char(92) << "[1ex]" << endl;}
	      else {outf << " &" << endl;}
	    }
	  
	  }
	  outf << char(92) << "hline" << endl;

	  // better: variable... -> doubled lines for plotting reasons etc.
	  /*
	  for (int i_b = 0; i_b < oset.extended_distribution[i_d].n_bins; i_b++){
	    outf << oset.extended_distribution[i_d].bin_edge[i_b] << " &" << endl;
	  */
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "].size() = " << scaleband_variable[i_sb].size() << endl;
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "][" << i_d << "].size() = " << scaleband_variable[i_sb][i_d].size() << endl;
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "][" << i_d << "][" << i_x << "].size() = " << scaleband_variable[i_sb][i_d][i_x].size() << endl;

	  for (int i_b = 0; i_b < scaleband_variable[i_sb][i_d][i_x].size() - 1; i_b++){
	    //	  for (int i_b = 0; i_b < scaleband_variable[i_sb][i_d][i_x].size(); i_b++){
	    outf << oset.E_CMS / 1000 << " " << char(92) << "nameobservable{" << scaleband_variable[i_sb][i_d][i_x][i_b] << "} &" << endl;


	    for (int j_o = 0; j_o < output_order.size(); j_o++){
	      int i_o = no_output_order[j_o];
	      if (i_o == -1){continue;}

	      
	      outf << "$" << output_result_deviation(scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], scaleband_central_deviation[i_sb][i_d][i_o][i_x][i_b], 1)
		   << "^{" << output_latex_percent(scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], scaleband_maximum_result[i_sb][i_d][i_o][i_x][i_b], 1) << "}"
		   << "_{" << output_latex_percent(scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], scaleband_minimum_result[i_sb][i_d][i_o][i_x][i_b], 1) << "}$";
	      /*
	      outf << setw(16) << setprecision(8) << scaleband_central_result[i_sb][i_d][i_o][i_x][i_b] << "("
		   << setw(16) << setprecision(8) << scaleband_central_deviation[i_sb][i_d][i_o][i_x][i_b] << ")"
		   << "^{" << (scaleband_maximum_result[i_sb][i_d][i_o][i_x][i_b] / scaleband_central_result[i_sb][i_d][i_o][i_x][i_b] - 1.) * 100. << char(92) << "%}"
		   << "_{" << (scaleband_minimum_result[i_sb][i_d][i_o][i_x][i_b] / scaleband_central_result[i_sb][i_d][i_o][i_x][i_b] - 1.) * 100. << char(92) << "%}";
	      */	      
	      if (j_o == output_order.size() - 1){
		outf << " " << char(92) << char(92) << endl;
		//		outf << char(92) << "hline" << endl;
	      }
	      else {outf << " &" << endl;}
	      //  if (i < minresult.size() - 1){outf << setw(16) << setprecision(8) << variable[i + 1] << setw(16) << setprecision(8) << minresult[i_b] << "   " << setw(16) << setprecision(8) << mindeviation[i_b] << endl;}
	    }
	    
	    outf << " &" << endl;
    	    for (int j_o = 0; j_o < output_order.size(); j_o++){
	      int i_o = no_output_order[j_o];
	      if (i_o == -1){continue;}

	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "].size() = " << scaleband_central_deviation[i_sb].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "].size() = " << scaleband_central_deviation[i_sb][i_d].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "][" << i_o << "].size() = " << scaleband_central_deviation[i_sb][i_d][i_o].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "][" << i_o << "][" << i_x << "].size() = " << scaleband_central_deviation[i_sb][i_d][i_o][i_x].size() << endl;
	      
	      if (output_order[j_o] == output_reference){
		//		outf << "$0$";
		outf << "$" << output_latex_percent(scaleband_central_result[i_sb][i_d][no_ref2][i_x][i_b], scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], 1) << "$";
		if (j_o == output_order.size() - 1){
		  outf << " " << char(92) << char(92) << endl;
		  outf << char(92) << "hline" << endl;
		}
		else {outf << " &" << endl;}
	      }
	      else {
		//		double temp_double = scaleband_central_result[i_sb][i_d][no_ref][i_x][i_b] + scaleband_central_result[i_sb][i_d][i_o][i_x][i_b];
		outf << "$" << output_latex_percent(scaleband_central_result[i_sb][i_d][no_ref][i_x][i_b], scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], 1) << "$";
	      
		if (j_o == output_order.size() - 1){
		  outf << " " << char(92) << char(92) << endl;
		  outf << char(92) << "hline" << endl;
		}
		else {outf << " &" << endl;}
		//  if (i < minresult.size() - 1){outf << setw(16) << setprecision(8) << variable[i + 1] << setw(16) << setprecision(8) << minresult[i_b] << "   " << setw(16) << setprecision(8) << mindeviation[i_b] << endl;}
	      }
	    }
	  }
	  outf << char(92) << "end{tabular}" << endl;
	  outf << char(92) << "end{center}" << endl;
	  outf << char(92) << "caption{" << char(92) << "captiontext}" << endl;
	  outf << char(92) << "end{table}" << endl;
	  outf.close();
	}
      }
    }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_generic::output_distribution_table_IS_splitting(){
  Logger logger("summary_generic::get_summary");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_sb = 0; i_sb < outpath_scaleband.size(); i_sb++){
    // create table output here (composition of NNLO corrections according to PDF channels) !!!

    for (int x_q = 0; x_q < osi_value_qTcut_distribution.size() + 1; x_q++){
      // only copy&paste: more elegant soultion for directory_qTcut !!!
      string directory_qTcut;
      string in_directory_qTcut;
      string out_directory_qTcut;
      if (x_q == osi_value_qTcut_distribution.size()){
	directory_qTcut = "";
	in_directory_qTcut = "";
	out_directory_qTcut = "";
      }
      else {
	stringstream qTcut_ss;
	qTcut_ss << "qTcut-" << osi_value_qTcut_distribution[x_q];
	directory_qTcut = "/" + qTcut_ss.str();
	in_directory_qTcut = "";
	out_directory_qTcut = directory_qTcut;
      }
      
      logger << LOG_DEBUG << "x_q = " << x_q << "   directory_qTcut = " << directory_qTcut << endl;
      for (int i_d = 0; i_d < oset.extended_distribution.size(); i_d++){
	if (!switch_output_distribution_table[i_d]){continue;}

	vector<string> name_plot;
	string temp_sdd = "." + oset.extended_distribution[i_d].xdistribution_name;
	
	name_plot.push_back("composition.plot" + temp_sdd + ".dat");
	name_plot.push_back("composition.norm" + temp_sdd + ".dat");
	
	if (i_d >= osi_dat.size()){
	  int i_ddd = i_d - osi_dat.size();
	  
	  // recombined distributions (essentially for validation)
	  stringstream name_rec;
	  name_rec << ".rec." << osi_dddat[i_ddd].distribution_2.xdistribution_name << ".from." << osi_dddat[i_ddd].name;
	  //	  string s0;
	  name_plot.push_back("composition.plot" + name_rec.str() + ".tex");
	  name_plot.push_back("composition.norm" + name_rec.str()  + ".tex");
	  
	  for (int i_b1 = 0; i_b1 < osi_dddat[i_ddd].distribution_1.n_bins; i_b1++){
	    stringstream name_split_bin;
	    name_split_bin << "_" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1] << "-" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1 + 1];
	    stringstream name_split_lt;
	    name_split_lt << "_lt" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1];
	    stringstream name_split_ge;
	    name_split_ge << "_ge" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1];
	    
	    name_plot.push_back("composition.norm.norm.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("composition.norm.plot.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("composition.plot.norm.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    name_plot.push_back("composition.plot.plot.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".tex");
	    
	    name_plot.push_back("composition.norm.two." + osi_dddat[i_ddd].name + name_split_lt.str() + ".tex");
	    name_plot.push_back("composition.norm.two." + osi_dddat[i_ddd].name + name_split_ge.str() + ".tex");
	    name_plot.push_back("composition.plot.two." + osi_dddat[i_ddd].name + name_split_lt.str() + ".tex");
	    name_plot.push_back("composition.plot.two." + osi_dddat[i_ddd].name + name_split_ge.str() + ".tex");
	  }
	}
	


	
	vector<string> output_order(6);
	output_order[0] = "dNNLO.QCD";
	output_order[1] = "dNNLO.QCD.qqxD";
	output_order[2] = "dNNLO.QCD.qqxND";
	output_order[3] = "dNNLO.QCD.qq";
	output_order[4] = "dNNLO.QCD.gq";
	output_order[5] = "dNNLO.QCD.gg";

	vector<int> no_output_order(output_order.size(), -1);
	for (int j_o = 0; j_o < output_order.size(); j_o++){ 
	  for (int i_o = 0; i_o < yorder.size(); i_o++){
	    if (output_order[j_o] == yorder[i_o].resultdirectory){no_output_order[j_o] = i_o; break;}
	  }
	}

	for (int j_o = 0; j_o < output_order.size(); j_o++){ 
	  logger << LOG_INFO << left
		 << setw(4) << j_o
		 << setw(20) << output_order[j_o]
		 << " -> "
		 << setw(4) << no_output_order[j_o] << endl;
	}

	int n_output_version = name_plot.size();
	for (int i_x = 0; i_x < n_output_version; i_x++){
	  logger << LOG_INFO << "i_x = " << i_x << endl;

	  string outfilename_latex_table = final_resultdirectory + "/" + outpath_scaleband[i_sb] + directory_qTcut + "/" + name_plot[i_x];

	  logger << LOG_INFO << "outfilename_latex_table = " <<  outfilename_latex_table << endl;

	  ofstream outf;
	  outf.open(outfilename_latex_table.c_str(), ofstream::out | ofstream::trunc);

	  outf << char(92) << "renewcommand" << char(92) << "arraystretch{1.5}" << endl;
	  outf << char(92) << "begin{table}" << endl;
	  outf << char(92) << "begin{center}" << endl;
	  outf << char(92) << "begin{tabular}{|c|";
	  for (int j_o = 0; j_o < output_order.size(); j_o++){
	    int i_o = no_output_order[j_o];
	    if (i_o == -1){continue;}
	    outf << "c|";
	    /*
	    if (j_o == 0){outf << "c|";}
	    else {outf << "cc|";}
	    */
	  }
	  outf << "c|"; //validation !!!
	  outf << "}" << endl;
	  outf << char(92) << "hline" << endl;
	  outf << "$" << char(92) << "sqrt{s}" << char(92) << ", (" << char(92) << "mathrm{TeV})$ &" << endl;
	  for (int j_o = 0; j_o < output_order.size(); j_o++){
	    int i_o = no_output_order[j_o];
	    if (i_o == -1){continue;}

	    if (output_order[j_o] == "dNNLO.QCD"){outf << "$d" << char(92) << "sigma_{" << char(92) << "mathrm{NNLO}}" << char(92) << ", (" << char(92) << "mathrm{" << oset.unit_distribution << "})$";}
	    /*
	    else if (output_order[j_o] == "dNNLO.QCD.qq"){outf << "$d" << char(92) << "sigma_{" << char(92) << "mathrm{NNLO}_{qq}}";}
	    else if (output_order[j_o] == "dNNLO.QCD.gq"){outf << "$d" << char(92) << "sigma_{" << char(92) << "mathrm{NNLO}_{gq}}";}
	    else if (output_order[j_o] == "dNNLO.QCD.gg"){outf << "$d" << char(92) << "sigma_{" << char(92) << "mathrm{NNLO}_{gg}}";}
	    */
	    else if (output_order[j_o] == "dNNLO.QCD.qqxD"){outf << "$d" << char(92) << "sigma^{q" << char(92) << "bar{q}}_{D}/d" << char(92) << "sigma$";}
	    else if (output_order[j_o] == "dNNLO.QCD.qqxND"){outf << "$d" << char(92) << "sigma^{q" << char(92) << "bar{q}}_{ND}/d" << char(92) << "sigma$";}
	    else if (output_order[j_o] == "dNNLO.QCD.qq"){outf << "$d" << char(92) << "sigma^{qq}/d" << char(92) << "sigma$";}
	    else if (output_order[j_o] == "dNNLO.QCD.gq"){outf << "$d" << char(92) << "sigma^{gq}/d" << char(92) << "sigma$";}
	    else if (output_order[j_o] == "dNNLO.QCD.gg"){outf << "$d" << char(92) << "sigma^{gg}/d" << char(92) << "sigma$";}
	    else {outf << "$" << char(92) << "sigma_{" << char(92) << "mathrm{output_order[j_o]}}$";}
	    /*
	    outf << char(92) << ", (" << char(92) << "mathrm{" << oset.unit_distribution << "})$";
	    if (output_order[j_o] == "dNNLO.QCD.qq" || output_order[j_o] == "dNNLO.QCD.gq" || output_order[j_o] == "dNNLO.QCD.gg"){outf << " & $(" << char(92) << "%)$";}
	    */
	    
	    if (j_o == output_order.size() - 1){outf << " & " << char(92) << char(92) << endl;} // validation !!!
	    //	    if (j_o == output_order.size() - 1){outf << " " << char(92) << char(92) << endl;}
	    else {outf << " &" << endl;}
	  }
	  outf << char(92) << "hline" << endl;

	  // better: variable... -> doubled lines for plotting reasons etc.
	  /*
	  for (int i_b = 0; i_b < oset.extended_distribution[i_d].n_bins; i_b++){
	    outf << oset.extended_distribution[i_d].bin_edge[i_b] << " &" << endl;
	  */
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "].size() = " << scaleband_variable[i_sb].size() << endl;
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "][" << i_d << "].size() = " << scaleband_variable[i_sb][i_d].size() << endl;
	  logger << LOG_DEBUG_VERBOSE << "scaleband_variable[" << i_sb << "][" << i_d << "][" << i_x << "].size() = " << scaleband_variable[i_sb][i_d][i_x].size() << endl;

	  for (int i_b = 0; i_b < scaleband_variable[i_sb][i_d][i_x].size() - 1; i_b++){
	    //	  for (int i_b = 0; i_b < scaleband_variable[i_sb][i_d][i_x].size(); i_b++){
	    outf << oset.E_CMS / 1000 << " " << char(92) << "nameobservable{" << scaleband_variable[i_sb][i_d][i_x][i_b] << "} &" << endl;

	    double checksum = 0.;
	    
	    for (int j_o = 0; j_o < output_order.size(); j_o++){
	      int i_o = no_output_order[j_o];
	      if (i_o == -1){continue;}

	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "].size() = " << scaleband_central_deviation[i_sb].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "].size() = " << scaleband_central_deviation[i_sb][i_d].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "][" << i_o << "].size() = " << scaleband_central_deviation[i_sb][i_d][i_o].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "scaleband_central_deviation[" << i_sb << "][" << i_d << "][" << i_o << "][" << i_x << "].size() = " << scaleband_central_deviation[i_sb][i_d][i_o][i_x].size() << endl;

	      if (j_o == 0){
		outf << "$" << output_result_deviation(scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], scaleband_central_deviation[i_sb][i_d][i_o][i_x][i_b], 1) << "$";
	      }
	      else if (j_o > 0){
		/*
		outf << "$" << output_result_deviation(scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], scaleband_central_deviation[i_sb][i_d][i_o][i_x][i_b], 1) << "$";
		*/
		checksum += scaleband_central_result[i_sb][i_d][i_o][i_x][i_b];
		double temp_double = scaleband_central_result[i_sb][i_d][no_output_order[0]][i_x][i_b] + scaleband_central_result[i_sb][i_d][i_o][i_x][i_b];
		/*
		outf << " & ";
		*/
		outf << "$" << output_latex_percent(scaleband_central_result[i_sb][i_d][no_output_order[0]][i_x][i_b], temp_double, 1) << "$";
	      }
	      /*
	      outf << "$" << output_result_deviation(scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], scaleband_central_deviation[i_sb][i_d][i_o][i_x][i_b], 1)
		   << "^{" << output_latex_percent(scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], scaleband_maximum_result[i_sb][i_d][i_o][i_x][i_b], 1) << "}"
		   << "_{" << output_latex_percent(scaleband_central_result[i_sb][i_d][i_o][i_x][i_b], scaleband_minimum_result[i_sb][i_d][i_o][i_x][i_b], 1) << "}$";
	      */
	      
	      if (j_o == output_order.size() - 1){
		
		outf << " & " << output_latex_percent(scaleband_central_result[i_sb][i_d][no_output_order[0]][i_x][i_b], checksum, 1) << " " << char(92) << char(92) << endl; // validation !!!
		//		outf << " " << char(92) << char(92) << "   % checksum: " << output_latex_percent(scaleband_central_result[i_sb][i_d][no_output_order[0]][i_x][i_b], checksum, 1) << endl;
		outf << char(92) << "hline" << endl;
	      }
	      else {outf << " &" << endl;}
	      //  if (i < minresult.size() - 1){outf << setw(16) << setprecision(8) << variable[i + 1] << setw(16) << setprecision(8) << minresult[i_b] << "   " << setw(16) << setprecision(8) << mindeviation[i_b] << endl;}
	    }
	  }
	  outf << char(92) << "end{tabular}" << endl;
	  outf << char(92) << "end{center}" << endl;
	  outf << char(92) << "caption{" << char(92) << "captiontext}" << endl;
	  outf << char(92) << "end{table}" << endl;
	  outf.close();
	}
      }
    }

  }
  
  logger << LOG_DEBUG << "finished" << endl;
}

