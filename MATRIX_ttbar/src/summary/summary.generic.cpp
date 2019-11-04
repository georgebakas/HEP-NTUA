#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"
////////////////////
//  constructors  //
////////////////////
summary_generic::summary_generic(){}

summary_generic::summary_generic(munich & xmunich){
  Logger logger("summary_generic::summary_generic");
  logger << LOG_DEBUG << "called" << endl;

  logger << LOG_DEBUG << "xmunich: order = " << xmunich.order << endl;
  logger << LOG_DEBUG << "xmunich: infilename = " << xmunich.infilename << endl;
  logger << LOG_DEBUG << "xmunich: infilename_scaleband = " << xmunich.infilename_scaleband << endl;

  csi = xmunich.csi;

  generic = &(xmunich.generic);
  isi = &(xmunich.isi);

  oset = observable_set(xmunich.isi, xmunich.csi);

  order = xmunich.order;
  infilename = xmunich.infilename;
  infilename_scaleband = xmunich.infilename_scaleband;

  logger << LOG_DEBUG << "order = " << order << endl;
  logger << LOG_DEBUG << "infilename = " << infilename << endl;
  logger << LOG_DEBUG << "infilename_scaleband = " << infilename_scaleband << endl;

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_generic::get_summary(){
  Logger logger("summary_generic::get_summary");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // default values:
  run_min_n_step = 10;
  run_factor_max = 2;
  run_min_n_event = 50000;
  run_max_time_per_job = 2 * 60 * 60;
  run_min_number_of_jobs_for_one_channel = 3;
  average_factor = 0;
  switch_generate_rundirectories = 0;
  no_qTcut_runtime_estimate = 0;
  deviation_tolerance_factor = 5.;

  switch_extrapolation_result = 2;
  switch_extrapolation_distribution = 2;
  
  min_qTcut_extrapolation = 0.05;
  max_qTcut_extrapolation = 0.5;
  error_extrapolation_range_chi2 = 4.;

  min_max_value_extrapolation_range = 0.25;
  min_n_qTcut_extrapolation_range = 10;

  switch_output_plot = 1;
  switch_output_result = 2;
  switch_output_overview = 3;
  // 0 - none
  // 1 - order
  // 2 - order + list
  // 3 - order + list + contribution
  // 4 - order + list + contribution + subprocess


  switch_output_subprocess = 0;
  // 0 - none
  // 1 - overview
  // 2 - overview + plot
  // 3 - overview + plot + result
  switch_output_contribution = 0;
  switch_output_list = 0;
  switch_output_order = 1;


  switch_output_table_order = 0;
  switch_output_table_Kfactor = 0;
  switch_output_table_crosssection_Kfactor = 0;
  switch_output_table_IS_splitting = 0;

  
  infix_name_moment.resize(oset.n_moments + 1);
  int start_i_m = 0; // temporary !!!
  for (int i_m = start_i_m; i_m < oset.n_moments + 1; i_m++){
    stringstream temp_ss;
    if (i_m < oset.n_moments){temp_ss << ".moment_" << i_m;}
    infix_name_moment[i_m] = temp_ss.str();
  }
  x_m = 0;

  // readin from specified input file:
  readin_combination_infile();




  system_execute(logger, "mkdir " + final_resultdirectory);
  
  /*
  if (osi_switch_CV){
    system_execute(logger, "mkdir " + final_resultdirectory + "/CV");
  }
  if (osi_switch_TSV){
    for (int i_s = 0; i_s < osi_n_extended_set_TSV; i_s++){
      system_execute(logger, "mkdir " + final_resultdirectory + "/" + osi_name_extended_set_TSV[i_s]);
    }
  }
  */



  for (int i_l = 0; i_l < list_contribution_file.size(); i_l++){logger << LOG_INFO << "list_contribution_file[" << i_l << "] = " << list_contribution_file[i_l] << endl;}

  resultdirectory.resize(list_contribution_file.size());
  oset.initialization_result();
  oset.initialization_CV();

  name_scale_variation_TSV.resize(7);
  name_scale_variation_TSV[0] = "complete";
  name_scale_variation_TSV[1] = "equal";
  name_scale_variation_TSV[2] = "ren";
  name_scale_variation_TSV[3] = "fact";
  name_scale_variation_TSV[4] = "antipodal";
  name_scale_variation_TSV[5] = "7-point";
  name_scale_variation_TSV[6] = "9-point";

  if (osi_switch_TSV){initialization_scaleset_TSV();}

  scalename_TSV.resize(osi_n_extended_set_TSV);
  for (int i_s = 0; i_s < osi_n_extended_set_TSV; i_s++){
    scalename_TSV[i_s].resize(osi_n_scale_ren_TSV[i_s], vector<string> (osi_n_scale_fact_TSV[i_s]));
    if (osi_switch_distribution_TSV[i_s] != 0){
      system_execute(logger, "mkdir " + final_resultdirectory + "/" + osi_name_extended_set_TSV[i_s]);
      for (int i_r = 0; i_r < osi_n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < osi_n_scale_fact_TSV[i_s]; i_f++){
	  stringstream scalepart;
	  scalepart << "scale." << i_r << "." << i_f;
	  scalename_TSV[i_s][i_r][i_f] = final_resultdirectory + "/" + osi_name_extended_set_TSV[i_s] + "/" + scalepart.str();
	  system_execute(logger, "mkdir " + scalename_TSV[i_s][i_r][i_f]);
	  logger << LOG_DEBUG_VERBOSE << "mkdir " << scalename_TSV[i_s][i_r][i_f] << endl;
	}
      }
    }
  }


  if (osi_switch_CV){
    system_execute(logger, "mkdir " + final_resultdirectory + "/CV");
  }

  if      (osi_switch_CV == 1){name_variation_CV = "equal";}
  else if (osi_switch_CV == 2){name_variation_CV = "ren";}
  else if (osi_switch_CV == 3){name_variation_CV = "fac";}
  else if (osi_switch_CV == 4){name_variation_CV = "antipodal";}
  else if (osi_switch_CV == 5){name_variation_CV = "7-point";}
  else if (osi_switch_CV == 6){name_variation_CV = "9-point";}

  scalename_CV.resize(osi_n_scales_CV);
  for (int i_s = 0; i_s < osi_n_scales_CV; i_s++){
    stringstream temp;
    temp << "scale." << i_s;// << "/";
    scalename_CV[i_s] = temp.str();
  }


  subgroup.push_back(csi.process_class); // ???

  // Initialization of xlist from 'list_contribution_file':
  initialization_summary_list();

  logger << LOG_INFO << "BEFORE   &xlist = " << &xlist << endl;
  for (int i_l = 0; i_l < xlist.size(); i_l++){
    logger << LOG_INFO << "BEFORE   &xlist[" << i_l << "] = " << &xlist[i_l] << endl;
    for (int i_c = 1; i_c < xlist[i_l].xcontribution.size(); i_c++){
      xlist[i_l].xcontribution[i_c].ylist = &xlist[i_l];
      for (int i_p = 1; i_p < xlist[i_l].xcontribution[i_c].xsubprocess.size(); i_p++){
	xlist[i_l].xcontribution[i_c].xsubprocess[i_p].ycontribution = &xlist[i_l].xcontribution[i_c];
      }
    }
  }
  for (int i_l = 0; i_l < xlist.size(); i_l++){
    logger << LOG_INFO << "AFTER   &xlist[" << i_l << "] = " << &xlist[i_l] << endl;
  }
  for (int i_l = 0; i_l < xlist.size(); i_l++){
    for (int i_c = 1; i_c < xlist[i_l].xcontribution.size(); i_c++){
      logger << LOG_INFO << "AFTER   &xlist[" << i_l << "].xcontribution[" << i_c << "] = " << &xlist[i_l].xcontribution[i_c] << endl;
    }
  }

  for (int i_l = 0; i_l < xlist.size(); i_l++){
    logger << LOG_INFO << "xlist[" << i_l << "]" << endl; 
    logger << LOG_INFO << "list      active_qTcut = " << xlist[i_l].active_qTcut << "   output_n_qTcut = " << setw(3) << xlist[i_l].output_n_qTcut << "   selection_n_qTcut = " << xlist[i_l].selection_n_qTcut << endl;
    for (int i_c = 0; i_c < xlist[i_l].xcontribution.size(); i_c++){
      logger << LOG_INFO << "i_c = " << i_c << "   active_qTcut = " << xlist[i_l].xcontribution[i_c].active_qTcut << "   output_n_qTcut = " << setw(3) << xlist[i_l].xcontribution[i_c].output_n_qTcut << "   selection_n_qTcut = " << xlist[i_l].xcontribution[i_c].selection_n_qTcut << endl;
    }
  }



  /*
  system_execute(logger, "mkdir " + final_resultdirectory);
  if (osi_switch_CV){system_execute(logger, "mkdir " + final_resultdirectory + "/CV");}
  if (osi_switch_TSV){for (int i_s = 0; i_s < osi_n_extended_set_TSV; i_s++){system_execute(logger, "mkdir " + final_resultdirectory + "/" + osi_name_extended_set_TSV[i_s]);}}
  */  

  /*
  for (int i_l = 0; i_l < list_contribution_file.size(); i_l++){
    for (int i_c = 1; i_c < xlist[i_l].xcontribution.size(); i_c++){
      logger << LOG_INFO << "BEFORE xlist[" << i_l << "].xcontribution[" << i_c << "].ylist->resultdirectory = " << xlist[i_l].xcontribution[i_c].ylist->resultdirectory << endl;
    }
  }
  */

  // order of linked libraries !!!
  double test_d = 1.;
  string test_s = time_hms_from_double(test_d);
  
  for (int i_o = 0; i_o < yorder.size(); i_o++){
    yorder[i_o].active_qTcut = 0;
    yorder[i_o].output_n_qTcut = 1;
    yorder[i_o].selection_n_qTcut = 1;
    for (int i_c = 0; i_c < yorder[i_o].combination.size(); i_c++){
      for (int i_l = 0; i_l < yorder[i_o].contribution_file.size(); i_l++){
	int x_l = mapping_contribution_file[yorder[i_o].contribution_file[i_l]];
	if (xlist[x_l].active_qTcut){
	  yorder[i_o].active_qTcut = xlist[x_l].active_qTcut;
	  yorder[i_o].output_n_qTcut = xlist[x_l].output_n_qTcut;
	  yorder[i_o].selection_n_qTcut = xlist[x_l].selection_n_qTcut;
	}
      }
    }
  }
  

  for (int i_l = 0; i_l < xlist.size(); i_l++){
    logger << LOG_INFO << "xlist[" << i_l << "] after order:" << endl; 
    logger << LOG_INFO << "list      active_qTcut = " << xlist[i_l].active_qTcut << "   output_n_qTcut = " << setw(3) << xlist[i_l].output_n_qTcut << "   selection_n_qTcut = " << xlist[i_l].selection_n_qTcut << endl;
    for (int i_c = 0; i_c < xlist[i_l].xcontribution.size(); i_c++){
      logger << LOG_INFO << "i_c = " << i_c << "   active_qTcut = " << xlist[i_l].xcontribution[i_c].active_qTcut << "   output_n_qTcut = " << setw(3) << xlist[i_l].xcontribution[i_c].output_n_qTcut << "   selection_n_qTcut = " << xlist[i_l].xcontribution[i_c].selection_n_qTcut << endl;
    }
  }


  for (int i_o = 0; i_o < yorder.size(); i_o++){
    logger << LOG_INFO << "yorder[" << i_o << "].resultdirectory = " << setw(16) << yorder[i_o].resultdirectory << "   active_qTcut = " << yorder[i_o].active_qTcut << "   output_n_qTcut = " << setw(3) << yorder[i_o].output_n_qTcut << "   selection_n_qTcut = " << yorder[i_o].selection_n_qTcut << endl;
  }


  /*
  for (int i_l = 0; i_l < list_contribution_file.size(); i_l++){
    for (int i_c = 1; i_c < xlist[i_l].xcontribution.size(); i_c++){
      logger << LOG_INFO << "BEFORE xlist[" << i_l << "].xcontribution[" << i_c << "].ylist->resultdirectory = " << xlist[i_l].xcontribution[i_c].ylist->resultdirectory << endl;
    }
  }
  */

  if (order == "result"){
    for (int i_l = 0; i_l < list_contribution_file.size(); i_l++){
      logger << LOG_DEBUG_VERBOSE << "remove_run:   list_contribution_file[" << i_l << "] = " << list_contribution_file[i_l] << endl;
      for (int i_c = 1; i_c < xlist[i_l].xcontribution.size(); i_c++){
	//	logger << LOG_INFO << "xlist[" << i_l << "].xcontribution[" << i_c << "].ylist->resultdirectory = " << xlist[i_l].xcontribution[i_c].ylist->resultdirectory << endl;
	xlist[i_l].xcontribution[i_c].readin_contribution_remove_run();
	//	logger << LOG_INFO << "xlist[" << i_l << "].xcontribution[" << i_c << "].ylist->resultdirectory = " << xlist[i_l].xcontribution[i_c].ylist->resultdirectory << endl;
      }
    }
    
    for (int i_l = 0; i_l < list_contribution_file.size(); i_l++){
      logger << LOG_DEBUG_VERBOSE << "collect_contribution:   list_contribution_file[" << i_l << "] = " << list_contribution_file[i_l] << endl;

      /*
      for (int i_c = 1; i_c < xlist[i_l].xcontribution.size(); i_c++){
	logger << LOG_INFO << "xlist[" << i_l << "].xcontribution[" << i_c << "].ylist->resultdirectory = " << xlist[i_l].xcontribution[i_c].ylist->resultdirectory << endl;
      }
      */

      if (osi_switch_CV){xlist[i_l].collect_contribution_result_CV();}
      if (osi_switch_TSV){xlist[i_l].collect_contribution_result_TSV();}
    }

    for (int i_l = 0; i_l < list_contribution_file.size(); i_l++){
      logger << LOG_DEBUG_VERBOSE << "readin_runtime:   list_contribution_file[" << i_l << "] = " << list_contribution_file[i_l] << endl;
      for (int i_c = 1; i_c < xlist[i_l].xcontribution.size(); i_c++){
	xlist[i_l].xcontribution[i_c].readin_runtime_contribution();

	for (int i_p = 1; i_p < xlist[i_l].xcontribution[i_c].xsubprocess.size(); i_p++){
	  if (osi_switch_CV){
	    vector<vector<vector<double> > > ().swap(xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_run_CV);
	    vector<vector<vector<double> > > ().swap(xlist[i_l].xcontribution[i_c].xsubprocess[i_p].deviation_run_CV);
	  }
	  if (osi_switch_TSV){
	    vector<vector<vector<vector<vector<vector<double> > > > > > ().swap(xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_run_TSV);
	    vector<vector<vector<vector<vector<vector<double> > > > > > ().swap(xlist[i_l].xcontribution[i_c].xsubprocess[i_p].deviation_run_TSV);
	  }
	}
      }
      xlist[i_l].calculate_runtime();
    }

    for (int i_o = 0; i_o < yorder.size(); i_o++){
      for (int i_l = 0; i_l < yorder[i_o].contribution_file.size(); i_l++){
	int x_l = mapping_contribution_file[yorder[i_o].contribution_file[i_l]];
	yorder[i_o].xlist.push_back(&xlist[x_l]); 
	//	yorder[i_o].xlist.push_back(xlist[x_l]); 
      }
    }
    
    if (osi_switch_CV){collect_order_result_CV();} // partially contains runtime evaluation !!!

    // old position:    determine_runtime_result_new();

    if (osi_switch_TSV){collect_order_result_TSV();}

    determine_runtime_result_new();

    // add scaleband output for cross sections !!!

    determine_runtime();
  }


  
  if (order == "distribution"){
    initialization_distribution();

    if (osi_switch_CV){
      for (int i_l = 0; i_l < list_contribution_file.size(); i_l++){
	// internal function of summary.generic !!!
	xlist[i_l].collect_contribution_distribution_CV();

      }
      // internal function of summary.generic !!!
      //      collect_order_distribution_CV(final_resultdirectory, mapping_contribution_file, subgroup, xlist, oset.extended_distribution, oset.fakeasymfactor, oset, yorder);
      //      collect_order_distribution_CV();
    }

    if (osi_switch_TSV){
      for (int i_l = 0; i_l < list_contribution_file.size(); i_l++){
	// internal function of summary.generic !!!
	xlist[i_l].collect_contribution_distribution_TSV();
      }
      // internal function of summary.generic !!!
      //      collect_order_distribution_TSV(final_resultdirectory, mapping_contribution_file, subgroup, xlist, oset.extended_distribution, oset.fakeasymfactor, oset, yorder);
      //      collect_order_distribution_TSV();
    }


    for (int i_o = 0; i_o < yorder.size(); i_o++){
      logger << LOG_DEBUG << "yorder[i_o = " << i_o << "].contribution_file.size() = " << yorder[i_o].contribution_file.size() << endl;
      for (int i_l = 0; i_l < yorder[i_o].contribution_file.size(); i_l++){
	int x_l = mapping_contribution_file[yorder[i_o].contribution_file[i_l]];
	logger << LOG_DEBUG << "yorder[i_o = " << i_o << "].xlist.push_back(xlist[x_l = " << x_l << "]);" << endl;
	yorder[i_o].xlist.push_back(&xlist[x_l]); 
	//	yorder[i_o].xlist.push_back(xlist[x_l]); 
      }
    }
    if (osi_switch_CV){collect_order_distribution_CV();}
    if (osi_switch_TSV){collect_order_distribution_TSV();}


    logger << LOG_INFO << "infilename_scaleband = " << infilename_scaleband << endl;
    logger << LOG_INFO << "clear xlist since not needed any longer !!!" << endl;
    for (int i_l = 0; i_l < xlist.size(); i_l++){
      xlist[i_l] = summary_list();
    }
    
    xlist.clear();

    logger << LOG_INFO << "clear xlist done." << endl;
    
    if (infilename_scaleband != ""){
      logger << LOG_DEBUG << "final_resultdirectory = " << final_resultdirectory << endl;
      // internal function of summary.generic !!!
      logger << LOG_INFO << "Determination of scaleband is started." << endl;
      determine_scaleband(); //final_resultdirectory, infilename_scaleband, oset, osi_unit_factor_distribution, oset.extended_distribution, oset.fakeasymfactor, yorder);
    }

    // add runtime estimate based on a selected distribution / a number of selected distributions
  }
  if (order == "scaleband"){
    initialization_distribution();

    if (infilename_scaleband != ""){
      logger << LOG_DEBUG << "final_resultdirectory = " << final_resultdirectory << endl;
      // internal function of summary.generic !!!
      logger << LOG_INFO << "Determination of scaleband is started." << endl;
      determine_scaleband(); //final_resultdirectory, infilename_scaleband, oset, osi_unit_factor_distribution, oset.extended_distribution, oset.fakeasymfactor, yorder);
    }
  }
  
}

void summary_generic::initialization_distribution(){
  Logger logger("summary_generic::initialization_distribution");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //    vector<xdistribution> extended_distribution; // !!! -> oset.extended_distribution
  //    vector<double> fakeasymfactor; // !!! -> oset.fakeasymfactor
  oset.determine_extended_distribution();
  logger << LOG_INFO << "output_selection_distribution.size() = " << output_selection_distribution.size() << endl;
  for (int i_s = 0; i_s < output_selection_distribution.size(); i_s++){
    logger << LOG_INFO << "output_selection_distribution[" << i_s << "] = " << output_selection_distribution[i_s] << endl;
    
  }
  logger << LOG_INFO << "switch_output_distribution.size() = " << switch_output_distribution.size() << endl;
  
  if (output_selection_distribution.size() == 0){
    switch_output_distribution.resize(oset.extended_distribution.size(), 1);
  }
  else {
    switch_output_distribution.resize(oset.extended_distribution.size(), 0);
    for (int i_d = 0; i_d < oset.extended_distribution.size(); i_d++){
      logger << LOG_INFO << "oset.extended_distribution[" << i_d << "].xdistribution_name = " << oset.extended_distribution[i_d].xdistribution_name << endl;
      //	if (!switch_output_distribution[i_d]){continue;}
      for (int i_s = 0; i_s < output_selection_distribution.size(); i_s++){
	if (oset.extended_distribution[i_d].xdistribution_name == output_selection_distribution[i_s]){switch_output_distribution[i_d] = 1; break;}
      }
      if (switch_output_distribution[i_d] == 1){
	logger << LOG_INFO << "switch_output_distribution[" << setw(2) << i_d << "] = " << switch_output_distribution[i_d] << endl;
      }
    }
  }
  
  logger << LOG_DEBUG << "oset.extended_distribution.size() = " << oset.extended_distribution.size() << endl;
  logger << LOG_DEBUG << "output_selection_distribution.size() = " << output_selection_distribution.size() << endl;
  for (int i_s = 0; i_s < output_selection_distribution.size(); i_s++){
    logger << LOG_DEBUG << "output_selection_distribution[" << i_s << "] = " << output_selection_distribution[i_s] << endl;
    
  }
  logger << LOG_DEBUG << "switch_output_distribution.size() = " << switch_output_distribution.size() << endl;
  for (int i_s = 0; i_s < switch_output_distribution.size(); i_s++){
    logger << LOG_DEBUG << "switch_output_distribution[" << i_s << "] = " << switch_output_distribution[i_s] << endl;
    
  }




  logger << LOG_INFO << "output_selection_distribution_table.size() = " << output_selection_distribution_table.size() << endl;
  for (int i_s = 0; i_s < output_selection_distribution_table.size(); i_s++){
    logger << LOG_INFO << "output_selection_distribution_table[" << i_s << "] = " << output_selection_distribution_table[i_s] << endl;
    
  }
  logger << LOG_INFO << "switch_output_distribution_table.size() = " << switch_output_distribution_table.size() << endl;
  
  if (output_selection_distribution_table.size() == 0){
    switch_output_distribution_table.resize(oset.extended_distribution.size(), 1);
  }
  else {
    switch_output_distribution_table.resize(oset.extended_distribution.size(), 0);
    for (int i_d = 0; i_d < oset.extended_distribution.size(); i_d++){
      logger << LOG_INFO << "oset.extended_distribution[" << i_d << "].xdistribution_name = " << oset.extended_distribution[i_d].xdistribution_name << endl;
      //	if (!switch_output_distribution_table[i_d]){continue;}
      for (int i_s = 0; i_s < output_selection_distribution_table.size(); i_s++){
	if (oset.extended_distribution[i_d].xdistribution_name == output_selection_distribution_table[i_s]){switch_output_distribution_table[i_d] = 1; break;}
      }
      if (switch_output_distribution_table[i_d] == 1){
	logger << LOG_INFO << "switch_output_distribution_table[" << setw(2) << i_d << "] = " << switch_output_distribution_table[i_d] << endl;
      }
    }
  }
  
  logger << LOG_DEBUG << "oset.extended_distribution.size() = " << oset.extended_distribution.size() << endl;
  logger << LOG_DEBUG << "output_selection_distribution_table.size() = " << output_selection_distribution_table.size() << endl;
  for (int i_s = 0; i_s < output_selection_distribution_table.size(); i_s++){
    logger << LOG_DEBUG << "output_selection_distribution_table[" << i_s << "] = " << output_selection_distribution_table[i_s] << endl;
    
  }









  

  // only used for asymmetries: could most likely be removed later...
  bin_edge.resize(oset.extended_distribution.size());
  for (int i_d = 0; i_d < oset.extended_distribution.size(); i_d++){
    if (!switch_output_distribution[i_d]){continue;}
    bin_edge[i_d] = oset.extended_distribution[i_d].bin_edge;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void summary_generic::initialization_scaleset_TSV(){
  Logger logger("summary_generic::initialization_scaleset_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_INFO << "output_selection_scaleset.size() = " << output_selection_scaleset.size() << endl;
  for (int j_s = 0; j_s < output_selection_scaleset.size(); j_s++){
    logger << LOG_INFO << "output_selection_scaleset[" << j_s << "] = " << output_selection_scaleset[j_s] << endl;
  }
  logger << LOG_INFO << "switch_output_scaleset.size() = " << switch_output_scaleset.size() << endl;
  
  if (output_selection_scaleset.size() == 0){
    switch_output_scaleset.resize(oset.n_extended_set_TSV, 1);
  }
  else {
    switch_output_scaleset.resize(oset.n_extended_set_TSV, 0);
    for (int i_s = 0; i_s < oset.n_extended_set_TSV; i_s++){
      logger << LOG_INFO << "oset.name_extended_set_TSV[" << i_s << "] = " << oset.name_extended_set_TSV[i_s] << endl;
      //	if (!switch_output_scaleset[i_s]){continue;}
      for (int j_s = 0; j_s < output_selection_scaleset.size(); j_s++){
	if (oset.name_extended_set_TSV[i_s] == output_selection_scaleset[j_s]){switch_output_scaleset[i_s] = 1; break;}
      }
      if (switch_output_scaleset[i_s] == 1){
	logger << LOG_INFO << "switch_output_scaleset[" << setw(2) << i_s << "] = " << switch_output_scaleset[i_s] << endl;
      }
    }
  }
  
  if (switch_output_scaleset[oset.no_reference_TSV] == 0){
    switch_output_scaleset[oset.no_reference_TSV] = 1;
    logger << LOG_INFO << "reference scale must be included in output: switch_output_scaleset[" << setw(2) << oset.no_reference_TSV << "] = " << switch_output_scaleset[oset.no_reference_TSV] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "oset.n_extended_set_TSV = " << oset.n_extended_set_TSV << endl;

  /*
  // only used for asymmetries: could most likely be removed later...
  bin_edge.resize(oset.n_extended_set_TSV);
  for (int i_d = 0; i_d < oset.n_extended_set_TSV; i_d++){
    if (!switch_output_scaleset[i_d]){continue;}
    bin_edge[i_d] = oset.extended_scaleset[i_d].bin_edge;
  }
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void summary_generic::initialization_summary_list(){
  Logger logger("summary_generic::initialization_summary_list");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  xlist.resize(list_contribution_file.size());

  for (int i_l = 0; i_l < list_contribution_file.size(); i_l++){
    logger << LOG_DEBUG << "list_contribution_file[" << i_l << "] = " << list_contribution_file[i_l] << endl;
    xlist[i_l] = summary_list(list_contribution_file[i_l], *this);

    /*
    summary_list temp_sl(list_contribution_file[i_l], *this);
    cout << "&temp_sl = " << &temp_sl << endl;
    xlist.push_back(temp_sl);
    */

    //  resultdirectory[i_l] -> e.g. LO.03, NLO.CS.13, etc.
    //  infix_order_contribution[i_l] -> e.g. born, born.i1, RA.QCD, etc.
    //  resultdirectory[i_l] = final_resultdirectory + "/" + resultdirectory[i_l]; 
    // ???
  }

  for (int i_l = 0; i_l < xlist.size(); i_l++){
    logger << LOG_DEBUG << "xlist[" << i_l << "] content" << endl;
    xlist[i_l].output_info();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void summary_generic::readin_combination_infile(){
  Logger logger("summary_generic::readin_combination_infile");
  logger << LOG_DEBUG << "called" << endl;
  vector<string> readin;
  char LineBuffer[128];
  logger << LOG_INFO << "infilename = " << infilename << endl;

  ifstream in_file(infilename.c_str());
  readin.clear();
  while (in_file.getline(LineBuffer, 1024)){readin.push_back(LineBuffer);}

  vector<string> user_variable;
  vector<string> user_value;
  int user_counter = -1;
  for (int i = 0; i < readin.size(); i++){
    if (readin[i].size() != 0){
      if (readin[i][0] != '/' && readin[i][0] != '#' && readin[i][0] != '%'){
	int start = 0;
	user_counter++;
	user_variable.push_back("");
	user_value.push_back("");
	for (int j = 0; j < readin[i].size(); j++){
	  if (start == 0 || start == 1){
	    if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 0){}
	    else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	      user_variable[user_counter].push_back(readin[i][j]);
	      if (start != 1){start = 1;}
	    }
	    else {start++;}
	  }
	  else if (start == 2){
	    if (readin[i][j] == '='){start++;}
	    else if ((readin[i][j] == ' ') || (readin[i][j] == char(9))){}
	    else {cout << readin[i][j] << "    Incorrect input in line " << i + 1 << endl; exit(1);}
	  }
	  else if (start == 3 || start == 4){
	    if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 3){}
	    else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	      user_value[user_counter].push_back(readin[i][j]);
	      if (start != 4){start = 4;}
	    }
	    else {start++;}
	  }
	  else {break;}
	}
	if (start == 0){
	  user_counter--;
	  user_variable.erase(user_variable.end() - 1, user_variable.end());
	  user_value.erase(user_value.end() - 1, user_value.end());
	}
      }
    }
  }

  int counter_directory = -1;
  int counter_combination = -1;


  phasespace_optimization.push_back("");


  vector<string> accuracy_normalization;

  //  cout << "order_combination.size() = " << order_combination.size() << endl;
  for (int i = 0; i < user_variable.size(); i++){
    logger << LOG_DEBUG << "i = " << setw(3) << i << "   " << user_variable[i] << endl;
    if (user_variable[i] == "final_resultdirectory"){
      final_resultdirectory = user_value[i];
    }
    else if (user_variable[i] == "switch_generate_rundirectories"){
      switch_generate_rundirectories = atoi(user_value[i].c_str());
    }
    else if (user_variable[i] == "run_min_n_step"){
      run_min_n_step = atoi(user_value[i].c_str());
    }
    else if (user_variable[i] == "run_factor_max"){
      run_factor_max = atoi(user_value[i].c_str());
    }
    else if (user_variable[i] == "run_min_n_event"){
      run_min_n_event = atoll(user_value[i].c_str());
    }
    else if (user_variable[i] == "run_max_time_per_job"){
      run_max_time_per_job = 3600 * atof(user_value[i].c_str());
    }
    else if (user_variable[i] == "run_min_number_of_jobs_for_one_channel"){
      run_min_number_of_jobs_for_one_channel = atoi(user_value[i].c_str());
    }
    else if (user_variable[i] == "no_qTcut_runtime_estimate"){
      no_qTcut_runtime_estimate = atoi(user_value[i].c_str());
    }
    else if (user_variable[i] == "directory_runtime_estimate"){
      directory_runtime_estimate = user_value[i].c_str();
    }
    else if (user_variable[i] == "deviation_tolerance_factor"){
      deviation_tolerance_factor = atof(user_value[i].c_str());
    }
    else if (user_variable[i] == "min_qTcut_extrapolation"){
      min_qTcut_extrapolation = atof(user_value[i].c_str());
    }
    else if (user_variable[i] == "switch_extrapolation_result"){
      switch_extrapolation_result = atof(user_value[i].c_str());
    }
    else if (user_variable[i] == "switch_extrapolation_distribution"){
      switch_extrapolation_distribution = atof(user_value[i].c_str());
    }
    else if (user_variable[i] == "max_qTcut_extrapolation"){
      max_qTcut_extrapolation = atof(user_value[i].c_str());
    }
    else if (user_variable[i] == "min_max_value_extrapolation_range"){
      min_max_value_extrapolation_range = atof(user_value[i].c_str());
    }
    else if (user_variable[i] == "min_n_qTcut_extrapolation_range"){
      min_n_qTcut_extrapolation_range = atoi(user_value[i].c_str());
    }
    else if (user_variable[i] == "error_extrapolation_range_chi2"){
      error_extrapolation_range_chi2 = atof(user_value[i].c_str());
    }
 
    else if (user_variable[i] == "average_factor"){
      average_factor = atoi(user_value[i].c_str());
      // average_factor decides, how the results are combined:
      // average_factor = 1: conservative
      // average_factor = 0: alternative
      // average_factor = n: hybrid
    }

    else if (user_variable[i] == "phasespace_optimization"){
      phasespace_optimization.push_back("." + user_value[i]);
    }

    else if (user_variable[i] == "resultdirectory"){
      yorder.push_back(summary_order(user_value[i], *this));
      counter_directory++;
      counter_combination = 0;
    }

    else if (user_variable[i] == "contribution_file"){
      yorder[counter_directory].contribution_file.push_back(user_value[i]);
      yorder[counter_directory].combination[counter_combination].push_back(yorder[counter_directory].contribution_file.size() - 1);
    }

    else if (user_variable[i] == "combination"){
      // to be made more sophisticated...
      yorder[counter_directory].combination.push_back(vector<int> ());
      counter_combination++;
      if (user_value[i] == "+"){yorder[counter_directory].combination_type.push_back(0);}
      else if (user_value[i] == "x"){yorder[counter_directory].combination_type.push_back(1);}
    }

    else if (user_variable[i] == "accuracy_relative"){
      yorder[counter_directory].accuracy_relative = atof(user_value[i].c_str());
    }

    else if (user_variable[i] == "output_selection_distribution"){
      output_selection_distribution.push_back(user_value[i]);
    }

    else if (user_variable[i] == "output_selection_scaleset"){
      output_selection_scaleset.push_back(user_value[i]);
    }

    else if (user_variable[i] == "switch_output_subprocess"){
      switch_output_subprocess = atoi(user_value[i].c_str());
    }

    else if (user_variable[i] == "switch_output_contribution"){
      switch_output_contribution = atoi(user_value[i].c_str());
    }

    else if (user_variable[i] == "switch_output_list"){
      switch_output_list = atoi(user_value[i].c_str());
    }

    else if (user_variable[i] == "switch_output_order"){
      switch_output_order = atoi(user_value[i].c_str());
    }

    else if (user_variable[i] == "switch_output_plot"){
      switch_output_plot = atoi(user_value[i].c_str());
    }

    else if (user_variable[i] == "switch_output_result"){
      switch_output_result = atoi(user_value[i].c_str());
    }

    else if (user_variable[i] == "switch_output_overview"){
      switch_output_overview = atoi(user_value[i].c_str());
    }

    else if (user_variable[i] == "switch_output_table_order"){
      switch_output_table_order = atoi(user_value[i].c_str());
    }

    else if (user_variable[i] == "switch_output_table_Kfactor"){
      switch_output_table_Kfactor = atoi(user_value[i].c_str());
    }

    else if (user_variable[i] == "switch_output_table_crosssection_Kfactor"){
      switch_output_table_crosssection_Kfactor = atoi(user_value[i].c_str());
    }

    else if (user_variable[i] == "switch_output_table_IS_splitting"){
      switch_output_table_IS_splitting = atoi(user_value[i].c_str());
    }

    else if (user_variable[i] == "output_selection_distribution_table"){
      output_selection_distribution_table.push_back(user_value[i]);
    }

    
    else if (user_variable[i] == "accuracy_normalization"){
      yorder[counter_directory].accuracy_normalization = user_value[i];
    }
  }

  if (phasespace_optimization.size() > 1){
    phasespace_optimization.push_back(".combined");
  }
  for (int i_z = 0; i_z < phasespace_optimization.size(); i_z++){
    logger << LOG_INFO << "phasespace_optimization[" << setw(2) << i_z << "] = " << phasespace_optimization[i_z] << endl;
  }
  
  for (int i_c = 0; i_c < yorder.size(); i_c++){
    if (yorder[i_c].accuracy_normalization == ""){yorder[i_c].accuracy_no_normalization = i_c;}
    else {
      for (int j_c = 0; j_c < yorder.size(); j_c++){
	if (yorder[i_c].accuracy_normalization == yorder[j_c].resultdirectory){yorder[i_c].accuracy_no_normalization = j_c; break;}
      }
    }
  }
  for (int i_c = 0; i_c < yorder.size(); i_c++){
    if (yorder[i_c].accuracy_relative == 0.){yorder[i_c].accuracy_relative = 0.001;}
  }



  for (int i_c = 0; i_c < yorder.size(); i_c++){
    for (int j_c = 0; j_c < yorder[i_c].contribution_file.size(); j_c++){
      int flag = list_contribution_file.size();
      for (int i_l = 0; i_l < list_contribution_file.size(); i_l++){
	if (yorder[i_c].contribution_file[j_c] == list_contribution_file[i_l]){flag = i_l; break;}
      }
      if (flag == list_contribution_file.size()){
	list_contribution_file.push_back(yorder[i_c].contribution_file[j_c]);
	mapping_contribution_file[yorder[i_c].contribution_file[j_c]] = flag;
      }
    }
  }

  logger << LOG_DEBUG << "final_resultdirectory = " << final_resultdirectory << endl;
  /*
  for (int i_l = 0; i_l < resultdirectory.size(); i_l++){
    logger << LOG_DEBUG << "resultdirectory[" << i_l << "] = " << resultdirectory[i_l] << endl;
  }

  for (int i_c = 0; i_c < contribution_file.size(); i_c++){
    for (int j_c = 0; j_c < contribution_file[i_c].size(); j_c++){
      logger << LOG_DEBUG << "contribution_file[" << i_c << "][" << j_c << "] = " << contribution_file[i_c][j_c] << endl;
    }
  }
  for (int i_o = 0; i_o < order_combination.size(); i_o++){
    for (int i_c = 0; i_c < order_combination[i_o].size(); i_c++){
      for (int j_c = 0; j_c < order_combination[i_o][i_c].size(); j_c++){
	logger << LOG_DEBUG << "order_combination[" << i_o << "][" << i_c << "][" << j_c << "] = " << order_combination[i_o][i_c][j_c] << endl;
      }
    }
  }
  */
  for (int i_l = 0; i_l < list_contribution_file.size(); i_l++){
    logger << LOG_DEBUG << "list_contribution_file[" << i_l << "] = " << list_contribution_file[i_l] << endl;
  }
  /*
  for (int i_l = 0; i_l < accuracy_relative.size(); i_l++){
    logger << LOG_DEBUG << "accuracy_relative[" << setw(2) << i_l << "] = " << setw(2) << accuracy_relative[i_l] << " " << setw(16) << accuracy_normalization[i_l] << " (" << accuracy_no_normalization[i_l] << ")" << endl;
  }
  */

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_generic::readin_infile_scaleband(){
  Logger logger("summary_generic::readin_infile_scaleband");
  logger << LOG_DEBUG << "called" << endl;

  string vs;
  string s0;
  vector<string> vs0;
  vector<int> vi0;

  string filename;
  string rundirectory;

  vector<string> readin;
  char LineBuffer[128];

  logger << LOG_DEBUG << "infilename_scaleband = " << infilename_scaleband << endl;
  ifstream in_file(infilename_scaleband.c_str());
  readin.clear();
  while (in_file.getline(LineBuffer, 1024)){readin.push_back(LineBuffer);}
  logger << LOG_DEBUG << "readin.size() = " << readin.size() << endl;
  vector<string> user_variable;
  vector<string> user_value;
  int user_counter = -1;
  for (int i = 0; i < readin.size(); i++){
    logger << LOG_DEBUG << "readin[" << setw(3) << i << "] = " << readin[i] << endl;

    if (readin[i].size() != 0){
      if (readin[i][0] != '/' && readin[i][0] != '#' && readin[i][0] != '%'){
	int start = 0;
	user_counter++;
	user_variable.push_back(s0);
	user_value.push_back(s0);
	for (int j = 0; j < readin[i].size(); j++){
	  if (start == 0 || start == 1){
	    if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 0){}
	    else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	      user_variable[user_counter].push_back(readin[i][j]);
	      if (start != 1){start = 1;}
	    }
	    else {start++;}
	  }
	  else if (start == 2){
	    if (readin[i][j] == '='){start++;}
	    else if ((readin[i][j] == ' ') || (readin[i][j] == char(9))){}
	    else {logger << LOG_FATAL << readin[i][j] << "    Incorrect input in line " << i + 1 << endl; exit(1);}
	  }
	  else if (start == 3 || start == 4){
	    if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 3){}
	    else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	      user_value[user_counter].push_back(readin[i][j]);
	      if (start != 4){start = 4;}
	    }
	    else {start++;}
	  }
	  else {break;}
	}
	if (start == 0){
	  user_counter--;
	  user_variable.erase(user_variable.end() - 1, user_variable.end());
	  user_value.erase(user_value.end() - 1, user_value.end());
	}
      }
    }
  }

  logger << LOG_DEBUG << "user_variable.size() = " << user_variable.size() << endl;
  logger << LOG_DEBUG << "user_value.size() = " << user_value.size() << endl;

  string outpath;
  int counter_contribution = -1;
  for (int i = 0; i < user_variable.size(); i++){
    logger << LOG_DEBUG << "user_variable[" << i << "] = " << user_variable[i] << endl;
    logger << LOG_DEBUG << "user_value[" << i << "] = " << user_value[i] << endl;

    if (user_variable[i] == "outpath"){
      counter_contribution++;
      outpath = user_value[i];
      //      system_execute(logger, "mkdir " + final_resultdirectory + "/" + outpath);
      inpath_scalecentral.push_back("");
      inpath_scalevariation.push_back(vector<string> (0));
      outpath_scalecentral.push_back(outpath + "/scale.central");
      outpath_scaleband.push_back(outpath + "/scale.band");
      outpath_scalemax.push_back(outpath + "/scale.max");
      outpath_scalemin.push_back(outpath + "/scale.min");
      outpath_scaleplusmax.push_back(outpath + "/scale.plusmax");
      outpath_scaleplusmin.push_back(outpath + "/scale.plusmin");
      outpath_scaleminusmax.push_back(outpath + "/scale.minusmax");
      outpath_scaleminusmin.push_back(outpath + "/scale.minusmin");
    }
    else if (user_variable[i] == "inpath_scalecentral"){inpath_scalecentral[counter_contribution] = user_value[i];}
    else if (user_variable[i] == "inpath_scalevariation"){inpath_scalevariation[counter_contribution].push_back(user_value[i]);}
    else if (user_variable[i] == "outpath_scalecentral"){outpath_scalecentral[counter_contribution] = user_value[i];}
    else if (user_variable[i] == "outpath_scaleband"){outpath_scaleband[counter_contribution] = user_value[i];}
    else if (user_variable[i] == "outpath_scalemax"){outpath_scalemax[counter_contribution] = user_value[i];}
    else if (user_variable[i] == "outpath_scalemin"){outpath_scalemin[counter_contribution] = user_value[i];}
    else if (user_variable[i] == "outpath_scaleplusmax"){outpath_scaleplusmax[counter_contribution] = user_value[i];}
    else if (user_variable[i] == "outpath_scaleplusmin"){outpath_scaleplusmin[counter_contribution] = user_value[i];}
    else if (user_variable[i] == "outpath_scaleminusmax"){outpath_scaleminusmax[counter_contribution] = user_value[i];}
    else if (user_variable[i] == "outpath_scaleminusmin"){outpath_scaleminusmin[counter_contribution] = user_value[i];}
    else {logger << LOG_WARN << "Input is not known." << endl;}
  }
    
  logger << LOG_DEBUG << "finished" << endl;
}




void summary_generic::determine_scaleband(){
  Logger logger("summary_generic::determine_scaleband");
  logger << LOG_DEBUG << "called" << endl;

  vector<string> subgroup;
  //  int plotmode = 0;

  logger << LOG_DEBUG << "final_resultdirectory = " << final_resultdirectory << endl;

  readin_infile_scaleband();

  logger << LOG_DEBUG_VERBOSE << "osi_value_qTcut_distribution.size() = " << osi_value_qTcut_distribution.size() << endl;

  scaleband_variable.resize(outpath_scaleband.size(), vector<vector<vector<double> > > (oset.extended_distribution.size(), vector<vector<double> > ()));
  scaleband_central_result.resize(outpath_scaleband.size(), vector<vector<vector<vector<double> > > > (oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ())));
  scaleband_central_deviation.resize(outpath_scaleband.size(), vector<vector<vector<vector<double> > > > (oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ())));
  scaleband_minimum_result.resize(outpath_scaleband.size(), vector<vector<vector<vector<double> > > > (oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ())));
  scaleband_minimum_deviation.resize(outpath_scaleband.size(), vector<vector<vector<vector<double> > > > (oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ())));
  scaleband_maximum_result.resize(outpath_scaleband.size(), vector<vector<vector<vector<double> > > > (oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ())));
  scaleband_maximum_deviation.resize(outpath_scaleband.size(), vector<vector<vector<vector<double> > > > (oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ())));
  
  for (int i_sb = 0; i_sb < outpath_scaleband.size(); i_sb++){
    logger << LOG_INFO << "Determination of scaleband for " << outpath_scaleband[i_sb] << endl;
    // needs to be shifted if used outside...
    /*
    vector<vector<vector<double> > > scaleband_variable(oset.extended_distribution.size(), vector<vector<double> > ());
    vector<vector<vector<vector<double> > > > scaleband_central_result(oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ()));
    vector<vector<vector<vector<double> > > > scaleband_central_deviation(oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ()));
    vector<vector<vector<vector<double> > > > scaleband_minimum_result(oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ()));
    vector<vector<vector<vector<double> > > > scaleband_minimum_deviation(oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ()));
    vector<vector<vector<vector<double> > > > scaleband_maximum_result(oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ()));
    vector<vector<vector<vector<double> > > > scaleband_maximum_deviation(oset.extended_distribution.size(), vector<vector<vector<double> > > (yorder.size(), vector<vector<double> > ()));
    */
    int x_v = -1;
    for (int i_v = 0; i_v < inpath_scalevariation[i_sb].size(); i_v++){
      if (inpath_scalecentral[i_sb] == inpath_scalevariation[i_sb][i_v]){x_v = i_v; break;}
    }
    if (x_v == -1){
      x_v = inpath_scalevariation[i_sb].size();
      inpath_scalevariation[i_sb].push_back(inpath_scalecentral[i_sb]);
    }
    
    system_execute(logger, "mkdir " + final_resultdirectory + "/" + outpath_scaleband[i_sb]);
    system_execute(logger, "mkdir " + final_resultdirectory + "/" + outpath_scalecentral[i_sb]);
    system_execute(logger, "mkdir " + final_resultdirectory + "/" + outpath_scalemax[i_sb]);
    system_execute(logger, "mkdir " + final_resultdirectory + "/" + outpath_scalemin[i_sb]);
   
    for (int x_q = 0; x_q < osi_value_qTcut_distribution.size() + 1; x_q++){
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
      
      if (x_q != osi_value_qTcut_distribution.size()){
	system_execute(logger, "mkdir " + final_resultdirectory + "/" + outpath_scaleband[i_sb] + directory_qTcut);
	system_execute(logger, "mkdir " + final_resultdirectory + "/" + outpath_scalecentral[i_sb] + directory_qTcut);
	system_execute(logger, "mkdir " + final_resultdirectory + "/" + outpath_scalemax[i_sb] + directory_qTcut);
	system_execute(logger, "mkdir " + final_resultdirectory + "/" + outpath_scalemin[i_sb] + directory_qTcut);
      }

     

      string infilename;
      
      char LineBuffer[1024];
      string s0;
      vector<string> vs0, vs1(1);
      logger << LOG_DEBUG << "oset.extended_distribution.size() = " << oset.extended_distribution.size() << endl;
      for (int i_d = 0; i_d < oset.extended_distribution.size(); i_d++){
	logger << LOG_INFO << "Determination of scaleband for " << outpath_scaleband[i_sb] << " : " << oset.extended_distribution[i_d].xdistribution_name << endl;

	logger << LOG_DEBUG << "switch_output_distribution[" << i_d << "] = " << switch_output_distribution[i_d] << endl;
	if (!switch_output_distribution[i_d]){continue;}
	logger << LOG_DEBUG << "oset.extended_distribution[" << i_d << "].xdistribution_name = " << oset.extended_distribution[i_d].xdistribution_name << endl;
	
	for (int i_o = 0; i_o < yorder.size(); i_o++){
	  logger << LOG_DEBUG << "(!yorder[i_o].active_qTcut && x_q < osi_value_qTcut_distribution.size()) = " << (!yorder[i_o].active_qTcut && x_q < osi_value_qTcut_distribution.size()) << endl;
	  logger << LOG_DEBUG << "osi_value_qTcut_distribution.size() = " << osi_value_qTcut_distribution.size() << endl;
	  logger << LOG_DEBUG << "yorder[" << i_o << "].active_qTcut = " << yorder[i_o].active_qTcut << endl;

	  if (!yorder[i_o].active_qTcut && x_q < osi_value_qTcut_distribution.size()){continue;}
	  
	  vector<string> name_plot;
	  string temp_sdd = "." + oset.extended_distribution[i_d].xdistribution_name + ".." + yorder[i_o].resultdirectory;
	  
	  name_plot.push_back("plot" + temp_sdd + ".dat");
	  name_plot.push_back("norm" + temp_sdd + ".dat");
	  
	  if (i_d >= osi_dat.size()){
	    int i_ddd = i_d - osi_dat.size();
	    
	    // recombined distributions (essentially for validation)
	    stringstream name_rec;
	    name_rec << ".rec." << osi_dddat[i_ddd].distribution_2.xdistribution_name << ".from." << osi_dddat[i_ddd].name << ".." << yorder[i_o].resultdirectory;
	    
	    name_plot.push_back("plot" + name_rec.str() + ".dat");
	    name_plot.push_back("norm" + name_rec.str()  + ".dat");
	    
	    for (int i_b1 = 0; i_b1 < osi_dddat[i_ddd].distribution_1.n_bins; i_b1++){
	      stringstream name_split_bin;
	      name_split_bin << "_" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1] << "-" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1 + 1];
	      stringstream name_split_lt;
	      name_split_lt << "_lt" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1];
	      stringstream name_split_ge;
	      name_split_ge << "_ge" << osi_dddat[i_ddd].distribution_1.bin_edge[i_b1];
	      
	      name_plot.push_back("norm.norm.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".." + yorder[i_o].resultdirectory + ".dat");
	      name_plot.push_back("norm.plot.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".." + yorder[i_o].resultdirectory + ".dat");
	      name_plot.push_back("plot.norm.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".." + yorder[i_o].resultdirectory + ".dat");
	      name_plot.push_back("plot.plot.split." + osi_dddat[i_ddd].name + name_split_bin.str() + ".." + yorder[i_o].resultdirectory + ".dat");
	      
	      name_plot.push_back("norm.two." + osi_dddat[i_ddd].name + name_split_lt.str() + ".." + yorder[i_o].resultdirectory + ".dat");
	      name_plot.push_back("norm.two." + osi_dddat[i_ddd].name + name_split_ge.str() + ".." + yorder[i_o].resultdirectory + ".dat");
	      name_plot.push_back("plot.two." + osi_dddat[i_ddd].name + name_split_lt.str() + ".." + yorder[i_o].resultdirectory + ".dat");
	      name_plot.push_back("plot.two." + osi_dddat[i_ddd].name + name_split_ge.str() + ".." + yorder[i_o].resultdirectory + ".dat");
	    }
	  }
	  
	  int n_output_version = name_plot.size();
	  logger << LOG_DEBUG << "n_output_version = " << n_output_version << endl;
	  
	  scaleband_variable[i_sb][i_d].resize(n_output_version);

	  logger << LOG_DEBUG << "scaleband_variable[" << i_sb << "].size() = " << scaleband_variable[i_sb].size() << endl;
	  logger << LOG_DEBUG << "scaleband_variable[" << i_sb << "][" << i_d << "].size() = " << scaleband_variable[i_sb][i_d].size() << endl;

	  scaleband_central_result[i_sb][i_d][i_o].resize(n_output_version);
	  scaleband_central_deviation[i_sb][i_d][i_o].resize(n_output_version);
	  scaleband_minimum_result[i_sb][i_d][i_o].resize(n_output_version);
	  scaleband_minimum_deviation[i_sb][i_d][i_o].resize(n_output_version);
	  scaleband_maximum_result[i_sb][i_d][i_o].resize(n_output_version);
	  scaleband_maximum_deviation[i_sb][i_d][i_o].resize(n_output_version);
	  
	  for (int i_x = 0; i_x < n_output_version; i_x++){
	    logger << LOG_DEBUG << "name_plot[" << setw(3) << i_x << "] = " << name_plot[i_x] << endl;
	    int switch_exist = 0;
	    
	    vector<vector<vector<double> > > xresult, xdeviation, xdeviation2;
	    vector<double> variable;
	    vector<double> centralresult;
	    vector<double> centraldeviation;
	    for (int i_v = 0; i_v < inpath_scalevariation[i_sb].size(); i_v++){
	      vector<string> readin;
	      vector<vector<string> > readin2;
	      
	      infilename = final_resultdirectory + "/" + inpath_scalevariation[i_sb][i_v] + directory_qTcut + "/" + name_plot[i_x];
	      logger << LOG_DEBUG_VERBOSE << "i_d = " << setw(2) << i_d << "   i_v = " << setw(2) << i_v << "   i_o = " << setw(2) << i_o << "   infilename = " << infilename << endl;
	      
	      ifstream inf(infilename.c_str());
	      readin = vs0;
	      while (inf.getline(LineBuffer, 1024)){readin.push_back(LineBuffer);}
	      logger << LOG_DEBUG << infilename << "   readin.size() = " << readin.size() << endl;
	      if (readin.size() == 0){continue;} 
	      else {switch_exist = 1;}
	      for (int i_b = 0; i_b < readin.size(); i_b++){
		if (readin[i_b][0] == char(35)){}
		else {
		  readin2.push_back(vs0);
		  int start = 0;
		  for (int j = 0; j < readin[i_b].size(); j++){
		    if (readin[i_b][j] != ' '){
		      if (start == 0){start = 1; readin2[readin2.size() - 1].push_back(s0); readin2[readin2.size() - 1][readin2[readin2.size() - 1].size() - 1] = readin[i_b][j];}
		       else {readin2[readin2.size() - 1][readin2[readin2.size() - 1].size() - 1] = readin2[readin2.size() - 1][readin2[readin2.size() - 1].size() - 1] + readin[i_b][j];}
		    }
		    else if ((readin[i_b][j] == ' ') && (start == 1)){start = 0;}
		  }
		}
		if (i_b == 0 || i_b == readin.size() - 1){logger << LOG_DEBUG_VERBOSE << "i_b = " << setw(2) << i_b << "   (" << readin.size() << " - " << readin2.size() << ")" << endl;}
	      }
	      
	      if (i_v == 0){variable.resize(readin2.size());}
	      vector<vector<double> > data(readin2.size());
	      vector<vector<double> > result(readin2.size()), deviation(readin2.size()), deviation2(readin2.size());
	      
	      
	      logger << LOG_DEBUG_VERBOSE << "i_v = " << i_v << "   readin2.size() = " << readin2.size() << endl;
	      
	      for (int i_b = 0; i_b < readin2.size(); i_b++){
		data[i_b].resize(readin2[i_b].size());
		result[i_b].resize((readin2[i_b].size() - 1) / 2);
		deviation[i_b].resize((readin2[i_b].size() - 1) / 2);
		deviation2[i_b].resize((readin2[i_b].size() - 1) / 2);
	      }
	      
	      for (int i_b = 0; i_b < readin2.size(); i_b++){
		for (int j = 0; j < readin2[i_b].size(); j++){
		  data[i_b][j] = atof(readin2[i_b][j].c_str());
		}
	      }
	      for (int i_b = 0; i_b < result.size(); i_b++){
		if (i_v == 0){variable[i_b] = data[i_b][0];}
		for (int j = 0; j < result[i_b].size(); j++){
		  result[i_b][j] = data[i_b][j * 2 + 1];
		  deviation[i_b][j] = data[i_b][j * 2 + 2];
		  deviation2[i_b][j] = pow(deviation[i_b][j], 2);
		}
	      }
	      
	      if (x_v == i_v){
		centralresult.resize(result.size());
		centraldeviation.resize(deviation.size());
		for (int i_b = 0; i_b < result.size(); i_b++){
		  centralresult[i_b] = result[i_b][0];
		  centraldeviation[i_b] = deviation[i_b][0];
		}
	      }
	      xresult.push_back(result);
	      xdeviation.push_back(deviation);
	      xdeviation2.push_back(deviation2);
	    }
	    if (!switch_exist){continue;}





	    vector<double> minresult(xresult[0].size(), 0.);
	    vector<double> maxresult(xresult[0].size(), 0.);
	    vector<double> mindeviation(xdeviation[0].size(), 0.);
	    vector<double> maxdeviation(xdeviation[0].size(), 0.);
	    for (int i_v = 0; i_v < xresult[0].size(); i_v++){
	      logger << LOG_DEBUG_VERBOSE << "i_v = " << i_v << "   xresult.size() = " << xresult.size() << endl;
	      
	      minresult[i_v] = 1.e99;
	      maxresult[i_v] = -1.e99;
	      for (int i_c = 0; i_c < xresult.size(); i_c++){
		if (xresult[i_c][i_v][0] < minresult[i_v]){
		  minresult[i_v] = xresult[i_c][i_v][0];
		  mindeviation[i_v] = xdeviation[i_c][i_v][0];
		}
		if (xresult[i_c][i_v][0] > maxresult[i_v]){
		  maxresult[i_v] = xresult[i_c][i_v][0];
		  maxdeviation[i_v] = xdeviation[i_c][i_v][0];
		}
	      }
	      //	if (plotmode == 1){if (minresult[i_b] < 0.){minresult[i_b] = 0.;}}
	      logger << LOG_DEBUG_VERBOSE << "i_v = " << i_v << "   xresult.size() = " << xresult.size() << " done." << endl;
	    }


	    scaleband_variable[i_sb][i_d][i_x] = variable;
	    scaleband_central_result[i_sb][i_d][i_o][i_x] = centralresult;
	    scaleband_central_deviation[i_sb][i_d][i_o][i_x] = centraldeviation;
	    scaleband_minimum_result[i_sb][i_d][i_o][i_x] = minresult;
	    scaleband_minimum_deviation[i_sb][i_d][i_o][i_x] = mindeviation;
	    scaleband_maximum_result[i_sb][i_d][i_o][i_x] = maxresult;
	    scaleband_maximum_deviation[i_sb][i_d][i_o][i_x] = maxdeviation;

	    
	    ofstream outf;
	    
	    string outfilename_scaleband_matrix = final_resultdirectory + "/" + outpath_scaleband[i_sb] + directory_qTcut + "/" + name_plot[i_x];
	    string outfilename_scaleband = final_resultdirectory + "/" + outpath_scaleband[i_sb] + directory_qTcut + "/fill" + name_plot[i_x];
	    string outfilename_scaleband_ref = final_resultdirectory + "/" + outpath_scaleband[i_sb] + directory_qTcut + "/reffill" + name_plot[i_x];
	    string outfilename_scalecentral = final_resultdirectory + "/" + outpath_scalecentral[i_sb] + directory_qTcut + "/" + name_plot[i_x];
	    string outfilename_scalemax = final_resultdirectory + "/" + outpath_scalemax[i_sb] + directory_qTcut + "/" + name_plot[i_x];
	    string outfilename_scalemin = final_resultdirectory + "/" + outpath_scalemin[i_sb] + directory_qTcut + "/" + name_plot[i_x];
	    
	    logger << LOG_DEBUG << "outfilename_scaleband_matrix = " << outfilename_scaleband_matrix << endl;
	    logger << LOG_DEBUG << "outfilename_scaleband = " << outfilename_scaleband << endl;
	    logger << LOG_DEBUG << "outfilename_scalecentral = " << outfilename_scalecentral << endl;
	    logger << LOG_DEBUG << "outfilename_scalemax = " << outfilename_scalemax << endl;
	    logger << LOG_DEBUG << "outfilename_scalemin = " << outfilename_scalemin << endl;
	    
	    /////////////////
	    //  scaleband  //
	    /////////////////
	    outf.open(outfilename_scaleband_matrix.c_str(), ofstream::out | ofstream::trunc);
	    for (int i_b = 0; i_b < centralresult.size(); i_b++){
	      outf << setw(16) << setprecision(8) << variable[i_b] << "   " 
		   << setw(16) << setprecision(8) << centralresult[i_b] << "   " << setw(16) << setprecision(8) << centraldeviation[i_b] 
		   << setw(16) << setprecision(8) << minresult[i_b] << "   " << setw(16) << setprecision(8) << mindeviation[i_b] 
		   << setw(16) << setprecision(8) << maxresult[i_b] << "   " << setw(16) << setprecision(8) << maxdeviation[i_b] 
		   << endl;
	    }
	    outf.close();
	    
	    // outfilename_scaleband_matrix contains all the output needed for using matrix plot scripts !!! 
	    
	    outf.open(outfilename_scaleband.c_str(), ofstream::out | ofstream::trunc);
	    for (int i_b = 0; i_b < minresult.size(); i_b++){

	      if (centralresult[i_b] < 0.){centralresult[i_b] = 0.;}
	      if (minresult[i_b] < 0.){minresult[i_b] = 0.;}
	      if (maxresult[i_b] < 0.){maxresult[i_b] = 0.;}
	      
	      outf << setw(16) << setprecision(8) << variable[i_b] << "   " 
		   << setw(16) << setprecision(8) << centralresult[i_b] << "   " << setw(16) << setprecision(8) << centraldeviation[i_b] 
		   << setw(16) << setprecision(8) << minresult[i_b] << "   " << setw(16) << setprecision(8) << mindeviation[i_b] 
		   << setw(16) << setprecision(8) << maxresult[i_b] << "   " << setw(16) << setprecision(8) << maxdeviation[i_b] 
		   << endl;
	      //	       outf << setw(16) << setprecision(8) << variable[i_b] << "   " << setw(16) << minresult[i_b] << setw(16) << maxresult[i_b] << endl;
	      
	      if (i_b < minresult.size() - 1){
		outf << setw(16) << setprecision(8) << variable[i_b + 1] << "   " 
		     << setw(16) << setprecision(8) << centralresult[i_b] << "   " << setw(16) << setprecision(8) << centraldeviation[i_b] 
		     << setw(16) << setprecision(8) << minresult[i_b] << "   " << setw(16) << setprecision(8) << mindeviation[i_b] 
		     << setw(16) << setprecision(8) << maxresult[i_b] << "   " << setw(16) << setprecision(8) << maxdeviation[i_b] 
		     << endl;
		//		 outf << setw(16) << setprecision(8) << variable[i_b + 1] << "   " << setw(16) << minresult[i_b] << setw(16) << maxresult[i_b] << endl;
	      }
	    }
	    outf.close();
	    
	    /////////////////////
	    //  scaleband-ref  //
	    /////////////////////
	    outf.open(outfilename_scaleband_ref.c_str(), ofstream::out | ofstream::trunc);
	    for (int i_b = 0; i_b < centralresult.size(); i_b++){
	      outf << setw(16) << setprecision(8) << variable[i_b] << setw(16) << setprecision(8) << centralresult[i_b] << endl;
	      if (i_b < centralresult.size() - 1){
		outf << setw(16) << setprecision(8) << variable[i_b + 1] << setw(16) << setprecision(8) << centralresult[i_b] << endl;
	      }
	    }
	    outf.close();
	    
	    ////////////////////
	    //  scalecentral  //
	    ////////////////////
	    /*
	    outf.open(outfilename_scalecentral.c_str(), ofstream::out | ofstream::trunc);
	    for (int i_b = 0; i_b < centralresult.size(); i_b++){
	      outf << setw(16) << setprecision(8) << variable[i_b] << setw(16) << setprecision(8) << centralresult[i_b] << "   " << setw(16) << setprecision(8) << centraldeviation[i_b] << endl;
	      //  if (i < centralresult.size() - 1){outf << setw(16) << setprecision(8) << variable[i + 1] << setw(16) << setprecision(8) << centralresult[i_b] << "   " << setw(16) << setprecision(8) << centraldeviation[i_b] << endl;}
	    }
	    outf.close();
	    */
	    outf.open(outfilename_scalecentral.c_str(), ofstream::out | ofstream::trunc);
	    for (int i_b = 0; i_b < centralresult.size() - 1; i_b++){
	      outf << setw(16) << setprecision(8) << variable[i_b] << "   " 
		   << setw(16) << setprecision(8) << variable[i_b + 1] << "   " 
		   << setw(16) << setprecision(8) << centralresult[i_b] << "   " << setw(16) << setprecision(8) << centraldeviation[i_b] 
		   << setw(16) << setprecision(8) << minresult[i_b] << "   " << setw(16) << setprecision(8) << mindeviation[i_b] 
		   << setw(16) << setprecision(8) << maxresult[i_b] << "   " << setw(16) << setprecision(8) << maxdeviation[i_b] 
		   << endl;
	    }
	    outf.close();

	    ////////////////
	    //  scalemax  //
	    ////////////////
	    outf.open(outfilename_scalemax.c_str(), ofstream::out | ofstream::trunc);
	    for (int i_b = 0; i_b < maxresult.size(); i_b++){
	      outf << setw(16) << setprecision(8) << variable[i_b] << setw(16) << setprecision(8) << maxresult[i_b] << "   " << setw(16) << setprecision(8) << maxdeviation[i_b] << endl;
	      //  if (i < maxresult.size() - 1){outf << setw(16) << setprecision(8) << variable[i + 1] << setw(16) << setprecision(8) << maxresult[i_b] << "   " << setw(16) << setprecision(8) << maxdeviation[i_b] << endl;}
	    }
	    outf.close();
	    
	    ////////////////
	    //  scalemin  //
	    ////////////////
	    outf.open(outfilename_scalemin.c_str(), ofstream::out | ofstream::trunc);
	    for (int i_b = 0; i_b < minresult.size(); i_b++){
	      outf << setw(16) << setprecision(8) << variable[i_b] << setw(16) << setprecision(8) << minresult[i_b] << "   " << setw(16) << setprecision(8) << mindeviation[i_b] << endl;
	      //  if (i < minresult.size() - 1){outf << setw(16) << setprecision(8) << variable[i + 1] << setw(16) << setprecision(8) << minresult[i_b] << "   " << setw(16) << setprecision(8) << mindeviation[i_b] << endl;}
	    }
	    outf.close();
	  }
	}

      }
    }

  }

  if (switch_output_table_order){output_distribution_table_order();}
  if (switch_output_table_Kfactor){output_distribution_table_Kfactor();}
  if (switch_output_table_crosssection_Kfactor){output_distribution_table_crosssection_Kfactor();}
  if (switch_output_table_IS_splitting){output_distribution_table_IS_splitting();}
 
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void summary_generic::determine_runtime(){
  Logger logger("summary_generic::determine_runtime");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  string filename;

  //////////////////////////////////////////////////////
  //  Output information about extrapolated runtimes  //
  //////////////////////////////////////////////////////

  double runtime_extrapolated_complete = 0.;
  for (int i_l = 0; i_l < xlist.size(); i_l++){
    logger << LOG_INFO << "xlist[" << i_l << "].xcontribution[0].xsubprocess[0].extrapolated_runtime = " << xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_runtime << endl;
  }
  for (int i_l = 0; i_l < xlist.size(); i_l++){runtime_extrapolated_complete += xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_runtime;}
  // check double-counting (of runtimes) around here !!!
  for (int i_l = 0; i_l < xlist.size(); i_l++){
    logger << LOG_INFO << "xlist[" << i_l << "].xcontribution[0].xsubprocess[0].extrapolated_runtime = " << xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_runtime << endl;
  }
  logger << LOG_INFO << "runtime_extrapolated_complete = " << runtime_extrapolated_complete << endl;

  ofstream outfile_extrapolated;
  filename = final_resultdirectory + "/info.extrapolated.runtime.txt";
  outfile_extrapolated.open(filename.c_str(), ofstream::out | ofstream::trunc);  
  outfile_extrapolated << "Overview over contributions   " << setprecision(8) << setw(16) << runtime_extrapolated_complete / 24 / 60 / 60 << " d" << endl;
  for (int i_s = 0; i_s < 27; i_s++){outfile_extrapolated << "=";} outfile_extrapolated << endl;
  outfile_extrapolated << endl;
  outfile_extrapolated << left << setw(25) << "folder" << setw(20) << "" << setw(21) << right << "runtime" << "   " << setw(12) << "" << "   " << setw(16) << "deviation" << setw(16) << "norm" << setw(16) << "norm_err" << endl;
  outfile_extrapolated << endl;

  outfile_extrapolated << left << setw(25) << "complete";
  outfile_extrapolated << left << setw(20) << "";
  outfile_extrapolated << time_hms_from_double(runtime_extrapolated_complete) << "   ";

  outfile_extrapolated << endl;
  for (int i_s = 0; i_s < 129; i_s++){outfile_extrapolated << "=";} outfile_extrapolated << endl;
  outfile_extrapolated << endl;

  for (int i_l = 0; i_l < xlist.size(); i_l++){
    for (int i_c = 0; i_c < xlist[i_l].xcontribution.size(); i_c++){
      //      logger << LOG_INFO << "xlist[" << i_l << "].xcontribution[" << i_c << "].xsubprocess[0].extrapolated_runtime = " << xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_runtime << endl;
      outfile_extrapolated << left << setw(25) << xlist[i_l].xcontribution[i_c].infix_path_contribution;
      outfile_extrapolated << left << setw(20) << xlist[i_l].xcontribution[i_c].xsubprocess[0].name;
      outfile_extrapolated << time_hms_from_double(xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_runtime) << "   ";
      outfile_extrapolated << left << setw(16) << "";
      outfile_extrapolated << right << setw(16) << setprecision(8) << xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_normalization_deviation;
      outfile_extrapolated << right << setw(16) << setprecision(8) << xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_normalization;
      outfile_extrapolated << right << setw(16) << setprecision(8) << xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_deviation;

  //      outfile_extrapolated << left << setw(25) << "folder" << setw(20) << "subprocess" << setw(21) << right << "runtime" << "   " << setw(12) << "n_event" << "   " << setw(16) << "deviation" << setw(16) << "norm" << setw(16) << "norm_err" << endl;

      outfile_extrapolated << endl;
      if (i_c == 0){for (int i_s = 0; i_s < 129; i_s++){outfile_extrapolated << "-";} outfile_extrapolated << endl;}
    }
    outfile_extrapolated << endl;
  }
  outfile_extrapolated << endl;
  outfile_extrapolated << "Overview over subprocesses" << endl;
  for (int i_s = 0; i_s < 26; i_s++){outfile_extrapolated << "=";} outfile_extrapolated << endl;
  outfile_extrapolated << endl;
  outfile_extrapolated << left << setw(25) << "folder" << setw(20) << "subprocess" << setw(21) << right << "runtime" << "   " << setw(12) << "n_event" << "   " << setw(16) << "deviation" << setw(16) << "norm" << setw(16) << "norm_err" << endl;
  outfile_extrapolated << endl;
  for (int i_l = 0; i_l < xlist.size(); i_l++){xlist[i_l].output_extrapolated_runtime(outfile_extrapolated);}
  outfile_extrapolated.close();



  ////////////////////////////////////////////////////////////////////////////////
  //  Output information about extrapolated runtimes - input for MATRIX script  //
  ////////////////////////////////////////////////////////////////////////////////

  ofstream outfile_readin_extrapolated;
  filename = "./runtime.dat";
  outfile_readin_extrapolated.open(filename.c_str(), ofstream::out | ofstream::trunc);  
  outfile_readin_extrapolated << left << setw(25) << "# folder" << setw(20) << "subprocess" << setw(12) << right << "runtime(s)" << "   " << setw(12) << "n_event" << "   " << setw(16) << "deviation" << setw(16) << "norm" << setw(16) << "norm_err" << endl;
  for (int i_l = 0; i_l < xlist.size(); i_l++){xlist[i_l].output_readin_extrapolated_runtime(outfile_readin_extrapolated);}
  outfile_readin_extrapolated.close();



  //////////////////////////////////////////////
  //  Output information about used runtimes  //
  //////////////////////////////////////////////

  double runtime_used_complete = 0.;
  for (int i_l = 0; i_l < xlist.size(); i_l++){runtime_used_complete += xlist[i_l].xcontribution[0].xsubprocess[0].used_runtime;}

  ofstream outfile_used;
  filename = final_resultdirectory + "/info.used.runtime.txt";
  outfile_used.open(filename.c_str(), ofstream::out | ofstream::trunc);  
  outfile_used << "Overview over contributions   " << setprecision(8) << setw(16) << runtime_used_complete / 24 / 60 / 60 << " d" << endl;
  for (int i_s = 0; i_s < 27; i_s++){outfile_used << "=";} outfile_used << endl;
  outfile_used << endl;
  outfile_used << left << setw(25) << "folder" << setw(20) << "" << setw(21) << right << "runtime" << endl;
  outfile_used << endl;

  outfile_used << left << setw(25) << "complete";
  outfile_used << left << setw(20) << "";
  outfile_used << time_hms_from_double(runtime_used_complete) << "   ";
  outfile_used << endl;
  for (int i_s = 0; i_s < 66; i_s++){outfile_used << "=";} outfile_used << endl;
  outfile_used << endl;

  for (int i_l = 0; i_l < xlist.size(); i_l++){
    for (int i_c = 0; i_c < xlist[i_l].xcontribution.size(); i_c++){
      stringstream temp_ext_runtime;
      outfile_used << left << setw(25) << xlist[i_l].xcontribution[i_c].infix_path_contribution;
      outfile_used << left << setw(20) << "";//xlist[i_l].xcontribution[i_c].xsubprocess[0].name;
      outfile_used << time_hms_from_double(xlist[i_l].xcontribution[i_c].xsubprocess[0].used_runtime) << "   ";
      outfile_used << endl;
      if (i_c == 0){for (int i_s = 0; i_s < 66; i_s++){outfile_used << "-";} outfile_used << endl;}
    }
    outfile_used << endl;
  }
  outfile_used << endl;
  outfile_used << "Overview over subprocesses" << endl;
  for (int i_s = 0; i_s < 26; i_s++){outfile_used << "=";} outfile_used << endl;
  outfile_used << endl;
  outfile_used << left << setw(25) << "folder" << setw(20) << "subprocess" << setw(16) << right << "runtime" << "   " << setw(12) << "n_event" << "   " << setw(16) << "deviation" << setw(16) << "norm" << setw(16) << "norm_err" << endl;
  outfile_used << endl;
  for (int i_l = 0; i_l < xlist.size(); i_l++){xlist[i_l].output_used_runtime(outfile_used);}
  outfile_used.close();

  /////////////////////////////////////////////////////////////////////
  //  Output for new run directories based on extrapolated runtimes  //
  /////////////////////////////////////////////////////////////////////
  if (switch_generate_rundirectories){

    int n_run_phasespace_optimization = phasespace_optimization.size();
    if (phasespace_optimization.size() > 1){n_run_phasespace_optimization = phasespace_optimization.size() - 1;}

    stringstream list_path;
    list_path << final_resultdirectory + "/newrundir";
    system_execute(logger, "mkdir " + list_path.str());

    string path_newrundir_result;
    path_newrundir_result =  final_resultdirectory + "/newrundir/result";
    system_execute(logger, "mkdir " + path_newrundir_result);
    
    list_path << "/start";
    for (int i_m = 0; i_m < n_run_phasespace_optimization; i_m++){
      system_execute(logger, "mkdir " + list_path.str() + phasespace_optimization[i_m]);
    }
    
    string path_newrundir_result_infile_runs;
    path_newrundir_result_infile_runs = final_resultdirectory + "/newrundir/result/infile.runs";
    for (int i_m = 0; i_m < phasespace_optimization.size(); i_m++){
      system_execute(logger, "mkdir " + path_newrundir_result_infile_runs + phasespace_optimization[i_m]);
    }


    vector<int> counter_complete_process(n_run_phasespace_optimization, 0);
    vector<vector<int> > counter_list_process(n_run_phasespace_optimization, vector<int> (xlist.size() , 0));
    vector<vector<vector<int> > > counter_contribution_process(n_run_phasespace_optimization, vector<vector<int> > (xlist.size()));

    for (int i_m = 0; i_m < n_run_phasespace_optimization; i_m++){
      ofstream out_runscript_complete;
      string runscript_complete_name = final_resultdirectory + "/newrundir" + "/MUNICH.complete" + phasespace_optimization[i_m] + ".start.sh";
      
      logger << LOG_DEBUG << "runscript_complete_name = " << runscript_complete_name << endl;
      out_runscript_complete.open(runscript_complete_name.c_str(), ofstream::out | ofstream::trunc);  
      
      ofstream out_processlist_complete;
      string runlist_complete_name = final_resultdirectory + "/newrundir" + "/list.run.subprocesses.complete" + phasespace_optimization[i_m] + ".txt";
      out_processlist_complete.open(runlist_complete_name.c_str(), ofstream::out | ofstream::trunc);  

      
      //      int counter_complete_process = 0;
      for (int i_l = 0; i_l < xlist.size(); i_l++){
	ofstream out_processlist_list;
	string name_processscript_list = final_resultdirectory + "/newrundir" + "/start" + phasespace_optimization[i_m] + "/list.run.subprocesses." + xlist[i_l].xcontribution[0].infix_order_contribution + ".txt";
	logger << LOG_DEBUG << "name_processscript_list = " << name_processscript_list << endl;
	out_processlist_list.open(name_processscript_list.c_str(), ofstream::out | ofstream::trunc);  
	
	ofstream out_runscript_list;
	string name_runscript_list = final_resultdirectory + "/newrundir/start" + phasespace_optimization[i_m] + "/MUNICH." + xlist[i_l].xcontribution[0].infix_order_contribution + ".start.sh";
	out_runscript_list.open(name_runscript_list.c_str(), ofstream::out | ofstream::trunc);  
	logger << LOG_DEBUG << "name_runscript_list = " << name_runscript_list << endl;
	
	//	int counter_list_process = 0;
	//	vector<int> counter_list_process(n_run_phasespace_optimization, 0);
	counter_contribution_process[i_m][i_l].resize(xlist[i_l].xcontribution.size(), 0);
	for (int i_c = 1; i_c < xlist[i_l].xcontribution.size(); i_c++){
	  //	  vector<int> counter_contribution_process(n_run_phasespace_optimization, 0);
	  //	  int counter_contribution_process = 0;
	  xlist[i_l].xcontribution[i_c].output_extrapolated_runtime_directories(i_m, counter_contribution_process[i_m][i_l][i_c]);
	  //out_runscript_list, out_processlist_list
	  out_processlist_list << right << setw(8) << "" << setw(30) << left << xlist[i_l].xcontribution[i_c].infix_order_contribution + phasespace_optimization[i_m] << "     " << right << setw(5) << counter_contribution_process[i_m][i_l][i_c] << " ( " << setw(5) << xlist[i_l].xcontribution[i_c].max_number_of_jobs_for_one_channel << " )" << endl;
	  out_runscript_list << "./start" + phasespace_optimization[i_m] + "/MUNICH." + xlist[i_l].xcontribution[i_c].infix_order_contribution + ".start.sh" << endl;
	  counter_list_process[i_m][i_l] += counter_contribution_process[i_m][i_l][i_c];
	}
	out_processlist_list << right << setw(8) << "========" << setw(30) << "==============================" << "=====" << "=====" << "===" << "=====" << "==" << endl;
	out_processlist_list << right << setw(8) << "" << setw(30) << left << xlist[i_l].xcontribution[0].infix_order_contribution + phasespace_optimization[i_m] << "     " << right << setw(5) << counter_list_process[i_m][i_l] << endl;
	out_runscript_list.close();
	system_execute(logger, "chmod 700 " + name_runscript_list);
	out_processlist_list.close();

	out_runscript_complete << "./start" + phasespace_optimization[i_m] + "/MUNICH." << xlist[i_l].xcontribution[0].infix_order_contribution << ".start.sh" << endl;
	out_processlist_complete << right << setw(8) << "" << setw(30) << left << xlist[i_l].xcontribution[0].infix_order_contribution + phasespace_optimization[i_m] << "     " << right << setw(5) << counter_list_process[i_m][i_l] << endl;
	counter_complete_process[i_m] += counter_list_process[i_m][i_l];
      }

      out_runscript_complete.close();
      system_execute(logger, "chmod 700 " + runscript_complete_name);
      out_processlist_complete << right << setw(8) << "========" << setw(30) << "==============================" << "=====" << "=====" << endl;
      out_processlist_complete << right << setw(8) << "" << setw(30) << left << "complete" << "     " << right << setw(5) << counter_complete_process[i_m] << endl;
      out_processlist_complete.close();
    }
      
      /*
      // old version
      // result/infile.runs/<xxx>.dat
      ofstream out_old_infile_result_sub;
      string old_infile_result_sub_name = path_newrundir_result_infile_runs + "/old." + xlist[i_l].xcontribution[0].infix_order_contribution + ".dat";
      out_old_infile_result_sub.open(old_infile_result_sub_name.c_str(), ofstream::out | ofstream::trunc);
      write_infile_string(out_old_infile_result_sub, "processname", xlist[i_l].processname);
      out_old_infile_result_sub << endl;
      write_infile_string(out_old_infile_result_sub, "resultdirectory", xlist[i_l].xcontribution[0].infix_order_contribution);
      write_infile_string(out_old_infile_result_sub, "type_perturbative_order", xlist[i_l].type_perturbative_order);
      if (xlist[i_l].type_subtraction_method == ""){write_infile_string(out_old_infile_result_sub, "subtraction_method", "---");}
      else {write_infile_string(out_old_infile_result_sub, "subtraction_method", xlist[i_l].type_subtraction_method);}
      //	  write_infile_string(out_old_infile_result_sub, "contribution_order", generic_parameter.advanced_directory[j_co][i_sm][j_o]);
      //  if (generic_parameter.advanced_directory[j_co][i_sm][j_o].substr(0, 1) == "a"){
      write_infile_int(out_old_infile_result_sub, "contribution_order_alpha_s", xlist[i_l].in_contribution_order_alpha_s);
      write_infile_int(out_old_infile_result_sub, "contribution_order_alpha_e", xlist[i_l].in_contribution_order_alpha_e);
      write_infile_int(out_old_infile_result_sub, "photon_induced", xlist[i_l].photon_induced);
      out_old_infile_result_sub << endl;

      for (int i_c = 1; i_c < xlist[i_l].xcontribution.size(); i_c++){
	write_infile_string(out_old_infile_result_sub, "type_contribution", xlist[i_l].xcontribution[i_c].type_contribution);
	if (xlist[i_l].xcontribution[i_c].type_correction == ""){write_infile_string(out_old_infile_result_sub, "type_correction", "---");}
	else {write_infile_string(out_old_infile_result_sub, "type_correction", xlist[i_l].xcontribution[i_c].type_correction);}
	write_infile_int(out_old_infile_result_sub, "interference", xlist[i_l].xcontribution[i_c].interference);
	for (int i_r = 0; i_r < xlist[i_l].xcontribution[i_c].max_number_of_jobs_for_one_channel; i_r++){
	  stringstream temp_ss;
	  temp_ss << xlist[i_l].xcontribution[i_c].infix_path_contribution << "/" << "run." << i_r;
	  write_infile_string(out_old_infile_result_sub, "directory", temp_ss.str());
	}
	out_old_infile_result_sub << endl;
      }
            
      out_old_infile_result_sub.close();
      */





      // new version
      // result/infile.runs/<xxx>.dat

    for (int i_m = 0; i_m < phasespace_optimization.size(); i_m++){
      for (int i_l = 0; i_l < xlist.size(); i_l++){
	ofstream out_infile_result_sub;
	string infile_result_sub_name = path_newrundir_result_infile_runs + phasespace_optimization[i_m] + "/" + xlist[i_l].xcontribution[0].infix_order_contribution + ".dat";
	out_infile_result_sub.open(infile_result_sub_name.c_str(), ofstream::out | ofstream::trunc);
	write_infile_string(out_infile_result_sub, "processname", xlist[i_l].processname);
	out_infile_result_sub << endl;
	write_infile_string(out_infile_result_sub, "resultdirectory", xlist[i_l].xcontribution[0].infix_order_contribution);
	write_infile_string(out_infile_result_sub, "type_perturbative_order", xlist[i_l].type_perturbative_order);
	if (xlist[i_l].type_subtraction_method == ""){write_infile_string(out_infile_result_sub, "subtraction_method", "---");}
	else {write_infile_string(out_infile_result_sub, "subtraction_method", xlist[i_l].type_subtraction_method);}
	//	  write_infile_string(out_infile_result_sub, "contribution_order", generic_parameter.advanced_directory[j_co][i_sm][j_o]);
	//  if (generic_parameter.advanced_directory[j_co][i_sm][j_o].substr(0, 1) == "a"){
	write_infile_int(out_infile_result_sub, "contribution_order_alpha_s", xlist[i_l].in_contribution_order_alpha_s);
	write_infile_int(out_infile_result_sub, "contribution_order_alpha_e", xlist[i_l].in_contribution_order_alpha_e);
	write_infile_int(out_infile_result_sub, "photon_induced", xlist[i_l].photon_induced);
	out_infile_result_sub << endl;
	
	for (int i_c = 1; i_c < xlist[i_l].xcontribution.size(); i_c++){
	  write_infile_string(out_infile_result_sub, "type_contribution", xlist[i_l].xcontribution[i_c].type_contribution);
	  if (xlist[i_l].xcontribution[i_c].type_correction == ""){write_infile_string(out_infile_result_sub, "type_correction", "---");}
	  else {write_infile_string(out_infile_result_sub, "type_correction", xlist[i_l].xcontribution[i_c].type_correction);}
	  write_infile_int(out_infile_result_sub, "interference", xlist[i_l].xcontribution[i_c].interference);
	  
	  if (phasespace_optimization[i_m] == ".combined"){
	    for (int j_z = 0; j_z < phasespace_optimization.size() - 1; j_z++){
	      stringstream temp_ss;
	      temp_ss << "start" << phasespace_optimization[j_z] << "/runlist." << xlist[i_l].xcontribution[i_c].infix_order_contribution << ".dat";
	      string temp_listname = "directory_list";
	      if (phasespace_optimization[j_z] != ""){
		temp_listname = temp_listname + " " + phasespace_optimization[j_z].substr(1, phasespace_optimization[j_z].size() - 1);
	      }
	      write_infile_string(out_infile_result_sub, temp_listname, temp_ss.str());

	    }
	  }
	  else {
	    stringstream temp_ss;
	    temp_ss << "start" << phasespace_optimization[i_m] << "/runlist." << xlist[i_l].xcontribution[i_c].infix_order_contribution << ".dat";
	    string temp_listname = "directory_list";
	    if (phasespace_optimization[i_m] != ""){
	      temp_listname = temp_listname + " " + phasespace_optimization[i_m].substr(1, phasespace_optimization[i_m].size() - 1);
	    }
	    write_infile_string(out_infile_result_sub, temp_listname, temp_ss.str());
	  }
	  
	  out_infile_result_sub << endl;
	}
                 
	out_infile_result_sub.close();
      }
     
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}







void summary_generic::collect_order_result_TSV(){
  Logger logger("summary_generic::collect_order_result_TSV");
  logger << LOG_DEBUG << "called" << endl;

  system_execute(logger, "mkdir " + final_resultdirectory);

  for (int i_o = 0; i_o < yorder.size(); i_o++){
    /*
    for (int i_l = 0; i_l < yorder[i_o].contribution_file.size(); i_l++){
      int x_l = mapping_contribution_file[yorder[i_o].contribution_file[i_l]];
      yorder[i_o].xlist.push_back(xlist[x_l]); 
    }
    */
    yorder[i_o].collect_result_TSV(subgroup);
  }
  /*
  for (int i_o = 0; i_o < yorder.size(); i_o++){
    yorder[i_o].output_result_TSV(final_resultdirectory, subgroup, oset, name_scale_variation_TSV);
    yorder[i_o].output_result_overview_TSV(final_resultdirectory, subgroup, oset, name_scale_variation_TSV);
    yorder[i_o].output_result_plot_TSV(final_resultdirectory, subgroup, oset, name_scale_variation_TSV);
    yorder[i_o].output_result_plot_qTcut_TSV(final_resultdirectory, subgroup, oset, name_scale_variation_TSV);
  }
  */

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_generic::collect_order_result_CV(){
  Logger logger("summary_generic::collect_order_result_CV");
  logger << LOG_DEBUG << "called" << endl;

  for (int i_o = 0; i_o < yorder.size(); i_o++){
    yorder[i_o].collect_result_CV(subgroup);
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_generic::determine_runtime_result_new(){
  Logger logger("summary_generic::determine_runtime_result_new");
  logger << LOG_DEBUG << "called" << endl;

  logger << LOG_DEBUG << "switch_reference = " << oset.switch_reference << endl;

  //  int x_m = 0;

  //  int x_q = 0;
  int x_q = oset.no_qTcut_reference_TSV;
  int x_s = oset.no_reference_TSV;
  int x_r = oset.no_scale_ren_reference_TSV;
  int x_f = oset.no_scale_fact_reference_TSV;

  logger << LOG_DEBUG << "default output from TSV (in summary routine):" << endl;
  logger << LOG_DEBUG << "name_reference_TSV          = " << oset.name_reference_TSV << endl;
  logger << LOG_DEBUG << "no_reference_TSV            = " << oset.no_reference_TSV << endl;
  logger << LOG_DEBUG << "no_scale_ren_reference_TSV  = " << oset.no_scale_ren_reference_TSV << endl;
  logger << LOG_DEBUG << "no_scale_fact_reference_TSV = " << oset.no_scale_fact_reference_TSV << endl;
  logger << LOG_DEBUG << "no_qTcut_reference_TSV      = " << oset.no_qTcut_reference_TSV << endl;

  // runtime extrapolation determination:

  for (int i_o = 0; i_o < yorder.size(); i_o++){
    yorder[i_o].error2_time = 0.;
    for (int i_l = 0; i_l < yorder[i_o].contribution_file.size(); i_l++){
      int x_l = mapping_contribution_file[yorder[i_o].contribution_file[i_l]];
      logger << LOG_DEBUG << "order_contribution_file[" << i_o << "][" << i_l << "] = " << yorder[i_o].contribution_file[i_l] << endl;
      logger << LOG_DEBUG << "x_l = mapping_contribution_file[order_contribution_file[" << i_o << "][" << i_l << "]] = " << mapping_contribution_file[yorder[i_o].contribution_file[i_l]] << endl;
      logger << LOG_DEBUG << "xlist[" << x_l << "].xcontribution[0].xsubprocess[0].error2_time = " << xlist[x_l].xcontribution[0].xsubprocess[0].error2_time << endl;
      yorder[i_o].error2_time += xlist[x_l].xcontribution[0].xsubprocess[0].error2_time;
    }
    logger << LOG_DEBUG << "yorder[" << i_o << "].error2_time = " << yorder[i_o].error2_time << endl;
  }

  logger << LOG_DEBUG << "oset.switch_reference = " << oset.switch_reference << endl;

  vector<double> order_cross_section(yorder.size());
  vector<double> order_relative_error(yorder.size());
  for (int i_o = 0; i_o < order_cross_section.size(); i_o++){
    logger << LOG_DEBUG << "yorder[" << i_o << "].accuracy_no_normalization = " << yorder[i_o].accuracy_no_normalization << endl;

    logger << LOG_DEBUG << "oset.switch_reference = " << oset.switch_reference << endl;
    if (oset.switch_reference == "TSV"){
      logger << LOG_DEBUG << "TSV" << endl;
      logger << LOG_DEBUG << "yorder[yorder[" << i_o << "].accuracy_no_normalization].result_TSV[0][0][" << x_s << "][" << x_r << "][" << x_f << "] = " << yorder[yorder[i_o].accuracy_no_normalization].result_TSV[0][0][x_s][x_r][x_f] << endl;
      order_cross_section[i_o] = yorder[yorder[i_o].accuracy_no_normalization].result_TSV[0][0][x_s][x_r][x_f];
    }
    else if (oset.switch_reference == "CV"){
      logger << LOG_DEBUG << "CV" << endl;
      logger << LOG_DEBUG << "yorder[yorder[" << i_o << "].accuracy_no_normalization].result_CV[0][0][" << osi_no_central_scale_CV << "] = " << yorder[yorder[i_o].accuracy_no_normalization].result_CV[0][0][osi_no_central_scale_CV] << endl;

      order_cross_section[i_o] = yorder[yorder[i_o].accuracy_no_normalization].result_CV[0][0][osi_no_central_scale_CV];
      ///      order_cross_section[i_o] = yorder[yorder[i_o].accuracy_no_normalization].result_CV[0][0][0][osi_no_central_scale_CV];
    }
    else {logger << LOG_FATAL << "Should not happen! No reference scale selected." << endl; exit(1);}

    /*
    if (osi_switch_TSV && oset.name_reference_TSV != ""){
      logger << LOG_DEBUG << "TSV" << endl;
      order_cross_section[i_o] = yorder[yorder[i_o].accuracy_no_normalization].result_TSV[0][0][x_s][x_r][x_f];
    }
    else {
      logger << LOG_DEBUG << "CV" << endl;
      order_cross_section[i_o] = yorder[yorder[i_o].accuracy_no_normalization].result_CV[0][0][0][osi_no_central_scale_CV];
    }
    */
    
    logger << LOG_INFO << "order_cross_section[" << i_o << "] = " << order_cross_section[i_o] << endl;
    order_relative_error[i_o] = yorder[i_o].accuracy_relative;
  }

  for (int i_l = 0; i_l < xlist.size(); i_l++){
    logger << LOG_INFO << "xlist[" << i_l << "]: " << xlist[i_l].xcontribution[0].infix_order_contribution << endl;
    vector<int> contribution_in_order;
    for (int i_o = 0; i_o < yorder.size(); i_o++){
      for (int j_l = 0; j_l < yorder[i_o].contribution_file.size(); j_l++){
	if (i_l == mapping_contribution_file[yorder[i_o].contribution_file[j_l]]){
	  contribution_in_order.push_back(i_o);
	  break;
	}
      }
    }
    for (int i_o = 0; i_o < contribution_in_order.size(); i_o++){
      logger << LOG_DEBUG << "contribution_in_order[" << i_o << "] = " << contribution_in_order[i_o] << endl;
    }
    int x_o = 0;

    double min_order_selection = 1.e99;
    for (int i_o = 0; i_o < contribution_in_order.size(); i_o++){
      double temp_order_selection = order_cross_section[contribution_in_order[i_o]] * order_relative_error[contribution_in_order[i_o]];
      logger << LOG_DEBUG << "order_selection[contribution_in_order[" << i_o << "] = " << contribution_in_order[i_o] << "] = " << setprecision(8) << setw(16) << temp_order_selection << "   order_cross_section = " << setprecision(8) << setw(16) << order_cross_section[contribution_in_order[i_o]] << "   order_relative_error = " << setprecision(8) << setw(16) << order_relative_error[contribution_in_order[i_o]] << endl;
      if (temp_order_selection < min_order_selection){
	min_order_selection = temp_order_selection;
	x_o = contribution_in_order[i_o];
      }
    }
    logger << LOG_DEBUG << "min_order_selection = " << min_order_selection << " @ x_o = " << x_o << endl;

    // reset to zero for xcontribution[0].xsubprocess[0] component !!!
    xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_runtime = 0.;
    
    xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_deviation = 0.;
    xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_normalization_deviation = 0.;
    for (int i_c = 1; i_c < xlist[i_l].xcontribution.size(); i_c++){
      xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_deviation = 0.;
      xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_normalization_deviation = 0.;

      logger << LOG_DEBUG << "xlist[" << i_l << "].xcontribution[" << i_c << "].type_contribution.type_correction = " << xlist[i_l].xcontribution[i_c].type_contribution << "." << xlist[i_l].xcontribution[i_c].type_correction << endl;

      for (int i_p = 1; i_p < xlist[i_l].xcontribution[i_c].xsubprocess.size(); i_p++){
	/*
	cout << "i_c = " << i_c << "   i_p = " << i_p << endl;
	cout << "xlist[" << i_l << "].extrapolated_runtime.size() = " << xlist[i_l].extrapolated_runtime.size() << endl;
	cout << "xlist[" << i_l << "].extrapolated_runtime[" << i_c << "].size() = " << xlist[i_l].extrapolated_runtime[i_c].size() << endl;
	cout << "xlist[" << i_l << "].error2_time.size() = " << xlist[i_l].error2_time.size() << endl;
	cout << "xlist[" << i_l << "].error2_time[" << i_c << "].size() = " << xlist[i_l].error2_time[i_c].size() << endl;
	cout << "xlist[" << i_l << "].extrapolated_n_event.size() = " << xlist[i_l].extrapolated_n_event.size() << endl;
	cout << "xlist[" << i_l << "].extrapolated_n_event[" << i_c << "].size() = " << xlist[i_l].extrapolated_n_event[i_c].size() << endl;
	cout << "xlist[" << i_l << "].extrapolated_sigma_deviation.size() = " << xlist[i_l].extrapolated_sigma_deviation.size() << endl;
	cout << "xlist[" << i_l << "].extrapolated_sigma_deviation[" << i_c << "].size() = " << xlist[i_l].extrapolated_sigma_deviation[i_c].size() << endl;
	*/
	int switch_no_contribution = 0;


	if (oset.switch_reference == "TSV"){
	  logger << LOG_INFO << "xlist[" << i_l << "].xcontribution[" << i_c << "].xsubprocess[" << i_p << "].result_TSV[0][0] = " << xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_TSV[0][x_q][x_s][x_r][x_f] << endl;
	  //	  logger << LOG_INFO << "xlist[" << i_l << "].xcontribution[" << i_c << "].xsubprocess[" << i_p << "].result_CV[0][0] = " << xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_TSV[0][x_q][x_s][x_r][x_f] << endl;
	  if (xlist[i_l].xcontribution[i_c].active_qTcut && no_qTcut_runtime_estimate){
	    x_q = no_qTcut_runtime_estimate;
	  }
	  if (xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_TSV[0][x_q][x_s][x_r][x_f] == 0.){switch_no_contribution = 1;}
	  logger << LOG_DEBUG << "xlist[" << i_l << "].xcontribution[" << i_c << "].xsubprocess[" << i_p << "].result_TSV[0][" << x_q << "][" << x_s << "][" << x_r << "][" << x_f << "] = " << xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_TSV[0][x_q][x_s][x_r][x_f] << endl;
	}
	else if (oset.switch_reference == "CV"){
	  // check why zero !!!
	  if (xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_CV[0][0] == 0.){switch_no_contribution = 1;}
	}
	/*
	if (osi_switch_CV){if (xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_CV[0][0] == 0.){switch_no_contribution = 1;}}
	if (osi_switch_TSV){
	  if (xlist[i_l].xcontribution[i_c].active_qTcut && no_qTcut_runtime_estimate){
	    x_q = no_qTcut_runtime_estimate;
	  }
	  if (xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_TSV[0][x_q][x_s][x_r][x_f] == 0.){switch_no_contribution = 1;}
	  logger << LOG_DEBUG << "xlist[" << i_l << "].xcontribution[" << i_c << "].xsubprocess[" << i_p << "].result_TSV[0][" << x_q << "][" << x_s << "][" << x_r << "][" << x_f << "] = " << xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_TSV[0][x_q][x_s][x_r][x_f] << endl;
	}
	*/

	//	if (osi_switch_TSV){if (xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_TSV[0][0][0][0][0] == 0.){switch_no_contribution = 1;}}
	
	if (osi_switch_CV){logger << LOG_DEBUG << "xlist[" << i_l << "].xcontribution[" << i_c << "].xsubprocess[" << i_p << "].result_CV[0][0] = " << xlist[i_l].xcontribution[i_c].xsubprocess[i_p].result_CV[0][0] << "   switch_no_contribution = " << switch_no_contribution << endl;}
	
	logger << LOG_DEBUG << "xlist[" << i_l << "].xcontribution[" << i_c << "].xsubprocess[" << i_p << "]:   switch_no_contribution = " << switch_no_contribution << endl;
	  
	if (switch_no_contribution){
	  xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_runtime = 0.;
	  xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_n_event = 0.;
	  xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_deviation = 0.;
	  xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization = 0.;
	  xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization_deviation = 0.;


	}
	else {
	  xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_runtime = yorder[x_o].error2_time * xlist[i_l].xcontribution[i_c].xsubprocess[i_p].error2_time / pow(order_relative_error[x_o] * order_cross_section[x_o], 2);

	  logger << LOG_DEBUG << "yorder[" << x_o << "].error2_time = " << yorder[x_o].error2_time << endl;
	  logger << LOG_DEBUG << "order_cross_section[" << x_o << "] = " << order_cross_section[x_o] << endl;
	  logger << LOG_DEBUG << "order_relative_error[" << x_o << "] = " << order_relative_error[x_o] << endl;

	  logger << LOG_DEBUG << "xlist[" << i_l << "].xcontribution[" << i_c << "].xsubprocess[" << i_p << "].error2_time = " << xlist[i_l].xcontribution[i_c].xsubprocess[i_p].error2_time << endl;

	  logger << LOG_DEBUG << "xlist[" << i_l << "].xcontribution[" << i_c << "].xsubprocess[" << i_p << "].extrapolated_runtime = " << xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_runtime << endl;

	  xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_runtime += xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_runtime;
	  
	  xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_n_event = xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_runtime / xlist[i_l].xcontribution[i_c].xsubprocess[i_p].time_per_event;
	  
	  xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_deviation = xlist[i_l].xcontribution[i_c].xsubprocess[i_p].error2_time / sqrt(xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_runtime);
	  xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_deviation += pow(xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_deviation, 2);


	  xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization = order_cross_section[x_o];
	  xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization_deviation = xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_deviation / order_cross_section[x_o];
	  xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_normalization_deviation += pow(xlist[i_l].xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization_deviation, 2);
	}
	// 	cout << "i_l = " << i_l << "   i_c = " << i_c << "   i_p = " << i_p << endl;
      }
      xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_runtime += xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_runtime;
      xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_normalization = order_cross_section[x_o];
      xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_normalization = order_cross_section[x_o];
      xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_deviation += xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_deviation;
      xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_deviation = sqrt(xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_deviation);
      xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_normalization_deviation += xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_normalization_deviation;
      xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_normalization_deviation = sqrt(xlist[i_l].xcontribution[i_c].xsubprocess[0].extrapolated_sigma_normalization_deviation);
      /*
      cout << "xlist[" << i_l << "].extrapolated_runtime.size() = " << xlist[i_l].extrapolated_runtime.size() << endl;
      cout << "xlist[" << i_l << "].extrapolated_runtime[0].size() = " << xlist[i_l].extrapolated_runtime[0].size() << endl;
      cout << "xlist[" << i_l << "].extrapolated_runtime[" << i_c << "].size() = " << xlist[i_l].extrapolated_runtime[i_c].size() << endl;
      */
      //      xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_normalization_deviation = sqrt(xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_normalization_deviation);
      /*
      cout << "xlist[" << i_l << "].extrapolated_n_event.size() = " << xlist[i_l].extrapolated_n_event.size() << endl;
      cout << "xlist[" << i_l << "].extrapolated_n_event[0].size() = " << xlist[i_l].extrapolated_n_event[0].size() << endl;
      */
      xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_n_event = 0.;
    }
    xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_normalization = order_cross_section[x_o];
    xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_deviation = sqrt(xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_deviation);
    xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_normalization_deviation = sqrt(xlist[i_l].xcontribution[0].xsubprocess[0].extrapolated_sigma_normalization_deviation);

    /*
    list_used_runtime_CV[i_l][0].resize(1, 0.);
    for (int i_c = 1; i_c < xlist[i_l].error2_time.size(); i_c++){
      //     list_used_runtime_CV[i_l][i_c].resize(xlist[i_l].error2_time[i_c].size(), 0.);
      for (int i_p = 1; i_p < list_used_runtime_CV[i_l][i_c].size(); i_p++){
	list_used_runtime_CV[i_l][i_c][0] += list_used_runtime_CV[i_l][i_c][i_p];
      }
      list_used_runtime_CV[i_l][0][0] += list_used_runtime_CV[i_l][i_c][0];
    }
    */
  }

  logger << LOG_DEBUG << "finished" << endl;
}




void summary_generic::collect_order_distribution_TSV(){
  Logger logger("summary_generic::collect_order_distribution_TSV");
  logger << LOG_DEBUG << "called" << endl;

  for (int i_o = 0; i_o < yorder.size(); i_o++){
    yorder[i_o].collect_distribution_TSV();
    //    yorder[i_o].collect_distribution_TSV(final_resultdirectory, subgroup, oset.extended_distribution, oset.fakeasymfactor, oset, scalename_TSV);
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}






void summary_generic::collect_order_distribution_CV(){
  Logger logger("summary_generic::collect_order_distribution_CV");
  logger << LOG_DEBUG << "called" << endl;

  for (int i_o = 0; i_o < yorder.size(); i_o++){
    yorder[i_o].collect_distribution_CV();
  }

  /*
  for (int i_o = 0; i_o < yorder.size(); i_o++){
    yorder[i_o].distribution_result_CV.resize(subgroup.size(), vector<vector<vector<double> > > (oset.extended_distribution.size()));
    yorder[i_o].distribution_deviation_CV.resize(subgroup.size(), vector<vector<vector<double> > > (oset.extended_distribution.size()));
    for (int i_g = 0; i_g < subgroup.size(); i_g++){
      for (int i_d = 0; i_d < oset.extended_distribution.size(); i_d++){
	if (!switch_output_distribution[i_d]){continue;}
	yorder[i_o].distribution_result_CV[i_g][i_d].resize(oset.extended_distribution[i_d].n_bins, vector<double> (osi_n_scales_CV, 0.));
	yorder[i_o].distribution_deviation_CV[i_g][i_d].resize(oset.extended_distribution[i_d].n_bins, vector<double> (osi_n_scales_CV, 0.));
      }
    }

    for (int i_g = 0; i_g < subgroup.size(); i_g++){
      for (int i_d = 0; i_d < oset.extended_distribution.size(); i_d++){
	if (!switch_output_distribution[i_d]){continue;}
	for (int i_b = 0; i_b < oset.extended_distribution[i_d].n_bins; i_b++){
	  for (int i_s = 0; i_s < osi_n_scales_CV; i_s++){
	    double dev = 0.;
	    for (int i_l = 0; i_l < yorder[i_o].contribution_file.size(); i_l++){
	      int x_l = mapping_contribution_file[yorder[i_o].contribution_file[i_l]];
	      yorder[i_o].distribution_result_CV[i_g][i_d][i_b][i_s] += xlist[x_l].xcontribution[0].distribution_result_CV[i_g][i_d][i_b][i_s];
	      dev += pow(xlist[x_l].xcontribution[0].distribution_deviation_CV[i_g][i_d][i_b][i_s], 2);
	    }
	    yorder[i_o].distribution_deviation_CV[i_g][i_d][i_b][i_s] = sqrt(dev);
	    stringstream temp_res;
	    temp_res << setw(23) << setprecision(15) << yorder[i_o].distribution_result_CV[i_g][i_d][i_b][i_s];
	    stringstream temp_dev;
	    temp_dev << setw(23) << setprecision(15) << yorder[i_o].distribution_deviation_CV[i_g][i_d][i_b][i_s];
	    logger << LOG_DEBUG_VERBOSE << "XXX   result_CV[" << i_o << "][" << i_g << "][" << i_d << "][" << i_b << "][" << i_s << "] = " << temp_res.str() << " +- " << temp_dev.str() << endl;
	  }
	}
      }
    }
  */
  for (int i_o = 0; i_o < 0; i_o++){
    exit(1);
    // to switch off this loop...
    //  for (int i_o = 0; i_o < yorder.size(); i_o++){
  
    if (switch_output_overview > 0){yorder[i_o].output_distribution_overview_CV();}

    for (int i_g = 0; i_g < subgroup.size(); i_g++){
      for (int i_d = 0; i_d < oset.extended_distribution.size(); i_d++){
	if (!switch_output_distribution[i_d]){continue;}
	for (int i_s = 0; i_s < osi_n_scales_CV; i_s++){
	  if ((oset.extended_distribution[i_d].xdistribution_name).substr(0, 4) == "asym"){
	  /*
	  double result_fw, result_bw, deviation_fw, deviation_bw;
	  ofstream out_asymfile;
	  string filename_asym = final_resultdirectory + "/CV/" + osi_directory_name_scale_CV[i_s] + "/" + oset.extended_distribution[i_d].xdistribution_name + "." + infix_order_contribution[0] + ".dat";
	  logger << LOG_DEBUG_VERBOSE << "filename_asym = " << filename_asym << endl;
	  out_asymfile.open(filename_asym.c_str(), ofstream::out | ofstream::trunc);  
	  if (oset.extended_distribution[i_d].xdistribution_type == "asym1" || oset.extended_distribution[i_d].xdistribution_type == "asym2"){
	    ofstream out_plotfile;
	    string filename_plot = final_resultdirectory + "/CV/" + osi_directory_name_scale_CV[i_s] + "/plot." + oset.extended_distribution[i_d].xdistribution_name + "." + infix_order_contributiony[0] + ".dat";
	    logger << LOG_DEBUG_VERBOSE << "filename_plot = " << filename_plot << endl;
	    out_plotfile.open(filename_plot.c_str(), ofstream::out | ofstream::trunc);  
	    // too slow // checksum << setw(16) << "interval" << " " << setw(16) << "backward cs." << setw(16) << "forward cs." << setw(16) << "cross section" << endl;
	    for (int i_b = 0; i_b < n_bin_distribution[i_d]; i_b++){
	      if (oset.extended_distribution[i_d].xdistribution_type == "asym1"){
		result_fw = yorder[i_o].distribution_result_CV[i_g][i_d][n_bin_distribution[i_d] + i_b][i_s] * oset.fakeasymfactor[i_d];
		result_bw = yorder[i_o].distribution_result_CV[i_g][i_d][i_b][i_s] * oset.fakeasymfactor[i_d];
		deviation_fw = yorder[i_o].distribution_deviation_CV[i_g][i_d][n_bin_distribution[i_d] + i_b][i_s] * oset.fakeasymfactor[i_d];
		deviation_bw = yorder[i_o].distribution_deviation_CV[i_g][i_d][i_b][i_s] * oset.fakeasymfactor[i_d];
	      }
	      else if (oset.extended_distribution[i_d].xdistribution_type == "asym2"){
		result_fw = yorder[i_o].distribution_result_CV[i_g][i_d][2 * i_b + 1][i_s] * oset.fakeasymfactor[i_d];
		result_bw = yorder[i_o].distribution_result_CV[i_g][i_d][2 * i_b + 0][i_s] * oset.fakeasymfactor[i_d];
		deviation_fw = yorder[i_o].distribution_deviation_CV[i_g][i_d][2 * i_b + 1][i_s] * oset.fakeasymfactor[i_d];
		deviation_bw = yorder[i_o].distribution_deviation_CV[i_g][i_d][2 * i_b + 0][i_s] * oset.fakeasymfactor[i_d];
	      }
	      out_asymfile << setw(25) << "backw. cross section: " << "   [" << setw(5) << setprecision(4) << bin_edge[i_d][i_b] << ";" << setw(5) << setprecision(4) << bin_edge[i_d][i_b + 1] << "]   " << setw(16) << setprecision(8) << result_bw * oset.extended_distribution[i_d].step << setw(16) << setprecision(8) << deviation_bw * oset.extended_distribution[i_d].step << endl;
	      out_asymfile << setw(25) << "forw. cross section: " << "   [" << setw(5) << setprecision(4) << bin_edge[i_d][i_b] << ";" << setw(5) << setprecision(4) << bin_edge[i_d][i_b + 1] << "]   " << setw(16) << setprecision(8) << result_fw * oset.extended_distribution[i_d].step << setw(16) << setprecision(8) << deviation_fw * oset.extended_distribution[i_d].step << endl;
	      out_asymfile << setw(25) << "forw.-backw. asymmetry: " << "   [" << setw(5) << setprecision(4) << bin_edge[i_d][i_b] << ";" << setw(5) << setprecision(4) << bin_edge[i_d][i_b + 1] << "]   " << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	      out_plotfile << setw(5) << setprecision(4) << bin_edge[i_d][i_b] << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	      double sum_result_fw = 0., sum_result_bw = 0.;
	      for (int c = i_b; c < n_bin_distribution[i_d]; c++){
		if (oset.extended_distribution[i_d].xdistribution_type == "asym1"){
		  sum_result_fw += yorder[i_o].distribution_result_CV[i_g][i_d][n_bin_distribution[i_d] + c][i_s] * oset.fakeasymfactor[i_d];
		  sum_result_bw += yorder[i_o].distribution_result_CV[i_g][i_d][c][i_s] * oset.fakeasymfactor[i_d];
		}
		else if (oset.extended_distribution[i_d].xdistribution_type == "asym2"){
		  sum_result_fw += yorder[i_o].distribution_result_CV[i_g][i_d][2 * c + 1][i_s] * oset.fakeasymfactor[i_d];
		  sum_result_bw += yorder[i_o].distribution_result_CV[i_g][i_d][2 * c + 0][i_s] * oset.fakeasymfactor[i_d];
		}
	      }
	    }
	    out_plotfile << setw(5) << setprecision(4) << oset.extended_distribution[i_d].end << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	    out_asymfile.close();
	    out_plotfile.close();
	  }
	  else {
	    result_fw = yorder[i_o].distribution_result_CV[i_g][i_d][1][i_s] * oset.fakeasymfactor[i_d];
	    result_bw = yorder[i_o].distribution_result_CV[i_g][i_d][0][i_s] * oset.fakeasymfactor[i_d];
	    deviation_fw = yorder[i_o].distribution_deviation_CV[i_g][i_d][1][i_s] * oset.fakeasymfactor[i_d];
	    deviation_bw = yorder[i_o].distribution_deviation_CV[i_g][i_d][0][i_s] * oset.fakeasymfactor[i_d];
	    out_asymfile << setw(25) << "backw. cross section: " << setw(16) << setprecision(8) << result_bw << setw(16) << setprecision(8) << deviation_bw << endl;
	    out_asymfile << setw(25) << "forw. cross section: " << setw(16) << setprecision(8) << result_fw << setw(16) << setprecision(8) << deviation_fw << endl;
	    out_asymfile << setw(25) << "forw.-backw. asymmetry: " << setw(16) << setprecision(8) << (result_fw - result_bw) / (result_fw + result_bw) << setw(16) << 2. * (deviation_fw * result_bw + deviation_bw * result_fw) / pow(result_bw + result_fw, 2) << endl;
	    out_asymfile.close();
	  }
	  */
	  }
	  else {
	    ofstream out_plotfile;
	    string filename_plot = final_resultdirectory + "/CV/" + osi_directory_name_scale_CV[i_s] + "/plot." + oset.extended_distribution[i_d].xdistribution_name + ".." + yorder[i_o].resultdirectory + ".dat";
	    logger << LOG_DEBUG_VERBOSE << "filename_plot = " << filename_plot << endl;
	    out_plotfile.open(filename_plot.c_str(), ofstream::out | ofstream::trunc);  
	    for (int i_b = 0; i_b < yorder[i_o].distribution_result_CV[i_g][i_d].size(); i_b++){

	      out_plotfile << setw(10) << oset.extended_distribution[i_d].bin_edge[i_b] << setw(16) << setprecision(8) << osi_unit_factor_distribution * yorder[i_o].distribution_result_CV[i_g][i_d][i_b][i_s] / oset.extended_distribution[i_d].bin_width[i_b] << setw(16) << setprecision(8) << osi_unit_factor_distribution * yorder[i_o].distribution_deviation_CV[i_g][i_d][i_b][i_s] / oset.extended_distribution[i_d].bin_width[i_b] << endl;
	      // norm ->	    out_plotfile << setw(10) << oset.extended_distribution[i_d].bin_edge[i_b] << setw(16) << setprecision(8) << osi_unit_factor_distribution * yorder[i_o].distribution_result_CV[i_g][i_d][i_b][i_s] << setw(16) << setprecision(8) << osi_unit_factor_distribution * yorder[i_o].distribution_deviation_CV[i_g][i_d][i_b][i_s] << endl;
	      double sum_result = 0.;
	      for (int j_b = i_b; j_b < oset.extended_distribution[i_d].n_bins; j_b++){sum_result += yorder[i_o].distribution_result_CV[i_g][i_d][j_b][i_s];}
	      //	      for (int j_b = i_b; j_b < oset.extended_distribution[i_d].n_bins; j_b++){sum_result += yorder[i_o].distribution_result_CV[i_g][i_d][j_b][i_s] / oset.extended_distribution[i_d].bin_width[i_b];}
	    }
	    out_plotfile << setw(10) << oset.extended_distribution[i_d].bin_edge[oset.extended_distribution[i_d].n_bins] << setw(16) << setprecision(8) << osi_unit_factor_distribution * yorder[i_o].distribution_result_CV[i_g][i_d][yorder[i_o].distribution_result_CV[i_g][i_d].size() - 1][i_s] / oset.extended_distribution[i_d].bin_width[oset.extended_distribution[i_d].n_bins - 1] << setw(16) << setprecision(8) << osi_unit_factor_distribution * yorder[i_o].distribution_deviation_CV[i_g][i_d][yorder[i_o].distribution_deviation_CV[i_g][i_d].size() - 1][i_s] / oset.extended_distribution[i_d].bin_width[oset.extended_distribution[i_d].n_bins - 1] << endl; // endpoint
	    // norm ->	    out_plotfile << setw(10) << oset.extended_distribution[i_d].bin_edge[oset.extended_distribution[i_d].n_bins] << setw(16) << setprecision(8) << osi_unit_factor_distribution * yorder[i_o].distribution_result_CV[i_g][i_d][yorder[i_o].distribution_result_CV[i_g][i_d].size() - 1][i_s] << setw(16) << setprecision(8) << osi_unit_factor_distribution * yorder[i_o].distribution_deviation_CV[i_g][i_d][yorder[i_o].distribution_deviation_CV[i_g][i_d].size() - 1][i_s] << endl; // endpoint
	    
	    out_plotfile.close();
	  }
	}
      }
    }



  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_generic::output_filename_TSV(string & filename_begin, string & filename_end){
  Logger logger("summary_generic::output_filename_TSV");
  logger << LOG_DEBUG << "called" << endl;

  filename_complete = filename_begin + "/complete/" + filename_end;
  filename_ren = filename_begin + "/ren/" + filename_end;
  filename_fact = filename_begin + "/fact/" + filename_end;
  filename_equal = filename_begin + "/equal/" + filename_end;
  filename_antipodal = filename_begin + "/antipodal/" + filename_end;

  logger << LOG_DEBUG << "finished" << endl;
}

