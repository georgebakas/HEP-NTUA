#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

////////////////////
//  constructors  //
////////////////////
summary_list::summary_list(){}

summary_list::summary_list(string & infilename, summary_generic & _generic){
  Logger logger("summary_list::summary_list");
  logger << LOG_DEBUG << "called" << endl;

  // only one 'summary_generic'
  ygeneric = &_generic;
  // only one 'observable_set'
  osi = &ygeneric->oset;

  xcontribution.resize(1, summary_contribution("", *this));

  logger << LOG_DEBUG << "ygeneric->oset.value_qTcut_distribution.size() = " << ygeneric->oset.value_qTcut_distribution.size() << endl;
  logger << LOG_DEBUG << "osi->value_qTcut_distribution.size() = " << osi->value_qTcut_distribution.size() << endl;

  in_contribution_order_alpha_s = 0;
  in_contribution_order_alpha_e = 0;
  photon_induced = 0;

  //  string vs;
  string s0;
  //  vector<string> vs0;
  //  vector<int> vi0;

  string filename;
  string rundirectory;
  vector<string> readin;

  char LineBuffer[128];
  //  vector<string> contribution_name;
  //  vector<vector<string> > contribution_directory;
  //  vector<int> contribution_variationmode;

  //  int QCD_order = 0;
  logger << LOG_DEBUG << "infilename = " << infilename << endl;
  ifstream in_file(infilename.c_str());
  readin.clear();
  while (in_file.getline(LineBuffer, 1024)){readin.push_back(LineBuffer);}
  logger << LOG_DEBUG << "readin.size() = " << readin.size() << endl;
  vector<string> user_variable;
  vector<string> user_value;
  vector<string> user_variable_additional;
  get_userinput_from_readin(user_variable, user_variable_additional, user_value, readin);
     
  int counter_contribution = 1;
  for (int i = 0; i < user_variable.size(); i++){if (user_variable[i] == "type_contribution"){counter_contribution++;}
    //    logger << LOG_DEBUG << "i = " << i << "   counter_contribution = " << counter_contribution << endl;
  }
  xcontribution.resize(counter_contribution);
  counter_contribution = 0;

  for (int i = 0; i < user_variable.size(); i++){
    if (user_variable[i] == "processname"){
      processname = user_value[i];
      xcontribution[counter_contribution].processname = user_value[i];
    }
    else if (user_variable[i] == "resultdirectory"){
      resultdirectory = user_value[i];
      xcontribution[counter_contribution].resultdirectory = user_value[i];
    }
    else if (user_variable[i] == "type_perturbative_order"){
      type_perturbative_order = user_value[i];
      xcontribution[counter_contribution].type_perturbative_order = user_value[i];
    }
    else if (user_variable[i] == "subtraction_method"){
      type_subtraction_method = user_value[i];
      xcontribution[counter_contribution].type_subtraction_method = user_value[i];
    }
    else if (user_variable[i] == "contribution_order_alpha_s"){
      in_contribution_order_alpha_s = atoi(user_value[i].c_str());
      xcontribution[counter_contribution].in_contribution_order_alpha_s = atoi(user_value[i].c_str());
    }
    else if (user_variable[i] == "contribution_order_alpha_e"){
      in_contribution_order_alpha_e = atoi(user_value[i].c_str());
      xcontribution[counter_contribution].in_contribution_order_alpha_e = atoi(user_value[i].c_str());
    }
    else if (user_variable[i] == "type_contribution"){
      counter_contribution++;
      xcontribution[counter_contribution] = summary_contribution(user_value[i], *this);
      //      logger << LOG_DEBUG << "xcontribution[" << counter_contribution << "]->ylist.resultdirectory = " << xcontribution[counter_contribution].ylist->resultdirectory << endl;
    }
    else if (user_variable[i] == "type_correction"){
      xcontribution[counter_contribution].type_correction = user_value[i];
    }
    else if (user_variable[i] == "interference"){
      xcontribution[counter_contribution].interference = atoi(user_value[i].c_str());
    }
    else if (user_variable[i] == "photon_induced"){
      photon_induced = atoi(user_value[i].c_str());
      xcontribution[counter_contribution].photon_induced = atoi(user_value[i].c_str());
    }
    else if (user_variable[i] == "directory" || user_variable[i] == "directory_list"){
      int x_n = -1;
      for (int i_n = 0; i_n < xcontribution[counter_contribution].extended_directory_name.size(); i_n++){
	if (user_variable_additional[i] == xcontribution[counter_contribution].extended_directory_name[i_n]){x_n = i_n; break;}
      }
      if (x_n == -1){
	xcontribution[counter_contribution].extended_directory_name.push_back(user_variable_additional[i]);
	xcontribution[counter_contribution].extended_directory.push_back(vector<string> ());
	x_n = xcontribution[counter_contribution].extended_directory_name.size() - 1;
      }
      if (user_variable[i] == "directory"){
	xcontribution[counter_contribution].extended_directory[x_n].push_back(user_value[i]);

	// directory and directory_extra should be later removed !!!
	xcontribution[counter_contribution].directory.push_back(user_value[i]);
	xcontribution[counter_contribution].directory_extra.push_back(0);
	
      }
      else if (user_variable[i] == "directory_list"){
	string infilepath = "../" + user_value[i];
	logger << LOG_INFO << "directory_list = " << user_value[i] << "   path = " << infilepath << endl;
	ifstream in_file(infilepath.c_str());
	vector<string> readin_list(0);
	while (in_file.getline(LineBuffer, 1024)){readin_list.push_back(LineBuffer);}
	for (int i_l = 0; i_l < readin_list.size(); i_l++){
	  logger << LOG_INFO << "directory_list: readin_list[" << setw(4) << right << i_l << "] = " << readin_list[i_l] << endl;
	  xcontribution[counter_contribution].extended_directory[x_n].push_back(readin_list[i_l]);

	  // directory and directory_extra should be later removed !!!
	  xcontribution[counter_contribution].directory.push_back(readin_list[i_l]);
	  xcontribution[counter_contribution].directory_extra.push_back(0);

	}	  
      }

      
    }
    // should be removed later...
    else if (user_variable[i] == "directory_extra"){
      xcontribution[counter_contribution].directory.push_back(user_value[i]);
      xcontribution[counter_contribution].directory_extra.push_back(1);
    }




    
  }
  logger << LOG_INFO << "resultdirectory = " << resultdirectory << endl;
  for (int i_c = 0; i_c < xcontribution.size(); i_c++){
    for (int i_n = 0; i_n < xcontribution[i_c].extended_directory_name.size(); i_n++){
      logger << LOG_INFO << "xcontribution[" << i_c << "].infix_order_contribution = " << xcontribution[i_c].infix_order_contribution << "   extended_directory_name[" << i_n << "] = " << xcontribution[i_c].extended_directory_name[i_n] << "   entries: " << xcontribution[i_c].extended_directory[i_n].size() << endl;
    }
    logger << LOG_INFO << "xcontribution[" << i_c << "].infix_order_contribution = " << xcontribution[i_c].infix_order_contribution << "     directory.size() = " << xcontribution[i_c].directory.size() << endl;
  }



  
  for (int i_c = 0; i_c < xcontribution.size(); i_c++){logger << LOG_DEBUG << "BEFORE   &xcontribution[" << i_c << "] = " << &xcontribution[i_c] << endl;}
  //  for (int i_x = 0; i_x < type_contribution.size(); i_x++){logger << LOG_DEBUG << "type_contribution[" << i_x << "] = " << type_contribution[i_x] << endl;}
  //  for (int i_x = 0; i_x < type_correction.size(); i_x++){logger << LOG_DEBUG << "type_correction[" << i_x << "] = " << type_correction[i_x] << endl;}
  //  for (int i_x = 0; i_x < type_subtraction_method.size(); i_x++){logger << LOG_DEBUG << "type_subtraction_method[" << i_x << "] = " << type_subtraction_method[i_x] << endl;}
  //  for (int i_x = 0; i_x < interference.size(); i_x++){logger << LOG_DEBUG << "interference[" << i_x << "] = " << interference[i_x] << endl;}
  //  for (int i_x = 0; i_x < type_contribution.size(); i_x++){logger << LOG_DEBUG << "type_contribution[" << i_x << "] = " << type_contribution[i_x] << endl;}
  //  for (int i_x = 0; i_x < type_contribution.size(); i_x++){logger << LOG_DEBUG << "type_contribution[" << i_x << "] = " << type_contribution[i_x] << endl;}

  stringstream temp_photon_induced;
  if (xcontribution[0].photon_induced){temp_photon_induced << "a";}
  stringstream temp_order;
  stringstream temp_path;
  temp_order << type_perturbative_order;
  temp_path << type_perturbative_order;
  if (type_subtraction_method != "---"){
    temp_order << "." << type_subtraction_method;
    temp_path << "." << type_subtraction_method;
  }
  temp_order << "." << temp_photon_induced.str() << in_contribution_order_alpha_s << in_contribution_order_alpha_e;
  temp_path << "/" << temp_photon_induced.str() << in_contribution_order_alpha_s << in_contribution_order_alpha_e;
  xcontribution[0].infix_order_contribution = temp_order.str();
  xcontribution[0].infix_path_contribution = temp_path.str();
  logger << LOG_DEBUG << "xcontribution[0].infix_order_contribution   = " << xcontribution[0].infix_order_contribution << endl;

  active_qTcut = 0;
  output_n_qTcut = 1;
  selection_n_qTcut = 1;
  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    if (xcontribution[i_c].type_contribution == "" ||
	xcontribution[i_c].type_contribution == "CT" ||
	xcontribution[i_c].type_contribution == "RT" ||
	xcontribution[i_c].type_contribution == "CT2" ||
	xcontribution[i_c].type_contribution == "RVA" ||
	xcontribution[i_c].type_contribution == "RCA" ||
	xcontribution[i_c].type_contribution == "RRA" ||
	xcontribution[i_c].type_contribution == "L2RT" ||
	xcontribution[i_c].type_contribution == "L2CT" ||
	xcontribution[i_c].type_contribution == "CJ" ||
	xcontribution[i_c].type_contribution == "RJ" ||
	xcontribution[i_c].type_contribution == "CJ2" ||
	xcontribution[i_c].type_contribution == "RVJ" ||
	xcontribution[i_c].type_contribution == "RCJ" ||
	xcontribution[i_c].type_contribution == "RRJ" ||
	xcontribution[i_c].type_contribution == "L2RJ" ||
	xcontribution[i_c].type_contribution == "L2CJ"){
      xcontribution[i_c].active_qTcut = 1;
      xcontribution[i_c].output_n_qTcut = osi->n_qTcut;
      xcontribution[i_c].selection_n_qTcut = osi->value_qTcut_distribution.size();
      active_qTcut = 1;
      output_n_qTcut = osi->n_qTcut;
      selection_n_qTcut = osi->value_qTcut_distribution.size();
    }
    else {
      xcontribution[i_c].active_qTcut = 0;
      xcontribution[i_c].output_n_qTcut = 1;
      xcontribution[i_c].selection_n_qTcut = 1;
    }
  }
  xcontribution[0].active_qTcut = active_qTcut;
  xcontribution[0].output_n_qTcut = output_n_qTcut;
  xcontribution[0].selection_n_qTcut = selection_n_qTcut;

  logger << LOG_INFO << "list      active_qTcut = " << active_qTcut << "   output_n_qTcut = " << setw(3) << output_n_qTcut << "   selection_n_qTcut = " << selection_n_qTcut << endl;
  for (int i_c = 0; i_c < xcontribution.size(); i_c++){
    logger << LOG_INFO << "i_c = " << i_c << "   active_qTcut = " << xcontribution[i_c].active_qTcut << "   output_n_qTcut = " << setw(3) << xcontribution[i_c].output_n_qTcut << "   selection_n_qTcut = " << xcontribution[i_c].selection_n_qTcut << endl;
  }

  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    stringstream temp_contribution;
    logger << LOG_DEBUG << "type_contribution[" << i_c << "] = " << xcontribution[i_c].type_contribution << endl;
    temp_contribution << xcontribution[i_c].type_contribution;
    logger << LOG_DEBUG << "type_correction[" << i_c << "] = " << xcontribution[i_c].type_correction << endl;
    if (xcontribution[i_c].type_correction != "---" && xcontribution[i_c].type_correction != ""){temp_contribution << "." << xcontribution[i_c].type_correction;}
    logger << LOG_DEBUG << "interference[" << i_c << "] = " << xcontribution[i_c].interference << endl;
    if (xcontribution[i_c].interference != 0){temp_contribution << ".i" << xcontribution[i_c].interference;}
    xcontribution[i_c].infix_order_contribution = temp_order.str() + "." + temp_contribution.str();
    xcontribution[i_c].infix_path_contribution = temp_path.str() + "/" + temp_contribution.str();
    xcontribution[i_c].infix_contribution = temp_contribution.str();
    logger << LOG_DEBUG << "infix_order_contribution[" << i_c << "] = " << xcontribution[i_c].infix_order_contribution << endl;
    logger << LOG_DEBUG << "infix_path_contribution [" << i_c << "] = " << xcontribution[i_c].infix_path_contribution << endl;
    logger << LOG_DEBUG << "infix_contribution      [" << i_c << "] = " << xcontribution[i_c].infix_contribution << endl;
    logger << LOG_DEBUG << "active_qTcut[" << i_c << "] = " << xcontribution[i_c].active_qTcut << endl;
    logger << LOG_DEBUG << "output_n_qTcut[" << i_c << "] = " << xcontribution[i_c].output_n_qTcut << endl;
    logger << LOG_DEBUG << "selection_n_qTcut[" << i_c << "] = " << xcontribution[i_c].selection_n_qTcut << endl;
    for (int i_q = 0; i_q < osi->value_qTcut_distribution.size(); i_q++){
      logger << LOG_DEBUG << "osi->value_qTcut_distribution[" << setw(2) << i_q << "] = " << osi->value_qTcut_distribution[i_q] << endl;
    }
  }
  
  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    xcontribution[i_c].average_factor = ygeneric->average_factor;
    if (xcontribution[i_c].type_contribution == "born"){ygeneric->generic->list_subprocess_born(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "VA" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_V_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "VA" && xcontribution[i_c].type_correction == "QEW"){ygeneric->generic->list_subprocess_V_QEW(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "VA" && xcontribution[i_c].type_correction == "MIX"){ygeneric->generic->list_subprocess_V_MIX(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "CA" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_C_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "CA" && xcontribution[i_c].type_correction == "QEW"){ygeneric->generic->list_subprocess_C_QEW(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "RA" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_R_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "RA" && xcontribution[i_c].type_correction == "QEW"){ygeneric->generic->list_subprocess_R_QEW(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "RA" && xcontribution[i_c].type_correction == "MIX"){ygeneric->generic->list_subprocess_R_MIX(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    
    else if (xcontribution[i_c].type_contribution == "VT" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_V_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "CT" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_C_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "RT" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_R_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    
    else if (xcontribution[i_c].type_contribution == "VJ" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_V_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "CJ" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_C_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "RJ" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_R_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    
    else if (xcontribution[i_c].type_contribution == "loop"){ygeneric->generic->list_subprocess_born(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    
    else if (xcontribution[i_c].type_contribution == "VT2" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_V2_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "CT2" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_C2_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    //    else if (xcontribution[i_c].type_contribution == "VT2" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_V2_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s - 1, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    //    else if (xcontribution[i_c].type_contribution == "CT2" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_C2_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s - 1, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "RVA" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_RV_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "RCA" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_RC_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "RRA" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_RR_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    
    else if (xcontribution[i_c].type_contribution == "VJ2" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_V2_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "CJ2" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_C2_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "RVJ" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_RV_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "RCJ" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_RC_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "RRJ" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_RR_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    
    else if (xcontribution[i_c].type_contribution == "L2I"){ygeneric->generic->list_subprocess_born(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    
    else if (xcontribution[i_c].type_contribution == "L2VA" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_V_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "L2CA" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_C_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "L2RA" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_R_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    
    else if (xcontribution[i_c].type_contribution == "L2VT" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_V_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "L2CT" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_C_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "L2RT" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_R_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    
    else if (xcontribution[i_c].type_contribution == "L2VJ" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_V_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "L2CJ" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_C_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    else if (xcontribution[i_c].type_contribution == "L2RJ" && xcontribution[i_c].type_correction == "QCD"){ygeneric->generic->list_subprocess_R_QCD(xcontribution[i_c].subprocess, xcontribution[i_c].subgroup_no_member, in_contribution_order_alpha_s, in_contribution_order_alpha_e, xcontribution[i_c].interference);}
    
    else {logger << LOG_ERROR << "ERROR:   type_contribution[" << i_c << "] = " << xcontribution[i_c].type_contribution <<  "type_correction[" << i_c << "] = " << xcontribution[i_c].type_correction << endl;}
    
    xcontribution[i_c].xsubprocess.resize(xcontribution[i_c].subprocess.size());
    for (int i_p = 1; i_p < xcontribution[i_c].xsubprocess.size(); i_p++){
      logger << LOG_INFO << "xcontribution[" << i_c << "].subprocess[" << i_p << "] = " << xcontribution[i_c].subprocess[i_p] << endl;
      xcontribution[i_c].xsubprocess[i_p] = summary_subprocess(xcontribution[i_c].subprocess[i_p], xcontribution[i_c]);
    }
  }
  
  // subprocess not yet filled here !!! to be checked !!!

  logger << LOG_DEBUG << "finished" << endl;
}

///////////////////////
//  access elements  //
///////////////////////

///////////////
//  methods  //
///////////////

void summary_list::output_info(){
  Logger logger("summary_list::output_info");
  logger << LOG_DEBUG << "called" << endl;
  
  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    logger << LOG_INFO << resultdirectory << "   xcontribution[" << i_c << "].infix_order_contribution[" << i_c << "] = " << xcontribution[i_c].infix_order_contribution << "   with   " << xcontribution[i_c].subprocess.size() - 1 << "   subprocesses" << endl;
    for (int i_p = 1; i_p < xcontribution[i_c].subprocess.size(); i_p++){
      logger << LOG_INFO << resultdirectory << "   xcontribution[" << i_c << "].subprocess[" << i_p << "] = " << xcontribution[i_c].subprocess[i_p] << endl;
    }
    logger.newLine(LOG_INFO);
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::calculate_runtime(){
  Logger logger("summary_list::calculate_runtime");
  logger << LOG_DEBUG << "called" << endl;

  logger << LOG_INFO << "xcontribution[0].xsubprocess.size() = " << xcontribution[0].xsubprocess.size() << endl;
  logger << LOG_INFO << "xcontribution[0].xsubprocess[0].error2_time = " << xcontribution[0].xsubprocess[0].error2_time << "   --- should be 0 here !!!" << endl;

  // reset to zero for xcontribution[0].xsubprocess[0] component !!!
  xcontribution[0].xsubprocess[0].used_runtime = 0.;
  xcontribution[0].xsubprocess[0].used_n_event = 0;
  
  xcontribution[0].xsubprocess[0].error2_time = 0.;
  
  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    xcontribution[0].xsubprocess[0].error2_time += xcontribution[i_c].xsubprocess[0].error2_time;
    logger << LOG_INFO << "SUMRUNTIME   contribution   " << setw(15) << xcontribution[0].xsubprocess[0].name << "   error2_time = " << xcontribution[0].xsubprocess[0].error2_time << endl;

    logger << LOG_INFO << "xcontribution[" << i_c << "].infix_order_contribution = " << xcontribution[i_c].infix_order_contribution << endl;
    for (int i_p = 0; i_p < xcontribution[i_c].subprocess.size(); i_p++){
      logger << LOG_INFO << "convergence[" << i_c << "][" << i_p << "] = " << setw(20) << xcontribution[i_c].subprocess[i_p] << "   list_error2_time = " << setw(20) << xcontribution[i_c].xsubprocess[i_p].error2_time << "   list_time_per_event = " <<  xcontribution[i_c].xsubprocess[i_p].time_per_event << endl;
    }
  }
  
  for (int i_c = 0; i_c < xcontribution.size(); i_c++){
    logger << LOG_INFO << "contribution convergence[" << i_c << "][0]: " << "   xcontribution[" << i_c << "].xsubprocess[0].error2_time = " << setw(20) << xcontribution[i_c].xsubprocess[0].error2_time << endl;
  }
  
  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    for (int i_p = 0; i_p < xcontribution[i_c].xsubprocess.size(); i_p++){
      logger << LOG_DEBUG << setw(20) << xcontribution[i_c].infix_order_contribution << "[" << i_c << "][" << setw(2) << i_p << "] = " << setw(20) << xcontribution[i_c].subprocess[i_p] << "   used_runtime = " << setprecision(3) << setw(10) << xcontribution[i_c].xsubprocess[i_p].used_runtime / 3600 << " h" << endl; 
    }
    logger << LOG_DEBUG << endl;
  }

  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    xcontribution[0].xsubprocess[0].used_runtime += xcontribution[i_c].xsubprocess[0].used_runtime;
    xcontribution[0].xsubprocess[0].used_n_event += xcontribution[i_c].xsubprocess[0].used_n_event;
  }

  logger << LOG_DEBUG << "      used_runtime = " << setprecision(3) << setw(10) << xcontribution[0].xsubprocess[0].used_runtime / 3600 << " h" << endl; 
  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    logger << LOG_DEBUG << "[" << i_c << "]   used_runtime = " << setprecision(3) << setw(10) << xcontribution[i_c].xsubprocess[0].used_runtime / 3600 << " h" << endl; 
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}





void summary_list::output_readin_extrapolated_runtime(ofstream & outfile_readin_extrapolated){
  Logger logger("summary_list::output_readin_extrapolated_runtime");
  logger << LOG_DEBUG << "called" << endl;

  /////////////////////////////////////////////////////////////
  //  extrapolated runtime - readin.extrapolated.runtime.dat //
  /////////////////////////////////////////////////////////////

  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    for (int i_p = 1; i_p < xcontribution[i_c].xsubprocess.size(); i_p++){
      outfile_readin_extrapolated << left << setw(25) << xcontribution[i_c].infix_path_contribution;
      outfile_readin_extrapolated << left << setw(20) << xcontribution[i_c].xsubprocess[i_p].name;
      outfile_readin_extrapolated << right << setw(12) << (long long)(xcontribution[i_c].xsubprocess[i_p].extrapolated_runtime) << "   ";
      outfile_readin_extrapolated << right << setw(12) << (long long)(xcontribution[i_c].xsubprocess[i_p].extrapolated_n_event) << "   ";
      outfile_readin_extrapolated << right << setw(15) << setprecision(8) << xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization_deviation;
      outfile_readin_extrapolated << right << setw(15) << setprecision(8) << xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization;
      outfile_readin_extrapolated << right << setw(15) << setprecision(8) << xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_deviation;
      outfile_readin_extrapolated << endl;
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::output_extrapolated_runtime(ofstream & outfile_extrapolated){
  Logger logger("summary_list::output_extrapolated_runtime");
  logger << LOG_DEBUG << "called" << endl;

  ///////////////////////////////////////////////////////////
  //  extrapolated runtime - info.extrapolated.runtime.txt //
  ///////////////////////////////////////////////////////////

  //  outfile_extrapolated << left << setw(25) << "# folder" << setw(20) << "subprocess" << setw(16) << right << "runtime" << "   " << setw(12) << "n_event" << "   " << left << setw(15) << "deviation" << setw(15) << "norm" << setw(15) << "norm_err" << endl;
  //  outfile_extrapolated << endl;


  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    for (int i_p = 0; i_p < xcontribution[i_c].xsubprocess.size(); i_p++){
      /*
      stringstream temp_ext_runtime;
      long long temp_h = (long long)(xcontribution[i_c].xsubprocess[i_p].extrapolated_runtime / 3600);
      long long temp_min = (long long)(xcontribution[i_c].xsubprocess[i_p].extrapolated_runtime / 60) - temp_h * 60;
      long long temp_sec = (long long)(xcontribution[i_c].xsubprocess[i_p].extrapolated_runtime) - temp_h * 3600 - temp_min * 60;
      temp_ext_runtime << setw(4) << temp_h << " h " << setw(2) << temp_min << " m " << setw(2) << temp_sec << " s";
      */
      outfile_extrapolated << left << setw(25) << xcontribution[i_c].infix_path_contribution;
      outfile_extrapolated << left << setw(20) << xcontribution[i_c].xsubprocess[i_p].name;
      outfile_extrapolated << time_hms_from_double(xcontribution[i_c].xsubprocess[i_p].extrapolated_runtime) << "   ";
      //	outfile_extrapolated << temp_ext_runtime.str() << "   ";
      if (i_p == 0){
	/*
	outfile_extrapolated << left << setw(25) << xcontribution[i_c].infix_path_contribution;
	outfile_extrapolated << left << setw(20) << xcontribution[i_c].xsubprocess[i_p].name;
	outfile_extrapolated << temp_ext_runtime.str() << "   ";
	*/
	outfile_extrapolated << endl;
	for (int i_s = 0; i_s < 129; i_s++){outfile_extrapolated << "-";} outfile_extrapolated << endl;
      }
      else {
	outfile_extrapolated << right << setw(12) << (long long)(xcontribution[i_c].xsubprocess[i_p].extrapolated_n_event) << "   ";
	outfile_extrapolated << right << setw(15) << setprecision(8) << xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization_deviation;
	outfile_extrapolated << right << setw(15) << setprecision(8) << xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization;
	outfile_extrapolated << right << setw(15) << setprecision(8) << xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_deviation;
	outfile_extrapolated << endl;
	
      //" << setprecision(3) << setw(10) << error2_time[i_c][i_p] << endl; 
	//	outfile_extrapolated << left << setw(20) << xcontribution[i_c].subprocess[i_p] << "   ext. runtime = " << temp_ext_runtime.str() << "   ext. n_event = " << setw(12) << (long long)(xcontribution[i_c].xsubprocess[i_p].extrapolated_n_event) << "   ext. abs. deviation = " << setprecision(8) << setw(15) << xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_deviation <<"   ext. rel. deviation = " << xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization_deviation <<  "   sigma_normalization = " << xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization << endl; // << "   constant = " << setprecision(3) << setw(10) << error2_time[i_c][i_p] << endl; 
      }
      //      outfile_extrapolated << endl;
      //extrapolated_runtime[i_c][i_p] / time_per_event[i_c][i_p]
      //error2_time[i_c][i_p] / sqrt(extrapolated_runtime[i_c][i_p])
    }
    outfile_extrapolated << endl;
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::output_used_runtime(ofstream & outfile_used){
  Logger logger("summary_list::output_used_runtime");
  logger << LOG_DEBUG << "called" << endl;

  ///////////////////////////////////////////
  //  used runtime - info.used.runtime.txt //
  ///////////////////////////////////////////

  for (int i_c = 1; i_c < xcontribution.size(); i_c++){
    for (int i_p = 0; i_p < xcontribution[i_c].xsubprocess.size(); i_p++){
      /*
      stringstream temp_used_runtime;
      long long temp_h = (long long)(xcontribution[i_c].xsubprocess[i_p].used_runtime / 3600);
      long long temp_min = (long long)(xcontribution[i_c].xsubprocess[i_p].used_runtime / 60) - temp_h * 60;
      long long temp_sec = (long long)(xcontribution[i_c].xsubprocess[i_p].used_runtime) - temp_h * 3600 - temp_min * 60;
      temp_used_runtime << setw(4) << temp_h << " h " << setw(2) << temp_min << " m " << setw(2) << temp_sec << " s";
      */
      outfile_used << left << setw(25) << xcontribution[i_c].infix_path_contribution;
      outfile_used << left << setw(20) << xcontribution[i_c].xsubprocess[i_p].name;
      //	outfile_used << temp_used_runtime.str() << "   ";
      outfile_used << time_hms_from_double(xcontribution[i_c].xsubprocess[i_p].used_runtime) << "   ";
      if (i_p == 0){
	/*
	outfile_used << left << setw(25) << xcontribution[i_c].infix_path_contribution;
	outfile_used << left << setw(20) << xcontribution[i_c].xsubprocess[i_p].name;
	//	outfile_used << temp_used_runtime.str() << "   ";
	outfile_used << time_hms_from_double(xcontribution[i_c].xsubprocess[i_p].used_runtime) << "   ";
	*/
	outfile_used << endl;
	for (int i_s = 0; i_s < 129; i_s++){outfile_used << "-";} outfile_used << endl;
      }
      else {
	outfile_used << right << setw(12) << (long long)(xcontribution[i_c].xsubprocess[i_p].used_n_event) << "   ";
	outfile_used << right << setw(15) << setprecision(8) << "";//xcontribution[i_c].xsubprocess[i_p].used_sigma_normalization_deviation;
	outfile_used << right << setw(15) << setprecision(8) << xcontribution[i_c].xsubprocess[i_p].extrapolated_sigma_normalization;
	outfile_used << right << setw(15) << setprecision(8) << "";//xcontribution[i_c].xsubprocess[i_p].used_sigma_deviation;
	outfile_used << right << setw(15) << setprecision(8) << xcontribution[i_c].xsubprocess[i_p].error2_time;
	outfile_used << endl;
      }
    }
    outfile_used << endl;
  }

  logger << LOG_DEBUG << "finished" << endl;
}




















// not yet implemented:
void summary_list::extrapolation_linear(vector<double> & xvalue, double & min_extrapolation, double & max_extrapolation, vector<double> & data_result, vector<double> & data_deviation2, vector<double> & X, vector<double> & dX, double & chi2_dof){
  Logger logger("summary_list::extrapolation_linear");
  logger << LOG_DEBUG << "called" << endl;

  vector<double> A(3, 0.);
  vector<double> B(2, 0.);
  vector<vector<double> > Ainv(3, vector<double> (3, 0.));
  double detA;
  vector<vector<double> > H(3, vector<double> (3, 0.));
  vector<vector<double> > Hinv(3, vector<double> (3, 0.));
  double detH;
  
  //  double a, b, c;
  //  double da, db, dc;
 
  //  vector<double> data_deviation2(xvalue.size());

  //  double min_extrapolation = 0.00;
  //  double min_extrapolation = 0.05;
  //  double max_extrapolation = 0.50;

  for (int i_q = 0; i_q < xvalue.size(); i_q++){
    if (xvalue[i_q] < min_extrapolation){continue;}
    if (xvalue[i_q] > max_extrapolation){break;}
    if (data_deviation2[i_q] != 0.){
      for (int i = 0; i < 3; i++){A[i] += pow(xvalue[i_q], i) / data_deviation2[i_q];}
      for (int i = 0; i < 2; i++){B[i] += pow(xvalue[i_q], i) * data_result[i_q] / data_deviation2[i_q];}
    }
  }
  
  Ainv[1][1] = A[0];
  Ainv[1][2] = -A[1];
  Ainv[2][2] = A[2];
  
  detA = A[0] * A[2] - A[1] * A[1];
  /*
  Ainv[1][1] = A[0] * A[2] - A[1] * A[1];
  Ainv[1][2] = A[1] * A[2] - A[0] * A[3];
  Ainv[1][3] = A[1] * A[3] - A[2] * A[2];
  
  Ainv[2][2] = A[0] * A[4] - A[2] * A[2];
  Ainv[2][3] = A[2] * A[3] - A[1] * A[4];
  
  Ainv[3][3] = A[2] * A[4] - A[3] * A[3];
  
  detA = A[0] * A[2] * A[4] + 2 * A[1] * A[2] * A[3] - A[0] * A[3] * A[3] - A[4] * A[1] * A[1] - A[2] * A[2] * A[2];
  */
  
  H[1][1] = 2 * A[2];
  H[1][2] = 2 * A[1];
  H[2][2] = 2 * A[0];
  
  detH = H[1][1] * H[2][2] - H[1][2] * H[1][2];
  /*
  H[1][1] = 2 * A[4];
  H[1][2] = 2 * A[3];
  H[1][3] = 2 * A[2];
  H[2][2] = 2 * A[2];
  H[2][3] = 2 * A[1];
  H[3][3] = 2 * A[0];

  detH = H[1][1] * H[2][2] * H[3][3] + 2 * H[1][2] * H[2][3] * H[1][3] - H[2][2] * H[1][3] * H[1][3] - H[3][3] * H[1][2] * H[1][2] - H[1][1] * H[2][3] * H[2][3];
  */  

  Hinv[1][1] = H[2][2] / detH;
  Hinv[2][2] = H[1][1] / detH;
  Hinv[1][2] = -H[1][2] / detH;

  /*
  Hinv[1][1] = (H[2][2] * H[3][3] - H[2][3] * H[2][3]) / detH;
  Hinv[1][2] = (H[1][3] * H[2][3] - H[1][2] * H[3][3]) / detH;
  Hinv[1][3] = (H[1][2] * H[2][3] - H[1][3] * H[2][2]) / detH;
  
  Hinv[2][2] = (H[1][1] * H[3][3] - H[1][3] * H[1][3]) / detH;
  Hinv[2][3] = (H[1][3] * H[1][2] - H[1][1] * H[2][3]) / detH;
  
  Hinv[3][3] = (H[1][1] * H[2][2] - H[1][2] * H[1][2]) / detH;
  */
  
  if (detA != 0.){
    
    X[1] = (Ainv[1][1] * B[1] + Ainv[1][2] * B[0]) / detA;
    X[0] = (Ainv[1][2] * B[1] + Ainv[2][2] * B[0]) / detA;
    /*
    X[2] = (Ainv[1][1] * B[2] + Ainv[1][2] * B[1] + Ainv[1][3] * B[0]) / detA;
    X[1] = (Ainv[1][2] * B[2] + Ainv[2][2] * B[1] + Ainv[2][3] * B[0]) / detA;
    X[0] = (Ainv[1][3] * B[2] + Ainv[2][3] * B[1] + Ainv[3][3] * B[0]) / detA;
    */
    
    dX[0] = sqrt(Hinv[1][1]);
    dX[1] = sqrt(Hinv[2][2]);
    /*
    dX[2] = Hinv[1][1];
    dX[1] = Hinv[2][2];
    dX[0] = Hinv[3][3];
    */
  }

  int n_data_point = 0;
  double chi2 = 0.;
  // chi2 determination:
  for (int i_q = 0; i_q < xvalue.size(); i_q++){
    if (xvalue[i_q] < min_extrapolation){continue;}
    if (xvalue[i_q] > max_extrapolation){break;}
    n_data_point++;
    //    chi2 += abs((data_result[i_q] - (X[0] + xvalue[i_q] * X[1])) / sqrt(data_deviation2[i_q]));
    chi2 += pow(data_result[i_q] - (X[0] + xvalue[i_q] * X[1]), 2) / data_deviation2[i_q];
  }
  chi2_dof = chi2 / (n_data_point - 2);



  /*
    else {
    //	start values:
    a = 0.;
    b = 0.;
    c = 0.;
    da = 0.;
    db = 0.;
    dc = 0.;
    }
  */

  logger << LOG_DEBUG << "finished" << endl;
}


void summary_list::extrapolation_quadratic(vector<double> & xvalue, double & min_extrapolation, double & max_extrapolation, vector<double> & data_result, vector<double> & data_deviation2, vector<double> & X, vector<double> & dX, double & chi2_dof){
  Logger logger("extrapolation_quadratic");
  logger << LOG_DEBUG << "called" << endl;

  vector<double> A(5, 0.);
  vector<double> B(3, 0.);
  vector<vector<double> > Ainv(4, vector<double> (4, 0.));
  double detA;
  vector<vector<double> > H(4, vector<double> (4, 0.));
  vector<vector<double> > Hinv(4, vector<double> (4, 0.));
  double detH;
  
  //  double a, b, c;
  //  double da, db, dc;
 
  //  vector<double> data_deviation2(xvalue.size());

  //  double min_extrapolation = 0.00;
  //  double min_extrapolation = 0.05;
  //  double max_extrapolation = 0.50;

  for (int i_q = 0; i_q < xvalue.size(); i_q++){
    if (xvalue[i_q] < min_extrapolation){continue;}
    if (xvalue[i_q] > max_extrapolation){break;}
    if (data_deviation2[i_q] != 0.){
      for (int i = 0; i < 5; i++){A[i] += pow(xvalue[i_q], i) / data_deviation2[i_q];}
      for (int i = 0; i < 3; i++){B[i] += pow(xvalue[i_q], i) * data_result[i_q] / data_deviation2[i_q];}
    }
  }
  
  Ainv[1][1] = A[0] * A[2] - A[1] * A[1];
  Ainv[1][2] = A[1] * A[2] - A[0] * A[3];
  Ainv[1][3] = A[1] * A[3] - A[2] * A[2];
  
  Ainv[2][2] = A[0] * A[4] - A[2] * A[2];
  Ainv[2][3] = A[2] * A[3] - A[1] * A[4];
  
  Ainv[3][3] = A[2] * A[4] - A[3] * A[3];
  
  detA = A[0] * A[2] * A[4] + 2 * A[1] * A[2] * A[3] - A[0] * A[3] * A[3] - A[4] * A[1] * A[1] - A[2] * A[2] * A[2];
  
  
  H[1][1] = 2 * A[4];
  H[1][2] = 2 * A[3];
  H[1][3] = 2 * A[2];
  H[2][2] = 2 * A[2];
  H[2][3] = 2 * A[1];
  H[3][3] = 2 * A[0];
  
  detH = H[1][1] * H[2][2] * H[3][3] + 2 * H[1][2] * H[2][3] * H[1][3] - H[2][2] * H[1][3] * H[1][3] - H[3][3] * H[1][2] * H[1][2] - H[1][1] * H[2][3] * H[2][3];
  
  Hinv[1][1] = (H[2][2] * H[3][3] - H[2][3] * H[2][3]) / detH;
  Hinv[1][2] = (H[1][3] * H[2][3] - H[1][2] * H[3][3]) / detH;
  Hinv[1][3] = (H[1][2] * H[2][3] - H[1][3] * H[2][2]) / detH;
  
  Hinv[2][2] = (H[1][1] * H[3][3] - H[1][3] * H[1][3]) / detH;
  Hinv[2][3] = (H[1][3] * H[1][2] - H[1][1] * H[2][3]) / detH;
  
  Hinv[3][3] = (H[1][1] * H[2][2] - H[1][2] * H[1][2]) / detH;
  
  
  if (detA != 0.){
    X[2] = (Ainv[1][1] * B[2] + Ainv[1][2] * B[1] + Ainv[1][3] * B[0]) / detA;
    X[1] = (Ainv[1][2] * B[2] + Ainv[2][2] * B[1] + Ainv[2][3] * B[0]) / detA;
    X[0] = (Ainv[1][3] * B[2] + Ainv[2][3] * B[1] + Ainv[3][3] * B[0]) / detA;
    
    dX[2] = Hinv[1][1];
    dX[1] = Hinv[2][2];
    dX[0] = Hinv[3][3];
  }

  int n_data_point = 0;
  double chi2 = 0.;
  // chi2 determination:
  for (int i_q = 0; i_q < xvalue.size(); i_q++){
    if (xvalue[i_q] < min_extrapolation){continue;}
    if (xvalue[i_q] > max_extrapolation){break;}
    n_data_point++;
    chi2 += pow(data_result[i_q] - (X[0] + xvalue[i_q] * X[1] + pow(xvalue[i_q], 2) * X[2]), 2) / data_deviation2[i_q];
    //    chi2 += abs((data_result[i_q] - (X[0] + xvalue[i_q] * X[1] + pow(xvalue[i_q], 2) * X[2])) / sqrt(data_deviation2[i_q]));
  }
  chi2_dof = chi2 / (n_data_point - 3);


  /*
    else {
    //	start values:
    a = 0.;
    b = 0.;
    c = 0.;
    da = 0.;
    db = 0.;
    dc = 0.;
    }
  */

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::extrapolation_TSV(int order_fit){
  Logger logger("summary_list::extrapolation_TSV");
  logger << LOG_DEBUG << "called" << endl;


  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
    if (i_m > 0){continue;} // no moments implemented so far !!!
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){

      if (!(osi->value_qTcut.size() == xcontribution[0].result_TSV[i_m][i_g].size())){cout << "osi->value_qTcut.size() = " << osi->value_qTcut.size() << " != " << xcontribution[0].result_TSV[i_m][i_g].size() << " = xcontribution[0].result_TSV[i_m][i_g].size()" << endl;}
      assert(osi->value_qTcut.size() == xcontribution[0].result_TSV[i_m][i_g].size());
      if (!(osi->value_qTcut.size() == xcontribution[0].deviation_TSV[i_m][i_g].size())){cout << "osi->value_qTcut.size() = " << osi->value_qTcut.size() << " != " << xcontribution[0].deviation_TSV[i_m][i_g].size() << " = xcontribution[0].deviation_TSV[i_m][i_g].size()" << endl;}
      assert(osi->value_qTcut.size() == xcontribution[0].deviation_TSV[i_m][i_g].size());
      
      static double xmin = ygeneric->min_qTcut_extrapolation;
      static double xmax = ygeneric->max_qTcut_extrapolation;
      
      static int x_s = osi->no_reference_TSV;
      static int x_r = osi->no_scale_ren_reference_TSV;
      static int x_f = osi->no_scale_fact_reference_TSV;
      logger << LOG_DEBUG << "x_s = " << x_s << endl;
      logger << LOG_DEBUG << "x_r = " << x_r << endl;
      logger << LOG_DEBUG << "x_f = " << x_f << endl;
      
      vector<double> all_chi2_dof(osi->value_qTcut.size(), 1.e99);
      vector<double> all_extrapolation_result(osi->value_qTcut.size(), 0.);
      vector<double> all_extrapolation_deviation(osi->value_qTcut.size(), 0.);
      vector<double> X(3, 0.);
      vector<double> dX(3, 0.);
      vector<double> data_result(osi->value_qTcut.size());
      vector<double> data_deviation2(osi->value_qTcut.size());
      
      vector<vector<vector<int> > > no_xmax_TSV(osi->n_extended_set_TSV);
      for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
	if (!ygeneric->switch_output_scaleset[i_s]){continue;}
	no_xmax_TSV[i_s].resize(osi->n_scale_ren_TSV[i_s], vector<int> (osi->n_scale_fact_TSV[i_s], 0));
      }
      
      for (int i_s = 0; i_s < result_TSV[i_m][i_g].size(); i_s++){
	for (int i_r = 0; i_r < result_TSV[i_m][i_g][i_s].size(); i_r++){
	  for (int i_f = 0; i_f < result_TSV[i_m][i_g][i_s][i_r].size(); i_f++){
	    for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q++){
	      data_result[i_q] = xcontribution[0].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f];
	      data_deviation2[i_q] = pow(xcontribution[0].deviation_TSV[i_m][i_g][i_q][i_s][i_r][i_f], 2);
	    }
	    int no_xmax = 0;
	    
	    logger << LOG_DEBUG << "Extrapolation at i_s = " << i_s << "   i_r = " << i_r << "   i_f = " << i_f << endl;
	    extrapolation_fit_range_determination(osi->value_qTcut, data_result, data_deviation2, all_extrapolation_result, all_extrapolation_deviation, all_chi2_dof, order_fit, no_xmax, xmin, xmax);
	    logger << LOG_DEBUG << "Extrapolation between  " << xmin << "  and  " << xmax << "  considered optimal." << endl;
	    deviation_extrapolation_TSV[i_m][i_g][i_s][i_r][i_f] = all_extrapolation_deviation[no_xmax];
	    no_xmax_TSV[i_s][i_r][i_f] = no_xmax;
	  }
	}
      }
      
  /*
  for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q++){
    data_result[i_q] = xcontribution[0].result_TSV[i_m][i_g][i_q][x_s][x_r][x_f];
    data_deviation2[i_q] = pow(xcontribution[0].deviation_TSV[i_m][i_g][i_q][x_s][x_r][x_f], 2);
  }
  int no_xmax = 0;

  extrapolation_fit_range_determination(osi->value_qTcut, data_result, data_deviation2, all_extrapolation_result, all_extrapolation_deviation, all_chi2_dof, order_fit, no_xmax, xmin, xmax);
  */

      /*
      xmax = ygeneric->max_qTcut_extrapolation;
      logger << LOG_DEBUG << "Fit performed with the standard value xmax = " << xmax << endl;
      */      
      
      
      for (int i_s = 0; i_s < result_TSV[i_m][i_g].size(); i_s++){
	for (int i_r = 0; i_r < result_TSV[i_m][i_g][i_s].size(); i_r++){
	  for (int i_f = 0; i_f < result_TSV[i_m][i_g][i_s][i_r].size(); i_f++){
	    xmax = osi->value_qTcut[no_xmax_TSV[i_s][i_r][i_f]];
	    vector<double> data_result(osi->value_qTcut.size());
	    vector<double> data_deviation2(osi->value_qTcut.size());
	    for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q++){
	      data_result[i_q] = xcontribution[0].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f];
	      data_deviation2[i_q] = pow(xcontribution[0].deviation_TSV[i_m][i_g][i_q][i_s][i_r][i_f], 2);
	    }
	    double chi2_dof = 0.;
	    vector<double> X(3, 0.);
	    vector<double> dX(3, 0.);
	    if	(order_fit == 1){extrapolation_linear(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	    else if (order_fit == 2){extrapolation_quadratic(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	    else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
	    result_TSV[i_m][i_g][i_s][i_r][i_f] = X[0];
	    //	out_chi2[i_s][i_r][i_f] = chi2_dof;

	    double temp_error_fitvalue_at_fixed_qTcut = 0.;
	    if	(order_fit == 1){
	      //	  y(qTcut) = X[0] + qTcut * X[1]
	      //	  temp_error_fitvalue_at_fixed_qTcut = abs(osi->value_qTcut[no_value_xmin] * X[1]) / 2;
	      temp_error_fitvalue_at_fixed_qTcut = abs(xmin * X[1]) / 2;
	    }
	    else if (order_fit == 2){
	      //	  y(qTcut) = X[0] + qTcut * X[1] + qTcut² * X[2]
	      //	  temp_error_fitvalue_at_fixed_qTcut = abs(osi->value_qTcut[no_value_xmin] * X[1] + pow(osi->value_qTcut[no_value_xmin], 2) * X[2]) / 2;
	      temp_error_fitvalue_at_fixed_qTcut = abs(xmin * X[1] + pow(xmin, 2) * X[2]) / 2;
	    }


	/*
	if (i_s == 0 && i_r == 1 && i_f == 1){
	  logger << LOG_DEBUG << "extrapolation    = " << setw(15) << setprecision(8) << result_TSV[i_m][i_g][i_s][i_r][i_f] << "   chi2_dof = " << chi2_dof << endl;
	  for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q++){
	    logger << LOG_DEBUG << "data_result[" << setw(3) << i_q << "] = " << setw(15) << setprecision(8) << data_result[i_q] << " +- " << setw(15) << setprecision(8) << sqrt(data_deviation2[i_q]) << endl;
	  }
	}
	*/
	
	    double extrapolation_plus = 0.;
	    for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q++){
	      data_result[i_q] = xcontribution[0].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f] + xcontribution[0].deviation_TSV[i_m][i_g][i_q][i_s][i_r][i_f];
	    }
	    X.resize(3, 0.);
	    dX.resize(3, 0.);
	    if	(order_fit == 1){extrapolation_linear(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	    else if (order_fit == 2){extrapolation_quadratic(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	    else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
	    extrapolation_plus = X[0];
	    
	    double extrapolation_minus = 0.;
	    for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q++){
	      data_result[i_q] = xcontribution[0].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f] - xcontribution[0].deviation_TSV[i_m][i_g][i_q][i_s][i_r][i_f];
	    }
	    X.resize(3, 0.);
	    dX.resize(3, 0.);
	    if	(order_fit == 1){extrapolation_linear(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	    else if (order_fit == 2){extrapolation_quadratic(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	    else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
	    extrapolation_minus = X[0];
	    
	    // larger error estimate for large uncetainty due to high lowest qTcut value:
	    
	    if (temp_error_fitvalue_at_fixed_qTcut > deviation_extrapolation_TSV[i_m][i_g][i_s][i_r][i_f]){
	      logger << LOG_DEBUG << resultdirectory << "   Extrapolation error is changed due to fixed lowest value:" << endl;
	      logger << LOG_DEBUG << "     old   deviation_extrapolation_CV[" << i_m << "][" << i_g << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << deviation_extrapolation_TSV[i_m][i_g][i_s][i_r][i_f] << endl;
	      deviation_extrapolation_TSV[i_m][i_g][i_s][i_r][i_f] = temp_error_fitvalue_at_fixed_qTcut;
	      logger << LOG_DEBUG << "     new   deviation_extrapolation_CV[" << i_m << "][" << i_g << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << deviation_extrapolation_TSV[i_m][i_g][i_s][i_r][i_f] << endl;
	    }
	    else {
	      logger << LOG_DEBUG << resultdirectory << "   Extrapolation error is not changed due to fixed lowest value:" << endl;
	      logger << LOG_DEBUG << "     old   deviation_extrapolation_CV[" << i_m << "][" << i_g << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << deviation_extrapolation_TSV[i_m][i_g][i_s][i_r][i_f] << endl;
	      logger << LOG_DEBUG << "    (new)  deviation_extrapolation_CV[" << i_m << "][" << i_g << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << temp_error_fitvalue_at_fixed_qTcut << endl;
	    }

	    deviation_statistics_TSV[i_m][i_g][i_s][i_r][i_f] = (extrapolation_plus - extrapolation_minus) / 2;
	    //  pure statistical error, "extrapolated" to qTcut->0.
	    //  A further error estimate from the extrapolation curve (in particular its actually random fit range) should be added !!!
	    
	    deviation_TSV[i_m][i_g][i_s][i_r][i_f] = sqrt(pow(deviation_statistics_TSV[i_m][i_g][i_s][i_r][i_f], 2) + pow(deviation_extrapolation_TSV[i_m][i_g][i_s][i_r][i_f], 2));
	    logger << LOG_DEBUG << "result_TSV[" << i_m << "][" << i_g << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << endl << setw(23) << setprecision(15) << result_TSV[i_m][i_g][i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << deviation_TSV[i_m][i_g][i_s][i_r][i_f] << "    [ " << setw(23) << setprecision(15) << deviation_statistics_TSV[i_m][i_g][i_s][i_r][i_f] << " (statistics) +- " << deviation_extrapolation_TSV[i_m][i_g][i_s][i_r][i_f] << " (extrapolation)" << " ]" << endl;


	    for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q += 5){
	      logger << LOG_DEBUG << setw(3) << i_q << "   qTcut = " << setw(10) << setprecision(5) << osi->value_qTcut[i_q] << "   result = " << setw(23) << setprecision(15) << xcontribution[0].result_TSV[i_m][i_g][i_q][i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << xcontribution[0].deviation_TSV[i_m][i_g][i_q][i_s][i_r][i_f] << endl;
	    }
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::extrapolation_CV(int order_fit){
  Logger logger("summary_list::extrapolation_CV");
  logger << LOG_DEBUG << "called" << endl;

  /*
  int no_value_xmin = 0;
  for (int i_q = 0; i_q < xvalue.size(); i_q++){
    if (xvalue[i_q] >= xmin){no_value_xmin = i_q; break;}
  }
  logger << LOG_DEBUG << "xmin = " << xmin << " at qTcut value no_value_xmin = " << no_value_xmin << endl;
  */  

  for (int i_m = 0; i_m < osi->n_moments + 1; i_m++){
    if (i_m > 0){continue;} // no moments implemented so far !!!
    for (int i_g = 0; i_g < ygeneric->subgroup.size(); i_g++){

      if (!(osi->value_qTcut.size() == xcontribution[0].result_CV[i_m][i_g].size())){cout << "osi->value_qTcut.size() = " << osi->value_qTcut.size() << " != " << xcontribution[0].result_CV[i_m][i_g].size() << " = xcontribution[0].result_CV[i_m][i_g].size()" << endl;}
      assert(osi->value_qTcut.size() == xcontribution[0].result_CV[i_m][i_g].size());
      if (!(osi->value_qTcut.size() == xcontribution[0].deviation_CV[i_m][i_g].size())){cout << "osi->value_qTcut.size() = " << osi->value_qTcut.size() << " != " << xcontribution[0].deviation_CV[i_m][i_g].size() << " = xcontribution[0].deviation_CV[i_m][i_g].size()" << endl;}
      assert(osi->value_qTcut.size() == xcontribution[0].deviation_CV[i_m][i_g].size());
      
      static double xmin = ygeneric->min_qTcut_extrapolation;
      static double xmax = ygeneric->max_qTcut_extrapolation;
      
      static int x_s = (osi->n_scales_CV - 1 ) / 2;
      logger << LOG_DEBUG << "x_s = " << x_s << endl;
      
      vector<double> all_chi2_dof(osi->value_qTcut.size(), 1.e99);
      vector<double> all_extrapolation_result(osi->value_qTcut.size(), 0.);
      vector<double> all_extrapolation_deviation(osi->value_qTcut.size(), 0.);
      vector<double> X(3, 0.);
      vector<double> dX(3, 0.);
      vector<double> data_result(osi->value_qTcut.size());
      vector<double> data_deviation2(osi->value_qTcut.size());
      
      vector<int> no_xmax_CV(osi->n_scales_CV, 0);

      for (int i_s = 0; i_s < result_CV[i_m][i_g].size(); i_s++){
	for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q++){
	  data_result[i_q] = xcontribution[0].result_CV[i_m][i_g][i_q][i_s];
	  data_deviation2[i_q] = pow(xcontribution[0].deviation_CV[i_m][i_g][i_q][i_s], 2);
	}
	int no_xmax = 0;
	logger << LOG_DEBUG << "Extrapolation at i_s = " << i_s << endl;
	extrapolation_fit_range_determination(osi->value_qTcut, data_result, data_deviation2, all_extrapolation_result, all_extrapolation_deviation, all_chi2_dof, order_fit, no_xmax, xmin, xmax);
	logger << LOG_DEBUG << "Extrapolation between  " << xmin << "  and  " << xmax << "  considered optimal." << endl;
	
	deviation_extrapolation_CV[i_m][i_g][i_s] = all_extrapolation_deviation[no_xmax];
	no_xmax_CV[i_s] = no_xmax;
      }
      
      for (int i_s = 0; i_s < result_CV[i_m][i_g].size(); i_s++){
	xmax = osi->value_qTcut[no_xmax_CV[i_s]];
	vector<double> data_result(osi->value_qTcut.size());
	vector<double> data_deviation2(osi->value_qTcut.size());
	for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q++){
	  data_result[i_q] = xcontribution[0].result_CV[i_m][i_g][i_q][i_s];
	  data_deviation2[i_q] = pow(xcontribution[0].deviation_CV[i_m][i_g][i_q][i_s], 2);
	}
	double chi2_dof = 0.;
	vector <double> X(3, 0.);
	vector <double> dX(3, 0.);
	if	(order_fit == 1){extrapolation_linear(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else if (order_fit == 2){extrapolation_quadratic(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
	result_CV[i_m][i_g][i_s] = X[0];
	//    out_chi2[i_s] = chi2_dof;
	
	double temp_error_fitvalue_at_fixed_qTcut = 0.;
	if	(order_fit == 1){
	  //	  y(qTcut) = X[0] + qTcut * X[1]
	  //	  temp_error_fitvalue_at_fixed_qTcut = abs(osi->value_qTcut[no_value_xmin] * X[1]) / 2;
	  temp_error_fitvalue_at_fixed_qTcut = abs(xmin * X[1]) / 2;
	}
	else if (order_fit == 2){
	  //	  y(qTcut) = X[0] + qTcut * X[1] + qTcut² * X[2]
	  //	  temp_error_fitvalue_at_fixed_qTcut = abs(osi->value_qTcut[no_value_xmin] * X[1] + pow(osi->value_qTcut[no_value_xmin], 2) * X[2]) / 2;
	  temp_error_fitvalue_at_fixed_qTcut = abs(xmin * X[1] + pow(xmin, 2) * X[2]) / 2;
	}

	/*
	  if (i_s == 0){
      logger << LOG_DEBUG << "extrapolation    = " << setw(15) << setprecision(8) << result_CV[i_m][i_g][i_s] << "   chi2_dof = " << chi2_dof << endl;
      for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q++){
	logger << LOG_DEBUG << "data_result[" << setw(3) << i_q << "] = " << setw(15) << setprecision(8) << data_result[i_q] << " +- " << setw(15) << setprecision(8) << sqrt(data_deviation2[i_q]) << endl;
      }
    }
    */

	double extrapolation_plus = 0.;
	for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q++){
	  data_result[i_q] = xcontribution[0].result_CV[i_m][i_g][i_q][i_s] + xcontribution[0].deviation_CV[i_m][i_g][i_q][i_s];
	}
	X.resize(3, 0.);
	dX.resize(3, 0.);
	if	(order_fit == 1){extrapolation_linear(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else if (order_fit == 2){extrapolation_quadratic(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
	extrapolation_plus = X[0];
	
	double extrapolation_minus = 0.;
	for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q++){
	  data_result[i_q] = xcontribution[0].result_CV[i_m][i_g][i_q][i_s] - xcontribution[0].deviation_CV[i_m][i_g][i_q][i_s];
	}
	X.resize(3, 0.);
	dX.resize(3, 0.);
	if	(order_fit == 1){extrapolation_linear(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else if (order_fit == 2){extrapolation_quadratic(osi->value_qTcut, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
	extrapolation_minus = X[0];


	// larger error estimate for large uncetainty due to high lowest qTcut value:

	if (temp_error_fitvalue_at_fixed_qTcut > deviation_extrapolation_CV[i_m][i_g][i_s]){
	  logger << LOG_DEBUG << resultdirectory << "   Extrapolation error is changed due to fixed lowest value:" << endl;
	  logger << LOG_DEBUG << "     old   deviation_extrapolation_CV[" << i_m << "][" << i_g << "][" << i_s << "] = " << deviation_extrapolation_CV[i_m][i_g][i_s] << endl;
	  deviation_extrapolation_CV[i_m][i_g][i_s] = temp_error_fitvalue_at_fixed_qTcut;
	  logger << LOG_DEBUG << "     new   deviation_extrapolation_CV[" << i_m << "][" << i_g << "][" << i_s << "] = " << deviation_extrapolation_CV[i_m][i_g][i_s] << endl;
	}
	else {
	  logger << LOG_DEBUG << resultdirectory << "   Extrapolation error is not changed due to fixed lowest value:" << endl;
	  logger << LOG_DEBUG << "     old   deviation_extrapolation_CV[" << i_m << "][" << i_g << "][" << i_s << "] = " << deviation_extrapolation_CV[i_m][i_g][i_s] << endl;
	  logger << LOG_DEBUG << "    (new)  deviation_extrapolation_CV[" << i_m << "][" << i_g << "][" << i_s << "] = " << temp_error_fitvalue_at_fixed_qTcut << endl;
	}

	deviation_statistics_CV[i_m][i_g][i_s] = (extrapolation_plus - extrapolation_minus) / 2;
	
	deviation_CV[i_m][i_g][i_s] = sqrt(pow(deviation_statistics_CV[i_m][i_g][i_s], 2) + pow(deviation_extrapolation_CV[i_m][i_g][i_s], 2));
	
	logger << LOG_INFO << "result_CV[" << i_m << "][" << i_g << "][" << i_s << "] = " << endl << setw(23) << setprecision(15) << result_CV[i_m][i_g][i_s] << " +- " << setw(23) << setprecision(15) << deviation_CV[i_m][i_g][i_s] << "    [ " << setw(23) << setprecision(15) << deviation_statistics_CV[i_m][i_g][i_s] << " (statistics) +- " << deviation_extrapolation_CV[i_m][i_g][i_s] << " (extrapolation)" << " ]" << endl;
	
	
	for (int i_q = 0; i_q < osi->value_qTcut.size(); i_q += 5){
	  logger << LOG_DEBUG << setw(3) << i_q << "   qTcut = " << setw(10) << setprecision(5) << osi->value_qTcut[i_q] << "   result = " << setw(23) << setprecision(15) << xcontribution[0].result_CV[i_m][i_g][i_q][i_s] << " +- " << setw(23) << setprecision(15) << xcontribution[0].deviation_CV[i_m][i_g][i_q][i_s] << endl;
	}
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::extrapolation_fit_range_determination(vector<double> & xvalue, vector<double> & data_result, vector<double> & data_deviation2, vector<double> & all_extrapolation_result, vector<double> & all_extrapolation_deviation, vector<double> & all_chi2_dof, int order_fit, int & no_xmax, double & xmin, double & xmax){
  Logger logger("extrapolation_fit_range_determination");
  logger << LOG_DEBUG << "called" << endl;

  vector<double> X(3, 0.);
  vector<double> dX(3, 0.);

  //  double min_max_value_extrapolation_range = 0.05;
  //  int min_n_qTcut_extrapolation_range = 10;

  int no_value_xmin = 0;
  for (int i_q = 0; i_q < xvalue.size(); i_q++){
    //    logger << LOG_INFO << "xmin = " << setw(23) << setprecision(15) << xmin << "   xvalue[" << setw(3) << i_q << "] = " << setw(23) << setprecision(15) << xvalue[i_q] << endl;  
    //    if (xvalue[i_q] == xmin){no_value_xmin = i_q; break;}
    if (xvalue[i_q] >= xmin){no_value_xmin = i_q; break;}
  }
  logger << LOG_DEBUG << "xmin = " << xmin << " at qTcut value no_value_xmin = " << no_value_xmin << endl;

  int start_value_xmax = 0;
  for (int i_q = 0; i_q < xvalue.size(); i_q++){if (xvalue[i_q] >= ygeneric->min_max_value_extrapolation_range){start_value_xmax = i_q; break;}}
  //  for (int i_q = 0; i_q < xvalue.size(); i_q++){if (xvalue[i_q] == ygeneric->min_max_value_extrapolation_range){start_value_xmax = i_q; break;}}
  logger << LOG_DEBUG << "start_value_xmax  after ygeneric->min_max_value_extrapolation_range = " << setw(5) << setprecision(3) << ygeneric->min_max_value_extrapolation_range << " : " << setw(3) << start_value_xmax << endl;

  if (start_value_xmax < no_value_xmin + ygeneric->min_n_qTcut_extrapolation_range){start_value_xmax = no_value_xmin + ygeneric->min_n_qTcut_extrapolation_range;}
  logger << LOG_DEBUG << "start_value_xmax  after ygeneric->min_n_qTcut_extrapolation_range   = " << setw(5) << ygeneric->min_n_qTcut_extrapolation_range << " : " << setw(3) << start_value_xmax << endl;

  int end_value_xmax = xvalue.size();
  for (int i_q = 0; i_q < xvalue.size(); i_q++){if (xvalue[i_q] >= ygeneric->max_qTcut_extrapolation){end_value_xmax = i_q + 1; break;}}
  logger << LOG_DEBUG << "end_value_xmax  after ygeneric->max_qTcut_extrapolation             = " << setw(5) << ygeneric->max_qTcut_extrapolation << " : " << setw(3) << end_value_xmax << endl;

  for (int j_q = start_value_xmax; j_q < end_value_xmax; j_q++){
    //  for (int j_q = start_value_xmax; j_q < xvalue.size(); j_q++){
    xmax = xvalue[j_q];
    
    if (order_fit == 1){extrapolation_linear(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, all_chi2_dof[j_q]);}
    else if (order_fit == 2){extrapolation_quadratic(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, all_chi2_dof[j_q]);}
    else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
    all_extrapolation_result[j_q] = X[0];

    logger << LOG_DEBUG << "ext. res. = " << setw(15) << setprecision(8) << all_extrapolation_result[j_q] << "   chi2_dof = " << setw(15) << setprecision(8) << all_chi2_dof[j_q] << "   [" << xmin << " ; " << xmax << "]" << endl;

  }
  double temp_min = 1.e99;
  for (int i_q = start_value_xmax; i_q < end_value_xmax; i_q++){
    //  for (int i_q = start_value_xmax; i_q < xvalue.size(); i_q++){
    if (all_chi2_dof[i_q] < temp_min){temp_min = all_chi2_dof[i_q]; no_xmax = i_q;}
  }
  xmax = xvalue[no_xmax];
  
  string order_fit_name;
  if (order_fit == 1){order_fit_name = "linear";}
  else if (order_fit == 2){order_fit_name = "quadratic";}
  logger << LOG_DEBUG << "Best chi2 value from " << setw(9) << order_fit_name << " fit at " << no_xmax << " :   extrapolation_result = " << setw(15) << setprecision(8) << all_extrapolation_result[no_xmax] << "   chi2 = " << all_chi2_dof[no_xmax] << endl;
  
  double best_fit_result = all_extrapolation_result[no_xmax];
  double least_chi2 = all_chi2_dof[no_xmax];
  for (int j_q = start_value_xmax; j_q < xvalue.size(); j_q++){
    //    logger << LOG_DEBUG << "extrapolation[" << setw(2) << j_q << "] / best_fit_result = " << setw(8) << setprecision(5) << showpoint << right << (all_extrapolation_result[j_q] / best_fit_result - 1.) * 100 << " %   at relative chi2 = " << setw(8) << setprecision(5) << noshowpoint << left << all_chi2_dof[j_q] / least_chi2 << endl;
  }

  double chi2_range = ygeneric->error_extrapolation_range_chi2;
  int no_max_deviation_in_chi2_range = 0;
  double max_deviation_in_chi2_range = 0.;
  for (int j_q = start_value_xmax; j_q < xvalue.size(); j_q++){
    if (all_chi2_dof[j_q] / least_chi2 < chi2_range){
      if (abs(all_extrapolation_result[j_q] - best_fit_result) > max_deviation_in_chi2_range){
	no_max_deviation_in_chi2_range = j_q;
	max_deviation_in_chi2_range = abs(all_extrapolation_result[j_q] - best_fit_result);
      }
    }
  }
  logger << LOG_DEBUG << "Maximum deviation taken into account in " << setw(9) << order_fit_name << " fit at " << no_max_deviation_in_chi2_range << " :   rel. chi2 = " << all_chi2_dof[no_max_deviation_in_chi2_range] / least_chi2 << "   deviation = " << max_deviation_in_chi2_range << endl;

  all_extrapolation_deviation[no_xmax] = max_deviation_in_chi2_range;

  logger << LOG_DEBUG << "finished" << endl;
}



// used in distributions so far:

void summary_list::extrapolation_TSV(vector<double> & xvalue, vector<vector<vector<vector<double> > > > & data_result_TSV, vector<vector<vector<vector<double> > > > & data_deviation_TSV, vector<vector<vector<double> > > & extrapolation_result_TSV, vector<vector<vector<double> > > & extrapolation_deviation_TSV, int order_fit){
  Logger logger("summary_list::extrapolation_TSV");
  logger << LOG_DEBUG << "called" << endl;

  if (!(xvalue.size() == data_result_TSV.size())){cout << "xvalue.size() = " << xvalue.size() << " != " << data_result_TSV.size() << " = data_result_TSV.size()" << endl;}
  assert(xvalue.size() == data_result_TSV.size());
  if (!(xvalue.size() == data_deviation_TSV.size())){cout << "xvalue.size() = " << xvalue.size() << " != " << data_deviation_TSV.size() << " = data_deviation_TSV.size()" << endl;}
  assert(xvalue.size() == data_deviation_TSV.size());

  static double xmin = ygeneric->min_qTcut_extrapolation;
  static double xmax = ygeneric->max_qTcut_extrapolation;

  static int x_s = osi->no_reference_TSV;
  static int x_r = osi->no_scale_ren_reference_TSV;
  static int x_f = osi->no_scale_fact_reference_TSV;
  logger << LOG_DEBUG << "x_s = " << x_s << endl;
  logger << LOG_DEBUG << "x_r = " << x_r << endl;
  logger << LOG_DEBUG << "x_f = " << x_f << endl;
 
  vector<vector<vector<double> > > new_extrapolation_deviation_TSV(osi->n_extended_set_TSV);
  for (int i_s = 0; i_s < osi->n_extended_set_TSV; i_s++){
    if (!ygeneric->switch_output_scaleset[i_s]){continue;}
    new_extrapolation_deviation_TSV[i_s].resize(osi->n_scale_ren_TSV[i_s], vector<double> (osi->n_scale_fact_TSV[i_s], 0.));
  }


  vector<double> all_chi2_dof(xvalue.size(), 1.e99);
  vector<double> all_extrapolation_result(xvalue.size(), 0.);
  vector<double> all_extrapolation_deviation(xvalue.size(), 0.);
  vector<double> X(3, 0.);
  vector<double> dX(3, 0.);
  vector<double> data_result(xvalue.size());
  vector<double> data_deviation2(xvalue.size());

  for (int i_s = 0; i_s < extrapolation_result_TSV.size(); i_s++){
    for (int i_r = 0; i_r < extrapolation_result_TSV[i_s].size(); i_r++){
      for (int i_f = 0; i_f < extrapolation_result_TSV[i_s][i_r].size(); i_f++){

	for (int i_q = 0; i_q < xvalue.size(); i_q++){
	  data_result[i_q] = data_result_TSV[i_q][i_s][i_r][i_f];
	  data_deviation2[i_q] = pow(data_deviation_TSV[i_q][i_s][i_r][i_f], 2);
	}
	int no_xmax = 0;
	
	extrapolation_fit_range_determination(xvalue, data_result, data_deviation2, all_extrapolation_result, all_extrapolation_deviation, all_chi2_dof, order_fit, no_xmax, xmin, xmax);
	new_extrapolation_deviation_TSV[i_s][i_r][i_f] = all_extrapolation_deviation[no_xmax];

      }
    }
  }

  /*
  for (int i_q = 0; i_q < xvalue.size(); i_q++){
    data_result[i_q] = data_result_TSV[i_q][x_s][x_r][x_f];
    data_deviation2[i_q] = pow(data_deviation_TSV[i_q][x_s][x_r][x_f], 2);
  }
  int no_xmax = 0;

  extrapolation_fit_range_determination(xvalue, data_result, data_deviation2, all_extrapolation_result, all_extrapolation_deviation, all_chi2_dof, order_fit, no_xmax, xmin, xmax);
  */

  xmax = ygeneric->max_qTcut_extrapolation;
  logger << LOG_DEBUG << "Fit performed with the standard value xmax = " << xmax << endl;



  for (int i_s = 0; i_s < extrapolation_result_TSV.size(); i_s++){
    for (int i_r = 0; i_r < extrapolation_result_TSV[i_s].size(); i_r++){
      for (int i_f = 0; i_f < extrapolation_result_TSV[i_s][i_r].size(); i_f++){

	vector<double> data_result(xvalue.size());
	vector<double> data_deviation2(xvalue.size());
	for (int i_q = 0; i_q < xvalue.size(); i_q++){
	  data_result[i_q] = data_result_TSV[i_q][i_s][i_r][i_f];
	  data_deviation2[i_q] = pow(data_deviation_TSV[i_q][i_s][i_r][i_f], 2);
	}
	double chi2_dof = 0.;
	vector<double> X(3, 0.);
	vector<double> dX(3, 0.);

	if	(order_fit == 1){extrapolation_linear(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else if (order_fit == 2){extrapolation_quadratic(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
	extrapolation_result_TSV[i_s][i_r][i_f] = X[0];
	//	out_chi2[i_s][i_r][i_f] = chi2_dof;

	/*
	if (i_s == 0 && i_r == 1 && i_f == 1){
	  logger << LOG_DEBUG << "extrapolation    = " << setw(15) << setprecision(8) << extrapolation_result_TSV[i_s][i_r][i_f] << "   chi2_dof = " << chi2_dof << endl;
	  for (int i_q = 0; i_q < xvalue.size(); i_q++){
	    logger << LOG_DEBUG << "data_result[" << setw(3) << i_q << "] = " << setw(15) << setprecision(8) << data_result[i_q] << " +- " << setw(15) << setprecision(8) << sqrt(data_deviation2[i_q]) << endl;
	  }
	}
	*/
	
	double extrapolation_plus = 0.;
	for (int i_q = 0; i_q < xvalue.size(); i_q++){
	  data_result[i_q] = data_result_TSV[i_q][i_s][i_r][i_f] + data_deviation_TSV[i_q][i_s][i_r][i_f];
	}
	X.resize(3, 0.);
	dX.resize(3, 0.);
	if	(order_fit == 1){extrapolation_linear(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else if (order_fit == 2){extrapolation_quadratic(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
	extrapolation_plus = X[0];

	double extrapolation_minus = 0.;
	for (int i_q = 0; i_q < xvalue.size(); i_q++){
	  data_result[i_q] = data_result_TSV[i_q][i_s][i_r][i_f] - data_deviation_TSV[i_q][i_s][i_r][i_f];
	}
	X.resize(3, 0.);
	dX.resize(3, 0.);
	if	(order_fit == 1){extrapolation_linear(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else if (order_fit == 2){extrapolation_quadratic(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
	else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
	extrapolation_minus = X[0];

	extrapolation_deviation_TSV[i_s][i_r][i_f] = (extrapolation_plus - extrapolation_minus) / 2;
	//  pure statistical error, "extrapolated" to qTcut->0.
	//  A further error estimate from the extrapolation curve (in particular its actually random fit range) should be added !!!

	logger << LOG_INFO << "extrapolation_result_TSV[" << i_s << "][" << i_r << "][" << i_f << "] = " << setw(23) << setprecision(15) << extrapolation_result_TSV[i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << extrapolation_deviation_TSV[i_s][i_r][i_f] << " (statistic) +- " << new_extrapolation_deviation_TSV[i_s][i_r][i_f] << " (extrapolation)" << endl;


	for (int i_q = 0; i_q < xvalue.size(); i_q += 5){
	  logger << LOG_DEBUG << setw(3) << i_q << "   qTcut = " << setw(10) << setprecision(5) << xvalue[i_q] << "   result = " << setw(23) << setprecision(15) << data_result_TSV[i_q][i_s][i_r][i_f] << " +- " << setw(23) << setprecision(15) << data_deviation_TSV[i_q][i_s][i_r][i_f] << endl;
	}
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void summary_list::extrapolation_CV(vector<double> & xvalue, vector<vector<double> > & data_result_CV, vector<vector<double> > & data_deviation_CV, vector<double> & extrapolation_result_CV, vector<double> & extrapolation_deviation_CV, int order_fit){
  Logger logger("summary_list::extrapolation_CV");
  logger << LOG_DEBUG << "called" << endl;

  if (!(xvalue.size() == data_result_CV.size())){cout << "xvalue.size() = " << xvalue.size() << " != " << data_result_CV.size() << " = data_result_CV.size()" << endl;}
  assert(xvalue.size() == data_result_CV.size());
  if (!(xvalue.size() == data_deviation_CV.size())){cout << "xvalue.size() = " << xvalue.size() << " != " << data_deviation_CV.size() << " = data_deviation_CV.size()" << endl;}
  assert(xvalue.size() == data_deviation_CV.size());

  static double xmin = ygeneric->min_qTcut_extrapolation;
  static double xmax = ygeneric->max_qTcut_extrapolation;

  static int x_s = (osi->n_scales_CV - 1 ) / 2;
  logger << LOG_DEBUG << "x_s = " << x_s << endl;

  vector<double> new_extrapolation_deviation_CV(osi->n_scales_CV, 0.);

  vector<double> all_chi2_dof(xvalue.size(), 1.e99);
  vector<double> all_extrapolation_result(xvalue.size(), 0.);
  vector<double> all_extrapolation_deviation(xvalue.size(), 0.);
  vector<double> X(3, 0.);
  vector<double> dX(3, 0.);
  vector<double> data_result(xvalue.size());
  vector<double> data_deviation2(xvalue.size());


  for (int i_s = 0; i_s < extrapolation_result_CV.size(); i_s++){
    for (int i_q = 0; i_q < xvalue.size(); i_q++){
      data_result[i_q] = data_result_CV[i_q][i_s];
      data_deviation2[i_q] = pow(data_deviation_CV[i_q][i_s], 2);
    }
    int no_xmax = 0;
    logger << LOG_DEBUG << "Extrapolation at i_s = " << i_s << endl;
    extrapolation_fit_range_determination(xvalue, data_result, data_deviation2, all_extrapolation_result, all_extrapolation_deviation, all_chi2_dof, order_fit, no_xmax, xmin, xmax);

    new_extrapolation_deviation_CV[i_s] = all_extrapolation_deviation[no_xmax];
  }
  /*
  for (int i_q = 0; i_q < xvalue.size(); i_q++){
    data_result[i_q] = data_result_CV[i_q][x_s];
    data_deviation2[i_q] = pow(data_deviation_CV[i_q][x_s], 2);
  }
  int no_xmax = 0;

  extrapolation_fit_range_determination(xvalue, data_result, data_deviation2, all_extrapolation_result, all_extrapolation_deviation, all_chi2_dof, order_fit, no_xmax, xmin, xmax);
  */

  xmax = ygeneric->max_qTcut_extrapolation;
  logger << LOG_DEBUG << "Fit performed with the standard value xmax = " << xmax << endl;



  for (int i_s = 0; i_s < extrapolation_result_CV.size(); i_s++){
    vector<double> data_result(xvalue.size());
    vector<double> data_deviation2(xvalue.size());
    for (int i_q = 0; i_q < xvalue.size(); i_q++){
      data_result[i_q] = data_result_CV[i_q][i_s];
      data_deviation2[i_q] = pow(data_deviation_CV[i_q][i_s], 2);
    }
    double chi2_dof = 0.;
    vector <double> X(3, 0.);
    vector <double> dX(3, 0.);
    if	(order_fit == 1){extrapolation_linear(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
    else if (order_fit == 2){extrapolation_quadratic(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
    else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
    extrapolation_result_CV[i_s] = X[0];
    //    out_chi2[i_s] = chi2_dof;

    /*
    if (i_s == 0){
      logger << LOG_DEBUG << "extrapolation    = " << setw(15) << setprecision(8) << extrapolation_result_CV[i_s] << "   chi2_dof = " << chi2_dof << endl;
      for (int i_q = 0; i_q < xvalue.size(); i_q++){
	logger << LOG_DEBUG << "data_result[" << setw(3) << i_q << "] = " << setw(15) << setprecision(8) << data_result[i_q] << " +- " << setw(15) << setprecision(8) << sqrt(data_deviation2[i_q]) << endl;
      }
    }
    */

    double extrapolation_plus = 0.;
    for (int i_q = 0; i_q < xvalue.size(); i_q++){
      data_result[i_q] = data_result_CV[i_q][i_s] + data_deviation_CV[i_q][i_s];
    }
    X.resize(3, 0.);
    dX.resize(3, 0.);
    if	(order_fit == 1){extrapolation_linear(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
    else if (order_fit == 2){extrapolation_quadratic(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
    else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
    extrapolation_plus = X[0];
    
    double extrapolation_minus = 0.;
    for (int i_q = 0; i_q < xvalue.size(); i_q++){
      data_result[i_q] = data_result_CV[i_q][i_s] - data_deviation_CV[i_q][i_s];
    }
    X.resize(3, 0.);
    dX.resize(3, 0.);
    if	(order_fit == 1){extrapolation_linear(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
    else if (order_fit == 2){extrapolation_quadratic(xvalue, xmin, xmax, data_result, data_deviation2, X, dX, chi2_dof);}
    else {logger << LOG_FATAL << "No fit order specified." << endl; exit(1);}
    extrapolation_minus = X[0];
    
    extrapolation_deviation_CV[i_s] = (extrapolation_plus - extrapolation_minus) / 2;
    
    logger << LOG_INFO << "extrapolation_result_CV[" << i_s << "] = " << setw(23) << setprecision(15) << extrapolation_result_CV[i_s] << " +- " << setw(23) << setprecision(15) << extrapolation_deviation_CV[i_s] << " (statistic) +- " << new_extrapolation_deviation_CV[i_s] << " (extrapolation)" << endl;
    for (int i_q = 0; i_q < xvalue.size(); i_q += 5){
      logger << LOG_DEBUG << setw(3) << i_q << "   qTcut = " << setw(10) << setprecision(5) << xvalue[i_q] << "   result = " << setw(23) << setprecision(15) << data_result_CV[i_q][i_s] << " +- " << setw(23) << setprecision(15) << data_deviation_CV[i_q][i_s] << endl;
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}


