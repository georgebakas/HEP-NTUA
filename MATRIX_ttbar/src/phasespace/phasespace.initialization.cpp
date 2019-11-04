#include "../include/classes.cxx"
void phasespace_set::initialization_mapping_parameter(inputparameter_set & isi){
  Logger logger("phasespace_set::initialization_parameter (isi)");
  logger << LOG_DEBUG << "started" << endl;

  //  MC_channel = 0;
  MC_channel_phasespace = 0;
  //  tau_channel = 0;
  //  x1x2_channel = 0;

  coll_choice = isi.coll_choice;
  //  smin_opt = isi.smin_opt;
  //  sqrtsmin_opt = isi.sqrtsmin_opt;

  E = isi.E;
  s_had = pow(2 * E, 2);

  switch_qTcut = isi.switch_qTcut;
  min_qTcut = isi.min_qTcut;
  

  //  mapping_cut_pT = isi.mapping_cut_pT;
  // !!! no input parameter !!!

  nu = isi.nu;
  nuxs = isi.nuxs;
  nuxt = isi.nuxt;
  exp_pdf = isi.exp_pdf;
  exp_pT = isi.exp_pT;
  exp_y = isi.exp_y;
  exp_ij_k_y = isi.exp_ij_k_y;
  exp_ij_k_z = isi.exp_ij_k_z;
  exp_ij_a_x = isi.exp_ij_a_x;
  exp_ij_a_z = isi.exp_ij_a_z;
  exp_ai_k_x = isi.exp_ai_k_x;
  exp_ai_k_u = isi.exp_ai_k_u;
  exp_ai_b_x = isi.exp_ai_b_x;
  exp_ai_b_v = isi.exp_ai_b_v;
  map_technical_s = isi.map_technical_s;
  map_technical_t = isi.map_technical_t;
  map_technical_x = isi.map_technical_x;
  mass0 = isi.mass0;

  cut2_pT = pow(min_qTcut, 2.);
  //  cut_pT_lep = pow(mapping_cut_pT[11], 2.);
  // !!! not use any longer ???

  //  g_onshell_decay = 1.;
  //  g_onshell_decay = isi.g_onshell_decay;
  // !!! set elsewhere !!!

  cut_technical = isi.cut_technical;

  switch_console_output_phasespace_issue = isi.switch_console_output_phasespace_issue;

  logger << LOG_DEBUG << "finished" << endl;
}

// !!! temporary !!!
void phasespace_set::initialization_minimum_tau(inputparameter_set & isi){
  Logger logger("phasespace_set::initialization_minimum_tau (isi)");
  logger << LOG_DEBUG << "started" << endl;

  ////////////////////////////////////////
  //  minimum CMS-energy determination  //
  ////////////////////////////////////////
  
  // !!! should be completely worked out anew and shifted elsewhere !!!
  // !!! shoud be here, but worked out completely anew !!!

  // extend to be used also with define_ET

  double minEjet = 0.;
  // minimum energy required for jet (including b-jet (and c-jet)) system
  //  if (isi.esi.pda[isi.esi.observed_object["ljet"]].n_observed_min > 0){
  //  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["ljet"]].n_observed_min; i++){minEljet.push_back(isi.esi.pda[isi.esi.observed_object["ljet"]].define_pT);}
  vector<double> minEljet(isi.esi.pda[isi.esi.observed_object["ljet"]].n_partonlevel, 0.);
  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["ljet"]].n_observed_min; i++){minEljet[i] = isi.esi.pda[isi.esi.observed_object["ljet"]].define_pT;}
  //  }
  //  vector<double> minEbjet;
  //  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel; i++){minEbjet.push_back(M[5]);}
  vector<double> minEbjet(isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel, M[5]);
  if (isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel < isi.esi.pda[isi.esi.observed_object["bjet"]].n_observed_min){} // cross-section contribution is zero
  if (isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel == isi.esi.pda[isi.esi.observed_object["bjet"]].n_observed_min){for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel; i++){minEbjet[i] = sqrt(M2[5] + pow(isi.esi.pda[isi.esi.observed_object["bjet"]].define_pT, 2));}}
  //  else if (isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel == isi.esi.pda[isi.esi.observed_object["bjet"]].n_observed_min){minEjet += isi.esi.pda[isi.esi.observed_object["bjet"]].n_observed_min * sqrt(M2[5] + pow(isi.esi.pda[isi.esi.observed_object["bjet"]].define_pT, 2));}
  else if ((isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel > isi.esi.pda[isi.esi.observed_object["bjet"]].n_observed_min) && (isi.esi.pda[isi.esi.observed_object["bjet"]].n_observed_min > 0)){
    int min_bs_per_bjet = isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel / isi.esi.pda[isi.esi.observed_object["bjet"]].n_observed_min;
    int n_more_bs_per_bjet = isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel - isi.esi.pda[isi.esi.observed_object["bjet"]].n_observed_min * min_bs_per_bjet;
    for (int i = 0; i < min_bs_per_bjet * (isi.esi.pda[isi.esi.observed_object["bjet"]].n_observed_min - n_more_bs_per_bjet); i++){
      if (i < min_bs_per_bjet * (isi.esi.pda[isi.esi.observed_object["bjet"]].n_observed_min - n_more_bs_per_bjet)){minEbjet[i] = sqrt(M2[5] + pow(isi.esi.pda[isi.esi.observed_object["bjet"]].define_pT / min_bs_per_bjet, 2));}
      else {minEbjet[i] = sqrt(M2[5] + pow(isi.esi.pda[isi.esi.observed_object["bjet"]].define_pT / (min_bs_per_bjet + 1), 2));}
    }
  }
  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["ljet"]].n_partonlevel; i++){logger << LOG_INFO << "minEljet[" << i << "] = " << minEljet[i] << endl;}
  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel; i++){logger << LOG_INFO << "minEbjet[" << i << "] = " << minEbjet[i] << endl;}
  minEjet = accumulate(minEljet.begin(), minEljet.end(), 0.) + accumulate(minEbjet.begin(), minEbjet.end(), 0.);
  logger << LOG_INFO << "minEjet = " << minEjet << endl;
  if (user.switch_value[user.switch_map["M_jetjet"]] == 1 &&
      isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel + isi.esi.pda[isi.esi.observed_object["ljet"]].n_partonlevel >= 2){
    if (minEjet < user.cut_value[user.cut_map["M_jetjet"]]){
      minEjet = user.cut_value[user.cut_map["M_jetjet"]];
    }
  }
  logger << LOG_INFO << "minEjet = " << minEjet << endl;


  double minElep = 0.;
  vector<double> minEe(isi.esi.pda[isi.esi.observed_object["e"]].n_partonlevel, 0.);
  vector<double> minEmu(isi.esi.pda[isi.esi.observed_object["mu"]].n_partonlevel, 0.);
  vector<double> minEtau(isi.esi.pda[isi.esi.observed_object["tau"]].n_partonlevel, 0.);
  vector<double> minEalep(isi.esi.pda[isi.esi.observed_object["lep"]].n_partonlevel, 0.);
  
//  cout << "osi_observed_object[lep] = " << osi_observed_object["lep"] << endl;
  cout << "isi.esi.pda[lep]].n_observed_min = " << isi.esi.pda[isi.esi.observed_object["lep"]].n_observed_min << endl;
  cout << "isi.esi.pda[lep]].n_partonlevel = " << isi.esi.pda[isi.esi.observed_object["lep"]].n_partonlevel << endl;
  cout << "minEalep.size() = " << minEalep.size() << endl;

  if (isi.esi.pda[isi.esi.observed_object["lep"]].n_partonlevel != 0){
    for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["lep"]].n_observed_min; i++){minEalep[i] = isi.esi.pda[isi.esi.observed_object["lep"]].define_pT;}
    for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["lep"]].n_partonlevel; i++){logger << LOG_INFO << "minEalep[" << i << "] = " << minEalep[i] << endl;}
    minElep = accumulate(minEalep.begin(), minEalep.end(), 0.);
    logger << LOG_INFO << "minElep = " << minElep << endl;
    
    if (user.switch_value[user.switch_map["M_Zrec"]] == 1 &&
	isi.esi.pda[isi.esi.observed_object["lep"]].n_partonlevel >= 2){
      if (minElep < user.cut_value[user.cut_map["min_M_Zrec"]]){
	minElep = user.cut_value[user.cut_map["min_M_Zrec"]];
      }
    }

    if (user.switch_value[user.switch_map["delta_M_Zrec_MZ "]] == 1 &&
	isi.esi.pda[isi.esi.observed_object["lep"]].n_partonlevel >= 2){
      if (minElep < M[23] - user.cut_value[user.cut_map["max_delta_M_Zrec_MZ"]]){
	minElep = M[23] - user.cut_value[user.cut_map["max_delta_M_Zrec_MZ"]];
      }
    }

    if (user.switch_value[user.switch_map["M_leplep"]] == 1 &&
	isi.esi.pda[isi.esi.observed_object["lep"]].n_partonlevel >= 2){
      if (minElep < user.cut_value[user.cut_map["min_M_leplep"]]){
	minElep = user.cut_value[user.cut_map["min_M_leplep"]];
      }
    }

    if (user.switch_value[user.switch_map["M_emep"]] == 1 &&
	isi.esi.pda[isi.esi.observed_object["e"]].n_partonlevel >= 2){
      if (minElep < user.cut_value[user.cut_map["min_M_emep"]]){
	minElep = user.cut_value[user.cut_map["min_M_emep"]];
      }
    }

    if (user.switch_value[user.switch_map["M_mummup"]] == 1 &&
	isi.esi.pda[isi.esi.observed_object["mu"]].n_partonlevel >= 2){
      if (minElep < user.cut_value[user.cut_map["min_M_mummup"]]){
	minElep = user.cut_value[user.cut_map["min_M_mummup"]];
      }
    }

    if (user.switch_value[user.switch_map["M_leplep"]] == 1 &&
	isi.esi.pda[isi.esi.observed_object["lep"]].n_partonlevel >= 4){
      if (minElep < 2 * user.cut_value[user.cut_map["min_M_leplep"]]){
	minElep = 2 * user.cut_value[user.cut_map["min_M_leplep"]];
      }
    }
  }
  logger << LOG_INFO << "minElep = " << minElep << endl;

  double minEwp = 0.;
  vector<double> minEawp(isi.esi.pda[isi.esi.observed_object["wp"]].n_partonlevel, 0.);
  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["wp"]].n_observed_min; i++){minEawp[i] = sqrt(M2[24] + pow(isi.esi.pda[isi.esi.observed_object["wp"]].define_pT, 2));}
  for (int i = isi.esi.pda[isi.esi.observed_object["wp"]].n_observed_min; i < isi.esi.pda[isi.esi.observed_object["wp"]].n_partonlevel; i++){minEawp[i] = M[24];}
  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["wp"]].n_partonlevel; i++){logger << LOG_INFO << "minEawp[" << i << "] = " << minEawp[i] << endl;}
  minEwp = accumulate(minEawp.begin(), minEawp.end(), 0.);
  logger << LOG_INFO << "minEwp = " << minEwp << endl;

  double minEwm = 0.;
  vector<double> minEawm(isi.esi.pda[isi.esi.observed_object["wm"]].n_partonlevel, 0.);
  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["wm"]].n_observed_min; i++){minEawm[i] = sqrt(M2[24] + pow(isi.esi.pda[isi.esi.observed_object["wm"]].define_pT, 2));}
  for (int i = isi.esi.pda[isi.esi.observed_object["wm"]].n_observed_min; i < isi.esi.pda[isi.esi.observed_object["wm"]].n_partonlevel; i++){minEawm[i] = M[24];}
  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["wm"]].n_partonlevel; i++){logger << LOG_INFO << "minEawm[" << i << "] = " << minEawm[i] << endl;}
  minEwm = accumulate(minEawm.begin(), minEawm.end(), 0.);
  logger << LOG_INFO << "minEwm = " << minEwm << endl;


  // !!! simplified so far !!!
  // no individual cuts on different leptons

  double minEphoton = 0.;
  vector<double> minEaphoton(isi.esi.pda[isi.esi.observed_object["photon"]].n_partonlevel, 0.);
  logger << LOG_INFO << "before minEphoton = " << minEphoton << endl;
//  cout << "osi_observed_object[photon] = " << osi_observed_object["photon"] << endl;
//  cout << "osi_n_observed_min.n_observed_min.size() = " <<  osi_n_observed_min.n_observed_min.size() << endl;
//  cout << "osi_n_partonlevel.n_partonlevel.size() = " << .n_partonlevel osi_n_partonlevel.size() << endl;
  cout << "minEaphoton.size() = " << minEaphoton.size() << endl;
  cout << "isi.esi.pda[photon]].n_observed_min = " << isi.esi.pda[isi.esi.observed_object["photon"]].n_observed_min << endl;
  cout << "isi.esi.pda[photon]].define_pT = " << isi.esi.pda[isi.esi.observed_object["photon"]].define_pT << endl;

  //  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["photon"]].n_observed_min; i++){
  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["photon"]].n_observed_min; i++){minEaphoton[i] = isi.esi.pda[isi.esi.observed_object["photon"]].define_pT;}
  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["photon"]].n_partonlevel; i++){logger << LOG_INFO << "minEaphoton[" << i << "] = " << minEaphoton[i] << endl;}
  //  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["photon"]].n_partonlevel; i++){minEaphoton[i] = isi.esi.pda[isi.esi.observed_object["photon"]].define_pT;}
  logger << LOG_INFO << "after minEphoton = " << minEphoton << endl;
  //  for (int i = 0; i < isi.esi.pda[isi.esi.observed_object["photon"]].n_partonlevel; i++){logger << LOG_INFO << "minEaphoton[" << i << "] = " << minEaphoton[i] << endl;}
  logger << LOG_INFO << "before minEphoton = " << minEphoton << endl;
  minEphoton = accumulate(minEaphoton.begin(), minEaphoton.end(), 0.);
  logger << LOG_INFO << "minEphoton = " << minEphoton << endl;

  double minEmissing = 0.;
  minEmissing = isi.esi.pda[isi.esi.observed_object["missing"]].define_pT;

  double add_masses = 0.;
  for (int i = 3; i < csi->type_parton[0].size(); i++){
    if (abs(csi->type_parton[0][i]) == 6){add_masses += M[6];}
    else if (abs(csi->type_parton[0][i]) == 23){add_masses += M[23];}
    //    else if (abs(csi->type_parton[0][i]) == 24){add_masses += M[24];}
    else if (abs(csi->type_parton[0][i]) == 25){add_masses += M[25];}
  }


  
  //  logger << LOG_DEBUG << "add_masses + minEjet + minElep + minEphoton + cut_pT_miss = " << add_masses + minEjet + minElep + minEphoton + cut_pT_miss << endl;
  //  double sqrts_min_tau_0 = add_masses + minEjet + minElep + minEphoton + minEwp + minEwm + cut_pT_miss;
  
  logger << LOG_INFO << "add_masses = " << add_masses << endl;
  logger << LOG_INFO << "minEjet = " << minEjet << endl;
  logger << LOG_INFO << "minElep = " << minElep << endl;
  logger << LOG_INFO << "minEphoton = " << minEphoton << endl;
  logger << LOG_INFO << "minEwp = " << minEwp << endl;
  logger << LOG_INFO << "minEwm = " << minEwm << endl;
  logger << LOG_INFO << "minEmissing = " << minEmissing << endl;
  double sqrts_min_tau_0 = add_masses + minEjet + minElep + minEphoton + minEwp + minEwm + minEmissing;
  logger << LOG_INFO << "sqrts_min_tau_0 = " << sqrts_min_tau_0 << endl;
  if (user.switch_value[user.switch_map["HT_all"]] == 1){
    if (sqrts_min_tau_0 < user.cut_value[user.cut_map["min_HT_all"]]){
      sqrts_min_tau_0 = user.cut_value[user.cut_map["min_HT_all"]];
    }
  }

  if (user.switch_value[user.switch_map["HT_jet"]] == 1){
    if (sqrts_min_tau_0 < user.cut_value[user.cut_map["min_HT_jet"]]){
      sqrts_min_tau_0 = user.cut_value[user.cut_map["min_HT_jet"]];
    }
  }

  if (user.switch_value[user.switch_map["pT_jet_1st"]] == 1){
    if (sqrts_min_tau_0 < user.cut_value[user.cut_map["min_pT_jet_1st"]]){
      sqrts_min_tau_0 = 2 * user.cut_value[user.cut_map["min_pT_jet_1st"]];
    }
  }

  if (user.switch_value[user.switch_map["pT_w"]] == 1){
    if (sqrts_min_tau_0 < user.cut_value[user.cut_map["min_pT_w"]]){
      sqrts_min_tau_0 = 2 * user.cut_value[user.cut_map["min_pT_w"]];
    }
  }

  if (user.switch_value[user.switch_map["optimize_sqrts_min"]] == 1){
    if (sqrts_min_tau_0 < user.double_value[user.double_map["optimize_sqrts_min"]]){
      sqrts_min_tau_0 = user.double_value[user.double_map["optimize_sqrts_min"]];
    }
  }

  logger << LOG_INFO << "sqrts_min_tau_0 = " << sqrts_min_tau_0 << "   (after HT-cut)" << endl;

  tau_0 = pow((sqrts_min_tau_0) / (2 * E), 2);
  logger << LOG_INFO << "tau_0 = " << tau_0 << endl;

  if (tau_0 == 0.){tau_0 = 1.E-07; logger << LOG_INFO << "tau_0 = " << tau_0 << endl;}
  //  if (tau_0 == 0.){tau_0 = 1.E-06; logger << LOG_INFO << "tau_0 = " << tau_0 << endl;}
  
  
  /////////////////////////////////////////////////
  //  minimum CMS-energy determination finished  //
  /////////////////////////////////////////////////

  /*
  tau_0 = isi.tau_0;
  */
  tau_0_s_had = tau_0 * s_had;
  
  logger << LOG_DEBUG << "finished" << endl;
}

// !!! should be completely worked out anew and shifted elsewhere !!!
void phasespace_set::initialization_minimum_phasespacecut(inputparameter_set & isi){
  Logger logger("phasespace_set::initialization_minimum_phasespacecut (isi)");
  logger << LOG_DEBUG << "started" << endl;

  mapping_cut_pT.resize(26, 0.);
  
  if (!csi->class_contribution_CS_real){
    // only applied in configuration where each parton must become a jet:
    if (isi.esi.pda[isi.esi.observed_object["jet"]].n_partonlevel == isi.esi.pda[isi.esi.observed_object["jet"]].n_observed_min){
      double temp_pT2_jet = pow(isi.esi.pda[isi.esi.observed_object["jet"]].define_pT, 2);
      mapping_cut_pT[0] = temp_pT2_jet;
      mapping_cut_pT[1] = temp_pT2_jet;
      mapping_cut_pT[2] = temp_pT2_jet;
      mapping_cut_pT[3] = temp_pT2_jet;
      mapping_cut_pT[4] = temp_pT2_jet;
      mapping_cut_pT[5] = temp_pT2_jet;
    }
    if (isi.esi.pda[isi.esi.observed_object["ljet"]].n_partonlevel == isi.esi.pda[isi.esi.observed_object["ljet"]].n_observed_min){
      double temp_pT2_ljet = pow( isi.esi.pda[isi.esi.observed_object["ljet"]].define_pT, 2);
      if (temp_pT2_ljet > mapping_cut_pT[0]){mapping_cut_pT[0] = temp_pT2_ljet;}
      if (temp_pT2_ljet > mapping_cut_pT[0]){mapping_cut_pT[1] = temp_pT2_ljet;}
      if (temp_pT2_ljet > mapping_cut_pT[0]){mapping_cut_pT[2] = temp_pT2_ljet;}
      if (temp_pT2_ljet > mapping_cut_pT[0]){mapping_cut_pT[3] = temp_pT2_ljet;}
      if (temp_pT2_ljet > mapping_cut_pT[0]){mapping_cut_pT[4] = temp_pT2_ljet;}
    }
    if (isi.esi.pda[isi.esi.observed_object["bjet"]].n_partonlevel == isi.esi.pda[isi.esi.observed_object["bjet"]].n_observed_min){
      double temp_pT2_bjet = pow(isi.esi.pda[isi.esi.observed_object["bjet"]].define_pT, 2);
      if (temp_pT2_bjet > mapping_cut_pT[5]){mapping_cut_pT[5] = temp_pT2_bjet;}
    }
  }
  if (!csi->class_contribution_CS_real){
    if (isi.esi.pda[isi.esi.observed_object["e"]].n_partonlevel == isi.esi.pda[isi.esi.observed_object["e"]].n_observed_min){mapping_cut_pT[11] = pow(isi.esi.pda[isi.esi.observed_object["e"]].define_pT, 2);}
    if (isi.esi.pda[isi.esi.observed_object["mu"]].n_partonlevel == isi.esi.pda[isi.esi.observed_object["mu"]].n_observed_min){mapping_cut_pT[13] = pow(isi.esi.pda[isi.esi.observed_object["mu"]].define_pT, 2);}
    if (isi.esi.pda[isi.esi.observed_object["tau"]].n_partonlevel == isi.esi.pda[isi.esi.observed_object["tau"]].n_observed_min){mapping_cut_pT[15] = pow(isi.esi.pda[isi.esi.observed_object["tau"]].define_pT, 2);}
    if (isi.esi.pda[isi.esi.observed_object["photon"]].n_partonlevel == isi.esi.pda[isi.esi.observed_object["photon"]].n_observed_min){mapping_cut_pT[22] = pow(isi.esi.pda[isi.esi.observed_object["photon"]].define_pT, 2);}
    if (isi.esi.n_parton_nu == 1){
      mapping_cut_pT[12] = pow(isi.esi.pda[isi.esi.observed_object["missing"]].define_pT, 2);
      mapping_cut_pT[14] = pow(isi.esi.pda[isi.esi.observed_object["missing"]].define_pT, 2);
      mapping_cut_pT[16] = pow(isi.esi.pda[isi.esi.observed_object["missing"]].define_pT, 2);
    }
  }

  for (int i = 0; i < mapping_cut_pT.size(); i++){logger << LOG_INFO << "mapping_cut_pT[" << setw(2) << i << "] = " << mapping_cut_pT[i] << endl;}

  /////////////////////////////////////////////////////
  //  minimum phasespace-cut determination finished  //
  /////////////////////////////////////////////////////

  logger << LOG_DEBUG << "finished" << endl;
}

void phasespace_set::initialization_contribution_order(contribution_set & _csi){
  Logger logger("phasespace_set::initialization_contribution_order (isi)");
  logger << LOG_DEBUG << "started" << endl;

  //  psi.csi & = csi;
  csi = &_csi;

  logger << LOG_DEBUG << "The same information is also contained in observable_set !!!" << endl;
  // Should be removed, as csi is accessible via psi:
  process_class = _csi.process_class;
  subprocess = _csi.subprocess;
  type_perturbative_order = _csi.type_perturbative_order;
  type_contribution = _csi.type_contribution;
  type_correction = _csi.type_correction;
  //

  if (csi->type_contribution == "born" || 
      csi->type_contribution == "L2I" || 
      csi->type_contribution == "loop"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    phasespace_order_interference.resize(1, csi->contribution_order_interference);
  }
  else if (csi->type_contribution == "CA" || 
	   csi->type_contribution == "L2CA"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    if (csi->type_correction == "QCD"){
      phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s - 1);
      phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
      phasespace_order_interference.resize(1, csi->contribution_order_interference);
    }
    else if (csi->type_correction == "QEW"){
      phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
      phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e - 1);
      phasespace_order_interference.resize(1, csi->contribution_order_interference);
    }
    else if (csi->type_correction == "MIX"){
      logger << LOG_FATAL << "Contribution CA.MIX is not defined!" << endl;
      exit(1);
    }
  }
  else if (csi->type_contribution == "VA" || 
	   csi->type_contribution == "L2VA"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    if (csi->type_correction == "QCD"){
      phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s - 1);
      phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
      phasespace_order_interference.resize(1, csi->contribution_order_interference);
    }
    else if (csi->type_correction == "QEW"){
      phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
      phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e - 1);
      phasespace_order_interference.resize(1, csi->contribution_order_interference);
    }
    else if (csi->type_correction == "MIX"){
      phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
      phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e - 1);
      phasespace_order_interference.resize(1, csi->contribution_order_interference);
      logger << LOG_INFO << "Phase-space is taken according to involved QEW correction! (could be improved later...)" << endl;
    }
  }
  else if (csi->type_contribution == "RA" || 
	   csi->type_contribution == "L2RA"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    phasespace_order_interference.resize(1, csi->contribution_order_interference);
    logger << LOG_INFO << "Information on dipole phase-spaces is supplied after dipole determination..." << endl;
  }

  else if (csi->type_contribution == "VT" || 
	   csi->type_contribution == "VJ" || 
	   csi->type_contribution == "L2VT" || 
	   csi->type_contribution == "L2VJ" || 
	   csi->type_contribution == "NLL_LO" || 
	   csi->type_contribution == "NLL_NLO" || 
	   csi->type_contribution == "NNLL_LO" || 
	   csi->type_contribution == "NNLL_NLO"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s - 1);
    phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    phasespace_order_interference.resize(1, csi->contribution_order_interference);
  }
  else if (csi->type_contribution == "CT" || 
	   csi->type_contribution == "CJ" || 
	   csi->type_contribution == "L2CT" || 
	   csi->type_contribution == "L2CJ"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s - 1);
    phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    phasespace_order_interference.resize(1, csi->contribution_order_interference);
  }
  else if (csi->type_contribution == "RT" || 
	   csi->type_contribution == "L2RT" ||
	   csi->type_contribution == "RJ" || 
	   csi->type_contribution == "L2RJ"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    phasespace_order_interference.resize(1, csi->contribution_order_interference);
  }
  else if (csi->type_contribution == "VT2" || 
	   csi->type_contribution == "VJ2" || 
	   csi->type_contribution == "NNLL_NNLO"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s - 2);
    phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    phasespace_order_interference.resize(1, csi->contribution_order_interference);
  }
  else if (csi->type_contribution == "CT2" ||
	   csi->type_contribution == "CJ2"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s - 2);
    phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    phasespace_order_interference.resize(1, csi->contribution_order_interference);
  }
  else if (csi->type_contribution == "RVA" ||
	   csi->type_contribution == "RVJ"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s - 1);
    phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    phasespace_order_interference.resize(1, csi->contribution_order_interference);
  }
  else if (csi->type_contribution == "RCA" ||
	   csi->type_contribution == "RCJ"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s - 1);
    phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    phasespace_order_interference.resize(1, csi->contribution_order_interference);
  }
  else if (csi->type_contribution == "RRA" ||
	   csi->type_contribution == "RRJ"){
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
    phasespace_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    phasespace_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    phasespace_order_interference.resize(1, csi->contribution_order_interference);
    logger << LOG_INFO << "Information on dipole phase-spaces is supplied after dipole determination..." << endl;
  }
  else {
    logger << LOG_FATAL << "Contribution " << csi->type_contribution << "." << csi->type_correction << " has not been defined yet!" << endl;
      exit(1);
  }

  /*
  // for all contributions !!!
    contribution_order_alpha_s.resize(1, csi->contribution_order_alpha_s);
    contribution_order_alpha_e.resize(1, csi->contribution_order_alpha_e);
    contribution_order_interference.resize(1, csi->contribution_order_interference);
  */
  /*
  contribution_order_alpha_s = csi->contribution_order_alpha_s;
  contribution_order_alpha_e = csi->contribution_order_alpha_e;
  contribution_order_interference = csi->contribution_order_interference;
  */
  logger << LOG_DEBUG << "before contribution output" << endl;

  logger << LOG_DEBUG << "process_class = " << process_class << endl;
  logger << LOG_DEBUG << "subprocess = " << subprocess << endl;
  logger << LOG_DEBUG << "type_perturbative_order = " << type_perturbative_order << endl;
  logger << LOG_DEBUG << "type_contribution = " << csi->type_contribution << endl;
  logger << LOG_DEBUG << "type_correction = " << type_correction << endl;
  
  logger << LOG_DEBUG << "csi->contribution_order_alpha_s = " << csi->contribution_order_alpha_s << endl;
  logger << LOG_DEBUG << "csi->contribution_order_alpha_e = " << csi->contribution_order_alpha_e << endl;
  logger << LOG_DEBUG << "csi->contribution_order_interference = " << csi->contribution_order_interference << endl;
  logger << LOG_DEBUG << "contribution_order_alpha_s = " << contribution_order_alpha_s[0] << endl;
  logger << LOG_DEBUG << "contribution_order_alpha_e = " << contribution_order_alpha_e[0] << endl;
  logger << LOG_DEBUG << "contribution_order_interference = " << contribution_order_interference[0] << endl;
  logger << LOG_DEBUG << "phasespace_order_alpha_s = " << phasespace_order_alpha_s[0] << endl;
  logger << LOG_DEBUG << "phasespace_order_alpha_e = " << phasespace_order_alpha_e[0] << endl;
  logger << LOG_DEBUG << "phasespace_order_interference = " << phasespace_order_interference[0] << endl;
  
  logger << LOG_DEBUG << "after contribution output" << endl;

  // strange location !!!
  no_map.resize(1);
  o_map.resize(1);
  no_prc.resize(1);
  o_prc.resize(1);
  MC_n_channel_phasespace.resize(1);
  MC_sum_channel_phasespace.resize(1);

  logger << LOG_DEBUG << "finished" << endl;
}


void phasespace_set::initialization_complete(inputparameter_set & isi){
  Logger logger("phasespace_set::initialization_complete (isi)");
  logger << LOG_DEBUG << "started" << endl;

  if (csi->type_contribution == "born" ||
      csi->type_contribution == "RT" ||
      csi->type_contribution == "L2I" ||
      csi->type_contribution == "loop" ||
      csi->type_contribution == "L2RT"){
    initialization_phasespace_born();
    generic->ax_psp_psp(0, *this);
    MC_n_channel = generic->determination_MCchannels_psp(0, *this);
  }
  if (csi->type_contribution == "CA" ||
      csi->type_contribution == "RCA" ||
      csi->type_contribution == "L2CA"){
    initialization_phasespace_born();
    generic->ax_psp_psp(0, *this);
    MC_n_channel = generic->determination_MCchannels_psp(0, *this);
  }
  if (csi->type_contribution == "VA" ||
      csi->type_contribution == "RVA" ||
      csi->type_contribution == "L2VA"){
    initialization_phasespace_born();
    generic->ax_psp_psp(0, *this);
    MC_n_channel = generic->determination_MCchannels_psp(0, *this);
    // Check if standard VA runs are affected !!!
    if (csi->type_correction == "QCD"){
      QCD_selection_phasespace_singularity(RA_singular_region, RA_singular_region_name, RA_singular_region_list, *this);
    }
    // QEW/MIX ???
  }
  if (csi->type_contribution == "CT" ||
      csi->type_contribution == "CT2" ||
      csi->type_contribution == "VT" ||
      csi->type_contribution == "VT2"){
    initialization_phasespace_born();
    generic->ax_psp_psp(0, *this);
    MC_n_channel = generic->determination_MCchannels_psp(0, *this);
  }
  
  initialization_mapping_parameter(isi);
  initialization_minimum_tau(isi);
  initialization_minimum_phasespacecut(isi);  
  //
  
  initialization_phasespace_subprocess();
  generic->optimize_minv_psp(*this);
  initialization_phasespace_subprocess_optimization();

  generic->ac_tau_psp_psp(0, tau_MC_map, *this);

  initialization_optimization(isi);
  initialization_optimization_grid();

  //new:
  if (csi->type_contribution == "RT"){initialization_fake_dipole_mapping_RT(*this, *generic);}
  if (csi->type_contribution == "L2RT"){initialization_fake_dipole_mapping_RT(*this, *generic);}
  if (csi->type_contribution == "RCA"){initialization_fake_dipole_mapping_RT(*this, *generic);}
  if (csi->type_contribution == "RVA"){initialization_fake_dipole_mapping_RT(*this, *generic);}

  initialization_MC();
  
  if (csi->type_contribution == "VT" ||
      csi->type_contribution == "VT2"){
    initialization_phasespace_IS_CX();
    if (do_resummation){initialization_phasespace_IS_QT();}
    // Switch needs to be part of psi or csi !!!
    //    if (osi.switch_old_qT_version){initialization_phasespace_IS_CT_from_CX();}
  }
  if (csi->type_contribution == "CT" ||
      csi->type_contribution == "CT2"){
    initialization_phasespace_IS_CX();
    initialization_phasespace_IS_QT();
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}


void phasespace_set::initialization_complete_RA(inputparameter_set & isi, vector<dipole_set> & dipole){
  Logger logger("phasespace_set::initialization_complete (isi)");
  logger << LOG_DEBUG << "started" << endl;


  initialization_mapping_parameter(isi);
  initialization_minimum_tau(isi);
  initialization_minimum_phasespacecut(isi);
  
  RA_dipole = &dipole;
  n_dipoles = (*RA_dipole).size();
  
  initialization_optimization(isi);

  initialization_phasespace_subprocess_dipole_RA();
  generic->optimize_minv_psp(*this);
  initialization_phasespace_subprocess_optimization_dipole_RA();

  initialization_phasespace_RA();

  // not here for dipole-less phase-spaces !!! ???
  //  generic->ac_tau_psp_psp(0, tau_MC_map, *this);

  initialization_optimization_grid();

  // no here for dipole-less phase-spaces !!! ???
  initialization_MC();

  logger << LOG_DEBUG << "finished" << endl;
}


void phasespace_set::initialization_optimization(inputparameter_set & isi){
  Logger logger("phasespace_set::initialization_optimization (isi)");
  logger << LOG_DEBUG << "started" << endl;

  n_alpha_steps = isi.n_alpha_steps;
  n_alpha_events = isi.n_alpha_events;
  n_alpha_epc = isi.n_alpha_epc;
  n_tau_steps = isi.n_tau_steps;
  n_tau_events = isi.n_tau_events;
  n_tau_bins = isi.n_tau_bins;
  n_x1x2_steps = isi.n_x1x2_steps;
  n_x1x2_events = isi.n_x1x2_events;
  n_x1x2_bins = isi.n_x1x2_bins;
  n_z1z2_steps = isi.n_z1z2_steps;
  n_z1z2_events = isi.n_z1z2_events;
  n_z1z2_bins = isi.n_z1z2_bins;
  n_IS_events = isi.n_IS_events;
  n_IS_events_factor = isi.n_IS_events_factor;
  n_IS_steps = isi.n_IS_steps;
  n_IS_gridsize = isi.n_IS_gridsize;
  n_IS_gridsize_p = isi.n_IS_gridsize_p;
  n_IS_gridsize_f = isi.n_IS_gridsize_f;
  n_IS_gridsize_t_t = isi.n_IS_gridsize_t_t;
  n_IS_gridsize_t_phi = isi.n_IS_gridsize_t_phi;
  n_IS_gridsize_d_cth = isi.n_IS_gridsize_d_cth;
  n_IS_gridsize_d_phi = isi.n_IS_gridsize_d_phi;
  n_IS_gridsize_xy = isi.n_IS_gridsize_xy;
  n_IS_gridsize_zuv = isi.n_IS_gridsize_zuv;
  n_events_max = isi.n_events_max;
  n_events_min = isi.n_events_min;
  n_step = isi.n_step;
  
  n_events_MC_opt = 0;
  opt_n_events_min = isi.n_events_min;

  string run_mode = isi.run_mode;

  logger << LOG_DEBUG << "Input value for phase-space optimization:" << endl;

  logger << LOG_DEBUG << "switch_MC = " << switch_MC << endl;
  logger << LOG_DEBUG << "switch_MC_tau = " << switch_MC_tau << endl;
  logger << LOG_DEBUG << "switch_MC_x_dipole = " << switch_MC_x_dipole << endl;
  logger << LOG_DEBUG << "switch_IS_MC = " << switch_IS_MC << endl;
  logger << LOG_DEBUG << "switch_IS_tau = " << switch_IS_tau << endl;
  logger << LOG_DEBUG << "switch_IS_x1x2 = " << switch_IS_x1x2 << endl;
  logger << LOG_DEBUG << "switch_IS_z1z2 = " << switch_IS_z1z2 << endl;

  if (run_mode == "grid"){
    if (switch_MC == 2){switch_MC = 1;}
    if (switch_MC_tau == 2){switch_MC_tau = 1;}
    if (switch_MC_x_dipole == 2){switch_MC_x_dipole = 1;}
    if (switch_IS_MC == 2){switch_IS_MC = 1;}
    if (switch_IS_tau == 2){switch_IS_tau = 1;}
    if (switch_IS_x1x2 == 2){switch_IS_x1x2 = 1;}
    if (switch_IS_z1z2 == 2){switch_IS_z1z2 = 1;}

    if (switch_step_mode_grid){i_step_mode = & i_gen;}
    else {i_step_mode = & i_acc;}
  }
  else if (run_mode == "gridacc"){
    if (switch_MC == 2){switch_MC = 1;}
    if (switch_MC_tau == 2){switch_MC_tau = 1;}
    if (switch_MC_x_dipole == 2){switch_MC_x_dipole = 1;}
    if (switch_IS_MC == 2){switch_IS_MC = 1;}
    if (switch_IS_tau == 2){switch_IS_tau = 1;}
    if (switch_IS_x1x2 == 2){switch_IS_x1x2 = 1;}
    if (switch_IS_z1z2 == 2){switch_IS_z1z2 = 1;}
    i_step_mode = & i_acc;
  }
  else if (run_mode == "grid2"){ // for possible 2nd grid optimization step
    if (switch_MC == 2){switch_MC = 3;}
    if (switch_MC_tau == 2){switch_MC_tau = 3;}
    if (switch_MC_x_dipole == 2){switch_MC_x_dipole = 3;}
    if (switch_IS_MC == 2){switch_IS_MC = 3;}
    if (switch_IS_tau == 2){switch_IS_tau = 3;}
    if (switch_IS_x1x2 == 2){switch_IS_x1x2 = 3;}
    if (switch_IS_z1z2 == 2){switch_IS_z1z2 = 3;}

    if (switch_step_mode_grid){i_step_mode = & i_gen;}
    else {i_step_mode = & i_acc;}
  }
  else {
    i_step_mode = & i_acc;
  }

  if (!csi->class_contribution_CS_collinear &&
      !csi->class_contribution_IRcut_implicit &&
      csi->type_contribution != "VT" &&
      csi->type_contribution != "L2VT" &&
      csi->type_contribution != "VT2" &&
      csi->type_contribution != "VJ" &&
      csi->type_contribution != "L2VJ" &&
      csi->type_contribution != "VJ2"){
    switch_IS_z1z2 = -1;
  }
  if (!csi->class_contribution_CS_real){
    switch_MC_x_dipole = -1;
  }
  
  logger << LOG_DEBUG << "Values after run_mode and contribution based modifications:" << endl;
  logger << LOG_DEBUG << "run_mode = " << run_mode << endl;
  logger << LOG_DEBUG << "switch_MC = " << switch_MC << endl;
  logger << LOG_DEBUG << "switch_MC_tau = " << switch_MC_tau << endl;
  logger << LOG_DEBUG << "switch_MC_x_dipole = " << switch_MC_x_dipole << endl;
  logger << LOG_DEBUG << "switch_IS_MC = " << switch_IS_MC << endl;
  logger << LOG_DEBUG << "switch_IS_tau = " << switch_IS_tau << endl;
  logger << LOG_DEBUG << "switch_IS_x1x2 = " << switch_IS_x1x2 << endl;
  logger << LOG_DEBUG << "switch_IS_z1z2 = " << switch_IS_z1z2 << endl;
								  //  logger << LOG_DEBUG << " = " <<  << endl;
				      
  // shifted here:
  logger << LOG_DEBUG << "Values set via input in file_parameter.dat:" << endl;
  logger << LOG_DEBUG << "n_alpha_events = " << isi.n_alpha_events << endl;
  logger << LOG_DEBUG << "n_IS_events = " << isi.n_IS_events << endl;
  logger << LOG_DEBUG << "n_tau_events = " << isi.n_tau_events << endl;
  logger << LOG_DEBUG << "n_x1x2_events = " << isi.n_x1x2_events << endl;
  logger << LOG_DEBUG << "n_step = " << isi.n_step << endl;

  // shifted here:
  switch_IS_mode_phasespace = isi.switch_IS_mode_phasespace;
  weight_in_contribution = isi.MCweight_in_contribution;
  weight_in_directory = isi.MCweight_in_directory;
  weight_min = isi.MCweight_min;
  weight_limit_min = isi.MCweight_limit_min;
  weight_limit_max = isi.MCweight_limit_max;

  // doubled !!!
  initialization_filename(isi);


  logger << LOG_DEBUG << "finished" << endl;
}


void phasespace_set::initialization_optimization_grid(){
  Logger logger("phasespace_set::initialization_optimization_grid");
  logger << LOG_DEBUG << "started" << endl;


  //  MC_n_channel = isi.MC_n_channel;
  logger << LOG_DEBUG << "MC_n_channel = " << MC_n_channel << endl;
  // " has been directly set before !!!" << endl;

  // No multi-channel optimization needed if MC_n_channel == 1:
  if (MC_n_channel == 1){
    switch_MC = -1;
    n_alpha_steps = 0;
    n_alpha_events = 0;
    n_alpha_epc = 0;
  }

  logger << LOG_DEBUG << "container_IS_name.size() = " << container_IS_name.size() << endl;
  for (int i_v = 0; i_v < container_IS_name.size(); i_v++){
    logger << LOG_DEBUG << "container_IS_name[" << i_v << "] = " << container_IS_name[i_v] << endl;
  }
  // No IS sampling on multi-channel variables needed if no such variables exist:
  if (container_IS_name.size() == 0){
    switch_IS_MC = -1;
    n_IS_events = 0;
    n_IS_events_factor = 0;
    n_IS_steps = 0;
  }

  logger << LOG_DEBUG << "tau_MC_map.size() = " << tau_MC_map.size() << endl;
  // No multi-channel optimization needed for tau if tau_MC_map.size() == 1:
  if (tau_MC_map.size() == 1){
    switch_MC_tau = -1;
  }
  logger << LOG_DEBUG << "After possible tau_MC_map.size() = 1 modification:" << endl;
  logger << LOG_DEBUG << "switch_MC_tau = " << switch_MC_tau << endl;

  /*
  // shifted:
  logger << LOG_DEBUG << "Values set via input in file_parameter.dat:" << endl;
  logger << LOG_DEBUG << "n_alpha_events = " << isi.n_alpha_events << endl;
  logger << LOG_DEBUG << "n_IS_events = " << isi.n_IS_events << endl;
  logger << LOG_DEBUG << "n_tau_events = " << isi.n_tau_events << endl;
  logger << LOG_DEBUG << "n_x1x2_events = " << isi.n_x1x2_events << endl;
  logger << LOG_DEBUG << "n_step = " << isi.n_step << endl;
  */
  
  if (n_alpha_events != 0 && n_alpha_events < n_step){
    //  if (n_alpha_events < n_step){
    n_alpha_events = n_step;
    logger << LOG_INFO << "n_alpha_events is set to n_step:   " << n_alpha_events << endl;
  }
  
  if (n_tau_events < n_step){
    n_tau_events = n_step;
    logger << LOG_INFO << "n_tau_events is set to n_step:   " << n_tau_events << endl;
  }
  
  if (n_x1x2_events < n_step){
    n_x1x2_events = n_step;
    logger << LOG_INFO << "n_x1x2_events is set to n_step:   " << n_x1x2_events << endl;
  }
  
  if (n_IS_events != 0 && n_IS_events < n_step){
    n_IS_events = n_step;
    logger << LOG_INFO << "n_IS_events is set to n_step:   " << n_IS_events << endl;
  }
 
  if (csi->class_contribution_CS_collinear){
    if (n_z1z2_events < n_step){
      n_z1z2_events = n_step;
      logger << LOG_INFO << "n_z1z2_events is set to n_step:   " << n_z1z2_events << endl;
    }
  }

  if (n_step != 0){
  //double temp = MC_n_channel * n_alpha_epc;
  logger << LOG_DEBUG << "n_alpha_events             = " << n_alpha_events << endl;
  logger << LOG_DEBUG << "MC_n_channel * n_alpha_epc = " << MC_n_channel * n_alpha_epc << endl;
  if (MC_n_channel * n_alpha_epc > n_alpha_events){
    n_alpha_events = (MC_n_channel * n_alpha_epc) - (MC_n_channel * n_alpha_epc) % n_step + 1 * n_step;
    logger << LOG_INFO << "n_alpha_events is increased due to n_alpha_epc:   " << n_alpha_events << endl;
  }
  
  logger << LOG_DEBUG << "n_IS_events                         = " << n_IS_events << endl;
  logger << LOG_DEBUG << "n_alpha_events * n_IS_events_factor = " << n_alpha_events * n_IS_events_factor << endl;
  if (n_alpha_events * n_IS_events_factor > n_IS_events){
    n_IS_events = n_alpha_events * n_IS_events_factor;
    logger << LOG_INFO << "n_IS_events is increased due to n_IS_events_factor:   " << n_IS_events << endl;
  }
  
  
  
  logger << LOG_DEBUG << "n_alpha_events = " << n_alpha_events << endl;
  logger << LOG_DEBUG << "n_IS_events = " << n_IS_events << endl;
  logger << LOG_DEBUG << "n_tau_events  = " << n_tau_events << endl;
  logger << LOG_DEBUG << "n_x1x2_events  = " << n_x1x2_events << endl;
  if (csi->class_contribution_CS_collinear){
    logger << LOG_DEBUG << "n_z1z2_events  = " << n_z1z2_events << endl;
  }




  if (switch_MC == 1 || switch_MC == 3){
    if (n_alpha_steps * n_alpha_events > n_events_MC_opt){n_events_MC_opt = n_alpha_steps * n_alpha_events;}
  }
  if (switch_IS_MC == 1 || switch_IS_MC == 3){
    if (n_IS_events * n_IS_steps > n_events_MC_opt){n_events_MC_opt = n_IS_events * n_IS_steps;}
  }


  if ((switch_MC == 1 || switch_MC == 3) || (switch_IS_MC == 1 || switch_IS_MC == 3)){
    n_tau_steps = n_tau_steps + n_events_MC_opt / n_tau_events;
    //  n_tau_events = n_tau_events + n_events_MC_opt;
    logger << LOG_INFO << "n_tau_steps is increased due to MC weight optimization:   " << n_tau_steps << endl;
    
    n_x1x2_steps = n_x1x2_steps + n_events_MC_opt / n_x1x2_events;
    //  n_x1x2_events = n_x1x2_events + n_events_MC_opt;
    logger << LOG_INFO << "n_x1x2_steps is increased due to MC weight optimization:   " << n_x1x2_steps << endl;
  

    if (csi->class_contribution_CS_collinear){
      n_z1z2_steps = n_z1z2_steps + n_events_MC_opt / n_z1z2_events;
      //  n_z1z2_events = n_z1z2_events + n_events_MC_opt;
      logger << LOG_INFO << "n_z1z2_steps is increased due to MC weight optimization:   " << n_z1z2_steps << endl;
    }
  }
  
  if ((switch_MC == 1 || switch_MC == 3) && (n_events_MC_opt > opt_n_events_min)){opt_n_events_min = n_events_MC_opt;}
  if ((switch_IS_MC == 1 || switch_IS_MC == 3) && (n_events_MC_opt > opt_n_events_min)){opt_n_events_min = n_events_MC_opt;}
  if ((switch_IS_tau == 1 || switch_IS_tau == 3) && (n_tau_steps * n_tau_events > opt_n_events_min)){opt_n_events_min = n_tau_steps * n_tau_events;}
  if ((switch_IS_x1x2 == 1 || switch_IS_x1x2 == 3) && (n_x1x2_steps * n_x1x2_events > opt_n_events_min)){opt_n_events_min = n_x1x2_steps * n_x1x2_events;}
  if (csi->class_contribution_CS_collinear){
    if ((switch_IS_z1z2 == 1 || switch_IS_z1z2 == 3) && (n_z1z2_steps * n_z1z2_events > opt_n_events_min)){opt_n_events_min = n_z1z2_steps * n_z1z2_events;}
  }
  logger << LOG_INFO << "opt_n_events_min:   " << opt_n_events_min << endl;
  //  switch_n_events_opt = isi.switch_n_events_opt;
  if (switch_n_events_opt == 0){}
  else if (switch_n_events_opt == 1){
    n_events_min = opt_n_events_min;
    logger << LOG_INFO << "n_events_min is increased so that optimization phase is always finished:   " << n_events_min << endl;
  }
  else if (switch_n_events_opt == 2){
    n_events_min = opt_n_events_min;
    n_events_max = opt_n_events_min;
    logger << LOG_INFO << "n_events_min is increased so that optimization phase is always finished:   " << n_events_min << endl;
    logger << LOG_INFO << "n_events_max is set to finish of optimization phase:   " << n_events_max << endl;
  }
  else if (switch_n_events_opt == 3){
    n_events_min = opt_n_events_min;
    //    n_warmup = opt_n_events_min;
    //    warmup = 1;
    logger << LOG_INFO << "n_events_min is increased so that optimization phase is always finished:   " << n_events_min << endl;
    //    logger << LOG_INFO << "n_warmup is set to finish of optimization phase:   " << n_warmup << endl;
  }
  else {
    logger << LOG_FATAL << "switch_n_events_opt ==" << switch_n_events_opt << " is deprecated!" << endl;
    assert(false);
  }
  
  logger.newLine(LOG_INFO);
  logger << LOG_INFO << "Run values:" << endl;
  logger.newLine(LOG_INFO);
  logger << LOG_INFO << "n_alpha_events   = " << setw(10) << n_alpha_events << " in " << "n_alpha_events = " << setw(10) << n_alpha_steps << " steps" << endl;
  logger << LOG_INFO << "n_IS_events      = " << setw(10) << n_IS_events << " in " << "n_IS_events    = " << setw(10) << n_IS_steps << " steps" << endl;
  logger << LOG_INFO << "n_tau_events     = " << setw(10) << n_tau_events << " in " << "n_tau_events   = " << setw(10) << n_tau_steps << " steps" << endl;
  logger << LOG_INFO << "n_x1x2_events    = " << setw(10) << n_x1x2_events << " in " << "n_x1x2_events  = " << setw(10) << n_x1x2_steps << " steps" << endl;
  if (csi->class_contribution_CS_collinear){
    logger << LOG_INFO << "n_z1z2_events    = " << setw(10) << n_z1z2_events << " in " << "n_z1z2_events  = " << setw(10) << n_z1z2_steps << " steps" << endl;
  }
  logger.newLine(LOG_INFO);
  logger << LOG_INFO << "opt_n_events_min = " << setw(10) << opt_n_events_min << endl;
  logger.newLine(LOG_INFO);
  
  }
  /*
  // shifted:
  switch_IS_mode_phasespace = isi.switch_IS_mode_phasespace;
  weight_in_contribution = isi.MCweight_in_contribution;
  weight_in_directory = isi.MCweight_in_directory;
  weight_min = isi.MCweight_min;
  weight_limit_min = isi.MCweight_limit_min;
  weight_limit_max = isi.MCweight_limit_max;

  // doubled !!!
  initialization_filename(isi);
  */
  
  MC_phasespace = multichannel_set("MC_phasespace", MC_n_channel, n_alpha_events, n_alpha_steps, weight_min, weight_limit_min, 0.001, switch_MC, filename_MCweight, filename_MCweight_in_contribution);

  logger << LOG_DEBUG << "finished" << endl;
}

void phasespace_set::initialization_filename(inputparameter_set & isi){
  Logger logger("phasespace_set::initialization_filename");
  logger << LOG_DEBUG << "started" << endl;

  string directory_MCweights_in_contribution = "../../../../" + isi.MCweight_in_directory + "/weights";
  // not clear yet !!!



  int switch_output_weights = 1; // temporary !!!

  //  assert(isi.csi.subprocess == csi->subprocess);

  string dir_MCweights = "weights";

  filename_MCweight = dir_MCweights + "/weights_MC_" + csi->subprocess + ".dat";
  filename_tauweight = dir_MCweights + "/weights_tau_" + csi->subprocess + ".dat";
  filename_x1x2weight = dir_MCweights + "/weights_x1x2_" + csi->subprocess + ".dat";
  filename_z1z2weight.resize(3);
  filename_z1z2weight[1] = dir_MCweights + "/weights_z1_" + csi->subprocess + ".dat";
  filename_z1z2weight[2] = dir_MCweights + "/weights_z2_" + csi->subprocess + ".dat";

  //  string dir_MCweights_in_contribution = "../../../../" + MCweight_in_directory + "/weights";

  filename_MCweight_in_contribution = directory_MCweights_in_contribution + "/weights_MC_" + csi->subprocess + ".dat";
  filename_tauweight_in_contribution = directory_MCweights_in_contribution + "/weights_tau_" + csi->subprocess + ".dat";
  filename_x1x2weight_in_contribution = directory_MCweights_in_contribution + "/weights_x1x2_" + csi->subprocess + ".dat";
  filename_z1z2weight_in_contribution.resize(3);
  filename_z1z2weight_in_contribution[1] = directory_MCweights_in_contribution + "/weights_z1_" + csi->subprocess + ".dat";
  filename_z1z2weight_in_contribution[2] = directory_MCweights_in_contribution + "/weights_z2_" + csi->subprocess + ".dat";

  if (switch_output_weights){system_execute(logger, "mkdir " + dir_MCweights);}

  logger << LOG_DEBUG << "filename_MCweight   = " << filename_MCweight << endl;
  logger << LOG_DEBUG << "filename_tauweight  = " << filename_tauweight << endl;
  logger << LOG_DEBUG << "filename_x1x2weight = " << filename_x1x2weight << endl;
  logger << LOG_DEBUG << "filename_MCweight_in_contribution   = " << filename_MCweight_in_contribution << endl;
  logger << LOG_DEBUG << "filename_tauweight_in_contribution  = " << filename_tauweight_in_contribution << endl;
  logger << LOG_DEBUG << "filename_x1x2weight_in_contribution = " << filename_x1x2weight_in_contribution << endl;

  logger << LOG_DEBUG << "filenames created" << endl;

  logger << LOG_DEBUG << "finished" << endl;
}



void phasespace_set::initialization_masses(model_set & msi){
  Logger logger("phasespace_set::initialization_masses");
  logger << LOG_DEBUG << "started" << endl;

  g_global_NWA = 1.; // temporarily !!! does not belong here !!!

  // Extension to improve mappings in loop-induced processes:

  M.resize(50, 0.);
  M2.resize(50, 0.);
  Gamma.resize(50, 0.);
  cM2.resize(50, 0.);
  map_Gamma.resize(50, 0.);
  reg_Gamma.resize(50, 0.);

  for (int i = 0; i < 26; i++){ 
    M[i] = msi.M[i];
    M2[i] = msi.M2[i];
    Gamma[i] = msi.Gamma[i];
    cM2[i] = msi.cM2[i];
    map_Gamma[i] = msi.map_Gamma[i];
    reg_Gamma[i] = msi.reg_Gamma[i];
  }

  for (int i = 1; i < 7; i++){
    if (M[i] != 0.){
      M[30 + i] = 2 * M[i];
      M2[30 + i] = pow(M2[i], 2);
      Gamma[30 + i] = 10 * Gamma[i];
      cM2[30 + i] = pow(M[30 + i], 2) - ri * M[30 + i] * Gamma[30 + i];
      map_Gamma[30 + i] = 10 * map_Gamma[i];
      map_Gamma[30 + i] = 10 * map_Gamma[i];

      M[40 + i] = 2 * M[i];
      M2[40 + i] = pow(M2[i], 2);
      Gamma[40 + i] = 10 * Gamma[i];
      cM2[40 + i] = pow(M[40 + i], 2) - ri * M[40 + i] * Gamma[40 + i];
      map_Gamma[40 + i] = 10 * map_Gamma[i];
      map_Gamma[40 + i] = 10 * map_Gamma[i];
    }
  }

  /*    
  // old version: should be identical for physical particles !!!
  M = msi.M;
  M2 = msi.M2;
  Gamma = msi.Gamma;
  cM2 = msi.cM2;
  map_Gamma = msi.map_Gamma;
  reg_Gamma = msi.reg_Gamma;
  */
  
  logger << LOG_DEBUG << "M.size() = " << M.size() << endl;
  logger << LOG_DEBUG << "M2.size() = " << M2.size() << endl;
  logger << LOG_DEBUG << "cM2.size() = " << cM2.size() << endl;
  logger << LOG_DEBUG << "Gamma.size() = " << Gamma.size() << endl;
  logger << LOG_DEBUG << "map_Gamma.size() = " << map_Gamma.size() << endl;
  logger << LOG_DEBUG << "reg_Gamma.size() = " << reg_Gamma.size() << endl;

  for (int i = 0; i < M.size(); i++){logger << LOG_DEBUG << "M[" << setw(2) << i << "]  = " << right << setprecision(15) << setw(23) << M[i] << "   Gamma[" << setw(2) << i << "]  = " << right << setprecision(15) << setw(23) << Gamma[i] << "   map = " << right << setprecision(8) << setw(8) << map_Gamma[i] << "   reg = " << right << setprecision(8) << setw(8) << reg_Gamma[i] << endl;}

  logger << LOG_DEBUG << "finished" << endl;
}



void phasespace_set::initialization_masses(vector<double> _M, vector<double> _M2, vector<double> _Gamma, vector<double_complex> _cM2, vector<double> _map_Gamma, vector<double> _reg_Gamma){
  Logger logger("phasespace_set::initialization_masses");
  logger << LOG_DEBUG << "started" << endl;

  g_global_NWA = 1.; // temporarily !!! does not belong here !!!

  M = _M;
  M2 = _M2;
  Gamma = _Gamma;
  cM2 = _cM2;
  map_Gamma = _map_Gamma;
  reg_Gamma = _reg_Gamma;

  logger << LOG_DEBUG << "M.size() = " << M.size() << endl;
  logger << LOG_DEBUG << "M2.size() = " << M2.size() << endl;
  logger << LOG_DEBUG << "cM2.size() = " << cM2.size() << endl;
  logger << LOG_DEBUG << "Gamma.size() = " << Gamma.size() << endl;
  logger << LOG_DEBUG << "map_Gamma.size() = " << map_Gamma.size() << endl;
  logger << LOG_DEBUG << "reg_Gamma.size() = " << reg_Gamma.size() << endl;

  for (int i = 0; i < 26; i++){logger << LOG_DEBUG << "M[" << setw(2) << i << "]  = " << right << setprecision(15) << setw(23) << M[i] << "   Gamma[" << setw(2) << i << "]  = " << right << setprecision(15) << setw(23) << Gamma[i] << "   map = " << right << setprecision(8) << setw(8) << map_Gamma[i] << "   reg = " << right << setprecision(8) << setw(8) << reg_Gamma[i] << endl;}

  logger << LOG_DEBUG << "finished" << endl;
}



void phasespace_set::initialization_phasespace_born(){
  Logger logger("phasespace_set::initialization_phasespace_born");
  logger << LOG_DEBUG << "started" << endl;

  n_dipoles = 1;

  container_IS_startvalue.resize(1, vector<int> (4));

  c_p.resize(1);
  c_f.resize(1);
  c_t.resize(1);
  c_d.resize(1);

  v_smin.resize(1);
  v_smax.resize(1);

  g_p.resize(1);
  g_f.resize(1);
  g_t.resize(1);
  g_d.resize(1);

  needed_v_smin.resize(1);
  needed_v_smax.resize(1);
  needed_g_p.resize(1);
  
  needed_g_f.resize(1);
  needed_g_t.resize(1);
  needed_g_d.resize(1);
  
  //  g_IS_all_RA.resize(1);

  // from initialization.particles.QEW.RA.cxx
  start_xbp_all.resize(1);
  start_xbs_all.resize(1);
  start_xbsqrts_all.resize(1);

  MC_optswitch = 1;

  logger << LOG_DEBUG << "finished" << endl;
}



void phasespace_set::initialization_phasespace_RA(){
  //void phasespace_set::initialization_phasespace_RA(int _n_dipoles, vector<dipole_set> _dipole){
  Logger logger("phasespace_set::initialization_phasespace_RA");
  logger << LOG_DEBUG << "started" << endl;

  logger << LOG_DEBUG << "before filling of dipole entries:" << endl;
  logger << LOG_DEBUG << "no_map.size() = " << no_map.size() << endl;
  logger << LOG_DEBUG << "o_map.size() = " << o_map.size() << endl;
  logger << LOG_DEBUG << "no_prc.size() = " << no_prc.size() << endl;
  logger << LOG_DEBUG << "o_prc.size() = " << o_prc.size() << endl;
  for (int i_a = 0; i_a < no_map.size(); i_a++){
    stringstream temp_map;
    temp_map << "no_map[" << i_a << "] = " << no_map[i_a] << "     ";
    temp_map << "o_map[" << i_a << "] = ";
    for (int i_p = 0; i_p < o_map[i_a].size(); i_p++){
      temp_map << o_map[i_a][i_p] << "  ";
    }
    logger << LOG_DEBUG << temp_map.str() << endl;
    stringstream temp_prc;
    temp_prc << "no_prc[" << i_a << "] = " << no_prc[i_a] << "     ";
    temp_prc << "o_prc[" << i_a << "] = ";
    for (int i_p = 0; i_p < o_prc[i_a].size(); i_p++){
      temp_prc << o_prc[i_a][i_p] << "  ";
    }
    logger << LOG_DEBUG << temp_prc.str() << endl;
  }
  logger << LOG_DEBUG << endl;
  
  //  n_dipoles = _n_dipoles;

  no_map.resize(n_dipoles);
  o_map.resize(n_dipoles);
  no_prc.resize(n_dipoles);
  o_prc.resize(n_dipoles);
  for (int i_a = 0; i_a < n_dipoles; i_a++){
    no_map[i_a] = (*RA_dipole)[i_a].no_map();
    o_map[i_a] = (*RA_dipole)[i_a].o_map();
    no_prc[i_a] = (*RA_dipole)[i_a].no_prc();
    o_prc[i_a] = (*RA_dipole)[i_a].o_prc();
  }
  logger << LOG_DEBUG << "after filling of dipole entries:" << endl;
  logger << LOG_DEBUG << "no_map.size() = " << no_map.size() << endl;
  logger << LOG_DEBUG << "o_map.size() = " << o_map.size() << endl;
  logger << LOG_DEBUG << "no_prc.size() = " << no_prc.size() << endl;
  logger << LOG_DEBUG << "o_prc.size() = " << o_prc.size() << endl;
  for (int i_a = 0; i_a < no_map.size(); i_a++){
    stringstream temp_map;
    temp_map << "no_map[" << i_a << "] = " << no_map[i_a] << "     ";
    temp_map << "o_map[" << i_a << "] = ";
    for (int i_p = 0; i_p < o_map[i_a].size(); i_p++){
      temp_map << o_map[i_a][i_p] << "  ";
    }
    logger << LOG_DEBUG << temp_map.str() << endl;
    stringstream temp_prc;
    temp_prc << "no_prc[" << i_a << "] = " << no_prc[i_a] << "     ";
    temp_prc << "o_prc[" << i_a << "] = ";
    for (int i_p = 0; i_p < o_prc[i_a].size(); i_p++){
      temp_prc << o_prc[i_a][i_p] << "  ";
    }
    logger << LOG_DEBUG << temp_prc.str() << endl;
  }

  logger << LOG_DEBUG << endl;

  container_IS_name.resize(0);
  container_IS_switch.resize(0);

  container_IS_startvalue.resize(n_dipoles, vector<int> (7));
  container_IS_startvalue[0].resize(4);
  c_p.resize(n_dipoles);
  c_f.resize(n_dipoles);
  c_t.resize(n_dipoles);
  c_d.resize(n_dipoles);

  v_smin.resize(n_dipoles);
  v_smax.resize(n_dipoles);

  g_p.resize(n_dipoles);
  g_f.resize(n_dipoles);
  g_t.resize(n_dipoles);
  g_d.resize(n_dipoles);

  needed_v_smin.resize(n_dipoles);
  needed_v_smax.resize(n_dipoles);
  needed_g_p.resize(n_dipoles);
  
  needed_g_f.resize(n_dipoles);
  needed_g_t.resize(n_dipoles);
  needed_g_d.resize(n_dipoles);
  
  //  g_IS_all_RA.resize(n_dipoles);

  // from initialization.particles.QEW.RA.cxx
  /*
  start_xbp_all.resize(n_dipoles);
  start_xbs_all.resize(n_dipoles);
  start_xbsqrts_all.resize(n_dipoles);
  */
  
  //  MC_x_dipole_mapping.resize(n_dipoles, vector<int> (1, 0));
  MC_x_dipole_mapping.resize(n_dipoles);

  //  dipole_sinx_min.resize(n_dipoles, 0.);
  dipole_x.resize(n_dipoles, 0.);

  no_random_dipole.resize(3);
  no_random_dipole[0] = 3 * csi->n_particle - 4 - 2;
  no_random_dipole[1] = 3 * csi->n_particle - 4 - 1;
  no_random_dipole[2] = 3 * csi->n_particle - 4;

  MC_optswitch = 0;

  //  logger << LOG_DEBUG << "MC_n_channel = " << MC_n_channel << endl;
  for (int i_a = 0; i_a < n_dipoles; i_a++){
    logger << LOG_DEBUG << "MC_n_channel_phasespace[" << setw(2) << i_a << "] = " << setw(4) << MC_n_channel_phasespace[i_a] << "   MC_sum_channel_phasespace[" << setw(2) << i_a << "] = " << setw(4) << MC_sum_channel_phasespace[i_a] << endl;
  }
  MC_n_channel = MC_sum_channel_phasespace[n_dipoles - 1];
  /*
  for (int i_a = 0; i_a < n_dipoles; i_a++){
    MC_n_channel += _dipole[i_a].n_channel();
    logger << LOG_DEBUG << "dipole[" << setw(2) << i_a << "].n_channel() = " << _dipole[i_a].n_channel() << endl;
  }
  */
  logger << LOG_DEBUG << "MC_n_channel = " << MC_n_channel << endl;

  // shifted from munich.integration.RA.QCD.cpp
  // if-statement should be redundant !!!
  if (csi->type_contribution == "RA" ||
      csi->type_contribution == "RRA" ||
      csi->type_contribution == "RRJ" ||
      csi->type_contribution == "L2RA"){
    logger << LOG_DEBUG << "start_xbs_all.size() = " << start_xbs_all.size() << endl;
    // Simplification in case there is no contributing dipole:
    if (n_dipoles == 1){switch_RS_mapping = 1;}
    for (int i_a = 1; i_a < n_dipoles; i_a++){
      logger << LOG_DEBUG << "start_xbs_all[" << i_a << "].size() = " << start_xbs_all[i_a].size() << endl;
      if (start_xbs_all[i_a][xb_max / 2 - 4]){switch_RS_mapping = 1;}
    }
  }
  logger << LOG_DEBUG << "switch_RS_mapping = " << switch_RS_mapping << endl;

  // Switch of dipole mapping with 2->1 kinematics
  // Otherwise, delta function need to be implemented for x-type dipole variables !!!
  if (switch_RS_mapping){MC_n_channel = MC_n_channel_phasespace[0];}

  generic->ax_psp_psp(0, *this);

  if (!switch_RS_mapping){
    //  if (!osi.switch_RS_mapping){
    for (int i_a = 1; i_a < n_dipoles; i_a++){generic->ax_psp_dipole(i_a, *this);}
    initialization_phasespace_subprocess_dipole_IS_RA();
    //    initialization_phasespace_subprocess_dipole_IS_RA(n_dipoles, _dipole);
  }

  // first if-statement should be redundant !!!
  if (csi->class_contribution_CS_real){
    if (!switch_RS_mapping){
      for (int i_a = 1; i_a < n_dipoles; i_a++){generic->ac_tau_psp_dipole(i_a, MC_x_dipole_mapping[i_a], *this);}
      if (csi->type_contribution == "RRA" ||
	  csi->type_contribution == "RRJ"){
	initialization_fake_dipole_mapping_RRA();
	//	initialization_fake_dipole_mapping_RRA(dipole, *this, generic);
      }
      initialization_phasespace_MC_x_dipole_RA();
    }
    else if (switch_RS_mapping){
      generic->ac_tau_psp_psp(0, MC_x_dipole_mapping[0], *this);
      tau_MC_map = MC_x_dipole_mapping[0];
    }
  }
  
  logger << LOG_DEBUG << "switch_RS_mapping = " << switch_RS_mapping << endl;
  logger << LOG_DEBUG << "MC_n_channel = " << MC_n_channel << endl;
  logger << LOG_DEBUG << "tau_MC_map.size() = " << tau_MC_map.size() << endl;

  
  logger << LOG_DEBUG << "finished" << endl;
}



void phasespace_set::initialization_phasespace_subprocess(){
  Logger logger("phasespace_set::initialization_phasespace_subprocess");
  logger << LOG_DEBUG << "started" << endl;

  x_pdf.resize(3);
  vector<double> pdf_factor(3);

  logger << LOG_DEBUG_VERBOSE << "start_xbp_all.size() = " << start_xbp_all.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "start_xbs_all.size() = " << start_xbs_all.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "start_xbsqrts_all.size() = " << start_xbsqrts_all.size() << endl;

  start_xbp_all[0].resize(xb_max);
  start_xbs_all[0].resize(xb_max);
  start_xbsqrts_all[0].resize(xb_max);


  if (coll_choice == 0){
    fourvector p1(E, 0, 0, E), p2(E, 0, 0, -E);
    start_xbp_all[0][0] = p1 + p2;
    start_xbp_all[0][1] = p1;
    start_xbp_all[0][2] = p2;
    //    s_hat = 4. * pow(E, 2);
    //    sqrt_s_hat = 2. * E;
    for (int i = 0; i < 3; i++){x_pdf[i] = 1.;}
  }

  //  double s_hadronic = 4. * pow(E, 2);
  logger << LOG_DEBUG << "definitions of particle momenta and masses" << endl;

  // **************************************************************************
  // *                                                                        *
  // *  definitions for automatized phase-space generation                    *
  // *                                                                        *
  // **************************************************************************

  //  int xb_max = intpow(2, csi->n_particle + 2);
  if (csi->n_particle > 0){
    for (int xi = 1; xi <= csi->n_particle + 2; xi++){
      start_xbs_all[0][intpow(2, xi - 1)] = M2[abs(csi->type_parton[0][xi])];
      start_xbsqrts_all[0][intpow(2, xi - 1)] = M[abs(csi->type_parton[0][xi])];
    }
    if (coll_choice == 0){
      start_xbs_all[0][0] = start_xbp_all[0][0].m2();
      start_xbsqrts_all[0][0] = sqrt(start_xbs_all[0][0]);
      start_xbp_all[0][xb_max - 4] = start_xbp_all[0][0];
      start_xbs_all[0][xb_max - 4] = start_xbs_all[0][0];
      start_xbsqrts_all[0][xb_max - 4] = start_xbsqrts_all[0][0];
      start_xbp_all[0][3] = start_xbp_all[0][0];
      start_xbs_all[0][3] = start_xbs_all[0][0];
      start_xbsqrts_all[0][3] = start_xbsqrts_all[0][0];
      pdf_factor[0] = 1.;
      pdf_factor[1] = 1.;
      pdf_factor[2] = 0.;
      //      tau_channel = 0;
    }
  }
  sqrtsmin_opt = start_xbsqrts_all; // = start_xbsqrts;
  smin_opt = start_xbs_all; // = start_xbs;

  /*
  vector<double> xbsmin_opt = start_xbs_all[0];
  vector<double> xbsqrtsmin_opt = start_xbsqrts_all[0];
  // new:
  vector<vector<double> > xbsqrtsmin_opt_all(n_ps, xbsqrtsmin_opt);
  vector<vector<double> > xbsmin_opt_all(n_ps, xbsmin_opt);
  // new.
  */
  
  logger << LOG_DEBUG << "finished" << endl;
}

void phasespace_set::initialization_phasespace_subprocess_optimization(){
  Logger logger("phasespace_set::initialization_phasespace_subprocess_optimization");
  logger << LOG_DEBUG << "started" << endl;

  for (int ix = 4; ix < sqrtsmin_opt[0].size(); ix += 4){
    if (sqrtsmin_opt[0][ix] == start_xbsqrts_all[0][ix]){}
    else if (sqrtsmin_opt[0][ix] < start_xbsqrts_all[0][ix]){sqrtsmin_opt[0][ix] = start_xbsqrts_all[0][ix];}
    else if (sqrtsmin_opt[0][ix] > start_xbsqrts_all[0][ix] && vectorbinary_from_binary(ix).size() == 1){
      cout << "phase-space cut is in contradiction to on-shell conditions of intermediate particle " << ix << "!" << endl;
      exit(1);
    }
  }

  for (int ix = 4; ix < sqrtsmin_opt[0].size(); ix += 4){
    if (vectorbinary_from_binary(ix).size() > 1){
      double temp_sqrts_opt = max_subset_from_binary(sqrtsmin_opt[0], ix);
      if (temp_sqrts_opt > sqrtsmin_opt[0][ix]){sqrtsmin_opt[0][ix] = temp_sqrts_opt;}
      smin_opt[0][ix] = pow(sqrtsmin_opt[0][ix], 2);
    }
  }

  logger << LOG_INFO << "OPT   list min_opt[0]" << endl;

  //cout << "smin_opt[0] before adding invariant-mass cuts:" << endl;
  for (int i = 4; i < smin_opt[0].size(); i += 4){
    //  logger << LOG_INFO << "sqrtsmin_opt[0][" << setw(3) << right << i << "] = " << setdl << sqrtsmin_opt[0][i] << "   smin_opt[0][" << setw(3) << right << i << "] = " << setdl << smin_opt[0][i] << endl;
    stringstream temp_ss;
    temp_ss << "sqrtsmin_opt[0][" << setw(3) << i << "] = " << setw(15) << setprecision(8) << sqrtsmin_opt[0][i] << "   smin_opt[0][" << setw(3) << i << "] = " << setw(15) << setprecision(8) << smin_opt[0][i] << endl;
    /*
    << "   ---   ";// << endl;
    vector<int> temp_vb = vectorbinary_from_binary(i);
    if (temp_vb.size() > 1){
      for (int j = 0; j < temp_vb.size(); j++){
	temp_ss << "[" << setw(3) << temp_vb[j] << "] = " << setw(15) << setprecision(8) << sqrtsmin_opt[0][temp_vb[j]] << "   ";
      }
    }
    temp_ss << endl;
    */
    if (sqrtsmin_opt[0][i] != 0.){logger << LOG_INFO << temp_ss.str();}
    else {logger << LOG_DEBUG << temp_ss.str();}
  }
  logger.newLine(LOG_INFO);
  
  for (int i = 0; i < start_xbp_all[0].size(); i++){
    if (start_xbs_all[0][i] != 0. || start_xbsqrts_all[0][i] != 0.){
      logger << LOG_INFO << "xbs[" << setw(3) << right << i << "] = " << start_xbs_all[0][i] << "   " << "xbsqrts[" << setw(3) << right << i << "] = " << start_xbsqrts_all[0][i] << "   " << endl;
    }
  }
  logger.newLine(LOG_INFO);

  logger << LOG_DEBUG << "tau_0 = " << tau_0 << "   -> sqrt(s^_min) = " << 2 * sqrt(tau_0) * E << endl;
  if (tau_0 < smin_opt[0][smin_opt[0].size() - 4] / (4 * pow(E, 2))){tau_0 = smin_opt[0][smin_opt[0].size() - 4] / (4 * pow(E, 2));}
  //  if (tau_0 < pow(cut_Mjetjet, 2) / (4 * pow(E, 2))){tau_0 = pow(cut_Mjetjet, 2) / (4 * pow(E, 2));}
//if (tau_0 < pow(100., 2) / (4 * pow(E, 2))){tau_0 = pow(100., 2) / (4 * pow(E, 2));}
logger << LOG_DEBUG << "tau_0 = " << tau_0 << "   -> sqrt(s^_min) = " << 2 * sqrt(tau_0) * E << endl;

  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "finished" << endl;
}

void phasespace_set::initialization_phasespace_subprocess_dipole_IS_RA(){
  //void phasespace_set::initialization_phasespace_subprocess_dipole_IS_RA(int _n_dipoles, vector<dipole_set> _dipole){
  Logger logger("phasespace_set::initialization_phasespace_subprocess_dipole_IS_RA");
  logger << LOG_DEBUG << "started" << endl;

  for (int i_a = 1; i_a < n_dipoles; i_a++){
    container_IS_startvalue[i_a][4] = container_IS_name.size();
    container_IS_name.push_back("x of dipole " + (*RA_dipole)[i_a].name());
    container_IS_startvalue[i_a][5] = container_IS_name.size();
    container_IS_name.push_back("z of dipole " + (*RA_dipole)[i_a].name());
    container_IS_startvalue[i_a][6] = container_IS_name.size();
    container_IS_name.push_back("phi of dipole " + (*RA_dipole)[i_a].name());
  }
  logger << LOG_DEBUG << "finished" << endl;
}

void phasespace_set::initialization_phasespace_subprocess_dipole_RA(){
  //void phasespace_set::initialization_phasespace_subprocess_dipole_RA(int _n_dipoles, vector<dipole_set> _dipole){
  Logger logger("phasespace_set::initialization_phasespace_subprocess_dipole_RA");
  logger << LOG_DEBUG << "started" << endl;

  /*
  for (int i_a = 1; i_a < n_dipoles; i_a++){
    container_IS_startvalue[i_a][4] = container_IS_name.size();
    container_IS_name.push_back("x of dipole " + (*RA_dipole)[i_a].name());
    container_IS_startvalue[i_a][5] = container_IS_name.size();
    container_IS_name.push_back("z of dipole " + (*RA_dipole)[i_a].name());
    container_IS_startvalue[i_a][6] = container_IS_name.size();
    container_IS_name.push_back("phi of dipole " + (*RA_dipole)[i_a].name());
  }
  */
  // from initialization.particles.QEW.RA.cxx
  // temporary - check if needed elsewhere !!!
  //  double s_hat;
  //  double sqrt_s_hat;
  //  int max_particles;
  //  int tau_channel;
  //  int coll_choice = 1;

  x_pdf.resize(3);
  vector<double> pdf_factor(3);
  vector<double> p_2(3 + csi->n_particle);
  vector<double> sqp_2(3 + csi->n_particle);
  /*  
  vector<vector<double> > p_2_RA(n_dipoles);
  vector<vector<double> > sqp_2_RA(n_dipoles);
  */
  vector<vector<vector<double> > > dx_s_RA(n_dipoles);
  vector<vector<vector<double> > > dx_sqrts_RA(n_dipoles);
  for (int i_a = 0; i_a < dx_s_RA.size(); i_a++){
    dx_s_RA[i_a].resize((*RA_dipole)[i_a].type_parton().size(), vector<double> (1));
    dx_sqrts_RA[i_a].resize((*RA_dipole)[i_a].type_parton().size(), vector<double> (1));
  }
  logger << LOG_DEBUG << "x_s_RA started" << endl;



  logger << LOG_DEBUG << "xb_max = " << xb_max << endl;

  start_xbp_all.resize(n_dipoles);
  start_xbs_all.resize(n_dipoles);
  start_xbsqrts_all.resize(n_dipoles);

  logger << LOG_DEBUG << "start_xbp_all.size() = " << start_xbp_all.size() << endl;
  start_xbp_all[0].resize(xb_max);
  start_xbs_all[0].resize(xb_max);
  start_xbsqrts_all[0].resize(xb_max);
//  logger << LOG_DEBUG << "   n_dipoles = " << n_dipoles << endl;
  for (int c = 1; c < n_dipoles; c++){
    start_xbp_all[c].resize(xb_max_dipoles);
    start_xbs_all[c].resize(xb_max_dipoles);
    start_xbsqrts_all[c].resize(xb_max_dipoles);
  }
  logger << LOG_DEBUG << "xb_max = " << xb_max << endl;



  p_2.resize(3 + csi->n_particle);
  sqp_2.resize(3 + csi->n_particle);
  for (int i = 1; i < csi->n_particle + 3; i++){
    logger << LOG_DEBUG << "i = " << i << endl;
    logger << LOG_DEBUG << "(*RA_dipole)[0].type_parton().size() = " << (*RA_dipole)[0].type_parton().size() << endl;
    logger << LOG_DEBUG << "(*RA_dipole)[0].type_parton()[" << i << "] = " << (*RA_dipole)[0].type_parton()[i] << endl;
    logger << LOG_DEBUG << "M.size() = " << M.size() << endl;
    logger << LOG_DEBUG << "M2.size() = " << M2.size() << endl;
    logger << LOG_DEBUG << "M[abs((*RA_dipole)[0].type_parton()[" << i << ")] = " << (*RA_dipole)[0].type_parton()[i] << "] = " << setw(23) << setprecision(15) << M[abs((*RA_dipole)[0].type_parton()[i])] << endl;
    logger << LOG_DEBUG << "M2[abs((*RA_dipole)[0].type_parton()[" << i << ")] = " << (*RA_dipole)[0].type_parton()[i] << "] = " << setw(23) << setprecision(15) << M2[abs((*RA_dipole)[0].type_parton()[i])] << endl;
    
    p_2[i] = M2[abs((*RA_dipole)[0].type_parton()[i])];
    sqp_2[i] = M[abs((*RA_dipole)[0].type_parton()[i])];
    //    p_2[i] = M2[abs(csi->type_parton[0][i])];
    //    sqp_2[i] = M[abs(csi->type_parton[0][i])];
  }
  sqp_2[0] = accumulate(sqp_2.begin(), sqp_2.end(), 0.); // ppwwbb
  p_2[0] = pow(sqp_2[0], 2);
   
  logger << LOG_DEBUG << "p_2 started" << endl;
    
  for (int i_a = 0; i_a < dx_s_RA.size(); i_a++){
    for (int i_b = 1; i_b < dx_s_RA[i_a].size(); i_b++){
      if (dx_s_RA[i_a][i_b].size() == 1){
	dx_s_RA[i_a][i_b][0] = M2[abs((*RA_dipole)[i_a].type_parton()[i_b])];
	dx_sqrts_RA[i_a][i_b][0] = M[abs((*RA_dipole)[i_a].type_parton()[i_b])];
      }
    }
  }
 




  logger << LOG_INFO << "definitions of particle momenta and masses" << endl;
// **************************************************************************
// *                                                                        *
// *  definitions for automatized phase-space generation                    *
// *                                                                        *
// **************************************************************************
  for (int i_a = 0; i_a < n_dipoles; i_a++){
    for (int xbi = 0; xbi < start_xbp_all[i_a].size(); xbi++){
      start_xbp_all[i_a][xbi] = nullvector;
      //      start_xbs_all[i_a][xbi] = sqrt(-1.);
      //      start_xbsqrts_all[i_a][xbi] = sqrt(-1.);
      start_xbs_all[i_a][xbi] = 0.;
      start_xbsqrts_all[i_a][xbi] = 0.;
    }
    int max_particles;
    if (i_a == 0){max_particles = csi->n_particle + 2 + 1;}
    else {max_particles = csi->n_particle + 2;}
    for (int xbi = 1; xbi < max_particles; xbi++){
      start_xbs_all[i_a][intpow(2, xbi - 1)] = dx_s_RA[i_a][xbi][0];
      start_xbsqrts_all[i_a][intpow(2, xbi - 1)] = dx_sqrts_RA[i_a][xbi][0];
      /*
      start_xbs_all[i_a][intpow(2, xbi - 1)] = p_2[xbi];
      start_xbsqrts_all[i_a][intpow(2, xbi - 1)] = sqp_2[xbi];
      */
    }
  }

  if (coll_choice == 0){
    logger << LOG_DEBUG << "coll_choice == 0" << endl;
    fourvector p1(E, 0, 0, E), p2(E, 0, 0, -E);
    //    s_hat = 4. * pow(E, 2);
    //    sqrt_s_hat = 2. * E;
    for (int i = 0; i < 3; i++){x_pdf[i] = 1.;}

    for (int i_a = 0; i_a < n_dipoles; i_a++){
      logger << LOG_DEBUG << i_a << endl;
      start_xbp_all[i_a][0] = p1 + p2;
      start_xbp_all[i_a][1] = p1;
      start_xbp_all[i_a][2] = p2;
      start_xbs_all[i_a][0] = start_xbp_all[0][0].m2();
      start_xbsqrts_all[i_a][0] = sqrt(start_xbs_all[i_a][0]);
      start_xbp_all[i_a][3] = start_xbp_all[i_a][0];
      start_xbs_all[i_a][3] = start_xbs_all[i_a][0];
      start_xbsqrts_all[i_a][3] = start_xbsqrts_all[i_a][0];
      //      for (int i_p = 0; i_p < 3; i_p++){start_xbp_all[i_a][i_p] = start_xbp_all[0][i_p];}
      if (i_a == 0){
	start_xbp_all[i_a][xb_max - 4] = start_xbp_all[i_a][0];
	start_xbs_all[i_a][xb_max - 4] = start_xbs_all[i_a][0];
	start_xbsqrts_all[i_a][xb_max - 4] = start_xbsqrts_all[i_a][0];
      }
      else {
	start_xbp_all[i_a][xb_max_dipoles - 4] = start_xbp_all[i_a][0];
	start_xbs_all[i_a][xb_max_dipoles - 4] = start_xbs_all[i_a][0];
	start_xbsqrts_all[i_a][xb_max_dipoles - 4] = start_xbsqrts_all[i_a][0];
      }
    }
    pdf_factor[0] = 1.;
    pdf_factor[1] = 1.;
    pdf_factor[2] = 0.;
    //    tau_channel = 0;
  }
  
  for (int i_a = 0; i_a < n_dipoles; i_a++){
    for (int i_x = 4; i_x < start_xbs_all[i_a].size(); i_x += 4){
      vector<int> temp_vb = vectorbinary_from_binary(i_x);
      logger << LOG_DEBUG << "start_xbsqrts_all[" << setw(2) << i_a << "][" << setw(3) << i_x << "] = " << left << setw(23) << setprecision(15) << start_xbsqrts_all[i_a][i_x] << "   start_xbs_all[" << setw(2) << i_a << "][" << setw(3) << i_x << "] = " << left << setw(23) << setprecision(15) << start_xbs_all[i_a][i_x] << endl;
      if (temp_vb.size() == 1 || start_xbs_all[i_a][i_x] != 0. || start_xbsqrts_all[i_a][i_x] != 0.){
	logger << LOG_DEBUG << "start_xbsqrts_all[" << setw(2) << i_a << "][" << setw(3) << i_x << "] = " << left << setw(23) << setprecision(15) << start_xbsqrts_all[i_a][i_x] << "   start_xbs_all[" << setw(2) << i_a << "][" << setw(3) << i_x << "] = " << left << setw(23) << setprecision(15) << start_xbs_all[i_a][i_x] << endl;
      }
    }
  }
 
  
  sqrtsmin_opt = start_xbsqrts_all; // = start_xbsqrts;
  smin_opt = start_xbs_all; // = start_xbs;

  logger << LOG_DEBUG << "finished" << endl;
}

void phasespace_set::initialization_phasespace_subprocess_optimization_dipole_RA(){
  //void phasespace_set::initialization_phasespace_subprocess_optimization_dipole_RA(int _n_dipoles, vector<dipole_set> _dipole){
  Logger logger("phasespace_set::initialization_phasespace_subprocess_dipole_optimization_RA");
  logger << LOG_DEBUG << "started" << endl;

  // temporary - check if needed elsewhere !!!
  //  vector<double> dipole_sinx_min(n_dipoles, 0.);

  for (int i_a = 0; i_a < n_dipoles; i_a++){
    for (int i_x = 4; i_x < sqrtsmin_opt[i_a].size(); i_x += 4){
      if (sqrtsmin_opt[i_a][i_x] == start_xbsqrts_all[i_a][i_x]){}
      else if (sqrtsmin_opt[i_a][i_x] < start_xbsqrts_all[i_a][i_x]){sqrtsmin_opt[i_a][i_x] = start_xbsqrts_all[i_a][i_x];} 
      else if (sqrtsmin_opt[i_a][i_x] > start_xbsqrts_all[i_a][i_x]){
	if (vectorbinary_from_binary(i_x).size() == 1){
	  logger << LOG_FATAL << "phase-space cut is in contradiction to on-shell conditions of intermediate particle " << i_x << "!" << endl;
	  exit(1);
	}
      }
    }
  }
  
  for (int i_a = 0; i_a < n_dipoles; i_a++){
    for (int i_x = 4; i_x < sqrtsmin_opt[i_a].size(); i_x += 4){
      if (vectorbinary_from_binary(i_x).size() > 1){
	//      double temp_sqrts = max_subset_from_binary(start_xbsqrts_all[i_a], i_x);
	double temp_sqrts_opt = max_subset_from_binary(sqrtsmin_opt[i_a], i_x);
	if (temp_sqrts_opt > sqrtsmin_opt[i_a][i_x]){sqrtsmin_opt[i_a][i_x] = temp_sqrts_opt;}
	smin_opt[i_a][i_x] = pow(sqrtsmin_opt[i_a][i_x], 2);
      }
    }
  }
  
  logger << LOG_DEBUG << "OPT   list min_opt[0]" << endl;

  
  
  ///  logger << LOG_INFO << "n_dipoles = " << n_dipoles << endl;
  for (int i_a = 0; i_a < n_dipoles; i_a++){
    ///    logger << LOG_INFO << "smin_opt[" << i_a << "].size() = " << smin_opt[i_a].size() << endl;
    for (int i = 4; i < smin_opt[i_a].size(); i += 4){
      ///      logger << LOG_INFO << "i = " << i << "   sqrtsmin_opt[i_a].size() = " << sqrtsmin_opt[i_a].size() << "   smin_opt[i_a].size() = " << smin_opt[i_a].size() << endl;
      stringstream temp_ss;
      temp_ss << "sqrtsmin_opt[" << i_a << "][" << setw(3) << i << "] = " << setw(15) << setprecision(8) << sqrtsmin_opt[i_a][i] << "   smin_opt[" << i_a << "][" << setw(3) << i << "] = " << setw(15) << setprecision(8) << smin_opt[i_a][i] << endl;
      /*
      temp_ss << "sqrtsmin_opt[" << setw(2) << i_a << "][" << setw(3) << i << "] = " << setw(15) << setprecision(8) << sqrtsmin_opt[i_a][i] << "   smin_opt[" << setw(2) << i_a << "][" << setw(3) << i << "] = " << setw(15) << setprecision(8) << smin_opt[i_a][i] << "   ---   " << endl;
      vector<int> temp_vb = vectorbinary_from_binary(i);
      if (temp_vb.size() > 1){
	for (int j = 0; j < temp_vb.size(); j++){
	  temp_ss << "[" << setw(2) << i_a << "][" << setw(3) << i << "] = " << setw(15) << setprecision(8) << sqrtsmin_opt[i_a][i] << "   ";
	}
      }
      temp_ss << endl;
      */
      ///      logger << LOG_INFO << "i = " << i << "   sqrtsmin_opt[i_a].size() = " << sqrtsmin_opt[i_a].size() << "   smin_opt[i_a].size() = " << smin_opt[i_a].size() << endl;
      if (sqrtsmin_opt[i_a][i] != 0.){logger << LOG_INFO << temp_ss.str();}
      else {logger << LOG_DEBUG << temp_ss.str();}
      ///       logger << LOG_INFO << "i = " << i << "   sqrtsmin_opt[i_a].size() = " << sqrtsmin_opt[i_a].size() << "   smin_opt[i_a].size() = " << smin_opt[i_a].size() << endl;
     //      logger << LOG_DEBUG << temp_ss.str();
    }
  }
  logger.newLine(LOG_DEBUG);
  
  
  for (int i_a = 0; i_a < n_dipoles; i_a++){
    for (int i_x = 4; i_x < smin_opt[i_a].size(); i_x += 4){
      vector<int> temp_vb = vectorbinary_from_binary(i_x);
      if (temp_vb.size() == 1 || start_xbs_all[i_a][i_x] != 0. || start_xbsqrts_all[i_a][i_x] != 0.){
	logger << LOG_DEBUG << "start_xbsqrts_all[" << setw(2) << i_a << "][" << setw(3) << i_x << "] = " << left << setw(23) << setprecision(15) << start_xbsqrts_all[i_a][i_x] << "   start_xbs_all[" << setw(2) << i_a << "][" << setw(3) << i_x << "] = " << left << setw(23) << setprecision(15) << start_xbs_all[i_a][i_x] << endl;
      }
    }
  }
  logger.newLine(LOG_DEBUG);
  
  vector<vector<vector<double> > > dx_s_RA(n_dipoles);
  vector<vector<vector<double> > > dx_sqrts_RA(n_dipoles);
  for (int i_a = 0; i_a < dx_s_RA.size(); i_a++){
    dx_s_RA[i_a].resize((*RA_dipole)[i_a].type_parton().size(), vector<double> (1));
    dx_sqrts_RA[i_a].resize((*RA_dipole)[i_a].type_parton().size(), vector<double> (1));
  }
  for (int i_a = 0; i_a < dx_s_RA.size(); i_a++){
    for (int i_b = 1; i_b < dx_s_RA[i_a].size(); i_b++){
      if (dx_s_RA[i_a][i_b].size() == 1){
	dx_s_RA[i_a][i_b][0] = M2[abs((*RA_dipole)[i_a].type_parton()[i_b])];
	dx_sqrts_RA[i_a][i_b][0] = M[abs((*RA_dipole)[i_a].type_parton()[i_b])];
      }
    }
  }

  
  dipole_sinx_min.resize(n_dipoles, 0.);
  //  dipole_x.resize(n_dipoles, 0.);

  for (int i_a = 0; i_a < dx_s_RA.size(); i_a++){
    //  cout << "i_a = " << i_a << "   dx_s_RA[" << i_a << "].size() = " << dx_s_RA[i_a].size() << endl;
    for (int ib = 3; ib < dx_s_RA[i_a].size(); ib++){
      dipole_sinx_min[i_a] += M[abs((*RA_dipole)[i_a].type_parton()[ib])];
      //    cout << "M[abs((*RA_dipole)[" << i_a << "].type_parton()[" << ib << "])] = " << M[abs((*RA_dipole)[i_a].type_parton()[ib])] << endl;
      //    cout << "(*RA_dipole)[" << i_a << "].type_parton()[" << ib << "] = " << (*RA_dipole)[i_a].type_parton()[ib] << endl;
    }
    dipole_sinx_min[i_a] = pow(dipole_sinx_min[i_a], 2);
    logger << LOG_DEBUG << "smin_opt (only from masses):   dipole_sinx_min[" << setw(2) << i_a << "] = " << dipole_sinx_min[i_a] << endl;
    dipole_sinx_min[i_a] = smin_opt[i_a][smin_opt[i_a].size() - 4];
    logger << LOG_DEBUG << "smin_opt (from smin_opt[" << i_a << "][" << smin_opt[i_a].size() - 4 << "]):    dipole_sinx_min[" << setw(2) << i_a << "] = " << dipole_sinx_min[i_a] << endl;
  }
    
  logger << LOG_DEBUG << "tau_0 = " << tau_0 << "   -> sqrt(s^_min) = " << 2 * sqrt(tau_0) * E << endl;
  
  if (tau_0 < smin_opt[0][smin_opt[0].size() - 4] / (4 * pow(E, 2))){tau_0 = smin_opt[0][smin_opt[0].size() - 4] / (4 * pow(E, 2));}
  
  logger << LOG_DEBUG << "tau_0 = " << tau_0 << "   -> sqrt(s^_min) = " << 2 * sqrt(tau_0) * E << endl;
  logger.newLine(LOG_DEBUG);

  // shifted elsewhere
  /*
  for (int j = 0; j < tau_MC_tau_gamma[0].size(); j++){
    double random_tau = (double(j + 1)) / double(tau_MC_tau_beta[0].size());
    tau_MC_tau_gamma[0][j] = h_propto_pot(random_tau, tau_0, exp_pdf);
    logger << LOG_DEBUG << "tau_MC_tau_gamma[0][" << setw(3) << j << "] = " << tau_MC_tau_gamma[0][j] << endl;
  }
  */  
  
  xbp_all = start_xbp_all;
  xbs_all = start_xbs_all;
  xbsqrts_all = start_xbsqrts_all;
  
  corrected_xbp_all = start_xbp_all;
  corrected_xbs_all = start_xbs_all;
  corrected_xbsqrts_all = start_xbsqrts_all;

  logger << LOG_DEBUG << "finished" << endl;
}



/*
void phasespace_set::initialization_phasespace_MC_tau(){
  Logger logger("phasespace_set::initialization_phasespace_MC_tau");
  logger << LOG_DEBUG << "started" << endl;
  logger << LOG_INFO << "2 * sqrt(tau_0) * E = " << 2 * sqrt(tau_0) * E << endl;

  logger << LOG_INFO << "tau_MC_map.size() = " << tau_MC_map.size() << endl;
  for (int i_m = tau_MC_map.size() - 1; i_m > 0; i_m--){
    logger << LOG_INFO << "tau_MC_map[" << i_m << "] = " << tau_MC_map[i_m] << endl;
    logger << LOG_INFO << "psi_M[" << abs(tau_MC_map[i_m]) << "] = " << M[abs(tau_MC_map[i_m])] << endl;
    logger << LOG_INFO << "psi_Gamma[" << abs(tau_MC_map[i_m]) << "] = " << Gamma[abs(tau_MC_map[i_m])] << endl;
    logger << LOG_INFO << "psi_map_Gamma[" << abs(tau_MC_map[i_m]) << "] = " << map_Gamma[abs(tau_MC_map[i_m])] << endl;
    if (map_Gamma[abs(tau_MC_map[i_m])] == 0.){tau_MC_map.erase(tau_MC_map.begin() + i_m);}
    else if (2 * sqrt(tau_0) * E > M[abs(tau_MC_map[i_m])]){tau_MC_map.erase(tau_MC_map.begin() + i_m);}
    // leads to problems if M_res < smin !!! needs to be investigated !!!
    //    else if (2 * sqrt(tau_0) * E > M[abs(tau_MC_map[i_m])] + 5 * map_Gamma[abs(tau_MC_map[i_m])]){tau_MC_map.erase(tau_MC_map.begin() + i_m);}
  }

  n_tau_MC_channel = tau_MC_map.size();
  tau_MC_tau_gamma.resize(n_tau_MC_channel, vector<double> (n_tau_bins));

  if (switch_MC_tau == -1){tau_MC_map.resize(1);}

  // !!! input !!!
  int switch_minimum_weight = 1;
  double limit_minimum_weight = 0.1;
  double reserved_minimum_weight = 0.1;
  //  double limit_minimum_weight = 1.e-4;
  //  double reserved_minimum_weight = 1.e-3;

  // MC_tau is initialized twice !!!
  MC_tau = multichannel_set("MC_tau", tau_MC_map.size(), n_alpha_events, n_alpha_steps, switch_minimum_weight, limit_minimum_weight, reserved_minimum_weight, switch_MC_tau, filename_MCweight);

  logger << LOG_DEBUG << "MC_tau.n_channel = " << MC_tau.n_channel << endl;
  for (int i_c = 0; i_c < MC_tau.n_channel; i_c++){
    logger << LOG_DEBUG << "tau_MC_map  [" << i_c << "] = " << tau_MC_map[i_c] << endl;
    logger << LOG_DEBUG << "MC_tau.alpha[" << i_c << "] = " << MC_tau.alpha[i_c] << endl;
    logger << LOG_DEBUG << "MC_tau.beta [" << i_c << "] = " << MC_tau.beta[i_c] << endl;
  }
  logger << LOG_DEBUG << "finished" << endl;
}
*/



void phasespace_set::initialization_phasespace_MC_tau(){
  Logger logger("phasespace_set::initialization_phasespace_MC_tau");
  logger << LOG_DEBUG << "started" << endl;

  // only non-output statements in initialization_phasespace_MC_tau (which are not repeated here):
  logger << LOG_DEBUG << "tau_MC_map.size() = " << tau_MC_map.size() << endl;
  for (int i_m = tau_MC_map.size() - 1; i_m > 0; i_m--){
    if (map_Gamma[abs(tau_MC_map[i_m])] == 0.){tau_MC_map.erase(tau_MC_map.begin() + i_m);}
    else if (2 * sqrt(tau_0) * E > M[abs(tau_MC_map[i_m])]){tau_MC_map.erase(tau_MC_map.begin() + i_m);}
  }
  logger << LOG_DEBUG << "tau_MC_map.size() = " << tau_MC_map.size() << endl;
  // initialization of initial-state multichannel (MC) related tau parameters
  // to be shifted elsewhere
  n_tau_MC_channel = tau_MC_map.size();
  tau_MC_tau_gamma.resize(n_tau_MC_channel, vector<double> (n_tau_bins));
  
  logger << LOG_DEBUG << "tau_MC_tau_gamma.size() = " << tau_MC_tau_gamma.size() << endl;
  for (int j = 0; j < tau_MC_tau_gamma[0].size(); j++){
    double random_tau = (double(j + 1)) / double(n_tau_bins);
    tau_MC_tau_gamma[0][j] = h_propto_pot(random_tau, tau_0, exp_pdf);
  }
  
  for (int j = 0; j < tau_MC_tau_gamma[0].size(); j++){
    logger << LOG_DEBUG << "tau_MC_tau_gamma[0][" << j << "]" << " = " << setprecision(20) << setw(28) << tau_MC_tau_gamma[0][j] << "   " << double2hexastr(tau_MC_tau_gamma[0][j]) << endl;
  }

  // to be shifted elsewhere - used so far somewhere...
  if (switch_MC_tau == -1){tau_MC_map.resize(1);}

  logger << LOG_DEBUG << "finished" << endl;
}

void phasespace_set::initialization_phasespace_MC_x_dipole_RA(){
  Logger logger("phasespace_set::initialization_phasespace_MC_x_dipole_RA");
  logger << LOG_DEBUG << "started" << endl;

  // determine tau_MC_map from MC_x_dipole_mapping (including 0!)
  for (int i_a = 0; i_a < n_dipoles; i_a++){

    for (int i_m = MC_x_dipole_mapping[i_a].size() - 1; i_m > 0; i_m--){
      logger << LOG_INFO << "MC_x_dipole_mapping[" << i_a << "][" << i_m << "] = " << MC_x_dipole_mapping[i_a][i_m] << endl;
      logger << LOG_INFO << "psi_M[" << abs(MC_x_dipole_mapping[i_a][i_m]) << "] = " << M[abs(MC_x_dipole_mapping[i_a][i_m])] << endl;
      logger << LOG_INFO << "psi_Gamma[" << abs(MC_x_dipole_mapping[i_a][i_m]) << "] = " << Gamma[abs(MC_x_dipole_mapping[i_a][i_m])] << endl;
      logger << LOG_INFO << "psi_map_Gamma[" << abs(MC_x_dipole_mapping[i_a][i_m]) << "] = " << map_Gamma[abs(MC_x_dipole_mapping[i_a][i_m])] << endl;
      if (map_Gamma[abs(MC_x_dipole_mapping[i_a][i_m])] == 0.){MC_x_dipole_mapping[i_a].erase(MC_x_dipole_mapping[i_a].begin() + i_m);}
      else if (2 * sqrt(tau_0) * E > M[abs(MC_x_dipole_mapping[i_a][i_m])]){MC_x_dipole_mapping[i_a].erase(MC_x_dipole_mapping[i_a].begin() + i_m);}
      // leads to problems if M_res < smin !!! needs to be investigated !!!
      //      else if (2 * sqrt(tau_0) * E > M[abs(MC_x_dipole_mapping[i_a][i_m])] + 5 * map_Gamma[abs(MC_x_dipole_mapping[i_a][i_m])]){MC_x_dipole_mapping[i_a].erase(MC_x_dipole_mapping[i_a].begin() + i_m);}
    }

    // to add MC_x_dipole mappings also to MC_tau
    for (int i_m = 0; i_m < MC_x_dipole_mapping[i_a].size(); i_m++){
      int flag = tau_MC_map.size();
      for (int j_m = 0; j_m < tau_MC_map.size(); j_m++){
	logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << "   i_m = " << i_m << "   j_m = " << j_m << endl;
 	if (MC_x_dipole_mapping[i_a][i_m] == tau_MC_map[j_m]){flag = j_m; break;}
      }
      if (flag == tau_MC_map.size()){tau_MC_map.push_back(MC_x_dipole_mapping[i_a][i_m]);}
    }
  }

  //  MC_x_dipole_mapping[0] = tau_MC_map;

  // output of tau_MC_map and MC_x_dipole_mapping
  stringstream temp_ss;
  temp_ss << "tau_MC_map = ";
  for (int i_m = 0; i_m < tau_MC_map.size(); i_m++){temp_ss << setw(4) << tau_MC_map[i_m];}
  logger << LOG_DEBUG << temp_ss.str() << endl;
  for (int i_a = 0; i_a < n_dipoles; i_a++){
    stringstream temp_ss;
    temp_ss << "MC_x_dipole_mapping[" << setw(2) << i_a << "] = ";
    for (int i_m = 0; i_m < MC_x_dipole_mapping[i_a].size(); i_m++){temp_ss << setw(4) << MC_x_dipole_mapping[i_a][i_m];}
    logger << LOG_DEBUG << temp_ss.str() << endl;
  }

  int switch_minimum_weight = 1;
  double limit_minimum_weight = 0.1;
  double reserved_minimum_weight = 0.1;
  //  double limit_minimum_weight = 0.01;
  //  double reserved_minimum_weight = 0.01;
  MC_x_dipole.resize(n_dipoles);
  for (int i_a = 1; i_a < n_dipoles; i_a++){
    if (switch_MC_x_dipole == -1){MC_x_dipole_mapping[i_a].resize(1);}

    stringstream temp_ss;
    temp_ss << "MC_x_dipole_" << i_a;

    logger << LOG_DEBUG << "temp_ss.str() = " << temp_ss.str() << endl;
    logger << LOG_DEBUG << "MC_x_dipole_mapping[" << i_a << "].size() = " << MC_x_dipole_mapping[i_a].size() << endl;
    logger << LOG_DEBUG << "n_alpha_events = " << n_alpha_events << endl;
    logger << LOG_DEBUG << "n_alpha_steps = " << n_alpha_steps << endl;
    logger << LOG_DEBUG << "switch_minimum_weight = " << switch_minimum_weight << endl;
    logger << LOG_DEBUG << "limit_minimum_weight = " << limit_minimum_weight << endl;
    logger << LOG_DEBUG << "reserved_minimum_weight = " << reserved_minimum_weight << endl;
    logger << LOG_DEBUG << "switch_MC_x_dipole = " << switch_MC_x_dipole << endl;
    logger << LOG_DEBUG << "filename_MCweight = " << filename_MCweight << endl;
    logger << LOG_DEBUG << "filename_MCweight_in_contribution = " << filename_MCweight_in_contribution << endl;
     
    MC_x_dipole[i_a] = multichannel_set(temp_ss.str(), MC_x_dipole_mapping[i_a].size(), n_alpha_events, n_alpha_steps, switch_minimum_weight, limit_minimum_weight, reserved_minimum_weight, switch_MC_x_dipole, filename_MCweight, filename_MCweight_in_contribution);
  }

  logger << LOG_DEBUG << "finished" << endl;
}

void phasespace_set::initialization_phasespace_IS(){
  Logger logger("phasespace_set::initialization_phasespace_IS");
  logger << LOG_DEBUG << "started" << endl;

  g_z_coll.resize(3); // could be deleted later... !!!

  z1z2_opt_end = -1; // not needed here - only in CA (or qTsubtraction)

  // should become input !!!
  int switch_minimum_weight = 1;
  double limit_minimum_weight = 0.1;
  double reserved_minimum_weight = 0.1;
  //  double limit_minimum_weight = 1.e-4;
  //  double reserved_minimum_weight = 1.e-3;


  /*
  // shifted to initialization_phasespace_RA !!!
  if (csi->class_contribution_CS_real){
    if (!switch_RS_mapping){
      for (int i_a = 1; i_a < n_dipoles; i_a++){generic->ac_tau_psp_dipole(i_a, MC_x_dipole_mapping[i_a], *this);}
      if (csi->type_contribution == "RRA" ||
	  csi->type_contribution == "RRJ"){
	initialization_fake_dipole_mapping_RRA();
	//	initialization_fake_dipole_mapping_RRA(dipole, *this, generic);
      }
      initialization_phasespace_MC_x_dipole_RA();
    }
    else if (switch_RS_mapping){
      tau_MC_map = MC_x_dipole_mapping[0];
    }
  }
  */
  
  /*
  // only non-output statements in initialization_phasespace_MC_tau (which are not repeated here):
  // shifted back to new initialization_phasespace_MC_tau function
  for (int i_m = tau_MC_map.size() - 1; i_m > 0; i_m--){
    if (map_Gamma[abs(tau_MC_map[i_m])] == 0.){tau_MC_map.erase(tau_MC_map.begin() + i_m);}
    else if (2 * sqrt(tau_0) * E > M[abs(tau_MC_map[i_m])]){tau_MC_map.erase(tau_MC_map.begin() + i_m);}
  }
  // initialization of initial-state multichannel (MC) related tau parameters
  // to be shifted elsewhere
  n_tau_MC_channel = tau_MC_map.size();
  tau_MC_tau_gamma.resize(n_tau_MC_channel, vector<double> (n_tau_bins));
  
  for (int j = 0; j < tau_MC_tau_gamma[0].size(); j++){
    double random_tau = (double(j + 1)) / double(n_tau_bins);
    tau_MC_tau_gamma[0][j] = h_propto_pot(random_tau, tau_0, exp_pdf);
  }

  // to be shifted elsewhere - used so far somewhere...
  if (switch_MC_tau == -1){tau_MC_map.resize(1);}
  */
    
  initialization_phasespace_MC_tau();
  
  // MC_tau is initialized twice !!! not any more...
  MC_tau = multichannel_set("MC_tau", tau_MC_map.size(), n_alpha_events, n_alpha_steps, switch_minimum_weight, limit_minimum_weight, reserved_minimum_weight, switch_MC_tau, filename_MCweight, filename_MCweight_in_contribution);

  // Check if tau is fixed, and if so, deactivate "IS_tau":
  if (start_xbs_all[0][xb_max - 4] != 0.){switch_IS_tau = -1;}
  // meaning undefined !!! switch_minimum_weight, limit_minimum_weight, reserved_minimum_weight
  IS_tau = importancesampling_set("IS_tau", n_tau_bins, n_tau_steps, n_tau_events, switch_minimum_weight, limit_minimum_weight, reserved_minimum_weight, switch_IS_tau, filename_tauweight, filename_tauweight_in_contribution);

  // meaning undefined !!! switch_minimum_weight, limit_minimum_weight, reserved_minimum_weight
  IS_x1x2 = importancesampling_set("IS_x1x2", n_x1x2_bins, n_x1x2_steps, n_x1x2_events, switch_minimum_weight, limit_minimum_weight, reserved_minimum_weight, switch_IS_x1x2, filename_x1x2weight, filename_x1x2weight_in_contribution);


  if (csi->class_contribution_CS_collinear){
    IS_z1z2.resize(3);
    for (int i_c = 1; i_c < 3; i_c++){
      stringstream temp_name_ss;
      temp_name_ss << "IS_z1z2[" << i_c << "]"; 
      string temp_name = temp_name_ss.str();
      IS_z1z2[i_c] = importancesampling_set(temp_name, n_z1z2_bins, n_z1z2_steps, n_z1z2_events, switch_minimum_weight, limit_minimum_weight, reserved_minimum_weight, switch_IS_z1z2, filename_z1z2weight[i_c], filename_z1z2weight_in_contribution[i_c]);
    }
    
    z_coll.resize(3);
    all_xz_coll_pdf.resize(3, vector<double> (3));
    g_z_coll.resize(3);
  }

  // determination of global variable 'active_optimization':
  // active_optimization:
  // active_optimization = -1: no optimization of any mapping (very unlikely)
  // active_optimization =  0: no active optimization (like in standard runs)
  // active_optimization =  1: at least one active optimization (like in grid runs)
  active_optimization = -1;
  if (MC_phasespace.active_optimization != -1){if (MC_phasespace.active_optimization > active_optimization){active_optimization = MC_phasespace.active_optimization;}}
  logger << LOG_DEBUG << "MC_phasespace.active_optimization = " << MC_phasespace.active_optimization << endl;
  if (MC_tau.active_optimization != -1){if (MC_tau.active_optimization > active_optimization){active_optimization = MC_tau.active_optimization;}}
  logger << LOG_DEBUG << "MC_tau.active_optimization = " << MC_tau.active_optimization << endl;
  if (csi->class_contribution_CS_real){
    for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
      //    for (int i_a = 1; i_a < n_dipoles; i_a++){
      logger << LOG_DEBUG << "MC_x_dipole[" << i_a << "].active_optimization = " << MC_x_dipole[i_a].active_optimization << endl;
      if (MC_x_dipole[i_a].active_optimization != -1){if (MC_x_dipole[i_a].active_optimization > active_optimization){active_optimization = MC_x_dipole[i_a].active_optimization;}}
    }
  }
  // IS_phasespace is treated in randommanager !!!
  logger << LOG_DEBUG << "IS_tau.active_optimization = " << IS_tau.active_optimization << endl;
  if (IS_tau.active_optimization != -1){if (IS_tau.active_optimization > active_optimization){active_optimization = IS_tau.active_optimization;}}
  logger << LOG_DEBUG << "IS_x1x2.active_optimization = " << IS_x1x2.active_optimization << endl;
  if (IS_x1x2.active_optimization != -1){if (IS_x1x2.active_optimization > active_optimization){active_optimization = IS_x1x2.active_optimization;}}
  if (csi->class_contribution_CS_collinear){
    for (int i_c = 1; i_c < 3; i_c++){
      logger << LOG_DEBUG << "IS_z1z2[" << i_c << "].active_optimization = " << IS_z1z2[i_c].active_optimization << endl;
      if (IS_z1z2[i_c].active_optimization != -1){if (IS_z1z2[i_c].active_optimization > active_optimization){active_optimization = IS_z1z2[i_c].active_optimization;}}
    }
  }  
  
  // determination of global variable 'end_optimization':


  logger << LOG_DEBUG << "finished" << endl;
}



// !!! could later be replaced by CX version (once collinear has been merged to multicollinear (ncollinear) !!!
void phasespace_set::initialization_phasespace_IS_CA(){
  Logger logger("phasespace_set::initialization_IS_CA");
  logger << LOG_DEBUG << "started" << endl;

  /*
  int switch_minimum_weight = 1;
  double limit_minimum_weight = 0.1;
  double reserved_minimum_weight = 0.1;
  IS_z1z2.resize(3);
  for (int ic = 1; ic < 3; ic++){
    stringstream temp_name_ss;
    temp_name_ss << "IS_z1z2[" << ic << "]"; 
    string temp_name = temp_name_ss.str();
    IS_z1z2[ic] = importancesampling_set(temp_name, n_z1z2_bins, n_z1z2_steps, n_z1z2_events, switch_minimum_weight, limit_minimum_weight, reserved_minimum_weight, switch_IS_z1z2, filename_z1z2weight[ic]);
  }

  z_coll.resize(3);
  all_xz_coll_pdf.resize(3, vector<double> (3));
  g_z_coll.resize(3);
  */
    
  logger << LOG_DEBUG << "finished" << endl;
}

void phasespace_set::initialization_phasespace_IS_CX(){
  Logger logger("phasespace_set::initialization_phasespace_IS_CX");
  logger << LOG_DEBUG << "started" << endl;


  QT_eps = 0;
  QT_random_IS = 0.;

  QT_random_z.resize(3);
  for (int i_z = 1; i_z < 3; i_z++){
    stringstream temp;
    temp << "z[" << i_z << "]";

    logger << LOG_INFO << "n_z1z2_events = " << n_z1z2_events << endl;
    logger << LOG_INFO << "switch_IS_z1z2 = " << switch_IS_z1z2 << endl;
    logger << LOG_INFO << "n_z1z2_steps = " << n_z1z2_steps << endl;
    logger << LOG_INFO << "n_z1z2_bins = " << n_z1z2_bins << endl;
    logger << LOG_INFO << "temp.str() = " << temp.str() << endl;
    
    QT_random_z[i_z] = new randomvariable(n_z1z2_events, switch_IS_z1z2, n_z1z2_steps, n_z1z2_bins, temp.str());
    ///    random_manager.register_variable(QT_random_z[i_z], true);
    random_manager.register_variable(QT_random_z[i_z]);
  }

  z_coll.resize(3);
  g_z_coll.resize(3);

  logger << LOG_DEBUG << "finished" << endl;
}

void phasespace_set::initialization_phasespace_IS_QT(){
  Logger logger("phasespace_set::initialization_phasespace_IS_QT");
  logger << LOG_DEBUG << "started" << endl;

  QT_eps = 0;
  QT_random_IS = 0.;

  QT_qt2 = 0.;
  QT_jacqt2 = 0.;
  QT_random_qt2 = 0.;

  QT_random_qt = new randomvariable(n_z1z2_events, switch_IS_z1z2, n_z1z2_steps, n_z1z2_bins, "qt2");
  ///  random_manager.register_variable(QT_random_qt, true);
  random_manager.register_variable(QT_random_qt);

  logger << LOG_DEBUG << "finished" << endl;
}




void phasespace_set::initialization_MC(){
  Logger logger("phasespace_set::initialization_MC");
  logger << LOG_DEBUG << "started" << endl;

  initialization_phasespace_IS();
  
  g_tot = 0.;
  g_MC = 0.;
  if (coll_choice){
    g_pdf = 0.;
    g_tau = 0.;
    g_x1x2 = 0.;
  }
  else {
    g_pdf = 1.;
    g_tau = 1.;
    g_x1x2 = 1.;
  }
  //  g_global_NWA = 1.;

  MC_g_IS_global = 1.;
  
  MC_sum_w_channel.resize(MC_n_channel, 0.);
  MC_g_IS_channel.resize(MC_n_channel, 0.);




  logger << LOG_DEBUG << "switch_IS_mode_phasespace = " << switch_IS_mode_phasespace << endl;
  if (switch_IS_mode_phasespace == 1 || switch_IS_mode_phasespace == 2){
    // create random variables for phase space integration with grids for each variable in each channel (MC_n_channel * (3n - 4))
    phasespace_randoms.resize(no_random * MC_n_channel);
    string tmp_name;
    ostringstream convert;
    for (int i = 0; i < no_random * MC_n_channel; i++) {
      convert.str(string());
      convert << "r" << i;
      tmp_name = convert.str();
      logger << LOG_DEBUG << "initialize randomvariable " << i << endl; 
      phasespace_randoms[i] = new randomvariable(n_IS_events, switch_IS_MC, n_IS_steps, n_IS_gridsize, tmp_name);
      ///      random_manager.register_variable(phasespace_randoms[i], true);
      random_manager.register_variable(phasespace_randoms[i]);
    }
    container_IS_switch.resize(container_IS_name.size(), 0);
    //    psi.phasespace_randoms=&phasespace_randoms;
  }
  else if (switch_IS_mode_phasespace == 3 || switch_IS_mode_phasespace == 4){
    // create random variables for phase space integration with one grid only for each mapping (propagator, t channel, etc.)
    phasespace_randoms.resize(container_IS_name.size());
    // later: only one file for all grids
    string tmp_name;
    ostringstream convert;
    logger << LOG_DEBUG_VERBOSE << "phasespace_randoms.size() = " << phasespace_randoms.size() << endl; 
    for (int i = 0; i < phasespace_randoms.size(); i++) {
      convert.str(string());
      convert << "yr" << i;
      //    tmp_name = convert.str();
      tmp_name = convert.str() + "_" + container_IS_name[i];
      logger << LOG_DEBUG_VERBOSE << "initialize randomvariable " << i << endl; 
      phasespace_randoms[i] = new randomvariable(n_IS_events, switch_IS_MC, n_IS_steps, n_IS_gridsize, tmp_name);

      //    phasespace_randoms[i] = new randomvariable(n_events_factor*n_z1z2_events,n_z1z2_steps,n_r_bins,tmp_name);
      logger << LOG_DEBUG_VERBOSE << container_IS_name[i] << "   n_IS_events = " << n_IS_events << "   n_IS_steps = " << n_IS_steps << "   n_IS_gridsize = " << n_IS_gridsize << endl;
      ///      random_manager.register_variable(phasespace_randoms[i], true);
      random_manager.register_variable(phasespace_randoms[i]);
    }
    container_IS_switch.resize(container_IS_name.size(), 1);
    //    psi.phasespace_randoms=&phasespace_randoms;
  }
  else {
    container_IS_switch.resize(container_IS_name.size(), 0);
  }
  logger << LOG_DEBUG << "definitions of integration variables" << endl;


  /// check if new implementation via multichannel_set::readin_MCweight_optimization() reproduces old result:
  /*  
  vector<double> temp_MC_phasespace_alpha = MC_phasespace.alpha;
  vector<double> temp_MC_tau_alpha = MC_tau.alpha;
  vector<vector<double> > temp_MC_x_dipole_alpha(MC_x_dipole.size());
  if (!csi->class_contribution_CS_real){
    for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
      temp_MC_x_dipole_alpha[i_a] = MC_x_dipole[i_a].alpha;
    }
  }
  */

  /*
  vector<double> temp_IS_x1x2_alpha = IS_x1x2.alpha;
  vector<double> temp_IS_tau_alpha = IS_tau.alpha;
  vector<double> temp_IS_x1x2_beta = IS_x1x2.beta;
  vector<double> temp_IS_tau_beta = IS_tau.beta;
  */

   
  // switch_IS_tau
  // switch_IS_tau = 0: No IS optimization of partonic CMS energy.
  // switch_IS_tau = 1: Basic tauweights are taken; optimization is done according to parameters ...
  // switch_IS_tau = 2: tauweights read in from 'filename_tauweight_in_contribution'; no further optimization
  // switch_IS_tau = 3: tauweights read in from 'filename_tauweight_in_contribution'; further optimization according to parameters ...
  if (switch_IS_tau == 2 || switch_IS_tau == 3){
    logger << LOG_INFO << "filename_tauweight_in_contribution = " << filename_tauweight_in_contribution << endl;
    IS_tau.readin_IS_optimization();
  }

  // switch_IS_x1x2
  // switch_IS_x1x2 = 0: No IS optimization of x1x2 mappings.
  // switch_IS_x1x2 = 1: Basic x1x2weights are taken; optimization is done according to parameters ...
  // switch_IS_x1x2 = 2: x1x2weights read in from 'filename_x1x2weight_in_contribution'; no further optimization
  // switch_IS_x1x2 = 3: x1x2weights read in from 'filename_x1x2weight_in_contribution'; further optimization according to parameters ...
  if (switch_IS_x1x2 == 2 || switch_IS_x1x2 == 3){
    logger << LOG_INFO << "filename_x1x2weight_in_contribution = " << filename_x1x2weight_in_contribution << endl;
    
    IS_x1x2.readin_IS_optimization();
  }
  logger << LOG_DEBUG << "initialization of x1x2weight optimization finished" << endl;
  


  /*
  /// check if new implementation via multichannel_set::readin_MCweight_optimization() reproduces old result:
  if (switch_IS_tau == 2 || switch_IS_tau == 3){
    if (temp_IS_tau_alpha != IS_tau.alpha || temp_IS_tau_beta != IS_tau.beta){logger << LOG_FATAL << "Input of IS_tau incorrect." << endl; exit(1);}
  }

  if (switch_IS_x1x2 == 2 || switch_IS_x1x2 == 3){
    if (temp_IS_x1x2_alpha != IS_x1x2.alpha || temp_IS_x1x2_beta != IS_x1x2.beta){logger << LOG_FATAL << "Input of IS_x1x2 incorrect." << endl; exit(1);}
  }
  // end
  */

  logger << LOG_DEBUG << "../include/initialization.MCweight.optimization.cxx   ended" << endl;

 
  logger << LOG_DEBUG << "finished" << endl;
}



//vector<dipole_set> & dipole, phasespace_set & psi, call_generic & generic
void phasespace_set::initialization_fake_dipole_mapping_RRA(){
  Logger logger("phasespace_set::initialization_fake_dipole_mapping_RRA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  int n_dipoles = dipole.size();
  logger << LOG_DEBUG << "fake_dipole to improve initial-state mappings started..." << endl;

  //    logger << LOG_DEBUG << "n_dipoles = " << n_dipoles << endl;
  logger << LOG_DEBUG << "REAL *RA_DIPOLE: (*RA_dipole).size() = " << (*RA_dipole).size() << endl;

    vector<vector<dipole_set> > fake_dipole((*RA_dipole).size());
    vector<vector<dipole_set> > fake_dipole_candidate((*RA_dipole).size());
    for (int i_a = 1; i_a < (*RA_dipole).size(); i_a++){
      phasespace_set fake_psi;// = psi;
      fake_psi.no_map.resize(1, 0);
      fake_psi.no_prc.resize(1, 0);
      fake_psi.o_map.resize(1, vector<int> (o_map.size() - 1));
      fake_psi.o_prc.resize(1, vector<int> (o_prc.size() - 1));
      fake_psi.phasespace_order_alpha_s.resize(1, phasespace_order_alpha_s[i_a]);
      fake_psi.phasespace_order_alpha_e.resize(1, phasespace_order_alpha_e[i_a]);
      fake_psi.phasespace_order_interference.resize(1, phasespace_order_interference[i_a]);
      fake_psi.MC_n_channel_phasespace.resize(1, MC_n_channel_phasespace[i_a]);
      fake_psi.MC_sum_channel_phasespace.resize(1, MC_sum_channel_phasespace[i_a]);

      fake_psi.M = M;
      /*
      fake_psi.no_map.erase(fake_psi.no_map.begin() + 1, fake_psi.no_map.end());
      fake_psi.no_prc.erase(fake_psi.no_map.begin() + 1, fake_psi.no_prc.end());
      fake_psi.o_map.erase(fake_psi.o_map.begin() + 1, fake_psi.o_map.end());
      fake_psi.o_prc.erase(fake_psi.o_prc.begin() + 1, fake_psi.o_prc.end());
      fake_psi.phasespace_order_alpha_s.erase(fake_psi.phasespace_order_alpha_s.begin() + 1, fake_psi.phasespace_order_alpha_s.end());
      fake_psi.phasespace_order_alpha_e.erase(fake_psi.phasespace_order_alpha_e.begin() + 1, fake_psi.phasespace_order_alpha_e.end());
      fake_psi.phasespace_order_interference.erase(fake_psi.phasespace_order_interference.begin() + 1, fake_psi.phasespace_order_interference.end());
      fake_psi.phasespace_order_alpha_s[0] = phasespace_order_alpha_s[i_a];
      fake_psi.phasespace_order_alpha_e[0] = phasespace_order_alpha_e[i_a];
      fake_psi.phasespace_order_interference[0] = phasespace_order_interference[i_a];
*/
      /*
      fake_psi.MC_n_channel_phasespace.erase(fake_psi.MC_n_channel_phasespace.begin() + 1, fake_psi.MC_n_channel_phasespace.end());
      fake_psi.MC_sum_channel_phasespace.erase(fake_psi.MC_sum_channel_phasespace.begin() + 1, fake_psi.MC_sum_channel_phasespace.end());
      */
      logger << LOG_DEBUG_VERBOSE << "phasespace_order_alpha_s[" << i_a << "] - 1 = " << phasespace_order_alpha_s[i_a] - 1 << endl;
      logger << LOG_DEBUG_VERBOSE << "phasespace_order_alpha_s[" << i_a << "]     = " << phasespace_order_alpha_e[i_a] << endl;
      logger << LOG_DEBUG_VERBOSE << "fake_psi.phasespace_order_alpha_s[" << 0 << "] = " << fake_psi.phasespace_order_alpha_s[0] << endl;
      logger << LOG_DEBUG_VERBOSE << "fake_psi.phasespace_order_alpha_s[" << 0 << "] = " << fake_psi.phasespace_order_alpha_e[0] << endl;
      logger << LOG_DEBUG_VERBOSE << "fake_psi.MC_n_channel_phasespace.size() = " << fake_psi.MC_n_channel_phasespace.size() << endl;
      logger << LOG_DEBUG_VERBOSE << "fake_psi.MC_sum_channel_phasespace.size() = " << fake_psi.MC_sum_channel_phasespace.size() << endl;


      vector<int> fake_type_parton = (*RA_dipole)[i_a].type_parton();
      vector<int> fake_basic_type_parton = (*RA_dipole)[i_a].basic_type_parton();
      fake_dipole[i_a].push_back(dipole_set((*RA_dipole)[i_a].name(), fake_type_parton, fake_basic_type_parton, symmetry_factor, no_map[i_a], o_map[i_a], no_prc[i_a], o_prc[i_a], 0));
      QCD_determine_dipoles(fake_dipole_candidate[i_a], fake_type_parton, fake_basic_type_parton);
      logger << LOG_DEBUG_VERBOSE << "REAL *RA_DIPOLE " << i_a << endl;



      QCD_selection_fake_dipoles(fake_dipole[i_a], fake_dipole_candidate[i_a], fake_psi.phasespace_order_alpha_s[0] - 1, fake_psi.phasespace_order_alpha_e[0], phasespace_order_interference[i_a], fake_psi.RA_singular_region, fake_psi.RA_singular_region_name, fake_psi.RA_singular_region_list, fake_psi, *generic, generic->determination_no_subprocess_doubledipole, generic->determination_MCchannels_doubledipole);



      //      QCD_selection_fake_dipoles(fake_dipole[i_a], fake_dipole_candidate[i_a], phasespace_order_alpha_s[i_a] - 1, phasespace_order_alpha_e[i_a], phasespace_order_interference[i_a], fake_psi.RA_singular_region, fake_psi.RA_singular_region_name, fake_psi.RA_singular_region_list, fake_psi, generic, generic->determination_no_subprocess_doubledipole, generic->determination_MCchannels_doubledipole);

      vector<vector<int> > fake_MC_x_dipole_mapping(fake_dipole[i_a].size());
      logger << LOG_DEBUG_VERBOSE << "fake_dipole[" << i_a << "].size() = " << fake_dipole[i_a].size() << endl;
      for (int j_a = 1; j_a < fake_dipole[i_a].size(); j_a++){
	logger << LOG_DEBUG_VERBOSE << "before: fake_MC_x_dipole_mapping[" << j_a << "].size() = " << fake_MC_x_dipole_mapping[j_a].size() << endl;
	generic->ac_tau_psp_doubledipole(j_a, fake_MC_x_dipole_mapping[j_a], fake_psi);
	logger << LOG_DEBUG_VERBOSE << "after:  fake_MC_x_dipole_mapping[" << j_a << "].size() = " << fake_MC_x_dipole_mapping[j_a].size() << endl;
      }

      logger << LOG_DEBUG_VERBOSE << "fake_dipole[" << i_a << "].size() = " << fake_dipole[i_a].size() << endl;
      for (int j_a = 1; j_a < fake_dipole[i_a].size(); j_a++){
	logger << LOG_DEBUG_VERBOSE << "fake_MC_x_dipole_mapping[" << j_a << "].size() = " << fake_MC_x_dipole_mapping[j_a].size() << endl;
	for (int i_m = 0; i_m < fake_MC_x_dipole_mapping[j_a].size(); i_m++){
	  int flag = MC_x_dipole_mapping[i_a].size();
	  for (int j_m = 0; j_m < MC_x_dipole_mapping[i_a].size(); j_m++){
	    logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << "   i_m = " << i_m << "   j_m = " << j_m << "   fake_MC_x_dipole_mapping[" << j_a << "][" << i_m << "] = " << fake_MC_x_dipole_mapping[j_a][i_m] << endl;
	    if (fake_MC_x_dipole_mapping[j_a][i_m] == MC_x_dipole_mapping[i_a][j_m]){flag = j_m; break;}
	  }
	  if (flag == MC_x_dipole_mapping[i_a].size()){MC_x_dipole_mapping[i_a].push_back(fake_MC_x_dipole_mapping[j_a][i_m]);}
	}
      }
    }

    for (int i_a = 1; i_a < (*RA_dipole).size(); i_a++){
      for (int j_a = 1; j_a < fake_dipole[i_a].size(); j_a++){
	logger << LOG_DEBUG << "fake_dipole[" << i_a << "][" << j_a << "].name() = " << (*RA_dipole)[i_a].name() << " - " << fake_dipole[i_a][j_a].name() << endl;
      }
      /*
      for (int j_m = 0; j_m < MC_x_dipole_mapping[i_a].size(); j_m++){
	logger << LOG_DEBUG << "MC_x_dipole_mapping[" << i_a << "][" << j_m << "] = " << MC_x_dipole_mapping[i_a][j_m] << endl;
      }
      */
    }


    for (int i_a = 0; i_a < (*RA_dipole).size(); i_a++){
      for (int j_m = 0; j_m < MC_x_dipole_mapping[i_a].size(); j_m++){
	logger << LOG_DEBUG << "MC_x_dipole_mapping[" << i_a << "][" << j_m << "] = " << MC_x_dipole_mapping[i_a][j_m] << endl;
      }
    }



    logger << LOG_DEBUG << "fake_dipole to improve initial-state mappings finished..." << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



