#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

void initialization_fake_dipole_mapping_RT(phasespace_set & psi, call_generic & generic){
  Logger logger("initialization_fake_dipole_mapping_RT");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  /*
  if (!switch_MC_tau){
    logger << LOG_DEBUG << "No fake_dipole mappings used..." << endl;
    return;
  }
  */

  logger << LOG_DEBUG << "fake_dipole to improve initial-state mappings started..." << endl;

  vector<dipole_set> fake_dipole;
  vector<dipole_set> fake_dipole_candidate;
  phasespace_set fake_psi = psi;
  //    fake_psi.no_map.erase(fake_psi.no_map.begin() + 1, fake_psi.no_map.end());
  
  //  vector<int> fake_type_parton = psi_type_parton;
  //  vector<int> fake_basic_type_parton = psi_basic_type_parton;
  vector<int> fake_type_parton = psi.csi->type_parton[0];
  vector<int> fake_basic_type_parton = psi.csi->basic_type_parton[0];
  fake_dipole.push_back(dipole_set("ME2_R", fake_type_parton, fake_basic_type_parton, psi_symmetry_factor, psi_no_map[0], psi_o_map[0], psi_no_prc[0], psi_o_prc[0], 0));
  QCD_determine_dipoles(fake_dipole_candidate, fake_type_parton, fake_basic_type_parton);
  QCD_selection_fake_dipoles(fake_dipole, fake_dipole_candidate, psi_phasespace_order_alpha_s[0] - 1, psi_phasespace_order_alpha_e[0], psi_phasespace_order_interference[0], fake_psi.RA_singular_region, fake_psi.RA_singular_region_name, fake_psi.RA_singular_region_list, fake_psi, generic, generic.determination_no_subprocess_dipole, generic.determination_MCchannels_dipole);
  
  vector<vector<int> > fake_MC_x_dipole_mapping(fake_dipole.size());
  for (int j_a = 1; j_a < fake_dipole.size(); j_a++){
    generic.ac_tau_psp_dipole(j_a, fake_MC_x_dipole_mapping[j_a], fake_psi);
  }
  
  // determine psi_tau_MC_map from fake_MC_x_dipole_mapping (including 0!)

  for (int i_a = 0; i_a < fake_dipole.size(); i_a++){
    //  for (int i_a = 0; i_a < n_dipoles; i_a++){
    for (int i_m = 0; i_m < fake_MC_x_dipole_mapping[i_a].size(); i_m++){
      int flag = psi_tau_MC_map.size();
      for (int j_m = 0; j_m < psi_tau_MC_map.size(); j_m++){
	logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << "   i_m = " << i_m << "   j_m = " << j_m << endl;
	if (fake_MC_x_dipole_mapping[i_a][i_m] == psi_tau_MC_map[j_m]){flag = j_m; break;}
      }
      if (flag == psi_tau_MC_map.size()){psi_tau_MC_map.push_back(fake_MC_x_dipole_mapping[i_a][i_m]);}
    }
  }
  
  /*
      for (int j_a = 1; j_a < fake_dipole.size(); j_a++){
	for (int i_m = 0; i_m < fake_MC_x_dipole_mapping[j_a].size(); i_m++){
	  int flag = psi_MC_x_dipole_mapping.size();
	  for (int j_m = 0; j_m < psi_MC_x_dipole_mapping[i_a].size(); j_m++){
	    logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << "   i_m = " << i_m << "   j_m = " << j_m << endl;
	    if (fake_MC_x_dipole_mapping[j_a][i_m] == psi_MC_x_dipole_mapping[i_a][j_m]){flag = j_m; break;}
	  }
	  if (flag == psi_MC_x_dipole_mapping[i_a].size()){psi_MC_x_dipole_mapping[i_a].push_back(fake_MC_x_dipole_mapping[j_a][i_m]);}
	}
      }
    }
  */

  for (int j_a = 1; j_a < fake_dipole.size(); j_a++){
    logger << LOG_DEBUG << "fake_dipole[" << j_a << "].name() = " << fake_dipole[j_a].name() << endl;
  }
      /*
	for (int j_m = 0; j_m < psi_MC_x_dipole_mapping[i_a].size(); j_m++){
	logger << LOG_DEBUG << "psi_MC_x_dipole_mapping[" << i_a << "][" << j_m << "] = " << psi_MC_x_dipole_mapping[i_a][j_m] << endl;
	}
      */
  
  
  for (int j_m = 0; j_m < psi_tau_MC_map.size(); j_m++){
    logger << LOG_DEBUG << "psi_MC_tau[" << j_m << "] = " << psi_tau_MC_map[j_m] << endl;
  }
  
  logger << LOG_DEBUG << "fake_dipole to improve initial-state mappings finished..." << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void initialization_fake_dipole_mapping_RRA(vector<dipole_set> & dipole, phasespace_set & psi, call_generic & generic){
  Logger logger("initialization_fake_dipole_mapping_RRA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  int n_dipoles = dipole.size();
  logger << LOG_DEBUG << "fake_dipole to improve initial-state mappings started..." << endl;

  //    logger << LOG_DEBUG << "psi_n_dipoles = " << psi_n_dipoles << endl;
  logger << LOG_DEBUG << "REAL DIPOLE: dipole.size() = " << dipole.size() << endl;

    vector<vector<dipole_set> > fake_dipole(dipole.size());
    vector<vector<dipole_set> > fake_dipole_candidate(dipole.size());
    for (int i_a = 1; i_a < dipole.size(); i_a++){
      phasespace_set fake_psi;// = psi;
      fake_psi.no_map.resize(1, 0);
      fake_psi.no_prc.resize(1, 0);
      fake_psi.o_map.resize(1, vector<int> (psi.o_map.size() - 1));
      fake_psi.o_prc.resize(1, vector<int> (psi.o_prc.size() - 1));
      fake_psi.phasespace_order_alpha_s.resize(1, psi_phasespace_order_alpha_s[i_a]);
      fake_psi.phasespace_order_alpha_e.resize(1, psi_phasespace_order_alpha_e[i_a]);
      fake_psi.phasespace_order_interference.resize(1, psi_phasespace_order_interference[i_a]);
      fake_psi.MC_n_channel_phasespace.resize(1, psi_MC_n_channel_phasespace[i_a]);
      fake_psi.MC_sum_channel_phasespace.resize(1, psi_MC_sum_channel_phasespace[i_a]);

      fake_psi.M = psi.M;
      /*
      fake_psi.no_map.erase(fake_psi.no_map.begin() + 1, fake_psi.no_map.end());
      fake_psi.no_prc.erase(fake_psi.no_map.begin() + 1, fake_psi.no_prc.end());
      fake_psi.o_map.erase(fake_psi.o_map.begin() + 1, fake_psi.o_map.end());
      fake_psi.o_prc.erase(fake_psi.o_prc.begin() + 1, fake_psi.o_prc.end());
      fake_psi.phasespace_order_alpha_s.erase(fake_psi.phasespace_order_alpha_s.begin() + 1, fake_psi.phasespace_order_alpha_s.end());
      fake_psi.phasespace_order_alpha_e.erase(fake_psi.phasespace_order_alpha_e.begin() + 1, fake_psi.phasespace_order_alpha_e.end());
      fake_psi.phasespace_order_interference.erase(fake_psi.phasespace_order_interference.begin() + 1, fake_psi.phasespace_order_interference.end());
      fake_psi.phasespace_order_alpha_s[0] = psi_phasespace_order_alpha_s[i_a];
      fake_psi.phasespace_order_alpha_e[0] = psi_phasespace_order_alpha_e[i_a];
      fake_psi.phasespace_order_interference[0] = psi_phasespace_order_interference[i_a];
*/
      /*
      fake_psi.MC_n_channel_phasespace.erase(fake_psi.MC_n_channel_phasespace.begin() + 1, fake_psi.MC_n_channel_phasespace.end());
      fake_psi.MC_sum_channel_phasespace.erase(fake_psi.MC_sum_channel_phasespace.begin() + 1, fake_psi.MC_sum_channel_phasespace.end());
      */
      logger << LOG_DEBUG_VERBOSE << "psi_phasespace_order_alpha_s[" << i_a << "] - 1 = " << psi_phasespace_order_alpha_s[i_a] - 1 << endl;
      logger << LOG_DEBUG_VERBOSE << "psi_phasespace_order_alpha_s[" << i_a << "]     = " << psi_phasespace_order_alpha_e[i_a] << endl;
      logger << LOG_DEBUG_VERBOSE << "fake_psi.phasespace_order_alpha_s[" << 0 << "] = " << fake_psi.phasespace_order_alpha_s[0] << endl;
      logger << LOG_DEBUG_VERBOSE << "fake_psi.phasespace_order_alpha_s[" << 0 << "] = " << fake_psi.phasespace_order_alpha_e[0] << endl;
      logger << LOG_DEBUG_VERBOSE << "fake_psi.MC_n_channel_phasespace.size() = " << fake_psi.MC_n_channel_phasespace.size() << endl;
      logger << LOG_DEBUG_VERBOSE << "fake_psi.MC_sum_channel_phasespace.size() = " << fake_psi.MC_sum_channel_phasespace.size() << endl;


      vector<int> fake_type_parton = dipole[i_a].type_parton();
      vector<int> fake_basic_type_parton = dipole[i_a].basic_type_parton();
      fake_dipole[i_a].push_back(dipole_set(dipole[i_a].name(), fake_type_parton, fake_basic_type_parton, psi_symmetry_factor, psi_no_map[i_a], psi_o_map[i_a], psi_no_prc[i_a], psi_o_prc[i_a], 0));
      QCD_determine_dipoles(fake_dipole_candidate[i_a], fake_type_parton, fake_basic_type_parton);
      logger << LOG_DEBUG_VERBOSE << "REAL DIPOLE " << i_a << endl;



      QCD_selection_fake_dipoles(fake_dipole[i_a], fake_dipole_candidate[i_a], fake_psi.phasespace_order_alpha_s[0] - 1, fake_psi.phasespace_order_alpha_e[0], psi_phasespace_order_interference[i_a], fake_psi.RA_singular_region, fake_psi.RA_singular_region_name, fake_psi.RA_singular_region_list, fake_psi, generic, generic.determination_no_subprocess_doubledipole, generic.determination_MCchannels_doubledipole);



      //      QCD_selection_fake_dipoles(fake_dipole[i_a], fake_dipole_candidate[i_a], psi_phasespace_order_alpha_s[i_a] - 1, psi_phasespace_order_alpha_e[i_a], psi_phasespace_order_interference[i_a], fake_psi.RA_singular_region, fake_psi.RA_singular_region_name, fake_psi.RA_singular_region_list, fake_psi, generic, generic.determination_no_subprocess_doubledipole, generic.determination_MCchannels_doubledipole);

      vector<vector<int> > fake_MC_x_dipole_mapping(fake_dipole[i_a].size());
      logger << LOG_DEBUG_VERBOSE << "fake_dipole[" << i_a << "].size() = " << fake_dipole[i_a].size() << endl;
      for (int j_a = 1; j_a < fake_dipole[i_a].size(); j_a++){
	logger << LOG_DEBUG_VERBOSE << "before: fake_MC_x_dipole_mapping[" << j_a << "].size() = " << fake_MC_x_dipole_mapping[j_a].size() << endl;
	generic.ac_tau_psp_doubledipole(j_a, fake_MC_x_dipole_mapping[j_a], fake_psi);
	logger << LOG_DEBUG_VERBOSE << "after:  fake_MC_x_dipole_mapping[" << j_a << "].size() = " << fake_MC_x_dipole_mapping[j_a].size() << endl;
      }

      logger << LOG_DEBUG_VERBOSE << "fake_dipole[" << i_a << "].size() = " << fake_dipole[i_a].size() << endl;
      for (int j_a = 1; j_a < fake_dipole[i_a].size(); j_a++){
	logger << LOG_DEBUG_VERBOSE << "fake_MC_x_dipole_mapping[" << j_a << "].size() = " << fake_MC_x_dipole_mapping[j_a].size() << endl;
	for (int i_m = 0; i_m < fake_MC_x_dipole_mapping[j_a].size(); i_m++){
	  int flag = psi_MC_x_dipole_mapping[i_a].size();
	  for (int j_m = 0; j_m < psi_MC_x_dipole_mapping[i_a].size(); j_m++){
	    logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << "   i_m = " << i_m << "   j_m = " << j_m << "   fake_MC_x_dipole_mapping[" << j_a << "][" << i_m << "] = " << fake_MC_x_dipole_mapping[j_a][i_m] << endl;
	    if (fake_MC_x_dipole_mapping[j_a][i_m] == psi_MC_x_dipole_mapping[i_a][j_m]){flag = j_m; break;}
	  }
	  if (flag == psi_MC_x_dipole_mapping[i_a].size()){psi_MC_x_dipole_mapping[i_a].push_back(fake_MC_x_dipole_mapping[j_a][i_m]);}
	}
      }
    }

    for (int i_a = 1; i_a < dipole.size(); i_a++){
      for (int j_a = 1; j_a < fake_dipole[i_a].size(); j_a++){
	logger << LOG_DEBUG << "fake_dipole[" << i_a << "][" << j_a << "].name() = " << dipole[i_a].name() << " - " << fake_dipole[i_a][j_a].name() << endl;
      }
      /*
      for (int j_m = 0; j_m < psi_MC_x_dipole_mapping[i_a].size(); j_m++){
	logger << LOG_DEBUG << "psi_MC_x_dipole_mapping[" << i_a << "][" << j_m << "] = " << psi_MC_x_dipole_mapping[i_a][j_m] << endl;
      }
      */
    }


    for (int i_a = 0; i_a < dipole.size(); i_a++){
      for (int j_m = 0; j_m < psi_MC_x_dipole_mapping[i_a].size(); j_m++){
	logger << LOG_DEBUG << "psi_MC_x_dipole_mapping[" << i_a << "][" << j_m << "] = " << psi_MC_x_dipole_mapping[i_a][j_m] << endl;
      }
    }



    logger << LOG_DEBUG << "fake_dipole to improve initial-state mappings finished..." << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



