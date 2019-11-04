#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

void observable_set::determine_ioperator_QCD(phasespace_set & psi){
  Logger logger("observable_set::CS_determine_ioperator_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static map<int, double> charge_particle;
  if (initialization){
    fill_charge_particle(charge_particle);
    initialization = 0;
  }
  vector<string> pa_name(csi->type_parton[0].size(), "");
  if (csi->type_parton[0][1] > -10 && csi->type_parton[0][1] < 10){pa_name[1] = "a";}
  if (csi->type_parton[0][2] > -10 && csi->type_parton[0][2] < 10){pa_name[2] = "b";}
  int count = 0;
  vector<string> alphabet(csi->type_parton[0].size() - 3, "");
  for (int i_p = 0; i_p < alphabet.size(); i_p++){alphabet[i_p] = char(105 + i_p);}
  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (csi->type_parton[0][i_p] > -10 && csi->type_parton[0][i_p] < 10){pa_name[i_p] = alphabet[count++];}}

  int temp_type_correction = 1;

  for (int temp_no_emitter = 1; temp_no_emitter < csi->type_parton[0].size(); temp_no_emitter++){
    for (int temp_no_spectator = 1; temp_no_spectator < csi->type_parton[0].size(); temp_no_spectator++){
      if (temp_no_emitter == temp_no_spectator){continue;}
      if (pa_name[temp_no_emitter] == "" || pa_name[temp_no_spectator] == ""){continue;}
      vector<int> temp_pair(2);
      temp_pair[0] = std::min(temp_no_emitter, temp_no_spectator);
      temp_pair[1] = std::max(temp_no_emitter, temp_no_spectator);

      double temp_charge_factor = 1.;
      
      int temp_type;
      if (csi->type_parton[0][temp_no_emitter] == 0){temp_type = 0;}
      else {temp_type = 1;}

      int temp_massive;
      if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 0;}
      else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 1;}
      else if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 2;}
      else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 3;}
      else {logger << LOG_WARN << "Should not happen!" << endl;}

      string temp_name;
      temp_name = "I^{" + pa_name[temp_no_emitter] + "," + pa_name[temp_no_spectator] + "}";

      int flag = (*VA_ioperator).size();
      for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
	if (temp_no_emitter == (*VA_ioperator)[i_a][0].no_emitter() && temp_type_correction == (*VA_ioperator)[i_a][0].type_correction()){
	  flag = i_a; 
	  break;
	}
      }
      if (flag == (*VA_ioperator).size()){(*VA_ioperator).push_back(vector<ioperator_set> ());}
      //      flag = (*VA_ioperator).size() - 1;
      //      }
      (*VA_ioperator)[flag].push_back(ioperator_set(temp_name, temp_type, temp_pair, psi_no_prc[0], csi->type_parton[0], temp_charge_factor, temp_no_emitter, temp_no_spectator, temp_type_correction, temp_massive));
    }
  }
  for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
    for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
      logger << LOG_DEBUG << "(*VA_ioperator)[" << i_a << "][" << j_a << "].name() = " << (*VA_ioperator)[i_a][j_a].name() << endl;
    }
  }


  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::calculate_ioperator_QCD(){
  Logger logger("observable_set::calculate_ioperator_QCD_CDST");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (massive_QCD){calculate_ioperator_QCD_CDST();}
  else {calculate_ioperator_QCD_CS();}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_ioperator_QCD_CS(){
  Logger logger("observable_set::calculate_ioperator_QS_CS");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static vector<vector<int> > pair;
  static vector<vector<double> > ioperator_pair_log(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > ioperator_pair_log2_2(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<int> type;
  static double prefactor_log[2][2];
  //  static double alpha_S_2pi = alpha_S * inv2pi;
  static vector<vector<int> > type_emitter((*VA_ioperator).size());
  static vector<vector<vector<int> > > ppair((*VA_ioperator).size());

  if (initialization){
    for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
      if ((*VA_ioperator)[i_a][0].type_correction() == 1){
	for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
	  int flag = -1;
	  for (int i_p = 0; i_p < pair.size(); i_p++){if ((*VA_ioperator)[i_a][j_a].pair() == pair[i_p]){flag = i_p; break;}}
	  if (flag == -1){pair.push_back((*VA_ioperator)[i_a][j_a].pair());}
	  flag = -1;
	  for (int i_t = 0; i_t < type.size(); i_t++){if ((*VA_ioperator)[i_a][j_a].type() == type[i_t]){flag = i_t; break;}}
	  if (flag == -1){type.push_back((*VA_ioperator)[i_a][j_a].type());}
	}
      }
    }
    for (int i_t = 0; i_t < type.size(); i_t++){
      if (switch_VI_bosonic_fermionic == 0){
	// both fermionic and bosonic contributions
	if (type[i_t] == 0){
	  prefactor_log[0][0] = VA_DeltaIR2 - pi2_3 + ((VA_DeltaIR1 + 1.) * gamma_g(N_f) + K_g(N_f)) / C_A;
	  // Gamma(1 + eps) = 1 - gamma eps + 1/12 * (6 gamma² + pi²) eps² + O(eps³)
	  // 1 / Gamma(1 - eps) = 1 - gamma eps + 1/12 * (6 gamma² - pi²) eps² + O(eps³)
	  // 1 / Gamma(1 - eps) = Gamma(1 + eps) - pi²/6 eps² + O(eps³)
	  // switch_polenrom = 0: (4pi)^eps * Gamma(1 + eps) * (c_0 + c_1 * 1/eps + c_2 * 1/eps²)
	  // switch_polenrom = 1: (4pi)^eps / Gamma(1 - eps) * (c_0 + c_1 * 1/eps + c_2 * 1/eps²)
	  if (switch_polenorm == 1){prefactor_log[0][0] += -pi2_6;} // BLHA -> COLI !!!
	  prefactor_log[0][1] = VA_DeltaIR1 + gamma_g(N_f) / C_A;
	}
	else if (type[i_t] == 1){
	  prefactor_log[1][0] = VA_DeltaIR2 - pi2_3 + ((VA_DeltaIR1 + 1.) * gamma_q + K_q) / C_F;
	  if (switch_polenorm == 1){prefactor_log[1][0] += -pi2_6;} // BLHA -> COLI !!!
	  prefactor_log[1][1] = VA_DeltaIR1 + gamma_q / C_F;
	}
      }
      else if (switch_VI_bosonic_fermionic == 1){
	// only bosonic contributions
	if (type[i_t] == 0){
	  prefactor_log[0][0] = VA_DeltaIR2 - pi2_3 + ((VA_DeltaIR1 + 1.) * gamma_g_bos + K_g_bos) / C_A;
	  if (switch_polenorm == 1){prefactor_log[0][0] += -pi2_6;} // BLHA -> COLI !!!
	  prefactor_log[0][1] = VA_DeltaIR1 + gamma_g_bos / C_A;
	}
	else if (type[i_t] == 1){
	  prefactor_log[1][0] = VA_DeltaIR2 - pi2_3 + ((VA_DeltaIR1 + 1.) * gamma_q + K_q) / C_F;
	  if (switch_polenorm == 1){prefactor_log[1][0] += -pi2_6;} // BLHA -> COLI !!!
	  prefactor_log[1][1] = VA_DeltaIR1 + gamma_q / C_F;
	}
      }
      else if (switch_VI_bosonic_fermionic == 2){
	// only fermionic contribution
	if (type[i_t] == 0){
	  prefactor_log[0][0] = ((VA_DeltaIR1 + 1.) * gamma_g_ferm(N_f) + K_g_ferm(N_f)) / C_A;
	  prefactor_log[0][1] = gamma_g_ferm(N_f) / C_A;
	}
	else if (type[i_t] == 1){
	  prefactor_log[1][0] = 0.;
	  prefactor_log[1][1] = 0.;
	}
      }
    }

    for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
      if ((*VA_ioperator)[i_a][0].type_correction() == 1){
	type_emitter[i_a].resize((*VA_ioperator)[i_a].size());
	ppair[i_a].resize((*VA_ioperator)[i_a].size());
	for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
	  type_emitter[i_a][j_a] = (*VA_ioperator)[i_a][j_a].type();
	  ppair[i_a][j_a] = (*VA_ioperator)[i_a][j_a].pair();
	}
      } 
    }

    for (int i_p = 0; i_p < pair.size(); i_p++){logger << LOG_DEBUG_VERBOSE << "pair[" << i_p << "] = " << "(" << pair[i_p][0] << ", " << pair[i_p][1] << ")" << endl;}
    for (int i_t = 0; i_t < type.size(); i_t++){logger << LOG_DEBUG_VERBOSE << "type[" << i_t << "] = " << type[i_t] << endl;}
    //    if (initialization){initialization = 0;}
    if (!VA_delta_flag){initialization = 0;}
  }

  double mu2IR = pow(var_mu_ren, 2);
  // modified to deal with mu_ren != mu_reg !!!
  if (user.double_value[user.double_map["mu_reg"]] != 0. && user.double_value[user.double_map["mu_reg"]] != -1.){mu2IR = pow(user.double_value[user.double_map["mu_reg"]], 2);}

  //  logger << LOG_DEBUG << "var_mu_ren = " << var_mu_ren << endl;
  //  logger << LOG_DEBUG << "mu2IR = " << mu2IR << endl;
  //  logger << LOG_DEBUG << "alpha_S_2pi = " << alpha_S_2pi << endl;

    for (int i_p = 0; i_p < pair.size(); i_p++){
    ioperator_pair_log[pair[i_p][0]][pair[i_p][1]] = log(mu2IR / (2. * p_parton[0][pair[i_p][0]] * p_parton[0][pair[i_p][1]]));
    ioperator_pair_log2_2[pair[i_p][0]][pair[i_p][1]] = pow(ioperator_pair_log[pair[i_p][0]][pair[i_p][1]], 2) / 2.;
  }

  for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
    for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
      if ((*VA_ioperator)[i_a][0].type_correction() == 1){
	VA_I_ME2_cf[i_a][j_a] = -alpha_S * inv2pi * (prefactor_log[type_emitter[i_a][j_a]][0] + prefactor_log[type_emitter[i_a][j_a]][1] * ioperator_pair_log[ppair[i_a][j_a][0]][ppair[i_a][j_a][1]] + ioperator_pair_log2_2[ppair[i_a][j_a][0]][ppair[i_a][j_a][1]]) * VA_ME2_cf[i_a][j_a];
	//	logger << LOG_DEBUG << "VA_I_ME2_cf[" << i_a << "][" << j_a << "] = " << VA_I_ME2_cf[i_a][j_a] << endl;
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_ioperator_QCD_CDST(){
  Logger logger("observable_set::calculate_ioperator_QCD_CDST");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static vector<vector<int> > pair;
  static vector<int> pair_massive;
  static vector<int> pair_code_emitter;
  //  static double prefactor_log[2][2];
  //  static double alpha_S_2pi = alpha_S * inv2pi;
  //  static int N_f = osi_N_f;

  static vector<vector<int> > type_emitter((*VA_ioperator).size());
  static vector<vector<int> > no_emitter((*VA_ioperator).size());
  static vector<vector<int> > no_spectator((*VA_ioperator).size());
  static vector<vector<vector<int> > > ppair((*VA_ioperator).size());

  static vector<double> pmass(csi->type_parton[0].size());
  static vector<double> pmass2(csi->type_parton[0].size());

  static vector<vector<double> > s_ik(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > Q2_ik(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > sqrtQ2_ik(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > mu2_i(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > rho2_i(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > v_ik(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > delta(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));

  static double temp_finite_massive_full_g_in;
  static double temp_finite_massive_full_g_out;
  static double temp_finite_massive_full_q;
  static vector<vector<double> > V_NS_massive_full(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<vector<double> > > V_S_massive_full(csi->type_parton[0].size(), vector<vector<double> > (csi->type_parton[0].size(), vector<double> (3)));
  static vector<vector<vector<double> > > finite_massive_full(csi->type_parton[0].size(), vector<vector<double> > (csi->type_parton[0].size(), vector<double> (2)));
  static vector<vector<double> > ioperator_pair_log(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > ioperator_pair_log2_2(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));

  static double Gamma_g;
  static double Gamma_g_out;
  static vector<vector<double> > Gamma_q(csi->type_parton[0].size(), vector<double> (2));

  if (initialization){
    for (int i_p = 0; i_p < csi->type_parton[0].size(); i_p++){
      pmass[i_p] = mass_parton[0][i_p];
      pmass2[i_p] = mass2_parton[0][i_p];
      if (pmass[i_p] == 0.){
	Gamma_q[i_p][0] = gamma_q * VA_DeltaIR1;
	Gamma_q[i_p][1] = 0.;
      }
      else if (pmass[i_p] != 0.){
	//    Gamma_q[i_p][0] = C_F * (VA_DeltaIR1 + .5 * log(m2_i / s_ik) - 2.);
	Gamma_q[i_p][1] = -.5 * C_F; // gamma_q = 1.5 * C_F !!! check if nothing goes wrong here !!! (explicit mu dependence)
      }
    }
    // Gamma_g not treated here -> depends on Q2aux
    
    for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
      if ((*VA_ioperator)[i_a][0].type_correction() == 1){
	for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
	  int flag = -1;
	  for (int i_p = 0; i_p < pair.size(); i_p++){if ((*VA_ioperator)[i_a][j_a].pair() == pair[i_p]){flag = i_p; break;}}
	  if (flag == -1){
	    pair.push_back((*VA_ioperator)[i_a][j_a].pair());
	    pair_massive.push_back((*VA_ioperator)[i_a][j_a].massive());
	    int temp_code_emitter;
	    if ((*VA_ioperator)[i_a][j_a].pair()[0] == (*VA_ioperator)[i_a][j_a].no_emitter()){temp_code_emitter = 1;}
	    else if ((*VA_ioperator)[i_a][j_a].pair()[1] == (*VA_ioperator)[i_a][j_a].no_emitter()){temp_code_emitter = 2;}
	    pair_code_emitter.push_back(temp_code_emitter);
	  }
	  else {pair_code_emitter[flag] = 3;}
	  /*
	    flag = 0;
	    for (int i_t = 0; i_t < type.size(); i_t++){if ((*VA_ioperator)[i_a][j_a].type() == type[i_t]){flag = 1; break;}}
	    if (flag == 0){type.push_back((*VA_ioperator)[i_a][j_a].type());}
	  */
	}
      }
    }

    for (int i_p = 0; i_p < pair.size(); i_p++){
      int ip_0 = pair[i_p][0];
      int ip_1 = pair[i_p][1];
      if (pair_massive[i_p] == 3){
	if (pair_code_emitter[i_p] != 2){V_S_massive_full[ip_0][ip_1][2] = 0.;}
	if (pair_code_emitter[i_p] != 1){V_S_massive_full[ip_1][ip_0][2] = 0.;}
      }
      else if (pair_massive[i_p] == 1 || pair_massive[i_p] == 2){
	if (pair_code_emitter[i_p] != 2){V_S_massive_full[ip_0][ip_1][2] = .5;}
	if (pair_code_emitter[i_p] != 1){V_S_massive_full[ip_1][ip_0][2] = .5;}
      }
      else if (pair_massive[i_p] == 0){
	// m_i == 0 (gluon/quark) and m_k == 0 (gluon/quark)
	if (pair_code_emitter[i_p] != 2){
	  V_S_massive_full[ip_0][ip_1][0] = VA_DeltaIR2;
	  if (switch_polenorm == 1){V_S_massive_full[ip_0][ip_1][0] += -pi2_6;} // BLHA -> COLI !!!
	  V_S_massive_full[ip_0][ip_1][1] = VA_DeltaIR1;
	  V_S_massive_full[ip_0][ip_1][2] = 1.;
	  if (csi->type_parton[0][ip_0] == 0 && ip_0 < 3){V_NS_massive_full[ip_0][ip_1] = 0.;}
	  else if (csi->type_parton[0][ip_0] != 0){V_NS_massive_full[ip_0][ip_1] = 0.;}
	}
	if (pair_code_emitter[i_p] != 1){
	  V_S_massive_full[ip_1][ip_0][0] = VA_DeltaIR2;
	  if (switch_polenorm == 1){V_S_massive_full[ip_1][ip_0][0] += -pi2_6;} // BLHA -> COLI !!!
	  V_S_massive_full[ip_1][ip_0][1] = VA_DeltaIR1;
	  V_S_massive_full[ip_1][ip_0][2] = 1.;
	  if (csi->type_parton[0][ip_1] == 0 && ip_1 < 3){V_NS_massive_full[ip_1][ip_0] = 0.;}
	  if (csi->type_parton[0][ip_1] != 0){V_NS_massive_full[ip_1][ip_0] = 0.;}
	}
      }
    
      if (pair_code_emitter[i_p] != 2){
	if (csi->type_parton[0][ip_0] == 0){
	  //	  finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_g + gamma_g(N_f) + K_g(N_f)) / C_A;   // Gamma_g depends on Q2_aux -> later
	  finite_massive_full[ip_0][ip_1][1] = + gamma_g(N_f) / C_A;
	}
	else if (csi->type_parton[0][ip_0] != 0){
	  if (pmass[ip_0] == 0){finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_q[ip_0][0] + gamma_q + K_q) / C_F;}
	  //	  else {finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_q[ip_0][0] + gamma_q + K_q) / C_F;}   // Gamma_q[0] depends on phase-space
	  finite_massive_full[ip_0][ip_1][1] = + (Gamma_q[ip_0][1] + gamma_q) / C_F;
	}
      }
      if (pair_code_emitter[i_p] != 1){
	if (csi->type_parton[0][ip_1] == 0){
	  //	  finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_g + gamma_g(N_f) + K_g(N_f)) / C_A;   // Gamma_g depends on Q2_aux -> later
	  finite_massive_full[ip_1][ip_0][1] = + gamma_g(N_f) / C_A;
	}
	else if (csi->type_parton[0][ip_1] != 0){
	  if (pmass[ip_1] == 0){finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_q[ip_1][0] + gamma_q + K_q) / C_F;}
	  //	  else {finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_q[ip_1][0] + gamma_q + K_q) / C_F;}   // Gamma_q[0] depends on phase-space
	  finite_massive_full[ip_1][ip_0][1] = + (Gamma_q[ip_1][1] + gamma_q) / C_F;
	}
      }
    }
 
    for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
      if ((*VA_ioperator)[i_a][0].type_correction() == 1){
	type_emitter[i_a].resize((*VA_ioperator)[i_a].size());
	no_emitter[i_a].resize((*VA_ioperator)[i_a].size());
	no_spectator[i_a].resize((*VA_ioperator)[i_a].size());
	ppair[i_a].resize((*VA_ioperator)[i_a].size());
	for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
	  type_emitter[i_a][j_a] = (*VA_ioperator)[i_a][j_a].type();
	  no_emitter[i_a][j_a] = (*VA_ioperator)[i_a][j_a].no_emitter();
	  no_spectator[i_a][j_a] = (*VA_ioperator)[i_a][j_a].no_spectator();
	  ppair[i_a][j_a] = (*VA_ioperator)[i_a][j_a].pair();
	}
      }
    }

    //    logger << LOG_DEBUG << "pair_massive[" << i_p << "] = " << pair_massive[i_p] << endl;
    for (int i_p = 0; i_p < pair.size(); i_p++){logger << LOG_DEBUG << "pair[" << i_p << "] = " << "(" << pair[i_p][0] << ", " << pair[i_p][1] << ")   massive = " << pair_massive[i_p] << "   code_emitter = " << pair_code_emitter[i_p] << endl;}
    //    for (int i_t = 0; i_t < type.size(); i_t++){logger << LOG_DEBUG << "type[" << i_t << "] = " << type[i_t] << endl;}
    //    initialization = 0;
    //    if (initialization){initialization = 0;}
    if (!VA_delta_flag){initialization = 0;}
  }

  //  double kappa = 2. / 3.;
  double mu2IR = pow(var_mu_ren, 2);
  // modified to deal with mu_ren != mu_reg !!!
  if (user.double_value[user.double_map["mu_reg"]] != 0. && user.double_value[user.double_map["mu_reg"]] != -1.){mu2IR = pow(user.double_value[user.double_map["mu_reg"]], 2);}

  double Q2_aux = mu2IR;

  // needed only if a gluon is around !!!
  Gamma_g = gamma_g(N_f) * VA_DeltaIR1;
  Gamma_g_out = Gamma_g;
  for (int i_p = 1; i_p < 7; i_p++){if (M[i_p] > 0.){Gamma_g_out += -2./3. * T_R * log(M2[i_p] / Q2_aux);}}
  temp_finite_massive_full_g_in = - pi2_3 + (Gamma_g + gamma_g(N_f) + K_g(N_f)) / C_A;
  temp_finite_massive_full_g_out = - pi2_3 + (Gamma_g_out + gamma_g(N_f) + K_g(N_f)) / C_A;
  // needed only if a (anti)quark is around !!!
  temp_finite_massive_full_q = - pi2_3 + (gamma_q + K_q) / C_F;


  for (int i_p = 0; i_p < pair.size(); i_p++){
    int ip_0 = pair[i_p][0];
    int ip_1 = pair[i_p][1];
    // always needed
    s_ik[ip_0][ip_1] = 2. * p_parton[0][ip_0] * p_parton[0][ip_1];
    ioperator_pair_log[ip_0][ip_1] = log(mu2IR / s_ik[ip_0][ip_1]);
    ioperator_pair_log2_2[ip_0][ip_1] = pow(ioperator_pair_log[ip_0][ip_1], 2) / 2.;
    
    if (pair_code_emitter[i_p] != 2){
      if (csi->type_parton[0][ip_0] == 0){
	if (ip_0 < 3){finite_massive_full[ip_0][ip_1][0] = temp_finite_massive_full_g_in;}
	else if (ip_0 > 2){finite_massive_full[ip_0][ip_1][0] = temp_finite_massive_full_g_out;}
      }
      else if (csi->type_parton[0][ip_0] != 0){
	if (pmass[ip_0] != 0){	
	  Gamma_q[ip_0][0] = C_F * (VA_DeltaIR1 + .5 * log(pmass2[ip_0] / s_ik[ip_0][ip_1]) - 2.);
	  finite_massive_full[ip_0][ip_1][0] = temp_finite_massive_full_q + Gamma_q[ip_0][0] / C_F;
	}
      }
    }   
    if (pair_code_emitter[i_p] != 1){
      if (csi->type_parton[0][ip_1] == 0){
	if (ip_1 < 3){finite_massive_full[ip_1][ip_0][0] = temp_finite_massive_full_g_in;}
	else if (ip_1 > 2){finite_massive_full[ip_1][ip_0][0] = temp_finite_massive_full_g_out;}
      }
      else if (csi->type_parton[0][ip_1] != 0){
	if (pmass[ip_1] != 0){	
	  Gamma_q[ip_1][0] = C_F * (VA_DeltaIR1 + .5 * log(pmass2[ip_1] / s_ik[ip_0][ip_1]) - 2.);
	  finite_massive_full[ip_1][ip_0][0] = temp_finite_massive_full_q + Gamma_q[ip_1][0] / C_F;
	}
      }
    }
      
      
    if (pair_massive[i_p] == 3){
    // m_i > 0 (quark) && m_k > 0 (quark) 
      Q2_ik[ip_0][ip_1] = pmass2[ip_0] + pmass2[ip_1] + s_ik[ip_0][ip_1];
      sqrtQ2_ik[ip_0][ip_1] = sqrt(Q2_ik[ip_0][ip_1]);
      mu2_i[ip_0][ip_1] = pmass2[ip_0] / Q2_ik[ip_0][ip_1];
      mu2_i[ip_1][ip_0] = pmass2[ip_1] / Q2_ik[ip_0][ip_1];
      double temp_r2jk = pmass2[ip_0] * pmass2[ip_1] / pow(s_ik[ip_0][ip_1], 2);
      if (temp_r2jk < 1.e-3){
	delta[ip_0][ip_1] = 2 * temp_r2jk + 2 * pow(temp_r2jk, 2) + 4 * pow(temp_r2jk, 3) + 10 * pow(temp_r2jk, 4) + 28 * pow(temp_r2jk, 5);
	v_ik[ip_0][ip_1] = 1. - delta[ip_0][ip_1];
      }
      else {
	v_ik[ip_0][ip_1] = sqrt(lambda(1., mu2_i[ip_0][ip_1], mu2_i[ip_1][ip_0])) / (1. - mu2_i[ip_0][ip_1] - mu2_i[ip_1][ip_0]);
	delta[ip_0][ip_1] = 1. - v_ik[ip_0][ip_1];
      }
      double rho2 = delta[ip_0][ip_1] / (1. + v_ik[ip_0][ip_1]);
      double logrho2 = log(rho2);
      double logrho = log(sqrt(rho2));
      rho2_i[ip_0][ip_1] = (1. - v_ik[ip_0][ip_1] + 2 * mu2_i[ip_0][ip_1] / (1. - mu2_i[ip_0][ip_1] - mu2_i[ip_1][ip_0])) / (1. + v_ik[ip_0][ip_1] + 2 * mu2_i[ip_0][ip_1] / (1. - mu2_i[ip_0][ip_1] - mu2_i[ip_1][ip_0]));
      rho2_i[ip_1][ip_0] = (1. - v_ik[ip_0][ip_1] + 2 * mu2_i[ip_1][ip_0] / (1. - mu2_i[ip_0][ip_1] - mu2_i[ip_1][ip_0])) / (1. + v_ik[ip_0][ip_1] + 2 * mu2_i[ip_1][ip_0] / (1. - mu2_i[ip_0][ip_1] - mu2_i[ip_1][ip_0]));
      double temp_V_S_massive_pair_0 = 1. / v_ik[ip_0][ip_1] * (logrho * VA_DeltaIR1 -.25 * pow(log(rho2_i[ip_0][ip_1]), 2) -.25 * pow(log(rho2_i[ip_1][ip_0]), 2) - pi2_6 + logrho * log(Q2_ik[ip_0][ip_1] / s_ik[ip_0][ip_1]));
      double temp_V_S_massive_pair_1 = logrho / v_ik[ip_0][ip_1];
      double temp_V_NS_massive_pair = gamma_q / C_F * log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) + 1. / v_ik[ip_0][ip_1] * (logrho2 * log(1. + rho2) + 2 * gsl_sf_dilog(rho2) - gsl_sf_dilog(1. - rho2_i[ip_0][ip_1]) - gsl_sf_dilog(1. - rho2_i[ip_1][ip_0]) - pi2_6);
      if (pair_code_emitter[i_p] != 2){
	V_S_massive_full[ip_0][ip_1][0] = temp_V_S_massive_pair_0;
	V_S_massive_full[ip_0][ip_1][1] = temp_V_S_massive_pair_1;
	//      V_S_massive_full[ip_0][ip_1][2] = 0.;
	V_NS_massive_full[ip_0][ip_1] = temp_V_NS_massive_pair + log((sqrtQ2_ik[ip_0][ip_1] - pmass[ip_1]) / sqrtQ2_ik[ip_0][ip_1]) - 2 * log((pow(sqrtQ2_ik[ip_0][ip_1] - pmass[ip_1], 2) - pmass2[ip_0]) / Q2_ik[ip_0][ip_1]) - 2 * pmass2[ip_0] / s_ik[ip_0][ip_1] * log(pmass[ip_0] / (sqrtQ2_ik[ip_0][ip_1] - pmass[ip_1])) - pmass[ip_1] / (sqrtQ2_ik[ip_0][ip_1] - pmass[ip_1]) + 2 * pmass[ip_1] * (2 * pmass[ip_1] - sqrtQ2_ik[ip_0][ip_1]) / s_ik[ip_0][ip_1] + pi2_2;
	/*
	if (csi->type_parton[0][ip_0] == 0){
	  finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_g + gamma_g(N_f) + K_g(N_f)) / C_A;
	}
	else if (csi->type_parton[0][ip_0] != 0){
	  Gamma_q[ip_0][0] = C_F * (VA_DeltaIR1 + .5 * log(pmass2[ip_0] / s_ik[ip_0][ip_1]) - 2.);
	  finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_q[ip_0][0] + gamma_q + K_q) / C_F;
	}
	*/
      }
      
      if (pair_code_emitter[i_p] != 1){
	V_S_massive_full[ip_1][ip_0][0] = temp_V_S_massive_pair_0;
	V_S_massive_full[ip_1][ip_0][1] = temp_V_S_massive_pair_1;
	//      V_S_massive_full[ip_0][ip_1][2] = 0.;
	V_NS_massive_full[ip_1][ip_0] = temp_V_NS_massive_pair + log((sqrtQ2_ik[ip_0][ip_1] - pmass[ip_0]) / sqrtQ2_ik[ip_0][ip_1]) - 2 * log((pow(sqrtQ2_ik[ip_0][ip_1] - pmass[ip_0], 2) - pmass2[ip_1]) / Q2_ik[ip_0][ip_1]) - 2 * pmass2[ip_1] / s_ik[ip_0][ip_1] * log(pmass[ip_1] / (sqrtQ2_ik[ip_0][ip_1] - pmass[ip_0])) - pmass[ip_0] / (sqrtQ2_ik[ip_0][ip_1] - pmass[ip_0]) + 2 * pmass[ip_0] * (2 * pmass[ip_0] - sqrtQ2_ik[ip_0][ip_1]) / s_ik[ip_0][ip_1] + pi2_2;
	/*
	if (csi->type_parton[0][ip_1] == 0){
	  finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_g + gamma_g(N_f) + K_g(N_f)) / C_A;
	}
	else if (csi->type_parton[0][ip_1] != 0){
	  Gamma_q[ip_1][0] = C_F * (VA_DeltaIR1 + .5 * log(pmass2[ip_1] / s_ik[ip_0][ip_1]) - 2.);
	  finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_q[ip_1][0] + gamma_q + K_q) / C_F;
	}
	*/
      }
    }
    else if (pair_massive[i_p] == 1 || pair_massive[i_p] == 2){
      double temp_mass2;
      if (pair_massive[i_p] == 1){temp_mass2 = pmass2[ip_0];}
      else if (pair_massive[i_p] == 2){temp_mass2 = pmass2[ip_1];}
      else {logger << LOG_WARN << "pair_massive[i_p = " << i_p << "] = " << pair_massive[i_p] << " --- should not happen!" << endl;}
      Q2_ik[ip_0][ip_1] = temp_mass2 + s_ik[ip_0][ip_1];
      double temp_log_s_Q2_ik = log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]);

      double temp_V_S_massive_pair_0 = .5 * VA_DeltaIR2 + .5 * log(temp_mass2 / s_ik[ip_0][ip_1]) * VA_DeltaIR1 - 0.25 * pow(log(temp_mass2 / s_ik[ip_0][ip_1]), 2) - pi2_12 - .5 * log(temp_mass2 / s_ik[ip_0][ip_1]) * temp_log_s_Q2_ik - .5 * log(temp_mass2 / Q2_ik[ip_0][ip_1]) * temp_log_s_Q2_ik;
      double temp_V_S_massive_pair_1 = .5 * VA_DeltaIR1 + .5 * log(temp_mass2 / s_ik[ip_0][ip_1]);
      if (switch_polenorm == 1){temp_V_S_massive_pair_0 += -.5 * pi2_6;} // BLHA -> COLI !!!
      if (pair_code_emitter[i_p] != 2){
	V_S_massive_full[ip_0][ip_1][0] = temp_V_S_massive_pair_0;
	V_S_massive_full[ip_0][ip_1][1] = temp_V_S_massive_pair_1;
      }
      if (pair_code_emitter[i_p] != 1){
	V_S_massive_full[ip_1][ip_0][0] = temp_V_S_massive_pair_0;
	V_S_massive_full[ip_1][ip_0][1] = temp_V_S_massive_pair_1;
      }
      if ((pair_massive[i_p] == 1 && pair_code_emitter[i_p] != 2) || (pair_massive[i_p] == 2 && pair_code_emitter[i_p] != 1)){
	// m_i > 0 (quark) and m_k == 0 (gluon/quark)
	//	  logger << LOG_DEBUG_VERBOSE << "m_i > 0 (quark) and m_k == 0 (gluon/quark):   i = " << dx_pa[ioperator_pair[ia][ib]][0] << "   k = " << dx_pa[ioperator_pair[ia][ic]][0] << endl;
	double temp_V_NS_massive_pair = gamma_q / C_F * log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) + pi2_6 - gsl_sf_dilog(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) - 2 * log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) - temp_mass2 / s_ik[ip_0][ip_1] * log(temp_mass2 / Q2_ik[ip_0][ip_1]);
	if (pair_massive[i_p] == 1 && pair_code_emitter[i_p] != 2){V_NS_massive_full[ip_0][ip_1] = temp_V_NS_massive_pair;}
	// massive emitter ip_0, massless spectator ip_1
	if (pair_code_emitter[i_p] != 2){
	/*
	  if (csi->type_parton[0][ip_0] == 0){
	    finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_g + gamma_g(N_f) + K_g(N_f)) / C_A;
	  }
	  else if (csi->type_parton[0][ip_0] != 0){
	    Gamma_q[ip_0][0] = C_F * (VA_DeltaIR1 + .5 * log(pmass2[ip_0] / s_ik[ip_0][ip_1]) - 2.);
	    finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_q[ip_0][0] + gamma_q + K_q) / C_F;
	  }
	*/
      	}
	if (pair_massive[i_p] == 2 && pair_code_emitter[i_p] != 1){V_NS_massive_full[ip_1][ip_0] = temp_V_NS_massive_pair;}
	// massive emitter ip_1, massless spectator ip_0
	if (pair_code_emitter[i_p] != 1){
	/*
	  if (csi->type_parton[0][ip_1] == 0){
	    finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_g + gamma_g(N_f) + K_g(N_f)) / C_A;
	  }
	  else if (csi->type_parton[0][ip_1] != 0){
	    Gamma_q[ip_1][0] = C_F * (VA_DeltaIR1 + .5 * log(pmass2[ip_1] / s_ik[ip_0][ip_1]) - 2.);
	    finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_q[ip_1][0] + gamma_q + K_q) / C_F;
	  }
	*/
	}
      }
      if ((pair_massive[i_p] == 1 && pair_code_emitter[i_p] != 1) || (pair_massive[i_p] == 2 && pair_code_emitter[i_p] != 2)){
	// m_i == 0 (gluon/quark) and m_k > 0 (quark)
	//	  logger << LOG_DEBUG_VERBOSE << "m_i == 0 (gluon/quark) and m_k > 0 (quark):   i = " << dx_pa[ioperator_pair[ia][ib]][0] << "   k = " << dx_pa[ioperator_pair[ia][ic]][0] << endl;
	double temp_mass;
	double temp_type;
	if (pair_massive[i_p] == 1){temp_mass = pmass[ip_0]; temp_type = csi->type_parton[0][ip_1];}
	else if (pair_massive[i_p] == 2){temp_mass = pmass[ip_1]; temp_type = csi->type_parton[0][ip_0];}
	sqrtQ2_ik[ip_0][ip_1] = sqrt(Q2_ik[ip_0][ip_1]);

	double temp_V_NS_massive_pair;
	if (temp_type == 0){
	  // m_i == 0 (gluon) and m_k > 0 (quark)
	  temp_V_NS_massive_pair = gamma_g(N_f) / C_A * (log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) - 2 * log((sqrtQ2_ik[ip_0][ip_1] - temp_mass) / sqrtQ2_ik[ip_0][ip_1]) - 2 * temp_mass / (sqrtQ2_ik[ip_0][ip_1] + temp_mass)) + pi2_6 - gsl_sf_dilog(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]); // ...
	  if ((pair_massive[i_p] == 1 && ip_1 > 2) || (pair_massive[i_p] == 2 && ip_0 > 2)){
	    // meaning of '7' has to be re-considered... replace by heaviest quark in real g->QQ~ splittings?
	    for (int i = 1; i < 7; i++){
	      if (M[i] > 0.){
		temp_V_NS_massive_pair += 2. / 3. * T_R / C_A * log(M2[i] / Q2_aux); // seems to be implemented only for kappa = 2/3 ??? !!!
		if (s_ik[ip_0][ip_1] > 4 * M[i] * (temp_mass + M[i])){
		  double rho1 = sqrt(1. - 4 * M2[i] / pow(sqrtQ2_ik[ip_0][ip_1] - temp_mass, 2));
		  temp_V_NS_massive_pair += 4. / 3. * T_R / C_A * (log((sqrtQ2_ik[ip_0][ip_1] - temp_mass) / sqrtQ2_ik[ip_0][ip_1]) + temp_mass * pow(rho1, 3) / (sqrtQ2_ik[ip_0][ip_1] + temp_mass) + log((1. + rho1) / 2.) - rho1 / 3. * (3. + pow(rho1, 2)) - .5 * log(M2[i] / Q2_ik[ip_0][ip_1]));
		}
	      }
	    }
	  }
	}
	else if (temp_type != 0){
	  // m_i == 0 (quark) and m_k > 0 (quark)
	  temp_V_NS_massive_pair = gamma_q / C_F * (log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) - 2 * log((sqrtQ2_ik[ip_0][ip_1] - temp_mass) / sqrtQ2_ik[ip_0][ip_1]) - 2 * temp_mass / (sqrtQ2_ik[ip_0][ip_1] + temp_mass)) + pi2_6 - gsl_sf_dilog(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]);
	}
	if (pair_massive[i_p] == 1 && pair_code_emitter[i_p] != 1){V_NS_massive_full[ip_1][ip_0] = temp_V_NS_massive_pair;} // massive emitter ip_1, massless spectator ip_0
	if (pair_massive[i_p] == 2 && pair_code_emitter[i_p] != 2){V_NS_massive_full[ip_0][ip_1] = temp_V_NS_massive_pair;} // massive emitter ip_0, massless spectator ip_1
      }
    }
    else if (pair_massive[i_p] == 0){
      // m_i == 0 (gluon/quark) and m_k == 0 (gluon/quark)
      //	  logger << LOG_DEBUG_VERBOSE << "m_i == 0 (gluon/quark) and m_k == 0 (gluon/quark):   i = " << dx_pa[ioperator_pair[ia][ib]][0] << "   k = " << dx_pa[ioperator_pair[ia][ic]][0] << endl;
      if ((csi->type_parton[0][ip_0] == 0 && ip_0 > 2 && pair_code_emitter[i_p] != 2) || (csi->type_parton[0][ip_1] == 0 && ip_1 > 2 && pair_code_emitter[i_p] != 1)){
	double temp_V_NS_massive = 0.; // ...
	//	double m_k = 0.; // ???
	for (int i = 1; i < 7; i++){
	  if (M[i] > 0.){
	    temp_V_NS_massive += 2. / 3. * T_R / C_A * log(M2[i] / Q2_aux);
	    //		logger << LOG_DEBUG_VERBOSE << "V_NS_massive(0., 0.)[" << i << ", N_F] = " << V_NS_massive << endl;
	    /*
	    // this must be commented out to reproduce the massless CS version !!!
	    if (s_ik[ip_0][ip_1] > 4 * M[i] * (m_k + M[i])){
	      double rho1 = sqrt(1. - 4 * M2[i] / pow(sqrtQ2_ik[ip_0][ip_1] - m_k, 2));
	      temp_V_NS_massive += 4. / 3. * T_R / C_A * (log((1. + rho1) / 2.) - rho1 / 3. * (3. + pow(rho1, 2)) - .5 * log(M2[i] / s_ik[ip_0][ip_1]));
	      //		  logger << LOG_DEBUG_VERBOSE << "V_NS_massive(0., 0.)[" << i << ", N_F^ik] = " << V_NS_massive << endl;
	    }
	    */
	  }
	}
	if (csi->type_parton[0][ip_0] == 0 && ip_0 > 2 && pair_code_emitter[i_p] != 2){V_NS_massive_full[ip_0][ip_1] = temp_V_NS_massive;}
	if (csi->type_parton[0][ip_1] == 0 && ip_1 > 2 && pair_code_emitter[i_p] != 1){V_NS_massive_full[ip_1][ip_0] = temp_V_NS_massive;}
      }
      if (csi->type_parton[0][ip_0] == 0 && pair_code_emitter[i_p] != 2){
	/*
	if (csi->type_parton[0][ip_0] == 0){
	  finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_g + gamma_g(N_f) + K_g(N_f)) / C_A;
	}
	else if (csi->type_parton[0][ip_0] != 0){
	  Gamma_q[ip_0][0] = C_F * (VA_DeltaIR1 + .5 * log(pmass2[ip_0] / s_ik[ip_0][ip_1]) - 2.);
	  finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_q[ip_0][0] + gamma_q + K_q) / C_F;
	}
	*/
      }

      if (csi->type_parton[0][ip_1] == 0 && pair_code_emitter[i_p] != 1){
	/*
	if (csi->type_parton[0][ip_1] == 0){
	  finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_g + gamma_g(N_f) + K_g(N_f)) / C_A;
	}
	else if (csi->type_parton[0][ip_1] != 0){
	  Gamma_q[ip_1][0] = C_F * (VA_DeltaIR1 + .5 * log(pmass2[ip_1] / s_ik[ip_0][ip_1]) - 2.);
	  finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_q[ip_1][0] + gamma_q + K_q) / C_F;
	}
	*/
      }
    }
  }

  for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
    if ((*VA_ioperator)[i_a][0].type_correction() == 1){
      for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
	int ip_0 = (*VA_ioperator)[i_a][j_a].pair()[0];
	int ip_1 = (*VA_ioperator)[i_a][j_a].pair()[1];
	int no_em = (*VA_ioperator)[i_a][j_a].no_emitter();
	int no_sp = (*VA_ioperator)[i_a][j_a].no_spectator();
	/*      
		prefactor_log[0] = V_S_massive_full[no_em][no_sp][0] + V_NS_massive_full[no_em][no_sp] + finite_massive_full[no_em][no_sp][0];
		prefactor_log[1] = V_S_massive_full[no_em][no_sp][1] + finite_massive_full[no_em][no_sp][1];
		prefactor_log[2] = V_S_massive_full[no_em][no_sp][2];
		//      prefactor_I_ME2_cs[i_a][j_a] = prefactor_log[0] + prefactor_log[1] * ioperator_pair_log[ia] + prefactor_log[2] * ioperator_pair_log2_2[ia];
		*/
	VA_I_ME2_cf[i_a][j_a] = -alpha_S * inv2pi * ((V_S_massive_full[no_em][no_sp][0] + V_NS_massive_full[no_em][no_sp] + finite_massive_full[no_em][no_sp][0]) 
						     + (V_S_massive_full[no_em][no_sp][1] + finite_massive_full[no_em][no_sp][1]) * ioperator_pair_log[ip_0][ip_1] 
						     + (V_S_massive_full[no_em][no_sp][2]) * ioperator_pair_log2_2[ip_0][ip_1]) * VA_ME2_cf[i_a][j_a];
	/*	
		logger << LOG_DEBUG_VERBOSE << (*VA_ioperator)[i_a][j_a].name() << endl;
		logger << LOG_DEBUG_VERBOSE << "VA_ME2_cf[" << i_a << "][" << j_a << "] = " << VA_ME2_cf[i_a][j_a] << endl;
		logger << LOG_DEBUG_VERBOSE << "    prefactor_log[" << no_em << "][" << no_sp << "][0]    = " << (V_S_massive_full[no_em][no_sp][0] + V_NS_massive_full[no_em][no_sp] + finite_massive_full[no_em][no_sp][0]) << endl;
		logger << LOG_DEBUG_VERBOSE << "    prefactor_log[" << no_em << "][" << no_sp << "][1]   = " << (V_S_massive_full[no_em][no_sp][1] + finite_massive_full[no_em][no_sp][1]) << endl;
		logger << LOG_DEBUG_VERBOSE << "    prefactor_log[" << no_em << "][" << no_sp << "][2]   = " << (V_S_massive_full[no_em][no_sp][2]) << endl;
		logger << LOG_DEBUG_VERBOSE << "   V_S_massive_full[" << no_em << "][" << no_sp << "][0] = " << V_S_massive_full[no_em][no_sp][0] << endl;
		logger << LOG_DEBUG_VERBOSE << "   V_S_massive_full[" << no_em << "][" << no_sp << "][1] = " << V_S_massive_full[no_em][no_sp][1] << endl;
		logger << LOG_DEBUG_VERBOSE << "   V_S_massive_full[" << no_em << "][" << no_sp << "][2] = " << V_S_massive_full[no_em][no_sp][2] << endl;
		logger << LOG_DEBUG_VERBOSE << "  V_NS_massive_full[" << no_em << "][" << no_sp << "]    = " << V_NS_massive_full[no_em][no_sp] << endl;
		logger << LOG_DEBUG_VERBOSE << "finite_massive_full[" << no_em << "][" << no_sp << "][0] = " << finite_massive_full[no_em][no_sp][0] << endl;
		logger << LOG_DEBUG_VERBOSE << "finite_massive_full[" << no_em << "][" << no_sp << "][1] = " << finite_massive_full[no_em][no_sp][1] << endl;
		logger << LOG_DEBUG_VERBOSE << "   ioperator_pair_log[" << ip_0 << "][" << ip_1 << "]    = " << ioperator_pair_log[ip_0][ip_1] << endl;
		logger << LOG_DEBUG_VERBOSE << "ioperator_pair_log2_2[" << ip_0 << "][" << ip_1 << "]    = " << ioperator_pair_log2_2[ip_0][ip_1] << endl;
		logger << LOG_DEBUG_VERBOSE << "VA_I_ME2_cf[" << i_a << "][" << j_a << "] = " << VA_I_ME2_cf[i_a][j_a] << endl;
*/
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::determine_ioperator_QEW(phasespace_set & psi){
  Logger logger("observable_set::determine_ioperator_QEW");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static map<int, double> charge_particle;
  if (initialization){
    fill_charge_particle(charge_particle);
    initialization = 0;
  }
  logger << LOG_DEBUG << "I-operator determination begin" << endl;
  vector<string> pa_name(csi->type_parton[0].size(), "");
  //  if (charge_particle[csi->type_parton[0][1]] != 0 || csi->type_parton[0][1] == 22){pa_name[1] = "a";}
  //  if (charge_particle[csi->type_parton[0][2]] != 0 || csi->type_parton[0][2] == 22){pa_name[2] = "b";}
  pa_name[1] = "a";
  pa_name[2] = "b";
  int count = 0;
  vector<string> alphabet(csi->type_parton[0].size() - 3, "");
  for (int i_p = 0; i_p < alphabet.size(); i_p++){alphabet[i_p] = char(105 + i_p);}
  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (charge_particle[csi->type_parton[0][i_p]] != 0 || csi->type_parton[0][i_p] == 22){pa_name[i_p] = alphabet[count++];}}

  int temp_type_correction = 2;

  for (int temp_no_emitter = 1; temp_no_emitter < csi->type_parton[0].size(); temp_no_emitter++){
    for (int temp_no_spectator = 1; temp_no_spectator < csi->type_parton[0].size(); temp_no_spectator++){
      if (temp_no_emitter == temp_no_spectator){continue;}
      if (pa_name[temp_no_emitter] == "" || pa_name[temp_no_spectator] == ""){continue;}
      vector<int> temp_pair(2);
      temp_pair[0] = std::min(temp_no_emitter, temp_no_spectator);
      temp_pair[1] = std::max(temp_no_emitter, temp_no_spectator);

      double temp_charge_factor = 1.;
      if (csi->type_parton[0][temp_no_emitter] != 22 && csi->type_parton[0][temp_no_spectator] != 22){
	temp_charge_factor = charge_particle[abs(csi->type_parton[0][temp_no_spectator])] / charge_particle[abs(csi->type_parton[0][temp_no_emitter])]; // Q²_ij / (Q_ij Q_k)
	//	temp_charge_factor = charge_particle[abs(csi->type_parton[0][temp_no_emitter])] * charge_particle[abs(csi->type_parton[0][temp_no_spectator])];
	if ((temp_no_emitter > 2 && csi->type_parton[0][temp_no_emitter] > 0) || (temp_no_emitter < 3 && csi->type_parton[0][temp_no_emitter] < 0)){temp_charge_factor = -temp_charge_factor;}
	if ((temp_no_spectator > 2 && csi->type_parton[0][temp_no_spectator] > 0) || (temp_no_spectator < 3 && csi->type_parton[0][temp_no_spectator] < 0)){temp_charge_factor = -temp_charge_factor;}
      }
      else if (csi->type_parton[0][temp_no_emitter] == 22){
	// For photon emitters, the following conventions are used:
	// initial-state emitter: kappa_ij,k = -1 for the other initial-state particle, 0 elsewhere
	// final-state emitter:   kappa_ij,k = -.5 for the two initial-state particles, 0 elsewhere
	if (temp_no_emitter < 3){
	  if (temp_no_spectator < 3){temp_charge_factor = -1.;} // kappa_ij,k = -1 for the other initial-state particle, 0 elsewhere
	  else {temp_charge_factor = 0.;}
	}
	else {
	  // temporary !!! deactivates photon-emitters !!! (not always correct !!!)
	  if (temp_no_spectator < 3){temp_charge_factor = -0.5;} // kappa_ij,k = -1 for the other initial-state particle, 0 elsewhere
	  else {temp_charge_factor = 0.;}
	}
	/*
	// old version (always inactive, see below):
	if (temp_no_spectator < 3){temp_charge_factor = -1.;} // kappa_ij,k = -1 for the other initial-state particle, 0 elsewhere
	else {temp_charge_factor = 0.;}
	*/
	// temporary !!! deactivates photon-emitters !!! (not always correct !!!)
	/////	temp_charge_factor = 0.;
      }
      else {temp_charge_factor = 0.;}
      
      int temp_type;
      if (csi->type_parton[0][temp_no_emitter] == 22){temp_type = 0;}
      else {temp_type = 1;}

      int temp_massive;
      if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 0;}
      else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 1;}
      else if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 2;}
      else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 3;}
      else {logger << LOG_WARN << "Should not happen!" << endl;}
      /*
      int temp_massive;
      if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 0;}
      else {temp_massive = 1;}
      */
      string temp_name;
      temp_name = "I^{" + pa_name[temp_no_emitter] + "," + pa_name[temp_no_spectator] + "}";

      int flag = (*VA_ioperator).size();
      for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
	if (temp_no_emitter == (*VA_ioperator)[i_a][0].no_emitter() && temp_type_correction == (*VA_ioperator)[i_a][0].type_correction()){
	  flag = i_a; 
	  break;
	}
      }
      //      int flag = (*VA_ioperator).size();
      //      for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){if (temp_no_emitter == (*VA_ioperator)[i_a][0].no_emitter()){flag = i_a; break;}}
      if (flag == (*VA_ioperator).size()){(*VA_ioperator).push_back(vector<ioperator_set> ());}
      //      flag = (*VA_ioperator).size() - 1;
      //      }
      (*VA_ioperator)[flag].push_back(ioperator_set(temp_name, temp_type, temp_pair, psi_no_prc[0], csi->type_parton[0], temp_charge_factor, temp_no_emitter, temp_no_spectator, temp_type_correction, temp_massive));
    }
  }
  logger << LOG_INFO << "Before selection of contributing I-operator 'dipoles':" << endl;
  logger.newLine(LOG_INFO);
  for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
    for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
      logger << LOG_INFO << "(*VA_ioperator)[" << i_a << "][" << j_a << "].name() = " << (*VA_ioperator)[i_a][j_a].name() << endl;
    }
  }

  // check if VA_ioperator-type gives a non-vanishing contribution:
  for (int i_a = (*VA_ioperator).size() - 1; i_a >= 0; i_a--){
    for (int j_a = (*VA_ioperator)[i_a].size() - 1; j_a >= 0; j_a--){
      // old version: contributions from photon splittings were deactivated in the I-operator (....type() == 0):
      //      if ((*VA_ioperator)[i_a][j_a].type() == 0 || (*VA_ioperator)[i_a][j_a].charge_factor() == 0.){
      // new veresion: arrange photon splittings via charge_factor of the respective I-operator contribution:
      if ((*VA_ioperator)[i_a][j_a].charge_factor() == 0.){
	(*VA_ioperator)[i_a].erase((*VA_ioperator)[i_a].begin() + j_a);
      }
    }
    if ((*VA_ioperator)[i_a].size() == 0){
      (*VA_ioperator).erase((*VA_ioperator).begin() + i_a);
    }
  }
  // pdf-selection effects are included later !!!

  logger << LOG_INFO << "After selection of contributing I-operator 'dipoles':" << endl;
  logger.newLine(LOG_INFO);
  for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
    for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
      logger << LOG_INFO << "(*VA_ioperator)[" << i_a << "][" << j_a << "].name() = " << (*VA_ioperator)[i_a][j_a].name() << endl;
    }
  }



  // Switch off selected ioperator contributions:

  for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
    logger << LOG_DEBUG << "kill:   i_a = " << i_a << endl;
    for (int j_a = (*VA_ioperator)[i_a].size() - 1; j_a >= 0 ; j_a--){
      //      logger << LOG_DEBUG << "kill:   i_a = " << i_a << "   j_a = " << j_a << endl;
      if (user.string_value[user.string_map["selection"]] == "ii"){
	// only initial-initial contribution:
	if ((*VA_ioperator)[i_a][j_a].no_emitter() > 2 || (*VA_ioperator)[i_a][j_a].no_spectator() > 2){
	  (*VA_ioperator)[i_a].erase((*VA_ioperator)[i_a].begin() + j_a);
	  logger << LOG_DEBUG << "kill:   (*VA_ioperator)[" << i_a << "][" << j_a << "].name() = " << (*VA_ioperator)[i_a][j_a].name() << endl;
	}
      }
      else if (user.string_value[user.string_map["selection"]] == "if"){
	// only initial-final contribution:
	if ((*VA_ioperator)[i_a][j_a].no_emitter() > 2 || (*VA_ioperator)[i_a][j_a].no_spectator() < 3){
	  (*VA_ioperator)[i_a].erase((*VA_ioperator)[i_a].begin() + j_a);
	  logger << LOG_DEBUG << "kill:   (*VA_ioperator)[" << i_a << "][" << j_a << "].name() = " << (*VA_ioperator)[i_a][j_a].name() << endl;
	}
      }
      else if (user.string_value[user.string_map["selection"]] == "fi"){
	// only final-initial contribution:
	if ((*VA_ioperator)[i_a][j_a].no_emitter() < 3 || (*VA_ioperator)[i_a][j_a].no_spectator() > 2){
	  (*VA_ioperator)[i_a].erase((*VA_ioperator)[i_a].begin() + j_a);
	  logger << LOG_DEBUG << "kill:   (*VA_ioperator)[" << i_a << "][" << j_a << "].name() = " << (*VA_ioperator)[i_a][j_a].name() << endl;
	}
      }
      else if (user.string_value[user.string_map["selection"]] == "ff"){
	// only final-final contribution:
	if ((*VA_ioperator)[i_a][j_a].no_emitter() < 3 || (*VA_ioperator)[i_a][j_a].no_spectator() < 3){
	  (*VA_ioperator)[i_a].erase((*VA_ioperator)[i_a].begin() + j_a);
	  logger << LOG_DEBUG << "kill:   (*VA_ioperator)[" << i_a << "][" << j_a << "].name() = " << (*VA_ioperator)[i_a][j_a].name() << endl;
	}
      }
    	
    }
  }


  for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
    for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
      logger << LOG_DEBUG << "(*VA_ioperator)[" << i_a << "][" << j_a << "].name() = " << (*VA_ioperator)[i_a][j_a].name() << endl;
    }
  }

  

  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::calculate_ioperator_QEW(){
  Logger logger("observable_set::calculate_ioperator_QEW");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (massive_QEW){calculate_ioperator_QEW_CDST();}
  else {calculate_ioperator_QEW_CS();}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_ioperator_QEW_CS(){
  Logger logger("observable_set::calculate_ioperator_QEW_CS");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static vector<vector<int> > pair;
  static vector<vector<double> > ioperator_pair_log(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > ioperator_pair_log2_2(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<int> type;
  static double prefactor_log[2][2];
  static double alpha_e_2pi = msi.alpha_e * inv2pi;
  static vector<vector<int> > type_emitter((*VA_ioperator).size());
  static vector<vector<double> > charge_factor((*VA_ioperator).size());
  static vector<vector<vector<int> > > ppair((*VA_ioperator).size());
  static map<int, double> charge_particle;
  static vector<double> charge2_particle(26);

  //  if (VA_delta_flag){initialization = 1;}

  if (initialization){
    fill_charge_particle(charge_particle);
    for (int i_p = 0; i_p < 26; i_p++){
      charge2_particle[i_p] = pow(charge_particle[i_p], 2);
    }
    for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
      if ((*VA_ioperator)[i_a][0].type_correction() == 2){
	for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
	  int flag = -1;
	  for (int i_p = 0; i_p < pair.size(); i_p++){if ((*VA_ioperator)[i_a][j_a].pair() == pair[i_p]){flag = i_p; break;}}
	  if (flag == -1){pair.push_back((*VA_ioperator)[i_a][j_a].pair());}
	  flag = -1;
	  for (int i_t = 0; i_t < type.size(); i_t++){if ((*VA_ioperator)[i_a][j_a].type() == type[i_t]){flag = i_t; break;}}
	  if (flag == -1){type.push_back((*VA_ioperator)[i_a][j_a].type());}
	}
      }
    }
    for (int i_t = 0; i_t < type.size(); i_t++){
      if (switch_VI_bosonic_fermionic == 0){
	// both fermionic and bosonic contributions
	if (type[i_t] == 0){
	  prefactor_log[0][0] = (VA_DeltaIR1 + 1.) * gamma_a(N_f) + K_a(N_f); // T²_i = C_A -> 0
	  prefactor_log[0][1] = gamma_a(N_f);
	}
	else if (type[i_t] == 1){
	  prefactor_log[1][0] = VA_DeltaIR2 - pi2_3 + ((VA_DeltaIR1 + 1.) * gamma_qew_q + K_qew_q);
	  if (switch_polenorm == 1){prefactor_log[1][0] += -pi2_6;} // BLHA -> COLI !!!
	  prefactor_log[1][1] = VA_DeltaIR1 + gamma_qew_q;
	}
      }
      else if (switch_VI_bosonic_fermionic == 1){
	// only bosonic contributions
	if (type[i_t] == 0){
	  prefactor_log[0][0] = 0.;
	  prefactor_log[0][1] = 0.;
	}
	else if (type[i_t] == 1){
	  prefactor_log[1][0] = VA_DeltaIR2 - pi2_3 + ((VA_DeltaIR1 + 1.) * gamma_qew_q + K_qew_q);
	  if (switch_polenorm == 1){prefactor_log[1][0] += -pi2_6;} // BLHA -> COLI !!!
	  prefactor_log[1][1] = VA_DeltaIR1 + gamma_qew_q;
	}
      }
      else if (switch_VI_bosonic_fermionic == 2){
	// only fermionic contribution
	if (type[i_t] == 0){
	  prefactor_log[0][0] = ((VA_DeltaIR1 + 1.) * gamma_a(N_f) + K_a(N_f));
	  prefactor_log[0][1] = gamma_a(N_f);
	}
	else if (type[i_t] == 1){
	  prefactor_log[1][0] = 0.;
	  prefactor_log[1][1] = 0.;
	}
      }
    }

    for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
      if ((*VA_ioperator)[i_a][0].type_correction() == 2){
	type_emitter[i_a].resize((*VA_ioperator)[i_a].size());
	charge_factor[i_a].resize((*VA_ioperator)[i_a].size());
	ppair[i_a].resize((*VA_ioperator)[i_a].size());
	for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
	  type_emitter[i_a][j_a] = (*VA_ioperator)[i_a][j_a].type();
	  charge_factor[i_a][j_a] = (*VA_ioperator)[i_a][j_a].charge_factor();
	  ppair[i_a][j_a] = (*VA_ioperator)[i_a][j_a].pair();
	}
      }
    }

    for (int i_p = 0; i_p < pair.size(); i_p++){logger << LOG_DEBUG_VERBOSE << "pair[" << i_p << "] = " << "(" << pair[i_p][0] << ", " << pair[i_p][1] << ")" << endl;}
    for (int i_t = 0; i_t < type.size(); i_t++){logger << LOG_DEBUG_VERBOSE << "type[" << i_t << "] = " << type[i_t] << endl;}
    //    initialization = 0;
    //    if (initialization){initialization = 0;}
    if (!VA_delta_flag){initialization = 0;}
  }

  // check if this causes the problem with I-operator: !!!
  double mu2IR = pow(var_mu_ren, 2);
  // modified to deal with mu_ren != mu_reg !!!
  if (user.double_value[user.double_map["mu_reg"]] != 0. && user.double_value[user.double_map["mu_reg"]] != -1.){mu2IR = pow(user.double_value[user.double_map["mu_reg"]], 2);}

  for (int i_p = 0; i_p < pair.size(); i_p++){
    ioperator_pair_log[pair[i_p][0]][pair[i_p][1]] = log(mu2IR / (2. * p_parton[0][pair[i_p][0]] * p_parton[0][pair[i_p][1]]));
    ioperator_pair_log2_2[pair[i_p][0]][pair[i_p][1]] = pow(ioperator_pair_log[pair[i_p][0]][pair[i_p][1]], 2) / 2.;
  }

  for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
    if ((*VA_ioperator)[i_a][0].type_correction() == 2){
      for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
	//	VA_I_ME2_cf[i_a][j_a] = -alpha_e_2pi * charge2_particle[abs(csi->type_parton[0][(*VA_ioperator)[i_a][j_a].no_emitter()])] * (prefactor_log[type_emitter[i_a][j_a]][0] + prefactor_log[type_emitter[i_a][j_a]][1] * ioperator_pair_log[ppair[i_a][j_a][0]][ppair[i_a][j_a][1]] + ioperator_pair_log2_2[ppair[i_a][j_a][0]][ppair[i_a][j_a][1]]) * VA_ME2_cf[i_a][j_a];
	
	if (csi->type_parton[0][(*VA_ioperator)[i_a][j_a].no_emitter()] == 22){ // check why the charge factor is 1 for external photons !!!
	  VA_I_ME2_cf[i_a][j_a] = -alpha_e_2pi * (prefactor_log[type_emitter[i_a][j_a]][0] + prefactor_log[type_emitter[i_a][j_a]][1] * ioperator_pair_log[ppair[i_a][j_a][0]][ppair[i_a][j_a][1]] + ioperator_pair_log2_2[ppair[i_a][j_a][0]][ppair[i_a][j_a][1]]) * VA_ME2_cf[i_a][j_a];
	}
	else { // What about charge_factor[i_a][j_a] ???
	  // original implementation:
	  // logger << LOG_DEBUG_VERBOSE << "charge2_particle[abs(csi->type_parton[0][(*VA_ioperator)[" << i_a << "][" << j_a << "].no_emitter()])] = " << charge2_particle[abs(csi->type_parton[0][(*VA_ioperator)[i_a][j_a].no_emitter()])] << endl;
	  VA_I_ME2_cf[i_a][j_a] = -alpha_e_2pi * charge2_particle[abs(csi->type_parton[0][(*VA_ioperator)[i_a][j_a].no_emitter()])] * (prefactor_log[type_emitter[i_a][j_a]][0] + prefactor_log[type_emitter[i_a][j_a]][1] * ioperator_pair_log[ppair[i_a][j_a][0]][ppair[i_a][j_a][1]] + ioperator_pair_log2_2[ppair[i_a][j_a][0]][ppair[i_a][j_a][1]]) * VA_ME2_cf[i_a][j_a];
	}
	
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




void observable_set::calculate_ioperator_QEW_CDST(){
  Logger logger("observable_set::calculate_ioperator_QEW_CDST");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static vector<vector<int> > pair;
  static vector<int> pair_massive;
  static vector<int> pair_code_emitter;
  //  static double prefactor_log[2][2];
  static double alpha_e_2pi = msi.alpha_e * inv2pi;
  //  static double alpha_S_2pi = alpha_S * inv2pi;

  //  static int N_f = osi_N_f;

  static vector<vector<int> > type_emitter((*VA_ioperator).size());
  static vector<vector<int> > no_emitter((*VA_ioperator).size());
  static vector<vector<int> > no_spectator((*VA_ioperator).size());
  static vector<vector<vector<int> > > ppair((*VA_ioperator).size());

  static vector<double> pmass(csi->type_parton[0].size());
  static vector<double> pmass2(csi->type_parton[0].size());

  static vector<vector<double> > s_ik(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > Q2_ik(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > sqrtQ2_ik(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > mu2_i(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > rho2_i(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > v_ik(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > delta(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));

  static vector<vector<double> > V_NS_massive_full(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<vector<double> > > V_S_massive_full(csi->type_parton[0].size(), vector<vector<double> > (csi->type_parton[0].size(), vector<double> (3)));
  static vector<vector<vector<double> > > finite_massive_full(csi->type_parton[0].size(), vector<vector<double> > (csi->type_parton[0].size(), vector<double> (2, 0.)));
  static vector<vector<double> > ioperator_pair_log(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));
  static vector<vector<double> > ioperator_pair_log2_2(csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size()));

  static double Gamma_a;
  static vector<vector<double> > Gamma_q(csi->type_parton[0].size(), vector<double> (2));

  if (initialization){
    for (int i_p = 0; i_p < csi->type_parton[0].size(); i_p++){
      pmass[i_p] = mass_parton[0][i_p];
      pmass2[i_p] = mass2_parton[0][i_p];
      if (pmass[i_p] == 0.){
	Gamma_q[i_p][0] = gamma_qew_q * VA_DeltaIR1;
	Gamma_q[i_p][1] = 0.;
      }
      else if (pmass[i_p] != 0.){
	//    Gamma_q[i_p][0] = (VA_DeltaIR1 + .5 * log(m2_i / s_ik) - 2.);
	Gamma_q[i_p][1] = -.5; // gamma_q = 1.5 !!! check if nothing goes wrong here !!! (explicit mu dependence)
      }
    }
    // Gamma_a not treated here -> depends on Q2aux

    for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
      if ((*VA_ioperator)[i_a][0].type_correction() == 2){
	for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
	  int flag = -1;
	  for (int i_p = 0; i_p < pair.size(); i_p++){if ((*VA_ioperator)[i_a][j_a].pair() == pair[i_p]){flag = i_p; break;}}
	  // !!! just modified !!! 
	  if (flag == -1){
	    pair.push_back((*VA_ioperator)[i_a][j_a].pair());
	    pair_massive.push_back((*VA_ioperator)[i_a][j_a].massive());
	    int temp_code_emitter;
	    if ((*VA_ioperator)[i_a][j_a].pair()[0] == (*VA_ioperator)[i_a][j_a].no_emitter()){temp_code_emitter = 1;}
	    else if ((*VA_ioperator)[i_a][j_a].pair()[1] == (*VA_ioperator)[i_a][j_a].no_emitter()){temp_code_emitter = 2;}
	    pair_code_emitter.push_back(temp_code_emitter);
	  }
	  else {pair_code_emitter[flag] = 3;}
	  /*
	    flag = 0;
	    for (int i_t = 0; i_t < type.size(); i_t++){if ((*VA_ioperator)[i_a][j_a].type() == type[i_t]){flag = 1; break;}}
	    if (flag == 0){type.push_back((*VA_ioperator)[i_a][j_a].type());}
	  */
	}
      }
    }

    for (int i_p = 0; i_p < pair.size(); i_p++){
      int ip_0 = pair[i_p][0];
      int ip_1 = pair[i_p][1];
      if (pair_massive[i_p] == 3){
	if (pair_code_emitter[i_p] != 2){V_S_massive_full[ip_0][ip_1][2] = 0.;}
	if (pair_code_emitter[i_p] != 1){V_S_massive_full[ip_1][ip_0][2] = 0.;}
      }
      else if (pair_massive[i_p] == 1 || pair_massive[i_p] == 2){
	if (pair_code_emitter[i_p] != 2){V_S_massive_full[ip_0][ip_1][2] = .5;}
	if (pair_code_emitter[i_p] != 1){V_S_massive_full[ip_1][ip_0][2] = .5;}
      }
      else if (pair_massive[i_p] == 0){
	// m_i == 0 (gluon/quark) and m_k == 0 (gluon/quark)
	if (pair_code_emitter[i_p] != 2){
	  V_S_massive_full[ip_0][ip_1][0] = VA_DeltaIR2;
	  if (switch_polenorm == 1){V_S_massive_full[ip_0][ip_1][0] += -pi2_6;} // BLHA -> COLI !!!
	  V_S_massive_full[ip_0][ip_1][1] = VA_DeltaIR1;
	  V_S_massive_full[ip_0][ip_1][2] = 1.;
	  if (csi->type_parton[0][ip_0] == 0 && ip_0 < 3){V_NS_massive_full[ip_0][ip_1] = 0.;}
	  else if (csi->type_parton[0][ip_0] != 0){V_NS_massive_full[ip_0][ip_1] = 0.;}
	}
	if (pair_code_emitter[i_p] != 1){
	  V_S_massive_full[ip_1][ip_0][0] = VA_DeltaIR2;
	  if (switch_polenorm == 1){V_S_massive_full[ip_1][ip_0][0] += -pi2_6;} // BLHA -> COLI !!!
	  V_S_massive_full[ip_1][ip_0][1] = VA_DeltaIR1;
	  V_S_massive_full[ip_1][ip_0][2] = 1.;
	  if (csi->type_parton[0][ip_1] == 0 && ip_1 < 3){V_NS_massive_full[ip_1][ip_0] = 0.;}
	  if (csi->type_parton[0][ip_1] != 0){V_NS_massive_full[ip_1][ip_0] = 0.;}
	}
      }


      ///      for (int i_p = 0; i_p < pair.size(); i_p++){
      logger << LOG_DEBUG << "	pair[i_p] = " << pair[i_p][0] << "  " << pair[i_p][1] << "   pair_massive[i_p] = " << pair_massive[i_p] << "   pair_code_emitter[" <<  i_p << "] = " << pair_code_emitter[i_p] << endl;
      ///      }

      
      if (pair_code_emitter[i_p] != 2){
	if (csi->type_parton[0][ip_0] == 0){
	  //	  finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_a + gamma_a(N_f) + K_a(N_f));   // Gamma_g depends on Q2_aux -> later
	  finite_massive_full[ip_0][ip_1][1] = + gamma_a(N_f);
	}
	else if (csi->type_parton[0][ip_0] != 0){
	  if (pmass[ip_0] == 0){
	    finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_q[ip_0][0] + gamma_qew_q + K_qew_q);
	    //	    logger << LOG_DEBUG_VERBOSE << "Gamma_q[" << ip_0 << "][0] = " << Gamma_q[ip_0][0] << endl;
	    //	    logger << LOG_DEBUG_VERBOSE << "gamma_qew_q = " << gamma_qew_q << endl;
	    //	    logger << LOG_DEBUG_VERBOSE << "K_qew_q = " << K_qew_q << endl;
	  }
	  //	  else {finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_q[ip_0][0] + gamma_q + K_q);}   // Gamma_q[0] depends on phase-space
	  finite_massive_full[ip_0][ip_1][1] = + (Gamma_q[ip_0][1] + gamma_qew_q);
	}
	//	logger << LOG_DEBUG_VERBOSE << "1st pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 2   finite_massive_full[" << ip_0 << "][" << ip_1 << "][0] = " << finite_massive_full[ip_0][ip_1][0] << endl;
	// 	logger << LOG_DEBUG_VERBOSE << "1st pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 2   finite_massive_full[" << ip_1 << "][" << ip_0 << "][0] = " << finite_massive_full[ip_1][ip_0][0] << endl;
      }
      if (pair_code_emitter[i_p] != 1){
	if (csi->type_parton[0][ip_1] == 0){
	  //	  finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_a + gamma_a(N_f) + K_a(N_f));   // Gamma_g depends on Q2_aux -> later
	  finite_massive_full[ip_1][ip_0][1] = + gamma_a(N_f);
	}
	else if (csi->type_parton[0][ip_1] != 0){
	  if (pmass[ip_1] == 0){finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_q[ip_1][0] + gamma_qew_q + K_qew_q);}
	  //	  else {finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_q[ip_1][0] + gamma_qew_q + K_qew_q);}   // Gamma_q[0] depends on phase-space
	  finite_massive_full[ip_1][ip_0][1] = + (Gamma_q[ip_1][1] + gamma_qew_q);
	}
	// 	logger << LOG_DEBUG_VERBOSE << "1st pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 1   finite_massive_full[" << ip_0 << "][" << ip_1 << "][0] = " << finite_massive_full[ip_0][ip_1][0] << endl;
	//	logger << LOG_DEBUG_VERBOSE << "1st pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 1   finite_massive_full[" << ip_1 << "][" << ip_0 << "][0] = " << finite_massive_full[ip_1][ip_0][0] << endl;
     }
    }
 
    for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
      if ((*VA_ioperator)[i_a][0].type_correction() == 1){
	type_emitter[i_a].resize((*VA_ioperator)[i_a].size());
	no_emitter[i_a].resize((*VA_ioperator)[i_a].size());
	no_spectator[i_a].resize((*VA_ioperator)[i_a].size());
	ppair[i_a].resize((*VA_ioperator)[i_a].size());
	for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
	  type_emitter[i_a][j_a] = (*VA_ioperator)[i_a][j_a].type();
	  no_emitter[i_a][j_a] = (*VA_ioperator)[i_a][j_a].no_emitter();
	  no_spectator[i_a][j_a] = (*VA_ioperator)[i_a][j_a].no_spectator();
	  ppair[i_a][j_a] = (*VA_ioperator)[i_a][j_a].pair();
	}
      }
    }

    //    logger << LOG_DEBUG << "pair_massive[" << i_p << "] = " << pair_massive[i_p] << endl;
    for (int i_p = 0; i_p < pair.size(); i_p++){logger << LOG_DEBUG << "pair[" << i_p << "] = " << "(" << pair[i_p][0] << ", " << pair[i_p][1] << ")   massive = " << pair_massive[i_p] << "   code_emitter = " << pair_code_emitter[i_p] << endl;}
    //    for (int i_t = 0; i_t < type.size(); i_t++){logger << LOG_DEBUG << "type[" << i_t << "] = " << type[i_t] << endl;}
    //    initialization = 0;
    //    if (initialization){initialization = 0;}
    if (!VA_delta_flag){initialization = 0;}
  }

  //  double kappa = 2. / 3.;
  double mu2IR = pow(var_mu_ren, 2);
  // modified to deal with mu_ren != mu_reg !!!
  if (user.double_value[user.double_map["mu_reg"]] != 0. && user.double_value[user.double_map["mu_reg"]] != -1.){mu2IR = pow(user.double_value[user.double_map["mu_reg"]], 2);}

  double Q2_aux = mu2IR;

  // needed only if a gluon is around !!!
  Gamma_a = gamma_a(N_f) * VA_DeltaIR1;
  for (int i_p = 1; i_p < 7; i_p++){if (M[i_p] > 0.){Gamma_a += -2./3. * N_c * log(M2[i_p] / Q2_aux);}} // T_R -> N_c




  for (int i_p = 0; i_p < pair.size(); i_p++){
    int ip_0 = pair[i_p][0];
    int ip_1 = pair[i_p][1];
    // always needed
    s_ik[ip_0][ip_1] = 2. * p_parton[0][ip_0] * p_parton[0][ip_1];
    ioperator_pair_log[ip_0][ip_1] = log(mu2IR / s_ik[ip_0][ip_1]);
    ioperator_pair_log2_2[ip_0][ip_1] = pow(ioperator_pair_log[ip_0][ip_1], 2) / 2.;
    
    
    
    if (pair_massive[i_p] == 3){
    // m_i > 0 (quark) && m_k > 0 (quark) 
      Q2_ik[ip_0][ip_1] = pmass2[ip_0] + pmass2[ip_1] + s_ik[ip_0][ip_1];
      sqrtQ2_ik[ip_0][ip_1] = sqrt(Q2_ik[ip_0][ip_1]);
      mu2_i[ip_0][ip_1] = pmass2[ip_0] / Q2_ik[ip_0][ip_1];
      mu2_i[ip_1][ip_0] = pmass2[ip_1] / Q2_ik[ip_0][ip_1];
      double temp_r2jk = pmass2[ip_0] * pmass2[ip_1] / pow(s_ik[ip_0][ip_1], 2);
      if (temp_r2jk < 1.e-3){
	delta[ip_0][ip_1] = 2 * temp_r2jk + 2 * pow(temp_r2jk, 2) + 4 * pow(temp_r2jk, 3) + 10 * pow(temp_r2jk, 4) + 28 * pow(temp_r2jk, 5);
	v_ik[ip_0][ip_1] = 1. - delta[ip_0][ip_1];
      }
      else {
	v_ik[ip_0][ip_1] = sqrt(lambda(1., mu2_i[ip_0][ip_1], mu2_i[ip_1][ip_0])) / (1. - mu2_i[ip_0][ip_1] - mu2_i[ip_1][ip_0]);
	delta[ip_0][ip_1] = 1. - v_ik[ip_0][ip_1];
      }
      double rho2 = delta[ip_0][ip_1] / (1. + v_ik[ip_0][ip_1]);
      double logrho2 = log(rho2);
      double logrho = log(sqrt(rho2));
      rho2_i[ip_0][ip_1] = (1. - v_ik[ip_0][ip_1] + 2 * mu2_i[ip_0][ip_1] / (1. - mu2_i[ip_0][ip_1] - mu2_i[ip_1][ip_0])) / (1. + v_ik[ip_0][ip_1] + 2 * mu2_i[ip_0][ip_1] / (1. - mu2_i[ip_0][ip_1] - mu2_i[ip_1][ip_0]));
      rho2_i[ip_1][ip_0] = (1. - v_ik[ip_0][ip_1] + 2 * mu2_i[ip_1][ip_0] / (1. - mu2_i[ip_0][ip_1] - mu2_i[ip_1][ip_0])) / (1. + v_ik[ip_0][ip_1] + 2 * mu2_i[ip_1][ip_0] / (1. - mu2_i[ip_0][ip_1] - mu2_i[ip_1][ip_0]));
      double temp_V_S_massive_pair_0 = 1. / v_ik[ip_0][ip_1] * (logrho * VA_DeltaIR1 -.25 * pow(log(rho2_i[ip_0][ip_1]), 2) -.25 * pow(log(rho2_i[ip_1][ip_0]), 2) - pi2_6 + logrho * log(Q2_ik[ip_0][ip_1] / s_ik[ip_0][ip_1]));
      double temp_V_S_massive_pair_1 = logrho / v_ik[ip_0][ip_1];
      double temp_V_NS_massive_pair = gamma_qew_q * log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) + 1. / v_ik[ip_0][ip_1] * (logrho2 * log(1. + rho2) + 2 * gsl_sf_dilog(rho2) - gsl_sf_dilog(1. - rho2_i[ip_0][ip_1]) - gsl_sf_dilog(1. - rho2_i[ip_1][ip_0]) - pi2_6);
      if (pair_code_emitter[i_p] != 2){
	V_S_massive_full[ip_0][ip_1][0] = temp_V_S_massive_pair_0;
	V_S_massive_full[ip_0][ip_1][1] = temp_V_S_massive_pair_1;
	//      V_S_massive_full[ip_0][ip_1][2] = 0.;
	V_NS_massive_full[ip_0][ip_1] = temp_V_NS_massive_pair + log((sqrtQ2_ik[ip_0][ip_1] - pmass[ip_1]) / sqrtQ2_ik[ip_0][ip_1]) - 2 * log((pow(sqrtQ2_ik[ip_0][ip_1] - pmass[ip_1], 2) - pmass2[ip_0]) / Q2_ik[ip_0][ip_1]) - 2 * pmass2[ip_0] / s_ik[ip_0][ip_1] * log(pmass[ip_0] / (sqrtQ2_ik[ip_0][ip_1] - pmass[ip_1])) - pmass[ip_1] / (sqrtQ2_ik[ip_0][ip_1] - pmass[ip_1]) + 2 * pmass[ip_1] * (2 * pmass[ip_1] - sqrtQ2_ik[ip_0][ip_1]) / s_ik[ip_0][ip_1] + pi2_2;
	if (csi->type_parton[0][ip_0] == 0){
	  finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_a + gamma_a(N_f) + K_a(N_f));
	}
	else if (csi->type_parton[0][ip_0] != 0){
	  Gamma_q[ip_0][0] = (VA_DeltaIR1 + .5 * log(pmass2[ip_0] / s_ik[ip_0][ip_1]) - 2.);
	  finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_q[ip_0][0] + gamma_qew_q + K_qew_q);
	}
	// 	logger << LOG_DEBUG_VERBOSE << "pair_massive[" << i_p << "] == 3   pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 2   finite_massive_full[" << ip_0 << "][" << ip_1 << "][0] = " << finite_massive_full[ip_0][ip_1][0] << endl;
	// 	logger << LOG_DEBUG_VERBOSE << "pair_massive[" << i_p << "] == 3   pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 2   finite_massive_full[" << ip_1 << "][" << ip_0 << "][0] = " << finite_massive_full[ip_1][ip_0][0] << endl;
      }
      
      if (pair_code_emitter[i_p] != 1){
	V_S_massive_full[ip_1][ip_0][0] = temp_V_S_massive_pair_0;
	V_S_massive_full[ip_1][ip_0][1] = temp_V_S_massive_pair_1;
	//      V_S_massive_full[ip_0][ip_1][2] = 0.;
	V_NS_massive_full[ip_1][ip_0] = temp_V_NS_massive_pair + log((sqrtQ2_ik[ip_0][ip_1] - pmass[ip_0]) / sqrtQ2_ik[ip_0][ip_1]) - 2 * log((pow(sqrtQ2_ik[ip_0][ip_1] - pmass[ip_0], 2) - pmass2[ip_1]) / Q2_ik[ip_0][ip_1]) - 2 * pmass2[ip_1] / s_ik[ip_0][ip_1] * log(pmass[ip_1] / (sqrtQ2_ik[ip_0][ip_1] - pmass[ip_0])) - pmass[ip_0] / (sqrtQ2_ik[ip_0][ip_1] - pmass[ip_0]) + 2 * pmass[ip_0] * (2 * pmass[ip_0] - sqrtQ2_ik[ip_0][ip_1]) / s_ik[ip_0][ip_1] + pi2_2;
	if (csi->type_parton[0][ip_1] == 0){
	  finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_a + gamma_a(N_f) + K_a(N_f));
	}
	else if (csi->type_parton[0][ip_1] != 0){
	  Gamma_q[ip_1][0] = (VA_DeltaIR1 + .5 * log(pmass2[ip_1] / s_ik[ip_0][ip_1]) - 2.);
	  finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_q[ip_0][0] + gamma_qew_q + K_qew_q);
	}
	//	logger << LOG_DEBUG_VERBOSE << "pair_massive[" << i_p << "] == 3   pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 1   finite_massive_full[" << ip_0 << "][" << ip_1 << "][0] = " << finite_massive_full[ip_0][ip_1][0] << endl;
	// 	logger << LOG_DEBUG_VERBOSE << "pair_massive[" << i_p << "] == 3   pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 1   finite_massive_full[" << ip_1 << "][" << ip_0 << "][0] = " << finite_massive_full[ip_1][ip_0][0] << endl;
      }
    }
    else if (pair_massive[i_p] == 1 || pair_massive[i_p] == 2){
      double temp_mass2;
      if (pair_massive[i_p] == 1){temp_mass2 = pmass2[ip_0];}
      else if (pair_massive[i_p] == 2){temp_mass2 = pmass2[ip_1];}
      else {logger << LOG_WARN << "Should not happen!" << endl;}

      // !!!
      if (temp_mass2 == 0. && pmass2[ip_0] != 0){temp_mass2 = pmass2[ip_0];}
      if (temp_mass2 == 0. && pmass2[ip_1] != 0){temp_mass2 = pmass2[ip_1];}
      // !!!

      Q2_ik[ip_0][ip_1] = temp_mass2 + s_ik[ip_0][ip_1];
      double temp_log_s_Q2_ik = log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]);
      logger << LOG_DEBUG << "s_ik[" << ip_0 << "][" << ip_1 << "] = " << s_ik[ip_0][ip_1] << endl;
      logger << LOG_DEBUG << "Q2_ik[" << ip_0 << "][" << ip_1 << "] = " << Q2_ik[ip_0][ip_1] << endl;
      logger << LOG_DEBUG << "temp_mass2 = " << temp_mass2 << endl;
      logger << LOG_DEBUG << "pmass2[" << ip_0 << "] = " << pmass2[ip_0] << endl;
      logger << LOG_DEBUG << "pmass2[" << ip_1 << "] = " << pmass2[ip_1] << endl;
      logger << LOG_DEBUG << "temp_log_s_Q2_ik = " << temp_log_s_Q2_ik << endl;
 
      double temp_V_S_massive_pair_0 = .5 * VA_DeltaIR2 + .5 * log(temp_mass2 / s_ik[ip_0][ip_1]) * VA_DeltaIR1 - 0.25 * pow(log(temp_mass2 / s_ik[ip_0][ip_1]), 2) - pi2_12 - .5 * log(temp_mass2 / s_ik[ip_0][ip_1]) * temp_log_s_Q2_ik - .5 * log(temp_mass2 / Q2_ik[ip_0][ip_1]) * temp_log_s_Q2_ik;
      double temp_V_S_massive_pair_1 = .5 * VA_DeltaIR1 + .5 * log(temp_mass2 / s_ik[ip_0][ip_1]);
      if (switch_polenorm == 1){temp_V_S_massive_pair_0 += -.5 * pi2_6;} // BLHA -> COLI !!!
      if (pair_code_emitter[i_p] != 2){
	V_S_massive_full[ip_0][ip_1][0] = temp_V_S_massive_pair_0;
	V_S_massive_full[ip_0][ip_1][1] = temp_V_S_massive_pair_1;
      }
      if (pair_code_emitter[i_p] != 1){
	V_S_massive_full[ip_1][ip_0][0] = temp_V_S_massive_pair_0;
	V_S_massive_full[ip_1][ip_0][1] = temp_V_S_massive_pair_1;
      }

      logger << LOG_DEBUG << "pair_massive[" << i_p << "] == 1 || pair_massive[" << i_p << "] == 2" << endl;
      logger << LOG_DEBUG << "temp_log_s_Q2_ik = " << temp_log_s_Q2_ik << endl;
      logger << LOG_DEBUG << "pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << endl;
      logger << LOG_DEBUG << "temp_V_S_massive_pair_0 = " << temp_V_S_massive_pair_0 << endl;
      logger << LOG_DEBUG << "temp_V_S_massive_pair_1 = " << temp_V_S_massive_pair_1 << endl;
      logger << LOG_DEBUG << "[" << i_p << "] = " << i_p << endl;
      
      
      if ((pair_massive[i_p] == 1 && pair_code_emitter[i_p] != 2) || (pair_massive[i_p] == 2 && pair_code_emitter[i_p] != 1)){
	// m_i > 0 (quark) and m_k == 0 (gluon/quark)
	//	  logger << LOG_DEBUG_VERBOSE << "m_i > 0 (quark) and m_k == 0 (gluon/quark):   i = " << dx_pa[ioperator_pair[ia][ib]][0] << "   k = " << dx_pa[ioperator_pair[ia][ic]][0] << endl;
	double temp_V_NS_massive_pair = gamma_qew_q * log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) + pi2_6 - gsl_sf_dilog(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) - 2 * log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) - temp_mass2 / s_ik[ip_0][ip_1] * log(temp_mass2 / Q2_ik[ip_0][ip_1]);
	if (pair_massive[i_p] == 1 && pair_code_emitter[i_p] != 2){V_NS_massive_full[ip_0][ip_1] = temp_V_NS_massive_pair;}
	// massive emitter ip_0, massless spectator ip_1
	if (pair_code_emitter[i_p] != 2){	  
	  if (csi->type_parton[0][ip_0] == 0){// ??? gluon ???
	    finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_a + gamma_a(N_f) + K_a(N_f));
	  }
	  else if (csi->type_parton[0][ip_0] != 0){
	    //	    logger << LOG_DEBUG_VERBOSE << "pmass2[" << ip_0 << "] = " << pmass2[ip_0] << endl;
	    //	    logger << LOG_DEBUG_VERBOSE << "pmass2[" << ip_1 << "] = " << pmass2[ip_1] << endl;
	    // !!! temp !!!
	    if (pmass[ip_0] != 0){	
	      Gamma_q[ip_0][0] = (VA_DeltaIR1 + .5 * log(pmass2[ip_0] / s_ik[ip_0][ip_1]) - 2.);
	      finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_q[ip_0][0] + gamma_qew_q + K_qew_q);
	    }
	  }
	  //	  logger << LOG_DEBUG_VERBOSE << "2 pair_massive[" << i_p << "] == " << pair_massive[i_p] << "   pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 2   finite_massive_full[" << ip_0 << "][" << ip_1 << "][0] = " << finite_massive_full[ip_0][ip_1][0] << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << "3 pair_massive[" << i_p << "] == " << pair_massive[i_p] << "   pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 2   finite_massive_full[" << ip_1 << "][" << ip_0 << "][0] = " << finite_massive_full[ip_1][ip_0][0] << endl;
     	}
	if (pair_massive[i_p] == 2 && pair_code_emitter[i_p] != 1){V_NS_massive_full[ip_1][ip_0] = temp_V_NS_massive_pair;}
	  // massive emitter ip_1, massless spectator ip_0
	if (pair_code_emitter[i_p] != 1){
	  if (csi->type_parton[0][ip_1] == 0){
	    finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_a + gamma_a(N_f) + K_a(N_f));
	  }
	  else if (csi->type_parton[0][ip_1] != 0){
	    //	    logger << LOG_DEBUG_VERBOSE << "pmass2[" << ip_0 << "] = " << pmass2[ip_0] << endl;
	    //	    logger << LOG_DEBUG_VERBOSE << "pmass2[" << ip_1 << "] = " << pmass2[ip_1] << endl;
	    // !!! temp !!!
	    if (pmass[ip_1] != 0){	
	      Gamma_q[ip_1][0] = (VA_DeltaIR1 + .5 * log(pmass2[ip_1] / s_ik[ip_0][ip_1]) - 2.);
	      finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_q[ip_1][0] + gamma_qew_q + K_qew_q);
	    }
	  }
	  //	  logger << LOG_DEBUG_VERBOSE << "4 pair_massive[" << i_p << "] == " << pair_massive[i_p] << "   pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 1   finite_massive_full[" << ip_0 << "][" << ip_1 << "][0] = " << finite_massive_full[ip_0][ip_1][0] << endl;
	  //	  logger << LOG_DEBUG_VERBOSE << "5 pair_massive[" << i_p << "] == " << pair_massive[i_p] << "   pair_code_emitter[" << i_p << "] = " << pair_code_emitter[i_p] << " != 1   finite_massive_full[" << ip_1 << "][" << ip_0 << "][0] = " << finite_massive_full[ip_1][ip_0][0] << endl;
	}
      }
      if ((pair_massive[i_p] == 1 && pair_code_emitter[i_p] != 1) || (pair_massive[i_p] == 2 && pair_code_emitter[i_p] != 2)){

	logger << LOG_DEBUG << " ((pair_massive[i_p] == 1 && pair_code_emitter[i_p] != 1) || (pair_massive[i_p] == 2 && pair_code_emitter[i_p] != 2))" << endl;
	
	// m_i == 0 (gluon/quark) and m_k > 0 (quark)
	//	  logger << LOG_DEBUG_VERBOSE << "m_i == 0 (gluon/quark) and m_k > 0 (quark):   i = " << dx_pa[ioperator_pair[ia][ib]][0] << "   k = " << dx_pa[ioperator_pair[ia][ic]][0] << endl;
	double temp_mass = 0.;
	double temp_type = 0.;
	if (pair_massive[i_p] == 1){temp_mass = pmass[ip_0]; temp_type = csi->type_parton[0][ip_1];}
	else if (pair_massive[i_p] == 2){temp_mass = pmass[ip_1]; temp_type = csi->type_parton[0][ip_0];}
	sqrtQ2_ik[ip_0][ip_1] = sqrt(Q2_ik[ip_0][ip_1]);

	double temp_V_NS_massive_pair;
	if (temp_type == 0){
	  // m_i == 0 (gluon) and m_k > 0 (quark)
	  temp_V_NS_massive_pair = gamma_a(N_f) * (log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) - 2 * log((sqrtQ2_ik[ip_0][ip_1] - temp_mass) / sqrtQ2_ik[ip_0][ip_1]) - 2 * temp_mass / (sqrtQ2_ik[ip_0][ip_1] + temp_mass)) + pi2_6 - gsl_sf_dilog(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]); // ...
	  if ((pair_massive[i_p] == 1 && ip_1 > 2) || (pair_massive[i_p] == 2 && ip_0 > 2)){
	    // meaning of '7' has to be re-considered... replace by heaviest quark in real g->QQ~ splittings?
	    for (int i = 1; i < 7; i++){
	      if (M[i] > 0.){
		temp_V_NS_massive_pair += 2. / 3. * N_c * log(M2[i] / Q2_aux); // T_R / C_A -> N_c
		if (s_ik[ip_0][ip_1] > 4 * M[i] * (temp_mass + M[i])){
		  double rho1 = sqrt(1. - 4 * M2[i] / pow(sqrtQ2_ik[ip_0][ip_1] - temp_mass, 2));
		  temp_V_NS_massive_pair += 4. / 3. * N_c * (log((sqrtQ2_ik[ip_0][ip_1] - temp_mass) / sqrtQ2_ik[ip_0][ip_1]) + temp_mass * pow(rho1, 3) / (sqrtQ2_ik[ip_0][ip_1] + temp_mass) + log((1. + rho1) / 2.) - rho1 / 3. * (3. + pow(rho1, 2)) - .5 * log(M2[i] / Q2_ik[ip_0][ip_1]));
		}
	      }
	    }
	  }
	}
	else if (temp_type != 0){
	  // m_i == 0 (quark) and m_k > 0 (quark)
	  temp_V_NS_massive_pair = gamma_qew_q * (log(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]) - 2 * log((sqrtQ2_ik[ip_0][ip_1] - temp_mass) / sqrtQ2_ik[ip_0][ip_1]) - 2 * temp_mass / (sqrtQ2_ik[ip_0][ip_1] + temp_mass)) + pi2_6 - gsl_sf_dilog(s_ik[ip_0][ip_1] / Q2_ik[ip_0][ip_1]);
	}
	if (pair_massive[i_p] == 1 && pair_code_emitter[i_p] != 1){V_NS_massive_full[ip_1][ip_0] = temp_V_NS_massive_pair;} // massive emitter ip_1, massless spectator ip_0
	if (pair_massive[i_p] == 2 && pair_code_emitter[i_p] != 2){V_NS_massive_full[ip_0][ip_1] = temp_V_NS_massive_pair;} // massive emitter ip_0, massless spectator ip_1
      }
    }
    else if (pair_massive[i_p] == 0){
      // m_i == 0 (gluon/quark) and m_k == 0 (gluon/quark)
      //	  logger << LOG_DEBUG_VERBOSE << "m_i == 0 (gluon/quark) and m_k == 0 (gluon/quark):   i = " << dx_pa[ioperator_pair[ia][ib]][0] << "   k = " << dx_pa[ioperator_pair[ia][ic]][0] << endl;
      if ((csi->type_parton[0][ip_0] == 0 && ip_0 > 2 && pair_code_emitter[i_p] != 2) || (csi->type_parton[0][ip_1] == 0 && ip_1 > 2 && pair_code_emitter[i_p] != 1)){
	double temp_V_NS_massive = 0.; // ...
	double m_k = 0.; // ???
	for (int i = 1; i < 7; i++){
	  if (M[i] > 0.){
	    temp_V_NS_massive += 2. / 3. * N_c * log(M2[i] / Q2_aux);
	    //		logger << LOG_DEBUG_VERBOSE << "V_NS_massive(0., 0.)[" << i << ", N_F] = " << V_NS_massive << endl;
	    if (s_ik[ip_0][ip_1] > 4 * M[i] * (m_k + M[i])){
	      double rho1 = sqrt(1. - 4 * M2[i] / pow(sqrtQ2_ik[ip_0][ip_1] - m_k, 2));
	      temp_V_NS_massive += 4. / 3. * N_c * (log((1. + rho1) / 2.) - rho1 / 3. * (3. + pow(rho1, 2)) - .5 * log(M2[i] / s_ik[ip_0][ip_1]));
	      //		  logger << LOG_DEBUG_VERBOSE << "V_NS_massive(0., 0.)[" << i << ", N_F^ik] = " << V_NS_massive << endl;
	    }
	  }
	}
	if (csi->type_parton[0][ip_0] == 0 && ip_0 > 2 && pair_code_emitter[i_p] != 2){V_NS_massive_full[ip_0][ip_1] = temp_V_NS_massive;}
	if (csi->type_parton[0][ip_1] == 0 && ip_1 > 2 && pair_code_emitter[i_p] != 1){V_NS_massive_full[ip_1][ip_0] = temp_V_NS_massive;}
      }
      if (csi->type_parton[0][ip_0] == 0 && pair_code_emitter[i_p] != 2){
	if (csi->type_parton[0][ip_0] == 0){
	  finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_a + gamma_a(N_f) + K_a(N_f));
	}
	else if (csi->type_parton[0][ip_0] != 0){
	  Gamma_q[ip_0][0] = (VA_DeltaIR1 + .5 * log(pmass2[ip_0] / s_ik[ip_0][ip_1]) - 2.);
	  finite_massive_full[ip_0][ip_1][0] = - pi2_3 + (Gamma_q[ip_0][0] + gamma_qew_q + K_qew_q);
	}
      }
      if (csi->type_parton[0][ip_1] == 0 && pair_code_emitter[i_p] != 1){
	if (csi->type_parton[0][ip_1] == 0){
	  finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_a + gamma_a(N_f) + K_a(N_f));
	}
	else if (csi->type_parton[0][ip_1] != 0){
	  Gamma_q[ip_1][0] = (VA_DeltaIR1 + .5 * log(pmass2[ip_1] / s_ik[ip_0][ip_1]) - 2.);
	  finite_massive_full[ip_1][ip_0][0] = - pi2_3 + (Gamma_q[ip_1][0] + gamma_qew_q + K_qew_q);
	}
      }
      //	logger << LOG_DEBUG_VERBOSE << "finite_massive_full[" << ip_0 << "][" << ip_1 << "][0] = " << finite_massive_full[ip_0][ip_1][0] << endl;
      // 	logger << LOG_DEBUG_VERBOSE << "finite_massive_full[" << ip_1 << "][" << ip_0 << "][0] = " << finite_massive_full[ip_1][ip_0][0] << endl;

    }
  }
  //  logger << LOG_DEBUG_VERBOSE << "xxx done!" << endl;

  for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
    for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
      if ((*VA_ioperator)[i_a][0].type_correction() == 2){
	int ip_0 = (*VA_ioperator)[i_a][j_a].pair()[0];
	int ip_1 = (*VA_ioperator)[i_a][j_a].pair()[1];
	int no_em = (*VA_ioperator)[i_a][j_a].no_emitter();
	int no_sp = (*VA_ioperator)[i_a][j_a].no_spectator();
	/*      
		prefactor_log[0] = V_S_massive_full[no_em][no_sp][0] + V_NS_massive_full[no_em][no_sp] + finite_massive_full[no_em][no_sp][0];
		prefactor_log[1] = V_S_massive_full[no_em][no_sp][1] + finite_massive_full[no_em][no_sp][1];
		prefactor_log[2] = V_S_massive_full[no_em][no_sp][2];
		//      prefactor_I_ME2_cs[i_a][j_a] = prefactor_log[0] + prefactor_log[1] * ioperator_pair_log[ia] + prefactor_log[2] * ioperator_pair_log2_2[ia];
      */
	VA_I_ME2_cf[i_a][j_a] = -alpha_e_2pi * ((V_S_massive_full[no_em][no_sp][0] + V_NS_massive_full[no_em][no_sp] + finite_massive_full[no_em][no_sp][0]) 
					     + (V_S_massive_full[no_em][no_sp][1] + finite_massive_full[no_em][no_sp][1]) * ioperator_pair_log[ip_0][ip_1] 
					     + (V_S_massive_full[no_em][no_sp][2]) * ioperator_pair_log2_2[ip_0][ip_1]) * VA_ME2_cf[i_a][j_a];
	
	  logger << LOG_DEBUG << (*VA_ioperator)[i_a][j_a].name() << endl;
	  logger << LOG_DEBUG << "VA_ME2_cf[" << i_a << "][" << j_a << "] = " << VA_ME2_cf[i_a][j_a] << endl;
	  logger << LOG_DEBUG << "    prefactor_log[" << no_em << "][" << no_sp << "][0]    = " << (V_S_massive_full[no_em][no_sp][0] + V_NS_massive_full[no_em][no_sp] + finite_massive_full[no_em][no_sp][0]) << endl;
	  logger << LOG_DEBUG << "    prefactor_log[" << no_em << "][" << no_sp << "][1]   = " << (V_S_massive_full[no_em][no_sp][1] + finite_massive_full[no_em][no_sp][1]) << endl;
	  logger << LOG_DEBUG << "    prefactor_log[" << no_em << "][" << no_sp << "][2]   = " << (V_S_massive_full[no_em][no_sp][2]) << endl;
	  logger << LOG_DEBUG << "   V_S_massive_full[" << no_em << "][" << no_sp << "][0] = " << V_S_massive_full[no_em][no_sp][0] << endl;
	  logger << LOG_DEBUG << "   V_S_massive_full[" << no_em << "][" << no_sp << "][1] = " << V_S_massive_full[no_em][no_sp][1] << endl;
	  logger << LOG_DEBUG << "   V_S_massive_full[" << no_em << "][" << no_sp << "][2] = " << V_S_massive_full[no_em][no_sp][2] << endl;
	  logger << LOG_DEBUG << "  V_NS_massive_full[" << no_em << "][" << no_sp << "]    = " << V_NS_massive_full[no_em][no_sp] << endl;
	  logger << LOG_DEBUG << "finite_massive_full[" << no_em << "][" << no_sp << "][0] = " << finite_massive_full[no_em][no_sp][0] << endl;
	  logger << LOG_DEBUG << "finite_massive_full[" << no_em << "][" << no_sp << "][1] = " << finite_massive_full[no_em][no_sp][1] << endl;
	  logger << LOG_DEBUG << "   ioperator_pair_log[" << ip_0 << "][" << ip_1 << "]    = " << ioperator_pair_log[ip_0][ip_1] << endl;
	  logger << LOG_DEBUG << "ioperator_pair_log2_2[" << ip_0 << "][" << ip_1 << "]    = " << ioperator_pair_log2_2[ip_0][ip_1] << endl;
	  logger << LOG_DEBUG << "VA_I_ME2_cf[" << i_a << "][" << j_a << "] = " << VA_I_ME2_cf[i_a][j_a] << endl;
	
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
