#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

void observable_set::determine_collinear_QCD(phasespace_set & psi){
  Logger logger("observable_set::determine_collinear_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  //  static map<int, double> charge_particle;
  if (initialization == 1){
    //    fill_charge_particle(charge_particle);
    initialization = 0;
  }
  logger << LOG_DEBUG << "new collinear-dipole determination" << endl << endl;
  for (int i1 = 0; i1 < combination_pdf.size(); i1++){
    stringstream temp_ss;
    temp_ss << " combination_pdf[" << setw(2) << i1 << "] = "; 
    for (int i = 0; i < 3; i++){temp_ss << setw(4) << combination_pdf[i1][i] << "   ";}
    logger << LOG_DEBUG << temp_ss.str() << endl;
  }
  logger.newLine(LOG_DEBUG);

  vector<string> pa_name(csi->type_parton[0].size(), "");
  if (csi->type_parton[0][1] >= -10 && csi->type_parton[0][1] <= 10){pa_name[1] = "a";}
  if (csi->type_parton[0][2] >= -10 && csi->type_parton[0][2] <= 10){pa_name[2] = "b";}
  //  int count = 0;
  //  vector<string> alphabet(csi->type_parton[0].size() - 3, "");
  //  for (int i_p = 0; i_p < alphabet.size(); i_p++){alphabet[i_p] = char(105 + i_p);}
  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (csi->type_parton[0][i_p] >= -10 && csi->type_parton[0][i_p] <= 10){pa_name[i_p] = char(105 + i_p - 3);}} 
  //  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (csi->type_parton[0][i_p] >= -10 && csi->type_parton[0][i_p] =< 10){pa_name[i_p] = alphabet[i_p - 3];}} 
  //  new labels: omit colourless particles in counting
  //  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (csi->type_parton[0][i_p] >= -10 && csi->type_parton[0][i_p] =< 10){pa_name[i_p] = alphabet[count++];}}

  stringstream temp_tp;
  stringstream temp_name;
  for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){
    temp_tp << setw(3) << csi->type_parton[0][i_p] << " ";
    temp_name << setw(3) << pa_name[i_p] << " ";
  }
  logger.newLine(LOG_INFO);
  logger << LOG_INFO << "Collinear QCD emission:" << endl;
  logger << LOG_INFO << "csi->type_parton[0] = " << temp_tp.str() << endl;
  logger << LOG_INFO << "pa_name             = " << temp_name.str() << endl;



  int temp_type_correction = 1;

  /*
  map <int,string> pname;
  fill_pname(pname);
  */

  for (int temp_no_emitter = 1; temp_no_emitter < 3; temp_no_emitter++){
    if (pa_name[temp_no_emitter] == ""){continue;}
    string temp_name;
    vector<string> temp_all_name(combination_pdf.size());
    int temp_type;
    vector<int> temp_in_collinear(3, 0);
    temp_in_collinear[temp_no_emitter] = 1;
    vector<vector<int> > temp_pdf_new;
    double temp_charge_factor = 1.;
    for (int i_t = 0; i_t < 2; i_t++){
      temp_pdf_new = combination_pdf;
      if (combination_pdf[0][temp_no_emitter] == 0){
	if (i_t == 0){
	  temp_type = 0; // hard process with g from g -> g g splitting
	  temp_in_collinear[0] = 1;
	  temp_name = "C_" + pa_name[temp_no_emitter] + "^{gg}";
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_all_name[i_x] = "C_" + pa_name[temp_no_emitter] + "^{gg}";}
	}
	else {
	  temp_type = 2; // hard process with a from q -> g q splitting
	  temp_in_collinear[0] = 0;
	  temp_name = "C_" + pa_name[temp_no_emitter] + "^{qg}";
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_pdf_new[i_x][temp_no_emitter] = 10;}
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_all_name[i_x] = "C_" + pa_name[temp_no_emitter] + "^{qg}";}
	}
      }
      else if (combination_pdf[0][temp_no_emitter] != 0){
	if (i_t == 0){
	  temp_type = 3; // hard process with q from g -> q qx splitting
	  temp_in_collinear[0] = 0;
	  temp_name = "C_" + pa_name[temp_no_emitter] + "^{g" + csi->name_particle[combination_pdf[0][temp_no_emitter]] + "}";
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_pdf_new[i_x][temp_no_emitter] = 0;}
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_all_name[i_x] = "C_" + pa_name[temp_no_emitter] + "^{g" + csi->name_particle[combination_pdf[i_x][temp_no_emitter]] + "}";}
	}
	else {
	  temp_type = 1; // hard process with q from q -> q g splitting
	  temp_in_collinear[0] = 1;
	  temp_name = "C_" + pa_name[temp_no_emitter] + "^{" + csi->name_particle[combination_pdf[0][temp_no_emitter]] + csi->name_particle[combination_pdf[0][temp_no_emitter]] + "}";
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_all_name[i_x] = "C_" + pa_name[temp_no_emitter] + "^{" + csi->name_particle[combination_pdf[i_x][temp_no_emitter]] + csi->name_particle[combination_pdf[i_x][temp_no_emitter]] + "}";}
	}
      }
      /*
      if (LHAPDF::hasPhoton() == 0 && (temp_type == 0 || temp_type == 3 || csi->type_parton[0][temp_no_emitter % 2 + 1] == 22)){
	logger << LOG_DEBUG << "No photon pdf available in pdf set " << LHAPDFname << endl;
	continue;
      }
      */
      int flag = (*CA_collinear).size();
      for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){if (temp_no_emitter == (*CA_collinear)[i_c][0].no_emitter() && temp_type == (*CA_collinear)[i_c][0].type()){flag = i_c; break;}}
      if (flag == (*CA_collinear).size()){(*CA_collinear).push_back(vector<collinear_set> ());}
      for (int temp_no_spectator = 0; temp_no_spectator < csi->type_parton[0].size(); temp_no_spectator++){
	if ((temp_no_spectator > 0 && pa_name[temp_no_spectator] == "") || temp_no_spectator == temp_no_emitter){continue;}
	vector<int> temp_pair(2);
	if (temp_no_emitter < temp_no_spectator){
	  temp_pair[0] = temp_no_emitter;
	  temp_pair[1] = temp_no_spectator;
	}
	else {
	  temp_pair[0] = temp_no_spectator;
	  temp_pair[1] = temp_no_emitter;
	}

	// already existing: mass_parton or so... ???
	int temp_massive;
	if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 0;}
	else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 1;}
	else if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 2;}
	else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 3;}
	else {cout << "Should not happen!" << endl;}


	vector<string> temp_all_name_spectator = temp_all_name;
	for (int i_x = 0; i_x < temp_all_name_spectator.size(); i_x++){temp_all_name_spectator[i_x] = temp_all_name_spectator[i_x] + "(" + pa_name[temp_no_spectator] + ")";}
	string temp_name_spectator = temp_name + "(" + pa_name[temp_no_spectator] + ")";
	(*CA_collinear)[flag].push_back(collinear_set(temp_name_spectator, temp_all_name_spectator, temp_type, temp_in_collinear, psi_no_prc[0], csi->type_parton[0], temp_pdf_new, temp_charge_factor, temp_no_emitter, temp_no_spectator, temp_pair, temp_type_correction, temp_massive));
      }
    }
  }

  logger.newLine(LOG_INFO);

  logger << LOG_INFO << "Before selection of contributing collinear 'dipoles':" << endl;
  logger.newLine(LOG_INFO);
  output_collinear();

  /*
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
      stringstream temp;
      if      ((*CA_collinear)[i_c][j_c].type_correction() == 1){temp << "QCD ";}
      else if ((*CA_collinear)[i_c][j_c].type_correction() == 2){temp << "QEW ";}
      temp << "collinear dipole " << i_c << ", " << j_c << ":   " << setw(15) << left << (*CA_collinear)[i_c][j_c].name() << ":   ";
      temp << "massive = " << (*CA_collinear)[i_c][j_c].massive() << "   ";
      temp << "psi_no_prc[0] = " << (*CA_collinear)[i_c][j_c].no_prc() << "   ";
      temp << "csi->type_parton[0] = ";
      for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){temp << (*CA_collinear)[i_c][j_c].type_parton()[i_p] << "   ";}
      temp << "charge_factor = " << setw(23) << setprecision(15) << (*CA_collinear)[i_c][j_c].charge_factor() << "   ";
      logger << LOG_DEBUG << temp.str() << endl;
      temp.str("");
      temp << setw(48) << "";
      temp << "type = " << (*CA_collinear)[i_c][j_c].type() << "   ";
      temp << "in_collinear = ";
      for (int i_em = 0; i_em < 3; i_em++){temp << (*CA_collinear)[i_c][j_c].in_collinear()[i_em] << "   ";}
      temp << "no_emitter = " << (*CA_collinear)[i_c][j_c].no_emitter() << "   ";
      temp << "no_spectator = " << (*CA_collinear)[i_c][j_c].no_spectator() << "   ";
      temp << "pair = ";
      for (int i_p = 0; i_p < 2; i_p++){temp << (*CA_collinear)[i_c][j_c].pair()[i_p] << "   ";}
      logger << LOG_DEBUG << temp.str() << endl;
      for (int i_x = 0; i_x < (*CA_collinear)[i_c][j_c].all_name().size(); i_x++){
	temp.str("");
	temp << setw(48) << "";
	temp << "pdf contributions: " << setw(15) << left << (*CA_collinear)[i_c][j_c].all_name()[i_x] << ":   ";
	for (int i_yxy = 0; i_yxy < 3; i_yxy++){temp << setw(2) << right << (*CA_collinear)[i_c][j_c].all_pdf()[i_x][i_yxy] << "  ";}
	logger << LOG_DEBUG << temp.str() << endl;
      }
    }
  }
  logger.newLine(LOG_DEBUG);
*/

  CA_dipole_splitting.resize(4, vector<int> (3, 0));

  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_em = 0; i_em < 3; i_em++){
      if ((*CA_collinear)[i_c][0].in_collinear()[i_em] == 1){
	CA_dipole_splitting[(*CA_collinear)[i_c][0].type()][i_em] = 1;
      }
    }
  }

  logger << LOG_DEBUG << "CA_dipole_splitting:" << endl;
  for (int i_em = 0; i_em < 3; i_em++){
    stringstream temp_ss;
    temp_ss << "i_em = " << i_em << ":   ";
    for (int i_t = 0; i_t < 4; i_t++){temp_ss << CA_dipole_splitting[i_t][i_em] << "   ";}
    logger << LOG_DEBUG << temp_ss.str() << endl;
  }
  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "new collinear dipoles determined " << endl << endl;


  CA_combination_pdf.resize((*CA_collinear).size(), vector<vector<vector<int> > > (combination_pdf.size()));
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      // CA_combination_pdf[i_c] contains a vector, which contains the usual osi_combination_pdf(3) with 0 -> direction (1, -1), 1 -> parton with x1 (in hadron 1/2 for +1/-1), 2 -> parton with x2 (in hadron 2/1 for +1/-1)
      if (((*CA_collinear)[i_c][0].all_pdf()[i_i][1] == 10) && ((*CA_collinear)[i_c][0].all_pdf()[i_i][2] == 10)){
	for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	  if (i_q == 0){continue;}
	  for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	    if (j_q == 0){continue;}
	    vector<int> new_temp_combination_pdf = (*CA_collinear)[i_c][0].all_pdf()[i_i];
	    new_temp_combination_pdf[1] = i_q;
	    new_temp_combination_pdf[2] = j_q;
	    CA_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	  }
	}
      }
      else if ((*CA_collinear)[i_c][0].all_pdf()[i_i][1] == 10){
	for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	  if (i_q == 0){continue;}
	  vector<int> new_temp_combination_pdf = (*CA_collinear)[i_c][0].all_pdf()[i_i];
	  new_temp_combination_pdf[1] = i_q;
	  CA_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	}
      }
      else if ((*CA_collinear)[i_c][0].all_pdf()[i_i][2] == 10){
	for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	  if (j_q == 0){continue;}
	  vector<int> new_temp_combination_pdf = (*CA_collinear)[i_c][0].all_pdf()[i_i];
	  new_temp_combination_pdf[2] = j_q;
	  CA_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	}
      }
      else {
	CA_combination_pdf[i_c][i_i].push_back((*CA_collinear)[i_c][0].all_pdf()[i_i]);
      }
    }
  }
  logger << LOG_DEBUG << "CA_combination_pdf determined " << endl << endl;

  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      logger << LOG_DEBUG << "pdf contributions: " << setw(15) << left << (*CA_collinear)[i_c][0].all_name()[i_i] << ":   " << endl;
      for (int i_x = 0; i_x < CA_combination_pdf[i_c][i_i].size(); i_x++){
	stringstream temp;
	temp.str("");
	temp << setw(48) << "";
	temp << "CA_combination_pdf[" << setw(2) << i_c << "][" << setw(2) << i_i << "][" << setw(2) << i_x << "] = "; 
	for (int i_y = 0; i_y < CA_combination_pdf[i_c][i_i][i_x].size(); i_y++){temp << setw(2) << right << CA_combination_pdf[i_c][i_i][i_x][i_y] << "  ";}
	logger << LOG_DEBUG << temp.str() << endl;
      }
    }
  }

  CA_Q2f.resize((*CA_collinear).size());
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    CA_Q2f[i_c].resize(CA_combination_pdf[i_c].size());
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      for (int i_x = 0; i_x < CA_combination_pdf[i_c][i_i].size(); i_x++){
	CA_Q2f[i_c][i_i].push_back(1.);
      }
    }
  }
  logger << LOG_DEBUG << "new collinear-dipole Q2f determined " << endl << endl;
  logger << LOG_DEBUG << "CA_Q2f.size() = " << CA_Q2f.size() << endl;
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    logger << LOG_DEBUG << "CA_Q2f[" << i_c << "].size() = " << CA_Q2f[i_c].size() << endl;
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      logger << LOG_DEBUG << "CA_Q2f[" << i_c << "][" << i_i << "].size() = " << CA_Q2f[i_c][i_i].size() << endl;
    }
  }


  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_collinear_QCD(){
  Logger logger("observable_set::calculate_collinear_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (massive_QCD){calculate_collinear_QCD_CDST();}
  else {calculate_collinear_QCD_CS();}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_collinear_QCD_CS(){
  Logger logger("observable_set::calculate_collinear_QCD_CS");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  // calculate all needed momentum-independent (splitting) functions
  static double Kbar[3][4] = {{0.}};
  static double Kt[3][4] = {{0.}};
  static double P[3][4] = {{0.}};
  static double Kbar_plus[3][4] = {{0.}};
  static double Kt_plus[3][4] = {{0.}};
  static double P_plus[3][4] = {{0.}};
  static double intKbar_plus[3][4] = {{0.}};
  static double intKt_plus[3][4] = {{0.}};
  static double intP_plus[3][4] = {{0.}};
  static double Kbar_delta[4] = {0.};
  static double Kt_delta[4] = {0.};
  static double P_delta[4] = {0.};

  static double alpha_S_2pi = alpha_S * inv2pi;
  static vector<vector<int> > pair;
  static vector<double> iT2_ap(3);
  static vector<vector<double> > gamma_i_T2_i((*CA_collinear).size());
  static vector<vector<double> > ln_papi(3, vector<double> (csi->type_parton[0].size()));
  static vector<vector<vector<vector<double> > > > CA_value_ln_muF_papi(CA_value_log_mu2_fact.size());
  static vector<vector<vector<vector<vector<double> > > > > value_dataP(CA_value_log_mu2_fact.size());

  static vector<vector<int> > type((*CA_collinear).size());
  static vector<vector<int> > no_emitter((*CA_collinear).size());
  static vector<vector<int> > no_spectator((*CA_collinear).size());
  static vector<vector<int> > collinear_singularity((*CA_collinear).size());
  static vector<vector<vector<int> > > ppair((*CA_collinear).size());

  if (initialization == 1){
    if (CA_dipole_splitting[0][1] == 1 || CA_dipole_splitting[0][2] == 1){    // g -> g (+g) splitting (0)
      Kbar_delta[0] = Kbar_gg_delta(N_f);
      Kt_delta[0] = Kt_gg_delta();
      P_delta[0] = P_gg_delta(N_f);
    }
    if (CA_dipole_splitting[1][1] == 1 || CA_dipole_splitting[1][2] == 1){    // q -> q (+g) splitting
      Kbar_delta[1] = Kbar_qq_delta();
      Kt_delta[1] = Kt_qq_delta();
      P_delta[1] = Pxx_qq_delta();
    }

    for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	if ((*CA_collinear)[i_c][j_c].no_spectator() == 0){continue;}
	int flag = -1;
	for (int i_p = 0; i_p < pair.size(); i_p++){
	  if ((*CA_collinear)[i_c][j_c].pair() == pair[i_p]){flag = i_p; break;}
	}
	if (flag == -1){pair.push_back((*CA_collinear)[i_c][j_c].pair());}
      }
    }

    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      CA_value_ln_muF_papi[sd].resize(CA_value_log_mu2_fact[sd].size(), vector<vector<double> > (3, vector<double> (csi->type_parton[0].size())));
    }
    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      value_dataP[sd].resize(CA_value_log_mu2_fact[sd].size(), vector<vector<vector<double> > > (3, vector<vector<double> > ((*CA_collinear).size())));
      for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	for (int i_x = 0; i_x < 3; i_x++){
	  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
	    value_dataP[sd][ss][i_x][i_c].resize((*CA_collinear)[i_c].size(), 0.);
	  }
	}
      }
    }

    for (int i_em = 1; i_em < 3; i_em++){
      if (csi->type_parton[0][i_em] == 0){iT2_ap[i_em] = 1. / C_A;}
      else {iT2_ap[i_em] = 1. / C_F;}
    }

    for (int i_c = 0; i_c < gamma_i_T2_i.size(); i_c++){ // i_c = -1 is not needed !!!
      gamma_i_T2_i[i_c].resize((*CA_collinear)[i_c].size());
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	if ((*CA_collinear)[i_c][j_c].type_parton()[(*CA_collinear)[i_c][j_c].no_spectator()] == 0){gamma_i_T2_i[i_c][j_c] = gamma_g(N_f) / C_A;}
	else {gamma_i_T2_i[i_c][j_c] = gamma_q / C_F;}
      }
    }

    for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
      type[i_c].resize((*CA_collinear)[i_c].size());
      no_emitter[i_c].resize((*CA_collinear)[i_c].size());
      no_spectator[i_c].resize((*CA_collinear)[i_c].size());
      collinear_singularity[i_c].resize((*CA_collinear)[i_c].size());
      ppair[i_c].resize((*CA_collinear)[i_c].size());
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	type[i_c][j_c] = (*CA_collinear)[i_c][j_c].type();
	no_emitter[i_c][j_c] = (*CA_collinear)[i_c][j_c].no_emitter();
	no_spectator[i_c][j_c] = (*CA_collinear)[i_c][j_c].no_spectator();
	collinear_singularity[i_c][j_c] = (*CA_collinear)[i_c][j_c].in_collinear()[0];
	ppair[i_c][j_c] = (*CA_collinear)[i_c][j_c].pair();
      }
    }

    for (int i_p = 0; i_p < pair.size(); i_p++){logger << LOG_DEBUG << "pair[" << i_p << "] = " << "(" << pair[i_p][0] << ", " << pair[i_p][1] << ")" << endl;}
    initialization = 0;
    logger << LOG_DEBUG_VERBOSE << "initialization finished!" << endl;
  }


  for (int i_em = 1; i_em < 3; i_em++){
    logger << LOG_DEBUG_VERBOSE << "z_coll[" << i_em << "] = " << z_coll[i_em] << "   CA_dipole_splitting[0][" << i_em << "] = " << CA_dipole_splitting[0][i_em] << endl;
    logger << LOG_DEBUG_VERBOSE << "z_coll[" << i_em << "] = " << z_coll[i_em] << "   CA_dipole_splitting[1][" << i_em << "] = " << CA_dipole_splitting[1][i_em] << endl;
    logger << LOG_DEBUG_VERBOSE << "z_coll[" << i_em << "] = " << z_coll[i_em] << "   CA_dipole_splitting[2][" << i_em << "] = " << CA_dipole_splitting[2][i_em] << endl;
    logger << LOG_DEBUG_VERBOSE << "z_coll[" << i_em << "] = " << z_coll[i_em] << "   CA_dipole_splitting[3][" << i_em << "] = " << CA_dipole_splitting[3][i_em] << endl;
    if (CA_dipole_splitting[0][i_em] == 1){      // g -> g (+g) splitting
      Kbar[i_em][0] = Kbar_gg(z_coll[i_em]);
      Kt[i_em][0] = Kt_gg(z_coll[i_em]);
      P[i_em][0] = P_gg(z_coll[i_em]);
      if (CA_dipole_splitting[0][0] == 1){
	Kbar_plus[i_em][0] = Kbar_gg_plus(z_coll[i_em]);
	Kt_plus[i_em][0] = Kt_gg_plus(z_coll[i_em]);
	P_plus[i_em][0] = P_gg_plus(z_coll[i_em]);
	intKbar_plus[i_em][0] = intKbar_gg_plus(x_pdf[i_em]);
	intKt_plus[i_em][0] = intKt_gg_plus(x_pdf[i_em]);
	intP_plus[i_em][0] = intP_gg_plus(x_pdf[i_em]);
      }
      else {logger << LOG_ERROR << "May not happen !!! g -> g splitting without irregular terms !!!" << endl;}
    }
    if (CA_dipole_splitting[1][i_em] == 1){      // q -> q (+g) splitting
      Kbar[i_em][1] = Kbar_qq(z_coll[i_em]);
      Kt[i_em][1] = Kt_qq(z_coll[i_em]);
      P[i_em][1] = Pxx_qq(z_coll[i_em]);
      if (CA_dipole_splitting[1][0] == 1){
	Kbar_plus[i_em][1] = Kbar_qq_plus(z_coll[i_em]);
	Kt_plus[i_em][1] = Kt_qq_plus(z_coll[i_em]);
	P_plus[i_em][1] = Pxx_qq_plus(z_coll[i_em]);
	intKbar_plus[i_em][1] = intKbar_qq_plus(x_pdf[i_em]);
	intKt_plus[i_em][1] = intKt_qq_plus(x_pdf[i_em]);
	intP_plus[i_em][1] = intPxx_qq_plus(x_pdf[i_em]);
      }
      else {logger << LOG_ERROR << "May not happen !!! q -> q splitting without irregular terms !!!" << endl;}
    }
    if (CA_dipole_splitting[2][i_em] == 1){      // q -> g (+q) splitting
      Kbar[i_em][2] = Kbar_qg(z_coll[i_em]);
      Kt[i_em][2] = Kt_qg(z_coll[i_em]);
      P[i_em][2] = P_qg(z_coll[i_em]);
    }
    if (CA_dipole_splitting[3][i_em] == 1){      // g -> q (+q~) splitting
      Kbar[i_em][3] = Kbar_gq(z_coll[i_em]);
      Kt[i_em][3] = Kt_gq(z_coll[i_em]);
      P[i_em][3] = P_gq(z_coll[i_em]);
    }
  }

  
    if (switch_TSV){
      for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	for (int v_xf = 0; v_xf < n_scale_dyn_fact[v_sf]; v_xf++){
	  logger << LOG_DEBUG << "value_central_logscale2_fact[" << v_sf << "] = " << value_central_logscale2_fact[v_sf] << endl;
	  logger << LOG_DEBUG << "value_relative_logscale2_fact[" << v_sf << "][" << v_xf << "] = " << value_relative_logscale2_fact[v_sf][v_xf] << endl;
	  logger << LOG_DEBUG << "relative_central_scale_fact_TSV[" << v_sf << "] = " << relative_central_scale_fact_TSV[v_sf] << endl;
	  logger << LOG_DEBUG << "relative_scale_fact_TSV[" << v_sf << "][" << v_xf << "] = " << relative_scale_fact_TSV[v_sf][v_xf] << endl;
	}
      }
    }
  
   for (int i_p = 0; i_p < pair.size(); i_p++){
    ln_papi[pair[i_p][0]][pair[i_p][1]] = log(2 * p_parton[0][pair[i_p][0]] * p_parton[0][pair[i_p][1]]);
    if (switch_TSV){
      for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	for (int v_xf = 0; v_xf < n_scale_dyn_fact[v_sf]; v_xf++){
	  value_logscale2_fact_papi[v_sf][v_xf][pair[i_p][0]][pair[i_p][1]] = value_central_logscale2_fact[v_sf] + value_relative_logscale2_fact[v_sf][v_xf] - ln_papi[pair[i_p][0]][pair[i_p][1]];
	  logger << LOG_DEBUG << "value_logscale2_fact_papi[" << v_sf << "][" << v_xf << "][" << pair[i_p][0] << "][" << pair[i_p][1] << "] = " << value_logscale2_fact_papi[v_sf][v_xf][pair[i_p][0]][pair[i_p][1]] << endl;
	}
      }
    }
    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	CA_value_ln_muF_papi[sd][ss][pair[i_p][0]][pair[i_p][1]] = CA_value_log_mu2_fact[sd][ss] - ln_papi[pair[i_p][0]][pair[i_p][1]];
      }
    }
  }

   
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
      // K terms 
      if (no_spectator[i_c][j_c] == 0){
	data_K[0][i_c][j_c] = (Kbar[no_emitter[i_c][j_c]][type[i_c][j_c]] + Kbar_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
	data_K[1][i_c][j_c] = (-Kbar_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];
	data_K[2][i_c][j_c] = (Kbar_delta[type[i_c][j_c]] - intKbar_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
      }
      else if (no_spectator[i_c][j_c] < 3){
	data_K[0][i_c][j_c] = -iT2_ap[no_emitter[i_c][j_c]] * (Kt[no_emitter[i_c][j_c]][type[i_c][j_c]] + Kt_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
	data_K[1][i_c][j_c] = -iT2_ap[no_emitter[i_c][j_c]] * (-Kt_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];
	data_K[2][i_c][j_c] = -iT2_ap[no_emitter[i_c][j_c]] * (Kt_delta[type[i_c][j_c]] - intKt_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
      }
      else if (no_spectator[i_c][j_c] > 2){
	if (collinear_singularity[i_c][j_c] == 1){
	  data_K[0][i_c][j_c] = gamma_i_T2_i[i_c][j_c] * (1. / (1. - z_coll[no_emitter[i_c][j_c]])) * CA_ME2_cf[i_c][j_c];
	  data_K[1][i_c][j_c] = -gamma_i_T2_i[i_c][j_c] * (1. / (1. - z_coll[no_emitter[i_c][j_c]])) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];
	  data_K[2][i_c][j_c] = gamma_i_T2_i[i_c][j_c] * (1. + log(1. - x_pdf[no_emitter[i_c][j_c]])) * CA_ME2_cf[i_c][j_c];
	}
      }
      if (no_spectator[i_c][j_c] != 0){
	// P terms 
	if (switch_TSV){
	  for (int v_sf = 0; v_sf < value_logscale2_fact_papi.size(); v_sf++){
	    for (int v_xf = 0; v_xf < value_logscale2_fact_papi[v_sf].size(); v_xf++){
	      value_data_P[v_sf][v_xf][0][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (P[no_emitter[i_c][j_c]][type[i_c][j_c]] + P_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	      value_data_P[v_sf][v_xf][1][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (-P_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	      value_data_P[v_sf][v_xf][2][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (P_delta[type[i_c][j_c]] - intP_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	    }
	  }
	}
	for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
	  for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	    value_dataP[sd][ss][0][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (P[no_emitter[i_c][j_c]][type[i_c][j_c]] + P_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	    value_dataP[sd][ss][1][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (-P_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	    value_dataP[sd][ss][2][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (P_delta[type[i_c][j_c]] - intP_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	  }
	}
      }
    }
    for (int i_x = 0; i_x < 3; i_x++){
      double temp_sumK = accumulate(data_K[i_x][i_c].begin(), data_K[i_x][i_c].end(), 0.);
      logger << LOG_DEBUG_VERBOSE << "i_x = " << i_x << "   temp_sumK = " << temp_sumK << endl;
      if (switch_TSV){
	for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	  for (int v_xf = 0; v_xf < value_logscale2_fact_papi[v_sf].size(); v_xf++){
	    if (switch_KP == 0){
	      // K + P terms
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_S_2pi * (temp_sumK + accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_S_2pi * (temp_sumK + accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	    }
	    else if (switch_KP == 1){
	      // P terms 
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_S_2pi * (accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_S_2pi * (accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	    }
	    else if (switch_KP == 2){
	      // K terms 
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_S_2pi * (temp_sumK);
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_S_2pi * (temp_sumK);
	    }
	  }
	}
      }

      for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
	for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	  if (switch_KP == 0){
	    // K + P term
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_S_2pi * (temp_sumK + accumulate(value_dataP[sd][ss][i_x][i_c].begin(), value_dataP[sd][ss][i_x][i_c].end(), 0.));
	  }
	  else if (switch_KP == 1){
	    // P terms 
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_S_2pi * (accumulate(value_dataP[sd][ss][i_x][i_c].begin(), value_dataP[sd][ss][i_x][i_c].end(), 0.));
	  }
	  else if (switch_KP == 2){
	    // K terms 
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_S_2pi * (temp_sumK);
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_collinear_QCD_CDST(){
  static Logger logger("observable_set::calculate_collinear_QCD_CDST");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  static int initialization = 1;
  // calculate all needed momentum-independent (splitting) functions
  // check if all elements are always zero!!!
  /*
  double Kbar[4][3] = {{0.}};
  double Kt[4][3] = {{0.}};
  double P[4][3] = {{0.}};
  double Kbar_plus[4][3] = {{0.}};
  double Kt_plus[4][3] = {{0.}};
  double P_plus[4][3] = {{0.}};
  double intKbar_plus[4][3] = {{0.}};
  double intKt_plus[4][3] = {{0.}};
  double intP_plus[4][3] = {{0.}};
  double Kbar_delta[4] = {0.};
  double Kt_delta[4] = {0.};
  double P_delta[4] = {0.};
  */
  // check if everything remains unchanged without setting all functions to 0 !!!
  // reconsider if different partonic channels are calculated simultaneously later !!!
  static double Kbar[4][3] = {{0.}};
  static double Kt[4][3] = {{0.}};
  static double P[4][3] = {{0.}};
  static double Kbar_plus[4][3] = {{0.}};
  static double Kt_plus[4][3] = {{0.}};
  static double P_plus[4][3] = {{0.}};
  static double intKbar_plus[4][3] = {{0.}};
  static double intKt_plus[4][3] = {{0.}};
  static double intP_plus[4][3] = {{0.}};
  static double Kbar_delta[4] = {0.};
  static double Kt_delta[4] = {0.};
  static double P_delta[4] = {0.};

  static double alpha_S_2pi = alpha_S * inv2pi;

  static vector<vector<int> > pair;
  static vector<double> iT2_ap(3);
  static vector<double> gamma_a_T2_ap(3);
  static vector<vector<double> > gamma_i_T2_i((*CA_collinear).size());
  static vector<vector<double> > ln_papi(3, vector<double> (csi->type_parton[0].size()));
  static vector<vector<vector<vector<double> > > > CA_value_ln_muF_papi(CA_value_log_mu2_fact.size());
  static vector<vector<vector<vector<vector<double> > > > > value_dataP(CA_value_log_mu2_fact.size());

  static vector<vector<int> > type((*CA_collinear).size());
  static vector<vector<int> > no_emitter((*CA_collinear).size());
  static vector<vector<int> > no_spectator((*CA_collinear).size());
  static vector<vector<int> > collinear_singularity((*CA_collinear).size());
  static vector<vector<vector<int> > > ppair((*CA_collinear).size());

  static int n_max_spectator = 0;
  if (initialization == 1){
    for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
      if ((*CA_collinear)[i_c].size() > n_max_spectator){n_max_spectator = (*CA_collinear)[i_c].size();}
    }
  }
  static vector<vector<vector<double> > > Kit(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > intKit_plus(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > Kit_plus_x(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > Kit_plus_1(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > Kit_plus_outside_x(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > Kit_plus_outside_1(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > Kit_delta(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<double> m_Q(csi->type_parton[0].size(), 0.);
  static vector<double> m2_Q(csi->type_parton[0].size(), 0.);
  static vector<vector<double> > sall_ja(3, vector<double> (csi->type_parton[0].size(), 0.));
  static vector<vector<double> > sall_ja_x(3, vector<double> (csi->type_parton[0].size(), 0.));
  static vector<vector<double> > mu2_Q(3, vector<double> (csi->type_parton[0].size(), 0.));
  static vector<vector<double> > mu2_Q_x(3, vector<double> (csi->type_parton[0].size(), 0.));

  if (initialization == 1){
    for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){
      if (i_p < 3 && M2[abs(csi->type_parton[0][i_p])] != 0.){cout << "Incoming massive partons are not supported!" << endl; int_end = 1;}
      if (M2[abs(csi->type_parton[0][i_p])] != 0.){
	//	m_Q[i_p] = M[abs(csi->type_parton[0][i_p])];
	//	m2_Q[i_p] = M2[abs(csi->type_parton[0][i_p])];
	m_Q[i_p] = mass_parton[0][i_p];
	m2_Q[i_p] = mass2_parton[0][i_p];
      }
      //      cout << "m2_Q[" << i_p << "] = " << m2_Q[i_p] << endl;
    }

    if (CA_dipole_splitting[0][1] == 1 || CA_dipole_splitting[0][2] == 1){    // g -> g (+g) splitting (0)
      Kbar_delta[0] = Kbar_gg_delta(N_f);
      Kt_delta[0] = Kt_gg_delta();
      P_delta[0] = P_gg_delta(N_f);
    }
    if (CA_dipole_splitting[1][1] == 1 || CA_dipole_splitting[1][2] == 1){    // q -> q (+g) splitting
      Kbar_delta[1] = Kbar_qq_delta();
      Kt_delta[1] = Kt_qq_delta();
      P_delta[1] = Pxx_qq_delta();
    }
  
    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      CA_value_ln_muF_papi[sd].resize(CA_value_log_mu2_fact[sd].size(), vector<vector<double> > (3, vector<double> (csi->type_parton[0].size())));
    }
    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      value_dataP[sd].resize(CA_value_log_mu2_fact[sd].size(), vector<vector<vector<double> > > (3, vector<vector<double> > ((*CA_collinear).size())));
      for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	for (int i_x = 0; i_x < 3; i_x++){
	  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
	    value_dataP[sd][ss][i_x][i_c].resize((*CA_collinear)[i_c].size(), 0.);
	  }
	}
      }
    }

    for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	if ((*CA_collinear)[i_c][j_c].no_spectator() == 0){continue;}
	int flag = -1;
	for (int i_p = 0; i_p < pair.size(); i_p++){
	  if ((*CA_collinear)[i_c][j_c].pair() == pair[i_p]){flag = i_p; break;}
	}
	if (flag == -1){pair.push_back((*CA_collinear)[i_c][j_c].pair());}
      }
    }
 
    for (int i_em = 1; i_em < 3; i_em++){
      if (csi->type_parton[0][i_em] == 0){iT2_ap[i_em] = 1. / C_A;}
      else {iT2_ap[i_em] = 1. / C_F;}
    }

    for (int i_em = 1; i_em < 3; i_em++){
      if (csi->type_parton[0][i_em] == 0){gamma_a_T2_ap[i_em] = gamma_g(N_f) / C_A;}
      else {gamma_a_T2_ap[i_em] = gamma_q / C_F;}
    }
    
    for (int i_c = 0; i_c < gamma_i_T2_i.size(); i_c++){ // i_c = -1 is not needed !!!
      gamma_i_T2_i[i_c].resize((*CA_collinear)[i_c].size());
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	if ((*CA_collinear)[i_c][j_c].type_parton()[(*CA_collinear)[i_c][j_c].no_spectator()] == 0){gamma_i_T2_i[i_c][j_c] = gamma_g(N_f) / C_A;}
	else {gamma_i_T2_i[i_c][j_c] = gamma_q / C_F;}
      }
    }

    for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
      type[i_c].resize((*CA_collinear)[i_c].size());
      no_emitter[i_c].resize((*CA_collinear)[i_c].size());
      no_spectator[i_c].resize((*CA_collinear)[i_c].size());
      collinear_singularity[i_c].resize((*CA_collinear)[i_c].size());
      ppair[i_c].resize((*CA_collinear)[i_c].size());
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	type[i_c][j_c] = (*CA_collinear)[i_c][j_c].type();
	no_emitter[i_c][j_c] = (*CA_collinear)[i_c][j_c].no_emitter();
	no_spectator[i_c][j_c] = (*CA_collinear)[i_c][j_c].no_spectator();
	collinear_singularity[i_c][j_c] = (*CA_collinear)[i_c][j_c].in_collinear()[0];
	ppair[i_c][j_c] = (*CA_collinear)[i_c][j_c].pair();
      }
    }


    for (int i_p = 0; i_p < pair.size(); i_p++){
      int pair_em = pair[i_p][0];
      int pair_sp = pair[i_p][1];

      if (m2_Q[pair_em] > 0.){return;}
      
      if (m2_Q[pair_sp] > 0.){
	for (int i_dt = 0; i_dt < 4; i_dt++){
	  Kit[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_x[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_1[i_dt][pair_em][pair_sp] = 0.;
	  intKit_plus[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_outside_x[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_outside_1[i_dt][pair_em][pair_sp] = 0.; 
	  Kit_delta[i_dt][pair_em][pair_sp] = 0.; 
	}
      }
    }

    for (int i_p = 0; i_p < pair.size(); i_p++){cout << "pair[" << i_p << "] = " << "(" << pair[i_p][0] << ", " << pair[i_p][1] << ")" << endl;}
    initialization = 0;
  }
  logger << LOG_DEBUG_VERBOSE << "initialization finished!" << endl;

  for (int i_em = 1; i_em < 3; i_em++){
    if (CA_dipole_splitting[0][i_em] == 1){      // g -> g (+g) splitting
      Kbar[0][i_em] = Kbar_gg(z_coll[i_em]);
      Kt[0][i_em] = Kt_gg(z_coll[i_em]);
      P[0][i_em] = P_gg(z_coll[i_em]);
      if (CA_dipole_splitting[0][0] == 1){
	Kbar_plus[0][i_em] = Kbar_gg_plus(z_coll[i_em]);
	Kt_plus[0][i_em] = Kt_gg_plus(z_coll[i_em]);
	P_plus[0][i_em] = P_gg_plus(z_coll[i_em]);
	intKbar_plus[0][i_em] = intKbar_gg_plus(x_pdf[i_em]);
	intKt_plus[0][i_em] = intKt_gg_plus(x_pdf[i_em]);
	intP_plus[0][i_em] = intP_gg_plus(x_pdf[i_em]);
      }
      else {cout << "May not happen !!! g -> g splitting without irregular terms !!!" << endl;}
    }
    if (CA_dipole_splitting[1][i_em] == 1){      // q -> q (+g) splitting
      Kbar[1][i_em] = Kbar_qq(z_coll[i_em]);
      Kt[1][i_em] = Kt_qq(z_coll[i_em]);
      P[1][i_em] = Pxx_qq(z_coll[i_em]);
      if (CA_dipole_splitting[1][0] == 1){
	Kbar_plus[1][i_em] = Kbar_qq_plus(z_coll[i_em]);
	Kt_plus[1][i_em] = Kt_qq_plus(z_coll[i_em]);
	P_plus[1][i_em] = Pxx_qq_plus(z_coll[i_em]);
	intKbar_plus[1][i_em] = intKbar_qq_plus(x_pdf[i_em]);
	intKt_plus[1][i_em] = intKt_qq_plus(x_pdf[i_em]);
	intP_plus[1][i_em] = intPxx_qq_plus(x_pdf[i_em]);
      }
      else {cout << "May not happen !!! q -> q splitting without irregular terms !!!" << endl;}
    }
    if (CA_dipole_splitting[2][i_em] == 1){      // q -> g (+q) splitting
      Kbar[2][i_em] = Kbar_qg(z_coll[i_em]);
      Kt[2][i_em] = Kt_qg(z_coll[i_em]);
      P[2][i_em] = P_qg(z_coll[i_em]);
    }
    if (CA_dipole_splitting[3][i_em] == 1){      // g -> q (+q~) splitting
      Kbar[3][i_em] = Kbar_gq(z_coll[i_em]);
      Kt[3][i_em] = Kt_gq(z_coll[i_em]);
      P[3][i_em] = P_gq(z_coll[i_em]);
    }
  }

  logger << LOG_DEBUG_VERBOSE << "splitting kernels finished!" << endl;


  for (int i_p = 0; i_p < pair.size(); i_p++){
    int pair_em = pair[i_p][0];
    int pair_sp = pair[i_p][1];
    // exception: 1 -- 2; however, in this case, for all relevant configurations the involved functions are symmetric!
    logger << LOG_DEBUG_VERBOSE << "p_parton[0].size() = " << p_parton[0].size() << endl;
    sall_ja[pair_em][pair_sp] = 2 * p_parton[0][pair_em] * p_parton[0][pair_sp];
    ln_papi[pair_em][pair_sp] = log(sall_ja[pair_em][pair_sp]);
    
    logger << LOG_DEBUG_VERBOSE << "m2_Q[" << pair_sp << "] = " << m2_Q[pair_sp] << endl;

    if (m2_Q[pair_sp] > 0.){
      // Only happens if no initial-initial dipole is discussed, i.e. pair_em and pair_sp really point at emitter and spectator, respectively.
      sall_ja_x[pair_em][pair_sp] = sall_ja[pair_em][pair_sp] / z_coll[pair_em];
      mu2_Q[pair_em][pair_sp] = m2_Q[pair_sp] / sall_ja[pair_em][pair_sp];
      mu2_Q_x[pair_em][pair_sp] = m2_Q[pair_sp] / sall_ja_x[pair_em][pair_sp];
      
      logger << LOG_DEBUG_VERBOSE << "i_p = " << i_p << "   pair_em = " << pair_em << "   pair_sp = " << pair_sp << "   z_coll[" << pair_em << "] = " << z_coll[pair_em] << endl;

      // only for massive quarks as spectators
      for (int i_dt = 0; i_dt < 4; i_dt++){
	if      (i_dt == 0 && (CA_dipole_splitting[i_dt][1] == 1 || CA_dipole_splitting[i_dt][2] == 1)){
	  // g -> g (+g) splitting
	  //	  cout << "g -> g (+g): dipole_phasespace[0][" << pair_sp + 6 << "] = " << dipole_phasespace[0][pair_sp + 6] << "   " << m_Q[pair_sp] << endl;
	  Kit[i_dt][pair_em][pair_sp] = 
	    - 2 * log(2. - z_coll[pair_em]) / (1. - z_coll[pair_em]) // from second term in (6.58) (-> K^qq_q) [included from (6.60)]
	    + 2 * m2_Q[pair_sp] / (z_coll[pair_em] * sall_ja_x[pair_em][pair_sp]) * log(m2_Q[pair_sp] / ((1. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp] + m2_Q[pair_sp])); // from (6.59) (-> K^qg_q) [included from (6.60)]
	  Kit_plus_x[i_dt][pair_em][pair_sp] = 
	    + 2 * log(1. - z_coll[pair_em]) / (1. - z_coll[pair_em]) // from first term in (6.58) 
	    + (1. - z_coll[pair_em]) / (2 * pow(1. - z_coll[pair_em] + mu2_Q_x[pair_em][pair_sp], 2)) // from first term in (5.58, J_gQ^a) included from (6.58) 
	    - 2 / (1. - z_coll[pair_em]) * (1. + log(1. - z_coll[pair_em] + mu2_Q_x[pair_em][pair_sp])) // from second term in (5.58, J_gQ^a) included from (6.58)
	    ;
	  Kit_plus_1[i_dt][pair_em][pair_sp] = 
	    + 2 * log(1. - z_coll[pair_em]) / (1. - z_coll[pair_em]) // from first term in (6.58) 
	    + (1. - z_coll[pair_em]) / (2 * pow(1. - z_coll[pair_em] + mu2_Q[pair_em][pair_sp], 2)) // from first term in (5.58, J_gQ^a) included from (6.58) 
	    - 2 / (1. - z_coll[pair_em]) * (1. + log(1. - z_coll[pair_em] + mu2_Q[pair_em][pair_sp])) // from second term in (5.58, J_gQ^a) included from (6.58)
	    ;
	  intKit_plus[i_dt][pair_em][pair_sp] = 
	    - pow(log(1. - x_pdf[pair_em]), 2) // from first term in (6.58)
	    + .5 * (- mu2_Q[pair_em][pair_sp] / (1. - x_pdf[pair_em] + mu2_Q[pair_em][pair_sp]) + mu2_Q[pair_em][pair_sp] / (1. + mu2_Q[pair_em][pair_sp]) - log((1. - x_pdf[pair_em] + mu2_Q[pair_em][pair_sp]) / (1. + mu2_Q[pair_em][pair_sp]))) // from first term in (5.58) included from (6.58)
	    + 2 * (gsl_sf_dilog(-1. / mu2_Q[pair_em][pair_sp]) - gsl_sf_dilog(-(1. - x_pdf[pair_em]) / mu2_Q[pair_em][pair_sp]) + log(1. - x_pdf[pair_em]) * (1. + log(mu2_Q[pair_em][pair_sp]))) // from second term in (5.58) included from (6.58)
	    ;
	  
	  // terms containing x-dependent pre-factor of (2/(1-z_coll[pair_em]))_+
	  Kit_plus_outside_x[i_dt][pair_em][pair_sp] = 
	    + log(((2. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp]) / ((2. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp] + m2_Q[pair_sp])) // from fourth term from (6.58)
	    + log(2. + mu2_Q_x[pair_em][pair_sp] - z_coll[pair_em]) // from third term in (5.58) included from (6.58)
	    ;
	  Kit_plus_outside_1[i_dt][pair_em][pair_sp] = 
	    + log(sall_ja[pair_em][pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp])) // from fourth term from (6.58)
	    + log(1. + mu2_Q[pair_em][pair_sp]) // from third term in (5.58) included from (6.58)
	    ;
	  Kit_delta[i_dt][pair_em][pair_sp] = -gamma_q / C_F // from fifth term from (6.58)
	    + mu2_Q[pair_em][pair_sp] * log(m2_Q[pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp])) // from sixth term from (6.58)
	    + .5 * m2_Q[pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp])// from seventh term from (6.58)
	    ;
	}
	else if (i_dt == 1 && (CA_dipole_splitting[i_dt][1] == 1 || CA_dipole_splitting[i_dt][2] == 1)){
	  // q -> q (+g) splitting
	  //	  cout << "q -> q (+g): dipole_phasespace[0][" << pair_sp + 6 << "] = " << dipole_phasespace[0][pair_sp + 6] << "   " << m_Q[pair_sp] << endl;
	  Kit[i_dt][pair_em][pair_sp] = -2 * log(2. - z_coll[pair_em]) / (1. - z_coll[pair_em]);
	  Kit_plus_x[i_dt][pair_em][pair_sp] = 0.
	    + 2 * log(1. - z_coll[pair_em]) / (1. - z_coll[pair_em]) // from first term in (6.58) 
	    + (1. - z_coll[pair_em]) / (2 * pow(1. - z_coll[pair_em] + mu2_Q_x[pair_em][pair_sp], 2)) // from first term in (5.58, J_gQ^a) included from (6.58) 
	    - 2 / (1. - z_coll[pair_em]) * (1. + log(1. - z_coll[pair_em] + mu2_Q_x[pair_em][pair_sp])) // from second term in (5.58, J_gQ^a) included from (6.58)
	    ;
	  
	  Kit_plus_1[i_dt][pair_em][pair_sp] = 0.
	    + 2 * log(1. - z_coll[pair_em]) / (1. - z_coll[pair_em]) // from first term in (6.58) 
	    + (1. - z_coll[pair_em]) / (2 * pow(1. - z_coll[pair_em] + mu2_Q[pair_em][pair_sp], 2)) // from first term in (5.58, J_gQ^a) included from (6.58) 
	    - 2 / (1. - z_coll[pair_em]) * (1. + log(1. - z_coll[pair_em] + mu2_Q[pair_em][pair_sp])) // from second term in (5.58, J_gQ^a) included from (6.58)
	    ;
	  
	  intKit_plus[i_dt][pair_em][pair_sp] = 0.
	    - pow(log(1. - x_pdf[pair_em]), 2) // from first term in (6.58)
	    + .5 * (- mu2_Q[pair_em][pair_sp] / (1. - x_pdf[pair_em] + mu2_Q[pair_em][pair_sp]) + mu2_Q[pair_em][pair_sp] / (1. + mu2_Q[pair_em][pair_sp]) - log((1. - x_pdf[pair_em] + mu2_Q[pair_em][pair_sp]) / (1. + mu2_Q[pair_em][pair_sp]))) // from first term in (5.58) included from (6.58)
	    + 2 * (gsl_sf_dilog(-1. / mu2_Q[pair_em][pair_sp]) - gsl_sf_dilog(-(1. - x_pdf[pair_em]) / mu2_Q[pair_em][pair_sp]) + log(1. - x_pdf[pair_em]) * (1. + log(mu2_Q[pair_em][pair_sp]))); // from second term in (5.58) included from (6.58)
	  
	  // terms containing x-dependent pre-factor of (2/(1-z_coll[pair_em]))_+
	  Kit_plus_outside_x[i_dt][pair_em][pair_sp] = 
	    + log(((2. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp]) / ((2. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp] + m2_Q[pair_sp])) // from fourth term from (6.58)
	    + log(2. + mu2_Q_x[pair_em][pair_sp] - z_coll[pair_em]) // from third term in (5.58) included from (6.58)
	    ;
	  Kit_plus_outside_1[i_dt][pair_em][pair_sp] = 
	    + log(sall_ja[pair_em][pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp])) // from fourth term from (6.58)
	    + log(1. + mu2_Q[pair_em][pair_sp]) // from third term in (5.58) included from (6.58)
	    ;
	  Kit_delta[i_dt][pair_em][pair_sp] = 
	    - gamma_q / C_F // from fifth term from (6.58)
	    + mu2_Q[pair_em][pair_sp] * log(m2_Q[pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp])) // from sixth term from (6.58)
	    + .5 * m2_Q[pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp]) // from seventh term from (6.58)
	    ;
	}
	else if (i_dt == 2 && (CA_dipole_splitting[i_dt][1] == 1 || CA_dipole_splitting[i_dt][2] == 1)){
	  // q -> g (+q) splitting
	  Kit[i_dt][pair_em][pair_sp] = 
	    2 * C_F / C_A * m2_Q[pair_sp] / (z_coll[pair_em] * sall_ja_x[pair_em][pair_sp]) * log(m2_Q[pair_sp] / ((1. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp] + m2_Q[pair_sp])); // from (6.59)
	  /*
	  Kit_plus_x[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_1[i_dt][pair_em][pair_sp] = 0.;
	  intKit_plus[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_outside_x[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_outside_1[i_dt][pair_em][pair_sp] = 0.; 
	  Kit_delta[i_dt][pair_em][pair_sp] = 0.; 
	  */
	}
	else if (i_dt == 3 && (CA_dipole_splitting[i_dt][1] == 1 || CA_dipole_splitting[i_dt][2] == 1)){
	  // g -> q (+g) splitting
	  //	  Kit[i_dt][pair_em][pair_sp] = 0.;
	  /*
	  Kit[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_x[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_1[i_dt][pair_em][pair_sp] = 0.;
	  intKit_plus[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_outside_x[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_outside_1[i_dt][pair_em][pair_sp] = 0.; 
	  Kit_delta[i_dt][pair_em][pair_sp] = 0.; 
	  */
	}
	/*
	else {
	  // ???
	  Kit[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_x[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_1[i_dt][pair_em][pair_sp] = 0.;
	  intKit_plus[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_outside_x[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_outside_1[i_dt][pair_em][pair_sp] = 0.; 
	  Kit_delta[i_dt][pair_em][pair_sp] = 0.; 
	}
	*/
      }
    }


    /*
    else if (dipole_phasespace[0][pair_sp + 6] > 2 && dx_pa[dipole_phasespace[0][pair_sp + 6]][0] == 0){
      // !!! does maybe not vanish in particular cases with outgoing gluons !!!
      // j == gluon case: additional terms for N_J^ja !!!
    }
    */

      
    logger << LOG_DEBUG_VERBOSE << "max_dyn_fact = " << max_dyn_fact << endl;
    logger << LOG_DEBUG_VERBOSE << "n_scale_dyn_fact.size() = " << n_scale_dyn_fact.size() << endl;

    if (switch_TSV){
      for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	//      logger << LOG_DEBUG_VERBOSE << "n_scale_dyn_fact[" << v_sf << "] = " << n_scale_dyn_fact[v_sf] << endl;
	logger << LOG_DEBUG << "value_central_logscale2_fact[" << v_sf << "] = " << value_central_logscale2_fact[v_sf] << endl;
	logger << LOG_DEBUG << "value_relative_logscale2_fact[" << v_sf << "].size() = " << value_relative_logscale2_fact[v_sf].size() << endl;

       	for (int v_xf = 0; v_xf < n_scale_dyn_fact[v_sf]; v_xf++){
	  logger << LOG_DEBUG << "value_relative_logscale2_fact[" << v_sf << "][" << v_xf << "] = " << value_relative_logscale2_fact[v_sf][v_xf] << endl;
	  logger << LOG_DEBUG << "value_absolute_logscale2_fact[" << v_sf << "][" << v_xf << "] = " << value_central_logscale2_fact[v_sf] + value_relative_logscale2_fact[v_sf][v_xf] << endl;
	  
	  logger << LOG_DEBUG << "value_logscale2_fact_papi.size() = " << value_logscale2_fact_papi.size() << endl;
	  logger << LOG_DEBUG << "value_logscale2_fact_papi[" << v_sf << "].size() = " << value_logscale2_fact_papi[v_sf].size() << endl;
	  logger << LOG_DEBUG << "value_logscale2_fact_papi[" << v_sf << "][" << v_xf << "].size() = " << value_logscale2_fact_papi[v_sf][v_xf].size() << endl;
	  
	  value_logscale2_fact_papi[v_sf][v_xf][pair_em][pair_sp] = value_central_logscale2_fact[v_sf] + value_relative_logscale2_fact[v_sf][v_xf] - ln_papi[pair_em][pair_sp];

	  logger << LOG_DEBUG << "value_logscale2_fact_papi[" << v_sf << "][" << v_xf << "][" << pair_em << "][" << pair_sp << "] = " << value_logscale2_fact_papi[v_sf][v_xf][pair_em][pair_sp] << endl;

	}
      }
    }

    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	CA_value_ln_muF_papi[sd][ss][pair_em][pair_sp] = CA_value_log_mu2_fact[sd][ss] - ln_papi[pair_em][pair_sp];
      }
    }
  }
  logger << LOG_DEBUG_VERBOSE << "pair finished!" << endl;

  /*
  for (int i_dt = 0; i_dt < 4; i_dt++){
    cout << "CA_dipole_splitting[" << i_dt << "] = ";
    for (int i_em = 0; i_em < 3; i_em++){
      cout << CA_dipole_splitting[i_dt][i_em] << "   ";
    }
    cout << endl;
  }

  for (int i_dt = 0; i_dt < 4; i_dt++){
    for (int i_p = 0; i_p < pair.size(); i_p++){
      int pair_em = pair[i_p][0];
      int pair_sp = pair[i_p][1];
      
      cout << "Kit                [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit[i_dt][pair_em][pair_sp] << endl;
      cout << "Kit_plus_x         [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit_plus_x[i_dt][pair_em][pair_sp] << endl;
      cout << "Kit_plus_1         [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit_plus_1[i_dt][pair_em][pair_sp] << endl;
      cout << "intKit_plus        [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << intKit_plus[i_dt][pair_em][pair_sp] << endl;
      cout << "Kit_plus_outside_x [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit_plus_outside_x[i_dt][pair_em][pair_sp] << endl;
      cout << "Kit_plus_outside_1 [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit_plus_outside_1[i_dt][pair_em][pair_sp] << endl;
      cout << "Kit_delta          [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit_delta[i_dt][pair_em][pair_sp] << endl;
      
    }
  }
    */

  /*
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
      cout << "CA_ME2_cf[" << i_c << "][" << j_c << "] = " << CA_ME2_cf[i_c][j_c] << endl;
    }
  }
  */

  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
      // K terms 
      if (no_spectator[i_c][j_c] == 0){
	data_K[0][i_c][j_c] = (Kbar[type[i_c][j_c]][no_emitter[i_c][j_c]] + Kbar_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
	data_K[1][i_c][j_c] = (-Kbar_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];
	data_K[2][i_c][j_c] = (Kbar_delta[type[i_c][j_c]] - intKbar_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
      }
      else if (no_spectator[i_c][j_c] < 3){
	data_K[0][i_c][j_c] = -iT2_ap[no_emitter[i_c][j_c]] * (Kt[type[i_c][j_c]][no_emitter[i_c][j_c]] + Kt_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];// reg
	data_K[1][i_c][j_c] = -iT2_ap[no_emitter[i_c][j_c]] * (-Kt_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];// plus
	data_K[2][i_c][j_c] = -iT2_ap[no_emitter[i_c][j_c]] * (Kt_delta[type[i_c][j_c]] - intKt_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];// delta
      }
      else if (no_spectator[i_c][j_c] > 2){
	//	cout << "m2_Q[no_spectator[" << i_c << "][" << j_c << "] = " << no_spectator[i_c][j_c] << "] = " << m2_Q[no_spectator[i_c][j_c]] << endl;
	if (m2_Q[no_spectator[i_c][j_c]] == 0.){
	  if (collinear_singularity[i_c][j_c] == 1){
	    data_K[0][i_c][j_c] = gamma_i_T2_i[i_c][j_c] * (1. / (1. - z_coll[no_emitter[i_c][j_c]])) * CA_ME2_cf[i_c][j_c];
	    data_K[1][i_c][j_c] = -gamma_i_T2_i[i_c][j_c] * (1. / (1. - z_coll[no_emitter[i_c][j_c]])) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];
	    data_K[2][i_c][j_c] = gamma_i_T2_i[i_c][j_c] * (1. + log(1. - x_pdf[no_emitter[i_c][j_c]])) * CA_ME2_cf[i_c][j_c];
	  }
	}
	else {
	  // -> Kit
	  //	  cout << "A data_K[0][" << i_c << "][" << j_c << "] = " << data_K[0][i_c][j_c] << endl;
	  data_K[0][i_c][j_c] = -(Kit[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] 
					+ Kit_plus_x[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] 
					+ 2. / (1. - z_coll[no_emitter[i_c][j_c]]) * Kit_plus_outside_x[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]
					) * CA_ME2_cf[i_c][j_c]; // reg
	  /*
	  cout << "  Kit[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]  = " << Kit[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]  << endl;
	  cout << "  Kit_plus_x[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]  = " << Kit_plus_x[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]  << endl;
	  cout << "  2. / (1. - z_coll[no_emitter[i_c][j_c]]) * Kit_plus_outside_x[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] = " << 2. / (1. - z_coll[no_emitter[i_c][j_c]]) * Kit_plus_outside_x[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] << endl;
	  cout << "  CA_ME2_cf[i_c][j_c] = " << CA_ME2_cf[i_c][j_c] << endl;
	  cout << "B data_K[0][" << i_c << "][" << j_c << "] = " << data_K[0][i_c][j_c] << endl;
	  cout << "A data_K[1][" << i_c << "][" << j_c << "] = " << data_K[1][i_c][j_c] << endl;
	  */
	  data_K[1][i_c][j_c] =  -(
					 - Kit_plus_1[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] 
					 - 2. / (1. - z_coll[no_emitter[i_c][j_c]]) * Kit_plus_outside_1[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]
					 ) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c]; // plus
	  /*
	  cout << "B data_K[1][" << i_c << "][" << j_c << "] = " << data_K[1][i_c][j_c] << endl;
	  cout << "A data_K[2][" << i_c << "][" << j_c << "] = " << data_K[2][i_c][j_c] << endl;
	  */
	  data_K[2][i_c][j_c] = -(Kit_delta[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] 
					- intKit_plus[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]
					- 2 * log(1. - x_pdf[no_emitter[i_c][j_c]]) * Kit_plus_outside_1[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]
					) * CA_ME2_cf[i_c][j_c]; // delta
	  //	  cout << "B data_K[2][" << i_c << "][" << j_c << "] = " << data_K[2][i_c][j_c] << endl;
	  
	  data_K[0][i_c][j_c] += -iT2_ap[no_emitter[i_c][j_c]] * P[type[i_c][j_c]][no_emitter[i_c][j_c]] * log(((1. - z_coll[no_emitter[i_c][j_c]]) * sall_ja_x[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]) / ((1. - z_coll[no_emitter[i_c][j_c]]) * sall_ja_x[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] + m2_Q[no_spectator[i_c][j_c]])) * CA_ME2_cf[i_c][j_c]; // from (6.55)
	  //	  cout << "C data_K[0][" << i_c << "][" << j_c << "] = " << data_K[0][i_c][j_c] << endl;
	  // no explizit plus term in (6.55)
	  //	  cout << "collinear_singularity[" << i_c << "][" << j_c << "] = " << collinear_singularity[i_c][j_c] << endl;
	  if (collinear_singularity[i_c][j_c] == 1){
	    // contributes only to irregular splittings 
	    //	    cout << "data_K[2][" << i_c << "][" << j_c << "] = " << data_K[2][i_c][j_c] << endl;
	    data_K[2][i_c][j_c] += -gamma_a_T2_ap[no_emitter[i_c][j_c]] * (log((sall_ja[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] - 2 * m_Q[no_spectator[i_c][j_c]] * sqrt(sall_ja[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] + m2_Q[no_spectator[i_c][j_c]]) + 2 * m2_Q[no_spectator[i_c][j_c]]) / sall_ja[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]) + 2 * m_Q[no_spectator[i_c][j_c]] / (sqrt(sall_ja[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] + m2_Q[no_spectator[i_c][j_c]]) + m_Q[no_spectator[i_c][j_c]])) * CA_ME2_cf[i_c][j_c];  // from (6.55)
	    //	    data_K[2][i_c][i_ca] += -gamma_a_T2_ap[i_em] * (log((sall_ja[i_em][i_cs] - 2 * m_j * sqrt(sall_ja[i_em][i_cs] + m2_j) + 2 * m2_j) / sall_ja[i_em][i_cs]) + 2 * m_j / (sqrt(sall_ja[i_em][i_cs] + m2_j) + m_j)) * CA_ME2cc[i_c][i_ca];  // from (6.55)
	    //	    cout << "data_K[2][" << i_c << "][" << j_c << "] = " << data_K[2][i_c][j_c] << endl;

	  }
	  //	  cout << "C data_K[2][" << i_c << "][" << j_c << "] = " << data_K[2][i_c][j_c] << endl;
	}
      }
      if (no_spectator[i_c][j_c] != 0){
	// P terms 
	if (switch_TSV){
	  for (int v_sf = 0; v_sf < value_logscale2_fact_papi.size(); v_sf++){
	    for (int v_xf = 0; v_xf < value_logscale2_fact_papi[v_sf].size(); v_xf++){
	      value_data_P[v_sf][v_xf][0][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (P[type[i_c][j_c]][no_emitter[i_c][j_c]] + P_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	      value_data_P[v_sf][v_xf][1][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (-P_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	      value_data_P[v_sf][v_xf][2][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (P_delta[type[i_c][j_c]] - intP_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	    }
	  }
	}
	for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
	  for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	    value_dataP[sd][ss][0][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (P[type[i_c][j_c]][no_emitter[i_c][j_c]] + P_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	    value_dataP[sd][ss][1][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (-P_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	    value_dataP[sd][ss][2][i_c][j_c] = iT2_ap[no_emitter[i_c][j_c]] * (P_delta[type[i_c][j_c]] - intP_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	  }
	}
      }
    }
    for (int i_x = 0; i_x < 3; i_x++){
      /*
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	cout << "data_K[" << i_x << "][" << i_c << "][" << j_c << "] = " << data_K[i_x][i_c][j_c] << endl;
      }
      */
      double temp_sumK = accumulate(data_K[i_x][i_c].begin(), data_K[i_x][i_c].end(), 0.);
      logger << LOG_DEBUG_VERBOSE << "i_x = " << i_x << "   temp_sumK = " << temp_sumK << endl;
      if (switch_TSV){
	for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	  for (int v_xf = 0; v_xf < value_logscale2_fact_papi[v_sf].size(); v_xf++){
	    if (switch_KP == 0){
	      // K + P terms
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_S_2pi * (temp_sumK + accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_S_2pi * (temp_sumK + accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	    }
	    else if (switch_KP == 1){
	      // P terms 
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_S_2pi * (accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_S_2pi * (accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	    }
	    else if (switch_KP == 2){
	      // K terms 
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_S_2pi * (temp_sumK);
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_S_2pi * (temp_sumK);
	    }
	  }
	}
      }
      //    }
      //    dataK = data_K;
      //    for (int i_x = 0; i_x < 3; i_x++){
      for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
	for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	  if (switch_KP == 0){
	    // K + P term
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_S_2pi * (temp_sumK + accumulate(value_dataP[sd][ss][i_x][i_c].begin(), value_dataP[sd][ss][i_x][i_c].end(), 0.));
	  }
	  else if (switch_KP == 1){
	    // P terms 
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_S_2pi * (accumulate(value_dataP[sd][ss][i_x][i_c].begin(), value_dataP[sd][ss][i_x][i_c].end(), 0.));
	  }
	  else if (switch_KP == 2){
	    // K terms 
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_S_2pi * (temp_sumK);
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




void observable_set::determine_collinear_QEW(phasespace_set & psi){
  Logger logger("observable_set::determine_collinear_QEW");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static map<int, double> charge_particle;
  static vector<double> charge2_particle(26);
  if (initialization == 1){
    fill_charge_particle(charge_particle);
    for (int i_p = 0; i_p < 26; i_p++){
      charge2_particle[i_p] = pow(charge_particle[i_p], 2);
    }
    initialization = 0;
  }
  logger << LOG_DEBUG << "new collinear-dipole determination" << endl << endl;
  for (int i1 = 0; i1 < combination_pdf.size(); i1++){
    stringstream temp_ss;
    temp_ss << " combination_pdf[" << setw(2) << i1 << "] = "; 
    for (int i = 0; i < 3; i++){temp_ss << setw(4) << combination_pdf[i1][i] << "   ";}
    logger << LOG_DEBUG << temp_ss.str() << endl;
  }
  logger.newLine(LOG_DEBUG);

  vector<string> pa_name(csi->type_parton[0].size(), "");
  // !!! Should remove gluons here !!!
  if (charge_particle[csi->type_parton[0][1]] != 0 || csi->type_parton[0][1] == 22){pa_name[1] = "a";}
  if (charge_particle[csi->type_parton[0][2]] != 0 || csi->type_parton[0][2] == 22){pa_name[2] = "b";}
  //  pa_name[1] = "a";
  //  pa_name[2] = "b";
  //  int count = 0;
  //  vector<string> alphabet(csi->type_parton[0].size() - 3, "");
  //  for (int i_p = 0; i_p < alphabet.size(); i_p++){alphabet[i_p] = char(105 + i_p);}
  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (charge_particle[csi->type_parton[0][i_p]] != 0 || csi->type_parton[0][i_p] == 22){pa_name[i_p] = char(105 + i_p - 3);}}
  //  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (charge_particle[csi->type_parton[0][i_p]] != 0 || csi->type_parton[0][i_p] == 22){pa_name[i_p] = alphabet[i_p - 3]; count++;}}
  //  for (int i_p = 0; i_p < alphabet.size(); i_p++){alphabet[i_p] = char(105 + i_p);}
  //  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (charge_particle[csi->type_parton[0][i_p]] != 0 || csi->type_parton[0][i_p] == 22){pa_name[i_p] = alphabet[count++];}}

  stringstream temp_tp;
  stringstream temp_name;
  for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){
    temp_tp << setw(3) << csi->type_parton[0][i_p] << " ";
    temp_name << setw(3) << pa_name[i_p] << " ";
  }
  logger.newLine(LOG_INFO);
  logger << LOG_INFO << "Collinear QED emission:" << endl;
  logger << LOG_INFO << "csi->type_parton[0] = " << temp_tp.str() << endl;
  logger << LOG_INFO << "pa_name             = " << temp_name.str() << endl;



  int temp_type_correction = 2;

  /*
  map <int,string> pname;
  fill_pname(pname);
  */

  for (int temp_no_emitter = 1; temp_no_emitter < 3; temp_no_emitter++){
    if (pa_name[temp_no_emitter] == ""){continue;}
    string temp_name;
    vector<string> temp_all_name(combination_pdf.size());
    int temp_type;
    vector<int> temp_in_collinear(3, 0);
    // to be checked !!!
    //    temp_in_collinear[temp_no_emitter] = 1;
    vector<vector<int> > temp_pdf_new;
    double temp_charge_factor = 1.;







    for (int i_t = 0; i_t < 2; i_t++){
      // to be checked !!!
      temp_in_collinear[temp_no_emitter] = 1;

      temp_pdf_new = combination_pdf;
      //      if (combination_pdf[0][temp_no_emitter] == 22){
      if (combination_pdf[0][temp_no_emitter] == 7){
	if (i_t == 0){
	  //	  cout << "temp_no_emitter = " << temp_no_emitter << "   i_t = " << i_t << endl;
	  //	  cout << "LHAPDF::hasPhoton() = " << LHAPDF::hasPhoton() << endl;
	  temp_type = 0; // hard process with a from a -> a a splitting (not in QEW)
	  temp_in_collinear[0] = 1;
	  // to be checked !!!
	  temp_in_collinear[temp_no_emitter] = 0;

	  temp_name = "C_" + pa_name[temp_no_emitter] + "^{aa}";
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_all_name[i_x] = "C_" + pa_name[temp_no_emitter] + "^{aa}";}
	}
	else {
	  temp_type = 2; // hard process with a from q -> q a splitting
	  temp_in_collinear[0] = 0;
	  temp_name = "C_" + pa_name[temp_no_emitter] + "^{qa}";
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_pdf_new[i_x][temp_no_emitter] = 10;}
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_all_name[i_x] = "C_" + pa_name[temp_no_emitter] + "^{qa}";}
	}
      }
      //      else if (combination_pdf[0][temp_no_emitter] != 22){
      else if (combination_pdf[0][temp_no_emitter] != 7){
	if (i_t == 0){
	  //	  cout << "temp_no_emitter = " << temp_no_emitter << "   i_t = " << i_t << endl;
	  //	  cout << "LHAPDF::hasPhoton() = " << LHAPDF::hasPhoton() << endl;
	  temp_type = 3; // hard process with q/qx from a -> q qx splitting
	  temp_in_collinear[0] = 0;
	  temp_name = "C_" + pa_name[temp_no_emitter] + "^{a" + csi->name_particle[combination_pdf[0][temp_no_emitter]] + "}";
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_pdf_new[i_x][temp_no_emitter] = 7;}
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_all_name[i_x] = "C_" + pa_name[temp_no_emitter] + "^{a" + csi->name_particle[combination_pdf[i_x][temp_no_emitter]] + "}";}
	}
	else {
	  temp_type = 1; // hard process with q from q -> q a splitting
	  temp_in_collinear[0] = 1;
	  temp_name = "C_" + pa_name[temp_no_emitter] + "^{" + csi->name_particle[combination_pdf[0][temp_no_emitter]] + csi->name_particle[combination_pdf[0][temp_no_emitter]] + "}";
	  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){temp_all_name[i_x] = "C_" + pa_name[temp_no_emitter] + "^{" + csi->name_particle[combination_pdf[i_x][temp_no_emitter]] + csi->name_particle[combination_pdf[i_x][temp_no_emitter]] + "}";}
	}
      }

      if ((LHAPDF::hasPhoton() == 0) && (temp_type == 0 || temp_type == 3 || csi->type_parton[0][temp_no_emitter % 2 + 1] == 22)){
	logger << LOG_DEBUG << "No photon pdf available in pdf set " << LHAPDFname << endl;
	continue;
      }
      /*
      if ((LHAPDF::hasPhoton() == 0) && (temp_type == 0 || temp_type == 3 || csi->type_parton[0][temp_no_emitter % 2 + 1] == 22)){
	logger << LOG_DEBUG << "No photon pdf available in pdf set " << LHAPDFname << endl;
	continue;
      }
      */
      int flag = (*CA_collinear).size();
      for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){if (temp_no_emitter == (*CA_collinear)[i_c][0].no_emitter() && temp_type == (*CA_collinear)[i_c][0].type()){flag = i_c; break;}}
      if (flag == (*CA_collinear).size()){(*CA_collinear).push_back(vector<collinear_set> ());}
      for (int temp_no_spectator = 0; temp_no_spectator < csi->type_parton[0].size(); temp_no_spectator++){
	if ((temp_no_spectator > 0 && pa_name[temp_no_spectator] == "") || temp_no_spectator == temp_no_emitter){continue;}
	vector<int> temp_pair(2);
	if (temp_no_emitter < temp_no_spectator){
	  temp_pair[0] = temp_no_emitter;
	  temp_pair[1] = temp_no_spectator;
	}
	else {
	  temp_pair[0] = temp_no_spectator;
	  temp_pair[1] = temp_no_emitter;
	}

	int temp_massive;
	if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 0;}
	else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 1;}
	else if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 2;}
	else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 3;}
	else {cout << "Should not happen!" << endl;}


	//    double temp_charge_factor;
	logger << LOG_DEBUG << "temp_no_emitter = " << temp_no_emitter << endl;
	logger << LOG_DEBUG << "csi->type_parton[0][temp_no_emitter   = " << temp_no_emitter << "] = " << csi->type_parton[0][temp_no_emitter] << endl;
	logger << LOG_DEBUG << "csi->type_parton[0][temp_no_spectator = " << temp_no_spectator << "] = " << csi->type_parton[0][temp_no_spectator] << endl;
	if (temp_no_spectator == 0){
	  temp_charge_factor = 1.;
	  //	  temp_charge_factor = pow(charge_particle[abs(csi->type_parton[0][temp_no_emitter])], 2);
	}
	else if (csi->type_parton[0][temp_no_emitter] != 22 && csi->type_parton[0][temp_no_spectator] != 22){
	  temp_charge_factor = charge_particle[abs(csi->type_parton[0][temp_no_spectator])] / charge_particle[abs(csi->type_parton[0][temp_no_emitter])]; // (Q_ij Q_k) / Q²_ij
	  // seems to be the other way round, which looks also more reasonable... // before: Q²_ij / (Q_ij Q_k)
	  //	  temp_charge_factor = charge_particle[abs(csi->type_parton[0][temp_no_emitter])] * charge_particle[abs(csi->type_parton[0][temp_no_spectator])];
	  if ((temp_no_emitter > 2 && csi->type_parton[0][temp_no_emitter] > 0) || (temp_no_emitter < 3 && csi->type_parton[0][temp_no_emitter] < 0)){temp_charge_factor = -temp_charge_factor;}
	  // * (-1) for outgoing quark / incoming antiquark as emitter (only the latter makes sense here !!!)
	  if ((temp_no_spectator > 2 && csi->type_parton[0][temp_no_spectator] > 0) || (temp_no_spectator < 3 && csi->type_parton[0][temp_no_spectator] < 0)){temp_charge_factor = -temp_charge_factor;}
	  // * (-1) for outgoing quark / incoming antiquark as spectator
	}
	else if (csi->type_parton[0][temp_no_emitter] == 22){
	  if (temp_no_spectator < 3){temp_charge_factor = -1.;} // kappa_ij,k = -1 for the other initial-state particle, 0 elsewhere
	  else {temp_charge_factor = 0.;}
	  //	  temp_charge_factor = 1.;
	  // charge factor should depend on the respective charge of the splitting quarks/antiquarks - not implemented in a first version in the pdf routines!!!
	  /*
	  temp_charge_factor = charge_particle[abs(csi->type_parton[0][dipole_candidate[i_a].no_R_emitter_1()])] * charge_particle[abs(csi->type_parton[0][dipole_candidate[i_a].no_R_emitter_2()])];
	  if ((dipole_candidate[i_a].no_R_emitter_1() > 2 && csi->type_parton[0][dipole_candidate[i_a].no_R_emitter_1()] > 0) || (dipole_candidate[i_a].no_R_emitter_1() < 3 && csi->type_parton[0][dipole_candidate[i_a].no_R_emitter_1()] < 0)){temp_charge_factor = -temp_charge_factor;}
	  if ((dipole_candidate[i_a].no_R_emitter_2() > 2 && csi->type_parton[0][dipole_candidate[i_a].no_R_emitter_2()] > 0) || (dipole_candidate[i_a].no_R_emitter_2() < 3 && csi->type_parton[0][dipole_candidate[i_a].no_R_emitter_2()] < 0)){temp_charge_factor = -temp_charge_factor;}
	  */
	}
	else {temp_charge_factor = 0.;}
	// case of outgoing photons not correctly implemented !!!

	logger << LOG_DEBUG << "temp_charge_factor = " << temp_charge_factor << endl;


	vector<string> temp_all_name_spectator = temp_all_name;
	for (int i_x = 0; i_x < temp_all_name_spectator.size(); i_x++){temp_all_name_spectator[i_x] = temp_all_name_spectator[i_x] + "(" + pa_name[temp_no_spectator] + ")";}
	string temp_name_spectator = temp_name + "(" + pa_name[temp_no_spectator] + ")";

	// check if CA_collinear-type gives a non-vanishing contribution:
	//	if (temp_type != 0 && temp_charge_factor != 0.){
	(*CA_collinear)[flag].push_back(collinear_set(temp_name_spectator, temp_all_name_spectator, temp_type, temp_in_collinear, psi_no_prc[0], csi->type_parton[0], temp_pdf_new, temp_charge_factor, temp_no_emitter, temp_no_spectator, temp_pair, temp_type_correction, temp_massive));
	  //	}
	// pdf-selection effects are included later !!!
      }
    }
  }

  //  exit(1);

  logger.newLine(LOG_INFO);
  logger << LOG_INFO << "Before selection of contributing collinear 'dipoles':" << endl;
  logger.newLine(LOG_INFO);

  output_collinear();

  logger.newLine(LOG_INFO);
  logger << LOG_INFO << "Before filling CA_combination_pdf for collinear 'dipoles':" << endl;
  logger.newLine(LOG_INFO);
  CA_combination_pdf.resize((*CA_collinear).size(), vector<vector<vector<int> > > (combination_pdf.size()));
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      // CA_combination_pdf[i_c] contains a vector, which contains the usual combination_pdf(3) with 0 -> direction (1, -1), 1 -> parton with x1 (in hadron 1/2 for +1/-1), 2 -> parton with x2 (in hadron 2/1 for +1/-1)
      if (((*CA_collinear)[i_c][0].all_pdf()[i_i][1] == 10) && ((*CA_collinear)[i_c][0].all_pdf()[i_i][2] == 10)){
	for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	  if (i_q == 0){continue;}
	  for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	    if (j_q == 0){continue;}
	    vector<int> new_temp_combination_pdf = (*CA_collinear)[i_c][0].all_pdf()[i_i];
	    new_temp_combination_pdf[1] = i_q;
	    new_temp_combination_pdf[2] = j_q;
	    CA_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	  }
	}
      }
      else if ((*CA_collinear)[i_c][0].all_pdf()[i_i][1] == 10){
	for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	  if (i_q == 0){continue;}
	  vector<int> new_temp_combination_pdf = (*CA_collinear)[i_c][0].all_pdf()[i_i];
	  new_temp_combination_pdf[1] = i_q;
	  CA_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	}
      }
      else if ((*CA_collinear)[i_c][0].all_pdf()[i_i][2] == 10){
	for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	  if (j_q == 0){continue;}
	  vector<int> new_temp_combination_pdf = (*CA_collinear)[i_c][0].all_pdf()[i_i];
	  new_temp_combination_pdf[2] = j_q;
	  CA_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	}
      }
      else {
	CA_combination_pdf[i_c][i_i].push_back((*CA_collinear)[i_c][0].all_pdf()[i_i]);
      }
    }
  }
  logger.newLine(LOG_INFO);
  logger << LOG_INFO << "After filling CA_combination_pdf for collinear 'dipoles'." << endl;



  logger << LOG_INFO << "Before filling CA_Q2f for collinear 'dipoles':" << endl;
  logger.newLine(LOG_INFO);
  CA_Q2f.resize((*CA_collinear).size());
  logger << LOG_INFO << "(*CA_collinear).size() = " << (*CA_collinear).size() << endl;
  logger << LOG_INFO << "CA_combination_pdf.size() = " << CA_combination_pdf.size() << endl;
  logger << LOG_INFO << "CA_Q2f.size() = " <<  CA_Q2f.size() << endl;
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    //    CA_Q2f[i_c].resize((*CA_collinear)[i_c].size());
    //    for (int i_i = 0; i_i < (*CA_collinear)[i_c].size(); i_i++){
    CA_Q2f[i_c].resize(CA_combination_pdf[i_c].size());
    logger << LOG_INFO << "(*CA_collinear)[" << i_c << "].size() = " << (*CA_collinear)[i_c].size() << endl;
    logger << LOG_INFO << "CA_combination_pdf[" << i_c << "].size() = " << CA_combination_pdf[i_c].size() << endl;
    logger << LOG_INFO << "CA_Q2f[" << i_c << "].size() = " << CA_Q2f[i_c].size() << endl;
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      logger << LOG_INFO << "CA_combination_pdf[" << i_c << "][" << i_i << "].size() = " << CA_combination_pdf[i_c][i_i].size() << endl;
      logger << LOG_INFO << "CA_Q2f[" << i_c << "][" << i_i << "].size() = " << CA_Q2f[i_c][i_i].size() << endl;
      for (int i_x = 0; i_x < CA_combination_pdf[i_c][i_i].size(); i_x++){
	//	logger << LOG_INFO << "CA_combination_pdf[" << i_c << "][" << i_i << "][" << i_x << "].size() = " << CA_combination_pdf[i_c][i_i][i_x].size() << endl;
	//(*CA_collinear)[i_c][i_x] -> (*CA_collinear)[i_c][i_i]
	logger << LOG_DEBUG << "i_c = " << i_c << "   i_i = " << i_i << "   i_x = " << i_x << endl;
	if ((*CA_collinear)[i_c][i_i].type_parton()[(*CA_collinear)[i_c][i_i].no_emitter()] == 22){  // emitter is a photon
	  if ((*CA_collinear)[i_c][i_i].type() == 0){  // hard process with a from a -> a a splitting
	    CA_Q2f[i_c][i_i].push_back(1.);
	    logger << LOG_INFO << "type = " << (*CA_collinear)[i_c][i_i].type() << "   CA_Q2f[" << i_c << "][" << i_i << "][" << CA_Q2f[i_c][i_i].size() - 1 << "] = " << CA_Q2f[i_c][i_i][CA_Q2f[i_c][i_i].size() - 1] << endl;
	  }
	  else if ((*CA_collinear)[i_c][i_i].type() == 2){  // hard process with a from q -> q a splitting
	    logger << LOG_INFO << "CA_combination_pdf[i_c = " << i_c << "][i_i = " << i_i << "][i_x = " << i_x << "][(*CA_collinear)[i_c = " << i_c << "][i_i = " << i_i << "].no_emitter() = " << (*CA_collinear)[i_c][i_i].no_emitter() << "] = " << CA_combination_pdf[i_c][i_i][i_x][(*CA_collinear)[i_c][i_i].no_emitter()] << endl;
	    CA_Q2f[i_c][i_i].push_back(charge2_particle[abs(CA_combination_pdf[i_c][i_i][i_x][(*CA_collinear)[i_c][i_i].no_emitter()])]);
	    logger << LOG_INFO << "type = " << (*CA_collinear)[i_c][i_i].type() << "   CA_Q2f[" << i_c << "][" << i_i << "][" << CA_Q2f[i_c][i_i].size() - 1 << "] = " << CA_Q2f[i_c][i_i][CA_Q2f[i_c][i_i].size() - 1] << endl;
	  }
	  else {logger << LOG_WARN << "Wrong splitting: may not happen!" << endl; exit(1);}
	}
	else { // emitter and spectator are charged particles
	  CA_Q2f[i_c][i_i].push_back(charge2_particle[abs((*CA_collinear)[i_c][i_i].type_parton()[(*CA_collinear)[i_c][i_i].no_emitter()])]);
	  logger << LOG_INFO << "type = " << (*CA_collinear)[i_c][i_i].type() << "   CA_Q2f[" << i_c << "][" << i_i << "][" << CA_Q2f[i_c][i_i].size() - 1 << "] = " << CA_Q2f[i_c][i_i][CA_Q2f[i_c][i_i].size() - 1] << endl;
	}
	logger << LOG_DEBUG << "i_c = " << i_c << "   i_i = " << i_i << "   i_x = " << i_x << endl;
      }
      logger << LOG_INFO << "CA_Q2f[" << i_c << "][" << i_i << "].size() = " << CA_Q2f[i_c][i_i].size() << endl;
    }
  }
  logger << LOG_INFO << "After filling CA_Q2f for collinear 'dipoles'." << endl;
  logger.newLine(LOG_INFO);



  logger << LOG_INFO << "Before checking if some CA_collinear gives  vanishing contribution due to charge_factor:" << endl;
  logger.newLine(LOG_INFO);
  logger << LOG_DEBUG << "(*CA_collinear).size() = " << (*CA_collinear).size() << endl;
  logger << LOG_DEBUG << "CA_Q2f.size() = " << CA_Q2f.size() << endl;
  for (int i_c = (*CA_collinear).size() - 1; i_c >=0; i_c--){
    logger << LOG_DEBUG << "(*CA_collinear)[" << i_c << "].size() = " << (*CA_collinear)[i_c].size() << endl;
    logger << LOG_DEBUG << "CA_Q2f[" << i_c << "].size() = " << CA_Q2f[i_c].size() << endl;
    for (int j_c = (*CA_collinear)[i_c].size() - 1; j_c >= 0; j_c--){
      if ((*CA_collinear)[i_c][j_c].charge_factor() == 0.){
	(*CA_collinear)[i_c].erase((*CA_collinear)[i_c].begin() + j_c);
      }
    }
    if ((*CA_collinear)[i_c].size() == 0){
      (*CA_collinear).erase((*CA_collinear).begin() + i_c);
      CA_Q2f.erase(CA_Q2f.begin() + i_c); 
    }
  }
  logger.newLine(LOG_INFO);
  logger << LOG_INFO << "After checking if some CA_collinear gives  vanishing contribution due to charge_factor." << endl;


  // pdf-selection effects are included later !!!

  logger.newLine(LOG_INFO);
  logger << LOG_INFO << "After selection of contributing collinear 'dipoles':" << endl;
  logger.newLine(LOG_INFO);
  ///  output_collinear();

  CA_dipole_splitting.resize(4, vector<int> (3, 0));
  
  //  logger << LOG_DEBUG_VERBOSE << "(*CA_collinear).size() = " << (*CA_collinear).size() << endl;
  //  logger << LOG_DEBUG_VERBOSE << "CA_dipole_splitting.size() = " << CA_dipole_splitting.size() << endl;
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    //    logger << LOG_DEBUG_VERBOSE << "(*CA_collinear)[" << i_c << "][0].in_collinear().size() = " << (*CA_collinear)[i_c][0].in_collinear().size() << endl;
    //    logger << LOG_DEBUG_VERBOSE << "CA_dipole_splitting[(*CA_collinear)[" << i_c << "].size() = " << (*CA_collinear)[i_c].size() << endl;
    //    logger << LOG_DEBUG_VERBOSE << "CA_dipole_splitting[(*CA_collinear)[i_c][0].type() = " << (*CA_collinear)[i_c][0].type() << "] = " << CA_dipole_splitting[(*CA_collinear)[i_c][0].type()].size() << endl;
    for (int i_em = 0; i_em < 3; i_em++){
      if ((*CA_collinear)[i_c][0].in_collinear()[i_em] == 1){
	CA_dipole_splitting[(*CA_collinear)[i_c][0].type()][i_em] = 1;
      }
    }
  }

  logger << LOG_DEBUG << "CA_dipole_splitting:" << endl;
  for (int i_em = 0; i_em < 3; i_em++){
    stringstream temp_ss;
    temp_ss << "i_em = " << i_em << ":   ";
    for (int i_dt = 0; i_dt < 4; i_dt++){temp_ss << CA_dipole_splitting[i_dt][i_em] << "   ";}
    logger << LOG_DEBUG << temp_ss.str() << endl;
  }
  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "new collinear dipoles determined " << endl << endl;

  /*
  CA_combination_pdf.resize((*CA_collinear).size(), vector<vector<vector<int> > > (combination_pdf.size()));
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      // CA_combination_pdf[i_c] contains a vector, which contains the usual combination_pdf(3) with 0 -> direction (1, -1), 1 -> parton with x1 (in hadron 1/2 for +1/-1), 2 -> parton with x2 (in hadron 2/1 for +1/-1)
      if (((*CA_collinear)[i_c][0].all_pdf()[i_i][1] == 10) && ((*CA_collinear)[i_c][0].all_pdf()[i_i][2] == 10)){
	for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	  if (i_q == 0){continue;}
	  for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	    if (j_q == 0){continue;}
	    vector<int> new_temp_combination_pdf = (*CA_collinear)[i_c][0].all_pdf()[i_i];
	    new_temp_combination_pdf[1] = i_q;
	    new_temp_combination_pdf[2] = j_q;
	    CA_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	  }
	}
      }
      else if ((*CA_collinear)[i_c][0].all_pdf()[i_i][1] == 10){
	for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	  if (i_q == 0){continue;}
	  vector<int> new_temp_combination_pdf = (*CA_collinear)[i_c][0].all_pdf()[i_i];
	  new_temp_combination_pdf[1] = i_q;
	  CA_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	}
      }
      else if ((*CA_collinear)[i_c][0].all_pdf()[i_i][2] == 10){
	for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	  if (j_q == 0){continue;}
	  vector<int> new_temp_combination_pdf = (*CA_collinear)[i_c][0].all_pdf()[i_i];
	  new_temp_combination_pdf[2] = j_q;
	  CA_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	}
      }
      else {
	CA_combination_pdf[i_c][i_i].push_back((*CA_collinear)[i_c][0].all_pdf()[i_i]);
      }
    }
  }
  logger << LOG_DEBUG << "new collinear-dipole pdfs determined " << endl << endl;
*/

  /* // old position
  CA_Q2f.resize((*CA_collinear).size());
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_x = 0; i_x < CA_combination_pdf[i_c][0].size(); i_x++){
      if ((*CA_collinear)[i_c][0].type_parton()[(*CA_collinear)[i_c][0].no_emitter()] == 22){  // emitter is a photon
	if ((*CA_collinear)[i_c][0].type() == 0){  // hard process with a from a -> a a splitting
	  CA_Q2f[i_c].push_back(1.);
 	}
	else if ((*CA_collinear)[i_c][0].type() == 2){  // hard process with a from q -> q a splitting
	  CA_Q2f[i_c].push_back(charge2_particle[abs(CA_combination_pdf[i_c][0][i_x][(*CA_collinear)[i_c][0].no_emitter()])]);
	}
	else {logger << LOG_WARN << "Wrong splitting: may not happen!" << endl; exit(1);}
      }
      else { // emitter and spectator are charged particles
	CA_Q2f[i_c].push_back(charge2_particle[abs((*CA_collinear)[i_c][0].type_parton()[(*CA_collinear)[i_c][0].no_emitter()])]);
      }
    }
  }
  logger << LOG_DEBUG << "new collinear-dipole Q2f determined " << endl << endl;
  */


  /*
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_x = 0; i_x < (*CA_collinear)[i_c].size(); i_x++){
      cout << "(*CA_collinear)[" << setw(2) << i_c << "][" << setw(2) << i_x << "].charge_factor() = " << (*CA_collinear)[i_c][i_x].charge_factor() << endl;
    }
  }

  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_x = 0; i_x < CA_combination_pdf[i_c][0].size(); i_x++){
      cout << "CA_Q2f[" << setw(2) << i_c << "][" << setw(2) << i_x << "] = " << CA_Q2f[i_c][i_x] << endl;
    }
  }
  */
  /*
  CA_Q2f.resize((*CA_collinear).size());
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    //    CA_Q2f[i_c].resize((*CA_collinear)[i_c].size());
    for (int i_x = 0; i_x < CA_combination_pdf[i_c][0].size(); i_x++){
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	if ((*CA_collinear)[i_c][j_c].no_spectator() == 0){
	  CA_Q2f[i_c][j_c].push_back(charge2_particle[CA_combination_pdf[i_c][0][i_x][(*CA_collinear)[i_c][j_c].type_parton()[(*CA_collinear)[i_c][j_c].no_emitter()]]]);
	}
	else {
	  if ((*CA_collinear)[i_c][j_c].type_parton()[(*CA_collinear)[i_c][j_c].no_emitter()] == 22){
	    if ((*CA_collinear)[i_c][j_c].no_spectator() < 3){
	      CA_Q2f[i_c][j_c].push_back(charge2_particle[(*CA_collinear)[i_c][0].all_pdf()[(*CA_collinear)[i_c][j_c].type_parton()[(*CA_collinear)[i_c][j_c].no_emitter()]]]);
	    }
	    else {
	      CA_Q2f[i_c][j_c].push_back(0.);
	    }
	  }
	  else { // emitter and spectator are charged particles
	    CA_Q2f[i_c][j_c].push_back(charge2_particle[(*CA_collinear)[i_c][j_c].type_parton()[(*CA_collinear)[i_c][j_c].no_emitter()]]);
	  }
	}
      }
    }
  }
  */


  logger << LOG_INFO << "CA_combination_pdf determined " << endl << endl;
  logger.newLine(LOG_INFO);
  output_collinear_pdf();




  logger << LOG_DEBUG << "user.string_value[user.string_map[selection]] = " << user.string_value[user.string_map["selection"]] << endl;
  
  for (int i_c = (*CA_collinear).size() - 1; i_c >=0; i_c--){
    for (int j_c = (*CA_collinear)[i_c].size() - 1; j_c >= 0; j_c--){

      if (user.string_value[user.string_map["selection"]] == "ii"){
	// only initial-initial contribution:
	if ((*CA_collinear)[i_c][j_c].no_emitter() > 2 || (*CA_collinear)[i_c][j_c].no_spectator() > 2){
	  (*CA_collinear)[i_c].erase((*CA_collinear)[i_c].begin() + j_c);
	  logger << LOG_DEBUG << "kill:   (*CA_collinear)[" << i_c << "][" << j_c << "].name() = " << (*CA_collinear)[i_c][j_c].name() << endl;
	}
      }
      else if (user.string_value[user.string_map["selection"]] == "if"){
	// only initial-final contribution:
	if ((*CA_collinear)[i_c][j_c].no_emitter() > 2 || (*CA_collinear)[i_c][j_c].no_spectator() < 3){
	  (*CA_collinear)[i_c].erase((*CA_collinear)[i_c].begin() + j_c);
	  logger << LOG_DEBUG << "kill:   (*CA_collinear)[" << i_c << "][" << j_c << "].name() = " << (*CA_collinear)[i_c][j_c].name() << endl;
	}
      }
      else if (user.string_value[user.string_map["selection"]] == "fi"){
	// only final-initial contribution:
	if ((*CA_collinear)[i_c][j_c].no_emitter() < 3 || (*CA_collinear)[i_c][j_c].no_spectator() > 2){
	  (*CA_collinear)[i_c].erase((*CA_collinear)[i_c].begin() + j_c);
	  logger << LOG_DEBUG << "kill:   (*CA_collinear)[" << i_c << "][" << j_c << "].name() = " << (*CA_collinear)[i_c][j_c].name() << endl;
	}
      }
      else if (user.string_value[user.string_map["selection"]] == "ff"){
	// only final-final contribution:
	if ((*CA_collinear)[i_c][j_c].no_emitter() < 3 || (*CA_collinear)[i_c][j_c].no_spectator() < 3){
	  (*CA_collinear)[i_c].erase((*CA_collinear)[i_c].begin() + j_c);
	  logger << LOG_DEBUG << "kill:   (*CA_collinear)[" << i_c << "][" << j_c << "].name() = " << (*CA_collinear)[i_c][j_c].name() << endl;
	}
      }
    }
  }
  
  logger.newLine(LOG_INFO);



  
  /*
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      logger << LOG_INFO << "pdf contributions: " << setw(15) << left << (*CA_collinear)[i_c][0].all_name()[i_i] << ":   " << endl;
      for (int i_x = 0; i_x < CA_combination_pdf[i_c][i_i].size(); i_x++){
	stringstream temp;
	temp.str("");
	temp << setw(19) << "";
	temp << "CA_Q2f[" << setw(2) << i_c << "][" << setw(2) << i_x << "] = " << setw(10) << setprecision(8) << CA_Q2f[i_c][i_x] << "     "; 
	temp << "CA_combination_pdf[" << setw(2) << i_c << "][" << setw(2) << i_i << "][" << setw(2) << i_x << "] = "; 
	for (int i_y = 0; i_y < CA_combination_pdf[i_c][i_i][i_x].size(); i_y++){temp << setw(2) << right << CA_combination_pdf[i_c][i_i][i_x][i_y] << "  ";}
	logger << LOG_INFO << temp.str() << endl;
      }
    }
    logger.newLine(LOG_INFO);
  }
  */



  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_collinear_QEW(){
  Logger logger("observable_set::calculate_collinear_QEW");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (massive_QEW){calculate_collinear_QEW_CDST();}
  else {calculate_collinear_QEW_CS();}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_collinear_QEW_CS(){
  Logger logger("observable_set::calculate_collinear_QEW_CS");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  static int initialization = 1;
  // calculate all needed momentum-independent (splitting) functions
  static double Kbar[3][4] = {{0.}};
  static double Kt[3][4] = {{0.}};
  static double P[3][4] = {{0.}};
  static double Kbar_plus[3][4] = {{0.}};
  static double Kt_plus[3][4] = {{0.}};
  static double P_plus[3][4] = {{0.}};
  static double intKbar_plus[3][4] = {{0.}};
  static double intKt_plus[3][4] = {{0.}};
  static double intP_plus[3][4] = {{0.}};
  static double Kbar_delta[4] = {0.};
  static double Kt_delta[4] = {0.};
  static double P_delta[4] = {0.};

  static double alpha_e_2pi = msi.alpha_e * inv2pi;
  static vector<vector<int> > pair;
  //  static vector<double> gamma_i((*CA_collinear).size());
  static vector<vector<double> > gamma_i((*CA_collinear).size());
  static vector<vector<double> > ln_papi(3, vector<double> (csi->type_parton[0].size()));
  static vector<vector<vector<vector<double> > > > CA_value_ln_muF_papi(CA_value_log_mu2_fact.size());
  //  static vector<vector<vector<double> > > dataK(3, vector<vector<double> > ((*CA_collinear).size(), vector<double> (1, 0.)));
  static vector<vector<vector<vector<vector<double> > > > > value_dataP(CA_value_log_mu2_fact.size());

  static vector<vector<int> > type((*CA_collinear).size());
  static vector<vector<int> > no_emitter((*CA_collinear).size());
  static vector<vector<int> > no_spectator((*CA_collinear).size());
  static vector<vector<int> > collinear_singularity((*CA_collinear).size());
  static vector<vector<vector<int> > > ppair((*CA_collinear).size());

  static map<int, double> charge_particle;
  static vector<double> charge2_particle(26);

  if (initialization == 1){
    fill_charge_particle(charge_particle);
    for (int i_pdg = 0; i_pdg < 26; i_pdg++){
      charge2_particle[i_pdg] = pow(charge_particle[i_pdg], 2);
    }
    
    for (int i_dt = 0; i_dt < CA_dipole_splitting.size(); i_dt++){
      for (int i_em = 0; i_em < CA_dipole_splitting.size(); i_em++){
	logger << LOG_DEBUG_VERBOSE << "CA_dipole_splitting[" << i_dt << "][" << i_em << "] = " << CA_dipole_splitting[i_dt][i_em] << endl;
      }
    }

    if (CA_dipole_splitting[0][0] == 1 || CA_dipole_splitting[0][1] == 1 || CA_dipole_splitting[0][2] == 1){    // a -> a (+a) splitting (0)
      // !!! massive bottom quarks !!!
      Kbar_delta[0] = Kbar_qew_aa_delta(N_f);  // 16/9*sum[Q²_f] // * (1. + 1. + 1. + N_c * (2 * 4./9. + 2 * 1./9.));
      Kt_delta[0] = Kt_qew_aa_delta();
      P_delta[0] = P_qew_aa_delta(N_f);  // -2/3*sum[Q²_f]
    }
    
    if (CA_dipole_splitting[1][0] == 1 || CA_dipole_splitting[1][1] == 1 || CA_dipole_splitting[1][2] == 1){    // q -> q (+a) splitting
      Kbar_delta[1] = Kbar_qew_qq_delta();
      Kt_delta[1] = Kt_qew_qq_delta();
      P_delta[1] = Pxx_qew_qq_delta();
    }

    for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	if ((*CA_collinear)[i_c][j_c].no_spectator() == 0){continue;}
	int flag = -1;
	for (int i_p = 0; i_p < pair.size(); i_p++){
	  if ((*CA_collinear)[i_c][j_c].pair() == pair[i_p]){flag = i_p; break;}
	}
	if (flag == -1){pair.push_back((*CA_collinear)[i_c][j_c].pair());}
      }
    }
 
    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      CA_value_ln_muF_papi[sd].resize(CA_value_log_mu2_fact[sd].size(), vector<vector<double> > (3, vector<double> (csi->type_parton[0].size())));
    }
    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      value_dataP[sd].resize(CA_value_log_mu2_fact[sd].size(), vector<vector<vector<double> > > (3, vector<vector<double> > ((*CA_collinear).size())));
      for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	for (int i_x = 0; i_x < 3; i_x++){
	  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
	    value_dataP[sd][ss][i_x][i_c].resize((*CA_collinear)[i_c].size(), 0.);
	  }
	}
      }
    }

    for (int i_c = 0; i_c < gamma_i.size(); i_c++){ // i_c = -1 is not needed !!!
      gamma_i[i_c].resize((*CA_collinear)[i_c].size());
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	if ((*CA_collinear)[i_c][j_c].type_parton()[(*CA_collinear)[i_c][j_c].no_spectator()] == 22){cout << "Does this ever happen ???" << endl; gamma_i[i_c][j_c] = gamma_qew_a(N_f);}
	else {gamma_i[i_c][j_c] = gamma_qew_q;}
      }
    }

    for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
      type[i_c].resize((*CA_collinear)[i_c].size());
      no_emitter[i_c].resize((*CA_collinear)[i_c].size());
      no_spectator[i_c].resize((*CA_collinear)[i_c].size());
      collinear_singularity[i_c].resize((*CA_collinear)[i_c].size());
      ppair[i_c].resize((*CA_collinear)[i_c].size());
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	type[i_c][j_c] = (*CA_collinear)[i_c][j_c].type();
	no_emitter[i_c][j_c] = (*CA_collinear)[i_c][j_c].no_emitter();
	no_spectator[i_c][j_c] = (*CA_collinear)[i_c][j_c].no_spectator();
	collinear_singularity[i_c][j_c] = (*CA_collinear)[i_c][j_c].in_collinear()[0];
	ppair[i_c][j_c] = (*CA_collinear)[i_c][j_c].pair();
      }
    }

    for (int i_p = 0; i_p < pair.size(); i_p++){logger << LOG_DEBUG << "pair[" << i_p << "] = " << "(" << pair[i_p][0] << ", " << pair[i_p][1] << ")" << endl;}
    initialization = 0;
    logger << LOG_DEBUG_VERBOSE << "initialization finished!" << endl;
  }


  for (int i_c = 1; i_c < 3; i_c++){
    
    if (CA_dipole_splitting[0][i_c] == 1){      // a -> a (+a) splitting
      Kbar[i_c][0] = Kbar_qew_aa(z_coll[i_c]);  // 0
      Kt[i_c][0] = Kt_qew_aa(z_coll[i_c]);  // 0
      P[i_c][0] = P_qew_aa(z_coll[i_c]);  // 0
      if (CA_dipole_splitting[0][0] == 1){
	Kbar_plus[i_c][0] = Kbar_qew_aa_plus(z_coll[i_c]);  // 0
	Kt_plus[i_c][0] = Kt_qew_aa_plus(z_coll[i_c]);  // 0
	P_plus[i_c][0] = P_qew_aa_plus(z_coll[i_c]);  // 0
	intKbar_plus[i_c][0] = intKbar_qew_aa_plus(x_pdf[i_c]);  // 0
	intKt_plus[i_c][0] = intKt_qew_aa_plus(x_pdf[i_c]);  // 0
	intP_plus[i_c][0] = intP_qew_aa_plus(x_pdf[i_c]);  // 0
      }
      else {logger << LOG_ERROR << "May not happen !!! a -> a splitting without irregular terms !!!" << endl;}
    }
    
    else if (CA_dipole_splitting[1][i_c] == 1){      // q -> q (+a) splitting
      Kbar[i_c][1] = Kbar_qew_qq(z_coll[i_c]);
      Kt[i_c][1] = Kt_qew_qq(z_coll[i_c]);
      P[i_c][1] = Pxx_qew_qq(z_coll[i_c]);
      if (CA_dipole_splitting[1][0] == 1){
	Kbar_plus[i_c][1] = Kbar_qew_qq_plus(z_coll[i_c]);
	Kt_plus[i_c][1] = Kt_qew_qq_plus(z_coll[i_c]);
	P_plus[i_c][1] = Pxx_qew_qq_plus(z_coll[i_c]);
	intKbar_plus[i_c][1] = intKbar_qew_qq_plus(x_pdf[i_c]);
	intKt_plus[i_c][1] = intKt_qew_qq_plus(x_pdf[i_c]);
	intP_plus[i_c][1] = intPxx_qew_qq_plus(x_pdf[i_c]);
      }
      else {logger << LOG_ERROR << "May not happen !!! q -> q splitting without irregular terms !!!" << endl;}
    }
    if (CA_dipole_splitting[2][i_c] == 1){      // q -> a (+q) splitting ??? not checked yet !!!
      Kbar[i_c][2] = Kbar_qew_qa(z_coll[i_c]);
      Kt[i_c][2] = Kt_qew_qa(z_coll[i_c]);
      P[i_c][2] = P_qew_qa(z_coll[i_c]);
    }
    if (CA_dipole_splitting[3][i_c] == 1){      // a -> q (+q~) splitting ??? not checked yet !!!
      Kbar[i_c][3] = Kbar_qew_aq(z_coll[i_c]);
      Kt[i_c][3] = Kt_qew_aq(z_coll[i_c]);
      P[i_c][3] = P_qew_aq(z_coll[i_c]);
      logger << LOG_DEBUG_VERBOSE << "Kbar[" << i_c << "][3] = " << Kbar[i_c][3] << endl;
      logger << LOG_DEBUG_VERBOSE << "Kt[" << i_c << "][3] = " << Kt[i_c][3] << endl;
      logger << LOG_DEBUG_VERBOSE << "P[" << i_c << "][3] = " << P[i_c][3] << endl;

    }
  }

  for (int i_p = 0; i_p < pair.size(); i_p++){
    ln_papi[pair[i_p][0]][pair[i_p][1]] = log(2 * p_parton[0][pair[i_p][0]] * p_parton[0][pair[i_p][1]]);
    //    logger << LOG_DEBUG_VERBOSE << "before   value_logscale2_fact_papi" << endl;
    if (switch_TSV){
      for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	for (int i_m = 0; i_m < n_scale_dyn_fact[v_sf]; i_m++){
	  value_logscale2_fact_papi[v_sf][i_m][pair[i_p][0]][pair[i_p][1]] = value_central_logscale2_fact[v_sf] + value_relative_logscale2_fact[v_sf][i_m] - ln_papi[pair[i_p][0]][pair[i_p][1]];
	  /*
	  //	  logger << LOG_DEBUG_VERBOSE << "value_scale_fact         [" << v_sf << "][" << i_m << "] = " << value_scale_fact[v_sf][i_m] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_central_logscale2_fact  [" << v_sf << "]    = " << value_central_logscale2_fact[v_sf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_relative_logscale2_fact [" << v_sf << "][" << i_m << "] = " << value_relative_logscale2_fact[v_sf][i_m] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_logscale2_fact          [" << v_sf << "][" << i_m << "] = " << value_central_logscale2_fact[v_sf] + value_relative_logscale2_fact[v_sf][i_m] << endl;
	  logger << LOG_DEBUG_VERBOSE << "pa = " << p_parton[0][pair[i_p][0]] << endl;
	  logger << LOG_DEBUG_VERBOSE << "pb = " << p_parton[0][pair[i_p][1]] << endl;
	  logger << LOG_DEBUG_VERBOSE << "2papb = " << 2 * p_parton[0][pair[i_p][0]] * p_parton[0][pair[i_p][1]] << endl;
	  logger << LOG_DEBUG_VERBOSE << "log(2papb) = " << ln_papi[pair[i_p][0]][pair[i_p][1]] << endl;
	  logger << LOG_DEBUG_VERBOSE << "log(muF²/2papb) = " << value_logscale2_fact_papi[v_sf][i_m][pair[i_p][0]][pair[i_p][1]] << endl;
	  */
	}
      }
    }
    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	CA_value_ln_muF_papi[sd][ss][pair[i_p][0]][pair[i_p][1]] = CA_value_log_mu2_fact[sd][ss] - ln_papi[pair[i_p][0]][pair[i_p][1]];
      }
    }
  }


  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
      //    int i_ca = 0;
      logger << LOG_DEBUG_VERBOSE << "(*CA_collinear)[" << i_c << "][" << j_c << "] = " << (*CA_collinear)[i_c][j_c].name() << endl;
      // K terms 
      if (no_spectator[i_c][j_c] == 0){ // for (anti-)quarks, factor Q²_f included in the matrix element CA_ME2_cf[i_c][j_c] (which is Born * Q²_f) ! ???
	data_K[0][i_c][j_c] = (Kbar[no_emitter[i_c][j_c]][type[i_c][j_c]] + Kbar_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
	data_K[1][i_c][j_c] = (-Kbar_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];
	data_K[2][i_c][j_c] = (Kbar_delta[type[i_c][j_c]] - intKbar_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
	logger << LOG_DEBUG_VERBOSE << "no_emitter[" << i_c << "][" << j_c << "] = " << no_emitter[i_c][j_c] << "   no_spectator[" << i_c << "][" << j_c << "] = " << no_spectator[i_c][j_c] << "   CS data_K[0][" << i_c << "][" << j_c << "] = " << data_K[0][i_c][j_c] << endl;
	logger << LOG_DEBUG_VERBOSE << "no_emitter[" << i_c << "][" << j_c << "] = " << no_emitter[i_c][j_c] << "   no_spectator[" << i_c << "][" << j_c << "] = " << no_spectator[i_c][j_c] << "   CS data_K[1][" << i_c << "][" << j_c << "] = " << data_K[1][i_c][j_c] << endl;
	logger << LOG_DEBUG_VERBOSE << "no_emitter[" << i_c << "][" << j_c << "] = " << no_emitter[i_c][j_c] << "   no_spectator[" << i_c << "][" << j_c << "] = " << no_spectator[i_c][j_c] << "   CS data_K[2][" << i_c << "][" << j_c << "] = " << data_K[2][i_c][j_c] << endl;
	logger << LOG_DEBUG_VERBOSE << "Kbar_delta[type[i_c][j_c] = " << type[i_c][j_c] << "] = " << Kbar_delta[type[i_c][j_c]] << endl;
	logger << LOG_DEBUG_VERBOSE << "intKbar_plus[no_emitter[i_c][j_c] = " << no_emitter[i_c][j_c] << "][type[i_c][j_c] = " << type[i_c][j_c] << "] = " << intKbar_plus[no_emitter[i_c][j_c]][type[i_c][j_c]] << endl;
      }
      else if (no_spectator[i_c][j_c] < 3){ // for (anti-)quarks, factor Q²_f always drops out due to T²_a->Q²_f in the denominators ???
	data_K[0][i_c][j_c] = -(Kt[no_emitter[i_c][j_c]][type[i_c][j_c]] + Kt_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
	data_K[1][i_c][j_c] = -(-Kt_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];
	data_K[2][i_c][j_c] = -(Kt_delta[type[i_c][j_c]] - intKt_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
	logger << LOG_DEBUG_VERBOSE << "no_emitter[" << i_c << "][" << j_c << "] = " << no_emitter[i_c][j_c] << "   no_spectator[" << i_c << "][" << j_c << "] = " << no_spectator[i_c][j_c] << "   CS data_K[0][" << i_c << "][" << j_c << "] = " << data_K[0][i_c][j_c] << endl;
	logger << LOG_DEBUG_VERBOSE << "Kt     [no_emitter[" << i_c << "][" << j_c << "] = " << no_emitter[i_c][j_c] << "][type[" << i_c << "][" << j_c << "] = " << type[i_c][j_c] << "] = " << Kt[no_emitter[i_c][j_c]][type[i_c][j_c]] << endl;
	logger << LOG_DEBUG_VERBOSE << "Kt_plus[no_emitter[" << i_c << "][" << j_c << "] = " << no_emitter[i_c][j_c] << "][type[" << i_c << "][" << j_c << "] = " << type[i_c][j_c] << "] = " << Kt_plus[no_emitter[i_c][j_c]][type[i_c][j_c]] << endl;
	logger << LOG_DEBUG_VERBOSE << "CA_ME2_cf[" << i_c << "][" << j_c << "] = " << CA_ME2_cf[i_c][j_c] << endl;
      }
      else if (no_spectator[i_c][j_c] > 2){ // for (anti-)quarks, factor Q²_f always drops out due to T²_a->Q²_f in the denominators ???
	if (collinear_singularity[i_c][j_c] == 1){
	  data_K[0][i_c][j_c] = gamma_i[i_c][j_c] * (1. / (1. - z_coll[no_emitter[i_c][j_c]])) * CA_ME2_cf[i_c][j_c];
	  data_K[1][i_c][j_c] = -gamma_i[i_c][j_c] * (1. / (1. - z_coll[no_emitter[i_c][j_c]])) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];
	  data_K[2][i_c][j_c] = gamma_i[i_c][j_c] * (1. + log(1. - x_pdf[no_emitter[i_c][j_c]])) * CA_ME2_cf[i_c][j_c];
 	}
	//      logger << LOG_DEBUG_VERBOSE << "no_emitter[" << i_c << "][" << j_c << "] = " << no_emitter[i_c][j_c] << "   no_spectator[" << i_c << "][" << j_c << "] = " << no_spectator[i_c][j_c] << "   CS data_K[0][" << i_c << "][" << j_c << "] = " << data_K[0][i_c][j_c] << endl;
      }
      
      //      logger << LOG_DEBUG_VERBOSE << "no_emitter[" << i_c << "][" << j_c << "] = " << no_emitter[i_c][j_c] << "   no_spectator[" << i_c << "][" << j_c << "] = " << no_spectator[i_c][j_c] << "   CS data_K[0][" << i_c << "][" << j_c << "] = " << data_K[0][i_c][j_c] << endl;
      //      logger << LOG_DEBUG_VERBOSE << "CS data_K[0][" << i_c << "][" << j_c << "] = " << data_K[0][i_c][j_c] << endl;
      //      logger << LOG_DEBUG_VERBOSE << "CS data_K[1][" << i_c << "][" << j_c << "] = " << data_K[1][i_c][j_c] << endl;
      //      logger << LOG_DEBUG_VERBOSE << "CS data_K[2][" << i_c << "][" << j_c << "] = " << data_K[2][i_c][j_c] << endl;
      
      if (no_spectator[i_c][j_c] != 0){ // for (anti-)quarks, factor Q²_f always drops out due to T²_a->Q²_f in the denominators
	// P terms 
	logger << LOG_DEBUG_VERBOSE << "before value_logscale2_fact_papi" << endl;
	if (switch_TSV){
	  for (int v_sf = 0; v_sf < value_logscale2_fact_papi.size(); v_sf++){
	    for (int v_xf = 0; v_xf < value_logscale2_fact_papi[v_sf].size(); v_xf++){
	      value_data_P[v_sf][v_xf][0][i_c][j_c] = (P[no_emitter[i_c][j_c]][type[i_c][j_c]] + P_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	      value_data_P[v_sf][v_xf][1][i_c][j_c] = (-P_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	      value_data_P[v_sf][v_xf][2][i_c][j_c] = (P_delta[type[i_c][j_c]] - intP_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	      /*
	      logger << LOG_DEBUG_VERBOSE << "log(muF²/2papb) = " << value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] << endl;
	      logger << LOG_DEBUG_VERBOSE << "P_delta[type[i_c][j_c]] = " << P_delta[type[i_c][j_c]] << endl;
	      logger << LOG_DEBUG_VERBOSE << "intP_plus[no_emitter[i_c][j_c]][type[i_c][j_c]] = " << intP_plus[no_emitter[i_c][j_c]][type[i_c][j_c]] << endl;
	      logger << LOG_DEBUG_VERBOSE << "CA_ME2_cf[i_c][j_c] = " << CA_ME2_cf[i_c][j_c] << endl;
	      logger << LOG_DEBUG_VERBOSE << "P_aa_delta * ME2_born * log(muF²/2papb) = " << value_data_P[v_sf][v_xf][2][i_c][j_c] << endl;
	      */
	    }
	  }
	}
	for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
	  for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	    value_dataP[sd][ss][0][i_c][j_c] = (P[no_emitter[i_c][j_c]][type[i_c][j_c]] + P_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	    value_dataP[sd][ss][1][i_c][j_c] = (-P_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	    value_dataP[sd][ss][2][i_c][j_c] = (P_delta[type[i_c][j_c]] - intP_plus[no_emitter[i_c][j_c]][type[i_c][j_c]]) * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	  }
	}
      }
    }
    logger << LOG_DEBUG_VERBOSE << "before value_ME2_KP calculation" << endl;
    for (int i_x = 0; i_x < 3; i_x++){
      double temp_sumK = accumulate(data_K[i_x][i_c].begin(), data_K[i_x][i_c].end(), 0.);
      if (switch_TSV){
	for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	  for (int v_xf = 0; v_xf < value_logscale2_fact_papi[v_sf].size(); v_xf++){
	    if (switch_KP == 0){
	      // K + P terms
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_e_2pi * (temp_sumK + accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_e_2pi * (temp_sumK + accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	    }
	    else if (switch_KP == 1){
	      // P terms 
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_e_2pi * (accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_e_2pi * (accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	    }
	    else if (switch_KP == 2){
	      // K terms 
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_e_2pi * (temp_sumK);
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_e_2pi * (temp_sumK);
	    }
	  }
	}
	logger << LOG_DEBUG_VERBOSE << "after value_ME2_KP calculation   i_x = " << i_x << endl;
      }

      //    }
      //    dataK = data_K;
      //    for (int i_x = 0; i_x < 3; i_x++){
      for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
	for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	  if (switch_KP == 0){
	    // K + P term
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_e_2pi * (temp_sumK + accumulate(value_dataP[sd][ss][i_x][i_c].begin(), value_dataP[sd][ss][i_x][i_c].end(), 0.));
	  }
	  else if (switch_KP == 1){
	    // P terms 
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_e_2pi * (accumulate(value_dataP[sd][ss][i_x][i_c].begin(), value_dataP[sd][ss][i_x][i_c].end(), 0.));
	  }
	  else if (switch_KP == 2){
	    // K terms 
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_e_2pi * (temp_sumK);
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_collinear_QEW_CDST(){
  static Logger logger("observable_set::calculate_collinear_QEW_CDST");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  static int initialization = 1;
  // calculate all needed momentum-independent (splitting) functions
  // check if all elements are always zero!!!
  /*
  double Kbar[4][3] = {{0.}};
  double Kt[4][3] = {{0.}};
  double P[4][3] = {{0.}};
  double Kbar_plus[4][3] = {{0.}};
  double Kt_plus[4][3] = {{0.}};
  double P_plus[4][3] = {{0.}};
  double intKbar_plus[4][3] = {{0.}};
  double intKt_plus[4][3] = {{0.}};
  double intP_plus[4][3] = {{0.}};
  double Kbar_delta[4] = {0.};
  double Kt_delta[4] = {0.};
  double P_delta[4] = {0.};
  */
  // check if everything remains unchanged without setting all functions to 0 !!!
  // reconsider if different partonic channels are calculated simultaneously later !!!
  static double Kbar[4][3] = {{0.}};
  static double Kt[4][3] = {{0.}};
  static double P[4][3] = {{0.}};
  static double Kbar_plus[4][3] = {{0.}};
  static double Kt_plus[4][3] = {{0.}};
  static double P_plus[4][3] = {{0.}};
  static double intKbar_plus[4][3] = {{0.}};
  static double intKt_plus[4][3] = {{0.}};
  static double intP_plus[4][3] = {{0.}};
  static double Kbar_delta[4] = {0.};
  static double Kt_delta[4] = {0.};
  static double P_delta[4] = {0.};

  static double alpha_e_2pi = msi.alpha_e * inv2pi;

  static vector<vector<int> > pair;
  /*
  static vector<double> iT2_ap(3);
  static vector<double> gamma_a_T2_ap(3);
  static vector<vector<double> > gamma_i_T2_i((*CA_collinear).size());
  */
  static vector<double> gamma_ax(3);
  static vector<vector<double> > gamma_i((*CA_collinear).size());

  static vector<vector<double> > ln_papi(3, vector<double> (csi->type_parton[0].size()));
  static vector<vector<vector<vector<double> > > > CA_value_ln_muF_papi(CA_value_log_mu2_fact.size());
  static vector<vector<vector<vector<vector<double> > > > > value_dataP(CA_value_log_mu2_fact.size());

  static vector<vector<int> > type((*CA_collinear).size());
  static vector<vector<int> > no_emitter((*CA_collinear).size());
  static vector<vector<int> > no_spectator((*CA_collinear).size());
  static vector<vector<int> > collinear_singularity((*CA_collinear).size());
  static vector<vector<vector<int> > > ppair((*CA_collinear).size());

  static int n_max_spectator = 0;
  if (initialization == 1){
    for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
      if ((*CA_collinear)[i_c].size() > n_max_spectator){n_max_spectator = (*CA_collinear)[i_c].size();}
    }
  }
  static vector<vector<vector<double> > > Kit(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > intKit_plus(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > Kit_plus_x(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > Kit_plus_1(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > Kit_plus_outside_x(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > Kit_plus_outside_1(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<vector<vector<double> > > Kit_delta(4, vector<vector<double> > (3, vector<double> (n_max_spectator, 0.)));
  static vector<double> m_Q(csi->type_parton[0].size(), 0.);
  static vector<double> m2_Q(csi->type_parton[0].size(), 0.);
  static vector<vector<double> > sall_ja(3, vector<double> (csi->type_parton[0].size(), 0.));
  static vector<vector<double> > sall_ja_x(3, vector<double> (csi->type_parton[0].size(), 0.));
  static vector<vector<double> > mu2_Q(3, vector<double> (csi->type_parton[0].size(), 0.));
  static vector<vector<double> > mu2_Q_x(3, vector<double> (csi->type_parton[0].size(), 0.));

  if (initialization == 1){
    for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){
      if (i_p < 3 && M2[abs(csi->type_parton[0][i_p])] != 0.){cout << "Incoming massive partons are not supported!" << endl; int_end = 1;}
      //      if (i_p < 3 && M2[abs(csi->type_parton[0][i_p])] != 0.){cout << "Incoming massive partons are not supported!" << endl; exit(1);}
      if (M2[abs(csi->type_parton[0][i_p])] != 0.){
	//	m_Q[i_p] = M[abs(csi->type_parton[0][i_p])];
	//	m2_Q[i_p] = M2[abs(csi->type_parton[0][i_p])];
	m_Q[i_p] = mass_parton[0][i_p];
	m2_Q[i_p] = mass2_parton[0][i_p];
      }
      //      cout << "m2_Q[" << i_p << "] = " << m2_Q[i_p] << endl;
    }
    /*
    if (CA_dipole_splitting[0][1] == 1 || CA_dipole_splitting[0][2] == 1){    // g -> g (+g) splitting (0)
      Kbar_delta[0] = Kbar_gg_delta(N_f);
      Kt_delta[0] = Kt_gg_delta();
      P_delta[0] = P_gg_delta(N_f);
    }
    */
    if (CA_dipole_splitting[1][1] == 1 || CA_dipole_splitting[1][2] == 1){    // q -> q (+g) splitting
      Kbar_delta[1] = Kbar_qq_delta() / C_F;
      Kt_delta[1] = Kt_qq_delta() / C_F;
      P_delta[1] = Pxx_qq_delta() / C_F;
    }
  
    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      CA_value_ln_muF_papi[sd].resize(CA_value_log_mu2_fact[sd].size(), vector<vector<double> > (3, vector<double> (csi->type_parton[0].size())));
    }
    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      value_dataP[sd].resize(CA_value_log_mu2_fact[sd].size(), vector<vector<vector<double> > > (3, vector<vector<double> > ((*CA_collinear).size())));
      for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	for (int i_x = 0; i_x < 3; i_x++){
	  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
	    value_dataP[sd][ss][i_x][i_c].resize((*CA_collinear)[i_c].size(), 0.);
	  }
	}
      }
    }

    for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	if ((*CA_collinear)[i_c][j_c].no_spectator() == 0){continue;}
	int flag = 0;
	for (int i_p = 0; i_p < pair.size(); i_p++){
	  if ((*CA_collinear)[i_c][j_c].pair() == pair[i_p]){flag = 1; break;}
	}
	if (flag == 0){pair.push_back((*CA_collinear)[i_c][j_c].pair());}
      }
    }
    /*
    for (int i_em = 1; i_em < 3; i_em++){
      if (csi->type_parton[0][i_em] == 0){iT2_ap[i_em] = 1. / C_A;}
      else {iT2_ap[i_em] = 1. / C_F;}
    }
    */
    /*
    for (int i_em = 1; i_em < 3; i_em++){
      if (csi->type_parton[0][i_em] == 0){gamma_a_T2_ap[i_em] = gamma_g(N_f) / C_A;}
      else {gamma_a_T2_ap[i_em] = gamma_q / C_F;}
    }
    for (int i_c = 0; i_c < gamma_i_T2_i.size(); i_c++){ // i_c = -1 is not needed !!!
      gamma_i_T2_i[i_c].resize((*CA_collinear)[i_c].size());
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	if ((*CA_collinear)[i_c][j_c].type_parton()[(*CA_collinear)[i_c][j_c].no_spectator()] == 0){gamma_i_T2_i[i_c][j_c] = gamma_g(N_f) / C_A;}
	else {gamma_i_T2_i[i_c][j_c] = gamma_q / C_F;}
      }
    }
    */
    for (int i_em = 1; i_em < 3; i_em++){
      if (csi->type_parton[0][i_em] == 22){gamma_ax[i_em] = gamma_a(N_f);}
      else {gamma_ax[i_em] = gamma_qew_q;}
    }
    for (int i_c = 0; i_c < gamma_i.size(); i_c++){ // i_c = -1 is not needed !!!
      gamma_i[i_c].resize((*CA_collinear)[i_c].size());
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	if ((*CA_collinear)[i_c][j_c].type_parton()[(*CA_collinear)[i_c][j_c].no_spectator()] == 22){gamma_i[i_c][j_c] = gamma_a(N_f);}
	else {gamma_i[i_c][j_c] = gamma_qew_q;}
      }
    }

    for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
      type[i_c].resize((*CA_collinear)[i_c].size());
      no_emitter[i_c].resize((*CA_collinear)[i_c].size());
      no_spectator[i_c].resize((*CA_collinear)[i_c].size());
      collinear_singularity[i_c].resize((*CA_collinear)[i_c].size());
      ppair[i_c].resize((*CA_collinear)[i_c].size());
      for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	type[i_c][j_c] = (*CA_collinear)[i_c][j_c].type();
	no_emitter[i_c][j_c] = (*CA_collinear)[i_c][j_c].no_emitter();
	no_spectator[i_c][j_c] = (*CA_collinear)[i_c][j_c].no_spectator();
	collinear_singularity[i_c][j_c] = (*CA_collinear)[i_c][j_c].in_collinear()[0];
	ppair[i_c][j_c] = (*CA_collinear)[i_c][j_c].pair();
      }
    }


    for (int i_p = 0; i_p < pair.size(); i_p++){
      int pair_em = pair[i_p][0];
      int pair_sp = pair[i_p][1];
      
      if (m2_Q[pair_sp] > 0.){
	for (int i_dt = 0; i_dt < 4; i_dt++){
	  Kit[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_x[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_1[i_dt][pair_em][pair_sp] = 0.;
	  intKit_plus[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_outside_x[i_dt][pair_em][pair_sp] = 0.;
	  Kit_plus_outside_1[i_dt][pair_em][pair_sp] = 0.; 
	  Kit_delta[i_dt][pair_em][pair_sp] = 0.; 
	}
      }
    }



    for (int i_p = 0; i_p < pair.size(); i_p++){cout << "pair[" << i_p << "] = " << "(" << pair[i_p][0] << ", " << pair[i_p][1] << ")" << endl;}
    initialization = 0;
  }
  logger << LOG_DEBUG_VERBOSE << "initialization finished!" << endl;
  for (int i_em = 1; i_em < 3; i_em++){
  /*
    if (CA_dipole_splitting[0][i_em] == 1){      // a -> a (+a) splitting
      Kbar[0][i_em] = Kbar_gg(z_coll[i_em]);
      Kt[0][i_em] = Kt_gg(z_coll[i_em]);
      P[0][i_em] = P_gg(z_coll[i_em]);
      if (CA_dipole_splitting[0][0] == 1){
	Kbar_plus[0][i_em] = Kbar_gg_plus(z_coll[i_em]);
	Kt_plus[0][i_em] = Kt_gg_plus(z_coll[i_em]);
	P_plus[0][i_em] = P_gg_plus(z_coll[i_em]);
	intKbar_plus[0][i_em] = intKbar_gg_plus(x_pdf[i_em]);
	intKt_plus[0][i_em] = intKt_gg_plus(x_pdf[i_em]);
	intP_plus[0][i_em] = intP_gg_plus(x_pdf[i_em]);
      }
      else {cout << "May not happen !!! g -> g splitting without irregular terms !!!" << endl;}
    }
  */
    if (CA_dipole_splitting[1][i_em] == 1){      // q -> q (+a) splitting
      Kbar[1][i_em] = Kbar_qq(z_coll[i_em]) / C_F;
      Kt[1][i_em] = Kt_qq(z_coll[i_em]) / C_F;
      P[1][i_em] = Pxx_qq(z_coll[i_em]) / C_F;
      if (CA_dipole_splitting[1][0] == 1){
	Kbar_plus[1][i_em] = Kbar_qq_plus(z_coll[i_em]) / C_F;
	Kt_plus[1][i_em] = Kt_qq_plus(z_coll[i_em]) / C_F;
	P_plus[1][i_em] = Pxx_qq_plus(z_coll[i_em]) / C_F;
	intKbar_plus[1][i_em] = intKbar_qq_plus(x_pdf[i_em]) / C_F;
	intKt_plus[1][i_em] = intKt_qq_plus(x_pdf[i_em]) / C_F;
	intP_plus[1][i_em] = intPxx_qq_plus(x_pdf[i_em]) / C_F;
      }
      else {cout << "May not happen !!! q -> q splitting without irregular terms !!!" << endl;}
    }
    if (CA_dipole_splitting[2][i_em] == 1){      // q -> a (+q) splitting
      Kbar[2][i_em] = Kbar_qg(z_coll[i_em]) / C_F;
      Kt[2][i_em] = Kt_qg(z_coll[i_em]) / C_F;
      P[2][i_em] = P_qg(z_coll[i_em]) / C_F;
    }
    if (CA_dipole_splitting[3][i_em] == 1){      // g -> q (+q~) splitting
      Kbar[3][i_em] = Kbar_gq(z_coll[i_em]) / T_R * N_c;
      Kt[3][i_em] = Kt_gq(z_coll[i_em]) / T_R * N_c;
      P[3][i_em] = P_gq(z_coll[i_em]) / T_R * N_c;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "splitting kernels finished!" << endl;

  for (int i_p = 0; i_p < pair.size(); i_p++){
    int pair_em = pair[i_p][0];
    int pair_sp = pair[i_p][1];
    // exception: 1 -- 2; however, in this case, for all relevant configurations the involved functions are symmetric!
    sall_ja[pair_em][pair_sp] = 2 * p_parton[0][pair_em] * p_parton[0][pair_sp];
    ln_papi[pair_em][pair_sp] = log(sall_ja[pair_em][pair_sp]);
    //    cout << "m2_Q[" << pair_sp << "] = " << m2_Q[pair_sp] << endl;

    if (m2_Q[pair_sp] > 0.){
      // Only happens if no initial-initial dipole is discussed, i.e. pair_em and pair_sp really point at emitter and spectator, respectively.
      sall_ja_x[pair_em][pair_sp] = sall_ja[pair_em][pair_sp] / z_coll[pair_em];
      mu2_Q[pair_em][pair_sp] = m2_Q[pair_sp] / sall_ja[pair_em][pair_sp];
      mu2_Q_x[pair_em][pair_sp] = m2_Q[pair_sp] / sall_ja_x[pair_em][pair_sp];
      //      cout << "i_p = " << i_p << "   pair_em = " << pair_em << "   pair_sp = " << pair_sp << "   z_coll[" << pair_em << "] = " << z_coll[pair_em] << endl;
      // only for massive quarks as spectators
      for (int i_dt = 0; i_dt < 4; i_dt++){
	/*
	if      (i_dt == 0 && (CA_dipole_splitting[i_dt][1] == 1 || CA_dipole_splitting[i_dt][2] == 1)){
	  // g -> g (+g) splitting
	  //	  cout << "g -> g (+g): dipole_phasespace[0][" << pair_sp + 6 << "] = " << dipole_phasespace[0][pair_sp + 6] << "   " << m_Q[pair_sp] << endl;
	  Kit[i_dt][pair_em][pair_sp] = 
	    - 2 * log(2. - z_coll[pair_em]) / (1. - z_coll[pair_em]) // from second term in (6.58) (-> K^qq_q) [included from (6.60)]
	    + 2 * m2_Q[pair_sp] / (z_coll[pair_em] * sall_ja_x[pair_em][pair_sp]) * log(m2_Q[pair_sp] / ((1. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp] + m2_Q[pair_sp])); // from (6.59) (-> K^qg_q) [included from (6.60)]
	  Kit_plus_x[i_dt][pair_em][pair_sp] = 
	    + 2 * log(1. - z_coll[pair_em]) / (1. - z_coll[pair_em]) // from first term in (6.58) 
	    + (1. - z_coll[pair_em]) / (2 * pow(1. - z_coll[pair_em] + mu2_Q_x[pair_em][pair_sp], 2)) // from first term in (5.58, J_gQ^a) included from (6.58) 
	    - 2 / (1. - z_coll[pair_em]) * (1. + log(1. - z_coll[pair_em] + mu2_Q_x[pair_em][pair_sp])) // from second term in (5.58, J_gQ^a) included from (6.58)
	    ;
	  Kit_plus_1[i_dt][pair_em][pair_sp] = 
	    + 2 * log(1. - z_coll[pair_em]) / (1. - z_coll[pair_em]) // from first term in (6.58) 
	    + (1. - z_coll[pair_em]) / (2 * pow(1. - z_coll[pair_em] + mu2_Q[pair_em][pair_sp], 2)) // from first term in (5.58, J_gQ^a) included from (6.58) 
	    - 2 / (1. - z_coll[pair_em]) * (1. + log(1. - z_coll[pair_em] + mu2_Q[pair_em][pair_sp])) // from second term in (5.58, J_gQ^a) included from (6.58)
	    ;
	  intKit_plus[i_dt][pair_em][pair_sp] = 
	    - pow(log(1. - x_pdf[pair_em]), 2) // from first term in (6.58)
	    + .5 * (- mu2_Q[pair_em][pair_sp] / (1. - x_pdf[pair_em] + mu2_Q[pair_em][pair_sp]) + mu2_Q[pair_em][pair_sp] / (1. + mu2_Q[pair_em][pair_sp]) - log((1. - x_pdf[pair_em] + mu2_Q[pair_em][pair_sp]) / (1. + mu2_Q[pair_em][pair_sp]))) // from first term in (5.58) included from (6.58)
	    + 2 * (gsl_sf_dilog(-1. / mu2_Q[pair_em][pair_sp]) - gsl_sf_dilog(-(1. - x_pdf[pair_em]) / mu2_Q[pair_em][pair_sp]) + log(1. - x_pdf[pair_em]) * (1. + log(mu2_Q[pair_em][pair_sp]))) // from second term in (5.58) included from (6.58)
	    ;
	  
	  // terms containing x-dependent pre-factor of (2/(1-z_coll[pair_em]))_+
	  Kit_plus_outside_x[i_dt][pair_em][pair_sp] = 
	    + log(((2. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp]) / ((2. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp] + m2_Q[pair_sp])) // from fourth term from (6.58)
	    + log(2. + mu2_Q_x[pair_em][pair_sp] - z_coll[pair_em]) // from third term in (5.58) included from (6.58)
	    ;
	  Kit_plus_outside_1[i_dt][pair_em][pair_sp] = 
	    + log(sall_ja[pair_em][pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp])) // from fourth term from (6.58)
	    + log(1. + mu2_Q[pair_em][pair_sp]) // from third term in (5.58) included from (6.58)
	    ;
	  Kit_delta[i_dt][pair_em][pair_sp] = -gamma_q / C_F // from fifth term from (6.58)
	    + mu2_Q[pair_em][pair_sp] * log(m2_Q[pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp])) // from sixth term from (6.58)
	    + .5 * m2_Q[pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp])// from seventh term from (6.58)
	    ;
	}
	*/
	if (i_dt == 1 && (CA_dipole_splitting[i_dt][1] == 1 || CA_dipole_splitting[i_dt][2] == 1)){
	  // q -> q (+g) splitting
	  //	  cout << "q -> q (+g): dipole_phasespace[0][" << pair_sp + 6 << "] = " << dipole_phasespace[0][pair_sp + 6] << "   " << m_Q[pair_sp] << endl;
	  Kit[i_dt][pair_em][pair_sp] = -2 * log(2. - z_coll[pair_em]) / (1. - z_coll[pair_em]);
	  Kit_plus_x[i_dt][pair_em][pair_sp] = 0.
	    + 2 * log(1. - z_coll[pair_em]) / (1. - z_coll[pair_em]) // from first term in (6.58) 
	    + (1. - z_coll[pair_em]) / (2 * pow(1. - z_coll[pair_em] + mu2_Q_x[pair_em][pair_sp], 2)) // from first term in (5.58, J_gQ^a) included from (6.58) 
	    - 2 / (1. - z_coll[pair_em]) * (1. + log(1. - z_coll[pair_em] + mu2_Q_x[pair_em][pair_sp])) // from second term in (5.58, J_gQ^a) included from (6.58)
	    ;
	  
	  Kit_plus_1[i_dt][pair_em][pair_sp] = 0.
	    + 2 * log(1. - z_coll[pair_em]) / (1. - z_coll[pair_em]) // from first term in (6.58) 
	    + (1. - z_coll[pair_em]) / (2 * pow(1. - z_coll[pair_em] + mu2_Q[pair_em][pair_sp], 2)) // from first term in (5.58, J_gQ^a) included from (6.58) 
	    - 2 / (1. - z_coll[pair_em]) * (1. + log(1. - z_coll[pair_em] + mu2_Q[pair_em][pair_sp])) // from second term in (5.58, J_gQ^a) included from (6.58)
	    ;
	  
	  intKit_plus[i_dt][pair_em][pair_sp] = 0.
	    - pow(log(1. - x_pdf[pair_em]), 2) // from first term in (6.58)
	    + .5 * (- mu2_Q[pair_em][pair_sp] / (1. - x_pdf[pair_em] + mu2_Q[pair_em][pair_sp]) + mu2_Q[pair_em][pair_sp] / (1. + mu2_Q[pair_em][pair_sp]) - log((1. - x_pdf[pair_em] + mu2_Q[pair_em][pair_sp]) / (1. + mu2_Q[pair_em][pair_sp]))) // from first term in (5.58) included from (6.58)
	    + 2 * (gsl_sf_dilog(-1. / mu2_Q[pair_em][pair_sp]) - gsl_sf_dilog(-(1. - x_pdf[pair_em]) / mu2_Q[pair_em][pair_sp]) + log(1. - x_pdf[pair_em]) * (1. + log(mu2_Q[pair_em][pair_sp]))); // from second term in (5.58) included from (6.58)
	  
	  // terms containing x-dependent pre-factor of (2/(1-z_coll[pair_em]))_+
	  Kit_plus_outside_x[i_dt][pair_em][pair_sp] = 
	    + log(((2. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp]) / ((2. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp] + m2_Q[pair_sp])) // from fourth term from (6.58)
	    + log(2. + mu2_Q_x[pair_em][pair_sp] - z_coll[pair_em]) // from third term in (5.58) included from (6.58)
	    ;
	  Kit_plus_outside_1[i_dt][pair_em][pair_sp] = 
	    + log(sall_ja[pair_em][pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp])) // from fourth term from (6.58)
	    + log(1. + mu2_Q[pair_em][pair_sp]) // from third term in (5.58) included from (6.58)
	    ;
	  Kit_delta[i_dt][pair_em][pair_sp] = 
	    - gamma_qew_q // from fifth term from (6.58)
	    + mu2_Q[pair_em][pair_sp] * log(m2_Q[pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp])) // from sixth term from (6.58)
	    + .5 * m2_Q[pair_sp] / (sall_ja[pair_em][pair_sp] + m2_Q[pair_sp]) // from seventh term from (6.58)
	    ;
	}
	else if (i_dt == 2 && (CA_dipole_splitting[i_dt][1] == 1 || CA_dipole_splitting[i_dt][2] == 1)){
	  // q -> g (+q) splitting
	  Kit[i_dt][pair_em][pair_sp] = 
	    2 * 0.5 * m2_Q[pair_sp] / (z_coll[pair_em] * sall_ja_x[pair_em][pair_sp]) * log(m2_Q[pair_sp] / ((1. - z_coll[pair_em]) * sall_ja_x[pair_em][pair_sp] + m2_Q[pair_sp])); // from (6.59) // C_F / C_A -> 0.5 ???
	}
	else if (i_dt == 3 && (CA_dipole_splitting[i_dt][1] == 1 || CA_dipole_splitting[i_dt][2] == 1)){
	  // g -> q (+g) splitting
	  //	  Kit[i_dt][pair_em][pair_sp] = 0.;
	}
      }
    }


    /*
    else if (dipole_phasespace[0][pair_sp + 6] > 2 && dx_pa[dipole_phasespace[0][pair_sp + 6]][0] == 0){
      // !!! does maybe not vanish in particular cases with outgoing gluons !!!
      // j == gluon case: additional terms for N_J^ja !!!
    }
    */

      
    if (switch_TSV){
      for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	for (int i_m = 0; i_m < n_scale_dyn_fact[v_sf]; i_m++){
	  value_logscale2_fact_papi[v_sf][i_m][pair_em][pair_sp] = value_central_logscale2_fact[v_sf] + value_relative_logscale2_fact[v_sf][i_m] - ln_papi[pair_em][pair_sp];
	}
      }
    }

    for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
      for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	CA_value_ln_muF_papi[sd][ss][pair_em][pair_sp] = CA_value_log_mu2_fact[sd][ss] - ln_papi[pair_em][pair_sp];
      }
    }
  }
  logger << LOG_DEBUG_VERBOSE << "pair finished!" << endl;

  /*
  for (int i_dt = 0; i_dt < 4; i_dt++){
    cout << "CA_dipole_splitting[" << i_dt << "] = ";
    for (int i_em = 0; i_em < 3; i_em++){
      cout << CA_dipole_splitting[i_dt][i_em] << "   ";
    }
    cout << endl;
  }

  for (int i_dt = 0; i_dt < 4; i_dt++){
    for (int i_p = 0; i_p < pair.size(); i_p++){
      int pair_em = pair[i_p][0];
      int pair_sp = pair[i_p][1];
      cout << "Kit                [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit[i_dt][pair_em][pair_sp] << endl;
      cout << "Kit_plus_x         [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit_plus_x[i_dt][pair_em][pair_sp] << endl;
      cout << "Kit_plus_1         [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit_plus_1[i_dt][pair_em][pair_sp] << endl;
      cout << "intKit_plus        [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << intKit_plus[i_dt][pair_em][pair_sp] << endl;
      cout << "Kit_plus_outside_x [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit_plus_outside_x[i_dt][pair_em][pair_sp] << endl;
      cout << "Kit_plus_outside_1 [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit_plus_outside_1[i_dt][pair_em][pair_sp] << endl;
      cout << "Kit_delta          [" << i_dt << "][" << pair_em << "][" << pair_sp << "] = " << setw(23) << setprecision(15) << Kit_delta[i_dt][pair_em][pair_sp] << endl;
    }
  }
  */  

  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
      // K terms 
      if (no_spectator[i_c][j_c] == 0){
	data_K[0][i_c][j_c] = (Kbar[type[i_c][j_c]][no_emitter[i_c][j_c]] + Kbar_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
	data_K[1][i_c][j_c] = (-Kbar_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];
	data_K[2][i_c][j_c] = (Kbar_delta[type[i_c][j_c]] - intKbar_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];
      }
      else if (no_spectator[i_c][j_c] < 3){
	data_K[0][i_c][j_c] = -(Kt[type[i_c][j_c]][no_emitter[i_c][j_c]] + Kt_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];// reg
	data_K[1][i_c][j_c] = -(-Kt_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];// plus
	data_K[2][i_c][j_c] = -(Kt_delta[type[i_c][j_c]] - intKt_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_ME2_cf[i_c][j_c];// delta
      }
      else if (no_spectator[i_c][j_c] > 2){
	if (m2_Q[no_spectator[i_c][j_c]] == 0.){
	  if (collinear_singularity[i_c][j_c] == 1){
	    data_K[0][i_c][j_c] = gamma_i[i_c][j_c] * (1. / (1. - z_coll[no_emitter[i_c][j_c]])) * CA_ME2_cf[i_c][j_c];
	    data_K[1][i_c][j_c] = -gamma_i[i_c][j_c] * (1. / (1. - z_coll[no_emitter[i_c][j_c]])) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c];
	    data_K[2][i_c][j_c] = gamma_i[i_c][j_c] * (1. + log(1. - x_pdf[no_emitter[i_c][j_c]])) * CA_ME2_cf[i_c][j_c];
	  }
	}
	else {
	  // -> Kit
	  data_K[0][i_c][j_c] = -(Kit[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] 
					+ Kit_plus_x[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] 
					+ 2. / (1. - z_coll[no_emitter[i_c][j_c]]) * Kit_plus_outside_x[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]
					) * CA_ME2_cf[i_c][j_c]; // reg
	  data_K[1][i_c][j_c] =  -(
					 - Kit_plus_1[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] 
					 - 2. / (1. - z_coll[no_emitter[i_c][j_c]]) * Kit_plus_outside_1[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]
					 ) * z_coll[no_emitter[i_c][j_c]] * CA_ME2_cf[i_c][j_c]; // plus
	  data_K[2][i_c][j_c] = -(Kit_delta[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] 
					- intKit_plus[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]
					- 2 * log(1. - x_pdf[no_emitter[i_c][j_c]]) * Kit_plus_outside_1[type[i_c][j_c]][no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]
					) * CA_ME2_cf[i_c][j_c]; // delta
	  
	  data_K[0][i_c][j_c] += -P[type[i_c][j_c]][no_emitter[i_c][j_c]] * log(((1. - z_coll[no_emitter[i_c][j_c]]) * sall_ja_x[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]) / ((1. - z_coll[no_emitter[i_c][j_c]]) * sall_ja_x[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] + m2_Q[no_spectator[i_c][j_c]])) * CA_ME2_cf[i_c][j_c]; // from (6.55)
	  // no explizit plus term in (6.55)
	  //	  cout << "collinear_singularity[" << i_c << "][" << j_c << "] = " << collinear_singularity[i_c][j_c] << endl;
	  if (collinear_singularity[i_c][j_c] == 1){
	    // contributes only to irregular splittings 
	    //	    cout << "data_K[2][" << i_c << "][" << j_c << "] = " << data_K[2][i_c][j_c] << endl;
	    data_K[2][i_c][j_c] += -gamma_ax[no_emitter[i_c][j_c]] * (log((sall_ja[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] - 2 * m_Q[no_spectator[i_c][j_c]] * sqrt(sall_ja[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] + m2_Q[no_spectator[i_c][j_c]]) + 2 * m2_Q[no_spectator[i_c][j_c]]) / sall_ja[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]]) + 2 * m_Q[no_spectator[i_c][j_c]] / (sqrt(sall_ja[no_emitter[i_c][j_c]][no_spectator[i_c][j_c]] + m2_Q[no_spectator[i_c][j_c]]) + m_Q[no_spectator[i_c][j_c]])) * CA_ME2_cf[i_c][j_c];  // from (6.55)
	    //	    data_K[2][i_c][i_ca] += -gamma_a_T2_ap[i_em] * (log((sall_ja[i_em][i_cs] - 2 * m_j * sqrt(sall_ja[i_em][i_cs] + m2_j) + 2 * m2_j) / sall_ja[i_em][i_cs]) + 2 * m_j / (sqrt(sall_ja[i_em][i_cs] + m2_j) + m_j)) * CA_ME2cc[i_c][i_ca];  // from (6.55)
	    //	    cout << "data_K[2][" << i_c << "][" << j_c << "] = " << data_K[2][i_c][j_c] << endl;

	  }
	}
      }
      if (no_spectator[i_c][j_c] != 0){
	// P terms 
	if (switch_TSV){
	  for (int v_sf = 0; v_sf < value_logscale2_fact_papi.size(); v_sf++){
	    for (int v_xf = 0; v_xf < value_logscale2_fact_papi[v_sf].size(); v_xf++){
	      value_data_P[v_sf][v_xf][0][i_c][j_c] = (P[type[i_c][j_c]][no_emitter[i_c][j_c]] + P_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	      value_data_P[v_sf][v_xf][1][i_c][j_c] = (-P_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	      value_data_P[v_sf][v_xf][2][i_c][j_c] = (P_delta[type[i_c][j_c]] - intP_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * value_logscale2_fact_papi[v_sf][v_xf][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	    }
	  }
	}
	for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
	  for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	    value_dataP[sd][ss][0][i_c][j_c] = (P[type[i_c][j_c]][no_emitter[i_c][j_c]] + P_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	    value_dataP[sd][ss][1][i_c][j_c] = (-P_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * z_coll[no_emitter[i_c][j_c]] * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	    value_dataP[sd][ss][2][i_c][j_c] = (P_delta[type[i_c][j_c]] - intP_plus[type[i_c][j_c]][no_emitter[i_c][j_c]]) * CA_value_ln_muF_papi[sd][ss][ppair[i_c][j_c][0]][ppair[i_c][j_c][1]] * CA_ME2_cf[i_c][j_c];
	  }
	}
      }
    }
    for (int i_x = 0; i_x < 3; i_x++){
      double temp_sumK = accumulate(data_K[i_x][i_c].begin(), data_K[i_x][i_c].end(), 0.);
      if (switch_TSV){
	for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	  for (int v_xf = 0; v_xf < value_logscale2_fact_papi[v_sf].size(); v_xf++){
	    if (switch_KP == 0){
	      // K + P terms
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_e_2pi * (temp_sumK + accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_e_2pi * (temp_sumK + accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	    }
	    else if (switch_KP == 1){
	      // P terms 
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_e_2pi * (accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_e_2pi * (accumulate(value_data_P[v_sf][v_xf][i_x][i_c].begin(), value_data_P[v_sf][v_xf][i_x][i_c].end(), 0.));
	    }
	    else if (switch_KP == 2){
	      // K terms 
	      value_ME2_KP[v_sf][v_xf][i_x][i_c] = alpha_e_2pi * (temp_sumK);
	      value_ME2term_fact[i_c][i_x][v_sf][v_xf] = alpha_e_2pi * (temp_sumK);
	    }
	  }
	}
      }
      //    }
      //    dataK = data_K;
      //    for (int i_x = 0; i_x < 3; i_x++){
      for (int sd = 0; sd < CA_value_ln_muF_papi.size(); sd++){
	for (int ss = 0; ss < CA_value_ln_muF_papi[sd].size(); ss++){
	  if (switch_KP == 0){
	    // K + P term
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_e_2pi * (temp_sumK + accumulate(value_dataP[sd][ss][i_x][i_c].begin(), value_dataP[sd][ss][i_x][i_c].end(), 0.));
	  }
	  else if (switch_KP == 1){
	    // P terms 
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_e_2pi * (accumulate(value_dataP[sd][ss][i_x][i_c].begin(), value_dataP[sd][ss][i_x][i_c].end(), 0.));
	  }
	  else if (switch_KP == 2){
	    // K terms 
	    CA_value_ME2_KP[sd][ss][i_x][i_c] = alpha_e_2pi * (temp_sumK);
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::output_collinear(){
  static Logger logger("observable_set::output_collinear");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
      stringstream temp;
      if      ((*CA_collinear)[i_c][j_c].type_correction() == 1){temp << "QCD ";}
      else if ((*CA_collinear)[i_c][j_c].type_correction() == 2){temp << "QEW ";}
      temp << "collinear dipole " << i_c << ", " << j_c << ":   " << setw(15) << left << (*CA_collinear)[i_c][j_c].name() << ":   ";
      temp << "csi->type_parton = ";
      for (int i_p = 1; i_p < (*CA_collinear)[i_c][j_c].type_parton().size(); i_p++){temp << (*CA_collinear)[i_c][j_c].type_parton()[i_p] << "   ";}
      //      for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){temp << (*CA_collinear)[i_c][j_c].type_parton()[i_p] << "   ";}
      temp << "no_prc = " << (*CA_collinear)[i_c][j_c].no_prc() << "   ";
      logger << LOG_INFO << temp.str() << endl;
      temp.str("");
      temp << setw(48) << "";
      temp << "charge_factor = " << setw(23) << setprecision(15) << (*CA_collinear)[i_c][j_c].charge_factor() << "   ";
      temp << "massive = " << (*CA_collinear)[i_c][j_c].massive() << "   ";
      logger << LOG_INFO << temp.str() << endl;
      temp.str("");
      temp << setw(48) << "";
      temp << "type = " << (*CA_collinear)[i_c][j_c].type() << "   ";
      temp << "in_collinear = ";
      for (int i_em = 0; i_em < 3; i_em++){temp << (*CA_collinear)[i_c][j_c].in_collinear()[i_em] << "   ";}
      logger << LOG_INFO << temp.str() << endl;
      temp.str("");
      temp << setw(48) << "";
      temp << "no_emitter = " << (*CA_collinear)[i_c][j_c].no_emitter() << "   ";
      temp << "no_spectator = " << (*CA_collinear)[i_c][j_c].no_spectator() << "   ";
      temp << "pair = ";
      for (int i_p = 0; i_p < 2; i_p++){temp << (*CA_collinear)[i_c][j_c].pair()[i_p] << "   ";}
      logger << LOG_INFO << temp.str() << endl;
      for (int i_x = 0; i_x < (*CA_collinear)[i_c][j_c].all_name().size(); i_x++){
	temp.str("");
	temp << setw(48) << "";
	temp << "pdf contributions: " << setw(15) << left << (*CA_collinear)[i_c][j_c].all_name()[i_x] << ":   ";
	for (int i_yxy = 0; i_yxy < 3; i_yxy++){temp << setw(2) << right << (*CA_collinear)[i_c][j_c].all_pdf()[i_x][i_yxy] << "  ";}
	logger << LOG_INFO << temp.str() << endl;
      }
    }
  }
  logger.newLine(LOG_INFO);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::output_collinear_pdf(){
  static Logger logger("observable_set::output_collinear_pdf");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      logger << LOG_INFO << "pdf contributions: " << setw(15) << left << (*CA_collinear)[i_c][0].all_name()[i_i] << ":   " << endl;
      for (int i_x = 0; i_x < CA_combination_pdf[i_c][i_i].size(); i_x++){
	stringstream temp;
	temp.str("");
	temp << setw(19) << "";
	temp << "CA_Q2f[" << setw(2) << i_c << "][" << setw(2) << i_i << "][" << setw(2) << i_x << "] = " << setw(10) << setprecision(8) << CA_Q2f[i_c][i_i][i_x] << "     "; 
	temp << "CA_combination_pdf[" << setw(2) << i_c << "][" << setw(2) << i_i << "][" << setw(2) << i_x << "] = "; 
	for (int i_y = 0; i_y < CA_combination_pdf[i_c][i_i][i_x].size(); i_y++){temp << setw(2) << right << CA_combination_pdf[i_c][i_i][i_x][i_y] << "  ";}
	logger << LOG_INFO << temp.str() << endl;
      }
    }
    logger.newLine(LOG_INFO);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
