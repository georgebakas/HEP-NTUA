#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
void observable_set::determine_CX_QCD(phasespace_set & psi, int n_emission){
  Logger logger("observable_set::determine_CX_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  /*
  static int initialization = 1;
  //  static map<int, double> charge_particle;
  if (initialization == 1){
    //    fill_charge_particle(charge_particle);
    initialization = 0;
  }
  */
  logger << LOG_DEBUG << "multicollinear determination" << endl;
  /*
  for (int i1 = 0; i1 < combination_pdf.size(); i1++){
    stringstream temp_ss;
    temp_ss << " combination_pdf[" << setw(2) << i1 << "] = "; 
    for (int i = 0; i < 3; i++){temp_ss << setw(4) << combination_pdf[i1][i] << "   ";}
    logger << LOG_DEBUG << temp_ss.str() << endl;
  }
  */
  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "csi->type_parton[0].size() = " << csi->type_parton[0].size() << endl;

  vector<string> pa_name(csi->type_parton[0].size(), "");
  logger << LOG_DEBUG << "pa_name.size() = " << pa_name.size() << endl;
  if (csi->type_parton[0][1] >= -10 && csi->type_parton[0][1] <= 10){pa_name[1] = "a";}
  if (csi->type_parton[0][2] >= -10 && csi->type_parton[0][2] <= 10){pa_name[2] = "b";}
  //  int count = 0;
  //  vector<string> alphabet(csi->type_parton[0].size() - 3, "");
  //  for (int i_p = 0; i_p < alphabet.size(); i_p++){alphabet[i_p] = char(105 + i_p);}
  for (int i_p = 3; i_p < pa_name.size(); i_p++){
    logger << LOG_DEBUG << "i_p = " << i_p << endl;
    if (csi->type_parton[0][i_p] >= -10 && csi->type_parton[0][i_p] <= 10){pa_name[i_p] = char(105 + i_p - 3);}
  } 
  //  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (csi->type_parton[0][i_p] >= -10 && csi->type_parton[0][i_p] =< 10){pa_name[i_p] = alphabet[i_p - 3];}} 
  //  new labels: omit colourless particles in counting
  //  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (csi->type_parton[0][i_p] >= -10 && csi->type_parton[0][i_p] =< 10){pa_name[i_p] = alphabet[count++];}}

  //  int temp_type_correction = 1;
  /*
  logger << LOG_DEBUG << "before pname = " << endl;

  map <int,string> pname;
  fill_pname(pname);
  */
  logger << LOG_DEBUG << "psi_no_prc[0] = " << psi_no_prc[0] << endl;
  logger << LOG_DEBUG << "csi->type_parton[0].size() = " << csi->type_parton[0].size() << endl;
  logger << LOG_DEBUG << "combination_pdf.size() = " << combination_pdf.size() << endl;

  multicollinear_set start_multicollinear(psi_no_prc[0], csi->type_parton[0], combination_pdf);
  //  vector<vector<multicollinear_set> > 
  multicollinear.resize(1);
  multicollinear[0].push_back(start_multicollinear);

  //string _name, vector<string> _all_name, vector<int> _type_splitting, vector<vector<int> > _in_collinear, int _no_prc, vector<int> _csi->type_parton, vector<vector<int> > _all_pdf, double _charge_factor, vector<int> _no_emitter, vector<int> _no_spectator, vector<vector<int> > _pair, vector<int> _type_correction, vector<int> _massive

  for (int i_e = 1; i_e < n_emission + 1; i_e++){
    logger << LOG_DEBUG << "i_e = " << i_e << endl;
    //    logger << LOG_DEBUG << "multicollinear[i_e - 1 = " << i_e - 1 << "].size() = " << multicollinear[i_e - 1].size() << endl;
    multicollinear.push_back(vector<multicollinear_set> (0));
    for (int i_c = 0; i_c < multicollinear[i_e - 1].size(); i_c++){
      logger << LOG_DEBUG << "i_c = " << i_c << endl;
      //      logger << LOG_DEBUG << "i_c = " << i_c << endl;
      //      logger << LOG_DEBUG << "multicollinear.size() = " << multicollinear.size() << endl;
      //      logger << LOG_DEBUG << "multicollinear.size() = " << multicollinear.size() << endl;
      //      logger << LOG_DEBUG << "multicollinear[i_e = " << i_e << "].size() = " << multicollinear[i_e].size() << endl;
      int x_c = 0;      
      for (int temp_no_emitter = 1; temp_no_emitter < 3; temp_no_emitter++){
	logger << LOG_DEBUG << "temp_no_emitter = " << temp_no_emitter << endl;
	if (pa_name[temp_no_emitter] == ""){continue;}
	if (temp_no_emitter < multicollinear[i_e - 1][i_c].no_emitter[i_e - 1]){continue;}
	for (int i_t = 0; i_t < 2; i_t++){
	  logger << LOG_DEBUG << "i_t = " << i_t << endl;

	  x_c = multicollinear[i_e].size();
	  multicollinear[i_e].push_back(multicollinear[i_e - 1][i_c]);
	  multicollinear[i_e][x_c].no_emitter.push_back(temp_no_emitter);
	  multicollinear[i_e][x_c].in_collinear.push_back(vector<int> (3, 0));
	  multicollinear[i_e][x_c].in_collinear[i_e][temp_no_emitter] = 1;
	  multicollinear[i_e][x_c].charge_factor.push_back(1.);
	  multicollinear[i_e][x_c].type_correction.push_back(1);
	  multicollinear[i_e][x_c].previous_splitting.push_back(i_c);

	  multicollinear[i_e][x_c].n_emission[0]++;
	  multicollinear[i_e][x_c].n_emission[temp_no_emitter]++;
	  multicollinear[i_e][x_c].emission[temp_no_emitter] = 1;

	  multicollinear[i_e][x_c].n_emission_leg_type[0][0]++;
	  multicollinear[i_e][x_c].n_emission_leg_type[0][1]++;
	  multicollinear[i_e][x_c].n_emission_leg_type[temp_no_emitter][0]++;
	  multicollinear[i_e][x_c].n_emission_leg_type[temp_no_emitter][1]++;

	  logger << LOG_DEBUG << "x_c = " << x_c << endl;
	  logger << LOG_DEBUG << "temp_no_emitter = " << temp_no_emitter << endl;
	  logger << LOG_DEBUG << "multicollinear.size() = " <<multicollinear.size()  << endl;
	  logger << LOG_DEBUG << "multicollinear[" << i_e << "].size() = " << multicollinear[i_e].size()  << endl;
	  logger << LOG_DEBUG << "multicollinear[" << i_e << "][" << x_c << "].pdf.size() = " << multicollinear[i_e][x_c].pdf.size()  << endl;
	  logger << LOG_DEBUG << "multicollinear[" << i_e << "][" << x_c << "].pdf[0].size() = " << multicollinear[i_e][x_c].pdf[0].size()  << endl;
	  logger << LOG_DEBUG << "multicollinear[" << i_e << "][" << x_c << "].pdf[0][" << temp_no_emitter << "] = " << multicollinear[i_e][x_c].pdf[0][temp_no_emitter] << endl;
	  if (multicollinear[i_e][x_c].pdf[0][temp_no_emitter] == 0){
	  logger << LOG_DEBUG << "   i_t = " << i_t << endl;
	    if (i_t == 0){
	      multicollinear[i_e][x_c].type_splitting.push_back(0); // hard process with g from g -> g g splitting
	      multicollinear[i_e][x_c].type_splitting_leg[temp_no_emitter].push_back(1); // hard process with g from g -> g g splitting
	      multicollinear[i_e][x_c].in_collinear[i_e][0] = 1;
	      multicollinear[i_e][x_c].name = multicollinear[i_e][x_c].name + "_" + pa_name[temp_no_emitter] + "^{gg}";
	      for (int i_x = 0; i_x < multicollinear[i_e][x_c].pdf.size(); i_x++){multicollinear[i_e][x_c].all_name[i_x] = multicollinear[i_e][x_c].all_name[i_x] + "_" + pa_name[temp_no_emitter] + "^{gg}";}
	    }
	    else {
	      multicollinear[i_e][x_c].type_splitting.push_back(2); // hard process with g from q -> g q splitting
	      multicollinear[i_e][x_c].type_splitting_leg[temp_no_emitter].push_back(3); // hard process with g from q -> g q splitting
	      multicollinear[i_e][x_c].in_collinear[i_e][0] = 0;
	      multicollinear[i_e][x_c].name = multicollinear[i_e][x_c].name + "_" + pa_name[temp_no_emitter] + "^{qg}";
	      for (int i_x = 0; i_x < multicollinear[i_e][x_c].pdf.size(); i_x++){multicollinear[i_e][x_c].all_name[i_x] = multicollinear[i_e][x_c].all_name[i_x] + "_" + pa_name[temp_no_emitter] + "^{qg}";}
	      for (int i_x = 0; i_x < multicollinear[i_e][x_c].pdf.size(); i_x++){multicollinear[i_e][x_c].pdf[i_x][temp_no_emitter] = 10;}
	    }
	  logger << LOG_DEBUG << "   i_t = " << i_t << endl;
	  }
	  else if (multicollinear[i_e][x_c].pdf[0][temp_no_emitter] != 0){
	    if (i_t == 0){
	      multicollinear[i_e][x_c].type_splitting.push_back(3); // hard process with q from g -> q qx splitting
	      multicollinear[i_e][x_c].type_splitting_leg[temp_no_emitter].push_back(4); // hard process with q from g -> q qx splitting
	      multicollinear[i_e][x_c].in_collinear[i_e][0] = 0;
	      multicollinear[i_e][x_c].name = multicollinear[i_e][x_c].name + "_" + pa_name[temp_no_emitter] + "^{g" + csi->name_particle[multicollinear[i_e][x_c].pdf[0][temp_no_emitter]] + "}";
	      for (int i_x = 0; i_x < multicollinear[i_e][x_c].pdf.size(); i_x++){multicollinear[i_e][x_c].all_name[i_x] = multicollinear[i_e][x_c].all_name[i_x] + "_" + pa_name[temp_no_emitter] + "^{g" + csi->name_particle[multicollinear[i_e][x_c].pdf[i_x][temp_no_emitter]] + "}";}
	      for (int i_x = 0; i_x < multicollinear[i_e][x_c].pdf.size(); i_x++){multicollinear[i_e][x_c].pdf[i_x][temp_no_emitter] = 0;}
	    }
	    else {
	      multicollinear[i_e][x_c].type_splitting.push_back(1); // hard process with q from q -> q g splitting
	      multicollinear[i_e][x_c].type_splitting_leg[temp_no_emitter].push_back(2); // hard process with q from q -> q g splitting
	      multicollinear[i_e][x_c].in_collinear[i_e][0] = 1;
	      multicollinear[i_e][x_c].name = multicollinear[i_e][x_c].name + "_" + pa_name[temp_no_emitter] + "^{" + csi->name_particle[multicollinear[i_e][x_c].pdf[0][temp_no_emitter]] + csi->name_particle[multicollinear[i_e][x_c].pdf[0][temp_no_emitter]] + "}";
	      for (int i_x = 0; i_x < multicollinear[i_e][x_c].pdf.size(); i_x++){multicollinear[i_e][x_c].all_name[i_x] = multicollinear[i_e][x_c].all_name[i_x] + "_" + pa_name[temp_no_emitter] + "^{" + csi->name_particle[multicollinear[i_e][x_c].pdf[i_x][temp_no_emitter]] + csi->name_particle[multicollinear[i_e][x_c].pdf[i_x][temp_no_emitter]] + "}";}
	    }
	  }
	  /*
	// already existing: mass_parton or so... ???
	int temp_massive;
	if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 0;}
	else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 1;}
	else if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 2;}
	else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 3;}
	else {cout << "Should not happen!" << endl;}
	multicollinear[i_e][x_c].massive.push_back(temp_massive);
	  */
	}

	// non-singlet contribution
	if (i_e < 2){continue;}
	if (temp_no_emitter != multicollinear[i_e - 1][i_c].no_emitter[i_e - 1]){continue;}
	if (multicollinear[i_e - 1][i_c].type_splitting[i_e - 1] != 3){continue;}

	x_c = multicollinear[i_e].size();
	logger << LOG_DEBUG << "x_c = " << x_c << endl;
	multicollinear[i_e].push_back(multicollinear[i_e - 1][i_c]);
	multicollinear[i_e][x_c].no_emitter.push_back(temp_no_emitter);
	multicollinear[i_e][x_c].in_collinear.push_back(vector<int> (3, 0));
	multicollinear[i_e][x_c].in_collinear[i_e][temp_no_emitter] = 1;
	multicollinear[i_e][x_c].charge_factor.push_back(1.);
	multicollinear[i_e][x_c].type_correction.push_back(1);
	multicollinear[i_e][x_c].previous_splitting.push_back(i_c);

	multicollinear[i_e][x_c].n_emission[0]++;
	multicollinear[i_e][x_c].n_emission[temp_no_emitter]++;
	multicollinear[i_e][x_c].emission[temp_no_emitter] = 1;

	multicollinear[i_e][x_c].type_splitting.push_back(4); // selected(a-q) splitting
	multicollinear[i_e][x_c].type_splitting_leg[temp_no_emitter].push_back(5); // selected(a-q) splitting
	multicollinear[i_e][x_c].in_collinear[i_e][0] = 0;
	multicollinear[i_e][x_c].name = multicollinear[i_e][x_c].name + "_" + pa_name[temp_no_emitter] + "^{" + csi->name_particle[-multicollinear[i_e - 2][multicollinear[i_e - 1][i_c].previous_splitting[i_e - 2]].pdf[0][temp_no_emitter]] + "g}";
	for (int i_x = 0; i_x < multicollinear[i_e][x_c].pdf.size(); i_x++){multicollinear[i_e][x_c].all_name[i_x] = multicollinear[i_e][x_c].all_name[i_x] + "_" + pa_name[temp_no_emitter] + "^{" + csi->name_particle[-multicollinear[i_e - 2][multicollinear[i_e - 1][i_c].previous_splitting[i_e - 2]].pdf[i_x][temp_no_emitter]] + "g}";}
	for (int i_x = 0; i_x < multicollinear[i_e][x_c].pdf.size(); i_x++){multicollinear[i_e][x_c].pdf[i_x][temp_no_emitter] = -multicollinear[i_e - 2][multicollinear[i_e - 1][i_c].previous_splitting[i_e - 2]].pdf[i_x][temp_no_emitter];}
	// END non-singlet contribution

      }

      /*
      int flag = (*CA_collinear).size();
      for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){if (temp_no_emitter == (*CA_collinear)[i_a][0].no_emitter() && temp_type == (*CA_collinear)[i_a][0].type()){flag = i_a; break;}}
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
      */
    }

    // sort out splitting sequence leading to identical pdf's:
    //	  logger << LOG_DEBUG << "i_e = " << i_e << << endl;

    if (i_e < 2){continue;}
    logger << LOG_DEBUG << "sort out splitting sequence leading to identical pdf's:" << endl;
    for (int i_c = 0; i_c < multicollinear[i_e].size(); i_c++){
      // hard process with q from  g -> q (qx) -> q (g) [42] is understood as part of  g -> g (g) -> q (qx) [14]
      if (multicollinear[i_e][i_c].type_splitting[i_e] == 0 && multicollinear[i_e][i_c].type_splitting[i_e - 1] == 3){
	//	for (int j_c = i_c + 1; j_c < multicollinear[i_e].size(); j_c++){
	for (int j_c = 0; j_c < multicollinear[i_e].size(); j_c++){
	  if (j_c == i_c){continue;}
	  /*
	  logger << LOG_DEBUG << "i_c = " << i_c << "   j_c = " << j_c << endl;
	  logger << LOG_DEBUG << setw(15) << "" << "multicollinear[" << i_e << "][" << j_c << "].type_splitting[" << i_e << "] = " << multicollinear[i_e][j_c].type_splitting[i_e] << endl;
	  logger << LOG_DEBUG << setw(15) << "" << "multicollinear[" << i_e << "][" << j_c << "].type_splitting[" << i_e - 1 << "] = " << multicollinear[i_e][j_c].type_splitting[i_e - 1] << endl;
	  logger << LOG_DEBUG << setw(15) << "" << "(multicollinear[i_e][i_c].no_emitter[i_e] == multicollinear[i_e][j_c].no_emitter[i_e]) = " << (multicollinear[i_e][i_c].no_emitter[i_e] == multicollinear[i_e][j_c].no_emitter[i_e]) << endl;
	  logger << LOG_DEBUG << setw(15) << "" << "(multicollinear[i_e][i_c].no_emitter[i_e - 1] == multicollinear[i_e][j_c].no_emitter[i_e - 1]) = " << (multicollinear[i_e][i_c].no_emitter[i_e - 1] == multicollinear[i_e][j_c].no_emitter[i_e - 1]) << endl;
	  */
	  if (multicollinear[i_e][j_c].type_splitting[i_e] == 3 && 
	      multicollinear[i_e][j_c].type_splitting[i_e - 1] == 1 && 
	      multicollinear[i_e][i_c].no_emitter[i_e] == multicollinear[i_e][j_c].no_emitter[i_e] && 
	      multicollinear[i_e][i_c].no_emitter[i_e - 1] == multicollinear[i_e][j_c].no_emitter[i_e - 1]){
	    //	    logger << LOG_DEBUG << "i_c = " << i_c << "   j_c = " << j_c << "   erase" << endl;

	    multicollinear[i_e].erase(multicollinear[i_e].begin() + j_c);
	    break;
	  }
	}
      }
      // hard process with q from  g from Q -> g (Q) -> g (g) [31] is understood as part of  Q -> Q (g) -> g (Q) [23]
      else if (multicollinear[i_e][i_c].type_splitting[i_e] == 1 && multicollinear[i_e][i_c].type_splitting[i_e - 1] == 2){
	//	for (int j_c = i_c + 1; j_c < multicollinear[i_e].size(); j_c++){
	for (int j_c = 0; j_c < multicollinear[i_e].size(); j_c++){
	  if (j_c == i_c){continue;}
	  if (multicollinear[i_e][j_c].type_splitting[i_e] == 2 && 
		   multicollinear[i_e][j_c].type_splitting[i_e - 1] == 0 && 
		   multicollinear[i_e][i_c].no_emitter[i_e] == multicollinear[i_e][j_c].no_emitter[i_e] && 
		   multicollinear[i_e][i_c].no_emitter[i_e - 1] == multicollinear[i_e][j_c].no_emitter[i_e - 1]){
	    //	    logger << LOG_DEBUG << "i_c = " << i_c << "   j_c = " << j_c << "   erase" << endl;

	    multicollinear[i_e].erase(multicollinear[i_e].begin() + j_c);
	    break;
	  }
	}
      }
      // hard process with g from g -> Q (Qx) -> g (Q) [43] is understood as part of  g -> g (g) -> g (g) [11]
      else if (multicollinear[i_e][i_c].type_splitting[i_e] == 0 && multicollinear[i_e][i_c].type_splitting[i_e - 1] == 0){
	//	for (int j_c = i_c + 1; j_c < multicollinear[i_e].size(); j_c++){
	for (int j_c = 0; j_c < multicollinear[i_e].size(); j_c++){
	  if (j_c == i_c){continue;}
	  if (multicollinear[i_e][j_c].type_splitting[i_e] == 3 && 
		   multicollinear[i_e][j_c].type_splitting[i_e - 1] == 2 && 
		   multicollinear[i_e][i_c].no_emitter[i_e] == multicollinear[i_e][j_c].no_emitter[i_e] && 
		   multicollinear[i_e][i_c].no_emitter[i_e - 1] == multicollinear[i_e][j_c].no_emitter[i_e - 1]){
	    //	    logger << LOG_DEBUG << "i_c = " << i_c << "   j_c = " << j_c << "   erase" << endl;

	    multicollinear[i_e].erase(multicollinear[i_e].begin() + j_c);
	    break;
	  }
	}
      }
    }
  }

  for (int i_e = 0; i_e < multicollinear.size(); i_e++){
    logger << LOG_DEBUG << "mc[i_e = " << i_e << "].size() = " << multicollinear[i_e].size() << endl;
    logger.newLine(LOG_DEBUG);
    for (int i_c = 0; i_c < multicollinear[i_e].size(); i_c++){
      //  for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
      //    for (int j_a = 0; j_a < (*CA_collinear)[i_a].size(); j_a++){
      stringstream temp;
      temp << "mc " << i_e << ", " << i_c << ":   " << setw(1 + 9 * i_e) << left << multicollinear[i_e][i_c].name << ":   ";
      temp << "csi->type_parton = (";
      for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){
	temp << setw(4) << right << multicollinear[i_e][i_c].type_parton[i_p];
	if (i_p == 2){temp << " ->";}
      }
      temp << ")   ";
      temp << "no_prc = " << multicollinear[i_e][i_c].no_prc << "   ";

      temp << "ps = " << multicollinear[i_e][i_c].previous_splitting[i_e] << "   ";
      temp << "n_emission = ";
      temp << "(";
      for (int i_z = 0; i_z < 3; i_z++){
	temp << multicollinear[i_e][i_c].n_emission[i_z];
	if (i_z < 2){temp << " ";}
      }
      temp << ")";
      logger << LOG_DEBUG << temp.str() << endl;
      temp.str("");

      for (int j_e = 1; j_e < i_e + 1; j_e++){
	temp << setw(16 + 9 * i_e) << "";
	temp << "emission " << j_e << " :   ";
	if      (multicollinear[i_e][i_c].type_correction[j_e] == 1){temp << "QCD   ";}
	else if (multicollinear[i_e][i_c].type_correction[j_e] == 2){temp << "QEW   ";}
	temp << "no_emitter = ";
	temp << multicollinear[i_e][i_c].no_emitter[j_e] << "   ";
	temp << "type_splitting = ";
	temp << multicollinear[i_e][i_c].type_splitting[j_e] << "   ";
	/*
	temp << "charge_factor = ";
	temp << setw(15) << setprecision(8) << multicollinear[i_e][i_c].charge_factor[j_e] << " ";
	temp << "  ";
	*/
	temp << "in_collinear = ";
	temp << "(";
	for (int i_z = 0; i_z < 3; i_z++){
	  temp << multicollinear[i_e][i_c].in_collinear[j_e][i_z];
	  if (i_z < 2){temp << " ";}
	}
	temp << ")";
	logger << LOG_DEBUG << temp.str() << endl;
	temp.str("");
      }




      /*
      temp << "massive = ";
      for (int j_e = 0; j_e < i_e + 1; j_e++){
	temp << multicollinear[i_e][i_c].massive[j_e] << " ";
      }
      temp << "  ";
      */
      /*
      temp << "charge_factor = ";
      for (int j_e = 0; j_e < i_e + 1; j_e++){
	temp << setw(15) << setprecision(8) << multicollinear[i_e][i_c].charge_factor[j_e] << " ";
      }
      logger << LOG_DEBUG << temp.str() << endl;
      temp.str("");
      temp << setw(48) << "";
      temp << "type_splitting = ";
      for (int j_e = 0; j_e < i_e + 1; j_e++){
	temp << multicollinear[i_e][i_c].type_splitting[j_e] << " ";
      }
      temp << "  ";
      temp << "in_collinear = ";
      for (int j_e = 0; j_e < i_e + 1; j_e++){
	temp << "(";
	for (int i_z = 0; i_z < 3; i_z++){
	  temp << multicollinear[i_e][i_c].in_collinear[j_e][i_z];
	  if (i_z < 2){temp << " ";}
	}
	temp << " ) ";
      }
      temp << "  ";
      temp << "no_emitter = ";
      for (int j_e = 0; j_e < i_e + 1; j_e++){
	temp << multicollinear[i_e][i_c].no_emitter[j_e] << "   ";
      }
      */
      /*
      temp << "no_spectator = " << multicollinear[i_e][i_c].no_spectator() << "   ";
      temp << "pair = ";
      for (int i_p = 0; i_p < 2; i_p++){temp << multicollinear[i_e][i_c].pair()[i_p] << "   ";}
      */
      //      logger << LOG_DEBUG << temp.str() << endl;
      for (int i_x = 0; i_x < multicollinear[i_e][i_c].all_name.size(); i_x++){
	temp.str("");
	temp << setw(11) << "";
	temp << setw(1 + 9 * i_e) << left << multicollinear[i_e][i_c].all_name[i_x] << ":   " << "pdf =     (";
	for (int i_z = 0; i_z < 3; i_z++){temp << setw(4) << right << multicollinear[i_e][i_c].pdf[i_x][i_z];}
	temp << ")";
	logger << LOG_DEBUG << temp.str() << endl;
      }
      logger.newLine(LOG_DEBUG);
    }
  }
  logger.newLine(LOG_DEBUG);

  determine_CX_ncollinear_QCD(psi, n_emission);

  perform_selection_content_pdf_list();

  output_CX_ncollinear_QCD();

  /*
  CX_combination_pdf.resize(multicollinear.size());

  //  CA_combination_pdf.resize(multicollinear[n_emission].size(), vector<vector<vector<int> > > (combination_pdf.size()));
  for (int i_e = 0; i_e < multicollinear.size(); i_e++){
    CX_combination_pdf[i_e].resize(multicollinear[i_e].size(), vector<vector<vector<int> > > (combination_pdf.size()));
    for (int i_c = 0; i_c < multicollinear[i_e].size(); i_c++){


      for (int i_i = 0; i_i < CX_combination_pdf[i_e][i_c].size(); i_i++){
	// CX_combination_pdf[i_e][i_c] contains a vector, which contains the usual osi_combination_pdf(3) with 0 -> direction (1, -1), 1 -> parton with x1 (in hadron 1/2 for +1/-1), 2 -> parton with x2 (in hadron 2/1 for +1/-1)
	if ((multicollinear[i_e][i_c].pdf[i_i][1] == 10) && (multicollinear[i_e][i_c].pdf[i_i][2] == 10)){
	  for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	    if (i_q == 0){continue;}
	    for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	      if (j_q == 0){continue;}
	      vector<int> new_temp_combination_pdf = multicollinear[i_e][i_c].pdf[i_i];
	      new_temp_combination_pdf[1] = i_q;
	      new_temp_combination_pdf[2] = j_q;
	      CX_combination_pdf[i_e][i_c][i_i].push_back(new_temp_combination_pdf);
	    }
	  }
	}
	else if (multicollinear[i_e][i_c].pdf[i_i][1] == 10){
	  for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	    if (i_q == 0){continue;}
	    vector<int> new_temp_combination_pdf = multicollinear[i_e][i_c].pdf[i_i];
	    new_temp_combination_pdf[1] = i_q;
	    CX_combination_pdf[i_e][i_c][i_i].push_back(new_temp_combination_pdf);
	  }
	}
	else if (multicollinear[i_e][i_c].pdf[i_i][2] == 10){
	  for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	    if (j_q == 0){continue;}
	    vector<int> new_temp_combination_pdf = multicollinear[i_e][i_c].pdf[i_i];
	    new_temp_combination_pdf[2] = j_q;
	    CX_combination_pdf[i_e][i_c][i_i].push_back(new_temp_combination_pdf);
	  }
	}
	else {
	  CX_combination_pdf[i_e][i_c][i_i].push_back(multicollinear[i_e][i_c].pdf[i_i]);
	}
      }
    }
    logger << LOG_DEBUG << "CX_combination_pdf[i_e] determined " << endl << endl;
    
    for (int i_c = 0; i_c < multicollinear[i_e].size(); i_c++){
      for (int i_i = 0; i_i < CX_combination_pdf[i_e][i_c].size(); i_i++){
	logger << LOG_DEBUG << "pdf contributions: " << setw(15) << left << multicollinear[i_e][i_c].all_name[i_i] << ":   " << endl;
	for (int i_z = 0; i_z < CX_combination_pdf[i_e][i_c][i_i].size(); i_z++){
	  stringstream temp;
	  temp.str("");
	  temp << setw(48) << "";
	  temp << "CX_combination_pdf[i_e][" << setw(2) << i_c << "][" << setw(2) << i_i << "][" << setw(2) << i_z << "] = "; 
	  for (int i_y = 0; i_y < CX_combination_pdf[i_e][i_c][i_i][i_z].size(); i_y++){temp << setw(2) << right << CX_combination_pdf[i_e][i_c][i_i][i_z][i_y] << "  ";}
	  logger << LOG_DEBUG << temp.str() << endl;
	}
      }
    }
  } 
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::output_CX_ncollinear_QCD(){
  Logger logger("observable_set::output_CX_ncollinear_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<string> pa_name(csi->type_parton[0].size(), "");
  if (csi->type_parton[0][1] >= -10 && csi->type_parton[0][1] <= 10){pa_name[1] = "a";}
  if (csi->type_parton[0][2] >= -10 && csi->type_parton[0][2] <= 10){pa_name[2] = "b";}
  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (csi->type_parton[0][i_p] >= -10 && csi->type_parton[0][i_p] <= 10){pa_name[i_p] = char(105 + i_p - 3);}} 
  /*
  map <int,string> pname;
  fill_pname(pname);
  */
  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "ncollinear.size() = " << ncollinear.size() << endl;
  logger.newLine(LOG_DEBUG);

  for (int i_c = 0; i_c < ncollinear.size(); i_c++){
    stringstream temp;
    int n_name = ncollinear[i_c].name.size();
    
    temp << "ncollinear " << setw(2) << i_c << ":   " << setw(n_name) << left << ncollinear[i_c].name << " :";
    //<< ncollinear[i_c].n_emission[0] << ", " 
    logger << LOG_DEBUG << temp.str() << endl;
    temp.str("");
    temp << setw(n_name + 22) << "";
    temp << "csi->type_parton = (";
    for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){
      temp << setw(4) << right << ncollinear[i_c].type_parton[i_p];
      if (i_p == 2){temp << " ->";}
    }
    temp << ")   ";
    temp << "no_prc = " << ncollinear[i_c].no_prc << "   ";
    
    //    temp << "ps = " << ncollinear[i_c].previous_splitting[ncollinear[i_c].n_emission[0]] << "   ";
    logger << LOG_DEBUG << temp.str() << endl;

    vector<string> name_emission(3);
    name_emission[0] = "n_emission_leg_type    "; 
    name_emission[1] = "n_emission_leg_type_QCD"; 
    name_emission[2] = "n_emission_leg_type_QEW"; 

    for (int i_t = 0; i_t < 3; i_t++){
      temp.str("");
      temp << setw(n_name + 22) << "";
      temp << name_emission[i_t] << " = ";
      temp << "(";
      for (int i_e = 0; i_e < 3; i_e++){
	temp << ncollinear[i_c].n_emission_leg_type[i_e][i_t];
	if (i_e < 2){temp << " ";}
      }
      temp << ")";
      if (i_t == 0){
	for (int i_e = 1; i_e < 3; i_e++){
	  //	  if (ncollinear[i_c].n_emission_leg_type[i_e][i_t] == 0){continue;}
	  temp << "   " << pa_name[i_e] << " -> " << ncollinear[i_c].type_splitting_full[i_e];
	}
      }
      logger << LOG_DEBUG << temp.str() << endl;
    }

    temp.str("");
    temp << setw(n_name + 22) << "";
    temp << "emission         -> no_pdf = " << setw(2) << right << ncollinear[i_c].no_pdf;
    logger << LOG_DEBUG << temp.str() << endl;

    for (int i_e = 0; i_e < 4; i_e++){
      temp.str("");
      temp << setw(n_name + 22) << "";
      if (ncollinear[i_c].no_endpoint[i_e] == -1){
	temp << "endpoint " << i_e << " -> " << setw(2) << "--" << " -> no_pdf = "  << setw(2) << "--";
      }
      else{
	temp << "endpoint " << i_e << " -> " << setw(2) << ncollinear[i_c].no_endpoint[i_e] << " -> no_pdf = " << setw(2) << right << ncollinear[ncollinear[i_c].no_endpoint[i_e]].no_pdf;
      } 
      //      temp << "endpoint " << i_e << " -> " << ncollinear[i_c].no_endpoint[i_e];
      logger << LOG_DEBUG << temp.str() << endl;
    }


    /*
    temp.str("");
    temp << setw(n_name + 20) << "";
    temp << "n_emission_QCD = ";
    temp << "(";
    for (int i_z = 0; i_z < 3; i_z++){
      temp << ncollinear[i_c].n_emission_leg_type[i_z][1];
      if (i_z < 2){temp << " ";}
    }
    temp << ")";
    logger << LOG_DEBUG << temp.str() << endl;

    temp.str("");
    temp << setw(n_name + 20) << "";
    temp << "n_emission_QEW = ";
    temp << "(";
    for (int i_z = 0; i_z < 3; i_z++){
      temp << ncollinear[i_c].n_emission_leg_type[i_z][2];
      if (i_z < 2){temp << " ";}
    }
    temp << ")";

    for (int i_e = 1; i_e < 3; i_e++){
      if (ncollinear[i_c].n_emission[i_e] == 0){continue;}
      temp << "   " << pa_name[i_e] << " -> " << ncollinear[i_c].type_splitting_full[i_e];
    }
    logger << LOG_DEBUG << temp.str() << endl;
    */

    /*
    temp.str("");

    for (int j_e = 1; j_e < ncollinear[i_c].n_emission[0] + 1; j_e++){
      temp << setw(16 + 9 * ncollinear[i_c].n_emission[0]) << "";
      temp << "emission " << j_e << " :   ";
      if      (ncollinear[i_c].type_correction[j_e] == 1){temp << "QCD   ";}
      else if (ncollinear[i_c].type_correction[j_e] == 2){temp << "QEW   ";}
      temp << "no_emitter = ";
      temp << ncollinear[i_c].no_emitter[j_e] << "   ";
      temp << "type_splitting = ";
      temp << ncollinear[i_c].type_splitting[j_e] << "   ";
      temp << "in_collinear = ";
      temp << "(";
      for (int i_z = 0; i_z < 3; i_z++){
	temp << ncollinear[i_c].in_collinear[j_e][i_z];
	if (i_z < 2){temp << " ";}
      }
      temp << ")";
      logger << LOG_DEBUG << temp.str() << endl;
      temp.str("");
    }
    */
    
    for (int i_x = 0; i_x < ncollinear[i_c].all_name.size(); i_x++){
      temp.str("");
      temp << setw(17) << "";

      temp << setw(n_name) << left << ncollinear[i_c].all_name[i_x] << " :   " << "pdf =     (";
      for (int i_z = 0; i_z < 3; i_z++){temp << setw(4) << right << ncollinear[i_c].pdf[i_x][i_z];}
      temp << ")";
      logger << LOG_DEBUG << temp.str() << endl;
    }
    
    logger.newLine(LOG_DEBUG);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_CX_ncollinear_QCD(phasespace_set & psi, int n_emission){
  Logger logger("observable_set::determine_CX_ncollinear_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<string> pa_name(csi->type_parton[0].size(), "");
  if (csi->type_parton[0][1] >= -10 && csi->type_parton[0][1] <= 10){pa_name[1] = "a";}
  if (csi->type_parton[0][2] >= -10 && csi->type_parton[0][2] <= 10){pa_name[2] = "b";}
  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (csi->type_parton[0][i_p] >= -10 && csi->type_parton[0][i_p] <= 10){pa_name[i_p] = char(105 + i_p - 3);}} 
  /*
  map <int,string> pname;
  fill_pname(pname);
  */

  for (int i_e = 0; i_e < multicollinear.size(); i_e += n_emission){
    for (int i_c = 0; i_c < multicollinear[i_e].size(); i_c++){
      ncollinear.push_back(multicollinear[i_e][i_c]);
    }
  }

  for (int i_c = 0; i_c < ncollinear.size(); i_c++){
    // determine more general name
    stringstream temp;
    temp << "C";
    for (int i_e = 1; i_e < 3; i_e++){
      if (ncollinear[i_c].n_emission[i_e] == 0){continue;}
	//      if (ncollinear[i_c].n_emission[i_e] > 0){
      temp << "_" << pa_name[i_e];
      if (ncollinear[i_c].n_emission[i_e] > 1){temp << ncollinear[i_c].n_emission[i_e];}
      temp << "^{";
      temp << csi->name_particle[ncollinear[i_c].type_parton[i_e]] << csi->name_particle[ncollinear[i_c].pdf[0][i_e]];
      temp << "}";
	//      }
    }

    // determine more general all_name
    for (int i_a = 0; i_a < ncollinear[i_c].pdf.size(); i_a++){
      stringstream temp_all;
      temp_all << "C";
      for (int i_e = 1; i_e < 3; i_e++){
	if (ncollinear[i_c].n_emission[i_e] == 0){continue;}
	//      if (ncollinear[i_c].n_emission[i_e] > 0){
	temp_all << "_" << pa_name[i_e];
	if (ncollinear[i_c].n_emission[i_e] > 1){temp_all << ncollinear[i_c].n_emission[i_e];}
	temp_all << "^{";
	temp_all << csi->name_particle[ncollinear[0].pdf[i_a][i_e]] << csi->name_particle[ncollinear[i_c].pdf[i_a][i_e]];
	temp_all << "}";
	//      }
      }
      ncollinear[i_c].all_name[i_a] = temp_all.str();
    }

    // determine leg emission:
    if (ncollinear[i_c].n_emission_leg_type[1][0] > 0 && ncollinear[i_c].n_emission_leg_type[2][0] > 0){ncollinear[i_c].leg_emission = 0;}
    else if (ncollinear[i_c].n_emission_leg_type[1][0] > 0 && ncollinear[i_c].n_emission_leg_type[2][0] == 0){ncollinear[i_c].leg_emission = 1;}
    else if (ncollinear[i_c].n_emission_leg_type[2][0] > 0 && ncollinear[i_c].n_emission_leg_type[1][0] == 0){ncollinear[i_c].leg_emission = 2;}
    else {ncollinear[i_c].leg_emission = -1;}

    // determine numbers of ncollinear entries with corresponding endpoint pdf's:
    for (int i_e = 1; i_e < 3; i_e++){
      int j_e = i_e % 2 + 1;
      //      logger << LOG_DEBUG << "i_e = " << i_e << "   j_e = " << j_e << endl;
      if (ncollinear[i_c].n_emission_leg_type[i_e][0] == 0){continue;}
      for (int j_c = 0; j_c < ncollinear.size(); j_c++){
	//	if (i_c == j_c){continue;}
	/*
	logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << i_c << "].no_endpoint[" << i_e << "] = " << ncollinear[i_c].no_endpoint[i_e] << endl;
	logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << j_c << "].no_endpoint[" << i_e << "] = " << ncollinear[j_c].no_endpoint[i_e] << endl;
	logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << i_c << "].no_endpoint[" << j_e << "] = " << ncollinear[i_c].no_endpoint[j_e] << endl;
	logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << j_c << "].no_endpoint[" << j_e << "] = " << ncollinear[j_c].no_endpoint[j_e] << endl;
	*/
	/*
	logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << i_c << "].n_emission_leg_type[" << i_e << "][0] = " << ncollinear[i_c].n_emission_leg_type[i_e][0] << endl;
	logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << j_c << "].n_emission_leg_type[" << i_e << "][0] = " << ncollinear[j_c].n_emission_leg_type[i_e][0] << endl;
	logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << i_c << "].n_emission_leg_type[" << j_e << "][0] = " << ncollinear[i_c].n_emission_leg_type[j_e][0] << endl;
	logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << j_c << "].n_emission_leg_type[" << j_e << "][0] = " << ncollinear[j_c].n_emission_leg_type[j_e][0] << endl;
	*/
	if (ncollinear[i_c].pdf == ncollinear[j_c].pdf){
	  /*
	  if (ncollinear[i_c].n_emission_leg_type[i_e][0] > 0 && ncollinear[j_c].n_emission_leg_type[i_e][0] == 0 && 
	      (ncollinear[i_c].n_emission_leg_type[j_e][0] > 0 && ncollinear[j_c].n_emission_leg_type[j_e][0] > 0) ||
	      (ncollinear[i_c].n_emission_leg_type[j_e][0] == 0 && ncollinear[j_c].n_emission_leg_type[j_e][0] == 0)){
	    ncollinear[i_c].no_endpoint[i_e] = j_c;
	    logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << i_c << "].no_endpoint[" << i_e << "] = " << ncollinear[i_c].no_endpoint[i_e] << endl;
	  }
	  
	  if (ncollinear[j_c].n_emission_leg_type[i_e][0] == 0 && ncollinear[j_c].n_emission_leg_type[j_e][0] == 0){
	    ncollinear[i_c].no_endpoint[0] = j_c; 
	    logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << i_c << "].no_endpoint[0] = " << ncollinear[i_c].no_endpoint[0] << endl;
	  }
	  */
	  logger << LOG_DEBUG_VERBOSE << "ncollinear[" << right << setw(2) << i_c << "].emission[" << i_e << "] = " << ncollinear[i_c].emission[i_e] << endl;
	  logger << LOG_DEBUG_VERBOSE << "ncollinear[" << right << setw(2) << i_c << "].emission[" << j_e << "] = " << ncollinear[i_c].emission[j_e] << endl;
	  logger << LOG_DEBUG_VERBOSE << "ncollinear[" << right << setw(2) << j_c << "].emission[" << i_e << "] = " << ncollinear[j_c].emission[i_e] << endl;
	  logger << LOG_DEBUG_VERBOSE << "ncollinear[" << right << setw(2) << j_c << "].emission[" << j_e << "] = " << ncollinear[j_c].emission[j_e] << endl;
	  if (ncollinear[i_c].emission[i_e] == ncollinear[j_c].emission[i_e] && 
	      ncollinear[i_c].emission[j_e] == ncollinear[j_c].emission[j_e]){
	    if (ncollinear[i_c].emission[i_e] == 0 && ncollinear[i_c].emission[j_e] == 0){ncollinear[i_c].no_endpoint[0] = j_c;}
	    if (ncollinear[i_c].emission[i_e] == 1 && ncollinear[i_c].emission[j_e] == 0){ncollinear[i_c].no_endpoint[j_e] = j_c;}
	    if (ncollinear[i_c].emission[i_e] == 1 && ncollinear[i_c].emission[j_e] == 1){ncollinear[i_c].no_endpoint[3] = j_c;}
	  }
	  if (ncollinear[i_c].emission[i_e] == 1 && ncollinear[j_c].emission[i_e] == 0){
	    if (ncollinear[i_c].emission[j_e] == 0 && ncollinear[j_c].emission[j_e] == 0){ncollinear[i_c].no_endpoint[0] = j_c;}
	    if (ncollinear[i_c].emission[j_e] == 1 && ncollinear[j_c].emission[j_e] == 1){ncollinear[i_c].no_endpoint[i_e] = j_c;}
	  }
	  if (ncollinear[i_c].emission[i_e] == 1 && ncollinear[j_c].emission[i_e] == 0 &&
	      ncollinear[i_c].emission[j_e] == 1 && ncollinear[j_c].emission[j_e] == 0){
	    ncollinear[i_c].no_endpoint[0] = j_c;
	  }

	  //	  if (ncollinear[j_c].emission[i_e] == 0 && ncollinear[j_c].emission[j_e] == 0){ncollinear[i_c].no_endpoint[0] = j_c;}
	  //	  if (ncollinear[i_c].emission[i_e] == 1 && ncollinear[j_c].emission[i_e] == 1 && ncollinear[j_c].emission[j_e] == 0){ncollinear[i_c].no_endpoint[j_e] = j_c;}
	  //	  if (ncollinear[j_c].emission[i_e] == 0 && ncollinear[j_c].emission[j_e] == 0){ncollinear[i_c].no_endpoint[0] = j_c;}
	  //	  if (ncollinear[j_c].emission[i_e] == 0 && ncollinear[j_c].emission[j_e] == 0){ncollinear[i_c].no_endpoint[0] = j_c;}



	  /*
	  if (ncollinear[i_c].n_emission_leg_type[i_e][0] > 0 && ncollinear[j_c].n_emission_leg_type[i_e][0] == 0 && 
	      (ncollinear[i_c].n_emission_leg_type[j_e][0] > 0 && ncollinear[j_c].n_emission_leg_type[j_e][0] > 0) ||
	      (ncollinear[i_c].n_emission_leg_type[j_e][0] == 0 && ncollinear[j_c].n_emission_leg_type[j_e][0] == 0)){
	    ncollinear[i_c].no_endpoint[i_e] = j_c;
	    logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << i_c << "].no_endpoint[" << i_e << "] = " << ncollinear[i_c].no_endpoint[i_e] << endl;
	  }
	  if (ncollinear[i_c].n_emission_leg_type[i_e][0] > 0 && ncollinear[j_c].n_emission_leg_type[i_e][0] == 0 && 
	      ncollinear[i_c].n_emission_leg_type[j_e][0] > 0 && ncollinear[j_c].n_emission_leg_type[j_e][0] == 0){
	    ncollinear[i_c].no_endpoint[0] = j_c; 
	    logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << i_c << "].no_endpoint[0] = " << ncollinear[i_c].no_endpoint[0] << endl;
	  }
	  */
	}
      }
    }

    // determine type_splitting_full (for each leg):
    // 14 -  0   +   0 - 14
    // 34 -  0   +   0 - 34
    // 54 -  0   +   0 - 54
    // 22 -  0   +   0 - 22
    //  2 -  2
    //  2 -  4   +   4 -  2
    //  4 -  4

    for (int i_e = 1; i_e < 3; i_e++){
      if (ncollinear[i_c].n_emission[i_e] == 0){continue;}
      for (int i_s = 0; i_s < ncollinear[i_c].type_splitting_leg[i_e].size(); i_s++){
	ncollinear[i_c].type_splitting_full[i_e] += pow(10., i_s) * (ncollinear[i_c].type_splitting_leg[i_e][i_s]);
      }
    }


    // determine order_emission
    if (ncollinear[i_c].leg_emission == 0){
      if (ncollinear[i_c].type_splitting_full[2] < ncollinear[i_c].type_splitting_full[1]){
	ncollinear[i_c].x_a = 2;
	ncollinear[i_c].x_b = 1;
      }
      else {
	ncollinear[i_c].x_a = 1;
	ncollinear[i_c].x_b = 2;
      }
    }
    else if (ncollinear[i_c].leg_emission > 0){
      ncollinear[i_c].x_a = ncollinear[i_c].leg_emission;
      ncollinear[i_c].x_b = ncollinear[i_c].x_a % 2 + 1;
    }



    logger << LOG_DEBUG << "ncollinear[" << right << setw(2) << i_c << "].name =    old: " << left << setw(20) << ncollinear[i_c].name << "    new: " << setw(20) << temp.str() << endl;
    ncollinear[i_c].name = temp.str();

    logger << LOG_DEBUG << left << setw(20) << ncollinear[i_c].name <<    "type_splitting_full[1] = " << setw(2) << ncollinear[i_c].type_splitting_full[1] << "   type_splitting_full[2] = " << setw(2) << ncollinear[i_c].type_splitting_full[2];
    logger.newLine(LOG_DEBUG);

    //      logger << LOG_DEBUG << "new name:   ncollinear[" << setw(2) << i_c << "].name = " << ncollinear[i_c].name << endl;
  }




  /*
  for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
    for (int i_z = 0; i_z < 3; i_z++){
      if ((*CA_collinear)[i_a][0].in_collinear()[i_z] == 1){CA_dipole_splitting[(*CA_collinear)[i_a][0].type()][i_z] = 1;}
    }
  }

  logger << LOG_DEBUG << "CA_dipole_splitting:" << endl;
  for (int i_z = 0; i_z < 3; i_z++){
    stringstream temp_ss;
    temp_ss << "i_z = " << i_z << ":   ";
    for (int i_t = 0; i_t < 4; i_t++){temp_ss << CA_dipole_splitting[i_t][i_z] << "   ";}
    logger << LOG_DEBUG << temp_ss.str() << endl;
  }
  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "new collinear dipoles determined " << endl << endl;
  */
  



  // intermediate step of CX_combination_pdf most likely not needed !!!

  /*
  CX_combination_pdf.resize(1, vector<vector<vector<vector<int> > > > (ncollinear.size(), vector<vector<vector<int> > > (combination_pdf.size())));

  for (int i_c = 0; i_c < ncollinear.size(); i_c++){
    for (int i_i = 0; i_i < CX_combination_pdf[0][i_c].size(); i_i++){
      // CX_combination_pdf[0][i_c] contains a vector, which contains the usual osi_combination_pdf(3) with 0 -> direction (1, -1), 1 -> parton with x1 (in hadron 1/2 for +1/-1), 2 -> parton with x2 (in hadron 2/1 for +1/-1)
      if ((ncollinear[i_c].pdf[i_i][1] == 10) && (ncollinear[i_c].pdf[i_i][2] == 10)){
	for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	  if (i_q == 0){continue;}
	  for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	    if (j_q == 0){continue;}
	    vector<int> new_temp_combination_pdf = ncollinear[i_c].pdf[i_i];
	    new_temp_combination_pdf[1] = i_q;
	    new_temp_combination_pdf[2] = j_q;
	    CX_combination_pdf[0][i_c][i_i].push_back(new_temp_combination_pdf);
	  }
	}
      }
      else if (ncollinear[i_c].pdf[i_i][1] == 10){
	for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	  if (i_q == 0){continue;}
	  vector<int> new_temp_combination_pdf = ncollinear[i_c].pdf[i_i];
	  new_temp_combination_pdf[1] = i_q;
	  CX_combination_pdf[0][i_c][i_i].push_back(new_temp_combination_pdf);
	}
      }
      else if (ncollinear[i_c].pdf[i_i][2] == 10){
	for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	  if (j_q == 0){continue;}
	  vector<int> new_temp_combination_pdf = ncollinear[i_c].pdf[i_i];
	  new_temp_combination_pdf[2] = j_q;
	  CX_combination_pdf[0][i_c][i_i].push_back(new_temp_combination_pdf);
	}
      }
      else {
	CX_combination_pdf[0][i_c][i_i].push_back(ncollinear[i_c].pdf[i_i]);
      }
    }
  }
  logger << LOG_DEBUG << "CX_combination_pdf[0] determined " << endl << endl;
  
  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "pdf contributions:   CX_combination_pdf" << endl;
  logger.newLine(LOG_DEBUG);
  for (int i_c = 0; i_c < ncollinear.size(); i_c++){
    for (int i_i = 0; i_i < CX_combination_pdf[0][i_c].size(); i_i++){
      logger << LOG_DEBUG << "pdf contributions: " << setw(15) << left << ncollinear[i_c].all_name[i_i] << ":   " << endl;
      for (int i_z = 0; i_z < CX_combination_pdf[0][i_c][i_i].size(); i_z++){
	stringstream temp;
	temp.str("");
	temp << setw(48) << "";
	temp << "CX_combination_pdf[0][" << setw(2) << i_c << "][" << setw(2) << i_i << "][" << setw(2) << i_z << "] = "; 
	for (int i_y = 0; i_y < CX_combination_pdf[0][i_c][i_i][i_z].size(); i_y++){temp << setw(2) << right << CX_combination_pdf[0][i_c][i_i][i_z][i_y] << "  ";}
	logger << LOG_DEBUG << temp.str() << endl;
      }
    }
  }

  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "CX_combination_pdf   determined." << endl;
  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "Apply 'perform_selection_pdf_content' on CX_combination_pdf" << endl;
  logger.newLine(LOG_DEBUG);
 



  // determination of list_combination_pdf (including list_combination_pdf_emission) from CX_combination_pdf:

  for (int i_c = 0; i_c < CX_combination_pdf[0].size(); i_c++){
    int flag = list_combination_pdf.size();
    for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
      if (CX_combination_pdf[0][i_c] == list_combination_pdf[i_l] && ncollinear[i_c].emission == list_combination_pdf_emission[i_l]){
	flag = i_l;
	break;
      }
    }
    ncollinear[i_c].no_pdf = flag;
    if (flag == list_combination_pdf.size()){
      list_combination_pdf.push_back(CX_combination_pdf[0][i_c]);
      list_combination_pdf_emission.push_back(ncollinear[i_c].emission);
    }
  }

*/



  // direct determination of list_combination_pdf from ncollinear (including list_combination_pdf_emission):

  list_combination_pdf.resize(ncollinear.size(), vector<vector<vector<int> > > (combination_pdf.size()));
  list_combination_pdf_emission.resize(ncollinear.size());

  for (int i_c = 0; i_c < ncollinear.size(); i_c++){
    list_combination_pdf_emission[i_c] = ncollinear[i_c].emission;
    ncollinear[i_c].no_pdf = i_c;
    for (int i_i = 0; i_i < list_combination_pdf[i_c].size(); i_i++){
      // list_combination_pdf[i_c] contains a vector, which contains the usual osi_combination_pdf(3) with 0 -> direction (1, -1), 1 -> parton with x1 (in hadron 1/2 for +1/-1), 2 -> parton with x2 (in hadron 2/1 for +1/-1)
      if ((ncollinear[i_c].pdf[i_i][1] == 10) && (ncollinear[i_c].pdf[i_i][2] == 10)){
	for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	  if (i_q == 0){continue;}
	  for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	    if (j_q == 0){continue;}
	    vector<int> new_temp_combination_pdf = ncollinear[i_c].pdf[i_i];
	    new_temp_combination_pdf[1] = i_q;
	    new_temp_combination_pdf[2] = j_q;
	    list_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	  }
	}
      }
      else if (ncollinear[i_c].pdf[i_i][1] == 10){
	for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	  if (i_q == 0){continue;}
	  vector<int> new_temp_combination_pdf = ncollinear[i_c].pdf[i_i];
	  new_temp_combination_pdf[1] = i_q;
	  list_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	}
      }
      else if (ncollinear[i_c].pdf[i_i][2] == 10){
	for (int j_q = -N_f_active; j_q < N_f_active + 1; j_q++){
	  if (j_q == 0){continue;}
	  vector<int> new_temp_combination_pdf = ncollinear[i_c].pdf[i_i];
	  new_temp_combination_pdf[2] = j_q;
	  list_combination_pdf[i_c][i_i].push_back(new_temp_combination_pdf);
	}
      }
      else {
	list_combination_pdf[i_c][i_i].push_back(ncollinear[i_c].pdf[i_i]);
      }
    }
  }
  logger << LOG_DEBUG << "list_combination_pdf determined " << endl << endl;






  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "pdf contributions:   list_combination_pdf, list_combination_pdf_emission" << endl;
  logger.newLine(LOG_DEBUG);
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    stringstream temp_em;
    for (int i_e = 1; i_e < 3; i_e++){
      temp_em << "   ";
      if (list_combination_pdf_emission[i_l][i_e] == 1){temp_em << "x_" << i_e << " / z_" << i_e;}
      else if (list_combination_pdf_emission[i_l][i_e] == 0){temp_em << "x_" << i_e;}
      else {logger << LOG_FATAL << temp_em.str() << " is not a valid value." << endl; exit(1);}
    }
    logger << LOG_DEBUG << "pdf contributions [" << setw(2) << i_l << "]   @" << temp_em.str() << endl;// << setw(15) << left << ncollinear[i_c].all_name[i_i] << ":   " << endl;

    for (int i_i = 0; i_i < list_combination_pdf[i_l].size(); i_i++){
       for (int i_z = 0; i_z < list_combination_pdf[i_l][i_i].size(); i_z++){
	stringstream temp;
	temp.str("");
	temp << setw(32) << "";
	temp << "list_combination_pdf[" << setw(2) << i_l << "][" << setw(2) << i_i << "][" << setw(2) << i_z << "] = "; 
	for (int i_y = 0; i_y < list_combination_pdf[i_l][i_i][i_z].size(); i_y++){temp << setw(2) << right << list_combination_pdf[i_l][i_i][i_z][i_y] << "  ";}
	logger << LOG_DEBUG << temp.str() << endl;
      }
    }
  }
  logger.newLine(LOG_DEBUG);

  //  output_CX_ncollinear_QCD();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




  // xz_pdf: 0 - 0 : empty - could be filled with tau = x1 * x2
  //         1 - 0 : x1 -> no collinear emission from 1st parton
  //         2 - 0 : x2 -> no collinear emission from 2nd parton
  //         0 - 1 : empty
  //         1 - 1 : x1 / z1 -> collinear emission from 1st parton
  //         2 - 1 : x2 / z2 -> collinear emission from 2nd parton

void observable_set::calculate_IS_CX(phasespace_set & psi){
  Logger logger("observable_set::calculate_IS_CX");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  logger << LOG_DEBUG_VERBOSE << "psi_x_pdf.size()  = " << psi_x_pdf.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "psi_z_coll.size() = " << psi_z_coll.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "x_pdf.size()  = " << x_pdf.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "z_coll.size() = " << z_coll.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "xz_pdf.size() = " << xz_pdf.size() << endl;
  for (int i_x = 0; i_x < 3; i_x++){
    logger << LOG_DEBUG_VERBOSE << "xz_pdf[" << i_x << "].size() = " << xz_pdf[i_x].size() << endl;
  }

  for (int i_x = 0; i_x < 3; i_x++){xz_pdf[i_x][0] = psi_x_pdf[i_x];}
  for (int i_x = 1; i_x < 3; i_x++){xz_pdf[i_x][1] = psi_x_pdf[i_x] / psi_z_coll[i_x];}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

