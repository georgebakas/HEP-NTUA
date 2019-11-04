#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
void QEW_determine_dipoles(vector<dipole_set> & QEW_dipole_candidate, vector<int> & type_parton, vector<int> & basic_type_parton){
  Logger logger("QEW_determine_dipoles");
  static int initialization = 1;
  static map<int, double> charge_particle;
  if (initialization == 1){
    fill_charge_particle(charge_particle);
    initialization = 0;
  }
  vector<string> pa_name(type_parton.size());
  if (charge_particle[type_parton[1]] != 0 || type_parton[1] == 22){pa_name[1] = "a";}
  if (charge_particle[type_parton[2]] != 0 || type_parton[2] == 22){pa_name[2] = "b";}
  //  pa_name[1] = "a";
  //  pa_name[2] = "b";
  /*
  int count = 0;
  vector<string> alphabet(type_parton.size() - 3, "");
  for (int i_p = 0; i_p < alphabet.size(); i_p++){alphabet[i_p] = char(105 + i_p);}
  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (charge_particle[type_parton[i_p]] != 0 || type_parton[i_p] == 22){pa_name[i_p] = alphabet[count++];}}
  */
  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (charge_particle[type_parton[i_p]] != 0 || type_parton[i_p] == 22){pa_name[i_p] = char(105 + i_p - 3);}}
  //  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (charge_particle[csi->type_parton[0][i_p]] != 0 || csi->type_parton[0][i_p] == 22){pa_name[i_p] = char(105 + i_p - 3);}}

  logger << LOG_DEBUG << "EW dipole determination " << endl << endl;
  int type_splitting;
  string temp_name;
  int no_A_emitter;
  int no_A_spectator;
  int type_correction = 2;

  logger << LOG_DEBUG << "final-state emitter, final-state spectator" << endl;
  for (int no_R_emitter_1 = 3; no_R_emitter_1 < pa_name.size(); no_R_emitter_1++){
    if (pa_name[no_R_emitter_1] == ""){continue;}
    for (int no_R_emitter_2 = 3; no_R_emitter_2 < pa_name.size(); no_R_emitter_2++){
      if (no_R_emitter_1 == no_R_emitter_2){continue;}
      if (pa_name[no_R_emitter_2] == ""){continue;}
      if      (type_parton[no_R_emitter_1] != 22 && type_parton[no_R_emitter_2] == 22){type_splitting = 2;}
      else if (type_parton[no_R_emitter_1] >   0 && type_parton[no_R_emitter_2] == -type_parton[no_R_emitter_1]){type_splitting = 1;}
      else {continue;}
      for (int no_R_spectator = 3; no_R_spectator < pa_name.size(); no_R_spectator++){
	if (pa_name[no_R_spectator] == ""){continue;}
	if (no_R_spectator == no_R_emitter_1){continue;}
	if (no_R_spectator == no_R_emitter_2){continue;}
	cout << "type_splitting = " << type_splitting << endl;
	cout << "pa_name[no_R_emitter_1 = " << no_R_emitter_1 << "] = " << pa_name[no_R_emitter_1] << "   type_parton[" << no_R_emitter_1 << "] = " << type_parton[no_R_emitter_1] << endl;
	cout << "pa_name[no_R_emitter_2 = " << no_R_emitter_2 << "] = " << pa_name[no_R_emitter_2] << "   type_parton[" << no_R_emitter_2 << "] = " << type_parton[no_R_emitter_2] << endl;
	cout << "pa_name[no_R_spectator = " << no_R_spectator << "] = " << pa_name[no_R_spectator] << "   type_parton[" << no_R_spectator << "] = " << type_parton[no_R_spectator] << endl;
	temp_name = "D_{" + pa_name[no_R_emitter_1] + pa_name[no_R_emitter_2] + "," + pa_name[no_R_spectator] + "}";
	vector<int> temp_type_parton = type_parton;
	if (type_splitting == 1){temp_type_parton[no_R_emitter_1] = 22;}
	temp_type_parton.erase(temp_type_parton.begin() + no_R_emitter_2);
	vector<int> temp_basic_type_parton = basic_type_parton;
	if (type_splitting == 1){temp_basic_type_parton[no_R_emitter_1] = 22;}
	temp_basic_type_parton.erase(temp_basic_type_parton.begin() + no_R_emitter_2);
	if (no_R_emitter_2 > no_R_emitter_1){no_A_emitter = no_R_emitter_1;}
	else {no_A_emitter = no_R_emitter_1 - 1;}
	if (no_R_emitter_2 > no_R_spectator){no_A_spectator = no_R_spectator;}
	else {no_A_spectator = no_R_spectator - 1;}
	//	else {no_A_emitter = no_R_spectator - 1;}
	cout << "pa_name[no_A_emitter   = " << no_A_emitter << "] = " << pa_name[no_A_emitter] << "   temp_type_parton[" << no_A_emitter << "] = " << temp_type_parton[no_A_emitter] << endl;
	cout << "pa_name[no_A_spectator = " << no_A_spectator << "] = " << pa_name[no_A_spectator] << "   temp_type_parton[" << no_A_spectator << "] = " << temp_type_parton[no_A_spectator] << endl;
	logger << LOG_DEBUG << temp_name << "   pa[R,e1:" << setw(2) << no_R_emitter_1 << "] = " << setw(3) << type_parton[no_R_emitter_1] << "   pa[R,e2:" << setw(2) << no_R_emitter_2 << "] = " << setw(3) << type_parton[no_R_emitter_2] << "   pa[R,sp:" << setw(2) << no_R_spectator << "] = " << setw(3) << type_parton[no_R_spectator] << "      pa[A,e:" << setw(2) << no_A_emitter << "] = " << setw(3) << temp_type_parton[no_A_emitter] << "   pa[A,sp:" << setw(2) << no_A_spectator << "] = " << setw(3) << temp_type_parton[no_A_spectator] << endl;

	QEW_dipole_candidate.push_back(dipole_set(temp_name, temp_type_parton, temp_basic_type_parton, 1, type_splitting, no_R_emitter_1, no_R_emitter_2, no_R_spectator, no_A_emitter, no_A_spectator, type_correction));
      }
    }
  }

  logger << LOG_DEBUG << "final-state emitter, initial-state spectator" << endl;
  for (int no_R_emitter_1 = 3; no_R_emitter_1 < pa_name.size(); no_R_emitter_1++){
    if (pa_name[no_R_emitter_1] == ""){continue;}
    for (int no_R_emitter_2 = 3; no_R_emitter_2 < pa_name.size(); no_R_emitter_2++){
      if (no_R_emitter_1 == no_R_emitter_2){continue;}
      if (pa_name[no_R_emitter_2] == ""){continue;}
      if      (type_parton[no_R_emitter_1] != 22 && type_parton[no_R_emitter_2] == 22){type_splitting = 2;}
      else if (type_parton[no_R_emitter_1] >   0 && type_parton[no_R_emitter_2] == -type_parton[no_R_emitter_1]){type_splitting = 1;}
      else {continue;}
      for (int no_R_spectator = 1; no_R_spectator < 3; no_R_spectator++){
	if (pa_name[no_R_spectator] == ""){continue;}
	temp_name = "D^" + pa_name[no_R_spectator] + "_{" + pa_name[no_R_emitter_1] + pa_name[no_R_emitter_2] + "}";
	vector<int> temp_type_parton = type_parton;
	if (type_splitting == 1){temp_type_parton[no_R_emitter_1] = 22;}
	temp_type_parton.erase(temp_type_parton.begin() + no_R_emitter_2);
	vector<int> temp_basic_type_parton = basic_type_parton;
	if (type_splitting == 1){temp_basic_type_parton[no_R_emitter_1] = 22;}
	temp_basic_type_parton.erase(temp_basic_type_parton.begin() + no_R_emitter_2);
	if (no_R_emitter_2 > no_R_emitter_1){no_A_emitter = no_R_emitter_1;}
	else {no_A_emitter = no_R_emitter_1 - 1;}
	no_A_spectator = no_R_spectator;
	logger << LOG_DEBUG << temp_name << "   pa[R,e1:" << setw(2) << no_R_emitter_1 << "] = " << setw(3) << type_parton[no_R_emitter_1] << "   pa[R,e2:" << setw(2) << no_R_emitter_2 << "] = " << setw(3) << type_parton[no_R_emitter_2] << "   pa[R,sp:" << setw(2) << no_R_spectator << "] = " << setw(3) << type_parton[no_R_spectator] << "      pa[A,e:" << setw(2) << no_A_emitter << "] = " << setw(3) << temp_type_parton[no_A_emitter] << "   pa[A,sp:" << setw(2) << no_A_spectator << "] = " << setw(3) << temp_type_parton[no_A_spectator] << endl;

	QEW_dipole_candidate.push_back(dipole_set(temp_name, temp_type_parton, temp_basic_type_parton, 2, type_splitting, no_R_emitter_1, no_R_emitter_2, no_R_spectator, no_A_emitter, no_A_spectator, type_correction));
      }
    }
  }
  
  logger << LOG_DEBUG << "initial-state emitter, final-state spectator" << endl;
  for (int no_R_emitter_1 = 1; no_R_emitter_1 < 3; no_R_emitter_1++){
    if (pa_name[no_R_emitter_1] == ""){continue;}
    for (int no_R_emitter_2 = 3; no_R_emitter_2 < pa_name.size(); no_R_emitter_2++){
      if (pa_name[no_R_emitter_2] == ""){continue;}
      if      (type_parton[no_R_emitter_1] != 22 && type_parton[no_R_emitter_2] == 22){type_splitting = 2;}
      else if (type_parton[no_R_emitter_1] == 22 && type_parton[no_R_emitter_2] != 22){type_splitting = 3;}
      else if (type_parton[no_R_emitter_1] != 22 && type_parton[no_R_emitter_2] == type_parton[no_R_emitter_1]){type_splitting = 1;}
      else {continue;}
      for (int no_R_spectator = 3; no_R_spectator < pa_name.size(); no_R_spectator++){
	if (pa_name[no_R_spectator] == ""){continue;}
	if (no_R_spectator == no_R_emitter_2){continue;}
	temp_name = "D^{" + pa_name[no_R_emitter_1] + pa_name[no_R_emitter_2] + "}_" + pa_name[no_R_spectator];
	vector<int> temp_type_parton = type_parton;
	if      (type_splitting == 3){temp_type_parton[no_R_emitter_1] = -type_parton[no_R_emitter_2];}
	else if (type_splitting == 1){temp_type_parton[no_R_emitter_1] = 22;}
	temp_type_parton.erase(temp_type_parton.begin() + no_R_emitter_2);
	vector<int> temp_basic_type_parton = basic_type_parton;
	if      (type_splitting == 3){temp_basic_type_parton[no_R_emitter_1] = -basic_type_parton[no_R_emitter_2];}
	else if (type_splitting == 1){temp_basic_type_parton[no_R_emitter_1] = 22;}
	temp_basic_type_parton.erase(temp_basic_type_parton.begin() + no_R_emitter_2);
	no_A_emitter = no_R_emitter_1;
	if (no_R_spectator < no_R_emitter_2){no_A_spectator = no_R_spectator;}
	else {no_A_spectator = no_R_spectator - 1;}
	logger << LOG_DEBUG << temp_name << "   pa[R,e1:" << setw(2) << no_R_emitter_1 << "] = " << setw(3) << type_parton[no_R_emitter_1] << "   pa[R,e2:" << setw(2) << no_R_emitter_2 << "] = " << setw(3) << type_parton[no_R_emitter_2] << "   pa[R,sp:" << setw(2) << no_R_spectator << "] = " << setw(3) << type_parton[no_R_spectator] << "      pa[A,e:" << setw(2) << no_A_emitter << "] = " << setw(3) << temp_type_parton[no_A_emitter] << "   pa[A,sp:" << setw(2) << no_A_spectator << "] = " << setw(3) << temp_type_parton[no_A_spectator] << endl;

	QEW_dipole_candidate.push_back(dipole_set(temp_name, temp_type_parton, temp_basic_type_parton, 3, type_splitting, no_R_emitter_1, no_R_emitter_2, no_R_spectator, no_A_emitter, no_A_spectator, type_correction));
     }
    }
  }

  logger << LOG_DEBUG << "initial-state emitter, initial-state spectator" << endl;
  for (int no_R_emitter_1 = 1; no_R_emitter_1 < 3; no_R_emitter_1++){
    if (pa_name[no_R_emitter_1] == ""){continue;}
    for (int no_R_emitter_2 = 3; no_R_emitter_2 < pa_name.size(); no_R_emitter_2++){
      if (pa_name[no_R_emitter_2] == ""){continue;}
      if      (type_parton[no_R_emitter_1] != 22 && type_parton[no_R_emitter_2] == 22){type_splitting = 2;}
      else if (type_parton[no_R_emitter_1] == 22 && type_parton[no_R_emitter_2] != 22){type_splitting = 3;}
      else if (type_parton[no_R_emitter_1] != 22 && type_parton[no_R_emitter_2] == type_parton[no_R_emitter_1]){type_splitting = 1;}
      else {continue;}
      for (int no_R_spectator = 1; no_R_spectator < 3; no_R_spectator++){
	if (pa_name[no_R_spectator] == ""){continue;}
	if (no_R_spectator == no_R_emitter_1){continue;}
	temp_name = "D^{" + pa_name[no_R_emitter_1] + pa_name[no_R_emitter_2] + "," + pa_name[no_R_spectator] + "}";
	vector<int> temp_type_parton = type_parton;
	if      (type_splitting == 3){temp_type_parton[no_R_emitter_1] = -type_parton[no_R_emitter_2];}
	else if (type_splitting == 1){temp_type_parton[no_R_emitter_1] = 22;}
	temp_type_parton.erase(temp_type_parton.begin() + no_R_emitter_2);
	vector<int> temp_basic_type_parton = basic_type_parton;
	if      (type_splitting == 3){temp_basic_type_parton[no_R_emitter_1] = -basic_type_parton[no_R_emitter_2];}
	else if (type_splitting == 1){temp_basic_type_parton[no_R_emitter_1] = 22;}
	temp_basic_type_parton.erase(temp_basic_type_parton.begin() + no_R_emitter_2);
	no_A_emitter = no_R_emitter_1;
	no_A_spectator = no_R_spectator;
	logger << LOG_DEBUG << temp_name << "   pa[R,e1:" << setw(2) << no_R_emitter_1 << "] = " << setw(3) << type_parton[no_R_emitter_1] << "   pa[R,e2:" << setw(2) << no_R_emitter_2 << "] = " << setw(3) << type_parton[no_R_emitter_2] << "   pa[R,sp:" << setw(2) << no_R_spectator << "] = " << setw(3) << type_parton[no_R_spectator] << "      pa[A,e:" << setw(2) << no_A_emitter << "] = " << setw(3) << temp_type_parton[no_A_emitter] << "   pa[A,sp:" << setw(2) << no_A_spectator << "] = " << setw(3) << temp_type_parton[no_A_spectator] << endl;

	QEW_dipole_candidate.push_back(dipole_set(temp_name, temp_type_parton, temp_basic_type_parton, 5, type_splitting, no_R_emitter_1, no_R_emitter_2, no_R_spectator, no_A_emitter, no_A_spectator, type_correction));
     }
    }
  }

  logger << LOG_DEBUG << "QEW_dipole_candidate.size() = " << QEW_dipole_candidate.size() << endl;
  for (int i_d = 0; i_d < QEW_dipole_candidate.size(); i_d++){
    logger << LOG_DEBUG << QEW_dipole_candidate[i_d].name() << endl;
  }
  logger << LOG_DEBUG << endl << "QEW dipoles determined " << endl << endl;
}



void QEW_selection_dipoles(vector<dipole_set> & dipole, vector<dipole_set> & dipole_candidate, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference, vector<vector<double> > & singular_region, vector<vector<string> > & singular_region_name, vector<vector<int> > & list_singular_regions, phasespace_set & psi, call_generic & generic){
  static Logger logger("QEW_selection_dipoles");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  static int initialization = 1;
  static vector<int> type_parton = dipole[0].type_parton();
  static map<int, double> charge_particle;
  if (initialization == 1){
    fill_charge_particle(charge_particle);
    initialization = 0;
  }
  singular_region.resize(type_parton.size(), vector<double> (type_parton.size()));
  singular_region_name.resize(type_parton.size(), vector<string> (type_parton.size()));
  vector<string> pa_name(type_parton.size(), "");
  if (charge_particle[type_parton[1]] != 0 || type_parton[1] == 22){pa_name[1] = "a";}
  if (charge_particle[type_parton[2]] != 0 || type_parton[2] == 22){pa_name[2] = "b";}
  int count = 0;
  vector<string> alphabet(type_parton.size() - 3, "");
  for (int i_p = 0; i_p < alphabet.size(); i_p++){alphabet[i_p] = char(105 + i_p);}
  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (charge_particle[type_parton[i_p]] != 0 || type_parton[i_p] == 22){pa_name[i_p] = alphabet[count++];}}
  logger << LOG_DEBUG << "QEW dipole selection from QEW dipole candidates started." << endl;
  //  int x_a = 0; // wrong if QCD and QEW dipoles are mixed !!!
  int x_a = dipole.size() - 1;
  for (int i_a = 0; i_a < dipole_candidate.size(); i_a++){
    vector<int> temp_o_map(dipole_candidate[i_a].type_parton().size());
    int temp_no_prc;
    vector<int> temp_o_prc;
    vector<int> temp_type_parton = dipole_candidate[i_a].type_parton();
    vector<int> temp_basic_type_parton = dipole_candidate[i_a].basic_type_parton();
    double temp_symmetry_factor;
    int temp_no_map;
    generic.determination_no_subprocess_dipole(temp_no_map, temp_o_map, temp_no_prc, temp_o_prc, temp_symmetry_factor, temp_basic_type_parton, basic_order_alpha_s, basic_order_alpha_e, basic_order_interference);
    //    int temp_no_map = generic.determination_no_subprocess_dipole(temp_o_map, temp_no_prc, temp_o_prc, temp_symmetry_factor, temp_type_parton, basic_order_alpha_s, basic_order_alpha_e, basic_order_interference);
    ///    int temp_no_map = generic.determination_no_subprocess_dipole(temp_o_map, temp_type_parton, basic_order_alpha_s, basic_order_alpha_e, basic_order_interference);
    if (temp_no_map == -1){continue;}
    double temp_charge_factor = 0.;
    if (dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_emitter()] != 22 && dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_spectator()] != 22){
      temp_charge_factor = charge_particle[abs(dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_emitter()])] * charge_particle[abs(dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_spectator()])];
      if ((dipole_candidate[i_a].no_A_emitter() > 2 && dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_emitter()] > 0) || (dipole_candidate[i_a].no_A_emitter() < 3 && dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_emitter()] < 0)){temp_charge_factor = -temp_charge_factor;}
      if ((dipole_candidate[i_a].no_A_spectator() > 2 && dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_spectator()] > 0) || (dipole_candidate[i_a].no_A_spectator() < 3 && dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_spectator()] < 0)){temp_charge_factor = -temp_charge_factor;}
    }
    else if (dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_emitter()] == 22){
      // charge_factor independent of spectator !?!
      temp_charge_factor = charge_particle[abs(type_parton[dipole_candidate[i_a].no_R_emitter_1()])] * charge_particle[abs(type_parton[dipole_candidate[i_a].no_R_emitter_2()])];
      if ((dipole_candidate[i_a].no_R_emitter_1() > 2 && type_parton[dipole_candidate[i_a].no_R_emitter_1()] > 0) || (dipole_candidate[i_a].no_R_emitter_1() < 3 && type_parton[dipole_candidate[i_a].no_R_emitter_1()] < 0)){temp_charge_factor = -temp_charge_factor;}
      if ((dipole_candidate[i_a].no_R_emitter_2() > 2 && type_parton[dipole_candidate[i_a].no_R_emitter_2()] > 0) || (dipole_candidate[i_a].no_R_emitter_2() < 3 && type_parton[dipole_candidate[i_a].no_R_emitter_2()] < 0)){temp_charge_factor = -temp_charge_factor;}
      // initial-state photon as emitter:
      // selection of kappa_ij,k such that kappa_ij,k = -1 for the other initial state: does not work for ga initial states so far !!!
      if (dipole_candidate[i_a].no_R_emitter_1() < 3){
	if (dipole_candidate[i_a].no_R_spectator() > 2){temp_charge_factor = 0.;}
      }
      cout << "temp_charge_factor = " << temp_charge_factor << endl;
    }
    else {
      temp_charge_factor = 0.;
    }
    if (temp_charge_factor == 0.){
      logger << LOG_DEBUG << "dipole_candidate[" << i_a << "] is not used because of spectator selection!" << endl;
      continue;
    }
    x_a++; // labels the present dipole
    /*
    vector<vector<su3generator> > T;
    vector<vector<su3structure> > f;
    vector<vector<su3delta> > delta;
    vector<vector<vector<double> > > dipole_colourmatrix(dipole_candidate.size());
    vector<vector<vector<int> > > dipole_spinorder(dipole_candidate.size());
    vector<vector<int> > dipole_fckm(dipole_candidate.size());
    vector<vector<int> > dipole_data(dipole_candidate.size(), vector<int> (3));
*/
    // could be replaced... (map == prc)
    ///    int temp_no_prc;
    ///    vector<int> temp_o_prc(dipole_candidate[i_a].type_parton().size());
    ///    generic.initialization_subprocess_dipole(temp_type_parton, temp_o_prc, temp_o_map, T, f, delta, dipole_colourmatrix[i_a], dipole_spinorder[i_a], dipole_fckm[i_a], dipole_data[i_a], temp_no_prc, temp_no_map, basic_order_alpha_s, basic_order_alpha_e, basic_order_interference);
    psi_no_map.push_back(temp_no_map);
    psi_o_map.push_back(temp_o_map);
    psi_no_prc.push_back(temp_no_prc);
    psi_o_prc.push_back(temp_o_prc);
    psi_phasespace_order_alpha_s.push_back(psi_phasespace_order_alpha_s[0]);
    psi_phasespace_order_alpha_e.push_back(psi_phasespace_order_alpha_e[0] - 1);
    psi_phasespace_order_interference.push_back(psi_phasespace_order_interference[0]);
    psi_MC_n_channel_phasespace.push_back(generic.determination_MCchannels_dipole(x_a, psi));
    psi_MC_sum_channel_phasespace.push_back(psi_MC_sum_channel_phasespace[x_a - 1] + psi_MC_n_channel_phasespace[x_a]);
    logger << LOG_DEBUG_VERBOSE << "psi_MC_n_channel_phasespace[x_a = " << x_a << "] = " << psi_MC_n_channel_phasespace[x_a] << endl;
    logger << LOG_DEBUG_VERBOSE << "psi_MC_sum_channel_phasespace[x_a = " << x_a << "] = " << psi_MC_n_channel_phasespace[x_a] << endl;
    //    int temp_n_channel = generic.determination_MCchannels_dipole(temp_no_map, basic_order_alpha_s, basic_order_alpha_e, basic_order_interference);
    if (temp_no_prc == 0){continue;}
    /*
    double temp_charge_factor = 0.;
    if (dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_emitter()] != 22 && dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_spectator()] != 22){
      temp_charge_factor = charge_particle[abs(dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_emitter()])] * charge_particle[abs(dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_spectator()])];
      if ((dipole_candidate[i_a].no_A_emitter() > 2 && dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_emitter()] > 0) || (dipole_candidate[i_a].no_A_emitter() < 3 && dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_emitter()] < 0)){temp_charge_factor = -temp_charge_factor;}
      if ((dipole_candidate[i_a].no_A_spectator() > 2 && dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_spectator()] > 0) || (dipole_candidate[i_a].no_A_spectator() < 3 && dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_spectator()] < 0)){temp_charge_factor = -temp_charge_factor;}
    }
    else if (dipole_candidate[i_a].type_parton()[dipole_candidate[i_a].no_A_emitter()] == 22){
      // charge_factor independent of spectator !?!
      temp_charge_factor = charge_particle[abs(type_parton[dipole_candidate[i_a].no_R_emitter_1()])] * charge_particle[abs(type_parton[dipole_candidate[i_a].no_R_emitter_2()])];
      if ((dipole_candidate[i_a].no_R_emitter_1() > 2 && type_parton[dipole_candidate[i_a].no_R_emitter_1()] > 0) || (dipole_candidate[i_a].no_R_emitter_1() < 3 && type_parton[dipole_candidate[i_a].no_R_emitter_1()] < 0)){temp_charge_factor = -temp_charge_factor;}
      if ((dipole_candidate[i_a].no_R_emitter_2() > 2 && type_parton[dipole_candidate[i_a].no_R_emitter_2()] > 0) || (dipole_candidate[i_a].no_R_emitter_2() < 3 && type_parton[dipole_candidate[i_a].no_R_emitter_2()] < 0)){temp_charge_factor = -temp_charge_factor;}
      cout << "temp_charge_factor = " << temp_charge_factor << endl;
    }
    else {
      temp_charge_factor = 0.;
    }
    */
    temp_symmetry_factor = temp_symmetry_factor / dipole[0].symmetry_factor();
    ///    double temp_symmetry_factor = dipole_data[i_a][1] / dipole[0].symmetry_factor();
    int temp_sum_channel = dipole[dipole.size() - 1].sum_channel() + psi_MC_n_channel_phasespace[x_a];
    int temp_massive;
    if (psi_M[abs(type_parton[dipole_candidate[i_a].no_R_emitter_1()])] == 0. && psi_M[abs(type_parton[dipole_candidate[i_a].no_R_emitter_2()])] == 0. && psi_M[abs(type_parton[dipole_candidate[i_a].no_R_spectator()])] == 0.){temp_massive = 0;}
    else {temp_massive = 1;}
    dipole.push_back(dipole_set(dipole_candidate[i_a], temp_no_map, temp_o_map, temp_no_prc, temp_o_prc, psi_MC_n_channel_phasespace[x_a], temp_sum_channel, temp_charge_factor, temp_symmetry_factor, temp_massive));

    int no_this_dipole = dipole.size() - 1;
    dipole[no_this_dipole].contribution_order_alpha_s = basic_order_alpha_s;
    dipole[no_this_dipole].contribution_order_alpha_e = basic_order_alpha_e;
    dipole[no_this_dipole].contribution_order_interference = basic_order_interference;

    vector<int> temp_singularity(2);
    temp_singularity[0] = dipole_candidate[i_a].no_R_emitter_1();
    temp_singularity[1] = dipole_candidate[i_a].no_R_emitter_2();
    sort(temp_singularity.begin(), temp_singularity.end());
    int new_region = 1;
    for (int i_x = 0; i_x < list_singular_regions.size(); i_x++){if (temp_singularity == list_singular_regions[i_x]){new_region = 0; break;}}
    if (new_region == 1){
      list_singular_regions.push_back(temp_singularity);
      singular_region_name[temp_singularity[0]][temp_singularity[1]] = "p" + pa_name[temp_singularity[0]] + "." + "p" + pa_name[temp_singularity[1]];
    }
  }
  logger << LOG_DEBUG << "QEW dipoles determined: " << dipole.size() - 1 << " dipoles contribute. " << endl;
  for (int i_a = 0; i_a < dipole.size(); i_a++){
    stringstream temp;
    temp << "dipole " << i_a << ":   " << setw(15) << left << dipole[i_a].name() << ":   ";
    temp << "   type_correction = " << dipole[i_a].type_correction();
    temp << "   massive = " << dipole[i_a].massive();
    temp << "   type_dipole = " << dipole[i_a].type_dipole();
    temp << "   type_splitting = " << dipole[i_a].type_splitting();
    logger << LOG_DEBUG << "     " << temp.str() << endl;
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
