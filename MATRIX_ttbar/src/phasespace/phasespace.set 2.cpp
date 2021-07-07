#include "../include/classes.cxx"

//#include "../include/definitions.observable.set.cxx"
////////////////////
//  constructors  //
////////////////////
phasespace_set::phasespace_set(){}
phasespace_set::phasespace_set(inputparameter_set & isi, contribution_set & _csi){
  Logger logger("phasespace_set::phasespace_set (isi)");
  logger << LOG_DEBUG << "called" << endl;

  csi = &_csi;
 
  coll_choice = isi.coll_choice;

  for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){
    logger << LOG_DEBUG << "type_parton[" << i_p << "] = " << csi->type_parton[0][i_p] << endl;
  }

  // !!! must be changed: type_parton is different in oss and pss !!! should be solved now...
  hcf = pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * csi->n_particle - 4));

  xb_max = intpow(2, csi->n_particle + 2);
  xb_max_dipoles = intpow(2, csi->n_particle + 2 - 1);
  xb_out_dipoles = xb_max_dipoles - 4;
  
  no_random = 3 * csi->n_particle - 4;
  r.resize(no_random + 1);

  // 'random_MC' should replace 'r' at some point !!!

  MC_opt_end = 0; // to avoid 'uninitialized' warning !!!
  end_optimization = 0; // to avoid 'uninitialized' warning !!!
  
  random_MC.resize(no_random + 1);
  
  random_MC_tau.resize(2);
  random_tau.resize(2);
  random_x12.resize(2);
  random_z.resize(3, vector<double> (2));

  RA_x_a = 0; // to avoid warning !!!

  initialization_random_number_generator(isi);
  
  logger << LOG_DEBUG << "isi.MCweight_in_directory = " << isi.MCweight_in_directory << endl;
  logger << LOG_DEBUG << "csi->subprocess = " << csi->subprocess << endl;
  string dir_MCweights_in_contribution = "../../../../" + isi.MCweight_in_directory + "/weights";
  random_manager = randommanager(sran, dir_MCweights_in_contribution, csi->subprocess);

  //  for (int i = 0; i < isi.zwahl; i++){x = ran(sran);}


  o_map.resize(1);
  o_map[0].resize(3 + csi->n_particle, 0);
  //  o_map.resize(3 + csi->n_particle, 0);

  o_prc.resize(1);
  o_prc[0].resize(3 + csi->n_particle, 0);

  switch_MC = isi.switch_MC;
  switch_MC_tau = isi.switch_MC_tau;
  switch_MC_x_dipole = isi.switch_MC_x_dipole;

  switch_IS_MC = isi.switch_IS_MC;
  switch_IS_tau = isi.switch_IS_tau;
  switch_IS_x1x2 = isi.switch_IS_x1x2;
  switch_IS_z1z2 = isi.switch_IS_z1z2;
  switch_n_events_opt = isi.switch_n_events_opt;

  switch_use_alpha_after_IS = isi.switch_use_alpha_after_IS;
  switch_step_mode_grid = isi.switch_step_mode_grid;

  switch_RS_mapping = isi.switch_RS_mapping;

  tau_0_num.resize(1, 1.);

  i_gen = 0;
  i_rej = 0;
  i_acc = 0;
  i_nan = 0;
  i_tec = 0;

  last_step_mode = 0;

  user = isi.user;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::initialization_random_number_generator(inputparameter_set & isi){
  Logger logger("phasespace_set::initialization_random_number_generator");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  //////////////////////////////////////////////////
  //  random-number initialization determination  //
  //////////////////////////////////////////////////

  double x = 0.;
  sran.resize(3);
  if (csi->subprocess != ""){
    // new implementation based on trivial prime number generator:
    int random_per_psp = 3 * csi->n_particle - 4 + 2 + 2;
    int random_sheet = isi.zwahl / random_per_psp;
    int max_n_prime_number = (random_sheet + 1) * 3;
    logger << LOG_DEBUG << "max_n_prime_number = " << max_n_prime_number << endl;
    vector<int> prime_number(1, 2);
    logger << LOG_DEBUG << "random_sheet no. 0" << endl;
    logger << LOG_DEBUG << "prime_number[" << setw(4) << prime_number.size() - 1 << "] = " << setw(10) << prime_number[prime_number.size() - 1] << "   sqrt = " << sqrt(prime_number[prime_number.size() - 1]) << endl;
    int number = 1;
    while (prime_number.size() < max_n_prime_number){
      number += 2;
      double max_number = sqrt(number);
      int flag = 0;
      for (int i_n = 0; i_n < prime_number.size(); i_n++){
	if (prime_number[i_n] > max_number){break;}
	if (number % prime_number[i_n] == 0){flag = 1; break;} 
      }
      if (!flag){
	if (prime_number.size() % 3 == 0){logger << LOG_DEBUG << "random_sheet no. " << prime_number.size() / 3 << endl;}
	prime_number.push_back(number);
	logger << LOG_DEBUG << "prime_number[" << setw(4) << prime_number.size() - 1 << "] = " << setw(10) << prime_number[prime_number.size() - 1] << "   sqrt = " << sqrt(prime_number[prime_number.size() - 1]) << endl;
      }
    }

    logger << LOG_DEBUG << "selected random_sheet: " << random_sheet << endl;
    sran[0] = sqrt(prime_number[prime_number.size() - 3]);
    sran[1] = sqrt(prime_number[prime_number.size() - 2]);
    sran[2] = sqrt(prime_number[prime_number.size() - 1]);

    for (int i_z = 0; i_z < 3; i_z++){sran[i_z] = sran[i_z] - int(sran[i_z]);}

    int zwahl_shift = isi.zwahl % random_per_psp;

    logger << LOG_DEBUG << "random_per_psp = " << random_per_psp << endl;
    logger << LOG_DEBUG << "isi.zwahl      = " << isi.zwahl << endl;
    logger << LOG_DEBUG << "random_sheet   = " << random_sheet << endl;
    
    for (int i = 0; i < zwahl_shift; i++){x = ran(sran);}
    logger << LOG_DEBUG << "x              = " << x << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::calculate_g_tot(){
  Logger logger("phasespace_set::calculate_g_tot");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  ///  logger << LOG_INFO << "g_MC = " << g_MC << endl;
  for (int j = 0; j < MC_phasespace.alpha.size(); j++){
    if (munich_isnan(MC_phasespace.alpha[j])){
      ///	logger << LOG_INFO << "MC_phasespace.g_channel[" << j << "] = " << MC_phasespace.g_channel[j] << "   MC_phasespace.alpha[" << j << "] = " << MC_phasespace.alpha[j] << endl;
      }
  }
  ///  logger << LOG_INFO << "g_pdf = " << g_pdf << endl;
  ///  logger << LOG_INFO << "g_tot = " << g_tot << endl;
  //  logger << LOG_INFO << "g_MC = " << g_MC << endl;
  
  g_MC = 0.;
  for (int j = 0; j < MC_phasespace.alpha.size(); j++){g_MC += MC_phasespace.g_channel[j] * MC_phasespace.alpha[j];}
  for (int j = 0; j < MC_phasespace.alpha.size(); j++){logger << LOG_DEBUG_VERBOSE << "MC_phasespace.alpha[" << j << "] = " << MC_phasespace.alpha[j] << endl;}
  for (int j = 0; j < MC_phasespace.g_channel.size(); j++){logger << LOG_DEBUG_VERBOSE << "MC_phasespace.g_channel[" << j << "] = " << MC_phasespace.g_channel[j] << endl;}
  
  if (g_global_NWA != 1.){g_MC *= g_global_NWA;} 
  if (coll_choice > 0){g_tot = g_MC * g_pdf;}
  else {g_tot = g_MC;}
  if (switch_IS_mode_phasespace == 1 || switch_IS_mode_phasespace == 3){g_tot = g_tot * MC_g_IS_global;}

  logger << LOG_DEBUG_VERBOSE << setw(20) << "MC_tau.channel" << " = " << MC_tau.channel << endl;
   logger << LOG_DEBUG_VERBOSE << setw(20) << "g_MC" << " = " << setprecision(20) << setw(28) << g_MC << "   " << double2hexastr(g_MC) << endl;
  logger << LOG_DEBUG_VERBOSE << setw(20) << "g_pdf" << " = " << setprecision(20) << setw(28) << g_pdf << "   " << double2hexastr(g_pdf) << endl;
  logger << LOG_DEBUG_VERBOSE << setw(20) << "g_tot" << " = " << setprecision(20) << setw(28) << g_tot << "   " << double2hexastr(g_tot) << endl;
  logger << LOG_DEBUG_VERBOSE << setw(20) << "xbs_all[0][0]" << " = " << setprecision(20) << setw(28) << xbs_all[0][0] << "   " << double2hexastr(xbs_all[0][0]) << endl;
  
  logger << LOG_DEBUG_VERBOSE << "g_MC = " << g_MC << endl;
  logger << LOG_DEBUG_VERBOSE << "g_pdf = " << g_pdf << endl;

  ps_factor = hcf / (g_tot * xbs_all[0][0]);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::calculate_IS(){
  Logger logger("phasespace_set::calculate_IS");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
  logger << LOG_DEBUG_VERBOSE << "i_gen = " << setw(12) << i_gen << "   i_acc = " << setw(12) << i_acc << "   i_rej = " << setw(12) << i_rej << "   i_nan = " << setw(12) << i_nan << "   i_tec = " << setw(12) << i_tec << endl;
  //  logger << LOG_DEBUG_VERBOSE << "n(gen. psp's -> i) = " << setw(12) << i << "   n(rej. psp's -> jets[0]) = " << setw(12) << jets[0] << "   n(acc. psp's -> jets[1]) = " << setw(12) << jets[1] << "   n(nan. psp's -> nancut) = " << setw(12) << nancut << endl;

  xbp_all = start_xbp_all;
  xbs_all = start_xbs_all;
  xbsqrts_all = start_xbsqrts_all;
  i_gen++;
  //  i++;

  if (coll_choice != 0){
    logger << LOG_DEBUG_VERBOSE << "start_xbs_all[0][xb_max - 4 = " << xb_max - 4 << "] = " << start_xbs_all[0][xb_max - 4] << endl;

    if (start_xbs_all[0][xb_max - 4] != 0.){
      if (switch_IS_x1x2){
	calculate_initial_tau_fixed_x1x2_IS();
      }
      else {
	calculate_initial_tau_fixed_x1x2();
      }
    }
    else {
      if (switch_IS_tau != 0 && switch_IS_x1x2 != 0){
	calculate_initial_tau_IS_x1x2_IS();
      }
      else if (switch_IS_tau != 0 && switch_IS_x1x2 == 0){
	calculate_initial_tau_IS_x1x2();
      }
      else if (switch_IS_tau == 0 && switch_IS_x1x2 != 0){
	calculate_initial_tau_x1x2_IS();
      }
      else {
	calculate_initial_tau_x1x2();
      }
    }

    boost = (x_pdf[1] - x_pdf[2]) / (x_pdf[1] + x_pdf[2]);
    g_pdf = g_tau * g_x1x2;
    logger << LOG_DEBUG_VERBOSE << "g_pdf = " << g_pdf << endl;
  }
  randomvector(sran, no_random, r);
  logger << LOG_DEBUG_VERBOSE << "r.size() = " << r.size() << endl;
  // remove MC_channel !!!
  if (r.size() > 0){for (int j = 0; j < MC_phasespace.beta.size(); j++){if (r[0] <= MC_phasespace.beta[j]){MC_phasespace.channel = j; break;}}}
  //  if (r.size() > 0){for (int j = 0; j < MC_phasespace.beta.size(); j++){if (r[0] <= MC_phasespace.beta[j]){MC_channel = j; MC_phasespace.channel = j; break;}}}

  //  MC_channel_phasespace = MC_channel;
  MC_g_IS_global = 1.;

  // more informative name for 'weight_IS', e.g. switch_IS_mode_phasespace ???
  logger << LOG_DEBUG_VERBOSE << "MC_channel_phasespace = " << MC_channel_phasespace << endl;
  if (switch_IS_mode_phasespace == 1 || switch_IS_mode_phasespace == 2){
    double g_IS_temp;
    for (int k=0; k < no_random; k++) {
      phasespace_randoms[MC_phasespace.channel * no_random + k]->get_random(r[k + 1], g_IS_temp);
      MC_g_IS_global *= g_IS_temp;
    }
  }

  if (o_map[0][1] == 2 && o_map[0][2] == 1){swap(xbp_all[0][1], xbp_all[0][2]);}

  /*
  for (int i_x = 0; i_x < xbp_all[0].size(); i_x++){
    //    logger << LOG_DEBUG << "start_xbp[" << setw(3) << i_x << "] = " << start_xbp_all[0][i_x] << " " << setw(23) << setprecision(15) << start_xbsqrts_all[0][i_x] << " " << setw(23) << setprecision(15) << start_xbs_all[0][i_x] << endl;
    //    logger << LOG_DEBUG << "xbp[" << setw(3) << i_x << "] = " << xbp_all[0][i_x] << " " << setw(23) << setprecision(15) << xbsqrts_all[0][i_x] << " " << setw(23) << setprecision(15) << xbs_all[0][i_x] << endl;
    if (xbs_all[0][i_x] != 0. || xbp_all[0][i_x] != nullvector){
      logger << LOG_DEBUG_VERBOSE << "xbp[" << setw(3) << i_x << "] = " << xbp_all[0][i_x] << " " << setw(23) << setprecision(15) << xbsqrts_all[0][i_x] << " " << setw(23) << setprecision(15) << xbs_all[0][i_x] << endl;
    }
  }
  */
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void phasespace_set::calculate_IS_RA(vector<dipole_set> & _dipole){
  Logger logger("phasespace_set::calculate_IS");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
  for (int i_a = 0; i_a < _dipole.size(); i_a++){
    if (MC_phasespace.channel < _dipole[i_a].sum_channel()){
      RA_x_a = i_a;
      if (i_a == 0){MC_channel_phasespace = MC_phasespace.channel;}
      else {MC_channel_phasespace = MC_phasespace.channel - MC_sum_channel_phasespace[i_a - 1];}
      break;
    }
  }
  int temp_zero = 0;
  if (RA_x_a > 0){temp_zero = _dipole[RA_x_a - 1].sum_channel();}
  logger << LOG_DEBUG_VERBOSE << "MC_phasespace.channel = " << MC_phasespace.channel << endl;
  logger << LOG_DEBUG_VERBOSE << "MC_channel_phasespace = " << MC_channel_phasespace << "   RA_x_a = " << RA_x_a << "   channel[" << RA_x_a << "] = " << MC_phasespace.channel - temp_zero << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


/*
void phasespace_set::calculate_IS_tau_x1x2_MC(){
  Logger logger("phasespace_set::calculate_IS_tau_x1x2_MC");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<double> z_tau_MC(1);
  randomvector(sran, 0, z_tau_MC);
  vector<double> z_tau(2);
  randomvector(sran, 1, z_tau);
  logger << LOG_DEBUG_VERBOSE << "begin tau_channel = " << tau_channel << endl;
  for (int j = 0; j < tau_MC_beta.size(); j++){if (z_tau_MC[0] <= tau_MC_beta[j]){tau_MC_channel = j; break;}}
  logger << LOG_DEBUG_VERBOSE << "tau_MC_channel = " << tau_MC_channel << endl;
  if (tau_MC_channel == 0){
    for (int j = 0; j < tau_beta.size(); j++){if (z_tau[0] <= tau_beta[j]){tau_channel = j; break;}}
    logger << LOG_DEBUG_VERBOSE << "after1 tau_channel = " << tau_channel << endl;
    double random_tau = (double(tau_channel) + z_tau[1]) / double(tau_beta.size());
    x_pdf[0] = h_propto_pot(random_tau, tau_0, exp_pdf);
  }
  else {
    if (tau_MC_map[tau_MC_channel] < 0){x_pdf[0] = c_propagator_Breit_Wigner(z_tau[0], -tau_MC_map[tau_MC_channel], tau_0_s_had, s_had) / s_had;}
    ///    else if (tau_MC_map[tau_MC_channel] > 0){x_pdf[0] = c_propagator_vanishing_width(z_tau[0], M2[tau_MC_map[tau_MC_channel]], tau_0_s_had, s_had, nuxs) / s_had;}
    for (int j = 0; j < tau_MC_tau_gamma[0].size(); j++){
      logger << LOG_DEBUG_VERBOSE << "tau_MC_tau_gamma[0][" << j << "] = " << tau_MC_tau_gamma[0][j] << endl;
      if (x_pdf[0] <= tau_MC_tau_gamma[0][j]){
	tau_channel = j; 
	logger << LOG_DEBUG_VERBOSE << "tau_channel = " << tau_channel << "   tau_min = " << setw(23) << setprecision(15) << h_propto_pot(tau_alpha[j], tau_0, exp_pdf) << " < " << setw(23) << setprecision(15) << x_pdf[0] << " < " << setw(23) << setprecision(15) << h_propto_pot(tau_alpha[j + 1], tau_0, exp_pdf) << " = tau_max" << endl;

	break;
      }
    }
    logger << LOG_DEBUG_VERBOSE << "after2 tau_channel = " << tau_channel << endl;
  }
  double tau_s_had = x_pdf[0] * s_had;
  logger << LOG_DEBUG_VERBOSE << "before g_tau_MC" << endl;
  vector<double> g_tau_MC(tau_MC_alpha.size());
  logger << LOG_DEBUG_VERBOSE << "g_tau_MC.size() = " << g_tau_MC.size() << endl;
  for (int j = 0; j < g_tau_MC.size(); j++){
    logger << LOG_DEBUG_VERBOSE << "j = " << j << endl;
    logger << LOG_DEBUG_VERBOSE << "tau_MC_map[" << j << "] = " << tau_MC_map[j] << endl;
    if (j == 0){
      logger << LOG_DEBUG_VERBOSE << "done tau_alpha.size() = " << tau_alpha.size() << endl;
      logger << LOG_DEBUG_VERBOSE << "done tau_channel = " << tau_channel << endl;

      g_tau_MC[j] = g_propto_pot(x_pdf[0], tau_0, exp_pdf) * (tau_alpha[tau_channel] * tau_alpha.size());
    }
    else {
      if (tau_MC_map[j] < 0){g_tau_MC[j] = s_had * g_propagator_Breit_Wigner(tau_s_had, -tau_MC_map[j], tau_0_s_had, s_had);}
      ///      else if (tau_MC_map[j] > 0){g_tau_MC[j] = s_had * g_propagator_vanishing_width(tau_s_had, M2[tau_MC_map[j]], tau_0_s_had, s_had, nuxs);}
    }
  }
  logger << LOG_DEBUG_VERBOSE << "done g_tau_MC.size() = " << g_tau_MC.size() << endl;

  g_tau = 0.;
  for (int i = 0; i < g_tau_MC.size(); i++){
    g_tau += tau_MC_alpha[i] * g_tau_MC[i];
  }
  logger << LOG_DEBUG_VERBOSE << "before z_x1x2" << endl;

  vector<double> z_x1x2(2);
  randomvector(sran, 1, z_x1x2);
  for (int j = 0; j < x1x2_beta.size(); j++){if (z_x1x2[0] <= x1x2_beta[j]){x1x2_channel = j; break;}}
  double random_x1x2 = (double(x1x2_channel) + z_x1x2[1]) / double(x1x2_beta.size());
  x_pdf[1] = exp(random_x1x2 * log(x_pdf[0]));
  x_pdf[2] = x_pdf[0] / x_pdf[1];
  g_x1x2 = 1. / (-log(x_pdf[0])) * (x1x2_alpha[x1x2_channel] * x1x2_alpha.size());
  double sqrttau = sqrt(x_pdf[0]);
  xbp_all[0][1] = fourvector(sqrttau * E, 0., 0., sqrttau * E);
  xbp_all[0][2] = fourvector(sqrttau * E, 0., 0., -sqrttau * E);
  xbp_all[0][0] = xbp_all[0][1] + xbp_all[0][2];
  xbs_all[0][0] = 4. * x_pdf[0] * pow(E, 2);
  xbsqrts_all[0][0] = sqrt(xbs_all[0][0]);
  //  int xb_max = xbp[0].size();
  xbp_all[0][xb_max - 4] = xbp_all[0][0];
  xbs_all[0][xb_max - 4] = xbs_all[0][0];
  xbsqrts_all[0][xb_max - 4] = xbsqrts_all[0][0];
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/

void phasespace_set::calculate_IS_CX(){
  Logger logger("phasespace_set::calculate_IS_CX");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  //  for (int i_x = 0; i_x < 3; i_x++){xz_pdf[i_x][0] = x_pdf[i_x];}
  for (int i_z = 1; i_z < 3; i_z++){
    QT_random_z[i_z]->get_random(QT_random_IS, QT_g_IS_);
    z_coll[i_z] = pow(x_pdf[i_z], QT_eps + (1. - QT_eps) * QT_random_IS);
    //    g_z_coll[i_z] = -1. / (log(x_pdf[i_z]) * (1. - QT_eps)) * QT_g_IS_;
    g_z_coll[i_z] = -1. / (log(x_pdf[i_z]) * (1. - QT_eps));
    g_pdf *= QT_g_IS_;
    //    cout << "z_coll[" << i_z << "] = " << z_coll[i_z] << endl;
    //    cout << "g_z_coll[" << i_z << "] = " << g_z_coll[i_z] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::calculate_IS_QT(){
  Logger logger("phasespace_set::calculate_IS_QT");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  QT_random_qt->get_random(QT_random_qt2,QT_g_IS_);
  //  cout << "qt: QT_random_qt2 = " << setw(23) << setprecision(15) << QT_random_qt2 << "   QT_g_IS_ = " << setw(23) << setprecision(15) << QT_g_IS_ << endl;
  g_pdf *= QT_g_IS_;

  if (type_contribution == "CT" ||
      type_contribution == "CT2"){
    // the (Bessel transformed) large logarithms vanish up to double precision accuracy for arguments larger than 40 -> use this as upper bound
    // 1.261 ~ b0*b0
    double maxqt2;
    double Qres_tmp=0.0, Qres_lower, Qres_upper;
    
    // Qres appears at two places here: as the reference value for the lower and for the upper qT cut
    // conceptually, these are partially unrelated
    // in fixed order computations, we have to use Qres as the lower cut off scale, because we precompute the qT integrals
    // in resummed computations, we explicitly integrate over qT, so this is unnecessary; in fact
    // it is more efficient to use mF as the lower reference scale, as qT/mF determine the universal limit
    // the upper cut off scale should alway be given by Qres, as otherwise qT/Qres can become huge, resulting in numerical
    // problems in the counterterm weight computation
    
    // Qres != mF does not work at the moment; has to be implemented in process specific cut files
    logger << LOG_DEBUG_VERBOSE << "do_resummation = " << do_resummation << endl;
    logger << LOG_DEBUG_VERBOSE << "Qres = " << Qres << endl;

    assert(Qres==0 || do_resummation==0);
    

    if (dynamical_Qres) {
      Qres_tmp = Qres_prefactor*xbsqrts_all[0][0];
    } else {
      Qres_tmp = Qres;
    }
    if (Qres_tmp == 0) {
      // default choice
      Qres_tmp = xbsqrts_all[0][0];
    }
    
    Qres_upper = Qres_tmp;
    
    if (do_resummation){
      Qres_lower = xbsqrts_all[0][0];
    } else {
      Qres_lower = Qres_tmp;
    }
    
    if (contribution_order_alpha_s[0] == 1){maxqt2 = 1.261 * 30. * 20. * Qres_upper * Qres_upper;}
    else {maxqt2 = 1.261 * 40. * 50. * Qres_upper * Qres_upper;}
    
    if (switch_qTcut == 1){
      double r0 = min_qTcut / 100.;
      //  version with cut on pT/Qres_cut
      QT_qt2 = r0 * r0 * pow(Qres_lower,2) * exp(log(maxqt2 / r0 / r0 / Qres_lower / Qres_lower) * QT_random_qt2);
//       cout << "qT=" << QT_qt2 << ", Qres=" << Qres_upper << ", " << Qres_lower << ", ratio=" << QT_qt2/Qres_upper << endl;
      QT_jacqt2 = log(maxqt2 / r0 / r0 / Qres_lower / Qres_lower) * QT_qt2;
      //  version with cut on pT/m_inv
      //    QT_qt2 = r0 * r0 * xbs_all[0][0] * exp(log(maxqt2 / r0 / r0 / xbs_all[0][0]) * QT_random_qt2);
      //    QT_jacqt2 = log(maxqt2 / r0 / r0 / xbs_all[0][0]) * QT_qt2;
    }
    else if (switch_qTcut == 2){
      // old version with cut on pT only (not relative to m_inv)
      //   qt2=min_qTcut*min_qTcut*exp(log(maxqt2/min_qTcut/min_qTcut)*random_qt2);
      //   jacqt2=log(maxqt2/min_qTcut/min_qTcut)*qt2;
      logger << LOG_FATAL << "old qT variation is deprecated" << endl;
      exit(0);
    }
  }

  else if (type_contribution == "NLL_LO" || type_contribution == "NLL_NLO" || type_contribution == "NNLL_LO" || type_contribution == "NNLL_NLO" || type_contribution == "NNLL_NNLO"){
    if (do_resummation) {
      min_qTcut = 0.1; // ??? why is that set manually (and overwrites the set value of min_qTcut) ???
      
      double maxqt2 = pow(2 * E, 2);
      maxqt2 = pow(E, 2); // ???

      if (lambda_qt2 == 0 || lambda_qt2 == 1) { // ??? meaning of lambda_qt2 ??? does it only choose the mapping of QT_qt2 ???
	QT_qt2 = min_qTcut * min_qTcut * exp(log(maxqt2 / min_qTcut / min_qTcut) * QT_random_qt2);
	QT_jacqt2 = log(maxqt2 / min_qTcut / min_qTcut) * QT_qt2;
      }

      //  random_qt2=0;

      //  double lambda=0.0355064;
      //  double C=-1.0/lambda*(exp(-lambda*sqrt(maxqt2))-exp(-lambda*min_qTcut));
      //  double y0=-1.0/lambda/C*exp(-lambda*min_qTcut);
      ////  cout << exp(lambda*sqrt(maxqt2)) << ", " << exp(lambda*min_qTcut) << endl;
      ////  cout << C << ", " << y0 << endl;
      //  double qt=-1.0/lambda*log(-lambda*C*(random_qt2+y0));
      //  double jacqt=C*exp(lambda*qt);
      //  qt2=qt*qt;
      //  jacqt2=jacqt*2*qt;

      else {
	double y0 = 1. / (pow(maxqt2 / min_qTcut / min_qTcut, 1 - lambda_qt2) - 1);
	double C = pow(min_qTcut, 2 * (1 - lambda_qt2)) / (1 - lambda_qt2) / y0;

	QT_qt2 = pow((1 - lambda_qt2) * C * (QT_random_qt2 + y0), 1. / (1 - lambda_qt2));
	QT_jacqt2 = C * pow(QT_qt2, lambda_qt2);

	//    cout << C << ", " << y0 << ", " << (1-lambda_qt2)*C*(random_qt2+y0) << endl;
	//    cout << qt2 << ", " << maxqt2 << endl;
      }

      //  maxqt2 = 1e6;
      //  qt2 = min_qTcut * min_qTcut + random_qt2 * (maxqt2 - min_qTcut * min_qTcut);
      //  jacqt2 = (maxqt2 - min_qTcut * min_qTcut);
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
/*
// original version: check if everything works without this one...
void phasespace_set::calculate_IS_QT(){
  Logger logger("phasespace_set::calculate_IS_QT");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  QT_random_qt->get_random(QT_random_qt2,QT_g_IS_);
  //  cout << "qt: QT_random_qt2 = " << setw(23) << setprecision(15) << QT_random_qt2 << "   QT_g_IS_ = " << setw(23) << setprecision(15) << QT_g_IS_ << endl;
  g_pdf *= QT_g_IS_;


  logger << LOG_DEBUG_VERBOSE << "type_correction = " << type_correction << "   type_contribution = " << type_contribution << endl;

  if (type_correction == "CT" || type_correction == "CT2"){
    // the (Bessel transformed) large logarithms vanish up to double precision accuracy for arguments larger than 40 -> use this as upper bound
    // 1.261 ~ b0*b0
    double maxqt2;
    double Qres;
    if (do_resummation){// resummation; there is (usually) no point in using Q_res as a reference scale for the cut. in fact, it only slows down convergence
      Qres = xbsqrts_all[0][0]; // should be the same as Qres = sqrt(xbs_all[0][0]); !!!
    }
    else { // f.o.; need to use correct lower limit for qT integral
      if (Qres_cut == 0 && !dynamical_Qres){Qres = xbsqrts_all[0][0];} // not specified -> assume standard f.o. computation // should be the same as Qres = sqrt(xbs_all[0][0]); !!!
      else if (dynamical_Qres){Qres = xbsqrts_all[0][0] * Qres_prefactor;} // should be the same as Qres = sqrt(xbs_all[0][0]) * Qres_prefactor; !!!
      else {Qres = Qres_cut;} // use specified fixed Qres
    }
    
    // ??? why is "maxqt2" multiplied by "Qres^2" here, even though maxqt2 only appears via (maxqt2 / Qres^2) ??? (or with Qres^2 -> xbs_all[0][0])
    if (contribution_order_alpha_s[0] == 1){maxqt2 = 1.261 * 30. * 20. * Qres * Qres;}
    else if (contribution_order_alpha_s[0] == 2){maxqt2 = 1.261 * 40. * 50. * Qres * Qres;}
    //  if (contribution_order_alpha_s[0] == 1){maxqt2 = 1.261 * 30. * 20. * xbs_all[0][0];}
    //  else if (contribution_order_alpha_s[0] == 2){maxqt2 = 1.261 * 40. * 50. * xbs_all[0][0];} // TODO: is this really necessary?
    
    if (switch_qTcut == 1){
      double r0 = min_qTcut / 100.;
      //  version with cut on pT/Qres_cut
      QT_qt2 = r0 * r0 * xbs_all[0][0] * exp(log(maxqt2 / r0 / r0 / Qres / Qres) * QT_random_qt2);
      QT_jacqt2 = log(maxqt2 / r0 / r0 / Qres / Qres) * QT_qt2;
      //  version with cut on pT/m_inv
      //    QT_qt2 = r0 * r0 * xbs_all[0][0] * exp(log(maxqt2 / r0 / r0 / xbs_all[0][0]) * QT_random_qt2);
      //    QT_jacqt2 = log(maxqt2 / r0 / r0 / xbs_all[0][0]) * QT_qt2;
    }
    else if (switch_qTcut == 2){
      // old version with cut on pT only (not relative to m_inv)
      //   qt2=min_qTcut*min_qTcut*exp(log(maxqt2/min_qTcut/min_qTcut)*random_qt2);
      //   jacqt2=log(maxqt2/min_qTcut/min_qTcut)*qt2;
      logger << LOG_FATAL << "old qT variation is deprecated" << endl;
      exit(0);
    }
  }

  else if (type_correction == "VT" || type_correction == "VT2"){
    if (do_resummation) {
      min_qTcut = 0.1; // ??? why is that set manually (and overwrites the set value of min_qTcut) ???
      
      double maxqt2 = pow(2 * E, 2);
      maxqt2 = pow(E, 2); // ???

      if (lambda_qt2 == 0 || lambda_qt2 == 1) { // ??? meaning of lambda_qt2 ??? does it only choose the mapping of QT_qt2 ???
	QT_qt2 = min_qTcut * min_qTcut * exp(log(maxqt2 / min_qTcut / min_qTcut) * QT_random_qt2);
	QT_jacqt2 = log(maxqt2 / min_qTcut / min_qTcut) * QT_qt2;
      }

      //  random_qt2=0;

      //  double lambda=0.0355064;
      //  double C=-1.0/lambda*(exp(-lambda*sqrt(maxqt2))-exp(-lambda*min_qTcut));
      //  double y0=-1.0/lambda/C*exp(-lambda*min_qTcut);
      ////  cout << exp(lambda*sqrt(maxqt2)) << ", " << exp(lambda*min_qTcut) << endl;
      ////  cout << C << ", " << y0 << endl;
      //  double qt=-1.0/lambda*log(-lambda*C*(random_qt2+y0));
      //  double jacqt=C*exp(lambda*qt);
      //  qt2=qt*qt;
      //  jacqt2=jacqt*2*qt;

      else {
	double y0 = 1. / (pow(maxqt2 / min_qTcut / min_qTcut, 1 - lambda_qt2) - 1);
	double C = pow(min_qTcut, 2 * (1 - lambda_qt2)) / (1 - lambda_qt2) / y0;

	QT_qt2 = pow((1 - lambda_qt2) * C * (QT_random_qt2 + y0), 1. / (1 - lambda_qt2));
	QT_jacqt2 = C * pow(QT_qt2, lambda_qt2);

	//    cout << C << ", " << y0 << ", " << (1-lambda_qt2)*C*(random_qt2+y0) << endl;
	//    cout << qt2 << ", " << maxqt2 << endl;
      }

      //  maxqt2 = 1e6;
      //  qt2 = min_qTcut * min_qTcut + random_qt2 * (maxqt2 - min_qTcut * min_qTcut);
      //  jacqt2 = (maxqt2 - min_qTcut * min_qTcut);
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/


void phasespace_set::determine_dipole_phasespace_RA(vector<dipole_set> & _dipole){
  Logger logger("phasespace_set::determine_dipole_phasespace_RA");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
  for (int xbi = 1; xbi < xbp_all[0].size(); xbi = xbi * 2){
    logger << LOG_DEBUG_VERBOSE << "xbp_all[" << 0 << "][" << xbi << "] = " << xbp_all[0][xbi] << endl;
  }

  for (int i_a = 1; i_a < _dipole.size(); i_a++){
    logger << LOG_DEBUG_VERBOSE << "construct phasespace of dipole no. " << i_a << ": " << _dipole[i_a].name() << endl; //"   " << _dipole[i_a].xy << endl;
    
    if (_dipole[i_a].type_dipole() == 1){
      if (_dipole[i_a].massive() == 0){phasespace_ij_k(i_a, _dipole);}
      else {phasespace_ij_k_massive(i_a, _dipole);}
    }
    if (_dipole[i_a].type_dipole() == 2){
      if (_dipole[i_a].massive() == 0){phasespace_ij_a(i_a, _dipole);}
      else {phasespace_ij_a_massive(i_a, _dipole);}
    }
    if (_dipole[i_a].type_dipole() == 3){phasespace_ai_k(i_a, _dipole);}
    if (_dipole[i_a].type_dipole() == 5){phasespace_ai_b(i_a, _dipole);}
    logger << LOG_DEBUG_VERBOSE << "construct phasespace of dipole no. " << i_a << ": " << _dipole[i_a].name() << "   " << _dipole[i_a].xy << endl;
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

#define p_i xbp_all[0][bp_i]
#define p_j xbp_all[0][bp_j]
#define p_k xbp_all[0][bp_k]
#define pt_ij xbp_all[i_a][bpt_ij]
#define pt_k xbp_all[i_a][bpt_k]
#define y_ij_k dipole[i_a].xy

void phasespace_set::phasespace_ij_k(int i_a, vector<dipole_set> & dipole){
  static Logger logger("phasespace_ij_k");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  int bp_i = dipole[i_a].binary_R_emitter_1();
  int bp_j = dipole[i_a].binary_R_emitter_2();
  int bp_k = dipole[i_a].binary_R_spectator();
  int bpt_ij = dipole[i_a].binary_A_emitter();
  int bpt_k = dipole[i_a].binary_A_spectator();
  logger << LOG_DEBUG_VERBOSE << "bp_i   = " << bp_i << endl;
  logger << LOG_DEBUG_VERBOSE << "bp_j   = " << bp_j << endl;
  logger << LOG_DEBUG_VERBOSE << "bp_k   = " << bp_k << endl;
  logger << LOG_DEBUG_VERBOSE << "bpt_ij = " << bpt_ij << endl;
  logger << LOG_DEBUG_VERBOSE << "bpt_k  = " << bpt_k << endl;
  logger << LOG_DEBUG_VERBOSE << "y_ij_k = " << y_ij_k << endl;
  y_ij_k = (p_i * p_j) / (p_i * p_j + p_j * p_k + p_k * p_i);
  pt_ij = p_i + p_j - y_ij_k / (1. - y_ij_k) * p_k;
  pt_k = 1. / (1. - y_ij_k) * p_k;

  for (int xbi = 1; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else if (xbi == bp_j){}
    else if (xbi == bp_k){}
    else if (xbi < bp_j){xbp_all[i_a][xbi] = xbp_all[0][xbi];}
    else if (xbi > bp_j){xbp_all[i_a][xbi / 2] = xbp_all[0][xbi];}
    else {cout << "Should not happen!" << endl;}
  }
  xbp_all[i_a][0] = xbp_all[i_a][1] + xbp_all[i_a][2];
  xbs_all[i_a][0] = xbs_all[0][0];
  xbsqrts_all[i_a][0] = xbsqrts_all[0][0];
  int xbn_all = xbp_all[i_a].size() - 4;
  xbp_all[i_a][xbn_all] = xbp_all[i_a][0];
  xbs_all[i_a][xbn_all] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xbn_all] = xbsqrts_all[i_a][0];
  /*
  for (int i = 0; i < osi_p_parton[i_a].size(); i++){
    int xbi = intpow(2, i - 1);
    osi_p_parton[i_a][i] = xbp_all[i_a][xbi];
  }
  */
  logger << LOG_DEBUG_VERBOSE << "phasespace_ij_k finished" << endl;
}
#undef p_i 
#undef p_j 
#undef p_k 
#undef pt_ij 
#undef pt_k 

#define p_i xbp_all[0][bp_i]
#define p_j xbp_all[0][bp_j]
#define p_a xbp_all[0][bp_a]
#define pt_ij xbp_all[i_a][bpt_ij]
#define pt_a xbp_all[i_a][bpt_a]
#define x_ij_a dipole[i_a].xy
void phasespace_set::phasespace_ij_a(int i_a, vector<dipole_set> & dipole){
  static Logger logger("phasespace_ij_a");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  int bp_i = dipole[i_a].binary_R_emitter_1();
  int bp_j = dipole[i_a].binary_R_emitter_2();
  int bp_a = dipole[i_a].binary_R_spectator(); // bp_b not needed !!!
  int bpt_ij = dipole[i_a].binary_A_emitter();
  int bpt_a = dipole[i_a].binary_A_spectator(); // always bp_a == bpt_a

  x_ij_a = 1. - (p_i * p_j) / ((p_i + p_j) * p_a);
  pt_ij = p_i + p_j - (1 - x_ij_a) * p_a;
  pt_a = x_ij_a * p_a;

  for (int xbi = 1; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else if (xbi == bp_j){}
    else if (xbi == bp_a){}
    else if (xbi < bp_j){xbp_all[i_a][xbi] = xbp_all[0][xbi];}
    else if (xbi > bp_j){xbp_all[i_a][xbi / 2] = xbp_all[0][xbi];}
    else {cout << "Should not happen!" << endl;}
  }

  xbp_all[i_a][0] = xbp_all[i_a][1] + xbp_all[i_a][2];
  xbs_all[i_a][0] = x_ij_a * xbs_all[0][0];
  xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);
  int xbn_all = xbp_all[i_a].size() - 4;
  xbp_all[i_a][xbn_all] = xbp_all[i_a][0];
  xbs_all[i_a][xbn_all] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xbn_all] = xbsqrts_all[i_a][0];
  /*
  for (int i = 0; i < osi_p_parton[i_a].size(); i++){
    int xbi = intpow(2, i - 1);
    osi_p_parton[i_a][i] = xbp_all[i_a][xbi];
  }
  */
  logger << LOG_DEBUG_VERBOSE << "phasespace_ij_k finished" << endl;
}
#undef p_i 
#undef p_j 
#undef p_a 
#undef pt_ij 
#undef pt_a 
#undef x_ij_a

#define p_a xbp_all[0][bp_a]
#define p_i xbp_all[0][bp_i]
#define p_k xbp_all[0][bp_k]
#define pt_ai xbp_all[i_a][bpt_ai]
#define pt_k xbp_all[i_a][bpt_k]
#define x_ik_a dipole[i_a].xy
/*
#define p_a xbp_all[0][dipole[i_a].binary_R_emitter_1()]
#define p_i xbp_all[0][dipole[i_a].binary_R_emitter_2()]
#define p_k xbp_all[0][dipole[i_a].binary_R_spectator()]
#define pt_ai xbp_all[i_a][dipole[i_a].binary_A_emitter()]
#define pt_k xbp_all[i_a][dipole[i_a].binary_R_spectator()]
*/
void phasespace_set::phasespace_ai_k(int i_a, vector<dipole_set> & dipole){
  static Logger logger("phasespace_ai_k");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  int bp_a = dipole[i_a].binary_R_emitter_1();
  int bp_i = dipole[i_a].binary_R_emitter_2();
  int bp_k = dipole[i_a].binary_R_spectator();
  int bpt_ai = dipole[i_a].binary_A_emitter();
  int bpt_k = dipole[i_a].binary_A_spectator();

  x_ik_a = 1. - (p_i * p_k) / ((p_i + p_k) * p_a);
  //double u_i = (p_i * p_a) / ((p_i + p_k) * p_a);
  pt_k = p_i + p_k - (1 - x_ik_a) * p_a;
  pt_ai = x_ik_a * p_a;

  for (int xbi = 1; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else if (xbi == bp_k){}
    else if (xbi == bp_a){}
    else if (xbi < bp_i){xbp_all[i_a][xbi] = xbp_all[0][xbi];}
    else if (xbi > bp_i){xbp_all[i_a][xbi / 2] = xbp_all[0][xbi];}
    else {cout << "Should not happen!" << endl;}
  }

  xbp_all[i_a][0] = xbp_all[i_a][1] + xbp_all[i_a][2];
  xbs_all[i_a][0] = x_ik_a * xbs_all[0][0];
  xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);
  int xbn_all = xbp_all[i_a].size() - 4;
  xbp_all[i_a][xbn_all] = xbp_all[i_a][0];
  xbs_all[i_a][xbn_all] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xbn_all] = xbsqrts_all[i_a][0];
  /*
  for (int i = 0; i < osi_p_parton[i_a].size(); i++){
    int xbi = intpow(2, i - 1);
    osi_p_parton[i_a][i] = xbp_all[i_a][xbi];
  }
  */
  logger << LOG_DEBUG_VERBOSE << "phasespace_ai_k finished" << endl;
}
#undef p_i 
#undef p_k 
#undef p_a 
#undef pt_ai 
#undef pt_k 
#undef x_ik_a 
/*
#define p_a xbp_all[0][dipole[i_a].binary_R_emitter_1()]
#define p_i xbp_all[0][dipole[i_a].binary_R_emitter_2()]
#define p_b xbp_all[0][dipole[i_a].binary_R_spectator()]
#define bpt_ai bp_a
#define bpt_b bp_b
#define pt_ai xbp_all[i_a][dipole[i_a].binary_A_emitter()]
#define pt_b xbp_all[i_a][dipole[i_a].binary_R_spectator()]
*/
#define bpt_ai bp_a
#define bpt_b bp_b
#define p_a xbp_all[0][bp_a]
#define p_i xbp_all[0][bp_i]
#define p_b xbp_all[0][bp_b]
#define pt_ai xbp_all[i_a][bpt_ai]
#define pt_b xbp_all[i_a][bpt_b]
#define x_i_ab dipole[i_a].xy
void phasespace_set::phasespace_ai_b(int i_a, vector<dipole_set> & dipole){
  static Logger logger("phasespace_ai_b");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  int bp_a = dipole[i_a].binary_R_emitter_1();
  int bp_i = dipole[i_a].binary_R_emitter_2();
  int bp_b = dipole[i_a].binary_R_spectator();
  //  int bpt_ai = dipole[i_a].binary_A_emitter(); // always bpt_ai == bp_a
  //  int bpt_b = dipole[i_a].binary_A_spectator(); // always bpt_b == bp_b
  x_i_ab = 1. - (p_i * (p_a + p_b)) / (p_a * p_b);
  //  double v_i = (p_i * p_a) / (p_a * p_b);
  pt_ai = x_i_ab * p_a;
  fourvector K = p_a + p_b - p_i;
  fourvector Kt = pt_ai + p_b;
  fourvector KKt = K + Kt;
  double KKt2 = KKt.m2();
  double K2 = K.m2();
  //  cout << "phasespace_ai_b   xbp_all.size() = " << xbp_all.size() << endl;
  //  cout << "phasespace_ai_b   xbp_all[0].size() = " << xbp_all[0].size() << endl;
  //  cout << "phasespace_ai_b   xbp_all[" << i_a << "].size() = " << xbp_all[i_a].size() << endl;

  for (int xbi = 1; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_a){}
    else if (xbi == bp_b){xbp_all[i_a][xbi] = xbp_all[0][xbi];}
    else if (xbi == bp_i){}
    else if (xbi < bp_i){xbp_all[i_a][xbi] = xbp_all[0][xbi] - ((2. * (xbp_all[0][xbi] * KKt)) / KKt2) * KKt + (2. * (xbp_all[0][xbi] * K) / K2) * Kt;}
    else if (xbi > bp_i){xbp_all[i_a][xbi / 2] = xbp_all[0][xbi] - ((2. * (xbp_all[0][xbi] * KKt)) / KKt2) * KKt + (2. * (xbp_all[0][xbi] * K) / K2) * Kt;}
    else {cout << "Should not happen!" << endl;}
  }
  xbp_all[i_a][0] = xbp_all[i_a][1] + xbp_all[i_a][2];
  xbs_all[i_a][0] = x_i_ab * xbs_all[0][0];
  xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);
  int xbn_all = xbp_all[i_a].size() - 4;

  xbp_all[i_a][xbn_all] = xbp_all[i_a][0];
  xbs_all[i_a][xbn_all] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xbn_all] = xbsqrts_all[i_a][0];
  /*
  for (int i = 0; i < osi_p_parton[i_a].size(); i++){
    int xbi = intpow(2, i - 1);
    osi_p_parton[i_a][i] = xbp_all[i_a][xbi];
  }
  */
  logger << LOG_DEBUG_VERBOSE << "phasespace_ai_b finished" << endl;
}
#undef p_a 
#undef p_i 
#undef p_b 
#undef pt_ai 
#undef pt_b 
#undef bpt_ai
#undef bpt_b
#undef x_i_ab








#define p_i xbp_all[0][bp_i]
#define p_j xbp_all[0][bp_j]
#define p_k xbp_all[0][bp_k]
#define pt_ij xbp_all[i_a][bpt_ij]
#define pt_k xbp_all[i_a][bpt_k]
#define m2_i xbs_all[0][bp_i]
#define m2_j xbs_all[0][bp_j]
#define m2_k xbs_all[0][bp_k]
#define m2_ij xbs_all[i_a][bpt_ij]
#define y_ij_k dipole[i_a].xy
void phasespace_set::phasespace_ij_k_massive(int i_a, vector<dipole_set> & dipole){
  static Logger logger("phasespace_ij_k_massive");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  int bp_i = dipole[i_a].binary_R_emitter_1();
  int bp_j = dipole[i_a].binary_R_emitter_2();
  int bp_k = dipole[i_a].binary_R_spectator();
  int bpt_ij = dipole[i_a].binary_A_emitter();
  int bpt_k = dipole[i_a].binary_A_spectator();
  /*
  logger << LOG_DEBUG_VERBOSE << "bp_i   = " << bp_i << endl;
  logger << LOG_DEBUG_VERBOSE << "bp_j   = " << bp_j << endl;
  logger << LOG_DEBUG_VERBOSE << "bp_k   = " << bp_k << endl;
  logger << LOG_DEBUG_VERBOSE << "bpt_ij = " << bpt_ij << endl;
  logger << LOG_DEBUG_VERBOSE << "bpt_k  = " << bpt_k << endl;
  logger << LOG_DEBUG_VERBOSE << "p_i   = " << xbp_all[0][bp_i] << endl;
  logger << LOG_DEBUG_VERBOSE << "p_j   = " << xbp_all[0][bp_j] << endl;
  logger << LOG_DEBUG_VERBOSE << "p_k   = " << xbp_all[0][bp_k] << endl;
  */
  fourvector Q = p_i + p_j + p_k;
  double Q2 = Q.m2();
  double pi_pj = p_i * p_j;
  /*
  logger << LOG_DEBUG_VERBOSE << "Q2    = " << Q2 << endl;
  logger << LOG_DEBUG_VERBOSE << "m2_i = " << m2_i << endl;
  logger << LOG_DEBUG_VERBOSE << "m2_j = " << m2_j << endl;
  logger << LOG_DEBUG_VERBOSE << "m2_ij = " << m2_ij << endl;
  logger << LOG_DEBUG_VERBOSE << "m2_k  = " << m2_k << endl;
  logger << LOG_DEBUG_VERBOSE << "pi_pj = " << pi_pj << endl;
  logger << LOG_DEBUG_VERBOSE << "lambda(Q2, m2_ij, m2_k) = " << lambda(Q2, m2_ij, m2_k) << endl;
  logger << LOG_DEBUG_VERBOSE << "lambda(Q2, m2_i + m2_j + 2 * pi_pj, m2_k) = " << lambda(Q2, m2_i + m2_j + 2 * pi_pj, m2_k) << endl;
  */
  pt_k = sqrt(lambda(Q2, m2_ij, m2_k) / lambda(Q2, m2_i + m2_j + 2 * pi_pj, m2_k)) * (p_k - ((Q * p_k) / Q2) * Q) + ((Q2 + m2_k - m2_ij) / (2 * Q2)) * Q;
  pt_ij = Q - pt_k;
  /*
  logger << LOG_DEBUG_VERBOSE << "Q     = " << Q << endl;
  logger << LOG_DEBUG_VERBOSE << "pt_k  = " << pt_k << endl;
  logger << LOG_DEBUG_VERBOSE << "pt_ij = " << pt_ij << endl;
  */
  // uncommented !!!
  y_ij_k = (p_i * p_j) / (p_i * p_j + p_j * p_k + p_k * p_i);
  //  logger << LOG_DEBUG_VERBOSE << "y_ij_k = " << y_ij_k << endl;
  /*
  pt_ij = p_i + p_j - y_ij_k / (1. - y_ij_k) * p_k;
  pt_k = 1. / (1. - y_ij_k) * p_k;
  */  

  for (int xbi = 1; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else if (xbi == bp_j){}
    else if (xbi == bp_k){}
    else if (xbi < bp_j){xbp_all[i_a][xbi] = xbp_all[0][xbi];}
    else if (xbi > bp_j){xbp_all[i_a][xbi / 2] = xbp_all[0][xbi];}
    else {cout << "Should not happen!" << endl;}
  }
  /*
  for (int xbi = 1; xbi < xbp_all[i_a].size(); xbi = xbi * 2){
    logger << LOG_DEBUG_VERBOSE << "xbp_all[" << i_a << "][" << xbi << "] = " << xbp_all[i_a][xbi] << endl;
  }
  */
  xbp_all[i_a][0] = xbp_all[i_a][1] + xbp_all[i_a][2];
  xbs_all[i_a][0] = xbs_all[0][0];
  xbsqrts_all[i_a][0] = xbsqrts_all[0][0];
  int xbn_all = xbp_all[i_a].size() - 4;
  xbp_all[i_a][xbn_all] = xbp_all[i_a][0];
  xbs_all[i_a][xbn_all] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xbn_all] = xbsqrts_all[i_a][0];
  /*
  for (int i = 0; i < osi_p_parton[i_a].size(); i++){
    int xbi = intpow(2, i - 1);
    osi_p_parton[i_a][i] = xbp_all[i_a][xbi];
  }
  */
  logger << LOG_DEBUG_VERBOSE << "phasespace_ij_k finished" << endl;
}
#undef m2_i
#undef m2_j
#undef m2_k
#undef m2_ij
#undef p_i 
#undef p_j 
#undef p_k 
#undef pt_ij 
#undef pt_k 
#undef y_ij_k

#define p_i xbp_all[0][bp_i]
#define p_j xbp_all[0][bp_j]
#define p_a xbp_all[0][bp_a]
#define pt_ij xbp_all[i_a][bpt_ij]
#define pt_a xbp_all[i_a][bpt_a]
#define m2_i xbs_all[0][bp_i]
#define m2_j xbs_all[0][bp_j]
#define m2_ij xbs_all[i_a][bpt_ij]
#define x_ij_a dipole[i_a].xy
void phasespace_set::phasespace_ij_a_massive(int i_a, vector<dipole_set> & dipole){
  static Logger logger("phasespace_ij_a_massive");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  int bp_i = dipole[i_a].binary_R_emitter_1();
  int bp_j = dipole[i_a].binary_R_emitter_2();
  int bp_a = dipole[i_a].binary_R_spectator(); // bp_b not needed !!!
  int bpt_ij = dipole[i_a].binary_A_emitter();
  int bpt_a = dipole[i_a].binary_A_spectator(); // always bp_a == bpt_a

  double pi_pj = p_i * p_j;
  double pi_pa = p_i * p_a;
  double pj_pa = p_j * p_a;

  double one_minus_x_ij_a = (pi_pj - .5 * (m2_ij - m2_i - m2_j)) / (pi_pa + pj_pa);
  x_ij_a = 1. - one_minus_x_ij_a;
  pt_ij = p_i + p_j - one_minus_x_ij_a * p_a;
  pt_a = x_ij_a * p_a;

  /*
  x_ij_a = 1. - (p_i * p_j) / ((p_i + p_j) * p_a);
  pt_ij = p_i + p_j - (1 - x_ij_a) * p_a;
  pt_a = x_ij_a * p_a;
  */
  for (int xbi = 1; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else if (xbi == bp_j){}
    else if (xbi == bp_a){}
    else if (xbi < bp_j){xbp_all[i_a][xbi] = xbp_all[0][xbi];}
    else if (xbi > bp_j){xbp_all[i_a][xbi / 2] = xbp_all[0][xbi];}
    else {cout << "Should not happen!" << endl;}
  }

  xbp_all[i_a][0] = xbp_all[i_a][1] + xbp_all[i_a][2];
  xbs_all[i_a][0] = x_ij_a * xbs_all[0][0];
  xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);
  int xbn_all = xbp_all[i_a].size() - 4;
  xbp_all[i_a][xbn_all] = xbp_all[i_a][0];
  xbs_all[i_a][xbn_all] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xbn_all] = xbsqrts_all[i_a][0];
  /*
  for (int i = 0; i < osi_p_parton[i_a].size(); i++){
    int xbi = intpow(2, i - 1);
    osi_p_parton[i_a][i] = xbp_all[i_a][xbi];
  }
  */
  logger << LOG_DEBUG_VERBOSE << "phasespace_ij_k finished" << endl;
}
#undef x_ij_a
#undef m2_i
#undef m2_j
#undef m2_ij
#undef p_i 
#undef p_j 
#undef p_a 
#undef pt_ij 
#undef pt_a 



void phasespace_set::output_check_tau_0(){
  Logger logger("phasespace_set::output_check_tau_0");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (x_pdf[0] < tau_0_num[tau_0_num.size() - 1]){
    tau_0_num.push_back(x_pdf[0]);
    logger << LOG_INFO << "tau_0 = " << setw(23) << setprecision(15) << tau_0 << "   min(" << setw(3) << tau_0_num.size() - 1 << "): " << setw(23) << setprecision(15) << tau_0_num[tau_0_num.size() - 1] << endl;
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




void phasespace_set::handling_cut_psp(){
  static Logger logger("phasespace_set::handling_cut_psp");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (MC_tau.end_optimization == 0){
    MC_tau.n_rej++;
    MC_tau.n_rej_channel[MC_tau.channel]++;
  }

  // count cutted events for weight optimization of dipole mappings
  if (RA_x_a != 0) {
    if ((switch_MC_x_dipole != -1) && csi->class_contribution_CS_real){
      /*
    if ((switch_MC_x_dipole != -1) &&(type_contribution == "RA" || 
					  type_contribution == "RRA")){
      */
      if (MC_x_dipole[RA_x_a].end_optimization == 0){
        MC_x_dipole[RA_x_a].n_rej++;
        MC_x_dipole[RA_x_a].n_rej_channel[MC_x_dipole[RA_x_a].channel]++;
      }
    }
  }
  /*
  if (tau_opt_end == 0){tau_cuts_channel[tau_channel]++;}
  if (x1x2_opt_end == 0){x1x2_cuts_channel[x1x2_channel]++;}
  ///  cuts_channel[MC_phasespace.channel]++;
  */
  
  if (!IS_tau.end_optimization){
    IS_tau.n_rej++;
    IS_tau.n_rej_channel[IS_tau.channel]++;
  }
  
  if (!IS_x1x2.end_optimization){
    IS_x1x2.n_rej++;
    IS_x1x2.n_rej_channel[IS_x1x2.channel]++;
  }

  i_rej++;
  MC_phasespace.n_rej_channel[MC_phasespace.channel]++;
  
  random_manager.increase_cut_counter();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::handling_techcut_psp(){
  static Logger logger("phasespace_set::handling_techcut_psp");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (MC_tau.end_optimization == 0){
    MC_tau.n_acc++;
    MC_tau.n_acc_channel[MC_tau.channel]++;
  }

  // count cutted events for weight optimization of dipole mappings
  if (RA_x_a != 0) {
    if ((switch_MC_x_dipole != -1) && csi->class_contribution_CS_real){
      if (MC_x_dipole[RA_x_a].end_optimization == 0){
        MC_x_dipole[RA_x_a].n_acc++;
        MC_x_dipole[RA_x_a].n_acc_channel[MC_x_dipole[RA_x_a].channel]++;
      }
    }
  }

  if (!IS_tau.end_optimization){
    IS_tau.n_acc++;
    IS_tau.n_acc_channel[IS_tau.channel]++;
  }
  
  if (!IS_x1x2.end_optimization){
    IS_x1x2.n_acc++;
    IS_x1x2.n_acc_channel[IS_x1x2.channel]++;
  }

  i_acc++;
  MC_phasespace.n_acc_channel[MC_phasespace.channel]++;
  
  random_manager.increase_counter_n_acc();
  //  random_manager.increase_cut_counter();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

//#include "OV.phasespace.set.cpp"
//#include "OV2.phasespace.set.cpp"
