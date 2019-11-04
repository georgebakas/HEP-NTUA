#include "../include/classes.cxx"
//#include "../include/definitions.phasespace.set.cxx"

void phasespace_set::ac_psp_RA_group_ij_k(vector<dipole_set> & dipole, int i_a){
  static Logger logger("phasespace_set::ac_psp_RA_group_ij_k");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double g_IS_temp;
  if (switch_IS_mode_phasespace == 3 || switch_IS_mode_phasespace == 4){
    for (int j = 2; j < 3; j++){ // manually switched off x and z IS
      if (container_IS_switch[container_IS_startvalue[i_a][4 + j]] == 1){phasespace_randoms[container_IS_startvalue[i_a][4 + j]]->get_random(r[no_random_dipole[j]], g_IS_temp);}
      MC_g_IS_global *= g_IS_temp;
    }
  }
  for (int i = 0; i < 3; i++){xbp_all[i_a][i] = xbp_all[0][i];}
  xbs_all[i_a][0] = xbs_all[0][0];
  xbsqrts_all[i_a][0] = xbsqrts_all[0][0];
  xbp_all[i_a][xb_out_dipoles] = xbp_all[i_a][0];
  xbs_all[i_a][xb_out_dipoles] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xb_out_dipoles] = xbsqrts_all[i_a][0];
  generic->ac_psp_dipole(i_a, MC_phasespace.channel - dipole[i_a - 1].sum_channel(), *this);
  int bp_i = dipole[i_a].binary_R_emitter_1();
  int bp_j = dipole[i_a].binary_R_emitter_2();
  int bp_k = dipole[i_a].binary_R_spectator();
  int bpt_ij = dipole[i_a].binary_A_emitter();
  int bpt_k = dipole[i_a].binary_A_spectator();
  fourvector temp_p_ijk = xbp_all[i_a][bpt_ij] + xbp_all[i_a][bpt_k];
  fourvector temp_p_k_boost = xbp_all[i_a][bpt_k].boost(temp_p_ijk);
  double sqrt_p2_ijk = temp_p_ijk.m();
  double y_ij_k = h_propto_pot(r[no_random_dipole[0]], 0., 1., exp_ij_k_y, map_technical_x);
  double one_minus_y_ij_k = 1. - y_ij_k;
  double one_minus_z_i = h_propto_pot(r[no_random_dipole[1]], 0., 1., exp_ij_k_z, map_technical_x);
  double z_i = 1. - one_minus_z_i;
  double phi_k_perp = c_phi(no_random_dipole[2]);
  double f2 = y_ij_k * z_i / one_minus_z_i;
  double f = sqrt(f2);
  fourvector khat_perp(-f2 * sqrt_p2_ijk, f * sqrt_p2_ijk * cos(phi_k_perp), f * sqrt_p2_ijk * sin(phi_k_perp), -f2 * sqrt_p2_ijk);
  fourvector k_perp = (khat_perp.rotateback(temp_p_k_boost)).boost(temp_p_ijk.Pinv());

  xbp_all[0][bp_i] = y_ij_k * (1. + z_i) * xbp_all[i_a][bpt_k] + z_i  * xbp_all[i_a][bpt_ij] + one_minus_z_i * k_perp;
  xbp_all[0][bp_j] = one_minus_z_i * (xbp_all[i_a][bpt_ij] - k_perp) - y_ij_k * z_i * xbp_all[i_a][bpt_k];
  xbp_all[0][bp_k] = one_minus_y_ij_k * xbp_all[i_a][bpt_k];
  for (int xbi = 4; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else if (xbi == bp_j){}
    else if (xbi == bp_k){}
    else if (xbi <  bp_j){xbp_all[0][xbi] = xbp_all[i_a][xbi];}
    else if (xbi >  bp_j){xbp_all[0][xbi] = xbp_all[i_a][xbi / 2];}
    else {logger << LOG_INFO << "Should not happen!" << endl;}
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::ag_psp_RA_group_ij_k(vector<dipole_set> & dipole, int i_a){
  static Logger logger("phasespace_set::ag_psp_RA_group_ij_k");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int bp_i = dipole[i_a].binary_R_emitter_1();
  int bp_j = dipole[i_a].binary_R_emitter_2();
  int bp_k = dipole[i_a].binary_R_spectator();
  int bpt_ij = dipole[i_a].binary_A_emitter();
  int bpt_k = dipole[i_a].binary_A_spectator();
  if (xbs_all[0][bp_i + bp_j] == 0.){xbs_all[0][bp_i + bp_j] = (xbp_all[0][bp_i] + xbp_all[0][bp_j]).m2();}
  double p_kp_ij = 2. * (xbp_all[0][bp_k] * (xbp_all[0][bp_i] + xbp_all[0][bp_j]));
  double y_ij_k = 1. / (1. + p_kp_ij / xbs_all[0][bp_i + bp_j]);
  //  double one_minus_y_ij_k = 1. - y_ij_k;
  double z_i = (xbp_all[0][bp_i] * xbp_all[0][bp_k]) / (xbp_all[0][bp_i] * xbp_all[0][bp_k] + xbp_all[0][bp_j] * xbp_all[0][bp_k]);
  double one_minus_z_i = (xbp_all[0][bp_j] * xbp_all[0][bp_k]) / (xbp_all[0][bp_i] * xbp_all[0][bp_k] + xbp_all[0][bp_j] * xbp_all[0][bp_k]);
  double g_alpha = g_propto_pot(y_ij_k, 0., 1., exp_ij_k_y, map_technical_x) * g_propto_pot(one_minus_z_i, 0., 1., exp_ij_k_z, map_technical_x) / (pi * (xbp_all[i_a][bpt_ij] * xbp_all[i_a][bpt_k]) * (1. - y_ij_k));
  if (switch_console_output_phasespace_issue){
    if (g_alpha < 0.){logger << LOG_INFO << "i_a = " << i_a << "   g_alpha = " << g_alpha << endl;}
  }
  if (xbs_all[i_a][0] == 0.){xbs_all[i_a][0] = xbs_all[0][0];}
  if (xbsqrts_all[i_a][0] == 0.){xbsqrts_all[i_a][0] = xbsqrts_all[0][0];}
  generic->ag_psp_dipole(i_a, dipole[i_a - 1].sum_channel(), *this);
  for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){MC_phasespace.g_channel[i] = MC_phasespace.g_channel[i] * g_alpha;}
  // invert mappings for VAMP implementation
  if (switch_IS_mode_phasespace == 2 || switch_IS_mode_phasespace == 4){
    //  static vector<double> inv_r(3);
    vector<double> inv_r(3);
    inv_r[0] = inv_propto_pot(y_ij_k, 0., 1., exp_ij_k_y, map_technical_x);
    inv_r[1] = inv_propto_pot(one_minus_z_i, 0., 1., exp_ij_k_z, map_technical_x);
    fourvector k_perp = xbp_all[0][bp_j] + y_ij_k * z_i * xbp_all[i_a][bpt_k];
    k_perp = -(k_perp - xbp_all[i_a][bpt_ij])/one_minus_z_i;
    fourvector temp_p_ijk = xbp_all[i_a][bpt_ij] + xbp_all[i_a][bpt_k];
    fourvector temp_p_k_boost = xbp_all[i_a][bpt_k].boost(temp_p_ijk);
    k_perp = (k_perp.boost(temp_p_ijk)).rotateback_inverse(temp_p_k_boost);
    double phi = atan2(k_perp.x2(),k_perp.x1());
    if (phi < 0){phi += 2 * M_PI;}
    inv_r[2] = inv_phi(phi);
    if (!(inv_r[2] >= 0) || inv_r[2] > 1){
      if (switch_console_output_phasespace_issue){
	logger << LOG_INFO << "inv_r[2] not in [0; 1] !   r = " << inv_r[2] <<"; k_perp.x,y,z=" << k_perp.x1() << ", " << k_perp.x2() << ", " << k_perp.x3() << endl;
      }
      if (!(inv_r[2] >= 0)){inv_r[2] = 0.;}
      if (inv_r[2] > 1){inv_r[2] = 1.;}
    }
    if (switch_IS_mode_phasespace == 2){
      for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){
        double g_IS_temp;
        double g_IS = 1.;
        for (int j = 0; j < 3; j++){
          phasespace_randoms[(i + 1) * 8 - 3 + j]->get_g_IS(inv_r[j], g_IS_temp);
          g_IS *= g_IS_temp;
        }
        MC_phasespace.g_channel[i] *= g_IS;
      }
    }
    else if (switch_IS_mode_phasespace == 4){
      double g_IS_temp;
      double g_IS = 1.;
      for (int j = 2; j < 3; j++){ // manually switched of x and z IS
        phasespace_randoms[container_IS_startvalue[i_a][4 + j]]->get_g_IS(inv_r[j], g_IS_temp);
        g_IS *= g_IS_temp;
      }
      for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){MC_phasespace.g_channel[i] *= g_IS;}
    }
    //  logger << LOG_INFO << "reconstr " << i_a << ": " << inv_r[0] << ", " << inv_r[1] << ", " << inv_r[2] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::ac_psp_RA_group_ij_a(vector<dipole_set> & dipole, double & sinx_ij_a_min, int i_a){
  static Logger logger("phasespace_set::ac_psp_RA_group_ij_a");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double g_IS_temp;
  if (switch_IS_mode_phasespace == 3 || switch_IS_mode_phasespace == 4){
    for (int j = 2; j < 3; j++){ // manually switched off x and z IS
      if (container_IS_switch[container_IS_startvalue[i_a][4 + j]] == 1){phasespace_randoms[container_IS_startvalue[i_a][4 + j]]->get_random(r[no_random_dipole[j]], g_IS_temp);}
      MC_g_IS_global *= g_IS_temp;
    }
  }
  int bp_i = dipole[i_a].binary_R_emitter_1();
  int bp_j = dipole[i_a].binary_R_emitter_2();
  int bp_a = dipole[i_a].binary_R_spectator();
  int bpt_ij = dipole[i_a].binary_A_emitter();
  //  int bpt_a = dipole[i_a].binary_A_spectator(); // always bp_a == bpt_a
  int bp_b = bp_a % 2 + 1;

  /*
  double one_minus_x_ij_a_min = 1. - sinx_ij_a_min / xbs_all[0][0];
  double one_minus_x_ij_a = h_propto_pot(r[no_random_dipole[0]], 0., one_minus_x_ij_a_min, exp_ij_a_x, map_technical_x);
  double x_ij_a = 1. - one_minus_x_ij_a;
  */
  // multi-channel, where channel 0 reproduces the original version, and channels >0 map potential IS resonances...
  double one_minus_x_ij_a_min = 1. - sinx_ij_a_min / xbs_all[0][0];
  double x_ij_a = 1.;
  double one_minus_x_ij_a = 0.;
  // random_MC_tau[1] -> increase number in calculate_intial_tau_IS_x1x2_IS !!!
  for (int i_c = 0; i_c < MC_x_dipole[i_a].beta.size(); i_c++){
    logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[" << i_a << "].beta[" << i_c << "] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].beta[i_c] << endl;
  }
  for (int i_c = 0; i_c < MC_x_dipole[i_a].beta.size(); i_c++){if (random_MC_tau[1] <= MC_x_dipole[i_a].beta[i_c]){MC_x_dipole[i_a].channel = i_c; break;}}
  logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[" << i_a << "].channel = " << MC_x_dipole[i_a].channel << endl;

  if (MC_x_dipole[i_a].channel == 0){
    one_minus_x_ij_a = h_propto_pot(r[no_random_dipole[0]], 0., one_minus_x_ij_a_min, exp_ij_a_x, map_technical_x);
    x_ij_a = 1. - one_minus_x_ij_a;
  }
  else {
    if (MC_x_dipole_mapping[i_a][MC_x_dipole[i_a].channel] < 0){
      x_ij_a = c_propagator_Breit_Wigner(r[no_random_dipole[0]], -MC_x_dipole_mapping[i_a][MC_x_dipole[i_a].channel], sinx_ij_a_min, xbs_all[0][0]) / xbs_all[0][0];
      one_minus_x_ij_a = 1. - x_ij_a;
      logger << LOG_DEBUG_VERBOSE << "s_dipole = " << setw(23) << setprecision(15) << sqrt(x_ij_a * xbs_all[0][0]) << "   " << "x_ij_a = " << setw(23) << setprecision(15) << x_ij_a << "   " << "xbs_all[0][0] = " << setw(23) << setprecision(15) << xbs_all[0][0] << endl;
    }
  }


  xbp_all[i_a][bp_a] = x_ij_a * xbp_all[0][bp_a];
  xbp_all[i_a][bp_b] = xbp_all[0][bp_b];
  xbp_all[i_a][0] = xbp_all[i_a][bp_a] + xbp_all[i_a][bp_b];
  xbs_all[i_a][0] = x_ij_a * xbs_all[0][0];
  xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);
  xbp_all[i_a][xb_out_dipoles] = xbp_all[i_a][0];
  xbs_all[i_a][xb_out_dipoles] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xb_out_dipoles] = xbsqrts_all[i_a][0];
  generic->ac_psp_dipole(i_a, MC_phasespace.channel - dipole[i_a - 1].sum_channel(), *this);
  fourvector temp_p_ija = xbp_all[i_a][bpt_ij] + xbp_all[i_a][bp_a];
  fourvector temp_p_a_boost = xbp_all[i_a][bp_a].boost(temp_p_ija);
  double sqrt_p2_ija = temp_p_ija.m();
  double one_minus_z_i = h_propto_pot(r[no_random_dipole[1]], 0., 1., exp_ij_a_z, map_technical_x);
  double z_i = 1. - one_minus_z_i;
  double phi_a_perp = c_phi(no_random_dipole[2]);
  double f2 = z_i * one_minus_x_ij_a / (one_minus_z_i * x_ij_a);
  double f = sqrt(f2);
  fourvector khat_perp(-f2 * sqrt_p2_ija, f * sqrt_p2_ija * cos(phi_a_perp), f * sqrt_p2_ija * sin(phi_a_perp), -f2 * sqrt_p2_ija);
  fourvector k_perp = (khat_perp.rotateback(temp_p_a_boost)).boost(temp_p_ija.Pinv());

  xbp_all[0][bp_i] = one_minus_x_ij_a / x_ij_a * (1. + z_i) * xbp_all[i_a][bp_a] + z_i * xbp_all[i_a][bpt_ij] + one_minus_z_i * k_perp;
  xbp_all[0][bp_j] = one_minus_z_i * (xbp_all[i_a][bpt_ij] - k_perp) - one_minus_x_ij_a / x_ij_a * z_i * xbp_all[i_a][bp_a];
  xbp_all[0][bp_a] = xbp_all[i_a][bp_a] / x_ij_a;
  for (int xbi = 4; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else if (xbi == bp_j){}
    else if (xbi < bp_j){xbp_all[0][xbi] = xbp_all[i_a][xbi];}
    else if (xbi > bp_j){xbp_all[0][xbi] = xbp_all[i_a][xbi / 2];}
    else {logger << LOG_INFO << "Should not happen!" << endl;}
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::ag_psp_RA_group_ij_a(double & sinx_ij_a_min, vector<dipole_set> & dipole, int i_a){
  static Logger logger("phasespace_set::ag_psp_RA_group_ij_a");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  //  static int xb_all_out = intpow(2, csi->n_particle + 1) - 4;
  int bp_i = dipole[i_a].binary_R_emitter_1();
  int bp_j = dipole[i_a].binary_R_emitter_2();
  int bp_a = dipole[i_a].binary_R_spectator(); // bp_b not needed !!!
  int bpt_ij = dipole[i_a].binary_A_emitter();
  //  int bpt_a = dipole[i_a].binary_A_spectator(); // always bp_a == bpt_a
  if (xbs_all[0][bp_i + bp_j] == 0.){xbs_all[0][bp_i + bp_j] = (xbp_all[0][bp_i] + xbp_all[0][bp_j]).m2();}
  double pa_p12 = 2. * (xbp_all[0][bp_a] * (xbp_all[0][bp_i] + xbp_all[0][bp_j]));

  double one_minus_x_ij_a_min = 1. - sinx_ij_a_min / xbs_all[0][0];
  double one_minus_x_ij_a = xbs_all[0][bp_i + bp_j] / pa_p12;
  double x_ij_a = 1. - one_minus_x_ij_a;

  // multi-channel, where channel 0 reproduces the original version, and channels >0 map potential IS resonances...
  double g_x_ij_a = g_propto_pot(one_minus_x_ij_a, 0., one_minus_x_ij_a_min, exp_ij_a_x, map_technical_x);
  double x_ij_a_s_part = x_ij_a * xbs_all[0][0];
  logger << LOG_DEBUG_VERBOSE << "x_ij_a_s_part = " << x_ij_a_s_part << "   sqrt -> " << sqrt(x_ij_a_s_part) << endl;
  for (int i_c = 0; i_c < MC_x_dipole[i_a].g_channel.size(); i_c++){
    if (i_c == 0){
      MC_x_dipole[i_a].g_channel[i_c] = g_propto_pot_mod(one_minus_x_ij_a, 0., one_minus_x_ij_a_min, exp_ij_a_x, map_technical_x);
    }
    else {
      if (MC_x_dipole_mapping[i_a][i_c] < 0){MC_x_dipole[i_a].g_channel[i_c] = xbs_all[0][0] * g_propagator_Breit_Wigner(x_ij_a_s_part, -MC_x_dipole_mapping[i_a][i_c], sinx_ij_a_min, xbs_all[0][0]);}
      logger << LOG_DEBUG_VERBOSE << "sinx_ij_a_min = " << setw(23) << setprecision(15) << sinx_ij_a_min << endl;
      logger << LOG_DEBUG_VERBOSE << "x_ij_a_s_part = " << setw(23) << setprecision(15) << x_ij_a_s_part << endl;
      logger << LOG_DEBUG_VERBOSE << "xbs_all[0][0] = " << setw(23) << setprecision(15) << xbs_all[0][0] << endl;
      logger << LOG_DEBUG_VERBOSE << "-MC_x_dipole_mapping[i_a][i_c] = " << setw(23) << setprecision(15) << -MC_x_dipole_mapping[i_a][MC_x_dipole[i_a].channel] << endl;
      logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[i_a].g_channel[i_c] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].g_channel[i_c] << endl;
      ///      else if (tau_MC_map[i_c] > 0){MC_x_dipole[i_a].g_channel[i_c] = s_had * g_propagator_vanishing_width(x_ij_a_s_part, M2[tau_MC_map[i_c]], tau_0_s_had, s_had, nuxs);}
    }
  }
  g_x_ij_a = 0.;
  for (int i = 0; i < MC_x_dipole[i_a].g_channel.size(); i++){
    logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[" << i_a << "].alpha[" << i << "] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].alpha[i] << "   " << "MC_x_dipole[" << i_a << "].g_channel[" << i << "] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].g_channel[i] << endl;
  }
  for (int i = 0; i < MC_x_dipole[i_a].g_channel.size(); i++){g_x_ij_a += MC_x_dipole[i_a].alpha[i] * MC_x_dipole[i_a].g_channel[i];}
  logger << LOG_DEBUG_VERBOSE << "g_x_ij_a = " << g_x_ij_a << endl;

  double z_i = (xbp_all[0][bp_i] * xbp_all[0][bp_a]) / (xbp_all[0][bp_i] * xbp_all[0][bp_a] + xbp_all[0][bp_j] * xbp_all[0][bp_a]);
  double one_minus_z_i = (xbp_all[0][bp_j] * xbp_all[0][bp_a]) / (xbp_all[0][bp_i] * xbp_all[0][bp_a] + xbp_all[0][bp_j] * xbp_all[0][bp_a]);
  double g_alpha = g_x_ij_a * g_propto_pot(one_minus_z_i, 0., 1., exp_ij_a_z, map_technical_x) / (pi * (xbp_all[i_a][bpt_ij] * xbp_all[0][bp_a]));
  //  double g_alpha = g_propto_pot(one_minus_x_ij_a, 0., one_minus_x_ij_a_min, exp_ij_a_x, map_technical_x) * g_propto_pot(one_minus_z_i, 0., 1., exp_ij_a_z, map_technical_x) / (pi * (xbp_all[i_a][bpt_ij] * xbp_all[0][bp_a]));
  if (switch_console_output_phasespace_issue){
    if (g_alpha < 0.){logger << LOG_INFO << "i_a = " << i_a << "   g_alpha = " << g_alpha << endl;}
  }
  if (xbs_all[i_a][0] == 0.){xbs_all[i_a][0] = x_ij_a * xbs_all[0][0];}
  if (xbsqrts_all[i_a][0] == 0.){xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);}
  generic->ag_psp_dipole(i_a, dipole[i_a - 1].sum_channel(), *this);
  for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){MC_phasespace.g_channel[i] = MC_phasespace.g_channel[i] * g_alpha;}
  // invert mappings for VAMP implementation
  if (switch_IS_mode_phasespace == 2 || switch_IS_mode_phasespace == 4){
    //  static vector<double> inv_r(3);
    vector<double> inv_r(3);
    inv_r[0] = inv_propto_pot(one_minus_x_ij_a, 0., one_minus_x_ij_a_min, exp_ij_a_x, map_technical_x);
    inv_r[1] = inv_propto_pot(one_minus_z_i, 0., 1., exp_ij_a_z, map_technical_x);
    fourvector k_perp = xbp_all[0][bp_j] + one_minus_x_ij_a / x_ij_a * z_i * xbp_all[i_a][bp_a];
    k_perp = -(k_perp - xbp_all[i_a][bpt_ij])/one_minus_z_i;
    fourvector temp_p_ija = xbp_all[i_a][bpt_ij] + xbp_all[i_a][bp_a];
    fourvector temp_p_a_boost = xbp_all[i_a][bp_a].boost(temp_p_ija);
    k_perp = (k_perp.boost(temp_p_ija)).rotateback_inverse(temp_p_a_boost);
    double phi = atan2(k_perp.x2(),k_perp.x1());
    if (phi < 0){phi += 2 * M_PI;}
    inv_r[2] = inv_phi(phi);
    if (!(inv_r[2] >= 0) || inv_r[2] > 1){
      if (switch_console_output_phasespace_issue){
	logger << LOG_INFO << "inv_r[2] not in [0; 1] !   r = " << inv_r[2] <<"; k_perp.x,y,z=" << k_perp.x1() << ", " << k_perp.x2() << ", " << k_perp.x3() << endl;
      }
      if (!(inv_r[2] >= 0)){inv_r[2] = 0.;}
      if (inv_r[2] > 1){inv_r[2] = 1.;}
    }
    if (switch_IS_mode_phasespace == 2){
      for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){
        double g_IS_temp;
        double g_IS = 1.;
        for (int j = 0; j < 3; j++){
          phasespace_randoms[(i + 1) * 8 - 3 + j]->get_g_IS(inv_r[j], g_IS_temp);
          g_IS *= g_IS_temp;
        }
        MC_phasespace.g_channel[i] *= g_IS;
      }
    }
    else if (switch_IS_mode_phasespace == 4){
      double g_IS_temp;
      double g_IS = 1.;
      for (int j = 2; j < 3; j++){ // manually switched of x and z IS
        phasespace_randoms[container_IS_startvalue[i_a][4 + j]]->get_g_IS(inv_r[j], g_IS_temp);
        g_IS *= g_IS_temp;
      }
      for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){MC_phasespace.g_channel[i] *= g_IS;}
    }
    //  logger << LOG_INFO << "reconstr " << i_a << ": " << inv_r[0] << ", " << inv_r[1] << ", " << inv_r[2] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::ac_psp_RA_group_ai_k(vector<dipole_set> & dipole, double & sinx_ik_a_min, int i_a){
  static Logger logger("phasespace_set::ac_psp_RA_group_ai_k");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double g_IS_temp;
  if (switch_IS_mode_phasespace == 3 || switch_IS_mode_phasespace == 4){
    for (int j = 2; j < 3; j++){ // manually switched of x and z IS
      if (container_IS_switch[container_IS_startvalue[i_a][4 + j]] == 1){phasespace_randoms[container_IS_startvalue[i_a][4 + j]]->get_random(r[no_random_dipole[j]], g_IS_temp);}
      MC_g_IS_global *= g_IS_temp;
    }
  }
  int bp_a = dipole[i_a].binary_R_emitter_1();
  int bp_i = dipole[i_a].binary_R_emitter_2();
  int bp_k = dipole[i_a].binary_R_spectator();
  //  int bpt_ai = dipole[i_a].binary_A_emitter();
  int bpt_k = dipole[i_a].binary_A_spectator();
  int bp_b = bp_a % 2 + 1;



  /*
  double one_minus_x_ik_a_min = 1. - sinx_ik_a_min / xbs_all[0][0];
  double one_minus_x_ik_a = h_propto_pot(r[no_random_dipole[0]], 0., one_minus_x_ik_a_min, exp_ai_k_x, map_technical_x);
  double x_ik_a = 1. - one_minus_x_ik_a;
  */
  // multi-channel, where channel 0 reproduces the original version, and channels >0 map potential IS resonances...
  double one_minus_x_ik_a_min = 1. - sinx_ik_a_min / xbs_all[0][0];
  double x_ik_a = 1.;
  double one_minus_x_ik_a = 0.;
  // random_MC_tau[1] -> increase number in calculate_intial_tau_IS_x1x2_IS !!!
  for (int i_c = 0; i_c < MC_x_dipole[i_a].beta.size(); i_c++){
    logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[" << i_a << "].beta[" << i_c << "] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].beta[i_c] << endl;
  }
  for (int i_c = 0; i_c < MC_x_dipole[i_a].beta.size(); i_c++){if (random_MC_tau[1] <= MC_x_dipole[i_a].beta[i_c]){MC_x_dipole[i_a].channel = i_c; break;}}
  logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[" << i_a << "].channel = " << MC_x_dipole[i_a].channel << endl;

  if (MC_x_dipole[i_a].channel == 0){
    one_minus_x_ik_a = h_propto_pot(r[no_random_dipole[0]], 0., one_minus_x_ik_a_min, exp_ai_k_x, map_technical_x);
    x_ik_a = 1. - one_minus_x_ik_a;
  }
  else {
    if (MC_x_dipole_mapping[i_a][MC_x_dipole[i_a].channel] < 0){
      x_ik_a = c_propagator_Breit_Wigner(r[no_random_dipole[0]], -MC_x_dipole_mapping[i_a][MC_x_dipole[i_a].channel], sinx_ik_a_min, xbs_all[0][0]) / xbs_all[0][0];
      one_minus_x_ik_a = 1. - x_ik_a;
      logger << LOG_DEBUG_VERBOSE << "s_dipole = " << setw(23) << setprecision(15) << sqrt(x_ik_a * xbs_all[0][0]) << "   " << "x_ik_a = " << setw(23) << setprecision(15) << x_ik_a << "   " << "xbs_all[0][0] = " << setw(23) << setprecision(15) << xbs_all[0][0] << endl;
    }
  }
  
  xbp_all[i_a][bp_a] = x_ik_a * xbp_all[0][bp_a];
  xbp_all[i_a][bp_b] = xbp_all[0][bp_b];
  xbp_all[i_a][0] = xbp_all[i_a][bp_a] + xbp_all[i_a][bp_b];
  xbs_all[i_a][0] = x_ik_a * xbs_all[0][0];
  xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);
  xbp_all[i_a][xb_out_dipoles] = xbp_all[i_a][0];
  xbs_all[i_a][xb_out_dipoles] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xb_out_dipoles] = xbsqrts_all[i_a][0];
  generic->ac_psp_dipole(i_a, MC_phasespace.channel - dipole[i_a - 1].sum_channel(), *this);
  fourvector temp_p_aik = xbp_all[i_a][bp_a] + xbp_all[i_a][bpt_k];
  fourvector temp_p_k_boost = xbp_all[i_a][bp_a].boost(temp_p_aik);
  double sqrt_p2_aik = temp_p_aik.m();
  double u_i = h_propto_pot(r[no_random_dipole[1]], 0., 1., exp_ai_k_u, map_technical_x);
  double one_minus_u_i = 1. - u_i;
  double phi_k_perp = c_phi(no_random_dipole[2]);
  double f2 = u_i * one_minus_x_ik_a / (one_minus_u_i * x_ik_a);
  double f = sqrt(f2);
  fourvector khat_perp(-f2 * sqrt_p2_aik, f * sqrt_p2_aik * cos(phi_k_perp), f * sqrt_p2_aik * sin(phi_k_perp), -f2 * sqrt_p2_aik);
  fourvector k_perp = (khat_perp.rotateback(temp_p_k_boost)).boost(temp_p_aik.Pinv());

  xbp_all[0][bp_i] = one_minus_x_ik_a / x_ik_a * (1. + u_i) * xbp_all[i_a][bp_a] + u_i * xbp_all[i_a][bpt_k] + one_minus_u_i * k_perp;
  xbp_all[0][bp_k] = one_minus_u_i * (xbp_all[i_a][bpt_k] - k_perp) - one_minus_x_ik_a / x_ik_a * u_i * xbp_all[i_a][bp_a];
  xbp_all[0][bp_a] = xbp_all[i_a][bp_a] / x_ik_a;
  for (int xbi = 4; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else if (xbi == bp_k){}
    else if (xbi < bp_i){xbp_all[0][xbi] = xbp_all[i_a][xbi];}
    else if (xbi > bp_i){xbp_all[0][xbi] = xbp_all[i_a][xbi / 2];}
    else {logger << LOG_INFO << "Should not happen!" << endl;}
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::ag_psp_RA_group_ai_k(double & sinx_ik_a_min, vector<dipole_set> & dipole, int i_a){
  static Logger logger("phasespace_set::ag_psp_RA_group_ai_k");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int bp_a = dipole[i_a].binary_R_emitter_1();
  int bp_i = dipole[i_a].binary_R_emitter_2();
  int bp_k = dipole[i_a].binary_R_spectator();
  //  int bpt_ai = dipole[i_a].binary_A_emitter();
  int bpt_k = dipole[i_a].binary_A_spectator();
  //  int bp_b = bp_a % 2 + 1;
  if (xbs_all[0][bp_i + bp_k] == 0.){xbs_all[0][bp_i + bp_k] = (xbp_all[0][bp_i] + xbp_all[0][bp_k]).m2();}
  double pa_p12 = 2. * (xbp_all[0][bp_a] * (xbp_all[0][bp_i] + xbp_all[0][bp_k]));

  double one_minus_x_ik_a_min = 1. - sinx_ik_a_min / xbs_all[0][0];
  double one_minus_x_ik_a = xbs_all[0][bp_i + bp_k] / pa_p12;
  double x_ik_a = 1. - one_minus_x_ik_a;

  // multi-channel, where channel 0 reproduces the original version, and channels >0 map potential IS resonances...
  double g_x_ik_a = g_propto_pot(one_minus_x_ik_a, 0., one_minus_x_ik_a_min, exp_ai_k_x, map_technical_x);
  double x_ik_a_s_part = x_ik_a * xbs_all[0][0];
  logger << LOG_DEBUG_VERBOSE << "x_ik_a_s_part = " << x_ik_a_s_part << "   sqrt -> " << sqrt(x_ik_a_s_part) << endl;
  for (int i_c = 0; i_c < MC_x_dipole[i_a].g_channel.size(); i_c++){
    if (i_c == 0){
      MC_x_dipole[i_a].g_channel[i_c] = g_propto_pot_mod(one_minus_x_ik_a, 0., one_minus_x_ik_a_min, exp_ai_k_x, map_technical_x);
    }
    else {
      if (MC_x_dipole_mapping[i_a][i_c] < 0){MC_x_dipole[i_a].g_channel[i_c] = xbs_all[0][0] * g_propagator_Breit_Wigner(x_ik_a_s_part, -MC_x_dipole_mapping[i_a][i_c], sinx_ik_a_min, xbs_all[0][0]);}
      logger << LOG_DEBUG_VERBOSE << "sinx_ik_a_min = " << setw(23) << setprecision(15) << sinx_ik_a_min << endl;
      logger << LOG_DEBUG_VERBOSE << "x_ik_a_s_part = " << setw(23) << setprecision(15) << x_ik_a_s_part << endl;
      logger << LOG_DEBUG_VERBOSE << "xbs_all[0][0] = " << setw(23) << setprecision(15) << xbs_all[0][0] << endl;
      logger << LOG_DEBUG_VERBOSE << "-MC_x_dipole_mapping[i_a][i_c] = " << setw(23) << setprecision(15) << -MC_x_dipole_mapping[i_a][MC_x_dipole[i_a].channel] << endl;
      logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[i_a].g_channel[i_c] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].g_channel[i_c] << endl;
      ///      else if (tau_MC_map[i_c] > 0){MC_x_dipole[i_a].g_channel[i_c] = s_had * g_propagator_vanishing_width(x_ik_a_s_part, M2[tau_MC_map[i_c]], tau_0_s_had, s_had, nuxs);}
    }
  }
  g_x_ik_a = 0.;
  for (int i = 0; i < MC_x_dipole[i_a].g_channel.size(); i++){
    logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[" << i_a << "].alpha[" << i << "] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].alpha[i] << "   " << "MC_x_dipole[" << i_a << "].g_channel[" << i << "] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].g_channel[i] << endl;
  }
  for (int i = 0; i < MC_x_dipole[i_a].g_channel.size(); i++){g_x_ik_a += MC_x_dipole[i_a].alpha[i] * MC_x_dipole[i_a].g_channel[i];}
  logger << LOG_DEBUG_VERBOSE << "g_x_ik_a = " << g_x_ik_a << endl;

  double u_i = (xbp_all[0][bp_i] * xbp_all[0][bp_a]) / (xbp_all[0][bp_i] * xbp_all[0][bp_a] + xbp_all[0][bp_k] * xbp_all[0][bp_a]);
  double one_minus_u_i = 1. - u_i;
  double g_alpha = g_x_ik_a * g_propto_pot(u_i, 0., 1., exp_ai_k_u, map_technical_x) / (pi * (xbp_all[i_a][bpt_k] * xbp_all[0][bp_a]));
  if (switch_console_output_phasespace_issue){
    if (g_alpha < 0.){logger << LOG_INFO << "i_a = " << i_a << "   g_alpha = " << g_alpha << endl;}
  }
  //  double g_alpha = g_propto_pot(one_minus_x_ik_a, 0., one_minus_x_ik_a_min, exp_ai_k_x, map_technical_x) * g_propto_pot(u_i, 0., 1., exp_ai_k_u, map_technical_x) / (pi * (xbp_all[i_a][bpt_k] * xbp_all[0][bp_a]));
  if (xbs_all[i_a][0] == 0.){xbs_all[i_a][0] = x_ik_a * xbs_all[0][0];}
  if (xbsqrts_all[i_a][0] == 0.){xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);}
  generic->ag_psp_dipole(i_a, dipole[i_a - 1].sum_channel(), *this);
  for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){MC_phasespace.g_channel[i] = MC_phasespace.g_channel[i] * g_alpha;}
  // invert mappings for VAMP implementation
  if (switch_IS_mode_phasespace == 2 || switch_IS_mode_phasespace == 4){
    //  static vector<double> inv_r(3);
    vector<double> inv_r(3);
    inv_r[0] = inv_propto_pot(one_minus_x_ik_a, 0., one_minus_x_ik_a_min, exp_ai_k_x, map_technical_x);
    // !! extension needed if IS for x_ik_a is activated !!!
    inv_r[1] = inv_propto_pot(u_i, 0., 1., exp_ai_k_u, map_technical_x);
    fourvector k_perp = xbp_all[0][bp_k] + one_minus_x_ik_a / x_ik_a * u_i * xbp_all[i_a][bp_a];
    k_perp = -(k_perp - xbp_all[i_a][bpt_k])/one_minus_u_i;
    fourvector temp_p_aik = xbp_all[i_a][bp_a] + xbp_all[i_a][bpt_k];
    fourvector temp_p_k_boost = xbp_all[i_a][bp_a].boost(temp_p_aik);
    k_perp = (k_perp.boost(temp_p_aik)).rotateback_inverse(temp_p_k_boost);
    double phi = atan2(k_perp.x2(),k_perp.x1());
    if (phi < 0){phi += 2 * M_PI;}
    inv_r[2] = inv_phi(phi);
    if (!(inv_r[2] >= 0) || inv_r[2] > 1){
      if (switch_console_output_phasespace_issue){
	logger << LOG_INFO << "inv_r[2] not in [0; 1] !   r = " << inv_r[2] <<"; k_perp.x,y,z=" << k_perp.x1() << ", " << k_perp.x2() << ", " << k_perp.x3() << endl;
      }
      if (!(inv_r[2] >= 0)){inv_r[2] = 0.;}
      if (inv_r[2] > 1){inv_r[2] = 1.;}
    }
    if (switch_IS_mode_phasespace == 2){
      for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){
        double g_IS_temp;
        double g_IS = 1.;
        for (int j = 0; j < 3; j++){
          phasespace_randoms[(i + 1) * 8 - 3 + j]->get_g_IS(inv_r[j], g_IS_temp);
          g_IS *= g_IS_temp;
        }
        MC_phasespace.g_channel[i] *= g_IS;
      }
    }
    else if (switch_IS_mode_phasespace == 4){
      double g_IS_temp;
      double g_IS = 1.;
      for (int j = 2; j < 3; j++){ // manually switched of x and z IS
        phasespace_randoms[container_IS_startvalue[i_a][4 + j]]->get_g_IS(inv_r[j], g_IS_temp);
        g_IS *= g_IS_temp;
      }
      for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){MC_phasespace.g_channel[i] *= g_IS;}
    }
    //  logger << LOG_INFO << "reconstr " << i_a << ": " << inv_r[0] << ", " << inv_r[1] << ", " << inv_r[2] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::ac_psp_RA_group_ai_b(vector<dipole_set> & dipole, double & sinx_i_ab_min, int i_a){
  static Logger logger("phasespace_set::ac_psp_RA_group_ai_b");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double g_IS_temp;
  if (switch_IS_mode_phasespace == 3 || switch_IS_mode_phasespace == 4){
    for (int j = 2; j < 3; j++){ // manually switched off x and z IS
      if (container_IS_switch[container_IS_startvalue[i_a][4 + j]] == 1){phasespace_randoms[container_IS_startvalue[i_a][4 + j]]->get_random(r[no_random_dipole[j]], g_IS_temp);}
      MC_g_IS_global *= g_IS_temp;
    }
  }
  int bp_a = dipole[i_a].binary_R_emitter_1();
  int bp_i = dipole[i_a].binary_R_emitter_2();
  int bp_b = dipole[i_a].binary_R_spectator();
  //  int bpt_ai = dipole[i_a].binary_A_emitter(); // always bpt_ai == bp_a
  //  int bpt_b = dipole[i_a].binary_A_spectator(); // always bpt_b == bp_b


  /*
  double one_minus_x_i_ab_min = 1. - sinx_i_ab_min / xbs_all[0][0];
  double one_minus_x_i_ab = h_propto_pot_mod(r[no_random_dipole[0]], 0., one_minus_x_i_ab_min, exp_ai_b_x, map_technical_x);
  double x_i_ab = 1. - one_minus_x_i_ab;
  */
  // multi-channel, where channel 0 reproduces the original version, and channels >0 map potential IS resonances...
  double one_minus_x_i_ab_min = 1. - sinx_i_ab_min / xbs_all[0][0];
  double x_i_ab = 1.;
  double one_minus_x_i_ab = 0.;
  // random_MC_tau[1] -> increase number in calculate_intial_tau_IS_x1x2_IS !!!
  for (int i_c = 0; i_c < MC_x_dipole[i_a].beta.size(); i_c++){
    logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[" << i_a << "].beta[" << i_c << "] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].beta[i_c] << endl;
  }
  for (int i_c = 0; i_c < MC_x_dipole[i_a].beta.size(); i_c++){if (random_MC_tau[1] <= MC_x_dipole[i_a].beta[i_c]){MC_x_dipole[i_a].channel = i_c; break;}}
  logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[" << i_a << "].channel = " << MC_x_dipole[i_a].channel << endl;

  if (MC_x_dipole[i_a].channel == 0){
    one_minus_x_i_ab = h_propto_pot_mod(r[no_random_dipole[0]], 0., one_minus_x_i_ab_min, exp_ai_b_x, map_technical_x);
    x_i_ab = 1. - one_minus_x_i_ab;
  }
  else {
    if (MC_x_dipole_mapping[i_a][MC_x_dipole[i_a].channel] < 0){
      //      double sinx_i_ab_min_s_part = sinx_i_ab_min / xbs_all[0][0];
      //      x_pdf[0]
      x_i_ab = c_propagator_Breit_Wigner(r[no_random_dipole[0]], -MC_x_dipole_mapping[i_a][MC_x_dipole[i_a].channel], sinx_i_ab_min, xbs_all[0][0]) / xbs_all[0][0];
      one_minus_x_i_ab = 1. - x_i_ab;
      logger << LOG_DEBUG_VERBOSE << "s_dipole = " << setw(23) << setprecision(15) << sqrt(x_i_ab * xbs_all[0][0]) << "   " << "x_i_ab = " << setw(23) << setprecision(15) << x_i_ab << "   " << "xbs_all[0][0] = " << setw(23) << setprecision(15) << xbs_all[0][0] << endl;
    }
  }
  logger << LOG_DEBUG_VERBOSE << "sqrt(sinx_i_ab_min) = " << sqrt(sinx_i_ab_min) << endl;
  logger << LOG_DEBUG_VERBOSE << "xbsqrts_all[0][0] = " << xbsqrts_all[0][0] << endl;
  logger << LOG_DEBUG_VERBOSE << "x_i_ab = " << x_i_ab << endl;
  logger << LOG_DEBUG_VERBOSE << "one_minus_x_i_ab = " << one_minus_x_i_ab << endl;

  if (x_i_ab != x_i_ab){logger << LOG_DEBUG_VERBOSE << "x_i_ab = " << x_i_ab << endl;}
  xbp_all[i_a][bp_a] = x_i_ab * xbp_all[0][bp_a];
  xbp_all[i_a][bp_b] = xbp_all[0][bp_b];
  xbp_all[i_a][0] = xbp_all[i_a][bp_a] + xbp_all[i_a][bp_b];
  xbs_all[i_a][0] = x_i_ab * xbs_all[0][0];
  xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);
  xbp_all[i_a][xb_out_dipoles] = xbp_all[i_a][0];
  xbs_all[i_a][xb_out_dipoles] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xb_out_dipoles] = xbsqrts_all[i_a][0];
  fourvector temp_p_aib = xbp_all[i_a][bp_a] + xbp_all[i_a][bp_b];
  fourvector temp_p_a_boost = xbp_all[i_a][bp_a].boost(temp_p_aib);
  double sqrt_p2_aib = temp_p_aib.m();
  double v_i = h_propto_pot_mod(r[no_random_dipole[1]], 0., one_minus_x_i_ab, exp_ai_b_v, map_technical_x);
logger << LOG_DEBUG_VERBOSE << "v_i = " << v_i << endl;
  if (v_i != v_i){logger << LOG_DEBUG_VERBOSE << "v_i = " << v_i << endl;}
  double one_minus_v_i = 1. - v_i;
  double phi_k_perp = c_phi(no_random_dipole[2]);
  if (phi_k_perp != phi_k_perp){logger << LOG_DEBUG_VERBOSE << "phi_k_perp = " << phi_k_perp << endl;}
  double f2 = v_i * (one_minus_x_i_ab - v_i) / x_i_ab;
  double f = sqrt(f2);
  fourvector khat_perp(0., f * sqrt_p2_aib * cos(phi_k_perp), f * sqrt_p2_aib * sin(phi_k_perp), 0.);
  fourvector k_perp = (khat_perp.rotateback(temp_p_a_boost)).boost(temp_p_aib.Pinv());
  fourvector K = (x_i_ab + v_i) / x_i_ab * xbp_all[i_a][bp_a] + one_minus_v_i * xbp_all[i_a][bp_b] - k_perp;
  if (K != K){logger << LOG_DEBUG_VERBOSE << "K = " << K << endl;}
  fourvector Kt = xbp_all[i_a][bp_a] + xbp_all[i_a][bp_b];
  if (Kt != Kt){logger << LOG_DEBUG_VERBOSE << "Kt = " << Kt << endl;}
  fourvector KKt = K + Kt;
  xbp_all[0][bp_i] = v_i * xbp_all[i_a][bp_b] + (one_minus_x_i_ab - v_i) / x_i_ab * xbp_all[i_a][bp_a] + k_perp;
  xbp_all[0][bp_a] = xbp_all[i_a][bp_a] / x_i_ab;
  xbp_all[0][bp_b] = xbp_all[i_a][bp_b];
  /*
  logger << LOG_DEBUG_VERBOSE << "xbp_all[" << i_a << "][" << bp_a << "] = " << xbp_all[i_a][bp_a] << endl;
  logger << LOG_DEBUG_VERBOSE << "xbp_all[" << i_a << "][" << bp_b << "] = " << xbp_all[i_a][bp_b] << endl;
  for (int xbi = 4; xbi < xbp_all[i_a].size(); xbi = xbi * 2){
    logger << LOG_DEBUG_VERBOSE << "xbp_all[" << i_a << "][" << xbi << "] = " << xbp_all[i_a][xbi] << endl;
  }
  logger << LOG_DEBUG_VERBOSE << "xb_all_out = " << xb_all_out << endl;
  logger << LOG_DEBUG_VERBOSE << "xbp_all[" << i_a << "][" << xb_all_out << "] = " << xbp_all[i_a][xb_all_out] << endl;
  logger << LOG_DEBUG_VERBOSE << "xbs_all[" << i_a << "][" << xb_all_out << "] = " << xbs_all[i_a][xb_all_out] << endl;
  */
  generic->ac_psp_dipole(i_a, MC_phasespace.channel - dipole[i_a - 1].sum_channel(), *this);
  /*
  logger << LOG_DEBUG_VERBOSE << "channel - dipole[" << i_a - 1 << "].sum_channel() = " << channel - dipole[i_a - 1].sum_channel() << endl;
  logger << LOG_DEBUG_VERBOSE << "xbp_all[" << i_a << "][" << bp_a << "] = " << xbp_all[i_a][bp_a] << endl;
  logger << LOG_DEBUG_VERBOSE << "xbp_all[" << i_a << "][" << bp_b << "] = " << xbp_all[i_a][bp_b] << endl;
  for (int xbi = 4; xbi < xbp_all[i_a].size(); xbi = xbi * 2){
    logger << LOG_DEBUG_VERBOSE << "xbp_all[" << i_a << "][" << xbi << "] = " << xbp_all[i_a][xbi] << endl;
  }
  */
  double prod_em_sp = xbp_all[i_a][bp_a] * xbp_all[i_a][bp_b];
  for (int xbi = 4; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else {
      int temp_xbi = 0;
      if      (xbi < bp_i){temp_xbi = xbi;}
      else if (xbi > bp_i){temp_xbi = xbi / 2;}
      double prod_em_xbi = xbp_all[i_a][temp_xbi] * xbp_all[i_a][bp_a];
      double prod_sp_xbi = xbp_all[i_a][temp_xbi] * xbp_all[i_a][bp_b];
      double prod_perp_xbi = xbp_all[i_a][temp_xbi] * k_perp;

      xbp_all[0][xbi] = xbp_all[i_a][temp_xbi]
        - ((2. * x_i_ab + v_i) * prod_em_xbi + (2. - v_i) * x_i_ab * prod_sp_xbi - x_i_ab * prod_perp_xbi) / ((4. * x_i_ab + v_i - x_i_ab * v_i) * prod_em_sp) * KKt
        + (prod_em_xbi + prod_sp_xbi) / prod_em_sp * K;

      if (xbp_all[0][xbi]  != xbp_all[0][xbi]){
	logger << LOG_DEBUG_VERBOSE << "i_acc = " << i_acc << "   xbi = " << xbi << "   xbp_all[" << i_a << "][" << temp_xbi << "] = " << xbp_all[i_a][temp_xbi] << endl;
	logger << LOG_DEBUG_VERBOSE << "K = " << K << endl;
	logger << LOG_DEBUG_VERBOSE << "Kt = " << Kt << endl;
	logger << LOG_DEBUG_VERBOSE << "x_i_ab = " << x_i_ab << endl;
	logger << LOG_DEBUG_VERBOSE << "v_i = " << v_i << endl;
	logger << LOG_DEBUG_VERBOSE << "phi_k_perp = " << phi_k_perp << endl;
	logger << LOG_DEBUG_VERBOSE << "prod_em_sp = " << prod_em_sp << endl;
	logger << LOG_DEBUG_VERBOSE << "prod_em_xbi = " << prod_em_xbi << endl;
	logger << LOG_DEBUG_VERBOSE << "prod_sp_xbi = " << prod_sp_xbi << endl;
	logger << LOG_DEBUG_VERBOSE << "xbp_all[" << i_a << "][" << bp_a << "] = " << xbp_all[i_a][bp_a] << endl;
	logger << LOG_DEBUG_VERBOSE << "xbp_all[" << i_a << "][" << bp_b << "] = " << xbp_all[i_a][bp_b] << endl;
	logger << LOG_DEBUG_VERBOSE << "xbp_all[" << i_a << "][" << temp_xbi << "] = " << xbp_all[i_a][temp_xbi] << endl;
	//	exit(1);
	/*
	logger << LOG_DEBUG_VERBOSE << "1st line: " << xbp_all[i_a][temp_xbi] << endl;
	logger << LOG_DEBUG_VERBOSE << "2nd line: " << - ((2. * x_i_ab + v_i) * prod_em_xbi + (2. - v_i) * x_i_ab * prod_sp_xbi - x_i_ab * prod_perp_xbi) / ((4. * x_i_ab + v_i - x_i_ab * v_i) * prod_em_sp) * KKt << endl;
	logger << LOG_DEBUG_VERBOSE << "3rd line: " << + (prod_em_xbi + prod_sp_xbi) / prod_em_sp * K << endl;
	*/
     }
    }
  }
  for (int xbi = 4; xbi < xbp_all[0].size(); xbi = xbi * 2){
    logger << LOG_DEBUG_VERBOSE << "xbp_all[" << 0 << "][" << xbi << "] = " << xbp_all[0][xbi] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::ag_psp_RA_group_ai_b(double & sinx_i_ab_min, vector<dipole_set> & dipole, int i_a){
  static Logger logger("phasespace_set::ag_psp_RA_group_ai_b");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int bp_a = dipole[i_a].binary_R_emitter_1();
  int bp_i = dipole[i_a].binary_R_emitter_2();
  int bp_b = dipole[i_a].binary_R_spectator();
  //  int bpt_ai = dipole[i_a].binary_A_emitter(); // always bpt_ai == bp_a
  //  int bpt_b = dipole[i_a].binary_A_spectator(); // always bpt_b == bp_b



  /*
  double one_minus_x_i_ab_min = 1. - sinx_i_ab_min / xbs_all[0][0];
  double one_minus_x_i_ab = 2 * xbp_all[0][bp_i].x0() / xbsqrts_all[0][0];
  double x_i_ab = 1. - one_minus_x_i_ab;
  double g_x_i_ab = g_propto_pot_mod(one_minus_x_i_ab, 0., one_minus_x_i_ab_min, exp_ai_b_x, map_technical_x)
  */
  double one_minus_x_i_ab_min = 1. - sinx_i_ab_min / xbs_all[0][0];
  double one_minus_x_i_ab = 2 * xbp_all[0][bp_i].x0() / xbsqrts_all[0][0];
  double x_i_ab = 1. - one_minus_x_i_ab;
  double g_x_i_ab;
  double x_i_ab_s_part = x_i_ab * xbs_all[0][0];






  
  logger << LOG_DEBUG_VERBOSE << "x_i_ab_s_part = " << x_i_ab_s_part << "   sqrt -> " << sqrt(x_i_ab_s_part) << endl;

  for (int i_c = 0; i_c < MC_x_dipole[i_a].g_channel.size(); i_c++){
    if (i_c == 0){
      MC_x_dipole[i_a].g_channel[i_c] = g_propto_pot_mod(one_minus_x_i_ab, 0., one_minus_x_i_ab_min, exp_ai_b_x, map_technical_x);
    }
    else {
      if (MC_x_dipole_mapping[i_a][i_c] < 0){MC_x_dipole[i_a].g_channel[i_c] = xbs_all[0][0] * g_propagator_Breit_Wigner(x_i_ab_s_part, -MC_x_dipole_mapping[i_a][i_c], sinx_i_ab_min, xbs_all[0][0]);}
      logger << LOG_DEBUG_VERBOSE << "sinx_i_ab_min = " << setw(23) << setprecision(15) << sinx_i_ab_min << endl;
      logger << LOG_DEBUG_VERBOSE << "x_i_ab_s_part = " << setw(23) << setprecision(15) << x_i_ab_s_part << endl;
      logger << LOG_DEBUG_VERBOSE << "xbs_all[0][0] = " << setw(23) << setprecision(15) << xbs_all[0][0] << endl;
      logger << LOG_DEBUG_VERBOSE << "MC_x_dipole_mapping.size() = " << MC_x_dipole_mapping.size() << endl;
      logger << LOG_DEBUG_VERBOSE << "MC_x_dipole_mapping[" << i_a << "].size() = " << MC_x_dipole_mapping[i_a].size() << endl;
      logger << LOG_DEBUG_VERBOSE << "MC_x_dipole.size() = " << MC_x_dipole.size() << endl;
      //      logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[" << i_a << "].channel = " << MC_x_dipole[i_a].channel << endl;
      //      logger << LOG_DEBUG_VERBOSE << "-MC_x_dipole_mapping[i_a][i_c] = " << setw(23) << setprecision(15) << -MC_x_dipole_mapping[i_a][MC_x_dipole[i_a].channel] << endl;
      //      logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[i_a].g_channel[i_c] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].g_channel[i_c] << endl;
      ///      else if (tau_MC_map[i_c] > 0){MC_x_dipole[i_a].g_channel[i_c] = s_had * g_propagator_vanishing_width(x_i_ab_s_part, M2[tau_MC_map[i_c]], tau_0_s_had, s_had, nuxs);}
    }
  }
  g_x_i_ab = 0.;
  for (int i = 0; i < MC_x_dipole[i_a].g_channel.size(); i++){
    logger << LOG_DEBUG_VERBOSE << "MC_x_dipole[" << i_a << "].alpha[" << i << "] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].alpha[i] << "   " << "MC_x_dipole[" << i_a << "].g_channel[" << i << "] = " << setw(23) << setprecision(15) << MC_x_dipole[i_a].g_channel[i] << endl;
  }
  for (int i = 0; i < MC_x_dipole[i_a].g_channel.size(); i++){g_x_i_ab += MC_x_dipole[i_a].alpha[i] * MC_x_dipole[i_a].g_channel[i];}
  logger << LOG_DEBUG_VERBOSE << "g_x_i_ab = " << g_x_i_ab << endl;



  double v_i = (xbp_all[0][bp_i] * xbp_all[0][bp_a]) / (xbp_all[0][bp_b] * xbp_all[0][bp_a]);
  double g_alpha = g_x_i_ab * g_propto_pot_mod(v_i, 0., one_minus_x_i_ab, exp_ai_b_v, map_technical_x) / (pi * (xbp_all[0][bp_b] * xbp_all[0][bp_a]));
  if (switch_console_output_phasespace_issue){
    if (g_alpha < 0.){logger << LOG_INFO << "i_a = " << i_a << "   g_alpha = " << g_alpha << endl;}
  }
  //  double g_alpha = g_propto_pot_mod(one_minus_x_i_ab, 0., one_minus_x_i_ab_min, exp_ai_b_x, map_technical_x) * g_propto_pot_mod(v_i, 0., one_minus_x_i_ab, exp_ai_b_v, map_technical_x) / (pi * (xbp_all[0][bp_b] * xbp_all[0][bp_a]));
  if (xbs_all[i_a][0] == 0.){xbs_all[i_a][0] = x_i_ab * xbs_all[0][0];}
  if (xbsqrts_all[i_a][0] == 0.){xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);}
  generic->ag_psp_dipole(i_a, dipole[i_a - 1].sum_channel(), *this);

  for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){MC_phasespace.g_channel[i] = MC_phasespace.g_channel[i] * g_alpha;}
  // invert mappings for VAMP implementation
  if (switch_IS_mode_phasespace == 2 || switch_IS_mode_phasespace == 4){
    //  static vector<double> inv_r(3);
    vector<double> inv_r(3);
    inv_r[0] = inv_propto_pot(one_minus_x_i_ab, 0., one_minus_x_i_ab_min, exp_ai_b_x, map_technical_x);
    inv_r[1] = inv_propto_pot(v_i, 0., one_minus_x_i_ab, exp_ai_b_v, map_technical_x);
    fourvector k_perp = xbp_all[0][bp_i] - v_i * xbp_all[i_a][bp_b] + (one_minus_x_i_ab - v_i) / x_i_ab * xbp_all[i_a][bp_a];
    fourvector temp_p_aib = xbp_all[i_a][bp_a] + xbp_all[i_a][bp_b];
    fourvector temp_p_a_boost = xbp_all[i_a][bp_a].boost(temp_p_aib);
    k_perp = (k_perp.boost(temp_p_aib)).rotateback_inverse(temp_p_a_boost);
    double phi = atan2(k_perp.x2(),k_perp.x1());
    if (phi < 0){phi += 2 * M_PI;}
    inv_r[2] = inv_phi(phi);
    if (!(inv_r[2] >= 0) || inv_r[2] > 1){
      if (switch_console_output_phasespace_issue){
	logger << LOG_INFO << "inv_r[2] not in [0; 1] !   r = " << inv_r[2] <<"; k_perp.x,y,z=" << k_perp.x1() << ", " << k_perp.x2() << ", " << k_perp.x3() << endl;
      }
      if (!(inv_r[2] >= 0)){inv_r[2] = 0.;}
      if (inv_r[2] > 1){inv_r[2] = 1.;}
    }
    if (switch_IS_mode_phasespace == 2){
      for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){
        double g_IS_temp;
        double g_IS = 1.;
        for (int j = 0; j < 3; j++){
          phasespace_randoms[(i + 1) * 8 - 3 + j]->get_g_IS(inv_r[j], g_IS_temp);
          g_IS *= g_IS_temp;
        }
        MC_phasespace.g_channel[i] *= g_IS;
      }
    }
    else if (switch_IS_mode_phasespace == 4){
      double g_IS_temp;
      double g_IS = 1.;
      for (int j = 2; j < 3; j++){ // manually switched of x and z IS
        phasespace_randoms[container_IS_startvalue[i_a][4 + j]]->get_g_IS(inv_r[j], g_IS_temp);
        g_IS *= g_IS_temp;
      }
      for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){MC_phasespace.g_channel[i] *= g_IS;}
     }
    //  logger << LOG_INFO << "reconstr " << i_a << ": " << inv_r[0] << ", " << inv_r[1] << ", " << inv_r[2] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



