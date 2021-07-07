#include "../include/classes.cxx"
//#include "../include/definitions.phasespace.set.cxx"

#define p_i xbp_all[0][bp_i]
#define p_j xbp_all[0][bp_j]
#define p_k xbp_all[0][bp_k]
#define m_i xbsqrts_all[0][bp_i]
#define m_j xbsqrts_all[0][bp_j]
#define m_k xbsqrts_all[0][bp_k]
#define m_ij xbsqrts_all[i_a][bpt_ij]
#define m2_i xbs_all[0][bp_i]
#define m2_j xbs_all[0][bp_j]
#define m2_k xbs_all[0][bp_k]
#define m2_ij xbs_all[i_a][bpt_ij]
void phasespace_set::ac_psp_RA_group_ij_k_massive(vector<dipole_set> & dipole, int i_a){
  static Logger logger("phasespace_set::ac_psp_RA_group_ij_k_massive");
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
  fourvector Q = xbp_all[i_a][bpt_ij] + xbp_all[i_a][bpt_k];
  double Q2 = Q.m2();
  double sqrtQ2 = sqrt(Q2);
  double mu_i = m_i / sqrtQ2;
  double mu_j = m_j / sqrtQ2;
  double mu_k = m_k / sqrtQ2;
  //  double mu_ij = m_ij / sqrtQ2;
  double mu2_i = m2_i / Q2;
  double mu2_j = m2_j / Q2;
  double mu2_k = m2_k / Q2;
  //  double mu2_ij = m2_ij / Q2;
  double one_minus_all_mu = 1. - mu2_i - mu2_j - mu2_k;
  double y_minus = (2 * mu_i * mu_j) / one_minus_all_mu;
  double y_plus = 1. - (2 * mu_k * (1. - mu_k)) / one_minus_all_mu;
  double y_ij_k = h_propto_pot_mod(r[no_random_dipole[0]], y_minus, y_plus, exp_ij_k_y, map_technical_x);
  double one_minus_y_ij_k = 1. - y_ij_k;
  double z_prefactor = (2 * mu2_i + one_minus_all_mu * y_ij_k) / (2 * (mu2_i + mu2_j + one_minus_all_mu * y_ij_k));
  double v_ij_k = sqrt(pow(2 * mu2_k + one_minus_all_mu * one_minus_y_ij_k, 2) - 4 * mu2_k) / (one_minus_all_mu * one_minus_y_ij_k);
  double v_ij_i = sqrt(pow(one_minus_all_mu * y_ij_k, 2) - 4 * mu2_i * mu2_j) / (one_minus_all_mu * y_ij_k + 2 * mu2_i);
  double z_minus = z_prefactor * (1 - v_ij_i * v_ij_k);
  double z_plus = z_prefactor * (1 + v_ij_i * v_ij_k);
  double one_minus_z_i = h_propto_pot_mod(r[no_random_dipole[1]], 1. - z_plus, 1. - z_minus, exp_ij_k_z, map_technical_x);
  double z_i = 1. - one_minus_z_i;
  double phi_k_perp = c_phi(no_random_dipole[2]);
  fourvector temp_p_k_boost = xbp_all[i_a][bpt_k].boost(Q);
  double s_ij = y_ij_k * (Q2 - m2_k - m2_i - m2_j) + m2_i + m2_j;
  double s2_ij = pow(s_ij, 2);
  double lambda_Q2_sij_mk2 = lambda(Q2, s_ij, m2_k);
  double sqrtlambda_Q2_sij_mk2 = sqrt(lambda_Q2_sij_mk2);
  fourvector Q_phat_k((Q2 - s_ij + m2_k) / (2 * sqrtQ2), 0., 0., sqrtlambda_Q2_sij_mk2 / (2 * sqrtQ2));
  fourvector Q_phat_ij((Q2 + s_ij - m2_k) / (2 * sqrtQ2), 0., 0., -sqrtlambda_Q2_sij_mk2 / (2 * sqrtQ2));
  double help_p_i_0 = (z_i * (Q2 - m2_k) + one_minus_z_i * s_ij + m2_i - m2_j) / (2 * sqrtQ2);
  double help_p_i_3 = (s_ij * (Q2 + m2_k) - z_i * pow(Q2 - m2_k, 2) - one_minus_z_i * s2_ij + (Q2 - s_ij + m2_k) * (m2_i - m2_j)) / (2 * sqrtQ2 * sqrtlambda_Q2_sij_mk2);
  double help_p_i_T;
  if (m2_i == 0. && m2_j == 0.){
    help_p_i_T = sqrt((s_ij * z_i * one_minus_z_i * pow(Q2 - s_ij - m2_k, 2) - s_ij * m2_k * (s_ij + 2 * (one_minus_z_i - z_i) * (m2_i - m2_j)) - m2_k * pow(m2_i - m2_j, 2)) / lambda_Q2_sij_mk2 - z_i * m2_j - one_minus_z_i * m2_i);
  }
  else {
    double delta_1 = y_ij_k * (Q2 - m2_k - m2_i - m2_j) / (m2_i + m2_j);
    help_p_i_T = sqrt((s_ij * z_i * one_minus_z_i * pow(Q2 - s_ij - m2_k, 2) + m2_k * (4 * m2_i * m2_j - (m2_i + m2_j) * (4 * (1 + delta_1) * (one_minus_z_i * m2_i + z_i * m2_j) + pow(delta_1, 2) * (m2_i + m2_j)))) / lambda_Q2_sij_mk2 - z_i * m2_j - one_minus_z_i * m2_i);
  }
  fourvector Q_phat_i(help_p_i_0, help_p_i_T * cos(phi_k_perp), help_p_i_T * sin(phi_k_perp), help_p_i_3);
  fourvector Q_phat_j = Q_phat_ij - Q_phat_i;

  for (int xbi = 4; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){xbp_all[0][bp_i] = (Q_phat_i.rotateback(temp_p_k_boost)).boost(Q.Pinv());}
    else if (xbi == bp_j){xbp_all[0][bp_j] = (Q_phat_j.rotateback(temp_p_k_boost)).boost(Q.Pinv());}
    else if (xbi == bp_k){xbp_all[0][bp_k] = (Q_phat_k.rotateback(temp_p_k_boost)).boost(Q.Pinv());}
    else if (xbi <  bp_j){xbp_all[0][xbi] = xbp_all[i_a][xbi];}
    else if (xbi >  bp_j){xbp_all[0][xbi] = xbp_all[i_a][xbi / 2];}
    else {logger << LOG_INFO << "Should not happen!" << endl;}
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::ag_psp_RA_group_ij_k_massive(vector<dipole_set> & dipole, int i_a){
  static Logger logger("phasespace_set::ag_psp_RA_group_ij_k_massive");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int bp_i = dipole[i_a].binary_R_emitter_1();
  int bp_j = dipole[i_a].binary_R_emitter_2();
  int bp_k = dipole[i_a].binary_R_spectator();
  int bpt_ij = dipole[i_a].binary_A_emitter();
  int bpt_k = dipole[i_a].binary_A_spectator();
  if (xbs_all[0][bp_i + bp_j] == 0.){xbs_all[0][bp_i + bp_j] = (xbp_all[0][bp_i] + xbp_all[0][bp_j]).m2();}
  double pi_pj = p_i * p_j;
  double pj_pk = p_j * p_k;
  double pi_pk = p_i * p_k;
  double y_ij_k = pi_pj / (pi_pj + pj_pk + pi_pk);
  double one_minus_y_ij_k = 1. - y_ij_k;
  double z_i = pi_pk / (pi_pk + pj_pk);
  double one_minus_z_i = pj_pk / (pi_pk + pj_pk);
  fourvector Q = xbp_all[0][bp_i] + xbp_all[0][bp_j] + xbp_all[0][bp_k];
  double Q2 = Q.m2();
  double sqrtQ2 = sqrt(Q2);
  double mu_i = m_i / sqrtQ2;
  double mu_j = m_j / sqrtQ2;
  double mu_k = m_k / sqrtQ2;
  //  double mu_ij = m_ij / sqrtQ2;
  double mu2_i = m2_i / Q2;
  double mu2_j = m2_j / Q2;
  double mu2_k = m2_k / Q2;
  double mu2_ij = m2_ij / Q2;
  double one_minus_all_mu = 1. - mu2_i - mu2_j - mu2_k;
  double y_minus = (2 * mu_i * mu_j) / one_minus_all_mu;
  double y_plus = 1. - (2 * mu_k * (1. - mu_k)) / one_minus_all_mu;
  double z_prefactor = (2 * mu2_i + one_minus_all_mu * y_ij_k) / (2 * (mu2_i + mu2_j + one_minus_all_mu * y_ij_k));
  double v_ij_k = sqrt(pow(2 * mu2_k + one_minus_all_mu * one_minus_y_ij_k, 2) - 4 * mu2_k) / (one_minus_all_mu * one_minus_y_ij_k);
  double v_ij_i = sqrt(pow(one_minus_all_mu * y_ij_k, 2) - 4 * mu2_i * mu2_j) / (one_minus_all_mu * y_ij_k + 2 * mu2_i);
  double z_minus = z_prefactor * (1 - v_ij_i * v_ij_k);
  double z_plus = z_prefactor * (1 + v_ij_i * v_ij_k);
  double g_alpha = g_propto_pot_mod(y_ij_k, y_minus, y_plus, exp_ij_k_y, map_technical_x) * g_propto_pot_mod(one_minus_z_i, 1. - z_plus, 1. - z_minus, exp_ij_k_z, map_technical_x) / (.5 * pi * Q2 * pow(one_minus_all_mu, 2) / sqrt(lambda(1., mu2_ij, mu2_k)) * (1. - y_ij_k));
  if (xbs_all[i_a][0] == 0.){xbs_all[i_a][0] = xbs_all[0][0];}
  if (xbsqrts_all[i_a][0] == 0.){xbsqrts_all[i_a][0] = xbsqrts_all[0][0];}
  generic->ag_psp_dipole(i_a, dipole[i_a - 1].sum_channel(), *this);
  for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){MC_phasespace.g_channel[i] = MC_phasespace.g_channel[i] * g_alpha;}
  // to be adapted to massive case !!!
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
#undef p_i 
#undef p_j 
#undef p_k 
#undef m_i
#undef m_j
#undef m_k
#undef m_ij
#undef m2_i
#undef m2_j
#undef m2_k
#undef m2_ij



#define p_i xbp_all[0][bp_i]
#define p_j xbp_all[0][bp_j]
#define p_a xbp_all[0][bp_a]
#define pt_ij xbp_all[i_a][bpt_ij]
#define m_i xbsqrts_all[0][bp_i]
#define m_j xbsqrts_all[0][bp_j]
#define m_ij xbsqrts_all[i_a][bpt_ij]
#define m2_i xbs_all[0][bp_i]
#define m2_j xbs_all[0][bp_j]
#define m2_ij xbs_all[i_a][bpt_ij]
void phasespace_set::ac_psp_RA_group_ij_a_massive(vector<dipole_set> & dipole, double & sinx_ij_a_min, int i_a){
  static Logger logger("phasespace_set::ac_psp_RA_group_ij_a_massive");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double g_IS_temp;
  if (switch_IS_mode_phasespace == 3 || switch_IS_mode_phasespace == 4){
    for (int j = 2; j < 3; j++){ // manually switched of x and z IS
      if (container_IS_switch[container_IS_startvalue[i_a][4 + j]] == 1){phasespace_randoms[container_IS_startvalue[i_a][4 + j]]->get_random(r[no_random_dipole[j]], g_IS_temp);}
      MC_g_IS_global *= g_IS_temp;
    }
  }
  int bp_a = dipole[i_a].binary_R_spectator();
  int bp_i = dipole[i_a].binary_R_emitter_1();
  int bp_j = dipole[i_a].binary_R_emitter_2();
  int bpt_ij = dipole[i_a].binary_A_emitter();
  //  int bpt_a = dipole[i_a].binary_A_spectator(); // always bp_a == bpt_a
  int bp_b = bp_a % 2 + 1;
  double xplus = 1.;
  //  double xplus = 1. + mu2_ij - pow(mu_i + mu_j, 2);   cannot be calculated here !!! check realization for m_i = m_j = m_Q, m_ij = 0 (g -> QQ~ splitting) 
  double one_minus_x_ij_a_min = 1. - sinx_ij_a_min / xbs_all[0][0];
  double one_minus_x_ij_a = h_propto_pot_mod(r[no_random_dipole[0]], 1. - xplus, one_minus_x_ij_a_min, exp_ij_a_x, map_technical_x);
  double x_ij_a = 1. - one_minus_x_ij_a;
  xbp_all[i_a][bp_a] = x_ij_a * xbp_all[0][bp_a];
  xbp_all[i_a][bp_b] = xbp_all[0][bp_b];
  xbp_all[i_a][0] = xbp_all[i_a][bp_a] + xbp_all[i_a][bp_b];
  xbs_all[i_a][0] = x_ij_a * xbs_all[0][0];
  xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);
  xbp_all[i_a][xb_out_dipoles] = xbp_all[i_a][0];
  xbs_all[i_a][xb_out_dipoles] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xb_out_dipoles] = xbsqrts_all[i_a][0];
  generic->ac_psp_dipole(i_a, MC_phasespace.channel - dipole[i_a - 1].sum_channel(), *this);
  double two_ptij_pa = 2 * pt_ij * p_a;
  //  double sqrt_two_ptij_pa = sqrt(two_ptij_pa);
  //  double mu_i = m_i / sqrt_two_ptij_pa;
  //  double mu_j = m_j / sqrt_two_ptij_pa;
  //  double mu_ij = m_ij / sqrt_two_ptij_pa;
  double mu2_i = m2_i / two_ptij_pa;
  double mu2_j = m2_j / two_ptij_pa;
  double mu2_ij = m2_ij / two_ptij_pa;
  double help_z_bound_symm_1 = one_minus_x_ij_a + mu2_ij - mu2_i - mu2_j;
  double help_z_bound_symm_2 = sqrt(pow(help_z_bound_symm_1, 2) - 4 * mu2_i * mu2_j);
  double help_z_bound_ij = one_minus_x_ij_a + mu2_ij + mu2_i - mu2_j;
  double zminus = (help_z_bound_ij - help_z_bound_symm_2) / (2 * (one_minus_x_ij_a + mu2_ij));
  double zplus = (help_z_bound_ij + help_z_bound_symm_2) / (2 * (one_minus_x_ij_a + mu2_ij));
  double one_minus_z_i = h_propto_pot_mod(r[no_random_dipole[1]], 1. - zplus, 1. - zminus, exp_ij_a_z, map_technical_x);//exp_ij_a());
  double z_i = 1. - one_minus_z_i;
  double phi_a_perp = c_phi(no_random_dipole[2]);
  fourvector Q = xbp_all[i_a][bpt_ij] + xbp_all[i_a][bp_a];
  double Q2 = Q.m2();
  double sqrtQ2 = sqrt(Q2);
  fourvector temp_p_a_boost = xbp_all[i_a][bp_a].boost(Q);
  // not used !   fourvector Q_phattilde_a((Q2 - m2_ij) / (2 * sqrtQ2), 0., 0., (Q2 - m2_ij) / (2 * sqrtQ2));
  // not used !   fourvector Q_phattilde_ij((Q2 + m2_ij) / (2 * sqrtQ2), 0., 0., -(Q2 - m2_ij) / (2 * sqrtQ2));
  // not used !   fourvector Q_phat_a = (1. / x_ij_a) * Q_phattilde_a;
  double help_p_i = one_minus_z_i * one_minus_x_ij_a * (Q2 - m2_ij) / (2 * x_ij_a * sqrtQ2) + (one_minus_z_i * m2_ij + m2_i - m2_j) / (2 * sqrtQ2);
  double help_p_i_T = sqrt((z_i * one_minus_z_i * one_minus_x_ij_a / x_ij_a) * (Q2 - m2_ij) + z_i * one_minus_z_i * m2_ij - one_minus_z_i * m2_i - z_i * m2_j);
  fourvector Q_phat_i(help_p_i + .5 * z_i * sqrtQ2, help_p_i_T * cos(phi_a_perp), help_p_i_T * sin(phi_a_perp), help_p_i - .5 * z_i * sqrtQ2);
  fourvector Q_phat_ij((Q2 - (one_minus_x_ij_a - x_ij_a) * m2_ij) / (2 * x_ij_a * sqrtQ2), 0., 0., (one_minus_x_ij_a - x_ij_a) * (Q2 - m2_ij) / (2 * x_ij_a * sqrtQ2));
  fourvector Q_phat_j = Q_phat_ij - Q_phat_i;
  for (int xbi = 4; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){xbp_all[0][bp_i] = (Q_phat_i.rotateback(temp_p_a_boost)).boost(Q.Pinv());}
    else if (xbi == bp_j){xbp_all[0][bp_j] = (Q_phat_j.rotateback(temp_p_a_boost)).boost(Q.Pinv());}
    else if (xbi < bp_j){xbp_all[0][xbi] = xbp_all[i_a][xbi];}
    else if (xbi > bp_j){xbp_all[0][xbi] = xbp_all[i_a][xbi / 2];}
    else {logger << LOG_INFO << "Should not happen!" << endl;}
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::ag_psp_RA_group_ij_a_massive(double & sinx_ij_a_min, vector<dipole_set> & dipole, int i_a){
  static Logger logger("phasespace_set::ag_psp_RA_group_ij_a_massive");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int bp_i = dipole[i_a].binary_R_emitter_1();
  int bp_j = dipole[i_a].binary_R_emitter_2();
  int bp_a = dipole[i_a].binary_R_spectator();
  int bpt_ij = dipole[i_a].binary_A_emitter();
  if (xbs_all[0][bp_i + bp_j] == 0.){xbs_all[0][bp_i + bp_j] = (xbp_all[0][bp_i] + xbp_all[0][bp_j]).m2();}
  double pi_pj = p_i * p_j;
  double pj_pa = p_j * p_a;
  double pi_pa = p_i * p_a;
  double two_ptij_pa = 2 * pt_ij * p_a;
  //  double sqrt_two_ptij_pa = sqrt(two_ptij_pa);
  //  double mu_i = m_i / sqrt_two_ptij_pa;
  //  double mu_j = m_j / sqrt_two_ptij_pa;
  //  double mu_ij = m_ij / sqrt_two_ptij_pa;
  double mu2_i = m2_i / two_ptij_pa;
  double mu2_j = m2_j / two_ptij_pa;
  double mu2_ij = m2_ij / two_ptij_pa;
  double xplus = 1.;
  //  double xplus = 1. + mu2_ij - pow(mu_i + mu_j, 2);   cannot be calculated here !!! check realization for m_i = m_j = m_Q, m_ij = 0 (g -> QQ~ splitting) 
  double one_minus_x_ij_a_min = 1. - sinx_ij_a_min / xbs_all[0][0];
  double one_minus_x_ij_a = (pi_pj - .5 * (m2_ij - m2_i - m2_j)) / (pi_pa + pj_pa);
  double x_ij_a = 1. - one_minus_x_ij_a;
  double z_i = pi_pa / (pi_pa + pj_pa);
  double z_j = pj_pa / (pi_pa + pj_pa);
  double one_minus_z_i = z_j;
  double help_z_bound_symm_1 = one_minus_x_ij_a + mu2_ij - mu2_i - mu2_j;
  double help_z_bound_symm_2 = sqrt(pow(help_z_bound_symm_1, 2) - 4 * mu2_i * mu2_j);
  double help_z_bound_ij = one_minus_x_ij_a + mu2_ij + mu2_i - mu2_j;
  double zminus = (help_z_bound_ij - help_z_bound_symm_2) / (2 * (one_minus_x_ij_a + mu2_ij));
  double zplus = (help_z_bound_ij + help_z_bound_symm_2) / (2 * (one_minus_x_ij_a + mu2_ij));
  double g_alpha = g_propto_pot_mod(one_minus_x_ij_a, 1. - xplus, one_minus_x_ij_a_min, exp_ij_a_x, map_technical_x) * g_propto_pot_mod(z_j, 1. - zplus, 1. - zminus, exp_ij_a_z, map_technical_x) / (pi * (xbp_all[i_a][bpt_ij] * xbp_all[0][bp_a]));
  if (xbs_all[i_a][0] == 0.){xbs_all[i_a][0] = x_ij_a * xbs_all[0][0];}
  if (xbsqrts_all[i_a][0] == 0.){xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);}
  generic->ag_psp_dipole(i_a, dipole[i_a - 1].sum_channel(), *this);
  for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){MC_phasespace.g_channel[i] = MC_phasespace.g_channel[i] * g_alpha;}
  // to be adapted to massive case !!!
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
#undef m_i
#undef m_j
#undef m_ij
#undef m2_i
#undef m2_j
#undef m2_ij



#define p_i xbp_all[0][bp_i]
#define p_k xbp_all[0][bp_k]
#define p_a xbp_all[0][bp_a]
#define pt_k xbp_all[i_a][bpt_k]
#define m_k xbsqrts_all[0][bp_k]
#define m2_k xbs_all[0][bp_k]
void phasespace_set::ac_psp_RA_group_ai_k_massive(vector<dipole_set> & dipole, double & sinx_ik_a_min, int i_a){
  static Logger logger("phasespace_set::ac_psp_RA_group_ai_k_massive");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double g_IS_temp;
  logger << LOG_DEBUG_VERBOSE << "1" << endl;

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
  double one_minus_x_ik_a_min = 1. - sinx_ik_a_min / xbs_all[0][0];
  double one_minus_x_ik_a = h_propto_pot_mod(r[no_random_dipole[0]], 0., one_minus_x_ik_a_min, exp_ai_k_x, map_technical_x);
  double x_ik_a = 1. - one_minus_x_ik_a;
  xbp_all[i_a][bp_a] = x_ik_a * xbp_all[0][bp_a];
  xbp_all[i_a][bp_b] = xbp_all[0][bp_b];
  xbp_all[i_a][0] = xbp_all[i_a][bp_a] + xbp_all[i_a][bp_b];
  xbs_all[i_a][0] = x_ik_a * xbs_all[0][0];
  xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);
  xbp_all[i_a][xb_out_dipoles] = xbp_all[i_a][0];
  xbs_all[i_a][xb_out_dipoles] = xbs_all[i_a][0];
  xbsqrts_all[i_a][xb_out_dipoles] = xbsqrts_all[i_a][0];
  logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << endl;
  logger << LOG_DEBUG_VERBOSE << "dipole[" << i_a << "].name() = " << dipole[i_a].name() << endl;
  generic->ac_psp_dipole(i_a, MC_phasespace.channel - dipole[i_a - 1].sum_channel(), *this);
  logger << LOG_DEBUG_VERBOSE << "3" << endl;
  double two_ptk_pa = 2 * pt_k * p_a;
  //  double sqrt_two_ptk_pa = sqrt(two_ptk_pa);
  //  double mu_k = m_k / sqrt_two_ptk_pa;
  double mu2_k = m2_k / two_ptk_pa;
  double uplus = one_minus_x_ik_a / (one_minus_x_ik_a + mu2_k);
  double u_i = h_propto_pot_mod(r[no_random_dipole[1]], 0., uplus, exp_ai_k_u, map_technical_x);//exp_ij_a());
  double one_minus_u_i = 1. - u_i;
  double phi_a_perp = c_phi(no_random_dipole[2]);
  fourvector Q = xbp_all[i_a][bpt_k] + xbp_all[i_a][bp_a];
  double Q2 = Q.m2();
  double sqrtQ2 = sqrt(Q2);
  fourvector temp_p_a_boost = xbp_all[i_a][bp_a].boost(Q);
  // not used !   fourvector Q_phattilde_ai((Q2 - m2_k) / (2 * sqrtQ2), 0., 0., (Q2 - m2_k) / (2 * sqrtQ2));
  // not used !   fourvector Q_phattilde_k((Q2 + m2_k) / (2 * sqrtQ2), 0., 0., -(Q2 - m2_k) / (2 * sqrtQ2));
  // not used !   fourvector Q_phat_a = (1. / x_ik_a) * Q_phattilde_ai;
  double help_p_i = (one_minus_u_i * one_minus_x_ik_a * (Q2 - m2_k) - u_i * x_ik_a * m2_k) / (2 * x_ik_a * sqrtQ2);
  double help_p_i_T = sqrt((u_i * one_minus_u_i * one_minus_x_ik_a / x_ik_a) * (Q2 - m2_k) - pow(u_i, 2) * m2_k);
  fourvector Q_phat_i(help_p_i + .5 * u_i * sqrtQ2, help_p_i_T * cos(phi_a_perp), help_p_i_T * sin(phi_a_perp), help_p_i - .5 * u_i * sqrtQ2);
  fourvector Q_phat_ik((Q2 - (one_minus_x_ik_a - x_ik_a) * m2_k) / (2 * x_ik_a * sqrtQ2), 0., 0., (one_minus_x_ik_a - x_ik_a) * (Q2 - m2_k) / (2 * x_ik_a * sqrtQ2));
  fourvector Q_phat_k = Q_phat_ik - Q_phat_i;
  for (int xbi = 4; xbi < xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){xbp_all[0][bp_i] = (Q_phat_i.rotateback(temp_p_a_boost)).boost(Q.Pinv());}
    else if (xbi == bp_k){xbp_all[0][bp_k] = (Q_phat_k.rotateback(temp_p_a_boost)).boost(Q.Pinv());}
    else if (xbi < bp_i){xbp_all[0][xbi] = xbp_all[i_a][xbi];}
    else if (xbi > bp_i){xbp_all[0][xbi] = xbp_all[i_a][xbi / 2];}
    else {logger << LOG_INFO << "Should not happen!" << endl;}
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::ag_psp_RA_group_ai_k_massive(double & sinx_ik_a_min, vector<dipole_set> & dipole, int i_a){
  static Logger logger("phasespace_set::ag_psp_RA_group_ai_k_massive");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int bp_a = dipole[i_a].binary_R_emitter_1(); // bp_b not needed !!!
  int bp_i = dipole[i_a].binary_R_emitter_2();
  int bp_k = dipole[i_a].binary_R_spectator();
  int bpt_k = dipole[i_a].binary_A_spectator();
  if (xbs_all[0][bp_i + bp_k] == 0.){xbs_all[0][bp_i + bp_k] = (xbp_all[0][bp_i] + xbp_all[0][bp_k]).m2();}
  double pi_pk = p_i * p_k;
  double pi_pa = p_i * p_a;
  double pk_pa = p_k * p_a;
  double two_ptk_pa = 2 * pt_k * p_a;
  //  double sqrt_two_ptk_pa = sqrt(two_ptk_pa);
  //  double mu_k = m_k / sqrt_two_ptk_pa;
  double mu2_k = m2_k / two_ptk_pa;
  double one_minus_x_ik_a_min = 1. - sinx_ik_a_min / xbs_all[0][0];
  double one_minus_x_ik_a = pi_pk / (pi_pa + pk_pa);
  double x_ik_a = 1. - one_minus_x_ik_a;
  double u_i = pi_pa / (pi_pa + pk_pa);
  double one_minus_u_i = 1. - u_i;
  double uplus = one_minus_x_ik_a / (one_minus_x_ik_a + mu2_k);
  double g_alpha = g_propto_pot_mod(one_minus_x_ik_a, 0., one_minus_x_ik_a_min, exp_ai_k_x, map_technical_x) * g_propto_pot_mod(u_i, 0., uplus, exp_ai_k_u, map_technical_x) / (pi * (xbp_all[i_a][bpt_k] * xbp_all[0][bp_a]));
  if (xbs_all[i_a][0] == 0.){xbs_all[i_a][0] = x_ik_a * xbs_all[0][0];}
  if (xbsqrts_all[i_a][0] == 0.){xbsqrts_all[i_a][0] = sqrt(xbs_all[i_a][0]);}
  generic->ag_psp_dipole(i_a, dipole[i_a - 1].sum_channel(), *this);
  for (int i = dipole[i_a - 1].sum_channel(); i < dipole[i_a].sum_channel(); i++){MC_phasespace.g_channel[i] = MC_phasespace.g_channel[i] * g_alpha;}
  // to be adapted to massive case !!!
  // invert mappings for VAMP implementation
  if (switch_IS_mode_phasespace == 2 || switch_IS_mode_phasespace == 4){
    //  static vector<double> inv_r(3);
    vector<double> inv_r(3);
    inv_r[0] = inv_propto_pot(one_minus_x_ik_a, 0., one_minus_x_ik_a_min, exp_ai_k_x, map_technical_x);
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
#undef p_i 
#undef p_a
#undef p_k 
#undef m_k
#undef m2_k
