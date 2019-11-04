#include "../include/classes.cxx"

void phasespace_set::calculate_initial_tau_IS_x1x2_IS(){
  Logger logger("phasespace_set::calculate_initial_tau_IS_x1x2_IS");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  randomvector(sran, 1, random_MC_tau);
  randomvector(sran, 1, random_tau);
  for (int i = 0; i < random_tau.size(); i++){
    logger << LOG_DEBUG_VERBOSE << "random_tau[" << i << "] = " << random_tau[i] << endl;
  }
  for (int j = 0; j < MC_tau.beta.size(); j++){if (random_MC_tau[0] <= MC_tau.beta[j]){MC_tau.channel = j; break;}}
  logger << LOG_DEBUG_VERBOSE << "MC_tau.channel = " << MC_tau.channel << endl;
  if (MC_tau.channel == 0){
    /*
    // old...
    for (int j = 0; j < tau_beta.size(); j++){if (random_tau[0] <= tau_beta[j]){tau_channel = j; break;}}
    double old_random_IS_tau = (double(tau_channel) + random_tau[1]) / double(tau_beta.size());
    // ...old
    */
    for (int j = 0; j < IS_tau.beta.size(); j++){if (random_tau[0] <= IS_tau.beta[j]){IS_tau.channel = j; break;}}
    double random_IS_tau = (double(IS_tau.channel) + random_tau[1]) / double(IS_tau.beta.size());
    /*
    assert(tau_channel == IS_tau.channel);
    assert(old_random_IS_tau == random_IS_tau);
    */
    x_pdf[0] = h_propto_pot(random_IS_tau, tau_0, exp_pdf);
  }
  else {
    x_pdf[0] = c_propagator_Breit_Wigner(random_tau[0], abs(tau_MC_map[MC_tau.channel]), tau_0_s_had, s_had) / s_had;
    logger << LOG_DEBUG_VERBOSE << setw(20) << "x_pdf[0]" << " = " << setprecision(20) << setw(28) << x_pdf[0] << "   " << double2hexastr(x_pdf[0]) << endl;
    /*
    for (int j = 0; j < tau_MC_tau_gamma[0].size(); j++){
      logger << LOG_DEBUG_VERBOSE << "tau_MC_tau_gamma[0][" << j << "]" << " = " << setprecision(20) << setw(28) << tau_MC_tau_gamma[0][j] << "   " << double2hexastr(tau_MC_tau_gamma[0][j]) << endl;
    }
    */
    
    for (int j = 0; j < tau_MC_tau_gamma[0].size(); j++){
      if (x_pdf[0] <= tau_MC_tau_gamma[0][j]){
	//	tau_channel = j;
	IS_tau.channel = j;
	break;
      }
    }
    logger << LOG_DEBUG_VERBOSE << "IS_tau.channel = " << IS_tau.channel << endl;
  }
  logger << LOG_DEBUG_VERBOSE << "x_pdf[0] = " << x_pdf[0] << endl;
  //  logger << LOG_DEBUG_VERBOSE << "tau_channel = " << tau_channel << endl;

  double tau_s_had = x_pdf[0] * s_had;
  // g_tau_MC could be replaced by MC_tau.g_channel
  //  vector<double> g_tau_MC(MC_tau.alpha.size());
  for (int j = 0; j < MC_tau.g_channel.size(); j++){
    if (j == 0){
      /*
      MC_tau.g_channel[j] = g_propto_pot(x_pdf[0], tau_0, exp_pdf) * (tau_alpha[tau_channel] * tau_alpha.size());
      double temp = MC_tau.g_channel[j];
      */
      MC_tau.g_channel[j] = g_propto_pot(x_pdf[0], tau_0, exp_pdf) * (IS_tau.alpha[IS_tau.channel] * IS_tau.alpha.size());
      //      assert(temp == MC_tau.g_channel[j]);
    }
    else {
      if (tau_MC_map[j] < 0){MC_tau.g_channel[j] = s_had * g_propagator_Breit_Wigner(tau_s_had, -tau_MC_map[j], tau_0_s_had, s_had);}
    }
  }
  g_tau = 0.;
  for (int i = 0; i < MC_tau.g_channel.size(); i++){g_tau += MC_tau.alpha[i] * MC_tau.g_channel[i];}

  for (int i = 0; i < MC_tau.g_channel.size(); i++){
    logger << LOG_DEBUG_VERBOSE << setw(20) << "MC_tau.alpha[" << i << "]" << " = " << setprecision(20) << setw(28) << MC_tau.alpha[i] << "   " << double2hexastr(MC_tau.alpha[i]) << endl;
  }
  logger << LOG_DEBUG_VERBOSE << "IS_tau.channel = " << IS_tau.channel << endl;
  logger << LOG_DEBUG_VERBOSE << setw(20) << "IS_tau.alpha[" << IS_tau.channel << "]" << " = " << setprecision(20) << setw(28) << IS_tau.alpha[IS_tau.channel] << "   " << double2hexastr(IS_tau.alpha[IS_tau.channel]) << endl;
  for (int i = 0; i < MC_tau.g_channel.size(); i++){
    logger << LOG_DEBUG_VERBOSE << setw(20) << "MC_tau.g_channel[" << i << "]" << " = " << setprecision(20) << setw(28) << MC_tau.g_channel[i] << "   " << double2hexastr(MC_tau.g_channel[i]) << endl;
  }

  logger << LOG_DEBUG_VERBOSE << setw(20) << "g_tau" << " = " << setprecision(20) << setw(28) << g_tau << "   " << double2hexastr(g_tau) << endl;

  //  vector<double> z_x1x2(2);
  randomvector(sran, 1, random_x12);
  // old...
  /*
  for (int j = 0; j < x1x2_beta.size(); j++){if (random_x12[0] <= x1x2_beta[j]){x1x2_channel = j; break;}}
  double old_random_x1x2 = (double(x1x2_channel) + random_x12[1]) / double(x1x2_beta.size());
  // ...old
  */
  for (int j = 0; j < IS_x1x2.beta.size(); j++){if (random_x12[0] <= IS_x1x2.beta[j]){IS_x1x2.channel = j; break;}}
  double random_x1x2 = (double(IS_x1x2.channel) + random_x12[1]) / double(IS_x1x2.beta.size());
  /*
  assert(x1x2_channel == IS_x1x2.channel);
  assert(old_random_x1x2 == random_x1x2);
  */
  double min_x_pdf_1 = GSL_MAX_DBL(x_pdf[0], 1.e-7);
  x_pdf[1] = exp(random_x1x2 * log(min_x_pdf_1));
  //  x_pdf[1] = exp(random_x1x2 * log(x_pdf[0]));
  x_pdf[2] = x_pdf[0] / x_pdf[1];
  /*
  g_x1x2 = 1. / (-log(min_x_pdf_1)) * (x1x2_alpha[x1x2_channel] * x1x2_alpha.size());
  //  g_x1x2 = 1. / (-log(x_pdf[0])) * (x1x2_alpha[x1x2_channel] * x1x2_alpha.size());
  double temp = g_x1x2;
  */
  g_x1x2 = 1. / (-log(min_x_pdf_1)) * (IS_x1x2.alpha[IS_x1x2.channel] * IS_x1x2.alpha.size());
  //  assert(temp == g_x1x2);
  
  logger << LOG_DEBUG_VERBOSE << setw(20) << "g_x1x2" << " = " << setprecision(20) << setw(28) << g_x1x2 << "   " << double2hexastr(g_x1x2) << endl;

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



void phasespace_set::calculate_initial_tau_IS_x1x2(){
  Logger logger("phasespace_set::calculate_initial_tau_IS_x1x2");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  randomvector(sran, 0, random_MC_tau);
  randomvector(sran, 1, random_tau);
  for (int j = 0; j < MC_tau.beta.size(); j++){if (random_MC_tau[0] <= MC_tau.beta[j]){MC_tau.channel = j; break;}}
  if (MC_tau.channel == 0){
    for (int j = 0; j < IS_tau.beta.size(); j++){if (random_tau[0] <= IS_tau.beta[j]){IS_tau.channel = j; break;}}
    double random_IS_tau = (double(IS_tau.channel) + random_tau[1]) / double(IS_tau.beta.size());
    x_pdf[0] = h_propto_pot(random_IS_tau, tau_0, exp_pdf);
  }
  else {
    x_pdf[0] = c_propagator_Breit_Wigner(random_tau[0], abs(tau_MC_map[MC_tau.channel]), tau_0_s_had, s_had) / s_had;
    for (int j = 0; j < tau_MC_tau_gamma[0].size(); j++){
      if (x_pdf[0] <= tau_MC_tau_gamma[0][j]){
	IS_tau.channel = j; 
	break;
      }
    }
  }
  double tau_s_had = x_pdf[0] * s_had;
//  vector<double> g_tau_MC(MC_tau.alpha.size());
  for (int j = 0; j < MC_tau.g_channel.size(); j++){
    if (j == 0){
      MC_tau.g_channel[j] = g_propto_pot(x_pdf[0], tau_0, exp_pdf) * (IS_tau.alpha[IS_tau.channel] * IS_tau.alpha.size());
    }
    else {
      if (tau_MC_map[j] < 0){MC_tau.g_channel[j] = s_had * g_propagator_Breit_Wigner(tau_s_had, -tau_MC_map[j], tau_0_s_had, s_had);}
      ///      else if (tau_MC_map[j] > 0){MC_tau.g_channel[j] = s_had * g_propagator_vanishing_width(tau_s_had, M2[tau_MC_map[j]], tau_0_s_had, s_had, nuxs);}
    }
  }
  g_tau = 0.;
  for (int i = 0; i < MC_tau.g_channel.size(); i++){g_tau += MC_tau.alpha[i] * MC_tau.g_channel[i];}

  //  vector<double> z_x1x2(2);
  randomvector(sran, 1, random_x12);
  x_pdf[1] = exp(random_x12[0] * log(x_pdf[0]));
  x_pdf[2] = x_pdf[0] / x_pdf[1];
  g_x1x2 = 1. / (-log(x_pdf[0]));

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



void phasespace_set::calculate_initial_tau_x1x2_IS(){
  Logger logger("phasespace_set::calculate_initial_tau_x1x2_IS");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // adapt to switch off tau IS mapping !!!
  randomvector(sran, 0, random_MC_tau);
  randomvector(sran, 1, random_tau);
  
  for (int j = 0; j < MC_tau.beta.size(); j++){if (random_MC_tau[0] <= MC_tau.beta[j]){MC_tau.channel = j; break;}}
  if (MC_tau.channel == 0){
    x_pdf[0] = h_propto_pot(random_tau[0], tau_0, exp_pdf);
  }
  else {
    x_pdf[0] = c_propagator_Breit_Wigner(random_tau[0], abs(tau_MC_map[MC_tau.channel]), tau_0_s_had, s_had) / s_had;
  }
  double tau_s_had = x_pdf[0] * s_had;
  for (int j = 0; j < MC_tau.g_channel.size(); j++){
    if (j == 0){
      MC_tau.g_channel[j] = g_propto_pot(x_pdf[0], tau_0, exp_pdf);
    }
    else {
      if (tau_MC_map[j] < 0){MC_tau.g_channel[j] = s_had * g_propagator_Breit_Wigner(tau_s_had, -tau_MC_map[j], tau_0_s_had, s_had);}
    }
  }
  g_tau = 0.;
  for (int i = 0; i < MC_tau.g_channel.size(); i++){g_tau += MC_tau.alpha[i] * MC_tau.g_channel[i];}

  //  vector<double> z_x1x2(2);
  randomvector(sran, 1, random_x12);
  for (int j = 0; j < IS_x1x2.beta.size(); j++){if (random_x12[0] <= IS_x1x2.beta[j]){IS_x1x2.channel = j; break;}}
  double random_x1x2 = (double(IS_x1x2.channel) + random_x12[1]) / double(IS_x1x2.beta.size());
  x_pdf[1] = exp(random_x1x2 * log(x_pdf[0]));
  x_pdf[2] = x_pdf[0] / x_pdf[1];
  g_x1x2 = 1. / (-log(x_pdf[0])) * (IS_x1x2.alpha[IS_x1x2.channel] * IS_x1x2.alpha.size());

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



void phasespace_set::calculate_initial_tau_x1x2(){
  Logger logger("phasespace_set::calculate_initial_tau_x1x2");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  vector<double> z_tau_MC(1);
  randomvector(sran, 0, random_MC_tau);
  //  vector<double> z_tau(2);
  randomvector(sran, 1, random_tau);
  for (int j = 0; j < MC_tau.beta.size(); j++){if (random_MC_tau[0] <= MC_tau.beta[j]){MC_tau.channel = j; break;}}
  if (MC_tau.channel == 0){
    x_pdf[0] = h_propto_pot(random_tau[0], tau_0, exp_pdf);
  }
  else {
    x_pdf[0] = c_propagator_Breit_Wigner(random_tau[0], abs(tau_MC_map[MC_tau.channel]), tau_0_s_had, s_had) / s_had;
    //    if (tau_MC_map[MC_tau.channel] < 0){x_pdf[0] = c_propagator_Breit_Wigner(random_tau[0], abs(tau_MC_map[MC_tau.channel]), tau_0_s_had, s_had) / s_had;}
    ///    else if (tau_MC_map[MC_tau.channel] > 0){x_pdf[0] = c_propagator_vanishing_width(random_tau[0], M2[tau_MC_map[MC_tau.channel]], tau_0_s_had, s_had, nuxs) / s_had;}
  }
  double tau_s_had = x_pdf[0] * s_had;
//  vector<double> g_tau_MC(MC_tau.alpha.size());
  for (int j = 0; j < MC_tau.g_channel.size(); j++){
    if (j == 0){
      MC_tau.g_channel[j] = g_propto_pot(x_pdf[0], tau_0, exp_pdf);
    }
    else {
      if (tau_MC_map[j] < 0){MC_tau.g_channel[j] = s_had * g_propagator_Breit_Wigner(tau_s_had, -tau_MC_map[j], tau_0_s_had, s_had);}
      ///      else if (tau_MC_map[j] > 0){MC_tau.g_channel[j] = s_had * g_propagator_vanishing_width(tau_s_had, M2[tau_MC_map[j]], tau_0_s_had, s_had, nuxs);}
    }
  }
  g_tau = 0.;
  for (int i = 0; i < MC_tau.g_channel.size(); i++){g_tau += MC_tau.alpha[i] * MC_tau.g_channel[i];}

  //  vector<double> z_x1x2(2);
  randomvector(sran, 1, random_x12);
  x_pdf[1] = exp(random_x12[0] * log(x_pdf[0]));
  x_pdf[2] = x_pdf[0] / x_pdf[1];
  g_x1x2 = 1. / (-log(x_pdf[0]));

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

void phasespace_set::calculate_initial_tau_fixed_x1x2_IS(){
  Logger logger("phasespace_set::calculate_initial_tau_fixed_x1x2_IS");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  x_pdf[0] = start_xbs_all[0][xb_max - 4] / s_had;
  g_tau = s_had;
  //  if (csi->n_particle == 1){g_tau = s_had;}
  //  else {g_tau = s_had / (pi * M[24] * map_Gamma[24]);} // 24 must be replaced by process-dependent number !!!

  randomvector(sran, 1, random_x12);
  for (int j = 0; j < IS_x1x2.beta.size(); j++){if (random_x12[0] <= IS_x1x2.beta[j]){IS_x1x2.channel = j; break;}}
  double random_x1x2 = (double(IS_x1x2.channel) + random_x12[1]) / double(IS_x1x2.beta.size());
  x_pdf[1] = exp(random_x1x2 * log(x_pdf[0]));
  x_pdf[2] = x_pdf[0] / x_pdf[1];
  g_x1x2 = 1. / (-log(x_pdf[0])) * (IS_x1x2.alpha[IS_x1x2.channel] * IS_x1x2.alpha.size());

  double sqrttau = sqrt(x_pdf[0]);
  xbp_all[0][1] = fourvector(sqrttau * E, 0., 0., sqrttau * E);
  xbp_all[0][2] = fourvector(sqrttau * E, 0., 0., -sqrttau * E);
  xbp_all[0][0] = xbp_all[0][1] + xbp_all[0][2];
  xbs_all[0][0] = 4. * x_pdf[0] * pow(E, 2);
  xbsqrts_all[0][0] = sqrt(xbs_all[0][0]);
  xbp_all[0][xb_max - 4] = xbp_all[0][0];
  xbs_all[0][xb_max - 4] = xbs_all[0][0];
  xbsqrts_all[0][xb_max - 4] = xbsqrts_all[0][0];

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void phasespace_set::calculate_initial_tau_fixed_x1x2(){
  Logger logger("phasespace_set::calculate_initial_tau_fixed_x1x2");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  x_pdf[0] = start_xbs_all[0][xb_max - 4] / s_had;
  g_tau = s_had;
  //  if (csi->n_particle == 1){g_tau = s_had;}
  //  else {g_tau = s_had / (pi * M[24] * map_Gamma[24]);} // 24 must be replaced by process-dependent number !!!

  randomvector(sran, 1, random_x12);
  x_pdf[1] = exp(random_x12[0] * log(x_pdf[0]));
  x_pdf[2] = x_pdf[0] / x_pdf[1];
  g_x1x2 = 1. / (-log(x_pdf[0]));

  double sqrttau = sqrt(x_pdf[0]);
  xbp_all[0][1] = fourvector(sqrttau * E, 0., 0., sqrttau * E);
  xbp_all[0][2] = fourvector(sqrttau * E, 0., 0., -sqrttau * E);
  xbp_all[0][0] = xbp_all[0][1] + xbp_all[0][2];
  xbs_all[0][0] = 4. * x_pdf[0] * pow(E, 2);
  xbsqrts_all[0][0] = sqrt(xbs_all[0][0]);
  xbp_all[0][xb_max - 4] = xbp_all[0][0];
  xbs_all[0][xb_max - 4] = xbs_all[0][0];
  xbsqrts_all[0][xb_max - 4] = xbsqrts_all[0][0];

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::calculate_initial_collinear_z1z2_IS(){
  Logger logger("phasespace_set::calculate_initial_collinear_z1z2_IS");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_z = 1; i_z < 3; i_z++){
    vector<double> z_z1z2(2);
    randomvector(sran, 1, z_z1z2);
    for (int j = 0; j < IS_z1z2[i_z].beta.size(); j++){
      if (z_z1z2[0] <= IS_z1z2[i_z].beta[j]){
	IS_z1z2[i_z].channel = j; 
	//	z1z2_channel[i_z] = j; 
	break;
      }
    }
    double random_z1z2 = (double(IS_z1z2[i_z].channel) + z_z1z2[1]) / double(IS_z1z2[i_z].beta.size());
    z_coll[i_z] = exp(random_z1z2 * log(x_pdf[i_z]));
    g_z_coll[i_z] = 1. / (-log(x_pdf[i_z])) * (IS_z1z2[i_z].alpha[IS_z1z2[i_z].channel] * IS_z1z2[i_z].alpha.size());
  }
  /*
  for (int i_z = 1; i_z < 3; i_z++){
    vector<double> z_z1z2(2);
    randomvector(sran, 1, z_z1z2);
    for (int j = 0; j < z1z2_beta[i_z].size(); j++){
      if (z_z1z2[0] <= z1z2_beta[i_z][j]){
	z1z2_channel[i_z] = j; 
	break;
      }
    }
    double random_z1z2 = (double(z1z2_channel[i_z]) + z_z1z2[1]) / double(z1z2_beta[i_z].size());
    z_coll[i_z] = exp(random_z1z2 * log(x_pdf[i_z]));
    g_z_coll[i_z] = 1. / (-log(x_pdf[i_z])) * (z1z2_alpha[i_z][z1z2_channel[i_z]] * z1z2_alpha[i_z].size());
  }
   */
  for (int i_x = 0; i_x < 3; i_x++){
    all_xz_coll_pdf[i_x][0] = x_pdf[i_x];
  }
  for (int i_z = 1; i_z < 3; i_z++){
    all_xz_coll_pdf[0][i_z] = z_coll[i_z];
    all_xz_coll_pdf[i_z][i_z] = x_pdf[i_z] / z_coll[i_z];
    all_xz_coll_pdf[i_z % 2 + 1][i_z] = x_pdf[0] / z_coll[i_z];
  }

  // all_xz_coll_pdf: 
  //  0 - 0 : tau = x1 * x2
  //  1 - 0 : x1 -> enters 1st pdf if ()_+ or delta distributions are present
  //  2 - 0 : x2 -> enters 2nd pdf if ()_+ or delta distributions are present
  //  0 - 1 : z1 -> collinear emission from a
  //  1 - 1 : x1 / z1 -> enters 1st pdf for regular and ()_+ distributions
  //  2 - 1 : (x1 / z1) * x2 -> devide LHAPDF value by this to get the pdf
  //  0 - 2 : z2 -> collinear emission from b
  //  1 - 2 : x1 * (x2 / z2) -> devide LHAPDF value by this to get the pdf
  //  2 - 2 : x2 / z2 -> enters 1st pdf for regular and ()_+ distributions

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
//#include "old.phasespace.mapping.initial.cpp"
