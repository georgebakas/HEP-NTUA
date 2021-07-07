#include "../include/classes.cxx"

void observable_set::calculate_dipole_Acc_QCD(int x_a, double & Dfactor){
  static Logger logger("observable_set::calculate_dipole_Acc_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (switch_OL == 1){

  static int n_momentum = 5 * (p_parton[x_a].size() - 1);
  for (int i = 1; i < p_parton[x_a].size(); i++){logger << LOG_DEBUG_VERBOSE << "p_parton[" << x_a << "][" << i << "] = " << p_parton[x_a][i] << endl;}
  logger << LOG_DEBUG_VERBOSE << "n_momentum = " << n_momentum << endl;
  static double *P;
  P = new double[n_momentum];
  for (int i = 1; i < p_parton[x_a].size(); i++){
    P[5 * (i - 1)]     = p_parton[x_a][i].x0();
    P[5 * (i - 1) + 1] = p_parton[x_a][i].x1();
    P[5 * (i - 1) + 2] = p_parton[x_a][i].x2();
    P[5 * (i - 1) + 3] = p_parton[x_a][i].x3();
    //    P[5 * (i - 1) + 4] = p_parton[x_a][i].m();
    P[5 * (i - 1) + 4] = 0.;
  }
  for (int i = 0; i < n_momentum; i++){logger << LOG_DEBUG_VERBOSE << "P[" << i << "] = " << P[i] << endl;}
  static int n_cc = (p_parton[x_a].size() - 1) * (p_parton[x_a].size() - 2) / 2;
  logger << LOG_DEBUG_VERBOSE << "n_cc = " << n_cc << endl;
  static double b_ME2;
  static double ewcc;
  static double *M2cc;
  M2cc = new double[n_cc];
  logger << LOG_DEBUG_VERBOSE << "(*RA_dipole)[x_a].process_id = " << (*RA_dipole)[x_a].process_id << endl;
  if ((csi->type_contribution == "RA" ||
       csi->type_contribution == "RRA") && user.string_value[user.string_map["model"]] != "Bornloop"){
    ol_evaluate_cc((*RA_dipole)[x_a].process_id, P, &b_ME2, M2cc, &ewcc);
  }
  else if (csi->type_contribution == "L2RA" || 
	   user.string_value[user.string_map["model"]] == "Bornloop"){
    ol_evaluate_cc2((*RA_dipole)[x_a].process_id, P, &b_ME2, M2cc, &ewcc);
  }

  logger << LOG_DEBUG_VERBOSE << "(*RA_dipole)[x_a].process_id = " << (*RA_dipole)[x_a].process_id << endl;
  for (int i = 0; i < n_cc; i++){logger << LOG_DEBUG_VERBOSE << "M2cc[" << i << "] = " << M2cc[i] << endl;}
  Dfactor = M2cc[(*RA_dipole)[x_a].no_BLHA_entry];
  logger << LOG_DEBUG_POINT << "OpenLoops:     Dfactor[" << x_a << "] = " << setw(23) << setprecision(15) << Dfactor << endl;
  logger << LOG_DEBUG_VERBOSE << "(*RA_dipole)[x_a].no_BLHA_entry = " << (*RA_dipole)[x_a].no_BLHA_entry << endl;
  delete [] M2cc;
  delete [] P;

  } 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_dipole_Asc_QCD(int x_a, fourvector & Vtensor, double & ME2_metric, double & ME2_vector){
  static Logger logger("observable_set::calculate_dipole_Asc_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (switch_OL == 1){

  static int n_momentum = 5 * (p_parton[x_a].size() - 1);
  static double *P;
  P = new double[n_momentum];
  for (int i = 1; i < p_parton[x_a].size(); i++){
    P[5 * (i - 1)]     = p_parton[x_a][i].x0();
    P[5 * (i - 1) + 1] = p_parton[x_a][i].x1();
    P[5 * (i - 1) + 2] = p_parton[x_a][i].x2();
    P[5 * (i - 1) + 3] = p_parton[x_a][i].x3();
    //    P[5 * (i - 1) + 4] = p_parton[x_a][i].m();
    P[5 * (i - 1) + 4] = 0.;
  }
  double *MOM;
  MOM = new double[4];
  MOM[0] = Vtensor.x0();
  MOM[1] = Vtensor.x1();
  MOM[2] = Vtensor.x2();
  MOM[3] = Vtensor.x3();
  static int n_cc = (p_parton[x_a].size() - 1) * (p_parton[x_a].size() - 2) / 2;
  logger << LOG_DEBUG_VERBOSE << "n_cc = " << n_cc << endl;
  static double b_ME2;
  static double ewcc;
  static double *M2cc;
  M2cc = new double[n_cc];

  if ((csi->type_contribution == "RA" ||
       csi->type_contribution == "RRA") && user.string_value[user.string_map["model"]] != "Bornloop"){
    ol_evaluate_cc((*RA_dipole)[x_a].process_id, P, &b_ME2, M2cc, &ewcc);
  }
  else if (csi->type_contribution == "L2RA" || 
	   user.string_value[user.string_map["model"]] == "Bornloop"){
    logger << LOG_DEBUG_VERBOSE << "before ol_evaluate_cc2" << endl;
    ol_evaluate_cc2((*RA_dipole)[x_a].process_id, P, &b_ME2, M2cc, &ewcc);
    logger << LOG_DEBUG_VERBOSE << "after ol_evaluate_cc2" << endl;
  }

  ME2_metric = M2cc[(*RA_dipole)[x_a].no_BLHA_entry];
  logger << LOG_DEBUG_POINT << "OpenLoops:     ME2_metric[" << x_a << "] = " << setw(23) << setprecision(15) << ME2_metric << endl;
  static int n_sc = p_parton[x_a].size() - 1;
  logger << LOG_DEBUG_VERBOSE << "n_sc = " << n_sc << endl;
  static double *M2sc;
  M2sc = new double[n_sc];

  if ((csi->type_contribution == "RA" ||
       csi->type_contribution == "RRA") && user.string_value[user.string_map["model"]] != "Bornloop"){
    ol_evaluate_sc((*RA_dipole)[x_a].process_id, P, (*RA_dipole)[x_a].no_A_emitter(), MOM, M2sc);
  }
  else if (csi->type_contribution == "L2RA" || 
	   user.string_value[user.string_map["model"]] == "Bornloop"){
    logger << LOG_DEBUG_VERBOSE << "before ol_evaluate_sc2" << endl;
    ol_evaluate_sc2((*RA_dipole)[x_a].process_id, P, (*RA_dipole)[x_a].no_A_emitter(), MOM, M2sc);
    logger << LOG_DEBUG_VERBOSE << "after ol_evaluate_sc2" << endl;
  }

  ME2_vector = M2sc[(*RA_dipole)[x_a].no_A_spectator() - 1];
  logger << LOG_DEBUG_POINT << "OpenLoops:     ME2_vector[" << x_a << "] = " << setw(23) << setprecision(15) << ME2_vector << endl;

  logger << LOG_DEBUG_VERBOSE << "Vtensor = " << Vtensor << "   M2sc[" << (*RA_dipole)[x_a].no_A_spectator() - 1 << "] = " << ME2_vector << endl;
  for (int i = 0; i < n_sc; i++){logger << LOG_DEBUG_VERBOSE << "M2sc[" << i << "] = " << M2sc[i] << endl;}

  delete [] M2sc;
  delete [] M2cc;
  delete [] MOM;
  delete [] P;

  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



#define p_i p_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define p_j p_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define p_k p_parton[0][(*RA_dipole)[x_a].no_R_spectator()]
#define pt_ij p_parton[x_a][(*RA_dipole)[x_a].no_A_emitter()]
#define pt_k p_parton[x_a][(*RA_dipole)[x_a].no_A_spectator()]
double observable_set::calculate_dipole_QCD_A_ij_k(int x_a){
  static Logger logger("observable_set::calculate_dipole_A_ij_k");
  logger << LOG_DEBUG_VERBOSE << "called with dipole no. " << x_a << endl;
  double ME2 = 0.;
  //  int dipole_no_A_emitter = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_emitter()];
  //  int dipole_no_A_spectator = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_spectator()];
  double pi_pj = p_i * p_j;
  double pj_pk = p_j * p_k;
  double pi_pk = p_i * p_k;
  double y_ij_k = pi_pj / (pi_pj + pi_pk + pj_pk);
  double one_minus_y_ij_k = 1. - y_ij_k;
  double z_i = pi_pk / (pi_pk + pj_pk);
  double z_j = pj_pk / (pi_pk + pj_pk);
  if ((*RA_dipole)[x_a].type_splitting() < 2){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = z_i * p_i - z_j * p_j;
    if ((*RA_dipole)[x_a].type_splitting() == 0){
      // emitter1 = gluon, emitter2 = gluon
      factor = -1. / (2 * pi_pj) * 16 * C_A * pi * alpha_S;
      factor_vector = 1. / pi_pj;
      factor_metric = 1. / (1. - z_i * one_minus_y_ij_k) + 1. / (1 - z_j * one_minus_y_ij_k) - 2.;
    }
    else if ((*RA_dipole)[x_a].type_splitting() == 1){
      // emitter1 = quark, emitter2 = antiquark
      factor = -1. / (2 * pi_pj) * 8 * T_R * pi * alpha_S;
      factor_vector = -2. / pi_pj;
      factor_metric = 1.;
    }

    calculate_dipole_Asc_QCD(x_a, Vtensor, ME2_metric, ME2_vector);

    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) / C_A * (*RA_dipole)[x_a].symmetry_factor();
  }
  else {
    // always (*RA_dipole)[x_a].type_splitting() == 2; emitter1 = (anti)quark, emitter2 = gluon
    double Dfactor = 0.;

    calculate_dipole_Acc_QCD(x_a, Dfactor);

    Dfactor = Dfactor / C_F * (*RA_dipole)[x_a].symmetry_factor();

    ME2 = -1. / (2. * pi_pj) * (8 * C_F * pi * alpha_S) * (2. / (1. - z_i * one_minus_y_ij_k) - (1. + z_i)) * Dfactor;
  }
  return -ME2;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

#define m2_i mass2_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define m2_j mass2_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define m2_k mass2_parton[0][(*RA_dipole)[x_a].no_R_spectator()]
#define m2_ij mass2_parton[x_a][(*RA_dipole)[x_a].no_A_emitter()]
double observable_set::calculate_dipole_QCD_A_ij_k_massive(int x_a){
  static Logger logger("observable_set::calculate_dipole_QCD_A_ij_k_massive");
  logger << LOG_DEBUG_VERBOSE << "called with dipole no. " << x_a << endl;
  logger << LOG_DEBUG_VERBOSE << "m2_i = " << m2_i << endl;
  logger << LOG_DEBUG_VERBOSE << "m2_j = " << m2_j << endl;
  logger << LOG_DEBUG_VERBOSE << "m2_k = " << m2_k << endl;
  logger << LOG_DEBUG_VERBOSE << "m2_ij = " << m2_ij << endl;
  double ME2 = 0.;
  //  int dipole_no_A_emitter = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_emitter()];
  //  int dipole_no_A_spectator = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_spectator()];
  double pi_pj = p_i * p_j;
  double pj_pk = p_j * p_k;
  double pi_pk = p_i * p_k;
  double y_ij_k = pi_pj / (pi_pj + pi_pk + pj_pk);
  double one_minus_y_ij_k = 1. - y_ij_k;
  double z_i = pi_pk / (pi_pk + pj_pk);
  double z_j = pj_pk / (pi_pk + pj_pk);
  fourvector Q = pt_ij + pt_k;
  double Q2 = Q.m2();
  double mu2_i = m2_i / Q2;
  double mu2_j = m2_j / Q2;
  double mu2_k = m2_k / Q2;
  double mu2_ij = m2_ij / Q2;
  double one_minus_all_mu = 1. - mu2_i - mu2_j - mu2_k;
  double z_prefactor = (2 * mu2_i + one_minus_all_mu * y_ij_k) / (2 * (mu2_i + mu2_j + one_minus_all_mu * y_ij_k));
  double v_ij_k = sqrt(pow(2 * mu2_k + one_minus_all_mu * one_minus_y_ij_k, 2) - 4 * mu2_k) / (one_minus_all_mu * one_minus_y_ij_k);
  double v_ij_i = sqrt(pow(one_minus_all_mu * y_ij_k, 2) - 4 * mu2_i * mu2_j) / (one_minus_all_mu * y_ij_k + 2 * mu2_i);
  double vt_ij_k = sqrt(lambda(1, mu2_ij, mu2_k)) / (1. - mu2_ij - mu2_k);
  double z_minus = z_prefactor * (1 - v_ij_i * v_ij_k);
  double z_plus = z_prefactor * (1 + v_ij_i * v_ij_k);
  double z_i_m = z_i - .5 * (1. - v_ij_k);
  double z_j_m = z_j - .5 * (1. - v_ij_k);
  double kappa = 2. / 3.;  // must match the I-operator implementation (which uses kappa = 2. / 3.)
  if ((*RA_dipole)[x_a].type_splitting() < 2){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = z_i_m * p_i - z_j_m * p_j;
    if ((*RA_dipole)[x_a].type_splitting() == 0){
      // emitter1 = gluon, emitter2 = gluon
      factor = -1. / (2 * pi_pj + m2_i + m2_j - m2_ij) * 16 * C_A * pi * alpha_S;
      factor_vector = 1. / (v_ij_k * pi_pj);
      factor_metric = 1. / (1. - z_i * one_minus_y_ij_k) + 1. / (1. - z_j * one_minus_y_ij_k) - (2. - kappa * z_plus * z_minus) / v_ij_k;
    }
    else if ((*RA_dipole)[x_a].type_splitting() == 1){
      // emitter1 = quark, emitter2 = antiquark
      factor = -1. / (2 * pi_pj + m2_i + m2_j - m2_ij) * 8 * T_R * pi * alpha_S / v_ij_k;
      factor_vector = -4 / (2 * pi_pj + m2_i + m2_j);
      if (m2_i != m2_j){cout << "m2_i = " << m2_i << " =/= " << m2_j << " = m2_j" << endl; exit(1);}
      factor_metric = 1. - 2 * kappa * (z_plus * z_minus - m2_i / (2 * pi_pj + m2_i + m2_j));
    }

    calculate_dipole_Asc_QCD(x_a, Vtensor, ME2_metric, ME2_vector);

    logger << LOG_DEBUG_VERBOSE << "dipole[" << x_a << "].symmetry_factor() = " << (*RA_dipole)[x_a].symmetry_factor() << endl;
    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) / C_A * (*RA_dipole)[x_a].symmetry_factor();
  }
  else {
    // always (*RA_dipole)[x_a].type_splitting() == 2; emitter1/2 = gluon/(anti)quark
    double Dfactor = 0.;

    calculate_dipole_Acc_QCD(x_a, Dfactor);

    logger << LOG_DEBUG_VERBOSE << "dipole[" << x_a << "].symmetry_factor() = " << (*RA_dipole)[x_a].symmetry_factor() << endl;
    Dfactor = Dfactor / C_F * (*RA_dipole)[x_a].symmetry_factor();

    if (csi->type_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()] == 0){
      // emitter1 = gluon, emitter2 = (anti-)quark
      ME2 = -1. / (2 * pi_pj + m2_i + m2_j - m2_ij) * (8 * C_F * pi * alpha_S) * (2. / (1. - z_j * one_minus_y_ij_k) - (vt_ij_k / v_ij_k) * (1. + z_j + m2_ij / pi_pj)) * Dfactor;
    }
    else {
      // emitter1 = (anti-)quark, emitter2 = gluon
      ME2 = -1. / (2 * pi_pj + m2_i + m2_j - m2_ij) * (8 * C_F * pi * alpha_S) * (2. / (1. - z_i * one_minus_y_ij_k) - (vt_ij_k / v_ij_k) * (1. + z_i + m2_ij / pi_pj)) * Dfactor;
    }
  }
  return -ME2;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
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

#define p_i p_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define p_j p_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define p_a p_parton[0][(*RA_dipole)[x_a].no_R_spectator()]
#define pt_ij p_parton[x_a][(*RA_dipole)[x_a].no_A_emitter()]
#define pt_a p_parton[x_a][(*RA_dipole)[x_a].no_A_spectator()]
double observable_set::calculate_dipole_QCD_A_ij_a(int x_a){
  static Logger logger("observable_set::calculate_dipole_A_ij_a");
  logger << LOG_DEBUG_VERBOSE << "called with dipole no. " << x_a << endl;
  double ME2 = 0.;
  //  int dipole_no_A_emitter = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_emitter()];
  //  int dipole_no_A_spectator = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_spectator()];
  double pi_pj = p_i * p_j;
  double pi_pa = p_i * p_a;
  double pj_pa = p_j * p_a;
  double one_minus_x_ij_a = pi_pj / (pi_pa + pj_pa);
  double x_ij_a = 1. - one_minus_x_ij_a;
  double z_i = pi_pa / (pi_pa + pj_pa);
  double z_j = pj_pa / (pi_pa + pj_pa);
  if ((*RA_dipole)[x_a].type_splitting() < 2){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = z_i * p_i - z_j * p_j;
    if ((*RA_dipole)[x_a].type_splitting() == 0){
      // emitter1 = gluon, emitter2 = gluon
      factor = -1. / (2 * pi_pj * x_ij_a) * 16 * C_A * pi * alpha_S;
      factor_vector = 1. / (pi_pj);
      factor_metric = 1. / (z_j + one_minus_x_ij_a) + 1. / (z_i + one_minus_x_ij_a) - 2.;
    }
    else if ((*RA_dipole)[x_a].type_splitting() == 1){
      // emitter1 = quark, emitter2 = antiquark
      factor = -1. / (2 * pi_pj * x_ij_a) * 8 * T_R * pi * alpha_S;
      factor_vector = -2. / pi_pj;
      factor_metric = 1.;
    }

    calculate_dipole_Asc_QCD(x_a, Vtensor, ME2_metric, ME2_vector);

    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) / C_A * (*RA_dipole)[x_a].symmetry_factor();
  }
  else {
    // always (*RA_dipole)[x_a].type_splitting() == 2; emitter1/2 = gluon/(anti)quark
    double Dfactor = 0.;

    calculate_dipole_Acc_QCD(x_a, Dfactor);

    Dfactor = Dfactor / C_F * (*RA_dipole)[x_a].symmetry_factor();

    ME2 = -1. / (2 * x_ij_a * pi_pj) * (8 * C_F * pi * alpha_S) * (2. / (z_j + one_minus_x_ij_a) - (1. + z_i)) * Dfactor;
  }
  return -ME2;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

#define m2_i mass2_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define m2_j mass2_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define m2_ij mass2_parton[x_a][(*RA_dipole)[x_a].no_A_emitter()]
double observable_set::calculate_dipole_QCD_A_ij_a_massive(int x_a){
  static Logger logger("observable_set::calculate_dipole_QCD_A_ij_a_massive");
  logger << LOG_DEBUG_VERBOSE << "called with dipole no. " << x_a << endl;
  double ME2 = 0.;
  //  int dipole_no_A_emitter = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_emitter()];
  //  int dipole_no_A_spectator = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_spectator()];
  double pi_pj = p_i * p_j;
  double pi_pa = p_i * p_a;
  double pj_pa = p_j * p_a;
  double one_minus_x_ij_a = (pi_pj - .5 * (m2_ij - m2_i - m2_j)) / (pi_pa + pj_pa);
  double x_ij_a = 1. - one_minus_x_ij_a;
  double z_i = pi_pa / (pi_pa + pj_pa);
  double z_j = pj_pa / (pi_pa + pj_pa);
  if ((*RA_dipole)[x_a].type_splitting() < 2){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = z_i * p_i - z_j * p_j;
    if ((*RA_dipole)[x_a].type_splitting() == 0){
      // emitter1 = gluon, emitter2 = gluon
      cout << "calculate_QCD_A_ij_a_massive: Should not happen! g -> gg splitting with initial-state spectator" << endl;
      //      exit(1);
    }
    else if ((*RA_dipole)[x_a].type_splitting() == 1){
      // emitter1 = quark, emitter2 = antiquark
      factor = -1. / ((2 * pi_pj + m2_i + m2_j - m2_ij) * x_ij_a) * 8 * T_R * pi * alpha_S;
      factor_vector = -4. / (2 * pi_pj + m2_i + m2_j);
      factor_metric = 1.;
    }

    calculate_dipole_Asc_QCD(x_a, Vtensor, ME2_metric, ME2_vector);

    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) / C_A * (*RA_dipole)[x_a].symmetry_factor();
  }
  else {
    // always (*RA_dipole)[x_a].type_splitting() == 2; emitter1/2 = gluon/(anti)quark
    double Dfactor = 0.;

    calculate_dipole_Acc_QCD(x_a, Dfactor);

    Dfactor = Dfactor / C_F * (*RA_dipole)[x_a].symmetry_factor();
    if (csi->type_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()] == 0){
      cout << "ij a: emitter1 = gluon, emitter2 = (anti-)quark" << endl;
      // emitter1 = gluon, emitter2 = (anti-)quark
      if (m2_j != m2_ij){exit(1);}
      ME2 = -1. / ((2 * pi_pj + m2_i + m2_j - m2_ij) * x_ij_a) * (8 * C_F * pi * alpha_S) * (2. / (one_minus_x_ij_a + z_i) - 1. - z_j - m2_ij / pi_pj) * Dfactor;
    }
    else {
      // emitter1 = (anti-)quark, emitter2 = gluon
      if (m2_i != m2_ij){exit(1);}
      ME2 = -1. / ((2 * pi_pj + m2_i + m2_j - m2_ij) * x_ij_a) * (8 * C_F * pi * alpha_S) * (2. / (one_minus_x_ij_a + z_j) - 1. - z_i - m2_ij / pi_pj) * Dfactor;
    }
  }
  return -ME2;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
#undef m2_i
#undef m2_j
#undef m2_ij

#undef p_i
#undef p_j
#undef p_a
#undef pt_ij
#undef pt_a

#define p_a p_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define p_i p_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define p_k p_parton[0][(*RA_dipole)[x_a].no_R_spectator()]
double observable_set::calculate_dipole_QCD_A_ai_k(int x_a){
  static Logger logger("observable_set::calculate_dipole_A_ai_k");
  logger << LOG_DEBUG_VERBOSE << "called with dipole no. " << x_a << endl;
  double ME2 = 0.;
  //  int dipole_no_A_emitter = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_emitter()];
  //  int dipole_no_A_spectator = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_spectator()];
  double pi_pk = p_i * p_k;
  double pi_pa = p_i * p_a;
  double pk_pa = p_k * p_a;
  double one_minus_x_ik_a = pi_pk / (pi_pa + pk_pa);
  double x_ik_a = 1. - one_minus_x_ik_a;
  double u_i = pi_pa / (pi_pa + pk_pa);
  double one_minus_u_i = 1. - u_i;
  if ((*RA_dipole)[x_a].type_splitting() < 2){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = p_i / u_i - p_k / one_minus_u_i;
    if ((*RA_dipole)[x_a].type_splitting() == 0){
      // emitter1 = gluon, emitter2 = gluon -> emitter12 = gluon
      factor = -1. / (2 * pi_pa * x_ik_a) * 16 * C_A * pi * alpha_S;
      factor_vector = (one_minus_x_ik_a * u_i * one_minus_u_i) / (x_ik_a * pi_pk);
      factor_metric = 1. / (one_minus_x_ik_a + u_i) - 1. + x_ik_a * one_minus_x_ik_a;
    }
    else if ((*RA_dipole)[x_a].type_splitting() == 1){
      // emitter1 = (anti)quark, emitter2 = (anti)quark -> emitter12 = gluon
      factor = -1. / (2 * pi_pa * x_ik_a) * 8 * C_F * pi * alpha_S;
      factor_vector = (2 * one_minus_x_ik_a * u_i * one_minus_u_i) / (x_ik_a * pi_pk);
      factor_metric = x_ik_a;
    }

    calculate_dipole_Asc_QCD(x_a, Vtensor, ME2_metric, ME2_vector);

    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) / C_A * (*RA_dipole)[x_a].symmetry_factor();
  }
  else {
    double Dfactor = 0.;

    calculate_dipole_Acc_QCD(x_a, Dfactor);

    Dfactor = Dfactor / C_F * (*RA_dipole)[x_a].symmetry_factor();
    if      ((*RA_dipole)[x_a].type_splitting() == 2){
      // emitter1 = (anti)quark, emitter2 = gluon -> emitter12 = (anti)quark
      ME2 = -1. / (2 * x_ik_a * pi_pa) * (8 * C_F * pi * alpha_S) * (2. / (one_minus_x_ik_a + u_i) - (1. + x_ik_a)) * Dfactor;
    }
    else if ((*RA_dipole)[x_a].type_splitting() == 3){
      // emitter1 = gluon, emitter2 = (anti)quark -> emitter12 = (anti)quark
      ME2 = -1. / (2 * x_ik_a * pi_pa) * (8 * T_R * pi * alpha_S) * (1. - 2 * x_ik_a * one_minus_x_ik_a) * Dfactor;
    }
  }
  return -ME2;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
#undef p_a
#undef p_i
#undef p_k

#define p_a p_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define p_i p_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define p_b p_parton[0][(*RA_dipole)[x_a].no_R_spectator()]
double observable_set::calculate_dipole_QCD_A_ai_b(int x_a){
  static Logger logger("observable_set::calculate_dipole_A_ai_b");
  logger << LOG_DEBUG_VERBOSE << "called with dipole no. " << x_a << endl;

  /*
  cout << "(*RA_dipole)[" << x_a << "].type_splitting() = " << (*RA_dipole)[x_a].type_splitting() << endl;
  cout << "(*RA_dipole)[" << x_a << "].no_A_emitter() = " << (*RA_dipole)[x_a].no_A_emitter() << endl;
  cout << "(*RA_dipole)[" << x_a << "].no_A_spectator() = " << (*RA_dipole)[x_a].no_A_spectator() << endl;
  cout << "(*RA_dipole)[" << x_a << "].o_prc().size() = " << (*RA_dipole)[x_a].o_prc().size() << endl;
  cout << "(*RA_dipole)[" << x_a << "].o_prc()[(*RA_dipole)[" << x_a << "].no_A_emitter()] = " << (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_emitter()] << endl;
  cout << "(*RA_dipole)[" << x_a << "].o_prc()[(*RA_dipole)[" << x_a << "].no_A_spectator()] = " << (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_spectator()] << endl;
  */

  double ME2 = 0.;
  //  int dipole_no_A_emitter = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_emitter()];
  //  int dipole_no_A_spectator = (*RA_dipole)[x_a].o_prc()[(*RA_dipole)[x_a].no_A_spectator()];
  double x_i_ab = 1. - (p_i * (p_a + p_b)) / (p_a * p_b);
  // define some z_i analogously !!!
  if ((*RA_dipole)[x_a].type_splitting() < 2){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = p_i - p_b * (p_i * p_a) / (p_b * p_a);
    if ((*RA_dipole)[x_a].type_splitting() == 0){
      // emitter1 = gluon, emitter2 = gluon -> emitter12 = gluon
      factor = -1. / (2. * x_i_ab * (p_a * p_i)) * 16. * C_A * pi * alpha_S;
      factor_vector = ((1. - x_i_ab) * (p_a * p_b)) / (x_i_ab * (p_i * p_b) * (p_i * p_a));
      factor_metric = x_i_ab / (1 - x_i_ab) + x_i_ab * (1. - x_i_ab);
    }
    else if ((*RA_dipole)[x_a].type_splitting() == 1){
      // emitter1 = (anti)quark, emitter2 = (anti)quark -> emitter12 = gluon
      factor = -1. / (2. * x_i_ab * (p_a * p_i)) * 8. * C_F * pi * alpha_S;
      factor_vector = (2. * (1. - x_i_ab) * (p_a * p_b)) / (x_i_ab * (p_a * p_i) * (p_i * p_b));
      factor_metric = x_i_ab;
    }

    calculate_dipole_Asc_QCD(x_a, Vtensor, ME2_metric, ME2_vector);

    /*
    cout << "factor_metric = " << factor_metric << endl;
    cout << "ME2_metric = " << ME2_metric << endl;
    cout << "factor_vector = " << factor_vector << endl;
    cout << "Vtensor.m2() = " << Vtensor.m2() << endl;
    cout << "ME2_vector = " << ME2_vector << endl;
    */
    /*
    // guessed factor 2 !!!
    if ((*RA_dipole)[x_a].type_splitting() == 0){
      ME2 = factor * (factor_metric * ME2_metric - 2 * Vtensor.m2() * factor_vector * ME2_vector) / C_A * (*RA_dipole)[x_a].symmetry_factor();
    }
    else if ((*RA_dipole)[x_a].type_splitting() == 1){
      ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) / C_A * (*RA_dipole)[x_a].symmetry_factor();
    }
    */
    //    cout << "osi_i_acc = " << setw(12) << "" << "   metric = " << factor * (factor_metric * ME2_metric) / C_A * (*RA_dipole)[x_a].symmetry_factor() << endl;
    //    cout << "            " << setw(12) << "" << "   vector = " << factor * (- Vtensor.m2() * factor_vector * ME2_vector) / C_A * (*RA_dipole)[x_a].symmetry_factor() << endl;



    // correct version - reactivate !!!
    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) / C_A * (*RA_dipole)[x_a].symmetry_factor();
  }
  else {
    double Dfactor = 0.;

    calculate_dipole_Acc_QCD(x_a, Dfactor);

    Dfactor = Dfactor / C_F * (*RA_dipole)[x_a].symmetry_factor();
    if      ((*RA_dipole)[x_a].type_splitting() == 2){
      // emitter1 = (anti)quark, emitter2 = gluon -> emitter12 = (anti)quark
      ME2 = -1. / (2. * x_i_ab * (p_a * p_i)) * (8. * C_F * pi * alpha_S) * (2. / (1. - x_i_ab) - (1. + x_i_ab)) * Dfactor;
    }
    else if ((*RA_dipole)[x_a].type_splitting() == 3){
      // emitter1 = gluon, emitter2 = (anti)quark -> emitter12 = (anti)quark
      ME2 = -1. / (2. * x_i_ab * (p_a * p_i)) * (8. * T_R * pi * alpha_S) * (1. - 2. * x_i_ab * (1. - x_i_ab)) * Dfactor;
    }
  }
  return -ME2;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
#undef p_a
#undef p_i
#undef p_b






void observable_set::calculate_dipole_Acc_QEW(int x_a, double & Dfactor){
  static Logger logger("observable_set::calculate_dipole_Acc_QEW");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (switch_OL == 1){

  static int n_momentum = 5 * (p_parton[x_a].size() - 1);
  static double *P;
  P = new double[n_momentum];
  for (int i = 1; i < p_parton[x_a].size(); i++){
    P[5 * (i - 1)]     = p_parton[x_a][i].x0();
    P[5 * (i - 1) + 1] = p_parton[x_a][i].x1();
    P[5 * (i - 1) + 2] = p_parton[x_a][i].x2();
    P[5 * (i - 1) + 3] = p_parton[x_a][i].x3();
    P[5 * (i - 1) + 4] = 0.;
    //    P[5 * (i - 1) + 4] = p_parton[x_a][i].m();
  }
  ol_evaluate_tree((*RA_dipole)[x_a].process_id, P, &Dfactor);
  logger << LOG_DEBUG_POINT << "OpenLoops:     Dfactor[" << x_a << "] = " << setw(23) << setprecision(15) << Dfactor << endl;
  delete [] P;

  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
void observable_set::calculate_dipole_Asc_QEW(int x_a, fourvector & Vtensor, double & ME2_metric, double & ME2_vector){
  static Logger logger("observable_set::calculate_dipole_Asc_QEW");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (switch_OL == 1){

  static int n_momentum = 5 * (p_parton[x_a].size() - 1);
  static double *P;
  P = new double[n_momentum];
  for (int i = 1; i < p_parton[x_a].size(); i++){
    P[5 * (i - 1)]     = p_parton[x_a][i].x0();
    P[5 * (i - 1) + 1] = p_parton[x_a][i].x1();
    P[5 * (i - 1) + 2] = p_parton[x_a][i].x2();
    P[5 * (i - 1) + 3] = p_parton[x_a][i].x3();
    P[5 * (i - 1) + 4] = p_parton[x_a][i].m();
  }
  double *MOM;
  MOM = new double[4];
  MOM[0] = Vtensor.x0();
  MOM[1] = Vtensor.x1();
  MOM[2] = Vtensor.x2();
  MOM[3] = Vtensor.x3();
  ol_evaluate_tree((*RA_dipole)[x_a].process_id, P, &ME2_metric);
  logger << LOG_DEBUG_POINT << "OpenLoops:     ME2_metric[" << x_a << "] = " << setw(23) << setprecision(15) << ME2_metric << endl;
  static int n_sc = p_parton[x_a].size() - 1;
  static double *M2sc;
  M2sc = new double[n_sc];
  ol_evaluate_sc((*RA_dipole)[x_a].process_id, P, (*RA_dipole)[x_a].no_A_emitter(), MOM, M2sc);
  ME2_vector = M2sc[0];
  logger << LOG_DEBUG_POINT << "OpenLoops:     ME2_vector[" << x_a << "] = " << setw(23) << setprecision(15) << ME2_vector << endl;
  /*
  logger << LOG_DEBUG << "Vtensor = " << Vtensor << "   M2sc[" << (*RA_dipole)[x_a].no_A_spectator() - 1 << "] = " << ME2_vector << endl;
  for (int i = 0; i < n_sc; i++){logger << LOG_DEBUG << "M2sc[" << i << "] = " << M2sc[i] << endl;}
  */
  //  ME2_vector = M2sc[(*RA_dipole)[x_a].no_A_spectator() - 1];
  // all M2sc empty !!!
  // check if correct colour-correlaton (none) is used !!!
  delete [] M2sc;
  delete [] MOM;
  delete [] P;

  }
}

#define p_i p_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define p_j p_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define p_k p_parton[0][(*RA_dipole)[x_a].no_R_spectator()]
#define pt_ij p_parton[x_a][(*RA_dipole)[x_a].no_A_emitter()]
#define pt_k p_parton[x_a][(*RA_dipole)[x_a].no_A_spectator()]
double observable_set::calculate_dipole_QEW_A_ij_k(int x_a){
  static Logger logger("observable_set::calculate_dipole_QEW_A_ij_k");
  logger << LOG_DEBUG_VERBOSE << "called    with dipole no. " << x_a << endl;
  double ME2 = 0.;
  double pi_pj = p_i * p_j;
  double pj_pk = p_j * p_k;
  double pi_pk = p_i * p_k;
  double y_ij_k = pi_pj / (pi_pj + pi_pk + pj_pk);
  double one_minus_y_ij_k = 1. - y_ij_k;
  double z_i = pi_pk / (pi_pk + pj_pk);
  double z_j = pj_pk / (pi_pk + pj_pk);
  if ((*RA_dipole)[x_a].type_splitting() == 1){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = z_i * p_i - z_j * p_j;
    // emitter1 = charged (anti)particle, emitter2 = charged (anti)particle
    factor = -1. / (2 * pi_pj) * 8 * N_c * pi * msi.alpha_e;
    //  T_R / C_A -> N_c ???
    // adapt to paper version... Q_f^2 etc. !!!
    // formula seems ok.
    factor_vector = -2. / pi_pj;
    factor_metric = 1.;

    calculate_dipole_Asc_QEW(x_a, Vtensor, ME2_metric, ME2_vector);

    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) * (*RA_dipole)[x_a].symmetry_factor() * (*RA_dipole)[x_a].charge_factor();
  }
  else {
    // always (*RA_dipole)[x_a].type_splitting() == 2; emitter1 = charged (anti)particle, emitter2 = photon
    double Dfactor;

    calculate_dipole_Acc_QEW(x_a, Dfactor);

    ME2 = -1. / (2. * pi_pj) * (8 * pi * msi.alpha_e) * (2. / (1. - z_i * one_minus_y_ij_k) - (1. + z_i)) * Dfactor * (*RA_dipole)[x_a].symmetry_factor() * (*RA_dipole)[x_a].charge_factor();
    // C_F -> Q_f^2, cancelled against T_ij^2 -> Q_f^2 in the denominator
  }
  return -ME2;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
#define m2_i mass2_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define m2_j mass2_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define m2_k mass2_parton[0][(*RA_dipole)[x_a].no_R_spectator()]
#define m2_ij mass2_parton[x_a][(*RA_dipole)[x_a].no_A_emitter()]
double observable_set::calculate_dipole_QEW_A_ij_k_massive(int x_a){
  static Logger logger("observable_set::calculate_dipole_QEW_A_ij_k_massive");
  logger << LOG_DEBUG_VERBOSE << "called    with dipole no. " << x_a << endl;
  double ME2 = 0.;
  double pi_pj = p_i * p_j;
  double pj_pk = p_j * p_k;
  double pi_pk = p_i * p_k;
  double y_ij_k = pi_pj / (pi_pj + pi_pk + pj_pk);
  double one_minus_y_ij_k = 1. - y_ij_k;
  double z_i = pi_pk / (pi_pk + pj_pk);
  double z_j = pj_pk / (pi_pk + pj_pk);
  fourvector Q = pt_ij + pt_k;
  double Q2 = Q.m2();
  double mu2_i = m2_i / Q2;
  double mu2_j = m2_j / Q2;
  double mu2_k = m2_k / Q2;
  double mu2_ij = m2_ij / Q2;
  double one_minus_all_mu = 1. - mu2_i - mu2_j - mu2_k;
  double z_prefactor = (2 * mu2_i + one_minus_all_mu * y_ij_k) / (2 * (mu2_i + mu2_j + one_minus_all_mu * y_ij_k));
  double v_ij_k = sqrt(pow(2 * mu2_k + one_minus_all_mu * one_minus_y_ij_k, 2) - 4 * mu2_k) / (one_minus_all_mu * one_minus_y_ij_k);
  double v_ij_i = sqrt(pow(one_minus_all_mu * y_ij_k, 2) - 4 * mu2_i * mu2_j) / (one_minus_all_mu * y_ij_k + 2 * mu2_i);
  double vt_ij_k = sqrt(lambda(1, mu2_ij, mu2_k)) / (1. - mu2_ij - mu2_k);
  double z_minus = z_prefactor * (1 - v_ij_i * v_ij_k);
  double z_plus = z_prefactor * (1 + v_ij_i * v_ij_k);
  double z_i_m = z_i - .5 * (1. - v_ij_k);
  double z_j_m = z_j - .5 * (1. - v_ij_k);
  double kappa = 2. / 3.;  // must match the I-operator implementation (which uses kappa = 2. / 3.)
  if ((*RA_dipole)[x_a].type_splitting() == 1){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = z_i_m * p_i - z_j_m * p_j;
    // emitter1 = charged particle, emitter2 = charged antiparticle
    factor = -1. / (2 * pi_pj + m2_i + m2_j - m2_ij) * 8 * T_R * pi * msi.alpha_e / v_ij_k;
    factor_vector = -4 / (2 * pi_pj + m2_i + m2_j);
    if (m2_i != m2_j){cout << "m2_i = " << m2_i << " =/= " << m2_j << " = m2_j" << endl; exit(1);}
    factor_metric = 1. - 2 * kappa * (z_plus * z_minus - m2_i / (2 * pi_pj + m2_i + m2_j));

    calculate_dipole_Asc_QEW(x_a, Vtensor, ME2_metric, ME2_vector);

    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) / C_A * (*RA_dipole)[x_a].symmetry_factor();
  }
  else {
    // always (*RA_dipole)[x_a].type_splitting() == 2; emitter1/2 = charged (anti-)particle/photon
    double Dfactor = 0.;

    calculate_dipole_Acc_QEW(x_a, Dfactor);

    Dfactor = Dfactor / C_F * (*RA_dipole)[x_a].symmetry_factor();

    // emitter1 = charged (anti-)particle, emitter2 = photon
    ME2 = -1. / (2 * pi_pj + m2_i + m2_j - m2_ij) * (8 * C_F * pi * msi.alpha_e) * (2. / (1. - z_i * one_minus_y_ij_k) - (vt_ij_k / v_ij_k) * (1. + z_i + m2_ij / pi_pj)) * Dfactor;
    // C_F -> Q_f^2, cancelled against T_ij^2 -> Q_f^2 in the denominator
    // has not yet been replaced here !!!
  }
  return -ME2 * (*RA_dipole)[x_a].charge_factor();
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
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

#define p_i p_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define p_j p_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define p_a p_parton[0][(*RA_dipole)[x_a].no_R_spectator()]
#define pt_ij p_parton[x_a][(*RA_dipole)[x_a].no_A_emitter()]
#define pt_a p_parton[x_a][(*RA_dipole)[x_a].no_A_spectator()]
double observable_set::calculate_dipole_QEW_A_ij_a(int x_a){
  static Logger logger("observable_set::calculate_dipole_QEW_A_ij_a");
  logger << LOG_DEBUG_VERBOSE << "called    with dipole no. " << x_a << endl;
  double ME2 = 0.;
  double pi_pj = p_i * p_j;
  double pi_pa = p_i * p_a;
  double pj_pa = p_j * p_a;
  double one_minus_x_ij_a = pi_pj / (pi_pa + pj_pa);
  double x_ij_a = 1. - one_minus_x_ij_a;
  double z_i = pi_pa / (pi_pa + pj_pa);
  double z_j = pj_pa / (pi_pa + pj_pa);
  if ((*RA_dipole)[x_a].type_splitting() == 1){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = z_i * p_i - z_j * p_j;
    // emitter1 = charged particle, emitter2 = charged antiparticle
    factor = -1. / (2 * pi_pj * x_ij_a) * 8 * N_c * pi * msi.alpha_e;
    //  T_R / C_A -> N_c ???
    // adapt to paper version... Q_f^2 etc. !!!
    factor_vector = -2. / pi_pj;
    factor_metric = 1.;

    calculate_dipole_Asc_QEW(x_a, Vtensor, ME2_metric, ME2_vector);

    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) * (*RA_dipole)[x_a].symmetry_factor() * (*RA_dipole)[x_a].charge_factor();
  }
  else {
    // always (*RA_dipole)[x_a].type_splitting() == 2; emitter1/2 = charged (anti)particle/photon
    double Dfactor;

    calculate_dipole_Acc_QEW(x_a, Dfactor);

    ME2 = -1. / (2 * x_ij_a * pi_pj) * (8 * pi * msi.alpha_e) * (2. / (z_j + one_minus_x_ij_a) - (1. + z_i)) * Dfactor * (*RA_dipole)[x_a].symmetry_factor() * (*RA_dipole)[x_a].charge_factor();
    // C_F -> Q_f^2, cancelled against Q_f^2 in the denominator
  }
  return -ME2;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
#define m2_i mass2_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define m2_j mass2_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define m2_ij mass2_parton[x_a][(*RA_dipole)[x_a].no_A_emitter()]
double observable_set::calculate_dipole_QEW_A_ij_a_massive(int x_a){
  static Logger logger("observable_set::calculate_dipole_QEW_A_ij_a_massive");
  logger << LOG_DEBUG_VERBOSE << "A_ij_a_massive   called    with dipole no. " << x_a << endl;
  double ME2 = 0.;
  double pi_pj = p_i * p_j;
  double pi_pa = p_i * p_a;
  double pj_pa = p_j * p_a;
  double one_minus_x_ij_a = (pi_pj - .5 * (m2_ij - m2_i - m2_j)) / (pi_pa + pj_pa);
  double x_ij_a = 1. - one_minus_x_ij_a;
  double z_i = pi_pa / (pi_pa + pj_pa);
  double z_j = pj_pa / (pi_pa + pj_pa);
  if ((*RA_dipole)[x_a].type_splitting() == 1){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = z_i * p_i - z_j * p_j;
    // emitter1 = quark, emitter2 = antiquark
    factor = -1. / ((2 * pi_pj + m2_i + m2_j - m2_ij) * x_ij_a) * 8 * T_R * pi * msi.alpha_e;
    factor_vector = -4. / (2 * pi_pj + m2_i + m2_j);
    factor_metric = 1.;

    calculate_dipole_Asc_QEW(x_a, Vtensor, ME2_metric, ME2_vector);

    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) / C_A * (*RA_dipole)[x_a].symmetry_factor();
  }
  else {
    // always (*RA_dipole)[x_a].type_splitting() == 2; emitter1/2 = charged (anti)particle/photon
    double Dfactor = 0.;

    calculate_dipole_Acc_QEW(x_a, Dfactor);

    Dfactor = Dfactor / C_F * (*RA_dipole)[x_a].symmetry_factor();
    // emitter1 = charged (anti)particle, emitter2 = photon
    if (m2_i != m2_ij){exit(1);}
    ME2 = -1. / ((2 * pi_pj + m2_i + m2_j - m2_ij) * x_ij_a) * (8 * C_F * pi * msi.alpha_e) * (2. / (one_minus_x_ij_a + z_j) - 1. - z_i - m2_ij / pi_pj) * Dfactor;
    // C_F -> Q_f^2, cancelled against T_ij^2 -> Q_f^2 in the denominator
    // has not yet been replaced here !!!
  }
  return -ME2 * (*RA_dipole)[x_a].charge_factor();
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
#undef m2_i
#undef m2_j
#undef m2_ij
#undef p_i
#undef p_j
#undef p_a
#undef pt_ij
#undef pt_a

#define p_a p_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define p_i p_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define p_k p_parton[0][(*RA_dipole)[x_a].no_R_spectator()]
double observable_set::calculate_dipole_QEW_A_ai_k(int x_a){
  static Logger logger("observable_set::calculate_dipole_QEW_A_ai_k");
  logger << LOG_DEBUG_VERBOSE << "called    with dipole no. " << x_a << endl;
  double ME2 = 0.;
  double pi_pk = p_i * p_k;
  double pi_pa = p_i * p_a;
  double pk_pa = p_k * p_a;
  double one_minus_x_ik_a = pi_pk / (pi_pa + pk_pa);
  double x_ik_a = 1. - one_minus_x_ik_a;
  double u_i = pi_pa / (pi_pa + pk_pa);
  double one_minus_u_i = 1. - u_i;
  if ((*RA_dipole)[x_a].type_splitting() == 1){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = p_i / u_i - p_k / one_minus_u_i;
    // emitter1 = charged (anti)particle, emitter2 = charged (anti)particle -> emitter12 = photon
    factor = -1. / (2 * pi_pa * x_ik_a) * 8 * pi * msi.alpha_e;
    //    factor = -1. / (2 * pi_pa * x_ik_a) * 8 * 0.5 * pi * msi.alpha_e;
    // C_F / C_A -> .5
    // adapt to paper !!! not used so far !!!
    factor_vector = (2 * one_minus_x_ik_a * u_i * one_minus_u_i) / (x_ik_a * pi_pk);
    factor_metric = x_ik_a;

    calculate_dipole_Asc_QEW(x_a, Vtensor, ME2_metric, ME2_vector);

    //    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) * (*RA_dipole)[x_a].symmetry_factor() * (*RA_dipole)[x_a].charge_factor();
    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) * (*RA_dipole)[x_a].symmetry_factor() * (*RA_dipole)[x_a].charge_factor() * 0.;
  }
  else {
    double Dfactor = 0.;

    calculate_dipole_Acc_QEW(x_a, Dfactor);

    if      ((*RA_dipole)[x_a].type_splitting() == 2){
      // emitter1 = charged (anti)particle, emitter2 = photon -> emitter12 = charged (anti)particle
      // adapt to paper !!! not used so far !!!
      ME2 = -1. / (2 * x_ik_a * pi_pa) * (8 * pi * msi.alpha_e) * (2. / (one_minus_x_ik_a + u_i) - (1. + x_ik_a)) * Dfactor * (*RA_dipole)[x_a].symmetry_factor() * (*RA_dipole)[x_a].charge_factor();
    }
    else if ((*RA_dipole)[x_a].type_splitting() == 3){
      // emitter1 = photon, emitter2 = charged (anti)particle -> emitter12 = charged (anti)particle
      ME2 = -1. / (2 * x_ik_a * pi_pa) * (8 * N_c * pi * msi.alpha_e) * (1. - 2 * x_ik_a * one_minus_x_ik_a) * Dfactor * (*RA_dipole)[x_a].symmetry_factor() * (*RA_dipole)[x_a].charge_factor();
      // T_R / C_F -> N_c
    // adapt to paper !!! not used so far !!!
    }
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return -ME2;
}
#undef p_a
#undef p_i
#undef p_k

#define p_a p_parton[0][(*RA_dipole)[x_a].no_R_emitter_1()]
#define p_i p_parton[0][(*RA_dipole)[x_a].no_R_emitter_2()]
#define p_b p_parton[0][(*RA_dipole)[x_a].no_R_spectator()]
double observable_set::calculate_dipole_QEW_A_ai_b(int x_a){
  static Logger logger("observable_set::calculate_dipole_QEW_A_ai_b");
  logger << LOG_DEBUG_VERBOSE << "called    with dipole no. " << x_a << endl;
  double ME2 = 0.;
  double x_i_ab = 1. - (p_i * (p_a + p_b)) / (p_a * p_b);
  // define some z_i analogously !!!
  if ((*RA_dipole)[x_a].type_splitting() == 1){
    double factor = 0.;
    double ME2_metric = 0.;
    double ME2_vector = 0.;
    double factor_metric = 0.;
    double factor_vector = 0.;
    fourvector Vtensor = p_i - p_b * (p_i * p_a) / (p_b * p_a);
    // emitter1 = charged (anti)particle, emitter2 = charged (anti)particle -> emitter12 = photon
    factor = -1. / (2. * x_i_ab * (p_a * p_i)) * 8. * pi * msi.alpha_e;
    //    factor = -1. / (2. * x_i_ab * (p_a * p_i)) * 8. * 0.5 * pi * msi.alpha_e;
    // C_F / C_A -> .5
    // adapt to paper !!! not used so far !!!
    factor_vector = (2. * (1. - x_i_ab) * (p_a * p_b)) / (x_i_ab * (p_a * p_i) * (p_i * p_b)); // manually introduced factor 0.25!!! -> should be fixed now !!!
    factor_metric = x_i_ab;

    calculate_dipole_Asc_QEW(x_a, Vtensor, ME2_metric, ME2_vector);

    ME2 = factor * (factor_metric * ME2_metric - Vtensor.m2() * factor_vector * ME2_vector) * (*RA_dipole)[x_a].symmetry_factor() * (*RA_dipole)[x_a].charge_factor();
  }
  else {
    double Dfactor = 0.;

    calculate_dipole_Acc_QEW(x_a, Dfactor);

    if      ((*RA_dipole)[x_a].type_splitting() == 2){
      // emitter1 = charged (anti)particle, emitter2 = photon -> emitter12 = charged (anti)particle
      ME2 = -1. / (2 * x_i_ab * (p_a * p_i)) * (8 * pi * msi.alpha_e) * (2. / (1. - x_i_ab) - (1. + x_i_ab)) * Dfactor * (*RA_dipole)[x_a].symmetry_factor() * (*RA_dipole)[x_a].charge_factor();
    // adapt to paper !!! not used so far !!!
    }
    else if ((*RA_dipole)[x_a].type_splitting() == 3){
      // emitter1 = photon, emitter2 = charged (anti)particle -> emitter12 = charged (anti)particle
      ME2 = -1. / (2 * x_i_ab * (p_a * p_i)) * (8 * N_c * pi * msi.alpha_e) * (1. - 2 * x_i_ab * (1. - x_i_ab)) * Dfactor * (*RA_dipole)[x_a].symmetry_factor() * (*RA_dipole)[x_a].charge_factor();
      // T_R / C_F -> N_c
    // adapt to paper !!! not used so far !!!
    }
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return -ME2;
}
#undef p_a
#undef p_i
#undef p_b
