#include "../include/classes.cxx"

#define p_i psi_xbp_all[0][bp_i]
#define p_j psi_xbp_all[0][bp_j]
#define p_k psi_xbp_all[0][bp_k]
#define pt_ij psi_xbp_all[i_a][bpt_ij]
#define pt_k psi_xbp_all[i_a][bpt_k]
#define m2_i psi_xbs_all[0][bp_i]
#define m2_j psi_xbs_all[0][bp_j]
#define m2_k psi_xbs_all[0][bp_k]
#define m2_ij psi_xbs_all[i_a][bpt_ij]
#define y_ij_k dipole[i_a].xy
void phasespace_ij_k_massive(vector<vector<fourvector> > & psi_xbp_all, vector<vector<double> > & psi_xbs_all, vector<vector<double> > & psi_xbsqrts_all, int i_a, vector<dipole_set> & dipole){
  static Logger logger("phasespace_ij_k_massive");
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

  fourvector Q = p_i + p_j + p_k;
  double Q2 = Q.m2();
  //  double sqrtQ2 = sqrt(Q2);
  double pi_pj = p_i * p_j;
  pt_k = sqrt(lambda(Q2, m2_ij, m2_k) / lambda(Q2, m2_i + m2_j + 2 * pi_pj, m2_k)) * (p_k - ((Q * p_k) / Q2) * Q) + ((Q2 + m2_k - m2_ij) / (2 * Q2)) * Q;
  pt_ij = Q - pt_k;

  // uncommented !!!
  y_ij_k = (p_i * p_j) / (p_i * p_j + p_j * p_k + p_k * p_i);
  //  pt_ij = p_i + p_j - y_ij_k / (1. - y_ij_k) * p_k;
  //  pt_k = 1. / (1. - y_ij_k) * p_k;
  
  for (int xbi = 1; xbi < psi_xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else if (xbi == bp_j){}
    else if (xbi == bp_k){}
    else if (xbi < bp_j){psi_xbp_all[i_a][xbi] = psi_xbp_all[0][xbi];}
    else if (xbi > bp_j){psi_xbp_all[i_a][xbi / 2] = psi_xbp_all[0][xbi];}
    else {cout << "Should not happen!" << endl;}
  }
  psi_xbp_all[i_a][0] = psi_xbp_all[i_a][1] + psi_xbp_all[i_a][2];
  psi_xbs_all[i_a][0] = psi_xbs_all[0][0];
  psi_xbsqrts_all[i_a][0] = psi_xbsqrts_all[0][0];
  int xbn_all = psi_xbp_all[i_a].size() - 4;
  psi_xbp_all[i_a][xbn_all] = psi_xbp_all[i_a][0];
  psi_xbs_all[i_a][xbn_all] = psi_xbs_all[i_a][0];
  psi_xbsqrts_all[i_a][xbn_all] = psi_xbsqrts_all[i_a][0];

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

#define p_i psi_xbp_all[0][bp_i]
#define p_j psi_xbp_all[0][bp_j]
#define p_a psi_xbp_all[0][bp_a]
#define pt_ij psi_xbp_all[i_a][bpt_ij]
#define pt_a psi_xbp_all[i_a][bpt_a]
#define m2_i psi_xbs_all[0][bp_i]
#define m2_j psi_xbs_all[0][bp_j]
#define m2_ij psi_xbs_all[i_a][bpt_ij]
#define x_ij_a dipole[i_a].xy
void phasespace_ij_a_massive(vector<vector<fourvector> > & psi_xbp_all, vector<vector<double> > & psi_xbs_all, vector<vector<double> > & psi_xbsqrts_all, int i_a, vector<dipole_set> & dipole){
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
  for (int xbi = 1; xbi < psi_xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else if (xbi == bp_j){}
    else if (xbi == bp_a){}
    else if (xbi < bp_j){psi_xbp_all[i_a][xbi] = psi_xbp_all[0][xbi];}
    else if (xbi > bp_j){psi_xbp_all[i_a][xbi / 2] = psi_xbp_all[0][xbi];}
    else {cout << "Should not happen!" << endl;}
  }

  psi_xbp_all[i_a][0] = psi_xbp_all[i_a][1] + psi_xbp_all[i_a][2];
  psi_xbs_all[i_a][0] = x_ij_a * psi_xbs_all[0][0];
  psi_xbsqrts_all[i_a][0] = sqrt(psi_xbs_all[i_a][0]);
  int xbn_all = psi_xbp_all[i_a].size() - 4;
  psi_xbp_all[i_a][xbn_all] = psi_xbp_all[i_a][0];
  psi_xbs_all[i_a][xbn_all] = psi_xbs_all[i_a][0];
  psi_xbsqrts_all[i_a][xbn_all] = psi_xbsqrts_all[i_a][0];

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
