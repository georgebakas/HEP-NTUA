#include "../include/classes.cxx"

#define p_i psi_xbp_all[0][bp_i]
#define p_j psi_xbp_all[0][bp_j]
#define p_k psi_xbp_all[0][bp_k]
#define pt_ij psi_xbp_all[i_a][bpt_ij]
#define pt_k psi_xbp_all[i_a][bpt_k]
#define y_ij_k dipole[i_a].xy
void phasespace_ij_k(vector<vector<fourvector> > & psi_xbp_all, vector<vector<double> > & psi_xbs_all, vector<vector<double> > & psi_xbsqrts_all, int i_a, vector<dipole_set> & dipole){
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
  /*
  for (int i = 0; i < osi_p_parton[i_a].size(); i++){
    int xbi = intpow(2, i - 1);
    osi_p_parton[i_a][i] = psi_xbp_all[i_a][xbi];
  }
  */
  logger << LOG_DEBUG_VERBOSE << "phasespace_ij_k finished" << endl;
}
#undef p_i 
#undef p_j 
#undef p_k 
#undef pt_ij 
#undef pt_k 

#define p_i psi_xbp_all[0][bp_i]
#define p_j psi_xbp_all[0][bp_j]
#define p_a psi_xbp_all[0][bp_a]
#define pt_ij psi_xbp_all[i_a][bpt_ij]
#define pt_a psi_xbp_all[i_a][bpt_a]
#define x_ij_a dipole[i_a].xy
void phasespace_ij_a(vector<vector<fourvector> > & psi_xbp_all, vector<vector<double> > & psi_xbs_all, vector<vector<double> > & psi_xbsqrts_all, int i_a, vector<dipole_set> & dipole){
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
  /*
  for (int i = 0; i < osi_p_parton[i_a].size(); i++){
    int xbi = intpow(2, i - 1);
    osi_p_parton[i_a][i] = psi_xbp_all[i_a][xbi];
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

#define p_a psi_xbp_all[0][bp_a]
#define p_i psi_xbp_all[0][bp_i]
#define p_k psi_xbp_all[0][bp_k]
#define pt_ai psi_xbp_all[i_a][bpt_ai]
#define pt_k psi_xbp_all[i_a][bpt_k]
#define x_ik_a dipole[i_a].xy
/*
#define p_a psi_xbp_all[0][dipole[i_a].binary_R_emitter_1()]
#define p_i psi_xbp_all[0][dipole[i_a].binary_R_emitter_2()]
#define p_k psi_xbp_all[0][dipole[i_a].binary_R_spectator()]
#define pt_ai psi_xbp_all[i_a][dipole[i_a].binary_A_emitter()]
#define pt_k psi_xbp_all[i_a][dipole[i_a].binary_R_spectator()]
*/
void phasespace_ai_k(vector<vector<fourvector> > & psi_xbp_all, vector<vector<double> > & psi_xbs_all, vector<vector<double> > & psi_xbsqrts_all, int i_a, vector<dipole_set> & dipole){
  static Logger logger("phasespace_ai_k");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  int bp_a = dipole[i_a].binary_R_emitter_1();
  int bp_i = dipole[i_a].binary_R_emitter_2();
  int bp_k = dipole[i_a].binary_R_spectator();
  int bpt_ai = dipole[i_a].binary_A_emitter();
  int bpt_k = dipole[i_a].binary_A_spectator();

  x_ik_a = 1. - (p_i * p_k) / ((p_i + p_k) * p_a);
  //  double u_i = (p_i * p_a) / ((p_i + p_k) * p_a);
  pt_k = p_i + p_k - (1 - x_ik_a) * p_a;
  pt_ai = x_ik_a * p_a;

  for (int xbi = 1; xbi < psi_xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_i){}
    else if (xbi == bp_k){}
    else if (xbi == bp_a){}
    else if (xbi < bp_i){psi_xbp_all[i_a][xbi] = psi_xbp_all[0][xbi];}
    else if (xbi > bp_i){psi_xbp_all[i_a][xbi / 2] = psi_xbp_all[0][xbi];}
    else {cout << "Should not happen!" << endl;}
  }

  psi_xbp_all[i_a][0] = psi_xbp_all[i_a][1] + psi_xbp_all[i_a][2];
  psi_xbs_all[i_a][0] = x_ik_a * psi_xbs_all[0][0];
  psi_xbsqrts_all[i_a][0] = sqrt(psi_xbs_all[i_a][0]);
  int xbn_all = psi_xbp_all[i_a].size() - 4;
  psi_xbp_all[i_a][xbn_all] = psi_xbp_all[i_a][0];
  psi_xbs_all[i_a][xbn_all] = psi_xbs_all[i_a][0];
  psi_xbsqrts_all[i_a][xbn_all] = psi_xbsqrts_all[i_a][0];
  /*
  for (int i = 0; i < osi_p_parton[i_a].size(); i++){
    int xbi = intpow(2, i - 1);
    osi_p_parton[i_a][i] = psi_xbp_all[i_a][xbi];
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
#define p_a psi_xbp_all[0][dipole[i_a].binary_R_emitter_1()]
#define p_i psi_xbp_all[0][dipole[i_a].binary_R_emitter_2()]
#define p_b psi_xbp_all[0][dipole[i_a].binary_R_spectator()]
#define bpt_ai bp_a
#define bpt_b bp_b
#define pt_ai psi_xbp_all[i_a][dipole[i_a].binary_A_emitter()]
#define pt_b psi_xbp_all[i_a][dipole[i_a].binary_R_spectator()]
*/
#define bpt_ai bp_a
#define bpt_b bp_b
#define p_a psi_xbp_all[0][bp_a]
#define p_i psi_xbp_all[0][bp_i]
#define p_b psi_xbp_all[0][bp_b]
#define pt_ai psi_xbp_all[i_a][bpt_ai]
#define pt_b psi_xbp_all[i_a][bpt_b]
#define x_i_ab dipole[i_a].xy
void phasespace_ai_b(vector<vector<fourvector> > & psi_xbp_all, vector<vector<double> > & psi_xbs_all, vector<vector<double> > & psi_xbsqrts_all, int i_a, vector<dipole_set> & dipole){
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
  //  cout << "phasespace_ai_b   psi_xbp_all.size() = " << psi_xbp_all.size() << endl;
  //  cout << "phasespace_ai_b   psi_xbp_all[0].size() = " << psi_xbp_all[0].size() << endl;
  //  cout << "phasespace_ai_b   psi_xbp_all[" << i_a << "].size() = " << psi_xbp_all[i_a].size() << endl;

  for (int xbi = 1; xbi < psi_xbp_all[0].size(); xbi = xbi * 2){
    if      (xbi == bp_a){}
    else if (xbi == bp_b){psi_xbp_all[i_a][xbi] = psi_xbp_all[0][xbi];}
    else if (xbi == bp_i){}
    else if (xbi < bp_i){psi_xbp_all[i_a][xbi] = psi_xbp_all[0][xbi] - ((2. * (psi_xbp_all[0][xbi] * KKt)) / KKt2) * KKt + (2. * (psi_xbp_all[0][xbi] * K) / K2) * Kt;}
    else if (xbi > bp_i){psi_xbp_all[i_a][xbi / 2] = psi_xbp_all[0][xbi] - ((2. * (psi_xbp_all[0][xbi] * KKt)) / KKt2) * KKt + (2. * (psi_xbp_all[0][xbi] * K) / K2) * Kt;}
    else {cout << "Should not happen!" << endl;}
  }
  psi_xbp_all[i_a][0] = psi_xbp_all[i_a][1] + psi_xbp_all[i_a][2];
  psi_xbs_all[i_a][0] = x_i_ab * psi_xbs_all[0][0];
  psi_xbsqrts_all[i_a][0] = sqrt(psi_xbs_all[i_a][0]);
  int xbn_all = psi_xbp_all[i_a].size() - 4;

  psi_xbp_all[i_a][xbn_all] = psi_xbp_all[i_a][0];
  psi_xbs_all[i_a][xbn_all] = psi_xbs_all[i_a][0];
  psi_xbsqrts_all[i_a][xbn_all] = psi_xbsqrts_all[i_a][0];
  /*
  for (int i = 0; i < osi_p_parton[i_a].size(); i++){
    int xbi = intpow(2, i - 1);
    osi_p_parton[i_a][i] = psi_xbp_all[i_a][xbi];
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

