// dipolesubtraction/CS.phasespace.dipoles.cpp
// !!! shift into phasespace_set !!!
void phasespace_ij_k(vector<vector<fourvector> > & xbp_RA, vector<vector<double> > & xbs_RA, vector<vector<double> > & xbsqrts_RA, int i_a, vector<dipole_set> & CS_dipole);
void phasespace_ij_a(vector<vector<fourvector> > & xbp_RA, vector<vector<double> > & xbs_RA, vector<vector<double> > & xbsqrts_RA, int i_a, vector<dipole_set> & CS_dipole);
void phasespace_ai_k(vector<vector<fourvector> > & xbp_RA, vector<vector<double> > & xbs_RA, vector<vector<double> > & xbsqrts_RA, int i_a, vector<dipole_set> & CS_dipole);
void phasespace_ai_b(vector<vector<fourvector> > & xbp_RA, vector<vector<double> > & xbs_RA, vector<vector<double> > & xbsqrts_RA, int i_a, vector<dipole_set> & CS_dipole);

// dipolesubtraction/CDST.phasespace.dipoles.cpp
// !!! shift into phasespace_set !!!
void phasespace_ij_k_massive(vector<vector<fourvector> > & xbp_RA, vector<vector<double> > & xbs_RA, vector<vector<double> > & xbsqrts_RA, int i_a, vector<dipole_set> & CDST_dipole);
void phasespace_ij_a_massive(vector<vector<fourvector> > & xbp_RA, vector<vector<double> > & xbs_RA, vector<vector<double> > & xbsqrts_RA, int i_a, vector<dipole_set> & CDST_dipole);



// dipolesubtraction/CS.phasespace.group.dipoles.cpp
// !!! shift into phasespace_set - if generic works there !!!
void ac_psp_RA_group_ij_k(vector<double> & r, int channel, vector<dipole_set> & dipole, int i_a, double & g_IS_channel, phasespace_set & pset, call_generic & generic);
void ag_psp_RA_group_ij_k(vector<dipole_set> & dipole, int i_a, phasespace_set & pset, call_generic & generic);
void ac_psp_RA_group_ij_a(vector<double> & r, int channel, vector<dipole_set> & dipole, double & sinx_ij_a_min, int i_a, double & g_IS_channel, phasespace_set & pset, call_generic & generic);
void ag_psp_RA_group_ij_a(double & sinx_ij_a_min, vector<dipole_set> & dipole, int i_a, phasespace_set & pset, call_generic & generic);
void ac_psp_RA_group_ai_k(vector<double> & r, int channel, vector<dipole_set> & dipole, double & sinx_ik_a_min, int i_a, double & g_IS_channel, phasespace_set & pset, call_generic & generic);
void ag_psp_RA_group_ai_k(double & sinx_ik_a_min, vector<dipole_set> & dipole, int i_a, phasespace_set & pset, call_generic & generic);
void ac_psp_RA_group_ai_b(vector<double> & r, int channel, vector<dipole_set> & dipole, double & sinx_i_ab_min, int i_a, double & g_IS_channel, phasespace_set & pset, call_generic & generic);
void ag_psp_RA_group_ai_b(double & sinx_i_ab_min, vector<dipole_set> & dipole, int i_a, phasespace_set & pset, call_generic & generic);

// dipolesubtraction/CDST.phasespace.group.dipoles.cpp
// !!! shift into phasespace_set - if generic works there !!!
void ac_psp_RA_group_ij_k_massive(vector<double> & r, int channel, vector<dipole_set> & dipole, int i_a, double & g_IS_channel, phasespace_set & pset, call_generic & generic);
void ag_psp_RA_group_ij_k_massive(vector<dipole_set> & dipole, int i_a, phasespace_set & pset, call_generic & generic);
void ac_psp_RA_group_ij_a_massive(vector<double> & r, int channel, vector<dipole_set> & dipole, double & sinx_ij_a_min, int i_a, double & g_IS_channel, phasespace_set & pset, call_generic & generic);
void ag_psp_RA_group_ij_a_massive(double & sinx_ij_a_min, vector<dipole_set> & dipole, int i_a, phasespace_set & pset, call_generic & generic);
void ac_psp_RA_group_ai_k_massive(vector<double> & r, int channel, vector<dipole_set> & dipole, double & sinx_ik_a_min, int i_a, double & g_IS_channel, phasespace_set & pset, call_generic & generic);
void ag_psp_RA_group_ai_k_massive(double & sinx_ik_a_min, vector<dipole_set> & dipole, int i_a, phasespace_set & pset, call_generic & generic);



// dipolesubtraction/QCD.CS.dipoles.cpp
void QCD_determine_dipoles(vector<dipole_set> & CS_dipole_candidate, vector<int> & type_parton, vector<int> & basic_type_parton);
void QCD_selection_dipoles(vector<dipole_set> & dipole, vector<dipole_set> & dipole_candidate, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference, vector<vector<double> > & singular_region, vector<vector<string> > & singular_region_name, vector<vector<int> > & list_singular_regions, phasespace_set & psi, call_generic & generic);
void QCD_selection_fake_dipoles(vector<dipole_set> & dipole, vector<dipole_set> & dipole_candidate, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference, vector<vector<double> > & RA_singular_region, vector<vector<string> > & RA_singular_region_name, vector<vector<int> > & RA_singular_region_list, phasespace_set & psi, call_generic & generic, void (*determination_no_subprocess_doubledipole)(int & no_map, vector<int> & o_map, int & no_prc, vector<int> & o_prc, double & factor_symmetry, vector<int> & tp, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference), int (*determination_MCchannels_doubledipole)(int no_ps, phasespace_set & psi));
void QCD_selection_phasespace_singularity(vector<vector<double> > & RA_singular_region, vector<vector<string> > & RA_singular_region_name, vector<vector<int> > & RA_singular_region_list, phasespace_set & psi);



// dipolesubtraction/QEW.CS.dipoles.cpp
void QEW_determine_dipoles(vector<dipole_set> & QEW_dipole_candidate, vector<int> & type_parton, vector<int> & basic_type_parton);
void QEW_selection_dipoles(vector<dipole_set> & dipole, vector<dipole_set> & dipole_candidate, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference, vector<vector<double> > & singular_region, vector<vector<string> > & singular_region_name, vector<vector<int> > & list_singular_regions, phasespace_set & psi, call_generic & generic);



// dipolesubtraction/dipolesubtraction.ME2.VA.cpp
void calculate_ME2_ioperator_VA_QCD(observable_set & oset);
void calculate_ME2_ioperator_VA_QEW(observable_set & oset);
void calculate_ME2_ioperator_VA_MIX(observable_set & oset);
void calculate_ME2_VA_QCD(observable_set & oset);
void calculate_ME2_VA_QEW(observable_set & oset);
void calculate_ME2_VA_MIX(observable_set & oset);
void calculate_ME2check_VA_QCD(observable_set & oset);
void calculate_ME2check_VA_QEW(observable_set & oset);
void calculate_ME2check_VA_MIX(observable_set & oset);



// dipolesubtraction/dipolesubtraction.ME2.CA.cpp
void calculate_ME2_CA_QCD(observable_set & oset);
void calculate_ME2_CA_QEW(observable_set & oset);
void calculate_ME2check_CA_QCD(observable_set & oset, phasespace_set & psi);
void calculate_ME2check_CA_QEW(observable_set & oset, phasespace_set & psi);



// dipolesubtraction/dipolesubtraction.ME2.RA.cpp
void calculate_ME2_RA_QCD(observable_set & oset);
void calculate_ME2_RA_QEW(observable_set & oset);
void calculate_ME2_RA_MIX(observable_set & oset);
void calculate_ME2check_RA_QCD(phasespace_set & psi, observable_set & oset);
void calculate_ME2check_RA_QEW(phasespace_set & psi, observable_set & oset);
void calculate_ME2check_RA_MIX(phasespace_set & psi, observable_set & oset);



// dipolesubtraction/splitting.cpp
double P_qg(double x);
double P_gq(double x);
double P_qq_reg(double x);
double P_qq(double x);
double P_qq_plus(double x);
double P_qq_delta();
double intP_qq_plus(double x);
double Pxx_qq(double x);
double Pxx_qq_reg(double x); // == Pxx_qq
double Pxx_qq_plus(double x);
double Pxx_qq_delta();
double intPxx_qq_plus(double x);
double P_gg(double x);
double P_gg_reg(double x); // == P_gg
double P_gg_plus(double x);
double P_gg_delta(int N_f);
double intP_gg_plus(double x);
double Kbar_qg(double x);
double Kbar_gq(double x);
double Kbar_qq(double x);
double Kbar_qq_plus(double x);
double Kbar_qq_delta();
double intKbar_qq_plus(double x);
double Kbar_gg(double x);
double Kbar_gg_plus(double x);
double Kbar_gg_delta(int N_f);
double intKbar_gg_plus(double x);
double Kt_qg(double x);
double Kt_gq(double x);
double Kt_qq(double x);
double Kt_qq_plus(double x);
double Kt_qq_delta();
double intKt_qq_plus(double x);
double Kt_gg(double x);
double Kt_gg_plus(double x);
double Kt_gg_delta();
double intKt_gg_plus(double x);
double gamma_g(int N_f);
double gamma_g_ferm(int N_f);
double K_g (int N_f);
double K_g_ferm (int N_f);
double gamma_a(int N_f);
double K_a (int N_f);
double gamma_qew_a(int N_f);
double K_qew_a(int N_f);
double P_qew_qa(double x);
double P_qew_aq(double x);
double P_qew_qq_reg(double x);
double P_qew_qq(double x);
double P_qew_qq_plus(double x);
double P_qew_qq_delta();
double intP_qew_qq_plus(double x);
double Pxx_qew_qq(double x);
double Pxx_qew_qq_reg(double x);
double Pxx_qew_qq_plus(double x);
double Pxx_qew_qq_delta();
double intPxx_qew_qq_plus(double x);
double P_qew_aa(double x);
double P_qew_aa_reg(double x);
double P_qew_aa_plus(double x);
double P_qew_aa_delta(int N_f);
double intP_qew_aa_plus(double x);
double Kbar_qew_qa(double x);
double Kbar_qew_aq(double x);
double Kbar_qew_qq(double x);
double Kbar_qew_qq_plus(double x);
double Kbar_qew_qq_delta();
double intKbar_qew_qq_plus(double x);
double Kbar_qew_aa(double x);
double Kbar_qew_aa_plus(double x);
double Kbar_qew_aa_delta(int N_f);
double intKbar_qew_aa_plus(double x);
double Kt_qew_qa(double x);
double Kt_qew_aq(double x);
double Kt_qew_qq(double x);
double Kt_qew_qq_plus(double x);
double Kt_qew_qq_delta();
double intKt_qew_qq_plus(double x);
double Kt_qew_aa(double x);
double Kt_qew_aa_plus(double x);
double Kt_qew_aa_delta();
double intKt_qew_aa_plus(double x);



