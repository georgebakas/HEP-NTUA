#ifndef MODEL_SET_H
#define MODEL_SET_H

#include <string>
#include <vector>
#include "logger.h"

using namespace std;

class model_set{
private:

public:
//
// constructors
//
  model_set();

  model_set(vector<string> & file_input, int N_f, int N_f_active, double scale_ren, int present_type_perturbative_order, vector<string> & collection_type_perturbative_order, vector<string> & contribution_LHAPDFname, vector<int> & contribution_LHAPDFsubset, int switch_alpha_CMS, int switch_costheta_real);

  void calculate_Gamma_W(int QCD_order);
  void calculate_Gamma_Z(int QCD_order);
  void calculate_Gamma_H(int QCD_order);
  void calculate_Gamma_t(int QCD_order);
  void determine_CKM_matrix();
  void determine_sintheta_real(vector<string> & det_M_gauge, vector<string> & det_theta_w);
  void determine_alpha_e_Gmu(vector<string> & det_EW_coupling, int switch_alpha_CMS);
  void determine_EWcouplings_real();
  void determine_EWcouplings_complex(int switch_costheta_real);
  void fill_particle_mass_vector();

  void get_simple_userinput_from_readin(vector<string> & user_variable, vector<string> & user_value, vector<string> & readin);
  void initialization_LHAPDF(string & this_lhapdfname, int this_lhapdfsubset);

  void output_model_file(int switch_costheta_real);
  
//
// access elements
//

  double empty_double;
  double_complex empty_double_complex;
  int empty_int;
  string empty_string;

  int ew_scheme;
  // ew_scheme  is corresponds to OpenLoops switch; it is handed over to OpenLoops as well
  // ew_scheme =  0: alpha(0) scheme
  // ew_scheme =  1: Gmu scheme
  // ew_scheme =  2: alpha(MZ) scheme

  int use_cms;
  // use_cms  is corresponds to OpenLoops switch; it is handed over to OpenLoops as well
  // use_cms =  0: real-mass scheme (only works without intermediate W/Z ersonances)
  // use_cms =  1: complex-mass scheme (alpha_e_Gmu from real M_W and cos_w)
  // use_cms =  2: complex-mass scheme (alpha_e_Gmu from absolute value of complex M_W and cos_w)

  int use_adapted_ew_coupling;
  // use_adapted_ew_coupling  switches to alpha_e_0 (alpha_e_EW_Gmu/MZ) for all Born-level
  // EW couplings that are (not) connected with outgoing identified photons: can be used
  // only for ew_scheme = 1, 2 (0).
  // use_adapted_ew_coupling = -1: use alpha_e as defined by ew_scheme everywhere
  // use_adapted_ew_coupling =  0: use alpha_e_0 for couplings associated with
  //                               identified outgoing photons (ew_scheme = 1. 2)
  // use_adapted_ew_coupling =  1: use alpha_e_Gmu for couplings not associated with
  //                               identified outgoing photons (ew_scheme = 0)
  // use_adapted_ew_coupling =  2: use alpha_e_MZ for couplings not associated with
  //                               identified outgoing photons (ew_scheme = 0)
  
  double M_W;
  double M_Z;
  double M_H;
  double M_t;
  double M_b;
  double M_c;
  double M_s;
  double M_u;
  double M_d;
  double M_e;
  double M_mu;
  double M_tau;
  double M_ve;
  double M_vm;
  double M_vt;

  double M2_W;
  double M2_Z;
  double M2_H;
  double M2_t;
  double M2_b;
  double M2_c;
  double M2_s;
  double M2_u;
  double M2_d;
  double M2_e;
  double M2_mu;
  double M2_tau;
  double M2_ve;
  double M2_vm;
  double M2_vt;

  double Gamma_W;
  double Gamma_Z;
  double Gamma_H;
  double Gamma_t;
  double Gamma_b;
  double Gamma_c;

  double map_Gamma_W;
  double map_Gamma_Z;
  double map_Gamma_H;
  double map_Gamma_t;
  double map_Gamma_b;
  double map_Gamma_c;
  
  double reg_Gamma_W;
  double reg_Gamma_Z;
  double reg_Gamma_H;
  double reg_Gamma_t;
  double reg_Gamma_b;
  double reg_Gamma_c;

  double_complex cM_W;
  double_complex cM_Z;
  double_complex cM_H;
  double_complex cM_t;
  double_complex cM_b;
  double_complex cM_c;

  double_complex cM2_W;
  double_complex cM2_Z;
  double_complex cM2_H;
  double_complex cM2_t;
  double_complex cM2_b;
  double_complex cM2_c;

  double G_F;

  double Gamma_Wlv;
  double Gamma_Wud;
  double Gamma_Zll;
  double Gamma_Zvv;
  double Gamma_Zdd;
  double Gamma_Zuu;
  double Gamma_tWb;

  double BR_Wlv;
  double BR_Wud;
  double BR_Zll;
  double BR_Zvv;
  double BR_Zdd;
  double BR_Zuu;
  double BR_tWb;

  double alpha_s;
  double g_s;

  double alpha_e;
  double e;
  double alpha_e_0;
  double alpha_e_MZ;
  double alpha_e_Gmu;

  double cos_w;
  double cos2_w;
  double sin2_w;
  double sin_w;

  double_complex ccos_w;
  double_complex ccos2_w;
  double_complex csin2_w;
  double_complex csin_w;

  string CKM_matrix;
  double theta_c;
  vector<vector<double_complex> > V_ckm;
  double V_du;
  double V_dc;
  double V_dt;
  double V_su;
  double V_sc;
  double V_st;
  double V_bu;
  double V_bc;
  double V_bt;

  vector<double> e_pow;
  vector<double> M;
  vector<double> M2;
  vector<double> Gamma;
  vector<double_complex> cM;
  vector<double_complex> cM2;
  vector<double> map_Gamma;
  vector<double> reg_Gamma;


  double Q_u;
  double Q_d;
  double Q_v;
  double Q_l;
  double Iw_u;
  double Iw_d;
  double Iw_n;
  double Iw_l;


  double Cplus_Zuu;
  double Cplus_Zdd;
  double Cminus_Zuu;
  double Cminus_Zdd;
  double Cplus_Znn;
  double Cplus_Zee;
  double Cminus_Znn;
  double Cminus_Zee;
  double C_ZWminusWplus;

  double Cplus_Auu;
  double Cplus_Add;
  double Cminus_Auu;
  double Cminus_Add;
  double Cplus_Ann;
  double Cplus_Aee;
  double Cminus_Ann;
  double Cminus_Aee;
  double C_AWminusWplus;

  //  double Cplus_W; == 0 in SM
  double Cminus_W;

  double_complex cCplus_Zuu;
  double_complex cCplus_Zdd;
  double_complex cCminus_Zuu;
  double_complex cCminus_Zdd;
  double_complex cCplus_Znn;
  double_complex cCplus_Zee;
  double_complex cCminus_Znn;
  double_complex cCminus_Zee;
  double_complex cC_ZWminusWplus;

  double_complex cCplus_Auu;
  double_complex cCplus_Add;
  double_complex cCminus_Auu;
  double_complex cCminus_Add;
  double_complex cCplus_Ann;
  double_complex cCplus_Aee;
  double_complex cCminus_Ann;
  double_complex cCminus_Aee;
  double_complex cC_AWminusWplus;

  double_complex cCminus_W;

  vector<int> vGamma_W_order;
  vector<int> valpha_S_W_order;
  vector<double> valpha_S_W;
  vector<double> valpha_S_W_scale;
  vector<double> vGamma_W;
  vector<double> vGamma_Wud;
  vector<double> vGamma_Wlv;
  vector<double> vBR_Wud;
  vector<double> vBR_Wlv;
  double alpha_S_W;

  vector<int> vGamma_Z_order;
  vector<int> valpha_S_Z_order;
  vector<double> valpha_S_Z;
  vector<double> valpha_S_Z_scale;
  vector<double> vGamma_Z;
  vector<double> vGamma_Zuu;
  vector<double> vGamma_Zdd;
  vector<double> vGamma_Zvv;
  vector<double> vGamma_Zll;
  vector<double> vBR_Zuu;
  vector<double> vBR_Zdd;
  vector<double> vBR_Zvv;
  vector<double> vBR_Zll;
  double alpha_S_Z;

  vector<double> vGamma_H;

  vector<int> vGamma_t_order;
  vector<int> valpha_S_t_order;
  vector<double> valpha_S_t;
  vector<double> valpha_S_t_scale;
  vector<double> vGamma_t;
  vector<double> vGamma_tWb;
  vector<double> vBR_tWb;
  double alpha_S_t;

  vector<string> contribution_LHAPDFname;
  vector<int> contribution_LHAPDFsubset;
};

#endif
