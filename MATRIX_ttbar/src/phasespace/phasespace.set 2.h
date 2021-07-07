#ifndef PHASESPACESET_H
#define PHASESPACESET_H

#include <string>
#include <vector>

#include "randommanager.h"
#include "../classes/header/logger.h"
#include "../include/ftypes.h"
#include "../classes/header/dddistribution.h"
#include "../classes/header/xdistribution.h"
#include "multichannel.set.h"
#include "importancesampling.set.h"


/*
// to be eliminated again:
  int n_alpha_phases;
  int n_tau_phases;
  int n_x1x2_phases;
  int n_z1z2_phases;
  int n_IS_phases;
  int n_warmup;
  int warmup;
*/

using namespace std;

class call_generic;

class phasespace_set{
private:

public:
////////////////////
//  constructors  //
////////////////////
  phasespace_set();
  phasespace_set(inputparameter_set & isi, contribution_set & _csi);

  void initialization_random_number_generator(inputparameter_set & isi);

  void initialization_mapping_parameter(inputparameter_set & isi);
  void initialization_minimum_tau(inputparameter_set & isi);
  void initialization_minimum_phasespacecut(inputparameter_set & isi);

  void initialization_complete(inputparameter_set & isi);
  void initialization_complete_RA(inputparameter_set & isi, vector<dipole_set> & dipole);

  void initialization_optimization(inputparameter_set & isi);
  void initialization_filename(inputparameter_set & isi);

  void initialization_optimization_grid();
  
  void initialization_contribution_order(contribution_set & csi);

  void initialization_masses(model_set & msi);
  void initialization_masses(vector<double> _M, vector<double> _M2, vector<double> _Gamma, vector<double_complex> _cM2, vector<double> _map_Gamma, vector<double> _reg_Gamma);

  void initialization_phasespace_born();
  void initialization_phasespace_subprocess();
  void initialization_phasespace_subprocess_optimization();

  void initialization_phasespace_RA();
  //  void initialization_phasespace_RA(int _n_dipoles, vector<dipole_set> _dipole);

  void initialization_phasespace_MC_tau();
  void initialization_phasespace_MC_x_dipole_RA();

  void initialization_phasespace_subprocess_dipole_IS_RA();
  void initialization_phasespace_subprocess_dipole_RA();
  void initialization_phasespace_subprocess_optimization_dipole_RA();
  //  void initialization_phasespace_subprocess_dipole_IS_RA(int _n_dipoles, vector<dipole_set> _dipole);
  //  void initialization_phasespace_subprocess_dipole_RA(int _n_dipoles, vector<dipole_set> _dipole);
  //  void initialization_phasespace_subprocess_optimization_dipole_RA(int _n_dipoles, vector<dipole_set> _dipole);

  void initialization_phasespace_IS_CA();
  void initialization_phasespace_IS();
  void initialization_phasespace_IS_CX();
  void initialization_phasespace_IS_QT();
  void initialization_phasespace_IS_CT_from_CX();


  void initialization_phasespace_IS_CT();
  void initialization_phasespace_IS_VT();

  void initialization_MC();
  //  void initialization_MC(randommanager _random_manager);


  void determine_dipole_phasespace_RA(vector<dipole_set> & _dipole);

  // dipolesubtraction/header/CS.phasespace.dipoles.h
  void phasespace_ij_k(int i_a, vector<dipole_set> & CS_dipole);
  void phasespace_ij_a(int i_a, vector<dipole_set> & CS_dipole);
  void phasespace_ai_k(int i_a, vector<dipole_set> & CS_dipole);
  void phasespace_ai_b(int i_a, vector<dipole_set> & CS_dipole);
  
  // dipolesubtraction/header/CDST.phasespace.dipoles.h
  void phasespace_ij_k_massive(int i_a, vector<dipole_set> & CDST_dipole);
  void phasespace_ij_a_massive(int i_a, vector<dipole_set> & CDST_dipole);



  void correct_phasespacepoint_real();
  void correct_phasespacepoint_dipole();




  void calculate_g_tot();


  void calculate_IS();
  void calculate_IS_RA(vector<dipole_set> & _dipole);
  void calculate_IS_tau_x1x2_MC();
  void calculate_initial_tau_IS_x1x2_IS();
  void calculate_initial_tau_IS_x1x2();
  void calculate_initial_tau_x1x2_IS();
  void calculate_initial_tau_x1x2();
  void calculate_initial_tau_fixed_x1x2_IS();
  void calculate_initial_tau_fixed_x1x2();

  void calculate_initial_collinear_z1z2_IS();

  void calculate_IS_CT();
  void calculate_IS_VT();

  void calculate_IS_QT_from_CX();
  void calculate_IS_CA_from_CT();
  void calculate_IS_CX_from_CT();
  void calculate_IS_CX_from_QT();

  void calculate_IS_CX();
  void calculate_IS_QT();


  double c_phi(int i_r);
  double inv_phi(double & phi);
  double c_costheta(int i_r);
  double inv_costheta(double & costheta);

  void c_propagator(int i_a, int m, int out, double smin, double smax, int i_r1);
  double g_propagator(int i_a, int m, int out, double smin, double smax);
  void inv_propagator(int i_a, int m, int out, double smin, double smax, double & r1);

  double c_propagator_Breit_Wigner(double & r1, int m, double & smin, double & smax);
  double g_propagator_Breit_Wigner(double & s, int m, double & smin, double & smax);
  double inv_propagator_Breit_Wigner(double s, int m, double & smin, double & smax);

  double c_propagator_narrow_width();
  double g_propagator_narrow_width(int m);
  double inv_propagator_narrow_width();

  double c_propagator_vanishing_width(double r, double mass2, double smin, double smax, double nu);
  double g_propagator_vanishing_width(double s, double mass2, double smin, double smax, double nu);
  double inv_propagator_vanishing_width(double s, double mass2, double smin, double smax, double nu);
  
  void c_timelikeinvariant(int i_a, int out, double smin, double smax, int i_r1);
  double g_timelikeinvariant(double smin, double smax);
  void inv_timelikeinvariant(int i_a, int out, double smin, double smax, double & r1);


  double c_t_propagator(double r, int m, double smin, double smax);
  double g_t_propagator(double s, int m, double smin, double smax);
  double inv_t_propagator(double s, int m, double smin, double smax);



  void c_decay(int i_a, int in, int out1, int out2, int i_r1, int i_r2);
  double g_decay(int i_a, int in, int out1, int out2);
  void inv_decay(int i_a, int in, int out1, int out2, double & r1, double & r2);

  void c_tchannel(int i_a, int m, int in, int in1, int in2, int out1, int out2, int i_r1, int i_r2);
  double g_tchannel(int i_a, int m, int in, int in1, int in2, int out1, int out2);
  void inv_tchannel(int i_a, int m, int in, int in1, int in2, int out1, int out2, double & r1, double & r2);

  void c_tchannel_opt(int i_a, int m, int out, int in1, int in2, int out1, int out2, int i_r1, int i_r2);
  double g_tchannel_opt(int i_a, int m, int out, int in1, int in2, int out1, int out2);
  void inv_tchannel_opt(int i_a, int m, int out, int in1, int in2, int out1, int out2, double & r1, double & r2);


  ///  void psp_MCweight_optimization(observable_set & oset);
  ///  void step_MCweight_optimization(observable_set & oset);
  ///  void result_MCweight_optimization(observable_set & oset, int & size_proc_generic, vector<int> & size_proceeding);
  ///  void step_tau_weight_optimization(observable_set & oset, int & size_proc_generic, vector<int> & size_proceeding);
  ///  void step_x1x2_weight_optimization(observable_set & oset, int & size_proc_generic, vector<int> & size_proceeding);
  ///  void step_z1z2_weight_optimization(observable_set & oset, int & size_proc_generic, vector<int> & size_proceeding);


  void psp_MCweight_optimization(double & integrand, double & this_psp_weight, double & this_psp_weight2, double & optimization_modifier);
  void psp_MCweight_optimization_CA(double & integrand, double & this_psp_weight, double & this_psp_weight2);

  void step_MCweight_optimization();
  void result_MCweight_optimization();
  void weight_output_MCweight_optimization();

  void step_tau_weight_optimization();
  void step_x1x2_weight_optimization();
  void step_z1z2_weight_optimization();

  void psp_optimization_complete(double & integrand, double & this_psp_weight, double & this_psp_weight2, double & optimization_modifier);
  void step_optimization_complete();
  void result_optimization_complete();
  void output_optimization_complete();

  void output_check_tau_0();

  void handling_cut_psp();
  void handling_techcut_psp();

  void ac_psp_RA_group_ij_k(vector<dipole_set> & dipole, int i_a);
  void ag_psp_RA_group_ij_k(vector<dipole_set> & dipole, int i_a);
  void ac_psp_RA_group_ij_a(vector<dipole_set> & dipole, double & sinx_ij_a_min, int i_a);
  void ag_psp_RA_group_ij_a(double & sinx_ij_a_min, vector<dipole_set> & dipole, int i_a);
  void ac_psp_RA_group_ai_k(vector<dipole_set> & dipole, double & sinx_ik_a_min, int i_a);
  void ag_psp_RA_group_ai_k(double & sinx_ik_a_min, vector<dipole_set> & dipole, int i_a);
  void ac_psp_RA_group_ai_b(vector<dipole_set> & dipole, double & sinx_i_ab_min, int i_a);
  void ag_psp_RA_group_ai_b(double & sinx_i_ab_min, vector<dipole_set> & dipole, int i_a);

  void ac_psp_RA_group_ij_k_massive(vector<dipole_set> & dipole, int i_a);
  void ag_psp_RA_group_ij_k_massive(vector<dipole_set> & dipole, int i_a);
  void ac_psp_RA_group_ij_a_massive(vector<dipole_set> & dipole, double & sinx_ij_a_min, int i_a);
  void ag_psp_RA_group_ij_a_massive(double & sinx_ij_a_min, vector<dipole_set> & dipole, int i_a);
  void ac_psp_RA_group_ai_k_massive(vector<dipole_set> & dipole, double & sinx_ik_a_min, int i_a);
  void ag_psp_RA_group_ai_k_massive(double & sinx_ik_a_min, vector<dipole_set> & dipole, int i_a);


  void initialization_fake_dipole_mapping_RRA();


  call_generic *generic;

  //  vector<steering_optimization> steering;

  //  vector<int> type_parton;
  //  vector<int> basic_type_parton;
  //  vector<vector<int> > type_parton;
  ///  int n_particle;

  /*  
      int no_map;
      vector<int> o_map;
  // most likely actually not needed:
  int no_prc;
  vector<int> o_prc;
  //
  */
  vector<int> no_map;
  vector<vector<int> > o_map;
  // most likely actually not needed:
  vector<int> no_prc;
  vector<vector<int> > o_prc;
  //



  string dir_MCweights;

  randommanager random_manager;
  vector<randomvariable*> phasespace_randoms;

  vector<randomvariable*> QT_random_z;
  randomvariable* QT_random_z1;
  randomvariable* QT_random_z2;
  randomvariable* QT_random_qt;

  vector<double> sran;


  // random numbers could become part of multichannel/importancesampling sets !!!
  vector<double> r;

  vector<double> random_MC;
  vector<double> random_MC_tau;
  vector<double> random_tau;
  vector<double> random_x12;
  vector<vector<double> > random_z;


  Logger *logger; // ???




  int switch_MC;
  int switch_MC_tau;
  int switch_MC_x_dipole;

  int switch_IS_MC;
  int switch_IS_tau;
  int switch_IS_x1x2;
  int switch_IS_z1z2;

  int switch_n_events_opt;

  int switch_qTcut;


  // MC parameters fixed for the whole run, maybe up to modifications depending on optimization procedure

  vector<vector<int> > container_IS_startvalue;
  vector<string> container_IS_name;
  vector<int> container_IS_switch;

  vector<vector<vector<int> > > c_p;
  vector<vector<vector<int> > > c_f;
  vector<vector<vector<int> > > c_t;
  vector<vector<vector<int> > > c_d;


  vector<vector<double> > v_smin;
  vector<vector<double> > v_smax;

  vector<vector<double> > g_p;
  vector<vector<double> > g_f;
  vector<vector<double> > g_t;
  vector<vector<double> > g_d;

  vector<vector<vector<int> > > needed_v_smin;
  vector<vector<vector<int> > > needed_v_smax;

  vector<vector<vector<int> > > needed_g_p;
  vector<vector<vector<int> > > needed_g_f;
  vector<vector<vector<int> > > needed_g_t;
  vector<vector<vector<int> > > needed_g_d;
  
  
  // variables changing at each phase-space point

  vector<vector<fourvector> > xbp_all;
  vector<vector<double> > xbs_all;
  vector<vector<double> > xbsqrts_all;



  // MC parameters fixed for the whole run

  int n_alpha_steps;
  int n_alpha_events;
  int n_alpha_epc;

  int n_tau_steps;
  int n_tau_events;
  int n_tau_bins;

  int n_x1x2_steps;
  int n_x1x2_events;
  int n_x1x2_bins;

  int n_z1z2_steps;
  int n_z1z2_events;
  int n_z1z2_bins;

  int n_IS_events;
  int n_IS_events_factor;
  int n_IS_steps;
  int n_IS_gridsize;
  int n_IS_gridsize_p;
  int n_IS_gridsize_f;
  int n_IS_gridsize_t_t;
  int n_IS_gridsize_t_phi;
  int n_IS_gridsize_d_cth;
  int n_IS_gridsize_d_phi;
  int n_IS_gridsize_xy;
  int n_IS_gridsize_zuv;

  int n_events_max;
  int n_events_min;
  int n_step;

  int n_events_MC_opt;
  int opt_n_events_min;

  int xb_max;
  int xb_max_dipoles;
  int xb_out_dipoles;
  int no_random;

  vector<int> no_random_dipole;

  vector<vector<fourvector> > start_xbp_all;
  vector<vector<double> > start_xbs_all;
  vector<vector<double> > start_xbsqrts_all;

  vector<vector<fourvector> > corrected_xbp_all;
  vector<vector<double> > corrected_xbs_all;
  vector<vector<double> > corrected_xbsqrts_all;


  int n_dipoles;
  vector<double> dipole_sinx_min;
  vector<double> dipole_x;

  double s_had;
  double tau_0_s_had;

  double nu;
  double nuxs;
  double nuxt;
  double exp_pdf;
  double exp_pT;
  double exp_y;
  double exp_ij_k_y;
  double exp_ij_k_z;
  double exp_ij_a_x;
  double exp_ij_a_z;
  double exp_ai_k_x;
  double exp_ai_k_u;
  double exp_ai_b_x;
  double exp_ai_b_v;
  double map_technical_s;
  double map_technical_t;
  double map_technical_x;
  double mass0;

  double cut2_pT;
  double cut_pT_lep;

  vector<vector<double> > smin_opt;
  vector<vector<double> > sqrtsmin_opt;

  double E;
  double tau_0;
  int coll_choice;
  double min_qTcut;

  vector<double> mapping_cut_pT;

  double g_onshell_decay;

  vector<double> tau_0_num;

  // MC parameters modified at each optimization step

  /////  int i_alpha_it;
  /////  vector<vector<double> > alpha_it;
  /////  vector<double> diff_w; // -> MC_diff_w_step
  /////  int Dmin;
  int MCweight_lc_min;
  double a_reserved_min;
  int switch_use_alpha_after_IS;
  int switch_step_mode_grid;
  int x_alpha_it_min;


  // model parameters, fixed for the whole run
  // could always be taken from model_set msi ???
  vector<double> M;
  vector<double> M2;
  vector<double> Gamma;
  vector<double_complex> cM2;
  vector<double> map_Gamma;
  vector<double> reg_Gamma;



  // input/output filenames, fixed for the whole run
  // could be replaced by a single file ???
  string filename_MCweight;
  string filename_tauweight;
  string filename_x1x2weight;
  vector<string> filename_z1z2weight;
  string filename_MCweight_in_contribution;
  string filename_tauweight_in_contribution;
  string filename_x1x2weight_in_contribution;
  vector<string> filename_z1z2weight_in_contribution;



  // MC optimization parameters

  int weight_opt;
  int switch_IS_mode_phasespace;
  string weight_in_contribution;
  string weight_in_directory;
  int weight_min;
  double weight_limit_min;
  double weight_limit_max;


  
  // MC parameters changing on runtime (not used so far !!!)

  //  vector<double> *r;
  int MC_optswitch;


  // number of channels within respective dipole
  vector<int> MC_n_channel_phasespace;
  // number of channels summed up to respective dipole
  vector<int> MC_sum_channel_phasespace;
  // number of channels in total
  int MC_n_channel;
  // selected channel within respective dipole (referring to MC_n_channel_phasespace[x_a])
  int MC_channel_phasespace;
  // selected channel in total (referring to MC_n_channel)
  int MC_channel;
  
  double hcf;
  double ps_factor;



  double g_tot;
  double g_global_NWA;

  double g_MC;
  double g_pdf;
  double g_tau;
  double g_x1x2;

  double MC_g_IS_global;

  /////  vector<double> temp_sum_w_channel;
  /////  vector<double> sum_w_channel;
  vector<double> MC_sum_w_channel;

  /////  vector<double> w_channel_av;

  /////  vector<double> MC_alpha;
  /////  vector<double> MC_beta;

  /////  vector<double> MC_g_channel;
  vector<double> MC_g_IS_channel;

  //  vector<int> cuts_channel;
  //  vector<int> MC_n_rej_channel;

  //  vector<int> counts_channel;
  //  vector<int> MC_n_acc_channel;

  int active_optimization;
  int end_optimization;

  int MC_opt_end;
  int IS_opt_end;
  /*///
  int tau_opt_end;
  int tau_channel;
  vector<double> tau_alpha;
  vector<double> tau_beta;
  vector<int> tau_counts_channel;
  vector<int> tau_cuts_channel;
  vector<double> sum_tau_channel_weight;
  vector<double> sum_tau_channel_weight2;
  *////
  int n_tau_MC_channel;
  /*///
  int tau_MC_channel;
  *////
  ///
  vector<int> tau_MC_map;
  ///
  /*///
  vector<double> tau_MC_alpha;
  vector<double> tau_MC_beta;
  vector<vector<int> > tau_MC_tau_channel;
  vector<vector<double> > tau_MC_tau_alpha;
  vector<vector<double> > tau_MC_tau_beta;
  *////
  ///
  vector<vector<double> > tau_MC_tau_gamma;

  multichannel_set MC_phasespace;

  // IS-multichannelling with optimization
  multichannel_set MC_tau;
  /*
  vector<double> tau_IS_MC_alpha;
  vector<double> tau_IS_MC_beta;
  vector<int> tau_IS_MC_n_acc_channel;
  vector<int> tau_IS_MC_n_rej_channel;
  vector<double> tau_IS_MC_g_channel;
  vector<double> tau_IS_MC_g_IS_channel;
  vector<double> tau_IS_MC_temp_sum_w_channel;
  vector<double> tau_IS_MC_sum_w_channel;
  vector<double> tau_IS_MC_MC_sum_w_channel;
  vector<double> tau_IS_MC_w_channel_av;
  */
  // IS-multichannelling inside dipole mappings with optimization
  vector<vector<int> > MC_x_dipole_mapping;
  vector<multichannel_set> MC_x_dipole;
  /*
  vector<vector<double> > tau_IS_MC_alpha_phasespace;
  vector<vector<double> > tau_IS_MC_beta_phasespace;
  vector<vector<int> > tau_IS_MC_n_acc_channel_phasespace;
  vector<vector<int> > tau_IS_MC_n_rej_channel_phasespace;
  vector<vector<double> > tau_IS_MC_g_channel_phasespace;
  vector<vector<double> > tau_IS_MC_g_IS_channel_phasespace;
  vector<vector<double> > tau_IS_MC_temp_sum_w_channel_phasespace;
  vector<vector<double> > tau_IS_MC_sum_w_channel_phasespace;
  vector<vector<double> > tau_IS_MC_MC_sum_w_channel_phasespace;
  vector<vector<double> > tau_IS_MC_w_channel_av_phasespace;
  */

  //  vector<importancesampling_set> IS_phasespace;
  
  importancesampling_set IS_tau;

  importancesampling_set IS_x1x2;

  vector<importancesampling_set> IS_z1z2;


  
  // not used so far
  /////  vector<vector<int> > tau_MC_tau_counts_channel;
  /////  vector<vector<int> > tau_MC_tau_cuts_channel;
  /////  vector<vector<double> > sum_tau_MC_tau_channel_weight;
  /////  vector<vector<double> > sum_tau_MC_tau_channel_weight2;
  // 1st entry: no. dipole (0 = real, >0 = dipoles)
  // 2nd entry: no. channel within each dipole (0 = usual mapping, >0 = BW IS mapping)
  //

  /*///
  int x1x2_opt_end;
  int x1x2_channel;
  ///
  vector<double> x1x2_alpha;
  vector<double> x1x2_beta;
  vector<int> x1x2_counts_channel;
  vector<int> x1x2_cuts_channel;
  vector<double> sum_x1x2_channel_weight;
  vector<double> sum_x1x2_channel_weight2;
  *////
  ///
  int z1z2_opt_end;
  /*///
  vector<int> z1z2_channel;
  ///
  vector<vector<double> > z1z2_alpha;
  vector<vector<double> > z1z2_beta;
  vector<vector<int> > z1z2_counts_channel;
  vector<vector<int> > z1z2_cuts_channel;
  vector<vector<double> > sum_z1z2_channel_weight;
  vector<vector<double> > sum_z1z2_channel_weight2;
  *////
  ///
  
  double boost;
  vector<double> x_pdf;
  // SK: Is that one used ??? clarify !!!
  vector<vector<double> > xz_pdf;
  
  vector<double> z_coll;
  vector<vector<double> > all_xz_coll_pdf;
  vector<double> g_z_coll;
  
  // QT definitions - should later be removed and unified with the ones used in CA !!!
  double QT_qt2;
  //  double QT_qt;
  //  double QT_y;
  //  double QT_Q;
  double QT_jacqt2;

  double Qres;
  bool dynamical_Qres;
  double Qres_prefactor;

  bool do_resummation;
  double lambda_qt2; // for generation of qT in a resummed computation

  // SK: QT_random_qt2 needed any longer ???
  double QT_random_qt2;
  // careful: contains z1,z2, NOT x1/z1,x2/z2
  vector<double> QT_zx_pdf;
  vector<double> QT_xz_pdf;
  // SK: zx_pdf and xz_pdf not used by me: we need to unify the notation here !!!
  //  vector<double> zx_pdf;
  //  vector<double> xz_pdf;
  //
  vector<double> zz_pdf;
  vector<double> z_pdf;
  vector<vector<int> > zall_pdf;
  double QT_g_z1;
  double QT_g_z2;
  double QT_eps;
  double QT_random_IS;
  double QT_g_IS_;
  //



  long long i_gen;
  long long i_rej;
  long long i_acc;
  long long i_nan;
  long long i_tec;

  long long last_step_mode;
  long long *i_step_mode;

  int random_seed;

  string process_class;
  string subprocess;
  string type_perturbative_order;
  string type_contribution;
  string type_correction;
  vector<int> contribution_order_alpha_s;
  vector<int> contribution_order_alpha_e;
  vector<int> contribution_order_interference;
  vector<int> phasespace_order_alpha_s;
  vector<int> phasespace_order_alpha_e;
  vector<int> phasespace_order_interference;



  vector<vector<double> > RA_singular_region;
  vector<vector<string> > RA_singular_region_name;
  vector<vector<int> > RA_singular_region_list;
  int RA_techcut;
  int RA_x_a;

  double cut_technical;

  double symmetry_factor;

  user_defined user;
  contribution_set *csi;

  //  int switch_console_output_tau_0;

  int switch_console_output_phasespace_issue;

  int switch_RS_mapping;
  
  vector<dipole_set> *RA_dipole;
};

#endif
