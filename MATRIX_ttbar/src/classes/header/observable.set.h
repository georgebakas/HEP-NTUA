#ifndef OBSERVABLESET_H
#define OBSERVABLESET_H

#include <map>
#include <string>
#include <vector>

#include "randommanager.h"
#include "logger.h"
#include "ftypes.h"
//include "../../include/ftypes.h"
#include "fourvector.h"
#include "particle.h"
#include "dddistribution.h"
#include "xdistribution.h"
#include "dipole.set.h"
#include "collinear.set.h"
#include "ioperator.set.h"
#include "phasespace.set.h"
#include "../../NJsubtraction/NJsubtraction.h"

using namespace std;

class munich;

class observable_set{
 private:

 public:
////////////////////
//  constructors  //
////////////////////
  observable_set();
  observable_set(inputparameter_set & isi, contribution_set & _csi);
  //  observable_set(inputparameter_set & isi, contribution_set & _csi, munich & _munich);
  observable_set(int _n_set_TSV);
  observable_set(int _switch_qTcut, int _n_qTcut, double _min_qTcut, double _step_qTcut, double _scale_ren, double _scale_fact, int _switch_TSV, int _n_set_TSV, vector<string> _name_set_TSV, vector<double> _central_scale_ren_TSV, vector<double> _central_scale_fact_TSV, vector<double> _relative_central_scale_ren_TSV, vector<double> _relative_central_scale_fact_TSV, vector<int> _n_scale_ren_TSV, vector<int> _n_scale_fact_TSV, vector<int> _factor_scale_ren_TSV, vector<int> _factor_scale_fact_TSV, vector<int> _dynamic_scale_ren_TSV, vector<int> _dynamic_scale_fact_TSV, vector<double> _min_qTcut_TSV, vector<double> _max_qTcut_TSV, vector<int> _switch_distribution_TSV, vector<int> _switch_moment_TSV, vector<int> _max_n_integrand_TSV, vector<double> _min_qTcut_distribution_TSV, vector<double> _max_qTcut_distribution_TSV, vector<int> _no_central_scale_ren_TSV, vector<int> _no_central_scale_fact_TSV, vector<vector<double> > _relative_scale_ren_TSV, vector<vector<double> > _relative_scale_fact_TSV, vector<string> _filename_integration_TSV, vector<string> _filename_result_TSV, vector<string> _filename_distribution_TSV, vector<string> _filename_moment_TSV, int _switch_VI, int _switch_KP, int _switch_CM, int _switch_OL, int _switch_yuk, int _order_y);

  void OpenLoops_testpoint_pptt_cc();

  void initialization_complete(inputparameter_set & isi, phasespace_set & psi);
  // process independent
  void initialization_path(inputparameter_set & isi);
  void initialization_unit(inputparameter_set & isi);
  double determine_unit_factor(string & unit);

  void initialization_object_event_selection(inputparameter_set & isi);
  void initialization_switches(inputparameter_set & isi);
  void initialization_qTcut(inputparameter_set & isi);
  void initialization_LHAPDF_parameters(inputparameter_set & isi);
  void initialization_integration_parameters(inputparameter_set & isi);
  void initialization_generic(inputparameter_set & isi);
  void initialization_OpenLoops_input(inputparameter_set & isi);
  void initialization_basic_CV(inputparameter_set & isi);
  void initialization_basic_TSV(inputparameter_set & isi);
  void initialization_CV();

  // process dependent
  void initialization_contribution(contribution_set & csi);
  void initialization_object_process();
  void initialization_filename();











  void initialization_masses(vector<double> _M, vector<double> _M2);


  void initialization_OpenLoops_input(vector<string> _OL_parameter, vector<string> _OL_value);
  void initialization_OpenLoops_parameter(phasespace_set & psi);
  void initialization_OpenLoops_process(phasespace_set & psi);
  int register_OL_subprocess(int i_a, int amptype);
  void testpoint_from_OL_rambo();
  
  void initialization_filename(int process_type);//vector<int> this_type_parton, 
  //  void initialization_filename(int process_type, vector<xdistribution> all_distributions, vector<dddistribution> all_dddistributions);//vector<int> this_type_parton, 

  void initialization_object_generic();


  

  void initialization_object_process(vector<vector<int> > type_parton);
  void initialization_object_process_part2();

  void initialization_runtime_partonlevel();
  void output_initialization_runtime_partonlevel();

  void determine_phasespace_object_partonlevel();

  void initialization_integration(phasespace_set & psi);
  void initialization_TSV();
  void initialization_distribution();
  void readin_file_distribution();
  void readin_file_dddistribution();
  void determine_extended_distribution();
  void get_userinput_extravalue_from_readin(vector<string> & user_variable, vector<string> & user_variable_additional, vector<string> & user_value, vector<string> & readin);
  void distribution_back_to_back_configuration(int i_a, int k_g);

  // new: shifted from routines.distribution.cpp (and modified)
  void determine_distribution_complete();
  void determine_single_distribution();
  void determine_double_distribution();
  void preparation_distribution_TSV();
  void determine_single_distribution_TSV();
  void determine_double_distribution_TSV();


  void initialization_distribution(vector<xdistribution> _dat, vector<dddistribution> _dddat);

  void initialization_scale_fact_CV();
  void initialization_scale_ren_CV();

  void initialization_result();

  void initialization_tree();
  void initialization_VA(); //vector<vector<ioperator_set> > & _ioperator
  void initialization_CA();
  void initialization_RA(vector<dipole_set> & _dipole);

  void initialization_specific_VA();
  void initialization_specific_CA();
  void initialization_specific_RA();

  void initialization_QT();
  void initialization_QT_coefficients();
  void initialization_specific_VT();
  void initialization_specific_QT();
  void initialization_specific_CT(phasespace_set & psi);

  void initialization_NJ();

  void initialization_CX_ncollinear();
  void initialization_CX_multicollinear();

  void calculate_IS_CX(phasespace_set & psi);


  void initialization_LHAPDF();

  // selection content pdf
  void determine_selection_content_pdf();
  void perform_selection_content_pdf();
  void perform_selection_content_pdf_collinear();
  void perform_selection_content_pdf_list();

  void calculate_pdf_LHAPDF_list_CV();
  void calculate_pdf_LHAPDF_scale_list_CV(int sd, int ss);

  void calculate_pdf_LHAPDF_list_TSV();
  void calculate_pdf_LHAPDF_scale_list_TSV(int i_a, int i_v, int i_m);


  void calculate_pdf_LHAPDF_CA_collinear_TSV(vector<vector<double> > & all_xz_coll_pdf);
  void calculate_pdf_LHAPDF_CA_collinear_CV(vector<vector<double> > & all_xz_coll_pdf);
  void calculate_pdf_LHAPDF_scale_CA_collinear_CV(vector<vector<double> > & all_xz_coll_pdf, double & mu_fact, vector<vector<vector<double> > > & all_pdf_factor);

  void calculate_pdf_LHAPDF_TSV();
  void calculate_pdf_LHAPDF_CV();
  void calculate_pdf_LHAPDF_RA_CV();
  void calculate_pdf_LHAPDF_scale_CV(double & mu_fact, vector<double> & pdf_factor);




  void determine_object_definition();
  void determine_object_partonlevel();
  void determine_equivalent_object();
  void determine_relevant_object();
  void determine_runtime_object();

  void output_determine_runtime_object();


  void determine_p_parton(phasespace_set & psi);

  void determine_scale();
  void determine_scale_CA();
  void determine_scale_RA(int i_a);

  void determine_integrand(phasespace_set & psi);
  void determine_integrand_VA(phasespace_set & psi);
  void determine_integrand_CA(phasespace_set & psi);
  void determine_integrand_RA(phasespace_set & psi);

  void determine_correlationoperator_QCD();

  void determine_splitting_tH1F(phasespace_set & psi);
  void determine_splitting_tH1F_tH1(phasespace_set & psi);
  void determine_splitting_tH1F_tH1_without_H1_delta(phasespace_set & psi);
  void determine_splitting_tH1_only_H1_delta(phasespace_set & psi);
  void determine_splitting_tgaga_tcga_tgamma2(phasespace_set & psi);
  void determine_splitting_tH2(phasespace_set & psi);


  void determine_splitting_tH1F(vector<vector<double> > & sum_coll, vector<vector<vector<double> > > & coll, phasespace_set & psi);
  void determine_splitting_tH1F_tH1(vector<vector<double> > & sum_coll, vector<vector<vector<double> > > & coll, phasespace_set & psi);
  void determine_splitting_tH1F_tH1_scale_dependent(double & H1_delta, vector<vector<double> > & sum_coll, vector<vector<vector<double> > > & coll, phasespace_set & psi);
  void determine_splitting_tgaga_tcga_tgamma2(double & A_F, vector<vector<double> > & sum_coll, vector<vector<vector<double> > > & coll, phasespace_set & psi);
  void determine_splitting_tgaga_tcga_tgamma2_scale_dependent(double & A_F, vector<vector<double> > & sum_coll, vector<vector<vector<double> > > & coll, phasespace_set & psi);
  void determine_splitting_tH2(vector<vector<double> > & sum_coll, vector<vector<vector<double> > > & coll, phasespace_set & psi);


  void calculate_A_F();
  void calculate_B2_A_F();

  void calculate_QT_CM_kinematics();
  
  void calculate_Ft1born();
  void calculate_Gamma1(); // to be replaced...
  void calculate_Gamma1born();
  void calculate_Gamma1loop();

  void calculate_Gamma2();
  void calculate_commutator_Gamma1_Ft1();
  void calculate_anticommutator_Gamma1_Ft1();
  void calculate_Gamma1_Ft1();
  void calculate_Ft1_Gamma1();
  void calculate_Gamma1_squared();

  double calculate_sigma12(double & pdf_factor_x1x2);
  double calculate_sigma11(double & pdf_factor_x1x2, double & tH1F, double & LQ);
  double calculate_sigma24(double & pdf_factor_x1x2);
  double calculate_sigma23(double & pdf_factor_x1x2, double & sig11);
  double calculate_sigma22(double & pdf_factor_x1x2, double & sig11, double & tH1, double & tH1F, double & H1full, double & tgaga, double & LR, double & LF, double & LQ);

  double calculate_sigma21(double & pdf_factor_x1x2, double & sig11, double & tH1, double & tH1F, double & H1full, double & tgaga, double & tcga, double & tgamma2, double & LR, double & LF, double & LQ, double & A_F);
  double calculate_sigma21(double & pdf_factor_x1x2, double & sig11, double & tH1_only_H1_delta, double & tH1, double & tH1F, double & H1full, double & tgaga, double & tcga, double & tgamma2, double & LR, double & LF, double & LQ, double & A_F);

  double calculate_H1(double pdf_factor_x1x2, double tH1, double tH1F, double LR, double LF, double LQ);
  double calculate_H1(double pdf_factor_x1x2, double tH1_only_H1_delta, double tH1, double tH1F, double LR, double LF, double LQ);
  double calculate_H2(double pdf_factor_x1x2, double tH1, double tH1F, double H1full, double sig11, double tgaga, double tcga, double tgamma2, double tH2, double LR, double LF, double LQ);

  void determine_NJconvolution_terms_1(phasespace_set &psi,vector<double> &conv_coeff_Pxx,vector<double> &conv_coeff_Ixx);

  void Njettiness_calculate_NJ_axes(int i_a);
  double Njettiness_calculate_NJ_axes_assigned_energy(fourvector & temp_jet);
  void Njettiness_calculate_NJ_axes_frame(int i_a);
  //  void Njettiness_calculate_NJ_axes_normalized(int i_a);
  //  void Njettiness_calculate_NJ_axes_energy(int i_a);
  //  void Njettiness_calculate_NJ_axes_Qi(int i_a);

  
  void determine_integrand_CX_ncollinear_CT(phasespace_set & psi);
  void determine_integrand_CX_ncollinear_CT2(phasespace_set & psi);
  void determine_integrand_CX_ncollinear_VT(phasespace_set & psi);
  void determine_integrand_CX_ncollinear_VT2(phasespace_set & psi);

  void determine_integrand_CX_ncollinear_CJ(phasespace_set & psi);
  void determine_integrand_CX_ncollinear_VJ(phasespace_set & psi);

  void determine_integrand_CX_multicollinear_CT(phasespace_set & psi);
  void determine_integrand_CX_CT2(phasespace_set & psi);

  void determine_integrand_CT(phasespace_set & psi);
  void determine_integrand_CT2(phasespace_set & psi);

  void determine_integrand_VT(phasespace_set & psi);
  void determine_integrand_VT2(phasespace_set & psi);

  void determine_psp_weight();
  void determine_psp_weight_TSV(phasespace_set & psi);

  void determine_psp_weight_RA(phasespace_set & psi);
  void determine_psp_weight_QT(phasespace_set & psi);

  // to be removed later ...
  void determine_psp_weight_VT(phasespace_set & psi);
  // to be removed later ...
  void determine_psp_weight_CT(phasespace_set & psi);
  void determine_psp_weight_CT2(phasespace_set & psi);


  // classes/observable.ME2.cpp
  void calculate_ME2_born();
  void calculate_ME2check_born();
  void check_vanishing_ME2_born();
  void calculate_ME2_loop();
  void calculate_ME2check_loop();



  void initialization_runtime();
  void determine_runtime(phasespace_set & psi);

  void header_integration(ofstream & out_integration, phasespace_set & psi);

  void output_check_value_scale();
  void output_check_running_alpha_S();

  void output_integration(phasespace_set & psi);
  void output_finalization_integration(phasespace_set & psi);
  void output_integration_not_needed(phasespace_set & psi);

  void output_integrand_maximum_psp(phasespace_set & psi);
  void output_integrand_maximum(phasespace_set & psi);
  void output_integrand_maximum_VA(phasespace_set & psi);
  void output_integrand_maximum_CA(phasespace_set & psi);//, vector<vector<collinear_set> > & _collinear
  void output_integrand_maximum_RA(phasespace_set & psi); //, vector<dipole_set> & _dipole

  void output_integrand_cancellation_check_RA(phasespace_set & psi);

  void output_integrand_maximum_VT2(phasespace_set & psi);

  void output_testpoint_born(ofstream & out_comparison);

  void output_testpoint_VA_ioperator(ofstream & out_comparison);
  void output_testpoint_VA_result(ofstream & out_comparison);
  void output_testpoint_VA_Delta(ofstream & out_comparison, int i, double & i_Delta, string & s_Delta);

  void output_testpoint_CA(phasespace_set & psi);

  void output_testpoint_RA(ofstream & out_comparison);


  void output_testpoint_VT_result(ofstream & out_comparison);


  void output_distribution_complete(phasespace_set & psi);
  void output_distribution(phasespace_set & psi);
  void output_dddistribution(phasespace_set & psi);
  void output_distribution_CV(phasespace_set & psi);
  void output_dddistribution_CV(phasespace_set & psi);
  void output_distribution_TSV(phasespace_set & psi);

  void output_pdf_comparison_CA_to_CT(phasespace_set & psi);
  void output_pdf_comparison_CX_to_CT(phasespace_set & psi);
  void output_pdf_comparison_CX_to_CT2(phasespace_set & psi);
  void output_pdf_comparison_CX_multicollinear_to_CT2(phasespace_set & psi);

  // only temporarily:
  double determine_integrand_CT(double qt2, double q2, double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb, int order, double A_F, double LR, double LF, double LQ, vector<double> & sigma_qT, phasespace_set & psi);


  void determine_ioperator_QCD(phasespace_set & psi);
  void calculate_ioperator_QCD();
  void calculate_ioperator_QCD_CS();
  void calculate_ioperator_QCD_CDST();
  void determine_ioperator_QEW(phasespace_set & psi);
  void calculate_ioperator_QEW();
  void calculate_ioperator_QEW_CS();
  void calculate_ioperator_QEW_CDST();

  void determine_collinear_QCD(phasespace_set & psi);
  void calculate_collinear_QCD();
  void calculate_collinear_QCD_CS();
  void calculate_collinear_QCD_CDST();
  void determine_collinear_QEW(phasespace_set & psi);
  void calculate_collinear_QEW();
  void calculate_collinear_QEW_CS();
  void calculate_collinear_QEW_CDST();
  void output_collinear();
  void output_collinear_pdf();


  void output_CX_ncollinear_QCD();

  void determine_CX_QCD(phasespace_set & psi, int n_emission);
  void determine_CX_ncollinear_QCD(phasespace_set & psi, int n_emission);

  void determine_CX_QCD_singlet(phasespace_set & psi, int n_emission);

  void calculate_dipole_Acc_QCD(int x_a, double & Dfactor);
  void calculate_dipole_Asc_QCD(int x_a, fourvector & Vtensor, double & ME2_metric, double & ME2_vector);
  double calculate_dipole_QCD_A_ij_k(int x_a);
  double calculate_dipole_QCD_A_ij_k_massive(int x_a);
  double calculate_dipole_QCD_A_ij_a(int x_a);
  double calculate_dipole_QCD_A_ij_a_massive(int x_a);
  double calculate_dipole_QCD_A_ai_k(int x_a);
  double calculate_dipole_QCD_A_ai_b(int x_a);
  void calculate_dipole_Acc_QEW(int x_a, double & Dfactor);
  void calculate_dipole_Asc_QEW(int x_a, fourvector & Vtensor, double & ME2_metric, double & ME2_vector);
  double calculate_dipole_QEW_A_ij_k(int x_a);
  double calculate_dipole_QEW_A_ij_k_massive(int x_a);
  double calculate_dipole_QEW_A_ij_a(int x_a);
  double calculate_dipole_QEW_A_ij_a_massive(int x_a);
  double calculate_dipole_QEW_A_ai_k(int x_a);
  double calculate_dipole_QEW_A_ai_b(int x_a);






  void output_momenta(ofstream & out_comparison);
  void output_momenta_phasespace(ofstream & out_comparison, int x_a);

  void output_testpoint_collinear(phasespace_set & psi);
  void output_testpoint(phasespace_set & psi);


  void calculate_intermediate_result(phasespace_set & psi);
  void calculate_Xsection(long long i, double & Xsection, double & Xsection_delta, double & sum_weights, double & sum_weights2, double & temp_sum_weights, double & temp_sum_weights2);

  void perform_integration_step_complete(phasespace_set & psi);

  void output_zero_contribution_complete(phasespace_set & psi);

  void output_step_integration(phasespace_set & psi);
  void output_step_integration_TSV(phasespace_set & psi);
  void output_step_result(phasespace_set & psi);
  void output_step_result_TSV(phasespace_set & psi);
  void output_step_distribution_TSV(phasespace_set & psi);

  void output_step_errorplot(phasespace_set & psi);
  void output_step_execution(phasespace_set & psi);

  void determine_techcut_RA(phasespace_set & psi);


  void perform_proceeding_in(phasespace_set & psi);
  void perform_proceeding_out(phasespace_set & psi);
  void perform_proceeding_check(phasespace_set & psi);

  void old_calculate_amplitude_doublevirtual();
  void old_calculate_amplitude_doublevirtual_Vgamma(int process, int parton_identity);

  // new
  void calculate_sigma_gg( double qt2, double q2, double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_qz2, double pdf_factor_z1q, double pdf_factor_qq, int order, double A_F,  double LR, double LF, vector<double> & sigma_qT );       
                                                              
  double calculate_virtual_gg(                    double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_qz2, double pdf_factor_z1q, double pdf_factor_qq, int order, double H1_delta, double H2_delta, double LR, double LF);
                                                             

  double calculate_sigma_qqbar(double qt2, double q2, double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb, int order, double A_F, double LR, double LF, double LQ, vector<double> & sigma_qT, phasespace_set & psi);

  double virtual_qqbar(double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb, int order, double H1_delta, double H2_delta, double LR, double LF, double LQ);

  //  double_complex F2_q(double LR);
  //  void calc_H_coefficients(double q2, double s, double A0, double A1, double A2, double &H1, double &H2, int order);

  //  string run_mode;

  NJsubtraction NJ; // initialize NJsubtraction class

  string filename_result;
  vector<string> filename_moment;
  string filename_integration;
  string filename_maxevent;
  string filename_comparison;
  string filename_gnuplot;
  string filename_makegnuplot;
  string filename_proceeding;
  string filename_proceeding_2;
  string filename_execution;
  vector<string> filename_distribution;
  vector<string> filename_dddistribution;

  string filename_distribution_all_CV;

  string filename_integration_CV;
  vector<vector<string> > filename_distribution_CV;
  vector<vector<string> > filename_dddistribution_CV;
  
  vector<string> filename_integration_TSV;
  vector<string> filename_distribution_TSV;
  vector<string> filename_result_TSV;
  vector<string> filename_moment_TSV;


  int switch_qTcut;
  int n_qTcut;
  double min_qTcut;
  double step_qTcut;
  double max_qTcut;
  string binning_qTcut;
  string selection_qTcut;

  vector<double> value_qTcut;
  
  string selection_qTcut_distribution;
  string selection_no_qTcut_distribution;
  vector<int> no_qTcut_distribution;
  vector<double> value_qTcut_distribution;

  string selection_qTcut_result;
  string selection_no_qTcut_result;
  vector<int> no_qTcut_result;
  vector<double> value_qTcut_result;

  string selection_qTcut_integration;
  string selection_no_qTcut_integration;
  vector<int> no_qTcut_integration;
  vector<double> value_qTcut_integration;

  vector<int> counter_killed_qTcut;
  vector<int> counter_acc_qTcut;


  int switch_NJcut;
  int switch_NJcut_axes;
  int switch_NJcut_axes_energy;
  int switch_NJcut_measure;

   
  double scale_ren;
  double scale_fact;

  int needed_scale2_ren;
  int needed_scale2_fact;

  vector<int> n_scale_TSV;
  int max_n_scale_ren_TSV;
  int max_n_scale_fact_TSV;
  int n_dynamic_scale_TSV;
  vector<int> switch_dynamic_scale_TSV;
  vector<int> dynamic_scale_TSV;
  vector<double> relative_central_scale_TSV;
  vector<int> factor_scale_TSV;
  vector<double> central_scale_TSV;
  vector<string> name_TSV;

  int switch_TSV;
  int n_set_TSV;
  vector<string> name_set_TSV;
  vector<double> central_scale_ren_TSV;
  vector<double> central_scale_fact_TSV;
  vector<double> relative_central_scale_ren_TSV;
  vector<double> relative_central_scale_fact_TSV;
  vector<int> n_scale_ren_TSV;
  vector<int> n_scale_fact_TSV;
  vector<int> factor_scale_ren_TSV; // integer ???
  vector<int> factor_scale_fact_TSV; // integer ???
  vector<int> dynamic_scale_ren_TSV;
  vector<int> dynamic_scale_fact_TSV;
  vector<double> min_qTcut_TSV;
  vector<double> max_qTcut_TSV;

  int switch_distribution_at_all_TSV;
  vector<int> switch_distribution_TSV;
  vector<int> max_n_integrand_TSV;
  vector<double> min_qTcut_distribution_TSV;
  vector<double> max_qTcut_distribution_TSV;



  vector<int> distribution_no_qTcut;
  vector<vector<int> > distribution_no_qTcut_phasespace;

  fourvector QT_Q;
  double QT_QT;
  double QT_sqrtQ2;


  vector<int> switch_moment_TSV;
 
  vector<int> no_central_scale_ren_TSV;
  vector<int> no_central_scale_fact_TSV;

  vector<vector<double> > relative_scale_ren_TSV;
  vector<vector<double> > relative_scale_fact_TSV;

  vector<string> name_diff_set_TSV;
  int n_diff_set_TSV;
  vector<string> name_diff_set_plus_TSV;
  vector<string> name_diff_set_minus_TSV;
  vector<int> no_diff_set_plus_TSV;
  vector<int> no_diff_set_minus_TSV;

  vector<string> name_extended_set_TSV;
  int n_extended_set_TSV;

  string switch_reference;

  string name_reference_TSV;
  int no_reference_TSV;
  int no_scale_ren_reference_TSV;
  int no_scale_fact_reference_TSV;
  int no_qTcut_reference_TSV;



  int n_particle;
  int process_type;
  double E_beam;
  double E_CMS;

  int n_ps;
  int n_pc;
  vector<int> relation_pc_ps;
  vector<int> cut_ps;
  int first_non_cut_ps;
  vector<int> change_cut;

  int process_id;

  // not yet "connected"...
  vector<dipole_set> * RA_dipole;
  vector<vector<ioperator_set> > * VA_ioperator;
  vector<vector<collinear_set> > * CA_collinear;
  //  vector<vector<collinear_set> > CA_collinear;

  vector<vector<multicollinear_set> > multicollinear;
  vector<multicollinear_set> ncollinear;


  int n_pz;
  int n_moments;

  vector<int> moment_symm;

  int QCD_order;
  int active_qTcut;
  int output_n_qTcut;

  int max_dyn_ren;
  int max_dyn_fact;
  int max_scale_dyn_ren;
  int max_scale_dyn_fact;
  vector<int> n_scale_dyn_ren;
  vector<int> n_scale_dyn_fact;

  int coll_choice;
  
  int N_f;
  int N_f_active;
  int N_quarks;

  string LHAPDFname;
  int LHAPDFsubset;
  int N_nondecoupled; // should be calculated from chosen PDF set





  double alpha_S;

  double rescaling_factor_alpha_e;

  
  
  int switch_VI;
  int switch_VI_bosonic_fermionic;
  int switch_KP;
  int switch_CM;
  int switch_OL;

  int switch_yuk;
  int order_y;

  int switch_H1gg;
  int switch_H2;

  int switch_polenorm;

  int switch_old_qT_version;

  int switch_RS;
  int switch_RS_mapping;

  int switch_testcut;
  int switch_resum;

  int switch_result;
  int switch_distribution;
  int switch_moment;

  int switch_output_execution;
  int switch_output_integration;
  int switch_output_maxevent;
  int switch_output_comparison;
  int switch_output_gnuplot;
  int switch_output_proceeding;
  int switch_output_weights;
  int switch_output_result;
  int switch_output_moment;
  int switch_output_distribution;

  int switch_output_cancellation_check;
  int switch_output_testpoint;
  int switch_output_cutinfo;

  int switch_console_output_runtime;
  int switch_console_output_tau_0;
  int switch_console_output_techcut_RA;
  int switch_console_output_phasespace_issue;
  int switch_console_output_ME2_issue;



  vector<vector<int> > no_value_ren_TSV;
  vector<vector<int> > no_value_fact_TSV;

  vector<int> value_no_central_scale_fact; // where is 'central' scale in each DS set
  vector<double> value_central_scale_fact; // where is 'central' scale in each DS set
  vector<double> value_central_logscale2_fact; // where is 'central' scale in each DS set

  vector<vector<double> > value_relative_scale_ren;
  vector<vector<double> > value_relative_scale2_ren;
  vector<vector<vector<double> > > value_alpha_S_TSV;
  vector<vector<vector<double> > > value_relative_factor_alpha_S;
  vector<vector<vector<double> > > value_scale_ren;
  vector<vector<vector<double> > > value_scale2_ren;

  vector<vector<double> > value_relative_scale_fact;
  vector<vector<double> > value_relative_scale2_fact;
  vector<vector<double> > value_relative_logscale2_fact;
  vector<vector<vector<double> > > value_scale_fact;
  vector<vector<vector<double> > > value_scale2_fact;

  vector<vector<vector<vector<vector<double> > > > > value_pdf_factor;
  // for debugging against Sherpa
  vector<vector<vector<vector<vector<double> > > > > value_pdf_factor_1;
  vector<vector<vector<vector<vector<double> > > > > value_pdf_factor_2;
  vector<vector<vector<vector<vector<vector<double> > > > > > value_pdf_factor_combination_1;
  vector<vector<vector<vector<vector<vector<double> > > > > > value_pdf_factor_combination_2;
  // end

  vector<vector<vector<vector<vector<double> > > > > value_list_pdf_factor_TSV;
  //  vector<vector<vector<double> > > CX_value_integrand_TSV_qTcut;
  ///  vector<vector<vector<vector<double> > > > value_integrand_qTcut_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > value_integrand_qTcut_TSV;
  vector<vector<vector<vector<double> > > > value_ren_integrand_qTcut_TSV;
  vector<vector<vector<vector<double> > > > value_fact_integrand_qTcut_TSV;
  //  vector<vector<vector<double> > > integrand_TSV_qTcut;
  vector<vector<vector<vector<double> > > > integrand_qTcut_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > ps_integrand_qTcut_TSV;




  vector<vector<double*> > pointer_relative_scale_ren;
  vector<vector<double*> > pointer_relative_scale2_ren;
  vector<vector<vector<double*> > > pointer_alpha_S_TSV;
  vector<vector<vector<double*> > > pointer_relative_factor_alpha_S;
  vector<vector<vector<double*> > > pointer_scale_ren;
  vector<vector<vector<double*> > > pointer_scale2_ren;

  vector<vector<double*> > pointer_relative_scale_fact;
  vector<vector<double*> > pointer_relative_scale2_fact;
  vector<vector<double*> > pointer_relative_logscale2_fact;
  vector<vector<vector<double*> > > pointer_scale_fact;
  vector<vector<vector<double*> > > pointer_scale2_fact;

  vector<vector<vector<vector<vector<double*> > > > > pointer_pdf_factor;
  // for debugging against Sherpa
  vector<vector<vector<vector<vector<double*> > > > > pointer_pdf_factor_1;
  vector<vector<vector<vector<vector<double*> > > > > pointer_pdf_factor_2;
  vector<vector<vector<vector<vector<vector<double*> > > > > > pointer_pdf_factor_combination_1;
  vector<vector<vector<vector<vector<vector<double*> > > > > > pointer_pdf_factor_combination_2;
  // end

  vector<vector<vector<vector<vector<double*> > > > > pointer_list_pdf_factor_TSV;

  vector<double> value_ME2term;
  vector<vector<vector<double> > > value_ME2term_ren;
  vector<vector<vector<vector<double> > > > value_ME2term_fact;

  vector<vector<vector<vector<vector<double*> > > > > pointer_ME2term;



  vector<vector<vector<vector<double> > > > value_logscale2_fact_papi; // ???
  vector<vector<vector<double> > > data_K;
  vector<vector<vector<vector<vector<double> > > > > value_data_P;
  vector<vector<vector<vector<double> > > > value_ME2_KP;


  vector<vector<vector<vector<vector<double> > > > > ps_integrand_TSV;
  vector<vector<vector<double> > > integrand_TSV;
  vector<vector<vector<double> > > sum_weight_TSV;
  vector<vector<vector<double> > > sum_weight2_TSV;
  vector<vector<vector<double> > > fullsum_weight_TSV;
  vector<vector<vector<double> > > fullsum_weight2_TSV;
  vector<vector<vector<vector<double> > > > sum_weight_qTcut_TSV;
  vector<vector<vector<vector<double> > > > sum_weight2_qTcut_TSV;
  vector<vector<vector<vector<double> > > > fullsum_weight_qTcut_TSV;
  vector<vector<vector<vector<double> > > > fullsum_weight2_qTcut_TSV;

  vector<vector<vector<vector<vector<vector<double> > > > > > ps_moment_TSV;
  vector<vector<vector<vector<double> > > > moment_TSV;
  vector<vector<vector<vector<double> > > > sum_moment_TSV;
  vector<vector<vector<vector<double> > > > sum_moment2_TSV;
  vector<vector<vector<vector<double> > > > fullsum_moment_TSV;
  vector<vector<vector<vector<double> > > > fullsum_moment2_TSV;
  vector<vector<vector<vector<vector<double> > > > > sum_moment_qTcut_TSV;
  vector<vector<vector<vector<vector<double> > > > > sum_moment2_qTcut_TSV;
  vector<vector<vector<vector<vector<double> > > > > fullsum_moment_qTcut_TSV;
  vector<vector<vector<vector<vector<double> > > > > fullsum_moment2_qTcut_TSV;

  

  vector<vector<vector<double> > > bin_count_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > bin_weight_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > bin_weight2_TSV;

  vector<vector<vector<double> > > change_weight_TSV;
  vector<vector<vector<double> > > change_weight2_TSV;
  vector<vector<vector<double> > > change_weight_2nd_TSV;
  vector<vector<vector<double> > > change_weight2_2nd_TSV;
  /*
  vector<vector<vector<vector<double> > > > change_bin_weight_TSV;
  vector<vector<vector<vector<double> > > > change_bin_weight2_TSV;
  vector<vector<vector<vector<double> > > > change_bin_weight_2nd_TSV;
  vector<vector<vector<vector<double> > > > change_bin_weight2_2nd_TSV;
  */
  vector<vector<vector<vector<vector<double> > > > > change_qTcut_bin_weight_TSV;
  vector<vector<vector<vector<vector<double> > > > > change_qTcut_bin_weight2_TSV;
  vector<vector<vector<vector<vector<double> > > > > change_qTcut_bin_weight_2nd_TSV;
  vector<vector<vector<vector<vector<double> > > > > change_qTcut_bin_weight2_2nd_TSV;


  vector<int> no_min_qTcut_TSV;
  vector<int> no_max_qTcut_TSV;
  vector<int> no_min_qTcut_distribution_TSV;
  vector<int> no_max_qTcut_distribution_TSV;
  



  map<string, int> access_object;

  event_set esi;

  int runtime_jet_algorithm; // ???
  vector<int> runtime_jet_recombination;
  vector<int> runtime_photon_recombination;
  vector<int> runtime_photon_isolation;
  vector<int> runtime_original;
  vector<int> runtime_missing;
  vector<int> runtime_object;
  vector<vector<int> > runtime_order;
  vector<vector<int> > runtime_order_inverse;

  vector<vector<int> > ps_n_partonlevel;
  vector<vector<int> > ps_relevant_n_partonlevel;

  vector<vector<int> > ps_runtime_jet_algorithm;
  vector<vector<int> > ps_runtime_photon;
  vector<vector<int> > ps_runtime_photon_recombination;
  vector<vector<int> > ps_runtime_photon_isolation;

  vector<vector<int> > ps_runtime_jet_recombination;

  vector<vector<int> > ps_runtime_original;
  vector<vector<int> > ps_runtime_missing;
  vector<vector<int> > ps_runtime_object;
  vector<vector<vector<int> > > ps_runtime_order;
  vector<vector<vector<int> > > ps_runtime_order_inverse;

  map<string, string> equivalent_object;
  map<int, int> equivalent_no_object;
  map<string, int> no_relevant_object;






  int frixione_isolation;
  double frixione_n;
  double frixione_epsilon;
  double frixione_fixed_ET_max;
  double frixione_delta_0;
  int frixione_jet_removal;

  int photon_recombination;
  int photon_R_definition;
  double photon_R;
  double photon_R2;
  double photon_E_threshold_ratio;
  int photon_jet_algorithm;

  int photon_photon_recombination;
  double photon_photon_recombination_R;

  int jet_algorithm;
  int jet_R_definition;
  double jet_R;
  double jet_R2;
  double parton_y_max;
  double parton_eta_max;


  string name_process;

  //  vector<vector<int> > type_parton;
  vector<vector<double> > mass_parton;
  vector<vector<double> > mass2_parton;

  vector<double> M;
  vector<double> M2;

  int massive_QCD;
  int massive_QEW;

  vector<vector<fourvector> > p_parton;
  vector<vector<vector<particle> > > particle_event;
  vector<vector<int> > n_object;

  //  vector<double> track_particle_x;

  vector<vector<fourvector> > start_p_parton;
  vector<vector<vector<particle> > > start_particle_event;
  vector<vector<int> > start_n_object;

  vector<int> recombination_history;

  // overlap with phasespace_set:
  vector<double> x_pdf;
  double boost;
  vector<double> z_coll;


  // not yet initialised !!!  
  vector<vector<int> > bin;
  vector<vector<int> > bin_max;

  vector<vector<long long> > bin_counts;
  
  // rename at some point:
  vector<xdistribution> dat;
  vector<dddistribution> dddat;
  vector<xdistribution> extended_distribution;

  vector<double> fakeasymfactor; // can most likely be removed !!!
  vector<vector<double> > bin_weight;
  vector<vector<double> > bin_weight2;
  vector<vector<vector<double> > > bin_weight_CV;
  vector<vector<vector<double> > > bin_weight2_CV;





  // reorganize completely:
  vector<vector<double> > moment;
  vector<vector<vector<double> > > directed_moment;

  double max_integrand;
  double sigma_normalization;
  double sigma_normalization_deviation;



  int switch_CV;
  int variation_mu_fact_CV;
  int variation_mu_ren_CV;
  int n_scales_CV;
  int variation_factor_CV;
  int no_central_scale_CV;
  double central_scale_CV;

  double prefactor_reference;
  double prefactor_CV;

  int dynamic_scale;
  int dynamic_scale_CV;


  vector<double> mu_fact;
  vector<double> mu_ren;

  vector<double> scale_ren_CV;
  vector<double> scale_fact_CV;
  vector<double> rel_scale_factor_CV;

  vector<vector<double> > mu_fact_CV;
  vector<vector<double> > mu_ren_CV;
  vector<vector<double> > alpha_S_CV;
  vector<vector<double> > rel_alpha_S_CV;

  vector<double> rel_scale_factor_ren_CV;
  vector<double> rel_scale_factor_fact_CV;

  vector<string> directory_name_scale_CV;


  double integrand;
  vector<double> integrand_CV;
  vector<vector<double> > integrand_qTcut_CV;

  vector<vector<double> > integrand_D;
  vector<vector<vector<double> > > integrand_D_CV;
  vector<vector<vector<vector<double> > > > integrand_D_qTcut_CV;

  double this_psp_weight;
  double this_psp_weight2;
  vector<vector<double> > this_psp_weight_CV;
  vector<vector<double> > this_psp_weight2_CV;
  double step_sum_weight;
  double step_sum_weight2;
  vector<vector<double> > step_sum_weight_CV;
  vector<vector<double> > step_sum_weight2_CV;
  double full_sum_weight;
  double full_sum_weight2;
  vector<vector<double> > full_sum_weight_CV;
  vector<vector<double> > full_sum_weight2_CV;


  vector<double> this_psp_moment;
  vector<double> this_psp_moment2;
  vector<double> step_sum_moment;
  vector<double> step_sum_moment2;
  vector<double> full_sum_moment; 
  vector<double> full_sum_moment2;

  vector<vector<vector<double> > > this_psp_moment_CV;
  vector<vector<vector<double> > > this_psp_moment2_CV;
  vector<vector<vector<double> > > step_sum_moment_CV;
  vector<vector<vector<double> > > step_sum_moment2_CV;
  vector<vector<vector<double> > > full_sum_moment_CV;
  vector<vector<vector<double> > > full_sum_moment2_CV;

  
  // only needed for RA contributions... really needed ???
  vector<vector<double> > var_A_integrand_cut;
  vector<vector<double> > var_R_integrand_cut;
  vector<vector<double> > var_R_integrand_cut_incl;

  vector<vector<vector<double> > > A_integrand_moment_cut;
  vector<vector<vector<double> > > R_integrand_moment_cut;
  vector<vector<vector<double> > > R_integrand_moment_cut_incl;


  int RA_x_a;
  double RA_techcut_integrand;








  int id_scales;

  int map_value_scale_fact;
  vector<int> map_value_scale_fact_CV;
  vector<double> value_mu_fact_central;
  vector<vector<double> > value_mu_fact_rel;
  vector<vector<double> > value_mu_fact;

  int map_value_scale_ren;
  vector<int> map_value_scale_ren_CV;
  vector<double> value_mu_ren_central;
  vector<vector<double> > value_mu_ren_rel;
  vector<vector<double> > value_mu_ren;
  vector<vector<double> > value_alpha_S;
  vector<vector<double> > value_factor_alpha_S;

// should be removed when ...value.. is completely installed ???

  double ME2;
  int check_vanishing_ME2_end;
  int flag_vanishing_ME2;
  int n_event_vanishing_ME2;

  double mu_central;
  double var_mu_ren;
  double var_mu_fact;
  double var_alpha_S_reference;
  double var_rel_alpha_S;
  vector<double> var_mu_ren_CV;
  vector<double> var_mu_fact_CV;
  vector<double> var_alpha_S_CV;
  vector<double> var_rel_alpha_S_CV;






  // constants relevant for qT subtraction

  //  int n_jet_born; // -> csi
  //  int n_particle_born; // -> csi
  vector<fourvector> NJ_q_axes;
  vector<fourvector> NJ_n_axes;
  vector<double> NJ_Ei;
  vector<fourvector> NJ_q_axes_frame;
  vector<fourvector> NJ_n_axes_frame;
  vector<double> NJ_Qi;
  vector<fourvector> NJ_p_parton_frame;





  // constants relevant for qT subtraction

  int order_alphas_born;
  int initial_channel;

  int initial_gg;
  int initial_qqx;

  int initial_pdf_gg;
  int initial_pdf_qqx;
  int initial_pdf_gq;
  int initial_pdf_rest;

  int initial_pdf_diag;

  int initial_diag;
  int initial_diag_gg;
  int initial_diag_qqx;

  int QT_initialstate_type;
  int QT_finalstate_massive_coloured;
  vector<correlationoperator> QT_correlationoperator;
  vector<double> QT_ME2_cf;

  vector<double> QT_ME2_loopcf;

  vector<double> QT_ME2_Born4cf;

  vector<double> I1_int;
  vector<double> I2_int;
  vector<double> I3_int;
  vector<double> I4_int;

  vector<double> LL0_int;
  vector<double> LL1_int;
  vector<double> LL2_int;
  vector<double> LL3_int;

  double beta0;
  double beta1;
  double Kappa;

  double A1g;
  double B1g;
  double A2g;
  double Delta2gg;
  double B2g;

  double D0gggg;
  double D1gggg;
  double Deltagggg;

  double H2ggD0;

  double A1q;
  double B1q;
  double A2q;
  double Delta2qq;
  double B2q;



  
  double gamma2Q;
  double gammacusp2;


  double q2;
  double v34;
  double log_v34;
  double dilogpT2m2;
  double logmT2m2;
  double L34;
  double gammacusp2_v;

  //  CC    Coefficients of D0 and D1 in P*P (as/pi normalization)
  double D0qqqq;
  double D1qqqq;
  //CC    Coefficients of delta(1-z) in P*P
  double Deltaqqqq;

  double H2qqD0;

  double A_F;
  double B2_A_F;


  
  double Ft1born;
  vector<double> Ft1born_4correlator;
  double Ft1loop;

  double Gamma1;
  double Gamma1born;
  double Gamma1loop;

  double Gamma2;
  double Gamma1squared;
  double commutator_Gammat1_Ft1;
  double anticommutator_Gammat1_Ft1;
  double Gammat1_Ft1;
  double Ft1_Gammat1;

  double Ft1;

  vector<vector<double> > coll_tH1F_contribution;
  vector<vector<double> > coll_tH1_contribution;
  vector<vector<double> > coll_tH1_only_H1_delta_contribution;
  vector<vector<double> > coll_tH1_without_H1_delta_contribution;
  vector<vector<double> > coll_tgaga_contribution;
  vector<vector<double> > coll_tcga_contribution;
  vector<vector<double> > coll_tgamma2_contribution;
  vector<vector<double> > coll_tH2_contribution;

  vector<double> coll_tH1F;
  vector<double> coll_tH1;
  vector<double> coll_tH1_only_H1_delta;
  vector<double> coll_tH1_without_H1_delta;
  vector<double> coll_tgaga;
  vector<double> coll_tcga;
  vector<double> coll_tgamma2;
  vector<double> coll_tH2;

  // central - CV - D
  vector<vector<double> > value_LR;
  vector<vector<double> > value_LF;
  vector<vector<double> > value_LQ;
 
  // central - CV
  vector<vector<vector<double > > > coll_tH1F_pdf;
  vector<vector<vector<double > > > coll_tH1_pdf;
  vector<vector<vector<double > > > coll_tH1_only_H1_delta_pdf;
  vector<vector<vector<double > > > coll_tH1_without_H1_delta_pdf;
  vector<vector<vector<double > > > coll_tgaga_pdf;
  vector<vector<vector<double > > > coll_tcga_pdf;
  vector<vector<vector<double > > > coll_tgamma2_pdf;
  vector<vector<vector<double > > > coll_tH2_pdf;

  vector<vector<double> > value_sig11;
  vector<vector<double> > value_tH1F;
  vector<vector<double> > value_tH1;
  vector<vector<double> > value_tH1_only_H1_delta;
  vector<vector<double> > value_tH1_without_H1_delta;

  vector<vector<double> > value_tgaga;
  vector<vector<double> > value_tcga;
  vector<vector<double> > value_tgamma2;

  vector<vector<double> > value_sig12;
  vector<vector<double> > value_sig23;
  vector<vector<double> > value_sig24;
  vector<vector<vector<vector<double> > > > value_sig21;
  vector<vector<vector<vector<double> > > > value_sig22;

  vector<vector<double> > value_tH2;

  vector<vector<vector<vector<double> > > > value_H1full;

  // D
  vector<vector<vector<vector<double> > > > coll_tH1F_pdf_D;
  vector<vector<vector<vector<double> > > > coll_tH1_pdf_D;
  vector<vector<vector<vector<double> > > > coll_tH1_only_H1_delta_pdf_D;
  vector<vector<vector<vector<double> > > > coll_tH1_without_H1_delta_pdf_D;
  vector<vector<vector<vector<double> > > > coll_tgaga_pdf_D;
  vector<vector<vector<vector<double> > > > coll_tcga_pdf_D;
  vector<vector<vector<vector<double> > > > coll_tgamma2_pdf_D;
  vector<vector<vector<vector<double> > > > coll_tH2_pdf_D;

  vector<vector<vector<double> > > value_sig11_D;
  vector<vector<vector<double> > > value_tH1F_D;
  vector<vector<vector<double> > > value_tH1_D;
  vector<vector<vector<double> > > value_tH1_only_H1_delta_D;
  vector<vector<vector<double> > > value_tH1_without_H1_delta_D;

  vector<vector<vector<double> > > value_tgaga_D;
  vector<vector<vector<double> > > value_tcga_D;
  vector<vector<vector<double> > > value_tgamma2_D;
 
  vector<vector<vector<double> > > value_sig12_D;
  vector<vector<vector<double> > > value_sig23_D;
  vector<vector<vector<double> > > value_sig24_D;
  vector<vector<vector<vector<vector<double> > > > > value_sig21_D;
  vector<vector<vector<vector<vector<double> > > > > value_sig22_D;

  vector<vector<vector<double> > > value_tH2_D;

  vector<vector<vector<vector<vector<double> > > > > value_H1full_D;

  // TSV
  vector<vector<vector<vector<double> > > > coll_tH1F_pdf_TSV;
  vector<vector<vector<vector<double> > > > coll_tH1_pdf_TSV;
  vector<vector<vector<vector<double> > > > coll_tH1_only_H1_delta_pdf_TSV;
  vector<vector<vector<vector<double> > > > coll_tH1_without_H1_delta_pdf_TSV;
  vector<vector<vector<vector<double> > > > coll_tgaga_pdf_TSV;
  vector<vector<vector<vector<double> > > > coll_tcga_pdf_TSV;
  vector<vector<vector<vector<double> > > > coll_tgamma2_pdf_TSV;
  vector<vector<vector<vector<double> > > > coll_tH2_pdf_TSV;

  vector<vector<double> > value_LR_TSV;
  vector<vector<double> > value_LF_TSV;
  vector<vector<double> > value_LQ_TSV;

  vector<vector<vector<double> > > value_sig11_TSV;
  vector<vector<vector<double> > > value_tH1F_TSV;
  vector<vector<vector<double> > > value_tH1_TSV; 
  vector<vector<vector<double> > > value_tH1_only_H1_delta_TSV; 
  vector<vector<vector<double> > > value_tH1_without_H1_delta_TSV; 
  
  vector<vector<vector<double> > > value_tgaga_TSV;
  vector<vector<vector<double> > > value_tcga_TSV;
  vector<vector<vector<double> > > value_tgamma2_TSV;
  
  vector<vector<vector<vector<vector<double> > > > > value_H1full_TSV; 
  
  vector<vector<vector<double> > > value_sig12_TSV;
  vector<vector<vector<double> > > value_sig23_TSV; 
  vector<vector<vector<double> > > value_sig24_TSV;
  vector<vector<vector<vector<vector<double> > > > > value_sig21_TSV;
  vector<vector<vector<vector<vector<double> > > > > value_sig22_TSV;


  vector<vector<vector<double> > > value_tH2_TSV;


  // used for qT subtraction 
  vector<double> pdf_factor;
  vector<double> QT_pdf_factor_z1x2;
  vector<double> QT_pdf_factor_x1z2;
  vector<double> QT_pdf_factor_gx2;
  vector<double> QT_pdf_factor_x1g;

  vector<double> QT_pdf_factor_gg;
  vector<double> QT_pdf_factor_qx2;
  vector<double> QT_pdf_factor_x1q;
  vector<double> QT_pdf_factor_z1z2;
  vector<double> QT_pdf_factor_gz2;
  vector<double> QT_pdf_factor_z1g;
  vector<double> QT_pdf_factor_qbx2;
  vector<double> QT_pdf_factor_x1qb;
  // new
  vector<double> QT_pdf_factor_qq;
  vector<double> QT_pdf_factor_qz2;
  vector<double> QT_pdf_factor_z1q;


  vector<vector<double> > pdf_factor_CV;
  vector<vector<double> > QT_pdf_factor_z1x2_CV;
  vector<vector<double> > QT_pdf_factor_x1z2_CV;
  vector<vector<double> > QT_pdf_factor_gx2_CV;
  vector<vector<double> > QT_pdf_factor_x1g_CV;
  vector<vector<double> > QT_pdf_factor_gg_CV;
  vector<vector<double> > QT_pdf_factor_qx2_CV;
  vector<vector<double> > QT_pdf_factor_x1q_CV;
  vector<vector<double> > QT_pdf_factor_z1z2_CV;
  vector<vector<double> > QT_pdf_factor_gz2_CV;
  vector<vector<double> > QT_pdf_factor_z1g_CV;
  vector<vector<double> > QT_pdf_factor_qbx2_CV;
  vector<vector<double> > QT_pdf_factor_x1qb_CV;
  // new
  vector<vector<double> > QT_pdf_factor_qq_CV;
  vector<vector<double> > QT_pdf_factor_qz2_CV;
  vector<vector<double> > QT_pdf_factor_z1q_CV;
  
  vector<vector<vector<double> > > QT_value_pdf_factor;
  vector<vector<vector<double> > > QT_value_pdf_factor_z1x2;
  vector<vector<vector<double> > > QT_value_pdf_factor_x1z2;
  vector<vector<vector<double> > > QT_value_pdf_factor_gx2;
  vector<vector<vector<double> > > QT_value_pdf_factor_x1g;
  vector<vector<vector<double> > > QT_value_pdf_factor_gg;
  vector<vector<vector<double> > > QT_value_pdf_factor_qx2;
  vector<vector<vector<double> > > QT_value_pdf_factor_x1q;
  vector<vector<vector<double> > > QT_value_pdf_factor_z1z2;
  vector<vector<vector<double> > > QT_value_pdf_factor_gz2;
  vector<vector<vector<double> > > QT_value_pdf_factor_z1g;
  vector<vector<vector<double> > > QT_value_pdf_factor_qbx2;
  vector<vector<vector<double> > > QT_value_pdf_factor_x1qb;
  // new
  vector<vector<vector<double> > > QT_value_pdf_factor_qq;
  vector<vector<vector<double> > > QT_value_pdf_factor_qz2;
  vector<vector<vector<double> > > QT_value_pdf_factor_z1q;

  //  vector<vector<vector<double> > > QT_value_sigma_qT;
  vector<vector<vector<double> > > QT_value_integrand_qTcut;
  vector<vector<vector<vector<double> > > > QT_value_integrand_qTcut_D;

  vector<double> QT_sigma_qT;
  vector<vector<double> > QT_sigma_qT_CV;
  
  vector<double> QT_sig_CV;
  vector<double> QT_virt_CV;

  double QT_A0;
  double QT_A1;
  double QT_A2;
  double QT_H1_delta;
  double QT_H2_delta;

  double QT_H0_doublespinflip;
  double QT_H0_DG;
  
  double QT_Qres;
  bool dynamical_Qres;
  double QT_Qres_prefactor;
  
  bool switch_resummation; 

  vector<double> RA_mu_ren;
  vector<double> RA_mu_fact;
  vector<double> RA_alpha_S_reference;
  vector<double> RA_rel_alpha_S;
  vector<vector<double> > RA_mu_ren_CV;
  vector<vector<double> > RA_mu_fact_CV;
  vector<vector<double> > RA_alpha_S_CV;
  vector<vector<double> > RA_rel_alpha_S_CV;
  vector<vector<double> > RA_pdf_factor;
  vector<vector<vector<double> > > RA_pdf_factor_CV;

  vector<vector<vector<double> > > RA_value_mu_fact;
  vector<vector<vector<double> > > RA_value_mu_ren;
  vector<vector<vector<double> > > RA_value_alpha_S;
  vector<vector<vector<double> > > RA_value_factor_alpha_S;

  vector<double> RA_ME2;
  vector<double> var_RA_ME2;
  vector<vector<double> > var_RA_ME2_CV;
  vector<vector<double> > RA_ME2_moment;
  vector<vector<vector<double> > > RA_ME2_moment_CV;

  vector<double> RA_integrand_moment;
  vector<vector<double> > RA_integrand_moment_CV;



  double VA_b_ME2;
  double VA_V_ME2;
  double VA_X_ME2;
  double VA_I_ME2;
  vector<double> VA_X_ME2_CV;
  vector<vector<double> > VA_X_ME2_vr_mr;

  double VA_DeltaUV;
  double VA_DeltaIR1;
  double VA_DeltaIR2;
  int VA_delta_flag;
  
  vector<vector<double> > VA_ME2_cf;
  vector<vector<double> > VA_I_ME2_cf;
  
  vector<vector<int> > CA_dipole_splitting;
  vector<vector<double> > CA_ME2_cf;

  vector<vector<double> > CA_value_log_mu2_fact;
  vector<vector<vector<vector<double> > > > CA_value_ME2_KP;
  vector<vector<vector<vector<vector<double> > > > > CA_value_pdf_factor;
  vector<vector<vector<double> > > CA_value_integrand;
  vector<vector<vector<vector<double> > > > CA_value_integrand_D;
  vector<vector<double> > CA_sum_value_integrand;
  vector<vector<vector<double> > > CA_sum_value_integrand_D;



  vector<vector<double> > xz_pdf;
  vector<vector<double> > xz_coll;
  vector<vector<double> > xz_factor;

  vector<vector<vector<vector<vector<double> > > > > CX_value_pdf_factor;
  vector<vector<vector<vector<vector<int> > > > > CX_combination_pdf;

  vector<vector<vector<vector<double> > > > value_list_pdf_factor;
  vector<vector<vector<vector<int> > > > list_combination_pdf;
  vector<vector<int> > list_combination_pdf_emission;


  vector<vector<vector<double> > > CX_value_integrand;
  vector<vector<vector<vector<double> > > > CX_value_integrand_D;
  vector<vector<double> > CX_sum_value_integrand;
  vector<vector<vector<double> > > CX_sum_value_integrand_D;
  vector<vector<vector<double> > > CX_value_integrand_qTcut;
  vector<vector<vector<vector<double> > > > CX_value_integrand_qTcut_D;

  vector<vector<vector<vector<vector<double> > > > > CX_value_integrand_RF;
  vector<vector<vector<vector<vector<vector<double> > > > > > CX_value_integrand_RF_D;
  vector<vector<vector<vector<double> > > > CX_sum_value_integrand_RF;
  vector<vector<vector<vector<vector<double> > > > > CX_sum_value_integrand_RF_D;
  vector<vector<vector<vector<vector<double> > > > > CX_value_integrand_RF_qTcut;
  vector<vector<vector<vector<vector<vector<double> > > > > > CX_value_integrand_RF_qTcut_D;

  munich *xmunich;

  contribution_set *csi;
  phasespace_set *psi;

  string process_class;
  string subprocess;
  string type_perturbative_order;
  string type_contribution;
  string type_correction;
  int contribution_order_alpha_s;
  int contribution_order_alpha_e;
  int contribution_order_interference;


  int int_end;
  
  clock_t start;
  clock_t end;
  int h;
  int min;
  int sec;
  int time_counter;
  int sec_import;

  model_set msi;

  //  parameter_model model;

  vector<string> OL_parameter;
  vector<string> OL_value;

  user_defined user;

  string path_to_main;


  double Xsection;
  double Xsection_delta;
  vector<vector<double> > Xsection_CV;
  vector<vector<double> > Xsection_delta_CV;
  vector<vector<vector<vector<double> > > > Xsection_TSV;
  vector<vector<vector<vector<double> > > > Xsection_delta_TSV;

  int temp_n_step;

  vector<vector<int> > combination_pdf;
  vector<vector<vector<vector<int> > > > CA_combination_pdf;
  //  vector<vector<double> > CA_QijQk_Q2ij;
  //  vector<vector<vector<double> > > CA_QijQk_Q2ij;
  vector<vector<vector<double> > > CA_Q2f;
  //  vector<vector<double> > CA_Q2f;
  //  vector<vector<vector<double> > > CA_Qk_Qij;

  //  string pdf_selection;
  //  string pdf_disable;
  vector<string> pdf_selection;
  vector<string> pdf_disable;
  //  int pdf_content_modify;
  vector<vector<vector<int> > > allowed_all_pdf;

  vector<string> jet_algorithm_selection;
  vector<string> jet_algorithm_disable;
  vector<int> jet_algorithm_list;

  vector<string> photon_recombination_selection;
  vector<string> photon_recombination_disable;
  vector<int> photon_recombination_list;

  string unit_calculation;
  string unit_result;
  string unit_distribution;

  double unit_factor_calculation;
  double unit_factor2_calculation;
  double unit_factor_result;
  double unit_factor_distribution;

  friend std::ostream & operator << (std::ostream &, const observable_set &);
};

#endif
