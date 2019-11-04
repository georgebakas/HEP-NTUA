#ifndef INPUTPARAMETERSET_H
#define INPUTPARAMETERSET_H

class inputparameter_set{
 private:

 public:
////////////////////
//  constructors  //
////////////////////
  inputparameter_set();
  inputparameter_set(string basic_process_class, string & subprocess);//, vector<string> & readin);

  void set_default_geneva(string run_mode);

  string output_level;

  //  observable_set oset;
  user_defined user;
  contribution_set csi;
  event_set esi;

  //////////////////////////
  //  unknown parameters  //
  //////////////////////////

  string path_MUNICH;

  string path_to_main;
  // ???

  int process_number;
  // ???

  int n_particle;
  // ??? needed here ???

  int test_output;
  // ???

  double tau_0;
  // should be evaluated later !!! actually no input !!!

  int ckm_choice;
  // not relevant by now - but maybe later...

  int select_contribution;
  // ???

  int n_moments;
  // ??? should be replaced by input file !!!

  string dir_directory;
  // ???

  //  vector<steering_optimization> steering;
  // !!! removed

  int no_contribution;
  // ??? meaning ???

  //  string filename_process;
  // !!! removed

  vector<string> impulsname;
  // ??? remove

  //  int n_madgraph;
  // removed !!!

  vector<string> name_fiducial_cut;

  //////////////////////////////////////////////////
  //  parameters determined during input read-in  //
  //////////////////////////////////////////////////

  int max_perturbative_QCD_order;
  int present_type_perturbative_order;
  int type_perturbative_order_counter;
  vector<string> collection_type_perturbative_order;
  int present_type_contribution;
  int type_contribution_counter;
  vector<string> collection_type_contribution;
  int present_type_correction;
  int type_correction_counter;
  vector<string> collection_type_correction;

  //////////////////////////////////////////////////////////////
  //  selection of process and contribution to be calculated  //
  //////////////////////////////////////////////////////////////

  string process_class;
  vector<string> decay;
  string type_perturbative_order;
  string type_contribution;
  string type_correction;
  int contribution_order_alpha_s;
  int contribution_order_alpha_e;
  int contribution_order_interference;
  // -> csi


  ////////////////////////////////////////////////
  //  switches to steer calculation of results  //
  ////////////////////////////////////////////////

  int switch_result;
  int switch_distribution;
  int switch_moment;

  ///////////////////////////////////////////////
  //  switches to steer output of calculation  //
  ///////////////////////////////////////////////

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

  int switch_testcut;  // ???

  int switch_console_output_runtime;
  int switch_console_output_tau_0;
  int switch_console_output_techcut_RA;
  int switch_console_output_phasespace_issue;
  int switch_console_output_ME2_issue;


  ////////////////////////////////////////////////////////////////////////
  //  unit of calculation output / result output / distribution output  //
  ////////////////////////////////////////////////////////////////////////

  string unit_calculation;
  string unit_result;
  string unit_distribution;

  ///////////////////////
  //  beam parameters  //
  ///////////////////////

  double E;
  int coll_choice;
  int pdf_set;
  // ??? only via LHAPDF ???

  //  int pdf_content_modify;
  // ???

  vector<string> pdf_selection;
  vector<string> pdf_disable;
  string LHAPDFname;
  int LHAPDFsubset;
  vector<string> contribution_LHAPDFname;
  vector<int> contribution_LHAPDFsubset;

  ////////////////////////////////////////
  //  jet and jet-algorithm parameters  //
  ////////////////////////////////////////

  int jet_algorithm;
  int jet_R_definition;
  double jet_R;
  double parton_y_max;
  double parton_eta_max;
  vector<string> jet_algorithm_selection;
  vector<string> jet_algorithm_disable;
  int N_f;
  int N_f_active;
  int N_quarks;
  int N_nondecoupled; 
  // !!! should be calculated from chosen PDF set

  ///////////////////////////////////////
  //  photon-recombination parameters  //
  ///////////////////////////////////////

  int photon_recombination;
  int photon_R_definition;
  double photon_R;
  double photon_E_threshold_ratio;
  int photon_jet_algorithm;
  vector<string> photon_recombination_selection;
  vector<string> photon_recombination_disable;
  int photon_photon_recombination;
  double photon_photon_recombination_R;

  ///////////////////////////////////
  //  photon-isolation parameters  //
  ///////////////////////////////////

  int frixione_isolation;
  double frixione_n;
  double frixione_epsilon;
  double frixione_fixed_ET_max;
  double frixione_delta_0;
  int frixione_jet_removal;

  /////////////////////////////////
  //  qT-subtraction parameters  //
  /////////////////////////////////

  int switch_qTcut;
  int n_qTcut;
  double min_qTcut;
  double step_qTcut;
  double max_qTcut;
  string binning_qTcut;
  string selection_qTcut;

  string selection_qTcut_distribution;
  string selection_no_qTcut_distribution;
  string selection_qTcut_result;
  string selection_no_qTcut_result;
  string selection_qTcut_integration;
  string selection_no_qTcut_integration;
  //  vector<double> value_qTcut;
  // !!! no input parameter !!!
  //  vector<int> no_qTcut_distribution;
  // !!! no input parameter !!!
  //  vector<double> value_qTcut_distribution;
  // !!! no input parameter !!!

  int switch_NJcut;
  int switch_NJcut_axes;
  int switch_NJcut_axes_energy;
  int switch_NJcut_measure;


  ////////////////////////////////////////////////////
  //  technical switches for selected contributions //
  ////////////////////////////////////////////////////

  int switch_KP;
  int switch_VI;
  int switch_VI_bosonic_fermionic;
  // ???

  int switch_H1gg;
  int switch_H2;

  int switch_polenorm; 
  // ??? should actually steer the OL settings on pole normalization
  int switch_old_qT_version; // 0 = new, 1 = old 

  int switch_RS;
  int switch_RS_mapping;

  int switch_CM;
  // ???

  int switch_OL;
  // ??? 


  int switch_yuk;
  int order_y;


  ///////////////////////////////
  //  scale-choice parameters  //
  ///////////////////////////////

  double scale_fact;
  double scale_ren;
  int dynamic_scale;
  double prefactor_reference;

  ///////////////////////////////////////
  //  scale variation parameters - CV  //
  ///////////////////////////////////////

  int switch_CV;
  int n_scales_CV;
  int dynamic_scale_CV;
  int variation_mu_ren_CV;
  // ???

  int variation_mu_fact_CV;
  // ???

  int variation_factor_CV;
  int central_scale_CV;
  // ???

  double prefactor_CV;


  //////////////////////////////////////
  //  weight-optimization parameters  //
  //////////////////////////////////////

  int switch_n_events_opt;

  ////////////////////////////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - multi-channel Monte Carlo for final-state variables  //
  ////////////////////////////////////////////////////////////////////////////////////////////

  int switch_MC;
  int n_alpha_steps;
  int n_alpha_events;
  int n_alpha_epc;
  int MCweight_min;
  double MCweight_limit_min;
  double MCweight_limit_max;
  int switch_use_alpha_after_IS;
  int switch_step_mode_grid;


  string MCweight_in_contribution;
  // ??? needed ??? maybe location of grid to be read in...

  string MCweight_in_directory;
  // ??? needed ???

  /////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - IS for final-state variables  //
  /////////////////////////////////////////////////////////////////////

  int switch_IS_MC;
  int switch_IS_mode_phasespace;
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
 
  /////////////////////////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - multi-channel Monte Carlo for CMS-energy mapping  //
  /////////////////////////////////////////////////////////////////////////////////////////

  int switch_MC_tau;

  //////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - IS for CMS-energy mapping  //
  //////////////////////////////////////////////////////////////////

  int switch_IS_tau;
  int n_tau_steps;
  int n_tau_events;
  int n_tau_bins;

  ////////////////////////////////////////////////////////////
  //  weight-optimization parameters - IS for x1x2 mapping  //
  ////////////////////////////////////////////////////////////

  int switch_IS_x1x2;
  int n_x1x2_steps;
  int n_x1x2_events;
  int n_x1x2_bins;

  ///////////////////////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - IS for z1 and z2 mappings (collinear emission)  //
  ///////////////////////////////////////////////////////////////////////////////////////

  int switch_IS_z1z2;
  int n_z1z2_steps;
  int n_z1z2_events;
  int n_z1z2_bins;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - multi-channel Monte Carlo for x-values in dipole mappings  //
  //////////////////////////////////////////////////////////////////////////////////////////////////

  int switch_MC_x_dipole;

  ////////////////////////////////////////////////////////////////////////////////
  //  weight-optimization parameters - IS for xy/zuv values in dipole mappings  //
  ////////////////////////////////////////////////////////////////////////////////

  int n_IS_gridsize_xy;
  int n_IS_gridsize_zuv;



  /////////////////////////////////////////
  //  phase-space generation parameters  //
  /////////////////////////////////////////

  double nu;
  // removed ???

  double nuxs;
  double nuxt;
  double exp_pdf;
  double exp_pT;
  // removed ???

  double exp_y;
  // removed ???

  double exp_ij_k_y;
  double exp_ij_k_z;
  double exp_ij_a_x;
  double exp_ij_a_z;
  double exp_ai_k_x;
  double exp_ai_k_u;
  double exp_ai_b_x;
  double exp_ai_b_v;

  ////////////////////////////////////////
  //  technical integration parameters  //
  ////////////////////////////////////////

  double mass0;
  double map_technical_s;
  double map_technical_t;
  double map_technical_x;
  double cut_technical;

  //////////////////////////////
  //  integration parameters  //
  //////////////////////////////

  double sigma_normalization;
  int zwahl;
  double sigma_normalization_deviation;
  int n_events_max;
  int n_events_min;
  int n_step;

  ////////////////////////////////////////////////////////////////////////////////////////
  //  parameters directly forwarded to OpenLoops - after the default settings are done  //
  ////////////////////////////////////////////////////////////////////////////////////////

  vector<string> OL_parameter;
  vector<string> OL_value;

  ////////////////////////////////////////
  //  scale variation parameters - TSV  //
  ////////////////////////////////////////

  int switch_TSV;

  vector<string> name_set_TSV;
  int n_set_TSV;

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

  vector<double> central_scale_TSV;
  vector<double> central_scale_ren_TSV;
  vector<double> central_scale_fact_TSV;
  vector<double> relative_central_scale_TSV;
  vector<double> relative_central_scale_ren_TSV;
  vector<double> relative_central_scale_fact_TSV;

  vector<int> n_scale_TSV;
  vector<int> n_scale_ren_TSV;
  vector<int> n_scale_fact_TSV;
  vector<int> factor_scale_TSV;
  vector<int> factor_scale_ren_TSV;
  vector<int> factor_scale_fact_TSV;

  vector<int> dynamic_scale_TSV;
  vector<int> dynamic_scale_ren_TSV;
  vector<int> dynamic_scale_fact_TSV;
  vector<double> min_qTcut_TSV;
  vector<double> max_qTcut_TSV;
  vector<int> switch_distribution_TSV;
  vector<double> min_qTcut_distribution_TSV;
  vector<double> max_qTcut_distribution_TSV;
  vector<int> switch_moment_TSV;

  vector<int> max_n_integrand_TSV;
  // !!! no input parameter !!!

  vector<int> no_central_scale_ren_TSV;
  // !!! no input parameter !!!

  vector<int> no_central_scale_fact_TSV;
  // !!! no input parameter !!!

  
  ///////////////////////////////////////////////
  //  selection of run_mode (grid, time, run)  //
  ///////////////////////////////////////////////

  string run_mode;


  void parameter_readin_file(string filename, vector<string> & readin, bool essential, int retry);
  void parameter_readin(string & subprocess, vector<string> & readin);
  void get_userinput_from_readin(vector<string> & user_variable, vector<string> & user_variable_additional, vector<string> & user_value, vector<string> & readin);

  friend std::ostream & operator << (std::ostream &, const inputparameter_set &);

};
#endif
