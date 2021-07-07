#ifndef SUMMARY_CONTRIBUTION_H
#define SUMMARY_CONTRIBUTION_H

//#include "summary.subprocess.h"

class summary_list; // check if reasonable !!!
class summary_generic;

class summary_contribution {
private:

public:

  ////////////////////
  //  constructors  //
  ////////////////////

  summary_contribution();
  summary_contribution(string _type_contribution, summary_list & _list);



  ///////////////////////
  //  access elements  //
  ///////////////////////

  // each 'summary_contribution' belongs to a specific 'summary_list'
  summary_list *ylist; // check if reasonable !!!
  // only one 'summary_generic'
  summary_generic *ygeneric;
  // only one 'observable_set'
  observable_set *osi;


  // same as in list...
  string processname;
  string resultdirectory; // not filled here !!!
  string type_perturbative_order;
  string type_subtraction_method;
  int in_contribution_order_alpha_s;
  int in_contribution_order_alpha_e;
  //

  string type_contribution;
  string type_correction;
  int interference;
  int photon_induced;

  // to later replace directory and directory_extra
  vector<string> extended_directory_name;
  vector<vector<string> > extended_directory;

  vector<string> directory;
  vector<int> directory_extra;
  string infix_contribution;
  string infix_order_contribution;
  string infix_path_contribution;

  int active_qTcut;
  int output_n_qTcut;
  int selection_n_qTcut;

  vector<string> subprocess;
  vector<summary_subprocess> xsubprocess;
  vector<vector<int > > subgroup_no_member;

  int max_number_of_jobs_for_one_channel;

  int average_factor;

  vector<vector<vector<vector<double> > > > result_CV;
  vector<vector<vector<vector<double> > > > deviation_CV;
  // 1st: i_m - n_moments + 1 (1 -> cross section only, >1 -> cross section + moments)
  // 2nd: i_g - n_subgroups (number of subprocess groups: 0 entry -> sum over all subgroups, 1...subgroup.size(): contributions from gg, gq, qq~ etc. initial states (as defined in process))
  // 3rd: i_q - osi_n_qTcut (number of qTcut values)
  // 4th: i_s - osi_n_scales_CV (number of scales evaluated for scale-variation plots in CV-files  )

  vector<vector<vector<vector<double> > > > distribution_result_CV;
  vector<vector<vector<vector<double> > > > distribution_deviation_CV;


  vector<vector<vector<vector<vector<vector<double> > > > > > result_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > deviation_TSV;
  // 1st: n_moments + 1 (1 -> cross section only, >1 -> cross section + moments)
  // 2nd: n_subgroups (number of subprocess groups: 0 entry -> sum over all subgroups, 1...subgroup.size(): contributions from gg, gq, qq~ etc. initial states (as defined in process))
  // 3rd: osi_n_qTcut (number of qTcut values)
  // 4th: n_scale_set_TSV (number of scale sets in TSV)
  // 5th: n_scales_ren_TSV (number of renormalization scales evaluated for scale-variation plots in TSV-files)
  // 6th: n_scales_fact_TSV (number of renormalization scales evaluated for scale-variation plots in TSV-files)

  ///  vector<double> n_bin_distribution_modasym;
  vector<int> n_bin_distribution_modasym;

  vector<vector<vector<vector<vector<vector<vector<double> > > > > > > distribution_result_TSV;
  vector<vector<vector<vector<vector<vector<vector<double> > > > > > > distribution_deviation_TSV;
  // 1st: i_g - subgroup.size() (number of subprocess groups: 0 entry -> sum over all subgroups, 1...subgroup.size(): contributions from gg, gq, qq~ etc. initial states (as defined in process files))
  // 2nd: i_d - n_distribution (number of (dd) distributions)
  // 3rd: i_b - n_bin (number of bins of the respective distribution)
  // 4th: x_q - selection_n_qTcut (number of qTcut values selected for distributions)
  // 5th: i_s - n_scale_set_TSV (number of scale sets in TSV)
  // 6th: i_r - n_scales_ren_TSV (number of renormalization scales evaluated for scale-variation plots in TSV-files)
  // 7th: i_f - n_scales_fact_TSV (number of renormalization scales evaluated for scale-variation plots in TSV-files)

  ///////////////
  //  methods  //
  ///////////////

  void readin_contribution_remove_run();
  void readin_runtime_contribution();

  //  summary/summary.contribution.result.cpp
  void readin_result_contribution_CV(string result_moment, int x_m);
  void readin_result_contribution_TSV(string result_moment, int x_m);

  //  summary/summary.contribution.output.result.cpp
  void output_result_overview_qTcut_TSV();
  void output_result_overview_qTcut_CV();
  void output_result_qTcut_TSV();
  void output_result_qTcut_CV();
  void output_result_plot_qTcut_TSV(); 
  void output_result_plot_qTcut_CV(); 

  void output_result_plot_CV();

  void output_subprocesses_result_plot_qTcut_TSV(); 
  void output_subprocesses_result_plot_qTcut_CV(); 

  void output_subprocesses_result_plot_CV();
  //  void output_contribution_result_plot_qTcut_CV(string result_moment, int x_m); // redundant (also in summary.list) !!!
  void output_result_check_CV();

  //  summary/summary.contribution.distribution.cpp
  //vector<vector<double> > & bin_edge
  void readin_distribution_contribution_CV();
  void readin_distribution_contribution_TSV();

  //  summary/summary.contribution.output.distribution.cpp
  void output_distribution_overview_qTcut_TSV();
  void output_distribution_overview_qTcut_CV();
  void output_distribution_norm_qTcut_TSV();
  void output_distribution_plot_qTcut_TSV();
  void output_distribution_plot_CV();

  void output_distribution_asymmetry_CV();

  void output_extrapolated_runtime_directories(int i_m, int & list_contribution_counter);
};
#endif
