#ifndef SUMMARY_LIST_H
#define SUMMARY_LIST_H

class summary_generic;

class summary_list {
private:

public:

////////////////////
//  constructors  //
////////////////////
  summary_list();
  summary_list(string & infilename, summary_generic & _generic);

///////////////////////
//  access elements  //
///////////////////////

  // only one 'summary_generic'
  summary_generic *ygeneric;
  // only one 'observable_set'
  observable_set *osi;

  string processname;
  string resultdirectory;
  string type_perturbative_order;
  string type_subtraction_method;
  int in_contribution_order_alpha_s;
  int in_contribution_order_alpha_e;
  int photon_induced;

  vector<summary_contribution> xcontribution;

  int active_qTcut;
  int output_n_qTcut;
  int selection_n_qTcut;

  vector<vector<vector<double> > > result_CV;
  vector<vector<vector<double> > > deviation_CV;
  vector<vector<vector<double> > > deviation_extrapolation_CV;
  vector<vector<vector<double> > > deviation_statistics_CV;

  vector<vector<vector<vector<vector<double> > > > > result_TSV;
  vector<vector<vector<vector<vector<double> > > > > deviation_TSV;
  vector<vector<vector<vector<vector<double> > > > > deviation_extrapolation_TSV;
  vector<vector<vector<vector<vector<double> > > > > deviation_statistics_TSV;
  // 1st: n_moments + 1 (1 -> cross section only, >1 -> cross section + moments)
  // 2nd: n_subgroups (number of subprocess groups: 0 entry -> sum over all subgroups, 1...subgroup.size(): contributions from gg, gq, qq~ etc. initial states (as defined in process))
  // 3rd: n_scale_set_TSV (number of scale sets in TSV)
  // 4th: n_scales_ren_TSV (number of renormalization scales evaluated for scale-variation plots in TSV-files)
  // 5th: n_scales_fact_TSV (number of renormalization scales evaluated for scale-variation plots in TSV-files)



  // missing: distribution_result_CV !!!
  // probably identical to xcontribution[0].distribution_result_CV



  vector<vector<vector<vector<vector<vector<double> > > > > > distribution_result_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > distribution_deviation_TSV;
  // 1st: i_g - subgroup.size() (number of subprocess groups: 0 entry -> sum over all subgroups, 1...subgroup.size(): contributions from gg, gq, qq~ etc. initial states (as defined in process files))
  // 2nd: i_d - n_distribution (number of (double-differential) distributions)
  // 3rd: i_b - n_bin (number of bins of the respective distribution)
  // 4th: i_s - n_scale_set_TSV (number of scale sets in TSV)
  // 5th: i_r - n_scales_ren_TSV (number of renormalization scales evaluated for scale-variation plots in TSV-files)
  // 6th: i_f - n_scales_fact_TSV (number of renormalization scales evaluated for scale-variation plots in TSV-files)

  vector<vector<vector<vector<vector<vector<vector<double> > > > > > > distribution_result_qTcut_TSV;
  vector<vector<vector<vector<vector<vector<vector<double> > > > > > > distribution_deviation_qTcut_TSV;
  // 1st: x_q - selection_n_qTcut (number of qTcut values selected for distributions)
  // 2nd: i_g - subgroup.size() (number of subprocess groups: 0 entry -> sum over all subgroups)
  // 3rd: i_d - n_distribution (number of (dd) distributions)
  // 4th: i_b - n_bin (number of bins of the respective distribution)
  // 5th: i_s - n_scale_set_TSV (number of scale sets in TSV)
  // 6th: i_r - n_scales_ren_TSV (number of renormalization scales evaluated in TSV files)
  // 7th: i_f - n_scales_fact_TSV (number of renormalization scales evaluated in TSV files)

  ///////////////
  //  methods  //
  ///////////////

  //  summary/summary.list.cpp

  void output_info();
  void calculate_runtime();

  //  void output_extrapolated_runtime_directories(ofstream & out_runscript_complete, ofstream & out_runlist_complete);

  void output_readin_extrapolated_runtime(ofstream & outfile_readin_extrapolated);
  void output_extrapolated_runtime(ofstream & outfile_extrapolated);
  void output_used_runtime(ofstream & outfile_used);

  void extrapolation_quadratic(vector<double> & xvalue, double & min_extrapolation, double & max_extrapolation, vector<double> & data_result, vector<double> & data_deviation2, vector<double> & X, vector<double> & dX, double & chi2_dof);
  void extrapolation_linear(vector<double> & xvalue, double & min_extrapolation, double & max_extrapolation, vector<double> & data_result, vector<double> & data_deviation2, vector<double> & X, vector<double> & dX, double & chi2_dof);

  void extrapolation_TSV(int order_fit);
  void extrapolation_CV(int order_fit);

  void extrapolation_TSV(vector<double> & xvalue, vector<vector<vector<vector<double> > > > & data_result_TSV, vector<vector<vector<vector<double> > > > & data_deviation_TSV, vector<vector<vector<double> > > & extrapolation_result_TSV, vector<vector<vector<double> > > & extrapolation_deviation_TSV, int order_fit);
  void extrapolation_CV(vector<double> & xvalue, vector<vector<double> > & data_result_CV, vector<vector<double> > & data_deviation_CV, vector<double> & extrapolation_result_CV, vector<double> & extrapolation_deviation_CV, int order_fit);

  void extrapolation_fit_range_determination(vector<double> & xvalue, vector<double> & data_result, vector<double> & data_deviation, vector<double> & all_extrapolation_result, vector<double> & all_extrapolation_deviation, vector<double> & all_chi2_dof, int order_fit, int & no_xmax, double & xmin, double & xmax);

  //  summary/summary.list.result.cpp
  void collect_contribution_result_CV();
  void collect_contribution_result_TSV();

  //  summary/summary.list.output.result.cpp
  void output_contribution_result_overview_qTcut_CV();
  void output_contribution_result_plot_CV();
  void output_contribution_result_plot_qTcut_CV();
  void output_contribution_result_CV();
  void output_contribution_result_qTcut_CV();

  void output_contribution_result_overview_qTcut_TSV();
  void output_contribution_result_plot_TSV();
  void output_contribution_result_plot_qTcut_TSV();
  void output_contribution_result_TSV();
  void output_contribution_result_qTcut_TSV();

  //  summary/summary.list.distribution.cpp
  void collect_contribution_distribution_CV();
  void collect_contribution_distribution_TSV();

  //  summary/summary.list.output.distribution.cpp
  void output_distribution_overview_TSV();
  void output_distribution_overview_CV();

  void output_contribution_distribution_norm_TSV();
  void output_contribution_distribution_plot_TSV();
  void output_contribution_distribution_norm_qTcut_TSV();
  void output_contribution_distribution_plot_qTcut_TSV();
  void output_contribution_distribution_check_TSV();

};
#endif
