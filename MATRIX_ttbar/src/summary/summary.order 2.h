#ifndef SUMMARY_ORDER_H
#define SUMMARY_ORDER_H

class summary_generic;

class summary_order {
private:

public:

////////////////////
//  constructors  //
////////////////////
  summary_order();
  summary_order(string _resultdirectory, summary_generic & _generic);



///////////////////////
//  access elements  //
///////////////////////

  // only one 'summary_generic'
  summary_generic *ygeneric;
  // only one 'observable_set'
  observable_set *osi;

  string processname;
  //  string resultdirectory;
  string type_perturbative_order;
  string type_subtraction_method;
  //  int in_contribution_order_alpha_s;
  //  int in_contribution_order_alpha_e;

  //  vector<summary_contribution> xcontribution;

  string resultdirectory;
  vector<string> contribution_file;
  vector<vector<int> > combination;
  vector<int> combination_type;
  double accuracy_relative;
  string accuracy_normalization;
  int accuracy_no_normalization;

  double error2_time;

  //  vector<summary_list> xlist;
  vector<summary_list*> xlist;

  int active_qTcut;
  int output_n_qTcut;
  int selection_n_qTcut;





  vector<vector<vector<double> > > result_CV;
  vector<vector<vector<double> > > deviation_CV;
  // 1st: i_g - subgroup.size() (number of subprocess groups: 0 entry -> sum over all subgroups)
  // 2nd: i_m - n_moment (number of moments)
  // 3rd: i_s - n_scales_CV (number of scales in CV)

  vector<vector<vector<vector<double> > > > result_qTcut_CV;
  vector<vector<vector<vector<double> > > > deviation_qTcut_CV;
  // 1st: i_g - subgroup.size() (number of subprocess groups: 0 entry -> sum over all subgroups)
  // 2nd: i_m - n_moment (number of moments)
  // 3rd: i_q - output_n_qTcut (number of qTcut values selected for results)
  // 4th: i_s - n_scales_CV (number of scales in CV)


  vector<vector<vector<vector<vector<double> > > > > result_TSV;
  vector<vector<vector<vector<vector<double> > > > > deviation_TSV;
  // 1st: i_g - subgroup.size() (number of subprocess groups: 0 entry -> sum over all subgroups)
  // 2nd: i_m - n_moment (number of moments)
  // 3rd: i_s - n_scale_set_TSV (number of scale sets in TSV)
  // 4th: i_r - n_scales_ren_TSV (number of renormalization scales in TSV files)
  // 5th: i_f - n_scales_fact_TSV (number of renormalization scales in TSV files)

  vector<vector<vector<vector<vector<vector<double> > > > > > result_qTcut_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > deviation_qTcut_TSV;
  // 1st: i_g - subgroup.size() (number of subprocess groups: 0 entry -> sum over all subgroups)
  // 2nd: i_m - n_moment (number of moments)
  // 3rd: i_q - output_n_qTcut (number of qTcut values selected for results)
  // 4th: i_s - n_scale_set_TSV (number of scale sets in TSV)
  // 5th: i_r - n_scales_ren_TSV (number of renormalization scales in TSV files)
  // 6th: i_f - n_scales_fact_TSV (number of renormalization scales in TSV files)



  vector<vector<vector<vector<double> > > > distribution_result_CV;
  vector<vector<vector<vector<double> > > > distribution_deviation_CV;

  vector<vector<vector<vector<vector<vector<double> > > > > > distribution_result_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > distribution_deviation_TSV;
  // 1st: i_g - subgroup.size() (number of subprocess groups: 0 entry -> sum over all subgroups)
  // 2nd: i_d - n_distribution (number of (dd) distributions)
  // 3rd: i_b - n_bin (number of bins of the respective distribution)
  // 4th: i_s - n_scale_set_TSV (number of scale sets in TSV)
  // 5th: i_r - n_scales_ren_TSV (number of renormalization scales in TSV files)
  // 6th: i_f - n_scales_fact_TSV (number of renormalization scales in TSV files)

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
  void collect_result_TSV(vector<string> & subgroup);
  void collect_result_CV(vector<string> & subgroup);

  void output_result_TSV();
  void output_result_overview_TSV();
  void output_result_plot_TSV();
  void output_result_plot_qTcut_TSV();

  void output_result_overview_CV();
  void output_result_CV();
  void output_result_plot_CV();
  void output_result_plot_qTcut_CV();

  //  summary/summary.order.distribution.cpp
  //  void collect_distribution_TSV(string & final_resultdirectory, vector<string> & subgroup, vector<xdistribution> & extended_distribution, vector<double> & fakeasymfactor, observable_set & oset, vector<vector<vector<string> > > & scalename_TSV);
  void collect_distribution_TSV();
  void collect_distribution_CV();

  //  summary/summary.order.output.distribution.cpp
  void output_distribution_overview_TSV();
  void output_distribution_overview_CV();

  void output_distribution_TSV();
  void output_distribution_qTcut_TSV();

  void output_sddistribution_TSV(string & identifier, vector<vector<vector<vector<double> > > > & this_distribution_result_TSV, vector<vector<vector<vector<double> > > > & this_distribution_deviation_TSV, int i_d, string & subdirectory);
  void output_dddistribution_reconstruct_sdd_TSV(string & identifier, vector<vector<vector<vector<double> > > > & this_distribution_result_TSV, vector<vector<vector<vector<double> > > > & this_distribution_deviation_TSV, int i_d, string & subdirectory);
  void output_dddistribution_split_in_first_sdd_TSV(string & identifier, vector<vector<vector<vector<double> > > > & this_distribution_result_TSV, vector<vector<vector<vector<double> > > > & this_distribution_deviation_TSV, int i_d, string & subdirectory);
  void output_dddistribution_ge_lt_in_first_sdd_TSV(string & identifier, vector<vector<vector<vector<double> > > > & this_distribution_result_TSV, vector<vector<vector<vector<double> > > > & this_distribution_deviation_TSV, int i_d, string & subdirectory);

  void output_distribution_CV();
  void output_sddistribution_CV(string & identifier, vector<vector<double> > & this_distribution_result_CV, vector<vector<double> > & this_distribution_deviation_CV, int i_d, string & subdirectory);

  void output_dddistribution_split_in_first_sdd_CV(string & identifier, vector<vector<double> > & this_distribution_result_CV, vector<vector<double> > & this_distribution_deviation_CV, int i_d, string & subdirectory);
};
#endif
