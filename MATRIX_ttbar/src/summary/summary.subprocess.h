#ifndef SUMMARY_SUBPROCESS_H
#define SUMMARY_SUBPROCESS_H

class summary_generic;
class summary_contribution;

class summary_subprocess {
private:

public:

  ////////////////////
  //  constructors  //
  ////////////////////

  summary_subprocess();
  summary_subprocess(string _name, summary_contribution & _contribution);
  summary_subprocess(string _name, int _n_seed, int _default_size_CV, summary_contribution & _contribution);

  ///////////////////////
  //  access elements  //
  ///////////////////////

  // each 'summary_subprocess' belongs to a specific 'summary_contribution'
  summary_contribution *ycontribution;
  // only one 'summary_generic'
  summary_generic *ygeneric;
  // only one 'observable_set'
  observable_set *osi;

  string name;
  int n_seed;
  int default_size_CV;

  int n_ps;
  
  double error2_time;
  double time_per_event;
  double used_runtime;
  long long used_n_event;

  double extrapolated_runtime;
  double extrapolated_n_event;
  double extrapolated_sigma_deviation;
  double extrapolated_sigma_normalization;
  double extrapolated_sigma_normalization_deviation;

  int max_number_of_jobs_for_one_channel;

  vector<double> result_run;
  vector<double> deviation_run;

  vector<vector<int> > remove_run_qTcut;
  vector<vector<int> > remove_run_qTcut_TSV;
  vector<vector<int> > remove_run_qTcut_CV;

  vector<double> error2_time_run;
  vector<double> time_per_event_run;
  vector<double> used_runtime_run;
  vector<long long> used_n_event_run;


  long long no_results;

  // only for reference result !!!
  vector<long long> N;
  long long N_all;


  double result;
  double deviation;
  double chi2;

  double Nd_event_CV;
  long long N_event_CV;
  vector<long long> N_valid_event_CV;
  // 1st: qTcut
  vector<vector<double> > result_CV;
  vector<vector<double> > deviation_CV;
  vector<vector<double> > chi2_CV;

 //  check meaning of n_parallel_runs_CV and related ones !!!
  vector<int> n_parallel_runs_reference_CV;
  vector<int> n_valid_parallel_runs_reference_CV;
  vector<vector<int> > n_parallel_runs_CV;
  vector<vector<int> > n_valid_parallel_runs_CV;

  vector<long long> N_event_run_CV;
  vector<vector<vector<double> > > result_run_CV;
  vector<vector<vector<double> > > deviation_run_CV;



  //  long long no_results_TSV;

  double Nd_event_TSV;
  long long N_event_TSV;
  vector<long long> N_valid_event_TSV;
  // 1st: qTcut
  vector<vector<vector<vector<vector<double> > > > > result_TSV;
  vector<vector<vector<vector<vector<double> > > > > deviation_TSV;
  vector<vector<vector<vector<vector<double> > > > > chi2_TSV;
  //  check meaning of n_parallel_runs_TSV !!!
  vector<int> n_parallel_runs_reference_TSV;
  vector<int> n_valid_parallel_runs_reference_TSV;
  vector<vector<vector<vector<vector<int> > > > > n_parallel_runs_TSV;
  vector<vector<vector<vector<vector<int> > > > > n_valid_parallel_runs_TSV;

  vector<long long> N_event_run_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > result_run_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > deviation_run_TSV;



  vector<vector<vector<double> > > distribution_result_CV;
  vector<vector<vector<double> > > distribution_deviation_CV;
  vector<vector<vector<double> > > distribution_chi2_CV;

  // in analogy to TSV:

  // total number of events for each [i_m] combination (summed over all [i_z], independent of [i_d][i_b][i_s]
  vector<long long> distribution_group_N_CV; // needed ???
  // binwise number of events for each [i_m][i_d][i_b][x_q] combination (summed over all [i_z], independent of [i_s]):
  vector<vector<vector<long long> > > distribution_group_N_binwise_CV;

  // combination of [i_m] combinations
  vector<vector<vector<vector<double> > > > distribution_group_result_CV;
  vector<vector<vector<vector<double> > > > distribution_group_deviation_CV;
  vector<vector<vector<vector<double> > > > distribution_group_chi2_CV;

  // total number of events for each [i_m][i_z] combination (independent of [i_d][i_b][i_s]):
  vector<vector<long long> > distribution_group_run_N_CV;
  // binwise number of events for each [i_m][i_z][i_d][i_b] combination (independent of [i_s]):
  vector<vector<vector<vector<long long> > > > distribution_group_run_N_binwise_CV;

  vector<vector<vector<vector<vector<double> > > > > distribution_group_run_result_CV;
  vector<vector<vector<vector<vector<double> > > > > distribution_group_run_deviation_CV;

  

  // total number of events from runs that have not been removed form a certain [i_m][i_d][i_b] result (removal determined separately for each [i_d][i_b] contribution):
  vector<vector<vector<long long> > > distribution_group_N_total_CV;
  // new: counts events per bin after removal of runs:
  vector<vector<vector<long long> > > distribution_group_N_total_binwise_CV;

  // counter of removed/nonzero runs (sum over all i_z) for certain [i_m][i_d][i_b]:
  vector<vector<vector<int> > > distribution_group_counter_removal_run_CV;
  vector<vector<vector<int> > > distribution_group_counter_nonzero_run_CV;

  //  switch for removal (0/1) of runs for certain [i_m][i_z][i_d][i_b]:
  vector<vector<vector<vector<int> > > > distribution_group_run_removal_run_CV;

  vector<vector<long long> > distribution_N_total_CV;
  // new: counts events per bin after removal of runs:
  vector<vector<long long> > distribution_N_total_binwise_CV;

  

  // TSV distribution

  // total number of contributing runs
  vector<long long> distribution_group_counter_run_TSV;

  // combination of [i_m] combinations
  vector<vector<vector<vector<vector<vector<double> > > > > > distribution_result_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > distribution_deviation_TSV;
  vector<vector<vector<vector<vector<vector<double> > > > > > distribution_chi2_TSV;
  
  // total number of events for each [i_m] combination (summed over all [i_z], independent of [i_d][i_b][x_q][i_s][i_r][i_f]
  vector<long long> distribution_group_N_TSV; // needed ???
  // binwise number of events for each [i_m][i_d][i_b][x_q] combination (summed over all [i_z], independent of [x_q][i_s][i_r][i_f]):
  vector<vector<vector<vector<long long> > > > distribution_group_N_binwise_TSV;

  vector<vector<vector<vector<vector<vector<vector<double> > > > > > > distribution_group_result_TSV;
  vector<vector<vector<vector<vector<vector<vector<double> > > > > > > distribution_group_deviation_TSV;
  vector<vector<vector<vector<vector<vector<vector<double> > > > > > > distribution_group_chi2_TSV;

  // total number of events for each [i_m][i_z] combination (independent of [i_d][i_b][x_q][i_s][i_r][i_f]):
  vector<vector<long long> > distribution_group_run_N_TSV;
  // binwise number of events for each [i_m][i_z][i_d][i_b][x_q] combination (independent of [i_s][i_r][i_f]):
  vector<vector<vector<vector<vector<long long> > > > > distribution_group_run_N_binwise_TSV;

  vector<vector<vector<vector<vector<vector<vector<vector<double> > > > > > > > distribution_group_run_result_TSV;
  vector<vector<vector<vector<vector<vector<vector<vector<double> > > > > > > > distribution_group_run_deviation_TSV;

  // total number of events from runs that have not been removed form a certain [i_m][i_d][i_b][x_q] result (removal determined separately for each [i_d][i_b][x_q] contribution):
  vector<vector<vector<vector<long long> > > > distribution_group_N_total_TSV;
  // new: counts events per bin after removal of runs:
  vector<vector<vector<vector<long long> > > > distribution_group_N_total_binwise_TSV;

  vector<vector<vector<long long> > > distribution_N_total_TSV;
  // new: counts events per bin after removal of runs:
  vector<vector<vector<long long> > > distribution_N_total_binwise_TSV;

  
  // counter of removed/nonzero runs (sum over all i_z) for certain [i_m][i_d][i_b][x_q]:
  vector<vector<vector<vector<int> > > > distribution_group_counter_removal_run_qTcut_TSV;
  vector<vector<vector<vector<int> > > > distribution_group_counter_nonzero_run_qTcut_TSV;

  //  switch for removal (0/1) of runs for certain [i_m][i_z][i_d][i_b][x_q]:
  vector<vector<vector<vector<vector<int> > > > > distribution_group_run_removal_run_qTcut_TSV;


  
  vector<vector<int> > type_parton;




  ///////////////
  //  methods  //
  ///////////////

  void readin_runtime(string & result_moment);

  void readin_result_CV(string & result_moment);
  void combination_result_CV(string & result_moment);
  void removal_exceptional_runs_result_CV();

  void combine_result_original();
  void combine_result_alternative();

  void combine_result_original_CV();
  void combine_result_alternative_CV();
  void combine_result_hybrid_CV();

  void readin_result_TSV(string & result_moment);
  void combination_result_TSV(string & result_moment);
  void removal_exceptional_runs_result_TSV();

  void combine_result_original_TSV();
  void combine_result_alternative_TSV();
  void combine_result_hybrid_TSV();

  void initialization_distribution_CV();
  void readin_distribution_CV();
  void combination_distribution_CV();
  void combine_distribution_conservative_CV(int i_d, int i_b, int i_s);

  //  void old_combination_distribution_CV();
  //  void old_combine_distribution_hybrid_CV(int x_z, int i_d, int i_b, int i_s);
  //  void old_combine_distribution_conservative_CV(int x_z, int i_d, int i_b, int i_s);

  //  void combine_distribution_original_CV();

  void initialization_distribution_TSV();
  void readin_distribution_TSV();
  void combination_distribution_TSV();
  void combine_distribution_hybrid_TSV(int x_q, int i_d, int i_b, int i_s, int i_r, int i_f);
  void combine_distribution_conservative_TSV(int x_q, int i_d, int i_b, int i_s, int i_r, int i_f);
  //  void combine_distribution_aggressive_TSV(int x_q, int i_d, int i_b, int i_s, int i_r, int i_f);
 
};
#endif
