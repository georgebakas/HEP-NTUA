#ifndef SUMMARY_GENERIC_H
#define SUMMARY_GENERIC_H

class summary_generic {
private:

public:

////////////////////
//  constructors  //
////////////////////
  summary_generic();
  summary_generic(munich & xmunich);
  summary_generic(string _order, string & _infilename, string & _infilename_scaleband, inputparameter_set & isi, call_generic & generic);

///////////////////////
//  access elements  //
///////////////////////

  inputparameter_set * isi;
  call_generic * generic;
  contribution_set csi;
  observable_set oset;

  string order;
  string infilename;
  string infilename_scaleband;


  string final_resultdirectory;
  vector<string> list_contribution_file;
  map<string, int> mapping_contribution_file;
  vector<summary_order> yorder;

  int run_min_n_step;
  int run_factor_max ;
  long long run_min_n_event;
  double run_max_time_per_job;
  int run_min_number_of_jobs_for_one_channel;
  int average_factor;
  int switch_generate_rundirectories;
  int no_qTcut_runtime_estimate;
  string directory_runtime_estimate;
  double deviation_tolerance_factor;

  int switch_extrapolation_result;
  int switch_extrapolation_distribution;

  double min_qTcut_extrapolation;
  double max_qTcut_extrapolation;
  double min_max_value_extrapolation_range;
  int min_n_qTcut_extrapolation_range;
  double error_extrapolation_range_chi2;

  vector<string> phasespace_optimization;

  vector<string> inpath_scalecentral;
  vector<vector<string> > inpath_scalevariation;
  vector<string> outpath_scaleband;
  vector<string> outpath_scalecentral;
  vector<string> outpath_scalemax;
  vector<string> outpath_scalemin;
  
  vector<string> outpath_scaleplusmax;
  vector<string> outpath_scaleplusmin;
  vector<string> outpath_scaleminusmax;
  vector<string> outpath_scaleminusmin;

  // only used for asymmetries: could most likely be removed later...
  vector<vector<double> > bin_edge;

  //  string processname; // needed ???
  vector<string> subgroup;
  vector<summary_list> xlist;
  vector<string> resultdirectory; // maybe a better name ???

  vector<string> name_scale_variation_TSV;
  vector<vector<vector<string> > > scalename_TSV;
  string name_variation_CV;
  vector<string> scalename_CV;

  vector<string> output_selection_distribution;
  vector<int> switch_output_distribution;

  vector<string> output_selection_scaleset;
  vector<int> switch_output_scaleset;

  int switch_output_subprocess;
  int switch_output_contribution;
  int switch_output_list;
  int switch_output_order;

  int switch_output_plot;
  int switch_output_result;
  int switch_output_overview;

  // new:
  int switch_output_table_order;
  int switch_output_table_Kfactor;
  int switch_output_table_crosssection_Kfactor;
  int switch_output_table_IS_splitting;

  vector<string> output_selection_distribution_table;
  vector<int> switch_output_distribution_table;


  vector<vector<vector<vector<double> > > > scaleband_variable;
  vector<vector<vector<vector<vector<double> > > > > scaleband_central_result;
  vector<vector<vector<vector<vector<double> > > > > scaleband_central_deviation;
  vector<vector<vector<vector<vector<double> > > > > scaleband_minimum_result;
  vector<vector<vector<vector<vector<double> > > > > scaleband_minimum_deviation;
  vector<vector<vector<vector<vector<double> > > > > scaleband_maximum_result;
  vector<vector<vector<vector<vector<double> > > > > scaleband_maximum_deviation;

  

  vector<string> infix_name_moment;
  int x_m;

  string filename_complete;
  string filename_ren;
  string filename_fact;
  string filename_equal;
  string filename_antipodal;
  
  string filename_CV;

  /*
  ofstream outfile_complete;
  ofstream outfile_ren;
  ofstream outfile_fact;
  ofstream outfile_equal;
  ofstream outfile_antipodal;

  ofstream outfile_CV;
  */

  void get_summary();

  void readin_combination_infile();
  void readin_infile_scaleband();
  void initialization_distribution();
  void initialization_scaleset_TSV();

  void initialization_summary_list();

  void collect_order_result_TSV();
  void collect_order_result_CV();

  void determine_runtime_result_new();

  void collect_order_distribution_TSV();
  void collect_order_distribution_CV();

  void determine_scaleband();
  void determine_runtime();

  void output_filename_TSV(string & filename_begin, string & filename_end);
  //  void output_TSV(stringstream & temp_ss);

  void output_distribution_table_order();
  void output_distribution_table_Kfactor();
  void output_distribution_table_crosssection_Kfactor();
  void output_distribution_table_IS_splitting();

};
#endif
