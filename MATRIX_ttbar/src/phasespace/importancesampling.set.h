#ifndef IMPORTANCESAMPLING_SET_H
#define IMPORTANCESAMPLING_SET_H

class importancesampling_set {
private:

public:
  // constant over full importancesampling
  string name;

  string filename;
  string filename_readin;

  int n_gridsize;
  //int n_x1x2_bins;
  int n_optimization_step;
  //int n_x1x2_steps;
  int n_event_per_step;
  //int n_x1x2_events;
  
  int switch_minimum_weight;
  double limit_minimum_weight;
  double reserved_minimum_weight;

  // used only inside optimization routines
  // could be private
  double a_reserved_min;

  // chanes for each phase-space point
  int channel;
  vector<double> g_channel;
  vector<double> g_IS_channel;

  // changes at optimization steps
  vector<double> alpha;
  //vector<double> x1x2_alpha;
  vector<double> beta;
  //vector<double> x1x2_beta;
  
  // needed for optimization phase
  int active_optimization;
  int end_optimization;
  //int x1x2_opt_end;

  int n_gen;
  int n_acc;
  int n_rej;
  int n_tec;
  int n_nan;

  vector<int> n_acc_channel;
  //vector<int> x1x2_counts_channel;
  vector<int> n_rej_channel;
  //vector<int> x1x2_cuts_channel;
  vector<double> sum_channel_weight;
  //vector<double> sum_x1x2_channel_weight;
  vector<double> sum_channel_weight2;
  //vector<double> sum_x1x2_channel_weight2;
  
  /*
  vector<double> w_channel_step_sum;
  vector<double> w_channel_full_sum;
  //  vector<double> MC_sum_w_channel;
  vector<double> w_channel_av;
  */
  
  // could be private
  /*
  int counter_minimum_weight;
  vector<int> no_channel_minimum_weight;
  */
  /*
  int i_alpha_it;
  vector<vector<double> > alpha_it;
  vector<double> diff_w; // -> MC_diff_w_step
  int x_minimum_diff_w;
  int x_alpha_it_min;
  */
  
////////////////////
//  constructors  //
////////////////////
  importancesampling_set();
  importancesampling_set(int _n_channel);
  importancesampling_set(string _name, int _n_gridsize, int _n_optimization_step, int _n_event_per_step, int _switch_minimum_weight, double _limit_minimum_weight, double _reserved_minimum_weight, int _active_optimization, int _end_optimization);
  importancesampling_set(string _name, int _n_gridsize, int _n_optimization_step, int _n_event_per_step, int _switch_minimum_weight, double _limit_minimum_weight, double _reserved_minimum_weight, int _switch_IS, string _filename, string _filename_readin);

///////////////
//  methods  //
///////////////
  void vegasgrid_average();
  void vegasgrid_calculation();
  void psp_IS_optimization(double & this_psp_weight, double & this_psp_weight2);
  void step_IS_optimization(int i_step_mode);
  void result_IS_optimization(int i_step_mode);
  void output_IS_optimization(int i_step_mode);
  void readin_IS_optimization();

///////////////////////
//  access elements  //
///////////////////////

};
#endif
