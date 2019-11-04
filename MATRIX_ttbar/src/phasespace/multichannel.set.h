#ifndef MULTICHANNEL_SET_H
#define MULTICHANNEL_SET_H

class multichannel_set {
private:

public:
  // constant over full multichannel
  string name;

  string filename;
  string filename_readin;
    
  int n_channel;
  int n_alpha_events;
  int n_alpha_steps;
  int switch_minimum_weight;
  double limit_minimum_weight;
  double reserved_minimum_weight;

  // used only inside optimization routines
  // could be private
  double a_reserved_min;

  // changes for each phase-space point
  int channel;
  vector<double> g_channel;
  vector<double> g_IS_channel;

  // changes at optimization steps
  vector<double> alpha;
  vector<double> beta;

  // needed for optimization phase
  //  int end_optimization;

  int active_optimization;
  int end_optimization;

  int n_gen;
  int n_acc;
  int n_rej;
  int n_tec;
  int n_nan;

  vector<int> n_acc_channel;
  vector<int> n_rej_channel;
  vector<double> w_channel_step_sum;
  vector<double> w_channel_full_sum;
  //  vector<double> MC_sum_w_channel;
  vector<double> w_channel_av;

  // could be private
  int counter_minimum_weight;
  vector<int> no_channel_minimum_weight;

  int i_alpha_it;
  vector<vector<double> > alpha_it;
  vector<double> diff_w; // -> MC_diff_w_step
  int x_minimum_diff_w;
  int x_alpha_it_min;

////////////////////
//  constructors  //
////////////////////
  multichannel_set();
  multichannel_set(int _n_channel);
  multichannel_set(string _name, int _n_channel, int _n_alpha_events, int _n_alpha_steps, int _switch_minimum_weight, double _limit_minimum_weight, double _reserved_minimum_weight, int _active_optimization, int _end_optimization);
  multichannel_set(string _name, int _n_channel, int _n_alpha_events, int _n_alpha_steps, int _switch_minimum_weight, double _limit_minimum_weight, double _reserved_minimum_weight, int _switch_MC, string _filename, string _filename_readin);
  
  void psp_MCweight_optimization(double & integrand, double & g_tot);
  void step_MCweight_optimization(int i_acc);
  void result_MCweight_optimization(int i_acc);
  void output_MCweight_optimization(int i_step_mode, string & filename_MCweight);
  void readin_MCweight_optimization();

};
#endif
