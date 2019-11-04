#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "logger.h"
#include "randommanager.h"

using std::vector;

class randommanager;

class randomvariable {
public:
  randomvariable(int _n_events, int _switch_IS_MC, int _n_steps, int _n_bins, const std::string &_name);
  //  void optimize(const double psp_weight);

  void vegasgrid_average();
  void psp_IS_optimization(double & psp_weight);
  //  void step_IS_optimization_temp();
  void step_IS_optimization(int events);
  void result_IS_optimization(int i_step_mode);
  
  void get_random(double &x, double &g_IS);
  void get_g_IS(double r, double &g_IS);
  void multiply_g_IS(double r, double & g_IS);
  //  double calc_variation();
  void proceeding_out(std::ofstream &out_proceeding);
  void proceeding_in(int &proc, vector<std::string> &readin);
  void check_proceeding_in(int &int_end, int &proc, vector<std::string> &readin);
  //  void proceeding_size(int &proc);
  //  void proceeding_size(vector<int> & save_proceeding);
  void proceeding_save();
  void save_weights(std::ofstream &out_weights);
  void readin_weights(vector<string> & readin);
  //  int readin_weights(std::vector<std::string> &readin);

  //   int do_optimization_step(int events);

  
  int channel;

  int n_events;
  int switch_IS_MC;
  int n_steps;
  int n_bins;
  //  int n_poweroftwo;

  //  int opt_end;
  bool used;

  bool imp_sampling;
  
  vector<double> s;
  vector<double> save_s;

  vector<double> alpha;
  vector<double> beta;
  vector<int> n_rej_channel;
  vector<int> n_acc_channel;
  vector<double> sum_channel_weight;
  vector<double> sum_channel_weight2;

  int active_optimization;
  int end_optimization;

  std::string name;
  std::string filename;

  randommanager *manager;

private:
  void init(int _n_events, int _switch_IS_MC, int _n_steps, int _n_bins, const std::string & _name);
  void get_random(double &x, double &g_IS, vector<double> &s);

  //  int save_opt_end;
  int save_end_optimization;

  vector<double> save_alpha;
  vector<double> save_beta;
  vector<int> save_n_acc_channel;
  vector<int> save_n_rej_channel;
  vector<double> save_sum_channel_weight;
  vector<double> save_sum_channel_weight2;
};

#endif
