#ifndef RANDOMMANAGER_H
#define RANDOMMANAGER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "logger.h"
#include "randomvariable.h"

using std::vector;

class randomvariable;

class randommanager{
public:
  randommanager();
  randommanager(vector<double> sran, std::string _weightdir, std::string _processname);
  ~randommanager();

  int active_optimization;
  int end_optimization;

  void register_variable(randomvariable* new_var);
  //  void register_variable(randomvariable* new_var, bool imp_sampling);
  
  void optimize(const double psp_weight, const double psp_weight2);
  void psp_IS_optimization(double & psp_weight, double & psp_weight2);
  
  void increase_cut_counter();
  void increase_counter_n_acc();

  void do_optimization_step(int events);
  void proceeding_out(std::ofstream &out_proceeding);
  void proceeding_in(int &proc, vector<std::string> &readin);
  void check_proceeding_in(int &int_end, int &proc, vector<std::string> &readin);
  //  void proceeding_size(int &proc);
  //  void proceeding_size(vector<int> & save_proceeding);
  void proceeding_save();
  void writeout_weights();

  void add_var_to_queue(randomvariable* variable);
  //  void nancut(bool count_as_cut);
  void nancut();

private:
  void readin_weights();

  vector<randomvariable*> random_variables;

  //  bool opt_end;
  vector<double> s;
  //  Logger logger;

  std::string weightdir;
  std::string processname;

  int save_end_optimization;

  vector<double> save_s;

  vector<randomvariable*> used_queue;
  int queue_counter;
  int queue_threshold;

  vector<std::string> readin;
};

#endif
