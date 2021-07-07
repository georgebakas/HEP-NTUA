#include "../include/classes.cxx"

////////////////////
//  constructors  //
////////////////////

importancesampling_set::importancesampling_set(){}

importancesampling_set::importancesampling_set(int _n_channel){}

importancesampling_set::importancesampling_set(string _name, int _n_gridsize, int _n_optimization_step, int _n_event_per_step, int _switch_minimum_weight, double _limit_minimum_weight, double _reserved_minimum_weight, int _active_optimization, int _end_optimization){
  Logger logger("importancesampling_set::importancesampling_set " + _name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  name = _name;
  n_gridsize = _n_gridsize;
  n_optimization_step = _n_optimization_step;
  n_event_per_step = _n_event_per_step;

  switch_minimum_weight = _switch_minimum_weight;
  limit_minimum_weight = _limit_minimum_weight;
  reserved_minimum_weight = _reserved_minimum_weight;

  active_optimization = _active_optimization;
  end_optimization = _end_optimization;


  alpha.resize(n_gridsize, 1. / double(n_gridsize));
  beta.resize(n_gridsize);
  n_acc_channel.resize(n_gridsize);
  n_rej_channel.resize(n_gridsize);
  sum_channel_weight.resize(n_gridsize);
  sum_channel_weight2.resize(n_gridsize);
  
  if (n_gridsize > 0){
    beta[0] = alpha[0];
    for (int i = 1; i < beta.size(); i++){beta[i] = beta[i - 1] + alpha[i];}
  }

  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

importancesampling_set::importancesampling_set(string _name, int _n_gridsize, int _n_optimization_step, int _n_event_per_step, int _switch_minimum_weight, double _limit_minimum_weight, double _reserved_minimum_weight, int _switch_IS, string _filename, string _filename_readin){
  Logger logger("importancesampling_set::importancesampling_set " + _name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  name = _name;
  n_gridsize = _n_gridsize;
  n_optimization_step = _n_optimization_step;
  n_event_per_step = _n_event_per_step;

  switch_minimum_weight = _switch_minimum_weight;
  limit_minimum_weight = _limit_minimum_weight;
  reserved_minimum_weight = _reserved_minimum_weight;

  if (_switch_IS == 0 || _switch_IS == 2){
    active_optimization = 0;
    end_optimization = 1;
  }
  else if (_switch_IS == 1 || _switch_IS == 3){
    active_optimization = 1;
    end_optimization = 0;
  }
  else if (_switch_IS == -1){
    active_optimization = -1;
    end_optimization = -1;
  }
  else {logger << LOG_FATAL << "No valid value of switch_IS chosen." << endl; exit(1);}
  
  filename = _filename;
  filename_readin = _filename_readin;

  
  alpha.resize(n_gridsize, 1. / double(n_gridsize));
  beta.resize(n_gridsize);
  n_acc_channel.resize(n_gridsize);
  n_rej_channel.resize(n_gridsize);
  sum_channel_weight.resize(n_gridsize);
  sum_channel_weight2.resize(n_gridsize);
  
  logger << LOG_INFO << "_switch_IS = " << _switch_IS << endl;
  /*
  if (_switch_IS == 2 || _switch_IS == 3){
    readin_IS_optimization();
  }
  */
  
  if (n_gridsize > 0){
    beta[0] = alpha[0];
    for (int i = 1; i < beta.size(); i++){beta[i] = beta[i - 1] + alpha[i];}
  }

  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



///////////////////////
//  access elements  //
///////////////////////



///////////////
//  methods  //
///////////////

void importancesampling_set::vegasgrid_average(){
  Logger logger("importancesampling_set::vegasgrid_average " + name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<double> alpha_new(alpha.size());

  alpha_new[0] = .75 * alpha[0] + .25 * alpha[1];
  for (int j = 1; j < alpha.size() - 1; j++){
    alpha_new[j] = .5 * alpha[j] + .25 * alpha[j - 1] + .25 * alpha[j + 1];
  }
  alpha_new[alpha.size() - 1] = .75 * alpha[alpha.size() - 1] + .25 * alpha[alpha.size() - 2];
  if (alpha_new.size() == 1){alpha_new[0] = 1.;}

  alpha = alpha_new;
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void importancesampling_set::vegasgrid_calculation(){
  Logger logger("importancesampling_set::vegasgrid_calculation " + name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  /*
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
void_calculation(vector<double> & alpha, vector<double> & beta, vector<int> & cuts_channel, vector<int> & counts_channel, vector<double> & sum_channel_weight, vector<double> & sum_channel_weight2){
  Logger logger("vegasgrid_calculation");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  */

  vector<double> averaged_channel_weight(alpha.size(), 0.);
  vector<double> averaged_channel_weight2(alpha.size(), 0.);
  for (int j = 0; j < alpha.size(); j++){
    if (n_acc_channel[j] + n_rej_channel[j] != 0){
      averaged_channel_weight[j] = sum_channel_weight[j] / (n_acc_channel[j] + n_rej_channel[j]);
      averaged_channel_weight2[j] = sum_channel_weight2[j] / (n_acc_channel[j] + n_rej_channel[j]);
    }
    else {
      averaged_channel_weight[j] = 0.;
      averaged_channel_weight2[j] = 0.;
    }
  }
  vector<double> value(alpha.size(), 0.);
  for (int j = 0; j < averaged_channel_weight.size(); j++){
    //    if (channel_weight[j].size() != 0){
    value[j] = averaged_channel_weight[j]; // original version !!!
    //    double temp = (averaged_channel_weight2[j] - pow(averaged_channel_weight[j], 2) / (n_acc_channel[j] + n_rej_channel[j])) / (n_acc_channel[j] + n_rej_channel[j] - 1);
    //    if (temp > 0.){value[j] = sqrt(temp);} // works quite well !!!
    //    else {value[j] = 0.;}
    //    double temp = averaged_channel_weight2[j] - (averaged_channel_weight[j], 2) / channel_weight[j].size();
    //    if (temp > 0.){value[j] = pow(temp, 0.25);} // works quite well !!!
    //    else {value[j] = 0.;}
    //      value[j] = averaged_channel_weight2[j] - (averaged_channel_weight[j], 2) / channel_weight[j].size(); // works worse than with sqrt(...)
    //    }
    //    else {value[j] = 0.;}
  }
  double sum_value = accumulate(value.begin(), value.end(), 0.);
  for (int j = 0; j < averaged_channel_weight.size(); j++){
    logger << LOG_DEBUG_VERBOSE <<"value[" << setw(3) << j << "] = " << setprecision(8) << setw(15) << value[j] << "   averaged_channel_weight[" << setw(3) << j << "] = " << setprecision(8) << setw(15) << averaged_channel_weight[j] << "   n_acc_events = " << n_acc_channel[j] << "   n_rej_events = " << n_rej_channel[j] << endl;
  }
  for (int j = 0; j < alpha.size(); j++){value[j] = value[j] / sum_value;}
  
  int MCweight_lc_min = 0;
  vector<int> MCweights_set_to_min(alpha.size());
  double MCweight_limit_min = 1.e-3 / alpha.size();
  for (int j = 0; j < value.size(); j++){
    if (value[j] < MCweight_limit_min){
      MCweight_lc_min++;
      MCweights_set_to_min[j] = 1;
      value[j] = 0.;
    }
  }
  sum_value = accumulate(value.begin(), value.end(), 0.);
  double normalization = (1. - MCweight_lc_min * MCweight_limit_min) / sum_value;
  for (int j = 0; j < value.size(); j++){
    if (MCweights_set_to_min[j] == 1){value[j] = MCweight_limit_min;}
    else {value[j] = value[j] * normalization;}
  }
  
  alpha = value;


  //  vector<double> alpha_new(alpha.size());
  vegasgrid_average();
  //    for (int j = 0; j < beta.size(); j++){cout << "alpha[" << setw(3) << j << "] = " << setprecision(8) << setw(15) << showpoint << alpha[j] << "   alpha_new[" << setw(3) << j << "] = " << setprecision(8) << setw(15) << showpoint << alpha_new[j] << endl;}
  //  alpha = alpha_new;

  beta.resize(alpha.size());
  beta[0] = alpha[0];
  for (int i = 1; i < beta.size(); i++){
    beta[i] = beta[i - 1] + alpha[i];
  }
  
  //  for (int j = 0; j < beta.size(); j++){cout << "alpha[" << setw(3) << j << "] = " << setprecision(8) << setw(15) << showpoint << alpha[j] << "   beta[" << setw(3) << j << "] = " << setprecision(8) << setw(15) << beta[j] << "   " << setprecision(8) << setw(15) << sqrt(h_propto_pot(double(j) / beta.size(), 0, psi_exp_pdf)) * (2 * psi_E) << "   n_accepted = " << setw(8) << channel_weight[j].size() << "   n_rejected = " << setw(8) << n_rej_channel[j] << "   n_total = " << setw(8) << channel_weight[j].size() + n_rej_channel[j] << endl;}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void importancesampling_set::psp_IS_optimization(double & this_psp_weight, double & this_psp_weight2){
  Logger logger("importancesampling_set::psp_IS_optimization " + name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (active_optimization != 1){return;}
  if (end_optimization != 0){return;}

  n_acc_channel[channel]++;
  sum_channel_weight[channel] += abs(this_psp_weight);
  sum_channel_weight2[channel] += this_psp_weight2;
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void importancesampling_set::step_IS_optimization(int i_step_mode){
  Logger logger("importancesampling_set::step_IS_optimization " + name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (end_optimization != 0){return;}
  if (active_optimization != 1){return;}
  if (i_step_mode % n_event_per_step != 0 || i_step_mode == 0){return;}

  vegasgrid_calculation();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void importancesampling_set::result_IS_optimization(int i_step_mode){
  Logger logger("importancesampling_set::result_IS_optimization " + name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (active_optimization != 1){return;}
  if (end_optimization != 0){return;}

  if (i_step_mode != n_event_per_step * n_optimization_step){return;}
  
  end_optimization = 1;
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void importancesampling_set::output_IS_optimization(int i_step_mode){
  Logger logger("importancesampling::output_IS_optimization " + name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (active_optimization == -1){return;}
  if (end_optimization != 1){return;}
  
  //  if (i_step_mode != n_event_per_step * n_optimization_step){return;}
  
  ofstream out_alpha;
  out_alpha.open(filename.c_str(), ofstream::out | ofstream::trunc);
  for (int j = 0; j < alpha.size(); j++){
    out_alpha << setprecision(15) << alpha[j] << endl;
  }
  out_alpha.close();

  end_optimization = 2;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void importancesampling_set::readin_IS_optimization(){
  Logger logger("importancesampling::readin_MCweight_optimization " + name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_INFO << "filename_readin = " << filename_readin << endl;
  ifstream in_MCweights(filename_readin.c_str());
  vector<string> readin;
  char LineBuffer[128];
  while (in_MCweights.getline(LineBuffer, 128)){readin.push_back(LineBuffer);}
  in_MCweights.close();
  
  int start_data = 0;
  string start_line = "% " + name;
  //  string start_line_temp = "# " + name;

  // to guarantee compatibility with old weight files
  if (readin.size() > 0){
    if (readin[0][0] != '#'){start_data = 1;}
  }
  // end
  
  vector<string> data;
  for (int i = 0; i < readin.size(); i++){
    if (readin[i] == start_line){
      logger << LOG_INFO << "Data readin started for " << name << " in line " << i << endl;
      start_data = 1;
    }
    else if (readin[i][0] == '%' && start_data == 1){
      logger << LOG_INFO << "Data readin finished for " << name << " in line " << i << endl;
      break;
    }
    else if (readin[i][0] == '#' && start_data == 1){
      logger << LOG_INFO << "Commented line ignored for " << name << " in line " << i << ": " << readin[i] << endl;
    }
    else if (start_data){
      data.push_back(readin[i]);
    }
    else {logger << LOG_DEBUG << "readin[" << setw(4) << i << "] = " << readin[i] << " is ignored." << endl;}
  }
  
  logger << LOG_INFO << "Number of lines in " << filename_readin << ": " << readin.size() << endl;
  logger << LOG_INFO << "Number of relevant input lines for " << name << ": " << data.size() << "   (n_gridsize = " << n_gridsize << ")" << endl;
  
  if (data.size() != n_gridsize){
    logger << LOG_INFO << "Number of input lines does not fit n_gridsize value." << endl;
  }

  if (n_gridsize > 0){
    for (int i_c = 0; i_c < n_gridsize; i_c++){
      alpha[i_c] = atof(data[i_c].c_str());
      logger << LOG_DEBUG << "alpha[" << setw(4) << i_c << "] = " << alpha[i_c] << endl;
    }
    
    beta[0] = alpha[0];
    for (int i_c = 1; i_c < n_gridsize; i_c++){
      beta[i_c] = beta[i_c - 1] + alpha[i_c];
    }
  }
  else {
    logger << LOG_INFO << "Should not happen: n_gridsize = " << n_gridsize << endl;
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
