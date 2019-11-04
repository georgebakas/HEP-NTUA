#include "../include/classes.cxx"
////////////////////
//  constructors  //
////////////////////
multichannel_set::multichannel_set(){}

multichannel_set::multichannel_set(int _n_channel){  // should not be used at all !!!
  n_channel = _n_channel;
  n_alpha_events = 100000;
  n_alpha_steps = 10;
  switch_minimum_weight = 1;
  limit_minimum_weight = 0.0001;
  reserved_minimum_weight = 0.001;

  active_optimization = 0;
  end_optimization = 0; // !!!
  //  reserved_minimum_weight = 0.001; // option to change via input ???

  alpha.resize(n_channel, 1. / double (n_channel));
  beta.resize(n_channel, 0.);
  beta[0] = alpha[0];
  for (int i = 1; i < beta.size(); i++){beta[i] = beta[i - 1] + alpha[i];}

  n_acc_channel.resize(n_channel, 0);
  n_rej_channel.resize(n_channel, 0);
  
  g_channel.resize(n_channel, 0.);
  g_IS_channel.resize(n_channel, 0.);
  
  w_channel_step_sum.resize(n_channel, 0.);
  w_channel_full_sum.resize(n_channel, 0.);
  //  MC_sum_w_channel.resize(n_channel, 0.);
  // ???
  w_channel_av.resize(n_channel, 0.);

  alpha_it.resize(n_alpha_steps, vector<double> (n_channel, 0.));
  diff_w.resize(n_alpha_steps, 0.);

  i_alpha_it = 0;
  x_alpha_it_min = 0;
}



multichannel_set::multichannel_set(string _name, int _n_channel, int _n_alpha_events, int _n_alpha_steps, int _switch_minimum_weight, double _limit_minimum_weight, double _reserved_minimum_weight, int _switch_MC, string _filename, string _filename_readin){
  Logger logger("multichannel_set::multichannel_set " + _name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  name = _name;
  n_channel = _n_channel;
  n_alpha_events = _n_alpha_events;
  n_alpha_steps = _n_alpha_steps;

  // switch_minimum_weight == 0: Weights below limit_minimum_weight are set to zero (after last iteration step).
  // switch_minimum_weight == 1: Weights below limit_minimum_weight are set to limit_minimum_weight.
  switch_minimum_weight = _switch_minimum_weight;

  // limit_minimum_weight = relative minimum weight, wrt. average weight (1. / n_channel)
  limit_minimum_weight = _limit_minimum_weight;
  
  ///  reserved_minimum_weight = _reserved_minimum_weight;
  // Could be a better choice:

  // reserved_minimum_weight = absolute minimum weight (reserved are in total: n_below_minimum / n_channel)
  reserved_minimum_weight = limit_minimum_weight / n_channel;
  
  a_reserved_min = reserved_minimum_weight;

  // input parameter "_reserved_minimum_weight" is not needed any more.
  // also "a_reserved_min" is not needed -> always identical to "reserved_minimum_weight"
  
  if (_switch_MC == 0 || _switch_MC == 2){
    active_optimization = 0;
    end_optimization = 1;
  }
  else if (_switch_MC == 1 || _switch_MC == 3){
    active_optimization = 1;
    end_optimization = 0;
  }
  else if (_switch_MC == -1){
    active_optimization = -1;
    end_optimization = -1;
  }
  else {logger << LOG_FATAL << "No valid value of switch_MC chosen." << endl; exit(1);}
  
  filename = _filename;
  filename_readin = _filename_readin;
  
  if (n_channel == 1){end_optimization = 1;} // check if output has to be created immediately !!!

  counter_minimum_weight = 0;
  alpha.resize(n_channel, 1. / double (n_channel));
  beta.resize(n_channel, 0.);
  beta[0] = alpha[0];
  for (int i = 1; i < beta.size(); i++){beta[i] = beta[i - 1] + alpha[i];}

  n_gen = 0;
  n_acc = 0;
  n_rej = 0;
  n_tec = 0;
  n_nan = 0;

  n_acc_channel.resize(n_channel, 0);
  n_rej_channel.resize(n_channel, 0);
  
  channel = 0;

  g_channel.resize(n_channel, 0.);
  g_IS_channel.resize(n_channel, 0.);
  
  w_channel_step_sum.resize(n_channel, 0.);
  w_channel_full_sum.resize(n_channel, 0.);
  //  MC_sum_w_channel.resize(n_channel, 0.);
  // ???
  w_channel_av.resize(n_channel, 0.);

  alpha_it.resize(n_alpha_steps, vector<double> (n_channel, 0.));
  diff_w.resize(n_alpha_steps, 0.);

  i_alpha_it = 0;
  x_alpha_it_min = 0;

  if (_switch_MC == 2 || _switch_MC == 3){
    readin_MCweight_optimization();
  }
  
  logger << LOG_DEBUG << setw(30) << "n_channel" << " = " << n_channel << endl;
  logger << LOG_DEBUG << setw(30) << "n_alpha_events" << " = " << n_alpha_events << endl;
  logger << LOG_DEBUG << setw(30) << "n_alpha_steps" << " = " << n_alpha_steps << endl;
  logger << LOG_DEBUG << setw(30) << "switch_minimum_weight" << " = " << switch_minimum_weight << endl;
  logger << LOG_DEBUG << setw(30) << "limit_minimum_weight" << " = " << limit_minimum_weight << endl;
  logger << LOG_DEBUG << setw(30) << "reserved_minimum_weight" << " = " << reserved_minimum_weight << endl;
  logger << LOG_DEBUG << setw(30) << "active_optimization" << " = " << active_optimization << endl;
  logger << LOG_DEBUG << setw(30) << "end_optimization" << " = " << end_optimization << endl;


  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

///////////////////////
//  access elements  //
///////////////////////

///////////////
//  methods  //
///////////////

void multichannel_set::psp_MCweight_optimization(double & integrand, double & g_tot){
  Logger logger("multichannel_set::psp_MCweight_optimization " + name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (active_optimization != 1){return;}
  if (end_optimization != 0){return;}

  n_acc++;
  n_acc_channel[channel]++;
  double w_help = pow(abs(integrand), 2) / g_tot;
  if (munich_isnan(w_help)){// || w_help == 0.){
    logger << LOG_DEBUG << "munich_isnan(w_help)" << endl;
    logger << LOG_DEBUG << "integrand = " << integrand << endl;
    logger << LOG_DEBUG << "g_tot = " << g_tot << endl;
    for (int i_c = 0; i_c < n_channel; i_c++){logger << LOG_DEBUG << "g_channel[" << i_c << "] = " << g_channel[i_c] << endl;}
  }
  else {
    for (int i_c = 0; i_c < n_channel; i_c++){
      if (munich_isnan(g_channel[i_c])){ // || munich_isnan(w_help)) { // doubled "w_help" !!!
	logger << LOG_ERROR << "weight in channel " << i_c << " of " << name << " is nan; " << g_channel[i_c] << ", " << w_help << endl;
      }
      else {
	w_channel_step_sum[i_c] += g_channel[i_c] * w_help;
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void multichannel_set::step_MCweight_optimization(int i_step_mode){
  Logger logger("multichannel_set::step_MCweight_optimization " + name);
  logger << LOG_DEBUG << "called" << endl;

  if (end_optimization != 0){return;}
  if (active_optimization != 1){return;}
  if (i_step_mode > n_alpha_events * n_alpha_steps){return;}
  if (i_step_mode % n_alpha_events != 0){return;}
  
  for (int i_c = 0; i_c < n_channel; i_c++){
    w_channel_full_sum[i_c] += w_channel_step_sum[i_c];
  }
  w_channel_step_sum = vector<double> (n_channel, 0.);

  vector<double> c_w_channel_sum(n_channel);
  double c_w_norm;

  logger << LOG_INFO << "i_alpha_it = " << i_alpha_it << "   x_alpha_it_min = " << x_alpha_it_min << endl;
  
  alpha_it[i_alpha_it] = alpha;

  // evaluation of quality measure diff_w -> decides if new alpha's should be determined:
  for (int i_c = 0; i_c < n_channel; i_c++){
    if (n_acc_channel[i_c] + n_rej_channel[i_c] != 0){// Is this check needed ???
      //w_channel_av[i_c] = w_channel_full_sum[i_c] / (n_rej + n_acc);
      //w_channel_av[i_c] = w_channel_full_sum[i_c] / (i_rej + i_acc);
      w_channel_av[i_c] = w_channel_full_sum[i_c] / n_alpha_events; // w_channel_av and w_channel_full_sum are essentially the same - could be devided by something channel dependent instead.
      //    logger << LOG_DEBUG << "w_channel_av[" << i_c << "] = " << w_channel_av[i_c] << endl;
    }
    else {w_channel_av[i_c] = 0.;}
    
    //  maybe try different power of w_channel_av[i_c] instead of .5 !!!
    c_w_channel_sum[i_c] = alpha[i_c] * sqrt(w_channel_av[i_c]);
    //    c_w_channel_sum[i_c] = alpha[i_c] *pow(w_channel_av[i_c], exponent_new_alpha_set);
  }
  w_channel_full_sum = vector<double> (n_channel, 0.);
  sort(w_channel_av.begin(), w_channel_av.end());
  diff_w[i_alpha_it] = w_channel_av[w_channel_av.size() - 1] - w_channel_av[0];

  c_w_norm = accumulate(c_w_channel_sum.begin(), c_w_channel_sum.end(), 0.);
  assert(!munich_isnan(c_w_norm));

  x_minimum_diff_w = x_alpha_it_min;
  logger << LOG_INFO << "diff_w[" << setw(2) << x_alpha_it_min << "] = " << diff_w[x_alpha_it_min] << "   alpha[0] = " << setprecision(15) << setw(23) << alpha_it[x_alpha_it_min][0] << endl;

  for (int i_s = x_alpha_it_min + 1; i_s < i_alpha_it; i_s++){  //  to use only alpha sets after previous IS step; before:    //  for (int i_s = 0; i_s < i_alpha_it; i_s++){
    logger << LOG_INFO << "diff_w[" << setw(2) << i_s << "] = " << setprecision(15) << setw(23) << diff_w[i_s] << "   alpha[0] = " << setprecision(15) << setw(23) << alpha_it[i_s][0] << endl;
    if (diff_w[i_s] < diff_w[x_minimum_diff_w]){
      x_minimum_diff_w = i_s;
    }
  }
  double diff_w_min = diff_w[x_minimum_diff_w];
  logger << LOG_INFO << "diff_w_min = " << setprecision(15) << setw(23) << diff_w_min << " = diff_w[" << setw(2) << x_minimum_diff_w << "]" << endl;
  logger << LOG_INFO << "diff_w[" << setw(2) << i_alpha_it << "] = " << setprecision(15) << setw(23) << diff_w[i_alpha_it] << "   alpha[0] = " << setprecision(15) << setw(23) << alpha_it[i_alpha_it][0] << " --- newly determined" << endl; //" --- newly determined" << endl;

  /*
  if (name == "MC_phasespace"){
    for (int i_s = 0; i_s < i_alpha_it + 1; i_s++){
      if (i_s == x_alpha_it_min){logger << LOG_INFO << "------------------------------------------------------------------------------------------------------------------" << endl;}
      logger << LOG_INFO << "diff_w[" << setw(2) << i_s << "]           = " << setprecision(15) << setw(23) << log10(diff_w[i_s]) << "     " << setprecision(15) << setw(23) << alpha_it[i_s][0] << endl;
      if (i_s == i_alpha_it){logger << LOG_INFO << "------------------------------------------------------------------------------------------------------------------" << endl;}
    }
    
    logger << LOG_INFO << setw(20) << "x_alpha_it_min" << " = " << x_alpha_it_min << endl;
    logger << LOG_INFO << setw(20) << "i_alpha_it" << " = " << i_alpha_it << endl;
    logger << LOG_INFO << setw(20) << "diff_w_min" << " = " << setprecision(15) << setw(23) << log10(diff_w_min) << endl;
    logger << LOG_INFO << setw(20) << "x_minimum_diff_w" << " = " << x_minimum_diff_w << endl;

    if ((i_alpha_it == 0) || (diff_w[i_alpha_it] <= diff_w_min)){logger << LOG_INFO << setw(30) << "" << "New set of alpha's (-> " << i_alpha_it + 1 << ") is used in next step." << endl;}
    else {logger << LOG_INFO << setw(30) << "" << "Set of alpha's no. " << x_minimum_diff_w << " is used in next step." << endl;}
    logger.newLine(LOG_INFO);
  }
  */
  
  logger << LOG_INFO << "limit_minimum_weight = " << limit_minimum_weight << endl;
  logger << LOG_INFO << "switch_minimum_weight = " << switch_minimum_weight << endl;
  logger << LOG_INFO << "reserved_minimum_weight = " << reserved_minimum_weight << endl;
  
  // check if i_alpha_it == x_alpha_it_min makes more sense here !!!
  if ((i_alpha_it == 0) || (diff_w[i_alpha_it] <= diff_w_min)){
    // determination of new alpha's for next iteration step:
    for (int i_c = 0; i_c < n_channel; i_c++){
      alpha[i_c] = c_w_channel_sum[i_c] / c_w_norm;
    }

    // new possibly iterative determination of weights after application of minimum_weight

    counter_minimum_weight = 0;
    no_channel_minimum_weight.clear();
    vector<int> new_no_channel_minimum_weight(1);
    int step_counter = 0;
    while (new_no_channel_minimum_weight.size() > 0){
      step_counter++;
      new_no_channel_minimum_weight.clear();
     
      for (int i_c = 0; i_c < n_channel; i_c++){
	if (alpha[i_c] <= reserved_minimum_weight){
	  int flag = 0;
	  for (int j_c = 0; j_c < no_channel_minimum_weight.size(); j_c++){
	    if (i_c == no_channel_minimum_weight[j_c]){flag = 1; break;}
	  }
	  if (!flag){	  
	    alpha[i_c] = 0.;
	    counter_minimum_weight++;
	    new_no_channel_minimum_weight.push_back(i_c);  // always: set weights to zero only in the result step !!!
	  }
	}
      }

      double a_norm;
      a_norm = accumulate(alpha.begin(), alpha.end(), 0.) / (1. - counter_minimum_weight * reserved_minimum_weight);
      for (int i_c = 0; i_c < n_channel; i_c++){
	if (alpha[i_c] != 0.){
	  alpha[i_c] = alpha[i_c] / a_norm;
	}
      }

      logger << LOG_INFO << "limit_minimum_weight  step no. " << setw(3) << step_counter << " : counter_minimum_weight = " << setw(3) << counter_minimum_weight << "   a_norm = " << a_norm << endl;
	
      for (int j_c = 0; j_c < new_no_channel_minimum_weight.size(); j_c++){
	no_channel_minimum_weight.push_back(new_no_channel_minimum_weight[j_c]);
      }
    }
    for (int i_m = 0; i_m < no_channel_minimum_weight.size(); i_m++){
      alpha[no_channel_minimum_weight[i_m]] = reserved_minimum_weight;
    }
  }
  else {
    // Try without '+1' (often falls back to original set, but might deliver better convergence) ???
    /*
    alpha = alpha_it[x_minimum_diff_w]; // why "+ 1" ??? check !!!
    beta[0] = alpha[0];
    for (int i = 1; i < beta.size(); i++){beta[i] = beta[i - 1] + alpha[i];}
    */
    alpha = alpha_it[x_minimum_diff_w + 1]; // why "+ 1" ??? check !!!
    //  alpha = alpha_it[no_diff_w_min + 1]; // why "+ 1" ??? check !!!

    // In particular for resumed runs:
    counter_minimum_weight = 0;
    no_channel_minimum_weight.clear();
    for (int i_c = 0; i_c < n_channel; i_c++){
      if (alpha[i_c] <= reserved_minimum_weight){
	counter_minimum_weight++;
	no_channel_minimum_weight.push_back(i_c);  // always: set weights to zero only in the result step !!!
      }
    }    
    
    logger << LOG_INFO << "optimization step doesn't seem to improve the set of weights:" << endl;
    logger << LOG_INFO << "alpha set no. " << x_minimum_diff_w + 1 << "  is used in the next step." << endl;
    logger << LOG_INFO << setw(40) << "counter_minimum_weight" << " = " << counter_minimum_weight << endl;
    logger << LOG_INFO << setw(40) << "no_channel_minimum_weight.size()" << " = " << no_channel_minimum_weight.size() << endl;
    for (int j_c = 0; j_c < no_channel_minimum_weight.size(); j_c++){
      logger << LOG_INFO << "minimum:   alpha[" << no_channel_minimum_weight[j_c] << "] = " << alpha[no_channel_minimum_weight[j_c]] << endl;
    }
  }
  
  beta[0] = alpha[0];
  for (int i = 1; i < beta.size(); i++){beta[i] = beta[i - 1] + alpha[i];}

  i_alpha_it++;
  
  for (int i_c = 0; i_c < n_channel; i_c++){
    if (alpha[i_c] == 0.){g_channel[i_c] = 0.;}
  }

  logger << LOG_INFO << "accumulate(alpha.begin(), alpha.end(), 0.) = " << accumulate(alpha.begin(), alpha.end(), 0.) << endl;
  logger << LOG_INFO << "alpha[" << 0 << "] = " << setprecision(15) << setw(23) << alpha[0] << "   beta[" << 0 << "] = " << setprecision(15) << setw(23) << beta[0] << endl;
  for (int i_c = 0; i_c < n_channel; i_c++){
    if (alpha[i_c] <= reserved_minimum_weight){
      logger << LOG_INFO << "alpha[" << i_c << "] = " << setprecision(15) << setw(23) << alpha[i_c] << "   beta[" << i_c << "] = " << setprecision(15) << setw(23) << beta[i_c] << endl;
    }
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}



void multichannel_set::result_MCweight_optimization(int i_step_mode){
  Logger logger("multichannel_set::result_MCweight_optimization " + name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_INFO << "i_step_mode = " << i_step_mode << endl;
  logger << LOG_INFO << "end_optimization = " << end_optimization << endl;

  if (active_optimization != 1){return;}
  if (end_optimization != 0){return;}
  if (i_step_mode != n_alpha_events * n_alpha_steps){return;}
  
  // changed from x_minimum_diff_w = 0, for (int j = 0... !!!
  //  x_minimum_diff_w = 0;
  //  x_minimum_diff_w = 1;
  // use only alpha sets after previous IS step:
  x_minimum_diff_w = x_alpha_it_min;
  
  //  for (int i_s = 0; i_s < diff_w.size(); i_s++){
  for (int i_s = x_alpha_it_min + 1; i_s < diff_w.size(); i_s++){  // to use only alpha sets after previous IS step; before:   //  for (int i_s = 1; i_s < diff_w.size(); i_s++){
    if (diff_w[i_s] < diff_w[x_minimum_diff_w]){
      x_minimum_diff_w = i_s;
    }
  }
  alpha = alpha_it[x_minimum_diff_w];
  logger << LOG_INFO << "alpha set no. " << x_minimum_diff_w << "  is used in the integration." << endl;

  logger << LOG_INFO << "limit_minimum_weight = " << limit_minimum_weight << endl;
  logger << LOG_INFO << "switch_minimum_weight = " << switch_minimum_weight << endl;
  logger << LOG_INFO << "reserved_minimum_weight = " << reserved_minimum_weight << endl;

  counter_minimum_weight = 0;
  no_channel_minimum_weight.clear();
  for (int i_c = 0; i_c < n_channel; i_c++){
    if (alpha[i_c] == reserved_minimum_weight){
      counter_minimum_weight++;
      no_channel_minimum_weight.push_back(i_c);
    }
  }

  if (switch_minimum_weight == 0){
    for (int j_c = 0; j_c < no_channel_minimum_weight.size(); j_c++){
      alpha[no_channel_minimum_weight[j_c]] = 0.;
    }
    double a_norm;
    a_norm = accumulate(alpha.begin(), alpha.end(), 0.);
    for (int i_c = 0; i_c < n_channel; i_c++){
      if (alpha[i_c] != 0.){
	alpha[i_c] = alpha[i_c] / a_norm;
      }
    }
  }

  beta[0] = alpha[0];
  for (int i = 1; i < beta.size(); i++){beta[i] = beta[i - 1] + alpha[i];}

  end_optimization = 1;

  logger << LOG_INFO << "counter_minimum_weight = " << counter_minimum_weight << endl;
  logger << LOG_INFO << "no_channel_minimum_weight.size() = " << no_channel_minimum_weight.size() << endl;
  logger << LOG_INFO << "accumulate(alpha.begin(), alpha.end(), 0.) = " << accumulate(alpha.begin(), alpha.end(), 0.) << endl;

  
  if (switch_minimum_weight == 0){
    logger << LOG_INFO << "The following alpha's are set to  0 :" << endl;
  }
  else if (switch_minimum_weight == 0){
    logger << LOG_INFO << "The following alpha's are set to  " << reserved_minimum_weight << " :" << endl;
  }
  for (int j_c = 0; j_c < no_channel_minimum_weight.size(); j_c++){
    logger << LOG_INFO << "alpha[" << no_channel_minimum_weight[j_c] << "] = " << alpha[no_channel_minimum_weight[j_c]] << endl;
  }
 
  /*
  for (int i_s = 0; i_s < alpha_it.size(); i_s++){
    logger << LOG_DEBUG << "diff_w[" << setw(2) << i_s << "] = " << setprecision(15) << setw(23) << diff_w[i_s] << endl;
  }
  */
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void multichannel_set::output_MCweight_optimization(int i_step_mode, string & filename_MCweight){
  Logger logger("multichannel::output_MCweight_optimization " + name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (active_optimization == -1){return;}
  if (end_optimization != 1){return;}
  //  if (i_step_mode != n_alpha_events * n_alpha_steps){return;}

  ofstream out_MCweights;
  out_MCweights.open(filename_MCweight.c_str(), ofstream::out | ofstream::app);
  out_MCweights << "% " << name << endl;
  for (int i_c = 0; i_c < n_channel; i_c++){out_MCweights << setprecision(16) << alpha[i_c] << endl;}
  out_MCweights.close();

  end_optimization = 2;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void multichannel_set::readin_MCweight_optimization(){
  Logger logger("multichannel::readin_MCweight_optimization " + name);
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
  for (int i = 0; i < readin.size(); i++){
    for (int j = 0; j < readin[i].size(); j++){
      if (readin[i][j] == '#'){readin[i][j] = '%';}
    }
  }
  logger << LOG_INFO << name << endl;
  logger << LOG_INFO << readin[0] << endl;
  //  if (name == "MC_phasespace"){
  if (readin[0] == "% MC"){
    readin[0] = "% MC_phasespace";
    logger << LOG_INFO << readin[0] << endl;
  }
  // end
  
  vector<string> data;
  for (int i = 0; i < readin.size(); i++){
    logger << LOG_INFO << "***" << readin[i] << "***" << endl;
    if (readin[i] == start_line){
      logger << LOG_INFO << "Data readin started for " << name << " in line " << i << endl;
      start_data = 1;
    }
    else if (readin[i][0] == '%' && start_data == 1){
      logger << LOG_INFO << "Data readin finished for " << name << " in line " << i << endl;
      //	readin.erase(readin.begin() + i);
      break;
    }
    else if (readin[i][0] == '#' && start_data == 1){
      logger << LOG_INFO << "Commented line ignored for " << name << " in line " << i << ": " << readin[i] << endl;
    }
    else if (start_data){
      //      logger << LOG_INFO << readin[i] << endl;
      data.push_back(readin[i]);
    }
    else {logger << LOG_DEBUG << "readin[" << setw(4) << i << "] = " << readin[i] << " is ignored." << endl;}
  }
  logger << LOG_INFO << "Number of lines in " << filename_readin << ": " << readin.size() << endl;
  logger << LOG_INFO << "Number of relevant input lines for " << name << ": " << data.size() << "   (n_channel = " << n_channel << ")" << endl;
  
  if (data.size() != n_channel){
    logger << LOG_INFO << "Number of input lines does not fit n_channel value." << endl;
  }

  for (int i_c = 0; i_c < n_channel; i_c++){
    alpha[i_c] = atof(data[i_c].c_str());
    logger << LOG_DEBUG << "alpha[" << setw(4) << i_c << "] = " << alpha[i_c] << endl;
  }

  beta[0] = alpha[0];
  for (int i_c = 1; i_c < beta.size(); i_c++){
    beta[i_c] = beta[i_c - 1] + alpha[i_c];
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
