#include "../include/classes.cxx"

double ran_temp(vector<double> & s){
  static Logger logger("ran_temp");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
  // TODO: maybe do something more clever here.. e.g. make s static and don't bother with the random seed anywhere else
  // problem: how to seed initially? simply incrementing zwahl times will lead to a periodicity if there are too many parallel computations

  assert(s[0]+s[1]+s[2]>0);
  assert(s[0]>=0);
  assert(s[0]<=1);
  assert(s[1]>=0);
  assert(s[1]<=1);
  assert(s[2]>=0);
  assert(s[2]<=1);

  s[0] = fmod(s[0] + s[1] + s[2], 1.);
  s[1] = fmod(s[0] + s[1] + s[2], 1.);
  s[2] = fmod(s[0] + s[1] + s[2], 1.);

//  return rand()/double(RAND_MAX);

  logger << LOG_DEBUG_VERBOSE << "finished - s[0] = " << s[0] << endl;
  return s[0];
}



randomvariable::randomvariable(int _n_events, int _switch_IS_MC, int _n_steps, int _n_bins, const string &_name) {
  static Logger logger("randomvariable::randomvariable");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  logger << LOG_DEBUG << "initialize " << _name << "(" << _n_events << ", " << _switch_IS_MC << ", " << _n_steps << ", " << _n_bins << ")" << endl;

  init(_n_events, _switch_IS_MC, _n_steps, _n_bins, _name);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void randomvariable::init(int _n_events, int _switch_IS_MC, int _n_steps, int _n_bins, const string & _name) {
  static Logger logger("randomvariable::init");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  n_events = _n_events;
  switch_IS_MC = _switch_IS_MC;
  n_steps = _n_steps;
  n_bins = _n_bins;
  name = _name;

  if (switch_IS_MC == 0 && n_bins == 0) {
    logger << LOG_DEBUG << "z1z2 input incomplete, but optimization turned off anyway" << endl;
    n_bins = 1;
  }

  assert(n_bins > 0);

  alpha.resize(n_bins, 1. / double(n_bins));
  beta.resize(n_bins);

  if (switch_IS_MC == 0 || switch_IS_MC == 2){end_optimization = 1;}
  else if (switch_IS_MC == 1 || switch_IS_MC == 3){end_optimization = 0;}
  //  if (switch_IS_MC == 0 || switch_IS_MC == 2){opt_end = 1;}
  //  else if (switch_IS_MC == 1 || switch_IS_MC == 3){opt_end = 0;}

  channel=0;

  // imp_sampling will be useful when complete random-number organization is done vai randommanager !!!
  // should then be used as an input paramter again !!!
  imp_sampling = true;

  used = false;

  beta[0] = alpha[0];
  for (int i_b = 1; i_b < n_bins; i_b++){beta[i_b] = beta[i_b - 1] + alpha[i_b];}

  n_acc_channel.resize(n_bins, 0);
  n_rej_channel.resize(n_bins, 0);
  sum_channel_weight.resize(n_bins, 0.);
  sum_channel_weight2.resize(n_bins, 0.);

  s.resize(3);
  s[0]=0;
  s[1]=0;
  s[2]=0;

  manager=NULL;

  proceeding_save();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void randomvariable::vegasgrid_average(){
  Logger logger("randomvariable::vegasgrid_average " + name);
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<double> alpha_new(n_bins);

  alpha_new[0] = .75 * alpha[0] + .25 * alpha[1];
  for (int j = 1; j < n_bins - 1; j++){
    alpha_new[j] = .5 * alpha[j] + .25 * alpha[j - 1] + .25 * alpha[j + 1];
  }
  alpha_new[n_bins - 1] = .75 * alpha[n_bins - 1] + .25 * alpha[n_bins - 2];
  if (alpha_new.size() == 1){alpha_new[0] = 1.;}

  alpha = alpha_new;
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void randomvariable::psp_IS_optimization(double & psp_weight) {
  static Logger logger("randomvariable::psp_IS_optimization");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  //  if (!imp_sampling || !used){return;}
  if (!imp_sampling){return;}
  if (!used){return;}
  if (end_optimization){return;}
  
  // imp_sampling ??? used ???
  // TODO: stop recording if opt_end==1? -> yes !!!

  n_acc_channel[channel]++;
  sum_channel_weight[channel] += fabs(psp_weight);
  //  alternative optimization -> check !!!
  //  sum_channel_weight2[channel] += pow(psp_weight, 2);

  used = false;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void randomvariable::step_IS_optimization(int i_step_mode){
  static Logger logger("randomvariable::step_IS_optimization");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  if (end_optimization){return;}

  if (i_step_mode % n_events == 0 && i_step_mode != 0){

    //    step_IS_optimization_temp();
    
    for (int j = 0; j < n_bins; j++){
      if (n_acc_channel[j] + n_rej_channel[j] != 0){
	alpha[j] = sum_channel_weight[j] / (n_acc_channel[j] + n_rej_channel[j]);
      }
      else {
	alpha[j] = 0.;
      }
    }
    
    double sum_value = accumulate(alpha.begin(), alpha.end(), 0.);
    
    if (sum_value == 0.){logger << LOG_DEBUG << "1st: sum_value = 0." << endl;}
    if (sum_value == 0.){
      for (int j = 0; j < n_bins; j++){
	alpha[j] = 1. / n_bins;
	sum_value = 1.;
      }
    }
    
    for (int j = 0; j < n_bins; j++){alpha[j] = alpha[j] / sum_value;}
    
    int MCweight_lc_min = 0;
    vector<int> MCweights_set_to_min(n_bins);
    
    double MCweight_limit_min = 1.e-2 / n_bins;
    sum_value = 0.;
    for (int j = 0; j < n_bins; j++){
      if (alpha[j] < MCweight_limit_min){
	MCweight_lc_min++;
	MCweights_set_to_min[j] = 1;
      }
      else {
	sum_value += alpha[j];
      }
    }
    if (sum_value == 0.){logger << LOG_DEBUG << "2nd: sum_value = 0." << endl;}

    double normalization = (1. - MCweight_lc_min * MCweight_limit_min) / sum_value;
    if (normalization != normalization){logger << LOG_WARN << "normalization != normalization" << endl;}
    if (normalization > 1.e200){logger << LOG_WARN << "normalization = +-inf" << endl;}
    
    for (int j = 0; j < n_bins; j++){
      if (MCweights_set_to_min[j] == 1){alpha[j] = MCweight_limit_min;}
      else {alpha[j] = alpha[j] * normalization;}
    }
    
    vegasgrid_average();
    //  vegasgrid_average_temp(alpha);
    
    beta[0] = alpha[0];
    for (int i = 1; i < n_bins; i++){beta[i] = beta[i - 1] + alpha[i];}
    
    //    vegasgrid_calculation_temp(alpha, beta, n_rej_channel, n_acc_channel, sum_channel_weight);
    
    //    logger << LOG_DEBUG << left << setw(50) << name << ": vegasgrid_calculation @ " << i_step_mode << " events" << "(" << accumulate(n_acc_channel.begin(),n_acc_channel.end(),0)+accumulate(n_rej_channel.begin(),n_rej_channel.end(),0) << " events used)" << "; variation is " << calc_variation() << endl;
    
    // shifted -> result_IS_optimization !!!
    result_IS_optimization(i_step_mode);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void randomvariable::result_IS_optimization(int i_step_mode) {
  static Logger logger("randomvariable::result_IS_optimization");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  if (end_optimization){return;}
  
  // Shifted here from step_IS_optimization:
  if (i_step_mode == n_events * n_steps){
    logger << LOG_DEBUG << "sampled random variable " << setw(50) << left << name << ":" << endl;
    for (int is = 0; is < n_bins; is++){
      logger << LOG_DEBUG << "alpha[" << right << setw(4) << is << "] = " << left << setw(15) << setprecision(8) << alpha[is] << "   beta[" << right << setw(4) << is << "] = " << left << setw(15) << setprecision(8) << beta[is] << endl;
    }
    end_optimization = 1;
    
    //      logger << LOG_DEBUG << setw(50) << left << name << ": variation is " << calc_variation() << endl;
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}













void randomvariable::get_random(double &x, double &g_IS, vector<double> &s) {
  static Logger logger("randomvariable::get_random(x, g_IS, s)");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
  if (used && !end_optimization) {
    logger << LOG_FATAL  << "randomvariable::get_random  (" << name << " used && !end_optimization) should never be called with used == true !" << endl;
    exit(1);
  }

  // At the moment, this simply overwrites the random number from psi.r[...] !!!
  // Later: only generate the random number here, but use global 's' !!!
  double random = ran_temp(s);

  // Later: If 'imp_sampling == false', this random numner is not modified by IS:
  if (!imp_sampling) {
    x = random;
    g_IS = 1.0;
    return;
  }

  // used refers only to IS sampling:
  used = true;

  // 'random_channel' is used to select IS bin:
  double random_channel = ran_temp(s);

  /*
  channel = 0;
  for (int i_b = n_bins / 2; i_b > 0; i_b /= 2){
    if (random_channel > beta[channel + i_b - 1]){channel += i_b;}
    
    if (n_bins != 64){
      //      logger << LOG_WARN << "i_b = " << i_b << "   channel = " << channel << "   random_channel = " << random_channel << " > " << beta[channel + i_b - 1] << " =  beta[" << channel + i_b - 1 << "]" << endl;
    }
    
    if (i_b == 3){
      logger << LOG_WARN << "i_b == 3 - before channel = " << channel << endl;
      for (int j_b = 0; j_b < 4; j_b++){
	logger << LOG_WARN << "random_channel = " << random_channel << " <= " << beta[channel + j_b] << " = beta[" << channel + j_b << "]" << endl;
	if (random_channel <= beta[channel + j_b]){
	  channel += j_b;
	  logger << LOG_WARN << "i_b == 3 - channel = " << channel << endl;
	  i_b = 0;
	  break;
	}
      }
    }
  }
  int temp_channel = channel;
  channel = 0;
  if (n_bins != 64){
    // More efficient: replace by nested intervals ???
    for (int i_b = 0; i_b < beta.size(); i_b++){
      if (random_channel <= beta[i_b]){
	channel = i_b;
	break;
      }
    }
    if (temp_channel != channel){
      logger << LOG_WARN << "temp_channel = " << temp_channel << " ! unequal = " << channel << " = channel   n_bins = " << n_bins << "   beta.size() = " << beta.size() << endl;
    }
  }
  */
  /*
  if (new_channel != channel){
    logger << LOG_WARN << "new_channel = " << new_channel << " != " << channel << " = channel" << endl;;
  }
  */
  channel = 0;
  if (n_bins == 64){
    for (int i_b = n_bins / 2; i_b > 0; i_b /= 2){
      if (random_channel > beta[channel + i_b - 1]){channel += i_b;}
    }
  }
  else {
    /*
    for (int i_b = 0; i_b < beta.size(); i_b++){
      if (random_channel <= beta[i_b]){
	channel = i_b;
	break;
      }
    }
    */
    /*
    int old_channel = 0;
    for (int i_b = 0; i_b < beta.size(); i_b++){
      if (random_channel <= beta[i_b]){
	old_channel = i_b;
	break;
      }
    }
    logger << LOG_DEBUG_VERBOSE << "old determination:  old_channel = " << old_channel << endl;
    */
    int extra = 0;
    for (int i_b = n_bins / 2; i_b > 0; i_b /= 2){
      //      logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << "   channel = " << channel << "   random_channel = " << random_channel << " > " << beta[channel + i_b - 1] << " =  beta[" << channel + i_b - 1 << "]" << endl;
      if (random_channel > beta[channel + i_b - 1]){channel += i_b;}
      //      logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << "   channel = " << channel << endl;
      if (i_b == 3){
	//	logger << LOG_DEBUG_VERBOSE << "i_b == 3 - before channel = " << channel << endl;
	for (int j_b = 0; j_b < 3 + extra; j_b++){
	  //	  logger << LOG_DEBUG_VERBOSE << "random_channel = " << random_channel << " <= " << beta[channel + j_b] << " = beta[" << channel + j_b << "]" << endl;
	  if (random_channel <= beta[channel + j_b]){
	    channel += j_b;
	    //	    logger << LOG_DEBUG_VERBOSE << "i_b == 3 - channel = " << channel << endl;
	    i_b = 0;
	    break;
	  }
	}
      }
      if (i_b % 2){extra++;}
    }
    /*
    logger << LOG_DEBUG_VERBOSE << "new determination:  channel = " << channel << endl;
    if (old_channel != channel){
      logger << LOG_WARN << "old_channel = " << old_channel << " !!! unequal !!! " << channel << " = channel" << endl;;
    }
    */
  }
  
  // 'x' replaces the flatly generated random number:
  //  x = (double(channel) + random) / double(beta.size());
  x = (double(channel) + random) / n_bins;

  // 'g_IS' is actually not needed here yet (only for 'psi.g_tot' determination) !!!
  // 'g_IS' is the additional contribution to the IS:
  g_IS = alpha[channel] * n_bins;

  //  logger << LOG_DEBUG_VERBOSE << name << ": " << channel << ", " << beta.size() << endl;
  //  logger << LOG_DEBUG_VERBOSE << setw(40) << name << " rch = " << setw(15) << setprecision(8) << random_channel << " r = " << setw(15) << setprecision(8) << random << " ch = " << setw(2) << channel << " alpha.size() = " << setw(2) << alpha.size() << "   a[ch] = " << setw(15) << setprecision(8) << alpha[channel] << "   g_IS = " << setw(15) << setprecision(8) << g_IS << endl;
  //  logger << LOG_DEBUG_VERBOSE << random_channel << ", " << random << ", " << channel << endl;
 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randomvariable::get_random(double &x, double &g_IS) {
  static Logger logger("randomvariable::get_random(x, g_IS)");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  // Only calls function with additional argument 's' (which is a class member anyway)
  // -> redunant ??? Should be merged !!!
  get_random(x, g_IS, s);

  if (!end_optimization) {
    // It tells manager that this variable has been used at this phase-space point (can be used later):
    if (manager != NULL){
      manager->add_var_to_queue(this);
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randomvariable::get_g_IS(double r, double & g_IS) {
  static Logger logger("randomvariable::get_g_IS");
  //  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  int temp_channel = (int)(r * n_bins);
  g_IS = alpha[temp_channel] * n_bins;
  
  if (r < 0. || r > 1. || munich_isnan(r)){
    temp_channel = 0;
    g_IS = sqrt(-1.);
  }
  
  
  //  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
/*
void randomvariable::get_g_IS(double r, double & g_IS) {
  //  static Logger logger("randomvariable::get_g_IS");
  //  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  int temp_channel = (int)(r * beta.size());
  //  int temp_channel = (int)(r * (double)beta.size());
  g_IS = alpha[temp_channel] * alpha.size();

  if (r < 0. || r > 1. || munich_isnan(r)){
    temp_channel = 0;
    g_IS = sqrt(-1.);
    //    logger << LOG_DEBUG << setw(40) << name << " r = " << setw(15) << setprecision(8) << r << " ch = " << setw(2) << channel << " alpha.size() = " << setw(2) << alpha.size() << "   a[ch] = " << setw(15) << setprecision(8) << alpha[channel] << "   g_IS = " << setw(15) << setprecision(8) << g_IS << endl;
  }
  
  //  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/
/*
void randomvariable::get_g_IS(double r, double & g_IS) {
  static Logger logger("randomvariable::get_g_IS");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  int temp_channel = (int)(r * (double)beta.size());
  g_IS = alpha[temp_channel] * alpha.size();

  if (r < 0. || r > 1. || munich_isnan(r)){
    temp_channel = 0;
    g_IS = sqrt(-1.);
    logger << LOG_DEBUG << setw(40) << name << " r = " << setw(15) << setprecision(8) << r << " ch = " << setw(2) << channel << " alpha.size() = " << setw(2) << alpha.size() << "   a[ch] = " << setw(15) << setprecision(8) << alpha[channel] << "   g_IS = " << setw(15) << setprecision(8) << g_IS << endl;
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/


void randomvariable::proceeding_out(ofstream &out_proceeding) {
  static Logger logger("randomvariable::proceeding_out");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  out_proceeding << "# " << name << endl;

  // TODO: handle imp_sampling!
  out_proceeding << save_end_optimization << endl;

  // Check what should happen dpendent on 'end_optimization' !!!
  // Maybe introduce active_optimization' to avoid output for production runs !!!
  if (end_optimization != -1){
    // There should actually not be a seed set for each variable !!! -> randommanager
    // SK: better only once for randommanager ???
    for (int i_r = 0; i_r < 3; i_r++){out_proceeding << double2hexastr(save_s[i_r]) << endl;}

    out_proceeding << "# " << name << " alpha (" << n_bins << ")" << endl;

    for (int i_b = 0; i_b < n_bins; i_b++){
      out_proceeding << double2hexastr(save_alpha[i_b]) << endl;
    }
  }

  if (end_optimization == 0){
    out_proceeding << "# " << name << " n_acc_channel - n_rej_channel - sum_channel_weight (3 x " << n_bins << ")" << endl;
    for (int i_b = 0; i_b < n_bins; i_b++){
      out_proceeding << save_n_acc_channel[i_b] << endl;
      out_proceeding << save_n_rej_channel[i_b] << endl;
      // Maybe switch to 'double2hexastr(...)' output !!!
      out_proceeding << double2hexastr(save_sum_channel_weight[i_b]) << endl;
      // ???
      //      out_proceeding << setprecision(16) << save_sum_channel_weight2[i_b] << endl;
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void randomvariable::proceeding_in(int &proc, vector<string> &readin) {
  static Logger logger("randomvariable::proceeding_in");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  end_optimization = atoi(readin[proc].c_str());
  proc += 1;
  if (end_optimization != -1){

    for (int i = 0; i < 3; i++){s[i] = hexastr2double(readin[proc + i]);}
    proc += 3;
    //    for (int i = 0; i < 3; i++){logger << LOG_DEBUG << "s[" << i << "] = " << s[i] << endl;}

    for (int i = 0; i < n_bins; i++){
      alpha[i] = hexastr2double(readin[proc + 1 * i]);
    }
    proc += n_bins;

    beta[0] = alpha[0];
    for (int i = 1; i < n_bins; i++){beta[i] = beta[i - 1] + alpha[i];}

    logger << LOG_DEBUG << name << endl;
    /*
    for (int i = 0; i < 3; i++){
      logger << LOG_DEBUG << "s[" << i << "] = " << s[i] << endl;
    }
    for (int i = 0; i < n_bins; i++){
      logger << LOG_DEBUG << "alpha[" << setw(4) << i << "] = " << setprecision(15) << setw(23) << alpha[i] << "   beta[" << setw(4) << i << "] = " << setprecision(15) << setw(23) << beta[i] << endl;
    }
    */
  }
  if (end_optimization == 0){
    for (int i = 0; i < n_bins; i++){
      n_acc_channel[i] = atoi(readin[proc + 3 * i].c_str());
      n_rej_channel[i] = atoi(readin[proc + 3 * i + 1].c_str());
      // Maybe switch to 'hexastr2double(...)' input (with corresponding input) !!!
      sum_channel_weight[i] = hexastr2double(readin[proc + 3 * i + 2]);
    }
    proc += 3 * n_bins;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randomvariable::check_proceeding_in(int & int_end, int & temp_check_size, vector<string> & readin){
  static Logger logger("randomvariable::check_proceeding_in");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // check if usage of 'int_end' makes sense here !!!
  
  end_optimization = atoi(readin[temp_check_size].c_str());
  temp_check_size += 1;
  logger << LOG_DEBUG << "end_optimization = " << end_optimization << endl;

  if (temp_check_size > readin.size()){
    int_end = 2;
    logger << LOG_DEBUG << "int_end = 2   temp_check_size = " << temp_check_size << " > " << readin.size() << " = readin.size()" << endl;
    return;
  }
  else {
    logger << LOG_DEBUG << "temp_check_size = " << temp_check_size << " <= " << readin.size() << " = readin.size()" << endl;
  }
  
  if (end_optimization != -1){
    //  Should be removed accordingly !!!
    //  s(3)
    temp_check_size += 3;

    if (temp_check_size > readin.size()){
      int_end = 2;
      logger << LOG_DEBUG << "int_end = 2   temp_check_size = " << temp_check_size << " > " << readin.size() << " = readin.size()" << endl;
      return;
    }
    else{
      logger << LOG_DEBUG << "temp_check_size = " << temp_check_size << " <= " << readin.size() << " = readin.size()" << endl;
    }
    
    temp_check_size += n_bins;
    //  alpha(n_bins)
    if (temp_check_size > readin.size()){
      int_end = 2;
      logger << LOG_DEBUG << "int_end = 2   temp_check_size = " << temp_check_size << " > " << readin.size() << " = readin.size()" << endl;
      return;
    }
    else {
      logger << LOG_DEBUG << "temp_check_size = " << temp_check_size << " <= " << readin.size() << " = readin.size()" << endl;
    }
  }
  
  if (end_optimization == 0){
    //  count_channel, n_rej_channel, sum_channel_weight (3 * n_bins)
    temp_check_size += 3 * n_bins;
    if (temp_check_size > readin.size()){
      int_end = 2;
      logger << LOG_DEBUG << "int_end = 2   temp_check_size = " << temp_check_size << " > " << readin.size() << " = readin.size()" << endl;
      return;
    }
    else{
      logger << LOG_DEBUG << "temp_check_size = " << temp_check_size << " <= " << readin.size() << " = readin.size()" << endl;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randomvariable::proceeding_save(){
  static Logger logger("randomvariable::proceeding_save");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  save_end_optimization = end_optimization;
  save_alpha = alpha;
  save_beta = beta;
  save_n_acc_channel = n_acc_channel;
  save_n_rej_channel = n_rej_channel;
  save_sum_channel_weight = sum_channel_weight;
  save_s = s;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randomvariable::save_weights(ofstream & out_weights){
  static Logger logger("randomvariable::save_weights");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  out_weights << "variable " << setw(50) << name << endl;

  for (int i_b = 0; i_b < n_bins; i_b++){
    // Why such low precision ???
    out_weights << setprecision(8) << alpha[i_b] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

// Change to void ... !!! Errors will anyway stop run.
void randomvariable::readin_weights(vector<string> &readin) {
  //int randomvariable::readin_weights(vector<string> &readin) {
  static Logger logger("randomvariable::readin_weights");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // nothing to be read in for switch_IS_MC != 2, 3:
  if (switch_IS_MC != 2 && switch_IS_MC != 3) {return;}
  //  if (switch_IS_MC != 2 && switch_IS_MC != 3) {return 0;}

  // check if correct !!!
  // find the correct section in the weights file
  int proc;
  for (proc=0; proc<readin.size(); proc++) {
    size_t pos=readin[proc].find("variable ");
    if (pos != string::npos) {
      pos += 9;
      string var_name=readin[proc].substr(pos,string::npos);
      int start_index = 0;
      for (int i = 0; i < var_name.size(); i++){
        if (var_name[i] != ' '){start_index = i; break;}
      }
      var_name=var_name.substr(start_index, var_name.size() - start_index);
      if (var_name == name) {
        logger << LOG_DEBUG_VERBOSE << "weights found for " << var_name << endl;
        break;
      }
    }
  }

  if (proc == readin.size()) {
    logger << LOG_ERROR << "could not read in weights for " << name << endl;
    assert(false);
    //    return 1;
  }

  logger << LOG_DEBUG_VERBOSE << "reading weights for " << readin[proc] << endl;
  proc++;

  for (int i_b = 0; i_b < n_bins; i_b++){
    alpha[i_b] = atof(readin[proc].c_str());
    proc++;
    logger << LOG_DEBUG_VERBOSE << "alpha[" << setw(4) << i_b << "] = " << setw(23) << setprecision(15) << alpha[i_b] << endl;
  }
  //  if (switch_IS_MC == 0 || switch_IS_MC == 2){logger << LOG_DEBUG << readin[proc] << " is not further optimized." << endl;}
  //  else if (switch_IS_MC == 1 || switch_IS_MC == 3){logger << LOG_DEBUG << readin[proc] << " is further optimized." << endl;}

  beta[0] = alpha[0];
  for (int i_b = 1; i_b < n_bins; i_b++){beta[i_b] = beta[i_b - 1] + alpha[i_b];}

  // Check for issues in read-in weights:
  if (fabs(beta[n_bins - 1] - 1) > 1e-6){
    for (int i_b = 0; i_b < n_bins; i_b++){
      logger << LOG_ERROR << "alpha[" << setw(4) << i_b << "] = " << alpha[i_b] << "   beta[" << setw(4) << i_b << "] = " << beta[i_b] << endl;
    }
    logger << LOG_FATAL << setw(50) << name << ": something went wrong in weight read-in, aborting.." << endl;
    exit(1);

    //    logger << LOG_DEBUG_VERBOSE << "finished   return value: 1" << endl;
    //    return 1;
  }
  /*
  logger << LOG_DEBUG_VERBOSE << "finished   return value: 0" << endl;
  return 0;
  */
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




/*
// redundant: should be replaced by step_IS_optimization
int randomvariable::do_optimization_step(int events) {
  static Logger logger("do_optimization_step");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  if (!imp_sampling)
    return 2;

  if (end_optimization == 0){
    //  if (opt_end == 0){
    if (events % n_events == 0 && events != 0){
      //logger << LOG_DEBUG << name << " vegasgrid_calculation @ " << events << " events" << "(" << accumulate(n_acc_channel.begin(),n_acc_channel.end(),0)+accumulate(n_rej_channel.begin(),n_rej_channel.end(),0) << " events used)" << endl;
      vegasgrid_calculation_temp(alpha, beta, n_rej_channel, n_acc_channel, sum_channel_weight);
      logger << LOG_DEBUG << left << setw(50) << name << ": vegasgrid_calculation @ " << events << " events" << "(" << accumulate(n_acc_channel.begin(),n_acc_channel.end(),0)+accumulate(n_rej_channel.begin(),n_rej_channel.end(),0) << " events used)" << "; variation is " << calc_variation() << endl;
      /*
      for (int is = 0; is < alpha.size(); is++){
	logger << LOG_DEBUG << "alpha[" << right << setw(4) << is << "] = " << left << setw(15) << setprecision(8) << alpha[is] << "   beta[" << right << setw(4) << is << "] = " << left << setw(15) << setprecision(8) << beta[is] << endl;
      }
*//*

      if (events == n_events * n_steps){
        logger << LOG_DEBUG << "sampled random variable " << setw(50) << left << name << ":" << endl;
        for (int is = 0; is < alpha.size(); is++){
	  logger << LOG_DEBUG << "alpha[" << right << setw(4) << is << "] = " << left << setw(15) << setprecision(8) << alpha[is] << "   beta[" << right << setw(4) << is << "] = " << left << setw(15) << setprecision(8) << beta[is] << endl;
        }
	//        opt_end = 1;
	end_optimization = 1;
	
        logger << LOG_DEBUG << setw(50) << left << name << ": variation is " << calc_variation() << endl;
        return 1;
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return 0;
}

*/


/*
void vegasgrid_average_temp(vector<double> & alpha, double averaging=1.0){
  double tmp=alpha[0],tmp2;
//  alpha[0] = 0.5*(averaging+2)/(1+averaging) * alpha[0] + 0.5/(1+averaging) * alpha[1];
//  for (int j = 1; j < alpha.size() - 1; j++){
////    alpha_new[j] = .5 * alpha[j] + .25 * alpha[j - 1] + .25 * alpha[j + 1];
//      tmp2=alpha[j];
//      alpha[j] = 1.0/(1+averaging) * alpha[j] + 0.5/(1+averaging) * tmp + 0.5/(1+averaging) * alpha[j + 1];
//      tmp=tmp2;
//  }
//  alpha[alpha.size() - 1] = 0.5*(averaging+2)/(1+averaging) * alpha[alpha.size() - 1] + 0.5/(1+averaging) * tmp;

  alpha[0] = .75 * alpha[0] + .25 * alpha[1];
  for (int j = 1; j < alpha.size() - 1; j++){
//    alpha_new[j] = .5 * alpha[j] + .25 * alpha[j - 1] + .25 * alpha[j + 1];
      tmp2=alpha[j];
      alpha[j] = .5 * alpha[j] + .25 * tmp + .25 * alpha[j + 1];
      tmp=tmp2;
  }
  alpha[alpha.size() - 1] = .75 * alpha[alpha.size() - 1] + .25 * tmp;
}

void vegasgrid_calculation_temp(vector<double> & alpha, vector<double> & beta, vector<int> & n_rej_channel, vector<int> & count_channel, vector<double> & sum_channel_weight){
  static Logger logger("vegasgrid_calculation_temp");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int j = 0; j < alpha.size(); j++){
    if (count_channel[j] + n_rej_channel[j] != 0){
      alpha[j] = sum_channel_weight[j] / (count_channel[j] + n_rej_channel[j]);
    }
    else {
      alpha[j] = 0.;
    }
  }

  double sum_value = std::accumulate(alpha.begin(), alpha.end(), 0.);

  if (sum_value == 0.){logger << LOG_DEBUG << "1st: sum_value = 0." << endl;}
  if (sum_value == 0.){
    for (int j = 0; j < alpha.size(); j++){
      alpha[j] = 1. / alpha.size();
      sum_value = 1.;
    }
  }

  for (int j = 0; j < alpha.size(); j++){alpha[j] = alpha[j] / sum_value;}

  int MCweight_lc_min = 0;
  vector<int> MCweights_set_to_min(alpha.size());

  double MCweight_limit_min = 1.e-2 / alpha.size();
  sum_value=0.0;
  for (int j = 0; j < alpha.size(); j++){
    if (alpha[j] < MCweight_limit_min){
      MCweight_lc_min++;
      MCweights_set_to_min[j] = 1;
    }
    else {
      sum_value+=alpha[j];
    }
  }
  if (sum_value == 0.){logger << LOG_DEBUG << "2nd: sum_value = 0." << endl;}

  double normalization = (1. - MCweight_lc_min * MCweight_limit_min) / sum_value;
  if (normalization != normalization){logger << LOG_WARN << "normalization != normalization" << endl;}
  if (normalization > 1.e200){logger << LOG_WARN << "normalization = +-inf" << endl;}

  for (int j = 0; j < alpha.size(); j++){
    if (MCweights_set_to_min[j] == 1){alpha[j] = MCweight_limit_min;}
    else {alpha[j] = alpha[j] * normalization;}
  }

  vegasgrid_average_temp(alpha);

  beta[0] = alpha[0];
  for (int i = 1; i < beta.size(); i++){beta[i] = beta[i - 1] + alpha[i];}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

*/

/*
void randomvariable::optimize(const double psp_weight) {
  if (!imp_sampling || !used)
    return;

  // TODO: stop recording if opt_end==1?

  n_acc_channel[channel]++;
  sum_channel_weight[channel] += fabs(psp_weight);
  //sum_channel_weight2[channel] += psp_weight*psp_weight;
  //cout << "RV: channel = " << channel << endl;
  // cout << "RV: n_acc_channel[" << setw(3) << channel << "] = " << n_acc_channel[channel] << endl;

  used=false;
}
*/


/*
void randomvariable::multiply_g_IS(double r, double & g_IS) {
  if (r < 0.0 || r > 1.0 || r != r){
    cout << left << setw(40) << name << "r = " << setw(15) << setprecision(8) << r << " alpha.size() = " << setw(2) << alpha.size() << "   a[" << setw(2) << channel << "] = " << setw(15) << setprecision(8) << alpha[channel] << "   g_IS = " << setw(15) << setprecision(8) << g_IS << endl;
  }
  
//  assert(r >= 0.);
//  assert(r <= 1.);
  
  int temp_channel = (int)(r * (double)beta.size());
  if (r < 0.0 || r > 1.0 || r != r){temp_channel = 0;}

  g_IS *= alpha[temp_channel] * alpha.size();
  if (r < 0.0 || r > 1.0 || r != r){g_IS = sqrt(-1.);}
}
*/


/*
double randomvariable::calc_variation() {
  double average=0;
  average=accumulate(alpha.begin(),alpha.end(),0.0);
  average/=alpha.size();

  double variation=0.0;
  for (int i=0; i<alpha.size(); i++) {
    variation+=pow(alpha[i]-average,2);
  }

  if (variation == 0.){variation = 1.;} // in this case, this channel has not been called at all! 0 might be problematic...

  // variation should be invariant under rescaling
  variation*=alpha.size();

  if (sqrt(variation) != sqrt(variation)){
    cout << "alpha.size() = " << alpha.size() << endl;
    for (int i_a = 0; i_a < alpha.size(); i_a++){
      cout << "alpha[" << i_a << "] = " << alpha[i_a] << endl;
    }
    cout << "variation = " << variation << endl;
    cout << "average 1st = " << accumulate(alpha.begin(),alpha.end(),0.0) << endl;
    cout << "average 2nd = " << accumulate(alpha.begin(),alpha.end(),0.0) / alpha.size() << endl;
  }


  return sqrt(variation);
}
*/



/*
void randomvariable::step_IS_optimization_temp(){
  static Logger logger("step_IS_optimization_temp");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;


  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/



/*
// Could be outdated !!!
void randomvariable::proceeding_size(int & proc) {
  static Logger logger("randomvariable::proceeding_size(int&)");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  proc += 1;
  if (end_optimization != -1){
    proc += 3;
    proc += n_bins;
  }
  if (end_optimization == 0){
    proc += 3 * n_bins;
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

// Could be outdated !!!
void randomvariable::proceeding_size(vector<int> & size_proceeding) {
  static Logger logger("randomvariable::proceeding_size(vector<int>&)");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  size_proceeding[0] += 1;
  if (end_optimization != -1){
    size_proceeding[0] += 3;
    size_proceeding[0] += n_bins;
  }
  if (end_optimization == 0){
    size_proceeding[5] += 3 * n_bins;
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/



