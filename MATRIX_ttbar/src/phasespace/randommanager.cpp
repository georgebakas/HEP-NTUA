#include "../include/classes.cxx"

randommanager::randommanager(){}
//randommanager::randommanager() : logger(Logger("random manager")) {}

randommanager::randommanager(vector<double> sran, string _weightdir, string _processname){
  //randommanager::randommanager(vector<double> sran, string _weightdir, string _processname) : logger(Logger("random manager")) {
  static Logger logger("randommanager::randommanager");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;


  s=sran;

  for (int i = 0; i < s.size(); i++){logger << LOG_DEBUG_VERBOSE << "s[" << i << "] = " << s[i] << endl;}

  // hack to get uncorrelated PRN
  s[0]/=2;
  s[1]/=3;
  s[2]/=4;

  for (int i = 0; i < s.size(); i++){logger << LOG_DEBUG_VERBOSE << "s[" << i << "] = " << s[i] << endl;}

  weightdir=_weightdir;
  processname=_processname;

  //  opt_end=true;
  // check if reasonable initialization
  end_optimization = 1;

  used_queue.resize(0);

  // Check meaning of queue_counter !!!
  queue_counter=0;

  // Check meaning of queue_threshold !!!
  queue_threshold=1e2;

  // Check order ??? proceeding before readin ???
  proceeding_save();

  logger << LOG_DEBUG << "weightdir = " << weightdir << endl;

  readin_weights();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

randommanager::~randommanager(){
  static Logger logger("randommanager::~randommanager");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_r = 0; i_r < random_variables.size(); i_r++) {
    //  FIXME: should the random manager be responsible for deleting the variables it manages?
    //  If so, we have to be careful that they actually still exists at this point
    //  delete random_variables[i_r];
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

///void randommanager::register_variable(randomvariable* new_var, bool imp_sampling) {
void randommanager::register_variable(randomvariable* new_var){
  static Logger logger("randommanager::register_variable");
  logger << LOG_DEBUG << "called" << endl;
  
  ostringstream convert;

  for (int i = 0; i < s.size(); i++){logger << LOG_DEBUG_VERBOSE << "s[" << i << "] = " << s[i] << endl;}

  if (new_var->name=="") {
    convert << "random_var" << random_variables.size();
    new_var->name = convert.str();
  }

  ///  new_var->imp_sampling=imp_sampling;


  //  SK TODO: completely remove the individual s(3) vectors to generate random numbers
  //  Should be done completely in a single one, within randommanager !!!
  // decorrelate random sequences
  // TODO: seed uniformly!
  new_var->s[0] = (new_var->s[0]+s[0])/4*(1+cos(random_variables.size()+2));
  new_var->s[1] = (new_var->s[1]+s[1])/4*(1+cos(random_variables.size()+2));
  new_var->s[2] = (new_var->s[2]+s[2])/4*(1+cos(random_variables.size()+2));

  for (int i = 0; i < new_var->s.size(); i++){logger << LOG_DEBUG_VERBOSE << "new_var->s[" << i << "] = " << new_var->s[i] << endl;}

  new_var->proceeding_save();

  new_var->manager=this;

  new_var->readin_weights(readin);

  random_variables.push_back(new_var);

  logger << LOG_DEBUG << "new random variable registered " << setw(50) << new_var->name << ": " << new_var->n_events << ", " << new_var->switch_IS_MC << ", " << new_var->n_steps << ", " << new_var->n_bins << endl;
  ///", " << new_var->imp_sampling << endl;

  //  opt_end=false;
  end_optimization = 0;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randommanager::psp_IS_optimization(double & psp_weight, double & psp_weight2) {
  static Logger logger("randommanager::psp_IS_optimization");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (end_optimization){return;}

  end_optimization = 1;
  // very confusing distinction !!! queue_threshold ???
  if (random_variables.size()>queue_threshold) {
    //    logger << LOG_DEBUG << "random_variables.size() > queue_threshold" << endl;
    for (int i=0; i<queue_counter; i++) {
      used_queue[i]->psp_IS_optimization(psp_weight);
     if (used_queue[i]->end_optimization == 0) {
        end_optimization = 0;
      }
    }
  }
  else {
    //    logger << LOG_DEBUG << "random_variables.size() <= queue_threshold" << endl;
    int first = 0;
    int last = random_variables.size() - 1;
    for (int i = first; i < last + 1; i++) {
      random_variables[i]->psp_IS_optimization(psp_weight);
      if (random_variables[i]->end_optimization == 0) {
        end_optimization = 0;
      }
    }
  }

  // ???
  queue_counter = 0;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void randommanager::increase_counter_n_acc(){
  static Logger logger("randommanager::increase_counter_n_acc");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (end_optimization){return;}

  if (random_variables.size()>queue_threshold) {
    for (int i=0; i<queue_counter; i++) {
      used_queue[i]->n_acc_channel[used_queue[i]->channel]++;
      if (used_queue[i]->end_optimization == 0) {
        used_queue[i]->used=false;
      }
    }
  }
  else {
    int first=0;
    int last=random_variables.size()-1;

    for (int i=first; i<last+1; i++) {
      if (random_variables[i]->used) {
        random_variables[i]->n_acc_channel[random_variables[i]->channel]++;
        random_variables[i]->used=false;
      }
    }
  }

  queue_counter=0;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randommanager::increase_cut_counter() {
  static Logger logger("randommanager::increase_cut_counter");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (end_optimization){return;}
  //  if (opt_end)
  //    return;

  if (random_variables.size()>queue_threshold) {
    for (int i=0; i<queue_counter; i++) {
      used_queue[i]->n_rej_channel[used_queue[i]->channel]++;
      //      if (used_queue[i]->opt_end == 0) {
      if (used_queue[i]->end_optimization == 0) {
        used_queue[i]->used=false;
      }
    }
  }
  else {
    int first=0;
    int last=random_variables.size()-1;

    for (int i=first; i<last+1; i++) {
      if (random_variables[i]->used) {
        random_variables[i]->n_rej_channel[random_variables[i]->channel]++;
        random_variables[i]->used=false;
      }
    }
  }

  queue_counter=0;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randommanager::do_optimization_step(int events) {
  static Logger logger("randommanager::do_optimization_step");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (end_optimization){return;}
  //  if (opt_end)
  //    return;

  //  int finished;
  //  bool writeout=false;
  //  cout << "events = " << events << endl;
  for (int i = 0; i < random_variables.size(); i++){
    random_variables[i]->step_IS_optimization(events);
    /*
    //    cout << "random_variables[i]->do_optimization_step(events) = " << random_variables[i]->do_optimization_step(events) << endl;
    finished=random_variables[i]->do_optimization_step(events);
    if (finished==1){
      writeout=true;
      //      end_optimization = 1;
    }
    */
  }
  end_optimization = 1;
  for (int i = 0; i < random_variables.size(); i++){
    if (random_variables[i]->end_optimization == 0){end_optimization = 0; break;}
  }

  if (end_optimization == 1){writeout_weights();}
  /*
  if (writeout) {
    writeout_weights();
//    //ofstream out_weights(("weights/weights_IS_MC_" + processname + ".dat").c_str());
//    ofstream out_weights(("weights/weights_IS_MC_" + processname + ".dat").c_str());
//
//    for (int i=0; i<random_variables.size(); i++) {
//      random_variables[i]->save_weights(out_weights);
//    }
  }
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randommanager::proceeding_out(ofstream &out_proceeding) {
  static Logger logger("randommanager::proceeding_out");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i = 0; i < 3; i++){out_proceeding << double2hexastr(save_s[i]) << endl;}

  for (int i_r = 0; i_r < random_variables.size(); i_r++){
    random_variables[i_r]->proceeding_out(out_proceeding);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randommanager::proceeding_in(int &proc, vector<string> &readin) {
  static Logger logger("randommanager::proceeding_in");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i = 0; i < 3; i++){s[i] = hexastr2double(readin[proc + i]);}
  proc += 3;

  for (int i_r = 0; i_r < random_variables.size(); i_r++){
    random_variables[i_r]->proceeding_in(proc, readin);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randommanager::check_proceeding_in(int &int_end, int &temp_check_size, vector<string> &readin) {
  static Logger logger("randommanager::check_proceeding_in");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG << "int_end = " << int_end << endl;
  logger << LOG_DEBUG << "temp_check_size = " << temp_check_size << endl;
  logger << LOG_DEBUG << "readin.size() = " << readin.size() << endl;

  //  s(3)
  temp_check_size += 3;
  if (temp_check_size > readin.size()){
    int_end = 2;
    logger << LOG_DEBUG << "int_end = 2   temp_check_size = " << temp_check_size << " > " << readin.size() << " = readin.size()" << endl;
    return;
  }
  else {
    logger << LOG_DEBUG << "temp_check_size = " << temp_check_size << " <= " << readin.size() << " = readin.size()" << endl;
  }

  logger << LOG_DEBUG << "random_variables.size() = " << random_variables.size() << endl;
  
  for (int i_r = 0; i_r < random_variables.size(); i_r++){
    logger << LOG_DEBUG << "random_variables[" << i_r << "] = " << random_variables[i_r]->name << "   n_bins = " << random_variables[i_r]->n_bins << endl;
    random_variables[i_r]->check_proceeding_in(int_end, temp_check_size, readin);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



/*
void randommanager::proceeding_size(int &proc){
  static Logger logger("randommanager::proceeding_size(int&)");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  proc += 3;
  for (int i_r = 0; i_r < random_variables.size(); i_r++){
    random_variables[i_r]->proceeding_size(proc);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void randommanager::proceeding_size(vector<int> & size_proceeding){
  static Logger logger("randommanager::proceeding_size(vector<int>&)");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  size_proceeding[0] += 3;
  for (int i_r = 0; i_r < random_variables.size(); i_r++){
    random_variables[i_r]->proceeding_size(size_proceeding);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/



void randommanager::proceeding_save() {
  static Logger logger("randommanager::proceeding_save");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  save_s = s;

  for (int i_r = 0; i_r < random_variables.size(); i_r++){
    random_variables[i_r]->proceeding_save();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void randommanager::readin_weights(){
  static Logger logger("randommanager::readin_weights");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

//  // check if we need to readin weights
//  bool bool_readin=false;
//  for (int i=0; i<random_variables.size(); i++) {
//    if (random_variables[i]->switch_IS_MC == 2 || random_variables[i]->switch_IS_MC == 3)
//      bool_readin=true;
//  }
//
//  if (!bool_readin) {
//    logger << LOG_DEBUG << "no weights need to be read" << endl;
//    return;
//  }

  // check if weights should actually be read in !!!
  
  string in_path = weightdir + "/weights_IS_MC_" + processname + ".dat";
  logger << LOG_DEBUG << "in_path = " << in_path << endl;
  ifstream readin_weights(in_path.c_str());

  if (readin_weights){
    char LineBuffer[128];
    while (readin_weights.getline(LineBuffer, 128)){readin.push_back(LineBuffer);}
    readin_weights.close();
  }
  else {
    //  if (!readin_weights) {
    logger << LOG_INFO << "weight file " << in_path << " could not be opened" << endl;
    readin_weights.close();
    return;
  }

  /*
  char LineBuffer[128];
  while (readin_weights.getline(LineBuffer, 128)) {readin.push_back(LineBuffer);}
  readin_weights.close();
  */
//  vector<bool> successfully_readin(random_variables.size(),false);
//
//  for (int proc=0; proc<readin.size(); proc++) {
//    size_t pos=readin[proc].find("variable ");
//    if (pos != string::npos) {
//      pos += 9;
//      string var_name=readin[proc].substr(pos,string::npos);
//      //      logger << LOG_DEBUG << "weight readin var_name = " << var_name << endl;
//      int start_index = 0;
//      for (int i = 0; i < var_name.size(); i++){
//        if (var_name[i] != ' '){start_index = i; break;}
//      }
//      var_name=var_name.substr(start_index, var_name.size() - start_index);
//      logger << LOG_INFO << "weights found for " << var_name << endl;
//
//      for (int i=0; i<random_variables.size(); i++) {
//	//	logger << LOG_DEBUG << "weight readin var_name = " << random_variables[i]->name << endl;
//      if (var_name==random_variables[i]->name) {
//        if (random_variables[i]->readin_weights(proc,readin)) {
//          logger << LOG_FATAL << "weight readin failed for variable " << random_variables[i]->name << ", check if ok!" << endl;
//	        exit(1);
//        } else {
//            successfully_readin[i]=true;
//          }
//        }
//      }
//    }
//  }
//
//  // verify that every variable with has successfully readin its weights
//  for (int i=0; i<random_variables.size(); i++) {
//    logger << LOG_DEBUG << "checking " << random_variables[i]->name << endl;
//    if ((random_variables[i]->switch_IS_MC == 2 || random_variables[i]->switch_IS_MC == 3) && !successfully_readin[i]) {
//      logger << LOG_FATAL << "weight readin failed for variable " << random_variables[i]->name << ", check if ok!" << endl;
//      //exit(1);
//    }
//  }
//


  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void randommanager::writeout_weights() {
  static Logger logger("randommanager::writeout_weights");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_weights(("weights/weights_IS_MC_" + processname + ".dat").c_str());
  for (int i_r = 0; i_r < random_variables.size(); i_r++){
    random_variables[i_r]->save_weights(out_weights);
  }
  end_optimization = 2;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void randommanager::add_var_to_queue(randomvariable* variable) {
  static Logger logger("randommanager::add_var_to_queue");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (end_optimization){return;}

  // check if this distinction makes sense !!!
  if (queue_counter >= used_queue.size()){
    used_queue.push_back(variable);
  }
  else {
    used_queue[queue_counter] = variable;
  }
  queue_counter++;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void randommanager::nancut() {
  static Logger logger("randommanager::nancut");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_r = 0; i_r < random_variables.size(); i_r++){
    if (random_variables[i_r]->used){
      random_variables[i_r]->used = false;
    }
  }
  queue_counter = 0;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}





/*
void randommanager::nancut(bool count_as_cut) {
  for (int i=0; i<random_variables.size(); i++) {
    if (random_variables[i]->used) {
      if (count_as_cut) {
        random_variables[i]->n_rej_channel[random_variables[i]->channel]++;
      }
      random_variables[i]->used=false;
    }
  }
  queue_counter=0;
}
*/


/*
void randommanager::optimize(const double psp_weight, const double psp_weight2) {
  static Logger logger("randommanager::optimize");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (end_optimization){return;}
  /*
  if (opt_end)
    return;
  
  opt_end=true;
*//*
  end_optimization = 1;

  // very confusing distinction !!!
  if (random_variables.size()>queue_threshold) {
    //    logger << LOG_DEBUG << "random_variables.size() > queue_threshold" << endl;
    for (int i=0; i<queue_counter; i++) {
      used_queue[i]->psp_IS_optimization(psp_weight);
      //      used_queue[i]->optimize(psp_weight);
      //      if (used_queue[i]->opt_end == 0) {
      //        opt_end=false;
      if (used_queue[i]->end_optimization == 0) {
        end_optimization = 0;
      }
    }
  }
  else {
    //    logger << LOG_DEBUG << "random_variables.size() <= queue_threshold" << endl;
    int first=0;
    int last=random_variables.size()-1;

    for (int i=first; i<last+1; i++) {
      random_variables[i]->psp_IS_optimization(psp_weight);
      //      random_variables[i]->optimize(psp_weight);
      //      if (random_variables[i]->opt_end == 0) {
      //        opt_end=false;
      if (random_variables[i]->end_optimization == 0) {
        end_optimization = 0;
      }
    }
  }

  queue_counter=0;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
  */
