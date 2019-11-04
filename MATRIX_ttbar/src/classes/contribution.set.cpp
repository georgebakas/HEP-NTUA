#include "../include/classes.cxx"

////////////////////
//  constructors  //
////////////////////
contribution_set::contribution_set(){
  switch_to_newscheme = 0;
  
  fill_code_particle();
  fill_name_particle();

  basic_process_class = "";
  process_class = "";
  decay.resize(0);
  type_perturbative_order = "";
  type_contribution = "";
  type_correction = "";
  contribution_order_alpha_s = 0;
  contribution_order_alpha_e = 0;
  contribution_order_interference = 0;

  basic_subprocess = "";
  subprocess = "";
  type_parton.resize(1);
  n_particle = 0;
  process_type = 0;

  order_alpha_s_born = 0;
  n_jet_born = 0;
  n_particle_born = 0;
  
  for (int j_l = 1; j_l < 7; j_l++){
    lepton_exchange[10 + j_l] = 10 + j_l;
    lepton_exchange[-10 - j_l] = -10 - j_l;
  }
  //  run_mode = "";

}


// This constructor does not seem to be used at all !!!
contribution_set::contribution_set(string _process_class, vector<string> _decay, string _type_perturbative_order, string _type_contribution, string _type_correction, int _order_alpha_s, int _order_alpha_e, int _order_interference, string _subprocess){
  Logger logger("contribution_set::contribution_set");
  logger << LOG_DEBUG << "started" << endl;
  //  logger << LOG_DEBUG << "The same information is also contained in observable_set !!!" << endl;

  //  basic_process_class = _basic_process_class;
  process_class = _process_class;
  decay = _decay;

  if (process_class[process_class.size() - 1] == 'X'){switch_to_newscheme = 0;}
  else {switch_to_newscheme = 1;}
  logger << LOG_DEBUG_VERBOSE << "switch_to_newscheme = " << switch_to_newscheme << endl;

  fill_code_particle();
  fill_name_particle();

  type_perturbative_order = _type_perturbative_order;
  type_contribution = _type_contribution;
  type_correction = _type_correction;

  contribution_order_alpha_s = _order_alpha_s;
  contribution_order_alpha_e = _order_alpha_e;
  contribution_order_interference = _order_interference;

  // ???
  //  basic_subprocess = _basic_subprocess;
  subprocess = _subprocess;

  //  type_parton.resize(1);
  //  subprocess_readin(process_class, subprocess, decay, process_type, n_particle, type_parton);
  readin_hadronic_process();
  readin_subprocess();

  for (int j_l = 1; j_l < 7; j_l++){
    lepton_exchange[10 + j_l] = 10 + j_l;
    lepton_exchange[-10 - j_l] = -10 - j_l;
  }

  logger << LOG_DEBUG << "basic_process_class = " << basic_process_class << endl;
  logger << LOG_DEBUG << "process_class = " << process_class << endl;
  for (int i_v = 0; i_v < decay.size(); i_v++){
    logger << LOG_DEBUG << "decay[" << i_v << "] = " << decay[i_v] << endl;
  }

  logger << LOG_DEBUG << "type_perturbative_order = " << type_perturbative_order << endl;
  logger << LOG_DEBUG << "type_contribution = " << type_contribution << endl;
  logger << LOG_DEBUG << "type_correction = " << type_correction << endl;
  
  logger << LOG_DEBUG << "contribution_order_alpha_s = " << contribution_order_alpha_s << endl;
  logger << LOG_DEBUG << "contribution_order_alpha_e = " << contribution_order_alpha_e << endl;
  logger << LOG_DEBUG << "contribution_order_interference = " << contribution_order_interference << endl;
  
  logger << LOG_DEBUG << "subprocess = " << subprocess << endl;


  determination_order_alpha_s_born();
  determination_class_contribution();

 
  logger << LOG_DEBUG << "finished" << endl;
}


void contribution_set::determination_order_alpha_s_born(){
  Logger logger("contribution_set::order_alpha_s_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (type_contribution == ""){logger << LOG_FATAL << "No contribution chosen!" << endl; exit(1);}
  
  if (type_contribution == "born"){
    order_alpha_s_born = contribution_order_alpha_s;
  } 
  else if (type_correction == "QCD" && 
	   (type_contribution == "RA" ||
	    type_contribution == "VA" ||
	    type_contribution == "CA" ||
	    type_contribution == "RT" ||
	    type_contribution == "VT" ||
	    type_contribution == "CT" ||
	    type_contribution == "RJ" ||
	    type_contribution == "VJ" ||
	    type_contribution == "CJ")){
    order_alpha_s_born = contribution_order_alpha_s - 1;
  }
  else if (type_correction == "QCD" && 
	   (type_contribution == "VT2" ||
	    type_contribution == "CT2" ||
	    type_contribution == "RVA" ||
	    type_contribution == "RCA" ||
	    type_contribution == "RRA" ||
	    type_contribution == "VJ2" ||
	    type_contribution == "CJ2" ||
	    type_contribution == "RVJ" ||
	    type_contribution == "RCJ" ||
	    type_contribution == "RRJ")){
    order_alpha_s_born = contribution_order_alpha_s - 2;
  }
  else if (//type_correction == "QCD" && 
	   (type_contribution == "loop" ||
	    type_contribution == "LI2")){
    order_alpha_s_born = contribution_order_alpha_s;
    // The loop-induced process is the "Born" process here !!!
    //    order_alpha_s_born = contribution_order_alpha_s - 2;
  }
  else if (//type_correction == "QCD" && 
	   (type_contribution == "L2RA" ||
	    type_contribution == "L2VA" ||
	    type_contribution == "L2CA" ||
	    type_contribution == "L2RT" ||
	    type_contribution == "L2VT" ||
	    type_contribution == "L2CT" ||
	    type_contribution == "L2RJ" ||
	    type_contribution == "L2VJ" ||
	    type_contribution == "L2CJ")){
    order_alpha_s_born = contribution_order_alpha_s - 1;
    // The loop-induced process is the "Born" process here !!!
    //    order_alpha_s_born = contribution_order_alpha_s - 3;
  }
  else {
  logger << LOG_INFO << "Should not happen !!!" << endl;
    order_alpha_s_born = contribution_order_alpha_s;
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}

void contribution_set::determination_class_contribution(){
  Logger logger("contribution_set::determination_class_contribution");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  class_contribution_born = 0;
  class_contribution_CS_collinear = 0;
  class_contribution_CS_virtual = 0;
  class_contribution_CS_real = 0;
  class_contribution_QT_collinear = 0;
  class_contribution_QT_virtual = 0;
  class_contribution_QT_real = 0;
  class_contribution_NJ_collinear = 0;
  class_contribution_NJ_virtual = 0;
  class_contribution_NJ_real = 0;
  class_contribution_collinear = 0;
  class_contribution_virtual = 0;
  class_contribution_real = 0;
  class_contribution_loopinduced = 0;
  class_contribution_IRcut_implicit = 0;
  class_contribution_IRcut_explicit = 0;
  class_contribution_CS = 0;
  class_contribution_QT = 0;
  class_contribution_NJ = 0;
  class_contribution_QT_CS = 0;
  class_contribution_NJ_CS = 0;

  if (type_contribution == ""){logger << LOG_FATAL << "No contribution chosen!" << endl; exit(1);}

  // class_contribution_... simplifies type_contribution-dependent statements for similar contributions
  
  if (type_contribution == "born"){
    class_contribution_born = 1;
  }
 
  if (type_contribution == "CT" ||
      type_contribution == "CT2" ||
      type_contribution == "CJ" ||
      type_contribution == "CJ2" ||
      type_contribution == "L2CT" ||
      type_contribution == "L2CJ"){
    class_contribution_IRcut_implicit = 1;
  }

  if (type_contribution == "RT" ||
      type_contribution == "RJ" ||
      type_contribution == "RVA" ||
      type_contribution == "RCA" ||
      type_contribution == "RRA" ||
      type_contribution == "RVJ" ||
      type_contribution == "RCJ" ||
      type_contribution == "RRJ"){
    class_contribution_IRcut_explicit = 1;
  }
  
  if (type_contribution == "VA" ||
      type_contribution == "RVA" ||
      type_contribution == "RVJ" ||
      type_contribution == "L2VA"){
    class_contribution_CS_virtual = 1;
  }
  if (type_contribution == "VT" ||
      type_contribution == "L2VT" ||
      type_contribution == "VT2"){
    class_contribution_QT_virtual = 1;
  }
  if (type_contribution == "VJ" ||
      type_contribution == "L2VJ" ||
      type_contribution == "VJ2"){
    class_contribution_NJ_virtual = 1;
  }

  if (type_contribution == "CA" ||
      type_contribution == "RCA" ||
      type_contribution == "RCJ" ||
      type_contribution == "L2CA"){
    class_contribution_CS_collinear = 1;
  }
  if (type_contribution == "CT" ||
      type_contribution == "L2CT" ||
      type_contribution == "CT2"){
    class_contribution_QT_collinear = 1;
  }
  if (type_contribution == "CJ" ||
      type_contribution == "L2CJ" ||
      type_contribution == "CJ2"){
    class_contribution_NJ_collinear = 1;
  }

  if (type_contribution == "RA" ||
      type_contribution == "RRA" ||
      type_contribution == "RRJ" ||
      type_contribution == "L2RA"){
    class_contribution_CS_real = 1;
  }
  
  if (type_contribution == "loop" ||
      type_contribution == "L2I" ||
      type_contribution == "L2VA" ||
      type_contribution == "L2CA" ||
      type_contribution == "L2RA" ||
      type_contribution == "L2VT" ||
      type_contribution == "L2CT" ||
      type_contribution == "L2RT" ||
      type_contribution == "L2VJ" ||
      type_contribution == "L2CJ" ||
      type_contribution == "L2RJ"){
    class_contribution_loopinduced = 1;
  }
  
  if (type_contribution == "VA" ||
      type_contribution == "CA" ||
      type_contribution == "RA" ||
      type_contribution == "L2VA" ||
      type_contribution == "L2CA" ||
      type_contribution == "L2RA"){
    class_contribution_CS = 1;
  }
  
  
  if (type_contribution == "VT" ||
      type_contribution == "CT" ||
      type_contribution == "RT" ||
      type_contribution == "L2VT" ||
      type_contribution == "L2CT" ||
      type_contribution == "L2RT"){
    class_contribution_QT = 1;
  }
  
  if (type_contribution == "VJ" ||
      type_contribution == "CJ" ||
      type_contribution == "RJ" ||
      type_contribution == "L2VJ" ||
      type_contribution == "L2CJ" ||
      type_contribution == "L2RJ"){
    class_contribution_NJ = 1;
  }
  
  if (type_contribution == "VT2" ||
      type_contribution == "CT2" ||
      type_contribution == "RVA" ||
      type_contribution == "RCA" ||
      type_contribution == "RRA"){
    class_contribution_QT_CS = 1;
  }
  
  if (type_contribution == "VJ2" ||
      type_contribution == "CJ2" ||
      type_contribution == "RVJ" ||
      type_contribution == "RCJ" ||
      type_contribution == "RRJ"){
    class_contribution_NJ_CS = 1;
  }

  logger << LOG_INFO << setw(35) << "class_contribution_born" << " = " << class_contribution_born << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_IRcut_implicit" << " = " << class_contribution_IRcut_implicit << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_IRcut_explicit" << " = " << class_contribution_IRcut_explicit << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_CS" << " = " << class_contribution_CS << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_CS_virtual" << " = " << class_contribution_CS_virtual << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_CS_collinear" << " = " << class_contribution_CS_collinear << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_CS_real" << " = " << class_contribution_CS_real << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_QT" << " = " << class_contribution_QT << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_QT_virtual" << " = " << class_contribution_QT_virtual << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_QT_collinear" << " = " << class_contribution_QT_collinear << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_QT_real" << " = " << class_contribution_QT_real << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_NJ" << " = " << class_contribution_NJ << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_NJ_virtual" << " = " << class_contribution_NJ_virtual << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_NJ_collinear" << " = " << class_contribution_NJ_collinear << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_NJ_real" << " = " << class_contribution_NJ_real << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_virtual" << " = " << class_contribution_virtual << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_collinear" << " = " << class_contribution_collinear << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_real" << " = " << class_contribution_real << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_loopinduced" << " = " << class_contribution_loopinduced << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_QT_CS" << " = " << class_contribution_QT_CS << endl;
  logger << LOG_INFO << setw(35) << "class_contribution_NJ_CS" << " = " << class_contribution_NJ_CS << endl;

  logger << LOG_DEBUG << "finished" << endl;
}



/*
void contribution_set::readin_subprocess(){
  Logger logger("contribution_set::readin_subprocess");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  basic_type_parton.resize(1);
  //  type_parton.resize(1);
  //  subprocess_readin(csi.process_class, csi.subprocess, csi.decay, csi.process_type, csi.n_particle, csi.type_parton);
  //  determine_basic_subprocess();

  // determine type_parton from basic_type_parton (possibly exchanged leptons):
  determine_subprocess();

  logger << LOG_DEBUG << "finished" << endl;
}
*/
/*
void contribution_set::determine_basic_subprocess(){
  Logger logger("contribution_set::determine_basic_subprocess");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  map<string,int> code_particle;
  //  fill_code_particle(code_particle);
  //  map<int,string> name_particle;
  //  fill_name_particle(name_particle);

  logger << LOG_DEBUG << "Partonic process:" << endl;
  vector<vector<string> > particles_name(2);
  //  Logger logger("process readin");
  readin_subprocess_from_name(subprocess, basic_type_parton[0]);
  //  determine_process(subprocess, particles_name, process_type, n_particle, basic_type_parton[0], code_particle);
  logger << LOG_DEBUG << "n_particle = " << n_particle << endl;
  for (int i_p = 1; i_p < basic_type_parton[0].size(); i_p++){logger << LOG_DEBUG << "basic_type_parton[0][" << setw(2) << i_p << "] = " << setw(5) << right << basic_type_parton[0][i_p] << "   " << setw(5) << left << name_particle[basic_type_parton[0][i_p]] << endl;}

  //  int dummy;
  logger << LOG_DEBUG << "Hadronic process:" << endl;
  //  int hadron_n_particle;
  vector<vector<int> > hadron_type_parton(1);
  vector<vector<string> > hadron_particles_name(2);
  logger << LOG_DEBUG << "process_class = " << process_class << endl;

  readin_subprocess_from_name(process_class, hadron_type_parton[0]);
//  determine_process(process_class, hadron_particles_name, dummy, hadron_n_particle, hadron_type_parton[0], code_particle);
//  logger << LOG_DEBUG << "hadron_n_particle = " << hadron_n_particle << endl;

  for (int i_p = 1; i_p < hadron_type_parton[0].size(); i_p++){logger << LOG_DEBUG << "hadron_type_parton[0][" << setw(2) << i_p << "] = " << right << setw(5) << hadron_type_parton[0][i_p] << "   " << setw(5) << left << name_particle[hadron_type_parton[0][i_p]] << endl;}

  logger << LOG_DEBUG << process_class << endl;
  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG << "process_type = " << process_type << endl;

  //  logger << LOG_DEBUG << "process file" << endl;
  logger << LOG_DEBUG << "process_class           = " << process_class << endl;
  logger << LOG_DEBUG << "n_particle              = " << n_particle << endl;
  //  cout << "process_number          = " << process_number << endl;
  if (process_type == 2){
    //  if (basic_type_parton[0][0] == 0){
    for (int i2 = 1; i2 < n_particle + 3; i2++){
      logger << LOG_DEBUG << "basic_type_parton[0][" << i2 << "]                  = " << basic_type_parton[0][i2] << endl;
    }
  }
  else if (process_type == 1){
    for (int i2 = 0; i2 < n_particle + 1; i2++){
      logger << LOG_DEBUG << "basic_type_parton[0][" << i2 << "]                  = " << basic_type_parton[0][i2] << endl;
    }
  }
  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/

void contribution_set::readin_hadronic_process(){
  Logger logger("contribution_set::readin_hadronic_process");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger.newLine(LOG_DEBUG);
    
  logger << LOG_DEBUG << "Basic hadronic process:" << endl;
  vector<vector<int> > basic_hadron_type_parton(1);
  vector<vector<string> > basic_hadron_particles_name(2);
  logger << LOG_DEBUG << "basic_process_class = " << basic_process_class << endl;
  readin_subprocess_from_name(basic_process_class, basic_hadron_type_parton[0], 0);
  logger << LOG_DEBUG << "n_particle = " << n_particle << "   (temporarily)" << endl;
  logger << LOG_DEBUG << "process_type = " << process_type << "   (temporarily)" << endl;
  for (int i_p = 1; i_p < basic_hadron_type_parton[0].size(); i_p++){
    logger << LOG_DEBUG << "basic_hadron_type_parton[0][" << setw(2) << i_p << "] = " << right << setw(5) << basic_hadron_type_parton[0][i_p] << "   " << setw(5) << left << name_particle[basic_hadron_type_parton[0][i_p]] << endl;
  }

  n_jet_born = 0;
  for (int i_s = 0; i_s < basic_process_class.size(); i_s++){
    if (basic_process_class[i_s] == 'j'){n_jet_born++;}
  }
  logger << LOG_INFO << "n_jet_born = " << n_jet_born << endl;

  n_photon_born = 0;
  for (int i_s = 0; i_s < basic_process_class.size(); i_s++){
    if (basic_process_class[i_s] == 'a'){n_photon_born++;}
  }
  logger << LOG_INFO << "n_photon_born = " << n_photon_born << endl;

  logger << LOG_INFO << "New determination based on basic_hadron_type_parton[0]:" << endl;
  n_jet_born = 0;
  for (int i_p = 3; i_p < basic_hadron_type_parton[0].size(); i_p++){
    if (basic_hadron_type_parton[0][i_p] == 100){n_jet_born++;}
  }
  logger << LOG_INFO << "n_jet_born = " << n_jet_born << endl;

  n_photon_born = 0;
  for (int i_p = 0; i_p < basic_process_class.size(); i_p++){
    if (basic_hadron_type_parton[0][i_p] == 22){n_photon_born++;}
  }
  logger << LOG_INFO << "n_photon_born = " << n_photon_born << endl;


  

  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG << "Hadronic process:" << endl;
  vector<vector<int> > hadron_type_parton(1);
  vector<vector<string> > hadron_particles_name(2);
  logger << LOG_DEBUG << "process_class = " << process_class << endl;
  readin_subprocess_from_name(process_class, hadron_type_parton[0], 2);
  logger << LOG_DEBUG << "n_particle = " << n_particle << "   (temporarily)" << endl;
  logger << LOG_DEBUG << "process_type = " << process_type << "   (temporarily)" << endl;
  for (int i_p = 1; i_p < hadron_type_parton[0].size(); i_p++){
    logger << LOG_DEBUG << "hadron_type_parton[0][" << setw(2) << i_p << "] = " << right << setw(5) << hadron_type_parton[0][i_p] << "   " << setw(5) << left << name_particle[hadron_type_parton[0][i_p]] << endl;
  }

  logger.newLine(LOG_DEBUG);
  
  map<int, int> lepton_counter;
  lepton_counter[11] = 1;
  lepton_counter[13] = 2;
  lepton_counter[15] = 3;
  lepton_counter[-11] = 4;
  lepton_counter[-13] = 5;
  lepton_counter[-15] = 6;
  lepton_counter[12] = 7;
  lepton_counter[14] = 8;
  lepton_counter[16] = 9;
  lepton_counter[-12] = 10;
  lepton_counter[-14] = 11;
  lepton_counter[-16] = 12;
  
  map<int, int> lepton_inverse;
  lepton_inverse[1] = 11;
  lepton_inverse[2] = 13;
  lepton_inverse[3] = 15;
  lepton_inverse[4] = -11;
  lepton_inverse[5] = -13;
  lepton_inverse[6] = -15;
  lepton_inverse[7] = 12;
  lepton_inverse[8] = 14;
  lepton_inverse[9] = 16;
  lepton_inverse[10] = -12;
  lepton_inverse[11] = -14;
  lepton_inverse[12] = -16;
  

  //  int use_basic_process = 0;
  
  vector<int> basic_process_counter(13, 0);
  for (int i_p = 1; i_p < basic_hadron_type_parton[0].size(); i_p++){
    basic_process_counter[lepton_counter[basic_hadron_type_parton[0][i_p]]]++;
  }
  /*
    for (int i_l = 0; i_l < 13; i_l++){
      cout << "basic_process_counter[" << setw(2) <<  i_l << "] = " << setw(2) << basic_process_counter[i_l] << endl;
    }
    cout << endl;
  */
  vector<int> process_exchange_counter(13, 0);
  vector<int> process_counter(13, 0);

  //  map<int, int> lepton_exchange;
  for (int i_l = 0; i_l < 7; i_l++){
    if (i_l == 0){
      for (int j_l = 1; j_l < 7; j_l++){
	lepton_exchange[10 + j_l] = 10 + j_l;
	lepton_exchange[-10 - j_l] = -10 - j_l;
      }
    }
    else if (i_l == 1){
      lepton_exchange[11] = 13;
      lepton_exchange[12] = 14;
      lepton_exchange[13] = 11;
      lepton_exchange[14] = 12;
      lepton_exchange[15] = 15;
      lepton_exchange[16] = 16;
      lepton_exchange[-11] = -13;
      lepton_exchange[-12] = -14;
      lepton_exchange[-13] = -11;
      lepton_exchange[-14] = -12;
      lepton_exchange[-15] = -15;
      lepton_exchange[-16] = -16;
    }
    else if (i_l == 2){
      lepton_exchange[11] = 15;
      lepton_exchange[12] = 16;
      lepton_exchange[13] = 13;
      lepton_exchange[14] = 14;
      lepton_exchange[15] = 11;
      lepton_exchange[16] = 12;
      lepton_exchange[-11] = -15;
      lepton_exchange[-12] = -16;
      lepton_exchange[-13] = -13;
      lepton_exchange[-14] = -14;
      lepton_exchange[-15] = -11;
      lepton_exchange[-16] = -12;
    }
    else if (i_l == 3){
      lepton_exchange[11] = 11;
      lepton_exchange[12] = 12;
      lepton_exchange[13] = 15;
      lepton_exchange[14] = 16;
      lepton_exchange[15] = 13;
      lepton_exchange[16] = 14;
      lepton_exchange[-11] = -11;
      lepton_exchange[-12] = -12;
      lepton_exchange[-13] = -15;
      lepton_exchange[-14] = -16;
      lepton_exchange[-15] = -13;
      lepton_exchange[-16] = -14;
    }
    else if (i_l == 4){
      lepton_exchange[11] = 13;
      lepton_exchange[12] = 14;
      lepton_exchange[13] = 15;
      lepton_exchange[14] = 16;
      lepton_exchange[15] = 11;
      lepton_exchange[16] = 12;
      lepton_exchange[-11] = -13;
      lepton_exchange[-12] = -14;
      lepton_exchange[-13] = -15;
      lepton_exchange[-14] = -16;
      lepton_exchange[-15] = -11;
      lepton_exchange[-16] = -12;
    }
    else if (i_l == 5){
      lepton_exchange[11] = 15;
      lepton_exchange[12] = 16;
      lepton_exchange[13] = 11;
      lepton_exchange[14] = 12;
      lepton_exchange[15] = 13;
      lepton_exchange[16] = 14;
      lepton_exchange[-11] = -15;
      lepton_exchange[-12] = -16;
      lepton_exchange[-13] = -11;
      lepton_exchange[-14] = -12;
      lepton_exchange[-15] = -13;
      lepton_exchange[-16] = -14;
    }
    else if (i_l == 6){
      logger << LOG_FATAL << "No allowed process specified." << endl;
      //      use_basic_process = 1;
      exit(1);
    }
    
    for (int j_l = 0; j_l < 13; j_l++){process_exchange_counter[j_l] = 0;}
    for (int j_l = 0; j_l < 13; j_l++){process_counter[j_l] = 0;}

    for (int i_p = 1; i_p < hadron_type_parton[0].size(); i_p++){
      process_exchange_counter[lepton_counter[lepton_exchange[hadron_type_parton[0][i_p]]]]++;
      process_counter[lepton_counter[hadron_type_parton[0][i_p]]]++;
    }
    
      for (int j_l = 0; j_l < 13; j_l++){
	logger << LOG_DEBUG << "process_exchange_counter[" << setw(2) <<  i_l << "][" << setw(2) <<  j_l << "] = " << setw(2) << process_counter[j_l] << " --- " << setw(2) << basic_process_counter[j_l] << " = basic_process_counter[" << setw(2) <<  j_l << "]" << endl;
      }
      logger.newLine(LOG_DEBUG);
    
    if (process_exchange_counter == basic_process_counter){
      logger << LOG_DEBUG << "Permutation " << i_l << " is chosen:" << endl;
      for (int j_l = 1; j_l < 7; j_l++){
	if (lepton_exchange[10 + j_l] != 10 + j_l){logger << LOG_DEBUG << "lepton_exchange[" << setw(3) << right << 10 + j_l << "] -> " << setw(3) << right << lepton_exchange[10 + j_l] << endl;}
	if (lepton_exchange[10 + j_l] != 10 + j_l){logger << LOG_DEBUG << "lepton_exchange[" << setw(3) << right << -10 - j_l << "] -> " << setw(3) << right << lepton_exchange[-10 - j_l] << endl;}
      }
      break;
    }
  }
  //  }
  



  int first_lepton = 1;
  string new_process_class = "";
  for (int i_p = 1; i_p < hadron_type_parton[0].size(); i_p++){
    if (abs(hadron_type_parton[0][i_p]) > 10 && abs(hadron_type_parton[0][i_p]) < 17){
      if (first_lepton){
	for (int i_l = 1; i_l < 13; i_l++){
	  for (int i_r = 0; i_r < process_counter[i_l]; i_r++){
	    new_process_class = new_process_class + name_particle[lepton_inverse[i_l]]; 
	  //    [lepton_counter[lepton_exchange[hadron_type_parton[0][i_p]]]]++;
	  }
	}
	first_lepton = 0;
      }
      else {
      }
    }
    else {
      new_process_class = new_process_class + name_particle[hadron_type_parton[0][i_p]];
    }
    if ((process_type == 1 && i_p == 1) || (process_type == 2 && i_p == 2)){new_process_class = new_process_class + "-";}
  }
  new_process_class = new_process_class + "+X";

  logger << LOG_DEBUG << "new_process_class = " << new_process_class << endl;
  process_class = new_process_class;

  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void contribution_set::readin_subprocess(){
  Logger logger("contribution_set::readin_subprocess");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

   
  logger << LOG_DEBUG << "Partonic process:" << endl;
  basic_type_parton.resize(1);
  vector<vector<string> > particles_name(2);
  //  Logger logger("process readin");
  readin_subprocess_from_name(subprocess, basic_type_parton[0], 1);
  //  determine_process(subprocess, particles_name, process_type, n_particle, basic_type_parton[0], code_particle);
  logger << LOG_DEBUG << "n_particle = " << n_particle << endl;
  logger << LOG_DEBUG << "process_type = " << process_type << endl;
  for (int i_p = 1; i_p < basic_type_parton[0].size(); i_p++){
    logger << LOG_DEBUG << "basic_type_parton[0][" << setw(2) << right << i_p << "] = " << setw(5) << right << basic_type_parton[0][i_p] << "   " << setw(5) << left << name_particle[basic_type_parton[0][i_p]] << endl;
  }
  
  logger.newLine(LOG_DEBUG);
  
  type_parton = basic_type_parton;
  
  for (int i_p = 1; i_p < type_parton[0].size(); i_p++){
    if (abs(type_parton[0][i_p]) > 10 && abs(type_parton[0][i_p]) < 17){
      type_parton[0][i_p] = lepton_exchange[type_parton[0][i_p]];
    } 
    logger << LOG_DEBUG << "type_parton[0][" << setw(2) << right << i_p << "] = " << right << setw(5) << type_parton[0][i_p] << "   " << setw(5) << left << name_particle[type_parton[0][i_p]] << endl;
  }
  
  /*
    map<int, string> datname;
    fill_datname(datname);
    
    string new_subprocess = "";
  if (process_type == 2){
    new_subprocess = datname[type_parton[0][1]] + datname[type_parton[0][2]] + "_";
    for (int i = 3; i < type_parton[0].size(); i++){
      new_subprocess = new_subprocess + datname[type_parton[0][i]];
    }
  }
  else if (process_type == 1){
    new_subprocess = datname[type_parton[0][0]] + "_";
    for (int i = 1; i < type_parton[0].size(); i++){
      new_subprocess = new_subprocess + datname[type_parton[0][i]];
    }
  }
  logger << LOG_DEBUG << "old:subprocess = " << subprocess << endl;
  logger << LOG_DEBUG << "new_subprocess = " << new_subprocess << endl;
  subprocess = new_subprocess;
  */
  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


// temp_switch == 0: Basic hadronic process (Born-level particle content)
// temp_switch == 1: Partonic process
// temp_switch == 2: Hadronic process (with possibly different lepton flavours)
void contribution_set::readin_subprocess_from_name(string & this_processname, vector<int> & this_type_parton, int temp_switch){
  Logger logger("contribution_set::readin_subprocess_from_name");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  /*
  if (this_processname[this_processname.size() - 1] == 'X'){switch_to_newscheme = 0;}
  else {switch_to_newscheme = 1;}
  logger << LOG_DEBUG_VERBOSE << "switch_to_newscheme = " << switch_to_newscheme << endl;
  */

  if (switch_to_newscheme){newscheme_readin_subprocess_from_name(this_processname, this_type_parton, temp_switch);}
  else {oldscheme_readin_subprocess_from_name(this_processname, this_type_parton, temp_switch);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void contribution_set::newscheme_readin_subprocess_from_name(string & this_processname, vector<int> & this_type_parton, int temp_switch){
  Logger logger("contribution_set::readin_subprocess_from_name");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  vector<vector<string> > particles_name(2);
  int particle_counter = 0;
  int exit = 0;
  logger << LOG_DEBUG_VERBOSE << "begin 'determine_process' for " << this_processname << endl;
  // Put the spaces that preserve from possible segmentation faults:
  this_processname = this_processname + "     ";
  for (int i = 0; i < this_processname.size(); i++){
    if (particles_name[0].size() == 2 && particles_name[1].size() == 0){particle_counter = 1;}
    logger << LOG_DEBUG_VERBOSE << i << "   " << this_processname[i] << "   " << particle_counter << endl;
    // particle counter switches from incoming to outgoing particles (could happen always after the second incoming one).
    if (this_processname[i] == '_' || this_processname[i] == '-'){particle_counter++;}
    else if (this_processname[i] == 'p'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("px"); i++;}
      else {particles_name[particle_counter].push_back("p");}
    }
    else if (this_processname[i] == 'j'){particles_name[particle_counter].push_back("j");}
    // generic fermions should not appear by now:
    else if (this_processname[i] == 'f'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("fx"); i++;}
      else {particles_name[particle_counter].push_back("f");}
    }

    else if (this_processname[i] == 'd'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("dx"); i++;}
      else {particles_name[particle_counter].push_back("d");}
    }
    else if (this_processname[i] == 'u'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("ux"); i++;}
      else {particles_name[particle_counter].push_back("u");}
    }
    else if (this_processname[i] == 's'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("sx"); i++;}
      else {particles_name[particle_counter].push_back("s");}
    }
    else if (this_processname[i] == 'c'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("cx"); i++;}
      else {particles_name[particle_counter].push_back("c");}
    }
    else if (this_processname[i] == 'b'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("bx"); i++;}
      else {particles_name[particle_counter].push_back("b");}
    }
    else if (this_processname[i] == 't'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("tx"); i++;}
      else {particles_name[particle_counter].push_back("t");}
    }
    /*
    // name for tau changed to y !!!
      else if (this_processname[i + 1] == 'a')
	if      (this_processname[i + 2] == 'm'){particles_name[particle_counter].push_back("tam"); i += 2;}
	else if (this_processname[i + 2] == 'p'){particles_name[particle_counter].push_back("tap"); i += 2;}
	else {exit = 1; break;}
      else {particles_name[particle_counter].push_back("t");}
    }
    */
    else if (this_processname[i] == 'e'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("ex"); i++;}
      else {particles_name[particle_counter].push_back("e");}
    }
    else if (this_processname[i] == 'm'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("mx"); i++;}
      else {particles_name[particle_counter].push_back("m");}
    }
    else if (this_processname[i] == 'y'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("yx"); i++;}
      else {particles_name[particle_counter].push_back("y");}
    }
    else if (this_processname[i] == 'n'){
      if      (this_processname[i + 1] == 'e'){
	if (this_processname[i + 2] == 'x'){particles_name[particle_counter].push_back("nex"); i += 2;}
	else {particles_name[particle_counter].push_back("ne"); i++;}
      }
      else if (this_processname[i + 1] == 'm'){
	if (this_processname[i + 2] == 'x'){particles_name[particle_counter].push_back("nmx"); i += 2;}
	else {particles_name[particle_counter].push_back("nm"); i++;}
      }
      else if (this_processname[i + 1] == 'y'){
	if (this_processname[i + 2] == 'x'){particles_name[particle_counter].push_back("nyx"); i += 2;}
	else {particles_name[particle_counter].push_back("nt"); i++;}
      }
      else {exit = 1; break;}
    }
    else if (this_processname[i] == 'g'){particles_name[particle_counter].push_back("g");}
    else if (this_processname[i] == 'w'){
      if (this_processname[i + 1] == 'x'){particles_name[particle_counter].push_back("wx"); i++;}
      else {particles_name[particle_counter].push_back("w");}
    }
    else if (this_processname[i] == 'z'){particles_name[particle_counter].push_back("z");}
    else if (this_processname[i] == 'h'){particles_name[particle_counter].push_back("h");}
    else if (this_processname[i] == 'a'){particles_name[particle_counter].push_back("a");}
    //    else if (this_processname[i] == ' '){}
    //    else if (this_processname[i] == '+'){}
    //    else if (this_processname[i] == 'X'){}
    else {exit = 1; break;}
  }
  if (exit){} // ??? use of exit ???

  logger << LOG_DEBUG_VERBOSE << "particles_name.size() = " << particles_name.size() << endl;

  if (temp_switch == 0){  // basic hadron level
    n_particle_born = particles_name[1].size();
    logger << LOG_INFO << "n_particle_born = " << n_particle_born << endl;
  }
  if (temp_switch == 1){  // parton level
    n_particle = particles_name[1].size();
    logger << LOG_INFO << "n_particle = " << n_particle << endl;
  }

  for (int i = 0; i < particles_name.size(); i++){
    for (int j = 0; j < particles_name[i].size(); j++){
      logger << LOG_DEBUG_VERBOSE << "particles_name[" << i << "][" << j << "] = " << particles_name[i][j] << endl;
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "particles_name[0].size() = " << particles_name[0].size() << endl;
  if (particles_name[0].size() == 1){}
  else if (particles_name[0].size() == 2){this_type_parton.push_back(0);}
  logger << LOG_DEBUG_VERBOSE << "particles_name[0].size() = " << particles_name[0].size() << endl;
  //  else {exit = 1; break;}
  process_type = particles_name[0].size();
  logger << LOG_DEBUG_VERBOSE << "process_type = " << process_type << endl;
  for (int i = 0; i < 2; i++){
    logger << LOG_DEBUG_VERBOSE << "particles_name[" << i << "].size() = " << particles_name[i].size() << endl;
    for (int j = 0; j < particles_name[i].size(); j++){
    logger << LOG_DEBUG_VERBOSE << i << "   " << j << endl;
      this_type_parton.push_back(code_particle[particles_name[i][j]]);
      logger << LOG_DEBUG_VERBOSE << "particles_name[" << right << setw(2) << i << "][" << setw(2) << j << "] = " << left << setw(5) << particles_name[i][j] << "   code_particle = " << right << setw(4) << code_particle[particles_name[i][j]] << endl;
    }
  }
  if (particles_name[0].size() == 1){for (int p = 0; p < this_type_parton.size(); p++){logger << LOG_DEBUG_VERBOSE << "this_type_parton[" << right << setw(2) << p << "] = " << right << this_type_parton[p] << endl;}}
  else if (particles_name[0].size() == 2){for (int p = 1; p < this_type_parton.size(); p++){logger << LOG_DEBUG_VERBOSE << "this_type_parton[" << right << setw(2) << p << "] = " << right << this_type_parton[p] << endl;}}
  // remove the spaces that preserve from possible segmentation faults:
  this_processname.erase(this_processname.end() - 5, this_processname.end());
  logger << LOG_DEBUG_VERBOSE << "end 'determine_process' for " << this_processname << endl;
}

void contribution_set::oldscheme_readin_subprocess_from_name(string & this_processname, vector<int> & this_type_parton, int temp_switch){
  Logger logger("oldscheme_readin_subprocess_from_name");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  vector<vector<string> > particles_name(2);
  int particle_counter = 0;
  int exit = 0;
  logger << LOG_DEBUG_VERBOSE << "begin 'determine_process' for " << this_processname << endl;
  this_processname = this_processname + "     ";
  for (int i = 0; i < this_processname.size(); i++){
    logger << LOG_DEBUG_VERBOSE << i << "   " << this_processname[i] << "   " << particle_counter << endl;
    if (this_processname[i] == '_' || this_processname[i] == '-'){particle_counter++;}
    else if (this_processname[i] == 'p'){
      if (this_processname[i + 1] == '~'){particles_name[particle_counter].push_back("p~"); i++;}
      else {particles_name[particle_counter].push_back("p");}
    }
    else if (this_processname[i] == 'j'){particles_name[particle_counter].push_back("j");}
    /*
    else if (this_processname[i] == 'j'){
      if (this_processname[i + 1] == 'e' && this_processname[i + 2] == 't'){particles_name[particle_counter].push_back("jet"); i += 2;}
      else {exit = 1; break;}
    }
    else if (this_processname[i] == '2'){
      if (this_processname[i + 1] == 'j' && this_processname[i + 2] == 'e' && this_processname[i + 3] == 't' && this_processname[i + 4] == 's'){
	for (int j = 0; j < 2; j++){particles_name[particle_counter].push_back("jet");}
	i += 4;
      }
      else {exit = 1; break;}
    }
    else if (this_processname[i] == '3'){
      if (this_processname[i + 1] == 'j' && this_processname[i + 2] == 'e' && this_processname[i + 3] == 't' && this_processname[i + 4] == 's'){
	for (int j = 0; j < 3; j++){particles_name[particle_counter].push_back("jet");}
	i += 4;
      }
      else {exit = 1; break;}
    }
    else if (this_processname[i] == '4'){
      if (this_processname[i + 1] == 'j' && this_processname[i + 2] == 'e' && this_processname[i + 3] == 't' && this_processname[i + 4] == 's'){
	for (int j = 0; j < 4; j++){particles_name[particle_counter].push_back("jet");}
	i += 4;
      }
      else {exit = 1; break;}
    }
    */
    else if (this_processname[i] == 'f'){
      if (this_processname[i + 1] == '~'){particles_name[particle_counter].push_back("f~"); i++;}
      else {particles_name[particle_counter].push_back("f");}
    }

    else if (this_processname[i] == 'd'){
      if (this_processname[i + 1] == '~'){particles_name[particle_counter].push_back("d~"); i++;}
      else {particles_name[particle_counter].push_back("d");}
    }
    else if (this_processname[i] == 'u'){
      if (this_processname[i + 1] == '~'){particles_name[particle_counter].push_back("u~"); i++;}
      else {particles_name[particle_counter].push_back("u");}
    }
    else if (this_processname[i] == 's'){
      if (this_processname[i + 1] == '~'){particles_name[particle_counter].push_back("s~"); i++;}
      else {particles_name[particle_counter].push_back("s");}
    }
    else if (this_processname[i] == 'c'){
      if (this_processname[i + 1] == '~'){particles_name[particle_counter].push_back("c~"); i++;}
      else {particles_name[particle_counter].push_back("c");}
    }
    else if (this_processname[i] == 'b'){
      if (this_processname[i + 1] == '~'){particles_name[particle_counter].push_back("b~"); i++;}
      else {particles_name[particle_counter].push_back("b");}
    }
    else if (this_processname[i] == 't'){
      if (this_processname[i + 1] == '~'){particles_name[particle_counter].push_back("t~"); i++;}
      else if (this_processname[i + 1] == 'a')
	if      (this_processname[i + 2] == 'm'){particles_name[particle_counter].push_back("tam"); i += 2;}
	else if (this_processname[i + 2] == 'p'){particles_name[particle_counter].push_back("tap"); i += 2;}
	else {exit = 1; break;}
      else {particles_name[particle_counter].push_back("t");}
    }
    else if (this_processname[i] == 'e'){
      if      (this_processname[i + 1] == 'm'){particles_name[particle_counter].push_back("em"); i++;}
      else if (this_processname[i + 1] == 'p'){particles_name[particle_counter].push_back("ep"); i++;}
      else {exit = 1; break;}
    }
    else if (this_processname[i] == 'm'){
      if      (this_processname[i + 1] == 'u' && this_processname[i + 2] == 'm'){particles_name[particle_counter].push_back("mum"); i += 2;}
      else if (this_processname[i + 1] == 'u' && this_processname[i + 2] == 'p'){particles_name[particle_counter].push_back("mup"); i += 2;}
      else {exit = 1; break;}
    }
    else if (this_processname[i] == 'v'){
      if      (this_processname[i + 1] == 'e'){
	if (this_processname[i + 2] == '~'){particles_name[particle_counter].push_back("ve~"); i += 2;}
	else {particles_name[particle_counter].push_back("ve"); i++;}
      }
      else if (this_processname[i + 1] == 'm'){
	if (this_processname[i + 2] == '~'){particles_name[particle_counter].push_back("vm~"); i += 2;}
	else {particles_name[particle_counter].push_back("vm"); i++;}
      }
      else if (this_processname[i + 1] == 't'){
	if (this_processname[i + 2] == '~'){particles_name[particle_counter].push_back("vt~"); i += 2;}
	else {particles_name[particle_counter].push_back("vt"); i++;}
      }
      else {exit = 1; break;}
    }
    else if (this_processname[i] == 'g'){particles_name[particle_counter].push_back("g");}
    else if (this_processname[i] == 'w'){
      if      (this_processname[i + 1] == 'm'){particles_name[particle_counter].push_back("wm"); i++;}
      else if (this_processname[i + 1] == 'p'){particles_name[particle_counter].push_back("wp"); i++;}
      else {exit = 1; break;}
    }
    else if (this_processname[i] == 'z'){particles_name[particle_counter].push_back("z");}
    else if (this_processname[i] == 'h'){particles_name[particle_counter].push_back("h");}
    else if (this_processname[i] == 'a'){particles_name[particle_counter].push_back("a");}
    else if (this_processname[i] == ' '){}
    else if (this_processname[i] == '+'){}
    else if (this_processname[i] == 'X'){}
    else {exit = 1; break;}
  }
  if (exit){} // ??? use of exit ???

  logger << LOG_DEBUG_VERBOSE << "particles_name.size() = " << particles_name.size() << endl;
  //  n_particle = particles_name[1].size();
  //  logger << LOG_DEBUG_VERBOSE << "n_particle = " << n_particle << endl;

  if (temp_switch == 0){  // basic hadron level
    n_particle_born = particles_name[1].size();
    logger << LOG_INFO << "n_particle_born = " << n_particle_born << endl;
  }
  if (temp_switch == 1){  // parton level
    n_particle = particles_name[1].size();
    logger << LOG_INFO << "n_particle = " << n_particle << endl;
  }

  for (int i = 0; i < particles_name.size(); i++){
    for (int j = 0; j < particles_name[i].size(); j++){
      logger << LOG_DEBUG_VERBOSE << "particles_name[" << i << "][" << j << "] = " << particles_name[i][j] << endl;
    }
  }
  logger << LOG_DEBUG_VERBOSE << "particles_name[0].size() = " << particles_name[0].size() << endl;
  if (particles_name[0].size() == 1){}
  else if (particles_name[0].size() == 2){this_type_parton.push_back(0);}
  logger << LOG_DEBUG_VERBOSE << "particles_name[0].size() = " << particles_name[0].size() << endl;
  //  else {exit = 1; break;}
  process_type = particles_name[0].size();
  logger << LOG_DEBUG_VERBOSE << "process_type = " << process_type << endl;
  for (int i = 0; i < 2; i++){
    logger << LOG_DEBUG_VERBOSE << "particles_name[" << i << "].size() = " << particles_name[i].size() << endl;
    for (int j = 0; j < particles_name[i].size(); j++){
    logger << LOG_DEBUG_VERBOSE << i << "   " << j << endl;
      this_type_parton.push_back(code_particle[particles_name[i][j]]);
      logger << LOG_DEBUG_VERBOSE << "particles_name[" << right << setw(2) << i << "][" << setw(2) << j << "] = " << left << setw(5) << particles_name[i][j] << "   code_particle = " << right << setw(4) << code_particle[particles_name[i][j]] << endl;
    }
  }
  if (particles_name[0].size() == 1){for (int p = 0; p < this_type_parton.size(); p++){logger << LOG_DEBUG_VERBOSE << "this_type_parton[" << right << setw(2) << p << "] = " << right << this_type_parton[p] << endl;}}
  else if (particles_name[0].size() == 2){for (int p = 1; p < this_type_parton.size(); p++){logger << LOG_DEBUG_VERBOSE << "this_type_parton[" << right << setw(2) << p << "] = " << right << this_type_parton[p] << endl;}}
  this_processname.erase(this_processname.end() - 5, this_processname.end());
  logger << LOG_DEBUG_VERBOSE << "end 'determine_process' for " << this_processname << endl;
}


void contribution_set::fill_code_particle(){
  Logger logger("fill_code_particle");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  if (switch_to_newscheme){newscheme_fill_code_particle();}
  else {oldscheme_fill_code_particle();}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void contribution_set::newscheme_fill_code_particle(){
  Logger logger("newscheme_fill_code_particle");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  code_particle["d"]   = 1;
  code_particle["u"]   = 2;
  code_particle["s"]   = 3;
  code_particle["c"]   = 4;
  code_particle["b"]   = 5;
  code_particle["t"]   = 6;
  code_particle["dx"]  = -1;
  code_particle["ux"]  = -2;
  code_particle["sx"]  = -3;
  code_particle["cx"]  = -4;
  code_particle["bx"]  = -5;
  code_particle["tx"]  = -6;
  code_particle["e"]   = 11;
  code_particle["m"]   = 13;
  code_particle["y"]   = 15;
  code_particle["ex"]  = -11;
  code_particle["mx"]  = -13;
  code_particle["yx"]  = -15;
  code_particle["ne"]  = 12;
  code_particle["nm"]  = 14;
  code_particle["nt"]  = 16;
  code_particle["nex"] = -12;
  code_particle["nmx"] = -14;
  code_particle["ntx"] = -16;
  code_particle["g"]   = 0;
  code_particle["a"]   = 22;
  code_particle["ax"]  = -22;
  code_particle["wx"]  = -24;
  code_particle["w"]   = 24;
  code_particle["z"]   = 23;
  code_particle["h"]   = 25;
  code_particle["p"]   = 101;
  code_particle["px"]  = -101;
  code_particle["j"]   = 100;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void contribution_set::oldscheme_fill_code_particle(){
  Logger logger("oldscheme_fill_code_particle");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  code_particle["d"]   = 1;
  code_particle["u"]   = 2;
  code_particle["s"]   = 3;
  code_particle["c"]   = 4;
  code_particle["b"]   = 5;
  code_particle["t"]   = 6;
  code_particle["d~"]  = -1;
  code_particle["u~"]  = -2;
  code_particle["s~"]  = -3;
  code_particle["c~"]  = -4;
  code_particle["b~"]  = -5;
  code_particle["t~"]  = -6;
  code_particle["em"]  = 11;
  code_particle["mum"] = 13;
  code_particle["tam"] = 15;
  code_particle["ep"]  = -11;
  code_particle["mup"] = -13;
  code_particle["tap"] = -15;
  code_particle["ve"]  = 12;
  code_particle["vm"]  = 14;
  code_particle["vt"]  = 16;
  code_particle["ve~"] = -12;
  code_particle["vm~"] = -14;
  code_particle["vt~"] = -16;
  code_particle["g"]   = 0;
  code_particle["a"]   = 22;
  code_particle["a~"]  = -22;
  code_particle["wp"]  = -24;
  code_particle["wm"]  = 24;
  code_particle["z"]   = 23;
  code_particle["h"]   = 25;
  code_particle["p"]  = 101;
  code_particle["p~"]  = -101;
  code_particle["j"]  = 100;
  //  code_particle["jet"]  = 100;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void contribution_set::fill_name_particle(){
  Logger logger("fill_name_particle");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  if (switch_to_newscheme){newscheme_fill_name_particle();}
  else {oldscheme_fill_name_particle();}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void contribution_set::newscheme_fill_name_particle(){
  Logger logger("newscheme_fill_name_particle");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  name_particle[  0] = "g";
  name_particle[ 22] = "a";
  name_particle[-22] = "ax";
  name_particle[-24] = "wx";
  name_particle[ 24] = "w";
  name_particle[ 23] = "z";
  name_particle[ 25] = "h";
  name_particle[  1] = "d";
  name_particle[  2] = "u";
  name_particle[  3] = "s";
  name_particle[  4] = "c";
  name_particle[  5] = "b";
  name_particle[  6] = "t";
  name_particle[ -1] = "dx";
  name_particle[ -2] = "ux";
  name_particle[ -3] = "sx";
  name_particle[ -4] = "cx";
  name_particle[ -5] = "bx";
  name_particle[ -6] = "tx";
  name_particle[ 11] = "e";
  name_particle[ 12] = "ne";
  name_particle[ 13] = "m";
  name_particle[ 14] = "nm";
  name_particle[ 15] = "y";
  name_particle[ 16] = "ny";
  name_particle[-11] = "ex";
  name_particle[-12] = "nex";
  name_particle[-13] = "mx";
  name_particle[-14] = "nmx";
  name_particle[-15] = "yx";
  name_particle[-16] = "nyx";
  name_particle[101] = "p";
  name_particle[-101] = "px";
  name_particle[100] = "j";

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void contribution_set::oldscheme_fill_name_particle(){
  Logger logger("oldscheme_fill_name_particle");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  name_particle[  0] = "g";
  name_particle[ 22] = "a";
  name_particle[-22] = "a~";
  name_particle[-24] = "wp";
  name_particle[ 24] = "wm";
  name_particle[ 23] = "z";
  name_particle[ 25] = "h";
  name_particle[  1] = "d";
  name_particle[  2] = "u";
  name_particle[  3] = "s";
  name_particle[  4] = "c";
  name_particle[  5] = "b";
  name_particle[  6] = "t";
  name_particle[ -1] = "d~";
  name_particle[ -2] = "u~";
  name_particle[ -3] = "s~";
  name_particle[ -4] = "c~";
  name_particle[ -5] = "b~";
  name_particle[ -6] = "t~";
  name_particle[ 11] = "em";
  name_particle[ 12] = "ve";
  name_particle[ 13] = "mum";
  name_particle[ 14] = "vm";
  name_particle[ 15] = "tam";
  name_particle[ 16] = "vt";
  name_particle[-11] = "ep";
  name_particle[-12] = "ve~";
  name_particle[-13] = "mup";
  name_particle[-14] = "vm~";
  name_particle[-15] = "tap";
  name_particle[-16] = "vt~";
  name_particle[101] = "p";
  name_particle[-101] = "p~";
  name_particle[100] = "j";
  //  name_particle[100] = "jet";

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
