#include "../include/classes.cxx"

////////////////////
//  constructors  //
////////////////////
define_particle_set::define_particle_set(){
 Logger logger("define_particle_set::define_particle_set");
  logger << LOG_DEBUG << "started" << endl;

  n_observed_min = 0;
  n_observed_max = 99;
  define_pT = 0.;
  define_ET = 0.;
  define_eta = 1.e99;
  define_y = 1.e99;

  n_partonlevel = 0;

  logger << LOG_DEBUG << "finished" << endl;
}
