#ifndef DEFINEPARTICLESET_H
#define DEFINEPARTICLESET_H

#include <map>
#include <string>
#include <vector>

using namespace std;

class define_particle_set{
 private:

 public:
////////////////////
//  constructors  //
////////////////////
  define_particle_set();

  int n_observed_min;
  int n_observed_max;
  double define_pT;
  double define_ET;
  double define_eta;
  double define_y;

  int n_partonlevel;

};
#endif
