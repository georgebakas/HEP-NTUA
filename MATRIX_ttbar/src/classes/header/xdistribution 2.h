#ifndef XDISTRIBUTION_H
#define XDISTRIBUTION_H

#include <string>

class observable_set;

// #include "observable.set.h"

using namespace std;

typedef enum {
    CUMULATIVE_NONE = 0,
    CUMULATIVE_LOWER,
    CUMULATIVE_UPPER,
    CUMULATIVE_BOTH
} TypeCumulative;

struct ParticleID {
  int type;
  int number;
  
  bool required_exists;
  bool required_definition;
  
  bool is_present;
  
  // number of particles of the same type
  int multiplicity;
  
  fourvector momentum;
};

class xdistribution {
 private:
  double performBinning(double value);
  void computeObservable(double &observable, double &observable_max);
  bool checkRequirements(int i_a);
  void fillMomenta(int i_a);
   
  observable_set *oset;
  vector<double> observable;
  vector<double> observable_max;
  vector<vector<ParticleID> > particles;
  vector<fourvector> reconstructedParticles;
  
 public:  
  // constructors
  xdistribution();
  xdistribution(string _name, string _xdistribution_type, vector<vector<int> > _all_particle, vector<vector<int> > _all_particle_group, vector<vector<int> > _required_subparticle, double _start, double _end, int _n_bins, double _step, string _edges, string _type_binning, observable_set *_oset=NULL);
  
  void determineBin();

  void initialization_distribution_bin();

//   int xdistribution_number; 
  vector<vector<int> > all_particle;
  vector<vector<int> > all_particle_group;
  vector<vector<int> > required_subparticle;
  vector<int> required_particle;
  double start;
  double end;
  int n_bins;
  double step;
  string xdistribution_name;
  string xdistribution_type;
  int symm;

  string edges;
  string type_binning;
  int type_cumulative;
  vector<double> bin_edge;
  vector<double> bin_width;
  
  vector<int> bin;
  //  bin   (contains the relevant bin number of the respective distribution (the smallest bin number in case of cumulated distributions)):
  //  1st index -> no_phasespace (i_a)
  //  entry     -> no_bin

  vector<int> bin_max;
  //  bin_max   contains the highest relevant bin number (i_b < dat[i_d].bin_max[i_a]) of the respective cumulated distribution (irrelevant for normal distributions):
  //  1st index -> no_phasespace (i_a)  
  //  entry     -> no_bin (+1)

  vector<int> distribution_no_bin;
  //  distribution_no_bin   (contains all bin numbers which are relevant for the present set of phase-space points):   
  //  1st index -> counter of relevant bin numbers   
  //  entry     -> no_bin

  vector<int> distribution_no_bin_max;
  //  distribution_no_bin_max   (contains all bin_max numbers which are relevant for the present set of phase-space points):
  //  1st index -> counter of relevant bin numbers
  //  entry     -> no_bin (+1)

  vector<int> distribution_no_bin_all;
  //  distribution_no_bin_all
  //  1st index -> counter of relevant bin/bin_max numbers   
  //  entry     -> no_bin (in ascending (?) order, appearing only once each)

  vector<vector<int> > distribution_no_bin_phasespace;
  //  distribution_no_bin_phasespace   each entry contains list with all contributing phase-spaces (which might still be cut at different qTcut values)
  //  1st index -> number of bin from  distribution_no_bin_all (bin accesible there)
  //  2nd index -> counter of phasespaces (i_a) with that bin value
  //  entry     -> no_phasespace (i_a)
  //  Is  distribution_no_bin_phasespace  really neeed ???

  vector<vector<int> > distribution_no_bin_max_phasespace;
  //  ???
  //  Is  distribution_no_bin_max_phasespace  really neeed ???

  vector<vector<vector<int> > > distribution_no_bin_no_qTcut_phasespace;
  //  distribution_no_bin_no_qTcut_phasespace
  //  1st entry -> bin value from  distribution_no_bin_all
  //  2nd entry -> cut value from  oset->distribution_no_qTcut
  //  3rd entry -> counter of no of phasespaces (i_a) with that bin and qTcut value
  //  entry     -> no of phasespace (i_a)


  TypeCumulative typeCumulative;

  friend std::ostream & operator << ( std::ostream &, const xdistribution &);
};

#endif
