#ifndef DDISTRIBUTION_H
#define DDISTRIBUTION_H

#include <string>

#include "xdistribution.h"

using namespace std;

class dddistribution {
 private:
  // coordinates, hidden
  /*
  string priv_name;
  int priv_d_1;
  int priv_d_2;
  int priv_n_bins;
  double priv_step;
  xdistribution priv_distribution_1;
  xdistribution priv_distribution_2;
  */

 public:  
  // constructors
  dddistribution();
  dddistribution(string _name, int _d_1, int _d_2, observable_set *_oset);
  //, vector<xdistribution> dat
  //  dddistribution(string name, int d_1, int d_2, xdistribution dat1, xdistribution dat2);

  void determineBin();
  void initialization_distribution_bin();

  observable_set *oset;

  string name;
  int d_1;
  int d_2;
  int n_bins;
  double step;
  xdistribution distribution_1;
  xdistribution distribution_2;
  //  xdistribution *distribution_1;
  //  xdistribution *distribution_2;

  int mirror_type;
  // -1 -> no contribution from rotated psp (symmetric IS - independent of symmetry of distribution)
  //  0 -> contribution from rotated psp (non-symmetric IS and fully non-symmetric distribution)
  // +1 -> contribution from rotated psp (non-symmetric IS and at least partially symmetric distribution)
  
  vector<int> bin;
  //  bin   (contains the relevant bin number of the respective distribution (the smallest bin number in case of cumulated distributions)):
  //  1st index -> no_phasespace (i_a)
  //  entry     -> no_bin

  vector<int> bin_max;
  //  bin_max   contains the highest relevant bin number (i_b < dat[i_d].bin_max[i_a]) of the respective cumulated distribution (irrelevant for normal distributions):
  //  1st index -> no_phasespace (i_a)  
  //  entry     -> no_bin (+1)


  
  vector<int> mirror_bin;
  //  bin   (contains the relevant bin number of the respective distribution (the smallest bin number in case of cumulated distributions)):
  //  1st index -> no_phasespace (i_a)
  //  entry     -> no_bin

  vector<int> mirror_bin_max;
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
  
  //  to be introduced !!!
  vector<int> distribution_no_mirror_bin;
  vector<int> distribution_no_mirror_bin_max;
  vector<int> distribution_no_mirror_bin_all;



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


  
  // not used yet, I guess... could be used to set several intervals in CUMULATIVE distributions
  vector<vector<int> > advanced_bin;
  vector<vector<int> > advanced_bin_max;



  friend std::ostream & operator << ( std::ostream &, const dddistribution &);
};

#endif
