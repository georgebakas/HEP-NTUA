#include "../include/classes.cxx"
//#include "../include/definitions.cxx"
#include "../include/definitions.observable.set.cxx"

// !!! cumulative distributions and dddistributions are probably not working properly with qT subtraction !!!

xdistribution get_fake_distribution_from_dddistribution(dddistribution & dddist){
  Logger logger("get_fake_distribution_from_dddistribution");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  xdistribution temp;

  temp.xdistribution_name = dddist.name;
  temp.n_bins = dddist.n_bins;

  temp.bin_edge.resize(temp.n_bins + 1);
  temp.bin_width.resize(temp.n_bins);
  for (int i_b1 = 0; i_b1 < dddist.distribution_1.n_bins; i_b1++){
    for (int i_b2 = 0; i_b2 < dddist.distribution_2.n_bins; i_b2++){
      int i_b = i_b1 * dddist.distribution_2.n_bins + i_b2;
      temp.bin_edge[i_b] = i_b1 * 1000 + i_b2;//dddist.distribution_1.bin_edge[i_b1] * dddist.distribution_2.bin_edge[i_b1];
      temp.bin_width[i_b] = dddist.distribution_1.bin_width[i_b1] * dddist.distribution_2.bin_width[i_b2];
    }
  }
  temp.bin_edge[temp.n_bins] = dddist.distribution_1.n_bins * 1000 + dddist.distribution_2.n_bins;

  return temp;
}



