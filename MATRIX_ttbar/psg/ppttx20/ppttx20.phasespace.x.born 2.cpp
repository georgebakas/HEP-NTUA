#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_optimize_minv_born(phasespace_set & psi){
  static Logger logger("ppttx20_optimize_minv_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

#include "user/specify.optimize.phasespace.minv.born.cxx"

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

int ppttx20_determination_MCchannels_born(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_determination_MCchannels_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int n_channel = 0;
  if      (psi_no_map[no_ps] == 0){n_channel = 3;}
  else if (psi_no_map[no_ps] == 1){n_channel = 3;}
  else if (psi_no_map[no_ps] == 2){n_channel = 1;}
  else if (psi_no_map[no_ps] == 3){n_channel = 1;}
  else if (psi_no_map[no_ps] == 4){n_channel = 1;}
  return n_channel;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_ac_tau_psp_born(int no_ps, vector<int> & tau_MC_map, phasespace_set & psi){
  static Logger logger("ppttx20_ac_tau_psp_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if      (psi.no_map[no_ps] == 0){}
  else if (psi.no_map[no_ps] == 1){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 2){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 3){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 4){
    tau_MC_map.push_back(0);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_ax_psp_born(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (psi_no_map[no_ps] == 0){}
  else if (psi_no_map[no_ps] == 1){ppttx20_ax_psp_200_gg_ttx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 2){ppttx20_ax_psp_200_ddx_ttx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 3){ppttx20_ax_psp_200_uux_ttx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 4){ppttx20_ax_psp_200_bbx_ttx(no_ps, psi);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
