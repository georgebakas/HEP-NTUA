#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_optimize_minv_real(phasespace_set & psi){
  static Logger logger("ppttx20_optimize_minv_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

#include "user/specify.optimize.phasespace.minv.real.cxx"

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

int ppttx20_determination_MCchannels_real(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_determination_MCchannels_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int n_channel = 0;
  if      (psi_no_map[no_ps] ==  0){n_channel = 18;}
  else if (psi_no_map[no_ps] ==  1){n_channel = 15;}
  else if (psi_no_map[no_ps] ==  2){n_channel = 5;}
  else if (psi_no_map[no_ps] ==  3){n_channel = 5;}
  else if (psi_no_map[no_ps] ==  4){n_channel = 5;}
  else if (psi_no_map[no_ps] ==  5){n_channel = 5;}
  else if (psi_no_map[no_ps] ==  6){n_channel = 5;}
  else if (psi_no_map[no_ps] ==  7){n_channel = 5;}
  else if (psi_no_map[no_ps] ==  8){n_channel = 5;}
  else if (psi_no_map[no_ps] ==  9){n_channel = 5;}
  else if (psi_no_map[no_ps] == 10){n_channel = 5;}
  return n_channel;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_ac_tau_psp_real(int no_ps, vector<int> & tau_MC_map, phasespace_set & psi){
  static Logger logger("ppttx20_ac_tau_psp_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if      (psi.no_map[no_ps] == 0){}
  else if (psi.no_map[no_ps] ==  1){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] ==  2){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] ==  3){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] ==  4){
    tau_MC_map.push_back(0);
    tau_MC_map.push_back(-5);
  }
  else if (psi.no_map[no_ps] ==  5){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] ==  6){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] ==  7){
    tau_MC_map.push_back(0);
    tau_MC_map.push_back(-5);
  }
  else if (psi.no_map[no_ps] ==  8){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] ==  9){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 10){
    tau_MC_map.push_back(0);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_ax_psp_real(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (psi_no_map[no_ps] == 0){}
  else if (psi_no_map[no_ps] ==  1){ppttx20_ax_psp_300_gg_ttxg(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  2){ppttx20_ax_psp_300_gd_ttxd(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  3){ppttx20_ax_psp_300_gu_ttxu(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  4){ppttx20_ax_psp_300_gb_ttxb(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  5){ppttx20_ax_psp_300_gdx_ttxdx(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  6){ppttx20_ax_psp_300_gux_ttxux(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  7){ppttx20_ax_psp_300_gbx_ttxbx(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  8){ppttx20_ax_psp_300_ddx_ttxg(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  9){ppttx20_ax_psp_300_uux_ttxg(no_ps, psi);}
  else if (psi_no_map[no_ps] == 10){ppttx20_ax_psp_300_bbx_ttxg(no_ps, psi);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
