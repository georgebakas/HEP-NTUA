#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_optimize_minv_doublereal(phasespace_set & psi){
  static Logger logger("ppttx20_optimize_minv_doublereal");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

#include "user/specify.optimize.phasespace.minv.doublereal.cxx"

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

int ppttx20_determination_MCchannels_doublereal(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_determination_MCchannels_doublereal");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int n_channel = 0;
  if      (psi_no_map[no_ps] ==  0){n_channel = 160;}
  else if (psi_no_map[no_ps] ==  1){n_channel = 105;}
  else if (psi_no_map[no_ps] ==  2){n_channel = 35;}
  else if (psi_no_map[no_ps] ==  3){n_channel = 35;}
  else if (psi_no_map[no_ps] ==  4){n_channel = 35;}
  else if (psi_no_map[no_ps] ==  5){n_channel = 35;}
  else if (psi_no_map[no_ps] ==  6){n_channel = 35;}
  else if (psi_no_map[no_ps] ==  7){n_channel = 35;}
  else if (psi_no_map[no_ps] ==  8){n_channel = 35;}
  else if (psi_no_map[no_ps] ==  9){n_channel = 35;}
  else if (psi_no_map[no_ps] == 10){n_channel = 35;}
  else if (psi_no_map[no_ps] == 11){n_channel = 14;}
  else if (psi_no_map[no_ps] == 12){n_channel = 7;}
  else if (psi_no_map[no_ps] == 13){n_channel = 7;}
  else if (psi_no_map[no_ps] == 14){n_channel = 7;}
  else if (psi_no_map[no_ps] == 16){n_channel = 7;}
  else if (psi_no_map[no_ps] == 17){n_channel = 35;}
  else if (psi_no_map[no_ps] == 18){n_channel = 14;}
  else if (psi_no_map[no_ps] == 19){n_channel = 7;}
  else if (psi_no_map[no_ps] == 20){n_channel = 7;}
  else if (psi_no_map[no_ps] == 21){n_channel = 7;}
  else if (psi_no_map[no_ps] == 22){n_channel = 7;}
  else if (psi_no_map[no_ps] == 23){n_channel = 7;}
  else if (psi_no_map[no_ps] == 25){n_channel = 7;}
  else if (psi_no_map[no_ps] == 27){n_channel = 7;}
  else if (psi_no_map[no_ps] == 28){n_channel = 7;}
  else if (psi_no_map[no_ps] == 29){n_channel = 14;}
  else if (psi_no_map[no_ps] == 30){n_channel = 7;}
  else if (psi_no_map[no_ps] == 31){n_channel = 7;}
  else if (psi_no_map[no_ps] == 32){n_channel = 7;}
  else if (psi_no_map[no_ps] == 34){n_channel = 35;}
  else if (psi_no_map[no_ps] == 35){n_channel = 7;}
  else if (psi_no_map[no_ps] == 36){n_channel = 14;}
  else if (psi_no_map[no_ps] == 37){n_channel = 7;}
  else if (psi_no_map[no_ps] == 38){n_channel = 7;}
  else if (psi_no_map[no_ps] == 39){n_channel = 7;}
  else if (psi_no_map[no_ps] == 40){n_channel = 7;}
  else if (psi_no_map[no_ps] == 42){n_channel = 7;}
  else if (psi_no_map[no_ps] == 43){n_channel = 7;}
  else if (psi_no_map[no_ps] == 44){n_channel = 14;}
  else if (psi_no_map[no_ps] == 45){n_channel = 7;}
  else if (psi_no_map[no_ps] == 46){n_channel = 7;}
  else if (psi_no_map[no_ps] == 47){n_channel = 35;}
  else if (psi_no_map[no_ps] == 48){n_channel = 7;}
  else if (psi_no_map[no_ps] == 49){n_channel = 7;}
  else if (psi_no_map[no_ps] == 50){n_channel = 14;}
  else if (psi_no_map[no_ps] == 51){n_channel = 14;}
  else if (psi_no_map[no_ps] == 52){n_channel = 7;}
  else if (psi_no_map[no_ps] == 53){n_channel = 7;}
  else if (psi_no_map[no_ps] == 54){n_channel = 7;}
  else if (psi_no_map[no_ps] == 56){n_channel = 7;}
  else if (psi_no_map[no_ps] == 57){n_channel = 14;}
  else if (psi_no_map[no_ps] == 58){n_channel = 7;}
  else if (psi_no_map[no_ps] == 59){n_channel = 7;}
  else if (psi_no_map[no_ps] == 60){n_channel = 14;}
  return n_channel;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_ac_tau_psp_doublereal(int no_ps, vector<int> & tau_MC_map, phasespace_set & psi){
  static Logger logger("ppttx20_ac_tau_psp_doublereal");
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
    tau_MC_map.push_back(-5);
  }
  else if (psi.no_map[no_ps] == 11){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 12){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 13){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 14){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 16){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 17){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 18){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 19){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 20){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 21){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 22){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 23){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 25){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 27){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 28){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 29){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 30){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 31){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 32){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 34){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 35){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 36){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 37){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 38){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 39){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 40){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 42){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 43){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 44){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 45){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 46){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 47){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 48){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 49){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 50){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 51){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 52){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 53){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 54){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 56){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 57){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 58){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 59){
    tau_MC_map.push_back(0);
  }
  else if (psi.no_map[no_ps] == 60){
    tau_MC_map.push_back(0);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_ax_psp_doublereal(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_doublereal");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (psi_no_map[no_ps] == 0){}
  else if (psi_no_map[no_ps] ==  1){ppttx20_ax_psp_400_gg_ttxgg(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  2){ppttx20_ax_psp_400_gg_ttxddx(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  3){ppttx20_ax_psp_400_gg_ttxuux(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  4){ppttx20_ax_psp_400_gg_ttxbbx(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  5){ppttx20_ax_psp_400_gd_ttxgd(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  6){ppttx20_ax_psp_400_gu_ttxgu(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  7){ppttx20_ax_psp_400_gb_ttxgb(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  8){ppttx20_ax_psp_400_gdx_ttxgdx(no_ps, psi);}
  else if (psi_no_map[no_ps] ==  9){ppttx20_ax_psp_400_gux_ttxgux(no_ps, psi);}
  else if (psi_no_map[no_ps] == 10){ppttx20_ax_psp_400_gbx_ttxgbx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 11){ppttx20_ax_psp_400_dd_ttxdd(no_ps, psi);}
  else if (psi_no_map[no_ps] == 12){ppttx20_ax_psp_400_du_ttxdu(no_ps, psi);}
  else if (psi_no_map[no_ps] == 13){ppttx20_ax_psp_400_ds_ttxds(no_ps, psi);}
  else if (psi_no_map[no_ps] == 14){ppttx20_ax_psp_400_dc_ttxdc(no_ps, psi);}
  else if (psi_no_map[no_ps] == 16){ppttx20_ax_psp_400_db_ttxdb(no_ps, psi);}
  else if (psi_no_map[no_ps] == 17){ppttx20_ax_psp_400_ddx_ttxgg(no_ps, psi);}
  else if (psi_no_map[no_ps] == 18){ppttx20_ax_psp_400_ddx_ttxddx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 19){ppttx20_ax_psp_400_ddx_ttxuux(no_ps, psi);}
  else if (psi_no_map[no_ps] == 20){ppttx20_ax_psp_400_ddx_ttxssx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 21){ppttx20_ax_psp_400_ddx_ttxccx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 22){ppttx20_ax_psp_400_ddx_ttxbbx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 23){ppttx20_ax_psp_400_dux_ttxdux(no_ps, psi);}
  else if (psi_no_map[no_ps] == 25){ppttx20_ax_psp_400_dsx_ttxdsx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 27){ppttx20_ax_psp_400_dcx_ttxdcx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 28){ppttx20_ax_psp_400_dbx_ttxdbx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 29){ppttx20_ax_psp_400_uu_ttxuu(no_ps, psi);}
  else if (psi_no_map[no_ps] == 30){ppttx20_ax_psp_400_uc_ttxuc(no_ps, psi);}
  else if (psi_no_map[no_ps] == 31){ppttx20_ax_psp_400_ub_ttxub(no_ps, psi);}
  else if (psi_no_map[no_ps] == 32){ppttx20_ax_psp_400_udx_ttxudx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 34){ppttx20_ax_psp_400_uux_ttxgg(no_ps, psi);}
  else if (psi_no_map[no_ps] == 35){ppttx20_ax_psp_400_uux_ttxddx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 36){ppttx20_ax_psp_400_uux_ttxuux(no_ps, psi);}
  else if (psi_no_map[no_ps] == 37){ppttx20_ax_psp_400_uux_ttxssx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 38){ppttx20_ax_psp_400_uux_ttxccx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 39){ppttx20_ax_psp_400_uux_ttxbbx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 40){ppttx20_ax_psp_400_usx_ttxusx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 42){ppttx20_ax_psp_400_ucx_ttxucx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 43){ppttx20_ax_psp_400_ubx_ttxubx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 44){ppttx20_ax_psp_400_bb_ttxbb(no_ps, psi);}
  else if (psi_no_map[no_ps] == 45){ppttx20_ax_psp_400_bdx_ttxbdx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 46){ppttx20_ax_psp_400_bux_ttxbux(no_ps, psi);}
  else if (psi_no_map[no_ps] == 47){ppttx20_ax_psp_400_bbx_ttxgg(no_ps, psi);}
  else if (psi_no_map[no_ps] == 48){ppttx20_ax_psp_400_bbx_ttxddx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 49){ppttx20_ax_psp_400_bbx_ttxuux(no_ps, psi);}
  else if (psi_no_map[no_ps] == 50){ppttx20_ax_psp_400_bbx_ttxbbx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 51){ppttx20_ax_psp_400_dxdx_ttxdxdx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 52){ppttx20_ax_psp_400_dxux_ttxdxux(no_ps, psi);}
  else if (psi_no_map[no_ps] == 53){ppttx20_ax_psp_400_dxsx_ttxdxsx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 54){ppttx20_ax_psp_400_dxcx_ttxdxcx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 56){ppttx20_ax_psp_400_dxbx_ttxdxbx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 57){ppttx20_ax_psp_400_uxux_ttxuxux(no_ps, psi);}
  else if (psi_no_map[no_ps] == 58){ppttx20_ax_psp_400_uxcx_ttxuxcx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 59){ppttx20_ax_psp_400_uxbx_ttxuxbx(no_ps, psi);}
  else if (psi_no_map[no_ps] == 60){ppttx20_ax_psp_400_bxbx_ttxbxbx(no_ps, psi);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
