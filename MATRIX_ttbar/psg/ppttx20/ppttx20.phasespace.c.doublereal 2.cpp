#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ac_psp_doublereal(int no_ps, int channel, phasespace_set & psi){
  static Logger logger("ppttx20_ac_psp_doublereal");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (psi_no_map[no_ps] == 0){}
  else if (psi_no_map[no_ps] ==  1){ppttx20_ac_psp_400_gg_ttxgg(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  2){ppttx20_ac_psp_400_gg_ttxddx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  3){ppttx20_ac_psp_400_gg_ttxuux(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  4){ppttx20_ac_psp_400_gg_ttxbbx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  5){ppttx20_ac_psp_400_gd_ttxgd(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  6){ppttx20_ac_psp_400_gu_ttxgu(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  7){ppttx20_ac_psp_400_gb_ttxgb(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  8){ppttx20_ac_psp_400_gdx_ttxgdx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  9){ppttx20_ac_psp_400_gux_ttxgux(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 10){ppttx20_ac_psp_400_gbx_ttxgbx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 11){ppttx20_ac_psp_400_dd_ttxdd(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 12){ppttx20_ac_psp_400_du_ttxdu(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 13){ppttx20_ac_psp_400_ds_ttxds(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 14){ppttx20_ac_psp_400_dc_ttxdc(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 16){ppttx20_ac_psp_400_db_ttxdb(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 17){ppttx20_ac_psp_400_ddx_ttxgg(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 18){ppttx20_ac_psp_400_ddx_ttxddx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 19){ppttx20_ac_psp_400_ddx_ttxuux(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 20){ppttx20_ac_psp_400_ddx_ttxssx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 21){ppttx20_ac_psp_400_ddx_ttxccx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 22){ppttx20_ac_psp_400_ddx_ttxbbx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 23){ppttx20_ac_psp_400_dux_ttxdux(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 25){ppttx20_ac_psp_400_dsx_ttxdsx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 27){ppttx20_ac_psp_400_dcx_ttxdcx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 28){ppttx20_ac_psp_400_dbx_ttxdbx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 29){ppttx20_ac_psp_400_uu_ttxuu(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 30){ppttx20_ac_psp_400_uc_ttxuc(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 31){ppttx20_ac_psp_400_ub_ttxub(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 32){ppttx20_ac_psp_400_udx_ttxudx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 34){ppttx20_ac_psp_400_uux_ttxgg(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 35){ppttx20_ac_psp_400_uux_ttxddx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 36){ppttx20_ac_psp_400_uux_ttxuux(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 37){ppttx20_ac_psp_400_uux_ttxssx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 38){ppttx20_ac_psp_400_uux_ttxccx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 39){ppttx20_ac_psp_400_uux_ttxbbx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 40){ppttx20_ac_psp_400_usx_ttxusx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 42){ppttx20_ac_psp_400_ucx_ttxucx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 43){ppttx20_ac_psp_400_ubx_ttxubx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 44){ppttx20_ac_psp_400_bb_ttxbb(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 45){ppttx20_ac_psp_400_bdx_ttxbdx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 46){ppttx20_ac_psp_400_bux_ttxbux(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 47){ppttx20_ac_psp_400_bbx_ttxgg(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 48){ppttx20_ac_psp_400_bbx_ttxddx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 49){ppttx20_ac_psp_400_bbx_ttxuux(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 50){ppttx20_ac_psp_400_bbx_ttxbbx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 51){ppttx20_ac_psp_400_dxdx_ttxdxdx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 52){ppttx20_ac_psp_400_dxux_ttxdxux(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 53){ppttx20_ac_psp_400_dxsx_ttxdxsx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 54){ppttx20_ac_psp_400_dxcx_ttxdxcx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 56){ppttx20_ac_psp_400_dxbx_ttxdxbx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 57){ppttx20_ac_psp_400_uxux_ttxuxux(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 58){ppttx20_ac_psp_400_uxcx_ttxuxcx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 59){ppttx20_ac_psp_400_uxbx_ttxuxbx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 60){ppttx20_ac_psp_400_bxbx_ttxbxbx(no_ps, channel, psi);} 

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
