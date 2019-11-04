#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ac_psp_real(int no_ps, int channel, phasespace_set & psi){
  static Logger logger("ppttx20_ac_psp_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (psi_no_map[no_ps] == 0){}
  else if (psi_no_map[no_ps] ==  1){ppttx20_ac_psp_300_gg_ttxg(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  2){ppttx20_ac_psp_300_gd_ttxd(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  3){ppttx20_ac_psp_300_gu_ttxu(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  4){ppttx20_ac_psp_300_gb_ttxb(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  5){ppttx20_ac_psp_300_gdx_ttxdx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  6){ppttx20_ac_psp_300_gux_ttxux(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  7){ppttx20_ac_psp_300_gbx_ttxbx(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  8){ppttx20_ac_psp_300_ddx_ttxg(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] ==  9){ppttx20_ac_psp_300_uux_ttxg(no_ps, channel, psi);} 
  else if (psi_no_map[no_ps] == 10){ppttx20_ac_psp_300_bbx_ttxg(no_ps, channel, psi);} 

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
