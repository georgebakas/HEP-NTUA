#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ag_psp_born(int no_ps, int zero, phasespace_set & psi){
  static Logger logger("ppttx20_ag_psp_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (psi_no_map[no_ps] == 0){}
  else if (psi_no_map[no_ps] == 1){ppttx20_ag_psp_200_gg_ttx(no_ps, zero, psi);}
  else if (psi_no_map[no_ps] == 2){ppttx20_ag_psp_200_ddx_ttx(no_ps, zero, psi);}
  else if (psi_no_map[no_ps] == 3){ppttx20_ag_psp_200_uux_ttx(no_ps, zero, psi);}
  else if (psi_no_map[no_ps] == 4){ppttx20_ag_psp_200_bbx_ttx(no_ps, zero, psi);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
