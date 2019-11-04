#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ac_psp_200_gg_ttx(int no_ps, int channel, phasespace_set & psi){
  static Logger logger("ppttx20_ac_psp_200_gg_ttx");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  double g_IS_temp;
  int ran = 1;
  for (int i_p = 0; i_p < psi_c_p[no_ps][channel].size(); i_p++){
  }
  ran += psi_c_p[no_ps][channel].size();

  for (int i_f = 0; i_f < psi_c_f[no_ps][channel].size(); i_f++){
  }
  ran += psi_c_f[no_ps][channel].size();

  for (int i_t = 0; i_t < psi_c_t[no_ps][channel].size(); i_t++){
    if      (psi_c_t[no_ps][channel][i_t] ==   0){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 0] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 0]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 0] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 1]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  12,   1,   2,   4,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==   1){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 1] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 2]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 1] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 3]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  12,   1,   2,   8,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
  }
  ran += 2 * psi_c_t[no_ps][channel].size();

  for (int i_d = 0; i_d < psi_c_d[no_ps][channel].size(); i_d++){
    if      (psi_c_d[no_ps][channel][i_d] ==   0){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 0] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 0]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 0] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 1]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  12,   4,   8, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
  }
  ran += 2 * psi_c_d[no_ps][channel].size();

  /*
  if      (channel ==    0){ac_channel_11(r, 1, 2, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6, psi);}
  else if (channel ==    1){ac_channel_11(r, 1, 2, 8, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6, psi);}
  else if (channel ==    2){ac_channel_2v(r, 1, 2, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0, psi);}
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
