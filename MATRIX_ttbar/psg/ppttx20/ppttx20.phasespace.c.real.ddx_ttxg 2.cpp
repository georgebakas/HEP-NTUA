#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ac_psp_300_ddx_ttxg(int no_ps, int channel, phasespace_set & psi){
  static Logger logger("ppttx20_ac_psp_300_ddx_ttxg");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  double g_IS_temp;
  int ran = 1;
  for (int i_p = 0; i_p < psi_c_p[no_ps][channel].size(); i_p++){
    if (psi_c_p[no_ps][channel][i_p] ==   0){
      psi_v_smin[no_ps][  0] = psi_smin_opt[no_ps][12];
      psi_v_smax[no_ps][  0] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 0] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 0]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   1){
      psi_v_smin[no_ps][  1] = psi_smin_opt[no_ps][24];
      psi_v_smax[no_ps][  1] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 1] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 1]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  24, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   2){
      psi_v_smin[no_ps][  2] = psi_smin_opt[no_ps][20];
      psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 2] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 2]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  20, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], ran + i_p);
    }
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
      psi.c_tchannel_opt(no_ps,   0,  28,   2,   1,  16,  12, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  28,   1,   2,  16,  12, ran + 2 * i_t, ran + 2 * i_t + 1);
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
    else if (psi_c_d[no_ps][channel][i_d] ==   1){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 1] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 2]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 1] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 3]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  28,   4,  24, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==   2){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 2] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 4]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 2] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 5]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  24,   8,  16, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==   3){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 3] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 6]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 3] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 7]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  28,  16,  12, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==   4){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 4] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 8]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 4] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 9]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  28,   8,  20, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==   5){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 5] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 10]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 5] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 11]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  20,   4,  16, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
  }
  ran += 2 * psi_c_d[no_ps][channel].size();

  /*
  if      (channel ==    0){ac_channel_2v1(r, 1, 2, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0, psi);}
  else if (channel ==    1){ac_channel_3yv(r, 1, 2, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0, psi);}
  else if (channel ==    2){ac_channel_3yv(r, 1, 2, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0, psi);}
  else if (channel ==    3){ac_channel_3yv(r, 1, 2, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0, psi);}
  else if (channel ==    4){ac_channel_12v(r, 1, 2, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0, psi);}
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
