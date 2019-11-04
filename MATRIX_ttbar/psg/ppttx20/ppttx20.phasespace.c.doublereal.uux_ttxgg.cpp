#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ac_psp_400_uux_ttxgg(int no_ps, int channel, phasespace_set & psi){
  static Logger logger("ppttx20_ac_psp_400_uux_ttxgg");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  double g_IS_temp;
  int ran = 1;
  for (int i_p = 0; i_p < psi_c_p[no_ps][channel].size(); i_p++){
    if (psi_c_p[no_ps][channel][i_p] ==   0){
      psi_v_smin[no_ps][  0] = psi_smin_opt[no_ps][12];
      psi_v_smax[no_ps][  0] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][48], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 0] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 0]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   1){
      psi_v_smin[no_ps][  1] = psi_smin_opt[no_ps][24];
      psi_v_smax[no_ps][  1] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][36], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 1] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 1]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  24, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   2){
      psi_v_smin[no_ps][  2] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][24], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 2] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 2]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   3){
      psi_v_smin[no_ps][  3] = psi_smin_opt[no_ps][40];
      psi_v_smax[no_ps][  3] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][20], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 3] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 3]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  40, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   4){
      psi_v_smin[no_ps][  4] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][40], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  4] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 4] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 4]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   5){
      psi_v_smin[no_ps][  5] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  4] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 5] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 5]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  4], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   6){
      psi_v_smin[no_ps][  6] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 6] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 6]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  2], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   7){
      psi_v_smin[no_ps][  7] = psi_smin_opt[no_ps][20];
      psi_v_smax[no_ps][  5] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][40], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 7] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 7]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  20, psi_v_smin[no_ps][  7], psi_v_smax[no_ps][  5], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   8){
      psi_v_smin[no_ps][  8] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][20], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 8] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 8]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  8], psi_v_smax[no_ps][  2], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   9){
      psi_v_smin[no_ps][  9] = psi_smin_opt[no_ps][36];
      psi_v_smax[no_ps][  6] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][24], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 9] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 9]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  36, psi_v_smin[no_ps][  9], psi_v_smax[no_ps][  6], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  10){
      psi_v_smin[no_ps][ 10] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][36], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  4] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 10] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 10]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  44, psi_v_smin[no_ps][ 10], psi_v_smax[no_ps][  4], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  11){
      psi_v_smin[no_ps][ 11] = psi_smin_opt[no_ps][48];
      psi_v_smax[no_ps][  7] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][12], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 11] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 11]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  48, psi_v_smin[no_ps][ 11], psi_v_smax[no_ps][  7], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  12){
      psi_v_smin[no_ps][ 12] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][56]);
      psi_v_smax[no_ps][  8] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 12] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 12]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 12], psi_v_smax[no_ps][  8], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  13){
      psi_v_smin[no_ps][ 13] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][24], 2), psi_smin_opt[no_ps][56]);
      psi_v_smax[no_ps][  8] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 13] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 13]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 13], psi_v_smax[no_ps][  8], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  14){
      psi_v_smin[no_ps][ 14] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][40], 2), psi_smin_opt[no_ps][56]);
      psi_v_smax[no_ps][  8] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 14] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 14]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 14], psi_v_smax[no_ps][  8], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  15){
      psi_v_smin[no_ps][ 15] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][52]);
      psi_v_smax[no_ps][  9] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 15] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 15]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  52, psi_v_smin[no_ps][ 15], psi_v_smax[no_ps][  9], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  16){
      psi_v_smin[no_ps][ 16] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][20], 2), psi_smin_opt[no_ps][52]);
      psi_v_smax[no_ps][  9] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 16] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 16]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  52, psi_v_smin[no_ps][ 16], psi_v_smax[no_ps][  9], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  17){
      psi_v_smin[no_ps][ 17] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][36], 2), psi_smin_opt[no_ps][52]);
      psi_v_smax[no_ps][  9] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 17] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 17]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  52, psi_v_smin[no_ps][ 17], psi_v_smax[no_ps][  9], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  18){
      psi_v_smin[no_ps][  0] = psi_smin_opt[no_ps][12];
      psi_v_smax[no_ps][ 10] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][48], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 18] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 18]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][ 10], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  19){
      psi_v_smin[no_ps][  9] = psi_smin_opt[no_ps][36];
      psi_v_smax[no_ps][ 11] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][24], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 19] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 19]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  36, psi_v_smin[no_ps][  9], psi_v_smax[no_ps][ 11], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  20){
      psi_v_smin[no_ps][  7] = psi_smin_opt[no_ps][20];
      psi_v_smax[no_ps][ 12] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][40], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 20] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 20]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  20, psi_v_smin[no_ps][  7], psi_v_smax[no_ps][ 12], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  21){
      psi_v_smin[no_ps][ 11] = psi_smin_opt[no_ps][48];
      psi_v_smax[no_ps][ 13] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][12], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 21] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 21]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  48, psi_v_smin[no_ps][ 11], psi_v_smax[no_ps][ 13], ran + i_p);
    }
  }
  ran += psi_c_p[no_ps][channel].size();

  for (int i_f = 0; i_f < psi_c_f[no_ps][channel].size(); i_f++){
    if (psi_c_f[no_ps][channel][i_f] ==   0){
      psi_v_smin[no_ps][  6] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 0] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 0]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  2], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   1){
      psi_v_smin[no_ps][  5] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  4] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 1] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 1]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  44, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  4], ran + i_f);
    }
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
      psi.c_tchannel_opt(no_ps,   0,  60,   2,   1,  32,  28, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  28,   1,  34,  12,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==   2){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 2] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 4]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 2] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 5]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  60,   2,   1,  16,  44, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==   3){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 3] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 6]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 3] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 7]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  44,   1,  18,  12,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==   4){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 4] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 8]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 4] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 9]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  60,   1,   2,  32,  28, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==   5){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 5] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 10]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 5] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 11]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  28,   2,  33,  16,  12, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==   6){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 6] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 12]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 6] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 13]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  60,   1,   2,  16,  44, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==   7){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 7] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 14]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 7] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 15]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  44,   2,  17,  32,  12, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==   8){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 8] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 16]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 8] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 17]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  44,   2,  17,  12,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==   9){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 9] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 18]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 9] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 19]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  28,   2,  33,  12,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  10){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 10] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 20]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 10] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 21]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  60,   1,   2,  12,  48, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  11){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 11] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 22]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 11] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 23]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  60,   1,   2,  48,  12, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_decay(no_ps,  44,   4,  40, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  40,   8,  32, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  44,  32,  12, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==   6){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 6] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 12]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 6] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 13]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  28,  16,  12, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==   7){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 7] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 14]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 7] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 15]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  28,   8,  20, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==   8){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 8] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 16]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 8] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 17]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  20,   4,  16, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==   9){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 9] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 18]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 9] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 19]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  44,   8,  36, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  10){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 10] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 20]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 10] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 21]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  36,   4,  32, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  11){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 11] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 22]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 11] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 23]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  60,   4,  56, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  12){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 12] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 24]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 12] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 25]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  56,   8,  48, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  13){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 13] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 26]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 13] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 27]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  48,  16,  32, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  14){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 14] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 28]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 14] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 29]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  56,  32,  24, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  15){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 15] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 30]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 15] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 31]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  60,  32,  28, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  16){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 16] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 32]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 16] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 33]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  56,  16,  40, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  17){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 17] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 34]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 17] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 35]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  60,  16,  44, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  18){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 18] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 36]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 18] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 37]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  60,   8,  52, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  19){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 19] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 38]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 19] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 39]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  52,   4,  48, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  20){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 20] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 40]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 20] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 41]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  52,  32,  20, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  21){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 21] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 42]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 21] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 43]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  52,  16,  36, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  22){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 22] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 44]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 22] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 45]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  60,  12,  48, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  23){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 23] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 46]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 23] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 47]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  60,  36,  24, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
    else if (psi_c_d[no_ps][channel][i_d] ==  24){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 24] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 48]->get_random(psi_r[ran + 2 * i_d], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][3] + 24] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][3] + 49]->get_random(psi_r[ran + 2 * i_d + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_decay(no_ps,  60,  20,  40, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
  }
  ran += 2 * psi_c_d[no_ps][channel].size();

  /*
  if      (channel ==    0){ac_channel_2v11(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==    1){ac_channel_2v11(r, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==    2){ac_channel_3yv1(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==    3){ac_channel_3yv1(r, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==    4){ac_channel_3yv1(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==    5){ac_channel_3yv1(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==    6){ac_channel_3yv1(r, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==    7){ac_channel_3yv1(r, 1, 2, 8, 4, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==    8){ac_channel_4yyv(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==    9){ac_channel_4yyv(r, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   10){ac_channel_4yyv(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   11){ac_channel_4yyv(r, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   12){ac_channel_4yyv(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   13){ac_channel_4yyv(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   14){ac_channel_4yyv(r, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   15){ac_channel_4yyv(r, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==   16){ac_channel_4yyv(r, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   17){ac_channel_4yyv(r, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   18){ac_channel_4yyv(r, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   19){ac_channel_4yyv(r, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   20){ac_channel_42vv(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   21){ac_channel_42vv(r, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   22){ac_channel_42vv(r, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   23){ac_channel_12v1(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   24){ac_channel_12v1(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   25){ac_channel_13yv(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   26){ac_channel_13yv(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   27){ac_channel_13yv(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   28){ac_channel_13yv(r, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   29){ac_channel_13yv(r, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   30){ac_channel_13yv(r, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   31){ac_channel_112v(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   32){ac_channel_112v(r, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   33){ac_channel_2v2v(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   34){ac_channel_2v2v(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
