#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ac_psp_400_gg_ttxgg(int no_ps, int channel, phasespace_set & psi){
  static Logger logger("ppttx20_ac_psp_400_gg_ttxgg");
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
      psi_v_smin[no_ps][  1] = psi_smin_opt[no_ps][36];
      psi_v_smax[no_ps][  1] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][24], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 1] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 1]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  36, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   2){
      psi_v_smin[no_ps][  2] = psi_smin_opt[no_ps][20];
      psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][40], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 2] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 2]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  20, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   3){
      psi_v_smin[no_ps][  3] = psi_smin_opt[no_ps][48];
      psi_v_smax[no_ps][  3] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][12], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 3] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 3]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   4){
      psi_v_smin[no_ps][  4] = psi_smin_opt[no_ps][40];
      psi_v_smax[no_ps][  4] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][20], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 4] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 4]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  40, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   5){
      psi_v_smin[no_ps][  5] = psi_smin_opt[no_ps][24];
      psi_v_smax[no_ps][  5] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][36], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 5] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 5]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  24, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  5], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   6){
      psi_v_smin[no_ps][  6] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][24], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  6] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 6] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 6]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  6], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   7){
      psi_v_smin[no_ps][  7] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][40], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  7] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 7] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 7]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  7], psi_v_smax[no_ps][  7], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   8){
      psi_v_smin[no_ps][  8] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  7] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 8] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 8]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  8], psi_v_smax[no_ps][  7], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   9){
      psi_v_smin[no_ps][  9] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][52]);
      psi_v_smax[no_ps][  8] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 9] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 9]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  52, psi_v_smin[no_ps][  9], psi_v_smax[no_ps][  8], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  10){
      psi_v_smin[no_ps][ 10] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][20], 2), psi_smin_opt[no_ps][52]);
      psi_v_smax[no_ps][  8] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 10] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 10]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  52, psi_v_smin[no_ps][ 10], psi_v_smax[no_ps][  8], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  11){
      psi_v_smin[no_ps][ 11] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  6] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 11] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 11]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  28, psi_v_smin[no_ps][ 11], psi_v_smax[no_ps][  6], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  12){
      psi_v_smin[no_ps][ 12] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][36], 2), psi_smin_opt[no_ps][52]);
      psi_v_smax[no_ps][  8] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 12] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 12]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  52, psi_v_smin[no_ps][ 12], psi_v_smax[no_ps][  8], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  13){
      psi_v_smin[no_ps][ 13] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][20], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  6] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 13] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 13]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  28, psi_v_smin[no_ps][ 13], psi_v_smax[no_ps][  6], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  14){
      psi_v_smin[no_ps][ 14] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][36], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  7] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 14] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 14]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  44, psi_v_smin[no_ps][ 14], psi_v_smax[no_ps][  7], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  15){
      psi_v_smin[no_ps][ 15] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][56]);
      psi_v_smax[no_ps][  9] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 15] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 15]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 15], psi_v_smax[no_ps][  9], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  16){
      psi_v_smin[no_ps][ 16] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][24], 2), psi_smin_opt[no_ps][56]);
      psi_v_smax[no_ps][  9] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 16] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 16]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 16], psi_v_smax[no_ps][  9], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  17){
      psi_v_smin[no_ps][ 17] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][40], 2), psi_smin_opt[no_ps][56]);
      psi_v_smax[no_ps][  9] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 17] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 17]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 17], psi_v_smax[no_ps][  9], ran + i_p);
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
      psi_v_smin[no_ps][  1] = psi_smin_opt[no_ps][36];
      psi_v_smax[no_ps][ 11] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][24], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 19] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 19]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  36, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][ 11], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  20){
      psi_v_smin[no_ps][  2] = psi_smin_opt[no_ps][20];
      psi_v_smax[no_ps][ 12] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][40], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 20] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 20]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  20, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][ 12], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  21){
      psi_v_smin[no_ps][  3] = psi_smin_opt[no_ps][48];
      psi_v_smax[no_ps][ 13] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][12], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 21] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 21]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][ 13], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  22){
      psi_v_smin[no_ps][  4] = psi_smin_opt[no_ps][40];
      psi_v_smax[no_ps][ 14] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][20], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 22] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 22]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  40, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][ 14], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==  23){
      psi_v_smin[no_ps][  5] = psi_smin_opt[no_ps][24];
      psi_v_smax[no_ps][ 15] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][36], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 23] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 23]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   6,  24, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][ 15], ran + i_p);
    }
  }
  ran += psi_c_p[no_ps][channel].size();

  for (int i_f = 0; i_f < psi_c_f[no_ps][channel].size(); i_f++){
    if (psi_c_f[no_ps][channel][i_f] ==   0){
      psi_v_smin[no_ps][  5] = psi_smin_opt[no_ps][24];
      psi_v_smax[no_ps][  5] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][36], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 0] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 0]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  24, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  5], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   1){
      psi_v_smin[no_ps][ 16] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][24], 2), psi_smin_opt[no_ps][56]);
      psi_v_smax[no_ps][  9] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 1] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 1]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  56, psi_v_smin[no_ps][ 16], psi_v_smax[no_ps][  9], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   2){
      psi_v_smin[no_ps][  4] = psi_smin_opt[no_ps][40];
      psi_v_smax[no_ps][  4] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][20], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 2] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 2]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  40, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   3){
      psi_v_smin[no_ps][ 17] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][40], 2), psi_smin_opt[no_ps][56]);
      psi_v_smax[no_ps][  9] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 3] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 3]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  56, psi_v_smin[no_ps][ 17], psi_v_smax[no_ps][  9], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   4){
      psi_v_smin[no_ps][  0] = psi_smin_opt[no_ps][12];
      psi_v_smax[no_ps][  0] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][48], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 4] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 4]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   5){
      psi_v_smin[no_ps][ 11] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  6] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 5] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 5]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][ 11], psi_v_smax[no_ps][  6], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   6){
      psi_v_smin[no_ps][  3] = psi_smin_opt[no_ps][48];
      psi_v_smax[no_ps][  3] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][12], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 6] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 6]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  48, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   7){
      psi_v_smin[no_ps][ 15] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][56]);
      psi_v_smax[no_ps][  9] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 7] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 7]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  56, psi_v_smin[no_ps][ 15], psi_v_smax[no_ps][  9], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   8){
      psi_v_smin[no_ps][  2] = psi_smin_opt[no_ps][20];
      psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][40], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 8] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 8]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  20, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   9){
      psi_v_smin[no_ps][ 13] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][20], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  6] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 9] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 9]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][ 13], psi_v_smax[no_ps][  6], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==  10){
      psi_v_smin[no_ps][  8] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  7] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 10] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 10]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  44, psi_v_smin[no_ps][  8], psi_v_smax[no_ps][  7], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==  11){
      psi_v_smin[no_ps][  1] = psi_smin_opt[no_ps][36];
      psi_v_smax[no_ps][  1] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][24], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 11] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 11]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  36, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==  12){
      psi_v_smin[no_ps][ 14] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][36], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  7] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 12] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 12]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  44, psi_v_smin[no_ps][ 14], psi_v_smax[no_ps][  7], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==  13){
      psi_v_smin[no_ps][ 10] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][20], 2), psi_smin_opt[no_ps][52]);
      psi_v_smax[no_ps][  8] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 13] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 13]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  52, psi_v_smin[no_ps][ 10], psi_v_smax[no_ps][  8], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==  14){
      psi_v_smin[no_ps][ 12] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][36], 2), psi_smin_opt[no_ps][52]);
      psi_v_smax[no_ps][  8] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 14] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 14]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  52, psi_v_smin[no_ps][ 12], psi_v_smax[no_ps][  8], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==  15){
      psi_v_smin[no_ps][  9] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][52]);
      psi_v_smax[no_ps][  8] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 15] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 15]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  52, psi_v_smin[no_ps][  9], psi_v_smax[no_ps][  8], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==  16){
      psi_v_smin[no_ps][  6] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][24], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  6] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 16] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 16]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  6], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==  17){
      psi_v_smin[no_ps][  7] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][40], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  7] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 17] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 17]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  44, psi_v_smin[no_ps][  7], psi_v_smax[no_ps][  7], ran + i_f);
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
      psi.c_tchannel_opt(no_ps,   6,  60,   1,   2,   4,  56, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  56,   2,   5,  32,  24, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel(no_ps,   0,  24,   5,  34,   8,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  56,   2,   5,  16,  40, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel(no_ps,   0,  40,   5,  18,   8,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel(no_ps,   6,  40,   5,  18,  32,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  60,   1,   2,  32,  28, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  28,   2,  33,  16,  12, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel(no_ps,   6,  12,  33,  18,   4,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel(no_ps,   6,  24,   5,  34,  16,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   6,  56,   2,   5,   8,  48, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel(no_ps,   6,  48,   5,  10,  16,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  12){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 12] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 24]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 12] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 25]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  48,   5,  10,  32,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  13){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 13] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 26]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 13] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 27]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  28,   2,  33,   8,  20, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  14){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 14] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 28]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 14] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 29]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  20,  33,  10,   4,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  15){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 15] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 30]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 15] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 31]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  60,   1,   2,  16,  44, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  16){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 16] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 32]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 16] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 33]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  44,   2,  17,  32,  12, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  17){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 17] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 34]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 17] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 35]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  12,  17,  34,   4,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  18){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 18] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 36]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 18] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 37]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  44,   2,  17,   8,  36, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  19){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 19] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 38]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 19] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 39]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  36,  17,  10,   4,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  20){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 20] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 40]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 20] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 41]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   0,  36,  17,  10,  32,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  21){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 21] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 42]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 21] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 43]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   0,  20,  33,  10,  16,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  22){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 22] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 44]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 22] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 45]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  60,   1,   2,   8,  52, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  23){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 23] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 46]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 23] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 47]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  52,   2,   9,  32,  20, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  24){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 24] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 48]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 24] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 49]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   0,  20,   9,  34,   4,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  25){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 25] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 50]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 25] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 51]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  52,   2,   9,  16,  36, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  26){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 26] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 52]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 26] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 53]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   0,  36,   9,  18,   4,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  27){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 27] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 54]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 27] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 55]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  36,   9,  18,  32,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  28){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 28] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 56]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 28] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 57]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  12,  33,  18,   8,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  29){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 29] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 58]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 29] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 59]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  20,   9,  34,  16,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  30){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 30] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 60]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 30] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 61]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  52,   2,   9,   4,  48, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  31){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 31] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 62]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 31] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 63]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  48,   9,   6,  16,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  32){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 32] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 64]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 32] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 65]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  48,   9,   6,  32,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  33){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 33] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 66]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 33] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 67]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  28,   2,  33,   4,  24, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  34){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 34] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 68]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 34] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 69]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  24,  33,   6,   8,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  35){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 35] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 70]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 35] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 71]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  12,  17,  34,   8,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  36){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 36] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 72]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 36] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 73]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  44,   2,  17,   4,  40, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  37){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 37] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 74]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 37] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 75]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   6,  40,  17,   6,   8,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  38){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 38] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 76]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 38] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 77]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   0,  40,  17,   6,  32,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  39){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 39] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 78]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 39] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 79]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel(no_ps,   0,  24,  33,   6,  16,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  40){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 40] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 80]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 40] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 81]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  60,   2,   1,  32,  28, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  41){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 41] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 82]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 41] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 83]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  28,   1,  34,  12,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  42){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 42] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 84]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 42] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 85]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  60,   2,   1,  16,  44, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  43){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 43] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 86]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 43] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 87]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  44,   1,  18,  12,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  44){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 44] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 88]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 44] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 89]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  44,   1,  18,  36,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  45){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 45] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 90]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 45] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 91]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  28,   1,  34,  20,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  46){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 46] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 92]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 46] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 93]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  60,   2,   1,   8,  52, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  47){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 47] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 94]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 47] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 95]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  52,   1,  10,  20,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  48){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 48] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 96]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 48] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 97]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  52,   1,  10,  36,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  49){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 49] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 98]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 49] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 99]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  52,   1,  10,  48,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  50){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 50] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 100]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 50] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 101]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  44,   1,  18,  40,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  51){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 51] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 102]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 51] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 103]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  28,   1,  34,  24,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  52){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 52] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 104]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 52] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 105]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  60,   2,   1,   4,  56, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  53){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 53] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 106]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 53] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 107]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  56,   1,   6,  24,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  54){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 54] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 108]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 54] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 109]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  56,   1,   6,  40,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  55){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 55] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 110]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 55] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 111]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  56,   1,   6,  48,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  56){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 56] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 112]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 56] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 113]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  56,   2,   5,  48,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  57){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 57] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 114]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 57] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 115]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  56,   2,   5,  24,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  58){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 58] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 116]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 58] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 117]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  28,   2,  33,  24,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  59){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 59] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 118]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 59] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 119]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  56,   2,   5,  40,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  60){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 60] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 120]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 60] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 121]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  44,   2,  17,  40,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  61){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 61] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 122]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 61] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 123]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  44,   2,  17,  12,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  62){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 62] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 124]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 62] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 125]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  28,   2,  33,  12,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  63){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 63] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 126]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 63] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 127]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  52,   2,   9,  48,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  64){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 64] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 128]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 64] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 129]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  52,   2,   9,  20,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  65){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 65] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 130]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 65] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 131]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  28,   2,  33,  20,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  66){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 66] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 132]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 66] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 133]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  52,   2,   9,  36,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  67){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 67] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 134]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 67] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 135]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  44,   2,  17,  36,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  68){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 68] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 136]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 68] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 137]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  60,   1,   2,  12,  48, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  69){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 69] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 138]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 69] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 139]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  60,   1,   2,  36,  24, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  70){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 70] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 140]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 70] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 141]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  60,   1,   2,  20,  40, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  71){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 71] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 142]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 71] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 143]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   0,  60,   1,   2,  48,  12, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  72){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 72] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 144]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 72] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 145]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  60,   1,   2,  40,  20, ran + 2 * i_t, ran + 2 * i_t + 1);
    }
    else if (psi_c_t[no_ps][channel][i_t] ==  73){
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 73] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 146]->get_random(psi_r[ran + 2 * i_t], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][2] + 73] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][2] + 147]->get_random(psi_r[ran + 2 * i_t + 1], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_tchannel_opt(no_ps,   6,  60,   1,   2,  24,  36, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_decay(no_ps,  36,   4,  32, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  20,   4,  16, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  48,  16,  32, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  24,   8,  16, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  28,   4,  24, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  44,   4,  40, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  44,  32,  12, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  52,   4,  48, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  52,  32,  20, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  28,  16,  12, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  52,  16,  36, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  28,   8,  20, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  44,   8,  36, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  56,   8,  48, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  56,  32,  24, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  56,  16,  40, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  60,   4,  56, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  60,  32,  28, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  60,  16,  44, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  60,   8,  52, ran + 2 * i_d, ran + 2 * i_d + 1);
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
  if      (channel ==    0){ac_channel_1111(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==    1){ac_channel_1111(r, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==    2){ac_channel_1111(r, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==    3){ac_channel_1111(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  else if (channel ==    4){ac_channel_1111(r, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==    5){ac_channel_1111(r, 1, 2, 4, 16, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==    6){ac_channel_1111(r, 1, 2, 4, 32, 16, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==    7){ac_channel_1111(r, 1, 2, 32, 4, 16, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  else if (channel ==    8){ac_channel_1111(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  else if (channel ==    9){ac_channel_1111(r, 1, 2, 16, 4, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  else if (channel ==   10){ac_channel_1111(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==   11){ac_channel_1111(r, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==   12){ac_channel_1111(r, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   13){ac_channel_1111(r, 1, 2, 8, 4, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   14){ac_channel_1111(r, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   15){ac_channel_1111(r, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  else if (channel ==   16){ac_channel_1111(r, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   17){ac_channel_1111(r, 1, 2, 8, 16, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   18){ac_channel_1111(r, 1, 2, 8, 32, 16, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   19){ac_channel_1111(r, 1, 2, 32, 8, 16, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  else if (channel ==   20){ac_channel_1111(r, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  else if (channel ==   21){ac_channel_1111(r, 1, 2, 16, 8, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  else if (channel ==   22){ac_channel_1111(r, 1, 2, 16, 32, 8, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==   23){ac_channel_1111(r, 1, 2, 32, 16, 8, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==   24){ac_channel_2v11(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   25){ac_channel_2v11(r, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   26){ac_channel_2v11(r, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   27){ac_channel_2v11(r, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   28){ac_channel_2v11(r, 1, 2, 4, 16, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   29){ac_channel_2v11(r, 1, 2, 4, 32, 16, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   30){ac_channel_2v11(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==   31){ac_channel_2v11(r, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   32){ac_channel_2v11(r, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   33){ac_channel_2v11(r, 1, 2, 8, 16, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   34){ac_channel_2v11(r, 1, 2, 8, 32, 16, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   35){ac_channel_2v11(r, 1, 2, 16, 32, 8, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==   36){ac_channel_3yv1(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   37){ac_channel_3yv1(r, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   38){ac_channel_3yv1(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   39){ac_channel_3yv1(r, 1, 2, 4, 16, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  else if (channel ==   40){ac_channel_3yv1(r, 1, 2, 32, 4, 16, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   41){ac_channel_3yv1(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   42){ac_channel_3yv1(r, 1, 2, 16, 4, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   43){ac_channel_3yv1(r, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   44){ac_channel_3yv1(r, 1, 2, 8, 4, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   45){ac_channel_3yv1(r, 1, 2, 8, 16, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  else if (channel ==   46){ac_channel_3yv1(r, 1, 2, 32, 8, 16, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   47){ac_channel_3yv1(r, 1, 2, 16, 8, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   48){ac_channel_4yyv(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==   49){ac_channel_4yyv(r, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   50){ac_channel_4yyv(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   51){ac_channel_4yyv(r, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   52){ac_channel_4yyv(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   53){ac_channel_4yyv(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   54){ac_channel_4yyv(r, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   55){ac_channel_4yyv(r, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==   56){ac_channel_4yyv(r, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   57){ac_channel_4yyv(r, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   58){ac_channel_4yyv(r, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   59){ac_channel_4yyv(r, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   60){ac_channel_42vv(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   61){ac_channel_42vv(r, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   62){ac_channel_42vv(r, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   63){ac_channel_12v1(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   64){ac_channel_12v1(r, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   65){ac_channel_12v1(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   66){ac_channel_12v1(r, 1, 2, 4, 16, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  else if (channel ==   67){ac_channel_12v1(r, 1, 2, 32, 4, 16, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   68){ac_channel_12v1(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   69){ac_channel_12v1(r, 1, 2, 16, 4, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   70){ac_channel_12v1(r, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   71){ac_channel_12v1(r, 1, 2, 8, 4, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  else if (channel ==   72){ac_channel_12v1(r, 1, 2, 8, 16, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  else if (channel ==   73){ac_channel_12v1(r, 1, 2, 32, 8, 16, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   74){ac_channel_12v1(r, 1, 2, 16, 8, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   75){ac_channel_13yv(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  else if (channel ==   76){ac_channel_13yv(r, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   77){ac_channel_13yv(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   78){ac_channel_13yv(r, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   79){ac_channel_13yv(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   80){ac_channel_13yv(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   81){ac_channel_13yv(r, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   82){ac_channel_13yv(r, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  else if (channel ==   83){ac_channel_13yv(r, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   84){ac_channel_13yv(r, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   85){ac_channel_13yv(r, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   86){ac_channel_13yv(r, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  else if (channel ==   87){ac_channel_112v(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==   88){ac_channel_112v(r, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   89){ac_channel_112v(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   90){ac_channel_112v(r, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   91){ac_channel_112v(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   92){ac_channel_112v(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   93){ac_channel_112v(r, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   94){ac_channel_112v(r, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  else if (channel ==   95){ac_channel_112v(r, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   96){ac_channel_112v(r, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   97){ac_channel_112v(r, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==   98){ac_channel_112v(r, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  else if (channel ==   99){ac_channel_2v2v(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==  100){ac_channel_2v2v(r, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==  101){ac_channel_2v2v(r, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==  102){ac_channel_2v2v(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==  103){ac_channel_2v2v(r, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  else if (channel ==  104){ac_channel_2v2v(r, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
