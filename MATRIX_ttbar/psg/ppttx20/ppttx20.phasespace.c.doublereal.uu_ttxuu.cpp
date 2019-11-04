#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ac_psp_400_uu_ttxuu(int no_ps, int channel, phasespace_set & psi){
  static Logger logger("ppttx20_ac_psp_400_uu_ttxuu");
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
      psi_v_smin[no_ps][  1] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  1] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 1] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 1]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], ran + i_p);
    }
    else if (psi_c_p[no_ps][channel][i_p] ==   2){
      psi_v_smin[no_ps][  2] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][0] + 2] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][0] + 2]->get_random(psi_r[ran + i_p], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], ran + i_p);
    }
  }
  ran += psi_c_p[no_ps][channel].size();

  for (int i_f = 0; i_f < psi_c_f[no_ps][channel].size(); i_f++){
    if (psi_c_f[no_ps][channel][i_f] ==   0){
      psi_v_smin[no_ps][  0] = psi_smin_opt[no_ps][12];
      psi_v_smax[no_ps][  0] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][48], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 0] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 0]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   1){
      psi_v_smin[no_ps][  2] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][28]);
      psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 1] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 1]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], ran + i_f);
    }
    else if (psi_c_f[no_ps][channel][i_f] ==   2){
      psi_v_smin[no_ps][  1] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][44]);
      psi_v_smax[no_ps][  1] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2);
      if (psi_container_IS_switch[psi_container_IS_startvalue[no_ps][1] + 2] == 1){
        psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][1] + 2]->get_random(psi_r[ran + i_f], g_IS_temp);
        psi_MC_g_IS_global *= g_IS_temp;
      }
      psi.c_timelikeinvariant(no_ps,  44, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], ran + i_f);
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
      psi.c_tchannel_opt(no_ps,   0,  60,   1,   2,  32,  28, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  28,   2,  33,  16,  12, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel(no_ps,   6,  12,  33,  18,   4,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  60,   1,   2,  16,  44, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  44,   2,  17,  32,  12, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel(no_ps,   6,  12,  17,  34,   4,   8, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel(no_ps,   6,  12,  33,  18,   8,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel(no_ps,   6,  12,  17,  34,   8,   4, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  60,   2,   1,  32,  28, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  28,   1,  34,  12,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  60,   2,   1,  16,  44, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  44,   1,  18,  12,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  44,   2,  17,  12,  32, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_tchannel_opt(no_ps,   0,  28,   2,  33,  12,  16, ran + 2 * i_t, ran + 2 * i_t + 1);
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
      psi.c_decay(no_ps,  44,  32,  12, ran + 2 * i_d, ran + 2 * i_d + 1);
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
      psi.c_decay(no_ps,  28,  16,  12, ran + 2 * i_d, ran + 2 * i_d + 1);
    }
  }
  ran += 2 * psi_c_d[no_ps][channel].size();

  /*
  if      (channel ==    0){ac_channel_1111(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  else if (channel ==    1){ac_channel_1111(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  else if (channel ==    2){ac_channel_1111(r, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  else if (channel ==    3){ac_channel_1111(r, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  else if (channel ==    4){ac_channel_2v11(r, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==    5){ac_channel_2v11(r, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==    6){ac_channel_3yv1(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==    7){ac_channel_3yv1(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==    8){ac_channel_12v1(r, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==    9){ac_channel_12v1(r, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   10){ac_channel_13yv(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   11){ac_channel_13yv(r, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   12){ac_channel_112v(r, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  else if (channel ==   13){ac_channel_112v(r, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
