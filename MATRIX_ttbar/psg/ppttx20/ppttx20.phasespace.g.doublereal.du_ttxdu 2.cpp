#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ag_psp_400_du_ttxdu(int no_ps, int zero, phasespace_set & psi){
  static Logger logger("ppttx20_ag_psp_400_du_ttxdu");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (psi_xbp[no_ps][12] == nullvector){psi_xbp[no_ps][12] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8];}
  if (psi_xbs[no_ps][12] == 0.){psi_xbs[no_ps][12] = psi_xbp[no_ps][12].m2();}
  if (psi_xbs[no_ps][12] < 0.){psi_xbs[no_ps][12] = 0.;}
  if (psi_xbsqrts[no_ps][12] == 0.){psi_xbsqrts[no_ps][12] = sqrt(psi_xbs[no_ps][12]);}
  if (psi_xbp[no_ps][17] == nullvector){psi_xbp[no_ps][17] = psi_xbp[no_ps][1] - psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][17] == 0.){psi_xbs[no_ps][17] = psi_xbp[no_ps][17].m2();}
  if (psi_xbs[no_ps][17] > 0.){psi_xbs[no_ps][17] = 0.;}
  if (psi_xbp[no_ps][34] == nullvector){psi_xbp[no_ps][34] = psi_xbp[no_ps][2] - psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][34] == 0.){psi_xbs[no_ps][34] = psi_xbp[no_ps][34].m2();}
  if (psi_xbs[no_ps][34] > 0.){psi_xbs[no_ps][34] = 0.;}
  if (psi_xbp[no_ps][13] == nullvector){psi_xbp[no_ps][13] = psi_xbp[no_ps][1] - psi_xbp[no_ps][4] - psi_xbp[no_ps][8];}
  if (psi_xbs[no_ps][13] == 0.){psi_xbs[no_ps][13] = psi_xbp[no_ps][13].m2();}
  if (psi_xbs[no_ps][13] > 0.){psi_xbs[no_ps][13] = 0.;}
  if (psi_xbp[no_ps][14] == nullvector){psi_xbp[no_ps][14] = psi_xbp[no_ps][2] - psi_xbp[no_ps][4] - psi_xbp[no_ps][8];}
  if (psi_xbs[no_ps][14] == 0.){psi_xbs[no_ps][14] = psi_xbp[no_ps][14].m2();}
  if (psi_xbs[no_ps][14] > 0.){psi_xbs[no_ps][14] = 0.;}
  if (psi_xbp[no_ps][21] == nullvector){psi_xbp[no_ps][21] = psi_xbp[no_ps][1] - psi_xbp[no_ps][4] - psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][21] == 0.){psi_xbs[no_ps][21] = psi_xbp[no_ps][21].m2();}
  if (psi_xbs[no_ps][21] > 0.){psi_xbs[no_ps][21] = 0.;}
  if (psi_xbp[no_ps][25] == nullvector){psi_xbp[no_ps][25] = psi_xbp[no_ps][1] - psi_xbp[no_ps][8] - psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][25] == 0.){psi_xbs[no_ps][25] = psi_xbp[no_ps][25].m2();}
  if (psi_xbs[no_ps][25] > 0.){psi_xbs[no_ps][25] = 0.;}
  if (psi_xbp[no_ps][28] == nullvector){psi_xbp[no_ps][28] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8] + psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][28] == 0.){psi_xbs[no_ps][28] = psi_xbp[no_ps][28].m2();}
  if (psi_xbs[no_ps][28] < 0.){psi_xbs[no_ps][28] = 0.;}
  if (psi_xbp[no_ps][44] == nullvector){psi_xbp[no_ps][44] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][44] == 0.){psi_xbs[no_ps][44] = psi_xbp[no_ps][44].m2();}
  if (psi_xbs[no_ps][44] < 0.){psi_xbs[no_ps][44] = 0.;}
  if (psi_xbp[no_ps][60] == nullvector){psi_xbp[no_ps][60] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8] + psi_xbp[no_ps][16] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][60] == 0.){psi_xbs[no_ps][60] = psi_xbp[no_ps][60].m2();}
  if (psi_xbs[no_ps][60] < 0.){psi_xbs[no_ps][60] = 0.;}

  psi_v_smin[no_ps][  0] = psi_smin_opt[no_ps][12]; // 7 times
  psi_v_smin[no_ps][  1] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][28]); // 2 times
  psi_v_smin[no_ps][  2] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][44]); // 5 times

  psi_v_smax[no_ps][  0] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][48], 2); // 7 times
  psi_v_smax[no_ps][  1] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2); // 2 times
  psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2); // 5 times

  //  vector<double> g_p(3);
  {psi_g_p[no_ps][  0] = psi.g_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0]); /* npp 2 */} // 5 times
  {psi_g_p[no_ps][  1] = psi.g_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1]); /* npp 3 */} // 1 times
  {psi_g_p[no_ps][  2] = psi.g_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2]); /* npp 3 */} // 1 times

  //  vector<double> g_f(3);
  {psi_g_f[no_ps][  0] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0]); /* nfp 2 */} // 2 times
  {psi_g_f[no_ps][  1] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2]); /* nfp 3 */} // 4 times
  {psi_g_f[no_ps][  2] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1]); /* nfp 3 */} // 1 times

  //  vector<double> g_t(7);
  {psi_g_t[no_ps][  0] = psi.g_tchannel_opt(no_ps,   0,  60,   1,   2,  16,  44); /* ntp 0 */} // 5 times
  {psi_g_t[no_ps][  1] = psi.g_tchannel_opt(no_ps,   0,  44,   2,  17,  32,  12); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][  2] = psi.g_tchannel(no_ps,   6,  12,  17,  34,   4,   8); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][  3] = psi.g_tchannel(no_ps,   6,  12,  17,  34,   8,   4); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][  4] = psi.g_tchannel_opt(no_ps,   0,  60,   2,   1,  32,  28); /* ntp 0 */} // 2 times
  {psi_g_t[no_ps][  5] = psi.g_tchannel_opt(no_ps,   0,  28,   1,  34,  12,  16); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][  6] = psi.g_tchannel_opt(no_ps,   0,  44,   2,  17,  12,  32); /* ntp 1 */} // 1 times

  //  vector<double> g_d(3);
  {psi_g_d[no_ps][  0] = psi.g_decay(no_ps,  12,   4,   8); /* ndp 2 */} // 5 times
  {psi_g_d[no_ps][  1] = psi.g_decay(no_ps,  28,  16,  12); /* ndp 3 */} // 1 times
  {psi_g_d[no_ps][  2] = psi.g_decay(no_ps,  44,  32,  12); /* ndp 3 */} // 1 times


  if (psi_weight_IS == 2 || psi_weight_IS == 4){
    double g_IS_temp;
    vector<vector<double> > inv_r(4);
    inv_r[0].resize(3);
    inv_r[1].resize(3);
    inv_r[2].resize(14);
    inv_r[3].resize(6);

    {psi.inv_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], inv_r[0][0]);} // 5 times
    {psi.inv_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], inv_r[0][1]);} // 1 times
    {psi.inv_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], inv_r[0][2]);} // 1 times

    {psi.inv_timelikeinvariant(no_ps,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], inv_r[1][0]);} // 2 times
    {psi.inv_timelikeinvariant(no_ps,  44, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], inv_r[1][1]);} // 4 times
    {psi.inv_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], inv_r[1][2]);} // 1 times

    {psi.inv_tchannel_opt(no_ps,   0,  60,   1,   2,  16,  44, inv_r[2][0], inv_r[2][1]);} // 5 times
    {psi.inv_tchannel_opt(no_ps,   0,  44,   2,  17,  32,  12, inv_r[2][2], inv_r[2][3]);} // 3 times
    {psi.inv_tchannel(no_ps,   6,  12,  17,  34,   4,   8, inv_r[2][4], inv_r[2][5]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  12,  17,  34,   8,   4, inv_r[2][6], inv_r[2][7]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   2,   1,  32,  28, inv_r[2][8], inv_r[2][9]);} // 2 times
    {psi.inv_tchannel_opt(no_ps,   0,  28,   1,  34,  12,  16, inv_r[2][10], inv_r[2][11]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  44,   2,  17,  12,  32, inv_r[2][12], inv_r[2][13]);} // 1 times

    {psi.inv_decay(no_ps,  12,   4,   8, inv_r[3][0], inv_r[3][1]);} // 5 times
    {psi.inv_decay(no_ps,  28,  16,  12, inv_r[3][2], inv_r[3][3]);} // 1 times
    {psi.inv_decay(no_ps,  44,  32,  12, inv_r[3][4], inv_r[3][5]);} // 1 times

    if (psi_weight_IS == 2){
      int n_random = 3 * psi_n_particle - 4;;
      vector<vector<double> > inv_channel_r(7, vector<double> (8));
      vector<double> g_IS(7, 1.);
      for (int i_c = 0; i_c < 7; i_c++){
        int ran = 0;
        if (psi_MC_alpha[zero + i_c] != 0.){
          for (int i_p = 0; i_p < psi_c_p[no_ps][i_c].size(); i_p++){
            inv_channel_r[i_c][i_p] = inv_r[0][psi_c_p[no_ps][i_c][i_p]];
          }
          ran += psi_c_p[no_ps][i_c].size();
          for (int i_f = 0; i_f < psi_c_f[no_ps][i_c].size(); i_f++){
            inv_channel_r[i_c][ran + i_f] = inv_r[1][psi_c_f[no_ps][i_c][i_f]];
          }
          ran += psi_c_f[no_ps][i_c].size();
          for (int i_t = 0; i_t < psi_c_t[no_ps][i_c].size(); i_t++){
            inv_channel_r[i_c][ran + 2 * i_t] = inv_r[2][2 * psi_c_t[no_ps][i_c][i_t]];
            inv_channel_r[i_c][ran + 2 * i_t + 1] = inv_r[2][2 * psi_c_t[no_ps][i_c][i_t] + 1];
          }
          ran += 2 * psi_c_t[no_ps][i_c].size();
          for (int i_d = 0; i_d < psi_c_d[no_ps][i_c].size(); i_d++){
            inv_channel_r[i_c][ran + 2 * i_d] = inv_r[3][2 * psi_c_d[no_ps][i_c][i_d]];
            inv_channel_r[i_c][ran + 2 * i_d + 1] = inv_r[3][2 * psi_c_d[no_ps][i_c][i_d] + 1];
          }
          ran += 2 * psi_c_d[no_ps][i_c].size();
          for (int i_r = 0; i_r < 8; i_r++){
            psi_phasespace_randoms[(zero + i_c) * n_random + i_r]->get_g_IS(inv_channel_r[i_c][i_r], g_IS_temp);
            g_IS[i_c] *= g_IS_temp;
          }
        }
      }
      if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = g_IS[  0] * psi_g_f[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_t[no_ps][  2];}
      if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = g_IS[  1] * psi_g_f[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_t[no_ps][  3];}
      if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = g_IS[  2] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][  2] * psi_g_t[no_ps][  4] * psi_g_t[no_ps][  5] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = g_IS[  3] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  1] * psi_g_t[no_ps][  4] * psi_g_d[no_ps][  1] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = g_IS[  4] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = g_IS[  5] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  2] * psi_g_t[no_ps][  0] * psi_g_d[no_ps][  2] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = g_IS[  6] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  6] * psi_g_d[no_ps][  0];}

    }
    else if (psi_weight_IS == 4){
      vector<vector<double> > g_IS(4);
      g_IS[0].resize(3, 1.);
      g_IS[1].resize(3, 1.);
      g_IS[2].resize(7, 1.);
      g_IS[3].resize(3, 1.);

      for (int i_m = 0; i_m < 4; i_m++){
        for (int i_r = 0; i_r < inv_r[i_m].size(); i_r++){
          psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][i_m] + i_r]->get_g_IS(inv_r[i_m][i_r], g_IS_temp);
          int i_g = i_r;
          if (i_m > 1){i_g = i_r / 2;}
          g_IS[i_m][i_g] *= g_IS_temp;
        }
      }
      if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = psi_g_f[no_ps][  0] * g_IS[1][  0] * psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  1] * g_IS[2][  1] * psi_g_t[no_ps][  2] * g_IS[2][  2];}
      if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = psi_g_f[no_ps][  0] * g_IS[1][  0] * psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  1] * g_IS[2][  1] * psi_g_t[no_ps][  3] * g_IS[2][  3];}
      if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][  2] * g_IS[1][  2] * psi_g_t[no_ps][  4] * g_IS[2][  4] * psi_g_t[no_ps][  5] * g_IS[2][  5] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_t[no_ps][  4] * g_IS[2][  4] * psi_g_d[no_ps][  1] * g_IS[3][  1] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  1] * g_IS[2][  1] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_d[no_ps][  2] * g_IS[3][  2] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
    }
  }

  else {
    if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = psi_g_f[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_t[no_ps][  2];}
    if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = psi_g_f[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_t[no_ps][  3];}
    if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][  2] * psi_g_t[no_ps][  4] * psi_g_t[no_ps][  5] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  1] * psi_g_t[no_ps][  4] * psi_g_d[no_ps][  1] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  2] * psi_g_t[no_ps][  0] * psi_g_d[no_ps][  2] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  6] * psi_g_d[no_ps][  0];}
  }

/*
  if (psi_MC_alpha[   0] != 0.){ag_channel_1111(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  if (psi_MC_alpha[   1] != 0.){ag_channel_1111(g, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  if (psi_MC_alpha[   2] != 0.){ag_channel_2v11(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[   3] != 0.){ag_channel_3yv1(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[   4] != 0.){ag_channel_12v1(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[   5] != 0.){ag_channel_13yv(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[   6] != 0.){ag_channel_112v(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  // channel    0: 1111(1, 2, 16, 4, 8, 32;   0,   0,   6)
  // channel    1: 1111(1, 2, 16, 8, 4, 32;   0,   0,   6)
  // channel    2: 2v11(1, 2, 4, 8, 16, 32;   0,   0,   0)
  // channel    3: 3yv1(1, 2, 16, 4, 8, 32;   0,   0,   0)
  // channel    4: 12v1(1, 2, 16, 4, 8, 32;   0,   0,   0)
  // channel    5: 13yv(1, 2, 16, 32, 4, 8;   0,   0,   0)
  // channel    6: 112v(1, 2, 16, 32, 4, 8;   0,   0,   0)
*/

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
