#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ag_psp_400_bbx_ttxuux(int no_ps, int zero, phasespace_set & psi){
  static Logger logger("ppttx20_ag_psp_400_bbx_ttxuux");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (psi_xbp[no_ps][12] == nullvector){psi_xbp[no_ps][12] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8];}
  if (psi_xbs[no_ps][12] == 0.){psi_xbs[no_ps][12] = psi_xbp[no_ps][12].m2();}
  if (psi_xbs[no_ps][12] < 0.){psi_xbs[no_ps][12] = 0.;}
  if (psi_xbsqrts[no_ps][12] == 0.){psi_xbsqrts[no_ps][12] = sqrt(psi_xbs[no_ps][12]);}
  if (psi_xbp[no_ps][48] == nullvector){psi_xbp[no_ps][48] = psi_xbp[no_ps][16] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][48] == 0.){psi_xbs[no_ps][48] = psi_xbp[no_ps][48].m2();}
  if (psi_xbs[no_ps][48] < 0.){psi_xbs[no_ps][48] = 0.;}
  if (psi_xbsqrts[no_ps][48] == 0.){psi_xbsqrts[no_ps][48] = sqrt(psi_xbs[no_ps][48]);}
  if (psi_xbp[no_ps][13] == nullvector){psi_xbp[no_ps][13] = psi_xbp[no_ps][1] - psi_xbp[no_ps][4] - psi_xbp[no_ps][8];}
  if (psi_xbs[no_ps][13] == 0.){psi_xbs[no_ps][13] = psi_xbp[no_ps][13].m2();}
  if (psi_xbs[no_ps][13] > 0.){psi_xbs[no_ps][13] = 0.;}
  if (psi_xbp[no_ps][28] == nullvector){psi_xbp[no_ps][28] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8] + psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][28] == 0.){psi_xbs[no_ps][28] = psi_xbp[no_ps][28].m2();}
  if (psi_xbs[no_ps][28] < 0.){psi_xbs[no_ps][28] = 0.;}
  if (psi_xbp[no_ps][44] == nullvector){psi_xbp[no_ps][44] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][44] == 0.){psi_xbs[no_ps][44] = psi_xbp[no_ps][44].m2();}
  if (psi_xbs[no_ps][44] < 0.){psi_xbs[no_ps][44] = 0.;}
  if (psi_xbp[no_ps][49] == nullvector){psi_xbp[no_ps][49] = psi_xbp[no_ps][1] - psi_xbp[no_ps][16] - psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][49] == 0.){psi_xbs[no_ps][49] = psi_xbp[no_ps][49].m2();}
  if (psi_xbs[no_ps][49] > 0.){psi_xbs[no_ps][49] = 0.;}
  if (psi_xbp[no_ps][52] == nullvector){psi_xbp[no_ps][52] = psi_xbp[no_ps][4] + psi_xbp[no_ps][16] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][52] == 0.){psi_xbs[no_ps][52] = psi_xbp[no_ps][52].m2();}
  if (psi_xbs[no_ps][52] < 0.){psi_xbs[no_ps][52] = 0.;}
  if (psi_xbp[no_ps][56] == nullvector){psi_xbp[no_ps][56] = psi_xbp[no_ps][8] + psi_xbp[no_ps][16] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][56] == 0.){psi_xbs[no_ps][56] = psi_xbp[no_ps][56].m2();}
  if (psi_xbs[no_ps][56] < 0.){psi_xbs[no_ps][56] = 0.;}
  if (psi_xbp[no_ps][60] == nullvector){psi_xbp[no_ps][60] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8] + psi_xbp[no_ps][16] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][60] == 0.){psi_xbs[no_ps][60] = psi_xbp[no_ps][60].m2();}
  if (psi_xbs[no_ps][60] < 0.){psi_xbs[no_ps][60] = 0.;}

  psi_v_smin[no_ps][  0] = psi_smin_opt[no_ps][48]; // 5 times
  psi_v_smin[no_ps][  1] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][56]); // 1 times
  psi_v_smin[no_ps][  2] = psi_smin_opt[no_ps][12]; // 5 times
  psi_v_smin[no_ps][  3] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][44]); // 1 times
  psi_v_smin[no_ps][  4] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][28]); // 1 times
  psi_v_smin[no_ps][  5] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][52]); // 1 times

  psi_v_smax[no_ps][  0] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][12], 2); // 4 times
  psi_v_smax[no_ps][  1] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2); // 1 times
  psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][48], 2); // 3 times
  psi_v_smax[no_ps][  3] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2); // 1 times
  psi_v_smax[no_ps][  4] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2); // 1 times
  psi_v_smax[no_ps][  5] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2); // 1 times
  psi_v_smax[no_ps][  6] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][48], 2); // 2 times
  psi_v_smax[no_ps][  7] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][12], 2); // 1 times

  //  vector<double> g_p(8);
  {psi_g_p[no_ps][  0] = psi.g_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0]); /* npp 2 */} // 4 times
  {psi_g_p[no_ps][  1] = psi.g_propagator(no_ps,   6,  56, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1]); /* npp 3 */} // 1 times
  {psi_g_p[no_ps][  2] = psi.g_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2]); /* npp 2 */} // 3 times
  {psi_g_p[no_ps][  3] = psi.g_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3]); /* npp 3 */} // 1 times
  {psi_g_p[no_ps][  4] = psi.g_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4]); /* npp 3 */} // 1 times
  {psi_g_p[no_ps][  5] = psi.g_propagator(no_ps,   6,  52, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  5]); /* npp 3 */} // 1 times
  {psi_g_p[no_ps][  6] = psi.g_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  6]); /* npp 2 */} // 2 times
  {psi_g_p[no_ps][  7] = psi.g_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  7]); /* npp 2 */} // 1 times

  //  vector<double> g_f(0);

  //  vector<double> g_t(2);
  {psi_g_t[no_ps][  0] = psi.g_tchannel_opt(no_ps,   5,  60,   1,   2,  12,  48); /* ntp 0 */} // 1 times
  {psi_g_t[no_ps][  1] = psi.g_tchannel_opt(no_ps,   5,  60,   1,   2,  48,  12); /* ntp 0 */} // 1 times

  //  vector<double> g_d(11);
  {psi_g_d[no_ps][  0] = psi.g_decay(no_ps,  60,   4,  56); /* ndp 4 */} // 1 times
  {psi_g_d[no_ps][  1] = psi.g_decay(no_ps,  56,   8,  48); /* ndp 3 */} // 1 times
  {psi_g_d[no_ps][  2] = psi.g_decay(no_ps,  48,  16,  32); /* ndp 2 */} // 5 times
  {psi_g_d[no_ps][  3] = psi.g_decay(no_ps,  60,  16,  44); /* ndp 4 */} // 1 times
  {psi_g_d[no_ps][  4] = psi.g_decay(no_ps,  44,  32,  12); /* ndp 3 */} // 1 times
  {psi_g_d[no_ps][  5] = psi.g_decay(no_ps,  12,   4,   8); /* ndp 2 */} // 5 times
  {psi_g_d[no_ps][  6] = psi.g_decay(no_ps,  60,  32,  28); /* ndp 4 */} // 1 times
  {psi_g_d[no_ps][  7] = psi.g_decay(no_ps,  28,  16,  12); /* ndp 3 */} // 1 times
  {psi_g_d[no_ps][  8] = psi.g_decay(no_ps,  60,   8,  52); /* ndp 4 */} // 1 times
  {psi_g_d[no_ps][  9] = psi.g_decay(no_ps,  52,   4,  48); /* ndp 3 */} // 1 times
  {psi_g_d[no_ps][ 10] = psi.g_decay(no_ps,  60,  12,  48); /* ndp 4 */} // 1 times


  if (psi_weight_IS == 2 || psi_weight_IS == 4){
    double g_IS_temp;
    vector<vector<double> > inv_r(4);
    inv_r[0].resize(8);
    inv_r[1].resize(0);
    inv_r[2].resize(4);
    inv_r[3].resize(22);

    {psi.inv_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], inv_r[0][0]);} // 4 times
    {psi.inv_propagator(no_ps,   6,  56, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], inv_r[0][1]);} // 1 times
    {psi.inv_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], inv_r[0][2]);} // 3 times
    {psi.inv_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3], inv_r[0][3]);} // 1 times
    {psi.inv_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4], inv_r[0][4]);} // 1 times
    {psi.inv_propagator(no_ps,   6,  52, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  5], inv_r[0][5]);} // 1 times
    {psi.inv_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  6], inv_r[0][6]);} // 2 times
    {psi.inv_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  7], inv_r[0][7]);} // 1 times


    {psi.inv_tchannel_opt(no_ps,   5,  60,   1,   2,  12,  48, inv_r[2][0], inv_r[2][1]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   5,  60,   1,   2,  48,  12, inv_r[2][2], inv_r[2][3]);} // 1 times

    {psi.inv_decay(no_ps,  60,   4,  56, inv_r[3][0], inv_r[3][1]);} // 1 times
    {psi.inv_decay(no_ps,  56,   8,  48, inv_r[3][2], inv_r[3][3]);} // 1 times
    {psi.inv_decay(no_ps,  48,  16,  32, inv_r[3][4], inv_r[3][5]);} // 5 times
    {psi.inv_decay(no_ps,  60,  16,  44, inv_r[3][6], inv_r[3][7]);} // 1 times
    {psi.inv_decay(no_ps,  44,  32,  12, inv_r[3][8], inv_r[3][9]);} // 1 times
    {psi.inv_decay(no_ps,  12,   4,   8, inv_r[3][10], inv_r[3][11]);} // 5 times
    {psi.inv_decay(no_ps,  60,  32,  28, inv_r[3][12], inv_r[3][13]);} // 1 times
    {psi.inv_decay(no_ps,  28,  16,  12, inv_r[3][14], inv_r[3][15]);} // 1 times
    {psi.inv_decay(no_ps,  60,   8,  52, inv_r[3][16], inv_r[3][17]);} // 1 times
    {psi.inv_decay(no_ps,  52,   4,  48, inv_r[3][18], inv_r[3][19]);} // 1 times
    {psi.inv_decay(no_ps,  60,  12,  48, inv_r[3][20], inv_r[3][21]);} // 1 times

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
      if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = g_IS[  0] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  1] * psi_g_d[no_ps][  0] * psi_g_d[no_ps][  1] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = g_IS[  1] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][  3] * psi_g_d[no_ps][  3] * psi_g_d[no_ps][  4] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = g_IS[  2] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][  4] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = g_IS[  3] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  5] * psi_g_d[no_ps][  8] * psi_g_d[no_ps][  9] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = g_IS[  4] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  6] * psi_g_d[no_ps][ 10] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = g_IS[  5] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  6] * psi_g_t[no_ps][  0] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = g_IS[  6] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][  7] * psi_g_t[no_ps][  1] * psi_g_d[no_ps][  2] * psi_g_d[no_ps][  5];}

    }
    else if (psi_weight_IS == 4){
      vector<vector<double> > g_IS(4);
      g_IS[0].resize(8, 1.);
      g_IS[1].resize(0, 1.);
      g_IS[2].resize(2, 1.);
      g_IS[3].resize(11, 1.);

      for (int i_m = 0; i_m < 4; i_m++){
        for (int i_r = 0; i_r < inv_r[i_m].size(); i_r++){
          psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][i_m] + i_r]->get_g_IS(inv_r[i_m][i_r], g_IS_temp);
          int i_g = i_r;
          if (i_m > 1){i_g = i_r / 2;}
          g_IS[i_m][i_g] *= g_IS_temp;
        }
      }
      if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_d[no_ps][  0] * g_IS[3][  0] * psi_g_d[no_ps][  1] * g_IS[3][  1] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_d[no_ps][  3] * g_IS[3][  3] * psi_g_d[no_ps][  4] * g_IS[3][  4] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_d[no_ps][  6] * g_IS[3][  6] * psi_g_d[no_ps][  7] * g_IS[3][  7] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_d[no_ps][  8] * g_IS[3][  8] * psi_g_d[no_ps][  9] * g_IS[3][  9] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  6] * g_IS[0][  6] * psi_g_d[no_ps][ 10] * g_IS[3][ 10] * psi_g_d[no_ps][  5] * g_IS[3][  5] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  6] * g_IS[0][  6] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_d[no_ps][  5] * g_IS[3][  5] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][  7] * g_IS[0][  7] * psi_g_t[no_ps][  1] * g_IS[2][  1] * psi_g_d[no_ps][  2] * g_IS[3][  2] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
    }
  }

  else {
    if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  1] * psi_g_d[no_ps][  0] * psi_g_d[no_ps][  1] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][  3] * psi_g_d[no_ps][  3] * psi_g_d[no_ps][  4] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][  4] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  5] * psi_g_d[no_ps][  8] * psi_g_d[no_ps][  9] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  6] * psi_g_d[no_ps][ 10] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  6] * psi_g_t[no_ps][  0] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][  7] * psi_g_t[no_ps][  1] * psi_g_d[no_ps][  2] * psi_g_d[no_ps][  5];}
  }

/*
  if (psi_MC_alpha[   0] != 0.){ag_channel_4yyv(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[   1] != 0.){ag_channel_4yyv(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[   2] != 0.){ag_channel_4yyv(g, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[   3] != 0.){ag_channel_4yyv(g, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[   4] != 0.){ag_channel_42vv(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[   5] != 0.){ag_channel_2v2v(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   5, psi);}
  if (psi_MC_alpha[   6] != 0.){ag_channel_2v2v(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   5, psi);}
  // channel    0: 4yyv(1, 2, 4, 8, 16, 32;   0,   6,   0)
  // channel    1: 4yyv(1, 2, 16, 32, 4, 8;   0,   0,   0)
  // channel    2: 4yyv(1, 2, 32, 16, 4, 8;   0,   0,   0)
  // channel    3: 4yyv(1, 2, 8, 4, 16, 32;   0,   6,   0)
  // channel    4: 42vv(1, 2, 4, 8, 16, 32;   0,   0,   0)
  // channel    5: 2v2v(1, 2, 4, 8, 16, 32;   0,   0,   5)
  // channel    6: 2v2v(1, 2, 16, 32, 4, 8;   0,   0,   5)
*/

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
