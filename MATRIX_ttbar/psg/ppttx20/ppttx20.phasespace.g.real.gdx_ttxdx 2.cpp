#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ag_psp_300_gdx_ttxdx(int no_ps, int zero, phasespace_set & psi){
  static Logger logger("ppttx20_ag_psp_300_gdx_ttxdx");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (psi_xbp[no_ps][5] == nullvector){psi_xbp[no_ps][5] = psi_xbp[no_ps][1] - psi_xbp[no_ps][4];}
  if (psi_xbs[no_ps][5] == 0.){psi_xbs[no_ps][5] = psi_xbp[no_ps][5].m2();}
  if (psi_xbs[no_ps][5] > 0.){psi_xbs[no_ps][5] = 0.;}
  if (psi_xbp[no_ps][9] == nullvector){psi_xbp[no_ps][9] = psi_xbp[no_ps][1] - psi_xbp[no_ps][8];}
  if (psi_xbs[no_ps][9] == 0.){psi_xbs[no_ps][9] = psi_xbp[no_ps][9].m2();}
  if (psi_xbs[no_ps][9] > 0.){psi_xbs[no_ps][9] = 0.;}
  if (psi_xbp[no_ps][12] == nullvector){psi_xbp[no_ps][12] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8];}
  if (psi_xbs[no_ps][12] == 0.){psi_xbs[no_ps][12] = psi_xbp[no_ps][12].m2();}
  if (psi_xbs[no_ps][12] < 0.){psi_xbs[no_ps][12] = 0.;}
  if (psi_xbp[no_ps][17] == nullvector){psi_xbp[no_ps][17] = psi_xbp[no_ps][1] - psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][17] == 0.){psi_xbs[no_ps][17] = psi_xbp[no_ps][17].m2();}
  if (psi_xbs[no_ps][17] > 0.){psi_xbs[no_ps][17] = 0.;}
  if (psi_xbp[no_ps][18] == nullvector){psi_xbp[no_ps][18] = psi_xbp[no_ps][2] - psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][18] == 0.){psi_xbs[no_ps][18] = psi_xbp[no_ps][18].m2();}
  if (psi_xbs[no_ps][18] > 0.){psi_xbs[no_ps][18] = 0.;}
  if (psi_xbp[no_ps][20] == nullvector){psi_xbp[no_ps][20] = psi_xbp[no_ps][4] + psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][20] == 0.){psi_xbs[no_ps][20] = psi_xbp[no_ps][20].m2();}
  if (psi_xbs[no_ps][20] < 0.){psi_xbs[no_ps][20] = 0.;}
  if (psi_xbp[no_ps][24] == nullvector){psi_xbp[no_ps][24] = psi_xbp[no_ps][8] + psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][24] == 0.){psi_xbs[no_ps][24] = psi_xbp[no_ps][24].m2();}
  if (psi_xbs[no_ps][24] < 0.){psi_xbs[no_ps][24] = 0.;}
  if (psi_xbp[no_ps][28] == nullvector){psi_xbp[no_ps][28] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8] + psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][28] == 0.){psi_xbs[no_ps][28] = psi_xbp[no_ps][28].m2();}
  if (psi_xbs[no_ps][28] < 0.){psi_xbs[no_ps][28] = 0.;}

  psi_v_smin[no_ps][  0] = psi_smin_opt[no_ps][12]; // 3 times
  psi_v_smin[no_ps][  1] = psi_smin_opt[no_ps][24]; // 1 times
  psi_v_smin[no_ps][  2] = psi_smin_opt[no_ps][20]; // 1 times

  psi_v_smax[no_ps][  0] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2); // 3 times
  psi_v_smax[no_ps][  1] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2); // 1 times
  psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2); // 1 times

  //  vector<double> g_p(1);
  {psi_g_p[no_ps][  0] = psi.g_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0]); /* npp 2 */} // 3 times

  //  vector<double> g_f(2);
  {psi_g_f[no_ps][  0] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1]); /* nfp 2 */} // 1 times
  {psi_g_f[no_ps][  1] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2]); /* nfp 2 */} // 1 times

  //  vector<double> g_t(6);
  {psi_g_t[no_ps][  0] = psi.g_tchannel_opt(no_ps,   6,  28,   1,   2,   4,  24); /* ntp 0 */} // 1 times
  {psi_g_t[no_ps][  1] = psi.g_tchannel_opt(no_ps,   0,  24,   2,   5,  16,   8); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][  2] = psi.g_tchannel_opt(no_ps,   6,  28,   1,   2,   8,  20); /* ntp 0 */} // 1 times
  {psi_g_t[no_ps][  3] = psi.g_tchannel_opt(no_ps,   0,  20,   2,   9,  16,   4); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][  4] = psi.g_tchannel_opt(no_ps,   0,  28,   2,   1,  16,  12); /* ntp 0 */} // 1 times
  {psi_g_t[no_ps][  5] = psi.g_tchannel_opt(no_ps,   0,  28,   1,   2,  16,  12); /* ntp 0 */} // 1 times

  //  vector<double> g_d(2);
  {psi_g_d[no_ps][  0] = psi.g_decay(no_ps,  12,   4,   8); /* ndp 2 */} // 3 times
  {psi_g_d[no_ps][  1] = psi.g_decay(no_ps,  28,  16,  12); /* ndp 3 */} // 1 times


  if (psi_weight_IS == 2 || psi_weight_IS == 4){
    double g_IS_temp;
    vector<vector<double> > inv_r(4);
    inv_r[0].resize(1);
    inv_r[1].resize(2);
    inv_r[2].resize(12);
    inv_r[3].resize(4);

    {psi.inv_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], inv_r[0][0]);} // 3 times

    {psi.inv_timelikeinvariant(no_ps,  24, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], inv_r[1][0]);} // 1 times
    {psi.inv_timelikeinvariant(no_ps,  20, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], inv_r[1][1]);} // 1 times

    {psi.inv_tchannel_opt(no_ps,   6,  28,   1,   2,   4,  24, inv_r[2][0], inv_r[2][1]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  24,   2,   5,  16,   8, inv_r[2][2], inv_r[2][3]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  28,   1,   2,   8,  20, inv_r[2][4], inv_r[2][5]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  20,   2,   9,  16,   4, inv_r[2][6], inv_r[2][7]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  28,   2,   1,  16,  12, inv_r[2][8], inv_r[2][9]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  28,   1,   2,  16,  12, inv_r[2][10], inv_r[2][11]);} // 1 times

    {psi.inv_decay(no_ps,  12,   4,   8, inv_r[3][0], inv_r[3][1]);} // 3 times
    {psi.inv_decay(no_ps,  28,  16,  12, inv_r[3][2], inv_r[3][3]);} // 1 times

    if (psi_weight_IS == 2){
      int n_random = 3 * psi_n_particle - 4;;
      vector<vector<double> > inv_channel_r(5, vector<double> (5));
      vector<double> g_IS(5, 1.);
      for (int i_c = 0; i_c < 5; i_c++){
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
          for (int i_r = 0; i_r < 5; i_r++){
            psi_phasespace_randoms[(zero + i_c) * n_random + i_r]->get_g_IS(inv_channel_r[i_c][i_r], g_IS_temp);
            g_IS[i_c] *= g_IS_temp;
          }
        }
      }
      if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = g_IS[  0] * psi_g_f[no_ps][  0] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1];}
      if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = g_IS[  1] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  2] * psi_g_t[no_ps][  3];}
      if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = g_IS[  2] * psi_g_p[no_ps][  0] * psi_g_t[no_ps][  4] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = g_IS[  3] * psi_g_p[no_ps][  0] * psi_g_d[no_ps][  1] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = g_IS[  4] * psi_g_p[no_ps][  0] * psi_g_t[no_ps][  5] * psi_g_d[no_ps][  0];}

    }
    else if (psi_weight_IS == 4){
      vector<vector<double> > g_IS(4);
      g_IS[0].resize(1, 1.);
      g_IS[1].resize(2, 1.);
      g_IS[2].resize(6, 1.);
      g_IS[3].resize(2, 1.);

      for (int i_m = 0; i_m < 4; i_m++){
        for (int i_r = 0; i_r < inv_r[i_m].size(); i_r++){
          psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][i_m] + i_r]->get_g_IS(inv_r[i_m][i_r], g_IS_temp);
          int i_g = i_r;
          if (i_m > 1){i_g = i_r / 2;}
          g_IS[i_m][i_g] *= g_IS_temp;
        }
      }
      if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = psi_g_f[no_ps][  0] * g_IS[1][  0] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  1] * g_IS[2][  1];}
      if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][  2] * g_IS[2][  2] * psi_g_t[no_ps][  3] * g_IS[2][  3];}
      if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_t[no_ps][  4] * g_IS[2][  4] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_d[no_ps][  1] * g_IS[3][  1] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_t[no_ps][  5] * g_IS[2][  5] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
    }
  }

  else {
    if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = psi_g_f[no_ps][  0] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1];}
    if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = psi_g_f[no_ps][  1] * psi_g_t[no_ps][  2] * psi_g_t[no_ps][  3];}
    if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = psi_g_p[no_ps][  0] * psi_g_t[no_ps][  4] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = psi_g_p[no_ps][  0] * psi_g_d[no_ps][  1] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = psi_g_p[no_ps][  0] * psi_g_t[no_ps][  5] * psi_g_d[no_ps][  0];}
  }

/*
  if (psi_MC_alpha[   0] != 0.){ag_channel_111(g, 1, 2, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0, psi);}
  if (psi_MC_alpha[   1] != 0.){ag_channel_111(g, 1, 2, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0, psi);}
  if (psi_MC_alpha[   2] != 0.){ag_channel_2v1(g, 1, 2, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0, psi);}
  if (psi_MC_alpha[   3] != 0.){ag_channel_3yv(g, 1, 2, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0, psi);}
  if (psi_MC_alpha[   4] != 0.){ag_channel_12v(g, 1, 2, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0, psi);}
  // channel    0: 111(1, 2, 4, 8, 16;   6,   0)
  // channel    1: 111(1, 2, 8, 4, 16;   6,   0)
  // channel    2: 2v1(1, 2, 4, 8, 16;   0,   0)
  // channel    3: 3yv(1, 2, 16, 4, 8;   0,   0)
  // channel    4: 12v(1, 2, 16, 4, 8;   0,   0)
*/

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
