#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ag_psp_200_uux_ttx(int no_ps, int zero, phasespace_set & psi){
  static Logger logger("ppttx20_ag_psp_200_uux_ttx");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (psi_xbp[no_ps][12] == nullvector){psi_xbp[no_ps][12] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8];}
  if (psi_xbs[no_ps][12] == 0.){psi_xbs[no_ps][12] = psi_xbp[no_ps][12].m2();}
  if (psi_xbs[no_ps][12] < 0.){psi_xbs[no_ps][12] = 0.;}



  //  vector<double> g_p(0);

  //  vector<double> g_f(0);

  //  vector<double> g_t(0);

  //  vector<double> g_d(1);
  {psi_g_d[no_ps][  0] = psi.g_decay(no_ps,  12,   4,   8); /* ndp 2 */} // 1 times


  if (psi_weight_IS == 2 || psi_weight_IS == 4){
    double g_IS_temp;
    vector<vector<double> > inv_r(4);
    inv_r[0].resize(0);
    inv_r[1].resize(0);
    inv_r[2].resize(0);
    inv_r[3].resize(2);




    {psi.inv_decay(no_ps,  12,   4,   8, inv_r[3][0], inv_r[3][1]);} // 1 times

    if (psi_weight_IS == 2){
      int n_random = 3 * psi_n_particle - 4;;
      vector<vector<double> > inv_channel_r(1, vector<double> (2));
      vector<double> g_IS(1, 1.);
      for (int i_c = 0; i_c < 1; i_c++){
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
          for (int i_r = 0; i_r < 2; i_r++){
            psi_phasespace_randoms[(zero + i_c) * n_random + i_r]->get_g_IS(inv_channel_r[i_c][i_r], g_IS_temp);
            g_IS[i_c] *= g_IS_temp;
          }
        }
      }
      if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = g_IS[  0] * psi_g_d[no_ps][  0];}

    }
    else if (psi_weight_IS == 4){
      vector<vector<double> > g_IS(4);
      g_IS[0].resize(0, 1.);
      g_IS[1].resize(0, 1.);
      g_IS[2].resize(0, 1.);
      g_IS[3].resize(1, 1.);

      for (int i_m = 0; i_m < 4; i_m++){
        for (int i_r = 0; i_r < inv_r[i_m].size(); i_r++){
          psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][i_m] + i_r]->get_g_IS(inv_r[i_m][i_r], g_IS_temp);
          int i_g = i_r;
          if (i_m > 1){i_g = i_r / 2;}
          g_IS[i_m][i_g] *= g_IS_temp;
        }
      }
      if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = psi_g_d[no_ps][  0] * g_IS[3][  0];}
    }
  }

  else {
    if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = psi_g_d[no_ps][  0];}
  }

/*
  if (psi_MC_alpha[   0] != 0.){ag_channel_2v(g, 1, 2, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0, psi);}
  // channel    0: 2v(1, 2, 4, 8;   0)
*/

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
