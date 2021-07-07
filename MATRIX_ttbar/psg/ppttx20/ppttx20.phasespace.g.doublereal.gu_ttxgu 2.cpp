#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ag_psp_400_gu_ttxgu(int no_ps, int zero, phasespace_set & psi){
  static Logger logger("ppttx20_ag_psp_400_gu_ttxgu");
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
  if (psi_xbsqrts[no_ps][12] == 0.){psi_xbsqrts[no_ps][12] = sqrt(psi_xbs[no_ps][12]);}
  if (psi_xbp[no_ps][17] == nullvector){psi_xbp[no_ps][17] = psi_xbp[no_ps][1] - psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][17] == 0.){psi_xbs[no_ps][17] = psi_xbp[no_ps][17].m2();}
  if (psi_xbs[no_ps][17] > 0.){psi_xbs[no_ps][17] = 0.;}
  if (psi_xbp[no_ps][18] == nullvector){psi_xbp[no_ps][18] = psi_xbp[no_ps][2] - psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][18] == 0.){psi_xbs[no_ps][18] = psi_xbp[no_ps][18].m2();}
  if (psi_xbs[no_ps][18] > 0.){psi_xbs[no_ps][18] = 0.;}
  if (psi_xbp[no_ps][20] == nullvector){psi_xbp[no_ps][20] = psi_xbp[no_ps][4] + psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][20] == 0.){psi_xbs[no_ps][20] = psi_xbp[no_ps][20].m2();}
  if (psi_xbs[no_ps][20] < 0.){psi_xbs[no_ps][20] = 0.;}
  if (psi_xbsqrts[no_ps][20] == 0.){psi_xbsqrts[no_ps][20] = sqrt(psi_xbs[no_ps][20]);}
  if (psi_xbp[no_ps][24] == nullvector){psi_xbp[no_ps][24] = psi_xbp[no_ps][8] + psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][24] == 0.){psi_xbs[no_ps][24] = psi_xbp[no_ps][24].m2();}
  if (psi_xbs[no_ps][24] < 0.){psi_xbs[no_ps][24] = 0.;}
  if (psi_xbsqrts[no_ps][24] == 0.){psi_xbsqrts[no_ps][24] = sqrt(psi_xbs[no_ps][24]);}
  if (psi_xbp[no_ps][33] == nullvector){psi_xbp[no_ps][33] = psi_xbp[no_ps][1] - psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][33] == 0.){psi_xbs[no_ps][33] = psi_xbp[no_ps][33].m2();}
  if (psi_xbs[no_ps][33] > 0.){psi_xbs[no_ps][33] = 0.;}
  if (psi_xbp[no_ps][34] == nullvector){psi_xbp[no_ps][34] = psi_xbp[no_ps][2] - psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][34] == 0.){psi_xbs[no_ps][34] = psi_xbp[no_ps][34].m2();}
  if (psi_xbs[no_ps][34] > 0.){psi_xbs[no_ps][34] = 0.;}
  if (psi_xbp[no_ps][36] == nullvector){psi_xbp[no_ps][36] = psi_xbp[no_ps][4] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][36] == 0.){psi_xbs[no_ps][36] = psi_xbp[no_ps][36].m2();}
  if (psi_xbs[no_ps][36] < 0.){psi_xbs[no_ps][36] = 0.;}
  if (psi_xbsqrts[no_ps][36] == 0.){psi_xbsqrts[no_ps][36] = sqrt(psi_xbs[no_ps][36]);}
  if (psi_xbp[no_ps][40] == nullvector){psi_xbp[no_ps][40] = psi_xbp[no_ps][8] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][40] == 0.){psi_xbs[no_ps][40] = psi_xbp[no_ps][40].m2();}
  if (psi_xbs[no_ps][40] < 0.){psi_xbs[no_ps][40] = 0.;}
  if (psi_xbsqrts[no_ps][40] == 0.){psi_xbsqrts[no_ps][40] = sqrt(psi_xbs[no_ps][40]);}
  if (psi_xbp[no_ps][48] == nullvector){psi_xbp[no_ps][48] = psi_xbp[no_ps][16] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][48] == 0.){psi_xbs[no_ps][48] = psi_xbp[no_ps][48].m2();}
  if (psi_xbs[no_ps][48] < 0.){psi_xbs[no_ps][48] = 0.;}
  if (psi_xbsqrts[no_ps][48] == 0.){psi_xbsqrts[no_ps][48] = sqrt(psi_xbs[no_ps][48]);}
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
  if (psi_xbp[no_ps][49] == nullvector){psi_xbp[no_ps][49] = psi_xbp[no_ps][1] - psi_xbp[no_ps][16] - psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][49] == 0.){psi_xbs[no_ps][49] = psi_xbp[no_ps][49].m2();}
  if (psi_xbs[no_ps][49] > 0.){psi_xbs[no_ps][49] = 0.;}
  if (psi_xbp[no_ps][50] == nullvector){psi_xbp[no_ps][50] = psi_xbp[no_ps][2] - psi_xbp[no_ps][16] - psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][50] == 0.){psi_xbs[no_ps][50] = psi_xbp[no_ps][50].m2();}
  if (psi_xbs[no_ps][50] > 0.){psi_xbs[no_ps][50] = 0.;}
  if (psi_xbp[no_ps][52] == nullvector){psi_xbp[no_ps][52] = psi_xbp[no_ps][4] + psi_xbp[no_ps][16] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][52] == 0.){psi_xbs[no_ps][52] = psi_xbp[no_ps][52].m2();}
  if (psi_xbs[no_ps][52] < 0.){psi_xbs[no_ps][52] = 0.;}
  if (psi_xbp[no_ps][56] == nullvector){psi_xbp[no_ps][56] = psi_xbp[no_ps][8] + psi_xbp[no_ps][16] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][56] == 0.){psi_xbs[no_ps][56] = psi_xbp[no_ps][56].m2();}
  if (psi_xbs[no_ps][56] < 0.){psi_xbs[no_ps][56] = 0.;}
  if (psi_xbp[no_ps][60] == nullvector){psi_xbp[no_ps][60] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8] + psi_xbp[no_ps][16] + psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][60] == 0.){psi_xbs[no_ps][60] = psi_xbp[no_ps][60].m2();}
  if (psi_xbs[no_ps][60] < 0.){psi_xbs[no_ps][60] = 0.;}

  psi_v_smin[no_ps][  0] = psi_smin_opt[no_ps][12]; // 17 times
  psi_v_smin[no_ps][  1] = psi_smin_opt[no_ps][20]; // 7 times
  psi_v_smin[no_ps][  2] = psi_smin_opt[no_ps][24]; // 7 times
  psi_v_smin[no_ps][  3] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][24], 2), psi_smin_opt[no_ps][28]); // 4 times
  psi_v_smin[no_ps][  4] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][44]); // 8 times
  psi_v_smin[no_ps][  5] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][28]); // 6 times
  psi_v_smin[no_ps][  6] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][20], 2), psi_smin_opt[no_ps][28]); // 4 times
  psi_v_smin[no_ps][  7] = psi_smin_opt[no_ps][48]; // 5 times
  psi_v_smin[no_ps][  8] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][24], 2), psi_smin_opt[no_ps][56]); // 3 times
  psi_v_smin[no_ps][  9] = psi_smin_opt[no_ps][40]; // 1 times
  psi_v_smin[no_ps][ 10] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][40], 2), psi_smin_opt[no_ps][56]); // 1 times
  psi_v_smin[no_ps][ 11] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][20], 2), psi_smin_opt[no_ps][52]); // 3 times
  psi_v_smin[no_ps][ 12] = psi_smin_opt[no_ps][36]; // 1 times
  psi_v_smin[no_ps][ 13] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][36], 2), psi_smin_opt[no_ps][52]); // 1 times
  psi_v_smin[no_ps][ 14] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][56]); // 1 times
  psi_v_smin[no_ps][ 15] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][52]); // 1 times

  psi_v_smax[no_ps][  0] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][48], 2); // 15 times
  psi_v_smax[no_ps][  1] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][40], 2); // 7 times
  psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][36], 2); // 7 times
  psi_v_smax[no_ps][  3] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2); // 14 times
  psi_v_smax[no_ps][  4] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2); // 8 times
  psi_v_smax[no_ps][  5] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][12], 2); // 4 times
  psi_v_smax[no_ps][  6] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][48], 2); // 2 times
  psi_v_smax[no_ps][  7] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][12], 2); // 1 times
  psi_v_smax[no_ps][  8] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2); // 5 times
  psi_v_smax[no_ps][  9] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][20], 2); // 1 times
  psi_v_smax[no_ps][ 10] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2); // 5 times
  psi_v_smax[no_ps][ 11] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][24], 2); // 1 times

  //  vector<double> g_p(10);
  {psi_g_p[no_ps][  0] = psi.g_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0]); /* npp 2 */} // 13 times
  {psi_g_p[no_ps][  1] = psi.g_propagator(no_ps,   6,  20, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1]); /* npp 2 */} // 5 times
  {psi_g_p[no_ps][  2] = psi.g_propagator(no_ps,   6,  24, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2]); /* npp 2 */} // 5 times
  {psi_g_p[no_ps][  3] = psi.g_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][  4] = psi.g_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][  5] = psi.g_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  3]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][  6] = psi.g_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  3]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][  7] = psi.g_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  7], psi_v_smax[no_ps][  5]); /* npp 2 */} // 4 times
  {psi_g_p[no_ps][  8] = psi.g_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  6]); /* npp 2 */} // 2 times
  {psi_g_p[no_ps][  9] = psi.g_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  7], psi_v_smax[no_ps][  7]); /* npp 2 */} // 1 times

  //  vector<double> g_f(15);
  {psi_g_f[no_ps][  0] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2]); /* nfp 2 */} // 2 times
  {psi_g_f[no_ps][  1] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  8], psi_v_smax[no_ps][  8]); /* nfp 3 */} // 3 times
  {psi_g_f[no_ps][  2] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  9], psi_v_smax[no_ps][  9]); /* nfp 2 */} // 1 times
  {psi_g_f[no_ps][  3] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 10], psi_v_smax[no_ps][  8]); /* nfp 3 */} // 1 times
  {psi_g_f[no_ps][  4] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0]); /* nfp 2 */} // 2 times
  {psi_g_f[no_ps][  5] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][  6] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1]); /* nfp 2 */} // 2 times
  {psi_g_f[no_ps][  7] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 11], psi_v_smax[no_ps][ 10]); /* nfp 3 */} // 3 times
  {psi_g_f[no_ps][  8] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 12], psi_v_smax[no_ps][ 11]); /* nfp 2 */} // 1 times
  {psi_g_f[no_ps][  9] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 13], psi_v_smax[no_ps][ 10]); /* nfp 3 */} // 1 times
  {psi_g_f[no_ps][ 10] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  3]); /* nfp 3 */} // 3 times
  {psi_g_f[no_ps][ 11] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  3]); /* nfp 3 */} // 1 times
  {psi_g_f[no_ps][ 12] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3]); /* nfp 3 */} // 1 times
  {psi_g_f[no_ps][ 13] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 14], psi_v_smax[no_ps][  8]); /* nfp 3 */} // 1 times
  {psi_g_f[no_ps][ 14] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 15], psi_v_smax[no_ps][ 10]); /* nfp 3 */} // 1 times

  //  vector<double> g_t(30);
  {psi_g_t[no_ps][  0] = psi.g_tchannel_opt(no_ps,   6,  60,   1,   2,   4,  56); /* ntp 0 */} // 5 times
  {psi_g_t[no_ps][  1] = psi.g_tchannel_opt(no_ps,   0,  56,   2,   5,  32,  24); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][  2] = psi.g_tchannel(no_ps,   0,  24,   5,  34,   8,  16); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][  3] = psi.g_tchannel_opt(no_ps,   0,  56,   2,   5,  16,  40); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][  4] = psi.g_tchannel(no_ps,   0,  40,   5,  18,   8,  32); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][  5] = psi.g_tchannel(no_ps,   6,  24,   5,  34,  16,   8); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][  6] = psi.g_tchannel_opt(no_ps,   0,  60,   1,   2,  16,  44); /* ntp 0 */} // 5 times
  {psi_g_t[no_ps][  7] = psi.g_tchannel_opt(no_ps,   0,  44,   2,  17,  32,  12); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][  8] = psi.g_tchannel(no_ps,   6,  12,  17,  34,   4,   8); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][  9] = psi.g_tchannel_opt(no_ps,   6,  60,   1,   2,   8,  52); /* ntp 0 */} // 5 times
  {psi_g_t[no_ps][ 10] = psi.g_tchannel_opt(no_ps,   0,  52,   2,   9,  32,  20); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][ 11] = psi.g_tchannel(no_ps,   0,  20,   9,  34,   4,  16); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 12] = psi.g_tchannel_opt(no_ps,   0,  52,   2,   9,  16,  36); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 13] = psi.g_tchannel(no_ps,   0,  36,   9,  18,   4,  32); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 14] = psi.g_tchannel(no_ps,   6,  20,   9,  34,  16,   4); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 15] = psi.g_tchannel(no_ps,   6,  12,  17,  34,   8,   4); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 16] = psi.g_tchannel_opt(no_ps,   0,  60,   2,   1,  32,  28); /* ntp 0 */} // 6 times
  {psi_g_t[no_ps][ 17] = psi.g_tchannel_opt(no_ps,   0,  28,   1,  34,  12,  16); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 18] = psi.g_tchannel_opt(no_ps,   0,  60,   2,   1,  16,  44); /* ntp 0 */} // 2 times
  {psi_g_t[no_ps][ 19] = psi.g_tchannel_opt(no_ps,   0,  44,   1,  18,  12,  32); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 20] = psi.g_tchannel_opt(no_ps,   6,  28,   1,  34,  20,   8); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 21] = psi.g_tchannel_opt(no_ps,   6,  28,   1,  34,  24,   4); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 22] = psi.g_tchannel_opt(no_ps,   0,  60,   1,   2,  32,  28); /* ntp 0 */} // 5 times
  {psi_g_t[no_ps][ 23] = psi.g_tchannel_opt(no_ps,   0,  28,   2,  33,  16,  12); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 24] = psi.g_tchannel_opt(no_ps,   0,  56,   2,   5,  48,   8); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 25] = psi.g_tchannel_opt(no_ps,   0,  44,   2,  17,  12,  32); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 26] = psi.g_tchannel_opt(no_ps,   0,  28,   2,  33,  12,  16); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 27] = psi.g_tchannel_opt(no_ps,   0,  52,   2,   9,  48,   4); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 28] = psi.g_tchannel_opt(no_ps,   0,  60,   1,   2,  12,  48); /* ntp 0 */} // 1 times
  {psi_g_t[no_ps][ 29] = psi.g_tchannel_opt(no_ps,   0,  60,   1,   2,  48,  12); /* ntp 0 */} // 1 times

  //  vector<double> g_d(11);
  {psi_g_d[no_ps][  0] = psi.g_decay(no_ps,  12,   4,   8); /* ndp 2 */} // 15 times
  {psi_g_d[no_ps][  1] = psi.g_decay(no_ps,  20,   4,  16); /* ndp 2 */} // 5 times
  {psi_g_d[no_ps][  2] = psi.g_decay(no_ps,  24,   8,  16); /* ndp 2 */} // 5 times
  {psi_g_d[no_ps][  3] = psi.g_decay(no_ps,  28,   4,  24); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][  4] = psi.g_decay(no_ps,  44,  32,  12); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][  5] = psi.g_decay(no_ps,  28,  16,  12); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][  6] = psi.g_decay(no_ps,  28,   8,  20); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][  7] = psi.g_decay(no_ps,  60,  32,  28); /* ndp 4 */} // 3 times
  {psi_g_d[no_ps][  8] = psi.g_decay(no_ps,  60,  16,  44); /* ndp 4 */} // 1 times
  {psi_g_d[no_ps][  9] = psi.g_decay(no_ps,  60,  12,  48); /* ndp 4 */} // 1 times
  {psi_g_d[no_ps][ 10] = psi.g_decay(no_ps,  48,  16,  32); /* ndp 2 */} // 5 times


  if (psi_weight_IS == 2 || psi_weight_IS == 4){
    double g_IS_temp;
    vector<vector<double> > inv_r(4);
    inv_r[0].resize(10);
    inv_r[1].resize(15);
    inv_r[2].resize(60);
    inv_r[3].resize(22);

    {psi.inv_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], inv_r[0][0]);} // 13 times
    {psi.inv_propagator(no_ps,   6,  20, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], inv_r[0][1]);} // 5 times
    {psi.inv_propagator(no_ps,   6,  24, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], inv_r[0][2]);} // 5 times
    {psi.inv_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3], inv_r[0][3]);} // 3 times
    {psi.inv_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4], inv_r[0][4]);} // 3 times
    {psi.inv_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  3], inv_r[0][5]);} // 3 times
    {psi.inv_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  3], inv_r[0][6]);} // 3 times
    {psi.inv_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  7], psi_v_smax[no_ps][  5], inv_r[0][7]);} // 4 times
    {psi.inv_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  6], inv_r[0][8]);} // 2 times
    {psi.inv_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  7], psi_v_smax[no_ps][  7], inv_r[0][9]);} // 1 times

    {psi.inv_timelikeinvariant(no_ps,  24, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], inv_r[1][0]);} // 2 times
    {psi.inv_timelikeinvariant(no_ps,  56, psi_v_smin[no_ps][  8], psi_v_smax[no_ps][  8], inv_r[1][1]);} // 3 times
    {psi.inv_timelikeinvariant(no_ps,  40, psi_v_smin[no_ps][  9], psi_v_smax[no_ps][  9], inv_r[1][2]);} // 1 times
    {psi.inv_timelikeinvariant(no_ps,  56, psi_v_smin[no_ps][ 10], psi_v_smax[no_ps][  8], inv_r[1][3]);} // 1 times
    {psi.inv_timelikeinvariant(no_ps,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], inv_r[1][4]);} // 2 times
    {psi.inv_timelikeinvariant(no_ps,  44, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4], inv_r[1][5]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  20, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], inv_r[1][6]);} // 2 times
    {psi.inv_timelikeinvariant(no_ps,  52, psi_v_smin[no_ps][ 11], psi_v_smax[no_ps][ 10], inv_r[1][7]);} // 3 times
    {psi.inv_timelikeinvariant(no_ps,  36, psi_v_smin[no_ps][ 12], psi_v_smax[no_ps][ 11], inv_r[1][8]);} // 1 times
    {psi.inv_timelikeinvariant(no_ps,  52, psi_v_smin[no_ps][ 13], psi_v_smax[no_ps][ 10], inv_r[1][9]);} // 1 times
    {psi.inv_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  3], inv_r[1][10]);} // 3 times
    {psi.inv_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  3], inv_r[1][11]);} // 1 times
    {psi.inv_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3], inv_r[1][12]);} // 1 times
    {psi.inv_timelikeinvariant(no_ps,  56, psi_v_smin[no_ps][ 14], psi_v_smax[no_ps][  8], inv_r[1][13]);} // 1 times
    {psi.inv_timelikeinvariant(no_ps,  52, psi_v_smin[no_ps][ 15], psi_v_smax[no_ps][ 10], inv_r[1][14]);} // 1 times

    {psi.inv_tchannel_opt(no_ps,   6,  60,   1,   2,   4,  56, inv_r[2][0], inv_r[2][1]);} // 5 times
    {psi.inv_tchannel_opt(no_ps,   0,  56,   2,   5,  32,  24, inv_r[2][2], inv_r[2][3]);} // 3 times
    {psi.inv_tchannel(no_ps,   0,  24,   5,  34,   8,  16, inv_r[2][4], inv_r[2][5]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  56,   2,   5,  16,  40, inv_r[2][6], inv_r[2][7]);} // 1 times
    {psi.inv_tchannel(no_ps,   0,  40,   5,  18,   8,  32, inv_r[2][8], inv_r[2][9]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  24,   5,  34,  16,   8, inv_r[2][10], inv_r[2][11]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   1,   2,  16,  44, inv_r[2][12], inv_r[2][13]);} // 5 times
    {psi.inv_tchannel_opt(no_ps,   0,  44,   2,  17,  32,  12, inv_r[2][14], inv_r[2][15]);} // 3 times
    {psi.inv_tchannel(no_ps,   6,  12,  17,  34,   4,   8, inv_r[2][16], inv_r[2][17]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  60,   1,   2,   8,  52, inv_r[2][18], inv_r[2][19]);} // 5 times
    {psi.inv_tchannel_opt(no_ps,   0,  52,   2,   9,  32,  20, inv_r[2][20], inv_r[2][21]);} // 3 times
    {psi.inv_tchannel(no_ps,   0,  20,   9,  34,   4,  16, inv_r[2][22], inv_r[2][23]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  52,   2,   9,  16,  36, inv_r[2][24], inv_r[2][25]);} // 1 times
    {psi.inv_tchannel(no_ps,   0,  36,   9,  18,   4,  32, inv_r[2][26], inv_r[2][27]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  20,   9,  34,  16,   4, inv_r[2][28], inv_r[2][29]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  12,  17,  34,   8,   4, inv_r[2][30], inv_r[2][31]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   2,   1,  32,  28, inv_r[2][32], inv_r[2][33]);} // 6 times
    {psi.inv_tchannel_opt(no_ps,   0,  28,   1,  34,  12,  16, inv_r[2][34], inv_r[2][35]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   2,   1,  16,  44, inv_r[2][36], inv_r[2][37]);} // 2 times
    {psi.inv_tchannel_opt(no_ps,   0,  44,   1,  18,  12,  32, inv_r[2][38], inv_r[2][39]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  28,   1,  34,  20,   8, inv_r[2][40], inv_r[2][41]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  28,   1,  34,  24,   4, inv_r[2][42], inv_r[2][43]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   1,   2,  32,  28, inv_r[2][44], inv_r[2][45]);} // 5 times
    {psi.inv_tchannel_opt(no_ps,   0,  28,   2,  33,  16,  12, inv_r[2][46], inv_r[2][47]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  56,   2,   5,  48,   8, inv_r[2][48], inv_r[2][49]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  44,   2,  17,  12,  32, inv_r[2][50], inv_r[2][51]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  28,   2,  33,  12,  16, inv_r[2][52], inv_r[2][53]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  52,   2,   9,  48,   4, inv_r[2][54], inv_r[2][55]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   1,   2,  12,  48, inv_r[2][56], inv_r[2][57]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   1,   2,  48,  12, inv_r[2][58], inv_r[2][59]);} // 1 times

    {psi.inv_decay(no_ps,  12,   4,   8, inv_r[3][0], inv_r[3][1]);} // 15 times
    {psi.inv_decay(no_ps,  20,   4,  16, inv_r[3][2], inv_r[3][3]);} // 5 times
    {psi.inv_decay(no_ps,  24,   8,  16, inv_r[3][4], inv_r[3][5]);} // 5 times
    {psi.inv_decay(no_ps,  28,   4,  24, inv_r[3][6], inv_r[3][7]);} // 3 times
    {psi.inv_decay(no_ps,  44,  32,  12, inv_r[3][8], inv_r[3][9]);} // 3 times
    {psi.inv_decay(no_ps,  28,  16,  12, inv_r[3][10], inv_r[3][11]);} // 3 times
    {psi.inv_decay(no_ps,  28,   8,  20, inv_r[3][12], inv_r[3][13]);} // 3 times
    {psi.inv_decay(no_ps,  60,  32,  28, inv_r[3][14], inv_r[3][15]);} // 3 times
    {psi.inv_decay(no_ps,  60,  16,  44, inv_r[3][16], inv_r[3][17]);} // 1 times
    {psi.inv_decay(no_ps,  60,  12,  48, inv_r[3][18], inv_r[3][19]);} // 1 times
    {psi.inv_decay(no_ps,  48,  16,  32, inv_r[3][20], inv_r[3][21]);} // 5 times

    if (psi_weight_IS == 2){
      int n_random = 3 * psi_n_particle - 4;;
      vector<vector<double> > inv_channel_r(35, vector<double> (8));
      vector<double> g_IS(35, 1.);
      for (int i_c = 0; i_c < 35; i_c++){
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
      if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = g_IS[  1] * psi_g_f[no_ps][  2] * psi_g_f[no_ps][  3] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  3] * psi_g_t[no_ps][  4];}
      if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = g_IS[  2] * psi_g_f[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_t[no_ps][  5];}
      if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = g_IS[  3] * psi_g_f[no_ps][  4] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_t[no_ps][  8];}
      if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = g_IS[  4] * psi_g_f[no_ps][  6] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  9] * psi_g_t[no_ps][ 10] * psi_g_t[no_ps][ 11];}
      if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = g_IS[  5] * psi_g_f[no_ps][  8] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][  9] * psi_g_t[no_ps][ 12] * psi_g_t[no_ps][ 13];}
      if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = g_IS[  6] * psi_g_f[no_ps][  6] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  9] * psi_g_t[no_ps][ 10] * psi_g_t[no_ps][ 14];}
      if (psi_MC_alpha[zero +   7] != 0.){psi_MC_g_channel[zero +   7] = g_IS[  7] * psi_g_f[no_ps][  4] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_t[no_ps][ 15];}
      if (psi_MC_alpha[zero +   8] != 0.){psi_MC_g_channel[zero +   8] = g_IS[  8] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 16] * psi_g_t[no_ps][ 17] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +   9] != 0.){psi_MC_g_channel[zero +   9] = g_IS[  9] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][ 18] * psi_g_t[no_ps][ 19] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  10] != 0.){psi_MC_g_channel[zero +  10] = g_IS[ 10] * psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 11] * psi_g_t[no_ps][ 16] * psi_g_t[no_ps][ 20] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  11] != 0.){psi_MC_g_channel[zero +  11] = g_IS[ 11] * psi_g_p[no_ps][  2] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 16] * psi_g_t[no_ps][ 21] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  12] != 0.){psi_MC_g_channel[zero +  12] = g_IS[ 12] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][  3] * psi_g_t[no_ps][ 16] * psi_g_d[no_ps][  3] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  13] != 0.){psi_MC_g_channel[zero +  13] = g_IS[ 13] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  4] * psi_g_t[no_ps][ 18] * psi_g_d[no_ps][  4] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  14] != 0.){psi_MC_g_channel[zero +  14] = g_IS[ 14] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  5] * psi_g_t[no_ps][ 16] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  15] != 0.){psi_MC_g_channel[zero +  15] = g_IS[ 15] * psi_g_p[no_ps][  1] * psi_g_p[no_ps][  6] * psi_g_t[no_ps][ 16] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  16] != 0.){psi_MC_g_channel[zero +  16] = g_IS[ 16] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][  3] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  3] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  17] != 0.){psi_MC_g_channel[zero +  17] = g_IS[ 17] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  4] * psi_g_d[no_ps][  8] * psi_g_d[no_ps][  4] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  18] != 0.){psi_MC_g_channel[zero +  18] = g_IS[ 18] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  5] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  19] != 0.){psi_MC_g_channel[zero +  19] = g_IS[ 19] * psi_g_p[no_ps][  1] * psi_g_p[no_ps][  6] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  20] != 0.){psi_MC_g_channel[zero +  20] = g_IS[ 20] * psi_g_p[no_ps][  7] * psi_g_p[no_ps][  8] * psi_g_d[no_ps][  9] * psi_g_d[no_ps][  0] * psi_g_d[no_ps][ 10];}
      if (psi_MC_alpha[zero +  21] != 0.){psi_MC_g_channel[zero +  21] = g_IS[ 21] * psi_g_p[no_ps][  2] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  22] != 0.){psi_MC_g_channel[zero +  22] = g_IS[ 22] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 23] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  23] != 0.){psi_MC_g_channel[zero +  23] = g_IS[ 23] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  24] != 0.){psi_MC_g_channel[zero +  24] = g_IS[ 24] * psi_g_p[no_ps][  1] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  9] * psi_g_t[no_ps][ 10] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  25] != 0.){psi_MC_g_channel[zero +  25] = g_IS[ 25] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][  3] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][  3] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  26] != 0.){psi_MC_g_channel[zero +  26] = g_IS[ 26] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  4] * psi_g_t[no_ps][  6] * psi_g_d[no_ps][  4] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  27] != 0.){psi_MC_g_channel[zero +  27] = g_IS[ 27] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  5] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  28] != 0.){psi_MC_g_channel[zero +  28] = g_IS[ 28] * psi_g_p[no_ps][  1] * psi_g_p[no_ps][  6] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  29] != 0.){psi_MC_g_channel[zero +  29] = g_IS[ 29] * psi_g_p[no_ps][  7] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 24] * psi_g_d[no_ps][ 10];}
      if (psi_MC_alpha[zero +  30] != 0.){psi_MC_g_channel[zero +  30] = g_IS[ 30] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 25] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  31] != 0.){psi_MC_g_channel[zero +  31] = g_IS[ 31] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 26] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  32] != 0.){psi_MC_g_channel[zero +  32] = g_IS[ 32] * psi_g_p[no_ps][  7] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][  9] * psi_g_t[no_ps][ 27] * psi_g_d[no_ps][ 10];}
      if (psi_MC_alpha[zero +  33] != 0.){psi_MC_g_channel[zero +  33] = g_IS[ 33] * psi_g_p[no_ps][  7] * psi_g_p[no_ps][  8] * psi_g_t[no_ps][ 28] * psi_g_d[no_ps][  0] * psi_g_d[no_ps][ 10];}
      if (psi_MC_alpha[zero +  34] != 0.){psi_MC_g_channel[zero +  34] = g_IS[ 34] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  9] * psi_g_t[no_ps][ 29] * psi_g_d[no_ps][ 10] * psi_g_d[no_ps][  0];}

    }
    else if (psi_weight_IS == 4){
      vector<vector<double> > g_IS(4);
      g_IS[0].resize(10, 1.);
      g_IS[1].resize(15, 1.);
      g_IS[2].resize(30, 1.);
      g_IS[3].resize(11, 1.);

      for (int i_m = 0; i_m < 4; i_m++){
        for (int i_r = 0; i_r < inv_r[i_m].size(); i_r++){
          psi_phasespace_randoms[psi_container_IS_startvalue[no_ps][i_m] + i_r]->get_g_IS(inv_r[i_m][i_r], g_IS_temp);
          int i_g = i_r;
          if (i_m > 1){i_g = i_r / 2;}
          g_IS[i_m][i_g] *= g_IS_temp;
        }
      }
      if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = psi_g_f[no_ps][  0] * g_IS[1][  0] * psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  1] * g_IS[2][  1] * psi_g_t[no_ps][  2] * g_IS[2][  2];}
      if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = psi_g_f[no_ps][  2] * g_IS[1][  2] * psi_g_f[no_ps][  3] * g_IS[1][  3] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  3] * g_IS[2][  3] * psi_g_t[no_ps][  4] * g_IS[2][  4];}
      if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = psi_g_f[no_ps][  0] * g_IS[1][  0] * psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  1] * g_IS[2][  1] * psi_g_t[no_ps][  5] * g_IS[2][  5];}
      if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = psi_g_f[no_ps][  4] * g_IS[1][  4] * psi_g_f[no_ps][  5] * g_IS[1][  5] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][  7] * g_IS[2][  7] * psi_g_t[no_ps][  8] * g_IS[2][  8];}
      if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = psi_g_f[no_ps][  6] * g_IS[1][  6] * psi_g_f[no_ps][  7] * g_IS[1][  7] * psi_g_t[no_ps][  9] * g_IS[2][  9] * psi_g_t[no_ps][ 10] * g_IS[2][ 10] * psi_g_t[no_ps][ 11] * g_IS[2][ 11];}
      if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = psi_g_f[no_ps][  8] * g_IS[1][  8] * psi_g_f[no_ps][  9] * g_IS[1][  9] * psi_g_t[no_ps][  9] * g_IS[2][  9] * psi_g_t[no_ps][ 12] * g_IS[2][ 12] * psi_g_t[no_ps][ 13] * g_IS[2][ 13];}
      if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = psi_g_f[no_ps][  6] * g_IS[1][  6] * psi_g_f[no_ps][  7] * g_IS[1][  7] * psi_g_t[no_ps][  9] * g_IS[2][  9] * psi_g_t[no_ps][ 10] * g_IS[2][ 10] * psi_g_t[no_ps][ 14] * g_IS[2][ 14];}
      if (psi_MC_alpha[zero +   7] != 0.){psi_MC_g_channel[zero +   7] = psi_g_f[no_ps][  4] * g_IS[1][  4] * psi_g_f[no_ps][  5] * g_IS[1][  5] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][  7] * g_IS[2][  7] * psi_g_t[no_ps][ 15] * g_IS[2][ 15];}
      if (psi_MC_alpha[zero +   8] != 0.){psi_MC_g_channel[zero +   8] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][ 10] * g_IS[1][ 10] * psi_g_t[no_ps][ 16] * g_IS[2][ 16] * psi_g_t[no_ps][ 17] * g_IS[2][ 17] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +   9] != 0.){psi_MC_g_channel[zero +   9] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][  5] * g_IS[1][  5] * psi_g_t[no_ps][ 18] * g_IS[2][ 18] * psi_g_t[no_ps][ 19] * g_IS[2][ 19] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  10] != 0.){psi_MC_g_channel[zero +  10] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_f[no_ps][ 11] * g_IS[1][ 11] * psi_g_t[no_ps][ 16] * g_IS[2][ 16] * psi_g_t[no_ps][ 20] * g_IS[2][ 20] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  11] != 0.){psi_MC_g_channel[zero +  11] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_f[no_ps][ 12] * g_IS[1][ 12] * psi_g_t[no_ps][ 16] * g_IS[2][ 16] * psi_g_t[no_ps][ 21] * g_IS[2][ 21] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  12] != 0.){psi_MC_g_channel[zero +  12] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_t[no_ps][ 16] * g_IS[2][ 16] * psi_g_d[no_ps][  3] * g_IS[3][  3] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  13] != 0.){psi_MC_g_channel[zero +  13] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_t[no_ps][ 18] * g_IS[2][ 18] * psi_g_d[no_ps][  4] * g_IS[3][  4] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  14] != 0.){psi_MC_g_channel[zero +  14] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_t[no_ps][ 16] * g_IS[2][ 16] * psi_g_d[no_ps][  5] * g_IS[3][  5] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  15] != 0.){psi_MC_g_channel[zero +  15] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_p[no_ps][  6] * g_IS[0][  6] * psi_g_t[no_ps][ 16] * g_IS[2][ 16] * psi_g_d[no_ps][  6] * g_IS[3][  6] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  16] != 0.){psi_MC_g_channel[zero +  16] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_d[no_ps][  7] * g_IS[3][  7] * psi_g_d[no_ps][  3] * g_IS[3][  3] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  17] != 0.){psi_MC_g_channel[zero +  17] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_d[no_ps][  8] * g_IS[3][  8] * psi_g_d[no_ps][  4] * g_IS[3][  4] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  18] != 0.){psi_MC_g_channel[zero +  18] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_d[no_ps][  7] * g_IS[3][  7] * psi_g_d[no_ps][  5] * g_IS[3][  5] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  19] != 0.){psi_MC_g_channel[zero +  19] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_p[no_ps][  6] * g_IS[0][  6] * psi_g_d[no_ps][  7] * g_IS[3][  7] * psi_g_d[no_ps][  6] * g_IS[3][  6] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  20] != 0.){psi_MC_g_channel[zero +  20] = psi_g_p[no_ps][  7] * g_IS[0][  7] * psi_g_p[no_ps][  8] * g_IS[0][  8] * psi_g_d[no_ps][  9] * g_IS[3][  9] * psi_g_d[no_ps][  0] * g_IS[3][  0] * psi_g_d[no_ps][ 10] * g_IS[3][ 10];}
      if (psi_MC_alpha[zero +  21] != 0.){psi_MC_g_channel[zero +  21] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  1] * g_IS[2][  1] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  22] != 0.){psi_MC_g_channel[zero +  22] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][ 10] * g_IS[1][ 10] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 23] * g_IS[2][ 23] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  23] != 0.){psi_MC_g_channel[zero +  23] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][  5] * g_IS[1][  5] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][  7] * g_IS[2][  7] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  24] != 0.){psi_MC_g_channel[zero +  24] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_f[no_ps][  7] * g_IS[1][  7] * psi_g_t[no_ps][  9] * g_IS[2][  9] * psi_g_t[no_ps][ 10] * g_IS[2][ 10] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  25] != 0.){psi_MC_g_channel[zero +  25] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_d[no_ps][  3] * g_IS[3][  3] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  26] != 0.){psi_MC_g_channel[zero +  26] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_d[no_ps][  4] * g_IS[3][  4] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  27] != 0.){psi_MC_g_channel[zero +  27] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_d[no_ps][  5] * g_IS[3][  5] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  28] != 0.){psi_MC_g_channel[zero +  28] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_p[no_ps][  6] * g_IS[0][  6] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_d[no_ps][  6] * g_IS[3][  6] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  29] != 0.){psi_MC_g_channel[zero +  29] = psi_g_p[no_ps][  7] * g_IS[0][  7] * psi_g_f[no_ps][ 13] * g_IS[1][ 13] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][ 24] * g_IS[2][ 24] * psi_g_d[no_ps][ 10] * g_IS[3][ 10];}
      if (psi_MC_alpha[zero +  30] != 0.){psi_MC_g_channel[zero +  30] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][  5] * g_IS[1][  5] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][ 25] * g_IS[2][ 25] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  31] != 0.){psi_MC_g_channel[zero +  31] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][ 10] * g_IS[1][ 10] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 26] * g_IS[2][ 26] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  32] != 0.){psi_MC_g_channel[zero +  32] = psi_g_p[no_ps][  7] * g_IS[0][  7] * psi_g_f[no_ps][ 14] * g_IS[1][ 14] * psi_g_t[no_ps][  9] * g_IS[2][  9] * psi_g_t[no_ps][ 27] * g_IS[2][ 27] * psi_g_d[no_ps][ 10] * g_IS[3][ 10];}
      if (psi_MC_alpha[zero +  33] != 0.){psi_MC_g_channel[zero +  33] = psi_g_p[no_ps][  7] * g_IS[0][  7] * psi_g_p[no_ps][  8] * g_IS[0][  8] * psi_g_t[no_ps][ 28] * g_IS[2][ 28] * psi_g_d[no_ps][  0] * g_IS[3][  0] * psi_g_d[no_ps][ 10] * g_IS[3][ 10];}
      if (psi_MC_alpha[zero +  34] != 0.){psi_MC_g_channel[zero +  34] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  9] * g_IS[0][  9] * psi_g_t[no_ps][ 29] * g_IS[2][ 29] * psi_g_d[no_ps][ 10] * g_IS[3][ 10] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
    }
  }

  else {
    if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = psi_g_f[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_t[no_ps][  2];}
    if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = psi_g_f[no_ps][  2] * psi_g_f[no_ps][  3] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  3] * psi_g_t[no_ps][  4];}
    if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = psi_g_f[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_t[no_ps][  5];}
    if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = psi_g_f[no_ps][  4] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_t[no_ps][  8];}
    if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = psi_g_f[no_ps][  6] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  9] * psi_g_t[no_ps][ 10] * psi_g_t[no_ps][ 11];}
    if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = psi_g_f[no_ps][  8] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][  9] * psi_g_t[no_ps][ 12] * psi_g_t[no_ps][ 13];}
    if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = psi_g_f[no_ps][  6] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  9] * psi_g_t[no_ps][ 10] * psi_g_t[no_ps][ 14];}
    if (psi_MC_alpha[zero +   7] != 0.){psi_MC_g_channel[zero +   7] = psi_g_f[no_ps][  4] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_t[no_ps][ 15];}
    if (psi_MC_alpha[zero +   8] != 0.){psi_MC_g_channel[zero +   8] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 16] * psi_g_t[no_ps][ 17] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +   9] != 0.){psi_MC_g_channel[zero +   9] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][ 18] * psi_g_t[no_ps][ 19] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  10] != 0.){psi_MC_g_channel[zero +  10] = psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 11] * psi_g_t[no_ps][ 16] * psi_g_t[no_ps][ 20] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  11] != 0.){psi_MC_g_channel[zero +  11] = psi_g_p[no_ps][  2] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 16] * psi_g_t[no_ps][ 21] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  12] != 0.){psi_MC_g_channel[zero +  12] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][  3] * psi_g_t[no_ps][ 16] * psi_g_d[no_ps][  3] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  13] != 0.){psi_MC_g_channel[zero +  13] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  4] * psi_g_t[no_ps][ 18] * psi_g_d[no_ps][  4] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  14] != 0.){psi_MC_g_channel[zero +  14] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  5] * psi_g_t[no_ps][ 16] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  15] != 0.){psi_MC_g_channel[zero +  15] = psi_g_p[no_ps][  1] * psi_g_p[no_ps][  6] * psi_g_t[no_ps][ 16] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  16] != 0.){psi_MC_g_channel[zero +  16] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][  3] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  3] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  17] != 0.){psi_MC_g_channel[zero +  17] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  4] * psi_g_d[no_ps][  8] * psi_g_d[no_ps][  4] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  18] != 0.){psi_MC_g_channel[zero +  18] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  5] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  19] != 0.){psi_MC_g_channel[zero +  19] = psi_g_p[no_ps][  1] * psi_g_p[no_ps][  6] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  20] != 0.){psi_MC_g_channel[zero +  20] = psi_g_p[no_ps][  7] * psi_g_p[no_ps][  8] * psi_g_d[no_ps][  9] * psi_g_d[no_ps][  0] * psi_g_d[no_ps][ 10];}
    if (psi_MC_alpha[zero +  21] != 0.){psi_MC_g_channel[zero +  21] = psi_g_p[no_ps][  2] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  22] != 0.){psi_MC_g_channel[zero +  22] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 23] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  23] != 0.){psi_MC_g_channel[zero +  23] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  24] != 0.){psi_MC_g_channel[zero +  24] = psi_g_p[no_ps][  1] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  9] * psi_g_t[no_ps][ 10] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  25] != 0.){psi_MC_g_channel[zero +  25] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][  3] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][  3] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  26] != 0.){psi_MC_g_channel[zero +  26] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  4] * psi_g_t[no_ps][  6] * psi_g_d[no_ps][  4] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  27] != 0.){psi_MC_g_channel[zero +  27] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  5] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  28] != 0.){psi_MC_g_channel[zero +  28] = psi_g_p[no_ps][  1] * psi_g_p[no_ps][  6] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  29] != 0.){psi_MC_g_channel[zero +  29] = psi_g_p[no_ps][  7] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 24] * psi_g_d[no_ps][ 10];}
    if (psi_MC_alpha[zero +  30] != 0.){psi_MC_g_channel[zero +  30] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 25] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  31] != 0.){psi_MC_g_channel[zero +  31] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 26] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  32] != 0.){psi_MC_g_channel[zero +  32] = psi_g_p[no_ps][  7] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][  9] * psi_g_t[no_ps][ 27] * psi_g_d[no_ps][ 10];}
    if (psi_MC_alpha[zero +  33] != 0.){psi_MC_g_channel[zero +  33] = psi_g_p[no_ps][  7] * psi_g_p[no_ps][  8] * psi_g_t[no_ps][ 28] * psi_g_d[no_ps][  0] * psi_g_d[no_ps][ 10];}
    if (psi_MC_alpha[zero +  34] != 0.){psi_MC_g_channel[zero +  34] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  9] * psi_g_t[no_ps][ 29] * psi_g_d[no_ps][ 10] * psi_g_d[no_ps][  0];}
  }

/*
  if (psi_MC_alpha[   0] != 0.){ag_channel_1111(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[   1] != 0.){ag_channel_1111(g, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[   2] != 0.){ag_channel_1111(g, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[   3] != 0.){ag_channel_1111(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  if (psi_MC_alpha[   4] != 0.){ag_channel_1111(g, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[   5] != 0.){ag_channel_1111(g, 1, 2, 8, 4, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[   6] != 0.){ag_channel_1111(g, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[   7] != 0.){ag_channel_1111(g, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  if (psi_MC_alpha[   8] != 0.){ag_channel_2v11(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[   9] != 0.){ag_channel_2v11(g, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  10] != 0.){ag_channel_2v11(g, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  11] != 0.){ag_channel_2v11(g, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  12] != 0.){ag_channel_3yv1(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  13] != 0.){ag_channel_3yv1(g, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  14] != 0.){ag_channel_3yv1(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  15] != 0.){ag_channel_3yv1(g, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  16] != 0.){ag_channel_4yyv(g, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  17] != 0.){ag_channel_4yyv(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  18] != 0.){ag_channel_4yyv(g, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  19] != 0.){ag_channel_4yyv(g, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  20] != 0.){ag_channel_42vv(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  21] != 0.){ag_channel_12v1(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  22] != 0.){ag_channel_12v1(g, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  23] != 0.){ag_channel_12v1(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  24] != 0.){ag_channel_12v1(g, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  25] != 0.){ag_channel_13yv(g, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  26] != 0.){ag_channel_13yv(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  27] != 0.){ag_channel_13yv(g, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  28] != 0.){ag_channel_13yv(g, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  29] != 0.){ag_channel_112v(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  30] != 0.){ag_channel_112v(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  31] != 0.){ag_channel_112v(g, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  32] != 0.){ag_channel_112v(g, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  33] != 0.){ag_channel_2v2v(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  34] != 0.){ag_channel_2v2v(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  // channel    0: 1111(1, 2, 4, 8, 16, 32;   6,   0,   0)
  // channel    1: 1111(1, 2, 4, 8, 32, 16;   6,   0,   0)
  // channel    2: 1111(1, 2, 4, 16, 8, 32;   6,   0,   6)
  // channel    3: 1111(1, 2, 16, 4, 8, 32;   0,   0,   6)
  // channel    4: 1111(1, 2, 8, 4, 16, 32;   6,   0,   0)
  // channel    5: 1111(1, 2, 8, 4, 32, 16;   6,   0,   0)
  // channel    6: 1111(1, 2, 8, 16, 4, 32;   6,   0,   6)
  // channel    7: 1111(1, 2, 16, 8, 4, 32;   0,   0,   6)
  // channel    8: 2v11(1, 2, 4, 8, 16, 32;   0,   0,   0)
  // channel    9: 2v11(1, 2, 4, 8, 32, 16;   0,   0,   0)
  // channel   10: 2v11(1, 2, 4, 16, 8, 32;   6,   0,   6)
  // channel   11: 2v11(1, 2, 8, 16, 4, 32;   6,   0,   6)
  // channel   12: 3yv1(1, 2, 4, 8, 16, 32;   6,   0,   0)
  // channel   13: 3yv1(1, 2, 32, 4, 8, 16;   0,   0,   0)
  // channel   14: 3yv1(1, 2, 16, 4, 8, 32;   0,   0,   0)
  // channel   15: 3yv1(1, 2, 8, 4, 16, 32;   6,   0,   0)
  // channel   16: 4yyv(1, 2, 32, 4, 8, 16;   6,   0,   0)
  // channel   17: 4yyv(1, 2, 16, 32, 4, 8;   0,   0,   0)
  // channel   18: 4yyv(1, 2, 32, 16, 4, 8;   0,   0,   0)
  // channel   19: 4yyv(1, 2, 32, 8, 4, 16;   6,   0,   0)
  // channel   20: 42vv(1, 2, 4, 8, 16, 32;   0,   0,   0)
  // channel   21: 12v1(1, 2, 4, 8, 16, 32;   6,   6,   0)
  // channel   22: 12v1(1, 2, 32, 4, 8, 16;   0,   0,   0)
  // channel   23: 12v1(1, 2, 16, 4, 8, 32;   0,   0,   0)
  // channel   24: 12v1(1, 2, 8, 4, 16, 32;   6,   6,   0)
  // channel   25: 13yv(1, 2, 32, 4, 8, 16;   6,   0,   0)
  // channel   26: 13yv(1, 2, 16, 32, 4, 8;   0,   0,   0)
  // channel   27: 13yv(1, 2, 32, 16, 4, 8;   0,   0,   0)
  // channel   28: 13yv(1, 2, 32, 8, 4, 16;   6,   0,   0)
  // channel   29: 112v(1, 2, 4, 8, 16, 32;   0,   6,   0)
  // channel   30: 112v(1, 2, 16, 32, 4, 8;   0,   0,   0)
  // channel   31: 112v(1, 2, 32, 16, 4, 8;   0,   0,   0)
  // channel   32: 112v(1, 2, 8, 4, 16, 32;   0,   6,   0)
  // channel   33: 2v2v(1, 2, 4, 8, 16, 32;   0,   0,   0)
  // channel   34: 2v2v(1, 2, 16, 32, 4, 8;   0,   0,   0)
*/

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
