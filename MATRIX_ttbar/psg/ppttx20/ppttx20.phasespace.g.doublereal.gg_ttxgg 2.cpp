#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ag_psp_400_gg_ttxgg(int no_ps, int zero, phasespace_set & psi){
  static Logger logger("ppttx20_ag_psp_400_gg_ttxgg");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (psi_xbp[no_ps][5] == nullvector){psi_xbp[no_ps][5] = psi_xbp[no_ps][1] - psi_xbp[no_ps][4];}
  if (psi_xbs[no_ps][5] == 0.){psi_xbs[no_ps][5] = psi_xbp[no_ps][5].m2();}
  if (psi_xbs[no_ps][5] > 0.){psi_xbs[no_ps][5] = 0.;}
  if (psi_xbp[no_ps][6] == nullvector){psi_xbp[no_ps][6] = psi_xbp[no_ps][2] - psi_xbp[no_ps][4];}
  if (psi_xbs[no_ps][6] == 0.){psi_xbs[no_ps][6] = psi_xbp[no_ps][6].m2();}
  if (psi_xbs[no_ps][6] > 0.){psi_xbs[no_ps][6] = 0.;}
  if (psi_xbp[no_ps][9] == nullvector){psi_xbp[no_ps][9] = psi_xbp[no_ps][1] - psi_xbp[no_ps][8];}
  if (psi_xbs[no_ps][9] == 0.){psi_xbs[no_ps][9] = psi_xbp[no_ps][9].m2();}
  if (psi_xbs[no_ps][9] > 0.){psi_xbs[no_ps][9] = 0.;}
  if (psi_xbp[no_ps][10] == nullvector){psi_xbp[no_ps][10] = psi_xbp[no_ps][2] - psi_xbp[no_ps][8];}
  if (psi_xbs[no_ps][10] == 0.){psi_xbs[no_ps][10] = psi_xbp[no_ps][10].m2();}
  if (psi_xbs[no_ps][10] > 0.){psi_xbs[no_ps][10] = 0.;}
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
  if (psi_xbp[no_ps][22] == nullvector){psi_xbp[no_ps][22] = psi_xbp[no_ps][2] - psi_xbp[no_ps][4] - psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][22] == 0.){psi_xbs[no_ps][22] = psi_xbp[no_ps][22].m2();}
  if (psi_xbs[no_ps][22] > 0.){psi_xbs[no_ps][22] = 0.;}
  if (psi_xbp[no_ps][25] == nullvector){psi_xbp[no_ps][25] = psi_xbp[no_ps][1] - psi_xbp[no_ps][8] - psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][25] == 0.){psi_xbs[no_ps][25] = psi_xbp[no_ps][25].m2();}
  if (psi_xbs[no_ps][25] > 0.){psi_xbs[no_ps][25] = 0.;}
  if (psi_xbp[no_ps][26] == nullvector){psi_xbp[no_ps][26] = psi_xbp[no_ps][2] - psi_xbp[no_ps][8] - psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][26] == 0.){psi_xbs[no_ps][26] = psi_xbp[no_ps][26].m2();}
  if (psi_xbs[no_ps][26] > 0.){psi_xbs[no_ps][26] = 0.;}
  if (psi_xbp[no_ps][28] == nullvector){psi_xbp[no_ps][28] = psi_xbp[no_ps][4] + psi_xbp[no_ps][8] + psi_xbp[no_ps][16];}
  if (psi_xbs[no_ps][28] == 0.){psi_xbs[no_ps][28] = psi_xbp[no_ps][28].m2();}
  if (psi_xbs[no_ps][28] < 0.){psi_xbs[no_ps][28] = 0.;}
  if (psi_xbp[no_ps][37] == nullvector){psi_xbp[no_ps][37] = psi_xbp[no_ps][1] - psi_xbp[no_ps][4] - psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][37] == 0.){psi_xbs[no_ps][37] = psi_xbp[no_ps][37].m2();}
  if (psi_xbs[no_ps][37] > 0.){psi_xbs[no_ps][37] = 0.;}
  if (psi_xbp[no_ps][38] == nullvector){psi_xbp[no_ps][38] = psi_xbp[no_ps][2] - psi_xbp[no_ps][4] - psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][38] == 0.){psi_xbs[no_ps][38] = psi_xbp[no_ps][38].m2();}
  if (psi_xbs[no_ps][38] > 0.){psi_xbs[no_ps][38] = 0.;}
  if (psi_xbp[no_ps][41] == nullvector){psi_xbp[no_ps][41] = psi_xbp[no_ps][1] - psi_xbp[no_ps][8] - psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][41] == 0.){psi_xbs[no_ps][41] = psi_xbp[no_ps][41].m2();}
  if (psi_xbs[no_ps][41] > 0.){psi_xbs[no_ps][41] = 0.;}
  if (psi_xbp[no_ps][42] == nullvector){psi_xbp[no_ps][42] = psi_xbp[no_ps][2] - psi_xbp[no_ps][8] - psi_xbp[no_ps][32];}
  if (psi_xbs[no_ps][42] == 0.){psi_xbs[no_ps][42] = psi_xbp[no_ps][42].m2();}
  if (psi_xbs[no_ps][42] > 0.){psi_xbs[no_ps][42] = 0.;}
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

  psi_v_smin[no_ps][  0] = psi_smin_opt[no_ps][12]; // 19 times
  psi_v_smin[no_ps][  1] = psi_smin_opt[no_ps][36]; // 19 times
  psi_v_smin[no_ps][  2] = psi_smin_opt[no_ps][20]; // 19 times
  psi_v_smin[no_ps][  3] = psi_smin_opt[no_ps][48]; // 19 times
  psi_v_smin[no_ps][  4] = psi_smin_opt[no_ps][40]; // 19 times
  psi_v_smin[no_ps][  5] = psi_smin_opt[no_ps][24]; // 19 times
  psi_v_smin[no_ps][  6] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][24], 2), psi_smin_opt[no_ps][28]); // 8 times
  psi_v_smin[no_ps][  7] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][40], 2), psi_smin_opt[no_ps][44]); // 8 times
  psi_v_smin[no_ps][  8] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][44]); // 8 times
  psi_v_smin[no_ps][  9] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][4] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][52]); // 8 times
  psi_v_smin[no_ps][ 10] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][20], 2), psi_smin_opt[no_ps][52]); // 8 times
  psi_v_smin[no_ps][ 11] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][12], 2), psi_smin_opt[no_ps][28]); // 8 times
  psi_v_smin[no_ps][ 12] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][36], 2), psi_smin_opt[no_ps][52]); // 8 times
  psi_v_smin[no_ps][ 13] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][20], 2), psi_smin_opt[no_ps][28]); // 8 times
  psi_v_smin[no_ps][ 14] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][36], 2), psi_smin_opt[no_ps][44]); // 8 times
  psi_v_smin[no_ps][ 15] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][8] + psi_xbsqrts[no_ps][48], 2), psi_smin_opt[no_ps][56]); // 8 times
  psi_v_smin[no_ps][ 16] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][32] + psi_xbsqrts[no_ps][24], 2), psi_smin_opt[no_ps][56]); // 8 times
  psi_v_smin[no_ps][ 17] = GSL_MAX_DBL(pow(psi_sqrtsmin_opt[no_ps][16] + psi_xbsqrts[no_ps][40], 2), psi_smin_opt[no_ps][56]); // 8 times

  psi_v_smax[no_ps][  0] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][48], 2); // 17 times
  psi_v_smax[no_ps][  1] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][24], 2); // 17 times
  psi_v_smax[no_ps][  2] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][40], 2); // 17 times
  psi_v_smax[no_ps][  3] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][12], 2); // 18 times
  psi_v_smax[no_ps][  4] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][20], 2); // 18 times
  psi_v_smax[no_ps][  5] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][36], 2); // 18 times
  psi_v_smax[no_ps][  6] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][32], 2); // 24 times
  psi_v_smax[no_ps][  7] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][16], 2); // 24 times
  psi_v_smax[no_ps][  8] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][8], 2); // 24 times
  psi_v_smax[no_ps][  9] = pow(psi_xbsqrts[no_ps][0] - psi_sqrtsmin_opt[no_ps][4], 2); // 24 times
  psi_v_smax[no_ps][ 10] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][48], 2); // 2 times
  psi_v_smax[no_ps][ 11] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][24], 2); // 2 times
  psi_v_smax[no_ps][ 12] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][40], 2); // 2 times
  psi_v_smax[no_ps][ 13] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][12], 2); // 1 times
  psi_v_smax[no_ps][ 14] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][20], 2); // 1 times
  psi_v_smax[no_ps][ 15] = pow(psi_xbsqrts[no_ps][0] - psi_xbsqrts[no_ps][36], 2); // 1 times

  //  vector<double> g_p(24);
  {psi_g_p[no_ps][  0] = psi.g_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0]); /* npp 2 */} // 13 times
  {psi_g_p[no_ps][  1] = psi.g_propagator(no_ps,   6,  36, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1]); /* npp 2 */} // 13 times
  {psi_g_p[no_ps][  2] = psi.g_propagator(no_ps,   6,  20, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2]); /* npp 2 */} // 13 times
  {psi_g_p[no_ps][  3] = psi.g_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3]); /* npp 2 */} // 14 times
  {psi_g_p[no_ps][  4] = psi.g_propagator(no_ps,   6,  40, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4]); /* npp 2 */} // 14 times
  {psi_g_p[no_ps][  5] = psi.g_propagator(no_ps,   6,  24, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  5]); /* npp 2 */} // 14 times
  {psi_g_p[no_ps][  6] = psi.g_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  6]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][  7] = psi.g_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  7], psi_v_smax[no_ps][  7]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][  8] = psi.g_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  8], psi_v_smax[no_ps][  7]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][  9] = psi.g_propagator(no_ps,   6,  52, psi_v_smin[no_ps][  9], psi_v_smax[no_ps][  8]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][ 10] = psi.g_propagator(no_ps,   6,  52, psi_v_smin[no_ps][ 10], psi_v_smax[no_ps][  8]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][ 11] = psi.g_propagator(no_ps,   0,  28, psi_v_smin[no_ps][ 11], psi_v_smax[no_ps][  6]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][ 12] = psi.g_propagator(no_ps,   6,  52, psi_v_smin[no_ps][ 12], psi_v_smax[no_ps][  8]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][ 13] = psi.g_propagator(no_ps,   0,  28, psi_v_smin[no_ps][ 13], psi_v_smax[no_ps][  6]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][ 14] = psi.g_propagator(no_ps,   0,  44, psi_v_smin[no_ps][ 14], psi_v_smax[no_ps][  7]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][ 15] = psi.g_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 15], psi_v_smax[no_ps][  9]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][ 16] = psi.g_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 16], psi_v_smax[no_ps][  9]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][ 17] = psi.g_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 17], psi_v_smax[no_ps][  9]); /* npp 3 */} // 3 times
  {psi_g_p[no_ps][ 18] = psi.g_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][ 10]); /* npp 2 */} // 2 times
  {psi_g_p[no_ps][ 19] = psi.g_propagator(no_ps,   6,  36, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][ 11]); /* npp 2 */} // 2 times
  {psi_g_p[no_ps][ 20] = psi.g_propagator(no_ps,   6,  20, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][ 12]); /* npp 2 */} // 2 times
  {psi_g_p[no_ps][ 21] = psi.g_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][ 13]); /* npp 2 */} // 1 times
  {psi_g_p[no_ps][ 22] = psi.g_propagator(no_ps,   6,  40, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][ 14]); /* npp 2 */} // 1 times
  {psi_g_p[no_ps][ 23] = psi.g_propagator(no_ps,   6,  24, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][ 15]); /* npp 2 */} // 1 times

  //  vector<double> g_f(18);
  {psi_g_f[no_ps][  0] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  5]); /* nfp 2 */} // 4 times
  {psi_g_f[no_ps][  1] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 16], psi_v_smax[no_ps][  9]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][  2] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4]); /* nfp 2 */} // 4 times
  {psi_g_f[no_ps][  3] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 17], psi_v_smax[no_ps][  9]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][  4] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0]); /* nfp 2 */} // 4 times
  {psi_g_f[no_ps][  5] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 11], psi_v_smax[no_ps][  6]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][  6] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3]); /* nfp 2 */} // 4 times
  {psi_g_f[no_ps][  7] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 15], psi_v_smax[no_ps][  9]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][  8] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2]); /* nfp 2 */} // 4 times
  {psi_g_f[no_ps][  9] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 13], psi_v_smax[no_ps][  6]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][ 10] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  8], psi_v_smax[no_ps][  7]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][ 11] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1]); /* nfp 2 */} // 4 times
  {psi_g_f[no_ps][ 12] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 14], psi_v_smax[no_ps][  7]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][ 13] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 10], psi_v_smax[no_ps][  8]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][ 14] = psi.g_timelikeinvariant(psi_v_smin[no_ps][ 12], psi_v_smax[no_ps][  8]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][ 15] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  9], psi_v_smax[no_ps][  8]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][ 16] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  6]); /* nfp 3 */} // 5 times
  {psi_g_f[no_ps][ 17] = psi.g_timelikeinvariant(psi_v_smin[no_ps][  7], psi_v_smax[no_ps][  7]); /* nfp 3 */} // 5 times

  //  vector<double> g_t(74);
  {psi_g_t[no_ps][  0] = psi.g_tchannel_opt(no_ps,   6,  60,   1,   2,   4,  56); /* ntp 0 */} // 15 times
  {psi_g_t[no_ps][  1] = psi.g_tchannel_opt(no_ps,   0,  56,   2,   5,  32,  24); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][  2] = psi.g_tchannel(no_ps,   0,  24,   5,  34,   8,  16); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][  3] = psi.g_tchannel_opt(no_ps,   0,  56,   2,   5,  16,  40); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][  4] = psi.g_tchannel(no_ps,   0,  40,   5,  18,   8,  32); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][  5] = psi.g_tchannel(no_ps,   6,  40,   5,  18,  32,   8); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][  6] = psi.g_tchannel_opt(no_ps,   0,  60,   1,   2,  32,  28); /* ntp 0 */} // 15 times
  {psi_g_t[no_ps][  7] = psi.g_tchannel_opt(no_ps,   0,  28,   2,  33,  16,  12); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][  8] = psi.g_tchannel(no_ps,   6,  12,  33,  18,   4,   8); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][  9] = psi.g_tchannel(no_ps,   6,  24,   5,  34,  16,   8); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 10] = psi.g_tchannel_opt(no_ps,   6,  56,   2,   5,   8,  48); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][ 11] = psi.g_tchannel(no_ps,   6,  48,   5,  10,  16,  32); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 12] = psi.g_tchannel(no_ps,   6,  48,   5,  10,  32,  16); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 13] = psi.g_tchannel_opt(no_ps,   6,  28,   2,  33,   8,  20); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][ 14] = psi.g_tchannel(no_ps,   6,  20,  33,  10,   4,  16); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 15] = psi.g_tchannel_opt(no_ps,   0,  60,   1,   2,  16,  44); /* ntp 0 */} // 15 times
  {psi_g_t[no_ps][ 16] = psi.g_tchannel_opt(no_ps,   0,  44,   2,  17,  32,  12); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][ 17] = psi.g_tchannel(no_ps,   6,  12,  17,  34,   4,   8); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 18] = psi.g_tchannel_opt(no_ps,   6,  44,   2,  17,   8,  36); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][ 19] = psi.g_tchannel(no_ps,   6,  36,  17,  10,   4,  32); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 20] = psi.g_tchannel(no_ps,   0,  36,  17,  10,  32,   4); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 21] = psi.g_tchannel(no_ps,   0,  20,  33,  10,  16,   4); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 22] = psi.g_tchannel_opt(no_ps,   6,  60,   1,   2,   8,  52); /* ntp 0 */} // 15 times
  {psi_g_t[no_ps][ 23] = psi.g_tchannel_opt(no_ps,   0,  52,   2,   9,  32,  20); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][ 24] = psi.g_tchannel(no_ps,   0,  20,   9,  34,   4,  16); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 25] = psi.g_tchannel_opt(no_ps,   0,  52,   2,   9,  16,  36); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][ 26] = psi.g_tchannel(no_ps,   0,  36,   9,  18,   4,  32); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 27] = psi.g_tchannel(no_ps,   6,  36,   9,  18,  32,   4); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 28] = psi.g_tchannel(no_ps,   6,  12,  33,  18,   8,   4); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 29] = psi.g_tchannel(no_ps,   6,  20,   9,  34,  16,   4); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 30] = psi.g_tchannel_opt(no_ps,   6,  52,   2,   9,   4,  48); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][ 31] = psi.g_tchannel(no_ps,   6,  48,   9,   6,  16,  32); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 32] = psi.g_tchannel(no_ps,   6,  48,   9,   6,  32,  16); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 33] = psi.g_tchannel_opt(no_ps,   6,  28,   2,  33,   4,  24); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][ 34] = psi.g_tchannel(no_ps,   6,  24,  33,   6,   8,  16); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 35] = psi.g_tchannel(no_ps,   6,  12,  17,  34,   8,   4); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 36] = psi.g_tchannel_opt(no_ps,   6,  44,   2,  17,   4,  40); /* ntp 1 */} // 3 times
  {psi_g_t[no_ps][ 37] = psi.g_tchannel(no_ps,   6,  40,  17,   6,   8,  32); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 38] = psi.g_tchannel(no_ps,   0,  40,  17,   6,  32,   8); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 39] = psi.g_tchannel(no_ps,   0,  24,  33,   6,  16,   8); /* ntp 2 */} // 1 times
  {psi_g_t[no_ps][ 40] = psi.g_tchannel_opt(no_ps,   0,  60,   2,   1,  32,  28); /* ntp 0 */} // 6 times
  {psi_g_t[no_ps][ 41] = psi.g_tchannel_opt(no_ps,   0,  28,   1,  34,  12,  16); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 42] = psi.g_tchannel_opt(no_ps,   0,  60,   2,   1,  16,  44); /* ntp 0 */} // 6 times
  {psi_g_t[no_ps][ 43] = psi.g_tchannel_opt(no_ps,   0,  44,   1,  18,  12,  32); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 44] = psi.g_tchannel_opt(no_ps,   6,  44,   1,  18,  36,   8); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 45] = psi.g_tchannel_opt(no_ps,   6,  28,   1,  34,  20,   8); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 46] = psi.g_tchannel_opt(no_ps,   6,  60,   2,   1,   8,  52); /* ntp 0 */} // 6 times
  {psi_g_t[no_ps][ 47] = psi.g_tchannel_opt(no_ps,   6,  52,   1,  10,  20,  32); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 48] = psi.g_tchannel_opt(no_ps,   6,  52,   1,  10,  36,  16); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 49] = psi.g_tchannel_opt(no_ps,   0,  52,   1,  10,  48,   4); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 50] = psi.g_tchannel_opt(no_ps,   6,  44,   1,  18,  40,   4); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 51] = psi.g_tchannel_opt(no_ps,   6,  28,   1,  34,  24,   4); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 52] = psi.g_tchannel_opt(no_ps,   6,  60,   2,   1,   4,  56); /* ntp 0 */} // 6 times
  {psi_g_t[no_ps][ 53] = psi.g_tchannel_opt(no_ps,   6,  56,   1,   6,  24,  32); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 54] = psi.g_tchannel_opt(no_ps,   6,  56,   1,   6,  40,  16); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 55] = psi.g_tchannel_opt(no_ps,   0,  56,   1,   6,  48,   8); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 56] = psi.g_tchannel_opt(no_ps,   0,  56,   2,   5,  48,   8); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 57] = psi.g_tchannel_opt(no_ps,   6,  56,   2,   5,  24,  32); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 58] = psi.g_tchannel_opt(no_ps,   6,  28,   2,  33,  24,   4); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 59] = psi.g_tchannel_opt(no_ps,   6,  56,   2,   5,  40,  16); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 60] = psi.g_tchannel_opt(no_ps,   6,  44,   2,  17,  40,   4); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 61] = psi.g_tchannel_opt(no_ps,   0,  44,   2,  17,  12,  32); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 62] = psi.g_tchannel_opt(no_ps,   0,  28,   2,  33,  12,  16); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 63] = psi.g_tchannel_opt(no_ps,   0,  52,   2,   9,  48,   4); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 64] = psi.g_tchannel_opt(no_ps,   6,  52,   2,   9,  20,  32); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 65] = psi.g_tchannel_opt(no_ps,   6,  28,   2,  33,  20,   8); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 66] = psi.g_tchannel_opt(no_ps,   6,  52,   2,   9,  36,  16); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 67] = psi.g_tchannel_opt(no_ps,   6,  44,   2,  17,  36,   8); /* ntp 1 */} // 1 times
  {psi_g_t[no_ps][ 68] = psi.g_tchannel_opt(no_ps,   0,  60,   1,   2,  12,  48); /* ntp 0 */} // 1 times
  {psi_g_t[no_ps][ 69] = psi.g_tchannel_opt(no_ps,   6,  60,   1,   2,  36,  24); /* ntp 0 */} // 1 times
  {psi_g_t[no_ps][ 70] = psi.g_tchannel_opt(no_ps,   6,  60,   1,   2,  20,  40); /* ntp 0 */} // 1 times
  {psi_g_t[no_ps][ 71] = psi.g_tchannel_opt(no_ps,   0,  60,   1,   2,  48,  12); /* ntp 0 */} // 1 times
  {psi_g_t[no_ps][ 72] = psi.g_tchannel_opt(no_ps,   6,  60,   1,   2,  40,  20); /* ntp 0 */} // 1 times
  {psi_g_t[no_ps][ 73] = psi.g_tchannel_opt(no_ps,   6,  60,   1,   2,  24,  36); /* ntp 0 */} // 1 times

  //  vector<double> g_d(25);
  {psi_g_d[no_ps][  0] = psi.g_decay(no_ps,  12,   4,   8); /* ndp 2 */} // 15 times
  {psi_g_d[no_ps][  1] = psi.g_decay(no_ps,  36,   4,  32); /* ndp 2 */} // 15 times
  {psi_g_d[no_ps][  2] = psi.g_decay(no_ps,  20,   4,  16); /* ndp 2 */} // 15 times
  {psi_g_d[no_ps][  3] = psi.g_decay(no_ps,  48,  16,  32); /* ndp 2 */} // 15 times
  {psi_g_d[no_ps][  4] = psi.g_decay(no_ps,  40,   8,  32); /* ndp 2 */} // 15 times
  {psi_g_d[no_ps][  5] = psi.g_decay(no_ps,  24,   8,  16); /* ndp 2 */} // 15 times
  {psi_g_d[no_ps][  6] = psi.g_decay(no_ps,  28,   4,  24); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][  7] = psi.g_decay(no_ps,  44,   4,  40); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][  8] = psi.g_decay(no_ps,  44,  32,  12); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][  9] = psi.g_decay(no_ps,  52,   4,  48); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][ 10] = psi.g_decay(no_ps,  52,  32,  20); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][ 11] = psi.g_decay(no_ps,  28,  16,  12); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][ 12] = psi.g_decay(no_ps,  52,  16,  36); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][ 13] = psi.g_decay(no_ps,  28,   8,  20); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][ 14] = psi.g_decay(no_ps,  44,   8,  36); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][ 15] = psi.g_decay(no_ps,  56,   8,  48); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][ 16] = psi.g_decay(no_ps,  56,  32,  24); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][ 17] = psi.g_decay(no_ps,  56,  16,  40); /* ndp 3 */} // 3 times
  {psi_g_d[no_ps][ 18] = psi.g_decay(no_ps,  60,   4,  56); /* ndp 4 */} // 3 times
  {psi_g_d[no_ps][ 19] = psi.g_decay(no_ps,  60,  32,  28); /* ndp 4 */} // 3 times
  {psi_g_d[no_ps][ 20] = psi.g_decay(no_ps,  60,  16,  44); /* ndp 4 */} // 3 times
  {psi_g_d[no_ps][ 21] = psi.g_decay(no_ps,  60,   8,  52); /* ndp 4 */} // 3 times
  {psi_g_d[no_ps][ 22] = psi.g_decay(no_ps,  60,  12,  48); /* ndp 4 */} // 1 times
  {psi_g_d[no_ps][ 23] = psi.g_decay(no_ps,  60,  36,  24); /* ndp 4 */} // 1 times
  {psi_g_d[no_ps][ 24] = psi.g_decay(no_ps,  60,  20,  40); /* ndp 4 */} // 1 times


  if (psi_weight_IS == 2 || psi_weight_IS == 4){
    double g_IS_temp;
    vector<vector<double> > inv_r(4);
    inv_r[0].resize(24);
    inv_r[1].resize(18);
    inv_r[2].resize(148);
    inv_r[3].resize(50);

    {psi.inv_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], inv_r[0][0]);} // 13 times
    {psi.inv_propagator(no_ps,   6,  36, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], inv_r[0][1]);} // 13 times
    {psi.inv_propagator(no_ps,   6,  20, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], inv_r[0][2]);} // 13 times
    {psi.inv_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3], inv_r[0][3]);} // 14 times
    {psi.inv_propagator(no_ps,   6,  40, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4], inv_r[0][4]);} // 14 times
    {psi.inv_propagator(no_ps,   6,  24, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  5], inv_r[0][5]);} // 14 times
    {psi.inv_propagator(no_ps,   0,  28, psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  6], inv_r[0][6]);} // 3 times
    {psi.inv_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  7], psi_v_smax[no_ps][  7], inv_r[0][7]);} // 3 times
    {psi.inv_propagator(no_ps,   0,  44, psi_v_smin[no_ps][  8], psi_v_smax[no_ps][  7], inv_r[0][8]);} // 3 times
    {psi.inv_propagator(no_ps,   6,  52, psi_v_smin[no_ps][  9], psi_v_smax[no_ps][  8], inv_r[0][9]);} // 3 times
    {psi.inv_propagator(no_ps,   6,  52, psi_v_smin[no_ps][ 10], psi_v_smax[no_ps][  8], inv_r[0][10]);} // 3 times
    {psi.inv_propagator(no_ps,   0,  28, psi_v_smin[no_ps][ 11], psi_v_smax[no_ps][  6], inv_r[0][11]);} // 3 times
    {psi.inv_propagator(no_ps,   6,  52, psi_v_smin[no_ps][ 12], psi_v_smax[no_ps][  8], inv_r[0][12]);} // 3 times
    {psi.inv_propagator(no_ps,   0,  28, psi_v_smin[no_ps][ 13], psi_v_smax[no_ps][  6], inv_r[0][13]);} // 3 times
    {psi.inv_propagator(no_ps,   0,  44, psi_v_smin[no_ps][ 14], psi_v_smax[no_ps][  7], inv_r[0][14]);} // 3 times
    {psi.inv_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 15], psi_v_smax[no_ps][  9], inv_r[0][15]);} // 3 times
    {psi.inv_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 16], psi_v_smax[no_ps][  9], inv_r[0][16]);} // 3 times
    {psi.inv_propagator(no_ps,   6,  56, psi_v_smin[no_ps][ 17], psi_v_smax[no_ps][  9], inv_r[0][17]);} // 3 times
    {psi.inv_propagator(no_ps,   0,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][ 10], inv_r[0][18]);} // 2 times
    {psi.inv_propagator(no_ps,   6,  36, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][ 11], inv_r[0][19]);} // 2 times
    {psi.inv_propagator(no_ps,   6,  20, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][ 12], inv_r[0][20]);} // 2 times
    {psi.inv_propagator(no_ps,   0,  48, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][ 13], inv_r[0][21]);} // 1 times
    {psi.inv_propagator(no_ps,   6,  40, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][ 14], inv_r[0][22]);} // 1 times
    {psi.inv_propagator(no_ps,   6,  24, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][ 15], inv_r[0][23]);} // 1 times

    {psi.inv_timelikeinvariant(no_ps,  24, psi_v_smin[no_ps][  5], psi_v_smax[no_ps][  5], inv_r[1][0]);} // 4 times
    {psi.inv_timelikeinvariant(no_ps,  56, psi_v_smin[no_ps][ 16], psi_v_smax[no_ps][  9], inv_r[1][1]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  40, psi_v_smin[no_ps][  4], psi_v_smax[no_ps][  4], inv_r[1][2]);} // 4 times
    {psi.inv_timelikeinvariant(no_ps,  56, psi_v_smin[no_ps][ 17], psi_v_smax[no_ps][  9], inv_r[1][3]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  12, psi_v_smin[no_ps][  0], psi_v_smax[no_ps][  0], inv_r[1][4]);} // 4 times
    {psi.inv_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][ 11], psi_v_smax[no_ps][  6], inv_r[1][5]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  48, psi_v_smin[no_ps][  3], psi_v_smax[no_ps][  3], inv_r[1][6]);} // 4 times
    {psi.inv_timelikeinvariant(no_ps,  56, psi_v_smin[no_ps][ 15], psi_v_smax[no_ps][  9], inv_r[1][7]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  20, psi_v_smin[no_ps][  2], psi_v_smax[no_ps][  2], inv_r[1][8]);} // 4 times
    {psi.inv_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][ 13], psi_v_smax[no_ps][  6], inv_r[1][9]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  44, psi_v_smin[no_ps][  8], psi_v_smax[no_ps][  7], inv_r[1][10]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  36, psi_v_smin[no_ps][  1], psi_v_smax[no_ps][  1], inv_r[1][11]);} // 4 times
    {psi.inv_timelikeinvariant(no_ps,  44, psi_v_smin[no_ps][ 14], psi_v_smax[no_ps][  7], inv_r[1][12]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  52, psi_v_smin[no_ps][ 10], psi_v_smax[no_ps][  8], inv_r[1][13]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  52, psi_v_smin[no_ps][ 12], psi_v_smax[no_ps][  8], inv_r[1][14]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  52, psi_v_smin[no_ps][  9], psi_v_smax[no_ps][  8], inv_r[1][15]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  28, psi_v_smin[no_ps][  6], psi_v_smax[no_ps][  6], inv_r[1][16]);} // 5 times
    {psi.inv_timelikeinvariant(no_ps,  44, psi_v_smin[no_ps][  7], psi_v_smax[no_ps][  7], inv_r[1][17]);} // 5 times

    {psi.inv_tchannel_opt(no_ps,   6,  60,   1,   2,   4,  56, inv_r[2][0], inv_r[2][1]);} // 15 times
    {psi.inv_tchannel_opt(no_ps,   0,  56,   2,   5,  32,  24, inv_r[2][2], inv_r[2][3]);} // 3 times
    {psi.inv_tchannel(no_ps,   0,  24,   5,  34,   8,  16, inv_r[2][4], inv_r[2][5]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  56,   2,   5,  16,  40, inv_r[2][6], inv_r[2][7]);} // 3 times
    {psi.inv_tchannel(no_ps,   0,  40,   5,  18,   8,  32, inv_r[2][8], inv_r[2][9]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  40,   5,  18,  32,   8, inv_r[2][10], inv_r[2][11]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   1,   2,  32,  28, inv_r[2][12], inv_r[2][13]);} // 15 times
    {psi.inv_tchannel_opt(no_ps,   0,  28,   2,  33,  16,  12, inv_r[2][14], inv_r[2][15]);} // 3 times
    {psi.inv_tchannel(no_ps,   6,  12,  33,  18,   4,   8, inv_r[2][16], inv_r[2][17]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  24,   5,  34,  16,   8, inv_r[2][18], inv_r[2][19]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  56,   2,   5,   8,  48, inv_r[2][20], inv_r[2][21]);} // 3 times
    {psi.inv_tchannel(no_ps,   6,  48,   5,  10,  16,  32, inv_r[2][22], inv_r[2][23]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  48,   5,  10,  32,  16, inv_r[2][24], inv_r[2][25]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  28,   2,  33,   8,  20, inv_r[2][26], inv_r[2][27]);} // 3 times
    {psi.inv_tchannel(no_ps,   6,  20,  33,  10,   4,  16, inv_r[2][28], inv_r[2][29]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   1,   2,  16,  44, inv_r[2][30], inv_r[2][31]);} // 15 times
    {psi.inv_tchannel_opt(no_ps,   0,  44,   2,  17,  32,  12, inv_r[2][32], inv_r[2][33]);} // 3 times
    {psi.inv_tchannel(no_ps,   6,  12,  17,  34,   4,   8, inv_r[2][34], inv_r[2][35]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  44,   2,  17,   8,  36, inv_r[2][36], inv_r[2][37]);} // 3 times
    {psi.inv_tchannel(no_ps,   6,  36,  17,  10,   4,  32, inv_r[2][38], inv_r[2][39]);} // 1 times
    {psi.inv_tchannel(no_ps,   0,  36,  17,  10,  32,   4, inv_r[2][40], inv_r[2][41]);} // 1 times
    {psi.inv_tchannel(no_ps,   0,  20,  33,  10,  16,   4, inv_r[2][42], inv_r[2][43]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  60,   1,   2,   8,  52, inv_r[2][44], inv_r[2][45]);} // 15 times
    {psi.inv_tchannel_opt(no_ps,   0,  52,   2,   9,  32,  20, inv_r[2][46], inv_r[2][47]);} // 3 times
    {psi.inv_tchannel(no_ps,   0,  20,   9,  34,   4,  16, inv_r[2][48], inv_r[2][49]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  52,   2,   9,  16,  36, inv_r[2][50], inv_r[2][51]);} // 3 times
    {psi.inv_tchannel(no_ps,   0,  36,   9,  18,   4,  32, inv_r[2][52], inv_r[2][53]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  36,   9,  18,  32,   4, inv_r[2][54], inv_r[2][55]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  12,  33,  18,   8,   4, inv_r[2][56], inv_r[2][57]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  20,   9,  34,  16,   4, inv_r[2][58], inv_r[2][59]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  52,   2,   9,   4,  48, inv_r[2][60], inv_r[2][61]);} // 3 times
    {psi.inv_tchannel(no_ps,   6,  48,   9,   6,  16,  32, inv_r[2][62], inv_r[2][63]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  48,   9,   6,  32,  16, inv_r[2][64], inv_r[2][65]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  28,   2,  33,   4,  24, inv_r[2][66], inv_r[2][67]);} // 3 times
    {psi.inv_tchannel(no_ps,   6,  24,  33,   6,   8,  16, inv_r[2][68], inv_r[2][69]);} // 1 times
    {psi.inv_tchannel(no_ps,   6,  12,  17,  34,   8,   4, inv_r[2][70], inv_r[2][71]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  44,   2,  17,   4,  40, inv_r[2][72], inv_r[2][73]);} // 3 times
    {psi.inv_tchannel(no_ps,   6,  40,  17,   6,   8,  32, inv_r[2][74], inv_r[2][75]);} // 1 times
    {psi.inv_tchannel(no_ps,   0,  40,  17,   6,  32,   8, inv_r[2][76], inv_r[2][77]);} // 1 times
    {psi.inv_tchannel(no_ps,   0,  24,  33,   6,  16,   8, inv_r[2][78], inv_r[2][79]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   2,   1,  32,  28, inv_r[2][80], inv_r[2][81]);} // 6 times
    {psi.inv_tchannel_opt(no_ps,   0,  28,   1,  34,  12,  16, inv_r[2][82], inv_r[2][83]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   2,   1,  16,  44, inv_r[2][84], inv_r[2][85]);} // 6 times
    {psi.inv_tchannel_opt(no_ps,   0,  44,   1,  18,  12,  32, inv_r[2][86], inv_r[2][87]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  44,   1,  18,  36,   8, inv_r[2][88], inv_r[2][89]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  28,   1,  34,  20,   8, inv_r[2][90], inv_r[2][91]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  60,   2,   1,   8,  52, inv_r[2][92], inv_r[2][93]);} // 6 times
    {psi.inv_tchannel_opt(no_ps,   6,  52,   1,  10,  20,  32, inv_r[2][94], inv_r[2][95]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  52,   1,  10,  36,  16, inv_r[2][96], inv_r[2][97]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  52,   1,  10,  48,   4, inv_r[2][98], inv_r[2][99]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  44,   1,  18,  40,   4, inv_r[2][100], inv_r[2][101]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  28,   1,  34,  24,   4, inv_r[2][102], inv_r[2][103]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  60,   2,   1,   4,  56, inv_r[2][104], inv_r[2][105]);} // 6 times
    {psi.inv_tchannel_opt(no_ps,   6,  56,   1,   6,  24,  32, inv_r[2][106], inv_r[2][107]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  56,   1,   6,  40,  16, inv_r[2][108], inv_r[2][109]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  56,   1,   6,  48,   8, inv_r[2][110], inv_r[2][111]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  56,   2,   5,  48,   8, inv_r[2][112], inv_r[2][113]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  56,   2,   5,  24,  32, inv_r[2][114], inv_r[2][115]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  28,   2,  33,  24,   4, inv_r[2][116], inv_r[2][117]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  56,   2,   5,  40,  16, inv_r[2][118], inv_r[2][119]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  44,   2,  17,  40,   4, inv_r[2][120], inv_r[2][121]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  44,   2,  17,  12,  32, inv_r[2][122], inv_r[2][123]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  28,   2,  33,  12,  16, inv_r[2][124], inv_r[2][125]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  52,   2,   9,  48,   4, inv_r[2][126], inv_r[2][127]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  52,   2,   9,  20,  32, inv_r[2][128], inv_r[2][129]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  28,   2,  33,  20,   8, inv_r[2][130], inv_r[2][131]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  52,   2,   9,  36,  16, inv_r[2][132], inv_r[2][133]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  44,   2,  17,  36,   8, inv_r[2][134], inv_r[2][135]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   1,   2,  12,  48, inv_r[2][136], inv_r[2][137]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  60,   1,   2,  36,  24, inv_r[2][138], inv_r[2][139]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  60,   1,   2,  20,  40, inv_r[2][140], inv_r[2][141]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   0,  60,   1,   2,  48,  12, inv_r[2][142], inv_r[2][143]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  60,   1,   2,  40,  20, inv_r[2][144], inv_r[2][145]);} // 1 times
    {psi.inv_tchannel_opt(no_ps,   6,  60,   1,   2,  24,  36, inv_r[2][146], inv_r[2][147]);} // 1 times

    {psi.inv_decay(no_ps,  12,   4,   8, inv_r[3][0], inv_r[3][1]);} // 15 times
    {psi.inv_decay(no_ps,  36,   4,  32, inv_r[3][2], inv_r[3][3]);} // 15 times
    {psi.inv_decay(no_ps,  20,   4,  16, inv_r[3][4], inv_r[3][5]);} // 15 times
    {psi.inv_decay(no_ps,  48,  16,  32, inv_r[3][6], inv_r[3][7]);} // 15 times
    {psi.inv_decay(no_ps,  40,   8,  32, inv_r[3][8], inv_r[3][9]);} // 15 times
    {psi.inv_decay(no_ps,  24,   8,  16, inv_r[3][10], inv_r[3][11]);} // 15 times
    {psi.inv_decay(no_ps,  28,   4,  24, inv_r[3][12], inv_r[3][13]);} // 3 times
    {psi.inv_decay(no_ps,  44,   4,  40, inv_r[3][14], inv_r[3][15]);} // 3 times
    {psi.inv_decay(no_ps,  44,  32,  12, inv_r[3][16], inv_r[3][17]);} // 3 times
    {psi.inv_decay(no_ps,  52,   4,  48, inv_r[3][18], inv_r[3][19]);} // 3 times
    {psi.inv_decay(no_ps,  52,  32,  20, inv_r[3][20], inv_r[3][21]);} // 3 times
    {psi.inv_decay(no_ps,  28,  16,  12, inv_r[3][22], inv_r[3][23]);} // 3 times
    {psi.inv_decay(no_ps,  52,  16,  36, inv_r[3][24], inv_r[3][25]);} // 3 times
    {psi.inv_decay(no_ps,  28,   8,  20, inv_r[3][26], inv_r[3][27]);} // 3 times
    {psi.inv_decay(no_ps,  44,   8,  36, inv_r[3][28], inv_r[3][29]);} // 3 times
    {psi.inv_decay(no_ps,  56,   8,  48, inv_r[3][30], inv_r[3][31]);} // 3 times
    {psi.inv_decay(no_ps,  56,  32,  24, inv_r[3][32], inv_r[3][33]);} // 3 times
    {psi.inv_decay(no_ps,  56,  16,  40, inv_r[3][34], inv_r[3][35]);} // 3 times
    {psi.inv_decay(no_ps,  60,   4,  56, inv_r[3][36], inv_r[3][37]);} // 3 times
    {psi.inv_decay(no_ps,  60,  32,  28, inv_r[3][38], inv_r[3][39]);} // 3 times
    {psi.inv_decay(no_ps,  60,  16,  44, inv_r[3][40], inv_r[3][41]);} // 3 times
    {psi.inv_decay(no_ps,  60,   8,  52, inv_r[3][42], inv_r[3][43]);} // 3 times
    {psi.inv_decay(no_ps,  60,  12,  48, inv_r[3][44], inv_r[3][45]);} // 1 times
    {psi.inv_decay(no_ps,  60,  36,  24, inv_r[3][46], inv_r[3][47]);} // 1 times
    {psi.inv_decay(no_ps,  60,  20,  40, inv_r[3][48], inv_r[3][49]);} // 1 times

    if (psi_weight_IS == 2){
      int n_random = 3 * psi_n_particle - 4;;
      vector<vector<double> > inv_channel_r(105, vector<double> (8));
      vector<double> g_IS(105, 1.);
      for (int i_c = 0; i_c < 105; i_c++){
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
      if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = g_IS[  2] * psi_g_f[no_ps][  2] * psi_g_f[no_ps][  3] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  3] * psi_g_t[no_ps][  5];}
      if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = g_IS[  3] * psi_g_f[no_ps][  4] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_t[no_ps][  8];}
      if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = g_IS[  4] * psi_g_f[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_t[no_ps][  9];}
      if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = g_IS[  5] * psi_g_f[no_ps][  6] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 10] * psi_g_t[no_ps][ 11];}
      if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = g_IS[  6] * psi_g_f[no_ps][  6] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 10] * psi_g_t[no_ps][ 12];}
      if (psi_MC_alpha[zero +   7] != 0.){psi_MC_g_channel[zero +   7] = g_IS[  7] * psi_g_f[no_ps][  8] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 13] * psi_g_t[no_ps][ 14];}
      if (psi_MC_alpha[zero +   8] != 0.){psi_MC_g_channel[zero +   8] = g_IS[  8] * psi_g_f[no_ps][  4] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 16] * psi_g_t[no_ps][ 17];}
      if (psi_MC_alpha[zero +   9] != 0.){psi_MC_g_channel[zero +   9] = g_IS[  9] * psi_g_f[no_ps][ 11] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 18] * psi_g_t[no_ps][ 19];}
      if (psi_MC_alpha[zero +  10] != 0.){psi_MC_g_channel[zero +  10] = g_IS[ 10] * psi_g_f[no_ps][ 11] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 18] * psi_g_t[no_ps][ 20];}
      if (psi_MC_alpha[zero +  11] != 0.){psi_MC_g_channel[zero +  11] = g_IS[ 11] * psi_g_f[no_ps][  8] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 13] * psi_g_t[no_ps][ 21];}
      if (psi_MC_alpha[zero +  12] != 0.){psi_MC_g_channel[zero +  12] = g_IS[ 12] * psi_g_f[no_ps][  8] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 23] * psi_g_t[no_ps][ 24];}
      if (psi_MC_alpha[zero +  13] != 0.){psi_MC_g_channel[zero +  13] = g_IS[ 13] * psi_g_f[no_ps][ 11] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 25] * psi_g_t[no_ps][ 26];}
      if (psi_MC_alpha[zero +  14] != 0.){psi_MC_g_channel[zero +  14] = g_IS[ 14] * psi_g_f[no_ps][ 11] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 25] * psi_g_t[no_ps][ 27];}
      if (psi_MC_alpha[zero +  15] != 0.){psi_MC_g_channel[zero +  15] = g_IS[ 15] * psi_g_f[no_ps][  4] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_t[no_ps][ 28];}
      if (psi_MC_alpha[zero +  16] != 0.){psi_MC_g_channel[zero +  16] = g_IS[ 16] * psi_g_f[no_ps][  8] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 23] * psi_g_t[no_ps][ 29];}
      if (psi_MC_alpha[zero +  17] != 0.){psi_MC_g_channel[zero +  17] = g_IS[ 17] * psi_g_f[no_ps][  6] * psi_g_f[no_ps][ 15] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 30] * psi_g_t[no_ps][ 31];}
      if (psi_MC_alpha[zero +  18] != 0.){psi_MC_g_channel[zero +  18] = g_IS[ 18] * psi_g_f[no_ps][  6] * psi_g_f[no_ps][ 15] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 30] * psi_g_t[no_ps][ 32];}
      if (psi_MC_alpha[zero +  19] != 0.){psi_MC_g_channel[zero +  19] = g_IS[ 19] * psi_g_f[no_ps][  0] * psi_g_f[no_ps][ 16] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 33] * psi_g_t[no_ps][ 34];}
      if (psi_MC_alpha[zero +  20] != 0.){psi_MC_g_channel[zero +  20] = g_IS[ 20] * psi_g_f[no_ps][  4] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 16] * psi_g_t[no_ps][ 35];}
      if (psi_MC_alpha[zero +  21] != 0.){psi_MC_g_channel[zero +  21] = g_IS[ 21] * psi_g_f[no_ps][  2] * psi_g_f[no_ps][ 17] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 36] * psi_g_t[no_ps][ 37];}
      if (psi_MC_alpha[zero +  22] != 0.){psi_MC_g_channel[zero +  22] = g_IS[ 22] * psi_g_f[no_ps][  2] * psi_g_f[no_ps][ 17] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 36] * psi_g_t[no_ps][ 38];}
      if (psi_MC_alpha[zero +  23] != 0.){psi_MC_g_channel[zero +  23] = g_IS[ 23] * psi_g_f[no_ps][  0] * psi_g_f[no_ps][ 16] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 33] * psi_g_t[no_ps][ 39];}
      if (psi_MC_alpha[zero +  24] != 0.){psi_MC_g_channel[zero +  24] = g_IS[ 24] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][ 40] * psi_g_t[no_ps][ 41] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  25] != 0.){psi_MC_g_channel[zero +  25] = g_IS[ 25] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 42] * psi_g_t[no_ps][ 43] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  26] != 0.){psi_MC_g_channel[zero +  26] = g_IS[ 26] * psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 42] * psi_g_t[no_ps][ 44] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  27] != 0.){psi_MC_g_channel[zero +  27] = g_IS[ 27] * psi_g_p[no_ps][  2] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][ 40] * psi_g_t[no_ps][ 45] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  28] != 0.){psi_MC_g_channel[zero +  28] = g_IS[ 28] * psi_g_p[no_ps][  2] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][ 46] * psi_g_t[no_ps][ 47] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  29] != 0.){psi_MC_g_channel[zero +  29] = g_IS[ 29] * psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][ 46] * psi_g_t[no_ps][ 48] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  30] != 0.){psi_MC_g_channel[zero +  30] = g_IS[ 30] * psi_g_p[no_ps][  3] * psi_g_f[no_ps][ 15] * psi_g_t[no_ps][ 46] * psi_g_t[no_ps][ 49] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  31] != 0.){psi_MC_g_channel[zero +  31] = g_IS[ 31] * psi_g_p[no_ps][  4] * psi_g_f[no_ps][ 17] * psi_g_t[no_ps][ 42] * psi_g_t[no_ps][ 50] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  32] != 0.){psi_MC_g_channel[zero +  32] = g_IS[ 32] * psi_g_p[no_ps][  5] * psi_g_f[no_ps][ 16] * psi_g_t[no_ps][ 40] * psi_g_t[no_ps][ 51] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  33] != 0.){psi_MC_g_channel[zero +  33] = g_IS[ 33] * psi_g_p[no_ps][  5] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][ 52] * psi_g_t[no_ps][ 53] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  34] != 0.){psi_MC_g_channel[zero +  34] = g_IS[ 34] * psi_g_p[no_ps][  4] * psi_g_f[no_ps][  3] * psi_g_t[no_ps][ 52] * psi_g_t[no_ps][ 54] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  35] != 0.){psi_MC_g_channel[zero +  35] = g_IS[ 35] * psi_g_p[no_ps][  3] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][ 52] * psi_g_t[no_ps][ 55] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  36] != 0.){psi_MC_g_channel[zero +  36] = g_IS[ 36] * psi_g_p[no_ps][  5] * psi_g_p[no_ps][  6] * psi_g_t[no_ps][ 40] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  37] != 0.){psi_MC_g_channel[zero +  37] = g_IS[ 37] * psi_g_p[no_ps][  4] * psi_g_p[no_ps][  7] * psi_g_t[no_ps][ 42] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  38] != 0.){psi_MC_g_channel[zero +  38] = g_IS[ 38] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  8] * psi_g_t[no_ps][ 42] * psi_g_d[no_ps][  8] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  39] != 0.){psi_MC_g_channel[zero +  39] = g_IS[ 39] * psi_g_p[no_ps][  3] * psi_g_p[no_ps][  9] * psi_g_t[no_ps][ 46] * psi_g_d[no_ps][  9] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  40] != 0.){psi_MC_g_channel[zero +  40] = g_IS[ 40] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 10] * psi_g_t[no_ps][ 46] * psi_g_d[no_ps][ 10] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  41] != 0.){psi_MC_g_channel[zero +  41] = g_IS[ 41] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][ 11] * psi_g_t[no_ps][ 40] * psi_g_d[no_ps][ 11] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  42] != 0.){psi_MC_g_channel[zero +  42] = g_IS[ 42] * psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 12] * psi_g_t[no_ps][ 46] * psi_g_d[no_ps][ 12] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  43] != 0.){psi_MC_g_channel[zero +  43] = g_IS[ 43] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 13] * psi_g_t[no_ps][ 40] * psi_g_d[no_ps][ 13] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  44] != 0.){psi_MC_g_channel[zero +  44] = g_IS[ 44] * psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 14] * psi_g_t[no_ps][ 42] * psi_g_d[no_ps][ 14] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  45] != 0.){psi_MC_g_channel[zero +  45] = g_IS[ 45] * psi_g_p[no_ps][  3] * psi_g_p[no_ps][ 15] * psi_g_t[no_ps][ 52] * psi_g_d[no_ps][ 15] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  46] != 0.){psi_MC_g_channel[zero +  46] = g_IS[ 46] * psi_g_p[no_ps][  5] * psi_g_p[no_ps][ 16] * psi_g_t[no_ps][ 52] * psi_g_d[no_ps][ 16] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  47] != 0.){psi_MC_g_channel[zero +  47] = g_IS[ 47] * psi_g_p[no_ps][  4] * psi_g_p[no_ps][ 17] * psi_g_t[no_ps][ 52] * psi_g_d[no_ps][ 17] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  48] != 0.){psi_MC_g_channel[zero +  48] = g_IS[ 48] * psi_g_p[no_ps][  3] * psi_g_p[no_ps][ 15] * psi_g_d[no_ps][ 18] * psi_g_d[no_ps][ 15] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  49] != 0.){psi_MC_g_channel[zero +  49] = g_IS[ 49] * psi_g_p[no_ps][  5] * psi_g_p[no_ps][ 16] * psi_g_d[no_ps][ 18] * psi_g_d[no_ps][ 16] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  50] != 0.){psi_MC_g_channel[zero +  50] = g_IS[ 50] * psi_g_p[no_ps][  5] * psi_g_p[no_ps][  6] * psi_g_d[no_ps][ 19] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  51] != 0.){psi_MC_g_channel[zero +  51] = g_IS[ 51] * psi_g_p[no_ps][  4] * psi_g_p[no_ps][ 17] * psi_g_d[no_ps][ 18] * psi_g_d[no_ps][ 17] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  52] != 0.){psi_MC_g_channel[zero +  52] = g_IS[ 52] * psi_g_p[no_ps][  4] * psi_g_p[no_ps][  7] * psi_g_d[no_ps][ 20] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  53] != 0.){psi_MC_g_channel[zero +  53] = g_IS[ 53] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  8] * psi_g_d[no_ps][ 20] * psi_g_d[no_ps][  8] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  54] != 0.){psi_MC_g_channel[zero +  54] = g_IS[ 54] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][ 11] * psi_g_d[no_ps][ 19] * psi_g_d[no_ps][ 11] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  55] != 0.){psi_MC_g_channel[zero +  55] = g_IS[ 55] * psi_g_p[no_ps][  3] * psi_g_p[no_ps][  9] * psi_g_d[no_ps][ 21] * psi_g_d[no_ps][  9] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  56] != 0.){psi_MC_g_channel[zero +  56] = g_IS[ 56] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 10] * psi_g_d[no_ps][ 21] * psi_g_d[no_ps][ 10] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  57] != 0.){psi_MC_g_channel[zero +  57] = g_IS[ 57] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 13] * psi_g_d[no_ps][ 19] * psi_g_d[no_ps][ 13] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  58] != 0.){psi_MC_g_channel[zero +  58] = g_IS[ 58] * psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 12] * psi_g_d[no_ps][ 21] * psi_g_d[no_ps][ 12] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  59] != 0.){psi_MC_g_channel[zero +  59] = g_IS[ 59] * psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 14] * psi_g_d[no_ps][ 20] * psi_g_d[no_ps][ 14] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  60] != 0.){psi_MC_g_channel[zero +  60] = g_IS[ 60] * psi_g_p[no_ps][  3] * psi_g_p[no_ps][ 18] * psi_g_d[no_ps][ 22] * psi_g_d[no_ps][  0] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  61] != 0.){psi_MC_g_channel[zero +  61] = g_IS[ 61] * psi_g_p[no_ps][  5] * psi_g_p[no_ps][ 19] * psi_g_d[no_ps][ 23] * psi_g_d[no_ps][  1] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  62] != 0.){psi_MC_g_channel[zero +  62] = g_IS[ 62] * psi_g_p[no_ps][  4] * psi_g_p[no_ps][ 20] * psi_g_d[no_ps][ 24] * psi_g_d[no_ps][  2] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  63] != 0.){psi_MC_g_channel[zero +  63] = g_IS[ 63] * psi_g_p[no_ps][  5] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  64] != 0.){psi_MC_g_channel[zero +  64] = g_IS[ 64] * psi_g_p[no_ps][  4] * psi_g_f[no_ps][  3] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  3] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  65] != 0.){psi_MC_g_channel[zero +  65] = g_IS[ 65] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  66] != 0.){psi_MC_g_channel[zero +  66] = g_IS[ 66] * psi_g_p[no_ps][  3] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 10] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  67] != 0.){psi_MC_g_channel[zero +  67] = g_IS[ 67] * psi_g_p[no_ps][  2] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 13] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  68] != 0.){psi_MC_g_channel[zero +  68] = g_IS[ 68] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 16] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  69] != 0.){psi_MC_g_channel[zero +  69] = g_IS[ 69] * psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 18] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  70] != 0.){psi_MC_g_channel[zero +  70] = g_IS[ 70] * psi_g_p[no_ps][  2] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 23] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  71] != 0.){psi_MC_g_channel[zero +  71] = g_IS[ 71] * psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 25] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  72] != 0.){psi_MC_g_channel[zero +  72] = g_IS[ 72] * psi_g_p[no_ps][  3] * psi_g_f[no_ps][ 15] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 30] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  73] != 0.){psi_MC_g_channel[zero +  73] = g_IS[ 73] * psi_g_p[no_ps][  5] * psi_g_f[no_ps][ 16] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 33] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  74] != 0.){psi_MC_g_channel[zero +  74] = g_IS[ 74] * psi_g_p[no_ps][  4] * psi_g_f[no_ps][ 17] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 36] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  75] != 0.){psi_MC_g_channel[zero +  75] = g_IS[ 75] * psi_g_p[no_ps][  3] * psi_g_p[no_ps][ 15] * psi_g_t[no_ps][  0] * psi_g_d[no_ps][ 15] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  76] != 0.){psi_MC_g_channel[zero +  76] = g_IS[ 76] * psi_g_p[no_ps][  5] * psi_g_p[no_ps][ 16] * psi_g_t[no_ps][  0] * psi_g_d[no_ps][ 16] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  77] != 0.){psi_MC_g_channel[zero +  77] = g_IS[ 77] * psi_g_p[no_ps][  5] * psi_g_p[no_ps][  6] * psi_g_t[no_ps][  6] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  78] != 0.){psi_MC_g_channel[zero +  78] = g_IS[ 78] * psi_g_p[no_ps][  4] * psi_g_p[no_ps][ 17] * psi_g_t[no_ps][  0] * psi_g_d[no_ps][ 17] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  79] != 0.){psi_MC_g_channel[zero +  79] = g_IS[ 79] * psi_g_p[no_ps][  4] * psi_g_p[no_ps][  7] * psi_g_t[no_ps][ 15] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  80] != 0.){psi_MC_g_channel[zero +  80] = g_IS[ 80] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][  8] * psi_g_t[no_ps][ 15] * psi_g_d[no_ps][  8] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  81] != 0.){psi_MC_g_channel[zero +  81] = g_IS[ 81] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][ 11] * psi_g_t[no_ps][  6] * psi_g_d[no_ps][ 11] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  82] != 0.){psi_MC_g_channel[zero +  82] = g_IS[ 82] * psi_g_p[no_ps][  3] * psi_g_p[no_ps][  9] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][  9] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  83] != 0.){psi_MC_g_channel[zero +  83] = g_IS[ 83] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 10] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][ 10] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  84] != 0.){psi_MC_g_channel[zero +  84] = g_IS[ 84] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 13] * psi_g_t[no_ps][  6] * psi_g_d[no_ps][ 13] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  85] != 0.){psi_MC_g_channel[zero +  85] = g_IS[ 85] * psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 12] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][ 12] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  86] != 0.){psi_MC_g_channel[zero +  86] = g_IS[ 86] * psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 14] * psi_g_t[no_ps][ 15] * psi_g_d[no_ps][ 14] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  87] != 0.){psi_MC_g_channel[zero +  87] = g_IS[ 87] * psi_g_p[no_ps][  3] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 56] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  88] != 0.){psi_MC_g_channel[zero +  88] = g_IS[ 88] * psi_g_p[no_ps][  5] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 57] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  89] != 0.){psi_MC_g_channel[zero +  89] = g_IS[ 89] * psi_g_p[no_ps][  5] * psi_g_f[no_ps][ 16] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 58] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero +  90] != 0.){psi_MC_g_channel[zero +  90] = g_IS[ 90] * psi_g_p[no_ps][  4] * psi_g_f[no_ps][  3] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 59] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  91] != 0.){psi_MC_g_channel[zero +  91] = g_IS[ 91] * psi_g_p[no_ps][  4] * psi_g_f[no_ps][ 17] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 60] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero +  92] != 0.){psi_MC_g_channel[zero +  92] = g_IS[ 92] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 61] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  93] != 0.){psi_MC_g_channel[zero +  93] = g_IS[ 93] * psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 62] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero +  94] != 0.){psi_MC_g_channel[zero +  94] = g_IS[ 94] * psi_g_p[no_ps][  3] * psi_g_f[no_ps][ 15] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 63] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero +  95] != 0.){psi_MC_g_channel[zero +  95] = g_IS[ 95] * psi_g_p[no_ps][  2] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 64] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  96] != 0.){psi_MC_g_channel[zero +  96] = g_IS[ 96] * psi_g_p[no_ps][  2] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 65] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero +  97] != 0.){psi_MC_g_channel[zero +  97] = g_IS[ 97] * psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 66] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  98] != 0.){psi_MC_g_channel[zero +  98] = g_IS[ 98] * psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 67] * psi_g_d[no_ps][  1];}
      if (psi_MC_alpha[zero +  99] != 0.){psi_MC_g_channel[zero +  99] = g_IS[ 99] * psi_g_p[no_ps][  3] * psi_g_p[no_ps][ 18] * psi_g_t[no_ps][ 68] * psi_g_d[no_ps][  0] * psi_g_d[no_ps][  3];}
      if (psi_MC_alpha[zero + 100] != 0.){psi_MC_g_channel[zero + 100] = g_IS[100] * psi_g_p[no_ps][  5] * psi_g_p[no_ps][ 19] * psi_g_t[no_ps][ 69] * psi_g_d[no_ps][  1] * psi_g_d[no_ps][  5];}
      if (psi_MC_alpha[zero + 101] != 0.){psi_MC_g_channel[zero + 101] = g_IS[101] * psi_g_p[no_ps][  4] * psi_g_p[no_ps][ 20] * psi_g_t[no_ps][ 70] * psi_g_d[no_ps][  2] * psi_g_d[no_ps][  4];}
      if (psi_MC_alpha[zero + 102] != 0.){psi_MC_g_channel[zero + 102] = g_IS[102] * psi_g_p[no_ps][  0] * psi_g_p[no_ps][ 21] * psi_g_t[no_ps][ 71] * psi_g_d[no_ps][  3] * psi_g_d[no_ps][  0];}
      if (psi_MC_alpha[zero + 103] != 0.){psi_MC_g_channel[zero + 103] = g_IS[103] * psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 22] * psi_g_t[no_ps][ 72] * psi_g_d[no_ps][  4] * psi_g_d[no_ps][  2];}
      if (psi_MC_alpha[zero + 104] != 0.){psi_MC_g_channel[zero + 104] = g_IS[104] * psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 23] * psi_g_t[no_ps][ 73] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  1];}

    }
    else if (psi_weight_IS == 4){
      vector<vector<double> > g_IS(4);
      g_IS[0].resize(24, 1.);
      g_IS[1].resize(18, 1.);
      g_IS[2].resize(74, 1.);
      g_IS[3].resize(25, 1.);

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
      if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = psi_g_f[no_ps][  2] * g_IS[1][  2] * psi_g_f[no_ps][  3] * g_IS[1][  3] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  3] * g_IS[2][  3] * psi_g_t[no_ps][  5] * g_IS[2][  5];}
      if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = psi_g_f[no_ps][  4] * g_IS[1][  4] * psi_g_f[no_ps][  5] * g_IS[1][  5] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][  7] * g_IS[2][  7] * psi_g_t[no_ps][  8] * g_IS[2][  8];}
      if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = psi_g_f[no_ps][  0] * g_IS[1][  0] * psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  1] * g_IS[2][  1] * psi_g_t[no_ps][  9] * g_IS[2][  9];}
      if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = psi_g_f[no_ps][  6] * g_IS[1][  6] * psi_g_f[no_ps][  7] * g_IS[1][  7] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][ 10] * g_IS[2][ 10] * psi_g_t[no_ps][ 11] * g_IS[2][ 11];}
      if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = psi_g_f[no_ps][  6] * g_IS[1][  6] * psi_g_f[no_ps][  7] * g_IS[1][  7] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][ 10] * g_IS[2][ 10] * psi_g_t[no_ps][ 12] * g_IS[2][ 12];}
      if (psi_MC_alpha[zero +   7] != 0.){psi_MC_g_channel[zero +   7] = psi_g_f[no_ps][  8] * g_IS[1][  8] * psi_g_f[no_ps][  9] * g_IS[1][  9] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][ 13] * g_IS[2][ 13] * psi_g_t[no_ps][ 14] * g_IS[2][ 14];}
      if (psi_MC_alpha[zero +   8] != 0.){psi_MC_g_channel[zero +   8] = psi_g_f[no_ps][  4] * g_IS[1][  4] * psi_g_f[no_ps][ 10] * g_IS[1][ 10] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 16] * g_IS[2][ 16] * psi_g_t[no_ps][ 17] * g_IS[2][ 17];}
      if (psi_MC_alpha[zero +   9] != 0.){psi_MC_g_channel[zero +   9] = psi_g_f[no_ps][ 11] * g_IS[1][ 11] * psi_g_f[no_ps][ 12] * g_IS[1][ 12] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 18] * g_IS[2][ 18] * psi_g_t[no_ps][ 19] * g_IS[2][ 19];}
      if (psi_MC_alpha[zero +  10] != 0.){psi_MC_g_channel[zero +  10] = psi_g_f[no_ps][ 11] * g_IS[1][ 11] * psi_g_f[no_ps][ 12] * g_IS[1][ 12] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 18] * g_IS[2][ 18] * psi_g_t[no_ps][ 20] * g_IS[2][ 20];}
      if (psi_MC_alpha[zero +  11] != 0.){psi_MC_g_channel[zero +  11] = psi_g_f[no_ps][  8] * g_IS[1][  8] * psi_g_f[no_ps][  9] * g_IS[1][  9] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][ 13] * g_IS[2][ 13] * psi_g_t[no_ps][ 21] * g_IS[2][ 21];}
      if (psi_MC_alpha[zero +  12] != 0.){psi_MC_g_channel[zero +  12] = psi_g_f[no_ps][  8] * g_IS[1][  8] * psi_g_f[no_ps][ 13] * g_IS[1][ 13] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 23] * g_IS[2][ 23] * psi_g_t[no_ps][ 24] * g_IS[2][ 24];}
      if (psi_MC_alpha[zero +  13] != 0.){psi_MC_g_channel[zero +  13] = psi_g_f[no_ps][ 11] * g_IS[1][ 11] * psi_g_f[no_ps][ 14] * g_IS[1][ 14] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 25] * g_IS[2][ 25] * psi_g_t[no_ps][ 26] * g_IS[2][ 26];}
      if (psi_MC_alpha[zero +  14] != 0.){psi_MC_g_channel[zero +  14] = psi_g_f[no_ps][ 11] * g_IS[1][ 11] * psi_g_f[no_ps][ 14] * g_IS[1][ 14] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 25] * g_IS[2][ 25] * psi_g_t[no_ps][ 27] * g_IS[2][ 27];}
      if (psi_MC_alpha[zero +  15] != 0.){psi_MC_g_channel[zero +  15] = psi_g_f[no_ps][  4] * g_IS[1][  4] * psi_g_f[no_ps][  5] * g_IS[1][  5] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][  7] * g_IS[2][  7] * psi_g_t[no_ps][ 28] * g_IS[2][ 28];}
      if (psi_MC_alpha[zero +  16] != 0.){psi_MC_g_channel[zero +  16] = psi_g_f[no_ps][  8] * g_IS[1][  8] * psi_g_f[no_ps][ 13] * g_IS[1][ 13] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 23] * g_IS[2][ 23] * psi_g_t[no_ps][ 29] * g_IS[2][ 29];}
      if (psi_MC_alpha[zero +  17] != 0.){psi_MC_g_channel[zero +  17] = psi_g_f[no_ps][  6] * g_IS[1][  6] * psi_g_f[no_ps][ 15] * g_IS[1][ 15] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 30] * g_IS[2][ 30] * psi_g_t[no_ps][ 31] * g_IS[2][ 31];}
      if (psi_MC_alpha[zero +  18] != 0.){psi_MC_g_channel[zero +  18] = psi_g_f[no_ps][  6] * g_IS[1][  6] * psi_g_f[no_ps][ 15] * g_IS[1][ 15] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 30] * g_IS[2][ 30] * psi_g_t[no_ps][ 32] * g_IS[2][ 32];}
      if (psi_MC_alpha[zero +  19] != 0.){psi_MC_g_channel[zero +  19] = psi_g_f[no_ps][  0] * g_IS[1][  0] * psi_g_f[no_ps][ 16] * g_IS[1][ 16] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][ 33] * g_IS[2][ 33] * psi_g_t[no_ps][ 34] * g_IS[2][ 34];}
      if (psi_MC_alpha[zero +  20] != 0.){psi_MC_g_channel[zero +  20] = psi_g_f[no_ps][  4] * g_IS[1][  4] * psi_g_f[no_ps][ 10] * g_IS[1][ 10] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 16] * g_IS[2][ 16] * psi_g_t[no_ps][ 35] * g_IS[2][ 35];}
      if (psi_MC_alpha[zero +  21] != 0.){psi_MC_g_channel[zero +  21] = psi_g_f[no_ps][  2] * g_IS[1][  2] * psi_g_f[no_ps][ 17] * g_IS[1][ 17] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 36] * g_IS[2][ 36] * psi_g_t[no_ps][ 37] * g_IS[2][ 37];}
      if (psi_MC_alpha[zero +  22] != 0.){psi_MC_g_channel[zero +  22] = psi_g_f[no_ps][  2] * g_IS[1][  2] * psi_g_f[no_ps][ 17] * g_IS[1][ 17] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 36] * g_IS[2][ 36] * psi_g_t[no_ps][ 38] * g_IS[2][ 38];}
      if (psi_MC_alpha[zero +  23] != 0.){psi_MC_g_channel[zero +  23] = psi_g_f[no_ps][  0] * g_IS[1][  0] * psi_g_f[no_ps][ 16] * g_IS[1][ 16] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][ 33] * g_IS[2][ 33] * psi_g_t[no_ps][ 39] * g_IS[2][ 39];}
      if (psi_MC_alpha[zero +  24] != 0.){psi_MC_g_channel[zero +  24] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][  5] * g_IS[1][  5] * psi_g_t[no_ps][ 40] * g_IS[2][ 40] * psi_g_t[no_ps][ 41] * g_IS[2][ 41] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  25] != 0.){psi_MC_g_channel[zero +  25] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][ 10] * g_IS[1][ 10] * psi_g_t[no_ps][ 42] * g_IS[2][ 42] * psi_g_t[no_ps][ 43] * g_IS[2][ 43] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  26] != 0.){psi_MC_g_channel[zero +  26] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_f[no_ps][ 12] * g_IS[1][ 12] * psi_g_t[no_ps][ 42] * g_IS[2][ 42] * psi_g_t[no_ps][ 44] * g_IS[2][ 44] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  27] != 0.){psi_MC_g_channel[zero +  27] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_f[no_ps][  9] * g_IS[1][  9] * psi_g_t[no_ps][ 40] * g_IS[2][ 40] * psi_g_t[no_ps][ 45] * g_IS[2][ 45] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  28] != 0.){psi_MC_g_channel[zero +  28] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_f[no_ps][ 13] * g_IS[1][ 13] * psi_g_t[no_ps][ 46] * g_IS[2][ 46] * psi_g_t[no_ps][ 47] * g_IS[2][ 47] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  29] != 0.){psi_MC_g_channel[zero +  29] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_f[no_ps][ 14] * g_IS[1][ 14] * psi_g_t[no_ps][ 46] * g_IS[2][ 46] * psi_g_t[no_ps][ 48] * g_IS[2][ 48] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  30] != 0.){psi_MC_g_channel[zero +  30] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_f[no_ps][ 15] * g_IS[1][ 15] * psi_g_t[no_ps][ 46] * g_IS[2][ 46] * psi_g_t[no_ps][ 49] * g_IS[2][ 49] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  31] != 0.){psi_MC_g_channel[zero +  31] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_f[no_ps][ 17] * g_IS[1][ 17] * psi_g_t[no_ps][ 42] * g_IS[2][ 42] * psi_g_t[no_ps][ 50] * g_IS[2][ 50] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  32] != 0.){psi_MC_g_channel[zero +  32] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_f[no_ps][ 16] * g_IS[1][ 16] * psi_g_t[no_ps][ 40] * g_IS[2][ 40] * psi_g_t[no_ps][ 51] * g_IS[2][ 51] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  33] != 0.){psi_MC_g_channel[zero +  33] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][ 52] * g_IS[2][ 52] * psi_g_t[no_ps][ 53] * g_IS[2][ 53] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  34] != 0.){psi_MC_g_channel[zero +  34] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_f[no_ps][  3] * g_IS[1][  3] * psi_g_t[no_ps][ 52] * g_IS[2][ 52] * psi_g_t[no_ps][ 54] * g_IS[2][ 54] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  35] != 0.){psi_MC_g_channel[zero +  35] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_f[no_ps][  7] * g_IS[1][  7] * psi_g_t[no_ps][ 52] * g_IS[2][ 52] * psi_g_t[no_ps][ 55] * g_IS[2][ 55] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  36] != 0.){psi_MC_g_channel[zero +  36] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_p[no_ps][  6] * g_IS[0][  6] * psi_g_t[no_ps][ 40] * g_IS[2][ 40] * psi_g_d[no_ps][  6] * g_IS[3][  6] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  37] != 0.){psi_MC_g_channel[zero +  37] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_p[no_ps][  7] * g_IS[0][  7] * psi_g_t[no_ps][ 42] * g_IS[2][ 42] * psi_g_d[no_ps][  7] * g_IS[3][  7] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  38] != 0.){psi_MC_g_channel[zero +  38] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  8] * g_IS[0][  8] * psi_g_t[no_ps][ 42] * g_IS[2][ 42] * psi_g_d[no_ps][  8] * g_IS[3][  8] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  39] != 0.){psi_MC_g_channel[zero +  39] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_p[no_ps][  9] * g_IS[0][  9] * psi_g_t[no_ps][ 46] * g_IS[2][ 46] * psi_g_d[no_ps][  9] * g_IS[3][  9] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  40] != 0.){psi_MC_g_channel[zero +  40] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][ 10] * g_IS[0][ 10] * psi_g_t[no_ps][ 46] * g_IS[2][ 46] * psi_g_d[no_ps][ 10] * g_IS[3][ 10] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  41] != 0.){psi_MC_g_channel[zero +  41] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][ 11] * g_IS[0][ 11] * psi_g_t[no_ps][ 40] * g_IS[2][ 40] * psi_g_d[no_ps][ 11] * g_IS[3][ 11] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  42] != 0.){psi_MC_g_channel[zero +  42] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_p[no_ps][ 12] * g_IS[0][ 12] * psi_g_t[no_ps][ 46] * g_IS[2][ 46] * psi_g_d[no_ps][ 12] * g_IS[3][ 12] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  43] != 0.){psi_MC_g_channel[zero +  43] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][ 13] * g_IS[0][ 13] * psi_g_t[no_ps][ 40] * g_IS[2][ 40] * psi_g_d[no_ps][ 13] * g_IS[3][ 13] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  44] != 0.){psi_MC_g_channel[zero +  44] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_p[no_ps][ 14] * g_IS[0][ 14] * psi_g_t[no_ps][ 42] * g_IS[2][ 42] * psi_g_d[no_ps][ 14] * g_IS[3][ 14] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  45] != 0.){psi_MC_g_channel[zero +  45] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_p[no_ps][ 15] * g_IS[0][ 15] * psi_g_t[no_ps][ 52] * g_IS[2][ 52] * psi_g_d[no_ps][ 15] * g_IS[3][ 15] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  46] != 0.){psi_MC_g_channel[zero +  46] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_p[no_ps][ 16] * g_IS[0][ 16] * psi_g_t[no_ps][ 52] * g_IS[2][ 52] * psi_g_d[no_ps][ 16] * g_IS[3][ 16] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  47] != 0.){psi_MC_g_channel[zero +  47] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_p[no_ps][ 17] * g_IS[0][ 17] * psi_g_t[no_ps][ 52] * g_IS[2][ 52] * psi_g_d[no_ps][ 17] * g_IS[3][ 17] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  48] != 0.){psi_MC_g_channel[zero +  48] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_p[no_ps][ 15] * g_IS[0][ 15] * psi_g_d[no_ps][ 18] * g_IS[3][ 18] * psi_g_d[no_ps][ 15] * g_IS[3][ 15] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  49] != 0.){psi_MC_g_channel[zero +  49] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_p[no_ps][ 16] * g_IS[0][ 16] * psi_g_d[no_ps][ 18] * g_IS[3][ 18] * psi_g_d[no_ps][ 16] * g_IS[3][ 16] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  50] != 0.){psi_MC_g_channel[zero +  50] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_p[no_ps][  6] * g_IS[0][  6] * psi_g_d[no_ps][ 19] * g_IS[3][ 19] * psi_g_d[no_ps][  6] * g_IS[3][  6] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  51] != 0.){psi_MC_g_channel[zero +  51] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_p[no_ps][ 17] * g_IS[0][ 17] * psi_g_d[no_ps][ 18] * g_IS[3][ 18] * psi_g_d[no_ps][ 17] * g_IS[3][ 17] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  52] != 0.){psi_MC_g_channel[zero +  52] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_p[no_ps][  7] * g_IS[0][  7] * psi_g_d[no_ps][ 20] * g_IS[3][ 20] * psi_g_d[no_ps][  7] * g_IS[3][  7] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  53] != 0.){psi_MC_g_channel[zero +  53] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  8] * g_IS[0][  8] * psi_g_d[no_ps][ 20] * g_IS[3][ 20] * psi_g_d[no_ps][  8] * g_IS[3][  8] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  54] != 0.){psi_MC_g_channel[zero +  54] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][ 11] * g_IS[0][ 11] * psi_g_d[no_ps][ 19] * g_IS[3][ 19] * psi_g_d[no_ps][ 11] * g_IS[3][ 11] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  55] != 0.){psi_MC_g_channel[zero +  55] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_p[no_ps][  9] * g_IS[0][  9] * psi_g_d[no_ps][ 21] * g_IS[3][ 21] * psi_g_d[no_ps][  9] * g_IS[3][  9] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  56] != 0.){psi_MC_g_channel[zero +  56] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][ 10] * g_IS[0][ 10] * psi_g_d[no_ps][ 21] * g_IS[3][ 21] * psi_g_d[no_ps][ 10] * g_IS[3][ 10] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  57] != 0.){psi_MC_g_channel[zero +  57] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][ 13] * g_IS[0][ 13] * psi_g_d[no_ps][ 19] * g_IS[3][ 19] * psi_g_d[no_ps][ 13] * g_IS[3][ 13] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  58] != 0.){psi_MC_g_channel[zero +  58] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_p[no_ps][ 12] * g_IS[0][ 12] * psi_g_d[no_ps][ 21] * g_IS[3][ 21] * psi_g_d[no_ps][ 12] * g_IS[3][ 12] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  59] != 0.){psi_MC_g_channel[zero +  59] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_p[no_ps][ 14] * g_IS[0][ 14] * psi_g_d[no_ps][ 20] * g_IS[3][ 20] * psi_g_d[no_ps][ 14] * g_IS[3][ 14] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  60] != 0.){psi_MC_g_channel[zero +  60] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_p[no_ps][ 18] * g_IS[0][ 18] * psi_g_d[no_ps][ 22] * g_IS[3][ 22] * psi_g_d[no_ps][  0] * g_IS[3][  0] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  61] != 0.){psi_MC_g_channel[zero +  61] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_p[no_ps][ 19] * g_IS[0][ 19] * psi_g_d[no_ps][ 23] * g_IS[3][ 23] * psi_g_d[no_ps][  1] * g_IS[3][  1] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  62] != 0.){psi_MC_g_channel[zero +  62] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_p[no_ps][ 20] * g_IS[0][ 20] * psi_g_d[no_ps][ 24] * g_IS[3][ 24] * psi_g_d[no_ps][  2] * g_IS[3][  2] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  63] != 0.){psi_MC_g_channel[zero +  63] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  1] * g_IS[2][  1] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  64] != 0.){psi_MC_g_channel[zero +  64] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_f[no_ps][  3] * g_IS[1][  3] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][  3] * g_IS[2][  3] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  65] != 0.){psi_MC_g_channel[zero +  65] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][  5] * g_IS[1][  5] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][  7] * g_IS[2][  7] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  66] != 0.){psi_MC_g_channel[zero +  66] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_f[no_ps][  7] * g_IS[1][  7] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][ 10] * g_IS[2][ 10] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  67] != 0.){psi_MC_g_channel[zero +  67] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_f[no_ps][  9] * g_IS[1][  9] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][ 13] * g_IS[2][ 13] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  68] != 0.){psi_MC_g_channel[zero +  68] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][ 10] * g_IS[1][ 10] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 16] * g_IS[2][ 16] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  69] != 0.){psi_MC_g_channel[zero +  69] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_f[no_ps][ 12] * g_IS[1][ 12] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 18] * g_IS[2][ 18] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  70] != 0.){psi_MC_g_channel[zero +  70] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_f[no_ps][ 13] * g_IS[1][ 13] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 23] * g_IS[2][ 23] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  71] != 0.){psi_MC_g_channel[zero +  71] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_f[no_ps][ 14] * g_IS[1][ 14] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 25] * g_IS[2][ 25] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  72] != 0.){psi_MC_g_channel[zero +  72] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_f[no_ps][ 15] * g_IS[1][ 15] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 30] * g_IS[2][ 30] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  73] != 0.){psi_MC_g_channel[zero +  73] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_f[no_ps][ 16] * g_IS[1][ 16] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][ 33] * g_IS[2][ 33] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  74] != 0.){psi_MC_g_channel[zero +  74] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_f[no_ps][ 17] * g_IS[1][ 17] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 36] * g_IS[2][ 36] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  75] != 0.){psi_MC_g_channel[zero +  75] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_p[no_ps][ 15] * g_IS[0][ 15] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_d[no_ps][ 15] * g_IS[3][ 15] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  76] != 0.){psi_MC_g_channel[zero +  76] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_p[no_ps][ 16] * g_IS[0][ 16] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_d[no_ps][ 16] * g_IS[3][ 16] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  77] != 0.){psi_MC_g_channel[zero +  77] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_p[no_ps][  6] * g_IS[0][  6] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_d[no_ps][  6] * g_IS[3][  6] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  78] != 0.){psi_MC_g_channel[zero +  78] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_p[no_ps][ 17] * g_IS[0][ 17] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_d[no_ps][ 17] * g_IS[3][ 17] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  79] != 0.){psi_MC_g_channel[zero +  79] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_p[no_ps][  7] * g_IS[0][  7] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_d[no_ps][  7] * g_IS[3][  7] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  80] != 0.){psi_MC_g_channel[zero +  80] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][  8] * g_IS[0][  8] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_d[no_ps][  8] * g_IS[3][  8] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  81] != 0.){psi_MC_g_channel[zero +  81] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][ 11] * g_IS[0][ 11] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_d[no_ps][ 11] * g_IS[3][ 11] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  82] != 0.){psi_MC_g_channel[zero +  82] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_p[no_ps][  9] * g_IS[0][  9] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_d[no_ps][  9] * g_IS[3][  9] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  83] != 0.){psi_MC_g_channel[zero +  83] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][ 10] * g_IS[0][ 10] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_d[no_ps][ 10] * g_IS[3][ 10] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  84] != 0.){psi_MC_g_channel[zero +  84] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][ 13] * g_IS[0][ 13] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_d[no_ps][ 13] * g_IS[3][ 13] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  85] != 0.){psi_MC_g_channel[zero +  85] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_p[no_ps][ 12] * g_IS[0][ 12] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_d[no_ps][ 12] * g_IS[3][ 12] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  86] != 0.){psi_MC_g_channel[zero +  86] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_p[no_ps][ 14] * g_IS[0][ 14] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_d[no_ps][ 14] * g_IS[3][ 14] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  87] != 0.){psi_MC_g_channel[zero +  87] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_f[no_ps][  7] * g_IS[1][  7] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][ 56] * g_IS[2][ 56] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  88] != 0.){psi_MC_g_channel[zero +  88] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_f[no_ps][  1] * g_IS[1][  1] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][ 57] * g_IS[2][ 57] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  89] != 0.){psi_MC_g_channel[zero +  89] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_f[no_ps][ 16] * g_IS[1][ 16] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][ 58] * g_IS[2][ 58] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero +  90] != 0.){psi_MC_g_channel[zero +  90] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_f[no_ps][  3] * g_IS[1][  3] * psi_g_t[no_ps][  0] * g_IS[2][  0] * psi_g_t[no_ps][ 59] * g_IS[2][ 59] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  91] != 0.){psi_MC_g_channel[zero +  91] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_f[no_ps][ 17] * g_IS[1][ 17] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 60] * g_IS[2][ 60] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero +  92] != 0.){psi_MC_g_channel[zero +  92] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][ 10] * g_IS[1][ 10] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 61] * g_IS[2][ 61] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  93] != 0.){psi_MC_g_channel[zero +  93] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_f[no_ps][  5] * g_IS[1][  5] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][ 62] * g_IS[2][ 62] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero +  94] != 0.){psi_MC_g_channel[zero +  94] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_f[no_ps][ 15] * g_IS[1][ 15] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 63] * g_IS[2][ 63] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero +  95] != 0.){psi_MC_g_channel[zero +  95] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_f[no_ps][ 13] * g_IS[1][ 13] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 64] * g_IS[2][ 64] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  96] != 0.){psi_MC_g_channel[zero +  96] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_f[no_ps][  9] * g_IS[1][  9] * psi_g_t[no_ps][  6] * g_IS[2][  6] * psi_g_t[no_ps][ 65] * g_IS[2][ 65] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero +  97] != 0.){psi_MC_g_channel[zero +  97] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_f[no_ps][ 14] * g_IS[1][ 14] * psi_g_t[no_ps][ 22] * g_IS[2][ 22] * psi_g_t[no_ps][ 66] * g_IS[2][ 66] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  98] != 0.){psi_MC_g_channel[zero +  98] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_f[no_ps][ 12] * g_IS[1][ 12] * psi_g_t[no_ps][ 15] * g_IS[2][ 15] * psi_g_t[no_ps][ 67] * g_IS[2][ 67] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
      if (psi_MC_alpha[zero +  99] != 0.){psi_MC_g_channel[zero +  99] = psi_g_p[no_ps][  3] * g_IS[0][  3] * psi_g_p[no_ps][ 18] * g_IS[0][ 18] * psi_g_t[no_ps][ 68] * g_IS[2][ 68] * psi_g_d[no_ps][  0] * g_IS[3][  0] * psi_g_d[no_ps][  3] * g_IS[3][  3];}
      if (psi_MC_alpha[zero + 100] != 0.){psi_MC_g_channel[zero + 100] = psi_g_p[no_ps][  5] * g_IS[0][  5] * psi_g_p[no_ps][ 19] * g_IS[0][ 19] * psi_g_t[no_ps][ 69] * g_IS[2][ 69] * psi_g_d[no_ps][  1] * g_IS[3][  1] * psi_g_d[no_ps][  5] * g_IS[3][  5];}
      if (psi_MC_alpha[zero + 101] != 0.){psi_MC_g_channel[zero + 101] = psi_g_p[no_ps][  4] * g_IS[0][  4] * psi_g_p[no_ps][ 20] * g_IS[0][ 20] * psi_g_t[no_ps][ 70] * g_IS[2][ 70] * psi_g_d[no_ps][  2] * g_IS[3][  2] * psi_g_d[no_ps][  4] * g_IS[3][  4];}
      if (psi_MC_alpha[zero + 102] != 0.){psi_MC_g_channel[zero + 102] = psi_g_p[no_ps][  0] * g_IS[0][  0] * psi_g_p[no_ps][ 21] * g_IS[0][ 21] * psi_g_t[no_ps][ 71] * g_IS[2][ 71] * psi_g_d[no_ps][  3] * g_IS[3][  3] * psi_g_d[no_ps][  0] * g_IS[3][  0];}
      if (psi_MC_alpha[zero + 103] != 0.){psi_MC_g_channel[zero + 103] = psi_g_p[no_ps][  2] * g_IS[0][  2] * psi_g_p[no_ps][ 22] * g_IS[0][ 22] * psi_g_t[no_ps][ 72] * g_IS[2][ 72] * psi_g_d[no_ps][  4] * g_IS[3][  4] * psi_g_d[no_ps][  2] * g_IS[3][  2];}
      if (psi_MC_alpha[zero + 104] != 0.){psi_MC_g_channel[zero + 104] = psi_g_p[no_ps][  1] * g_IS[0][  1] * psi_g_p[no_ps][ 23] * g_IS[0][ 23] * psi_g_t[no_ps][ 73] * g_IS[2][ 73] * psi_g_d[no_ps][  5] * g_IS[3][  5] * psi_g_d[no_ps][  1] * g_IS[3][  1];}
    }
  }

  else {
    if (psi_MC_alpha[zero +   0] != 0.){psi_MC_g_channel[zero +   0] = psi_g_f[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_t[no_ps][  2];}
    if (psi_MC_alpha[zero +   1] != 0.){psi_MC_g_channel[zero +   1] = psi_g_f[no_ps][  2] * psi_g_f[no_ps][  3] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  3] * psi_g_t[no_ps][  4];}
    if (psi_MC_alpha[zero +   2] != 0.){psi_MC_g_channel[zero +   2] = psi_g_f[no_ps][  2] * psi_g_f[no_ps][  3] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  3] * psi_g_t[no_ps][  5];}
    if (psi_MC_alpha[zero +   3] != 0.){psi_MC_g_channel[zero +   3] = psi_g_f[no_ps][  4] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_t[no_ps][  8];}
    if (psi_MC_alpha[zero +   4] != 0.){psi_MC_g_channel[zero +   4] = psi_g_f[no_ps][  0] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_t[no_ps][  9];}
    if (psi_MC_alpha[zero +   5] != 0.){psi_MC_g_channel[zero +   5] = psi_g_f[no_ps][  6] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 10] * psi_g_t[no_ps][ 11];}
    if (psi_MC_alpha[zero +   6] != 0.){psi_MC_g_channel[zero +   6] = psi_g_f[no_ps][  6] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 10] * psi_g_t[no_ps][ 12];}
    if (psi_MC_alpha[zero +   7] != 0.){psi_MC_g_channel[zero +   7] = psi_g_f[no_ps][  8] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 13] * psi_g_t[no_ps][ 14];}
    if (psi_MC_alpha[zero +   8] != 0.){psi_MC_g_channel[zero +   8] = psi_g_f[no_ps][  4] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 16] * psi_g_t[no_ps][ 17];}
    if (psi_MC_alpha[zero +   9] != 0.){psi_MC_g_channel[zero +   9] = psi_g_f[no_ps][ 11] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 18] * psi_g_t[no_ps][ 19];}
    if (psi_MC_alpha[zero +  10] != 0.){psi_MC_g_channel[zero +  10] = psi_g_f[no_ps][ 11] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 18] * psi_g_t[no_ps][ 20];}
    if (psi_MC_alpha[zero +  11] != 0.){psi_MC_g_channel[zero +  11] = psi_g_f[no_ps][  8] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 13] * psi_g_t[no_ps][ 21];}
    if (psi_MC_alpha[zero +  12] != 0.){psi_MC_g_channel[zero +  12] = psi_g_f[no_ps][  8] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 23] * psi_g_t[no_ps][ 24];}
    if (psi_MC_alpha[zero +  13] != 0.){psi_MC_g_channel[zero +  13] = psi_g_f[no_ps][ 11] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 25] * psi_g_t[no_ps][ 26];}
    if (psi_MC_alpha[zero +  14] != 0.){psi_MC_g_channel[zero +  14] = psi_g_f[no_ps][ 11] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 25] * psi_g_t[no_ps][ 27];}
    if (psi_MC_alpha[zero +  15] != 0.){psi_MC_g_channel[zero +  15] = psi_g_f[no_ps][  4] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_t[no_ps][ 28];}
    if (psi_MC_alpha[zero +  16] != 0.){psi_MC_g_channel[zero +  16] = psi_g_f[no_ps][  8] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 23] * psi_g_t[no_ps][ 29];}
    if (psi_MC_alpha[zero +  17] != 0.){psi_MC_g_channel[zero +  17] = psi_g_f[no_ps][  6] * psi_g_f[no_ps][ 15] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 30] * psi_g_t[no_ps][ 31];}
    if (psi_MC_alpha[zero +  18] != 0.){psi_MC_g_channel[zero +  18] = psi_g_f[no_ps][  6] * psi_g_f[no_ps][ 15] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 30] * psi_g_t[no_ps][ 32];}
    if (psi_MC_alpha[zero +  19] != 0.){psi_MC_g_channel[zero +  19] = psi_g_f[no_ps][  0] * psi_g_f[no_ps][ 16] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 33] * psi_g_t[no_ps][ 34];}
    if (psi_MC_alpha[zero +  20] != 0.){psi_MC_g_channel[zero +  20] = psi_g_f[no_ps][  4] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 16] * psi_g_t[no_ps][ 35];}
    if (psi_MC_alpha[zero +  21] != 0.){psi_MC_g_channel[zero +  21] = psi_g_f[no_ps][  2] * psi_g_f[no_ps][ 17] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 36] * psi_g_t[no_ps][ 37];}
    if (psi_MC_alpha[zero +  22] != 0.){psi_MC_g_channel[zero +  22] = psi_g_f[no_ps][  2] * psi_g_f[no_ps][ 17] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 36] * psi_g_t[no_ps][ 38];}
    if (psi_MC_alpha[zero +  23] != 0.){psi_MC_g_channel[zero +  23] = psi_g_f[no_ps][  0] * psi_g_f[no_ps][ 16] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 33] * psi_g_t[no_ps][ 39];}
    if (psi_MC_alpha[zero +  24] != 0.){psi_MC_g_channel[zero +  24] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][ 40] * psi_g_t[no_ps][ 41] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  25] != 0.){psi_MC_g_channel[zero +  25] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 42] * psi_g_t[no_ps][ 43] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  26] != 0.){psi_MC_g_channel[zero +  26] = psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 42] * psi_g_t[no_ps][ 44] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  27] != 0.){psi_MC_g_channel[zero +  27] = psi_g_p[no_ps][  2] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][ 40] * psi_g_t[no_ps][ 45] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  28] != 0.){psi_MC_g_channel[zero +  28] = psi_g_p[no_ps][  2] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][ 46] * psi_g_t[no_ps][ 47] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  29] != 0.){psi_MC_g_channel[zero +  29] = psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][ 46] * psi_g_t[no_ps][ 48] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  30] != 0.){psi_MC_g_channel[zero +  30] = psi_g_p[no_ps][  3] * psi_g_f[no_ps][ 15] * psi_g_t[no_ps][ 46] * psi_g_t[no_ps][ 49] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  31] != 0.){psi_MC_g_channel[zero +  31] = psi_g_p[no_ps][  4] * psi_g_f[no_ps][ 17] * psi_g_t[no_ps][ 42] * psi_g_t[no_ps][ 50] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  32] != 0.){psi_MC_g_channel[zero +  32] = psi_g_p[no_ps][  5] * psi_g_f[no_ps][ 16] * psi_g_t[no_ps][ 40] * psi_g_t[no_ps][ 51] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  33] != 0.){psi_MC_g_channel[zero +  33] = psi_g_p[no_ps][  5] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][ 52] * psi_g_t[no_ps][ 53] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  34] != 0.){psi_MC_g_channel[zero +  34] = psi_g_p[no_ps][  4] * psi_g_f[no_ps][  3] * psi_g_t[no_ps][ 52] * psi_g_t[no_ps][ 54] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  35] != 0.){psi_MC_g_channel[zero +  35] = psi_g_p[no_ps][  3] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][ 52] * psi_g_t[no_ps][ 55] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  36] != 0.){psi_MC_g_channel[zero +  36] = psi_g_p[no_ps][  5] * psi_g_p[no_ps][  6] * psi_g_t[no_ps][ 40] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  37] != 0.){psi_MC_g_channel[zero +  37] = psi_g_p[no_ps][  4] * psi_g_p[no_ps][  7] * psi_g_t[no_ps][ 42] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  38] != 0.){psi_MC_g_channel[zero +  38] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  8] * psi_g_t[no_ps][ 42] * psi_g_d[no_ps][  8] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  39] != 0.){psi_MC_g_channel[zero +  39] = psi_g_p[no_ps][  3] * psi_g_p[no_ps][  9] * psi_g_t[no_ps][ 46] * psi_g_d[no_ps][  9] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  40] != 0.){psi_MC_g_channel[zero +  40] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 10] * psi_g_t[no_ps][ 46] * psi_g_d[no_ps][ 10] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  41] != 0.){psi_MC_g_channel[zero +  41] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][ 11] * psi_g_t[no_ps][ 40] * psi_g_d[no_ps][ 11] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  42] != 0.){psi_MC_g_channel[zero +  42] = psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 12] * psi_g_t[no_ps][ 46] * psi_g_d[no_ps][ 12] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  43] != 0.){psi_MC_g_channel[zero +  43] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 13] * psi_g_t[no_ps][ 40] * psi_g_d[no_ps][ 13] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  44] != 0.){psi_MC_g_channel[zero +  44] = psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 14] * psi_g_t[no_ps][ 42] * psi_g_d[no_ps][ 14] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  45] != 0.){psi_MC_g_channel[zero +  45] = psi_g_p[no_ps][  3] * psi_g_p[no_ps][ 15] * psi_g_t[no_ps][ 52] * psi_g_d[no_ps][ 15] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  46] != 0.){psi_MC_g_channel[zero +  46] = psi_g_p[no_ps][  5] * psi_g_p[no_ps][ 16] * psi_g_t[no_ps][ 52] * psi_g_d[no_ps][ 16] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  47] != 0.){psi_MC_g_channel[zero +  47] = psi_g_p[no_ps][  4] * psi_g_p[no_ps][ 17] * psi_g_t[no_ps][ 52] * psi_g_d[no_ps][ 17] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  48] != 0.){psi_MC_g_channel[zero +  48] = psi_g_p[no_ps][  3] * psi_g_p[no_ps][ 15] * psi_g_d[no_ps][ 18] * psi_g_d[no_ps][ 15] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  49] != 0.){psi_MC_g_channel[zero +  49] = psi_g_p[no_ps][  5] * psi_g_p[no_ps][ 16] * psi_g_d[no_ps][ 18] * psi_g_d[no_ps][ 16] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  50] != 0.){psi_MC_g_channel[zero +  50] = psi_g_p[no_ps][  5] * psi_g_p[no_ps][  6] * psi_g_d[no_ps][ 19] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  51] != 0.){psi_MC_g_channel[zero +  51] = psi_g_p[no_ps][  4] * psi_g_p[no_ps][ 17] * psi_g_d[no_ps][ 18] * psi_g_d[no_ps][ 17] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  52] != 0.){psi_MC_g_channel[zero +  52] = psi_g_p[no_ps][  4] * psi_g_p[no_ps][  7] * psi_g_d[no_ps][ 20] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  53] != 0.){psi_MC_g_channel[zero +  53] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  8] * psi_g_d[no_ps][ 20] * psi_g_d[no_ps][  8] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  54] != 0.){psi_MC_g_channel[zero +  54] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][ 11] * psi_g_d[no_ps][ 19] * psi_g_d[no_ps][ 11] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  55] != 0.){psi_MC_g_channel[zero +  55] = psi_g_p[no_ps][  3] * psi_g_p[no_ps][  9] * psi_g_d[no_ps][ 21] * psi_g_d[no_ps][  9] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  56] != 0.){psi_MC_g_channel[zero +  56] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 10] * psi_g_d[no_ps][ 21] * psi_g_d[no_ps][ 10] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  57] != 0.){psi_MC_g_channel[zero +  57] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 13] * psi_g_d[no_ps][ 19] * psi_g_d[no_ps][ 13] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  58] != 0.){psi_MC_g_channel[zero +  58] = psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 12] * psi_g_d[no_ps][ 21] * psi_g_d[no_ps][ 12] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  59] != 0.){psi_MC_g_channel[zero +  59] = psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 14] * psi_g_d[no_ps][ 20] * psi_g_d[no_ps][ 14] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  60] != 0.){psi_MC_g_channel[zero +  60] = psi_g_p[no_ps][  3] * psi_g_p[no_ps][ 18] * psi_g_d[no_ps][ 22] * psi_g_d[no_ps][  0] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  61] != 0.){psi_MC_g_channel[zero +  61] = psi_g_p[no_ps][  5] * psi_g_p[no_ps][ 19] * psi_g_d[no_ps][ 23] * psi_g_d[no_ps][  1] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  62] != 0.){psi_MC_g_channel[zero +  62] = psi_g_p[no_ps][  4] * psi_g_p[no_ps][ 20] * psi_g_d[no_ps][ 24] * psi_g_d[no_ps][  2] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  63] != 0.){psi_MC_g_channel[zero +  63] = psi_g_p[no_ps][  5] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  1] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  64] != 0.){psi_MC_g_channel[zero +  64] = psi_g_p[no_ps][  4] * psi_g_f[no_ps][  3] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][  3] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  65] != 0.){psi_MC_g_channel[zero +  65] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][  7] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  66] != 0.){psi_MC_g_channel[zero +  66] = psi_g_p[no_ps][  3] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 10] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  67] != 0.){psi_MC_g_channel[zero +  67] = psi_g_p[no_ps][  2] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 13] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  68] != 0.){psi_MC_g_channel[zero +  68] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 16] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  69] != 0.){psi_MC_g_channel[zero +  69] = psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 18] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  70] != 0.){psi_MC_g_channel[zero +  70] = psi_g_p[no_ps][  2] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 23] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  71] != 0.){psi_MC_g_channel[zero +  71] = psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 25] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  72] != 0.){psi_MC_g_channel[zero +  72] = psi_g_p[no_ps][  3] * psi_g_f[no_ps][ 15] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 30] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  73] != 0.){psi_MC_g_channel[zero +  73] = psi_g_p[no_ps][  5] * psi_g_f[no_ps][ 16] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 33] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  74] != 0.){psi_MC_g_channel[zero +  74] = psi_g_p[no_ps][  4] * psi_g_f[no_ps][ 17] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 36] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  75] != 0.){psi_MC_g_channel[zero +  75] = psi_g_p[no_ps][  3] * psi_g_p[no_ps][ 15] * psi_g_t[no_ps][  0] * psi_g_d[no_ps][ 15] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  76] != 0.){psi_MC_g_channel[zero +  76] = psi_g_p[no_ps][  5] * psi_g_p[no_ps][ 16] * psi_g_t[no_ps][  0] * psi_g_d[no_ps][ 16] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  77] != 0.){psi_MC_g_channel[zero +  77] = psi_g_p[no_ps][  5] * psi_g_p[no_ps][  6] * psi_g_t[no_ps][  6] * psi_g_d[no_ps][  6] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  78] != 0.){psi_MC_g_channel[zero +  78] = psi_g_p[no_ps][  4] * psi_g_p[no_ps][ 17] * psi_g_t[no_ps][  0] * psi_g_d[no_ps][ 17] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  79] != 0.){psi_MC_g_channel[zero +  79] = psi_g_p[no_ps][  4] * psi_g_p[no_ps][  7] * psi_g_t[no_ps][ 15] * psi_g_d[no_ps][  7] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  80] != 0.){psi_MC_g_channel[zero +  80] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][  8] * psi_g_t[no_ps][ 15] * psi_g_d[no_ps][  8] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  81] != 0.){psi_MC_g_channel[zero +  81] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][ 11] * psi_g_t[no_ps][  6] * psi_g_d[no_ps][ 11] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  82] != 0.){psi_MC_g_channel[zero +  82] = psi_g_p[no_ps][  3] * psi_g_p[no_ps][  9] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][  9] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  83] != 0.){psi_MC_g_channel[zero +  83] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 10] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][ 10] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  84] != 0.){psi_MC_g_channel[zero +  84] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 13] * psi_g_t[no_ps][  6] * psi_g_d[no_ps][ 13] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  85] != 0.){psi_MC_g_channel[zero +  85] = psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 12] * psi_g_t[no_ps][ 22] * psi_g_d[no_ps][ 12] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  86] != 0.){psi_MC_g_channel[zero +  86] = psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 14] * psi_g_t[no_ps][ 15] * psi_g_d[no_ps][ 14] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  87] != 0.){psi_MC_g_channel[zero +  87] = psi_g_p[no_ps][  3] * psi_g_f[no_ps][  7] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 56] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  88] != 0.){psi_MC_g_channel[zero +  88] = psi_g_p[no_ps][  5] * psi_g_f[no_ps][  1] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 57] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  89] != 0.){psi_MC_g_channel[zero +  89] = psi_g_p[no_ps][  5] * psi_g_f[no_ps][ 16] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 58] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero +  90] != 0.){psi_MC_g_channel[zero +  90] = psi_g_p[no_ps][  4] * psi_g_f[no_ps][  3] * psi_g_t[no_ps][  0] * psi_g_t[no_ps][ 59] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  91] != 0.){psi_MC_g_channel[zero +  91] = psi_g_p[no_ps][  4] * psi_g_f[no_ps][ 17] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 60] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero +  92] != 0.){psi_MC_g_channel[zero +  92] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][ 10] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 61] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  93] != 0.){psi_MC_g_channel[zero +  93] = psi_g_p[no_ps][  0] * psi_g_f[no_ps][  5] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 62] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero +  94] != 0.){psi_MC_g_channel[zero +  94] = psi_g_p[no_ps][  3] * psi_g_f[no_ps][ 15] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 63] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero +  95] != 0.){psi_MC_g_channel[zero +  95] = psi_g_p[no_ps][  2] * psi_g_f[no_ps][ 13] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 64] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  96] != 0.){psi_MC_g_channel[zero +  96] = psi_g_p[no_ps][  2] * psi_g_f[no_ps][  9] * psi_g_t[no_ps][  6] * psi_g_t[no_ps][ 65] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero +  97] != 0.){psi_MC_g_channel[zero +  97] = psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 14] * psi_g_t[no_ps][ 22] * psi_g_t[no_ps][ 66] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  98] != 0.){psi_MC_g_channel[zero +  98] = psi_g_p[no_ps][  1] * psi_g_f[no_ps][ 12] * psi_g_t[no_ps][ 15] * psi_g_t[no_ps][ 67] * psi_g_d[no_ps][  1];}
    if (psi_MC_alpha[zero +  99] != 0.){psi_MC_g_channel[zero +  99] = psi_g_p[no_ps][  3] * psi_g_p[no_ps][ 18] * psi_g_t[no_ps][ 68] * psi_g_d[no_ps][  0] * psi_g_d[no_ps][  3];}
    if (psi_MC_alpha[zero + 100] != 0.){psi_MC_g_channel[zero + 100] = psi_g_p[no_ps][  5] * psi_g_p[no_ps][ 19] * psi_g_t[no_ps][ 69] * psi_g_d[no_ps][  1] * psi_g_d[no_ps][  5];}
    if (psi_MC_alpha[zero + 101] != 0.){psi_MC_g_channel[zero + 101] = psi_g_p[no_ps][  4] * psi_g_p[no_ps][ 20] * psi_g_t[no_ps][ 70] * psi_g_d[no_ps][  2] * psi_g_d[no_ps][  4];}
    if (psi_MC_alpha[zero + 102] != 0.){psi_MC_g_channel[zero + 102] = psi_g_p[no_ps][  0] * psi_g_p[no_ps][ 21] * psi_g_t[no_ps][ 71] * psi_g_d[no_ps][  3] * psi_g_d[no_ps][  0];}
    if (psi_MC_alpha[zero + 103] != 0.){psi_MC_g_channel[zero + 103] = psi_g_p[no_ps][  2] * psi_g_p[no_ps][ 22] * psi_g_t[no_ps][ 72] * psi_g_d[no_ps][  4] * psi_g_d[no_ps][  2];}
    if (psi_MC_alpha[zero + 104] != 0.){psi_MC_g_channel[zero + 104] = psi_g_p[no_ps][  1] * psi_g_p[no_ps][ 23] * psi_g_t[no_ps][ 73] * psi_g_d[no_ps][  5] * psi_g_d[no_ps][  1];}
  }

/*
  if (psi_MC_alpha[   0] != 0.){ag_channel_1111(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[   1] != 0.){ag_channel_1111(g, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[   2] != 0.){ag_channel_1111(g, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[   3] != 0.){ag_channel_1111(g, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  if (psi_MC_alpha[   4] != 0.){ag_channel_1111(g, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[   5] != 0.){ag_channel_1111(g, 1, 2, 4, 16, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[   6] != 0.){ag_channel_1111(g, 1, 2, 4, 32, 16, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[   7] != 0.){ag_channel_1111(g, 1, 2, 32, 4, 16, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  if (psi_MC_alpha[   8] != 0.){ag_channel_1111(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  if (psi_MC_alpha[   9] != 0.){ag_channel_1111(g, 1, 2, 16, 4, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  if (psi_MC_alpha[  10] != 0.){ag_channel_1111(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  11] != 0.){ag_channel_1111(g, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  12] != 0.){ag_channel_1111(g, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  13] != 0.){ag_channel_1111(g, 1, 2, 8, 4, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  14] != 0.){ag_channel_1111(g, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  15] != 0.){ag_channel_1111(g, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  if (psi_MC_alpha[  16] != 0.){ag_channel_1111(g, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  17] != 0.){ag_channel_1111(g, 1, 2, 8, 16, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  18] != 0.){ag_channel_1111(g, 1, 2, 8, 32, 16, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  19] != 0.){ag_channel_1111(g, 1, 2, 32, 8, 16, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  if (psi_MC_alpha[  20] != 0.){ag_channel_1111(g, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   6, psi);}
  if (psi_MC_alpha[  21] != 0.){ag_channel_1111(g, 1, 2, 16, 8, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  if (psi_MC_alpha[  22] != 0.){ag_channel_1111(g, 1, 2, 16, 32, 8, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  23] != 0.){ag_channel_1111(g, 1, 2, 32, 16, 8, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  24] != 0.){ag_channel_2v11(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  25] != 0.){ag_channel_2v11(g, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  26] != 0.){ag_channel_2v11(g, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  27] != 0.){ag_channel_2v11(g, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  28] != 0.){ag_channel_2v11(g, 1, 2, 4, 16, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  29] != 0.){ag_channel_2v11(g, 1, 2, 4, 32, 16, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  30] != 0.){ag_channel_2v11(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  31] != 0.){ag_channel_2v11(g, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  32] != 0.){ag_channel_2v11(g, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  33] != 0.){ag_channel_2v11(g, 1, 2, 8, 16, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  34] != 0.){ag_channel_2v11(g, 1, 2, 8, 32, 16, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  35] != 0.){ag_channel_2v11(g, 1, 2, 16, 32, 8, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  36] != 0.){ag_channel_3yv1(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  37] != 0.){ag_channel_3yv1(g, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  38] != 0.){ag_channel_3yv1(g, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  39] != 0.){ag_channel_3yv1(g, 1, 2, 4, 16, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  if (psi_MC_alpha[  40] != 0.){ag_channel_3yv1(g, 1, 2, 32, 4, 16, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  41] != 0.){ag_channel_3yv1(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  42] != 0.){ag_channel_3yv1(g, 1, 2, 16, 4, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  43] != 0.){ag_channel_3yv1(g, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  44] != 0.){ag_channel_3yv1(g, 1, 2, 8, 4, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  45] != 0.){ag_channel_3yv1(g, 1, 2, 8, 16, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  if (psi_MC_alpha[  46] != 0.){ag_channel_3yv1(g, 1, 2, 32, 8, 16, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  47] != 0.){ag_channel_3yv1(g, 1, 2, 16, 8, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  48] != 0.){ag_channel_4yyv(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  49] != 0.){ag_channel_4yyv(g, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  50] != 0.){ag_channel_4yyv(g, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  51] != 0.){ag_channel_4yyv(g, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  52] != 0.){ag_channel_4yyv(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  53] != 0.){ag_channel_4yyv(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  54] != 0.){ag_channel_4yyv(g, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  55] != 0.){ag_channel_4yyv(g, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  56] != 0.){ag_channel_4yyv(g, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  57] != 0.){ag_channel_4yyv(g, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  58] != 0.){ag_channel_4yyv(g, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  59] != 0.){ag_channel_4yyv(g, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  60] != 0.){ag_channel_42vv(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  61] != 0.){ag_channel_42vv(g, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  62] != 0.){ag_channel_42vv(g, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  63] != 0.){ag_channel_12v1(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  64] != 0.){ag_channel_12v1(g, 1, 2, 4, 8, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  65] != 0.){ag_channel_12v1(g, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  66] != 0.){ag_channel_12v1(g, 1, 2, 4, 16, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  if (psi_MC_alpha[  67] != 0.){ag_channel_12v1(g, 1, 2, 32, 4, 16, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  68] != 0.){ag_channel_12v1(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  69] != 0.){ag_channel_12v1(g, 1, 2, 16, 4, 32, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  70] != 0.){ag_channel_12v1(g, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  71] != 0.){ag_channel_12v1(g, 1, 2, 8, 4, 32, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   0, psi);}
  if (psi_MC_alpha[  72] != 0.){ag_channel_12v1(g, 1, 2, 8, 16, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  if (psi_MC_alpha[  73] != 0.){ag_channel_12v1(g, 1, 2, 32, 8, 16, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  74] != 0.){ag_channel_12v1(g, 1, 2, 16, 8, 32, 4, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  75] != 0.){ag_channel_13yv(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  if (psi_MC_alpha[  76] != 0.){ag_channel_13yv(g, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  77] != 0.){ag_channel_13yv(g, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  78] != 0.){ag_channel_13yv(g, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  79] != 0.){ag_channel_13yv(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  80] != 0.){ag_channel_13yv(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  81] != 0.){ag_channel_13yv(g, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  82] != 0.){ag_channel_13yv(g, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   6, psi);}
  if (psi_MC_alpha[  83] != 0.){ag_channel_13yv(g, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  84] != 0.){ag_channel_13yv(g, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  85] != 0.){ag_channel_13yv(g, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  86] != 0.){ag_channel_13yv(g, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   0, psi);}
  if (psi_MC_alpha[  87] != 0.){ag_channel_112v(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  88] != 0.){ag_channel_112v(g, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  89] != 0.){ag_channel_112v(g, 1, 2, 32, 4, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  90] != 0.){ag_channel_112v(g, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  91] != 0.){ag_channel_112v(g, 1, 2, 16, 4, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  92] != 0.){ag_channel_112v(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  93] != 0.){ag_channel_112v(g, 1, 2, 32, 16, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[  94] != 0.){ag_channel_112v(g, 1, 2, 8, 4, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   6,   0, psi);}
  if (psi_MC_alpha[  95] != 0.){ag_channel_112v(g, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  96] != 0.){ag_channel_112v(g, 1, 2, 32, 8, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  97] != 0.){ag_channel_112v(g, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[  98] != 0.){ag_channel_112v(g, 1, 2, 16, 8, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   0,   6, psi);}
  if (psi_MC_alpha[  99] != 0.){ag_channel_2v2v(g, 1, 2, 4, 8, 16, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[ 100] != 0.){ag_channel_2v2v(g, 1, 2, 4, 32, 8, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[ 101] != 0.){ag_channel_2v2v(g, 1, 2, 4, 16, 8, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[ 102] != 0.){ag_channel_2v2v(g, 1, 2, 16, 32, 4, 8, p, s, sqrts, smin_opt, sqrtsmin_opt,   0,   0,   0, psi);}
  if (psi_MC_alpha[ 103] != 0.){ag_channel_2v2v(g, 1, 2, 8, 32, 4, 16, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  if (psi_MC_alpha[ 104] != 0.){ag_channel_2v2v(g, 1, 2, 8, 16, 4, 32, p, s, sqrts, smin_opt, sqrtsmin_opt,   6,   6,   6, psi);}
  // channel    0: 1111(1, 2, 4, 8, 16, 32;   6,   0,   0)
  // channel    1: 1111(1, 2, 4, 8, 32, 16;   6,   0,   0)
  // channel    2: 1111(1, 2, 4, 32, 8, 16;   6,   0,   6)
  // channel    3: 1111(1, 2, 32, 4, 8, 16;   0,   0,   6)
  // channel    4: 1111(1, 2, 4, 16, 8, 32;   6,   0,   6)
  // channel    5: 1111(1, 2, 4, 16, 32, 8;   6,   6,   6)
  // channel    6: 1111(1, 2, 4, 32, 16, 8;   6,   6,   6)
  // channel    7: 1111(1, 2, 32, 4, 16, 8;   0,   6,   6)
  // channel    8: 1111(1, 2, 16, 4, 8, 32;   0,   0,   6)
  // channel    9: 1111(1, 2, 16, 4, 32, 8;   0,   6,   6)
  // channel   10: 1111(1, 2, 16, 32, 4, 8;   0,   6,   0)
  // channel   11: 1111(1, 2, 32, 16, 4, 8;   0,   6,   0)
  // channel   12: 1111(1, 2, 8, 4, 16, 32;   6,   0,   0)
  // channel   13: 1111(1, 2, 8, 4, 32, 16;   6,   0,   0)
  // channel   14: 1111(1, 2, 8, 32, 4, 16;   6,   0,   6)
  // channel   15: 1111(1, 2, 32, 8, 4, 16;   0,   0,   6)
  // channel   16: 1111(1, 2, 8, 16, 4, 32;   6,   0,   6)
  // channel   17: 1111(1, 2, 8, 16, 32, 4;   6,   6,   6)
  // channel   18: 1111(1, 2, 8, 32, 16, 4;   6,   6,   6)
  // channel   19: 1111(1, 2, 32, 8, 16, 4;   0,   6,   6)
  // channel   20: 1111(1, 2, 16, 8, 4, 32;   0,   0,   6)
  // channel   21: 1111(1, 2, 16, 8, 32, 4;   0,   6,   6)
  // channel   22: 1111(1, 2, 16, 32, 8, 4;   0,   6,   0)
  // channel   23: 1111(1, 2, 32, 16, 8, 4;   0,   6,   0)
  // channel   24: 2v11(1, 2, 4, 8, 16, 32;   0,   0,   0)
  // channel   25: 2v11(1, 2, 4, 8, 32, 16;   0,   0,   0)
  // channel   26: 2v11(1, 2, 4, 32, 8, 16;   6,   0,   6)
  // channel   27: 2v11(1, 2, 4, 16, 8, 32;   6,   0,   6)
  // channel   28: 2v11(1, 2, 4, 16, 32, 8;   6,   6,   6)
  // channel   29: 2v11(1, 2, 4, 32, 16, 8;   6,   6,   6)
  // channel   30: 2v11(1, 2, 16, 32, 4, 8;   0,   6,   0)
  // channel   31: 2v11(1, 2, 8, 32, 4, 16;   6,   0,   6)
  // channel   32: 2v11(1, 2, 8, 16, 4, 32;   6,   0,   6)
  // channel   33: 2v11(1, 2, 8, 16, 32, 4;   6,   6,   6)
  // channel   34: 2v11(1, 2, 8, 32, 16, 4;   6,   6,   6)
  // channel   35: 2v11(1, 2, 16, 32, 8, 4;   0,   6,   0)
  // channel   36: 3yv1(1, 2, 4, 8, 16, 32;   6,   0,   0)
  // channel   37: 3yv1(1, 2, 4, 8, 32, 16;   6,   0,   0)
  // channel   38: 3yv1(1, 2, 32, 4, 8, 16;   0,   0,   0)
  // channel   39: 3yv1(1, 2, 4, 16, 32, 8;   0,   6,   6)
  // channel   40: 3yv1(1, 2, 32, 4, 16, 8;   6,   6,   6)
  // channel   41: 3yv1(1, 2, 16, 4, 8, 32;   0,   0,   0)
  // channel   42: 3yv1(1, 2, 16, 4, 32, 8;   6,   6,   6)
  // channel   43: 3yv1(1, 2, 8, 4, 16, 32;   6,   0,   0)
  // channel   44: 3yv1(1, 2, 8, 4, 32, 16;   6,   0,   0)
  // channel   45: 3yv1(1, 2, 8, 16, 32, 4;   0,   6,   6)
  // channel   46: 3yv1(1, 2, 32, 8, 16, 4;   6,   6,   6)
  // channel   47: 3yv1(1, 2, 16, 8, 32, 4;   6,   6,   6)
  // channel   48: 4yyv(1, 2, 4, 8, 16, 32;   0,   6,   0)
  // channel   49: 4yyv(1, 2, 4, 32, 8, 16;   6,   6,   0)
  // channel   50: 4yyv(1, 2, 32, 4, 8, 16;   6,   0,   0)
  // channel   51: 4yyv(1, 2, 4, 16, 8, 32;   6,   6,   0)
  // channel   52: 4yyv(1, 2, 16, 4, 8, 32;   6,   0,   0)
  // channel   53: 4yyv(1, 2, 16, 32, 4, 8;   0,   0,   0)
  // channel   54: 4yyv(1, 2, 32, 16, 4, 8;   0,   0,   0)
  // channel   55: 4yyv(1, 2, 8, 4, 16, 32;   0,   6,   0)
  // channel   56: 4yyv(1, 2, 8, 32, 4, 16;   6,   6,   0)
  // channel   57: 4yyv(1, 2, 32, 8, 4, 16;   6,   0,   0)
  // channel   58: 4yyv(1, 2, 8, 16, 4, 32;   6,   6,   0)
  // channel   59: 4yyv(1, 2, 16, 8, 4, 32;   6,   0,   0)
  // channel   60: 42vv(1, 2, 4, 8, 16, 32;   0,   0,   0)
  // channel   61: 42vv(1, 2, 4, 32, 8, 16;   6,   6,   0)
  // channel   62: 42vv(1, 2, 4, 16, 8, 32;   6,   6,   0)
  // channel   63: 12v1(1, 2, 4, 8, 16, 32;   6,   6,   0)
  // channel   64: 12v1(1, 2, 4, 8, 32, 16;   6,   6,   0)
  // channel   65: 12v1(1, 2, 32, 4, 8, 16;   0,   0,   0)
  // channel   66: 12v1(1, 2, 4, 16, 32, 8;   0,   6,   6)
  // channel   67: 12v1(1, 2, 32, 4, 16, 8;   6,   0,   6)
  // channel   68: 12v1(1, 2, 16, 4, 8, 32;   0,   0,   0)
  // channel   69: 12v1(1, 2, 16, 4, 32, 8;   6,   0,   6)
  // channel   70: 12v1(1, 2, 8, 4, 16, 32;   6,   6,   0)
  // channel   71: 12v1(1, 2, 8, 4, 32, 16;   6,   6,   0)
  // channel   72: 12v1(1, 2, 8, 16, 32, 4;   0,   6,   6)
  // channel   73: 12v1(1, 2, 32, 8, 16, 4;   6,   0,   6)
  // channel   74: 12v1(1, 2, 16, 8, 32, 4;   6,   0,   6)
  // channel   75: 13yv(1, 2, 4, 8, 16, 32;   0,   6,   6)
  // channel   76: 13yv(1, 2, 4, 32, 8, 16;   6,   6,   6)
  // channel   77: 13yv(1, 2, 32, 4, 8, 16;   6,   0,   0)
  // channel   78: 13yv(1, 2, 4, 16, 8, 32;   6,   6,   6)
  // channel   79: 13yv(1, 2, 16, 4, 8, 32;   6,   0,   0)
  // channel   80: 13yv(1, 2, 16, 32, 4, 8;   0,   0,   0)
  // channel   81: 13yv(1, 2, 32, 16, 4, 8;   0,   0,   0)
  // channel   82: 13yv(1, 2, 8, 4, 16, 32;   0,   6,   6)
  // channel   83: 13yv(1, 2, 8, 32, 4, 16;   6,   6,   6)
  // channel   84: 13yv(1, 2, 32, 8, 4, 16;   6,   0,   0)
  // channel   85: 13yv(1, 2, 8, 16, 4, 32;   6,   6,   6)
  // channel   86: 13yv(1, 2, 16, 8, 4, 32;   6,   0,   0)
  // channel   87: 112v(1, 2, 4, 8, 16, 32;   0,   6,   0)
  // channel   88: 112v(1, 2, 4, 32, 8, 16;   6,   6,   6)
  // channel   89: 112v(1, 2, 32, 4, 8, 16;   6,   0,   6)
  // channel   90: 112v(1, 2, 4, 16, 8, 32;   6,   6,   6)
  // channel   91: 112v(1, 2, 16, 4, 8, 32;   6,   0,   6)
  // channel   92: 112v(1, 2, 16, 32, 4, 8;   0,   0,   0)
  // channel   93: 112v(1, 2, 32, 16, 4, 8;   0,   0,   0)
  // channel   94: 112v(1, 2, 8, 4, 16, 32;   0,   6,   0)
  // channel   95: 112v(1, 2, 8, 32, 4, 16;   6,   6,   6)
  // channel   96: 112v(1, 2, 32, 8, 4, 16;   6,   0,   6)
  // channel   97: 112v(1, 2, 8, 16, 4, 32;   6,   6,   6)
  // channel   98: 112v(1, 2, 16, 8, 4, 32;   6,   0,   6)
  // channel   99: 2v2v(1, 2, 4, 8, 16, 32;   0,   0,   0)
  // channel  100: 2v2v(1, 2, 4, 32, 8, 16;   6,   6,   6)
  // channel  101: 2v2v(1, 2, 4, 16, 8, 32;   6,   6,   6)
  // channel  102: 2v2v(1, 2, 16, 32, 4, 8;   0,   0,   0)
  // channel  103: 2v2v(1, 2, 8, 32, 4, 16;   6,   6,   6)
  // channel  104: 2v2v(1, 2, 8, 16, 4, 32;   6,   6,   6)
*/

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
