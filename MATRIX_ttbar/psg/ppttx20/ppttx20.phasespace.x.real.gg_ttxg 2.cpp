#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ax_psp_300_gg_ttxg(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_300_gg_ttxg");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  psi_v_smin[no_ps].resize(3);
  psi_needed_v_smin[no_ps].resize(3);
  int needed_v_smin_0[] = {   6,   10,   13,    2,    5}; //   5
  psi_needed_v_smin[no_ps][  0] = get_vector_from_array_int(needed_v_smin_0, 5);
  int needed_v_smin_1[] = {   7,   11,   14,    3,    4}; //   5
  psi_needed_v_smin[no_ps][  1] = get_vector_from_array_int(needed_v_smin_1, 5);
  int needed_v_smin_2[] = {   8,    9,   12,    0,    1}; //   5
  psi_needed_v_smin[no_ps][  2] = get_vector_from_array_int(needed_v_smin_2, 5);

  psi_v_smax[no_ps].resize(3);
  psi_needed_v_smax[no_ps].resize(3);
  int needed_v_smax_0[] = {   6,   10,   13,    2,    5}; //   5
  psi_needed_v_smax[no_ps][  0] = get_vector_from_array_int(needed_v_smax_0, 5);
  int needed_v_smax_1[] = {   7,   11,   14,    3,    4}; //   5
  psi_needed_v_smax[no_ps][  1] = get_vector_from_array_int(needed_v_smax_1, 5);
  int needed_v_smax_2[] = {   8,    9,   12,    0,    1}; //   5
  psi_needed_v_smax[no_ps][  2] = get_vector_from_array_int(needed_v_smax_2, 5);

  psi_container_IS_startvalue[no_ps][0] = psi_container_IS_name.size();
  psi_g_p[no_ps].resize(3);
  psi_needed_g_p[no_ps].resize(3);
  int needed_g_p_0[] = {   6,   10,   13}; //   3
  psi_needed_g_p[no_ps][0] = get_vector_from_array_int(needed_g_p_0, 3);
  psi_container_IS_name.push_back("propagator_0_smin0_smax0");
  int needed_g_p_1[] = {   7,   11,   14}; //   3
  psi_needed_g_p[no_ps][1] = get_vector_from_array_int(needed_g_p_1, 3);
  psi_container_IS_name.push_back("propagator_6_smin1_smax1");
  int needed_g_p_2[] = {   8,    9,   12}; //   3
  psi_needed_g_p[no_ps][2] = get_vector_from_array_int(needed_g_p_2, 3);
  psi_container_IS_name.push_back("propagator_6_smin2_smax2");

  psi_container_IS_startvalue[no_ps][1] = psi_container_IS_name.size();
  psi_g_f[no_ps].resize(3);
  psi_needed_g_f[no_ps].resize(3);
  int needed_g_f_0[] = {   0,    1}; //   2
  psi_needed_g_f[no_ps][0] = get_vector_from_array_int(needed_g_f_0, 2);
  psi_container_IS_name.push_back("timelikeinvariant_smin2_smax2");
  int needed_g_f_1[] = {   2,    5}; //   2
  psi_needed_g_f[no_ps][1] = get_vector_from_array_int(needed_g_f_1, 2);
  psi_container_IS_name.push_back("timelikeinvariant_smin0_smax0");
  int needed_g_f_2[] = {   3,    4}; //   2
  psi_needed_g_f[no_ps][2] = get_vector_from_array_int(needed_g_f_2, 2);
  psi_container_IS_name.push_back("timelikeinvariant_smin1_smax1");

  psi_container_IS_startvalue[no_ps][2] = psi_container_IS_name.size();
  psi_g_t[no_ps].resize(12);
  psi_needed_g_t[no_ps].resize(12);
  int needed_g_t_0[] = {   0,    1,   12}; //   3
  psi_needed_g_t[no_ps][0] = get_vector_from_array_int(needed_g_t_0, 3);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__1_2__4_24");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__1_2__4_24");
  int needed_g_t_1[] = {   0}; //   1
  psi_needed_g_t[no_ps][1] = get_vector_from_array_int(needed_g_t_1, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_5__16_8");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_5__16_8");
  int needed_g_t_2[] = {   1}; //   1
  psi_needed_g_t[no_ps][2] = get_vector_from_array_int(needed_g_t_2, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_5__8_16");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_5__8_16");
  int needed_g_t_3[] = {   2,    5,   13}; //   3
  psi_needed_g_t[no_ps][3] = get_vector_from_array_int(needed_g_t_3, 3);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__16_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__16_12");
  int needed_g_t_4[] = {   2}; //   1
  psi_needed_g_t[no_ps][4] = get_vector_from_array_int(needed_g_t_4, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_17__8_4");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_17__8_4");
  int needed_g_t_5[] = {   3,    4,   14}; //   3
  psi_needed_g_t[no_ps][5] = get_vector_from_array_int(needed_g_t_5, 3);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__1_2__8_20");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__1_2__8_20");
  int needed_g_t_6[] = {   3}; //   1
  psi_needed_g_t[no_ps][6] = get_vector_from_array_int(needed_g_t_6, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_9__16_4");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_9__16_4");
  int needed_g_t_7[] = {   4}; //   1
  psi_needed_g_t[no_ps][7] = get_vector_from_array_int(needed_g_t_7, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_9__4_16");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_9__4_16");
  int needed_g_t_8[] = {   5}; //   1
  psi_needed_g_t[no_ps][8] = get_vector_from_array_int(needed_g_t_8, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_17__4_8");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_17__4_8");
  int needed_g_t_9[] = {   6}; //   1
  psi_needed_g_t[no_ps][9] = get_vector_from_array_int(needed_g_t_9, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_1__16_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_1__16_12");
  int needed_g_t_10[] = {   7}; //   1
  psi_needed_g_t[no_ps][10] = get_vector_from_array_int(needed_g_t_10, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_1__8_20");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_1__8_20");
  int needed_g_t_11[] = {   8}; //   1
  psi_needed_g_t[no_ps][11] = get_vector_from_array_int(needed_g_t_11, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_1__4_24");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_1__4_24");

  psi_container_IS_startvalue[no_ps][3] = psi_container_IS_name.size();
  psi_g_d[no_ps].resize(6);
  psi_needed_g_d[no_ps].resize(6);
  int needed_g_d_0[] = {   6,   10,   13}; //   3
  psi_needed_g_d[no_ps][0] = get_vector_from_array_int(needed_g_d_0, 3);
  psi_container_IS_name.push_back("costheta_of_decay_12__4_8");
  psi_container_IS_name.push_back("phi_of_decay_12__4_8");
  int needed_g_d_1[] = {   7,   11,   14}; //   3
  psi_needed_g_d[no_ps][1] = get_vector_from_array_int(needed_g_d_1, 3);
  psi_container_IS_name.push_back("costheta_of_decay_20__4_16");
  psi_container_IS_name.push_back("phi_of_decay_20__4_16");
  int needed_g_d_2[] = {   8,    9,   12}; //   3
  psi_needed_g_d[no_ps][2] = get_vector_from_array_int(needed_g_d_2, 3);
  psi_container_IS_name.push_back("costheta_of_decay_24__8_16");
  psi_container_IS_name.push_back("phi_of_decay_24__8_16");
  int needed_g_d_3[] = {   9}; //   1
  psi_needed_g_d[no_ps][3] = get_vector_from_array_int(needed_g_d_3, 1);
  psi_container_IS_name.push_back("costheta_of_decay_28__4_24");
  psi_container_IS_name.push_back("phi_of_decay_28__4_24");
  int needed_g_d_4[] = {  10}; //   1
  psi_needed_g_d[no_ps][4] = get_vector_from_array_int(needed_g_d_4, 1);
  psi_container_IS_name.push_back("costheta_of_decay_28__16_12");
  psi_container_IS_name.push_back("phi_of_decay_28__16_12");
  int needed_g_d_5[] = {  11}; //   1
  psi_needed_g_d[no_ps][5] = get_vector_from_array_int(needed_g_d_5, 1);
  psi_container_IS_name.push_back("costheta_of_decay_28__8_20");
  psi_container_IS_name.push_back("phi_of_decay_28__8_20");

  psi_c_p[no_ps].resize(15);
  psi_c_f[no_ps].resize(15);
  psi_c_t[no_ps].resize(15);
  psi_c_d[no_ps].resize(15);

  int c_p_0[] = {}; //   0
  psi_c_p[no_ps][0] = get_vector_from_array_int(c_p_0, 0);
  int c_f_0[] = {  0}; //   1
  psi_c_f[no_ps][0] = get_vector_from_array_int(c_f_0, 1);
  int c_t_0[] = {  0,   1}; //   2
  psi_c_t[no_ps][0] = get_vector_from_array_int(c_t_0, 2);
  int c_d_0[] = {}; //   0
  psi_c_d[no_ps][0] = get_vector_from_array_int(c_d_0, 0);

  int c_p_1[] = {}; //   0
  psi_c_p[no_ps][1] = get_vector_from_array_int(c_p_1, 0);
  int c_f_1[] = {  0}; //   1
  psi_c_f[no_ps][1] = get_vector_from_array_int(c_f_1, 1);
  int c_t_1[] = {  0,   2}; //   2
  psi_c_t[no_ps][1] = get_vector_from_array_int(c_t_1, 2);
  int c_d_1[] = {}; //   0
  psi_c_d[no_ps][1] = get_vector_from_array_int(c_d_1, 0);

  int c_p_2[] = {}; //   0
  psi_c_p[no_ps][2] = get_vector_from_array_int(c_p_2, 0);
  int c_f_2[] = {  1}; //   1
  psi_c_f[no_ps][2] = get_vector_from_array_int(c_f_2, 1);
  int c_t_2[] = {  3,   4}; //   2
  psi_c_t[no_ps][2] = get_vector_from_array_int(c_t_2, 2);
  int c_d_2[] = {}; //   0
  psi_c_d[no_ps][2] = get_vector_from_array_int(c_d_2, 0);

  int c_p_3[] = {}; //   0
  psi_c_p[no_ps][3] = get_vector_from_array_int(c_p_3, 0);
  int c_f_3[] = {  2}; //   1
  psi_c_f[no_ps][3] = get_vector_from_array_int(c_f_3, 1);
  int c_t_3[] = {  5,   6}; //   2
  psi_c_t[no_ps][3] = get_vector_from_array_int(c_t_3, 2);
  int c_d_3[] = {}; //   0
  psi_c_d[no_ps][3] = get_vector_from_array_int(c_d_3, 0);

  int c_p_4[] = {}; //   0
  psi_c_p[no_ps][4] = get_vector_from_array_int(c_p_4, 0);
  int c_f_4[] = {  2}; //   1
  psi_c_f[no_ps][4] = get_vector_from_array_int(c_f_4, 1);
  int c_t_4[] = {  5,   7}; //   2
  psi_c_t[no_ps][4] = get_vector_from_array_int(c_t_4, 2);
  int c_d_4[] = {}; //   0
  psi_c_d[no_ps][4] = get_vector_from_array_int(c_d_4, 0);

  int c_p_5[] = {}; //   0
  psi_c_p[no_ps][5] = get_vector_from_array_int(c_p_5, 0);
  int c_f_5[] = {  1}; //   1
  psi_c_f[no_ps][5] = get_vector_from_array_int(c_f_5, 1);
  int c_t_5[] = {  3,   8}; //   2
  psi_c_t[no_ps][5] = get_vector_from_array_int(c_t_5, 2);
  int c_d_5[] = {}; //   0
  psi_c_d[no_ps][5] = get_vector_from_array_int(c_d_5, 0);

  int c_p_6[] = {  0}; //   1
  psi_c_p[no_ps][6] = get_vector_from_array_int(c_p_6, 1);
  int c_f_6[] = {}; //   0
  psi_c_f[no_ps][6] = get_vector_from_array_int(c_f_6, 0);
  int c_t_6[] = {  9}; //   1
  psi_c_t[no_ps][6] = get_vector_from_array_int(c_t_6, 1);
  int c_d_6[] = {  0}; //   1
  psi_c_d[no_ps][6] = get_vector_from_array_int(c_d_6, 1);

  int c_p_7[] = {  1}; //   1
  psi_c_p[no_ps][7] = get_vector_from_array_int(c_p_7, 1);
  int c_f_7[] = {}; //   0
  psi_c_f[no_ps][7] = get_vector_from_array_int(c_f_7, 0);
  int c_t_7[] = { 10}; //   1
  psi_c_t[no_ps][7] = get_vector_from_array_int(c_t_7, 1);
  int c_d_7[] = {  1}; //   1
  psi_c_d[no_ps][7] = get_vector_from_array_int(c_d_7, 1);

  int c_p_8[] = {  2}; //   1
  psi_c_p[no_ps][8] = get_vector_from_array_int(c_p_8, 1);
  int c_f_8[] = {}; //   0
  psi_c_f[no_ps][8] = get_vector_from_array_int(c_f_8, 0);
  int c_t_8[] = { 11}; //   1
  psi_c_t[no_ps][8] = get_vector_from_array_int(c_t_8, 1);
  int c_d_8[] = {  2}; //   1
  psi_c_d[no_ps][8] = get_vector_from_array_int(c_d_8, 1);

  int c_p_9[] = {  2}; //   1
  psi_c_p[no_ps][9] = get_vector_from_array_int(c_p_9, 1);
  int c_f_9[] = {}; //   0
  psi_c_f[no_ps][9] = get_vector_from_array_int(c_f_9, 0);
  int c_t_9[] = {}; //   0
  psi_c_t[no_ps][9] = get_vector_from_array_int(c_t_9, 0);
  int c_d_9[] = {  3,   2}; //   2
  psi_c_d[no_ps][9] = get_vector_from_array_int(c_d_9, 2);

  int c_p_10[] = {  0}; //   1
  psi_c_p[no_ps][10] = get_vector_from_array_int(c_p_10, 1);
  int c_f_10[] = {}; //   0
  psi_c_f[no_ps][10] = get_vector_from_array_int(c_f_10, 0);
  int c_t_10[] = {}; //   0
  psi_c_t[no_ps][10] = get_vector_from_array_int(c_t_10, 0);
  int c_d_10[] = {  4,   0}; //   2
  psi_c_d[no_ps][10] = get_vector_from_array_int(c_d_10, 2);

  int c_p_11[] = {  1}; //   1
  psi_c_p[no_ps][11] = get_vector_from_array_int(c_p_11, 1);
  int c_f_11[] = {}; //   0
  psi_c_f[no_ps][11] = get_vector_from_array_int(c_f_11, 0);
  int c_t_11[] = {}; //   0
  psi_c_t[no_ps][11] = get_vector_from_array_int(c_t_11, 0);
  int c_d_11[] = {  5,   1}; //   2
  psi_c_d[no_ps][11] = get_vector_from_array_int(c_d_11, 2);

  int c_p_12[] = {  2}; //   1
  psi_c_p[no_ps][12] = get_vector_from_array_int(c_p_12, 1);
  int c_f_12[] = {}; //   0
  psi_c_f[no_ps][12] = get_vector_from_array_int(c_f_12, 0);
  int c_t_12[] = {  0}; //   1
  psi_c_t[no_ps][12] = get_vector_from_array_int(c_t_12, 1);
  int c_d_12[] = {  2}; //   1
  psi_c_d[no_ps][12] = get_vector_from_array_int(c_d_12, 1);

  int c_p_13[] = {  0}; //   1
  psi_c_p[no_ps][13] = get_vector_from_array_int(c_p_13, 1);
  int c_f_13[] = {}; //   0
  psi_c_f[no_ps][13] = get_vector_from_array_int(c_f_13, 0);
  int c_t_13[] = {  3}; //   1
  psi_c_t[no_ps][13] = get_vector_from_array_int(c_t_13, 1);
  int c_d_13[] = {  0}; //   1
  psi_c_d[no_ps][13] = get_vector_from_array_int(c_d_13, 1);

  int c_p_14[] = {  1}; //   1
  psi_c_p[no_ps][14] = get_vector_from_array_int(c_p_14, 1);
  int c_f_14[] = {}; //   0
  psi_c_f[no_ps][14] = get_vector_from_array_int(c_f_14, 0);
  int c_t_14[] = {  5}; //   1
  psi_c_t[no_ps][14] = get_vector_from_array_int(c_t_14, 1);
  int c_d_14[] = {  1}; //   1
  psi_c_d[no_ps][14] = get_vector_from_array_int(c_d_14, 1);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
