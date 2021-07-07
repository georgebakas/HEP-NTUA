#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ax_psp_400_bxbx_ttxbxbx(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_400_bxbx_ttxbxbx");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  psi_v_smin[no_ps].resize(3);
  psi_needed_v_smin[no_ps].resize(3);
  int needed_v_smin_0[] = {   4,    5,    6,    7,    8,    9,   10,   11,   12,   13,    0,    1,    2,    3}; //  14
  psi_needed_v_smin[no_ps][  0] = get_vector_from_array_int(needed_v_smin_0, 14);
  int needed_v_smin_1[] = {   6,   10,    1,    3,    5,    9,   12}; //   7
  psi_needed_v_smin[no_ps][  1] = get_vector_from_array_int(needed_v_smin_1, 7);
  int needed_v_smin_2[] = {   7,   11,    0,    2,    4,    8,   13}; //   7
  psi_needed_v_smin[no_ps][  2] = get_vector_from_array_int(needed_v_smin_2, 7);

  psi_v_smax[no_ps].resize(3);
  psi_needed_v_smax[no_ps].resize(3);
  int needed_v_smax_0[] = {   4,    5,    6,    7,    8,    9,   10,   11,   12,   13,    0,    1,    2,    3}; //  14
  psi_needed_v_smax[no_ps][  0] = get_vector_from_array_int(needed_v_smax_0, 14);
  int needed_v_smax_1[] = {   6,   10,    1,    3,    5,    9,   12}; //   7
  psi_needed_v_smax[no_ps][  1] = get_vector_from_array_int(needed_v_smax_1, 7);
  int needed_v_smax_2[] = {   7,   11,    0,    2,    4,    8,   13}; //   7
  psi_needed_v_smax[no_ps][  2] = get_vector_from_array_int(needed_v_smax_2, 7);

  psi_container_IS_startvalue[no_ps][0] = psi_container_IS_name.size();
  psi_g_p[no_ps].resize(3);
  psi_needed_g_p[no_ps].resize(3);
  int needed_g_p_0[] = {   4,    5,    6,    7,    8,    9,   10,   11,   12,   13}; //  10
  psi_needed_g_p[no_ps][0] = get_vector_from_array_int(needed_g_p_0, 10);
  psi_container_IS_name.push_back("propagator_0_smin0_smax0");
  int needed_g_p_1[] = {   6,   10}; //   2
  psi_needed_g_p[no_ps][1] = get_vector_from_array_int(needed_g_p_1, 2);
  psi_container_IS_name.push_back("propagator_5_smin1_smax1");
  int needed_g_p_2[] = {   7,   11}; //   2
  psi_needed_g_p[no_ps][2] = get_vector_from_array_int(needed_g_p_2, 2);
  psi_container_IS_name.push_back("propagator_5_smin2_smax2");

  psi_container_IS_startvalue[no_ps][1] = psi_container_IS_name.size();
  psi_g_f[no_ps].resize(3);
  psi_needed_g_f[no_ps].resize(3);
  int needed_g_f_0[] = {   0,    1,    2,    3}; //   4
  psi_needed_g_f[no_ps][0] = get_vector_from_array_int(needed_g_f_0, 4);
  psi_container_IS_name.push_back("timelikeinvariant_smin0_smax0");
  int needed_g_f_1[] = {   0,    2,    4,    8,   13}; //   5
  psi_needed_g_f[no_ps][1] = get_vector_from_array_int(needed_g_f_1, 5);
  psi_container_IS_name.push_back("timelikeinvariant_smin2_smax2");
  int needed_g_f_2[] = {   1,    3,    5,    9,   12}; //   5
  psi_needed_g_f[no_ps][2] = get_vector_from_array_int(needed_g_f_2, 5);
  psi_container_IS_name.push_back("timelikeinvariant_smin1_smax1");

  psi_container_IS_startvalue[no_ps][2] = psi_container_IS_name.size();
  psi_g_t[no_ps].resize(14);
  psi_needed_g_t[no_ps].resize(14);
  int needed_g_t_0[] = {   0,    2,    8,   11,   13}; //   5
  psi_needed_g_t[no_ps][0] = get_vector_from_array_int(needed_g_t_0, 5);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__32_28");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__32_28");
  int needed_g_t_1[] = {   0,    2,    8}; //   3
  psi_needed_g_t[no_ps][1] = get_vector_from_array_int(needed_g_t_1, 3);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_33__16_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_33__16_12");
  int needed_g_t_2[] = {   0}; //   1
  psi_needed_g_t[no_ps][2] = get_vector_from_array_int(needed_g_t_2, 1);
  psi_container_IS_name.push_back("t_of_tchannel_6__33_18__4_8");
  psi_container_IS_name.push_back("phi_of_tchannel_6__33_18__4_8");
  int needed_g_t_3[] = {   1,    3,    9,   10,   12}; //   5
  psi_needed_g_t[no_ps][3] = get_vector_from_array_int(needed_g_t_3, 5);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__16_44");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__16_44");
  int needed_g_t_4[] = {   1,    3,    9}; //   3
  psi_needed_g_t[no_ps][4] = get_vector_from_array_int(needed_g_t_4, 3);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_17__32_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_17__32_12");
  int needed_g_t_5[] = {   1}; //   1
  psi_needed_g_t[no_ps][5] = get_vector_from_array_int(needed_g_t_5, 1);
  psi_container_IS_name.push_back("t_of_tchannel_6__17_34__4_8");
  psi_container_IS_name.push_back("phi_of_tchannel_6__17_34__4_8");
  int needed_g_t_6[] = {   2}; //   1
  psi_needed_g_t[no_ps][6] = get_vector_from_array_int(needed_g_t_6, 1);
  psi_container_IS_name.push_back("t_of_tchannel_6__33_18__8_4");
  psi_container_IS_name.push_back("phi_of_tchannel_6__33_18__8_4");
  int needed_g_t_7[] = {   3}; //   1
  psi_needed_g_t[no_ps][7] = get_vector_from_array_int(needed_g_t_7, 1);
  psi_container_IS_name.push_back("t_of_tchannel_6__17_34__8_4");
  psi_container_IS_name.push_back("phi_of_tchannel_6__17_34__8_4");
  int needed_g_t_8[] = {   4,    7}; //   2
  psi_needed_g_t[no_ps][8] = get_vector_from_array_int(needed_g_t_8, 2);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_1__32_28");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_1__32_28");
  int needed_g_t_9[] = {   4}; //   1
  psi_needed_g_t[no_ps][9] = get_vector_from_array_int(needed_g_t_9, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__1_34__12_16");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__1_34__12_16");
  int needed_g_t_10[] = {   5,    6}; //   2
  psi_needed_g_t[no_ps][10] = get_vector_from_array_int(needed_g_t_10, 2);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_1__16_44");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_1__16_44");
  int needed_g_t_11[] = {   5}; //   1
  psi_needed_g_t[no_ps][11] = get_vector_from_array_int(needed_g_t_11, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__1_18__12_32");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__1_18__12_32");
  int needed_g_t_12[] = {  12}; //   1
  psi_needed_g_t[no_ps][12] = get_vector_from_array_int(needed_g_t_12, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__2_17__12_32");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__2_17__12_32");
  int needed_g_t_13[] = {  13}; //   1
  psi_needed_g_t[no_ps][13] = get_vector_from_array_int(needed_g_t_13, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__2_33__12_16");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__2_33__12_16");

  psi_container_IS_startvalue[no_ps][3] = psi_container_IS_name.size();
  psi_g_d[no_ps].resize(3);
  psi_needed_g_d[no_ps].resize(3);
  int needed_g_d_0[] = {   4,    5,    6,    7,    8,    9,   10,   11,   12,   13}; //  10
  psi_needed_g_d[no_ps][0] = get_vector_from_array_int(needed_g_d_0, 10);
  psi_container_IS_name.push_back("costheta_of_decay_12__4_8");
  psi_container_IS_name.push_back("phi_of_decay_12__4_8");
  int needed_g_d_1[] = {   6,   10}; //   2
  psi_needed_g_d[no_ps][1] = get_vector_from_array_int(needed_g_d_1, 2);
  psi_container_IS_name.push_back("costheta_of_decay_44__32_12");
  psi_container_IS_name.push_back("phi_of_decay_44__32_12");
  int needed_g_d_2[] = {   7,   11}; //   2
  psi_needed_g_d[no_ps][2] = get_vector_from_array_int(needed_g_d_2, 2);
  psi_container_IS_name.push_back("costheta_of_decay_28__16_12");
  psi_container_IS_name.push_back("phi_of_decay_28__16_12");

  psi_c_p[no_ps].resize(14);
  psi_c_f[no_ps].resize(14);
  psi_c_t[no_ps].resize(14);
  psi_c_d[no_ps].resize(14);

  int c_p_0[] = {}; //   0
  psi_c_p[no_ps][0] = get_vector_from_array_int(c_p_0, 0);
  int c_f_0[] = {  0,   1}; //   2
  psi_c_f[no_ps][0] = get_vector_from_array_int(c_f_0, 2);
  int c_t_0[] = {  0,   1,   2}; //   3
  psi_c_t[no_ps][0] = get_vector_from_array_int(c_t_0, 3);
  int c_d_0[] = {}; //   0
  psi_c_d[no_ps][0] = get_vector_from_array_int(c_d_0, 0);

  int c_p_1[] = {}; //   0
  psi_c_p[no_ps][1] = get_vector_from_array_int(c_p_1, 0);
  int c_f_1[] = {  0,   2}; //   2
  psi_c_f[no_ps][1] = get_vector_from_array_int(c_f_1, 2);
  int c_t_1[] = {  3,   4,   5}; //   3
  psi_c_t[no_ps][1] = get_vector_from_array_int(c_t_1, 3);
  int c_d_1[] = {}; //   0
  psi_c_d[no_ps][1] = get_vector_from_array_int(c_d_1, 0);

  int c_p_2[] = {}; //   0
  psi_c_p[no_ps][2] = get_vector_from_array_int(c_p_2, 0);
  int c_f_2[] = {  0,   1}; //   2
  psi_c_f[no_ps][2] = get_vector_from_array_int(c_f_2, 2);
  int c_t_2[] = {  0,   1,   6}; //   3
  psi_c_t[no_ps][2] = get_vector_from_array_int(c_t_2, 3);
  int c_d_2[] = {}; //   0
  psi_c_d[no_ps][2] = get_vector_from_array_int(c_d_2, 0);

  int c_p_3[] = {}; //   0
  psi_c_p[no_ps][3] = get_vector_from_array_int(c_p_3, 0);
  int c_f_3[] = {  0,   2}; //   2
  psi_c_f[no_ps][3] = get_vector_from_array_int(c_f_3, 2);
  int c_t_3[] = {  3,   4,   7}; //   3
  psi_c_t[no_ps][3] = get_vector_from_array_int(c_t_3, 3);
  int c_d_3[] = {}; //   0
  psi_c_d[no_ps][3] = get_vector_from_array_int(c_d_3, 0);

  int c_p_4[] = {  0}; //   1
  psi_c_p[no_ps][4] = get_vector_from_array_int(c_p_4, 1);
  int c_f_4[] = {  1}; //   1
  psi_c_f[no_ps][4] = get_vector_from_array_int(c_f_4, 1);
  int c_t_4[] = {  8,   9}; //   2
  psi_c_t[no_ps][4] = get_vector_from_array_int(c_t_4, 2);
  int c_d_4[] = {  0}; //   1
  psi_c_d[no_ps][4] = get_vector_from_array_int(c_d_4, 1);

  int c_p_5[] = {  0}; //   1
  psi_c_p[no_ps][5] = get_vector_from_array_int(c_p_5, 1);
  int c_f_5[] = {  2}; //   1
  psi_c_f[no_ps][5] = get_vector_from_array_int(c_f_5, 1);
  int c_t_5[] = { 10,  11}; //   2
  psi_c_t[no_ps][5] = get_vector_from_array_int(c_t_5, 2);
  int c_d_5[] = {  0}; //   1
  psi_c_d[no_ps][5] = get_vector_from_array_int(c_d_5, 1);

  int c_p_6[] = {  0,   1}; //   2
  psi_c_p[no_ps][6] = get_vector_from_array_int(c_p_6, 2);
  int c_f_6[] = {}; //   0
  psi_c_f[no_ps][6] = get_vector_from_array_int(c_f_6, 0);
  int c_t_6[] = { 10}; //   1
  psi_c_t[no_ps][6] = get_vector_from_array_int(c_t_6, 1);
  int c_d_6[] = {  1,   0}; //   2
  psi_c_d[no_ps][6] = get_vector_from_array_int(c_d_6, 2);

  int c_p_7[] = {  0,   2}; //   2
  psi_c_p[no_ps][7] = get_vector_from_array_int(c_p_7, 2);
  int c_f_7[] = {}; //   0
  psi_c_f[no_ps][7] = get_vector_from_array_int(c_f_7, 0);
  int c_t_7[] = {  8}; //   1
  psi_c_t[no_ps][7] = get_vector_from_array_int(c_t_7, 1);
  int c_d_7[] = {  2,   0}; //   2
  psi_c_d[no_ps][7] = get_vector_from_array_int(c_d_7, 2);

  int c_p_8[] = {  0}; //   1
  psi_c_p[no_ps][8] = get_vector_from_array_int(c_p_8, 1);
  int c_f_8[] = {  1}; //   1
  psi_c_f[no_ps][8] = get_vector_from_array_int(c_f_8, 1);
  int c_t_8[] = {  0,   1}; //   2
  psi_c_t[no_ps][8] = get_vector_from_array_int(c_t_8, 2);
  int c_d_8[] = {  0}; //   1
  psi_c_d[no_ps][8] = get_vector_from_array_int(c_d_8, 1);

  int c_p_9[] = {  0}; //   1
  psi_c_p[no_ps][9] = get_vector_from_array_int(c_p_9, 1);
  int c_f_9[] = {  2}; //   1
  psi_c_f[no_ps][9] = get_vector_from_array_int(c_f_9, 1);
  int c_t_9[] = {  3,   4}; //   2
  psi_c_t[no_ps][9] = get_vector_from_array_int(c_t_9, 2);
  int c_d_9[] = {  0}; //   1
  psi_c_d[no_ps][9] = get_vector_from_array_int(c_d_9, 1);

  int c_p_10[] = {  0,   1}; //   2
  psi_c_p[no_ps][10] = get_vector_from_array_int(c_p_10, 2);
  int c_f_10[] = {}; //   0
  psi_c_f[no_ps][10] = get_vector_from_array_int(c_f_10, 0);
  int c_t_10[] = {  3}; //   1
  psi_c_t[no_ps][10] = get_vector_from_array_int(c_t_10, 1);
  int c_d_10[] = {  1,   0}; //   2
  psi_c_d[no_ps][10] = get_vector_from_array_int(c_d_10, 2);

  int c_p_11[] = {  0,   2}; //   2
  psi_c_p[no_ps][11] = get_vector_from_array_int(c_p_11, 2);
  int c_f_11[] = {}; //   0
  psi_c_f[no_ps][11] = get_vector_from_array_int(c_f_11, 0);
  int c_t_11[] = {  0}; //   1
  psi_c_t[no_ps][11] = get_vector_from_array_int(c_t_11, 1);
  int c_d_11[] = {  2,   0}; //   2
  psi_c_d[no_ps][11] = get_vector_from_array_int(c_d_11, 2);

  int c_p_12[] = {  0}; //   1
  psi_c_p[no_ps][12] = get_vector_from_array_int(c_p_12, 1);
  int c_f_12[] = {  2}; //   1
  psi_c_f[no_ps][12] = get_vector_from_array_int(c_f_12, 1);
  int c_t_12[] = {  3,  12}; //   2
  psi_c_t[no_ps][12] = get_vector_from_array_int(c_t_12, 2);
  int c_d_12[] = {  0}; //   1
  psi_c_d[no_ps][12] = get_vector_from_array_int(c_d_12, 1);

  int c_p_13[] = {  0}; //   1
  psi_c_p[no_ps][13] = get_vector_from_array_int(c_p_13, 1);
  int c_f_13[] = {  1}; //   1
  psi_c_f[no_ps][13] = get_vector_from_array_int(c_f_13, 1);
  int c_t_13[] = {  0,  13}; //   2
  psi_c_t[no_ps][13] = get_vector_from_array_int(c_t_13, 2);
  int c_d_13[] = {  0}; //   1
  psi_c_d[no_ps][13] = get_vector_from_array_int(c_d_13, 1);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
