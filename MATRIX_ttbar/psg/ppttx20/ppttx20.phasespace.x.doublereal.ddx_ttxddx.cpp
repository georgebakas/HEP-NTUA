#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ax_psp_400_ddx_ttxddx(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_400_ddx_ttxddx");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  psi_v_smin[no_ps].resize(6);
  psi_needed_v_smin[no_ps].resize(6);
  int needed_v_smin_0[] = {   2,    3,    5,    6,    8,    9,   10,   11,   12,   13,    0,    1}; //  12
  psi_needed_v_smin[no_ps][  0] = get_vector_from_array_int(needed_v_smin_0, 12);
  int needed_v_smin_1[] = {   3,    6,    2}; //   3
  psi_needed_v_smin[no_ps][  1] = get_vector_from_array_int(needed_v_smin_1, 3);
  int needed_v_smin_2[] = {   4,    7,    8,   12,   13}; //   5
  psi_needed_v_smin[no_ps][  2] = get_vector_from_array_int(needed_v_smin_2, 5);
  int needed_v_smin_3[] = {   4}; //   1
  psi_needed_v_smin[no_ps][  3] = get_vector_from_array_int(needed_v_smin_3, 1);
  int needed_v_smin_4[] = {   5,   10,    0,    1,    9,   11}; //   6
  psi_needed_v_smin[no_ps][  4] = get_vector_from_array_int(needed_v_smin_4, 6);
  int needed_v_smin_5[] = {   7}; //   1
  psi_needed_v_smin[no_ps][  5] = get_vector_from_array_int(needed_v_smin_5, 1);

  psi_v_smax[no_ps].resize(8);
  psi_needed_v_smax[no_ps].resize(8);
  int needed_v_smax_0[] = {   2,    3,    5,    6,    9,   10,   11,   13,    0,    1}; //  10
  psi_needed_v_smax[no_ps][  0] = get_vector_from_array_int(needed_v_smax_0, 10);
  int needed_v_smax_1[] = {   3,    6,    2}; //   3
  psi_needed_v_smax[no_ps][  1] = get_vector_from_array_int(needed_v_smax_1, 3);
  int needed_v_smax_2[] = {   4,    7,    8,   12}; //   4
  psi_needed_v_smax[no_ps][  2] = get_vector_from_array_int(needed_v_smax_2, 4);
  int needed_v_smax_3[] = {   4}; //   1
  psi_needed_v_smax[no_ps][  3] = get_vector_from_array_int(needed_v_smax_3, 1);
  int needed_v_smax_4[] = {   5,   10,    0,    1,    9,   11}; //   6
  psi_needed_v_smax[no_ps][  4] = get_vector_from_array_int(needed_v_smax_4, 6);
  int needed_v_smax_5[] = {   7}; //   1
  psi_needed_v_smax[no_ps][  5] = get_vector_from_array_int(needed_v_smax_5, 1);
  int needed_v_smax_6[] = {   8,   12}; //   2
  psi_needed_v_smax[no_ps][  6] = get_vector_from_array_int(needed_v_smax_6, 2);
  int needed_v_smax_7[] = {  13}; //   1
  psi_needed_v_smax[no_ps][  7] = get_vector_from_array_int(needed_v_smax_7, 1);

  psi_container_IS_startvalue[no_ps][0] = psi_container_IS_name.size();
  psi_g_p[no_ps].resize(8);
  psi_needed_g_p[no_ps].resize(8);
  int needed_g_p_0[] = {   2,    3,    5,    6,    9,   10,   11,   13}; //   8
  psi_needed_g_p[no_ps][0] = get_vector_from_array_int(needed_g_p_0, 8);
  psi_container_IS_name.push_back("propagator_0_smin0_smax0");
  int needed_g_p_1[] = {   3,    6}; //   2
  psi_needed_g_p[no_ps][1] = get_vector_from_array_int(needed_g_p_1, 2);
  psi_container_IS_name.push_back("propagator_0_smin1_smax1");
  int needed_g_p_2[] = {   4,    7,    8,   12}; //   4
  psi_needed_g_p[no_ps][2] = get_vector_from_array_int(needed_g_p_2, 4);
  psi_container_IS_name.push_back("propagator_0_smin2_smax2");
  int needed_g_p_3[] = {   4}; //   1
  psi_needed_g_p[no_ps][3] = get_vector_from_array_int(needed_g_p_3, 1);
  psi_container_IS_name.push_back("propagator_6_smin3_smax3");
  int needed_g_p_4[] = {   5,   10}; //   2
  psi_needed_g_p[no_ps][4] = get_vector_from_array_int(needed_g_p_4, 2);
  psi_container_IS_name.push_back("propagator_0_smin4_smax4");
  int needed_g_p_5[] = {   7}; //   1
  psi_needed_g_p[no_ps][5] = get_vector_from_array_int(needed_g_p_5, 1);
  psi_container_IS_name.push_back("propagator_6_smin5_smax5");
  int needed_g_p_6[] = {   8,   12}; //   2
  psi_needed_g_p[no_ps][6] = get_vector_from_array_int(needed_g_p_6, 2);
  psi_container_IS_name.push_back("propagator_0_smin0_smax6");
  int needed_g_p_7[] = {  13}; //   1
  psi_needed_g_p[no_ps][7] = get_vector_from_array_int(needed_g_p_7, 1);
  psi_container_IS_name.push_back("propagator_0_smin2_smax7");

  psi_container_IS_startvalue[no_ps][1] = psi_container_IS_name.size();
  psi_g_f[no_ps].resize(3);
  psi_needed_g_f[no_ps].resize(3);
  int needed_g_f_0[] = {   0,    1}; //   2
  psi_needed_g_f[no_ps][0] = get_vector_from_array_int(needed_g_f_0, 2);
  psi_container_IS_name.push_back("timelikeinvariant_smin0_smax0");
  int needed_g_f_1[] = {   0,    1,    9,   11}; //   4
  psi_needed_g_f[no_ps][1] = get_vector_from_array_int(needed_g_f_1, 4);
  psi_container_IS_name.push_back("timelikeinvariant_smin4_smax4");
  int needed_g_f_2[] = {   2}; //   1
  psi_needed_g_f[no_ps][2] = get_vector_from_array_int(needed_g_f_2, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin1_smax1");

  psi_container_IS_startvalue[no_ps][2] = psi_container_IS_name.size();
  psi_g_t[no_ps].resize(9);
  psi_needed_g_t[no_ps].resize(9);
  int needed_g_t_0[] = {   0,    1,    9,   10,   11}; //   5
  psi_needed_g_t[no_ps][0] = get_vector_from_array_int(needed_g_t_0, 5);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__16_44");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__16_44");
  int needed_g_t_1[] = {   0,    1,    9}; //   3
  psi_needed_g_t[no_ps][1] = get_vector_from_array_int(needed_g_t_1, 3);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_17__32_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_17__32_12");
  int needed_g_t_2[] = {   0}; //   1
  psi_needed_g_t[no_ps][2] = get_vector_from_array_int(needed_g_t_2, 1);
  psi_container_IS_name.push_back("t_of_tchannel_6__17_34__4_8");
  psi_container_IS_name.push_back("phi_of_tchannel_6__17_34__4_8");
  int needed_g_t_3[] = {   1}; //   1
  psi_needed_g_t[no_ps][3] = get_vector_from_array_int(needed_g_t_3, 1);
  psi_container_IS_name.push_back("t_of_tchannel_6__17_34__8_4");
  psi_container_IS_name.push_back("phi_of_tchannel_6__17_34__8_4");
  int needed_g_t_4[] = {   2,    3}; //   2
  psi_needed_g_t[no_ps][4] = get_vector_from_array_int(needed_g_t_4, 2);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_1__32_28");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_1__32_28");
  int needed_g_t_5[] = {   2}; //   1
  psi_needed_g_t[no_ps][5] = get_vector_from_array_int(needed_g_t_5, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_34__12_16");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_34__12_16");
  int needed_g_t_6[] = {  11}; //   1
  psi_needed_g_t[no_ps][6] = get_vector_from_array_int(needed_g_t_6, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_17__12_32");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_17__12_32");
  int needed_g_t_7[] = {  12}; //   1
  psi_needed_g_t[no_ps][7] = get_vector_from_array_int(needed_g_t_7, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__12_48");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__12_48");
  int needed_g_t_8[] = {  13}; //   1
  psi_needed_g_t[no_ps][8] = get_vector_from_array_int(needed_g_t_8, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__48_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__48_12");

  psi_container_IS_startvalue[no_ps][3] = psi_container_IS_name.size();
  psi_g_d[no_ps].resize(11);
  psi_needed_g_d[no_ps].resize(11);
  int needed_g_d_0[] = {   2,    3,    5,    6,    8,    9,   10,   11,   12,   13}; //  10
  psi_needed_g_d[no_ps][0] = get_vector_from_array_int(needed_g_d_0, 10);
  psi_container_IS_name.push_back("costheta_of_decay_12__4_8");
  psi_container_IS_name.push_back("phi_of_decay_12__4_8");
  int needed_g_d_1[] = {   3,    6}; //   2
  psi_needed_g_d[no_ps][1] = get_vector_from_array_int(needed_g_d_1, 2);
  psi_container_IS_name.push_back("costheta_of_decay_28__16_12");
  psi_container_IS_name.push_back("phi_of_decay_28__16_12");
  int needed_g_d_2[] = {   4}; //   1
  psi_needed_g_d[no_ps][2] = get_vector_from_array_int(needed_g_d_2, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__4_56");
  psi_container_IS_name.push_back("phi_of_decay_60__4_56");
  int needed_g_d_3[] = {   4}; //   1
  psi_needed_g_d[no_ps][3] = get_vector_from_array_int(needed_g_d_3, 1);
  psi_container_IS_name.push_back("costheta_of_decay_56__8_48");
  psi_container_IS_name.push_back("phi_of_decay_56__8_48");
  int needed_g_d_4[] = {   4,    7,    8,   12,   13}; //   5
  psi_needed_g_d[no_ps][4] = get_vector_from_array_int(needed_g_d_4, 5);
  psi_container_IS_name.push_back("costheta_of_decay_48__16_32");
  psi_container_IS_name.push_back("phi_of_decay_48__16_32");
  int needed_g_d_5[] = {   5}; //   1
  psi_needed_g_d[no_ps][5] = get_vector_from_array_int(needed_g_d_5, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__16_44");
  psi_container_IS_name.push_back("phi_of_decay_60__16_44");
  int needed_g_d_6[] = {   5,   10}; //   2
  psi_needed_g_d[no_ps][6] = get_vector_from_array_int(needed_g_d_6, 2);
  psi_container_IS_name.push_back("costheta_of_decay_44__32_12");
  psi_container_IS_name.push_back("phi_of_decay_44__32_12");
  int needed_g_d_7[] = {   6}; //   1
  psi_needed_g_d[no_ps][7] = get_vector_from_array_int(needed_g_d_7, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__32_28");
  psi_container_IS_name.push_back("phi_of_decay_60__32_28");
  int needed_g_d_8[] = {   7}; //   1
  psi_needed_g_d[no_ps][8] = get_vector_from_array_int(needed_g_d_8, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__8_52");
  psi_container_IS_name.push_back("phi_of_decay_60__8_52");
  int needed_g_d_9[] = {   7}; //   1
  psi_needed_g_d[no_ps][9] = get_vector_from_array_int(needed_g_d_9, 1);
  psi_container_IS_name.push_back("costheta_of_decay_52__4_48");
  psi_container_IS_name.push_back("phi_of_decay_52__4_48");
  int needed_g_d_10[] = {   8}; //   1
  psi_needed_g_d[no_ps][10] = get_vector_from_array_int(needed_g_d_10, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__12_48");
  psi_container_IS_name.push_back("phi_of_decay_60__12_48");

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
  int c_f_1[] = {  0,   1}; //   2
  psi_c_f[no_ps][1] = get_vector_from_array_int(c_f_1, 2);
  int c_t_1[] = {  0,   1,   3}; //   3
  psi_c_t[no_ps][1] = get_vector_from_array_int(c_t_1, 3);
  int c_d_1[] = {}; //   0
  psi_c_d[no_ps][1] = get_vector_from_array_int(c_d_1, 0);

  int c_p_2[] = {  0}; //   1
  psi_c_p[no_ps][2] = get_vector_from_array_int(c_p_2, 1);
  int c_f_2[] = {  2}; //   1
  psi_c_f[no_ps][2] = get_vector_from_array_int(c_f_2, 1);
  int c_t_2[] = {  4,   5}; //   2
  psi_c_t[no_ps][2] = get_vector_from_array_int(c_t_2, 2);
  int c_d_2[] = {  0}; //   1
  psi_c_d[no_ps][2] = get_vector_from_array_int(c_d_2, 1);

  int c_p_3[] = {  0,   1}; //   2
  psi_c_p[no_ps][3] = get_vector_from_array_int(c_p_3, 2);
  int c_f_3[] = {}; //   0
  psi_c_f[no_ps][3] = get_vector_from_array_int(c_f_3, 0);
  int c_t_3[] = {  4}; //   1
  psi_c_t[no_ps][3] = get_vector_from_array_int(c_t_3, 1);
  int c_d_3[] = {  1,   0}; //   2
  psi_c_d[no_ps][3] = get_vector_from_array_int(c_d_3, 2);

  int c_p_4[] = {  2,   3}; //   2
  psi_c_p[no_ps][4] = get_vector_from_array_int(c_p_4, 2);
  int c_f_4[] = {}; //   0
  psi_c_f[no_ps][4] = get_vector_from_array_int(c_f_4, 0);
  int c_t_4[] = {}; //   0
  psi_c_t[no_ps][4] = get_vector_from_array_int(c_t_4, 0);
  int c_d_4[] = {  2,   3,   4}; //   3
  psi_c_d[no_ps][4] = get_vector_from_array_int(c_d_4, 3);

  int c_p_5[] = {  0,   4}; //   2
  psi_c_p[no_ps][5] = get_vector_from_array_int(c_p_5, 2);
  int c_f_5[] = {}; //   0
  psi_c_f[no_ps][5] = get_vector_from_array_int(c_f_5, 0);
  int c_t_5[] = {}; //   0
  psi_c_t[no_ps][5] = get_vector_from_array_int(c_t_5, 0);
  int c_d_5[] = {  5,   6,   0}; //   3
  psi_c_d[no_ps][5] = get_vector_from_array_int(c_d_5, 3);

  int c_p_6[] = {  0,   1}; //   2
  psi_c_p[no_ps][6] = get_vector_from_array_int(c_p_6, 2);
  int c_f_6[] = {}; //   0
  psi_c_f[no_ps][6] = get_vector_from_array_int(c_f_6, 0);
  int c_t_6[] = {}; //   0
  psi_c_t[no_ps][6] = get_vector_from_array_int(c_t_6, 0);
  int c_d_6[] = {  7,   1,   0}; //   3
  psi_c_d[no_ps][6] = get_vector_from_array_int(c_d_6, 3);

  int c_p_7[] = {  2,   5}; //   2
  psi_c_p[no_ps][7] = get_vector_from_array_int(c_p_7, 2);
  int c_f_7[] = {}; //   0
  psi_c_f[no_ps][7] = get_vector_from_array_int(c_f_7, 0);
  int c_t_7[] = {}; //   0
  psi_c_t[no_ps][7] = get_vector_from_array_int(c_t_7, 0);
  int c_d_7[] = {  8,   9,   4}; //   3
  psi_c_d[no_ps][7] = get_vector_from_array_int(c_d_7, 3);

  int c_p_8[] = {  2,   6}; //   2
  psi_c_p[no_ps][8] = get_vector_from_array_int(c_p_8, 2);
  int c_f_8[] = {}; //   0
  psi_c_f[no_ps][8] = get_vector_from_array_int(c_f_8, 0);
  int c_t_8[] = {}; //   0
  psi_c_t[no_ps][8] = get_vector_from_array_int(c_t_8, 0);
  int c_d_8[] = { 10,   0,   4}; //   3
  psi_c_d[no_ps][8] = get_vector_from_array_int(c_d_8, 3);

  int c_p_9[] = {  0}; //   1
  psi_c_p[no_ps][9] = get_vector_from_array_int(c_p_9, 1);
  int c_f_9[] = {  1}; //   1
  psi_c_f[no_ps][9] = get_vector_from_array_int(c_f_9, 1);
  int c_t_9[] = {  0,   1}; //   2
  psi_c_t[no_ps][9] = get_vector_from_array_int(c_t_9, 2);
  int c_d_9[] = {  0}; //   1
  psi_c_d[no_ps][9] = get_vector_from_array_int(c_d_9, 1);

  int c_p_10[] = {  0,   4}; //   2
  psi_c_p[no_ps][10] = get_vector_from_array_int(c_p_10, 2);
  int c_f_10[] = {}; //   0
  psi_c_f[no_ps][10] = get_vector_from_array_int(c_f_10, 0);
  int c_t_10[] = {  0}; //   1
  psi_c_t[no_ps][10] = get_vector_from_array_int(c_t_10, 1);
  int c_d_10[] = {  6,   0}; //   2
  psi_c_d[no_ps][10] = get_vector_from_array_int(c_d_10, 2);

  int c_p_11[] = {  0}; //   1
  psi_c_p[no_ps][11] = get_vector_from_array_int(c_p_11, 1);
  int c_f_11[] = {  1}; //   1
  psi_c_f[no_ps][11] = get_vector_from_array_int(c_f_11, 1);
  int c_t_11[] = {  0,   6}; //   2
  psi_c_t[no_ps][11] = get_vector_from_array_int(c_t_11, 2);
  int c_d_11[] = {  0}; //   1
  psi_c_d[no_ps][11] = get_vector_from_array_int(c_d_11, 1);

  int c_p_12[] = {  2,   6}; //   2
  psi_c_p[no_ps][12] = get_vector_from_array_int(c_p_12, 2);
  int c_f_12[] = {}; //   0
  psi_c_f[no_ps][12] = get_vector_from_array_int(c_f_12, 0);
  int c_t_12[] = {  7}; //   1
  psi_c_t[no_ps][12] = get_vector_from_array_int(c_t_12, 1);
  int c_d_12[] = {  0,   4}; //   2
  psi_c_d[no_ps][12] = get_vector_from_array_int(c_d_12, 2);

  int c_p_13[] = {  0,   7}; //   2
  psi_c_p[no_ps][13] = get_vector_from_array_int(c_p_13, 2);
  int c_f_13[] = {}; //   0
  psi_c_f[no_ps][13] = get_vector_from_array_int(c_f_13, 0);
  int c_t_13[] = {  8}; //   1
  psi_c_t[no_ps][13] = get_vector_from_array_int(c_t_13, 1);
  int c_d_13[] = {  4,   0}; //   2
  psi_c_d[no_ps][13] = get_vector_from_array_int(c_d_13, 2);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
