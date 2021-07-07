#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ax_psp_400_gg_ttxuux(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_400_gg_ttxuux");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  psi_v_smin[no_ps].resize(18);
  psi_needed_v_smin[no_ps].resize(18);
  int needed_v_smin_0[] = {   8,    9,   12,   14,   17,   18,   20,   21,   23,   26,   27,   30,   31,   33,   34}; //  15
  psi_needed_v_smin[no_ps][  0] = get_vector_from_array_int(needed_v_smin_0, 15);
  int needed_v_smin_1[] = {  10,   11,   13,   15,   16,   19,   20,   22,   24,   25,   28,   29,   32,   33,   34}; //  15
  psi_needed_v_smin[no_ps][  1] = get_vector_from_array_int(needed_v_smin_1, 15);
  int needed_v_smin_2[] = {  12,   17,   26,    9,   23,   30}; //   6
  psi_needed_v_smin[no_ps][  2] = get_vector_from_array_int(needed_v_smin_2, 6);
  int needed_v_smin_3[] = {  13,   19,   28,   10,   24,   32}; //   6
  psi_needed_v_smin[no_ps][  3] = get_vector_from_array_int(needed_v_smin_3, 6);
  int needed_v_smin_4[] = {  14,   18,   27,    8,   21,   31}; //   6
  psi_needed_v_smin[no_ps][  4] = get_vector_from_array_int(needed_v_smin_4, 6);
  int needed_v_smin_5[] = {  15,   16,   25,   11,   22,   29}; //   6
  psi_needed_v_smin[no_ps][  5] = get_vector_from_array_int(needed_v_smin_5, 6);
  int needed_v_smin_6[] = {   0,    7}; //   2
  psi_needed_v_smin[no_ps][  6] = get_vector_from_array_int(needed_v_smin_6, 2);
  int needed_v_smin_7[] = {   0}; //   1
  psi_needed_v_smin[no_ps][  7] = get_vector_from_array_int(needed_v_smin_7, 1);
  int needed_v_smin_8[] = {   1,    6}; //   2
  psi_needed_v_smin[no_ps][  8] = get_vector_from_array_int(needed_v_smin_8, 2);
  int needed_v_smin_9[] = {   1}; //   1
  psi_needed_v_smin[no_ps][  9] = get_vector_from_array_int(needed_v_smin_9, 1);
  int needed_v_smin_10[] = {   2,    5}; //   2
  psi_needed_v_smin[no_ps][ 10] = get_vector_from_array_int(needed_v_smin_10, 2);
  int needed_v_smin_11[] = {   2}; //   1
  psi_needed_v_smin[no_ps][ 11] = get_vector_from_array_int(needed_v_smin_11, 1);
  int needed_v_smin_12[] = {   3,    4}; //   2
  psi_needed_v_smin[no_ps][ 12] = get_vector_from_array_int(needed_v_smin_12, 2);
  int needed_v_smin_13[] = {   3}; //   1
  psi_needed_v_smin[no_ps][ 13] = get_vector_from_array_int(needed_v_smin_13, 1);
  int needed_v_smin_14[] = {   4}; //   1
  psi_needed_v_smin[no_ps][ 14] = get_vector_from_array_int(needed_v_smin_14, 1);
  int needed_v_smin_15[] = {   5}; //   1
  psi_needed_v_smin[no_ps][ 15] = get_vector_from_array_int(needed_v_smin_15, 1);
  int needed_v_smin_16[] = {   6}; //   1
  psi_needed_v_smin[no_ps][ 16] = get_vector_from_array_int(needed_v_smin_16, 1);
  int needed_v_smin_17[] = {   7}; //   1
  psi_needed_v_smin[no_ps][ 17] = get_vector_from_array_int(needed_v_smin_17, 1);

  psi_v_smax[no_ps].resize(12);
  psi_needed_v_smax[no_ps].resize(12);
  int needed_v_smax_0[] = {   8,    9,   12,   14,   17,   18,   21,   23,   26,   27,   30,   31,   34}; //  13
  psi_needed_v_smax[no_ps][  0] = get_vector_from_array_int(needed_v_smax_0, 13);
  int needed_v_smax_1[] = {  10,   11,   13,   15,   16,   19,   20,   22,   24,   25,   28,   29,   32,   33}; //  14
  psi_needed_v_smax[no_ps][  1] = get_vector_from_array_int(needed_v_smax_1, 14);
  int needed_v_smax_2[] = {  12,   17,   26,    2,    6,    9,   23,   30}; //   8
  psi_needed_v_smax[no_ps][  2] = get_vector_from_array_int(needed_v_smax_2, 8);
  int needed_v_smax_3[] = {  13,   19,   28,    4,    5,   10,   24,   32}; //   8
  psi_needed_v_smax[no_ps][  3] = get_vector_from_array_int(needed_v_smax_3, 8);
  int needed_v_smax_4[] = {  14,   18,   27,    3,    7,    8,   21,   31}; //   8
  psi_needed_v_smax[no_ps][  4] = get_vector_from_array_int(needed_v_smax_4, 8);
  int needed_v_smax_5[] = {  15,   16,   25,    0,    1,   11,   22,   29}; //   8
  psi_needed_v_smax[no_ps][  5] = get_vector_from_array_int(needed_v_smax_5, 8);
  int needed_v_smax_6[] = {  20,   33}; //   2
  psi_needed_v_smax[no_ps][  6] = get_vector_from_array_int(needed_v_smax_6, 2);
  int needed_v_smax_7[] = {  34}; //   1
  psi_needed_v_smax[no_ps][  7] = get_vector_from_array_int(needed_v_smax_7, 1);
  int needed_v_smax_8[] = {   0,    7}; //   2
  psi_needed_v_smax[no_ps][  8] = get_vector_from_array_int(needed_v_smax_8, 2);
  int needed_v_smax_9[] = {   1,    6}; //   2
  psi_needed_v_smax[no_ps][  9] = get_vector_from_array_int(needed_v_smax_9, 2);
  int needed_v_smax_10[] = {   2,    5}; //   2
  psi_needed_v_smax[no_ps][ 10] = get_vector_from_array_int(needed_v_smax_10, 2);
  int needed_v_smax_11[] = {   3,    4}; //   2
  psi_needed_v_smax[no_ps][ 11] = get_vector_from_array_int(needed_v_smax_11, 2);

  psi_container_IS_startvalue[no_ps][0] = psi_container_IS_name.size();
  psi_g_p[no_ps].resize(8);
  psi_needed_g_p[no_ps].resize(8);
  int needed_g_p_0[] = {   8,    9,   12,   14,   17,   18,   21,   23,   26,   27,   30,   31,   34}; //  13
  psi_needed_g_p[no_ps][0] = get_vector_from_array_int(needed_g_p_0, 13);
  psi_container_IS_name.push_back("propagator_0_smin0_smax0");
  int needed_g_p_1[] = {  10,   11,   13,   15,   16,   19,   20,   22,   24,   25,   28,   29,   32,   33}; //  14
  psi_needed_g_p[no_ps][1] = get_vector_from_array_int(needed_g_p_1, 14);
  psi_container_IS_name.push_back("propagator_0_smin1_smax1");
  int needed_g_p_2[] = {  12,   17,   26}; //   3
  psi_needed_g_p[no_ps][2] = get_vector_from_array_int(needed_g_p_2, 3);
  psi_container_IS_name.push_back("propagator_0_smin2_smax2");
  int needed_g_p_3[] = {  13,   19,   28}; //   3
  psi_needed_g_p[no_ps][3] = get_vector_from_array_int(needed_g_p_3, 3);
  psi_container_IS_name.push_back("propagator_6_smin3_smax3");
  int needed_g_p_4[] = {  14,   18,   27}; //   3
  psi_needed_g_p[no_ps][4] = get_vector_from_array_int(needed_g_p_4, 3);
  psi_container_IS_name.push_back("propagator_0_smin4_smax4");
  int needed_g_p_5[] = {  15,   16,   25}; //   3
  psi_needed_g_p[no_ps][5] = get_vector_from_array_int(needed_g_p_5, 3);
  psi_container_IS_name.push_back("propagator_6_smin5_smax5");
  int needed_g_p_6[] = {  20,   33}; //   2
  psi_needed_g_p[no_ps][6] = get_vector_from_array_int(needed_g_p_6, 2);
  psi_container_IS_name.push_back("propagator_0_smin0_smax6");
  int needed_g_p_7[] = {  34}; //   1
  psi_needed_g_p[no_ps][7] = get_vector_from_array_int(needed_g_p_7, 1);
  psi_container_IS_name.push_back("propagator_0_smin1_smax7");

  psi_container_IS_startvalue[no_ps][1] = psi_container_IS_name.size();
  psi_g_f[no_ps].resize(16);
  psi_needed_g_f[no_ps].resize(16);
  int needed_g_f_0[] = {   0,    7}; //   2
  psi_needed_g_f[no_ps][0] = get_vector_from_array_int(needed_g_f_0, 2);
  psi_container_IS_name.push_back("timelikeinvariant_smin6_smax8");
  int needed_g_f_1[] = {   0}; //   1
  psi_needed_g_f[no_ps][1] = get_vector_from_array_int(needed_g_f_1, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin7_smax5");
  int needed_g_f_2[] = {   1,    6}; //   2
  psi_needed_g_f[no_ps][2] = get_vector_from_array_int(needed_g_f_2, 2);
  psi_container_IS_name.push_back("timelikeinvariant_smin8_smax9");
  int needed_g_f_3[] = {   1}; //   1
  psi_needed_g_f[no_ps][3] = get_vector_from_array_int(needed_g_f_3, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin9_smax5");
  int needed_g_f_4[] = {   2,    5}; //   2
  psi_needed_g_f[no_ps][4] = get_vector_from_array_int(needed_g_f_4, 2);
  psi_container_IS_name.push_back("timelikeinvariant_smin10_smax10");
  int needed_g_f_5[] = {   2}; //   1
  psi_needed_g_f[no_ps][5] = get_vector_from_array_int(needed_g_f_5, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin11_smax2");
  int needed_g_f_6[] = {   3,    4}; //   2
  psi_needed_g_f[no_ps][6] = get_vector_from_array_int(needed_g_f_6, 2);
  psi_container_IS_name.push_back("timelikeinvariant_smin12_smax11");
  int needed_g_f_7[] = {   3}; //   1
  psi_needed_g_f[no_ps][7] = get_vector_from_array_int(needed_g_f_7, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin13_smax4");
  int needed_g_f_8[] = {   4}; //   1
  psi_needed_g_f[no_ps][8] = get_vector_from_array_int(needed_g_f_8, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin14_smax3");
  int needed_g_f_9[] = {   5}; //   1
  psi_needed_g_f[no_ps][9] = get_vector_from_array_int(needed_g_f_9, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin15_smax3");
  int needed_g_f_10[] = {   6}; //   1
  psi_needed_g_f[no_ps][10] = get_vector_from_array_int(needed_g_f_10, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin16_smax2");
  int needed_g_f_11[] = {   7}; //   1
  psi_needed_g_f[no_ps][11] = get_vector_from_array_int(needed_g_f_11, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin17_smax4");
  int needed_g_f_12[] = {   8,   21,   31}; //   3
  psi_needed_g_f[no_ps][12] = get_vector_from_array_int(needed_g_f_12, 3);
  psi_container_IS_name.push_back("timelikeinvariant_smin4_smax4");
  int needed_g_f_13[] = {   9,   23,   30}; //   3
  psi_needed_g_f[no_ps][13] = get_vector_from_array_int(needed_g_f_13, 3);
  psi_container_IS_name.push_back("timelikeinvariant_smin2_smax2");
  int needed_g_f_14[] = {  10,   24,   32}; //   3
  psi_needed_g_f[no_ps][14] = get_vector_from_array_int(needed_g_f_14, 3);
  psi_container_IS_name.push_back("timelikeinvariant_smin3_smax3");
  int needed_g_f_15[] = {  11,   22,   29}; //   3
  psi_needed_g_f[no_ps][15] = get_vector_from_array_int(needed_g_f_15, 3);
  psi_container_IS_name.push_back("timelikeinvariant_smin5_smax5");

  psi_container_IS_startvalue[no_ps][2] = psi_container_IS_name.size();
  psi_g_t[no_ps].resize(38);
  psi_needed_g_t[no_ps].resize(38);
  int needed_g_t_0[] = {   0,    1,   22,   25,   29}; //   5
  psi_needed_g_t[no_ps][0] = get_vector_from_array_int(needed_g_t_0, 5);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__1_2__4_56");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__1_2__4_56");
  int needed_g_t_1[] = {   0}; //   1
  psi_needed_g_t[no_ps][1] = get_vector_from_array_int(needed_g_t_1, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_5__32_24");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_5__32_24");
  int needed_g_t_2[] = {   0}; //   1
  psi_needed_g_t[no_ps][2] = get_vector_from_array_int(needed_g_t_2, 1);
  psi_container_IS_name.push_back("t_of_tchannel_0__5_34__8_16");
  psi_container_IS_name.push_back("phi_of_tchannel_0__5_34__8_16");
  int needed_g_t_3[] = {   1}; //   1
  psi_needed_g_t[no_ps][3] = get_vector_from_array_int(needed_g_t_3, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_5__16_40");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_5__16_40");
  int needed_g_t_4[] = {   1}; //   1
  psi_needed_g_t[no_ps][4] = get_vector_from_array_int(needed_g_t_4, 1);
  psi_container_IS_name.push_back("t_of_tchannel_0__5_18__8_32");
  psi_container_IS_name.push_back("phi_of_tchannel_0__5_18__8_32");
  int needed_g_t_5[] = {   2,    6,   23,   26,   30}; //   5
  psi_needed_g_t[no_ps][5] = get_vector_from_array_int(needed_g_t_5, 5);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__16_44");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__16_44");
  int needed_g_t_6[] = {   2}; //   1
  psi_needed_g_t[no_ps][6] = get_vector_from_array_int(needed_g_t_6, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_17__8_36");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_17__8_36");
  int needed_g_t_7[] = {   2}; //   1
  psi_needed_g_t[no_ps][7] = get_vector_from_array_int(needed_g_t_7, 1);
  psi_container_IS_name.push_back("t_of_tchannel_0__17_10__32_4");
  psi_container_IS_name.push_back("phi_of_tchannel_0__17_10__32_4");
  int needed_g_t_8[] = {   3,    7,   21,   27,   31}; //   5
  psi_needed_g_t[no_ps][8] = get_vector_from_array_int(needed_g_t_8, 5);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__32_28");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__32_28");
  int needed_g_t_9[] = {   3}; //   1
  psi_needed_g_t[no_ps][9] = get_vector_from_array_int(needed_g_t_9, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_33__8_20");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_33__8_20");
  int needed_g_t_10[] = {   3}; //   1
  psi_needed_g_t[no_ps][10] = get_vector_from_array_int(needed_g_t_10, 1);
  psi_container_IS_name.push_back("t_of_tchannel_0__33_10__16_4");
  psi_container_IS_name.push_back("phi_of_tchannel_0__33_10__16_4");
  int needed_g_t_11[] = {   4,    5,   24,   28,   32}; //   5
  psi_needed_g_t[no_ps][11] = get_vector_from_array_int(needed_g_t_11, 5);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__1_2__8_52");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__1_2__8_52");
  int needed_g_t_12[] = {   4}; //   1
  psi_needed_g_t[no_ps][12] = get_vector_from_array_int(needed_g_t_12, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_9__32_20");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_9__32_20");
  int needed_g_t_13[] = {   4}; //   1
  psi_needed_g_t[no_ps][13] = get_vector_from_array_int(needed_g_t_13, 1);
  psi_container_IS_name.push_back("t_of_tchannel_0__9_34__4_16");
  psi_container_IS_name.push_back("phi_of_tchannel_0__9_34__4_16");
  int needed_g_t_14[] = {   5}; //   1
  psi_needed_g_t[no_ps][14] = get_vector_from_array_int(needed_g_t_14, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_9__16_36");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_9__16_36");
  int needed_g_t_15[] = {   5}; //   1
  psi_needed_g_t[no_ps][15] = get_vector_from_array_int(needed_g_t_15, 1);
  psi_container_IS_name.push_back("t_of_tchannel_0__9_18__4_32");
  psi_container_IS_name.push_back("phi_of_tchannel_0__9_18__4_32");
  int needed_g_t_16[] = {   6}; //   1
  psi_needed_g_t[no_ps][16] = get_vector_from_array_int(needed_g_t_16, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_17__4_40");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_17__4_40");
  int needed_g_t_17[] = {   6}; //   1
  psi_needed_g_t[no_ps][17] = get_vector_from_array_int(needed_g_t_17, 1);
  psi_container_IS_name.push_back("t_of_tchannel_0__17_6__32_8");
  psi_container_IS_name.push_back("phi_of_tchannel_0__17_6__32_8");
  int needed_g_t_18[] = {   7}; //   1
  psi_needed_g_t[no_ps][18] = get_vector_from_array_int(needed_g_t_18, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_33__4_24");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_33__4_24");
  int needed_g_t_19[] = {   7}; //   1
  psi_needed_g_t[no_ps][19] = get_vector_from_array_int(needed_g_t_19, 1);
  psi_container_IS_name.push_back("t_of_tchannel_0__33_6__16_8");
  psi_container_IS_name.push_back("phi_of_tchannel_0__33_6__16_8");
  int needed_g_t_20[] = {   8,   14}; //   2
  psi_needed_g_t[no_ps][20] = get_vector_from_array_int(needed_g_t_20, 2);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_1__32_28");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_1__32_28");
  int needed_g_t_21[] = {   8}; //   1
  psi_needed_g_t[no_ps][21] = get_vector_from_array_int(needed_g_t_21, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_34__12_16");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_34__12_16");
  int needed_g_t_22[] = {   9,   12}; //   2
  psi_needed_g_t[no_ps][22] = get_vector_from_array_int(needed_g_t_22, 2);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_1__16_44");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_1__16_44");
  int needed_g_t_23[] = {   9}; //   1
  psi_needed_g_t[no_ps][23] = get_vector_from_array_int(needed_g_t_23, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_18__12_32");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_18__12_32");
  int needed_g_t_24[] = {  10,   13}; //   2
  psi_needed_g_t[no_ps][24] = get_vector_from_array_int(needed_g_t_24, 2);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_1__8_52");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_1__8_52");
  int needed_g_t_25[] = {  10}; //   1
  psi_needed_g_t[no_ps][25] = get_vector_from_array_int(needed_g_t_25, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_10__48_4");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_10__48_4");
  int needed_g_t_26[] = {  11,   15}; //   2
  psi_needed_g_t[no_ps][26] = get_vector_from_array_int(needed_g_t_26, 2);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_1__4_56");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_1__4_56");
  int needed_g_t_27[] = {  11}; //   1
  psi_needed_g_t[no_ps][27] = get_vector_from_array_int(needed_g_t_27, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_6__48_8");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_6__48_8");
  int needed_g_t_28[] = {  21}; //   1
  psi_needed_g_t[no_ps][28] = get_vector_from_array_int(needed_g_t_28, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_33__16_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_33__16_12");
  int needed_g_t_29[] = {  22}; //   1
  psi_needed_g_t[no_ps][29] = get_vector_from_array_int(needed_g_t_29, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_5__8_48");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_5__8_48");
  int needed_g_t_30[] = {  23}; //   1
  psi_needed_g_t[no_ps][30] = get_vector_from_array_int(needed_g_t_30, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_17__32_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_17__32_12");
  int needed_g_t_31[] = {  24}; //   1
  psi_needed_g_t[no_ps][31] = get_vector_from_array_int(needed_g_t_31, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__2_9__4_48");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__2_9__4_48");
  int needed_g_t_32[] = {  29}; //   1
  psi_needed_g_t[no_ps][32] = get_vector_from_array_int(needed_g_t_32, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_5__48_8");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_5__48_8");
  int needed_g_t_33[] = {  30}; //   1
  psi_needed_g_t[no_ps][33] = get_vector_from_array_int(needed_g_t_33, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_17__12_32");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_17__12_32");
  int needed_g_t_34[] = {  31}; //   1
  psi_needed_g_t[no_ps][34] = get_vector_from_array_int(needed_g_t_34, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_33__12_16");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_33__12_16");
  int needed_g_t_35[] = {  32}; //   1
  psi_needed_g_t[no_ps][35] = get_vector_from_array_int(needed_g_t_35, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_9__48_4");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_9__48_4");
  int needed_g_t_36[] = {  33}; //   1
  psi_needed_g_t[no_ps][36] = get_vector_from_array_int(needed_g_t_36, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__12_48");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__12_48");
  int needed_g_t_37[] = {  34}; //   1
  psi_needed_g_t[no_ps][37] = get_vector_from_array_int(needed_g_t_37, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__48_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__48_12");

  psi_container_IS_startvalue[no_ps][3] = psi_container_IS_name.size();
  psi_g_d[no_ps].resize(11);
  psi_needed_g_d[no_ps].resize(11);
  int needed_g_d_0[] = {   8,    9,   12,   14,   17,   18,   20,   21,   23,   26,   27,   30,   31,   33,   34}; //  15
  psi_needed_g_d[no_ps][0] = get_vector_from_array_int(needed_g_d_0, 15);
  psi_container_IS_name.push_back("costheta_of_decay_12__4_8");
  psi_container_IS_name.push_back("phi_of_decay_12__4_8");
  int needed_g_d_1[] = {  10,   11,   13,   15,   16,   19,   20,   22,   24,   25,   28,   29,   32,   33,   34}; //  15
  psi_needed_g_d[no_ps][1] = get_vector_from_array_int(needed_g_d_1, 15);
  psi_container_IS_name.push_back("costheta_of_decay_48__16_32");
  psi_container_IS_name.push_back("phi_of_decay_48__16_32");
  int needed_g_d_2[] = {  12,   17,   26}; //   3
  psi_needed_g_d[no_ps][2] = get_vector_from_array_int(needed_g_d_2, 3);
  psi_container_IS_name.push_back("costheta_of_decay_44__32_12");
  psi_container_IS_name.push_back("phi_of_decay_44__32_12");
  int needed_g_d_3[] = {  13,   19,   28}; //   3
  psi_needed_g_d[no_ps][3] = get_vector_from_array_int(needed_g_d_3, 3);
  psi_container_IS_name.push_back("costheta_of_decay_52__4_48");
  psi_container_IS_name.push_back("phi_of_decay_52__4_48");
  int needed_g_d_4[] = {  14,   18,   27}; //   3
  psi_needed_g_d[no_ps][4] = get_vector_from_array_int(needed_g_d_4, 3);
  psi_container_IS_name.push_back("costheta_of_decay_28__16_12");
  psi_container_IS_name.push_back("phi_of_decay_28__16_12");
  int needed_g_d_5[] = {  15,   16,   25}; //   3
  psi_needed_g_d[no_ps][5] = get_vector_from_array_int(needed_g_d_5, 3);
  psi_container_IS_name.push_back("costheta_of_decay_56__8_48");
  psi_container_IS_name.push_back("phi_of_decay_56__8_48");
  int needed_g_d_6[] = {  16}; //   1
  psi_needed_g_d[no_ps][6] = get_vector_from_array_int(needed_g_d_6, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__4_56");
  psi_container_IS_name.push_back("phi_of_decay_60__4_56");
  int needed_g_d_7[] = {  17}; //   1
  psi_needed_g_d[no_ps][7] = get_vector_from_array_int(needed_g_d_7, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__16_44");
  psi_container_IS_name.push_back("phi_of_decay_60__16_44");
  int needed_g_d_8[] = {  18}; //   1
  psi_needed_g_d[no_ps][8] = get_vector_from_array_int(needed_g_d_8, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__32_28");
  psi_container_IS_name.push_back("phi_of_decay_60__32_28");
  int needed_g_d_9[] = {  19}; //   1
  psi_needed_g_d[no_ps][9] = get_vector_from_array_int(needed_g_d_9, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__8_52");
  psi_container_IS_name.push_back("phi_of_decay_60__8_52");
  int needed_g_d_10[] = {  20}; //   1
  psi_needed_g_d[no_ps][10] = get_vector_from_array_int(needed_g_d_10, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__12_48");
  psi_container_IS_name.push_back("phi_of_decay_60__12_48");

  psi_c_p[no_ps].resize(35);
  psi_c_f[no_ps].resize(35);
  psi_c_t[no_ps].resize(35);
  psi_c_d[no_ps].resize(35);

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
  int c_f_1[] = {  2,   3}; //   2
  psi_c_f[no_ps][1] = get_vector_from_array_int(c_f_1, 2);
  int c_t_1[] = {  0,   3,   4}; //   3
  psi_c_t[no_ps][1] = get_vector_from_array_int(c_t_1, 3);
  int c_d_1[] = {}; //   0
  psi_c_d[no_ps][1] = get_vector_from_array_int(c_d_1, 0);

  int c_p_2[] = {}; //   0
  psi_c_p[no_ps][2] = get_vector_from_array_int(c_p_2, 0);
  int c_f_2[] = {  4,   5}; //   2
  psi_c_f[no_ps][2] = get_vector_from_array_int(c_f_2, 2);
  int c_t_2[] = {  5,   6,   7}; //   3
  psi_c_t[no_ps][2] = get_vector_from_array_int(c_t_2, 3);
  int c_d_2[] = {}; //   0
  psi_c_d[no_ps][2] = get_vector_from_array_int(c_d_2, 0);

  int c_p_3[] = {}; //   0
  psi_c_p[no_ps][3] = get_vector_from_array_int(c_p_3, 0);
  int c_f_3[] = {  6,   7}; //   2
  psi_c_f[no_ps][3] = get_vector_from_array_int(c_f_3, 2);
  int c_t_3[] = {  8,   9,  10}; //   3
  psi_c_t[no_ps][3] = get_vector_from_array_int(c_t_3, 3);
  int c_d_3[] = {}; //   0
  psi_c_d[no_ps][3] = get_vector_from_array_int(c_d_3, 0);

  int c_p_4[] = {}; //   0
  psi_c_p[no_ps][4] = get_vector_from_array_int(c_p_4, 0);
  int c_f_4[] = {  6,   8}; //   2
  psi_c_f[no_ps][4] = get_vector_from_array_int(c_f_4, 2);
  int c_t_4[] = { 11,  12,  13}; //   3
  psi_c_t[no_ps][4] = get_vector_from_array_int(c_t_4, 3);
  int c_d_4[] = {}; //   0
  psi_c_d[no_ps][4] = get_vector_from_array_int(c_d_4, 0);

  int c_p_5[] = {}; //   0
  psi_c_p[no_ps][5] = get_vector_from_array_int(c_p_5, 0);
  int c_f_5[] = {  4,   9}; //   2
  psi_c_f[no_ps][5] = get_vector_from_array_int(c_f_5, 2);
  int c_t_5[] = { 11,  14,  15}; //   3
  psi_c_t[no_ps][5] = get_vector_from_array_int(c_t_5, 3);
  int c_d_5[] = {}; //   0
  psi_c_d[no_ps][5] = get_vector_from_array_int(c_d_5, 0);

  int c_p_6[] = {}; //   0
  psi_c_p[no_ps][6] = get_vector_from_array_int(c_p_6, 0);
  int c_f_6[] = {  2,  10}; //   2
  psi_c_f[no_ps][6] = get_vector_from_array_int(c_f_6, 2);
  int c_t_6[] = {  5,  16,  17}; //   3
  psi_c_t[no_ps][6] = get_vector_from_array_int(c_t_6, 3);
  int c_d_6[] = {}; //   0
  psi_c_d[no_ps][6] = get_vector_from_array_int(c_d_6, 0);

  int c_p_7[] = {}; //   0
  psi_c_p[no_ps][7] = get_vector_from_array_int(c_p_7, 0);
  int c_f_7[] = {  0,  11}; //   2
  psi_c_f[no_ps][7] = get_vector_from_array_int(c_f_7, 2);
  int c_t_7[] = {  8,  18,  19}; //   3
  psi_c_t[no_ps][7] = get_vector_from_array_int(c_t_7, 3);
  int c_d_7[] = {}; //   0
  psi_c_d[no_ps][7] = get_vector_from_array_int(c_d_7, 0);

  int c_p_8[] = {  0}; //   1
  psi_c_p[no_ps][8] = get_vector_from_array_int(c_p_8, 1);
  int c_f_8[] = { 12}; //   1
  psi_c_f[no_ps][8] = get_vector_from_array_int(c_f_8, 1);
  int c_t_8[] = { 20,  21}; //   2
  psi_c_t[no_ps][8] = get_vector_from_array_int(c_t_8, 2);
  int c_d_8[] = {  0}; //   1
  psi_c_d[no_ps][8] = get_vector_from_array_int(c_d_8, 1);

  int c_p_9[] = {  0}; //   1
  psi_c_p[no_ps][9] = get_vector_from_array_int(c_p_9, 1);
  int c_f_9[] = { 13}; //   1
  psi_c_f[no_ps][9] = get_vector_from_array_int(c_f_9, 1);
  int c_t_9[] = { 22,  23}; //   2
  psi_c_t[no_ps][9] = get_vector_from_array_int(c_t_9, 2);
  int c_d_9[] = {  0}; //   1
  psi_c_d[no_ps][9] = get_vector_from_array_int(c_d_9, 1);

  int c_p_10[] = {  1}; //   1
  psi_c_p[no_ps][10] = get_vector_from_array_int(c_p_10, 1);
  int c_f_10[] = { 14}; //   1
  psi_c_f[no_ps][10] = get_vector_from_array_int(c_f_10, 1);
  int c_t_10[] = { 24,  25}; //   2
  psi_c_t[no_ps][10] = get_vector_from_array_int(c_t_10, 2);
  int c_d_10[] = {  1}; //   1
  psi_c_d[no_ps][10] = get_vector_from_array_int(c_d_10, 1);

  int c_p_11[] = {  1}; //   1
  psi_c_p[no_ps][11] = get_vector_from_array_int(c_p_11, 1);
  int c_f_11[] = { 15}; //   1
  psi_c_f[no_ps][11] = get_vector_from_array_int(c_f_11, 1);
  int c_t_11[] = { 26,  27}; //   2
  psi_c_t[no_ps][11] = get_vector_from_array_int(c_t_11, 2);
  int c_d_11[] = {  1}; //   1
  psi_c_d[no_ps][11] = get_vector_from_array_int(c_d_11, 1);

  int c_p_12[] = {  0,   2}; //   2
  psi_c_p[no_ps][12] = get_vector_from_array_int(c_p_12, 2);
  int c_f_12[] = {}; //   0
  psi_c_f[no_ps][12] = get_vector_from_array_int(c_f_12, 0);
  int c_t_12[] = { 22}; //   1
  psi_c_t[no_ps][12] = get_vector_from_array_int(c_t_12, 1);
  int c_d_12[] = {  2,   0}; //   2
  psi_c_d[no_ps][12] = get_vector_from_array_int(c_d_12, 2);

  int c_p_13[] = {  1,   3}; //   2
  psi_c_p[no_ps][13] = get_vector_from_array_int(c_p_13, 2);
  int c_f_13[] = {}; //   0
  psi_c_f[no_ps][13] = get_vector_from_array_int(c_f_13, 0);
  int c_t_13[] = { 24}; //   1
  psi_c_t[no_ps][13] = get_vector_from_array_int(c_t_13, 1);
  int c_d_13[] = {  3,   1}; //   2
  psi_c_d[no_ps][13] = get_vector_from_array_int(c_d_13, 2);

  int c_p_14[] = {  0,   4}; //   2
  psi_c_p[no_ps][14] = get_vector_from_array_int(c_p_14, 2);
  int c_f_14[] = {}; //   0
  psi_c_f[no_ps][14] = get_vector_from_array_int(c_f_14, 0);
  int c_t_14[] = { 20}; //   1
  psi_c_t[no_ps][14] = get_vector_from_array_int(c_t_14, 1);
  int c_d_14[] = {  4,   0}; //   2
  psi_c_d[no_ps][14] = get_vector_from_array_int(c_d_14, 2);

  int c_p_15[] = {  1,   5}; //   2
  psi_c_p[no_ps][15] = get_vector_from_array_int(c_p_15, 2);
  int c_f_15[] = {}; //   0
  psi_c_f[no_ps][15] = get_vector_from_array_int(c_f_15, 0);
  int c_t_15[] = { 26}; //   1
  psi_c_t[no_ps][15] = get_vector_from_array_int(c_t_15, 1);
  int c_d_15[] = {  5,   1}; //   2
  psi_c_d[no_ps][15] = get_vector_from_array_int(c_d_15, 2);

  int c_p_16[] = {  1,   5}; //   2
  psi_c_p[no_ps][16] = get_vector_from_array_int(c_p_16, 2);
  int c_f_16[] = {}; //   0
  psi_c_f[no_ps][16] = get_vector_from_array_int(c_f_16, 0);
  int c_t_16[] = {}; //   0
  psi_c_t[no_ps][16] = get_vector_from_array_int(c_t_16, 0);
  int c_d_16[] = {  6,   5,   1}; //   3
  psi_c_d[no_ps][16] = get_vector_from_array_int(c_d_16, 3);

  int c_p_17[] = {  0,   2}; //   2
  psi_c_p[no_ps][17] = get_vector_from_array_int(c_p_17, 2);
  int c_f_17[] = {}; //   0
  psi_c_f[no_ps][17] = get_vector_from_array_int(c_f_17, 0);
  int c_t_17[] = {}; //   0
  psi_c_t[no_ps][17] = get_vector_from_array_int(c_t_17, 0);
  int c_d_17[] = {  7,   2,   0}; //   3
  psi_c_d[no_ps][17] = get_vector_from_array_int(c_d_17, 3);

  int c_p_18[] = {  0,   4}; //   2
  psi_c_p[no_ps][18] = get_vector_from_array_int(c_p_18, 2);
  int c_f_18[] = {}; //   0
  psi_c_f[no_ps][18] = get_vector_from_array_int(c_f_18, 0);
  int c_t_18[] = {}; //   0
  psi_c_t[no_ps][18] = get_vector_from_array_int(c_t_18, 0);
  int c_d_18[] = {  8,   4,   0}; //   3
  psi_c_d[no_ps][18] = get_vector_from_array_int(c_d_18, 3);

  int c_p_19[] = {  1,   3}; //   2
  psi_c_p[no_ps][19] = get_vector_from_array_int(c_p_19, 2);
  int c_f_19[] = {}; //   0
  psi_c_f[no_ps][19] = get_vector_from_array_int(c_f_19, 0);
  int c_t_19[] = {}; //   0
  psi_c_t[no_ps][19] = get_vector_from_array_int(c_t_19, 0);
  int c_d_19[] = {  9,   3,   1}; //   3
  psi_c_d[no_ps][19] = get_vector_from_array_int(c_d_19, 3);

  int c_p_20[] = {  1,   6}; //   2
  psi_c_p[no_ps][20] = get_vector_from_array_int(c_p_20, 2);
  int c_f_20[] = {}; //   0
  psi_c_f[no_ps][20] = get_vector_from_array_int(c_f_20, 0);
  int c_t_20[] = {}; //   0
  psi_c_t[no_ps][20] = get_vector_from_array_int(c_t_20, 0);
  int c_d_20[] = { 10,   0,   1}; //   3
  psi_c_d[no_ps][20] = get_vector_from_array_int(c_d_20, 3);

  int c_p_21[] = {  0}; //   1
  psi_c_p[no_ps][21] = get_vector_from_array_int(c_p_21, 1);
  int c_f_21[] = { 12}; //   1
  psi_c_f[no_ps][21] = get_vector_from_array_int(c_f_21, 1);
  int c_t_21[] = {  8,  28}; //   2
  psi_c_t[no_ps][21] = get_vector_from_array_int(c_t_21, 2);
  int c_d_21[] = {  0}; //   1
  psi_c_d[no_ps][21] = get_vector_from_array_int(c_d_21, 1);

  int c_p_22[] = {  1}; //   1
  psi_c_p[no_ps][22] = get_vector_from_array_int(c_p_22, 1);
  int c_f_22[] = { 15}; //   1
  psi_c_f[no_ps][22] = get_vector_from_array_int(c_f_22, 1);
  int c_t_22[] = {  0,  29}; //   2
  psi_c_t[no_ps][22] = get_vector_from_array_int(c_t_22, 2);
  int c_d_22[] = {  1}; //   1
  psi_c_d[no_ps][22] = get_vector_from_array_int(c_d_22, 1);

  int c_p_23[] = {  0}; //   1
  psi_c_p[no_ps][23] = get_vector_from_array_int(c_p_23, 1);
  int c_f_23[] = { 13}; //   1
  psi_c_f[no_ps][23] = get_vector_from_array_int(c_f_23, 1);
  int c_t_23[] = {  5,  30}; //   2
  psi_c_t[no_ps][23] = get_vector_from_array_int(c_t_23, 2);
  int c_d_23[] = {  0}; //   1
  psi_c_d[no_ps][23] = get_vector_from_array_int(c_d_23, 1);

  int c_p_24[] = {  1}; //   1
  psi_c_p[no_ps][24] = get_vector_from_array_int(c_p_24, 1);
  int c_f_24[] = { 14}; //   1
  psi_c_f[no_ps][24] = get_vector_from_array_int(c_f_24, 1);
  int c_t_24[] = { 11,  31}; //   2
  psi_c_t[no_ps][24] = get_vector_from_array_int(c_t_24, 2);
  int c_d_24[] = {  1}; //   1
  psi_c_d[no_ps][24] = get_vector_from_array_int(c_d_24, 1);

  int c_p_25[] = {  1,   5}; //   2
  psi_c_p[no_ps][25] = get_vector_from_array_int(c_p_25, 2);
  int c_f_25[] = {}; //   0
  psi_c_f[no_ps][25] = get_vector_from_array_int(c_f_25, 0);
  int c_t_25[] = {  0}; //   1
  psi_c_t[no_ps][25] = get_vector_from_array_int(c_t_25, 1);
  int c_d_25[] = {  5,   1}; //   2
  psi_c_d[no_ps][25] = get_vector_from_array_int(c_d_25, 2);

  int c_p_26[] = {  0,   2}; //   2
  psi_c_p[no_ps][26] = get_vector_from_array_int(c_p_26, 2);
  int c_f_26[] = {}; //   0
  psi_c_f[no_ps][26] = get_vector_from_array_int(c_f_26, 0);
  int c_t_26[] = {  5}; //   1
  psi_c_t[no_ps][26] = get_vector_from_array_int(c_t_26, 1);
  int c_d_26[] = {  2,   0}; //   2
  psi_c_d[no_ps][26] = get_vector_from_array_int(c_d_26, 2);

  int c_p_27[] = {  0,   4}; //   2
  psi_c_p[no_ps][27] = get_vector_from_array_int(c_p_27, 2);
  int c_f_27[] = {}; //   0
  psi_c_f[no_ps][27] = get_vector_from_array_int(c_f_27, 0);
  int c_t_27[] = {  8}; //   1
  psi_c_t[no_ps][27] = get_vector_from_array_int(c_t_27, 1);
  int c_d_27[] = {  4,   0}; //   2
  psi_c_d[no_ps][27] = get_vector_from_array_int(c_d_27, 2);

  int c_p_28[] = {  1,   3}; //   2
  psi_c_p[no_ps][28] = get_vector_from_array_int(c_p_28, 2);
  int c_f_28[] = {}; //   0
  psi_c_f[no_ps][28] = get_vector_from_array_int(c_f_28, 0);
  int c_t_28[] = { 11}; //   1
  psi_c_t[no_ps][28] = get_vector_from_array_int(c_t_28, 1);
  int c_d_28[] = {  3,   1}; //   2
  psi_c_d[no_ps][28] = get_vector_from_array_int(c_d_28, 2);

  int c_p_29[] = {  1}; //   1
  psi_c_p[no_ps][29] = get_vector_from_array_int(c_p_29, 1);
  int c_f_29[] = { 15}; //   1
  psi_c_f[no_ps][29] = get_vector_from_array_int(c_f_29, 1);
  int c_t_29[] = {  0,  32}; //   2
  psi_c_t[no_ps][29] = get_vector_from_array_int(c_t_29, 2);
  int c_d_29[] = {  1}; //   1
  psi_c_d[no_ps][29] = get_vector_from_array_int(c_d_29, 1);

  int c_p_30[] = {  0}; //   1
  psi_c_p[no_ps][30] = get_vector_from_array_int(c_p_30, 1);
  int c_f_30[] = { 13}; //   1
  psi_c_f[no_ps][30] = get_vector_from_array_int(c_f_30, 1);
  int c_t_30[] = {  5,  33}; //   2
  psi_c_t[no_ps][30] = get_vector_from_array_int(c_t_30, 2);
  int c_d_30[] = {  0}; //   1
  psi_c_d[no_ps][30] = get_vector_from_array_int(c_d_30, 1);

  int c_p_31[] = {  0}; //   1
  psi_c_p[no_ps][31] = get_vector_from_array_int(c_p_31, 1);
  int c_f_31[] = { 12}; //   1
  psi_c_f[no_ps][31] = get_vector_from_array_int(c_f_31, 1);
  int c_t_31[] = {  8,  34}; //   2
  psi_c_t[no_ps][31] = get_vector_from_array_int(c_t_31, 2);
  int c_d_31[] = {  0}; //   1
  psi_c_d[no_ps][31] = get_vector_from_array_int(c_d_31, 1);

  int c_p_32[] = {  1}; //   1
  psi_c_p[no_ps][32] = get_vector_from_array_int(c_p_32, 1);
  int c_f_32[] = { 14}; //   1
  psi_c_f[no_ps][32] = get_vector_from_array_int(c_f_32, 1);
  int c_t_32[] = { 11,  35}; //   2
  psi_c_t[no_ps][32] = get_vector_from_array_int(c_t_32, 2);
  int c_d_32[] = {  1}; //   1
  psi_c_d[no_ps][32] = get_vector_from_array_int(c_d_32, 1);

  int c_p_33[] = {  1,   6}; //   2
  psi_c_p[no_ps][33] = get_vector_from_array_int(c_p_33, 2);
  int c_f_33[] = {}; //   0
  psi_c_f[no_ps][33] = get_vector_from_array_int(c_f_33, 0);
  int c_t_33[] = { 36}; //   1
  psi_c_t[no_ps][33] = get_vector_from_array_int(c_t_33, 1);
  int c_d_33[] = {  0,   1}; //   2
  psi_c_d[no_ps][33] = get_vector_from_array_int(c_d_33, 2);

  int c_p_34[] = {  0,   7}; //   2
  psi_c_p[no_ps][34] = get_vector_from_array_int(c_p_34, 2);
  int c_f_34[] = {}; //   0
  psi_c_f[no_ps][34] = get_vector_from_array_int(c_f_34, 0);
  int c_t_34[] = { 37}; //   1
  psi_c_t[no_ps][34] = get_vector_from_array_int(c_t_34, 1);
  int c_d_34[] = {  1,   0}; //   2
  psi_c_d[no_ps][34] = get_vector_from_array_int(c_d_34, 2);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
