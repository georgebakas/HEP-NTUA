#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ax_psp_400_bbx_ttxgg(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_400_bbx_ttxgg");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  psi_v_smin[no_ps].resize(18);
  psi_needed_v_smin[no_ps].resize(18);
  int needed_v_smin_0[] = {   0,    1,    4,    5,   13,   14,   20,   23,   24,   27,   28,   31,   32,   33,   34}; //  15
  psi_needed_v_smin[no_ps][  0] = get_vector_from_array_int(needed_v_smin_0, 15);
  int needed_v_smin_1[] = {   2,    9,   10,   21,   25}; //   5
  psi_needed_v_smin[no_ps][  1] = get_vector_from_array_int(needed_v_smin_1, 5);
  int needed_v_smin_2[] = {   2,   10,   25}; //   3
  psi_needed_v_smin[no_ps][  2] = get_vector_from_array_int(needed_v_smin_2, 3);
  int needed_v_smin_3[] = {   3,   11,   12,   22,   26}; //   5
  psi_needed_v_smin[no_ps][  3] = get_vector_from_array_int(needed_v_smin_3, 5);
  int needed_v_smin_4[] = {   3,   12,   26}; //   3
  psi_needed_v_smin[no_ps][  4] = get_vector_from_array_int(needed_v_smin_4, 3);
  int needed_v_smin_5[] = {   4,   13,   27,    1,   24,   31}; //   6
  psi_needed_v_smin[no_ps][  5] = get_vector_from_array_int(needed_v_smin_5, 6);
  int needed_v_smin_6[] = {   5,   14,   28,    0,   23,   32}; //   6
  psi_needed_v_smin[no_ps][  6] = get_vector_from_array_int(needed_v_smin_6, 6);
  int needed_v_smin_7[] = {   6,   16,   17,   22,   29}; //   5
  psi_needed_v_smin[no_ps][  7] = get_vector_from_array_int(needed_v_smin_7, 5);
  int needed_v_smin_8[] = {   6,   17,   29}; //   3
  psi_needed_v_smin[no_ps][  8] = get_vector_from_array_int(needed_v_smin_8, 3);
  int needed_v_smin_9[] = {   7,   18,   19,   21,   30}; //   5
  psi_needed_v_smin[no_ps][  9] = get_vector_from_array_int(needed_v_smin_9, 5);
  int needed_v_smin_10[] = {   7,   19,   30}; //   3
  psi_needed_v_smin[no_ps][ 10] = get_vector_from_array_int(needed_v_smin_10, 3);
  int needed_v_smin_11[] = {   8,   15,   20,   33,   34}; //   5
  psi_needed_v_smin[no_ps][ 11] = get_vector_from_array_int(needed_v_smin_11, 5);
  int needed_v_smin_12[] = {   8}; //   1
  psi_needed_v_smin[no_ps][ 12] = get_vector_from_array_int(needed_v_smin_12, 1);
  int needed_v_smin_13[] = {   9}; //   1
  psi_needed_v_smin[no_ps][ 13] = get_vector_from_array_int(needed_v_smin_13, 1);
  int needed_v_smin_14[] = {  11}; //   1
  psi_needed_v_smin[no_ps][ 14] = get_vector_from_array_int(needed_v_smin_14, 1);
  int needed_v_smin_15[] = {  15}; //   1
  psi_needed_v_smin[no_ps][ 15] = get_vector_from_array_int(needed_v_smin_15, 1);
  int needed_v_smin_16[] = {  16}; //   1
  psi_needed_v_smin[no_ps][ 16] = get_vector_from_array_int(needed_v_smin_16, 1);
  int needed_v_smin_17[] = {  18}; //   1
  psi_needed_v_smin[no_ps][ 17] = get_vector_from_array_int(needed_v_smin_17, 1);

  psi_v_smax[no_ps].resize(14);
  psi_needed_v_smax[no_ps].resize(14);
  int needed_v_smax_0[] = {   0,    1,    4,    5,   13,   14,   23,   24,   27,   28,   31,   32,   34}; //  13
  psi_needed_v_smax[no_ps][  0] = get_vector_from_array_int(needed_v_smax_0, 13);
  int needed_v_smax_1[] = {   2,    9,   10,   21,   25}; //   5
  psi_needed_v_smax[no_ps][  1] = get_vector_from_array_int(needed_v_smax_1, 5);
  int needed_v_smax_2[] = {   2,    5,    6,   10,   14,   17,   25,   28,   29,    0,   23,   32}; //  12
  psi_needed_v_smax[no_ps][  2] = get_vector_from_array_int(needed_v_smax_2, 12);
  int needed_v_smax_3[] = {   3,   11,   12,   22,   26}; //   5
  psi_needed_v_smax[no_ps][  3] = get_vector_from_array_int(needed_v_smax_3, 5);
  int needed_v_smax_4[] = {   3,    4,    7,   12,   13,   19,   26,   27,   30,    1,   24,   31}; //  12
  psi_needed_v_smax[no_ps][  4] = get_vector_from_array_int(needed_v_smax_4, 12);
  int needed_v_smax_5[] = {   6,   16,   17,   29}; //   4
  psi_needed_v_smax[no_ps][  5] = get_vector_from_array_int(needed_v_smax_5, 4);
  int needed_v_smax_6[] = {   7,   18,   19,   30}; //   4
  psi_needed_v_smax[no_ps][  6] = get_vector_from_array_int(needed_v_smax_6, 4);
  int needed_v_smax_7[] = {   8,   15,   20,   33}; //   4
  psi_needed_v_smax[no_ps][  7] = get_vector_from_array_int(needed_v_smax_7, 4);
  int needed_v_smax_8[] = {   8,    9,   11}; //   3
  psi_needed_v_smax[no_ps][  8] = get_vector_from_array_int(needed_v_smax_8, 3);
  int needed_v_smax_9[] = {  15,   16,   18}; //   3
  psi_needed_v_smax[no_ps][  9] = get_vector_from_array_int(needed_v_smax_9, 3);
  int needed_v_smax_10[] = {  20,   33}; //   2
  psi_needed_v_smax[no_ps][ 10] = get_vector_from_array_int(needed_v_smax_10, 2);
  int needed_v_smax_11[] = {  21}; //   1
  psi_needed_v_smax[no_ps][ 11] = get_vector_from_array_int(needed_v_smax_11, 1);
  int needed_v_smax_12[] = {  22}; //   1
  psi_needed_v_smax[no_ps][ 12] = get_vector_from_array_int(needed_v_smax_12, 1);
  int needed_v_smax_13[] = {  34}; //   1
  psi_needed_v_smax[no_ps][ 13] = get_vector_from_array_int(needed_v_smax_13, 1);

  psi_container_IS_startvalue[no_ps][0] = psi_container_IS_name.size();
  psi_g_p[no_ps].resize(22);
  psi_needed_g_p[no_ps].resize(22);
  int needed_g_p_0[] = {   0,    1,    4,    5,   13,   14,   23,   24,   27,   28,   31,   32,   34}; //  13
  psi_needed_g_p[no_ps][0] = get_vector_from_array_int(needed_g_p_0, 13);
  psi_container_IS_name.push_back("propagator_0_smin0_smax0");
  int needed_g_p_1[] = {   2,    9,   10,   21,   25}; //   5
  psi_needed_g_p[no_ps][1] = get_vector_from_array_int(needed_g_p_1, 5);
  psi_container_IS_name.push_back("propagator_6_smin1_smax1");
  int needed_g_p_2[] = {   2,   10,   25}; //   3
  psi_needed_g_p[no_ps][2] = get_vector_from_array_int(needed_g_p_2, 3);
  psi_container_IS_name.push_back("propagator_0_smin2_smax2");
  int needed_g_p_3[] = {   3,   11,   12,   22,   26}; //   5
  psi_needed_g_p[no_ps][3] = get_vector_from_array_int(needed_g_p_3, 5);
  psi_container_IS_name.push_back("propagator_6_smin3_smax3");
  int needed_g_p_4[] = {   3,   12,   26}; //   3
  psi_needed_g_p[no_ps][4] = get_vector_from_array_int(needed_g_p_4, 3);
  psi_container_IS_name.push_back("propagator_0_smin4_smax4");
  int needed_g_p_5[] = {   4,   13,   27}; //   3
  psi_needed_g_p[no_ps][5] = get_vector_from_array_int(needed_g_p_5, 3);
  psi_container_IS_name.push_back("propagator_0_smin5_smax4");
  int needed_g_p_6[] = {   5,   14,   28}; //   3
  psi_needed_g_p[no_ps][6] = get_vector_from_array_int(needed_g_p_6, 3);
  psi_container_IS_name.push_back("propagator_0_smin6_smax2");
  int needed_g_p_7[] = {   6,   16,   17,   29}; //   4
  psi_needed_g_p[no_ps][7] = get_vector_from_array_int(needed_g_p_7, 4);
  psi_container_IS_name.push_back("propagator_6_smin7_smax5");
  int needed_g_p_8[] = {   6,   17,   29}; //   3
  psi_needed_g_p[no_ps][8] = get_vector_from_array_int(needed_g_p_8, 3);
  psi_container_IS_name.push_back("propagator_0_smin8_smax2");
  int needed_g_p_9[] = {   7,   18,   19,   30}; //   4
  psi_needed_g_p[no_ps][9] = get_vector_from_array_int(needed_g_p_9, 4);
  psi_container_IS_name.push_back("propagator_6_smin9_smax6");
  int needed_g_p_10[] = {   7,   19,   30}; //   3
  psi_needed_g_p[no_ps][10] = get_vector_from_array_int(needed_g_p_10, 3);
  psi_container_IS_name.push_back("propagator_0_smin10_smax4");
  int needed_g_p_11[] = {   8,   15,   20,   33}; //   4
  psi_needed_g_p[no_ps][11] = get_vector_from_array_int(needed_g_p_11, 4);
  psi_container_IS_name.push_back("propagator_0_smin11_smax7");
  int needed_g_p_12[] = {   8}; //   1
  psi_needed_g_p[no_ps][12] = get_vector_from_array_int(needed_g_p_12, 1);
  psi_container_IS_name.push_back("propagator_6_smin12_smax8");
  int needed_g_p_13[] = {   9}; //   1
  psi_needed_g_p[no_ps][13] = get_vector_from_array_int(needed_g_p_13, 1);
  psi_container_IS_name.push_back("propagator_6_smin13_smax8");
  int needed_g_p_14[] = {  11}; //   1
  psi_needed_g_p[no_ps][14] = get_vector_from_array_int(needed_g_p_14, 1);
  psi_container_IS_name.push_back("propagator_6_smin14_smax8");
  int needed_g_p_15[] = {  15}; //   1
  psi_needed_g_p[no_ps][15] = get_vector_from_array_int(needed_g_p_15, 1);
  psi_container_IS_name.push_back("propagator_6_smin15_smax9");
  int needed_g_p_16[] = {  16}; //   1
  psi_needed_g_p[no_ps][16] = get_vector_from_array_int(needed_g_p_16, 1);
  psi_container_IS_name.push_back("propagator_6_smin16_smax9");
  int needed_g_p_17[] = {  18}; //   1
  psi_needed_g_p[no_ps][17] = get_vector_from_array_int(needed_g_p_17, 1);
  psi_container_IS_name.push_back("propagator_6_smin17_smax9");
  int needed_g_p_18[] = {  20,   33}; //   2
  psi_needed_g_p[no_ps][18] = get_vector_from_array_int(needed_g_p_18, 2);
  psi_container_IS_name.push_back("propagator_0_smin0_smax10");
  int needed_g_p_19[] = {  21}; //   1
  psi_needed_g_p[no_ps][19] = get_vector_from_array_int(needed_g_p_19, 1);
  psi_container_IS_name.push_back("propagator_6_smin9_smax11");
  int needed_g_p_20[] = {  22}; //   1
  psi_needed_g_p[no_ps][20] = get_vector_from_array_int(needed_g_p_20, 1);
  psi_container_IS_name.push_back("propagator_6_smin7_smax12");
  int needed_g_p_21[] = {  34}; //   1
  psi_needed_g_p[no_ps][21] = get_vector_from_array_int(needed_g_p_21, 1);
  psi_container_IS_name.push_back("propagator_0_smin11_smax13");

  psi_container_IS_startvalue[no_ps][1] = psi_container_IS_name.size();
  psi_g_f[no_ps].resize(2);
  psi_needed_g_f[no_ps].resize(2);
  int needed_g_f_0[] = {   0,   23,   32}; //   3
  psi_needed_g_f[no_ps][0] = get_vector_from_array_int(needed_g_f_0, 3);
  psi_container_IS_name.push_back("timelikeinvariant_smin6_smax2");
  int needed_g_f_1[] = {   1,   24,   31}; //   3
  psi_needed_g_f[no_ps][1] = get_vector_from_array_int(needed_g_f_1, 3);
  psi_container_IS_name.push_back("timelikeinvariant_smin5_smax4");

  psi_container_IS_startvalue[no_ps][2] = psi_container_IS_name.size();
  psi_g_t[no_ps].resize(12);
  psi_needed_g_t[no_ps].resize(12);
  int needed_g_t_0[] = {   0,    2,    5,    6}; //   4
  psi_needed_g_t[no_ps][0] = get_vector_from_array_int(needed_g_t_0, 4);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__2_1__32_28");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__2_1__32_28");
  int needed_g_t_1[] = {   0}; //   1
  psi_needed_g_t[no_ps][1] = get_vector_from_array_int(needed_g_t_1, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__1_34__12_16");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__1_34__12_16");
  int needed_g_t_2[] = {   1,    3,    4,    7}; //   4
  psi_needed_g_t[no_ps][2] = get_vector_from_array_int(needed_g_t_2, 4);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__2_1__16_44");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__2_1__16_44");
  int needed_g_t_3[] = {   1}; //   1
  psi_needed_g_t[no_ps][3] = get_vector_from_array_int(needed_g_t_3, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__1_18__12_32");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__1_18__12_32");
  int needed_g_t_4[] = {  23,   25,   28,   29,   32}; //   5
  psi_needed_g_t[no_ps][4] = get_vector_from_array_int(needed_g_t_4, 5);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__1_2__32_28");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__1_2__32_28");
  int needed_g_t_5[] = {  23}; //   1
  psi_needed_g_t[no_ps][5] = get_vector_from_array_int(needed_g_t_5, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__2_33__16_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__2_33__16_12");
  int needed_g_t_6[] = {  24,   26,   27,   30,   31}; //   5
  psi_needed_g_t[no_ps][6] = get_vector_from_array_int(needed_g_t_6, 5);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__1_2__16_44");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__1_2__16_44");
  int needed_g_t_7[] = {  24}; //   1
  psi_needed_g_t[no_ps][7] = get_vector_from_array_int(needed_g_t_7, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__2_17__32_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__2_17__32_12");
  int needed_g_t_8[] = {  31}; //   1
  psi_needed_g_t[no_ps][8] = get_vector_from_array_int(needed_g_t_8, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__2_17__12_32");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__2_17__12_32");
  int needed_g_t_9[] = {  32}; //   1
  psi_needed_g_t[no_ps][9] = get_vector_from_array_int(needed_g_t_9, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__2_33__12_16");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__2_33__12_16");
  int needed_g_t_10[] = {  33}; //   1
  psi_needed_g_t[no_ps][10] = get_vector_from_array_int(needed_g_t_10, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__1_2__12_48");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__1_2__12_48");
  int needed_g_t_11[] = {  34}; //   1
  psi_needed_g_t[no_ps][11] = get_vector_from_array_int(needed_g_t_11, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__1_2__48_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__1_2__48_12");

  psi_container_IS_startvalue[no_ps][3] = psi_container_IS_name.size();
  psi_g_d[no_ps].resize(25);
  psi_needed_g_d[no_ps].resize(25);
  int needed_g_d_0[] = {   0,    1,    4,    5,   13,   14,   20,   23,   24,   27,   28,   31,   32,   33,   34}; //  15
  psi_needed_g_d[no_ps][0] = get_vector_from_array_int(needed_g_d_0, 15);
  psi_container_IS_name.push_back("costheta_of_decay_12__4_8");
  psi_container_IS_name.push_back("phi_of_decay_12__4_8");
  int needed_g_d_1[] = {   2,   10,   25}; //   3
  psi_needed_g_d[no_ps][1] = get_vector_from_array_int(needed_g_d_1, 3);
  psi_container_IS_name.push_back("costheta_of_decay_28__4_24");
  psi_container_IS_name.push_back("phi_of_decay_28__4_24");
  int needed_g_d_2[] = {   2,    9,   10,   21,   25}; //   5
  psi_needed_g_d[no_ps][2] = get_vector_from_array_int(needed_g_d_2, 5);
  psi_container_IS_name.push_back("costheta_of_decay_24__8_16");
  psi_container_IS_name.push_back("phi_of_decay_24__8_16");
  int needed_g_d_3[] = {   3,   12,   26}; //   3
  psi_needed_g_d[no_ps][3] = get_vector_from_array_int(needed_g_d_3, 3);
  psi_container_IS_name.push_back("costheta_of_decay_44__4_40");
  psi_container_IS_name.push_back("phi_of_decay_44__4_40");
  int needed_g_d_4[] = {   3,   11,   12,   22,   26}; //   5
  psi_needed_g_d[no_ps][4] = get_vector_from_array_int(needed_g_d_4, 5);
  psi_container_IS_name.push_back("costheta_of_decay_40__8_32");
  psi_container_IS_name.push_back("phi_of_decay_40__8_32");
  int needed_g_d_5[] = {   4,   13,   27}; //   3
  psi_needed_g_d[no_ps][5] = get_vector_from_array_int(needed_g_d_5, 3);
  psi_container_IS_name.push_back("costheta_of_decay_44__32_12");
  psi_container_IS_name.push_back("phi_of_decay_44__32_12");
  int needed_g_d_6[] = {   5,   14,   28}; //   3
  psi_needed_g_d[no_ps][6] = get_vector_from_array_int(needed_g_d_6, 3);
  psi_container_IS_name.push_back("costheta_of_decay_28__16_12");
  psi_container_IS_name.push_back("phi_of_decay_28__16_12");
  int needed_g_d_7[] = {   6,   17,   29}; //   3
  psi_needed_g_d[no_ps][7] = get_vector_from_array_int(needed_g_d_7, 3);
  psi_container_IS_name.push_back("costheta_of_decay_28__8_20");
  psi_container_IS_name.push_back("phi_of_decay_28__8_20");
  int needed_g_d_8[] = {   6,   16,   17,   22,   29}; //   5
  psi_needed_g_d[no_ps][8] = get_vector_from_array_int(needed_g_d_8, 5);
  psi_container_IS_name.push_back("costheta_of_decay_20__4_16");
  psi_container_IS_name.push_back("phi_of_decay_20__4_16");
  int needed_g_d_9[] = {   7,   19,   30}; //   3
  psi_needed_g_d[no_ps][9] = get_vector_from_array_int(needed_g_d_9, 3);
  psi_container_IS_name.push_back("costheta_of_decay_44__8_36");
  psi_container_IS_name.push_back("phi_of_decay_44__8_36");
  int needed_g_d_10[] = {   7,   18,   19,   21,   30}; //   5
  psi_needed_g_d[no_ps][10] = get_vector_from_array_int(needed_g_d_10, 5);
  psi_container_IS_name.push_back("costheta_of_decay_36__4_32");
  psi_container_IS_name.push_back("phi_of_decay_36__4_32");
  int needed_g_d_11[] = {   8,    9,   11}; //   3
  psi_needed_g_d[no_ps][11] = get_vector_from_array_int(needed_g_d_11, 3);
  psi_container_IS_name.push_back("costheta_of_decay_60__4_56");
  psi_container_IS_name.push_back("phi_of_decay_60__4_56");
  int needed_g_d_12[] = {   8}; //   1
  psi_needed_g_d[no_ps][12] = get_vector_from_array_int(needed_g_d_12, 1);
  psi_container_IS_name.push_back("costheta_of_decay_56__8_48");
  psi_container_IS_name.push_back("phi_of_decay_56__8_48");
  int needed_g_d_13[] = {   8,   15,   20,   33,   34}; //   5
  psi_needed_g_d[no_ps][13] = get_vector_from_array_int(needed_g_d_13, 5);
  psi_container_IS_name.push_back("costheta_of_decay_48__16_32");
  psi_container_IS_name.push_back("phi_of_decay_48__16_32");
  int needed_g_d_14[] = {   9}; //   1
  psi_needed_g_d[no_ps][14] = get_vector_from_array_int(needed_g_d_14, 1);
  psi_container_IS_name.push_back("costheta_of_decay_56__32_24");
  psi_container_IS_name.push_back("phi_of_decay_56__32_24");
  int needed_g_d_15[] = {  10,   14,   17}; //   3
  psi_needed_g_d[no_ps][15] = get_vector_from_array_int(needed_g_d_15, 3);
  psi_container_IS_name.push_back("costheta_of_decay_60__32_28");
  psi_container_IS_name.push_back("phi_of_decay_60__32_28");
  int needed_g_d_16[] = {  11}; //   1
  psi_needed_g_d[no_ps][16] = get_vector_from_array_int(needed_g_d_16, 1);
  psi_container_IS_name.push_back("costheta_of_decay_56__16_40");
  psi_container_IS_name.push_back("phi_of_decay_56__16_40");
  int needed_g_d_17[] = {  12,   13,   19}; //   3
  psi_needed_g_d[no_ps][17] = get_vector_from_array_int(needed_g_d_17, 3);
  psi_container_IS_name.push_back("costheta_of_decay_60__16_44");
  psi_container_IS_name.push_back("phi_of_decay_60__16_44");
  int needed_g_d_18[] = {  15,   16,   18}; //   3
  psi_needed_g_d[no_ps][18] = get_vector_from_array_int(needed_g_d_18, 3);
  psi_container_IS_name.push_back("costheta_of_decay_60__8_52");
  psi_container_IS_name.push_back("phi_of_decay_60__8_52");
  int needed_g_d_19[] = {  15}; //   1
  psi_needed_g_d[no_ps][19] = get_vector_from_array_int(needed_g_d_19, 1);
  psi_container_IS_name.push_back("costheta_of_decay_52__4_48");
  psi_container_IS_name.push_back("phi_of_decay_52__4_48");
  int needed_g_d_20[] = {  16}; //   1
  psi_needed_g_d[no_ps][20] = get_vector_from_array_int(needed_g_d_20, 1);
  psi_container_IS_name.push_back("costheta_of_decay_52__32_20");
  psi_container_IS_name.push_back("phi_of_decay_52__32_20");
  int needed_g_d_21[] = {  18}; //   1
  psi_needed_g_d[no_ps][21] = get_vector_from_array_int(needed_g_d_21, 1);
  psi_container_IS_name.push_back("costheta_of_decay_52__16_36");
  psi_container_IS_name.push_back("phi_of_decay_52__16_36");
  int needed_g_d_22[] = {  20}; //   1
  psi_needed_g_d[no_ps][22] = get_vector_from_array_int(needed_g_d_22, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__12_48");
  psi_container_IS_name.push_back("phi_of_decay_60__12_48");
  int needed_g_d_23[] = {  21}; //   1
  psi_needed_g_d[no_ps][23] = get_vector_from_array_int(needed_g_d_23, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__36_24");
  psi_container_IS_name.push_back("phi_of_decay_60__36_24");
  int needed_g_d_24[] = {  22}; //   1
  psi_needed_g_d[no_ps][24] = get_vector_from_array_int(needed_g_d_24, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__20_40");
  psi_container_IS_name.push_back("phi_of_decay_60__20_40");

  psi_c_p[no_ps].resize(35);
  psi_c_f[no_ps].resize(35);
  psi_c_t[no_ps].resize(35);
  psi_c_d[no_ps].resize(35);

  int c_p_0[] = {  0}; //   1
  psi_c_p[no_ps][0] = get_vector_from_array_int(c_p_0, 1);
  int c_f_0[] = {  0}; //   1
  psi_c_f[no_ps][0] = get_vector_from_array_int(c_f_0, 1);
  int c_t_0[] = {  0,   1}; //   2
  psi_c_t[no_ps][0] = get_vector_from_array_int(c_t_0, 2);
  int c_d_0[] = {  0}; //   1
  psi_c_d[no_ps][0] = get_vector_from_array_int(c_d_0, 1);

  int c_p_1[] = {  0}; //   1
  psi_c_p[no_ps][1] = get_vector_from_array_int(c_p_1, 1);
  int c_f_1[] = {  1}; //   1
  psi_c_f[no_ps][1] = get_vector_from_array_int(c_f_1, 1);
  int c_t_1[] = {  2,   3}; //   2
  psi_c_t[no_ps][1] = get_vector_from_array_int(c_t_1, 2);
  int c_d_1[] = {  0}; //   1
  psi_c_d[no_ps][1] = get_vector_from_array_int(c_d_1, 1);

  int c_p_2[] = {  1,   2}; //   2
  psi_c_p[no_ps][2] = get_vector_from_array_int(c_p_2, 2);
  int c_f_2[] = {}; //   0
  psi_c_f[no_ps][2] = get_vector_from_array_int(c_f_2, 0);
  int c_t_2[] = {  0}; //   1
  psi_c_t[no_ps][2] = get_vector_from_array_int(c_t_2, 1);
  int c_d_2[] = {  1,   2}; //   2
  psi_c_d[no_ps][2] = get_vector_from_array_int(c_d_2, 2);

  int c_p_3[] = {  3,   4}; //   2
  psi_c_p[no_ps][3] = get_vector_from_array_int(c_p_3, 2);
  int c_f_3[] = {}; //   0
  psi_c_f[no_ps][3] = get_vector_from_array_int(c_f_3, 0);
  int c_t_3[] = {  2}; //   1
  psi_c_t[no_ps][3] = get_vector_from_array_int(c_t_3, 1);
  int c_d_3[] = {  3,   4}; //   2
  psi_c_d[no_ps][3] = get_vector_from_array_int(c_d_3, 2);

  int c_p_4[] = {  0,   5}; //   2
  psi_c_p[no_ps][4] = get_vector_from_array_int(c_p_4, 2);
  int c_f_4[] = {}; //   0
  psi_c_f[no_ps][4] = get_vector_from_array_int(c_f_4, 0);
  int c_t_4[] = {  2}; //   1
  psi_c_t[no_ps][4] = get_vector_from_array_int(c_t_4, 1);
  int c_d_4[] = {  5,   0}; //   2
  psi_c_d[no_ps][4] = get_vector_from_array_int(c_d_4, 2);

  int c_p_5[] = {  0,   6}; //   2
  psi_c_p[no_ps][5] = get_vector_from_array_int(c_p_5, 2);
  int c_f_5[] = {}; //   0
  psi_c_f[no_ps][5] = get_vector_from_array_int(c_f_5, 0);
  int c_t_5[] = {  0}; //   1
  psi_c_t[no_ps][5] = get_vector_from_array_int(c_t_5, 1);
  int c_d_5[] = {  6,   0}; //   2
  psi_c_d[no_ps][5] = get_vector_from_array_int(c_d_5, 2);

  int c_p_6[] = {  7,   8}; //   2
  psi_c_p[no_ps][6] = get_vector_from_array_int(c_p_6, 2);
  int c_f_6[] = {}; //   0
  psi_c_f[no_ps][6] = get_vector_from_array_int(c_f_6, 0);
  int c_t_6[] = {  0}; //   1
  psi_c_t[no_ps][6] = get_vector_from_array_int(c_t_6, 1);
  int c_d_6[] = {  7,   8}; //   2
  psi_c_d[no_ps][6] = get_vector_from_array_int(c_d_6, 2);

  int c_p_7[] = {  9,  10}; //   2
  psi_c_p[no_ps][7] = get_vector_from_array_int(c_p_7, 2);
  int c_f_7[] = {}; //   0
  psi_c_f[no_ps][7] = get_vector_from_array_int(c_f_7, 0);
  int c_t_7[] = {  2}; //   1
  psi_c_t[no_ps][7] = get_vector_from_array_int(c_t_7, 1);
  int c_d_7[] = {  9,  10}; //   2
  psi_c_d[no_ps][7] = get_vector_from_array_int(c_d_7, 2);

  int c_p_8[] = { 11,  12}; //   2
  psi_c_p[no_ps][8] = get_vector_from_array_int(c_p_8, 2);
  int c_f_8[] = {}; //   0
  psi_c_f[no_ps][8] = get_vector_from_array_int(c_f_8, 0);
  int c_t_8[] = {}; //   0
  psi_c_t[no_ps][8] = get_vector_from_array_int(c_t_8, 0);
  int c_d_8[] = { 11,  12,  13}; //   3
  psi_c_d[no_ps][8] = get_vector_from_array_int(c_d_8, 3);

  int c_p_9[] = {  1,  13}; //   2
  psi_c_p[no_ps][9] = get_vector_from_array_int(c_p_9, 2);
  int c_f_9[] = {}; //   0
  psi_c_f[no_ps][9] = get_vector_from_array_int(c_f_9, 0);
  int c_t_9[] = {}; //   0
  psi_c_t[no_ps][9] = get_vector_from_array_int(c_t_9, 0);
  int c_d_9[] = { 11,  14,   2}; //   3
  psi_c_d[no_ps][9] = get_vector_from_array_int(c_d_9, 3);

  int c_p_10[] = {  1,   2}; //   2
  psi_c_p[no_ps][10] = get_vector_from_array_int(c_p_10, 2);
  int c_f_10[] = {}; //   0
  psi_c_f[no_ps][10] = get_vector_from_array_int(c_f_10, 0);
  int c_t_10[] = {}; //   0
  psi_c_t[no_ps][10] = get_vector_from_array_int(c_t_10, 0);
  int c_d_10[] = { 15,   1,   2}; //   3
  psi_c_d[no_ps][10] = get_vector_from_array_int(c_d_10, 3);

  int c_p_11[] = {  3,  14}; //   2
  psi_c_p[no_ps][11] = get_vector_from_array_int(c_p_11, 2);
  int c_f_11[] = {}; //   0
  psi_c_f[no_ps][11] = get_vector_from_array_int(c_f_11, 0);
  int c_t_11[] = {}; //   0
  psi_c_t[no_ps][11] = get_vector_from_array_int(c_t_11, 0);
  int c_d_11[] = { 11,  16,   4}; //   3
  psi_c_d[no_ps][11] = get_vector_from_array_int(c_d_11, 3);

  int c_p_12[] = {  3,   4}; //   2
  psi_c_p[no_ps][12] = get_vector_from_array_int(c_p_12, 2);
  int c_f_12[] = {}; //   0
  psi_c_f[no_ps][12] = get_vector_from_array_int(c_f_12, 0);
  int c_t_12[] = {}; //   0
  psi_c_t[no_ps][12] = get_vector_from_array_int(c_t_12, 0);
  int c_d_12[] = { 17,   3,   4}; //   3
  psi_c_d[no_ps][12] = get_vector_from_array_int(c_d_12, 3);

  int c_p_13[] = {  0,   5}; //   2
  psi_c_p[no_ps][13] = get_vector_from_array_int(c_p_13, 2);
  int c_f_13[] = {}; //   0
  psi_c_f[no_ps][13] = get_vector_from_array_int(c_f_13, 0);
  int c_t_13[] = {}; //   0
  psi_c_t[no_ps][13] = get_vector_from_array_int(c_t_13, 0);
  int c_d_13[] = { 17,   5,   0}; //   3
  psi_c_d[no_ps][13] = get_vector_from_array_int(c_d_13, 3);

  int c_p_14[] = {  0,   6}; //   2
  psi_c_p[no_ps][14] = get_vector_from_array_int(c_p_14, 2);
  int c_f_14[] = {}; //   0
  psi_c_f[no_ps][14] = get_vector_from_array_int(c_f_14, 0);
  int c_t_14[] = {}; //   0
  psi_c_t[no_ps][14] = get_vector_from_array_int(c_t_14, 0);
  int c_d_14[] = { 15,   6,   0}; //   3
  psi_c_d[no_ps][14] = get_vector_from_array_int(c_d_14, 3);

  int c_p_15[] = { 11,  15}; //   2
  psi_c_p[no_ps][15] = get_vector_from_array_int(c_p_15, 2);
  int c_f_15[] = {}; //   0
  psi_c_f[no_ps][15] = get_vector_from_array_int(c_f_15, 0);
  int c_t_15[] = {}; //   0
  psi_c_t[no_ps][15] = get_vector_from_array_int(c_t_15, 0);
  int c_d_15[] = { 18,  19,  13}; //   3
  psi_c_d[no_ps][15] = get_vector_from_array_int(c_d_15, 3);

  int c_p_16[] = {  7,  16}; //   2
  psi_c_p[no_ps][16] = get_vector_from_array_int(c_p_16, 2);
  int c_f_16[] = {}; //   0
  psi_c_f[no_ps][16] = get_vector_from_array_int(c_f_16, 0);
  int c_t_16[] = {}; //   0
  psi_c_t[no_ps][16] = get_vector_from_array_int(c_t_16, 0);
  int c_d_16[] = { 18,  20,   8}; //   3
  psi_c_d[no_ps][16] = get_vector_from_array_int(c_d_16, 3);

  int c_p_17[] = {  7,   8}; //   2
  psi_c_p[no_ps][17] = get_vector_from_array_int(c_p_17, 2);
  int c_f_17[] = {}; //   0
  psi_c_f[no_ps][17] = get_vector_from_array_int(c_f_17, 0);
  int c_t_17[] = {}; //   0
  psi_c_t[no_ps][17] = get_vector_from_array_int(c_t_17, 0);
  int c_d_17[] = { 15,   7,   8}; //   3
  psi_c_d[no_ps][17] = get_vector_from_array_int(c_d_17, 3);

  int c_p_18[] = {  9,  17}; //   2
  psi_c_p[no_ps][18] = get_vector_from_array_int(c_p_18, 2);
  int c_f_18[] = {}; //   0
  psi_c_f[no_ps][18] = get_vector_from_array_int(c_f_18, 0);
  int c_t_18[] = {}; //   0
  psi_c_t[no_ps][18] = get_vector_from_array_int(c_t_18, 0);
  int c_d_18[] = { 18,  21,  10}; //   3
  psi_c_d[no_ps][18] = get_vector_from_array_int(c_d_18, 3);

  int c_p_19[] = {  9,  10}; //   2
  psi_c_p[no_ps][19] = get_vector_from_array_int(c_p_19, 2);
  int c_f_19[] = {}; //   0
  psi_c_f[no_ps][19] = get_vector_from_array_int(c_f_19, 0);
  int c_t_19[] = {}; //   0
  psi_c_t[no_ps][19] = get_vector_from_array_int(c_t_19, 0);
  int c_d_19[] = { 17,   9,  10}; //   3
  psi_c_d[no_ps][19] = get_vector_from_array_int(c_d_19, 3);

  int c_p_20[] = { 11,  18}; //   2
  psi_c_p[no_ps][20] = get_vector_from_array_int(c_p_20, 2);
  int c_f_20[] = {}; //   0
  psi_c_f[no_ps][20] = get_vector_from_array_int(c_f_20, 0);
  int c_t_20[] = {}; //   0
  psi_c_t[no_ps][20] = get_vector_from_array_int(c_t_20, 0);
  int c_d_20[] = { 22,   0,  13}; //   3
  psi_c_d[no_ps][20] = get_vector_from_array_int(c_d_20, 3);

  int c_p_21[] = {  1,  19}; //   2
  psi_c_p[no_ps][21] = get_vector_from_array_int(c_p_21, 2);
  int c_f_21[] = {}; //   0
  psi_c_f[no_ps][21] = get_vector_from_array_int(c_f_21, 0);
  int c_t_21[] = {}; //   0
  psi_c_t[no_ps][21] = get_vector_from_array_int(c_t_21, 0);
  int c_d_21[] = { 23,  10,   2}; //   3
  psi_c_d[no_ps][21] = get_vector_from_array_int(c_d_21, 3);

  int c_p_22[] = {  3,  20}; //   2
  psi_c_p[no_ps][22] = get_vector_from_array_int(c_p_22, 2);
  int c_f_22[] = {}; //   0
  psi_c_f[no_ps][22] = get_vector_from_array_int(c_f_22, 0);
  int c_t_22[] = {}; //   0
  psi_c_t[no_ps][22] = get_vector_from_array_int(c_t_22, 0);
  int c_d_22[] = { 24,   8,   4}; //   3
  psi_c_d[no_ps][22] = get_vector_from_array_int(c_d_22, 3);

  int c_p_23[] = {  0}; //   1
  psi_c_p[no_ps][23] = get_vector_from_array_int(c_p_23, 1);
  int c_f_23[] = {  0}; //   1
  psi_c_f[no_ps][23] = get_vector_from_array_int(c_f_23, 1);
  int c_t_23[] = {  4,   5}; //   2
  psi_c_t[no_ps][23] = get_vector_from_array_int(c_t_23, 2);
  int c_d_23[] = {  0}; //   1
  psi_c_d[no_ps][23] = get_vector_from_array_int(c_d_23, 1);

  int c_p_24[] = {  0}; //   1
  psi_c_p[no_ps][24] = get_vector_from_array_int(c_p_24, 1);
  int c_f_24[] = {  1}; //   1
  psi_c_f[no_ps][24] = get_vector_from_array_int(c_f_24, 1);
  int c_t_24[] = {  6,   7}; //   2
  psi_c_t[no_ps][24] = get_vector_from_array_int(c_t_24, 2);
  int c_d_24[] = {  0}; //   1
  psi_c_d[no_ps][24] = get_vector_from_array_int(c_d_24, 1);

  int c_p_25[] = {  1,   2}; //   2
  psi_c_p[no_ps][25] = get_vector_from_array_int(c_p_25, 2);
  int c_f_25[] = {}; //   0
  psi_c_f[no_ps][25] = get_vector_from_array_int(c_f_25, 0);
  int c_t_25[] = {  4}; //   1
  psi_c_t[no_ps][25] = get_vector_from_array_int(c_t_25, 1);
  int c_d_25[] = {  1,   2}; //   2
  psi_c_d[no_ps][25] = get_vector_from_array_int(c_d_25, 2);

  int c_p_26[] = {  3,   4}; //   2
  psi_c_p[no_ps][26] = get_vector_from_array_int(c_p_26, 2);
  int c_f_26[] = {}; //   0
  psi_c_f[no_ps][26] = get_vector_from_array_int(c_f_26, 0);
  int c_t_26[] = {  6}; //   1
  psi_c_t[no_ps][26] = get_vector_from_array_int(c_t_26, 1);
  int c_d_26[] = {  3,   4}; //   2
  psi_c_d[no_ps][26] = get_vector_from_array_int(c_d_26, 2);

  int c_p_27[] = {  0,   5}; //   2
  psi_c_p[no_ps][27] = get_vector_from_array_int(c_p_27, 2);
  int c_f_27[] = {}; //   0
  psi_c_f[no_ps][27] = get_vector_from_array_int(c_f_27, 0);
  int c_t_27[] = {  6}; //   1
  psi_c_t[no_ps][27] = get_vector_from_array_int(c_t_27, 1);
  int c_d_27[] = {  5,   0}; //   2
  psi_c_d[no_ps][27] = get_vector_from_array_int(c_d_27, 2);

  int c_p_28[] = {  0,   6}; //   2
  psi_c_p[no_ps][28] = get_vector_from_array_int(c_p_28, 2);
  int c_f_28[] = {}; //   0
  psi_c_f[no_ps][28] = get_vector_from_array_int(c_f_28, 0);
  int c_t_28[] = {  4}; //   1
  psi_c_t[no_ps][28] = get_vector_from_array_int(c_t_28, 1);
  int c_d_28[] = {  6,   0}; //   2
  psi_c_d[no_ps][28] = get_vector_from_array_int(c_d_28, 2);

  int c_p_29[] = {  7,   8}; //   2
  psi_c_p[no_ps][29] = get_vector_from_array_int(c_p_29, 2);
  int c_f_29[] = {}; //   0
  psi_c_f[no_ps][29] = get_vector_from_array_int(c_f_29, 0);
  int c_t_29[] = {  4}; //   1
  psi_c_t[no_ps][29] = get_vector_from_array_int(c_t_29, 1);
  int c_d_29[] = {  7,   8}; //   2
  psi_c_d[no_ps][29] = get_vector_from_array_int(c_d_29, 2);

  int c_p_30[] = {  9,  10}; //   2
  psi_c_p[no_ps][30] = get_vector_from_array_int(c_p_30, 2);
  int c_f_30[] = {}; //   0
  psi_c_f[no_ps][30] = get_vector_from_array_int(c_f_30, 0);
  int c_t_30[] = {  6}; //   1
  psi_c_t[no_ps][30] = get_vector_from_array_int(c_t_30, 1);
  int c_d_30[] = {  9,  10}; //   2
  psi_c_d[no_ps][30] = get_vector_from_array_int(c_d_30, 2);

  int c_p_31[] = {  0}; //   1
  psi_c_p[no_ps][31] = get_vector_from_array_int(c_p_31, 1);
  int c_f_31[] = {  1}; //   1
  psi_c_f[no_ps][31] = get_vector_from_array_int(c_f_31, 1);
  int c_t_31[] = {  6,   8}; //   2
  psi_c_t[no_ps][31] = get_vector_from_array_int(c_t_31, 2);
  int c_d_31[] = {  0}; //   1
  psi_c_d[no_ps][31] = get_vector_from_array_int(c_d_31, 1);

  int c_p_32[] = {  0}; //   1
  psi_c_p[no_ps][32] = get_vector_from_array_int(c_p_32, 1);
  int c_f_32[] = {  0}; //   1
  psi_c_f[no_ps][32] = get_vector_from_array_int(c_f_32, 1);
  int c_t_32[] = {  4,   9}; //   2
  psi_c_t[no_ps][32] = get_vector_from_array_int(c_t_32, 2);
  int c_d_32[] = {  0}; //   1
  psi_c_d[no_ps][32] = get_vector_from_array_int(c_d_32, 1);

  int c_p_33[] = { 11,  18}; //   2
  psi_c_p[no_ps][33] = get_vector_from_array_int(c_p_33, 2);
  int c_f_33[] = {}; //   0
  psi_c_f[no_ps][33] = get_vector_from_array_int(c_f_33, 0);
  int c_t_33[] = { 10}; //   1
  psi_c_t[no_ps][33] = get_vector_from_array_int(c_t_33, 1);
  int c_d_33[] = {  0,  13}; //   2
  psi_c_d[no_ps][33] = get_vector_from_array_int(c_d_33, 2);

  int c_p_34[] = {  0,  21}; //   2
  psi_c_p[no_ps][34] = get_vector_from_array_int(c_p_34, 2);
  int c_f_34[] = {}; //   0
  psi_c_f[no_ps][34] = get_vector_from_array_int(c_f_34, 0);
  int c_t_34[] = { 11}; //   1
  psi_c_t[no_ps][34] = get_vector_from_array_int(c_t_34, 1);
  int c_d_34[] = { 13,   0}; //   2
  psi_c_d[no_ps][34] = get_vector_from_array_int(c_d_34, 2);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
