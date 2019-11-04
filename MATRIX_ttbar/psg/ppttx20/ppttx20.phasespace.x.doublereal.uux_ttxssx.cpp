#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ax_psp_400_uux_ttxssx(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_400_uux_ttxssx");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  psi_v_smin[no_ps].resize(6);
  psi_needed_v_smin[no_ps].resize(6);
  int needed_v_smin_0[] = {   0,    3,    4,    5,    6}; //   5
  psi_needed_v_smin[no_ps][  0] = get_vector_from_array_int(needed_v_smin_0, 5);
  int needed_v_smin_1[] = {   0}; //   1
  psi_needed_v_smin[no_ps][  1] = get_vector_from_array_int(needed_v_smin_1, 1);
  int needed_v_smin_2[] = {   1,    2,    4,    5,    6}; //   5
  psi_needed_v_smin[no_ps][  2] = get_vector_from_array_int(needed_v_smin_2, 5);
  int needed_v_smin_3[] = {   1}; //   1
  psi_needed_v_smin[no_ps][  3] = get_vector_from_array_int(needed_v_smin_3, 1);
  int needed_v_smin_4[] = {   2}; //   1
  psi_needed_v_smin[no_ps][  4] = get_vector_from_array_int(needed_v_smin_4, 1);
  int needed_v_smin_5[] = {   3}; //   1
  psi_needed_v_smin[no_ps][  5] = get_vector_from_array_int(needed_v_smin_5, 1);

  psi_v_smax[no_ps].resize(8);
  psi_needed_v_smax[no_ps].resize(8);
  int needed_v_smax_0[] = {   0,    3,    4,    5}; //   4
  psi_needed_v_smax[no_ps][  0] = get_vector_from_array_int(needed_v_smax_0, 4);
  int needed_v_smax_1[] = {   0}; //   1
  psi_needed_v_smax[no_ps][  1] = get_vector_from_array_int(needed_v_smax_1, 1);
  int needed_v_smax_2[] = {   1,    2,    6}; //   3
  psi_needed_v_smax[no_ps][  2] = get_vector_from_array_int(needed_v_smax_2, 3);
  int needed_v_smax_3[] = {   1}; //   1
  psi_needed_v_smax[no_ps][  3] = get_vector_from_array_int(needed_v_smax_3, 1);
  int needed_v_smax_4[] = {   2}; //   1
  psi_needed_v_smax[no_ps][  4] = get_vector_from_array_int(needed_v_smax_4, 1);
  int needed_v_smax_5[] = {   3}; //   1
  psi_needed_v_smax[no_ps][  5] = get_vector_from_array_int(needed_v_smax_5, 1);
  int needed_v_smax_6[] = {   4,    5}; //   2
  psi_needed_v_smax[no_ps][  6] = get_vector_from_array_int(needed_v_smax_6, 2);
  int needed_v_smax_7[] = {   6}; //   1
  psi_needed_v_smax[no_ps][  7] = get_vector_from_array_int(needed_v_smax_7, 1);

  psi_container_IS_startvalue[no_ps][0] = psi_container_IS_name.size();
  psi_g_p[no_ps].resize(8);
  psi_needed_g_p[no_ps].resize(8);
  int needed_g_p_0[] = {   0,    3,    4,    5}; //   4
  psi_needed_g_p[no_ps][0] = get_vector_from_array_int(needed_g_p_0, 4);
  psi_container_IS_name.push_back("propagator_0_smin0_smax0");
  int needed_g_p_1[] = {   0}; //   1
  psi_needed_g_p[no_ps][1] = get_vector_from_array_int(needed_g_p_1, 1);
  psi_container_IS_name.push_back("propagator_6_smin1_smax1");
  int needed_g_p_2[] = {   1,    2,    6}; //   3
  psi_needed_g_p[no_ps][2] = get_vector_from_array_int(needed_g_p_2, 3);
  psi_container_IS_name.push_back("propagator_0_smin2_smax2");
  int needed_g_p_3[] = {   1}; //   1
  psi_needed_g_p[no_ps][3] = get_vector_from_array_int(needed_g_p_3, 1);
  psi_container_IS_name.push_back("propagator_0_smin3_smax3");
  int needed_g_p_4[] = {   2}; //   1
  psi_needed_g_p[no_ps][4] = get_vector_from_array_int(needed_g_p_4, 1);
  psi_container_IS_name.push_back("propagator_0_smin4_smax4");
  int needed_g_p_5[] = {   3}; //   1
  psi_needed_g_p[no_ps][5] = get_vector_from_array_int(needed_g_p_5, 1);
  psi_container_IS_name.push_back("propagator_6_smin5_smax5");
  int needed_g_p_6[] = {   4,    5}; //   2
  psi_needed_g_p[no_ps][6] = get_vector_from_array_int(needed_g_p_6, 2);
  psi_container_IS_name.push_back("propagator_0_smin2_smax6");
  int needed_g_p_7[] = {   6}; //   1
  psi_needed_g_p[no_ps][7] = get_vector_from_array_int(needed_g_p_7, 1);
  psi_container_IS_name.push_back("propagator_0_smin0_smax7");

  psi_container_IS_startvalue[no_ps][1] = psi_container_IS_name.size();
  psi_g_f[no_ps].resize(0);
  psi_needed_g_f[no_ps].resize(0);

  psi_container_IS_startvalue[no_ps][2] = psi_container_IS_name.size();
  psi_g_t[no_ps].resize(2);
  psi_needed_g_t[no_ps].resize(2);
  int needed_g_t_0[] = {   5}; //   1
  psi_needed_g_t[no_ps][0] = get_vector_from_array_int(needed_g_t_0, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__12_48");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__12_48");
  int needed_g_t_1[] = {   6}; //   1
  psi_needed_g_t[no_ps][1] = get_vector_from_array_int(needed_g_t_1, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__48_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__48_12");

  psi_container_IS_startvalue[no_ps][3] = psi_container_IS_name.size();
  psi_g_d[no_ps].resize(11);
  psi_needed_g_d[no_ps].resize(11);
  int needed_g_d_0[] = {   0}; //   1
  psi_needed_g_d[no_ps][0] = get_vector_from_array_int(needed_g_d_0, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__4_56");
  psi_container_IS_name.push_back("phi_of_decay_60__4_56");
  int needed_g_d_1[] = {   0}; //   1
  psi_needed_g_d[no_ps][1] = get_vector_from_array_int(needed_g_d_1, 1);
  psi_container_IS_name.push_back("costheta_of_decay_56__8_48");
  psi_container_IS_name.push_back("phi_of_decay_56__8_48");
  int needed_g_d_2[] = {   0,    3,    4,    5,    6}; //   5
  psi_needed_g_d[no_ps][2] = get_vector_from_array_int(needed_g_d_2, 5);
  psi_container_IS_name.push_back("costheta_of_decay_48__16_32");
  psi_container_IS_name.push_back("phi_of_decay_48__16_32");
  int needed_g_d_3[] = {   1}; //   1
  psi_needed_g_d[no_ps][3] = get_vector_from_array_int(needed_g_d_3, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__16_44");
  psi_container_IS_name.push_back("phi_of_decay_60__16_44");
  int needed_g_d_4[] = {   1}; //   1
  psi_needed_g_d[no_ps][4] = get_vector_from_array_int(needed_g_d_4, 1);
  psi_container_IS_name.push_back("costheta_of_decay_44__32_12");
  psi_container_IS_name.push_back("phi_of_decay_44__32_12");
  int needed_g_d_5[] = {   1,    2,    4,    5,    6}; //   5
  psi_needed_g_d[no_ps][5] = get_vector_from_array_int(needed_g_d_5, 5);
  psi_container_IS_name.push_back("costheta_of_decay_12__4_8");
  psi_container_IS_name.push_back("phi_of_decay_12__4_8");
  int needed_g_d_6[] = {   2}; //   1
  psi_needed_g_d[no_ps][6] = get_vector_from_array_int(needed_g_d_6, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__32_28");
  psi_container_IS_name.push_back("phi_of_decay_60__32_28");
  int needed_g_d_7[] = {   2}; //   1
  psi_needed_g_d[no_ps][7] = get_vector_from_array_int(needed_g_d_7, 1);
  psi_container_IS_name.push_back("costheta_of_decay_28__16_12");
  psi_container_IS_name.push_back("phi_of_decay_28__16_12");
  int needed_g_d_8[] = {   3}; //   1
  psi_needed_g_d[no_ps][8] = get_vector_from_array_int(needed_g_d_8, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__8_52");
  psi_container_IS_name.push_back("phi_of_decay_60__8_52");
  int needed_g_d_9[] = {   3}; //   1
  psi_needed_g_d[no_ps][9] = get_vector_from_array_int(needed_g_d_9, 1);
  psi_container_IS_name.push_back("costheta_of_decay_52__4_48");
  psi_container_IS_name.push_back("phi_of_decay_52__4_48");
  int needed_g_d_10[] = {   4}; //   1
  psi_needed_g_d[no_ps][10] = get_vector_from_array_int(needed_g_d_10, 1);
  psi_container_IS_name.push_back("costheta_of_decay_60__12_48");
  psi_container_IS_name.push_back("phi_of_decay_60__12_48");

  psi_c_p[no_ps].resize(7);
  psi_c_f[no_ps].resize(7);
  psi_c_t[no_ps].resize(7);
  psi_c_d[no_ps].resize(7);

  int c_p_0[] = {  0,   1}; //   2
  psi_c_p[no_ps][0] = get_vector_from_array_int(c_p_0, 2);
  int c_f_0[] = {}; //   0
  psi_c_f[no_ps][0] = get_vector_from_array_int(c_f_0, 0);
  int c_t_0[] = {}; //   0
  psi_c_t[no_ps][0] = get_vector_from_array_int(c_t_0, 0);
  int c_d_0[] = {  0,   1,   2}; //   3
  psi_c_d[no_ps][0] = get_vector_from_array_int(c_d_0, 3);

  int c_p_1[] = {  2,   3}; //   2
  psi_c_p[no_ps][1] = get_vector_from_array_int(c_p_1, 2);
  int c_f_1[] = {}; //   0
  psi_c_f[no_ps][1] = get_vector_from_array_int(c_f_1, 0);
  int c_t_1[] = {}; //   0
  psi_c_t[no_ps][1] = get_vector_from_array_int(c_t_1, 0);
  int c_d_1[] = {  3,   4,   5}; //   3
  psi_c_d[no_ps][1] = get_vector_from_array_int(c_d_1, 3);

  int c_p_2[] = {  2,   4}; //   2
  psi_c_p[no_ps][2] = get_vector_from_array_int(c_p_2, 2);
  int c_f_2[] = {}; //   0
  psi_c_f[no_ps][2] = get_vector_from_array_int(c_f_2, 0);
  int c_t_2[] = {}; //   0
  psi_c_t[no_ps][2] = get_vector_from_array_int(c_t_2, 0);
  int c_d_2[] = {  6,   7,   5}; //   3
  psi_c_d[no_ps][2] = get_vector_from_array_int(c_d_2, 3);

  int c_p_3[] = {  0,   5}; //   2
  psi_c_p[no_ps][3] = get_vector_from_array_int(c_p_3, 2);
  int c_f_3[] = {}; //   0
  psi_c_f[no_ps][3] = get_vector_from_array_int(c_f_3, 0);
  int c_t_3[] = {}; //   0
  psi_c_t[no_ps][3] = get_vector_from_array_int(c_t_3, 0);
  int c_d_3[] = {  8,   9,   2}; //   3
  psi_c_d[no_ps][3] = get_vector_from_array_int(c_d_3, 3);

  int c_p_4[] = {  0,   6}; //   2
  psi_c_p[no_ps][4] = get_vector_from_array_int(c_p_4, 2);
  int c_f_4[] = {}; //   0
  psi_c_f[no_ps][4] = get_vector_from_array_int(c_f_4, 0);
  int c_t_4[] = {}; //   0
  psi_c_t[no_ps][4] = get_vector_from_array_int(c_t_4, 0);
  int c_d_4[] = { 10,   5,   2}; //   3
  psi_c_d[no_ps][4] = get_vector_from_array_int(c_d_4, 3);

  int c_p_5[] = {  0,   6}; //   2
  psi_c_p[no_ps][5] = get_vector_from_array_int(c_p_5, 2);
  int c_f_5[] = {}; //   0
  psi_c_f[no_ps][5] = get_vector_from_array_int(c_f_5, 0);
  int c_t_5[] = {  0}; //   1
  psi_c_t[no_ps][5] = get_vector_from_array_int(c_t_5, 1);
  int c_d_5[] = {  5,   2}; //   2
  psi_c_d[no_ps][5] = get_vector_from_array_int(c_d_5, 2);

  int c_p_6[] = {  2,   7}; //   2
  psi_c_p[no_ps][6] = get_vector_from_array_int(c_p_6, 2);
  int c_f_6[] = {}; //   0
  psi_c_f[no_ps][6] = get_vector_from_array_int(c_f_6, 0);
  int c_t_6[] = {  1}; //   1
  psi_c_t[no_ps][6] = get_vector_from_array_int(c_t_6, 1);
  int c_d_6[] = {  2,   5}; //   2
  psi_c_d[no_ps][6] = get_vector_from_array_int(c_d_6, 2);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
