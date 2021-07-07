#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ax_psp_300_gb_ttxb(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_300_gb_ttxb");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  psi_v_smin[no_ps].resize(3);
  psi_needed_v_smin[no_ps].resize(3);
  int needed_v_smin_0[] = {   2,    3,    4}; //   3
  psi_needed_v_smin[no_ps][  0] = get_vector_from_array_int(needed_v_smin_0, 3);
  int needed_v_smin_1[] = {   0}; //   1
  psi_needed_v_smin[no_ps][  1] = get_vector_from_array_int(needed_v_smin_1, 1);
  int needed_v_smin_2[] = {   1}; //   1
  psi_needed_v_smin[no_ps][  2] = get_vector_from_array_int(needed_v_smin_2, 1);

  psi_v_smax[no_ps].resize(3);
  psi_needed_v_smax[no_ps].resize(3);
  int needed_v_smax_0[] = {   2,    3,    4}; //   3
  psi_needed_v_smax[no_ps][  0] = get_vector_from_array_int(needed_v_smax_0, 3);
  int needed_v_smax_1[] = {   0}; //   1
  psi_needed_v_smax[no_ps][  1] = get_vector_from_array_int(needed_v_smax_1, 1);
  int needed_v_smax_2[] = {   1}; //   1
  psi_needed_v_smax[no_ps][  2] = get_vector_from_array_int(needed_v_smax_2, 1);

  psi_container_IS_startvalue[no_ps][0] = psi_container_IS_name.size();
  psi_g_p[no_ps].resize(1);
  psi_needed_g_p[no_ps].resize(1);
  int needed_g_p_0[] = {   2,    3,    4}; //   3
  psi_needed_g_p[no_ps][0] = get_vector_from_array_int(needed_g_p_0, 3);
  psi_container_IS_name.push_back("propagator_0_smin0_smax0");

  psi_container_IS_startvalue[no_ps][1] = psi_container_IS_name.size();
  psi_g_f[no_ps].resize(2);
  psi_needed_g_f[no_ps].resize(2);
  int needed_g_f_0[] = {   0}; //   1
  psi_needed_g_f[no_ps][0] = get_vector_from_array_int(needed_g_f_0, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin1_smax1");
  int needed_g_f_1[] = {   1}; //   1
  psi_needed_g_f[no_ps][1] = get_vector_from_array_int(needed_g_f_1, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin2_smax2");

  psi_container_IS_startvalue[no_ps][2] = psi_container_IS_name.size();
  psi_g_t[no_ps].resize(6);
  psi_needed_g_t[no_ps].resize(6);
  int needed_g_t_0[] = {   0}; //   1
  psi_needed_g_t[no_ps][0] = get_vector_from_array_int(needed_g_t_0, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__1_2__4_24");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__1_2__4_24");
  int needed_g_t_1[] = {   0}; //   1
  psi_needed_g_t[no_ps][1] = get_vector_from_array_int(needed_g_t_1, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_5__16_8");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_5__16_8");
  int needed_g_t_2[] = {   1}; //   1
  psi_needed_g_t[no_ps][2] = get_vector_from_array_int(needed_g_t_2, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_6__1_2__8_20");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_6__1_2__8_20");
  int needed_g_t_3[] = {   1}; //   1
  psi_needed_g_t[no_ps][3] = get_vector_from_array_int(needed_g_t_3, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_9__16_4");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_9__16_4");
  int needed_g_t_4[] = {   2}; //   1
  psi_needed_g_t[no_ps][4] = get_vector_from_array_int(needed_g_t_4, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_1__16_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_1__16_12");
  int needed_g_t_5[] = {   4}; //   1
  psi_needed_g_t[no_ps][5] = get_vector_from_array_int(needed_g_t_5, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__1_2__16_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__1_2__16_12");

  psi_container_IS_startvalue[no_ps][3] = psi_container_IS_name.size();
  psi_g_d[no_ps].resize(2);
  psi_needed_g_d[no_ps].resize(2);
  int needed_g_d_0[] = {   2,    3,    4}; //   3
  psi_needed_g_d[no_ps][0] = get_vector_from_array_int(needed_g_d_0, 3);
  psi_container_IS_name.push_back("costheta_of_decay_12__4_8");
  psi_container_IS_name.push_back("phi_of_decay_12__4_8");
  int needed_g_d_1[] = {   3}; //   1
  psi_needed_g_d[no_ps][1] = get_vector_from_array_int(needed_g_d_1, 1);
  psi_container_IS_name.push_back("costheta_of_decay_28__16_12");
  psi_container_IS_name.push_back("phi_of_decay_28__16_12");

  psi_c_p[no_ps].resize(5);
  psi_c_f[no_ps].resize(5);
  psi_c_t[no_ps].resize(5);
  psi_c_d[no_ps].resize(5);

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
  int c_f_1[] = {  1}; //   1
  psi_c_f[no_ps][1] = get_vector_from_array_int(c_f_1, 1);
  int c_t_1[] = {  2,   3}; //   2
  psi_c_t[no_ps][1] = get_vector_from_array_int(c_t_1, 2);
  int c_d_1[] = {}; //   0
  psi_c_d[no_ps][1] = get_vector_from_array_int(c_d_1, 0);

  int c_p_2[] = {  0}; //   1
  psi_c_p[no_ps][2] = get_vector_from_array_int(c_p_2, 1);
  int c_f_2[] = {}; //   0
  psi_c_f[no_ps][2] = get_vector_from_array_int(c_f_2, 0);
  int c_t_2[] = {  4}; //   1
  psi_c_t[no_ps][2] = get_vector_from_array_int(c_t_2, 1);
  int c_d_2[] = {  0}; //   1
  psi_c_d[no_ps][2] = get_vector_from_array_int(c_d_2, 1);

  int c_p_3[] = {  0}; //   1
  psi_c_p[no_ps][3] = get_vector_from_array_int(c_p_3, 1);
  int c_f_3[] = {}; //   0
  psi_c_f[no_ps][3] = get_vector_from_array_int(c_f_3, 0);
  int c_t_3[] = {}; //   0
  psi_c_t[no_ps][3] = get_vector_from_array_int(c_t_3, 0);
  int c_d_3[] = {  1,   0}; //   2
  psi_c_d[no_ps][3] = get_vector_from_array_int(c_d_3, 2);

  int c_p_4[] = {  0}; //   1
  psi_c_p[no_ps][4] = get_vector_from_array_int(c_p_4, 1);
  int c_f_4[] = {}; //   0
  psi_c_f[no_ps][4] = get_vector_from_array_int(c_f_4, 0);
  int c_t_4[] = {  5}; //   1
  psi_c_t[no_ps][4] = get_vector_from_array_int(c_t_4, 1);
  int c_d_4[] = {  0}; //   1
  psi_c_d[no_ps][4] = get_vector_from_array_int(c_d_4, 1);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
