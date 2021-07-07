#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ax_psp_300_bbx_ttxg(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_300_bbx_ttxg");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  psi_v_smin[no_ps].resize(3);
  psi_needed_v_smin[no_ps].resize(3);
  int needed_v_smin_0[] = {   0,    2,    4}; //   3
  psi_needed_v_smin[no_ps][  0] = get_vector_from_array_int(needed_v_smin_0, 3);
  int needed_v_smin_1[] = {   1}; //   1
  psi_needed_v_smin[no_ps][  1] = get_vector_from_array_int(needed_v_smin_1, 1);
  int needed_v_smin_2[] = {   3}; //   1
  psi_needed_v_smin[no_ps][  2] = get_vector_from_array_int(needed_v_smin_2, 1);

  psi_v_smax[no_ps].resize(3);
  psi_needed_v_smax[no_ps].resize(3);
  int needed_v_smax_0[] = {   0,    2,    4}; //   3
  psi_needed_v_smax[no_ps][  0] = get_vector_from_array_int(needed_v_smax_0, 3);
  int needed_v_smax_1[] = {   1}; //   1
  psi_needed_v_smax[no_ps][  1] = get_vector_from_array_int(needed_v_smax_1, 1);
  int needed_v_smax_2[] = {   3}; //   1
  psi_needed_v_smax[no_ps][  2] = get_vector_from_array_int(needed_v_smax_2, 1);

  psi_container_IS_startvalue[no_ps][0] = psi_container_IS_name.size();
  psi_g_p[no_ps].resize(3);
  psi_needed_g_p[no_ps].resize(3);
  int needed_g_p_0[] = {   0,    2,    4}; //   3
  psi_needed_g_p[no_ps][0] = get_vector_from_array_int(needed_g_p_0, 3);
  psi_container_IS_name.push_back("propagator_0_smin0_smax0");
  int needed_g_p_1[] = {   1}; //   1
  psi_needed_g_p[no_ps][1] = get_vector_from_array_int(needed_g_p_1, 1);
  psi_container_IS_name.push_back("propagator_6_smin1_smax1");
  int needed_g_p_2[] = {   3}; //   1
  psi_needed_g_p[no_ps][2] = get_vector_from_array_int(needed_g_p_2, 1);
  psi_container_IS_name.push_back("propagator_6_smin2_smax2");

  psi_container_IS_startvalue[no_ps][1] = psi_container_IS_name.size();
  psi_g_f[no_ps].resize(0);
  psi_needed_g_f[no_ps].resize(0);

  psi_container_IS_startvalue[no_ps][2] = psi_container_IS_name.size();
  psi_g_t[no_ps].resize(2);
  psi_needed_g_t[no_ps].resize(2);
  int needed_g_t_0[] = {   0}; //   1
  psi_needed_g_t[no_ps][0] = get_vector_from_array_int(needed_g_t_0, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__2_1__16_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__2_1__16_12");
  int needed_g_t_1[] = {   4}; //   1
  psi_needed_g_t[no_ps][1] = get_vector_from_array_int(needed_g_t_1, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_5__1_2__16_12");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_5__1_2__16_12");

  psi_container_IS_startvalue[no_ps][3] = psi_container_IS_name.size();
  psi_g_d[no_ps].resize(6);
  psi_needed_g_d[no_ps].resize(6);
  int needed_g_d_0[] = {   0,    2,    4}; //   3
  psi_needed_g_d[no_ps][0] = get_vector_from_array_int(needed_g_d_0, 3);
  psi_container_IS_name.push_back("costheta_of_decay_12__4_8");
  psi_container_IS_name.push_back("phi_of_decay_12__4_8");
  int needed_g_d_1[] = {   1}; //   1
  psi_needed_g_d[no_ps][1] = get_vector_from_array_int(needed_g_d_1, 1);
  psi_container_IS_name.push_back("costheta_of_decay_28__4_24");
  psi_container_IS_name.push_back("phi_of_decay_28__4_24");
  int needed_g_d_2[] = {   1}; //   1
  psi_needed_g_d[no_ps][2] = get_vector_from_array_int(needed_g_d_2, 1);
  psi_container_IS_name.push_back("costheta_of_decay_24__8_16");
  psi_container_IS_name.push_back("phi_of_decay_24__8_16");
  int needed_g_d_3[] = {   2}; //   1
  psi_needed_g_d[no_ps][3] = get_vector_from_array_int(needed_g_d_3, 1);
  psi_container_IS_name.push_back("costheta_of_decay_28__16_12");
  psi_container_IS_name.push_back("phi_of_decay_28__16_12");
  int needed_g_d_4[] = {   3}; //   1
  psi_needed_g_d[no_ps][4] = get_vector_from_array_int(needed_g_d_4, 1);
  psi_container_IS_name.push_back("costheta_of_decay_28__8_20");
  psi_container_IS_name.push_back("phi_of_decay_28__8_20");
  int needed_g_d_5[] = {   3}; //   1
  psi_needed_g_d[no_ps][5] = get_vector_from_array_int(needed_g_d_5, 1);
  psi_container_IS_name.push_back("costheta_of_decay_20__4_16");
  psi_container_IS_name.push_back("phi_of_decay_20__4_16");

  psi_c_p[no_ps].resize(5);
  psi_c_f[no_ps].resize(5);
  psi_c_t[no_ps].resize(5);
  psi_c_d[no_ps].resize(5);

  int c_p_0[] = {  0}; //   1
  psi_c_p[no_ps][0] = get_vector_from_array_int(c_p_0, 1);
  int c_f_0[] = {}; //   0
  psi_c_f[no_ps][0] = get_vector_from_array_int(c_f_0, 0);
  int c_t_0[] = {  0}; //   1
  psi_c_t[no_ps][0] = get_vector_from_array_int(c_t_0, 1);
  int c_d_0[] = {  0}; //   1
  psi_c_d[no_ps][0] = get_vector_from_array_int(c_d_0, 1);

  int c_p_1[] = {  1}; //   1
  psi_c_p[no_ps][1] = get_vector_from_array_int(c_p_1, 1);
  int c_f_1[] = {}; //   0
  psi_c_f[no_ps][1] = get_vector_from_array_int(c_f_1, 0);
  int c_t_1[] = {}; //   0
  psi_c_t[no_ps][1] = get_vector_from_array_int(c_t_1, 0);
  int c_d_1[] = {  1,   2}; //   2
  psi_c_d[no_ps][1] = get_vector_from_array_int(c_d_1, 2);

  int c_p_2[] = {  0}; //   1
  psi_c_p[no_ps][2] = get_vector_from_array_int(c_p_2, 1);
  int c_f_2[] = {}; //   0
  psi_c_f[no_ps][2] = get_vector_from_array_int(c_f_2, 0);
  int c_t_2[] = {}; //   0
  psi_c_t[no_ps][2] = get_vector_from_array_int(c_t_2, 0);
  int c_d_2[] = {  3,   0}; //   2
  psi_c_d[no_ps][2] = get_vector_from_array_int(c_d_2, 2);

  int c_p_3[] = {  2}; //   1
  psi_c_p[no_ps][3] = get_vector_from_array_int(c_p_3, 1);
  int c_f_3[] = {}; //   0
  psi_c_f[no_ps][3] = get_vector_from_array_int(c_f_3, 0);
  int c_t_3[] = {}; //   0
  psi_c_t[no_ps][3] = get_vector_from_array_int(c_t_3, 0);
  int c_d_3[] = {  4,   5}; //   2
  psi_c_d[no_ps][3] = get_vector_from_array_int(c_d_3, 2);

  int c_p_4[] = {  0}; //   1
  psi_c_p[no_ps][4] = get_vector_from_array_int(c_p_4, 1);
  int c_f_4[] = {}; //   0
  psi_c_f[no_ps][4] = get_vector_from_array_int(c_f_4, 0);
  int c_t_4[] = {  1}; //   1
  psi_c_t[no_ps][4] = get_vector_from_array_int(c_t_4, 1);
  int c_d_4[] = {  0}; //   1
  psi_c_d[no_ps][4] = get_vector_from_array_int(c_d_4, 1);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
