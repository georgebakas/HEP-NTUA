#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ax_psp_400_dxsx_ttxdxsx(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_400_dxsx_ttxdxsx");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  psi_v_smin[no_ps].resize(3);
  psi_needed_v_smin[no_ps].resize(3);
  int needed_v_smin_0[] = {   2,    3,    4,    5,    6,    0,    1}; //   7
  psi_needed_v_smin[no_ps][  0] = get_vector_from_array_int(needed_v_smin_0, 7);
  int needed_v_smin_1[] = {   3,    2}; //   2
  psi_needed_v_smin[no_ps][  1] = get_vector_from_array_int(needed_v_smin_1, 2);
  int needed_v_smin_2[] = {   5,    0,    1,    4,    6}; //   5
  psi_needed_v_smin[no_ps][  2] = get_vector_from_array_int(needed_v_smin_2, 5);

  psi_v_smax[no_ps].resize(3);
  psi_needed_v_smax[no_ps].resize(3);
  int needed_v_smax_0[] = {   2,    3,    4,    5,    6,    0,    1}; //   7
  psi_needed_v_smax[no_ps][  0] = get_vector_from_array_int(needed_v_smax_0, 7);
  int needed_v_smax_1[] = {   3,    2}; //   2
  psi_needed_v_smax[no_ps][  1] = get_vector_from_array_int(needed_v_smax_1, 2);
  int needed_v_smax_2[] = {   5,    0,    1,    4,    6}; //   5
  psi_needed_v_smax[no_ps][  2] = get_vector_from_array_int(needed_v_smax_2, 5);

  psi_container_IS_startvalue[no_ps][0] = psi_container_IS_name.size();
  psi_g_p[no_ps].resize(3);
  psi_needed_g_p[no_ps].resize(3);
  int needed_g_p_0[] = {   2,    3,    4,    5,    6}; //   5
  psi_needed_g_p[no_ps][0] = get_vector_from_array_int(needed_g_p_0, 5);
  psi_container_IS_name.push_back("propagator_0_smin0_smax0");
  int needed_g_p_1[] = {   3}; //   1
  psi_needed_g_p[no_ps][1] = get_vector_from_array_int(needed_g_p_1, 1);
  psi_container_IS_name.push_back("propagator_0_smin1_smax1");
  int needed_g_p_2[] = {   5}; //   1
  psi_needed_g_p[no_ps][2] = get_vector_from_array_int(needed_g_p_2, 1);
  psi_container_IS_name.push_back("propagator_0_smin2_smax2");

  psi_container_IS_startvalue[no_ps][1] = psi_container_IS_name.size();
  psi_g_f[no_ps].resize(3);
  psi_needed_g_f[no_ps].resize(3);
  int needed_g_f_0[] = {   0,    1}; //   2
  psi_needed_g_f[no_ps][0] = get_vector_from_array_int(needed_g_f_0, 2);
  psi_container_IS_name.push_back("timelikeinvariant_smin0_smax0");
  int needed_g_f_1[] = {   0,    1,    4,    6}; //   4
  psi_needed_g_f[no_ps][1] = get_vector_from_array_int(needed_g_f_1, 4);
  psi_container_IS_name.push_back("timelikeinvariant_smin2_smax2");
  int needed_g_f_2[] = {   2}; //   1
  psi_needed_g_f[no_ps][2] = get_vector_from_array_int(needed_g_f_2, 1);
  psi_container_IS_name.push_back("timelikeinvariant_smin1_smax1");

  psi_container_IS_startvalue[no_ps][2] = psi_container_IS_name.size();
  psi_g_t[no_ps].resize(7);
  psi_needed_g_t[no_ps].resize(7);
  int needed_g_t_0[] = {   0,    1,    4,    5,    6}; //   5
  psi_needed_g_t[no_ps][0] = get_vector_from_array_int(needed_g_t_0, 5);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__1_2__16_44");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__1_2__16_44");
  int needed_g_t_1[] = {   0,    1,    4}; //   3
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
  int needed_g_t_6[] = {   6}; //   1
  psi_needed_g_t[no_ps][6] = get_vector_from_array_int(needed_g_t_6, 1);
  psi_container_IS_name.push_back("t_of_tchannel_opt_0__2_17__12_32");
  psi_container_IS_name.push_back("phi_of_tchannel_opt_0__2_17__12_32");

  psi_container_IS_startvalue[no_ps][3] = psi_container_IS_name.size();
  psi_g_d[no_ps].resize(3);
  psi_needed_g_d[no_ps].resize(3);
  int needed_g_d_0[] = {   2,    3,    4,    5,    6}; //   5
  psi_needed_g_d[no_ps][0] = get_vector_from_array_int(needed_g_d_0, 5);
  psi_container_IS_name.push_back("costheta_of_decay_12__4_8");
  psi_container_IS_name.push_back("phi_of_decay_12__4_8");
  int needed_g_d_1[] = {   3}; //   1
  psi_needed_g_d[no_ps][1] = get_vector_from_array_int(needed_g_d_1, 1);
  psi_container_IS_name.push_back("costheta_of_decay_28__16_12");
  psi_container_IS_name.push_back("phi_of_decay_28__16_12");
  int needed_g_d_2[] = {   5}; //   1
  psi_needed_g_d[no_ps][2] = get_vector_from_array_int(needed_g_d_2, 1);
  psi_container_IS_name.push_back("costheta_of_decay_44__32_12");
  psi_container_IS_name.push_back("phi_of_decay_44__32_12");

  psi_c_p[no_ps].resize(7);
  psi_c_f[no_ps].resize(7);
  psi_c_t[no_ps].resize(7);
  psi_c_d[no_ps].resize(7);

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

  int c_p_4[] = {  0}; //   1
  psi_c_p[no_ps][4] = get_vector_from_array_int(c_p_4, 1);
  int c_f_4[] = {  1}; //   1
  psi_c_f[no_ps][4] = get_vector_from_array_int(c_f_4, 1);
  int c_t_4[] = {  0,   1}; //   2
  psi_c_t[no_ps][4] = get_vector_from_array_int(c_t_4, 2);
  int c_d_4[] = {  0}; //   1
  psi_c_d[no_ps][4] = get_vector_from_array_int(c_d_4, 1);

  int c_p_5[] = {  0,   2}; //   2
  psi_c_p[no_ps][5] = get_vector_from_array_int(c_p_5, 2);
  int c_f_5[] = {}; //   0
  psi_c_f[no_ps][5] = get_vector_from_array_int(c_f_5, 0);
  int c_t_5[] = {  0}; //   1
  psi_c_t[no_ps][5] = get_vector_from_array_int(c_t_5, 1);
  int c_d_5[] = {  2,   0}; //   2
  psi_c_d[no_ps][5] = get_vector_from_array_int(c_d_5, 2);

  int c_p_6[] = {  0}; //   1
  psi_c_p[no_ps][6] = get_vector_from_array_int(c_p_6, 1);
  int c_f_6[] = {  1}; //   1
  psi_c_f[no_ps][6] = get_vector_from_array_int(c_f_6, 1);
  int c_t_6[] = {  0,   6}; //   2
  psi_c_t[no_ps][6] = get_vector_from_array_int(c_t_6, 2);
  int c_d_6[] = {  0}; //   1
  psi_c_d[no_ps][6] = get_vector_from_array_int(c_d_6, 1);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
