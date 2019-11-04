#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../ppttx20/header/ppttx20.phasespace.h"
void ppttx20_ax_psp_200_ddx_ttx(int no_ps, phasespace_set & psi){
  static Logger logger("ppttx20_ax_psp_200_ddx_ttx");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  psi_v_smin[no_ps].resize(0);
  psi_needed_v_smin[no_ps].resize(0);

  psi_v_smax[no_ps].resize(0);
  psi_needed_v_smax[no_ps].resize(0);

  psi_container_IS_startvalue[no_ps][0] = psi_container_IS_name.size();
  psi_g_p[no_ps].resize(0);
  psi_needed_g_p[no_ps].resize(0);

  psi_container_IS_startvalue[no_ps][1] = psi_container_IS_name.size();
  psi_g_f[no_ps].resize(0);
  psi_needed_g_f[no_ps].resize(0);

  psi_container_IS_startvalue[no_ps][2] = psi_container_IS_name.size();
  psi_g_t[no_ps].resize(0);
  psi_needed_g_t[no_ps].resize(0);

  psi_container_IS_startvalue[no_ps][3] = psi_container_IS_name.size();
  psi_g_d[no_ps].resize(1);
  psi_needed_g_d[no_ps].resize(1);
  int needed_g_d_0[] = {   0}; //   1
  psi_needed_g_d[no_ps][0] = get_vector_from_array_int(needed_g_d_0, 1);
  psi_container_IS_name.push_back("costheta_of_decay_12__4_8");
  psi_container_IS_name.push_back("phi_of_decay_12__4_8");

  psi_c_p[no_ps].resize(1);
  psi_c_f[no_ps].resize(1);
  psi_c_t[no_ps].resize(1);
  psi_c_d[no_ps].resize(1);

  int c_p_0[] = {}; //   0
  psi_c_p[no_ps][0] = get_vector_from_array_int(c_p_0, 0);
  int c_f_0[] = {}; //   0
  psi_c_f[no_ps][0] = get_vector_from_array_int(c_f_0, 0);
  int c_t_0[] = {}; //   0
  psi_c_t[no_ps][0] = get_vector_from_array_int(c_t_0, 0);
  int c_d_0[] = {  0}; //   1
  psi_c_d[no_ps][0] = get_vector_from_array_int(c_d_0, 1);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
