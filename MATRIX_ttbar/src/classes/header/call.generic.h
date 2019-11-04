#ifndef CALLGENERIC_H
#define CALLGENERIC_H

#include <string>
#include <vector>

#include "randommanager.h"
#include "logger.h"
#include "../../include/ftypes.h"
#include "fourvector.h"
#include "particle.h"
//#include "phasespace.set.h"
#include "dddistribution.h"
#include "xdistribution.h"

using namespace std;

class call_generic{
 private:

 public:
////////////////////
//  constructors  //
////////////////////
  call_generic();

  void (*calculate_dynamic_scale)(int i_a, observable_set & oset);
  void (*calculate_dynamic_scale_RA)(int i_a, observable_set & oset);
  void (*calculate_dynamic_scale_TSV)(int i_a, observable_set & oset);
  void (*moments)(observable_set & oset);
  void (*particles)(int i_a, observable_set & oset);
  void (*cuts)(int i_a, observable_set & oset);
  void (*cuts_test)(int i_a, observable_set & oset);
  void (*determination_subprocess_psp)(int i_a, phasespace_set & psi);
  //  void (*determination_subprocess_psp)(int i_a, vector<int> & tp, phasespace_set & psi);
  void (*combination_subprocess_psp)(int i_a, phasespace_set & psi, observable_set & oset);
  void (*optimize_minv_psp)(phasespace_set & psi);
  int (*determination_MCchannels_psp)(int no_ps, phasespace_set & psi);
  void (*ac_tau_psp_psp)(int no_ps, vector<int> & tau_MC_map, phasespace_set & psi);
  void (*ax_psp_psp)(int no_ps, phasespace_set & psi);
  void (*ac_psp_psp)(int no_ps, int MC_channel, phasespace_set & psi);
  void (*ag_psp_psp)(int no_ps, int zero, phasespace_set & psi);
  void (*phasespacepoint_psp)(observable_set & oset);

  // only required for RA contributions
  void (*determination_no_subprocess_dipole)(int & no_map, vector<int> & o_map, int & no_prc, vector<int> & o_prc, double & factor_symmetry, vector<int> & tp, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference);
  void (*determination_subprocess_dipole)(int i_a, phasespace_set & psi);
  //  void (*determination_subprocess_dipole)(int i_a, vector<int> & tp, phasespace_set & psi);

  int (*determination_MCchannels_dipole)(int no_ps, phasespace_set & psi);
  void (*ac_tau_psp_dipole)(int no_ps, vector<int> & tau_MC_map, phasespace_set & psi);
  void (*ax_psp_dipole)(int no_ps, phasespace_set & psi);
  void (*ac_psp_dipole)(int no_ps, int MC_channel, phasespace_set & psi);
  void (*ag_psp_dipole)(int no_ps, int zero, phasespace_set & psi);

  void (*determination_no_subprocess_doubledipole)(int & no_map, vector<int> & o_map, int & no_prc, vector<int> & o_prc, double & factor_symmetry, vector<int> & tp, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference);
  void (*determination_subprocess_doubledipole)(int i_a, phasespace_set & psi);
  //  void (*determination_subprocess_doubledipole)(int i_a, vector<int> & tp, phasespace_set & psi);
  int (*determination_MCchannels_doubledipole)(int no_ps, phasespace_set & psi);
  void (*ac_tau_psp_doubledipole)(int no_ps, vector<int> & tau_MC_map, phasespace_set & psi);
  void (*ax_psp_doubledipole)(int no_ps, phasespace_set & psi);
  void (*ac_psp_doubledipole)(int no_ps, int MC_channel, phasespace_set & psi);
  void (*ag_psp_doubledipole)(int no_ps, int zero, phasespace_set & psi);
  
  /*
  // could be combined at some point...
  void (*QCD_calculate_collinear)(vector<vector<collinear_set> > & CS_collinear, vector<double> & x_pdf, vector<double> & z_coll, vector<vector<double> > & ME2_cf, vector<vector<double> > & CA_value_log_mu2_fact, vector<vector<vector<vector<double> > > > & CA_value_ME2_KP, vector<int> & type_parton, vector<vector<int> > & dipole_splitting, parameter_set & pset, observable_set & oset);
  void (*QEW_calculate_collinear)(vector<vector<collinear_set> > & CS_collinear, vector<double> & x_pdf, vector<double> & z_coll, vector<vector<double> > & ME2_cf, vector<vector<double> > & CA_value_log_mu2_fact, vector<vector<vector<vector<double> > > > & CA_value_ME2_KP, vector<int> & type_parton, vector<vector<int> > & dipole_splitting, parameter_set & pset, observable_set & oset);
  */

  void (*calculate_H2)(observable_set & oset);
  void (*calculate_H1gg)(observable_set & oset);

  void (*list_subprocess_born)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);

  // CS/CDST subtraction at NLO
  void (*list_subprocess_V_QCD)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_C_QCD)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_R_QCD)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);

  // CS/CDST subtraction - extension to QEW corrections at NLO
  void (*list_subprocess_V_QEW)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_C_QEW)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_R_QEW)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);

  void (*list_subprocess_V_MIX)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_R_MIX)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);

  // QT subtraction at NNLO QCD
  void (*list_subprocess_loop)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);

  void (*list_subprocess_V2_QCD)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_C2_QCD)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_RV_QCD)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_RC_QCD)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_RR_QCD)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);

  // QT subtraction - extension to mixed QCD-EW in NNLO
  void (*list_subprocess_V2_QEW)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_C2_QEW)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_RV_QEW)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_RC_QEW)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_RR_QEW)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);

  // QT subtraction - extension to EW in NNLO
  void (*list_subprocess_V2_MIX)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_C2_MIX)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_RV_MIX)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
  void (*list_subprocess_RR_MIX)(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
};

#endif
