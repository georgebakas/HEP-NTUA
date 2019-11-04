#include "ppttx20.header.cxx"
#include "../include/definitions.observable.set.cxx"

#define USERSWITCH(__switch__) osi_user_switch_value[osi_user_switch_map[__switch__]]
#define USERINT(__int__) osi_user_int_value[osi_user_int_map[__int__]]
#define USERDOUBLE(__double__) osi_user_double_value[osi_user_double_map[__double__]]
#define USERSTRING(__string__) osi_user_string_value[osi_user_string_map[__string__]]
#define USERCUT(__cut__) osi_user_cut_value[osi_user_cut_map[__cut__]]
#define USERPARTICLE(__particle__) user_particle[osi_user_particle_map[__particle__]]
#define PARTICLE(__particle__) osi_particle_event[osi_access_object[__particle__]][i_a]
#define NUMBER(__particle__) osi_n_object[osi_access_object[__particle__]][i_a]

void ppttx20_particles(int i_a, observable_set & oset){
  static Logger logger("ppttx20_particles");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static stringstream info_cut;
  if (osi_switch_output_cutinfo){
    info_cut.clear();
    info_cut.str("");
  }

  vector<vector<particle> > user_particle(oset.user.particle_name.size());

#include "user/specify.particles.cxx"

  for (int i_g = 0; i_g < oset.user.particle_name.size(); i_g++){
    int k_g = osi_access_object[oset.user.particle_name[i_g]];
    for (int i_p = 0; i_p < user_particle[i_g].size(); i_p++){
      if (abs(user_particle[i_g][i_p].rapidity) < osi_esi.pda[k_g].define_y && 
          abs(user_particle[i_g][i_p].eta) < osi_esi.pda[k_g].define_eta && 
          user_particle[i_g][i_p].pT > osi_esi.pda[k_g].define_pT && 
          user_particle[i_g][i_p].ET > osi_esi.pda[k_g].define_ET){
        osi_n_object[k_g][i_a]++;
        osi_particle_event[k_g][i_a].push_back(user_particle[i_g][i_p]);
      }
    }
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   " << osi_esi.pda[k_g].n_observed_min << " <= osi_n_object[" << setw(2) << k_g << "][" << setw(2) << i_a << "] = " << osi_n_object[k_g][i_a] << " <= " << osi_esi.pda[k_g].n_observed_max << endl;}
    if ((osi_n_object[k_g][i_a] < osi_esi.pda[k_g].n_observed_min) || (osi_n_object[k_g][i_a] > osi_esi.pda[k_g].n_observed_max)){
      osi_cut_ps[i_a] = -1; 
      if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   cut after user_particle[" << i_g << "] " << oset.user.particle_name[i_g] << endl;}
      logger << LOG_DEBUG_VERBOSE << "event at [i_a = " << i_a << "] cut after user_particle[" << i_g << "] " << oset.user.particle_name[i_g] << endl; 
      if (osi_switch_output_cutinfo){logger << LOG_DEBUG << endl << info_cut.str(); }
      return;
    }
    if (osi_n_object[k_g][i_a] > 1){sort(osi_particle_event[k_g][i_a].begin(), osi_particle_event[k_g][i_a].end(), greaterBypT);}
  }

  if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   event selection of user-defined particles passed and osi_cut_ps[" << i_a << "] set   (osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << ")." << endl;}
  if (osi_switch_output_cutinfo){logger << LOG_DEBUG << endl << info_cut.str();}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_cuts(int i_a, observable_set & oset){
  static Logger logger("ppttx20_cuts");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static stringstream info_cut;
  if (osi_switch_output_cutinfo){
    info_cut.clear();
    info_cut.str("");
    for (int i_o = 0; i_o < osi_esi.object_list_selection.size(); i_o++){
      int j_o = osi_equivalent_no_object[i_o]; 
      // i_o -> number in relevant_object_list
      // j_o -> number in object_list
      info_cut << "[" << setw(2) << i_a << "]   osi_esi.object_list_selection[" << setw(2) << i_o << "] = " << setw(6) << osi_esi.object_list_selection[i_o] << "   osi_n_object[" << setw(2) << j_o << "][" << setw(2) << i_a << "] = " << osi_n_object[j_o][i_a] << "   osi_particle_event[" << setw(2) << j_o << "][" << setw(2) << i_a << "].size() = " << osi_particle_event[j_o][i_a].size() << endl;
      for (int i_p = 0; i_p < osi_n_object[j_o][i_a]; i_p++){
        info_cut << "       osi_particle_event[" << setw(2) << j_o << "][" << setw(2) << i_a << "][" << setw(2) << i_p << "] = " << osi_particle_event[j_o][i_a][i_p].momentum << endl;
      }
    }
  }
#include "user/specify.cuts.cxx"

if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   event selection of user-defined cuts passed and osi_cut_ps[" << i_a << "] set   (osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << ")." << endl;}
if (osi_switch_output_cutinfo){logger << LOG_DEBUG << endl << info_cut.str();}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_cuts_test(int i_a, observable_set & oset){
  static Logger logger("ppttx20_cuts_test");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_p = 3; i_p < osi_p_parton[i_a].size(); i_p++){osi_particle_event[0][i_a][i_p] = particle(osi_p_parton[i_a][i_p].zboost(-osi_boost));}
#include "user/specify.cuts.test.cxx"

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_moments(observable_set & oset){

#include "user/specify.moments.cxx"

}

void ppttx20_calculate_dynamic_scale(int i_a, observable_set & oset){
  static Logger logger("ppttx20_calculate_dynamic_scale");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
#include "user/specify.prepare.scales.cxx"
  for (int sd = 1; sd < osi_value_mu_ren_rel.size(); sd++){
    if (osi_value_mu_ren_rel[sd].size() != 0){
      double temp_mu_central = 1.;
#include "user/specify.scales.cxx"
      for (int ss = 0; ss < osi_value_mu_ren_rel[sd].size(); ss++){
        osi_value_mu_ren[sd][ss] = temp_mu_central * osi_value_mu_ren_rel[sd][ss];
        osi_value_alpha_S[sd][ss] = LHAPDF::alphasPDF(osi_value_mu_ren[sd][ss]);
        osi_value_factor_alpha_S[sd][ss] = pow(osi_value_alpha_S[sd][ss] / osi_alpha_S, osi_contribution_order_alpha_s);
      }
    }
  }
  if (osi_id_scales == 1){
    osi_value_mu_fact = osi_value_mu_ren;
  }
  else {
    for (int sd = 1; sd < osi_value_mu_fact_rel.size(); sd++){
      if (osi_value_mu_fact_rel[sd].size() != 0){
        double temp_mu_central = 1.;
#include "user/specify.scales.cxx"
        for (int ss = 0; ss < osi_value_mu_fact_rel[sd].size(); ss++){
          osi_value_mu_fact[sd][ss] = temp_mu_central * osi_value_mu_fact_rel[sd][ss];
        }
      }
    }
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_calculate_dynamic_scale_RA(int i_a, observable_set & oset){
  static Logger logger("ppttx20_calculate_dynamic_scale_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
#include "user/specify.prepare.scales.cxx"
  for (int sd = 1; sd < osi_value_mu_ren_rel.size(); sd++){
    if (osi_value_mu_ren_rel[sd].size() != 0){
      double temp_mu_central = 1.;
#include "user/specify.scales.cxx"
      for (int ss = 0; ss < osi_value_mu_ren_rel[sd].size(); ss++){
        osi_RA_value_mu_ren[i_a][sd][ss] = temp_mu_central * osi_value_mu_ren_rel[sd][ss];
        osi_RA_value_alpha_S[i_a][sd][ss] = LHAPDF::alphasPDF(osi_RA_value_mu_ren[i_a][sd][ss]);
        osi_RA_value_factor_alpha_S[i_a][sd][ss] = pow(osi_RA_value_alpha_S[i_a][sd][ss] / osi_alpha_S, osi_contribution_order_alpha_s);
      }
    }
  }
  if (osi_id_scales == 1){
    osi_RA_value_mu_fact[i_a] = osi_RA_value_mu_ren[i_a];
  }
  else {
    for (int sd = 1; sd < osi_value_mu_fact_rel.size(); sd++){
      if (osi_value_mu_fact_rel[sd].size() != 0){
        double temp_mu_central = 1.;
#include "user/specify.scales.cxx"
        for (int ss = 0; ss < osi_value_mu_fact_rel[sd].size(); ss++){
          osi_RA_value_mu_fact[i_a][sd][ss] = temp_mu_central * osi_value_mu_fact_rel[sd][ss];
        }
      }
    }
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_calculate_dynamic_scale_TSV(int i_a, observable_set & oset){
  static Logger logger("ppttx20_calculate_dynamic_scale_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
#include "user/specify.prepare.scales.cxx"
  for (int sd = 1; sd < osi_max_dyn_ren + 1; sd++){
    double temp_mu_central = 1.;
#include "user/specify.scales.cxx"
    for (int ss = 0; ss < osi_n_scale_dyn_ren[sd]; ss++){
      osi_value_scale_ren[i_a][sd][ss] = temp_mu_central * osi_value_relative_scale_ren[sd][ss];
      if (osi_needed_scale2_ren){osi_value_scale2_ren[i_a][sd][ss] = pow(osi_value_scale_ren[i_a][sd][ss], 2);}
      osi_value_alpha_S_TSV[i_a][sd][ss] = LHAPDF::alphasPDF(osi_value_scale_ren[i_a][sd][ss]);
      osi_value_relative_factor_alpha_S[i_a][sd][ss] = pow(osi_value_alpha_S_TSV[i_a][sd][ss] / osi_alpha_S, osi_contribution_order_alpha_s);
    }
  }
  for (int sd = 1; sd < osi_max_dyn_fact + 1; sd++){
    double temp_mu_central = 1.;
#include "user/specify.scales.cxx"
    osi_value_central_scale_fact[sd] = temp_mu_central;
    if (osi_needed_scale2_fact){osi_value_central_logscale2_fact[sd] = 2 * log(temp_mu_central);}
    for (int ss = 0; ss < osi_n_scale_dyn_fact[sd]; ss++){
      osi_value_scale_fact[i_a][sd][ss] = osi_value_central_scale_fact[sd] * osi_value_relative_scale_fact[sd][ss];
      if (osi_needed_scale2_fact){osi_value_scale2_fact[i_a][sd][ss] = pow(osi_value_scale_fact[i_a][sd][ss], 2);}
    }
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void ppttx20_phasespacepoint_born(observable_set & oset){
  static Logger logger("ppttx20_phasespacepoint_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

#include "user/testpoint.born.cxx"

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
} 

void ppttx20_phasespacepoint_collinear(observable_set & oset){
  static Logger logger("ppttx20_phasespacepoint_collinear");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

#include "user/testpoint.collinear.cxx"

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
} 

void ppttx20_phasespacepoint_real(observable_set & oset){
  static Logger logger("ppttx20_phasespacepoint_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

#include "user/testpoint.real.cxx"

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
} 

void ppttx20_phasespacepoint_realcollinear(observable_set & oset){
  static Logger logger("ppttx20_phasespacepoint_realcollinear");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

#include "user/testpoint.realcollinear.cxx"

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
} 

void ppttx20_phasespacepoint_doublereal(observable_set & oset){
  static Logger logger("ppttx20_phasespacepoint_doublereal");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

#include "user/testpoint.doublereal.cxx"

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
} 
