#include "../include/classes.cxx"
void munich::integration_RA_MIX(){
  static Logger logger("munich::integration_RA_MIX");
  logger << LOG_INFO << "called" << endl;

  generic.determination_subprocess_psp(0, psi);
  generic.combination_subprocess_psp(0, psi, osi);

  // could be modified later ??? merge with MC-mappings of dipoles ???
  psi.MC_n_channel_phasespace[0] = generic.determination_MCchannels_psp(0, psi);
  dipole.push_back(dipole_set("ME2_R", csi.type_parton[0], csi.basic_type_parton[0], psi.symmetry_factor, psi.no_map[0], psi.o_map[0], psi.no_prc[0], psi.o_prc[0], psi.MC_n_channel_phasespace[0]));
  QCD_determine_dipoles(dipole_candidate, csi.type_parton[0], csi.basic_type_parton[0]);
  QCD_selection_dipoles(dipole, dipole_candidate, psi.phasespace_order_alpha_s[0] - 1, psi.phasespace_order_alpha_e[0], psi.phasespace_order_interference[0], psi.RA_singular_region, psi.RA_singular_region_name, psi.RA_singular_region_list, psi, generic);
  QEW_determine_dipoles(dipole_candidate, csi.type_parton[0], csi.basic_type_parton[0]);
  QEW_selection_dipoles(dipole, dipole_candidate, psi.phasespace_order_alpha_s[0], psi.phasespace_order_alpha_e[0] - 1, psi.phasespace_order_interference[0], psi.RA_singular_region, psi.RA_singular_region_name, psi.RA_singular_region_list, psi, generic);

  // Should be replaced !!!
  static int n_dipoles = dipole.size();
  osi.initialization_complete(isi, psi);

  psi.initialization_complete_RA(isi, dipole);

  osi.initialization_integration(psi);

  observable_set save_osi = osi;
  phasespace_set save_psi = psi;
  runresumption_set rsi(osi, psi, save_osi, save_psi);

  calculate_ME2check_RA_MIX(psi, osi);
  osi.max_integrand = osi.sigma_normalization; // ???

  if (psi.n_events_max == 0){osi.int_end = 1;}
  if (osi.user.int_value[osi.user.int_map["rescaling_exponent"]] != 0){psi.hcf = psi.hcf * pow(osi.user.double_value[osi.user.double_map["rescaling_factor"]], osi.user.int_value[osi.user.int_map["rescaling_exponent"]]);}

  osi.initialization_runtime();
  if (osi.int_end){osi.output_zero_contribution_complete(psi);}
  while (osi.int_end == 0){
    rsi.perform_iteration_step();
    if (osi.int_end == 1){break;}

    psi.calculate_IS();

    osi.determine_p_parton(psi);
    // replace by proper nan-check !!!
    if (psi.xbp_all[0] != psi.xbp_all[0]){errorhandling_c_psp_initial(); continue;}
    psi.calculate_IS_RA(dipole);
    if (psi.RA_x_a == 0){generic.ac_psp_psp(0, psi.MC_phasespace.channel, psi);}
    else {  //  if !(osi.switch_RS_mapping){}
      if      (dipole[psi.RA_x_a].type_dipole() == 1 && dipole[psi.RA_x_a].massive() == 1){psi.ac_psp_RA_group_ij_k_massive(dipole, psi.RA_x_a);}
      else if (dipole[psi.RA_x_a].type_dipole() == 2 && dipole[psi.RA_x_a].massive() == 1){psi.ac_psp_RA_group_ij_a_massive(dipole, psi.dipole_sinx_min[psi.RA_x_a], psi.RA_x_a);}
      else if (dipole[psi.RA_x_a].type_dipole() == 3 && dipole[psi.RA_x_a].massive() == 1){psi.ac_psp_RA_group_ai_k_massive(dipole, psi.dipole_sinx_min[psi.RA_x_a], psi.RA_x_a);}
      else if (dipole[psi.RA_x_a].type_dipole() == 1 && dipole[psi.RA_x_a].massive() == 0){psi.ac_psp_RA_group_ij_k(dipole, psi.RA_x_a);}
      else if (dipole[psi.RA_x_a].type_dipole() == 2 && dipole[psi.RA_x_a].massive() == 0){psi.ac_psp_RA_group_ij_a(dipole, psi.dipole_sinx_min[psi.RA_x_a], psi.RA_x_a);}
      else if (dipole[psi.RA_x_a].type_dipole() == 3 && dipole[psi.RA_x_a].massive() == 0){psi.ac_psp_RA_group_ai_k(dipole, psi.dipole_sinx_min[psi.RA_x_a], psi.RA_x_a);}
      else if (dipole[psi.RA_x_a].type_dipole() == 5){psi.ac_psp_RA_group_ai_b(dipole, psi.dipole_sinx_min[psi.RA_x_a], psi.RA_x_a);}
    }
    psi.correct_phasespacepoint_real();
    for (int i_p = 0; i_p < osi.p_parton[0].size(); i_p++){osi.p_parton[0][psi.o_map[0][i_p]] = psi.xbp_all[0][intpow(2, i_p - 1)];}
    // replace by proper nan-check !!!
    if (psi.xbp_all[0] != psi.xbp_all[0]){errorhandling_c_xbpsp(); continue;}
    for (int ib = 1; ib < osi.p_parton[0].size(); ib++){logger << LOG_DEBUG_VERBOSE << "osi.p_parton[0][" << ib << "][0] = " << psi.xbp_all[0][intpow(2, ib - 1)] << "   " << sqrt(abs(psi.xbp_all[0][intpow(2, ib - 1)].m2())) << endl;}
    for (int ib = 0; ib < osi.p_parton[0].size(); ib++){logger << LOG_DEBUG_VERBOSE << "osi.p_parton[0][" << ib << "] = " << osi.p_parton[0][ib] << "   " << osi.p_parton[0][ib].m2() << "   " << sqrt(abs(osi.p_parton[0][ib].m2())) << endl;}

    psi.determine_dipole_phasespace_RA(dipole);
    for (int i_a = 1; i_a < n_dipoles; i_a++){for (int i_p = 0; i_p < osi.p_parton[i_a].size(); i_p++){osi.p_parton[i_a][i_p] = psi.xbp_all[i_a][intpow(2, i_p - 1)];}}
    psi.correct_phasespacepoint_dipole(); // Which phase-space is acutally corrected here (psi -- osi) ???

    // replace by proper nan-check !!!
    if (osi.p_parton[0] != osi.p_parton[0]){errorhandling_c_psp(); continue;}
    perform_event_selection(osi, generic);
    if (osi.switch_RS == 1){for (int i_a = 1; i_a < n_dipoles; i_a++){osi.cut_ps[i_a] = -1;}}
    else if (osi.switch_RS == 2){osi.cut_ps[0] = -1;}
    osi.first_non_cut_ps = -1;
    for (int i_a = 0; i_a < osi.cut_ps.size(); i_a++){if (osi.cut_ps[i_a] > -1){osi.first_non_cut_ps = i_a; break;}}
    if (osi.first_non_cut_ps == -1){psi.handling_cut_psp(); continue;}

    if (osi.switch_console_output_tau_0){psi.output_check_tau_0();}
    if (osi.switch_output_testpoint){osi.output_testpoint(psi);}
    ///    if (psi.i_acc == psi.n_events_max){osi.int_end = 1;}
    generic.ag_psp_psp(0, 0, psi);
    if (psi.switch_RS_mapping == 0){
      //    if (osi.switch_RS_mapping == 0){
      for (int i_a = 1; i_a < n_dipoles; i_a++){
        if      (dipole[i_a].type_dipole() == 1 && dipole[i_a].massive() == 1){psi.ag_psp_RA_group_ij_k_massive(dipole, i_a);}
        else if (dipole[i_a].type_dipole() == 2 && dipole[i_a].massive() == 1){psi.ag_psp_RA_group_ij_a_massive(psi.dipole_sinx_min[i_a], dipole, i_a);}
        else if (dipole[i_a].type_dipole() == 3 && dipole[i_a].massive() == 1){psi.ag_psp_RA_group_ai_k_massive(psi.dipole_sinx_min[i_a], dipole, i_a);}
        else if (dipole[i_a].type_dipole() == 1 && dipole[i_a].massive() == 0){psi.ag_psp_RA_group_ij_k(dipole, i_a);}
        else if (dipole[i_a].type_dipole() == 2 && dipole[i_a].massive() == 0){psi.ag_psp_RA_group_ij_a(psi.dipole_sinx_min[i_a], dipole, i_a);}
        else if (dipole[i_a].type_dipole() == 3 && dipole[i_a].massive() == 0){psi.ag_psp_RA_group_ai_k(psi.dipole_sinx_min[i_a], dipole, i_a);}
        else if (dipole[i_a].type_dipole() == 5){psi.ag_psp_RA_group_ai_b(psi.dipole_sinx_min[i_a], dipole, i_a);}
      }
    }
    psi.calculate_g_tot();
    if (munich_isnan(psi.g_tot)){errorhandling_RA_gtot(); continue;} 

    calculate_ME2_RA_MIX(osi);

    double temp_sum_RA_ME2 = accumulate(osi.RA_ME2.begin(), osi.RA_ME2.end(), 0.);
    if (munich_isinf(temp_sum_RA_ME2)){errorhandling_RA_me2(3); continue;}
    if (munich_isnan(temp_sum_RA_ME2)){errorhandling_RA_me2(0); continue;}
    if (!(osi.check_vanishing_ME2_end)){handling_vanishing_me2(); if (osi.flag_vanishing_ME2){continue;}}

    if (psi.g_tot == 0.){logger << LOG_INFO << "psi.g_tot == 0." << endl; errorhandling_RA_me2(1); continue;} // does that really happen?
    if (munich_isinf(psi.g_tot)){psi.handling_cut_psp(); continue;}

    for (int i_a = 0; i_a < n_dipoles; i_a++){
      if (osi.cut_ps[i_a] != -1){
        generic.calculate_dynamic_scale_RA(i_a, osi);
        generic.calculate_dynamic_scale_TSV(i_a, osi);
	osi.determine_scale_RA(i_a);
      }
    }

    // add proper errorhandling, also for mu_fact !!!
    // replace by proper nan-check !!!
    if (osi.switch_CV && osi.RA_value_factor_alpha_S != osi.RA_value_factor_alpha_S){errorhandling_alpha_S(); continue;}
    if (osi.switch_TSV && osi.value_relative_factor_alpha_S != osi.value_relative_factor_alpha_S){errorhandling_alpha_S(); continue;}

    if (osi.switch_moment){generic.moments(osi);}

    osi.calculate_pdf_LHAPDF_RA_CV();
    //    calculate_pdf_LHAPDF_RA_CV(osi);
    osi.calculate_pdf_LHAPDF_TSV();
    //    calculate_pdf_LHAPDF_TSV(osi);

    osi.determine_integrand_RA(psi);

    // introduce proper nan-check !!!
    if (osi.var_RA_ME2 != osi.var_RA_ME2){errorhandling_RA_me2(1); continue;}
    if (osi.switch_CV && osi.var_RA_ME2_CV != osi.var_RA_ME2_CV){errorhandling_RA_me2(1); continue;}
    osi.determine_techcut_RA(psi);
    if (psi.RA_techcut){psi.RA_techcut = 0; psi.handling_cut_psp(); continue;}

    // introduce proper nan-check !!!
    if (osi.var_RA_ME2_CV != osi.var_RA_ME2_CV){errorhandling_RA_me2(2); continue;}
    if (osi.switch_output_cancellation_check){osi.output_integrand_cancellation_check_RA(psi);}
    if (osi.switch_output_maxevent){osi.output_integrand_maximum_RA(psi);}
    //    logger << LOG_INFO << setw(8) << psi.i_acc << "   RA_ME2 = " << setprecision(15) << setw(23) << abs(accumulate(osi.RA_ME2.begin(), osi.RA_ME2.end(), 0.)) << "   g_tot = " << setprecision(15) << psi.g_tot << endl;
    osi.determine_psp_weight_TSV(psi);
    osi.determine_psp_weight_RA(psi);

    if (osi.switch_distribution){osi.determine_distribution_complete();}
    static double optimization_modifier = 1.;
    if (osi.user.switch_value[osi.user.switch_map["optimization_modifier"]]){
      //      optimization_modifier = pow(osi.particle_event[osi.access_object["Vrec"]][0][0].pT, osi.user.switch_value[osi.user.switch_map["optimization_modifier"]]);
      optimization_modifier = pow(psi.x_pdf[0], osi.user.switch_value[osi.user.switch_map["optimization_modifier"]]);
      logger << LOG_DEBUG << "optimization_modifier = " << optimization_modifier << endl;
    }   

    psi.psp_MCweight_optimization(osi.integrand, osi.this_psp_weight, osi.this_psp_weight2, optimization_modifier);
    psi.i_acc++;
    osi.determine_runtime(psi);
  }
  osi.output_finalization_integration(psi);
}
