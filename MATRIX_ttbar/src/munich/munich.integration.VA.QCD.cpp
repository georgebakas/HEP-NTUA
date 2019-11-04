#include "../include/classes.cxx"
void munich::integration_VA_QCD(){
  static Logger logger("munich::integration_VA_QCD");
  logger << LOG_INFO << "called" << endl;

  generic.determination_subprocess_psp(0, psi);
  generic.combination_subprocess_psp(0, psi, osi);

  osi.initialization_complete(isi, psi);

  psi.initialization_complete(isi);

  osi.initialization_integration(psi);

  observable_set save_osi = osi;
  phasespace_set save_psi = psi;
  runresumption_set rsi(osi, psi, save_osi, save_psi);

  calculate_ME2check_VA_QCD(osi);
 
  if (psi.n_events_max == 0){osi.int_end = 1;}
  if (osi.user.int_value[osi.user.int_map["rescaling_exponent"]] != 0){psi.hcf = psi.hcf * pow(osi.user.double_value[osi.user.double_map["rescaling_factor"]], osi.user.int_value[osi.user.int_map["rescaling_exponent"]]);}

  if (osi.int_end){osi.output_zero_contribution_complete(psi);}
  osi.initialization_runtime();
  while (osi.int_end == 0){
    rsi.perform_iteration_step();
    if (osi.int_end == 1){break;}

    psi.calculate_IS();
    generic.ac_psp_psp(0, psi.MC_phasespace.channel, psi);
    osi.determine_p_parton(psi);
    // replace by proper nan-check !!!
    if (osi.p_parton[0] != osi.p_parton[0]){errorhandling_c_psp(); continue;}
    perform_event_selection(osi, generic);
    logger << LOG_DEBUG_VERBOSE << "psi.i_gen = " << psi.i_gen << "   psi.i_acc = " << psi.i_acc << "   psi.i_rej = " << psi.i_rej << "   osi.cut_ps[0] = " << osi.cut_ps[0] << endl;
    if (osi.cut_ps[0] == -1){psi.handling_cut_psp(); continue;}

    if (osi.switch_console_output_tau_0){psi.output_check_tau_0();}
    generic.ag_psp_psp(0, 0, psi);
    psi.calculate_g_tot();
    if (munich_isnan(psi.g_tot) || munich_isinf(psi.g_tot)){errorhandling_gtot(); continue;} // isinf case ???

    generic.calculate_dynamic_scale(0, osi);
    generic.calculate_dynamic_scale_TSV(0, osi);
    osi.determine_scale();

    // to avoid issues with COLLIER in phase-space divergent integration
    /*
    bool temp_cut_technical = false;
    for (int sr = 0; sr < psi.RA_singular_region_list.size(); sr++){
      int x1 = psi.RA_singular_region_list[sr][0];
      int x2 = psi.RA_singular_region_list[sr][1];
      psi.RA_singular_region[x1][x2] = osi.p_parton[0][x1] * osi.p_parton[0][x2] / psi.xbs_all[0][0];
      logger << LOG_DEBUG_VERBOSE << "psi.i_gen = " << psi.i_gen << "   psi.i_acc = " << psi.i_acc << "   psi.i_rej = " << psi.i_rej << "    psi.RA_singular_region[" << x1 << "][" << x2 << "] = " << psi.RA_singular_region[x1][x2] << "   psi.cut_technical = " << psi.cut_technical << endl;
      if (psi.RA_singular_region[x1][x2] < psi.cut_technical){temp_cut_technical = true;}
    }
    if (temp_cut_technical){errorhandling_cut_technical(); continue;}
    */

    ///**///
    calculate_ME2_VA_QCD(osi);
    ///**///    osi.VA_V_ME2 = 1.;
    
    if (munich_isnan(osi.VA_V_ME2) || munich_isnan(osi.VA_X_ME2) || munich_isnan(osi.VA_I_ME2)){errorhandling_me2(); continue;}
    if (osi.type_contribution == "RVA"){osi.counter_acc_qTcut[osi.cut_ps[0]]++;}
    if ((osi.type_contribution == "VA" || osi.type_contribution == "RVA") && osi.VA_V_ME2 == 0. && osi.VA_X_ME2 == 0. && osi.switch_VI != 2){errorhandling_OL(); continue;} // maybe better error output in case of RVA contribution !!!

    if (osi.switch_moment){generic.moments(osi);}

    ///**///
    osi.calculate_pdf_LHAPDF_CV();
    ///**///
    osi.calculate_pdf_LHAPDF_TSV();

    osi.determine_integrand_VA(psi);
    osi.determine_psp_weight_TSV(psi);
    osi.determine_psp_weight();
    if (osi.switch_output_maxevent){osi.output_integrand_maximum_VA(psi);}
    if (osi.switch_distribution){osi.determine_distribution_complete();}
    static double optimization_modifier = 1.;
    if (osi.user.switch_value[osi.user.switch_map["optimization_modifier"]]){
      optimization_modifier = pow(osi.particle_event[osi.access_object["Vrec"]][0][0].pT, osi.user.switch_value[osi.user.switch_map["optimization_modifier"]]);
      logger << LOG_DEBUG << "optimization_modifier = " << optimization_modifier << endl;
    }   

    psi.psp_MCweight_optimization(osi.integrand, osi.this_psp_weight, osi.this_psp_weight2, optimization_modifier);
    psi.i_acc++;
    logger << LOG_DEBUG_VERBOSE << "psi.i_gen = " << psi.i_gen << "   psi.i_acc = " << psi.i_acc << "   psi.i_rej = " << psi.i_rej << "   osi.VA_V_ME2 = " << osi.VA_V_ME2 << endl;
    osi.determine_runtime(psi);
  }
  osi.output_finalization_integration(psi);
}
