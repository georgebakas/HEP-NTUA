#include "../include/classes.cxx"
void munich::integration_VA_MIX(){
  static Logger logger("munich::integration_VA_MIX");
  logger << LOG_INFO << "called" << endl;

  generic.determination_subprocess_psp(0, psi);
  generic.combination_subprocess_psp(0, psi, osi);

  osi.initialization_complete(isi, psi);
  
  psi.initialization_complete(isi);

  osi.initialization_integration(psi);

  observable_set save_osi = osi;
  phasespace_set save_psi = psi;
  runresumption_set rsi(osi, psi, save_osi, save_psi);

  generic.phasespacepoint_psp(osi);
  calculate_ME2check_VA_MIX(osi);

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
    if (osi.cut_ps[0] == -1){psi.handling_cut_psp(); continue;}

    if (osi.switch_console_output_tau_0){psi.output_check_tau_0();}
    generic.ag_psp_psp(0, 0, psi);
    psi.calculate_g_tot();
    if (munich_isnan(psi.g_tot) || munich_isinf(psi.g_tot)){errorhandling_gtot(); continue;} // isinf case ??? 

    generic.calculate_dynamic_scale(0, osi);
    generic.calculate_dynamic_scale_TSV(0, osi);
    osi.determine_scale();

    calculate_ME2_VA_MIX(osi);

    if (munich_isnan(osi.VA_V_ME2) || munich_isnan(osi.VA_X_ME2) || munich_isnan(osi.VA_I_ME2)){errorhandling_me2(); continue;}
    if (osi.VA_V_ME2 == 0. && osi.switch_VI != 2){errorhandling_OL(); continue;}

    if (osi.switch_moment){generic.moments(osi);}

    osi.calculate_pdf_LHAPDF_CV();
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
    osi.determine_runtime(psi);
  }
  osi.output_finalization_integration(psi);
}
