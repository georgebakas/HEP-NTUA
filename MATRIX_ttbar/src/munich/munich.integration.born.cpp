#include "../include/classes.cxx"
void munich::integration_born(){
  static Logger logger("munich::integration_born");
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
  osi.calculate_ME2check_born();

  osi.output_check_running_alpha_S();
  //  integrate_topwidth(osi, psi);

  //  if (psi.n_events_max == 0){osi.int_end = 1;}
  if (osi.user.int_value[osi.user.int_map["rescaling_exponent"]] != 0){psi.hcf = psi.hcf * pow(osi.user.double_value[osi.user.double_map["rescaling_factor"]], osi.user.int_value[osi.user.int_map["rescaling_exponent"]]);}

  if (osi.int_end){osi.output_zero_contribution_complete(psi);}
  osi.initialization_runtime();
  while (osi.int_end == 0){
    rsi.perform_iteration_step();
    if (osi.int_end == 1){break;}
    //    if (*psi.i_step_mode > 0 && *psi.i_step_mode % psi.n_step == 0){rsi.perform_iteration_step();}
    psi.calculate_IS();
    generic.ac_psp_psp(0, psi.MC_phasespace.channel, psi);
    osi.determine_p_parton(psi);
    if (osi.p_parton != osi.p_parton){errorhandling_c_psp(); continue;}
    osi.cut_ps[0] = 0;
    perform_event_selection(osi, generic);
    if (osi.cut_ps[0] == -1){psi.handling_cut_psp(); continue;}
    if (osi.switch_console_output_tau_0){psi.output_check_tau_0();}
    if (osi.switch_output_testpoint){osi.output_testpoint(psi);}
    generic.ag_psp_psp(0, 0, psi);
    psi.calculate_g_tot();
    if (munich_isnan(psi.g_tot) || munich_isinf(psi.g_tot)){errorhandling_gtot(); continue;} // isinf case ???
    ///**///
    osi.calculate_ME2_born();
    
    if (munich_isnan(osi.ME2)){errorhandling_me2(); continue;}
    // new:
    if (csi.type_contribution == "L2RT" ||
	csi.type_contribution == "L2RJ"){
      osi.counter_acc_qTcut[osi.cut_ps[0]]++;
      if (osi.ME2 == 0.){
	osi.counter_killed_qTcut[osi.cut_ps[0]]++;
	psi.i_tec++;
	psi.handling_techcut_psp();
	continue;
      }
    }
    /*
    // old: would typically impose a bias for all qTcut values !!!
    if (csi.type_contribution == "L2RT"){osi.counter_acc_qTcut[osi.cut_ps[0]]++;}
    if (csi.type_contribution == "L2RT" && osi.ME2 == 0.){errorhandling_OL(); continue;}
    */
    if (!(osi.check_vanishing_ME2_end) && csi.type_contribution == "born"){handling_vanishing_me2(); if (osi.flag_vanishing_ME2){continue;}}

    generic.calculate_dynamic_scale(0, osi);
    generic.calculate_dynamic_scale_TSV(0, osi);
    osi.determine_scale();

    if (osi.switch_moment){generic.moments(osi);}

    ///**///
    osi.calculate_pdf_LHAPDF_CV();
    ///**///
    osi.calculate_pdf_LHAPDF_TSV();

    if (munich_isnan(osi.pdf_factor[0])){errorhandling_pdf(); continue;}
    if (munich_isnan(osi.var_rel_alpha_S)){errorhandling_alpha_S(); continue;}

    osi.determine_integrand(psi);
    osi.determine_psp_weight_TSV(psi);
    osi.determine_psp_weight();
    if (osi.switch_output_maxevent){osi.output_integrand_maximum(psi);}
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
