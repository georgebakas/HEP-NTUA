#include "../include/classes.cxx"
void munich::integration_CA_QEW(){
  static Logger logger("munich::integration_CA_QEW");
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
  perform_event_selection(osi, generic);
  osi.cut_ps[0] = 0; // to enforce the point to be evaluated !!!
  generic.calculate_dynamic_scale_TSV(0, osi);
  calculate_ME2check_CA_QEW(osi, psi);

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
    psi.calculate_initial_collinear_z1z2_IS();
    osi.z_coll = psi.z_coll;
    if (osi.switch_output_testpoint){osi.output_testpoint(psi);}
    generic.ag_psp_psp(0, 0, psi);
    psi.calculate_g_tot();
    if (munich_isnan(psi.g_tot) || munich_isinf(psi.g_tot)){errorhandling_gtot(); continue;} 

    calculate_ME2_CA_QEW(osi);

    // replace by proper nan-check !!!
    if (osi.CA_ME2_cf != osi.CA_ME2_cf){errorhandling_collinear_me2(); continue;}
    if (munich_isinf(psi.g_tot)){psi.handling_cut_psp(); continue;}

    generic.calculate_dynamic_scale(0, osi);
    generic.calculate_dynamic_scale_TSV(0, osi);
    osi.determine_scale();
    osi.determine_scale_CA();

    if (osi.switch_moment){generic.moments(osi);}

    osi.calculate_pdf_LHAPDF_CA_collinear_CV(psi.all_xz_coll_pdf);
    if (osi.switch_TSV){osi.calculate_pdf_LHAPDF_CA_collinear_TSV(psi.all_xz_coll_pdf);}
    osi.calculate_collinear_QEW();
    osi.determine_integrand_CA(psi);
    if (osi.switch_output_maxevent){osi.output_integrand_maximum_CA(psi);}
    if (munich_isnan(osi.integrand)){errorhandling_collinear_me2(); continue;}
    osi.determine_psp_weight_TSV(psi);
    osi.determine_psp_weight();
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
