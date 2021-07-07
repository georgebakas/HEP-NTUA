#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../include/definitions.observable.set.cxx"
#include "../qTsubtraction/header/more.h"
void munich::integration_VJ_QCD(){
  static Logger logger("munich::integration_VJ_QCD");
  logger << LOG_INFO << "called" << endl;

  psi.initialization_phasespace_born();

  generic.determination_subprocess_psp(0, psi);
  generic.combination_subprocess_psp(0, psi, osi);

  generic.ax_psp_psp(0, psi);
  psi_MC_n_channel = generic.determination_MCchannels_psp(0, psi);
  psi.initialization_mapping_parameter(isi);
  psi.initialization_minimum_tau(isi);
  psi.initialization_minimum_phasespacecut(isi);
  /*
  psi.initialization_optimization(isi);
  ///  psi.initialization_filename(isi);

  // generation of optimized xbs limits (including invariant-mass cuts)
  psi.initialization_phasespace_subprocess();
  generic.optimize_minv_psp(psi);
  psi.initialization_phasespace_subprocess_optimization();
  */
    
  psi.initialization_complete(isi);

  // generation of MC_tau mappings
  //  generic.ac_tau_psp_psp(0, psi_tau_MC_map, psi);
  ///  psi.initialization_phasespace_MC_tau();

  ///  psi.initialization_phasespace_IS();
  psi.initialization_MC();

  // check position of the following phasespace_IS calls !!!
  osi.determine_CX_QCD(psi, 1);

  psi.initialization_phasespace_IS_CX();

  osi.initialization_NJ();
  osi.initialization_TSV();
  osi.initialization_CV();
  osi.initialization_CX_ncollinear();

  osi.initialization_runtime_partonlevel();
  osi.initialization_distribution();
  osi.initialization_integration(psi);

  osi.initialization_QT_coefficients(); // needed for NJ subtraction ??? --> so far I have only put initialization of VA_X_ME2_vr_mr put here

  observable_set save_osi = osi;
  phasespace_set save_psi = psi;
  runresumption_set rsi(osi, psi, save_osi, save_psi);

  osi.initialization_LHAPDF();
  osi.initialization_OpenLoops_process(psi);
  generic.phasespacepoint_psp(osi);
  calculate_ME2check_VJ_QCD(osi);

  if (osi.mass_parton[0][1] > 0. || osi.mass_parton[0][2] > 0.){osi.int_end = 1;}
  if (psi_n_events_max == 0){osi.int_end = 1;}
  if (osi.user.int_value[osi.user.int_map["rescaling_exponent"]] != 0){psi_hcf = psi_hcf * pow(osi.user.double_value[osi.user.double_map["rescaling_factor"]], osi.user.int_value[osi.user.int_map["rescaling_exponent"]]);}

  if (osi.int_end){osi.output_zero_contribution_complete(psi);}
  osi.initialization_runtime();
  while (osi.int_end == 0){
    rsi.perform_iteration_step();
    if (osi.int_end == 1){break;}

    psi.calculate_IS();
    
    psi.calculate_IS_CX();

    osi.calculate_IS_CX(psi);
    generic.ac_psp_psp(0, psi.MC_phasespace.channel, psi);
    osi.determine_p_parton(psi);
    // replace by proper nan-check !!!
    if (osi.p_parton[0] != osi.p_parton[0]){errorhandling_c_psp(); continue;}

    perform_event_selection(osi, generic);
    if (osi.cut_ps[0] == -1){psi.handling_cut_psp(); continue;}

    if (osi.switch_console_output_tau_0){psi.output_check_tau_0();}
    if (osi.switch_output_testpoint){osi.output_testpoint(psi);}
    ///    if (psi_i_acc == psi_n_events_max){osi.int_end = 1;}
    generic.ag_psp_psp(0, 0, psi);
    psi.calculate_g_tot();
    if (munich_isnan(psi_g_tot) || munich_isinf(psi_g_tot)){errorhandling_gtot(); continue;}

    generic.calculate_dynamic_scale(0, osi);
    generic.calculate_dynamic_scale_TSV(0, osi);
    osi.determine_scale();

    calculate_ME2_VJ_QCD(osi); // change order -> to before scale determination !!!
    osi.ME2 = osi.VA_V_ME2 + osi.VA_X_ME2;
    if (munich_isnan(osi.ME2)){errorhandling_me2(); continue;}

    if (osi.switch_moment){generic.moments(osi);}

    osi.calculate_pdf_LHAPDF_list_CV();
    if (osi.switch_TSV){osi.calculate_pdf_LHAPDF_list_TSV();}
    osi.determine_integrand_CX_ncollinear_VJ(psi);
    osi.determine_psp_weight_TSV(psi);
    osi.determine_psp_weight_QT(psi);

    if (osi.switch_output_maxevent){osi.output_integrand_maximum_VA(psi);}
    if (osi.switch_distribution){osi.determine_distribution_complete();}
    //    if (osi.switch_distribution){determine_distribution_complete(osi);}
    static double optimization_modifier = 1.;
    if (osi.user.switch_value[osi.user.switch_map["optimization_modifier"]]){
      optimization_modifier = pow(osi.particle_event[osi.access_object["Vrec"]][0][0].pT, osi.user.switch_value[osi.user.switch_map["optimization_modifier"]]);
      logger << LOG_DEBUG << "optimization_modifier = " << optimization_modifier << endl;
    }   

    psi.psp_MCweight_optimization(osi.integrand, osi.this_psp_weight, osi.this_psp_weight2, optimization_modifier);
    psi_i_acc++;
    ///    if (psi_i_acc % psi_n_step == 0){rsi.perform_iteration_step();}
    osi.determine_runtime(psi);
  }
  osi.output_finalization_integration(psi);
}
