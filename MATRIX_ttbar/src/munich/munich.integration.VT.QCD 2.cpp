#include "../include/classes.cxx"
#include "../qTsubtraction/header/more.h"
void munich::integration_VT_QCD(){
  static Logger logger("munich::integration_VT_QCD");
  logger << LOG_INFO << "called" << endl;

  osi.QT_Qres = osi.user.cut_value[osi.user.cut_map["Qres"]];
  osi.switch_resum = osi.user.switch_value[osi.user.switch_map["do_resummation"]];
  osi.dynamical_Qres = osi.user.switch_value[osi.user.switch_map["dynamical_Qres"]];
  osi.QT_Qres_prefactor = osi.user.cut_value[osi.user.cut_map["Qres_prefactor"]];
  if (osi.QT_Qres_prefactor==0) {
    osi.QT_Qres_prefactor=1;
  }
  psi.do_resummation = osi.switch_resum;
  psi.Qres = osi.QT_Qres;
  psi.dynamical_Qres = osi.dynamical_Qres;
  psi.Qres_prefactor = osi.QT_Qres_prefactor;

  generic.determination_subprocess_psp(0, psi);
  generic.combination_subprocess_psp(0, psi, osi);

  osi.initialization_complete(isi, psi);

  psi.initialization_complete(isi);

  if (osi.switch_old_qT_version){psi.initialization_phasespace_IS_CT_from_CX();}

  osi.initialization_integration(psi);
  
  observable_set save_osi = osi;
  phasespace_set save_psi = psi;
  runresumption_set rsi(osi, psi, save_osi, save_psi);

  generic.phasespacepoint_psp(osi);
  calculate_ME2check_VT_QCD(osi, generic);

  if (osi.mass_parton[0][1] > 0. || osi.mass_parton[0][2] > 0.){osi.int_end = 1;}
  if (psi.n_events_max == 0){osi.int_end = 1;}
  if (osi.user.int_value[osi.user.int_map["rescaling_exponent"]] != 0){psi.hcf = psi.hcf * pow(osi.user.double_value[osi.user.double_map["rescaling_factor"]], osi.user.int_value[osi.user.int_map["rescaling_exponent"]]);}

  if (osi.int_end){osi.output_zero_contribution_complete(psi);}
  osi.initialization_runtime();
  while (osi.int_end == 0){
    rsi.perform_iteration_step();
    if (osi.int_end == 1){break;}

    psi.calculate_IS();
    
    // new version:
    if (!osi.switch_old_qT_version) {
      psi.calculate_IS_CX();
    }
    // old version:
    else if (osi.switch_old_qT_version){
      psi.calculate_IS_VT();
    }
    if (psi.do_resummation){psi.calculate_IS_QT();}

    osi.calculate_IS_CX(psi);
    generic.ac_psp_psp(0, psi.MC_phasespace.channel, psi);
    osi.determine_p_parton(psi);
    // replace by proper nan-check !!!
    if (osi.p_parton[0] != osi.p_parton[0]){errorhandling_c_psp(); continue;}

#ifdef MORE
    if (osi.switch_resummation){
      double Q = sqrt(psi.xbs_all[0][0]);
      double y = log(sqrt(psi.x_pdf[1] / psi.x_pdf[2]));
      performQTBoost_cms(psi.QT_qt2, Q, y, osi.p_parton);
    }
#endif

    perform_event_selection(osi, generic);
    if (osi.cut_ps[0] == -1){psi.handling_cut_psp(); continue;}

    if (osi.dynamical_Qres) {
      osi.QT_Qres = osi.QT_Qres_prefactor*sqrt(psi.xbs_all[0][0]);
    }

    if (osi.switch_console_output_tau_0){psi.output_check_tau_0();}
    if (osi.switch_output_testpoint){osi.output_testpoint(psi);}
    ///    if (psi.i_acc == psi.n_events_max){osi.int_end = 1;}
    generic.ag_psp_psp(0, 0, psi);
    psi.calculate_g_tot();
    if (munich_isnan(psi.g_tot) || munich_isinf(psi.g_tot)){errorhandling_gtot(); continue;}

    generic.calculate_dynamic_scale(0, osi);
    generic.calculate_dynamic_scale_TSV(0, osi);
    osi.determine_scale();

    calculate_ME2_VT_QCD(osi, generic); // change order -> to before scale determination !!!
    osi.ME2 = osi.VA_V_ME2 + osi.VA_X_ME2;
    if (munich_isnan(osi.ME2)){errorhandling_me2(); continue;}

    if (osi.switch_moment){generic.moments(osi);}

    // old version:
    if (osi.switch_old_qT_version) {
      calculate_pdf_LHAPDF_QT(osi.combination_pdf, psi, osi, psi.contribution_order_alpha_s[0]);
      osi.determine_integrand_VT(psi);
      osi.determine_psp_weight_VT(psi);
    }
    // new version:
    if (!osi.switch_old_qT_version){
      osi.calculate_pdf_LHAPDF_list_CV();
      if (osi.switch_TSV){osi.calculate_pdf_LHAPDF_list_TSV();}
      osi.determine_integrand_CX_ncollinear_VT(psi);
      osi.determine_psp_weight_TSV(psi);
      osi.determine_psp_weight_QT(psi);
    }

    if (osi.switch_output_maxevent){osi.output_integrand_maximum_VA(psi);}
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
