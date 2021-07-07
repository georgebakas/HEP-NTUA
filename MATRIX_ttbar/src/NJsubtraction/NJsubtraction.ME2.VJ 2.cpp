#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"
#include "../include/definitions.phasespace.set.cxx"

// {{{ calculate_ME2_VJ_QCD(observable_set & oset)
void calculate_ME2_VJ_QCD(observable_set & oset){
  static Logger logger("calculate_ME2_VJ_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double one = 1;
  static int n_momentum = 5 * (osi_n_particle + 2);
  double *P;
  P = new double[n_momentum];
  for (int i = 1; i < osi_p_parton[0].size(); i++){
    P[5 * (i - 1)]     = osi_p_parton[0][i].x0();
    P[5 * (i - 1) + 1] = osi_p_parton[0][i].x1();
    P[5 * (i - 1) + 2] = osi_p_parton[0][i].x2();
    P[5 * (i - 1) + 3] = osi_p_parton[0][i].x3();
    P[5 * (i - 1) + 4] = osi_p_parton[0][i].m();
  }

  static double acc;
  static double M2L0;
  double *M2L1;
  double *IRL1;
  M2L1 = new double[3];
  IRL1 = new double[3];
  double *M2L2;
  double *IRL2;
  M2L2 = new double[5];
  IRL2 = new double[5];
  
  static char * OL_mu_ren = stch("muren");
  static char * OL_mu_reg = stch("mureg");
  static char * pole_uv = stch("pole_uv");
  static char * pole_ir1 = stch("pole_ir1");
  static char * pole_ir2 = stch("pole_ir2");
  static char * fact_uv = stch("fact_uv");
  static char * fact_ir = stch("fact_ir");

  ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
  ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
  ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);

  ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
  ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
  ol_setparameter_double(fact_uv, one);
  ol_setparameter_double(fact_ir, one);

  ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);
  osi_VA_V_ME2 = M2L1[0];

  for (int i_vr =  0; i_vr < oset.value_mu_ren.size(); i_vr++){
    for (int i_mr = 0; i_mr < oset.value_mu_ren[i_vr].size(); i_mr++){
      oset.VA_X_ME2_vr_mr[i_vr][i_mr] = 0.;
      double inv_factor_CV = oset.value_mu_ren[i_vr][i_mr] / oset.var_mu_ren;
      ol_setparameter_double(OL_mu_ren, oset.value_mu_ren[i_vr][i_mr]);
      ol_setparameter_double(OL_mu_reg, oset.value_mu_ren[i_vr][i_mr]);
      ol_setparameter_double(fact_uv, inv_factor_CV);
      ol_setparameter_double(fact_ir, inv_factor_CV);
      ol_evaluate_ct(1, P, &M2L0, &oset.VA_X_ME2_vr_mr[i_vr][i_mr]);
    }
  }
  // if (osi_switch_TSV){
  //   for (int i_v = 0; i_v < osi_max_dyn_ren + 1; i_v++){
  //     for (int i_r = 0; i_r < osi_n_scale_dyn_ren[i_v]; i_r++){
  // 	double inv_factor_TSV = osi_value_scale_ren[0][i_v][i_r] / osi_var_mu_ren;
  // 	ol_setparameter_double(OL_mu_ren, osi_value_scale_ren[0][i_v][i_r]);
  // 	ol_setparameter_double(OL_mu_reg, osi_value_scale_ren[0][i_v][i_r]);
  // 	ol_setparameter_double(fact_uv, inv_factor_TSV);
  // 	ol_setparameter_double(fact_ir, inv_factor_TSV);
  // 	ol_evaluate_ct(1, P, &M2L0, &osi_value_ME2term_ren[0][i_v][i_r]);
  // 	osi_value_ME2term_ren[0][i_v][i_r] += osi_VA_V_ME2 + osi_VA_I_ME2;
  // 	logger << LOG_DEBUG_VERBOSE << "after ol_evaluate_ct i_v = " << i_v << "   i_r = " << i_r << endl;
  //     }
  //   }
  // }
  
  delete [] P;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ calculate_ME2check_VJ_QCD(observable_set & oset)
void calculate_ME2check_VJ_QCD(observable_set & oset){
  static Logger logger("calculate_ME2check_VJ_QCD");
  // logger << LOG_DEBUG_VERBOSE << "called" << endl;

  osi_VA_delta_flag = 1;
  string sDelta;

  osi_VA_b_ME2 = 0.;
  osi_VA_V_ME2 = 0.;
  osi_VA_X_ME2 = 0.;
  osi_VA_I_ME2 = 0.;

  osi_QT_A0 = 0.;
  osi_QT_A1 = 0.;
  osi_QT_A2 = 0.;
  osi_QT_H1_delta = 0.;
  osi_QT_H2_delta = 0.;
  for (int i_s = 0; i_s < osi_n_scales_CV; i_s++){osi_VA_X_ME2_CV[i_s] = 0.;}

  if (osi_p_parton[0][0].x0() != 0.){
    ofstream out_comparison;
    out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
    osi_VA_DeltaUV = 0.;
    osi_VA_DeltaIR1 = 0.;
    osi_VA_DeltaIR2 = 0.;
    static char * OL_mu_ren = stch("muren");
    static char * OL_mu_reg = stch("mureg");
    static char * pole_uv = stch("pole_uv");
    static char * pole_ir1 = stch("pole_ir1");
    static char * pole_ir2 = stch("pole_ir2");
    ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
    ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
    ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
    ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
    ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
    calculate_ME2_VJ_QCD(oset);
    out_comparison << "Absolute results: " << endl << endl;
    oset.output_testpoint_VA_result(out_comparison);
    out_comparison << endl;
    out_comparison << "Particle momenta: " << endl << endl;
    output_momenta(out_comparison, oset);
    out_comparison.close();
  }
 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ calculate_ME2_VJ2_QCD(observable_set & oset, call_generic & generic)
void calculate_ME2_VJ2_QCD(observable_set & oset, call_generic & generic){
  static Logger logger("calculate_ME2_VJ2_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // static double acc;
  // static int n_momentum = 5 * (osi_n_particle + 2);
  // static double *P;
  // P = new double[n_momentum];
  // for (int i = 1; i < osi_p_parton[0].size(); i++){
  //   P[5 * (i - 1)]     = osi_p_parton[0][i].x0();
  //   P[5 * (i - 1) + 1] = osi_p_parton[0][i].x1();
  //   P[5 * (i - 1) + 2] = osi_p_parton[0][i].x2();
  //   P[5 * (i - 1) + 3] = osi_p_parton[0][i].x3();
  //   P[5 * (i - 1) + 4] = osi_p_parton[0][i].m();
  // }

  // double *M2L1;
  // double *IRL1;
  // M2L1 = new double[3];
  // IRL1 = new double[3];
  // double *M2L2;
  // double *IRL2;
  // M2L2 = new double[5];
  // IRL2 = new double[5];

  // double mu_Q = (osi_p_parton[0][1] + osi_p_parton[0][2]).m();
  // static double one = 1;

  // static char * OL_mu_ren = stch("muren");
  // static char * OL_mu_reg = stch("mureg");
  //   //  static char * renscale = stch("renscale");
  // static char * pole_uv = stch("pole_uv");
  // static char * pole_ir1 = stch("pole_ir1");
  // static char * pole_ir2 = stch("pole_ir2");
  // //  ol_setparameter_double(renscale, osi_var_mu_ren);
  // ///  ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
  // ///  ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
  // ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
  // ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
  // ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
  
  // int polenorm_0=0;
  // static char * polenorm = stch("polenorm");
  // ol_setparameter_int(polenorm, polenorm_0);
  
  // //  ol_setparameter_double(stch("renscale"), mu_Q);
  // ol_setparameter_double(OL_mu_ren, mu_Q);
  // ol_setparameter_double(OL_mu_reg, mu_Q);
  // static char * fact_uv = stch("fact_uv");
  // static char * fact_ir = stch("fact_ir");
  // ol_setparameter_double(fact_uv, one);
  // ol_setparameter_double(fact_ir, one);
  // //  ol_setparameter_double(stch("fact_uv"), one);
  // //  ol_setparameter_double(stch("fact_ir"), one);
  
  // //  ol_evaluate_loop2(1, P, M2L2, &acc);
  // ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);

  // // previous version !!! check if correct !!!
  // //  osi_VA_V_ME2 = M2L1[0];

  // static double M2L0;
  // ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2);
  // osi_VA_V_ME2 = M2L1[0]+osi_VA_X_ME2;
  
  // generic.calculate_H2(oset);
  
  // //  logger << LOG_DEBUG_VERBOSE << "ratios = " << osi_QT_A0 / osi_VA_b_ME2 << ", " << osi_QT_A1 / M2L1[0] << endl;
  // logger << LOG_DEBUG_VERBOSE << "ratios = " << osi_QT_A0 / osi_VA_b_ME2 << ", " << osi_QT_A1 / osi_VA_V_ME2 << endl;
  
  // delete [] M2L1;
  // delete [] IRL1;
  // delete [] M2L2;
  // delete [] IRL2;
  // delete [] P;

  // logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ calculate_ME2check_VJ2_QCD(observable_set & oset, call_generic & generic)
void calculate_ME2check_VJ2_QCD(observable_set & oset, call_generic & generic){
  static Logger logger("calculate_ME2check_VJ2_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // osi_VA_delta_flag = 1;
  // string sDelta;

  // osi_VA_b_ME2 = 0.;
  // osi_VA_V_ME2 = 0.;
  // osi_VA_X_ME2 = 0.;
  // osi_VA_I_ME2 = 0.;

  // osi_QT_A0 = 0.;
  // osi_QT_A1 = 0.;
  // osi_QT_A2 = 0.;
  // osi_QT_H1_delta = 0.;
  // osi_QT_H2_delta = 0.;
  // for (int i_s = 0; i_s < osi_n_scales_CV; i_s++){osi_VA_X_ME2_CV[i_s] = 0.;}

  // if (osi_p_parton[0][0].x0() != 0.){
  //   ofstream out_comparison;
  //   out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
  //   osi_VA_DeltaUV = 0.;
  //   osi_VA_DeltaIR1 = 0.;
  //   osi_VA_DeltaIR2 = 0.;
  //   static char * OL_mu_ren = stch("muren");
  //   static char * OL_mu_reg = stch("mureg");
  //   static char * pole_uv = stch("pole_uv");
  //   static char * pole_ir1 = stch("pole_ir1");
  //   static char * pole_ir2 = stch("pole_ir2");
  //   ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
  //   ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
  //   ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
  //   ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
  //   ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
  //   calculate_ME2_VJ2_QCD(oset, generic);
  //   out_comparison << "Absolute results: " << endl << endl;
  //   oset.output_testpoint_VA_result(out_comparison);
  //   out_comparison << endl;
  //   out_comparison << "Particle momenta: " << endl << endl;
  //   output_momenta(out_comparison, oset);
  //   out_comparison.close();
  // }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
