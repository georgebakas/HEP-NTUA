#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

// {{{ observable_set::determine_integrand_CX_ncollinear_CJ(phasespace_set & psi)
void observable_set::determine_integrand_CX_ncollinear_CJ(phasespace_set & psi){
  Logger logger("observable_set::determine_integrand_CX_ncollinear_CJ");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  // static int order = 1;
  // double TauNoverQ = sqrt(psi_QT_qt2 / QT_Qres / QT_Qres); // FIXME: comput N-jettiness TauN here
  // double LL1, LL2, LL3, LL4;

  // first: do implementation with binning in TauN, later try to integrate logs directly
  // LL0 = NJ.compute_NJ_log(TauNoverQ, 0);
  // LL1 = NJ.compute_NJ_log(TauNoverQ, 1);
  // LL2 = NJ.compute_NJ_log(TauNoverQ, 2);
  // LL3 = NJ.compute_NJ_log(TauNoverQ, 3);

  double global_factor = alpha_S / 2. / pi; // note: different expansion than in qT subtraction!

  vector<double> conv_coeff_Pxx,conv_coeff_Ixx;
  conv_coeff_Pxx.resize(list_combination_pdf.size());
  conv_coeff_Ixx.resize(list_combination_pdf.size());
  // compute terms convoluted over addition z integration with PDFs for subtraction terms
  determine_NJconvolution_terms_1(psi,conv_coeff_Pxx,conv_coeff_Ixx);

  // double LQ = log(psi_xbs_all[0][0] / pow(QT_Qres, 2));
  // if (QT_Qres == 0){LQ = 0.;}

  double Q_a,Q_b; // FIXME: have to set it here correctly !!!

  vector<vector<double> > convolutions_Pxx,convolutions_Ixx,value_C11,value_C10;
  convolutions_Pxx.resize(value_mu_fact.size());
  convolutions_Ixx.resize(value_mu_fact.size());
  value_C11.resize(value_mu_fact.size());
  value_C10.resize(value_mu_fact.size());
  vector<vector<vector<vector<double> > > > value_C10_muR;
  value_C10_muR.resize(value_mu_ren.size(), vector<vector<vector<double> > > (value_mu_fact.size()));

  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    convolutions_Pxx[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    convolutions_Ixx[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    // coefficients without born amplitude; multiplied with LO squared amplitude below
    value_C11[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    value_C10[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
      value_C10_muR[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<double> (value_mu_fact[i_vf].size(), double (ncollinear.size())));
    }
  }

  Q_a = 1.; // FIXME: have to set it here correctly !!!
  Q_b = 1.; // FIXME: have to set it here correctly !!!
  NJ.sab = (p_parton[0][1]+p_parton[0][2]).m2(); // (p1+p2)^2
  int jetindex = 0;  // FIXME: remove hardcoded
  if (jetindex) {
    NJ.sa1 = (p_parton[0][1]-p_parton[0][jetindex]).m2(); // (p1-pj)^2
    NJ.s1b = (p_parton[0][2]-p_parton[0][jetindex]).m2(); // (p2-pj)^2
  }

  // loop over factorization scales and Sigma terms for subtraction
  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){ // loop over muF scale variation
    for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){ // loop over different muF scales?
      for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){ // loops over PDF contributions to compute final convolutions; 
                                                                   // done here because coefficient muF independent, but PDF factors muF dependend
  	convolutions_Pxx[i_vf][i_mf] += value_list_pdf_factor[i_vf][i_mf][i_l][0] * conv_coeff_Pxx[i_l];
  	convolutions_Ixx[i_vf][i_mf] += value_list_pdf_factor[i_vf][i_mf][i_l][0] * conv_coeff_Ixx[i_l];
      }
      // compute muR independent parts
      value_C11[i_vf][i_mf] = NJ.calculate_C11(value_list_pdf_factor[i_vf][i_mf][0][0]);
      value_C10[i_vf][i_mf] = NJ.calculate_C10(value_list_pdf_factor[i_vf][i_mf][0][0],convolutions_Pxx[i_vf][i_mf],convolutions_Ixx[i_vf][i_mf]);

      for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	  //compute muR dependent parts separately
	  value_C10_muR[i_vr][i_vf][i_mr][i_mf] = NJ.calculate_C10_muR(value_list_pdf_factor[i_vf][i_mf][0][0],convolutions_Pxx[i_vf][i_mf],convolutions_Ixx[i_vf][i_mf],Q_a,Q_b,value_mu_ren[i_vr][i_mr]);
	}
      }


      // WHAT IS i_y ???
      // WHY IS THERE AN EXTRA INTEGRANT FOR DISTRIBUTIONS ???


      // PUT BACK FOR DISTRIBUTIONS!

      // if (switch_distribution){
      //        for (int i_y = 0; i_y < 3; i_y++){
      //          for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
      //            coll_tH1F_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tH1F[i_l];
      //          }
      //          value_tH1F_D[i_vf][i_mf][i_y] = accumulate(coll_tH1F_pdf_D[i_vf][i_mf][i_y].begin(), coll_tH1F_pdf_D[i_vf][i_mf][i_y].end(), 0.);
      //          value_sig12_D[i_vf][i_mf][i_y] = calculate_sigma12(value_list_pdf_factor[i_vf][i_mf][0][i_y]);
      //          value_sig11_D[i_vf][i_mf][i_y] = calculate_sigma11(value_list_pdf_factor[i_vf][i_mf][0][i_y], value_tH1F_D[i_vf][i_mf][i_y], LQ);
      //        }
      // }

      // PUT BACK FOR DISTRIBUTIONS!
    }
  }

  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
      logger << LOG_DEBUG_VERBOSE << "value_C11    [" << i_vf << "]   [" << i_mf << "]    = " << value_C11[i_vf][i_mf] << endl;
      logger << LOG_DEBUG_VERBOSE << "value_C10    [" << i_vf << "]   [" << i_mf << "]    = " << value_C10[i_vf][i_mf] << endl;
      logger.newLine(LOG_DEBUG_VERBOSE);
    }
  }


  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
      for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	  for (int i_q = 0; i_q < n_qTcut; i_q++){
	    CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf][i_mr][i_mf] = global_factor * (value_C11[i_vf][i_mf] * LL1_int[i_q] + (value_C10[i_vf][i_mf]+value_C10_muR[i_vr][i_vf][i_mr][i_mf]) * LL0_int[i_q]);
	    logger << LOG_DEBUG_VERBOSE << "CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf][i_mr][i_mf] = " << CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf][i_mr][i_mf] << endl;
	    logger << LOG_DEBUG_VERBOSE << "single factors of it: " << global_factor << " " << value_C11[i_vf][i_mf]  << " " << LL1_int[i_q] << " " << value_C10[i_vf][i_mf] << " " << value_C10_muR[i_vr][i_vf][i_mr][i_mf] << " " <<  LL0_int[i_q] << endl;
	    // PUT BACK FOR DISTRIBUTIONS!
	    // if (switch_distribution){
	    //   for (int i_y = 0; i_y < 3; i_y++){
	    //        CX_value_integrand_RF_qTcut_D[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = global_factor * (value_sig12_D[i_vf][i_mf][i_y] * I2_int[i_q] + value_sig11_D[i_vf][i_mf][i_y] * I1_int[i_q]);
	    //   }
	    // }
	    // PUT BACK FOR DISTRIBUTIONS!
	  }
	}
      }
    }
  }


    // WHAT HAPPENS HERE ???

  integrand = -var_rel_alpha_S * psi_ps_factor * ME2 * CX_value_integrand_RF_qTcut[0][dynamic_scale][dynamic_scale][map_value_scale_ren][map_value_scale_fact];

  if (switch_CV){
    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      integrand_CV[i_s] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * CX_value_integrand_RF_qTcut[0][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]];
      for (int i_q = 0; i_q < n_qTcut; i_q++){
        integrand_qTcut_CV[i_q][i_s] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * CX_value_integrand_RF_qTcut[i_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]];
      }
    }
  }
  
  // if (switch_distribution){
  //   int x_q = 0;
  //   for (int j = 0; j < 3; j++){
  //     integrand_D[j][0] = -var_rel_alpha_S * psi_ps_factor * ME2 * CX_value_integrand_RF_qTcut_D[x_q][dynamic_scale][map_value_scale_fact][j];
  //     logger << LOG_DEBUG_VERBOSE << "new: integrand_D[" << j << "][0] = " << integrand_D[j][0] << endl;
  //     if (switch_CV){
  //       for (int i_s = 0; i_s < n_scales_CV; i_s++){
  //         integrand_D_CV[i_s][j][0] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * CX_value_integrand_qTcut_D[x_q][dynamic_scale_CV][map_value_scale_fact_CV[i_s]][j];
  //         for (int i_q = 0; i_q < n_qTcut; i_q++){
  //           integrand_D_qTcut_CV[i_q][i_s][j][0] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * CX_value_integrand_RF_qTcut_D[i_q][dynamic_scale_CV][map_value_scale_fact_CV[i_s]][j];
  //         }
  //       }
  //     }
  //   }
  //  }


  // WHY IS THIS DONE AGAIN for TSV ???

  // put back later
  // if (switch_TSV){
  //   for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
  //     for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
  //       for (int i_y = 0; i_y < 3; i_y++){
  //         for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
  //           coll_tH1F_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tH1F[i_l];
  //         }
  //         value_tH1F_TSV[i_vf][i_mf][i_y] = accumulate(coll_tH1F_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tH1F_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
  //         value_sig12_TSV[i_vf][i_mf][i_y] = calculate_sigma12(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y]);
  //         value_sig11_TSV[i_vf][i_mf][i_y] = calculate_sigma11(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y], value_tH1F_TSV[i_vf][i_mf][i_y], LQ);

  //         for (int i_q = 0; i_q < n_qTcut; i_q++){
  //           if (switch_resummation){
  //             if (i_q <= cut_ps[0]){
  //               value_fact_integrand_qTcut_TSV[i_q][i_vf][i_mf][i_y] = -ME2 * global_factor * (value_sig12_TSV[i_vf][i_mf][i_y] * LL2 + value_sig11_TSV[i_vf][i_mf][i_y] * LL1) / pow(QT_Qres,2) * psi_QT_jacqt2;
  //             }
  //             else {value_fact_integrand_qTcut_TSV[i_q][i_vf][i_mf][i_y] = 0.;}
  //           }
  //           else {
  //             value_fact_integrand_qTcut_TSV[i_q][i_vf][i_mf][i_y] = -ME2 * global_factor * (value_sig12_TSV[i_vf][i_mf][i_y] * I2_int[i_q] + value_sig11_TSV[i_vf][i_mf][i_y] * I1_int[i_q]);
  //           }
  //         }
  //       }
  //     }
  //   }

  // WHY IS THIS DONE AGAIN for TSV ???

  //   for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
  //     for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
  //       logger << LOG_DEBUG_VERBOSE << "value_sig11_TSV    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig11_TSV[i_vf][i_mf][0] << endl;
  //       logger << LOG_DEBUG_VERBOSE << "value_sig12_TSV    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig12_TSV[i_vf][i_mf][0] << endl;
  //       logger << LOG_DEBUG_VERBOSE << "value_tH1F_TSV     [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1F_TSV[i_vf][i_mf][0] << endl;
  //     }
  //   }
  // }
  // put back later

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ observable_set::determine_NJconvolution_terms_1(phasespace_set &psi,vector<double> &conv_coeff_Pxx,vector<double> &conv_coeff_Ixx)
void observable_set::determine_NJconvolution_terms_1(phasespace_set &psi,vector<double> &conv_coeff_Pxx,vector<double> &conv_coeff_Ixx){
  // computes convolution coefficients without the PDF factor
  static Logger logger("observable_set::determine_NJconvolution_terms_1");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  // first empty previous convolution results
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    conv_coeff_Pxx[i_l] = 0.;
    conv_coeff_Ixx[i_l] = 0.;
  }
  for (int i_c = 1; i_c < ncollinear.size(); i_c++){
    int x_e = ncollinear[i_c].x_a;
    int y_e = ncollinear[i_c].x_b;
    int no_zx = ncollinear[ncollinear[i_c].no_endpoint[y_e]].no_pdf;
    int no_xx = ncollinear[ncollinear[i_c].no_endpoint[0]].no_pdf;

    logger << LOG_DEBUG_VERBOSE << "ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << endl;
    if (x_e == 0){logger << LOG_FATAL << "Should not happen at 1st order!!!" << endl; exit(1);}
    if (ncollinear[i_c].type_splitting_full[x_e] == 0){logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = 0 -> no splitting type defined !!!" << endl; exit(1);}
    // depending on type of splitting (ncollinear[i_c]) compute different convolutions
    // type_splitting_full (leg-wise):
    // x_e - y_e     ---  ---

    //   1 -   0   +   0 -  1   // hard process with g from g -> g (g) splitting
    //   2 -   0   +   0 -  2   // hard process with q from q -> q (g) splitting
    //   3 -   0   +   0 -  3   // hard process with g from q -> g (q) splitting
    //   4 -   0   +   0 -  4   // hard process with q from g -> q (qx) splitting       

    else if (ncollinear[i_c].type_splitting_full[x_e] == 1){ // g -> g (g) splitting
      // compute convolutions with splitting function: P(z1) x f(x1/z1) x f(x2)  (PDF factors multiplied later)
      conv_coeff_Pxx[no_zx] += NJ.P1gg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + NJ.P1gg_z(psi_z_coll[x_e]) / psi_g_z_coll[x_e];
// conv_coeff_Pxx[no_zx] += (3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] + Pggreg(psi_z_coll[x_e]) / psi_g_z_coll[x_e])*2.;
      logger << LOG_DEBUG_VERBOSE << "Pgg: compare these two, they should be equal:" << endl;
      logger << LOG_DEBUG_VERBOSE << NJ.P1gg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + NJ.P1gg_z(psi_z_coll[x_e]) / psi_g_z_coll[x_e] << endl;
      logger << LOG_DEBUG_VERBOSE << (3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] + Pggreg(psi_z_coll[x_e]) / psi_g_z_coll[x_e])*2. << endl;
      // -log( x1 ) * (pdf_factor_z1x2 * Pggreg(z1)
      // -log( x1 ) * (pdf_factor_z1x2)*3/(1-z1)

      // compute delta piece of distributions with splitting function: P(z1) x f(x1) x f(x2)  (PDF factors multiplied later)
      conv_coeff_Pxx[no_xx] += -psi_z_coll[x_e] * NJ.P1gg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.P1gg_int(x_pdf[x_e]) + NJ.P1gg_delta(psi_z_coll[x_e]);
// conv_coeff_Pxx[no_xx] += (-psi_z_coll[x_e] * 3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(x_pdf[x_e]) + beta0)*2.;
      logger << LOG_DEBUG_VERBOSE << "Pgg: compare these two, they should be equal:" << endl;
      logger << LOG_DEBUG_VERBOSE << -psi_z_coll[x_e] * NJ.P1gg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.P1gg_int(x_pdf[x_e]) + NJ.P1gg_delta(psi_z_coll[x_e]) << endl;
      logger << LOG_DEBUG_VERBOSE << (-psi_z_coll[x_e] * 3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(x_pdf[x_e]) + beta0)*2. << endl;
      // -log( x1 ) * (- pdf_factor_x1x2 * z1 * 3 / (1 - z1)) - 3 * D0int(x1) * pdf_factor_x1x2
      // + beta0 * pdf_factor_x1x2 (half of this (because added for each leg)!)
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_P1_contribution[no_zx = " << no_zx << "] += " << conv_coeff_Pxx[no_zx] << endl;
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_P1_contribution[no_xx = " << no_xx << "] += " << conv_coeff_Pxx[no_xx] << endl;
 
      // now the same for matching functions of 0-jettiness Beam functions
      conv_coeff_Ixx[no_zx] += NJ.I1gg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + NJ.I1gg_z(psi_z_coll[x_e]) / psi_g_z_coll[x_e];
      conv_coeff_Ixx[no_xx] += -psi_z_coll[x_e] * NJ.I1gg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.I1gg_int(x_pdf[x_e]) + NJ.I1gg_delta(psi_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 2){ // hard process with q from q -> q (g) splitting
      // compute convolutions with splitting function: P(z1) x f(x1/z1) x f(x2)  (PDF factors multiplied later)
      conv_coeff_Pxx[no_zx] += NJ.P1qq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + NJ.P1qq_z(psi_z_coll[x_e]) / psi_g_z_coll[x_e];
// conv_coeff_Pxx[no_zx] += (Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e])*2.;
      logger << LOG_DEBUG_VERBOSE << "Pqq: compare these two, they should be equal:" << endl;
      logger << LOG_DEBUG_VERBOSE << -psi_z_coll[x_e] * NJ.P1qq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.P1qq_int(x_pdf[x_e]) + NJ.P1qq_delta(psi_z_coll[x_e]) << endl;
      logger << LOG_DEBUG_VERBOSE << (Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e])*2. << endl;
      // - g_z1 * (pdf_factor_z1x2) * Pqq(z1)  // from H1F

      // compute delta piece of distributions with splitting function: P(z1) x f(x1) x f(x2)  (PDF factors multiplied later)
      conv_coeff_Pxx[no_xx] += -psi_z_coll[x_e] * NJ.P1qq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.P1qq_int(x_pdf[x_e]) + NJ.P1qq_delta(psi_z_coll[x_e]);
// conv_coeff_Pxx[no_xx] += ((-psi_z_coll[x_e] * Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - Pqqint(x_pdf[x_e])))*2.;
      logger << LOG_DEBUG_VERBOSE << "Pqq: compare these two, they should be equal:" << endl;
      logger << LOG_DEBUG_VERBOSE << -psi_z_coll[x_e] * NJ.P1qq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.P1qq_int(x_pdf[x_e]) + NJ.P1qq_delta(psi_z_coll[x_e]) << endl;
      logger << LOG_DEBUG_VERBOSE << ((-psi_z_coll[x_e] * Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - Pqqint(x_pdf[x_e])))*2. << endl;
      // - g_z1 * (- pdf_factor_x1x2 * z1) * Pqq(z1) - pdf_factor_x1x2 * Pqqint(x1)  // from H1F
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_P1_contribution[no_zx = " << no_zx << "] += " << conv_coeff_Pxx[no_zx] << endl;
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_P1_contribution[no_xx = " << no_xx << "] += " << conv_coeff_Pxx[no_xx] << endl;

      // now the same for matching functions of 0-jettiness Beam functions
      conv_coeff_Ixx[no_zx] += NJ.I1qq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + NJ.I1qq_z(psi_z_coll[x_e]) / psi_g_z_coll[x_e];
      conv_coeff_Ixx[no_xx] += -psi_z_coll[x_e] * NJ.I1qq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.I1qq_int(x_pdf[x_e]) + NJ.I1qq_delta(psi_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 3){ // hard process with g from q -> g (q) splitting
      // compute convolutions with splitting function: P(z1) x f(x1/z1) x f(x2)  (PDF factors multiplied later)
      conv_coeff_Pxx[no_zx] += NJ.P1gq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + NJ.P1gq_z(psi_z_coll[x_e]) / psi_g_z_coll[x_e];
// conv_coeff_Pxx[no_zx] += (Pgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e])*2.;
      logger << LOG_DEBUG_VERBOSE << "Pgq: compare these two, they should be equal:" << endl;
      logger << LOG_DEBUG_VERBOSE << NJ.P1gq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + NJ.P1gq_z(psi_z_coll[x_e]) / psi_g_z_coll[x_e] << endl;
      logger << LOG_DEBUG_VERBOSE << (Pgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e])*2. << endl;
      // -log( x1 ) * (pdf_factor_qx2 * Pgq(z1) )
      logger << LOG_DEBUG_VERBOSE << "q->g(q):  coll_P1_contribution[no_zx = " << no_zx << "] += " << conv_coeff_Pxx[no_zx] << endl;

      conv_coeff_Pxx[no_xx] += -psi_z_coll[x_e] * NJ.P1gq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.P1gq_int(x_pdf[x_e]) + NJ.P1gq_delta(psi_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "Pgq: compare these two, they should be equal:" << endl;
      logger << LOG_DEBUG_VERBOSE << -psi_z_coll[x_e] * NJ.P1gq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.P1gq_int(x_pdf[x_e]) + NJ.P1gq_delta(psi_z_coll[x_e]) << endl;
      logger << LOG_DEBUG_VERBOSE << 0.*2. << endl;

      // now the same for matching functions of 0-jettiness Beam functions
      conv_coeff_Ixx[no_zx] += NJ.I1gq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + NJ.I1gq_z(psi_z_coll[x_e]) / psi_g_z_coll[x_e];
      conv_coeff_Ixx[no_xx] += -psi_z_coll[x_e] * NJ.I1gq_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.I1gq_int(x_pdf[x_e]) + NJ.I1gq_delta(psi_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 4){ // hard process with q from g -> q (qx) splitting
      // compute convolutions with splitting function: P(z1) x f(x1/z1) x f(x2)  (PDF factors multiplied later)
      conv_coeff_Pxx[no_zx] += NJ.P1qg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + NJ.P1qg_z(psi_z_coll[x_e]) / psi_g_z_coll[x_e];
// conv_coeff_Pxx[no_zx] += (Pqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e])*2.;
      logger << LOG_DEBUG_VERBOSE << "Pqg: compare these two, they should be equal:" << endl;
      logger << LOG_DEBUG_VERBOSE << NJ.P1qg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + NJ.P1qg_z(psi_z_coll[x_e]) / psi_g_z_coll[x_e] << endl;
      logger << LOG_DEBUG_VERBOSE << (Pqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e])*2. << endl;
      // - g_z1 * Pqg(z1) * (pdf_factor_gx2)  // from H1F
      logger << LOG_DEBUG_VERBOSE << "g->q(qx): coll_P1_contribution[no_zx = " << no_zx << "] += " << conv_coeff_Pxx[no_zx] << endl;

      conv_coeff_Pxx[no_xx] += -psi_z_coll[x_e] * NJ.P1qg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.P1qg_int(x_pdf[x_e]) + NJ.P1qg_delta(psi_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "Pqg: compare these two, they should be equal:" << endl;
      logger << LOG_DEBUG_VERBOSE << -psi_z_coll[x_e] * NJ.P1qg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.P1qg_int(x_pdf[x_e]) + NJ.P1qg_delta(psi_z_coll[x_e]) << endl;
      logger << LOG_DEBUG_VERBOSE << 0.*2. << endl;

      // now the same for matching functions of 0-jettiness Beam functions
      conv_coeff_Ixx[no_zx] += NJ.I1qg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + NJ.I1qg_z(psi_z_coll[x_e]) / psi_g_z_coll[x_e];
      conv_coeff_Ixx[no_xx] += -psi_z_coll[x_e] * NJ.I1qg_plus(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - NJ.I1qg_int(x_pdf[x_e]) + NJ.I1qg_delta(psi_z_coll[x_e]);
    }
    else {logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << " -> selected splitting type not allowed !!!" << endl; exit(1);}
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
