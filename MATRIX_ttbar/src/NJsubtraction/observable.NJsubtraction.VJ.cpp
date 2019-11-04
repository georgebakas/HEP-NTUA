#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

// {{{ observable_set::determine_integrand_CX_ncollinear_VJ(phasespace_set & psi)
void observable_set::determine_integrand_CX_ncollinear_VJ(phasespace_set & psi){
  Logger logger("observable_set::determine_integrand_CX_ncollinear_VJ");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
     
  double global_factor = alpha_S / 2. / pi; // note: different expansion than in qT subtraction!

  vector<double> conv_coeff_Pxx,conv_coeff_Ixx;
  conv_coeff_Pxx.resize(list_combination_pdf.size());
  conv_coeff_Ixx.resize(list_combination_pdf.size());
  // compute terms convoluted over addition z integration with PDFs for subtraction terms
  determine_NJconvolution_terms_1(psi,conv_coeff_Pxx,conv_coeff_Ixx);

  // double LQ = log(psi_xbs_all[0][0] / pow(QT_Qres, 2));
  // if (QT_Qres == 0){LQ = 0.;}

  vector<vector<double> > convolutions_Pxx,convolutions_Ixx,value_C1m1;
  convolutions_Pxx.resize(value_mu_fact.size());
  convolutions_Ixx.resize(value_mu_fact.size());
  value_C1m1.resize(value_mu_fact.size());
  vector<vector<vector<vector<double> > > > value_C1m1_muR;
  value_C1m1_muR.resize(value_mu_ren.size(), vector<vector<vector<double> > > (value_mu_fact.size()));

  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    convolutions_Pxx[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    convolutions_Ixx[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    // coefficients without born amplitude; multiplied with LO squared amplitude below
    value_C1m1[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
      value_C1m1_muR[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<double> (value_mu_fact[i_vf].size(), double (ncollinear.size())));
    }
  }

  // loop over factorization scales and Sigma terms for subtraction
  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){ // loop over muF scale variation
    for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){ // loop over different muF scales?
      // compute logarithms ln(Q/muF)
      for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){ // loops over PDF contributions to compute final convolutions; 
                                                                   // done here because coefficient muF independent, but PDF factors muF dependend
  	convolutions_Pxx[i_vf][i_mf] += value_list_pdf_factor[i_vf][i_mf][i_l][0] * conv_coeff_Pxx[i_l];
  	convolutions_Ixx[i_vf][i_mf] += value_list_pdf_factor[i_vf][i_mf][i_l][0] * conv_coeff_Ixx[i_l];
      }
      value_C1m1[i_vf][i_mf] = NJ.calculate_C1m1(value_list_pdf_factor[i_vf][i_mf][0][0],convolutions_Pxx[i_vf][i_mf],convolutions_Ixx[i_vf][i_mf]);
      for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	  double HardF = (VA_V_ME2+VA_X_ME2_vr_mr[i_vr][i_mr])/VA_b_ME2; // FIXME: compute hard function correctly
	  double lnRF = log(value_mu_ren[i_vr][i_mr] / value_mu_fact[i_vf][i_mf]);
	  double Q_a = 1.; // FIXME: have to set it here correctly !!!
	  double Q_b = 1.; // FIXME: have to set it here correctly !!!
	  //compute muR dependent parts separately
	  value_C1m1_muR[i_vr][i_vf][i_mr][i_mf] = NJ.calculate_C1m1_muR(value_list_pdf_factor[i_vf][i_mf][0][0],HardF,convolutions_Pxx[i_vf][i_mf],convolutions_Ixx[i_vf][i_mf],lnRF,Q_a,Q_b,value_mu_ren[i_vr][i_mr]);
	  for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	    CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf][i_mr][i_mf] = global_factor * (value_C1m1[i_vf][i_mf]+value_C1m1_muR[i_vr][i_vf][i_mr][i_mf]);
	  }
	}
      }


      // if (switch_distribution){
      //   for (int i_y = 0; i_y < 3; i_y++){
      // 	  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
      // 	    coll_tH1F_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tH1F[i_l];
      // 	    coll_tH1_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tH1[i_l];
      // 	  }
      // 	  value_tH1F_D[i_vf][i_mf][i_y] = accumulate(coll_tH1F_pdf_D[i_vf][i_mf][i_y].begin(), coll_tH1F_pdf_D[i_vf][i_mf][i_y].end(), 0.);
      // 	  value_tH1_D[i_vf][i_mf][i_y] = accumulate(coll_tH1_pdf_D[i_vf][i_mf][i_y].begin(), coll_tH1_pdf_D[i_vf][i_mf][i_y].end(), 0.);
      //   }
      // }
    }
  }

  // for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
  //   for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
  //     for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
  // 	for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
  // 	  logger << LOG_DEBUG_VERBOSE << "QT_H1_delta                 = " << QT_H1_delta << endl;
  // 	  logger << LOG_DEBUG_VERBOSE << "value_LR          [" << i_vr << "]   [" << i_mr << "] = " << value_LR[i_vr][i_mr] << endl;
  // 	  logger << LOG_DEBUG_VERBOSE << "value_LF       [" << i_vf << "]   [" << i_mf << "]    = " << value_LF[i_vf][i_mf] << endl;
  // 	  logger << LOG_DEBUG_VERBOSE << "value_tH1      [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1[i_vf][i_mf] << endl;
  // 	  logger << LOG_DEBUG_VERBOSE << "value_tH1F     [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1F[i_vf][i_mf] << endl;
  // 	}
  //     }
  //   }
  // }


  int x_q = 0;
  integrand = var_rel_alpha_S * psi_ps_factor * VA_b_ME2 * CX_value_integrand_RF_qTcut[x_q][dynamic_scale][dynamic_scale][map_value_scale_ren][map_value_scale_fact];
  if (switch_CV){
    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      integrand_CV[i_s] = var_rel_alpha_S_CV[i_s] * psi_ps_factor * VA_b_ME2 * CX_value_integrand_RF_qTcut[x_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]];
      for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	integrand_qTcut_CV[i_q][i_s] = var_rel_alpha_S_CV[i_s] * psi_ps_factor * VA_b_ME2 * CX_value_integrand_RF_qTcut[i_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]];
      }
    }
  }
    
  // if (switch_distribution){
  //   for (int i_y = 0; i_y < 3; i_y++){
  //     integrand_D[i_y][0] = var_rel_alpha_S * psi_ps_factor * VA_b_ME2 * CX_value_integrand_RF_qTcut_D[x_q][dynamic_scale][dynamic_scale][map_value_scale_ren][map_value_scale_fact][i_y];
  //     if (switch_CV){
  // 	for (int i_s = 0; i_s < n_scales_CV; i_s++){
  // 	  integrand_D_CV[i_s][i_y][0] = var_rel_alpha_S_CV[i_s] * psi_ps_factor * VA_b_ME2 * CX_value_integrand_RF_qTcut_D[x_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]][i_y];
  // 	  for (int i_q = 0; i_q < output_n_qTcut; i_q++){
  // 	    integrand_D_qTcut_CV[i_q][i_s][i_y][0] = var_rel_alpha_S_CV[i_s] * psi_ps_factor * VA_b_ME2 * CX_value_integrand_RF_qTcut_D[i_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]][i_y];
  // 	  }
  // 	}
  //     }
  //   }
  // }



  // if (switch_TSV){
  //   for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
  //     for (int i_mr = 0; i_mr < n_scale_dyn_ren[i_vr]; i_mr++){
  // 	value_LR_TSV[i_vr][i_mr] = log(psi_xbs_all[0][0] / pow(value_scale_ren[0][i_vr][i_mr], 2));
  //     }
  //   }
  //   for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
  //     for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
  // 	value_LF_TSV[i_vf][i_mf] = log(psi_xbs_all[0][0] / pow(value_scale_fact[0][i_vf][i_mf], 2));
  // 	for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
  // 	  for (int i_y = 0; i_y < 3; i_y++){
  // 	    coll_tH1F_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tH1F[i_l];
  // 	    coll_tH1_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tH1[i_l];
  // 	  }
  // 	}
  // 	for (int i_y = 0; i_y < 3; i_y++){
  // 	  value_tH1F_TSV[i_vf][i_mf][i_y] = accumulate(coll_tH1F_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tH1F_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
  // 	  value_tH1_TSV[i_vf][i_mf][i_y] = accumulate(coll_tH1_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tH1_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
  // 	}
  //     }
  //   }

  //   for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
  //     for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
  // 	for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
  // 	  for (int i_mr = 0; i_mr < n_scale_dyn_ren[i_vr]; i_mr++){
  // 	    logger << LOG_DEBUG_VERBOSE << "QT_H1_delta                     = " << QT_H1_delta << endl;
  // 	    logger << LOG_DEBUG_VERBOSE << "value_LR_TSV          [" << i_vr << "]   [" << i_mr << "] = " << value_LR_TSV[i_vr][i_mr] << endl;
  // 	    logger << LOG_DEBUG_VERBOSE << "value_LF_TSV       [" << i_vf << "]   [" << i_mf << "]    = " << value_LF_TSV[i_vf][i_mf] << endl;
  // 	    logger << LOG_DEBUG_VERBOSE << "value_tH1_TSV      [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1_TSV[i_vf][i_mf][0] << endl;
  // 	    logger << LOG_DEBUG_VERBOSE << "value_tH1F_TSV     [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1F_TSV[i_vf][i_mf][0] << endl;
  // 	    logger.newLine(LOG_DEBUG_VERBOSE);
  //  	  }
  // 	}
  //     }
  //   }


  //   for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
  //     for (int i_mr = 0; i_mr < n_scale_dyn_ren[i_vr]; i_mr++){
  // 	for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
  // 	  for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
  // 	    for (int i_y = 0; i_y < 3; i_y++){
  // 	      for (int i_q = 0; i_q < output_n_qTcut; i_q++){
  // 		value_integrand_qTcut_TSV[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = global_factor * VA_b_ME2 * calculate_H1(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y], value_tH1_TSV[i_vf][i_mf][i_y], value_tH1F_TSV[i_vf][i_mf][i_y], value_LR_TSV[i_vr][i_mr], value_LF_TSV[i_vf][i_mf], LQ);
  // 	      }
  // 	    }
  // 	  }
  // 	}
  //     }
  //   }
  // }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
