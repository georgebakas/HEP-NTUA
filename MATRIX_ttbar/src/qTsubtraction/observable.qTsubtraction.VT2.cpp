#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

void observable_set::determine_integrand_CX_ncollinear_VT2(phasespace_set & psi){
  Logger logger("observable_set::determine_integrand_CX_ncollinear_VT2");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
#ifdef MORE
  static int initialization = 1;
  static int gg_initiated = 0;
  static int charged = 0;
  if (initialization){
    string initial_state = name_process.substr(0,3);
    if (initial_state == "uu~") {
      gg_initiated = -1;
    }
    else if (initial_state == "dd~") {
      gg_initiated = -2;
    }
    else if (initial_state == "bb~") {
      gg_initiated = -12;
    }
    else if (initial_state == "gg_") {
      gg_initiated = 1;
    }
    else if (initial_state == "du~") {
      charged = 1;
      gg_initiated = -1;
    }
    else if (initial_state == "ud~") {
      charged = 1;
      gg_initiated = -2;
    }
    else {
      assert(false && "initial state for resummation not known");
    }
    initialization = 0;
  }
#endif
   
  //  CX_value_integrand_qTcut  must be redefined !!!

  double global_factor = pow(alpha_S, 2) / pi2;

  if (QT_finalstate_massive_coloured){
    //    calculate_Gamma1();
    if (initial_pdf_diag || initial_pdf_gq){
      // new: 20190126
      calculate_Gamma1born();
    }
  }

  calculate_A_F();
  //  A_F = calculate_A_F();
  // only with resummation ???
  calculate_B2_A_F();
  //  B2_A_F = calculate_B2_A_F();

  determine_splitting_tH1F_tH1(psi);
  ///  should be replaced by the following two contributions:
  ///  determine_splitting_tH1F_tH1_without_H1_delta(psi);
  //// -> coll_tH1F_contribution[no_xx/no_zx][...]
  ///  determine_splitting_tH1_only_H1_delta(psi);
  //// -> coll_tH1_only_H1_delta_contribution[no_xx][...]
  ////  new  coll_tH1  is accompanied by  coll_tH1_only_H1_delta  and does not contain respective terms any longer !!!

  determine_splitting_tgaga_tcga_tgamma2(psi);
  determine_splitting_tH2(psi);
  ///  contains H1 contributions that might have to be split...

  double LQ = log(psi_xbs_all[0][0] / pow(QT_Qres, 2));
  if (QT_Qres == 0){LQ = 0.;}

  for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
    for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
      value_LR[i_vr][i_mr] = log(psi_xbs_all[0][0] / pow(value_mu_ren[i_vr][i_mr], 2));
    }
  }

  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
      value_LF[i_vf][i_mf] = log(psi_xbs_all[0][0] / pow(value_mu_fact[i_vf][i_mf], 2));
      for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){

	if (i_vf == 0 && i_mf == 0){logger << LOG_DEBUG_POINT << "value_list_pdf_factor[0][0][" << i_l << "][0] = " << value_list_pdf_factor[i_vf][i_mf][i_l][0] << endl;}
	  
	coll_tH1F_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tH1F[i_l];
	coll_tH1_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tH1[i_l];
	//// to be introduced: additional term and coll_tH1 without H1 terms !!!
	////	coll_tH1_only_H1_delta_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tH1_only_H1_delta[i_l];

	coll_tgaga_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tgaga[i_l];
	coll_tcga_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tcga[i_l];
	coll_tgamma2_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tgamma2[i_l];
	coll_tH2_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tH2[i_l];
	//// modification in coll_tH2 needed ??? Split off H1 terms and replace by loopcc terms ???
      }
      value_tH1F[i_vf][i_mf] = accumulate(coll_tH1F_pdf[i_vf][i_mf].begin(), coll_tH1F_pdf[i_vf][i_mf].end(), 0.);
      value_tH1[i_vf][i_mf] = accumulate(coll_tH1_pdf[i_vf][i_mf].begin(), coll_tH1_pdf[i_vf][i_mf].end(), 0.);
      //// to be introduced: additional term and coll_tH1 without H1 terms !!!
      ////      value_tH1_only_H1_delta[i_vf][i_mf] = accumulate(coll_tH1_only_H1_delta_pdf[i_vf][i_mf].begin(), coll_tH1_only_H1_delta_pdf[i_vf][i_mf].end(), 0.);

      value_sig11[i_vf][i_mf] = calculate_sigma11(value_list_pdf_factor[i_vf][i_mf][0][0], value_tH1F[i_vf][i_mf], LQ);
      for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	  logger << LOG_DEBUG_POINT << "i_vf = " << i_vf << "   i_mf = " << i_mf << "   i_vr = " << i_vr << "   i_mr = " << i_mr << endl;
	  value_H1full[i_vr][i_vf][i_mr][i_mf] = calculate_H1(value_list_pdf_factor[i_vf][i_mf][0][0], value_tH1[i_vf][i_mf], value_tH1F[i_vf][i_mf], value_LR[i_vr][i_mr], value_LF[i_vf][i_mf], LQ);
	  //// to be introduced: New  calculate_H1  call based on tH1 and tH1_only_H1_delta !!!
	}
      }
      value_tgaga[i_vf][i_mf] = accumulate(coll_tgaga_pdf[i_vf][i_mf].begin(), coll_tgaga_pdf[i_vf][i_mf].end(), 0.);
      value_tcga[i_vf][i_mf] = accumulate(coll_tcga_pdf[i_vf][i_mf].begin(), coll_tcga_pdf[i_vf][i_mf].end(), 0.);
      value_tgamma2[i_vf][i_mf] = accumulate(coll_tgamma2_pdf[i_vf][i_mf].begin(), coll_tgamma2_pdf[i_vf][i_mf].end(), 0.);
      value_tH2[i_vf][i_mf] = accumulate(coll_tH2_pdf[i_vf][i_mf].begin(), coll_tH2_pdf[i_vf][i_mf].end(), 0.);
      // The following term does not exist in H2_gg, but in H2_qqbar !!! -> must be shifted to tH2(qqx) !!!
      /*
      value_tH2[i_vf][i_mf] += QT_H1_delta * (value_tH1[i_vf][i_mf] - QT_H1_delta * value_list_pdf_factor[i_vf][i_mf][0][0]);
      if (i_vf == 1 && i_mf == 1){
	logger.newLine(LOG_DEBUG);
	logger << LOG_DEBUG << "tH2 at [" << i_vf << "][" << i_mf << "]:" << endl;
	logger << LOG_DEBUG << "value_tH2[" << i_vf << "][" << i_mf << "] = " << accumulate(coll_tH2_pdf[i_vf][i_mf].begin(), coll_tH2_pdf[i_vf][i_mf].end(), 0.) << endl;
	logger << LOG_DEBUG << "tH2 extra term (involving H1²) = " << QT_H1_delta * (value_tH1[i_vf][i_mf] - QT_H1_delta * value_list_pdf_factor[i_vf][i_mf][0][0]) << endl;
      }
      */

      logger << LOG_DEBUG_VERBOSE << "switch_distribution = " << switch_distribution << endl;
      if (switch_distribution){
	for (int i_y = 0; i_y < 3; i_y++){
	  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
	    coll_tH1F_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tH1F[i_l];
	    coll_tH1_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tH1[i_l];
	    coll_tgaga_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tgaga[i_l];
	    coll_tcga_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tcga[i_l];
	    coll_tgamma2_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tgamma2[i_l];
	    coll_tH2_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tH2[i_l];
	  }
	  value_tH1F_D[i_vf][i_mf][i_y] = accumulate(coll_tH1F_pdf_D[i_vf][i_mf][i_y].begin(), coll_tH1F_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  value_tH1_D[i_vf][i_mf][i_y] = accumulate(coll_tH1_pdf_D[i_vf][i_mf][i_y].begin(), coll_tH1_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  value_sig11_D[i_vf][i_mf][i_y] = calculate_sigma11(value_list_pdf_factor[i_vf][i_mf][0][i_y], value_tH1F_D[i_vf][i_mf][i_y], LQ);
	  for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	    for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	      value_H1full_D[i_vr][i_vf][i_mr][i_mf][i_y] = calculate_H1(value_list_pdf_factor[i_vf][i_mf][0][i_y], value_tH1_D[i_vf][i_mf][i_y], value_tH1F_D[i_vf][i_mf][i_y], value_LR[i_vr][i_mr], value_LF[i_vf][i_mf], LQ);
	    }
	  }
	  value_tgaga_D[i_vf][i_mf][i_y] = accumulate(coll_tgaga_pdf_D[i_vf][i_mf][i_y].begin(), coll_tgaga_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  value_tcga_D[i_vf][i_mf][i_y] = accumulate(coll_tcga_pdf_D[i_vf][i_mf][i_y].begin(), coll_tcga_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  value_tgamma2_D[i_vf][i_mf][i_y] = accumulate(coll_tgamma2_pdf_D[i_vf][i_mf][i_y].begin(), coll_tgamma2_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  value_tH2_D[i_vf][i_mf][i_y] = accumulate(coll_tH2_pdf_D[i_vf][i_mf][i_y].begin(), coll_tH2_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  /*
	  // The following term does not exist in H2_gg, but in H2_qqbar !!! -> must be shifted to tH2(qqx) !!!
	  value_tH2_D[i_vf][i_mf][i_y] += QT_H1_delta * (value_tH1_D[i_vf][i_mf][i_y] - QT_H1_delta * value_list_pdf_factor[i_vf][i_mf][0][i_y]);
	  */
	}
      }
    }
  }

  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
      for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	  if (i_vr == 1 && i_mr == 1 && i_vf == 1 && i_mf == 1){
	    logger << LOG_DEBUG << "QT_H1_delta                 = " << QT_H1_delta << endl;
	    logger << LOG_DEBUG << "QT_H2_delta                 = " << QT_H2_delta << endl;
	    logger << LOG_DEBUG << "A_F                         = " << A_F << endl;
	    logger << LOG_DEBUG << "value_LR          [" << i_vr << "]   [" << i_mr << "] = " << value_LR[i_vr][i_mr] << endl;
	    logger << LOG_DEBUG << "value_LF       [" << i_vf << "]   [" << i_mf << "]    = " << value_LF[i_vf][i_mf] << endl;
	    logger << LOG_DEBUG << "value_sig11    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig11[i_vf][i_mf] << endl;
	    logger << LOG_DEBUG << "value_tH1      [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1[i_vf][i_mf] << endl;
	    logger << LOG_DEBUG << "value_tH1F     [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1F[i_vf][i_mf] << endl;
	    logger << LOG_DEBUG << "value_H1full   [" << i_vr << "][" << i_vf << "][" << i_mr << "][" << i_mf << "] = " << value_H1full[i_vr][i_vf][i_mr][i_mf] << endl;
	    logger << LOG_DEBUG << "value_tgaga    [" << i_vf << "]   [" << i_mf << "]    = " << value_tgaga[i_vf][i_mf] << endl;
	    logger << LOG_DEBUG << "value_tcga     [" << i_vf << "]   [" << i_mf << "]    = " << value_tcga[i_vf][i_mf] << endl;
	    logger << LOG_DEBUG << "value_tgamma2  [" << i_vf << "]   [" << i_mf << "]    = " << value_tgamma2[i_vf][i_mf] << endl;
	    logger << LOG_DEBUG << "value_tH2      [" << i_vf << "]   [" << i_mf << "]    = " << value_tH2[i_vf][i_mf] << endl;
	    logger.newLine(LOG_DEBUG);
	  }
	  logger << LOG_DEBUG_VERBOSE << "QT_H1_delta                 = " << QT_H1_delta << endl;
	  logger << LOG_DEBUG_VERBOSE << "QT_H2_delta                 = " << QT_H2_delta << endl;
	  logger << LOG_DEBUG_VERBOSE << "A_F                         = " << A_F << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_LR          [" << i_vr << "]   [" << i_mr << "] = " << value_LR[i_vr][i_mr] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_LF       [" << i_vf << "]   [" << i_mf << "]    = " << value_LF[i_vf][i_mf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_sig11    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig11[i_vf][i_mf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_tH1      [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1[i_vf][i_mf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_tH1F     [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1F[i_vf][i_mf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_H1full   [" << i_vr << "][" << i_vf << "][" << i_mr << "][" << i_mf << "] = " << value_H1full[i_vr][i_vf][i_mr][i_mf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_tgaga    [" << i_vf << "]   [" << i_mf << "]    = " << value_tgaga[i_vf][i_mf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_tcga     [" << i_vf << "]   [" << i_mf << "]    = " << value_tcga[i_vf][i_mf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_tgamma2  [" << i_vf << "]   [" << i_mf << "]    = " << value_tgamma2[i_vf][i_mf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_tH2      [" << i_vf << "]   [" << i_mf << "]    = " << value_tH2[i_vf][i_mf] << endl;
	  logger.newLine(LOG_DEBUG_VERBOSE);
	}
      }
    }
  }



  if (switch_resummation){
    // !!!
#ifdef MORE
    double Q = sqrt(psi_xbs_all[0][0]);
    double y = log(sqrt(psi_x_pdf[1] / psi_x_pdf[2]));
    double muF = var_mu_fact;
    double muRen = var_mu_ren;
    int channel = 0;
    double sqrt_shad = 2 * psi_E;
    double qt = sqrt(psi_QT_qt2);
    //  performQTBoost(psi_QT_qt2, Q, y, particle_event);
    logger << LOG_DEBUG << "Q=" << Q << ", " << "y=" << y << ", " << "qt=" << qt << ", x1=" << x_pdf[1] << ", x2=" << x_pdf[2] << endl;

    double temp = 0.;
    if (abs(y) > -0.5 * log(pow(somescale_.startscale / sqrt_shad, 2))){cout << "y too large, skipping point" << endl;}
    else {
      resumm_(&temp, &gg_initiated, &channel, &Q, &y, &qt, &QT_Qres, &muF, &muRen, &sqrt_shad, &VA_b_ME2, &QT_H1_delta, &QT_H2_delta);
      if (munich_isnan(temp)){temp = 0.; cout << "CX_value_integrand_RF_qTcut[0][i_vr][i_vf][i_mr][i_mf] is nan, setting to zero!" << endl;}
    }
    for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
      for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
	  for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
	    CX_value_integrand_RF_qTcut[0][i_vr][i_vf][i_mr][i_mf] = temp;
	  }
	}
      }
    }

    integrand = psi_ps_factor * rescaling_factor_alpha_e * CX_value_integrand_RF_qTcut[0][dynamic_scale][dynamic_scale][map_value_scale_ren][map_value_scale_fact] / x_pdf[0] * psi_QT_jacqt2 / 2 / qt;
    if (switch_CV){
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	integrand_CV[i_s] = psi_ps_factor * rescaling_factor_alpha_e * CX_value_integrand_RF_qTcut[0][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]] / x_pdf[0] * psi_QT_jacqt2 / 2 / qt;
      }
    }

    if (switch_distribution){
      for (int i_y = 0; i_y < 3; i_y++){
        integrand_D[i_y][0] = psi_ps_factor * rescaling_factor_alpha_e * CX_value_integrand_RF_qTcut[0][dynamic_scale][dynamic_scale][map_value_scale_ren][map_value_scale_fact] / x_pdf[0] * psi_QT_jacqt2 / 2 / qt;
	if (i_y > 0){integrand_D[i_y][0] /= 2;}  //         TODO: asymmetric distributions not implemented !!! (refers to issues pointed out before...)
 
        if (switch_CV){
          for (int i_s = 0; i_s < n_scales_CV; i_s++){
            integrand_D_CV[i_s][i_y][0] = psi_ps_factor * rescaling_factor_alpha_e * CX_value_integrand_RF_qTcut[0][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]] / x_pdf[0] * psi_QT_jacqt2 / 2 / qt;
	    if (i_y > 0){integrand_D_CV[i_s][i_y][0] /= 2;}  //         TODO: asymmetric distributions not implemented !!! (refers to issues pointed out before...)
         }
        }
      }
    }
#endif
  }
  else {
    for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
      for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
	for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	  for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	      CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf][i_mr][i_mf] = global_factor * calculate_H2(value_list_pdf_factor[i_vf][i_mf][0][0], value_tH1[i_vf][i_mf], value_tH1F[i_vf][i_mf], value_H1full[i_vr][i_vf][i_mr][i_mf], value_sig11[i_vf][i_mf], value_tgaga[i_vf][i_mf], value_tcga[i_vf][i_mf], value_tgamma2[i_vf][i_mf], value_tH2[i_vf][i_mf], value_LR[i_vr][i_mr], value_LF[i_vf][i_mf], LQ);
	      if (switch_distribution){
		for (int i_y = 0; i_y < 3; i_y++){
		  CX_value_integrand_RF_qTcut_D[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = global_factor * calculate_H2(value_list_pdf_factor[i_vf][i_mf][0][i_y], value_tH1_D[i_vf][i_mf][i_y], value_tH1F_D[i_vf][i_mf][i_y], value_H1full_D[i_vr][i_vf][i_mr][i_mf][i_y], value_sig11_D[i_vf][i_mf][i_y], value_tgaga_D[i_vf][i_mf][i_y], value_tcga_D[i_vf][i_mf][i_y], value_tgamma2_D[i_vf][i_mf][i_y], value_tH2_D[i_vf][i_mf][i_y], value_LR[i_vr][i_mr], value_LF[i_vf][i_mf], LQ);
		}
	      }
	    }
	  }
	}
      }
    }

    int x_q = 0;
    integrand = var_rel_alpha_S * psi_ps_factor * rescaling_factor_alpha_e * VA_b_ME2 * CX_value_integrand_RF_qTcut[x_q][dynamic_scale][dynamic_scale][map_value_scale_ren][map_value_scale_fact];
    if (switch_CV){
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	integrand_CV[i_s] = var_rel_alpha_S_CV[i_s] * psi_ps_factor * rescaling_factor_alpha_e * VA_b_ME2 * CX_value_integrand_RF_qTcut[x_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]];
	for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	  integrand_qTcut_CV[i_q][i_s] = var_rel_alpha_S_CV[i_s] * psi_ps_factor * rescaling_factor_alpha_e * VA_b_ME2 * CX_value_integrand_RF_qTcut[i_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]];
	}
      }
    }
    
    if (switch_distribution){
      for (int i_y = 0; i_y < 3; i_y++){
	integrand_D[i_y][0] = var_rel_alpha_S * psi_ps_factor * rescaling_factor_alpha_e * VA_b_ME2 * CX_value_integrand_RF_qTcut_D[x_q][dynamic_scale][dynamic_scale][map_value_scale_ren][map_value_scale_fact][i_y];
	if (switch_CV){
	  for (int i_s = 0; i_s < n_scales_CV; i_s++){
	    integrand_D_CV[i_s][i_y][0] = var_rel_alpha_S_CV[i_s] * psi_ps_factor * rescaling_factor_alpha_e * VA_b_ME2 * CX_value_integrand_RF_qTcut_D[x_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]][i_y];
	    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	      integrand_D_qTcut_CV[i_q][i_s][i_y][0] = var_rel_alpha_S_CV[i_s] * psi_ps_factor * rescaling_factor_alpha_e * VA_b_ME2 * CX_value_integrand_RF_qTcut_D[i_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]][i_y];
	    }
	  }
	}
      }
    }
  }



  if (switch_TSV){
    for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
      for (int i_mr = 0; i_mr < n_scale_dyn_ren[i_vr]; i_mr++){
	value_LR_TSV[i_vr][i_mr] = log(psi_xbs_all[0][0] / pow(value_scale_ren[0][i_vr][i_mr], 2));
      }
    }
    for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
      for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
	value_LF_TSV[i_vf][i_mf] = log(psi_xbs_all[0][0] / pow(value_scale_fact[0][i_vf][i_mf], 2));
	for (int i_y = 0; i_y < 3; i_y++){
	  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
	    if (i_vf == 0 && i_mf == 0 && i_y == 0){
	      logger << LOG_DEBUG << "value_list_pdf_factor_TSV[0][0][" << i_l << "][0] = " << value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][0] << endl;
	      logger << LOG_DEBUG << "   coll_tH1F[" << i_l << "] = " << coll_tH1F[i_l] << endl;
	      logger << LOG_DEBUG << "    coll_tH1[" << i_l << "] = " << coll_tH1[i_l] << endl;
	      logger << LOG_DEBUG << "  coll_tgaga[" << i_l << "] = " << coll_tgaga[i_l] << endl;
	      logger << LOG_DEBUG << "   coll_tcga[" << i_l << "] = " << coll_tcga[i_l] << endl;
	      logger << LOG_DEBUG << "coll_tgamma2[" << i_l << "] = " << coll_tgamma2[i_l] << endl;
	      logger << LOG_DEBUG << "    coll_tH2[" << i_l << "] = " << coll_tH2[i_l] << endl;
	    }
	  
	    coll_tH1F_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tH1F[i_l];
	    coll_tH1_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tH1[i_l];
	    coll_tgaga_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tgaga[i_l];
	    coll_tcga_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tcga[i_l];
	    coll_tgamma2_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tgamma2[i_l];
	    coll_tH2_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tH2[i_l];
  
	    if (i_vf == 0 && i_mf == 0 && i_y == 0){
	      logger << LOG_DEBUG_POINT << "value_list_pdf_factor_TSV[0][0][" << i_l << "][0] = " << value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][0] << endl;
	      logger << LOG_DEBUG_POINT << "   coll_tH1F_pdf_TSV[" << i_vf << "][" << i_mf << "][" << i_y << "][" << i_l << "] = " << coll_tH1F_pdf_TSV[i_vf][i_mf][i_y][i_l] << endl;
	      logger << LOG_DEBUG_POINT << "    coll_tH1_pdf_TSV[" << i_vf << "][" << i_mf << "][" << i_y << "][" << i_l << "] = " << coll_tH1_pdf_TSV[i_vf][i_mf][i_y][i_l] << endl;
	      logger << LOG_DEBUG_POINT << "  coll_tgaga_pdf_TSV[" << i_vf << "][" << i_mf << "][" << i_y << "][" << i_l << "] = " << coll_tgaga_pdf_TSV[i_vf][i_mf][i_y][i_l] << endl;
	      logger << LOG_DEBUG_POINT << "   coll_tcga_pdf_TSV[" << i_vf << "][" << i_mf << "][" << i_y << "][" << i_l << "] = " << coll_tcga_pdf_TSV[i_vf][i_mf][i_y][i_l] << endl;
	      logger << LOG_DEBUG_POINT << "coll_tgamma2_pdf_TSV[" << i_vf << "][" << i_mf << "][" << i_y << "][" << i_l << "] = " << coll_tgamma2_pdf_TSV[i_vf][i_mf][i_y][i_l] << endl;
	      logger << LOG_DEBUG_POINT << "    coll_tH2_pdf_TSV[" << i_vf << "][" << i_mf << "][" << i_y << "][" << i_l << "] = " << coll_tH2_pdf_TSV[i_vf][i_mf][i_y][i_l] << endl;
	    }
	  }
	  value_tH1F_TSV[i_vf][i_mf][i_y] = accumulate(coll_tH1F_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tH1F_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
	  value_tH1_TSV[i_vf][i_mf][i_y] = accumulate(coll_tH1_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tH1_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
	  value_sig11_TSV[i_vf][i_mf][i_y] = calculate_sigma11(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y], value_tH1F_TSV[i_vf][i_mf][i_y], LQ);
	  for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
	    for (int i_mr = 0; i_mr < n_scale_dyn_ren[i_vr]; i_mr++){
	      logger << LOG_DEBUG_POINT << "TSV   i_y = " << i_y << "   i_vf = " << i_vf << "   i_mf = " << i_mf << "   i_vr = " << i_vr << "   i_mr = " << i_mr << endl;
	      value_H1full_TSV[i_vr][i_vf][i_mr][i_mf][i_y] = calculate_H1(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y], value_tH1_TSV[i_vf][i_mf][i_y], value_tH1F_TSV[i_vf][i_mf][i_y], value_LR_TSV[i_vr][i_mr], value_LF_TSV[i_vf][i_mf], LQ);
	    }
	  }
	  value_tgaga_TSV[i_vf][i_mf][i_y] = accumulate(coll_tgaga_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tgaga_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
	  value_tcga_TSV[i_vf][i_mf][i_y] = accumulate(coll_tcga_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tcga_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
	  value_tgamma2_TSV[i_vf][i_mf][i_y] = accumulate(coll_tgamma2_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tgamma2_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
	  value_tH2_TSV[i_vf][i_mf][i_y] = accumulate(coll_tH2_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tH2_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);

	  /*
	  // The following term does not exist in H2_gg, but in H2_qqbar !!! -> must be shifted to tH2(qqx) !!!
	  value_tH2_TSV[i_vf][i_mf][i_y] += QT_H1_delta * (value_tH1_TSV[i_vf][i_mf][i_y] - QT_H1_delta * value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y]);
	  if (i_vf == 1 && i_mf == 1){
	    logger.newLine(LOG_DEBUG);
	    logger << LOG_DEBUG << "tH2_TSV at [" << i_vf << "][" << i_mf << "][" << i_y << "]:" << endl;
	    logger << LOG_DEBUG << "value_tH2_TSV[" << i_vf << "][" << i_mf << "][" << i_y << "] = " << accumulate(coll_tH2_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tH2_pdf_TSV[i_vf][i_mf][i_y].end(), 0.) << endl;
	    logger << LOG_DEBUG << "tH2 extra term (involving H1²) = " << QT_H1_delta * (value_tH1_TSV[i_vf][i_mf][i_y] - QT_H1_delta * value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y]) << endl;
	  }
	  */
	}
      }
    }

    for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
      for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
	for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
	  for (int i_mr = 0; i_mr < n_scale_dyn_ren[i_vr]; i_mr++){
	    if (i_vr == 1 && i_mr == 1 && i_vf == 1 && i_mf == 1){
	      logger << LOG_DEBUG << "QT_H1_delta                     = " << QT_H1_delta << endl;
	      logger << LOG_DEBUG << "QT_H2_delta                     = " << QT_H2_delta << endl;
	      logger << LOG_DEBUG << "A_F                             = " << A_F << endl;
	      logger << LOG_DEBUG << "value_LR_TSV          [" << i_vr << "]   [" << i_mr << "] = " << value_LR_TSV[i_vr][i_mr] << endl;
	      logger << LOG_DEBUG << "value_LF_TSV       [" << i_vf << "]   [" << i_mf << "]    = " << value_LF_TSV[i_vf][i_mf] << endl;
	      logger << LOG_DEBUG << "value_sig11_TSV    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig11_TSV[i_vf][i_mf][0] << endl;
	      logger << LOG_DEBUG << "value_tH1_TSV      [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1_TSV[i_vf][i_mf][0] << endl;
	      logger << LOG_DEBUG << "value_tH1F_TSV     [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1F_TSV[i_vf][i_mf][0] << endl;
	      logger << LOG_DEBUG << "value_H1full_TSV   [" << i_vr << "][" << i_vf << "][" << i_mr << "][" << i_mf << "] = " << value_H1full_TSV[i_vr][i_vf][i_mr][i_mf][0] << endl;
	      logger << LOG_DEBUG << "value_tgaga_TSV    [" << i_vf << "]   [" << i_mf << "]    = " << value_tgaga_TSV[i_vf][i_mf][0] << endl;
	      logger << LOG_DEBUG << "value_tcga_TSV     [" << i_vf << "]   [" << i_mf << "]    = " << value_tcga_TSV[i_vf][i_mf][0] << endl;
	      logger << LOG_DEBUG << "value_tgamma2_TSV  [" << i_vf << "]   [" << i_mf << "]    = " << value_tgamma2_TSV[i_vf][i_mf][0] << endl;
	      logger << LOG_DEBUG << "value_tH2_TSV      [" << i_vf << "]   [" << i_mf << "]    = " << value_tH2_TSV[i_vf][i_mf][0] << endl;
	      logger.newLine(LOG_DEBUG);
	    }
	    logger << LOG_DEBUG_VERBOSE << "QT_H1_delta                     = " << QT_H1_delta << endl;
	    logger << LOG_DEBUG_VERBOSE << "QT_H2_delta                     = " << QT_H2_delta << endl;
	    logger << LOG_DEBUG_VERBOSE << "A_F                             = " << A_F << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_LR_TSV          [" << i_vr << "]   [" << i_mr << "] = " << value_LR_TSV[i_vr][i_mr] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_LF_TSV       [" << i_vf << "]   [" << i_mf << "]    = " << value_LF_TSV[i_vf][i_mf] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_sig11_TSV    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig11_TSV[i_vf][i_mf][0] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_tH1_TSV      [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1_TSV[i_vf][i_mf][0] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_tH1F_TSV     [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1F_TSV[i_vf][i_mf][0] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_H1full_TSV   [" << i_vr << "][" << i_vf << "][" << i_mr << "][" << i_mf << "] = " << value_H1full_TSV[i_vr][i_vf][i_mr][i_mf][0] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_tgaga_TSV    [" << i_vf << "]   [" << i_mf << "]    = " << value_tgaga_TSV[i_vf][i_mf][0] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_tcga_TSV     [" << i_vf << "]   [" << i_mf << "]    = " << value_tcga_TSV[i_vf][i_mf][0] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_tgamma2_TSV  [" << i_vf << "]   [" << i_mf << "]    = " << value_tgamma2_TSV[i_vf][i_mf][0] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_tH2_TSV      [" << i_vf << "]   [" << i_mf << "]    = " << value_tH2_TSV[i_vf][i_mf][0] << endl;
	    logger.newLine(LOG_DEBUG_VERBOSE);
   	  }
	}
      }
    }


    if (switch_resummation){
#ifdef MORE
      double Q = sqrt(psi_xbs_all[0][0]);
      double y = log(sqrt(psi_x_pdf[1] / psi_x_pdf[2]));
      double muF = var_mu_fact;
      double muRen = var_mu_ren;
      //      int gg_initiated;
      int channel = 0;
      double sqrt_shad = 2 * psi_E;
      double qt = sqrt(psi_QT_qt2);
        
      double temp = 0.;
      if (abs(y) > -0.5 * log(pow(somescale_.startscale / sqrt_shad, 2))){cout << "y too large, skipping point" << endl;}
      else {
	resumm_(&temp, &gg_initiated, &channel, &Q, &y, &qt, &QT_Qres, &muF, &muRen, &sqrt_shad, &VA_b_ME2, &QT_H1_delta, &QT_H2_delta);
	if (munich_isnan(temp)){temp = 0.; cout << "value_integrand_qTcut_TSV[i_q][i_vr][i_vf][i_mr][i_mf][i_y] is nan, setting to zero!" << endl;}
      }
      for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
	for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
	  for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
	    for (int i_mr = 0; i_mr < n_scale_dyn_ren[i_vr]; i_mr++){
	      for (int i_y = 0; i_y < 3; i_y++){
		// check i_q dependence !!!
		for (int i_q = 0; i_q < output_n_qTcut; i_q++){
		  // Most likely, psi_ps_factor * rescaling_factor_alpha_e needs to be removed !!! What about global_factor, or VA_b_ME2?
		  value_integrand_qTcut_TSV[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = psi_ps_factor * rescaling_factor_alpha_e * temp / x_pdf[0] * psi_QT_jacqt2 / 2 / qt;
		}
	      }
	    }
	  }
	}
      }
#endif
    }
    else {
      for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
	for (int i_mr = 0; i_mr < n_scale_dyn_ren[i_vr]; i_mr++){
	  for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
	    for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
	      for (int i_y = 0; i_y < 3; i_y++){
		for (int i_q = 0; i_q < output_n_qTcut; i_q++){
		  logger << LOG_DEBUG_VERBOSE << "H2_TSV at [" << i_vr << "][" << i_vf << "][" << i_mr << "][" << i_mf << "][" << i_y << "][" << i_q << "]:" << endl;
		  value_integrand_qTcut_TSV[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = global_factor * VA_b_ME2 * calculate_H2(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y], value_tH1_TSV[i_vf][i_mf][i_y], value_tH1F_TSV[i_vf][i_mf][i_y], value_H1full_TSV[i_vr][i_vf][i_mr][i_mf][i_y], value_sig11_TSV[i_vf][i_mf][i_y], value_tgaga_TSV[i_vf][i_mf][i_y], value_tcga_TSV[i_vf][i_mf][i_y], value_tgamma2_TSV[i_vf][i_mf][i_y], value_tH2_TSV[i_vf][i_mf][i_y], value_LR_TSV[i_vr][i_mr], value_LF_TSV[i_vf][i_mf], LQ);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
