#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

void observable_set::determine_integrand_CX_ncollinear_CT2(phasespace_set & psi){
  Logger logger("observable_set::determine_integrand_CX_ncollinear_CT2");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  static int order = 2;

  double global_factor = pow(alpha_S, 2) / pi2;
  double qtoverq = sqrt(psi_QT_qt2 / QT_Qres / QT_Qres);
  double LL1, LL2, LL3, LL4;
  if (switch_resummation){Itilde(qtoverq, order, LL1, LL2, LL3, LL4);}
  // in FO computations, we do the qT integrals once and for all in the beginning
  // in resummed computations, we have to bin the "artificial" qT, so that is not easily possible -> revert back to old implementation  

  if (QT_finalstate_massive_coloured){calculate_Gamma1();}

  calculate_A_F();
  calculate_B2_A_F();
  //  A_F = calculate_A_F();
  determine_splitting_tH1F_tH1(psi);
  determine_splitting_tgaga_tcga_tgamma2(psi);

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
	//	coll_tcga[i_l] = 0.; // !!!
	coll_tH1F_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tH1F[i_l];
	coll_tH1_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tH1[i_l];
	coll_tgaga_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tgaga[i_l];
	coll_tcga_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tcga[i_l];
	coll_tgamma2_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tgamma2[i_l];
      }
      value_tH1F[i_vf][i_mf] = accumulate(coll_tH1F_pdf[i_vf][i_mf].begin(), coll_tH1F_pdf[i_vf][i_mf].end(), 0.);
      value_tH1[i_vf][i_mf] = accumulate(coll_tH1_pdf[i_vf][i_mf].begin(), coll_tH1_pdf[i_vf][i_mf].end(), 0.);
      value_tgaga[i_vf][i_mf] = accumulate(coll_tgaga_pdf[i_vf][i_mf].begin(), coll_tgaga_pdf[i_vf][i_mf].end(), 0.);
      value_tcga[i_vf][i_mf] = accumulate(coll_tcga_pdf[i_vf][i_mf].begin(), coll_tcga_pdf[i_vf][i_mf].end(), 0.);
      value_tgamma2[i_vf][i_mf] = accumulate(coll_tgamma2_pdf[i_vf][i_mf].begin(), coll_tgamma2_pdf[i_vf][i_mf].end(), 0.);
      value_sig11[i_vf][i_mf] = calculate_sigma11(value_list_pdf_factor[i_vf][i_mf][0][0], value_tH1F[i_vf][i_mf], LQ);
      value_sig24[i_vf][i_mf] = calculate_sigma24(value_list_pdf_factor[i_vf][i_mf][0][0]);
      value_sig23[i_vf][i_mf] = calculate_sigma23(value_list_pdf_factor[i_vf][i_mf][0][0], value_sig11[i_vf][i_mf]);

      for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	  value_H1full[i_vr][i_vf][i_mr][i_mf] = calculate_H1(value_list_pdf_factor[i_vf][i_mf][0][0], value_tH1[i_vf][i_mf], value_tH1F[i_vf][i_mf], value_LR[i_vr][i_mr], value_LF[i_vf][i_mf], LQ);
	  value_sig22[i_vr][i_vf][i_mr][i_mf] = calculate_sigma22(value_list_pdf_factor[i_vf][i_mf][0][0], value_sig11[i_vf][i_mf], value_tH1[i_vf][i_mf], value_tH1F[i_vf][i_mf], value_H1full[i_vr][i_vf][i_mr][i_mf], value_tgaga[i_vf][i_mf], value_LR[i_vr][i_mr], value_LF[i_vf][i_mf], LQ);
	  value_sig21[i_vr][i_vf][i_mr][i_mf] = calculate_sigma21(value_list_pdf_factor[i_vf][i_mf][0][0], value_sig11[i_vf][i_mf], value_tH1[i_vf][i_mf], value_tH1F[i_vf][i_mf], value_H1full[i_vr][i_vf][i_mr][i_mf], value_tgaga[i_vf][i_mf], value_tcga[i_vf][i_mf], value_tgamma2[i_vf][i_mf], value_LR[i_vr][i_mr], value_LF[i_vf][i_mf], LQ, A_F);
	}
      }

      if (switch_distribution){
	for (int i_y = 0; i_y < 3; i_y++){
	  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
	    coll_tH1F_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tH1F[i_l];
	    coll_tH1_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tH1[i_l];
	    coll_tgaga_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tgaga[i_l];
	    coll_tcga_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tcga[i_l];
	    coll_tgamma2_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tgamma2[i_l];
	  }
	  value_tH1F_D[i_vf][i_mf][i_y] = accumulate(coll_tH1F_pdf_D[i_vf][i_mf][i_y].begin(), coll_tH1F_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  value_tH1_D[i_vf][i_mf][i_y] = accumulate(coll_tH1_pdf_D[i_vf][i_mf][i_y].begin(), coll_tH1_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  value_tgaga_D[i_vf][i_mf][i_y] = accumulate(coll_tgaga_pdf_D[i_vf][i_mf][i_y].begin(), coll_tgaga_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  value_tcga_D[i_vf][i_mf][i_y] = accumulate(coll_tcga_pdf_D[i_vf][i_mf][i_y].begin(), coll_tcga_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  value_tgamma2_D[i_vf][i_mf][i_y] = accumulate(coll_tgamma2_pdf_D[i_vf][i_mf][i_y].begin(), coll_tgamma2_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  value_sig11_D[i_vf][i_mf][i_y] = calculate_sigma11(value_list_pdf_factor[i_vf][i_mf][0][i_y], value_tH1F_D[i_vf][i_mf][i_y], LQ);
	  value_sig24_D[i_vf][i_mf][i_y] = calculate_sigma24(value_list_pdf_factor[i_vf][i_mf][0][i_y]);
	  value_sig23_D[i_vf][i_mf][i_y] = calculate_sigma23(value_list_pdf_factor[i_vf][i_mf][0][i_y], value_sig11_D[i_vf][i_mf][i_y]);

	  for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	    for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	      value_H1full_D[i_vr][i_vf][i_mr][i_mf][i_y] = calculate_H1(value_list_pdf_factor[i_vf][i_mf][0][i_y], value_tH1_D[i_vf][i_mf][i_y], value_tH1F_D[i_vf][i_mf][i_y], value_LR[i_vr][i_mr], value_LF[i_vf][i_mf], LQ);
	      value_sig22_D[i_vr][i_vf][i_mr][i_mf][i_y] = calculate_sigma22(value_list_pdf_factor[i_vf][i_mf][0][i_y], value_sig11_D[i_vf][i_mf][i_y], value_tH1_D[i_vf][i_mf][i_y], value_tH1F_D[i_vf][i_mf][i_y], value_H1full_D[i_vr][i_vf][i_mr][i_mf][i_y], value_tgaga_D[i_vf][i_mf][i_y], value_LR[i_vr][i_mr], value_LF[i_vf][i_mf], LQ);
	      value_sig21_D[i_vr][i_vf][i_mr][i_mf][i_y] = calculate_sigma21(value_list_pdf_factor[i_vf][i_mf][0][i_y], value_sig11_D[i_vf][i_mf][i_y], value_tH1_D[i_vf][i_mf][i_y], value_tH1F_D[i_vf][i_mf][i_y], value_H1full_D[i_vr][i_vf][i_mr][i_mf][i_y], value_tgaga_D[i_vf][i_mf][i_y], value_tcga_D[i_vf][i_mf][i_y], value_tgamma2_D[i_vf][i_mf][i_y], value_LR[i_vr][i_mr], value_LF[i_vf][i_mf], LQ, A_F);
	    }
	  }
	}
      }
    }
  }
      
  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
      for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	  logger << LOG_DEBUG_VERBOSE << "QT_H1_delta                 = " << QT_H1_delta << endl;
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
	  logger << LOG_DEBUG_VERBOSE << "value_sig21    [" << i_vr << "][" << i_vf << "][" << i_mr << "][" << i_mf << "] = " << value_sig21[i_vr][i_vf][i_mr][i_mf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_sig22    [" << i_vr << "][" << i_vf << "][" << i_mr << "][" << i_mf << "] = " << value_sig22[i_vr][i_vf][i_mr][i_mf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_sig23    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig23[i_vf][i_mf] << endl;
	  logger << LOG_DEBUG_VERBOSE << "value_sig24    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig24[i_vf][i_mf] << endl;
	  logger.newLine(LOG_DEBUG_VERBOSE);
	}
      }
    }
  }

  if (switch_resummation){
    for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
      for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
	for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	  for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	    for (int i_q = 0; i_q < n_qTcut; i_q++){
	      if (i_q <= cut_ps[0]){
		CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf][i_mr][i_mf] = global_factor 
		  * (value_sig21[i_vr][i_vf][i_mr][i_mf] * LL1 + 
		     value_sig22[i_vr][i_vf][i_mr][i_mf] * LL2 + 
		     value_sig23[i_vf][i_mf] * LL3 + 
		     value_sig24[i_vf][i_mf] * LL4) / pow(QT_Qres,2) * psi_QT_jacqt2;
		if (switch_distribution){
		  for (int i_y = 0; i_y < 3; i_y++){
		    CX_value_integrand_RF_qTcut_D[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = global_factor 
		      * (value_sig21_D[i_vr][i_vf][i_mr][i_mf][i_y] * LL1 + 
			 value_sig22_D[i_vr][i_vf][i_mr][i_mf][i_y] * LL2 + 
			 value_sig23_D[i_vf][i_mf][i_y] * LL3 + 
			 value_sig24_D[i_vf][i_mf][i_y] * LL4) / pow(QT_Qres,2) * psi_QT_jacqt2;
		  }
		}
	      }
	      else {
		CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf][i_mr][i_mf] = 0.;
		if (switch_distribution){
		  for (int i_y = 0; i_y < 3; i_y++){
		    CX_value_integrand_RF_qTcut_D[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = 0.;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  else {
    for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
      for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
	for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	  for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	    for (int i_q = 0; i_q < n_qTcut; i_q++){
	      CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf][i_mr][i_mf] = global_factor 
		* (value_sig21[i_vr][i_vf][i_mr][i_mf] * I1_int[i_q] + 
		   value_sig22[i_vr][i_vf][i_mr][i_mf] * I2_int[i_q] + 
		   value_sig23[i_vf][i_mf] * I3_int[i_q] + 
		   value_sig24[i_vf][i_mf] * I4_int[i_q]);
	      if (switch_distribution){
		for (int i_y = 0; i_y < 3; i_y++){
		  CX_value_integrand_RF_qTcut_D[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = global_factor 
		    * (value_sig21_D[i_vr][i_vf][i_mr][i_mf][i_y] * I1_int[i_q] + 
		       value_sig22_D[i_vr][i_vf][i_mr][i_mf][i_y] * I2_int[i_q] + 
		       value_sig23_D[i_vf][i_mf][i_y] * I3_int[i_q] + 
		       value_sig24_D[i_vf][i_mf][i_y] * I4_int[i_q]);
		}
	      }
	    }
	  }
        }
      }
    }
  }
  /*
  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
      for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	for (int i_mr = 0; i_mr < value_mu_ren[i_vr].size(); i_mr++){
	  for (int i_q = 0; i_q < n_qTcut; i_q++){
	    if (switch_resummation){
	      if (i_q <= cut_ps[0]){
		CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf][i_mr][i_mf] = global_factor 
		  * (value_sig21[i_vr][i_vf][i_mr][i_mf] * LL1 + 
		     value_sig22[i_vr][i_vf][i_mr][i_mf] * LL2 + 
		     value_sig23[i_vf][i_mf] * LL3 + 
		     value_sig24[i_vf][i_mf] * LL4) / pow(QT_Qres,2) * psi_QT_jacqt2;
		if (switch_distribution){
		  for (int i_y = 0; i_y < 3; i_y++){
		    CX_value_integrand_RF_qTcut_D[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = global_factor 
		      * (value_sig21_D[i_vr][i_vf][i_mr][i_mf][i_y] * LL1 + 
			 value_sig22_D[i_vr][i_vf][i_mr][i_mf][i_y] * LL2 + 
			 value_sig23_D[i_vf][i_mf][i_y] * LL3 + 
			 value_sig24_D[i_vf][i_mf][i_y] * LL4) / pow(QT_Qres,2) * psi_QT_jacqt2;
		  }
		}
	      }
	      else {
		CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf][i_mr][i_mf] = 0.;
		if (switch_distribution){
		  for (int i_y = 0; i_y < 3; i_y++){
		    CX_value_integrand_RF_qTcut_D[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = 0.;
		  }
		}
	      }
	    }

	    else {
	      CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf][i_mr][i_mf] = global_factor 
		* (value_sig21[i_vr][i_vf][i_mr][i_mf] * I1_int[i_q] + 
		   value_sig22[i_vr][i_vf][i_mr][i_mf] * I2_int[i_q] + 
		   value_sig23[i_vf][i_mf] * I3_int[i_q] + 
		   value_sig24[i_vf][i_mf] * I4_int[i_q]);
	      if (switch_distribution){
		for (int i_y = 0; i_y < 3; i_y++){
		  CX_value_integrand_RF_qTcut_D[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = global_factor 
		    * (value_sig21_D[i_vr][i_vf][i_mr][i_mf][i_y] * I1_int[i_q] + 
		       value_sig22_D[i_vr][i_vf][i_mr][i_mf][i_y] * I2_int[i_q] + 
		       value_sig23_D[i_vf][i_mf][i_y] * I3_int[i_q] + 
		       value_sig24_D[i_vf][i_mf][i_y] * I4_int[i_q]);
		}
	      }
	    }
	  }
        }
      }
    }
  }
  */

  integrand = -var_rel_alpha_S * psi_ps_factor * ME2 * CX_value_integrand_RF_qTcut[0][dynamic_scale][dynamic_scale][map_value_scale_ren][map_value_scale_fact];
  if (switch_CV){
    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      integrand_CV[i_s] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * CX_value_integrand_RF_qTcut[0][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]];
      for (int i_q = 0; i_q < n_qTcut; i_q++){
	integrand_qTcut_CV[i_q][i_s] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * CX_value_integrand_RF_qTcut[i_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]];
      }
    }
  }
  
  if (switch_distribution){
    int x_q = 0;
    for (int j = 0; j < 3; j++){
      integrand_D[j][0] = -var_rel_alpha_S * psi_ps_factor * ME2 * CX_value_integrand_RF_qTcut_D[x_q][dynamic_scale][dynamic_scale][map_value_scale_ren][map_value_scale_fact][j];
      if (switch_CV){
	for (int i_s = 0; i_s < n_scales_CV; i_s++){
	  integrand_D_CV[i_s][j][0] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * CX_value_integrand_RF_qTcut_D[x_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]][j];
	  for (int i_q = 0; i_q < n_qTcut; i_q++){
	    integrand_D_qTcut_CV[i_q][i_s][j][0] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * CX_value_integrand_RF_qTcut_D[i_q][dynamic_scale_CV][dynamic_scale_CV][map_value_scale_ren_CV[i_s]][map_value_scale_fact_CV[i_s]][j];
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
	    coll_tH1F_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tH1F[i_l];
	    coll_tH1_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tH1[i_l];
	    coll_tgaga_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tgaga[i_l];
	    coll_tcga_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tcga[i_l];
	    coll_tgamma2_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tgamma2[i_l];
	    coll_tH2_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tH2[i_l];
	  }
	  value_tH1F_TSV[i_vf][i_mf][i_y] = accumulate(coll_tH1F_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tH1F_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
	  value_tH1_TSV[i_vf][i_mf][i_y] = accumulate(coll_tH1_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tH1_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
	  value_tgaga_TSV[i_vf][i_mf][i_y] = accumulate(coll_tgaga_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tgaga_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
	  value_tcga_TSV[i_vf][i_mf][i_y] = accumulate(coll_tcga_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tcga_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
	  value_tgamma2_TSV[i_vf][i_mf][i_y] = accumulate(coll_tgamma2_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tgamma2_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);	
	  value_sig11_TSV[i_vf][i_mf][i_y] = calculate_sigma11(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y], value_tH1F_TSV[i_vf][i_mf][i_y], LQ);
	  value_sig23_TSV[i_vf][i_mf][i_y] = calculate_sigma23(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y], value_sig11_TSV[i_vf][i_mf][i_y]);
	  value_sig24_TSV[i_vf][i_mf][i_y] = calculate_sigma24(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y]);
	  for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
	    for (int i_mr = 0; i_mr < n_scale_dyn_ren[i_vr]; i_mr++){
	      value_H1full_TSV[i_vr][i_vf][i_mr][i_mf][i_y] = calculate_H1(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y], value_tH1_TSV[i_vf][i_mf][i_y], value_tH1F_TSV[i_vf][i_mf][i_y], value_LR_TSV[i_vr][i_mr], value_LF_TSV[i_vf][i_mf], LQ);
	      value_sig22_TSV[i_vr][i_vf][i_mr][i_mf][i_y] = calculate_sigma22(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y], value_sig11_TSV[i_vf][i_mf][i_y], value_tH1_TSV[i_vf][i_mf][i_y], value_tH1F_TSV[i_vf][i_mf][i_y], value_H1full_TSV[i_vr][i_vf][i_mr][i_mf][i_y], value_tgaga_TSV[i_vf][i_mf][i_y], value_LR_TSV[i_vr][i_mr], value_LF_TSV[i_vf][i_mf], LQ);
	      logger << LOG_DEBUG_VERBOSE << "sigma21_TSV at [" << i_vr << "][" << i_vf << "][" << i_mr << "][" << i_mf << "][" << i_y << "]:" << endl;
	      value_sig21_TSV[i_vr][i_vf][i_mr][i_mf][i_y] = calculate_sigma21(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y], value_sig11_TSV[i_vf][i_mf][i_y], value_tH1_TSV[i_vf][i_mf][i_y], value_tH1F_TSV[i_vf][i_mf][i_y], value_H1full_TSV[i_vr][i_vf][i_mr][i_mf][i_y], value_tgaga_TSV[i_vf][i_mf][i_y], value_tcga_TSV[i_vf][i_mf][i_y], value_tgamma2_TSV[i_vf][i_mf][i_y], value_LR_TSV[i_vr][i_mr], value_LF_TSV[i_vf][i_mf], LQ, A_F);
	      for (int i_q = 0; i_q < n_qTcut; i_q++){
		if (switch_resummation){
		  if (i_q <= cut_ps[0]){
		    value_integrand_qTcut_TSV[i_q][i_vr][i_vf][i_mr][i_mf][i_y] 
		      = -ME2 * global_factor * (value_sig21_TSV[i_vr][i_vf][i_mr][i_mf][i_y] * LL1 + 
						value_sig22_TSV[i_vr][i_vf][i_mr][i_mf][i_y] * LL2 + 
						value_sig23_TSV[i_vf][i_mf][i_y] * LL3 + 
						value_sig24_TSV[i_vf][i_mf][i_y] * LL4) / pow(QT_Qres,2) * psi_QT_jacqt2;
		  }
		  else {value_integrand_qTcut_TSV[i_q][i_vr][i_vf][i_mr][i_mf][i_y] = 0.;}
		}
		else {
		  value_integrand_qTcut_TSV[i_q][i_vr][i_vf][i_mr][i_mf][i_y] 
		    = -ME2 * global_factor * (value_sig21_TSV[i_vr][i_vf][i_mr][i_mf][i_y] * I1_int[i_q] + 
					      value_sig22_TSV[i_vr][i_vf][i_mr][i_mf][i_y] * I2_int[i_q] + 
					      value_sig23_TSV[i_vf][i_mf][i_y] * I3_int[i_q] + 
					      value_sig24_TSV[i_vf][i_mf][i_y] * I4_int[i_q]);
		}
	      }
	    logger << LOG_DEBUG_VERBOSE << "value_integrand_qTcut_TSV[0][" << i_vr << "][" << i_vf << "][" << i_mr << "][" << i_mf << "][0] = " << value_integrand_qTcut_TSV[0][i_vr][i_vf][i_mr][i_mf][0] / (-ME2) << endl;
	    }
	  }
	}
      }
    }

    for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
      for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
	for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
	  for (int i_mr = 0; i_mr < n_scale_dyn_ren[i_vr]; i_mr++){
	    logger << LOG_DEBUG_VERBOSE << "QT_H1_delta                     = " << QT_H1_delta << endl;
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
	    logger << LOG_DEBUG_VERBOSE << "value_sig21_TSV    [" << i_vr << "][" << i_vf << "][" << i_mr << "][" << i_mf << "] = " << value_sig21_TSV[i_vr][i_vf][i_mr][i_mf][0] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_sig22_TSV    [" << i_vr << "][" << i_vf << "][" << i_mr << "][" << i_mf << "] = " << value_sig22_TSV[i_vr][i_vf][i_mr][i_mf][0] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_sig23_TSV    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig23_TSV[i_vf][i_mf][0] << endl;
	    logger << LOG_DEBUG_VERBOSE << "value_sig24_TSV    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig24_TSV[i_vf][i_mf][0] << endl;
	    logger.newLine(LOG_DEBUG_VERBOSE);
   	  }
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




