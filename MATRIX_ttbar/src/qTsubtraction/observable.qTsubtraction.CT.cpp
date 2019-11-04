#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

void observable_set::determine_integrand_CX_ncollinear_CT(phasespace_set & psi){
  Logger logger("observable_set::determine_integrand_CX_ncollinear_CT");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  static int order = 1;

  double qtoverq = sqrt(psi_QT_qt2 / QT_Qres / QT_Qres);
  double LL1, LL2, LL3, LL4;
  if (switch_resummation){Itilde(qtoverq, order, LL1, LL2, LL3, LL4);}
  // in FO computations, we do the qT integrals once and for all in the beginning
  // in resummed computations, we have to bin the "artificial" qT, so that is not easily possible -> revert back to old implementation  

  double global_factor = alpha_S / pi;

  if (QT_finalstate_massive_coloured){
    //    calculate_Gamma1();
    if (initial_pdf_diag){calculate_Gamma1born();}
  }

  determine_splitting_tH1F(psi);
  double LQ = log(psi_xbs_all[0][0] / pow(QT_Qres, 2));
  if (QT_Qres == 0){LQ = 0.;}




  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
      for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
	coll_tH1F_pdf[i_vf][i_mf][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][0] * coll_tH1F[i_l];
	logger << LOG_DEBUG_VERBOSE << "coll_tH1F[" << i_l << "] = " << coll_tH1F[i_l] << endl;
	logger << LOG_DEBUG_VERBOSE << "value_list_pdf_factor[" << i_vf << "][" << i_mf << "][" << i_l << "][0] = " << value_list_pdf_factor[i_vf][i_mf][i_l][0] << endl;
      }
      value_tH1F[i_vf][i_mf] = accumulate(coll_tH1F_pdf[i_vf][i_mf].begin(), coll_tH1F_pdf[i_vf][i_mf].end(), 0.);
      value_sig12[i_vf][i_mf] = calculate_sigma12(value_list_pdf_factor[i_vf][i_mf][0][0]);
      value_sig11[i_vf][i_mf] = calculate_sigma11(value_list_pdf_factor[i_vf][i_mf][0][0], value_tH1F[i_vf][i_mf], LQ);

      if (switch_distribution){
	for (int i_y = 0; i_y < 3; i_y++){
	  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
	    coll_tH1F_pdf_D[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor[i_vf][i_mf][i_l][i_y] * coll_tH1F[i_l];
	  }
	  value_tH1F_D[i_vf][i_mf][i_y] = accumulate(coll_tH1F_pdf_D[i_vf][i_mf][i_y].begin(), coll_tH1F_pdf_D[i_vf][i_mf][i_y].end(), 0.);
	  value_sig12_D[i_vf][i_mf][i_y] = calculate_sigma12(value_list_pdf_factor[i_vf][i_mf][0][i_y]);
	  value_sig11_D[i_vf][i_mf][i_y] = calculate_sigma11(value_list_pdf_factor[i_vf][i_mf][0][i_y], value_tH1F_D[i_vf][i_mf][i_y], LQ);
	}
      }
    }
  }

  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
      logger << LOG_DEBUG_VERBOSE << "value_sig11    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig11[i_vf][i_mf] << endl;
      logger << LOG_DEBUG_VERBOSE << "value_sig12    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig12[i_vf][i_mf] << endl;
      logger << LOG_DEBUG_VERBOSE << "value_tH1F     [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1F[i_vf][i_mf] << endl;
      logger.newLine(LOG_DEBUG_VERBOSE);
    }
  }


  if (switch_resummation){
    for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
      for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
	for (int i_q = 0; i_q < n_qTcut; i_q++){
	  if (i_q <= cut_ps[0]){
	    CX_value_integrand_qTcut[i_q][i_vf][i_mf] = global_factor * (value_sig12[i_vf][i_mf] * LL2 + value_sig11[i_vf][i_mf] * LL1) / pow(QT_Qres,2) * psi_QT_jacqt2;
	  }
	  else {CX_value_integrand_qTcut[i_q][i_vf][i_mf] = 0.;}
	  if (switch_distribution){
	    for (int i_y = 0; i_y < 3; i_y++){
              if (i_q <= cut_ps[0]){
		CX_value_integrand_qTcut_D[i_q][i_vf][i_mf][i_y] = global_factor * (value_sig12_D[i_vf][i_mf][i_y] * LL2 + value_sig11_D[i_vf][i_mf][i_y] * LL1) / pow(QT_Qres,2) * psi_QT_jacqt2;
              }
              else {CX_value_integrand_qTcut_D[i_q][i_vf][i_mf][i_y] = 0.;}
	    }
	  }
	}
      }
    }
  }
  else {
    for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
      for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
	for (int i_q = 0; i_q < n_qTcut; i_q++){
	  CX_value_integrand_qTcut[i_q][i_vf][i_mf] = global_factor * (value_sig12[i_vf][i_mf] * I2_int[i_q] + value_sig11[i_vf][i_mf] * I1_int[i_q]);
	  if (switch_distribution){
	    for (int i_y = 0; i_y < 3; i_y++){
              CX_value_integrand_qTcut_D[i_q][i_vf][i_mf][i_y] = global_factor * (value_sig12_D[i_vf][i_mf][i_y] * I2_int[i_q] + value_sig11_D[i_vf][i_mf][i_y] * I1_int[i_q]);

	    }
	  }
	}
      }
    }
  }
      
  /*
  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    for (int i_mf = 0; i_mf < value_mu_fact[i_vf].size(); i_mf++){
      for (int i_q = 0; i_q < n_qTcut; i_q++){
	if (switch_resummation){
	  if (i_q <= cut_ps[0]){
	    CX_value_integrand_qTcut[i_q][i_vf][i_mf] = global_factor * (value_sig12[i_vf][i_mf] * LL2 + value_sig11[i_vf][i_mf] * LL1) / pow(QT_Qres,2) * psi_QT_jacqt2;
	  }
	  else {CX_value_integrand_qTcut[i_q][i_vf][i_mf] = 0.;}
	}
	else {
	  CX_value_integrand_qTcut[i_q][i_vf][i_mf] = global_factor * (value_sig12[i_vf][i_mf] * I2_int[i_q] + value_sig11[i_vf][i_mf] * I1_int[i_q]);
	}
	logger << LOG_DEBUG_VERBOSE << "CX_value_integrand_qTcut[" << i_q << "][" << i_vf << "][" << i_mf << "]  = " << CX_value_integrand_qTcut[i_q][i_vf][i_mf] << endl;
      }
 
      if (switch_distribution){
        for (int i_y = 0; i_y < 3; i_y++){
	  for (int i_q = 0; i_q < n_qTcut; i_q++){
            if (switch_resummation){
              if (i_q <= cut_ps[0]){
		CX_value_integrand_qTcut_D[i_q][i_vf][i_mf][i_y] = global_factor * (value_sig12_D[i_vf][i_mf][i_y] * LL2 + value_sig11_D[i_vf][i_mf][i_y] * LL1) / pow(QT_Qres,2) * psi_QT_jacqt2;
              }
              else {CX_value_integrand_qTcut_D[i_q][i_vf][i_mf][i_y] = 0.;}
            }
            else {
              CX_value_integrand_qTcut_D[i_q][i_vf][i_mf][i_y] = global_factor * (value_sig12_D[i_vf][i_mf][i_y] * I2_int[i_q] + value_sig11_D[i_vf][i_mf][i_y] * I1_int[i_q]);
            }
          }
        }
      }
    }
  }
  */

  // var_rel_alpha_S_CV[i_s] == value_factor_alpha_S[dynamic_scale_CV][map_value_scale_ren_CV[i_s]]
  // var_rel_alpha_S == value_factor_alpha_S[dynamic_scale][map_value_scale_ren]


  integrand = -var_rel_alpha_S * psi_ps_factor * rescaling_factor_alpha_e * ME2 * CX_value_integrand_qTcut[0][dynamic_scale][map_value_scale_fact];

  if (switch_CV){
    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      integrand_CV[i_s] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * rescaling_factor_alpha_e * ME2 * CX_value_integrand_qTcut[0][dynamic_scale_CV][map_value_scale_fact_CV[i_s]];
      for (int i_q = 0; i_q < n_qTcut; i_q++){
	integrand_qTcut_CV[i_q][i_s] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * rescaling_factor_alpha_e * ME2 * CX_value_integrand_qTcut[i_q][dynamic_scale_CV][map_value_scale_fact_CV[i_s]];
      }
    }
  }
  
  if (switch_distribution){
    int x_q = 0;
    for (int j = 0; j < 3; j++){
      integrand_D[j][0] = -var_rel_alpha_S * psi_ps_factor * rescaling_factor_alpha_e * ME2 * CX_value_integrand_qTcut_D[x_q][dynamic_scale][map_value_scale_fact][j];
      logger << LOG_DEBUG_VERBOSE << "new: integrand_D[" << j << "][0] = " << integrand_D[j][0] << endl;
      if (switch_CV){
	for (int i_s = 0; i_s < n_scales_CV; i_s++){
	  integrand_D_CV[i_s][j][0] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * rescaling_factor_alpha_e * ME2 * CX_value_integrand_qTcut_D[x_q][dynamic_scale_CV][map_value_scale_fact_CV[i_s]][j];
	  for (int i_q = 0; i_q < n_qTcut; i_q++){
	    integrand_D_qTcut_CV[i_q][i_s][j][0] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * rescaling_factor_alpha_e * ME2 * CX_value_integrand_qTcut_D[i_q][dynamic_scale_CV][map_value_scale_fact_CV[i_s]][j];
	  }
	}
      }
    }
  }



  if (switch_TSV){
    for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
      for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
	for (int i_y = 0; i_y < 3; i_y++){
	  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
	    coll_tH1F_pdf_TSV[i_vf][i_mf][i_y][i_l] = value_list_pdf_factor_TSV[0][i_vf][i_mf][i_l][i_y] * coll_tH1F[i_l];
	  }
	  value_tH1F_TSV[i_vf][i_mf][i_y] = accumulate(coll_tH1F_pdf_TSV[i_vf][i_mf][i_y].begin(), coll_tH1F_pdf_TSV[i_vf][i_mf][i_y].end(), 0.);
	  value_sig12_TSV[i_vf][i_mf][i_y] = calculate_sigma12(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y]);
	  value_sig11_TSV[i_vf][i_mf][i_y] = calculate_sigma11(value_list_pdf_factor_TSV[0][i_vf][i_mf][0][i_y], value_tH1F_TSV[i_vf][i_mf][i_y], LQ);

	  for (int i_q = 0; i_q < n_qTcut; i_q++){
	    if (switch_resummation){
	      if (i_q <= cut_ps[0]){
		value_fact_integrand_qTcut_TSV[i_q][i_vf][i_mf][i_y] = -ME2 * global_factor * (value_sig12_TSV[i_vf][i_mf][i_y] * LL2 + value_sig11_TSV[i_vf][i_mf][i_y] * LL1) / pow(QT_Qres,2) * psi_QT_jacqt2;
	      }
	      else {value_fact_integrand_qTcut_TSV[i_q][i_vf][i_mf][i_y] = 0.;}
	    }
	    else {
	      value_fact_integrand_qTcut_TSV[i_q][i_vf][i_mf][i_y] = -ME2 * global_factor * (value_sig12_TSV[i_vf][i_mf][i_y] * I2_int[i_q] + value_sig11_TSV[i_vf][i_mf][i_y] * I1_int[i_q]);
	    }
	  }
	}
      }
    }

    for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
      for (int i_mf = 0; i_mf < n_scale_dyn_fact[i_vf]; i_mf++){
	logger << LOG_DEBUG_VERBOSE << "value_sig11_TSV    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig11_TSV[i_vf][i_mf][0] << endl;
	logger << LOG_DEBUG_VERBOSE << "value_sig12_TSV    [" << i_vf << "]   [" << i_mf << "]    = " << value_sig12_TSV[i_vf][i_mf][0] << endl;
	logger << LOG_DEBUG_VERBOSE << "value_tH1F_TSV     [" << i_vf << "]   [" << i_mf << "]    = " << value_tH1F_TSV[i_vf][i_mf][0] << endl;
      }
    }
  }

  // shifted elsewhere ???

  //#ifdef MORE
  //  if (switch_resummation){
  //    double Q = sqrt(psi_xbs_all[0][0]);
  //    double y = log(sqrt(psi_x_pdf[1] / psi_x_pdf[2]));
  //    performQTBoost(psi_QT_qt2, Q, y, particle_event);
  //  }
  //#endif

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

