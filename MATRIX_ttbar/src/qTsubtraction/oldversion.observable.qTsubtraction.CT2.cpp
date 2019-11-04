#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

// !!! analogous adaptation as in CT !!! ???
void observable_set::determine_psp_weight_CT2(phasespace_set & psi){
  Logger logger("observable_set::determine_psp_weight_CT2");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  double LR = log(psi_xbs_all[0][0] / var_mu_ren  / var_mu_ren );
  double LF = log(psi_xbs_all[0][0] / var_mu_fact / var_mu_fact);
  double LQ = log(psi_xbs_all[0][0] / QT_Qres     / QT_Qres    );
  if (QT_Qres == 0) LQ = 0;

  double A_F = 0.;
  double sig = 0.;
  if(name_process[0] == 'g')
    A_F = 2 * (QT_H1_delta - pi2_6 * C_A);
  else 
    A_F = 2 * (QT_H1_delta - pi2_6 * C_F);


  if(name_process[0] == 'g')
    {
      calculate_sigma_gg( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_x_pdf[1], psi_x_pdf[2], pdf_factor[0], QT_pdf_factor_z1x2[0], QT_pdf_factor_x1z2[0],  QT_pdf_factor_qx2[0], QT_pdf_factor_x1q[0], QT_pdf_factor_z1z2[0], QT_pdf_factor_qz2[0], QT_pdf_factor_z1q[0], QT_pdf_factor_qq[0], psi_contribution_order_alpha_s[0], A_F, LR, LF, QT_sigma_qT);
      logger << LOG_DEBUG_VERBOSE << "sig = " << sig << endl;
    }
  else
    {
      logger << LOG_DEBUG_VERBOSE << "calculate_sigma_qqbar" << endl;
      sig = calculate_sigma_qqbar( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1, psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor[0], QT_pdf_factor_z1x2[0], QT_pdf_factor_x1z2[0], QT_pdf_factor_gx2[0], QT_pdf_factor_x1g[0], QT_pdf_factor_gg[0], QT_pdf_factor_qx2[0], QT_pdf_factor_x1q[0], QT_pdf_factor_z1z2[0], QT_pdf_factor_gz2[0], QT_pdf_factor_z1g[0], QT_pdf_factor_qbx2[0], QT_pdf_factor_x1qb[0], psi_contribution_order_alpha_s[0], A_F, LR, LF, LQ, QT_sigma_qT, psi );
    }

//cout << "sig=" << sig << endl;
//sig=calculate_sigma_qqbar(psi_QT_qt2,psi_xbs_all[0][0],psi_zz_pdf[1],psi_zz_pdf[2],psi_QT_g_z1,psi_QT_g_z2,psi_x_pdf[1],psi_x_pdf[2],pdf_factor[0],QT_pdf_factor_z1x2[0],QT_pdf_factor_x1z2[0],QT_pdf_factor_gx2[0],QT_pdf_factor_x1g[0],QT_pdf_factor_gg[0],0,0,0,0,0,0,0,psi_contribution_order_alpha_s[0],alpha_S,A_F,LR,LF,LQ);
//cout << "-> " << sig << endl << endl;

  integrand = -var_rel_alpha_S * psi_ps_factor * ME2 * QT_sigma_qT[0];

  this_psp_weight = integrand;
  this_psp_weight2  = pow(this_psp_weight, 2);
  step_sum_weight  += this_psp_weight;
  step_sum_weight2 += this_psp_weight2;
  
  if (switch_CV != 0)
    {
      for (int s = 0; s < n_scales_CV; s++)
	{
	  LR  = log(psi_xbs_all[0][0]/value_mu_ren[ dynamic_scale_CV][map_value_scale_ren_CV[ s]]/value_mu_ren[ dynamic_scale_CV][map_value_scale_ren_CV[ s]]);
	  LF  = log(psi_xbs_all[0][0]/value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]]/value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]]);

	  if(name_process[0] == 'g')
	    {
	      calculate_sigma_gg( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2],  psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][0], QT_pdf_factor_z1x2_CV[s][0], QT_pdf_factor_x1z2_CV[s][0], QT_pdf_factor_qx2_CV[s][0], QT_pdf_factor_x1q_CV[s][0], QT_pdf_factor_z1z2_CV[s][0], QT_pdf_factor_qz2_CV[s][0], QT_pdf_factor_z1q_CV[s][0], QT_pdf_factor_qq_CV[s][0], psi_contribution_order_alpha_s[0], A_F, LR, LF, QT_sigma_qT_CV[s]);  
	    }
	  else
	    {
	      logger << LOG_DEBUG_VERBOSE << "calculate_sigma_qqbar s = " << s << endl;
	      QT_sig_CV[s] = calculate_sigma_qqbar( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1, psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][0], QT_pdf_factor_z1x2_CV[s][0], QT_pdf_factor_x1z2_CV[s][0], QT_pdf_factor_gx2_CV[s][0], QT_pdf_factor_x1g_CV[s][0], QT_pdf_factor_gg_CV[s][0], QT_pdf_factor_qx2_CV[s][0], QT_pdf_factor_x1q_CV[s][0], QT_pdf_factor_z1z2_CV[s][0], QT_pdf_factor_gz2_CV[s][0], QT_pdf_factor_z1g_CV[s][0], QT_pdf_factor_qbx2_CV[s][0], QT_pdf_factor_x1qb_CV[s][0], psi_contribution_order_alpha_s[0], A_F, LR, LF, LQ, QT_sigma_qT_CV[s], psi );
	    }
	  for (int c = 0; c < n_qTcut; c++)
	    {
//        if (c <= cut_ps[0]){
//        cout << s << ", " << c << "," << cut << endl;
//        temp_psp_weight_CV[c][s] = -var_rel_alpha_S_CV[s] * psi_ps_factor * ME2* QT_sig_CV[s] * psi_QT_jacqt2;
//        temp_psp_weight_CV[c][s]=0;
//        for (int l=0; l<sigma.size(); l++) {
//          if (c <= cut_qt[l]){
//            temp_psp_weight_CV[c][s] -= var_rel_alpha_S_CV[s] * psi_ps_factor * ME2 * sigma_CV[s][l] * psi_QT_jacqt2[l];
//          }
//        }
	      this_psp_weight_CV[ c][s]  = -var_rel_alpha_S_CV[s] * psi_ps_factor * ME2 * QT_sigma_qT_CV[s][c];
	      this_psp_weight2_CV[c][s]  = pow(this_psp_weight_CV[c][s], 2);
	      step_sum_weight_CV[ c][s] += this_psp_weight_CV[ c][s];
	      step_sum_weight2_CV[c][s] += this_psp_weight2_CV[c][s];
	    }
//      cout << "alphaS=" << alpha_S << ", LR=" << LR << ", LF=" << LF << endl;
//      cout << s << ": " << pdf_factor_CV[s][0] << ", " << QT_pdf_factor_z1x2_CV[s][0] << ", " << QT_pdf_factor_x1z2_CV[s][0] << ", " << QT_pdf_factor_gx2_CV[s][0] << ", " << QT_pdf_factor_x1g_CV[s][0] << endl;
//      cout << s << ": VME2=" << VA_V_ME2 << ", VA_X_ME2_CV[s]=" << VA_X_ME2_CV[s] << ", AF=" << A_F << ", var_rel_alpha_S_CV[s]=" << sqrt(var_rel_alpha_S_CV[s]) << ", pspw=" << temp_psp_weight_CV[0][s] << endl;
	}
    }

  //------------------------------------------------------------------------------------------------------------------------------
  
  if (switch_distribution != 0)
    {

      integrand_D[0][0] = -var_rel_alpha_S * psi_ps_factor * ME2 * QT_sigma_qT[0];
      if (switch_CV != 0)
	{
	  for (int s = 0; s < n_scales_CV; s++)
	    {
	      integrand_D_CV[s][0][0] = -var_rel_alpha_S_CV[s] * psi_ps_factor * ME2 * QT_sigma_qT_CV[s][0];
	    }
	}
      for (int j = 1; j < 3; j++)
	{
	  LR = log(psi_xbs_all[0][0]/var_mu_ren/var_mu_ren);
	  LF = log(psi_xbs_all[0][0]/var_mu_fact/var_mu_fact);	  
	  
	  if(name_process[0] == 'g')
	    {
	      calculate_sigma_gg( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_x_pdf[1], psi_x_pdf[2], pdf_factor[j], QT_pdf_factor_z1x2[j], QT_pdf_factor_x1z2[j],  QT_pdf_factor_qx2[j], QT_pdf_factor_x1q[j], QT_pdf_factor_z1z2[j], QT_pdf_factor_qz2[j], QT_pdf_factor_z1q[j], QT_pdf_factor_qq[j], psi_contribution_order_alpha_s[0], A_F, LR, LF, QT_sigma_qT );
	    }
	  else
	    {
	      logger << LOG_DEBUG_VERBOSE << "calculate_sigma_qqbar   j = " << j << endl;
	      sig = calculate_sigma_qqbar( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1, psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor[j], QT_pdf_factor_z1x2[j], QT_pdf_factor_x1z2[j], QT_pdf_factor_gx2[j], QT_pdf_factor_x1g[j], QT_pdf_factor_gg[j], QT_pdf_factor_qx2[j], QT_pdf_factor_x1q[j], QT_pdf_factor_z1z2[j], QT_pdf_factor_gz2[j], QT_pdf_factor_z1g[j], QT_pdf_factor_qbx2[j], QT_pdf_factor_x1qb[j], psi_contribution_order_alpha_s[0], A_F, LR, LF, LQ, QT_sigma_qT, psi );
	    }

	  integrand_D[j][0] = -var_rel_alpha_S * psi_ps_factor * ME2 * QT_sigma_qT[0];
	  if (switch_CV != 0)
	    {
	      for (int s = 0; s < n_scales_CV; s++)
		{
		  LR = log(psi_xbs_all[0][0]/value_mu_ren[dynamic_scale_CV][map_value_scale_ren_CV[s]]/value_mu_ren[dynamic_scale_CV][map_value_scale_ren_CV[s]]);
		  LF = log(psi_xbs_all[0][0]/value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]]/value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]]);
		  
		  if(name_process[0] == 'g')
		    {
		      calculate_sigma_gg( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][j], QT_pdf_factor_z1x2_CV[s][j], QT_pdf_factor_x1z2_CV[s][j], QT_pdf_factor_qx2_CV[s][j], QT_pdf_factor_x1q_CV[s][j], QT_pdf_factor_z1z2_CV[s][j], QT_pdf_factor_qz2_CV[s][j], QT_pdf_factor_z1q_CV[s][j], QT_pdf_factor_qq_CV[s][j], psi_contribution_order_alpha_s[0], A_F, LR, LF, QT_sigma_qT_CV[s]);
		    }
		  else
		    {
		      logger << LOG_DEBUG_VERBOSE << "calculate_sigma_qqbar   j = " << j << "   s = " << s << endl;
		      QT_sig_CV[s] = calculate_sigma_qqbar( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1, psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][j], QT_pdf_factor_z1x2_CV[s][j], QT_pdf_factor_x1z2_CV[s][j], QT_pdf_factor_gx2_CV[s][j], QT_pdf_factor_x1g_CV[s][j], QT_pdf_factor_gg_CV[s][j], QT_pdf_factor_qx2_CV[s][j], QT_pdf_factor_x1q_CV[s][j], QT_pdf_factor_z1z2_CV[s][j], QT_pdf_factor_gz2_CV[s][j], QT_pdf_factor_z1g_CV[s][j], QT_pdf_factor_qbx2_CV[s][j], QT_pdf_factor_x1qb_CV[s][j], psi_contribution_order_alpha_s[0], A_F, LR, LF, LQ, QT_sigma_qT_CV[s], psi );
		    }
		  
		  integrand_D_CV[s][j][0] = -var_rel_alpha_S_CV[s] * psi_ps_factor * ME2 * QT_sigma_qT_CV[s][0];
		}
	    }
	}
    }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



