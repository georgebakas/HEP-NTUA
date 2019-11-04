#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

void observable_set::determine_psp_weight_CT(phasespace_set & psi)
 {
   Logger logger("observable_set::determine_psp_weight_CT");
   logger << LOG_DEBUG_VERBOSE << "started" << endl;
   
   double LR = log(psi_xbs_all[0][0] / var_mu_ren  / var_mu_ren);
   double LF = log(psi_xbs_all[0][0] / var_mu_fact / var_mu_fact);
   double LQ = log(psi_xbs_all[0][0] / QT_Qres     / QT_Qres);
   if (QT_Qres == 0) LQ = 0;
   double A_F = 0, sig = 0;
   
   // the counterterm
   if(name_process[0] == 'g')
     {
       calculate_sigma_gg( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_x_pdf[1], psi_x_pdf[2], pdf_factor[0], QT_pdf_factor_z1x2[0], QT_pdf_factor_x1z2[0],  QT_pdf_factor_qx2[0], QT_pdf_factor_x1q[0], QT_pdf_factor_z1z2[0], QT_pdf_factor_qz2[0], QT_pdf_factor_z1q[0], QT_pdf_factor_qq[0], psi_contribution_order_alpha_s[0], A_F, LR, LF, QT_sigma_qT);
     }
   else
     {
       sig = calculate_sigma_qqbar(psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1, psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor[0], QT_pdf_factor_z1x2[0], QT_pdf_factor_x1z2[0],  QT_pdf_factor_gx2[0], QT_pdf_factor_x1g[0], QT_pdf_factor_gg[0], QT_pdf_factor_qx2[0], QT_pdf_factor_x1q[0], QT_pdf_factor_z1z2[0], QT_pdf_factor_gz2[0], QT_pdf_factor_z1g[0], QT_pdf_factor_qbx2[0], QT_pdf_factor_x1qb[0], psi_contribution_order_alpha_s[0], A_F, LR, LF, LQ, QT_sigma_qT, psi);
     }

   if (switch_resummation == 1) 
     integrand = -var_rel_alpha_S * psi_ps_factor * ME2 * sig * psi_QT_jacqt2;
   else 
     integrand = -var_rel_alpha_S * psi_ps_factor * ME2 * QT_sigma_qT[0];
  
   logger << LOG_DEBUG << "old: var_rel_alpha_S = " << var_rel_alpha_S << endl;
   logger << LOG_DEBUG << "old: psi_ps_factor = " << psi_ps_factor << endl;
   logger << LOG_DEBUG << "old: ME2 = " << ME2 << endl; 
   logger << LOG_DEBUG << "old: QT_sigma_qT[0] = " << QT_sigma_qT[0] << endl;


   this_psp_weight   = integrand;
   this_psp_weight2  = pow(this_psp_weight, 2);
   step_sum_weight  += this_psp_weight;
   step_sum_weight2 += this_psp_weight2;   
   
   if (switch_CV != 0)
     {
       for (int s = 0; s < n_scales_CV; s++)
	 {
	   LR = log(psi_xbs_all[0][0] / value_mu_ren[dynamic_scale_CV ][map_value_scale_ren_CV[ s]] / value_mu_ren[ dynamic_scale_CV][map_value_scale_ren_CV[ s]]);
	   LF = log(psi_xbs_all[0][0] / value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]] / value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]]);
	   
	   if(name_process[0] == 'g')
	     {
	       calculate_sigma_gg( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][0], QT_pdf_factor_z1x2_CV[s][0], QT_pdf_factor_x1z2_CV[s][0], QT_pdf_factor_qx2_CV[s][0], QT_pdf_factor_x1q_CV[s][0], QT_pdf_factor_z1z2_CV[s][0], QT_pdf_factor_qz2_CV[s][0], QT_pdf_factor_z1q_CV[s][0], QT_pdf_factor_qq_CV[s][0], psi_contribution_order_alpha_s[0], A_F, LR, LF, QT_sigma_qT_CV[s] );
	     }
	   else
	     {
	       QT_sig_CV[s] = calculate_sigma_qqbar( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1, psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][0], QT_pdf_factor_z1x2_CV[s][0], QT_pdf_factor_x1z2_CV[s][0], QT_pdf_factor_gx2_CV[s][0], QT_pdf_factor_x1g_CV[s][0], QT_pdf_factor_gg[0], QT_pdf_factor_qx2[0], QT_pdf_factor_x1q[0], QT_pdf_factor_z1z2[0], QT_pdf_factor_gz2[0], QT_pdf_factor_z1g[0], QT_pdf_factor_qbx2[0], QT_pdf_factor_x1qb[0], psi_contribution_order_alpha_s[0], A_F, LR, LF, LQ, QT_sigma_qT_CV[s], psi );    
	     }
	   
	   for (int c = 0; c < n_qTcut; c++)
	     {
	       if (switch_resummation == 1) 
		 {
		   this_psp_weight_CV[c][s] = 0;
		   if (c <= cut_ps[0])
		     this_psp_weight_CV[c][s] = -var_rel_alpha_S_CV[s] * psi_ps_factor * ME2 * QT_sig_CV[s] * psi_QT_jacqt2;
		 } 
	       else 
		 {
		   this_psp_weight_CV[c][s] = -var_rel_alpha_S_CV[s] * psi_ps_factor * ME2 * QT_sigma_qT_CV[s][c];
		 }
	       this_psp_weight2_CV[c][s]  = pow(this_psp_weight_CV[c][s], 2);
	       step_sum_weight_CV[ c][s] += this_psp_weight_CV[ c][s];
	       step_sum_weight2_CV[c][s] += this_psp_weight2_CV[c][s];
	     }
	 }
     }

   //-------------------------------------------------------------------------------------------------
  
  if (switch_distribution != 0)
    {    
      //  integrand_D[0][0] = -var_rel_alpha_S * psi_ps_factor * ME2 * sig * psi_QT_jacqt2;
      //  integrand_D[0][0] = -var_rel_alpha_S * psi_ps_factor * ME2 * QT_sigma_qT[0];
      if (switch_resummation == 1) 
	{
	  integrand_D[0][0] = -var_rel_alpha_S * psi_ps_factor * ME2 * sig * psi_QT_jacqt2;
	} 
      else 
	{
	  integrand_D[0][0] = -var_rel_alpha_S * psi_ps_factor * ME2 *  QT_sigma_qT[0];
	}
      
      if (switch_CV != 0)
	{
	  for (int s = 0; s < n_scales_CV; s++)
	    {
	      //integrand_D_CV[i_s][0][0] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * QT_sig_CV[i_s] * psi_QT_jacqt2;
	      //      integrand_D_CV[i_s][0][0] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * QT_sigma_qT_CV[i_s][0];
	      if (switch_resummation == 1) 
		{
		  integrand_D_CV[s][0][0] = -var_rel_alpha_S_CV[s] * psi_ps_factor * ME2 * QT_sig_CV[s] * psi_QT_jacqt2;
		} 
	      else 
		{
		  integrand_D_CV[s][0][0] = -var_rel_alpha_S_CV[s] * psi_ps_factor * ME2 * QT_sigma_qT_CV[s][0];
		}
	    }
	}
      
      for (int j = 1; j < 3; j++)
	{
	  // why - and why here ???
	  double LR = log(psi_xbs_all[0][0] / var_mu_ren  / var_mu_ren);
	  double LF = log(psi_xbs_all[0][0] / var_mu_fact / var_mu_fact);
	  //
	  if(name_process[0] == 'g')
	    {
	      calculate_sigma_gg( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_x_pdf[1], psi_x_pdf[2], pdf_factor[j], QT_pdf_factor_z1x2[j], QT_pdf_factor_x1z2[j],  QT_pdf_factor_qx2[j], QT_pdf_factor_x1q[j], QT_pdf_factor_z1z2[j], QT_pdf_factor_qz2[j], QT_pdf_factor_z1q[j], QT_pdf_factor_qq[j], psi_contribution_order_alpha_s[0], A_F, LR, LF, QT_sigma_qT);
	    }
	  else
	    {
	      sig = calculate_sigma_qqbar( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1, psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor[j], QT_pdf_factor_z1x2[j], QT_pdf_factor_x1z2[j], QT_pdf_factor_gx2[j], QT_pdf_factor_x1g[j], QT_pdf_factor_gg[j],   QT_pdf_factor_qx2[j], QT_pdf_factor_x1q[j], QT_pdf_factor_z1z2[j], QT_pdf_factor_gz2[j], QT_pdf_factor_z1g[j], QT_pdf_factor_qbx2[j], QT_pdf_factor_x1qb[j], psi_contribution_order_alpha_s[0], A_F, LR, LF, LQ, QT_sigma_qT, psi);
	    }
	  //    integrand_D[j][0] = -var_rel_alpha_S * psi_ps_factor * ME2 * sig * psi_QT_jacqt2;
	  //    integrand_D[j][0] = -var_rel_alpha_S * psi_ps_factor * ME2 * QT_sigma_qT[0];
	  if (switch_resummation == 1) 
	    {
	      integrand_D[j][0] = -var_rel_alpha_S * psi_ps_factor * ME2 * sig * psi_QT_jacqt2;
	    } 
	  else 
	    {
	      integrand_D[j][0] = -var_rel_alpha_S * psi_ps_factor * ME2 * QT_sigma_qT[0];
	    }
	  if (switch_CV != 0)
	    {
	      for (int s = 0; s < n_scales_CV; s++)
		{
		  LR = log(psi_xbs_all[0][0] / value_mu_ren[dynamic_scale_CV][map_value_scale_ren_CV[s]] / value_mu_ren[dynamic_scale_CV][map_value_scale_ren_CV[s]]);
		  LF = log(psi_xbs_all[0][0] / value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]] / value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]]);
		  
		  if(name_process[0] == 'g')
		    {
		      calculate_sigma_gg(                   psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2],                           psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][j], QT_pdf_factor_z1x2_CV[s][j], QT_pdf_factor_x1z2_CV[s][j], QT_pdf_factor_qx2_CV[s][j], QT_pdf_factor_x1q_CV[s][j], QT_pdf_factor_z1z2_CV[s][j], QT_pdf_factor_qz2_CV[s][j], QT_pdf_factor_z1q_CV[s][j], QT_pdf_factor_qq_CV[s][j], psi_contribution_order_alpha_s[0], A_F, LR, LF, QT_sigma_qT_CV[s] );
		    }
		  else
		    {
		      QT_sig_CV[s] = calculate_sigma_qqbar( psi_QT_qt2, psi_xbs_all[0][0], psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1, psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][j], QT_pdf_factor_z1x2_CV[s][j], QT_pdf_factor_x1z2_CV[s][j], QT_pdf_factor_gx2_CV[s][j], QT_pdf_factor_x1g_CV[s][j], QT_pdf_factor_gg_CV[s][j], QT_pdf_factor_qx2_CV[s][j], QT_pdf_factor_x1q_CV[s][j], QT_pdf_factor_z1z2_CV[s][j], QT_pdf_factor_gz2_CV[s][j], QT_pdf_factor_z1g_CV[s][j], QT_pdf_factor_qbx2_CV[s][j], QT_pdf_factor_x1qb_CV[s][j], psi_contribution_order_alpha_s[0], A_F, LR, LF, LQ, QT_sigma_qT_CV[s], psi);
		    }

		  //        integrand_D_CV[i_s][j][0] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * QT_sig_CV[i_s] * psi_QT_jacqt2;
		  //        integrand_D_CV[i_s][j][0] = -var_rel_alpha_S_CV[i_s] * psi_ps_factor * ME2 * QT_sigma_qT_CV[i_s][0];
		  if (switch_resummation == 1) 
		    {
		      integrand_D_CV[s][j][0] = -var_rel_alpha_S_CV[s] * psi_ps_factor * ME2 * QT_sig_CV[s] * psi_QT_jacqt2;
		    } 
		  else 
		    {
		      integrand_D_CV[s][j][0] = -var_rel_alpha_S_CV[s] * psi_ps_factor * ME2 * QT_sigma_qT_CV[s][0];
		    }
		}
	    }
	}
  
//#ifdef MORE
//    if (switch_resummation==1) {
//      double Q = sqrt(psi_xbs_all[0][0]);
//      double y = log(sqrt(psi_x_pdf[1]/psi_x_pdf[2]));
//      
//      performQTBoost(psi_QT_qt2,Q,y,particle_event);
//    }
//#endif
    }

//if (switch_moment)
//  for (int nm = 0; nm < moment.size(); nm++){
//    for (int ia = 0; ia < moment[nm].size(); ia++){
//      ///
//      this_psp_moment[nm] = this_psp_weight * moment[nm][ia];
//      this_psp_moment2[nm] = pow(this_psp_moment[nm], 2.);
//      step_sum_moment[nm] += this_psp_moment[nm];
//      step_sum_moment2[nm] += this_psp_moment2[nm];
//      ///
//      //      psp_moment[nm][psp_counter] = psp_weight[psp_counter] * moment[nm][ia];
//      //      psp_momentdelta[nm][psp_counter] = pow(psp_moment[nm][psp_counter], 2.);
//      if (switch_CV != 0){
//	for (int c = 0; c < n_qTcut; c++){
//	  for (int i_s = 0; i_s < n_scales_CV; i_s++){
//	    if (c <= cut_ps[0]){
//	      ///
//	      this_psp_moment_CV[nm][c][i_s] = this_psp_weight_CV[c][i_s] * moment[nm][ia];
//	      this_psp_moment2_CV[nm][c][i_s] = pow(this_psp_moment_CV[nm][c][i_s], 2.);
//	      step_sum_moment_CV[nm][c][i_s] += this_psp_moment_CV[nm][c][i_s];
//	      step_sum_moment2_CV[nm][c][i_s] += this_psp_moment2_CV[nm][c][i_s];
//	      ///
//	      //	      psp_moment_CV[nm][c][i_s][psp_counter] = psp_weight_CV[c][i_s][psp_counter] * moment[nm][ia];
//	      //	      psp_momentdelta_CV[nm][c][i_s][psp_counter] = pow(psp_moment_CV[nm][c][i_s][psp_counter], 2.);
//	    }
//	  }
//	}
//      }
//    }
//  }
// }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
 }



