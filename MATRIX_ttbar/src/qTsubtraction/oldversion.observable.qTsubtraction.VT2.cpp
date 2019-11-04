#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

void observable_set::determine_integrand_VT2(phasespace_set & psi){
  Logger logger("observable_set::determine_integrand_VT2");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  // see notes, Eq. (12)-(16)
  double LR = log(psi_xbs_all[0][0] / pow(var_mu_ren, 2));
  double LF = log(psi_xbs_all[0][0] / pow(var_mu_fact, 2));
  double LQ = log(psi_xbs_all[0][0] / pow(QT_Qres, 2));
  if (QT_Qres == 0) LQ = 0;
  // careful: LR, LF defined with a minus sign

  if (psi_contribution_order_alpha_s[0] == 1 && QT_H1_delta==0) 
    {
      assert(false);
      QT_H2_delta = 0;
    } 
  else if (psi_contribution_order_alpha_s[0]==2) 
    {      
      logger << LOG_DEBUG_VERBOSE << "QT_H1_delta = " << QT_H1_delta << endl;
      logger << LOG_DEBUG_VERBOSE << "QT_H2_delta = " << QT_H2_delta << endl;
    /*
    if (QT_H2_delta==0) { // temporary, until all ME2.VVG2 are updated
     logger << LOG_DEBUG_VERBOSE << "QT_H2_delta == 0 should not happen !!!" << endl;
     //     calculate_amplitude(p_parton[0],QT_A0,QT_A1,QT_A2,2);
     calculate_amplitude_doublevirtual();

     logger << LOG_DEBUG_VERBOSE << "QT_A0 = " << QT_A0 << endl;
     if (QT_A0==0) {
	QT_A0=VA_b_ME2;
	QT_A1=ME2;
	A_F=f2pi/alpha_S*(VA_V_ME2 + VA_X_ME2)/VA_b_ME2-2.0/3*C_F*pi2-2.0*gamma_q*LR+C_F*LR*LR;
	QT_H1_delta = pi2_6*C_F+0.5*A_F;
	logger << LOG_DEBUG_VERBOSE << "QT_A0 = " << QT_A0 << endl;
	logger << LOG_DEBUG_VERBOSE << "QT_A1 = " << QT_A1 << endl;
	logger << LOG_DEBUG_VERBOSE << "A_F = " << A_F << endl;
	logger << LOG_DEBUG_VERBOSE << "QT_H1_delta = " << QT_H1_delta << endl;
	//        logger << LOG_DEBUG_VERBOSE << "old: " << VA_b_ME2 << ", " << QT_H1_delta << endl;
	//  exit(0);
      } 
      else {
	double q2=2*p_parton[0][3]*p_parton[0][4];
	// FIXME: only works for Vgamma
	calc_H_coefficients(q2,psi_xbs_all[0][0],QT_A0,QT_A1,QT_A2,QT_H1_delta,QT_H2_delta,2);
	
	// adapt alpha_S/Pi normalisation
	QT_H1_delta /= 2;
	QT_H2_delta /= 4;
	
	QT_H1_delta /= (alpha_S/2/pi*QT_A0);
	QT_H2_delta /= (pow(alpha_S/2/pi,2)*QT_A0);
	
	A_F=f2pi/alpha_S*(VA_V_ME2 + VA_X_ME2)/VA_b_ME2-2.0/3*C_F*pi2-2.0*gamma_q*LR+C_F*LR*LR;
	//      logger << LOG_DEBUG_VERBOSE << "H1: " << (pi2_6*C_F+0.5*A_F)*VA_b_ME2 << ", " << QT_H1_delta*QT_A0 << ", " <<  (pi2_6*C_F+0.5*A_F)/QT_H1_delta*VA_b_ME2/QT_A0 << endl;
      }
    }
*/
      logger << LOG_DEBUG_VERBOSE << "VA_b_ME2 = " << VA_b_ME2 << endl;
      VA_b_ME2 = QT_A0;
      logger << LOG_DEBUG_VERBOSE << "VA_b_ME2 = " << VA_b_ME2 << endl;
    }
    
  //  logger << LOG_DEBUG_VERBOSE << "A0/1/2 = " << QT_A0 << ", " << QT_A1 << ", " << QT_A2 << endl;
  //  logger << LOG_DEBUG_VERBOSE << "H1/2 = " << QT_H1_delta << ", " << QT_H2_delta << endl;
  
  double virt = 0.;
  
  if (switch_resummation==0) 
    {
      if(name_process[0] == 'g')
	{
	  virt = calculate_virtual_gg( psi_zz_pdf[1], psi_zz_pdf[2], psi_x_pdf[1], psi_x_pdf[2], pdf_factor[0], QT_pdf_factor_z1x2[0], QT_pdf_factor_x1z2[0], QT_pdf_factor_qx2[0], QT_pdf_factor_x1q[0], QT_pdf_factor_z1z2[0], QT_pdf_factor_qz2[0], QT_pdf_factor_z1q[0], QT_pdf_factor_qq[0], psi_contribution_order_alpha_s[0], QT_H1_delta, QT_H2_delta,  LR, LF );
	}
      else
	{
	  virt = virtual_qqbar(  psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1, psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor[0], QT_pdf_factor_z1x2[0], QT_pdf_factor_x1z2[0], QT_pdf_factor_gx2[0], QT_pdf_factor_x1g[0], QT_pdf_factor_gg[0],QT_pdf_factor_qx2[0], QT_pdf_factor_x1q[0], QT_pdf_factor_z1z2[0], QT_pdf_factor_gz2[0], QT_pdf_factor_z1g[0], QT_pdf_factor_qbx2[0], QT_pdf_factor_x1qb[0], psi_contribution_order_alpha_s[0], QT_H1_delta, QT_H2_delta, LR, LF, LQ);
	}
      integrand = psi_ps_factor * VA_b_ME2 * var_rel_alpha_S * virt;
    } 
  else 
    {
#ifdef MORE
    double Q = sqrt(psi_xbs_all[0][0]);
    double y = log(sqrt(psi_x_pdf[1]/psi_x_pdf[2]));
    double muF = var_mu_fact;
    double muRen = var_mu_ren;
    int gg_initiated;
    int channel = 0;

    double sqrt_shad = 2*psi_E;
    double qt = sqrt(psi_QT_qt2);

    //performQTBoost(psi_QT_qt2,Q,y,particle_event);

    logger << LOG_DEBUG << "Q=" << Q << ", " << "y=" << y << ", " << "qt=" << qt << ", x1=" << x_pdf[1] << ", x2=" << x_pdf[2] << endl;

    if (abs(y)>-0.5*log(pow(somescale_.startscale/sqrt_shad,2))) {
      cout << "y too large, skipping point" << endl;
      virt=0.0;
      integrand=0.0;
    } else { 
      virt=0;
      string initial_state = name_process.substr(0,4);
      if (initial_state == "uu~_") {
        gg_initiated=-1;
        resumm_(&virt,&gg_initiated,&channel,&Q,&y,&qt,&QT_Qres,&muF,&muRen,&sqrt_shad,&VA_b_ME2,&QT_H1_delta,&QT_H2_delta);
      } else if (initial_state == "dd~_") {
        gg_initiated=-2;
        resumm_(&virt,&gg_initiated,&channel,&Q,&y,&qt,&QT_Qres,&muF,&muRen,&sqrt_shad,&VA_b_ME2,&QT_H1_delta,&QT_H2_delta);
      } else {
        assert(false);
      }
      if (isnan(virt)) {
        cout << "virt is nan, setting to zero!" << endl;
        virt=0;
      }

      integrand = psi_ps_factor * virt/x_pdf[0] * psi_QT_jacqt2/2/qt;
      logger << LOG_DEBUG << "event " << psi_i_gen << ": virt=" << virt << ", integrand=" << psi_ps_factor * virt/x_pdf[0] * psi_QT_jacqt2/2/qt << endl;
    }
#endif

      /*
#ifdef MORE
    double Q = sqrt(psi_xbs_all[0][0]);
    double y = log(sqrt(psi_x_pdf[1]/psi_x_pdf[2]));
    double muF = var_mu_fact;
    double muRen = var_mu_ren;
    int gg_initiated;
    int channel = 0;

    double sqrt_shad = 2*psi_E;
    double qt = sqrt(psi_QT_qt2);

    performQTBoost(psi_QT_qt2,Q,y,particle_event);

    cout << "Q=" << Q << ", " << "y=" << y << ", " << "qt=" << qt << ", x1=" << x_pdf[1] << ", x2=" << x_pdf[2] << endl;

    if (abs(y)>-0.5*log(pow(somescale_.startscale/sqrt_shad,2))) {
      cout << "y too large, skipping point" << endl;
      virt=0.0;
      integrand=0.0;
    } else { 
      virt=0;
      // FIXME
      if (name_process == "uu~_zz" || name_process == "uu~_wpwm") {
        gg_initiated=-1;
        resumm_(&virt,&gg_initiated,&channel,&Q,&y,&qt,&QT_Qres,&muF,&muRen,&sqrt_shad,&VA_b_ME2,&QT_H1_delta,&QT_H2_delta);
      } else if (name_process == "dd~_zz" || name_process == "dd~_wpwm") {
        gg_initiated=-2;
        resumm_(&virt,&gg_initiated,&channel,&Q,&y,&qt,&QT_Qres,&muF,&muRen,&sqrt_shad,&VA_b_ME2,&QT_H1_delta,&QT_H2_delta);
      } else {
        assert(false);
      }
      if (isnan(virt)) {
        cout << "virt is nan, setting to zero!" << endl;
        virt=0;
      }

      integrand = psi_ps_factor * virt/x_pdf[0] * psi_QT_jacqt2/2/qt;
      cout << "event " << "?" << ": virt=" << virt << ", integrand=" << psi_ps_factor * virt/x_pdf[0] * psi_QT_jacqt2/2/qt << endl;
    }
#endif
      */
    }

  if (switch_CV != 0)
    {
      for (int s = 0; s < n_scales_CV; s++)
	{
	  LR = log(psi_xbs_all[0][0]/value_mu_ren[dynamic_scale_CV][map_value_scale_ren_CV[s]]  /value_mu_ren[dynamic_scale_CV][map_value_scale_ren_CV[s]]);
	  LF = log(psi_xbs_all[0][0]/value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]]/value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]]);
	  logger << LOG_DEBUG_VERBOSE << "s = " << s << endl;
	  
	  //    cout << s << ": " << value_mu_ren[dynamic_scale_CV][map_value_scale_ren_CV[s]] << ", " << value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]] << endl;
	
	  if(name_process[0] == 'g')
	    {            
	      QT_virt_CV[s] = calculate_virtual_gg( psi_zz_pdf[1], psi_zz_pdf[2], psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][0], QT_pdf_factor_z1x2_CV[s][0], QT_pdf_factor_x1z2_CV[s][0], QT_pdf_factor_qx2_CV[s][0], QT_pdf_factor_x1q_CV[s][0], QT_pdf_factor_z1z2_CV[s][0], QT_pdf_factor_qz2_CV[s][0], QT_pdf_factor_z1q_CV[s][0], QT_pdf_factor_qq_CV[s][0],  psi_contribution_order_alpha_s[0], QT_H1_delta, QT_H2_delta,  LR, LF );
	    }
          else
	    {
	      QT_virt_CV[s] = virtual_qqbar( psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1, psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][0], QT_pdf_factor_z1x2_CV[s][0], QT_pdf_factor_x1z2_CV[s][0], QT_pdf_factor_gx2_CV[s][0], QT_pdf_factor_x1g_CV[s][0], QT_pdf_factor_gg_CV[s][0], QT_pdf_factor_qx2_CV[s][0], QT_pdf_factor_x1q_CV[s][0], QT_pdf_factor_z1z2_CV[s][0], QT_pdf_factor_gz2_CV[s][0], QT_pdf_factor_z1g_CV[s][0], QT_pdf_factor_qbx2_CV[s][0], QT_pdf_factor_x1qb_CV[s][0], psi_contribution_order_alpha_s[0], QT_H1_delta, QT_H2_delta, LR, LF, LQ );
	    }
	}
    }


  //-------------------------------------------------------------------------------------------------------------------
  
  if (switch_distribution != 0)
    {
      integrand_D[0][0] = var_rel_alpha_S * psi_ps_factor * VA_b_ME2 * virt;
      if (switch_CV != 0)
	{
	  for (int s = 0; s < n_scales_CV; s++)
	    {
	      integrand_D_CV[s][0][0] = var_rel_alpha_S_CV[s] * psi_ps_factor * VA_b_ME2 * QT_virt_CV[s];
	    }
	}
      for (int j = 1; j < 3; j++)
	{
	  LR = log(psi_xbs_all[0][0]/var_mu_ren/var_mu_ren); 
	  LF = log(psi_xbs_all[0][0]/var_mu_fact/var_mu_fact); 
  
	  if(name_process[0] == 'g')
	    {
	      virt = calculate_virtual_gg( psi_zz_pdf[1], psi_zz_pdf[2], psi_x_pdf[1], psi_x_pdf[2], pdf_factor[j], QT_pdf_factor_z1x2[j], QT_pdf_factor_x1z2[j], QT_pdf_factor_qx2[j], QT_pdf_factor_x1q[j], QT_pdf_factor_z1z2[j], QT_pdf_factor_qz2[j], QT_pdf_factor_z1q[j], QT_pdf_factor_qq[j], psi_contribution_order_alpha_s[0], QT_H1_delta, QT_H2_delta, LR, LF );

	    }
	  else
	    {
	      virt = virtual_qqbar( psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1,psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor[j], QT_pdf_factor_z1x2[j], QT_pdf_factor_x1z2[j], QT_pdf_factor_gx2[j], QT_pdf_factor_x1g[j], QT_pdf_factor_gg[j], QT_pdf_factor_qx2[j], QT_pdf_factor_x1q[j], QT_pdf_factor_z1z2[j], QT_pdf_factor_gz2[j], QT_pdf_factor_z1g[j], QT_pdf_factor_qbx2[j], QT_pdf_factor_x1qb[j], psi_contribution_order_alpha_s[0], QT_H1_delta, QT_H2_delta, LR, LF, LQ );
	    }
	  
	  integrand_D[j][0] = var_rel_alpha_S * psi_ps_factor * VA_b_ME2 * virt;

	  if (switch_CV != 0)
	    {
	      for (int s = 0; s < n_scales_CV; s++)
		{
		  LR = log(psi_xbs_all[0][0]/value_mu_ren[dynamic_scale_CV][map_value_scale_ren_CV[s]]/value_mu_ren[dynamic_scale_CV][map_value_scale_ren_CV[s]]);
		  LF = log(psi_xbs_all[0][0]/value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]]/value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]]);

		  double virt_temp = 0.;
		  if(name_process[0] == 'g')
		    {
		      virt_temp = calculate_virtual_gg( psi_zz_pdf[1], psi_zz_pdf[2], psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][j], QT_pdf_factor_z1x2_CV[s][j], QT_pdf_factor_x1z2_CV[s][j], QT_pdf_factor_qx2_CV[s][j], QT_pdf_factor_x1q_CV[s][j], QT_pdf_factor_z1z2_CV[s][j], QT_pdf_factor_qz2_CV[s][j], QT_pdf_factor_z1q_CV[s][j], QT_pdf_factor_qq_CV[s][j],  psi_contribution_order_alpha_s[0], QT_H1_delta, QT_H2_delta,  LR, LF );
		    }
		  else
		    {
		      virt_temp = virtual_qqbar( psi_zz_pdf[1], psi_zz_pdf[2], psi_QT_g_z1, psi_QT_g_z2, psi_x_pdf[1], psi_x_pdf[2], pdf_factor_CV[s][j], QT_pdf_factor_z1x2_CV[s][j], QT_pdf_factor_x1z2_CV[s][j], QT_pdf_factor_gx2_CV[s][j], QT_pdf_factor_x1g_CV[s][j], QT_pdf_factor_gg_CV[s][j], QT_pdf_factor_qx2_CV[s][j], QT_pdf_factor_x1q_CV[s][j], QT_pdf_factor_z1z2_CV[s][j], QT_pdf_factor_gz2_CV[s][j], QT_pdf_factor_z1g_CV[s][j], QT_pdf_factor_qbx2_CV[s][j], QT_pdf_factor_x1qb_CV[s][j], psi_contribution_order_alpha_s[0], QT_H1_delta, QT_H2_delta, LR, LF, LQ );
		    }
	  
		  integrand_D_CV[s][j][0] = var_rel_alpha_S_CV[s] * psi_ps_factor * VA_b_ME2 * virt_temp;
		}
	    }
	}
    }
  
  QT_A0 = 0.; 
  QT_A1 = 0.; 
  QT_A2 = 0.;
  QT_H1_delta = 0.;
  QT_H2_delta = 0.;
  
  //  exit(0);
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
 }



