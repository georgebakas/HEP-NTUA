#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"
#include "../include/definitions.phasespace.set.cxx"

void calculate_pdf_LHAPDF_QT(vector<vector<int> > & all_pdf, phasespace_set & psi, observable_set & oset, int order){
  static Logger logger("calculate_pdf_LHAPDF_QT");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  
  for (int sd = 0; sd < osi_value_mu_fact.size(); sd++){
    for (int ss = 0; ss < osi_value_mu_fact[sd].size(); ss++){
      //osi_value_mu_fact[sd][ss], osi_QT_value_pdf_factor[sd][ss], osi_QT_value_pdf_factor_z1x2[sd][ss], osi_QT_value_pdf_factor_x1z2[sd][ss], osi_QT_value_pdf_factor_gx2[sd][ss], osi_QT_value_pdf_factor_x1g[sd][ss], osi_QT_value_pdf_factor_gg[sd][ss], osi_QT_value_pdf_factor_qx2[sd][ss], osi_QT_value_pdf_factor_x1q[sd][ss], osi_QT_value_pdf_factor_z1z2[sd][ss], osi_QT_value_pdf_factor_gz2[sd][ss], osi_QT_value_pdf_factor_z1g[sd][ss], osi_QT_value_pdf_factor_qbx2[sd][ss], osi_QT_value_pdf_factor_x1qb[sd][ss], 
      if (psi_coll_choice == 1){calculate_pdf_LHAPDF_LHC_QT(all_pdf, sd, ss, psi, oset, order);}
      else if (psi_coll_choice == 2){calculate_pdf_LHAPDF_Tevatron_QT(all_pdf, sd, ss, psi, oset, order);}
    }
  }

  osi_pdf_factor         = osi_QT_value_pdf_factor[     osi_dynamic_scale][osi_map_value_scale_fact];  
  osi_QT_pdf_factor_z1x2 = osi_QT_value_pdf_factor_z1x2[osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_x1z2 = osi_QT_value_pdf_factor_x1z2[osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_gx2  = osi_QT_value_pdf_factor_gx2[ osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_x1g  = osi_QT_value_pdf_factor_x1g[ osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_gg   = osi_QT_value_pdf_factor_gg[  osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_qx2  = osi_QT_value_pdf_factor_qx2[ osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_x1q  = osi_QT_value_pdf_factor_x1q[ osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_z1z2 = osi_QT_value_pdf_factor_z1z2[osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_gz2  = osi_QT_value_pdf_factor_gz2[ osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_z1g  = osi_QT_value_pdf_factor_z1g[ osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_qbx2 = osi_QT_value_pdf_factor_qbx2[osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_x1qb = osi_QT_value_pdf_factor_x1qb[osi_dynamic_scale][osi_map_value_scale_fact];
  // new                                                
  osi_QT_pdf_factor_qq   = osi_QT_value_pdf_factor_qq[  osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_qz2  = osi_QT_value_pdf_factor_qz2[ osi_dynamic_scale][osi_map_value_scale_fact];
  osi_QT_pdf_factor_z1q  = osi_QT_value_pdf_factor_z1q[ osi_dynamic_scale][osi_map_value_scale_fact];


  if (osi_switch_CV){
      for (int i_s = 0; i_s < osi_n_scales_CV; i_s++){
	  osi_pdf_factor_CV[        i_s] = osi_QT_value_pdf_factor[     osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];	  
	  osi_QT_pdf_factor_z1x2_CV[i_s] = osi_QT_value_pdf_factor_z1x2[osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  osi_QT_pdf_factor_x1z2_CV[i_s] = osi_QT_value_pdf_factor_x1z2[osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  osi_QT_pdf_factor_gx2_CV[ i_s] = osi_QT_value_pdf_factor_gx2[ osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  osi_QT_pdf_factor_x1g_CV[ i_s] = osi_QT_value_pdf_factor_x1g[ osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  osi_QT_pdf_factor_gg_CV[  i_s] = osi_QT_value_pdf_factor_gg[  osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  osi_QT_pdf_factor_qx2_CV[ i_s] = osi_QT_value_pdf_factor_qx2[ osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  osi_QT_pdf_factor_x1q_CV[ i_s] = osi_QT_value_pdf_factor_x1q[ osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  osi_QT_pdf_factor_z1z2_CV[i_s] = osi_QT_value_pdf_factor_z1z2[osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  osi_QT_pdf_factor_gz2_CV[ i_s] = osi_QT_value_pdf_factor_gz2[ osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  osi_QT_pdf_factor_z1g_CV[ i_s] = osi_QT_value_pdf_factor_z1g[ osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  osi_QT_pdf_factor_qbx2_CV[i_s] = osi_QT_value_pdf_factor_qbx2[osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  osi_QT_pdf_factor_x1qb_CV[i_s] = osi_QT_value_pdf_factor_x1qb[osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	  // new                                                          
	  osi_QT_pdf_factor_qq_CV[  i_s] = osi_QT_value_pdf_factor_qq[  osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
          osi_QT_pdf_factor_qz2_CV[ i_s] = osi_QT_value_pdf_factor_qz2[ osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
          osi_QT_pdf_factor_z1q_CV[ i_s] = osi_QT_value_pdf_factor_z1q[ osi_dynamic_scale_CV][osi_map_value_scale_fact_CV[i_s]];
	}
    }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void calculate_pdf_LHAPDF_LHC_QT(vector<vector<int> > & all_pdf, int sd, int ss, phasespace_set & psi, observable_set & oset, int order) {
  static Logger logger("calculate_pdf_LHAPDF_LHC_QT");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<vector<double> > lhapdf_result(3);
  vector<vector<double> > lhapdf_result_z(3);

  for (int j = 1; j < 3; j++)
    {
      //    cout << "psi_x_pdf[" << j << "] = " << psi_x_pdf[j] << "osi_value_mu_fact[" << sd << "][" << ss << "] = " << osi_value_mu_fact[sd][ss] << endl;
      lhapdf_result[j] = LHAPDF::xfx(psi_x_pdf[j], osi_value_mu_fact[sd][ss]); // pdfs @ x_i
      lhapdf_result_z[j] = LHAPDF::xfx(psi_z_pdf[j], osi_value_mu_fact[sd][ss]); // pdfs @ x_i / z_i
      //      modify_pdf_content(lhapdf_result[j], oset);
      //      modify_pdf_content(lhapdf_result_z[j], oset);
      for (int i_p = 0; i_p < lhapdf_result[j].size(); i_p++)
	logger << LOG_DEBUG_VERBOSE << "psi_x_pdf[" << j << "] = " << setw(23) << setprecision(15) << psi_x_pdf[j] << "   lhapdf_result[" << j << "][" << i_p << "] = " << lhapdf_result[j][i_p] << endl;
      for (int i_p = 0; i_p < lhapdf_result_z[j].size(); i_p++)
	logger << LOG_DEBUG_VERBOSE << "psi_z_pdf[" << j << "] = " << setw(23) << setprecision(15) << psi_z_pdf[j] << "   lhapdf_result_z[" << j << "][" << i_p << "] = " << lhapdf_result_z[j][i_p] << endl;
    }

  logger << LOG_DEBUG_VERBOSE << "osi_QT_value_pdf_factor.size() = " << osi_QT_value_pdf_factor.size() << endl;
  
  for (int i = 0; i < 3; i++) {
      osi_QT_value_pdf_factor[     sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_z1x2[sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_x1z2[sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_gx2[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_x1g[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_gg[  sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_qx2[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_x1q[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_z1z2[sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_gz2[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_z1g[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_qbx2[sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_x1qb[sd][ss][i] = 0.;
      // new                       
      osi_QT_value_pdf_factor_qq[  sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_qz2[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_z1q[ sd][ss][i] = 0.;
    }

  //  for (int j = 0; j < 3; j++){cout << "osi_QT_value_pdf_factor["<< sd << "][" << ss << "][" << j << "] = " << osi_QT_value_pdf_factor[sd][ss][j] << endl;}
  for (int i = 0; i < all_pdf.size(); i++)
    {
      //    logger << LOG_DEBUG_VERBOSE << "i = " << i << "   osi_QT_value_pdf_factor.size() = " << osi_QT_value_pdf_factor.size() << endl;
      if (all_pdf[i][0] ==  1)
	{
	  osi_QT_value_pdf_factor[     sd][ss][1] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result[  2][6 + all_pdf[i][2]];
	  osi_QT_value_pdf_factor_z1x2[sd][ss][1] += lhapdf_result_z[1][6 + all_pdf[i][1]] * lhapdf_result[  2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_x1z2[sd][ss][1] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_gx2[ sd][ss][1] += lhapdf_result_z[1][6                ] * lhapdf_result[  2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_x1g[ sd][ss][1] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6                ];
	  osi_QT_value_pdf_factor_gg[  sd][ss][1] += lhapdf_result_z[1][6                ] * lhapdf_result_z[2][6                ];
	  osi_QT_value_pdf_factor_z1z2[sd][ss][1] += lhapdf_result_z[1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6 + all_pdf[i][2]];
	  osi_QT_value_pdf_factor_gz2[ sd][ss][1] += lhapdf_result_z[1][6                ] * lhapdf_result_z[2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_z1g[ sd][ss][1] += lhapdf_result_z[1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6                ];
          osi_QT_value_pdf_factor_qbx2[sd][ss][1] += lhapdf_result_z[1][6 - all_pdf[i][1]] * lhapdf_result[  2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_x1qb[sd][ss][1] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6 - all_pdf[i][2]];

	  for (int l=0; l<osi_N_f_active; l++)
            {
              osi_QT_value_pdf_factor_qx2[sd][ss][1] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) * lhapdf_result[   2][6 + all_pdf[i][2]];
              osi_QT_value_pdf_factor_x1q[sd][ss][1] +=  lhapdf_result[  1][6 + all_pdf[i][1]]              * (lhapdf_result_z[2][5-l] + lhapdf_result_z[2][7+l]); 
	      // new
	      osi_QT_value_pdf_factor_qz2[sd][ss][1] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) * lhapdf_result_z[ 2][6 + all_pdf[i][2]];
	      osi_QT_value_pdf_factor_z1q[sd][ss][1] +=  lhapdf_result_z[  1][6 + all_pdf[i][1]]            * (lhapdf_result_z[2][5-l] + lhapdf_result_z[2][7+l]);
	      for (int k=0; k<osi_N_f_active; k++)
		osi_QT_value_pdf_factor_qq[sd][ss][1] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) * (lhapdf_result_z[2][5-k] + lhapdf_result_z[2][7+k]);
            }
	}

      else if (all_pdf[i][0] == -1)
	{
	  osi_QT_value_pdf_factor[     sd][ss][2] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result[  2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_z1x2[sd][ss][2] += lhapdf_result_z[1][6 + all_pdf[i][1]] * lhapdf_result[  2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_x1z2[sd][ss][2] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_gx2[ sd][ss][2] += lhapdf_result_z[1][6                ] * lhapdf_result[  2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_x1g[ sd][ss][2] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6                ];
	  osi_QT_value_pdf_factor_gg[  sd][ss][2] += lhapdf_result_z[1][6                ] * lhapdf_result_z[2][6                ];
	  osi_QT_value_pdf_factor_z1z2[sd][ss][2] += lhapdf_result_z[1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6 + all_pdf[i][2]];
	  osi_QT_value_pdf_factor_gz2[ sd][ss][2] += lhapdf_result_z[1][6                ] * lhapdf_result_z[2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_z1g[ sd][ss][2] += lhapdf_result_z[1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6                ];
          osi_QT_value_pdf_factor_qbx2[sd][ss][2] += lhapdf_result_z[1][6 - all_pdf[i][1]] * lhapdf_result[  2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_x1qb[sd][ss][2] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6 - all_pdf[i][2]];
	  
	  for (int l=0; l<osi_N_f_active; l++)
            {
              osi_QT_value_pdf_factor_qx2[sd][ss][2] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) *  lhapdf_result[  2][6 + all_pdf[i][2]];
              osi_QT_value_pdf_factor_x1q[sd][ss][2] +=  lhapdf_result[  1][6 + all_pdf[i][1]]              * (lhapdf_result_z[2][5-l] + lhapdf_result_z[2][7+l]);
	      // new
	      osi_QT_value_pdf_factor_qz2[sd][ss][2] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) *  lhapdf_result_z[2][6 + all_pdf[i][2]];
	      osi_QT_value_pdf_factor_z1q[sd][ss][2] +=  lhapdf_result_z[1][6 + all_pdf[i][1]]              * (lhapdf_result_z[2][5-l] + lhapdf_result_z[2][7+l]);
	      for (int k=0; k<osi_N_f_active; k++)
                osi_QT_value_pdf_factor_qq[sd][ss][2] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) * (lhapdf_result_z[2][5-k] + lhapdf_result_z[2][7+k]);
	    }
	}
      
      else 
	{cout << "No valid pdf selection!" << endl; exit(1);}
    }

  osi_QT_value_pdf_factor[     sd][ss][0] = osi_QT_value_pdf_factor[     sd][ss][1] + osi_QT_value_pdf_factor[     sd][ss][2];  
  //  for (int j = 0; j < 3; j++){cout << "with x-factors: osi_QT_value_pdf_factor["<< sd << "][" << ss << "][" << j << "] = " << osi_QT_value_pdf_factor[sd][ss][j] << endl;}
  osi_QT_value_pdf_factor_z1x2[sd][ss][0] = osi_QT_value_pdf_factor_z1x2[sd][ss][1] + osi_QT_value_pdf_factor_z1x2[sd][ss][2];
  osi_QT_value_pdf_factor_x1z2[sd][ss][0] = osi_QT_value_pdf_factor_x1z2[sd][ss][1] + osi_QT_value_pdf_factor_x1z2[sd][ss][2];
  osi_QT_value_pdf_factor_gx2[ sd][ss][0] = osi_QT_value_pdf_factor_gx2[ sd][ss][1] + osi_QT_value_pdf_factor_gx2[ sd][ss][2];
  osi_QT_value_pdf_factor_x1g[ sd][ss][0] = osi_QT_value_pdf_factor_x1g[ sd][ss][1] + osi_QT_value_pdf_factor_x1g[ sd][ss][2];
  osi_QT_value_pdf_factor_gg[  sd][ss][0] = osi_QT_value_pdf_factor_gg[  sd][ss][1] + osi_QT_value_pdf_factor_gg[  sd][ss][2];
  osi_QT_value_pdf_factor_z1z2[sd][ss][0] = osi_QT_value_pdf_factor_z1z2[sd][ss][1] + osi_QT_value_pdf_factor_z1z2[sd][ss][2];
  osi_QT_value_pdf_factor_gz2[ sd][ss][0] = osi_QT_value_pdf_factor_gz2[ sd][ss][1] + osi_QT_value_pdf_factor_gz2[ sd][ss][2];
  osi_QT_value_pdf_factor_z1g[ sd][ss][0] = osi_QT_value_pdf_factor_z1g[ sd][ss][1] + osi_QT_value_pdf_factor_z1g[ sd][ss][2];
  osi_QT_value_pdf_factor_qbx2[sd][ss][0] = osi_QT_value_pdf_factor_qbx2[sd][ss][1] + osi_QT_value_pdf_factor_qbx2[sd][ss][2];
  osi_QT_value_pdf_factor_x1qb[sd][ss][0] = osi_QT_value_pdf_factor_x1qb[sd][ss][1] + osi_QT_value_pdf_factor_x1qb[sd][ss][2];
  osi_QT_value_pdf_factor_qx2[ sd][ss][0] = osi_QT_value_pdf_factor_qx2[ sd][ss][1] + osi_QT_value_pdf_factor_qx2[ sd][ss][2];
  osi_QT_value_pdf_factor_x1q[ sd][ss][0] = osi_QT_value_pdf_factor_x1q[ sd][ss][1] + osi_QT_value_pdf_factor_x1q[ sd][ss][2];
  // new                                  
  osi_QT_value_pdf_factor_qz2[ sd][ss][0] = osi_QT_value_pdf_factor_qz2[ sd][ss][1] + osi_QT_value_pdf_factor_qz2[ sd][ss][2];
  osi_QT_value_pdf_factor_z1q[ sd][ss][0] = osi_QT_value_pdf_factor_z1q[ sd][ss][1] + osi_QT_value_pdf_factor_z1q[ sd][ss][2];
  osi_QT_value_pdf_factor_qq[  sd][ss][0] = osi_QT_value_pdf_factor_qq[  sd][ss][1] + osi_QT_value_pdf_factor_qq[  sd][ss][2];

  for (int i = 0; i < 3; i++)
    {
      osi_QT_value_pdf_factor[     sd][ss][i] = osi_QT_value_pdf_factor[     sd][ss][i] /  psi_x_pdf[0];
      osi_QT_value_pdf_factor_z1x2[sd][ss][i] = osi_QT_value_pdf_factor_z1x2[sd][ss][i] / (psi_z_pdf[1]*psi_x_pdf[2]);
      osi_QT_value_pdf_factor_x1z2[sd][ss][i] = osi_QT_value_pdf_factor_x1z2[sd][ss][i] / (psi_x_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_gx2[ sd][ss][i] = osi_QT_value_pdf_factor_gx2[ sd][ss][i] / (psi_z_pdf[1]*psi_x_pdf[2]);
      osi_QT_value_pdf_factor_x1g[ sd][ss][i] = osi_QT_value_pdf_factor_x1g[ sd][ss][i] / (psi_x_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_gg[  sd][ss][i] = osi_QT_value_pdf_factor_gg[  sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_z1z2[sd][ss][i] = osi_QT_value_pdf_factor_z1z2[sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]); 
      osi_QT_value_pdf_factor_gz2[ sd][ss][i] = osi_QT_value_pdf_factor_gz2[ sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_z1g[ sd][ss][i] = osi_QT_value_pdf_factor_z1g[ sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_qbx2[sd][ss][i] = osi_QT_value_pdf_factor_qbx2[sd][ss][i] / (psi_z_pdf[1]*psi_x_pdf[2]);
      osi_QT_value_pdf_factor_x1qb[sd][ss][i] = osi_QT_value_pdf_factor_x1qb[sd][ss][i] / (psi_x_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_qx2[ sd][ss][i] = osi_QT_value_pdf_factor_qx2[ sd][ss][i] / (psi_z_pdf[1]*psi_x_pdf[2]);
      osi_QT_value_pdf_factor_x1q[ sd][ss][i] = osi_QT_value_pdf_factor_x1q[ sd][ss][i] / (psi_x_pdf[1]*psi_z_pdf[2]);
      // new
      osi_QT_value_pdf_factor_qz2[ sd][ss][i] = osi_QT_value_pdf_factor_qz2[ sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_z1q[ sd][ss][i] = osi_QT_value_pdf_factor_z1q[ sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_qq[  sd][ss][i] = osi_QT_value_pdf_factor_qq[  sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
    }

  //  cout << pdf_factor[0] << ", " << pdf_factor_gg[0] << ", " << pdf_factor_qx2[0] << ", " << pdf_factor_x1q[0] << ", " << pdf_factor_z1z2[0] << ", " << pdf_factor_gz2[0] << ", " <<  pdf_factor_z1g[0] << ", " << pdf_factor_qbx2[0] << ", " << pdf_factor_x1qb[0] << endl;  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}





void calculate_pdf_LHAPDF_Tevatron_QT(vector<vector<int> > & all_pdf, int sd, int ss, phasespace_set & psi, observable_set & oset, int order) {
  static Logger logger("calculate_pdf_LHAPDF_Tevatron_QT");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<vector<double> > lhapdf_result(3);
  vector<vector<double> > lhapdf_result_z(3);

  for (int j = 1; j < 3; j++)
    {
      //    cout << "psi_x_pdf[" << j << "] = " << psi_x_pdf[j] << "osi_value_mu_fact[" << sd << "][" << ss << "] = " << osi_value_mu_fact[sd][ss] << endl;
      lhapdf_result[j] = LHAPDF::xfx(psi_x_pdf[j], osi_value_mu_fact[sd][ss]); // pdfs @ x_i
      lhapdf_result_z[j] = LHAPDF::xfx(psi_z_pdf[j], osi_value_mu_fact[sd][ss]); // pdfs @ x_i / z_i
      //      modify_pdf_content(lhapdf_result[j], oset);
      //      modify_pdf_content(lhapdf_result_z[j], oset);
      for (int i_p = 0; i_p < lhapdf_result[j].size(); i_p++)
	logger << LOG_DEBUG_VERBOSE << "psi_x_pdf[" << j << "] = " << setw(23) << setprecision(15) << psi_x_pdf[j] << "   lhapdf_result[" << j << "][" << i_p << "] = " << lhapdf_result[j][i_p] << endl;
      for (int i_p = 0; i_p < lhapdf_result_z[j].size(); i_p++)
	logger << LOG_DEBUG_VERBOSE << "psi_z_pdf[" << j << "] = " << setw(23) << setprecision(15) << psi_z_pdf[j] << "   lhapdf_result_z[" << j << "][" << i_p << "] = " << lhapdf_result_z[j][i_p] << endl;
    }

  logger << LOG_DEBUG_VERBOSE << "osi_QT_value_pdf_factor.size() = " << osi_QT_value_pdf_factor.size() << endl;
  
  for (int i = 0; i < 3; i++) {
      osi_QT_value_pdf_factor[     sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_z1x2[sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_x1z2[sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_gx2[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_x1g[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_gg[  sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_qx2[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_x1q[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_z1z2[sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_gz2[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_z1g[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_qbx2[sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_x1qb[sd][ss][i] = 0.;
      // new                       
      osi_QT_value_pdf_factor_qq[  sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_qz2[ sd][ss][i] = 0.;
      osi_QT_value_pdf_factor_z1q[ sd][ss][i] = 0.;
    }

  //  for (int j = 0; j < 3; j++){cout << "osi_QT_value_pdf_factor["<< sd << "][" << ss << "][" << j << "] = " << osi_QT_value_pdf_factor[sd][ss][j] << endl;}
  for (int i = 0; i < all_pdf.size(); i++)
    {
      //    logger << LOG_DEBUG_VERBOSE << "i = " << i << "   osi_QT_value_pdf_factor.size() = " << osi_QT_value_pdf_factor.size() << endl;
      if (all_pdf[i][0] ==  1)
	{
	  osi_QT_value_pdf_factor[     sd][ss][1] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result[  2][6 - all_pdf[i][2]];
	  osi_QT_value_pdf_factor_z1x2[sd][ss][1] += lhapdf_result_z[1][6 + all_pdf[i][1]] * lhapdf_result[  2][6 - all_pdf[i][2]];
      osi_QT_value_pdf_factor_x1z2[sd][ss][1] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6 - all_pdf[i][2]];
      osi_QT_value_pdf_factor_gx2[ sd][ss][1] += lhapdf_result_z[1][6                ] * lhapdf_result[  2][6 - all_pdf[i][2]];
      osi_QT_value_pdf_factor_x1g[ sd][ss][1] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6                ];
	  osi_QT_value_pdf_factor_gg[  sd][ss][1] += lhapdf_result_z[1][6                ] * lhapdf_result_z[2][6                ];
	  osi_QT_value_pdf_factor_z1z2[sd][ss][1] += lhapdf_result_z[1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6 - all_pdf[i][2]];
	  osi_QT_value_pdf_factor_gz2[ sd][ss][1] += lhapdf_result_z[1][6                ] * lhapdf_result_z[2][6 - all_pdf[i][2]];
      osi_QT_value_pdf_factor_z1g[ sd][ss][1] += lhapdf_result_z[1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6                ];
      osi_QT_value_pdf_factor_qbx2[sd][ss][1] += lhapdf_result_z[1][6 - all_pdf[i][1]] * lhapdf_result[  2][6 - all_pdf[i][2]];
      osi_QT_value_pdf_factor_x1qb[sd][ss][1] += lhapdf_result[  1][6 + all_pdf[i][1]] * lhapdf_result_z[2][6 + all_pdf[i][2]];

	  for (int l=0; l<osi_N_f_active; l++)
            {
              osi_QT_value_pdf_factor_qx2[sd][ss][1] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) * lhapdf_result[   2][6 - all_pdf[i][2]];
              osi_QT_value_pdf_factor_x1q[sd][ss][1] +=  lhapdf_result[  1][6 + all_pdf[i][1]]              * (lhapdf_result_z[2][5-l] + lhapdf_result_z[2][7+l]); 
	      // new
	      osi_QT_value_pdf_factor_qz2[sd][ss][1] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) * lhapdf_result_z[ 2][6 - all_pdf[i][2]];
	      osi_QT_value_pdf_factor_z1q[sd][ss][1] +=  lhapdf_result_z[  1][6 + all_pdf[i][1]]            * (lhapdf_result_z[2][5-l] + lhapdf_result_z[2][7+l]);
	      for (int k=0; k<osi_N_f_active; k++)
		osi_QT_value_pdf_factor_qq[sd][ss][1] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) * (lhapdf_result_z[2][5-k] + lhapdf_result_z[2][7+k]);
            }
	}

      else if (all_pdf[i][0] == -1)
	{
	  osi_QT_value_pdf_factor[     sd][ss][2] += lhapdf_result[  1][6 - all_pdf[i][1]] * lhapdf_result[  2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_z1x2[sd][ss][2] += lhapdf_result_z[1][6 - all_pdf[i][1]] * lhapdf_result[  2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_x1z2[sd][ss][2] += lhapdf_result[  1][6 - all_pdf[i][1]] * lhapdf_result_z[2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_gx2[ sd][ss][2] += lhapdf_result_z[1][6                ] * lhapdf_result[  2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_x1g[ sd][ss][2] += lhapdf_result[  1][6 - all_pdf[i][1]] * lhapdf_result_z[2][6                ];
	  osi_QT_value_pdf_factor_gg[  sd][ss][2] += lhapdf_result_z[1][6                ] * lhapdf_result_z[2][6                ];
	  osi_QT_value_pdf_factor_z1z2[sd][ss][2] += lhapdf_result_z[1][6 - all_pdf[i][1]] * lhapdf_result_z[2][6 + all_pdf[i][2]];
	  osi_QT_value_pdf_factor_gz2[ sd][ss][2] += lhapdf_result_z[1][6                ] * lhapdf_result_z[2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_z1g[ sd][ss][2] += lhapdf_result_z[1][6 - all_pdf[i][1]] * lhapdf_result_z[2][6                ];
          osi_QT_value_pdf_factor_qbx2[sd][ss][2] += lhapdf_result_z[1][6 + all_pdf[i][1]] * lhapdf_result[  2][6 + all_pdf[i][2]];
          osi_QT_value_pdf_factor_x1qb[sd][ss][2] += lhapdf_result[  1][6 - all_pdf[i][1]] * lhapdf_result_z[2][6 - all_pdf[i][2]];
	  
	  for (int l=0; l<osi_N_f_active; l++)
            {
              osi_QT_value_pdf_factor_qx2[sd][ss][2] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) *  lhapdf_result[  2][6 + all_pdf[i][2]];
              osi_QT_value_pdf_factor_x1q[sd][ss][2] +=  lhapdf_result[  1][6 - all_pdf[i][1]]              * (lhapdf_result_z[2][5-l] + lhapdf_result_z[2][7+l]);
	      // new
	      osi_QT_value_pdf_factor_qz2[sd][ss][2] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) *  lhapdf_result_z[2][6 + all_pdf[i][2]];
	      osi_QT_value_pdf_factor_z1q[sd][ss][2] +=  lhapdf_result_z[1][6 - all_pdf[i][1]]              * (lhapdf_result_z[2][5-l] + lhapdf_result_z[2][7+l]);
	      for (int k=0; k<osi_N_f_active; k++)
                osi_QT_value_pdf_factor_qq[sd][ss][2] += (lhapdf_result_z[1][5-l] + lhapdf_result_z[1][7+l]) * (lhapdf_result_z[2][5-k] + lhapdf_result_z[2][7+k]);
	    }
	}
      
      else 
	{cout << "No valid pdf selection!" << endl; exit(1);}
    }

  osi_QT_value_pdf_factor[     sd][ss][0] = osi_QT_value_pdf_factor[     sd][ss][1] + osi_QT_value_pdf_factor[     sd][ss][2];  
  //  for (int j = 0; j < 3; j++){cout << "with x-factors: osi_QT_value_pdf_factor["<< sd << "][" << ss << "][" << j << "] = " << osi_QT_value_pdf_factor[sd][ss][j] << endl;}
  osi_QT_value_pdf_factor_z1x2[sd][ss][0] = osi_QT_value_pdf_factor_z1x2[sd][ss][1] + osi_QT_value_pdf_factor_z1x2[sd][ss][2];
  osi_QT_value_pdf_factor_x1z2[sd][ss][0] = osi_QT_value_pdf_factor_x1z2[sd][ss][1] + osi_QT_value_pdf_factor_x1z2[sd][ss][2];
  osi_QT_value_pdf_factor_gx2[ sd][ss][0] = osi_QT_value_pdf_factor_gx2[ sd][ss][1] + osi_QT_value_pdf_factor_gx2[ sd][ss][2];
  osi_QT_value_pdf_factor_x1g[ sd][ss][0] = osi_QT_value_pdf_factor_x1g[ sd][ss][1] + osi_QT_value_pdf_factor_x1g[ sd][ss][2];
  osi_QT_value_pdf_factor_gg[  sd][ss][0] = osi_QT_value_pdf_factor_gg[  sd][ss][1] + osi_QT_value_pdf_factor_gg[  sd][ss][2];
  osi_QT_value_pdf_factor_z1z2[sd][ss][0] = osi_QT_value_pdf_factor_z1z2[sd][ss][1] + osi_QT_value_pdf_factor_z1z2[sd][ss][2];
  osi_QT_value_pdf_factor_gz2[ sd][ss][0] = osi_QT_value_pdf_factor_gz2[ sd][ss][1] + osi_QT_value_pdf_factor_gz2[ sd][ss][2];
  osi_QT_value_pdf_factor_z1g[ sd][ss][0] = osi_QT_value_pdf_factor_z1g[ sd][ss][1] + osi_QT_value_pdf_factor_z1g[ sd][ss][2];
  osi_QT_value_pdf_factor_qbx2[sd][ss][0] = osi_QT_value_pdf_factor_qbx2[sd][ss][1] + osi_QT_value_pdf_factor_qbx2[sd][ss][2];
  osi_QT_value_pdf_factor_x1qb[sd][ss][0] = osi_QT_value_pdf_factor_x1qb[sd][ss][1] + osi_QT_value_pdf_factor_x1qb[sd][ss][2];
  osi_QT_value_pdf_factor_qx2[ sd][ss][0] = osi_QT_value_pdf_factor_qx2[ sd][ss][1] + osi_QT_value_pdf_factor_qx2[ sd][ss][2];
  osi_QT_value_pdf_factor_x1q[ sd][ss][0] = osi_QT_value_pdf_factor_x1q[ sd][ss][1] + osi_QT_value_pdf_factor_x1q[ sd][ss][2];
  // new                                  
  osi_QT_value_pdf_factor_qz2[ sd][ss][0] = osi_QT_value_pdf_factor_qz2[ sd][ss][1] + osi_QT_value_pdf_factor_qz2[ sd][ss][2];
  osi_QT_value_pdf_factor_z1q[ sd][ss][0] = osi_QT_value_pdf_factor_z1q[ sd][ss][1] + osi_QT_value_pdf_factor_z1q[ sd][ss][2];
  osi_QT_value_pdf_factor_qq[  sd][ss][0] = osi_QT_value_pdf_factor_qq[  sd][ss][1] + osi_QT_value_pdf_factor_qq[  sd][ss][2];

  for (int i = 0; i < 3; i++)
    {
      osi_QT_value_pdf_factor[     sd][ss][i] = osi_QT_value_pdf_factor[     sd][ss][i] /  psi_x_pdf[0];
      osi_QT_value_pdf_factor_z1x2[sd][ss][i] = osi_QT_value_pdf_factor_z1x2[sd][ss][i] / (psi_z_pdf[1]*psi_x_pdf[2]);
      osi_QT_value_pdf_factor_x1z2[sd][ss][i] = osi_QT_value_pdf_factor_x1z2[sd][ss][i] / (psi_x_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_gx2[ sd][ss][i] = osi_QT_value_pdf_factor_gx2[ sd][ss][i] / (psi_z_pdf[1]*psi_x_pdf[2]);
      osi_QT_value_pdf_factor_x1g[ sd][ss][i] = osi_QT_value_pdf_factor_x1g[ sd][ss][i] / (psi_x_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_gg[  sd][ss][i] = osi_QT_value_pdf_factor_gg[  sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_z1z2[sd][ss][i] = osi_QT_value_pdf_factor_z1z2[sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]); 
      osi_QT_value_pdf_factor_gz2[ sd][ss][i] = osi_QT_value_pdf_factor_gz2[ sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_z1g[ sd][ss][i] = osi_QT_value_pdf_factor_z1g[ sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_qbx2[sd][ss][i] = osi_QT_value_pdf_factor_qbx2[sd][ss][i] / (psi_z_pdf[1]*psi_x_pdf[2]);
      osi_QT_value_pdf_factor_x1qb[sd][ss][i] = osi_QT_value_pdf_factor_x1qb[sd][ss][i] / (psi_x_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_qx2[ sd][ss][i] = osi_QT_value_pdf_factor_qx2[ sd][ss][i] / (psi_z_pdf[1]*psi_x_pdf[2]);
      osi_QT_value_pdf_factor_x1q[ sd][ss][i] = osi_QT_value_pdf_factor_x1q[ sd][ss][i] / (psi_x_pdf[1]*psi_z_pdf[2]);
      // new
      osi_QT_value_pdf_factor_qz2[ sd][ss][i] = osi_QT_value_pdf_factor_qz2[ sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_z1q[ sd][ss][i] = osi_QT_value_pdf_factor_z1q[ sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
      osi_QT_value_pdf_factor_qq[  sd][ss][i] = osi_QT_value_pdf_factor_qq[  sd][ss][i] / (psi_z_pdf[1]*psi_z_pdf[2]);
    }


  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


