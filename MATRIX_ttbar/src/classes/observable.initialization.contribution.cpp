#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"
#include "../NJsubtraction/NJsubtraction.h"

void observable_set::initialization_result(){
  Logger logger("observable_set::initialization_result");
  logger << LOG_DEBUG << "started" << endl;

  n_pc = 1;
  n_ps = 1;
  n_pz = 1;

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_contribution(contribution_set & _csi){
  Logger logger("observable_set::initialization_contribution");
  logger << LOG_DEBUG << "started" << endl;

  //  csi = &_csi;

  // Could contain the calls for initialization_tree, VA, etc.
  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_tree(){
  Logger logger("observable_set::initialization_tree");
  logger << LOG_DEBUG << "started" << endl;

  n_pc = 1;
  n_ps = 1;
  n_pz = 1;

  logger << LOG_DEBUG << "finished" << endl;
}

// introduce something similar for dipoles in order to eliminate such dipoles that cannot contribute due to their parton content (happens only if photon are not considered as jets)!!!

void observable_set::initialization_VA(){
  Logger logger("observable_set::initialization_VA");
  logger << LOG_DEBUG << "started" << endl;

  n_pc = 1;
  n_ps = 1;
  n_pz = 1;

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_specific_VA(){
  Logger logger("observable_set::initialization_specific_VA");
  logger << LOG_DEBUG << "started" << endl;

  VA_ME2_cf.resize((*VA_ioperator).size());
  VA_I_ME2_cf.resize((*VA_ioperator).size());
  for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
    VA_ME2_cf[i_a].resize((*VA_ioperator)[i_a].size(), 0.);
    VA_I_ME2_cf[i_a].resize((*VA_ioperator)[i_a].size(), 0.);
  }
  
  VA_b_ME2 = 0.;
  VA_V_ME2 = 0.;
  VA_X_ME2 = 0.;
  VA_I_ME2 = 0.;
  VA_X_ME2_CV.resize(n_scales_CV, 0.);
  
  VA_DeltaUV = 0.;
  VA_DeltaIR1 = 0.;
  VA_DeltaIR2 = 0.;
  VA_delta_flag = 0;


  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_CA(){
  Logger logger("observable_set::initialization_CA");
  logger << LOG_DEBUG << "started" << endl;

  n_pc = (*CA_collinear).size(); // might be changed later if pdf_selection or pdf_disable is used !
  n_ps = 1;
  n_pz = 3;

  x_pdf.resize(3);
  z_coll.resize(3);

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_specific_CA(){//vector<vector<collinear_set> > _collinear, vector<int> _type_parton
  Logger logger("observable_set::initialization_specific_CA");
  logger << LOG_DEBUG << "started" << endl;

  //  CA_dipole_splitting.resize(4, vector<int> (3, 0));

  CA_ME2_cf.resize((*CA_collinear).size());
  for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
    CA_ME2_cf[i_a].resize((*CA_collinear)[i_a].size(), 0.);
  }

  CA_value_log_mu2_fact = value_mu_fact;
  if (value_mu_fact[0].size() != 0){
    for (int ss = 0; ss < value_mu_fact[0].size(); ss++){
      CA_value_log_mu2_fact[0][ss] = 2 * log(value_mu_fact[0][ss]);
    }
  }
  
  logger << LOG_DEBUG << "value_mu_fact.size() = " << value_mu_fact.size() << endl;
  if (value_mu_fact[0].size() != 0){
    logger << LOG_DEBUG << "value_mu_fact[0].size() = " << value_mu_fact[0].size() << endl;
    logger << LOG_DEBUG << "value_mu_fact_rel.size() = " << value_mu_fact_rel.size() << endl;
    for (int j = 0; j < value_mu_fact_rel[0].size(); j++){
      logger << LOG_DEBUG << "value_mu_fact[" << setw(2) << 0 << "][" << setw(2) << j << "] = " << setw(25) << value_mu_fact[0][j] << "   CA_value_mu2_fact[" << setw(2) << 0 << "][" << setw(2) << j << "] = " << setw(25) << CA_value_log_mu2_fact[0][j] <<  endl;
    }
  }
  logger << LOG_DEBUG << "before CA_value_ME2_KP" << endl;

  CA_value_ME2_KP.resize(value_mu_fact.size());
  CA_value_pdf_factor.resize(value_mu_fact.size());
  
  CA_value_integrand.resize(value_mu_fact.size());
  CA_value_integrand_D.resize(value_mu_fact.size());
  
  CA_sum_value_integrand.resize(value_mu_fact.size());
  CA_sum_value_integrand_D.resize(value_mu_fact.size());
  
  logger << LOG_DEBUG << "before CA_value_ME2_KP[sd]" << endl;
  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    logger << LOG_DEBUG << "sd = " << sd << endl;
    CA_value_ME2_KP[sd].resize(value_mu_fact[sd].size(), vector<vector<double> > (3, vector<double> ((*CA_collinear).size(), 0.)));
    CA_value_pdf_factor[sd].resize(value_mu_fact[sd].size(), vector<vector<vector<double> > > ((*CA_collinear).size(), vector<vector<double> > (3, vector<double> (3, 0.))));
    
    CA_value_integrand[sd].resize(value_mu_fact[sd].size(), vector<double> ((*CA_collinear).size()));
    CA_value_integrand_D[sd].resize(value_mu_fact[sd].size(), vector<vector<double> > (3, vector<double> ((*CA_collinear).size())));
    
    CA_sum_value_integrand[sd].resize(value_mu_fact[sd].size());
    CA_sum_value_integrand_D[sd].resize(value_mu_fact[sd].size(), vector<double> (3));
    logger << LOG_DEBUG << "sd = " << sd << endl;
  }
  logger << LOG_DEBUG << "after CA_value_ME2_KP[sd]" << endl;




  data_K.resize(3, vector<vector<double> > ((*CA_collinear).size()));
  for (int i_x = 0; i_x < 3; i_x++){
    for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
      data_K[i_x][i_a].resize((*CA_collinear)[i_a].size(), 0.);
    }
  }

  logger << LOG_DEBUG << "after data_K" << endl;

  value_logscale2_fact_papi.resize(max_dyn_fact + 1);
  value_data_P.resize(max_dyn_fact + 1);
  value_ME2_KP.resize(max_dyn_fact + 1);
  
  if (switch_TSV){
    logger << LOG_DEBUG << "after value_ME2_KP.resize" << endl;
    value_central_scale_fact[0] = scale_fact / prefactor_reference;
    // !!!    value_central_scale_fact[0] = scale_fact;
    value_central_logscale2_fact[0] = 2 * log(value_central_scale_fact[0]);
    // !!!    value_central_logscale2_fact[0] = 2 * log(scale_fact);
    logger << LOG_DEBUG << "after value_central_scale_fact[0]" << endl;
    
    for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){  //
      value_logscale2_fact_papi[v_sf].resize(n_scale_dyn_fact[v_sf], vector<vector<double> > (csi->type_parton[0].size(), vector<double> (csi->type_parton[0].size(),  0.)));
      //    value_logscale2_fact_papi[v_sf].resize(n_scale_dyn_fact[v_sf], vector<vector<double> > (_type_parton.size(), vector<double> (_type_parton.size(),  0.)));
      value_data_P[v_sf].resize(n_scale_dyn_fact[v_sf]);
      value_ME2_KP[v_sf].resize(n_scale_dyn_fact[v_sf]);
      for (int i_f = 0; i_f < value_relative_scale_fact[v_sf].size(); i_f++){
	value_data_P[v_sf][i_f].resize(3, vector<vector<double> > ((*CA_collinear).size())); 
	for (int i_x = 0; i_x < 3; i_x++){
	  for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
	    value_data_P[v_sf][i_f][i_x][i_a].resize((*CA_collinear)[i_a].size(), 0.); 
	  }
	}
	value_ME2_KP[v_sf][i_f].resize(3, vector<double> ((*CA_collinear).size(), 0.));   
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_CX_ncollinear(){//vector<vector<collinear_set> > & _collinear, vector<int> _type_parton
  Logger logger("observable_set::initialization_CX");
  logger << LOG_DEBUG << "started" << endl;

  int switch_collinear = ncollinear.size();
  logger << LOG_DEBUG << "switch_collinear = " << switch_collinear << endl;
  logger << LOG_DEBUG << "std::min(2, switch_collinear) = " << std::min(2, switch_collinear) << endl;

  x_pdf.resize(3);
  z_coll.resize(3);

  xz_pdf.resize(3, vector<double> (std::min(2, switch_collinear), 0.));

  CX_value_pdf_factor.resize(value_mu_fact.size());
  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    CX_value_pdf_factor[sd].resize(value_mu_fact[sd].size(), vector<vector<vector<double> > > (1, vector<vector<double> > (ncollinear.size(), vector<double> (3, 0.))));
  }

  logger << LOG_DEBUG << "switch_TSV = " << switch_TSV << endl;
  if (switch_TSV){
    //    value_list_pdf_factor_TSV[i_a][v_sf][i_m][i_l][direction]
    logger << LOG_DEBUG << "max_dyn_fact = " << max_dyn_fact << endl;
    logger << LOG_DEBUG << "n_ps = " << n_ps << endl;
    logger << LOG_DEBUG << "value_list_pdf_factor_TSV.size() = " << value_list_pdf_factor_TSV.size() << endl;
    value_list_pdf_factor_TSV.resize(n_ps, vector<vector<vector<vector<double> > > > (max_dyn_fact + 1));
    logger << LOG_DEBUG << "value_list_pdf_factor_TSV.size() = " << value_list_pdf_factor_TSV.size() << endl;

    for (int i_a = 0; i_a < n_ps; i_a++){
      for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	logger << LOG_DEBUG << "n_scale_dyn_fact.size() = " << n_scale_dyn_fact.size() << endl;
	logger << LOG_DEBUG << "n_scale_dyn_fact[" << v_sf << "] = " << n_scale_dyn_fact[v_sf] << endl;

	value_list_pdf_factor_TSV[i_a][v_sf].resize(n_scale_dyn_fact[v_sf]);
	//	value_list_pdf_factor_TSV[i_a][v_sf].resize(n_scale_dyn_fact.size());
	for (int i_m = 0; i_m < n_scale_dyn_fact[v_sf]; i_m++){
	  value_list_pdf_factor_TSV[i_a][v_sf][i_m].resize(list_combination_pdf.size(), vector<double> (3, 0.));
	}
      }
    }

    logger << LOG_DEBUG << "n_set_TSV = " << n_set_TSV << endl;
    pointer_list_pdf_factor_TSV.resize(n_ps, vector<vector<vector<vector<double*> > > > (n_set_TSV));
    logger << LOG_DEBUG << "pointer_list_pdf_factor_TSV.size() = " << pointer_list_pdf_factor_TSV.size() << endl;
    ///  max_n_scale_fact_TSV ???
    for (int i_a = 0; i_a < n_ps; i_a++){
      logger << LOG_DEBUG << "pointer_list_pdf_factor_TSV[" << i_a << "].size() = " << pointer_list_pdf_factor_TSV[i_a].size() << endl;
      for (int i_s = 0; i_s < n_set_TSV; i_s++){
	logger << LOG_DEBUG << "pointer_list_pdf_factor_TSV[" << i_a << "][" << i_s << "].size() = " << pointer_list_pdf_factor_TSV[i_a][i_s].size() << endl;
	logger << LOG_DEBUG << "n_scale_fact_TSV[" << i_s << "] = " << n_scale_fact_TSV[i_s] << endl;
	logger << LOG_DEBUG << "list_combination_pdf.size() = " << list_combination_pdf.size() << endl;
	pointer_list_pdf_factor_TSV[i_a][i_s].resize(n_scale_fact_TSV[i_s], vector<vector<double*> > (list_combination_pdf.size(), vector<double*> (3)));
	for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
	    for (int i_h = 0; i_h < 3; i_h++){
	      logger << LOG_DEBUG << "before   i_a = " << i_a << "   i_s = " << i_s << "   i_f = " << i_f << "   i_l = " << i_l << "   i_h = " << i_h << endl;
	      logger << LOG_DEBUG << "no_value_fact_TSV[" << i_s << "][" << i_f << "] = " << no_value_fact_TSV[i_s][i_f] << endl;

	      logger << LOG_DEBUG << "pointer_list_pdf_factor_TSV[" << i_a << "][" << i_s << "][" << i_f << "].size() = " << pointer_list_pdf_factor_TSV[i_a][i_s][i_f].size() << endl;
	      logger << LOG_DEBUG << "pointer_list_pdf_factor_TSV[" << i_a << "][" << i_s << "][" << i_f << "][" << i_l << "].size() = " << pointer_list_pdf_factor_TSV[i_a][i_s][i_f][i_l].size() << endl;
	      pointer_list_pdf_factor_TSV[i_a][i_s][i_f][i_l][i_h] = &value_list_pdf_factor_TSV[i_a][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_l][i_h];
	      logger << LOG_DEBUG << "after    i_a = " << i_a << "   i_s = " << i_s << "   i_f = " << i_f << "   i_l = " << i_l << "   i_h = " << i_h << endl;
	    }
	  }
	}
      }
    }
    logger << LOG_DEBUG << "n_set_TSV = " << n_set_TSV << endl;

    /*
    CX_value_integrand_TSV.resize(value_mu_fact.size());
    CX_value_integrand_TSV_D.resize(value_mu_fact.size());
    
    CX_sum_value_integrand_TSV.resize(value_mu_fact.size());
    CX_sum_value_integrand_TSV_D.resize(value_mu_fact.size());
    
    for (int sd = 0; sd < value_mu_fact.size(); sd++){
      CX_value_integrand_TSV[sd].resize(value_mu_fact[sd].size(), vector<double> (ncollinear.size()));
      CX_value_integrand_TSV_D[sd].resize(value_mu_fact[sd].size(), vector<vector<double> > (3, vector<double> (ncollinear.size(), 0.)));
      
      CX_sum_value_integrand_TSV[sd].resize(value_mu_fact[sd].size());
      CX_sum_value_integrand_TSV_D[sd].resize(value_mu_fact[sd].size(), vector<double> (3));
    }
    */
    //    CX_value_integrand_TSV_qTcut.resize(n_qTcut, vector<vector<double> > (max_dyn_fact + 1));


    value_ren_integrand_qTcut_TSV.resize(n_qTcut, vector<vector<vector<double> > > (max_dyn_ren + 1));
    for (int i_q = 0; i_q < n_qTcut; i_q++){
      for (int v_sr = 0; v_sr < max_dyn_ren + 1; v_sr++){
	//	CX_value_integrand_TSV_qTcut[i_q][v_sr].resize(n_scale_dyn_ren[v_sr]);
	value_ren_integrand_qTcut_TSV[i_q][v_sr].resize(n_scale_dyn_ren[v_sr], vector<double> (3));
      }
    }

    value_fact_integrand_qTcut_TSV.resize(n_qTcut, vector<vector<vector<double> > > (max_dyn_fact + 1));
    for (int i_q = 0; i_q < n_qTcut; i_q++){
      for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	//	CX_value_integrand_TSV_qTcut[i_q][v_sf].resize(n_scale_dyn_fact[v_sf]);
	value_fact_integrand_qTcut_TSV[i_q][v_sf].resize(n_scale_dyn_fact[v_sf], vector<double> (3));
      }
    }
    value_integrand_qTcut_TSV.resize(n_qTcut, vector<vector<vector<vector<vector<double> > > > > (max_dyn_ren + 1, vector<vector<vector<vector<double> > > > (max_dyn_fact + 1)));
    for (int i_q = 0; i_q < n_qTcut; i_q++){
      for (int v_sr = 0; v_sr < max_dyn_fact + 1; v_sr++){
	for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
	  //	CX_value_integrand_TSV_qTcut[i_q][i_v].resize(n_scale_dyn_fact[i_v]);
	  value_integrand_qTcut_TSV[i_q][v_sr][v_sf].resize(n_scale_dyn_ren[v_sr], vector<vector<double> > (n_scale_dyn_fact[v_sf], vector<double> (3, 0.)));
	}
      }
    }

    /*
    //    integrand_TSV_qTcut.resize(n_qTcut, vector<vector<double> > (n_set_TSV));
    integrand_qTcut_TSV.resize(n_qTcut, vector<vector<vector<double> > > (n_set_TSV));
    for (int i_q = 0; i_q < n_qTcut; i_q++){
      for (int i_s = 0; i_s < n_set_TSV; i_s++){
	//	integrand_TSV_qTcut[i_q][i_s].resize(n_scale_fact_TSV[i_s]);
	integrand_qTcut_TSV[i_q][i_s].resize(n_scale_fact_TSV[i_s], vector<double> (3));
      }
    }
    */

  }


  value_list_pdf_factor.resize(value_mu_fact.size());
  for (int sd = 0; sd < value_mu_fact.size(); sd++){value_list_pdf_factor[sd].resize(value_mu_fact[sd].size(), vector<vector<double> > (list_combination_pdf.size(), vector<double> (3, 0.)));}
  
  CX_value_integrand.resize(value_mu_fact.size());
  CX_value_integrand_D.resize(value_mu_fact.size());
  
  CX_sum_value_integrand.resize(value_mu_fact.size());
  CX_sum_value_integrand_D.resize(value_mu_fact.size());
  
  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    CX_value_integrand[sd].resize(value_mu_fact[sd].size(), vector<double> (ncollinear.size()));
    CX_value_integrand_D[sd].resize(value_mu_fact[sd].size(), vector<vector<double> > (3, vector<double> (ncollinear.size(), 0.)));
    
    CX_sum_value_integrand[sd].resize(value_mu_fact[sd].size());
    CX_sum_value_integrand_D[sd].resize(value_mu_fact[sd].size(), vector<double> (3));
  }
  
  CX_value_integrand_qTcut.resize(n_qTcut, vector<vector<double> > (value_mu_fact.size()));
  CX_value_integrand_qTcut_D.resize(n_qTcut, vector<vector<vector<double> > > (value_mu_fact.size()));
  for (int i_q = 0; i_q < n_qTcut; i_q++){
    for (int sd = 0; sd < value_mu_fact.size(); sd++){
      CX_value_integrand_qTcut[i_q][sd].resize(value_mu_fact[sd].size());
      CX_value_integrand_qTcut_D[i_q][sd].resize(value_mu_fact[sd].size(), vector<double> (3));
    }
  }


  // new tretment of mu_ren -- mu_fact variations in CV (closer to TSV)

  CX_value_integrand_RF.resize(value_mu_ren.size(), vector<vector<vector<vector<double> > > > (value_mu_fact.size()));
  CX_value_integrand_RF_D.resize(value_mu_ren.size(), vector<vector<vector<vector<vector<double> > > > > (value_mu_fact.size()));
  
  CX_sum_value_integrand_RF.resize(value_mu_ren.size(), vector<vector<vector<double> > > (value_mu_fact.size()));
  CX_sum_value_integrand_RF_D.resize(value_mu_ren.size(), vector<vector<vector<vector<double> > > > (value_mu_fact.size()));
  
  for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
    for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
      CX_value_integrand_RF[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<vector<double> > (value_mu_fact[i_vf].size(), vector<double> (ncollinear.size())));
      CX_value_integrand_RF_D[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<vector<vector<double> > > (value_mu_fact[i_vf].size(), vector<vector<double> > (3, vector<double> (ncollinear.size(), 0.))));
      
      CX_sum_value_integrand_RF[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<double> (value_mu_fact[i_vf].size(), 0.));
      CX_sum_value_integrand_RF_D[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<vector<double> > (value_mu_fact[i_vf].size(), vector<double> (3, 0.)));
    }
  }
  
  CX_value_integrand_RF_qTcut.resize(n_qTcut, vector<vector<vector<vector<double> > > > (value_mu_ren.size(), vector<vector<vector<double> > > (value_mu_fact.size())));
  CX_value_integrand_RF_qTcut_D.resize(n_qTcut, vector<vector<vector<vector<vector<double> > > > > (value_mu_ren.size(), vector<vector<vector<vector<double> > > > (value_mu_fact.size())));
  for (int i_q = 0; i_q < n_qTcut; i_q++){
    for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
      for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
	CX_value_integrand_RF_qTcut[i_q][i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<double> (value_mu_fact[i_vf].size(), 0.));
	CX_value_integrand_RF_qTcut_D[i_q][i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<vector<double> > (value_mu_fact[i_vf].size(), vector<double> (3, 0.)));
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_CX_multicollinear(){//vector<vector<collinear_set> > & _collinear, vector<int> _type_parton
  Logger logger("observable_set::initialization_CX");
  logger << LOG_DEBUG << "started" << endl;

  int switch_collinear = multicollinear.size();
  logger << LOG_DEBUG << "switch_collinear = " << switch_collinear << endl;
  logger << LOG_DEBUG << "std::min(2, switch_collinear) = " << std::min(2, switch_collinear) << endl;

  x_pdf.resize(3);
  z_coll.resize(3);

  xz_pdf.resize(3, vector<double> (std::min(2, switch_collinear), 0.));

  CX_value_pdf_factor.resize(value_mu_fact.size());
  
  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    CX_value_pdf_factor[sd].resize(value_mu_fact[sd].size(), vector<vector<vector<double> > > (multicollinear.size()));
    for (int ss = 0; ss < value_mu_fact[sd].size(); ss++){
      for (int i_e = 0; i_e < multicollinear.size(); i_e++){
	CX_value_pdf_factor[sd][ss][i_e].resize(multicollinear[i_e].size(), vector<double> (3, 0.));
      }
    }
  }
  
  CX_value_integrand.resize(value_mu_fact.size());
  CX_value_integrand_D.resize(value_mu_fact.size());
  
  CX_sum_value_integrand.resize(value_mu_fact.size());
  CX_sum_value_integrand_D.resize(value_mu_fact.size());
  
  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    CX_value_integrand[sd].resize(value_mu_fact[sd].size(), vector<double> (multicollinear[multicollinear.size() - 1].size()));
    CX_value_integrand_D[sd].resize(value_mu_fact[sd].size(), vector<vector<double> > (3, vector<double> (multicollinear[multicollinear.size() - 1].size())));
    
    CX_sum_value_integrand[sd].resize(value_mu_fact[sd].size());
    CX_sum_value_integrand_D[sd].resize(value_mu_fact[sd].size(), vector<double> (3));
  }
  
  CX_value_integrand_qTcut.resize(n_qTcut, vector<vector<double> > (value_mu_fact.size()));
  CX_value_integrand_qTcut_D.resize(n_qTcut, vector<vector<vector<double> > > (value_mu_fact.size()));
  for (int i_q = 0; i_q < n_qTcut; i_q++){
    for (int sd = 0; sd < value_mu_fact.size(); sd++){
      CX_value_integrand_qTcut[i_q][sd].resize(value_mu_fact[sd].size());
      CX_value_integrand_qTcut_D[i_q][sd].resize(value_mu_fact[sd].size(), vector<double> (3));
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void observable_set::initialization_RA(vector<dipole_set> & _dipole){
  Logger logger("observable_set::initialization_RA");
  logger << LOG_DEBUG << "started" << endl;

  n_ps = _dipole.size();
  n_pc = _dipole.size();
  n_pz = 1;

  RA_dipole =& _dipole;

  csi->type_parton.resize(n_ps);
  csi->basic_type_parton.resize(n_ps);
  for (int i_a = 0; i_a < (*RA_dipole).size(); i_a++){
    csi->type_parton[i_a] = (*RA_dipole)[i_a].type_parton();
    csi->basic_type_parton[i_a] = (*RA_dipole)[i_a].basic_type_parton();
    logger << LOG_DEBUG_VERBOSE << "csi->type_parton[" << i_a << "].size() = " << csi->type_parton[i_a].size() << endl;
    for (int i_p = 1; i_p < csi->type_parton[i_a].size(); i_p++){
      logger << LOG_DEBUG_VERBOSE << "csi->type_parton[" << i_a << "][" << i_p << "] = " << csi->type_parton[i_a][i_p] << endl;
    }
  }


  /*
  type_parton.resize(n_ps);
  for (int i_a = 0; i_a < (*RA_dipole).size(); i_a++){
    type_parton[i_a] = (*RA_dipole)[i_a].type_parton();
  }
  */

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_specific_RA(){
  Logger logger("observable_set::initialization_specific_RA");
  logger << LOG_DEBUG << "started" << endl;

  // only for RA contributions !!! (needed in case of n_qTcut dependence ???)
  var_A_integrand_cut.resize(n_qTcut, vector<double> (n_scales_CV));
  var_R_integrand_cut.resize(n_qTcut, vector<double> (n_scales_CV));
  var_R_integrand_cut_incl.resize(n_qTcut, vector<double> (n_scales_CV));

  A_integrand_moment_cut.resize(n_qTcut, vector<vector<double> > (n_scales_CV, vector<double> (n_moments)));
  R_integrand_moment_cut.resize(n_qTcut, vector<vector<double> > (n_scales_CV, vector<double> (n_moments)));
  R_integrand_moment_cut_incl.resize(n_qTcut, vector<vector<double> > (n_scales_CV, vector<double> (n_moments)));

  // only used for RA !!!
  RA_mu_ren.resize(n_ps);
  RA_mu_fact.resize(n_ps);
  RA_alpha_S_reference.resize(n_ps);
  RA_rel_alpha_S.resize(n_ps);
  RA_mu_ren_CV.resize(n_ps, vector<double> (n_scales_CV));
  RA_mu_fact_CV.resize(n_ps, vector<double> (n_scales_CV));
  RA_alpha_S_CV.resize(n_ps, vector<double> (n_scales_CV));
  RA_rel_alpha_S_CV.resize(n_ps, vector<double> (n_scales_CV));
  for (int i_a = 0; i_a < n_ps; i_a++){
    if (dynamic_scale == 0){  
      RA_mu_fact[i_a] = var_mu_fact;
      RA_mu_ren[i_a] = var_mu_ren;
      RA_alpha_S_reference[i_a] = var_alpha_S_reference;
      RA_rel_alpha_S[i_a] = var_rel_alpha_S;
    }
    if (dynamic_scale_CV == 0){
      RA_mu_fact_CV[i_a] = var_mu_fact_CV;
      RA_mu_ren_CV[i_a] = var_mu_ren_CV;
      RA_alpha_S_CV[i_a] = var_alpha_S_CV;
      RA_rel_alpha_S_CV[i_a] = var_rel_alpha_S_CV;
    }
  }

  RA_pdf_factor.resize(n_ps, pdf_factor);
  logger << LOG_DEBUG_VERBOSE << "n_ps = " << n_ps << endl;
  logger << LOG_DEBUG_VERBOSE << "RA_pdf_factor.size() = " << RA_pdf_factor.size() << endl;

  RA_pdf_factor_CV.resize(n_ps, pdf_factor_CV);

  RA_value_mu_fact.resize(n_ps, value_mu_fact);
  RA_value_mu_ren.resize(n_ps, value_mu_ren);
  RA_value_alpha_S.resize(n_ps, value_alpha_S);
  RA_value_factor_alpha_S.resize(n_ps, value_factor_alpha_S);


  RA_ME2.resize(n_ps, 0.);

  var_RA_ME2.resize(n_ps, 0.);
  var_RA_ME2_CV.resize(n_scales_CV, var_RA_ME2);

  RA_ME2_moment.resize(n_moments, vector<double> (n_ps, 0.));
  RA_ME2_moment_CV.resize(n_scales_CV, RA_ME2_moment);

  RA_integrand_moment.resize(n_moments);
  RA_integrand_moment_CV.resize(n_scales_CV, RA_integrand_moment);


  /*
  for (int i_a = 0; i_a < (*RA_dipole).size(); i_a++){
    csi->type_parton[i_a] = (*RA_dipole)[i_a].type_parton();
  }
  */
  for (int i_a = 1; i_a < n_ps; i_a++){
    p_parton[i_a].erase(p_parton[i_a].end() - 1, p_parton[i_a].end());
    //    csi->type_parton[i_a].erase(csi->type_parton[i_a].end() - 1, csi->type_parton[i_a].end());
    mass_parton[i_a].erase(mass_parton[i_a].end() - 1, mass_parton[i_a].end());
    mass2_parton[i_a].erase(mass2_parton[i_a].end() - 1, mass2_parton[i_a].end());
  }
  for (int i_a = 1; i_a < n_ps; i_a++){
    for (int i_p = 1; i_p < csi->type_parton[i_a].size(); i_p++){
      mass_parton[i_a][i_p] = M[abs(csi->type_parton[i_a][i_p])];
      mass2_parton[i_a][i_p] = M2[abs(csi->type_parton[i_a][i_p])];
    }
  }

  start_p_parton = p_parton;

  for (int i_a = 0; i_a < n_ps; i_a++){
    logger << LOG_DEBUG_VERBOSE << "p_parton[" << i_a << "].size() = " << p_parton[i_a].size() << endl;
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void observable_set::initialization_QT(){
  Logger logger("observable_set::initialization_QT");
  logger << LOG_DEBUG << "started" << endl;

  n_pc = 1;
  n_ps = 1;
  n_pz = 1;

  I1_int.resize(n_qTcut);
  I2_int.resize(n_qTcut);
  I3_int.resize(n_qTcut);
  I4_int.resize(n_qTcut);
  computeItildaIntegrals(min_qTcut / 100, step_qTcut / 100, n_qTcut, I1_int, I2_int, I3_int, I4_int);

  cout << "Itilde" << endl;





  // coefficients needed in several qT-subtraction contributions:
  beta0 = (33. - 2 * N_f) / 12;
  beta1 = (153. - 19 * N_f) / 24;

  Kappa = 67. / 6 - pi2_2 - 5. / 9 * N_f;

  A1g = C_A; //3.;
  B1g = -2 * beta0;
  //  Kappa = 67. / 6 - M_PI * M_PI / 2 - 5. / 9 * N_f;
  A2g = 0.5 * A1g * Kappa;


  D0gggg = 6 * beta0;
  D1gggg = 18.;
  Deltagggg = beta0 * beta0 - 3. * pi2_2;
  //  Deltagggg = beta0 * beta0 - 3. / 2 * M_PI * M_PI;

  //  Delta2gg = (9 * (8. / 3 + 3 * zeta3) - 2. / 3 * N_f - 2. * N_f) / 4;
  Delta2gg = (9 * (8. / 3 + 3 * zeta3) - 2. / 3 * N_f - 2 * N_f) / 4;
  B2g = -2 * Delta2gg + .5 * beta0 * (2. / 3 * 3 * M_PI * M_PI);// + A_F);

  H2ggD0 = -101. / 3 + 14. * N_f / 9 + 63. / 2 * zeta3;


  A1q = C_F; //4. / 3;
  B1q = -2;

  A2q = 0.5 * A1q * Kappa;

  Delta2qq = (16. / 9 * (3. / 8 - pi2_2 + 6 * zeta3) + 4 * (17. / 24 + 11. * pi2_18 - 3 * zeta3) - 2. / 3 * N_f * (1. / 6 + 2 * pi2_9)) / 4;
  //  Delta2qq = (16. / 9 * (3. / 8 - M_PI * M_PI / 2 + 6 * zeta3) + 4 * (17. / 24 + 11. * M_PI * M_PI / 18 - 3 * zeta3) - 2. / 3 * N_f * (1. / 6 + 2 * M_PI * M_PI / 9)) / 4;
  B2q = -2 * Delta2qq + 0.5 * beta0 * (2. / 3 * 4. / 3 * M_PI * M_PI);// + A_F);
  // B2q_A_F = B2q + .5 * beta0 * A_F (psp-dependent !!!)

  //  CC    Coefficients of D0 and D1 in P*P (as/pi normalization)
  D0qqqq = 8. / 3;
  D1qqqq = 32. / 9;

  //CC    Coefficients of delta(1-z) in P*P
  Deltaqqqq = 4. / 9 * (9. / 4 - 2 * pi2_3);
  //  Deltaqqqq = 4. / 9 * (9. / 4 - 2 * M_PI * M_PI / 3);

  H2qqD0 = -404. / 27 + (56. * N_f) / 81 + 14. * zeta3;
  //  cout << "1. H2qqD0 = " << H2qqD0 << endl;

  gamma2Q = (1. / 16.) * (C_F * C_A * (2 * pi2_3 - 98. / 9 - 4 * zeta3) + 40. / 9 * T_R * C_F * N_f);

  gammacusp2 = (67. / 36 - pi2_12) * C_A - 5. / 9 * T_R * N_f;




  // could be specified in terms of CT, VT, CT2, VT2 contributions ...
  A_F = 0;
  cout << "list_combination_pdf.size() = " << list_combination_pdf.size() << endl;

  coll_tH1F_contribution.resize(list_combination_pdf.size());
  coll_tH1_contribution.resize(list_combination_pdf.size());
  coll_tH1_only_H1_delta_contribution.resize(list_combination_pdf.size());
  coll_tgaga_contribution.resize(list_combination_pdf.size());
  coll_tcga_contribution.resize(list_combination_pdf.size());
  coll_tgamma2_contribution.resize(list_combination_pdf.size());
  coll_tH2_contribution.resize(list_combination_pdf.size());

  coll_tH1F.resize(list_combination_pdf.size(), 0.);
  coll_tH1.resize(list_combination_pdf.size(), 0.);
  coll_tH1_only_H1_delta.resize(list_combination_pdf.size(), 0.);
  coll_tgaga.resize(list_combination_pdf.size(), 0.);
  coll_tcga.resize(list_combination_pdf.size(), 0.);
  coll_tgamma2.resize(list_combination_pdf.size(), 0.);
  coll_tH2.resize(list_combination_pdf.size(), 0.);


  /*
  // TSV
  cout << "max_dyn_ren = " << max_dyn_ren << endl;
  cout << "max_dyn_fact = " << max_dyn_fact << endl;

  value_LR_TSV.resize(max_dyn_ren + 1);
  value_LF_TSV.resize(max_dyn_fact + 1);
  value_LQ_TSV.resize(1);

  value_sig11_TSV.resize(max_dyn_fact + 1);
  value_tH1F_TSV.resize(max_dyn_fact + 1);
  value_tH1_TSV.resize(max_dyn_fact + 1); 
  value_tH1_only_H1_delta_TSV.resize(max_dyn_fact + 1); 
  
  value_tgaga_TSV.resize(max_dyn_fact + 1);
  value_tcga_TSV.resize(max_dyn_fact + 1);
  value_tgamma2_TSV.resize(max_dyn_fact + 1);
  
  value_H1full_TSV.resize(max_dyn_ren + 1, vector<vector<vector<vector<double> > > > (max_dyn_fact + 1)); 
  

  value_tH2_TSV.resize(max_dyn_fact + 1);
  
  for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
    cout << "n_scale_dyn_ren[" << i_vr << "] = " << n_scale_dyn_ren[i_vr] << endl;
    value_LR_TSV[i_vr].resize(n_scale_dyn_ren[i_vr], 0.);
  }
  
  for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
    value_LF_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], 0.);
    
    value_tH1F_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
    value_tH1_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
    value_tH1_only_H1_delta_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
    value_sig11_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
    value_tgaga_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
    value_tcga_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
    value_tgamma2_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
    
    value_tH2_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
  }
  
  for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
    for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
      value_H1full_TSV[i_vr][i_vf].resize(n_scale_dyn_ren[i_vr], vector<vector<double> > (n_scale_dyn_fact[i_vf], vector<double> (3, 0.)));
    }
  }
  */


  //Maybe shift also I1,2,3,4_int here ???
  /* 
  //  from CT
  static int order = 1;
  static vector<double> I1_int;
  static vector<double> I2_int;
  static vector<double> I3_int;
  static vector<double> I4_int;
  if (initialization == 1){
    I1_int.resize(n_qTcut);
    I2_int.resize(n_qTcut);
    I3_int.resize(n_qTcut);
    I4_int.resize(n_qTcut);
    if (!switch_resummation) {
      computeItildaIntegrals(min_qTcut / 100, step_qTcut / 100, n_qTcut, I1_int, I2_int, I3_int, I4_int);
    }

    initialization = 0;
  }
  //  from CT2
  static int order = 2;
  static vector<double> I1_int;
  static vector<double> I2_int;
  static vector<double> I3_int;
  static vector<double> I4_int;
  if (initialization == 1){
    I1_int.resize(n_qTcut);
    I2_int.resize(n_qTcut);
    I3_int.resize(n_qTcut);
    I4_int.resize(n_qTcut);
    computeItildaIntegrals(min_qTcut / 100, step_qTcut / 100, n_qTcut, I1_int, I2_int, I3_int, I4_int);

    initialization = 0;
  }
  */

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_NJ(){
  Logger logger("observable_set::initialization_NJ");
  logger << LOG_DEBUG << "started" << endl;

  NJ.xi = 1.; // xi: N-jettiness subtraction scheme variable
  NJ.nf = N_f; // xi: N-jettiness subtraction scheme variable
  NJ.beta0 = (33. - 2 * N_f) / 12;
  NJ.beta1 = (153. - 19 * N_f) / 24;

  // WHY ??? 
  n_pc = 1;
  n_ps = 1;
  n_pz = 1;
  // WHY ??? 

  // compute integrated N-jettiness amplitudes once and for all
  LL0_int.resize(n_qTcut);
  LL1_int.resize(n_qTcut);
  LL2_int.resize(n_qTcut);
  LL3_int.resize(n_qTcut);

  NJ.computeLLnIntegrals(min_qTcut / 100, step_qTcut / 100, n_qTcut, LL0_int, LL1_int, LL2_int, LL3_int);

  logger << LOG_DEBUG << "LL0_int[0] = " << LL0_int[0] << endl;
  logger << LOG_DEBUG << "LL1_int[0] = " << LL1_int[0] << endl;
  logger << LOG_DEBUG << "LL2_int[0] = " << LL2_int[0] << endl;
  logger << LOG_DEBUG << "LL3_int[0] = " << LL3_int[0] << endl;

  logger << LOG_DEBUG << "finished" << endl;
}


void observable_set::initialization_QT_coefficients(){
  Logger logger("observable_set::initialization_QT_coefficients");
  logger << LOG_DEBUG << "started" << endl;

  // better nomix NJ and QT subtractions here ???
  // added for NJsubtraction
  VA_X_ME2_vr_mr.resize(value_mu_ren.size());
  for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
    VA_X_ME2_vr_mr[i_vr].resize(value_mu_ren[i_vr].size(), 0.);
  }


  QT_initialstate_type = 0;
  string initial_state = name_process.substr(0,4);
  if (csi->type_parton[0][1] == 0 && csi->type_parton[0][2] == 0){
    //  if (initial_state == "gg_") {
    QT_initialstate_type = 0;
  }
  else if (initial_state == "uu~_") {
    QT_initialstate_type = 1;
  }
  else if (initial_state == "dd~_") {
    QT_initialstate_type = 2;
  }
  else {
    //    assert(false);
  }

  vector<string> pa_name(csi->type_parton[0].size(), "");
  if (csi->type_parton[0][1] >= -10 && csi->type_parton[0][1] <= 10){pa_name[1] = "a";}
  if (csi->type_parton[0][2] >= -10 && csi->type_parton[0][2] <= 10){pa_name[2] = "b";}
  for (int i_p = 3; i_p < pa_name.size(); i_p++){
     if (csi->type_parton[0][i_p] >= -10 && csi->type_parton[0][i_p] <= 10){pa_name[i_p] = char(105 + i_p - 3);}
  } 

  QT_finalstate_massive_coloured  = 0;
  for (int i_p = 3; i_p < pa_name.size(); i_p++){
    if (csi->type_parton[0][i_p] >= -10 && csi->type_parton[0][i_p] <= 10 && M[abs(csi->type_parton[0][i_p])] >= 0.){
      QT_finalstate_massive_coloured = 1;
    }
  } 

  if (QT_finalstate_massive_coloured){
    determine_correlationoperator_QCD();
    QT_ME2_cf.resize(QT_correlationoperator.size(), 0.);
    for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
      logger << LOG_INFO << "QT_correlationoperator[" << i_c << "]:   emitter = " << QT_correlationoperator[i_c].no_emitter << "   spectator = " << QT_correlationoperator[i_c].no_spectator << "   no_BLHA_entry = " << QT_correlationoperator[i_c].no_BLHA_entry << endl;
    }

    // only needed for qqx and gg channels
    QT_ME2_loopcf.resize(QT_correlationoperator.size(), 0.);
    for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
      logger << LOG_INFO << "loop x Born:   QT_correlationoperator[" << i_c << "]:   emitter = " << QT_correlationoperator[i_c].no_emitter << "   spectator = " << QT_correlationoperator[i_c].no_spectator << "   no_BLHA_entry = " << QT_correlationoperator[i_c].no_BLHA_entry << endl;
    }
    /*
    QT_ME2_Born4cf.resize(QT_correlationoperator.size(), 0.);
    for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
      logger << LOG_INFO << "loop x Born:   QT_correlationoperator[" << i_c << "]:   emitter = " << QT_correlationoperator[i_c].no_emitter << "   spectator = " << QT_correlationoperator[i_c].no_spectator << "   no_BLHA_entry = " << QT_correlationoperator[i_c].no_BLHA_entry << endl;
    }
    */
    Ft1born_4correlator.resize(QT_correlationoperator.size(), 0.);


  }




 
  // central - CV
  value_LR.resize(value_mu_ren.size());
  value_LF.resize(value_mu_fact.size());
  //  value_LQ.resize(1);
 
  value_sig11.resize(value_mu_fact.size());
  value_tH1F.resize(value_mu_fact.size());
  value_tH1.resize(value_mu_fact.size());
  value_tH1_only_H1_delta.resize(value_mu_fact.size());
  value_H1full.resize(value_mu_ren.size(), vector<vector<vector<double> > > (value_mu_fact.size()));
  value_tgaga.resize(value_mu_fact.size());
  value_tcga.resize(value_mu_fact.size());
  value_tgamma2.resize(value_mu_fact.size());

  value_sig12.resize(value_mu_fact.size());
  value_sig23.resize(value_mu_fact.size());
  value_sig24.resize(value_mu_fact.size());
  value_sig21.resize(value_mu_ren.size(), vector<vector<vector<double> > > (value_mu_fact.size()));
  value_sig22.resize(value_mu_ren.size(), vector<vector<vector<double> > > (value_mu_fact.size()));

  value_tH2.resize(value_mu_fact.size());

  coll_tH1F_pdf.resize(value_mu_fact.size());
  coll_tH1_pdf.resize(value_mu_fact.size());
  coll_tH1_only_H1_delta_pdf.resize(value_mu_fact.size());
  coll_tgaga_pdf.resize(value_mu_fact.size());
  coll_tcga_pdf.resize(value_mu_fact.size());
  coll_tgamma2_pdf.resize(value_mu_fact.size());
  coll_tH2_pdf.resize(value_mu_fact.size());

  /*
  for (int i_vq = 0; i_vq < 1; i_vq++){
    value_LQ[i_vq].resize(1, 0.);
  }
  */

  for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
    value_LR[i_vr].resize(value_mu_ren[i_vr].size(), 0.);
  }

  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    value_LF[i_vf].resize(value_mu_fact[i_vf].size(), 0.);

    value_tH1F[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    value_tH1[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    value_tH1_only_H1_delta[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
 
    value_sig11[i_vf].resize(value_mu_fact[i_vf].size(), 0.);

    value_tgaga[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    value_tcga[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    value_tgamma2[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    //  }

    //  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    value_tH2[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    //  }
    coll_tH1F_pdf[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (list_combination_pdf.size(), 0.));
    coll_tH1_pdf[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (list_combination_pdf.size(), 0.));
    coll_tH1_only_H1_delta_pdf[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (list_combination_pdf.size(), 0.));
    coll_tgaga_pdf[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (list_combination_pdf.size(), 0.));
    coll_tcga_pdf[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (list_combination_pdf.size(), 0.));
    coll_tgamma2_pdf[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (list_combination_pdf.size(), 0.));

    value_sig12[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    value_sig23[i_vf].resize(value_mu_fact[i_vf].size(), 0.);
    value_sig24[i_vf].resize(value_mu_fact[i_vf].size(), 0.);

    coll_tH2_pdf[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (list_combination_pdf.size(), 0.));

    //  for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
    for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
      value_H1full[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<double> (value_mu_fact[i_vf].size(), 0.));

      value_sig21[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<double> (value_mu_fact[i_vf].size(), 0.));
      value_sig22[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<double> (value_mu_fact[i_vf].size(), 0.));
    }
  }

  if (switch_distribution){
    // D
    value_sig11_D.resize(value_mu_fact.size());
    value_tH1F_D.resize(value_mu_fact.size());
    value_tH1_D.resize(value_mu_fact.size());
    value_tH1_only_H1_delta_D.resize(value_mu_fact.size());
    value_H1full_D.resize(value_mu_ren.size(), vector<vector<vector<vector<double> > > > (value_mu_fact.size()));
    value_tgaga_D.resize(value_mu_fact.size());
    value_tcga_D.resize(value_mu_fact.size());
    value_tgamma2_D.resize(value_mu_fact.size());

    value_tH2_D.resize(value_mu_fact.size());

    coll_tH1F_pdf_D.resize(value_mu_fact.size());
    coll_tH1_pdf_D.resize(value_mu_fact.size());
    coll_tH1_only_H1_delta_pdf_D.resize(value_mu_fact.size());
    coll_tgaga_pdf_D.resize(value_mu_fact.size());
    coll_tcga_pdf_D.resize(value_mu_fact.size());
    coll_tgamma2_pdf_D.resize(value_mu_fact.size());
    coll_tH2_pdf_D.resize(value_mu_fact.size());

    value_sig12_D.resize(value_mu_fact.size());
    value_sig23_D.resize(value_mu_fact.size());
    value_sig24_D.resize(value_mu_fact.size());

    value_sig21_D.resize(value_mu_ren.size(), vector<vector<vector<vector<double> > > > (value_mu_fact.size()));
    value_sig22_D.resize(value_mu_ren.size(), vector<vector<vector<vector<double> > > > (value_mu_fact.size()));

    for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
      value_tH1F_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (3, 0.));
      value_tH1_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (3, 0.));
      value_tH1_only_H1_delta_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (3, 0.));
      
      value_sig11_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (3, 0.));
      
      value_tgaga_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (3, 0.));
      value_tcga_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (3, 0.));
      value_tgamma2_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (3, 0.));
      //      for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
      value_tH2_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (3, 0.));

      coll_tH1F_pdf_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
      coll_tH1_pdf_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
      coll_tH1_only_H1_delta_pdf_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
      coll_tgaga_pdf_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
      coll_tcga_pdf_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
      coll_tgamma2_pdf_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));

      value_sig12_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (3, 0.));
      value_sig23_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (3, 0.));
      value_sig24_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<double> (3, 0.));

      coll_tH2_pdf_D[i_vf].resize(value_mu_fact[i_vf].size(), vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));

	//      }
      //      for (int i_vf = 0; i_vf < value_mu_fact.size(); i_vf++){
      for (int i_vr = 0; i_vr < value_mu_ren.size(); i_vr++){
	value_H1full_D[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<vector<double> > (value_mu_fact[i_vf].size(), vector<double> (3, 0.)));

	value_sig21_D[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<vector<double> > (value_mu_fact[i_vf].size(), vector<double> (3, 0.)));
	value_sig22_D[i_vr][i_vf].resize(value_mu_ren[i_vr].size(), vector<vector<double> > (value_mu_fact[i_vf].size(), vector<double> (3, 0.)));
      }
    }
  }

  if (switch_TSV){
    logger << LOG_DEBUG << "max_dyn_ren = " << max_dyn_ren << endl;
    logger << LOG_DEBUG << "max_dyn_fact = " << max_dyn_fact << endl;
    
    value_LR_TSV.resize(max_dyn_ren + 1);
    value_LF_TSV.resize(max_dyn_fact + 1);
    value_LQ_TSV.resize(1);

    value_sig11_TSV.resize(max_dyn_fact + 1);
    value_tH1F_TSV.resize(max_dyn_fact + 1);
    value_tH1_TSV.resize(max_dyn_fact + 1); 
    value_tH1_only_H1_delta_TSV.resize(max_dyn_fact + 1); 
    
    value_tgaga_TSV.resize(max_dyn_fact + 1);
    value_tcga_TSV.resize(max_dyn_fact + 1);
    value_tgamma2_TSV.resize(max_dyn_fact + 1);
    
    value_H1full_TSV.resize(max_dyn_ren + 1, vector<vector<vector<vector<double> > > > (max_dyn_fact + 1)); 
    
    value_tH2_TSV.resize(max_dyn_fact + 1);
    
    coll_tH1F_pdf_TSV.resize(max_dyn_fact + 1);
    coll_tH1_pdf_TSV.resize(max_dyn_fact + 1);
    coll_tH1_only_H1_delta_pdf_TSV.resize(max_dyn_fact + 1);
    coll_tgaga_pdf_TSV.resize(max_dyn_fact + 1);
    coll_tcga_pdf_TSV.resize(max_dyn_fact + 1);
    coll_tgamma2_pdf_TSV.resize(max_dyn_fact + 1);

    value_sig12_TSV.resize(max_dyn_fact + 1);
    value_sig23_TSV.resize(max_dyn_fact + 1);
    value_sig24_TSV.resize(max_dyn_fact + 1);

    value_sig21_TSV.resize(max_dyn_ren + 1, vector<vector<vector<vector<double> > > > (max_dyn_fact + 1));
    value_sig22_TSV.resize(max_dyn_ren + 1, vector<vector<vector<vector<double> > > > (max_dyn_fact + 1));

    coll_tH2_pdf_TSV.resize(max_dyn_fact + 1);
    
    for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
      logger << LOG_DEBUG << "n_scale_dyn_ren[" << i_vr << "] = " << n_scale_dyn_ren[i_vr] << endl;
      value_LR_TSV[i_vr].resize(n_scale_dyn_ren[i_vr], 0.);
    }
    
    for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
      value_LF_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], 0.);
      
      value_tH1F_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
      value_tH1_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
      value_tH1_only_H1_delta_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
      value_sig11_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));

      value_tgaga_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
      value_tcga_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
      value_tgamma2_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
      
      value_sig12_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
      value_sig23_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
      value_sig24_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));



      value_tH2_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<double> (3, 0.));
      
      coll_tH1F_pdf_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
      coll_tH1_pdf_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
      coll_tH1_only_H1_delta_pdf_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
      coll_tgaga_pdf_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
      coll_tcga_pdf_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
      coll_tgamma2_pdf_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
      coll_tH2_pdf_TSV[i_vf].resize(n_scale_dyn_fact[i_vf], vector<vector<double> > (3, vector<double> (list_combination_pdf.size(), 0.)));
    }
    
    for (int i_vr = 0; i_vr < max_dyn_ren + 1; i_vr++){
      for (int i_vf = 0; i_vf < max_dyn_fact + 1; i_vf++){
	value_H1full_TSV[i_vr][i_vf].resize(n_scale_dyn_ren[i_vr], vector<vector<double> > (n_scale_dyn_fact[i_vf], vector<double> (3, 0.)));

	value_sig21_TSV[i_vr][i_vf].resize(n_scale_dyn_ren[i_vr], vector<vector<double> > (n_scale_dyn_fact[i_vf], vector<double> (3, 0.)));
	value_sig22_TSV[i_vr][i_vf].resize(n_scale_dyn_ren[i_vr], vector<vector<double> > (n_scale_dyn_fact[i_vf], vector<double> (3, 0.)));

      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_specific_VT(){
  Logger logger("observable_set::initialization_specific_VT");
  logger << LOG_DEBUG << "started" << endl;

  VA_b_ME2 = 0.;
  VA_V_ME2 = 0.;
  VA_X_ME2 = 0.;
  VA_X_ME2_CV.resize(n_scales_CV, 0.);
  
  VA_DeltaUV = 0.;
  VA_DeltaIR1 = 0.;
  VA_DeltaIR2 = 0.;
  VA_delta_flag = 0;

  if (switch_resummation) {
#ifdef MORE
    if (type_contribution == "NLL_LO" || type_contribution == "NLL_NLO") {
      resorder_.resorder=1;
    } else {
      resorder_.resorder=2;
    }
    if (type_contribution == "NLL_LO" || type_contribution == "NNLL_LO") {
      resorder_.pertorder=0;
    } else if (type_contribution == "NLL_NLO" || type_contribution == "NNLL_NLO") {
      resorder_.pertorder=1;
    } else {
      resorder_.pertorder=2;
    }
    // FIXME: not available here, probably should introduce new switch
    //resorder_.approxNNLL=user_switch_value()[user_switch_map["approx_NNLL"]];
    if (resorder_.approxNNLL) {
      cout << "computing approximated NNLL" << endl;
    }

    resumchan_.pdf_content_modify=-1;//pdf_content_modify;

    lopowers_.pcf=0;
    lopowers_.kcf=0;

    flaginitres_.doinitres=1;

    char test[64]="MSTW2008nnlo90cl.LHgrid               "; //FIXME
    for (int i=0; i<LHAPDFname.size(); i++) {
      test[i]=LHAPDFname.c_str()[i];
  //    cout << i << ": " << test[i] << endl;
    }
    for (int i=LHAPDFname.size(); i<63; i++) {
      test[i]=' ';
  //    cout << i << ": " << test[i] << endl;
    }
    char fstring[64];
    ConvertToFortran(fstring,sizeof test,test);

    for (int i=0; i<64; i++) {
      pdfhres_.name[i]=fstring[i];
  //    cout << i << ": " << pdfhres_.name[i] << endl;
    }

    pdfhres_.member=1;
    pdfhres_.ih1=1;
    pdfhres_.ih2=1;

    nf_.nf = N_f;

    pdfset_();

    // Constants and resummation coefficients are initialized
  //  init_const_and_coeff_();

    // should not enter anywhere
    masses_.zmass=0;//M_Z;
    masses_.mtop=0;//M_t;
    masses_.mbot=0;//M_b;

    somescale_.startscale=var_mu_fact;
#else
    assert(false);
#endif
  }

  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::initialization_specific_CT(phasespace_set & psi){
  Logger logger("observable_set::initialization_specific_CT");

  QT_Qres = user.cut_value[user.cut_map["Qres"]];
  switch_resummation = user.switch_value[user.switch_map["do_resummation"]];
  dynamical_Qres = user.switch_value[user.switch_map["dynamical_Qres"]];
  QT_Qres_prefactor = user.cut_value[user.cut_map["Qres_prefactor"]];
  if (QT_Qres_prefactor==0) {
    QT_Qres_prefactor=1;
  }
  psi.do_resummation = switch_resummation;
  psi.Qres = QT_Qres;
  psi.dynamical_Qres = dynamical_Qres;
  psi.Qres_prefactor = QT_Qres_prefactor;

  logger << LOG_DEBUG << "finished" << endl;
}
