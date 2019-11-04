#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

void observable_set::initialization_specific_QT(){
  Logger logger("observable_set::initialization_specific_QT");
  logger << LOG_DEBUG << "started" << endl;

  QT_pdf_factor_z1x2.resize(3, 0.);
  QT_pdf_factor_x1z2.resize(3, 0.);
  QT_pdf_factor_gx2.resize( 3, 0.);
  QT_pdf_factor_x1g.resize( 3, 0.);
  QT_pdf_factor_gg.resize(  3, 0.);
  QT_pdf_factor_qx2.resize( 3, 0.);
  QT_pdf_factor_x1q.resize( 3, 0.);
  QT_pdf_factor_z1z2.resize(3, 0.);
  QT_pdf_factor_gz2.resize( 3, 0.);
  QT_pdf_factor_z1g.resize( 3, 0.);
  QT_pdf_factor_qbx2.resize(3, 0.);
  QT_pdf_factor_x1qb.resize(3, 0.);
  // new
  QT_pdf_factor_qq.resize(  3, 0.);
  QT_pdf_factor_qz2.resize( 3, 0.);
  QT_pdf_factor_z1q.resize( 3, 0.);


  QT_pdf_factor_z1x2_CV.resize(n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_x1z2_CV.resize(n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_gx2_CV.resize( n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_x1g_CV.resize( n_scales_CV, vector<double> (3, 0.)); 
  QT_pdf_factor_gg_CV.resize(  n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_qx2_CV.resize( n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_x1q_CV.resize( n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_z1z2_CV.resize(n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_gz2_CV.resize( n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_z1g_CV.resize( n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_qbx2_CV.resize(n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_x1qb_CV.resize(n_scales_CV, vector<double> (3, 0.));
  // new                       
  QT_pdf_factor_qq_CV.resize(  n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_qz2_CV.resize( n_scales_CV, vector<double> (3, 0.));
  QT_pdf_factor_z1q_CV.resize( n_scales_CV, vector<double> (3, 0.));

  QT_sig_CV.resize(     n_scales_CV, 0.);
  QT_virt_CV.resize(    n_scales_CV, 0.);
  QT_sigma_qT.resize(   n_qTcut    , 0.);
  QT_sigma_qT_CV.resize(n_scales_CV, vector<double>(n_qTcut, 0.));

  QT_A0       = 0.;
  QT_A1       = 0.;
  QT_A2       = 0.;
  QT_H1_delta = 0;
  QT_H2_delta = 0;


  QT_value_pdf_factor.resize(     value_mu_fact.size());
  QT_value_pdf_factor_z1x2.resize(value_mu_fact.size());
  QT_value_pdf_factor_x1z2.resize(value_mu_fact.size());
  QT_value_pdf_factor_gx2.resize( value_mu_fact.size());
  QT_value_pdf_factor_x1g.resize( value_mu_fact.size());
  QT_value_pdf_factor_gg.resize(  value_mu_fact.size());
  QT_value_pdf_factor_qx2.resize( value_mu_fact.size());
  QT_value_pdf_factor_x1q.resize( value_mu_fact.size());
  QT_value_pdf_factor_z1z2.resize(value_mu_fact.size());
  QT_value_pdf_factor_gz2.resize( value_mu_fact.size());
  QT_value_pdf_factor_z1g.resize( value_mu_fact.size());
  QT_value_pdf_factor_qbx2.resize(value_mu_fact.size());
  QT_value_pdf_factor_x1qb.resize(value_mu_fact.size());
  // new
  QT_value_pdf_factor_qq.resize(  value_mu_fact.size());
  QT_value_pdf_factor_qz2.resize( value_mu_fact.size());
  QT_value_pdf_factor_z1q.resize( value_mu_fact.size());

  for (int sd = 0; sd < value_mu_fact.size(); sd++)
    {
      QT_value_pdf_factor[     sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));      
      QT_value_pdf_factor_z1x2[sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_x1z2[sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_gx2[ sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_x1g[ sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_gg[  sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_qx2[ sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_x1q[ sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_z1z2[sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_gz2[ sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_z1g[ sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_qbx2[sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_x1qb[sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      // new                   
      QT_value_pdf_factor_qq[  sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_qz2[ sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));
      QT_value_pdf_factor_z1q[ sd].resize(value_mu_fact[sd].size(), vector<double> (3, 0.));

      for (int ss = 0; ss < value_mu_fact[sd].size(); ss++)
	{
	  if (coll_choice == 0)
	    {
	      QT_value_pdf_factor[sd][ss][1] = 1.; 
	      QT_value_pdf_factor[sd][ss][2] = 0.; 
	      QT_value_pdf_factor[sd][ss][0] = QT_value_pdf_factor[sd][ss][1];
	    }
	}
    }

  QT_value_integrand_qTcut.resize(n_qTcut, vector<vector<double> > (value_mu_fact.size()));
  QT_value_integrand_qTcut_D.resize(n_qTcut, vector<vector<vector<double> > > (value_mu_fact.size()));
  for (int i_q = 0; i_q < n_qTcut; i_q++){
    for (int sd = 0; sd < value_mu_fact.size(); sd++){
      QT_value_integrand_qTcut[i_q][sd].resize(value_mu_fact[sd].size());
      QT_value_integrand_qTcut_D[i_q][sd].resize(value_mu_fact[sd].size(), vector<double> (3));
    }
  }
  
  logger << LOG_DEBUG << "finished" << endl;
}

