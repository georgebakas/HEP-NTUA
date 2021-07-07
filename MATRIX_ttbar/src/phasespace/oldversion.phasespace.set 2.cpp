#include "../include/classes.cxx"

void phasespace_set::calculate_IS_CT(){
  Logger logger("phasespace_set::calculate_IS_CT");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  //  cout << "QT_g_IS_ = " << setw(23) << setprecision(15) << QT_g_IS_ << endl;
  //  cout << "QT_eps = " << setw(23) << setprecision(15) << QT_eps << endl;
  
  //  double QT_random_qt2;
  
  QT_random_z1->get_random(QT_random_IS,QT_g_IS_);
  zz_pdf[1]=pow(x_pdf[1],QT_eps+(1-QT_eps)*QT_random_IS);
  QT_g_z1=log(x_pdf[1])*(1-QT_eps);
  //  cout << "z1: QT_random_IS = " << setw(23) << setprecision(15) << QT_random_IS << "   zz_pdf[1] = " << setw(23) << setprecision(15) << zz_pdf[1] << "   QT_g_z1 = " << setw(23) << setprecision(15) << QT_g_z1 << "   QT_g_IS_ = " << setw(23) << setprecision(15) << QT_g_IS_ << endl;

  //  cout << "g_pdf = " << g_pdf << endl;

  g_pdf*=QT_g_IS_;
  
  QT_random_z2->get_random(QT_random_IS,QT_g_IS_);
  zz_pdf[2]=pow(x_pdf[2],QT_eps+(1-QT_eps)*QT_random_IS);
  QT_g_z2=log(x_pdf[2])*(1-QT_eps);
  //  cout << "z2: QT_random_IS = " << setw(23) << setprecision(15) << QT_random_IS << "   zz_pdf[2] = " << setw(23) << setprecision(15) << zz_pdf[2] << "   QT_g_z2 = " << setw(23) << setprecision(15) << QT_g_z2 << "   QT_g_IS_ = " << setw(23) << setprecision(15) << QT_g_IS_ << endl;
  g_pdf*=QT_g_IS_;
  
  QT_random_qt->get_random(QT_random_qt2,QT_g_IS_);
  //  cout << "qt: QT_random_qt2 = " << setw(23) << setprecision(15) << QT_random_qt2 << "   QT_g_IS_ = " << setw(23) << setprecision(15) << QT_g_IS_ << endl;
  g_pdf*=QT_g_IS_;
  
  // for pdf calculation
  z_pdf[1]=x_pdf[1]/zz_pdf[1];
  z_pdf[2]=x_pdf[2]/zz_pdf[2];
  z_pdf[0]=z_pdf[1]*z_pdf[2];
  
  // the (Bessel transformed) large logarithms vanish up to double precision accuracy for arguments larger than 40 -> use this as upper bound
  // 1.261 ~ b0*b0
  double maxqt2 = 0.;
  double Qres_tmp=0.0, Qres_lower, Qres_upper;
  
  // Qres appears at two places here: as the reference value for the lower and for the upper qT cut
  // conceptually, these are partially unrelated
  // in fixed order computations, we have to use Qres as the lower cut off scale, because we precompute the qT integrals
  // in resummed computations, we explicitly integrate over qT, so this is unnecessary; in fact
  // it is more efficient to use mF as the lower reference scale, as qT/mF determine the universal limit
  // the upper cut off scale should alway be given by Qres, as otherwise qT/Qres can become huge, resulting in numerical
  // problems in the counterterm weight computation
  
  // Qres != mF does not work at the moment; has to be implemented in process specific cut files
    assert(Qres==0 || do_resummation);
    
    if (dynamical_Qres) {
      Qres_tmp = Qres_prefactor*xbsqrts_all[0][0];
    } else {
      Qres_tmp = Qres;
    }
    if (Qres_tmp == 0) {
      // default choice
      Qres_tmp = xbsqrts_all[0][0];
    }
    
    Qres_upper = Qres_tmp;
    
    if (do_resummation){
      Qres_lower = xbsqrts_all[0][0];
    } else {
      Qres_lower = Qres_tmp;
    }
    
    if (contribution_order_alpha_s[0] == 1){maxqt2 = 1.261 * 30. * 20. * Qres_upper * Qres_upper;}
    else if (contribution_order_alpha_s[0] == 2){maxqt2 = 1.261 * 40. * 50. * Qres_upper * Qres_upper;}
    
    if (switch_qTcut == 1){
      double r0 = min_qTcut / 100.;
      //  version with cut on pT/Qres_cut
      QT_qt2 = r0 * r0 * pow(Qres_lower,2) * exp(log(maxqt2 / r0 / r0 / Qres_lower / Qres_lower) * QT_random_qt2);
//       cout << "qT=" << QT_qt2 << ", Qres=" << Qres_upper << ", " << Qres_lower << ", ratio=" << QT_qt2/Qres_upper << endl;
      QT_jacqt2 = log(maxqt2 / r0 / r0 / Qres_lower / Qres_lower) * QT_qt2;
      //  version with cut on pT/m_inv
      //    QT_qt2 = r0 * r0 * xbs_all[0][0] * exp(log(maxqt2 / r0 / r0 / xbs_all[0][0]) * QT_random_qt2);
      //    QT_jacqt2 = log(maxqt2 / r0 / r0 / xbs_all[0][0]) * QT_qt2;
    }
    else if (switch_qTcut == 2){
      // old version with cut on pT only (not relative to m_inv)
      //   qt2=min_qTcut*min_qTcut*exp(log(maxqt2/min_qTcut/min_qTcut)*random_qt2);
      //   jacqt2=log(maxqt2/min_qTcut/min_qTcut)*qt2;
      logger << LOG_FATAL << "old qT variation is deprecated" << endl;
      exit(0);
    }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
/*
// original version: check if everything works without this one...
void phasespace_set::calculate_IS_CT(){
  Logger logger("phasespace_set::calculate_IS_CT");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  //  cout << "QT_g_IS_ = " << setw(23) << setprecision(15) << QT_g_IS_ << endl;
  //  cout << "QT_eps = " << setw(23) << setprecision(15) << QT_eps << endl;
  
  //  double QT_random_qt2;
  
  QT_random_z1->get_random(QT_random_IS,QT_g_IS_);
  zz_pdf[1]=pow(x_pdf[1],QT_eps+(1-QT_eps)*QT_random_IS);
  QT_g_z1=log(x_pdf[1])*(1-QT_eps);
  //  cout << "z1: QT_random_IS = " << setw(23) << setprecision(15) << QT_random_IS << "   zz_pdf[1] = " << setw(23) << setprecision(15) << zz_pdf[1] << "   QT_g_z1 = " << setw(23) << setprecision(15) << QT_g_z1 << "   QT_g_IS_ = " << setw(23) << setprecision(15) << QT_g_IS_ << endl;
  g_pdf*=QT_g_IS_;
  
  QT_random_z2->get_random(QT_random_IS,QT_g_IS_);
  zz_pdf[2]=pow(x_pdf[2],QT_eps+(1-QT_eps)*QT_random_IS);
  QT_g_z2=log(x_pdf[2])*(1-QT_eps);
  //  cout << "z2: QT_random_IS = " << setw(23) << setprecision(15) << QT_random_IS << "   zz_pdf[2] = " << setw(23) << setprecision(15) << zz_pdf[2] << "   QT_g_z2 = " << setw(23) << setprecision(15) << QT_g_z2 << "   QT_g_IS_ = " << setw(23) << setprecision(15) << QT_g_IS_ << endl;
  g_pdf*=QT_g_IS_;
  
  QT_random_qt->get_random(QT_random_qt2,QT_g_IS_);
  //  cout << "qt: QT_random_qt2 = " << setw(23) << setprecision(15) << QT_random_qt2 << "   QT_g_IS_ = " << setw(23) << setprecision(15) << QT_g_IS_ << endl;
  g_pdf*=QT_g_IS_;
  
  // for pdf calculation
  z_pdf[1]=x_pdf[1]/zz_pdf[1];
  z_pdf[2]=x_pdf[2]/zz_pdf[2];
  z_pdf[0]=z_pdf[1]*z_pdf[2];
  
  // the (Bessel transformed) large logarithms vanish up to double precision accuracy for arguments larger than 40 -> use this as upper bound
  // 1.261 ~ b0*b0
  double maxqt2;
  double Qres;
  
  if (!do_resummation) { // f.o.; need to use correct lower limit for qT integral
    if (Qres_cut==0 && !dynamical_Qres) { // not specified -> assume standard f.o. computation
      Qres=sqrt(xbs_all[0][0]);
    } else if (dynamical_Qres) {
      Qres=sqrt(xbs_all[0][0])*Qres_prefactor;
    } else { // use specified fixed Qres
      Qres=Qres_cut;
    }
  }
  else { // resummation; there is (usually) no point in using Q_res as a reference scale for the cut. in fact, it only slows down convergence
    Qres=sqrt(xbs_all[0][0]);
  }

  if (contribution_order_alpha_s[0] == 1){maxqt2 = 1.261 * 30. * 20. * Qres * Qres;}
  else if (contribution_order_alpha_s[0] == 2){maxqt2 = 1.261 * 40. * 50. * Qres * Qres;}
  //  if (contribution_order_alpha_s[0] == 1){maxqt2 = 1.261 * 30. * 20. * xbs_all[0][0];}
  //  else if (contribution_order_alpha_s[0] == 2){maxqt2 = 1.261 * 40. * 50. * xbs_all[0][0];} // TODO: is this really necessary?
  //  cout << "contribution_order_alpha_s[0] = " << contribution_order_alpha_s[0] << endl;
  //  cout << "min_qTcut = " << min_qTcut << endl;

  if (switch_qTcut == 1){
    double r0 = min_qTcut / 100.;
    //  version with cut on pT/Qres_cut
    QT_qt2 = r0 * r0 * xbs_all[0][0] * exp(log(maxqt2 / r0 / r0 / Qres / Qres) * QT_random_qt2);
    QT_jacqt2 = log(maxqt2 / r0 / r0 / Qres / Qres) * QT_qt2;
    //  version with cut on pT/m_inv
    //    QT_qt2 = r0 * r0 * xbs_all[0][0] * exp(log(maxqt2 / r0 / r0 / xbs_all[0][0]) * QT_random_qt2);
    //    QT_jacqt2 = log(maxqt2 / r0 / r0 / xbs_all[0][0]) * QT_qt2;
  }
  else if (switch_qTcut == 2){
   // old version with cut on pT only (not relative to m_inv)
   //   qt2=min_qTcut*min_qTcut*exp(log(maxqt2/min_qTcut/min_qTcut)*random_qt2);
   //   jacqt2=log(maxqt2/min_qTcut/min_qTcut)*qt2;
    logger << LOG_FATAL << "old qT variation is deprecated" << endl;
    exit(0);
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/


void phasespace_set::calculate_IS_CX_from_QT(){
  Logger logger("phasespace_set::calculate_IS_CX");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
// for comparison with old result only !!!    

  for (int i_z = 1; i_z < 3; i_z++){
    z_coll[i_z] = zz_pdf[i_z];
  }
  // for pdf calculation
  z_pdf[1]=x_pdf[1] / zz_pdf[1];
  z_pdf[2]=x_pdf[2] / zz_pdf[2];
  z_pdf[0]=z_pdf[1] * z_pdf[2];
 
  g_z_coll[1] = -1. / QT_g_z1;
  g_z_coll[2] = -1. / QT_g_z2;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::calculate_IS_QT_from_CX(){
  Logger logger("phasespace_set::calculate_IS_CX");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
// for comparison with old result only !!!    

  for (int i_z = 1; i_z < 3; i_z++){
    zz_pdf[i_z] = z_coll[i_z];
  }
  // for pdf calculation
  z_pdf[1]=x_pdf[1] / zz_pdf[1];
  z_pdf[2]=x_pdf[2] / zz_pdf[2];
  z_pdf[0]=z_pdf[1] * z_pdf[2];
 
  QT_g_z1 = -1. / g_z_coll[1];
  QT_g_z2 = -1. / g_z_coll[2];

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::calculate_IS_CA_from_CT(){
  Logger logger("phasespace_set::calculate_IS_CA_from_CT");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_z = 1; i_z < 3; i_z++){
    z_coll[i_z] = zz_pdf[i_z];
  }

  g_z_coll[1] = -1. / QT_g_z1;
  g_z_coll[2] = -1. / QT_g_z2;

  for (int i_z = 0; i_z < 3; i_z++){
    all_xz_coll_pdf[i_z][0] = x_pdf[i_z];
  }
  for (int i_c = 1; i_c < 3; i_c++){
    all_xz_coll_pdf[0][i_c] = z_coll[i_c];
    all_xz_coll_pdf[i_c][i_c] = x_pdf[i_c] / z_coll[i_c];
    all_xz_coll_pdf[i_c % 2 + 1][i_c] = x_pdf[0] / z_coll[i_c];
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::calculate_IS_VT(){
  Logger logger("phasespace_set::calculate_IS_VT");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  QT_random_z1->get_random(QT_random_IS,QT_g_IS_);
  zz_pdf[1]=pow(x_pdf[1],QT_eps+(1-QT_eps)*QT_random_IS);
  QT_g_z1=log(x_pdf[1])*(1-QT_eps);
  g_pdf*=QT_g_IS_;
  
  QT_random_z2->get_random(QT_random_IS,QT_g_IS_);
  zz_pdf[2]=pow(x_pdf[2],QT_eps+(1-QT_eps)*QT_random_IS);
  QT_g_z2=log(x_pdf[2])*(1-QT_eps);
  g_pdf*=QT_g_IS_;

  // in order to (hopefully) circumvent seg fault...
  z_coll.resize(3);
  g_z_coll.resize(3);
  
  // for pdf calculation
  z_pdf[1]=x_pdf[1]/zz_pdf[1];
  z_pdf[2]=x_pdf[2]/zz_pdf[2];
  z_pdf[0]=z_pdf[1]*z_pdf[2];
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

/*
void phasespace_set::calculate_IS_VT(){
  Logger logger("phasespace_set::calculate_IS_VT");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  QT_random_z1->get_random(QT_random_IS,QT_g_IS_);
  zz_pdf[1]=pow(x_pdf[1],QT_eps+(1-QT_eps)*QT_random_IS);
  QT_g_z1=log(x_pdf[1])*(1-QT_eps);
  g_pdf*=QT_g_IS_;
  
  QT_random_z2->get_random(QT_random_IS,QT_g_IS_);
  zz_pdf[2]=pow(x_pdf[2],QT_eps+(1-QT_eps)*QT_random_IS);
  QT_g_z2=log(x_pdf[2])*(1-QT_eps);
  g_pdf*=QT_g_IS_;
  /*
  QT_random_qt->get_random(QT_random_qt2,QT_g_IS_);
  g_pdf*=QT_g_IS_;
*//*
  // !!! needs to be transported into new version !!!
  if (do_resummation) {
    min_qTcut=0.1;
    double QT_random_qt2;
    
    QT_random_qt->get_random(QT_random_qt2,QT_g_IS_);
    g_pdf*=QT_g_IS_;

    double maxqt2=pow(2*E,2);
    maxqt2=pow(E,2);

    if (lambda_qt2==0 || lambda_qt2==1) {
      QT_qt2=min_qTcut*min_qTcut*exp(log(maxqt2/min_qTcut/min_qTcut)*QT_random_qt2);
      QT_jacqt2=log(maxqt2/min_qTcut/min_qTcut)*QT_qt2;
    }

  //  random_qt2=0;

  //  double lambda=0.0355064;
  //  double C=-1.0/lambda*(exp(-lambda*sqrt(maxqt2))-exp(-lambda*min_qTcut));
  //  double y0=-1.0/lambda/C*exp(-lambda*min_qTcut);
  ////  cout << exp(lambda*sqrt(maxqt2)) << ", " << exp(lambda*min_qTcut) << endl;
  ////  cout << C << ", " << y0 << endl;
  //  double qt=-1.0/lambda*log(-lambda*C*(random_qt2+y0));
  //  double jacqt=C*exp(lambda*qt);
  //  qt2=qt*qt;
  //  jacqt2=jacqt*2*qt;
    else {
      double y0=1.0/(pow(maxqt2/min_qTcut/min_qTcut,1-lambda_qt2)-1);
      double C=pow(min_qTcut,2*(1-lambda_qt2))/(1-lambda_qt2)/y0;

      QT_qt2=pow((1-lambda_qt2)*C*(QT_random_qt2+y0),1.0/(1-lambda_qt2));
      QT_jacqt2=C*pow(QT_qt2,lambda_qt2);

  //    cout << C << ", " << y0 << ", " << (1-lambda_qt2)*C*(random_qt2+y0) << endl;
  //    cout << qt2 << ", " << maxqt2 << endl;
    }

  //  maxqt2 = 1e6;
  //  qt2 = min_qTcut*min_qTcut+random_qt2*(maxqt2-min_qTcut*min_qTcut);
  //  jacqt2 = (maxqt2-min_qTcut*min_qTcut);

  }

  // in order to (hopefully) circumvent seg fault...
  z_coll.resize(3);
  g_z_coll.resize(3);
  
  // for pdf calculation
  z_pdf[1]=x_pdf[1]/zz_pdf[1];
  z_pdf[2]=x_pdf[2]/zz_pdf[2];
  z_pdf[0]=z_pdf[1]*z_pdf[2];
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/

