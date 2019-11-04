#include "classes.cxx"
#include <header/pptt40.amplitude.doublevirtual.h>
#include "definitions.observable.set.cxx"

extern "C" {
  void init_interpolation_t3t4_();
  void init_interpolation_tttt_();
  double h1_tt_(int *i1, int *i2, double* beta, double* cost);
  double h2_tt_(int *i1, int *i2, double* beta, double* cost);
  double h2_tttt_(int *i1, int *i2, int *i3, int *i4, double* beta, double* cost);
  double virtgg_(double* beta, double* cost);
  double virtqq_(double* beta, double* cost);
  void initialization_momenta_(double* P, double* mu_Q, double* m_t);
  void fourcorrelators_(double* B4, int *i1, int *i2, int *i3, int *i4, int *channel);
  double_complex hpl2_(int *i1, int *i2, double_complex *expB34);
  double_complex hpl3_(int *i1, int *i2, int *i3, double_complex *expB34);
}


void pptt40_calculate_H2(observable_set & oset){
  static Logger logger("pptt40_calculate_H2");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static double CA = C_A;
  static double nf = double(oset.N_f);
  // for the moment only implemented for nf=5
  if (oset.N_f != 5) {logger << LOG_FATAL << "H2 only implemented for nf = 5." << endl; exit(1);}
  static double b0 = 11. / 2 - nf / 3;
  static double b1 = 51. / 2 - 19. / 6 * nf;

  static int temp_channel = 0;
  static double Ci = 0.;
  //  double virtxx = 0.;
  static double gama1 = 0.;
  static double gama2 = 0.;
  if (initialization){
    init_interpolation_t3t4_();
    init_interpolation_tttt_();

    if (oset.initial_gg){
      temp_channel = 1;
      Ci = C_A;
      //    virtxx = virtgg_(&beta, &cost);
      //    logger << LOG_DEBUG_POINT << "virtgg = " << setw(23) << setprecision(17) << showpoint << virtxx << endl;
      gama1 = -b0 / 2;
      gama2 = (2 * C_F * nf + (C_A * (256. / 27 - (2 * pi2) / 9) * nf) / 2
	       + C_A*C_A * (-692. / 27 + (11. * pi2) / 18 + 2 * zeta3)) / 16;
    }
    else if (oset.initial_qqx){
      temp_channel = 2;
      Ci = C_F;
      //      virtxx = virtqq_(&beta, &cost);
      //      logger << LOG_DEBUG_POINT << "virtqq = " << setw(23) << setprecision(17) << showpoint << virtxx << endl;
      gama1 = (-3 * C_F) / 4;
      gama2 = ((C_F * (130. / 27 + (2 * pi2) / 3) * nf ) / 2
	       + C_F*C_F * (-3. / 2 + 2 * pi2 - 24 * zeta3)
	       + C_A * C_F * (-961. / 54 - (11 * pi2) / 6 + 26 * zeta3)) / 16;
    }
    initialization = 0;
  }

  //  double ME2_L_L = 0.;
  //  double ME2_L2_B = 0.;

  /*
  double *M2L2;
  M2L2 = new double[5];
  for (int i_e = 0; i_e < 5; i_e++){M2L2[i_e] = 0.;}
  */
  double M1M1 = 0.;
  double M0M0 = 0.;
  double M0M0_3434 = 0.;
  double M2M0 = 0.;
  vector<vector<vector<vector<double> > > > M0M04c(5, vector<vector<vector<double> > > (5, vector<vector<double> > (5, vector<double> (5, 0.))));
  vector<vector<double> > M0M0cc(5, vector<double> (5, 0.));
  vector<vector<double> > M1M0cc(5, vector<double> (5, 0.));
  vector<vector<double> > imM1M0cc(5, vector<double> (5, 0.));
  double M1M0 = 0.;
  double imM1M0 = 0.;

  double m_HQ = oset.msi.M[oset.csi->type_parton[0][3]];
  double m2_HQ = m_HQ*m_HQ;


  logger << LOG_DEBUG_VERBOSE << "oset.switch_OL = " << oset.switch_OL << endl;
  if (oset.switch_OL){

  static double one = 1.;
  static int n_momentum = 5 * (osi_n_particle + 2);
  double *P;
  P = new double[n_momentum];
  for (int i = 1; i < osi_p_parton[0].size(); i++){
    P[5 * (i - 1)]     = osi_p_parton[0][i].x0();
    P[5 * (i - 1) + 1] = osi_p_parton[0][i].x1();
    P[5 * (i - 1) + 2] = osi_p_parton[0][i].x2();
    P[5 * (i - 1) + 3] = osi_p_parton[0][i].x3();
    P[5 * (i - 1) + 4] = osi_p_parton[0][i].m();
  }

  int CT_on = 1; // modification of ...ME2_VA and OpenLoops calls needed !!!
  ol_setparameter_int(stch("ct_on"), CT_on); // modification of ...ME2_VA and OpenLoops calls needed !!!

  static double acc;

  static char * OL_mu_ren = stch("muren");
  static char * OL_mu_reg = stch("mureg");
  static char * pole_uv = stch("pole_uv");
  static char * pole_ir1 = stch("pole_ir1");
  static char * pole_ir2 = stch("pole_ir2");
  ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
  ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
  ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);

  double mu_Q = (osi_p_parton[0][1] + osi_p_parton[0][2]).m();
  ol_setparameter_double(OL_mu_ren, mu_Q);
  ol_setparameter_double(OL_mu_reg, mu_Q);
  static char * fact_uv = stch("fact_uv");
  static char * fact_ir = stch("fact_ir");
  ol_setparameter_double(fact_uv, one);
  ol_setparameter_double(fact_ir, one);

  initialization_momenta_(P, &mu_Q, &m_HQ);

  double dummy = 0.;
  double m0_OL = 0.;
  double m2cc[16];
  double *m2l1;
  m2l1 = new double[3];
  double m2loopcc[16];
  double *m2l1_im;
  m2l1_im = new double[3];
  double m2loopcc_im[16];
  double *m2l2;
  m2l2 = new double[5];
  /*
  double *m2cc;
  m2cc = new double[5][5];
  */

  // Needed to initialize momenta !!! Should be separated from loopxBorn evaluation !!!
  /*
  double m_t = osi_msi.M_t;
  double temp_H1 = 0.;
  double temp_H1_T34 = 0.;
  double temp_H1_T13 = 0.;
  double temp_H1_T23 = 0.;
  deltaext_(P, &mu_Q, &m_t, &temp_H1, &temp_H1_T34, &temp_H1_T13, &temp_H1_T23, &temp_channel);
  temp_H1 = temp_H1 * pow(osi_alpha_S, 2);
  temp_H1_T34 = temp_H1_T34 * pow(osi_alpha_S, 2);
  temp_H1_T13 = temp_H1_T13 * pow(osi_alpha_S, 2);
  //  temp_H1_T23 = temp_H1_T23 * pow(osi_alpha_S, 2);
  // temp_H1_T23 is not calculated directly -> colour conservation !!!
  temp_H1_T23 = -(temp_H1_T34 + temp_H1_T13 + C_F * temp_H1);
  // alpha_S running is included here (subtract via cc Born), which is again done later in observable.qTsubtraction.cpp routines.
  logger << LOG_DEBUG_POINT << "Hayk:   osi_QT_H1_delta_Hayk     = " << temp_H1 * 2 << endl;
  logger << LOG_DEBUG_POINT << "Hayk:   osi_QT_H1_T34_delta_Hayk = " << temp_H1_T34 * 2 << endl;
  logger << LOG_DEBUG_POINT << "Hayk:   osi_QT_H1_T13_delta_Hayk = " << temp_H1_T13 * 2 << endl;
  logger << LOG_DEBUG_POINT << "Hayk:   osi_QT_H1_T23_delta_Hayk = " << temp_H1_T23 * 2 << endl;
  logger << LOG_DEBUG_POINT << "Hayk:   T13_T23_T34_T33_Hayk     = " << temp_H1_T34 + temp_H1_T13 + temp_H1_T23 + C_F * temp_H1 << endl;
  */

  ol_evaluate_ccmatrix(1, P, &m0_OL, m2cc, &acc);
  osi_VA_b_ME2 = m0_OL;

  /////  ol_evaluate_loop(1, P, dummy, m2l1, acc);
  ol_evaluate_loopccmatrix(1, P, &dummy, m2l1, m2loopcc, &acc);
  osi_VA_V_ME2 = m2l1[0];
  /*
  logger << LOG_DEBUG_POINT << "m2l1[0]    = " << m2l1[0] << endl;
  for (int i = 0; i < 16; i++){
    logger << LOG_DEBUG_POINT << "m2loopcc[" << i % 4 + 1 << "][" << i / 4 + 1 << "] = " << m2loopcc[i] << endl;
  }
  */

  /////  ol_evaluate_loop(3, P, dummy, m2l1_im, acc);
  ///// !!! 1 -> 3 (imaginary part)
  ol_evaluate_loopccmatrix(3, P, &dummy, m2l1_im, m2loopcc_im, &acc);
  //  ol_evaluate_loopccmatrix(1, P, &dummy, m2l1_im, m2loopcc_im, &acc);

  logger << LOG_DEBUG_POINT << "m2l1_im[0]    = " << m2l1_im[0] << endl;
  for (int i = 0; i < 16; i++){
    logger << LOG_DEBUG_POINT << "m2loopcc_im[" << i % 4 + 1 << "][" << i / 4 + 1 << "] = " << m2loopcc_im[i] << endl;
  }


  ol_evaluate_loop2(2, P, m2l2, &acc);



  static int n_cc = (osi_n_particle + 2) * (osi_n_particle + 1) / 2;
  static double ewcc = 0.;

  double *xxxm2l1_im;
  xxxm2l1_im = new double[3];
  static double *xxxM2loopcc_im;
  xxxM2loopcc_im = new double[n_cc];
  ol_evaluate_loopcc(3, P, &dummy, xxxm2l1_im, xxxM2loopcc_im, &ewcc);

  logger << LOG_DEBUG_POINT << "xxxm2l1_im[0] = " << xxxm2l1_im[0] << endl;
  for (int i = 0; i < n_cc; i++){
    logger << LOG_DEBUG_POINT << "xxxM2loopcc_im[" << i << "] = " << xxxM2loopcc_im[i] << endl;
  }

  for (int i = 1; i < 5; i++){
    for (int j = 1; j < i + 1; j++){
      if (i == j && i < 3){
	imM1M0cc[i][j] = Ci * xxxm2l1_im[0];
      }
      else if (i == j && i > 2){
	imM1M0cc[i][j] = C_F * xxxm2l1_im[0];
      }
      else{
	imM1M0cc[i][j] = xxxM2loopcc_im[(i-1)*(i-2)/2+(j-1)];
	imM1M0cc[j][i] = xxxM2loopcc_im[(i-1)*(i-2)/2+(j-1)];
      }
    }
  }


  double *xxxm2l1;
  xxxm2l1 = new double[3];
  static double *xxxM2loopcc;
  xxxM2loopcc = new double[n_cc];
  ol_evaluate_loopcc(1, P, &dummy, xxxm2l1, xxxM2loopcc, &ewcc);
  /*
  logger << LOG_DEBUG_POINT << "xxxm2l1[0] = " << xxxm2l1[0] << endl;
  for (int i = 0; i < n_cc; i++){
    logger << LOG_DEBUG_POINT << "xxxM2loopcc[" << i << "] = " << xxxM2loopcc[i] << endl;
   }
  */
  for (int i = 1; i < 5; i++){
    for (int j = 1; j < i + 1; j++){
      if (i == j && i < 3){
	M1M0cc[i][j] = Ci * xxxm2l1[0];
      }
      else if (i == j && i > 2){
	M1M0cc[i][j] = C_F * xxxm2l1[0];
      }
      else{
	M1M0cc[i][j] = xxxM2loopcc[(i-1)*(i-2)/2+(j-1)];
	M1M0cc[j][i] = xxxM2loopcc[(i-1)*(i-2)/2+(j-1)];
      }
    }
  }



  /*
  12
  13
  23
  14
  24
  34
  */

  M0M0 = m0_OL;
  M1M1 = m2l2[0];
  logger << LOG_DEBUG_POINT << "M1M1 = " << M1M1 << endl;

  /*
  for (int i = 0; i < 16; i++){
    logger << LOG_DEBUG_POINT << "M0M0cc[" << i % 4 + 1 << "][" << i / 4 + 1 << "] = m2cc[" << setw(2) << i << "] = " << m2cc[i] << endl;
  }
  */
  for (int i = 0; i < 16; i++){
    M0M0cc[i % 4 + 1][i / 4 + 1] = m2cc[i];
    // in order to switch back to the loopccmatrix-interface
    /////    M1M0cc[i % 4 + 1][i / 4 + 1] = m2loopcc[i];
    /////    imM1M0cc[i % 4 + 1][i / 4 + 1] = m2loopcc_im[i];
    // only for tests:
    //    imM1M0cc[i % 4 + 1][i / 4 + 1] = M0M0cc[i % 4 + 1][i / 4 + 1] ;
    //    imM1M0cc[i % 4 + 1][i / 4 + 1] = 0.;
  }




  //C The M1*M0 + cc I can take from OpenLoops, also the cc ones:
  M1M0 = m2l1[0] / (osi_alpha_S * inv2pi);
  imM1M0 = xxxm2l1_im[0] / (osi_alpha_S * inv2pi);
  // only for tests:
  //  imM1M0 = M1M0;
  //  imM1M0 = 0.;
  logger << LOG_DEBUG_POINT << "norm   imM1M0 = " << imM1M0 << endl;
  for (int i = 1; i < 5; i++){
    for (int j = 1; j < 5; j++){
      M1M0cc[i][j] /= (osi_alpha_S * inv2pi);
      imM1M0cc[i][j] /= (osi_alpha_S * inv2pi);
    }
  }
  //C The M1*M1 I also take from OL
  M1M1 = M1M1 / pow(osi_alpha_S * inv2pi, 2);
  logger << LOG_DEBUG_POINT << "norm   M1M1 = " << M1M1 << endl;
  for (int i = 0; i < 16; i++){
    logger << LOG_DEBUG_POINT << "norm   M1M0cc[" << i % 4 + 1 << "][" << i / 4 + 1 << "] = " << M1M0cc[i % 4 + 1][i / 4 + 1] << endl;
  }






  //  osi_VA_b_ME2 = 0.;

  /*
  double var_mu_ren = (osi_p_parton[0][1] + osi_p_parton[0][2]).m();
  static char * renscale = stch("renscale");
  ol_setparameter_double(renscale, var_mu_ren);
  */
  /*
  logger.newLine(LOG_DEBUG_POINT);
    logger << LOG_DEBUG_VERBOSE << "Before   ol_evaluate_loop2(1, P, M2L2, &acc);" << endl;
    ol_evaluate_loop2(1, P, M2L2, &acc);
    logger << LOG_DEBUG_VERBOSE << "After   ol_evaluate_loop2(1, P, M2L2, &acc);" << endl;
    logger << LOG_DEBUG_POINT << "ol_evaluate_loop" << endl;

    osi_VA_b_ME2 = M2L2[0];

    for (int i_e = 0; i_e < 5; i_e++){
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_L2I (eps^" << i_e << ") = " << setw(23) << setprecision(15) << M2L2[i_e] << endl;
    }
    logger.newLine(LOG_DEBUG_POINT);
    delete [] M2L2;


    double *M2L1;
    double *IRL1;
    M2L1 = new double[3];
    IRL1 = new double[3];
    //    double *M2L2;
    double *IRL2;
    M2L2 = new double[5];
    IRL2 = new double[5];
    logger << LOG_DEBUG_POINT << "ol_evaluate_full" << endl;
    ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_B           = " << setw(23) << setprecision(15) << osi_VA_b_ME2 << endl;
      logger.newLine(LOG_DEBUG_POINT);
    for (int i_e = 0; i_e < 3; i_e++){
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_B*V (eps^" << i_e << ") = " << setw(23) << setprecision(15) << M2L1[i_e] << endl;
    }
      logger.newLine(LOG_DEBUG_POINT);
    for (int i_e = 0; i_e < 5; i_e++){
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_L2I (eps^" << i_e << ") = " << setw(23) << setprecision(15) << M2L2[i_e] << endl;
    }
      logger.newLine(LOG_DEBUG_POINT);

    delete [] M2L1;
    delete [] IRL1;
    delete [] M2L2;
    delete [] IRL2;
*/
    CT_on = 0; // modification of ...ME2_VA and OpenLoops calls needed !!!
    ol_setparameter_int(stch("ct_on"), CT_on); // modification of ...ME2_VA and OpenLoops calls needed !!!
  }


  double sn = (osi_p_parton[0][1] + osi_p_parton[0][2]).m2();
  //  double sn = 2 * osi_p_parton[0][1] * osi_p_parton[0][2];
  double tn = (osi_p_parton[0][1] - osi_p_parton[0][3]).m2();
  //  2d0*dot(p,1,3)+mt**2.d0;
  double un = (osi_p_parton[0][1] - osi_p_parton[0][4]).m2();
    //2d0*dot(p,2,3)+mt**2.d0;
  double beta = sqrt(1. - 4 * m2_HQ / sn);
  //  double beta = dsqrt(1.d0-4.d0*mt**2.d0/sn);
  //  double cost = (tn-un)/sn/beta;
  double cost = (tn - un) / sn / beta;
  logger << LOG_DEBUG_POINT << "beta = " << setw(23) << setprecision(17) << showpoint << beta << endl;
  logger << LOG_DEBUG_POINT << "cost = " << setw(23) << setprecision(17) << showpoint << cost << endl;

  /*
  //  only once per run:
  double CA = C_A;
  double nf = double(oset.N_f);
  init_interpolation_t3t4_();
  init_interpolation_tttt_();
  */
  int i1 = 1;
  int i2 = 2;
  int i3 = 3;
  int i4 = 4;
  //  double test = 0.;
  /*
  beta = .5;
  cost = .5;

  beta = sqrt(1.-0.5);
  cost = 0.974679434480896;

  //  test = h1_TT_(&i1, &i2, &beta, &costheta);
  //  test = h1_TT_(&i1, &i2, &beta, &costheta);
  //  test = h1_TTTT_(&i1, &i2, &i3, &i4, &beta, &costheta);
  //  test = virtgg_(&beta, &cost);
  //  cout << "test = " << test << endl;
  //  test = xxxvirtgg_(&beta, &cost);
  //  cout << "h1_tt_(&i3,4) = " << h1_tt_(&beta,&cost) << endl;
  //  double x1 = h1_tt_(&i3,&i4,&beta,&cost);
  cout << "h1_tt_(&i3,4) = " << setw(23) << setprecision(17) << showpoint << h1_tt_(&i3,&i4,&beta,&cost) << endl;
  cout << "h1_tt_(&i3,3) = " << setw(23) << setprecision(17) << showpoint << h1_tt_(&i3,&i3,&beta,&cost) << endl;
  cout << "h2_tt_(&i3,4) = " << setw(23) << setprecision(17) << showpoint << h2_tt_(&i3,&i4,&beta,&cost) << endl;
  cout << "h2_tt_(&i3,1) = " << setw(23) << setprecision(17) << showpoint << h2_tt_(&i3,&i1,&beta,&cost) << endl;
  cout << "h2_tt_(&i4,2) = " << setw(23) << setprecision(17) << showpoint << h2_tt_(&i4,&i2,&beta,&cost) << endl;
  cout << "h2_tt_(&i4,1) = " << setw(23) << setprecision(17) << showpoint << h2_tt_(&i4,&i1,&beta,&cost) << endl;
  cout << "h2_tt_(&i3,2) = " << setw(23) << setprecision(17) << showpoint << h2_tt_(&i3,&i2,&beta,&cost) << endl;
  cout << "h2_tttt_(&i3,&i4,&i3,4) = " << setw(23) << setprecision(17) << showpoint << h2_tttt_(&i3,&i4,&i3,&i4,&beta,&cost) << endl;
  cout << "h2_tttt_(&i3,&i4,&i3,3) = " << setw(23) << setprecision(17) << showpoint << h2_tttt_(&i3,&i4,&i3,&i3,&beta,&cost) << endl;
  cout << "h2_tttt_(&i3,&i3,&i3,3) = " << setw(23) << setprecision(17) << showpoint << h2_tttt_(&i3,&i3,&i3,&i3,&beta,&cost) << endl;
  for (int i = 2; i < 21; i++){
    beta = 1. / i;
    cout << setw(23) << setprecision(17) << showpoint << beta << "   " << setw(23) << setprecision(17) << showpoint << virtgg_(&beta, &cost) << "   " << setw(23) << setprecision(17) << showpoint << virtqq_(&beta, &cost) << endl;
  }
  */
  //C I evaluate the M2*M0 + cc. Following equation 3.6 of
  //C 1312.6279, this can be obtained from the grid as

  double virtxx = 0.;
  if (oset.initial_gg){
    virtxx = virtgg_(&beta, &cost);
    logger << LOG_DEBUG_POINT << "virtgg = " << setw(23) << setprecision(17) << showpoint << virtxx << endl;
  }
  else if (oset.initial_qqx){
    virtxx = virtqq_(&beta, &cost);
    logger << LOG_DEBUG_POINT << "virtqq = " << setw(23) << setprecision(17) << showpoint << virtxx << endl;
  }

  M2M0 = virtxx * 16 * pi / (beta * (1. - pow(beta, 2)));
  /////    M2M0 = virtgg_(&beta, &cost) * 16 * pi / (beta * (1. - pow(beta, 2)));
  /////    M2M0 = virtgg_(&beta, &cost) * 4096 * pi / (beta * (1. - pow(beta, 2)));
  //C where the expansion for M is done in powers of as/(2*pi).
  //C To be consistent with the rest, I add the overall power
  //C of as, which otherwise is missing:
  M2M0 = M2M0 * pow(4 * pi * osi_alpha_S, 2);
  //    logger << LOG_DEBUG_POINT << "4 * pi * osi_alpha_S = " << setw(23) << setprecision(17) << showpoint << 4 * pi * osi_alpha_S << endl;

  //C I add factors from avergage of initial state d.o.f
  /////    double facgg = 1. / (8 * 8 * 2 * 2);
  /////    M2M0 = M2M0 * facgg;

  double temp_B4;
  fourcorrelators_(&temp_B4, &i3, &i4, &i3, &i4, &temp_channel);
  // C Finally, the tree-level 4 parton correlator 3434 I take from Hayk
  M0M0_3434 = temp_B4 * M0M0;
  logger << LOG_DEBUG_POINT << "fourcorrelators_(&temp_B4, " << i3 << ", " << i4 << ", " << i3 << ", " << i4 << ", " << temp_channel << "); -> temp_B4 = " << temp_B4 << endl;

  /////    M0M0 = m0_OL;
  //C The convention from Massimiliano is different from the one in
  //C OpenLoops, so I need to fix that:

  double fac = Ci * pi2_12;
  for (int i = 1; i < 5; i++){
    for (int j = 1; j < 5; j++){
      M1M0cc[i][j] = M1M0cc[i][j] + 2 * fac * M0M0cc[i][j];
    }
  }
  M1M1 = M1M1 + fac * M1M0 + pow(fac, 2) * M0M0;
  M1M0 = M1M0 + 2 * fac * M0M0;

  logger << LOG_DEBUG_POINT << "norm   M0M0 = " << setw(23) << setprecision(17) << showpoint << M0M0 << endl;
  logger << LOG_DEBUG_POINT << "norm   M0M0_3434 = " << setw(23) << setprecision(17) << showpoint << M0M0_3434 << endl;
  logger << LOG_DEBUG_POINT << "norm   imM1M0 = " << setw(23) << setprecision(17) << showpoint << imM1M0 << endl;
  logger << LOG_DEBUG_POINT << "norm   M2M0 = " << setw(23) << setprecision(17) << showpoint << M2M0 << endl;
  logger << LOG_DEBUG_POINT << "shifted   M1M0 = " << setw(23) << setprecision(17) << showpoint << M1M0 << endl;
  logger << LOG_DEBUG_POINT << "shifted   M1M1 = " << setw(23) << setprecision(17) << showpoint << M1M1 << endl;
  for (int i = 0; i < 16; i++){logger << LOG_DEBUG_POINT << "M0M0cc[" << i % 4 + 1 << "][" << i / 4 + 1 << "] = " << M0M0cc[i % 4 + 1][i / 4 + 1] << endl;}
  for (int i = 0; i < 16; i++){logger << LOG_DEBUG_POINT << "shifted   M1M0cc[" << i % 4 + 1 << "][" << i / 4 + 1 << "] = " << M1M0cc[i % 4 + 1][i / 4 + 1] << endl;}
  for (int i = 0; i < 16; i++){logger << LOG_DEBUG_POINT << "norm   imM1M0cc[" << i % 4 + 1 << "][" << i / 4 + 1 << "] = " << imM1M0cc[i % 4 + 1][i / 4 + 1] << endl;}

  osi_QT_H2_delta =
    M1M1
    + M2M0
    + (M1M0 * (-(Ci * pi2) / 3 + 4 * C_F * h1_tt_(&i3, &i3, &beta, &cost))) / 2
    + M0M0 * ((Ci*Ci * pi4) / 72
	      + (Ci * (C_A * (4856. - 603. * pi2 + 18. * pi4 - 2772. * zeta3)
		       + nf * (-656. + 90. * pi2 + 504. * zeta3)
		       - 216. * C_F * pi2 * h1_tt_(&i3, &i3, &beta, &cost))) / 648
	      + C_F*C_F * (pow(h1_tt_(&i3, &i3, &beta, &cost), 2)
			   + 2 * h2_tttt_(&i3, &i3, &i3, &i3, &beta, &cost)))
    + (pow(h1_tt_(&i3,&i4, &beta, &cost), 2)
       + 2 * h2_tttt_(&i3, &i4, &i3, &i4, &beta, &cost)) * M0M0_3434
    + 2 * h2_tt_(&i3, &i1, &beta, &cost) * M0M0cc[1][3]
    + 2 * h2_tt_(&i4, &i1, &beta, &cost) * M0M0cc[1][4]
    + 2 * h2_tt_(&i3, &i2, &beta, &cost) * M0M0cc[2][3]
    + 2 * h2_tt_(&i4, &i2, &beta, &cost) * M0M0cc[2][4]
    + (-(Ci * pi2 * h1_tt_(&i3, &i4, &beta, &cost)) / 3
       + 2 * (h2_tt_(&i3, &i4, &beta, &cost)
	      + C_F * (h1_tt_(&i3, &i3, &beta, &cost) * h1_tt_(&i3, &i4, &beta, &cost)
		       + h2_tttt_(&i3, &i4, &i3, &i3, &beta, &cost)))) * M0M0cc[3][4]
    + 2 * h1_tt_(&i3, &i4, &beta, &cost) * M1M0cc[3][4];

  logger << LOG_DEBUG_POINT << "h1_tt_(" << i3 << ", " << i3 << ") = " << setw(23) << setprecision(17) << showpoint << h1_tt_(&i3,&i3,&beta,&cost) << endl;
  logger << LOG_DEBUG_POINT << "h1_tt_(" << i3 << ", " << i4 << ") = " << setw(23) << setprecision(17) << showpoint << h1_tt_(&i3,&i4,&beta,&cost) << endl;
  logger << LOG_DEBUG_POINT << "h2_tt_(" << i3 << ", " << i1 << ") = " << setw(23) << setprecision(17) << showpoint << h2_tt_(&i3,&i1,&beta,&cost) << endl;
  logger << LOG_DEBUG_POINT << "h2_tt_(" << i4 << ", " << i1 << ") = " << setw(23) << setprecision(17) << showpoint << h2_tt_(&i4,&i1,&beta,&cost) << endl;
  logger << LOG_DEBUG_POINT << "h2_tt_(" << i3 << ", " << i2 << ") = " << setw(23) << setprecision(17) << showpoint << h2_tt_(&i3,&i2,&beta,&cost) << endl;
  logger << LOG_DEBUG_POINT << "h2_tt_(" << i4 << ", " << i2 << ") = " << setw(23) << setprecision(17) << showpoint << h2_tt_(&i4,&i2,&beta,&cost) << endl;
  logger << LOG_DEBUG_POINT << "h2_tt_(" << i3 << ", " << i4 << ") = " << setw(23) << setprecision(17) << showpoint << h2_tt_(&i3,&i4,&beta,&cost) << endl;
  logger << LOG_DEBUG_POINT << "h2_tttt_(" << i3 << ", " << i3 << ", " << i3 << ", " << i3 << ") = " << setw(23) << setprecision(17) << showpoint << h2_tttt_(&i3,&i3,&i3,&i3,&beta,&cost) << endl;
  logger << LOG_DEBUG_POINT << "h2_tttt_(" << i3 << ", " << i4 << ", " << i3 << ", " << i4 << ") = " << setw(23) << setprecision(17) << showpoint << h2_tttt_(&i3,&i4,&i3,&i4,&beta,&cost) << endl;
  logger << LOG_DEBUG_POINT << "h2_tttt_(" << i3 << ", " << i4 << ", " << i3 << ", " << i3 << ") = " << setw(23) << setprecision(17) << showpoint << h2_tttt_(&i3,&i4,&i3,&i3,&beta,&cost) << endl;

  logger << LOG_DEBUG_POINT << "osi_QT_H2_delta(1) = " << osi_QT_H2_delta << endl;

  //C I add the term for the scale mismatch between mt and M in 2-loop virtuals
  double logM2mt2 = log(sn / m2_HQ);

  osi_QT_H2_delta += 2 * logM2mt2 * (b1 * M0M0 + b0 * M1M0 - b0*b0 * M0M0 * logM2mt2);

  logger << LOG_DEBUG_POINT << "osi_QT_H2_delta(2) = " << osi_QT_H2_delta << endl;

  //C I add now the additional term due to the fact that the subtraction
  //C is performed for mu=mt instead of mu=M
  //C I start by defining gama1g and gama2g and K

  /////    double gama1 = -b0 / 2;
  /////    double gama2 = (2 * C_F * nf + (C_A * (256. / 27 - (2 * pi2) / 9) * nf) / 2
  /////		      + C_A*C_A * (-692. / 27 + (11. * pi2) / 18 + 2 * zeta3)) / 16;
  double KK = (67. / 18 - pi2_6) * C_A - 5. / 9 * nf;
  double gamaN = (20. * C_F * nf) / 9 + C_A * C_F * (-98. / 9 + (2 * pi2) / 3 - 4 * zeta3);

  //C I define v and b34
  ///      double complex b34,expb34,li2b34,li3b34,HPL2,HPL3,CothB34,aux1C,
  ///     &               aux2C,M,Mcc(4,4),aux3C

    double vv = sqrt(1. - pow(m2_HQ / (sn / 2 - m2_HQ), 2));
		     //dot(p,3,4)**2.d0);
    double_complex b34 = ri * pi + 1. / 2 * log((1. - vv) / (1. + vv));
    double_complex expB34 = exp(-2. * b34);
    int i0 = 0;
    double_complex li2b34 = hpl2_(&i0, &i1, &expB34);
    double_complex li3b34 = hpl3_(&i0, &i0, &i1, &expB34);
    double_complex CothB34 = (exp(-b34) + exp(b34)) / (-exp(-b34) + exp(b34));

    // C I write now all together the terms
    // C -4*b0*logM2mt2*Re(<M0|U1|M0>)+2*Re(<M0|U1|M1>)

    // C I define an auxiliary matrix element:
    // C TODO CHECK IMAGINARY PART
    double_complex M = -4. * b0 * logM2mt2 * M0M0 + M1M0 + ri * imM1M0;
    vector<vector<double_complex> > Mcc(5, vector<double_complex> (5, 0.));
    for (int i = 1; i < 5; i++){
      for (int j = 1; j < 5; j++){
        Mcc[i][j] = -4. * b0 * logM2mt2 * M0M0cc[i][j] + M1M0cc[i][j] + ri * imM1M0cc[i][j];
      }
    }

    double_complex aux1C = (logM2mt2
			    * (-2 * C_F * M + 4 * gama1 * M
			       + Ci * logM2mt2 * M - 2 * Ci * ri * M * pi
			       - 2. * b34 * CothB34 * Mcc[3][4])) / 2.;
    for (int i = 3; i < 5; i++){
      for (int j = 1; j < 3; j++){
	aux1C += ((logM2mt2 * (-(ri * pi)
			       + log((m2_HQ * sn) / (4. * pow(osi_p_parton[0][i] * osi_p_parton[0][j], 2)))) * Mcc[i][j]) / 2.);
      }
    }

    logger << LOG_DEBUG_POINT << "aux1C = " << aux1C << endl;
    logger << LOG_DEBUG_POINT << "real(aux1C) = " << real(aux1C) << endl;

    osi_QT_H2_delta += real(aux1C);

    logger << LOG_DEBUG_POINT << "osi_QT_H2_delta(3) = " << osi_QT_H2_delta << endl;

    //C I write now the term 2*Re(<M0|U2-U1^2|M0>)
    //    double temp_B4;
    for (int i = 3; i < 5; i++){
      for (int j = 1; j < 3; j++){
	temp_B4 = 0.;
	fourcorrelators_(&temp_B4, &j, &i, &i3, &i4, &temp_channel);
	logger << LOG_DEBUG_POINT << "fourcorrelators_(&temp_B4, " << j << ", " << i << ", " << i3 << ", " << i4 << ", " << temp_channel << "); -> temp_B4 = " << temp_B4 << endl;
        M0M04c[i][j][3][4] = temp_B4 * M0M0;
	for (int k = 3; k < 5; k++){
	  for (int l = 1; l < 3; l++){
	    temp_B4 = 0.;
	    fourcorrelators_(&temp_B4, &j, &i, &l, &k, &temp_channel);
	logger << LOG_DEBUG_POINT << "fourcorrelators_(&temp_B4, " << j << ", " << i << ", " << l << ", " << k << ", " << temp_channel << "); -> temp_B4 = " << temp_B4 << endl;
	    M0M04c[i][j][k][l] = temp_B4 * M0M0;
	  }
	}
      }
    }

    double_complex aux2C = (logM2mt2
			    * (-36. * b34*b34 * CothB34*CothB34 * logM2mt2*M0M0_3434
			       + 6. * M0M0 * (48. * gama2 + 3. * gamaN - 6. * b0 * C_F * logM2mt2
					      + 12. * b0 * gama1 * logM2mt2 + 6. * Ci * KK * logM2mt2
					      + 4. * b0 * Ci * logM2mt2 * logM2mt2 - 12. * Ci * ri * KK * pi
					      - 6. * b0 * Ci * ri * logM2mt2 * pi)
			       - 9. * logM2mt2 * M0M0 * pow(-2. * C_F + 4. * gama1 + Ci * (logM2mt2 - 2. * ri * pi), 2)
			       - 4. * (-18. * b34*b34 * C_A * (-1. + CothB34)
				       + 6. * pow(b34, 3) * C_A * (-1. + CothB34) * CothB34
				       + b34 * CothB34 * (9. * b0 * logM2mt2 + 18. * C_F * logM2mt2
							  - 36. * gama1 * logM2mt2 - 9. * Ci * logM2mt2 * logM2mt2
							  + 18. * Ci * ri * logM2mt2 * pi
							  + C_A * (67. - 6. * pi2
								   + 3. * CothB34 * (6. * li2b34 + pi2)) - 10. * nf)
				       + 3. * C_A * (pi2 + CothB34 * (6. * li2b34 - pi2)
						     + 6. * CothB34*CothB34 * (li3b34 - zeta3) + 6. * zeta3)
				       - 36. * b34 * C_A * CothB34 * log(1. - expB34)) * M0M0cc[3][4])) / 36.;

    for (int i = 3; i < 5; i++){
      for (int j = 1; j < 3; j++){
	aux2C +=
	  (logM2mt2 * (pi * (-9. * b0 * ri * logM2mt2 - 18. * C_F * ri * logM2mt2
			     + 36. * gama1 * ri * logM2mt2 + 9. * Ci * ri * logM2mt2 * logM2mt2
			     + 18. * Ci * logM2mt2 * pi + C_A * ri * (-67. + 3. * pi2)
			     + 10. * ri * nf)
		       + (9. * b0 * logM2mt2 + 18. * C_F * logM2mt2
			  - 36. * gama1 * logM2mt2 - 9. * Ci * logM2mt2 * logM2mt2
			  + 18. * Ci * ri * logM2mt2 * pi + C_A * (67. - 3. * pi2)
			  - 10. * nf) * log((m2_HQ * sn) / (4. * pow(osi_p_parton[0][i] * osi_p_parton[0][j], 2))))
	   * M0M0cc[i][j]) / 18.
	  + b34 * CothB34 * logM2mt2*logM2mt2
	  * (-(ri * pi) + log((m2_HQ * sn) / (4. * pow(osi_p_parton[0][i] * osi_p_parton[0][j], 2))))
	  * M0M04c[i][j][3][4];

	for (int k = 3; k < 5; k++){
	  for (int l = 1; l < 3; l++){
	    aux2C += (
		      (logM2mt2*logM2mt2
		       * (pi2 + log((m2_HQ * sn) / (4. * pow(osi_p_parton[0][i] * osi_p_parton[0][j], 2)))
			  * (2. * ri * pi - log((m2_HQ * sn) / (4. * pow(osi_p_parton[0][k] * osi_p_parton[0][l], 2)))))) / 4.) * M0M04c[i][j][k][l];
	  }
	}
      }
    }

    logger << LOG_DEBUG_POINT << "aux2C = " << aux2C << endl;
    logger << LOG_DEBUG_POINT << "real(aux2C) = " << real(aux2C) << endl;

    osi_QT_H2_delta += real(aux2C);

    logger << LOG_DEBUG_POINT << "osi_QT_H2_delta = " << osi_QT_H2_delta << endl;




  osi_QT_A0 = osi_VA_b_ME2;
  // Check if this if statements change the "rest" contribution !!!
  // Better: Avoid calculation of the ingredients if not finally needed for the chosen channels !!!
  //  if (oset.initial_diag || oset.initial_pdf_gq){
  osi_QT_A1 = osi_VA_V_ME2;
  osi_QT_H1_delta = osi_QT_A1 / 2 / (osi_alpha_S * inv2pi * osi_QT_A0);
  //  }
  //  else {osi_QT_H1_delta = 0.;}
  //  if (oset.initial_diag){
  //  osi_QT_A2 = ME2_L_L + ME2_L2_B;
  osi_QT_H2_delta = osi_QT_H2_delta / 4 / (osi_QT_A0);
  logger << LOG_DEBUG_POINT << "osi_QT_H2_delta / 4 / (osi_QT_A0) = " << osi_QT_H2_delta << endl;

    //  osi_QT_H2_delta = osi_QT_A2 / 4 / (osi_QT_A0);
  //pow(osi_alpha_S * inv2pi, 2) *

  //  }
  //  else {osi_QT_H2_delta = 0.;}


  //  osi_QT_H1_delta = 0.;
  //  osi_QT_H2_delta = 0.;

  // Catani scheme -> qT scheme

  // adapt alpha_S/Pi normalisation
  //  osi_QT_H1_delta /= 2;
  //  osi_QT_H2_delta /= 4;

  //  osi_QT_H1_delta /= (inv2pi * osi_alpha_S * osi_QT_A0);
  //  osi_QT_H2_delta /= (pow(inv2pi * osi_alpha_S, 2) * osi_QT_A0);


  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
