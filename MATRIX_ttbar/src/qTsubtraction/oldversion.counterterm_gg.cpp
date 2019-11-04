#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../include/definitions.observable.set.cxx"


double sigma12_gg(double pdf_factor_x1x2, double A1g)
 {
   return -0.5*A1g*pdf_factor_x1x2;
 }


double sigma11_gg(double pdf_factor_x1x2, double B1g, double tH1F)
 {
   return - (pdf_factor_x1x2*B1g + tH1F);
 }

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

double sigma24_gg(double pdf_factor_x1x2, double A1g)
 {
   // sig24=(A1g)**2/8*tdelta
   return 1./8 * A1g * A1g * pdf_factor_x1x2;
 }


double sigma23_gg(double pdf_factor_x1x2, double A1g, double sig11, double beta0)
 {
   // sig23=-beta0*A1g/3*tdelta-0.5d0*A1g*sig11                                                
   return -A1g * (beta0/3.0*pdf_factor_x1x2 + 0.5*sig11);
 }


double sigma22_gg(double pdf_factor_x1x2, double A1g, double B1g, double A2g, double beta0, double sig11, double tH1F, double H1full, double tgaga,  double LR, int nf)
 {
   return   -0.5  *A1g*(    H1full        - beta0*LR*pdf_factor_x1x2) - 0.5*(A2g*pdf_factor_x1x2 + (B1g-beta0)*sig11) + 0.5*B1g*tH1F     + 0.5  *tgaga;
 // sig22 = -0.5d0*A1g*((tH1st+LF*tH1stF) - beta0*LR) - 0.5*(A2g*tdelta          + (B1g-beta0)*sig11) + 0.5d0*B1g*tH1stF + 0.5d0*tgaga        
 }


double sigma21_gg( double pdf_factor_x1x2, double B1g, double A_F, double beta0, double sig11, double tH1, double tH1F, double H1full, double tgaga, double tcga, double tgamma2, double LR, double LF,  int nf )
 {
   Logger logger("sigma21_gg");
   logger << LOG_DEBUG_VERBOSE << "started" << endl;

   double Delta2gg = (9*(8.0/3+3*zeta3)-2.0/3*nf-2*nf)/4.0;
   //  const double B2g = 1.0/16*( (-64.0/3-24*zeta3)*9 + 16*nf + 16.0/3*nf ) + M_PI*M_PI*M_PI*beta0/2 ;
   double B2g = -2*Delta2gg + 0.5*beta0*(2.0/3*3*M_PI*M_PI + A_F);

   //   cout << "Delta2gg = " << Delta2gg << endl;
   //   cout << "B2g = " << B2g << endl;

   logger << LOG_DEBUG_VERBOSE << "-beta0*LR*sig11 = " << -beta0*LR*sig11 << endl;
   logger << LOG_DEBUG_VERBOSE << "-B1g*(tH1+LF*tH1F) = " << - B1g*(tH1+LF*tH1F) << endl;
   logger << LOG_DEBUG_VERBOSE << "-LF*tgaga = " << - LF*tgaga << endl;
   logger << LOG_DEBUG_VERBOSE << "-B2g*pdf_factor_x1x2 = " << - B2g*pdf_factor_x1x2 << endl;
   logger << LOG_DEBUG_VERBOSE << "+beta0*tH1 = " << + beta0*tH1 << endl;
   logger << LOG_DEBUG_VERBOSE << "-tcga = " << - tcga << endl;
   logger << LOG_DEBUG_VERBOSE << "-tgamma2 = " << - tgamma2 << endl;
   logger << LOG_DEBUG_VERBOSE << "-C1ggdelta(A_F)*tH1F = " << - C1ggdelta(A_F)*tH1F << endl;
   logger << LOG_DEBUG_VERBOSE << "-2*Delta2gg*pdf_factor_x1x2 = " << - 2*Delta2gg*pdf_factor_x1x2 << endl;
   logger << LOG_DEBUG_VERBOSE << "+2*beta0*LR*(B1g*pdf_factor_x1x2+tH1F) = " << +2*beta0*LR*(B1g*pdf_factor_x1x2+tH1F) << endl;


   return  -beta0*LR*sig11 - B1g*(tH1+LF*tH1F) - LF*tgaga - B2g*pdf_factor_x1x2 + beta0*tH1   - tcga - tgamma2 - C1ggdelta(A_F)*tH1F   - 2*Delta2gg*pdf_factor_x1x2+2*beta0*LR*(B1g*pdf_factor_x1x2+tH1F);


// sig21 = -beta0*LR*sig11 - B1g*H1full        - LF*tgaga - B2g*tdelta          + beta0*tH1st - tcga - tgamma2 - C1ggdelta     *tH1stF - 2*Delta2gg*tdelta 

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
 }


void calcIntermediateTerms_gg( double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_qz2, double pdf_factor_z1q, double pdf_factor_qq, int nf, double A_F, double beta0, double Kappa, double &tgaga, double &tcga, double &tgamma2 )
 {
  Logger logger("calcIntermediateTerms_gg");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
   /*        
     diff10 = - 3*[ dlog(xx10) * (fx1p(0)-fx10(0)*z1)/(1-z1) + D0int(xx10)*fx10(0) ] + beta0*fx10(0) - dlog(xx10)*fx1p(0)*Pggreg(z1)      
     diff20 = - 3*[ dlog(xx20) * (fx2p(0)-fx20(0)*z2)/(1-z2) + D0int(xx20)*fx20(0) ] + beta0*fx20(0) - dlog(xx20)*fx2p(0)*Pggreg(z2)

     diff1f = - dlog(xx10)*Pgq(z1) * (fx1p(j)+fx1p(-j))
     diff2f = - dlog(xx20)*Pgq(z2) * (fx2p(j)+fx2p(-j))

     Second part: gamma*gamma terms   
     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ + Pijjk(z) + Deltaijjk delta(1-z)

     tgaga = tgaga + 2*msqb * (diff10*diff20 + diff1f*diff2f + diff10*diff2f + diff1f*diff20) 
   */

   double D0gggg = 6*beta0, D1gggg = 18.0, Deltagggg = beta0*beta0-3.0/2*M_PI*M_PI;

   tgaga = 0;

// tgaga = tgaga - dlog(xx10)*( (fx1p(0)*fx20(0)-fx10(0)*fx20(0)*xx10**beta)*(D0gggg/(1-z1)+D1gggg*dlog(1-z1)/(1-z1)) + fx1p(0)*fx20(0)*(Pggggreg(z1)+Pgqqg(z1)) )      + (Deltagggg-D0gggg*D0int(xx10)-D1gggg*D1int(xx10))*fx10(0)*fx20(0)
   tgaga = tgaga -  log( x1 )*( (pdf_factor_z1x2-pdf_factor_x1x2* z1 )*(D0gggg/(1-z1)+D1gggg* log(1-z1)/(1-z1)) + pdf_factor_z1x2*(Pggggreg(z1, beta0)+Pgqqg(z1, nf)) ) + (Deltagggg-D0gggg*D0int( x1 )-D1gggg*D1int( x1 ))*pdf_factor_x1x2;
// tgaga = tgaga - dlog(xx20)*( (fx10(0)*fx2p(0)-fx10(0)*fx20(0)*xx20**alfa)*(D0gggg/(1-z2)+D1gggg*dlog(1-z2)/(1-z2)) + fx10(0)*fx2p(0)*(Pggggreg(z2)+Pgqqg(z2)) )      + (Deltagggg-D0gggg*D0int(xx20)-D1gggg*D1int(xx20))*fx10(0)*fx20(0)
   tgaga = tgaga -  log( x2 )*( (pdf_factor_x1z2-pdf_factor_x1x2* z2 )*(D0gggg/(1-z2)+D1gggg* log(1-z2)/(1-z2)) + pdf_factor_x1z2*(Pggggreg(z2, beta0)+Pgqqg(z2, nf)) ) + (Deltagggg-D0gggg*D0int( x2 )-D1gggg*D1int( x2 ))*pdf_factor_x1x2;
   
// tgaga = tgaga - dlog(xx10) * (Pgqqq(z1)+Pgggq(z1))        * (fx1p(j)+fx1p(-j)) * fx20(0)*msqb   
   tgaga = tgaga -  log( x1 ) * (Pgqqq(z1)+Pgggq(z1, beta0)) * pdf_factor_qx2;
// tgaga = tgaga - dlog(xx20) * (Pgqqq(z2)+Pgggq(z2))        * (fx2p(j)+fx2p(-j)) * fx10(0)*msqb  
   tgaga = tgaga -  log( x2 ) * (Pgqqq(z2)+Pgggq(z2, beta0)) * pdf_factor_x1q;

   double tx1 = 3.0*log(x1)*z1/(1-z1) - 3*D0int(x1) + beta0,  tz1 = -log(x1)*( 3.0/(1-z1) + Pggreg(z1) ) ;
   double tx2 = 3.0*log(x2)*z2/(1-z2) - 3*D0int(x2) + beta0,  tz2 = -log(x2)*( 3.0/(1-z2) + Pggreg(z2) ) ;
   double diff10_diff20 = pdf_factor_x1x2*tx1*tx2 + pdf_factor_x1z2*tx1*tz2 + pdf_factor_z1x2*tz1*tx2 + pdf_factor_z1z2*tz1*tz2 ;
   double diff10_diff2f =                 -( pdf_factor_x1q*tx1 + pdf_factor_z1q*tz1 )*log(x2)*Pgq(z2);
   double diff1f_diff20 = -log(x1)*Pgq(z1)*( pdf_factor_qx2*tx2 + pdf_factor_qz2*tz2 );
   double diff1f_diff2f =  log(x1)*Pgq(z1)                                            *log(x2)*Pgq(z2)*pdf_factor_qq;

   tgaga = tgaga + 2*(diff10_diff20 + diff10_diff2f + diff1f_diff20 + diff1f_diff2f);

   //------------------------------------------------------------------------------------------------------------------------------------------------

   /*
      diffg10 = -dlog(xx10) * ((fx1p(0)-fx10(0)*xx10**beta)*3d0/(1-z1) + Pggreg(z1)*fx1p(0)) + fx10(0)*(beta0-3*D0int(xx10))
      diffg1f = -dlog(xx10)*(fx1p(j)+fx1p(-j))*Pgq(z1)
      diffc20 = C1ggdelta*fx20(0)       
      diffc2f = -dlog(xx20)*(fx2p(j)+fx2p(-j))*Cgq(z2) 
      
      diffc10 = C1ggdelta*fx10(0)
      diffc1f = -dlog(xx10)*(fx1p(j)+fx1p(-j))*Cgq(z1)
      diffg20 = -dlog(xx20) * ((fx2p(0)-fx20(0)*xx20**alfa)*3d0/(1-z2) + Pggreg(z2)*fx2p(0)) + fx20(0)*(beta0-3*D0int(xx20))
      diffg2f = -dlog(xx20)*(fx2p(j)+fx2p(-j))*Pgq(z2) 
   
c    gamma first leg, C second leg
      tcga = tcga + msqb*(diffg10*diffc20 + diffg1f*diffc2f + diffg10*diffc2f + diffg1f*diffc20) = (diffg10 + diffg1f) * (diffc20 + diffc2f)
c    gamma second leg, C first leg
      tcga = tcga + msqb*(diffc10*diffg20 + diffc1f*diffg2f + diffc1f*diffg20 + diffc10*diffg2f) = (diffc10 + diffc1f) * (diffg20 + diffg2f)
   */

   tcga = 0;
   
// tcga = tcga - CgqPqg(z1)    *dlog(xx10)*flgq* fx1p(0)*fx20(0) *msqb
   tcga = tcga - CgqPqg(z1, nf)* log( x1 )     * pdf_factor_z1x2;
   logger << LOG_DEBUG_VERBOSE << "tcga: z1x2:   " << - CgqPqg(z1, nf)* log( x1 ) << endl;

// tcga = tcga - CgqPqg(z2)    *dlog(xx20)*flgq* fx2p(0)*fx10(0) *msqb
   tcga = tcga - CgqPqg(z2, nf)* log( x2 )     * pdf_factor_x1z2;
   logger << LOG_DEBUG_VERBOSE << "tcga: x1z2:   " << - CgqPqg(z2, nf)* log( x2 ) << endl;
 
// tcga = tcga - CgqPqq(z1)*dlog(xx10)*(fx1p(j)+fx1p(-j))*fx20(0)*msqb
   tcga = tcga - CgqPqq(z1)* log( x1 )*pdf_factor_qx2;
   logger << LOG_DEBUG_VERBOSE << "tcga: qx2:   " << - CgqPqq(z1)* log( x1 ) << endl;

// tcga = tcga - CgqPqq(z2)*dlog(xx20)*(fx2p(j)+fx2p(-j))*fx10(0)*msqb
   tcga = tcga - CgqPqq(z2)* log( x2 )*pdf_factor_x1q;
   logger << LOG_DEBUG_VERBOSE << "tcga: x1q:   " << - CgqPqq(z2)* log( x2 ) << endl;

   double diffg10_diffc20 = ( pdf_factor_z1x2*tz1 + pdf_factor_x1x2*tx1 )*C1ggdelta(A_F) ;
   logger << LOG_DEBUG_VERBOSE << "tcga: C1ggdelta(A_F):   " << C1ggdelta(A_F) << endl;
   logger << LOG_DEBUG_VERBOSE << "tcga: z1x2:   " << tz1*C1ggdelta(A_F) << endl;
   logger << LOG_DEBUG_VERBOSE << "tcga: x1x2:   " << tx1*C1ggdelta(A_F) << endl;

   double diffg1f_diffc2f = pdf_factor_qq*log(x1)*log(x2)*Pgq(z1)*Cgq(z2) ;
   logger << LOG_DEBUG_VERBOSE << "tcga: qq:   " << log(x1)*log(x2)*Pgq(z1)*Cgq(z2) << endl;

   double diffg10_diffc2f = -log(x2)*Cgq(z2)*( pdf_factor_z1q*tz1 + pdf_factor_x1q*tx1 ) ;
   logger << LOG_DEBUG_VERBOSE << "tcga: z1q:   " << -log(x2)*Cgq(z2)*tz1 << endl;
   logger << LOG_DEBUG_VERBOSE << "tcga: x1q:   " << -log(x2)*Cgq(z2)*tx1 << endl;

   double diffg1f_diffc20 = -log(x1)*Pgq(z1)*C1ggdelta(A_F)*pdf_factor_qx2;
   logger << LOG_DEBUG_VERBOSE << "tcga: C1ggdelta(A_F):   " << C1ggdelta(A_F) << endl;
   ///   logger << LOG_DEBUG_VERBOSE << "tcga; qx2 missing in new implementation: !!!" << endl;
   logger << LOG_DEBUG_VERBOSE << "tcga: qx2:   " << -log(x1)*Pgq(z1)*C1ggdelta(A_F) << endl;

   tcga = tcga + diffg10_diffc20 + diffg1f_diffc2f + diffg10_diffc2f + diffg1f_diffc20;
   
   double diffc10_diffg20 = C1ggdelta(A_F)*( pdf_factor_x1z2*tz2 + pdf_factor_x1x2*tx2 ) ;
   logger << LOG_DEBUG_VERBOSE << "tcga: C1ggdelta(A_F):   " << C1ggdelta(A_F) << endl;
   logger << LOG_DEBUG_VERBOSE << "tcga: x1z2:   " << tz2*C1ggdelta(A_F) << endl;
   logger << LOG_DEBUG_VERBOSE << "tcga: x1x2:   " << tx2*C1ggdelta(A_F) << endl;

   double diffc1f_diffg2f = pdf_factor_qq*log(x1)*log(x2)*Cgq(z1)*Pgq(z2) ;
   logger << LOG_DEBUG_VERBOSE << "tcga: qq:   " << log(x1)*log(x2)*Cgq(z1)*Pgq(z2) << endl;

   double diffc1f_diffg20 = -log(x1)*Cgq(z1)*(pdf_factor_qz2*tz2 + pdf_factor_qx2*tx2) ;
   logger << LOG_DEBUG_VERBOSE << "tcga: qz2:   " << -log(x1)*Cgq(z1)*tz2 << endl;
   logger << LOG_DEBUG_VERBOSE << "tcga: qx2:   " << -log(x1)*Cgq(z1)*tx2 << endl;

   double diffc10_diffg2f = -C1ggdelta(A_F)*log(x2)*Pgq(z2)*pdf_factor_x1q ;
   logger << LOG_DEBUG_VERBOSE << "tcga: C1ggdelta(A_F):   " << C1ggdelta(A_F) << endl;
   ///   logger << LOG_DEBUG_VERBOSE << "tcga; x1q missing in new implementation: !!!" << endl;
   logger << LOG_DEBUG_VERBOSE << "tcga: x1q:   " << -log(x2)*Pgq(z2)*C1ggdelta(A_F) << endl;



   tcga = tcga + diffc10_diffg20 + diffc1f_diffg2f + diffc1f_diffg20 + diffc10_diffg2f;

   //---------------------------------------------------------------------------------------------------------------------------------------------------

   tgamma2 = 0;
   
//  gamma2: diagonal part
//   First leg
// tgamma2 = tgamma2 - ( 1.5d0*Kappa*( dlog(xx10)*(fx1p(0)*fx20(0)-fx10(0)*fx20(0)*xx10**beta)/(1-z1) + D0int(xx10)*fx10(0)*fx20(0)) + dlog(xx10)*P2gg(z1)    *fx1p(0)*fx20(0) ) *flgq*msqb
   tgamma2 = tgamma2 - ( 1.5  *Kappa*(  log( x1 )*(pdf_factor_z1x2-pdf_factor_x1x2* z1       )/(1-z1) + D0int( x1 )*pdf_factor_x1x2) +  log( x1 )*P2gg(z1, nf)*pdf_factor_z1x2 );  

//   Second leg
// tgamma2 = tgamma2 - ( 1.5d0*Kappa*( dlog(xx20)*(fx10(0)*fx2p(0)-fx10(0)*fx20(0)*xx20**alfa)/(1-z2) + D0int(xx20)*fx10(0)*fx20(0)) + dlog(xx20)*P2gg(z2)    *fx10(0)*fx2p(0) ) *flgq*msqb 
   tgamma2 = tgamma2 - ( 1.5  *Kappa*(  log( x2 )*(pdf_factor_x1z2-pdf_factor_x1x2* z2       )/(1-z2) + D0int( x2 )*pdf_factor_x1x2) +  log( x2 )*P2gg(z2, nf)*pdf_factor_x1z2 );

//  gamma2: qg channel
//   First leg
// tgamma2 = tgamma2 - dlog(xx10)*P2gq(z1)    * (fx1p(j)+fx1p(-j))*fx20(0) *msqb 
   tgamma2 = tgamma2 -  log( x1 )*P2gq(z1, nf)* pdf_factor_qx2;
   
//  Second leg
// tgamma2 = tgamma2 - dlog(xx20)*P2gq(z2)    * fx10(0)*(fx2p(j)+fx2p(-j)) *msqb
   tgamma2 = tgamma2 -  log( x2 )*P2gq(z2, nf)* pdf_factor_x1q;

   logger << LOG_DEBUG_VERBOSE << "tgaga   = " << tgaga << endl;
   logger << LOG_DEBUG_VERBOSE << "tcga    = " << tcga << endl;
   logger << LOG_DEBUG_VERBOSE << "tgamma2 = " << tgamma2 << endl;
   //   tcga = 0.; // !!!

   logger << LOG_DEBUG_VERBOSE << "finished" << endl;
 }

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

void observable_set::calculate_sigma_gg( double qt2, double q2, double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_qz2, double pdf_factor_z1q, double pdf_factor_qq, int order, double A_F, double LR, double LF, vector<double> & sigma_qT )
 {
  Logger logger("observable_set::calculate_sigma_gg");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

   static int init = 1;
   static vector<double> I1_int, I2_int, I3_int, I4_int;
   if(init == 1)
     {
       I1_int.resize(n_qTcut); // n_qTcut = 200;
       I2_int.resize(n_qTcut);
       I3_int.resize(n_qTcut);
       I4_int.resize(n_qTcut);
       computeItildaIntegrals(min_qTcut/100, step_qTcut/100, n_qTcut, I1_int, I2_int, I3_int, I4_int);
       init = 0;
     }

   double beta0 = (33.0-2.0*N_f)/12 ;
   double tH1F  = H1F_gg(z1, z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_qx2, pdf_factor_x1q, beta0);
   
   // first order resummation coefficients
   double A1g = 3.0;
   double B1g = -2*beta0;
   
   double sig12 = sigma12_gg( pdf_factor_x1x2, A1g);
   double sig11 = sigma11_gg( pdf_factor_x1x2, B1g, tH1F);
   
   if (order == 3) {
     for (int i=0; i<n_qTcut; i++) {
       sigma_qT[i] = alpha_S/M_PI*(sig12*I2_int[i] + sig11*I1_int[i]);
     }

     //	if (sd == 0 && ss == 1){
	  logger << LOG_DEBUG << "global_factor  = " << alpha_S/M_PI << endl;
	  logger << LOG_DEBUG << "I2_int[0] = " << I2_int[0] << endl;
	  logger << LOG_DEBUG << "sig12 = " << sig12 << endl;
	  logger << LOG_DEBUG << "I1_int[0] = " << I1_int[0] << endl;
	  logger << LOG_DEBUG << "sig11 = " << sig11 << endl;
	  //	}

     return;
   }

   double Kappa = 67.0/6 - M_PI*M_PI/2 - 5.0/9*N_f;
   double A2g   = 0.5 * A1g * Kappa;

   double C_A      = 3;
   double H1_delta = M_PI*M_PI/6 * C_A + 0.5 * A_F;

//   H1_delta=0;
//   A_F=-2*M_PI*M_PI/6 * C_A;

   double tH1      = H1_M_gg(z1, z2, x1, x2, pdf_factor_x1x2, pdf_factor_qx2, pdf_factor_x1q, H1_delta );
   double H1full   = H1_gg(pdf_factor_x1x2, tH1, tH1F, LR, LF, beta0);

   double tgaga, tcga, tgamma2;
   calcIntermediateTerms_gg(z1, z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_qx2, pdf_factor_x1q, pdf_factor_z1z2, pdf_factor_qz2, pdf_factor_z1q, pdf_factor_qq, N_f, A_F, beta0, Kappa, tgaga, tcga, tgamma2 );
   
   //   tgaga = 0.; // !!!
   //   tcga = 0.; // !!!
   //   tgamma2 = 0.; // !!!


   double sig24 = sigma24_gg( pdf_factor_x1x2, A1g );
   double sig23 = sigma23_gg( pdf_factor_x1x2, A1g, sig11, beta0 );
   double sig22 = sigma22_gg( pdf_factor_x1x2, A1g, B1g, A2g, beta0, sig11, tH1F, H1full, tgaga, LR, N_f);
   double sig21 = sigma21_gg( pdf_factor_x1x2, B1g, A_F, beta0, sig11, tH1, tH1F, H1full, tgaga, tcga, tgamma2, LR, LF, N_f);

   for (int i=0; i<n_qTcut; i++) {
     sigma_qT[i] = pow(alpha_S/M_PI,2)*(sig24*I4_int[i] + sig23*I3_int[i] + sig22*I2_int[i] + sig21*I1_int[i]);
     //      cout << LOG_DEBUG << "QT:   sigma_qT[" << setw(3) << i << "] = " << sigma_qT[i] << endl;
   }
     
   //logger << LOG_DEBUG_VERBOSE << "z: " << z1 << ", " << z2 << endl;
   //logger << LOG_DEBUG_VERBOSE << "x: " << x1 << ", " << x2 << endl;
      logger << LOG_DEBUG_VERBOSE << "H1_delta = " << H1_delta << endl; 
      logger << LOG_DEBUG_VERBOSE << "A_F      = " << A_F      << endl; 
      logger << LOG_DEBUG_VERBOSE << "LR       = " << LR       << endl; 
      logger << LOG_DEBUG_VERBOSE << "LF       = " << LF       << endl; 
      logger << LOG_DEBUG_VERBOSE << "sig11    = " << sig11    << endl; 
      logger << LOG_DEBUG_VERBOSE << "tH1      = " << tH1      << endl;
      logger << LOG_DEBUG_VERBOSE << "tH1F     = " << tH1F     << endl;
      logger << LOG_DEBUG_VERBOSE << "H1full   = " << H1full   << endl;
      logger << LOG_DEBUG_VERBOSE << "tgaga    = " << tgaga    << endl;
      logger << LOG_DEBUG_VERBOSE << "tcga     = " << tcga     << endl;
      logger << LOG_DEBUG_VERBOSE << "tgamma2  = " << tgamma2  << endl;
      logger << LOG_DEBUG_VERBOSE << "sig21    = " << sig21    << endl; 
      logger << LOG_DEBUG_VERBOSE << "sig22    = " << sig22    << endl;
      logger << LOG_DEBUG_VERBOSE << "sig23    = " << sig23    << endl;
      logger << LOG_DEBUG_VERBOSE << "sig24    = " << sig24    << endl;
   //   exit(0); 
 }


double H1F_gg( double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_qx2, double pdf_factor_x1q, double beta0 )
 {
   double tH1F = 0.0;
   
   // tH1stF = - dlog(xx10)*( (fx1p(0)*fx20(0) - fx10(0)*fx20(0)*xx10**beta)*3/(1-z1) + fx1p(0)*fx20(0)*Pggreg(z1) + (fx1p(j) + fx1p(-j))*fx20(0) * Pgq(z1) ) - 3*D0int(xx10)*fx10(0)*fx20(0)   -                                
   tH1F = tH1F -  log( x1 )*( (pdf_factor_z1x2 - pdf_factor_x1x2*  z1      )*3/(1-z1) + pdf_factor_z1x2*Pggreg(z1) +  pdf_factor_qx2              * Pgq(z1) ) - 3*D0int(x1)  *pdf_factor_x1x2;   

   
   //    cout << "gg xx leg 1: pdf_factor_x1x2   " << -  log( x1 )*( (-  z1      )*3/(1-z1) ) - 3*D0int(x1)  + beta0 << endl;
   //    cout << "gg zx leg 1: pdf_factor_z1x2   " << -  log( x1 )*( 3/(1-z1) + Pggreg(z1)) << endl;
   //    cout << "qg zx leg 1: pdf_factor_qx2   " << -  log( x1 )*(Pgq(z1) )  << endl;
   

   //          - dlog(xx20)*( (fx10(0)*fx2p(0) - fx10(0)*fx20(0)*xx20**alfa)*3/(1-z2) + fx10(0)*fx2p(0)*Pggreg(z2) + (fx2p(j) + fx2p(-j))*fx10(0) * Pgq(z2) ) - 3*D0int(xx20)*fx10(0)*fx20(0)   +                     
    tH1F = tH1F -  log( x2 )*( (pdf_factor_x1z2 - pdf_factor_x1x2*  z2      )*3/(1-z2) + pdf_factor_x1z2*Pggreg(z2) +  pdf_factor_x1q              * Pgq(z2) ) - 3*D0int(x2)  *pdf_factor_x1x2;

   //          + 2*beta0*fx10(0)*fx20(0)                                                                                                                                                                                 
    //    cout << "gg xx leg 2: pdf_factor_x1x2   " << -  log( x2 )*( (-  z2      )*3/(1-z2) ) - 3*D0int(x2)  + beta0 << endl;
    //    cout << "gg zx leg 2: pdf_factor_z1x2   " << -  log( x2 )*( 3/(1-z2) + Pggreg(z2)) << endl;
    //    cout << "qg zx leg 2: pdf_factor_qx2   " << -  log( x2 )*(Pgq(z2) )  << endl;

   tH1F = tH1F + 2*beta0*pdf_factor_x1x2;

   return tH1F;
 }


double H1_M_gg( double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_qx2, double pdf_factor_x1q, double H1_delta )
 {
   double tH1 = 0;

   // tH1st=tH1st + 2*C1ggdelta*fx10(0)*fx20(0)*msqb
   tH1 = tH1      + H1_delta   *pdf_factor_x1x2;
   // tH1st=tH1st - dlog(xx10)*(fx1p(j)+fx1p(-j))*fx20(0)*Cgq(z1)*msqb
   tH1 = tH1      -  log( x1 )*pdf_factor_qx2            *Cgq(z1);
   // tH1st=tH1st - dlog(xx20)*(fx2p(j)+fx2p(-j))*fx10(0)*Cgq(z2)*msqb
   tH1 = tH1      -  log( x2 )*pdf_factor_x1q            *Cgq(z2);
   
//    cout << "C1delta=" << 0.5*H1_delta << endl;
   
   return tH1;
 }


double H1_gg( double pdf_factor_x1x2, double tH1, double tH1F, double LR, double LF, double beta0 )
 {
   return  tH1 + LF*tH1F - 2*beta0*LR*pdf_factor_x1x2;
 }


double H2_M_gg( double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2,  double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_qz2, double pdf_factor_z1q, double pdf_factor_qq, int nf, double H1_delta, double H2_delta )
{
  double H2ggD0 = -101.0/3 + 14.0*nf/9 + 63.0/2*zeta3 ;
  double tH2    = 0;
  
  //  gg channel  
  tH2 = tH2 + H2_delta*pdf_factor_x1x2;
  
  tH2 = tH2 - 0.5*( log(x1)*(pdf_factor_z1x2-pdf_factor_x1x2*z1)*H2ggD0/(1-z1) + H2ggD0*D0int(x1)*pdf_factor_x1x2 );
  tH2 = tH2 - 0.5*( log(x2)*(pdf_factor_x1z2-pdf_factor_x1x2*z2)*H2ggD0/(1-z2) + H2ggD0*D0int(x2)*pdf_factor_x1x2 );

  tH2 = tH2 - log(x1)*pdf_factor_z1x2*C2ggreg(z1,nf); 
  tH2 = tH2 - log(x2)*pdf_factor_x1z2*C2ggreg(z2,nf);

  tH2 = tH2 + log(x1)*log(x2)*Ggg(z1)*Ggg(z2)*pdf_factor_z1z2;  // tH2sp = msqb*dlog(xx10)*dlog(xx20)*fx1p(0)*fx2p(0)*Ggg(z1)*Ggg(z2)*flgq
  
  //  qg channel  
  tH2 = tH2 - log(x1)*pdf_factor_qx2*C2gq(z1,nf);               // tH2st = tH2st+(fx1p(j)+fx1p(-j))*H2stgq(z1)*(-dlog(xx10))*fx20(0)*msqb
  tH2 = tH2 - log(x2)*pdf_factor_x1q*C2gq(z2,nf);               // tH2st = tH2st+(fx2p(j)+fx2p(-j))*H2stgq(z2)*(-dlog(xx20))*fx10(0)*msqb
  tH2 = tH2 - log(x1)*H1_delta*pdf_factor_qx2*Cgq(z1);          // doesn't appear explicitly in HHNLO
  tH2 = tH2 - log(x2)*H1_delta*pdf_factor_x1q*Cgq(z2);          // doesn't appear explicitly in HHNLO 
  tH2 = tH2 + log(x1)*log(x2)*Ggg(z1)*Ggq(z2)*pdf_factor_z1q;   // tH2sp=tH2sp-msqb*dlog(xx10)*fx1p(0)*Ggg(z1)*spgq2
  tH2 = tH2 + log(x1)*log(x2)*Ggq(z1)*Ggg(z2)*pdf_factor_qz2;   // tH2sp=tH2sp-msqb*dlog(xx20)*fx2p(0)*Ggg(z2)*spgq1
  
  //  qq channel (C1C1 and G1G1)
  tH2 = tH2 + log(x1)*log(x2)*Cgq(z1)*Cgq(z2)*pdf_factor_qq;    // tH2st = tH2st+msqb*diffc1f*diffc2f*flgq    
  tH2 = tH2 + log(x1)*log(x2)*Ggq(z1)*Ggq(z2)*pdf_factor_qq;    // tH2sp = tH2sp+msqb*spgq1*spgq2*flgq
  
  return tH2;
}


double H2_gg( double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_qz2, double pdf_factor_z1q, double pdf_factor_qq, int nf,  double H1_delta, double H2_delta,  double LR, double LF, double tH1, double tH1F )
 {
  Logger logger("H2_gg");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

   double beta0 = (33.0  - 2.0 *nf)/12;
   double beta1 = (153.0 - 19.0*nf)/24;
   double Kappa = 67.0/6 - M_PI*M_PI/2 - 5.0/9*nf;
   double Delta2gg = (9*(8.0/3+3*zeta3)-2.0/3*nf-2*nf)/4.0;

   double C_A = 3;

   double tgaga, tcga, tgamma2;
   double A_F = 2 * (H1_delta - M_PI*M_PI/6 * C_A);;
   calcIntermediateTerms_gg(z1, z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_qx2, pdf_factor_x1q, pdf_factor_z1z2, pdf_factor_qz2, pdf_factor_z1q, pdf_factor_qq, nf, A_F, beta0, Kappa, tgaga, tcga, tgamma2);

   double H = 0;

   H = H + H2_M_gg( z1, z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_qx2, pdf_factor_x1q, pdf_factor_z1z2, pdf_factor_qz2, pdf_factor_z1q, pdf_factor_qq, nf, H1_delta, H2_delta );



   // 20170523: In H2_qqx, there is an extra term:
   // tH2 += H1_delta*(tH1-H1_delta*pdf_factor_x1x2);
   // Is that term missing here ???



   // add scale dependence at NNLO  
   // xmsq = xmsq + asopi**2*( 0.5d0*beta0*LF**2*tH1stF + tgamma2*LF - 2*beta0*LR*(tH1st+LF*tH1stF-2*beta0*LR*tdelta        ) - beta0*LF*LR*tH1stF + LF*tcga - beta0*LR*tH1st+0.5d0*LF**2*tgaga - 2*(0.5d0*(beta0*LR)**2 +beta1*LR)*tdelta )
   H      = H    +          ( 0.5  *beta0*LF*LF*tH1F   + tgamma2*LF - 2*beta0*LR*(tH1  +LF*tH1F  -2*beta0*LR*pdf_factor_x1x2) - beta0*LF*LR*tH1F   + LF*tcga - beta0*LR*tH1  +0.5  *LF*LF*tgaga - 2*(0.5  *(beta0*LR)*(beta0*LR)+beta1*LR)*pdf_factor_x1x2) ;

  // Include missing delta term from C*gamma (no factor 2 here !)
  // xmsq = xmsq + asopi**2*(LF*C1ggdelta*tH1stF)
   H      = H    +           LF*C1ggdelta(A_F)*tH1F;

  // Include missing term from contact term in 2 loop AP
  // xmsq = xmsq + asopi**2*(2*Delta2gg*tdelta)*LF
   H      = H    +           2*Delta2gg*pdf_factor_x1x2*LF;
   

    logger << LOG_DEBUG_VERBOSE << setw(50) << "LF" << " = " << LF << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "LR" << " = " << LR << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "H2_M_gg" << " = " << H2_M_gg( z1, z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_qx2, pdf_factor_x1q, pdf_factor_z1z2, pdf_factor_qz2, pdf_factor_z1q, pdf_factor_qq, nf, H1_delta, H2_delta ) << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "+ 0.5*beta0*LF*LF*tH1F" << " = " << 0.5  *beta0*LF*LF*tH1F << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "+ tgamma2*LF" << " = " << tgamma2*LF << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "- 2*beta0*LR*(tH1  +LF*tH1F  -2*beta0*LR*pdf_factor_x1x2)" << " = " << - 2*beta0*LR*(tH1  +LF*tH1F  -2*beta0*LR*pdf_factor_x1x2) << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "- beta0*LF*LR*tH1F" << " = " << - beta0*LF*LR*tH1F << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "+ LF*tcga" << " = " << + LF*tcga << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "- beta0*LR*tH1" << " = " << - beta0*LR*tH1 << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "+ 0.5  *LF*LF*tgaga" << " = " << +0.5  *LF*LF*tgaga << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "- 2*(0.5  *(beta0*LR)*(beta0*LR)+beta1*LR)*pdf_factor_x1x2" << " = " << - 2*(0.5  *(beta0*LR)*(beta0*LR)+beta1*LR)*pdf_factor_x1x2 << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "LF*C1ggdelta(A_F)*tH1F" << " = " << LF*C1ggdelta(A_F)*tH1F << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "2*Delta2gg*pdf_factor_x1x2*LF" << " = " << 2*Delta2gg*pdf_factor_x1x2*LF << endl;


   return H;
 }


double observable_set::calculate_virtual_gg( double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_qz2, double pdf_factor_z1q, double pdf_factor_qq, int order, double H1_delta, double H2_delta, double LR, double LF ) 
 {                                                                                                                                                                                                                        
   double beta0 = (33.0-2.0*N_f)/12;                                                                         
   if(order == 2)       
     return pdf_factor_x1x2;
   
   double tH1   = H1_M_gg(z1, z2, x1, x2, pdf_factor_x1x2, pdf_factor_qx2, pdf_factor_x1q, H1_delta );
   double tH1F  = H1F_gg( z1, z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_qx2, pdf_factor_x1q, beta0);     

   if(order == 3)
     return alpha_S/M_PI * H1_gg( pdf_factor_x1x2, tH1, tH1F, LR, LF, beta0 );
    
   double nf = (double)(N_f);

   /*
   cout<<"z1 = "<<z1<<"; z2 = "<<z2<<endl;                                          
   cout<<"x1 = "<<x1<<"; x2 = "<<x2<<endl;                                                                                                                                                                                              
   cout<<"beta = "<<log(z1)/log(x1)<<endl;                                                                                                                                                                                              
   cout<<"alfa = "<<log(z2)/log(x2)<<endl;                                                                                                                                                                                         
   cout<<"LF = "<<LF<<endl;
   cout<<"LR = "<<LR<<endl;
   
   double x = pow(alpha_S/M_PI, 2) * H2_gg( z1, z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_qx2, pdf_factor_x1q, pdf_factor_z1z2, pdf_factor_qz2, pdf_factor_z1q, pdf_factor_qq, nf, H1_delta, H2_delta, LR, LF, tH1, tH1F );
   
   cout<<"x = "<<x<<endl;
   
    exit(0); 
   */
   return pow(alpha_S/M_PI,2) * H2_gg( z1, z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_qx2, pdf_factor_x1q, pdf_factor_z1z2, pdf_factor_qz2, pdf_factor_z1q, pdf_factor_qq, nf, H1_delta, H2_delta, LR, LF, tH1, tH1F ); 

 }
