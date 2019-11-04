#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../include/definitions.observable.set.cxx"

double A1q=4.0/3;
double B1q=-2;

double sigma12_qqbar(double pdf_factor_x1x2) {
  return -0.5 * A1q * pdf_factor_x1x2;
}

double sigma11_qqbar(double pdf_factor_x1x2, double tH1F, double LQ) {
  return -B1q * pdf_factor_x1x2 - tH1F - pdf_factor_x1x2 * A1q * LQ;
}

double sigma24_qqbar(double pdf_factor_x1x2) {
  return 1. / 8 * A1q * A1q * pdf_factor_x1x2;
}

double sigma23_qqbar(double pdf_factor_x1x2, double sig11, double beta0) {
  return -beta0 * A1q / 3 * pdf_factor_x1x2 - 0.5 * A1q * sig11;
}

double sigma22_qqbar(double pdf_factor_x1x2, double beta0, double sig11, double tH1, double tH1F, double H1full, double tgaga, double LR, double LF, double LQ, int nf) {
  double Kappa=67.0/6-(M_PI*M_PI)/2-5.0/9*nf;
  double A2q=0.5*A1q*Kappa;

  double sig22=0.5*(beta0*A1q*LR-A2q)*pdf_factor_x1x2-0.5*A1q*H1full-0.5*(B1q-beta0)*sig11+0.5*B1q*tH1F+0.5*tgaga;

  // Q dependence
  sig22 -= 0.5*A1q*LQ*beta0*pdf_factor_x1x2;
  sig22 += 0.5*A1q*LQ*tH1F;
  sig22 -= 0.5*A1q*LQ*sig11;

  return sig22;

//        sig22=0.5d0*(beta0*A1q*LR-A2q)*tdelta
//     &     -0.5d0*A1q*(tH1st+LF*tH1stF)
//     &     -0.5d0*(B1q-beta0)*sig11
//     &     +0.5d0*B1q*tH1stF
//     &     +0.5d0*tgaga
}

double sigma21_qqbar(double pdf_factor_x1x2, double beta0, double sig11, double tH1, double tH1F, double H1full, double tgaga, double tcga, double tgamma2, double LR, double LF, double LQ, double A_F, int nf) {
//  const double B2q=4.0/9*(M_PI*M_PI-3.0/4-12.0*zeta3)+(11.0/9*M_PI*M_PI-193.0/12+6*zeta3)+nf/6.0*(17.0/3-4.0/9*M_PI*M_PI);
//  B2q=4d0/9*(pi**2-3d0/4-12*zeta3)+(11d0/9*pi**2-193d0/12+6*zeta3)
//     & +nf/6d0*(17d0/3-4d0/9*pi**2)
  const double Delta2qq=(16.0/9*(3.0/8-M_PI*M_PI/2+6*zeta3)+4*(17.0/24+11.0*M_PI*M_PI/18-3*zeta3)-2.0/3*nf*(1.0/6+2*M_PI*M_PI/9))/4;
  const double B2q=-2*Delta2qq+0.5*beta0*(2.0/3*4.0/3*M_PI*M_PI+A_F);
//Delta2qq=16d0/9*(3d0/8-pi**2/2+6*zeta3)
//     &   +4*(17d0/24+11d0*pi**2/18-3*zeta3)-2d0/3*nf*(1d0/6+2*pi**2/9d0)
//      Delta2qq=Delta2qq/4d0

  double sig21 = -beta0*LR*sig11-B1q*H1full-LF*tgaga-B2q*pdf_factor_x1x2+beta0*tH1-tcga-tgamma2-C1qqdelta(A_F)*tH1F-2*Delta2qq*pdf_factor_x1x2;

  // Q dependence
  double Kappa=67.0/6-(M_PI*M_PI)/2-5.0/9*nf;
  double A2q=0.5*A1q*Kappa;
  sig21 += LQ*sig11*beta0;
  sig21 -= LQ*H1full*A1q;
  sig21 += LQ*(tgaga+(B1q+0.5*A1q*LQ)*tH1F);
  sig21 -= LQ*A2q*pdf_factor_x1x2;

//  sig21 += LQ*(sig11*(beta0-B1q-0.5*A1q*LQ) - A2q*pdf_factor_x1x2 + A1q*(2*beta0*LR*pdf_factor_x1x2-tH1+(LQ-LF)*tH1F) + tgaga);

  return sig21;


//        sig21=-beta0*LR*sig11-B1q*(tH1st+LF*tH1stF)
//     &     -LF*tgaga-B2q*tdelta+beta0*tH1st-tcga-tgamma2
//
//c     Include missing delta term from C*gamma (no factor 2 here !)
//
//      sig21=sig21-C1qqdelta*tH1stF
//
//
//C     Include missing term from contact term in 2 loop AP
//
//      sig21=sig21-2*Delta2qq*tdelta
}

void calcIntermediateTerms(double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_z1z2, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb, int nf, double A_F, double beta0, double &tgaga, double &tcga, double &tgamma2) {
  double Kappa = 67.0/6 - (M_PI*M_PI)/2 - 5.0/9*nf;

//  CC    Coefficients of D0 and D1 in P*P (as/pi normalization)
  double D0qqqq = 8.0/3;
  double D1qqqq = 32.0/9;

//CC    Coefficients of delta(1-z) in P*P
  double Deltaqqqq = 4.0/9*(9.0/4 - 2*M_PI*M_PI/3);

  tgaga = 0.0;

// notes, Eq. (1) (checked twice)
//  double diffg1fg2f=pdf_factor_x1x2*(g_z1*g_z2*Pqq(z1)*Pqq(z2)*z1*z2+Pqqint(x1)*Pqqint(x2)-g_z1*Pqq(z1)*Pqqint(x2)*z1-g_z2*Pqq(z2)*Pqqint(x1)*z2);

//  (checked twice)
//  diffg1fg2f+=pdf_factor_z1x2*(-g_z1*g_z2*Pqq(z1)*Pqq(z2)*z2+g_z1*Pqq(z1)*Pqqint(x2));
//  diffg1fg2f+=pdf_factor_x1z2*(-g_z1*g_z2*Pqq(z1)*Pqq(z2)*z1+g_z2*Pqq(z2)*Pqqint(x1));
//  diffg1fg2f+=pdf_factor_z1z2*(g_z1*g_z2*Pqq(z1)*Pqq(z2));

  double diffg1fg2f = pdf_factor_x1x2*Pqqint(x1)*Pqqint(x2);
  diffg1fg2f += g_z1*Pqq(z1)*Pqqint(x2)* (-pdf_factor_x1x2*z1+pdf_factor_z1x2);
  diffg1fg2f += g_z2*Pqq(z2)*Pqqint(x1)* (-pdf_factor_x1x2*z2+pdf_factor_x1z2);
  diffg1fg2f += g_z1*g_z2*Pqq(z1)*Pqq(z2)* (pdf_factor_x1x2*z1*z2 - pdf_factor_z1x2*z2 - pdf_factor_x1z2*z1 + pdf_factor_z1z2);

  // Eq. (2) (checked twice)
  double diffg10g20 = pdf_factor_gg*g_z1*g_z2*Pqg(z1)*Pqg(z2);

  //Eq. (3) (checked twice)
  double diffg10g2f = g_z1*Pqg(z1)* ( pdf_factor_gx2*(-g_z2*z2*Pqq(z2)+Pqqint(x2)) + pdf_factor_gz2*g_z2*Pqq(z2) );
  double diffg1fg20 = g_z2*Pqg(z2)* ( pdf_factor_x1g*(-g_z1*z1*Pqq(z1)+Pqqint(x1)) + pdf_factor_z1g*g_z1*Pqq(z1) );

//  diffg1f = -dlog(xx10)*(fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1) - Pqqint(xx10)*fx10(j)
//  diffg10 = -dlog(xx10)*fx1p(0)*Pqg(z1)
//  diffg2f = -dlog(xx20)*(fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2) - Pqqint(xx20)*fx20(k)
//  diffg20 = -dlog(xx20)*fx2p(0)*Pqg(z2)

  // (checked twice)
  tgaga = 2.0*(diffg10g20 + diffg1fg2f + diffg10g2f + diffg1fg20);

//  tgaga=tgaga + 2*(flgq*diffg10*diffg20 + flgq*diffg1f*diffg2 + diffg10*diffg2f + diffg1f*diffg20)*msqc(j,k)

//CC     Second part: gamma*gamma terms

//C     First leg (checked twice)
//       diff1 = -dlog(xx10)*(flgq*(fx1p(j)         -    fx10(j)     *xx10**beta)*(D0qqqq/(1-z1)+D1qqqq*dlog(1-z1)/(1-z1)) + fx1p(j)        *Pqqqq(z1)*flgq + fx1p(0)       *(Pqqqg(z1)+Pqggg(z1)      ) ) 
  double diff1 = -g_z1      *(     (pdf_factor_z1x2 - pdf_factor_x1x2* z1       )*(D0qqqq/(1-z1)+D1qqqq* log(1-z1)/(1-z1)) + pdf_factor_z1x2*Pqqqq(z1)      + pdf_factor_gx2*(Pqqqg(z1)+Pqggg(z1,beta0)) );

  // +                     ( Deltaqqqq - D0qqqq*D0int(xx10) - D1qqqq*D1int(xx10))*fx10(j)*flgq
  diff1 += pdf_factor_x1x2*( Deltaqqqq - D0qqqq*D0int( x1 ) - D1qqqq*D1int( x1 ) );


//C    Second leg (checked twice)
  double diff2 = -g_z2* ( (pdf_factor_x1z2-pdf_factor_x1x2*z2)*(D0qqqq/(1-z2)+D1qqqq*log(1-z2)/(1-z2)) + pdf_factor_x1z2*Pqqqq(z2) + pdf_factor_x1g*(Pqqqg(z2)+Pqggg(z2,beta0)) );
  diff2 += pdf_factor_x1x2*( Deltaqqqq - D0qqqq*D0int(x2) - D1qqqq*D1int(x2) );

// diff2 = -dlog(xx20)*(flgq*(fx2p(k)-fx20(k)*xx20**alfa)*(D0qqqq/(1-z2)+D1qqqq*dlog(1-z2)/(1-z2)) + fx2p(k)*Pqqqq(z2)*flgq+fx2p(0)*(Pqqqg(z2)+Pqggg(z2))) + (Deltaqqqq-D0qqqq*D0int(xx20)-D1qqqq*D1int(xx20))*fx20(k)*flgq

  tgaga = tgaga + diff1 + diff2;

//C     Include Pqggq (checked twice)
  tgaga = tgaga - g_z1*pdf_factor_qx2*Pqggq(z1);
  tgaga = tgaga - g_z2*pdf_factor_x1q*Pqggq(z2);

//      do l=1,nf
//      diff1 = diff1 -dlog(xx10)*(fx1p(l)+fx1p(-l))*Pqggq(z1)*flgq
//      diff2 = diff2 -dlog(xx20)*(fx2p(l)+fx2p(-l))*Pqggq(z2)*flgq
//      enddo

//      tgaga=tgaga+diff1*fx20(k)*msqc(j,k)
//      tgaga=tgaga+diff2*fx10(j)*msqc(j,k)

//---------------------------------------------------------------------------------------------------

  tcga = 0.0;
  // (checked twice)
  double diffg10c20 = pdf_factor_gg*g_z1*g_z2*Pqg(z1)*Cqg(z2);
  // Eq. (7) (checked twice)
  double diffg1fc2f = pdf_factor_x1x2*C1qqdelta(A_F)*(g_z1*z1*Pqq(z1)-Pqqint(x1));
  diffg1fc2f += pdf_factor_z1x2*(-g_z1*Pqq(z1)*C1qqdelta(A_F));
  diffg1fc2f += pdf_factor_z1z2*(g_z1*g_z2*Cqq(z2)*Pqq(z1));
  diffg1fc2f += pdf_factor_x1z2*g_z2*Cqq(z2)*(-g_z1*Pqq(z1)*z1+Pqqint(x1));

  // (checked twice)
  double diffg10c2f = pdf_factor_gz2*g_z1*g_z2*Pqg(z1)*Cqq(z2);
  diffg10c2f += pdf_factor_gx2*(-g_z1*C1qqdelta(A_F)*Pqg(z1));

  // (checked twice)
  double diffg1fc20 = pdf_factor_x1g*g_z2*Cqg(z2)*(-g_z1*Pqq(z1)*z1+Pqqint(x1));
  diffg1fc20 += pdf_factor_z1g*g_z1*g_z2*Pqq(z1)*Cqg(z2);

  tcga = tcga + diffg10c20 + diffg1fc2f + diffg10c2f + diffg1fc20;

  //  cout << setw(19) << "contribution: " << "tcga 22 =   " << diffg1fc2f << endl;
  //  cout << setw(19) << "contribution: " << "tcga 24 =   " << diffg1fc20 << endl;
  //  cout << setw(19) << "contribution: " << "tcga 42 =   " << diffg10c2f << endl;
  //  cout << setw(19) << "contribution: " << "tcga 44 =   " << diffg10c20 << endl;

  //  C    Start  (C+C)*(gamma+gamma) term

//c    gamma first leg, C second leg

//      diffc2f = -dlog(xx20)*fx2p(k)*Cqq(z2) + C1qqdelta*fx20(k)

//      diffc20 = -dlog(xx20)*fx2p(0)*Cqg(z2)

//      tcga = tcga + msqc(j,k)*(flgq*diffg10*diffc20+flgq*diffg1f*diffc2f+diffg10*diffc2f+diffg1f*diffc20)

  // (checked twice)
  double diffc10g20 = pdf_factor_gg*g_z1*g_z2*Cqg(z1)*Pqg(z2);

  // Eq. (5)
    // (checked twice)
//  double diffc1fg2f=pdf_factor_x1x2*C1qqdelta(A_F)*(g_z2*z2*Pqq(z2)-Pqqint(x2));
//  diffc1fg2f+=pdf_factor_z1x2*g_z1*Cqq(z1)*(-g_z2*z2*Pqq(z2)+Pqqint(x2));
//  diffc1fg2f+=pdf_factor_z1z2*g_z1*g_z2*Cqq(z1)*Pqq(z2);
//  diffc1fg2f+=pdf_factor_x1z2*(-C1qqdelta(A_F)*g_z2*Pqq(z2));

  double diffc1fg2f = C1qqdelta(A_F)*g_z2*Pqq(z2) * (pdf_factor_x1x2*z2 - pdf_factor_x1z2);
  diffc1fg2f += g_z1*Cqq(z1)*g_z2*Pqq(z2) * (pdf_factor_z1z2 - pdf_factor_z1x2*z2);
  diffc1fg2f += Pqqint(x2)*(pdf_factor_z1x2*g_z1*Cqq(z1) - pdf_factor_x1x2*C1qqdelta(A_F));

  // (checked twice)
  double diffc10g2f = pdf_factor_gx2*g_z1*Cqg(z1) * (-g_z2*z2*Pqq(z2) + Pqqint(x2));
  diffc10g2f += pdf_factor_gz2*g_z1*g_z2*Cqg(z1)*Pqq(z2);

  // (checked twice)
  double diffc1fg20 = pdf_factor_z1g*g_z1*g_z2*Cqq(z1)*Pqg(z2);
  diffc1fg20 += pdf_factor_x1g*(-g_z2*C1qqdelta(A_F)*Pqg(z2));

  tcga = tcga + diffc10g20 + diffc1fg2f + diffc10g2f + diffc1fg20;

  //  cout << setw(19) << "contribution: " << "tcga 22 =   " << diffc1fg2f << endl;
  //  cout << setw(19) << "contribution: " << "tcga 24 =   " << diffc1fg20 << endl;
  //  cout << setw(19) << "contribution: " << "tcga 42 =   " << diffc10g2f << endl;
  //  cout << setw(19) << "contribution: " << "tcga 44 =   " << diffc10g20 << endl;

  //  cout << setw(19) << "contribution: " << "tcga 22 sum =   " << diffc1fg2f + diffg1fc2f << endl;
  //  cout << setw(19) << "contribution: " << "tcga 24 sum =   " << diffc1fg20 + diffg1fc20 << endl;
  //  cout << setw(19) << "contribution: " << "tcga 42 sum =   " << diffc10g2f + diffg10c2f << endl;
  //  cout << setw(19) << "contribution: " << "tcga 44 sum =   " << diffc10g20 + diffg10c20 << endl;
  /*
  cout << setw(19) << "contribution: " << "tcga 22 sum1 =   " << pdf_factor_z1z2*(g_z1*g_z2*Cqq(z2)*Pqq(z1)) + g_z1*Cqq(z1)*g_z2*Pqq(z2)*(pdf_factor_z1z2) << endl;
  cout << setw(19) << "contribution: " << "tcga 22 sum2 =   " << pdf_factor_z1x2*(-g_z1*Pqq(z1)*C1qqdelta(A_F)) + g_z1*Cqq(z1)*g_z2*Pqq(z2)*(-pdf_factor_z1x2*z2) + Pqqint(x2)*(pdf_factor_z1x2*g_z1*Cqq(z1)) << endl;
  cout << setw(19) << "contribution: " << "tcga 22 sum3 =   " << pdf_factor_x1z2*g_z2*Cqq(z2)*(-g_z1*Pqq(z1)*z1+Pqqint(x1)) + C1qqdelta(A_F)*g_z2*Pqq(z2)*(-pdf_factor_x1z2) << endl;
  cout << setw(19) << "contribution: " << "tcga 22 sum4 =   " << pdf_factor_x1x2*C1qqdelta(A_F)*(g_z1*z1*Pqq(z1)-Pqqint(x1)) + C1qqdelta(A_F)*g_z2*Pqq(z2)*(pdf_factor_x1x2*z2) + Pqqint(x2)*(-pdf_factor_x1x2*C1qqdelta(A_F)) << endl;
  
  cout << setw(19) << "contribution: " << "tcga 24 sum1 =   " << pdf_factor_z1g*g_z1*g_z2*Cqq(z1)*Pqg(z2) + pdf_factor_z1g*g_z1*g_z2*Pqq(z1)*Cqg(z2) << endl;
  cout << setw(19) << "contribution: " << "tcga 24 sum0 =   " << pdf_factor_x1g*(-g_z2*C1qqdelta(A_F)*Pqg(z2)) + pdf_factor_x1g*g_z2*Cqg(z2)*(-g_z1*Pqq(z1)*z1+Pqqint(x1)) << endl;

  cout << setw(19) << "contribution: " << "tcga 42 sum1 =   " << pdf_factor_gx2*g_z1*Cqg(z1)*(-g_z2*z2*Pqq(z2)+Pqqint(x2)) + pdf_factor_gx2*(-g_z1*C1qqdelta(A_F)*Pqg(z1)) << endl;
  cout << setw(19) << "contribution: " << "tcga 42 sum0 =   " << pdf_factor_gz2*g_z1*g_z2*Cqg(z1)*Pqq(z2) + pdf_factor_gz2*g_z1*g_z2*Pqg(z1)*Cqq(z2) << endl;
  */
//c    C first leg, gamma second leg

//      diffc1f=-dlog(xx10)*fx1p(j)*Cqq(z1)+C1qqdelta*fx10(j)

//      diffc10=-dlog(xx10)*fx1p(0)*Cqg(z1)

//      tcga=tcga+msqc(j,k)*
//     # (flgq*diffc10*diffg20+flgq*diffc1f*diffg2f
//     #          +diffc10*diffg2f+diffc1f*diffg20)

  // (checked twice)
  tcga = tcga - g_z1 * (pdf_factor_z1x2*CqqPqq(z1)+pdf_factor_gx2*(CqqPqg(z1)+CqgPgg(z1,beta0)));
  tcga = tcga - g_z2 * (pdf_factor_x1z2*CqqPqq(z2)+pdf_factor_x1g*(CqqPqg(z2)+CqgPgg(z2,beta0)));

  //  cout << setw(19) << "contribution: " << "tcga qg =   " << -g_z1*(pdf_factor_gx2*(CqqPqg(z1)+CqgPgg(z1,beta0))) << endl;
  //  cout << setw(19) << "contribution: " << "tcga qq =   " << -g_z1*(pdf_factor_z1x2*CqqPqq(z1)) << endl;
  //  cout << setw(19) << "contribution: " << "tcga qg =   " << -g_z2*(pdf_factor_x1g*(CqqPqg(z2)+CqgPgg(z2,beta0))) << endl;
  //  cout << setw(19) << "contribution: " << "tcga qq =   " << -g_z2*(pdf_factor_x1z2*CqqPqq(z2)) << endl;

//c    C*gamma: first leg (ignore delta term in Cqq: taken into account with tH1stF)

//      tcga = tcga
//     &     +(fx1p(j)*CqqPqq(z1)*flgq+fx1p(0)*(CqqPqg(z1)+CqgPgg(z1)))
//     &     *(-dlog(xx10))*fx20(k)*msqc(j,k)

//c    C*gamma: second leg (ignore delta term in Cqq: taken into account with tH1stF)

//      tcga = tcga
//     &     +(fx2p(k)*CqqPqq(z2)*flgq+fx2p(0)*(CqqPqg(z2)+CqgPgg(z2)))
//     &     *(-dlog(xx20))*fx10(j)*msqc(j,k)

  // (checked twice)
  tcga = tcga - g_z1*pdf_factor_qx2*CqgPgq(z1);
  tcga = tcga - g_z2*pdf_factor_x1q*CqgPgq(z2);
  //  cout << setw(19) << "contribution: " << "tcga qQ =   " << -g_z1*pdf_factor_qx2*CqgPgq(z1) << endl;
  //  cout << setw(19) << "contribution: " << "tcga qQ =   " << -g_z2*pdf_factor_x1q*CqgPgq(z2) << endl;

//c    Add Cqg*Pgq contribution

//      do l=1,nf
//      tcga=tcga+(fx1p(l)+fx1p(-l))*CqgPgq(z1)
//     &           *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
//      tcga=tcga+(fx2p(l)+fx2p(-l))*CqgPgq(z2)
//     &           *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
//      enddo

//CC  Start 2-loop AP

//-------------------------------------------------------------------------------------------------------------

//C   Gluon + pure singlet
  // (checked twice)
  tgamma2 = -g_z1*P2qg(z1)*pdf_factor_gx2;
  tgamma2 = tgamma2 - g_z2*P2qg( z2)*pdf_factor_x1g;
  tgamma2 = tgamma2 - g_z1*P2qqS(z1)*pdf_factor_qx2;
  tgamma2 = tgamma2 - g_z2*P2qqS(z2)*pdf_factor_x1q;

//      do l=-nf,nf
//      if(l.eq.0) then
//      tgamma2=tgamma2+fx1p(0)*P2qg(z1)
//     & *(-dlog(xx10))*fx20(k)*msqc(j,k)
//      tgamma2=tgamma2+fx2p(0)*P2qg(z2)
//     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
//      else
//      tgamma2=tgamma2+fx1p(l)*P2qqS(z1)
//     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
//      tgamma2=tgamma2+fx2p(l)*P2qqS(z2)
//     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
//      endif
//      enddo

//C   P2qq non-singlet: regular part
  // (checked twice)
  tgamma2 = tgamma2 - g_z1*P2qqV(z1,nf)*pdf_factor_z1x2;
  tgamma2 = tgamma2 - g_z2*P2qqV(z2,nf)*pdf_factor_x1z2;

//      tgamma2=tgamma2+fx1p(j)*P2qqV(z1)
//     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
//      tgamma2=tgamma2+fx2p(k)*P2qqV(z2)
//     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq

//C   P2qq non-singlet: 1/(1-z)_+
// (checked twice)
//  tgamma2=tgamma2+2.0/3*Kappa*(pdf_factor_x1x2*(g_z1*z1/(1-z1)-D0int(x1)) - g_z1*pdf_factor_z1x2/(1-z1));
//  tgamma2=tgamma2+2.0/3*Kappa*(pdf_factor_x1x2*(g_z2*z2/(1-z2)-D0int(x2)) - g_z2*pdf_factor_x1z2/(1-z2));

  tgamma2 = tgamma2 + 2.0/3*Kappa* ( g_z1/(1-z1)*(pdf_factor_x1x2*z1-pdf_factor_z1x2) - D0int(x1)*pdf_factor_x1x2 );
  tgamma2 = tgamma2 + 2.0/3*Kappa* ( g_z2/(1-z2)*(pdf_factor_x1x2*z2-pdf_factor_x1z2) - D0int(x2)*pdf_factor_x1x2 );

//      diff=-dlog(xx10)
//     &  *(fx1p(j)-fx10(j)*z1)/(1-z1)
//     &  - D0int(xx10)*fx10(j)
//
//      tgamma2=tgamma2+2d0/3*Kappa*diff*fx20(k)*msqc(j,k)*flgq

//      diff=-dlog(xx20)
//     &  *(fx2p(k)-fx20(k)*z2)/(1-z2)
//     &  - D0int(xx20)*fx20(k)

//      tgamma2=tgamma2+2d0/3*Kappa*diff*fx10(j)*msqc(j,k)*flgq

//C   P2qqb non singlet
  // (checked twice)
  tgamma2 = tgamma2 - g_z1*P2qqbV(z1)*pdf_factor_qbx2;
  tgamma2 = tgamma2 - g_z2*P2qqbV(z2)*pdf_factor_x1qb;

//      tgamma2=tgamma2+fx1p(-j)*P2qqbV(z1)
//     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

//      tgamma2=tgamma2+fx2p(-k)*P2qqbV(z2)
//     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
}

//#include "old.counterterm.cpp"


//-------------------------------------------------------------------------------------------------------------------------------------------------                                                                         


double observable_set::calculate_sigma_qqbar(double qt2, double q2, double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb, int order, double A_F, double LR, double LF, double LQ, vector<double> & sigma_qT, phasespace_set & psi){
  Logger logger("observable_set::sigma_qqbar");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
  static int init = 1;
  static vector<double> I1_int, I2_int, I3_int, I4_int;
  if (init == 1) 
    {
      I1_int.resize(n_qTcut);
      I2_int.resize(n_qTcut);
      I3_int.resize(n_qTcut);
      I4_int.resize(n_qTcut);
      computeItildaIntegrals(min_qTcut/100, step_qTcut/100, n_qTcut, I1_int, I2_int, I3_int, I4_int);      
      init = 0;
    }

  double qtoverq = sqrt(qt2/q2);

  logger << LOG_DEBUG_VERBOSE << "before Itilde qtoverq = " << qtoverq << endl;
  logger << LOG_DEBUG_VERBOSE << "before Itilde qt2 = " << qt2 << endl;
  logger << LOG_DEBUG_VERBOSE << "before Itilde q2 = " << q2 << endl;
  logger << LOG_DEBUG_VERBOSE << "before Itilde order = " << order << endl;

  double LL1, LL2, LL3, LL4;

  // Itilde(qtoverq,order,LL1,LL2,LL3,LL4);
  // in FO computations, we do the qT integrals once and for all in the beginning
  // in resummed computations, we have to bin the "artificial" qT, so that is not easily possible -> revert back to old implementation  
  if (switch_resummation){Itilde(qtoverq, order, LL1, LL2, LL3, LL4);}

  logger << LOG_DEBUG_VERBOSE << "after Itilde LL1 = " << LL1 << endl;
  logger << LOG_DEBUG_VERBOSE << "after Itilde LL2 = " << LL2 << endl;
  logger << LOG_DEBUG_VERBOSE << "after Itilde LL3 = " << LL3 << endl;
  logger << LOG_DEBUG_VERBOSE << "after Itilde LL4 = " << LL4 << endl;

  double sig = 0.;

  double tH1F = H1F(z1, z2, g_z1, g_z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_gx2, pdf_factor_x1g);
  logger << LOG_DEBUG_VERBOSE << "after tH1F " << endl;

  double sig12 = sigma12_qqbar(pdf_factor_x1x2);
  double sig11 = sigma11_qqbar(pdf_factor_x1x2, tH1F, LQ);
  logger << LOG_DEBUG_VERBOSE << "after sig11 " << endl;

  if (order == 1){
    sig = alpha_S/M_PI*(sig12*LL2 + sig11*LL1);
      
    for (int i=0; i<n_qTcut; i++){
      sigma_qT[i] = alpha_S/M_PI*(sig12*I2_int[i] + sig11*I1_int[i]);
    }

    logger << LOG_DEBUG_VERBOSE << "finished sig/q2 = " << sig/q2 << endl;
    return sig/q2;
  }
  
  double tgaga = 0.0, tcga = 0.0, tgamma2 = 0.0;
  
  const double CF = 4.0/3;
  double H1_delta = M_PI*M_PI/6*CF + 0.5*A_F;

  
  logger << LOG_DEBUG_VERBOSE << "QT:   z1       = " << z1       << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   z2       = " << z2       << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   g_z1     = " << g_z1     << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   g_z2     = " << g_z2     << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   A_F      = " << A_F      << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   LR       = " << LR       << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   LF       = " << LF       << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   LQ       = " << LQ       << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   H1_delta = " << H1_delta << endl; 

  double tH1    = H1_M_qqbar(z1, z2, g_z1, g_z2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_gx2, pdf_factor_x1g, H1_delta);
  double H1full = H1_qqbar(pdf_factor_x1x2, tH1, tH1F, LR, LF, LQ);

  double beta0 = (33.0-2.0*N_f)/12;

  calcIntermediateTerms(z1, z2, g_z1, g_z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_gx2, pdf_factor_x1g, pdf_factor_z1z2, pdf_factor_gg, pdf_factor_qx2, pdf_factor_x1q, pdf_factor_gz2, pdf_factor_z1g, pdf_factor_qbx2, pdf_factor_x1qb, N_f, A_F, beta0, tgaga, tcga, tgamma2);

//  cout << "H1=" << tH1 << ", H1F=" << tH1F << ", tgaga=" << tgaga << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   LR       = " << LR       << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   LF       = " << LF       << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   A_F      = " << A_F      << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   sig11   = " << sig11 << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tH1     = " << tH1 << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tH1F    = " << tH1F << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   H1full  = " << H1full << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tgaga   = " << tgaga << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tcga    = " << tcga << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tgamma2 = " << tgamma2 << endl;

  double sig24 = sigma24_qqbar(pdf_factor_x1x2);
  double sig23 = sigma23_qqbar(pdf_factor_x1x2, sig11, beta0);
  double sig22 = sigma22_qqbar(pdf_factor_x1x2, beta0, sig11, tH1, tH1F, H1full, tgaga, LR, LF, LQ, N_f);
  double sig21 = sigma21_qqbar(pdf_factor_x1x2, beta0, sig11, tH1, tH1F, H1full, tgaga, tcga, tgamma2, LR, LF, LQ, A_F, N_f);

  /*
  logger << LOG_DEBUG << "QT:   sig21   = " << sig21 << endl;
  logger << LOG_DEBUG << "QT:   sig22   = " << sig22 << endl;
  logger << LOG_DEBUG << "QT:   sig23   = " << sig23 << endl;
  logger << LOG_DEBUG << "QT:   sig24   = " << sig24 << endl;
  */

  sig = pow(alpha_S/M_PI,2)*(sig24*LL4+sig23*LL3+sig22*LL2+sig21*LL1);

  for (int i=0; i<n_qTcut; i++){
    sigma_qT[i] = pow(alpha_S/M_PI,2)*(sig24*I4_int[i]+sig23*I3_int[i]+sig22*I2_int[i]+sig21*I1_int[i]);
    logger << LOG_DEBUG_VERBOSE << "QT:   sigma_qT[" << setw(3) << i << "] = " << sigma_qT[i] << endl;
  }
  //  cout << "sigma_qT[0] = " << sigma_qT[0] << endl;
  //  sig=pow(alpha_S/M_PI,2)*LL3;

  //  cout << sig << endl;
  //  sig2=sig24*LL4+sig23*LL3+sig22*LL2+sig21*LL1

  //  cout << sig12 << ", " << sig11 << ", " << sig24 << ", " << sig23 << ", " << sig22 << ", " << sig21 << endl;
  //  cout << LL1 << ", " << LL2 << ", " << LL3 << ", " << LL4 << endl << endl;

  logger << LOG_DEBUG_VERBOSE << "finished   sig/q2 = " << sig/q2 << endl;
  return sig/q2;
}

double H1_qqbar(double pdf_factor_x1x2, double tH1, double tH1F, double LR, double LF, double LQ) {
  double H = tH1 + (LF - LQ) * tH1F - LQ * pdf_factor_x1x2 * (B1q + 0.5 * A1q * LQ);

  return H;
}

double H1_M_qqbar(double z1, double z2, double g_z1, double g_z2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double H1_delta) {
//  C     H1st delta term
//double tH1st = 2*C1qqdelta(A_F)*pdf_factor_x1x2;
  double tH1st = H1_delta        *pdf_factor_x1x2;

  /*
  cout << "H1st delta term     = " << H1_delta*pdf_factor_x1x2 << endl;
  cout << "H1st delta term / 2 = " << H1_delta*pdf_factor_x1x2 / 2 << endl;
  */

//     H1st: non delta terms, first leg
  tH1st = tH1st - g_z1* ( Cqq(z1)*pdf_factor_z1x2 + Cqg(z1)*pdf_factor_gx2 );
  /*
  cout << "H1st  Cqq(z1)  term     = " << Cqq(z1) << endl;
  cout << "H1st  Cqg(z1)  term     = " << Cqg(z1) << endl;
  cout << "H1st  Cqq(z1)*(-g_z1)  term     = " << Cqq(z1)*(-g_z1) << endl;
  cout << "H1st  Cqg(z1)*(-g_z1)  term     = " << Cqg(z1)*(-g_z1) << endl;
  cout << "H1st  pdf  qq1  term     = " << (pdf_factor_z1x2) << endl;
  cout << "H1st  pdf  qg1  term     = " << (pdf_factor_gx2) << endl;
  cout << "H1st  qq1  term     = " << (Cqq(z1)*pdf_factor_z1x2)*(-g_z1) << endl;
  cout << "H1st  qg1  term     = " << (Cqg(z1)*pdf_factor_gx2)*(-g_z1) << endl;
  */  
//      tH1st=tH1st+(fx1p(j)*Cqq(z1)*flgq+fx1p(0)*Cqg(z1))*(-dlog(xx10))*fx20(k)*msqc(j,k)

//     H1st: non delta terms, second leg
  tH1st = tH1st - g_z2* ( Cqq(z2)*pdf_factor_x1z2 + Cqg(z2)*pdf_factor_x1g );
  /*
  cout << "H1st  Cqq(z2)  term     = " << Cqq(z2) << endl;
  cout << "H1st  Cqg(z2)  term     = " << Cqg(z2) << endl;
  cout << "H1st  Cqq(z2)*(-g_z2)  term     = " << Cqq(z2)*(-g_z2) << endl;
  cout << "H1st  Cqg(z2)*(-g_z2)  term     = " << Cqg(z2)*(-g_z2) << endl;
  cout << "H1st  pdf  qq2  term     = " << (pdf_factor_x1z2) << endl;
  cout << "H1st  pdf  qg2  term     = " << (pdf_factor_x1g) << endl;
  cout << "H1st  qq2  term     = " << (Cqq(z2)*pdf_factor_x1z2)*(-g_z2) << endl;
  cout << "H1st  qg2  term     = " << (Cqg(z2)*pdf_factor_x1g)*(-g_z2) << endl;
  */
//      tH1st=tH1st+(fx2p(k)*Cqq(z2)*flgq+fx2p(0)*Cqg(z2))*(-dlog(xx20))*fx10(j)*msqc(j,k)

  return tH1st;
}

double H1F(double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g) {
//C     H1st: muf dependence
//c     gammaqq and gammaqg: first leg
  double tH1stF =  - g_z1      *( (pdf_factor_z1x2 - pdf_factor_x1x2*z1)*Pqq(z1)      + Pqg(z1)*pdf_factor_gx2)            - pdf_factor_x1x2*Pqqint(x1);
//       tH1stF =  - dlog(xx10)*( (fx1p(j)         - fx10(j)*xx10**beta)*Pqq(z1)*flgq + fx1p(0)*Pqg(z1))*fx20(k)*msqc(j,k) - Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)*flgq
                       
//c     gammaqq and gammaqg: second leg
   tH1stF = tH1stF - g_z2      *( (pdf_factor_x1z2 - pdf_factor_x1x2*z2)*Pqq(z2)      + Pqg(z2)*pdf_factor_x1g)            - pdf_factor_x1x2*Pqqint(x2);
// tH1stF = tH1stF - dlog(xx20)*( (fx2p(k)         - fx20(k)*xx20**alfa)*Pqq(z2)*flgq + fx2p(0)*Pqg(z2))*fx10(j)*msqc(j,k) - Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)*flgq

  return tH1stF;
}



double H2_M_qqbar(double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb, double nf, double H2_delta) {
  const double H2qqD0 = -404.0/27+(56*nf)/81+14*zeta3;

  //  CC    H2st gg contribution
// checked once
  double tH2 =          g_z1*g_z2*Cqg(z1)*Cqg(z2)*pdf_factor_gg;
//  cout << setw(50) << "H2st gg contribution" << setw(20) << "pdf_factor_gg = " << g_z1*g_z2*Cqg(z1)*Cqg(z2)*pdf_factor_gg << endl;
//      tH2st = tH2st + fx1p(0)*Cqg(z1)*(-dlog(xx10)) * fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)*flgq

//CC    H2st qqbar contribution from C1*C1 (without delta term)

//C     regular*regular
// checked once
  tH2=tH2+g_z1*g_z2*Cqq(z1)*Cqq(z2)*pdf_factor_z1z2;
//  cout << setw(50) << "H2st qqbar contribution from C1*C1 (without delta term)" << setw(20) << "pdf_factor_z1z2 = " << +g_z1*g_z2*Cqq(z1)*Cqq(z2)*pdf_factor_z1z2 << endl;
//      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10)) * fx2p(k)*Cqq(z2)*(-dlog(xx20))*msqc(j,k)*flgq

//C     regular-delta
// checked once
//  tH2=tH2-g_z1*C1qqdelta(A_F)*Cqq(z1)*pdf_factor_z1x2;
//  tH2=tH2-g_z2*C1qqdelta(A_F)*Cqq(z2)*pdf_factor_x1z2;
//  tH2=tH2-0.5*g_z1*H1_delta*Cqq(z1)*pdf_factor_z1x2;
//  tH2=tH2-0.5*g_z2*H1_delta*Cqq(z2)*pdf_factor_x1z2;

//      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
//     &            fx20(k)*C1qqdelta*msqc(j,k)*flgq
//
//      tH2st=tH2st+fx2p(k)*Cqq(z2)*(-dlog(xx20))*
//     &            fx10(j)*C1qqdelta*msqc(j,k)*flgq


//CC    H2st qg contribution from C1*C1

//C     regular*regular
// checked once
  tH2=tH2+g_z1*g_z2*Cqg(z1)*Cqq(z2)*pdf_factor_gz2;
//  cout << setw(50) << "H2st qg contribution from C1*C1" << setw(20) << "pdf_factor_gz2 = " << +g_z1*g_z2*Cqg(z1)*Cqq(z2)*pdf_factor_gz2 << endl;
  tH2=tH2+g_z1*g_z2*Cqg(z2)*Cqq(z1)*pdf_factor_z1g;
//  cout << setw(50) << "H2st qg contribution from C1*C1" << setw(20) << "pdf_factor_z1g = " << +g_z1*g_z2*Cqg(z2)*Cqq(z1)*pdf_factor_z1g << endl;

//      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10)) * fx2p(k)*Cqq(z2)*(-dlog(xx20))*msqc(j,k)
//
//      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10)) * fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)


//C     regular-delta
// checked once
//  tH2=tH2-g_z1*Cqg(z1)*C1qqdelta(A_F)*pdf_factor_gx2;
//  tH2=tH2-g_z2*Cqg(z2)*C1qqdelta(A_F)*pdf_factor_x1g;
//  tH2=tH2-0.5*g_z1*Cqg(z1)*H1_delta*pdf_factor_gx2;
//  tH2=tH2-0.5*g_z2*Cqg(z2)*H1_delta*pdf_factor_x1g;

//      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
//     &            fx20(k)*C1qqdelta*msqc(j,k)
//
//      tH2st=tH2st+fx2p(0)*Cqg(z2)*(-dlog(xx20))*
//     &            fx10(j)*C1qqdelta*msqc(j,k)

//CC    H2st qqbar channel: D0(z), first leg
// checked once
  tH2=tH2+0.5*g_z1*H2qqD0/(1-z1)*(-pdf_factor_z1x2+z1*pdf_factor_x1x2);
  //  cout << setw(50) << "H2st qqbar channel: D0(z), first leg" << setw(20) << " = " << +0.5*g_z1*H2qqD0/(1-z1)*(-pdf_factor_z1x2+z1*pdf_factor_x1x2) << endl;
//  cout << setw(50) << "H2st qqbar channel: D0(z), first leg" << setw(20) << "pdf_factor_x1x2 = " << +0.5*g_z1*H2qqD0/(1-z1)*(+z1*pdf_factor_x1x2) << "   " << pdf_factor_x1x2 << endl;
//  cout << setw(50) << "H2st qqbar channel: D0(z), first leg" << setw(20) << "pdf_factor_z1x2 = " << +0.5*g_z1*H2qqD0/(1-z1)*(-pdf_factor_z1x2) << "   " << pdf_factor_z1x2 << endl;
  tH2=tH2-0.5*pdf_factor_x1x2*H2qqD0*D0int(x1);
//  cout << setw(50) << "H2st qqbar channel: D0(z), first leg" << setw(20) << "pdf_factor_x1x2 = " << -0.5*pdf_factor_x1x2*H2qqD0*D0int(x1) << "   " << pdf_factor_x1x2 << endl;

//      diff=-dlog(xx10)*(fx1p(j)-fx10(j)*xx10**beta)*H2qqD0/(1-z1)
//
//      tH2st=tH2st+0.5d0*diff*fx20(k)*msqc(j,k)*flgq
//      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx10)
//     &            *fx10(j)*fx20(k)*msqc(j,k)*flgq

//CC    H2st, qqbar channel: D0(z), second leg
// checked once
  tH2=tH2+0.5*g_z2*H2qqD0/(1-z2)*(-pdf_factor_x1z2+z2*pdf_factor_x1x2);
  //  cout << setw(50) << "H2st, qqbar channel: D0(z), second leg" << setw(20) << " = " << +0.5*g_z2*H2qqD0/(1-z2)*(-pdf_factor_x1z2+z2*pdf_factor_x1x2) << endl;
//  cout << setw(50) << "H2st, qqbar channel: D0(z), second leg" << setw(20) << "pdf_factor_x1x2 = " << +0.5*g_z2*H2qqD0/(1-z2)*(z2*pdf_factor_x1x2) << "   " << pdf_factor_x1x2 << endl;
//  cout << setw(50) << "H2st, qqbar channel: D0(z), second leg" << setw(20) << "pdf_factor_x1z2 = " << +0.5*g_z2*H2qqD0/(1-z2)*(-pdf_factor_x1z2) << "   " << pdf_factor_x1z2 << endl;
  tH2=tH2-0.5*pdf_factor_x1x2*H2qqD0*D0int(x2);
//  cout << setw(50) << "H2st, qqbar channel: D0(z), second leg" << setw(20) << "pdf_factor_x1x2 = " << -0.5*pdf_factor_x1x2*H2qqD0*D0int(x2) << "   " << pdf_factor_x1x2 << endl;

//      diff=-dlog(xx20)*(fx2p(k)-fx20(k)*xx20**alfa)*H2qqD0/(1-z2)
//
//      tH2st=tH2st+0.5d0*diff*fx10(j)*msqc(j,k)*flgq
//      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx20)
//     &            *fx10(j)*fx20(k)*msqc(j,k)*flgq

//CC    C2qq, regular part, first leg
// checked once
  tH2=tH2-g_z1*C2qqreg(z1,nf)*pdf_factor_z1x2;
//  cout << setw(50) << "C2qq, regular part, first leg = " << setw(20) << "pdf_factor_z1x2 = " << -g_z1*C2qqreg(z1,nf)*pdf_factor_z1x2 << endl;

//      tH2st=tH2st+fx1p(j)*C2qqreg(z1)
//     &                   *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

//CC    C2qq, regular part, second leg
// checked once
  tH2=tH2-g_z2*C2qqreg(z2,nf)*pdf_factor_x1z2;
//  cout << setw(50) << "C2qq, regular part, second leg" << setw(20) << "pdf_factor_x1z2 = " << -g_z2*C2qqreg(z2,nf)*pdf_factor_x1z2 << endl;

//      tH2st=tH2st+fx2p(k)*C2qqreg(z2)
//     &                   *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq

//CC    C2qg, first leg
// checked once
  tH2=tH2-g_z1*C2qg(z1)*pdf_factor_gx2;
//  cout << setw(50) << "C2qg, first leg" << setw(20) << "pdf_factor_gx2 = " << -g_z1*C2qg(z1)*pdf_factor_gx2 << endl;

//      tH2st=tH2st+fx1p(0)*C2qg(z1)*(-dlog(xx10))*fx20(k)*msqc(j,k)

//CC    C2qg, second leg
// checked once
  tH2=tH2-g_z2*C2qg(z2)*pdf_factor_x1g;
//  cout << setw(50) << "C2qg, second leg" << setw(20) << "pdf_factor_x1g = " << -g_z2*C2qg(z2)*pdf_factor_x1g << endl;

//      tH2st=tH2st+fx2p(0)*C2qg(z2)*(-dlog(xx20))*fx10(j)*msqc(j,k)

//CC    Cqqbar contribution: first leg
// checked once
  tH2=tH2-g_z1*C2qqb(z1,nf)*pdf_factor_qbx2;
//  cout << setw(50) << "Cqqbar contribution: first leg" << setw(20) << "pdf_factor_qbx2 = " << -g_z1*C2qqb(z1,nf)*pdf_factor_qbx2 << endl;

//      tH2st=tH2st+fx1p(-j)*C2qqb(z1)
//     &                    *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

//CC    Cqqbar contribution: second leg
// checked once
  tH2=tH2-g_z2*C2qqb(z2,nf)*pdf_factor_x1qb;
//  cout << setw(50) << "Cqqbar contribution: second leg" << setw(20) << "pdf_factor_x1qb = " << -g_z2*C2qqb(z2,nf)*pdf_factor_x1qb << endl;

//      tH2st=tH2st+fx2p(-k)*C2qqb(z2) *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq

//      do a=1,nf

//CC    Cqqp contribution: first leg
  // checked once
  tH2=tH2-g_z1*C2qqp(z1,nf)*pdf_factor_qx2;
  tH2=tH2+g_z1*C2qqp(z1,nf)*(pdf_factor_qbx2+pdf_factor_z1x2);

//      if(a.ne.abs(j)) then
//       tH2st=tH2st+(fx1p(a)+fx1p(-a))*
//     &       C2qqp(z1)*(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
//      endif

//CC    Cqqp contribution: second leg
// checked once
  tH2=tH2-g_z2*C2qqp(z2,nf)*pdf_factor_x1q;
  tH2=tH2+g_z2*C2qqp(z2,nf)*(pdf_factor_x1qb+pdf_factor_x1z2);

//      if(a.ne.abs(k)) then
//       tH2st=tH2st+(fx2p(a)+fx2p(-a))*
//     &       C2qqp(z2)*(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
//      endif

//      enddo

  // add delta contribution
  tH2=tH2+H2_delta*pdf_factor_x1x2;

  return tH2;
}

double H2_qqbar(double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb,double tH1, double tH1F, double nf, double H1_delta, double H2_delta, double LR, double LF, double LQ) {
  Logger logger("H2_qqbar");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  double H=0;

  double sig11=sigma11_qqbar(pdf_factor_x1x2, tH1F, LQ);

  double beta0=(33.0-2.0*nf)/12;
  const double CF=4.0/3;

  double tgaga=0.0,tcga=0.0,tgamma2=0.0;
  double A_F=2*(H1_delta-M_PI*M_PI/6*CF);
  calcIntermediateTerms(z1,z2,g_z1,g_z2,x1,x2,pdf_factor_x1x2,pdf_factor_z1x2,pdf_factor_x1z2,pdf_factor_gx2,pdf_factor_x1g,pdf_factor_z1z2,pdf_factor_gg, pdf_factor_qx2,pdf_factor_x1q,pdf_factor_gz2,pdf_factor_z1g,pdf_factor_qbx2,pdf_factor_x1qb,nf,A_F,beta0,tgaga,tcga,tgamma2);

  double tH2=H2_M_qqbar(z1, z2, g_z1, g_z2, x1, x2, pdf_factor_x1x2,pdf_factor_z1x2,pdf_factor_x1z2,pdf_factor_gx2,pdf_factor_x1g,pdf_factor_gg, pdf_factor_qx2,pdf_factor_x1q,pdf_factor_z1z2,pdf_factor_gz2,pdf_factor_z1g,pdf_factor_qbx2,pdf_factor_x1qb,nf, H2_delta);

  logger << LOG_DEBUG_VERBOSE << "QT:   tH2      = " << tH2 << endl;
  tH2 += H1_delta*(tH1-H1_delta*pdf_factor_x1x2);
  logger << LOG_DEBUG_VERBOSE << "QT:   tH2 (incl. H1²) = " << tH2 << endl;

  logger << LOG_DEBUG_VERBOSE << "QT:   LR              = " << LR       << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   LF              = " << LF       << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   A_F             = " << A_F      << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   sig11           = " << sig11 << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tH1             = " << tH1 << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tH1F            = " << tH1F << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tgaga           = " << tgaga << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tcga            = " << tcga << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tgamma2         = " << tgamma2 << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tH2             = " << tH2 << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   H1_delta        = " << H1_delta << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   H2_delta        = " << H2_delta << endl;


  H=tH2;

//  xmsq=xmsq+asopi*(tH1st+LF*tH1stF) + asopi**2*(tdelta*H2qqdelta+tH2st)

//CC     add scale dependence at NNLO
  double H1full=H1_qqbar(pdf_factor_x1x2,tH1,tH1F,LR,LF,LQ);
  logger << LOG_DEBUG_VERBOSE << "QT:   H1full   = " << H1full << endl;

  H=H+0.5*beta0*LF*LF*tH1F+tgamma2*LF-beta0*LR*H1full+(LF-LQ)*tcga+0.5*(LF-LQ)*(LF-LQ)*tgaga;

  // missing delta terms. Note that tcga does contain the delta term, so we only need to add it once
  H=H+(LF-LQ)*C1qqdelta(A_F)*tH1F;

  const double Delta2qq=(16.0/9*(3.0/8-M_PI*M_PI/2+6*zeta3)+4*(17.0/24+11.0*M_PI*M_PI/18-3*zeta3)-2.0/3*nf*(1.0/6+2*M_PI*M_PI/9))/4;
  H=H+LF*2*Delta2qq*pdf_factor_x1x2;

  // add Q dependence
  double Kappa=67.0/6-(M_PI*M_PI)/2-5.0/9*nf;
  double A2q=0.5*A1q*Kappa;
  H=H+(1.0/6*A1q*beta0*LQ*LQ*LQ*pdf_factor_x1x2+0.5*LQ*LQ*(A2q*pdf_factor_x1x2+beta0*sig11));

  const double B2q=-2*Delta2qq+0.5*beta0*(2.0/3*4.0/3*M_PI*M_PI+A_F);
  H=H-LQ*((B2q+LQ*A2q)*pdf_factor_x1x2-beta0*tH1+tgamma2+2*Delta2qq*pdf_factor_x1x2);

  H=H-0.5*LQ*(B1q+0.5*A1q*LQ)*(H1full+tH1);

  H=H+0.5*LQ*LQ*(B1q+0.5*A1q*LQ)*(B1q+0.5*A1q*LQ)*pdf_factor_x1x2;

  H=H-0.5*LQ*(LF-LQ)*(B1q+0.5*A1q*LQ)*tH1F;
  
  return H;
}

//--------------------------------------------------------------------------------------------------------------------------------------------

double observable_set::virtual_qqbar(double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb, int order, double H1_delta, double H2_delta, double LR, double LF, double LQ) {
  Logger logger("observable_set::virtual_qqbar");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  //  cout << "virtual_qqbar: " << order << ", " << LQ << endl;

  if (order == 0)
    return pdf_factor_x1x2;

  double tH1  = H1_M_qqbar(z1, z2, g_z1, g_z2,         pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_gx2, pdf_factor_x1g, H1_delta);
  double tH1F = H1F (z1, z2, g_z1, g_z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_gx2, pdf_factor_x1g);

  if (order == 1)     
    return alpha_S/M_PI * H1_qqbar( pdf_factor_x1x2, tH1, tH1F, LR, LF, LQ );
       
  double nf=(double)(N_f);

  logger << LOG_DEBUG_VERBOSE << "QT:   LR       = " << LR       << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   LF       = " << LF       << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tH1      = " << tH1 << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   tH1F     = " << tH1F << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   H1_delta = " << H1_delta << endl;
  logger << LOG_DEBUG_VERBOSE << "QT:   H2_delta = " << H1_delta << endl;
  
  return pow(alpha_S/M_PI,2)*H2_qqbar(z1, z2, g_z1, g_z2, x1, x2, pdf_factor_x1x2, pdf_factor_z1x2, pdf_factor_x1z2, pdf_factor_gx2, pdf_factor_x1g, pdf_factor_gg, pdf_factor_qx2, pdf_factor_x1q, pdf_factor_z1z2, pdf_factor_gz2, pdf_factor_z1g, pdf_factor_qbx2, pdf_factor_x1qb, tH1, tH1F, nf, H1_delta, H2_delta, LR, LF, LQ);
}
