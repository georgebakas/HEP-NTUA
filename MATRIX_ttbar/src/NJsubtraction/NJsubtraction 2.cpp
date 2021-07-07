#include "NJsubtraction.h"
#include <cmath>
#include <assert.h>

// {{{ NJsubtraction::NJsubtraction ()
//NJsubtraction::NJsubtraction () : Gamma_0(2),Gamma_1(2),gammaB0(2),gammaB1(2),gammaJ0(2),gammaJ1(2),Ccoef(2),T2(2),TiTj(2,vector<vector<double> >(2,vector<double>(3))) {
NJsubtraction::NJsubtraction(){
  Gamma_0.resize(2);
  Gamma_1.resize(2);
  gammaB0.resize(2);
  gammaB1.resize(2);
  gammaJ0.resize(2);
  gammaJ1.resize(2);
  ///  Ccoef.resize(2);
  Ccoef.resize(3);
  T2.resize(2);
  TiTj.resize(2, vector<vector<double> >(2, vector<double> (3)));
  
  xi = 1.; // N-jettiness subtraction scheme variable
  // initialize constant coefficients for N-jettiness
  pi = M_PI; // FIXME: put correct value !!!
  cA = 3.;
  cF = 4./3.;
  tR = 1./2.;
  zeta2 = pi*pi/6.;
  zeta3 = 1.; // FIXME: put correct value !!!
  beta0 = -11./3. * cA + 2./3. * nf; // FIXME
  beta1 = 1.; // FIXME: put correct value !!!
  Gamma0 = 4.;
  Gamma1 = 4./3. * (cA * (4-pi*pi) + 5 * beta0); // = 4. * (cA * (67./9. - pi*pi/3) - 20./9. * tR * nf);
  // coefficients in Beam function definition
  Gamma_0[0] = Gamma0 * cA; // i = 0: gluon
  Gamma_0[1] = Gamma0 * cF; // i = 1: quark
  Gamma_1[0] = Gamma1 * cA; // i = 0: gluon
  Gamma_1[1] = Gamma1 * cF; // i = 1: quark
  gammaB0[0] = 2. * beta0;  // i = 0: gluon
  gammaB0[1] = 6. * cF;     // i = 1: quark
  gammaB1[0] = cA * (cA * (182./9-32.*zeta3) + beta0 * (94./9.-2.*pi*pi/3.)) + 2. * beta1; // i = 0: gluon
  gammaB1[1] = cF * (cA * (146./9.-80.*zeta3) + cF * (3.-4.*pi*pi+48.*zeta3) + beta0 * (121./9.+2.*pi*pi/3.)); // i = 1: quark
  // coefficients in Jet function definition
  gammaJ0[0] = gammaB0[0]; // same as for Beam functions
  gammaJ0[1] = gammaB0[1]; // same as for Beam functions
  gammaJ1[0] = gammaB1[0]; // same as for Beam functions
  gammaJ1[1] = gammaB1[1]; // same as for Beam functions
  // coefficients in soft function
  Ccoef[0] = cA; // i = 0: gluon
  Ccoef[1] = cF; // i = 1: quark
  Ccoef[2] = 0.; // to turn off contribution
  // compute colour algebra
  T2[0] = cA; // squared colour operator gluon
  T2[1] = cF; // squared colour operator quark
  // remember: 0=gluon, 1=quark, 2=nothing
  for (int i = 0; i<=1; i++){
    TiTj[i][i][2] = - T2[i];
  }
  for (int i = 0; i<=1; i++){
    for (int j = 0; j<=1; j++){
      for (int k = 0; k<=1; k++){
  	TiTj[i][j][k] = (T2[k] - T2[i] - T2[j])/2.;
      }
    }
  }
}
// }}}
// {{{ NJsubtraction::SoftF_1(int a, int b, int i)
double NJsubtraction::SoftF_1(int a, int b, int i) { // a,b,i corresponds to sum over partons (up to 1-jettiness)
// compute soft function
  if (a != b) {
    assert(false && "ERROR: a and b initial states must be equal, either both quark (=1) or gluon (=0).");
  }
  double SoftF_1 = - 2. * Gamma0 * (Ccoef[a]+Ccoef[b]+Ccoef[i]);
  return SoftF_1;
}
// }}}
// {{{ NJsubtraction::SoftF_0(int a, int b, int i, double sab, double sa1, double s1b))
double NJsubtraction::SoftF_0(int a, int b, int i, double sab, double sa1, double s1b) { // a,b,i corresponds to sum over partons (up to 1-jettiness)
// compute soft function
  if (a != b) {
    assert(false && "ERROR: a and b initial states must be equal, either both quark (=1) or gluon (=0).");
  }
  // compute invariants from momenta
  double SoftF_0 = - Gamma0 * 2. * TiTj[a][b][i] * log(sab); // T_a * T_b * ln(s_ab) + T_b * T_a * ln(s_ba)
  if (i != 2) {
    SoftF_0 += - Gamma0 * 2. * TiTj[a][b][i] * log(sa1); // T_a * T_1 * ln(s_a1) + T_1 * T_a * ln(s_1a)
    SoftF_0 += - Gamma0 * 2. * TiTj[a][b][i] * log(s1b); // T_1 * T_b * ln(s_1b) + T_b * T_1 * ln(s_b1)
  }
  return SoftF_0;
}
// }}}
// {{{ NJsubtraction::SoftF_m1(int a, int b, int i, double sab, double sa1, double s1b))
double NJsubtraction::SoftF_m1(int a, int b, int i, double sab, double sa1, double s1b) { // a,b,i corresponds to sum over partons (up to 1-jettiness)
// compute soft function
  if (a != b) {
    assert(false && "ERROR: a and b initial states must be equal, either both quark (=1) or gluon (=0).");
  }
  // compute invariants from momenta
  double SUMoverIijm = 0.; // not needed for 0-jettiness; FIXME: implement for 1-jettiness
  double SoftF_m1 = 2. * TiTj[a][b][i] * (pow(log(sab),2) + zeta2 + 4. * SUMoverIijm); // T_a * T_b * (...) + T_b * T_a * (...)
  if (i != 2) {
    SoftF_m1 += 2. * TiTj[a][b][i] * (pow(log(sa1),2) + zeta2 + 4. * SUMoverIijm); // T_a * T_1 * (...) + T_1 * T_a * (...)
    SoftF_m1 += 2. * TiTj[a][b][i] * (pow(log(s1b),2) + zeta2 + 4. * SUMoverIijm); // T_1 * T_b * (...) + T_b * T_1 * (...)
  }
  return SoftF_m1;
}
// }}}
// {{{ NJsubtraction::calculate_C11(double pdf_factor_x1x2)
double NJsubtraction::calculate_C11(double pdf_factor_x1x2) {
// coefficient computed without LO squared amplitude
  int i = 1; // quark contribution
  int a = 1; // parton a = quark
  int b = 1; // parton b = quark
  int j = 2; // parton j (first non initial parton) does not exist (= 2)
  double C11 = 2. * Gamma_0[i] * pdf_factor_x1x2; // B_a^(1) * f_b + f_a * B_b^(1) = 2 * Gamma0[i] * f_a * f_b
  C11 =  C11 + SoftF_1(a,b,j)  * pdf_factor_x1x2; // S^(1)_n * f_a * f_b
  return C11;
}
// }}}
// {{{ NJsubtraction::calculate_C10(double pdf_factor_x1x2,double convolutions_Pxx,double convolutions_Ixx))
double NJsubtraction::calculate_C10(double pdf_factor_x1x2,double convolutions_Pxx,double convolutions_Ixx) {// FIXME: pass over momenta or sab invariants !!!
// coefficient computed without LO squared amplitude
  double sab = 1.; // FIXME: remove hardcoded !!!
  double sa1 = 1.; // FIXME: remove hardcoded !!!
  double s1b = 1.; // FIXME: remove hardcoded !!!
  int i = 1; // quark contribution
  int a = 1; // parton a = quark
  int b = 1; // parton b = quark
  int j = 2; // parton j (first non initial parton) does not exist (= 2)
  double C10 = 2. * (-gammaB0[i]/2. * pdf_factor_x1x2 + 2. * convolutions_Pxx); // B_a^(1) * f_b + f_a * B_b^(1) = 2 * (...)
// put in extra muR part:  C10 =  C10 + Gamma_0[i] * log(lambda_a) + Gamma_0[i] * log(lambda_b); // lambda dependend contribution in eq. (A.16) of 1505:04794
  C10 =  C10 + SoftF_0(a,b,j,sab,sa1,s1b) * pdf_factor_x1x2; // S^(1)_n * f_a * f_b
  return C10;
}
// }}}
// {{{ NJsubtraction::calculate_C10_muR(double pdf_factor_x1x2,double convolutions_Pxx,double convolutions_Ixx,double Q_a,double Q_b,double muR))
double NJsubtraction::calculate_C10_muR(double pdf_factor_x1x2,double convolutions_Pxx,double convolutions_Ixx,double Q_a,double Q_b,double muR) {// FIXME: pass over momenta or sab invariants !!!
// coefficient computed without LO squared amplitude
  int i = 1; // quark contribution
  int a = 1; // parton a = quark
  int b = 1; // parton b = quark
  int j = 2; // parton j (first non initial parton) does not exist (= 2)
  double Llambda_a = log(Q_a/muR * xi);
  double Llambda_b = log(Q_b/muR * xi);
  double LlambdaR  = log(xi/muR);
  double jetfunction = 0.; // not needed for 0-jettiness; FIXME: implement for 1-jettiness
  double C10 = jetfunction;
// put in muR independent part:  C10 =  C10 + 2. * (-gammaB0[i]/2. * pdf_factor_x1x2 + 2. * convolutions_Pxx); // B_a^(1) * f_b + f_a * B_b^(1) = 2 * (...)
  C10 =  C10 + Gamma_0[i] * (Llambda_a+Llambda_b) * pdf_factor_x1x2; // lambda dependend contribution in eq. (A.16) of 1505:04794
// put in muR independent part:  C10 =  C10 + SoftF_0(a,b,j,sab,sa1,s1b) * pdf_factor_x1x2; // S^(1)_n * f_a * f_b
  C10 =  C10 + SoftF_1(a,b,j) * LlambdaR * pdf_factor_x1x2; // lambda dependend contribution in eq. (A.24) of 1505:04794
  return C10;
}
// }}}
// {{{ NJsubtraction::calculate_C1m1(double pdf_factor_x1x2,double convolutions_Pxx,double convolutions_Ixx)
double NJsubtraction::calculate_C1m1(double pdf_factor_x1x2,double convolutions_Pxx,double convolutions_Ixx) {// FIXME: pass over momenta or sab invariants !!!
// coefficient computed without LO squared amplitude
  double sab = 1.; // FIXME: remove hardcoded !!!
  double sa1 = 1.; // FIXME: remove hardcoded !!!
  double s1b = 1.; // FIXME: remove hardcoded !!!
  //  int i = 1; // quark contribution
  int a = 1; // parton a = quark
  int b = 1; // parton b = quark
  int j = 2; // parton j (first non initial parton) does not exist (= 2)
  double C1m1 = 4. * convolutions_Ixx; // B_a^(1) * f_b + f_a * B_b^(1) --> 2 * 2*Ixx
  C1m1 = C1m1 + SoftF_m1(a,b,j,sab,sa1,s1b) * pdf_factor_x1x2; // S^(1)_n * f_a * f_b
  return C1m1;
}
// }}}
// {{{ NJsubtraction::calculate_C1m1_muR(double pdf_factor_x1x2,double convolutions_Pxx,double convolutions_Ixx,double,double lnRF Q_aOVERmuR,double Q_bOVERmuR)
double NJsubtraction::calculate_C1m1_muR(double pdf_factor_x1x2,double HardF1,double convolutions_Pxx,double convolutions_Ixx,double lnRF,double Q_a,double Q_b,double muR) {// FIXME: pass over momenta or sab invariants !!!
// coefficient computed without LO squared amplitude
  double sab = 1.; // FIXME: remove hardcoded !!!
  double sa1 = 1.; // FIXME: remove hardcoded !!!
  double s1b = 1.; // FIXME: remove hardcoded !!!
  int i = 1; // quark contribution
  int a = 1; // parton a = quark
  int b = 1; // parton b = quark
  int j = 2; // parton j (first non initial parton) does not exist (= 2)
  double Llambda_a = log(Q_a/muR * xi);
  double Llambda_b = log(Q_b/muR * xi);
  double LlambdaR  = log(xi/muR);
  double jetfunction = 0.; // not needed for 0-jettiness; FIXME: implement for 1-jettiness
  double C1m1 = HardF1 * pdf_factor_x1x2;
  C1m1 = C1m1 + jetfunction;
  // put in muR independent part: C1m1 = C1m1 + 4. * convolutions_Ixx // B_a^(1) * f_b + f_a * B_b^(1) --> 2 * 2*Ixx
  C1m1 = C1m1 + lnRF * 4. * convolutions_Pxx; // B_a^(1) * f_b + f_a * B_b^(1) --> 2 * 2*Pxx
  C1m1 = C1m1 + (-gammaB0[i]/2. * pdf_factor_x1x2 + 2. * convolutions_Pxx) * (Llambda_a+Llambda_b); // lambda dependend contribution in eq. (A.16) of 1505:04794
  C1m1 = C1m1 + Gamma_0[i] * (pow(Llambda_a,2)/2.+pow(Llambda_b,2)/2.) * pdf_factor_x1x2; // lambda dependend contribution in eq. (A.16) of 1505:04794
  // put in muR independent part: C1m1 =  C1m1 + SoftF_0(a,b,j,sab,sa1,s1b) * pdf_factor_x1x2; // S^(1)_n * f_a * f_b
  C1m1 = C1m1 + SoftF_0(a,b,j,sab,sa1,s1b) * LlambdaR * pdf_factor_x1x2; // S^(1)_n * f_a * f_b
  C1m1 = C1m1 + SoftF_1(a,b,j) * pow(LlambdaR,2)/2. * pdf_factor_x1x2; // lambda dependend contribution in eq. (A.24) of 1505:047945D
  return C1m1;
}
// }}}
// {{{ NJsubtraction::computeLLnIntegrals(double x_min,double x_steps,int n_steps,vector<double> &LL0_int,vector<double> &LL1_int,vector<double> &LL2_int,vector<double> &LL3_int)
void NJsubtraction::computeLLnIntegrals(double x_min,double x_steps,int n_steps,vector<double> &LL0_int,vector<double> &LL1_int,vector<double> &LL2_int,vector<double> &LL3_int) {
  for (int i=0; i<n_steps; i++) {
    double x_cut=x_min+i*x_steps;
    LL0_int[i] = log(x_cut);
    LL1_int[i] = pow(log(x_cut),2)/2.;
    LL2_int[i] = pow(log(x_cut),3)/3.;
    LL3_int[i] = pow(log(x_cut),4)/4.;
  }
}
// }}}
// {{{ NJsubtraction::compute_NJ_log(double TauNoverQ,int order,int n)
double NJsubtraction::compute_NJ_log(double TauNoverQ,int order,int n) {
// compute the differntial logs at each order
  if (TauNoverQ >= 0) {
    if (n==0) {
      return 1./TauNoverQ;
    }
    else {
      return pow(log(TauNoverQ),n)/TauNoverQ;
    }
  }
  else { // theta function in definition of eq. (2.8) in 1505:04794,
         // but can this ever happen? shouldn't there be an assert if TauNoverQ < 0?
    return 0.;
  }
}
// }}}
// {{{ NJsubtraction::Plusdd(int n,double z)
double NJsubtraction::Plusdd(int n,double x) {
  double Plusdd = 1/x;
  if (n > 0) {
    Plusdd *= pow(log(x),n);
  }
  return Plusdd;
}
// }}}
// {{{ NJsubtraction::I1qq_plus(double z)
double NJsubtraction::I1qq_plus(double z) {
  double I1qq_plus = cF * 2. * Plusdd(1,1-z);
  return I1qq_plus;
}
// }}}
// {{{ NJsubtraction::I1qq_z(double z)
double NJsubtraction::I1qq_z(double z) {
  double I1qq_z = cF * (-log(1-z)*(1+z)+1-z-(1+z*z)/(1-z)*log(z)); // first term: left over from plus distribution
  return I1qq_z;
}
// }}}
// {{{ NJsubtraction::I1qq_delta(double z)
double NJsubtraction::I1qq_delta(double z) {
  double I1qq_delta = -zeta2;
  return I1qq_delta;
}
// }}}
// {{{ NJsubtraction::I1qq_int(double z)
double NJsubtraction::I1qq_int(double z) {
  double I1qq_int = 0.;
  return I1qq_int;
}
// }}}
// {{{ NJsubtraction::I1qg_z(double z)
double NJsubtraction::I1qg_z(double z) {
  double I1qg_z = P1qg_z(z) * (log((1-z)/z)-1) + 1; // last term: actually theta(1-z), but this can never happen, can it?
  return I1qg_z;
}
// }}}
// {{{ NJsubtraction::I1gg_plus(double z)
// double NJsubtraction::I1gg_plus(double z) {
//   double I1gg_plus = 0.;
//   return I1gg_plus;
// }
// // }}}
// // {{{ NJsubtraction::I1gg_z(double z)
// double NJsubtraction::I1gg_z(double z) {
//   double I1gg_z = 0.;
//   return I1gg_z;
// }
// // }}}
// // {{{ NJsubtraction::I1gg_delta(double z)
// double NJsubtraction::I1gg_delta(double z) {
//   double I1gg_delta = 0.;
//   return I1gg_delta;
// }
// // }}}
// // {{{ NJsubtraction::I1gg_int(double z)
// double NJsubtraction::I1gg_int(double z) {
//   double I1gg_int = 0.;
//   return I1gg_int;
// }
// }}}
// {{{ NJsubtraction::P1qq_plus(double z)
double NJsubtraction::P1qq_plus(double z) {
  double P1qq_plus = cF * 2. * Plusdd(0,1-z);
  return P1qq_plus;
}
// }}}
// {{{ NJsubtraction::P1qq_z(double z)
double NJsubtraction::P1qq_z(double z) {
  double P1qq_z = cF * (-1-z);
  return P1qq_z;
}
// }}}
// {{{ NJsubtraction::P1qq_delta(double z)
double NJsubtraction::P1qq_delta(double z) {
  double P1qq_delta = cF * 3./2.;
  return P1qq_delta;
}
// }}}
// {{{ NJsubtraction::P1qq_int(double z)
double NJsubtraction::P1qq_int(double z) {
  double P1qq_int = cF * 2. * (-log(1-z));
  return P1qq_int;
}
// }}}
// {{{ NJsubtraction::P1qg_z(double z)
double NJsubtraction::P1qg_z(double z) {
  double P1qg_z = tR * (1-2.*z*(1-z));
  return P1qg_z;
}
// }}}
// {{{ NJsubtraction::P1gq_z(double z)
double NJsubtraction::P1gq_z(double z) {
  double P1gq_z = cF * (1+(1-z)*(1-z))/z;
  return P1gq_z;
}
// }}}
// {{{ NJsubtraction::P1gg_plus(double z)
double NJsubtraction::P1gg_plus(double z) {
  double P1gg_plus = 2. * cA * Plusdd(0,1-z);
  return P1gg_plus;
}
// }}}
// {{{ NJsubtraction::P1gg_z(double z)
double NJsubtraction::P1gg_z(double z) {
  double P1gg_z = 2. * cA * (1/z-2+z+z*z);
  return P1gg_z;
}
// }}}
// {{{ NJsubtraction::P1gg_delta(double z)
double NJsubtraction::P1gg_delta(double z) {
  double P1gg_delta = cA * beta0;
  return P1gg_delta;
}
// }}}
// {{{ NJsubtraction::P1gg_int(double z)
double NJsubtraction::P1gg_int(double z) {
  double P1gg_int = 2. * cA * (-log(1-z));
  return P1gg_int;
}
// }}}
