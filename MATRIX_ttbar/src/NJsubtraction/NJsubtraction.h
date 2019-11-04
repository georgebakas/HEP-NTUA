#ifndef NJSUBRACTION_H
#define NJSUBRACTION_H

#include <iostream>
#include <vector>

using namespace std;

class NJsubtraction {
 private:
  // Beam function coeficients
  vector<double> Gamma_0, Gamma_1,gammaB0,gammaB1,gammaJ0,gammaJ1,Ccoef,T2; // i = 0: gluon; i = 1: quark
  vector<vector<vector<double> > > TiTj;
  double pi,cA,cF,tR;
  double Gamma0,Gamma1;
  double Plusdd(int,double);
 public:
  double xi; // N-jettiness subtraction scheme
  double sab,sa1,s1b; // set from outside
  double nf,beta0,beta1,zeta2,zeta3;
  NJsubtraction(); // default constructor
  double compute_NJ_log(double,int,int);
  // soft function
  double SoftF_1(int,int,int); // commputable up to
  double SoftF_0(int,int,int,double,double sa1=1.,double s1b=1.); // commputable up to
  double SoftF_m1(int,int,int,double,double sa1=1.,double s1b=1.); // commputable up to
  // coefficient functions for subtraction
  double calculate_C11(double);
  double calculate_C10(double,double,double);
  double calculate_C10_muR(double,double,double,double,double,double);
  double calculate_C1m1(double,double,double);
  double calculate_C1m1_muR(double,double,double,double,double,double,double,double);
  // Ixx(z) functions inside beam function;
  // split as plus-distribution part (no additional z dependence), z-dependent part and delta(1-z) piece
  // all functions return zero except if their component exists
  // additional the integral of each function from 0 to z with z < 1
  double I1qq_plus(double z);
  double I1qq_z(double z);
  double I1qq_delta(double z);
  double I1qq_int(double z);
  double I1qg_plus(double z) {return 0.;};
  double I1qg_z(double z);
  double I1qg_delta(double z) {return 0.;};
  double I1qg_int(double z) {return 0.;};
  double I1gq_plus(double z) {return 0.;};
  double I1gq_z(double z) {return 0.;};
  double I1gq_delta(double z) {return 0.;};
  double I1gq_int(double z) {return 0.;};
  double I1gg_plus(double z) {return 0.;};
  double I1gg_z(double z) {return 0.;};
  double I1gg_delta(double z) {return 0.;};
  double I1gg_int(double z) {return 0.;};
  // same notation for all first-order splitting functions with alpha_S/2pi normalization
  double P1qq_plus(double z);
  double P1qq_z(double z);
  double P1qq_delta(double z);
  double P1qq_int(double z);
  double P1qg_plus(double z) {return 0.;};
  double P1qg_z(double z);
  double P1qg_delta(double z) {return 0.;};
  double P1qg_int(double z) {return 0.;};
  double P1gq_plus(double z) {return 0.;};
  double P1gq_z(double z);
  double P1gq_delta(double z) {return 0.;};
  double P1gq_int(double z){return 0.;};
  double P1gg_plus(double z);
  double P1gg_z(double z);
  double P1gg_delta(double z);
  double P1gg_int(double z);
  // Integrals of N-jettiness logarithms
  void computeLLnIntegrals(double,double,int,vector<double> &,vector<double> &,vector<double> &,vector<double> &);
};

#endif // NJSUBTRACTION_H
