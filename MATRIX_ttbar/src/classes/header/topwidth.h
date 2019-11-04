#ifndef TOPWIDTH_H
#define TOPWIDTH_H

using namespace std;

class topwidth{
 private:

 public:
////////////////////
//  constructors  //
////////////////////
  topwidth();
  topwidth(model_set & msi);

  model_set *msi;
  phasespace_set *psi;
  observable_set *osi;

  double CE_f(double & omega2, double & beta2);
  double CE_z_m(double & omega);
  double CE_P0(double & omega2, double & z);
  double CE_P3(double & omega2, double & z);
  double CE_W0(double & omega2, double & z);
  double CE_Pplus(double & omega2, double & z);
  double CE_Pminus(double & omega2, double & z);
  double CE_Yp(double & omega2, double & z);
  double CE_Wplus(double & omega2, double & z);
  double CE_Wminus(double & omega2, double & z);
  double CE_Yw(double & omega2, double & z);
  double CE_Gamma_0(double & Gamma_t_infinity, double & omega2, double & beta2);
  double CE_Gamma_1(double & Gamma_t_infinity, double & omega2, double & beta2);
  double c_propagator_Breit_Wigner(double & r, int m, double smin, double smax);
  double g_propagator_Breit_Wigner(double & s, int m, double smin, double smax);
  double integrate_topwidth(double & alpha_S_t, int order_alpha_S, double relative_accuracy);

};
#endif
