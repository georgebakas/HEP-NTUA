#include "../include/classes.cxx"
double P_qg(double x){return C_F * (1. + pow(1. - x, 2)) / x;}

double P_gq(double x){return T_R * (pow(x, 2) + pow(1. - x, 2));}

// CS convention
double P_qq_reg(double x){return -C_F * (1. + x);}
double P_qq(double x){return 0.;}
double P_qq_plus(double x){return C_F * (1. + pow(x, 2)) / (1. - x);}
double P_qq_delta(){return 0.;}
double intP_qq_plus(double x){return C_F * (-.5 * pow(x, 2) - x - 2 * log(1. - x));}

/*
// CS convention
double Pxx_qq_reg(double x){return -C_F * (1. + x);}
double Pxx_qq(double x){return 0.;}
double Pxx_qq_plus(double x){return C_F * (1. + pow(x, 2)) / (1. - x);}
double Pxx_qq_delta(){return 0.;}
double intPxx_qq_plus(double x){return C_F * (-.5 * pow(x, 2) - x - 2 * log(1. - x));}
*/

// CDST convention
double Pxx_qq(double x){return -C_F * (1. + x);}
double Pxx_qq_reg(double x){return -C_F * (1. + x);}
double Pxx_qq_plus(double x){return 2 * C_F * (1. / (1. - x));}
double Pxx_qq_delta(){return 3. / 2. * C_F;}
double intPxx_qq_plus(double x){return -2 * C_F * log(1. - x);}


double P_gg(double x){return 2 * C_A * ((1. - x) / x - 1. + x * (1. - x));}
double P_gg_reg(double x){return 2 * C_A * ((1. - x) / x - 1. + x * (1. - x));}
double P_gg_plus(double x){return 2 * C_A / (1. - x);}
double P_gg_delta(int N_f){return (11. / 6. * C_A - 2. / 3. * N_f * T_R);}
double intP_gg_plus(double x){return -2 * C_A * log(1. - x);}

double Kbar_qg(double x){return P_qg(x) * log((1. - x) / x) + C_F * x;}

double Kbar_gq(double x){return P_gq(x) * log((1. - x) / x) + 2 * T_R * x * (1. - x);}

double Kbar_qq(double x){return C_F * (-(1. + x) * log((1. - x) / x) + (1. - x));}
double Kbar_qq_plus(double x){return C_F * ((2 / (1. - x)) * log((1. - x) / x));}
double Kbar_qq_delta(){return -C_F * (5. - pi2);}
double intKbar_qq_plus(double x){return C_F * (-pow(log(1. - x), 2) - 2 * gsl_sf_dilog(1. - x) + pi2_3);}

double Kbar_gg(double x){return 2 * C_A * ((1. - x) / x - 1. + x * (1. - x)) * log((1. - x) / x);}
double Kbar_gg_plus(double x){return C_A * ((2 / (1. - x)) * log((1. - x) / x));}
double Kbar_gg_delta(int N_f){return -(C_A * ((50. / 9.) - pi2) - T_R * N_f * (16. / 9.));}
double intKbar_gg_plus(double x){return C_A * (-pow(log(1. - x), 2) - 2 * gsl_sf_dilog(1. - x) + pi2_3);}

double Kt_qg(double x){return P_qg(x) * log(1. - x);}

double Kt_gq(double x){return P_gq(x) * log(1. - x);}

double Kt_qq(double x){return P_qq_reg(x) * log(1. - x);}
double Kt_qq_plus(double x){return C_F * (2 / (1. - x)) * log(1. - x);}
double Kt_qq_delta(){return C_F * (-pi2 / 3.);}
double intKt_qq_plus(double x){return -C_F * pow(log(1. - x), 2);}

double Kt_gg(double x){return P_gg_reg(x) * log(1. - x);}
double Kt_gg_plus(double x){return C_A * (2 / (1. - x)) * log(1. - x);}
double Kt_gg_delta(){return C_A * (-pi2 / 3.);}
double intKt_gg_plus(double x){return -C_A * pow(log(1. - x), 2);}


double gamma_g(int N_f){return (11. / 6.) * C_A - (2. / 3.) * T_R * N_f;}
double gamma_g_ferm(int N_f){return - (2. / 3.) * T_R * N_f;}

double K_g (int N_f){return (67. / 18. - pi2_6) * C_A - (10. / 9.) * T_R * N_f;}
double K_g_ferm (int N_f){return  - (10. / 9.) * T_R * N_f;}



double gamma_a(int N_f){
  double sum_f_charge2 = 1. + 1. + 1.; // leptons
  sum_f_charge2 += N_c * (2 * (4. / 9.) + 2 * (1. / 9.)); // 2 light quark generations
  if (N_f > 4){sum_f_charge2 += N_c * (1. / 9.);} // bottom quark
  if (N_f > 5){sum_f_charge2 += N_c * (4. / 9.);} // top quark
  return -(2. / 3.) * sum_f_charge2;
} // C_A -> 0, T_R -> sum_f N_c,f * Q²_f   ??? T_R -> 1 before...
//double gamma_a(int N_f){return -(2. / 3.) * N_f * N_c;} // C_A -> 0, T_R -> sum_f N_c,f * Q²_f   ??? T_R -> 1 before...

double gamma_qew_a(int N_f){
  double sum_f_charge2 = 1. + 1. + 1.; // leptons
  sum_f_charge2 += N_c * (2 * (4. / 9.) + 2 * (1. / 9.)); // 2 light quark generations
  if (N_f > 4){sum_f_charge2 += N_c * (1. / 9.);} // bottom quark
  if (N_f > 5){sum_f_charge2 += N_c * (4. / 9.);} // top quark
  return -(2. / 3.) * sum_f_charge2;
} // C_A -> 0, T_R -> sum_f N_c,f * Q²_f   ??? T_R -> 1 before...
//double gamma_qew_a(int N_f){return -(2. / 3.) * N_f * N_c;} // C_A -> 0, T_R -> sum_f N_c,f * Q²_f   ??? T_R -> 1 before...


double K_a(int N_f){
  double sum_f_charge2 = 1. + 1. + 1.; // leptons
  sum_f_charge2 += N_c * (2 * (4. / 9.) + 2 * (1. / 9.)); // 2 light quark generations
  if (N_f > 4){sum_f_charge2 += N_c * (1. / 9.);} // bottom quark
  if (N_f > 5){sum_f_charge2 += N_c * (4. / 9.);} // top quark
  return -(10. / 9.) * sum_f_charge2;
} // C_A -> 0, T_R -> sum_f N_c,f * Q²_f   ??? T_R -> 1 before...
//double K_a(int N_f){return -(10. / 9.) * N_f * N_c;} // C_A -> 0, T_R -> sum_f N_c,f * Q²_f   ??? T_R -> 1 before...

double K_qew_a(int N_f){
  double sum_f_charge2 = 1. + 1. + 1.; // leptons
  sum_f_charge2 += N_c * (2 * (4. / 9.) + 2 * (1. / 9.)); // 2 light quark generations
  if (N_f > 4){sum_f_charge2 += N_c * (1. / 9.);} // bottom quark
  if (N_f > 5){sum_f_charge2 += N_c * (4. / 9.);} // top quark
  return -(10. / 9.) * sum_f_charge2;
} // C_A -> 0, T_R -> sum_f N_c,f * Q²_f   ??? T_R -> 1 before...
//double K_qew_a(int N_f){return -(10. / 9.) * N_f * N_c;} // C_A -> 0, T_R -> sum_f N_c,f * Q²_f   ??? T_R -> 1 before...



double P_qew_qa(double x){return (1. + pow(1. - x, 2)) / x;} // C_F -> Q²_f

double P_qew_aq(double x){return N_c * (pow(x, 2) + pow(1. - x, 2));} // T_R -> N_c,f * Q²_f

// CS convention
double P_qew_qq_reg(double x){return -(1. + x);} // C_F -> Q²_f
double P_qew_qq(double x){return 0.;} // C_F -> Q²_f
double P_qew_qq_plus(double x){return (1. + pow(x, 2)) / (1. - x);} // C_F -> Q²_f
double P_qew_qq_delta(){return 0.;} // C_F -> Q²_f
double intP_qew_qq_plus(double x){return -.5 * pow(x, 2) - x - 2 * log(1. - x);} // C_F -> Q²_f


// CS convention
double Pxx_qew_qq_reg(double x){return -(1. + x);} // C_F -> Q²_f
double Pxx_qew_qq(double x){return 0.;} // C_F -> Q²_f
double Pxx_qew_qq_plus(double x){return (1. + pow(x, 2)) / (1. - x);} // C_F -> Q²_f
double Pxx_qew_qq_delta(){return 0.;} // C_F -> Q²_f
double intPxx_qew_qq_plus(double x){return -.5 * pow(x, 2) - x - 2 * log(1. - x);} // C_F -> Q²_f

/*
// CDST convention
double Pxx_qew_qq(double x){return -(1. + x);} // C_F -> Q²_f
double Pxx_qew_qq_reg(double x){return -(1. + x);} // C_F -> Q²_f
double Pxx_qew_qq_plus(double x){return 2. / (1. - x);} // C_F -> Q²_f
double Pxx_qew_qq_delta(){return 3. / 2.;} // C_F -> Q²_f
double intPxx_qew_qq_plus(double x){return -2 * log(1. - x);} // C_F -> Q²_f
*/

double P_qew_aa(double x){return 0.;} // C_A -> 0
//double P_qew_aa(double x){return 2 * 0. * ((1. - x) / x - 1. + x * (1. - x));} // C_A -> 0

double P_qew_aa_reg(double x){return 0.;} // C_A -> 0
//double P_qew_aa_reg(double x){return 2 * 0. * ((1. - x) / x - 1. + x * (1. - x));} // C_A -> 0

double P_qew_aa_plus(double x){return 0.;} // C_A -> 0
//double P_qew_aa_plus(double x){return 2 * 0. / (1. - x);} // C_A -> 0

double P_qew_aa_delta(int N_f){
  double sum_f_charge2 = 1. + 1. + 1.; // leptons
  sum_f_charge2 += N_c * (2 * (4. / 9.) + 2 * (1. / 9.)); // 2 light quark generations
  if (N_f > 4){sum_f_charge2 += N_c * (1. / 9.);} // bottom quark
  if (N_f > 5){sum_f_charge2 += N_c * (4. / 9.);} // top quark
  return (-2. / 3.) * sum_f_charge2;
} // C_A -> 0, T_R -> sum_f N_c,f * Q²_f
//double P_qew_aa_delta(int N_f){return (11. / 6. * 0. - 2. / 3. * N_f * N_c);} // C_A -> 0, T_R -> sum_f N_c,f * Q²_f

double intP_qew_aa_plus(double x){return 0.;} // C_A -> 0
//double intP_qew_aa_plus(double x){return -2 * 0. * log(1. - x);} // C_A -> 0



double Kbar_qew_qa(double x){return P_qew_qa(x) * log((1. - x) / x) + x;} // C_F -> Q²_f

double Kbar_qew_aq(double x){return P_qew_aq(x) * log((1. - x) / x) + 2 * N_c * x * (1. - x);} // T_R -> N_c,f * Q²_f

double Kbar_qew_qq(double x){return (-(1. + x) * log((1. - x) / x) + (1. - x));} // C_F -> Q²_f
double Kbar_qew_qq_plus(double x){return ((2 / (1. - x)) * log((1. - x) / x));} // C_F -> Q²_f
double Kbar_qew_qq_delta(){return -(5. - pi2);} // C_F -> Q²_f
double intKbar_qew_qq_plus(double x){return (-pow(log(1. - x), 2) - 2 * gsl_sf_dilog(1. - x) + pi2_3);} // C_F -> Q²_f

double Kbar_qew_aa(double x){return 0.;}
//double Kbar_qew_aa(double x){return 2 * 0. * ((1. - x) / x - 1. + x * (1. - x)) * log((1. - x) / x);} // C_A -> 0

double Kbar_qew_aa_plus(double x){return 0.;}
//double Kbar_qew_aa_plus(double x){return 0. * ((2 / (1. - x)) * log((1. - x) / x));} // C_A -> 0

double Kbar_qew_aa_delta(int N_f){
  double sum_f_charge2 = 1. + 1. + 1.; // leptons
  sum_f_charge2 += N_c * (2 * (4. / 9.) + 2 * (1. / 9.)); // 2 light quark generations
  if (N_f > 4){sum_f_charge2 += N_c * (1. / 9.);} // bottom quark
  if (N_f > 5){sum_f_charge2 += N_c * (4. / 9.);} // top quark
  cout << "sum_f_charge2 = " << sum_f_charge2 << endl;
  return (16. / 9.) * sum_f_charge2;
} // C_A -> 0, T_R * N_f -> sum_f N_c,f * Q²_f
//double Kbar_qew_aa_delta(int N_f){return -(0. * ((50. / 9.) - pi2) - N_c * N_f * (16. / 9.));} // C_A -> 0, T_R -> sum_f N_c,f * Q²_f

double intKbar_qew_aa_plus(double x){return 0.;}
//double intKbar_qew_aa_plus(double x){return 0. * (-pow(log(1. - x), 2) - 2 * gsl_sf_dilog(1. - x) + pi2_3);} // C_A -> 0



double Kt_qew_qa(double x){return P_qew_qa(x) * log(1. - x);} // C_F -> Q²_f

double Kt_qew_aq(double x){return P_qew_aq(x) * log(1. - x);} // T_R -> N_c,f * Q²_f

double Kt_qew_qq(double x){return P_qew_qq_reg(x) * log(1. - x);} // C_F -> Q²_f
double Kt_qew_qq_plus(double x){return (2 / (1. - x)) * log(1. - x);} // C_F -> Q²_f
double Kt_qew_qq_delta(){return (-pi2 / 3.);} // C_F -> Q²_f
double intKt_qew_qq_plus(double x){return -pow(log(1. - x), 2);} // C_F -> Q²_f

double Kt_qew_aa(double x){return P_qew_aa_reg(x) * log(1. - x);} // C_A -> 0, T_R -> N_c,f * Q²_f
double Kt_qew_aa_plus(double x){return 0. * (2 / (1. - x)) * log(1. - x);} // C_A -> 0
double Kt_qew_aa_delta(){return 0. * (-pi2 / 3.);} // C_A -> 0
double intKt_qew_aa_plus(double x){return -0. * pow(log(1. - x), 2);} // C_A -> 0
