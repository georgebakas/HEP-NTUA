namespace mathconst
{
  const double pi = 4. * atan(1.);
  const double pi_48 = pi / 48.;
  const double f2pi = 2. * pi;
  const double f4pi = 4. * pi;
  const double pi2 = pow(pi, 2);
  const double pi2_2 = pi2 / 2.;
  const double pi2_3 = pi2 / 3.;
  const double pi2_6 = pi2 / 6.;
  const double pi2_9 = pi2 / 9.;
  const double pi2_12 = pi2 / 12.;
  const double pi2_18 = pi2 / 18.;
  const double pi2_24 = pi2 / 24.;
  const double pi2_48 = pi2 / 48.;
  const double pi2_72 = pi2 / 72.;
  const double pi2_144 = pi2 / 144.;
  const double pi4 = pow(pi, 4);
  const double pi4_90 = pi4 / 90.;
  const double pi4_180 = pi4 / 180.;
  const double pi4_480 = pi4 / 480.;
  const double sqrt2 = sqrt(2.);
  const double inv2pi = 1. / (2. * pi);
  const double inv3pi = 1. / (3. * pi);
  const double inv4pi = 1. / (4. * pi);
  const double inv6pi = 1. / (6. * pi);
  const double inv12pi = 1. / (12. * pi);
  const double Eulergamma = 0.5772156649015328606065121;
  const double b0_Eulergamma = 2 * exp(-Eulergamma);
  const double zeta2 = pi2_6;
  const double zeta3 = 1.20205690315959429;
  const double zeta4 = pi4_90;
  const double PolyGamma_21 = -2 * zeta3;
  const double_complex ri(0, 1.);
};
namespace physconst
{
  using namespace mathconst;
  const fourvector nullvector(0., 0., 0., 0.);

  // colour algebra
  const double infinity = 1.e99;
  const int N_c = 3;
  const double dN_c=(double)N_c;
  const double C_F = 4. / 3.;
  const double C_A = double(N_c);//3.;
  const double T_R = 1. / 2.;

  const double gamma_q = (3. / 2.) * C_F;
  const double K_q = (7. / 2. - pi2_6) * C_F;
  const double gamma_g_bos = (11. / 6.) * C_A;
  const double K_g_bos = (67. / 18. - pi2_6) * C_A;

  const double gamma_qew_q = 3. / 2.; // C_F -> Q²_f
  const double K_qew_q = 7. / 2. - pi2_6; // C_F -> Q²_f

  //  const double gamma_g = (11. / 6.) * C_A - (2. / 3.) * T_R * N_f;
  //  const double K_g = (67. / 18. - pi2_6) * C_A - (10. / 9.) * T_R * N_f;
  /*
  const double psf2 = 1. / (2. * pow(2. * pi, 2.));
  const double psf3 = 1. / (2. * pow(2. * pi, 5.));
  // constant phase-space factor:  (hbar * c [GeV cm]) ^ 2 * fbarn [cm^2] --> result in fbarn

  const double hcf2 = pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 2.));
  const double hcf3 = pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 5.));
  const double hcf4 = pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 8.));
  const double hcf5 = pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 11.));
  const double hcf6 = pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 14.));
  const double hcf7 = pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 17.));

  const double hcf_n_particle[21] = {1., 
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 1 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 2 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 3 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 4 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 5 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 6 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 7 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 8 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 9 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 10 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 11 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 12 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 13 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 14 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 15 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 16 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 17 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 18 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 19 - 4)),
				     pow(0.197326968E-13, 2.) * 1.e+39 / (2. * pow(2. * pi, 3 * 20 - 4))};
  */
  /*
  const double Q_u = 2. / 3.;
  const double Q_d = -1. / 3.;
  const double Q_n = 0.;
  const double Q_l = -1.;
  const double Iw_u = 1. / 2.;
  const double Iw_d = -1. / 2.;
  const double Iw_n = 1. / 2.;
  const double Iw_l = -1. / 2.;
  const double C_gamma_ee = 1.;
  const double C_gammaWminusWplus = -1.;
  */
  const double C_qA = -2.;
  const double C_qF = 1. / 2.;
  const double C_gF = -9. / 2.;
};
