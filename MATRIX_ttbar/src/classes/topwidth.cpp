#include "../include/classes.cxx"
//#include "../include/definitions.cxx"
#include "../include/definitions.observable.set.cxx"
#include "../include/definitions.phasespace.set.cxx"

topwidth::topwidth(){
  Logger logger("topwidth::topwidth");
  logger << LOG_DEBUG << "started" << endl;

  logger << LOG_DEBUG << "finished" << endl;
}


// probably better without phasespace_set and observable_set !!!
//phasespace_set & _psi, observable_set & _osi
topwidth::topwidth(model_set & _msi){
  Logger logger("topwidth::topwidth");
  logger << LOG_DEBUG << "started" << endl;

  msi = &_msi;
  //  alpha_S_t = _alpha_S_t;
  //  order_alpha_S = _order_alpha_S;
  //  relative_accuracy = _relative_accuracy;

  //  psi = &_psi;
  //  osi = &_osi;

  logger << LOG_DEBUG << "finished" << endl;
}



double topwidth::CE_f(double & omega2, double & beta2){
  return pow(1. - beta2, 2) + omega2 * (1. + beta2) - 2 * pow(omega2, 2);
}
double topwidth::CE_z_m(double & omega){
  return pow(1. - omega, 2);
}
double topwidth::CE_P0(double & omega2, double & z){
  return .5 * (1. - omega2 + z);
}
double topwidth::CE_P3(double & omega2, double & z){
  return .5 * sqrt(lambda(1., omega2, z));
}
double topwidth::CE_W0(double & omega2, double & z){
  return .5 * (1. + omega2 - z);
}
double topwidth::CE_Pplus(double & omega2, double & z){
  return .5 * (1. - omega2 + z) + .5 * sqrt(lambda(1., omega2, z));
}
double topwidth::CE_Pminus(double & omega2, double & z){
  return .5 * (1. - omega2 + z) - .5 * sqrt(lambda(1., omega2, z));
}
double topwidth::CE_Yp(double & omega2, double & z){
  return .5 * log((.5 * (1. - omega2 + z) + .5 * sqrt(lambda(1., omega2, z))) / (.5 * (1. - omega2 + z) - .5 * sqrt(lambda(1., omega2, z))));
}
double topwidth::CE_Wplus(double & omega2, double & z){
  return .5 * (1. + omega2 - z) + .5 * sqrt(lambda(1., omega2, z));
}
double topwidth::CE_Wminus(double & omega2, double & z){
  return .5 * (1. + omega2 - z) - .5 * sqrt(lambda(1., omega2, z));
}
double topwidth::CE_Yw(double & omega2, double & z){
  return .5 * log((.5 * (1. + omega2 - z) + .5 * sqrt(lambda(1., omega2, z))) / (.5 * (1. + omega2 - z) - .5 * sqrt(lambda(1., omega2, z))));
}

double topwidth::CE_Gamma_0(double & Gamma_t_infinity, double & omega2, double & beta2){
  return Gamma_t_infinity * 2 * CE_P3(omega2, beta2) * CE_f(omega2, beta2); 
}
double topwidth::CE_Gamma_1(double & Gamma_t_infinity, double & omega2, double & beta2){
  if (beta2 == 0.){
    return Gamma_t_infinity * C_F / (2 * pi) 
      * (4 * pow(1. - omega2, 2) * (1. + 2 * omega2) * (gsl_sf_dilog(1. - omega2) - pi2_3 + .5 * log(1. - omega2) * log(omega2))
	 - 2 * omega2 * (1. + omega2) * (1. - 2 * omega2) * log(omega2)
	 - pow(1. - omega2, 2) * (4 * omega2 + 5.) * log(1. - omega2)
	 + .5 * (1. - omega2) * (5. + 9 * omega2 - 6 * pow(omega2, 2))
	 );
  }
  else {
    return Gamma_t_infinity * C_F / (2 * pi) 
    * (8 * CE_f(omega2, beta2) * CE_P0(omega2, beta2) 
       * (gsl_sf_dilog(1. - CE_Pminus(omega2, beta2)) 
	  - gsl_sf_dilog(1. - CE_Pplus(omega2, beta2))
	  - 2 * gsl_sf_dilog(1. - CE_Pminus(omega2, beta2) / CE_Pplus(omega2, beta2))
	  + CE_Yp(omega2, beta2) * log((4 * pow(CE_P3(omega2, beta2), 2)) / (pow(CE_Pplus(omega2, beta2), 2) * CE_Wplus(omega2, beta2)))
	  + CE_Yw(omega2, beta2) * log(CE_Pplus(omega2, beta2)))
       + 4 * (1 - beta2) * (pow(1 - beta2, 2) + omega2 * (1 + beta2) - 4 * pow(omega2, 2)) * CE_Yw(omega2, beta2)
       + (3. - beta2 + 11 * pow(beta2, 2) - pow(beta2, 3) + omega2 * (6. - 12 * beta2 + 2 * pow(beta2, 2)) - pow(omega2, 2) * (21. + 5 * beta2) + 12 * pow(omega2, 3)) * CE_Yp(omega2, beta2) 
       + 8 * CE_f(omega2, beta2) * CE_P3(omega2, beta2) * log(sqrt(omega2) / (4 * pow(CE_P3(omega2, beta2), 2))) 
       + 6 * (1. - 4 * beta2 + 3 * pow(beta2, 2) + omega2 * (3. + beta2) - 4 * pow(omega2, 2)) * CE_P3(omega2, beta2) * log(sqrt(beta2))
       + (5. - 22 * beta2 + 5 * pow(beta2, 2) + 9 * omega2 * (1. + beta2) - 6 * pow(omega2, 2)) * CE_P3(omega2, beta2)
       );
  }
}



double topwidth::c_propagator_Breit_Wigner(double & r, int m, double smin, double smax){
  Logger logger("c_propagator_Breit_Wigner");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // see Bronstein, Sec. 2.8.7, Eq. (2.7)
  double arg1=(smin - msi->M2[m]) / (msi->M[m] * msi->Gamma[m]);
  double arg2=(smax - msi->M2[m]) / (msi->M[m] * msi->Gamma[m]);
  double delta_y;
  if (arg1*arg2>-1)
    delta_y=atan((arg2-arg1)/(1+arg1*arg2));
  else if (arg2>0)
    delta_y=M_PI+atan((arg2-arg1)/(1+arg1*arg2));
  else
    delta_y=-M_PI+atan((arg2-arg1)/(1+arg1*arg2));

  // see Bronstein, Sec. 2.7.2.3, Eq. (2.87)
  double s = msi->M[m] * msi->Gamma[m] * (arg1+tan(r*delta_y))/(1-arg1*tan(r*delta_y)) + msi->M2[m];

  logger << LOG_DEBUG_VERBOSE << "finished   s = " << s << endl;
  return s;
}



double topwidth::g_propagator_Breit_Wigner(double & s, int m, double smin, double smax){
  Logger logger("g_propagator_Breit_Wigner");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // see Bronstein, Sec. 2.8.7, Eq. (2.7)
  double arg1=(smin - msi->M2[m]) / (msi->M[m] * msi->Gamma[m]);
  double arg2=(smax - msi->M2[m]) / (msi->M[m] * msi->Gamma[m]);
  double delta_y;
  if (arg1*arg2>-1)
    delta_y=atan((arg2-arg1)/(1+arg1*arg2));
  else if (arg2>0)
    delta_y=M_PI+atan((arg2-arg1)/(1+arg1*arg2));
  else
    delta_y=-M_PI+atan((arg2-arg1)/(1+arg1*arg2));

  double gs = (msi->M[m] * msi->Gamma[m]) / (delta_y * (pow(s - msi->M2[m], 2) + msi->M2[m] * pow(msi->Gamma[m], 2)));

  logger << LOG_DEBUG_VERBOSE << "finished   gs = " << gs << endl;
  return gs;
}



double topwidth::integrate_topwidth(double & alpha_S_t, int order_alpha_S, double relative_accuracy){
  Logger logger("integrate_topwidth");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  //  cout << "Gamma_t = " << Gamma_t << endl;
  cout << "msi->M_W  = " << msi->M_W << endl;
  cout << "msi->M_t  = " << msi->M_t << endl;
  cout << "msi->M_b  = " << msi->M_b << endl;
  cout << "G_F       = " << msi->G_F << endl;
  cout << "Gamma_W   = " << msi->Gamma_W << endl;
  cout << "alpha_S_t = " << alpha_S_t << endl;
  cout << "msi->M2_W = " << msi->M2_W << endl;
  cout << "msi->M2_t = " << msi->M2_t << endl;
  cout << "msi->M2_b = " << msi->M2_b << endl;

  //  double CE_beta = msi->M_b / msi->M_t;
  double CE_beta2 = msi->M2_b / msi->M2_t;
  cout << "CE_beta2 = " << CE_beta2 << endl;
  //  double CE_xi = M2_t / M2_W;
  //  double CE_gamma = Gamma_W / msi->M_W;
  double CE_omega2;
  //  double CE_z;

  double Gamma_t_infinity = msi->G_F * pow(msi->M_t, 3) / (8 * pi * sqrt(2.));
  cout << "Gamma_t_infinity = " << Gamma_t_infinity << endl;

  CE_omega2 = msi->M2_W / msi->M2_t;

  cout << "CE_omega2 = " << CE_omega2 << endl;

  double Gamma_0_NW = CE_Gamma_0(Gamma_t_infinity, CE_omega2, CE_beta2);

  cout << "Gamma_t_0(M2_W) = " << Gamma_0_NW << endl;

  double Gamma_1_NW = CE_Gamma_1(Gamma_t_infinity, CE_omega2, CE_beta2);
  cout << "Gamma_t_1(M2_W) = " << Gamma_1_NW << endl;
  cout << "Gamma_t_1(M2_W) / Gamma_t_0(M2_W)= " << Gamma_1_NW/ Gamma_0_NW << endl;
  cout << endl;
  cout << "Gamma_t_LO,NW = " << Gamma_0_NW << endl;
  cout << "Gamma_t_NLO,NW = " << Gamma_0_NW + alpha_S_t * Gamma_1_NW << endl;

  cout << endl;
  double Gamma_0_BW = 0.;
  double Gamma_1_BW = 0.;
  double integrand = 0.;
  double integrand_0_BW = 0.;
  double integrand_1_BW = 0.;
  double smax = pow(msi->M_t - msi->M_b, 2);
  double p2W = 0.;
  double g_p2W = 0.;

  double sum_int_Gamma = 0.;
  double sum2_int_Gamma = 0.;
  double sum_int_Gamma_0_BW = 0.;
  double sum2_int_Gamma_0_BW = 0.;
  double sum_int_Gamma_1_BW = 0.;
  double sum2_int_Gamma_1_BW = 0.;

  double result_Gamma_t = 0.;
  double result_Gamma_t_delta = 0.;
  double result_Gamma_0_BW = 0.;
  double result_Gamma_0_BW_delta = 0.;
  //  double result_Gamma_1_BW;
  //  double result_Gamma_1_BW_delta;

  long long i = 0;
  cout << "smax = " << smax << endl;
  //  cout << "psi_n_step = " << psi_n_step << endl;
  long long n_step = 100000;
  long long n_events_max = 1000000000;

  cout << "n_events_max = " << n_events_max << endl;
  double dn_events_max;
  if (n_events_max < 100){
    dn_events_max = pow(10., n_events_max);
  }
  else {
    dn_events_max = n_events_max;
  }
  cout << "dn_events_max = " << dn_events_max << endl;
  cout << endl;
  //  no_random = 0;
  //  vector<double> r(no_random + 1);

  int no_random = 0;
  vector<double> random;
  random.resize(no_random + 1);
  vector<double> sran(3);
  sran[0] = sqrt(2.) - 1.;
  sran[1] = sqrt(3.) - 1.;
  sran[2] = sqrt(5.) - 2.;

  long long ll_point = 0;
  double di = 0.;
  int int_end = 0;
  while (int_end == 0){
    i++;
    di++;
    ll_point++;
    logger << LOG_DEBUG << "i = " << i << "   di = " << di << endl;

    logger << LOG_DEBUG << "random.size() = " << random.size() << endl;
    logger << LOG_DEBUG << "sran.size() = " << sran.size() << endl;
    randomvector(sran, 0, random);
    logger << LOG_DEBUG << "random.size() = " << random.size() << endl;
    logger << LOG_DEBUG << "sran.size() = " << sran.size() << endl;
    //    cout << "random[0] = " << random[0] << endl;

    //double topwidth::c_propagator_Breit_Wigner(double & r, int m, double & smin, double & smax)
    logger << LOG_DEBUG << "random[0] = " << random[0] << endl;
    logger << LOG_DEBUG << "msi->Gamma[24] = " << msi->Gamma[24] << endl;

    p2W = c_propagator_Breit_Wigner(random[0], 24, 0., smax);
//    p2W = h_propagator(random[0], -24, 0., smax, psi);
    logger << LOG_DEBUG << "p2W = " << p2W << endl;

    CE_omega2 = p2W / msi->M2_t;
    //    cout << "CE_omega2 = " << CE_omega2 << endl;

    g_p2W = g_propagator_Breit_Wigner(p2W, 24, 0., smax);
    logger << LOG_DEBUG << "g_p2W = " << g_p2W << endl;
    //    g_p2W = gs_propagator(p2W, -24, 0., smax, psi);
    //    cout << "g_p2W = " << g_p2W << endl;
    Gamma_0_BW = CE_Gamma_0(Gamma_t_infinity, CE_omega2, CE_beta2);
    Gamma_1_BW = CE_Gamma_1(Gamma_t_infinity, CE_omega2, CE_beta2);
    
    integrand = (msi->M_W * msi->Gamma_W / pi) / (pow(p2W - msi->M2_W, 2) + msi->M2_W * pow(msi->Gamma_W, 2)) * (Gamma_0_BW + alpha_S_t * Gamma_1_BW);
    integrand_0_BW = (msi->M_W * msi->Gamma_W / pi) / (pow(p2W - msi->M2_W, 2) + msi->M2_W * pow(msi->Gamma_W, 2)) * Gamma_0_BW;
    integrand_1_BW = (msi->M_W * msi->Gamma_W / pi) / (pow(p2W - msi->M2_W, 2) + msi->M2_W * pow(msi->Gamma_W, 2)) * Gamma_1_BW;
    //    cout << "integrand = " << integrand << endl;

    sum_int_Gamma += integrand / g_p2W;
    sum2_int_Gamma += pow(integrand / g_p2W, 2);

    sum_int_Gamma_0_BW += integrand_0_BW / g_p2W;
    sum2_int_Gamma_0_BW += pow(integrand_0_BW / g_p2W, 2);

    sum_int_Gamma_1_BW += integrand_1_BW / g_p2W;
    sum2_int_Gamma_1_BW += pow(integrand_1_BW / g_p2W, 2);

    //    if (i > 0 && (i % psi_n_step) == 0){
    if (ll_point == n_step){
      ll_point = 0;
      result_Gamma_t = sum_int_Gamma / di;
      result_Gamma_t_delta = sqrt((sum2_int_Gamma - pow(sum_int_Gamma, 2) / di)) / (di - 1.);

      result_Gamma_0_BW = sum_int_Gamma_0_BW / i;
      result_Gamma_0_BW_delta = sqrt((sum2_int_Gamma_0_BW - pow(sum_int_Gamma_0_BW, 2) / di)) / (di - 1.);

      //      result_Gamma_1_BW = sum_int_Gamma_1_BW / i;
      //      result_Gamma_1_BW_delta = sqrt((sum2_int_Gamma_1_BW - pow(sum_int_Gamma_1_BW, 2) / di)) / (di - 1.);
      cout << showpoint;
      int n_points_log = int(log10(di) + 1);
      cout << setw(16) << setprecision(n_points_log) << noshowpoint << di << "   Gamma_t_LO  = " << setw(12) << setprecision(8) << showpoint << result_Gamma_0_BW << " +/- " << setw(12) << setprecision(8) << result_Gamma_0_BW_delta << endl;
      cout << setw(16) << setprecision(n_points_log) << noshowpoint << di << "   Gamma_t_NLO = " << setw(12) << setprecision(8) << showpoint << result_Gamma_t << " +/- " << setw(12) << setprecision(8) << result_Gamma_t_delta << endl;
      if (order_alpha_S == 0 && result_Gamma_0_BW_delta / result_Gamma_0_BW < relative_accuracy){int_end = 1;}
      else if (order_alpha_S == 1 && result_Gamma_t_delta / result_Gamma_t < relative_accuracy){int_end = 1;}
    }

    if (di >= dn_events_max){int_end = 1;}
  }

  
  if (order_alpha_S == 0){
    logger << LOG_INFO << "Gamma_tWb = " << result_Gamma_0_BW << endl;
    logger << LOG_DEBUG_VERBOSE << "finished" << endl;
    return result_Gamma_0_BW;
  }
  else if (order_alpha_S == 1){
    logger << LOG_INFO << "Gamma_tWb = " << result_Gamma_t<< endl;
    logger << LOG_DEBUG_VERBOSE << "finished" << endl;
    return result_Gamma_t;
  }
  else {logger << LOG_FATAL << "Width could not be calculated at the requested order." << endl; return 0.;}
}

