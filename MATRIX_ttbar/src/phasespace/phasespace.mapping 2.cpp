#include "../include/classes.cxx"

/////////////////////
//  angle mapping  //
/////////////////////

double phasespace_set::c_phi(int i_r){
  static Logger logger("phasespace_set::c_phi");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double phi = 2 * pi * r[i_r];

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << phi << endl;
  return phi;
}

double phasespace_set::inv_phi(double & phi){
  static Logger logger("phasespace_set::inv_phi");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double r_phi = phi / (2 * pi);

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << r_phi << endl;
  return r_phi;
}

double phasespace_set::c_costheta(int i_r){
  static Logger logger("phasespace_set::c_costheta");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double costheta = 2 * r[i_r] - 1.;

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << costheta << endl;
  return costheta;
}

double phasespace_set::inv_costheta(double & costheta){
  static Logger logger("phasespace_set::inv_costheta");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double r_costheta = (costheta + 1.) / 2;

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << r_costheta << endl;
  return r_costheta;
}



//////////////////////////
//  propagator mapping  //
//////////////////////////

void phasespace_set::c_propagator(int i_a, int m, int out, double smin, double smax, int i_r1){
  static Logger logger("phasespace_set::c_propagator");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG_VERBOSE << "s[" << out << "] with m = " << m << "    map_Gamma = " << map_Gamma[abs(m)] << "    M2 = " << M2[abs(m)] << "   smin = " << smin << "   smax = " << smax << endl;

  // Just for compatibility; could be removed now, I guess !!!
  m = abs(m);

  if (map_Gamma[m] == 0){
    double mass2;
    if (M2[m] == 0.){mass2 = mass0;}
    else {mass2 = M2[m] + mass0;}
    xbs_all[i_a][out] = c_propagator_vanishing_width(r[i_r1], mass2, smin, smax, nuxs);
  }
  else {
    xbs_all[i_a][out] = c_propagator_Breit_Wigner(r[i_r1], m, smin, smax);
    //    xbs_all[i_a][out] = c_propagator_Breit_Wigner(r[i_r1], m, smin, smax, *this);
    //    xbs_all[i_a][out] = c_propagator_Breit_Wigner(r[i_r1], m, smin, smax, psi);
  }
  xbsqrts_all[i_a][out] = sqrt(xbs_all[i_a][out]);

  if (munich_isnan(xbsqrts_all[i_a][out])) {
    logger << LOG_DEBUG << "xbsqrts_all[i_a][out] is NaN (r[" << i_r1 << "]=" << r[i_r1] << ", m=" << m << ", out=" << out << "smin=" << smin << ", smax=" << smax << ", map_Gamma[m]=" << map_Gamma[m] << ", M2[m]=" << M2[m] << ", map_Gamma[m]=" << map_Gamma[m] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

double phasespace_set::g_propagator(int i_a, int m, int out, double smin, double smax){
  Logger logger("phasespace_set::g_propagator");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG_VERBOSE << "s[" << out << "] with m = " << m << "    map_Gamma = " << map_Gamma[abs(m)] << "    M2 = " << M2[abs(m)] << endl;

  // !!! What should happen in this case ???
  if (switch_console_output_phasespace_issue){
    if (smax <= smin){
      logger << LOG_ERROR << "smin = " << smin << "   smax = " << smax << ", m = " << m << ", out = " << out << endl;
      for (int i = 0; i < xbs_all[i_a].size(); i++){
	if (xbs_all[i_a][i] != 0.){
	  logger << LOG_ERROR << "xbs_all[" << i_a << "][" << i << "] = " << xbs_all[i_a][i] << endl;
	}
      }
    }
  }

  double gs;

  // Just for compatibility; could be removed now, I guess !!!
  m = abs(m);

  if (map_Gamma[m] == 0){
    double mass2;
    if (M2[m] == 0.){mass2 = mass0;}
    else {mass2 = M2[m] + mass0;}
    gs = g_propagator_vanishing_width(xbs_all[i_a][out], mass2, smin, smax, nuxs);
  }
  else {
    gs = g_propagator_Breit_Wigner(xbs_all[i_a][out], m, smin, smax);
    //    gs = g_propagator_Breit_Wigner(xbs_all[i_a][out], m, smin, smax, *this);
    //    gs = g_propagator_Breit_Wigner(xbs_all[i_a][out], m, smin, smax, psi);
  }
  //  if (gs != gs){cout << "smin = " << smin << "   s[" << out << "] = " << xbs_all[i_a][out] << "   smax = " << smax << "   gs = " << gs << "   gtot = nan from s-channel mapping" << endl;}
   //  if (gs != gs){cout << "smin = " << smin << "   smax = " << smax << "   gtot = nan from s-channel mapping" << endl;}

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << gs << endl;
  if (gs < 0.){logger << LOG_DEBUG << "gs < 0.   smin = " << smin << "   s[" << out << "] = " << xbs_all[i_a][out] << "   smax = " << smax << "   gs = " << gs << "   M[" << m << "] = " << M[m] << "   map_Gamma[" << m << "] = " << map_Gamma[m] << endl;}
  return gs;
}

void phasespace_set::inv_propagator(int i_a, int m, int out, double smin, double smax, double & r1){
  static Logger logger("phasespace_set::inv_propagator");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // Just for compatibility; could be removed now, I guess !!!
  m = abs(m);

  if (map_Gamma[m] == 0){
    double mass2;
    if (M2[m] == 0.){mass2 = mass0;}
    else {mass2 = M2[m] + mass0;}
    r1 = inv_propagator_vanishing_width(xbs_all[i_a][out], mass2, smin, smax, nuxs);
  }
  else {
    r1 = inv_propagator_Breit_Wigner(xbs_all[i_a][out], m, smin, smax);
    //    r1 = inv_propagator_Breit_Wigner(xbs_all[i_a][out], m, smin, smax, *this);
    //    r1 = inv_propagator_Breit_Wigner(xbs_all[i_a][out], m, smin, smax, psi);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



///////////////////////////////////////
//  Breit-Wigner propagator mapping  //
///////////////////////////////////////

double phasespace_set::c_propagator_Breit_Wigner(double & r1, int m, double & smin, double & smax){
  Logger logger("phasespace_set::c_propagator_Breit_Wigner");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // double y1 = atan((smin - M2[m]) / (M[m] * map_Gamma[m]));
  // double y2 = atan((smax - M2[m]) / (M[m] * map_Gamma[m]));
  // see Bronstein, Sec. 2.8.7, Eq. (2.7)
  double arg1 = (smin - M2[m]) / (M[m] * map_Gamma[m]);
  double arg2 = (smax - M2[m]) / (M[m] * map_Gamma[m]);
  double delta_y;
  if (arg1 * arg2 > -1){delta_y = atan((arg2 - arg1) / (1 + arg1 * arg2));}
  else if (arg2 > 0){delta_y = M_PI + atan((arg2 - arg1) / (1. + arg1 * arg2));}
  else {delta_y = -M_PI + atan((arg2 - arg1) / (1 + arg1 * arg2));}
  // see Bronstein, Sec. 2.7.2.3, Eq. (2.87)
  double s = M[m] * map_Gamma[m] * (arg1 + tan(r1 * delta_y)) / (1 - arg1 * tan(r1 * delta_y)) + M2[m];

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << s << endl;
  return s;
}

double phasespace_set::g_propagator_Breit_Wigner(double & s, int m, double & smin, double & smax){
  Logger logger("phasespace_set::g_propagator_Breit_Wigner");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // double y1 = atan((smin - M2[m]) / (M[m] * map_Gamma[m]));
  // double y2 = atan((smax - M2[m]) / (M[m] * map_Gamma[m]));
  // see Bronstein, Sec. 2.8.7, Eq. (2.7)
  double arg1 = (smin - M2[m]) / (M[m] * map_Gamma[m]);
  double arg2 = (smax - M2[m]) / (M[m] * map_Gamma[m]);
  double delta_y;
  if (arg1 * arg2 > -1){delta_y = atan((arg2 - arg1)/(1 + arg1 * arg2));}
  else if (arg2 > 0){delta_y = M_PI + atan((arg2 - arg1) / (1 + arg1 * arg2));}
  else {delta_y = -M_PI + atan((arg2 - arg1) / (1 + arg1 * arg2));}
  // double gs = (M[m] * map_Gamma[m]) / ((y2 - y1) * (pow(s - M2[m], 2) + M2[m] * pow(map_Gamma[m], 2)));
  double gs = (M[m] * map_Gamma[m]) / (delta_y * (pow(s - M2[m], 2) + M2[m] * pow(map_Gamma[m], 2)));

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << gs << endl;
  if (gs < 0.){logger << LOG_DEBUG << "gs < 0.   smin = " << smin << "   s = " << s << "   smax = " << smax << "   gs = " << gs << "   M[" << m << "] = " << M[m] << "   map_Gamma[" << m << "] = " << map_Gamma[m] << endl;}
  return gs;
}

double phasespace_set::inv_propagator_Breit_Wigner(double s, int m, double & smin, double & smax){
  Logger logger("phasespace_set::inv_propagator_Breit_Wigner");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // TODO: micro optimize as above
  double arg1 = atan((smin - M2[m]) / (M[m] * map_Gamma[m]));
  double arg2 = atan((smax - M2[m]) / (M[m] * map_Gamma[m]));

  double r = (atan((s - M2[m]) / (M[m] * map_Gamma[m])) - arg1) / (arg2 - arg1);
  if (munich_isnan(r)) {
    logger << LOG_DEBUG << "r is NaN in inv_propagator_Breit_Wigner!" << endl;
    logger << LOG_DEBUG << "    r = " << setw(23) << setprecision(15) << r << ";            smin = " << setw(23) << setprecision(15) << smin << ";      s = " << setw(23) << setprecision(15) << s << ";   smax = " << setw(23) << setprecision(15) << smax << endl;
    logger << LOG_DEBUG << "M[" << setw(2) << m << "] = " << setw(23) << setprecision(15) << M[m] << ";   map_Gamma[" << setw(2) << m << "] = " << setw(23) << setprecision(15) << map_Gamma[m] << ";   arg1 = " << setw(23) << setprecision(15) << arg1 << ";   arg2 = " << setw(23) << setprecision(15) << arg2 << endl;
    if (s <= smin){r = map_technical_s;}
    else if (s >= smax){r = 1. - map_technical_s;}
    else {
      if (switch_console_output_phasespace_issue){
        logger << LOG_ERROR << "r remains NaN in inv_propagator_Breit_Wigner!" << endl;
        logger << LOG_ERROR << "   smin = " << smin << "   s = " << s << "   smax = " << smax << endl;
	logger << LOG_ERROR << "   arg1 = " << arg1 << "   arg2 = " << arg2 << endl;
	logger << LOG_ERROR << "   M[" << m << "] = " << M[m] << "   M2[" << m << "] = " << M2[m] << "   map_Gamma[" << m << "] = " << map_Gamma[m] << endl;
      }
      r = map_technical_s;
      // Otherwise program might stop !!! Reason for "nan" should be clarified !!!
    }
  }
  else if (r <= 0.) {
    logger << LOG_DEBUG << "r <= 0 in inv_propagator_Breit_Wigner!" << endl;
    logger << LOG_DEBUG << "    r = " << setw(23) << setprecision(15) << r << ";            smin = " << setw(23) << setprecision(15) << smin << ";      s = " << setw(23) << setprecision(15) << s << ";   smax = " << setw(23) << setprecision(15) << smax << endl;
    logger << LOG_DEBUG << "M[" << setw(2) << m << "] = " << setw(23) << setprecision(15) << M[m] << ";   map_Gamma[" << setw(2) << m << "] = " << setw(23) << setprecision(15) << map_Gamma[m] << ";   arg1 = " << setw(23) << setprecision(15) << arg1 << ";   arg2 = " << setw(23) << setprecision(15) << arg2 << endl;
    r = map_technical_s;
  }
  else if (r >= 1.) {
    logger << LOG_DEBUG << "r >= 1 in inv_propagator_Breit_Wigner!" << endl;
    logger << LOG_DEBUG << "    r = " << setw(23) << setprecision(15) << r << ";            smin = " << setw(23) << setprecision(15) << smin << ";      s = " << setw(23) << setprecision(15) << s << ";   smax = " << setw(23) << setprecision(15) << smax << endl;
    logger << LOG_DEBUG << "M[" << setw(2) << m << "] = " << setw(23) << setprecision(15) << M[m] << ";   map_Gamma[" << setw(2) << m << "] = " << setw(23) << setprecision(15) << map_Gamma[m] << ";   arg1 = " << setw(23) << setprecision(15) << arg1 << ";   arg2 = " << setw(23) << setprecision(15) << arg2 << endl;
    r = 1. - map_technical_s;
  }

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << r << endl;
  return r;
}



///////////////////////////////////////
//  narrow-width propagator mapping  //
///////////////////////////////////////

double phasespace_set::c_propagator_narrow_width(){return 0.;}

double phasespace_set::g_propagator_narrow_width(int m){
  Logger logger("phasespace_set::g_propagator_narrow_width");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double gs = 1. / (pi * M[m] * map_Gamma[m]);

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << gs << endl;
  if (gs < 0.){logger << LOG_DEBUG << "gs < 0.   gs = " << gs << "   M[" << m << "] = " << M[m] << "   map_Gamma[" << m << "] = " << map_Gamma[m] << endl;}
  return gs;
}

double phasespace_set::inv_propagator_narrow_width(){return 0.;}



//////////////////////////////////////////
//  vanishing-width propagator mapping  //
//////////////////////////////////////////

double phasespace_set::c_propagator_vanishing_width(double r, double mass2, double smin, double smax, double nu){
  Logger logger("c_propagator_vanishing_width");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double s;
  if (nu != 1.){s = pow(r * pow(smax - mass2, 1. - nu) + (1. - r) * pow(smin - mass2, 1. - nu), (1. / (1. - nu))) + mass2; }
  else {s = exp(r * log(smax - mass2) + (1. - r) * log(smin - mass2)) + mass2;}
  logger << LOG_DEBUG_VERBOSE << "s = " << s << endl;
  //  if (s != s){exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << s << endl;
  return s;
}

double phasespace_set::g_propagator_vanishing_width(double s, double mass2, double smin, double smax, double nu){
  Logger logger("g_propagator_vanishing_width");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG_VERBOSE << "smin = " << smin << "  s = " << s << "   smax = " << smax << "   mass2 = " << mass2 << "   smin - mass0 = " << smin - mass2 << " ^(1-nu = " << nu << ") = " << pow(smin - mass2, 1. - nu) << endl;
  double gs;
  if (smax <= smin) {
    if (switch_console_output_phasespace_issue){
      logger << LOG_ERROR << "smin = " << smin << "  s = " << s << "   smax = " << smax << "   mass2 = " << mass2 << "   smin - mass0 = " << smin-mass2 << " ^(1-nu = " << nu << ") = " << pow(smin-mass2,1-nu) << endl;
      logger << LOG_ERROR << "gs is set to infinity." << endl;
    }
    gs = 1. / 0.; //1.e99;
  }
  else {
    if (nu != 1.){gs = (1. - nu) / ((pow(smax - mass2, 1. - nu) - pow(smin - mass2, 1. - nu)) * pow(s - mass2, nu));}
    else {gs = 1. / ((log(smax - mass2) - log(smin - mass2)) * (s - mass2));}
  }
  //  if (gs != gs){exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << gs << endl;
  if (gs < 0.){logger << LOG_DEBUG << "gs < 0.   smin = " << smin << "   s = " << s << "   smax = " << smax << "   gs = " << gs << "   mass2 = " << mass2 << endl;}
  return gs;
}

double phasespace_set::inv_propagator_vanishing_width(double s, double mass2, double smin, double smax, double nu){
  Logger logger("inv_propagator_vanishing_width");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double r;
  if (nu != 1.) {
    double sminm0 = pow(smin-mass2,1-nu);
    r = (pow(s - mass2, 1. - nu) - sminm0) / (pow(smax - mass2, 1. - nu) - sminm0);
  } 
  else {r = log((s - mass2) / (smin - mass2)) / log((smax - mass2) / (smin - mass2));}

  static double temp_map_technical_r = 1.e-8;
  if (munich_isnan(r)) {
    logger << LOG_DEBUG << "smin = " << smin << endl;
    logger << LOG_DEBUG << "sminm0 = " << pow(smin-mass2,1-nu) << endl;
    logger << LOG_DEBUG << "smax = " << smax << endl;
    logger << LOG_DEBUG << "mass2= " << mass2 << endl;
    logger << LOG_DEBUG << "s    = " << s << endl;
    logger << LOG_DEBUG << "r is NaN in inv_propagator_vanishing_width!" << endl;
    logger << LOG_DEBUG << "r=" << r <<"; smin, s, smax=" << smin << ", " << s << ", " << smax << endl;
    //    r=0.5;
    if (s <= smin){r = temp_map_technical_r;}
    else if (s >= smax){r = 1. - temp_map_technical_r;}
    else {
      if (switch_console_output_phasespace_issue){
        logger << LOG_ERROR << "r remains NaN in inv_propagator_vanishing_width!" << endl;
      }
      r = temp_map_technical_r;
      // Otherwise program might stop !!! Reason for "nan" should be clarified !!!
    }
  }
  else if (r <= 0.){
    logger << LOG_DEBUG << "r<0 in inv_propagator_vanishing_width!" << endl;
    logger << LOG_DEBUG << "r=" << r <<"; smin, s, smax=" << smin << ", " << s << ", " << smax << endl;
    r = temp_map_technical_r;
  }
  if (r >= 1.){
    logger << LOG_DEBUG << "r>1 in inv_propagator_vanishing_width!" << endl;
    logger << LOG_DEBUG << "r=" << r <<"; smin, s, smax=" << smin << ", " << s << ", " << smax << endl;
    r = 1. - temp_map_technical_r;
  }

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << r << endl;
  return r;
}



////////////////////////////////////////
//  no-propagator-equivalent mapping  //
////////////////////////////////////////

void phasespace_set::c_timelikeinvariant(int i_a, int out, double smin, double smax, int i_r1){
  static Logger logger("phasespace_set::c_timelikeinvariant");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  xbs_all[i_a][out] = r[i_r1] * smax + (1. - r[i_r1]) * smin;
  xbsqrts_all[i_a][out] = sqrt(xbs_all[i_a][out]);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

double phasespace_set::g_timelikeinvariant(double smin, double smax){
  static Logger logger("phasespace_set::g_timelikeinvariant");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double gs = 1. / (smax - smin);

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << gs << endl;
  if (gs < 0.){logger << LOG_DEBUG << "gs < 0.   smin = " << smin << "   smax = " << smax << "   gs = " << gs << endl;}
  return gs;
}

void phasespace_set::inv_timelikeinvariant(int i_a, int out, double smin, double smax, double & r1){
  static Logger logger("inv_timelikeinvariant");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  r1 = (xbs_all[i_a][out] - smin) / (smax - smin);

  static double temp_map_technical_r = 1.e-8;
  if (munich_isnan(r1)) {
    logger << LOG_DEBUG << "r1 is NaN in inv_timelikeinvariant!" << endl;
    logger << LOG_DEBUG << "r1 = " << r1 <<"; smin, s, smax=" << smin << ", " << xbs_all[i_a][out] << ", " << smax << endl;
    //    r=0.5;
    if (xbs_all[i_a][out] <= smin){r1 = temp_map_technical_r;}
    else if (xbs_all[i_a][out] >= smax){r1 = 1. - temp_map_technical_r;}
    else {
      if (switch_console_output_phasespace_issue){
        logger << LOG_ERROR << "r1 remains NaN in inv_timelikeinvariant!" << endl;
      }
      r1 = temp_map_technical_r;
      // Otherwise program might stop !!! Reason for "nan" should be clarified !!!
    }
  }
  else if (r1 <= 0.){
    logger << LOG_DEBUG << "r1 < 0 in inv_timelikeinvariant!" << endl;
    logger << LOG_DEBUG << "r1 = " << r1 <<"; smin, s, smax=" << smin << ", " << xbs_all[i_a][out] << ", " << smin << endl;
    r1 = temp_map_technical_r;
  }
  if (r1 >= 1.){
    logger << LOG_DEBUG << "r1 > 1 in inv_timelikeinvariant!" << endl;
    logger << LOG_DEBUG << "r1 = " << r1 <<"; smin, s, smax=" << smin << ", " << xbs_all[i_a][out] << ", " << smin << endl;
    r1 = 1. - temp_map_technical_r;
  }
}



////////////////////////////////////
//  t-channel propagator mapping  //
////////////////////////////////////

double phasespace_set::c_t_propagator(double r, int m, double smin, double smax){
  static Logger logger("c_t_propagator");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double s, mass2, nuexp;
  if (M2[m] == 0.){
    mass2 = mass0;
    nuexp = nuxt;
  }
  else {
    mass2 = -M2[m];// + mass0;
    nuexp = nuxt;
  }
  if (nuexp != 1.){s = pow(r * pow(smax - mass2, 1 - nuexp) + (1 - r) * pow(smin - mass2, 1. - nuexp), (1. / (1. - nuexp))) + mass2;}
  else {s = exp(r * log(smax - mass2) + (1. - r) * log(smin - mass2)) + mass2;}
  //  if (s != s){cout << "tmin = " << -smin << "   tmax = " << -smax << "   p = nan from t-channel mapping" << endl;}

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << s << endl;
  if (s < 0.){logger << LOG_DEBUG << "s < 0.   smin = " << smin << "   s = " << s << "   smax = " << smax << endl;}

  return s;
}
double phasespace_set::g_t_propagator(double s, int m, double smin, double smax){
  static Logger logger("g_t_propagator");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double gs, mass2, nuexp;
  if (M2[m] == 0.){
    mass2 = mass0;
    nuexp = nuxt;
  }
  else {
    mass2 = -M2[m];// + mass0;
    nuexp = nuxt;
  }
  if (nuexp != 1.){gs = (1. - nuexp) / ((pow(smax - mass2, 1. - nuexp) - pow(smin - mass2, 1. - nuexp)) * pow(s - mass2, nuexp));}
  else {gs = 1. / ((log(smax - mass2) - log(smin - mass2)) * (s - mass2));}
  //  if (gs != gs){cout << "tmin = " << -smin << "   tmax = " << -smax << "   gtot = nan from t-channel mapping" << endl;}

  if (munich_isnan(gs)) {
    logger << LOG_DEBUG << "gs == nan in g_t_propagator!   -tmax = " << setprecision(15) << setw(23) << smin << "   -t = " << setprecision(15) << setw(23) << s << "   -tmin = " << setprecision(15) << setw(23) << smax << endl;
  }


  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << gs << endl;
  if (gs > 0.){logger << LOG_DEBUG << "gs > 0.   smin = " << smin << "   s = " << s << "   smax = " << smax << "   gs = " << gs << "   gtot = nan from s-channel mapping" << endl;}
  return gs;
}

double phasespace_set::inv_t_propagator(double s, int m, double smin, double smax){
  static Logger logger("inv_t_propagator");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double r, mass2, nuexp;
  if (M2[m] == 0.){
    mass2 = mass0;
    nuexp = nuxt;
  }
  else {
    mass2 = -M2[m];// + mass0;
    nuexp = nuxt;
  }
  if (nuexp != 1.) {
    double sminm0 = pow(smin - mass2, 1. - nuexp);
    r = (pow(s - mass2, 1. - nuexp) - sminm0) / (pow(smax - mass2, 1. - nuexp) - sminm0);
  } 
  else {r = log((s - mass2) / (smin - mass2)) / log((smax - mass2) / (smin - mass2));}

  static double temp_map_technical_r = 1.e-8;
  if (munich_isnan(r)) {
    logger << LOG_DEBUG << "r is NaN in inv_t_propagator!  r = " << setprecision(15) << setw(23) << r << "   smin = " << setprecision(15) << setw(23) << smin << "   s = " << setprecision(15) << setw(23) << s << "   smax = " << setprecision(15) << setw(23) << smax << endl;
    //    r=0.5;
    if (s <= smin){r = temp_map_technical_r;}
    else if (s >= smax){r = 1. - temp_map_technical_r;}
    else {
      if (switch_console_output_phasespace_issue){
	//	r = 0.5;
	logger << LOG_ERROR << "r remains NaN in inv_t_propagator!" << endl;
      }
      r = temp_map_technical_r;
      // Otherwise program might stop !!! Reason for "nan" should be clarified !!!
    }
  }
  else if (r <= 0.){
    //    logger << LOG_DEBUG << "r<0 in inv_t_propagator!" << endl;
    logger << LOG_DEBUG << "r < 0 in inv_t_propagator!  r = " << setprecision(15) << setw(23) << r << "   -tmax = " << setprecision(15) << setw(23) << smin << "   -t = " << setprecision(15) << setw(23) << s << "   -tmin = " << setprecision(15) << setw(23) << smax << endl;
    r = temp_map_technical_r;
  }
  else if (r >= 1.){
    logger << LOG_DEBUG << "r > 1 in inv_t_propagator!  r = " << setprecision(15) << setw(23) << r << "   -tmax = " << setprecision(15) << setw(23) << smin << "   -t = " << setprecision(15) << setw(23) << s << "   -tmin = " << setprecision(15) << setw(23) << smax << endl;
    r = 1. - temp_map_technical_r;
  }

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << r << endl;
  return r;
}



/////////////////////
//  decay mapping  //
/////////////////////

void phasespace_set::c_decay(int i_a, int in, int out1, int out2, int i_r1, int i_r2){
  Logger logger("phasespace_set::c_decay");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  logger << LOG_DEBUG << "c_decay: This routines is used !!!" << endl; 

  double costheta = c_costheta(i_r1);
  //  double costheta = c_costheta(r[i_r1]);
  double phi = c_phi(i_r2);
  //  double phi = c_phi(r[i_r2]);
  fourvector k10(((xbs_all[i_a][in] + xbs_all[i_a][out1] - xbs_all[i_a][out2]) / (2 * xbsqrts_all[i_a][in])), 0., 0., (sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][out1], xbs_all[i_a][out2])) / (2 * xbsqrts_all[i_a][in])));
  xbp_all[i_a][out1] = k10.rot(phi, costheta);
  xbp_all[i_a][out1] = xbp_all[i_a][out1].boost(xbp_all[i_a][in].Pinv(), xbs_all[i_a][in]);
  xbp_all[i_a][out2] = xbp_all[i_a][in] - xbp_all[i_a][out1];

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

double phasespace_set::g_decay(int i_a, int in, int out1, int out2){
  Logger logger("phasespace_set::g_decay");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  logger << LOG_DEBUG << "g_decay: This routines is used !!!" << endl; 

  double gd = (2 * xbs_all[i_a][in]) / (pi * sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][out1], xbs_all[i_a][out2])));

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << gd << endl;
  if (gd < 0.){logger << LOG_DEBUG << "gd < 0.   in = " << in << "   out1 = " << out1 << "   out2 = " << out2 << "   xbs_all[" << i_a << "][" << in << "] = " << xbs_all[i_a][in] << endl;}
  return gd;
}

void phasespace_set::inv_decay(int i_a, int in, int out1, int out2, double & r1, double & r2){
  static Logger logger("phasespace_set::inv_decay");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  fourvector k10 = xbp_all[i_a][out1].boost(xbp_all[i_a][in], xbs_all[i_a][in]);

  //  logger << LOG_DEBUG << "inv_decay: This routines is used !!!" << endl; 

  double phi = 0.;
  double costheta = 0.;
  if (k10.x1() != 0 || k10.x2() != 0) {
    phi = atan2(k10.x2(), k10.x1());
    if (phi < 0.){phi += 2 * pi;}
  }
  else {
    logger << LOG_DEBUG << "x1=x2=0 in inv_decay" << endl;
    phi = atan2(1., 1.);
  }

  //  costheta=k10.x3()/(sqrt(lambda(s[in], s[out1], s[out2])) / (2 * sqrt(s[in])));
  costheta = k10.x3() / k10.r();

  r2 = inv_phi(phi);
  r1 = inv_costheta(costheta);

  static double temp_map_technical_r = 1.e-8;

  if (munich_isnan(r1)) {
    logger << LOG_DEBUG << "r1 is NaN in inv_decay!" << endl;
    logger << LOG_DEBUG << "r1 = " << r1 <<"; k10.x, k10.y = " << k10.x1() << ", " << k10.x2() << endl;
    r1 = 0.5; // why???
  }
  else if (r1 <= 0.){
    logger << LOG_DEBUG << "r1 < 0 in inv_decay!" << endl;
    logger << LOG_DEBUG << "r1 = " << r1 <<"; k10.x, k10.y = " << k10.x1() << ", " << k10.x2() << endl;
    r1 = temp_map_technical_r;
  }
  else if (r1 >= 1.){
    logger << LOG_DEBUG << "r1 > 1 in inv_decay!" << endl;
    logger << LOG_DEBUG << "r1 = " << r1 << "; k10.x, k10.y=" << k10.x1() << ", " << k10.x2() << endl;
    r1 = 1. - temp_map_technical_r;
  }

  if (munich_isnan(r2)){
    logger << LOG_DEBUG << "r2 is NaN in inv_decay!" << endl;
    logger << LOG_DEBUG << "r2 = " << r2 <<"; k10.z, r = " << k10.x3() << ", " << k10.r() << endl;
    logger << LOG_DEBUG << xbp_all[i_a][in] << endl;
    r2 = 0.5;
  }
  // important typo before:   !!! check if relevant (most likely not used yet) !!!
  //  if (r2 < 1){
  else if (r2 <= 0.){
    logger << LOG_DEBUG << "r2 < 0 in inv_decay!" << endl;
    logger << LOG_DEBUG << "r2 = " << r2 <<"; k10.z, r = " << k10.x3() << ", " << k10.r() << endl;
    logger << LOG_DEBUG << xbp_all[i_a][in] << endl;
    r2 = temp_map_technical_r;
  }
  else if (r2 >= 1.){
    logger << LOG_DEBUG << "r2 > 1 in inv_decay!" << endl;
    logger << LOG_DEBUG << "r2 = " << r2 <<"; k10.z, r=" << k10.x3() << ", " << k10.r() << endl;
    logger << LOG_DEBUG << xbp_all[i_a][in] << endl;
    r2 = 1. - temp_map_technical_r;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



/////////////////////////
//  t-channel mapping  //
/////////////////////////

void phasespace_set::c_tchannel(int i_a, int m, int in, int in1, int in2, int out1, int out2, int i_r1, int i_r2){
  static Logger logger("phasespace_set::c_tchannel");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double term1 = (xbs_all[i_a][in] + xbs_all[i_a][out1] - xbs_all[i_a][out2]) * (xbs_all[i_a][in] + xbs_all[i_a][in1] - xbs_all[i_a][in2]);
  double term2 = sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][out1], xbs_all[i_a][out2])) * sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][in1], xbs_all[i_a][in2]));
  double term3 = (2 * sqrt(xbs_all[i_a][in]));
  double tmin = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 + term2) / (2 * xbs_all[i_a][in]);
  double tmax = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 - term2) / (2 * xbs_all[i_a][in]);
  xbs_all[i_a][in1 + out1] = -c_t_propagator(r[i_r1], m, -tmin, -tmax);
  //xbs_all[i_a][in1 + out1] = -c_t_propagator(r[i_r1], m, -tmin, -tmax, *this);
  double phi = c_phi(i_r2);
  //  double phi = c_phi(r[i_r2]);
  double costheta = ((2 * xbs_all[i_a][in] * (xbs_all[i_a][in1 + out1] - xbs_all[i_a][out1] - xbs_all[i_a][in1])) + term1) / (term2);
  fourvector p10(((xbs_all[i_a][in] + xbs_all[i_a][in1] - xbs_all[i_a][in2]) / term3), 0, 0, (sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][in1], xbs_all[i_a][in2])) / term3));
  fourvector p20(((xbs_all[i_a][in] - xbs_all[i_a][in1] + xbs_all[i_a][in2]) / term3), 0, 0, -(sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][in1], xbs_all[i_a][in2])) / term3));
  fourvector k10(((xbs_all[i_a][in] + xbs_all[i_a][out1] - xbs_all[i_a][out2]) / term3), 0, 0, (sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][out1], xbs_all[i_a][out2])) / term3));
  k10 = k10.rot3(phi, costheta);
  fourvector k20 = p10 + p20 - k10;
  fourvector p1b = xbp_all[i_a][in1].boost(xbp_all[i_a][in], xbs_all[i_a][in]);
  xbp_all[i_a][out1] = (k10.rot(p1b)).boost(xbp_all[i_a][in].Pinv(), xbs_all[i_a][in]);
  xbp_all[i_a][out2] = (k20.rot(p1b)).boost(xbp_all[i_a][in].Pinv(), xbs_all[i_a][in]);
  xbp_all[i_a][in1 + out1] = xbp_all[i_a][in1] - xbp_all[i_a][out1];

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

double phasespace_set::phasespace_set::g_tchannel(int i_a, int m, int in, int in1, int in2, int out1, int out2){
  static Logger logger("phasespace_set::g_tchannel");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double term1 = (xbs_all[i_a][in] + xbs_all[i_a][out1] - xbs_all[i_a][out2]) * (xbs_all[i_a][in] + xbs_all[i_a][in1] - xbs_all[i_a][in2]);
  double term2 = sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][out1], xbs_all[i_a][out2])) * sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][in1], xbs_all[i_a][in2]));
  double tmin = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 + term2) / (2 * xbs_all[i_a][in]);
  double tmax = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 - term2) / (2 * xbs_all[i_a][in]);
  if (tmin > 0.){logger << LOG_DEBUG << "tmin = " << tmin << endl;}
  if (tmax > 0.){logger << LOG_DEBUG << "tmax = " << tmax << endl; tmax = 0.;}
  double gp = -g_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax) * (2 * sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][in1], xbs_all[i_a][in2]))) / pi;
  //  double gp = -g_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax, *this) * (2 * sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][in1], xbs_all[i_a][in2]))) / pi;

  logger << LOG_DEBUG_VERBOSE << "finished - return " << setw(23) << setprecision(15) << gp << endl;
  if (gp < 0.){logger << LOG_DEBUG << "gp < 0." << "   in = " << in << "   in1 = " << in1 << "   in2 = " << in2 << "   out1 = " << out1 << "   out2 = " << out2 << endl; gp = 0.;}
  return gp;
}

void phasespace_set::inv_tchannel(int i_a, int m, int in, int in1, int in2, int out1, int out2, double & r1, double & r2){
  static Logger logger("phasespace_set::inv_tchannel");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double term1 = (xbs_all[i_a][in] + xbs_all[i_a][out1] - xbs_all[i_a][out2]) * (xbs_all[i_a][in] + xbs_all[i_a][in1] - xbs_all[i_a][in2]);
  double term2 = sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][out1], xbs_all[i_a][out2])) * sqrt(lambda(xbs_all[i_a][in], xbs_all[i_a][in1], xbs_all[i_a][in2]));
  //  double term3 = (2 * sqrt(xbs_all[i_a][in]));
  double tmin = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 + term2) / (2 * xbs_all[i_a][in]);
  double tmax = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 - term2) / (2 * xbs_all[i_a][in]);

  double phi = 0.;
  fourvector p1b = xbp_all[i_a][in1].boost(xbp_all[i_a][in], xbs_all[i_a][in]);

  // invert second part
  fourvector k10 = xbp_all[i_a][out1].boost(xbp_all[i_a][in], xbs_all[i_a][in]);
  k10 = k10.rot_inverse(p1b);

  if (k10.x1()!=0 || k10.x2()!=0) {
    phi = atan2(k10.x2(),k10.x1());
    if (phi<0)
      phi+=2*M_PI;
  }
  else
    logger << LOG_DEBUG << "x1==x20 in inv_tchannel" << endl;

  r2=phi/(2*M_PI);

  if (!(r2>=0)) {
    logger << LOG_DEBUG << "r2<0 in inv_tchannel!" << endl;
    logger << LOG_DEBUG << "r2=" << r2 <<"; k10.x, k10.y=" << k10.x1() << ", " << k10.x2() << endl;
    logger << LOG_DEBUG << xbp_all[i_a][in] << endl;
    r2=0;
  }
  if (r2>1) {
    logger << LOG_DEBUG << "r2>1 in inv_tchannel!" << endl;
    logger << LOG_DEBUG << "r2=" << r2 <<"; k10.x, k10.y=" << k10.x1() << ", " << k10.x2() << endl;
    logger << LOG_DEBUG << xbp_all[i_a][in] << endl;
    r2=1;
  }

  // invert first part
  r1 = inv_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax);
  //  r1 = inv_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax, *this);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::c_tchannel_opt(int i_a, int m, int out, int in1, int in2, int out1, int out2, int i_r1, int i_r2){
  static Logger logger("phasespace_set::c_tchannel_opt");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double term1 = (xbs_all[i_a][out] + xbs_all[i_a][out1] - xbs_all[i_a][out2]) * (xbs_all[i_a][out] + xbs_all[i_a][in1] - xbs_all[i_a][in2]);
  double term2 = sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][out1], xbs_all[i_a][out2])) * sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2]));
  double term3 = (2 * sqrt(xbs_all[i_a][out])); // !!! why not xbsqrts_all[i_a][out] ???
  //  double term3 = (2 * xbsqrts_all[i_a][out]); // !!! why not xbsqrts_all[i_a][out] ??? - check if calculated ...

  double tmin = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 + term2) / (2 * xbs_all[i_a][out]);
  double tmax = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 - term2) / (2 * xbs_all[i_a][out]);
  double tmax_cut = 0.;
  logger << LOG_DEBUG_VERBOSE << "MC_optswitch = " << MC_optswitch << endl;
  logger << LOG_DEBUG_VERBOSE << "tmin = " << tmin << "   tmax = " << tmax << endl;
  if (MC_optswitch > 0 && xbs_all[i_a][out1] == 0.){
    vector<int> temp_out_particle = vectorint_from_binary(out1);
    if (temp_out_particle.size() == 1){
      vector<int> temp_smin_out = vectorbinary_from_binary(out - out1); // all outgoing ones, but out1
      double smin_out = 0.;
      for (int i = 0; i < temp_smin_out.size(); i++){smin_out += xbsqrts_all[i_a][temp_smin_out[i]];}
      smin_out = pow(smin_out, 2);
      tmax_cut = .5 *(xbs_all[i_a][0] - smin_out) * (sqrt(1. - (4. * xbs_all[i_a][0] * mapping_cut_pT[abs(csi->type_parton[0][temp_out_particle[0]])]) / pow((xbs_all[i_a][0] - smin_out), 2.)) - 1.);
    }
  }
  else if (MC_optswitch > 0 && xbs_all[i_a][out1] > 0.){
    vector<int> temp_out_particle = vectorint_from_binary(out1);
    if (temp_out_particle.size() == 1){
      vector<int> temp_smin_out = vectorbinary_from_binary(out - out1); // all outgoing ones, but out1
      double smin_out = 0.;
      for (int i = 0; i < temp_smin_out.size(); i++){smin_out += xbsqrts_all[i_a][temp_smin_out[i]];}
      smin_out = pow(smin_out, 2);

      tmax_cut = M2[abs(csi->type_parton[0][temp_out_particle[0]])] + .5 * (xbs_all[i_a][0] + M2[abs(csi->type_parton[0][temp_out_particle[0]])] - smin_out) * (sqrt(1. - (4. * xbs_all[i_a][0] * (M2[abs(csi->type_parton[0][temp_out_particle[0]])] + mapping_cut_pT[abs(csi->type_parton[0][temp_out_particle[0]])]) / pow((xbs_all[i_a][0] + M2[abs(csi->type_parton[0][temp_out_particle[0]])] - smin_out), 2.))) - 1.);
    }
  }

  if (MC_optswitch < 2){
    if ((MC_optswitch == 1) && (tmax_cut < tmax) && (tmin < tmax_cut)){tmax = tmax_cut;}
    xbs_all[i_a][in1 + out1] = -c_t_propagator(r[i_r1], m, -tmin, -tmax);
    //  xbs_all[i_a][in1 + out1] = -c_t_propagator(r[i_r1], m, -tmin, -tmax, *this);
  }
  else {
    if ((tmax_cut < tmax) && (tmin < tmax_cut)){
      logger << LOG_DEBUG << "(tmax_cut < tmax) && (tmin < tmax_cut):   " << "(" << tmax_cut << " < " << tmax << ") && (" << tmin << " < " << tmax_cut << ")" << endl;
      double r0 = 0.95;
      if (r[i_r1] < r0){
	xbs_all[i_a][in1 + out1] = -c_t_propagator(r[i_r1] / r0, m, -tmin, -tmax_cut);
	//	xbs_all[i_a][in1 + out1] = -c_t_propagator(r[i_r1] / r0, m, -tmin, -tmax_cut, *this);
      }
      else {

      }
    }
    else {
      xbs_all[i_a][in1 + out1] = -c_t_propagator(r[i_r1], m, -tmin, -tmax);
      //      xbs_all[i_a][in1 + out1] = -c_t_propagator(r[i_r1], m, -tmin, -tmax, *this);
    }
  }
  double phi = c_phi(i_r2);
  //  double phi = c_phi(r[i_r2]);
  double costheta = ((2 * xbs_all[i_a][out] * (xbs_all[i_a][in1 + out1] - xbs_all[i_a][out1] - xbs_all[i_a][in1])) + term1) / (term2);
  fourvector p10(((xbs_all[i_a][out] + xbs_all[i_a][in1] - xbs_all[i_a][in2]) / term3), 0, 0, (sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2])) / term3));
  fourvector p20(((xbs_all[i_a][out] - xbs_all[i_a][in1] + xbs_all[i_a][in2]) / term3), 0, 0, -(sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2])) / term3));
  fourvector k10(((xbs_all[i_a][out] + xbs_all[i_a][out1] - xbs_all[i_a][out2]) / term3), 0, 0, (sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][out1], xbs_all[i_a][out2])) / term3));
  k10 = k10.rot3(phi, costheta);
  fourvector k20 = p10 + p20 - k10;
  fourvector p1b = xbp_all[i_a][in1].boost(xbp_all[i_a][out], xbs_all[i_a][out]);
  xbp_all[i_a][out1] = (k10.rot(p1b)).boost(xbp_all[i_a][out].Pinv(), xbs_all[i_a][out]);
  xbp_all[i_a][out2] = (k20.rot(p1b)).boost(xbp_all[i_a][out].Pinv(), xbs_all[i_a][out]);
  xbp_all[i_a][in1 + out1] = xbp_all[i_a][in1] - xbp_all[i_a][out1];

  logger << LOG_DEBUG_VERBOSE << "xbp_all[i_a][" << in1 << "] = " << xbp_all[i_a][in1] << endl;
  logger << LOG_DEBUG_VERBOSE << "xbp_all[i_a][" << in2 << "] = " << xbp_all[i_a][in2] << endl;
  logger << LOG_DEBUG_VERBOSE << "xbp_all[i_a][" << out1 << "] = " << xbp_all[i_a][out1] << endl;
  logger << LOG_DEBUG_VERBOSE << "xbp_all[i_a][" << out2 << "] = " << xbp_all[i_a][out2] << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

double phasespace_set::g_tchannel_opt(int i_a, int m, int out, int in1, int in2, int out1, int out2){
  static Logger logger("phasespace_set::g_tchannel_opt");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  double gp = 0.;
  double term1 = (xbs_all[i_a][out] + xbs_all[i_a][out1] - xbs_all[i_a][out2]) * (xbs_all[i_a][out] + xbs_all[i_a][in1] - xbs_all[i_a][in2]);
  double term2 = sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][out1], xbs_all[i_a][out2])) * sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2]));
  double tmin = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 + term2) / (2 * xbs_all[i_a][out]);
  double tmax = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 - term2) / (2 * xbs_all[i_a][out]);
  if (tmin > 0.){logger << LOG_DEBUG << "tmin = " << tmin << endl;}
  if (tmax > 0.){logger << LOG_DEBUG << "tmax = " << tmax << endl; tmax = 0.;}
  double tmax_cut = 0.;

  if (MC_optswitch > 0 && xbs_all[i_a][out1] == 0.){
    vector<int> temp_smin_out = vectorbinary_from_binary(out - out1); // all outgoing ones, but out1
    double smin_out = 0.;
    for (int i = 0; i < temp_smin_out.size(); i++){smin_out += xbsqrts_all[i_a][temp_smin_out[i]];}
    smin_out = pow(smin_out, 2);
    vector<int> temp_out_particle = vectorint_from_binary(out1);
    if (temp_out_particle.size() == 1){
      tmax_cut = .5 *(xbs_all[i_a][0] - smin_out) * (sqrt(1. - (4. * xbs_all[i_a][0] * mapping_cut_pT[abs(csi->type_parton[0][temp_out_particle[0]])]) / pow((xbs_all[i_a][0] - smin_out), 2.)) - 1.);
    }
  }
  else if (MC_optswitch > 0 && xbs_all[i_a][out1] > 0.){
    vector<int> temp_smin_out = vectorbinary_from_binary(out - out1); // all outgoing ones, but out1
    double smin_out = 0.;
    for (int i = 0; i < temp_smin_out.size(); i++){smin_out += xbsqrts_all[i_a][temp_smin_out[i]];}
    smin_out = pow(smin_out, 2);
    vector<int> temp_out_particle = vectorint_from_binary(out1);
    if (temp_out_particle.size() == 1){
      tmax_cut = M2[abs(csi->type_parton[0][temp_out_particle[0]])] + .5 * (xbs_all[i_a][0] + M2[abs(csi->type_parton[0][temp_out_particle[0]])] - smin_out) * (sqrt(1. - (4. * xbs_all[i_a][0] * (M2[abs(csi->type_parton[0][temp_out_particle[0]])] + mapping_cut_pT[abs(csi->type_parton[0][temp_out_particle[0]])]) / pow((xbs_all[i_a][0] + M2[abs(csi->type_parton[0][temp_out_particle[0]])] - smin_out), 2.))) - 1.);
    }
  }

  if (MC_optswitch < 2){
    if ((MC_optswitch == 1) && (tmax_cut < tmax) && (tmin < tmax_cut)){tmax = tmax_cut;}
    gp = -g_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax) * (2 * sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2]))) / pi;
    //    gp = -g_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax, *this) * (2 * sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2]))) / pi;
  }
  else {
    if ((tmax_cut < tmax) && (tmin < tmax_cut)){

    }
    else {
      gp = -g_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax) * (2 * sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2]))) / pi;
      //      gp = -g_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax, *this) * (2 * sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2]))) / pi;
    }
  }

  if (gp != gp){
    logger << LOG_DEBUG << "error output: g_t != g_t" << endl;
    for (int ip = 1; ip < xbp_all[i_a].size(); ip*=2){
      logger << LOG_DEBUG << "xbp_all[i_a][" << setw(3) << ip << "] = " << xbp_all[i_a][ip] << endl;
    }
    logger << LOG_DEBUG << "m = " << m << endl;
    logger << LOG_DEBUG << "out = " << out << endl;
    logger << LOG_DEBUG << "in1 = " << in1 << endl;
    logger << LOG_DEBUG << "in2 = " << in2 << endl;
    logger << LOG_DEBUG << "out1 = " << out1 << endl;
    logger << LOG_DEBUG << "out2 = " << out2 << endl;
    logger << LOG_DEBUG << "term1 = " << term1 << endl;
    logger << LOG_DEBUG << "term2 = " << term2 << endl;
    logger << LOG_DEBUG << "tmin = " << tmin << endl;
    logger << LOG_DEBUG << "tmax = " << tmax << endl;
    logger << LOG_DEBUG << "tmax_cut = " << tmax_cut << endl;
    logger << LOG_DEBUG << "-xbs_all[i_a][in1 + out1] = " << -xbs_all[i_a][in1 + out1] << endl;
    logger << LOG_DEBUG << "g_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax) = " << g_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax) << endl;
    //    logger << LOG_DEBUG << "g_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax, *this) = " << g_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax, *this) << endl;
    logger << LOG_DEBUG << "lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2]) = " << lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2]) << endl;
    logger << LOG_DEBUG << "sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2])) = " << sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2])) << endl;
    logger << LOG_DEBUG << "lambda(xbs_all[i_a][out], xbs_all[i_a][out1], xbs_all[i_a][out2]) = " << lambda(xbs_all[i_a][out], xbs_all[i_a][out1], xbs_all[i_a][out2]) << endl;
    logger << LOG_DEBUG << "sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][out1], xbs_all[i_a][out2])) = " << sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][out1], xbs_all[i_a][out2])) << endl;
    logger << LOG_DEBUG << "xbs_all[i_a][out = " << out << "] = " << xbs_all[i_a][out] << endl;
    logger << LOG_DEBUG << "xbs_all[i_a][out1 = " << out1 << "] = " << xbs_all[i_a][out1] << endl;
    logger << LOG_DEBUG << "xbs_all[i_a][out2 = " << out2 << "] = " << xbs_all[i_a][out2] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  if (gp < 0.){logger << LOG_DEBUG << "gp < 0." << "   out = " << out << "   in1 = " << in1 << "   in2 = " << in2 << "   out1 = " << out1 << "   out2 = " << out2 << endl; gp = 0.;}
  return gp;
}

void phasespace_set::inv_tchannel_opt(int i_a, int m, int out, int in1, int in2, int out1, int out2, double & r1, double & r2){
  static Logger logger("phasespace_set::inv_tchannel_opt");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  double term1 = (xbs_all[i_a][out] + xbs_all[i_a][out1] - xbs_all[i_a][out2]) * (xbs_all[i_a][out] + xbs_all[i_a][in1] - xbs_all[i_a][in2]);
  double term2 = sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][out1], xbs_all[i_a][out2])) * sqrt(lambda(xbs_all[i_a][out], xbs_all[i_a][in1], xbs_all[i_a][in2]));
  //  double term3 = (2 * sqrt(xbs_all[i_a][out])); // !!! why not xbsqrts_all[i_a][out] ???
  //  double term3 = (2 * xbsqrts_all[i_a][out]); // !!! why not xbsqrts_all[i_a][out] ??? - check if calculated ...
  double tmin = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 + term2) / (2 * xbs_all[i_a][out]);
  double tmax = xbs_all[i_a][out1] + xbs_all[i_a][in1] - (term1 - term2) / (2 * xbs_all[i_a][out]);
  double tmax_cut = 0.;
  if (MC_optswitch > 0 && xbs_all[i_a][out1] == 0.){
    vector<int> temp_out_particle = vectorint_from_binary(out1);
    //    int temp_out_particle_size;
    //    int temp_out_particle_index = calc_strange_number(out1,temp_out_particle_size);
    if (temp_out_particle.size() == 1){
      //    if (temp_out_particle_size == 1) {
      vector<int> temp_smin_out = vectorbinary_from_binary(out - out1); // all outgoing ones, but out1
      double smin_out = 0.;
      //      smin_out=calc_strange_sum(out - out1);
      for (int i = 0; i < temp_smin_out.size(); i++){smin_out += xbsqrts_all[i_a][temp_smin_out[i]];}// cout << i << "   " << temp_smin_out[i] << "   " << xbsqrts_all[i_a][temp_smin_out[i]] << endl;}
      smin_out = pow(smin_out, 2);
      tmax_cut = .5 *(xbs_all[i_a][0] - smin_out) * (sqrt(1. - (4. * xbs_all[i_a][0] * mapping_cut_pT[abs(csi->type_parton[0][temp_out_particle[0]])]) / pow((xbs_all[i_a][0] - smin_out), 2.)) - 1.);
    }
  }
  else if (MC_optswitch > 0 && xbs_all[i_a][out1] > 0.){
    vector<int> temp_out_particle = vectorint_from_binary(out1);
    //    int temp_out_particle_size;
    //    int temp_out_particle_index=calc_strange_number(out1,temp_out_particle_size);
    if (temp_out_particle.size() == 1){
      //    if (temp_out_particle_size == 1) {
      vector<int> temp_smin_out = vectorbinary_from_binary(out - out1); // all outgoing ones, but out1
      double smin_out = 0.;
      for (int i = 0; i < temp_smin_out.size(); i++){smin_out += xbsqrts_all[i_a][temp_smin_out[i]];}
      smin_out = pow(smin_out, 2);

      tmax_cut = M2[abs(csi->type_parton[0][temp_out_particle[0]])] + .5 * (xbs_all[i_a][0] + M2[abs(csi->type_parton[0][temp_out_particle[0]])] - smin_out) * (sqrt(1. - (4. * xbs_all[i_a][0] * (M2[abs(csi->type_parton[0][temp_out_particle[0]])] + mapping_cut_pT[abs(csi->type_parton[0][temp_out_particle[0]])]) / pow((xbs_all[i_a][0] + M2[abs(csi->type_parton[0][temp_out_particle[0]])] - smin_out), 2.))) - 1.);
    }
  }

  double phi = 0.;
  fourvector p1b = xbp_all[i_a][in1].boost(xbp_all[i_a][out], xbs_all[i_a][out]);

  // invert second part
  fourvector k10 = xbp_all[i_a][out1].boost(xbp_all[i_a][out], xbs_all[i_a][out]);
  k10 = k10.rot_inverse(p1b);

  if (k10.x1()!=0 || k10.x2()!=0) {
    phi = atan2(k10.x2(),k10.x1());
    if (phi<0)
      phi+=2*M_PI;
  }
  else {
    logger << LOG_DEBUG << "x1=x2=0 in inv_tchannel_opt!" << endl;
    phi=0;
  }

  r2=phi/(2*M_PI);
  if (!(r2>=0)) {
    logger << LOG_DEBUG << "r2<0 in inv_tchannel_opt!" << endl;
    logger << LOG_DEBUG << "r2=" << r2 <<"; k10.x, k10.y=" << k10.x1() << ", " << k10.x2() << endl;
    r2=0;
  }
  if ((r2>1)) {
    logger << LOG_DEBUG << "r2>1 in inv_tchannel_opt!" << endl;
    logger << LOG_DEBUG << "r2=" << r2 <<"; k10.x, k10.y=" << k10.x1() << ", " << k10.x2() << endl;
    r2=1;
  }

  // invert first part
  if (MC_optswitch < 2){
    if ((MC_optswitch == 1) && (tmax_cut < tmax) && (tmin < tmax_cut)){tmax = tmax_cut;}
    r1 = inv_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax);
    //    r1 = inv_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax, *this);
  }
  else {
    logger << LOG_FATAL << "What the heck??" << endl;
    exit(0);
    if ((tmax_cut < tmax) && (tmin < tmax_cut)){
      logger << LOG_DEBUG << "(tmax_cut < tmax) && (tmin < tmax_cut): " << "(" << tmax_cut << " < " << tmax << ") && (" << tmin << " < " << tmax_cut << ")" << endl;
      double r0 = 0.95;
//      if (r1 < r0){
	r1 = r0 * c_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax_cut);
	//	r1 = r0*c_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax_cut, *this);
//      }
//      else {
//        // FIXME: ???
//      }
    }
    else {
      r1 = inv_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax);
      //      r1 = inv_t_propagator(-xbs_all[i_a][in1 + out1], m, -tmin, -tmax, *this);
    }
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

