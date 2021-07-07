#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
// **************************************************************************
// *                                                                        *
// *  various mappings for x in [0;1] -> h in [smin;smax]                   *
// *                                                                        *
// **************************************************************************
double h_pot_min(double r, double smin, double smax, double n){
  if (n == 0.){return smax * exp(r * log(smin / smax));}
  else{return (smax - smin) * pow(r, n) + smin;}
}
double g_pot_min(double s, double smin, double smax, double n){
  if (n == 0.){return 1. / (-log(smin / smax) * s);}
  else{return 1. / (pow(s - smin, (1. - 1. / n)) * pow(smax - smin, 1. / n) * n);}
}
double h_pot_max(double r, double smin, double smax, double n){
  if (n == 0.){return smax + smin - smax * exp(r * log(smin / smax));}
  //  else{return (smax - smin) * (1. - pow(r, n)) + smin;}
  else{return (smax - smin) * (1. - pow(1. - r, n)) + smin;}
}
double g_pot_max(double s, double smin, double smax, double n){
  if (n == 0.){return 1. / (log(smin / smax) * (s - (smin + smax)));}
  else{return 1. / (n * pow(smax - s, 1. - 1. / n) * pow(smax - smin, 1. / n));}
}



double h_propto_pot(double r, double smin, double n){
  if (n == 1.){return smin * exp(-log(smin) * r);}
  else{return pow(pow(smin, 1. - n) + (1. -  pow(smin, 1. - n)) * r, 1. / (1. - n));}
}
double g_propto_pot(double s, double smin, double n){
  if (n == 1.){return 1. / (-log(smin) * s);}
  else{return 1. / (((1. - pow(smin, 1. - n)) / (1. - n)) * pow(s, n));}
}



double h_propto_pot(double r, double smin, double smax, double nu, double cut_technical){
  static Logger logger("h_propto_pot");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double s;
  if (nu != 1.){s = pow(r * pow(smax + cut_technical, 1. - nu) + (1. - r) * pow(smin + cut_technical, 1. - nu), (1. / (1. - nu))) - cut_technical;}
  else {s = exp(r * log(smax + cut_technical) + (1. - r) * log(smin + cut_technical)) - cut_technical;}
  return s;
}

double g_propto_pot(double s, double smin, double smax, double nu, double cut_technical){
  static Logger logger("g_propto_pot");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double gs;
  if (nu != 1.){gs = (1. - nu) / ((pow(smax + cut_technical, 1. - nu) - pow(smin + cut_technical, 1. - nu)) * pow(s + cut_technical, nu));}
  else {gs = 1. / ((log(smax + cut_technical) - log(smin + cut_technical)) * (s + cut_technical));}
  return gs;
}

double inv_propto_pot(double s, double smin, double smax, double nu, double cut_technical){
  static Logger logger("inv_propto_pot");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // temporary !!!
  int switch_console_output_phasespace_issue = 0;

  double r;
  if (nu != 1.) {
    double sminm0 = pow(smin + cut_technical, 1. - nu);
    r = (pow(s + cut_technical, 1. - nu) - sminm0) / (pow(smax + cut_technical, 1. - nu) - sminm0);
  } 
  else {r = log((s + cut_technical) / (smin + cut_technical)) / log((smax + cut_technical) / (smin + cut_technical));}

  if (munich_isnan(r)) {
    if (switch_console_output_phasespace_issue){
      cout << "r is NaN in inv_propto_pot!" << endl;
      cout << "r = " << r <<"; smin, s, smax=" << smin << ", " << s << ", " << smax << endl;
    }
    //    r = 0.5;
    if (s <= smin){r = cut_technical;}
    else if (s >= smax){r = 1. - cut_technical;}
  }

  else if (r < 0.){
    if (switch_console_output_phasespace_issue){
      cout << "r < 0 in inv_propto_pot!" << endl;
      cout << "r = " << r <<"; smin, s, smax=" << smin << ", " << s << ", " << smax << endl;
    }
    r = cut_technical;
  }
  else if (r > 1.){
    if (switch_console_output_phasespace_issue){
      cout << "r > 1 in inv_propto_pot!" << endl;
      cout << "r = " << r <<"; smin, s, smax=" << smin << ", " << s << ", " << smax << endl;
    }
    r = 1. - cut_technical;
  }

  return r;
}

// check difference to h_propto_pot_mod !!! (in particular in nu == 1 case) !!!

double h_propto_pot_mod(double r, double smin, double smax, double nu, double cut_technical){
  static Logger logger("h_propto_pot_mod");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double s;
  if (nu != 1.){s = pow(r * pow(smax + cut_technical, 1. - nu) + (1. - r) * pow(smin + cut_technical, 1. - nu), (1. / (1. - nu))) - cut_technical;}
  else {s = smin + (smax - smin) * (exp(r * log(1. + cut_technical) + (1. - r) * log(cut_technical)) - cut_technical);}
  return s;
}

double g_propto_pot_mod(double s, double smin, double smax, double nu, double cut_technical){
  static Logger logger("g_propto_pot_mod");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double gs;
  if (nu != 1.){gs = (1. - nu) / ((pow(smax + cut_technical, 1. - nu) - pow(smin + cut_technical, 1. - nu)) * pow(s + cut_technical, nu));}
  else{gs = 1. / ((log(1. + cut_technical) - log(cut_technical)) * ((s - smin) + cut_technical * (smax - smin)));}
  return gs;
}


// probably unused at the moment !!!
double h_propto_pot_symm(double r, double smin, double smax, double nu, double cut_technical){
  static Logger logger("h_propto_pot_symm");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // temporary !!!
  int switch_console_output_phasespace_issue = 0;

  double s;
  double smin_smax_2 = .5 * (smax - smin);
  if (nu != 1.){
    if (r <= 0.5){
      s = smin + smin_smax_2 * pow(2 * r * pow(1. + cut_technical, 1. - nu) + (1. - 2 * r) * pow(cut_technical, 1. - nu), (1. / (1. - nu))) - cut_technical;
    }
    else {
      s = smax - smin_smax_2 * pow(2 * (1. - r) * pow(1. + cut_technical, 1. - nu) + (2 * r - 1.) * pow(cut_technical, 1. - nu), (1. / (1. - nu))) - cut_technical;
    }
  }
  else{
    if (r <= 0.5){
      s = smin + smin_smax_2 * (exp(2 * r * log(1. + cut_technical) + (1. - 2 * r) * log(cut_technical)) - cut_technical);
    }
    else {
      s = smax - smin_smax_2 * (exp(2 * (1. - r) * log(1. + cut_technical) + (2 * r - 1.) * log(cut_technical)) - cut_technical);
    }
    if (switch_console_output_phasespace_issue){
      if (s < 0. || s > 1.){cout << "s = " << setw(25) << setprecision(16) << s << "    sminsmax_2 = " << .5 * (smax + smin) << endl;}
    }
  }
  return s;
}
double g_propto_pot_symm(double s, double smin, double smax, double nu, double cut_technical){
  static Logger logger("g_propto_pot_symm");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double gs;
  double sminsmax_2 = .5 * (smax + smin);
  double smin_smax_2 = .5 * (smax - smin);
  if (nu != 1.){
    if (s <= sminsmax_2){
      gs = (1. - nu) / (2 * pow(smin_smax_2, 1. - nu) * (pow(1. + cut_technical, 1. - nu) - pow(cut_technical, 1. - nu)) * pow(s - smin + cut_technical * smin_smax_2, nu));
    }
    else {
      gs = (1. - nu) / (2 * pow(smin_smax_2, 1. - nu) * (pow(1. + cut_technical, 1. - nu) - pow(cut_technical, 1. - nu)) * pow(smax - s + cut_technical * smin_smax_2, nu));
    }
  }
  else {
    if (s <= sminsmax_2){
      gs = 1. / (2 * (log(1. + cut_technical) - log(cut_technical)) * ((s - smin) + cut_technical * smin_smax_2));
    }
    else {
      gs = 1. / (2 * (log(1. + cut_technical) - log(cut_technical)) * ((smax - s) + cut_technical * smin_smax_2));
    }
  }
  return gs;
}



// **************************************************************************
// *                                                                        *
// *  mapping for angles                                                    *
// *                                                                        *
// **************************************************************************
/*
double c_phi(double r){
  return 2 * pi * r;
}

double inv_phi(double phi){
  return phi / (2 * pi);
}

double c_costheta(double r){
  return 2 * r - 1.;
}

double inv_costheta(double costheta){
  return (costheta + 1) / 2;
}
*/


//#include "old.phasespace.generator.cpp"
