#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

//  important issues:
// - check which version in H2 is correct (potentially doubled term in VT2 contribution, which is already included via H1full -> relevant if LQ =/= 0) !!!
// - check if 'scale-dependent' parts are really scale dependent !!!
// - check missing LQ contributions in gg channel !!!

// {{{ observable_set::calculate_A_F()
void observable_set::calculate_A_F(){
  static Logger logger("observable_set::calculate_A_F");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (initial_gg){A_F = 2 * (QT_H1_delta - pi2_6 * C_A);}
  else if (initial_qqx){A_F = 2 * (QT_H1_delta - pi2_6 * C_F);}
  else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ observable_set::calculate_B2_A_F()
void observable_set::calculate_B2_A_F(){
  static Logger logger("observable_set::calculate_B2_A_F");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (initial_gg){B2_A_F = B2g + .5 * beta0 * A_F;}
  else if (initial_qqx){B2_A_F = B2q + .5 * beta0 * A_F;}
  else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ observable_set::calculate_Gamma2()
void observable_set::calculate_Gamma2(){
  static Logger logger("observable_set::calculate_Gamma2");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  if (QT_finalstate_massive_coloured){

    double gammacusp2_v = 0.;

    if (initial_pdf_diag){
      //    if (initial_pdf_gg){ // check which one is needed here !!!
      //    if ((name_process[0] == 'g' && name_process[1] == 'g') || (name_process[0] != 'g' && name_process[1] != 'g')){
      logger << LOG_FATAL << "VT or CT2 term for gg and qqx channel is not validated (and not completely implemented) yet !!!" << endl;

      Gamma2 = 0.;

      double q2 = (p_parton[0][3] + p_parton[0][4]).m2();
      double v34 = sqrt(1. - pow(msi.M2_t / (p_parton[0][3] * p_parton[0][4]), 2));
      double log_v34 = log((1. + v34) / (1. - v34));
      for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
	if (QT_correlationoperator[i_c].type_combination == 1){
	  //	  Gamma2 += .5 * QT_ME2_cf[i_c];
	}
	else if (QT_correlationoperator[i_c].type_combination == 2){
	  Gamma2 += .5 * 2 * gammacusp2_v * QT_ME2_cf[i_c];
	}
	else if (QT_correlationoperator[i_c].type_combination == 3){
	  Gamma2 += .5 * log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * msi.M2_t)) * QT_ME2_cf[i_c] * gammacusp2;
	}
 	else if (QT_correlationoperator[i_c].type_combination == 0){
	  Gamma2 += .5 * (-4 * gamma2Q); // multiplied by 1 (born / born)
	}
	// three-correlator don't contribute for ttbar
      }

      // Ft1 is also needed elsewhere and should thus better be calculated only once !!!
      double Ft1 = 0.;
      //      double v34 = sqrt(1. - pow(2 * msi.M2_t / (2 * (p_parton[0][3] * p_parton[0][4])), 2));
      double frac_v34 = (1. + v34) / (1. - v34);
      //      double log_v34 = log(frac_v34);
      double y34 = p_parton[0][3].rapidity() - p_parton[0][4].rapidity();
      double dilogpT2m2 = gsl_sf_dilog(-p_parton[0][3].pT2() / msi.M2_t);
      double logmT2m2 = log(p_parton[0][3].ET2() / msi.M2_t);
      double L34 = log_v34 * logmT2m2
	- 2 * gsl_sf_dilog(2 * v34 / (1. + v34))
	- .25 * pow(log_v34, 2)
	+ 2 * (gsl_sf_dilog(1. - exp(y34) / sqrt(frac_v34)) + gsl_sf_dilog(1. - exp(-y34) / sqrt(frac_v34)) + .5 * pow(y34, 2));
      for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
	if (QT_correlationoperator[i_c].type_combination == 1){
	  Ft1 += (logmT2m2 + dilogpT2m2) * QT_ME2_cf[i_c];
	}
	else if (QT_correlationoperator[i_c].type_combination == 2){
	  Ft1 += (2 * dilogpT2m2 + (1. / v34) * L34) * QT_ME2_cf[i_c];
	}
      }
      //      QT_H1_delta += .5 * Ft1;
      //  H1 -> H1 - It1  (It1 = 2 Re < M0 | It1 | M0 >) 
    
      // [Gamma_t^1, F_t^1] -> only in gg channel !!!    

      Gamma2 += .5 * commutator_Gammat1_Ft1;

      Gamma2 += .5 * pi * beta0 * Ft1;


      exit(1);
    }
    else if (initial_pdf_qqx){
      logger << LOG_FATAL << "VT or CT2 term for qqx channel is not implemented yet !!!" << endl;
      exit(1);
    }
    else if (initial_pdf_gq){
      logger << LOG_FATAL << "Gamma²_t does not contribute in the gq channel !!!" << endl;
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ observable_set::calculate_commutator_Gamma1_Ft1()
void observable_set::calculate_commutator_Gamma1_Ft1(){
  static Logger logger("observable_set::calculate_commutator_Gamma1_Ft1");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  commutator_Gammat1_Ft1 = 0.;

  double q2 = (p_parton[0][3] + p_parton[0][4]).m2();
  double v34 = sqrt(1. - pow(2 * msi.M2_t / (2 * (p_parton[0][3] * p_parton[0][4])), 2));
  double frac_v34 = (1. + v34) / (1. - v34);
  double log_v34 = log(frac_v34);
  double y34 = p_parton[0][3].rapidity() - p_parton[0][4].rapidity();
  double dilogpT2m2 = gsl_sf_dilog(-p_parton[0][3].pT2() / msi.M2_t);
  double logmT2m2 = log(p_parton[0][3].ET2() / msi.M2_t);
  double L34 = log_v34 * logmT2m2
    - 2 * gsl_sf_dilog(2 * v34 / (1. + v34))
    - .25 * pow(log_v34, 2)
    + 2 * (gsl_sf_dilog(1. - exp(y34) / sqrt(frac_v34)) + gsl_sf_dilog(1. - exp(-y34) / sqrt(frac_v34)) + .5 * pow(y34, 2));

  for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
    if (QT_correlationoperator[i_c].type_combination == 999){// some new type, probably in a new correlation-operator class !!!
      commutator_Gammat1_Ft1 += (-.25 * log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * msi.M2_t)))
	* (2 * dilogpT2m2 + (1. / v34) * L34) 
	* QT_ME2_cf[i_c]; //  * QT_ME2_cf[i_c] must be the 4-colour-correlator commutator !!!
    }
  }
  // T_
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ observable_set::calculate_Gamma1_Ft1()
void observable_set::calculate_Gamma1_Ft1(){
  static Logger logger("observable_set::calculate_Gamma1_Ft1");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  Gammat1_Ft1 = 0.;
  /*
  double q2 = (p_parton[0][3] + p_parton[0][4]).m2();
  double dilogpT2m2 = gsl_sf_dilog(-p_parton[0][3].pT2() / msi.M2_t);
  double v34 = sqrt(1. - pow(2 * msi.M2_t / (2 * (p_parton[0][3] * p_parton[0][4])), 2));
  Gammat1_Ft1 += (-.25 * log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * msi.M2_t)))
    * (2 * dilogpT2m2 + (1. / v34) * L34) 
    * QT_ME2_cf[i_c]; //  * QT_ME2_cf[i_c] must be the 4-colour-correlator commutator !!!
  // T_
  */
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ observable_set::calculate_Gamma1_squared()
void observable_set::calculate_Gamma1_squared(){
  static Logger logger("observable_set::calculate_Gamma1_squared");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  if (QT_finalstate_massive_coloured){
    // 4 * Re (-1/4) (-1/4) * {...}

    if (initial_pdf_diag){
      //    if (initial_pdf_gg || initial_pdf_qqx){
      Gamma1squared = 0.;
      double q2 = (p_parton[0][3] + p_parton[0][4]).m2();
      double v34 = sqrt(1. - pow(msi.M2_t / (p_parton[0][3] * p_parton[0][4]), 2));
      double log_v34 = log((1. + v34) / (1. - v34));
      for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
	if (QT_correlationoperator[i_c].type_combination == 1){
	  Gamma1 += .5 * QT_ME2_cf[i_c];
	}
	else if (QT_correlationoperator[i_c].type_combination == 2){
	  Gamma1 += .5 * (1. / v34) * log_v34 * QT_ME2_cf[i_c];
	}
	else if (QT_correlationoperator[i_c].type_combination == 3){
	  Gamma1 += .5 * log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * msi.M2_t)) * QT_ME2_cf[i_c];
	}
      }
      //    }
      // Gammat = -2 * Gammat1(Eg.33, arXiv:1408.4564) !!! -> 2 * Re(Gammat1)
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ observable_set::calculate_Gamma1()
void observable_set::calculate_Gamma1(){
  static Logger logger("observable_set::calculate_Gamma1");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  if (QT_finalstate_massive_coloured){
    Gamma1 = 0.;
    /*
    if (initial_pdf_gg){
      logger << LOG_FATAL << "VT or CT2 term for gg channel is not implemented yet !!!" << endl;
      exit(1);
    }
    if (initial_pdf_qqx){
      logger << LOG_FATAL << "VT or CT2 term for qqx channel is not implemented yet !!!" << endl;
      exit(1);
    }
    // only for gq channel !!!
    if (initial_pdf_gq){
    */
      double q2 = (p_parton[0][3] + p_parton[0][4]).m2();
      double v34 = sqrt(1. - pow(msi.M2_t / (p_parton[0][3] * p_parton[0][4]), 2));
      double log_v34 = log((1. + v34) / (1. - v34));
      for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
	if (QT_correlationoperator[i_c].type_combination == 1){
	  Gamma1 += .5 * QT_ME2_cf[i_c];
	}
	else if (QT_correlationoperator[i_c].type_combination == 2){
	  Gamma1 += .5 * (1. / v34) * log_v34 * QT_ME2_cf[i_c];
	}
	else if (QT_correlationoperator[i_c].type_combination == 3){
	  Gamma1 += .5 * log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * msi.M2_t)) * QT_ME2_cf[i_c];
	}
      }
      /*
    }
      */
    // Gammat = -2 * Gammat1(Eg.33, arXiv:1408.4564) !!! -> 2 * Re(Gammat1)
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ observable_set::calculate_sigma12(double & pdf_factor_x1x2)
double observable_set::calculate_sigma12(double & pdf_factor_x1x2){
  static Logger logger("observable_set::calculate_sigma12");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double sig12 = 0.;
  if (initial_gg){sig12 = -0.5 * A1g * pdf_factor_x1x2;}
  else if (initial_qqx){sig12 = -0.5 * A1q * pdf_factor_x1x2;}
  else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig12;
}
// }}}
// {{{ observable_set::calculate_sigma11(double & pdf_factor_x1x2, double & tH1F, double & LQ)
double observable_set::calculate_sigma11(double & pdf_factor_x1x2, double & tH1F, double & LQ){
  static Logger logger("observable_set::calculate_sigma11");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double sig11 = 0.;
  if (initial_gg){
    sig11 = 
      - pdf_factor_x1x2 * B1g 
      - tH1F;
    // LQ contribution not implemented yet !!!
    // - pdf_factor_x1x2 * A1g * LQ; // ???
  }
  else if (initial_qqx){
    sig11 = 
      - B1q * pdf_factor_x1x2 
      - tH1F 
      - pdf_factor_x1x2 * A1q * LQ;
  }

  if (QT_finalstate_massive_coloured){
    sig11 += pdf_factor_x1x2 * Gamma1;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig11;
}
// }}}
// {{{ observable_set::calculate_sigma24(double & pdf_factor_x1x2)
double observable_set::calculate_sigma24(double & pdf_factor_x1x2){
  static Logger logger("observable_set::calculate_sigma24");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double sig24 = 0.;
  if (initial_gg){sig24 = 1. / 8 * A1g * A1g * pdf_factor_x1x2;}
  else if (initial_qqx){sig24 = 1. / 8 * A1q * A1q * pdf_factor_x1x2;}
  else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig24;
}
// }}}
// {{{ observable_set::calculate_sigma23(double & pdf_factor_x1x2, double & sig11)
double observable_set::calculate_sigma23(double & pdf_factor_x1x2, double & sig11){
  static Logger logger("observable_set::calculate_sigma23");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double sig23 = 0.;
  if (initial_gg){
    sig23 = 
      - A1g * beta0 / 3. * pdf_factor_x1x2 
      - A1g * 0.5 * sig11;
  }
  else if (initial_qqx){
    sig23 = 
      - beta0 * A1q / 3 * pdf_factor_x1x2 
      - 0.5 * A1q * sig11;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig23;
}
// }}}
// {{{ observable_set::calculate_sigma22(double & pdf_factor_x1x2, double & sig11, double & tH1, double & tH1F, double & H1full, double & tgaga, double & LR, double & LF, double & LQ)
double observable_set::calculate_sigma22(double & pdf_factor_x1x2, double & sig11, double & tH1, double & tH1F, double & H1full, double & tgaga, double & LR, double & LF, double & LQ){
  static Logger logger("observable_set::calculate_sigma22");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double sig22 = 0.;
  //  if (name_process[0] == 'g'){
  if (initial_gg){
    sig22 = 
      - 0.5 * A1g * H1full
      + 0.5 * (beta0 * A1g * LR - A2g) * pdf_factor_x1x2
      - 0.5 * (B1g - beta0) * sig11
      + 0.5 * B1g * tH1F
      + 0.5 * tgaga;
    // LQ contribution missing !!!
  }
  //  else {
  else if (initial_qqx){
   sig22 = 
      - 0.5 * A1q * H1full 
      + 0.5 * (beta0 * A1q * LR - A2q) * pdf_factor_x1x2 
      - 0.5 * (B1q - beta0) * sig11 
      + 0.5 * B1q * tH1F 
      + 0.5 * tgaga;
    
    // Q dependence
    sig22 -= 0.5 * A1q * LQ * beta0 * pdf_factor_x1x2;
    sig22 += 0.5 * A1q * LQ * tH1F;
    sig22 -= 0.5 * A1q * LQ * sig11;
  }

  if (QT_finalstate_massive_coloured){
    /* check selection !!!
    if (initial_gg && initial_pdf_gg){
      sig22 -= 0.5 * Gamma1 * B1g * pdf_factor_x1x2;
      // Gamma1squared ???
    }
    if (initial_qqx && initial_pdf_qqx){
      sig22 -= 0.5 * Gamma1 * B1q * pdf_factor_x1x2;
      sig22 += 0.5 * Gamma1squared * pdf_factor_x1x2;
    }
    */
    if (initial_pdf_gq){
      // only contribution needed for gq channel (more colour correlators in diagonal channels) !!!
      sig22 -= Gamma1 * tH1F;  
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig22;
}
// }}}
// {{{ observable_set::calculate_sigma21(double & pdf_factor_x1x2, double & sig11, double & tH1, double & tH1F, double & H1full, double & tgaga, double & tcga, double & tgamma2, double & LR, double & LF, double & LQ, double & A_F)
double observable_set::calculate_sigma21(double & pdf_factor_x1x2, double & sig11, double & tH1, double & tH1F, double & H1full, double & tgaga, double & tcga, double & tgamma2, double & LR, double & LF, double & LQ, double & A_F){
  static Logger logger("observable_set::calculate_sigma21");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double B2_A_F = 0.;
  static double sig21 = 0.;
  //  if (name_process[0] == 'g'){
  if (initial_gg){
    // B2_A_F recalculated - function 'calculate_B2_A_F()' could be used !!!
    B2_A_F = B2g + .5 * beta0 * A_F;
    // B2g must be B2(g)_A_F in the following formula (misleading name in counterterm_gg)
    sig21 = 
      - beta0 * LR * sig11 
      - B1g * H1full 
      - LF * tgaga 
      - B2_A_F * pdf_factor_x1x2 
      + beta0 * tH1
      - tcga 
      - tgamma2 
      - C1ggdelta(A_F) * tH1F 
      - 2 * Delta2gg * pdf_factor_x1x2
      + 2 * beta0 * LR * tH1F;
    //20170522:
    // + 2 * beta0 * LR * (B1g * pdf_factor_x1x2  seems doubled
    // (contained in  -B1g * H1full)
    // org:      + 2 * beta0 * LR * (B1g * pdf_factor_x1x2 + tH1F);
    // new:      + 2 * beta0 * LR * tH1F;

    logger << LOG_DEBUG_VERBOSE << "- beta0 * LR * sig11 = " << - beta0 * LR * sig11 << endl; //ok
    logger << LOG_DEBUG_VERBOSE << "- B1g * H1full = " << - B1g * H1full << endl;//different !!!
    logger << LOG_DEBUG_VERBOSE << "- LF * tgaga = " << - LF * tgaga << endl; //ok
    logger << LOG_DEBUG_VERBOSE << "- B2_A_F * pdf_factor_x1x2 = " << - B2_A_F * pdf_factor_x1x2 << endl; //ok
    logger << LOG_DEBUG_VERBOSE << "+ beta0 * tH1 = " << + beta0 * tH1 << endl; //ok
    logger << LOG_DEBUG_VERBOSE << "- tcga = " << - tcga << endl; //ok
    logger << LOG_DEBUG_VERBOSE << "- tgamma2 = " << - tgamma2 << endl; //ok
    logger << LOG_DEBUG_VERBOSE << "- C1ggdelta(A_F) * tH1F = " << - C1ggdelta(A_F) * tH1F << endl; //ok
    logger << LOG_DEBUG_VERBOSE << "- 2 * Delta2gg * pdf_factor_x1x2 = " << - 2 * Delta2gg * pdf_factor_x1x2 << endl; //ok
    logger << LOG_DEBUG_VERBOSE << "+ 2 * beta0 * LR * (B1g * pdf_factor_x1x2 + tH1F) = " << + 2 * beta0 * LR * (B1g * pdf_factor_x1x2 + tH1F) << endl; //ok

  }
  //  else {
  else if (initial_qqx){
    // B2_A_F recalculated - function 'calculate_B2_A_F()' could be used !!!
    B2_A_F = B2q + .5 * beta0 * A_F;
    sig21 = 
      - beta0 * LR * sig11 
      - B1q * H1full 
      - LF * tgaga 
      - B2_A_F * pdf_factor_x1x2 
      + beta0 * tH1 
      - tcga 
      - tgamma2 
      - C1qqdelta(A_F) * tH1F 
      - 2 * Delta2qq * pdf_factor_x1x2;
    
    // Q dependence // check if identical !!!
    sig21 += LQ * sig11 * beta0;
    sig21 -= LQ * H1full * A1q;
    sig21 += LQ * (tgaga + (B1q + 0.5 * A1q * LQ) * tH1F);
    sig21 -= LQ * A2q * pdf_factor_x1x2;
    //  sig21 += LQ*(sig11*(beta0-B1q-0.5*A1q*LQ) - A2q*pdf_factor_x1x2 + A1q*(2*beta0*LR*pdf_factor_x1x2-tH1+(LQ-LF)*tH1F) + tgaga);
  }

  if (QT_finalstate_massive_coloured){
    //if (name_process[0] == 'g'){
    if (initial_gg){
      sig21 += Gamma1 * H1full ;
    }
    //    else {
    else if (initial_qqx){
      sig21 += Gamma1 * H1full;
      // 20170522: check if this is correct !!!
      // 2 * beta0 * LR * B1q * pdf_factor_x1x2  term could be doubled !!!
      // (contained already in H1full ???)
      sig21 += 2 * beta0 * LR * (B1q * pdf_factor_x1x2 + tH1F);  // alpha_S² at LO, not particularly because of mass - but it appears only in ttx production so far !!!
    }
    else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}
  }
  // needs to be refined for qqx and gg channel !!!
  /*
  if (QT_finalstate_massive_coloured){
    if ((name_process[0] == 'g' && name_process[1] != 'g') || (name_process[0] != 'g' && name_process[1] == 'g')){
      sig21 += Gamma1 * H1full;
      //    }
    else if (name_process[0] == 'g' && name_process[1] == 'g'){
      sig21 += Gamma2 * pdf_factor_x1x2;
      sig21 -= 2 * beta0 * LR * Gamma1 * pdf_factor_x1x2;
      //      sig21 += Gamma1_H1full; // colour-correlated version !!!
    }
    else if (name_process[0] != 'g' && name_process[1] != 'g'){
      sig21 += Gamma2 * pdf_factor_x1x2;
      sig21 -= 2 * beta0 * LR * Gamma1 * pdf_factor_x1x2;
      //      sig21 += Gamma1_H1full; // colour-correlated version !!!
      sig21 += 2 * beta0 * LR * (B1q * pdf_factor_x1x2 + tH1F);  // alpha_S² at LO, not particularly because of mass - but it appears only in ttx production so far !!!
    }
    else {
      // not defined !!!
      //      exit(1);
    }
    }
*/

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig21;
}
// }}}
// no LR dependence !!!
// {{{ observable_set::calculate_H1(double pdf_factor_x1x2, double tH1, double tH1F, double LR, double LF, double LQ)
double observable_set::calculate_H1(double pdf_factor_x1x2, double tH1, double tH1F, double LR, double LF, double LQ){
  static Logger logger("observable_set::calculate_H1");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double result_H1 = 0.;
  if (name_process[0] == 'g'){
    result_H1 = 
      + tH1 
      + LF * tH1F 
      - 2 * beta0 * LR * pdf_factor_x1x2;
    // LQ-dependent term missing !!!
    logger << LOG_DEBUG_VERBOSE << setw(40) << "result_H1" << " = " << result_H1 << endl;
    logger << LOG_DEBUG_VERBOSE << setw(40) << "tH1" << " = " << tH1 << endl;
    logger << LOG_DEBUG_VERBOSE << setw(40) << "LF * tH1F" << " = " << LF * tH1F << endl;
    logger << LOG_DEBUG_VERBOSE << setw(40) << "- 2 * beta0 * LR * pdf_factor_x1x2" << " = " << - 2 * beta0 * LR * pdf_factor_x1x2 << endl;

  }
  else {
    result_H1 = 
      + tH1 
      + (LF - LQ) * tH1F 
      - LQ * pdf_factor_x1x2 * (B1q + 0.5 * A1q * LQ);
  }

  if (QT_finalstate_massive_coloured){
    if (name_process[0] == 'g'){
    }
    else {
      result_H1 += - 2 * beta0 * LR * pdf_factor_x1x2;  // alpha_S² at LO, not particularly because of mass - but it appears only in ttx production so far !!!
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return result_H1;
}
// }}}
// {{{ observable_set::calculate_H2(double pdf_factor_x1x2, double tH1, double tH1F, double H1full, double sig11, double tgaga, double tcga, double tgamma2, double tH2, double LR, double LF, double LQ)
double observable_set::calculate_H2(double pdf_factor_x1x2, double tH1, double tH1F, double H1full, double sig11, double tgaga, double tcga, double tgamma2, double tH2, double LR, double LF, double LQ){
  static Logger logger("observable_set::calculate_H2");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG_VERBOSE << "name_process[0] = " << name_process[0] << endl;

  static double result_H2 = 0.;
  if (name_process[0] == 'g'){
    result_H2 = 1.  //global_factor 
      * (tH2 // equivalent in qqx
	 + 0.5 * beta0 * LF * LF * tH1F // equivalent in qqx
	 + tgamma2 * LF // equivalent in qqx
	 // 20170523: replace the following term by H1full, analogously to qqbar case:
	 //	 - 2 * beta0 * LR * (tH1  
	 //			     + LF * tH1F 
	 //			     - 2 * beta0 * LR * pdf_factor_x1x2) // to be analyzed...
	 - 2 * beta0 * LR * H1full // equivalent in qqx
	 // until here
	 - beta0 * LF * LR * tH1F // LF * LR term in qqx hidden somewhere... 3 times in total (including the part in H1full ???)
	 + LF * tcga // equivalent in qqx
	 - beta0 * LR * tH1 // to be analyzed... 3 times in total (including the part in H1full ???)
	 + 0.5 * LF * LF * tgaga // equivalent in qqx
	 - 2 * (0.5 * (beta0 * LR) * (beta0 * LR) + beta1 * LR) * pdf_factor_x1x2 // to be analyzed... beta0²*LR² term -3 times in total (including the part in H1full ???)
	 // Include missing delta term from C*gamma (no factor 2 here !)
	 // xmsq = xmsq + asopi**2*(LF*C1ggdelta*tH1stF)
	 + LF * C1ggdelta(A_F) * tH1F // equivalent in qqx
	 // Include missing term from contact term in 2 loop AP
	 // xmsq = xmsq + asopi**2*(2*Delta2gg*tdelta)*LF
	 + 2 * Delta2gg * pdf_factor_x1x2 * LF); // equivalent in qqx

    logger << LOG_DEBUG_VERBOSE << setw(50) << "LF" << " = " << LF << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "LR" << " = " << LR << endl;
    logger << LOG_DEBUG_VERBOSE << setw(50) << "tH2" << " = " << tH2 << endl; // different from H2_M_gg !!!
    logger << LOG_DEBUG_VERBOSE << setw(50) << "+ 0.5 * beta0 * LF * LF * tH1F" << " = " << + 0.5 * beta0 * LF * LF * tH1F << endl; //ok
    logger << LOG_DEBUG_VERBOSE << setw(50) << "+ tgamma2 * LF" << " = " << + tgamma2 * LF << endl; //ok
    logger << LOG_DEBUG_VERBOSE << setw(50) << "- 2 * beta0 * LR * (tH1 + LF * tH1F - 2 * beta0 * LR * pdf_factor_x1x2)" << " = " << - 2 * beta0 * LR * (tH1 + LF * tH1F - 2 * beta0 * LR * pdf_factor_x1x2) << endl; //ok
    logger << LOG_DEBUG_VERBOSE << setw(50) << "- 2 * beta0 * LR * H1full" << " = " << - 2 * beta0 * LR * H1full << endl; // to replace the previous term
    logger << LOG_DEBUG_VERBOSE << setw(50) << "- beta0 * LF * LR * tH1F" << " = " << - beta0 * LF * LR * tH1F << endl; //ok
    logger << LOG_DEBUG_VERBOSE << setw(50) << "+ LF * tcga" << " = " << + LF * tcga << endl; //ok
    logger << LOG_DEBUG_VERBOSE << setw(50) << "- beta0 * LR * tH1" << " = " << - beta0 * LR * tH1 << endl; //ok

    logger << LOG_DEBUG_VERBOSE << setw(50) << "+ 0.5 * LF * LF * tgaga" << " = " << + 0.5 * LF * LF * tgaga << endl; //ok
    logger << LOG_DEBUG_VERBOSE << setw(50) << "- 2 * (0.5 * (beta0 * LR) * (beta0 * LR) + beta1 * LR) * pdf_factor_x1x2" << " = " << - 2 * (0.5 * (beta0 * LR) * (beta0 * LR) + beta1 * LR) * pdf_factor_x1x2 << endl; //ok
    logger << LOG_DEBUG_VERBOSE << setw(50) << "+ LF * C1ggdelta(A_F) * tH1F" << " = " << + LF * C1ggdelta(A_F) * tH1F << endl; //ok
    logger << LOG_DEBUG_VERBOSE << setw(50) << "+ 2 * Delta2gg * pdf_factor_x1x2 * LF" << " = " << + 2 * Delta2gg * pdf_factor_x1x2 * LF << endl; //ok
    //    logger << LOG_DEBUG_VERBOSE << setw(50) << "" << " = " <<  << endl;
    //    logger << LOG_DEBUG_VERBOSE << setw(50) << "" << " = " <<  << endl;
  }
  else {
    result_H2 = 1.  //global_factor 
      * (tH2
	 + 0.5 * beta0 * pow(LF, 2) * tH1F 
	 + tgamma2 * LF 
	 - beta0 * LR * H1full
	 + (LF - LQ) * tcga 
	 + 0.5 * (LF - LQ) * (LF - LQ) * tgaga
	 // missing delta terms. Note that tcga does contain the delta term, so we only need to add it once
	 + (LF - LQ) * C1qqdelta(A_F) * tH1F
	 + LF * 2 * Delta2qq * pdf_factor_x1x2
	 + (1. / 6 * A1q * beta0 * pow(LQ, 3) * pdf_factor_x1x2 
	    + 0.5 * pow(LQ, 2) * (A2q * pdf_factor_x1x2 
				  + beta0 * sig11))
	 - LQ * ((B2_A_F + LQ * A2q) * pdf_factor_x1x2 
		 - beta0 * tH1
		 + tgamma2 
		 + 2 * Delta2qq * pdf_factor_x1x2)
	 - 0.5 * LQ * (B1q + 0.5 * A1q * LQ) * (H1full + tH1)
	 //  should already be included via H1full. check!!!
	 //  The term was excluded in the non-distribution calculation (already included via H1full) and should thus probably be erased also at the other places, as done now !!!
	 // 	       + 0.5 * pow(LQ, 2) * (B1q + 0.5 * A1q * LQ) * (B1q + 0.5 * A1q * LQ) * pdf_factor_x1x2
	 - 0.5 * LQ * (LF - LQ) * (B1q + 0.5 * A1q * LQ) * tH1F);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return result_H2;
}
// }}}

// !!! resummation changes must be imported in complete file !!!

#define no_zx ncollinear[ncollinear[i_c].no_endpoint[y_e]].no_pdf
#define no_xx ncollinear[ncollinear[i_c].no_endpoint[0]].no_pdf
#define x_e ncollinear[i_c].x_a
#define y_e ncollinear[i_c].x_b
// {{{ observable_set::determine_splitting_tH1F(phasespace_set & psi)
void observable_set::determine_splitting_tH1F(phasespace_set & psi){
  static Logger logger("observable_set::determine_splitting_tH1F");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH1F_contribution[i_l].clear();}
  //  tH1F: prefactors of splittings, to be multiplied with pdf factors
  for (int i_c = 1; i_c < ncollinear.size(); i_c++){
    logger << LOG_DEBUG_VERBOSE << "ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << endl;
    if (x_e == 0){logger << LOG_FATAL << "Should not happen at 1st order!!!" << endl; exit(1);}

    // type_splitting_full (leg-wise):
    // x_e - y_e     ---  ---

    //   1 -   0   +   0 -  1   // hard process with g from g -> g (g) splitting
    //   2 -   0   +   0 -  2   // hard process with q from q -> q (g) splitting
    //   3 -   0   +   0 -  3   // hard process with g from q -> g (q) splitting
    //   4 -   0   +   0 -  4   // hard process with q from g -> q (qx) splitting       
    
    if (ncollinear[i_c].type_splitting_full[x_e] == 0){logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = 0 -> no splitting type defined !!!" << endl; exit(1);}
    else if (ncollinear[i_c].type_splitting_full[x_e] == 1){ // hard process with g from g -> g (g) splitting
      coll_tH1F_contribution[no_zx].push_back(3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] + Pggreg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      // -log( x1 ) * (pdf_factor_z1x2 * Pggreg(z1)
      // -log( x1 ) * (pdf_factor_z1x2)*3/(1-z1)
      coll_tH1F_contribution[no_xx].push_back(-psi_z_coll[x_e] * 3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(x_pdf[x_e]) + beta0);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1F_contribution[no_xx = " << no_xx << "] += " << coll_tH1F_contribution[no_xx][coll_tH1F_contribution[no_xx].size() - 1] << endl;
      // -log( x1 ) * (- pdf_factor_x1x2 * z1 * 3 / (1 - z1)) - 3 * D0int(x1) * pdf_factor_x1x2
      // + beta0 * pdf_factor_x1x2 (half of this (because added for each leg)!)
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 2){ // hard process with q from q -> q (g) splitting
      coll_tH1F_contribution[no_zx].push_back(Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      // - g_z1 * (pdf_factor_z1x2) * Pqq(z1)  // from H1F
      coll_tH1F_contribution[no_xx].push_back((-psi_z_coll[x_e] * Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - Pqqint(x_pdf[x_e])));
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1F_contribution[no_xx = " << no_xx << "] += " << coll_tH1F_contribution[no_xx][coll_tH1F_contribution[no_xx].size() - 1] << endl;
      // - g_z1 * (- pdf_factor_x1x2 * z1) * Pqq(z1) - pdf_factor_x1x2 * Pqqint(x1)  // from H1F 
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 3){ // hard process with g from q -> g (q) splitting
      coll_tH1F_contribution[no_zx].push_back(Pgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->g(q):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      // -log( x1 ) * (pdf_factor_qx2 * Pgq(z1) )
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 4){ // hard process with q from g -> q (qx) splitting
      coll_tH1F_contribution[no_zx].push_back(Pqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->q(qx): coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      // - g_z1 * Pqg(z1) * (pdf_factor_gx2)  // from H1F
    }
    else {logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << " -> selected splitting type not allowed !!!" << endl; exit(1);}
  }
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    coll_tH1F[i_l] = accumulate(coll_tH1F_contribution[i_l].begin(), coll_tH1F_contribution[i_l].end(), 0.);
    logger << LOG_DEBUG_VERBOSE << "coll_tH1F[" << i_l << "] = " << coll_tH1F[i_l] << endl;
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ observable_set::determine_splitting_tH1F_tH1(phasespace_set & psi)
void observable_set::determine_splitting_tH1F_tH1(phasespace_set & psi){
  static Logger logger("observable_set::determine_splitting_tH1F_tH1");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH1F_contribution[i_l].clear();}
  //  tH1F: prefactors of splittings, to be multiplied with pdf factors
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH1_contribution[i_l].clear();}
  //  tH1: prefactors of splittings, to be multiplied with pdf factors

  for (int i_c = 1; i_c < ncollinear.size(); i_c++){
    logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << "   x_e = " << x_e << endl;
    if (x_e == 0){logger << LOG_DEBUG_VERBOSE << "ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << " does not contribute here." << endl; continue;}

    // type_splitting_full (leg-wise):
    // x_e - y_e     ---  ---

    //  14 -   0   +   0 - 14   // hard process with q from g -> g (g)  -> q (qx) splitting (-> 4  in sig11 calculation)
    // (42 -   0   +   0 - 42   // hard process with q from g -> q (q)  -> q (g)  splitting (part of  14 -  0   +   0 - 14))
    //  11 -   0   +   0 - 11   // hard process with g from g -> g (g)  -> g (g)  splitting (-> 1  in sig11 calculation)
    // (43 -   0   +   0 - 43   // hard process with g from g -> q (qx) -> g (Q)  splitting (part of  11 -  0   +   0 - 11))
    //  23 -   0   +   0 - 23   // hard process with g from q -> q (g)  -> g (Q)  splitting (-> 3  in sig11 calculation)
    // (31 -   0   +   0 - 31   // hard process with g from q -> g (q)  -> g (g)  splitting (part of  23 -  0   +   0 - 23))
    //  22 -   0   +   0 - 22   // hard process with q from q -> q (g)  -> q (g)  splitting (-> 2  in sig11 calculation)

    //  34 -   0   +   0 - 34   // hard process with q from q -> g (q)  -> q (qx) splitting (----  no equivalent)
    //  54 -   0   +   0 - 54   // hard process with q from qx -> ...   -> q      splitting (non-singlet)    (----  no equivalent)
    //   2 -   2                // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) q from q -> q (g)  splitting
    //   2 -   4   +   4 -  2   // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   4 -   4                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting
    //   1 -   1                // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from g -> g (g)  splitting
    //   1 -   3   +   3 -  1   // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   3 -   3                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting

    if (ncollinear[i_c].type_splitting_full[x_e] == 0){logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = 0 -> no splitting type defined !!!" << endl; exit(1);}
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 11 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with q from g -> g (g) splitting
      coll_tH1F_contribution[no_zx].push_back(3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] + Pggreg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      coll_tH1F_contribution[no_xx].push_back(-psi_z_coll[x_e] * 3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(x_pdf[x_e]) + beta0);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1F_contribution[no_xx = " << no_xx << "] += " << coll_tH1F_contribution[no_xx][coll_tH1F_contribution[no_xx].size() - 1] << endl;
      coll_tH1_contribution[no_xx].push_back(.5 * QT_H1_delta);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1_contribution[no_xx = " << no_xx << "] += " << coll_tH1_contribution[no_xx][coll_tH1_contribution[no_xx].size() - 1] << endl;
    }
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 23 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 3 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with g from q -> g (q) splitting
      coll_tH1F_contribution[no_zx].push_back(Pgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->g(q):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      coll_tH1_contribution[no_zx].push_back(Cgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->g(q):  coll_tH1_contribution[no_zx = " << no_zx << "] += " << coll_tH1_contribution[no_zx][coll_tH1_contribution[no_zx].size() - 1] << endl;
    }
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 22 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with q from q -> q (g) splitting
      coll_tH1F_contribution[no_zx].push_back(Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      coll_tH1F_contribution[no_xx].push_back((-psi_z_coll[x_e] * Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - Pqqint(x_pdf[x_e])));
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1F_contribution[no_xx = " << no_xx << "] += " << coll_tH1F_contribution[no_xx][coll_tH1F_contribution[no_xx].size() - 1] << endl;
      coll_tH1_contribution[no_zx].push_back(Cqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1_contribution[no_zx = " << no_zx << "] += " << coll_tH1_contribution[no_zx][coll_tH1_contribution[no_zx].size() - 1] << endl;
      coll_tH1_contribution[no_xx].push_back(.5 * QT_H1_delta);
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1_contribution[no_xx = " << no_xx << "] += " << coll_tH1_contribution[no_xx][coll_tH1_contribution[no_xx].size() - 1] << endl;
    }
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 14 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 4 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with q from g -> q (qx) splitting
      coll_tH1F_contribution[no_zx].push_back(Pqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->q(qx): coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      coll_tH1_contribution[no_zx].push_back(Cqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->q(qx): coll_tH1_contribution[no_zx = " << no_zx << "] += " << coll_tH1_contribution[no_zx][coll_tH1_contribution[no_zx].size() - 1] << endl;
    }
    else {
      //  logger << LOG_DEBUG_VERBOSE << "ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << " does not contribute here." << endl; 
      continue;
    }
  }
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    coll_tH1F[i_l] = accumulate(coll_tH1F_contribution[i_l].begin(), coll_tH1F_contribution[i_l].end(), 0.);
    coll_tH1[i_l] = accumulate(coll_tH1_contribution[i_l].begin(), coll_tH1_contribution[i_l].end(), 0.);
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
#undef no_zx 
#undef no_xx
#undef x_e
#undef y_e

#define no_zz ncollinear[ncollinear[i_c].no_endpoint[3]].no_pdf
#define no_xz ncollinear[ncollinear[i_c].no_endpoint[x_e]].no_pdf
#define no_zx ncollinear[ncollinear[i_c].no_endpoint[y_e]].no_pdf
#define no_xx ncollinear[ncollinear[i_c].no_endpoint[0]].no_pdf
#define x_e ncollinear[i_c].x_a
#define y_e ncollinear[i_c].x_b
// {{{ observable_set::determine_splitting_tgaga_tcga_tgamma2(phasespace_set & psi)
void observable_set::determine_splitting_tgaga_tcga_tgamma2(phasespace_set & psi){
  static Logger logger("observable_set::determine_splitting_tgaga_tcga_tgamma2");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tgaga_contribution[i_l].clear();}
  //  tgaga: prefactors of splittings, to be multiplied with pdf factors
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tcga_contribution[i_l].clear();}
  //  tcga: prefactors of splittings, to be multiplied with pdf factors
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tgamma2_contribution[i_l].clear();}
  //  tgamma2: prefactors of splittings, to be multiplied with pdf factors

  for (int i_c = 1; i_c < ncollinear.size(); i_c++){

    // type_splitting_full (leg-wise):
    // x_e - y_e     ---  ---

    //  14 -   0   +   0 - 14   // hard process with q from g -> g (g)  -> q (qx) splitting (-> 4  in sig11 calculation)
    // (42 -   0   +   0 - 42   // hard process with q from g -> q (q)  -> q (g)  splitting (part of  14 -  0   +   0 - 14))
    //  11 -   0   +   0 - 11   // hard process with g from g -> g (g)  -> g (g)  splitting (-> 1  in sig11 calculation)
    // (43 -   0   +   0 - 43   // hard process with g from g -> q (qx) -> g (Q)  splitting (part of  11 -  0   +   0 - 11))
    //  23 -   0   +   0 - 23   // hard process with g from q -> q (g)  -> g (Q)  splitting (-> 3  in sig11 calculation)
    // (31 -   0   +   0 - 31   // hard process with g from q -> g (q)  -> g (g)  splitting (part of  23 -  0   +   0 - 23))
    //  22 -   0   +   0 - 22   // hard process with q from q -> q (g)  -> q (g)  splitting (-> 2  in sig11 calculation)

    //  34 -   0   +   0 - 34   // hard process with q from q -> g (q)  -> q (qx) splitting (----  no equivalent)
    //  54 -   0   +   0 - 54   // hard process with q from qx -> ...   -> q      splitting (non-singlet)    (----  no equivalent)
    //   2 -   2                // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) q from q -> q (g)  splitting
    //   2 -   4   +   4 -  2   // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   4 -   4                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting
    //   1 -   1                // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from g -> g (g)  splitting
    //   1 -   3   +   3 -  1   // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   3 -   3                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting

    if (ncollinear[i_c].type_splitting_full[x_e] == 0){logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = 0 -> no splitting type defined !!!" << endl; exit(1);}


    //    tgaga
    // 11 - 0
    //   -  log( x1 )*( (pdf_factor_z1x2-pdf_factor_x1x2* z1 )*(D0gggg/(1-z1)+D1gggg* log(1-z1)/(1-z1)) + pdf_factor_z1x2*(Pggggreg(z1, beta0)+Pgqqg(z1, nf)) ) + (Deltagggg-D0gggg*D0int( x1 )-D1gggg*D1int( x1 ))*pdf_factor_x1x2;
    // second leg (automatically included)   -  log( x2 )*( (pdf_factor_x1z2-pdf_factor_x1x2* z2 )*(D0gggg/(1-z2)+D1gggg* log(1-z2)/(1-z2)) + pdf_factor_x1z2*(Pggggreg(z2, beta0)+Pgqqg(z2, nf)) ) + (Deltagggg-D0gggg*D0int( x2 )-D1gggg*D1int( x2 ))*pdf_factor_x1x2;
    // 23 - 0
    //   -  log( x1 ) * (Pgqqq(z1)+Pgggq(z1, beta0)) * pdf_factor_qx2;
    // second leg (automatically included)   -  log( x2 ) * (Pgqqq(z2)+Pgggq(z2, beta0)) * pdf_factor_x1q;
    //   double tx1 = 3.0*log(x1)*z1/(1-z1) - 3*D0int(x1) + beta0,  tz1 = -log(x1)*( 3.0/(1-z1) + Pggreg(z1) ) ;
    //   double tx2 = 3.0*log(x2)*z2/(1-z2) - 3*D0int(x2) + beta0,  tz2 = -log(x2)*( 3.0/(1-z2) + Pggreg(z2) ) ;
    // 1 - 1
    //   2 * (pdf_factor_x1x2*tx1*tx2 + pdf_factor_x1z2*tx1*tz2 + pdf_factor_z1x2*tz1*tx2 + pdf_factor_z1z2*tz1*tz2);
    // 1 - 3
    //   -2 * ( pdf_factor_x1q*tx1 + pdf_factor_z1q*tz1 )*log(x2)*Pgq(z2);
    // second leg (automatically included)      -2 * log(x1)*Pgq(z1)*( pdf_factor_qx2*tx2 + pdf_factor_qz2*tz2 );
    // 3 - 3
    //   2 * log(x1)*Pgq(z1)*log(x2)*Pgq(z2)*pdf_factor_qq;
    

    //   tcga
    // 11 - 0
    //   tcga = tcga - CgqPqg(z1, nf)* log( x1 )     * pdf_factor_z1x2;
    // second leg (automatically included)   tcga = tcga - CgqPqg(z2, nf)* log( x2 )     * pdf_factor_x1z2;
    // 23 - 0
    //   tcga = tcga - CgqPqq(z1)* log( x1 )*pdf_factor_qx2;
    // second leg (automatically included)   tcga = tcga - CgqPqq(z2)* log( x2 )*pdf_factor_x1q;
    // 11 - 0
    //   double diffg10_diffc20 = ( pdf_factor_z1x2*tz1 + pdf_factor_x1x2*tx1 )*C1ggdelta(A_F) ;
    // 3 - 3
    //   double diffg1f_diffc2f = pdf_factor_qq * log(x1)*log(x2)*Pgq(z1)*Cgq(z2) ;
    // 1 - 3
    //   double diffg10_diffc2f = -log(x2)*Cgq(z2)*( pdf_factor_z1q*tz1 + pdf_factor_x1q*tx1 ) ;
    // 23 - 0
    //   double diffg1f_diffc20 = -log(x1)*Pgq(z1)*C1ggdelta(A_F)*pdf_factor_qx2;
    // 0 - 11
    // second leg (automatically included)   double diffc10_diffg20 = C1ggdelta(A_F)*( pdf_factor_x1z2*tz2 + pdf_factor_x1x2*tx2 ) ;
    // 3 - 3
    //   double diffc1f_diffg2f = pdf_factor_qq*log(x1)*log(x2)*Cgq(z1)*Pgq(z2) ;
    // 3 - 1
    // second leg (automatically included)   double diffc1f_diffg20 = -log(x1)*Cgq(z1)*(pdf_factor_qz2*tz2 + pdf_factor_qx2*tx2) ;
    // 0 - 23 
    // second leg (automatically included)   double diffc10_diffg2f = -C1ggdelta(A_F)*log(x2)*Pgq(z2)*pdf_factor_x1q ;
    //   tcga = tcga + diffg10_diffc20 + diffg1f_diffc2f + diffg10_diffc2f + diffg1f_diffc20;
    //   tcga = tcga + diffc10_diffg20 + diffc1f_diffg2f + diffc1f_diffg20 + diffc10_diffg2f;


    //   tgamma2
    //  gamma2: diagonal part
    //   First leg
    // 11 - 0
    //    tgamma2 = tgamma2 - ( 1.5  *Kappa*(  log( x1 )*(pdf_factor_z1x2-pdf_factor_x1x2* z1       )/(1-z1) + D0int( x1 )*pdf_factor_x1x2) +  log( x1 )*P2gg(z1, nf)*pdf_factor_z1x2 );  
    //   Second leg
    // second leg (automatically included)   tgamma2 = tgamma2 - ( 1.5  *Kappa*(  log( x2 )*(pdf_factor_x1z2-pdf_factor_x1x2* z2       )/(1-z2) + D0int( x2 )*pdf_factor_x1x2) +  log( x2 )*P2gg(z2, nf)*pdf_factor_x1z2 );
    //  gamma2: qg channel
    //   First leg
    // 23 - 0
    //   tgamma2 = tgamma2 -  log( x1 )*P2gq(z1, nf)* pdf_factor_qx2;
    //  Second leg
    // second leg (automatically included)   tgamma2 = tgamma2 -  log( x2 )*P2gq(z2, nf)* pdf_factor_x1q;


    else if (ncollinear[i_c].type_splitting_full[x_e] == 11 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with g from g -> g (g) -> g (g)  or  g -> Q (Qx) -> g (Q) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with g from g -> g (g) -> g (g)  or  g -> Q (Qx) -> g (Q) splitting" << endl;
      double temp_Dxgggg = (D0gggg / (1. - psi_z_coll[x_e]) + D1gggg * log(1. - psi_z_coll[x_e]) / (1. - psi_z_coll[x_e]));
      coll_tgaga_contribution[no_zx].push_back(temp_Dxgggg / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgaga  [no_zx = " << setw(2) << no_zx << "] += " << coll_tgaga_contribution[no_zx][coll_tgaga_contribution[no_zx].size() - 1] << endl;
      coll_tgaga_contribution[no_zx].push_back((Pggggreg(psi_z_coll[x_e], beta0) + Pgqqg(psi_z_coll[x_e], N_f)) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(gg/qqx):  tgaga  [no_zx = " << setw(2) << no_zx << "] += " << coll_tgaga_contribution[no_zx][coll_tgaga_contribution[no_zx].size() - 1] << endl;
      coll_tgaga_contribution[no_xx].push_back(-psi_z_coll[x_e] * temp_Dxgggg / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgaga  [no_xx = " << setw(2) << no_xx << "] += " << coll_tgaga_contribution[no_xx][coll_tgaga_contribution[no_xx].size() - 1] << endl;
      coll_tgaga_contribution[no_xx].push_back(Deltagggg - D0gggg * D0int(psi_x_pdf[x_e]) - D1gggg * D1int(psi_x_pdf[x_e]));
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgaga  [no_xx = " << setw(2) << no_xx << "] += " << coll_tgaga_contribution[no_xx][coll_tgaga_contribution[no_xx].size() - 1] << endl;

      coll_tcga_contribution[no_zx].push_back(CgqPqg(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(qqx):     tcga   [no_zx = " << setw(2) << no_zx << "] += " << coll_tcga_contribution[no_zx][coll_tcga_contribution[no_zx].size() - 1] << endl;

      double temp_C1ggdelta_A_F = C1ggdelta(A_F);
      double tx1 = -3 * psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(psi_x_pdf[x_e]) + beta0;
      double tz1 = (3 / (1. - psi_z_coll[x_e]) + Pggreg(psi_z_coll[x_e])) / psi_g_z_coll[x_e];
      coll_tcga_contribution[no_zx].push_back(tz1 * temp_C1ggdelta_A_F); // scale dependent
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tcga   [no_zx = " << setw(2) << no_zx << "] += " << coll_tcga_contribution[no_zx][coll_tcga_contribution[no_zx].size() - 1] << endl;
      coll_tcga_contribution[no_xx].push_back(tx1 * temp_C1ggdelta_A_F); // scale dependent
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tcga   [no_xx = " << setw(2) << no_xx << "] += " << coll_tcga_contribution[no_xx][coll_tcga_contribution[no_xx].size() - 1] << endl;

      coll_tgamma2_contribution[no_zx].push_back(1.5 * Kappa / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] );
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgamma2[no_zx = " << setw(2) << no_zx << "] += " << coll_tgamma2_contribution[no_zx][coll_tgamma2_contribution[no_zx].size() - 1] << endl;
      coll_tgamma2_contribution[no_zx].push_back(P2gg(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgamma2[no_zx = " << setw(2) << no_zx << "] += " << coll_tgamma2_contribution[no_zx][coll_tgamma2_contribution[no_zx].size() - 1] << endl;
      coll_tgamma2_contribution[no_xx].push_back(-1.5 * Kappa * psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgamma2[no_xx = " << setw(2) << no_xx << "] += " << coll_tgamma2_contribution[no_xx][coll_tgamma2_contribution[no_xx].size() - 1] << endl;
      coll_tgamma2_contribution[no_xx].push_back(-1.5 * Kappa * D0int(psi_x_pdf[x_e]));
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgamma2[no_xx = " << setw(2) << no_xx << "] += " << coll_tgamma2_contribution[no_xx][coll_tgamma2_contribution[no_xx].size() - 1] << endl;
    }

    else if (ncollinear[i_c].type_splitting_full[x_e] == 23 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with g from Q -> Q (g) -> g (g)  or  Q -> g (Q) -> g (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with g from Q -> Q (g) -> g (g)  or  Q -> g (Q) -> g (g) splitting" << endl;
      coll_tgaga_contribution[no_zx].push_back((Pgqqq(psi_z_coll[x_e]) + Pgggq(psi_z_coll[x_e], beta0)) / psi_g_z_coll[x_e]);
 
      coll_tcga_contribution[no_zx].push_back(CgqPqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      // reactivate this term: !!!
      coll_tcga_contribution[no_zx].push_back(Pgq(psi_z_coll[x_e]) * C1ggdelta(A_F) / psi_g_z_coll[x_e]); // scale dependent
      // until here.

      coll_tgamma2_contribution[no_zx].push_back(P2gq(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 1){ 
      // hard process with both legs x_e and y_e with g from g -> g (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << " with g from g -> g (g) splitting" << endl;

      double tx1 = -3 * psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(psi_x_pdf[x_e]) + beta0;
      double tz1 = (3 / (1. - psi_z_coll[x_e]) + Pggreg(psi_z_coll[x_e])) / psi_g_z_coll[x_e];
      double tx2 = -3 * psi_z_coll[y_e] / (1. - psi_z_coll[y_e]) / psi_g_z_coll[y_e] - 3 * D0int(psi_x_pdf[y_e]) + beta0;
      double tz2 = (3 / (1. - psi_z_coll[y_e]) + Pggreg(psi_z_coll[y_e])) / psi_g_z_coll[y_e];

      coll_tgaga_contribution[no_zz].push_back(2 * tz1 * tz2);
      coll_tgaga_contribution[no_zx].push_back(2 * tz1 * tx2);
      coll_tgaga_contribution[no_xz].push_back(2 * tx1 * tz2);
      coll_tgaga_contribution[no_xx].push_back(2 * tx1 * tx2);
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 3 && ncollinear[i_c].type_splitting_full[y_e] == 3){
      // hard process with both legs x_e and y_e with g from q -> g (q) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << "with g from q -> g (q) splitting" << endl;
      coll_tgaga_contribution[no_zz].push_back(2 * Pgq(psi_z_coll[x_e]) * Pgq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);

      coll_tcga_contribution[no_zz].push_back(Pgq(psi_z_coll[x_e]) * Cgq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_zz].push_back(Cgq(psi_z_coll[x_e]) * Pgq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);

      /*
      coll_tgamma2_contribution[no_zz].push_back();
      */
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 3){
      // hard process with leg x_e with g from g -> g (g) splitting and leg y_e with g from q -> g (q) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with leg " << x_e << " with g from g -> g (g) splitting and leg " << y_e << " with g from q -> g (q) splitting" << endl;
      double tx1 = -3 * psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(psi_x_pdf[x_e]) + beta0;
      double tz1 = (3 / (1. - psi_z_coll[x_e]) + Pggreg(psi_z_coll[x_e])) / psi_g_z_coll[x_e];
      coll_tgaga_contribution[no_zz].push_back(-2 * tz1 * log(psi_x_pdf[y_e]) * Pgq(psi_z_coll[y_e]));
      coll_tgaga_contribution[no_xz].push_back(-2 * tx1 * log(psi_x_pdf[y_e]) * Pgq(psi_z_coll[y_e]));

      coll_tcga_contribution[no_zz].push_back(Cgq(psi_z_coll[y_e]) * tz1 / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_xz].push_back(Cgq(psi_z_coll[y_e]) * tx1 / psi_g_z_coll[y_e]);

      /*
      coll_tgamma2_contribution[no_zz].push_back();
      coll_tgamma2_contribution[no_xz].push_back();
      */
    }

    else if (ncollinear[i_c].type_splitting_full[x_e] == 14 && ncollinear[i_c].type_splitting_full[y_e] == 0){ 
      // hard process with q from g -> g (g) -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from g -> g (g) -> q (qx) splitting" << endl;
      coll_tgaga_contribution[no_zx].push_back(1. / psi_g_z_coll[x_e] * (Pqqqg(psi_z_coll[x_e]) + Pqggg(psi_z_coll[x_e], beta0)));
      coll_tcga_contribution[no_zx].push_back((CqqPqg(psi_z_coll[x_e]) + CqgPgg(psi_z_coll[x_e], beta0)) / psi_g_z_coll[x_e]);
      coll_tgamma2_contribution[no_zx].push_back(P2qg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 34 && ncollinear[i_c].type_splitting_full[y_e] == 0){ 
      // hard process with q from Q -> g (Q) -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from Q -> g (Q) -> q (qx) splitting" << endl;
      coll_tgaga_contribution[no_zx].push_back(1. / psi_g_z_coll[x_e] * Pqggq(psi_z_coll[x_e]));
      coll_tcga_contribution[no_zx].push_back(CqgPgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tgamma2_contribution[no_zx].push_back(P2qqS(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 22 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with q from q -> q (g) -> q (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from q -> q (g) -> q (g) splitting" << endl;
      double temp_Dqqqq = (D0qqqq / (1. - psi_z_coll[x_e]) + D1qqqq * log(1. - psi_z_coll[x_e]) / (1. - psi_z_coll[x_e])) / psi_g_z_coll[x_e];
      coll_tgaga_contribution[no_zx].push_back(temp_Dqqqq);
      coll_tgaga_contribution[no_zx].push_back(Pqqqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tgaga_contribution[no_xx].push_back(-psi_z_coll[x_e] * temp_Dqqqq);
      coll_tgaga_contribution[no_xx].push_back((Deltaqqqq - D0qqqq * D0int(psi_x_pdf[x_e]) - D1qqqq * D1int(psi_x_pdf[x_e])));
      coll_tcga_contribution[no_zx].push_back(CqqPqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tgamma2_contribution[no_zx].push_back(P2qqV(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      coll_tgamma2_contribution[no_zx].push_back(2.0 / 3 * Kappa * 1. / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tgamma2_contribution[no_xx].push_back(-2.0 / 3 * Kappa * (psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] + D0int(psi_x_pdf[x_e])));
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 54 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with q from qx -> ... -> q  splitting (non-singlet)
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from qx -> ... -> q  splitting (non-singlet)" << endl;
      coll_tgamma2_contribution[no_zx].push_back(P2qqbV(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 2){ 
      // hard process with both legs x_e and y_e with q from q -> q (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << " with q from q -> q (g) splitting" << endl;
      double temp_Pqq_zz = 2 * Pqq(psi_z_coll[x_e]) * Pqq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e];
      coll_tgaga_contribution[no_zz].push_back(temp_Pqq_zz);
      coll_tgaga_contribution[no_xz].push_back(-2 * Pqq(psi_z_coll[y_e]) * Pqqint(psi_x_pdf[x_e]) / psi_g_z_coll[y_e]);
      coll_tgaga_contribution[no_xz].push_back(-temp_Pqq_zz * psi_z_coll[x_e]);
      coll_tgaga_contribution[no_zx].push_back(-2 * Pqq(psi_z_coll[x_e]) * Pqqint(psi_x_pdf[y_e]) / psi_g_z_coll[x_e]);
      coll_tgaga_contribution[no_zx].push_back(-temp_Pqq_zz * psi_z_coll[y_e]);
      coll_tgaga_contribution[no_xx].push_back(2 * Pqqint(psi_x_pdf[x_e]) * Pqqint(psi_x_pdf[y_e]));
      coll_tgaga_contribution[no_xx].push_back(2 * Pqq(psi_z_coll[x_e]) * Pqqint(psi_x_pdf[y_e]) / psi_g_z_coll[x_e] * psi_z_coll[x_e]);
      coll_tgaga_contribution[no_xx].push_back(2 * Pqq(psi_z_coll[y_e]) * Pqqint(psi_x_pdf[x_e]) / psi_g_z_coll[y_e] * psi_z_coll[y_e] );
      coll_tgaga_contribution[no_xx].push_back(temp_Pqq_zz * psi_z_coll[x_e] * psi_z_coll[y_e]);
      double temp_C1qqdelta_A_F = C1qqdelta(A_F);
      //      double temp_C1qqdelta_A_F = C1qqdelta(value_A_F[sd][ss]);
      coll_tcga_contribution[no_zz].push_back((Cqq(psi_z_coll[x_e]) * Pqq(psi_z_coll[y_e]) + Pqq(psi_z_coll[x_e]) * Cqq(psi_z_coll[y_e])) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_xz].push_back(temp_C1qqdelta_A_F * Pqq(psi_z_coll[y_e]) / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_xz].push_back(-Cqq(psi_z_coll[y_e]) / psi_g_z_coll[y_e] * (Pqq(psi_z_coll[x_e]) * psi_z_coll[x_e] / psi_g_z_coll[x_e] + Pqqint(psi_x_pdf[x_e])));
      coll_tcga_contribution[no_zx].push_back(Pqq(psi_z_coll[x_e]) * temp_C1qqdelta_A_F / psi_g_z_coll[x_e]);
      coll_tcga_contribution[no_zx].push_back(-Cqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e] * (Pqq(psi_z_coll[y_e]) * psi_z_coll[y_e] / psi_g_z_coll[y_e] + Pqqint(psi_x_pdf[y_e])));
      coll_tcga_contribution[no_xx].push_back(-temp_C1qqdelta_A_F * (psi_z_coll[x_e] * Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + Pqqint(psi_x_pdf[x_e])));
      coll_tcga_contribution[no_xx].push_back(-temp_C1qqdelta_A_F * (psi_z_coll[y_e] * Pqq(psi_z_coll[y_e]) / psi_g_z_coll[y_e] + Pqqint(psi_x_pdf[y_e])));


    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 4 && ncollinear[i_c].type_splitting_full[y_e] == 4){
      // hard process with both legs x_e and y_e with q from g -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << "with q from g -> q (qx) splitting" << endl;
      coll_tgaga_contribution[no_zz].push_back(2 * Pqg(psi_z_coll[x_e]) * Pqg(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_zz].push_back((Pqg(psi_z_coll[x_e]) * Cqg(psi_z_coll[y_e]) + Cqg(psi_z_coll[x_e]) * Pqg(psi_z_coll[y_e])) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 4){
      // hard process with leg x_e with q from q -> q (g) splitting and leg y_e with q from g -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with leg " << x_e << " with q from q -> q (g) splitting and leg " << y_e << " with q from g -> q (qx) splitting" << endl;
      coll_tgaga_contribution[no_zz].push_back(2 * Pqg(psi_z_coll[y_e]) / psi_g_z_coll[y_e] * Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tgaga_contribution[no_xz].push_back(-2 * Pqg(psi_z_coll[y_e]) / psi_g_z_coll[y_e] * (Pqq(psi_z_coll[x_e]) * psi_z_coll[x_e] / psi_g_z_coll[x_e] + Pqqint(psi_x_pdf[x_e])));
      coll_tcga_contribution[no_zz].push_back((Pqq(psi_z_coll[x_e]) * Cqg(psi_z_coll[y_e]) + Cqq(psi_z_coll[x_e]) * Pqg(psi_z_coll[y_e])) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_xz].push_back(-(Pqq(psi_z_coll[x_e]) * psi_z_coll[x_e] / psi_g_z_coll[x_e] + Pqqint(psi_x_pdf[x_e])) * Cqg(psi_z_coll[y_e]) / psi_g_z_coll[y_e]);
      //      coll_tcga_contribution[no_xz].push_back(C1qqdelta(value_A_F[sd][ss]) * Pqg(psi_z_coll[y_e]) / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_xz].push_back(C1qqdelta(A_F) * Pqg(psi_z_coll[y_e]) / psi_g_z_coll[y_e]);
    }
    else {logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << "   ncollinear[" << i_c << "].type_splitting_full[" << y_e << "] = " << ncollinear[i_c].type_splitting_full[y_e] << " -> selected splitting type not allowed !!!" << endl; exit(1);}
  }
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    for (int i_c = 0; i_c < coll_tcga_contribution[i_l].size(); i_c++){
      logger << LOG_DEBUG_VERBOSE << "tcga[" << i_l << "][" << i_c << "] = " << coll_tcga_contribution[i_l][i_c] << endl;
    }

    coll_tgaga[i_l] = accumulate(coll_tgaga_contribution[i_l].begin(), coll_tgaga_contribution[i_l].end(), 0.);
    coll_tcga[i_l] = accumulate(coll_tcga_contribution[i_l].begin(), coll_tcga_contribution[i_l].end(), 0.);
    coll_tgamma2[i_l] = accumulate(coll_tgamma2_contribution[i_l].begin(), coll_tgamma2_contribution[i_l].end(), 0.);
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
// {{{ observable_set::determine_splitting_tH2(phasespace_set & psi)
void observable_set::determine_splitting_tH2(phasespace_set & psi){
  static Logger logger("observable_set::determine_splitting_tgaga_tcga_tgamma2");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH2_contribution[i_l].clear();}

  ///  static double H2qqD0 = -404. / 27 + (56. * N_f) / 81 + 14. * zeta3;

  //  coll: prefactors of splittings, to be multiplied with pdf factors
  //  coll[0] -> tH2
  for (int i_c = 1; i_c < ncollinear.size(); i_c++){
    // type_splitting_full (leg-wise):
    // x_e y_e      --- ---

    // 14 -  0   +   0 - 14   // hard process with q from g -> g (g) -> q (qx) splitting             (-> 4  in sig11 calculation)
    //(42 -  0   +   0 - 42   // hard process with q from g -> q (qx) -> q (g) splitting             (part of  14 -  0   +   0 - 14))

    // 34 -  0   +   0 - 34   // hard process with q from Q -> g (Q) -> q (qx) splitting             (----  no equivalent
    // 54 -  0   +   0 - 54   // hard process with q from qx -> ... -> q  splitting (non-singlet)    (----  no equivalent

    // 11 -  0   +   0 - 11   // hard process with g from g -> g (g) -> g (g) splitting              (-> 1  in sig11 calculation)
    //(43 -  0   +   0 - 43   // hard process with g from g -> Q (Qx) -> g (Q) splitting             (part of  11 -  0   +   0 - 11))

    // 23 -  0   +   0 - 23   // hard process with g from Q -> Q (g) -> g (Q) splitting              (-> 3  in sig11 calculation)
    //(31 -  0   +   0 - 31   // hard process with g from Q -> g (Q) -> g (g) splitting              (part of  23 -  0   +   0 - 23))

    // 22 -  0   +   0 - 22   // hard process with q from q -> q (g) -> q (g) splitting              (-> 2  in sig11 calculation)

    //  2 -  2                // hard process with (leg 1) q from q -> q (g) splitting and (leg 2) q from q -> q (g) splitting
    //  2 -  4   +   4 -  2   // hard process with (leg 1) q from q -> q (g) splitting and (leg 2) g from q -> g (q) splitting
    //  4 -  4                // hard process with (leg 1) g from q -> g (q) splitting and (leg 2) g from q -> g (q) splitting





    //  gg channel  
    // 11 - 0 (with a factor .5) or 1 - 1 (without !!!)
    // analogously to qqbar channel: 1 - 1
    //  tH2 = tH2 + H2_delta*pdf_factor_x1x2;
    // 11 - 0
    //  tH2 = tH2 - 0.5*( log(x1)*(pdf_factor_z1x2-pdf_factor_x1x2*z1)*H2ggD0/(1-z1) + H2ggD0*D0int(x1)*pdf_factor_x1x2 );
    // second leg (automatically included)  tH2 = tH2 - 0.5*( log(x2)*(pdf_factor_x1z2-pdf_factor_x1x2*z2)*H2ggD0/(1-z2) + H2ggD0*D0int(x2)*pdf_factor_x1x2 );
    // 11 - 0
    //  tH2 = tH2 - log(x1)*pdf_factor_z1x2*C2ggreg(z1,nf); 
    // second leg (automatically included)  tH2 = tH2 - log(x2)*pdf_factor_x1z2*C2ggreg(z2,nf);
    // 1 - 1
    //  tH2 = tH2 + log(x1)*log(x2)*Ggg(z1)*Ggg(z2)*pdf_factor_z1z2;
    //  qg channel  
    // 23 - 0
    //  tH2 = tH2 - log(x1)*pdf_factor_qx2*C2gq(z1,nf);
    // second leg (automatically included)  tH2 = tH2 - log(x2)*pdf_factor_x1q*C2gq(z2,nf);
    //  tH2 = tH2 - log(x1)*H1_delta*pdf_factor_qx2*Cgq(z1);          // doesn't appear explicitly in HHNLO
    // second leg (automatically included)  tH2 = tH2 - log(x2)*H1_delta*pdf_factor_x1q*Cgq(z2);          // doesn't appear explicitly in HHNLO 
    // 1 - 3 (both needed)
    //  tH2 = tH2 + log(x1)*log(x2)*Ggg(z1)*Ggq(z2)*pdf_factor_z1q;
    //  tH2 = tH2 + log(x1)*log(x2)*Ggq(z1)*Ggg(z2)*pdf_factor_qz2;
    //  qq channel (C1C1 and G1G1)
    // 3 - 3
    //  tH2 = tH2 + log(x1)*log(x2)*Cgq(z1)*Cgq(z2)*pdf_factor_qq;    
    //  tH2 = tH2 + log(x1)*log(x2)*Ggq(z1)*Ggq(z2)*pdf_factor_qq;


    //    cout << ncollinear[i_c].type_splitting_full[x_e] << "   " << ncollinear[i_c].type_splitting_full[y_e] << "   no_xx = " << no_xx << "   no_zx = " << no_zx << "   no_xz = " << no_xz << "   no_zz = " << no_zz << endl;   


    if (ncollinear[i_c].type_splitting_full[x_e] == 0){logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = 0 -> no splitting type defined !!!" << endl; exit(1);}
    else if (ncollinear[i_c].type_splitting_full[x_e] == 11 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with g from g -> g (g) -> g (g)  or  g -> Q (Qx) -> g (Q) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with g from g -> g (g) -> g (g)  or  g -> Q (Qx) -> g (Q) splitting" << endl;
      coll_tH2_contribution[no_zx].push_back(.5 * H2ggD0 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_zx].push_back(C2ggreg(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_xx].push_back(-.5 * psi_z_coll[x_e] * H2ggD0 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_xx].push_back(-.5* H2ggD0 * D0int(psi_x_pdf[x_e]));
    }

    else if (ncollinear[i_c].type_splitting_full[x_e] == 23 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with g from Q -> Q (g) -> g (g)  or  Q -> g (Q) -> g (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with g from Q -> Q (g) -> g (g)  or  Q -> g (Q) -> g (g) splitting" << endl;
      coll_tH2_contribution[no_zx].push_back(C2gq(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_zx].push_back(QT_H1_delta * Cgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 1){ 
      // hard process with both legs x_e and y_e with g from g -> g (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << " with g from g -> g (g) splitting" << endl;
      coll_tH2_contribution[no_zz].push_back(Ggg(psi_z_coll[x_e]) * Ggg(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tH2_contribution[no_xx].push_back(QT_H2_delta);
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 3 && ncollinear[i_c].type_splitting_full[y_e] == 3){
      // hard process with both legs x_e and y_e with g from q -> g (q) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << "with g from q -> g (q) splitting" << endl;
      coll_tH2_contribution[no_zz].push_back(Cgq(psi_z_coll[x_e]) * Cgq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tH2_contribution[no_zz].push_back(Ggq(psi_z_coll[x_e]) * Ggq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 3){
      // hard process with leg x_e with g from g -> g (g) splitting and leg y_e with g from q -> g (q) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with leg " << x_e << " with g from g -> g (g) splitting and leg " << y_e << " with g from q -> g (q) splitting" << endl;
      coll_tH2_contribution[no_zz].push_back(Ggg(psi_z_coll[x_e]) * Ggq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      // first leg is gg, second is gq (called twice!!!)      coll_tH2_contribution[no_zz].push_back(Ggq(psi_z_coll[x_e]) * Ggg(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }



    else if (ncollinear[i_c].type_splitting_full[x_e] == 14 && ncollinear[i_c].type_splitting_full[y_e] == 0){ 
      // hard process with q from g -> g (g) -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from g -> g (g) -> q (qx) splitting" << endl;
      coll_tH2_contribution[no_zx].push_back(C2qg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      // no corresponding term (23-0) proportional to QT_H1_delta !!!
      // could be the one that is added outside of tH2 ???
      // new term, shifted here from outside this function:
      coll_tH2_contribution[no_zx].push_back(QT_H1_delta * Cqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 34 && ncollinear[i_c].type_splitting_full[y_e] == 0){ 
      // hard process with q from Q -> g (Q) -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from Q -> g (Q) -> q (qx) splitting" << endl;
      // this contribution should be subtracted again for the specified quarks q and qx, which are treated in 22 and 54, respectively !!!
      coll_tH2_contribution[no_zx].push_back(C2qqp(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 22 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with q from q -> q (g) -> q (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from q -> q (g) -> q (g) splitting" << endl;
      coll_tH2_contribution[no_zx].push_back(.5 * H2qqD0 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_zx].push_back(C2qqreg(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      // this contribution should be subtracted again for the specified quarks q and qx, which are treated in 22 and 54, respectively (from 34)!!!
      coll_tH2_contribution[no_zx].push_back(-C2qqp(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_xx].push_back(-.5 * H2qqD0 * psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_xx].push_back(-.5 * H2qqD0 * D0int(psi_x_pdf[x_e]));
      // new term, shifted here from outside this function:
      coll_tH2_contribution[no_zx].push_back(QT_H1_delta * Cqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 54 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with q from qx -> ... -> q  splitting (non-singlet)
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from qx -> ... -> q  splitting (non-singlet)" << endl;
      coll_tH2_contribution[no_zx].push_back(C2qqb(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      // this contribution should be subtracted again for the specified quarks q and qx, which are treated in 22 and 54, respectively (from 34)!!!
      coll_tH2_contribution[no_zx].push_back(-C2qqp(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 2){ 
      // hard process with both legs x_e and y_e with q from q -> q (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << " with q from q -> q (g) splitting" << endl;
      coll_tH2_contribution[no_zz].push_back(Cqq(psi_z_coll[x_e]) * Cqq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tH2_contribution[no_xx].push_back(QT_H2_delta);
    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 4 && ncollinear[i_c].type_splitting_full[y_e] == 4){
      // hard process with both legs x_e and y_e with q from g -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << "with q from g -> q (qx) splitting" << endl;
      coll_tH2_contribution[no_zz].push_back(Cqg(psi_z_coll[x_e]) * Cqg(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 4){
      // hard process with leg x_e with q from q -> q (g) splitting and leg y_e with q from g -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with leg " << x_e << " with q from q -> q (g) splitting and leg " << y_e << " with q from g -> q (qx) splitting" << endl;
      coll_tH2_contribution[no_zz].push_back(Cqq(psi_z_coll[x_e]) * Cqg(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
    }
    else {logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << "   ncollinear[" << i_c << "].type_splitting_full[" << y_e << "] = " << ncollinear[i_c].type_splitting_full[y_e] << " -> selected splitting type not allowed !!!" << endl; exit(1);}
  }
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    coll_tH2[i_l] = accumulate(coll_tH2_contribution[i_l].begin(), coll_tH2_contribution[i_l].end(), 0.);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
#undef no_zz 
#undef no_xz
#undef no_zx 
#undef no_xx
#undef x_e
#undef y_e
 // {{{ observable_set::determine_psp_weight_QT(phasespace_set & psi)
void observable_set::determine_psp_weight_QT(phasespace_set & psi){
  static Logger logger("observable_set::determine_psp_weight_QT");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  this_psp_weight = integrand;
  this_psp_weight2 = pow(this_psp_weight, 2.);
  step_sum_weight += this_psp_weight;
  step_sum_weight2 += this_psp_weight2;

  if (switch_CV){
    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	this_psp_weight_CV[i_q][i_s] = integrand_qTcut_CV[i_q][i_s];
	this_psp_weight2_CV[i_q][i_s] = pow(this_psp_weight_CV[i_q][i_s], 2.);
	step_sum_weight_CV[i_q][i_s] += this_psp_weight_CV[i_q][i_s];
	step_sum_weight2_CV[i_q][i_s] += this_psp_weight2_CV[i_q][i_s];
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}

//-------------------------------------------------------------------------------------------------------------------------------  
//-------------------------------------------------------------------------------------------------------------------------------
//#include "OV.observable.qTsubtraction.cpp"

// {{{ observable_set::determine_correlationoperator_QCD()
void observable_set::determine_correlationoperator_QCD(){
  static Logger logger("observable_set::CS_determine_correlationoperator_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static map<int, double> charge_particle;
  if (initialization == 1){
    fill_charge_particle(charge_particle);
    initialization = 0;
  }
  vector<string> pa_name(csi->type_parton[0].size(), "");
  if (csi->type_parton[0][1] > -10 && csi->type_parton[0][1] < 10){pa_name[1] = "a";}
  if (csi->type_parton[0][2] > -10 && csi->type_parton[0][2] < 10){pa_name[2] = "b";}
  int count = 0;
  vector<string> alphabet(csi->type_parton[0].size() - 3, "");
  for (int i_p = 0; i_p < alphabet.size(); i_p++){alphabet[i_p] = char(105 + i_p);}
  for (int i_p = 3; i_p < pa_name.size(); i_p++){if (csi->type_parton[0][i_p] > -10 && csi->type_parton[0][i_p] < 10){pa_name[i_p] = alphabet[count++];}}

  int temp_type_correction = 1;

  for (int temp_no_emitter = 1; temp_no_emitter < csi->type_parton[0].size(); temp_no_emitter++){
    for (int temp_no_spectator = temp_no_emitter; temp_no_spectator < csi->type_parton[0].size(); temp_no_spectator++){
      logger << LOG_DEBUG << "temp_no_emitter = " << temp_no_emitter << "   temp_no_spectator = " << temp_no_spectator << endl;
      //    for (int temp_no_spectator = 1; temp_no_spectator < csi->type_parton[0].size(); temp_no_spectator++){
      //      if (temp_no_emitter == temp_no_spectator){continue;}
      if (pa_name[temp_no_emitter] == "" || pa_name[temp_no_spectator] == ""){continue;}
      //      if (temp_no_emitter < 3 && temp_no_spectator < 3){continue;}
      //      if (csi.type_contribution == "VT" && (temp_no_emitter < 3 || temp_no_spectator < 3)){continue;}

      vector<int> temp_pair(2);
      temp_pair[0] = std::min(temp_no_emitter, temp_no_spectator);
      temp_pair[1] = std::max(temp_no_emitter, temp_no_spectator);

      int temp_type;
      if (csi->type_parton[0][temp_no_emitter] == 0){temp_type = 0;}
      else {temp_type = 1;}

      int temp_massive;
      if (mass_parton[0][temp_no_emitter] == 0. && mass_parton[0][temp_no_spectator] == 0.){temp_massive = 0;}
      else if (mass_parton[0][temp_no_emitter] != 0. && mass_parton[0][temp_no_spectator] == 0.){temp_massive = 1;}
      else if (mass_parton[0][temp_no_emitter] == 0. && mass_parton[0][temp_no_spectator] != 0.){temp_massive = 2;}
      else if (mass_parton[0][temp_no_emitter] != 0. && mass_parton[0][temp_no_spectator] != 0.){temp_massive = 3;}
      else {
	cout << "Should not happen!" << endl;
	cout << "mass_parton[0][temp_no_emitter = " << temp_no_emitter << "] = " << mass_parton[0][temp_no_emitter] << endl;
	cout << "mass_parton[0][temp_no_spectator = " << temp_no_spectator << "] = " << mass_parton[0][temp_no_spectator] << endl;
      }

      int temp_type_combination = 0;
      if ((temp_no_emitter > 2 && temp_no_spectator > 2) && temp_no_emitter == temp_no_spectator){temp_type_combination = 1;}
      else if (temp_no_emitter > 2 && temp_no_spectator > 2){temp_type_combination = 2;}
      else if (temp_no_emitter < 3 && temp_no_spectator > 2){temp_type_combination = 3;}
      //      else if ((temp_no_emitter < 3 && temp_no_spectator < 3) && temp_no_emitter == temp_no_spectator){temp_type_combination = 4;}
      else if ((temp_no_emitter < 3 && temp_no_spectator < 3)){temp_type_combination = 5;}
      else {logger << LOG_INFO << "Combination should not appear." << endl; exit(1);}

      double temp_charge_factor = 1.;
      if (temp_type_combination == 1){// || temp_type_combination == 4){
	if (csi->type_parton[0][temp_no_emitter] == 0){temp_charge_factor = C_A;}
	else {temp_charge_factor = C_F;}
      }

      string temp_name;
      temp_name = "C^{" + pa_name[temp_no_emitter] + "," + pa_name[temp_no_spectator] + "}";
      QT_correlationoperator.push_back(correlationoperator(temp_name, temp_type, temp_pair, 0, csi->type_parton[0], temp_charge_factor, temp_no_emitter, temp_no_spectator, temp_type_correction, temp_type_combination, temp_massive));
    }
  }

  for (int i_a = 0; i_a < QT_correlationoperator.size(); i_a++){
    logger << LOG_DEBUG << "QT_correlationoperator[" << i_a << "].name() = " << QT_correlationoperator[i_a].name << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
// }}}
