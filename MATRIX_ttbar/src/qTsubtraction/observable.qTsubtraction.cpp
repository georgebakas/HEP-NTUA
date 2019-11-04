#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

//  important issues:
// - check which version in H2 is correct (potentially doubled term in VT2 contribution, which is already included via H1full -> relevant if LQ =/= 0) !!!
// - check if 'scale-dependent' parts are really scale dependent !!!
// - check missing LQ contributions in gg channel !!!

void observable_set::calculate_A_F(){
  static Logger logger("observable_set::calculate_A_F");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (initial_gg){
    A_F = 2 * (QT_H1_delta - pi2_6 * C_A);
  }
  else if (initial_qqx){
    A_F = 2 * (QT_H1_delta - pi2_6 * C_F);
  }
  else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_B2_A_F(){
  static Logger logger("observable_set::calculate_B2_A_F");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (initial_gg){
    B2_A_F = B2g + .5 * beta0 * A_F;
  }
  else if (initial_qqx){
    B2_A_F = B2q + .5 * beta0 * A_F;
  }
  else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

double observable_set::calculate_sigma12(double & pdf_factor_x1x2){
  static Logger logger("observable_set::calculate_sigma12");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double sig12 = 0.;
  if (initial_gg){
    sig12 = -0.5 * A1g * pdf_factor_x1x2;
  }
  else if (initial_qqx){
    sig12 = -0.5 * A1q * pdf_factor_x1x2;
  }
  else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig12;
}

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
  else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}

  if (QT_finalstate_massive_coloured){
    sig11 += pdf_factor_x1x2 * Gamma1born;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig11;
}

double observable_set::calculate_sigma24(double & pdf_factor_x1x2){
  static Logger logger("observable_set::calculate_sigma24");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double sig24 = 0.;
  if (initial_gg){
    sig24 = 1. / 8 * A1g * A1g * pdf_factor_x1x2;
  }
  else if (initial_qqx){
    sig24 = 1. / 8 * A1q * A1q * pdf_factor_x1x2;
  }
  else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig24;
}

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
  else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig23;
}

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
  else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}

  if (QT_finalstate_massive_coloured){
    // new implementation (should also work for simultaneous calculation of channels):
    //    if (initial_diag_gg || initial_diag_qqx || initial_pdf_gq){ // <-  should be identical to the next version
    if (initial_diag || initial_pdf_gq){
      sig22 -= Gamma1born * tH1F;
    }
    if (initial_diag_gg){
      sig22 += - 0.5 * B1g * Gamma1born * pdf_factor_x1x2;
    }
    if (initial_diag_qqx){
      sig22 += - 0.5 * B1q * Gamma1born * pdf_factor_x1x2;
    }
    if (initial_diag){
      sig22 += + 0.5 * Gamma1squared * pdf_factor_x1x2;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig22;
}

// new version with additional argument (as for calculate_H1 ...) !!!
double observable_set::calculate_sigma21(double & pdf_factor_x1x2, double & sig11, double & tH1_only_H1_delta, double & tH1_without_H1_delta, double & tH1F, double & H1full, double & tgaga, double & tcga, double & tgamma2, double & LR, double & LF, double & LQ, double & A_F){
  static Logger logger("observable_set::calculate_sigma21");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double B2_A_F = 0.;
  static double sig21 = 0.;
  if (initial_gg){
    // B2_A_F recalculated - function 'calculate_B2_A_F()' could be used !!!
    B2_A_F = B2g + .5 * beta0 * A_F;
    // B2g must be B2(g)_A_F in the following formula (misleading name in counterterm_gg)
    sig21 = 
      - beta0 * LR * sig11 
      - B1g * H1full 
      - LF * tgaga  
      - B2_A_F * pdf_factor_x1x2 
      + beta0 * (tH1_without_H1_delta + tH1_only_H1_delta)
      - tcga 
      - tgamma2 
      - C1ggdelta(A_F) * tH1F 
      - 2 * Delta2gg * pdf_factor_x1x2
      + csi->order_alpha_s_born * beta0 * LR * tH1F;
    // replaces the following by using the alpha_S-order of the underlying Born process
    //      + 2 * beta0 * LR * tH1F;
    //20170522:
    // + 2 * beta0 * LR * (B1g * pdf_factor_x1x2  seems doubled
    // (contained in  -B1g * H1full)
    // org:      + 2 * beta0 * LR * (B1g * pdf_factor_x1x2 + tH1F);
    // new:      + 2 * beta0 * LR * tH1F;
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
      + beta0 * (tH1_without_H1_delta + tH1_only_H1_delta) 
      - tcga 
      - tgamma2 
      - C1qqdelta(A_F) * tH1F 
      - 2 * Delta2qq * pdf_factor_x1x2
      + csi->order_alpha_s_born * beta0 * LR * tH1F;
     // replaces 'QT_finalstate_massive_coloured' implementation by using the alpha_S-order of the underlying Born process
   
    // Q dependence // check if identical !!!
    sig21 += LQ * sig11 * beta0;
    sig21 -= LQ * H1full * A1q;
    sig21 += LQ * (tgaga + (B1q + 0.5 * A1q * LQ) * tH1F);
    sig21 -= LQ * A2q * pdf_factor_x1x2;
    //  sig21 += LQ*(sig11*(beta0-B1q-0.5*A1q*LQ) - A2q*pdf_factor_x1x2 + A1q*(2*beta0*LR*pdf_factor_x1x2-(tH1_without_H1_delta + tH1_only_H1_delta)+(LQ-LF)*tH1F) + tgaga);
  }

  if (QT_finalstate_massive_coloured){
    /*
    if (initial_qqx){
      // old comment: check if this part might be missing for qqx - solved ???
      sig21 += 2 * beta0 * LR * tH1F;  // alpha_S² at LO, not particularly because of mass - but it appears only in ttx production so far !!! ???
    }
    */
    // new implementation (should also work for simultaneous calculation of channels):
    if (initial_diag || initial_pdf_gq){
      sig21 += Gamma1born * (tH1_without_H1_delta + LF * tH1F);
      // old comment: only in the non-diagonal channel, where H1full does not contain a virtual amplitudes - solved after splitting off the H1full contribution ???
    }
    if (initial_diag){
      sig21 += 
	- Gamma1born * 2 * beta0 * LR * pdf_factor_x1x2
	+ Gamma1loop * pdf_factor_x1x2
	+ Gamma2 * pdf_factor_x1x2;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return sig21;
}


// new version, where QT_H1_delta part is treated separately (old version further exists in order to keep the VT implementation unchanged for this moment...) !!!
// no LR dependence !!! ???
double observable_set::calculate_H1(double pdf_factor_x1x2, double tH1_only_H1_delta, double tH1_without_H1_delta, double tH1F, double LR, double LF, double LQ){
  static Logger logger("observable_set::calculate_H1");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double result_H1 = 0.;
  //  if (name_process[0] == 'g'){
  if (initial_gg){
    result_H1 = 
      + tH1_only_H1_delta 
      // must be without H1_delta term !!! !!!
      + tH1_without_H1_delta 
      + LF * tH1F 
      - 2 * beta0 * LR * pdf_factor_x1x2;
    // LQ-dependent term missing !!!
    logger << LOG_DEBUG_VERBOSE << setw(40) << "result_H1" << " = " << result_H1 << endl;
    logger << LOG_DEBUG_VERBOSE << setw(40) << "tH1" << " = " << tH1_only_H1_delta + tH1_without_H1_delta << endl;
    logger << LOG_DEBUG_VERBOSE << setw(40) << "tH1_only_H1_delta" << " = " << tH1_only_H1_delta << endl;
    logger << LOG_DEBUG_VERBOSE << setw(40) << "tH1_without_H1_delta" << " = " << tH1_without_H1_delta << endl;
    logger << LOG_DEBUG_VERBOSE << setw(40) << "LF * tH1F" << " = " << LF * tH1F << endl;
    logger << LOG_DEBUG_VERBOSE << setw(40) << "- 2 * beta0 * LR * pdf_factor_x1x2" << " = " << - 2 * beta0 * LR * pdf_factor_x1x2 << endl;

  }
  //  else {
  else if (initial_qqx){
    result_H1 = 
      + tH1_only_H1_delta 
      + tH1_without_H1_delta 
      + (LF - LQ) * tH1F 
      - LQ * pdf_factor_x1x2 * (B1q + 0.5 * A1q * LQ);
    logger << LOG_DEBUG_VERBOSE << setw(40) << "result_H1" << " = " << result_H1 << endl;
    logger << LOG_DEBUG_VERBOSE << setw(40) << "tH1" << " = " << tH1_only_H1_delta + tH1_without_H1_delta << endl;
    logger << LOG_DEBUG_VERBOSE << setw(40) << "tH1_only_H1_delta" << " = " << tH1_only_H1_delta << endl;
    logger << LOG_DEBUG_VERBOSE << setw(40) << "tH1_without_H1_delta" << " = " << tH1_without_H1_delta << endl;
  }
  else {logger << LOG_DEBUG_VERBOSE << "No reasonable core process chosen." << endl;}

  if (QT_finalstate_massive_coloured){
    //    if (name_process[0] == 'g'){
    if (initial_gg){
    }
    else if (initial_qqx){
      //    else {
      result_H1 += - 2 * beta0 * LR * pdf_factor_x1x2;  // alpha_S² at LO, not particularly because of mass - but it appears only in ttx production so far !!!
    }
    else {logger << LOG_DEBUG_VERBOSE << "No reasonable core process chosen." << endl;}
  }
    logger << LOG_DEBUG_POINT << "separated H1 function: result_H1 = " << result_H1 << endl;
    logger << LOG_DEBUG_POINT << "separated H1 function: beta0 = " << beta0 << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return result_H1;
}

// no LR dependence !!!
double observable_set::calculate_H1(double pdf_factor_x1x2, double tH1, double tH1F, double LR, double LF, double LQ){
  static Logger logger("observable_set::calculate_H1 (combined: only one tH1)");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double result_H1 = 0.;
  //  if (name_process[0] == 'g'){
  if (initial_gg){
    result_H1 = 
      + tH1 
      + LF * tH1F 
      - csi->order_alpha_s_born * beta0 * LR * pdf_factor_x1x2;
    // replaces the following by using the alpha_S-order of the underlying Born process
    //      - 2 * beta0 * LR * pdf_factor_x1x2;
    // LQ-dependent term missing !!!
  }
  //  else {
  else if (initial_qqx){
    result_H1 = 
      + tH1 
      + (LF - LQ) * tH1F
      - csi->order_alpha_s_born * beta0 * LR * pdf_factor_x1x2
    // replaces 'QT_finalstate_massive_coloured' implementation by using the alpha_S-order of the underlying Born process
      - LQ * pdf_factor_x1x2 * (B1q + 0.5 * A1q * LQ);
  }
  
  else {logger << LOG_DEBUG_VERBOSE << "No reasonable core process chosen." << endl;}

  logger << LOG_DEBUG_POINT << setw(40) << "result_H1" << " = " << result_H1 << endl;
  logger << LOG_DEBUG_POINT << setw(40) << "tH1" << " = " << tH1 << endl;
  logger << LOG_DEBUG_POINT << setw(40) << "LF * tH1F" << " = " << LF * tH1F << endl;
  logger << LOG_DEBUG_POINT << setw(40) << "- 2 * beta0 * LR * pdf_factor_x1x2" << " = " << - 2 * beta0 * LR * pdf_factor_x1x2 << endl;
  logger << LOG_DEBUG_POINT << setw(40) << "csi->order_alpha_s_born = " << csi->order_alpha_s_born << endl;
  logger << LOG_DEBUG_POINT << setw(40) << "beta0" << " = " << beta0 << endl;
  logger << LOG_DEBUG_POINT << setw(40) << "LR" << " = " << LR << endl;
  logger << LOG_DEBUG_POINT << setw(40) << "pdf_factor_x1x2" << " = " << pdf_factor_x1x2 << endl;
  
  logger << LOG_DEBUG_POINT << "tH1 = " << tH1 << endl;
  logger << LOG_DEBUG_POINT << "tH1F = " << tH1F << endl;
  //   logger << LOG_DEBUG_POINT << " = " <<  << endl;

  /*
  if (QT_finalstate_massive_coloured){
    if (initial_gg){
      //    if (name_process[0] == 'g'){
    }
    //    else {
    else if (initial_qqx){
      result_H1 += - 2 * beta0 * LR * pdf_factor_x1x2;  // alpha_S² at LO, not particularly because of mass - but it appears only in ttx production so far !!!
    }
    else {logger << LOG_DEBUG_VERBOSE << "No reasonable core process chosen." << endl;}
  }
  */
  logger << LOG_DEBUG_POINT << "result_H1 = " << result_H1 << endl;
  logger << LOG_DEBUG_POINT << "beta0 = " << beta0 << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return result_H1;
}



double observable_set::calculate_H2(double pdf_factor_x1x2, double tH1, double tH1F, double H1full, double sig11, double tgaga, double tcga, double tgamma2, double tH2, double LR, double LF, double LQ){
  static Logger logger("observable_set::calculate_H2");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG_VERBOSE << "name_process[0] = " << name_process[0] << endl;

  static double result_H2 = 0.;
  //  if (name_process[0] == 'g'){
  if (initial_gg){
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
	 + LF * 2 * Delta2gg * pdf_factor_x1x2); // equivalent in qqx

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
  else if (initial_qqx){
    //  else {
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

	 // Manually added (terms missing because they are zero in colourless case due to csi->order_alpha_s_born = 0) !!!
	 // Insert analogous term to gg channel here:
	 - .5 * csi->order_alpha_s_born * beta0 * LF * LR * tH1F
	 // adds missing beta0 dependence of LF * C1qqdelta(A_F) * tH1F ???
	 - .5 * csi->order_alpha_s_born * beta0 * LR * H1full
	 - .5 * csi->order_alpha_s_born * beta0 * LR * tH1
	 - csi->order_alpha_s_born * (0.5 * (beta0 * LR) * (beta0 * LR) + beta1 * LR) * pdf_factor_x1x2
	 // until here.
	 
	 ///
	 /// below here: only relevant for resummation:
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
  else {logger << LOG_DEBUG_VERBOSE << "No reasonable core process chosen." << endl;}


  if (QT_finalstate_massive_coloured){
    logger << LOG_DEBUG_VERBOSE << "Soft function in H2 not implemented yet." << endl;
    //  logger << LOG_FATAL << "H2 not implemented yet." << endl; exit(1);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return result_H2;
}



// !!! resummation changes must be imported in complete file !!!










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

//-------------------------------------------------------------------------------------------------------------------------------  
//-------------------------------------------------------------------------------------------------------------------------------
//#include "OV.observable.qTsubtraction.cpp"

