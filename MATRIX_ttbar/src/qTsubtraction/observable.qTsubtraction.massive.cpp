#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

extern "C" {
  void fourcorrelators_(double* B4, int *i1, int *i2, int *i3, int *i4, int *channel);
  double callmyli3_(double* x, double* result);
}

void observable_set::determine_correlationoperator_QCD(){
 double m_HQ = msi.M[csi->type_parton[0][3]];
  double m2_HQ = m_HQ*m_HQ;
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
      if      (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 0;}
      else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] == 0.){temp_massive = 1;}
      else if (M[abs(csi->type_parton[0][temp_no_emitter])] == 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 2;}
      else if (M[abs(csi->type_parton[0][temp_no_emitter])] != 0. && M[abs(csi->type_parton[0][temp_no_spectator])] != 0.){temp_massive = 3;}
      else {
	logger << LOG_WARN << "Should not happen!" << endl;
	logger << LOG_WARN << "mass_parton[0][temp_no_emitter = " << temp_no_emitter << "] = " << M[abs(csi->type_parton[0][temp_no_emitter])] << endl;
	logger << LOG_WARN << "mass_parton[0][temp_no_spectator = " << temp_no_spectator << "] = " << M[abs(csi->type_parton[0][temp_no_spectator])] << endl;
      }
//
      int temp_type_combination = 0;
      if      ((temp_no_emitter > 2 && temp_no_spectator > 2) && temp_no_emitter == temp_no_spectator){temp_type_combination = 1;}
      else if (temp_no_emitter > 2 && temp_no_spectator > 2){temp_type_combination = 2;}
      else if (temp_no_emitter < 3 && temp_no_spectator > 2){temp_type_combination = 3;}
      else if ((temp_no_emitter < 3 && temp_no_spectator < 3) && temp_no_emitter == temp_no_spectator){temp_type_combination = 4;}
      else if ((temp_no_emitter < 3 && temp_no_spectator < 3)){temp_type_combination = 5;}
      else {logger << LOG_WARN << "Combination should not appear." << endl; exit(1);}

      double temp_charge_factor = 1.;
      if (temp_type_combination == 1 || temp_type_combination == 4){
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



// Introduce a function where all required kinematic quantities are calculated (stored in observable_set) !!!


void observable_set::calculate_QT_CM_kinematics(){
 double m_HQ = msi.M[csi->type_parton[0][3]];
  double m2_HQ = m_HQ*m_HQ;
  /*
  // VT

  // CT
  if (initial_pdf_diag){
    calculate_Gamma1born();
  }

  // VT2
  if (initial_pdf_diag || initial_pdf_gq){
    calculate_Gamma1born();
  }
  
  // CT2
  if (initial_diag){
    calculate_Gamma1loop();
    if (initial_diag_gg){calculate_commutator_Gamma1_Ft1();}
    calculate_Gamma1_squared();
    calculate_Gamma2();  // calls  calculate_Ft1born() , uses  Ft1born  and  commutator_Gammat1_Ft1
  }
  if (initial_pdf_diag || initial_pdf_gq){
    calculate_Gamma1born();
  }
  */
  
  q2 = (p_parton[0][3] + p_parton[0][4]).m2();
  // calculate_Gamma1born
  // calculate_Gamma1loop
  // calculate_commutator_Gamma1_Ft1
  // calculate_Gamma2
  // calculate_Gamma1_squared
  v34 = sqrt(1. - pow(2 * m2_HQ / (2 * (p_parton[0][3] * p_parton[0][4])), 2));
  // calculate_Ft1born
  // calculate_Gamma1born
  // calculate_Gamma1loop
  // calculate_commutator_Gamma1_Ft1
  // calculate_Gamma2
  // calculate_Gamma1_squared
  
  double frac_v34 = (1. + v34) / (1. - v34);
  // calculate_Ft1born
  // calculate_commutator_Gamma1_Ft1
  // 
  // 
  // 

  log_v34 = log(frac_v34);
  // calculate_Ft1born
  // calculate_Gamma1born
  // calculate_Gamma1loop
  // calculate_commutator_Gamma1_Ft1
  // calculate_Gamma2
  // calculate_Gamma1_squared

  dilogpT2m2 = gsl_sf_dilog(-p_parton[0][3].pT2() / m2_HQ);
  // calculate_Ft1born
  // calculate_commutator_Gamma1_Ft1
  // 
  // 
  // 
  logmT2m2 = log(p_parton[0][3].ET2() / m2_HQ);
  // calculate_Ft1born
  // calculate_commutator_Gamma1_Ft1
  // 
  // 
  // 
  double y34 = p_parton[0][3].rapidity() - p_parton[0][4].rapidity();
  L34 = log_v34 * logmT2m2
    - 2 * gsl_sf_dilog(2 * v34 / (1. + v34))
    - .25 * pow(log_v34, 2)
    + 2 * (gsl_sf_dilog(1. - exp(y34) / sqrt(frac_v34)) + gsl_sf_dilog(1. - exp(-y34) / sqrt(frac_v34)) + .5 * pow(y34, 2));
  // calculate_Ft1born
  // calculate_commutator_Gamma1_Ft1
  // 
  // 
  //

  //  gammacusp2_v
  double beta34 = .5 * log_v34;
  double arg_1mv34_1pv34 = (1. - v34) / (1. + v34);
  double dilog1mv34_1pv34 = gsl_sf_dilog(arg_1mv34_1pv34);
  double trilog1mv34_1pv34 = 0.;
  callmyli3_(&arg_1mv34_1pv34, &trilog1mv34_1pv34);
  double log_2v34_1pv34 = log(2 * v34 / (1. + v34));
  gammacusp2_v = C_A / 2. * (-5 * zeta2 + zeta3 + pow(beta34, 2)
			     + pow(1. / v34, 2) * (trilog1mv34_1pv34 + beta34 * dilog1mv34_1pv34 + pow(beta34, 3) / 3 - (5. / 6) * pi2 * beta34 - zeta3)
			     + (1. / v34) * (dilog1mv34_1pv34 - 2. * beta34 * log_2v34_1pv34 + (5. / 6) * pi2 * beta34 - pow(beta34, 2) + (5. / 6) * pi2 - pow(beta34, 3) / 3));
  // calculate_Gamma2

  

}





void observable_set::calculate_Ft1born(){
 double m_HQ = msi.M[csi->type_parton[0][3]];
  double m2_HQ = m_HQ*m_HQ; 
 static Logger logger("observable_set::calculate_Ft1born");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  double v34 = sqrt(1. - pow(2 * m2_HQ / (2 * (p_parton[0][3] * p_parton[0][4])), 2));
  double frac_v34 = (1. + v34) / (1. - v34);
  double log_v34 = log(frac_v34);
  double y34 = p_parton[0][3].rapidity() - p_parton[0][4].rapidity();
  double dilogpT2m2 = gsl_sf_dilog(-p_parton[0][3].pT2() / m2_HQ);
  double logmT2m2 = log(p_parton[0][3].ET2() / m2_HQ);
  double L34 = log_v34 * logmT2m2
    - 2 * gsl_sf_dilog(2 * v34 / (1. + v34))
    - .25 * pow(log_v34, 2)
    + 2 * (gsl_sf_dilog(1. - exp(y34) / sqrt(frac_v34)) + gsl_sf_dilog(1. - exp(-y34) / sqrt(frac_v34)) + .5 * pow(y34, 2));

  Ft1born = 0.;
  for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
    if (QT_correlationoperator[i_c].type_combination == 1){
      Ft1born += (logmT2m2 + dilogpT2m2) * QT_ME2_cf[i_c];
    }
    else if (QT_correlationoperator[i_c].type_combination == 2){
      Ft1born += (2 * dilogpT2m2 + (1. / v34) * L34) * QT_ME2_cf[i_c];
    }
  }

  logger << LOG_DEBUG_POINT << "Ft1born = " << Ft1born << endl;
  
  if (type_contribution == "CT2" && initial_diag){
    for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
      Ft1born_4correlator[i_c] = 0.;
      int i1 = QT_correlationoperator[i_c].no_emitter;
      int i2 = QT_correlationoperator[i_c].no_spectator;
      
      for (int j_c = 0; j_c < QT_correlationoperator.size(); j_c++){
	double temp_B4 = 0.;
	double temp_B4hc = 0.;
	int j1 = QT_correlationoperator[j_c].no_emitter;
	int j2 = QT_correlationoperator[j_c].no_spectator;
	if (QT_correlationoperator[j_c].type_combination == 1 || QT_correlationoperator[j_c].type_combination == 2){
	  fourcorrelators_(&temp_B4, &i1, &i2, &j1, &j2, &initial_channel);
	  fourcorrelators_(&temp_B4hc, &j1, &j2, &i1, &i2, &initial_channel);
	  logger << LOG_DEBUG_POINT << "fourcorrelators (" << i_c << " - " << j_c << ") = temp_B4   (" << setw(2) << i1 << ", " << setw(2) << i2 << ", " << setw(2) << j1 << ", " << setw(2) << j2 << ") = " << setw(23) << setprecision(15) << temp_B4 << endl; 
	  logger << LOG_DEBUG_POINT << "fourcorrelators (" << i_c << " - " << j_c << ") = temp_B4hc (" << setw(2) << j1 << ", " << setw(2) << j2 << ", " << setw(2) << i1 << ", " << setw(2) << i2 << ") = " << setw(23) << setprecision(15) << temp_B4hc << endl;
	  ///	  temp_B4 = 0.;
	  ///	  temp_B4hc = 0.;
	}
	if (QT_correlationoperator[j_c].type_combination == 1){
	  Ft1born_4correlator[i_c] += (logmT2m2 + dilogpT2m2) * .5 * (temp_B4 + temp_B4hc);
	  logger << LOG_DEBUG_POINT << "Ft1born_4correlator[" << i_c << "][" << j_c << "] = " << (logmT2m2 + dilogpT2m2) * .5 * (temp_B4 + temp_B4hc) << endl;
	}
	else if (QT_correlationoperator[j_c].type_combination == 2){
	  Ft1born_4correlator[i_c] += (2 * dilogpT2m2 + (1. / v34) * L34) * .5 * (temp_B4 + temp_B4hc);
	  logger << LOG_DEBUG_POINT << "Ft1born_4correlator[" << i_c << "][" << j_c << "] = " << (2 * dilogpT2m2 + (1. / v34) * L34) * .5 * (temp_B4 + temp_B4hc) << endl;
	}
      }
      logger << LOG_DEBUG_POINT << "Ft1born_4correlator[" << i_c << "] = " << Ft1born_4correlator[i_c] << endl;
      // Just to check the impact of the contribution !!!
      //      Ft1born_4correlator[i_c] = 0.;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::calculate_Gamma1born(){

  static Logger logger("observable_set::calculate_Gamma1born");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

//Substitution of all the msi.M2_t with m2_HQ
  double m_HQ = msi.M[csi->type_parton[0][3]];
  double m2_HQ = m_HQ*m_HQ;

  if (QT_finalstate_massive_coloured){
    double q2 = (p_parton[0][3] + p_parton[0][4]).m2();

    double v34 = sqrt(1. - pow(m2_HQ / (p_parton[0][3] * p_parton[0][4]), 2));
    double log_v34 = log((1. + v34) / (1. - v34));
    Gamma1born = 0.;
    for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
      if (QT_correlationoperator[i_c].type_combination == 1){
	Gamma1born += .5 * QT_ME2_cf[i_c];
      }
      else if (QT_correlationoperator[i_c].type_combination == 2){
	Gamma1born += .5 * (1. / v34) * log_v34 * QT_ME2_cf[i_c];
      }
      else if (QT_correlationoperator[i_c].type_combination == 3){
	Gamma1born += .5 * log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * m2_HQ)) * QT_ME2_cf[i_c];
      }
    }
    //Gammat = -2 * Gammat1(Eg.33, arXiv:1408.4564) !!! -> 2 * Re(Gammat1)
 }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::calculate_Gamma1loop(){
 double m_HQ = msi.M[csi->type_parton[0][3]];
  double m2_HQ = m_HQ*m_HQ;
  static Logger logger("observable_set::calculate_Gamma1loop");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  if (QT_finalstate_massive_coloured){
    Gamma1loop = 0.;
    double q2 = (p_parton[0][3] + p_parton[0][4]).m2();
    double v34 = sqrt(1. - pow(m2_HQ / (p_parton[0][3] * p_parton[0][4]), 2));
    double log_v34 = log((1. + v34) / (1. - v34));
    for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
      if (QT_correlationoperator[i_c].type_combination == 1){
	Gamma1loop += .5 * QT_ME2_loopcf[i_c];
	logger << LOG_DEBUG_POINT << "QT_ME2_loopcf[" << i_c << "] = " << setw(23) << setprecision(15) << QT_ME2_loopcf[i_c] << endl;
      }
      else if (QT_correlationoperator[i_c].type_combination == 2){
	Gamma1loop += .5 * (1. / v34) * log_v34 * QT_ME2_loopcf[i_c];
	logger << LOG_DEBUG_POINT << "QT_ME2_loopcf[" << i_c << "] = " << setw(23) << setprecision(15) << QT_ME2_loopcf[i_c] << endl;
      }
      else if (QT_correlationoperator[i_c].type_combination == 3){
	Gamma1loop += .5 * log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * m2_HQ)) * QT_ME2_loopcf[i_c];
	logger << LOG_DEBUG_POINT << "QT_ME2_loopcf[" << i_c << "] = " << setw(23) << setprecision(15) << QT_ME2_loopcf[i_c] << endl;
      }
      else {
	logger << LOG_DEBUG_POINT << "QT_ME2_loopcf[" << i_c << "] = " << setw(23) << setprecision(15) << QT_ME2_loopcf[i_c] << "   (not used)" << endl;
      }
    }
    // Gammat = -2 * Gammat1(Eg.33, arXiv:1408.4564) !!! -> 2 * Re(Gammat1)
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::calculate_commutator_Gamma1_Ft1(){
 double m_HQ = msi.M[csi->type_parton[0][3]];
  double m2_HQ = m_HQ*m_HQ;
  static Logger logger("observable_set::calculate_commutator_Gamma1_Ft1");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  int temp_channel = 0;
  if (initial_gg){temp_channel = 1;}
  else if (initial_qqx){temp_channel = 2;}

  double q2 = (p_parton[0][3] + p_parton[0][4]).m2();
  double v34 = sqrt(1. - pow(2 * m2_HQ / (2 * (p_parton[0][3] * p_parton[0][4])), 2));
  double frac_v34 = (1. + v34) / (1. - v34);
  double log_v34 = log(frac_v34);
  double y34 = p_parton[0][3].rapidity() - p_parton[0][4].rapidity();
  double dilogpT2m2 = gsl_sf_dilog(-p_parton[0][3].pT2() / m2_HQ);
  double logmT2m2 = log(p_parton[0][3].ET2() / m2_HQ);
  double L34 = log_v34 * logmT2m2
    - 2 * gsl_sf_dilog(2 * v34 / (1. + v34))
    - .25 * pow(log_v34, 2)
    + 2 * (gsl_sf_dilog(1. - exp(y34) / sqrt(frac_v34)) + gsl_sf_dilog(1. - exp(-y34) / sqrt(frac_v34)) + .5 * pow(y34, 2));

  commutator_Gammat1_Ft1 = 0.;
  for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
    int i1 = QT_correlationoperator[i_c].no_emitter;
    int i2 = QT_correlationoperator[i_c].no_spectator;
    for (int j_c = 0; j_c < QT_correlationoperator.size(); j_c++){
      int j1 = QT_correlationoperator[j_c].no_emitter;
      int j2 = QT_correlationoperator[j_c].no_spectator;
      double temp_B4 = 0.;
      double temp_B4ro = 0.;
      fourcorrelators_(&temp_B4, &i1, &i2, &j1, &j2, &temp_channel);
      logger << LOG_DEBUG_POINT << "fourcorrelators (" << i_c << " - " << j_c << ") = temp_B4   (" << setw(2) << i1 << ", " << setw(2) << i2 << ", " << setw(2) << j1 << ", " << setw(2) << j2 << ") = " << setw(23) << setprecision(15) << temp_B4 << endl; 
      fourcorrelators_(&temp_B4ro, &j1, &j2, &i1, &i2, &temp_channel);
      logger << LOG_DEBUG_POINT << "fourcorrelators (" << i_c << " - " << j_c << ") = temp_B4ro (" << setw(2) << j1 << ", " << setw(2) << j2 << ", " << setw(2) << i1 << ", " << setw(2) << i2 << ") = " << setw(23) << setprecision(15) << temp_B4ro << endl;
      //	  temp_B4 = 0.;
      //	  temp_B4ro = 0.;
      //	(temp_B4 - temp_B4ro)  is the 4-colour-correlator commutator !!!
      
	if (QT_correlationoperator[i_c].type_combination == 1 && QT_correlationoperator[j_c].type_combination == 1){
	  commutator_Gammat1_Ft1 += .5 * (logmT2m2 + dilogpT2m2) * (temp_B4 - temp_B4ro);
	  logger << LOG_DEBUG_VERBOSE 
		 << " type_contribution[" << i_c << "] = " << QT_correlationoperator[i_c].type_combination << "   "
		 << " type_contribution[" << j_c << "] = " << QT_correlationoperator[j_c].type_combination << "   "
		 << .5 * (logmT2m2 + dilogpT2m2) * (temp_B4 - temp_B4ro)
		 << endl;
	}
	else if (QT_correlationoperator[i_c].type_combination == 1 && QT_correlationoperator[j_c].type_combination == 2){
	  commutator_Gammat1_Ft1 += .5 * (2 * dilogpT2m2 + (1. / v34) * L34) * (temp_B4 - temp_B4ro);
	  logger << LOG_DEBUG_VERBOSE 
		 << " type_contribution[" << i_c << "] = " << QT_correlationoperator[i_c].type_combination << "   "
		 << " type_contribution[" << j_c << "] = " << QT_correlationoperator[j_c].type_combination << "   "
		 << .5 * (2 * dilogpT2m2 + (1. / v34) * L34) * (temp_B4 - temp_B4ro)
		 << endl;
	}
	//	else if (QT_correlationoperator[i_c].type_combination == 1 && QT_correlationoperator[j_c].type_combination == 3){
	// already covered in 3 - 1 ???
	//	}
	else if (QT_correlationoperator[i_c].type_combination == 2 && QT_correlationoperator[j_c].type_combination == 1){
	  commutator_Gammat1_Ft1 += .5 * (2 * dilogpT2m2 + (1. / v34) * L34) * (temp_B4 - temp_B4ro);
	}
	else if (QT_correlationoperator[i_c].type_combination == 2 && QT_correlationoperator[j_c].type_combination == 2){
	  commutator_Gammat1_Ft1 += .5 * (1. / v34) * log_v34 * (2 * dilogpT2m2 + (1. / v34) * L34) * (temp_B4 - temp_B4ro);
	  logger << LOG_DEBUG_VERBOSE 
		 << " type_contribution[" << i_c << "] = " << QT_correlationoperator[i_c].type_combination << "   "
		 << " type_contribution[" << j_c << "] = " << QT_correlationoperator[j_c].type_combination << "   "
		 << .5 * (1. / v34) * log_v34 * (2 * dilogpT2m2 + (1. / v34) * L34) * (temp_B4 - temp_B4ro)
		 << endl;
	}
	//	else if (QT_correlationoperator[i_c].type_combination == 2 && QT_correlationoperator[j_c].type_combination == 3){
	// already covered in 3 - 2 ???
	//	}
	else if (QT_correlationoperator[i_c].type_combination == 3 && QT_correlationoperator[j_c].type_combination == 1){
	  commutator_Gammat1_Ft1 += .5 * (logmT2m2 + dilogpT2m2)
	    * (log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * m2_HQ))) * (temp_B4 - temp_B4ro);
	  logger << LOG_DEBUG_VERBOSE 
		 << " type_contribution[" << i_c << "] = " << QT_correlationoperator[i_c].type_combination << "   "
		 << " type_contribution[" << j_c << "] = " << QT_correlationoperator[j_c].type_combination << "   "
		 << .5 * (logmT2m2 + dilogpT2m2)
	    * (log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * m2_HQ))) * (temp_B4 - temp_B4ro)
		 << endl;
	}
	else if (QT_correlationoperator[i_c].type_combination == 3 && QT_correlationoperator[j_c].type_combination == 2){
	  commutator_Gammat1_Ft1 += .5 * (2 * dilogpT2m2 + (1. / v34) * L34)
	    * (log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * m2_HQ))) * (temp_B4 - temp_B4ro);
	  logger << LOG_DEBUG_VERBOSE
		 << " type_contribution[" << i_c << "] = " << QT_correlationoperator[i_c].type_combination << "   "
		 << " type_contribution[" << j_c << "] = " << QT_correlationoperator[j_c].type_combination << "   "
		 << .5 * (2 * dilogpT2m2 + (1. / v34) * L34)
	    * (log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * m2_HQ))) * (temp_B4 - temp_B4ro)
		 << endl;
	}
	//	else if (QT_correlationoperator[i_c].type_combination == 3 && QT_correlationoperator[j_c].type_combination == 3){
	// should not exist !!!
	//	}
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::calculate_Gamma2(){
 double m_HQ = msi.M[csi->type_parton[0][3]];
  double m2_HQ = m_HQ*m_HQ; 
 static Logger logger("observable_set::calculate_Gamma2");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  if (QT_finalstate_massive_coloured){

      double v34 = sqrt(1. - pow(m2_HQ / (p_parton[0][3] * p_parton[0][4]), 2));
      double log_v34 = log((1. + v34) / (1. - v34));
      double beta34 = .5 * log_v34;
      double arg_1mv34_1pv34 = (1. - v34) / (1. + v34);
      double dilog1mv34_1pv34 = gsl_sf_dilog(arg_1mv34_1pv34);
      //     double trilog1mv34_1pv34 = myLI3_(&arg_1mv34_1pv34);
      double trilog1mv34_1pv34 = 0.;
      callmyli3_(&arg_1mv34_1pv34, &trilog1mv34_1pv34);
      double log_2v34_1pv34 = log(2 * v34 / (1. + v34));
      gammacusp2_v = C_A / 2. * (-5 * zeta2 + zeta3 + pow(beta34, 2)
				 + pow(1. / v34, 2) * (trilog1mv34_1pv34 + beta34 * dilog1mv34_1pv34 + pow(beta34, 3) / 3 - (5. / 6) * pi2 * beta34 - zeta3)
				 + (1. / v34) * (dilog1mv34_1pv34 - 2. * beta34 * log_2v34_1pv34 + (5. / 6) * pi2 * beta34 - pow(beta34, 2) + (5. / 6) * pi2 - pow(beta34, 3) / 3));

    if (initial_pdf_diag){
      //    if (initial_pdf_gg){ // check which one is needed here !!!
      //    if ((name_process[0] == 'g' && name_process[1] == 'g') || (name_process[0] != 'g' && name_process[1] != 'g')){
      ///      logger << LOG_FATAL << "VT or CT2 term for gg and qqx channel is not validated (and not completely implemented) yet !!!" << endl;
      /*
      Gamma2q=0.5d0*gammacusp2q*(
------------------     .      1d0/v34*dlog((1d0+v34)/(1d0-v34))
------------------     .      *Tqq(3,4)
------------------     .      +dlog(4d0*dot(ptrans,1,3)**2/(q2*mt**2))*Tqq(1,3)
------------------     .      +dlog(4d0*dot(ptrans,1,4)**2/(q2*mt**2))*Tqq(1,4)
------------------     .      +dlog(4d0*dot(ptrans,2,3)**2/(q2*mt**2))*Tqq(2,3)
------------------     .      +dlog(4d0*dot(ptrans,2,4)**2/(q2*mt**2))*Tqq(2,4))
------------------     .         +gammacuspprime2q*Tqq(3,4)
------------------     .         -2d0*gammaQ1
------------------     .         +0.5d0*beta0*F1q
      */


      Gamma2 = 0.;

      double q2 = (p_parton[0][3] + p_parton[0][4]).m2();
      //      double v34 = sqrt(1. - pow(m2_HQ / (p_parton[0][3] * p_parton[0][4]), 2));
      //      double log_v34 = log((1. + v34) / (1. - v34));
      for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
	//	if (QT_correlationoperator[i_c].type_combination == 1){
	//	  Gamma2 += .5 * QT_ME2_cf[i_c];
	//	}
	if (QT_correlationoperator[i_c].type_combination == 2){
	  Gamma2 += .5 * gammacusp2 * (1. / v34) * log_v34 * QT_ME2_cf[i_c];
	  Gamma2 += .5 * 2 * gammacusp2_v * QT_ME2_cf[i_c];
	}
	else if (QT_correlationoperator[i_c].type_combination == 3){
	  Gamma2 += .5 * gammacusp2 * log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * m2_HQ)) * QT_ME2_cf[i_c];
	}
	// 	else if (QT_correlationoperator[i_c].type_combination == 0){
	//	}
	// three-correlator don't contribute for ttbar
      }
      Gamma2 += .5 * (-4 * gamma2Q); // multiplied by 1 (born / born)

      calculate_Ft1born();
      // Ft1 is also needed elsewhere and should thus better be calculated only once !!! Maybe already calculated beforehand whenever relevant (in qTsubtraction.ME2.cpp) ???
      Gamma2 += .5 * beta0 * Ft1born;
      //      Gamma2 += .5 * pi * beta0 * Ft1born; // pi if beta0 is defined with a pi in the denominator !!!

      // [Gamma_t^1, F_t^1] -> only in gg channel !!!
      // Maybe put  if (initial_diag_gg){ }  ???
      Gamma2 += .5 * commutator_Gammat1_Ft1;

      //      exit(1);
    }
    else if (initial_pdf_qqx){
      // Should never happen because of  if (initial_pdf_diag){}  ...  else if (initial_pdf_qqx){}  ???
      logger << LOG_FATAL << "VT or CT2 term for qqx channel is not implemented yet !!!" << endl;
      exit(1);
    }
    else if (initial_pdf_gq){
      logger << LOG_FATAL << "Gamma²_t does not contribute in the gq channel !!!" << endl;
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_Gamma1_squared(){
 double m_HQ = msi.M[csi->type_parton[0][3]];
  double m2_HQ = m_HQ*m_HQ;
  static Logger logger("observable_set::calculate_Gamma1_squared");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  if (QT_finalstate_massive_coloured){
    // 4 * Re (-1/4) (-1/4) * {...}
    
    int temp_channel = 0;
    if (initial_gg){temp_channel = 1;}
    else if (initial_qqx){temp_channel = 2;}

    double q2 = (p_parton[0][3] + p_parton[0][4]).m2();
    double v34 = sqrt(1. - pow(m2_HQ / (p_parton[0][3] * p_parton[0][4]), 2));
    double log_v34 = log((1. + v34) / (1. - v34));

    Gamma1squared = 0.;
    for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
      int i1 = QT_correlationoperator[i_c].no_emitter;
      int i2 = QT_correlationoperator[i_c].no_spectator;
      for (int j_c = 0; j_c < QT_correlationoperator.size(); j_c++){
	int j1 = QT_correlationoperator[j_c].no_emitter;
	int j2 = QT_correlationoperator[j_c].no_spectator;
	double temp_B4 = 0.;
	fourcorrelators_(&temp_B4, &i1, &i2, &j1, &j2, &temp_channel);
	  logger << LOG_DEBUG_POINT << "fourcorrelators (" << i_c << " - " << j_c << ") = temp_B4   (" << setw(2) << j1 << ", " << setw(2) << j2 << ", " << setw(2) << i1 << ", " << setw(2) << i2 << ") = " << setw(23) << setprecision(15) << temp_B4 << endl;
	  ///	  temp_B4 = 0.;
	
	if (QT_correlationoperator[i_c].type_combination == 1 && QT_correlationoperator[j_c].type_combination == 1){
	  Gamma1squared += .25 * temp_B4; // e.g. QT_ME2_cfcf[i_c][j_c];
	}
	else if (QT_correlationoperator[i_c].type_combination == 1 && QT_correlationoperator[j_c].type_combination == 2){
	  Gamma1squared += .25 * temp_B4 
	    * (1. / v34) * log_v34;
	}
	else if (QT_correlationoperator[i_c].type_combination == 1 && QT_correlationoperator[j_c].type_combination == 3){
	  Gamma1squared += .25 * temp_B4
	    * log(pow(2 * (p_parton[0][QT_correlationoperator[j_c].no_emitter] * p_parton[0][QT_correlationoperator[j_c].no_spectator]), 2) / (q2 * m2_HQ));
	}
	else if (QT_correlationoperator[i_c].type_combination == 2 && QT_correlationoperator[j_c].type_combination == 1){
	  Gamma1squared += .25 * temp_B4
	    * (1. / v34) * log_v34;
	}
	else if (QT_correlationoperator[i_c].type_combination == 2 && QT_correlationoperator[j_c].type_combination == 2){
	  Gamma1squared += .25 * temp_B4
	    * pow((1. / v34) * log_v34, 2);
	}
	else if (QT_correlationoperator[i_c].type_combination == 2 && QT_correlationoperator[j_c].type_combination == 3){
	  Gamma1squared += .25 * temp_B4
	    * (1. / v34) * log_v34 
	    * log(pow(2 * (p_parton[0][QT_correlationoperator[j_c].no_emitter] * p_parton[0][QT_correlationoperator[j_c].no_spectator]), 2) / (q2 * m2_HQ));
	}
	else if (QT_correlationoperator[i_c].type_combination == 3 && QT_correlationoperator[j_c].type_combination == 1){
	  Gamma1squared += .25 * temp_B4
	    * log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * m2_HQ));
	}
	else if (QT_correlationoperator[i_c].type_combination == 3 && QT_correlationoperator[j_c].type_combination == 2){
	  Gamma1squared += .25 * temp_B4
	    * (1. / v34) * log_v34 
	    * log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * m2_HQ));
	}
	else if (QT_correlationoperator[i_c].type_combination == 3 && QT_correlationoperator[j_c].type_combination == 3){
	  Gamma1squared += .25 * temp_B4
	    * log(pow(2 * (p_parton[0][QT_correlationoperator[i_c].no_emitter] * p_parton[0][QT_correlationoperator[i_c].no_spectator]), 2) / (q2 * m2_HQ)) 
	    * log(pow(2 * (p_parton[0][QT_correlationoperator[j_c].no_emitter] * p_parton[0][QT_correlationoperator[j_c].no_spectator]), 2) / (q2 * m2_HQ));
	}
      }
    }

    // ???
    if (initial_pdf_diag){
      //    if (initial_pdf_gg || initial_pdf_qqx){
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}







