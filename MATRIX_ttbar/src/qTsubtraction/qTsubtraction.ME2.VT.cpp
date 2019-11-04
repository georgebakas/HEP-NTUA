#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"
#include "../include/definitions.phasespace.set.cxx"

extern "C" {
  void xtest_();
  void deltaext_(double* P, double* mu_Q, double* m_HQ, double* H1, double* H1_T34, double* H1_T13, double* H1_T23, int *channel);
  //  void callmsqGGav_(double* P, double* msqGGav);
  void callmsqggav_(double* P, double* m_HQ, double* msqGGav);
  void callmsqdgav_(double* P, double* m_HQ, double* msqDGav);
  void fourcorrelators_(double* B4, int *i1, int *i2, int *i3, int *i4, int *channel);
}

void calculate_ME2_VT_QCD(observable_set & oset, call_generic & generic){
  static Logger logger("calculate_ME2_VT_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (oset.switch_OL){

  static double one = 1;
  static int n_momentum = 5 * (osi_n_particle + 2);
  double *P;
  P = new double[n_momentum];
  for (int i = 1; i < osi_p_parton[0].size(); i++){
    P[5 * (i - 1)]     = osi_p_parton[0][i].x0();
    P[5 * (i - 1) + 1] = osi_p_parton[0][i].x1();
    P[5 * (i - 1) + 2] = osi_p_parton[0][i].x2();
    P[5 * (i - 1) + 3] = osi_p_parton[0][i].x3();
    P[5 * (i - 1) + 4] = osi_p_parton[0][i].m();
  }

  if ((osi_type_contribution == "VT" ||
       osi_type_contribution == "CT2" ||
       osi_type_contribution == "VJ" ||
       osi_type_contribution == "CJ2") && osi_user_string_value[osi_user_string_map["model"]] != "Bornloop"){

    static double acc;
    static double M2L0;
    double *M2L1;
    double *IRL1;
    M2L1 = new double[3];
    IRL1 = new double[3];
    double *M2L2;
    double *IRL2;
    M2L2 = new double[5];
    IRL2 = new double[5];
    
    static char * OL_mu_ren = stch("muren");
    static char * OL_mu_reg = stch("mureg");
    static char * pole_uv = stch("pole_uv");
    static char * pole_ir1 = stch("pole_ir1");
    static char * pole_ir2 = stch("pole_ir2");
    ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
    ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
    ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
    
    double mu_Q = (osi_p_parton[0][1] + osi_p_parton[0][2]).m();
    ol_setparameter_double(OL_mu_ren, mu_Q);
    ol_setparameter_double(OL_mu_reg, mu_Q);
    static char * fact_uv = stch("fact_uv");
    static char * fact_ir = stch("fact_ir");
    ol_setparameter_double(fact_uv, one);
    ol_setparameter_double(fact_ir, one);
    
    
    logger << LOG_DEBUG_VERBOSE << "ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);" << endl;
    ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);
    osi_VA_V_ME2 = M2L1[0];
    logger << LOG_DEBUG_VERBOSE << "ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2);" << endl;
    ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2);
    logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_B = " << setw(23) << setprecision(15) << osi_VA_b_ME2 << endl; 
    logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_VX = " << setw(23) << setprecision(15) << osi_VA_V_ME2 + osi_VA_X_ME2 << endl; 

    if (osi_QT_finalstate_massive_coloured && (osi_type_contribution == "VT" || osi_type_contribution == "CT2")){

      double temp_H1 = 0.;
      double temp_H1_T34 = 0.;
      double temp_H1_T13 = 0.;
      double temp_H1_T23 = 0.;
      
      if (osi_type_contribution == "CT2" && oset.initial_diag){
	//    if (osi_type_contribution == "CT2" && ((oset.initial_gg && oset.initial_pdf_gg) || (oset.initial_qqx && oset.initial_pdf_qqx))){
	// Linked to Hayk's code: not yet connected to the rest of the code !!!

	logger << LOG_DEBUG_VERBOSE << "new psp:" << endl;
	int temp_channel = 0;
	if (oset.initial_gg){temp_channel = 1;}
	else if (oset.initial_qqx){temp_channel = 2;}

        // Replacing to generic heavy quark mass
	//double m_t = osi_msi.M_t;
        double m_HQ = osi_msi.M[oset.csi->type_parton[0][3]];

	// Here, mu_ren = mu_Q is set.
	deltaext_(P, &mu_Q, &m_HQ, &temp_H1, &temp_H1_T34, &temp_H1_T13, &temp_H1_T23, &temp_channel);
	//    xtest_();
	//    logger << LOG_DEBUG_VERBOSE << "temp_H1 = " << temp_H1 << endl;
	
	// temp_H1/H1_T34/H1_T13/H1_T23 are not normalized to Born here !!!
	temp_H1 = temp_H1 * pow(osi_alpha_S, 2);
	temp_H1_T34 = temp_H1_T34 * pow(osi_alpha_S, 2);
	temp_H1_T13 = temp_H1_T13 * pow(osi_alpha_S, 2);
	temp_H1_T23 = temp_H1_T23 * pow(osi_alpha_S, 2);
	// alpha_S running is included here (subtract via cc Born), which is again done later in observable.qTsubtraction.cpp routines.
	logger << LOG_DEBUG_POINT << "Hayk:   osi_QT_H1_delta_Hayk     = " << temp_H1 << endl;
	logger << LOG_DEBUG_POINT << "Hayk:   osi_QT_H1_T34_delta_Hayk = " << temp_H1_T34 << endl;
	logger << LOG_DEBUG_POINT << "Hayk:   osi_QT_H1_T13_delta_Hayk = " << temp_H1_T13 << endl;
	logger << LOG_DEBUG_POINT << "Hayk:   osi_QT_H1_T23_delta_Hayk = " << temp_H1_T23 << endl;
	logger << LOG_DEBUG_POINT << "Hayk:   T13_T23_T34_T33_Hayk     = " << temp_H1_T34 + temp_H1_T13 + temp_H1_T23 + C_F * temp_H1 << endl;
	
	double temp_B4 = 0.;
	int i1 = 1;
	int i2 = 3;
	int i3 = 2;
	int i4 = 3;
	
	fourcorrelators_(&temp_B4, &i1, &i2, &i3, &i4, &temp_channel);
	//    temp_H1 = temp_H1 * pow(osi_alpha_S, 2);
	logger << LOG_DEBUG_POINT << "Hayk:   B4(" << i1 << i2 << i3 << i4 << ") = " << temp_B4 << endl;
	// temp_B4 is normalized to the squared Born ampitude.
      }
      
      //#include "colourstrippedexperiment.qTsubtraction.ME2.VT.cpp"
      
      // different distinction needed: the partonic process is always gg or qqx !!!
      // e.g. check pdf contributions and determine from this if the channels 0->gg, 1->qqx, 2->gq, 3->rest are active !!!
      //    if ((oset.name_process[0] == 'g' && oset.name_process[1] != 'g') || (oset.name_process[0] == 'g' && oset.name_process[1] != 'g')){



      // The colour correlators  M2cc -> osi_QT_ME2_cf[i_c]  are needed in all cases (massive and massless --- not in the massless case, I guess ...)
      // The colour correlators  M2cc -> osi_QT_ME2_cf[i_c]  are needed in  VT  and CT2  contributions:
      static int n_cc = (osi_n_particle + 2) * (osi_n_particle + 1) / 2;
      static double ewcc;
      static double *M2cc;
      M2cc = new double[n_cc];
      ol_evaluate_cc(1, P, &osi_ME2, M2cc, &ewcc);
      
      for (int i_c = 0; i_c < n_cc; i_c++){
	logger << LOG_DEBUG_VERBOSE << "M2cc[" << i_c << "] = " << setw(23) << setprecision(15) << M2cc[i_c] << endl;
      }
      
      for (int i_c = 0; i_c < osi_QT_correlationoperator.size(); i_c++){
	//      if (osi_QT_correlationoperator[i_c].type_combination == 1 || osi_QT_correlationoperator[i_c].type_combination == 4){osi_QT_ME2_cf[i_c] = osi_QT_correlationoperator[i_c].charge_factor;}
	if (osi_QT_correlationoperator[i_c].type_combination == 1 || 
	    osi_QT_correlationoperator[i_c].type_combination == 4){
	  osi_QT_ME2_cf[i_c] = osi_QT_correlationoperator[i_c].charge_factor;
	}
	//      else if (osi_QT_correlationoperator[i_c].type_combination == 4){osi_QT_ME2_cf[i_c] = osi_QT_correlationoperator[i_c].charge_factor;}
	else {
	  osi_QT_ME2_cf[i_c] = M2cc[osi_QT_correlationoperator[i_c].no_BLHA_entry] / osi_ME2;
	}
	
	logger << LOG_DEBUG_POINT << "OpenLoops:  QT_ME2_cf[" << i_c << "] (em = " << osi_QT_correlationoperator[i_c].no_emitter << ", sp = " << osi_QT_correlationoperator[i_c].no_spectator << ") = " << setw(23) << setprecision(15) << osi_QT_ME2_cf[i_c] << "   charge_factor = " << setw(23) << setprecision(15) << osi_QT_correlationoperator[i_c].charge_factor << "   no_BLHA_entry = " << osi_QT_correlationoperator[i_c].no_BLHA_entry << endl;
      }
      //    delete [] M2cc;


      if (oset.QT_finalstate_massive_coloured && osi_type_contribution == "CT2" && oset.initial_diag){
	//    if (oset.QT_finalstate_massive_coloured && osi_type_contribution == "CT2" && ((oset.initial_gg && oset.initial_pdf_gg) || (oset.initial_qqx && oset.initial_pdf_qqx))){
	// The colour correlators  M2loopcc -> osi_QT_ME2_loopcf[i_c]  are needed in only in the diagonal qqx and gg channelss (massive):
	static double *M2loopcc;
	M2loopcc = new double[n_cc];

	//	double *M2L1;
	//	M2L1 = new double[3];

	double ME2loop = 0.;
	double ME2tree = 0.;

	static char * OL_ct_on = stch("ct_on");
	ol_setparameter_int(OL_ct_on, 1);
	ol_evaluate_loopcc(1, P, &ME2tree, M2L1, M2loopcc, &ewcc);
	//&ME2loop
	ME2loop = M2L1[0];
	//	ol_evaluate_loopcc(2, P, &ME2loop, M2loopcc, &ewcc);
	ol_setparameter_int(OL_ct_on, 0);
	
	logger << LOG_DEBUG_POINT << "OpenLoops:  ME2tree = " << setw(23) << setprecision(15) << osi_ME2 << "   ME2loop     = " << setw(23) << setprecision(15) << ME2tree << endl;
	logger << LOG_DEBUG_POINT << "OpenLoops:  ME2     = " << setw(23) << setprecision(15) << osi_ME2 << "   ME2loop     = " << setw(23) << setprecision(15) << ME2loop << endl;
	for (int i_c = 0; i_c < n_cc; i_c++){
	  logger << LOG_DEBUG_POINT << "OpenLoops:  M2cc[" << i_c << "] = " << setw(23) << setprecision(15) << M2cc[i_c] << "   M2loopcc[" << i_c << "] = " << setw(23) << setprecision(15) << M2loopcc[i_c] << endl;
	}
	// Check does not seem to be needed any longer...
	/*	
	logger << LOG_DEBUG_POINT << "check if ol_evaluate_full still gives the correct result:" << endl;
	logger << LOG_DEBUG_POINT << "(osi_VA_V_ME2 + osi_VA_X_ME2) * C_F / (alpha_S / pi) = " << C_F * (osi_VA_V_ME2 + osi_VA_X_ME2) / (osi_alpha_S / pi) << endl;
	ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);
	osi_VA_V_ME2 = M2L1[0];
	ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2);
	logger << LOG_DEBUG_POINT << "(osi_VA_V_ME2 + osi_VA_X_ME2) * C_F / (alpha_S / pi) = " << C_F * (osi_VA_V_ME2 + osi_VA_X_ME2) / (osi_alpha_S / pi) << endl;
	
	logger << LOG_DEBUG_POINT << "check if ol_evaluate_cc still gives the correct result:" << endl;
	delete [] M2cc;
	M2cc = new double[n_cc];
	ol_evaluate_cc(1, P, &osi_ME2, M2cc, &ewcc);
*/
	/*
	for (int i_c = 0; i_c < n_cc; i_c++){
	  logger << LOG_DEBUG_POINT << "M2cc[" << i_c << "] = " << setw(23) << setprecision(15) << M2cc[i_c] << endl;
	}
	
	// end checks.

	logger.newLine(LOG_DEBUG_POINT);
	logger << LOG_DEBUG_POINT << "Born: osi_ME2 * C_F = " << C_F * osi_ME2 << endl;
	logger << LOG_DEBUG_POINT << "ME2loop * C_F / (alpha_S / pi) = " << C_F * ME2loop / (osi_alpha_S / pi) << endl;
	ME2loop = (osi_VA_V_ME2 + osi_VA_X_ME2);
	logger << LOG_DEBUG_POINT << "(osi_VA_V_ME2 + osi_VA_X_ME2) * C_F / (alpha_S / pi) = " << C_F * (osi_VA_V_ME2 + osi_VA_X_ME2) / (osi_alpha_S / pi) << endl;
	*/
	logger << LOG_DEBUG_POINT << "Do not overwrite  osi_QT_ME2_loopcf  results with Hayk's ones!" << endl;
	
	for (int i_c = 0; i_c < osi_QT_correlationoperator.size(); i_c++){
	  if (osi_QT_correlationoperator[i_c].type_combination == 1 || 
	      osi_QT_correlationoperator[i_c].type_combination == 4){
	    osi_QT_ME2_loopcf[i_c] = osi_QT_correlationoperator[i_c].charge_factor * ME2loop / (osi_alpha_S / pi) / osi_ME2;
	    logger << LOG_DEBUG_POINT << "QT_ME2_loopcf[" << i_c << "] (em = " << osi_QT_correlationoperator[i_c].no_emitter << ", sp = " << osi_QT_correlationoperator[i_c].no_spectator << ") = " << setw(23) << setprecision(15) << osi_QT_ME2_loopcf[i_c] * osi_VA_b_ME2 << "   no cc from OpenLoops" << endl;
	  }
	  else {
	    osi_QT_ME2_loopcf[i_c] = M2loopcc[osi_QT_correlationoperator[i_c].no_BLHA_entry] / (osi_alpha_S / pi) / osi_ME2;

	    logger << LOG_DEBUG_POINT << "QT_ME2_loopcf[" << i_c << "] (em = " << osi_QT_correlationoperator[i_c].no_emitter << ", sp = " << osi_QT_correlationoperator[i_c].no_spectator << ") = " << setw(23) << setprecision(15) << osi_QT_ME2_loopcf[i_c] * osi_VA_b_ME2 << "   from OpenLoops" << endl;

	    /*
	    double temp = osi_QT_ME2_loopcf[i_c];
	    // overwrite amplitudes with Hayk's ones !!!
	    if (i_c == 1){osi_QT_ME2_loopcf[i_c] = (temp_H1_T34 + (osi_QT_correlationoperator[7].charge_factor - osi_QT_correlationoperator[0].charge_factor) * ME2loop / (osi_alpha_S / pi)) / osi_ME2;}
	    if (i_c == 2){osi_QT_ME2_loopcf[i_c] = temp_H1_T13 / osi_ME2;}
	    if (i_c == 3){osi_QT_ME2_loopcf[i_c] = -(temp_H1_T13 + temp_H1_T34 + osi_QT_correlationoperator[7].charge_factor * ME2loop / (osi_alpha_S / pi)) / osi_ME2;}
	    //	  if (i_c == 3){osi_QT_ME2_loopcf[i_c] = temp_H1_T23 / osi_ME2;}
	    if (i_c == 5){osi_QT_ME2_loopcf[i_c] = -(temp_H1_T13 + temp_H1_T34 + osi_QT_correlationoperator[7].charge_factor * ME2loop / (osi_alpha_S / pi)) / osi_ME2;}
	    //	  if (i_c == 5){osi_QT_ME2_loopcf[i_c] = temp_H1_T23 / osi_ME2;}
	    if (i_c == 6){osi_QT_ME2_loopcf[i_c] = temp_H1_T13 / osi_ME2;}
	    if (i_c == 8){osi_QT_ME2_loopcf[i_c] = temp_H1_T34 / osi_ME2;}
	    // Those results seem to be correct !!!
	    
	    logger << LOG_DEBUG_POINT << "QT_ME2_loopcf[" << i_c << "] (em = " << osi_QT_correlationoperator[i_c].no_emitter << ", sp = " << osi_QT_correlationoperator[i_c].no_spectator << ") = " << setw(23) << setprecision(15) << osi_QT_ME2_loopcf[i_c] * osi_VA_b_ME2 << "   from Hayk" << endl;
	    logger << LOG_DEBUG_POINT << "QT_ME2_loopcf[" << i_c << "] (em = " << osi_QT_correlationoperator[i_c].no_emitter << ", sp = " << osi_QT_correlationoperator[i_c].no_spectator << ") = " << setw(23) << setprecision(15) << osi_QT_ME2_loopcf[i_c] / temp << "   ratio" << endl;
	    */
	  }
	  // All  osi_QT_ME2_loopcf[i_c]  are normalized to Born !!!

	  //DEBUG_VERBOSE !!!
	  logger << LOG_DEBUG_VERBOSE << "QT_ME2_loopcf[" << i_c << "] (em = " << osi_QT_correlationoperator[i_c].no_emitter << ", sp = " << osi_QT_correlationoperator[i_c].no_spectator << ") = " << setw(23) << setprecision(15) 
		 << osi_QT_ME2_loopcf[i_c] * osi_VA_b_ME2 // / (osi_alpha_S / pi) //
		 << "   charge_factor = " << setw(23) << setprecision(15) << osi_QT_correlationoperator[i_c].charge_factor 
		 << "   no_BLHA_entry = " << osi_QT_correlationoperator[i_c].no_BLHA_entry << endl;
	}
	delete [] M2loopcc;
	
	//      logger << LOG_DEBUG_VERBOSE << "shifted osi_QT_H1_T34_delta_Hayk = " << temp_H1_T34 + 2 * oset.beta0 * log((osi_p_parton[0][1] + osi_p_parton[0][2]).m2() / pow(osi_var_mu_ren, 2)) * osi_QT_ME2_cf[8] * osi_VA_b_ME2 << endl;
	//      logger << LOG_DEBUG_VERBOSE << "shifted osi_QT_H1_T13_delta_Hayk = " << temp_H1_T13 + 2 * oset.beta0 * log((osi_p_parton[0][1] + osi_p_parton[0][2]).m2() / pow(osi_var_mu_ren, 2)) * osi_QT_ME2_cf[2] * osi_VA_b_ME2 << endl;
	//      logger << LOG_DEBUG_VERBOSE << "shifted osi_QT_H1_T23_delta_Hayk = " << temp_H1_T23 + 2 * oset.beta0 * log((osi_p_parton[0][1] + osi_p_parton[0][2]).m2() / pow(osi_var_mu_ren, 2)) * osi_QT_ME2_cf[5] * osi_VA_b_ME2 << endl;
	
	logger << LOG_DEBUG_POINT << "OL_checksum (1) = " << osi_QT_ME2_loopcf[0] + osi_QT_ME2_loopcf[1] + osi_QT_ME2_loopcf[2] + osi_QT_ME2_loopcf[3] << endl;
	logger << LOG_DEBUG_POINT << "OL_checksum (2) = " << osi_QT_ME2_loopcf[1] + osi_QT_ME2_loopcf[4] + osi_QT_ME2_loopcf[5] + osi_QT_ME2_loopcf[6] << endl;
	logger << LOG_DEBUG_POINT << "OL_checksum (3) = " << osi_QT_ME2_loopcf[2] + osi_QT_ME2_loopcf[5] + osi_QT_ME2_loopcf[7] + osi_QT_ME2_loopcf[8] << endl;
	logger << LOG_DEBUG_POINT << "OL_checksum (4) = " << osi_QT_ME2_loopcf[3] + osi_QT_ME2_loopcf[6] + osi_QT_ME2_loopcf[8] + osi_QT_ME2_loopcf[9] << endl;
	
      }
      delete [] M2cc;
    }
    delete [] M2L2;
    delete [] IRL2;
    delete [] M2L1;
    delete [] IRL1;

    // Check if this if statement changes the "rest" contribution !!!
    if (oset.initial_diag || oset.initial_pdf_gq){
    
    //  osi_VA_V_ME2 = 2 Re < M0 | M1 >  etc.
    osi_QT_H1_delta = (osi_VA_V_ME2 + osi_VA_X_ME2) / osi_VA_b_ME2 / (osi_alpha_S / pi);
    // same translation for osi_QT_H1_delta as in ppllll24_calculate_H2.
    //  (alpha_S / pi) normalization of H1
    if (oset.switch_polenorm == 1){
      if (oset.initial_gg){osi_QT_H1_delta -= pi2_6 * C_A;}
      else if (oset.initial_qqx){osi_QT_H1_delta -= pi2_6 * C_F;}
      else {logger << LOG_FATAL << "Wrong process specified." << endl; exit(1);}
    }
    
    }
    else {
      osi_QT_H1_delta = 0.;
    }

    /*
    // begin temporary: check how the mu_ren dependence is treated:
    logger << LOG_DEBUG << "V + X (mu = Q)  = " << osi_VA_V_ME2 + osi_VA_X_ME2 << endl;
    if (oset.psi->xbs_all.size() != 0){
      osi_VA_X_ME2 = osi_VA_X_ME2 - oset.csi->order_alpha_s_born * oset.beta0 * log(oset.psi->xbs_all[0][0] / pow(osi_var_mu_ren, 2)) * osi_VA_b_ME2 * (osi_alpha_S / pi);
      logger << LOG_DEBUG << "V + X (mu = Q)x = " << osi_VA_V_ME2 + osi_VA_X_ME2 << endl;
    }
    M2L1 = new double[3];
    IRL1 = new double[3];
    M2L2 = new double[5];
    IRL2 = new double[5];
    //  ol_setparameter_double(OL_mu_ren, mu_Q);
    ol_setparameter_double(OL_mu_reg, mu_Q);
    ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
    //    ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
    ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);
    osi_VA_V_ME2 = M2L1[0];
    ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2);
    logger << LOG_DEBUG << "V + X (mu = MZ) = " << osi_VA_V_ME2 + osi_VA_X_ME2 << endl;
    // end temporary
    */


  


  
  if (osi_QT_finalstate_massive_coloured){
    oset.calculate_Ft1born();
    logger << LOG_DEBUG_VERBOSE << "osi_QT_H1_delta_OL   = " << osi_QT_H1_delta * osi_VA_b_ME2 << endl;
    //  H1 -> H1 - It1  (It1 = 2 Re < M0 | It1 | M0 >) 
    logger << LOG_DEBUG_POINT << "osi_QT_H1_delta                     = " << osi_QT_H1_delta << endl;
    osi_QT_H1_delta += .5 * oset.Ft1born;
    logger << LOG_DEBUG_POINT << "osi_QT_H1_delta + .5 * oset.Ft1born = " << osi_QT_H1_delta << endl;
    logger << LOG_DEBUG_POINT << "                  .5 * oset.Ft1born = " << .5 * oset.Ft1born << endl;

    if (oset.QT_finalstate_massive_coloured && osi_type_contribution == "CT2" && oset.initial_diag){
      for (int i_c = 0; i_c < osi_QT_correlationoperator.size(); i_c++){
	osi_QT_ME2_loopcf[i_c] += .5 * oset.Ft1born_4correlator[i_c];
	// bug in previous version !!!
	//	osi_QT_ME2_loopcf[i_c] += .5 * oset.Ft1born_4correlator[i_c] * osi_VA_b_ME2;
	////	osi_QT_ME2_loopcf[i_c] += .5 * temp_Ft1born * osi_VA_b_ME2;
	//	QT_ME2_loopcf_shifted[i_c] = osi_QT_ME2_loopcf[i_c] + .5 * temp_Ft1born * osi_VA_b_ME2;
	//	  - 2 * oset.beta0 * log((osi_p_parton[0][1] + osi_p_parton[0][2]).m2() / pow(osi_var_mu_ren, 2)) * osi_QT_ME2_cf[i_c] * osi_VA_b_ME2;
	
	logger << LOG_DEBUG_POINT << "shifted: osi_QT_ME2_loopcf[" << i_c << "] (em = " << osi_QT_correlationoperator[i_c].no_emitter << ", sp = " << osi_QT_correlationoperator[i_c].no_spectator << ") = " << osi_QT_ME2_loopcf[i_c] << endl;
      }
      
      logger << LOG_DEBUG_POINT << "shifted: OL_checksum (1) = " << osi_QT_ME2_loopcf[0] + osi_QT_ME2_loopcf[1] + osi_QT_ME2_loopcf[2] + osi_QT_ME2_loopcf[3] << endl;
      logger << LOG_DEBUG_POINT << "shifted: OL_checksum (2) = " << osi_QT_ME2_loopcf[1] + osi_QT_ME2_loopcf[4] + osi_QT_ME2_loopcf[5] + osi_QT_ME2_loopcf[6] << endl;
      logger << LOG_DEBUG_POINT << "shifted: OL_checksum (3) = " << osi_QT_ME2_loopcf[2] + osi_QT_ME2_loopcf[5] + osi_QT_ME2_loopcf[7] + osi_QT_ME2_loopcf[8] << endl;
      logger << LOG_DEBUG_POINT << "shifted: OL_checksum (4) = " << osi_QT_ME2_loopcf[3] + osi_QT_ME2_loopcf[6] + osi_QT_ME2_loopcf[8] + osi_QT_ME2_loopcf[9] << endl;
      
    }
  }

  }
  else if (osi_type_contribution == "L2VT" || 
	   osi_user_string_value[osi_user_string_map["model"]] == "Bornloop"){

    static double acc;
    double *M2L2;
    M2L2 = new double[5];
  
    static char * OL_mu_ren = stch("muren");
    static char * OL_mu_reg = stch("mureg");
    static char * pole_uv = stch("pole_uv");
    static char * pole_ir1 = stch("pole_ir1");
    static char * pole_ir2 = stch("pole_ir2");
    ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
    ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
    ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
  
    double mu_Q = (osi_p_parton[0][1] + osi_p_parton[0][2]).m();
    ol_setparameter_double(OL_mu_ren, mu_Q);
    ol_setparameter_double(OL_mu_reg, mu_Q);
    static char * fact_uv = stch("fact_uv");
    static char * fact_ir = stch("fact_ir");
    ol_setparameter_double(fact_uv, one);
    ol_setparameter_double(fact_ir, one);

    osi_VA_b_ME2 = 0.;

    double var_mu_ren= (osi_p_parton[0][1] + osi_p_parton[0][2]).m();
    static char * renscale = stch("renscale");
    ol_setparameter_double(renscale, var_mu_ren);
    ol_evaluate_loop2(1, P, M2L2, &acc);
    osi_VA_b_ME2 = M2L2[0];
    logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_L2I = " << setw(23) << setprecision(15) << osi_VA_b_ME2 << endl; 
 
    if (oset.initial_diag){
      if (osi_switch_H1gg){
	osi_QT_A0 = osi_VA_b_ME2;
	osi_QT_A1 = 0.;
	osi_QT_H1_delta = 0.;
      }
      else {
	generic.calculate_H1gg(oset);
	// osi_QT_A0 determined in calculate_H1gg does not contain massive-quark loops.
	// The number of quark flavours in the loop is set by osi.N_f .
      
	//  logger << LOG_DEBUG_VERBOSE << "ratios = " << osi_QT_A0 / osi_VA_b_ME2 << ", " << osi_QT_A1 / M2L1[0] << endl;
	logger << LOG_DEBUG_POINT << "ratio:   born(VVamp) / b_ME2(OL) = " << setw(23) << setprecision(15) << osi_QT_A0 / osi_VA_b_ME2 << endl;

	// A0 and b_ME2 should be identical if same flavour schemes are used.
	// To avoid mismatch between L2VT and L2VA, set:
	// reweighting 2-loop amplitude with mt-dependence (by commenting the following line):
	// osi_QT_H1_delta = osi_QT_H1_delta / osi_VA_b_ME2 * osi_QT_A0;
	// osi_QT_A0 = osi_VA_b_ME2;
	// Changed to OpenLoops result in order to check if the CS-QT difference can be explained from this !!!
      
	logger << LOG_DEBUG_VERBOSE << "born(VVamp) = " << setw(23) << setprecision(15) << osi_QT_A0 << "   born(OL) = " << setw(23) << setprecision(15) << osi_VA_b_ME2 << ", " << "1-loop VVamp = " << setw(23) << setprecision(15) << osi_QT_A1 << endl;
      }
      // cout << "stopping to check" << endl;
      // assert(false);
    }
    else {
      osi_QT_A0 = osi_VA_b_ME2;
      osi_QT_A1 = 0.;
      osi_QT_H1_delta = 0.;
    }

    
    delete [] M2L2;
  }
  
  delete [] P;

  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void calculate_ME2check_VT_QCD(observable_set & oset, call_generic & generic){
  static Logger logger("calculate_ME2check_VT_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  osi_VA_delta_flag = 1;
  string sDelta;

  osi_VA_b_ME2 = 0.;
  osi_VA_V_ME2 = 0.;
  osi_VA_X_ME2 = 0.;
  osi_VA_I_ME2 = 0.;

  osi_QT_A0 = 0.;
  osi_QT_A1 = 0.;
  osi_QT_A2 = 0.;
  osi_QT_H1_delta = 0.;
  osi_QT_H2_delta = 0.;
  for (int i_s = 0; i_s < osi_n_scales_CV; i_s++){osi_VA_X_ME2_CV[i_s] = 0.;}

  if (osi_p_parton[0][0].x0() != 0.){
    ofstream out_comparison;
    out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
    osi_VA_DeltaUV = 0.;
    osi_VA_DeltaIR1 = 0.;
    osi_VA_DeltaIR2 = 0.;
    static char * OL_mu_ren = stch("muren");
    static char * OL_mu_reg = stch("mureg");
    static char * pole_uv = stch("pole_uv");
    static char * pole_ir1 = stch("pole_ir1");
    static char * pole_ir2 = stch("pole_ir2");
    ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
    ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
    ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
    ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
    ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
    calculate_ME2_VT_QCD(oset, generic);

    if (oset.switch_OL){OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));}
    
    out_comparison << "Absolute results: " << endl << endl;
    oset.output_testpoint_VA_result(out_comparison);
    out_comparison << endl;
    out_comparison << "Particle momenta: " << endl << endl;
    output_momenta(out_comparison, oset);
    out_comparison.close();
  }
  //  OLP_PrintParameter(stch(osi_filename_olparameters.c_str()));
 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void calculate_ME2_VT2_QCD(observable_set & oset, call_generic & generic){
  static Logger logger("calculate_ME2_VT2_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (oset.switch_OL){

  static double acc;
  static int n_momentum = 5 * (osi_n_particle + 2);
  static double *P;
  P = new double[n_momentum];
  for (int i = 1; i < osi_p_parton[0].size(); i++){
    P[5 * (i - 1)]     = osi_p_parton[0][i].x0();
    P[5 * (i - 1) + 1] = osi_p_parton[0][i].x1();
    P[5 * (i - 1) + 2] = osi_p_parton[0][i].x2();
    P[5 * (i - 1) + 3] = osi_p_parton[0][i].x3();
    P[5 * (i - 1) + 4] = osi_p_parton[0][i].m();
  }

  // new (20190124):
  if (osi_QT_finalstate_massive_coloured){
    if (oset.initial_pdf_diag || oset.initial_pdf_gq){
      static int n_cc = (osi_n_particle + 2) * (osi_n_particle + 1) / 2;
      static double ewcc;
      static double *M2cc;
      M2cc = new double[n_cc];
      ol_evaluate_cc(1, P, &osi_ME2, M2cc, &ewcc);
      
      for (int i_c = 0; i_c < n_cc; i_c++){
	logger << LOG_DEBUG_VERBOSE << "M2cc[" << i_c << "] = " << setw(23) << setprecision(15) << M2cc[i_c] << endl;
      }

      for (int i_c = 0; i_c < osi_QT_correlationoperator.size(); i_c++){
	//      if (osi_QT_correlationoperator[i_c].type_combination == 1 || osi_QT_correlationoperator[i_c].type_combination == 4){osi_QT_ME2_cf[i_c] = osi_QT_correlationoperator[i_c].charge_factor;}
	if (osi_QT_correlationoperator[i_c].type_combination == 1 || 
	    osi_QT_correlationoperator[i_c].type_combination == 4){
	  osi_QT_ME2_cf[i_c] = osi_QT_correlationoperator[i_c].charge_factor;
	}
	//      else if (osi_QT_correlationoperator[i_c].type_combination == 4){osi_QT_ME2_cf[i_c] = osi_QT_correlationoperator[i_c].charge_factor;}
	else {
	  osi_QT_ME2_cf[i_c] = M2cc[osi_QT_correlationoperator[i_c].no_BLHA_entry] / osi_ME2;
	}
	
	logger << LOG_DEBUG_POINT << "OpenLoops:  QT_ME2_cf[" << i_c << "] (em = " << osi_QT_correlationoperator[i_c].no_emitter << ", sp = " << osi_QT_correlationoperator[i_c].no_spectator << ") = " << setw(23) << setprecision(15) << osi_QT_ME2_cf[i_c] << "   charge_factor = " << setw(23) << setprecision(15) << osi_QT_correlationoperator[i_c].charge_factor << "   no_BLHA_entry = " << osi_QT_correlationoperator[i_c].no_BLHA_entry << endl;
      }
      //    delete [] M2cc;
    }
  }
  
  double *M2L1;
  double *IRL1;
  M2L1 = new double[3];
  IRL1 = new double[3];
  double *M2L2;
  double *IRL2;
  M2L2 = new double[5];
  IRL2 = new double[5];

  double mu_Q = (osi_p_parton[0][1] + osi_p_parton[0][2]).m();
  static double one = 1;

  static char * OL_mu_ren = stch("muren");
  static char * OL_mu_reg = stch("mureg");
    //  static char * renscale = stch("renscale");
  static char * pole_uv = stch("pole_uv");
  static char * pole_ir1 = stch("pole_ir1");
  static char * pole_ir2 = stch("pole_ir2");
  //  ol_setparameter_double(renscale, osi_var_mu_ren);
  ///  ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
  ///  ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
  ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
  ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
  ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
  
  int polenorm_0=0;
  static char * polenorm = stch("polenorm");
  ol_setparameter_int(polenorm, polenorm_0);
  
  //  ol_setparameter_double(stch("renscale"), mu_Q);
  ol_setparameter_double(OL_mu_ren, mu_Q);
  ol_setparameter_double(OL_mu_reg, mu_Q);
  static char * fact_uv = stch("fact_uv");
  static char * fact_ir = stch("fact_ir");
  ol_setparameter_double(fact_uv, one);
  ol_setparameter_double(fact_ir, one);
  //  ol_setparameter_double(stch("fact_uv"), one);
  //  ol_setparameter_double(stch("fact_ir"), one);
  
  //  ol_evaluate_loop2(1, P, M2L2, &acc);
  ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);

  // previous version !!! check if correct !!!
  //  osi_VA_V_ME2 = M2L1[0];

  static double M2L0;
  ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2);
  osi_VA_V_ME2 = M2L1[0] + osi_VA_X_ME2;
  // ??? One could simply use the sum as usually later...



  // double-spin-flip contribution for generic gg-initiated processes:
  if (oset.initial_gg){
    oset.QT_H0_doublespinflip = 0.;

    // Replacing to generic heavy quark mass
    //double m_t = osi_msi.M_t;
    double m_HQ = osi_msi.M[oset.csi->type_parton[0][3]];
    

    callmsqggav_(P, &m_HQ, &oset.QT_H0_doublespinflip);
    //    callmsqGGav_(P, &QT_H0_doublespinflip);
    logger << LOG_DEBUG_POINT << "ME2_B_doublespinflip = " <<  oset.QT_H0_doublespinflip << "   ME2_B = " << osi_VA_b_ME2 << endl;
    oset.QT_H0_doublespinflip = oset.QT_H0_doublespinflip / osi_VA_b_ME2 * pow(osi_alpha_S, 2);
    logger << LOG_DEBUG_POINT << "oset.QT_H0_doublespinflip = " <<  oset.QT_H0_doublespinflip << endl;

    oset.QT_H0_DG = 0.;
    callmsqdgav_(P, &m_HQ, &oset.QT_H0_DG);
    //    callmsqGGav_(P, &QT_H0_doublespinflip);
    logger << LOG_DEBUG_POINT << "ME2_B_oset.QT_H0_DG = " <<  oset.QT_H0_DG << "   ME2_B = " << osi_VA_b_ME2 << endl;
    oset.QT_H0_DG = oset.QT_H0_DG / osi_VA_b_ME2 * pow(osi_alpha_S, 2);
    logger << LOG_DEBUG_POINT << "oset.QT_H0_DG = " << oset.QT_H0_DG << endl;


  }
  
  delete [] M2L1;
  delete [] IRL1;
  delete [] M2L2;
  delete [] IRL2;
  delete [] P;

  }
  
  if (osi_switch_H2){
    osi_QT_A0 = osi_VA_b_ME2;
    osi_QT_A1 = osi_VA_V_ME2;
    osi_QT_H1_delta = osi_QT_A1 / 2 / (osi_alpha_S * inv2pi * osi_QT_A0);
    osi_QT_A2 = 0.;
    osi_QT_H2_delta = 0.;
    //osi_QT_H2_delta = osi_QT_A2 / 4 / (pow(osi_alpha_S * inv2pi, 2) * osi_QT_A0);
  }
  else {
    generic.calculate_H2(oset);
  
    //  logger << LOG_DEBUG_VERBOSE << "ratios = " << osi_QT_A0 / osi_VA_b_ME2 << ", " << osi_QT_A1 / M2L1[0] << endl;
    ///    logger << LOG_DEBUG_VERBOSE << "ratios = " << osi_QT_A0 / osi_VA_b_ME2 << ", " << osi_QT_A1 / osi_VA_V_ME2 << endl;
    ///    cerr << "polylog   A2_1 = " << setw(15) << setprecision(8) << osi_QT_A2 / osi_QT_A1
    ///	 << "     s_part = " << setw(15) << setprecision(8) << osi_p_parton[0][0].m() << endl;
      /*
    cerr << "polylog   A2_1 = " << setw(15) << setprecision(8) << osi_QT_A2 / osi_QT_A1
	 << "     A1_0 = " << setw(15) << setprecision(8) << osi_QT_A1 / osi_QT_A0
	 << "     A20_11 = " << setw(15) << setprecision(8) << osi_QT_A2 * osi_QT_A0 / osi_QT_A1 / osi_QT_A1 << endl;
    cerr << "A2 = " << setw(23) << setprecision(15) << osi_QT_A2
	 << "     A1 = " << setw(23) << setprecision(15) << osi_QT_A1
	 << "     A0 = " << setw(23) << setprecision(15) << osi_QT_A0 << endl;
    cerr << "ratios (A0, A1) = " << setw(23) << setprecision(15) << osi_QT_A0 / osi_VA_b_ME2 << ", " << setw(23) << setprecision(15) << osi_QT_A1 / osi_VA_V_ME2 << endl;
      */
    // needed because of polylog issue in GiNaC (via VVamp) !!!
    /*
    if (abs(osi_QT_A2 / osi_QT_A1) > 5.){
      osi_QT_A2 = 0.;
      osi_QT_H2_delta = 0.;
      (oset.psi->i_nan)++;
    }
    */    
  }

  // new (20190124):
  if (osi_QT_finalstate_massive_coloured){
    if (oset.initial_pdf_diag || oset.initial_pdf_gq){
      
      oset.calculate_Ft1born();
      logger << LOG_DEBUG_VERBOSE << "osi_QT_H1_delta_OL   = " << osi_QT_H1_delta * osi_VA_b_ME2 << endl;
      //  H1 -> H1 - It1  (It1 = 2 Re < M0 | It1 | M0 >) 
      logger << LOG_DEBUG_POINT << "osi_QT_H1_delta                     = " << osi_QT_H1_delta << endl;
      osi_QT_H1_delta += .5 * oset.Ft1born;
      logger << LOG_DEBUG_POINT << "osi_QT_H1_delta + .5 * oset.Ft1born = " << osi_QT_H1_delta << endl;
      logger << LOG_DEBUG_POINT << "                  .5 * oset.Ft1born = " << .5 * oset.Ft1born << endl;
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void calculate_ME2check_VT2_QCD(observable_set & oset, call_generic & generic){
  static Logger logger("calculate_ME2check_VT2_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  osi_VA_delta_flag = 1;
  string sDelta;

  osi_VA_b_ME2 = 0.;
  osi_VA_V_ME2 = 0.;
  osi_VA_X_ME2 = 0.;
  osi_VA_I_ME2 = 0.;

  osi_QT_A0 = 0.;
  osi_QT_A1 = 0.;
  osi_QT_A2 = 0.;
  osi_QT_H1_delta = 0.;
  osi_QT_H2_delta = 0.;
  for (int i_s = 0; i_s < osi_n_scales_CV; i_s++){osi_VA_X_ME2_CV[i_s] = 0.;}

  if (osi_p_parton[0][0].x0() != 0.){
    ofstream out_comparison;
    out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
    osi_VA_DeltaUV = 0.;
    osi_VA_DeltaIR1 = 0.;
    osi_VA_DeltaIR2 = 0.;
    static char * OL_mu_ren = stch("muren");
    static char * OL_mu_reg = stch("mureg");
    static char * pole_uv = stch("pole_uv");
    static char * pole_ir1 = stch("pole_ir1");
    static char * pole_ir2 = stch("pole_ir2");
    ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
    ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
    ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
    ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
    ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
    calculate_ME2_VT2_QCD(oset, generic);

    if (oset.switch_OL){OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));}
    
    out_comparison << "Absolute results: " << endl << endl;
    oset.output_testpoint_VA_result(out_comparison);
    out_comparison << endl;
    out_comparison << "Particle momenta: " << endl << endl;
    output_momenta(out_comparison, oset);
    out_comparison.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
