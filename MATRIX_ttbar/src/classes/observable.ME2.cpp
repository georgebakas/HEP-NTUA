#include "../include/classes.cxx"
//#include "../include/definitions.cxx"
//#include "../include/definitions.observable.set.cxx"

void observable_set::calculate_ME2_born(){
  static Logger logger("observable_set::calculate_ME2_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (switch_OL == 1){

  static int n_momentum = 5 * (n_particle + 2);
  static double *P;
  P = new double[n_momentum];
  for (int i = 1; i < p_parton[0].size(); i++){
    P[5 * (i - 1)]     = p_parton[0][i].x0();
    P[5 * (i - 1) + 1] = p_parton[0][i].x1();
    P[5 * (i - 1) + 2] = p_parton[0][i].x2();
    P[5 * (i - 1) + 3] = p_parton[0][i].x3();
    P[5 * (i - 1) + 4] = p_parton[0][i].m();
  }

  if ((csi->type_contribution == "born" ||
       csi->type_contribution == "RT" ||
       csi->type_contribution == "RJ" ||
       csi->type_contribution == "CT" ||
       csi->type_contribution == "CJ") && user.string_value[user.string_map["model"]] != "Bornloop"){
    
    if (csi->type_contribution == "CT" && QT_finalstate_massive_coloured){
      static int n_cc = (n_particle + 2) * (n_particle + 1) / 2;
      static double ewcc;
      static double *M2cc;
      M2cc = new double[n_cc];
      ol_evaluate_cc(1, P, &value_ME2term[0], M2cc, &ewcc);
      ME2 = value_ME2term[0];
      for (int i_c = 0; i_c < QT_correlationoperator.size(); i_c++){
	if (QT_correlationoperator[i_c].type_combination == 1){QT_ME2_cf[i_c] = QT_correlationoperator[i_c].charge_factor;}
	else {QT_ME2_cf[i_c] = M2cc[QT_correlationoperator[i_c].no_BLHA_entry] / ME2;} 
	logger << LOG_DEBUG_VERBOSE << "QT_ME2_cf[" << i_c << "] = " << QT_ME2_cf[i_c] << "   charge_factor = " << QT_correlationoperator[i_c].charge_factor << "   no_BLHA_entry = " << QT_correlationoperator[i_c].no_BLHA_entry << endl;
      }
      delete [] M2cc;
    }
    else {
      ol_evaluate_tree(1, P, &value_ME2term[0]);
      ME2 = value_ME2term[0];
      
      if (!(check_vanishing_ME2_end)){
	if (ME2 == 0.){flag_vanishing_ME2 = 1;}
	else {
	  //      check_vanishing_ME2_born(oset);
	  static double b_ME2;
	  static int n_cc = (n_particle + 2) * (n_particle + 1) / 2;
	  static double ewcc;
	  static double *M2cc;
	  M2cc = new double[n_cc];
	  //  ol_evaluate_tree(1, P, &b_ME2);
	  flag_vanishing_ME2 = 0;
	  ol_evaluate_cc(1, P, &b_ME2, M2cc, &ewcc);
	  logger << LOG_DEBUG_VERBOSE << "b_ME2 = " << b_ME2 << endl;
	  for (int i_c = 0; i_c < n_cc; i_c++){
	    //	    logger << LOG_DEBUG_VERBOSE << "M2cc[" << setw(2) << i_c << "] = " << M2cc[i_c] << endl;
	    logger << LOG_DEBUG << "OpenLoops:  M2cc[" << i_c << "] = " << M2cc[i_c] << endl;
	    if (abs(M2cc[i_c]) > 1.e12 * abs(b_ME2)){
	      logger << LOG_DEBUG << "b_ME2 = 0. due to numerical cancellations!" << endl;
	      flag_vanishing_ME2 = 1;
	      break;
	    }
	  }
	  delete [] M2cc;
	}
      }

      logger << LOG_DEBUG << "ME2_OL  = " << ME2 <<endl;
    }
  }
  else if (csi->type_contribution == "loop" ||
	   csi->type_contribution == "L2I" ||
	   csi->type_contribution == "L2CT" ||
	   csi->type_contribution == "L2CJ" ||
	   csi->type_contribution == "L2RT" || 
	   csi->type_contribution == "L2RJ" || 
	   user.string_value[user.string_map["model"]] == "Bornloop"){
    double *M2L2;
    M2L2 = new double[5];
    static double one = 1;
    static double acc;
    static char * renscale = stch("renscale");
    static char * fact_uv = stch("fact_uv");
    static char * fact_ir = stch("fact_ir");
    ol_setparameter_double(renscale, var_mu_ren);
    ol_setparameter_double(fact_uv, one);
    ol_setparameter_double(fact_ir, one);
    ol_evaluate_loop2(1, P, M2L2, &acc);
    ME2 = M2L2[0];
    logger << LOG_DEBUG_VERBOSE << "M2L2[0] = " << M2L2[0] << endl;
    //  VA_b_ME2 = M2L2[0];
    value_ME2term[0] = M2L2[0];
    delete [] M2L2;

    logger << LOG_DEBUG << "ME2_OL  = " << ME2 <<endl;
  }

  delete [] P;

  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::calculate_ME2check_born(){
  static Logger logger("observable_set::calculate_ME2check_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  
  if (switch_output_comparison){
    
    xmunich->generic.phasespacepoint_psp(*this);

    if (p_parton[0][0].x0() == 0.){
      if (switch_OL){testpoint_from_OL_rambo();}
    }

    perform_event_selection(*this, xmunich->generic);
    if (cut_ps[0] == -1){
      cut_ps[0] = 0;
    }
    else {
      // 1 - Testpoint is calculated at basic fixed scale (prefactor * scale_ren).
      // 2 - Testpoint is calculated at output scale.
      xmunich->generic.calculate_dynamic_scale(0, *this);
      xmunich->generic.calculate_dynamic_scale_TSV(0, *this);
      determine_scale();
      /////  }
    }
    {
      //	if (p_parton[0][0].x0() != 0.){
      ofstream out_comparison;
      out_comparison.open(filename_comparison.c_str(), ofstream::out | ofstream::app);  
      
      calculate_ME2_born();

      output_testpoint_born(out_comparison);
      
      out_comparison << endl << "Particle momenta: " << endl << endl;
      output_momenta(out_comparison);
      
      int i_a = 0;
      out_comparison << endl;
      for (int sd = 1; sd < max_dyn_ren + 1; sd++){
	for (int ss = 0; ss < n_scale_dyn_ren[sd]; ss++){
	  out_comparison << "value_scale_ren[" << i_a << "][" << sd << "][" << ss << "] = " << value_scale_ren[i_a][sd][ss] << endl;
	}
      }
      for (int sd = 1; sd < max_dyn_fact + 1; sd++){
	for (int ss = 0; ss < n_scale_dyn_fact[sd]; ss++){
	  out_comparison << "value_scale_fact[" << i_a << "][" << sd << "][" << ss << "] = " << value_scale_fact[i_a][sd][ss] << endl;
	}
      }
      
      out_comparison.close();
    }
    
  }
  //    OLP_PrintParameter(stch("log/olparameters." + name_process + ".txt"));
  if (switch_OL){OLP_PrintParameter(stch("log/olparameters." + name_process + ".txt"));}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::check_vanishing_ME2_born(){
  static Logger logger("observable_set::check_vanishing_ME2_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int n_momentum = 5 * (n_particle + 2);
  static double *P;
  P = new double[n_momentum];
  for (int i = 1; i < p_parton[0].size(); i++){
    P[5 * (i - 1)]     = p_parton[0][i].x0();
    P[5 * (i - 1) + 1] = p_parton[0][i].x1();
    P[5 * (i - 1) + 2] = p_parton[0][i].x2();
    P[5 * (i - 1) + 3] = p_parton[0][i].x3();
    P[5 * (i - 1) + 4] = p_parton[0][i].m();
  }
  static double b_ME2;
  static int n_cc = (n_particle + 2) * (n_particle + 1) / 2;
  static double ewcc;
  static double *M2cc;
  M2cc = new double[n_cc];
  //  ol_evaluate_tree(1, P, &b_ME2);
  ol_evaluate_cc(1, P, &b_ME2, M2cc, &ewcc);
  logger << LOG_DEBUG << "b_ME2 = " << b_ME2 << endl;
  for (int i_c = 0; i_c < n_cc; i_c++){
    logger << LOG_DEBUG << "M2cc[" << setw(2) << i_c << "] = " << M2cc[i_c] << endl;
    if (abs(M2cc[i_c]) > 1.e12 * abs(b_ME2)){
      logger << LOG_DEBUG << "b_ME2 = 0. due to numerical cancellations!" << endl;
      int_end = 1;

      exit(0);
    }
  }

  delete [] M2cc;
  delete [] P;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::calculate_ME2_loop(){
  static Logger logger("observable_set::calculate_ME2_loop");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double *M2L2;
  M2L2 = new double[5];
  static double one = 1;
  static double acc;
  static int n_momentum = 5 * (n_particle + 2);
  static double *P;
  P = new double[n_momentum];
  for (int i = 1; i < p_parton[0].size(); i++){
    P[5 * (i - 1)]     = p_parton[0][i].x0();
    P[5 * (i - 1) + 1] = p_parton[0][i].x1();
    P[5 * (i - 1) + 2] = p_parton[0][i].x2();
    P[5 * (i - 1) + 3] = p_parton[0][i].x3();
    P[5 * (i - 1) + 4] = p_parton[0][i].m();
  }
  static char * renscale = stch("renscale");
  static char * fact_uv = stch("fact_uv");
  static char * fact_ir = stch("fact_ir");
  ol_setparameter_double(renscale, var_mu_ren);
  ol_setparameter_double(fact_uv, one);
  ol_setparameter_double(fact_ir, one);
  ol_evaluate_loop2(1, P, M2L2, &acc);
  ME2 = M2L2[0];
  //  VA_b_ME2 = M2L2[0];
  value_ME2term[0] = M2L2[0];
  delete [] M2L2;
  delete [] P;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


/*
void observable_set::calculate_ME2check_loop(){
  static Logger logger("observable_set::calculate_ME2check_loop");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  VA_delta_flag = 1;
  string sDelta;
  VA_b_ME2 = 0.;
  VA_V_ME2 = 0.;
  VA_X_ME2 = 0.;
  VA_I_ME2 = 0.;
  if (p_parton[0][0].x0() != 0.){
    ofstream out_comparison;
    out_comparison.open(filename_comparison.c_str(), ofstream::out | ofstream::app);  
    VA_DeltaUV = 0.;
    VA_DeltaIR1 = 0.;
    VA_DeltaIR2 = 0.;
    static char * renscale = stch("renscale");
    static char * pole_uv = stch("pole_uv");
    static char * pole_ir1 = stch("pole_ir1");
    static char * pole_ir2 = stch("pole_ir2");
    ol_setparameter_double(renscale, var_mu_ren);
    ol_setparameter_double(pole_uv, VA_DeltaUV);
    ol_setparameter_double(pole_ir1, VA_DeltaIR1);
    ol_setparameter_double(pole_ir2, VA_DeltaIR2);
    calculate_ME2_loop();
    //    double ME2;
    out_comparison << "Absolute results: " << endl << endl;
    output_result_VA(out_comparison, VA_b_ME2, VA_V_ME2, VA_X_ME2, VA_I_ME2);
    out_comparison << endl;
    out_comparison << "Particle momenta: " << endl << endl;
    output_momenta(out_comparison);
    out_comparison << "On-shell-projected particle momenta: " << endl << endl;
    output_momenta(out_comparison);
    out_comparison.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/

void observable_set::calculate_ME2check_loop(){
  static Logger logger("observable_set::calculate_ME2check_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (p_parton[0][0].x0() != 0.){
    ofstream out_comparison;
    out_comparison.open(filename_comparison.c_str(), ofstream::out | ofstream::app);  

    calculate_ME2_loop();
    if (switch_OL){OLP_PrintParameter(stch("log/olparameters." + name_process + ".txt"));}
    //    OLP_PrintParameter(stch("log/olparameters." + name_process + ".txt"));

    output_testpoint_born(out_comparison);

    out_comparison << endl << "Particle momenta: " << endl << endl;
    output_momenta(out_comparison);
    out_comparison.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

/*
void calculate_ME2check_loop(double & mu_ren, vector<double> & mu_ren_CV, double & rel_alpha_S, vector<double> & rel_alpha_S_CV, observable_set & oset, int & delta_flag){
  static Logger logger("calculate_ME2check_loop");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  delta_flag = 1;
  string sDelta;
  double b_ME2 = 0., V_ME2 = 0., X_ME2 = 0., I_ME2 = 0.;
  vector<double> X_ME2_CV(osi_n_scales_CV);
  if (osi_p_parton[0][0].x0() != 0.){
    ofstream out_comparison;
    out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
    double DeltaUV = 0.;
    double DeltaIR1 = 0.;
    double DeltaIR2 = 0.;
    static char * renscale = stch("renscale");
    static char * pole_uv = stch("pole_uv");
    static char * pole_ir1 = stch("pole_ir1");
    static char * pole_ir2 = stch("pole_ir2");
    ol_setparameter_double(renscale, mu_ren);
    ol_setparameter_double(pole_uv, DeltaUV);
    ol_setparameter_double(pole_ir1, DeltaIR1);
    ol_setparameter_double(pole_ir2, DeltaIR2);
    calculate_ME2_loop(mu_ren, mu_ren_CV, rel_alpha_S, rel_alpha_S_CV, V_ME2, DeltaUV, DeltaIR1, DeltaIR2, delta_flag, oset);
    out_comparison << "Absolute results: " << endl << endl;
    output_result_VA(out_comparison, b_ME2, V_ME2, X_ME2, I_ME2);
    out_comparison << endl;
    out_comparison << "Particle momenta: " << endl << endl;
    output_momenta(out_comparison, oset);
    out_comparison << "On-shell-projected particle momenta: " << endl << endl;
    output_momenta(out_comparison, oset);
    out_comparison.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/
