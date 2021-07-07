#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"
#include "../include/definitions.phasespace.set.cxx"

// different type_corrections could easily be merged !!!

void calculate_ME2_ioperator_VA_QCD(observable_set & oset){
  static Logger logger("calculate_ME2_ioperator_VA_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static vector<double> I_ME2_emitter(osi_VA_ioperator.size());
  if (initialization == 1){
    for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
      osi_VA_ME2_cf[i_a].resize(osi_VA_ioperator[i_a].size());
      osi_VA_I_ME2_cf[i_a].resize(osi_VA_ioperator[i_a].size());
    }
    initialization = 0;
  }

  if (oset.switch_OL == 1){

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

  if ((osi_type_contribution == "VA" ||
       osi_type_contribution == "RVA" ||
       osi_type_contribution == "RVJ") && osi_user_string_value[osi_user_string_map["model"]] != "Bornloop"){
    ol_evaluate_tree(1, P, &osi_VA_b_ME2);
    static int n_cc = (osi_n_particle + 2) * (osi_n_particle + 1) / 2;
    static double ewcc = 0.;
    static double *M2cc;
    M2cc = new double[n_cc];
    ol_evaluate_cc(1, P, &osi_VA_b_ME2, M2cc, &ewcc);
    for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
      for (int j_a = 0; j_a < osi_VA_ioperator[i_a].size(); j_a++){
	osi_VA_ME2_cf[i_a][j_a] = M2cc[osi_VA_ioperator[i_a][j_a].no_BLHA_entry];
	logger << LOG_DEBUG_POINT << "OpenLoops:  VA_ME2_cf[" << i_a << "][" << j_a << "] = " << setw(23) << setprecision(15) << osi_VA_ME2_cf[i_a][j_a] << endl;
      }
    }
    delete [] M2cc;
  }
  else if (osi_type_contribution == "L2VA" || 
	   osi_user_string_value[osi_user_string_map["model"]] == "Bornloop"){
    static int n_cc = (osi_n_particle + 2) * (osi_n_particle + 1) / 2;
    static double ewcc = 0.;
    static double *M2cc;
    M2cc = new double[n_cc];
    ol_evaluate_cc2(1, P, &osi_VA_b_ME2, M2cc, &ewcc);
    for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
      for (int j_a = 0; j_a < osi_VA_ioperator[i_a].size(); j_a++){
	osi_VA_ME2_cf[i_a][j_a] = M2cc[osi_VA_ioperator[i_a][j_a].no_BLHA_entry];
	logger << LOG_DEBUG_POINT << "OpenLoops:  VA_ME2_cf[" << i_a << "][" << j_a << "] = " << setw(23) << setprecision(15) << osi_VA_ME2_cf[i_a][j_a] << endl;
      }
    }
    delete [] M2cc;
  }

  delete [] P;

  }

  if (osi_massive_QCD){oset.calculate_ioperator_QCD_CDST();}
  else {oset.calculate_ioperator_QCD_CS();}
  
  for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
    I_ME2_emitter[i_a] = accumulate(osi_VA_I_ME2_cf[i_a].begin(), osi_VA_I_ME2_cf[i_a].end(), 0.);
  }
  osi_VA_I_ME2 = accumulate(I_ME2_emitter.begin(), I_ME2_emitter.end(), 0.);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void calculate_ME2_ioperator_VA_QEW(observable_set & oset){
  static Logger logger("calculate_ME2_ioperator_VA_QEW");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static vector<double> I_ME2_emitter(osi_VA_ioperator.size());
  if (initialization == 1){
    for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
      osi_VA_ME2_cf[i_a].resize(osi_VA_ioperator[i_a].size());
      osi_VA_I_ME2_cf[i_a].resize(osi_VA_ioperator[i_a].size());
    }
    initialization = 0;
  }

  if (oset.switch_OL == 1){

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

  ol_evaluate_tree(1, P, &osi_VA_b_ME2);

  delete [] P;

  for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
    for (int j_a = 0; j_a < osi_VA_ioperator[i_a].size(); j_a++){
      osi_VA_ME2_cf[i_a][j_a] = osi_VA_ioperator[i_a][j_a].charge_factor() * osi_VA_b_ME2;
    }
  }

  }

  if (osi_massive_QEW){oset.calculate_ioperator_QEW_CDST();}
  else {oset.calculate_ioperator_QEW_CS();}

  for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
    I_ME2_emitter[i_a] = accumulate(osi_VA_I_ME2_cf[i_a].begin(), osi_VA_I_ME2_cf[i_a].end(), 0.);
  }
  osi_VA_I_ME2 = accumulate(I_ME2_emitter.begin(), I_ME2_emitter.end(), 0.);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void calculate_ME2_ioperator_VA_MIX(observable_set & oset){
  static Logger logger("calculate_ME2_ioperator_VA_MIX");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static vector<double> I_ME2_emitter(osi_VA_ioperator.size());
  if (initialization == 1){
    for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
      osi_VA_ME2_cf[i_a].resize(osi_VA_ioperator[i_a].size());
      osi_VA_I_ME2_cf[i_a].resize(osi_VA_ioperator[i_a].size());
    }
    initialization = 0;
  }
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

  ////  ol_evaluate_tree(1, P, &osi_VA_b_ME2);
  static int n_cc = (osi_n_particle + 2) * (osi_n_particle + 1) / 2;
  static double ewcc = 0.;
  static double *M2cc;
  M2cc = new double[n_cc];
  //  ol_evaluate_cc(1, P, M2cc);
  ////  ol_evaluate_cc(1, P, &osi_VA_b_ME2, M2cc, &ewcc);

  int flag_QCD = 0;
  int flag_QEW = 0;
  for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
    for (int j_a = 0; j_a < osi_VA_ioperator[i_a].size(); j_a++){
      if (flag_QCD == 0){
	if (osi_VA_ioperator[i_a][j_a].type_correction() == 1){
	  double dummy_b_ME2;
	  ol_evaluate_cc(osi_VA_ioperator[i_a][j_a].process_id, P, &dummy_b_ME2, M2cc, &ewcc);
	  flag_QCD = 1;
	}
      }
      if (flag_QEW == 0){
	if (osi_VA_ioperator[i_a][j_a].type_correction() == 2){
	  ol_evaluate_tree(osi_VA_ioperator[i_a][j_a].process_id, P, &osi_VA_b_ME2);
	  flag_QEW = 1;
	}
      }
    }
  }

  for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
    if (osi_VA_ioperator[i_a][0].type_correction() == 1){
      for (int j_a = 0; j_a < osi_VA_ioperator[i_a].size(); j_a++){
        osi_VA_ME2_cf[i_a][j_a] = M2cc[osi_VA_ioperator[i_a][j_a].no_BLHA_entry];
      }
    }
    else if (osi_VA_ioperator[i_a][0].type_correction() == 2){
      for (int j_a = 0; j_a < osi_VA_ioperator[i_a].size(); j_a++){
        osi_VA_ME2_cf[i_a][j_a] = osi_VA_ioperator[i_a][j_a].charge_factor() * osi_VA_b_ME2;
      }
    }
  }
  delete [] M2cc;
  delete [] P;
  
  /*
  for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
    for (int j_a = 0; j_a < osi_VA_ioperator[i_a].size(); j_a++){
      logger << LOG_DEBUG_VERBOSE << "osi_VA_ME2_cf[" << i_a << "][" << j_a << "] = " << setw(23) << setprecision(15) << osi_VA_ME2_cf[i_a][j_a] << endl;
    }
  }
  */

  if (osi_massive_QCD){oset.calculate_ioperator_QCD_CDST();}
  else {oset.calculate_ioperator_QCD_CS();}
  
  if (osi_massive_QEW){oset.calculate_ioperator_QEW_CDST();}
  else {oset.calculate_ioperator_QEW_CS();}

  for (int i_a = 0; i_a < osi_VA_ioperator.size(); i_a++){
    I_ME2_emitter[i_a] = accumulate(osi_VA_I_ME2_cf[i_a].begin(), osi_VA_I_ME2_cf[i_a].end(), 0.);
  }
  osi_VA_I_ME2 = accumulate(I_ME2_emitter.begin(), I_ME2_emitter.end(), 0.);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void calculate_ME2_VA_QCD(observable_set & oset){
  static Logger logger("calculate_ME2_VA_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi_switch_VI == 0 || osi_switch_VI == 1){

    if (oset.switch_OL == 1){

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

    if ((osi_type_contribution == "VA" ||
	 osi_type_contribution == "RVA" ||
	 osi_type_contribution == "RVJ") &&
	osi_user_string_value[osi_user_string_map["model"]] != "Bornloop"){

      static char * OL_mu_ren = stch("muren");
      static char * OL_mu_reg = stch("mureg");
      static char * fact_uv = stch("fact_uv");
      static char * fact_ir = stch("fact_ir");
      ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
      ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
      ol_setparameter_double(fact_uv, one);
      ol_setparameter_double(fact_ir, one);

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

      /*
      // Another option:  Set  CT_on = 1 , set  osi_VA_V_ME2 = M2L1[0] - osi_VA_X_ME2;  after CT evaluation (remaining mu_ren-dependent counterterms as before):
      int CT_on = 1; 
      ol_setparameter_int(stch("ct_on"), CT_on); // modification of ...ME2_VA and OpenLoops calls needed !!!
      */    
      ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);
      osi_VA_V_ME2 = M2L1[0];
      ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2);

      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_B  = " << setw(23) << setprecision(15) << osi_VA_b_ME2 << endl; 
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_VX = " << setw(23) << setprecision(15) << osi_VA_V_ME2 + osi_VA_X_ME2 << endl;

      /*
      if (oset.switch_output_comparison){
	int R2_on = 0;
	ol_setparameter_int(stch("r2_on"), R2_on);
	double V_ME2_D4 = 0.;
	ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);
	V_ME2_D4 = M2L1[0];
	logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_D4 = " << setw(23) << setprecision(15) << V_ME2_D4 << endl; 
	double V_ME2_R2 = osi_VA_V_ME2 - V_ME2_D4;
	logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_R2 = " << setw(23) << setprecision(15) << V_ME2_R2 << endl; 
	R2_on = 1;
	ol_setparameter_int(stch("r2_on"), R2_on);
      }
      */
	
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_V  = " << setw(23) << setprecision(15) << osi_VA_V_ME2 << endl; 
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_CT = " << setw(23) << setprecision(15) << osi_VA_X_ME2 << endl; 
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_VX = " << setw(23) << setprecision(15) << osi_VA_V_ME2 + osi_VA_X_ME2 << endl;
 
      if (osi_switch_VI == 0){
	if (osi_VA_DeltaIR1 == 0. && osi_VA_DeltaIR2 == 0.){osi_VA_I_ME2 = IRL1[0];}
	else {osi_VA_I_ME2 = IRL1[0] + osi_VA_DeltaIR1 * IRL1[1] + osi_VA_DeltaIR2 * IRL1[2];}
	//  osi_VA_I_ME2 = i_DeltaIR1 * IR1 + osi_VA_DeltaIR2 * IR2; // !!! I-operator switched off !!!
      }
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_I  = " << setw(23) << setprecision(15) << osi_VA_I_ME2 << endl; 

      if (osi_switch_CV){
	for (int s = 0; s < osi_n_scales_CV; s++){
	  double inv_factor_CV = osi_var_mu_ren_CV[s] / osi_var_mu_ren;
	  ol_setparameter_double(OL_mu_ren, osi_var_mu_ren_CV[s]);
	  ol_setparameter_double(OL_mu_reg, osi_var_mu_ren_CV[s]);
	  ol_setparameter_double(fact_uv, inv_factor_CV);
	  ol_setparameter_double(fact_ir, inv_factor_CV);
	  ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2_CV[s]);
	  logger << LOG_DEBUG_POINT << "OpenLoops:  VA_X_ME2_CV[" << s << "] = " << setw(23) << setprecision(15) << osi_VA_X_ME2_CV[s] + osi_VA_V_ME2 << endl;	 
	  //	  logger << LOG_DEBUG_POINT << "OpenLoops:  VA_X_ME2_CV[" << s << "] = " << osi_VA_X_ME2_CV[s] << endl; 
	  logger << LOG_DEBUG_VERBOSE << "after ol_evaluate_ct s = " << s << endl;
	}
      }
      if (osi_switch_TSV){
	for (int i_v = 0; i_v < osi_max_dyn_ren + 1; i_v++){
	  for (int i_r = 0; i_r < osi_n_scale_dyn_ren[i_v]; i_r++){
	    double inv_factor_TSV = osi_value_scale_ren[0][i_v][i_r] / osi_var_mu_ren;
	    ol_setparameter_double(OL_mu_ren, osi_value_scale_ren[0][i_v][i_r]);
	    ol_setparameter_double(OL_mu_reg, osi_value_scale_ren[0][i_v][i_r]);
	    ol_setparameter_double(fact_uv, inv_factor_TSV);
	    ol_setparameter_double(fact_ir, inv_factor_TSV);
	    ol_evaluate_ct(1, P, &M2L0, &osi_value_ME2term_ren[0][i_v][i_r]);
	    //	    logger << LOG_DEBUG_POINT << "OpenLoops:  VA_X_ME2_TSV[" << i_v << "][" << i_r << "] = " << osi_value_ME2term_ren[0][i_v][i_r] << endl; 
	    logger << LOG_DEBUG_POINT << "OpenLoops:  VA_VX_ME2_TSV[" << i_v << "][" << i_r << "] = " << osi_value_ME2term_ren[0][i_v][i_r] + osi_VA_V_ME2 << endl; 
	    osi_value_ME2term_ren[0][i_v][i_r] += osi_VA_V_ME2 + osi_VA_I_ME2;
	    logger << LOG_DEBUG_VERBOSE << "after ol_evaluate_ct i_v = " << i_v << "   i_r = " << i_r << endl;
	  }
	}
      }
      
      delete [] M2L2;
      delete [] IRL2;
      delete [] M2L1;
      delete [] IRL1;
    }
    
    else if (osi_type_contribution == "L2VA" || 
	     osi_user_string_value[osi_user_string_map["model"]] == "Bornloop"){
      static double acc;
      //      static double M2L0;
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
      ol_evaluate_loop2(1, P, M2L2, &acc);
      
      logger << LOG_DEBUG_VERBOSE << "M2L2[0] = " << M2L2[0] << endl;
      
      osi_VA_b_ME2 = M2L2[0];
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_L2I = " << setw(23) << setprecision(15) << osi_VA_b_ME2 << endl; 

      if (osi_switch_H1gg){
	osi_QT_A0 = 0.;
	osi_QT_A1 = 0.;
	osi_QT_H1_delta = 0.;
      }
      else {
	oset.xmunich->generic.calculate_H1gg(oset);
      }

      if (osi_QT_A0 == 0. && osi_QT_A1 == 0. && osi_QT_H1_delta == 0.){
	osi_VA_V_ME2 = 0.;
      }
      else {
	logger << LOG_DEBUG << "ratios = " << osi_QT_A0 / osi_VA_b_ME2 << endl;
	logger << LOG_DEBUG_VERBOSE << "ratios = " << osi_QT_A0 / osi_VA_b_ME2 << endl;
	
	logger << LOG_DEBUG_VERBOSE << "born VVamp = " << setprecision(15) << setw(23) << osi_QT_A0 << " ,   "
	       << "1-loop VVamp = " << setprecision(15) << setw(23) << osi_QT_A1 << endl;
	logger << LOG_DEBUG_VERBOSE << "born OL    = " << setprecision(15) << setw(23) << osi_VA_b_ME2 << " ,   "
	       << "             = " << setprecision(15) << setw(23) << osi_VA_V_ME2 << endl;

	// reweighting 2-loop amplitude with mt-dependence (from osi_VA_b_ME2):
	osi_VA_V_ME2 = osi_QT_A1 * osi_VA_b_ME2 / osi_QT_A0;

	// no mt-dependence (from osi_VA_b_ME2) in 2-loop amplitude:
	//      osi_VA_V_ME2 = osi_QT_A1;

	// osi.QT_A1 is only the two-loop amplitude (times one-loop...),
	// without reweighting with any loop² result
      }
      
      // temporary solution because mu_reg does not exist as a standard parameter -> use mu_reg = mu_Q = s^ in I-operator !!!
      double save_osi_var_mu_ren = oset.var_mu_ren;
      oset.var_mu_ren = mu_Q;
      calculate_ME2_ioperator_VA_QCD(oset);
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_L2II  = " << setw(23) << setprecision(15) << osi_VA_I_ME2 << "   (corrected for var_rel_alpha_S)" << endl;
      oset.var_mu_ren = save_osi_var_mu_ren;
 
      // beta0 usually not initialized in CS subtraction !!!
      static double beta0 = (33. - 2 * oset.N_f) / 12;
      // scale variation (from mu_ren = s^ = mu_Q to the usual scales...) !!!
      osi_VA_X_ME2 = -oset.csi->order_alpha_s_born * beta0 * log(pow(mu_Q / osi_var_mu_ren, 2)) * osi_VA_b_ME2 * (osi_alpha_S / pi);
      if (osi_switch_CV){
	for (int i_s = 0; i_s < osi_n_scales_CV; i_s++){
	  osi_VA_X_ME2_CV[i_s] = -oset.csi->order_alpha_s_born * beta0 * log(pow(mu_Q / osi_var_mu_ren_CV[i_s], 2)) * osi_VA_b_ME2 * (osi_alpha_S / pi);
	}
      }
      
      if (osi_switch_TSV){
	for (int i_v = 0; i_v < osi_max_dyn_ren + 1; i_v++){
	  for (int i_r = 0; i_r < osi_n_scale_dyn_ren[i_v]; i_r++){
	    osi_value_ME2term_ren[0][i_v][i_r] = osi_VA_V_ME2 + osi_VA_I_ME2 - oset.csi->order_alpha_s_born * beta0 * log(pow(mu_Q / osi_value_scale_ren[0][i_v][i_r], 2)) * osi_VA_b_ME2 * (osi_alpha_S / pi);
	  }
	}
      }

      logger << LOG_DEBUG_VERBOSE
	     << "   V = " << setprecision(15) << setw(23) << osi_VA_V_ME2
	     << "   X = " << setprecision(15) << setw(23) << osi_VA_X_ME2
	     << "   I = " << setprecision(15) << setw(23) << osi_VA_I_ME2 << endl;
     
    }
    delete [] P;

    }
    
  }
  else if (osi_switch_VI == 2){
    logger << LOG_DEBUG_VERBOSE << "OL: osi_VA_I_ME2 = " << osi_VA_I_ME2 << endl;
    //  osi_VA_I_ME2 = 0.;
    //  osi_VA_V_ME2 = 0.;
    //  osi_VA_X_ME2 = 0.;
    //  for (int s = 0; s < osi_n_scales_CV; s++){osi_VA_X_ME2_CV[s] = 0.;}
    //  I-operator evaluation -> osi_VA_I_ME2
    calculate_ME2_ioperator_VA_QCD(oset);
    
    for (int i_v = 0; i_v < osi_max_dyn_ren + 1; i_v++){
      for (int i_r = 0; i_r < osi_n_scale_dyn_ren[i_v]; i_r++){
        osi_value_ME2term_ren[0][i_v][i_r] = osi_VA_I_ME2;
      }
    }
    logger << LOG_DEBUG_VERBOSE << "SK: osi_VA_I_ME2 = " << osi_VA_I_ME2 << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void calculate_ME2_VA_QEW(observable_set & oset){
  static Logger logger("calculate_ME2_VA_QEW");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi_switch_VI == 0 || osi_switch_VI == 1){

    if (oset.switch_OL == 1){

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
    static char * OL_mu_ren = stch("muren");
    static char * OL_mu_reg = stch("mureg");
    static char * fact_uv = stch("fact_uv");
    static char * fact_ir = stch("fact_ir");
	
    ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
    //  ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
    logger << LOG_DEBUG << "osi_user_double_value[osi_user_double_map[mureg]] = " << osi_user_double_value[osi_user_double_map["mu_reg"]] << endl;
    if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
    else if (osi_user_double_value[osi_user_double_map["mu_reg"]] == -1.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
    else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}
    logger << LOG_DEBUG_VERBOSE << "osi_user_double_value[osi_user_double_map[mureg]] = " << osi_user_double_value[osi_user_double_map["mu_reg"]] << endl;
	
    ol_setparameter_double(fact_uv, one);
    ol_setparameter_double(fact_ir, one);

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
    logger << LOG_DEBUG_VERBOSE << "before" << endl;

    /////
    /////    int CT_on = 1; // splitting into V and X does not work properly in present OpenLoops version !!!
    /////    ol_setparameter_int(stch("ct_on"), CT_on);
    /*// for ol_evaluate_ct check !!!
    int CT_on = 0;
    ol_setparameter_int(stch("ct_on"), CT_on);
    *///
   
    /*
    int CT_on = 1; // splitting into V and X does not work properly in present OpenLoops version !!!
    ol_setparameter_int(stch("ct_on"), CT_on);
    ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);
    *//*
    logger << LOG_INFO << "w/ CT   M2L1[0]  = " << M2L1[0] << endl;
    logger << LOG_INFO << "w/ CT   M2L1[1]  = " << M2L1[1] << endl;
    logger << LOG_INFO << "w/ CT   M2L1[2]  = " << M2L1[2] << endl;
    if (osi_VA_DeltaIR1 == 0. && osi_VA_DeltaIR2 == 0.){osi_VA_I_ME2 = IRL1[0];}
    else {osi_VA_I_ME2 = IRL1[0] + osi_VA_DeltaIR1 * IRL1[1] + osi_VA_DeltaIR2 * IRL1[2];}
    logger << LOG_INFO << "V+X_ME2+I_ME2 = " << M2L1[0] + osi_VA_I_ME2 << endl;
    logger.newLine(LOG_INFO);
    ///
    *//*
    osi_VA_V_ME2 = M2L1[0];
    osi_VA_X_ME2 = 0.;
    
    CT_on = 0;
    ol_setparameter_int(stch("ct_on"), CT_on);
    */
    
    // actual version (3 lines):
    ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);
    osi_VA_V_ME2 = M2L1[0];
    ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2);
    
    
    /*
    logger << LOG_INFO << "w/o CT  M2L1[0]  = " << M2L1[0] << endl;
    logger << LOG_INFO << "w/o CT  M2L1[1]  = " << M2L1[1] << endl;
    logger << LOG_INFO << "w/o CT  M2L1[2]  = " << M2L1[2] << endl;
    // I temporary here !!!
    if (osi_VA_DeltaIR1 == 0. && osi_VA_DeltaIR2 == 0.){osi_VA_I_ME2 = IRL1[0];}
    else {osi_VA_I_ME2 = IRL1[0] + osi_VA_DeltaIR1 * IRL1[1] + osi_VA_DeltaIR2 * IRL1[2];}
    logger << LOG_INFO << "osi_VA_V_ME2  = " << osi_VA_V_ME2 << endl;
    logger << LOG_INFO << "osi_VA_X_ME2  = " << osi_VA_X_ME2 << endl;
    logger << LOG_INFO << "osi_VA_I_ME2  = " << osi_VA_I_ME2 << endl;
    logger << LOG_INFO << "V+X+I_ME2     = " << osi_VA_V_ME2 + osi_VA_X_ME2 + osi_VA_I_ME2 << endl;
    */
    
    /////    osi_VA_X_ME2 = 0.;
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_B  = " << setw(23) << setprecision(15) << osi_VA_b_ME2 << endl; 
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_VX = " << setw(23) << setprecision(15) << osi_VA_V_ME2 + osi_VA_X_ME2 << endl;
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_V  = " << setw(23) << setprecision(15) << osi_VA_V_ME2 << endl; 
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_CT = " << setw(23) << setprecision(15) << osi_VA_X_ME2 << endl; 
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_VX = " << setw(23) << setprecision(15) << osi_VA_V_ME2 + osi_VA_X_ME2 << endl;
 

      if (osi_switch_VI == 0){
      // I-operator in OpenLoops not defined correctly for external photons:
      // check if still correct without external photons !!!
      
      /*
      calculate_ME2_ioperator_VA_QEW(oset);
      logger << LOG_INFO << setw(10) << oset.psi->i_acc << "   I-operator (SK)   = " << setprecision(15) << osi_VA_I_ME2 << endl;
      osi_VA_I_ME2 = 0.;
      */
      
      //      if (osi_VA_delta_flag){
      ///            if (0 == 1){
	if (oset.user.switch_map["I_Munich"]){
	  calculate_ME2_ioperator_VA_QEW(oset);
	  logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_I  = " << setw(23) << setprecision(15) << osi_VA_I_ME2 << "   (from Munich)" << endl;
	}
	else {
	  // original version (2 lines):
	  if (osi_VA_DeltaIR1 == 0. && osi_VA_DeltaIR2 == 0.){osi_VA_I_ME2 = IRL1[0];}
	  else {osi_VA_I_ME2 = IRL1[0] + osi_VA_DeltaIR1 * IRL1[1] + osi_VA_DeltaIR2 * IRL1[2];}
	  logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_I  = " << setw(23) << setprecision(15) << osi_VA_I_ME2 << "   (from OpenLoops)" << endl;
	}

	
	//  osi_VA_I_ME2 = i_DeltaIR1 * IR1 + osi_VA_DeltaIR2 * IR2; // !!! I-operator switched off !!!

	/*
	// alternative version using MUNICH I-operator
	calculate_ME2_ioperator_VA_QEW(oset);
	*/
	
	/////	logger << LOG_DEBUG << setw(10) << oset.psi->i_acc << "   I-operator (OL)   = " << setprecision(15) << osi_VA_I_ME2 << endl;
	/*///
      }
      else {
	calculate_ME2_ioperator_VA_QEW(oset);
	logger << LOG_DEBUG << setw(10) << oset.psi->i_acc << "   I-operator (SK)   = " << setprecision(15) << osi_VA_I_ME2 << endl;
      }
	*/	///
	// !!!

      /* // commented out now: should only check if poles fit with SK I-operator
	// Shift V by the difference between I_OL and I_SK
	double I_ME2_OL = osi_VA_I_ME2;
	osi_VA_I_ME2 = 0.;
	calculate_ME2_ioperator_VA_QEW(oset);
	logger << LOG_DEBUG << setw(10) << oset.psi->i_acc << "   I-operator (SK)   = " << setprecision(15) << osi_VA_I_ME2 << endl;
	osi_VA_V_ME2 += I_ME2_OL - osi_VA_I_ME2;
      */ //
	/*
      osi_VA_I_ME2 = 0.;
      calculate_ME2_ioperator_VA_QEW(oset);
      logger << LOG_INFO << setw(10) << oset.psi->i_acc << "   I-operator (SK)   = " << setprecision(15) << osi_VA_I_ME2 << endl;
      */
	/*
	logger << LOG_INFO << "osi_VA_delta_flag = " << osi_VA_delta_flag << endl;
	logger << LOG_INFO << "osi_VA_DeltaUV    = " << osi_VA_DeltaUV << endl;
	logger << LOG_INFO << "osi_VA_DeltaIR1   = " << osi_VA_DeltaIR1 << endl;
	logger << LOG_INFO << "osi_VA_DeltaIR2   = " << osi_VA_DeltaIR2 << endl;
	logger << LOG_INFO << "I-operator (OL)   = " << setprecision(15) << osi_VA_I_ME2 << endl;
	osi_VA_I_ME2 = 0.;
	calculate_ME2_ioperator_VA_QEW(oset);
	logger << LOG_INFO << "I-operator (SK)   = " << setprecision(15) << osi_VA_I_ME2 << endl;
	*/
      }
      else if (osi_switch_VI == 1){
	osi_VA_I_ME2 = 0.;
      }
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_I  = " << setw(23) << setprecision(15) << osi_VA_I_ME2 << endl; 

    //    cout << right << setw(5) << "I1" << " = " << setprecision(15) << setw(23) << left << IRL1[1] << endl;
    //    cout << right << setw(5) << "I2" << " = " << setprecision(15) << setw(23) << left << IRL1[2] << endl;

      // Shifted to above:
    ///    ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2);
    /////    osi_VA_X_ME2 = 0.;
    if (osi_switch_CV){
      for (int s = 0; s < osi_n_scales_CV; s++){
        double inv_factor_CV = osi_var_mu_ren_CV[s] / osi_var_mu_ren;
        ol_setparameter_double(OL_mu_ren, osi_var_mu_ren_CV[s]);
	if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren_CV[s]);}
	else if (osi_user_double_value[osi_user_double_map["mu_reg"]] == -1.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
	else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}

        ol_setparameter_double(fact_uv, inv_factor_CV);
        ol_setparameter_double(fact_ir, inv_factor_CV);
	ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2_CV[s]);
	logger << LOG_DEBUG_POINT << "OpenLoops:  VA_X_ME2_CV[" << s << "] = " << setw(23) << setprecision(15) << osi_VA_X_ME2_CV[s] + osi_VA_V_ME2 << endl;	 
	/////	osi_VA_X_ME2_CV[s] = 0.;
	logger << LOG_DEBUG_VERBOSE << "after ol_evaluate_ct s = " << s << endl;
      }
    }
    if (osi_switch_TSV){
      for (int i_v = 0; i_v < osi_max_dyn_ren + 1; i_v++){
	for (int i_r = 0; i_r < osi_n_scale_dyn_ren[i_v]; i_r++){
	  double inv_factor_TSV = osi_value_scale_ren[0][i_v][i_r] / osi_var_mu_ren;
	  ol_setparameter_double(OL_mu_ren, osi_value_scale_ren[0][i_v][i_r]);
	  if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_value_scale_ren[0][i_v][i_r]);}
	  else if (osi_user_double_value[osi_user_double_map["mu_reg"]] == -1.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
	  else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}
	  
	  ol_setparameter_double(fact_uv, inv_factor_TSV);
	  ol_setparameter_double(fact_ir, inv_factor_TSV);
	  ol_evaluate_ct(1, P, &M2L0, &osi_value_ME2term_ren[0][i_v][i_r]);
	  logger << LOG_DEBUG_POINT << "OpenLoops:  VA_VX_ME2_TSV[" << i_v << "][" << i_r << "] = " << osi_value_ME2term_ren[0][i_v][i_r] + osi_VA_V_ME2 << endl; 
	  osi_value_ME2term_ren[0][i_v][i_r] += osi_VA_V_ME2 + osi_VA_I_ME2;
	  /////	  osi_value_ME2term_ren[0][i_v][i_r] = osi_VA_V_ME2 + osi_VA_I_ME2;
	  logger << LOG_DEBUG_VERBOSE << "after ol_evaluate_ct i_v = " << i_v << "   i_r = " << i_r << endl;
	}
      }
    }
    delete [] M2L2;
    delete [] IRL2;
    delete [] M2L1;
    delete [] IRL1;
    delete [] P;
    }
  }
  else if (osi_switch_VI == 2){
    logger << LOG_DEBUG_VERBOSE << "OL: osi_VA_I_ME2 = " << osi_VA_I_ME2 << endl;
    //  osi_VA_I_ME2 = 0.;
    //  osi_VA_V_ME2 = 0.;
    //  osi_VA_X_ME2 = 0.;
    //  for (int s = 0; s < osi_n_scales_CV; s++){osi_VA_X_ME2_CV[s] = 0.;}
    //  I-operator evaluation -> osi_VA_I_ME2
    calculate_ME2_ioperator_VA_QEW(oset);
    for (int i_v = 0; i_v < osi_max_dyn_ren + 1; i_v++){
      for (int i_r = 0; i_r < osi_n_scale_dyn_ren[i_v]; i_r++){
        osi_value_ME2term_ren[0][i_v][i_r] = osi_VA_I_ME2;
      }
    }
    logger << LOG_DEBUG_VERBOSE << "SK: osi_VA_I_ME2 = " << osi_VA_I_ME2 << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void calculate_ME2_VA_MIX(observable_set & oset){
  static Logger logger("calculate_ME2_VA_MIX");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi_switch_VI == 0 || osi_switch_VI == 1){
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
    //    static char * renscale = stch("renscale");
    static char * OL_mu_ren = stch("muren");
    static char * OL_mu_reg = stch("mureg");
    static char * fact_uv = stch("fact_uv");
    static char * fact_ir = stch("fact_ir");
    ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
    if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
    else if (osi_user_double_value[osi_user_double_map["mu_reg"]] == -1.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
    else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}

    ol_setparameter_double(fact_uv, one);
    ol_setparameter_double(fact_ir, one);

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
    ol_evaluate_full(1, P, &osi_VA_b_ME2, M2L1, IRL1, M2L2, IRL2, &acc);
 
    osi_VA_V_ME2 = M2L1[0];
    if (osi_switch_VI == 0){
      if (osi_VA_DeltaIR1 == 0. && osi_VA_DeltaIR2 == 0.){osi_VA_I_ME2 = IRL1[0];}
      else {osi_VA_I_ME2 = IRL1[0] + osi_VA_DeltaIR1 * IRL1[1] + osi_VA_DeltaIR2 * IRL1[2];}
      //  osi_VA_I_ME2 = i_DeltaIR1 * IR1 + osi_VA_DeltaIR2 * IR2; // !!! I-operator switched off !!!
    }

    ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2);
    if (osi_switch_CV){
      for (int s = 0; s < osi_n_scales_CV; s++){
        double inv_factor_CV = osi_var_mu_ren_CV[s] / osi_var_mu_ren;
        ol_setparameter_double(OL_mu_ren, osi_var_mu_ren_CV[s]);
	if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren_CV[s]);}
	else if (osi_user_double_value[osi_user_double_map["mu_reg"]] == -1.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
	else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}

        ol_setparameter_double(fact_uv, inv_factor_CV);
        ol_setparameter_double(fact_ir, inv_factor_CV);
	ol_evaluate_ct(1, P, &M2L0, &osi_VA_X_ME2_CV[s]);
	logger << LOG_DEBUG_VERBOSE << "after ol_evaluate_ct s = " << s << endl;
      }
    }
    if (osi_switch_TSV){
      for (int i_v = 0; i_v < osi_max_dyn_ren + 1; i_v++){
	for (int i_r = 0; i_r < osi_n_scale_dyn_ren[i_v]; i_r++){
	  double inv_factor_TSV = osi_value_scale_ren[0][i_v][i_r] / osi_var_mu_ren;
	  ol_setparameter_double(OL_mu_ren, osi_value_scale_ren[0][i_v][i_r]);
	  if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_value_scale_ren[0][i_v][i_r]);}
	  else if (osi_user_double_value[osi_user_double_map["mu_reg"]] == -1.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
	  else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}
	  ol_setparameter_double(fact_uv, inv_factor_TSV);
	  ol_setparameter_double(fact_ir, inv_factor_TSV);
	  ol_evaluate_ct(1, P, &M2L0, &osi_value_ME2term_ren[0][i_v][i_r]);
	  osi_value_ME2term_ren[0][i_v][i_r] += osi_VA_V_ME2 + osi_VA_I_ME2;
	  logger << LOG_DEBUG_VERBOSE << "after ol_evaluate_ct i_v = " << i_v << "   i_r = " << i_r << endl;
	}
      }
    }
    delete [] M2L2;
    delete [] IRL2;
    delete [] M2L1;
    delete [] IRL1;
    delete [] P;
  }
  else if (osi_switch_VI == 2){
    logger << LOG_DEBUG_VERBOSE << "OL: osi_VA_I_ME2 = " << osi_VA_I_ME2 << endl;
    //  osi_VA_I_ME2 = 0.;
    //  osi_VA_V_ME2 = 0.;
    //  osi_VA_X_ME2 = 0.;
    //  for (int s = 0; s < osi_n_scales_CV; s++){osi_VA_X_ME2_CV[s] = 0.;}
    //  I-operator evaluation -> osi_VA_I_ME2
    calculate_ME2_ioperator_VA_MIX(oset);
    for (int i_v = 0; i_v < osi_max_dyn_ren + 1; i_v++){
      for (int i_r = 0; i_r < osi_n_scale_dyn_ren[i_v]; i_r++){
        osi_value_ME2term_ren[0][i_v][i_r] = osi_VA_I_ME2;
      }
    }
    logger << LOG_DEBUG_VERBOSE << "SK: osi_VA_I_ME2 = " << osi_VA_I_ME2 << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void calculate_ME2check_VA_QCD(observable_set & oset){
  static Logger logger("calculate_ME2check_VA_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (oset.switch_output_comparison){

    oset.xmunich->generic.phasespacepoint_psp(oset);

    if (osi_p_parton[0][0].x0() == 0.){
      if (oset.switch_OL){oset.testpoint_from_OL_rambo();}
    }
   
    if (osi_p_parton[0][0].x0() != 0.){
    
      ofstream out_comparison;
      logger << LOG_DEBUG_VERBOSE << "osi_filename_comparison = " << osi_filename_comparison << endl;
      out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  

      perform_event_selection(oset, oset.xmunich->generic);
      
      if (oset.cut_ps[0] == -1){
	out_comparison << "Phase-space 0 is cut. -> Default scale is used." << endl;
	for (int sd = 1; sd < oset.max_dyn_ren + 1; sd++){
	  for (int ss = 0; ss < oset.n_scale_dyn_ren[sd]; ss++){
	    oset.value_scale_ren[0][sd][ss] = 100.;
	  }
	}
	for (int sd = 1; sd < oset.max_dyn_fact + 1; sd++){
	  for (int ss = 0; ss < oset.n_scale_dyn_fact[sd]; ss++){
	    oset.value_scale_fact[0][sd][ss] = 100.;
	  }
	}
	oset.cut_ps[0] = 0;
      }
      else {
	// Maybe not valid any longer:
	// 1 - Testpoint is calculated at basic fixed scale (prefactor * scale_ren).
	// 2 - Testpoint is calculated at output scale.
	oset.xmunich->generic.calculate_dynamic_scale(0, oset);
	oset.xmunich->generic.calculate_dynamic_scale_TSV(0, oset);
	oset.determine_scale();
	/////  }
      }

      static char * OL_mu_ren = stch("muren");
      static char * OL_mu_reg = stch("mureg");
      
      static char * pole_uv = stch("pole_uv");
      static char * pole_ir1 = stch("pole_ir1");
      static char * pole_ir2 = stch("pole_ir2");
      //  static char * me_cache = stch("me_cache");
      
      static char * OL_alpha_s = stch("alpha_s");
      
      osi_VA_delta_flag = 1;
      string s_Delta;
      osi_VA_b_ME2 = 0.;
      osi_VA_V_ME2 = 0.;
      osi_VA_X_ME2 = 0.;
      osi_VA_I_ME2 = 0.;
      osi_VA_X_ME2_CV.resize(osi_n_scales_CV, 0.);
      
    /*
    if (osi_p_parton[0][0].x0() != 0.){
    */
      if (oset.switch_output_comparison == 2){
	// With OpenLoops, MUNICH does not reset the value of alpha_S for different scales, but adds the corresponding relative factors.
	// For the test-point output, the value of alpha_S is set in OpenLoops here:
	//      osi_alpha_S = osi_alpha_S * osi_var_rel_alpha_S;
	//      osi_alpha_S = osi_alpha_S * pow(osi_var_rel_alpha_S, double(oset.csi->contribution_order_alpha_s - 1) / oset.csi->contribution_order_alpha_s); 
	//      osi_alpha_S = oset.var_alpha_S_reference; 
	//      ol_setparameter_double(OL_alpha_s, osi_alpha_S);
	ol_setparameter_double(OL_alpha_s, oset.var_alpha_S_reference);
      }
      /*
      ofstream out_comparison;
      logger << LOG_DEBUG_VERBOSE << "osi_filename_comparison = " << osi_filename_comparison << endl;
      out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
    */

      double i_Delta;
      osi_VA_DeltaUV = 0.;
      osi_VA_DeltaIR1 = 0.;
      osi_VA_DeltaIR2 = 0.;
      
      logger << LOG_DEBUG << "osi_VA_DeltaUV  = " << osi_VA_DeltaUV << endl;
      logger << LOG_DEBUG << "osi_VA_DeltaIR1 = " << osi_VA_DeltaIR1 << endl;
      logger << LOG_DEBUG << "osi_VA_DeltaIR2 = " << osi_VA_DeltaIR2 << endl;
      
      if (oset.switch_OL){
      ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
      ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
      ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
      ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
      ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
      //  ol_setparameter_double(me_cache, 0);
      }
    
      if (oset.switch_OL){OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));}

      calculate_ME2_VA_QCD(oset);

      out_comparison << "Settings: " << endl << endl;
      out_comparison << setw(12) << "  mu_reg  = " << setprecision(15) << setw(23) << osi_var_mu_ren << "  ( = mu_ren by default)" << endl;
      out_comparison << setw(12) << "  mu_ren  = " << setprecision(15) << setw(23) << osi_var_mu_ren << "  " << "alpha_S(mu_ren) = " << osi_alpha_S << endl;
      out_comparison << setw(12) << "  mu_fact = " << setprecision(15) << setw(23) << osi_var_mu_fact << endl;
      out_comparison << endl;
      
      out_comparison << "Absolute results: " << endl << endl;
      oset.output_testpoint_VA_result(out_comparison);
      double OL_I_ME2 = osi_VA_I_ME2;
      if (osi_switch_VI != 2){
	osi_VA_I_ME2 = 0.;
	calculate_ME2_ioperator_VA_QCD(oset);
	oset.output_testpoint_VA_ioperator(out_comparison);
	osi_VA_I_ME2 = OL_I_ME2;
	//    double SK_I_ME2 = osi_VA_I_ME2;
	//    oset.output_testpoint_VA_ioperator(out_comparison);
      }
      osi_VA_V_ME2 = osi_VA_V_ME2 / osi_VA_b_ME2;
      osi_VA_X_ME2 = osi_VA_X_ME2 / osi_VA_b_ME2;
      osi_VA_I_ME2 = osi_VA_I_ME2 / osi_VA_b_ME2;
      //    osi_VA_b_ME2 = osi_VA_b_ME2;
      out_comparison << "Relative results (corrections devided by ME2_born, ME2_born divided by coupling constants): " << endl << endl;
      oset.output_testpoint_VA_result(out_comparison);
      out_comparison << endl;
      out_comparison << "Particle momenta: " << endl << endl;
      output_momenta(out_comparison, oset);


      out_comparison << "Numerical check of (UV and IR) finiteness: " << endl << endl;

      int temp_switch_check_IRneqUV = 1;

      string set_OL_model = "";
      for (int i_o = 0; i_o < osi_OL_parameter.size(); i_o++){
	if (osi_OL_parameter[i_o] == "model"){set_OL_model = osi_OL_value[i_o];}
      }
      if (set_OL_model == "heft"){temp_switch_check_IRneqUV = 0;}
      
      if (temp_switch_check_IRneqUV){
      s_Delta = "Delta_UV";
      for (int i = 0; i < 3; i++){
	i_Delta = (double(i) - 1.) * 1.;
	osi_VA_DeltaUV = i_Delta;
	
	logger << LOG_DEBUG << "osi_VA_DeltaUV  = " << osi_VA_DeltaUV << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR1 = " << osi_VA_DeltaIR1 << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR2 = " << osi_VA_DeltaIR2 << endl;
	
	if (oset.switch_OL){
	ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
	ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
	ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
	ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
	ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
	//  ol_setparameter_double(me_cache, 0);
	}

	calculate_ME2_VA_QCD(oset);
	oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
      }
      osi_VA_DeltaUV = 0.;
      out_comparison << endl;

      s_Delta = "Delta_IR_1";
      for (int i = 0; i < 3; i++){
	i_Delta = (double(i) - 1.) * 1.;
	osi_VA_DeltaIR1 = i_Delta;
	
	logger << LOG_DEBUG << "osi_VA_DeltaUV  = " << osi_VA_DeltaUV << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR1 = " << osi_VA_DeltaIR1 << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR2 = " << osi_VA_DeltaIR2 << endl;
	
	if (oset.switch_OL){
	ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
	ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
	ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
	ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
	ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
	//  ol_setparameter_double(me_cache, 0);
	}

	calculate_ME2_VA_QCD(oset);
	oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
      }
      osi_VA_DeltaIR1 = 0.;
      out_comparison << endl;
      }
    
    
      s_Delta = "Delta_IR_2";
      for (int i = 0; i < 3; i++){
	i_Delta = (double(i) - 1.) * 1.;
	
	osi_VA_DeltaIR2 = i_Delta;
	
	logger << LOG_DEBUG << "osi_VA_DeltaUV  = " << osi_VA_DeltaUV << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR1 = " << osi_VA_DeltaIR1 << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR2 = " << osi_VA_DeltaIR2 << endl;
	
	if (oset.switch_OL){
	ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
	ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
	ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
	ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
	ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
	//  ol_setparameter_double(me_cache, 0);
	}

	calculate_ME2_VA_QCD(oset);
	oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
      }
      osi_VA_DeltaIR2 = 0.;
      out_comparison << endl;
      
      s_Delta = "Delta_UV = Delta_IR_1";
      for (int i = 0; i < 3; i++){
	i_Delta = (double(i) - 1.) * 1.;
	osi_VA_DeltaUV = i_Delta;
	osi_VA_DeltaIR1 = i_Delta;
	
	logger << LOG_DEBUG << "osi_VA_DeltaUV  = " << osi_VA_DeltaUV << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR1 = " << osi_VA_DeltaIR1 << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR2 = " << osi_VA_DeltaIR2 << endl;
	
	if (oset.switch_OL){
	ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
	ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
	ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
	ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
	ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
	//  ol_setparameter_double(me_cache, 0);
	}

	calculate_ME2_VA_QCD(oset);
	oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
      }
      osi_VA_DeltaUV = 0.;
      osi_VA_DeltaIR1 = 0.;
      out_comparison << endl;
  

      logger << LOG_DEBUG << "osi_VA_DeltaUV  = " << osi_VA_DeltaUV << endl;
      logger << LOG_DEBUG << "osi_VA_DeltaIR1 = " << osi_VA_DeltaIR1 << endl;
      logger << LOG_DEBUG << "osi_VA_DeltaIR2 = " << osi_VA_DeltaIR2 << endl;
      
      if (oset.switch_OL){
      ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
      ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
      ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
      ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
      ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
      //  ol_setparameter_double(me_cache, 0);
      }
      
      calculate_ME2_VA_QCD(oset);
  
      osi_VA_delta_flag = 0;
      out_comparison.close();

      if (oset.switch_output_comparison == 2){
	// With OpenLoops, MUNICH does not reset the value of alpha_S for different scales, but adds the corresponding relative factors.
	// For the test-point output, the value of alpha_S is set in OpenLoops and reset to its standard value here:
	//      osi_alpha_S = osi_alpha_S / osi_var_rel_alpha_S;
	///      osi_alpha_S = osi_alpha_S / pow(osi_var_rel_alpha_S, double(oset.csi->contribution_order_alpha_s - 1) / oset.csi->contribution_order_alpha_s); 
	ol_setparameter_double(OL_alpha_s, osi_alpha_S);
      }
    }
    else {
      if (oset.switch_OL){OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));}
    }
    
    oset.switch_output_comparison = 0;
  }
  else {
    if (oset.switch_OL){OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));}
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void calculate_ME2check_VA_QEW(observable_set & oset){
  static Logger logger("calculate_ME2check_VA_QEW");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (oset.switch_output_comparison){

    oset.xmunich->generic.phasespacepoint_psp(oset);
    
    if (osi_p_parton[0][0].x0() == 0.){
      if (oset.switch_OL){oset.testpoint_from_OL_rambo();}
    }
   
    if (osi_p_parton[0][0].x0() != 0.){
      ofstream out_comparison;
      out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
      
      perform_event_selection(oset, oset.xmunich->generic);
      if (oset.cut_ps[0] == -1){
	out_comparison << "Phase-space 0 is cut. -> Default scale is used." << endl;
	for (int sd = 1; sd < oset.max_dyn_ren + 1; sd++){
	  for (int ss = 0; ss < oset.n_scale_dyn_ren[sd]; ss++){
	    oset.value_scale_ren[0][sd][ss] = 100.;
	  }
	}
	for (int sd = 1; sd < oset.max_dyn_fact + 1; sd++){
	  for (int ss = 0; ss < oset.n_scale_dyn_fact[sd]; ss++){
	    oset.value_scale_fact[0][sd][ss] = 100.;
	  }
	}
	osi_var_mu_ren = 100.;
	osi_var_mu_fact = 100.;
	
	oset.cut_ps[0] = 0;
      }
      else {
	// Maybe not valid any longer:
	// 1 - Testpoint is calculated at basic fixed scale (prefactor * scale_ren).
	// 2 - Testpoint is calculated at output scale.
	oset.xmunich->generic.calculate_dynamic_scale(0, oset);
	oset.xmunich->generic.calculate_dynamic_scale_TSV(0, oset);
	oset.determine_scale();
	/////  }
      }
      
      static char * OL_mu_ren = stch("muren");
      static char * OL_mu_reg = stch("mureg");
      
      static char * pole_uv = stch("pole_uv");
      static char * pole_ir1 = stch("pole_ir1");
      static char * pole_ir2 = stch("pole_ir2");
      //  static char * me_cache = stch("me_cache");
      
      static char * OL_alpha_s = stch("alpha_s");
      
      osi_VA_delta_flag = 1;
      string s_Delta;
      osi_VA_b_ME2 = 0.;
      osi_VA_V_ME2 = 0.;
      osi_VA_X_ME2 = 0.;
      osi_VA_I_ME2 = 0.;
      osi_VA_X_ME2_CV.resize(osi_n_scales_CV, 0.);
    
    /*
    if (osi_p_parton[0][0].x0() != 0.){
      ofstream out_comparison;
      out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
    */
    
      double coefficient_B = 0.;
      double coefficient_VF = 0.;
      double coefficient_V1 = 0.;
      double coefficient_V2 = 0.;
      
      if (oset.switch_output_comparison == 2){
	ol_setparameter_double(OL_alpha_s, oset.var_alpha_S_reference);
      }

      double i_Delta;
      osi_VA_DeltaUV = 0.;
      osi_VA_DeltaIR1 = 0.;
      osi_VA_DeltaIR2 = 0.;

      ///    ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
      ///    if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
      ///    else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}

      if (oset.switch_OL){
      ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
      ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
      ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
      ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
      ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
      //  ol_setparameter_double(me_cache, 0);
      }

    
      if (oset.switch_OL){OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));}

      calculate_ME2_VA_QEW(oset);

      coefficient_B = osi_VA_b_ME2;
      coefficient_VF = osi_VA_V_ME2 + osi_VA_X_ME2;

      out_comparison << "Settings: " << endl << endl;
      out_comparison << setw(12) << "  mu_reg  = " << setprecision(15) << setw(23) << osi_var_mu_ren << "  ( = mu_ren by default)" << endl;
      out_comparison << setw(12) << "  mu_ren  = " << setprecision(15) << setw(23) << osi_var_mu_ren << "  " << "alpha_S(mu_ren) = " << osi_alpha_S << endl;
      out_comparison << setw(12) << "  mu_fact = " << setprecision(15) << setw(23) << osi_var_mu_fact << endl;
      out_comparison << endl;
    
      out_comparison << "Absolute results: " << endl << endl;
      oset.output_testpoint_VA_result(out_comparison);
      double OL_I_ME2 = osi_VA_I_ME2;
      if (osi_switch_VI != 2){
	osi_VA_I_ME2 = 0.;
	calculate_ME2_ioperator_VA_QEW(oset);
	//      double SK_I_ME2 = osi_VA_I_ME2;
	oset.output_testpoint_VA_ioperator(out_comparison);
	osi_VA_I_ME2 = OL_I_ME2;
      }
      /* // not needed, thus output switched off
	 osi_VA_V_ME2 = osi_VA_V_ME2 / osi_VA_b_ME2;
	 osi_VA_X_ME2 = osi_VA_X_ME2 / osi_VA_b_ME2;
	 osi_VA_I_ME2 = osi_VA_I_ME2 / osi_VA_b_ME2;
	 //    osi_VA_b_ME2 = osi_VA_b_ME2;
	 out_comparison << "Relative results (corrections devided by ME2_born, ME2_born divided by coupling constants): " << endl << endl;
	 oset.output_testpoint_VA_result(out_comparison);
	 out_comparison << endl;
      */
      out_comparison << "Particle momenta: " << endl << endl;
      output_momenta(out_comparison, oset);
      
      out_comparison << "Numerical check of (UV and IR) finiteness: " << endl << endl;

      int temp_switch_check_IRneqUV = 1;

      if (temp_switch_check_IRneqUV){
      s_Delta = "Delta_UV";
      for (int i = 0; i < 3; i++){
	i_Delta = (double(i) - 1.) * 1.;
	osi_VA_DeltaUV = i_Delta;
	
	logger << LOG_DEBUG << "osi_VA_DeltaUV  = " << osi_VA_DeltaUV << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR1 = " << osi_VA_DeltaIR1 << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR2 = " << osi_VA_DeltaIR2 << endl;
	
	if (oset.switch_OL){
	  ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
	  ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
	  ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
	  ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
	  ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
	  //  ol_setparameter_double(me_cache, 0);
	}
	
	calculate_ME2_VA_QEW(oset);
	// use SK I-operator (needed for external photons by now) !!!
	///  calculate_ME2_ioperator_VA_QEW(oset);
	oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
      }
      osi_VA_DeltaUV = 0.;
      out_comparison << endl;
      
      s_Delta = "Delta_IR_1";
      for (int i = 0; i < 3; i++){
	i_Delta = (double(i) - 1.) * 1.;
	osi_VA_DeltaIR1 = i_Delta;
	
	logger << LOG_DEBUG << "osi_VA_DeltaUV  = " << osi_VA_DeltaUV << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR1 = " << osi_VA_DeltaIR1 << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR2 = " << osi_VA_DeltaIR2 << endl;
	
	if (oset.switch_OL){
	ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
	ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
	ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
	ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
	ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
	//  ol_setparameter_double(me_cache, 0);
	}

	calculate_ME2_VA_QEW(oset);
	// use SK I-operator (needed for external photons by now) !!!
	///  calculate_ME2_ioperator_VA_QEW(oset);
	oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
      }
      coefficient_V1 = osi_VA_V_ME2 + osi_VA_X_ME2 - coefficient_VF;
      osi_VA_DeltaIR1 = 0.;
      out_comparison << endl;
      }
    
      s_Delta = "Delta_IR_2";
      for (int i = 0; i < 3; i++){
	i_Delta = (double(i) - 1.) * 1.;
      
	osi_VA_DeltaIR2 = i_Delta;
	
	logger << LOG_DEBUG << "osi_VA_DeltaUV  = " << osi_VA_DeltaUV << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR1 = " << osi_VA_DeltaIR1 << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR2 = " << osi_VA_DeltaIR2 << endl;
	
	if (oset.switch_OL){
	ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
	ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
	ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
	ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
	ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
	//  ol_setparameter_double(me_cache, 0);
	}

	calculate_ME2_VA_QEW(oset);
	// use SK I-operator (needed for external photons by now) !!!
	///  calculate_ME2_ioperator_VA_QEW(oset);
	oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
      }
      coefficient_V2 = osi_VA_V_ME2 + osi_VA_X_ME2 - coefficient_VF;
    
      osi_VA_DeltaIR2 = 0.;
      out_comparison << endl;
      s_Delta = "Delta_UV = Delta_IR_1";
      for (int i = 0; i < 3; i++){
	i_Delta = (double(i) - 1.) * 1.;
	osi_VA_DeltaUV = i_Delta;
	osi_VA_DeltaIR1 = i_Delta;
	
	logger << LOG_DEBUG << "osi_VA_DeltaUV  = " << osi_VA_DeltaUV << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR1 = " << osi_VA_DeltaIR1 << endl;
	logger << LOG_DEBUG << "osi_VA_DeltaIR2 = " << osi_VA_DeltaIR2 << endl;
	
	if (oset.switch_OL){
	ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
	ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
	ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
	ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
	ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
	//  ol_setparameter_double(me_cache, 0);
	}

	calculate_ME2_VA_QEW(oset);
	// use SK I-operator (needed for external photons by now) !!!
	///  calculate_ME2_ioperator_VA_QEW(oset);
	oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
      }
      if (!temp_switch_check_IRneqUV){
	coefficient_V1 = osi_VA_V_ME2 + osi_VA_X_ME2 - coefficient_VF;
      }
      osi_VA_DeltaUV = 0.;
      osi_VA_DeltaIR1 = 0.;
      out_comparison << endl;

      logger << LOG_DEBUG << "osi_VA_DeltaUV  = " << osi_VA_DeltaUV << endl;
      logger << LOG_DEBUG << "osi_VA_DeltaIR1 = " << osi_VA_DeltaIR1 << endl;
      logger << LOG_DEBUG << "osi_VA_DeltaIR2 = " << osi_VA_DeltaIR2 << endl;

      if (oset.switch_OL){
      ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
      ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);
      ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
      ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
      ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
      //  ol_setparameter_double(me_cache, 0);
      }

      calculate_ME2_VA_QEW(oset);
      // in order to reset Delta-dependent SK I-operator settings !!!
      ///  calculate_ME2_ioperator_VA_QEW(oset);
      osi_VA_delta_flag = 0;
    
      out_comparison << endl;
      out_comparison << endl;
      out_comparison << "Coefficients:" << endl;
      out_comparison << endl;
      out_comparison << right << setw(5) << "B" << " = " << setprecision(15) << setw(23) << left << coefficient_B << endl;
      out_comparison << right << setw(5) << "VF" << " = " << setprecision(15) << setw(23) << left << coefficient_VF << endl;
      out_comparison << right << setw(5) << "V1" << " = " << setprecision(15) << setw(23) << left << coefficient_V1 << endl;
      out_comparison << right << setw(5) << "V2" << " = " << setprecision(15) << setw(23) << left << coefficient_V2 << endl;
      out_comparison << endl;
      out_comparison << endl;
      
      out_comparison << "#" << right
		     << setprecision(6) << setw(11) << "lam"
		     << setprecision(15) << setw(23) << "B"
		     << setprecision(15) << setw(23) << "VF"
		     << setprecision(15) << setw(23) << "V1"
		     << setprecision(15) << setw(23) << "V2"
		     << setprecision(15) << setw(23) << "VF/(B*alpha)"
		     << setprecision(15) << setw(23) << "V1/(B*alpha)"
		     << setprecision(15) << setw(23) << "V2/(B*alpha)"
		     << endl;
      
      out_comparison << right
		     << setprecision(7) << setw(12) << oset.msi.alpha_e / 0.0075552541674291547
		     << setprecision(15) << setw(23) << coefficient_B
		     << setprecision(15) << setw(23) << coefficient_VF
		     << setprecision(15) << setw(23) << coefficient_V1
		     << setprecision(15) << setw(23) << coefficient_V2 
		     << setprecision(15) << setw(23) << coefficient_VF / coefficient_B / oset.msi.alpha_e
		     << setprecision(15) << setw(23) << coefficient_V1 / coefficient_B / oset.msi.alpha_e
		     << setprecision(15) << setw(23) << coefficient_V2 / coefficient_B / oset.msi.alpha_e
		     << endl;
      
      out_comparison.close();
    
      if (oset.switch_output_comparison == 2){
	ol_setparameter_double(OL_alpha_s, osi_alpha_S);
      }
    }
    else {
      if (oset.switch_OL){OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));}
    }
    
    oset.switch_output_comparison = 0;
  }
  else {
    if (oset.switch_OL){OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));}
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void calculate_ME2check_VA_MIX(observable_set & oset){
  static Logger logger("calculate_ME2check_VA_MIX");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  
  //    static char * renscale = stch("renscale");
  static char * OL_mu_ren = stch("muren");
  static char * OL_mu_reg = stch("mureg");

  static char * pole_uv = stch("pole_uv");
  static char * pole_ir1 = stch("pole_ir1");
  static char * pole_ir2 = stch("pole_ir2");
  //  static char * me_cache = stch("me_cache");

  static char * stability_mode = stch("stability_mode");
  static char * redlib1 = stch("redlib1");
  static char * redlib2 = stch("redlib2");

  osi_VA_delta_flag = 1;
  string s_Delta;
  osi_VA_b_ME2 = 0.;
  osi_VA_V_ME2 = 0.;
  osi_VA_X_ME2 = 0.;
  osi_VA_I_ME2 = 0.;
  osi_VA_X_ME2_CV.resize(osi_n_scales_CV, 0.);
  if (osi_p_parton[0][0].x0() != 0.){
    ofstream out_comparison;
    out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
    double i_Delta;
    osi_VA_DeltaUV = 0.;
    osi_VA_DeltaIR1 = 0.;
    osi_VA_DeltaIR2 = 0.;

    ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
    if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
    else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}

    ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
    ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
    ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
    //  ol_setparameter_double(me_cache, 0);
    calculate_ME2_VA_MIX(oset);
    //    OLP_PrintParameter(stch("olparameters.txt"));
    if (oset.switch_OL){OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));}
    //    double ME2 = 0;
    out_comparison << "Absolute results: " << endl << endl;
    oset.output_testpoint_VA_result(out_comparison);
    double OL_I_ME2 = osi_VA_I_ME2;
    if (osi_switch_VI != 2){
      osi_VA_I_ME2 = 0.;
      calculate_ME2_ioperator_VA_MIX(oset);
     //      double SK_I_ME2 = osi_VA_I_ME2;
      oset.output_testpoint_VA_ioperator(out_comparison);
      osi_VA_I_ME2 = OL_I_ME2;
    }
    osi_VA_V_ME2 = osi_VA_V_ME2 / osi_VA_b_ME2;
    osi_VA_X_ME2 = osi_VA_X_ME2 / osi_VA_b_ME2;
    osi_VA_I_ME2 = osi_VA_I_ME2 / osi_VA_b_ME2;
    //    osi_VA_b_ME2 = osi_VA_b_ME2;
    out_comparison << "Relative results (corrections devided by ME2_born, ME2_born divided by coupling constants): " << endl << endl;
    oset.output_testpoint_VA_result(out_comparison);
    out_comparison << endl;
    out_comparison << "Particle momenta: " << endl << endl;
    output_momenta(out_comparison, oset);
    out_comparison << "On-shell-projected particle momenta: " << endl << endl;
    output_momenta(out_comparison, oset);
    out_comparison << "Numerical check of (UV and IR) finiteness: " << endl << endl;

    // Delta_UV =/= Delta_IR not allowed in CutTOols !!!
    int set_OL_stability_mode;
    ol_getparameter_int(stability_mode, &set_OL_stability_mode);
    int temp_OL_stability_mode = set_OL_stability_mode;
    int set_OL_redlib1;
    ol_getparameter_int(redlib1, &set_OL_redlib1);
    int temp_OL_redlib1 = set_OL_redlib1;
    int set_OL_redlib2;
    ol_getparameter_int(redlib2, &set_OL_redlib2);
    int temp_OL_redlib2 = set_OL_redlib2;
    //    char* set_OL_model;
    //    ol_getparameter_string("model", set_OL_model);
    string set_OL_model;
    for (int i_o = 0; i_o < osi_OL_parameter.size(); i_o++){
      if (osi_OL_parameter[i_o] == "model"){set_OL_model = osi_OL_value[i_o];}
    }
    string temp_OL_model = set_OL_model;

    logger << LOG_DEBUG << "set_OL_stability_mode  = " << set_OL_stability_mode << endl;
    logger << LOG_DEBUG << "set_OL_redlib1         = " << set_OL_redlib1 << endl;
    logger << LOG_DEBUG << "set_OL_redlib2         = " << set_OL_redlib2 << endl;
    logger << LOG_DEBUG << "set_OL_model           = " << set_OL_model << endl;

    int temp_switch_check_IRneqUV = 1;

    if (set_OL_model == "heft"){temp_switch_check_IRneqUV = 0;}
    else if (set_OL_stability_mode < 20){temp_switch_check_IRneqUV = 0;}
    else {
      if (set_OL_stability_mode == 23){
	temp_OL_stability_mode = 21;
	if      (set_OL_redlib1 == 5 && set_OL_redlib2 == 1){temp_OL_redlib1 = 7;}
	else if (set_OL_redlib1 == 5 && set_OL_redlib2 == 7){temp_OL_redlib1 = 1;}
	else if (set_OL_redlib2 == 5 && set_OL_redlib1 == 1){temp_OL_redlib2 = 7;}
	else if (set_OL_redlib2 == 5 && set_OL_redlib1 == 7){temp_OL_redlib2 = 1;}

	ol_setparameter_int(stch(stability_mode), temp_OL_stability_mode);
	ol_setparameter_int(stch(redlib1), temp_OL_redlib1);
	ol_setparameter_int(stch(redlib2), temp_OL_redlib2);
      }
      logger << LOG_DEBUG << "temp_OL_stability_mode = " << temp_OL_stability_mode << endl;
      logger << LOG_DEBUG << "temp_OL_redlib1        = " << temp_OL_redlib1 << endl;
      logger << LOG_DEBUG << "temp_OL_redlib2        = " << temp_OL_redlib2 << endl;
      logger << LOG_DEBUG << "temp_OL_model          = " << temp_OL_model << endl;
    }

    if (temp_switch_check_IRneqUV){
    s_Delta = "Delta_UV";
    for (int i = 0; i < 3; i++){
      i_Delta = (double(i) - 1.) * 1.;
      osi_VA_DeltaUV = i_Delta;

      ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
      if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
      else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}

      ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
      ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
      ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
      //  ol_setparameter_double(me_cache, 0);

      calculate_ME2_VA_MIX(oset);
      oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
    }
    osi_VA_DeltaUV = 0.;
    out_comparison << endl;

    s_Delta = "Delta_IR_1";
    for (int i = 0; i < 3; i++){
      i_Delta = (double(i) - 1.) * 1.;
      osi_VA_DeltaIR1 = i_Delta;

      ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
      if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
      else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}

      ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
      ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
      ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
      //  ol_setparameter_double(me_cache, 0);
      calculate_ME2_VA_MIX(oset);
      oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
    }
    osi_VA_DeltaIR1 = 0.;
    out_comparison << endl;
    }

    if (set_OL_stability_mode != temp_OL_stability_mode){
      ol_setparameter_int(stability_mode, set_OL_stability_mode);
      ol_setparameter_int(redlib1, set_OL_redlib1);
      ol_setparameter_int(redlib2, set_OL_redlib2);
    }

    s_Delta = "Delta_IR_2";
    for (int i = 0; i < 3; i++){
      i_Delta = (double(i) - 1.) * 1.;
      osi_VA_DeltaIR2 = i_Delta;

      ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
      if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
      else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}

      ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
      ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
      ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
      //  ol_setparameter_double(me_cache, 0);
      calculate_ME2_VA_MIX(oset);
      oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
    }
    osi_VA_DeltaIR2 = 0.;
    out_comparison << endl;
    s_Delta = "Delta_UV = Delta_IR_1";
    for (int i = 0; i < 3; i++){
      i_Delta = (double(i) - 1.) * 1.;
      osi_VA_DeltaUV = i_Delta;
      osi_VA_DeltaIR1 = i_Delta;

      ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
      if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
      else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}

      ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
      ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
      ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
      //  ol_setparameter_double(me_cache, 0);
      calculate_ME2_VA_MIX(oset);
      oset.output_testpoint_VA_Delta(out_comparison, i, i_Delta, s_Delta);
    }
    osi_VA_DeltaUV = 0.;
    osi_VA_DeltaIR1 = 0.;
    out_comparison << endl;

    ol_setparameter_double(OL_mu_ren, osi_var_mu_ren);
    if (osi_user_double_value[osi_user_double_map["mu_reg"]] == 0.){ol_setparameter_double(OL_mu_reg, osi_var_mu_ren);}
    else {ol_setparameter_double(OL_mu_reg, osi_user_double_value[osi_user_double_map["mu_reg"]]);}

    ol_setparameter_double(pole_uv, osi_VA_DeltaUV);
    ol_setparameter_double(pole_ir1, osi_VA_DeltaIR1);
    ol_setparameter_double(pole_ir2, osi_VA_DeltaIR2);
    //  ol_setparameter_double(me_cache, 0);
    calculate_ME2_VA_MIX(oset);
    osi_VA_delta_flag = 0;
    out_comparison.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




