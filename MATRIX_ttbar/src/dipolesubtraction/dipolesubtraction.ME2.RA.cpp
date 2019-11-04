#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"
#include "../include/definitions.phasespace.set.cxx"

// different type_corrections could maybe be merged !!!

void calculate_ME2_RA_QCD(observable_set & oset){
  static Logger logger("calculate_ME2_QCD_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi_cut_ps[0] > -1){

    if (oset.switch_OL == 1){
    
    static int n_momentum = 5 * (osi_n_particle + 2);
    double *P;
    P = new double[n_momentum];
    for (int i = 1; i < osi_p_parton[0].size(); i++){
      logger << LOG_DEBUG_VERBOSE << "i = " << i << endl;
      P[5 * (i - 1)]     = osi_p_parton[0][i].x0();
      P[5 * (i - 1) + 1] = osi_p_parton[0][i].x1();
      P[5 * (i - 1) + 2] = osi_p_parton[0][i].x2();
      P[5 * (i - 1) + 3] = osi_p_parton[0][i].x3();
      P[5 * (i - 1) + 4] = osi_p_parton[0][i].m();
    }
    //    logger << LOG_DEBUG << "osi_type_contribution = " << osi_type_contribution << "   model = " << osi_user_string_value[osi_user_string_map["model"]] << endl;

    if ((osi_type_contribution == "RA" ||
	 osi_type_contribution == "RRA") && osi_user_string_value[osi_user_string_map["model"]] != "Bornloop"){
      ol_evaluate_tree(osi_RA_dipole[0].process_id, P, &osi_value_ME2term[0]);
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_R  = " << setw(23) << setprecision(15) << osi_value_ME2term[0] << endl; 

      if (!(osi_check_vanishing_ME2_end)){
	static double b_ME2;
	static int n_cc = (osi_n_particle + 2) * (osi_n_particle + 1) / 2;
	static double ewcc;
	static double *M2cc;
	M2cc = new double[n_cc];
	osi_flag_vanishing_ME2 = 0;
	ol_evaluate_cc(1, P, &b_ME2, M2cc, &ewcc);
	logger << LOG_DEBUG_VERBOSE << "b_ME2 = " << b_ME2 << endl;
	for (int i_c = 0; i_c < n_cc; i_c++){
	  logger << LOG_DEBUG_VERBOSE << "M2cc[" << setw(2) << i_c << "] = " << M2cc[i_c] << endl;
	  if (abs(M2cc[i_c]) > 1.e12 * abs(b_ME2)){
	    logger << LOG_DEBUG << "b_ME2 = 0. due to numerical cancellations!" << endl;
	    osi_flag_vanishing_ME2 = 1;
	    break;
	  }
	}
	delete [] M2cc;
      }
    }
    else if (osi_type_contribution == "L2RA" || 
	       osi_user_string_value[osi_user_string_map["model"]] == "Bornloop"){
      logger << LOG_DEBUG << "L2RA" << endl;
      double *M2L2;
      M2L2 = new double[5];
      static double one = 1;
      static double acc;
      static char * renscale = stch("renscale");
      static char * fact_uv = stch("fact_uv");
      static char * fact_ir = stch("fact_ir");
      ol_setparameter_double(renscale, osi_var_mu_ren);
      ol_setparameter_double(fact_uv, one);
      ol_setparameter_double(fact_ir, one);
      logger << LOG_DEBUG_VERBOSE << "osi_RA_dipole[0].process_id = " << osi_RA_dipole[0].process_id << endl;
      ol_evaluate_loop2(osi_RA_dipole[0].process_id, P, M2L2, &acc);
      osi_value_ME2term[0] = M2L2[0];
      logger << LOG_DEBUG_POINT << "OpenLoops:  ME2_R  = " << setw(23) << setprecision(15) << osi_value_ME2term[0] << endl; 

      // This is probably not helpful for loop-induced processes (unstable ME'2 might be set to 0):
      /*
      osi_flag_vanishing_ME2 = 0;
      if (osi_value_ME2term[0] == 0.){
	osi_flag_vanishing_ME2 = 1;
	logger << LOG_DEBUG << "b_ME2 = 0. due to numerical cancellations!" << endl;
      }
      */
      
      delete [] M2L2;
    }
    else {
      logger << LOG_FATAL << "Should not happen!" << endl;
    }
    delete [] P;

    }
  }
  else {osi_value_ME2term[0] = 0.;}

  for (int i_a = 1; i_a < osi_n_ps; i_a++){
    if (osi_cut_ps[i_a] >= 0){
      if (osi_RA_dipole[i_a].massive() == 0){      
	if      (osi_RA_dipole[i_a].type_dipole() == 1){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ij_k(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 2){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ij_a(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 3){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ai_k(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 5){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ai_b(i_a);}
	else {cout << "Should not happen!" << endl;}
      }
      else {
	if      (osi_RA_dipole[i_a].type_dipole() == 1){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ij_k_massive(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 2){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ij_a_massive(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 3){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ai_k(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 5){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ai_b(i_a);}
      }
    }
    else {osi_value_ME2term[i_a] = 0.;}
  }
  osi_RA_ME2 = osi_value_ME2term;
  /*
  for (int i_a = 0; i_a < osi_n_ps; i_a++){
    logger << LOG_DEBUG << "osi_cut_ps[" << setw(2) << i_a << "] = " << setw(3) << osi_cut_ps[i_a] << "   osi_RA_ME2[" << setw(2) << i_a << "] = " << setw(23) << setw(15) << osi_RA_ME2[i_a] << endl;
  }
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void calculate_ME2check_RA_QCD(phasespace_set & psi, observable_set & oset){
  static Logger logger("calculate_ME2check_QCD_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (oset.switch_output_comparison){
    
    oset.xmunich->generic.phasespacepoint_psp(oset);

    if (osi_p_parton[0][0].x0() == 0.){
      if (oset.switch_OL){oset.testpoint_from_OL_rambo();}
    }

    if (osi_p_parton[0][0].x0() != 0.){
      ofstream out_comparison;
      out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
      for (int i_p = 0; i_p < osi_p_parton[0].size(); i_p++){psi_xbp_all[0][intpow(2, i_p - 1)] = osi_p_parton[0][i_p];}
      for (int i_p = 0; i_p < osi_p_parton[0].size(); i_p++){logger << LOG_DEBUG_VERBOSE << "osi_p_parton[0][" << i_p << "] = " << osi_p_parton[0][i_p] << endl;}
      psi.determine_dipole_phasespace_RA(osi_RA_dipole);

      for (int i_a = 1; i_a < osi_n_ps; i_a++){
	for (int i_p = 0; i_p < osi_p_parton[i_a].size(); i_p++){
	  osi_p_parton[i_a][i_p] = psi_xbp_all[i_a][intpow(2, i_p - 1)];
	}
      }

    perform_event_selection(oset, oset.xmunich->generic);
    for (int i_a = 0; i_a < osi_n_ps; i_a++){
      if (oset.cut_ps[i_a] == -1){
	out_comparison << "Phase-space " << i_a << " is cut." << endl;
	for (int sd = 1; sd < oset.max_dyn_ren + 1; sd++){
	  for (int ss = 0; ss < oset.n_scale_dyn_ren[sd]; ss++){
	    oset.value_scale_ren[i_a][sd][ss] = 100.;
	  }
	}
	for (int sd = 1; sd < oset.max_dyn_fact + 1; sd++){
	  for (int ss = 0; ss < oset.n_scale_dyn_fact[sd]; ss++){
	    oset.value_scale_fact[i_a][sd][ss] = 100.;
	  }
	}
	
	oset.cut_ps[i_a] = 0;
	// scales need to be set somehow !!!
      }
      else {
        oset.xmunich->generic.calculate_dynamic_scale_RA(i_a, oset);
        oset.xmunich->generic.calculate_dynamic_scale_TSV(i_a, oset);
	oset.determine_scale_RA(i_a);
      }
    }

    calculate_ME2_RA_QCD(oset);
      
    oset.output_testpoint_RA(out_comparison);
      
    out_comparison << "Corresponding phase-space points:" << endl << endl;
    for (int i_a = 0; i_a < osi_n_ps; i_a++){
      out_comparison << setw(12) << osi_RA_dipole[i_a].name() << endl << endl;
      vector<vector<int> > xdx_pa = osi_RA_dipole[i_a].dx_pa();
      oset.output_momenta_phasespace(out_comparison, i_a);
    }
      
    /*
      // scale output:

      out_comparison << "osi_particle_event[osi_access_object[bjet]][0].size() = " << osi_particle_event[osi_access_object["bjet"]][0].size() << endl;
      out_comparison << "osi_n_object[osi_access_object[bjet]][0] = " << osi_n_object[osi_access_object["bjet"]][0] << endl;

      out_comparison << "osi_particle_event[osi_access_object[ljet]][0].size() = " << osi_particle_event[osi_access_object["ljet"]][0].size() << endl;
      out_comparison << "osi_n_object[osi_access_object[ljet]][0] = " << osi_n_object[osi_access_object["ljet"]][0] << endl;

      out_comparison << "osi_particle_event[osi_access_object[jet]][0].size() = " << osi_particle_event[osi_access_object["jet"]][0].size() << endl;
      out_comparison << "osi_n_object[osi_access_object[jet]][0] = " << osi_n_object[osi_access_object["jet"]][0] << endl;

      for (int i_p = 0; i_p < osi_particle_event[osi_access_object["bjet"]][0].size(); i_p++){
	out_comparison << "bjet[" << i_p << "]: pT = " << osi_particle_event[osi_access_object["bjet"]][0][i_p].pT << "   eta = " << osi_particle_event[osi_access_object["bjet"]][0][i_p].eta << "   m = " << osi_particle_event[osi_access_object["bjet"]][0][i_p].m << endl;
      }
    */
    
    for (int i_a = 0; i_a < osi_n_ps; i_a++){
      out_comparison << endl;
      logger << LOG_DEBUG_VERBOSE << "oset.max_dyn_ren = " << oset.max_dyn_ren << endl;
      logger << LOG_DEBUG_VERBOSE << "oset.n_scale_dyn_ren.size() = " << oset.n_scale_dyn_ren.size() << endl;
      for (int sd = 1; sd < oset.max_dyn_ren + 1; sd++){
	logger << LOG_DEBUG_VERBOSE << "oset.n_scale_dyn_ren[" << sd << "] = " << oset.n_scale_dyn_ren[sd] << endl;
	for (int ss = 0; ss < oset.n_scale_dyn_ren[sd]; ss++){
	  out_comparison << "value_scale_ren[" << i_a << "][" << sd << "][" << ss << "] = " << oset.value_scale_ren[i_a][sd][ss] << endl;
	}
      }
      logger << LOG_DEBUG_VERBOSE << "oset.max_dyn_fact = " << oset.max_dyn_fact << endl;
      logger << LOG_DEBUG_VERBOSE << "oset.n_scale_fact_ren.size() = " << oset.n_scale_dyn_fact.size() << endl;
      for (int sd = 1; sd < oset.max_dyn_fact + 1; sd++){
	logger << LOG_DEBUG_VERBOSE << "oset.n_scale_dyn_fact[" << sd << "] = " << oset.n_scale_dyn_fact[sd] << endl;
	for (int ss = 0; ss < oset.n_scale_dyn_fact[sd]; ss++){
	  out_comparison << "value_scale_fact[" << i_a << "][" << sd << "][" << ss << "] = " << oset.value_scale_fact[i_a][sd][ss] << endl;
	}
      }
    }
    
    
    out_comparison.close();
    }
  }
  
  OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




void calculate_ME2_RA_QEW(observable_set & oset){
  static Logger logger("calculate_ME2_QEW_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi_cut_ps[0] > -1){

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
    ol_evaluate_tree(osi_RA_dipole[0].process_id, P, &osi_value_ME2term[0]);
    logger << LOG_DEBUG_POINT << "OpenLoops:     R_ME2 = " << setw(23) << setprecision(15) << osi_value_ME2term[0] << endl;
    delete [] P;

    }
  }
  else {osi_value_ME2term[0] = 0.;}

  for (int i_a = 1; i_a < osi_n_ps; i_a++){
    if (osi_cut_ps[i_a] >= 0){
      if (osi_RA_dipole[i_a].massive() == 0){
	if      (osi_RA_dipole[i_a].type_dipole() == 1){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ij_k(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 2){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ij_a(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 3){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ai_k(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 5){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ai_b(i_a);}
	else {cout << "Should not happen!" << endl;}
      }
      else {
	if      (osi_RA_dipole[i_a].type_dipole() == 1){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ij_k_massive(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 2){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ij_a_massive(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 3){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ai_k(i_a);}
	else if (osi_RA_dipole[i_a].type_dipole() == 5){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ai_b(i_a);}
	else {cout << "Should not happen!" << endl;}
      }
    }
    else {osi_value_ME2term[i_a] = 0.;}
  }
  osi_RA_ME2 = osi_value_ME2term;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void calculate_ME2check_RA_QEW(phasespace_set & psi, observable_set & oset){
  static Logger logger("calculate_ME2check_QEW_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (oset.switch_output_comparison){
    
    oset.xmunich->generic.phasespacepoint_psp(oset);
    
    if (osi_p_parton[0][0].x0() == 0.){
      if (oset.switch_OL){oset.testpoint_from_OL_rambo();}
    }

    if (osi_p_parton[0][0].x0() != 0.){
      ofstream out_comparison;
      out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
      for (int i_p = 0; i_p < osi_p_parton[0].size(); i_p++){psi_xbp_all[0][intpow(2, i_p - 1)] = osi_p_parton[0][i_p];}
      psi.determine_dipole_phasespace_RA(osi_RA_dipole);

      for (int i_a = 1; i_a < osi_n_ps; i_a++){
	for (int i_p = 0; i_p < osi_p_parton[i_a].size(); i_p++){
	  osi_p_parton[i_a][i_p] = psi_xbp_all[i_a][intpow(2, i_p - 1)];
	}
      }
    
      perform_event_selection(oset, oset.xmunich->generic);
      for (int i_a = 0; i_a < osi_n_ps; i_a++){
	if (oset.cut_ps[i_a] == -1){
	  out_comparison << "Phase-space " << i_a << " is cut." << endl;
	  oset.cut_ps[i_a] = 0;
	  // scales need to be set somehow !!!
	}
	else {
	  oset.xmunich->generic.calculate_dynamic_scale_RA(i_a, oset);
	  oset.xmunich->generic.calculate_dynamic_scale_TSV(i_a, oset);
	  oset.determine_scale_RA(i_a);
	}
      }
    
      calculate_ME2_RA_QEW(oset);

      oset.output_testpoint_RA(out_comparison);

      out_comparison << "Corresponding phase-space points:" << endl << endl;
      for (int i_a = 0; i_a < osi_n_ps; i_a++){
	out_comparison << setw(12) << osi_RA_dipole[i_a].name() << endl << endl;
	vector<vector<int> > xdx_pa = osi_RA_dipole[i_a].dx_pa();
	oset.output_momenta_phasespace(out_comparison, i_a);
      }
      
      out_comparison.close();
    }
  }
  
  OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void calculate_ME2_RA_MIX(observable_set & oset){
  static Logger logger("calculate_ME2_MIX_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi_cut_ps[0] > -1){
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
    ol_evaluate_tree(osi_RA_dipole[0].process_id, P, &osi_value_ME2term[0]);
    delete [] P;
  }
  else {osi_value_ME2term[0] = 0.;}

  for (int i_a = 1; i_a < osi_n_ps; i_a++){
    if (osi_cut_ps[i_a] >= 0){
      if      (osi_RA_dipole[i_a].type_correction() == 1){
	if (osi_RA_dipole[i_a].massive() == 0){      
	  if      (osi_RA_dipole[i_a].type_dipole() == 1){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ij_k(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 2){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ij_a(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 3){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ai_k(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 5){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ai_b(i_a);}
	  else {cout << "Should not happen!" << endl;}
	}
	else {
	  if      (osi_RA_dipole[i_a].type_dipole() == 1){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ij_k_massive(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 2){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ij_a_massive(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 3){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ai_k(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 5){osi_value_ME2term[i_a] = oset.calculate_dipole_QCD_A_ai_b(i_a);}
	}
      }
      else if (osi_RA_dipole[i_a].type_correction() == 2){
	if (osi_RA_dipole[i_a].massive() == 0){
	  if      (osi_RA_dipole[i_a].type_dipole() == 1){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ij_k(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 2){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ij_a(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 3){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ai_k(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 5){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ai_b(i_a);}
	  else {cout << "Should not happen!" << endl;}
	}
	else {
	  if      (osi_RA_dipole[i_a].type_dipole() == 1){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ij_k_massive(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 2){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ij_a_massive(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 3){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ai_k(i_a);}
	  else if (osi_RA_dipole[i_a].type_dipole() == 5){osi_value_ME2term[i_a] = oset.calculate_dipole_QEW_A_ai_b(i_a);}
	  else {cout << "Should not happen!" << endl;}
	}
      }
    }
    else {osi_value_ME2term[i_a] = 0.;}
  }
  osi_RA_ME2 = osi_value_ME2term;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void calculate_ME2check_RA_MIX(phasespace_set & psi, observable_set & oset){
  static Logger logger("calculate_ME2check_MIX_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (oset.switch_output_comparison){
    
    oset.xmunich->generic.phasespacepoint_psp(oset);
    
    if (osi_p_parton[0][0].x0() == 0.){
      if (oset.switch_OL){oset.testpoint_from_OL_rambo();}
    }

    if (osi_p_parton[0][0].x0() != 0.){
      ofstream out_comparison;
      out_comparison.open(osi_filename_comparison.c_str(), ofstream::out | ofstream::app);  
      for (int i_p = 0; i_p < osi_p_parton[0].size(); i_p++){psi_xbp_all[0][intpow(2, i_p - 1)] = osi_p_parton[0][i_p];}
      psi.determine_dipole_phasespace_RA(osi_RA_dipole);
      
      for (int i_a = 1; i_a < osi_n_ps; i_a++){
	for (int i_p = 0; i_p < osi_p_parton[i_a].size(); i_p++){
	  osi_p_parton[i_a][i_p] = psi_xbp_all[i_a][intpow(2, i_p - 1)];
	}
      }
      
      perform_event_selection(oset, oset.xmunich->generic);
      for (int i_a = 0; i_a < osi_n_ps; i_a++){
	if (oset.cut_ps[i_a] == -1){
	  out_comparison << "Phase-space " << i_a << " is cut." << endl;
	  oset.cut_ps[i_a] = 0;
	  // scales need to be set somehow !!!
	}
	else {
	  oset.xmunich->generic.calculate_dynamic_scale_RA(i_a, oset);
	  oset.xmunich->generic.calculate_dynamic_scale_TSV(i_a, oset);
	  oset.determine_scale_RA(i_a);
	}
      }
    
      calculate_ME2_RA_MIX(oset);
      
      oset.output_testpoint_RA(out_comparison);
      
      out_comparison << "Corresponding phase-space points:" << endl << endl;
      for (int i_a = 0; i_a < osi_n_ps; i_a++){
	out_comparison << setw(12) << osi_RA_dipole[i_a].name() << endl << endl;
	vector<vector<int> > xdx_pa = osi_RA_dipole[i_a].dx_pa();
	oset.output_momenta_phasespace(out_comparison, i_a);
      }
      
      out_comparison.close();
    }
  }

  OLP_PrintParameter(stch("log/olparameters." + osi_name_process + ".txt"));

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

