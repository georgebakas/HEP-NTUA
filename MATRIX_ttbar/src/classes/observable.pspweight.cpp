#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

// Relevant for TSV variation:

void observable_set::determine_psp_weight_TSV(phasespace_set & psi){
  Logger logger("observable_set::determine_psp_weight_TSV");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  if (!switch_TSV){return;}

  logger << LOG_DEBUG_VERBOSE << "type_contribution = " << type_contribution << endl;

  static int x_q = 0;
  for (int i_s = 0; i_s < n_set_TSV; i_s++){
    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	for (int i_h = 0; i_h < max_n_integrand_TSV[i_s]; i_h++){
	  if (type_contribution == "CT" || 
	      type_contribution == "CJ" || 
	      type_contribution == "L2CT" || 
	      type_contribution == "L2CJ"){ // CT, L2CT (intrinsic qTcut dependence -> not affected by cut_ps)
	    ps_integrand_TSV[i_s][i_r][i_f][i_h][0] = psi_ps_factor * rescaling_factor_alpha_e 
	      * value_fact_integrand_qTcut_TSV[x_q][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h]
	      * *(pointer_relative_factor_alpha_S)[0][i_s][i_r];
	    for (int i_q = 0; i_q < n_qTcut; i_q++){
	      ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][i_h][0] = psi_ps_factor * rescaling_factor_alpha_e 
		* value_fact_integrand_qTcut_TSV[i_q][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h]
		* *(pointer_relative_factor_alpha_S)[0][i_s][i_r];
	    }
	    logger << LOG_DEBUG_VERBOSE << "value_fact_integrand_qTcut_TSV[" << x_q << "][" << dynamic_scale_fact_TSV[i_s] << "][" << no_value_fact_TSV[i_s][i_f] << "][" << i_h << "] = " << value_fact_integrand_qTcut_TSV[x_q][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h] << endl;
	  }
	  // split from VT/VT2/CT2 !!! -> simplifications possible...
	  else if (type_contribution == "CT2" || 
		   type_contribution == "CJ2"){ // CT2 (intrinsic qTcut dependence -> not affected by cut_ps)
	    ps_integrand_TSV[i_s][i_r][i_f][i_h][0] = psi_ps_factor * rescaling_factor_alpha_e 
	      * value_integrand_qTcut_TSV[x_q][dynamic_scale_ren_TSV[i_s]][dynamic_scale_fact_TSV[i_s]][no_value_ren_TSV[i_s][i_r]][no_value_fact_TSV[i_s][i_f]][i_h]
	      * *(pointer_relative_factor_alpha_S)[0][i_s][i_r];
	    // output_n_qTcut  or  n_qTcut ??? (difference between CT2 and the others...)
	    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	      ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][i_h][0] = psi_ps_factor * rescaling_factor_alpha_e 
		* value_integrand_qTcut_TSV[i_q][dynamic_scale_ren_TSV[i_s]][dynamic_scale_fact_TSV[i_s]][no_value_ren_TSV[i_s][i_r]][no_value_fact_TSV[i_s][i_f]][i_h]
		* *(pointer_relative_factor_alpha_S)[0][i_s][i_r];
	    }
	    logger << LOG_DEBUG_VERBOSE << "value_integrand_qTcut_TSV[" << x_q << "][" << dynamic_scale_ren_TSV[i_s] << "][" << dynamic_scale_fact_TSV[i_s] << "][" << no_value_ren_TSV[i_s][i_r] << "][" << no_value_fact_TSV[i_s][i_f] << "][" << i_h << "] = " << value_integrand_qTcut_TSV[x_q][dynamic_scale_ren_TSV[i_s]][dynamic_scale_fact_TSV[i_s]][no_value_ren_TSV[i_s][i_r]][no_value_fact_TSV[i_s][i_f]][i_h] << endl;
	  }
	  // split from VT/VT2/CT2 !!! -> simplifications possible...
	  else if (type_contribution == "VT" || 
		   type_contribution == "VT2" || 
		   type_contribution == "VJ" || 
		   type_contribution == "VJ2" || 
		   type_contribution == "L2VT" || 
		   type_contribution == "L2VJ"){ // VT, VT2, L2VT (no qTcut dependence -> affected only by cut_ps)
	    ps_integrand_TSV[i_s][i_r][i_f][i_h][0] = psi_ps_factor * rescaling_factor_alpha_e 
	      * value_integrand_qTcut_TSV[x_q][dynamic_scale_ren_TSV[i_s]][dynamic_scale_fact_TSV[i_s]][no_value_ren_TSV[i_s][i_r]][no_value_fact_TSV[i_s][i_f]][i_h]
	      * *(pointer_relative_factor_alpha_S)[0][i_s][i_r];
	    // output_n_qTcut  or  n_qTcut ??? (difference between CT2 and the others...)
	    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	      ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][i_h][0] = psi_ps_factor * rescaling_factor_alpha_e 
		* value_integrand_qTcut_TSV[i_q][dynamic_scale_ren_TSV[i_s]][dynamic_scale_fact_TSV[i_s]][no_value_ren_TSV[i_s][i_r]][no_value_fact_TSV[i_s][i_f]][i_h]
		* *(pointer_relative_factor_alpha_S)[0][i_s][i_r];
	    }
	    logger << LOG_DEBUG_VERBOSE << "value_integrand_qTcut_TSV[" << x_q << "][" << dynamic_scale_ren_TSV[i_s] << "][" << dynamic_scale_fact_TSV[i_s] << "][" << no_value_ren_TSV[i_s][i_r] << "][" << no_value_fact_TSV[i_s][i_f] << "][" << i_h << "] = " << value_integrand_qTcut_TSV[x_q][dynamic_scale_ren_TSV[i_s]][dynamic_scale_fact_TSV[i_s]][no_value_ren_TSV[i_s][i_r]][no_value_fact_TSV[i_s][i_f]][i_h] << endl;
	  }
	  else if (n_pz == 1){ // Born, VA, RA, loop, RVA, RRA (no qTcut dependence -> affected only by cut_ps)
	    for (int i_a = 0; i_a < n_ps; i_a++){
	      ps_integrand_TSV[i_s][i_r][i_f][i_h][i_a] = psi_ps_factor * rescaling_factor_alpha_e 
		* *(pointer_ME2term)[i_a][0][i_s][i_r][i_f]
		* *(pointer_pdf_factor)[i_a][0][i_s][i_f][i_h]
		* *(pointer_relative_factor_alpha_S)[i_a][i_s][i_r];
	      // logger << LOG_DEBUG_VERBOSE << "*(pointer_ME2term)[i_a][0][i_s][i_r][i_f] = " << *(pointer_ME2term)[i_a][0][i_s][i_r][i_f] << endl;
	      // logger << LOG_DEBUG_VERBOSE << "*(pointer_pdf_factor)[i_a][0][i_s][i_f][i_h] = " << *(pointer_pdf_factor)[i_a][0][i_s][i_f][i_h] << endl;
	      // logger << LOG_DEBUG_VERBOSE << "*(pointer_relative_factor_alpha_S)[i_a][i_s][i_r] = " << *(pointer_relative_factor_alpha_S)[i_a][i_s][i_r] << endl;
	    }
	  }
	  else { // CA, RCA (no qTcut dependence -> affected only by cut_ps)
	    vector<double> temp_c(n_pc, 0.);
	    for (int i_c = 0; i_c < n_pc; i_c++){
	      for (int i_z = 0; i_z < n_pz; i_z++){
		if (i_z == 0 && (*CA_collinear)[i_c][0].in_collinear()[i_z] == 1){
		  temp_c[i_c] += *(pointer_pdf_factor)[i_c][i_z][i_s][i_f][i_h] 
		    * (  *(pointer_ME2term)[i_c][1][i_s][i_r][i_f] / psi_g_z_coll[(*CA_collinear)[i_c][0].no_emitter()] 
			 + *(pointer_ME2term)[i_c][2][i_s][i_r][i_f]);
		}
		else if ((*CA_collinear)[i_c][0].in_collinear()[i_z] == 1){ // always i_z == dipole_phasespace[i_c][4] !!!
		  temp_c[i_c] += *(pointer_pdf_factor)[i_c][i_z][i_s][i_f][i_h] 
		    * *(pointer_ME2term)[i_c][0][i_s][i_r][i_f] / psi_g_z_coll[i_z];
		}
	      }
	    }
	    ps_integrand_TSV[i_s][i_r][i_f][i_h][0] = psi_ps_factor * rescaling_factor_alpha_e * accumulate(temp_c.begin(), temp_c.end(), 0.) * *(pointer_relative_factor_alpha_S)[0][i_s][i_r];
	  }

	  logger << LOG_DEBUG_VERBOSE << "ps_integrand_TSV[" << i_s << "][" << i_r << "][" << i_f << "][" << i_h << "][0] = " << ps_integrand_TSV[i_s][i_r][i_f][i_h][0] << endl;
	  logger << LOG_DEBUG_VERBOSE << "ps_integrand_qTcut_TSV[" << x_q << "][" << i_s << "][" << i_r << "][" << i_f << "][" << i_h << "][0] = " << ps_integrand_qTcut_TSV[x_q][i_s][i_r][i_f][i_h][0] << endl;
	}

	if (n_ps == 1 && n_pc == 1){integrand_TSV[i_s][i_r][i_f] = ps_integrand_TSV[i_s][i_r][i_f][0][0];}
	// What happens here for CT/CT2 ? Result ignored ???
	else {integrand_TSV[i_s][i_r][i_f] = accumulate(ps_integrand_TSV[i_s][i_r][i_f][0].begin(), ps_integrand_TSV[i_s][i_r][i_f][0].end(), 0.);}

	logger << LOG_DEBUG_VERBOSE << "integrand_TSV[" << i_s << "][" << i_r << "][" << i_f << "] = " << integrand_TSV[i_s][i_r][i_f] << endl;

	if (type_contribution == "CT" || 
	    type_contribution == "CJ" || 
	    type_contribution == "CT2" || 
	    type_contribution == "CJ2" || 
	    type_contribution == "L2CT" || 
	    type_contribution == "L2CJ"){
	  for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	    integrand_qTcut_TSV[i_q][i_s][i_r][i_f] = ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][0][0];
	  }
	  logger << LOG_DEBUG_VERBOSE << "integrand_qTcut_TSV[" << x_q << "][" << i_s << "][" << i_r << "][" << i_f << "] = " << integrand_qTcut_TSV[x_q][i_s][i_r][i_f] << endl;
	}
      }
    }
  }



  // scale-difference extension start
  // (filled with differences)
  // to be filled / extended in size:
  // ps_integrand_TSV
  // integrand_TSV
  // ps_integrand_qTcut_TSV (VT, CT, VT2, CT2, L2CT, L2VT)
  // integrand_qTcut_TSV (CT, CT2, L2CT)

  for (int i_s = 0; i_s < n_diff_set_TSV; i_s++){
    int i_es = n_set_TSV + i_s;

    for (int i_r = 0; i_r < n_scale_ren_TSV[i_es]; i_r++){
      for (int i_f = 0; i_f < n_scale_fact_TSV[i_es]; i_f++){
	for (int i_h = 0; i_h < max_n_integrand_TSV[i_es]; i_h++){
	  logger << LOG_DEBUG_VERBOSE << "i_h = " << i_h << endl;

	  ps_integrand_TSV[i_es][i_r][i_f][i_h][0] = ps_integrand_TSV[no_diff_set_plus_TSV[i_s]][i_r][i_f][i_h][0] - ps_integrand_TSV[no_diff_set_minus_TSV[i_s]][i_r][i_f][i_h][0];
	  logger << LOG_DEBUG_VERBOSE << "ps_integrand_TSV[" << i_es << "][" << i_r << "][" << i_f << "][" << i_h << "][0] = " << ps_integrand_TSV[i_es][i_r][i_f][i_h][0] << endl;

	  integrand_TSV[i_es][i_r][i_f] = integrand_TSV[no_diff_set_plus_TSV[i_s]][i_r][i_f] - integrand_TSV[no_diff_set_minus_TSV[i_s]][i_r][i_f];
	  logger << LOG_DEBUG_VERBOSE << "integrand_TSV[" << i_es << "][" << i_r << "][" << i_f << "] = " << integrand_TSV[i_es][i_r][i_f] << endl;

	  if (type_contribution == "VT" || 
	      type_contribution == "CT" || 
	      type_contribution == "VT2" || 
	      type_contribution == "CT2" || 
	      type_contribution == "VJ" || 
	      type_contribution == "CJ" || 
	      type_contribution == "VJ2" || 
	      type_contribution == "CJ2" || 
	      type_contribution == "L2CT"||
	      type_contribution == "L2VT"){
	    for (int i_q = 0; i_q < n_qTcut; i_q++){
	      ps_integrand_qTcut_TSV[i_q][i_es][i_r][i_f][i_h][0] = ps_integrand_qTcut_TSV[i_q][no_diff_set_plus_TSV[i_s]][i_r][i_f][i_h][0] - ps_integrand_qTcut_TSV[i_q][no_diff_set_minus_TSV[i_s]][i_r][i_f][i_h][0];
	      logger << LOG_DEBUG_VERBOSE << "ps_integrand_qTcut_TSV[" << x_q << "][" << i_es << "][" << i_r << "][" << i_f << "][" << i_h << "][0] = " << ps_integrand_qTcut_TSV[x_q][i_es][i_r][i_f][i_h][0] << endl;
	    }
	  }

	  if (type_contribution == "CT" || 
	      type_contribution == "CT2" || 
	      type_contribution == "CJ" || 
	      type_contribution == "CJ2" || 
	      type_contribution == "L2CT" || 
	      type_contribution == "L2CJ"){
	    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	      integrand_qTcut_TSV[i_q][i_es][i_r][i_f] = integrand_qTcut_TSV[i_q][no_diff_set_plus_TSV[i_s]][i_r][i_f] - integrand_qTcut_TSV[i_q][no_diff_set_minus_TSV[i_s]][i_r][i_f];
	      logger << LOG_DEBUG_VERBOSE << "integrand_qTcut_TSV[" << x_q << "][" << i_es << "][" << i_r << "][" << i_f << "] = " << integrand_qTcut_TSV[x_q][i_es][i_r][i_f] << endl;
	    }
	  }

	}
      }
    }
  }

  // to be filled / extended in size:
  // sum_weight_TSV
  // sum_weight2_TSV
  // ps_moment_TSV
  // moment_TSV
  // sum_weight_qTcut_TSV
  // sum_weight2_qTcut_TSV
  // sum_moment_qTcut_TSV
  // sum_moment2_qTcut_TSV
  // 
  // end



  logger << LOG_DEBUG_VERBOSE << "active_qTcut = " << active_qTcut << endl;
  logger << LOG_DEBUG_VERBOSE << "type_contribution = " << type_contribution << endl;

  if (!active_qTcut){ // Born, VA, CA, RA, VT, VT2, loop
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	  sum_weight_TSV[i_s][i_r][i_f] += integrand_TSV[i_s][i_r][i_f];
	  sum_weight2_TSV[i_s][i_r][i_f] += pow(integrand_TSV[i_s][i_r][i_f], 2);
	}
      }
      if (switch_moment_TSV[i_s] > 0){
	for (int i_m = 0; i_m < n_moments; i_m++){
	  for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	      for (int i_h = 0; i_h < max_n_integrand_TSV[i_s]; i_h++){
		for (int i_a = 0; i_a < n_ps; i_a++){
		  ps_moment_TSV[i_s][i_m][i_r][i_f][i_h][i_a] = ps_integrand_TSV[i_s][i_r][i_f][i_h][i_a] * directed_moment[i_m][i_h][i_a];
		}
	      }
	      if (n_ps == 1){
		// if the moment is 'directed' in the sense of being different when the psp is rotated, the 2nd version must be used !!!
		if (moment_symm[i_m] == 0){moment_TSV[i_s][i_m][i_r][i_f] = ps_moment_TSV[i_s][i_m][i_r][i_f][0][0];}
		else if (moment_symm[i_m] == 1){moment_TSV[i_s][i_m][i_r][i_f] = ps_moment_TSV[i_s][i_m][i_r][i_f][1][0] + ps_moment_TSV[i_s][i_m][i_r][i_f][2][0];}
	      }
	      else {
		vector<double> temp_moment(max_n_integrand_TSV[i_s]);
		if (moment_symm[i_m] == 0){
		  moment_TSV[i_s][i_m][i_r][i_f] = accumulate(ps_moment_TSV[i_s][i_m][i_r][i_f][0].begin(), ps_moment_TSV[i_s][i_m][i_r][i_f][0].end(), 0.);
		}
		else if (moment_symm[i_m] == 1){
		  for (int i_h = 0; i_h < max_n_integrand_TSV[i_s]; i_h++){
		    temp_moment[i_h] = accumulate(ps_moment_TSV[i_s][i_m][i_r][i_f][i_h].begin(), ps_moment_TSV[i_s][i_m][i_r][i_f][i_h].end(), 0.);
		  }
		  moment_TSV[i_s][i_m][i_r][i_f] = accumulate(temp_moment.begin(), temp_moment.end(), 0.);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  else if (active_qTcut && (type_contribution == "CT" || 
			    type_contribution == "CT2" || 
			    type_contribution == "CJ" || 
			    type_contribution == "CJ2" || 
			    type_contribution == "L2CT" || 
			    type_contribution == "L2CJ")){
    for (int i_c = 0; i_c < n_qTcut; i_c++){
      for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	    sum_weight_qTcut_TSV[i_s][i_c][i_r][i_f] += integrand_qTcut_TSV[i_c][i_s][i_r][i_f];
	    sum_weight2_qTcut_TSV[i_s][i_c][i_r][i_f] += pow(integrand_qTcut_TSV[i_c][i_s][i_r][i_f], 2);;
	    logger << LOG_DEBUG_VERBOSE << "sum_weight_qTcut_TSV[" << i_s << "][" << i_c << "][" << i_r << "][" << i_f << "] = " << sum_weight_qTcut_TSV[i_s][i_c][i_r][i_f] << endl;
	  }
	}
      }
    }
  }
  else { // active_qTcut, but entering only via cut_ps[i_a]: RT, RVA, RCA, RRA
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      if (n_ps == 1){ // only one phasespace: RT, RVA, RCA
	vector<vector<double> > temp_weight2(n_scale_ren_TSV[i_s], vector<double> (n_scale_ren_TSV[i_s], 0.));
	for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	    temp_weight2[i_r][i_f] += pow(integrand_TSV[i_s][i_r][i_f], 2);
	  }
	}
	for (int i_c = 0; i_c < n_qTcut; i_c++){
	  if (i_c <= cut_ps[0]){
	    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		sum_weight_qTcut_TSV[i_s][i_c][i_r][i_f] += integrand_TSV[i_s][i_r][i_f];
		sum_weight2_qTcut_TSV[i_s][i_c][i_r][i_f] += temp_weight2[i_r][i_f];
	      }
	    }
	  }
	}
	if (switch_moment_TSV[i_s] > 0){
	  for (int i_m = 0; i_m < n_moments; i_m++){
	    vector<vector<double> > temp_moment(n_scale_ren_TSV[i_s], vector<double> (n_scale_ren_TSV[i_s], 0.));
	    vector<vector<double> > temp_moment2(n_scale_ren_TSV[i_s], vector<double> (n_scale_ren_TSV[i_s], 0.));
	    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		vector<double> temp_moment_directed(3, 0.);
		for (int i_h = 0; i_h < max_n_integrand_TSV[i_s]; i_h++){
		  temp_moment_directed[i_h] = ps_integrand_TSV[i_s][i_r][i_f][i_h][0] * directed_moment[i_m][i_h][0];
		}
		if (moment_symm[i_m] == 0){temp_moment[i_r][i_f] = temp_moment_directed[0];}
		else if (moment_symm[i_m] == 1){temp_moment[i_r][i_f] = temp_moment_directed[1] + temp_moment_directed[2];}
		temp_moment2[i_r][i_f] = pow(temp_moment[i_r][i_f], 2);
	      }
	    }
	    for (int i_c = 0; i_c < n_qTcut; i_c++){
	      if (i_c <= cut_ps[0]){
		for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		    sum_moment_qTcut_TSV[i_s][i_m][i_c][i_r][i_f] += temp_moment[i_r][i_f];
		    sum_moment2_qTcut_TSV[i_s][i_m][i_c][i_r][i_f] += temp_moment2[i_r][i_f];
		  }
		}
	      }
	    }
	  }
	}
      }
      else { // more than one phasespace: RRA
	vector<vector<double> > temp_weight(n_scale_ren_TSV[i_s], vector<double> (n_scale_ren_TSV[i_s], 0.));
	vector<vector<double> > temp_weight2(n_scale_ren_TSV[i_s], vector<double> (n_scale_ren_TSV[i_s], 0.));
	for (int i_c = 0; i_c < n_qTcut; i_c++){
	  if (i_c == 0 || change_cut[i_c] == 1){
	    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		temp_weight[i_r][i_f] = 0.;
		for (int i_a = 0; i_a < n_ps; i_a++){
		  if (i_c <= cut_ps[i_a]){
		    temp_weight[i_r][i_f] += ps_integrand_TSV[i_s][i_r][i_f][0][i_a];
		  }
		}
		temp_weight2[i_r][i_f] = pow(temp_weight[i_r][i_f], 2);
	      }
	    }
	    //	    change_cut[i_c] = 0; // would kill any second independent scale-variation calculation !!! shifted to later...
	  }
	  
	  for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	      sum_weight_qTcut_TSV[i_s][i_c][i_r][i_f] += temp_weight[i_r][i_f];
	      sum_weight2_qTcut_TSV[i_s][i_c][i_r][i_f] += temp_weight2[i_r][i_f];
	    }
	  }
	}

	if (switch_moment_TSV[i_s] > 0){
	  for (int i_m = 0; i_m < n_moments; i_m++){
	    vector<vector<double> > temp_moment(n_scale_ren_TSV[i_s], vector<double> (n_scale_ren_TSV[i_s], 0.));
	    vector<vector<double> > temp_moment2(n_scale_ren_TSV[i_s], vector<double> (n_scale_ren_TSV[i_s], 0.));
	    for (int i_c = 0; i_c < n_qTcut; i_c++){
	      if (i_c == 0 || change_cut[i_c] == 1){
		for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		    vector<double> temp_moment_directed(3, 0.);
		    for (int i_a = 0; i_a < n_ps; i_a++){
		      if (i_c <= cut_ps[i_a]){
			for (int i_h = 0; i_h < max_n_integrand_TSV[i_s]; i_h++){
			  temp_moment_directed[i_h] += ps_integrand_TSV[i_s][i_r][i_f][i_h][i_a] * directed_moment[i_m][i_h][i_a];
			}
		      }
		    }
		    if (moment_symm[i_m] == 0){temp_moment[i_r][i_f] = temp_moment_directed[0];}
		    else if (moment_symm[i_m] == 1){temp_moment[i_r][i_f] = temp_moment_directed[1] + temp_moment_directed[2];}
		    temp_moment2[i_r][i_f] = pow(temp_moment[i_r][i_f], 2);
		  }
		}
	      }
	      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		  sum_moment_qTcut_TSV[i_s][i_m][i_c][i_r][i_f] += temp_moment[i_r][i_f];
		  sum_moment2_qTcut_TSV[i_s][i_m][i_c][i_r][i_f] += temp_moment2[i_r][i_f];
		}
	      }
	    }
	  }
	}
      }
    }

    // reset change_cut for the next phase-space point:
    for (int i_c = 0; i_c < n_qTcut; i_c++){if (change_cut[i_c] == 1){change_cut[i_c] = 0;}}
  }

  //  if (active_qTcut){logger << LOG_DEBUG << setw(10) << left << psi_i_acc << "   " << setw(10) << left << psi_i_gen << "   " << sum_weight_qTcut_TSV[0][0][0][0] << endl;}
  //  else {logger << LOG_DEBUG << setw(10) << left << psi_i_acc << "   " << setw(10) << left << psi_i_gen << "   " << sum_weight_TSV[0][0][0] << endl;}

  logger << LOG_DEBUG_VERBOSE << "no_reference_TSV            = " << no_reference_TSV << endl;
  logger << LOG_DEBUG_VERBOSE << "no_scale_ren_reference_TSV  = " << no_scale_ren_reference_TSV << endl;
  logger << LOG_DEBUG_VERBOSE << "no_scale_fact_reference_TSV = " << no_scale_fact_reference_TSV << endl;
  logger << LOG_DEBUG_VERBOSE << "no_qTcut_reference_TSV      = " << no_qTcut_reference_TSV << endl;

  if (no_reference_TSV > -1){
    logger << LOG_DEBUG_VERBOSE << "before: integrand = " << integrand << endl;
    /*
    // to be refined...
    if (n_ps == 1 && n_pc == 1){ // born, VA, CA
      integrand = integrand_TSV[no_reference_TSV][no_scale_ren_reference_TSV][no_scale_fact_reference_TSV];
    }
    else { // RA - actually identical here
      integrand = integrand_TSV[no_reference_TSV][no_scale_ren_reference_TSV][no_scale_fact_reference_TSV];
    }
    */

    if (type_contribution == "CT" || 
	type_contribution == "CT2" || 
	type_contribution == "CJ" || 
	type_contribution == "CJ2" || 
	type_contribution == "L2CT" || 
	type_contribution == "L2CJ"){
      integrand = integrand_qTcut_TSV[no_qTcut_reference_TSV][no_reference_TSV][no_scale_ren_reference_TSV][no_scale_fact_reference_TSV];
    }
    else {
      integrand = integrand_TSV[no_reference_TSV][no_scale_ren_reference_TSV][no_scale_fact_reference_TSV];
    }
    logger << LOG_DEBUG_VERBOSE << "after: integrand  = " << integrand << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}





// Relevant for CV variation and reference result:

void observable_set::determine_scale(){
  Logger logger("observable_set::determine_scale");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (dynamic_scale != 0){
    var_mu_fact = value_mu_fact[dynamic_scale][map_value_scale_fact];
    var_mu_ren = value_mu_ren[dynamic_scale][map_value_scale_ren];
    var_alpha_S_reference = value_alpha_S[dynamic_scale][map_value_scale_ren];
    var_rel_alpha_S = value_factor_alpha_S[dynamic_scale][map_value_scale_ren];
    logger << LOG_DEBUG_VERBOSE << "var_mu_fact = " << var_mu_fact << endl;
    logger << LOG_DEBUG_VERBOSE << "var_mu_ren = " << var_mu_ren << endl;
    logger << LOG_DEBUG_VERBOSE << "var_rel_alpha_S = " << var_rel_alpha_S << endl;
  }
  
  if (switch_CV != 0){
    if (dynamic_scale_CV != 0){
      for (int s = 0; s < n_scales_CV; s++){
	var_mu_fact_CV[s] = value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]];
	var_mu_ren_CV[s] = value_mu_ren[dynamic_scale_CV][map_value_scale_ren_CV[s]];
	var_alpha_S_CV[s] = value_alpha_S[dynamic_scale_CV][map_value_scale_ren_CV[s]];
	var_rel_alpha_S_CV[s] = value_factor_alpha_S[dynamic_scale_CV][map_value_scale_ren_CV[s]];
	logger << LOG_DEBUG_VERBOSE << "value_mu_fact[dynamic_scale_CV = " << dynamic_scale_CV << "][map_value_scale_fact_CV[s = " << s << "] = " << map_value_scale_fact_CV[s] << "] = " << setw(15) << value_mu_fact[dynamic_scale_CV][map_value_scale_fact_CV[s]] << endl;
	logger << LOG_DEBUG_VERBOSE << "var_mu_fact[" << s << "] = " << setw(8) << var_mu_fact_CV[s] << "   var_mu_ren[" << s << "] = " << setw(8) << var_mu_ren_CV[s] << "   var_rel_alpha_S[" << s << "] = " << setw(15) << var_rel_alpha_S_CV[s] << endl;
      }
    }
  }
  
  // only LOG_DEBUG_VERBOSE output:
  if (switch_TSV){
    for (int i_v = 0; i_v < max_dyn_ren + 1; i_v++){
      logger << LOG_DEBUG_VERBOSE << left << "DS: " << i_v << right << setw(25) << "rel_scale_ren" << setw(24) << "rel_scale_ren" << "²" << setw(25) << "rel_factor_alpha_S" << setw(25) << "scale_ren" << setw(24) << "scale_ren" << "²" << endl;
      for (int i_m = 0; i_m < value_relative_scale_ren[i_v].size(); i_m++){
	logger << LOG_DEBUG_VERBOSE << setw(5) << "" << setw(25) << setprecision(15) << value_relative_scale_ren[i_v][i_m] << setw(25) << setprecision(15) << value_relative_scale2_ren[i_v][i_m] << setw(25) << setprecision(15) << value_relative_factor_alpha_S[0][i_v][i_m] << setw(25) << setprecision(15) << value_scale_ren[0][i_v][i_m] << setw(25) << setprecision(15) << value_scale2_ren[0][i_v][i_m] << endl;
    }
      logger.newLine(LOG_DEBUG_VERBOSE);
    }
    logger.newLine(LOG_DEBUG_VERBOSE);
    for (int i_v = 0; i_v < max_dyn_fact + 1; i_v++){
      logger << LOG_DEBUG_VERBOSE << "DS: " << i_v << setw(25) << "rel_scale_fact" << setw(24) << "rel_scale_fact" << "²" << setw(23) << "log(rel_scale_fact" << "²)" << setw(25) << "scale_fact" << setw(24) << "scale_fact" << "²" << endl;
      for (int i_m = 0; i_m < value_relative_scale_fact[i_v].size(); i_m++){
	logger << LOG_DEBUG_VERBOSE << setw(5) << "" << setw(25) << setprecision(15) << value_relative_scale_fact[i_v][i_m] << setw(25) << setprecision(15) << value_relative_scale2_fact[i_v][i_m] << setw(25) << setprecision(15) << value_relative_logscale2_fact[i_v][i_m] << setw(25) << setprecision(15) << value_scale_fact[0][i_v][i_m] << setw(25) << setprecision(15) << value_scale2_fact[0][i_v][i_m] << endl;
      }
      logger.newLine(LOG_DEBUG_VERBOSE);
    }
    logger.newLine(LOG_DEBUG_VERBOSE);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_scale_CA(){
  Logger logger("observable_set::determine_scale_CA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int sd = 1; sd < value_mu_fact_rel.size(); sd++){
    if (value_mu_fact_rel[sd].size() != 0){
      for (int ss = 0; ss < value_mu_fact_rel[sd].size(); ss++){
	CA_value_log_mu2_fact[sd][ss] = 2 * log(value_mu_fact[sd][ss]);
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_scale_RA(int i_a){
  Logger logger("observable_set::determine_scale_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (dynamic_scale != 0){
    RA_mu_fact[i_a] = RA_value_mu_fact[i_a][dynamic_scale][map_value_scale_fact];
    RA_mu_ren[i_a] = RA_value_mu_ren[i_a][dynamic_scale][map_value_scale_ren];
    RA_alpha_S_reference[i_a]  = RA_value_alpha_S[i_a][dynamic_scale][map_value_scale_ren];
    RA_rel_alpha_S[i_a]  = RA_value_factor_alpha_S[i_a][dynamic_scale][map_value_scale_ren];
  }

  if (switch_CV != 0){
    if (dynamic_scale_CV != 0){
      for (int s = 0; s < n_scales_CV; s++){
	RA_mu_fact_CV[i_a][s] = RA_value_mu_fact[i_a][dynamic_scale_CV][map_value_scale_fact_CV[s]];
	RA_mu_ren_CV[i_a][s] = RA_value_mu_ren[i_a][dynamic_scale_CV][map_value_scale_ren_CV[s]];
	RA_alpha_S_CV[i_a][s] = RA_value_alpha_S[i_a][dynamic_scale_CV][map_value_scale_ren_CV[s]];
	RA_rel_alpha_S_CV[i_a][s] = RA_value_factor_alpha_S[i_a][dynamic_scale_CV][map_value_scale_ren_CV[s]];
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_integrand(phasespace_set & psi){
  Logger logger("observable_set::determine_integrand");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  logger << LOG_DEBUG_VERBOSE << "ME2 = " << ME2 << endl;
  integrand = psi_ps_factor * rescaling_factor_alpha_e * pdf_factor[0] * var_rel_alpha_S * ME2;

  logger << LOG_DEBUG_VERBOSE << "integrand = " << integrand << "   pdf_factor[0] = " << pdf_factor[0] << "   var_rel_alpha_S = " << var_rel_alpha_S << "   psi_ps_factor = " << psi_ps_factor << "   rescaling_factor_alpha_e = " << rescaling_factor_alpha_e << "   ME2 = " << ME2 << endl;

  logger << LOG_DEBUG_VERBOSE << setw(20) << "integrand" << " = " << setprecision(20) << setw(28) << integrand << "   " << double2hexastr(integrand) << endl;
  logger << LOG_DEBUG_VERBOSE << setw(20) << "pdf_factor[0]" << " = " << setprecision(20) << setw(28) << pdf_factor[0] << "   " << double2hexastr(pdf_factor[0]) << endl;
  logger << LOG_DEBUG_VERBOSE << setw(20) << "var_rel_alpha_S" << " = " << setprecision(20) << setw(28) << var_rel_alpha_S << "   " << double2hexastr(var_rel_alpha_S) << endl;
  logger << LOG_DEBUG_VERBOSE << setw(20) << "psi_ps_factor * rescaling_factor_alpha_e" << " = " << setprecision(20) << setw(28) << psi_ps_factor * rescaling_factor_alpha_e << "   " << double2hexastr(psi_ps_factor * rescaling_factor_alpha_e) << endl;
  logger << LOG_DEBUG_VERBOSE << setw(20) << "ME2" << " = " << setprecision(20) << setw(28) << ME2 << "   " << double2hexastr(ME2) << endl;

  if (switch_CV != 0){
    for (int s = 0; s < n_scales_CV; s++){
      integrand_CV[s] = psi_ps_factor * rescaling_factor_alpha_e * pdf_factor_CV[s][0] * var_rel_alpha_S_CV[s] * ME2;
    }
  }
    
  if (switch_distribution != 0){
    for (int id = 0; id < 3; id++){
      integrand_D[id][0] = psi_ps_factor * rescaling_factor_alpha_e * pdf_factor[id] * var_rel_alpha_S * ME2;
      if (switch_CV != 0){
	for (int s = 0; s < n_scales_CV; s++){
	  integrand_D_CV[s][id][0] = psi_ps_factor * rescaling_factor_alpha_e * pdf_factor_CV[s][id] * var_rel_alpha_S_CV[s] * ME2;
	}
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_integrand_VA(phasespace_set & psi){
  Logger logger("observable_set::determined_integrand_VA");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  // !!! reference does not work here !!!
  
  integrand = psi_ps_factor * rescaling_factor_alpha_e * pdf_factor[0] * var_rel_alpha_S * (VA_V_ME2 + VA_X_ME2 + VA_I_ME2);
  if (switch_CV != 0){
    for (int s = 0; s < n_scales_CV; s++){
      integrand_CV[s] = psi_ps_factor * rescaling_factor_alpha_e * pdf_factor_CV[s][0] * var_rel_alpha_S_CV[s] * (VA_V_ME2 + VA_X_ME2_CV[s] + VA_I_ME2);
    }
  }

  if (switch_distribution != 0){
    for (int id = 0; id < 3; id++){
      integrand_D[id][0] = psi_ps_factor * rescaling_factor_alpha_e * pdf_factor[id] * var_rel_alpha_S * (VA_V_ME2 + VA_X_ME2 + VA_I_ME2);
      if (switch_CV != 0){
	for (int s = 0; s < n_scales_CV; s++){
	  integrand_D_CV[s][id][0] = psi_ps_factor * rescaling_factor_alpha_e * pdf_factor_CV[s][id] * var_rel_alpha_S_CV[s] * (VA_V_ME2 + VA_X_ME2_CV[s] + VA_I_ME2);
	}
      }
    }
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_integrand_CA(phasespace_set & psi){//, vector<vector<collinear_set> > _collinear
  Logger logger("observable_set::determine_integrand_CA");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  //  for (int sd = 0; sd < value_mu_ren.size(); sd++){
  //    for (int ss = 0; ss < value_mu_ren[sd].size(); ss++){
  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    for (int ss = 0; ss < value_mu_fact[sd].size(); ss++){
      for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
	logger << LOG_DEBUG_VERBOSE << "sd = " << sd << "   ss = " << ss << "   i_a = " << i_a << endl;
	double temp_z = CA_value_ME2_KP[sd][ss][0][i_a] / psi_g_z_coll[(*CA_collinear)[i_a][0].no_emitter()];
	//	double temp_z = psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[sd][ss] * CA_value_ME2_KP[sd][ss][0][i_a] / psi_g_z_coll[(*CA_collinear)[i_a][0].no_emitter()];
	CA_value_integrand[sd][ss][i_a] = temp_z * CA_value_pdf_factor[sd][ss][i_a][(*CA_collinear)[i_a][0].no_emitter()][0];
	if (switch_distribution != 0){
	  for (int id = 0; id < 3; id++){
	    CA_value_integrand_D[sd][ss][id][i_a] = temp_z * CA_value_pdf_factor[sd][ss][i_a][(*CA_collinear)[i_a][0].no_emitter()][id];
	  }
	}
	if ((*CA_collinear)[i_a][0].in_collinear()[0] == 1){
	  double temp_1 = (CA_value_ME2_KP[sd][ss][1][i_a] / psi_g_z_coll[(*CA_collinear)[i_a][0].no_emitter()] + CA_value_ME2_KP[sd][ss][2][i_a]);
	  //	  double temp_1 = psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[sd][ss] * (CA_value_ME2_KP[sd][ss][1][i_a] / psi_g_z_coll[(*CA_collinear)[i_a][0].no_emitter()] + CA_value_ME2_KP[sd][ss][2][i_a]);
	  CA_value_integrand[sd][ss][i_a] += temp_1 * CA_value_pdf_factor[sd][ss][i_a][0][0];
	  if (switch_distribution != 0){
	    for (int id = 0; id < 3; id ++){
	      CA_value_integrand_D[sd][ss][id][i_a] += temp_1 * CA_value_pdf_factor[sd][ss][i_a][0][id];
	    }
	  }
	logger << LOG_DEBUG_VERBOSE << "sd = " << sd << "   ss = " << ss << "   i_a = " << i_a << endl;
	}
      }
      // could be multiplied only here:
      // psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[sd][ss]
      CA_sum_value_integrand[sd][ss] = accumulate(CA_value_integrand[sd][ss].begin(), CA_value_integrand[sd][ss].end(), 0.);
      if (switch_distribution != 0){
	for (int iy = 0; iy < 3; iy ++){
	  CA_sum_value_integrand_D[sd][ss][iy] = accumulate(CA_value_integrand_D[sd][ss][iy].begin(), CA_value_integrand_D[sd][ss][iy].end(), 0.);
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "after CA_sum_value_integrand evaluation" << endl;
    
  integrand = psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[dynamic_scale][map_value_scale_ren] * CA_sum_value_integrand[dynamic_scale][map_value_scale_fact];

  logger << LOG_DEBUG_VERBOSE << "after integrand evaluation" << endl;

  if (munich_isnan(integrand) || munich_isinf(integrand)) {
    logger << LOG_ERROR << "integrand is nan/inf" << endl;
    int sd = dynamic_scale;
    int ss = map_value_scale_fact;
    for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
      double temp_1 = psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[sd][ss] * (CA_value_ME2_KP[sd][ss][1][i_a] / psi_g_z_coll[(*CA_collinear)[i_a][0].no_emitter()] + CA_value_ME2_KP[sd][ss][2][i_a]);
      logger << LOG_ERROR << i_a << ": temp_1=" << temp_1 << endl;
      logger << LOG_ERROR << value_factor_alpha_S[sd][ss] << ", " << CA_value_ME2_KP[sd][ss][1][i_a] << ", " << psi_g_z_coll[(*CA_collinear)[i_a][0].no_emitter()] << ", " << CA_value_ME2_KP[sd][ss][2][i_a] << endl;
      logger << LOG_ERROR << "PDF factor=" << CA_value_pdf_factor[sd][ss][i_a][0][0] << endl;
    }
  }

  if (switch_CV){
    for (int s = 0; s < n_scales_CV; s++){
      integrand_CV[s] = psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[dynamic_scale_CV][map_value_scale_ren_CV[s]] * CA_sum_value_integrand[dynamic_scale_CV][map_value_scale_fact_CV[s]];
      logger << LOG_DEBUG_VERBOSE << "CA_sum_value_integrand[dynamic_scale_CV = " << dynamic_scale_CV << "][map_value_scale_fact_CV[" << s << "] = " << map_value_scale_fact_CV[s] << "] = " << CA_sum_value_integrand[dynamic_scale_CV][map_value_scale_fact_CV[s]] << endl;
    }
  }

  if (switch_distribution != 0){
    for (int iy = 0; iy < 3; iy++){
      integrand_D[iy][0] = psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[dynamic_scale][map_value_scale_ren] * CA_sum_value_integrand_D[dynamic_scale][map_value_scale_fact][iy];
      if (switch_CV != 0){
	for (int s = 0; s < n_scales_CV; s++){
	  integrand_D_CV[s][iy][0] = psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[dynamic_scale_CV][map_value_scale_ren_CV[s]] * CA_sum_value_integrand_D[dynamic_scale_CV][map_value_scale_fact_CV[s]][iy];
	}
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



/*
void observable_set::determine_integrand_CA(phasespace_set & psi){//, vector<vector<collinear_set> > _collinear
  Logger logger("observable_set::determine_integrand_CA");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int sd = 0; sd < value_mu_ren.size(); sd++){
    for (int ss = 0; ss < value_mu_ren[sd].size(); ss++){
      for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
	double temp_z = psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[sd][ss] * CA_value_ME2_KP[sd][ss][0][i_a] / psi_g_z_coll[(*CA_collinear)[i_a][0].no_emitter()];
	CA_value_integrand[sd][ss][i_a] = temp_z * CA_value_pdf_factor[sd][ss][i_a][(*CA_collinear)[i_a][0].no_emitter()][0];
	if (switch_distribution != 0){
	  for (int id = 0; id < 3; id++){
	    CA_value_integrand_D[sd][ss][id][i_a] = temp_z * CA_value_pdf_factor[sd][ss][i_a][(*CA_collinear)[i_a][0].no_emitter()][id];
	  }
	}
	if ((*CA_collinear)[i_a][0].in_collinear()[0] == 1){
	  double temp_1 = psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[sd][ss] * (CA_value_ME2_KP[sd][ss][1][i_a] / psi_g_z_coll[(*CA_collinear)[i_a][0].no_emitter()] + CA_value_ME2_KP[sd][ss][2][i_a]);
	  CA_value_integrand[sd][ss][i_a] += temp_1 * CA_value_pdf_factor[sd][ss][i_a][0][0];
	  if (switch_distribution != 0){
	    for (int id = 0; id < 3; id ++){
	      CA_value_integrand_D[sd][ss][id][i_a] += temp_1 * CA_value_pdf_factor[sd][ss][i_a][0][id];
	    }
	  }
	}
      }
      // could be multiplied only here:
      // psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[sd][ss]
      CA_sum_value_integrand[sd][ss] = accumulate(CA_value_integrand[sd][ss].begin(), CA_value_integrand[sd][ss].end(), 0.);
      if (switch_distribution != 0){
	for (int iy = 0; iy < 3; iy ++){
	  CA_sum_value_integrand_D[sd][ss][iy] = accumulate(CA_value_integrand_D[sd][ss][iy].begin(), CA_value_integrand_D[sd][ss][iy].end(), 0.);
	}
      }
    }
  }
    
  integrand = CA_sum_value_integrand[dynamic_scale][map_value_scale_fact];

  if (munich_isnan(integrand) || munich_isinf(integrand)) {
    logger << LOG_ERROR << "integrand is nan/inf" << endl;
    int sd = dynamic_scale;
    int ss = map_value_scale_fact;
    for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
      double temp_1 = psi_ps_factor * rescaling_factor_alpha_e * value_factor_alpha_S[sd][ss] * (CA_value_ME2_KP[sd][ss][1][i_a] / psi_g_z_coll[(*CA_collinear)[i_a][0].no_emitter()] + CA_value_ME2_KP[sd][ss][2][i_a]);
      logger << LOG_ERROR << i_a << ": temp_1=" << temp_1 << endl;
      logger << LOG_ERROR << value_factor_alpha_S[sd][ss] << ", " << CA_value_ME2_KP[sd][ss][1][i_a] << ", " << psi_g_z_coll[(*CA_collinear)[i_a][0].no_emitter()] << ", " << CA_value_ME2_KP[sd][ss][2][i_a] << endl;
      logger << LOG_ERROR << "PDF factor=" << CA_value_pdf_factor[sd][ss][i_a][0][0] << endl;
    }
  }

  if (switch_CV){
    for (int s = 0; s < n_scales_CV; s++){
      integrand_CV[s] = CA_sum_value_integrand[dynamic_scale_CV][map_value_scale_fact_CV[s]];
      logger << LOG_DEBUG << "CA_sum_value_integrand[dynamic_scale_CV = " << dynamic_scale_CV << "][map_value_scale_fact_CV[" << s << "] = " << map_value_scale_fact_CV[s] << "] = " << CA_sum_value_integrand[dynamic_scale_CV][map_value_scale_fact_CV[s]] << endl;
    }
  }

  if (switch_distribution != 0){
    for (int iy = 0; iy < 3; iy++){
      integrand_D[iy][0] = CA_sum_value_integrand_D[dynamic_scale][map_value_scale_fact][iy];
      if (switch_CV != 0){
	for (int s = 0; s < n_scales_CV; s++){
	  integrand_D_CV[s][iy][0] = CA_sum_value_integrand_D[dynamic_scale_CV][map_value_scale_fact_CV[s]][iy];
	}
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/


void observable_set::determine_integrand_RA(phasespace_set & psi){
  Logger logger("observable_set::determined_integrand_RA");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  //  psi_ps_factor * rescaling_factor_alpha_e = psi_hcf / (psi_g_tot * psi_xbs_all[0][0]);
  for (int i_a = 0; i_a < n_ps; i_a++){
    //  for (int i_a = 0; i_a < n_dipoles; i_a++){
    //    logger << LOG_INFO << "RA_rel_alpha_S[" << i_a << "] = " << RA_rel_alpha_S[i_a] << endl;
    //    logger << LOG_INFO << "RA_pdf_factor[" << i_a << "][0] = " << RA_pdf_factor[i_a][0] << endl;
    var_RA_ME2[i_a] = RA_pdf_factor[i_a][0] * RA_rel_alpha_S[i_a] * RA_ME2[i_a];
    logger << LOG_DEBUG_VERBOSE << "RA_pdf_factor[" << i_a << "][0] = " << setw(23) << setprecision(15) << RA_pdf_factor[i_a][0] << "   RA_rel_alpha_S[" << i_a << "] = " << setw(23) << setprecision(15) << RA_rel_alpha_S[i_a] << "   RA_ME2[" << i_a << "] = " << setw(23) << setprecision(15) << RA_ME2[i_a] << endl;
    if (switch_CV != 0){
      for (int s = 0; s < n_scales_CV; s++){
	var_RA_ME2_CV[s][i_a] = RA_pdf_factor_CV[i_a][s][0] * RA_rel_alpha_S_CV[i_a][s] * RA_ME2[i_a];
      }
    }
    if (switch_distribution != 0){
      for (int j = 0; j < 3; j++){
	integrand_D[j][i_a] = RA_pdf_factor[i_a][j] * RA_rel_alpha_S[i_a] * RA_ME2[i_a] * psi_ps_factor * rescaling_factor_alpha_e;
	if (switch_CV != 0){
	  for (int s = 0; s < n_scales_CV; s++){
	    integrand_D_CV[s][j][i_a] = RA_pdf_factor_CV[i_a][s][j] * RA_rel_alpha_S_CV[i_a][s] * RA_ME2[i_a] * psi_ps_factor * rescaling_factor_alpha_e;
	  }
	}
      }
    }
  }

  integrand = accumulate(var_RA_ME2.begin(), var_RA_ME2.end(), 0.) * psi_ps_factor * rescaling_factor_alpha_e;

  logger << LOG_DEBUG_VERBOSE << "integrand = " << setw(24) << setprecision(15) << integrand << endl;

  if (switch_CV != 0){
    for (int s = 0; s < n_scales_CV; s++){
      integrand_CV[s] = accumulate(var_RA_ME2_CV[s].begin(), var_RA_ME2_CV[s].end(), 0.) * psi_ps_factor * rescaling_factor_alpha_e;
    }
  }

  if (n_moments != 0){
    for (int nm = 0; nm < moment.size(); nm++){
      for (int i_a = 0; i_a < moment[nm].size(); i_a++){
	RA_ME2_moment[nm][i_a] = var_RA_ME2[i_a] * moment[nm][i_a];
	if (switch_CV != 0){
	  for (int s = 0; s < n_scales_CV; s++){
	    RA_ME2_moment_CV[s][nm][i_a] = var_RA_ME2_CV[s][i_a] * moment[nm][i_a];
	  }
	}
      }
      RA_integrand_moment[nm] = accumulate(RA_ME2_moment[nm].begin(), RA_ME2_moment[nm].end(), 0.) * psi_ps_factor * rescaling_factor_alpha_e;
      
      if (switch_CV != 0){
	for (int s = 0; s < n_scales_CV; s++){
	  RA_integrand_moment_CV[s][nm] = accumulate(RA_ME2_moment_CV[s][nm].begin(), RA_ME2_moment_CV[s][nm].end(), 0.) * psi_ps_factor * rescaling_factor_alpha_e;
	}
      }
    }
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_psp_weight(){
  Logger logger("observable_set::determine_psp_weight");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  logger << LOG_DEBUG_VERBOSE << "psp_weight = " << setprecision(20) << setw(28) << integrand << "   " << double2hexastr(integrand) << endl;
  
  this_psp_weight = integrand;
  this_psp_weight2 = pow(this_psp_weight, 2.);
  step_sum_weight += this_psp_weight;
  step_sum_weight2 += this_psp_weight2;

  logger << LOG_DEBUG_VERBOSE << "step_sum_weight = " << step_sum_weight << "   this_psp_weight = " << this_psp_weight << endl;

  if (switch_CV != 0){
    for (int c = 0; c < n_qTcut; c++){
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	if (c <= cut_ps[0]){
	  this_psp_weight_CV[c][i_s] = integrand_CV[i_s];
	  this_psp_weight2_CV[c][i_s] = pow(this_psp_weight_CV[c][i_s], 2.);
	  step_sum_weight_CV[c][i_s] += this_psp_weight_CV[c][i_s];
	  step_sum_weight2_CV[c][i_s] += this_psp_weight2_CV[c][i_s];
	}
      }
    }
  }

  if (n_moments != 0){
    for (int i_m = 0; i_m < moment.size(); i_m++){
      for (int i_a = 0; i_a < moment[i_m].size(); i_a++){
	this_psp_moment[i_m] = this_psp_weight * moment[i_m][i_a];
	this_psp_moment2[i_m] = pow(this_psp_moment[i_m], 2.);
	step_sum_moment[i_m] += this_psp_moment[i_m];
	step_sum_moment2[i_m] += this_psp_moment2[i_m];
	if (switch_CV != 0){
	  for (int c = 0; c < n_qTcut; c++){
	    for (int i_s = 0; i_s < n_scales_CV; i_s++){
	      if (c <= cut_ps[0]){
		this_psp_moment_CV[i_m][c][i_s] = this_psp_weight_CV[c][i_s] * moment[i_m][i_a];
		this_psp_moment2_CV[i_m][c][i_s] = pow(this_psp_moment_CV[i_m][c][i_s], 2.);
		step_sum_moment_CV[i_m][c][i_s] += this_psp_moment_CV[i_m][c][i_s];
		step_sum_moment2_CV[i_m][c][i_s] += this_psp_moment2_CV[i_m][c][i_s];
	      }
	    }
	  }
	}
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_psp_weight_RA(phasespace_set & psi){
  Logger logger("observable_set::determine_psp_weight_RA");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  this_psp_weight = integrand;
  this_psp_weight2 = pow(this_psp_weight, 2.);
  step_sum_weight += this_psp_weight;
  step_sum_weight2 += this_psp_weight2;

  if (switch_CV){
    for (int c = 0; c < n_qTcut; c++){
      for (int s = 0; s < n_scales_CV; s++){
	var_R_integrand_cut[c][s] = 0.;
	var_A_integrand_cut[c][s] = 0.;
      }
    }
    for (int s = 0; s < n_scales_CV; s++){
      for (int ia = 1; ia < n_ps; ia++){
	for (int c = 0; c < n_qTcut; c++){
	  if (c <= cut_ps[ia])
	    var_A_integrand_cut[c][s] += var_RA_ME2_CV[s][ia];
	}
      }
      for (int c = 0; c < n_qTcut; c++){
	if (c <= cut_ps[0])
	  var_R_integrand_cut[c][s] = var_RA_ME2_CV[s][0];
      }
    }
    for (int c = 0; c < n_qTcut; c++){
      for (int s = 0; s < n_scales_CV; s++){
	this_psp_weight_CV[c][s] = (var_R_integrand_cut[c][s] + var_A_integrand_cut[c][s]) * psi_ps_factor * rescaling_factor_alpha_e;// sign !!!
	this_psp_weight2_CV[c][s] = pow(this_psp_weight_CV[c][s], 2);
	step_sum_weight_CV[c][s] += this_psp_weight_CV[c][s];
	step_sum_weight2_CV[c][s] += this_psp_weight2_CV[c][s];
      }
    }
  }
  /*
  if (psi_RA_techcut == 0){
    this_psp_weight = integrand;
    this_psp_weight2 = pow(this_psp_weight, 2.);
    step_sum_weight += this_psp_weight;
    step_sum_weight2 += this_psp_weight2;
  }
  else {
    this_psp_weight = 0.;
    this_psp_weight2 = 0.;
    this_psp_weight = 0.;
    this_psp_weight2 = 0.;
  }
  if (switch_CV){
    if (psi_RA_techcut == 0){
      for (int c = 0; c < n_qTcut; c++){
	for (int s = 0; s < n_scales_CV; s++){
	  var_R_integrand_cut[c][s] = 0.;
	  var_A_integrand_cut[c][s] = 0.;
	  //	var_R_integrand_cut_incl[c][s] = 0.;
	}
      }
      for (int s = 0; s < n_scales_CV; s++){
	for (int ia = 1; ia < n_ps; ia++){
	  for (int c = 0; c < n_qTcut; c++){
	    if (c <= cut_ps[ia])
	      var_A_integrand_cut[c][s] += var_RA_ME2_CV[s][ia];
	  }
	}
	for (int c = 0; c < n_qTcut; c++){
	  if (c <= cut_ps[0])
	    var_R_integrand_cut[c][s] = var_RA_ME2_CV[s][0];
	}
      }
      for (int c = 0; c < n_qTcut; c++){
	for (int s = 0; s < n_scales_CV; s++){
	  this_psp_weight_CV[c][s] = (var_R_integrand_cut[c][s] + var_A_integrand_cut[c][s]) * psi_ps_factor * rescaling_factor_alpha_e;// sign !!!
	  this_psp_weight2_CV[c][s] = pow(this_psp_weight_CV[c][s], 2);
	  step_sum_weight_CV[c][s] += this_psp_weight_CV[c][s];
	  step_sum_weight2_CV[c][s] += this_psp_weight2_CV[c][s];
	}
      }
    }
    else{
      for (int c = 0; c < n_qTcut; c++){
	for (int s = 0; s < n_scales_CV; s++){
	  this_psp_weight_CV[c][s] = 0.;
	  this_psp_weight2_CV[c][s] = 0.;
	}
      }
    }
  }
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



