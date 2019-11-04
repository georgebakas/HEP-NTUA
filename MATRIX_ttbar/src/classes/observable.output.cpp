#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
void observable_set::header_integration(ofstream & out_integration, phasespace_set & psi){
  Logger logger("observable_set::header_integration");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  out_integration << setprecision(15);
  out_integration << left << setw(30) << "process" << " = " << csi->process_class << endl;
  if (csi->process_class != csi->basic_process_class){
    out_integration << left << setw(30) << "basic process" << " = " << csi->basic_process_class << endl;
  }
  if (coll_choice == 0){out_integration << setw(30) << "collision type" << " = " << "partonic" << endl;}
  if (coll_choice == 1){out_integration << setw(30) << "collision type" << " = " << "proton-proton" << endl;}
  if (coll_choice == 2){out_integration << setw(30) << "collision type" << " = " << "proton-antiproton" << endl;}
  out_integration << left << setw(30) << "beam energy" << " = " << psi_E << " GeV" << endl;
  out_integration << left << setw(30) << "partonic subprocess" << " = " << csi->subprocess << endl;

  out_integration << endl;
  out_integration << "Particle masses: " << endl;
  out_integration << left << setw(30) << "M_W" << " = " << msi.M_W << " GeV" << endl;
  out_integration << left << setw(30) << "M_Z" << " = " << msi.M_Z << " GeV" << endl;
  out_integration << left << setw(30) << "M_H" << " = " << msi.M_H << " GeV" << endl;

  out_integration << left << setw(30) << "M_d" << " = " << msi.M_d << " GeV" << endl;
  out_integration << left << setw(30) << "M_u" << " = " << msi.M_u << " GeV" << endl;
  out_integration << left << setw(30) << "M_s" << " = " << msi.M_s << " GeV" << endl;
  out_integration << left << setw(30) << "M_c" << " = " << msi.M_c << " GeV" << endl;
  out_integration << left << setw(30) << "M_b" << " = " << msi.M_b << " GeV" << endl;
  out_integration << left << setw(30) << "M_t" << " = " << msi.M_t << " GeV" << endl;
  out_integration << endl;
  out_integration << "Particle widths: " << endl;
  out_integration << left << setw(30) << "Gamma_W" << " = " << msi.Gamma_W << " GeV" << endl;
  out_integration << left << setw(30) << "Gamma_Z" << " = " << msi.Gamma_Z << " GeV" << endl;
  out_integration << left << setw(30) << "Gamma_H" << " = " << msi.Gamma_H << " GeV" << endl;
  out_integration << left << setw(30) << "Gamma_t" << " = " << msi.Gamma_t << " GeV" << endl;
  out_integration << endl;
  out_integration << "Branching ratios: " << endl;
  out_integration << left << setw(30) << "BR_Wlv" << " = " << msi.BR_Wlv << endl;
  out_integration << left << setw(30) << "BR_Zll" << " = " << msi.BR_Zll << endl;
  out_integration << left << setw(30) << "BR_tWb" << " = " << msi.BR_tWb << endl;
  out_integration << endl;
  out_integration << "Electroweak Standard Model parameters: " << endl;
  //  out_integration << left << setw(30) << "1 / alpha_e" << " = " << 1 / msi.alpha_e << endl;
  //  out_integration << left << setw(30) << "alpha_e" << " = " << msi.alpha_e << endl;

  if (msi.ew_scheme == 0 || msi.use_adapted_ew_coupling == 0){
    out_integration << left << setw(30) << "1 / alpha_e_0" << " = " << setw(30) << 1 / msi.alpha_e_0 << setw(5) << setw(30) << "alpha_e_0" << " = " << msi.alpha_e_0 << endl;
  }
  if (msi.ew_scheme == 1 || msi.ew_scheme == -1 || msi.use_adapted_ew_coupling == 1){
    out_integration << left << setw(30) << "G_F" << " = " << msi.G_F << endl;
    out_integration << left << setw(30) << "1 / alpha_e_Gmu" << " = " << setw(30) << 1 / msi.alpha_e_Gmu << setw(5) << setw(30) << "alpha_e_Gmu" << " = " << msi.alpha_e_Gmu << endl;
  }
  if (msi.ew_scheme == 2 || msi.use_adapted_ew_coupling == 2){
    out_integration << left << setw(30) << "1 / alpha_e_MZ" << " = " << setw(30) << 1 / msi.alpha_e_MZ << setw(5) << setw(30) << "alpha_e_MZ" << " = " << msi.alpha_e_MZ << endl;
  }
  out_integration << endl;
  
  out_integration << left << setw(30) << "ew_scheme" << " = " << msi.ew_scheme;
  if (msi.ew_scheme == 0){out_integration << setw(24) << " (alpha(0) scheme)" << "->   alpha_e = alpha_e_0 = " << msi.alpha_e << endl;}
  else if (msi.ew_scheme == 1 || msi.ew_scheme == -1){out_integration << setw(24) << " (G_mu scheme)" << "->   alpha_e = alpha_e_Gmu = " << msi.alpha_e << endl;}
  else if (msi.ew_scheme == 2){out_integration << setw(24) << " (alpha(M_Z) scheme)" << "->   alpha_e = alpha_e_MZ = " << msi.alpha_e << endl;}
  else {out_integration << " (No allowed scheme specified!)" << endl; exit(1);}
  //  out_integration << left << setw(30) << "1 / alpha_e" << " = " << setw(30) << 1 / msi.alpha_e << setw(5) << setw(30) << "alpha_e" << " = " << msi.alpha_e << endl;
  if (msi.use_adapted_ew_coupling != -1){
    out_integration << left << setw(30) << "use_adapted_ew_coupling" << " = " << msi.use_adapted_ew_coupling << "   ->   ";
    if (msi.use_adapted_ew_coupling == 0){
      if (msi.ew_scheme == 1 || msi.ew_scheme == -1){
	out_integration << "Amplitudes are rescaled by a factor of (alpha_e_0 / alpha_e_Gmu) ^ " << csi->n_photon_born << " ." << endl;
      }
      else if (msi.ew_scheme == 2){
	out_integration << "Amplitudes are rescaled by a factor of (alpha_e_0 / alpha_e_MZ) ^ " << csi->n_photon_born << " ." << endl;
      }
      else {exit(1);}
    }
    else {
      int temp_exp = 0;
      if (csi->type_correction == "QEW" || csi->type_correction == "MIX"){temp_exp = csi->contribution_order_alpha_e - 1 - csi->n_photon_born;}
      else {temp_exp = csi->contribution_order_alpha_e - csi->n_photon_born;}
      if (msi.use_adapted_ew_coupling == 1){
	out_integration << "Amplitudes are rescaled by a factor of (alpha_e_Gmu / alpha_e_0) ^ " << temp_exp << " ." << endl;
      }
      else if (msi.use_adapted_ew_coupling == 2){
	out_integration << "Amplitudes are rescaled by a factor of (alpha_e_MZ / alpha_e_0) ^ " << temp_exp << " ." << endl;
	
      }
    }
  }
  //  out_integration << left << setw(30) << "use_adapted_ew_coupling" << " = " << msi.use_adapted_ew_coupling << endl;

  out_integration << left << setw(30) << "e" << " = " << msi.e_pow[1] << endl;
  out_integration << left << setw(30) << "cos_w" << " = " << msi.cos_w << endl;
  out_integration << left << setw(30) << "sin_w" << " = " << msi.sin_w << endl;
  out_integration << left << setw(30) << "cos2_w" << " = " << msi.cos2_w << endl;
  out_integration << left << setw(30) << "sin2_w" << " = " << msi.sin2_w << endl;
  //  out_integration << endl;

  /*
  if (ckm_choice == 0){
    out_integration << left << setw(30) << "V_ckm" << " = " << "trivial" << endl;
  }
  else {
    for (int i1 = 1; i1 < 4; i1++){
      if (i1 == 2){out_integration << left << setw(30) << "V_ckm" << " = ";}
      else {out_integration << left << setw(30) << "";}
      for (int i2 = 1; i2 < 4; i2++){
	out_integration << setprecision(8) << setw(20) << msi.V_ckm[i1][i2];
      }
      out_integration << endl;
    }
  }
  */
  out_integration << endl;


  out_integration << "QCD parameters: " << endl;
  if (coll_choice != 0){
   out_integration << left << setw(30) << "LHAPDF set" << " = "  << LHAPDFname << endl;
   out_integration << left << setw(30) << "LHAPDF subset" << " = "  << LHAPDFsubset << endl;
   out_integration << left << setw(30) << "factorization scale" << " = " << scale_fact << " GeV" << endl;
   out_integration << endl;
  }
  out_integration << left << setw(30) << "renormalization scale" << " = " << scale_ren << " GeV" << endl;
  out_integration << endl;
  out_integration << left << setw(30) << "N_f" << " = " << N_f << endl;
  out_integration << left << setw(30) << "N_f_active" << " = " << N_f_active << endl;
  stringstream temp_alpha_S;
  temp_alpha_S << "alpha_S(mu" << " = " << setprecision(6) << mu_ren[1] << ")";
  out_integration << left << setw(30) << temp_alpha_S.str() << " = " << setprecision(15) << alpha_S << endl;
  //  alpha_S = LHAPDF::alphasPDF(mu_ren[1]);
  stringstream temp_g_S;
  temp_g_S << "g_S(mu" << " = " << setprecision(6) << mu_ren[1] << ")";
  out_integration << left << setw(30) << temp_g_S.str() << " = " << setprecision(15) << g_from_alpha(alpha_S) << endl;
  out_integration << endl;

  out_integration << "Definitions for event selection: " << endl;
  out_integration << setw(5) << "no" << setw(20) << "object" << setw(20) << "n_observed_min" << setw(20) << "n_observed_max" << setw(20) << "define_pT_min" << setw(20) << "define_ET_min" << setw(20) << "define_|eta|_max"  << setw(20) << "define_|y|_max"<< endl;
  //  out_integration << setw(5) << right << "no" << setw(10) << "object" << setw(10) << "n_min" << setw(10) << "n_max" << setw(10) << "pT" << setw(10) << "ET" << setw(10) << "|eta|"  << setw(10) << "|y|"<< endl;

  for (int i = 1; i < esi.object_list.size(); i++){
    out_integration << setw(5) << i << setw(20) << esi.object_list[i] << setw(20) << esi.pda[i].n_observed_min << setw(20) << esi.pda[i].n_observed_max << setw(20) << esi.pda[i].define_pT << setw(20) << esi.pda[i].define_ET << setw(20) << esi.pda[i].define_eta << setw(20) << esi.pda[i].define_y << endl;
    //    out_integration << setw(5) << right << i << setw(10) << esi.object_list[i] << setw(10) << esi.pda[i].n_observed_min << setw(10) << esi.pda[i].n_observed_max << setw(10) << esi.pda[i].define_pT << setw(10) << esi.pda[i].define_ET << setw(10) << esi.pda[i].define_eta << setw(10) << esi.pda[i].define_y << endl;

  }
  out_integration << endl;
  
  out_integration << "Relevant object definitions for event selection: -> oset" << endl;
  
  out_integration << setw(5) << "no" << setw(20) << "object" << setw(20) << "n_observed_min" << setw(20) << "n_observed_max" << setw(20) << "define_pT_min" << setw(20) << "define_ET_min" << setw(20) << "define_|eta|_max"  << setw(20) << "define_|y|_max"<< endl;
  ///  out_integration << setw(5) << right << "no" << setw(10) << "object" << setw(10) << "access" << setw(10) << "n_min" << setw(10) << "n_max" << setw(10) << "pT" << setw(10) << "ET" << setw(10) << "|eta|"  << setw(10) << "|y|"<< endl;

  for (int i = 1; i < (esi.object_list_selection).size(); i++){
    out_integration << setw(5) << i << setw(20) << esi.object_list_selection[i] << setw(20) << esi.pds[i].n_observed_min << setw(20) << esi.pds[i].n_observed_max << setw(20) << esi.pds[i].define_pT << setw(20) << esi.pds[i].define_ET << setw(20) << esi.pds[i].define_eta << setw(20) << esi.pds[i].define_y << endl;
    ///    out_integration << setw(5) << right << i << setw(10) << esi.object_list_selection[i] << setw(10) << access_object[esi.object_list_selection[i]] << setw(10) << esi.pds[i].n_observed_min << setw(10) << esi.pds[i].n_observed_max << setw(10) << esi.pds[i].define_pT << setw(10) << esi.pds[i].define_ET << setw(10) << esi.pds[i].define_eta << setw(10) << esi.pds[i].define_y << endl;
  }

  /*
  for (int i = 1; i < (relevant_object_list).size(); i++){
    out_integration << setw(5) << i << setw(20) << relevant_object_list[i] << setw(20) << relevant_n_observed_min[i] << setw(20) << relevant_n_observed_max[i] << setw(20) << relevant_define_pT[i] << setw(20) << relevant_define_ET[i] << setw(20) << relevant_define_eta[i] << setw(20) << relevant_define_y[i] << endl;
  }
  */  
  out_integration << endl;

  out_integration << "Additional cuts: " << endl;
  out_integration << endl;
  out_integration << "Integration parameters: " << endl;
  out_integration << setw(30) << "n_events_max" << " = " << psi_n_events_max << endl;
  out_integration << setw(30) << "n_events_min" << " = " << psi_n_events_min << endl;
  out_integration << setw(30) << "n_step" << " = " << psi_n_step << endl;
  out_integration << endl;
  out_integration << setw(30) << "mass0" << " = " << psi_mass0 << endl;
  out_integration << setw(30) << "nuxs" << " = " << psi_nuxs << endl;
  out_integration << setw(30) << "nuxt" << " = " << psi_nuxt << endl;
  out_integration << setw(30) << "cut_technical" << " = " << psi_cut_technical << endl;
  out_integration << setw(30) << "exp_pdf" << " = " << psi_exp_pdf << endl;
  if (csi->type_contribution == "RA" || csi->type_contribution == "RRA"){
    out_integration << setw(30) << "exp_ij_k_y" << " = " << psi_exp_ij_k_y << endl;
    out_integration << setw(30) << "exp_ij_k_z" << " = " << psi_exp_ij_k_z << endl;
    out_integration << setw(30) << "exp_ij_a_x" << " = " << psi_exp_ij_a_x << endl;
    out_integration << setw(30) << "exp_ij_a_z" << " = " << psi_exp_ij_a_z << endl;
    out_integration << setw(30) << "exp_ai_k_x" << " = " << psi_exp_ai_k_x << endl;
    out_integration << setw(30) << "exp_ai_k_u" << " = " << psi_exp_ai_k_u << endl;
    out_integration << setw(30) << "exp_ai_b_x" << " = " << psi_exp_ai_b_x << endl;
    out_integration << setw(30) << "exp_ai_b_v" << " = " << psi_exp_ai_b_v << endl;
  }
  out_integration << endl;
  out_integration << endl;

  out_integration << setw(30) << "switch_MC" << " = " << psi_switch_MC << endl;
  out_integration << setw(30) << "MC_n_channel" << " = " << psi_MC_n_channel << endl;
  out_integration << setw(30) << "n_alpha_events" << " = " << psi_n_alpha_events << endl;
  out_integration << setw(30) << "n_alpha_steps" << " = " << psi_n_alpha_steps << endl;
  out_integration << setw(30) << "MCweight_min" << " = " << psi_MCweight_min << endl;
  out_integration << setw(30) << "MCweight_limit_min" << " = " << psi_MCweight_limit_min << endl;
  out_integration << setw(30) << "MCweight_limit_max" << " = " << psi_MCweight_limit_max << endl;
  out_integration << endl;

  if (csi->type_contribution == "RA" || csi->type_contribution == "RRA"){
     //  out_integration << setw(30) << "MC_n_channel" << " = " << psi_MC_n_channel << endl;
    for (int i_a = 0; i_a < psi_MC_n_channel_phasespace.size(); i_a++){
      out_integration << setw(21) << "MC_n_channel_phasespace[" << setw(3) << right << i_a << "]   = " << left << psi_MC_n_channel_phasespace[i_a] << endl;
    }
    out_integration << endl;
   }

  out_integration << setw(30) << "switch_MC_tau" << " = " << psi_switch_MC_tau << endl;
  for (int i_c = 0; i_c < psi_tau_MC_map.size(); i_c++){
    stringstream temp;
    temp <<    "tau_MC_map[" << i_c << "]";
    out_integration << setw(30) << temp.str() << " = " << psi_tau_MC_map[i_c] << endl;
  }
  out_integration << endl;

  if (csi->type_contribution == "RA" || csi->type_contribution == "RRA"){
    //  out_integration << setw(30) << "MC_n_channel" << " = " << psi_MC_n_channel << endl;
    out_integration << setw(30) << "switch_MC_x_dipole" << " = " << psi_switch_MC_x_dipole << endl;
    for (int i_a = 1; i_a < psi_MC_x_dipole_mapping.size(); i_a++){
      for (int i_m = 0; i_m < psi_MC_x_dipole_mapping[i_a].size(); i_m++){
	out_integration << setw(20) << "MC_x_dipole_mapping[" << setw(3) << right << i_a << "][" << setw(1) << right << i_m << "]    = " << left << psi_MC_x_dipole_mapping[i_a][i_m] << endl;
      }
    }
    out_integration << endl;
  }



  out_integration << setw(30) << "switch_IS_MC" << " = " << psi_switch_IS_MC << endl;
  out_integration << setw(30) << "weight_IS" << " = " << psi_weight_IS << endl;
  out_integration << setw(30) << "n_IS_gridsize" << " = " << psi_n_IS_gridsize << endl;
  out_integration << setw(30) << "n_IS_gridsize_p" << " = " << psi_n_IS_gridsize_p << endl;
  out_integration << setw(30) << "n_IS_gridsize_f" << " = " << psi_n_IS_gridsize_f << endl;
  out_integration << setw(30) << "n_IS_gridsize_t_t" << " = " << psi_n_IS_gridsize_t_t << endl;
  out_integration << setw(30) << "n_IS_gridsize_t_phi" << " = " << psi_n_IS_gridsize_t_phi << endl;
  out_integration << setw(30) << "n_IS_gridsize_d_cth" << " = " << psi_n_IS_gridsize_d_cth << endl;
  out_integration << setw(30) << "n_IS_gridsize_d_phi" << " = " << psi_n_IS_gridsize_d_phi << endl;
  out_integration << endl;
  out_integration << setw(30) << "switch_IS_tau" << " = " << psi_switch_IS_tau << endl;
  out_integration << setw(30) << "n_tau_bins" << " = " << psi_n_tau_bins << endl;
  out_integration << setw(30) << "n_tau_events" << " = " << psi_n_tau_events << endl;
  out_integration << setw(30) << "n_tau_steps" << " = " << psi_n_tau_steps << endl;
  out_integration << endl;




  out_integration << setw(30) << "switch_IS_x1x2" << " = " << psi_switch_IS_x1x2 << endl;
  out_integration << setw(30) << "n_x1x2_bins" << " = " << psi_n_x1x2_bins << endl;
  out_integration << setw(30) << "n_x1x2_events" << " = " << psi_n_x1x2_events << endl;
  out_integration << setw(30) << "n_x1x2_steps" << " = " << psi_n_x1x2_steps << endl;
  out_integration << endl;
  if (csi->type_contribution == "CA" || 
      csi->type_contribution == "RCA" || 
      csi->type_contribution == "RCJ" || 
      csi->type_contribution == "CT" || 
      csi->type_contribution == "CJ" || 
      csi->type_contribution == "CT2" || 
      csi->type_contribution == "CJ2" || 
      csi->type_contribution == "L2CT" || 
      csi->type_contribution == "L2CJ"){
    out_integration << setw(30) << "switch_IS_z1z2" << " = " << psi_switch_IS_z1z2 << endl;
    out_integration << setw(30) << "n_z1z2_bins" << " = " << psi_n_z1z2_bins << endl;
    out_integration << setw(30) << "n_z1z2_events" << " = " << psi_n_z1z2_events << endl;
    out_integration << setw(30) << "n_z1z2_steps" << " = " << psi_n_z1z2_steps << endl;
    out_integration << endl;
  }

  out_integration << endl;

  out_integration << "Relevant objects:" << endl;
  
  out_integration << setw(10) << "no" << setw(20) << "object" << setw(20) << "n_observed_min" << setw(20) << "n_observed_max" << setw(20) << "define_pT" << setw(20) << "define_ET" << setw(20) << "define_eta"  << setw(20) << "define_y"<< endl;
  for (int i = 0; i < (esi.object_list_selection).size(); i++){
    out_integration << setw(5) << i << setw(20) << esi.object_list_selection[i] << setw(20) << esi.pds[i].n_observed_min << setw(20) << esi.pds[i].n_observed_max << setw(20) << esi.pds[i].define_pT << setw(20) << esi.pds[i].define_ET << setw(20) << esi.pds[i].define_eta << setw(20) << esi.pds[i].define_y << endl;
  }

  out_integration << endl;

  out_integration << "Relevant objects:   runtime_jet_recombination   collects partons that enter jet recombination." << endl;
  for (int i = 0; i < runtime_jet_recombination.size(); i++){
    out_integration << "runtime_jet_recombination[" << setw(3) << i << "] = " << setw(3) << runtime_jet_recombination[i] << "[" << setw(10) <<  left << esi.object_list_selection[runtime_jet_recombination[i]] << "]" << endl;
  }
  out_integration << endl;
  
  out_integration << "Relevant objects:   runtime_photon_isolation   collects partons that could become isolated photons a la Frixione." << endl;
  for (int i = 0; i < runtime_photon_isolation.size(); i++){
    out_integration << "runtime_photon_isolation[" << setw(3) << i << "] = " << setw(3) << runtime_photon_isolation[i] << "[" << setw(10) <<   left << esi.object_list_selection[runtime_photon_isolation[i]] << "]" << endl;
  }
  out_integration << endl;

  /*
  //  runtime_photon_recombination is never used !!!
  out_integration << "Relevant objects:   runtime_photon_recombination   collects partons that enter photon recombination." << endl;
  for (int i = 0; i < runtime_photon_recombination.size(); i++){
    out_integration << "runtime_photon_recombination[" << setw(3) << i << "] = " << setw(3) << runtime_photon_recombination[i] << "[" << setw(10) << left << esi.object_list_selection[runtime_photon_recombination[i]] << "]" << endl;
  }
  out_integration << endl;
  */
  /*
  out_integration << "Relevant objects:   runtime_missing   collects partons that enter event selection with their parton-level momenta." << endl;
  for (int i = 0; i < runtime_missing.size(); i++){
    out_integration << "runtime_missing[" << setw(3) << i << "] = " << setw(3) << runtime_missing[i] << "[" << setw(10) << left << esi.object_list_selection[runtime_missing[i]] << "]" << endl;
  }
  out_integration << endl;
  */

  out_integration << "Relevant objects:   runtime_original   collects partons that enter event selection with their parton-level momenta." << endl;
  for (int i = 0; i < runtime_original.size(); i++){
    out_integration << "runtime_original[" << setw(3) << i << "] = " << setw(3) << runtime_original[i] << "[" << setw(10) << left << esi.object_list_selection[runtime_original[i]] << "]" << endl;
  }
  out_integration << endl;
  
  out_integration << "Relevant objects:   runtime_order   collects partons that belong to the respective object class." << endl;

  for (int i = 0; i < esi.object_list_selection.size(); i++){
    stringstream sstemp;
    for (int j = 0; j < runtime_order[i].size(); j++){sstemp << setw(3) << runtime_order[i][j] << "[" << setw(6) << csi->type_parton[0][runtime_order[i][j]] << "]";}// << setw(3) << runtime_order[i][csi->type_parton[i]] << " ["
    out_integration << setw(10) << esi.object_list_selection[i] << " -> " << sstemp.str() << endl;
  }
  out_integration << endl;

  out_integration << "Relevant objects:   runtime_order_inverse   collects object classes a parton belongs to." << endl;

  // !!! csi->type_parton has n_ps entries (vector<vector<int> >) !!!
  for (int i = 3; i < csi->type_parton[0].size(); i++){
    stringstream sstemp;
    for (int j = 0; j < runtime_order_inverse[i].size(); j++){
      sstemp << setw(10) << esi.object_list_selection[runtime_order_inverse[i][j]] << " [" << setw(3) << runtime_order_inverse[i][j] << "]";
    }
    out_integration << "pa[" << setw(3) << i << "] = " << setw(5) << csi->type_parton[0][i] << " -> " << sstemp.str() << endl;
  }
  out_integration << endl;
  
  out_integration << endl;
  /*
  //  old version:
  out_integration << "   generated (nan/tech)         cut    counted    Xsection contribution +- absolute error        dev/XScont   dev/XS@norm" << endl;
  out_integration << "*************************************************************************************************************************" << endl;
  */
  out_integration << "   generated [nan]    accepted (techcut)    rejected    Xsection contribution +- absolute error        dev/XScont   dev/XS@norm" << endl;
  out_integration << "*******************************************************************************************************************************" << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::output_check_value_scale(){
  Logger logger("observable_set::output_check_value_scale");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  /*
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
  */
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::output_check_running_alpha_S(){
  Logger logger("observable_set::output_check_running_alpha_S");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
 
  ofstream out_check_alpha_S;
  string filename_check_alpha_S = "alphaS." + LHAPDFname + ".dat";
  out_check_alpha_S.open(filename_check_alpha_S.c_str(), ofstream::out | ofstream::trunc);
  out_check_alpha_S << "# " << left << setw(40) << LHAPDFname << "getNf = " << LHAPDF::getNf() << "     N_nondecoupled = " << N_nondecoupled << endl;
  int decade = 3;
  int step = 1000;
  for (int i_d = 0; i_d < decade; i_d++){
    int temp_step = step;
    if (i_d == decade - 1){temp_step++;}
    for (int i_s = 0; i_s < temp_step; i_s++){
      double temp_mu_ren = pow(10., double(i_d) + double(i_s) / step);
      double temp_alpha_S = LHAPDF::alphasPDF(temp_mu_ren);
      out_check_alpha_S << setprecision(15) << setw(23) << temp_mu_ren << setprecision(15) << setw(23) << temp_alpha_S << endl;
    }
  }
  out_check_alpha_S.close();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::output_finalization_integration(phasespace_set & psi){
  Logger logger("observable_set::output_finalization_integration");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_INFO << "evaluation time = " << h << "h " << min << "min " << sec << "sec     n_gen = " << psi_i_gen << "   n_acc = " << psi_i_acc << endl;
  logger << LOG_INFO << "evaluation time per generated event = "  << showpoint << setprecision(4) << double(3600 * h + 60 * min + sec) / psi_i_gen << " / sec" << endl;
  logger << LOG_INFO << "evaluation time per accepted event  = "  << showpoint << setprecision(4) << double(3600 * h + 60 * min + sec) / psi_i_acc << " / sec" << endl;

  ofstream out_integration;
  out_integration.open(filename_integration.c_str(), ofstream::out | ofstream::app);
  static string stars = "***********************************************";
  out_integration << stars << " final result " << stars << endl;
  out_integration.close();

  ofstream out_execution;
  out_execution.open(filename_execution.c_str(), ofstream::out | ofstream::trunc);
  out_execution << "final result" << endl;
  out_execution.close();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

// access to the already calculated 'value' scales:
/*
for (int i_a = 0; i_a < osi_n_ps; i_a++){} // loop over all phasespaces -> i_a
for (int i_s = 0; i_s < osi_n_set_TSV; i_s++){} // // loop over all scaleset -> i_s
for (int i_m = 0; i_m < osi_n_scale_dyn_ren[osi_dynamic_scale_ren_TSV[i_s]].size(); i_m++){} // loop over ren. scales for scaleset (i_s) -> i_m 
for (int i_m = 0; i_m < osi_n_scale_dyn_fact[osi_dynamic_scale_fact_TSV[i_s]].size(); i_m++){} // loop over fact. scales for scaleset (i_s) -> i_m 

osi_value_relative_scale_ren[osi_dynamic_scale_ren_TSV[i_s]][i_m]
osi_value_relative_scale2_ren[osi_dynamic_scale_ren_TSV[i_s]][i_m]
osi_value_relative_factor_alpha_S[i_a][osi_dynamic_scale_ren_TSV[i_s]][i_m]
osi_value_scale_ren[i_a][osi_dynamic_scale_ren_TSV[i_s]][i_m]
osi_value_scale2_ren[i_a][osi_dynamic_scale_ren_TSV[i_s]][i_m]

osi_value_relative_scale_fact[osi_dynamic_scale_fact_TSV[i_s]][i_m]
osi_value_relative_scale2_fact[osi_dynamic_scale_fact_TSV[i_s]][i_m]
osi_value_relative_logscale2_fact[osi_dynamic_scale_fact_TSV[i_s]][i_m]
osi_value_scale_fact[i_a][osi_dynamic_scale_fact_TSV[i_s]][i_m]
osi_value_scale2_fact[i_a][osi_dynamic_scale_fact_TSV[i_s]][i_m]
*/


/*
for (int i_s = 0; i_s < osi_n_set_TSV; i_s++){
  int temp_DS = osi_dynamic_scale_ren_TSV[i_s];
  for (int i_r = 0; i_r < osi_n_scale_ren_TSV[i_s]; i_r++){
    if (i_s == 1 && (i_r == 2 || i_r == 4 || i_r == 6)){
      logger << LOG_DEBUG_VERBOSE << right << setw(40) << "oset_value_scale_ren[0][" << temp_DS << "][" << osi_no_value_ren_TSV[i_s][i_r] << "] = " << setw(24) << setprecision(15) << (*oset.value_scale_ren)[0][temp_DS][osi_no_value_ren_TSV[i_s][i_r]] << "   =   " << setw(24) << setprecision(15) << osi_value_scale_ren[0][temp_DS][osi_no_value_ren_TSV[i_s][i_r]] << endl;
      logger << LOG_DEBUG_VERBOSE << right << setw(40) << "oset_pointer_scale_ren[0][" << i_s << "][i_r] = " <<  setw(24) << setprecision(15) << *(*oset.pointer_scale_ren)[0][i_s][i_r] << "   =   " << setw(24) << setprecision(15) << *osi_pointer_scale_ren[0][i_s][i_r] << endl;
    }
  }
 }
*/

void observable_set::output_integrand_maximum_psp(phasespace_set & psi){
  Logger logger("observable_set::output_integrand_maximum_psp");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_maxevent;
  out_maxevent.open(filename_maxevent.c_str(), ofstream::out | ofstream::app);  
  
  for (int i_p = 0; i_p < p_parton[0].size(); i_p++){
    out_maxevent << "p_parton[0][" << i_p << "] = " << p_parton[0][i_p] << "   pT = " << p_parton[0][i_p].pT() << endl;
  }
  out_maxevent << endl;
  /*
  // write out all relevant pairings to spot badly mapped resonances
  if (csi->n_particle > 1){
    out_maxevent << "Possible cancidates for 2-particle resonances:" << endl;
    for (int i_p = 3; i_p < p_parton[0].size(); i_p++){
      for (int j_p = i_p + 1; j_p < p_parton[0].size(); j_p++){
	out_maxevent << "p_parton[" << i_p << " + " << j_p << "] = " << p_parton[0][i_p] + p_parton[0][j_p] << ", sqrt(s)= "  << (p_parton[0][i_p] + p_parton[0][j_p]).m() << endl;
      }
    }
    out_maxevent << endl;
  }
  if (csi->n_particle > 2){
    out_maxevent << "Possible cancidates for 3-particle resonances:" << endl;
    for (int i_p = 3; i_p < p_parton[0].size(); i_p++){
      for (int j_p = i_p + 1; j_p < p_parton[0].size(); j_p++){
	for (int k_p = j_p + 1; k_p < p_parton[0].size(); k_p++){
	  out_maxevent << "p_parton[" << i_p << " + " << j_p << " + " << k_p << "] = " << p_parton[0][i_p] + p_parton[0][j_p] + p_parton[0][k_p] << ", sqrt(s)= "  << (p_parton[0][i_p] + p_parton[0][j_p] + p_parton[0][k_p]).m() << endl;
	}
      }
    }
    out_maxevent << endl;
  }
  if (csi->n_particle > 3){
    out_maxevent << "Possible cancidates for 4-particle resonances:" << endl;
    for (int i_p = 3; i_p < p_parton[0].size(); i_p++){
      for (int j_p = i_p + 1; j_p < p_parton[0].size(); j_p++){
	for (int k_p = j_p + 1; k_p < p_parton[0].size(); k_p++){
	  for (int l_p = k_p + 1; l_p < p_parton[0].size(); l_p++){
	    out_maxevent << "p_parton[0][" << i_p << " + " << j_p << " + " << k_p << " + " << l_p << "] = " << p_parton[0][i_p] + p_parton[0][j_p] + p_parton[0][k_p] + p_parton[0][l_p] << ", sqrt(s)= "  << (p_parton[0][i_p] + p_parton[0][j_p] + p_parton[0][k_p] + p_parton[0][l_p]).m() << endl;
	  }
	}
      }
    }
    out_maxevent << endl;
  }
  */
  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i_u = 2; i_u < p_parton[i_a].size() - 2; i_u ++){
      //    for (int i_u = 2; i_u < csi->n_particle + 1; i_u ++){
      out_maxevent << "Possible cancidates for " << i_u << "-particle resonances for phase-space " << i_a << ":" << endl;
      for (int i_x = 4; i_x < psi_xbp_all[i_a].size(); i_x += 4){
	vector<int> xbno_contained_particle = vectorint_from_binary(i_x);
	//      if (xbno_contained_particle.size() > 1){
	if (xbno_contained_particle.size() == i_u){
	  stringstream tempss;
	  tempss << setw(3) << i_x << " =^= ";
	  for (int i_xi = 0; i_xi < xbno_contained_particle.size(); i_xi++){
	    tempss << setw(3) << xbno_contained_particle[i_xi];
	    if (i_xi < xbno_contained_particle.size() - 1){tempss << " + ";}
	  }
	  if (psi_xbp_all[i_a][i_x] != nullvector){out_maxevent << "sqrt(s[" << tempss.str() << "]) = " << psi_xbp_all[i_a][i_x].m() << endl;}
	}
      }
      out_maxevent << endl;
    }
  }

  out_maxevent << "sqrt(s_part)    = " << (p_parton[0][1] + p_parton[0][2]).m() << endl;
  out_maxevent << "tau             = " << psi_x_pdf[0] << endl;
  out_maxevent << "g_tot           = " << setw(25) << setprecision(16) << psi_g_tot << endl;
  out_maxevent << endl;
  out_maxevent << "psi.MC_tau.channel  = " << setw(5) << psi_MC_tau.channel << "   (absolute probability = " << setw(25) << setprecision(16) << psi_MC_tau.alpha[psi_MC_tau.channel] << endl;
  ///  out_maxevent << "tau_MC_channel  = " << setw(5) << psi_tau_MC_channel << "   (absolute probability = " << setw(25) << setprecision(16) << psi_tau_MC_alpha[psi_tau_MC_channel] << endl;
  //psi_tau_MC_tau_gamma[0][psi_tau_channel] << endl;
  out_maxevent << "psi.MC_phasespace.channel      = " << setw(5) << psi.MC_phasespace.channel << "   (absolute probability = " << setw(25) << setprecision(16) << psi.MC_phasespace.alpha[psi.MC_phasespace.channel] << endl;

  logger << LOG_DEBUG_VERBOSE << "psi.IS_tau.active_optimization = " << psi.IS_tau.active_optimization << endl;
  if (psi.IS_tau.active_optimization != -1){
    logger << LOG_DEBUG_VERBOSE << "psi.IS_tau.channel = " << psi.IS_tau.channel << endl;
    out_maxevent << "tau_channel     = " << setw(5) << psi.IS_tau.channel << "   (absolute probability = " << setw(25) << setprecision(16) << psi.IS_tau.alpha[psi.IS_tau.channel] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "psi.IS_x1x2.active_optimization = " << psi.IS_x1x2.active_optimization << endl;
  if (psi.IS_x1x2.active_optimization != -1){
    logger << LOG_DEBUG_VERBOSE << "psi.IS_x1x2.channel = " << psi.IS_x1x2.channel << endl;
    out_maxevent << "x1x2_channel    = " << setw(5) << psi.IS_x1x2.channel << "   (absolute probability = " << setw(25) << setprecision(16) << psi.IS_x1x2.alpha[psi.IS_x1x2.channel] << endl;
  }
  /*
  out_maxevent << "tau_channel     = " << setw(5) << psi_tau_channel << "   (absolute probability = " << setw(25) << setprecision(16) << psi_tau_alpha[psi_tau_channel] << endl;
  out_maxevent << "x1x2_channel    = " << setw(5) << psi_x1x2_channel << "   (absolute probability = " << setw(25) << setprecision(16) << psi_x1x2_alpha[psi_x1x2_channel] << endl;
  */  /*
    out_maxevent << "psi.MC_phasespace.channel   = " << setw(4) << psi.MC_phasespace.channel << endl;
    out_maxevent << "psi_tau_channel  = " << setw(4) << psi_tau_channel << "   psi_tau_alpha[" << setw(4) << psi_tau_channel << "] = " << psi_tau_alpha[psi_tau_channel] << endl;
    out_maxevent << "psi_x1x2_channel = " << setw(4) << psi_x1x2_channel << "  psi_x1x2_alpha[" << setw(4) << psi_x1x2_channel << "] = " << psi_x1x2_alpha[psi_x1x2_channel] << endl;
    */
    /*
      for (int ig = 0; ig < MC_phasespace.alpha.size(); ig++){
      out_maxevent << "dMC_alpha[" << setw(3) << ig + 1 << "] = " << setw(25) << setprecision(16) << MC_phasespace.alpha[ig] << "   dg[" << setw(3) << ig + 1 << "] = " << setw(25) << setprecision(16) << g[ig] << endl;
      }
    */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




void observable_set::output_testpoint_VA_ioperator(ofstream & out_comparison){
  Logger logger("observable_set::output_testpoint_VA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  out_comparison << endl;
  out_comparison << setw(12) << " SK_ME2_I = " << setprecision(15) << setw(23) << VA_I_ME2 << endl;
  out_comparison << endl;
  for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
    for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
      out_comparison << left << setw(26) << (*VA_ioperator)[i_a][j_a].name() << ": " << "I_ME2_cf[i_a = " << setw(2) << right << i_a << "][j_a = " << setw(2) << right << j_a << "] = " << setw(23) << setprecision(15) << VA_I_ME2_cf[i_a][j_a] << "   ME2_cf[" << setw(2) << right << i_a << "][" << setw(2) << right << j_a << "] = " << setw(23) << setprecision(15) << VA_ME2_cf[i_a][j_a] << endl;
    }
  }
  out_comparison << endl;
  out_comparison << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_testpoint_VA_result(ofstream & out_comparison){
  Logger logger("observable_set::output_testpoint_VA_result");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  out_comparison << setw(12) << " ME2_born = " << setprecision(15) << setw(23) << VA_b_ME2 << endl;
  out_comparison << setw(12) << "ME2_V+X+I = " << setprecision(15) << setw(23) << (VA_V_ME2 + VA_X_ME2 + VA_I_ME2) << endl;
  out_comparison << setw(12) << "  ME2_V+X = " << setprecision(15) << setw(23) << (VA_V_ME2 + VA_X_ME2) << endl;
  out_comparison << endl;
  out_comparison << setw(12) << "    ME2_V = " << setprecision(15) << setw(23) << VA_V_ME2 << endl;
  out_comparison << setw(12) << "    ME2_X = " << setprecision(15) << setw(23) << VA_X_ME2 << endl;
  out_comparison << setw(12) << "    ME2_I = " << setprecision(15) << setw(23) << VA_I_ME2 << endl;
  out_comparison << endl;
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_testpoint_VA_Delta(ofstream & out_comparison, int i, double & i_Delta, string & s_Delta){
  Logger logger("observable_set::output_testpoint_VA_Delta");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (i == 0){out_comparison << setw(23) << s_Delta << setw(23) << "ME2_V+X+I" << setw(23) << "ME2_V" << setw(23) << "ME2_X" << setw(23) << "ME2_I" << endl;}
  out_comparison << setprecision(5) << noshowpoint << setw(23) << i_Delta << setprecision(15) << setw(23) << showpoint << VA_V_ME2 + VA_X_ME2 + VA_I_ME2 << setprecision(15) << setw(23) << showpoint << VA_V_ME2 << setprecision(15) << setw(23) << showpoint << VA_X_ME2 << setprecision(15) << setw(23) << showpoint << VA_I_ME2 << endl;
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_testpoint_CA(phasespace_set & psi){
  Logger logger("observable_set::output_testpoint_CA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if ((*CA_collinear).size() == 0){return;}

  ofstream out_comparison;
  out_comparison.open(filename_comparison.c_str(), ofstream::out | ofstream::app);
  logger << LOG_DEBUG << "(*CA_collinear).size() = " << (*CA_collinear).size() << endl;

  int temp_symmetry_factor = 1;
  for (int i_p = 1; i_p < 3; i_p++){
    if (csi->type_parton[0][i_p] == 0){temp_symmetry_factor *= 2 * 8;}
    else if (abs(csi->type_parton[0][i_p]) < 10){temp_symmetry_factor *= 2 * 3;}
    else if (csi->type_parton[0][i_p] == 22){temp_symmetry_factor *= 2 * 1;}
  }

  out_comparison << "Colour-correlated matrix elements evaluated at alpha_S(" << scale_ren << "): " << endl << endl;

  logger << LOG_DEBUG << "(*CA_collinear).size() = " << (*CA_collinear).size() << endl;
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    logger << LOG_DEBUG << "(*CA_collinear)[" << i_c << "].size() = " << (*CA_collinear)[i_c].size() << endl;
    for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
      out_comparison << "Collinear " << setw(2) << right << i_c << ":   ";
      if (j_c == 0){
	out_comparison << left << setw(26) << (*CA_collinear)[i_c][j_c].name() << ": " << "ME2_cf[i_c = " << setw(2) << right << i_c << "][j_c = " << setw(2) << right << j_c << "] = " << setw(23) << setprecision(15) << CA_ME2_cf[i_c][j_c] / (*CA_collinear)[i_c][j_c].charge_factor() << endl;
      }
      else {
	out_comparison << left << setw(26) << (*CA_collinear)[i_c][j_c].name() << ": " << "ME2_cf[i_c = " << setw(2) << right << i_c << "][j_c = " << setw(2) << right << j_c << "] = " << setw(23) << setprecision(15) << CA_ME2_cf[i_c][j_c] << "   ratio = " << CA_ME2_cf[i_c][j_c] / (CA_ME2_cf[i_c][0] / (*CA_collinear)[i_c][0].charge_factor()) << endl;
      }

      // in this version, CA_ME2_cf[i_c][0] contains a factor Q²_f !!!
      //      out_comparison << left << setw(26) << (*CA_collinear)[i_c][j_c].name() << ": " << "ME2_cf[i_c = " << setw(2) << right << i_c << "][j_c = " << setw(2) << right << j_c << "] = " << setw(23) << setprecision(15) << CA_ME2_cf[i_c][j_c] << "   ratio = " << CA_ME2_cf[i_c][j_c] / CA_ME2_cf[i_c][0] << endl;
      //      out_comparison << left << setw(26) << (*CA_collinear)[i_c][j_c].name() << ": " << "ME2_cf[i_c = " << setw(2) << right << i_c << "][j_c = " << setw(2) << right << j_c << "] * 36 = " << setw(23) << setprecision(15) << CA_ME2_cf[i_c][j_c] * 36. << endl;
    }
    out_comparison << endl;

    out_comparison << "Colour-correlated matrix elements including correct (relative) factors of powers of alpha_S(mu_R): " << endl << endl;

    double rel_born_exponent = 1.;
    if (type_correction == "QCD"){rel_born_exponent = 1. - 1. / double(contribution_order_alpha_s);}
    else if (type_correction == "QEW"){rel_born_exponent = 1.;}

    logger << LOG_DEBUG << "n_set_TSV = " << n_set_TSV << endl;

    for (int i_s = 0; i_s < n_set_TSV; i_s++){
      //      int v_sf = dynamic_scale_fact_TSV[i_s];
      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	//	int v_sr = dynamic_scale_ren_TSV[i_s];
	out_comparison << "alpha_S(" << *(pointer_scale_ren)[0][i_s][i_r] << ") = " << setw(23) << setprecision(15) << LHAPDF::alphasPDF(*(pointer_scale_ren)[0][i_s][i_r]) << endl << endl;
	//	out_comparison << "alpha_S(" << value_scale_ren[0][v_sr][v_xr] << ") = " << setw(23) << setprecision(15) << LHAPDF::alphasPDF(value_scale_ren[0][v_sr][v_xr]) << endl << endl;

	for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
	  stringstream sstemp;
	  sstemp << "mu_R = " << setw(15) << setprecision(8) << *(pointer_scale_ren)[0][i_s][i_r] << "   ";
	  sstemp << left << setw(26) << (*CA_collinear)[i_c][j_c].name() << ": " << "ME2_cf[i_c = " << setw(2) << right << i_c << "][j_c = " << setw(2) << right << j_c << "] = ";
	  if (j_c == 0){
	    sstemp << setw(23) << setprecision(15) << CA_ME2_cf[i_c][j_c] / (*CA_collinear)[i_c][j_c].charge_factor() * pow(*(pointer_relative_factor_alpha_S)[0][i_s][i_r], rel_born_exponent);
	    sstemp << "  , ... * " << temp_symmetry_factor << " = " << setw(23) << setprecision(15) << CA_ME2_cf[i_c][j_c] / (*CA_collinear)[i_c][j_c].charge_factor() * pow(*(pointer_relative_factor_alpha_S)[0][i_s][i_r], rel_born_exponent) * temp_symmetry_factor;
	  }
	  else {
	    sstemp << setw(23) << setprecision(15) << CA_ME2_cf[i_c][j_c] * pow(*(pointer_relative_factor_alpha_S)[0][i_s][i_r], rel_born_exponent);
	    sstemp << "  , ... * " << temp_symmetry_factor << " = " << setw(23) << setprecision(15) << CA_ME2_cf[i_c][j_c] * pow(*(pointer_relative_factor_alpha_S)[0][i_s][i_r], rel_born_exponent) * temp_symmetry_factor;
	  }
	  // << "   ratio = " << CA_ME2_cf[i_c][j_c] / (CA_ME2_cf[i_c][0] / (*CA_collinear)[i_c][0].charge_factor())
	  out_comparison << sstemp.str() << endl;
	  if (*(pointer_relative_factor_alpha_S)[0][i_s][i_r] != 1.){
	    out_comparison << "pow(" << *(pointer_relative_factor_alpha_S)[0][i_s][i_r] << ", " << rel_born_exponent << ") = " << pow(*(pointer_relative_factor_alpha_S)[0][i_s][i_r], rel_born_exponent) << endl;
	  }
	}
	out_comparison << endl;
      }
    }

    stringstream sstempmode;
    if (switch_KP == 0){sstempmode << "KP terms";}
    else if (switch_KP == 1){sstempmode << "P terms";}
    else if (switch_KP == 2){sstempmode << "K terms";}

    out_comparison << sstempmode.str() << " evaluated at alpha_S(" << scale_ren << "): " << endl << endl;

    for (int i_s = 0; i_s < n_set_TSV; i_s++){
      logger << LOG_DEBUG_VERBOSE << "begin: i_s = " << i_s << endl;
      //    for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
      int v_sf = dynamic_scale_fact_TSV[i_s];
      logger << LOG_DEBUG_VERBOSE << "value_logscale2_fact_papi[" << v_sf << "].size() = " << value_logscale2_fact_papi[v_sf].size() << endl;
      logger << LOG_DEBUG_VERBOSE << "value_ME2_KP[" << v_sf << "].size() = " << value_ME2_KP[v_sf].size() << endl;
      //      assert(value_logscale2_fact_papi[v_sf].size() == value_scale_fact[0][v_sf].size());
      for (int v_xf = 0; v_xf < value_scale_fact[0][v_sf].size(); v_xf++){
	for (int i_x = 0; i_x < 3; i_x++){
	  if (value_ME2_KP[v_sf][v_xf][i_x][i_c] != 0.){
	    stringstream sstemp;
	    sstemp << "mu_F = " << setw(15) << setprecision(8) << value_scale_fact[0][v_sf][v_xf] << "   ";
	    if      (i_x == 0){sstemp << left << setw(28) << "regular, regular []_+     : ";}
	    else if (i_x == 1){sstemp << left << setw(28) << "endpoint []_+             : ";}
	    else if (i_x == 2){sstemp << left << setw(28) << "endpoint, integrated []_+ : ";}
	    //	    logger << LOG_DEBUG_VERBOSE << "*(pointer_ME2term)[" << i_c << "][" << i_x << "][" << i_s << "][" << i_r << "].size() = " << (pointer_ME2term)[i_c][i_x][i_s][i_r].size() << endl;  // different meaning of i_r and v_xf for i_s and v_sf in general !!!
	    //	    sstemp << setw(23) << setprecision(15) << *(pointer_ME2term)[i_c][i_x][i_s][i_r][i_f];
	    logger << LOG_DEBUG_VERBOSE << "value_ME2term_fact[" << i_c << "][" << i_x << "][" << v_sf << "].size() = " << value_ME2term_fact[i_c][i_x][v_sf].size() << endl;

	    sstemp << setw(23) << setprecision(15) << value_ME2term_fact[i_c][i_x][v_sf][v_xf];
	    out_comparison << sstemp.str() << endl;
	  }
	logger << LOG_DEBUG_VERBOSE << "end: v_xf = " << v_xf << endl;
	}
	out_comparison << endl;
      }
      logger << LOG_DEBUG_VERBOSE << "end: i_s = " << i_s << endl;
    }


    
    out_comparison << sstempmode.str() << " including correct (relative) factors of powers of alpha_S(mu_R): " << endl << endl;
    
    //    for (int v_sf = 0; v_sf < max_dyn_fact + 1; v_sf++){
    for (int i_s = 0; i_s < n_set_TSV; i_s++){
      int v_sf = dynamic_scale_fact_TSV[i_s];
      int v_sr = dynamic_scale_ren_TSV[i_s];
      for (int v_xr = 0; v_xr < value_scale_ren[0][v_sr].size(); v_xr++){
	//      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	//	assert (value_logscale2_fact_papi[v_sf].size() == value_scale_fact[0][v_sf].size());
	//	logger << LOG_DEBUG << "value_logscale2_fact_papi[" << v_sf << "].size() = " << value_logscale2_fact_papi[v_sf].size() << " =?= " << value_scale_fact[0][v_sf].size() << " = value_scale_fact[0][" << v_sf << "].size()" << endl;
	for (int v_xf = 0; v_xf < value_scale_fact[0][v_sf].size(); v_xf++){
	  for (int i_x = 0; i_x < 3; i_x++){
	    if (value_ME2_KP[v_sf][v_xf][i_x][i_c] != 0.){
	      stringstream sstemp;
	      sstemp << "v_sf = " << setw(2) << v_sf << "   " << "v_xf = " << setw(2) << v_xf << "   " << "v_sr = " << setw(2) << v_sr << "   " << "   " << "v_xr = " << setw(2) << v_xr << "   ";
	      sstemp << "mu_R = " << setw(15) << setprecision(8) << value_scale_ren[0][v_sr][v_xr] << "   ";
	      sstemp << "mu_F = " << setw(15) << setprecision(8) << value_scale_fact[0][v_sf][v_xf] << "   ";
	      if      (i_x == 0){sstemp << left << setw(28) << "regular, regular []_+     : ";}
	      else if (i_x == 1){sstemp << left << setw(28) << "endpoint []_+             : ";}
	      else if (i_x == 2){sstemp << left << setw(28) << "endpoint, integrated []_+ : ";}
	      //	      sstemp << setw(23) << setprecision(15) << value_ME2_KP[v_sf][v_xf][i_x][i_c] * *(pointer_relative_factor_alpha_S)[0][i_s][i_r];
	      //	      sstemp << setw(23) << setprecision(15) << *(pointer_ME2term)[i_c][i_x][i_s][i_r][i_f] * *(pointer_relative_factor_alpha_S)[0][i_s][i_r];
	      //	      sstemp << setw(23) << setprecision(15) << value_ME2term_fact[i_c][i_x][v_sf][v_xf] * *(pointer_relative_factor_alpha_S)[0][i_s][i_r];
	      sstemp << setw(23) << setprecision(15) << value_ME2term_fact[i_c][i_x][v_sf][v_xf] * value_relative_factor_alpha_S[0][v_sr][v_xr];

	      logger << LOG_DEBUG << sstemp.str() << endl;
      
	      out_comparison << sstemp.str() << endl;
	    }
	  }
	  out_comparison << endl;
	}
      }
    }
    out_comparison << endl << endl;
  }
  logger << LOG_DEBUG << "after (*CA_collinear).size() loop " << endl;

  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
 
    }
  }

  vector<double> g_z_coll(3.,0);
  for (int i_z = 1; i_z < 3; i_z++){g_z_coll[i_z] = 1. / (1. - x_pdf[i_z]);}



  vector<vector<double> > all_xz_coll_pdf(3, vector<double> (3));
  for (int i_x = 0; i_x < 3; i_x++){all_xz_coll_pdf[i_x][0] = x_pdf[i_x];}
  for (int i_z = 1; i_z < 3; i_z++){
    all_xz_coll_pdf[0][i_z] = z_coll[i_z];
    all_xz_coll_pdf[i_z][i_z] = x_pdf[i_z] / z_coll[i_z];
    all_xz_coll_pdf[i_z % 2 + 1][i_z] = x_pdf[0] / z_coll[i_z];
  }

  /*
  cout << "after all_xz_coll_pdf" << endl;
  for (int i_x = 0; i_x < 3; i_x++){
    for (int i_z = 0; i_z < 3; i_z++){
      cout << "all_xz_coll_pdf[" << i_x << "][" << i_z << "] = " << all_xz_coll_pdf[i_x][i_z] << endl;
    }
  }
  */
  //  if (switch_TSV){calculate_pdf_LHAPDF_CA_collinear_TSV((*CA_collinear), all_xz_coll_pdf, this);}

  if (switch_TSV){

    calculate_pdf_LHAPDF_CA_collinear_TSV(all_xz_coll_pdf);

    /*
    logger << LOG_DEBUG << "CHECK" << endl;
    for (int i_c = 0; i_c < n_pc; i_c++){
      for (int i_z = 0; i_z < n_pz; i_z++){
	for (int i_s = 0; i_s < n_set_TSV; i_s++){
	  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	    for (int i_h = 0; i_h < 3; i_h++){
	      for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
		logger << LOG_DEBUG << "value_pdf_factor_combination_1[" << i_c << "][" << i_z << "][" << dynamic_scale_fact_TSV[i_s] << "][" << no_value_fact_TSV[i_s][i_f] << "][" << i_h << "][" << i_i << "] = " << value_pdf_factor_combination_1[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h][i_i] << endl;
		logger << LOG_DEBUG << "value_pdf_factor_combination_2[" << i_c << "][" << i_z << "][" << dynamic_scale_fact_TSV[i_s] << "][" << no_value_fact_TSV[i_s][i_f] << "][" << i_h << "][" << i_i << "] = " << value_pdf_factor_combination_2[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h][i_i] << endl;
	      }
	    }
	  }
	}
      }
    }
    */

    out_comparison << endl << endl << endl; 
    out_comparison << "Comparison to Sherpa K+P terms" << endl << endl;

    static double alpha_S_2pi = alpha_S * inv2pi;
    static double alpha_e_2pi = msi.alpha_e * inv2pi;
    double alpha_x_2pi = 0.;
    if ((*CA_collinear)[0][0].type_correction() == 1){alpha_x_2pi = alpha_S_2pi;}
    else if ((*CA_collinear)[0][0].type_correction() == 2){alpha_x_2pi = alpha_e_2pi;}
    
    
    for (int i_s = 0; i_s < n_set_TSV; i_s++){
      int v_sr = dynamic_scale_ren_TSV[i_s];
      int v_sf = dynamic_scale_fact_TSV[i_s];
      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	int v_xr = no_value_ren_TSV[i_s][i_r];
	for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	  int v_xf = no_value_fact_TSV[i_s][i_f];
	  /////	for (int v_xf = 0; v_xf < value_scale_fact[0][v_sf].size(); v_xf++){
	  vector<double> temp_Sherpa_KPterms_pdf((*CA_collinear).size(), 0.);
	  vector<double> temp_Sherpa_KPterms((*CA_collinear).size(), 0.);
	  vector<double> temp_Sherpa_KPterms_delta((*CA_collinear).size(), 0.);
	  vector<double> temp_Sherpa_KPterms_endpoint((*CA_collinear).size(), 0.);
	  vector<double> temp_Sherpa_KPterms_regular((*CA_collinear).size(), 0.);
	  vector<double> temp_sum_c((*CA_collinear).size(), 0.);
	  for (int i_c = 0; i_c < n_pc; i_c++){
	    
	    //	out_comparison << endl << "value_scale_ren[" << 0 << "][" << v_sr << "][" << i_r << "] = " << value_scale_ren[0][v_sr][i_r] << endl;
	    double temp_c = 0.;
	    double temp_c_z = 0.;
	    for (int i_z = 0; i_z < n_pz; i_z++){
	      for (int i_h = 1; i_h < 3; i_h++){
		if (i_z == 0 && (*CA_collinear)[i_c][0].in_collinear()[i_z] == 1){
		  temp_c += *(pointer_pdf_factor)[i_c][i_z][i_s][i_f][i_h] // added again for different usage
		    * (  *(pointer_ME2term)[i_c][1][i_s][i_r][i_f] / g_z_coll[(*CA_collinear)[i_c][0].no_emitter()] 
			 + *(pointer_ME2term)[i_c][2][i_s][i_r][i_f])
		    ;
		}
		else if ((*CA_collinear)[i_c][0].in_collinear()[i_z] == 1){ // always i_z == dipole_phasespace[i_c][4] !!!
		  temp_c_z += *(pointer_pdf_factor)[i_c][i_z][i_s][i_f][i_h]  // added again for different usage
		    * *(pointer_ME2term)[i_c][0][i_s][i_r][i_f] / g_z_coll[i_z]
		    ;
		}
	      }
	    }
	    temp_sum_c[i_c] = (temp_c + temp_c_z) * *(pointer_relative_factor_alpha_S)[0][i_s][i_r];
	    
	    stringstream sstemp;
	    sstemp << "Collinear " << setw(2) << right << i_c << ":   " << left << setw(26) << (*CA_collinear)[i_c][0].name();
	    sstemp << "mu_R = " << setw(15) << setprecision(8) << value_scale_ren[0][v_sr][v_xr] << "   ";
	    sstemp << "mu_F = " << setw(15) << setprecision(8) << value_scale_fact[0][v_sf][v_xf] << "   ";
	    sstemp << endl << endl;
	    int tablevel0 = 0;
	    int tablevel1 = tablevel0 + 56;
	    int tablevel2 = tablevel1 + 50;
	    
	    for (int i_z = 0; i_z < 3; i_z++){
	      if (i_z == 0){
		
		// endpoint []_+(1)   and   delta   (including []_+(int))
		
		stringstream name_Sherpa_1;
		name_Sherpa_1 << "kpc";
		if ((*CA_collinear)[i_c][0].no_emitter() == 1){name_Sherpa_1 << "a";}
		if ((*CA_collinear)[i_c][0].no_emitter() == 2){name_Sherpa_1 << "b";}
		if ((*CA_collinear)[i_c][0].type() == 0 || (*CA_collinear)[i_c][0].type() == 3){name_Sherpa_1 << "[2]";}
		if ((*CA_collinear)[i_c][0].type() == 1 || (*CA_collinear)[i_c][0].type() == 2){name_Sherpa_1 << "[0]";}
		
		stringstream name_Sherpa_pdf_1;
		name_Sherpa_pdf_1 << "f";
		if ((*CA_collinear)[i_c][0].no_emitter() == 1){name_Sherpa_pdf_1 << "a";}
		if ((*CA_collinear)[i_c][0].no_emitter() == 2){name_Sherpa_pdf_1 << "b";}
		if ((*CA_collinear)[i_c][0].type() == 0 || (*CA_collinear)[i_c][0].type() == 3){name_Sherpa_pdf_1 << "g";}
		if ((*CA_collinear)[i_c][0].type() == 1 || (*CA_collinear)[i_c][0].type() == 2){name_Sherpa_pdf_1 << "q";}
		name_Sherpa_pdf_1 << " * f";
		if ((*CA_collinear)[i_c][0].no_emitter() == 1){name_Sherpa_pdf_1 << "b";}
		if ((*CA_collinear)[i_c][0].no_emitter() == 2){name_Sherpa_pdf_1 << "a";}
		
		//  m_kpca[0]*faq
		//  m_kpca[2]*fag
		//  m_kpcb[0]*fbq
		//  m_kpcb[2]*fbg
		
		sstemp << left;
		
		double temp_Munich_KP_endpoint = 0.;
		double temp_Munich_kpc_endpoint = 0.;
		double temp_Munich_KP_delta = 0.;	  
		double temp_Munich_kpc_delta = 0.;
		
		if (*(pointer_ME2term)[i_c][1][i_s][i_r][i_f] != 0){
		  
		  // endpoint []_+
		  
		  temp_Munich_KP_endpoint = *(pointer_ME2term)[i_c][1][i_s][i_r][i_f] 
		    * *(pointer_relative_factor_alpha_S)[0][i_s][i_r] 
		    //		  		  * z_coll[(*CA_collinear)[i_c][0].no_emitter()] // !!! factor shifted between different contributions !!!
		    * temp_symmetry_factor // !!! factor shifted between different contributions !!!
		    ;
		  temp_Munich_kpc_endpoint = temp_Munich_KP_endpoint / g_z_coll[(*CA_collinear)[i_c][0].no_emitter()]//;
		    / z_coll[(*CA_collinear)[i_c][0].no_emitter()]; // !!! check if correct here !!! related to g_z_coll !!! not here to compare with Sherpa !!
		  temp_Sherpa_KPterms_endpoint[i_c] += temp_Munich_kpc_endpoint;
		  ///		  temp_Sherpa_KPterms[i_c] += temp_Munich_kpc_endpoint;
		  sstemp << setw(tablevel0) << "" << left << name_Sherpa_1.str() << setw(20) << "_([]_+(1))" << " = " << setw(23) << setprecision(15) << temp_Munich_kpc_endpoint 
			 << " = " << "K+P_factor / g_z_coll" << endl;
		  //		       << " = " << "K+P_factor * pdf_factor / g_z_coll" << endl;
		  
		  // K+P_factor   []_+(1)
		  
		  string sum_KPterm = "";
		  stringstream split_KPterm;
		  int counter = 0;
		  double check_sum_KPterm = 0.;
		  for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
		    string name_Kterm;
		    string name_Pterm;
		    if ((*CA_collinear)[i_c][j_c].no_spectator() == 0){name_Kterm = "Kbar";}
		    else if ((*CA_collinear)[i_c][j_c].no_spectator() < 3){
		      stringstream temp;
		      temp << (*CA_collinear)[i_c][j_c].no_spectator();
		      name_Kterm = "Kt(" + temp.str() + ")";
		      name_Pterm = "P(" + temp.str() + ")";
		    }
		    else if ((*CA_collinear)[i_c][j_c].no_spectator() > 2){
		      stringstream temp;
		      temp << (*CA_collinear)[i_c][j_c].no_spectator();
		      name_Kterm = "gamma(" + temp.str() + ")";
		      name_Pterm = "P(" + temp.str() + ")";
		    }
		    if (data_K[1][i_c][j_c] != 0.){
		      if (counter > 0){sum_KPterm = sum_KPterm + " + ";}
		      counter++;
		      double temp_Sherpa_K_endpoint = 
			alpha_x_2pi 
			* data_K[1][i_c][j_c] 
			/ z_coll[(*CA_collinear)[i_c][0].no_emitter()] // !!! check if correct here !!! related to g_z_coll !!! not here to compare with Sherpa !!
			* *(pointer_relative_factor_alpha_S)[0][i_s][i_r] 
			//		      * z_coll[(*CA_collinear)[i_c][0].no_emitter()] // !!! factor shifted between different contributions (contained in pdf_1/2) !!!
			* temp_symmetry_factor // !!! factor shifted between different contributions !!!
			;
		      check_sum_KPterm += temp_Sherpa_K_endpoint;
		      sum_KPterm = sum_KPterm + name_Kterm;
		      split_KPterm << setw(tablevel2) << "" << left << setw(12) << name_Kterm << " = " << setw(23) << setprecision(15) << temp_Sherpa_K_endpoint << endl;
		    }
		    if (value_data_P[v_sf][v_xf][1][i_c][j_c] != 0.){
		      if (counter > 0){sum_KPterm = sum_KPterm + " + ";}
		      counter++;
		      double temp_Sherpa_P_endpoint = 
			alpha_x_2pi 
			* value_data_P[v_sf][v_xf][1][i_c][j_c] 
			/ z_coll[(*CA_collinear)[i_c][0].no_emitter()] // !!! check if correct here !!! related to g_z_coll !!! not here to compare with Sherpa !!
			* *(pointer_relative_factor_alpha_S)[0][i_s][i_r] 
			//		      * z_coll[(*CA_collinear)[i_c][0].no_emitter()] // !!! factor shifted between different contributions (contained in pdf_1/2) !!!
			* temp_symmetry_factor // !!! factor shifted between different contributions !!!
			;
		      check_sum_KPterm += temp_Sherpa_P_endpoint;
		      sum_KPterm = sum_KPterm + name_Pterm;
		      split_KPterm << setw(tablevel2) << "" << left << setw(12) << name_Pterm << " = " << setw(23) << setprecision(15) << temp_Sherpa_P_endpoint << endl;
		    }
		  }
		  sstemp << setw(tablevel1) << "" << setw(12) << "K+P_factor" << " = " << setw(23) << setprecision(15) << temp_Munich_KP_endpoint
			 << setw(9) << "" << " = " << sum_KPterm << endl << split_KPterm.str(); 
		  sstemp << setw(tablevel1) << "" << setw(12) << "   check" << " = " << setw(23) << setprecision(15) << check_sum_KPterm << endl;
		  sstemp << endl;
		  
		  // 1. / g_z_coll   []_+(1)
		  
		  sstemp << setw(tablevel1) << "" << setw(12) << "1/g_z_coll" << " = " << setw(23) << setprecision(15) << 
		    1. / g_z_coll[(*CA_collinear)[i_c][0].no_emitter()] << endl;
		  
		  // check   []_+(1)
		  
		  sstemp << setw(tablevel0) << "" << left << setw(27) << "   check" << " = " << setw(23) << setprecision(15) << 
		    check_sum_KPterm
		    / g_z_coll[(*CA_collinear)[i_c][0].no_emitter()]
		    //		  * temp_pdf_factor
			 << endl;
		  
		  sstemp << endl;
		}
		
		if (*(pointer_ME2term)[i_c][2][i_s][i_r][i_f] != 0){
		  
		  // delta + []_+(int)
		  
		  temp_Munich_KP_delta = *(pointer_ME2term)[i_c][2][i_s][i_r][i_f] 
		    * *(pointer_relative_factor_alpha_S)[0][i_s][i_r] 
		    * temp_symmetry_factor // !!! factor shifted between different contributions !!!
		    ;
		  temp_Munich_kpc_delta = temp_Munich_KP_delta;
		  temp_Sherpa_KPterms_delta[i_c] += temp_Munich_kpc_delta;	  
		  ///		  temp_Sherpa_KPterms[i_c] += temp_Munich_kpc_delta;	  
		  sstemp << setw(tablevel0) << "" << name_Sherpa_1.str() << setw(20) << "_(delta + []_+(int))" << " = " << setw(23) << setprecision(15) << temp_Munich_kpc_delta
			 << " = " << "K+P_factor" << endl;
		  //		       << " = " << "K+P_factor * pdf_factor" << endl;
		  
		  // K+P_factor delta + []_+(int)
		  
		  string sum_KPterm = "";
		  stringstream split_KPterm;
		  int counter = 0;
		  double check_sum_KPterm = 0.;
		  for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
		    string name_Kterm;
		    string name_Pterm;
		    if ((*CA_collinear)[i_c][j_c].no_spectator() == 0.){name_Kterm = "Kbar";} // && data_K[2][i_c][j_c] != 0
		    else if ((*CA_collinear)[i_c][j_c].no_spectator() < 3){
		      stringstream temp;
		      temp << (*CA_collinear)[i_c][j_c].no_spectator();
		      name_Kterm = "Kt(" + temp.str() + ")";
		      name_Pterm = "P(" + temp.str() + ")";
		    }
		    else if ((*CA_collinear)[i_c][j_c].no_spectator() > 2){
		      stringstream temp;
		      temp << (*CA_collinear)[i_c][j_c].no_spectator();
		      name_Kterm = "gamma(" + temp.str() + ")";
		      name_Pterm = "P(" + temp.str() + ")";
		    }
		    if (data_K[2][i_c][j_c] != 0.){
		      if (counter > 0){sum_KPterm = sum_KPterm + " + ";}
		      counter++;
		      double temp_Sherpa_K_delta = 
			alpha_x_2pi 
			* data_K[2][i_c][j_c] 
			* *(pointer_relative_factor_alpha_S)[0][i_s][i_r] 
			* temp_symmetry_factor
			;
		      check_sum_KPterm += temp_Sherpa_K_delta;
		      sum_KPterm = sum_KPterm + name_Kterm;
		      split_KPterm << setw(tablevel2) << "" << left << setw(12) << name_Kterm << " = " << setw(23) << setprecision(15) << temp_Sherpa_K_delta << endl;
		    }
		    if (value_data_P[v_sf][v_xf][2][i_c][j_c] != 0.){
		      if (counter > 0){sum_KPterm = sum_KPterm + " + ";}
		      counter++;
		      double temp_Sherpa_P_delta = 
			alpha_x_2pi 
			* value_data_P[v_sf][v_xf][2][i_c][j_c] 
			* *(pointer_relative_factor_alpha_S)[0][i_s][i_r] 
			* temp_symmetry_factor
			;
		      logger << LOG_DEBUG << "alpha_x_2pi  = " << alpha_x_2pi  << endl;
		      logger << LOG_DEBUG << "value_data_P[" << v_sf << "][" << v_xf << "][2][i_c][j_c]  = " << value_data_P[v_sf][v_xf][2][i_c][j_c]  << endl;
		      logger << LOG_DEBUG << "*(pointer_relative_factor_alpha_S)[0][i_s][i_r]  = " << *(pointer_relative_factor_alpha_S)[0][i_s][i_r]  << endl;
		      logger << LOG_DEBUG << "temp_symmetry_factor = " << temp_symmetry_factor << endl;
		      logger << LOG_DEBUG << "temp_Sherpa_P_delta = " <<  temp_Sherpa_P_delta << endl;
		      check_sum_KPterm += temp_Sherpa_P_delta;
		      sum_KPterm = sum_KPterm + name_Pterm;
		      split_KPterm << setw(tablevel2) << "" << left << setw(12) << name_Pterm << " = " << setw(23) << setprecision(15) << temp_Sherpa_P_delta << endl;
		    }
		  }
		  sstemp << setw(tablevel1) << "" << setw(12) << "K+P_factor" << " = " << setw(23) << setprecision(15) << temp_Munich_KP_delta
			 << " = " << setw(9) << "" << sum_KPterm << endl << split_KPterm.str(); 
		  sstemp << setw(tablevel1) << "" << setw(12) << "   check" << " = " << setw(23) << setprecision(15) << check_sum_KPterm << endl;
		  sstemp << endl;
		  
		  sstemp << setw(tablevel0) << "" << left << setw(27) << "   check" << " = " << setw(23) << setprecision(15) << 
		    check_sum_KPterm // / z_coll[(*CA_collinear)[i_c][0].no_emitter()]
		    //		  * temp_pdf_factor
		    //		  * *(pointer_pdf_factor_1)[i_c][i_z][i_s][i_f][0] * *(pointer_pdf_factor_2)[i_c][i_z][i_s][i_f][0]
			 << endl;
		  
		  sstemp << endl;
		}
		
		if (*(pointer_ME2term)[i_c][1][i_s][i_r][i_f] != 0. || *(pointer_ME2term)[i_c][2][i_s][i_r][i_f] != 0.){
		  
		  // pdf_factor   []_+(1)   and   delta + []_+(int)
		  
		  sstemp << setw(27) << name_Sherpa_1.str() << " = " << setw(23) << setprecision(15) << temp_Munich_kpc_endpoint + temp_Munich_kpc_delta;
		  
		  sstemp << " = ";
		  if (*(pointer_ME2term)[i_c][1][i_s][i_r][i_f] != 0.){sstemp << name_Sherpa_1.str() << "_([]_+(1))"; }
		  if (*(pointer_ME2term)[i_c][1][i_s][i_r][i_f] != 0. && *(pointer_ME2term)[i_c][2][i_s][i_r][i_f] != 0.){sstemp << " + ";}
		  if (*(pointer_ME2term)[i_c][2][i_s][i_r][i_f] != 0.){sstemp << name_Sherpa_1.str() << "_(delta + []_+(int))";}
		  sstemp << endl;
		  sstemp << endl;
		  
		  //		sstemp << setw(tablevel1) << "" << setw(12) << "pdf_factor" << " = " << setw(23) << setprecision(15) << 
		  sstemp << setw(tablevel0) << "" << setw(27) << name_Sherpa_pdf_1.str() << " = " << setw(23) << setprecision(15) << 
		    *(pointer_pdf_factor)[i_c][i_z][i_s][i_f][0] << " = " << "sum_i[pdf(a)_i * pdf(b)_i]" << endl;
		  for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
		    if (*(pointer_pdf_factor_combination_1)[i_c][i_z][i_s][i_f][0][i_i] != 0. && *(pointer_pdf_factor_combination_2)[i_c][i_z][i_s][i_f][0][i_i] != 0.){
		      //		    sstemp << setw(tablevel2) << "" << "   " << (*CA_collinear)[i_c][0].all_name()[i_i] << endl;
		      sstemp << setw(tablevel1) << "" << "pdf(a)*pdf(b)_{" << setw(12) << (*CA_collinear)[i_c][0].all_name()[i_i] << "} = " << setw(23) << setprecision(15) << 
			*(pointer_pdf_factor_combination_1)[i_c][i_z][i_s][i_f][0][i_i] * *(pointer_pdf_factor_combination_2)[i_c][i_z][i_s][i_f][0][i_i] << endl;
		      sstemp << setw(tablevel2) << "" << "pdf(a)_" << i_i << " = " << setw(23) << setprecision(15) << 
			*(pointer_pdf_factor_combination_1)[i_c][i_z][i_s][i_f][0][i_i] << endl;
		      sstemp << setw(tablevel2) << "" << "pdf(b)_" << i_i << " = " << setw(23) << setprecision(15) << 
			*(pointer_pdf_factor_combination_2)[i_c][i_z][i_s][i_f][0][i_i] << endl;
		    }
		  }
		  double temp_pdf_factor = 0.;
		  for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
		    temp_pdf_factor += *(pointer_pdf_factor_combination_1)[i_c][i_z][i_s][i_f][0][i_i] * *(pointer_pdf_factor_combination_2)[i_c][i_z][i_s][i_f][0][i_i];
		  }
		  temp_Sherpa_KPterms_pdf[i_c] += (temp_Sherpa_KPterms_delta[i_c] + temp_Sherpa_KPterms_endpoint[i_c]) * temp_pdf_factor;

		  sstemp << setw(tablevel0) << "" << setw(27) << "   check" << " = " << setw(23) << setprecision(15) << 
		    temp_pdf_factor << endl;
		  
		  sstemp << endl;
		 
		  sstemp << endl;
		}
	      }
	      
	      
	      else {
		
		// regular + []_+(z)
		
		if (i_z != (*CA_collinear)[i_c][0].no_emitter()){continue;}
		if (*(pointer_ME2term)[i_c][0][i_s][i_r][i_f] != 0){
		  stringstream name_Sherpa_z;
		  name_Sherpa_z << "kpc";
		  if ((*CA_collinear)[i_c][0].no_emitter() == 1){name_Sherpa_z << "a";}
		  if ((*CA_collinear)[i_c][0].no_emitter() == 2){name_Sherpa_z << "b";}
		  if ((*CA_collinear)[i_c][0].type() == 0 || (*CA_collinear)[i_c][0].type() == 3){name_Sherpa_z << "[3]";}
		  if ((*CA_collinear)[i_c][0].type() == 1 || (*CA_collinear)[i_c][0].type() == 2){name_Sherpa_z << "[1]";}
		  
		  stringstream name_Sherpa_pdf_z;
		  name_Sherpa_pdf_z << "f";
		  if ((*CA_collinear)[i_c][0].no_emitter() == 1){name_Sherpa_pdf_z << "a";}
		  if ((*CA_collinear)[i_c][0].no_emitter() == 2){name_Sherpa_pdf_z << "b";}
		  if ((*CA_collinear)[i_c][0].type() == 0 || (*CA_collinear)[i_c][0].type() == 3){name_Sherpa_pdf_z << "gx";}
		  if ((*CA_collinear)[i_c][0].type() == 1 || (*CA_collinear)[i_c][0].type() == 2){name_Sherpa_pdf_z << "qx";}
		  name_Sherpa_pdf_z << " * f";
		  if ((*CA_collinear)[i_c][0].no_emitter() == 1){name_Sherpa_pdf_z << "b";}
		  if ((*CA_collinear)[i_c][0].no_emitter() == 2){name_Sherpa_pdf_z << "a";}
		  
		  //  m_kpca[0]*faqx
		  //  m_kpca[2]*fagx
		  //  m_kpcb[0]*fbqx
		  //  m_kpcb[2]*fbgx
		  
		  double temp_Munich_KP_regular = 0.;	  
		  double temp_Munich_kpc_regular = 0.;
		  
		  temp_Munich_KP_regular = *(pointer_ME2term)[i_c][0][i_s][i_r][i_f]
		    * *(pointer_relative_factor_alpha_S)[0][i_s][i_r]
		    //		  * *(pointer_pdf_factor)[i_c][i_z][i_s][i_f][0] 
		    * temp_symmetry_factor // !!! factor shifted between different contributions !!!
		    ;
		  temp_Munich_kpc_regular = temp_Munich_KP_regular / g_z_coll[(*CA_collinear)[i_c][0].no_emitter()];
		  temp_Sherpa_KPterms_regular[i_c] += temp_Munich_kpc_regular;
		  ///		  temp_Sherpa_KPterms[i_c] += temp_Munich_kpc_regular;
		  
		  //		sstemp << setw(tablevel0) << "" << left << setw(11) << "reg + []_+(z)" << " = " << setw(23) << setprecision(15) << 
		  sstemp << setw(tablevel0) << "" << left << name_Sherpa_z.str() << setw(20) << "_(reg + []_+(z))" << " = " << setw(23) << setprecision(15) << temp_Munich_kpc_regular
		    		  * z_coll[(*CA_collinear)[i_c][0].no_emitter()] // !!! factor shifted here from PDF's !!! compensated by g_z_coll !!! not here to compare with Sherpa !! -> reactivated !!!
			 << " = " << "K+P_factor / g_z_coll * z_coll" << endl; // includes reactived factor !!!
		  //			 << " = " << "K+P_factor / g_z_coll" << endl;
		  //		       << " = " << "K+P_factor * pdf_factor / g_z_coll" << endl;
		  
		  // K+P_factor   regular + []_+(z)
		  
		  string sum_KPterm = "";
		  stringstream split_KPterm;
		  int counter = 0;
		  double check_sum_KPterm = 0.;
		  for (int j_c = 0; j_c < (*CA_collinear)[i_c].size(); j_c++){
		    string name_Kterm;
		    string name_Pterm;
		    if ((*CA_collinear)[i_c][j_c].no_spectator() == 0){name_Kterm = "Kbar"; }
		    else if ((*CA_collinear)[i_c][j_c].no_spectator() < 3){
		      stringstream temp;
		      temp << (*CA_collinear)[i_c][j_c].no_spectator();
		      name_Kterm = "Kt(" + temp.str() + ")";
		      name_Pterm = "P(" + temp.str() + ")";
		    }
		    else if ((*CA_collinear)[i_c][j_c].no_spectator() > 2){
		      stringstream temp;
		      temp << (*CA_collinear)[i_c][j_c].no_spectator();
		      name_Kterm = "gamma(" + temp.str() + ")";
		      name_Pterm = "P(" + temp.str() + ")";
		    }
		    if (data_K[0][i_c][j_c] != 0.){
		      if (counter > 0){sum_KPterm = sum_KPterm + " + ";}
		      counter++;
		      double temp_Sherpa_K_regular = 
			alpha_x_2pi 
			* data_K[0][i_c][j_c] 
			* *(pointer_relative_factor_alpha_S)[0][i_s][i_r] 
			* z_coll[(*CA_collinear)[i_c][0].no_emitter()] // !!! factor shifted here from PDF's !!! compensated by g_z_coll !!! not here to compare with Sherpa !! -> reactivated !!!
			* temp_symmetry_factor // !!! factor shifted between different contributions !!!
			;
		      check_sum_KPterm += temp_Sherpa_K_regular;
		      sum_KPterm = sum_KPterm + name_Kterm;
		      split_KPterm << setw(tablevel2) << "" << left << setw(12) << name_Kterm << " = " << setw(23) << setprecision(15) << temp_Sherpa_K_regular << endl;
		    }
		    if (value_data_P[v_sf][v_xf][0][i_c][j_c] != 0.){
		      if (counter > 0){sum_KPterm = sum_KPterm + " + ";}
		      counter++;
		      double temp_Sherpa_P_regular = 
			alpha_x_2pi 
			* value_data_P[v_sf][v_xf][0][i_c][j_c] 
			* *(pointer_relative_factor_alpha_S)[0][i_s][i_r] 
			* z_coll[(*CA_collinear)[i_c][0].no_emitter()] // !!! factor shifted here from PDF's !!! compensated by g_z_coll !!! not here to compare with Sherpa !! -> reactivated !!!
			* temp_symmetry_factor // !!! factor shifted between different contributions !!!
			;
		      check_sum_KPterm += temp_Sherpa_P_regular;
		      sum_KPterm = sum_KPterm + name_Pterm;
		      split_KPterm << setw(tablevel2) << "" << left << setw(12) << name_Pterm << " = " << setw(23) << setprecision(15) << temp_Sherpa_P_regular << endl;
		    }
		  }
		  sstemp << setw(tablevel1) << "" << setw(12) << "K+P_factor" << " = " << setw(23) << setprecision(15) << temp_Munich_KP_regular
		    * z_coll[(*CA_collinear)[i_c][0].no_emitter()] // !!! factor shifted here from PDF's !!! compensated by g_z_coll !!! not here to compare with Sherpa !! -> reactivated !!!
			 << setw(9) << "" << " = " << sum_KPterm << endl << split_KPterm.str(); 
		  sstemp << setw(tablevel1) << "" << setw(12) << "   check" << " = " << setw(23) << setprecision(15) << check_sum_KPterm << endl;
		  sstemp << endl;
		  
		  // 1. / g_z_coll   regular + []_+(z)
		  
		  sstemp << setw(tablevel1) << "" << setw(12) << "1/g_z_coll" << " = " << setw(23) << setprecision(15) << 
		    1. / g_z_coll[(*CA_collinear)[i_c][0].no_emitter()] << endl;

		  // z_coll   regular + []_+(z) 
		  sstemp << setw(tablevel1) << "" << setw(12) << "z_coll" << " = " << setw(23) << setprecision(15) << 
		    z_coll[(*CA_collinear)[i_c][0].no_emitter()] << endl; // -> reactivated factor !!!
		  
		  // check   regular + []_+(z)
		  
		  sstemp << setw(tablevel0) << "" << left << setw(27) << "   check" << " = " << setw(23) << setprecision(15) << 
		    check_sum_KPterm
		    //		  * temp_pdf_factor
		    / g_z_coll[(*CA_collinear)[i_c][0].no_emitter()]
			 << endl;
		  sstemp << endl;
		  
		  sstemp << setw(27) << name_Sherpa_z.str() << " = " << setw(23) << setprecision(15) << temp_Munich_kpc_regular * z_coll[(*CA_collinear)[i_c][0].no_emitter()]; // !!! re-appearance maybe here ??? not checked !!! -> could be included into some previous definition (re-activation) !!!
		  sstemp << " = ";
		  if (*(pointer_ME2term)[i_c][0][i_s][i_r][i_f] != 0.){sstemp << name_Sherpa_z.str() << "_(reg + []_+(z))"; }
		  sstemp << endl;
		  sstemp << endl;
		  
		  
		  
		  
		  
		  // pdf_factor   regular + []_+(z)
		  
		  sstemp << setw(tablevel0) << "" << setw(27) << name_Sherpa_pdf_z.str() << " = " << setw(23) << setprecision(15) << 
		    //		sstemp << setw(tablevel1) << "" << setw(12) << "pdf_factor" << " = " << setw(23) << setprecision(15) << 
		    *(pointer_pdf_factor)[i_c][i_z][i_s][i_f][0]
		    / z_coll[(*CA_collinear)[i_c][0].no_emitter()] // !!! factor shifted between different contributions !!! but where does it re-appear ???
			 << " = " << "sum_i[pdf(a)_i * pdf(b)_i]" << endl;
		  for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
		    if (*(pointer_pdf_factor_combination_1)[i_c][i_z][i_s][i_f][0][i_i] != 0. && *(pointer_pdf_factor_combination_2)[i_c][i_z][i_s][i_f][0][i_i] != 0.){
		      //		    sstemp << setw(tablevel2) << "" << "   " << (*CA_collinear)[i_c][0].all_name()[i_i] << endl;
		      sstemp << setw(tablevel1) << "" << "pdf(a)*pdf(b)_{" << setw(12) << (*CA_collinear)[i_c][0].all_name()[i_i] << "} = " << setw(23) << setprecision(15) << 
			*(pointer_pdf_factor_combination_1)[i_c][i_z][i_s][i_f][0][i_i] * *(pointer_pdf_factor_combination_2)[i_c][i_z][i_s][i_f][0][i_i] << endl;
		      sstemp << setw(tablevel2) << "" << "pdf(a)_" << i_i << " = " << setw(23) << setprecision(15) << 
			*(pointer_pdf_factor_combination_1)[i_c][i_z][i_s][i_f][0][i_i] << endl;
		      sstemp << setw(tablevel2) << "" << "pdf(b)_" << i_i << " = " << setw(23) << setprecision(15) << 
			*(pointer_pdf_factor_combination_2)[i_c][i_z][i_s][i_f][0][i_i] << endl;
		    }
		    logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << "   i_z = " << i_z << "   i_s = " << i_s << "   i_f = " << i_f << "   i_i = " << i_i << endl;
		  }
		  double temp_pdf_factor = 0.;
		  for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
		    temp_pdf_factor += *(pointer_pdf_factor_combination_1)[i_c][i_z][i_s][i_f][0][i_i] * *(pointer_pdf_factor_combination_2)[i_c][i_z][i_s][i_f][0][i_i];
		  }
		  temp_Sherpa_KPterms_pdf[i_c] += (temp_Sherpa_KPterms_regular[i_c]) * temp_pdf_factor * z_coll[(*CA_collinear)[i_c][0].no_emitter()] // !!! re-appearance maybe here ??? not checked !!! -> could be included into some previous definition (re-activation) !!!
;
		  sstemp << setw(tablevel0) << "" << setw(27) << "   check" << " = " << setw(23) << setprecision(15) << 
		    temp_pdf_factor << endl;
		  sstemp << endl;
		  
		}
	      }
	    }
	    out_comparison << sstemp.str() << endl;
	    
	  }
	  //	ps_integrand_TSV[i_s][i_r][i_f][i_h][0] = psi_ps_factor * rescaling_factor_alpha_e * accumulate(temp_c.begin(), temp_c.end(), 0.) * *(pointer_relative_factor_alpha_S)[0][i_s][i_r];


    
	  out_comparison << left << endl;
	  for (int i_c = 0; i_c < n_pc; i_c++){
	    out_comparison << "temp_Sherpa_KPterms[" << setw(12) << (*CA_collinear)[i_c][0].name() << "] = " << left << setw(23) << setprecision(15) << temp_Sherpa_KPterms_pdf[i_c] << endl;
	    ///      out_comparison << "temp_Sherpa_KPterms[" << setw(12) << (*CA_collinear)[i_c][0].name() << "] = " << left << setw(23) << setprecision(15) << temp_Sherpa_KPterms[i_c] << endl;
	    out_comparison << "   check             " << setw(12) << "" << " = " << setw(23) << setprecision(15) << 
	      temp_sum_c[i_c] 
	      * temp_symmetry_factor // !!! factor shifted between different contributions !!!
			   << endl;
	    out_comparison << endl;
	  }
	  
	  out_comparison << left << endl;
	  for (int i_c = 0; i_c < n_pc; i_c++){
	    out_comparison << "result[" << setw(12) << (*CA_collinear)[i_c][0].name() << "]" << " = " << left << setw(23) << setprecision(15) << psi.hcf / p_parton[0][0].m2() * temp_Sherpa_KPterms_pdf[i_c] << " = " << "temp_Sherpa_KPterms[" << setw(12) << (*CA_collinear)[i_c][0].name() << "] * hcf * flux_factor * symmetry_factor" << endl;
	    ///      out_comparison << "result[" << setw(12) << (*CA_collinear)[i_c][0].name() << "]" << " = " << left << setw(23) << setprecision(15) << psi.hcf / p_parton[0][0].m2() * temp_Sherpa_KPterms[i_c] << " = " << "temp_Sherpa_KPterms[" << setw(12) << (*CA_collinear)[i_c][0].name() << "] * hcf * flux_factor * symmetry_factor" << endl;
	    out_comparison << left << setw(49) << "" << "temp_Sherpa_KPterms[" << setw(12) << (*CA_collinear)[i_c][0].name() << "] = " << left << setw(23) << setprecision(15) << temp_Sherpa_KPterms_pdf[i_c] << endl;
	    ///      out_comparison << left << setw(49) << "" << "temp_Sherpa_KPterms[" << setw(12) << (*CA_collinear)[i_c][0].name() << "] = " << left << setw(23) << setprecision(15) << temp_Sherpa_KPterms[i_c] << endl;
	    out_comparison << setw(49) << "" << setw(33) << "hcf" << " = " << setw(23) << setprecision(15) << psi.hcf << endl;
	    out_comparison << setw(49) << "" << setw(33) << "flux_factor" << " = " << setw(23) << setprecision(15) << 1. / p_parton[0][0].m2() << endl;
	    out_comparison << setw(49) << "" << setw(33) << "symmetry_factor" << " = " << setw(23) << setprecision(15) << 1. / temp_symmetry_factor << endl;
	    out_comparison << endl;
	  }
	  
	  out_comparison << left << endl;
	  out_comparison << setw(15) << "hcf" << " = " << setw(23) << setprecision(15) << psi.hcf << " = " << "(hbar * c [GeV cm]) ^ 2 * fbarn [cm^2] / 2 / (2pi) ^ (3 * " << csi->n_particle << " - 4)" << endl; //x_pdf[0] * 4. * pow(E, 2) <<
	  /*
	    out_comparison << "psi.xbs_all.size() = " << psi.xbs_all.size() << endl;
	    out_comparison << "psi.xbs_all[0].size() = " << psi.xbs_all[0].size() << endl;
	    out_comparison << endl << "flux factor 1 / s^ = 1 / " << psi.xbs_all[0][0] << " = " << 1. / psi.xbs_all[0][0] << endl << endl; //x_pdf[0] * 4. * pow(E, 2) <<
	  */
	  out_comparison << setw(15) << "flux_factor" << " = " << setw(23) << setprecision(15) << 1. / p_parton[0][0].m2() << " = " << "1 / s^" << " = 1 / " << p_parton[0][0].m2() << endl; //x_pdf[0] * 4. * pow(E, 2) <<
	  out_comparison << setw(15) << "symmetry_factor" << " = " << setw(23) << setprecision(15) << 1. / temp_symmetry_factor << " = 1 / " << temp_symmetry_factor << endl; //x_pdf[0] * 4. * pow(E, 2) <<
	  out_comparison << endl;
	  out_comparison << "not included: phase-space weight (apart from factor from linear collinear-emission mapping, where applicable)" << endl;
	  
	  out_comparison << endl;
	  out_comparison << endl;
	  out_comparison << endl;
    	}
      }
    }
  }




  logger << LOG_DEBUG << "after calculate_pdf_LHAPDF_CA_collinear" << endl;

  out_comparison << endl << "Particle momenta in CMS frame: " << endl << endl;
  output_momenta(out_comparison);

  boost = (x_pdf[1] - x_pdf[2]) / (x_pdf[1] + x_pdf[2]);
  out_comparison << "x_pdf[1] = " << x_pdf[1] << endl;
  out_comparison << "x_pdf[2] = " << x_pdf[2] << endl;
  out_comparison << "boost    = " << boost << endl;

  for (int ib = 1; ib < p_parton[0].size(); ib++){
    p_parton[0][ib] = p_parton[0][ib].zboost(-boost);
  }
  out_comparison << endl << "Particle momenta in LAB frame: " << endl << endl;
  output_momenta(out_comparison);
  out_comparison << endl;

  out_comparison << endl << "Momentum fractions (for pdfs): " << endl << endl;
  for (int i_x = 1; i_x < 3; i_x++){
    out_comparison << "x_pdf[" << i_x << "] = " << setw(23) << showpoint << setprecision(15) << x_pdf[i_x] << endl;
  }
  out_comparison << endl;
  for (int i_x = 1; i_x < 3; i_x++){
    out_comparison << "x_pdf[" << i_x << "] / z_coll[" << i_x << "] = " << setw(23) << showpoint << setprecision(15) << x_pdf[i_x] / z_coll[i_x] << endl;
  }

  out_comparison << endl << "Collinear radiation fractions: " << endl << endl;
  for (int i_z = 1; i_z < 3; i_z++){
    out_comparison << "z_coll[" << i_z << "] = " << setw(23) << showpoint << setprecision(15) << z_coll[i_z] << endl;
  }

  //  out_comparison << endl << "flux factor 1 / s^ = 1 / " << p_parton[0][0].m2() << " = " << 1. / p_parton[0][0].m2() << endl << endl; //x_pdf[0] * 4. * pow(E, 2) <<

  out_comparison.close();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::output_testpoint_RA(ofstream & out_comparison){
  Logger logger("observable_set::output_testpoint_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  out_comparison << "Real-correction matrix element and dipole terms:" << endl << endl;
  for (int i_a = 0; i_a < (*RA_dipole).size(); i_a++){
    out_comparison << setw(12) << (*RA_dipole)[i_a].name() << " = " << setw(23) << setprecision(15) << RA_ME2[i_a] << endl;
    if (i_a == 0){out_comparison << endl;}
  }
  out_comparison << endl;
  out_comparison << "RA    = " << accumulate(RA_ME2.begin(), RA_ME2.end(), 0.) << endl;
  out_comparison << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}





void observable_set::output_testpoint_VT_result(ofstream & out_comparison){
  Logger logger("observable_set::output_testpoint_VA_result");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  out_comparison << setw(12) << " ME2_born = " << setprecision(15) << setw(23) << VA_b_ME2 << endl;
  out_comparison << setw(12) << "  ME2_V+X = " << setprecision(15) << setw(23) << (VA_V_ME2 + VA_X_ME2) << endl;
  out_comparison << endl;
  out_comparison << setw(12) << "    ME2_V = " << setprecision(15) << setw(23) << VA_V_ME2 << endl;
  out_comparison << setw(12) << "    ME2_X = " << setprecision(15) << setw(23) << VA_X_ME2 << endl;
  out_comparison << endl;
  out_comparison << setw(12) << " H1_delta = " << setprecision(15) << setw(23) << QT_H1_delta << endl;
  out_comparison << endl;
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}






void observable_set::output_integrand_maximum(phasespace_set & psi){
  Logger logger("observable_set::output_integrand_maximum");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (abs(this_psp_weight) > max_integrand){
    ofstream out_maxevent;
    out_maxevent.open(filename_maxevent.c_str(), ofstream::out | ofstream::app);  
    max_integrand = abs(this_psp_weight);
    out_maxevent << right << setw(10) << psi_i_gen << right << setw(10) << psi_i_acc << "   max_integrand = " << this_psp_weight << "   " << "max_integrand / sigma_normalization = " << max_integrand / sigma_normalization << endl;
    output_integrand_maximum_psp(psi);
    out_maxevent << endl;
    out_maxevent << "integrand = " << setw(25) << setprecision(16) << integrand << endl;
    out_maxevent << "ME2       = " << setw(25) << setprecision(16) << ME2 << endl;

    out_maxevent << endl;
    out_maxevent.close();
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_integrand_maximum_VA(phasespace_set & psi){
  Logger logger("observable_set::output_integrand_maximum_VA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (abs(this_psp_weight) > max_integrand){
    ofstream out_maxevent;
    out_maxevent.open(filename_maxevent.c_str(), ofstream::out | ofstream::app);

    max_integrand = abs(this_psp_weight);
    out_maxevent << right << setw(10) << psi_i_gen << right << setw(10) << psi_i_acc << "   max_integrand = " << this_psp_weight << "   " << "max_integrand / sigma_normalization = " << max_integrand / sigma_normalization << endl;

    output_integrand_maximum_psp(psi);

    out_maxevent << "integrand = " << setw(25) << setprecision(16) << integrand << endl;
    out_maxevent << "ME2       = " << setw(25) << setprecision(16) << ME2 << endl;
    out_maxevent << "V_ME2     = " << setw(25) << setprecision(16) << VA_V_ME2 << endl;
    out_maxevent << "X_ME2     = " << setw(25) << setprecision(16) << VA_X_ME2 << endl;
    out_maxevent << "I_ME2     = " << setw(25) << setprecision(16) << VA_I_ME2 << endl;
    out_maxevent << "b_ME2     = " << setw(25) << setprecision(16) << VA_b_ME2 << endl;
    out_maxevent << "(V+X+I)_ME2 / b_ME2     = " << setw(25) << setprecision(16) << (VA_V_ME2 + VA_I_ME2 + VA_X_ME2) / VA_b_ME2 << endl;
    out_maxevent << endl;    

    out_maxevent.close();
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::output_integrand_maximum_CA(phasespace_set & psi){//, vector<vector<collinear_set> > & _collinear
  Logger logger("observable_set::output_integrand_maximum_CA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  
  if (fabs(integrand) > fabs(max_integrand)){
    ofstream out_maxevent;
    max_integrand = fabs(integrand);
    out_maxevent.open(filename_maxevent.c_str(), ofstream::out | ofstream::app);  
    
    out_maxevent << right << setw(10) << psi_i_gen << right << setw(10) << psi_i_acc << "   max_integrand = " << this_psp_weight << "   " << "max_integrand / sigma_normalization = " << max_integrand / sigma_normalization << endl;

    output_integrand_maximum_psp(psi);

    for (int ic = 1; ic < 3; ic++){
      out_maxevent << "z_coll[" << ic << "] = " << right << setw(25) << setprecision(16) << showpoint << psi_z_coll[ic] << endl;
    }

    for (int sd = 0; sd < value_mu_fact.size(); sd++){
      for (int ss = 0; ss < value_mu_fact[sd].size(); ss++){
	for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
	  out_maxevent << "CA_value_pdf_factor[" << sd << "][" << ss << "][" << i_a << "][" << (*CA_collinear)[i_a][0].no_emitter() << "][0] = " << right << setw(25) << setprecision(16) << showpoint << CA_value_pdf_factor[sd][ss][i_a][(*CA_collinear)[i_a][0].no_emitter()][0] << endl;
	  out_maxevent << "CA_value_pdf_factor[" << sd << "][" << ss << "][" << i_a << "][0][0] = " << right << setw(25) << setprecision(16) << showpoint << CA_value_pdf_factor[sd][ss][i_a][0][0] << endl;
	}
      }
    }
    out_maxevent << endl;
    out_maxevent.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::output_integrand_maximum_RA(phasespace_set & psi){//, vector<dipole_set> & _dipole
  Logger logger("observable_set::output_integrand_maximum_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (integrand > max_integrand){
    ofstream out_maxevent;
    out_maxevent.open(filename_maxevent.c_str(), ofstream::out | ofstream::app);  

    max_integrand = abs(-accumulate(var_RA_ME2.begin() + 1, var_RA_ME2.end(), -var_RA_ME2[0]) * psi_ps_factor * rescaling_factor_alpha_e);

    out_maxevent << right << setw(10) << psi_i_gen << right << setw(10) << psi_i_acc << "   max_integrand = " << -accumulate(var_RA_ME2.begin() + 1, var_RA_ME2.end(), -var_RA_ME2[0]) * psi_ps_factor * rescaling_factor_alpha_e << "   " << "max_integrand / sigma_normalization = " << max_integrand / sigma_normalization << endl;
    output_integrand_maximum_psp(psi);

    /*
    for (int ib = 0; ib < p_parton[0].size(); ib++){out_maxevent << "p_parton[0][" << ib << "] = " << p_parton[0][ib] << endl;}

  // write out all relevant pairings to spot badly mapped resonances
  out_maxevent << endl;
  for (int ib1 = 3; ib1 < p_parton[0].size(); ib1++){
    for (int ib2 = ib1+1; ib2 < p_parton[0].size(); ib2++){
      out_maxevent << "p_parton[0][" << ib1 << "]+p_parton[0][" << ib2 << "] = " << p_parton[0][ib1]+p_parton[0][ib2] << ", sqrt(s)= "  << (p_parton[0][ib1]+p_parton[0][ib2]).m() << endl;
    }
  }
  for (int ib1 = 3; ib1 < p_parton[0].size(); ib1++){
    for (int ib2 = ib1+1; ib2 < p_parton[0].size(); ib2++){
      for (int ib3 = ib2+1; ib3 < p_parton[0].size(); ib3++){
      out_maxevent << "p_parton[0][" << ib1 << "]+p_parton[0][" << ib2 << "]+p_parton[0][" << ib3 << "] = " << p_parton[0][ib1]+p_parton[0][ib2]+p_parton[0][ib3] << ", sqrt(s)= "  << (p_parton[0][ib1]+p_parton[0][ib2]+p_parton[0][ib3]).m() << endl;
      }
    }
  }
  out_maxevent << endl;
    */
    int ccount = 0;
    for (int i_a = 1; i_a < cut_ps.size(); i_a++){if (cut_ps[i_a] == -1){ccount++;}}
    out_maxevent << "cut dipoles: " << setw(2) << ccount << endl;
    out_maxevent << left << setw(8) << "" << "   " << right << setw(25) << setprecision(16) << showpoint << "" << "   " << right << setw(25) << setprecision(16) << showpoint << "A_ME2" << "   " << right << setw(25) << setprecision(16) << showpoint << "A_ME2 / R_ME2" << endl; // "   " << right << setw(25) << setprecision(16) << showpoint << "m_top" << "   " << right << setw(25) << setprecision(16) << showpoint << "m_atop" << "   " << endl;
  logger << LOG_DEBUG_VERBOSE << "cut_ps.size() = " << cut_ps.size() << endl;

  for (int i_a = 0; i_a < cut_ps.size(); i_a++){
    out_maxevent << left << setw(8) << (*RA_dipole)[i_a].name() << "   " << right << setw(25) << setprecision(16) << showpoint << (*RA_dipole)[i_a].xy << "   " << right << setw(25) << setprecision(16) << showpoint << RA_ME2[i_a] << "   " << right << setw(25) << setprecision(16) << showpoint << RA_ME2[i_a] / RA_ME2[0] << "   " << endl; //right << setw(25) << setprecision(16) << showpoint << (p_parton[i_a][3] + p_parton[i_a][5]).m() << "   " << right << setw(25) << setprecision(16) << showpoint << (p_parton[i_a][4] + p_parton[i_a][6]).m() << endl;
  }
  out_maxevent << endl;
  out_maxevent << left << setw(8) << "Sum(dip) = " << right << setw(25) << setprecision(16) << showpoint << 0 << "   " << right << setw(25) << setprecision(16) << showpoint << accumulate(RA_ME2.begin() + 1, RA_ME2.end(), 0.) << "   " << right << setw(25) << setprecision(16) << showpoint << accumulate(RA_ME2.begin() + 1, RA_ME2.end(), 0.) / RA_ME2[0] << "   " << endl; //right << setw(25) << setprecision(16) << showpoint << (d  
  //  if (ccount != 0){
    out_maxevent << "Number of cut dipoles: " << ccount << endl;
    //    out_maxevent << left << setw(8) << "" << "   " << right << setw(25) << setprecision(16) << showpoint << "m_bb" << right << setw(25) << setprecision(16) << showpoint << "p_T(b)" << right << setw(25) << setprecision(16) << showpoint << "p_T(b~)" << endl;	  
    for (int i_a = 0; i_a < cut_ps.size(); i_a++){
      for (int ib = 0; ib < p_parton[i_a].size(); ib++){out_maxevent << "p_parton[" << i_a << "][" << ib << "] = " << p_parton[i_a][ib] << endl;}
      //      out_maxevent << left << setw(8) << dipole_name[i_a] << "   " << right << setw(25) << setprecision(16) << showpoint << (p_parton[i_a][7] + p_parton[i_a][8]).m() << right << setw(25) << setprecision(16) << showpoint << p_parton[i_a][7].pT() << right << setw(25) << setprecision(16) << showpoint << (p_parton[i_a][8]).pT() << endl;
    }
    //  }
  //  else {
    for (int i_a = 1; i_a < cut_ps.size(); i_a++){
      if      ((*RA_dipole)[i_a].type_dipole() == 1){
	double y_ij_k = (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()]) / (p_parton[0][(*RA_dipole)[i_a].no_R_spectator()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] + p_parton[0][(*RA_dipole)[i_a].no_R_spectator()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()] + p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()]);
	double z_i = (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]) / (p_parton[0][(*RA_dipole)[i_a].no_R_spectator()] * (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] + p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()]));
	double z_j = (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()] * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]) / (p_parton[0][(*RA_dipole)[i_a].no_R_spectator()] * (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] + p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()]));
	out_maxevent << setw(2) << i_a << "   y_" << (*RA_dipole)[i_a].no_R_emitter_1() << (*RA_dipole)[i_a].no_R_emitter_2() << "," << (*RA_dipole)[i_a].no_R_spectator() << " = " << setw(25) << setprecision(16) << y_ij_k << "   z_" << (*RA_dipole)[i_a].no_R_emitter_1() << " = " << setw(25) << setprecision(16) << z_i << "      z_" << (*RA_dipole)[i_a].no_R_emitter_2() << " = " << setw(25) << setprecision(16) << z_j << endl;
	out_maxevent << setw(2) << i_a << "   y_" << (*RA_dipole)[i_a].no_R_emitter_1() << (*RA_dipole)[i_a].no_R_emitter_2() << "," << (*RA_dipole)[i_a].no_R_spectator() << " = " << setw(25) << setprecision(16) << (*RA_dipole)[i_a].xy << endl;//"   z_" << (*RA_dipole)[i_a].no_R_emitter_1() << " = " << setw(25) << setprecision(16) << z_i << "      z_" << (*RA_dipole)[i_a].no_R_emitter_2() << " = " << setw(25) << setprecision(16) << z_j << endl;
      }
      else if ((*RA_dipole)[i_a].type_dipole() == 2){
	double x_ij_a = (p_parton[0][(*RA_dipole)[i_a].no_R_spectator()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] + p_parton[0][(*RA_dipole)[i_a].no_R_spectator()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()] - p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()]) / ((p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] + p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()]) * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]);
	double z_i = (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]) / (p_parton[0][(*RA_dipole)[i_a].no_R_spectator()] * (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] + p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()]));
	double z_j = (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()] * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]) / (p_parton[0][(*RA_dipole)[i_a].no_R_spectator()] * (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] + p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()]));
	out_maxevent << setw(2) << i_a << "   x_" << (*RA_dipole)[i_a].no_R_emitter_1() << (*RA_dipole)[i_a].no_R_emitter_2() << "," << (*RA_dipole)[i_a].no_R_spectator() << " = " << setw(25) << setprecision(16) << x_ij_a << "   z_" << (*RA_dipole)[i_a].no_R_emitter_1() << " = " << setw(25) << setprecision(16) << z_i << "      z_" << (*RA_dipole)[i_a].no_R_emitter_2() << " = " << setw(25) << setprecision(16) << z_j << endl;
	out_maxevent << setw(2) << i_a << "   x_" << (*RA_dipole)[i_a].no_R_emitter_1() << (*RA_dipole)[i_a].no_R_emitter_2() << "," << (*RA_dipole)[i_a].no_R_spectator() << " = " << setw(25) << setprecision(16) << (*RA_dipole)[i_a].xy << endl;//"   z_" << (*RA_dipole)[i_a].no_R_emitter_1() << " = " << setw(25) << setprecision(16) << z_i << "      z_" << (*RA_dipole)[i_a].no_R_emitter_2() << " = " << setw(25) << setprecision(16) << z_j << endl;
      }
      else if ((*RA_dipole)[i_a].type_dipole() == 3){
	double x_ik_a = (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()] + p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()] - p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()] * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]) / ((p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()] + p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]) * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()]);
	double u_i = (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()]) / (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()] + p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]));
	out_maxevent << setw(2) << i_a << "   x_" << (*RA_dipole)[i_a].no_R_emitter_2() << (*RA_dipole)[i_a].no_R_spectator() << "," << (*RA_dipole)[i_a].no_R_emitter_1() << " = " << setw(25) << setprecision(16) << x_ik_a << "   u_" << (*RA_dipole)[i_a].no_R_spectator() << " = " << setw(25) << setprecision(16) << u_i << endl;
	out_maxevent << setw(2) << i_a << "   x_" << (*RA_dipole)[i_a].no_R_emitter_2() << (*RA_dipole)[i_a].no_R_spectator() << "," << (*RA_dipole)[i_a].no_R_emitter_1() << " = " << setw(25) << setprecision(16) << (*RA_dipole)[i_a].xy << endl;//"   u_" << (*RA_dipole)[i_a].no_R_spectator() << " = " << setw(25) << setprecision(16) << u_i << endl;
      }
      else if ((*RA_dipole)[i_a].type_dipole() == 5){
	double x_i_ab = (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()] - p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] - p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()] * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]) / (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]);
	out_maxevent << setw(2) << i_a << "   x_" << (*RA_dipole)[i_a].no_R_spectator() << "," << (*RA_dipole)[i_a].no_R_emitter_1() << (*RA_dipole)[i_a].no_R_emitter_2() << " = " << setw(25) << setprecision(16) << x_i_ab << "   v_" << (*RA_dipole)[i_a].no_R_emitter_2() << " = " << setw(25) << setprecision(16) << (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()]) / (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]) << endl;
	out_maxevent << setw(2) << i_a << "   x_" << (*RA_dipole)[i_a].no_R_spectator() << "," << (*RA_dipole)[i_a].no_R_emitter_1() << (*RA_dipole)[i_a].no_R_emitter_2() << " = " << setw(25) << setprecision(16) << (*RA_dipole)[i_a].xy << endl;//"   v_" << (*RA_dipole)[i_a].no_R_emitter_2() << " = " << setw(25) << setprecision(16) << (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_emitter_2()]) / (p_parton[0][(*RA_dipole)[i_a].no_R_emitter_1()] * p_parton[0][(*RA_dipole)[i_a].no_R_spectator()]) << endl;
      }
    }
    //  }
  for (int sr = 0; sr < psi_RA_singular_region_list.size(); sr++){
    int x1 = psi_RA_singular_region_list[sr][0];
    int x2 = psi_RA_singular_region_list[sr][1];
    psi_RA_singular_region[x1][x2] = p_parton[0][x1] * p_parton[0][x2] / psi_xbs_all[0][0];
    out_maxevent << psi_RA_singular_region_name[x1][x2] << "/^s = "<< right << setw(25) << setprecision(16) << showpoint << psi_RA_singular_region[x1][x2] << endl; //"   A/R = " << right << setw(25) << setprecision(16) << showpoint << accumulate(RA_ME2.begin() + 1, RA_ME2.end(), 0.) / RA_ME2[0] << endl;
  }
  out_maxevent << right << setw(14) << psi_i_gen << right << setw(12) << psi_i_acc << "   maxI = " << setw(15) << setprecision(8) << -accumulate(var_RA_ME2.begin() + 1, var_RA_ME2.end(), -var_RA_ME2[0]) * psi_ps_factor * rescaling_factor_alpha_e << "   " << "maxI/LO = " << setw(15) << setprecision(8) << max_integrand / sigma_normalization << "   A/R = " << setw(25) << setprecision(16) << accumulate(RA_ME2.begin() + 1, RA_ME2.end(), 0.) / RA_ME2[0] << "   n_cut-dipoles: " << setw(2) << ccount << endl;
  out_maxevent << endl << endl << endl;
  //  cout << right << setw(14) << psi_i_gen << right << setw(12) << psi_i_acc << "   maxI = " << setw(15) << setprecision(8) << -accumulate(var_RA_ME2.begin() + 1, var_RA_ME2.end(), -var_RA_ME2[0]) * psi_ps_factor * rescaling_factor_alpha_e << "   " << "maxI/LO = " << setw(15) << setprecision(8) << max_integrand / sigma_normalization << "   A/R = " << setw(25) << setprecision(16) << accumulate(RA_ME2.begin() + 1, RA_ME2.end(), 0.) / RA_ME2[0] << "   n_cut-dipoles: " << setw(2) << ccount << endl;
  out_maxevent.close();
 }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::output_integrand_cancellation_check_RA(phasespace_set & psi){//, vector<dipole_set> & _dipole
  Logger logger("observable_set::output_integrand_cancellation_check_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int flag = 0;
  for (int sr = 0; sr < psi_RA_singular_region_list.size(); sr++){
    int x1 = psi_RA_singular_region_list[sr][0];
    int x2 = psi_RA_singular_region_list[sr][1];
    psi_RA_singular_region[x1][x2] = p_parton[0][x1] * p_parton[0][x2] / psi_xbs_all[0][0];
    if (psi_RA_singular_region[x1][x2] < psi_cut_technical * 1.e1){
      //      if (x1 == 2 || x2 == 2){
      flag++;
      //      }
    }
  }
  if (flag){
    int ccount = 0;
    for (int i_a = 1; i_a < cut_ps.size(); i_a++){if (cut_ps[i_a] == -1){ccount++;}}
    //    cout << "cut dipoles: " << setw(2) << ccount << endl;
    //    cout << left << setw(8) << "" << "   " << right << setw(25) << setprecision(16) << showpoint << "" << "   " << right << setw(25) << setprecision(16) << showpoint << "A_ME2" << "   " << right << setw(25) << setprecision(16) << showpoint << "A_ME2 / R_ME2" << endl;

    cout << right << setw(14) << psi_i_gen << right << setw(12) << psi_i_acc << "   maxI = " << setw(15) << setprecision(8) << -accumulate(var_RA_ME2.begin() + 1, var_RA_ME2.end(), -var_RA_ME2[0]) * psi_ps_factor * rescaling_factor_alpha_e << "   " << "maxI/LO = " << setw(15) << setprecision(8) << max_integrand / sigma_normalization << "   A/R = " << setw(25) << setprecision(16) << accumulate(RA_ME2.begin() + 1, RA_ME2.end(), 0.) / RA_ME2[0] << "   n_cut-dipoles: " << setw(2) << ccount << endl;
    for (int ib = 0; ib < p_parton[0].size(); ib++){cout << "p_parton[" << 0 << "][" << ib << "] = " << p_parton[0][ib] << endl;}
    for (int i_a = 0; i_a < cut_ps.size(); i_a++){
      cout << left << setw(8) << (*RA_dipole)[i_a].name() << "   " << right << setw(25) << setprecision(16) << showpoint << (*RA_dipole)[i_a].xy << "   " << right << setw(25) << setprecision(16) << showpoint << RA_ME2[i_a] << "   " << right << setw(25) << setprecision(16) << showpoint << RA_ME2[i_a] / RA_ME2[0] << "   " << endl;
    }
    for (int sr = 0; sr < psi_RA_singular_region_list.size(); sr++){
      int x1 = psi_RA_singular_region_list[sr][0];
      int x2 = psi_RA_singular_region_list[sr][1];
      cout << psi_RA_singular_region_name[x1][x2] << "/^s = "<< right << setw(25) << setprecision(16) << showpoint << psi_RA_singular_region[x1][x2] << endl;
    }
    
    cout << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::output_integrand_maximum_VT2(phasespace_set & psi){
  Logger logger("observable_set::output_integrand_maximum_VT2");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // TODO: put into perform.maxintegrand.checks
  if (fabs(this_psp_weight) > max_integrand) {
    max_integrand = fabs(this_psp_weight);
    cout << psi_i_acc << ", maxint: " << max_integrand << endl;
    cout << psi_zz_pdf[1] << ", " << psi_zz_pdf[2] << endl;
  }
  //  cout << QT_H1_delta << ", " << QT_H2_delta << endl;
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::output_testpoint_born(ofstream & out_comparison){
  Logger logger("observable_set::output_testpoint_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  stringstream scouplings_ss;
  scouplings_ss << "(alpha_s^" << contribution_order_alpha_s << " alpha_e^" << contribution_order_alpha_e << ")";
  string scouplings = scouplings_ss.str();
  double couplings = pow(alpha_S, contribution_order_alpha_s) * pow(msi.alpha_e, contribution_order_alpha_e);

  out_comparison << "ME2(production): " << endl;
  out_comparison << "ME2_born = " << setprecision(15) << showpoint << setw(25) << value_ME2term[0] << endl << endl;
  out_comparison << "ME2(production), devided by coupling constants: " << endl;
  out_comparison << "ME2_born / " << scouplings << " = " << setprecision(15) << showpoint << setw(25) << value_ME2term[0] / couplings << endl;


  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::output_testpoint_collinear(phasespace_set & psi){
  Logger logger("observable_set::output_testpoint_collinear");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;



  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::output_testpoint(phasespace_set & psi){
  Logger logger("observable_set::output_testpoint");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<int> list_external_mass;
  for (int i_p = 1; i_p < csi->type_parton[0].size(); i_p++){
    if (mass_parton[0][i_p] != 0.){
      int flag = -1;
      for (int i_l = 0; i_l < list_external_mass.size(); i_l++){
	if (abs(csi->type_parton[0][i_p]) == list_external_mass[i_l]){flag = i_l; break;}
      }
      if (flag == -1){list_external_mass.push_back(abs(csi->type_parton[0][i_p]));}
    }
  }


  int nd = 0;
  for (int i_a = 0; i_a < cut_ps.size(); i_a++){if (cut_ps[i_a] > -1){nd++;}}
    nd = cut_ps.size(); // would print every phase-space point
  if (nd == cut_ps.size()){
    //    int tab_txt = 1;
    if (list_external_mass.size() > 0){
      cout << "if (";
      for (int i_l = 0; i_l < list_external_mass.size(); i_l++){
	cout << "osi_M[" << list_external_mass[i_l] << "] == " << setprecision(6) << psi_M[list_external_mass[i_l]];
	if (i_l < list_external_mass.size() - 1){cout << " && ";}
      }
      cout << "){" << endl;
    }
    else {cout << "{" << endl;}

    if (type_contribution == "CA" || 
	type_contribution == "RCA"){
      for (int j = 0; j < 3; j++){
	cout << "  osi_x_pdf[" << j << "] = " << right << setw(25) << setprecision(16) << psi_x_pdf[j] << ";" << endl;
      }
      cout << endl;
      for (int j = 1; j < 3; j++){
	cout << "  osi_z_coll[" << j << "] = " << right << setw(25) << setprecision(16) << psi_z_coll[j] << ";" << endl;
      }
      cout << endl;
    }

    for (int i_p = 0; i_p < p_parton[0].size(); i_p++){
      cout << "  osi_p_parton[0][" << i_p << "] = fourvector(" << right << setw(25) << setprecision(16) << p_parton[0][i_p].x0() << "," << right << setw(25) << setprecision(16) << p_parton[0][i_p].x1() << "," << right << setw(25) << setprecision(16) << p_parton[0][i_p].x2() << "," << right << setw(25) << setprecision(16) << p_parton[0][i_p].x3() << ");" << endl;
    }
    cout << "}" << endl;
    //    if (list_external_mass.size() > 0){cout << "}" << endl;}
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}





void observable_set::perform_integration_step_complete(phasespace_set & psi){
  Logger logger("observable_set::perform_integration_step_complete");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (temp_n_step == 0){temp_n_step = psi_n_step;}

  //  calculate_intermediate_result(psi);

  if (switch_output_integration){
    output_step_integration(psi);
    output_step_integration_TSV(psi);
  }
  if (switch_output_result){
    output_step_result(psi);
    output_step_result_TSV(psi);
  }

  if (switch_output_gnuplot){output_step_errorplot(psi);}
  if (switch_output_execution){output_step_execution(psi);}

  //  cout << "psi_i_tec = " << psi_i_tec << endl;
  
  if (switch_output_moment){
    for (int nm = 0; nm < n_moments; nm++){
      //      new_output_integration(psi_i_gen, psi_i_rej, psi_i_acc, psi_i_nan, int_end, full_sum_moment[nm], full_sum_moment2[nm], step_sum_moment[nm], step_sum_moment2[nm], full_sum_moment_CV[nm], full_sum_moment2_CV[nm], step_sum_moment_CV[nm], step_sum_moment2_CV[nm], filename_moment[nm], 0, temp_n_step, psi, oset);
    }
  }
  
  if (psi_i_acc == 10 * temp_n_step){
    temp_n_step = 10 * temp_n_step;
    cout << "temp_n_step = " << temp_n_step << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_zero_contribution_complete(phasespace_set & psi){
  Logger logger("observable_set::output_zero_contribution_complete");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  /*
  if (psi_i_gen==1) {
    // avoid runs with only one event, to avoid division by zero when computing the sample variance
    psi_i_gen++;
  }
  */

  // needs to be updated !!!
  /*
  psi.weight_output_MCweight_optimization();
  psi_random_manager.writeout_weights();
  output_weight_vegas(psi_filename_tauweight, psi_tau_alpha);
  output_weight_vegas(psi_filename_x1x2weight, psi_x1x2_alpha);
  */
  logger << LOG_INFO << "before psi.output_optimization_complete" << endl; 
  logger << LOG_INFO << "psi.MC_opt_end = " << psi.MC_opt_end << endl;
  logger << LOG_INFO << "psi.active_optimization = " << psi.active_optimization << endl;
  logger << LOG_INFO << "psi.end_optimization    = " << psi.end_optimization << endl;

  logger << LOG_INFO << "psi.MC_phasespace.active_optimization = " << psi.MC_phasespace.active_optimization << endl;
  logger << LOG_INFO << "psi.MC_phasespace.end_optimization    = " << psi.MC_phasespace.end_optimization << endl;
  logger << LOG_INFO << "psi.MC_tau.active_optimization = " << psi.MC_tau.active_optimization << endl;
  logger << LOG_INFO << "psi.MC_tau.end_optimization    = " << psi.MC_tau.end_optimization << endl;
  if (csi->class_contribution_CS_real){
    for (int i_a = 1; i_a < psi.MC_x_dipole.size(); i_a++){
      logger << LOG_INFO << "psi.MC_x_dipole[" << i_a << "].active_optimization = " << psi.MC_x_dipole[i_a].active_optimization << endl;
      logger << LOG_INFO << "psi.MC_x_dipole[" << i_a << "].end_optimization    = " << psi.MC_x_dipole[i_a].end_optimization << endl;
    }
  }
  logger << LOG_INFO << "psi.IS_tau.active_optimization = " << psi.IS_tau.active_optimization << endl;
  logger << LOG_INFO << "psi.IS_tau.end_optimization    = " << psi.IS_tau.end_optimization << endl;
  logger << LOG_INFO << "psi.IS_x1x2.active_optimization = " << psi.IS_x1x2.active_optimization << endl;
  logger << LOG_INFO << "psi.IS_x1x2.end_optimization    = " << psi.IS_x1x2.end_optimization << endl;
  if (csi->class_contribution_CS_collinear){
    for (int i_z = 1; i_z < 3; i_z++){
       logger << LOG_INFO << "psi.IS_z1z2[" << i_z << "].active_optimization = " << psi.IS_z1z2[i_z].active_optimization << endl;
      logger << LOG_INFO << "psi.IS_z1z2[" << i_z << "].end_optimization    = " << psi.IS_z1z2[i_z].end_optimization << endl;
    }
  }

  
  /*
  logger << LOG_INFO << " = " <<  << endl;
  logger << LOG_INFO << " = " <<  << endl;
  */
  psi.MC_opt_end = 1;
  psi.end_optimization = 1;
  if (psi.MC_phasespace.active_optimization){psi.MC_phasespace.end_optimization = 1;}
  if (psi.MC_tau.active_optimization){psi.MC_tau.end_optimization = 1;}
  if (csi->class_contribution_CS_real){
    for (int i_a = 1; i_a < psi.MC_x_dipole.size(); i_a++){
      if (psi.MC_x_dipole[i_a].active_optimization){psi.MC_x_dipole[i_a].end_optimization = 1;}
    }
  }
  if (psi.IS_tau.active_optimization){psi.IS_tau.end_optimization = 1;}
  if (psi.IS_x1x2.active_optimization){psi.IS_x1x2.end_optimization = 1;}
  if (csi->class_contribution_CS_collinear){
    for (int i_z = 1; i_z < 3; i_z++){
      if (psi.IS_z1z2[i_z].active_optimization){psi.IS_z1z2[i_z].end_optimization = 1;}
    }
  }
  
  psi.random_manager.end_optimization = 1;
  
  logger << LOG_INFO << "before psi.output_optimization_complete" << endl;
  logger << LOG_INFO << "psi.MC_opt_end = " << psi.MC_opt_end << endl;
  logger << LOG_INFO << "psi.active_optimization = " << psi.active_optimization << endl;
  logger << LOG_INFO << "psi.end_optimization    = " << psi.end_optimization << endl;

  logger << LOG_INFO << "psi.MC_phasespace.active_optimization = " << psi.MC_phasespace.active_optimization << endl;
  logger << LOG_INFO << "psi.MC_phasespace.end_optimization    = " << psi.MC_phasespace.end_optimization << endl;
  logger << LOG_INFO << "psi.MC_tau.active_optimization = " << psi.MC_tau.active_optimization << endl;
  logger << LOG_INFO << "psi.MC_tau.end_optimization    = " << psi.MC_tau.end_optimization << endl;
  if (csi->class_contribution_CS_real){
    for (int i_a = 1; i_a < psi.MC_x_dipole.size(); i_a++){
      logger << LOG_INFO << "psi.MC_x_dipole[" << i_a << "].active_optimization = " << psi.MC_x_dipole[i_a].active_optimization << endl;
      logger << LOG_INFO << "psi.MC_x_dipole[" << i_a << "].end_optimization    = " << psi.MC_x_dipole[i_a].end_optimization << endl;
    }
  }
  logger << LOG_INFO << "psi.IS_tau.active_optimization = " << psi.IS_tau.active_optimization << endl;
  logger << LOG_INFO << "psi.IS_tau.end_optimization    = " << psi.IS_tau.end_optimization << endl;
  logger << LOG_INFO << "psi.IS_x1x2.active_optimization = " << psi.IS_x1x2.active_optimization << endl;
  logger << LOG_INFO << "psi.IS_x1x2.end_optimization    = " << psi.IS_x1x2.end_optimization << endl;
  if (csi->class_contribution_CS_collinear){
    for (int i_z = 1; i_z < 3; i_z++){
       logger << LOG_INFO << "psi.IS_z1z2[" << i_z << "].active_optimization = " << psi.IS_z1z2[i_z].active_optimization << endl;
      logger << LOG_INFO << "psi.IS_z1z2[" << i_z << "].end_optimization    = " << psi.IS_z1z2[i_z].end_optimization << endl;
    }
  }



  
  //  psi.result_optimization_complete();
  psi.output_optimization_complete();

  logger << LOG_INFO << "after psi.output_optimization_complete" << endl; 
  
  //      output_step_integration(psi);
  if (temp_n_step == 0){temp_n_step = psi_n_step;}
  if (switch_output_integration){
    output_step_integration(psi);
    output_step_integration_TSV(psi);
  }

  if (switch_output_result){
    output_step_result(psi);
    output_step_result_TSV(psi);
  }
  if (switch_output_execution){output_step_execution(psi);}
  
  if (switch_output_distribution){output_distribution_complete(psi);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_step_integration_TSV(phasespace_set & psi){
  Logger logger("observable_set::output_step_integration_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
    ofstream out_integration;
    out_integration.open((filename_integration_TSV[i_s]).c_str(), ofstream::out | ofstream::app);

    int next_no_qTcut = 0;
    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
      if (i_q != no_qTcut_integration[next_no_qTcut]){continue;}

      logger << LOG_DEBUG << "i_q = " << i_q << "   no_qTcut_integration[next_no_qTcut = " << next_no_qTcut << "] = " << no_qTcut_integration[next_no_qTcut] << endl;

      out_integration << endl;
      if (active_qTcut == 0){
	out_integration << left << setw(21) << "qTcut - independent";
      }
      else {
	out_integration << left << setw(8) << "qTcut = " << setprecision(3) << setw(6) << showpoint << value_qTcut[i_q] << setw(7) << " GeV   ";
      }
      if (i_q == 0){
	out_integration << "n(gen) = " << setw(12) << psi_i_gen << "   n(acc) = " << setw(12) << psi_i_acc << "   n(rej) = " << setw(12) << psi_i_rej << "   n(nan) = " << setw(12) << psi_i_nan << "   n(tec) = " << setw(12) << psi_i_tec;
      }
      out_integration << endl;
      out_integration << endl;

      double abs_central_Xsection = abs(Xsection_TSV[i_s][i_q][no_central_scale_ren_TSV[i_s]][no_central_scale_fact_TSV[i_s]]);
      
      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	  double delta_NLO_LO = abs(Xsection_delta_TSV[i_s][i_q][i_r][i_f] / sigma_normalization);
	  
	  int shift;
	  string unit;
	  double output_Xsection_delta = 0.;
	  double output_Xsection = 0.;
	  if (abs_central_Xsection >= 1.e3){output_Xsection = Xsection_TSV[i_s][i_q][i_r][i_f] / 1.e3; output_Xsection_delta = Xsection_delta_TSV[i_s][i_q][i_r][i_f] / 1.e3; unit = "pb";}// shift = -3;}
	  else if (abs_central_Xsection >= 1.e0){output_Xsection = Xsection_TSV[i_s][i_q][i_r][i_f]; output_Xsection_delta = Xsection_delta_TSV[i_s][i_q][i_r][i_f]; unit = "fb"; shift = 0;}
	  else if (abs_central_Xsection >= 1.e-3){output_Xsection = Xsection_TSV[i_s][i_q][i_r][i_f] * 1.e3; output_Xsection_delta = Xsection_delta_TSV[i_s][i_q][i_r][i_f] * 1.e3; unit = "ab";}// shift = 3;}
	  else if (abs_central_Xsection >= 1.e-6){output_Xsection = Xsection_TSV[i_s][i_q][i_r][i_f] * 1.e6; output_Xsection_delta = Xsection_delta_TSV[i_s][i_q][i_r][i_f] * 1.e6; unit = "zb";}// shift = 6;}
	  else if (abs_central_Xsection >= 1.e-9){output_Xsection = Xsection_TSV[i_s][i_q][i_r][i_f] * 1.e9; output_Xsection_delta = Xsection_delta_TSV[i_s][i_q][i_r][i_f] * 1.e9; unit = "yb";}// shift = 9;}
	  /*
	  double output_Xsection = 0.;
	  if (abs(Xsection) >= 1.e3){output_Xsection = Xsection_TSV[i_s][i_q][i_r][i_f] / 1.e3; output_Xsection_delta = Xsection_delta_TSV[i_s][i_q][i_r][i_f] / 1.e3; unit = "pb";}// shift = -3;}
	  else if (abs(Xsection) >= 1.e0){output_Xsection = Xsection_TSV[i_s][i_q][i_r][i_f]; output_Xsection_delta = Xsection_delta_TSV[i_s][i_q][i_r][i_f]; unit = "fb"; shift = 0;}
	  else if (abs(Xsection) >= 1.e-3){output_Xsection = Xsection_TSV[i_s][i_q][i_r][i_f] * 1.e3; output_Xsection_delta = Xsection_delta_TSV[i_s][i_q][i_r][i_f] * 1.e3; unit = "ab";}// shift = 3;}
	  else if (abs(Xsection) >= 1.e-6){output_Xsection = Xsection_TSV[i_s][i_q][i_r][i_f] * 1.e6; output_Xsection_delta = Xsection_delta_TSV[i_s][i_q][i_r][i_f] * 1.e6; unit = "zb";}// shift = 6;}
	  else if (abs(Xsection) >= 1.e-9){output_Xsection = Xsection_TSV[i_s][i_q][i_r][i_f] * 1.e9; output_Xsection_delta = Xsection_delta_TSV[i_s][i_q][i_r][i_f] * 1.e9; unit = "yb";}// shift = 9;}
	  */
	  
	  double dlog_delta_NLO_LO = log10(abs(delta_NLO_LO));
	  int log_delta_NLO_LO;
	  log_delta_NLO_LO = int(dlog_delta_NLO_LO);
	  if (dlog_delta_NLO_LO < 0.){log_delta_NLO_LO--;}
	  
	  double dlog_Xsection = log10(abs(output_Xsection));
	  int log_Xsection;
	  log_Xsection = int(dlog_Xsection);
	  if (dlog_Xsection < 0.){log_Xsection--;}
	  
	  double dlog_Xsection_delta = log10(abs(output_Xsection_delta));
	  int log_Xsection_delta;
	  log_Xsection_delta = int(dlog_Xsection_delta);
	  if (dlog_Xsection_delta < 0.){log_Xsection_delta--;}
	
	  double dlog_rel = log10(abs(output_Xsection_delta / output_Xsection));
	  int log_rel;
	  log_rel = int(dlog_rel);
	  if (dlog_rel < 0.){log_rel--;}
	  
	  double dlog_rel_mu_R = log10(relative_scale_ren_TSV[i_s][i_r]);
	  int log_rel_mu_R;
	  log_rel_mu_R = int(dlog_rel_mu_R);
	  if (dlog_rel_mu_R < 0.){log_rel_mu_R--;}
	  //	  out_integration << log_rel_mu_R << endl;
	  
	  double dlog_rel_mu_F = log10(relative_scale_fact_TSV[i_s][i_f]);
	  int log_rel_mu_F; 
	  log_rel_mu_F = int(dlog_rel_mu_F);
	  if (dlog_rel_mu_F < 0.){log_rel_mu_F--;}
	  //	  out_integration << log_rel_mu_F << endl;
	  
	  //	if ((i_q == 0) && (i_r == 0) && (i_f == 0)){}
	  out_integration << setw(10) << "";
	  out_integration << setw(11) << "rel_mu_R = " << setprecision(6 + log_rel_mu_R) << setw(8) << showpoint << relative_scale_ren_TSV[i_s][i_r] << "   ";
	  out_integration << setw(11) << "rel_mu_F = " << setprecision(6 + log_rel_mu_F) << setw(8) << showpoint << relative_scale_fact_TSV[i_s][i_f] << "   ";
	  
	  out_integration << "cXS = " << setw(10) << setprecision(7 + log_Xsection + shift) << showpoint << right << output_Xsection << " +- " << setw(10) << setprecision(7 + log_Xsection_delta + shift) << output_Xsection_delta << " " << unit << "   d(cXS) / cXS = " << left << setw(11) << setprecision(6) << abs(output_Xsection_delta / output_Xsection) << "   d(cXS) / XS = " << setw(11) << setprecision(6) << delta_NLO_LO << endl;
	  /*
	  out_integration << "cXS = " << setw(10) << setprecision(7 + log_Xsection + shift) << showpoint << right << output_Xsection << " +- " << setw(10) << setprecision(7 + log_Xsection_delta + shift) << output_Xsection_delta << " " << unit << "   d(cXS) / cXS = " << left << setw(11) << setprecision(6) << abs(Xsection_delta / Xsection) << "   d(cXS) / XS = " << setw(11) << setprecision(6) << delta_NLO_LO << endl;
	  */
	    }
      }
      if (next_no_qTcut < no_qTcut_integration.size() - 1){next_no_qTcut++;}

    }
    out_integration.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_step_result_TSV(phasespace_set & psi){
  Logger logger("observable_set::output_step_result_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
    ofstream out_result;
    out_result.open((filename_result_TSV[i_s]).c_str(), ofstream::out | ofstream::trunc);
    out_result << setw(12) << left << psi_i_gen << endl;
    for (int i_q = 0; i_q < output_n_qTcut; i_q++){

      // Add feature to print only selected values later !!!

      if (active_qTcut == 0){
	out_result << "#   qTcut - independent" << endl;
      }
      else {
	out_result << "#   qTcut[" << setw(3) << i_q << "] = " << setw(15) << setprecision(8) << value_qTcut[i_q] << endl;
      }
      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	  out_result << "  " << setw(3) << i_r << " " << setw(15) << setprecision(8) << relative_scale_ren_TSV[i_s][i_r];
	  out_result << "  " << setw(3) << i_f << " " << setw(15) << setprecision(8) << relative_scale_fact_TSV[i_s][i_f];
	  out_result << "  " << setw(24) << setprecision(15) << unit_factor_calculation * Xsection_TSV[i_s][i_q][i_r][i_f] << setw(24) << setprecision(15) << unit_factor_calculation * Xsection_delta_TSV[i_s][i_q][i_r][i_f];
	  out_result << endl;
	}
      }
    }
    out_result.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

//not used: !!!
void observable_set::output_step_distribution_TSV(phasespace_set & psi){
  Logger logger("observable_set::output_step_distribution_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_distribution;
  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
    if (switch_distribution_TSV[i_s] > 0){
      out_distribution.open((filename_distribution_TSV[i_s]).c_str(), ofstream::out | ofstream::trunc);
      // Replaced 'psi_i_gen' by 'n_ps * psi_i_gen' !!!
      logger << LOG_DEBUG_VERBOSE << "n_ps = " << n_ps << endl;
      logger << LOG_DEBUG_VERBOSE << "psi_i_gen = " << psi_i_gen << endl;
      logger << LOG_DEBUG_VERBOSE << "n_ps * psi_i_gen = " << n_ps * psi_i_gen << endl;
      out_distribution << setw(12) << left << n_ps * psi_i_gen << endl;
      //      out_distribution << setw(12) << left << psi_i_gen << endl;
      //      out_distribution << setw(12) << left << psi_i_gen << endl;
      
      for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	//	if (active_qTcut > 0 || i_q == 0){
	if (active_qTcut == 0){out_distribution << "#   no qTcut - independent" << endl;}
	else {out_distribution << "#   qTcut[" << setw(3) << i_q << "] = " << setw(15) << setprecision(8) << value_qTcut[i_q] << endl;}
	for (int i_d = 0; i_d < dat.size(); i_d++){
	  out_distribution << "#   dat.name[" << setw(3) << i_d << "] = " << dat[i_d].xdistribution_name << endl;
	  for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	    for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	      out_distribution << "#   ";
	      out_distribution << "relative_scale_ren[" << setw(3) << i_r << "] = " << setw(15) << setprecision(8) << relative_scale_ren_TSV[i_s][i_r];
	      out_distribution << "   relative_scale_fact[" << setw(3) << i_f << "] = " << setw(15) << setprecision(8) << relative_scale_fact_TSV[i_s][i_f];
	      out_distribution << endl;
	      for (int i_b = 0; i_b < dat[i_d].n_bins; i_b++){
		out_distribution << setw(12) << bin_count_TSV[i_q][i_d][i_b] << "  ";
		out_distribution << setw(24) << setprecision(15) << unit_factor_calculation * bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][i_b] << "  ";
		out_distribution << setw(24) << setprecision(15) << unit_factor2_calculation * bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][i_b];
		out_distribution << endl;
	      }
	    }
	  }
	}
      }
      out_distribution.close();
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::output_step_integration(phasespace_set & psi){
  Logger logger("observable_set::output_step_integration");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_integration;
  out_integration.open(filename_integration.c_str(), ofstream::out | ofstream::app);
  //  static string stars = "***********************************************";

  //  output_integration_one_line(out_integration, psi_i_gen, psi_i_rej, psi_i_acc, psi_i_nan, Xsection, Xsection_delta, 0, oset);
  double delta_NLO_LO = abs(Xsection_delta / sigma_normalization);
  //  int log_Xsection_delta = int(log10(Xsection_delta));
  int log_rel = int(log10(abs(Xsection_delta / Xsection)));
  //  if (d == 0){
  /*
  //  old version:
  out_integration << setw(12) << psi_i_gen << " (" << setw(8) << psi_i_nan << ") " << " ";
  out_integration << setw(10) << psi_i_rej << " ";
  out_integration << setw(10) << psi_i_acc << " ";
  */
  out_integration << setw(12) << psi_i_gen << " [" << setw(3) << psi_i_nan << "] " << " ";
  out_integration << setw(10) << psi_i_acc << " (" << setw(7) << psi_i_tec << ") " << " ";
  out_integration << setw(10) << psi_i_rej << " ";
  //  }
  //  else {out_integration << setw(57);}
  out_integration << "  sigma = " << setw(14) << setprecision(8) << showpoint << unit_factor_calculation * Xsection << " +- " << setw(14) << setprecision(8 + log_rel) << unit_factor_calculation * Xsection_delta << " " << unit_calculation << "    " << setw(11) << setprecision(7 + log_rel) << abs(Xsection_delta / Xsection);
  out_integration << "   " << setw(11) << setprecision(6) << delta_NLO_LO << endl;
  /*
  if (int_end == 1){
    //    cout << "psi_n_events_min = " << psi_n_events_min << endl;
    //    cout << "psi_n_events_max = " << psi_n_events_max << endl;
    //    cout << "psi_i_acc = " << psi_i_acc << endl;
    //    cout << "abs(Xsection_delta / sigma_normalization) = " << abs(Xsection_delta / sigma_normalization) << endl;
    //    cout << "sigma_normalization_deviation = " << sigma_normalization_deviation << endl;
    out_integration << stars << " final result " << stars << endl;
  }
  */
  out_integration.close();

  /*
  // !!!
  static int initialization = 1;
  static vector<double> rel_scale_factor_ren_CV(n_scales_CV);
  static vector<double> rel_scale_factor_fact_CV(n_scales_CV);
  if (initialization == 1){
     for (int s = 0; s < n_scales_CV; s++){
      if (n_scales_CV > 1){
	if (variation_mu_ren_CV == 0){rel_scale_factor_ren_CV[s] = 1.;}
	else if (variation_mu_ren_CV == 1){rel_scale_factor_ren_CV[s] = pow(10., log10(double(variation_factor_CV)) * double(2. * s / (n_scales_CV - 1) - 1));}
	else if (variation_mu_ren_CV == -1){rel_scale_factor_ren_CV[s] = pow(10., -log10(double(variation_factor_CV)) * double(2. * s / (n_scales_CV - 1) - 1));}
	if (variation_mu_fact_CV == 0){rel_scale_factor_fact_CV[s] = 1.;}
	else if (variation_mu_fact_CV == 1){rel_scale_factor_fact_CV[s] = pow(10., log10(double(variation_factor_CV)) * double(2. * s / (n_scales_CV - 1) - 1));}
	else if (variation_mu_fact_CV == -1){rel_scale_factor_fact_CV[s] = pow(10., -log10(double(variation_factor_CV)) * double(2. * s / (n_scales_CV - 1) - 1));}
      }
      else {
	rel_scale_factor_ren_CV[s] = 1.;
	rel_scale_factor_fact_CV[s] = 1.;
      }
    }
    initialization = 0;
  }
  */

  if (switch_CV){

    out_integration.open((filename_integration_CV).c_str(), ofstream::out | ofstream::app);
    for (int i_c = 0; i_c < output_n_qTcut; i_c++){

      /*
    int next_no_qTcut = 0;
    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
      if (i_q != no_qTcut_distribution[next_no_qTcut]){continue;}
      */

      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << endl;
	logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
	logger << LOG_DEBUG_VERBOSE << "Xsection_CV.size() = " << Xsection_CV.size() << endl;
	logger << LOG_DEBUG_VERBOSE << "Xsection_delta_CV.size() = " << Xsection_delta_CV.size() << endl;
	logger << LOG_DEBUG_VERBOSE << "Xsection_CV[" << i_c << "].size() = " << Xsection_CV[i_c].size() << endl;
	logger << LOG_DEBUG_VERBOSE << "Xsection_delta_CV[" << i_c << "].size() = " << Xsection_delta_CV[i_c].size() << endl;
	logger << LOG_DEBUG_VERBOSE << "temp_n_step = " << temp_n_step << endl;
	logger << LOG_DEBUG_VERBOSE << "(psi_i_acc % temp_n_step) = " << (psi_i_acc % temp_n_step) << endl;
	logger << LOG_DEBUG_VERBOSE << "psi_i_acc == psi_n_events_max = " << (psi_i_acc == psi_n_events_max) << endl;
	logger << LOG_DEBUG_VERBOSE << "int_end == 1 = " << (int_end == 1) << endl;
	if ((psi_i_acc % temp_n_step) == 0 || psi_i_acc == psi_n_events_max || int_end == 1){
	  //	  output_integration_single_CV(out_integration, psi_i_gen, psi_i_rej, psi_i_acc, psi_i_nan, Xsection, Xsection_delta, 0, i_c, i_s, oset);
	  double delta_NLO_LO = abs(Xsection_delta_CV[i_c][i_s] / sigma_normalization);
	  //	  int log_Xsection_delta = int(log10(Xsection_delta_CV[i_c][i_s]));
	  int log_rel = int(log10(abs(Xsection_delta_CV[i_c][i_s] / Xsection_CV[i_c][i_s])));
	  if (i_c == 0 && i_s == 0){
	    out_integration << "n(gen) = " << setw(12) << psi_i_gen << "   n(acc) = " << setw(12) << psi_i_acc << "   n(rej) = " << setw(12) << psi_i_rej << "   n(nan) = " << setw(12) << psi_i_nan << "   n(tec) = " << setw(12) << psi_i_tec << endl;
	/*
	    out_integration << setw(8) << psi_i_gen << " (" << setw(5) << psi_i_nan << ") " << " ";
	    out_integration << setw(8) << psi_i_rej << " ";
	    out_integration << setw(8) << psi_i_acc << " ";
	    out_integration << endl;
	*/
	  }
	  if (i_s == 0){
	  logger << LOG_DEBUG_VERBOSE << "value_qTcut.size() = " << value_qTcut.size() << endl;
	    if (active_qTcut == 1){out_integration << setw(6) << "qTcut = " << setprecision(3) << setw(6) << showpoint << value_qTcut[i_c] << " GeV   ";}
	    else {out_integration << setw(6) << "qTcut - independent  ";}
	  }
	  else {out_integration << setw(21) << "";}
	  
	  if (dynamic_scale_CV == 0){
	    out_integration << setw(9) << "mu_ren = " << setprecision(5) << setw(8) << showpoint << mu_ren_CV[1][i_s] << " GeV   " << setw(10) << "mu_fact = " << setprecision(5) << setw(8) << showpoint << mu_fact_CV[1][i_s] << " GeV ";
	  }
	  else {
	    out_integration << setw(9) << "mu_ren = " << setprecision(5) << setw(8) << showpoint << rel_scale_factor_ren_CV[i_s] << " DS" << dynamic_scale_CV << setw(10) << "    mu_fact = " << setprecision(5) << setw(8) << showpoint << rel_scale_factor_fact_CV[i_s] << " DS" << dynamic_scale_CV;
	  }

	  out_integration << "  sigma = " << setw(15) << setprecision(10) << showpoint << unit_factor_calculation * Xsection_CV[i_c][i_s] << " +- " << setw(15) << setprecision(10 + log_rel) << unit_factor_calculation * Xsection_delta_CV[i_c][i_s] << " " << unit_calculation << "    " << setw(11) << setprecision(7 + log_rel) << abs(Xsection_delta_CV[i_c][i_s] / Xsection_CV[i_c][i_s]) << "   " << setw(11) << setprecision(6) << delta_NLO_LO << endl;

	  if ((i_c == output_n_qTcut - 1) && (i_s == n_scales_CV - 1)){out_integration << endl;}
	}
      }
    }
    //    if (int_end == 1){out_integration << stars << " final result " << stars << endl;}
    out_integration.close();
  }

  //  for (int i_s = 0; i_s < n_set_TSV; i_s++){output_integration_TSV(i_gen, i_rej, psi_i_acc, i_nan, oset, i_s);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_step_result(phasespace_set & psi){
  Logger logger("observable_set::output_step_result");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_result;
  out_result.open(filename_result.c_str(), ofstream::out | ofstream::trunc);
  out_result << psi_i_gen << endl << setprecision(15) << unit_factor_calculation * Xsection << endl << setprecision(15) << unit_factor_calculation * Xsection_delta << endl;
  if (switch_CV){
    //    for (int i_c = 0; i_c < n_qTcut; i_c++){
    for (int i_c = 0; i_c < output_n_qTcut; i_c++){
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	out_result << setprecision(15) << unit_factor_calculation * Xsection_CV[i_c][i_s] << endl << setprecision(15) << unit_factor_calculation * Xsection_delta_CV[i_c][i_s] << endl;
      }
    }
  }

  long long temp_time = 3600 * h + 60 * min + sec;
  out_result << temp_time << endl;
  out_result << setprecision(15) << temp_time / double(psi_i_acc) << endl;
  out_result << setprecision(15) << pow(Xsection_delta, 2) * temp_time << endl;
  out_result << setprecision(15) << pow(Xsection_delta, 2) * double(psi_i_acc) << endl;

  out_result.close();

  //  for (int i_s = 0; i_s < n_set_TSV; i_s++){output_integration_TSV(i_gen, i_rej, psi_i_acc, i_nan, oset, i_s);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_step_errorplot(phasespace_set & psi){
  Logger logger("observable_set::output_step_errorplot");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_gnuplot;
  out_gnuplot.open(filename_gnuplot.c_str(), ofstream::out | ofstream::app);  
  out_gnuplot << setw(10) << psi_i_acc << setw(25) << Xsection << setw(25) << Xsection_delta << setw(25) << abs(Xsection_delta / sigma_normalization) << endl;
  //  if (psi_i_acc == temp_n_step || psi_i_acc % (10 * temp_n_step) == 0 || psi_i_acc == psi_n_events_max || int_end == 1){
  double delta = abs(Xsection_delta / sigma_normalization);
  output_gnuplotfile(filename_makegnuplot, psi_i_acc, delta, psi);
    //  }
  out_gnuplot.close();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_step_execution(phasespace_set & psi){
  Logger logger("observable_set::output_step_execution");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_execution;
  if (*psi.i_step_mode == psi_n_step){
    //  else if (psi_i_acc == psi_n_step){
    out_execution.open(filename_execution.c_str(), ofstream::out | ofstream::trunc);
    out_execution << "started" << endl;
    out_execution.close();
  }

  if (int_end == 1){
    out_execution.open(filename_execution.c_str(), ofstream::out | ofstream::trunc);
    out_execution << "final result" << endl;
    out_execution.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_momenta(ofstream & out_comparison){
  Logger logger("observable_set::output_momenta");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int ib = 1; ib < p_parton[0].size(); ib++){
    out_comparison << "   p[" << setw(3) << csi->type_parton[0][ib] << "] = (" <<
      setprecision(15) << setw(19) << p_parton[0][ib].x0() << "; " <<
      setprecision(15) << setw(19) << p_parton[0][ib].x1() << ", " <<
      setprecision(15) << setw(19) << p_parton[0][ib].x2() << ", " <<
      setprecision(15) << setw(19) << p_parton[0][ib].x3() << ")" <<
      "   " << sqrt(abs(p_parton[0][ib].m2())) << endl;
  }
  out_comparison << endl;
  fourvector check;
  for (int ib = 1; ib < p_parton[0].size(); ib++){
    if (ib < 3){check = check + p_parton[0][ib];}
    else{check = check - p_parton[0][ib];}
  }
  out_comparison << "momentum conservation check = (" <<
    setprecision(15) << setw(19) << check.x0() << "; " <<
    setprecision(15) << setw(19) << check.x1() << ", " <<
    setprecision(15) << setw(19) << check.x2() << ", " <<
    setprecision(15) << setw(19) << check.x3() << ")" <<
    "   " << sqrt(abs(check.m2())) << endl;
  out_comparison << "sqrt(abs((p(in) - p(out))^2 / p(in)^2)) = " << sqrt(abs(check.m2() / p_parton[0][0].m2())) << endl;
  out_comparison << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_momenta_phasespace(ofstream & out_comparison, int x_a){
  Logger logger("observable_set::output_momenta");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int ib = 1; ib < p_parton[x_a].size(); ib++){
    out_comparison << "   p[" << setw(3) << csi->type_parton[x_a][ib] << "] = (" <<
      setprecision(15) << setw(19) << p_parton[x_a][ib].x0() << "; " <<
      setprecision(15) << setw(19) << p_parton[x_a][ib].x1() << ", " <<
      setprecision(15) << setw(19) << p_parton[x_a][ib].x2() << ", " <<
      setprecision(15) << setw(19) << p_parton[x_a][ib].x3() << ")" <<
      "   " << sqrt(abs(p_parton[x_a][ib].m2())) << endl;
  }
  out_comparison << endl;
  fourvector check;
  for (int ib = 1; ib < p_parton[x_a].size(); ib++){
    if (ib < 3){check = check + p_parton[x_a][ib];}
    else{check = check - p_parton[x_a][ib];}
  }
  out_comparison << "momentum conservation check = (" <<
    setprecision(15) << setw(19) << check.x0() << "; " <<
    setprecision(15) << setw(19) << check.x1() << ", " <<
    setprecision(15) << setw(19) << check.x2() << ", " <<
    setprecision(15) << setw(19) << check.x3() << ")" <<
    "   " << sqrt(abs(check.m2())) << endl;
  out_comparison << "sqrt(abs((p(in) - p(out))^2 / p(in)^2)) = " << sqrt(abs(check.m2() / p_parton[x_a][0].m2())) << endl;
  out_comparison << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::output_distribution_complete(phasespace_set & psi){
  Logger logger("observable_set::output_distribution_complete");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG_VERBOSE << "switch_distribution = " << switch_distribution << endl;
  if (switch_distribution){
    //    output_distribution(psi);
    //    output_dddistribution(psi);
    logger << LOG_DEBUG_VERBOSE << "switch_CV = " << switch_CV << endl;
    if (switch_CV){
      output_distribution_CV(psi);
      output_dddistribution_CV(psi);
    }
    output_distribution_TSV(psi);
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_distribution(phasespace_set & psi){
  Logger logger("observable_set::output_distribution");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_distribution;
  logger << LOG_DEBUG_VERBOSE << "XXX dat.size() = " << dat.size() << endl;
  for (int i_d = 0; i_d < dat.size(); i_d++){
  logger << LOG_DEBUG_VERBOSE << "XXX filename_distribution[" << i_d << "] = " << filename_distribution[i_d] << endl;
    out_distribution.open(filename_distribution[i_d].c_str(), ofstream::out | ofstream::trunc);
    out_distribution << psi_i_gen << endl;
    for (int i_b = 0; i_b < dat[i_d].n_bins; i_b++){
      out_distribution << bin_counts[i_d][i_b] << endl;
      out_distribution << setprecision(15) << unit_factor_calculation * bin_weight[i_d][i_b] << endl;
      out_distribution << setprecision(15) << unit_factor2_calculation * bin_weight2[i_d][i_b] << endl;
    }
    out_distribution.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_dddistribution(phasespace_set & psi){
  Logger logger("observable_set::output_dddistribution");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_dddistribution;

  for (int i_d = 0; i_d < dddat.size(); i_d++){
  logger << LOG_DEBUG_VERBOSE << "XXX filename_dddistribution[" << i_d << "] = " << filename_dddistribution[i_d] << endl;
    out_dddistribution.open(filename_dddistribution[i_d].c_str(), ofstream::out | ofstream::trunc);
    out_dddistribution << psi_i_gen << endl;
    for (int i_b = 0; i_b < dddat[i_d].n_bins; i_b++){
      out_dddistribution << bin_counts[dat.size() + i_d][i_b] << endl;
      out_dddistribution << setprecision(15) << unit_factor_calculation * bin_weight[dat.size() + i_d][i_b] << endl;
      out_dddistribution << setprecision(15) << unit_factor2_calculation * bin_weight2[dat.size() + i_d][i_b] << endl;
    }
    out_dddistribution.close();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_distribution_CV(phasespace_set & psi){
  Logger logger("observable_set::output_distribution_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_distribution;
  //        logger << LOG_DEBUG << "unit_factor_calculation = " << unit_factor_calculation << endl;
  //        logger << LOG_DEBUG << "unit_factor2_calculation = " << unit_factor2_calculation << endl;

  for (int i_s = 0; i_s < n_scales_CV; i_s++){
    logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
    logger << LOG_DEBUG_VERBOSE << "dat.size() = " << dat.size() << endl;
    for (int i_d = 0; i_d < dat.size(); i_d++){
      logger << LOG_DEBUG_VERBOSE << "filename_distribution_CV[" << i_s << "][" << i_d << "] = " << filename_distribution_CV[i_s][i_d] << endl;
      out_distribution.open(filename_distribution_CV[i_s][i_d].c_str(), ofstream::out | ofstream::trunc);
      //      out_distribution << "# " << dat[i_d].distname() << endl;
      out_distribution << psi_i_gen << endl;
      for (int i_b = 0; i_b < dat[i_d].n_bins; i_b++){
	out_distribution << bin_counts[i_d][i_b] << endl;
	out_distribution << setprecision(15) << unit_factor_calculation * bin_weight_CV[i_s][i_d][i_b] << endl;
	out_distribution << setprecision(15) << unit_factor2_calculation * bin_weight2_CV[i_s][i_d][i_b] << endl;
      }
      out_distribution.close();
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_dddistribution_CV(phasespace_set & psi){
  Logger logger("observable_set::output_dddistribution_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_dddistribution;

  for (int i_s = 0; i_s < n_scales_CV; i_s++){
    logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
    for (int i_d = 0; i_d < dddat.size(); i_d++){
      logger << LOG_DEBUG_VERBOSE << "filename_dddistribution_CV[" << i_s << "][" << i_d << "] = " << filename_dddistribution_CV[i_s][i_d] << endl;
      out_dddistribution.open(filename_dddistribution_CV[i_s][i_d].c_str(), ofstream::out | ofstream::trunc);
      out_dddistribution << psi_i_gen << endl;
      for (int i_b = 0; i_b < dddat[i_d].n_bins; i_b++){
	out_dddistribution << bin_counts[dat.size() + i_d][i_b] << endl;
	out_dddistribution << setprecision(15) << unit_factor_calculation * bin_weight_CV[i_s][dat.size() + i_d][i_b] << endl;
	out_dddistribution << setprecision(15) << unit_factor2_calculation * bin_weight2_CV[i_s][dat.size() + i_d][i_b] << endl;
      }
      out_dddistribution.close();
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_distribution_TSV(phasespace_set & psi){
  Logger logger("observable_set::output_distribution_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_distribution;
  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
    if (switch_distribution_TSV[i_s] > 0){
      out_distribution.open((filename_distribution_TSV[i_s]).c_str(), ofstream::out | ofstream::trunc);
      // Replaced 'psi_i_gen' by 'n_ps * psi_i_gen' !!!
      //      out_distribution << setw(12) << left << n_ps * psi_i_gen << endl;
      out_distribution << setw(12) << left << psi_i_gen << endl;
      out_distribution << setw(12) << left << n_ps << endl;
       int next_no_qTcut = 0;
      for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	// output only for values selected by 'selection_qTcut_distribution':
	if (i_q != no_qTcut_distribution[next_no_qTcut]){continue;}
	logger << LOG_DEBUG << "i_q = " << i_q << "   no_qTcut_distribution[next_no_qTcut = " << next_no_qTcut << "] = " << no_qTcut_distribution[next_no_qTcut] << endl;
	
	  if (active_qTcut == 0){out_distribution << "#   no qTcut" << endl;}
	  else {out_distribution << "#   qTcut[" << setw(3) << i_q << "] = " << setw(15) << setprecision(8) << value_qTcut[i_q] << endl;}
	  for (int i_d = 0; i_d < dat.size(); i_d++){
	    out_distribution << "#   dat.name[" << setw(3) << i_d << "] = " << dat[i_d].xdistribution_name << endl;
	    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		out_distribution << "#   ";
		out_distribution << "relative_scale_ren[" << setw(3) << i_r << "] = " << setw(15) << setprecision(8) << relative_scale_ren_TSV[i_s][i_r];
		out_distribution << "   relative_scale_fact[" << setw(3) << i_f << "] = " << setw(15) << setprecision(8) << relative_scale_fact_TSV[i_s][i_f];
		out_distribution << endl;
		/*
		cout << "bin_count_TSV.size() = " << bin_count_TSV.size() << endl;
		cout << "bin_count_TSV[" << i_q << "].size() = " << bin_count_TSV[i_q].size() << endl;
		cout << "bin_count_TSV[" << i_q << "][" << i_d << "].size() = " << bin_count_TSV[i_q][i_d].size() << endl;
		*/
		for (int i_b = 0; i_b < dat[i_d].n_bins; i_b++){
		  out_distribution << setw(12) << bin_count_TSV[i_q][i_d][i_b] << "  ";
		  out_distribution << setw(24) << setprecision(15) << unit_factor_calculation * bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][i_b] << "  ";
		  out_distribution << setw(24) << setprecision(15) << unit_factor2_calculation * bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][i_b];
		  out_distribution << endl;
		}
	      }
	    }
	    //	  }
	  }

	  for (int i_d = 0; i_d < dddat.size(); i_d++){
	    int i_ddd = dat.size() + i_d;
	    out_distribution << "#   dddat.name[" << setw(3) << i_d << " - " << setw(3) << i_ddd << "] = " << dddat[i_d].name << endl;
	    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		out_distribution << "#   ";
		out_distribution << "relative_scale_ren[" << setw(3) << i_r << "] = " << setw(15) << setprecision(8) << relative_scale_ren_TSV[i_s][i_r];
		out_distribution << "   relative_scale_fact[" << setw(3) << i_f << "] = " << setw(15) << setprecision(8) << relative_scale_fact_TSV[i_s][i_f];
		out_distribution << endl;
		/*
		cout << "bin_count_TSV.size() = " << bin_count_TSV.size() << endl;
		cout << "bin_count_TSV[" << i_q << "].size() = " << bin_count_TSV[i_q].size() << endl;
		cout << "bin_count_TSV[" << i_q << "][" << i_d << "].size() = " << bin_count_TSV[i_q][i_ddd].size() << endl;
		*/
		for (int i_b = 0; i_b < dddat[i_d].n_bins; i_b++){
		  out_distribution << setw(12) << bin_count_TSV[i_q][i_ddd][i_b] << "  ";
		  out_distribution << setw(24) << setprecision(15) << unit_factor_calculation * bin_weight_TSV[i_s][i_q][i_r][i_f][i_ddd][i_b] << "  ";
		  out_distribution << setw(24) << setprecision(15) << unit_factor2_calculation * bin_weight2_TSV[i_s][i_q][i_r][i_f][i_ddd][i_b];
		  out_distribution << endl;
		}
	      }
	    }
	    //	  }
	  }
	  if (next_no_qTcut < no_qTcut_distribution.size() - 1){next_no_qTcut++;}
      }
      out_distribution.close();
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::output_pdf_comparison_CA_to_CT(phasespace_set & psi){
  Logger logger("observable_set::output_pdf_comparison_CA_to_CT");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    for (int ss = 0; ss < value_mu_fact[sd].size(); ss++){
      if (sd != 1 || ss != 1){continue;}
      for (int i_x = 0; i_x < 3; i_x++){
	//      if (i_x != 0){continue;}
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "]      [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] z1x2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_z1x2[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1z2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1z2[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] g_x2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_gx2[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1_g [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1g[sd][ss][i_x] << endl;
	logger.newLine(LOG_DEBUG);

	for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
	  for (int i_z = 0; i_z < 3; i_z++){
	    if (CA_value_pdf_factor[sd][ss][i_a][i_z][i_x] == 0.){continue;}
	    logger << LOG_DEBUG << "CA_value_pdf_factor[" << sd << "][" << ss << "][" << i_a << "][" << i_z << "][" << i_x << "] = " << setw(23) << setprecision(15) << CA_value_pdf_factor[sd][ss][i_a][i_z][i_x] << "   " << setw(15) << (*CA_collinear)[i_a][0].name() << "   " << (*CA_collinear)[i_a][0].type() << endl;
	  }
	}
	logger.newLine(LOG_DEBUG);
      }	
    }
  }

  
  for (int i_x = 0; i_x < 3; i_x++){
     if (i_x != 0){continue;}
    logger << LOG_DEBUG << "pdf_factor[" << i_x << "] = " << pdf_factor[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_z1x2[" << i_x << "] = " << QT_pdf_factor_z1x2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_x1z2[" << i_x << "] = " << QT_pdf_factor_x1z2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_gx2[" << i_x << "] = " << QT_pdf_factor_gx2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_x1g[" << i_x << "] = " << QT_pdf_factor_x1g[i_x] << endl;
    logger.newLine(LOG_DEBUG);
  }/*
    if (switch_CV){
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	logger << LOG_DEBUG << "pdf_factor_CV[" << i_s << "][" << i_x << "] = " << pdf_factor_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_z1x2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_z1x2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_x1z2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_x1z2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_gx2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_gx2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_x1g_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_x1g_CV[i_s][i_x] << endl;
	logger.newLine(LOG_DEBUG);
      }
    }
    logger.newLine(LOG_DEBUG);
  }
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
void observable_set::output_pdf_comparison_CX_to_CT(phasespace_set & psi){
  Logger logger("observable_set::output_pdf_comparison_CA_to_CT");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    for (int ss = 0; ss < value_mu_fact[sd].size(); ss++){
      if (sd != 1 || ss != 1){continue;}
      for (int i_x = 0; i_x < 3; i_x++){
	//      if (i_x != 0){continue;}

	for (int i_e = 0; i_e < multicollinear.size(); i_e++){
	  if (i_e == 0){
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "]      [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor[sd][ss][i_x] << endl;
	  }
	  if (i_e == 1){
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] z1x2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_z1x2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1z2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1z2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] g_x2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_gx2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1_g [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1g[sd][ss][i_x] << endl;
	  }
	  if (i_e == 2){
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] g__g [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_gg[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] q_x2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_qx2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1_q [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1q[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] z1z2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_z1z2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] g_z2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_gz2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] z1_g [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_z1g[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] qbx2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_qbx2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1qb [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1qb[sd][ss][i_x] << endl;
	  }

	  for (int i_c = 0; i_c < multicollinear[i_e].size(); i_c++){
	    if (CX_value_pdf_factor[sd][ss][i_e][i_c][i_x] == 0.){continue;}
	    stringstream temp;
	    temp << "CX_value_pdf_factor[" << sd << "][" << ss << "][" << i_e << "][" << setw(2) << i_c << "][" << i_x << "] = " << setw(23) << setprecision(15) << CX_value_pdf_factor[sd][ss][i_e][i_c][i_x] << "   " << setw(1+ 9 * i_e) << multicollinear[i_e][i_c].name << "  ";
	    for (int i_i = 1; i_i < multicollinear[i_e][i_c].type_splitting.size(); i_i++){temp  << " " << multicollinear[i_e][i_c].type_splitting[i_i];}
	    logger << LOG_DEBUG << temp.str() << endl;
	  }
	  logger.newLine(LOG_DEBUG);
	}
      }	
    }
  }

  
  for (int i_x = 0; i_x < 3; i_x++){
     if (i_x != 0){continue;}
    logger << LOG_DEBUG << "pdf_factor[" << i_x << "] = " << pdf_factor[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_z1x2[" << i_x << "] = " << QT_pdf_factor_z1x2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_x1z2[" << i_x << "] = " << QT_pdf_factor_x1z2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_gx2[" << i_x << "] = " << QT_pdf_factor_gx2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_x1g[" << i_x << "] = " << QT_pdf_factor_x1g[i_x] << endl;
    logger.newLine(LOG_DEBUG);
  }/*
    if (switch_CV){
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	logger << LOG_DEBUG << "pdf_factor_CV[" << i_s << "][" << i_x << "] = " << pdf_factor_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_z1x2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_z1x2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_x1z2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_x1z2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_gx2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_gx2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_x1g_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_x1g_CV[i_s][i_x] << endl;
	logger.newLine(LOG_DEBUG);
      }
    }
    logger.newLine(LOG_DEBUG);
  }
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
void observable_set::output_pdf_comparison_CX_to_CT2(phasespace_set & psi){
  Logger logger("observable_set::output_pdf_comparison_CA_to_CT");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    for (int ss = 0; ss < value_mu_fact[sd].size(); ss++){
      if (sd != 1 || ss != 1){continue;}
      for (int i_x = 0; i_x < 3; i_x++){
	//      if (i_x != 0){continue;}
	//	for (int i_e = 0; i_e < multicollinear.size(); i_e++){
	//	  if (i_e == 0){
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "]      [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor[sd][ss][i_x] << endl;
	//	  }
	//	  if (i_e == 1){
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] z1x2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_z1x2[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1z2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1z2[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] g_x2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_gx2[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1_g [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1g[sd][ss][i_x] << endl;
	//	  }
	//	  if (i_e == 2){
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] g__g [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_gg[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] q_x2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_qx2[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1_q [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1q[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] z1z2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_z1z2[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] g_z2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_gz2[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] z1_g [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_z1g[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] qbx2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_qbx2[sd][ss][i_x] << endl;
	logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1qb [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1qb[sd][ss][i_x] << endl;
	//	  }
	
	// needs to be adapted for ncollinear instead of multicollinear !!!
	/*
	  for (int i_c = 0; i_c < multicollinear[i_e].size(); i_c++){
	    if (CX_value_pdf_factor[sd][ss][i_e][i_c][i_x] == 0.){continue;}
	    stringstream temp;
	    temp << "CX_value_pdf_factor[" << sd << "][" << ss << "][" << i_e << "][" << setw(2) << i_c << "][" << i_x << "] = " << setw(23) << setprecision(15) << CX_value_pdf_factor[sd][ss][i_e][i_c][i_x] << "   " << setw(1+ 9 * i_e) << multicollinear[i_e][i_c].name << "  ";
	    for (int i_i = 1; i_i < multicollinear[i_e][i_c].type_splitting.size(); i_i++){temp  << " " << multicollinear[i_e][i_c].type_splitting[i_i];}
	    logger << LOG_DEBUG << temp.str() << endl;
	  }
	*/
	for (int i_c = 0; i_c < ncollinear.size(); i_c++){
	  if (CX_value_pdf_factor[sd][ss][0][i_c][i_x] == 0.){continue;}
	  stringstream temp;
	  temp << "CX_value_pdf_factor[" << sd << "][" << ss << "][" << 0 << "][" << setw(2) << i_c << "][" << i_x << "] = " << setw(23) << setprecision(15) << CX_value_pdf_factor[sd][ss][0][i_c][i_x] << "   " << setw(15) << ncollinear[i_c].name << "  ";
	  for (int i_i = 1; i_i < ncollinear[i_c].type_splitting.size(); i_i++){temp  << " " << ncollinear[i_c].type_splitting[i_i];}
	  logger << LOG_DEBUG << temp.str() << endl;
	}
	logger.newLine(LOG_DEBUG);
      }	
    }
  }

  
  for (int i_x = 0; i_x < 3; i_x++){
     if (i_x != 0){continue;}
    logger << LOG_DEBUG << "pdf_factor[" << i_x << "] = " << pdf_factor[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_z1x2[" << i_x << "] = " << QT_pdf_factor_z1x2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_x1z2[" << i_x << "] = " << QT_pdf_factor_x1z2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_gx2[" << i_x << "] = " << QT_pdf_factor_gx2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_x1g[" << i_x << "] = " << QT_pdf_factor_x1g[i_x] << endl;
    logger.newLine(LOG_DEBUG);
  }
  /*
    if (switch_CV){
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	logger << LOG_DEBUG << "pdf_factor_CV[" << i_s << "][" << i_x << "] = " << pdf_factor_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_z1x2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_z1x2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_x1z2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_x1z2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_gx2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_gx2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_x1g_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_x1g_CV[i_s][i_x] << endl;
	logger.newLine(LOG_DEBUG);
      }
    }
    logger.newLine(LOG_DEBUG);
  }
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::output_pdf_comparison_CX_multicollinear_to_CT2(phasespace_set & psi){
  Logger logger("observable_set::output_pdf_comparison_CA_to_CT");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    for (int ss = 0; ss < value_mu_fact[sd].size(); ss++){
      if (sd != 1 || ss != 1){continue;}
      for (int i_x = 0; i_x < 3; i_x++){
	//      if (i_x != 0){continue;}

	for (int i_e = 0; i_e < multicollinear.size(); i_e++){
	  if (i_e == 0){
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "]      [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor[sd][ss][i_x] << endl;
	  }
	  if (i_e == 1){
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] z1x2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_z1x2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1z2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1z2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] g_x2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_gx2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1_g [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1g[sd][ss][i_x] << endl;
	  }
	  if (i_e == 2){
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] g__g [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_gg[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] q_x2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_qx2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1_q [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1q[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] z1z2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_z1z2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] g_z2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_gz2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] z1_g [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_z1g[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] qbx2 [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_qbx2[sd][ss][i_x] << endl;
	    logger << LOG_DEBUG << "QT_value_pdf_factor[" << sd << "][" << ss << "] x1qb [" << i_x << "] = " << setw(23) << setprecision(15) << QT_value_pdf_factor_x1qb[sd][ss][i_x] << endl;
	  }

	  for (int i_c = 0; i_c < multicollinear[i_e].size(); i_c++){
	    if (CX_value_pdf_factor[sd][ss][i_e][i_c][i_x] == 0.){continue;}
	    stringstream temp;
	    temp << "CX_value_pdf_factor[" << sd << "][" << ss << "][" << i_e << "][" << setw(2) << i_c << "][" << i_x << "] = " << setw(23) << setprecision(15) << CX_value_pdf_factor[sd][ss][i_e][i_c][i_x] << "   " << setw(1+ 9 * i_e) << multicollinear[i_e][i_c].name << "  ";
	    for (int i_i = 1; i_i < multicollinear[i_e][i_c].type_splitting.size(); i_i++){temp  << " " << multicollinear[i_e][i_c].type_splitting[i_i];}
	    logger << LOG_DEBUG << temp.str() << endl;
	  }

	  logger.newLine(LOG_DEBUG);
	}
      }	
    }
  }

  
  for (int i_x = 0; i_x < 3; i_x++){
     if (i_x != 0){continue;}
    logger << LOG_DEBUG << "pdf_factor[" << i_x << "] = " << pdf_factor[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_z1x2[" << i_x << "] = " << QT_pdf_factor_z1x2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_x1z2[" << i_x << "] = " << QT_pdf_factor_x1z2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_gx2[" << i_x << "] = " << QT_pdf_factor_gx2[i_x] << endl;
    logger << LOG_DEBUG << "QT_pdf_factor_x1g[" << i_x << "] = " << QT_pdf_factor_x1g[i_x] << endl;
    logger.newLine(LOG_DEBUG);
  }/*
    if (switch_CV){
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	logger << LOG_DEBUG << "pdf_factor_CV[" << i_s << "][" << i_x << "] = " << pdf_factor_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_z1x2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_z1x2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_x1z2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_x1z2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_gx2_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_gx2_CV[i_s][i_x] << endl;
	logger << LOG_DEBUG << "QT_pdf_factor_x1g_CV[" << i_s << "][" << i_x << "] = " << QT_pdf_factor_x1g_CV[i_s][i_x] << endl;
	logger.newLine(LOG_DEBUG);
      }
    }
    logger.newLine(LOG_DEBUG);
  }
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
