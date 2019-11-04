#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

////////////////////
//  constructors  //
////////////////////
observable_set::observable_set(){
  int_end = 0;
  max_dyn_ren = 0;
  max_dyn_fact = 0;
}



observable_set::observable_set(int _n_set_TSV){
  n_set_TSV = _n_set_TSV;
  int_end = 0;  
  max_dyn_ren = 0;
  max_dyn_fact = 0;
}



observable_set::observable_set(inputparameter_set & isi, contribution_set & _csi){
  Logger logger("observable_set::observable_set");
  logger << LOG_DEBUG << "called" << endl;

  csi = &_csi;

  logger << LOG_INFO << setw(35) << "isi.switch_TSV" << " = " << isi.switch_TSV << endl;
  int_end = 0;
  max_dyn_ren = 0;
  max_dyn_fact = 0;

  //  run_mode = isi.run_mode;

  switch_reference = isi.switch_reference;

  // fake parameters which can be removed later:
  //  int pdf_set = 0;

  /////////////////////////////////////////////////////////
  //  set observable_set values from inputparameter_set  //
  /////////////////////////////////////////////////////////


  ////////////////////////////
  //  contribution_set csi  //
  ////////////////////////////

  //  csi = isi.csi;

  /////////////////////////
  //  beam / CMS energy  //
  /////////////////////////

  E_beam = isi.E;
  E_CMS = 2 * E_beam;

  //////////////////////////////////////////////////////////////
  //  selection of process and contribution to be calculated  //
  //////////////////////////////////////////////////////////////

  process_class = isi.csi.process_class;
  subprocess = isi.csi.subprocess;
  type_perturbative_order = isi.csi.type_perturbative_order;
  type_contribution = isi.csi.type_contribution;
  type_correction = isi.csi.type_correction;
  contribution_order_alpha_s = isi.csi.contribution_order_alpha_s;
  contribution_order_alpha_e = isi.csi.contribution_order_alpha_e;
  contribution_order_interference = isi.csi.contribution_order_interference;


  /*
  assert(subprocess == csi->subprocess && "subprocess == csi->subprocess");
  assert(type_perturbative_order == csi->type_perturbative_order && "type_perturbative_order == csi->type_perturbative_order");
  assert(type_contribution == csi->type_contribution && "type_contribution == csi->type_contribution");
  assert(type_correction == csi->type_correction && "type_correction == csi->type_correction");
  assert(contribution_order_alpha_s == csi->contribution_order_alpha_s && "contribution_order_alpha_s == csi->contribution_order_alpha_s");
  assert(contribution_order_alpha_e == csi->contribution_order_alpha_e && "contribution_order_alpha_e == csi->contribution_order_alpha_e");
  assert(contribution_order_interference == csi->contribution_order_interference && "contribution_order_interference == csi->contribution_order_interference");
  */
  /*
  if (csi->type_contribution == "born"){
    order_alphas_born = csi->contribution_order_alpha_s;
  } 
  else if (csi->type_correction == "QCD" && 
	   (csi->type_contribution == "RA" ||
	    csi->type_contribution == "VA" ||
	    csi->type_contribution == "CA" ||
	    csi->type_contribution == "RT" ||
	    csi->type_contribution == "VT" ||
	    csi->type_contribution == "CT")){
    order_alphas_born = csi->contribution_order_alpha_s - 1;
  }
  else if (csi->type_correction == "QCD" && 
	   (csi->type_contribution == "VT2" ||
	    csi->type_contribution == "CT2" ||
	    csi->type_contribution == "RVA" ||
	    csi->type_contribution == "RCA" ||
	    csi->type_contribution == "RRA")){
    order_alphas_born = csi->contribution_order_alpha_s - 2;
  }
  else if (//csi->type_correction == "QCD" && 
	   (csi->type_contribution == "loop" ||
	    csi->type_contribution == "LI2")){
    order_alphas_born = csi->contribution_order_alpha_s - 2;
  }
  else if (//csi->type_correction == "QCD" && 
	   (csi->type_contribution == "L2RA" ||
	    csi->type_contribution == "L2VA" ||
	    csi->type_contribution == "L2CA" ||
	    csi->type_contribution == "L2RT" ||
	    csi->type_contribution == "L2VT" ||
	    csi->type_contribution == "L2CT")){
    order_alphas_born = csi->contribution_order_alpha_s - 3;
  }
  else {
    order_alphas_born = csi->contribution_order_alpha_s;
  }
  */
  
  order_alphas_born = csi->order_alpha_s_born;

  //  order_alphas_born is probably not very useful for the N-jettiness subtraction.
  //  Replace (or add?) parameter based on number of Born-level jets:


  NJ_q_axes.resize(3 + csi->n_jet_born);
  NJ_n_axes.resize(3 + csi->n_jet_born);
  NJ_Ei.resize(3 + csi->n_jet_born);

  NJ_q_axes_frame.resize(3 + csi->n_jet_born);
  NJ_n_axes_frame.resize(3 + csi->n_jet_born);
  NJ_Qi.resize(3 + csi->n_jet_born);

  /*
  int n_jet_born = 0;
  for (int i_s = 0; i_s < csi->process_class.size(); i_s++){
    if (csi->process_class[i_s] == 'j'){n_jet_born++;}
  }
  */
  /*
  NJ_q_axes.resize(3 + csi->n_jet_born);
  NJ_q_axes[1] = fourvector(1., 0., 0., 1.);
  NJ_q_axes[2] = fourvector(1., 0., 0., -1.);
  */
  // only if needed ???
  
  logger << LOG_INFO << "csi->order_alpha_s_born = " << setw(23) << setprecision(15) << csi->order_alpha_s_born << endl;
  logger << LOG_INFO << "order_alphas_born = " << order_alphas_born << endl;
  logger << LOG_INFO << "csi->n_jet_born (no of Born-level jets) = " << csi->n_jet_born << endl;



  // should finally not be needed any more...
  // replaces parts of "initialization_contribution();"
  // still needed: n_particle, type_parton

  /////////////////////
  //  event_set esi  //
  /////////////////////

  esi = isi.esi;

  /////////////////////
  //  user_set user  //
  /////////////////////

  user = isi.user;

  initialization_path(isi);
  initialization_unit(isi);


  // ->|  exported to inputparameter_set until here !!!


// **************************************************************************
// *                                                                        *
// *  determination of perturbative QCD order of present contribution       *
// *                                                                        *
// **************************************************************************

  initialization_LHAPDF_parameters(isi);
  initialization_LHAPDF();

  logger << LOG_DEBUG << "LHAPDF has been initialized." << endl;

  int model_perturbative_order = 0;
  vector<string> model_contribution_LHAPDFname(isi.max_perturbative_QCD_order + 1);
  vector<int> model_contribution_LHAPDFsubset(isi.max_perturbative_QCD_order + 1);
  for (int i_s = 0; i_s < isi.type_perturbative_order_counter + 1; i_s++){
    if (isi.collection_type_perturbative_order[i_s] == "LO"){
      model_perturbative_order = 0;
      model_contribution_LHAPDFname[0] = isi.contribution_LHAPDFname[i_s];
      model_contribution_LHAPDFsubset[0] = isi.contribution_LHAPDFsubset[i_s];
    }
    if (isi.collection_type_perturbative_order[i_s] == "NLO"){
      model_perturbative_order = 1;
      model_contribution_LHAPDFname[1] = isi.contribution_LHAPDFname[i_s];
      model_contribution_LHAPDFsubset[1] = isi.contribution_LHAPDFsubset[i_s];
    }
    if (isi.collection_type_perturbative_order[i_s] == "NNLO"){
      model_perturbative_order = 2;
      model_contribution_LHAPDFname[2] = isi.contribution_LHAPDFname[i_s];
      model_contribution_LHAPDFsubset[2] = isi.contribution_LHAPDFsubset[i_s];
    }
  }
  for (int i_s = 0; i_s < model_contribution_LHAPDFname.size(); i_s++){
    if (LHAPDFname == model_contribution_LHAPDFname[i_s]){model_perturbative_order = i_s; break;}
  }
  logger << LOG_DEBUG << "model_perturbative_order = " << model_perturbative_order << endl;
  for (int i_s = 0; i_s < model_contribution_LHAPDFname.size(); i_s++){
    logger << LOG_DEBUG << "model_contribution_LHAPDFname[" << i_s << " = " << setw(25) << model_contribution_LHAPDFname[i_s] << "   " << model_contribution_LHAPDFsubset[i_s] << endl;
  }
  // Only needed for model ???




  int switch_alpha_CMS = 0;
  logger << LOG_DEBUG << "switch_alpha_CMS = " << switch_alpha_CMS << endl;
  logger << LOG_DEBUG << "user.switch_map[switch_alpha_CMS] = " << user.switch_map["switch_alpha_CMS"] << endl;
  logger << LOG_DEBUG << "user.switch_value.size() = " << user.switch_value.size() << endl;
  logger << LOG_DEBUG << "user.switch_value[user.switch_map[switch_alpha_CMS]] = " << user.switch_value[user.switch_map["switch_alpha_CMS"]] << endl;

  switch_alpha_CMS = user.switch_value[user.switch_map["switch_alpha_CMS"]];

  int switch_cosw_real = user.switch_value[user.switch_map["switch_cosw_real"]];


  
  //  model_set msi;
  vector<string> file_input;
  //isi.max_perturbative_QCD_order

  logger << LOG_DEBUG_VERBOSE << "isi.present_type_perturbative_order = " << isi.present_type_perturbative_order << endl;

  msi = model_set(file_input, N_f, N_f_active, scale_ren, model_perturbative_order, isi.collection_type_perturbative_order, model_contribution_LHAPDFname, model_contribution_LHAPDFsubset, switch_alpha_CMS, switch_cosw_real);

  //  msi = model_set(file_input, N_f, N_f_active, scale_ren, isi.present_type_perturbative_order, isi.collection_type_perturbative_order, isi.contribution_LHAPDFname, isi.contribution_LHAPDFsubset, switch_alpha_CMS, switch_cosw_real);
  //  msi = model_set(file_input, N_f, N_f_active, isi.max_perturbative_QCD_order, scale_ren, isi.present_type_perturbative_order, isi.collection_type_perturbative_order, isi.contribution_LHAPDFname, isi.contribution_LHAPDFsubset, switch_alpha_CMS, switch_cosw_real);

  //  model = readin_model(N_f, N_f_active, model_perturbative_order, isi.max_perturbative_QCD_order, model_contribution_LHAPDFname, model_contribution_LHAPDFsubset, switch_alpha_CMS, switch_cosw_real);

  //  parameter_model pmod = readin_model(N_f, N_f_active, model_perturbative_order, isi.max_perturbative_QCD_order, model_contribution_LHAPDFname, model_contribution_LHAPDFsubset, switch_alpha_CMS, switch_cosw_real);


  // alpha_e rescaling factor because of use_adapted_ew_coupling:
  logger << LOG_INFO << "msi.ew_scheme = " << msi.ew_scheme << endl;
  logger << LOG_INFO << "msi.use_cms = " << msi.use_cms << endl;
  logger << LOG_INFO << "msi.use_adapted_ew_coupling = " << msi.use_adapted_ew_coupling << endl;
  logger << LOG_INFO << "n_photon_born = " << csi->n_photon_born << endl;
  logger << LOG_INFO << "csi->contribution_order_alpha_e = " << csi->contribution_order_alpha_e << endl;
  logger << LOG_INFO << "msi.alpha_e = " << setw(23) << setprecision(15) << msi.alpha_e << endl;
  logger << LOG_INFO << "msi.alpha_e_0 = " << setw(23) << setprecision(15) << msi.alpha_e_0 << endl;
  logger << LOG_INFO << "msi.alpha_e_Gmu = " << setw(23) << setprecision(15) << msi.alpha_e_Gmu << endl;
  logger << LOG_INFO << "msi.alpha_e_MZ = " << setw(23) << setprecision(15) << msi.alpha_e_MZ << endl;
 
  if (msi.use_adapted_ew_coupling == -1){
    rescaling_factor_alpha_e = 1.;
  }
  else if (msi.use_adapted_ew_coupling == 0){
    // i.e. ew_scheme = 1 or 2:
    logger << LOG_INFO << "msi.alpha_e = " << setw(23) << setprecision(15) << msi.alpha_e << endl;
    logger << LOG_INFO << "msi.alpha_e_0 = " << setw(23) << setprecision(15) << msi.alpha_e_0 << endl;
    logger << LOG_INFO << "csi->n_photon_born = " << csi->n_photon_born << endl;
    rescaling_factor_alpha_e = pow(msi.alpha_e_0 / msi.alpha_e, csi->n_photon_born);
  }
  else if (msi.use_adapted_ew_coupling == 1 || msi.use_adapted_ew_coupling == 2){
      // check if always correct for "MIX" contributions !!!
    logger << LOG_INFO << "msi.alpha_e = " << setw(23) << setprecision(15) << msi.alpha_e << endl;
    if (msi.use_adapted_ew_coupling == 1){logger << LOG_INFO << setw(23) << setprecision(15) << "msi.alpha_e_Gmu = " << msi.alpha_e_Gmu << endl;}
    else if (msi.use_adapted_ew_coupling == 2){logger << LOG_INFO << "msi.alpha_e_MZ = " << setw(23) << setprecision(15) << msi.alpha_e_MZ << endl;}
    logger << LOG_INFO << "csi->n_photon_born = " << csi->n_photon_born << endl;
    if (csi->type_correction == "QEW" || csi->type_correction == "MIX"){
      logger << LOG_INFO << "csi->contribution_order_alpha_e(" << csi->contribution_order_alpha_e << ") - 1 - csi->n_photon_born(" << csi->n_photon_born << ") = " << csi->contribution_order_alpha_e - 1 - csi->n_photon_born << endl;
      if (msi.use_adapted_ew_coupling == 1){rescaling_factor_alpha_e = pow(msi.alpha_e_Gmu / msi.alpha_e, csi->contribution_order_alpha_e - 1 - csi->n_photon_born);}
      else if (msi.use_adapted_ew_coupling == 2){rescaling_factor_alpha_e = pow(msi.alpha_e_MZ / msi.alpha_e, csi->contribution_order_alpha_e - 1 - csi->n_photon_born);}
    }
    else {
      logger << LOG_INFO << "csi->contribution_order_alpha_e(" << csi->contribution_order_alpha_e << ") - csi->n_photon_born(" << csi->n_photon_born << ") = " << csi->contribution_order_alpha_e - csi->n_photon_born << endl;
      if (msi.use_adapted_ew_coupling == 1){rescaling_factor_alpha_e = pow(msi.alpha_e_Gmu / msi.alpha_e, csi->contribution_order_alpha_e - csi->n_photon_born);}
      else if (msi.use_adapted_ew_coupling == 2){rescaling_factor_alpha_e = pow(msi.alpha_e_MZ / msi.alpha_e, csi->contribution_order_alpha_e - csi->n_photon_born);}
    }
  }

  // rescaling_factor_alpha_e needs to be introduced in  NJ  contributions as well !!!

  logger << LOG_INFO << "rescaling_factor_alpha_e = " << setw(23) << setprecision(15) << rescaling_factor_alpha_e << endl;

  
  
  logger << LOG_INFO << "alpha_S at various scales for " << LHAPDFname << ":" << endl;
  logger << LOG_INFO << "alpha_S(M_W = " << setprecision(5) << setw(7) << msi.M_W << " GeV)) = " << setw(23) << setprecision(15) << LHAPDF::alphasPDF(msi.M_W) << endl;
  logger << LOG_INFO << "alpha_S(M_Z = " << setprecision(6) << setw(7) << msi.M_Z << " GeV)) = " << setw(23) << setprecision(15) << LHAPDF::alphasPDF(msi.M_Z) << endl;
  logger << LOG_INFO << "alpha_S(M_t = " << setprecision(4) << setw(7) << msi.M_t << " GeV)) = " << setw(23) << setprecision(15) << LHAPDF::alphasPDF(msi.M_t) << endl;
  if (msi.M_b != 0.){logger << LOG_INFO << "alpha_S(M_b = " << setprecision(3) << setw(7) << msi.M_b << " GeV)) = " << setw(23) << setprecision(15) << LHAPDF::alphasPDF(msi.M_b) << endl;}
  logger << LOG_INFO << "alpha_S(1 TeV)   = " << setw(23) << setprecision(15) << LHAPDF::alphasPDF(1000.) << endl;
  /*
  // check if content of msi and model (and 'pmod') are identical !!! If so, remove pmod, model and parameter.model.cpp !!!
  assert(model.M_W() == msi.M_W);
  assert(model.M_Z() == msi.M_Z);
  assert(model.M_H() == msi.M_H);
  assert(model.M_t() == msi.M_t);
  assert(model.M_b() == msi.M_b);
  assert(model.M_c() == msi.M_c);
  assert(model.M_s() == msi.M_s);
  assert(model.M_u() == msi.M_u);
  assert(model.M_d() == msi.M_d);
  assert(model.M_e() == msi.M_e);
  assert(model.M_mu() == msi.M_mu);
  assert(model.M_tau() == msi.M_tau);
  assert(model.M_ve() == msi.M_ve);
  assert(model.M_vm() == msi.M_vm);
  assert(model.M_vt() == msi.M_vt);
 
  assert(model.M2_W() == msi.M2_W);
  assert(model.M2_Z() == msi.M2_Z);
  assert(model.M2_H() == msi.M2_H);
  assert(model.M2_t() == msi.M2_t);
  assert(model.M2_b() == msi.M2_b);
  assert(model.M2_c() == msi.M2_c);
  assert(model.M2_s() == msi.M2_s);
  assert(model.M2_u() == msi.M2_u);
  assert(model.M2_d() == msi.M2_d);
  assert(model.M2_e() == msi.M2_e);
  assert(model.M2_mu() == msi.M2_mu);
  assert(model.M2_tau() == msi.M2_tau);
  assert(model.M2_ve() == msi.M2_ve);
  assert(model.M2_vm() == msi.M2_vm);
  assert(model.M2_vt() == msi.M2_vt);
 
  assert(model.Gamma_W() == msi.Gamma_W);
  assert(model.Gamma_Z() == msi.Gamma_Z);
  assert(model.Gamma_H() == msi.Gamma_H);
  assert(model.Gamma_t() == msi.Gamma_t);
  //  assert(model.Gamma_b() == msi.Gamma_b);
  //  assert(model.Gamma_c() == msi.Gamma_c);

  logger << LOG_INFO << left << setw(30) << "msi.map_Gamma_W" << " = " << msi.map_Gamma_W << " GeV" << endl;
  logger << LOG_INFO << left << setw(30) << "model.map_Gamma_W()" << " = " << model.map_Gamma_W() << " GeV" << endl;
  logger << LOG_INFO << left << setw(30) << "msi.map_Gamma_Z" << " = " << msi.map_Gamma_Z << " GeV" << endl;
  logger << LOG_INFO << left << setw(30) << "model.map_Gamma_Z()" << " = " << model.map_Gamma_Z() << " GeV" << endl;
  logger << LOG_INFO << left << setw(30) << "msi.map_Gamma_H" << " = " << msi.map_Gamma_H << " GeV" << endl;
  logger << LOG_INFO << left << setw(30) << "model.map_Gamma_H()" << " = " << model.map_Gamma_H() << " GeV" << endl;
  logger << LOG_INFO << left << setw(30) << "msi.map_Gamma_t" << " = " << msi.map_Gamma_t << " GeV" << endl;
  logger << LOG_INFO << left << setw(30) << "model.map_Gamma_t()" << " = " << model.map_Gamma_t() << " GeV" << endl;

  
  //  assert(model.map_Gamma_W() == msi.map_Gamma_W);
  //  assert(model.map_Gamma_Z() == msi.map_Gamma_Z);
  //  assert(model.map_Gamma_H() == msi.map_Gamma_H);
  //  assert(model.map_Gamma_t() == msi.map_Gamma_t);
  //  //  assert(model.map_Gamma_b() == msi.map_Gamma_b);
  //  //  assert(model.map_Gamma_c() == msi.map_Gamma_c);

  //  assert(model.reg_Gamma_W() == msi.reg_Gamma_W);
  //  assert(model.reg_Gamma_Z() == msi.reg_Gamma_Z);
  //  assert(model.reg_Gamma_H() == msi.reg_Gamma_H);
  //  assert(model.reg_Gamma_t() == msi.reg_Gamma_t);
  //  //  assert(model.reg_Gamma_b() == msi.reg_Gamma_b);
  //  //  assert(model.reg_Gamma_c() == msi.reg_Gamma_c);
  

  assert(model.cM_W() == msi.cM_W);
  assert(model.cM_Z() == msi.cM_Z);
  assert(model.cM_H() == msi.cM_H);
  assert(model.cM_t() == msi.cM_t);
  //  assert(model.cM_b() == msi.cM_b);
  //  assert(model.cM_c() == msi.cM_c);

  assert(model.cM2_W() == msi.cM2_W);
  assert(model.cM2_Z() == msi.cM2_Z);
  assert(model.cM2_H() == msi.cM2_H);
  assert(model.cM2_t() == msi.cM2_t);
  //  assert(model.cM2_b() == msi.cM2_b);
  //  assert(model.cM2_c() == msi.cM2_c);

  assert(model.G_F() == msi.G_F);

  assert(model.Gamma_Wlv() == msi.Gamma_Wlv);
  assert(model.Gamma_Wud() == msi.Gamma_Wud);
  assert(model.Gamma_Zll() == msi.Gamma_Zll);
  assert(model.Gamma_Zvv() == msi.Gamma_Zvv);
  assert(model.Gamma_Zuu() == msi.Gamma_Zuu);
  assert(model.Gamma_Zdd() == msi.Gamma_Zdd);

  assert(model.BR_Wlv() == msi.BR_Wlv);
  assert(model.BR_Wud() == msi.BR_Wud);
  assert(model.BR_Zll() == msi.BR_Zll);
  assert(model.BR_Zvv() == msi.BR_Zvv);
  assert(model.BR_Zuu() == msi.BR_Zuu);
  assert(model.BR_Zdd() == msi.BR_Zdd);

  //  assert(model.alpha_S() == msi.alpha_s);
  //  assert(model.g_S() == msi.g_s);
  assert(model.alpha_e() == msi.alpha_e);
  assert(model.e() == msi.e);

  assert(model.e_pow() == msi.e_pow);
  assert(model.V_ckm() == msi.V_ckm);
  assert(model.M() == msi.M);
  assert(model.M2() == msi.M2);
  assert(model.Gamma() == msi.Gamma);
  assert(model.cM() == msi.cM);
  assert(model.cM2() == msi.cM2);
  assert(model.map_Gamma() == msi.map_Gamma);
  assert(model.reg_Gamma() == msi.reg_Gamma);
  */


  // process-dependent part (type_parton, process_type)
  // process_type = 0: no process
  // process_type = 1: particles decay
  // process_type = 2: collision of 2 particles
  //  vector<vector<int> > type_parton(1);

  logger << LOG_DEBUG << "process_class = " << process_class << endl;
  logger << LOG_DEBUG << "isi.csi.subprocess = " << isi.csi.subprocess << endl;

  //  type_parton = isi.csi.type_parton;
  n_particle = isi.csi.n_particle;
  process_type = isi.csi.process_type;

  /*
  logger << LOG_DEBUG << "before" << endl;
  for (int i_x = 0; i_x < csi.type_parton.size(); i_x++){
    for (int i_p = 0; i_p < csi.type_parton[i_x].size(); i_p++){
      logger << LOG_DEBUG << "csi.type_parton[" << i_x << "][" << i_p << "] = " << csi.type_parton[i_x][i_p] << endl;
    }
  }

  type_parton.resize(1);

//  subprocess_readin(csi.process_class, csi.subprocess, csi.decay, process_type, n_particle, type_parton);
  */

  logger << LOG_DEBUG << "after" << endl;
  for (int i_x = 0; i_x < csi->type_parton.size(); i_x++){
    for (int i_p = 0; i_p < csi->type_parton[i_x].size(); i_p++){
      logger << LOG_DEBUG << "csi->type_parton[" << i_x << "][" << i_p << "] = " << csi->type_parton[i_x][i_p] << endl;
    }
  }

  initialization_object_event_selection(isi);
  //  fills (in particular):
  //  photon_recombination_list -> contains PDG labels of partons that will enter a jet algorithm
  //  jet_algorithm_list        -> contains PDG labels of partonsthat might be 'dressed' with photons 
  //  (photon_recombination_selection, photon_recombination_disable, jet_algorithm_selection, jet_algorithm_disable   are NOT used any longer afterwards.)
  //  sets all basic parameters for jet algorithm, photon recombination and Frixione isolation.
  //  nothing process-specific so far.

  determine_object_definition();
  //  initialization_object_process();
  //  determine_object_definition();
  //  define_pT/ET/eta/y and n_observed_min/max are adapted between different objects (like lep and e/mu/tau etc.).
  //  check if this works correctly !!!

  esi.determine_n_partonlevel(csi->type_parton);
  //  fill n_partonlevel: uses process information !!!
  //  very trivial right now, no use of jet_algorithm_list, etc., no splitting between different phasespaces.
  // is this actually needed here ??? More sophisticated determination happens later !!!
  //  initialization_object_process();


  isi.esi = esi;  // !!! Should be done the other way round !!!


  
  initialization_masses(msi.M, msi.M2);

  logger << LOG_DEBUG << "switch_distribution = " << isi.switch_distribution << endl;

  // this is needed because initialization_TSV is (in the routines functions) called before initialization_CV, but seems to be responsible for computing the central alphaS value. if we want this at the scales multiplied by the prefactor, initialization_generic needs to know about it..
  //  prefactor_CV = isi.prefactor_CV;
  // shifted elsewhere...

  initialization_switches(isi);
  initialization_qTcut(isi);
  initialization_integration_parameters(isi);
  initialization_basic_CV(isi);
  initialization_basic_TSV(isi);
  initialization_filename();
  initialization_OpenLoops_input(isi);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




void observable_set::Njettiness_calculate_NJ_axes(int i_a){
  static Logger logger("observable_set::Njettiness_calculate_NJ_axes");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // NJ_q_axes are always defined in the hadronic frame with three options:
  // switch_NJcut_axes == 0:   axes defined via jet algorithm (according to input)
  // switch_NJcut_axes == 1:   axes defined via partitioning using 1-jettiness
  
  NJ_q_axes[1] = fourvector(x_pdf[1] * E_beam, 0., 0., x_pdf[1] * E_beam);
  NJ_q_axes[2] = fourvector(x_pdf[2] * E_beam, 0., 0., -x_pdf[2] * E_beam);

  NJ_n_axes[1] = fourvector(1., 0., 0., 1.);
  NJ_n_axes[2] = fourvector(1., 0., 0., -1.);

  NJ_Ei[1] = x_pdf[1] * E_beam;
  NJ_Ei[2] = x_pdf[2] * E_beam;


  if (switch_NJcut_axes == 0){
    for (int i_j = 0; i_j < csi->n_jet_born; i_j++){
      // use first n_jet_born results from the applied jet algorithm here for now...
      // particle_event -> hadronic centre-of-mass system
      fourvector temp_jet = particle_event[access_object["jet"]][i_a][i_j].momentum;
      NJ_n_axes[3 + i_j] = (1. / temp_jet.r()) * fourvector(temp_jet.r(), temp_jet.x1(), temp_jet.x2(), temp_jet.x3());
      // use different normalizations of (light-like) momentum (Eq. (2.6)):
      NJ_Ei[3 + i_j] = Njettiness_calculate_NJ_axes_assigned_energy(temp_jet);
      NJ_q_axes[3 + i_j] = NJ_Ei[3 + i_j] * NJ_n_axes[3 + i_j];
    }
  }
  else if (switch_NJcut_axes == 1){
   // 
    int temp_order = ps_runtime_jet_algorithm[i_a].size() - csi->n_jet_born;
    logger << LOG_DEBUG << "temp_order = " << temp_order << endl;
    if (temp_order == 1){// NLO
      if (csi->n_jet_born == 1){// only for 1-jettiness
	vector<double> tau_1_NLO(3, 0.);
	vector<fourvector> parton(ps_runtime_jet_algorithm[i_a].size());
	for (int i_p = 0; i_p < ps_runtime_jet_algorithm[i_a].size(); i_p++){
	  parton[i_p] = particle_event[0][i_a][ps_runtime_jet_algorithm[i_a][i_p]].momentum;
	}
	tau_1_NLO[0] = parton[0].r() - abs(parton[0].x3());
	tau_1_NLO[1] = parton[1].r() - abs(parton[1].x3());
	tau_1_NLO[2] = parton[0].r() + parton[1].r() - (parton[0] + parton[1]).r();
	for (int i = 0; i < 3; i++){
	  logger << LOG_DEBUG << "tau_1_NLO[" << i << "] = " << tau_1_NLO[i] << endl;
	}

	int min_index = -1;
	double min = 1.e99;
	for (int i_x = 0; i_x < tau_1_NLO.size(); i_x++){
	  if (tau_1_NLO[i_x] < min){min = tau_1_NLO[i_x]; min_index = i_x;}
	}

	for (int i_j = 0; i_j < csi->n_jet_born; i_j++){
	  fourvector temp_jet;
	  if (min_index == 0 || min_index == 1){temp_jet = parton[min_index];}
	  else {temp_jet = parton[0] + parton[1];}
	  NJ_n_axes[3 + i_j] = (1. / temp_jet.r()) * fourvector(temp_jet.r(), temp_jet.x1(), temp_jet.x2(), temp_jet.x3());
	  // use different normalizations of (light-like) momentum (Eq. (2.6)):
	  NJ_Ei[3 + i_j] = Njettiness_calculate_NJ_axes_assigned_energy(temp_jet);
	  NJ_q_axes[3 + i_j] = NJ_Ei[3 + i_j] * NJ_n_axes[3 + i_j];
	  /*
	  // use different normalizations of (light-like) momentum (Eq. (2.6)):
	  double normalization_E = Njettiness_calculate_NJ_axes_assigned_energy(temp_jet);
	  */
	  //	NJ_q_axes[3] = normalization_E * fourvector(temp_jet.r(), temp_jet.x1(), temp_jet.x2(), temp_jet.x3());
	}
      }
      else {
	// generic case could be implemented...
      }
    }
    else if (temp_order == 2){// NNLO

    }   
  }
  else {
    logger << LOG_FATAL << "switch_NJcut_axes not defined!" << endl; exit(1);
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



double observable_set::Njettiness_calculate_NJ_axes_assigned_energy(fourvector & temp_jet){
  static Logger logger("observable_set::Njettiness_calculate_NJ_axes_assigned_energy");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // NJ_Ei are always defined in the hadronic frame with three options:
  // switch_NJcut_axes_energy == 0:   E_i = p_i^0
  // switch_NJcut_axes_energy == 1:   E_i = |p_i|
  // switch_NJcut_axes_energy == 2:   E_i = (p_i^0 + |p_i|) / 2

  if (switch_NJcut_axes_energy == 0){return temp_jet.x0();}
  else if (switch_NJcut_axes_energy == 1){return temp_jet.r();}
  else if (switch_NJcut_axes_energy == 2){return .5 * (temp_jet.x0() + temp_jet.r());}
  else {logger << LOG_FATAL << "switch_NJcut_axes_energy is not properly set!" << endl; exit(1);}
}




void observable_set::Njettiness_calculate_NJ_axes_frame(int i_a){
  Logger logger("observable_set::Njettiness_calculate_NJ_axes_frame");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // NJ_measures are typically defined in special reference frames (exception: 0):
  // switch_NJcut_measure == 0:   invariant-mass measure - Qi = Q
  // switch_NJcut_measure == 1:   geometric measure: hadronic frame (no boost)
  // switch_NJcut_measure == 2:   geometric measure: Born frame (Y)
  // switch_NJcut_measure == 3:   geometric measure: 1-jettiness-axis frame (Y_1)
  // switch_NJcut_measure == 4:   geometric measure: leptonic frame (Y_L)

  NJ_p_parton_frame.resize(ps_runtime_jet_algorithm[i_a].size());
  
  if (switch_NJcut_measure == 0){
    // hadronic frame
    double temp_Q = sqrt(x_pdf[1] * x_pdf[2] * pow(E_CMS, 2));
    for (int i_j = 1; i_j < NJ_Qi.size(); i_j++){
      NJ_q_axes_frame[i_j] = NJ_q_axes[i_j];
      NJ_n_axes_frame[i_j] = NJ_n_axes[i_j];
      NJ_Qi[i_j] = temp_Q;
    }
    for (int i_p = 0; i_p < ps_runtime_jet_algorithm[i_a].size(); i_p++){
      NJ_p_parton_frame[i_p] = particle_event[0][i_a][ps_runtime_jet_algorithm[i_a][i_p]].momentum;
    }
  }
  else if (switch_NJcut_measure == 1){
    // hadronic frame
    for (int i_j = 1; i_j < NJ_Qi.size(); i_j++){
      NJ_q_axes_frame[i_j] = NJ_q_axes[i_j];
      NJ_n_axes_frame[i_j] = NJ_n_axes[i_j];
      NJ_Qi[i_j] = 2 * NJ_Ei[i_j];
    }
    for (int i_p = 0; i_p < ps_runtime_jet_algorithm[i_a].size(); i_p++){
      NJ_p_parton_frame[i_p] = particle_event[0][i_a][ps_runtime_jet_algorithm[i_a][i_p]].momentum;
    }
  }
  else if (switch_NJcut_measure == 2){
    // Born frame (Y)
    double boost_Born = (x_pdf[1] - x_pdf[2]) / (x_pdf[1] + x_pdf[2]);
    if (boost != boost_Born){logger << LOG_DEBUG << "boost = " << boost << " != " << boost_Born << " = boost_Born" << endl;}

    logger.newLine(LOG_DEBUG);
    for (int i_p = 1; i_p < p_parton[i_a].size(); i_p++){
      logger << LOG_DEBUG << "p(Y frame)[" << i_p << "] = " << particle_event[0][i_a][i_p].momentum.zboost(boost_Born) << endl;
    }

    for (int i_j = 1; i_j < NJ_Qi.size(); i_j++){
      NJ_q_axes_frame[i_j] = NJ_q_axes[i_j].zboost(boost_Born);
      NJ_n_axes_frame[i_j] = NJ_n_axes[i_j].zboost(boost_Born);
      NJ_Qi[i_j] = 2 * NJ_q_axes_frame[i_j].x0();
    }
    for (int i_p = 0; i_p < ps_runtime_jet_algorithm[i_a].size(); i_p++){
      NJ_p_parton_frame[i_p] = particle_event[0][i_a][ps_runtime_jet_algorithm[i_a][i_p]].momentum.zboost(boost_Born);
    }
    /*
    p_parton[i_a][ps_runtime_jet_algorithm[i_a][i_p]] == NJ_p_parton_frame[i_p]
    double temp_Y = .5 * log(x_pdf[1] / x_pdf[2]);
    2 * NJ_q_axes_frame[1].x0() == 2 * NJ_q_axes[1].x0() * exp(-temp_Y)
    2 * NJ_q_axes_frame[2].x0() == 2 * NJ_q_axes[2].x0() * exp(+temp_Y)
    double temp_Q2 = x_pdf[1] * x_pdf[2] * pow(E_CMS, 2);
    double temp_Q = sqrt(temp_Q2);
    */
  }
  else if (switch_NJcut_measure == 3){
    // 1-jettiness-axis frame (Y_1)
    // only defined here for 1-jettiness case
    fourvector p_1 = NJ_q_axes[3];
    double boost_1 = p_1.x3() / p_1.x0();

    logger.newLine(LOG_DEBUG);
    for (int i_p = 1; i_p < p_parton[i_a].size(); i_p++){
      logger << LOG_DEBUG << "p(Y_1 frame)[" << i_p << "] = " << particle_event[0][i_a][i_p].momentum.zboost(boost_1) << endl;
    }

    for (int i_j = 1; i_j < NJ_Qi.size(); i_j++){
      NJ_q_axes_frame[i_j] = NJ_q_axes[i_j].zboost(boost_1);
      NJ_n_axes_frame[i_j] = NJ_n_axes[i_j].zboost(boost_1);
      NJ_Qi[i_j] = 2 * NJ_q_axes_frame[i_j].x0();
    }
    for (int i_p = 0; i_p < ps_runtime_jet_algorithm[i_a].size(); i_p++){
      NJ_p_parton_frame[i_p] = particle_event[0][i_a][ps_runtime_jet_algorithm[i_a][i_p]].momentum.zboost(boost_1);
    }
  }
  else if (switch_NJcut_measure == 4){
    // leptonic frame (Y_L)
    fourvector p_L(0., 0., 0., 0.);
    for (int i_p = 3; i_p < 3 + csi->n_particle_born - csi->n_jet_born; i_p++){
      p_L = p_L + particle_event[0][i_a][i_p].momentum;
    }
    double boost_L = p_L.x3() / p_L.x0();

    logger.newLine(LOG_DEBUG);
    for (int i_p = 1; i_p < p_parton[i_a].size(); i_p++){
      logger << LOG_DEBUG << "p(Y_L frame)[" << i_p << "] = " << particle_event[0][i_a][i_p].momentum.zboost(boost_L) << endl;
    }
    
    for (int i_j = 1; i_j < NJ_Qi.size(); i_j++){
      NJ_q_axes_frame[i_j] = NJ_q_axes[i_j].zboost(boost_L);
      NJ_n_axes_frame[i_j] = NJ_n_axes[i_j].zboost(boost_L);
      NJ_Qi[i_j] = 2 * NJ_q_axes_frame[i_j].x0();
    }
    for (int i_p = 0; i_p < ps_runtime_jet_algorithm[i_a].size(); i_p++){
      NJ_p_parton_frame[i_p] = particle_event[0][i_a][ps_runtime_jet_algorithm[i_a][i_p]].momentum.zboost(boost_L);
    }
  }

  logger.newLine(LOG_DEBUG);
  for (int i_q = 1; i_q < NJ_q_axes.size(); i_q++){
    logger << LOG_DEBUG << "NJ_q_axes_frame[" << i_q << "] = " << NJ_q_axes_frame[i_q] << endl;
  }
  logger.newLine(LOG_DEBUG);
  for (int i_j = 1; i_j < NJ_Qi.size(); i_j++){
    logger << LOG_DEBUG << "NJ_Qi[" << i_j << "] = " << NJ_Qi[i_j] << endl;
  }
  logger.newLine(LOG_DEBUG);
  for (int i_q = 1; i_q < NJ_q_axes.size(); i_q++){
    logger << LOG_DEBUG << "NJ_n_axes_frame[" << i_q << "] = " << NJ_n_axes_frame[i_q] << endl;
  }
  logger.newLine(LOG_DEBUG);
  for (int i_p = 0; i_p < ps_runtime_jet_algorithm[i_a].size(); i_p++){
    logger << LOG_DEBUG << "NJ_p_parton_frame[" << i_p << "] = " << NJ_p_parton_frame[i_p] << endl;
  }
 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}






void observable_set::calculate_intermediate_result(phasespace_set & psi){
  Logger logger("observable_set::calculate_intermediate_result");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  calculate_Xsection(psi_i_gen, Xsection, Xsection_delta, full_sum_weight, full_sum_weight2, step_sum_weight, step_sum_weight2);

  logger << LOG_DEBUG << "full_sum_weight = " << full_sum_weight << endl;

  if (switch_CV){
    for (int i_c = 0; i_c < n_qTcut; i_c++){
      for (int i_s = 0; i_s < n_scales_CV; i_s++){
	calculate_Xsection(psi_i_gen, Xsection_CV[i_c][i_s], Xsection_delta_CV[i_c][i_s], full_sum_weight_CV[i_c][i_s], full_sum_weight2_CV[i_c][i_s], step_sum_weight_CV[i_c][i_s], step_sum_weight2_CV[i_c][i_s]);
      }
    }
  }

  if (switch_TSV){
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	  for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	    if (active_qTcut){
	      logger << LOG_DEBUG_VERBOSE << "    sum_weight_qTcut_TSV[" << i_s << "][" << i_q << "][" << i_r << "][" << i_f << "] = " << sum_weight_qTcut_TSV[i_s][i_q][i_r][i_f] << endl;
	      calculate_Xsection(psi_i_gen, Xsection_TSV[i_s][i_q][i_r][i_f], Xsection_delta_TSV[i_s][i_q][i_r][i_f], fullsum_weight_qTcut_TSV[i_s][i_q][i_r][i_f], fullsum_weight2_qTcut_TSV[i_s][i_q][i_r][i_f], sum_weight_qTcut_TSV[i_s][i_q][i_r][i_f], sum_weight2_qTcut_TSV[i_s][i_q][i_r][i_f]);
	      logger << LOG_DEBUG_VERBOSE << "fullsum_weight_qTcut_TSV[" << i_s << "][" << i_q << "][" << i_r << "][" << i_f << "] = " << fullsum_weight_qTcut_TSV[i_s][i_q][i_r][i_f] << endl;
	      logger << LOG_DEBUG_VERBOSE << "      Xsection_TSV[" << i_s << "][" << i_q << "][" << i_r << "][" << i_f << "] = " << Xsection_TSV[i_s][i_q][i_r][i_f] << endl;
	    }
	    else {
	      logger << LOG_DEBUG_VERBOSE << "    sum_weight_TSV[" << i_s << "][" << i_r << "][" << i_f << "] = " << sum_weight_TSV[i_s][i_r][i_f] << endl;
	      calculate_Xsection(psi_i_gen, Xsection_TSV[i_s][i_q][i_r][i_f], Xsection_delta_TSV[i_s][i_q][i_r][i_f], fullsum_weight_TSV[i_s][i_r][i_f], fullsum_weight2_TSV[i_s][i_r][i_f], sum_weight_TSV[i_s][i_r][i_f], sum_weight2_TSV[i_s][i_r][i_f]);
	      logger << LOG_DEBUG_VERBOSE << "fullsum_weight_TSV[" << i_s << "][" << i_r << "][" << i_f << "] = " << fullsum_weight_TSV[i_s][i_r][i_f] << endl;
	      logger << LOG_DEBUG_VERBOSE << "      Xsection_TSV[" << i_s << "][" << i_q << "][" << i_r << "][" << i_f << "] = " << Xsection_TSV[i_s][i_q][i_r][i_f] << endl;
	    }
	  }
	  /*
	  // Should be identical !!!
	  if (active_qTcut){
	    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
	      calculate_Xsection(psi_i_gen, Xsection_TSV[i_s][i_q][i_r][i_f], Xsection_delta_TSV[i_s][i_q][i_r][i_f], fullsum_weight_qTcut_TSV[i_s][i_q][i_r][i_f], fullsum_weight2_qTcut_TSV[i_s][i_q][i_r][i_f], sum_weight_qTcut_TSV[i_s][i_q][i_r][i_f], sum_weight2_qTcut_TSV[i_s][i_q][i_r][i_f]);
	    }
	  }
	  else {
	    calculate_Xsection(psi_i_gen, Xsection_TSV[i_s][0][i_r][i_f], Xsection_delta_TSV[i_s][0][i_r][i_f], fullsum_weight_TSV[i_s][i_r][i_f], fullsum_weight2_TSV[i_s][i_r][i_f], sum_weight_TSV[i_s][i_r][i_f], sum_weight2_TSV[i_s][i_r][i_f]);
	  }
	  */
	}
      }
    }
  }

  if ((psi_i_acc == psi_n_events_max) || (psi_i_acc >= psi_n_events_min && abs(Xsection_delta / sigma_normalization) < sigma_normalization_deviation)){int_end = 1;}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::calculate_Xsection(long long i_gen, double & this_Xsection, double & this_Xsection_delta, double & this_sum_weights, double & this_sum_weights2, double & this_temp_sum_weights, double & this_temp_sum_weights2){
  Logger logger("observable_set::calculate_Xsection");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  this_sum_weights += this_temp_sum_weights;
  this_sum_weights2 += this_temp_sum_weights2;
  this_temp_sum_weights = 0.;
  this_temp_sum_weights2 = 0.;
  this_Xsection_delta = sqrt((this_sum_weights2 - pow(this_sum_weights, 2) / i_gen)) / (i_gen - 1);
  this_Xsection = this_sum_weights / i_gen;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_p_parton(phasespace_set & psi){
  Logger logger("observable_set::determine_p_parton");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  p_parton = start_p_parton;
  for (int xi = 0; xi <= csi->n_particle + 2; xi++){p_parton[0][psi_o_map[0][xi]] = psi_xbp_all[0][intpow(2, xi - 1)];}
  x_pdf = psi_x_pdf;
  boost = psi_boost;
  /*
  logger << LOG_DEBUG_VERBOSE << "type_contribution = " << type_contribution << endl;
  if (type_contribution == "CA" || type_contribution == "RCA"){z_coll = psi_z_coll;}
  */
  for (int xi = 0; xi <= csi->n_particle + 2; xi++){logger << LOG_DEBUG_VERBOSE << "p_parton[0][" << psi_o_map[0][xi] << "] = " << p_parton[0][psi_o_map[0][xi]] << endl;}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::initialization_runtime(){
  Logger logger("observable_set::initialization_runtime");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  start = clock();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::determine_runtime(phasespace_set & psi){
  Logger logger("observable_set::determine_runtime");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (start > clock()){start = clock(); time_counter = 0;}
  if (start < 0 && clock() > 0){start = clock(); time_counter = 0;}
  
  //  int seconds_passed = clock()/CLOCKS_PER_SEC - start/CLOCKS_PER_SEC;
  int seconds_passed = clock()/CLOCKS_PER_SEC - start/CLOCKS_PER_SEC + sec_import;
  
  int old_h = h;
  int old_min = min;
  
  if (seconds_passed > 3600) {
    h = int(seconds_passed / 3600);
    min = 0;
    seconds_passed -= h * 3600;
  }
  if (seconds_passed > 60) {
    min = int(seconds_passed / 60);
    seconds_passed -= min * 60;
  }
  sec = seconds_passed;
  //  sec = seconds_passed + sec_import;
  
  /*
  static long long i_acc_old = 0;
  static long long i_gen_old = 0;
  double efficiency = 0;

  if (min > old_min || h > old_h){
    //  if ((min > old_min && h == 0) || h > old_h){
    if (i_gen_old != psi_i_gen) {
      efficiency = double(psi_i_acc-i_acc_old)/(psi_i_gen-i_gen_old);

      i_acc_old = psi_i_acc;
      i_gen_old = psi_i_gen;
    }
  }
  */

  if ((switch_console_output_runtime == 1 && ((min > old_min && h == 0) || h > old_h)) ||
      (switch_console_output_runtime == 2 && (min > old_min || h > old_h))){
    // complete efficiency, not only that one of the last interval.
    double efficiency = double(psi_i_acc) / psi_i_gen;

    logger << LOG_INFO << setw(3) << right << h << " hours" << setw(3) << right << min << " minutes needed for " << psi_i_acc << " accepted (" << psi_i_gen << " generated) events (efficiency = " << efficiency << ")." << endl;
  }
  /*
  if (h > old_h) {
    logger << LOG_INFO << setw(3) << right << h << " hours" << setw(3) << right << min << " minutes needed for " << psi_i_acc << " accepted (" << psi_i_gen << " generated) events (efficiency = " << efficiency << ")." << endl;
  }
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::determine_techcut_RA(phasespace_set & psi){
  Logger logger("observable_set::determine_techcut_RA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  psi_RA_techcut = 0;
  static int num_printouts = 0;
  const int max_printouts = 1000;
  static double severity_threshold = 1.e-12;
  static double impact_threshold = 1.e-3;

  // check singular regions (p_i * p_j / s_hat) < cut_technical
  for (int sr = 0; sr < psi_RA_singular_region_list.size(); sr++){
    int x1 = psi_RA_singular_region_list[sr][0];
    int x2 = psi_RA_singular_region_list[sr][1];
    psi_RA_singular_region[x1][x2] = p_parton[0][x1] * p_parton[0][x2] / psi_xbs_all[0][0];

    if (psi_RA_singular_region[x1][x2] < psi_cut_technical){
      int ccount = 0;
      for (int i_a = 1; i_a < cut_ps.size(); i_a++){if (cut_ps[i_a] == -1){ccount++;}}

      double severity = abs(1. + accumulate(RA_ME2.begin() + 1, RA_ME2.end(), 0.) / RA_ME2[0]);
      //      double severity = 1.;

      // suppress warning if point is (relatively) harmless and output_level is low
      //      if (csi->type_contribution != "L2RA" &&
      //	  (((severity > severity_threshold || Log::getLogThreshold() <= LOG_DEBUG) && 
      if ((((severity > severity_threshold || Log::getLogThreshold() <= LOG_DEBUG) && 
	    (abs(integrand / sigma_normalization) > impact_threshold || sigma_normalization == 1.) && 
	    num_printouts < max_printouts) || 
	   Log::getLogThreshold() <= LOG_DEBUG_VERBOSE)){
        logger << LOG_DEBUG << endl;
        logger << LOG_DEBUG << "applying technical cut near dipole singularity" << endl;
        logger << LOG_DEBUG << "cut_technical = " << psi_cut_technical << endl;
        for (int ib = 0; ib < p_parton[0].size(); ib++){logger << LOG_DEBUG << "p_parton[0][" << ib << "] = " << p_parton[0][ib] << "   " << p_parton[0][ib].m2() << "   " << sqrt(abs(p_parton[0][ib].m2())) << endl;}
        logger << LOG_DEBUG << "contributing dipoles:   " << endl;
        for (int i_a = 0; i_a < cut_ps.size(); i_a++){
          if (cut_ps[i_a] != -1){logger << LOG_DEBUG << setw(8) << (*RA_dipole)[i_a].name() << "   " << "RA_ME2[" << i_a << "] = " << RA_ME2[i_a] << endl;}
        }

	if (switch_console_output_techcut_RA){
	  num_printouts++;
	  logger << LOG_INFO 
		 << right << setw(12) << psi_i_gen 
		 << right << setw(10) << psi_i_acc 
		 << "   (" 
		 << right << setw(2) << cut_ps[0] << "/" 
		 << right << setw(2) << ccount 
		 << "   " 
		 << psi_RA_singular_region_name[x1][x2] << "/^s = " 
		 << right << setw(15) << setprecision(8) << showpoint << psi_RA_singular_region[x1][x2] 
		 << "):   int/LO = " 
		 << right << setw(15) << setprecision(8) << showpoint << integrand / sigma_normalization 
		 << "   A/R = " 
		 << left << setw(15) << setprecision(8) << accumulate(RA_ME2.begin() + 1, RA_ME2.end(), 0.) / RA_ME2[0] 
		 << endl;
	}
      }

      /*
      // RA_techcut_integrand could be used to check the dependence on the parameter cut_technical !!!
      if (psi_RA_singular_region[x1][x2] > psi_cut_technical / 10.){RA_techcut_integrand = integrand;}
      else {RA_techcut_integrand = 0.;}
      */
      /*
      if (abs(integrand / sigma_normalization) > 1. || psi_RA_singular_region[x1][x2] < psi_cut_technical / 100.) {
	psi_RA_techcut = 1; 
	psi_i_tec++;
	integrand = 0.;
      }
      */
      // new:
      // All points with (p_i * p_j / s_hat) < cut_technical are technically cut.
      // This might require an adaptation of the corresponding values...
      psi_RA_techcut = 1; 
      psi_i_tec++;
      integrand = 0.;
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



///////////////////////
//  access elements  //
///////////////////////

ostream & operator << (ostream & s, const observable_set & osi){
  s << "observable_set independent stuff:" << endl;
  s << setw(38) << right << "switch_qTcut" << setw(5) << "" << "=" << setw(5) << "" << osi.switch_qTcut << endl;
  s << setw(38) << right << "n_qTcut" << setw(5) << "" << "=" << setw(5) << "" << osi.n_qTcut << endl;
  s << setw(38) << right << "min_qTcut" << setw(5) << "" << "=" << setw(5) << "" << osi.min_qTcut << endl;
  s << setw(38) << right << "step_qTcut" << setw(5) << "" << "=" << setw(5) << "" << osi.step_qTcut << endl;
  s << setw(38) << right << "switch_distribution_at_all_TSV" << setw(5) << "" << "=" << setw(5) << "" << osi.switch_distribution_at_all_TSV << endl;

    s << "observable_set output:   n_set_TSV = " << osi.n_set_TSV << endl << endl;

  if (osi.switch_TSV == 0){
    s << "observable_set output:   TSV has been switched off." << endl << endl;
  }
  else {
    s << "observable_set output:   n_set_TSV = " << osi.n_set_TSV << endl << endl;
    for (int i_s = 0; i_s < osi.n_set_TSV; i_s++){
      s << setw(35) << right << "name_set_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.name_set_TSV[i_s] << endl;
      s << setw(35) << right << "dynamic_scale_ren_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.dynamic_scale_ren_TSV[i_s] << endl;
      s << setw(35) << right << "central_scale_ren_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.central_scale_ren_TSV[i_s] << endl;
      s << setw(35) << right << "relative_central_scale_ren_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.relative_central_scale_ren_TSV[i_s] << endl;
      s << setw(35) << right << "n_scale_ren_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.n_scale_ren_TSV[i_s] << endl;
      s << setw(35) << right << "factor_scale_ren_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.factor_scale_ren_TSV[i_s] << endl;
      s << setw(35) << right << "dynamic_scale_fact_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.dynamic_scale_fact_TSV[i_s] << endl;
      s << setw(35) << right << "central_scale_fact_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.central_scale_fact_TSV[i_s] << endl;
      s << setw(35) << right << "relative_central_scale_fact_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.relative_central_scale_fact_TSV[i_s] << endl;
      s << setw(35) << right << "n_scale_fact_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.n_scale_fact_TSV[i_s] << endl;
      s << setw(35) << right << "factor_scale_fact_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.factor_scale_fact_TSV[i_s] << endl;
      s << setw(35) << right << "min_qTcut_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.min_qTcut_TSV[i_s] << endl;
      s << setw(35) << right << "max_qTcut_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.max_qTcut_TSV[i_s] << endl;
      s << setw(35) << right << "switch_distribution_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.switch_distribution_TSV[i_s] << endl;
      s << setw(35) << right << "max_n_integrand_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.max_n_integrand_TSV[i_s] << endl;
      s << setw(35) << right << "min_qTcut_distribution_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.min_qTcut_distribution_TSV[i_s] << endl;
      s << setw(35) << right << "max_qTcut_distribution_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.max_qTcut_distribution_TSV[i_s] << endl;
      s << setw(35) << right << "switch_moment_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.switch_moment_TSV[i_s] << endl;
      s << setw(35) << right << "filename_integration_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.filename_integration_TSV[i_s] << endl;
      s << setw(35) << right << "filename_result_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.filename_result_TSV[i_s] << endl;
      s << setw(35) << right << "filename_distribution_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.filename_distribution_TSV[i_s] << endl;
      s << setw(35) << right << "filename_moment_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.filename_moment_TSV[i_s] << endl;
      s << endl;
    }
    for (int i_s = 0; i_s < osi.n_extended_set_TSV; i_s++){
      s << setw(35) << right << "name_extended_set_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.name_extended_set_TSV[i_s] << endl;
      s << setw(35) << right << "n_scale_ren_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.n_scale_ren_TSV[i_s] << endl;
      s << setw(35) << right << "n_scale_fact_TSV[" << setw(2) << i_s << "]" << setw(5) << "" << "=" << setw(5) << "" << osi.n_scale_fact_TSV[i_s] << endl;
    }

    s << "observable_set output:   n_set_TSV = " << osi.n_set_TSV << endl << endl;
  }
  return s;
}


