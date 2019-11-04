#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

void observable_set::initialization_OpenLoops_input(inputparameter_set & isi){
  Logger logger("observable_set::initialization_OpenLoops_input (isi)");
  logger << LOG_DEBUG << "called" << endl;

  OL_parameter = isi.OL_parameter;
  OL_value = isi.OL_value;

  logger << LOG_DEBUG << "finished" << endl;
}



void observable_set::initialization_OpenLoops_parameter(phasespace_set & psi){
  static Logger logger("observable_set::initialization_OpenLoops_parameter");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ////////////////////////////////////
  //  N_nondecoupled determination  //
  ////////////////////////////////////

  //  string LHAPDFversion = LHAPDF::version();
#ifdef LHAPDF5
  string LHAPDFversion = "5";
#else
  string LHAPDFversion = "6";
#endif

  logger << LOG_INFO << "LHAPDFversion = " << LHAPDFversion << endl;
  logger << LOG_INFO << "LHAPDFversion.substr(0, 1) = " << LHAPDFversion.substr(0, 1) << endl;

  int LHAPDF_N_nondecoupled = LHAPDF::getNf();

  if (LHAPDFversion.substr(0, 1) == "6"){
    logger << LOG_INFO << "LHAPDF version 6 (" << LHAPDFversion << ") is used." << endl;
    N_nondecoupled = LHAPDF_N_nondecoupled;
    if (LHAPDFname.substr(0, 10) == "NNPDF21_lo"){N_nondecoupled = 6;}
  }
  else if (LHAPDFversion.substr(0, 1) == "5"){
    logger << LOG_INFO << "LHAPDF version 5 (" << LHAPDFversion << ") is used." << endl;

    logger << LOG_INFO << "LHAPDFname = " << LHAPDFname << endl;
    logger << LOG_INFO << "LHAPDF_N_nondecoupled = " << LHAPDF_N_nondecoupled << endl;
    
    if (LHAPDFname.substr(0, 7) == "NNPDF21"){
      logger << LOG_DEBUG << "LHAPDFname = " << LHAPDFname << endl;
      N_nondecoupled = 6;
      for (int i_s = LHAPDFname.size() - 3; i_s >=0; i_s--){
	string temp_Nf = LHAPDFname.substr(i_s, 3);
	logger << LOG_DEBUG << "LHAPDFname.substr(" << i_s << ", 3) = " << temp_Nf << endl;
	if (temp_Nf == "NF3"){N_nondecoupled = 3;}
	else if (temp_Nf == "NF4"){N_nondecoupled = 4;}
	else if (temp_Nf == "NF5"){N_nondecoupled = 5;}
      }
      //      logger << LOG_DEBUG << "N_nondecoupled = " << N_nondecoupled << endl;
    }

    else if (LHAPDFname.substr(0, 7) == "NNPDF23"){
      logger << LOG_DEBUG << "LHAPDFname = " << LHAPDFname << endl;
      N_nondecoupled = 6;
      for (int i_s = LHAPDFname.size() - 3; i_s >=0; i_s--){
	string temp_Nf = LHAPDFname.substr(i_s, 3);
	logger << LOG_DEBUG << "LHAPDFname.substr(" << i_s << ", 3) = " << temp_Nf << endl;
	if (temp_Nf == "NF3"){N_nondecoupled = 3;}
	else if (temp_Nf == "NF4"){N_nondecoupled = 4;}
	else if (temp_Nf == "NF5"){N_nondecoupled = 5;}
      }
      //      logger << LOG_DEBUG << "N_nondecoupled = " << N_nondecoupled << endl;
    }

    else if (LHAPDFname.substr(0, 6) == "NNPDF3"){
      logger << LOG_DEBUG << "LHAPDFname = " << LHAPDFname << endl;
      N_nondecoupled = 5;
      for (int i_s = LHAPDFname.size() - 4; i_s >=0; i_s--){
	string temp_Nf = LHAPDFname.substr(i_s, 4);
	logger << LOG_DEBUG << "LHAPDFname.substr(" << i_s << ", 4) = " << temp_Nf << endl;
	if (temp_Nf == "nf_3"){N_nondecoupled = 3;}
	else if (temp_Nf == "nf_4"){N_nondecoupled = 4;}
	else if (temp_Nf == "nf_6"){N_nondecoupled = 6;}
      }
      //      logger << LOG_DEBUG << "N_nondecoupled = " << N_nondecoupled << endl;
    }

    else if (LHAPDFname.substr(0, 4) == "CT14"){
      logger << LOG_DEBUG << "LHAPDFname = " << LHAPDFname << endl;
      N_nondecoupled = 5;
      for (int i_s = LHAPDFname.size() - 3; i_s >=0; i_s--){
	string temp_Nf = LHAPDFname.substr(i_s, 3);
	logger << LOG_DEBUG << "LHAPDFname.substr(" << i_s << ", 3) = " << temp_Nf << endl;
	if (temp_Nf == "NF3"){N_nondecoupled = 3;}
	else if (temp_Nf == "NF4"){N_nondecoupled = 4;}
	else if (temp_Nf == "NF6"){N_nondecoupled = 6;}
      }
      //      logger << LOG_DEBUG << "N_nondecoupled = " << N_nondecoupled << endl;
    }

    else if (LHAPDFname.substr(0, 8) == "MMHT2014"){
      logger << LOG_DEBUG << "LHAPDFname = " << LHAPDFname << endl;
      N_nondecoupled = 5;
      for (int i_s = LHAPDFname.size() - 3; i_s >=0; i_s--){
	string temp_Nf = LHAPDFname.substr(i_s, 3);
	logger << LOG_DEBUG << "LHAPDFname.substr(" << i_s << ", 3) = " << temp_Nf << endl;
	if (temp_Nf == "nf3"){N_nondecoupled = 3;}
	else if (temp_Nf == "nf4"){N_nondecoupled = 4;}
	else if (temp_Nf == "nf6"){N_nondecoupled = 6;}
      }
      //      logger << LOG_DEBUG << "N_nondecoupled = " << N_nondecoupled << endl;
    }

    else if (LHAPDFname.substr(0, 9) == "PDF4LHC15"){
      logger << LOG_DEBUG << "LHAPDFname = " << LHAPDFname << endl;
      N_nondecoupled = 5;
      for (int i_s = LHAPDFname.size() - 3; i_s >=0; i_s--){
	string temp_Nf = LHAPDFname.substr(i_s, 3);
	logger << LOG_DEBUG << "LHAPDFname.substr(" << i_s << ", 3) = " << temp_Nf << endl;
	if (temp_Nf == "nf3"){N_nondecoupled = 3;}
	else if (temp_Nf == "nf4"){N_nondecoupled = 4;}
	else if (temp_Nf == "nf6"){N_nondecoupled = 6;}
      }
      //      logger << LOG_DEBUG << "N_nondecoupled = " << N_nondecoupled << endl;
    }

    else if (LHAPDFname.substr(0, 6) == "LUXqed"){
      logger << LOG_DEBUG << "LHAPDFname = " << LHAPDFname << endl;
      N_nondecoupled = 5;
      for (int i_s = LHAPDFname.size() - 3; i_s >=0; i_s--){
	string temp_Nf = LHAPDFname.substr(i_s, 3);
	logger << LOG_DEBUG << "LHAPDFname.substr(" << i_s << ", 3) = " << temp_Nf << endl;
	if (temp_Nf == "nf3"){N_nondecoupled = 3;}
	else if (temp_Nf == "nf4"){N_nondecoupled = 4;}
	else if (temp_Nf == "nf6"){N_nondecoupled = 6;}
      }
      //      logger << LOG_DEBUG << "N_nondecoupled = " << N_nondecoupled << endl;
    }

    else {
      N_nondecoupled = LHAPDF_N_nondecoupled;
    }
  }

  if (LHAPDFname.substr(0, 8) == "MSTW2008" || LHAPDFname.substr(0, 8) == "MMHT2014"){
    for (int i_s = LHAPDFname.size() - 6; i_s >=0; i_s--){
      string temp_Nf = LHAPDFname.substr(i_s, 6);
      logger << LOG_DEBUG << "LHAPDFname.substr(" << i_s << ", 6) = " << temp_Nf << endl;
      if (temp_Nf == "nf4as5"){N_nondecoupled = 5;}
    }
  }
  
  if (N_nondecoupled != LHAPDF_N_nondecoupled){
    logger << LOG_WARN << "N_nondecoupled has been changed wrt. LHAPDF_N_nondecoupled in LHAPDF " << LHAPDFversion << endl;
  }

  logger << LOG_INFO << "N_nondecoupled = " << N_nondecoupled << endl;


  // generic version to replace mass/width initialization by particle name
  //  for (int i_p = 1; i_p < psi_M.size(); i_p++){
  for (int i_p = 1; i_p < 26; i_p++){
    if ((i_p > 6 && i_p < 11) || (i_p > 16 && i_p < 21) || i_p == 12 || i_p == 14 || i_p == 16 || i_p == 21 || i_p == 22){continue;}
    stringstream temp_mass_ss;
    temp_mass_ss << "mass(" << i_p << ")";
    //    string temp_mass_s = temp_mass_ss.str();
    ol_setparameter_double(stch(temp_mass_ss.str()), psi_M[i_p]);
  }

  //  ol_setparameter_double(stch("yuk(5)"), psi_M[5]);



  //  for (int i_p = 1; i_p < psi_Gamma.size(); i_p++){
  for (int i_p = 1; i_p < 26; i_p++){
    if ((i_p > 6 && i_p < 11) || (i_p > 16 && i_p < 21) || i_p == 12 || i_p == 14 || i_p == 16 || i_p == 21 || i_p == 22){continue;}
    stringstream temp_width_ss;
    temp_width_ss << "width(" << i_p << ")";
    //    string temp_width_s = temp_width_ss.str();
    ol_setparameter_double(stch(temp_width_ss.str()), psi_Gamma[i_p]);
  }

  ol_setparameter_int(stch("n_quarks"), N_quarks);
  /*
  int minnf_alphasrun = p_pdf->info().get_entry_as<int>("NumFlavors");
  logger << LOG_INFO << "manual:    N_nondecoupled = " << N_nondecoupled << endl;
  logger << LOG_INFO << "automatic: minnf_alphasrun = " << minnf_alphasrun << endl;
  assert(N_nondecoupled == minnf_alphasrun);
  */
  logger << LOG_INFO << "N_nondecoupled = " << N_nondecoupled << endl;
  ol_setparameter_int(stch("nq_nondecoupled"), N_nondecoupled); // should be calculated from chosen PDF set
  
  logger << LOG_INFO << "scale_ren = " << scale_ren << endl;
  ol_setparameter_double(stch("renscale"), scale_ren);

  ol_setparameter_int(stch("polenorm"), switch_polenorm);
  logger << LOG_INFO << "switch_polenorm = " << switch_polenorm << endl;
  //  ol_setparameter_int(stch("flavour_mapping"), 0);
  //  ol_setparameter_int(stch("ew_renorm_scheme"), 1);
  ol_setparameter_int(stch("verbose"), 2);

  logger << LOG_INFO << "msi.ew_scheme = " << msi.ew_scheme << endl;
  ol_setparameter_int(stch("ew_scheme"), msi.ew_scheme);

  /*
  logger << LOG_INFO << "msi.ew_renorm_scheme = " << msi.ew_scheme << " (per default equal to ew_scheme)" << endl;
  ol_setparameter_int(stch("ew_renorm_scheme"), msi.ew_scheme);
  */
  
  logger << LOG_INFO << "msi.use_cms = " << msi.use_cms << endl;
  ol_setparameter_int(stch("use_cms"), msi.use_cms);

  logger << LOG_INFO << "Gmu = " << msi.G_F << endl;
  ol_setparameter_double(stch("Gmu"), msi.G_F);

  logger << LOG_INFO << "alpha_S = " << alpha_S << endl;
  ol_setparameter_double(stch("alpha_s"), alpha_S);


  // CKM matrix
  if (msi.CKM_matrix != "trivial"){
    ol_setparameter_int(stch("ckmorder"), 1);
    ol_setparameter_double(stch("VCKMdu"), msi.V_du);
    ol_setparameter_double(stch("VCKMsu"), msi.V_su);
    ol_setparameter_double(stch("VCKMbu"), msi.V_bu);
    ol_setparameter_double(stch("VCKMdc"), msi.V_dc);
    ol_setparameter_double(stch("VCKMsc"), msi.V_sc);
    ol_setparameter_double(stch("VCKMbc"), msi.V_bc);
    ol_setparameter_double(stch("VCKMdt"), msi.V_dt);
    ol_setparameter_double(stch("VCKMst"), msi.V_st);
    ol_setparameter_double(stch("VCKMbt"), msi.V_bt);
  }
  else {
    ol_setparameter_int(stch("ckmorder"), 0);
  }


  int last_switch = 1;
  ol_setparameter_int(stch("last_switch"), last_switch);

  ///  int amp_switch = 1;
  ///  ol_setparameter_int(stch("redlib1"), amp_switch);
  ///  //  cout << "redlib1 is set to " << amp_switch << "." << endl;

  ///  int amp_switch_rescue = 7;
  ///  ol_setparameter_int(stch("redlib2"), amp_switch_rescue);
  ///  //  cout << "redlib2 is set to " << amp_switch_rescue << "." << endl;
  
  int use_coli_cache = 1;
  ol_setparameter_int(stch("use_coli_cache"), use_coli_cache);

  int out_symmetry = 1;
  ol_setparameter_int(stch("out_symmetry"), out_symmetry);

  int leading_colour = 0;
  ol_setparameter_int(stch("leading_colour"), leading_colour);

  int n_log = 1;
  //  int n_log = n_step;
  ol_setparameter_int(stch("stability_log"), n_log);
  
  ///  int xtest;
  ///  ol_getparameter_int(stch("stability_mode"), &xtest);
  ///  logger << LOG_DEBUG << "default: stability_mode = " << xtest << endl;

  ///  int stability_mode = 23;
  ///  ol_setparameter_int(stch("stability_mode"), stability_mode);


  // ???
  double pole1_UV = 0.;
  ol_setparameter_double(stch("pole_uv"), pole1_UV);
  
  double pole1_IR = 0.;
  ol_setparameter_double(stch("pole_ir1"), pole1_IR);

  double pole2_IR = 0.;
  ol_setparameter_double(stch("pole_ir2"), pole2_IR);
  

  int CT_on = 0; // modification of ...ME2_VA and OpenLoops calls needed !!!
  ol_setparameter_int(stch("ct_on"), CT_on); // modification of ...ME2_VA and OpenLoops calls needed !!!

  int R2_on = 1;
  ol_setparameter_int(stch("r2_on"), R2_on);

  int IR_on = 1;
  ol_setparameter_int(stch("ir_on"), IR_on);

  // ???
  double fact_uv = 1.;
  ol_setparameter_double(stch("fact_uv"), fact_uv);

  double fact_ir = 1.;
  ol_setparameter_double(stch("fact_ir"), fact_ir);
  
  
  int polecheck = 0;
  ol_setparameter_int(stch("polecheck"), polecheck);
  
  
  for (int i_ol = 0; i_ol < OL_parameter.size(); i_ol++){
    logger << LOG_DEBUG << "OL_parameter[" << i_ol << "] = " << OL_parameter[i_ol] << " = " << OL_value[i_ol] << endl;
    
    string temp_parameter = OL_parameter[i_ol];
    string temp_value = OL_value[i_ol];
    char * char_parameter = new char[temp_parameter.size() + 1];
    std::copy(temp_parameter.begin(), temp_parameter.end(), char_parameter);
    char_parameter[temp_parameter.size()] = '\0';
    char * char_value = new char[temp_value.size() + 1];
    std::copy(temp_value.begin(), temp_value.end(), char_value);
    char_value[temp_value.size()] = '\0';
    ol_setparameter_string(char_parameter, char_value);
    
    //    ol_setparameter_string(stch(OL_parameter[i_ol]), stch(OL_value[i_ol]));
  }

  if (msi.ew_scheme == -1){
    ol_setparameter_double(stch("alpha_QED"), msi.alpha_e_Gmu);
  }
  if (msi.ew_scheme == 0){
    ol_setparameter_double(stch("alpha_QED_0"), msi.alpha_e_0);
  }
  if (msi.ew_scheme == 1){
    ol_setparameter_double(stch("alpha_QED"), msi.alpha_e_Gmu);
  }
  if (msi.ew_scheme == 2){
    ol_setparameter_double(stch("alpha_QED_MZ"), msi.alpha_e_MZ);
  }

  /*
  ol_setparameter_double(stch("alpha_QED_0"), msi.alpha_e_0);
  ol_setparameter_double(stch("alpha_QED"), msi.alpha_e_Gmu);
  ol_setparameter_double(stch("alpha_QED_MZ"), msi.alpha_e_MZ);
  */
  
  ol_start();
 
  double temp_alpha = 0.;
  ol_getparameter_double(stch("alpha_QED"), &temp_alpha);
  logger << LOG_INFO << "OpenLoops:   msi.alpha_e     = " << setw(23) << setprecision(15) << temp_alpha << "   1 / msi.alpha_e     = " << setw(23) << setprecision(15) << 1. / temp_alpha << endl;
  ol_getparameter_double(stch("alpha_QED_0"), &temp_alpha);
  logger << LOG_INFO << "OpenLoops:   msi.alpha_e_0   = " << setw(23) << setprecision(15) << temp_alpha << "   1 / msi.alpha_e_0   = " << setw(23) << setprecision(15) << 1. / temp_alpha << endl;
  ol_getparameter_double(stch("alpha_QED_MZ"), &temp_alpha);
  logger << LOG_INFO << "OpenLoops:   msi.alpha_e_MZ  = " << setw(23) << setprecision(15) << temp_alpha << "   1 / msi.alpha_e_MZ  = " << setw(23) << setprecision(15) << 1. / temp_alpha << endl;


  /*
  logger << LOG_INFO << "MUNICH:  msi.alpha_e     = " << setw(23) << setprecision(15) << msi.alpha_e << "   1 / msi.alpha_e     = " << setw(23) << setprecision(15) << 1. / msi.alpha_e << endl;
  ol_getparameter_double(stch("alpha_QED"), &temp_alpha);
  
  logger << LOG_INFO << "Before:  msi.alpha_e     = " << setw(23) << setprecision(15) << temp_alpha << "   1 / msi.alpha_e     = " << setw(23) << setprecision(15) << 1. / temp_alpha << endl;
  ol_setparameter_double(stch("alpha_QED"), msi.alpha_e);
  ol_getparameter_double(stch("alpha_QED"), &temp_alpha);
  logger << LOG_INFO << "After:   msi.alpha_e     = " << setw(23) << setprecision(15) << temp_alpha << "   1 / msi.alpha_e     = " << setw(23) << setprecision(15) << 1. / temp_alpha << endl;
  
  if (msi.ew_scheme == 1){
    logger << LOG_INFO << "MUNICH:  msi.alpha_e_Gmu = " << setw(23) << setprecision(15) << msi.alpha_e_Gmu << "   1 / msi.alpha_e_Gmu = " << setw(23) << setprecision(15) << 1. / msi.alpha_e_Gmu << endl;
    ol_getparameter_double(stch("alpha_QED"), &temp_alpha);
    logger << LOG_INFO << "Before:  msi.alpha_e_Gmu = " << setw(23) << setprecision(15) << temp_alpha << "   1 / msi.alpha_e_Gmu = " << setw(23) << setprecision(15) << 1. / temp_alpha << endl;
    ol_setparameter_double(stch("alpha_QED"), msi.alpha_e_Gmu);
    ol_getparameter_double(stch("alpha_QED"), &temp_alpha);
    logger << LOG_INFO << "After:   msi.alpha_e_Gmu = " << setw(23) << setprecision(15) << temp_alpha << "   1 / msi.alpha_e_Gmu = " << setw(23) << setprecision(15) << 1. / temp_alpha << endl;
  }
  
  logger << LOG_INFO << "MUNICH:  msi.alpha_e_0   = " << setw(23) << setprecision(15) << msi.alpha_e_0 << "   1 / msi.alpha_e_0   = " << setw(23) << setprecision(15) << 1. / msi.alpha_e_0 << endl;
  ol_getparameter_double(stch("alpha_QED"), &temp_alpha);
  logger << LOG_INFO << "Before:  msi.alpha_e_0   = " << setw(23) << setprecision(15) << temp_alpha << "   1 / msi.alpha_e_0   = " << setw(23) << setprecision(15) << 1. / temp_alpha << endl;
  ol_setparameter_double(stch("alpha_QED_0"), msi.alpha_e_0);
  ol_getparameter_double(stch("alpha_QED_0"), &temp_alpha);
  logger << LOG_INFO << "After:   msi.alpha_e_0   = " << setw(23) << setprecision(15) << temp_alpha << "   1 / msi.alpha_e_0   = " << setw(23) << setprecision(15) << 1. / temp_alpha << endl;
  
  logger << LOG_INFO << "MUNICH:  msi.alpha_e_MZ  = " << setw(23) << setprecision(15) << msi.alpha_e_MZ << "   1 / msi.alpha_e_MZ  = " << setw(23) << setprecision(15) << 1. / msi.alpha_e_MZ << endl;
  ol_getparameter_double(stch("alpha_QED_MZ"), &temp_alpha);
  logger << LOG_INFO << "Before:  msi.alpha_e_MZ  = " << setw(23) << setprecision(15) << temp_alpha << "   1 / msi.alpha_e_MZ  = " << setw(23) << setprecision(15) << 1. / temp_alpha << endl;
  ol_setparameter_double(stch("alpha_QED_MZ"), msi.alpha_e_MZ);
  ol_getparameter_double(stch("alpha_QED_MZ"), &temp_alpha);
  logger << LOG_INFO << "After:   msi.alpha_e_MZ  = " << setw(23) << setprecision(15) << temp_alpha << "   1 / msi.alpha_e_MZ  = " << setw(23) << setprecision(15) << 1. / temp_alpha << endl;
  */

  
  /*
  ol_getparameter_int("stability_mode", &xtest);
  logger << LOG_DEBUG << "input:   stability_mode = " << xtest << endl;
  */
  /*
  ol_getparameter_int(stch("polenorm"), &test_polenorm);
  cout << "test_polenorm = " << test_polenorm << endl;
  */
  //  OLP_PrintParameter("olparameters.txt");
  //  ol_setparameter_int(stch("parameters_verbose"), 1);
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::initialization_OpenLoops_process(phasespace_set & psi){
  static Logger logger("observable_set::initialization_OpenLoops_process");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  initialization_OpenLoops_parameter(psi);

  logger << LOG_DEBUG << "csi->type_perturbative_order         = " << csi->type_perturbative_order << endl;
  logger << LOG_DEBUG << "csi->type_contribution               = " << csi->type_contribution << endl;
  logger << LOG_DEBUG << "csi->type_correction                 = " << csi->type_correction << endl;
  logger << LOG_DEBUG << "csi->contribution_order_alpha_s      = " << csi->contribution_order_alpha_s << endl;
  logger << LOG_DEBUG << "csi->contribution_order_alpha_e      = " << csi->contribution_order_alpha_e << endl;
  logger << LOG_DEBUG << "csi->contribution_order_interference = " << csi->contribution_order_interference << endl;

  process_id = -1;
  //  int process_id = -1;
  if (csi->type_contribution == "born" || 
      csi->type_contribution == "RT" || 
      csi->type_contribution == "RJ"){
    int type_amplitude = 1;
    if (user.string_value[user.string_map["model"]] == "Bornloop"){type_amplitude = 12;}
    //  select_OL_born_mode(QCD_order, QEW_order);
    logger << LOG_DEBUG << "csi->contribution_order_interference  = " << csi->contribution_order_interference << endl;
    if (csi->contribution_order_interference == 0){
      ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
      logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
      process_id = register_OL_subprocess(0, type_amplitude);
      if (process_id == -1){
	ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	process_id = register_OL_subprocess(0, type_amplitude);
      }
    }
    else if (csi->contribution_order_interference == 1){
      /*
	ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	process_id = register_OL_subprocess(0, 11);
	logger << LOG_DEBUG << "process_id (QCD corr. to int)    = " << process_id << endl;
	// usual initialization for QCD processes, where QCD virtual amplitudes are available
	*/
      //      if (process_id == -1){
      ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
      logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
      process_id = register_OL_subprocess(0, 1);
      logger << LOG_DEBUG << "process_id (bare born)           = " << process_id << endl;
      
      ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s + 1); // ????
      logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s + 1 << " (only to check if ME2 == 0.)" << endl;
      process_id = register_OL_subprocess(0, 11);
      logger << LOG_DEBUG << "process_id (colour correlations) = " << process_id << endl;
      //      }
    }
    else {
      logger << LOG_FATAL << "to be solved " << endl;
      exit(1);
    }
    //    if (process_id == -1 && csi->contribution_order_interference == 1){
  }
  else if (csi->type_contribution == "loop" || 
	   csi->type_contribution == "L2I" || 
	   csi->type_contribution == "L2RT" || 
	   csi->type_contribution == "L2RJ" || 
	   csi->type_contribution == "L2VT" || 
	   csi->type_contribution == "L2VJ" || 
	   csi->type_contribution == "L2VA" || 
	   csi->type_contribution == "L2CT" || 
	   csi->type_contribution == "L2CJ"){
    logger << LOG_DEBUG << "csi->type_contribution = " << csi->type_contribution << endl;
    //  select_OL_born_mode(QCD_order, QEW_order);
    //  ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
    logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
    process_id = register_OL_subprocess(0, 12);
    logger << LOG_DEBUG << "process_id = " << process_id << endl;
  }
  else if (csi->type_contribution == "CA" || 
	   csi->type_contribution == "RCA" || 
	   csi->type_contribution == "RCJ"){
    int type_amplitude = 11;
    if (user.string_value[user.string_map["model"]] == "Bornloop"){type_amplitude = 12;}
    if (csi->type_correction == "QCD"){
      //   ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
      if (csi->contribution_order_interference == 0){
	ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	process_id = register_OL_subprocess(0, type_amplitude); // temporary, to force the use of correct colour correlations !!!
      }
      else if (csi->contribution_order_interference == 1){
	ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	process_id = register_OL_subprocess(0, type_amplitude); // temporary, to force the use of correct colour correlations !!!
	logger << LOG_DEBUG << "process_id (colour correlations) = " << process_id << endl;
	
	ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	process_id = register_OL_subprocess(0, 1); // no colour correlations needed in EW corrections !!! adapt for loop-induced processes !!!
	logger << LOG_DEBUG << "process_id (bare born)           = " << process_id << endl;
	// if (process_id == 2) different libraries are needed for different Born contributions
      }
      else {
	logger << LOG_FATAL << "to be solved " << endl;
	exit(1);
      }
    }
    else if (csi->type_correction == "QEW"){
      ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
      //      ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e - 1);
      logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
      process_id = register_OL_subprocess(0, 1); // temporary, to force the use of correct colour correlations !!!
    }
    else {logger << LOG_FATAL << "Wrong correction type: " << csi->type_contribution << " - " << csi->type_correction << " does not exist." << endl;}
    //    process_id = register_OL_subprocess(0, 1);
    //    process_id = register_OL_subprocess(0, 2);
    //    process_id = register_OL_subprocess(0, 11); // temporary, to force the use of correct colour correlations !!!
  }

  else if (csi->type_contribution == "L2CA"){
    int type_amplitude = 12;
    if (csi->type_correction == "QCD"){
      //   ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
      if (csi->contribution_order_interference == 0){
	ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	process_id = register_OL_subprocess(0, type_amplitude); // temporary, to force the use of correct colour correlations !!!
      }
      /*
      else if (csi->contribution_order_interference == 1){
	ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	process_id = register_OL_subprocess(0, type_amplitude); // temporary, to force the use of correct colour correlations !!!
	logger << LOG_DEBUG << "process_id (colour correlations) = " << process_id << endl;
	
	ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	process_id = register_OL_subprocess(0, type_amplitude); // no colour correlations needed in EW corrections !!! adapt for loop-induced processes !!!
	logger << LOG_DEBUG << "process_id (bare born)           = " << process_id << endl;
	// if (process_id == 2) different libraries are needed for different Born contributions
      }
      */
      else {
	logger << LOG_FATAL << "to be solved " << endl;
	exit(1);
      }
    }
    /*
    else if (csi->type_correction == "QEW"){
      ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
      logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
      process_id = register_OL_subprocess(0, 1); // temporary, to force the use of correct colour correlations !!!
    }
    */
    else {logger << LOG_FATAL << "Wrong correction type: " << csi->type_contribution << " - " << csi->type_correction << " does not exist." << endl;}
  }


  else if (csi->type_contribution == "VA" || 
	   csi->type_contribution == "RVA" || 
	   csi->type_contribution == "RVJ"){
    if (csi->type_correction == "QCD"){
      //   ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
      ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
      logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
      process_id = register_OL_subprocess(0, 11);
    }
    else if (csi->type_correction == "QEW"){
      /////      ol_setparameter_int(stch("ew_renorm"), 1); // needed ???
      ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
      logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
      process_id = register_OL_subprocess(0, 11);
    }
    else if (csi->type_correction == "MIX"){
      /////      ol_setparameter_int(stch("ew_renorm"), 1); // needed ???
      ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
      logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
      process_id = register_OL_subprocess(0, 11);
      if (process_id == -1){
	ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	process_id = register_OL_subprocess(0, 11);
      }
      logger << LOG_DEBUG << "Virtual amplitude initialized.   process_id = " << process_id << endl;
      
      logger << LOG_DEBUG << "(*VA_ioperator).size() = " << setw(15) << (*VA_ioperator).size() << endl;
      
      for (int i_a = 0; i_a < (*VA_ioperator).size(); i_a++){
	logger << LOG_DEBUG << "(*VA_ioperator)[" << i_a << "].size() = " << setw(15) << (*VA_ioperator)[i_a].size() << endl;
	for (int j_a = 0; j_a < (*VA_ioperator)[i_a].size(); j_a++){
	  logger << LOG_DEBUG << "(*VA_ioperator)[" << i_a << "][" << j_a << "].name() = " << setw(15) << (*VA_ioperator)[i_a][j_a].name() << "   type_correction = " << (*VA_ioperator)[i_a][j_a].type_correction() << "   to be processed..." << endl;
	  if ((*VA_ioperator)[i_a][j_a].type_correction() == 1){
	    ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	    logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	    (*VA_ioperator)[i_a][j_a].process_id = register_OL_subprocess(0, 11);
	    if ((*VA_ioperator)[i_a][j_a].process_id == -1){
	      ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	      logger << LOG_DEBUG << "order_ew = " << csi->contribution_order_alpha_e << endl;
	      (*VA_ioperator)[i_a][j_a].process_id = register_OL_subprocess(0, 11);
	    }
	  }
	  if ((*VA_ioperator)[i_a][j_a].type_correction() == 2){
	    ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	    logger << LOG_DEBUG << "order_qcd = " << csi->contribution_order_alpha_s << endl;
	    (*VA_ioperator)[i_a][j_a].process_id = register_OL_subprocess(0, 1);
	  }
	  
	  logger << LOG_DEBUG << "(*VA_ioperator)[" << i_a << "][" << j_a << "].name() = " << setw(15) << (*VA_ioperator)[i_a][j_a].name() << "   type_correction = " << (*VA_ioperator)[i_a][j_a].type_correction() << "   process_id = " << (*VA_ioperator)[i_a][j_a].process_id << endl;
	}
      }
      
      //      ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
      /*
	if (csi->contribution_order_interference == 0){
	ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	}
	else if (csi->contribution_order_interference == 1){
	ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	}
      */
    }
    else {logger << LOG_FATAL << "Wrong correction type: " << csi->type_contribution << " - " << csi->type_correction << " does not exist." << endl;}
    //    process_id = register_OL_subprocess(0, 11);
  }
  
  
  
  else if (csi->type_contribution == "RA" || 
	   csi->type_contribution == "RRA" || 
	   csi->type_contribution == "RRJ"){
    logger << LOG_DEBUG_VERBOSE << "(R)RA:   (*RA_dipole).size() = " << (*RA_dipole).size() << endl;
    int type_amplitude = 1;
    if (user.string_value[user.string_map["model"]] == "Bornloop"){type_amplitude = 12;}
    for (int i_a = 0; i_a < (*RA_dipole).size(); i_a++){
      logger << LOG_DEBUG_VERBOSE << "(R)RA:   i_a = " << i_a << endl;
      if (i_a == 0){
	ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	(*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude);
	logger << LOG_DEBUG << "(*RA_dipole)[0].process_id = " << (*RA_dipole)[0].process_id << endl;
	logger << LOG_DEBUG << "type_amplitude = " << type_amplitude << endl;
	if ((*RA_dipole)[i_a].process_id == -1){
	  ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	  logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	  (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude);
	}
	/*
	  ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	  (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, 1);
	*/
      }
      else {
	if (csi->type_correction == "QCD"){
	  type_amplitude = 11;
	  if (csi->contribution_order_interference == 1){
	    ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	    logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	    (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude); // temporary, to force the use of correct colour correlations !!!
	  }
	  else {
	    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	    logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	    (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude); // temporary, to force the use of correct colour correlations !!!
	  }
	}
	else if (csi->type_correction == "QEW"){
	  ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	  logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	  (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude); // temporary, only Born needed here !!!
	}
	else if (csi->type_correction == "MIX"){
	  if ((*RA_dipole)[i_a].type_correction() == 1){
	    type_amplitude = 11;
	    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	    logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	    //	    (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, 1); // temporary, to force the use of correct colour correlations !!!
	    (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude); // temporary, to force the use of correct colour correlations !!!
	  }
	  else if ((*RA_dipole)[i_a].type_correction() == 2){
	    ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	    logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	    (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude); // temporary, only Born needed here !!!
	  }
	  else {logger << LOG_FATAL << "Wrong correction type: " << (*RA_dipole)[i_a].type_correction() << " does not exist." << endl;}
	  //	  (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, 1);
	}
	//	  (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, 11); // temporary, to force the use of correct colour correlations !!!
	//      select_OL_born_mode(QCD_order - 1, QEW_order);
      }
      //      (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, 2);
      //      (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, 1);
      logger << LOG_DEBUG << "dipole[" << i_a << "].process_id  = " << (*RA_dipole)[i_a].process_id << endl;
      if ((*RA_dipole)[i_a].process_id == -1){logger << LOG_FATAL << "Requested amplitudes are not available." << endl;}
    }
    logger.newLine(LOG_DEBUG);
  }


  else if (csi->type_contribution == "L2RA"){
    logger << LOG_DEBUG_VERBOSE << "(L2RA:   (*RA_dipole).size() = " << (*RA_dipole).size() << endl;
    int type_amplitude = 12;
    for (int i_a = 0; i_a < (*RA_dipole).size(); i_a++){
      logger << LOG_DEBUG_VERBOSE << "L2RA:   i_a = " << i_a << endl;
      if (i_a == 0){
	ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	(*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude);
	logger << LOG_DEBUG << "(*RA_dipole)[0].process_id = " << (*RA_dipole)[0].process_id << endl;
	logger << LOG_DEBUG << "type_amplitude = " << type_amplitude << endl;
	if ((*RA_dipole)[i_a].process_id == -1){
	  ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	  logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	  (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude);
	}
      }
      else {
	if (csi->type_correction == "QCD"){
	  if (csi->contribution_order_interference == 1){
	    ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	    logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	    (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude); // temporary, to force the use of correct colour correlations !!!
	  }
	  else {
	    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	    logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	    (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude); // temporary, to force the use of correct colour correlations !!!
	  }
	}
	/*
	else if (csi->type_correction == "QEW"){
	  ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	  logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	  (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude); // temporary, only Born needed here !!!
	}
	else if (csi->type_correction == "MIX"){
	  if ((*RA_dipole)[i_a].type_correction() == 1){
	    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
	    logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
	    //	    (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, 1); // temporary, to force the use of correct colour correlations !!!
	    (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude); // temporary, to force the use of correct colour correlations !!!
	  }
	  else if ((*RA_dipole)[i_a].type_correction() == 2){
	    ol_setparameter_int(stch("order_qcd"), csi->contribution_order_alpha_s);
	    logger << LOG_DEBUG << "order_qcd   set to   " << csi->contribution_order_alpha_s << endl;
	    (*RA_dipole)[i_a].process_id = register_OL_subprocess(i_a, type_amplitude); // temporary, only Born needed here !!!
	  }
	  else {logger << LOG_FATAL << "Wrong correction type: " << (*RA_dipole)[i_a].type_correction() << " does not exist." << endl;}
	}
	*/
      }
      logger << LOG_DEBUG << "dipole[" << i_a << "].process_id  = " << (*RA_dipole)[i_a].process_id << endl;
      if ((*RA_dipole)[i_a].process_id == -1){logger << LOG_FATAL << "Requested amplitudes are not available." << endl;}
    }
    logger.newLine(LOG_DEBUG);
  }




  
  else if (csi->type_contribution == "CT" ||
	   csi->type_contribution == "CJ"){
    int type_amplitude = 1;
    if (user.string_value[user.string_map["model"]] == "Bornloop"){type_amplitude = 12;}
    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
    logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
    process_id = register_OL_subprocess(0, type_amplitude);
  }
  else if (csi->type_contribution == "CT2" ||
	   csi->type_contribution == "CJ2"){
    //  ol_setparameter_int(stch("polenorm"), 1); // polenorm prescribed by H^(1) construction
    if (QT_finalstate_massive_coloured){
      ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
      logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
      process_id = register_OL_subprocess(0, 11);
      // loop-squared amplitudes are needed to run the double-virtual result to m_inv (why here in CT2 ???):
      process_id = register_OL_subprocess(0, 12);
    }
    else {
      ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
      logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
      process_id = register_OL_subprocess(0, 11);
    }
  }
  else if (csi->type_contribution == "VT" ||
	   csi->type_contribution == "VJ"){
    int type_amplitude = 11;
    if (user.string_value[user.string_map["model"]] == "Bornloop"){type_amplitude = 12;}
    //  ol_setparameter_int(stch("polenorm"), 1); // polenorm prescribed by H^(1) construction
    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
    process_id = register_OL_subprocess(0, type_amplitude);
  }
  else if (csi->type_contribution == "VT2" ||
	   csi->type_contribution == "VJ2"){
    //  ol_setparameter_int(stch("polenorm"), 1); // polenorm prescribed by H^(1) construction
    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
    logger << LOG_DEBUG << "order_ew   set to   " << csi->contribution_order_alpha_e << endl;
    if (QT_finalstate_massive_coloured){
      process_id = register_OL_subprocess(0, 11);
      // loop-squared amplitudes are needed to run the double-virtual result to m_inv:
      int process_id2 = register_OL_subprocess(0, 12);
      if (process_id2 == -1){process_id = -1;}
       // imaginary parts of loop*Born amplitudes (and colour correlators) are needed to run the double-virtual result to m_inv:
     ol_setparameter_string(stch("approx"), "imag");
      int process_id3 = register_OL_subprocess(0, 11);
      if (process_id3 == -1){process_id = -1;}
    }
    else {
      process_id = register_OL_subprocess(0, 11);
    }
  }
  
  else if (csi->type_contribution == "NLL_LO"){
    //  ol_setparameter_int(stch("polenorm"), 1); // polenorm prescribed by H^(1) construction
    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
    process_id = register_OL_subprocess(0, 11);
  }
  else if (csi->type_contribution == "NLL_NLO"){
    //  ol_setparameter_int(stch("polenorm"), 1); // polenorm prescribed by H^(1) construction
    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
    process_id = register_OL_subprocess(0, 11);
  }
  else if (csi->type_contribution == "NNLL_LO"){
    //  ol_setparameter_int(stch("polenorm"), 1); // polenorm prescribed by H^(1) construction
    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
    process_id = register_OL_subprocess(0, 11);
  }
  else if (csi->type_contribution == "NNLL_NLO"){
    //  ol_setparameter_int(stch("polenorm"), 1); // polenorm prescribed by H^(1) construction
    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
    process_id = register_OL_subprocess(0, 11);
  }
  else if (csi->type_contribution == "NNLL_NNLO"){
    //  ol_setparameter_int(stch("polenorm"), 1); // polenorm prescribed by H^(1) construction
    ol_setparameter_int(stch("order_ew"), csi->contribution_order_alpha_e);
    process_id = register_OL_subprocess(0, 11);
  }
  
  if (csi->type_contribution != "RA" && 
      csi->type_contribution != "RRA" && 
      csi->type_contribution != "RRJ" && 
      csi->type_contribution != "L2RA"){
    logger << LOG_DEBUG << "process_id = " << process_id << endl;
    if (process_id == -1){logger << LOG_FATAL << "Requested amplitudes are not available." << endl;}
  }
  
  ol_start();
  //  OLP_PrintParameter("olparameters.txt");
  //  logger << LOG_DEBUG << "1st testpoint" << endl;
  //  OpenLoops_testpoint_pptt_cc();
  //  logger << LOG_DEBUG << "2nd testpoint" << endl;
  //  OpenLoops_testpoint_pptt_cc();
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



int observable_set::register_OL_subprocess(int i_a, int amptype){
  static Logger logger("register_OL_subprocess");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  stringstream temp_processname_ss;
  for (int i_p = 1; i_p < csi->type_parton[i_a].size(); i_p++){
    if (csi->type_parton[i_a][i_p] == 0){temp_processname_ss << 21;}
    // pdf convention:
    // W+ = +24 !!!
    // W- = -24 !!!
    else if (csi->type_parton[i_a][i_p] == 24){temp_processname_ss << -24;}
    else if (csi->type_parton[i_a][i_p] == -24){temp_processname_ss << 24;}

    /// only off-shell Photonen !!! needs to be controlled elsewhere...
    else if (csi->type_parton[i_a][i_p] == 22){temp_processname_ss << 22;}
    else if (csi->type_parton[i_a][i_p] == -22){temp_processname_ss << 22;}
    //    else if (csi->type_parton[i_a][i_p] == 22){temp_processname_ss << -22;}
    ///
    
    else {temp_processname_ss << csi->type_parton[i_a][i_p];}
    if (i_p == 2){temp_processname_ss << " -> ";}
    else if (i_p == csi->type_parton[i_a].size() - 1){}
    else {temp_processname_ss << " ";}
  }
  string temp_processname = temp_processname_ss.str();

  char * char_processname = new char[temp_processname.size() + 1];
  std::copy(temp_processname.begin(), temp_processname.end(), char_processname);
  char_processname[temp_processname.size()] = '\0';

  logger << LOG_DEBUG << "before registration:   processname = " << setw(20) << char_processname << "   amptype = " << amptype << endl;

  int no_reg = ol_register_process(char_processname, amptype);

  logger << LOG_DEBUG << " after registration:   processname = " << setw(20) << char_processname << "   amptype = " << amptype << "   id = " << no_reg << endl;

  if (no_reg == -1){
   stringstream temp_symm_processname_ss;
   vector<int> symm_type_parton = csi->type_parton[i_a];
    for (int i_p = 1; i_p < symm_type_parton.size(); i_p++){
      if   (symm_type_parton[i_p] == 1){symm_type_parton[i_p] = 3;}
      else if (symm_type_parton[i_p] == 2){symm_type_parton[i_p] = 4;}
      else if (symm_type_parton[i_p] == 3){symm_type_parton[i_p] = 1;}
      else if (symm_type_parton[i_p] == 4){symm_type_parton[i_p] = 2;}
      else if (symm_type_parton[i_p] == -1){symm_type_parton[i_p] = -3;}
      else if (symm_type_parton[i_p] == -2){symm_type_parton[i_p] = -4;}
      else if (symm_type_parton[i_p] == -3){symm_type_parton[i_p] = -1;}
      else if (symm_type_parton[i_p] == -4){symm_type_parton[i_p] = -2;}
    }
    for (int i_p = 1; i_p < symm_type_parton.size(); i_p++){
      if (symm_type_parton[i_p] == 0){temp_symm_processname_ss << 21;}
      // pdf convention:
      // W+ = +24 !!!
      // W- = -24 !!!
      else if (symm_type_parton[i_p] == 24){temp_symm_processname_ss << -24;}
      else if (symm_type_parton[i_p] == -24){temp_symm_processname_ss << 24;}

      // temporary !!! Needs to distinguish between ax (-22) and a (22) !!!
      else if (symm_type_parton[i_p] == 22){temp_symm_processname_ss << -22;}

      else {temp_symm_processname_ss << symm_type_parton[i_p];}
      if (i_p == 2){temp_symm_processname_ss << " -> ";}
      else if (i_p == symm_type_parton.size() - 1){}
      else {temp_symm_processname_ss << " ";}
    }
    string temp_symm_processname = temp_symm_processname_ss.str();
    char * char_symm_processname = new char[temp_symm_processname.size() + 1];
    std::copy(temp_symm_processname.begin(), temp_symm_processname.end(), char_symm_processname);
    char_symm_processname[temp_symm_processname.size()] = '\0';

    logger << LOG_DEBUG << "before registration:   symm_processname = " << setw(20) << char_symm_processname << "   amptype = " << amptype << endl;

    no_reg = ol_register_process(char_symm_processname, amptype);

    logger << LOG_DEBUG << " after registration:   symm_processname = " << setw(20) << char_symm_processname << "   amptype = " << amptype << "   id = " << no_reg << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return no_reg;
}

void observable_set::testpoint_from_OL_rambo(){
  static Logger logger("testpoint_from_OL_rambo");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  double energy = oset.E_CMS;
  static int n_momentum = 5 * (n_particle + 2);
  logger << LOG_DEBUG_VERBOSE << "n_momentum = " << n_momentum << endl;
  double *P;
  P = new double[n_momentum];
  ol_phase_space_point(1, E_CMS, P);
  for (int i = 1; i < p_parton[0].size(); i++){
    p_parton[0][i] = fourvector(P[5 * (i - 1)], P[5 * (i - 1) + 1], P[5 * (i - 1) + 2], P[5 * (i - 1) + 3]);
  }
  p_parton[0][0] = p_parton[0][1] + p_parton[0][2];
  delete [] P;
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
