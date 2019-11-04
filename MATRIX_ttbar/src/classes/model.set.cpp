#include "../include/classes.cxx"
//
// constructors
//
model_set::model_set(){};

model_set::model_set(vector<string> & file_input, int N_f, int N_f_active, double scale_ren, int present_type_perturbative_order, vector<string> & collection_type_perturbative_order, vector<string> & _contribution_LHAPDFname, vector<int> & _contribution_LHAPDFsubset, int switch_alpha_rescaling_exp, int switch_costheta_real){

  Logger logger("model_set::model_set - input file");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  logger << LOG_DEBUG << "" << endl;
  //  logger << LOG_DEBUG << "max_perturbative_QCD_order = " << max_perturbative_QCD_order << endl;
  int n_perturbative_order = collection_type_perturbative_order.size();
  logger << LOG_DEBUG << "n_perturbative_order = " << n_perturbative_order << endl;

  contribution_LHAPDFname = _contribution_LHAPDFname;
  contribution_LHAPDFsubset = _contribution_LHAPDFsubset;
  
  char LineBuffer[256];
  logger << LOG_DEBUG << "reading parameters from 'file_model.dat'" << endl;

  empty_double = -1.e100;
  empty_double_complex = double_complex(-1.e100, -1.e100);
  empty_int = -32000;
  empty_string = "empty";

  V_du = empty_double;
  V_su = empty_double;
  V_dc = empty_double;
  V_sc = empty_double;
  V_bu = empty_double;
  V_bc = empty_double;
  V_dt = empty_double;
  V_st = empty_double;
  V_bt = empty_double;

  Q_u = empty_double;
  Q_d = empty_double;
  Q_v = empty_double;
  Q_l = empty_double;
  Iw_u = empty_double;
  Iw_d = empty_double;
  Iw_n = empty_double;
  Iw_l = empty_double;
  
  ew_scheme = empty_int;
  use_cms = empty_int;
  use_adapted_ew_coupling = empty_int;
  
  G_F = empty_double;
  theta_c = empty_double;

  M_W = empty_double;
  M_Z = empty_double;
  M_t = empty_double;
  M_b = empty_double;
  M_c = empty_double;
  M_s = empty_double;
  M_u = empty_double;
  M_d = empty_double;
  M_e = empty_double;
  M_mu = empty_double;
  M_tau = empty_double;
  M_ve = empty_double;
  M_vm = empty_double;
  M_vt = empty_double;
  M_H = empty_double;

  reg_Gamma_Z = empty_double;
  reg_Gamma_W = empty_double;
  reg_Gamma_H = empty_double;
  reg_Gamma_t = empty_double;

  map_Gamma_Z = empty_double;
  map_Gamma_W = empty_double;
  map_Gamma_H = empty_double;
  map_Gamma_t = empty_double;

  alpha_e = empty_double;
  e = empty_double;

  alpha_e_0 = empty_double;
  alpha_e_MZ = empty_double;
  alpha_e_Gmu = empty_double;

  cos_w = empty_double;
  cos2_w = empty_double;
  sin2_w = empty_double;
  sin_w = empty_double;

  Cplus_Zuu = empty_double;
  Cplus_Zdd = empty_double;
  Cminus_Zuu = empty_double;
  Cminus_Zdd = empty_double;
  Cplus_Znn = empty_double;
  Cplus_Zee = empty_double;
  Cminus_Znn = empty_double;
  Cminus_Zee = empty_double;
  C_ZWminusWplus = empty_double;

  Cplus_Auu = empty_double;
  Cplus_Add = empty_double;
  Cminus_Auu = empty_double;
  Cminus_Add = empty_double;
  Cplus_Ann = empty_double;
  Cplus_Aee = empty_double;
  Cminus_Ann = empty_double;
  Cminus_Aee = empty_double;
  C_AWminusWplus = empty_double;
  Cminus_W = empty_double;

  ccos_w = empty_double_complex;
  ccos2_w = empty_double_complex;
  csin2_w = empty_double_complex;
  csin_w = empty_double_complex;
  cCplus_Zuu = empty_double_complex;
  cCplus_Zdd = empty_double_complex;
  cCminus_Zuu = empty_double_complex;
  cCminus_Zdd = empty_double_complex;
  cCplus_Znn = empty_double_complex;
  cCplus_Zee = empty_double_complex;
  cCminus_Znn = empty_double_complex;
  cCminus_Zee = empty_double_complex;

  cC_ZWminusWplus = empty_double_complex;
  
  //  cCplus_Auu = empty_double_complex;
  //  cCplus_Add = empty_double_complex;
  //  cCminus_Auu = empty_double_complex;
  //  cCminus_Add = empty_double_complex;
  //  cCplus_Ann = empty_double_complex;
  //  cCplus_Aee = empty_double_complex;
  //  cCminus_Ann = empty_double_complex;
  //  cCminus_Aee = empty_double_complex;
  //  cC_AWminusWplus = empty_double_complex;
  
  cCminus_W = empty_double_complex;

  CKM_matrix = empty_string;

  //  conversion_to_pole_mass = 0;
  
  vGamma_W_order.resize(n_perturbative_order, empty_int);
  valpha_S_W_order.resize(n_perturbative_order, empty_int);
  valpha_S_W.resize(n_perturbative_order, empty_double);
  valpha_S_W_scale.resize(n_perturbative_order, empty_double);
  vGamma_W.resize(n_perturbative_order, empty_double);
  vGamma_Wud.resize(n_perturbative_order, empty_double);
  vGamma_Wlv.resize(n_perturbative_order, empty_double);
  vBR_Wud.resize(n_perturbative_order, empty_double);
  vBR_Wlv.resize(n_perturbative_order, empty_double);
  //  double alpha_S_W;

  vGamma_Z_order.resize(n_perturbative_order, empty_int);
  valpha_S_Z_order.resize(n_perturbative_order, empty_int);
  valpha_S_Z.resize(n_perturbative_order, empty_double);
  valpha_S_Z_scale.resize(n_perturbative_order, empty_double);
  vGamma_Z.resize(n_perturbative_order, empty_double);
  vGamma_Zuu.resize(n_perturbative_order, empty_double);
  vGamma_Zdd.resize(n_perturbative_order, empty_double);
  vGamma_Zvv.resize(n_perturbative_order, empty_double);
  vGamma_Zll.resize(n_perturbative_order, empty_double);
  vBR_Zuu.resize(n_perturbative_order, empty_double);
  vBR_Zdd.resize(n_perturbative_order, empty_double);
  vBR_Zvv.resize(n_perturbative_order, empty_double);
  vBR_Zll.resize(n_perturbative_order, empty_double);
  //  double alpha_S_Z;

  vGamma_H.resize(n_perturbative_order, empty_double);


  vGamma_t_order.resize(n_perturbative_order, empty_int);
  valpha_S_t_order.resize(n_perturbative_order, empty_int);
  valpha_S_t.resize(n_perturbative_order, empty_double);
  valpha_S_t_scale.resize(n_perturbative_order, empty_double);
  vGamma_t.resize(n_perturbative_order, empty_double);
  vGamma_tWb.resize(n_perturbative_order, empty_double);
  vBR_tWb.resize(n_perturbative_order, empty_double);


  //  double alpha_S_t;

  vector<string> readin;

  // replace the following by loop over 'file_input' !!!

  string filename;
  filename = "../file_model.dat";
  ifstream in_model_default_org(filename.c_str());
  while (in_model_default_org.getline(LineBuffer, 256)){readin.push_back(LineBuffer);}
  in_model_default_org.close();

  filename = "../../../../file_model.dat";
  ifstream in_model_default_gen(filename.c_str());
  while (in_model_default_gen.getline(LineBuffer, 256)){readin.push_back(LineBuffer);}
  in_model_default_gen.close();

  filename = "../../../file_model.dat";
  ifstream in_model_default_posm(filename.c_str());
  while (in_model_default_posm.getline(LineBuffer, 256)){readin.push_back(LineBuffer);}
  in_model_default_posm.close();

  filename = "../../file_model.dat";
  ifstream in_model_default_co(filename.c_str());
  while (in_model_default_co.getline(LineBuffer, 256)){readin.push_back(LineBuffer);}
  in_model_default_co.close();

  filename = "../file_model.dat";
  ifstream in_model_default_cc(filename.c_str());
  while (in_model_default_cc.getline(LineBuffer, 256)){readin.push_back(LineBuffer);}
  in_model_default_cc.close();

  filename = "file_model.dat";
  ifstream in_model_default_con(filename.c_str());
  while (in_model_default_con.getline(LineBuffer, 256)){readin.push_back(LineBuffer);}
  in_model_default_con.close();

  if (readin.size() == 0) {
    logger << LOG_ERROR << "file_model.dat could not be read in" << endl;
    assert(readin.size() > 0);
  }

  vector<string> user_variable;
  vector<string> user_value;
  get_simple_userinput_from_readin(user_variable, user_value, readin);


  // ??? meaning of QCD_order here ???
  int QCD_order = -1;

  vector<string> manipulate_variable;
  vector<string> manipulate_variable_counter;
  vector<string> manipulate_value;
  vector<string> det_theta_w;
  vector<string> det_M_gauge;
  vector<string> det_EW_coupling;
  for (int i = 0; i < user_variable.size(); i++){

    if (user_variable[i] == "ew_scheme"){ew_scheme = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "use_cms"){use_cms = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "use_adapted_ew_coupling"){use_adapted_ew_coupling = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "G_F"){G_F = atof(user_value[i].c_str());
      det_EW_coupling.push_back(user_variable[i]);}
    else if (user_variable[i] == "sin_w"){sin_w = atof(user_value[i].c_str());
      det_theta_w.push_back(user_variable[i]);}
    else if (user_variable[i] == "sin2_w"){sin2_w = atof(user_value[i].c_str());
      det_theta_w.push_back(user_variable[i]);}
    else if (user_variable[i] == "cos2_w"){cos2_w = atof(user_value[i].c_str());
      det_theta_w.push_back(user_variable[i]);}
    else if (user_variable[i] == "cos_w"){cos_w = atof(user_value[i].c_str());
      det_theta_w.push_back(user_variable[i]);}
    // alpha_e  is not allowed as a direct input parameter: use
    // - alpha_e_0 (with ew_scheme = 0),
    // - alpha_e_Gmu (with ew_scheme = 1),
    // - alpha_e_MZ (with ew_scheme = 2).
    //    else if (user_variable[i] == "alpha_e"){alpha_e = atof(user_value[i].c_str());
    //      det_EW_coupling.push_back(user_variable[i]);}
    //    else if (user_variable[i] == "e"){alpha_e_0 = pow(atof(user_value[i].c_str()), 2) / (4. * pi);}
    //    else if (user_variable[i] == "e"){e = atof(user_value[i].c_str());
    //      det_EW_coupling.push_back(user_variable[i]);}
    else if (user_variable[i] == "M_W"){M_W = atof(user_value[i].c_str()); 
      det_M_gauge.push_back(user_variable[i]);}
    else if (user_variable[i] == "M_Z"){M_Z = atof(user_value[i].c_str());
      det_M_gauge.push_back(user_variable[i]);}
    else if (user_variable[i] == "alpha_e_0"){alpha_e_0 = atof(user_value[i].c_str());}
    else if (user_variable[i] == "1/alpha_e_0"){alpha_e_0 = 1. / atof(user_value[i].c_str());}
    else if (user_variable[i] == "e_0"){alpha_e_0 = pow(atof(user_value[i].c_str()), 2) / (4. * pi);}
    else if (user_variable[i] == "alpha_e_MZ"){alpha_e_MZ = atof(user_value[i].c_str());}
    else if (user_variable[i] == "1/alpha_e_MZ"){alpha_e_MZ = 1. / atof(user_value[i].c_str());}
    else if (user_variable[i] == "e_MZ"){alpha_e_MZ = pow(atof(user_value[i].c_str()), 2) / (4. * pi);}
    else if (user_variable[i] == "alpha_e_Gmu"){alpha_e_Gmu = atof(user_value[i].c_str());}
    else if (user_variable[i] == "1/alpha_e_Gmu"){alpha_e_Gmu = 1. / atof(user_value[i].c_str());}
    else if (user_variable[i] == "e_Gmu"){alpha_e_Gmu = pow(atof(user_value[i].c_str()), 2) / (4. * pi);}

    else if (user_variable[i] == "M_H"){M_H = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_t"){M_t = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_b"){M_b = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_c"){M_c = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_s"){M_s = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_u"){M_u = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_d"){M_d = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_e"){M_e = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_mu"){M_mu = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_tau"){M_tau = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_ve"){M_ve = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_vm"){M_vm = atof(user_value[i].c_str());}
    else if (user_variable[i] == "M_vt"){M_vt = atof(user_value[i].c_str());}

    else if (user_variable[i] == "reg_Gamma_Z"){reg_Gamma_Z = atof(user_value[i].c_str());} 
    else if (user_variable[i] == "reg_Gamma_W"){reg_Gamma_W = atof(user_value[i].c_str());}
    else if (user_variable[i] == "reg_Gamma_H"){reg_Gamma_H = atof(user_value[i].c_str());} 
    else if (user_variable[i] == "reg_Gamma_t"){reg_Gamma_t = atof(user_value[i].c_str());}

    else if (user_variable[i] == "map_Gamma_Z"){map_Gamma_Z = atof(user_value[i].c_str());}
    else if (user_variable[i] == "map_Gamma_W"){map_Gamma_W = atof(user_value[i].c_str());}
    else if (user_variable[i] == "map_Gamma_H"){map_Gamma_H = atof(user_value[i].c_str());}
    else if (user_variable[i] == "map_Gamma_t"){map_Gamma_t = atof(user_value[i].c_str());}

    else if (user_variable[i] == "CKM_matrix"){CKM_matrix = user_value[i];}
    else if (user_variable[i] == "theta_c"){theta_c = atof(user_value[i].c_str());}

    else if (user_variable[i] == "V_du"){V_du = atof(user_value[i].c_str());}
    else if (user_variable[i] == "V_su"){V_su = atof(user_value[i].c_str());}
    else if (user_variable[i] == "V_bu"){V_bu = atof(user_value[i].c_str());}
    else if (user_variable[i] == "V_dc"){V_dc = atof(user_value[i].c_str());}
    else if (user_variable[i] == "V_sc"){V_sc = atof(user_value[i].c_str());}
    else if (user_variable[i] == "V_bc"){V_bc = atof(user_value[i].c_str());}
    else if (user_variable[i] == "V_dt"){V_dt = atof(user_value[i].c_str());}
    else if (user_variable[i] == "V_st"){V_st = atof(user_value[i].c_str());}
    else if (user_variable[i] == "V_bt"){V_bt = atof(user_value[i].c_str());}

    else if (user_variable[i] == "Q_u"){Q_u = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Q_d"){Q_d = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Q_v"){Q_v = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Q_l"){Q_l = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Iw_u"){Iw_u = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Iw_d"){Iw_d = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Iw_n"){Iw_n = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Iw_l"){Iw_l = atof(user_value[i].c_str());}

    else if (user_variable[i] == "Cplus_Zuu"){Cplus_Zuu = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cplus_Zdd"){Cplus_Zdd = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cminus_Zuu"){Cminus_Zuu = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cminus_Zdd"){Cminus_Zdd = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cplus_Znn"){Cplus_Znn = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cplus_Zee"){Cplus_Zee = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cminus_Znn"){Cminus_Znn = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cminus_Zee"){Cminus_Zee = atof(user_value[i].c_str());}
    else if (user_variable[i] == "C_ZWminusWplus"){C_ZWminusWplus = atof(user_value[i].c_str());}

    else if (user_variable[i] == "Cplus_Auu"){Cplus_Auu = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cplus_Add"){Cplus_Add = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cminus_Auu"){Cminus_Auu = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cminus_Add"){Cminus_Add = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cplus_Ann"){Cplus_Ann = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cplus_Aee"){Cplus_Aee = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cminus_Ann"){Cminus_Ann = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cminus_Aee"){Cminus_Aee = atof(user_value[i].c_str());}
    else if (user_variable[i] == "C_AWminusWplus"){C_AWminusWplus = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Cminus_W"){Cminus_W = atof(user_value[i].c_str());}

    //    else if (user_variable[i] == "conversion_to_pole_mass"){conversion_to_pole_mass = atoi(user_value[i].c_str());}





    else if (user_variable[i] == "contribution"){
      QCD_order = -1;
      /*
      for (int i_o = 0; i_o < collection_type_perturbative_order.size(); i_o++){
	if (user_value[i] == collection_type_perturbative_order[i_o]){QCD_order = i_o; break;}
      }
      if (QCD_order == -1){logger << LOG_FATAL << "No pdfset specified for this order!" << endl; exit(1);}
      */

      if      (user_value[i] == "LO"){QCD_order = 0;}
      else if (user_value[i] == "NLO"){QCD_order = 1;}
      else if (user_value[i] == "NNLO"){QCD_order = 2;}// ???
      else if (user_value[i] == "NNNLO"){QCD_order = 3;}// ???
            
      else {logger << LOG_DEBUG << "QCD_order = " << QCD_order << endl;}
    }

    else if (user_variable[i] == "Gamma_W"){vGamma_W[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Gamma_Wlv"){vGamma_Wlv[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Gamma_Wud"){vGamma_Wud[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "BR_Wlv"){vBR_Wlv[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "BR_Wud"){vBR_Wud[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Gamma_W_order"){vGamma_W_order[QCD_order] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "alpha_S_W_order"){valpha_S_W_order[QCD_order] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "alpha_S_W_scale"){valpha_S_W_scale[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "alpha_S_W"){valpha_S_W[QCD_order] = atof(user_value[i].c_str());}

    else if (user_variable[i] == "Gamma_Z"){vGamma_Z[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Gamma_Zll"){vGamma_Zll[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Gamma_Zvv"){vGamma_Zvv[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Gamma_Zdd"){vGamma_Zdd[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Gamma_Zuu"){vGamma_Zuu[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "BR_Zll"){vBR_Zll[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "BR_Zvv"){vBR_Zvv[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "BR_Zdd"){vBR_Zdd[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "BR_Zuu"){vBR_Zuu[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Gamma_Z_order"){vGamma_Z_order[QCD_order] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "alpha_S_Z_order"){valpha_S_Z_order[QCD_order] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "alpha_S_Z_scale"){valpha_S_Z_scale[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "alpha_S_Z"){valpha_S_Z[QCD_order] = atof(user_value[i].c_str());}

    else if (user_variable[i] == "Gamma_H"){vGamma_H[QCD_order] = atof(user_value[i].c_str());}

    else if (user_variable[i] == "Gamma_t"){
      vGamma_t[QCD_order] = atof(user_value[i].c_str());
      logger << LOG_DEBUG_VERBOSE << "vGamma_t[" << QCD_order << "] = " << vGamma_t[QCD_order] << endl;
    }
    else if (user_variable[i] == "Gamma_tWb"){vGamma_tWb[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "BR_tWb"){vBR_tWb[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "Gamma_t_order"){vGamma_t_order[QCD_order] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "alpha_S_t_order"){valpha_S_t_order[QCD_order] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "alpha_S_t_scale"){valpha_S_t_scale[QCD_order] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "alpha_S_t"){valpha_S_t[QCD_order] = atof(user_value[i].c_str());}
    
  

    //    else if (user_variable[i] == "M_Z"){M_Z = atof(user_value[i].c_str());}
    else if (user_variable[i] == "N_f_active"){N_f_active = atoi(user_value[i].c_str());}
  }

  /////////////////////////////////////
  //  end of file_model.dat read-in  //
  /////////////////////////////////////

  logger << LOG_DEBUG << "manipulate_variable.size() = " << manipulate_variable.size() << endl;
  logger << LOG_DEBUG << "manipulate_value.size() = " << manipulate_value.size() << endl;


  logger << LOG_DEBUG << "The following SM parameters and masses were not specified and thus set to the following standard values:" << endl;

  if (M_W == empty_double){M_W = 80.399;}
  if (M_Z == empty_double){M_Z = 91.1876;}
  if (M_H == empty_double){M_H = 126.;}
  if (M_t == empty_double){M_t = 172.0E+00;}
  if (M_b == empty_double){M_b = 4.7E+00;}
  if (M_c == empty_double){M_c = 0.;}
  if (M_s == empty_double){M_s = 0.;}
  if (M_u == empty_double){M_u = 0.;}
  if (M_d == empty_double){M_d = 0.;}
  if (M_e == empty_double){M_e = 0.;}
  if (M_mu == empty_double){M_mu = 0.;}
  if (M_tau == empty_double){M_tau = 0.;}
  if (M_ve == empty_double){M_ve = 0.;}
  if (M_vm == empty_double){M_vm = 0.;}
  if (M_vt == empty_double){M_vt = 0.;}

  /*
  double M_W_pole = M_W;
  if (conversion_to_pole_mass){
    double M_W_pole = M_W / sqrt(1. + pow(Gamma_W / M_W, 2));
    
  }
  */


  // default values:
  if (ew_scheme == empty_int){ew_scheme = 1;}
  if (use_cms == empty_int){use_cms = 1;}
  if (use_adapted_ew_coupling == empty_int){use_adapted_ew_coupling = -1;}

  // Here one could reactivate the input parameter  alpha_e / 1/alpha_e !


  determine_sintheta_real(det_M_gauge, det_theta_w);
  determine_EWcouplings_real();
  
  // Calculation of Gmu-scheme couplings and of Gamma_W/Z based on them:
  int save_use_cms = use_cms;
  if (use_cms == 1 || use_cms == 3){use_cms = 2;} // otherwise Gamma_W/Z and alpha_e_Gmu evaluation becomes circular
  determine_alpha_e_Gmu(det_EW_coupling, switch_alpha_rescaling_exp);
  calculate_Gamma_W(present_type_perturbative_order);
  calculate_Gamma_Z(present_type_perturbative_order);
  use_cms = save_use_cms;
  if (use_cms == 1 || use_cms == 3){determine_alpha_e_Gmu(det_EW_coupling, switch_alpha_rescaling_exp);}
  calculate_Gamma_H(present_type_perturbative_order);
  calculate_Gamma_t(present_type_perturbative_order);
  fill_particle_mass_vector();

  if (ew_scheme == 0){
    if (alpha_e_0 == empty_double){alpha_e_0 = 1. / 137.03599907399999;}
    alpha_e = alpha_e_0;
    logger << LOG_DEBUG << "alpha_e  is set to chosen value of  alpha_e_0 = " << setprecision(15) << setw(23) << alpha_e_0 << " = 1 / " << setprecision(15) << setw(23) <<  alpha_e_0 << " ." << endl;
  }
  else if (ew_scheme == 1 || ew_scheme == -1){
    alpha_e = alpha_e_Gmu;
    if (use_cms == 3){use_cms = 2;}
    logger << LOG_DEBUG << "alpha_e  is set to chosen value of  alpha_e_Gmu = " << setprecision(15) << setw(23) << alpha_e_Gmu << " = 1 / " << setprecision(15) << setw(23) <<  alpha_e_Gmu << " ." << endl;
  }
  /*
  else if (ew_scheme == -1){
    alpha_e = alpha_e_Gmu;
    logger << LOG_DEBUG << "alpha_e  is set to chosen value of  alpha_e_Gmu = " << setprecision(15) << setw(23) << alpha_e_Gmu << " = 1 / " << setprecision(15) << setw(23) <<  alpha_e_Gmu << " ." << endl;
  }
  */
  else if (ew_scheme == 2){
    if (alpha_e_MZ == empty_double){alpha_e_MZ = 1. / 128.;}
    alpha_e = alpha_e_MZ;
    logger << LOG_DEBUG << "alpha_e  is set to chosen value of  alpha_e_MZ = " << setprecision(15) << setw(23) << alpha_e_MZ << " = 1 / " << setprecision(15) << setw(23) << alpha_e_MZ << " ." << endl;
  }
  else {logger << LOG_ERROR << "Invalid value:  ew_scheme = " << ew_scheme << endl; exit(1);}

  e = sqrt(4 * pi * alpha_e);
  determine_EWcouplings_complex(switch_costheta_real);

  logger << LOG_INFO << "ew_scheme = " << ew_scheme << endl;
  logger << LOG_INFO << "use_adapted_ew_coupling = " << use_adapted_ew_coupling << endl;
  if (ew_scheme == use_adapted_ew_coupling){use_adapted_ew_coupling = -1;}  // No modification of alpha_e factor required
  logger << LOG_INFO << "use_adapted_ew_coupling = " << use_adapted_ew_coupling << endl;
   
  /*
  logger << LOG_DEBUG_VERBOSE << "alpha_e     = " << alpha_e << endl;
  logger << LOG_DEBUG_VERBOSE << "e           = " << e << endl;
  logger << LOG_DEBUG_VERBOSE << "G_F         = " << G_F << endl;
  logger << LOG_DEBUG_VERBOSE << "det_EW_coupling.size() = " << det_EW_coupling.size() << endl;

  determine_alpha_e_Gmu(det_EW_coupling, switch_alpha_rescaling_exp);
  */

  // At present, OpenLoops does support non-trivial CKM matrices only for selected amplitudes.
  // CKM_matrix must be set to "trivial" for the remaining ones.
  determine_CKM_matrix();

  e_pow.resize(11);
  for (int i = 0; i < e_pow.size(); i++){e_pow[i] = pow(e, i);}



  /*
  logger << LOG_INFO << "alpha_e = " << setprecision(15) << alpha_e << endl;
  if (use_cms == 1){alpha_e = abs(sqrt2 * pow(csin_w, 2) * pow(cM_W, 2) * G_F / pi);}
  logger << LOG_INFO << "use_cms = " << use_cms << endl;
  logger << LOG_INFO << "alpha_e = " << setprecision(15) << alpha_e << endl;
  */



  output_model_file(switch_costheta_real);
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void model_set::calculate_Gamma_W(int QCD_order){
  Logger logger("model_set::calculate_Gamma_W");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // determine at which order in perturbation theory (QCD) Gamma_W is calculated
  //  logger << LOG_DEBUG << "pdf_set = " << pdf_set << endl;
  if (vGamma_W_order[QCD_order] == empty_int){vGamma_W_order[QCD_order] = QCD_order;}
  int switch_W = vGamma_W_order[QCD_order];
  logger << LOG_DEBUG << "switch_W = " << switch_W << endl;
  logger << LOG_DEBUG << "QCD_order = " << QCD_order << endl;
  // determine alpha_S_W used in Gamma_W
  if (valpha_S_W[QCD_order] == empty_double){
    if (valpha_S_W_scale[QCD_order] == empty_double){valpha_S_W_scale[QCD_order] = M_W;}
    if (valpha_S_W_order[QCD_order] == empty_int){valpha_S_W_order[QCD_order] = switch_W;}
    logger << LOG_DEBUG << "LHAPDF:   QCD_order = " << QCD_order << endl;
    logger << LOG_DEBUG << "LHAPDF:   valpha_S_W_order[QCD_order = " << QCD_order << "] = " << valpha_S_W_order[QCD_order] << endl;
    logger << LOG_DEBUG << "LHAPDF:   contribution_LHAPDFsubset[valpha_S_W_order[QCD_order = " << QCD_order << "] = " << valpha_S_W_order[QCD_order] << "] = " << contribution_LHAPDFsubset[valpha_S_W_order[QCD_order]] << endl;
    initialization_LHAPDF(contribution_LHAPDFname[valpha_S_W_order[QCD_order]], contribution_LHAPDFsubset[valpha_S_W_order[QCD_order]]);
    alpha_S_W = LHAPDF::alphasPDF(valpha_S_W_scale[QCD_order]);
    logger << LOG_DEBUG << "LHAPDF: alpha_S_W = " << setw(25) << setprecision(16) << alpha_S_W << endl;
  }
  else {alpha_S_W = valpha_S_W[QCD_order];}

  // determine partial width at the chosen perturbative order with the chosen alpha_S
  if (switch_W == 0){
    if (vGamma_Wud[QCD_order] == empty_double){vGamma_Wud[QCD_order] = sqrt2 / (12. * pi) * G_F * pow(M_W, 3) * 3.;}
    if (vGamma_Wlv[QCD_order] == empty_double){vGamma_Wlv[QCD_order] = sqrt2 / (12. * pi) * G_F * pow(M_W, 3);}
  }
  else if (switch_W == 1){
    if (vGamma_Wud[QCD_order] == empty_double){vGamma_Wud[QCD_order] = sqrt2 / (12. * pi) * G_F * pow(M_W, 3) * 3. * (1. + alpha_S_W / pi);}
    if (vGamma_Wlv[QCD_order] == empty_double){vGamma_Wlv[QCD_order] = sqrt2 / (12. * pi) * G_F * pow(M_W, 3);}
  }
  else if (switch_W == 2){
    // temporarily !!! should be replaced by NNLO QCD value !!!
    if (vGamma_Wud[QCD_order] == empty_double){vGamma_Wud[QCD_order] = sqrt2 / (12. * pi) * G_F * pow(M_W, 3) * 3. * (1. + alpha_S_W / pi);}
    if (vGamma_Wlv[QCD_order] == empty_double){vGamma_Wlv[QCD_order] = sqrt2 / (12. * pi) * G_F * pow(M_W, 3);}
  }
  // determine total width and branching ratios from the partial widths
  if (vGamma_W[QCD_order] == empty_double){vGamma_W[QCD_order] = 2 * vGamma_Wud[QCD_order] + 3 * vGamma_Wlv[QCD_order];}
  if (vBR_Wlv[QCD_order] == empty_double){vBR_Wlv[QCD_order] = vGamma_Wlv[QCD_order] / vGamma_W[QCD_order];}
  if (vBR_Wud[QCD_order] == empty_double){vBR_Wud[QCD_order] = vGamma_Wud[QCD_order] / vGamma_W[QCD_order];}

  Gamma_W = vGamma_W[QCD_order];
  Gamma_Wlv = vGamma_Wlv[QCD_order];
  Gamma_Wud = vGamma_Wud[QCD_order];
  BR_Wlv = vBR_Wlv[QCD_order];
  BR_Wud = vBR_Wud[QCD_order];

  cM2_W = pow(M_W, 2) - ri * M_W * Gamma_W;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void model_set::calculate_Gamma_Z(int QCD_order){
  Logger logger("model_set::calculate_Gamma_Z");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // determine at which order in perturbation theory (QCD) Gamma_Z is calculated
  //  logger << LOG_DEBUG << "pdf_set = " << pdf_set << endl;
  if (vGamma_Z_order[QCD_order] == empty_int){vGamma_Z_order[QCD_order] = QCD_order;}
  int switch_Z = vGamma_Z_order[QCD_order];
  logger << LOG_DEBUG << "switch_Z = " << switch_Z << endl;
  // determine alpha_S_Z used in Gamma_Z
  if (valpha_S_Z[QCD_order] == empty_double){
    if (valpha_S_Z_scale[QCD_order] == empty_double){valpha_S_Z_scale[QCD_order] = M_Z;}
    if (valpha_S_Z_order[QCD_order] == empty_int){valpha_S_Z_order[QCD_order] = switch_Z;}
    //    alpha_S_Z = calc_alpha_S(pdf_set, N_f, valpha_S_Z_order[QCD_order], valpha_S_Z_scale[QCD_order]);
    //    cout << "old: alpha_S_Z = " << setw(25) << setprecision(16) << alpha_S_Z << endl;
    logger << LOG_DEBUG << "LHAPDF:   QCD_order = " << QCD_order << endl;
    logger << LOG_DEBUG << "LHAPDF:   valpha_S_Z_order[QCD_order = " << QCD_order << "] = " << valpha_S_Z_order[QCD_order] << endl;
    logger << LOG_DEBUG << "LHAPDF:   contribution_LHAPDFsubset[valpha_S_Z_order[QCD_order = " << QCD_order << "] = " << valpha_S_Z_order[QCD_order] << "] = " << contribution_LHAPDFsubset[valpha_S_Z_order[QCD_order]] << endl;
    initialization_LHAPDF(contribution_LHAPDFname[valpha_S_Z_order[QCD_order]], contribution_LHAPDFsubset[valpha_S_Z_order[QCD_order]]);
    alpha_S_Z = LHAPDF::alphasPDF(valpha_S_Z_scale[QCD_order]);
    logger << LOG_DEBUG << "LHAPDF:   alpha_S_Z = " << setw(25) << setprecision(16) << alpha_S_Z << endl;
  }
  else {alpha_S_Z = valpha_S_Z[QCD_order];}

  // determine partial width at the chosen perturbative order with the chosen alpha_S
  if (switch_Z == 0){
    if (vGamma_Zll[QCD_order] == empty_double){vGamma_Zll[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Zee * Cplus_Zee + Cminus_Zee * Cminus_Zee);}
    if (vGamma_Zvv[QCD_order] == empty_double){vGamma_Zvv[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Znn * Cplus_Znn + Cminus_Znn * Cminus_Znn);}
    if (vGamma_Zdd[QCD_order] == empty_double){vGamma_Zdd[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Zdd * Cplus_Zdd + Cminus_Zdd * Cminus_Zdd) * 3.;}
    if (vGamma_Zuu[QCD_order] == empty_double){vGamma_Zuu[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Zuu * Cplus_Zuu + Cminus_Zuu * Cminus_Zuu) * 3.;}
  }
  else if (switch_Z == 1){
    if (vGamma_Zll[QCD_order] == empty_double){vGamma_Zll[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Zee * Cplus_Zee + Cminus_Zee * Cminus_Zee);}
    if (vGamma_Zvv[QCD_order] == empty_double){vGamma_Zvv[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Znn * Cplus_Znn + Cminus_Znn * Cminus_Znn);}
    if (vGamma_Zdd[QCD_order] == empty_double){vGamma_Zdd[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Zdd * Cplus_Zdd + Cminus_Zdd * Cminus_Zdd) * 3. * (1. + alpha_S_Z / pi);}
    if (vGamma_Zuu[QCD_order] == empty_double){vGamma_Zuu[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Zuu * Cplus_Zuu + Cminus_Zuu * Cminus_Zuu) * 3. * (1. + alpha_S_Z / pi);}
  }
  else if (switch_Z == 2){
    // temporarily !!! should be replaced by NNLO QCD value !!!
    if (vGamma_Zll[QCD_order] == empty_double){vGamma_Zll[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Zee * Cplus_Zee + Cminus_Zee * Cminus_Zee);}
    if (vGamma_Zvv[QCD_order] == empty_double){vGamma_Zvv[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Znn * Cplus_Znn + Cminus_Znn * Cminus_Znn);}
    if (vGamma_Zdd[QCD_order] == empty_double){vGamma_Zdd[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Zdd * Cplus_Zdd + Cminus_Zdd * Cminus_Zdd) * 3. * (1. + alpha_S_Z / pi);}
    if (vGamma_Zuu[QCD_order] == empty_double){vGamma_Zuu[QCD_order] = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * pow(M_W, 2) * (Cplus_Zuu * Cplus_Zuu + Cminus_Zuu * Cminus_Zuu) * 3. * (1. + alpha_S_Z / pi);}
  }
  // determine total width and branching ratios from the partial widths
  if (vGamma_Z[QCD_order] == empty_double){vGamma_Z[QCD_order] = 2 * vGamma_Zuu[QCD_order] + 3 * vGamma_Zdd[QCD_order] + 3 * vGamma_Zvv[QCD_order] + 3 * vGamma_Zll[QCD_order];}
  if (vBR_Zll[QCD_order] == empty_double){vBR_Zll[QCD_order] = vGamma_Zll[QCD_order] / vGamma_Z[QCD_order];}
  if (vBR_Zvv[QCD_order] == empty_double){vBR_Zvv[QCD_order] = vGamma_Zvv[QCD_order] / vGamma_Z[QCD_order];}
  if (vBR_Zdd[QCD_order] == empty_double){vBR_Zdd[QCD_order] = vGamma_Zdd[QCD_order] / vGamma_Z[QCD_order];}
  if (vBR_Zuu[QCD_order] == empty_double){vBR_Zuu[QCD_order] = vGamma_Zuu[QCD_order] / vGamma_Z[QCD_order];}

  Gamma_Z = vGamma_Z[QCD_order];
  Gamma_Zuu = vGamma_Zuu[QCD_order];
  Gamma_Zdd = vGamma_Zdd[QCD_order];
  Gamma_Zvv = vGamma_Zvv[QCD_order];
  Gamma_Zll = vGamma_Zll[QCD_order];
  BR_Zdd = vBR_Zdd[QCD_order];
  BR_Zuu = vBR_Zuu[QCD_order];
  BR_Zll = vBR_Zll[QCD_order];
  BR_Zvv = vBR_Zvv[QCD_order];

  cM2_Z = pow(M_Z, 2) - ri * M_Z * Gamma_Z;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void model_set::calculate_Gamma_H(int QCD_order){
  Logger logger("model_set::calculate_Gamma_H");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (vGamma_H[QCD_order] == empty_double){vGamma_H[QCD_order] = 0.017;}
  Gamma_H = vGamma_H[QCD_order];

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void model_set::calculate_Gamma_t(int QCD_order){
  Logger logger("model_set::calculate_Gamma_t");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (vGamma_t[QCD_order] == empty_double && vGamma_tWb[QCD_order] == empty_double){
    // determine at which order in perturbation theory (QCD) Gamma_t is calculated
    if (vGamma_t_order[QCD_order] == empty_int){vGamma_t_order[QCD_order] = QCD_order;}
    int switch_t = vGamma_t_order[QCD_order];
    logger << LOG_DEBUG << "switch_t = " << switch_t << endl;
    logger << LOG_DEBUG << "vGamma_t_order[" << QCD_order << "] = " << vGamma_t_order[QCD_order] << endl;
    logger << LOG_DEBUG << "valpha_S_t_order[" << QCD_order << "] = " << valpha_S_t_order[QCD_order] << endl;
    logger << LOG_DEBUG << "QCD_order = " << QCD_order << endl;
    // determine alpha_S_t used in Gamma_t
    if (valpha_S_t[QCD_order] == empty_double){
      if (valpha_S_t_scale[QCD_order] == empty_double){valpha_S_t_scale[QCD_order] = M_t;}
      if (valpha_S_t_order[QCD_order] == empty_int){valpha_S_t_order[QCD_order] = switch_t;}
      logger << LOG_DEBUG << "LHAPDF:   QCD_order = " << QCD_order << endl;
      logger << LOG_DEBUG << "LHAPDF:   valpha_S_t_order[QCD_order = " << QCD_order << "] = " << valpha_S_t_order[QCD_order] << endl;
      logger << LOG_DEBUG << "LHAPDF:   contribution_LHAPDFsubset[valpha_S_t_order[QCD_order = " << QCD_order << "] = " << valpha_S_t_order[QCD_order] << "] = " << contribution_LHAPDFsubset[valpha_S_t_order[QCD_order]] << endl;
      initialization_LHAPDF(contribution_LHAPDFname[valpha_S_t_order[QCD_order]], contribution_LHAPDFsubset[valpha_S_t_order[QCD_order]]);
      alpha_S_t = LHAPDF::alphasPDF(valpha_S_t_scale[QCD_order]);
      logger << LOG_DEBUG << "LHAPDF:   alpha_S_t = " << setw(25) << setprecision(16) << alpha_S_t << endl;
    }
    else {alpha_S_t = valpha_S_t[QCD_order];}
    // determine partial width at the chosen perturbative order with the chosen alpha_S
    if (switch_t == 0){
      // LO decay of top quark
      if (vGamma_tWb[QCD_order] == empty_double){
	if (Gamma_W == 0.){
	  if (M_b == 0.){
	    // LO calculation of top-quark width with a stable W boson and a massless bottom quark
	    vGamma_tWb[QCD_order] = pow(e, 2) * pow(Cminus_W, 2) * pow(pow(M_t, 2) - pow(M_W, 2), 2) / (32. * pi * pow(M_t, 3)) * (2 + pow(M_t, 2) / pow(M_W, 2));
	  }
	  else {
	    // LO calculation of top-quark width with a stable W boson and a massive bottom quark
	    vGamma_tWb[QCD_order] = alpha_e_Gmu * pow(M_t, 3) / (16 * sin2_w * pow(M_W, 2))
	      * sqrt(lambda(1., pow(M_W / M_t, 2), pow(M_b / M_t, 2)))
	      * (pow(1. - pow(M_b / M_t, 2), 2) +  pow(M_W / M_t, 2) * (1. + pow(M_b / M_t, 2)) - 2 * pow(M_W / M_t, 4));
	  }
	}
	else {
	  // LO calculation of top-quark width with a leptonically decaying W boson (corrected for BR) and a massless/massive bottom quark ->numerical calculation
	  topwidth topwidth(*this);
	  vGamma_tWb[QCD_order] = topwidth.integrate_topwidth(alpha_S_t, vGamma_t_order[QCD_order], 0.00001);
	}
      }
    }
    else if (switch_t == 1){
      // NLO QCD decay of top quark
      if (vGamma_tWb[QCD_order] == empty_double){
	if (Gamma_W == 0.){
	  if (M_b == 0.){
	    // NLO QCD calculation of top-quark width with a stable W boson and a massless bottom quark
	    double r2 = pow(M_W / M_t, 2);
	    double Gamma_t_1 = -inv2pi * C_F  * (2. / 3. * pi2 + 4. * gsl_sf_dilog(r2) - 1.5 - 2. * log(r2 / (1. - r2)) + 2. * log(r2) * log (1. - r2) - 4. / (3. * (1. - r2)) + (22. - 34. * r2) / (9. * pow((1. - r2), 2)) * log(r2) + (3. + 27. * log(1. - r2) - 4. * log(r2)) / (9. * (1. + 2. * r2)));
	    logger << LOG_DEBUG << "Gamma_t_1   = " << setw(23) << setprecision(16) << Gamma_t_1 << endl;
	    
	    double Gamma_t_1b = -inv2pi * C_F  * ((pi2 + 2. * gsl_sf_dilog(r2) - 2. * gsl_sf_dilog(1. - r2)) + (4. * r2 * (1. - r2 - 2. * pow(r2, 2)) * log(r2) + 2. * pow(1. - r2, 2) * (5. + 4. * r2) * log(1. - r2) - (1. - r2) * (5. + 9. * r2 - 6. * pow(r2, 2))) / (2. * pow(1. - r2, 2) * (1. + 2. * r2)))
	      ;
	    logger << LOG_DEBUG << "Gamma_t_1b  = " << setw(23) << setprecision(16) << Gamma_t_1b << endl;
	    
	    vGamma_tWb[QCD_order] = pow(e, 2) * pow(Cminus_W, 2) * pow(pow(M_t, 2) - pow(M_W, 2), 2) / (32. * pi * pow(M_t, 3)) * (2 + pow(M_t, 2) / pow(M_W, 2)) * (1. + alpha_S_t * Gamma_t_1);
	  }
	  else {
	    // NLO QCD calculation of top-quark width with a stable W boson and a massive bottom quark
	    // !!! not yet calculated ???
	  }
	}
	else {
	  // NLO QCD calculation of top-quark width with a leptonically decaying W boson (corrected for BR) and a massless/massive bottom quark ->numerical calculation
	  topwidth topwidth(*this);
	  vGamma_tWb[QCD_order] = topwidth.integrate_topwidth(alpha_S_t, vGamma_t_order[QCD_order], 0.00001);
	}
      } 
      logger << LOG_DEBUG << "vGamma_tWb[QCD_order] = " << vGamma_tWb[QCD_order] << endl;
    }
    else if (switch_t == 2){
      // temporarily !!! should be replaced by NNLO QCD value !!!
      if (vGamma_tWb[QCD_order] == empty_double){
	double r2 = pow(M_W / M_t, 2);
	//      double_complex x = 1.;//epsdilog(r2);gsl_sf_dilog(r2)
	double Gamma_t_1 = -inv2pi * C_F  * (2. / 3. * pi2 + 4. * gsl_sf_dilog(r2) - 1.5 - 2. * log(r2 / (1. - r2)) + 2. * log(r2) * log (1. - r2) - 4. / (3. * (1. - r2)) + (22. - 34. * r2) / (9. * pow((1. - r2), 2)) * log(r2) + (3. + 27. * log(1. - r2) - 4. * log(r2)) / (9. * (1. + 2. * r2)));
	logger << LOG_DEBUG << "Gamma_t_1   = " << setw(23) << setprecision(16) << Gamma_t_1 << endl;
	
	double Gamma_t_1b = -inv2pi * C_F  * ((pi2 + 2. * gsl_sf_dilog(r2) - 2. * gsl_sf_dilog(1. - r2)) + (4. * r2 * (1. - r2 - 2. * pow(r2, 2)) * log(r2) + 2. * pow(1. - r2, 2) * (5. + 4. * r2) * log(1. - r2) - (1. - r2) * (5. + 9. * r2 - 6. * pow(r2, 2))) / (2. * pow(1. - r2, 2) * (1. + 2. * r2)))
	  ;
	logger << LOG_DEBUG << "Gamma_t_1b  = " << setw(23) << setprecision(16) << Gamma_t_1b << endl;
	
	vGamma_tWb[QCD_order] = pow(e, 2) * pow(Cminus_W, 2) * pow(pow(M_t, 2) - pow(M_W, 2), 2) / (32. * pi * pow(M_t, 3)) * (2 + pow(M_t, 2) / pow(M_W, 2)) * (1. + alpha_S_t * Gamma_t_1);} // !!! not yet calculated
      logger << LOG_DEBUG << "vGamma_tWb[QCD_order] = " << vGamma_tWb[QCD_order] << endl;
    }
    vGamma_t[QCD_order] = vGamma_tWb[QCD_order];
  }
  else if (vGamma_t[QCD_order] == empty_double && vGamma_tWb[QCD_order] != empty_double){
    vGamma_t[QCD_order] = vGamma_tWb[QCD_order];
  }
  else if (vGamma_t[QCD_order] != empty_double && vGamma_tWb[QCD_order] == empty_double){
    vGamma_tWb[QCD_order] = vGamma_t[QCD_order];
  }
  else {
    logger << LOG_INFO << "vGamma_t[" << QCD_order << "]   = " << vGamma_t[QCD_order] << endl;
    logger << LOG_INFO << "vGamma_tWb[" << QCD_order << "] = " << vGamma_tWb[QCD_order] << endl;
    logger << LOG_INFO << "vGamma_t[" << QCD_order << "] and vGamma_tWb[" << QCD_order << "] are independently defined." << endl;
  }

  // determine total width and branching ratios from the partial widths
  // name 'vGamma_t' is misleading: it is not always a 1->2p decay !!!
  if (vGamma_t[QCD_order] == empty_double){vGamma_t[QCD_order] = vGamma_tWb[QCD_order];}
  if (vBR_tWb[QCD_order] == empty_double){vBR_tWb[QCD_order] = vGamma_tWb[QCD_order] / vGamma_t[QCD_order];}

  Gamma_t = vGamma_t[QCD_order];
  Gamma_tWb = vGamma_tWb[QCD_order];

  BR_tWb = vBR_tWb[QCD_order];



  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

  //  double Gamma_t = 1.4;
  //  double Gamma_t = 1.4654916649759706; // !!! only LO !!!
  //  double Gamma_t = 1.3101250320515745; // !!! only NLO !!!
  //  double Gamma_t = e2 * pow(Cminus_W, 2) * pow(M2_t - M2_W, 2) / (32. * pi * pow(M_t, 3)) * (2 + M2_t / M2_W);
  //  double BR_Wlv = sqrt2 / (12. * pi) * G_F * pow(M_W, 3) / Gamma_W;
  //  double BR_Zee = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2.) * M2_W * (Cplus_Zee * Cplus_Zee + Cminus_Zee * Cminus_Zee) / Gamma_Z;
  //  double BR_tWb = 1.;

  //  double BR_Wlv = sqrt2 / (12. * pi) * G_F * pow(M_W, 3) / Gamma_W;
  //  double BR_Zee = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2.) * M2_W * (Cplus_Zee * Cplus_Zee + Cminus_Zee * Cminus_Zee) / Gamma_Z;
  //  double BR_tWb = 1.;

  //  double alpha_S_MZ = 0.1176;
  //  double alpha_S_MW = alpha_S_from_Lambda_QCD(N_f, temp_loop_order, temp_Lambda_QCD, M_W);
  //  double alpha_S_MZ = alpha_S_from_Lambda_QCD(N_f, temp_loop_order, temp_Lambda_QCD, M_Z);
  //  double Gamma_W = (3 * 2. + 3) * sqrt2 / (12. * pi) * G_F * pow(M_W, 3) * (1. + (2. * alpha_S_MZ) / (3. * pi)); // !!!
  //  double Gamma_W = 2.099736097449861;
  //  double vGamma_W = (3 * 2. + 3) * sqrt2 / (12. * pi) * G_F * pow(M_W, 3) * (1. + (2. * alpha_S_MZ) / (3. * pi)); // !!!
  //  double Gamma_Z = sqrt2 / (6. * pi) * G_F * M_Z * pow(sin_w, 2) * M2_W
  //    * (  6. * (Cplus_Zuu * Cplus_Zuu + Cminus_Zuu * Cminus_Zuu) * (1. + alpha_S_MZ / pi)
  //	 + 9. * (Cplus_Zdd * Cplus_Zdd + Cminus_Zdd * Cminus_Zdd) * (1. + alpha_S_MZ / pi) 
  //	 + 3. * (Cplus_Zee * Cplus_Zee + Cminus_Zee * Cminus_Zee) 
  //	 + 3. * (Cplus_Znn * Cplus_Znn + Cminus_Znn * Cminus_Znn));
  //  inr perturbative_order;



void model_set::determine_CKM_matrix(){
  Logger logger("model_set::determine_CKM_matrix");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (CKM_matrix == empty_string){CKM_matrix = "individual";}

  if (CKM_matrix == "Cabibbo"){
    if (theta_c == empty_double){theta_c = 0.227;}
    V_du = cos(theta_c);
    V_su = sin(theta_c);
    V_bu = 0.;
    V_dc = -sin(theta_c);
    V_sc = cos(theta_c);
    V_bc = 0.;
    V_dt = 0.;
    V_st = 0.;
    V_bt = 1.;
  }

  else if (CKM_matrix == "trivial"){
    V_du = 1.;
    V_su = 0.;
    V_bu = 0.;
    V_dc = 0.;
    V_sc = 1.;
    V_bc = 0.;
    V_dt = 0.;
    V_st = 0.;
    V_bt = 1.;
  }

  else if (CKM_matrix == "individual"){
    if (V_du == empty_double){V_du = 1.;}
    if (V_su == empty_double){V_su = 0.;}
    if (V_bu == empty_double){V_bu = 0.;}
    if (V_dc == empty_double){V_dc = 0.;}
    if (V_sc == empty_double){V_sc = 1.;}
    if (V_bc == empty_double){V_bc = 0.;}
    if (V_dt == empty_double){V_dt = 0.;}
    if (V_st == empty_double){V_st = 0.;}
    if (V_bt == empty_double){V_bt = 1.;}
  }

  V_ckm = vector<vector<double_complex> > (4, vector<double_complex> (4));
  V_ckm[0][0] = 1.;
  V_ckm[0][1] = 0.;
  V_ckm[0][2] = 0.;
  V_ckm[0][3] = 0.;
  V_ckm[1][0] = 0.;
  V_ckm[1][1] = V_du;
  V_ckm[1][2] = V_su;
  V_ckm[1][3] = V_bu;
  V_ckm[2][0] = 0.;
  V_ckm[2][1] = V_dc;
  V_ckm[2][2] = V_sc;
  V_ckm[2][3] = V_bc;
  V_ckm[3][0] = 0.;
  V_ckm[3][1] = V_dt;
  V_ckm[3][2] = V_st;
  V_ckm[3][3] = V_bt;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void model_set::determine_sintheta_real(vector<string> & det_M_gauge, vector<string> & det_theta_w){
  Logger logger("model_set::determine_couplings_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // Electroweak schemes to determine correlated parameters
  if (det_M_gauge.size() == 2){
    if (det_theta_w.size() > 0){cout << "theta_w is overdetermined. It's calculated via cos_w = M_W / M_Z." << endl;}
    cos_w = M_W / M_Z;
    cos2_w = pow(cos_w, 2);
    sin2_w = 1. - cos2_w;
    sin_w = sqrt(sin2_w);
  }
  else if (det_M_gauge.size() == 1){
    if (sin_w == empty_double && sin2_w == empty_double && cos2_w == empty_double && cos_w == empty_double){
      cout << "theta_w is underdetermined. Use default for ";
      if (det_M_gauge[det_M_gauge.size() - 1] == "M_Z"){cout << "M_W ";}
      else if (det_M_gauge[det_M_gauge.size() - 1] == "M_W"){cout << "M_Z ";}
      cout << "and calculate theta_w via cos_w = M_W / M_Z." << endl;
      cos_w = M_W / M_Z;
      cos2_w = pow(cos_w, 2);
      sin2_w = 1. - cos2_w;
      sin_w = sqrt(sin2_w);
    }
    else {
      if (det_theta_w.size() > 1){
	cout << "theta_w is overdetermined. Use " << det_theta_w[det_theta_w.size() - 1] << " to determine theta_w." << endl;
	if      (det_theta_w[det_theta_w.size() - 1] == "sin_w"){sin2_w = empty_double; cos2_w = empty_double; cos_w = empty_double;}
	else if (det_theta_w[det_theta_w.size() - 1] == "sin2_w"){sin_w = empty_double; cos2_w = empty_double; cos_w = empty_double;}
	else if (det_theta_w[det_theta_w.size() - 1] == "cos2_w"){sin_w = empty_double; sin2_w = empty_double; cos_w = empty_double;}
	else if (det_theta_w[det_theta_w.size() - 1] == "cos_w"){sin_w = empty_double; sin2_w = empty_double; cos2_w = empty_double;}
      }
      if (sin_w != empty_double){
	sin2_w = pow(sin_w, 2);
	cos2_w = 1. - sin2_w;
	cos_w = sqrt(cos2_w);
      }
      else if (sin2_w != empty_double){
	cos2_w = 1. - sin2_w;
	cos_w = sqrt(cos2_w);
	sin_w = sqrt(sin2_w);
      }
      else if (cos2_w != empty_double){
	sin2_w = 1. - cos2_w;
	cos_w = sqrt(cos2_w);
	sin_w = sqrt(sin2_w);
      }
      else if (cos_w != empty_double){
	cos2_w = pow(cos_w, 2);
	sin2_w = 1. - cos2_w;
	sin_w = sqrt(sin2_w);
      }
      if      (det_M_gauge[det_M_gauge.size() - 1] == "M_Z"){M_W = M_Z * cos_w;}
      else if (det_M_gauge[det_M_gauge.size() - 1] == "M_W"){M_Z = M_W / cos_w;}
    }
  }
  else if (det_M_gauge.size() == 0){
    if (sin_w == empty_double && sin2_w == empty_double && cos2_w == empty_double && cos_w == empty_double){
      cout << "theta_w is underdetermined. Use default for M_Z and M_W and calculate theta_w via cos_w = M_W / M_Z." << endl;
      cos_w = M_W / M_Z;
      cos2_w = pow(cos_w, 2);
      sin2_w = 1. - cos2_w;
      sin_w = sqrt(sin2_w);
    }
    else {
      if (det_theta_w.size() > 1){
	cout << "theta_w is overdetermined. Use " << det_theta_w[det_theta_w.size() - 1] << " to determine theta_w." << endl;
	if      (det_theta_w[det_theta_w.size() - 1] == "sin_w"){sin2_w = empty_double; cos2_w = empty_double; cos_w = empty_double;}
	else if (det_theta_w[det_theta_w.size() - 1] == "sin2_w"){sin_w = empty_double; cos2_w = empty_double; cos_w = empty_double;}
	else if (det_theta_w[det_theta_w.size() - 1] == "cos2_w"){sin_w = empty_double; sin2_w = empty_double; cos_w = empty_double;}
	else if (det_theta_w[det_theta_w.size() - 1] == "cos_w"){sin_w = empty_double; sin2_w = empty_double; cos2_w = empty_double;}
      }
      if (sin_w != empty_double){
	sin2_w = pow(sin_w, 2);
	cos2_w = 1. - sin2_w;
	cos_w = sqrt(cos2_w);
      }
      else if (sin2_w != empty_double){
	cos2_w = 1. - sin2_w;
	cos_w = sqrt(cos2_w);
	sin_w = sqrt(sin2_w);
      }
      else if (cos2_w != empty_double){
	sin2_w = 1. - cos2_w;
	cos_w = sqrt(cos2_w);
	sin_w = sqrt(sin2_w);
      }
      else if (cos_w != empty_double){
	cos2_w = pow(cos_w, 2);
	sin2_w = 1. - cos2_w;
	sin_w = sqrt(sin2_w);
      }
      M_W = M_Z * cos_w;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void model_set::determine_alpha_e_Gmu(vector<string> & det_EW_coupling, int switch_alpha_rescaling_exp){
  Logger logger("model_set::determine_alpha_e_Gmu");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static double firstcall_G_F = G_F;
  static double firstcall_alpha_e_Gmu = alpha_e_Gmu;
  
  if (firstcall_G_F != empty_double && firstcall_alpha_e_Gmu != empty_double){
    logger << LOG_ERROR << "EW coupling is overdetermined. Use " << det_EW_coupling[det_EW_coupling.size() - 1] << " to determine EW_coupling." << endl;
  }
  else if (firstcall_G_F == empty_double && firstcall_alpha_e_Gmu == empty_double){
    G_F = 1.16637E-05;
    logger << LOG_ERROR << "EW coupling is underdetermined. Use default for G_F to determine EW coupling: G_F = " << G_F << endl;
  }

  if (firstcall_alpha_e_Gmu == empty_double){
    if (use_cms == 0 || use_cms == 2){alpha_e_Gmu = sqrt2 * pow(M_W, 2) * (1. - pow(M_W, 2) / pow(M_Z, 2)) * G_F / pi;}
    else if (use_cms == 1){alpha_e_Gmu = sqrt2 * abs(cM2_W * (1. - cM2_W / cM2_Z)) * G_F / pi;}
    //    else if (use_cms == 3){alpha_e_Gmu = sqrt2 * real(cM2_W * (1. - cM2_W / cM2_Z)) * G_F / pi;}
    else if (use_cms == 3){alpha_e_Gmu = sqrt2 * pow(real(sqrt(cM2_W * (1. - cM2_W / cM2_Z))), 2) * G_F / pi;}

    //    if (use_cms == 0){alpha_e_Gmu = sqrt2 * pow(sin_w, 2) * pow(M_W, 2) * G_F / pi;}
    //    else if (use_cms == 1){alpha_e_Gmu = abs(sqrt2 * pow(csin_w, 2) * pow(cM_W, 2) * G_F / pi);}
    //    else if (use_cms == 2){alpha_e_Gmu = sqrt2 * pow(sin_w, 2) * pow(M_W, 2) * G_F / pi;}
  }
  else if (firstcall_G_F == empty_double){
    if (use_cms == 0 || use_cms == 2){G_F = alpha_e_Gmu * pi / (sqrt2 * pow(M_W, 2) * (1. - pow(M_W, 2) / pow(M_Z, 2)));}
    else if (use_cms == 1){G_F = alpha_e_Gmu * pi / (sqrt2 * abs(cM2_W * (1. - cM2_W / cM2_Z)));}
  }
  
  // Has nothing to do with CMS: was used for alpha-rescaling test. !!!
  if (switch_alpha_rescaling_exp < 0){
    alpha_e_Gmu = alpha_e_Gmu * pow(10, switch_alpha_rescaling_exp);
    logger << LOG_INFO << "alpha_e_Gmu is rescaled by a factor of " << pow(10, switch_alpha_rescaling_exp) << endl;
    logger << LOG_INFO << "alpha_e_Gmu = " << setprecision(15) << alpha_e_Gmu << endl;
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void model_set::determine_EWcouplings_real(){
  Logger logger("model_set::determine_couplings_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (Q_u == empty_double){Q_u = 2. / 3.;}
  if (Q_d == empty_double){Q_d = -1. / 3.;}
  if (Q_v == empty_double){Q_v = 0.;}
  if (Q_l == empty_double){Q_l = -1.;}
  if (Iw_u == empty_double){Iw_u = 1. / 2.;}
  if (Iw_d == empty_double){Iw_d = -1. / 2.;}
  if (Iw_n == empty_double){Iw_n = 1. / 2.;}
  if (Iw_l == empty_double){Iw_l = -1. / 2.;}

  if (Cplus_Zuu == empty_double){Cplus_Zuu = -Q_u * sin_w / cos_w;}
  if (Cplus_Zdd == empty_double){Cplus_Zdd = -Q_d * sin_w / cos_w;}
  if (Cminus_Zuu == empty_double){Cminus_Zuu = -Q_u * sin_w / cos_w + Iw_u / (sin_w * cos_w);}
  if (Cminus_Zdd == empty_double){Cminus_Zdd = -Q_d * sin_w / cos_w + Iw_d / (sin_w * cos_w);}
  if (Cplus_Znn == empty_double){Cplus_Znn = -Q_v * sin_w / cos_w;}
  if (Cplus_Zee == empty_double){Cplus_Zee = -Q_l * sin_w / cos_w;}
  if (Cminus_Znn == empty_double){Cminus_Znn = -Q_v * sin_w / cos_w + Iw_n / (sin_w * cos_w);}
  if (Cminus_Zee == empty_double){Cminus_Zee = -Q_l * sin_w / cos_w + Iw_l / (sin_w * cos_w);}
  if (C_ZWminusWplus == empty_double){C_ZWminusWplus = cos_w / sin_w;}

  if (Cplus_Auu == empty_double){Cplus_Auu = -Q_u;}
  if (Cplus_Add == empty_double){Cplus_Add = -Q_d;}
  if (Cminus_Auu == empty_double){Cminus_Auu = -Q_u;}
  if (Cminus_Add == empty_double){Cminus_Add = -Q_d;}
  if (Cplus_Ann == empty_double){Cplus_Ann = 0.;}
  if (Cplus_Aee == empty_double){Cplus_Aee = 1.;}
  if (Cminus_Ann == empty_double){Cminus_Ann = 0.;}
  if (Cminus_Aee == empty_double){Cminus_Aee = 1.;}
  if (C_AWminusWplus == empty_double){C_AWminusWplus = -1.;}
  if (Cminus_W == empty_double){Cminus_W = 1. / (sqrt2 * sin_w);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void model_set::determine_EWcouplings_complex(int switch_costheta_real){
  Logger logger("model_set::determine_couplings_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;


  if (switch_costheta_real){
    /*
    M2_W = pow(M_W, 2);
    M2_Z = pow(M_Z, 2);
    if (ccos2_w == empty_double_complex){ccos2_w = M2_W / M2_Z;}
    */
    if (ccos2_w == empty_double_complex){ccos2_w = cos2_w;}
    if (csin2_w == empty_double_complex){csin2_w = 1. - ccos2_w;}
    if (ccos_w == empty_double_complex){ccos_w = sqrt(ccos2_w);}
    if (csin_w == empty_double_complex){csin_w = sqrt(csin2_w);}
  }
  else {
    cM2_W = M2_W - ri * M_W * Gamma_W;
    cM2_Z = M2_Z - ri * M_Z * Gamma_Z;
    if (ccos2_w == empty_double_complex){ccos2_w = cM2_W / cM2_Z;}
    if (csin2_w == empty_double_complex){csin2_w = 1. - ccos2_w;}
    if (ccos_w == empty_double_complex){ccos_w = sqrt(ccos2_w);}
    if (csin_w == empty_double_complex){csin_w = sqrt(csin2_w);}
  }

  if (cCplus_Zuu == empty_double_complex){cCplus_Zuu = -Q_u * csin_w / ccos_w;}
  if (cCplus_Zdd == empty_double_complex){cCplus_Zdd = -Q_d * csin_w / ccos_w;}
  if (cCminus_Zuu == empty_double_complex){cCminus_Zuu = -Q_u * csin_w / ccos_w + Iw_u / (csin_w * ccos_w);}
  if (cCminus_Zdd == empty_double_complex){cCminus_Zdd = -Q_d * csin_w / ccos_w + Iw_d / (csin_w * ccos_w);}
  if (cCplus_Znn == empty_double_complex){cCplus_Znn = -Q_v * csin_w / ccos_w;}
  if (cCplus_Zee == empty_double_complex){cCplus_Zee = -Q_l * csin_w / ccos_w;}
  if (cCminus_Znn == empty_double_complex){cCminus_Znn = -Q_v * csin_w / ccos_w + Iw_n / (csin_w * ccos_w);}
  if (cCminus_Zee == empty_double_complex){cCminus_Zee = -Q_l * csin_w / ccos_w + Iw_l / (csin_w * ccos_w);}
  if (cC_ZWminusWplus == empty_double_complex){cC_ZWminusWplus = ccos_w / csin_w;}
  if (cCminus_W == empty_double_complex){cCminus_W = 1. / (sqrt2 * csin_w);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void model_set::fill_particle_mass_vector(){
  Logger logger("model_set::fill_particle_mass_vector");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG << "   Gamma_W = " << Gamma_W << endl;
  if (reg_Gamma_W == empty_double){reg_Gamma_W = 0.;}
  if (Gamma_W != 0. && reg_Gamma_W == 0.){reg_Gamma_W = Gamma_W;}
  logger << LOG_DEBUG << "reg_Gamma_W = " << reg_Gamma_W << endl;
  if      (Gamma_W == 0. && reg_Gamma_W != 0.){map_Gamma_W = reg_Gamma_W;}
  else if (Gamma_W == 0. && reg_Gamma_W == 0.){map_Gamma_W = 0.;}
  else if (Gamma_W != 0. && reg_Gamma_W != 0.){map_Gamma_W = Gamma_W;}
  //  else if (Gamma_W != 0. && reg_Gamma_W == 0.){map_Gamma_W = Gamma_W;}
  logger << LOG_DEBUG << "map_Gamma_W = " << map_Gamma_W << endl;

  logger << LOG_DEBUG << "   Gamma_Z = " << Gamma_Z << endl;
  if (map_Gamma_Z == empty_double){
    if (reg_Gamma_Z == empty_double){reg_Gamma_Z = 0.;}
    if (Gamma_Z != 0. && reg_Gamma_Z == 0.){reg_Gamma_Z = Gamma_Z;}
    logger << LOG_DEBUG << "reg_Gamma_Z = " << reg_Gamma_Z << endl;
    if      (Gamma_Z == 0. && reg_Gamma_Z != 0.){map_Gamma_Z = reg_Gamma_Z;}
    else if (Gamma_Z == 0. && reg_Gamma_Z == 0.){map_Gamma_Z = 0.;}
    else if (Gamma_Z != 0. && reg_Gamma_Z != 0.){map_Gamma_Z = Gamma_Z;}
    //  else if (Gamma_Z != 0. && reg_Gamma_Z == 0.){map_Gamma_Z = Gamma_Z;}
    logger << LOG_DEBUG << "map_Gamma_Z = " << map_Gamma_Z << endl;
  }
  
  logger << LOG_DEBUG << "   Gamma_H = " << Gamma_H << endl;
  if (reg_Gamma_H == empty_double){reg_Gamma_H = 0.;}
  if (Gamma_H != 0. && reg_Gamma_H == 0.){reg_Gamma_H = Gamma_H;}
  logger << LOG_DEBUG << "reg_Gamma_H = " << reg_Gamma_H << endl;
  if      (Gamma_H == 0. && reg_Gamma_H != 0.){map_Gamma_H = reg_Gamma_H;}
  else if (Gamma_H == 0. && reg_Gamma_H == 0.){map_Gamma_H = 0.;}
  else if (Gamma_H != 0. && reg_Gamma_H != 0.){map_Gamma_H = Gamma_H;}
  //  else if (Gamma_H != 0. && reg_Gamma_H == 0.){map_Gamma_H = Gamma_H;}
  logger << LOG_DEBUG << "map_Gamma_H = " << map_Gamma_H << endl;

  logger << LOG_DEBUG << "   Gamma_t = " << Gamma_t << endl;
  if (reg_Gamma_t == empty_double){reg_Gamma_t = 0.;}
  if (Gamma_t != 0. && reg_Gamma_t == 0.){reg_Gamma_t = Gamma_t;}
  logger << LOG_DEBUG << "reg_Gamma_t = " << reg_Gamma_t << endl;
  if      (Gamma_t == 0. && reg_Gamma_t != 0.){map_Gamma_t = reg_Gamma_t;}
  else if (Gamma_t == 0. && reg_Gamma_t == 0.){map_Gamma_t = 0.;}
  else if (Gamma_t != 0. && reg_Gamma_t != 0.){map_Gamma_t = Gamma_t;}
  //  else if (Gamma_t != 0. && reg_Gamma_t == 0.){map_Gamma_t = Gamma_t;}
  logger << LOG_DEBUG << "map_Gamma_t = " << map_Gamma_t << endl;



  M2_W = pow(M_W, 2);
  M2_Z = pow(M_Z, 2);
  M2_H = pow(M_H, 2);
  M2_t = pow(M_t, 2);
  M2_b = pow(M_b, 2);
  M2_c = pow(M_c, 2);
  M2_s = pow(M_s, 2);
  M2_u = pow(M_u, 2);
  M2_d = pow(M_d, 2);
  M2_e = pow(M_e, 2);
  M2_mu = pow(M_mu, 2);
  M2_tau = pow(M_tau, 2);
  M2_ve = pow(M_ve, 2);
  M2_vm = pow(M_vm, 2);
  M2_vt = pow(M_vt, 2);

  cM2_W = M2_W - ri * M_W * Gamma_W;
  cM2_Z = M2_Z - ri * M_Z * Gamma_Z;
  cM2_H = M2_H - ri * M_H * Gamma_H;
  cM2_t = M2_t - ri * M_t * Gamma_t;

  cM_W = sqrt(cM2_W);
  cM_Z = sqrt(cM2_Z);
  cM_H = sqrt(cM2_H);
  cM_t = sqrt(cM2_t);



  M.resize(26, 0.);
  M[5] = M_b;
  M[6] = M_t;
  M[15] = M_tau;
  M[24] = M_W;
  M[23] = M_Z;
  M[25] = M_H;

  M2.resize(26, 0.);
  M2[5] = M2_b;
  M2[6] = M2_t;
  M2[15] = M2_tau;
  M2[24] = M2_W;
  M2[23] = M2_Z;
  M2[25] = M2_H;

  Gamma.resize(26, 0.);
  //  Gamma[5] = Gamma_b;
  Gamma[6] = Gamma_t;
  Gamma[24] = Gamma_W;
  Gamma[23] = Gamma_Z;
  Gamma[25] = Gamma_H;


  //  vector<double> 
  map_Gamma.resize(26, 0.);
  //  map_Gamma[5] = map_Gamma_b;
  map_Gamma[6] = map_Gamma_t;
  map_Gamma[24] = map_Gamma_W;
  map_Gamma[23] = map_Gamma_Z;
  map_Gamma[25] = map_Gamma_H;

  //  vector<double> 
  reg_Gamma.resize(26, 0.);
  //  reg_Gamma[5] = reg_Gamma_b;
  reg_Gamma[6] = reg_Gamma_t;
  reg_Gamma[24] = reg_Gamma_W;
  reg_Gamma[23] = reg_Gamma_Z;
  reg_Gamma[25] = reg_Gamma_H;

  //  vector<double_complex> 
  cM.resize(26, 0.);
  //  cM[5] = cM_b;
  cM[6] = cM_t;
  cM[24] = cM_W;
  cM[23] = cM_Z;
  cM[25] = cM_H;
  //  vector<double_complex> 
  cM2.resize(26, 0.);
  //  cM2[5] = cM2_b;
  cM2[6] = cM2_t;
  cM2[24] = cM2_W;
  cM2[23] = cM2_Z;
  cM2[25] = cM2_H;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void model_set::get_simple_userinput_from_readin(vector<string> & user_variable, vector<string> & user_value, vector<string> & readin){
  Logger logger("model_set::get_simple_userinput_from_readin");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  string s0 = "";
  string readindata;
  for (int i = 0; i < readin.size(); i++){
    readindata = readin[i][0];
    if (readindata != "/" && readindata != "#" && readindata != "%"){
      int start = 0;
      user_variable.push_back(s0);
      user_value.push_back(s0);
      for (int j = 0; j < readin[i].size(); j++){
	if (start == 0 || start == 1){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 0){}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    user_variable[user_variable.size() - 1].push_back(readin[i][j]);
	    if (start != 1){start = 1;}
	  }
	  else {start++;}
	}
	else if (start == 2){
	  if (readin[i][j] == '='){start++;}
	  else if ((readin[i][j] == ' ') || (readin[i][j] == char(9))){}
	  else {logger << LOG_WARN << "Incorrect input in line " << i  << endl;
	    user_variable.erase(user_variable.end(), user_variable.end());
	    user_value.erase(user_value.end(), user_value.end());
	    break;
	  }
	}
	else if (start == 3 || start == 4){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 3){}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    user_value[user_value.size() - 1].push_back(readin[i][j]);
	    if (start != 4){start = 4;}
	  }
	  else {start++;}
	}
	else {break;}
      }
    }
  }

  for (int i = 0; i < user_variable.size(); i++){
    logger << LOG_DEBUG << "user_variable[" << setw(3) << i << "] = " << setw(30) << user_variable[i] << " = " << setw(30) << user_value[i] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void model_set::initialization_LHAPDF(string & this_lhapdfname, int this_lhapdfsubset){
  static Logger logger("model_set::initialization_LHAPDF");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  const int SUBSET = this_lhapdfsubset;

  cout << "this_lhapdfname.length() = " << this_lhapdfname.length() << endl;
  if (this_lhapdfname.length()>=6) {
    if (this_lhapdfname.substr(this_lhapdfname.size() - 6, 6) == ".LHpdf"){
      const string NAME = this_lhapdfname.substr(0, this_lhapdfname.size() - 6);
      LHAPDF::initPDFSet(NAME, LHAPDF::LHPDF, SUBSET);
      //    LHAPDF::initPDFSet(NAME, LHAPDF::LHPDF, SUBSET);
      cout << "XXX initialization done here!!! this_lhapdfname.length()>=6 .LHpdf" << endl;
      return;
    }
  }
  if (this_lhapdfname.length()>=7) {
    if (this_lhapdfname.substr(this_lhapdfname.size() - 7, 7) == ".LHgrid"){
      const string NAME = this_lhapdfname.substr(0, this_lhapdfname.size() - 7);
      LHAPDF::initPDFSet(NAME, LHAPDF::LHGRID, SUBSET);
      cout << "XXX initialization done here!!! this_lhapdfname.length()>=7 .LHgrid" << endl;
      return;
    }
  }

  const string NAME = this_lhapdfname;
  LHAPDF::initPDFSet(NAME, LHAPDF::LHGRID, SUBSET);
  logger << LOG_DEBUG_VERBOSE << "LHAPDF::initPDFSet(NAME = " << NAME << ", LHAPDF::LHGRID, SUBSET = " << SUBSET << ");" << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void model_set::output_model_file(int switch_costheta_real){
  Logger logger("model_set::fill_particle_mass_vector");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG << "model file" << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "M_W   = " << right << setprecision(16) << setw(23) << M_W << endl;
  logger << LOG_DEBUG << "M_Z   = " << right << setprecision(16) << setw(23) << M_Z << endl;
  logger << LOG_DEBUG << "M_H   = " << right << setprecision(16) << setw(23) << M_H << endl;
  logger << LOG_DEBUG << "M_t   = " << right << setprecision(16) << setw(23) << M_t << endl;
  logger << LOG_DEBUG << "M_b   = " << right << setprecision(16) << setw(23) << M_b << endl;
  logger << LOG_DEBUG << "M_c   = " << right << setprecision(16) << setw(23) << M_c << endl;
  logger << LOG_DEBUG << "M_s   = " << right << setprecision(16) << setw(23) << M_s << endl;
  logger << LOG_DEBUG << "M_u   = " << right << setprecision(16) << setw(23) << M_u << endl;
  logger << LOG_DEBUG << "M_d   = " << right << setprecision(16) << setw(23) << M_d << endl;
  logger << LOG_DEBUG << "M_e   = " << right << setprecision(16) << setw(23) << M_e << endl;
  logger << LOG_DEBUG << "M_mu  = " << right << setprecision(16) << setw(23) << M_mu << endl;
  logger << LOG_DEBUG << "M_tau = " << right << setprecision(16) << setw(23) << M_tau << endl;
  logger << LOG_DEBUG << "M_ve  = " << right << setprecision(16) << setw(23) << M_ve << endl;
  logger << LOG_DEBUG << "M_vm  = " << right << setprecision(16) << setw(23) << M_vm << endl;
  logger << LOG_DEBUG << "M_vt  = " << right << setprecision(16) << setw(23) << M_vt << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "Gamma_W   = " << setprecision(16) << setw(23) << Gamma_W << endl;
  logger << LOG_DEBUG << "Gamma_Wlv = " << setprecision(16) << setw(23) << Gamma_Wlv << endl;
  logger << LOG_DEBUG << "Gamma_Wud = " << setprecision(16) << setw(23) << Gamma_Wud << endl;
  logger << LOG_DEBUG << "BR_Wlv    = " << setprecision(16) << setw(23) << BR_Wlv << endl;
  logger << LOG_DEBUG << "BR_Wud    = " << setprecision(16) << setw(23) << BR_Wud << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "Gamma_Z   = " << setprecision(16) << setw(23) << Gamma_Z << endl;
  logger << LOG_DEBUG << "Gamma_Zll = " << setprecision(16) << setw(23) << Gamma_Zll << endl;
  logger << LOG_DEBUG << "Gamma_Zvv = " << setprecision(16) << setw(23) << Gamma_Zvv << endl;
  logger << LOG_DEBUG << "Gamma_Zdd = " << setprecision(16) << setw(23) << Gamma_Zdd << endl;
  logger << LOG_DEBUG << "Gamma_Zuu = " << setprecision(16) << setw(23) << Gamma_Zuu << endl;
  logger << LOG_DEBUG << "BR_Zll    = " << setprecision(16) << setw(23) << BR_Zll << endl;
  logger << LOG_DEBUG << "BR_Zvv    = " << setprecision(16) << setw(23) << BR_Zvv << endl;
  logger << LOG_DEBUG << "BR_Zdd    = " << setprecision(16) << setw(23) << BR_Zdd << endl;
  logger << LOG_DEBUG << "BR_Zuu    = " << setprecision(16) << setw(23) << BR_Zuu << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "Gamma_H   = " << setprecision(16) << setw(23) << Gamma_H << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "Gamma_t   = " << setprecision(16) << setw(23) << Gamma_t << endl;
  logger << LOG_DEBUG << "Gamma_tWb = " << setprecision(16) << setw(23) << Gamma_tWb << endl;
  logger << LOG_DEBUG << "BR_tWb    = " << setprecision(16) << setw(23) << BR_tWb << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "cM2_W   = " << right << setprecision(16) << setw(23) << cM2_W << endl;
  logger << LOG_DEBUG << "cM2_Z   = " << right << setprecision(16) << setw(23) << cM2_Z << endl;
  logger << LOG_DEBUG << "cM2_H   = " << right << setprecision(16) << setw(23) << cM2_H << endl;
  logger << LOG_DEBUG << "cM2_t   = " << right << setprecision(16) << setw(23) << cM2_t << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "cM_W   = " << right << setprecision(16) << setw(23) << cM_W << endl;
  logger << LOG_DEBUG << "cM_Z   = " << right << setprecision(16) << setw(23) << cM_Z << endl;
  logger << LOG_DEBUG << "cM_H   = " << right << setprecision(16) << setw(23) << cM_H << endl;
  logger << LOG_DEBUG << "cM_t   = " << right << setprecision(16) << setw(23) << cM_t << endl;
  logger << LOG_DEBUG << endl;
  stringstream temp_ss_u;
  temp_ss_u << setw(8) << "" << "( V_du  V_su  V_bu )   (";
  for (int i = 1; i < 4; i++){temp_ss_u << setprecision(16) << setw(25) << V_ckm[1][i];}
  temp_ss_u << " )";
  logger << LOG_DEBUG << temp_ss_u.str() << endl;
  stringstream temp_ss_c;
  temp_ss_c << setw(8) << "V_ckm = " << "( V_dc  V_sc  V_bc ) = (";
  for (int i = 1; i < 4; i++){temp_ss_c << setprecision(16) << setw(25) << V_ckm[2][i];}
  temp_ss_c << " )";
  logger << LOG_DEBUG << temp_ss_c.str() << endl;
  stringstream temp_ss_t;
  temp_ss_t << setw(8) << "" << "( V_dt  V_st  V_bt )   (";
  for (int i = 1; i < 4; i++){temp_ss_t << setprecision(16) << setw(25) << V_ckm[3][i];}
  temp_ss_t << " )";
  logger << LOG_DEBUG << temp_ss_t.str() << endl;
  logger << LOG_DEBUG << endl;
  stringstream name_ew_scheme;
  if (ew_scheme == 0){name_ew_scheme << "(alpha_qed_0)";}
  else if (ew_scheme == 1){name_ew_scheme << "(Gmu)";}
  else if (ew_scheme == 2){name_ew_scheme << "(alpha_qed_mz)";}
  logger << LOG_DEBUG << "ew_scheme   = " << ew_scheme << " " << name_ew_scheme.str() << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "alpha_e_0   = " << setprecision(16) << setw(23) << alpha_e_0 << endl;
  logger << LOG_DEBUG << "alpha_e_Gmu = " << setprecision(16) << setw(23) << alpha_e_Gmu << endl;
  logger << LOG_DEBUG << "alpha_e_MZ  = " << setprecision(16) << setw(23) << alpha_e_MZ << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "alpha_e     = " << setprecision(16) << setw(23) << alpha_e << " " << name_ew_scheme.str() << endl;
  logger << LOG_DEBUG << "e       = " << setprecision(16) << setw(23) << e << endl;
  if (ew_scheme == 1){
    logger << LOG_DEBUG << "G_F     = " << setprecision(16) << setw(23) << G_F << endl;
  }
  logger << LOG_DEBUG << endl;
  if (switch_costheta_real == 0){
    logger << LOG_DEBUG << "cos_w   = " << setprecision(16) << setw(23) << cos_w << endl;
    logger << LOG_DEBUG << "cos2_w  = " << setprecision(16) << setw(23) << cos2_w << endl;
    logger << LOG_DEBUG << "sin2_w  = " << setprecision(16) << setw(23) << sin2_w << endl;
    logger << LOG_DEBUG << "sin_w   = " << setprecision(16) << setw(23) << sin_w << endl;
  }
  else if (switch_costheta_real == 1){
    logger << LOG_DEBUG << "ccos_w   = " << setprecision(16) << setw(23) << ccos_w << endl;
    logger << LOG_DEBUG << "ccos2_w  = " << setprecision(16) << setw(23) << ccos2_w << endl;
    logger << LOG_DEBUG << "csin2_w  = " << setprecision(16) << setw(23) << csin2_w << endl;
    logger << LOG_DEBUG << "csin_w   = " << setprecision(16) << setw(23) << csin_w << endl;
  }
  else {
    logger << LOG_ERROR << "Invalid value:  switch_costheta_real = " << switch_costheta_real << endl;
    exit(1);
  }
  //  logger << LOG_DEBUG << "alpha_e = " << setprecision(16) << setw(23) << alpha_e << endl;
  logger << LOG_DEBUG << endl;
  if (switch_costheta_real == 0){
    logger << LOG_DEBUG << "Cplus_Zuu = " << Cplus_Zuu << endl;
    logger << LOG_DEBUG << "Cplus_Zdd = " << Cplus_Zdd << endl;
    logger << LOG_DEBUG << "Cminus_Zuu = " << Cminus_Zuu << endl;
    logger << LOG_DEBUG << "Cminus_Zdd = " << Cminus_Zdd << endl;
    logger << LOG_DEBUG << "Cplus_Znn = " << Cplus_Znn << endl;
    logger << LOG_DEBUG << "Cplus_Zee = " << Cplus_Zee << endl;
    logger << LOG_DEBUG << "Cminus_Znn = " << Cminus_Znn << endl;
    logger << LOG_DEBUG << "Cminus_Zee = " << Cminus_Zee << endl;
    logger << LOG_DEBUG << "C_ZWminusWplus = " << C_ZWminusWplus << endl;
  }
  else if (switch_costheta_real == 1){
    logger << LOG_DEBUG << "cCplus_Zuu = " << cCplus_Zuu << endl;
    logger << LOG_DEBUG << "cCplus_Zdd = " << cCplus_Zdd << endl;
    logger << LOG_DEBUG << "cCminus_Zuu = " << cCminus_Zuu << endl;
    logger << LOG_DEBUG << "cCminus_Zdd = " << cCminus_Zdd << endl;
    logger << LOG_DEBUG << "cCplus_Znn = " << cCplus_Znn << endl;
    logger << LOG_DEBUG << "cCplus_Zee = " << cCplus_Zee << endl;
    logger << LOG_DEBUG << "cCminus_Znn = " << cCminus_Znn << endl;
    logger << LOG_DEBUG << "cCminus_Zee = " << cCminus_Zee << endl;
    logger << LOG_DEBUG << "cC_ZWminusWplus = " << cC_ZWminusWplus << endl;
    logger << LOG_DEBUG << "cCminus_W = " << cCminus_W << endl;
  }
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "Cplus_Auu = " << Cplus_Auu << endl;
  logger << LOG_DEBUG << "Cplus_Add = " << Cplus_Add << endl;
  logger << LOG_DEBUG << "Cminus_Auu = " << Cminus_Auu << endl;
  logger << LOG_DEBUG << "Cminus_Add = " << Cminus_Add << endl;
  logger << LOG_DEBUG << "Cplus_Ann = " << Cplus_Ann << endl;
  logger << LOG_DEBUG << "Cplus_Aee = " << Cplus_Aee << endl;
  logger << LOG_DEBUG << "Cminus_Ann = " << Cminus_Ann << endl;
  logger << LOG_DEBUG << "Cminus_Aee = " << Cminus_Aee << endl;
  logger << LOG_DEBUG << "C_AWminusWplus = " << C_AWminusWplus << endl;
  logger << LOG_DEBUG << "Cminus_W = " << Cminus_W << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "map_Gamma_W   = " << setprecision(16) << setw(23) << map_Gamma_W << endl;
  logger << LOG_DEBUG << "map_Gamma_Z   = " << setprecision(16) << setw(23) << map_Gamma_Z << endl;
  logger << LOG_DEBUG << "map_Gamma_H   = " << setprecision(16) << setw(23) << map_Gamma_H << endl;
  logger << LOG_DEBUG << "map_Gamma_t   = " << setprecision(16) << setw(23) << map_Gamma_t << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG << "reg_Gamma_W   = " << setprecision(16) << setw(23) << reg_Gamma_W << endl;
  logger << LOG_DEBUG << "reg_Gamma_Z   = " << setprecision(16) << setw(23) << reg_Gamma_Z << endl;
  logger << LOG_DEBUG << "reg_Gamma_H   = " << setprecision(16) << setw(23) << reg_Gamma_H << endl;
  logger << LOG_DEBUG << "reg_Gamma_t   = " << setprecision(16) << setw(23) << reg_Gamma_t << endl;
  logger << LOG_DEBUG << endl;
  logger << LOG_DEBUG_VERBOSE << "msi.BR_Wlv = " << BR_Wlv << endl;
  logger << LOG_DEBUG << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
