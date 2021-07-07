#include "../include/classes.cxx"
////////////////////
//  constructors  //
////////////////////
munich::munich(){}

munich::munich(int argc, char *argv[], string basic_process_class){
  static Logger logger("munich::munich");
  logger << LOG_INFO << "called" << endl;

  Log::setLogThreshold(LOG_WARN);

  logger << LOG_DEBUG << "argc = " << argc << endl;
  if (argc == 0){logger << LOG_ERROR << "No input specified." << endl; exit(1);}
  else if (argc < 3){
    if (argc == 1){
      cin >> subprocess;
      logger << LOG_DEBUG << "console input = " << subprocess << endl;
    }
    else if (argc == 2){
      logger << LOG_DEBUG << "argv[1] = " << argv[1] << endl;
      subprocess = argv[1];
    }

    //    vector<string> readin;
    //    parameter_readin(0, subprocess, readin);
    //


    logger << LOG_DEBUG << "basic_process_class = " << basic_process_class << endl;
    isi = inputparameter_set(basic_process_class, subprocess);//, readin);
    
    csi = isi.csi;
    logger << LOG_DEBUG << "csi.basic_process_class = " << csi.basic_process_class << endl;
    logger << LOG_DEBUG << "csi.process_class = " << csi.process_class << endl;
    csi.determination_order_alpha_s_born();
    csi.determination_class_contribution();

    
    esi = isi.esi;
    user = isi.user;
    
    osi = observable_set(isi, csi);
    msi = osi.msi;
    
    
    psi = phasespace_set(isi, csi);
    psi.generic = &generic;

    osi.psi = &psi;
    osi.xmunich = this;
    
    logger.newLine(LOG_INFO);
    logger << LOG_INFO << "Settings from input files:" << endl << endl << isi << endl;
    //    cout << isi << endl;
  
  }
  else if (argc > 2){
    //    vector<string> readin;
    //    parameter_readin(readin);
    isi = inputparameter_set(basic_process_class, subprocess);//, readin);

    order = argv[1];
    if (argc == 2){
      cout << "argv[1] = " << argv[1] << endl;
      if (order == "result"){cout << "calculate all results" << endl;}
      else if (order == "distribution"){cout << "calculate all distributions" << endl;}
      infilename = "infile.result.all";
    }
    else if (argc > 2){
      cout << "argv[1] = " << argv[1] << endl;
      cout << "argv[2] = " << argv[2] << endl;
      if (order == "result"){cout << "calculate the results specified in " << argv[2] << endl;}
      else if (order == "distribution"){cout << "calculate the results specified in " << argv[2] << endl;}
      else if (order == "scaleband"){cout << "calculate the results specified in " << argv[2] << endl;}
      infilename = argv[2];
    }
    //    string 
    infilename_scaleband = "";
    if (argc == 4){
      cout << "argv[3] = " << argv[3] << endl;
      infilename_scaleband = argv[3];
      cout << "infilename_scaleband = " << infilename_scaleband << endl;
    }
    
    osi = observable_set(isi, csi);
    //    ysi = summary_generic(*this);
  }

  logger << LOG_INFO << "finished" << endl;
}

void munich::run_initialization(){
  static Logger logger("munich::run_initialization");
  logger << LOG_INFO << "called" << endl;

  // Should be done via msi in the future !!!
  //  psi.initialization_masses(osi.M, osi.M2, osi.msi.Gamma, osi.msi.cM2, osi.msi.map_Gamma, osi.msi.reg_Gamma);
  psi.initialization_masses(msi);

  psi.initialization_contribution_order(csi);


  logger << LOG_INFO << "finished" << endl;
}

void munich::run_integration(){
  static Logger logger("munich::run_integration");
  logger << LOG_INFO << "called" << endl;
  logger << LOG_INFO << "csi.type_contribution = " << csi.type_contribution << endl;
  logger << LOG_INFO << "csi.type_correction = " << csi.type_correction << endl;

  if (csi.type_contribution == "born"){integration_born();}
  if (csi.type_contribution == "loop"){integration_born();}
  if (csi.type_contribution == "L2I"){integration_born();}

  if (csi.type_contribution == "VA" && csi.type_correction == "QCD"){integration_VA_QCD();}
  if (csi.type_contribution == "CA" && csi.type_correction == "QCD"){integration_CA_QCD();}
  if (csi.type_contribution == "RA" && csi.type_correction == "QCD"){integration_RA_QCD();}
  
  if (csi.type_contribution == "L2VA" && csi.type_correction == "QCD"){integration_VA_QCD();}
  if (csi.type_contribution == "L2CA" && csi.type_correction == "QCD"){integration_CA_QCD();}
  if (csi.type_contribution == "L2RA" && csi.type_correction == "QCD"){integration_RA_QCD();}
  
  if (csi.type_contribution == "VA" && csi.type_correction == "QEW"){integration_VA_QEW();}
  if (csi.type_contribution == "CA" && csi.type_correction == "QEW"){integration_CA_QEW();}
  if (csi.type_contribution == "RA" && csi.type_correction == "QEW"){integration_RA_QEW();}
  if (csi.type_contribution == "VA" && csi.type_correction == "MIX"){integration_VA_MIX();}
  if (csi.type_contribution == "RA" && csi.type_correction == "MIX"){integration_RA_MIX();}

  if (csi.type_contribution == "VT" && csi.type_correction == "QCD"){integration_VT_QCD();}
  if (csi.type_contribution == "CT" && csi.type_correction == "QCD"){integration_CT_QCD();}
  if (csi.type_contribution == "RT" && csi.type_correction == "QCD"){integration_born();}

  if (csi.type_contribution == "L2VT" && csi.type_correction == "QCD"){integration_VT_QCD();}
  if (csi.type_contribution == "L2CT" && csi.type_correction == "QCD"){integration_CT_QCD();}
  if (csi.type_contribution == "L2RT" && csi.type_correction == "QCD"){integration_born();}

  if (csi.type_contribution == "VT2" && csi.type_correction == "QCD"){integration_VT2_QCD();}
  if (csi.type_contribution == "CT2" && csi.type_correction == "QCD"){integration_CT2_QCD();}
  if (csi.type_contribution == "RVA" && csi.type_correction == "QCD"){integration_VA_QCD();}
  if (csi.type_contribution == "RCA" && csi.type_correction == "QCD"){integration_CA_QCD();}
  if (csi.type_contribution == "RRA" && csi.type_correction == "QCD"){integration_RA_QCD();}


  if (csi.type_contribution == "NLL_LO" && csi.type_correction == "QCD"){integration_VT_QCD();}
  if (csi.type_contribution == "NLL_NLO" && csi.type_correction == "QCD"){integration_VT_QCD();}
  if (csi.type_contribution == "NNLL_LO" && csi.type_correction == "QCD"){integration_VT_QCD();}
  if (csi.type_contribution == "NNLL_NLO" && csi.type_correction == "QCD"){integration_VT_QCD();}
  if (csi.type_contribution == "NNLL_NNLO" && csi.type_correction == "QCD"){integration_VT2_QCD();}

  
  if (csi.type_contribution == "VJ" && csi.type_correction == "QCD"){integration_VJ_QCD();}
  if (csi.type_contribution == "CJ" && csi.type_correction == "QCD"){integration_CJ_QCD();}
  if (csi.type_contribution == "RJ" && csi.type_correction == "QCD"){integration_born();}

  //  if (csi.type_contribution == "VJ2" && csi.type_correction == "QCD"){integration_VJ2_QCD();}
  //  if (csi.type_contribution == "CJ2" && csi.type_correction == "QCD"){integration_CJ2_QCD();}
  if (csi.type_contribution == "RVJ" && csi.type_correction == "QCD"){integration_VA_QCD();}
  if (csi.type_contribution == "RCJ" && csi.type_correction == "QCD"){integration_CA_QCD();}
  if (csi.type_contribution == "RRJ" && csi.type_correction == "QCD"){integration_RA_QCD();}

  
  logger << LOG_INFO << "finished" << endl;
}

void munich::get_summary(){
  static Logger logger("munich::get_summary");
  logger << LOG_INFO << "called" << endl;

  ysi = summary_generic(*this);

  logger << LOG_INFO << "ysi.order = " << ysi.order << endl;
  logger << LOG_INFO << "munich: order = " << order << endl;
  logger << LOG_INFO << "ysi.infilename = " << ysi.infilename << endl;
  logger << LOG_INFO << "munich: infilename = " << infilename << endl;
  logger << LOG_INFO << "ysi.infilename_scaleband = " << ysi.infilename_scaleband << endl;
  logger << LOG_INFO << "munich: infilename_scaleband = " << infilename_scaleband << endl;

  ysi.get_summary();
  
  logger << LOG_INFO << "finished" << endl;
}



void munich::calculate_p_parton(){
  Logger logger("observable_set::calculate_p_parton");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  osi.p_parton = osi.start_p_parton;
  for (int xi = 0; xi <= csi.n_particle + 2; xi++){
    osi.p_parton[0][psi.o_map[0][xi]] = psi.xbp_all[0][intpow(2, xi - 1)];
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
