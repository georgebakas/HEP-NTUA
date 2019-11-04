#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"
#include "../include/definitions.phasespace.set.cxx"

// different type_corrections could maybe be merged !!!

void calculate_ME2_CA_QCD(observable_set & oset){
  static Logger logger("calculate_ME2_CA_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

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
  static double b_ME2;
  static int n_cc = (osi_n_particle + 2) * (osi_n_particle + 1) / 2;
  static double ewcc;
  static double *M2cc;
  M2cc = new double[n_cc];
  //  ol_evaluate_tree(1, P, &b_ME2);
  if ((osi_type_contribution == "CA" ||
       osi_type_contribution == "RCA") && osi_user_string_value[osi_user_string_map["model"]] != "Bornloop"){
    ol_evaluate_cc(1, P, &b_ME2, M2cc, &ewcc);
    // What is this part made for ???
    /*
    if (osi_process_id == 2){
      logger << LOG_DEBUG_VERBOSE << "b_ME2 (from cc)   = " << b_ME2 << endl;
      ol_evaluate_tree(2, P, &b_ME2);
      logger << LOG_DEBUG_VERBOSE << "b_ME2 (from tree) = " << b_ME2 << endl;
    }
    */
    // ???
  }
  else if (osi_type_contribution == "L2CA" || 
	   osi_user_string_value[osi_user_string_map["model"]] == "Bornloop"){
    ol_evaluate_cc2(1, P, &b_ME2, M2cc, &ewcc);
  }

  for (int i_a = 0; i_a < osi_CA_collinear.size(); i_a++){
    for (int j_a = 0; j_a < osi_CA_collinear[i_a].size(); j_a++){
      if (j_a == 0){osi_CA_ME2_cf[i_a][j_a] = b_ME2;}
      else {osi_CA_ME2_cf[i_a][j_a] = M2cc[osi_CA_collinear[i_a][j_a].no_BLHA_entry];}
      logger << LOG_DEBUG_POINT << "OpenLoops:  CA_ME2[" << i_a << "][" << j_a << "] = " << osi_CA_ME2_cf[i_a][j_a] << endl;
    }
  }
  delete [] M2cc;
  delete [] P;

  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void calculate_ME2_CA_QEW(observable_set & oset){
  static Logger logger("calculate_ME2_CA_QEW");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

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
  static double b_ME2;
  ol_evaluate_tree(1, P, &b_ME2);
  for (int i_a = 0; i_a < osi_CA_collinear.size(); i_a++){
    for (int j_a = 0; j_a < osi_CA_collinear[i_a].size(); j_a++){
      // better shift charge_factor to pdfs !!! (could be wrong in sums over all (anti-)quarks otherwise) !!!
      osi_CA_ME2_cf[i_a][j_a] = osi_CA_collinear[i_a][j_a].charge_factor() * b_ME2;
      logger << LOG_DEBUG_POINT << "OpenLoops:  CA_ME2[" << i_a << "][" << j_a << "] = " << osi_CA_ME2_cf[i_a][j_a] << endl;
    }
  }
  delete [] P;

  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void calculate_ME2check_CA_QCD(observable_set & oset, phasespace_set & psi){
  static Logger logger("calculate_ME2check_CA_QCD");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi_p_parton[0][0].x0() != 0.){
    calculate_ME2_CA_QCD(oset);

    if (osi_massive_QCD){oset.calculate_collinear_QCD_CDST();}
    else {oset.calculate_collinear_QCD_CS();}

    oset.output_testpoint_CA(psi);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void calculate_ME2check_CA_QEW(observable_set & oset, phasespace_set & psi){
  static Logger logger("calculate_ME2check_CA_QEW");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi_p_parton[0][0].x0() != 0.){
    calculate_ME2_CA_QEW(oset);

    if (osi_massive_QEW){oset.calculate_collinear_QEW_CDST();}
    else {oset.calculate_collinear_QEW_CS();}

    oset.output_testpoint_CA(psi);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

