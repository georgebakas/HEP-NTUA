#include "../include/classes.cxx"
//#include "../include/definitions.phasespace.set.cxx"
//#include "../include/definitions.observable.set.cxx"

/*
void munich::handling_cut_psp(){
  static Logger logger("munich::handling_cut_psp");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (psi.MC_tau.end_optimization == 0){
    psi.MC_tau.n_rej++;
    psi.MC_tau.n_rej_channel[psi.MC_tau.channel]++;
  }

  // count cutted events for weight optimization of dipole mappings
  if (psi.RA_x_a != 0) {
    if ((psi.switch_MC_x_dipole != -1) &&(psi.type_contribution == "RA" || 
					  psi.type_contribution == "RRA")){
      if (psi.MC_x_dipole[psi.RA_x_a].end_optimization == 0){
        psi.MC_x_dipole[psi.RA_x_a].n_rej++;
        psi.MC_x_dipole[psi.RA_x_a].n_rej_channel[psi.MC_x_dipole[psi.RA_x_a].channel]++;
      }
    }
  }
  if (psi.tau_opt_end == 0){psi.tau_cuts_channel[psi.tau_channel]++;}
  if (psi.x1x2_opt_end == 0){psi.x1x2_cuts_channel[psi.x1x2_channel]++;}
  psi.i_rej++;
  psi.cuts_channel[psi.MC_phasespace.channel]++;
  psi.random_manager.increase_cut_counter();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/

void munich::handling_vanishing_me2(){
  static Logger logger("munich::handling_vanishing_me2");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG_VERBOSE << "psi.i_tec = " << psi.i_tec << endl;
  logger << LOG_DEBUG_VERBOSE << "psi.i_acc = " << psi.i_acc << endl;
  logger << LOG_DEBUG_VERBOSE << "psi.i_gen = " << psi.i_gen << endl;
  logger << LOG_DEBUG_VERBOSE << "osi.n_event_vanishing_ME2 = " << osi.n_event_vanishing_ME2 << endl;
  logger << LOG_DEBUG_VERBOSE << "osi.check_vanishing_ME2_end = " << osi.check_vanishing_ME2_end << endl;
  logger << LOG_DEBUG_VERBOSE << "osi.flag_vanishing_ME2 = " << osi.flag_vanishing_ME2 << endl;
    
  if (osi.flag_vanishing_ME2 == 1){
    if (osi.switch_console_output_ME2_issue){
      logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej <<  "   " << setw(12) << psi.i_tec << "   ME2 is numerically zero!" << endl;
    }
    psi.i_tec++;
    psi.i_acc++;
  }
  if (psi.i_acc >= osi.n_event_vanishing_ME2){
    if (psi.i_acc == psi.i_tec){
      osi.int_end = 1;
      osi.output_zero_contribution_complete(psi);
    }
    else {
      osi.check_vanishing_ME2_end = 1; 
      if (osi.switch_console_output_ME2_issue){
	logger << LOG_DEBUG << "Matrix element does not vanish!" << endl;
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void munich::errorhandling_c_psp(){
  static Logger logger("munich::errorhandling_c_psp");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi.switch_console_output_phasespace_issue){
    logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "     p != p       @ MC_phasespace.channel  = " << psi.MC_phasespace.channel << endl;
    for (int i_a = 0; i_a < osi.p_parton.size(); i_a++){
      for (int i_p = 0; i_p < osi.p_parton[i_a].size(); i_p++){
	logger << LOG_DEBUG << "osi.p_parton[" << setw(3) << i_a << "][" << setw(3) << i_p << "] = " << osi.p_parton[i_a][i_p] << endl;
      }
    }
  }

  psi.i_nan++;
  psi.random_manager.nancut();
  psi.i_gen--;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void munich::errorhandling_me2(){
  static Logger logger("munich::errorhandling_me2");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi.switch_console_output_ME2_issue){
    logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "   ME2 != ME2     @ channel  = " << psi.MC_phasespace.channel << endl;
    for (int i_a = 0; i_a < osi.p_parton.size(); i_a++){
      for (int i_p = 0; i_p < osi.p_parton[i_a].size(); i_p++){
	logger << LOG_DEBUG << "osi.p_parton[" << setw(3) << i_p << "] = " << osi.p_parton[i_a][i_p] << endl;
      }
    }
  }

  psi.i_nan++;
  psi.random_manager.nancut();
  psi.i_gen--;
}

void munich::errorhandling_alpha_S(){
  static Logger logger("munich::errorhandling_alpha_S");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi.switch_console_output_ME2_issue){
    logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "   alpha_S != alpha_S     @ channel  = " << psi.MC_phasespace.channel << endl;
    //  logger << LOG_INFO << "osi.var_rel_alpha_S != osi.var_rel_alpha_S" << endl;   for (int i_a = 0; i_a < osi.p_parton.size(); i_a++){
    for (int i_a = 0; i_a < osi.p_parton.size(); i_a++){
      for (int i_p = 0; i_p < osi.p_parton[i_a].size(); i_p++){
	logger << LOG_DEBUG << "osi.p_parton[" << setw(3) << i_p << "] = " << osi.p_parton[i_a][i_p] << endl;
      }
    }
  }

  psi.i_nan++;
  psi.random_manager.nancut();
  psi.i_gen--;
}

void munich::errorhandling_pdf(){
  static Logger logger("munich::errorhandling_pdf");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi.switch_console_output_ME2_issue){
    logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "   pdf_factor != pdf_factor     @ channel  = " << psi.MC_phasespace.channel << endl;
    //  logger << LOG_INFO << "osi.pdf_factor[0] != osi.pdf_factor[0]" << endl;
    for (int i_p = 0; i_p < osi.x_pdf.size(); i_p++){
      logger << LOG_DEBUG << "osi.x_pdf[" << setw(3) << i_p << "] = " << osi.x_pdf[i_p] << endl;
    }
    for (int i_c = 0; i_c < osi.pdf_factor.size(); i_c++){
      logger << LOG_DEBUG << "osi.pdf_factor[" << setw(3) << i_c << "] = " << osi.pdf_factor[i_c] << endl;
    }
    for (int i_a = 0; i_a < osi.p_parton.size(); i_a++){
      for (int i_p = 0; i_p < osi.p_parton[i_a].size(); i_p++){
	logger << LOG_DEBUG << "osi.p_parton[" << setw(3) << i_p << "] = " << osi.p_parton[i_a][i_p] << endl;
      }
    }
  }

  psi.i_nan++;
  psi.random_manager.nancut();
  psi.i_gen--;
}

void munich::errorhandling_gtot(){
  static Logger logger("munich::errorhandling_gtot");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi.switch_console_output_phasespace_issue){
    logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "  gtot != gtot    @ channel  = " << psi.MC_phasespace.channel << endl;
    for (int i_p = 0; i_p < osi.p_parton[0].size(); i_p++){
      logger << LOG_DEBUG << "osi.p_parton[" << setw(3) << i_p << "] = " << osi.p_parton[0][i_p] << endl;
    }
    logger << LOG_DEBUG << "g_MC = " << psi.g_MC << endl;
    logger << LOG_DEBUG << "g_pdf = " << psi.g_pdf << endl;
    logger << LOG_DEBUG << "psi.MC_g_IS_global = " << psi.MC_g_IS_global << endl;
    if (psi.g_MC != psi.g_MC){
      for (int i_c = 0; i_c < psi.MC_n_channel; i_c++){
	if (psi.MC_phasespace.g_channel[i_c] != psi.MC_phasespace.g_channel[i_c]){
	  logger << LOG_DEBUG << "psi.MC_phasespace.g_channel[" << setw(4) << i_c << "] = " << psi.MC_phasespace.g_channel[i_c] << endl;
	}
      }
    }
  }

  if (psi.g_MC != psi.g_MC){
    for (int i_c = 0; i_c < psi.MC_n_channel; i_c++){
      if (psi.MC_phasespace.g_channel[i_c] != psi.MC_phasespace.g_channel[i_c]){
	psi.MC_phasespace.g_channel[i_c] = 0.;
      }
    }
  }

  psi.i_nan++;
  psi.random_manager.nancut();
  psi.i_gen--;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void munich::errorhandling_OL(){
  static Logger logger("munich::errorhandling_OL");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (csi.type_contribution == "RVA" ||
      csi.type_contribution == "L2RT" ||
      csi.type_contribution == "L2RJ"){
    osi.counter_killed_qTcut[osi.cut_ps[0]]++;
  }
  
  if (osi.switch_console_output_ME2_issue){
    if (csi.type_contribution == "RVA" ||
	csi.type_contribution == "L2RT" ||
	csi.type_contribution == "L2RJ"){
      logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "   killed by OpenLoops   cut_qT = " << osi.cut_ps[0] << endl;
      //      osi.counter_killed_qTcut[osi.cut_ps[0]]++;
    }
    else {
      logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "   killed by OpenLoops" << endl;
    }
    //  logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   killed by OpenLoops" << endl;
    for (int i_a = 0; i_a < osi.p_parton.size(); i_a++){
      for (int i_p = 0; i_p < osi.p_parton[i_a].size(); i_p++){
	logger << LOG_DEBUG << "osi.p_parton[" << setw(3) << i_p << "] = " << osi.p_parton[i_a][i_p] << endl;
      }
    }
  }

  // to guarantee resuming runs with identical result (???)
  psi.i_nan++;
  psi.random_manager.nancut();
  psi.i_gen--;
  //
  psi.i_tec++;
  psi.random_manager.increase_cut_counter();
  // questionable !!! Why would this be called here ???
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

// not used at present !!!
void munich::errorhandling_cut_technical(){
  static Logger logger("munich::errorhandling_cut_technical");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi.switch_console_output_phasespace_issue){
    logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "   treated due to cut_technical" << endl;
  }

  bool cut_event = false;

  for (int sr = 0; sr < psi.RA_singular_region_list.size(); sr++){
    int x1 = psi.RA_singular_region_list[sr][0];
    int x2 = psi.RA_singular_region_list[sr][1];
    //    psi.RA_singular_region[x1][x2] = osi.p_parton[0][x1] * osi.p_parton[0][x2] / psi.xbs_all[0][0];
    if (psi.RA_singular_region[x1][x2] < psi.cut_technical){
      cut_event = true;
      if (osi.switch_console_output_phasespace_issue){
	logger << LOG_WARN << right << setw(10) << psi.i_acc << "   int/LO = " << setsr << osi.integrand / osi.sigma_normalization << "    " << psi.RA_singular_region_name[x1][x2] << "/^s = "<< setdr << psi.RA_singular_region[x1][x2] << endl;
      }
    }
  }

  if (cut_event == true){
    psi.i_rej++;
    psi.random_manager.increase_cut_counter();
    if (osi.switch_console_output_phasespace_issue){
      logger << LOG_WARN << "Event counted as cut event, not as (nan/tech)!" << endl;
    }
  }
  else {
    psi.i_nan++;
    psi.random_manager.nancut();
    psi.i_gen--;
    if (osi.switch_console_output_phasespace_issue){
      logger << LOG_WARN << "Event counted as nan event (removed from i_gen)!" << endl;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void munich::errorhandling_c_xbpsp(){
  static Logger logger("munich::errorhandling_c_xbpsp");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi.switch_console_output_phasespace_issue){
    logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "     p != p       @ channel  = " << psi.MC_phasespace.channel << endl;
    for (int i_a = 0; i_a < psi.xbp_all.size(); i_a++){
      for (int i_p = 0; i_p < psi.xbp_all[i_a].size(); i_p++){
	if (psi.xbs_all[i_a][i_p] != 0.){logger << LOG_DEBUG << "psi.xbs_all[" << setw(3) << i_a << "][" << setw(3) << i_p << "] = " << psi.xbs_all[i_a][i_p] << endl;}
	if (psi.xbsqrts_all[i_a][i_p] != 0.){logger << LOG_DEBUG << "psi.xbsqrts_all[" << setw(3) << i_a << "][" << setw(3) << i_p << "] = " << psi.xbsqrts_all[i_a][i_p] << endl;}
	if (psi.xbp_all[i_a][i_p] != nullvector){logger << LOG_DEBUG << "psi.xbp_all[" << setw(3) << i_a << "][" << setw(3) << i_p << "] = " << psi.xbp_all[i_a][i_p] << endl;}
      }
    }
  }

  psi.i_nan++;
  psi.random_manager.nancut();
  psi.i_gen--;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void munich::errorhandling_c_psp_initial(){
  static Logger logger("munich::errorhandling_c_initial");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi.switch_console_output_phasespace_issue){
    logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "     p != p       @ intial momenta generation" << endl;
    for (int i_p = 0; i_p < 3; i_p++){
      logger << LOG_DEBUG << "osi.p_parton[0][" << setw(3) << i_p << "] = " << osi.p_parton[0][i_p] << endl;
    }
  }

  psi.i_nan++;
  psi.random_manager.nancut();
  psi.i_gen--;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void munich::errorhandling_RA_me2(int xswitch){
  static Logger logger("munich::errorhandling_RA_me2");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi.switch_console_output_ME2_issue){
    logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "   MC_phasespace.channel = " << setw(4) << psi.MC_phasespace.channel << "   @ " << dipole[psi.RA_x_a].name() << "   integrand = " << setsr << osi.integrand << "   gtot = " << setsr << psi.g_tot << endl;
    if (osi.switch_console_output_ME2_issue){
      if (xswitch == 0){logger << LOG_WARN << "munich_isnan(sum_RA_ME2)" << endl;}
      if (xswitch == 1){logger << LOG_WARN << "gtot == 0 (for numerical reasons)" << endl;} // why would this happen?
      if (xswitch == 2){logger << LOG_WARN << "RA_ME2_CV != RA_ME2_CV (alpha_S or pdf's are nan)" << endl;}
      if (xswitch == 3){logger << LOG_WARN << "munich_isinf(sum_RA_ME2)" << endl;}
    }
  }

  bool cut_event = false;
  for (int sr = 0; sr < psi.RA_singular_region_list.size(); sr++){
    int x1 = psi.RA_singular_region_list[sr][0];
    int x2 = psi.RA_singular_region_list[sr][1];
    psi.RA_singular_region[x1][x2] = osi.p_parton[0][x1] * osi.p_parton[0][x2] / psi.xbs_all[0][0];
    if (psi.RA_singular_region[x1][x2] < psi.cut_technical){
      cut_event = true;
      if (osi.switch_console_output_ME2_issue){
	logger << LOG_WARN << right << setw(10) << psi.i_acc << "   int/LO = " << setsr << osi.integrand / osi.sigma_normalization << "    " << psi.RA_singular_region_name[x1][x2] << "/^s = "<< setdr << psi.RA_singular_region[x1][x2] << "   A/R = " << setdr << accumulate(osi.RA_ME2.begin() + 1, osi.RA_ME2.end(), 0.) / osi.RA_ME2[0] << endl;
      }
    }
  }

  if (cut_event){
    psi.i_rej++;
    psi.random_manager.increase_cut_counter();
    if (osi.switch_console_output_ME2_issue){
      logger << LOG_WARN << "Event counted as cut event, not as (nan/tech)!" << endl;
    }
  }
  else {
    psi.i_nan++;
    psi.random_manager.nancut();
    psi.i_gen--;
    if (osi.switch_console_output_ME2_issue){
      logger << LOG_WARN << "Event counted as nan event (removed from i_gen)!" << endl;
    }
  }

  // reformulated !!!
  /*
  bool cut_event=false;
  for (int sr = 0; sr < psi.RA_singular_region_list.size(); sr++){
    if (psi.RA_singular_region[psi.RA_singular_region_list[sr][0]][psi.RA_singular_region_list[sr][1]] < psi.cut_technical){
      psi.i_nan--;
      psi.i_gen++;
      psi.i_rej++;
      logger << LOG_WARN << "Event counted as cut event, not as (nan/tech)!" << endl;
      cut_event=true;
      break;
    }
  }

  //	  for (int ib = 0; ib < osi.p_parton[0].size(); ib++){cout << "osi.p_parton[0][" << ib << "] = " << osi.p_parton[0][ib] << setw(25) << osi.p_parton[0][ib].m2() << setw(25) << sqrt(abs(osi.p_parton[0][ib].m2())) << endl;}

  if (cut_event)
    psi.random_manager.increase_cut_counter();
  else
    psi.random_manager.nancut();
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void munich::errorhandling_RA_gtot(){
  static Logger logger("munich::errorhandling_RA_gtot");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi.switch_console_output_phasespace_issue){
    logger << LOG_WARN << "i_gen = " << setw(12) << psi.i_gen << "   i_acc = " << setw(12) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "   MC_phasespace.channel = " << setw(4) << psi.MC_phasespace.channel << "   @ " << dipole[psi.RA_x_a].name() << "   integrand = " << setsr << osi.integrand << "   gtot = " << setsr << psi.g_tot << endl;

    // better use independend n_channel/sum_channel variables !!!
    for (int i_a = 0; i_a < dipole.size(); i_a++){
      for (int jc = 0; jc < dipole[i_a].n_channel(); jc++){
	int temp_zero = 0;
	if (i_a > 0){temp_zero = dipole[i_a - 1].sum_channel();}
	logger << LOG_DEBUG << setw(10) << dipole[i_a].name() << ":   psi.MC_phasespace.g_channel[" << setw(3) << temp_zero << " + " << setw(3) << jc << "] = " << psi.MC_phasespace.g_channel[temp_zero + jc] << endl;
      }
    }
  }

  bool cut_event = false;
  for (int sr = 0; sr < psi.RA_singular_region_list.size(); sr++){
    int x1 = psi.RA_singular_region_list[sr][0];
    int x2 = psi.RA_singular_region_list[sr][1];
    psi.RA_singular_region[x1][x2] = osi.p_parton[0][x1] * osi.p_parton[0][x2] / psi.xbs_all[0][0];
    if (psi.RA_singular_region[x1][x2] < psi.cut_technical){
      cut_event = true;
      if (osi.switch_console_output_phasespace_issue){
	logger << LOG_WARN << right << setw(10) << psi.i_acc << "   int/LO = " << setsr << osi.integrand / osi.sigma_normalization << "    " << psi.RA_singular_region_name[x1][x2] << "/^s = "<< setdr << psi.RA_singular_region[x1][x2] << "   A/R = " << setdr << accumulate(osi.RA_ME2.begin() + 1, osi.RA_ME2.end(), 0.) / osi.RA_ME2[0] << endl;
      }
    }
  }

  if (cut_event){
    psi.i_rej++;
    psi.random_manager.increase_cut_counter();
    if (osi.switch_console_output_ME2_issue){
      logger << LOG_WARN << "Event counted as cut event, not as (nan/tech)!" << endl;
    }
  }
  else {
    psi.i_nan++;
    psi.random_manager.nancut();
    psi.i_gen--;
    if (osi.switch_console_output_ME2_issue){
      logger << LOG_WARN << "Event counted as nan event (removed from i_gen)!" << endl;
    }
  }


  /*
  bool cut_event = false;
  for (int sr = 0; sr < psi.RA_singular_region_list.size(); sr++){
    if (psi.RA_singular_region[psi.RA_singular_region_list[sr][0]][psi.RA_singular_region_list[sr][1]] < psi.cut_technical){
      psi.i_nan--;
      psi.i_gen++;
      psi.i_rej++;
      if (osi.switch_console_output_phasespace_issue){
	logger << LOG_WARN << "Event counted as cut event, not as (nan/tech)!" << endl;
      }
      cut_event=true;
      break;
    }
  }
  //	  for (int ib = 0; ib < osi.p_parton[0].size(); ib++){cout << "osi.p_parton[0][" << ib << "] = " << osi.p_parton[0][ib] << setw(25) << osi.p_parton[0][ib].m2() << setw(25) << sqrt(abs(osi.p_parton[0][ib].m2())) << endl;}
  psi.i_nan++;
  psi.i_gen--;
  if (cut_event)
    psi.random_manager.increase_cut_counter();
  else
    psi.random_manager.nancut();
  */

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void munich::errorhandling_collinear_me2(){
  static Logger logger("munich::errorhandling_collinear_me2");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi.switch_console_output_ME2_issue){
    logger << LOG_WARN << psi.i_gen << setw(10) << psi.i_acc << "   i_rej = " << setw(12) << psi.i_rej << "   (" << setw(5) << psi.i_nan << ")" << "   MC_phasespace.channel = " << setw(4) << psi.MC_phasespace.channel << "   integrand = " << setprecision(16) << setw(25) << osi.integrand << "   gtot = " << setprecision(16) << setw(25) << psi.g_tot << endl;
    logger << LOG_DEBUG << psi.i_nan << " nancuts   (ME2/gtot)" << endl;
    for (int i_p = 1; i_p < osi.p_parton[0].size(); i_p++){
      logger << LOG_DEBUG << "osi.p_parton[0][" << i_p << "] = " << osi.p_parton[0][i_p] << endl;
    }
    logger << LOG_DEBUG << "gtot = " << psi.g_tot << endl;
    for (int i_c = 1; i_c < 3; i_c++){
      logger << LOG_DEBUG << "g_z_coll[" << i_c << "] = " << psi.g_z_coll[i_c] << endl;
    }
    for (int i_a = 0; i_a < collinear.size(); i_a++){
      for (int j_a = 0; j_a < collinear[i_a].size(); j_a++){
	logger << LOG_DEBUG << "CA_ME2_cf[" << i_a << "][" << j_a << "] = " << osi.CA_ME2_cf[i_a][j_a] << endl;
      }
    }
    logger << LOG_DEBUG << "integrand = " << osi.integrand << endl;
    logger << LOG_DEBUG << endl;
  }

  psi.i_nan++;
  psi.random_manager.nancut();
  psi.i_gen--;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


