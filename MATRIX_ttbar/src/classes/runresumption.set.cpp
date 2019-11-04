#include "../include/classes.cxx"

////////////////////
//  constructors  //
////////////////////
runresumption_set::runresumption_set(){
  size_proc_generic = 0;
  size_proceeding.resize(6);
}
runresumption_set::runresumption_set(observable_set & _osi, phasespace_set & _psi, observable_set & _save_osi, phasespace_set & _save_psi){

  size_proc_generic = 0;
  size_proceeding.resize(6);

  osi = &_osi;
  psi = &_psi;

  resumption_osi = &_save_osi;
  resumption_psi = &_save_psi;

  (*resumption_psi).random_manager.proceeding_save();
}

void runresumption_set::perform_proceeding_step(){
  Logger logger("runresumption_set::perform_proceeding_step");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  *resumption_osi = *osi;
  *resumption_psi = *psi;
  (*resumption_psi).random_manager.proceeding_save();
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void runresumption_set::perform_proceeding_in(){
  Logger logger("runresumption_set::perform_proceeding_in");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
  (*osi).perform_proceeding_in(*psi);
    
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void runresumption_set::perform_proceeding_out(){
  Logger logger("runresumption_set::perform_proceeding_out");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
  // Should depend on value of switch_proceeding when used !!!
  (*resumption_osi).perform_proceeding_out(*resumption_psi);
  //  (*osi).perform_proceeding_out(*psi);
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
 
}
void runresumption_set::perform_proceeding_check(){
  Logger logger("runresumption_set::perform_proceeding_check");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
  (*osi).perform_proceeding_check(*psi);
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void runresumption_set::perform_iteration_step(){
  Logger logger("runresumption_set::perform_iteration_step");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  ///  static int last_step_mode = 0;
  logger << LOG_DEBUG_VERBOSE << "*psi->i_step_mode = " << setw(12) << *psi->i_step_mode << "   "
	 << setw(12) << psi->i_gen << " [" << setw(3) << psi->i_nan << "] " << " "
	 << setw(10) << psi->i_acc << " (" << setw(7) << psi->i_tec << ") " << " "
	 << setw(10) << psi->i_rej << " "
	 << endl;
  
  if (*psi->i_step_mode > psi->last_step_mode && *psi->i_step_mode % psi->n_step == 0){
    if (osi->switch_output_proceeding == 2){perform_proceeding_out();}
    osi->calculate_intermediate_result(*psi);

    logger << LOG_DEBUG_VERBOSE << "osi->perform_integration_step_complete" << endl;
    osi->perform_integration_step_complete(*psi);
    if (osi->switch_output_distribution){osi->output_distribution_complete(*psi);}
    
    if (osi->csi->type_contribution == "RVA" ||
	osi->csi->type_contribution == "L2RT" ||
	osi->csi->type_contribution == "L2RJ"){
      logger << LOG_INFO << "Information about killed phase-space points at " << psi->i_acc << " events:" << endl;
      for (int i_q = 0; i_q < osi->n_qTcut; i_q++){
	if (osi->counter_acc_qTcut[i_q] > 0){
	  logger << LOG_INFO << "counter_killed_qTcut[" << setw(3) << i_q << "] = " << setw(8) << osi->counter_killed_qTcut[i_q] << "   counter_acc_qTcut[" << setw(3) << i_q << "] = " << setw(8) << osi->counter_acc_qTcut[i_q] << "   ratio = " << double(osi->counter_killed_qTcut[i_q]) / double(osi->counter_acc_qTcut[i_q]) << endl;
	}
      }
    }

    // introduce global 'active_optimization' variable to avoid weight output in normal (non-grid) runs !!!
    
    psi->step_optimization_complete();
    // to be implemented !!! better distinction for normal (non-grid) runs !!!
    // maybe psi->end_optimization = 0, 1, 2 ???

    psi->result_optimization_complete();

    logger << LOG_DEBUG_VERBOSE << "psi->MC_phasespace.end_optimization = " << psi->MC_phasespace.end_optimization << endl;

    output_weight_optimization(*psi, *osi, size_proc_generic, size_proceeding);

    logger << LOG_DEBUG_VERBOSE << "psi->MC_phasespace.end_optimization = " << psi->MC_phasespace.end_optimization << endl;

    psi->output_optimization_complete();

    
    /*
    if (!psi->end_optimization){psi->result_optimization_complete();}
    if (psi->end_optimization){psi->output_optimization_complete();}
    */
  /*
  psi->step_MCweight_optimization();
  ///  psi->random_manager.do_optimization_step(psi->i_acc);
  ///  psi->random_manager.do_optimization_step(psi->i_gen);
  psi->random_manager.do_optimization_step(*psi->i_step_mode);
  psi->result_MCweight_optimization();
  psi->step_tau_weight_optimization();
  psi->step_x1x2_weight_optimization();
  if (osi->csi->type_contribution == "CA" || 
      osi->csi->type_contribution == "RCA" || 
      osi->csi->type_contribution == "L2CA"){
    psi->step_z1z2_weight_optimization();
    // for CT, CT2 this is handled differently.
  }
  */
    
    ///    output_weight_optimization(*psi, *osi, size_proc_generic, size_proceeding);
    if (osi->switch_output_proceeding){perform_proceeding_step();}
    if (osi->switch_output_proceeding == 1){perform_proceeding_out();}

    psi->last_step_mode = *psi->i_step_mode;
  }

  if (*psi->i_step_mode == psi->n_events_max){osi->int_end = 1;}
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
