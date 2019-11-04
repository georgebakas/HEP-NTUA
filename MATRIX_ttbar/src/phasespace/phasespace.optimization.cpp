#include "../include/classes.cxx"

//#include "old.phasespace.optimization.cpp"

void phasespace_set::psp_MCweight_optimization(double & integrand, double & this_psp_weight, double & this_psp_weight2, double & optimization_modifier){
  Logger logger("phasespace_set::psp_MCweight_optimization");
  
  double mod_integrand = integrand * optimization_modifier;
  double mod_this_psp_weight = this_psp_weight * optimization_modifier;
  double mod_this_psp_weight2 = this_psp_weight2 * pow(optimization_modifier, 2);
  logger << LOG_DEBUG << "integrand = " <<integrand<< endl;
  // MC_phasespace - multichannel for final-state phasespace
  ///  MC_phasespace.channel = MC_channel;
  MC_phasespace.psp_MCweight_optimization(mod_integrand, g_tot);
  ///  counts_channel = MC_phasespace.n_acc_channel;
   
  // MC_tau - multichannel for partonic CMS energy tau
  MC_tau.psp_MCweight_optimization(mod_integrand, g_tot);


  //  logger << LOG_INFO << "RA_x_a = " << RA_x_a << endl;
  // RA_x_a is the dipole whose channel has been chosen for phase-space generation.

  if (csi->class_contribution_CS_real){
    // MC_x_dipole[i_a] - multichannel for partonic CMS energy in dipole kinematics
    //    if (RA_x_a > 0){MC_x_dipole[RA_x_a].psp_MCweight_optimization(mod_integrand, g_tot);}
    for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
      // check if this should be done only for the used channel RA_x_a !!!
      if (MC_x_dipole[i_a].active_optimization != -1){
	MC_x_dipole[i_a].psp_MCweight_optimization(mod_integrand, g_tot);
      }
    }
  }

  if (!IS_tau.end_optimization){
    IS_tau.psp_IS_optimization(mod_this_psp_weight, mod_this_psp_weight2);
  }
  if (!IS_x1x2.end_optimization){
    IS_x1x2.psp_IS_optimization(mod_this_psp_weight, mod_this_psp_weight2);
  }
  /*
  if (tau_opt_end == 0){
    tau_counts_channel[tau_channel]++;
    sum_tau_channel_weight[tau_channel] += abs(mod_this_psp_weight);
    sum_tau_channel_weight2[tau_channel] += mod_this_psp_weight2;
  }
  
  if (x1x2_opt_end == 0){
    x1x2_counts_channel[x1x2_channel]++;
    sum_x1x2_channel_weight[x1x2_channel] += abs(mod_this_psp_weight);
    sum_x1x2_channel_weight2[x1x2_channel] += mod_this_psp_weight2;
  }
  */
  
  random_manager.psp_IS_optimization(mod_this_psp_weight, mod_this_psp_weight2);
  //  random_manager.optimize(mod_this_psp_weight, mod_this_psp_weight2);

  if (munich_isnan(mod_integrand)){logger << LOG_FATAL << "unidentified error - mod_integrand nan - run stopped!" << endl; exit(1);}
  if (munich_isinf(mod_integrand)){logger << LOG_FATAL << "unidentified error - mod_integrand inf - run stopped!" << endl; exit(1);}
  if (munich_isnan(g_tot)){logger << LOG_FATAL << "unidentified error - gtot nan - run stopped!" << endl; exit(1);}
  //  if (munich_isinf(g_tot)){logger << LOG_DEBUG << "unidentified error - gtot inf - run stopped!" << endl; exit(1);}

  if (csi->class_contribution_CS_collinear){
    for (int i_z = 1; i_z < 3; i_z++){
      if (!IS_z1z2[i_z].end_optimization){
	IS_z1z2[i_z].psp_IS_optimization(this_psp_weight, this_psp_weight2);
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::step_MCweight_optimization(){
  Logger logger("phasespace_set::step_MCweight_optimization");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // MC_phasespace - multichannel for final-state phasespace
  MC_phasespace.step_MCweight_optimization(*i_step_mode);

  // MC_tau - multichannel for partonic CMS energy tau
  MC_tau.step_MCweight_optimization(*i_step_mode);

  if (csi->class_contribution_CS_real){
    // MC_x_dipole[i_a] - multichannel for partonic CMS energy in dipole kinematics
    //    for (int i_a = 1; i_a < n_dipoles; i_a++){
    for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
      logger << LOG_DEBUG << "MC_x_dipole[" << i_a << "].n_channel = " << MC_x_dipole[i_a].n_channel << endl;
      MC_x_dipole[i_a].step_MCweight_optimization(*i_step_mode);
    }
  }
 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::result_MCweight_optimization(){
  Logger logger("phasespace_set::result_MCweight_optimization");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  /*/////
  if (MC_opt_end != 0){return;}
  if (*i_step_mode != n_alpha_events * n_alpha_steps){return;}

  ///  MC_opt_end = 1;

  //  string dummy = "";

  // MC_phasespace - multichannel for final-state phasespace
  MC_phasespace.result_MCweight_optimization(*i_step_mode);
  //, i_rej, dummy);

  
  // not active for some strange reason: should be reactivated !!!
  // test reactivation !!!
  // MC_tau - multichannel for partonic CMS energy tau
  MC_tau.result_MCweight_optimization(*i_step_mode);

  if (csi->class_contribution_CS_real){
    // MC_x_dipole[i_a] - multichannel for partonic CMS energy in dipole kinematics
    for (int i_a = 1; i_a < n_dipoles; i_a++){
      MC_x_dipole[i_a].result_MCweight_optimization(*i_step_mode);
    }
  }
  */

  
  /*
  // All multichannel optimizations must be synchronized:
  MC_opt_end = MC_phasespace.end_optimization;

  // reactivate result_... for MC_tau and MC_x_dipole[i_a] ???
  
  logger << LOG_DEBUG_VERBOSE << "before weight_output" << endl;

  //  weight_output_MCweight_optimization();
  */
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


/// might be needed in 'observable_set::initialization_integration' and 'observable_set::output_zero_contribution_complete' !!!
void phasespace_set::weight_output_MCweight_optimization(){
  Logger logger("phasespace_set::weight_output_MCweight_optimization");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  /*
  //  ofstream out_MCweights(filename_MCweight.c_str());
  ofstream out_MCweights;
  //  out_integration.open(filename_MCweight.c_str(), ofstream::out | ofstream::app);
  out_MCweights.open(filename_MCweight.c_str(), ofstream::out | ofstream::trunc);
  out_MCweights << "# " << "MC" << endl;
  for (int j = 0; j < MC_alpha.size(); j++){
    out_MCweights << setprecision(16) << MC_alpha[j] << endl;
  }
  out_MCweights.close();
  MC_opt_end = 1; // is done twice !!! needed ???
  ///  size_proceeding[2] = 0;
  ///  size_proc_generic = accumulate(size_proceeding.begin(), size_proceeding.end(), 0);
  */

  // everything above could be done for "MC_phasespace" !!!


  ofstream out_MCweights;
  out_MCweights.open(filename_MCweight.c_str(), ofstream::out | ofstream::trunc);
  out_MCweights.close();

  if (switch_MC != -1){
    MC_phasespace.output_MCweight_optimization(*i_step_mode, filename_MCweight);
  }


  
  // Where does the "result_MCweight_optimization" happen for MC_tau and MC_x_dipole[i_a] ??? !!!

  // MC_tau - multichannel for partonic CMS energy tau
  //  MC_tau.result_MCweight_optimization(i_acc, i_rej, filename_MCweight);
  if (switch_MC_tau != -1){
    MC_tau.output_MCweight_optimization(*i_step_mode, filename_MCweight);
  }

  if ((switch_MC_x_dipole != -1) && csi->class_contribution_CS_real){
    // MC_x_dipole[i_a] - multichannel for partonic CMS energy in dipole kinematics
    //    for (int i_a = 1; i_a < n_dipoles; i_a++){
    for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
      MC_x_dipole[i_a].output_MCweight_optimization(*i_step_mode, filename_MCweight);
      //      MC_x_dipole[i_a].result_MCweight_optimization(i_acc, i_rej, filename_MCweight);
    }
  }
 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




void phasespace_set::step_optimization_complete(){
  Logger logger("phasespace_set::step_optimization_complete");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // MC_phasespace - multichannel for final-state phasespace
  MC_phasespace.step_MCweight_optimization(*i_step_mode);

  // MC_tau - multichannel for partonic CMS energy tau
  MC_tau.step_MCweight_optimization(*i_step_mode);

  if (csi->class_contribution_CS_real){
    // MC_x_dipole[i_a] - multichannel for partonic CMS energy in dipole kinematics
    for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
      MC_x_dipole[i_a].step_MCweight_optimization(*i_step_mode);
    }    
  }

  // MC_opt_end should be replaced by something general like 'end_optimization'  !!!
  //  if (!MC_opt_end){step_MCweight_optimization();}

  random_manager.do_optimization_step(*i_step_mode);

  IS_tau.step_IS_optimization(*i_step_mode);
  
  IS_x1x2.step_IS_optimization(*i_step_mode);
  
  if (csi->class_contribution_CS_collinear){
    for (int i_z = 1; i_z < 3; i_z++){
      IS_z1z2[i_z].step_IS_optimization(*i_step_mode);
    }
  }

  // What about VT and VT2 ???
  
  if (switch_use_alpha_after_IS){
    // choose best alpha set only after previous IS optimization step:
    if (switch_IS_MC && !MC_opt_end && *i_step_mode % n_IS_events == 0){

      MC_phasespace.x_alpha_it_min = MC_phasespace.i_alpha_it;
      
      if (MC_tau.active_optimization){
	MC_tau.x_alpha_it_min = MC_tau.i_alpha_it;
      }

      if ((switch_MC_x_dipole != -1) && csi->class_contribution_CS_real){
	for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
	  if (MC_x_dipole[i_a].active_optimization){
	    MC_x_dipole[i_a].x_alpha_it_min = MC_x_dipole[i_a].i_alpha_it;
	  }
	}
      }
    }
  }
    
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


    
void phasespace_set::psp_optimization_complete(double & integrand, double & this_psp_weight, double & this_psp_weight2, double & optimization_modifier){
  Logger logger("phasespace_set::psp_optimization_complete");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  double mod_integrand = integrand * optimization_modifier;
  double mod_this_psp_weight = this_psp_weight * optimization_modifier;
  double mod_this_psp_weight2 = this_psp_weight2 * pow(optimization_modifier, 2);

  // MC_phasespace - multichannel for final-state phasespace
  MC_phasespace.psp_MCweight_optimization(mod_integrand, g_tot);
   
  // MC_tau - multichannel for partonic CMS energy tau
  MC_tau.psp_MCweight_optimization(mod_integrand, g_tot);

  //  logger << LOG_INFO << "RA_x_a = " << RA_x_a << endl;
  // RA_x_a is the dipole whose channel has been chosen for phase-space generation.

  if (csi->class_contribution_CS_real){
    // MC_x_dipole[i_a] - multichannel for partonic CMS energy in dipole kinematics
    //    if (RA_x_a > 0){MC_x_dipole[RA_x_a].psp_MCweight_optimization(mod_integrand, g_tot);}
    for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
      // check if this should be done only for the used channel RA_x_a !!!
      MC_x_dipole[i_a].psp_MCweight_optimization(mod_integrand, g_tot);
    }
  }

  random_manager.psp_IS_optimization(mod_this_psp_weight, mod_this_psp_weight2);
  //  random_manager.optimize(mod_this_psp_weight, mod_this_psp_weight2);

  IS_tau.psp_IS_optimization(mod_this_psp_weight, mod_this_psp_weight2);

  IS_x1x2.psp_IS_optimization(mod_this_psp_weight, mod_this_psp_weight2);
  
  if (csi->class_contribution_CS_collinear){
    for (int i_z = 1; i_z < 3; i_z++){
      IS_z1z2[i_z].psp_IS_optimization(mod_this_psp_weight, mod_this_psp_weight2);
    }
  }

  if (munich_isnan(mod_integrand)){logger << LOG_FATAL << "unidentified error - mod_integrand nan - run stopped!" << endl; exit(1);}
  if (munich_isinf(mod_integrand)){logger << LOG_FATAL << "unidentified error - mod_integrand inf - run stopped!" << endl; exit(1);}
  if (munich_isnan(g_tot)){logger << LOG_FATAL << "unidentified error - gtot nan - run stopped!" << endl; exit(1);}
  //  if (munich_isinf(g_tot)){logger << LOG_DEBUG << "unidentified error - gtot inf - run stopped!" << endl; exit(1);}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::result_optimization_complete(){
  Logger logger("phasespace_set::result_optimization_complete");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (MC_opt_end == 0){
    // MC_phasespace - multichannel for final-state phasespace
    MC_phasespace.result_MCweight_optimization(*i_step_mode);
    
    // was not active for some strange reason: should be reactivated !!!
    
    // MC_tau - multichannel for partonic CMS energy tau
    MC_tau.result_MCweight_optimization(*i_step_mode);
    
    if (csi->class_contribution_CS_real){
      logger << LOG_DEBUG << "MC_x_dipole.size() = " << MC_x_dipole.size() << endl;
      for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
	//      for (int i_a = 1; i_a < n_dipoles; i_a++){
	if (MC_x_dipole[i_a].active_optimization != -1){
	  // MC_x_dipole[i_a] - multichannel for partonic CMS energy in dipole kinematics
	  MC_x_dipole[i_a].result_MCweight_optimization(*i_step_mode);
	}
      }
    }
    
    // All multichannel optimizations are synchronized (nothing else tested...)

    MC_opt_end = 1;
    if (MC_phasespace.active_optimization == 1 && MC_phasespace.end_optimization != 1){MC_opt_end = 0;}
    if (MC_tau.active_optimization == 1 && MC_tau.end_optimization != 1){MC_opt_end = 0;}
    if (csi->class_contribution_CS_real){
      for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
	//      for (int i_a = 1; i_a < n_dipoles; i_a++){
	if (MC_x_dipole[i_a].active_optimization == 1 && MC_x_dipole[i_a].end_optimization != 1){MC_opt_end = 0;}
      }
    }
  }

  IS_tau.result_IS_optimization(*i_step_mode);

  IS_x1x2.result_IS_optimization(*i_step_mode);

  if (csi->class_contribution_CS_collinear){
    for (int i_z = 1; i_z < 3; i_z++){
      IS_z1z2[i_z].result_IS_optimization(*i_step_mode);
    }
  }

  int temp_IS_opt_end = 1;
  if (IS_tau.active_optimization == 1 && IS_tau.end_optimization != 1){temp_IS_opt_end = 0;}
  if (IS_x1x2.active_optimization == 1 && IS_x1x2.end_optimization != 1){temp_IS_opt_end = 0;}
  if (csi->class_contribution_CS_collinear){
    for (int i_z = 1; i_z < 3; i_z++){
      if (IS_z1z2[i_z].active_optimization == 1 && IS_z1z2[i_z].end_optimization != 1){temp_IS_opt_end = 0;}
    }
  }
  
  if (end_optimization == 0){
    end_optimization = 1;
    if (MC_opt_end != 1){end_optimization = 0;}
    if (temp_IS_opt_end != 1){end_optimization = 0;}
  }
  
  /*
  end_optimization = 1;
  
  if (!MC_phasespace.end_optimization){end_optimization = 0;}
  if (!MC_tau.end_optimization){end_optimization = 0;}
  if (csi->class_contribution_CS_real){
    for (int i_a = 1; i_a < n_dipoles; i_a++){
      if (!MC_x_dipole[i_a].end_optimization){end_optimization = 0;}
    }
  }
  // add IS stuff here !!!  
  */
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void phasespace_set::output_optimization_complete(){
  Logger logger("phasespace_set::output_optimization_complete");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (MC_opt_end == 1){
    ofstream out_MCweights;
    out_MCweights.open(filename_MCweight.c_str(), ofstream::out | ofstream::trunc);
    out_MCweights.close();
    
    // MC_phasespace - multichannel for final-state phasespace
    MC_phasespace.output_MCweight_optimization(*i_step_mode, filename_MCweight);
    
    // MC_tau - multichannel for partonic CMS energy tau
    MC_tau.output_MCweight_optimization(*i_step_mode, filename_MCweight);
    
    if (csi->class_contribution_CS_real){
      for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
	if (MC_x_dipole[i_a].active_optimization != -1){
	  // MC_x_dipole[i_a] - multichannel for partonic CMS energy in dipole kinematics
	  MC_x_dipole[i_a].output_MCweight_optimization(*i_step_mode, filename_MCweight);
	}
      }
    }
    
    MC_opt_end = 2;
    if (MC_phasespace.active_optimization == 1 && MC_phasespace.end_optimization != 2){MC_opt_end = 1;}
    if (MC_tau.active_optimization == 1 && MC_tau.end_optimization != 2){MC_opt_end = 1;}
    if (csi->class_contribution_CS_real){
      //      for (int i_a = 1; i_a < n_dipoles; i_a++){
      for (int i_a = 1; i_a < MC_x_dipole.size(); i_a++){
	if (MC_x_dipole[i_a].active_optimization == 1 && MC_x_dipole[i_a].end_optimization != 2){MC_opt_end = 1;}
      }
    }
  }
  
  // IS_tau - importance sampling for partonic CMS energy tau
  IS_tau.output_IS_optimization(*i_step_mode);

  // IS_x1x2 - importance sampling for partonic momentum fractions x1 and x2
  IS_x1x2.output_IS_optimization(*i_step_mode);

  if (csi->class_contribution_CS_collinear){
    for (int i_z = 1; i_z < 3; i_z++){
      // IS_z1z2[i_z] - importance sampling for collinear emission from incoming parton i_z
      IS_z1z2[i_z].output_IS_optimization(*i_step_mode);
    }
  }
  
  // random_manager
  if (random_manager.end_optimization){random_manager.writeout_weights();}
  
  
  int temp_IS_opt_end = 2;
  if (IS_tau.active_optimization == 1 && IS_tau.end_optimization != 2){temp_IS_opt_end = 1;}
  if (IS_x1x2.active_optimization == 1 && IS_x1x2.end_optimization != 2){temp_IS_opt_end = 1;}
  if (csi->class_contribution_CS_collinear){
    for (int i_z = 1; i_z < 3; i_z++){
      if (IS_z1z2[i_z].active_optimization == 1 && IS_z1z2[i_z].end_optimization != 2){temp_IS_opt_end = 1;}
    }
  }
  
  if (end_optimization == 1){
    end_optimization = 2;
    if (MC_opt_end != 2){end_optimization = 1;}
    if (temp_IS_opt_end != 2){end_optimization = 1;}
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


