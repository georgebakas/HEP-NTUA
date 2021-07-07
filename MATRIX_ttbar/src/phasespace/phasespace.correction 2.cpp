#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

void phasespace_set::correct_phasespacepoint_real(){
  Logger logger("phasespace_set::correct_phasespacepoint_real");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  /*
  corrected_xbp_all = start_xbp_all;
  corrected_xbs_all = start_xbs_all;
  corrected_xbsqrts_all = start_xbsqrts_all;
  */

  fourvector momentum_correction_violation(-xbsqrts_all[0][0], 0., 0., 0.);
  for (int xbi = 4; xbi < xbp_all[0].size(); xbi = xbi * 2){momentum_correction_violation = momentum_correction_violation + xbp_all[0][xbi];}
  //cout << "momentum_correction_violation = " << momentum_correction_violation << endl;
  momentum_correction_violation = momentum_correction_violation / csi->n_particle;
  
  if (momentum_correction_violation.r() > 1.e-16 * xbsqrts_all[0][0]){
    //    logger << LOG_DEBUG << "i_acc = " << setw(8) << i_acc << "   i_gen = " << setw(8) << i_gen << "   momentum corrected:   " << momentum_correction_violation.r() << endl;
    logger << LOG_DEBUG_VERBOSE << "i_acc = " << setw(8) << right << i_acc << "   i_gen = " << setw(8) << right << i_gen << "   momentum corrected:   " << left << setprecision(15) << setw(23) << momentum_correction_violation.r() << endl;

    corrected_xbp_all[0] = start_xbp_all[0];
    corrected_xbs_all[0] = start_xbs_all[0];
    corrected_xbsqrts_all[0] = start_xbsqrts_all[0];

    //cout << "momentum_correction_violation = " << momentum_correction_violation << endl;
    for (int xbi = 4; xbi < xbp_all[0].size(); xbi = xbi * 2){
      fourvector temp = xbp_all[0][xbi] - momentum_correction_violation;
      corrected_xbp_all[0][xbi] = fourvector(sqrt(xbs_all[0][xbi] + temp.r2()), temp.x1(), temp.x2(), temp.x3());
    }
    double corrected_E = 0.;
    double corrected_E_beam = 0.;
    
    
    logger << LOG_DEBUG_VERBOSE << "before psp correction:" << endl;
    for (int xbi = 0; xbi < xbp_all[0].size(); xbi++){
      logger << LOG_DEBUG_VERBOSE << "xbp_all[0][" << setw(3) << xbi << "] = " << xbp_all[0][xbi] << "   xbs_all[0][" << setw(3) << xbi << "] = " << xbs_all[0][xbi] << endl;
    }
    
    for (int xbi = 4; xbi < xbp_all[0].size(); xbi = xbi * 2){corrected_E += xbp_all[0][xbi].x0();}
    
    int binary_all_out = intpow(2, csi->n_particle + 2) - 4;
    
    //cout << "corrected_E = " << corrected_E << endl;
    corrected_xbp_all[0][0] = fourvector(corrected_E, 0., 0., 0.);
    corrected_xbp_all[0][binary_all_out] = corrected_xbp_all[0][0];
    
    corrected_xbs_all[0][0] = pow(corrected_E, 2);
    corrected_xbs_all[0][binary_all_out] = pow(corrected_E, 2);
    corrected_xbsqrts_all[0][0] = corrected_E;
    corrected_xbsqrts_all[0][binary_all_out] = corrected_E;
    
    corrected_E_beam = corrected_E / 2;
    if (o_map[0][1] == 1 && o_map[0][2] == 2){
      corrected_xbp_all[0][1] = fourvector(corrected_E_beam, 0., 0., corrected_E_beam);
      corrected_xbp_all[0][2] = fourvector(corrected_E_beam, 0., 0., -corrected_E_beam);
    }
    else if (o_map[0][1] == 2 && o_map[0][2] == 1){
      corrected_xbp_all[0][1] = fourvector(corrected_E_beam, 0., 0., -corrected_E_beam);
      corrected_xbp_all[0][2] = fourvector(corrected_E_beam, 0., 0., corrected_E_beam);
    }
    
    
    xbp_all[0] = corrected_xbp_all[0];
    xbs_all[0] = corrected_xbs_all[0];
    xbsqrts_all[0] = corrected_xbsqrts_all[0];
    /*
    xbp_all = corrected_xbp_all;
    xbs_all = corrected_xbs_all;
    xbsqrts_all = corrected_xbsqrts_all;
    */
    
    logger << LOG_DEBUG_VERBOSE << "after psp correction:" << endl;
    for (int xbi = 0; xbi < xbp_all[0].size(); xbi++){
      logger << LOG_DEBUG_VERBOSE << "xbp_all[0][" << setw(3) << xbi << "] = " << xbp_all[0][xbi] << "   xbs_all[0][" << setw(3) << xbi << "] = " << xbs_all[0][xbi] << endl;
    }
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void phasespace_set::correct_phasespacepoint_dipole(){
  Logger logger("phasespace_set::correct_phasespacepoint_dipole");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_a = 1; i_a < n_dipoles; i_a++){
    
    //  fourvector momentum_correction_violation(-xbsqrts_all[i_a][0], 0., 0., 0.);
    fourvector momentum_correction_violation(0., 0., 0., 0.);
    for (int xbi = 1; xbi < 4; xbi = xbi * 2){
      momentum_correction_violation = momentum_correction_violation - xbp_all[i_a][xbi];
    }
    for (int xbi = 4; xbi < xbp_all[i_a].size(); xbi = xbi * 2){
      momentum_correction_violation = momentum_correction_violation + xbp_all[i_a][xbi];
    }
    //  cout << "momentum_correction_violation = " << momentum_correction_violation << endl;
    momentum_correction_violation = momentum_correction_violation / (csi->n_particle - 1);
    
    if (momentum_correction_violation.r() > 1.e-16 * xbsqrts_all[i_a][0]){
      logger << LOG_DEBUG_VERBOSE << "i_acc = " << setw(8) << right << i_acc << "   i_gen = " << setw(8) << right << i_gen << "   momentum corrected at i_a = " << i_a << ":   " << left << setprecision(15) << setw(23) << momentum_correction_violation.r() << endl;

      corrected_xbp_all[i_a] = start_xbp_all[i_a];
      //      corrected_xbs_all[i_a] = start_xbs_all[i_a];
      //      corrected_xbsqrts_all[i_a] = start_xbsqrts_all[i_a];

      //  cout << "momentum_correction_violation / " << (csi->n_particle - 1) << " = " << momentum_correction_violation << endl;
      for (int xbi = 4; xbi < xbp_all[i_a].size(); xbi = xbi * 2){
	fourvector temp = xbp_all[i_a][xbi] - momentum_correction_violation;
	corrected_xbp_all[i_a][xbi] = fourvector(sqrt(xbs_all[i_a][xbi] + temp.r2()), temp.x1(), temp.x2(), temp.x3());
      }
      for (int xbi = 4; xbi < xbp_all[i_a].size(); xbi = xbi * 2){
	double prec = abs(1 - sqrt(xbs_all[i_a][xbi] + xbp_all[i_a][xbi].r2()) / xbp_all[i_a][xbi].x0());
	double prec_new = abs(1 - sqrt(xbs_all[i_a][xbi] + corrected_xbp_all[i_a][xbi].r2()) / corrected_xbp_all[i_a][xbi].x0());
	if (prec_new > 1.e-9){
	  logger << LOG_DEBUG << "prec     = " << prec << endl;
	  logger << LOG_DEBUG << "prec_new = " << prec_new << "   " << corrected_xbp_all[i_a][xbi].r2() << "   " << sqrt(xbs_all[i_a][xbi] + corrected_xbp_all[i_a][xbi].r2()) << "   " << corrected_xbp_all[i_a][xbi].x0() << endl;
	  logger << LOG_DEBUG << "On-shell condition badly satisfied for xbp_all[" << i_a << "][" << xbi << "], namely to       " << right << setw(25) << setprecision(16) << showpoint << -log10(prec) << " digits." << endl;
	  logger << LOG_DEBUG << "                                                               new:   " << right << setw(25) << setprecision(16) << showpoint << -log10(prec_new) << " digits." << endl;
	  logger << LOG_DEBUG << "          xbp_all[" << i_a << "][" << setw(3) << xbi << "] = " << xbp_all[i_a][xbi] << "   m = " << sqrt(abs(xbp_all[i_a][xbi].m2())) << "   xbs_all[" << i_a << "][" << setw(3) << xbi << "] = " << xbsqrts_all[i_a][xbi] << endl;
	  logger << LOG_DEBUG << "corrected_xbp_all[" << i_a << "][" << setw(3) << xbi << "] = " << corrected_xbp_all[i_a][xbi] << "   m = " << sqrt(abs(corrected_xbp_all[i_a][xbi].m2())) << endl;
	}
	if (prec_new > 1.e-9){logger << LOG_DEBUG << "On-shell condition still badly satisfied for corrected_xbp_all[" << i_a << "][" << xbi << "], namely to " << right << setw(25) << setprecision(16) << showpoint << -log10(prec_new) << " digits." << endl;}
      }
      for (int xbi = 4; xbi < xbp_all[i_a].size(); xbi = xbi * 2){momentum_correction_violation = momentum_correction_violation + xbp_all[i_a][xbi];}
      for (int xbi = 4; xbi < xbp_all[i_a].size(); xbi = xbi * 2){xbp_all[i_a][xbi] = corrected_xbp_all[i_a][xbi];}
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
