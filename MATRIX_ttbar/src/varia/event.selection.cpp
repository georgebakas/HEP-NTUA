#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

void perform_event_selection(observable_set & oset, call_generic & generic){
  static Logger logger("perform_event_selection");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //      if (osi_switch_output_cutinfo){
  //	logger.newLine(LOG_DEBUG);
  //	logger.newLine(LOG_DEBUG);
  ///  logger.newLine(LOG_DEBUG);
  logger.newLine(LOG_DEBUG_POINT);
  logger << LOG_DEBUG // << "[" << setw(2) << i_a << "]"
	 << " New PSP: " 
	 << "   i_gen = " << setw(12) << oset.psi->i_gen 
	 << "   i_acc = " << setw(12) << oset.psi->i_acc 
	 << "   i_rej = " << setw(12) << oset.psi->i_rej 
	 << "   i_nan = " << setw(4) << oset.psi->i_nan 
	 << "   i_tec = " << setw(4) << oset.psi->i_tec 
	 << endl;
  //  logger.newLine(LOG_DEBUG);

  /*
  if (oset.csi->type_contribution == "VT2"){
    cerr << " New PSP: " 
	 << "   i_gen = " << setw(12) << oset.psi->i_gen 
	 << "   i_acc = " << setw(12) << oset.psi->i_acc 
	 << "   i_rej = " << setw(12) << oset.psi->i_rej 
	 << "   i_nan = " << setw(4) << oset.psi->i_nan 
	 << "   i_tec = " << setw(4) << oset.psi->i_tec 
	 << endl;
  }
  */
  
  stringstream temp_psp;
  temp_psp << "Final-state parton-level momenta (partonic CMS frame):" << endl;
  for (int i_p = 1; i_p < osi_p_parton[0].size(); i_p++){
    temp_psp << "type_parton[" << setw(2) << 0 << "][" << i_p << "] = " 
	     << setw(3) << osi_type_parton[0][i_p] 
	     << " -> " 
	     << osi_p_parton[0][i_p] << endl; 
  }
  logger << LOG_DEBUG_POINT << temp_psp.str() << endl;
  /*
  temp_psp << endl;
  // not yet filled here !!!
  temp_psp << "Final-state parton-level momenta (hadronic CMS = LAB frame):" << endl;
  for (int i_p = 3; i_p < osi_p_parton[0].size(); i_p++){
    temp_psp << "type_parton[" << setw(2) << 0 << "][" << i_p << "] = " 
	     << setw(3) << osi_type_parton[0][i_p] 
	     << " -> " 
	     << osi_particle_event[0][0][i_p].momentum << endl; 
    // << " m = " << osi_particle_event[0][0][i_p].m << " pT = " << osi_particle_event[0][0][i_p].pT << " eta = " << osi_particle_event[0][0][i_p].eta << endl;
  }
  */

  //  cerr << temp_psp.str() << endl;
  
  //  logger.newLine(LOG_DEBUG_POINT);
  //      }


  if (osi_switch_testcut == 0){
    for (int i_a = 0; i_a < osi_cut_ps.size(); i_a++){
      // should not be needed, as set to correct value later:
      osi_cut_ps[i_a] = 0;

      logger << LOG_DEBUG_VERBOSE << "before: osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << endl;
      perform_event_selection_phasespace(i_a, oset);
      ///
      if (osi_switch_output_cutinfo){logger << LOG_DEBUG_VERBOSE << "[" << setw(2) << i_a << "]" << "   after generic event selection, before user-defined particles   (osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << ")." << endl;}

      if (osi_cut_ps[i_a] == -1){continue;}

      if (oset.user.particle_name.size() > 0){generic.particles(i_a, oset);}
      ///
      if (osi_switch_output_cutinfo){logger << LOG_DEBUG << "[" << setw(2) << i_a << "]" << "   after user-defined particles, before individual cuts   (osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << ")." << endl;}

      if (osi_cut_ps[i_a] == -1){continue;}

      // application of general cuts set via fiducialcut interface:
      logger << LOG_DEBUG_POINT << "osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << "   before fiducialcut." << endl;
      for (int i_f = 0; i_f < oset.esi.fiducial_cut.size(); i_f++){
	oset.esi.fiducial_cut[i_f].apply_fiducialcut(i_a);
	if (osi_cut_ps[i_a] == -1){continue;}
      }

      //      if (osi_cut_ps[i_a] == -1){continue;}

      if (osi_switch_output_cutinfo){logger << LOG_DEBUG << "[" << setw(2) << i_a << "]" << "   after generic particles, before individual event selection   (osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << ")." << endl;}
      generic.cuts(i_a, oset);
      ///
      if (osi_switch_output_cutinfo){logger << LOG_DEBUG << "[" << setw(2) << i_a << "]" << "   after individual cuts   (osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << ")." << endl;}
      logger << LOG_DEBUG_VERBOSE << "after:  osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << endl;
    }
  }
  else if (osi_switch_testcut == 1){
    generic.cuts_test(0, oset);
  }


  for (int i_a = 0; i_a < osi_cut_ps.size(); i_a++){
    logger << LOG_DEBUG_POINT << "osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << endl;
  }
  
}



void perform_event_selection_phasespace(int i_a, observable_set & oset){
  static Logger logger("perform_event_selection_phasespace");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  static int initialization = 1;

  logger << LOG_DEBUG_VERBOSE << "initialization = " << initialization << endl;
  static int n_original = 0;

  // temporary N-jettiness implementation
  //  static int Njettiness = 0;
  //  static vector<fourvector> NJ_axes;
  // -> observable.set.cpp
  
  if (initialization == 1){
    // re-extract the number of particles that enter the event selection without modifications
    logger << LOG_DEBUG_VERBOSE << "initialization = " << initialization << endl;
    vector<int> temp_original;
    for (int i_g = 0; i_g < osi_runtime_original.size(); i_g++){
      logger << LOG_DEBUG_VERBOSE << "i_g = " << i_g << endl;
      int j_g = osi_runtime_original[i_g];
      for (int i_o = 0; i_o < osi_runtime_order[j_g].size(); i_o++){
        int i_p = osi_runtime_order[j_g][i_o];
        int flag = 0;
        for (int i = 0; i < temp_original.size(); i++){if (i_p == temp_original[i]){flag = 1; break;}}
        if (flag == 0){temp_original.push_back(i_p);}
      }
    }
    n_original = temp_original.size();
    for (int i = 0; i < temp_original.size(); i++){logger << LOG_DEBUG_VERBOSE << "temp_original[" << i << "] = " << temp_original[i] << endl;}
    logger << LOG_DEBUG_VERBOSE << "n_original = " << n_original << endl;

    /* -> observable.set.cpp
    // temporary N-jettiness implementation
    if (osi_switch_qTcut == 3){
      Njettiness = oset.order_alphas_born;
      /*
      if (osi_type_contribution == "RT"){Njettiness = oset.order_alphas_born - 1;}
      else if (osi_type_contribution == "L2RT"){Njettiness = oset.order_alphas_born - 3;} // not generic !!!
      else if (osi_type_contribution == "RRA"){Njettiness = oset.order_alphas_born - 2;}
      else if (osi_type_contribution == "RCA"){Njettiness = oset.order_alphas_born - 1;}
      else if (osi_type_contribution == "RVA"){Njettiness = oset.order_alphas_born + 1 - 1;} // not generic !!!
    *//*
     
      NJ_n_axes.resize(3 + Njettiness); // not generic !!! there should be n_jet_born or so... !!!
      NJ_n_axes[1] = fourvector(1., 0., 0., 1.);
      NJ_n_axes[2] = fourvector(1., 0., 0., -1.);
    }
  */
    
    initialization = 0;
  }

  
  if (osi_type_contribution == "RRA" ||
      osi_type_contribution == "RCA" ||
      osi_type_contribution == "RVA" ||
      osi_type_contribution == "RT" ||
      osi_type_contribution == "L2RT"){
    //    fourvector Q;
    /*
    int n_born_particles;
    if (osi_type_contribution == "RRA" && i_a == 0) {
      n_born_particles = osi_p_parton[i_a].size()-5;
    } else {
      n_born_particles = osi_p_parton[i_a].size()-4;
    }

    if (n_born_particles != oset.csi->n_particle_born){logger << LOG_FATAL << "n_born_particles = " << n_born_particles << " != " << oset.csi->n_particle_born << " = oset.csi->n_particle_born" << endl; exit(1);}

	
    for (int i_p = 3; i_p < 3+n_born_particles; i_p++){Q = Q + osi_p_parton[i_a][i_p];}
    */
    oset.QT_Q = fourvector();
    for (int i_p = 3; i_p < 3 + oset.csi->n_particle_born; i_p++){oset.QT_Q = oset.QT_Q + osi_p_parton[i_a][i_p];}
    oset.QT_QT = oset.QT_Q.pT();
    if (osi_switch_qTcut == 1){
      oset.QT_sqrtQ2 = oset.QT_Q.m();
      if (oset.QT_QT / oset.QT_sqrtQ2 < osi_min_qTcut / 100.){osi_cut_ps[i_a] = -1; return;}
      else {
	if (oset.binning_qTcut == "linear"){
	  osi_cut_ps[i_a] = GSL_MIN_INT(int(((oset.QT_QT / oset.QT_sqrtQ2 - osi_min_qTcut / 100.) * 100) / osi_step_qTcut), osi_n_qTcut - 1);
	  /*
	  // temporary check:
	  double temp_no_cut = 0;
	  double temp_value_cut = oset.QT_QT / oset.QT_sqrtQ2 * 100;
	  for (int i_q = osi_n_qTcut - 1; i_q >= 0; i_q--){
	    if (temp_value_cut > oset.value_qTcut[i_q]){temp_no_cut = i_q; break;}
	  }
	  if (osi_cut_ps[i_a] != temp_no_cut){logger << LOG_ERROR << "osi_cut_ps[i_a] = " << osi_cut_ps[i_a] << " != " << temp_no_cut << " = temp_no_cut" << endl;}
	  */
	}
	else {
	  double temp_value_cut = oset.QT_QT / oset.QT_sqrtQ2 * 100;
	  for (int i_q = osi_n_qTcut - 1; i_q >= 0; i_q--){
	    if (temp_value_cut > oset.value_qTcut[i_q]){osi_cut_ps[i_a] = i_q; break;}
	  }
	}
	  
      }
    }
    else if (osi_switch_qTcut == 2){
      if (oset.QT_QT < osi_min_qTcut){osi_cut_ps[i_a] = -1; return;}
      else {
	if (oset.binning_qTcut == "linear"){
	  osi_cut_ps[i_a] = GSL_MIN_INT(int((oset.QT_QT - osi_min_qTcut) / osi_step_qTcut), osi_n_qTcut - 1);
	//        osi_cut_ps[i_a] = GSL_MIN_INT(int((QTpartonic - osi_min_qTcut) / osi_step_qTcut), osi_n_qTcut);
	}
	else {
	  double temp_value_cut = oset.QT_QT;
	  for (int i_q = osi_n_qTcut - 1; i_q >= 0; i_q--){
	    if (temp_value_cut > oset.value_qTcut[i_q]){osi_cut_ps[i_a] = i_q; break;}
	  }
	}
      }
    }
    /*
    for (int i_p = 3; i_p < 3+oset.csi->n_particle_born; i_p++){Q = Q + osi_p_parton[i_a][i_p];}
    double QTpartonic = Q.pT();
    if (osi_switch_qTcut == 1){
      double Qpartonic = Q.m();
      if (QTpartonic / Qpartonic < osi_min_qTcut / 100.){osi_cut_ps[i_a] = -1; return;}
      else {
        osi_cut_ps[i_a] = GSL_MIN_INT(int(((QTpartonic / Qpartonic - osi_min_qTcut / 100.) * 100) / osi_step_qTcut), osi_n_qTcut - 1);
      }
    }
    else if (osi_switch_qTcut == 2){
      if (QTpartonic < osi_min_qTcut){osi_cut_ps[i_a] = -1; return;}
      else {
        osi_cut_ps[i_a] = GSL_MIN_INT(int((QTpartonic - osi_min_qTcut) / osi_step_qTcut), osi_n_qTcut - 1);
	//        osi_cut_ps[i_a] = GSL_MIN_INT(int((QTpartonic - osi_min_qTcut) / osi_step_qTcut), osi_n_qTcut);
      }
    }
    */
    
  }



  for (int i_p = 1; i_p < 3; i_p++){osi_particle_event[0][i_a][i_p] = particle(osi_p_parton[i_a][i_p].zboost(-osi_boost));}
  
  for (int i_p = 3; i_p < osi_p_parton[i_a].size(); i_p++){osi_particle_event[0][i_a][i_p] = particle(osi_p_parton[i_a][i_p].zboost(-osi_boost));}
  for (int i_p = 3; i_p < osi_p_parton[i_a].size(); i_p++){logger << LOG_DEBUG_VERBOSE << "osi_particle_event[0][" << i_a << "][" << i_p << "] = " << osi_particle_event[0][i_a][i_p] << endl;}

  



  

  logger << LOG_DEBUG_VERBOSE << "before: osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << endl;

  //  static int switch_info_cut = 1;
  static stringstream info_cut;
  if (osi_switch_output_cutinfo){
    info_cut.clear();
    info_cut.str("");
    //  info_cut << "phasespace i_a = " << i_a << endl;
    info_cut << endl;
    info_cut << "[" << setw(2) << i_a << "]   generic event selection" << endl;
    info_cut << "[" << setw(2) << i_a << "]   parton-level momenta:" << endl;
    for (int i_p = 1; i_p < osi_p_parton[i_a].size(); i_p++){
      info_cut << "[" << setw(2) << i_a << "]   type_parton[" << i_p << "] = " << setw(3) << osi_type_parton[i_a][i_p] << " -> " << osi_particle_event[0][i_a][i_p].momentum << endl;//<< " pT = " << osi_particle_event[0][i_a][i_p].pT << " eta = " << osi_particle_event[0][i_a][i_p].eta << endl;
    }
    info_cut << endl;
    //  logger << LOG_DEBUG_VERBOSE << info_cut.str();
  }
  
  for (int i_o = 1; i_o < osi_esi.object_list_selection.size(); i_o++){
    //  for (int i_o = 1; i_o < osi_relevant_object_list.size(); i_o++){
    int j_o = osi_equivalent_no_object[i_o]; 
    //    logger << LOG_INFO << "i_o = " << i_o << "   j_o = " << j_o << "   osi_equivalent_no_object.size() = " << osi_equivalent_no_object.size() << endl;
     // i_o -> number in relevant_object_list
    // j_o -> number in object_list
    //    int j_o = osi_no_relevant_object[osi_relevant_object_list[i_o]]; // i_o -> number in relevant_object_list; j_o -> i_o
    osi_particle_event[j_o][i_a].clear();
    osi_n_object[j_o][i_a] = 0;
    //    logger << LOG_INFO << "osi_n_object[" << j_o << "][" << i_a << "] = " << osi_n_object[j_o][i_a] << "   osi_particle_event[" << j_o << "][" << i_a << "].size() = " << osi_particle_event[j_o][i_a].size() << endl;
    //    logger << LOG_INFO << "osi_n_object[" << j_o << "][" << i_a << "] = " << osi_n_object[j_o][i_a] << endl;
    logger << LOG_DEBUG_VERBOSE << "osi_particle_event[" << j_o << "][" << i_a << "].size() = " << osi_particle_event[j_o][i_a].size() << endl;
    //        logger << LOG_INFO << "osi_n_object[" << j_o << "][" << i_a << "] = " << osi_n_object[j_o][i_a] << "   osi_particle_event[" << j_o << "][" << i_a << "].size() = " << osi_particle_event[j_o][i_a].size() << endl;
  }
  
  logger << LOG_DEBUG_VERBOSE << "photon_recombination" << endl;
  logger << LOG_DEBUG_VERBOSE << "osi_photon_recombination = " << osi_photon_recombination << endl;

  ///  vector<int> no_unrecombined_photon;
  static vector<int> no_unrecombined_photon;
  no_unrecombined_photon.clear();

  //  for (int i_p = 3; i_p < osi_p_parton[i_a].size(); i_p++){logger << LOG_DEBUG_VERBOSE << "[" << setw(2) << i_a << "]   type_parton[" << i_p << "] = " << setw(3) << osi_type_parton[i_a][i_p] << " -> " << osi_particle_event[0][i_a][i_p].momentum << endl;


  ////////////////////////////
  //  photon recombination  //
  ////////////////////////////

  //  Recombinable particles are defined in  photon_recombination_list  (defined via  photon_recombination_selection  and  photon_recombination_disable ).
  //  If a particle is dressed with a photon, the combined momentum is stored at the position of its naked momentum, and the momentum of the recombined photon is empty.
  if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   photon_recombination: " << "   photon_recombination = " << osi_photon_recombination << endl;}

  if (osi_photon_recombination){
    photon_recombination(no_unrecombined_photon, i_a, oset);
    /*
    for (int i_c = 0; i_c < osi_ps_runtime_photon_recombination[i_a].size(); i_c++){
      //      logger << LOG_INFO << "osi_ps_runtime_photon_recombination[" << i_a << "][" << i_c << "] = " << osi_ps_runtime_photon_recombination[i_a][i_c] << "   pT2 = " << osi_particle_event[0][i_a][osi_ps_runtime_photon_recombination[i_a][i_c]].pT2  << endl;
      if (osi_particle_event[0][i_a][osi_ps_runtime_photon_recombination[i_a][i_c]].pT2 == 0.){
	logger << LOG_INFO << "i_a = " << i_a << "   Lepton " << osi_ps_runtime_photon_recombination[i_a][i_c] << " is removed in recombination procedure." << "   pT2 = " << osi_particle_event[0][i_a][osi_ps_runtime_photon_recombination[i_a][i_c]].pT2 << endl;
	
      }
    }
    */
  }
  else {no_unrecombined_photon = osi_ps_runtime_photon[i_a];}
  //  Afterwards,  no_unrecombined_photon  contains positions of un-recombined photons in  osi_particle_event[0][i_a]  (all photons if no dressing is done).

  if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   no_unrecombined_photon.size() = " << no_unrecombined_photon.size() << endl;}
  if (osi_switch_output_cutinfo){
    if (no_unrecombined_photon.size() < osi_ps_runtime_photon[i_a].size()){
      for (int i_p = 3; i_p < osi_p_parton[i_a].size(); i_p++){
	info_cut << "[" << setw(2) << i_a << "]   type_parton[" << i_p << "] = " << setw(3) << osi_type_parton[i_a][i_p] << " -> " << osi_particle_event[0][i_a][i_p].momentum << endl;//<< " pT = " << osi_particle_event[0][i_a][i_p].pT << " eta = " << osi_particle_event[0][i_a][i_p].eta << endl;
      }
    }
  }
  logger << LOG_DEBUG_VERBOSE << "no_unrecombined_photon.size() = " << no_unrecombined_photon.size() << endl;



  ////////////////////////
  //  runtime original  //
  ////////////////////////

  //  runtime_original   actually needed only once because "original" particles don't change on reduced phase spaces -> check evaluation of runtime_original !!!
  //  runtime_original   collects partons that enter event selection with their parton-level momenta (after photon dressing if applicable).
  if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   runtime_original: " << "   osi_runtime_original.size() = " << osi_runtime_original.size() << endl;}

  logger << LOG_DEBUG_VERBOSE << "runtime_original" << "   osi_runtime_original.size() = " << osi_runtime_original.size() << endl;
  for (int i_g = 0; i_g < osi_runtime_original.size(); i_g++){
    // j_g -> running number (directing to number in relevant_object_list) of objects contained in runtime_original
    int j_g = osi_runtime_original[i_g];  // j_g -> number in relevant_object_list
    int k_g = osi_equivalent_no_object[j_g];  // k_g -> number in object_list
    // later...   int k_g = j_g;
    
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   osi_esi.object_list_selection[runtime_original[" << i_g << "] = " << j_g << "] = " << osi_esi.object_list_selection[j_g] << endl;}

    logger << LOG_DEBUG_VERBOSE << "osi_runtime_original: i_g = " << i_g << "   j_g = " << j_g << " (" << osi_esi.object_list_selection[j_g] << ")   k_g = " << k_g << " (" << osi_esi.object_list[k_g] << ")" << endl;
    logger << LOG_DEBUG_VERBOSE << "osi_n_object[" << k_g << "][" << i_a << "] = " << osi_n_object[k_g][i_a] << endl;
    
    for (int i_o = 0; i_o < osi_runtime_order[j_g].size(); i_o++){
      int i_p = osi_runtime_order[j_g][i_o];  // i_p -> position in  osi_particle_event[0][i_a]  of parton belong to the respective object class
      if (osi_switch_output_cutinfo){info_cut << "runtime_original: " << i_g << "   osi_particle_event[0][" << i_a << "][" << i_p << "] = " << osi_particle_event[0][i_a][i_p] << endl;}
      if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << setprecision(5)
					      << "   y: " << setw(10) << abs(osi_particle_event[0][i_a][i_p].rapidity) << " < " << setw(8) << osi_esi.pds[j_g].define_y 
					      << "   eta:  " << setw(10) << abs(osi_particle_event[0][i_a][i_p].eta) << " < " << setw(8) << osi_esi.pds[j_g].define_eta 
					      << "   pT:  " << setw(10) << abs(osi_particle_event[0][i_a][i_p].pT) << " > " << setw(8) << osi_esi.pds[j_g].define_pT 
					      << "   ET:  " << setw(10) << abs(osi_particle_event[0][i_a][i_p].ET) << " > " << setw(8) << osi_esi.pds[j_g].define_ET 
					      << "   ---   " 
					      << (abs(osi_particle_event[0][i_a][i_p].rapidity) < osi_esi.pds[j_g].define_y && 
						  abs(osi_particle_event[0][i_a][i_p].eta) < osi_esi.pds[j_g].define_eta && 
						  osi_particle_event[0][i_a][i_p].pT >= osi_esi.pds[j_g].define_pT && 
						  osi_particle_event[0][i_a][i_p].ET >= osi_esi.pds[j_g].define_ET) << endl;}
      
      if (abs(osi_particle_event[0][i_a][i_p].rapidity) < osi_esi.pds[j_g].define_y && 
	  abs(osi_particle_event[0][i_a][i_p].eta) < osi_esi.pds[j_g].define_eta && 
	  osi_particle_event[0][i_a][i_p].pT >= osi_esi.pds[j_g].define_pT && 
	  osi_particle_event[0][i_a][i_p].ET >= osi_esi.pds[j_g].define_ET){
	osi_n_object[k_g][i_a]++;
	osi_particle_event[k_g][i_a].push_back(osi_particle_event[0][i_a][i_p]);
      }
    }
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   " << osi_esi.pds[j_g].n_observed_min << " <= osi_n_object[" << setw(2) <<k_g << "][" << setw(2) <<i_a << "] = " << osi_n_object[k_g][i_a] << " <= " << osi_esi.pds[j_g].n_observed_max << endl;}
	
    if ((osi_n_object[k_g][i_a] < osi_esi.pds[j_g].n_observed_min) || (osi_n_object[k_g][i_a] > osi_esi.pds[j_g].n_observed_max)){
      osi_cut_ps[i_a] = -1; 
      if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   cut after runtime_original" << endl;}
      logger << LOG_DEBUG_VERBOSE << "event at [i_a = " << i_a << "] cut after runtime_original" << endl; 
      if (osi_switch_output_cutinfo){logger << LOG_DEBUG << endl << info_cut.str(); }
      return;
    }
    if (osi_n_object[k_g][i_a] > 1){sort(osi_particle_event[k_g][i_a].begin(), osi_particle_event[k_g][i_a].end(), greaterBypT);}
  }

  if (osi_switch_output_cutinfo){
    info_cut << endl;
    info_cut << "Summary after runtime_original:" << endl;
    for (int i_g = 0; i_g < osi_runtime_original.size(); i_g++){
      // j_g -> running number (directing to number in relevant_object_list) of objects contained in runtime_original
      int j_g = osi_runtime_original[i_g];  // j_g -> number in relevant_object_list
      int k_g = osi_equivalent_no_object[j_g];  // k_g -> number in object_list
      //      info_cut << "[" << setw(2) << i_a << "]   osi_esi.object_list_selection[runtime_original[" << i_g << "] = " << j_g << "] = " << osi_esi.object_list_selection[j_g] << endl;
      info_cut << "[" << setw(2) << i_a << "]   osi_particle_event[" << setw(2) << k_g << " -> " << setw(10) << osi_esi.object_list_selection[j_g] << "][" << setw(2) << i_a << "].size() = " << osi_particle_event[k_g][i_a].size() << endl;
    }
    info_cut << endl;
  }


  ///////////////////////
  //  runtime missing  //
  ///////////////////////

  if (oset.ps_runtime_missing[i_a].size() > 0 ){
    //      runtime_missing   actually needed only once because "invisible" particles don't change on reduced phase spaces -> check evaluation of runtime_missing !!!
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   runtime_missing: " << "   oset.ps_runtime_missing[" << i_a << "].size() = " << oset.ps_runtime_missing[i_a].size() << endl;}
    //    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   runtime_missing: " << "   osi_runtime_missing.size() = " << osi_runtime_missing.size() << endl;}
    logger << LOG_DEBUG_VERBOSE << "runtime_missing" << "   osi_runtime_missing.size() = " << osi_runtime_missing.size() << endl;
    //    logger << LOG_INFO << "runtime_missing" << "   osi_runtime_missing.size() = " << osi_runtime_missing.size() << endl;
    //    logger << LOG_INFO << "runtime_missing" << "   oset.ps_runtime_missing[" << i_a << "].size() = " << oset.ps_runtime_missing[i_a].size() << endl;
    
    int j_g = osi_esi.observed_object_selection["missing"];  // j_g -> number in relevant_object_list
    int k_g = osi_esi.observed_object["missing"];  //  k_g -> number in object_list
    // later...   int k_g = j_g;
    //  cout << "j_g = " << j_g << "   k_g = " << k_g << endl;
    fourvector p_missing;
    for (int i_p = 0; i_p < oset.ps_runtime_missing[i_a].size(); i_p++){
      //    cout << "
      //      logger << LOG_INFO << "osi_particle_event[0][" << i_a << "][" << oset.ps_runtime_missing[i_a][i_p] << "].momentum = " << osi_particle_event[0][i_a][oset.ps_runtime_missing[i_a][i_p]].momentum << endl;
      p_missing = p_missing + osi_particle_event[0][i_a][oset.ps_runtime_missing[i_a][i_p]].momentum;
    }
    osi_particle_event[k_g][i_a].push_back(particle(p_missing));
    osi_n_object[k_g][i_a] = 1;
    
    //    logger << LOG_INFO << "osi_particle_event[" << k_g << "][" << i_a << "][0] = " << osi_particle_event[k_g][i_a][0].momentum << "   pT = " << osi_particle_event[k_g][i_a][0].pT << endl;
    
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << setprecision(5)
					    << "   pT:  " << setw(10) << osi_particle_event[k_g][i_a][0].pT << " > " << setw(8) << osi_esi.pds[j_g].define_pT 
					    << "   ET:  " << setw(10) << osi_particle_event[k_g][i_a][0].ET << " > " << setw(8) << osi_esi.pds[j_g].define_ET 
					    << "   ---   " 
					    << (osi_particle_event[k_g][i_a][0].pT < osi_esi.pds[j_g].define_pT ||
						osi_particle_event[k_g][i_a][0].ET < osi_esi.pds[j_g].define_ET) << endl;}
    
    if (osi_particle_event[k_g][i_a][0].pT < osi_esi.pds[j_g].define_pT || osi_particle_event[k_g][i_a][0].ET < osi_esi.pds[j_g].define_ET){
      osi_cut_ps[i_a] = -1; 
      if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   cut after runtime_missing" << endl;}
      logger << LOG_DEBUG_VERBOSE << "event at [i_a = " << i_a << "] cut after runtime_missing" << endl; 
      if (osi_switch_output_cutinfo){logger << LOG_DEBUG << endl << info_cut.str(); }
      return;
    }
    
    if (osi_switch_output_cutinfo){
      info_cut << endl;
      info_cut << "Summary after runtime_missing:" << endl;
      //      for (int i_g = 0; i_g < osi_runtime_missing.size(); i_g++){
      for (int i_g = 0; i_g < oset.ps_runtime_missing[i_a].size(); i_g++){
	info_cut << "[" << setw(2) << i_a << "]   osi_particle_event[" << setw(2) << k_g << " -> " << setw(10) << osi_esi.object_list_selection[j_g] << "][" << setw(2) << i_a << "].size() = " << osi_particle_event[k_g][i_a].size() << endl;
      }
      info_cut << endl;
    }
    //    logger << LOG_INFO << "[" << setw(2) << i_a << "]   osi_particle_event[" << setw(2) << k_g << " -> " << setw(10) << osi_esi.object_list_selection[j_g] << "][" << setw(2) << i_a << "].size() = " << osi_particle_event[k_g][i_a].size() << endl;
  }

  /*
  //      runtime_missing   actually needed only once because "invisible" particles don't change on reduced phase spaces -> check evaluation of runtime_missing !!!
  if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   runtime_missing: " << "   osi_runtime_missing.size() = " << osi_runtime_missing.size() << endl;}
  logger << LOG_DEBUG_VERBOSE << "runtime_missing" << "   osi_runtime_missing.size() = " << osi_runtime_missing.size() << endl;

  for (int i_g = 0; i_g < osi_runtime_missing.size(); i_g++){
    int j_g = osi_runtime_missing[i_g];       // j_g -> number in relevant_object_list
    int k_g = osi_equivalent_no_object[j_g];  // k_g -> number in object_list
    // later...   int k_g = j_g;

    fourvector p_missing;
    for (int i_o = 0; i_o < osi_runtime_order[j_g].size(); i_o++){
      int i_p = osi_runtime_order[j_g][i_o];
      p_missing = p_missing + osi_particle_event[0][i_a][i_p].momentum;
    }
    osi_particle_event[k_g][i_a].push_back(particle(p_missing));
        
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << setprecision(5)
					    << "   pT:  " << setw(10) << osi_particle_event[k_g][i_a][0].pT << " > " << setw(8) << osi_esi.pds[j_g].define_pT 
					    << "   ET:  " << setw(10) << osi_particle_event[k_g][i_a][0].ET << " > " << setw(8) << osi_esi.pds[j_g].define_ET 
					    << "   ---   " 
					    << (osi_particle_event[k_g][i_a][0].pT < osi_esi.pds[j_g].define_pT ||
						osi_particle_event[k_g][i_a][0].ET < osi_esi.pds[j_g].define_ET) << endl;}
    
    
    if (osi_particle_event[k_g][i_a][0].pT < osi_esi.pds[j_g].define_pT || osi_particle_event[k_g][i_a][0].ET < osi_esi.pds[j_g].define_ET){
      osi_cut_ps[i_a] = -1; 
      if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   cut after runtime_missing" << endl;}
      logger << LOG_DEBUG_VERBOSE << "event at [i_a = " << i_a << "] cut after runtime_missing" << endl; 
      if (osi_switch_output_cutinfo){logger << LOG_DEBUG << endl << info_cut.str(); }
      return;
    }
  }
      
  if (osi_switch_output_cutinfo){
    info_cut << endl;
    info_cut << "Summary after runtime_missing:" << endl;
    for (int i_g = 0; i_g < osi_runtime_missing.size(); i_g++){
      // j_g -> running number (directing to number in relevant_object_list) of objects contained in runtime_missing
      int j_g = osi_runtime_missing[i_g];  // j_g -> number in relevant_object_list
      int k_g = osi_equivalent_no_object[j_g];  // k_g -> number in object_list
      //      info_cut << "[" << setw(2) << i_a << "]   osi_esi.object_list_selection[runtime_missing[" << i_g << "] = " << j_g << "] = " << osi_esi.object_list_selection[j_g] << endl;
      info_cut << "[" << setw(2) << i_a << "]   osi_particle_event[" << setw(2) << k_g << " -> " << setw(10) << osi_esi.object_list_selection[j_g] << "][" << setw(2) << i_a << "].size() = " << osi_particle_event[k_g][i_a].size() << endl;
    }
    info_cut << endl;
  }

  */
      
  //////////////////////////////////
  //  determination of protojets  //
  //////////////////////////////////

  // particle_protojet  does not contain photons here - added later after Frixione isolation.
  if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   jet_algorithm = " << osi_jet_algorithm << endl;}
  logger << LOG_DEBUG_VERBOSE << "jet algorithm: select protojet candidates" << endl;

  vector<particle> particle_protojet;
  vector<vector<int> > protojet_flavour;
  vector<vector<int> > protojet_parton_origin;

  logger << LOG_DEBUG_VERBOSE << "osi_ps_runtime_jet_algorithm.size() = " << osi_ps_runtime_jet_algorithm.size() << endl;
  for (int i_l = 0; i_l < osi_ps_runtime_jet_algorithm[i_a].size(); i_l++){
    logger << LOG_DEBUG_VERBOSE << "osi_ps_runtime_jet_algorithm[" << i_a << "].size() = " << osi_ps_runtime_jet_algorithm[i_a].size() << endl;
    int i_p = osi_ps_runtime_jet_algorithm[i_a][i_l];
    particle_protojet.push_back(osi_particle_event[0][i_a][i_p]);

    protojet_parton_origin.push_back(vector<int> (1, i_p));
    logger << LOG_DEBUG_VERBOSE << "osi_type_parton[" << i_a << "].size() = " << osi_type_parton[i_a].size() << endl;
    if (abs(osi_type_parton[i_a][i_p]) < 5){protojet_flavour.push_back(vector<int> (1, 0));}  // treatment of light jets -> flavour numbers between -4 and +4  ->  0 (ljet)
    else {protojet_flavour.push_back(vector<int> (1, osi_type_parton[i_a][i_p]));}
    logger << LOG_DEBUG_VERBOSE << "osi_ps_runtime_jet_algorithm[" << i_a << "].size() = " << osi_ps_runtime_jet_algorithm[i_a].size() << endl;
  }
  //  particle_protojet  contains all particles relevant for Frixione isolation.

  for (int i_p = 0; i_p < particle_protojet.size(); i_p++){
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   particle_protojet[" << i_p << "]   pT = " << particle_protojet[i_p].pT << "   ET = " << particle_protojet[i_p].ET << "   flavour = " << protojet_flavour[i_p][0] << "   origin = " << protojet_parton_origin[i_p][0] << endl;}
  }
  logger << LOG_DEBUG_VERBOSE << "phasespace " << i_a << " " << particle_protojet.size() << ": protojets" << endl;
  for (int i_p = 0; i_p < particle_protojet.size(); i_p++){
    logger << LOG_DEBUG_VERBOSE << "particle_protojet[" << i_p << "] = " << particle_protojet[i_p] << endl << "   flavour = " << protojet_flavour[i_p][0] << "   origin = " << protojet_parton_origin[i_p][0] << endl;
  }
  


  //////////////////////////
  //  Frixione isolation  //
  //////////////////////////

  if (osi_frixione_isolation){
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   photon_isolation" << endl;}
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   photon_isolation:   runtime_photon_isolation.size() = " << osi_runtime_photon_isolation.size() << endl;}
    logger << LOG_DEBUG_VERBOSE << "runtime_photon_isolation" << endl;
    logger << LOG_DEBUG_VERBOSE << "osi_runtime_photon_isolation.size() = " << osi_runtime_photon_isolation.size() << endl;
    
    vector<particle> particle_protojet_photon;
    vector<vector<int> > protojet_flavour_photon;
    vector<vector<int> > protojet_parton_origin_photon;
    
    for (int i_g = 0; i_g < osi_runtime_photon_isolation.size(); i_g++){
      int j_g = osi_runtime_photon_isolation[i_g];
      int k_g = osi_equivalent_no_object[j_g];
      logger << LOG_DEBUG_VERBOSE << "no_unrecombined_photon.size() = " << no_unrecombined_photon.size() << endl;
      for (int i_u = 0; i_u < no_unrecombined_photon.size(); i_u++){
	int i_p = no_unrecombined_photon[i_u];
	if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << setprecision(5)
						<< "   y: " << setw(10) << abs(osi_particle_event[0][i_a][i_p].rapidity) << " < " << setw(8) << osi_esi.pds[j_g].define_y 
						<< "   eta:  " << setw(10) << abs(osi_particle_event[0][i_a][i_p].eta) << " < " << setw(8) << osi_esi.pds[j_g].define_eta 
						<< "   pT:  " << setw(10) << abs(osi_particle_event[0][i_a][i_p].pT) << " > " << setw(8) << osi_esi.pds[j_g].define_pT 
						<< "   ET:  " << setw(10) << abs(osi_particle_event[0][i_a][i_p].ET) << " > " << setw(8) << osi_esi.pds[j_g].define_ET 
						<< "   ---   " 
						<< (abs(osi_particle_event[0][i_a][i_p].rapidity) < osi_esi.pds[j_g].define_y && 
						    abs(osi_particle_event[0][i_a][i_p].eta) < osi_esi.pds[j_g].define_eta && 
						    osi_particle_event[0][i_a][i_p].pT >= osi_esi.pds[j_g].define_pT && 
						    osi_particle_event[0][i_a][i_p].ET >= osi_esi.pds[j_g].define_ET) << endl;}
	
	if (abs(osi_particle_event[0][i_a][i_p].rapidity) < osi_esi.pds[j_g].define_y && 
	    abs(osi_particle_event[0][i_a][i_p].eta) < osi_esi.pds[j_g].define_eta && 
	    osi_particle_event[0][i_a][i_p].pT >= osi_esi.pds[j_g].define_pT && 
	    osi_particle_event[0][i_a][i_p].ET >= osi_esi.pds[j_g].define_ET){
	  int temp_number_photon = osi_n_object[k_g][i_a];  //  temp_number_photon = number of identified photons before the one considered now (i_p)
	  logger << LOG_DEBUG_VERBOSE << "temp_number_photon = " << temp_number_photon << endl;
	  frixione_isolation(osi_n_object[k_g][i_a], osi_particle_event[k_g][i_a], osi_particle_event[0][i_a][i_p], particle_protojet, oset);
	  // temporary: particle_protojet_photon etc. could be directly filled...
	  if (temp_number_photon == osi_n_object[k_g][i_a] && osi_photon_jet_algorithm){  //  Happens if the considered photon is not identified. Non-identified photons will enter the jet algorithm afterwards.
	    logger << LOG_DEBUG_VERBOSE << "Add un-identified photons to jet candidates" << endl;
	    particle_protojet_photon.push_back(osi_particle_event[0][i_a][i_p]);
	    protojet_parton_origin_photon.push_back(vector<int> (1, i_p));
	    protojet_flavour_photon.push_back(vector<int> (1, 22)); // photon
	    logger << LOG_DEBUG_VERBOSE << "particle_protojet_photon.size() = " << particle_protojet_photon.size() << endl;
	  }
	  logger << LOG_DEBUG_VERBOSE << "osi_particle_event[" << k_g << "][" << i_a << "].size() = " << osi_particle_event[k_g][i_a].size() << endl;
	  logger << LOG_DEBUG_VERBOSE << "osi_n_object[" << k_g << "][" << i_a << "] = " << osi_n_object[k_g][i_a] << endl;
	}
	else if (osi_photon_jet_algorithm){  //  Non-identified photons will enter the jet algorithm afterwards.
	  particle_protojet_photon.push_back(osi_particle_event[0][i_a][i_p]);
	  protojet_parton_origin_photon.push_back(vector<int> (1, i_p));
	  protojet_flavour_photon.push_back(vector<int> (1, 22)); // photon
	  logger << LOG_DEBUG_VERBOSE << "particle_protojet_photon.size() = " << particle_protojet_photon.size() << endl;
	}
      }

      // Not fully generic if photon--photon recombination is applied!
      
      if (!osi_photon_photon_recombination){
	if ((osi_n_object[k_g][i_a] < osi_esi.pds[j_g].n_observed_min) || (osi_n_object[k_g][i_a] > osi_esi.pds[j_g].n_observed_max)){
	  osi_cut_ps[i_a] = -1; 
	  if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   cut after runtime_photon_isolation" << endl;}
	  logger << LOG_DEBUG_VERBOSE << "event at [i_a = " << i_a << "] cut after runtime_photon_isolation" << endl; 
	  if (osi_switch_output_cutinfo){logger << LOG_DEBUG << endl << info_cut.str(); }
	  return;
	}
	// shouldn't this be osi_particle_event[k_g][i_a].size() instead of osi_n_object[k_g][i_a] ???
	//      if (osi_n_object[k_g][i_a] > 1){sort(osi_particle_event[k_g][i_a].begin(), osi_particle_event[k_g][i_a].end(), greaterBypT);}
	// let's try:
	if (osi_particle_event[k_g][i_a].size() > 1){sort(osi_particle_event[k_g][i_a].begin(), osi_particle_event[k_g][i_a].end(), greaterBypT);}
	// Should be identical.
      }
      
      // added to let non-identified photons enter the jet algorithm
      if (osi_photon_jet_algorithm){  //  Non-identified photons are added to protojets.
	for (int i_p = 0; i_p < particle_protojet_photon.size(); i_p++){
	  particle_protojet.push_back(particle_protojet_photon[i_p]);
	  protojet_parton_origin.push_back(protojet_parton_origin_photon[i_p]);
	  protojet_flavour.push_back(protojet_flavour_photon[i_p]);
	}
      }
      //
      else {
	// otherwise add to "photon" particle vector - without increasing photon number ???
      }
    }
    /*
    if (i_a == 0){
      logger << LOG_INFO << "photon_photon_recombination input (isolated photon: " << osi_particle_event[osi_access_object["photon"]][i_a].size() << " - rest: " << particle_protojet_photon.size() << "): " << osi_particle_event[osi_access_object["photon"]][i_a].size() + particle_protojet_photon.size() << endl;
    }
    */
    /*
    if (i_a == 0){
      for (int i_p = 3; i_p < osi_p_parton[i_a].size(); i_p++){
	logger << LOG_INFO << "[" << setw(2) << i_a << "]   type_parton[" << i_p << "] = " << setw(3) << osi_type_parton[i_a][i_p] << " -> " << osi_particle_event[0][i_a][i_p].momentum << endl;//<< " pT = " << osi_particle_event[0][i_a][i_p].pT << " eta = " << osi_particle_event[0][i_a][i_p].eta << endl;
      }


      logger << LOG_INFO << "photon_photon_recombination input (isolated photon: " << osi_particle_event[osi_access_object["photon"]][i_a].size() << " - rest: " << particle_protojet_photon.size() << "):" << endl;
      //      logger << LOG_INFO << "photon.size() = " << osi_particle_event[osi_access_object["photon"]][i_a].size() << endl;
      for (int i_p = 0; i_p < osi_particle_event[osi_access_object["photon"]][i_a].size(); i_p++){
	logger << LOG_INFO << "photon[" << i_p << "] = " << osi_particle_event[osi_access_object["photon"]][i_a][i_p].momentum << endl;
      }
      //      logger << LOG_INFO << "particle_protojet_photon.size() = " << particle_protojet_photon.size() << endl;
      for (int i_p = 0; i_p < particle_protojet_photon.size(); i_p++){
	logger << LOG_INFO << "particle_protojet_photon[" << i_p << "] = " << particle_protojet_photon[i_p].momentum << endl;
      }
      logger.newLine(LOG_INFO);
    }
    */
  }
  else {
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   No Frixione isolation performed" << endl;}
    logger << LOG_DEBUG_VERBOSE << "No Frixione isolation is done." << endl;

    //    logger << LOG_INFO << "i_a = " << i_a << "   no_unrecombined_photon.size() = " << no_unrecombined_photon.size() << endl;

    if (osi_photon_jet_algorithm){  //  All photons (i.e. un-recombined photons) are added to protojets if this option is chosen.
      for (int i_u = 0; i_u < no_unrecombined_photon.size(); i_u++){
	int i_p = no_unrecombined_photon[i_u];
	particle_protojet.push_back(osi_particle_event[0][i_a][i_p]);
	protojet_parton_origin.push_back(vector<int> (1, i_p));
	protojet_flavour.push_back(vector<int> (1, osi_type_parton[i_a][i_p]));
      }
    }
  }

  
  
  //  logger << LOG_INFO << "i_a = " << i_a << "   particle_protojet.size() = " << particle_protojet.size() << endl;

  if (osi_switch_output_cutinfo){
    info_cut << endl;
    info_cut << "Summary after Frixione isolation:" << endl;
    for (int i_g = 0; i_g < osi_runtime_photon_isolation.size(); i_g++){
      // j_g -> running number (directing to number in relevant_object_list) of objects contained in runtime_photon_isolation
      int j_g = osi_runtime_photon_isolation[i_g];  // j_g -> number in relevant_object_list
      int k_g = osi_equivalent_no_object[j_g];  // k_g -> number in object_list
      //      info_cut << "[" << setw(2) << i_a << "]   osi_esi.object_list_selection[runtime_photon_isolation[" << i_g << "] = " << j_g << "] = " << osi_esi.object_list_selection[j_g] << endl;
      info_cut << "[" << setw(2) << i_a << "]   osi_particle_event[" << setw(2) << k_g << " -> " << setw(10) << osi_esi.object_list_selection[j_g] << "][" << setw(2) << i_a << "].size() = " << osi_particle_event[k_g][i_a].size() << endl;
    }
    info_cut << "Protojets after Frixione isolation:" << endl;
    for (int i_p = 0; i_p < particle_protojet.size(); i_p++){
      info_cut << "protojet[" << i_p << "] = " << particle_protojet[i_p].momentum << "   flavour = " << protojet_flavour[i_p][0] << endl;
    }
    info_cut << endl;
  }


  //////////////////////////////////////////////////
  //  elimination of protojets too close to beam  //
  //////////////////////////////////////////////////

  if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   osi_parton_y_max = " << osi_parton_y_max << "   osi_parton_eta_max = " << osi_parton_eta_max << endl;}
  logger << LOG_DEBUG_VERBOSE << "jet algorithm: exclude protojet candidates" << endl;

  // exclude protojets which are too close to the beam (in terms of y or eta)
  if (osi_parton_y_max != 0. || osi_parton_eta_max != 0.){
    for (int j = particle_protojet.size() - 1; j >= 0; j--){
      if (abs(particle_protojet[j].rapidity) > osi_parton_y_max || abs(particle_protojet[j].eta) > osi_parton_eta_max){
	particle_protojet.erase(particle_protojet.begin() + j);
	protojet_flavour.erase(protojet_flavour.begin() + j);
	protojet_parton_origin.erase(protojet_parton_origin.begin() + j);
      }
    }
  }

  //#include "user/specify.extra.exclusion.jet.algorithm.cxx"

  /////////////////////
  //  jet algorithm  //
  /////////////////////

  if (osi_switch_output_cutinfo){
    for (int i_p = 0; i_p < particle_protojet.size(); i_p++){
      info_cut << "[" << setw(2) << i_a << "]" << "   particle_protojet[" << i_p << "]   pT = " << particle_protojet[i_p].pT << "   ET = " << particle_protojet[i_p].ET << "   flavour = " << protojet_flavour[i_p][0] << "   origin = " << protojet_parton_origin[i_p][0] << endl;
    }
  }
      
  logger << LOG_DEBUG_VERBOSE << "jet algorithm: combine protojets to jets" << endl;

  vector<particle> particle_jet;
  vector<vector<int> > jet_flavour;
  vector<vector<int> > jet_parton_origin;
  
  jet_algorithm_flavour(particle_protojet, protojet_flavour, protojet_parton_origin, particle_jet, jet_flavour, jet_parton_origin, oset);

  if (osi_switch_output_cutinfo){
    for (int i_p = 0; i_p < particle_jet.size(); i_p++){
      stringstream jet_flavour_ss;
      jet_flavour_ss << "( ";
      for (int i_q = 0; i_q < jet_flavour[i_p].size(); i_q++){
	jet_flavour_ss << jet_flavour[i_p][i_q];
	if (i_q < jet_flavour[i_p].size() - 1){jet_flavour_ss << " ";}
      }
      jet_flavour_ss << " )";
      
      stringstream jet_parton_origin_ss;
      jet_parton_origin_ss << "( ";
      for (int i_q = 0; i_q < jet_parton_origin[i_p].size(); i_q++){
	jet_parton_origin_ss << jet_parton_origin[i_p][i_q];
	if (i_q < jet_parton_origin[i_p].size() - 1){jet_parton_origin_ss << " ";}
      }
      jet_parton_origin_ss << " )";
      
      info_cut << "[" << setw(2) << i_a << "]" << "   particle_jet[" << i_p << "]   pT = " << particle_jet[i_p].pT << "   ET = " << particle_jet[i_p].ET << "   flavour = " << jet_flavour_ss.str() << "   origin = " << jet_parton_origin_ss.str() << endl;
      //      info_cut << "[" << setw(2) << i_a << "]" << "   particle_jet[" << i_p << "]   pT = " << particle_jet[i_p].pT << "   ET = " << particle_jet[i_p].ET << "   flavour = " << jet_flavour[i_p][0] << "   origin = " << jet_parton_origin[i_p][0] << endl;
    }
  }

  //  logger << LOG_INFO << "i_a = " << i_a << "   particle_jet.size() = " << particle_jet.size() << endl;

  //#include "user/specify.add.excluded.jet.cxx"



  ////////////////////////////////////////////////////////////
  //  Frixione isolation - removal of jets close to photon  //
  ////////////////////////////////////////////////////////////

  // sort out jets which are too close to an 'isolated photon' according to Frixione
  // might need to be adapted for photon_recombination & frixione_isolation !!!
  if (osi_frixione_isolation && osi_frixione_jet_removal){ // up to date ??? osi_equivalent_no_object[osi_no_relevant_object[...]] constructions !!!
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   osi_runtime_photon_isolation.size() = " << osi_runtime_photon_isolation.size() << endl;}
    
    // temporarily uncomment to see effect !!!
    if (osi_runtime_photon_isolation.size() != 0){//osi_no_relevant_object["photon"]
      logger << LOG_DEBUG_VERBOSE << "photon-isolation caused jet removal" << endl;
       
      for (int i_phot = 0; i_phot < osi_n_object[osi_access_object["photon"]][i_a]; i_phot++){
	for (int i_jet = particle_jet.size() - 1; i_jet >= 0; i_jet--){
	  logger << LOG_DEBUG_VERBOSE << "sqrt(R2_eta(osi_particle_event[osi_equivalent_no_object[osi_no_relevant_object[photon]]][i_a][i_phot], particle_jet[i_jet])) = " << sqrt(R2_eta(osi_particle_event[osi_access_object["photon"]][i_a][i_phot], particle_jet[i_jet])) << "   frixione_delta_0 = " << osi_frixione_delta_0 << endl;
	  if (sqrt(R2_eta(osi_particle_event[osi_access_object["photon"]][i_a][i_phot], particle_jet[i_jet])) < osi_frixione_delta_0){
	    logger << LOG_DEBUG_VERBOSE << "particle_jet[" << i_jet << "] = " << particle_jet[i_jet] << "   removed.   frixione_delta_0 = " << osi_frixione_delta_0 << endl;
	    particle_jet.erase(particle_jet.begin() + i_jet);
	    jet_flavour.erase(jet_flavour.begin() + i_jet);
	  }
	}
      }
    }
    
    for (int i_p = 0; i_p < particle_jet.size(); i_p++){
      if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   particle_jet[" << i_p << "]   pT = " << particle_jet[i_p].pT << "   ET = " << particle_jet[i_p].ET << "   flavour = " << jet_flavour[i_p][0] << "   origin = " << jet_parton_origin[i_p][0] << endl;}
    }
  }
  

  /////////////////////////////////////////////////
  //  jet algorithm - interpretation of outcome  //
  /////////////////////////////////////////////////

  logger << LOG_DEBUG_VERBOSE << "particle_jet_id" << endl;
  //  vector<particle_id> particle_jet_id(particle_jet.size());
  static vector<particle_id> particle_jet_id;
  particle_jet_id.clear();
  particle_jet_id.resize(particle_jet.size());

  for (int i_j = 0; i_j < particle_jet.size(); i_j++){
    particle_jet_id[i_j].xparticle = &particle_jet[i_j];
    particle_jet_id[i_j].content = jet_flavour[i_j];

    // IR-safe defintion - b and bx are combined to light jet: switch needed !!!
    particle_jet_id[i_j].identity = accumulate(particle_jet_id[i_j].content.begin(), particle_jet_id[i_j].content.end(), 0);
    if (oset.msi.M_b > 0. && particle_jet_id[i_j].identity == 0){
      // no IR-safe defintion in case of massless b quarks - b and bx are combined to bottom jet: switch needed !!!
      for (int i_q = 0; i_q < particle_jet_id[i_j].content.size(); i_q++){
	if (abs(particle_jet_id[i_j].content[i_q]) == 5){particle_jet_id[i_j].identity = 5; break;}
      }
    }
    
    // modification compared to IR-safe version of flavoured jet combination:
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   particle_jet_id[" << i_j << "].content.size() = " << particle_jet_id[i_j].content.size() << endl;}
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   particle_jet_id[" << i_j << "].identity = " << particle_jet_id[i_j].identity << endl;}
	
    // valid only for W+jets !!! does, however, not look wrong for other processes with photon recombination !!!
    // with this (E_phot -- E_had) criterion, it could happen that non-isolated photon collect enough pT to be considered a photon,
    // which, however, at present has no impact because the (properly isolated) photons are counted right after the Frixione isolation...
    // For processes without photons in the jet algorithm, nothing should happen here !!!
    
    // deactivate: this could allow to treat photons as jets:
    //    if (particle_jet_id[i_j].content.size() == 1){continue;}
	
    if (particle_jet_id[i_j].identity == 0){} // light jet or b-bx jet
    else if (particle_jet_id[i_j].identity == 5 || // bjet_b
	     particle_jet_id[i_j].identity == -5){continue;} // bjet_bx
    else {
      double E_phot = 0.;
      double E_had = 0.;
      int count_phot = 0;
      for (int i_l = 0; i_l < particle_jet_id[i_j].content.size(); i_l++){
	if (particle_jet_id[i_j].content[i_l] == 22){
	  E_phot += (osi_particle_event[0][i_a][jet_parton_origin[i_j][i_l]].momentum).x0();
	  count_phot++;
	}
	else {
	  E_had += (osi_particle_event[0][i_a][jet_parton_origin[i_j][i_l]].momentum).x0();
	}
      }

      if (osi_switch_output_cutinfo){
	info_cut << "E_phot = " << E_phot << endl;
	info_cut << "E_had = " << E_had << endl;
      }

      if (E_phot / (E_phot + E_had) >= osi_photon_E_threshold_ratio){
	//     if (E_phot / (E_phot + E_had) > osi_photon_E_threshold_ratio){
	// bug !!!      if (E_phot / (E_had + E_had) > osi_photon_E_threshold_ratio){
	particle_jet_id[i_j].identity = 22;
      }
      else {
	particle_jet_id[i_j].identity = particle_jet_id[i_j].identity - count_phot * 22;
      }
      if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   particle_jet_id[" << i_j << "].identity = " << particle_jet_id[i_j].identity << endl;}
      logger << LOG_DEBUG_VERBOSE << "E_phot = " << E_phot << endl;
      logger << LOG_DEBUG_VERBOSE << "E_had = " << E_had << endl;
      logger << LOG_DEBUG_VERBOSE << "particle_jet_id[" << i_j << "].identity = " << particle_jet_id[i_j].identity << endl;
    }

    //#include "user/specify.jet.recombination.cxx"
  }
      
  sort(particle_jet_id.begin(), particle_jet_id.end(), idgreaterBypT);

  

  if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   runtime_jet_recombination:   osi_runtime_jet_recombination.size() = " << osi_runtime_jet_recombination.size() << endl;}
      
  logger << LOG_DEBUG_VERBOSE << "runtime_jet_recombination" << endl;
  for (int i_g = 0; i_g < osi_runtime_jet_recombination.size(); i_g++){
    int j_g = osi_runtime_jet_recombination[i_g];
    int k_g = osi_equivalent_no_object[j_g];
    logger << LOG_DEBUG_VERBOSE << "runtime_jet_recombination: i_g = " << i_g << "   j_g = " << j_g << "   k_g = " << k_g << endl;
    //    logger << LOG_DEBUG_VERBOSE << "osi_relevant_define_y[" << j_g << "] = " << osi_esi.pds[j_g].define_y << endl;
    //    logger << LOG_DEBUG_VERBOSE << "osi_relevant_define_eta[" << j_g << "] = " << osi_esi.pds[j_g].define_eta << endl;
    //    logger << LOG_DEBUG_VERBOSE << "osi_relevant_define_pT[" << j_g << "] = " << osi_esi.pds[j_g].define_pT << endl;
    //    logger << LOG_DEBUG_VERBOSE << "osi_relevant_define_ET[" << j_g << "] = " << osi_esi.pds[j_g].define_ET << endl;
    if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]   osi_esi.object_list_selection[runtime_jet_recombination[" << i_g << "] = " << j_g << "] = " << osi_esi.object_list_selection[j_g] << endl;}
    
    logger << LOG_DEBUG_VERBOSE << "particle_jet_id.size() = " << particle_jet_id.size() << endl;
    for (int i_j = 0; i_j < particle_jet_id.size(); i_j++){
      logger << LOG_DEBUG_VERBOSE << "(*particle_jet_id[" << i_j << "].xparticle) = " << (*particle_jet_id[i_j].xparticle) << endl;
      /*
	logger << LOG_DEBUG_VERBOSE << "osi_no_relevant_object[photon] = " << osi_no_relevant_object["photon"] << endl;
	logger << LOG_DEBUG_VERBOSE << "k_g = " << k_g << endl;
	logger << LOG_DEBUG_VERBOSE << "j_g = " << j_g << endl;
	logger << LOG_DEBUG_VERBOSE << "particle_jet_id[" << i_j << "].identity = " <<particle_jet_id[i_j].identity  << endl;
      */
      // k_g (relevant_object) -> j_g (object)
      if      (j_g == osi_no_relevant_object["jet"] && particle_jet_id[i_j].identity == 22){continue;}
      else if (j_g == osi_no_relevant_object["photon"] && particle_jet_id[i_j].identity != 22){continue;}
      else if (j_g == osi_no_relevant_object["ljet"] && particle_jet_id[i_j].identity != 0){continue;}
      else if (j_g == osi_no_relevant_object["bjet"] && abs(particle_jet_id[i_j].identity) != 5){continue;}
      else if (j_g == osi_no_relevant_object["bjet_b"] && particle_jet_id[i_j].identity != +5){continue;}
      else if (j_g == osi_no_relevant_object["bjet_bx"] && particle_jet_id[i_j].identity != -5){continue;}
      
      if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << setprecision(5)
					      << "   y: " << setw(10) << abs((*particle_jet_id[i_j].xparticle).rapidity) << " < " << setw(8) << osi_esi.pds[j_g].define_y 
					      << "   eta:  " << setw(10) << abs((*particle_jet_id[i_j].xparticle).eta) << " < " << setw(8) << osi_esi.pds[j_g].define_eta 
					      << "   pT:  " << setw(10) << abs((*particle_jet_id[i_j].xparticle).pT) << " > " << setw(8) << osi_esi.pds[j_g].define_pT 
					      << "   ET:  " << setw(10) << abs((*particle_jet_id[i_j].xparticle).ET) << " > " << setw(8) << osi_esi.pds[j_g].define_ET 
					      << "   ---   " 
					      << (abs((*particle_jet_id[i_j].xparticle).rapidity) < osi_esi.pds[j_g].define_y && 
						  abs((*particle_jet_id[i_j].xparticle).eta) < osi_esi.pds[j_g].define_eta)
					      << "   ---   " 
					      << ((*particle_jet_id[i_j].xparticle).pT >= osi_esi.pds[j_g].define_pT && 
						  (*particle_jet_id[i_j].xparticle).ET >= osi_esi.pds[j_g].define_ET) << endl;}
      
      
      if (abs((*particle_jet_id[i_j].xparticle).rapidity) < osi_esi.pds[j_g].define_y && 
	  abs((*particle_jet_id[i_j].xparticle).eta) < osi_esi.pds[j_g].define_eta){ 
	if ((*particle_jet_id[i_j].xparticle).pT >= osi_esi.pds[j_g].define_pT &&
	    (*particle_jet_id[i_j].xparticle).ET >= osi_esi.pds[j_g].define_ET){
	  // not counted as an identified photon if Frixione isolation is used !!!
	  if (j_g == osi_no_relevant_object["photon"] && osi_frixione_isolation){}
	  else {
	    osi_n_object[k_g][i_a]++;
	  }
	  logger << LOG_DEBUG_VERBOSE << "osi_n_object[" << k_g << "][" << i_a << "] = " << osi_n_object[k_g][i_a] << endl;
	}
	osi_particle_event[k_g][i_a].push_back(*particle_jet_id[i_j].xparticle);
	logger << LOG_DEBUG_VERBOSE << "osi_particle_event[" << k_g << "][" << i_a << "].size() = " << osi_particle_event[k_g][i_a].size() << endl;
      }
    }
    logger << LOG_DEBUG_VERBOSE << "osi_esi.pds[" << j_g << "].n_observed_min = " << osi_esi.pds[j_g].n_observed_min << endl;
    logger << LOG_DEBUG_VERBOSE << "osi_esi.pds[" << j_g << "].n_observed_max = " << osi_esi.pds[j_g].n_observed_max << endl;

    if (!(j_g == osi_no_relevant_object["photon"] && osi_frixione_isolation && osi_photon_photon_recombination)){
      if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   " << osi_esi.pds[j_g].n_observed_min << " <= osi_n_object[" << setw(2) <<k_g << "][" << setw(2) <<i_a << "] = " << osi_n_object[k_g][i_a] << " <= " << osi_esi.pds[j_g].n_observed_max << endl;}
      
      if ((osi_n_object[k_g][i_a] < osi_esi.pds[j_g].n_observed_min) || (osi_n_object[k_g][i_a] > osi_esi.pds[j_g].n_observed_max)){
	osi_cut_ps[i_a] = -1; 
	if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   cut after runtime_jet_recombination" << endl;}
	logger << LOG_DEBUG_VERBOSE << "event at [i_a = " << i_a << "] cut after runtime_jet_recombination" << endl; 
      if (osi_switch_output_cutinfo){logger << LOG_DEBUG << endl << info_cut.str(); }
      return;
      }
    }
  }


    // Only applied for two photons: (combine with osi_photon_E_threshold_ratio = 1.!)
  static int photon_photon_counter = 0;
  if (osi_frixione_isolation && osi_photon_photon_recombination){
    int j_g = osi_no_relevant_object["photon"];
    int k_g = osi_access_object["photon"];
    int stop_photon_photon_recombination = 0;
    if (osi_particle_event[k_g][i_a].size() < 2){stop_photon_photon_recombination = 1;}
    
    int change_photon_photon_recombination = 0;
    //    if (k_g != osi_access_object["photon"]){logger << LOG_DEBUG << "k_g = " << k_g << " != " << osi_access_object["photon"] << " = osi_access_object[photon]" << endl;}
    
    while (!stop_photon_photon_recombination){
      double distance = osi_photon_photon_recombination_R;
      int i_rec = -1;
      int j_rec = -1;
      for (int i_p = 0; i_p < osi_particle_event[k_g][i_a].size(); i_p++){
	for (int j_p = i_p + 1; j_p < osi_particle_event[k_g][i_a].size(); j_p++){
	  double temp_distance = 0.;
	  if (osi_photon_R_definition == 0){temp_distance = R2_eta(osi_particle_event[k_g][i_a][i_p], osi_particle_event[k_g][i_a][j_p]);}
	  else if (osi_photon_R_definition == 1){temp_distance = R2_rapidity(osi_particle_event[k_g][i_a][i_p], osi_particle_event[k_g][i_a][j_p]);}
	  //	  logger << LOG_INFO << "photon_photon_distance = " << temp_distance << endl;
	  if (temp_distance < distance){
	    //	    logger << LOG_INFO << "photon_photon_distance = " << temp_distance << " < " << distance << endl;
	    distance = temp_distance;
	    i_rec = i_p;
	    j_rec = j_p;	      
	  }
	}
      }
      if (i_rec != -1 && j_rec != -1){
	change_photon_photon_recombination = 1;
	logger << LOG_DEBUG << "photon_photon_recombination applied (" << ++photon_photon_counter << "): distance = " << distance << endl;
	logger << LOG_DEBUG << "osi_particle_event[" << k_g << "][" << i_a << "][" << i_rec << "] = " << osi_particle_event[k_g][i_a][i_rec].momentum << endl;
	logger << LOG_DEBUG << "osi_particle_event[" << k_g << "][" << i_a << "][" << j_rec << "] = " << osi_particle_event[k_g][i_a][j_rec].momentum << endl;
	osi_particle_event[k_g][i_a][i_rec] = osi_particle_event[k_g][i_a][i_rec] + osi_particle_event[k_g][i_a][j_rec];
	osi_particle_event[k_g][i_a].erase(osi_particle_event[k_g][i_a].begin() + j_rec);
	if (j_rec < osi_n_object[k_g][i_a]){osi_n_object[k_g][i_a]--;}
      }
      else {
	stop_photon_photon_recombination = 1;
	// end while loop (to be added) !!!
      }
    }
    //      int old_n_photon = osi_n_object[k_g][i_a];
    if (change_photon_photon_recombination){
      osi_n_object[k_g][i_a] = 0;
      for (int i_p = 0; i_p < osi_particle_event[k_g][i_a].size(); i_p++){
	if (abs(osi_particle_event[0][i_a][i_p].rapidity) < osi_esi.pds[j_g].define_y && 
	    abs(osi_particle_event[0][i_a][i_p].eta) < osi_esi.pds[j_g].define_eta && 
	    osi_particle_event[0][i_a][i_p].pT >= osi_esi.pds[j_g].define_pT && 
	    osi_particle_event[0][i_a][i_p].ET >= osi_esi.pds[j_g].define_ET){
	  osi_n_object[k_g][i_a]++;
	}
      }
    }
    
    if ((osi_n_object[k_g][i_a] < osi_esi.pds[j_g].n_observed_min) || (osi_n_object[k_g][i_a] > osi_esi.pds[j_g].n_observed_max)){
      osi_cut_ps[i_a] = -1; 
      if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   cut after photon_photon_recombination" << endl;}
      logger << LOG_DEBUG_VERBOSE << "event at [i_a = " << i_a << "] cut after photon_photon_recombination" << endl; 
      if (osi_switch_output_cutinfo){logger << LOG_DEBUG << endl << info_cut.str(); }
      return;
    }
  }

    

  
  if (osi_switch_output_cutinfo){
    info_cut << endl;
    info_cut << "Summary after runtime_jet_recombination:" << endl;
    for (int i_g = 0; i_g < osi_runtime_jet_recombination.size(); i_g++){
      // j_g -> running number (directing to number in relevant_object_list) of objects contained in runtime_jet_recombination
      int j_g = osi_runtime_jet_recombination[i_g];  // j_g -> number in relevant_object_list
      int k_g = osi_equivalent_no_object[j_g];  // k_g -> number in object_list
      //      info_cut << "[" << setw(2) << i_a << "]   osi_esi.object_list_selection[runtime_jet_recombination[" << i_g << "] = " << j_g << "] = " << osi_esi.object_list_selection[j_g] << endl;
      info_cut << "[" << setw(2) << i_a << "]   osi_particle_event[" << setw(2) << k_g << " -> " << setw(10) << osi_esi.object_list_selection[j_g] << "][" << setw(2) << i_a << "].size() = " << osi_particle_event[k_g][i_a].size() << endl;
    }
    info_cut << endl;
  }

  ///


  // Old N-jettiness implementation shifted from here !!!

  


  if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   generic event selection passed   (osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << ")." << endl;}
  //      if (osi_switch_output_cutinfo){logger << LOG_DEBUG << endl << info_cut.str(); }
  

  logger << LOG_DEBUG_VERBOSE << "after:  osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << endl;
  logger << LOG_DEBUG_VERBOSE << "osi_switch_qTcut = " << osi_switch_qTcut << endl;
  logger << LOG_DEBUG_VERBOSE << "osi_output_n_qTcut = " << osi_output_n_qTcut << endl;
  if (osi_switch_qTcut == 0){osi_cut_ps[i_a] = 0;}
  else if (osi_switch_qTcut == 1 && osi_output_n_qTcut == 1){osi_cut_ps[i_a] = 0;}
  else if (osi_switch_qTcut == 3 && osi_output_n_qTcut == 1){osi_cut_ps[i_a] = 0;}
  else if (osi_type_contribution == "CT" || 
	   osi_type_contribution == "CJ" || 
	   osi_type_contribution == "CT2" || 
	   osi_type_contribution == "CJ2" || 
	   osi_type_contribution == "L2CT" || 
	   osi_type_contribution == "L2CJ"){
    // Contributions with a Born phase-space, but qTcut dependence:
    osi_cut_ps[i_a] = osi_n_qTcut - 1;  // check if '- 1' is correct !!!
  }
  else {
    logger << LOG_DEBUG_VERBOSE << "osi_switch_qTcut = " << osi_switch_qTcut << endl;
    logger << LOG_DEBUG_VERBOSE << "osi_output_n_qTcut = " << osi_output_n_qTcut << endl;
    logger << LOG_DEBUG_VERBOSE << "osi_n_qTcut - 1 = " << osi_n_qTcut - 1 << endl;
    // osi_cut_ps[i_a] should remain as it is -> No need to re-evaluate !!!
    //    osi_cut_ps[i_a] = osi_n_qTcut - 1;  // check if '- 1' is correct !!!
  }
  

  if (osi_switch_output_cutinfo){info_cut << "[" << setw(2) << i_a << "]" << "   generic event selection passed and osi_cut_ps[" << i_a << "] set   (osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << ")." << endl;}
  if (osi_switch_output_cutinfo){logger << LOG_DEBUG << endl << info_cut.str();}
  
  /*
    for (int i_o = 1; i_o < osi_esi.object_list_selection.size(); i_o++){// i_o -> number in relevant_object_list
    int j_o = osi_equivalent_no_object[i_o];  // j_o -> number in object_list
    //    int j_o = osi_no_relevant_object[osi_esi.object_list_selection[i_o]]; // i_o -> number in relevant_object_list; j_o -> i_o
    logger << LOG_DEBUG_VERBOSE << "osi_esi.object_list_selection[" << i_o << "] = " << setw(10) << osi_esi.object_list_selection[i_o] << "   osi_n_object[" << j_o << "][" << i_a << "] = " << osi_n_object[j_o][i_a] << "   osi_particle_event[" << j_o << "][" << i_a << "].size() = " << osi_particle_event[j_o][i_a].size() << endl;
    for (int i_p = 0; i_p < osi_particle_event[j_o][i_a].size(); i_p++){
    logger << LOG_DEBUG_VERBOSE << "          osi_particle_event[" << j_o << "][" << i_a << "][" << i_p << "] = " << osi_particle_event[j_o][i_a][i_p].momentum << endl; 
    }
    }
  */

  // N-jettiness implementation (typically uses results from jet algorithm, thus located here.)
  // changed to simulate NJ implementation in QT environment:
  if (osi_switch_qTcut == 3){
    /*
    if (osi_type_contribution == "RRA" ||
	osi_type_contribution == "RCA" ||
	osi_type_contribution == "RVA" ||
	osi_type_contribution == "RT" ||
	osi_type_contribution == "L2RT"){
    */
      
  if (osi_type_contribution == "RRJ" ||
      osi_type_contribution == "RCJ" ||
      osi_type_contribution == "RVJ" ||
      osi_type_contribution == "RJ" ||
      osi_type_contribution == "L2RJ"){
      
      //      logger.newLine(LOG_DEBUG);
      for (int i_p = 1; i_p < osi_p_parton[i_a].size(); i_p++){
	logger << LOG_DEBUG << "p(hadronic frame)[" << i_p << "] = " << osi_particle_event[0][i_a][i_p].momentum << endl;
      }
      
      oset.Njettiness_calculate_NJ_axes(i_a);

      logger.newLine(LOG_DEBUG);
      for (int i_j = 1; i_j < oset.NJ_q_axes.size(); i_j++){
	logger << LOG_DEBUG << "NJ_q_axes[" << i_j << "] = " << oset.NJ_q_axes[i_j] << endl;
	//	logger << LOG_DEBUG << "       [" << i_j << "] = " << oset.NJ_Ei[i_j] * oset.NJ_n_axes[i_j] << endl;
      }
      logger.newLine(LOG_DEBUG);
      for (int i_j = 1; i_j < oset.NJ_n_axes.size(); i_j++){
	logger << LOG_DEBUG << "NJ_n_axes[" << i_j << "] = " << oset.NJ_n_axes[i_j] << endl;
      }

      oset.Njettiness_calculate_NJ_axes_frame(i_a);

      // new definition in respective rest frame:
      double tau_jettiness = 0.;
      for (int i_p = 0; i_p < osi_ps_runtime_jet_algorithm[i_a].size(); i_p++){
	vector<double> delta_tau(3 + oset.csi->n_jet_born);
	for (int i_q = 1; i_q < oset.NJ_q_axes.size(); i_q++){
	  delta_tau[i_q] = abs(2 * oset.NJ_q_axes_frame[i_q] * oset.NJ_p_parton_frame[i_p] / oset.NJ_Qi[i_q]);
	}
	tau_jettiness += *min_element(delta_tau.begin() + 1, delta_tau.end());
      }

      // determination of cut_ps[i_a] for taucut scan:
      if (tau_jettiness < osi_min_qTcut){osi_cut_ps[i_a] = -1; return;}
      else {
	if (oset.binning_qTcut == "linear"){
	  osi_cut_ps[i_a] = GSL_MIN_INT(int((tau_jettiness - osi_min_qTcut) / osi_step_qTcut), osi_n_qTcut - 1);
	}
	else {
	  double temp_value_cut = tau_jettiness;
	  for (int i_q = osi_n_qTcut - 1; i_q >= 0; i_q--){
	    if (temp_value_cut > oset.value_qTcut[i_q]){osi_cut_ps[i_a] = i_q; break;}
	  }
	}
      }
    }
  }  


    
  // temporary N-jettiness implementation (still in the qT-subtraction framework, including contribution names etc.)
  if (osi_switch_qTcut == -3){// to deactivate this for the moment !!!
    //  if (osi_switch_qTcut == 3){
    if (osi_type_contribution == "RRA" ||
	osi_type_contribution == "RCA" ||
	osi_type_contribution == "RVA" ||
	osi_type_contribution == "RT" ||
	osi_type_contribution == "L2RT"){
      /*
      NJ_n_axes[1] = fourvector(1., 0., 0., 1.);
      NJ_n_axes[2] = fourvector(1., 0., 0., -1.);
      */
      double tau_jettiness = 0.;
      for (int i_j = 0; i_j < oset.csi->n_jet_born; i_j++){
	fourvector temp_jet = osi_particle_event[osi_access_object["jet"]][i_a][i_j].momentum;
	double temp_jet_R = temp_jet.r();
	oset.NJ_n_axes[3 + i_j] = fourvector(1., temp_jet.x1() / temp_jet_R, temp_jet.x2() / temp_jet_R, temp_jet.x3() / temp_jet_R);
      }
      // all available QCD partons
      stringstream temp;
      for (int i_p = 0; i_p < osi_ps_runtime_jet_algorithm[i_a].size(); i_p++){	temp << setw(5) << osi_ps_runtime_jet_algorithm[i_a][i_p];}
      //      logger << LOG_INFO << "oset.csi->n_jet_born = " << oset.csi->n_jet_born << "   osi_ps_runtime_jet_algorithm[" << i_a << "].size() = " << osi_ps_runtime_jet_algorithm[i_a].size() << " --- " << temp.str() << endl;
      
      // calculate N-jettiness tau
      for (int i_p = 0; i_p < osi_ps_runtime_jet_algorithm[i_a].size(); i_p++){
	vector<double> delta_tau(3 + oset.csi->n_jet_born);
	for (int i_j = 1; i_j < oset.NJ_n_axes.size(); i_j++){
	  //osi_particle_event[0][i_a][i_p]
	  delta_tau[i_j] = abs(oset.NJ_n_axes[i_j] * osi_particle_event[0][i_a][osi_ps_runtime_jet_algorithm[i_a][i_p]].momentum);
	  //	  delta_tau[i_j] = oset.NJ_n_axes[i_j] * osi_p_parton[i_a][osi_ps_runtime_jet_algorithm[i_a][i_p]];
	  /*
	  if (delta_tau[i_j] < 0.){
	    logger << LOG_INFO << "oset.NJ_n_axes[" << i_j << "] * osi_particle_event[0][" << i_a << "][" << osi_ps_runtime_jet_algorithm[i_a][i_p] << "] is negative: delta_tau[" << i_j << "] = " << delta_tau[i_j] << endl;
	    delta_tau[i_j] = abs(delta_tau[i_j]);
	  }
	  */
	}
	//	logger << LOG_INFO << "i_p = " << i_p << "   tau_jettiness = " << tau_jettiness << endl;
	/*
	double min = 1.e99;
	for (int i_j = 0; i_j < oset.NJ_n_axes.size(); i_j++){
	  if (delta_tau[i_j] < min){min = delta_tau[i_j];}
	}
	logger << LOG_INFO << "min = " << min << " =?= " << *min_element(delta_tau.begin(), delta_tau.end()) << " = *min_element(delta_tau.begin(), delta_tau.end())" << endl;
	*/

	tau_jettiness += *min_element(delta_tau.begin() + 1, delta_tau.end());
	//	tau_jettiness = tau_jettiness + *min_element(delta_tau.begin(), delta_tau.end());
	//	logger << LOG_INFO << "i_p = " << i_p << "   tau_jettiness = " << tau_jettiness << endl;
      }

      for (int i_j = 1; i_j < oset.NJ_n_axes.size(); i_j++){
	logger << LOG_DEBUG << "oset.NJ_n_axes[" << i_j << "] = " << oset.NJ_n_axes[i_j] << endl;//"   delta_tau[" << i_j << "] = " << delta_tau[i_j] << endl;
      }
      logger << LOG_DEBUG << "oset.binning_qTcut = " << oset.binning_qTcut << endl;
     
      //      logger << LOG_INFO << "oset.csi->n_jet_born = " << oset.csi->n_jet_born << "   tau_jettiness = " << tau_jettiness << "   osi_min_qTcut = " << osi_min_qTcut << endl;
      if (tau_jettiness < osi_min_qTcut){osi_cut_ps[i_a] = -1; return;}
      else {
	if (oset.binning_qTcut == "linear"){
	  osi_cut_ps[i_a] = GSL_MIN_INT(int((tau_jettiness - osi_min_qTcut) / osi_step_qTcut), osi_n_qTcut - 1);
	}
	else {
	  double temp_value_cut = tau_jettiness;
	  for (int i_q = osi_n_qTcut - 1; i_q >= 0; i_q--){
	    if (temp_value_cut > oset.value_qTcut[i_q]){osi_cut_ps[i_a] = i_q; break;}
	  }
	}
      ///	osi_cut_ps[i_a] = GSL_MIN_INT(int((tau_jettiness - osi_min_qTcut) / osi_step_qTcut), osi_n_qTcut - 1);
      }
      logger << LOG_DEBUG << "osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << endl;

      //      logger << LOG_INFO << "osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << endl;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "after:  osi_cut_ps[" << i_a << "] = " << osi_cut_ps[i_a] << endl;
      
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;

}
