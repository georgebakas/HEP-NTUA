#include "../include/classes.cxx"


/*
double observable_set::Njettiness_calculate_NJ_axes_assigned_energy(fourvector & temp_jet){
  double axis_energy = 1.;
  if (switch_NJcut_axes_energy == 1){
    //    axis_energy = temp_jet.x0() / temp_jet.r();
    axis_energy = temp_jet.x0();
  }    
  else if (switch_NJcut_axes_energy == 2){
    //    axis_energy = 1.;
    axis_energy = temp_jet.r();
  }
  else if (switch_NJcut_axes_energy == 3){
    //    axis_energy = .5 * (temp_jet.x0() + temp_jet.r()) / temp_jet.r();
    axis_energy = .5 * (temp_jet.x0() + temp_jet.r());
  }
  return axis_energy;
}
*/

/*
  void observable_set::Njettiness_calculate_NJ_axes_normalized(int i_a){
  Logger logger("observable_set::Njettiness_calculate_NJ_axes_normalized");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // NJ_axes_normalized are always defined in the hadronic frame with two options:
  // switch_NJcut_axes == 1:   axes defined via jet algorithm (according to input)
  // switch_NJcut_axes == 2:   axes defined via partitioning using 1-jettiness

  NJ_axes_normalized.resize(3 + csi->n_jet_born);
  // NJ_axes (1, 2) -> light-like momenta of the beam directions (normalized)
  // NJ_axes (3, ... 2 + n) -> light-like momenta, normalized through deviding the 3-momentum of the jet by its r() (-> independent of possible jet mass)
  NJ_axes_normalized[1] = fourvector(1., 0., 0., 1.);
  NJ_axes_normalized[2] = fourvector(1., 0., 0., -1.);

  if (switch_NJcut_axes == 1){
    for (int i_j = 0; i_j < csi->n_jet_born; i_j++){
      // use first n_jet_born results from the applied jet algorithm here for now...
      // particle_event -> hadronic centre-of-mass system
      // particle_event[access_object["jet"]][i_a] -> jet momenta in hadronic centre-of-mass system
      fourvector temp_jet = particle_event[access_object["jet"]][i_a][i_j].momentum;
      double temp_jet_R = temp_jet.r();
      NJ_axes_normalized[3 + i_j] = fourvector(1., temp_jet.x1() / temp_jet_R, temp_jet.x2() / temp_jet_R, temp_jet.x3() / temp_jet_R);
    }
  }
  else if (switch_NJcut_axes == 2){
   // 
    int temp_order = ps_runtime_jet_algorithm[i_a].size() - csi->n_jet_born;
    if (temp_order == 1){// NLO
      if (csi->n_jet_born == 1){
	vector<double> tau_1_NLO(3, 0.);
	vector<fourvector> parton(ps_runtime_jet_algorithm[i_a].size());
	for (int i_p = 0; i_p < ps_runtime_jet_algorithm[i_a].size(); i_p++){
	  parton[i_p] = particle_event[0][i_a][ps_runtime_jet_algorithm[i_a][i_p]].momentum;
	}
	tau_1_NLO[0] = parton[0].r() - abs(parton[0].x3());
	tau_1_NLO[1] = parton[1].r() - abs(parton[1].x3());
	tau_1_NLO[2] = parton[0].r() + parton[1].r() - (parton[1] + parton[2]).r();
	int min_index = -1;
	double min = 1.e99;
	for (int i_x = 0; i_x < tau_1_NLO.size(); i_x++){
	  if (tau_1_NLO[i_x] < min){min = tau_1_NLO[i_x]; min_index = i_x;}
	}

	fourvector temp_jet;
	if (min_index == 0 || min_index == 1){temp_jet = parton[min_index];}
	else {temp_jet = parton[0] + parton[1];}
	
	double temp_jet_R = temp_jet.r();
	NJ_axes_normalized[3] = fourvector(1., temp_jet.x1() / temp_jet_R, temp_jet.x2() / temp_jet_R, temp_jet.x3() / temp_jet_R);
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
*/

/*
void observable_set::Njettiness_calculate_NJ_axes_energy(int i_a){
  Logger logger("observable_set::Njettiness_calculate_NJ_axes_energy");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // NJ_axes_energy are always defined in the hadronic frame with three options:
  // switch_NJcut_axes_energy == 1:   E_i = p_i^0
  // switch_NJcut_axes_energy == 2:   E_i = |p_i|
  // switch_NJcut_axes_energy == 3:   E_i = (p_i^0 + |p_i|) / 2
  
  NJ_axes_energy.resize(3 + csi->n_jet_born);

  NJ_axes_energy[1] = x_pdf[1] * E_beam;
  NJ_axes_energy[2] = x_pdf[2] * E_beam;
  // not generic !!! Axes need to be defined first !!!
  for (int i_j = 0; i_j < csi->n_jet_born; i_j++){
    if (switch_NJcut_axes_energy == 1){
      NJ_axes_energy[3 + i_j] = particle_event[access_object["jet"]][i_a][i_j].momentum.x0();
    }    
    else if (switch_NJcut_axes_energy == 2){
      NJ_axes_energy[3 + i_j] = particle_event[access_object["jet"]][i_a][i_j].momentum.r();
    }
    else if (switch_NJcut_axes_energy == 3){
      NJ_axes_energy[3 + i_j] = .5 * (particle_event[access_object["jet"]][i_a][i_j].momentum.x0() + particle_event[access_object["jet"]][i_a][i_j].momentum.r());
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/

/*
void observable_set::Njettiness_calculate_NJ_axes_Qi(int i_a){
  Logger logger("observable_set::Njettiness_calculate_NJ_axes_Qi");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_j = 1; i_j < NJ_Qi.size(); i_j++){
    if (switch_NJcut_measure == 1){
      // hadronic frame
      NJ_Qi[i_j] = 2 * NJ_axes_energy[i_j];
    }    
    else if (switch_NJcut_measure == 2){
      // Born frame
      NJ_Qi[3 + i_j] = 2 * particle_event[access_object["jet"]][i_a][i_j].momentum.r();
    }
    else if (switch_NJcut_measure == 3){
       // 1-jettiness-axis frame
     NJ_Qi[3 + i_j] = 2 * .5 * (particle_event[access_object["jet"]][i_a][i_j].momentum.x0() + particle_event[access_object["jet"]][i_a][i_j].momentum.r());
    }
    else if (switch_NJcut_measure == 4){
      // leptonic frame
    }
  }
 
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/



