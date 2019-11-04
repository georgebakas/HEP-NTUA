#include "../include/classes.cxx"
xdistribution::xdistribution(){
  xdistribution_name = "";
  xdistribution_type = "";
  all_particle.push_back(vector<int> (0));
  all_particle_group.push_back(vector<int> (0));
  required_subparticle.push_back(vector<int> (0));
  required_particle.push_back(0);
  n_bins = 0;
  start = 0.;
  end = 0.;
  step = 0.;
  string edges = "";
  symm = 0;
}

void xdistribution::determineBin() {
  Logger logger("xdistribution::determineBin");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  observable.resize(oset->n_ps);
  observable_max.resize(oset->n_ps);
  bin.resize(oset->n_ps);
  bin_max.resize(oset->n_ps);
  
  for (int i_a = 0; i_a < oset->n_ps; i_a++) {
    if (oset->cut_ps[i_a] == -1) {
      // special cases should only be needed for backward compatibility; can be removed once old implementation has been removed
      if (xdistribution_type!="XS" && xdistribution_type!="multiplicity") {
//         bin[i_a]=0;
        continue;
      }
    }
    
    // make sure that all required particles exist, otherwise don't bin
    if (!checkRequirements(i_a)) {
      bin[i_a] = -1;
      bin_max[i_a] = -1;
      
      if (xdistribution_type == "pTveto" || xdistribution_type == "pTvetoexcl") {
        // nothing needs to be vetoed!
        bin[i_a]=0;
        bin_max[i_a] = n_bins;
      }
      
      continue;
    }
    
    fillMomenta(i_a);
    int save_RA_x_a = oset->RA_x_a;
    oset->RA_x_a = i_a;
    computeObservable(observable[i_a], observable_max[i_a]);
    oset->RA_x_a = save_RA_x_a;

    // why double ???
    double _bin = -1;
    double _bin_max = -1;


    if (typeCumulative == CUMULATIVE_NONE){
      _bin = performBinning(observable[i_a]);
      if (_bin < 0){_bin = -1;}
      else if (_bin >= n_bins){_bin = -1;}
    }
    else if (typeCumulative == CUMULATIVE_LOWER){
      _bin = 1 + performBinning(observable[i_a]);
      if (_bin >= n_bins){_bin = -1;}
      if (_bin >= 0){_bin_max = n_bins;}
      else {_bin = -1; _bin_max = -1;}
    }
    else if (typeCumulative == CUMULATIVE_UPPER) {
      _bin_max = 1+performBinning(observable_max[i_a]);
      _bin_max = min(_bin_max,double(n_bins));
      if (_bin_max >= 0) {
        _bin = 0;
      } else {
        _bin = -1;
        _bin_max = -1;
      }
    } else if (typeCumulative == CUMULATIVE_BOTH) {
      _bin = 1+performBinning(observable[i_a]);
      _bin_max = 1+performBinning(observable_max[i_a]);
      if (_bin >= n_bins) {
        _bin = -1;
      }
      if (_bin >= 0 && _bin_max >= 0) {
        _bin_max = min(_bin_max,double(n_bins));
      } else {
        _bin = -1;
        _bin_max = -1;
      }
    }
    
    bin[i_a] = int(_bin);
    bin_max[i_a] = int(_bin_max);
      
    if (bin[i_a]<-1) {
      cout << bin[i_a] << ", " << _bin << endl;
    }

     //cout << xdistribution_name << ", " << i_a << ": " << bin[i_a] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

bool xdistribution::checkRequirements(int i_a) {
  // special cases
  if (xdistribution_type == "multiplicity")
    return true;
  for (int group=0; group<particles.size(); group++) {
    for (int particle=0; particle<particles[group].size(); particle++) {
      if (particles[group][particle].required_definition) {
        if (particles[group][particle].number > oset->n_object[particles[group][particle].type][i_a]) {
          return false;
        }
      } else if (particles[group][particle].required_exists) {
        if (particles[group][particle].number > oset->particle_event[particles[group][particle].type][i_a].size()) {
          return false;
        }
      }
    }
  }
  
  return true;
}

void xdistribution::fillMomenta(int i_a) {
  for (int group=0; group<particles.size(); group++) {
    reconstructedParticles[group] = nullvector;
    for (int particle=0; particle<particles[group].size(); particle++) {
      if (particles[group][particle].number <= oset->particle_event[particles[group][particle].type][i_a].size() && particles[group][particle].number>=1) {
        particles[group][particle].momentum = oset->particle_event[particles[group][particle].type][i_a][particles[group][particle].number-1].momentum;
        particles[group][particle].is_present = true;
      } else {
        particles[group][particle].momentum = fourvector(0.0,0.0,0.0,0.0);
        particles[group][particle].is_present = false;
      }
      particles[group][particle].multiplicity = oset->n_object[particles[group][particle].type][i_a];
      
      reconstructedParticles[group] = reconstructedParticles[group] + particles[group][particle].momentum;
    }
  }
}

double xdistribution::performBinning(double value) {
   //cout << xdistribution_name << ": startpoint=" << start << ", val=" << value << ", endpoint=" << end << endl;

  double bin = -1;

  if (type_binning == "linear"){
    bin = (value - start) / step;
  }
  else if (type_binning == "logarithmic"){
    bin = (log10(value) - start) / step;
  }
  else if (type_binning == "irregular"){
    if (value < bin_edge[0] || value > bin_edge[n_bins]){
      return -1;
    }
    for (int i_b = 1; i_b < n_bins + 1; i_b++){
      //       cout << "val=" << value << ", edge=" << edges[i_b] << endl;
      if (value < bin_edge[i_b]){
        bin = i_b - 1;
        break;
      }
    }
    //    return -1;
  }
  else {
    assert(false);
    return -1;
  }

  if (munich_isnan(bin) || munich_isinf(bin)) {
    cout << xdistribution_name << ": bin is inf/nan, bin=" << bin << ", value=" << value << endl;
    bin = -1;
  }

  /*
  // check for absurdly large bin values
  if (xdistribution_name == "absetamaxvetopTlep20_lep"){
    cout << xdistribution_name << "   bin = " << bin << endl;
    cout << xdistribution_name << "   (double)std::numeric_limits<int>::max() = " << (double)std::numeric_limits<int>::max() << endl;
    cout << xdistribution_name << "   (double)std::numeric_limits<int>::min() = " << (double)std::numeric_limits<int>::min() << endl;
    //    cout << xdistribution_name << "   bin = " << bin << endl;
  }
  if (bin > (double)std::numeric_limits<int>::max() || bin < (double)std::numeric_limits<int>::min()) {
    //    logger << LOG_DEBUG_VERBOSE << xdistribution_name << ": bin too large, bin=" << bin << ", value=" << value << endl;
    bin = -1;
  }
  */
  
  return bin;
}

void xdistribution::computeObservable(double &observable, double &observable_max) {
  Logger logger("xdistribution::computeObservable");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

//  else if (xdistribution_type == "ymax"){xdistribution_number =  8;}
//  else if (xdistribution_type == "ymin"){xdistribution_number =  9;}

 
  observable = 0.0;
  observable_max = 0.0;
  
  // To be added somwhere:
  // 'reconstructedParticle[i-1] == particle i' in distrubtion.dat
  // .eta() is class function to compute pseudo-rapidity

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //////////////////////////////////////////////////////////
  //  Short Tuturial for creating your own distributions  //
  //////////////////////////////////////////////////////////
  //
  // There can be several groups 'particle i' (i=1,2,3,...) defined in distribution.dat
  // The behavior when several particle groups are defined depends on the observable
  //
  // There can also be several particles in one group, eg:
  // particle 1 = lep 1  # hardest lepton
  // particle 1 = lep 2  # second-hardest lepton
  // The momentum of group 'particle 1' is then computed as vectorial sum of the 4-momenta:
  // p(particle 1) = p(lep 1) + p(lep 2)
  //
  // The reconstructed (vectorial) sum of each particle group is accessible below by:
  // reconstructedParticles[group]
  // where
  // reconstructedParticles[group] == 'particle group+1' in distribution.dat, ie:
  // reconstructedParticles[0] == 'particle 1'
  // reconstructedParticles[1] == 'particle 2'
  // reconstructedParticles[2] == 'particle 3'
  // ... 
  // 
  //
  // 4-momenta (fourvector class)
  // ----------------------------
  //
  // The C++ type of 'reconstructedParticles' is 'fourvector', ie a 4-momentum.
  // 4-momenta can be summed (use 'fourvector myfourvector' to define a new fourvector), eg:
  // fourvector sum
  // fourvector electron = reconstructedParticles[0]
  // fourvector muon     = reconstructedParticles[1]
  // sum = electron + muon
  //
  // or direclty:
  // fourvector sum
  // sum = reconstructedParticles[0] + reconstructedParticles[1]
  // 
  //
  // Define observables
  // ------------------
  //
  // Some observables are already predefined through functions in the class and can be directly accessed, eg:
  // reconstructedParticles[0].pT()   gives the transverse momentum pT
  // reconstructedParticles[0].eta()  gives the pseudo-rapidity eta
  // ...
  // see the class definition (src-MUNICH/classes/fourvector.cpp) and header (src-MUNICH/classes/header/fourvector.h)
  // for more further pre-defined observables.
  //
  // ALL (other) observables can be computed directly from the 4-momenta:
  // E = reconstructedParticles[0].x0()
  // x = reconstructedParticles[0].x1()
  // y = reconstructedParticles[0].x2()
  // z = reconstructedParticles[0].x3()
  //
  //
  //
  // Define new distribution
  // -----------------------
  //
  // Follow the pre-defined distributions below and define your own 'xdistribution_type'. Copy a an 'else if(...){...}'
  // block, choose a unique name of your distribition and define your 'observable', eg:
  //
  // ...
  // else if (xdistribution_type == "unique_distribution_name") {
  //   assert(particles.size()==1 && "ERROR in xdistribution_type unique_distribution_name: only one particle group allowed");
  //   E = reconstructedParticles[0].x0()
  //   x = reconstructedParticles[0].x1()
  //   y = reconstructedParticles[0].x2()
  //   z = reconstructedParticles[0].x3()
  //   observable = E+x+y+z
  // }
  // ...
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  if (xdistribution_type == "XS"){
    observable = 0.;
  } 
  else if (xdistribution_type == "multiplicity"){
    observable = particles[0][0].multiplicity;
  } 
  else if (xdistribution_type == "multiplicity+"){
    observable_max = particles[0][0].multiplicity;
    // FIXME: more general version not ported. still needed?
  }
  else if (xdistribution_type == "muR"){
    // bins at the renormalization scale
    assert(particles.size()==0);
    observable = oset->var_mu_ren;
  } 
  else if (xdistribution_type == "muF"){
    // bins at the factorization scale
    assert(particles.size()==0);
    observable = oset->var_mu_fact;
  } 

  ///////////////////////////////////////////////////////////////////
  //  cumulative taucut distribution of particle/sum of particles  //
  ///////////////////////////////////////////////////////////////////
  
  else if (xdistribution_type == "taucut"){
    // temporary N-jettiness implementation (never applied so far, I guess !!!)
    static int Njettiness = 0;
    static vector<fourvector> NJ_axes;
    Njettiness = oset->csi->n_jet_born;
    //    Njettiness = oset->order_alphas_born;
    NJ_axes.resize(3 + Njettiness); // not generic !!! there should be n_jet_born or so... !!!
    NJ_axes[1] = fourvector(1., 0., 0., 1.);
    NJ_axes[2] = fourvector(1., 0., 0., -1.);
    double tau_jettiness = 0.;
    for (int i_j = 0; i_j < Njettiness; i_j++){
      fourvector temp_jet = oset->particle_event[oset->access_object["jet"]][oset->RA_x_a][i_j].momentum;
      double temp_jet_R = temp_jet.r();
      NJ_axes[3 + i_j] = fourvector(1., temp_jet.x1() / temp_jet_R, temp_jet.x2() / temp_jet_R, temp_jet.x3() / temp_jet_R);
    }
    // all available QCD partons
    stringstream temp;
    logger << LOG_DEBUG_VERBOSE << "oset->ps_runtime_jet_algorithm[oset->RA_x_a].size() = " << oset->ps_runtime_jet_algorithm[oset->RA_x_a].size() << endl;

    for (int i_p = 0; i_p < oset->ps_runtime_jet_algorithm[oset->RA_x_a].size(); i_p++){temp << setw(5) << oset->ps_runtime_jet_algorithm[oset->RA_x_a][i_p];}
    logger << LOG_DEBUG_VERBOSE << "Njettiness = " << Njettiness << "   oset->ps_runtime_jet_algorithm[" << oset->RA_x_a << "].size() = " << oset->ps_runtime_jet_algorithm[oset->RA_x_a].size() << " --- " << temp.str() << endl;
    for (int i_j = 1; i_j < NJ_axes.size(); i_j++){
      logger << LOG_DEBUG_VERBOSE << "NJ_axes[" << i_j << "] = " << NJ_axes[i_j] << endl;
    }    
    // calculate N-jettiness tau
    for (int i_p = 0; i_p < oset->ps_runtime_jet_algorithm[oset->RA_x_a].size(); i_p++){
      vector<double> delta_tau(3 + Njettiness);
      for (int i_j = 1; i_j < NJ_axes.size(); i_j++){
	delta_tau[i_j] = abs(NJ_axes[i_j] * oset->particle_event[0][oset->RA_x_a][oset->ps_runtime_jet_algorithm[oset->RA_x_a][i_p]].momentum);
	//	  delta_tau[i_j] = NJ_axes[i_j] * oset->p_parton[oset->RA_x_a][oset->ps_runtime_jet_algorithm[oset->RA_x_a][i_p]];
	logger << LOG_DEBUG_VERBOSE << "i_p = " << i_p << "   delta_tau[" << i_j << "] = " << delta_tau[i_j] << endl;
	
	if (delta_tau[i_j] < 0.){
	  //	  logger << LOG_DEBUG_VERBOSE << "NJ_axes[" << i_j << "] * oset->particle_event[0][" << oset->RA_x_a << "][" << oset->ps_runtime_jet_algorithm[oset->RA_x_a][i_p] << "] is negative: delta_tau[" << i_j << "] = " << delta_tau[i_j] << endl;

	  delta_tau[i_j] = abs(delta_tau[i_j]);
	}
      }
     
      double min = 1.e99;
      for (int i_j = 1; i_j < NJ_axes.size(); i_j++){
	if (delta_tau[i_j] < min){min = delta_tau[i_j];}
      }
      //      logger << LOG_DEBUG_VERBOSE << "min = " << min << " =?= " << *min_element(delta_tau.begin(), delta_tau.end()) << " = *min_element(delta_tau.begin(), delta_tau.end())" << endl;
      logger << LOG_DEBUG_VERBOSE << "i_p = " << i_p << "   tau_jettiness += " << *min_element(delta_tau.begin(), delta_tau.end()) << endl;
      tau_jettiness += *min_element(delta_tau.begin() + 1, delta_tau.end());
      //	tau_jettiness = tau_jettiness + *min_element(delta_tau.begin(), delta_tau.end());
    }
    logger << LOG_DEBUG_VERBOSE << "tau_jettiness = " << tau_jettiness << endl;
    observable_max += tau_jettiness;
  }



  /////////////////////////////////////////
  //  transverse-momentum distributions  //
  /////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////////////////////
  //  transverse-momentum distribution of particle group (scalar sum when several groups)  //
  //  pT = sum_i pT(particle i);  pT(particle i) = sqrt(p_x1^2 + p_x2^2)                   //
  ///////////////////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "pT") {
    observable = 0.;
    // loop over all particle groups (group=0,1,2,..., 'particle group+1' distribution.dat)  
    for (int group=0; group<particles.size(); group++) {
      observable += reconstructedParticles[group].pT(); // scalar sum of pT of all particle groups
    }
  } 
  ///////////////////////////////////////
  //  pTmax distribution of particles  //
  ///////////////////////////////////////
  else if (xdistribution_type == "pTmax") {
    double max_pT=0;
    for (int group=0; group<particles.size(); group++) {
      double tmp_pT=reconstructedParticles[group].pT();
      if (tmp_pT>max_pT) {
        max_pT=tmp_pT;
      }
    }
    observable = max_pT;
  } 
  ///////////////////////////////////////
  //  pTmin distribution of particles  //
  ///////////////////////////////////////
  else if (xdistribution_type == "pTmin") {
    double min_pT=1e99;
    for (int group=0; group<particles.size(); group++) {
      double tmp_pT=reconstructedParticles[group].pT();
      if (tmp_pT<min_pT) {
        min_pT=tmp_pT;
      }
    }
    observable = min_pT;
  }
  ////////////////////////////////////////////////////
  //  pT ratio of particle/sum of particles         //
  ////////////////////////////////////////////////////
  else if (xdistribution_type == "pTratio") {
    assert(particles.size()==2);
  
    observable = reconstructedParticles[0].pT()/reconstructedParticles[1].pT();
  }
  //////////////////////////////////////
  //  projected ET_miss distribution  //
  //////////////////////////////////////
  else if (xdistribution_type == "ETproj") {
    assert(particles.size()>1);

    // pT from neutrinos is the first reconstructed particle
    double phi_missing = reconstructedParticles[0].phi();
    
    // now find minimum delta phi between this and the other particles
    double dphi_min=0.5*pi;
    for (int group=1; group<particles.size(); group++) {
      double dphi_tmp = min(f2pi-abs(phi_missing-reconstructedParticles[group].phi()),abs(phi_missing-reconstructedParticles[group].phi()));
      if (dphi_tmp<dphi_min) {
        dphi_min = dphi_tmp;
      }
    }
    observable = reconstructedParticles[0].pT()*sin(dphi_min);
  }
  //////////////////////////////////////////////////////////
  //  pTveto 'distribution' of particle/sum of particles  //
  //////////////////////////////////////////////////////////
  else if (xdistribution_type == "pTveto"){
    assert(particles.size()==1);
    observable = reconstructedParticles[0].pT();
  } 
  //////////////////////////////////////////////////////////
  //  'exclusive' pTveto 'distribution' of particle/sum of particles  //
  //////////////////////////////////////////////////////////
  else if (xdistribution_type == "pTvetoexcl"){
    assert(particles.size()==2);
    observable = reconstructedParticles[0].pT();
    observable_max = reconstructedParticles[1].pT();
  } 
  //////////////////////////////////////////////////////////
  //  pT-cut 'distribution' of particle/sum of particles  //
  //////////////////////////////////////////////////////////
  else if (xdistribution_type == "pTcut"){
    for (int group=0; group<particles.size(); group++) {
      observable_max += reconstructedParticles[group].pT();
    }
  } 



  //  new cumulative distributiona for HE-LHC and FCC benchmarks
  
  /////////////////////////////////////////////////////////////
  //  pTmin-cut 'distribution' of particle/sum of particles  //
  /////////////////////////////////////////////////////////////
  else if (xdistribution_type == "pTmincut"){
    observable_max = 1.e99;
    for (int group = 0; group < particles.size(); group++){
      if (observable_max > reconstructedParticles[group].pT()){observable_max = reconstructedParticles[group].pT();}
    }
  } 

  ////////////////////////////////////////////////////////////////
  //  maximal |pseudo-rapidity| veto distribution of particle:  //
  ////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "absetamaxveto"){
    //    assert(particles.size() == 1 && "ERROR in xdistribution_type 'abseta': exactly one particle required");
    observable = 0.;
    for (int group = 0; group < particles.size(); group++){
      if (observable < abs(reconstructedParticles[group].eta())){observable = abs(reconstructedParticles[group].eta());}
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //  maximal |pseudo-rapidity| veto distribution of particle (extra condition: pTlep > 20 GeV):  //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "absetamaxvetopTlep20"){
    //    assert(particles.size() == 1 && "ERROR in xdistribution_type 'abseta': exactly one particle required");
    int cut_condition = 0;
    if (oset->particle_event[oset->access_object["nua"]][oset->RA_x_a].size() > 0){
      if (oset->particle_event[oset->access_object["missing"]][oset->RA_x_a][0].pT < 20.){cut_condition = 1;}
    }

    if (!cut_condition){
      for (int i_p = 0; i_p < oset->particle_event[oset->access_object["lep"]][oset->RA_x_a].size(); i_p++){
	if (oset->particle_event[oset->access_object["lep"]][oset->RA_x_a][i_p].pT < 20.){cut_condition = 1; break;}
      }
    }
    
    if (cut_condition){
      observable = 1.e99;
    }
    else {
      observable = 0.;
      for (int group = 0; group < particles.size(); group++){
	if (observable < abs(reconstructedParticles[group].eta())){observable = abs(reconstructedParticles[group].eta());}
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //  maximal |pseudo-rapidity| veto distribution of particle (extra condition: pTlep > 100 GeV):  //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "absetamaxvetopTlep100"){
    //    assert(particles.size() == 1 && "ERROR in xdistribution_type 'abseta': exactly one particle required");
    int cut_condition = 0;
    if (oset->particle_event[oset->access_object["nua"]][oset->RA_x_a].size() > 0){
      if (oset->particle_event[oset->access_object["missing"]][oset->RA_x_a][0].pT < 100.){cut_condition = 1;}
    }
    if (!cut_condition){
      for (int i_p = 0; i_p < oset->particle_event[oset->access_object["lep"]][oset->RA_x_a].size(); i_p++){
	if (oset->particle_event[oset->access_object["lep"]][oset->RA_x_a][i_p].pT < 100.){cut_condition = 1; break;}
      }
    }
    
    if (cut_condition){
      observable = 1.e99;
    }
    else {
      observable = 0.;
      for (int group = 0; group < particles.size(); group++){
	if (observable < abs(reconstructedParticles[group].eta())){observable = abs(reconstructedParticles[group].eta());}
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //  cumulative invariant-mass 'distribution' of particle/sum of particles  //
  /////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "mcutpTlep20" || xdistribution_type == "McutpTlep20"){
    int cut_condition = 0;
    if (oset->particle_event[oset->access_object["nua"]][oset->RA_x_a].size() > 0){
      if (oset->particle_event[oset->access_object["missing"]][oset->RA_x_a][0].pT < 20.){cut_condition = 1;}
    }
    if (!cut_condition){
      for (int i_p = 0; i_p < oset->particle_event[oset->access_object["lep"]][oset->RA_x_a].size(); i_p++){
	if (oset->particle_event[oset->access_object["lep"]][oset->RA_x_a][i_p].pT < 20.){cut_condition = 1; break;}
      }
    }
    
    if (cut_condition){
      observable_max = -1.;
    }
    else {
      for (int i_g = 0; i_g < particles.size(); i_g++) {
	observable_max += reconstructedParticles[i_g].m();
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //  cumulative invariant-mass 'distribution' of particle/sum of particles  //
  /////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "mcutpTlep100" || xdistribution_type == "McutpTlep100"){
    int cut_condition = 0;
    if (oset->particle_event[oset->access_object["nua"]][oset->RA_x_a].size() > 0){
      if (oset->particle_event[oset->access_object["missing"]][oset->RA_x_a][0].pT < 100.){cut_condition = 1;}
    }
    if (!cut_condition){
      for (int i_p = 0; i_p < oset->particle_event[oset->access_object["lep"]][oset->RA_x_a].size(); i_p++){
	if (oset->particle_event[oset->access_object["lep"]][oset->RA_x_a][i_p].pT < 100.){cut_condition = 1; break;}
      }
    }
    
    if (cut_condition){
      observable_max = -1.;
    }
    else {
      for (int i_g = 0; i_g < particles.size(); i_g++) {
	observable_max += reconstructedParticles[i_g].m();
      }
    }
  }



  ///////////////////////////////////////
  //  transverse-energy distributions  //
  ///////////////////////////////////////


  ///////////////////////////////////////////////////////////
  //  distribution in the (scalar) sum of ET of particles  //
  ///////////////////////////////////////////////////////////
  else if (xdistribution_type == "ET" || xdistribution_type == "HT") {
    double ET_scalar_sum=0;
    for (int group=0; group<particles.size(); group++) {
      ET_scalar_sum += reconstructedParticles[group].ET();
    }
    observable = ET_scalar_sum;
  } 
  /*
  // SK: Should be identical to ET definition:
  ///////////////////////////////////////////////////////////////////////////////
  //  transverse-mass distribution of particle group (sum when several groups) //
  //  HT = sum_i mT(particle i);  mT(particle i) = sqrt(m^2+pT^2)              //
  ///////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "HT") {
    double m_i, pT_i, mT_i;
    observable = 0.;
    // loop over all particle groups (group=0,1,2,..., 'particle group+1' distribution.dat)  
    for (int group=0; group<particles.size(); group++) {
      m_i  = reconstructedParticles[group].m();
      pT_i = reconstructedParticles[group].pT();
      mT_i = sqrt(m_i*m_i+pT_i*pT_i);
      observable += mT_i; // scalar sum of mT of all particle groups
    }
  }
  */


  /////////////////////////////////////
  //  transverse-mass distributions  //
  /////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  transverse-mass distribution of one, two or more particle groups                                         //
  //  one particle group (theoretical definition):                                                             //
  //  mT(particle 1) = sqrt(m^2+pT^2)                                                                          //
  //    -> identical to ET                                                                                     //
  //  two particle groups (experimental definition, used when neutrinos plus one other particle involved):     //
  //  mT = sqrt{[ET(particle 1) + ET(particle 2)]^2 - [pTvec(particle 1)+pTvec(particle 2)]^2};                //
  //  ET(particle i) = sqrt(m_i^2+ pT_i^2)                                                                     //
  //  more than two particle groups (more than one non-neutrino involved):                                     //
  //  mT = sqrt{[ET(particle 1) + sum_i>1 ET(particle i)]^2 - [pTvec(particle 1)+sum_i pTvec(particle i)]^2};  //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "mT" || xdistribution_type == "MT") {
    // use theoretical definition, when only one particle group
    if (particles.size() == 1){
      observable = reconstructedParticles[0].ET(); // theoretical definition of tranverse mass
    }
    // use experimental definition, when two particle group
    else if (particles.size() == 2){
      // all neutrinos should go in the first reconstructed particle
      // experimental definition of tranverse mass when neutrinos and one other particle are involved
      //      observable = sqrt(pow(reconstructedParticles[0].ET() + reconstructedParticles[1].ET(), 2) - (reconstructedParticles[0] + reconstructedParticles[1]).pT2());
      // different definition in MATRIX: !!!
      observable = sqrt(pow(reconstructedParticles[0].pT() + reconstructedParticles[1].pT(), 2) - (reconstructedParticles[0] + reconstructedParticles[1]).pT2());
    }
    else if (particles.size() > 2){
      // all neutrinos should go into the first reconstructed particle
      // experimental definition of tranverse mass when neutrinos and more than one non-neutrino are involved
      fourvector sum_momentum;
      double sum_ET = 0.;
      for (int i_g = 0; i_g < particles.size(); i_g++) {
	//	sum_ET += reconstructedParticles[i_g].ET();
	// different definition in MATRIX: !!!
	sum_ET += reconstructedParticles[i_g].pT();
	sum_momentum = sum_momentum + reconstructedParticles[i_g];
      }
      observable = sqrt(pow(sum_ET, 2) - sum_momentum.pT2());
    }
    /*
    // alternative definition for 2 involved massless particles (E == ET == |pT|):
    double dphi=abs(reconstructedParticles[0].phi()*reconstructedParticles[1].phi());
    observable = sqrt(2*reconstructedParticles[0].x0()*reconstructedParticles[1].x0()*(1.0-cos(dphi)));
    */
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  transverse-mass distribution of particle/sum of particles                                               //
  //  mTCMS = sqrt(2 pT(particle 1) pT(sum_i>1 particle i) [1-cos(deltaPhi(particle1, sum_i>1 particle i))]   //
  //  for two particles same as mT, but only valid for massless particles; for more than two particles pT of  //
  //  vectorial sum instead of scalar sum of pT's (probably wrong!? FIXME)                                    //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "mTCMS" || xdistribution_type == "MTCMS") {
    assert(particles.size()>1  && "ERROR in xdistribution_type mTCMS: at least two particle groups required");
    
    // all neutrinos should go in the first reconstructed particle

    // sum momenta of all particle group after the first one
    fourvector temp_fvl;
    for (int group=1; group<particles.size(); group++) {
      temp_fvl = temp_fvl + reconstructedParticles[group];
    }
    // compute difference in phi between (sum of) neutrino(s) and (sum of) other particle(s)
    double temp_dphi = abs(temp_fvl.phi() - reconstructedParticles[0].phi());
    if (temp_dphi > pi){temp_dphi = 2 * pi - temp_dphi;}

    // compute transverse mass as defined by experiments when neutrinos involved for massless particles
    observable = sqrt(2 * temp_fvl.pT() * reconstructedParticles[0].pT() * (1. - cos(temp_dphi)));
    // double observable2 = sqrt(pow(reconstructedParticles[0].ET() + reconstructedParticles[1].ET(), 2) - (reconstructedParticles[0] + reconstructedParticles[1]).pT2());

    // if (abs(observable - observable2)>0.00001){

    //   cout << observable << endl;
    //   cout << observable2 << endl;
    //   cout << "--------------------------" << endl;
    // }
  }
  //
  // MW: pretty sure this is wrong (pT of vectorial sum, instead of scalar sum of pT's)
  // ==> THIS SHOULD NOT BE USED ANYMORE, mT should not include the correct formula
  //
  else if (xdistribution_type == "mTATLAS" || xdistribution_type == "MTATLAS") {
    fourvector temp_fvm = reconstructedParticles[0];
    double temp_pTm = temp_fvm.pT();
    
    fourvector temp_fvl;
    for (int group=1; group<particles.size(); group++) {
      temp_fvl = temp_fvl + reconstructedParticles[group];
    }
    observable = sqrt(pow(temp_fvl.ET() + temp_pTm, 2) - (temp_fvl + temp_fvm).pT2());
  } 
  ////////////////////////////////////////////////////////////////////////
  //  transverse cluster mass distribution of paticle/sum of particles  //
  ////////////////////////////////////////////////////////////////////////
  //
  // MW: what is this needed for??? --> remove? (FIXME)
  //
  else if (xdistribution_type == "mTcluster" || xdistribution_type == "MTcluster") {
    assert(particles.size()==2);
    //    assert(particles.size()<=2);
    //    if (particles.size()==1) {
    // note: in contrast to old implementation, we *do* sum over the particle momenta
    //      observable = reconstructedParticles[0].mT();
    //    } else {
    observable = sqrt(pow(sqrt(reconstructedParticles[0].m2() + reconstructedParticles[0].pT2()) + reconstructedParticles[1].pT(), 2) - (reconstructedParticles[0] + reconstructedParticles[1]).pT2());
    //    }
  }
  else if (xdistribution_type == "mTZZ") {
    assert(particles.size()==2);
    //    assert(particles.size()<=2);
    //    if (particles.size()==1) {
    // note: in contrast to old implementation, we *do* sum over the particle momenta
    //      observable = reconstructedParticles[0].mT();
    //    } else {
    double mZ = oset->M[23];
    observable = sqrt(pow(sqrt(reconstructedParticles[0].pT2()+mZ*mZ) + sqrt(reconstructedParticles[1].pT2()+mZ*mZ), 2) - (reconstructedParticles[0] + reconstructedParticles[1]).pT2());
    //    }
  }


  ////////////////////////////////////
  //  invariant-mass distributions  //
  ////////////////////////////////////


  /////////////////////////////////////////////////
  //  invariant-mass 'distribution' of particle  //
  /////////////////////////////////////////////////
  else if (xdistribution_type == "m" || xdistribution_type == "M") {
    assert(particles.size() == 1 && "ERROR in xdistribution_type 'm/M': exactly one particle required");
    observable = reconstructedParticles[0].m();
  } 
  //////////////////////////////////////////////////////////////////////////
  //  invariant-mass-difference distribution of two particles: m_1 - m_2  //
  //  m(particle i) = sqrt(p^2)                                           //
  //////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "dm" || xdistribution_type == "dM") {
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'dm/dM': exactly two particles required");
    observable = reconstructedParticles[0].m() - reconstructedParticles[1].m();
  }
  //////////////////////////////////////////////////////////////////////////////
  //  |invariant-mass-difference| distribution of two particles: |m_1 - m_2|  //
  //////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "absdm" || xdistribution_type == "absdM") {
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'absdm/absdM': exactly two particles required");
    observable = abs(reconstructedParticles[0].m() - reconstructedParticles[1].m());
  }
  ////////////////////////////////////////////////////////////
  //  maximum invariant-mass distribution of two particles  //
  ////////////////////////////////////////////////////////////
  else if (xdistribution_type == "mmax" || xdistribution_type == "Mmax") {
    // could be extended to >2 particles
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'mmax/Mmax': exactly two particles required");
    observable = max(reconstructedParticles[0].m(),reconstructedParticles[1].m());
  }
  ////////////////////////////////////////////////////////////
  //  minimum invariant-mass distribution of two particles  //
  ////////////////////////////////////////////////////////////
  else if (xdistribution_type == "mmin" || xdistribution_type == "Mmin") {
    // could be extended to >2 particles
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'mmin/Mmin': exactly two particles required");
    observable = min(reconstructedParticles[0].m(),reconstructedParticles[1].m());
  } 
  /////////////////////////////////////////////////////////////////////////////
  //  cumulative invariant-mass 'distribution' of particle/sum of particles  //
  /////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "mcut" || xdistribution_type == "Mcut"){
    for (int i_g = 0; i_g < particles.size(); i_g++) {
      observable_max += reconstructedParticles[i_g].m();
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////
  //  cumulative |invariant-mass-difference| distribution of particle/sum of particles  //
  ////////////////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "absdmcut" || xdistribution_type == "absdMcut"){
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'absdmcut/absdMcut': exactly two particles required");
    observable_max = abs(reconstructedParticles[0].m() - reconstructedParticles[1].m());
  }


  ///////////////////////////////////////
  //  (pseudo-)rapidity distributions  //
  ///////////////////////////////////////


  ////////////////////////////////////////////
  //  rapidity distribution of particle: y  //
  ////////////////////////////////////////////
  else if (xdistribution_type == "y") {
    assert(particles.size() == 1 && "ERROR in xdistribution_type 'y': exactly one particle required");
    observable = reconstructedParticles[0].rapidity();
  } 
  //////////////////////////////////////////////////////////////////////////
  //  pseudo-rapidity distribution of particle: eta (= artanh(p_x3/|p|))  //
  //////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "eta") {
    assert(particles.size() == 1 && "ERROR in xdistribution_type 'eta': exactly one particle required");
    observable = reconstructedParticles[0].eta();
  }
  ////////////////////////////////////////////////
  //  |rapidity| distribution of particle: |y|  //
  ////////////////////////////////////////////////
  else if (xdistribution_type == "absy") {
    assert(particles.size() == 1 && "ERROR in xdistribution_type 'absy': exactly one particle required");
    observable = abs(reconstructedParticles[0].rapidity());
  }
  /////////////////////////////////////////////////////////
  //  |pseudo-rapidity| distribution of particle: |eta|  //
  /////////////////////////////////////////////////////////
  else if (xdistribution_type == "abseta") {
    assert(particles.size() == 1 && "ERROR in xdistribution_type 'abseta': exactly one particle required");
    observable = abs(reconstructedParticles[0].eta());
  }
  //////////////////////////////////////////////////////////////////////
  //  rapidity-difference distribution of two particles: (y_1 - y_2)  //
  //////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "dy") {
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'dy': exactly two particles required");
    observable = reconstructedParticles[0].rapidity() - reconstructedParticles[1].rapidity();
  }
  /////////////////////////////////////////////////////////////////////////////////
  //  pseudo-rapidity-difference distribution of two particles: (eta_1 - eta_2)  //
  /////////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "deta") {
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'deta': exactly two particles required");
    observable = reconstructedParticles[0].eta() - reconstructedParticles[1].eta();
  }
  else if (xdistribution_type == "dabsy") {
    ////////////////////////////////////////////////////////////////////////
  //  |rapidity|-difference distribution of two particles: |y_1| - |y_2|  //
    ////////////////////////////////////////////////////////////////////////
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'dabsy': exactly two particles required");
    observable = abs(reconstructedParticles[0].rapidity()) - abs(reconstructedParticles[1].rapidity());
  }
  /////////////////////////////////////////////////////////////////////////////////////
  //  |pseudo-rapidity|-difference distribution of two particles: |eta_1| - |eta_2|  //
  /////////////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "dabseta") {
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'dabseta': exactly two particles required");
    observable = abs(reconstructedParticles[0].eta()) - abs(reconstructedParticles[1].eta());
  }
  ////////////////////////////////////////////////////////////////////////
  //  |rapidity-difference| distribution of two particles: |y_1 - y_2|  //
  ////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "absdy") {
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'absdy': exactly two particles required");
    observable = abs(reconstructedParticles[0].rapidity() - reconstructedParticles[1].rapidity());
  }
  ///////////////////////////////////////////////////////////////////////////////////
  //  |pseudo-rapidity-difference| distribution of two particles: |eta_1 - eta_2|  //
  ///////////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "absdeta") {
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'absdeta': exactly two particles required");
    observable = abs(reconstructedParticles[0].eta() - reconstructedParticles[1].eta());
  } 
  else if (xdistribution_type == "absdabsy") {
    //////////////////////////////////////////////////////////////////////////////
    //  ||rapidity|-difference| distribution of two particles: ||y_1| - |y_2||  //
    //////////////////////////////////////////////////////////////////////////////
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'absdabsy': exactly two particles required");
    observable = abs(abs(reconstructedParticles[0].rapidity()) - abs(reconstructedParticles[1].rapidity()));
  } 
  else if (xdistribution_type == "absdabseta") {
    /////////////////////////////////////////////////////////////////////////////////////////
    //  ||pseudo-rapidity|-difference| distribution of two particles: ||eta_1| - |eta_2||  //
    /////////////////////////////////////////////////////////////////////////////////////////
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'dabsdeta': exactly two particles required");
    observable = abs(abs(reconstructedParticles[0].eta()) - abs(reconstructedParticles[1].eta()));
  } 
  ///////////////////////////////////////////////////
  //  maximum-|rapidity| distribution of particle  //
  ///////////////////////////////////////////////////
  else if (xdistribution_type == "absymax") {
    // FIXME: we might want to change this such that it computes the maximum eta of the *reconstructed* particles
    assert(particles.size() == 1 && "ERROR in xdistribution_type 'absymax': exactly one particle required");
    double max_y=0;
    for (int particle=0; particle<particles[0].size(); particle++) {
      double tmp_y=abs(particles[0][particle].momentum.rapidity());
      if (tmp_y>max_y) {
        max_y=tmp_y;
      }
    }
    observable = max_y;
  } 
  //////////////////////////////////////////////////////////
  //  maximum-|pseudo-rapidity| distribution of particle  //
  //////////////////////////////////////////////////////////
  else if (xdistribution_type == "absetamax") {
    // different definition in MATRIX !!!
    // FIXME: we might want to change this such that it computes the maximum eta of the *reconstructed* particles
    /*
    assert(particles.size() == 1 && "ERROR in xdistribution_type 'absetamax': exactly one particle required");
    double max_eta=0;
    for (int particle=0; particle<particles[0].size(); particle++) {
      double tmp_eta=abs(particles[0][particle].momentum.eta());
      if (tmp_eta>max_eta) {
        max_eta=tmp_eta;
      }
    }
    observable = max_eta;
    */
    observable = 0.;
    for (int group = 0; group < particles.size(); group++){
      if (observable < abs(reconstructedParticles[group].eta())){observable = abs(reconstructedParticles[group].eta());}
    }

  } 
  ///////////////////////////////////////////////////
  //  minimum-|rapidity| distribution of particle  //
  ///////////////////////////////////////////////////
  else if (xdistribution_type == "absymin") {
    // FIXME: we might want to change this such that it computes the maximum eta of the *reconstructed* particles
    assert(particles.size() == 1 && "ERROR in xdistribution_type 'absymin': exactly one particle required");
    
    double min_y=1.e99;
    for (int particle=0; particle<particles[0].size(); particle++) {
      double tmp_y=abs(particles[0][particle].momentum.rapidity());
      if (tmp_y<min_y) {
        min_y=tmp_y;
      }
    }
    observable = min_y;
  }
  //////////////////////////////////////////////////////////
  //  minimum-|pseudo-rapidity| distribution of particle  //
  //////////////////////////////////////////////////////////
  else if (xdistribution_type == "absetamin") {
    // FIXME: we might want to change this such that it computes the maximum eta of the *reconstructed* particles
    assert(particles.size() == 1 && "ERROR in xdistribution_type 'absetamin': exactly one particle required");
    double min_eta=1.e99;
    for (int particle=0; particle<particles[0].size(); particle++) {
      double tmp_eta=abs(particles[0][particle].momentum.eta());
      if (tmp_eta<min_eta) {
        min_eta=tmp_eta;
      }
    }
    observable = min_eta;
  } 


  ///////////////////////////////////////////
  //  angular (and similar) distributions  //
  ///////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////
  //  transverse-angle distribution of particle a                  //
  //  if particle_2nd != 0, transverse-angle difference of a & b   //
  ///////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "phi") {
    assert((particles.size() == 1 || particles.size() == 2) && "ERROR in xdistribution_type 'phi': exactly one or two particles required");
    if (particles.size() == 1) {
      observable = reconstructedParticles[0].phi();
    } else {
      double dphi_abs = abs(reconstructedParticles[0].phi()-reconstructedParticles[1].phi());
      observable = min(f2pi-dphi_abs,dphi_abs);
    }
  }
  /////////////////////////////////////////////////
  //  costheta distribution of particle          //
  //  if particle_2nd != 0, angle between 1 & 2  //
  /////////////////////////////////////////////////
  else if (xdistribution_type == "costheta") {
    assert((particles.size() == 1 || particles.size() == 2) && "ERROR in xdistribution_type 'costheta': exactly one or two particles required");
    if (particles.size()==1) {
      observable = reconstructedParticles[0].costheta();
    } else {
      observable = cosangle(reconstructedParticles[0],reconstructedParticles[1]);
    }
  }
  ///////////////////////////////////////////////////////////////////
  //  transverse-angle--rapidity plane distribution of particle a  //
  ///////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "dR") {
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'dR': exactly two particles required");
    double temp_dphi = abs(reconstructedParticles[0].phi() - reconstructedParticles[1].phi());
    if (temp_dphi > pi){temp_dphi = 2 * pi - temp_dphi;}
    double dy = reconstructedParticles[0].rapidity() - reconstructedParticles[1].rapidity();
    observable = sqrt(pow(dy, 2) + pow(temp_dphi, 2));
  } 
  //////////////////////////////////////////////////////////////////////////
  //  transverse-angle--pseudo-rapidity plane distribution of particle a  //
  //////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "dReta") {
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'dReta': exactly two particles required");
    double temp_dphi = abs(reconstructedParticles[0].phi() - reconstructedParticles[1].phi());
    if (temp_dphi > pi){temp_dphi = 2 * pi - temp_dphi;}
    double deta = reconstructedParticles[0].eta() - reconstructedParticles[1].eta();
    observable = sqrt(pow(deta, 2) + pow(temp_dphi, 2));
  } 
  ///////////////////////////////////////////////////////////////////
  //      costhetastar variable for diphoton and Drell-Yan         //
  //       definition as in Collins and Soper, PRD 16 (1977) 2219  //
  ///////////////////////////////////////////////////////////////////   
  else if (xdistribution_type == "costhetastar") {
    assert(particles.size() == 2 && "ERROR in xdistribution_type 'costhetastar': exactly two particles required");
    double pt3 = reconstructedParticles[0].pT();
    double pt4 = reconstructedParticles[1].pT();
    double eta3 = reconstructedParticles[0].eta();
    double eta4 = reconstructedParticles[1].eta();
    fourvector p34 = reconstructedParticles[0]+reconstructedParticles[1];
    double m34 = p34.m();
    double pt34 = p34.pT(); 
    observable = abs(2 * pt3 * pt4 * sinh(eta3 - eta4) / m34 / sqrt(pow(m34, 2) + pow(pt34, 2)));
  }


  ///////////////////////////////////////////////////////////////////////////////////////
  //  special distributions, which should become obsolete with user-defined particles  //
  ///////////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////
  //  delta phi of reconstructed Z's according to the CMS pairing prescription  //
  //  FIXME: should be replaced by a user defined distribution                  //
  ////////////////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "phi_pair_CMS") {
    assert(particles.size()==2);
    assert(particles[0].size()==2);
    assert(particles[1].size()==2);
  
    double dMZ1 = abs((particles[0][0].momentum+particles[1][0].momentum).m()-oset->M[23]);
    double dMZ2 = abs((particles[0][0].momentum+particles[1][1].momentum).m()-oset->M[23]);
    double dMZ3 = abs((particles[0][1].momentum+particles[1][0].momentum).m()-oset->M[23]);
    double dMZ4 = abs((particles[0][1].momentum+particles[1][1].momentum).m()-oset->M[23]);
    
    fourvector Z1,Z2;
    
    double dMZmin = min(dMZ1,min(dMZ2,min(dMZ3,dMZ4)));
    
    if (dMZ1 == dMZmin) {
      Z1 = particles[0][0].momentum+particles[1][0].momentum;
      Z2 = particles[0][1].momentum+particles[1][1].momentum;
    } else if (dMZ2 == dMZmin) {
      Z1 = particles[0][0].momentum+particles[1][1].momentum;
      Z2 = particles[0][1].momentum+particles[1][0].momentum;
    } else if (dMZ3 == dMZmin) {
      Z1 = particles[0][1].momentum+particles[1][0].momentum;
      Z2 = particles[0][0].momentum+particles[1][1].momentum;
    } else {
      assert(dMZ4==dMZmin);
      Z1 = particles[0][1].momentum+particles[1][1].momentum;
      Z2 = particles[0][0].momentum+particles[1][0].momentum;
    }
      
    observable = min(f2pi-abs(Z1.phi()-Z2.phi()),abs(Z1.phi()-Z2.phi()));
  }
  ///////////////////////////////////////////////////////////////////
  //  invariant mass of the subleading Z, according to the CMS pairing prescription //
  //  FIXME: should be replaced by a user defined distribution   //
  ///////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "m34_CMS") {

    // each particle contains a list of Z boson candidates
    double dMZ_min = 1e99;
    int Z_group=-1;
    int first_Z=-1;

    for (int group=0; group<particles.size(); group++) {
      for (int particle=0; particle<particles[group].size()/2; particle++) {
        double dMZ = abs((particles[group][2*particle].momentum+particles[group][2*particle+1].momentum).m()-oset->M[23]);
        if (dMZ < dMZ_min) {
          Z_group = group;
          first_Z = particle;
          dMZ_min = dMZ;
        }
      }
    }
    fourvector Z2 = particles[Z_group][(2*first_Z+2)%4].momentum+particles[Z_group][(2*first_Z+1+2)%4].momentum;
  
    observable = Z2.m();
  }
  /////////////////////////////////////////////////////////////////////
  //  costheta of the leading Z, evaluated in the rest frame of the  //
  //  four leptons, according to the CMS pairing prescription        //
  //  FIXME: should be replaced by a user defined distribution       //
  /////////////////////////////////////////////////////////////////////
  else if (xdistribution_type == "costheta12_CMS") {
    
    // each particle contains a list of Z boson candidates
    double dMZ_min = 1e99;
    int Z_group=-1;
    int first_Z=-1;

    for (int group=0; group<particles.size(); group++) {
      for (int particle=0; particle<particles[group].size()/2; particle++) {
        double dMZ = abs((particles[group][2*particle].momentum+particles[group][2*particle+1].momentum).m()-oset->M[23]);
        if (dMZ < dMZ_min) {
          Z_group = group;
          first_Z = particle;
          dMZ_min = dMZ;
        }
      }
    }
    
    fourvector Z1 = particles[Z_group][2*first_Z].momentum+particles[Z_group][2*first_Z+1].momentum;
    fourvector Z2 = particles[Z_group][(2*first_Z+2)%4].momentum+particles[Z_group][(2*first_Z+1+2)%4].momentum;

    fourvector Z_boosted = Z1.boost(Z1+Z2);
  
    observable = abs(Z_boosted.costheta());
  }

  else {
    cout << "distribution not implemented: " << xdistribution_type << endl;
    assert(false);
  }
}

xdistribution::xdistribution(string _xdistribution_name, string _xdistribution_type, vector<vector<int> > _all_particle, vector<vector<int> > _all_particle_group, vector<vector<int> > _required_subparticle, double _start, double _end, int _n_bins, double _step, string _edges, string _type_binning, observable_set *_oset){
  Logger logger("xdistribution::xdistribution");
  logger << LOG_DEBUG << "started" << endl;

  xdistribution_name = _xdistribution_name;
  xdistribution_type = _xdistribution_type;
  
  oset = _oset;
  
  for (int i_p=0; i_p<_all_particle.size(); i_p++) {
    vector<ParticleID> new_group;
    for (int i_s=0; i_s<_all_particle[i_p].size(); i_s++) {
      logger << LOG_DEBUG << "new particle " << i_p << ", " << i_s << endl;
      ParticleID new_particle;
      new_particle.type=_all_particle_group[i_p][i_s];
      new_particle.number=_all_particle[i_p][i_s];
      logger << LOG_DEBUG << new_particle.type << ", " << new_particle.number << endl;
      if (_required_subparticle[i_p][i_s] == 2) {
        new_particle.required_definition = true;
      } else {
         new_particle.required_definition = false;
      }
      if (_required_subparticle[i_p][i_s] == 1) {
        new_particle.required_exists = true;
      } else {
         new_particle.required_exists = false;
      }
      if (_required_subparticle[i_p][i_s] == -1) {
        // FIXME: what does this mean?
        assert(false);
      }
      new_group.push_back(new_particle);
    }
    particles.push_back(new_group);
  }
  
  reconstructedParticles.resize(particles.size());
  
  
  
  if (xdistribution_type == "multiplicity+" || 
      xdistribution_type == "pTcut" || 
      xdistribution_type == "mcut" || 
      xdistribution_type == "Mcut" || 
      xdistribution_type == "absdmcut" || 
      xdistribution_type == "absdMcut" || 
      xdistribution_type == "taucut" || 
      xdistribution_type == "pTmincut" || 
      xdistribution_type == "mcutpTlep20" || 
      xdistribution_type == "McutpTlep20" || 
      xdistribution_type == "mcutpTlep100" || 
      xdistribution_type == "McutpTlep100") {
    typeCumulative = CUMULATIVE_UPPER;
    // observable >= 'left bin edge x' -> bins from 0 to 'observable_max' are filled (check if +1 or not !!!).
  }
  else if (xdistribution_type == "pTveto" || 
	   xdistribution_type == "absetamaxveto" || 
	   xdistribution_type == "absetamaxvetopTlep20" || 
	   xdistribution_type == "absetamaxvetopTlep100") {
    typeCumulative = CUMULATIVE_LOWER;
  }
  else if (xdistribution_type == "pTvetoexcl") {
    typeCumulative = CUMULATIVE_BOTH;
  }
  else {
    typeCumulative = CUMULATIVE_NONE;
    // 'left bin edge x' <= observable < 'left bin edge x+1' -> bin 'observable = x' is filled.
  }
  
//   if (xdistribution_type == "pT"){xdistribution_number =  1;}
//   else if (xdistribution_type == "eta"){xdistribution_number =  2;}
//   else if (xdistribution_type == "pTmax"){xdistribution_number =  2;}
//   else if (xdistribution_type == "pTmin"){xdistribution_number =  2;}
//   else if (xdistribution_type == "m"){xdistribution_number =  3;}
//   else if (xdistribution_type == "dm"){xdistribution_number =  87;}
//   else if (xdistribution_type == "absdm"){xdistribution_number =  88;}
//   else if (xdistribution_type == "mmin"){xdistribution_number =  52;}
//   else if (xdistribution_type == "phi"){xdistribution_number =  4;}
//   else if (xdistribution_type == "phi_pair_CMS"){xdistribution_number =  49;}
//   else if (xdistribution_type == "dpT_djet_pT"){xdistribution_number =  50;}
//   else if (xdistribution_type == "costheta"){xdistribution_number =  6;}
//   else if (xdistribution_type == "pTratio"){xdistribution_number =  12;}
//   else if (xdistribution_type == "absetamax"){xdistribution_number =  13;}
//   else if (xdistribution_type == "absetamin"){xdistribution_number =  14;}
//   else if (xdistribution_type == "ET"){xdistribution_number =  71;}
//   else if (xdistribution_type == "ETproj"){xdistribution_number =  72;}
//   else if (xdistribution_type == "mT"){xdistribution_number =  15;}
//   else if (xdistribution_type == "mTATLAS"){xdistribution_number =  62;}
//   else if (xdistribution_type == "mTCMS"){xdistribution_number =  63;}
//   else if (xdistribution_type == "mTcluster"){xdistribution_number =  16;}
//   else if (xdistribution_type == "y"){xdistribution_number =  17;}
//   else if (xdistribution_type == "ymax"){xdistribution_number =  8;}
//   else if (xdistribution_type == "ymin"){xdistribution_number =  9;}
//   else if (xdistribution_type == "deta"){xdistribution_number =  18;}
//   else if (xdistribution_type == "dabseta"){xdistribution_number =  19;}
//   else if (xdistribution_type == "absdeta"){xdistribution_number =  85;}
//   else if (xdistribution_type == "absdabseta"){xdistribution_number =  86;}
//   else if (xdistribution_type == "dy"){xdistribution_number =  20;}
//   else if (xdistribution_type == "dabsy"){xdistribution_number =  21;}
//   else if (xdistribution_type == "absdy"){xdistribution_number =  26;}
//   else if (xdistribution_type == "absdabsy"){xdistribution_number =  27;}
//   else if (xdistribution_type == "absymax"){xdistribution_number =  24;}
//   else if (xdistribution_type == "absymin"){xdistribution_number =  25;}
//   else if (xdistribution_type == "dR"){xdistribution_number =  30;}
//   else if (xdistribution_type == "dReta"){xdistribution_number =  33;}
//   else if (xdistribution_type == "XS"){xdistribution_number =  98;}
//   else if (xdistribution_type == "multiplicity"){xdistribution_number =  99;}
//   else if (xdistribution_type == "multiplicity+"){xdistribution_number =  299;}
//   else if (xdistribution_type == "pTveto"){xdistribution_number =  101;}
//   else if (xdistribution_type == "pTvetoexcl"){xdistribution_number =  301;}
//   else if (xdistribution_type == "pTcut"){xdistribution_number =  201;}
//   else if (xdistribution_type == "mcut"){xdistribution_number =  202;}
//   else if (xdistribution_type == "absdmcut"){xdistribution_number =  203;}
//   else {xdistribution_number =  atoi(xdistribution_type.c_str());}
//   
//   if (xdistribution_number>300) {
//     typeCumulative = CUMULATIVE_BOTH;
//   } else if (xdistribution_number>200) {
//     typeCumulative = CUMULATIVE_UPPER;
//   } else if (xdistribution_number>100) {
//     typeCumulative = CUMULATIVE_LOWER;
//   } else {
//     typeCumulative = NORMAL;
//   } 

  all_particle = _all_particle;
  all_particle_group = _all_particle_group;
  required_subparticle = _required_subparticle;
  logger << LOG_DEBUG << "xdistribution_name = " << _xdistribution_name << endl;
  logger << LOG_DEBUG << "type_binning = " << _type_binning << endl;
  logger << LOG_DEBUG << "start = " << _start << endl;
  logger << LOG_DEBUG << "end = " << _end << endl;
  logger << LOG_DEBUG << "step = " << _step << endl;
  logger << LOG_DEBUG << "n_bins = " << _n_bins << endl;
  logger << LOG_DEBUG << "required_subparticle.size() = " << required_subparticle.size() << endl;
  logger << LOG_DEBUG << "required_particle.size() = " << required_particle.size() << endl;
  required_particle.resize(required_subparticle.size());
  for (int i_p = 0; i_p < required_particle.size(); i_p++){
    int flag = 0;
    logger << LOG_DEBUG << "required_subparticle[" << i_p << "].size() = " << required_subparticle[i_p].size() << endl;
    for (int i_s = 0; i_s < required_subparticle[i_p].size(); i_s++){
      if (required_subparticle[i_p][i_s] >= 0){flag = GSL_MIN(1, required_subparticle[i_p][i_s]); break;}
    }
    if (flag > 0){required_particle[i_p] = flag;}
  logger << LOG_DEBUG << "required_particle[" << i_p << "] = " << required_particle[i_p] << endl;
  }
  logger << LOG_DEBUG << "done!" << endl;
  type_binning = _type_binning;
  logger << LOG_DEBUG << "type_binning = " << type_binning << endl;

  int type_input = 0;
  if (_n_bins == 0){type_input = 1;}
  else if (_step == 0.){type_input = 4;}
  else if (_end == 0.){type_input = 3;}
  else if (_start == 0.){type_input = 2;}

  if (type_binning == "linear"){
    logger << LOG_DEBUG << "linear: type_input = " << type_input << endl;
    if (type_input == 0){
      n_bins = _n_bins;
      start = _start;
      end = _end;
      step = _step;
      if ((end - start) != (n_bins * step)){logger << LOG_FATAL << "ill-defined distribution " << xdistribution_name << " (" << xdistribution_type << "!" << endl; exit(1);}
    }
    else if (type_input == 1){
      // no n_bins
      start = _start;
      end = _end;
      step = _step;
      assert(step>0);
      n_bins = int((end - start) / step);
    }
    else if (type_input == 2){
      // no start
      n_bins = _n_bins;
      step = _step;
      end = _end;
      start = end - n_bins * step;
    }
    else if (type_input == 3){
      // no end
      n_bins = _n_bins;
      start = _start;
      step = _step;
      end = start + n_bins * step;
    }
    else if (type_input == 4){
      // no step
      n_bins = _n_bins;
      assert(n_bins>0);
      start = _start;
      end = _end;
      step = (end - start) / n_bins;
    }

    bin_edge.resize(n_bins + 1);
    bin_width.resize(n_bins);
    for (int i_b = 0; i_b < n_bins; i_b++){
      bin_edge[i_b] = start + i_b * step;
      bin_width[i_b] = step;
    }
    bin_edge[n_bins] = end;

    logger << LOG_DEBUG << "start = " << start << endl;
    logger << LOG_DEBUG << "end = " << end << endl;
    logger << LOG_DEBUG << "step = " << step << endl;
    logger << LOG_DEBUG << "n_bins = " << n_bins << endl;

  }
  else if (_type_binning == "logarithmic"){
    if (type_input != 4){logger << LOG_FATAL << "ill-defined distribution; equally-sized bins are not compatible with logarithmic binning " << xdistribution_name << " (" << xdistribution_type << "!" << endl; exit(1);}
    n_bins = _n_bins;
    if (_start <= 0 || _end <= 0){logger << LOG_FATAL << "ill-defined distribution; start/end are not compatible with logarithmic binning " << xdistribution_name << " (" << xdistribution_type << "!" << endl; exit(1);}
    start = log10(_start);
    end = log10(_end);
    step = (end - start) / n_bins;
    bin_edge.resize(n_bins + 1);
    bin_width.resize(n_bins);
    for (int i_b = 0; i_b < n_bins; i_b++){
      bin_edge[i_b] = pow(10., start + i_b * step);
      bin_width[i_b] = bin_edge[i_b] * (pow(10., step) - 1);
      //      bin_edge[i_b] = start + i_b * step;
      //      bin_width[i_b] = pow(10., bin_edge[i_b]) * (pow(10., step) - 1);
      logger << LOG_DEBUG << "bin_edge[" << i_b << "] = " << setprecision(8) << setw(15) << bin_edge[i_b] << "   " << setprecision(8) << setw(15) << bin_width[i_b] << endl;
    }
    bin_edge[n_bins] = pow(10., end);
  }
  else if (_type_binning == "irregular"){
    edges = _edges;
    logger << LOG_DEBUG << "irregular: edges = " << edges << endl;
    vector<string> vs_edge(1);
    for (int i_b = 0; i_b < edges.size(); i_b++){
      logger << LOG_DEBUG << "i_b = " << i_b << "   vs_edge.size() = " << vs_edge.size() << endl;
      if (edges[i_b] == ':'){vs_edge.push_back("");}
      else if (edges[i_b] != ':'){vs_edge[vs_edge.size() - 1].push_back(edges[i_b]);}
    }
    n_bins = vs_edge.size() - 1;
    for (int i_b = 0; i_b < n_bins; i_b++){logger << LOG_DEBUG << "vs_edge[" << i_b << "] = " << setprecision(8) << setw(15) << vs_edge[i_b] << endl;}
    bin_edge.resize(n_bins + 1);
    bin_width.resize(n_bins);
    for (int i_b = 0; i_b < n_bins + 1; i_b++){bin_edge[i_b] = atof(vs_edge[i_b].c_str());}
    for (int i_b = 0; i_b < n_bins; i_b++){bin_width[i_b] = bin_edge[i_b + 1] - bin_edge[i_b];}
    for (int i_b = 0; i_b < n_bins; i_b++){if (bin_width[i_b] < 0.){logger << LOG_FATAL << "ill-defined distribution; irregular bins are not correctly ordered: " << xdistribution_name << " (" << xdistribution_type << "!" << endl; exit(1);}}

    step = 1.;
    logger << LOG_DEBUG << xdistribution_name << " (" << xdistribution_type << ")" << endl;     
    for (int i_b = 0; i_b < n_bins; i_b++){logger << LOG_DEBUG << "bin_edge[" << i_b << "] = " << setprecision(8) << setw(15) << bin_edge[i_b] << "   " << setprecision(8) << setw(15) << bin_width[i_b] << endl;}

    start = bin_edge[0];
    end = bin_edge[n_bins];
  }


  if      ((xdistribution_type == "pT") ||
	   (xdistribution_type == "pTratio") || 
	   (xdistribution_type == "pTveto") || 
	   (xdistribution_type == "pTvetoexcl") || 
	   (xdistribution_type == "pTcut") ||
	   (xdistribution_type == "absy") ||
	   (xdistribution_type == "abseta") ||
	   (xdistribution_type == "dabsy") ||
	   (xdistribution_type == "dabseta") ||
	   (xdistribution_type == "absdy") ||
	   (xdistribution_type == "absdeta") ||
	   (xdistribution_type == "absdabsy") ||
	   (xdistribution_type == "absdabseta") ||
	   (xdistribution_type == "absymax") ||
	   (xdistribution_type == "absetamax") ||
	   (xdistribution_type == "absymin") ||
	   (xdistribution_type == "absetamin") ||
	   (xdistribution_type == "m") ||
	   (xdistribution_type == "mmin") ||
	   (xdistribution_type == "mmax") ||
	   (xdistribution_type == "mcut") ||
	   (xdistribution_type == "dm") ||
	   (xdistribution_type == "absdm") ||
	   (xdistribution_type == "absdmcut") ||
	   (xdistribution_type == "dR") ||
	   (xdistribution_type == "ET") ||
	   (xdistribution_type == "ETproj") ||
	   (xdistribution_type == "mT") ||
	   (xdistribution_type == "mTATLAS") ||
	   (xdistribution_type == "mTCMS") ||
	   (xdistribution_type == "m34_CMS") ||
	   (xdistribution_type == "mTcluster") ||
	   (xdistribution_type == "XS") ||
	   (xdistribution_type == "taucut") ||
	   (xdistribution_type == "muR") ||
	   (xdistribution_type == "muF") ||
	   (xdistribution_type == "multiplicity") || 
	   (xdistribution_type == "multiplicity+") || 
	   xdistribution_type == "pTmincut" || 
	   xdistribution_type == "mcutpTlep20" || 
	   xdistribution_type == "McutpTlep20" || 
	   xdistribution_type == "mcutpTlep100" || 
	   xdistribution_type == "McutpTlep100" || 
	   xdistribution_type == "absetamaxveto" || 
	   xdistribution_type == "absetamaxvetopTlep20" || 
	   xdistribution_type == "absetamaxvetopTlep100"){symm = 0;}
  else if ((xdistribution_type == "y") || 
	   (xdistribution_type == "eta") || 
	   (xdistribution_type == "dy") || 
	   (xdistribution_type == "deta" || 
	   (xdistribution_type == "etadiff"))){
    // make sure distribution is defined symmetrically
    assert(start == -end);
    symm = 1;
  }
  else if ((xdistribution_type == "phi") ||
	   (xdistribution_type == "phi_pair_CMS") ||
	   (xdistribution_type == "costheta12_CMS") ||
	   (xdistribution_type == "dpT_djet_pT") ||
	   (xdistribution_type == "costheta")){
    if (all_particle.size() == 1){
      // make sure distribution is defined symmetrically
      assert(start == -end);
      symm = 1; // ???
    }
    //    if (particle_2nd.size() == 0){symm = 1;}
    else {symm = 0;}
  }

  logger << LOG_DEBUG << "finished" << endl;
}

ostream & operator << (ostream &s, const xdistribution &sp){
  s << setw(25) << sp.xdistribution_name << "   ";// << endl;
  s << "type: " << setw(15) << sp.xdistribution_type << "   ";// << endl; 
  //  s << endl << setw(15) << "";
  s << "interval: " << setw(5) << sp.start << " - " << setw(4) << sp.step << "(" << setw(3) << sp.n_bins << ") - " << setw(5) << sp.end << "   ";
  s << "particles: ";
  for (int i_p = 0; i_p < (sp.all_particle).size(); i_p++){
    s << "[";
    for (int i_s = 0; i_s < (sp.all_particle[i_p]).size(); i_s++){
      s << sp.all_particle[i_p][i_s] << "(" << sp.all_particle_group[i_p][i_s] << ", ";
      if (sp.required_subparticle[i_p][i_s] == -1){s << "-";}
      else if (sp.required_subparticle[i_p][i_s] == 0){s << "0";}
      else if (sp.required_subparticle[i_p][i_s] == 1){s << "+";}
      s << ")";
      if (i_s < sp.all_particle[i_p].size() - 1){s << ", ";}
    }
    s << ", ";
    if (sp.required_particle[i_p] == 0){s << "-";}
    //    if (sp.required_particle[i_p] == -1){s << "-";}
    //    else if (sp.required_particle[i_p] == 0){s << "0";}
    else if (sp.required_particle[i_p] == 1){s << "+";}
    s << "]";

    if (i_p < sp.all_particle.size() - 1){s << " - ";}
  }
  s << "   ";// << endl;
  /*
  s << "1st: ";
  for (int i = 0; i < (sp.particle).size(); i++){
    s << sp.particle[i];
    if (i < sp.particle.size() - 1){s << ", ";}
  }
  s << "   ";// << endl;
  if (sp.particle_2nd.size() != 0){
    s << "2nd: = ";
    for (int i = 0; i < (sp.particle_2nd).size(); i++){
      s << sp.particle_2nd[i];
      if (i < sp.particle_2nd.size() - 1){s << ", ";}
    }
    //    s << endl;
  }
  */
  //  s << endl;
  /*
  s << endl << setw(52) << "";
  s << "n_bin: " << setw(4) << sp.n_bins << "   ";// << endl; 
  s << "start: " << setw(5) << sp.start << "   ";// << endl; 
  s << "end: " << setw(5) << sp.end << "   ";// << endl; 
  s << "width: " << setw(3) << sp.step << "   ";// << endl; 
  s << "symm: " << setw(2) << sp.symm << "   ";// << endl; 
  */
  return s;
  /*
  s << "xdistribution name  = " << setw(15) << sp.xdistribution_name << "   ";// << endl;
  
  s << "xdistribution type  = " << setw(15) << sp.xdistribution_type << "   ";// << endl; 
  s << "number of bins     = " << setw(5) << sp.n_bins << "   ";// << endl; 
  s << "start of first bin = " << setw(5) << sp.start << "   ";// << endl; 
  s << "end of last bin    = " << setw(5) << sp.end << "   ";// << endl; 
  s << "bin width          = " << setw(5) << sp.step << "   ";// << endl; 
  s << "symmetry           = " << setw(5) << sp.symm << "   ";// << endl; 
  s << "particles          = ";
  for (int i = 0; i < (sp.particle).size(); i++){
    s << sp.particle[i];
    if (i < sp.particle.size() - 1){s << ", ";}
  }
  s << "   ";// << endl;
  if (sp.particle_2nd.size() != 0){
    s << "particles 2nd      = ";
    for (int i = 0; i < (sp.particle_2nd).size(); i++){
      s << sp.particle_2nd[i];
      if (i < sp.particle_2nd.size() - 1){s << ", ";}
    }
    s << endl;
  }
  return s;
  */
}



void xdistribution::initialization_distribution_bin(){
  static Logger logger("xdistribution::initialization_distribution_bin");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  distribution_no_bin_phasespace.clear();
  distribution_no_bin_no_qTcut_phasespace.clear();

  if (oset->n_ps == 1){
    logger << LOG_DEBUG_VERBOSE << xdistribution_name << "   CASE   oset->n_ps == 1" << endl;

    if (typeCumulative != CUMULATIVE_NONE){
      logger << LOG_DEBUG_VERBOSE << xdistribution_name << "   CASE   oset->n_ps == 1 && typeCumulative != CUMULATIVE_NONE" << endl;
      distribution_no_bin = bin; // could still be -1
      distribution_no_bin_max = bin_max; // could still be -1
      if (bin[0] != -1 && bin_max[0] != -1){
	distribution_no_bin_phasespace = vector<vector<int> > (1, vector<int> (1, 0));
	///	distribution_no_bin_phasespace.resize(1, vector<int> (1, 0)); // -> (1, 0) only 1 phasespace, namely 0; resize(1, ...) -> only 1 bin
	/////	    distribution_no_qTcut_phasespace.resize(1, vector<int> (1, 0)); // -> (1, 0) only 1 phasespace value, namely 0; resize(1, ...) -> only 1 qTcut value

	distribution_no_bin_no_qTcut_phasespace = vector<vector<vector<int> > > (1, vector<vector<int> > (1, vector<int> (1, 0)));
	///	distribution_no_bin_no_qTcut_phasespace.resize(1, vector<vector<int> > (1, vector<int> (1, 0)));
      }
    }
      
    else if (typeCumulative == CUMULATIVE_NONE){
      logger << LOG_DEBUG_VERBOSE << xdistribution_name << "   CASE   oset->n_ps == 1 && typeCumulative == CUMULATIVE_NONE" << endl;

      if (bin[0] != -1){
	distribution_no_bin = bin; // could still be -1
	// why this -> ???    distribution_no_bin.resize(1, bin[0]);
	// not needed:   distribution_no_bin_max = bin_max; // could still be -1
	distribution_no_bin_phasespace = vector<vector<int> > (1, vector<int> (1, 0));
	///	distribution_no_bin_phasespace.resize(1, vector<int> (1, 0)); // -> (1, 0) only 1 phasespace, namely 0; resize(1, ...) -> only 1 bin
	distribution_no_bin_no_qTcut_phasespace = vector<vector<vector<int> > > (1, vector<vector<int> > (1, vector<int> (1, 0)));
	///	distribution_no_bin_no_qTcut_phasespace.resize(1, vector<vector<int> > (1, vector<int> (1, 0)));
      }
    }

    else {
      logger << LOG_FATAL << "Wrong  typeCumulative  set!" << endl; exit(1);
    }
  }

  else if (oset->n_ps > 1){
    logger << LOG_DEBUG_VERBOSE << xdistribution_name << "   CASE   oset->n_ps > 1" << endl;
    
    if (typeCumulative != CUMULATIVE_NONE){
      logger << LOG_DEBUG_VERBOSE << xdistribution_name << "   CASE   oset->n_ps > 1 && typeCumulative != CUMULATIVE_NONE" << endl;

      // starting point: distribution_no_bin   contains a vector<int> (oset->n_ps) with the bin values for each phase-space as entries 
      distribution_no_bin = bin;
      distribution_no_bin_max = bin_max;
      
      sort(distribution_no_bin.begin(), distribution_no_bin.end());
      for (int i_b = oset->n_ps - 1; i_b > 0; i_b--){
	if (distribution_no_bin[i_b] == distribution_no_bin[i_b - 1]){distribution_no_bin.erase(distribution_no_bin.begin() + i_b);}
      }
      if (distribution_no_bin[0] == -1){distribution_no_bin.erase(distribution_no_bin.begin());}

      sort(distribution_no_bin_max.begin(), distribution_no_bin_max.end());
      for (int i_b = oset->n_ps - 1; i_b > 0; i_b--){
	if (distribution_no_bin_max[i_b] == distribution_no_bin_max[i_b - 1]){distribution_no_bin_max.erase(distribution_no_bin_max.begin() + i_b);}
      }
      if (distribution_no_bin_max[0] == -1){distribution_no_bin_max.erase(distribution_no_bin_max.begin());}

      distribution_no_bin_all = distribution_no_bin;
      for (int i_b = 0; i_b < distribution_no_bin_max.size(); i_b++){distribution_no_bin_all.push_back(distribution_no_bin_max[i_b]);}
      sort(distribution_no_bin_all.begin(), distribution_no_bin_all.end());
      for (int i_b = distribution_no_bin_all.size() - 1; i_b > 0; i_b--){
	if (distribution_no_bin_all[i_b] == distribution_no_bin_all[i_b - 1]){distribution_no_bin_all.erase(distribution_no_bin_all.begin() + i_b);}
      }
      
      stringstream out_dnb;
      for (int i_b = 0; i_b < distribution_no_bin.size(); i_b++){out_dnb << distribution_no_bin[i_b] << "   ";}
      logger << LOG_DEBUG << "out_dnb [" << xdistribution_name << "] = " << out_dnb.str() << endl;
      stringstream out_dnbx;
      for (int i_b = 0; i_b < distribution_no_bin_max.size(); i_b++){out_dnbx << distribution_no_bin_max[i_b] << "   ";}
      logger << LOG_DEBUG << "out_dnbx[" << xdistribution_name << "] = " << out_dnbx.str() << endl;
      stringstream out_dnba;
      for (int i_b = 0; i_b < distribution_no_bin_all.size(); i_b++){out_dnba << distribution_no_bin_all[i_b] << "   ";}
      logger << LOG_DEBUG << "out_dnba[" << xdistribution_name << "] = " << out_dnba.str() << endl;
      
      distribution_no_bin_phasespace.resize(distribution_no_bin_all.size());
      for (int i_a = 0; i_a < oset->n_ps; i_a++){
	for (int i_b = 0; i_b < distribution_no_bin_all.size(); i_b++){
	  // needs modification for multi-differential distributions:
	  if (distribution_no_bin_all[i_b] >= bin[i_a] && distribution_no_bin_all[i_b] < bin_max[i_a]){distribution_no_bin_phasespace[i_b].push_back(i_a);}
	  if (distribution_no_bin_all[i_b] == bin_max[i_a]){break;}
	}
      }
      
      for (int x_b = 0; x_b < distribution_no_bin.size(); x_b++){
	stringstream out_dnbp;
	for (int x_a = 0; x_a < distribution_no_bin_phasespace[x_b].size(); x_a++){
	  out_dnbp << distribution_no_bin_phasespace[x_b][x_a] << " ";
	}
	logger << LOG_DEBUG << "out_dnb_phasespace [" << xdistribution_name << "][" << x_b << "] @ " << setw(3) << distribution_no_bin_all[x_b] << " = " << out_dnbp.str() << endl;
      }
      
      if (distribution_no_bin_all.size() != 0){
	distribution_no_bin_no_qTcut_phasespace.resize(distribution_no_bin_all.size() - 1, vector<vector<int> > (oset->distribution_no_qTcut.size()));
	for (int i_b = 0; i_b < distribution_no_bin_all.size() - 1; i_b++){ // -1, because last entry contains no dipoles numbers any more
	  for (int i_q = 0; i_q < oset->distribution_no_qTcut.size(); i_q++){
	    for (int i_ba = 0; i_ba < distribution_no_bin_phasespace[i_b].size(); i_ba++){
	      for (int i_qa = 0; i_qa < oset->distribution_no_qTcut_phasespace[i_q].size(); i_qa++){
		if (distribution_no_bin_phasespace[i_b][i_ba] == oset->distribution_no_qTcut_phasespace[i_q][i_qa]){
		  distribution_no_bin_no_qTcut_phasespace[i_b][i_q].push_back(distribution_no_bin_phasespace[i_b][i_ba]);
		}
	      }
	    }
	  }
	}
      }
    }

    else if (typeCumulative == CUMULATIVE_NONE){
      logger << LOG_DEBUG_VERBOSE << xdistribution_name << "   CASE   oset->n_ps > 1 && typeCumulative == CUMULATIVE_NONE" << endl;

      distribution_no_bin = bin;
      sort(distribution_no_bin.begin(), distribution_no_bin.end());

      // check why n_ps - 1 !!!
      for (int i_b = oset->n_ps - 1; i_b > 0; i_b--){
	if (distribution_no_bin[i_b] == distribution_no_bin[i_b - 1]){distribution_no_bin.erase(distribution_no_bin.begin() + i_b);}
      }
      if (distribution_no_bin[0] == -1){distribution_no_bin.erase(distribution_no_bin.begin());}
      
      distribution_no_bin_phasespace.resize(distribution_no_bin.size());
      for (int i_a = 0; i_a < oset->n_ps; i_a++){
	for (int i_b = 0; i_b < distribution_no_bin.size(); i_b++){
	  if (distribution_no_bin[i_b] == bin[i_a]){distribution_no_bin_phasespace[i_b].push_back(i_a); break;}
	}
      }

      distribution_no_bin_no_qTcut_phasespace.resize(distribution_no_bin.size(), vector<vector<int> > (oset->distribution_no_qTcut.size()));
      for (int i_b = 0; i_b < distribution_no_bin.size(); i_b++){
	for (int i_q = 0; i_q < oset->distribution_no_qTcut.size(); i_q++){
	  for (int i_ba = 0; i_ba < distribution_no_bin_phasespace[i_b].size(); i_ba++){
	    for (int i_qa = 0; i_qa < oset->distribution_no_qTcut_phasespace[i_q].size(); i_qa++){
	      if (distribution_no_bin_phasespace[i_b][i_ba] == oset->distribution_no_qTcut_phasespace[i_q][i_qa]){
		distribution_no_bin_no_qTcut_phasespace[i_b][i_q].push_back(distribution_no_bin_phasespace[i_b][i_ba]);
	      }
	    }
	  }
	}
      }
      
    }
  }





  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
