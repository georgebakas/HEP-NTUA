#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"

// Eq. (3.4) without E_gamma prefactor
double frixione_discr(double delta, observable_set & oset) {
  if (osi_frixione_n == 1){return osi_frixione_epsilon * (1. - cos(delta)) / (1. - cos(osi_frixione_delta_0));}
  else {return osi_frixione_epsilon * pow((1. - cos(delta)) / (1. - cos(osi_frixione_delta_0)), osi_frixione_n);}
}

double frixione_discr(double delta, double delta_dynamic, observable_set & oset) {
  if (osi_frixione_n == 1){return osi_frixione_epsilon * (1. - cos(delta)) / (1. - cos(delta_dynamic));}
  else {return osi_frixione_epsilon * pow((1. - cos(delta)) / (1. - cos(delta_dynamic)), osi_frixione_n);}
}

double frixione_discriminant_R2(double delta, double delta_dynamic, observable_set & oset) {
  //replace [(1- cos R)/(1- cos R_0)]^n  =>  [(R^2/R_0^2)]^n  in Frixione's isolation formula (this matters when R>1, i.e. at small pT)  
  if (osi_frixione_n == 1){return osi_frixione_epsilon * (pow(delta, 2)) / (pow(delta_dynamic, 2));}
  else {return osi_frixione_epsilon * pow((pow(delta, 2)) / (pow(delta_dynamic, 2)), osi_frixione_n);}
}

struct frix_parton {
  int index;
  double delta;
  bool operator<(const frix_parton& rhs) const{return delta<rhs.delta;}
};

void photon_recombination(vector<int> & no_unrecombined_photon, int i_a, observable_set & oset){
  static Logger logger("photon_recombination");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi_ps_runtime_photon[i_a].size() == 0){return;} // obviously, no recombination needed in that case !

  vector<int> no_rec_photon;
  vector<int> no_rec_charged_particle;
  no_rec_charged_particle = osi_ps_runtime_photon_recombination[i_a];
  no_rec_photon = osi_ps_runtime_photon[i_a];

  /*
  for (int i_p = 3; i_p < osi_p_parton[i_a].size(); i_p++){
    if ((osi_type_parton[i_a][i_p] > 0 && osi_type_parton[i_a][i_p] < 7) || 
	(osi_type_parton[i_a][i_p] < 0 && osi_type_parton[i_a][i_p] > -7) || 
	osi_type_parton[i_a][i_p] == 11 || osi_type_parton[i_a][i_p] == -11 || 
	osi_type_parton[i_a][i_p] == 13 || osi_type_parton[i_a][i_p] == -13 || 
	osi_type_parton[i_a][i_p] == 15 || osi_type_parton[i_a][i_p] == -15){no_rec_charged_particle.push_back(i_p);}
    // || osi_type_parton[i_a][i_p] == 24 || osi_type_parton[i_a][i_p] == -24
    else if (osi_type_parton[i_a][i_p] == 22){no_rec_photon.push_back(i_p);}
  }

  if (no_rec_photon.size() == 0){return;} // obviously, no recombination needed in that case !
  */

  logger << LOG_DEBUG_VERBOSE << "no_rec_charged_particle.size() = " << no_rec_charged_particle.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "no_rec_photon.size()           = " << no_rec_photon.size() << endl;

  for (int i_n = 0; i_n < no_rec_charged_particle.size(); i_n++){
    logger << LOG_DEBUG_VERBOSE << "i_a = " << setw(2) << i_a << "   no_rec_charged_particle[" << i_n << "] = " << no_rec_charged_particle[i_n] << endl;
  }
  for (int i_n = 0; i_n < no_rec_photon.size(); i_n++){
    logger << LOG_DEBUG_VERBOSE << "i_a = " << setw(2) << i_a << "   no_rec_photon[" << i_n << "] = " << no_rec_photon[i_n] << endl;
  }
  
  vector<vector<double> > distance(no_rec_photon.size(), vector<double> (no_rec_charged_particle.size(), 0.));
  for (int i_r = 0; i_r < no_rec_photon.size(); i_r++){
    for (int i_c = 0; i_c < no_rec_charged_particle.size(); i_c++){
      if (osi_photon_R_definition == 0){distance[i_r][i_c] = R2_eta(osi_particle_event[0][i_a][no_rec_photon[i_r]], osi_particle_event[0][i_a][no_rec_charged_particle[i_c]]);}
      else if (osi_photon_R_definition == 1){distance[i_r][i_c] = R2_rapidity(osi_particle_event[0][i_a][no_rec_photon[i_r]], osi_particle_event[0][i_a][no_rec_charged_particle[i_c]]);}
      logger << LOG_DEBUG_VERBOSE << "distance[" << i_r << "][" << i_c << "] = " << distance[i_r][i_c] << endl;
    }
  }

  while (no_rec_photon.size() > 0){
    int min1_rec_photon = 0;
    int min1_rec_charged_particle = 0;
    double min1_distance = 1.e99;

    for (int i_r = 0; i_r < no_rec_photon.size(); i_r++){
      logger << LOG_DEBUG_VERBOSE << "phot.: [0][" << i_a << "][" << no_rec_photon[i_r] << "] = " << osi_particle_event[0][i_a][no_rec_photon[i_r]].momentum << endl;
      for (int i_c = 0; i_c < no_rec_charged_particle.size(); i_c++){
	logger << LOG_DEBUG_VERBOSE << "ch.p.: [0][" << i_a << "][" << no_rec_charged_particle[i_c] << "] = " << osi_particle_event[0][i_a][no_rec_charged_particle[i_c]].momentum << "   distance[" << i_r << "][" << i_c << "] = " << setprecision(5) << distance[i_r][i_c] << " photon_R2 = " << osi_photon_R2 << endl;
	if (distance[i_r][i_c] < min1_distance){
	  min1_rec_photon = i_r;
	  min1_rec_charged_particle = i_c;
	  min1_distance = distance[i_r][i_c];

	}
      }
    }
    if (min1_distance < osi_photon_R2){
      osi_particle_event[0][i_a][no_rec_charged_particle[min1_rec_charged_particle]] = osi_particle_event[0][i_a][no_rec_charged_particle[min1_rec_charged_particle]] + osi_particle_event[0][i_a][no_rec_photon[min1_rec_photon]];
      osi_particle_event[0][i_a][no_rec_photon[min1_rec_photon]] = particle();

      logger << LOG_DEBUG_VERBOSE << "osi_particle_event[0][" << i_a << "][" << no_rec_charged_particle[min1_rec_charged_particle] << "] = " << osi_particle_event[0][i_a][no_rec_charged_particle[min1_rec_charged_particle]].momentum << endl;
	    //      no_unrecombined_photon.push_back(no_rec_photon[min1_rec_photon]); // why should this happen ???
      no_rec_photon.erase(no_rec_photon.begin() + min1_rec_photon, no_rec_photon.begin() + min1_rec_photon + 1);
      distance.erase(distance.begin() + min1_rec_photon, distance.begin() + min1_rec_photon + 1);
      logger << LOG_DEBUG_VERBOSE << "Photon " << min1_rec_charged_particle << " (" << no_rec_charged_particle[min1_rec_charged_particle] << ") recombined." << endl;
      for (int i_r = 0; i_r < no_rec_photon.size(); i_r++){
	if (osi_photon_R_definition == 0){distance[i_r][min1_rec_charged_particle] = R2_eta(osi_particle_event[0][i_a][no_rec_photon[i_r]], osi_particle_event[0][i_a][no_rec_charged_particle[min1_rec_charged_particle]]);}
	else if (osi_photon_R_definition == 1){distance[i_r][min1_rec_charged_particle] = R2_rapidity(osi_particle_event[0][i_a][no_rec_photon[i_r]], osi_particle_event[0][i_a][no_rec_charged_particle[min1_rec_charged_particle]]);}
	//      logger << LOG_DEBUG_VERBOSE << "Photon " << min1_rec_charged_particle << " (" << no_rec_charged_particle[min1_rec_charged_particle] << ") recombined." << endl;
	//	logger << LOG_DEBUG_VERBOSE << "new_distance[" << i_r << "][" << min1_rec_charged_particle << "] = " << distance[i_r][min1_rec_charged_particle] << endl;
		//
      //	distance[i_r][min1_rec_charged_particle] = R2_rapidity(osi_particle_event[0][i_a][no_rec_photon[i_r]], osi_particle_event[0][i_a][no_rec_charged_particle[min1_rec_charged_particle]]);
      }
    }
    else {
      //      logger << LOG_DEBUG_VERBOSE << "No combination done!" << endl;
      for (int i_r = 0; i_r < no_rec_photon.size(); i_r++){
	no_unrecombined_photon.push_back(no_rec_photon[i_r]);
	//	osi_particle_event[0][i_a][no_rec_photon[i_r]] = fourvector();
	//	logger << LOG_DEBUG_VERBOSE << "no_rec_photon[i_r] = " << no_rec_photon[i_r] << endl;
	no_rec_photon.erase(no_rec_photon.begin() + i_r);

	//	logger << LOG_DEBUG_VERBOSE << "osi_particle_event[0][i_a][no_rec_photon[i_r]] = " << osi_particle_event[0][i_a][no_rec_photon[i_r]].momentum() << endl;
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "no_rec_photon.size() = " << no_rec_photon.size() << endl;
  for (int i_r = 0; i_r < no_rec_photon.size(); i_r++){
    logger << LOG_DEBUG_VERBOSE << "no_rec_photon[i_r] = " << no_rec_photon[i_r] << endl;
  }
  logger << LOG_DEBUG_VERBOSE << "no_unrecombined_photon.size() = " << no_unrecombined_photon.size() << endl;
  for (int i_r = 0; i_r < no_unrecombined_photon.size(); i_r++){
    logger << LOG_DEBUG_VERBOSE << "no_unrecombined_photon[i_r] = " << no_unrecombined_photon[i_r] << endl;
  }
  logger << LOG_DEBUG_VERBOSE << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

// Frixione photon isolation routine with updated interface. 
// old: Adds 'momentum' to vector 'isolated_photon' and 1 to number_photon
void frixione_isolation(int & number_photon, vector<particle> & isolated_photon, particle & photon, vector<particle> & protojet, observable_set & oset){
  static Logger logger("frixione_isolation");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  vector<frix_parton> frixione_order_parton;
  //  frix_parton tmp_parton;
  static vector<frix_parton> frixione_order_parton;
  frixione_order_parton.clear();
  static frix_parton tmp_parton;

  //  static double frixione_delta2_0 = pow(osi_frixione_delta_0,2);
  static double frixione_delta2_0 = pow(osi_frixione_delta_0, 2);
  static double delta_dynamic = 0.;
  if (osi_frixione_isolation == 3 || osi_frixione_isolation == 4 || osi_frixione_isolation == 5){
    delta_dynamic = osi_msi.M_Z / photon.pT / sqrt(osi_frixione_epsilon); // FIXME -> M_Z as a fixed value by now
    frixione_delta2_0 = pow(delta_dynamic, 2);
  }
  if (osi_frixione_isolation == 6){
    delta_dynamic = osi_msi.M_Z / photon.pT / sqrt(osi_frixione_epsilon) * sqrt(2.); // FIXME -> M_Z as a fixed value by now
    frixione_delta2_0 = pow(delta_dynamic, 2);
  }
  if (osi_frixione_isolation == 7){
    delta_dynamic = osi_msi.M_Z / photon.pT / sqrt(osi_frixione_epsilon) * 2.; // FIXME -> M_Z as a fixed value by now
    frixione_delta2_0 = pow(delta_dynamic, 2);
  }
  if (osi_frixione_isolation == 8){
    delta_dynamic = osi_msi.M_Z / photon.pT / sqrt(osi_frixione_epsilon) * 2. * sqrt(2.); // FIXME -> M_Z as a fixed value by now
    frixione_delta2_0 = pow(delta_dynamic, 2);
  }
  if (osi_frixione_isolation == 9){
    delta_dynamic = osi_msi.M_Z / photon.pT / sqrt(osi_frixione_epsilon) * 4.; // FIXME -> M_Z as a fixed value by now
    frixione_delta2_0 = pow(delta_dynamic, 2);
  }
  if (osi_frixione_isolation == 10){
    if (photon.pT < 2 * osi_msi.M_Z){delta_dynamic = osi_frixione_delta_0;}
    else {delta_dynamic = osi_msi.M_Z / photon.pT;}
    frixione_delta2_0 = pow(delta_dynamic, 2);
  }
  if (osi_frixione_isolation == 11){
    delta_dynamic = osi_msi.M_Z / photon.pT / sqrt(osi_frixione_epsilon);
    if (delta_dynamic > osi_frixione_delta_0){delta_dynamic = osi_frixione_delta_0;}
    frixione_delta2_0 = pow(delta_dynamic, 2);
  }

  //    cout << "delta_dynamic = " << delta_dynamic << "   osi_msi.M_Z = " << osi_msi.M_Z << "   photon.pT = " << photon.pT << endl;

  // identify partons inside the cone
  logger << LOG_DEBUG_VERBOSE << "photon = " << photon.momentum << endl;
  for (int i = 0; i < protojet.size(); i++){
    //    logger << LOG_DEBUG_VERBOSE << "protojet[" << i << "] = " << protojet[i].momentum << endl;
    //    logger << LOG_DEBUG_VERBOSE << "delta_2 = " << delta_2 << endl;
    //    logger << LOG_DEBUG_VERBOSE << "frixione_delta2_0 = " << frixione_delta2_0 << endl;
    double delta_2 = 0.;
    if (osi_frixione_isolation == 5 || 
	osi_frixione_isolation == 6 || 
	osi_frixione_isolation == 7 || 
	osi_frixione_isolation == 8 || 
	osi_frixione_isolation == 9){delta_2 = R2_coshrapidity_cosphi(photon.momentum, protojet[i].momentum);}
    else {delta_2 = R2_eta(photon, protojet[i]);}

    if (delta_2 < frixione_delta2_0) {
      tmp_parton.index = i;
      tmp_parton.delta = sqrt(delta_2);//sqrt(R2_eta(photon, protojet[i]));
      frixione_order_parton.push_back(tmp_parton);
    }
  }
  sort(frixione_order_parton.begin(), frixione_order_parton.end());
  logger << LOG_DEBUG_VERBOSE << "frixione_order_parton.size() = " << frixione_order_parton.size() << endl;

  int isolation = 0;
  double E_i_sum = 0.;
  for (int i = 0; i < frixione_order_parton.size(); i++){
    double H_i=0.0;
    if      (osi_frixione_isolation == 1){H_i = photon.pT * frixione_discr(frixione_order_parton[i].delta, oset);} // standard definition
    else if (osi_frixione_isolation == 2){H_i = 10 * frixione_discr(frixione_order_parton[i].delta, oset);} // FIXME -> 10GeV as a fixed limit by now
    else if (osi_frixione_isolation == 3){H_i = photon.pT * frixione_discr(frixione_order_parton[i].delta, delta_dynamic, oset);} // dynamic cone (photon.pT dependent) 
    else if (osi_frixione_isolation == 4){H_i = photon.pT * frixione_discriminant_R2(frixione_order_parton[i].delta, delta_dynamic, oset);} // delta definition via (R²/R_0²)^n !!!
    else if (osi_frixione_isolation == 5){H_i = photon.pT * frixione_discriminant_R2(frixione_order_parton[i].delta, delta_dynamic, oset);} // delta definition via (R²/R_0²)^n !!!
    else if (osi_frixione_isolation == 6){H_i = photon.pT * frixione_discriminant_R2(frixione_order_parton[i].delta, delta_dynamic, oset);} // delta definition via (R²/R_0²)^n !!!
    else if (osi_frixione_isolation == 7){H_i = photon.pT * frixione_discriminant_R2(frixione_order_parton[i].delta, delta_dynamic, oset);} // delta definition via (R²/R_0²)^n !!!
    else if (osi_frixione_isolation == 8){H_i = photon.pT * frixione_discriminant_R2(frixione_order_parton[i].delta, delta_dynamic, oset);} // delta definition via (R²/R_0²)^n !!!
    else if (osi_frixione_isolation == 9){H_i = photon.pT * frixione_discriminant_R2(frixione_order_parton[i].delta, delta_dynamic, oset);} // delta definition via (R²/R_0²)^n !!!
    else if (osi_frixione_isolation == 10){H_i = photon.pT * frixione_discr(frixione_order_parton[i].delta, delta_dynamic, oset);} // dynamic cone (photon.pT dependent) 
    else if (osi_frixione_isolation == 11){H_i = photon.pT * frixione_discr(frixione_order_parton[i].delta, delta_dynamic, oset);} // dynamic cone (photon.pT dependent) 
    else {assert(false);}
   
    E_i_sum += protojet[frixione_order_parton[i].index].pT;
    if (E_i_sum > H_i){
      /*
      if (i == protojet.size() - 1){
	cout << "WARNING !!!   " << "i = " << i << "   R_ia = " << setprecision(8) << setw(16) << frixione_order_parton[i].delta << "   R_0 = " << setprecision(8) << setw(16) << delta_dynamic << "   pTA = " << setprecision(8) << setw(16) << photon.pT << "   discr = " << setprecision(8) << setw(16) << H_i / photon.pT << endl;
	cout << "photon       = " << photon.momentum << "   photon.pT = " << photon.pT << endl;
	cout << "frixione_order_parton.size() = " << frixione_order_parton.size() << endl;
	for (int j = 0; j < frixione_order_parton.size(); j++){
	  cout << "protojet[" << frixione_order_parton[j].index << "] = " << protojet[frixione_order_parton[j].index].momentum << "   protojet[" << frixione_order_parton[j].index << "].pT = " << protojet[frixione_order_parton[j].index].pT << "   frixione_order_parton[" << j << "].delta = " << frixione_order_parton[j].delta << endl;
	}
	cout << endl;
      }
      */
      isolation = -1; break;
    }
    
  }
  if (isolation == 0){
    number_photon++;
    isolated_photon.push_back(photon);
    logger << LOG_DEBUG_VERBOSE << "number_photon = " << number_photon << endl;
    logger << LOG_DEBUG_VERBOSE << "isolated_photon = " << isolated_photon[isolated_photon.size() - 1] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
