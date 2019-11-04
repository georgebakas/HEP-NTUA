#include "../include/classes.cxx"
#include "../include/definitions.observable.set.cxx"
void jet_algorithm_flavour(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin, observable_set & oset){
  static Logger logger("jet_algorithm_flavour");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  logger << LOG_DEBUG_VERBOSE << "osi_jet_algorithm = " << osi_jet_algorithm << endl;
  if      (osi_jet_algorithm == 0){
    jet.swap(protojet);
    jet_flavour.swap(protojet_flavour);
    jet_parton_origin.swap(protojet_parton_origin);
    /*
    jet = protojet;
    jet_flavour = protojet_flavour;
    jet_parton_origin = protojet_parton_origin;
    protojet.erase(protojet.begin(), protojet.end());
    protojet_flavour.erase(protojet_flavour.begin(), protojet_flavour.end());
    protojet_parton_origin.erase(protojet_parton_origin.begin(), protojet_parton_origin.end());
    */
  }
  else if (osi_jet_algorithm == 1){sc_Ellis_Soper_flavour(protojet, protojet_flavour, protojet_parton_origin, jet, jet_flavour, jet_parton_origin, oset);}
  else if (osi_jet_algorithm == 2){kTrun2_flavour(protojet, protojet_flavour, protojet_parton_origin, jet, jet_flavour, jet_parton_origin, oset);}
  else if (osi_jet_algorithm == 3){antikT_flavour(protojet, protojet_flavour, protojet_parton_origin, jet, jet_flavour, jet_parton_origin, oset);}
  else if (osi_jet_algorithm == 4){sc_Ellis_Soper_mod_flavour(protojet, protojet_flavour, protojet_parton_origin, jet, jet_flavour, jet_parton_origin, oset);}
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void sc_Ellis_Soper_flavour(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin, observable_set & oset){
  int min_ot = 1;
  int min1 = 0;
  int min2 = 0;
  double d_min1 = 1.e99;//, delta_phi;
  //  int min_ot, min1, min2;
  //  double d_min;//, delta_phi;
  vector<vector<double> > d_two(protojet.size(), vector<double> (protojet.size()));
  while (protojet.size() > 0){
    d_min1 = 1.e99;
    min_ot = 1;
    for (int i = 0; i < protojet.size(); i++){
      if (protojet[i].pT2 < d_min1){d_min1 = protojet[i].pT2; min1 = i;}
    }
    for (int i = 0; i < protojet.size(); i++){
      for (int j = i + 1; j < protojet.size(); j++){
	if (protojet[i].pT2 < protojet[j].pT2){d_two[i][j] = protojet[i].pT2;}
	else{d_two[i][j] = protojet[j].pT2;}
	if (osi_jet_R_definition == 0){d_two[i][j] = d_two[i][j] * R2_eta(protojet[i], protojet[j]) / osi_jet_R2;}
	else if (osi_jet_R_definition == 1){d_two[i][j] = d_two[i][j] * R2_rapidity(protojet[i], protojet[j]) / osi_jet_R2;}
	if (d_two[i][j] < d_min1){d_min1 = d_two[i][j]; min1 = i; min2 = j; min_ot = 2;}
      }
    }
    if (min_ot == 1){
      transfer_protojet_to_jet(protojet, protojet_flavour, protojet_parton_origin, min1, jet, jet_flavour, jet_parton_origin);
    }
    else {
      double phi_i = protojet[min1].phi;
      double phi_j = protojet[min2].phi;
      double ET_i = protojet[min1].pT;
      double ET_j = protojet[min2].pT;
      double ET = ET_i + ET_j;
      double deltaphi = abs(phi_i - phi_j);
      double phi;
      if (deltaphi <= pi){phi = (ET_i * phi_i + ET_j * phi_j) / ET;}
      else {
	if (phi_i > pi){phi_j = phi_j + f2pi;}
	else           {phi_i = phi_i + f2pi;}
	phi = (ET_i * phi_i + ET_j * phi_j) / ET;
	if (phi > f2pi){phi -= f2pi;}
      }
      /*
      double eta = (ET_i * protojet[min1].eta + ET_j * protojet[min2].eta) / ET;
      double pz = ET / tan(2 * atan(exp(-eta)));
      particle new_protojet(fourvector (sqrt(pow(ET, 2) + pow(pz, 2)), cos(phi) * ET, sin(phi) * ET, pz));
      combine_protojet(protojet, protojet_flavour, protojet_parton_origin, min1, min2, new_protojet);
      */
      combine_protojet(protojet, protojet_flavour, protojet_parton_origin, min1, min2);
      // not the correct combination procedure !!!
     }
  }
}

void sc_Ellis_Soper_mod_flavour(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin, observable_set & oset){
  //  int min_ot, min1, min2;
  //  double d_min1;//, delta_phi;
  int min_ot = 1;
  int min1 = 0;
  int min2 = 0;
  double d_min1 = 1.e99;//, delta_phi;
  vector<vector<double> > d_two(protojet.size(), vector<double> (protojet.size()));
  while (protojet.size() > 0){
    d_min1 = 1.e99;
    min_ot = 1;
    for (int i = 0; i < protojet.size(); i++){
      if (protojet[i].pT2 < d_min1){d_min1 = protojet[i].pT2; min1 = i;}
    }
    for (int i = 0; i < protojet.size(); i++){
      for (int j = i + 1; j < protojet.size(); j++){
	if (protojet[i].pT2 < protojet[j].pT2){d_two[i][j] = protojet[i].pT2;}
	else{d_two[i][j] = protojet[j].pT2;}
	if (osi_jet_R_definition == 0){d_two[i][j] = d_two[i][j] * R2_eta(protojet[i], protojet[j]) / osi_jet_R2;}
	else if (osi_jet_R_definition == 1){d_two[i][j] = d_two[i][j] * R2_rapidity(protojet[i], protojet[j]) / osi_jet_R2;}
	if (d_two[i][j] < d_min1){d_min1 = d_two[i][j]; min1 = i; min2 = j; min_ot = 2;}
      }
    }
    if (min_ot == 1){
      transfer_protojet_to_jet(protojet, protojet_flavour, protojet_parton_origin, min1, jet, jet_flavour, jet_parton_origin);
    }
    else {
      //, new_protojet      particle new_protojet = protojet[min1] + protojet[min2];
      combine_protojet(protojet, protojet_flavour, protojet_parton_origin, min1, min2);
    }
  }
}



void kTrun2_flavour(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin, observable_set & oset){
  int min_ot = 1;
  int min1 = 0;
  int min2 = 0;
  double d_min1 = 1.e99;//, delta_phi;
  vector<vector<double> > d_two(protojet.size(), vector<double> (protojet.size()));
  while (protojet.size() > 0){
    d_min1 = 1.e99;
    min_ot = 1;
    for (int i = 0; i < protojet.size(); i++){
      if (protojet[i].pT2 < d_min1){d_min1 = protojet[i].pT2; min1 = i;}
    }
    for (int i = 0; i < protojet.size(); i++){
      for (int j = i + 1; j < protojet.size(); j++){
	if (protojet[i].pT2 < protojet[j].pT2){d_two[i][j] = protojet[i].pT2;}
	else{d_two[i][j] = protojet[j].pT2;}
	if (osi_jet_R_definition == 0){d_two[i][j] = d_two[i][j] * R2_eta(protojet[i], protojet[j]) / osi_jet_R2;}
	else if (osi_jet_R_definition == 1){d_two[i][j] = d_two[i][j] * R2_rapidity(protojet[i], protojet[j]) / osi_jet_R2;}
	if (d_two[i][j] < d_min1){d_min1 = d_two[i][j]; min1 = i; min2 = j; min_ot = 2;}
      }
    }
    if (min_ot == 1){
      transfer_protojet_to_jet(protojet, protojet_flavour, protojet_parton_origin, min1, jet, jet_flavour, jet_parton_origin);
    }
    else{
      //, new_protojet      particle new_protojet = protojet[min1] + protojet[min2];
      combine_protojet(protojet, protojet_flavour, protojet_parton_origin, min1, min2);
    }
  }
}

void antikT_flavour(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin, observable_set & oset){
  static Logger logger("antikT_flavour");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  int min_ot = 1;
  int min1 = 0;
  int min2 = 0;
  double d_min1 = 1.e99;//, delta_phi;
  vector<vector<double> > d_two(protojet.size(), vector<double> (protojet.size()));
  while (protojet.size() > 0){
    d_min1 = 1.e99;
    min_ot = 1;
    min1 = 0;
    for (int i = 0; i < protojet.size(); i++){
      if (1. / protojet[i].pT2 < d_min1){d_min1 = 1. / protojet[i].pT2; min1 = i;}
    }
    for (int i = 0; i < protojet.size(); i++){
      for (int j = i + 1; j < protojet.size(); j++){
	if (protojet[i].pT2 > protojet[j].pT2){d_two[i][j] = 1. / protojet[i].pT2;}
	else {d_two[i][j] = 1. / protojet[j].pT2;}
	if (osi_jet_R_definition == 0){d_two[i][j] = d_two[i][j] * R2_eta(protojet[i], protojet[j]) / osi_jet_R2;}
	else if (osi_jet_R_definition == 1){d_two[i][j] = d_two[i][j] * R2_rapidity(protojet[i], protojet[j]) / osi_jet_R2;}
	if (d_two[i][j] < d_min1){d_min1 = d_two[i][j]; min1 = i; min2 = j; min_ot = 2;}
      }
    }
    if (min_ot == 1){
      transfer_protojet_to_jet(protojet, protojet_flavour, protojet_parton_origin, min1, jet, jet_flavour, jet_parton_origin);
    }
    else {
      //, new_protojet      particle new_protojet = protojet[min1] + protojet[min2];
      combine_protojet(protojet, protojet_flavour, protojet_parton_origin, min1, min2);
    }
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void transfer_protojet_to_jet(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, int pj1, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin){
  jet.push_back(protojet[pj1]);
  jet_flavour.push_back(protojet_flavour[pj1]);
  jet_parton_origin.push_back(protojet_parton_origin[pj1]);
  protojet.erase(protojet.begin() + pj1);
  protojet_flavour.erase(protojet_flavour.begin() + pj1);
  protojet_parton_origin.erase(protojet_parton_origin.begin() + pj1);
}

void combine_protojet(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, int pj1, int pj2){
  //  protojet.push_back(new_protojet); //, particle & new_protojet
  protojet.push_back(protojet[pj1] + protojet[pj2]);
  vector<int> temp_flavour(0);
  for (int i_p = 0; i_p < protojet_flavour[pj1].size(); i_p++){temp_flavour.push_back(protojet_flavour[pj1][i_p]);}
  for (int i_p = 0; i_p < protojet_flavour[pj2].size(); i_p++){temp_flavour.push_back(protojet_flavour[pj2][i_p]);}
  
  vector<int> temp_parton_origin(0);
  for (int i_p = 0; i_p < protojet_parton_origin[pj1].size(); i_p++){temp_parton_origin.push_back(protojet_parton_origin[pj1][i_p]);}
  for (int i_p = 0; i_p < protojet_parton_origin[pj2].size(); i_p++){temp_parton_origin.push_back(protojet_parton_origin[pj2][i_p]);}
  
  protojet_flavour.push_back(temp_flavour);
  protojet_parton_origin.push_back(temp_parton_origin);

  protojet.erase(protojet.begin() + pj2);
  protojet.erase(protojet.begin() + pj1);
  
  protojet_flavour.erase(protojet_flavour.begin() + pj2);
  protojet_flavour.erase(protojet_flavour.begin() + pj1);

  protojet_parton_origin.erase(protojet_parton_origin.begin() + pj2);
  protojet_parton_origin.erase(protojet_parton_origin.begin() + pj1);
}





