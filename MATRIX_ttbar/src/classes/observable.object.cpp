#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

void observable_set::determine_phasespace_object_partonlevel(){
  static  Logger logger("observable_set::determine_phasespace_object_partonlevel");
  logger << LOG_DEBUG << "called" << endl;

  for (int i_a = 0; i_a < n_ps; i_a++){
    logger << LOG_DEBUG << "csi->type_parton[" << i_a << "].size() = " << csi->type_parton[i_a].size() << endl;
    for (int i_p = 3; i_p < csi->type_parton[i_a].size(); i_p++){
      logger << LOG_DEBUG << "i_a = " << i_a << "   i_p = " << i_p << endl;
      if (abs(csi->type_parton[i_a][i_p]) == 22){ps_n_partonlevel[i_a][esi.observed_object["photon"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 11 || abs(csi->type_parton[i_a][i_p]) == 13 || abs(csi->type_parton[i_a][i_p]) == 15){ps_n_partonlevel[i_a][esi.observed_object["lep"]]++;}
      if (csi->type_parton[i_a][i_p] == 11 || csi->type_parton[i_a][i_p] == 13 || csi->type_parton[i_a][i_p] == 15){ps_n_partonlevel[i_a][esi.observed_object["lm"]]++;}
      if (csi->type_parton[i_a][i_p] == -11 || csi->type_parton[i_a][i_p] == -13 || csi->type_parton[i_a][i_p] == -15){ps_n_partonlevel[i_a][esi.observed_object["lp"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 11){ps_n_partonlevel[i_a][esi.observed_object["e"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 13){ps_n_partonlevel[i_a][esi.observed_object["mu"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 15){ps_n_partonlevel[i_a][esi.observed_object["tau"]]++;}
      if (csi->type_parton[i_a][i_p] == 11){ps_n_partonlevel[i_a][esi.observed_object["em"]]++;}
      if (csi->type_parton[i_a][i_p] == 13){ps_n_partonlevel[i_a][esi.observed_object["mum"]]++;}
      if (csi->type_parton[i_a][i_p] == 15){ps_n_partonlevel[i_a][esi.observed_object["taum"]]++;}
      if (csi->type_parton[i_a][i_p] == -11){ps_n_partonlevel[i_a][esi.observed_object["ep"]]++;}
      if (csi->type_parton[i_a][i_p] == -13){ps_n_partonlevel[i_a][esi.observed_object["mup"]]++;}
      if (csi->type_parton[i_a][i_p] == -15){ps_n_partonlevel[i_a][esi.observed_object["taup"]]++;}
      //	if (abs(csi->type_parton[i_a][i_p]) == 12 || abs(csi->type_parton[i_a][i_p]) == 14 || abs(csi->type_parton[i_a][i_p]) == 16){ps_n_partonlevel[i_a][esi.observed_object["nu"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 12 || abs(csi->type_parton[i_a][i_p]) == 14 || abs(csi->type_parton[i_a][i_p]) == 16){ps_n_partonlevel[i_a][esi.observed_object["nua"]]++;}
      if (csi->type_parton[i_a][i_p] == 12 || csi->type_parton[i_a][i_p] == 14 || csi->type_parton[i_a][i_p] == 16){ps_n_partonlevel[i_a][esi.observed_object["nu"]]++;}
      if (csi->type_parton[i_a][i_p] == -12 || csi->type_parton[i_a][i_p] == -14 || csi->type_parton[i_a][i_p] == -16){ps_n_partonlevel[i_a][esi.observed_object["nux"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 12){ps_n_partonlevel[i_a][esi.observed_object["nea"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 14){ps_n_partonlevel[i_a][esi.observed_object["nma"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 16){ps_n_partonlevel[i_a][esi.observed_object["nta"]]++;}
      if (csi->type_parton[i_a][i_p] == 12){ps_n_partonlevel[i_a][esi.observed_object["ne"]]++;}
      if (csi->type_parton[i_a][i_p] == 14){ps_n_partonlevel[i_a][esi.observed_object["nm"]]++;}
      if (csi->type_parton[i_a][i_p] == 16){ps_n_partonlevel[i_a][esi.observed_object["nt"]]++;}
      if (csi->type_parton[i_a][i_p] == -12){ps_n_partonlevel[i_a][esi.observed_object["nex"]]++;}
      if (csi->type_parton[i_a][i_p] == -14){ps_n_partonlevel[i_a][esi.observed_object["nmx"]]++;}
      if (csi->type_parton[i_a][i_p] == -16){ps_n_partonlevel[i_a][esi.observed_object["ntx"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 12 || abs(csi->type_parton[i_a][i_p]) == 14 || abs(csi->type_parton[i_a][i_p]) == 16){ps_n_partonlevel[i_a][esi.observed_object["missing"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 12 || abs(csi->type_parton[i_a][i_p]) == 14 || abs(csi->type_parton[i_a][i_p]) == 16){esi.n_parton_nu++;}
      if (abs(csi->type_parton[i_a][i_p]) < 5){ps_n_partonlevel[i_a][esi.observed_object["ljet"]]++;}
      if (csi->type_parton[i_a][i_p] == 5){ps_n_partonlevel[i_a][esi.observed_object["bjet_b"]]++;}
      if (csi->type_parton[i_a][i_p] == -5){ps_n_partonlevel[i_a][esi.observed_object["bjet_bx"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 5){ps_n_partonlevel[i_a][esi.observed_object["bjet"]]++;}
      if (csi->type_parton[i_a][i_p] == 6){ps_n_partonlevel[i_a][esi.observed_object["top"]]++;}
      if (csi->type_parton[i_a][i_p] == -6){ps_n_partonlevel[i_a][esi.observed_object["atop"]]++;}
      if (abs(csi->type_parton[i_a][i_p]) == 6){ps_n_partonlevel[i_a][esi.observed_object["tjet"]]++;}
      /*
      if (abs(csi->type_parton[i_a][i_p]) < 7){
	int flag = 0;
	for (int i_j = 0; i_j < jet_algorithm_list.size(); i_j++){
	  if (jet_algorithm_list[i_j] == csi->type_parton[i_a][i_p]){flag++; break;}
	}
	if (flag){
	  ps_n_partonlevel[i_a][esi.observed_object["jet"]]++;
	}
      }
      */
      if (csi->type_parton[i_a][i_p] == 24){ps_n_partonlevel[i_a][esi.observed_object["wm"]]++;}
      if (csi->type_parton[i_a][i_p] == -24){ps_n_partonlevel[i_a][esi.observed_object["wp"]]++;}
      if (csi->type_parton[i_a][i_p] == 23){ps_n_partonlevel[i_a][esi.observed_object["z"]]++;}
      if (csi->type_parton[i_a][i_p] == 25){ps_n_partonlevel[i_a][esi.observed_object["h"]]++;}
    }

    // all potential jets are collected in ps_n_partonlevel[i_a]:
    /*
    logger << LOG_INFO << "ps_runtime_jet_algorithm.size() = " << ps_runtime_jet_algorithm.size() << endl;
    logger << LOG_INFO << "ps_runtime_jet_algorithm[" << i_a << "].size() = " << ps_runtime_jet_algorithm[i_a].size() << endl;
    logger << LOG_INFO << "esi.observed_object[jet]" << esi.observed_object["jet"] << endl;
    logger << LOG_INFO << "ps_n_partonlevel[" << i_a << "].size() = " << ps_n_partonlevel[i_a].size() << endl;
    */

    ps_n_partonlevel[i_a][esi.observed_object["jet"]] = ps_runtime_jet_algorithm[i_a].size();
    if (photon_jet_algorithm){ps_n_partonlevel[i_a][esi.observed_object["jet"]] += ps_runtime_photon[i_a].size();}
  }

  ps_runtime_order_inverse.resize(n_ps, vector<vector<int> > (esi.object_list.size()));
  for (int i_a = 0; i_a < n_ps; i_a++){
    logger << LOG_DEBUG << "csi->type_parton[" << i_a << "].size() = " << csi->type_parton[i_a].size() << endl;
    for (int i_p = 3; i_p < csi->type_parton[i_a].size(); i_p++){

      if (csi->type_parton[i_a][i_p] == 0){
	ps_runtime_order_inverse[i_a][esi.observed_object["jet"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["ljet"]].push_back(i_p);
      }
      if (abs(csi->type_parton[i_a][i_p]) > 0 && abs(csi->type_parton[i_a][i_p]) < 5){
	ps_runtime_order_inverse[i_a][esi.observed_object["jet"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["ljet"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == 5){
	ps_runtime_order_inverse[i_a][esi.observed_object["jet"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["ljet"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["bjet"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["bjet_b"]].push_back(i_p);
	// ljet...
      }
      if (csi->type_parton[i_a][i_p] == -5){
	ps_runtime_order_inverse[i_a][esi.observed_object["jet"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["ljet"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["bjet"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["bjet_bx"]].push_back(i_p);
	// ljet...
      }
      if (csi->type_parton[i_a][i_p] == 6){
	ps_runtime_order_inverse[i_a][esi.observed_object["tjet"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["top"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == -6){
	ps_runtime_order_inverse[i_a][esi.observed_object["tjet"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["atop"]].push_back(i_p);
      }

      if (csi->type_parton[i_a][i_p] == 11){
	ps_runtime_order_inverse[i_a][esi.observed_object["lep"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["lm"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["e"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["em"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == -11){
	ps_runtime_order_inverse[i_a][esi.observed_object["lep"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["lp"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["e"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["ep"]].push_back(i_p);
      }

      if (csi->type_parton[i_a][i_p] == 13){
	ps_runtime_order_inverse[i_a][esi.observed_object["lep"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["lm"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["mu"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["mum"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == -13){
	ps_runtime_order_inverse[i_a][esi.observed_object["lep"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["lp"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["mu"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["mup"]].push_back(i_p);
      }

      if (csi->type_parton[i_a][i_p] == 15){
	ps_runtime_order_inverse[i_a][esi.observed_object["lep"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["lm"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["tau"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["taum"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == -15){
	ps_runtime_order_inverse[i_a][esi.observed_object["lep"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["lp"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["tau"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["taup"]].push_back(i_p);
      }

      if (csi->type_parton[i_a][i_p] == 12){
	ps_runtime_order_inverse[i_a][esi.observed_object["nua"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nu"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nea"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["ne"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["missing"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == -12){
	ps_runtime_order_inverse[i_a][esi.observed_object["nua"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nux"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nea"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nex"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["missing"]].push_back(i_p);
      }

      if (csi->type_parton[i_a][i_p] == 14){
	ps_runtime_order_inverse[i_a][esi.observed_object["nua"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nu"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nma"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nm"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["missing"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == -14){
	ps_runtime_order_inverse[i_a][esi.observed_object["nua"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nux"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nma"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nmx"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["missing"]].push_back(i_p);
      }

      if (csi->type_parton[i_a][i_p] == 16){
	ps_runtime_order_inverse[i_a][esi.observed_object["nua"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nu"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nta"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nt"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["missing"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == -16){
	ps_runtime_order_inverse[i_a][esi.observed_object["nua"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nux"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["nta"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["ntx"]].push_back(i_p);
	ps_runtime_order_inverse[i_a][esi.observed_object["missing"]].push_back(i_p);
      }

      if (csi->type_parton[i_a][i_p] == 22){
	ps_runtime_order_inverse[i_a][esi.observed_object["photon"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == 23){
	ps_runtime_order_inverse[i_a][esi.observed_object["z"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == 24){
	ps_runtime_order_inverse[i_a][esi.observed_object["wm"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == -24){
	ps_runtime_order_inverse[i_a][esi.observed_object["wp"]].push_back(i_p);
      }
      if (csi->type_parton[i_a][i_p] == 25){
	ps_runtime_order_inverse[i_a][esi.observed_object["h"]].push_back(i_p);
      }
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void observable_set::determine_object_definition(){
  static  Logger logger("observable_set::determine_object_definition");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG << endl << "Object definition after definition (fully inclusive): " << endl;
  logger << LOG_DEBUG << setw(5) << "no" << setw(20) << "object" << setw(20) << "n_observed_min" << setw(20) << "n_observed_max" << setw(20) << "define_pT" << setw(20) << "define_ET" << setw(20) << "define_eta"  << setw(20) << "define_y"<< endl;
  for (int i = 1; i < esi.object_list.size(); i++){
    //    logger << LOG_DEBUG << setw(5) << i << setw(20) << esi.object_list[i] << setw(20) << n_observed_min[i] << setw(20) << n_observed_max[i] << setw(20) << define_pT[i] << setw(20) << define_ET[i] << setw(20) << define_eta[i] << setw(20) << define_y[i] << endl;
    logger << LOG_DEBUG << setw(5) << i << setw(20) << esi.object_list[i] << setw(20) << esi.pda[i].n_observed_min << setw(20) << esi.pda[i].n_observed_max << setw(20) << esi.pda[i].define_pT << setw(20) << esi.pda[i].define_ET << setw(20) << esi.pda[i].define_eta << setw(20) << esi.pda[i].define_y << endl;
  }
  logger.newLine(LOG_DEBUG);

  // adapt define_pT (tjet-top/atop) !

  if (esi.pda[esi.observed_object["top"]].define_pT == 0. && esi.pda[esi.observed_object["tjet"]].define_pT > 0.){esi.pda[esi.observed_object["top"]].define_pT = esi.pda[esi.observed_object["tjet"]].define_pT;}
  if (esi.pda[esi.observed_object["atop"]].define_pT == 0. && esi.pda[esi.observed_object["tjet"]].define_pT > 0.){esi.pda[esi.observed_object["atop"]].define_pT = esi.pda[esi.observed_object["tjet"]].define_pT;}
  if (esi.pda[esi.observed_object["tjet"]].define_pT < esi.pda[esi.observed_object["top"]].define_pT &&
      esi.pda[esi.observed_object["tjet"]].define_pT < esi.pda[esi.observed_object["atop"]].define_pT){
    if (esi.pda[esi.observed_object["top"]].define_pT >= esi.pda[esi.observed_object["atop"]].define_pT){esi.pda[esi.observed_object["tjet"]].define_pT = esi.pda[esi.observed_object["atop"]].define_pT;}
    else {esi.pda[esi.observed_object["tjet"]].define_pT = esi.pda[esi.observed_object["top"]].define_pT;}
  }
  
  // adapt define_pT (bjet-bjet_b/bjet_bx) !
  
  if (esi.pda[esi.observed_object["bjet_b"]].define_pT == 0. && esi.pda[esi.observed_object["bjet"]].define_pT > 0.){esi.pda[esi.observed_object["bjet_b"]].define_pT = esi.pda[esi.observed_object["bjet"]].define_pT;}
  if (esi.pda[esi.observed_object["bjet_bx"]].define_pT == 0. && esi.pda[esi.observed_object["bjet"]].define_pT > 0.){esi.pda[esi.observed_object["bjet_bx"]].define_pT = esi.pda[esi.observed_object["bjet"]].define_pT;}
  if (esi.pda[esi.observed_object["bjet"]].define_pT < esi.pda[esi.observed_object["bjet_b"]].define_pT &&
      esi.pda[esi.observed_object["bjet"]].define_pT < esi.pda[esi.observed_object["bjet_bx"]].define_pT){
    if (esi.pda[esi.observed_object["bjet_b"]].define_pT >= esi.pda[esi.observed_object["bjet_bx"]].define_pT){esi.pda[esi.observed_object["bjet"]].define_pT = esi.pda[esi.observed_object["bjet_bx"]].define_pT;}
    else {esi.pda[esi.observed_object["bjet"]].define_pT = esi.pda[esi.observed_object["bjet_b"]].define_pT;}
  }
  /*
  // adapt define_pT (jet-ljet/bjet) !
  if (esi.pda[esi.observed_object["ljet"]].define_pT == 0. && esi.pda[esi.observed_object["jet"]].define_pT > 0.){esi.pda[esi.observed_object["ljet"]].define_pT = esi.pda[esi.observed_object["jet"]].define_pT;}
  if (esi.pda[esi.observed_object["bjet"]].define_pT == 0. && esi.pda[esi.observed_object["jet"]].define_pT > 0.){esi.pda[esi.observed_object["bjet"]].define_pT = esi.pda[esi.observed_object["jet"]].define_pT;}

  if (esi.pda[esi.observed_object["jet"]].define_pT == 0. && esi.pda[esi.observed_object["bjet"]].define_pT == 0. && esi.pda[esi.observed_object["ljet"]].define_pT > 0.){esi.pda[esi.observed_object["jet"]].define_pT = esi.pda[esi.observed_object["ljet"]].define_pT;}
  if (esi.pda[esi.observed_object["jet"]].define_pT == 0. && esi.pda[esi.observed_object["ljet"]].define_pT == 0. && esi.pda[esi.observed_object["bjet"]].define_pT > 0.){esi.pda[esi.observed_object["jet"]].define_pT = esi.pda[esi.observed_object["bjet"]].define_pT;}

  if (esi.pda[esi.observed_object["jet"]].define_pT < esi.pda[esi.observed_object["ljet"]].define_pT &&
      esi.pda[esi.observed_object["jet"]].define_pT < esi.pda[esi.observed_object["bjet"]].define_pT){
    if (esi.pda[esi.observed_object["ljet"]].define_pT >= esi.pda[esi.observed_object["bjet"]].define_pT){esi.pda[esi.observed_object["jet"]].define_pT = esi.pda[esi.observed_object["bjet"]].define_pT;}
    else {esi.pda[esi.observed_object["jet"]].define_pT = esi.pda[esi.observed_object["ljet"]].define_pT;}
  }
  */
  // adapt define_pT (lep-lm/lp) !
  if (esi.pda[esi.observed_object["lm"]].define_pT == 0. && esi.pda[esi.observed_object["lep"]].define_pT > 0.){esi.pda[esi.observed_object["lm"]].define_pT = esi.pda[esi.observed_object["lep"]].define_pT;}
  if (esi.pda[esi.observed_object["lp"]].define_pT == 0. && esi.pda[esi.observed_object["lep"]].define_pT > 0.){esi.pda[esi.observed_object["lp"]].define_pT = esi.pda[esi.observed_object["lep"]].define_pT;}
  if (esi.pda[esi.observed_object["lep"]].define_pT < esi.pda[esi.observed_object["lm"]].define_pT &&
      esi.pda[esi.observed_object["lep"]].define_pT < esi.pda[esi.observed_object["lp"]].define_pT){
    if (esi.pda[esi.observed_object["lm"]].define_pT >= esi.pda[esi.observed_object["lp"]].define_pT){esi.pda[esi.observed_object["lep"]].define_pT = esi.pda[esi.observed_object["lp"]].define_pT;}
    else {esi.pda[esi.observed_object["lep"]].define_pT = esi.pda[esi.observed_object["lm"]].define_pT;}
  }

  // adapt define_pT (lep-e/mu/tau) !
  if (esi.pda[esi.observed_object["e"]].define_pT == 0. && esi.pda[esi.observed_object["lep"]].define_pT > 0.){esi.pda[esi.observed_object["e"]].define_pT = esi.pda[esi.observed_object["lep"]].define_pT;}
  if (esi.pda[esi.observed_object["mu"]].define_pT == 0. && esi.pda[esi.observed_object["lep"]].define_pT > 0.){esi.pda[esi.observed_object["mu"]].define_pT = esi.pda[esi.observed_object["lep"]].define_pT;}
  if (esi.pda[esi.observed_object["tau"]].define_pT == 0. && esi.pda[esi.observed_object["lep"]].define_pT > 0.){esi.pda[esi.observed_object["tau"]].define_pT = esi.pda[esi.observed_object["lep"]].define_pT;}
  if (esi.pda[esi.observed_object["lep"]].define_pT < esi.pda[esi.observed_object["e"]].define_pT &&
      esi.pda[esi.observed_object["lep"]].define_pT < esi.pda[esi.observed_object["mu"]].define_pT &&
      esi.pda[esi.observed_object["lep"]].define_pT < esi.pda[esi.observed_object["tau"]].define_pT){
    if (esi.pda[esi.observed_object["tau"]].define_pT >= esi.pda[esi.observed_object["e"]].define_pT &&
	esi.pda[esi.observed_object["mu"]].define_pT >= esi.pda[esi.observed_object["e"]].define_pT){esi.pda[esi.observed_object["lep"]].define_pT = esi.pda[esi.observed_object["e"]].define_pT;}
    else if (esi.pda[esi.observed_object["tau"]].define_pT >= esi.pda[esi.observed_object["mu"]].define_pT){esi.pda[esi.observed_object["lep"]].define_pT = esi.pda[esi.observed_object["mu"]].define_pT;}
    else {esi.pda[esi.observed_object["lep"]].define_pT = esi.pda[esi.observed_object["tau"]].define_pT;}
  }

  // adapt define_pT (e-em/ep) !
  if (esi.pda[esi.observed_object["em"]].define_pT == 0. && esi.pda[esi.observed_object["e"]].define_pT > 0.){esi.pda[esi.observed_object["em"]].define_pT = esi.pda[esi.observed_object["e"]].define_pT;}
  if (esi.pda[esi.observed_object["ep"]].define_pT == 0. && esi.pda[esi.observed_object["e"]].define_pT > 0.){esi.pda[esi.observed_object["ep"]].define_pT = esi.pda[esi.observed_object["e"]].define_pT;}
  if (esi.pda[esi.observed_object["e"]].define_pT < esi.pda[esi.observed_object["em"]].define_pT &&
      esi.pda[esi.observed_object["e"]].define_pT < esi.pda[esi.observed_object["ep"]].define_pT){
    if (esi.pda[esi.observed_object["em"]].define_pT > esi.pda[esi.observed_object["ep"]].define_pT){esi.pda[esi.observed_object["e"]].define_pT = esi.pda[esi.observed_object["ep"]].define_pT;}
    else {esi.pda[esi.observed_object["e"]].define_pT = esi.pda[esi.observed_object["em"]].define_pT;}
  }

  // adapt define_pT (mu-mum/mup) !
  if (esi.pda[esi.observed_object["mum"]].define_pT == 0. && esi.pda[esi.observed_object["mu"]].define_pT > 0.){esi.pda[esi.observed_object["mum"]].define_pT = esi.pda[esi.observed_object["mu"]].define_pT;}
  if (esi.pda[esi.observed_object["mup"]].define_pT == 0. && esi.pda[esi.observed_object["mu"]].define_pT > 0.){esi.pda[esi.observed_object["mup"]].define_pT = esi.pda[esi.observed_object["mu"]].define_pT;}
  if (esi.pda[esi.observed_object["mu"]].define_pT < esi.pda[esi.observed_object["mum"]].define_pT &&
      esi.pda[esi.observed_object["mu"]].define_pT < esi.pda[esi.observed_object["mup"]].define_pT){
    if (esi.pda[esi.observed_object["mum"]].define_pT > esi.pda[esi.observed_object["mup"]].define_pT){esi.pda[esi.observed_object["mu"]].define_pT = esi.pda[esi.observed_object["mup"]].define_pT;}
    else {esi.pda[esi.observed_object["mu"]].define_pT = esi.pda[esi.observed_object["mum"]].define_pT;}
  }

  // adapt define_pT (tau-taum/taup) !
  if (esi.pda[esi.observed_object["taum"]].define_pT == 0. && esi.pda[esi.observed_object["tau"]].define_pT > 0.){esi.pda[esi.observed_object["taum"]].define_pT = esi.pda[esi.observed_object["tau"]].define_pT;}
  if (esi.pda[esi.observed_object["taup"]].define_pT == 0. && esi.pda[esi.observed_object["tau"]].define_pT > 0.){esi.pda[esi.observed_object["taup"]].define_pT = esi.pda[esi.observed_object["tau"]].define_pT;}
  if (esi.pda[esi.observed_object["tau"]].define_pT < esi.pda[esi.observed_object["taum"]].define_pT &&
      esi.pda[esi.observed_object["tau"]].define_pT < esi.pda[esi.observed_object["taup"]].define_pT){
    if (esi.pda[esi.observed_object["taum"]].define_pT > esi.pda[esi.observed_object["taup"]].define_pT){esi.pda[esi.observed_object["tau"]].define_pT = esi.pda[esi.observed_object["taup"]].define_pT;}
    else {esi.pda[esi.observed_object["tau"]].define_pT = esi.pda[esi.observed_object["taum"]].define_pT;}
  }

  // adapt define_pT (nua-nu/nux) !
  if (esi.pda[esi.observed_object["nu"]].define_pT == 0. && esi.pda[esi.observed_object["nua"]].define_pT > 0.){esi.pda[esi.observed_object["nu"]].define_pT = esi.pda[esi.observed_object["nua"]].define_pT;}
  if (esi.pda[esi.observed_object["nux"]].define_pT == 0. && esi.pda[esi.observed_object["nua"]].define_pT > 0.){esi.pda[esi.observed_object["nux"]].define_pT = esi.pda[esi.observed_object["nua"]].define_pT;}
  if (esi.pda[esi.observed_object["nua"]].define_pT < esi.pda[esi.observed_object["nu"]].define_pT &&
      esi.pda[esi.observed_object["nua"]].define_pT < esi.pda[esi.observed_object["nux"]].define_pT){
    if (esi.pda[esi.observed_object["nu"]].define_pT >= esi.pda[esi.observed_object["nux"]].define_pT){esi.pda[esi.observed_object["nua"]].define_pT = esi.pda[esi.observed_object["nux"]].define_pT;}
    else {esi.pda[esi.observed_object["nua"]].define_pT = esi.pda[esi.observed_object["nu"]].define_pT;}
  }

  // adapt define_pT (nua-nea/nma/nta) !
  if (esi.pda[esi.observed_object["nea"]].define_pT == 0. && esi.pda[esi.observed_object["nua"]].define_pT > 0.){esi.pda[esi.observed_object["nea"]].define_pT = esi.pda[esi.observed_object["nua"]].define_pT;}
  if (esi.pda[esi.observed_object["nma"]].define_pT == 0. && esi.pda[esi.observed_object["nua"]].define_pT > 0.){esi.pda[esi.observed_object["nma"]].define_pT = esi.pda[esi.observed_object["nua"]].define_pT;}
  if (esi.pda[esi.observed_object["nta"]].define_pT == 0. && esi.pda[esi.observed_object["nua"]].define_pT > 0.){esi.pda[esi.observed_object["nta"]].define_pT = esi.pda[esi.observed_object["nua"]].define_pT;}
  if (esi.pda[esi.observed_object["nua"]].define_pT < esi.pda[esi.observed_object["nea"]].define_pT &&
      esi.pda[esi.observed_object["nua"]].define_pT < esi.pda[esi.observed_object["nma"]].define_pT &&
      esi.pda[esi.observed_object["nua"]].define_pT < esi.pda[esi.observed_object["nta"]].define_pT){
    if (esi.pda[esi.observed_object["nta"]].define_pT >= esi.pda[esi.observed_object["nea"]].define_pT &&
	esi.pda[esi.observed_object["nma"]].define_pT >= esi.pda[esi.observed_object["nea"]].define_pT){esi.pda[esi.observed_object["nua"]].define_pT = esi.pda[esi.observed_object["nea"]].define_pT;}
    else if (esi.pda[esi.observed_object["nta"]].define_pT >= esi.pda[esi.observed_object["nma"]].define_pT){esi.pda[esi.observed_object["nua"]].define_pT = esi.pda[esi.observed_object["nma"]].define_pT;}
    else {esi.pda[esi.observed_object["nua"]].define_pT = esi.pda[esi.observed_object["nta"]].define_pT;}
  }

  // adapt define_pT (nea-ne/nex) !
  if (esi.pda[esi.observed_object["ne"]].define_pT == 0. && esi.pda[esi.observed_object["nea"]].define_pT > 0.){esi.pda[esi.observed_object["ne"]].define_pT = esi.pda[esi.observed_object["nea"]].define_pT;}
  if (esi.pda[esi.observed_object["nex"]].define_pT == 0. && esi.pda[esi.observed_object["nea"]].define_pT > 0.){esi.pda[esi.observed_object["nex"]].define_pT = esi.pda[esi.observed_object["nea"]].define_pT;}
  if (esi.pda[esi.observed_object["nea"]].define_pT < esi.pda[esi.observed_object["ne"]].define_pT &&
      esi.pda[esi.observed_object["nea"]].define_pT < esi.pda[esi.observed_object["nex"]].define_pT){
    if (esi.pda[esi.observed_object["ne"]].define_pT > esi.pda[esi.observed_object["nex"]].define_pT){esi.pda[esi.observed_object["nea"]].define_pT = esi.pda[esi.observed_object["nex"]].define_pT;}
    else {esi.pda[esi.observed_object["nea"]].define_pT = esi.pda[esi.observed_object["ne"]].define_pT;}
  }

  // adapt define_pT (nma-nm/nmx) !
  if (esi.pda[esi.observed_object["nm"]].define_pT == 0. && esi.pda[esi.observed_object["nma"]].define_pT > 0.){esi.pda[esi.observed_object["nm"]].define_pT = esi.pda[esi.observed_object["nma"]].define_pT;}
  if (esi.pda[esi.observed_object["nmx"]].define_pT == 0. && esi.pda[esi.observed_object["nma"]].define_pT > 0.){esi.pda[esi.observed_object["nmx"]].define_pT = esi.pda[esi.observed_object["nma"]].define_pT;}
  if (esi.pda[esi.observed_object["nma"]].define_pT < esi.pda[esi.observed_object["nm"]].define_pT &&
      esi.pda[esi.observed_object["nma"]].define_pT < esi.pda[esi.observed_object["nmx"]].define_pT){
    if (esi.pda[esi.observed_object["nm"]].define_pT > esi.pda[esi.observed_object["nmx"]].define_pT){esi.pda[esi.observed_object["nma"]].define_pT = esi.pda[esi.observed_object["nmx"]].define_pT;}
    else {esi.pda[esi.observed_object["nma"]].define_pT = esi.pda[esi.observed_object["nm"]].define_pT;}
  }

  // adapt define_pT (nta-nt/ntx) !
  if (esi.pda[esi.observed_object["nt"]].define_pT == 0. && esi.pda[esi.observed_object["nta"]].define_pT > 0.){esi.pda[esi.observed_object["nt"]].define_pT = esi.pda[esi.observed_object["nta"]].define_pT;}
  if (esi.pda[esi.observed_object["ntx"]].define_pT == 0. && esi.pda[esi.observed_object["nta"]].define_pT > 0.){esi.pda[esi.observed_object["ntx"]].define_pT = esi.pda[esi.observed_object["nta"]].define_pT;}
  if (esi.pda[esi.observed_object["nta"]].define_pT < esi.pda[esi.observed_object["nt"]].define_pT &&
      esi.pda[esi.observed_object["nta"]].define_pT < esi.pda[esi.observed_object["ntx"]].define_pT){
    if (esi.pda[esi.observed_object["nt"]].define_pT > esi.pda[esi.observed_object["ntx"]].define_pT){esi.pda[esi.observed_object["nta"]].define_pT = esi.pda[esi.observed_object["ntx"]].define_pT;}
    else {esi.pda[esi.observed_object["nta"]].define_pT = esi.pda[esi.observed_object["nt"]].define_pT;}
  }



  // adapt define_ET (tjet-top/atop) !
  if (esi.pda[esi.observed_object["top"]].define_ET == 0. && esi.pda[esi.observed_object["tjet"]].define_ET > 0.){esi.pda[esi.observed_object["top"]].define_ET = esi.pda[esi.observed_object["tjet"]].define_ET;}
  if (esi.pda[esi.observed_object["atop"]].define_ET == 0. && esi.pda[esi.observed_object["tjet"]].define_ET > 0.){esi.pda[esi.observed_object["atop"]].define_ET = esi.pda[esi.observed_object["tjet"]].define_ET;}
  if (esi.pda[esi.observed_object["tjet"]].define_ET < esi.pda[esi.observed_object["top"]].define_ET &&
      esi.pda[esi.observed_object["tjet"]].define_ET < esi.pda[esi.observed_object["atop"]].define_ET){
    if (esi.pda[esi.observed_object["top"]].define_ET >= esi.pda[esi.observed_object["atop"]].define_ET){esi.pda[esi.observed_object["tjet"]].define_ET = esi.pda[esi.observed_object["atop"]].define_ET;}
    else {esi.pda[esi.observed_object["tjet"]].define_ET = esi.pda[esi.observed_object["top"]].define_ET;}
  }
  
  // adapt define_ET (bjet-bjet_b/bjet_bx) !
  if (esi.pda[esi.observed_object["bjet_b"]].define_ET == 0. && esi.pda[esi.observed_object["bjet"]].define_ET > 0.){esi.pda[esi.observed_object["bjet_b"]].define_ET = esi.pda[esi.observed_object["bjet"]].define_ET;}
  if (esi.pda[esi.observed_object["bjet_bx"]].define_ET == 0. && esi.pda[esi.observed_object["bjet"]].define_ET > 0.){esi.pda[esi.observed_object["bjet_bx"]].define_ET = esi.pda[esi.observed_object["bjet"]].define_ET;}
  if (esi.pda[esi.observed_object["bjet"]].define_ET < esi.pda[esi.observed_object["bjet_b"]].define_ET &&
      esi.pda[esi.observed_object["bjet"]].define_ET < esi.pda[esi.observed_object["bjet_bx"]].define_ET){
    if (esi.pda[esi.observed_object["bjet_b"]].define_ET >= esi.pda[esi.observed_object["bjet_bx"]].define_ET){esi.pda[esi.observed_object["bjet"]].define_ET = esi.pda[esi.observed_object["bjet_bx"]].define_ET;}
    else {esi.pda[esi.observed_object["bjet"]].define_ET = esi.pda[esi.observed_object["bjet_b"]].define_ET;}
  }
  /*
  // adapt define_ET (jet-ljet/bjet) !
  if (esi.pda[esi.observed_object["ljet"]].define_ET == 0. && esi.pda[esi.observed_object["jet"]].define_ET > 0.){esi.pda[esi.observed_object["ljet"]].define_ET = esi.pda[esi.observed_object["jet"]].define_ET;}
  if (esi.pda[esi.observed_object["bjet"]].define_ET == 0. && esi.pda[esi.observed_object["jet"]].define_ET > 0.){esi.pda[esi.observed_object["bjet"]].define_ET = esi.pda[esi.observed_object["jet"]].define_ET;}
  if (esi.pda[esi.observed_object["jet"]].define_ET == 0. && esi.pda[esi.observed_object["bjet"]].define_ET == 0. && esi.pda[esi.observed_object["ljet"]].define_ET > 0.){esi.pda[esi.observed_object["jet"]].define_ET = esi.pda[esi.observed_object["ljet"]].define_ET;}
  if (esi.pda[esi.observed_object["jet"]].define_ET == 0. && esi.pda[esi.observed_object["ljet"]].define_ET == 0. && esi.pda[esi.observed_object["bjet"]].define_ET > 0.){esi.pda[esi.observed_object["jet"]].define_ET = esi.pda[esi.observed_object["bjet"]].define_ET;}

  if (esi.pda[esi.observed_object["jet"]].define_ET < esi.pda[esi.observed_object["ljet"]].define_ET &&
      esi.pda[esi.observed_object["jet"]].define_ET < esi.pda[esi.observed_object["bjet"]].define_ET){
    if (esi.pda[esi.observed_object["ljet"]].define_ET >= esi.pda[esi.observed_object["bjet"]].define_ET){esi.pda[esi.observed_object["jet"]].define_ET = esi.pda[esi.observed_object["bjet"]].define_ET;}
    else {esi.pda[esi.observed_object["jet"]].define_ET = esi.pda[esi.observed_object["ljet"]].define_ET;}
  }
  */
  // adapt define_ET (lep-lm/lp) !
  if (esi.pda[esi.observed_object["lm"]].define_ET == 0. && esi.pda[esi.observed_object["lep"]].define_ET > 0.){esi.pda[esi.observed_object["lm"]].define_ET = esi.pda[esi.observed_object["lep"]].define_ET;}
  if (esi.pda[esi.observed_object["lp"]].define_ET == 0. && esi.pda[esi.observed_object["lep"]].define_ET > 0.){esi.pda[esi.observed_object["lp"]].define_ET = esi.pda[esi.observed_object["lep"]].define_ET;}
  if (esi.pda[esi.observed_object["lep"]].define_ET < esi.pda[esi.observed_object["lm"]].define_ET &&
      esi.pda[esi.observed_object["lep"]].define_ET < esi.pda[esi.observed_object["lp"]].define_ET){
    if (esi.pda[esi.observed_object["lm"]].define_ET >= esi.pda[esi.observed_object["lp"]].define_ET){esi.pda[esi.observed_object["lep"]].define_ET = esi.pda[esi.observed_object["lp"]].define_ET;}
    else {esi.pda[esi.observed_object["lep"]].define_ET = esi.pda[esi.observed_object["lm"]].define_ET;}
  }

  // adapt define_ET (lep-e/mu/tau) !
  if (esi.pda[esi.observed_object["e"]].define_ET == 0. && esi.pda[esi.observed_object["lep"]].define_ET > 0.){esi.pda[esi.observed_object["e"]].define_ET = esi.pda[esi.observed_object["lep"]].define_ET;}
  if (esi.pda[esi.observed_object["mu"]].define_ET == 0. && esi.pda[esi.observed_object["lep"]].define_ET > 0.){esi.pda[esi.observed_object["mu"]].define_ET = esi.pda[esi.observed_object["lep"]].define_ET;}
  if (esi.pda[esi.observed_object["tau"]].define_ET == 0. && esi.pda[esi.observed_object["lep"]].define_ET > 0.){esi.pda[esi.observed_object["tau"]].define_ET = esi.pda[esi.observed_object["lep"]].define_ET;}
  if (esi.pda[esi.observed_object["lep"]].define_ET < esi.pda[esi.observed_object["e"]].define_ET &&
      esi.pda[esi.observed_object["lep"]].define_ET < esi.pda[esi.observed_object["mu"]].define_ET &&
      esi.pda[esi.observed_object["lep"]].define_ET < esi.pda[esi.observed_object["tau"]].define_ET){
    if (esi.pda[esi.observed_object["tau"]].define_ET >= esi.pda[esi.observed_object["e"]].define_ET &&
	esi.pda[esi.observed_object["mu"]].define_ET >= esi.pda[esi.observed_object["e"]].define_ET){esi.pda[esi.observed_object["lep"]].define_ET = esi.pda[esi.observed_object["e"]].define_ET;}
    else if (esi.pda[esi.observed_object["tau"]].define_ET >= esi.pda[esi.observed_object["mu"]].define_ET){esi.pda[esi.observed_object["lep"]].define_ET = esi.pda[esi.observed_object["mu"]].define_ET;}
    else {esi.pda[esi.observed_object["lep"]].define_ET = esi.pda[esi.observed_object["tau"]].define_ET;}
  }

  // adapt define_ET (e-em/ep) !
  if (esi.pda[esi.observed_object["em"]].define_ET == 0. && esi.pda[esi.observed_object["e"]].define_ET > 0.){esi.pda[esi.observed_object["em"]].define_ET = esi.pda[esi.observed_object["e"]].define_ET;}
  if (esi.pda[esi.observed_object["ep"]].define_ET == 0. && esi.pda[esi.observed_object["e"]].define_ET > 0.){esi.pda[esi.observed_object["ep"]].define_ET = esi.pda[esi.observed_object["e"]].define_ET;}
  if (esi.pda[esi.observed_object["e"]].define_ET < esi.pda[esi.observed_object["em"]].define_ET &&
      esi.pda[esi.observed_object["e"]].define_ET < esi.pda[esi.observed_object["ep"]].define_ET){
    if (esi.pda[esi.observed_object["em"]].define_ET > esi.pda[esi.observed_object["ep"]].define_ET){esi.pda[esi.observed_object["e"]].define_ET = esi.pda[esi.observed_object["ep"]].define_ET;}
    else {esi.pda[esi.observed_object["e"]].define_ET = esi.pda[esi.observed_object["em"]].define_ET;}
  }

  // adapt define_ET (mu-mum/mup) !
  if (esi.pda[esi.observed_object["mum"]].define_ET == 0. && esi.pda[esi.observed_object["mu"]].define_ET > 0.){esi.pda[esi.observed_object["mum"]].define_ET = esi.pda[esi.observed_object["mu"]].define_ET;}
  if (esi.pda[esi.observed_object["mup"]].define_ET == 0. && esi.pda[esi.observed_object["mu"]].define_ET > 0.){esi.pda[esi.observed_object["mup"]].define_ET = esi.pda[esi.observed_object["mu"]].define_ET;}
  if (esi.pda[esi.observed_object["mu"]].define_ET < esi.pda[esi.observed_object["mum"]].define_ET &&
      esi.pda[esi.observed_object["mu"]].define_ET < esi.pda[esi.observed_object["mup"]].define_ET){
    if (esi.pda[esi.observed_object["mum"]].define_ET > esi.pda[esi.observed_object["mup"]].define_ET){esi.pda[esi.observed_object["mu"]].define_ET = esi.pda[esi.observed_object["mup"]].define_ET;}
    else {esi.pda[esi.observed_object["mu"]].define_ET = esi.pda[esi.observed_object["mum"]].define_ET;}
  }

  // adapt define_ET (tau-taum/taup) !
  if (esi.pda[esi.observed_object["taum"]].define_ET == 0. && esi.pda[esi.observed_object["tau"]].define_ET > 0.){esi.pda[esi.observed_object["taum"]].define_ET = esi.pda[esi.observed_object["tau"]].define_ET;}
  if (esi.pda[esi.observed_object["taup"]].define_ET == 0. && esi.pda[esi.observed_object["tau"]].define_ET > 0.){esi.pda[esi.observed_object["taup"]].define_ET = esi.pda[esi.observed_object["tau"]].define_ET;}
  if (esi.pda[esi.observed_object["tau"]].define_ET < esi.pda[esi.observed_object["taum"]].define_ET &&
      esi.pda[esi.observed_object["tau"]].define_ET < esi.pda[esi.observed_object["taup"]].define_ET){
    if (esi.pda[esi.observed_object["taum"]].define_ET > esi.pda[esi.observed_object["taup"]].define_ET){esi.pda[esi.observed_object["tau"]].define_ET = esi.pda[esi.observed_object["taup"]].define_ET;}
    else {esi.pda[esi.observed_object["tau"]].define_ET = esi.pda[esi.observed_object["taum"]].define_ET;}
  }

  // adapt define_ET (nua-nu/nux) !
  if (esi.pda[esi.observed_object["nu"]].define_ET == 0. && esi.pda[esi.observed_object["nua"]].define_ET > 0.){esi.pda[esi.observed_object["nu"]].define_ET = esi.pda[esi.observed_object["nua"]].define_ET;}
  if (esi.pda[esi.observed_object["nux"]].define_ET == 0. && esi.pda[esi.observed_object["nua"]].define_ET > 0.){esi.pda[esi.observed_object["nux"]].define_ET = esi.pda[esi.observed_object["nua"]].define_ET;}
  if (esi.pda[esi.observed_object["nua"]].define_ET < esi.pda[esi.observed_object["nu"]].define_ET &&
      esi.pda[esi.observed_object["nua"]].define_ET < esi.pda[esi.observed_object["nux"]].define_ET){
    if (esi.pda[esi.observed_object["nu"]].define_ET >= esi.pda[esi.observed_object["nux"]].define_ET){esi.pda[esi.observed_object["nua"]].define_ET = esi.pda[esi.observed_object["nux"]].define_ET;}
    else {esi.pda[esi.observed_object["nua"]].define_ET = esi.pda[esi.observed_object["nu"]].define_ET;}
  }

  // adapt define_ET (nua-nea/nma/nta) !
  if (esi.pda[esi.observed_object["nea"]].define_ET == 0. && esi.pda[esi.observed_object["nua"]].define_ET > 0.){esi.pda[esi.observed_object["nea"]].define_ET = esi.pda[esi.observed_object["nua"]].define_ET;}
  if (esi.pda[esi.observed_object["nma"]].define_ET == 0. && esi.pda[esi.observed_object["nua"]].define_ET > 0.){esi.pda[esi.observed_object["nma"]].define_ET = esi.pda[esi.observed_object["nua"]].define_ET;}
  if (esi.pda[esi.observed_object["nta"]].define_ET == 0. && esi.pda[esi.observed_object["nua"]].define_ET > 0.){esi.pda[esi.observed_object["nta"]].define_ET = esi.pda[esi.observed_object["nua"]].define_ET;}
  if (esi.pda[esi.observed_object["nua"]].define_ET < esi.pda[esi.observed_object["nea"]].define_ET &&
      esi.pda[esi.observed_object["nua"]].define_ET < esi.pda[esi.observed_object["nma"]].define_ET &&
      esi.pda[esi.observed_object["nua"]].define_ET < esi.pda[esi.observed_object["nta"]].define_ET){
    if (esi.pda[esi.observed_object["nta"]].define_ET >= esi.pda[esi.observed_object["nea"]].define_ET &&
	esi.pda[esi.observed_object["nma"]].define_ET >= esi.pda[esi.observed_object["nea"]].define_ET){esi.pda[esi.observed_object["nua"]].define_ET = esi.pda[esi.observed_object["nea"]].define_ET;}
    else if (esi.pda[esi.observed_object["nta"]].define_ET >= esi.pda[esi.observed_object["nma"]].define_ET){esi.pda[esi.observed_object["nua"]].define_ET = esi.pda[esi.observed_object["nma"]].define_ET;}
    else {esi.pda[esi.observed_object["nua"]].define_ET = esi.pda[esi.observed_object["nta"]].define_ET;}
  }

  // adapt define_ET (nea-ne/nex) !
  if (esi.pda[esi.observed_object["ne"]].define_ET == 0. && esi.pda[esi.observed_object["nea"]].define_ET > 0.){esi.pda[esi.observed_object["ne"]].define_ET = esi.pda[esi.observed_object["nea"]].define_ET;}
  if (esi.pda[esi.observed_object["nex"]].define_ET == 0. && esi.pda[esi.observed_object["nea"]].define_ET > 0.){esi.pda[esi.observed_object["nex"]].define_ET = esi.pda[esi.observed_object["nea"]].define_ET;}
  if (esi.pda[esi.observed_object["nea"]].define_ET < esi.pda[esi.observed_object["ne"]].define_ET &&
      esi.pda[esi.observed_object["nea"]].define_ET < esi.pda[esi.observed_object["nex"]].define_ET){
    if (esi.pda[esi.observed_object["ne"]].define_ET > esi.pda[esi.observed_object["nex"]].define_ET){esi.pda[esi.observed_object["nea"]].define_ET = esi.pda[esi.observed_object["nex"]].define_ET;}
    else {esi.pda[esi.observed_object["nea"]].define_ET = esi.pda[esi.observed_object["ne"]].define_ET;}
  }

  // adapt define_ET (nma-nm/nmx) !
  if (esi.pda[esi.observed_object["nm"]].define_ET == 0. && esi.pda[esi.observed_object["nma"]].define_ET > 0.){esi.pda[esi.observed_object["nm"]].define_ET = esi.pda[esi.observed_object["nma"]].define_ET;}
  if (esi.pda[esi.observed_object["nmx"]].define_ET == 0. && esi.pda[esi.observed_object["nma"]].define_ET > 0.){esi.pda[esi.observed_object["nmx"]].define_ET = esi.pda[esi.observed_object["nma"]].define_ET;}
  if (esi.pda[esi.observed_object["nma"]].define_ET < esi.pda[esi.observed_object["nm"]].define_ET &&
      esi.pda[esi.observed_object["nma"]].define_ET < esi.pda[esi.observed_object["nmx"]].define_ET){
    if (esi.pda[esi.observed_object["nm"]].define_ET > esi.pda[esi.observed_object["nmx"]].define_ET){esi.pda[esi.observed_object["nma"]].define_ET = esi.pda[esi.observed_object["nmx"]].define_ET;}
    else {esi.pda[esi.observed_object["nma"]].define_ET = esi.pda[esi.observed_object["nm"]].define_ET;}
  }

  // adapt define_ET (nta-nt/ntx) !
  if (esi.pda[esi.observed_object["nt"]].define_ET == 0. && esi.pda[esi.observed_object["nta"]].define_ET > 0.){esi.pda[esi.observed_object["nt"]].define_ET = esi.pda[esi.observed_object["nta"]].define_ET;}
  if (esi.pda[esi.observed_object["ntx"]].define_ET == 0. && esi.pda[esi.observed_object["nta"]].define_ET > 0.){esi.pda[esi.observed_object["ntx"]].define_ET = esi.pda[esi.observed_object["nta"]].define_ET;}
  if (esi.pda[esi.observed_object["nta"]].define_ET < esi.pda[esi.observed_object["nt"]].define_ET &&
      esi.pda[esi.observed_object["nta"]].define_ET < esi.pda[esi.observed_object["ntx"]].define_ET){
    if (esi.pda[esi.observed_object["nt"]].define_ET > esi.pda[esi.observed_object["ntx"]].define_ET){esi.pda[esi.observed_object["nta"]].define_ET = esi.pda[esi.observed_object["ntx"]].define_ET;}
    else {esi.pda[esi.observed_object["nta"]].define_ET = esi.pda[esi.observed_object["nt"]].define_ET;}
  }



  // FIXME: shouldn't that be 'max'??
  //define_eta[esi.observed_object["jet"]] = min(define_eta[esi.observed_object["jet"]],min(define_eta[esi.observed_object["ljet"]],define_eta[esi.observed_object["bjet"]]));



  // adapt define_eta (tjet-top/atop) !
  if (esi.pda[esi.observed_object["top"]].define_eta == 1.e99 && esi.pda[esi.observed_object["tjet"]].define_eta < 1.e99){esi.pda[esi.observed_object["top"]].define_eta = esi.pda[esi.observed_object["tjet"]].define_eta;}
  if (esi.pda[esi.observed_object["atop"]].define_eta == 1.e99 && esi.pda[esi.observed_object["tjet"]].define_eta < 1.e99){esi.pda[esi.observed_object["atop"]].define_eta = esi.pda[esi.observed_object["tjet"]].define_eta;}
  if (esi.pda[esi.observed_object["tjet"]].define_eta < esi.pda[esi.observed_object["top"]].define_eta &&
      esi.pda[esi.observed_object["tjet"]].define_eta < esi.pda[esi.observed_object["atop"]].define_eta){
    if (esi.pda[esi.observed_object["top"]].define_eta >= esi.pda[esi.observed_object["atop"]].define_eta){esi.pda[esi.observed_object["tjet"]].define_eta = esi.pda[esi.observed_object["atop"]].define_eta;}
    else {esi.pda[esi.observed_object["tjet"]].define_eta = esi.pda[esi.observed_object["top"]].define_eta;}
  }

  // adapt define_eta (bjet-bjet_b/bjet_bx) !
  if (esi.pda[esi.observed_object["bjet_b"]].define_eta == 1.e99 && esi.pda[esi.observed_object["bjet"]].define_eta < 1.e99){esi.pda[esi.observed_object["bjet_b"]].define_eta = esi.pda[esi.observed_object["bjet"]].define_eta;}
  if (esi.pda[esi.observed_object["bjet_bx"]].define_eta == 1.e99 && esi.pda[esi.observed_object["bjet"]].define_eta < 1.e99){esi.pda[esi.observed_object["bjet_bx"]].define_eta = esi.pda[esi.observed_object["bjet"]].define_eta;}
  if (esi.pda[esi.observed_object["bjet"]].define_eta < esi.pda[esi.observed_object["bjet_b"]].define_eta &&
      esi.pda[esi.observed_object["bjet"]].define_eta < esi.pda[esi.observed_object["bjet_bx"]].define_eta){
    if (esi.pda[esi.observed_object["bjet_b"]].define_eta >= esi.pda[esi.observed_object["bjet_bx"]].define_eta){esi.pda[esi.observed_object["bjet"]].define_eta = esi.pda[esi.observed_object["bjet_bx"]].define_eta;}
    else {esi.pda[esi.observed_object["bjet"]].define_eta = esi.pda[esi.observed_object["bjet_b"]].define_eta;}
  }
  /*
  // adapt define_eta (jet-ljet/bjet) !
  if (esi.pda[esi.observed_object["ljet"]].define_eta == 1.e99 && esi.pda[esi.observed_object["jet"]].define_eta < 1.e99){esi.pda[esi.observed_object["ljet"]].define_eta = esi.pda[esi.observed_object["jet"]].define_eta;}
  if (esi.pda[esi.observed_object["bjet"]].define_eta == 1.e99 && esi.pda[esi.observed_object["jet"]].define_eta < 1.e99){esi.pda[esi.observed_object["bjet"]].define_eta = esi.pda[esi.observed_object["jet"]].define_eta;}

  if (esi.pda[esi.observed_object["jet"]].define_eta == 0. && esi.pda[esi.observed_object["bjet"]].define_eta == 0. && esi.pda[esi.observed_object["ljet"]].define_eta > 0.){esi.pda[esi.observed_object["jet"]].define_eta = esi.pda[esi.observed_object["ljet"]].define_eta;}
  if (esi.pda[esi.observed_object["jet"]].define_eta == 0. && esi.pda[esi.observed_object["ljet"]].define_eta == 0. && esi.pda[esi.observed_object["bjet"]].define_eta > 0.){esi.pda[esi.observed_object["jet"]].define_eta = esi.pda[esi.observed_object["bjet"]].define_eta;}

  if (esi.pda[esi.observed_object["jet"]].define_eta < esi.pda[esi.observed_object["ljet"]].define_eta &&
      esi.pda[esi.observed_object["jet"]].define_eta < esi.pda[esi.observed_object["bjet"]].define_eta){
    if (esi.pda[esi.observed_object["ljet"]].define_eta >= esi.pda[esi.observed_object["bjet"]].define_eta){esi.pda[esi.observed_object["jet"]].define_eta = esi.pda[esi.observed_object["bjet"]].define_eta;}
    else {esi.pda[esi.observed_object["jet"]].define_eta = esi.pda[esi.observed_object["ljet"]].define_eta;}
  }
  */
  // adapt define_eta (lep-lm/lp) !
  if (esi.pda[esi.observed_object["lm"]].define_eta == 1.e99 && esi.pda[esi.observed_object["lep"]].define_eta < 1.e99){esi.pda[esi.observed_object["lm"]].define_eta = esi.pda[esi.observed_object["lep"]].define_eta;}
  if (esi.pda[esi.observed_object["lp"]].define_eta == 1.e99 && esi.pda[esi.observed_object["lep"]].define_eta < 1.e99){esi.pda[esi.observed_object["lp"]].define_eta = esi.pda[esi.observed_object["lep"]].define_eta;}
  if (esi.pda[esi.observed_object["lep"]].define_eta < esi.pda[esi.observed_object["lm"]].define_eta &&
      esi.pda[esi.observed_object["lep"]].define_eta < esi.pda[esi.observed_object["lp"]].define_eta){
    if (esi.pda[esi.observed_object["lm"]].define_eta >= esi.pda[esi.observed_object["lp"]].define_eta){esi.pda[esi.observed_object["lep"]].define_eta = esi.pda[esi.observed_object["lp"]].define_eta;}
    else {esi.pda[esi.observed_object["lep"]].define_eta = esi.pda[esi.observed_object["lm"]].define_eta;}
  }

  // adapt define_eta (lep-e/mu/tau) !
  if (esi.pda[esi.observed_object["e"]].define_eta == 1.e99 && esi.pda[esi.observed_object["lep"]].define_eta < 1.e99){esi.pda[esi.observed_object["e"]].define_eta = esi.pda[esi.observed_object["lep"]].define_eta;}
  if (esi.pda[esi.observed_object["mu"]].define_eta == 1.e99 && esi.pda[esi.observed_object["lep"]].define_eta < 1.e99){esi.pda[esi.observed_object["mu"]].define_eta = esi.pda[esi.observed_object["lep"]].define_eta;}
  if (esi.pda[esi.observed_object["tau"]].define_eta == 1.e99 && esi.pda[esi.observed_object["lep"]].define_eta < 1.e99){esi.pda[esi.observed_object["tau"]].define_eta = esi.pda[esi.observed_object["lep"]].define_eta;}
  if (esi.pda[esi.observed_object["lep"]].define_eta < esi.pda[esi.observed_object["e"]].define_eta &&
      esi.pda[esi.observed_object["lep"]].define_eta < esi.pda[esi.observed_object["mu"]].define_eta &&
      esi.pda[esi.observed_object["lep"]].define_eta < esi.pda[esi.observed_object["tau"]].define_eta){
    if (esi.pda[esi.observed_object["tau"]].define_eta >= esi.pda[esi.observed_object["e"]].define_eta &&
	esi.pda[esi.observed_object["mu"]].define_eta >= esi.pda[esi.observed_object["e"]].define_eta){esi.pda[esi.observed_object["lep"]].define_eta = esi.pda[esi.observed_object["e"]].define_eta;}
    else if (esi.pda[esi.observed_object["tau"]].define_eta >= esi.pda[esi.observed_object["mu"]].define_eta){esi.pda[esi.observed_object["lep"]].define_eta = esi.pda[esi.observed_object["mu"]].define_eta;}
    else {esi.pda[esi.observed_object["lep"]].define_eta = esi.pda[esi.observed_object["tau"]].define_eta;}
  }

  // adapt define_eta (e-em/ep) !
  if (esi.pda[esi.observed_object["em"]].define_eta == 1.e99 && esi.pda[esi.observed_object["e"]].define_eta < 1.e99){esi.pda[esi.observed_object["em"]].define_eta = esi.pda[esi.observed_object["e"]].define_eta;}
  if (esi.pda[esi.observed_object["ep"]].define_eta == 1.e99 && esi.pda[esi.observed_object["e"]].define_eta < 1.e99){esi.pda[esi.observed_object["ep"]].define_eta = esi.pda[esi.observed_object["e"]].define_eta;}
  if (esi.pda[esi.observed_object["e"]].define_eta < esi.pda[esi.observed_object["em"]].define_eta &&
      esi.pda[esi.observed_object["e"]].define_eta < esi.pda[esi.observed_object["ep"]].define_eta){
    if (esi.pda[esi.observed_object["em"]].define_eta > esi.pda[esi.observed_object["ep"]].define_eta){esi.pda[esi.observed_object["e"]].define_eta = esi.pda[esi.observed_object["ep"]].define_eta;}
    else {esi.pda[esi.observed_object["e"]].define_eta = esi.pda[esi.observed_object["em"]].define_eta;}
  }

  // adapt define_eta (mu-mum/mup) !
  if (esi.pda[esi.observed_object["mum"]].define_eta == 1.e99 && esi.pda[esi.observed_object["mu"]].define_eta < 1.e99){esi.pda[esi.observed_object["mum"]].define_eta = esi.pda[esi.observed_object["mu"]].define_eta;}
  if (esi.pda[esi.observed_object["mup"]].define_eta == 1.e99 && esi.pda[esi.observed_object["mu"]].define_eta < 1.e99){esi.pda[esi.observed_object["mup"]].define_eta = esi.pda[esi.observed_object["mu"]].define_eta;}
  if (esi.pda[esi.observed_object["mu"]].define_eta < esi.pda[esi.observed_object["mum"]].define_eta &&
      esi.pda[esi.observed_object["mu"]].define_eta < esi.pda[esi.observed_object["mup"]].define_eta){
    if (esi.pda[esi.observed_object["mum"]].define_eta > esi.pda[esi.observed_object["mup"]].define_eta){esi.pda[esi.observed_object["mu"]].define_eta = esi.pda[esi.observed_object["mup"]].define_eta;}
    else {esi.pda[esi.observed_object["mu"]].define_eta = esi.pda[esi.observed_object["mum"]].define_eta;}
  }

  // adapt define_eta (tau-taum/taup) !
  if (esi.pda[esi.observed_object["taum"]].define_eta == 1.e99 && esi.pda[esi.observed_object["tau"]].define_eta < 1.e99){esi.pda[esi.observed_object["taum"]].define_eta = esi.pda[esi.observed_object["tau"]].define_eta;}
  if (esi.pda[esi.observed_object["taup"]].define_eta == 1.e99 && esi.pda[esi.observed_object["tau"]].define_eta < 1.e99){esi.pda[esi.observed_object["taup"]].define_eta = esi.pda[esi.observed_object["tau"]].define_eta;}
  if (esi.pda[esi.observed_object["tau"]].define_eta < esi.pda[esi.observed_object["taum"]].define_eta &&
      esi.pda[esi.observed_object["tau"]].define_eta < esi.pda[esi.observed_object["taup"]].define_eta){
    if (esi.pda[esi.observed_object["taum"]].define_eta > esi.pda[esi.observed_object["taup"]].define_eta){esi.pda[esi.observed_object["tau"]].define_eta = esi.pda[esi.observed_object["taup"]].define_eta;}
    else {esi.pda[esi.observed_object["tau"]].define_eta = esi.pda[esi.observed_object["taum"]].define_eta;}
  }

  // adapt define_eta (nua-nu/nux) !
  if (esi.pda[esi.observed_object["nu"]].define_eta == 1.e99 && esi.pda[esi.observed_object["nua"]].define_eta < 1.e99){esi.pda[esi.observed_object["nu"]].define_eta = esi.pda[esi.observed_object["nua"]].define_eta;}
  if (esi.pda[esi.observed_object["nux"]].define_eta == 1.e99 && esi.pda[esi.observed_object["nua"]].define_eta < 1.e99){esi.pda[esi.observed_object["nux"]].define_eta = esi.pda[esi.observed_object["nua"]].define_eta;}
  if (esi.pda[esi.observed_object["nua"]].define_eta < esi.pda[esi.observed_object["nu"]].define_eta &&
      esi.pda[esi.observed_object["nua"]].define_eta < esi.pda[esi.observed_object["nux"]].define_eta){
    if (esi.pda[esi.observed_object["nu"]].define_eta >= esi.pda[esi.observed_object["nux"]].define_eta){esi.pda[esi.observed_object["nua"]].define_eta = esi.pda[esi.observed_object["nux"]].define_eta;}
    else {esi.pda[esi.observed_object["nua"]].define_eta = esi.pda[esi.observed_object["nu"]].define_eta;}
  }

  // adapt define_eta (nua-nea/nma/nta) !
  if (esi.pda[esi.observed_object["nea"]].define_eta == 1.e99 && esi.pda[esi.observed_object["nua"]].define_eta < 1.e99){esi.pda[esi.observed_object["nea"]].define_eta = esi.pda[esi.observed_object["nua"]].define_eta;}
  if (esi.pda[esi.observed_object["nma"]].define_eta == 1.e99 && esi.pda[esi.observed_object["nua"]].define_eta < 1.e99){esi.pda[esi.observed_object["nma"]].define_eta = esi.pda[esi.observed_object["nua"]].define_eta;}
  if (esi.pda[esi.observed_object["nta"]].define_eta == 1.e99 && esi.pda[esi.observed_object["nua"]].define_eta < 1.e99){esi.pda[esi.observed_object["nta"]].define_eta = esi.pda[esi.observed_object["nua"]].define_eta;}
  if (esi.pda[esi.observed_object["nua"]].define_eta < esi.pda[esi.observed_object["nea"]].define_eta &&
      esi.pda[esi.observed_object["nua"]].define_eta < esi.pda[esi.observed_object["nma"]].define_eta &&
      esi.pda[esi.observed_object["nua"]].define_eta < esi.pda[esi.observed_object["nta"]].define_eta){
    if (esi.pda[esi.observed_object["nta"]].define_eta >= esi.pda[esi.observed_object["nea"]].define_eta &&
	esi.pda[esi.observed_object["nma"]].define_eta >= esi.pda[esi.observed_object["nea"]].define_eta){esi.pda[esi.observed_object["nua"]].define_eta = esi.pda[esi.observed_object["nea"]].define_eta;}
    else if (esi.pda[esi.observed_object["nta"]].define_eta >= esi.pda[esi.observed_object["nma"]].define_eta){esi.pda[esi.observed_object["nua"]].define_eta = esi.pda[esi.observed_object["nma"]].define_eta;}
    else {esi.pda[esi.observed_object["nua"]].define_eta = esi.pda[esi.observed_object["nta"]].define_eta;}
  }

  // adapt define_eta (nea-ne/nex) !
  if (esi.pda[esi.observed_object["ne"]].define_eta == 1.e99 && esi.pda[esi.observed_object["nea"]].define_eta < 1.e99){esi.pda[esi.observed_object["ne"]].define_eta = esi.pda[esi.observed_object["nea"]].define_eta;}
  if (esi.pda[esi.observed_object["nex"]].define_eta == 1.e99 && esi.pda[esi.observed_object["nea"]].define_eta < 1.e99){esi.pda[esi.observed_object["nex"]].define_eta = esi.pda[esi.observed_object["nea"]].define_eta;}
  if (esi.pda[esi.observed_object["nea"]].define_eta < esi.pda[esi.observed_object["ne"]].define_eta &&
      esi.pda[esi.observed_object["nea"]].define_eta < esi.pda[esi.observed_object["nex"]].define_eta){
    if (esi.pda[esi.observed_object["ne"]].define_eta > esi.pda[esi.observed_object["nex"]].define_eta){esi.pda[esi.observed_object["nea"]].define_eta = esi.pda[esi.observed_object["nex"]].define_eta;}
    else {esi.pda[esi.observed_object["nea"]].define_eta = esi.pda[esi.observed_object["ne"]].define_eta;}
  }

  // adapt define_eta (nma-nm/nmx) !
  if (esi.pda[esi.observed_object["nm"]].define_eta == 1.e99 && esi.pda[esi.observed_object["nma"]].define_eta < 1.e99){esi.pda[esi.observed_object["nm"]].define_eta = esi.pda[esi.observed_object["nma"]].define_eta;}
  if (esi.pda[esi.observed_object["nmx"]].define_eta == 1.e99 && esi.pda[esi.observed_object["nma"]].define_eta < 1.e99){esi.pda[esi.observed_object["nmx"]].define_eta = esi.pda[esi.observed_object["nma"]].define_eta;}
  if (esi.pda[esi.observed_object["nma"]].define_eta < esi.pda[esi.observed_object["nm"]].define_eta &&
      esi.pda[esi.observed_object["nma"]].define_eta < esi.pda[esi.observed_object["nmx"]].define_eta){
    if (esi.pda[esi.observed_object["nm"]].define_eta > esi.pda[esi.observed_object["nmx"]].define_eta){esi.pda[esi.observed_object["nma"]].define_eta = esi.pda[esi.observed_object["nmx"]].define_eta;}
    else {esi.pda[esi.observed_object["nma"]].define_eta = esi.pda[esi.observed_object["nm"]].define_eta;}
  }

  // adapt define_eta (nta-nt/ntx) !
  if (esi.pda[esi.observed_object["nt"]].define_eta == 1.e99 && esi.pda[esi.observed_object["nta"]].define_eta < 1.e99){esi.pda[esi.observed_object["nt"]].define_eta = esi.pda[esi.observed_object["nta"]].define_eta;}
  if (esi.pda[esi.observed_object["ntx"]].define_eta == 1.e99 && esi.pda[esi.observed_object["nta"]].define_eta < 1.e99){esi.pda[esi.observed_object["ntx"]].define_eta = esi.pda[esi.observed_object["nta"]].define_eta;}
  if (esi.pda[esi.observed_object["nta"]].define_eta < esi.pda[esi.observed_object["nt"]].define_eta &&
      esi.pda[esi.observed_object["nta"]].define_eta < esi.pda[esi.observed_object["ntx"]].define_eta){
    if (esi.pda[esi.observed_object["nt"]].define_eta > esi.pda[esi.observed_object["ntx"]].define_eta){esi.pda[esi.observed_object["nta"]].define_eta = esi.pda[esi.observed_object["ntx"]].define_eta;}
    else {esi.pda[esi.observed_object["nta"]].define_eta = esi.pda[esi.observed_object["nt"]].define_eta;}
  }


  // adapt define_y (tjet-top/atop) !
  if (esi.pda[esi.observed_object["top"]].define_y == 1.e99 && esi.pda[esi.observed_object["tjet"]].define_y < 1.e99){esi.pda[esi.observed_object["top"]].define_y = esi.pda[esi.observed_object["tjet"]].define_y;}
  if (esi.pda[esi.observed_object["atop"]].define_y == 1.e99 && esi.pda[esi.observed_object["tjet"]].define_y < 1.e99){esi.pda[esi.observed_object["atop"]].define_y = esi.pda[esi.observed_object["tjet"]].define_y;}
  if (esi.pda[esi.observed_object["tjet"]].define_y < esi.pda[esi.observed_object["top"]].define_y &&
      esi.pda[esi.observed_object["tjet"]].define_y < esi.pda[esi.observed_object["atop"]].define_y){
    if (esi.pda[esi.observed_object["top"]].define_y >= esi.pda[esi.observed_object["atop"]].define_y){esi.pda[esi.observed_object["tjet"]].define_y = esi.pda[esi.observed_object["atop"]].define_y;}
    else {esi.pda[esi.observed_object["tjet"]].define_y = esi.pda[esi.observed_object["top"]].define_y;}
  }

  // adapt define_y (bjet-bjet_b/bjet_bx) !
  if (esi.pda[esi.observed_object["bjet_b"]].define_y == 1.e99 && esi.pda[esi.observed_object["bjet"]].define_y < 1.e99){esi.pda[esi.observed_object["bjet_b"]].define_y = esi.pda[esi.observed_object["bjet"]].define_y;}
  if (esi.pda[esi.observed_object["bjet_bx"]].define_y == 1.e99 && esi.pda[esi.observed_object["bjet"]].define_y < 1.e99){esi.pda[esi.observed_object["bjet_bx"]].define_y = esi.pda[esi.observed_object["bjet"]].define_y;}
  if (esi.pda[esi.observed_object["bjet"]].define_y < esi.pda[esi.observed_object["bjet_b"]].define_y &&
      esi.pda[esi.observed_object["bjet"]].define_y < esi.pda[esi.observed_object["bjet_bx"]].define_y){
    if (esi.pda[esi.observed_object["bjet_b"]].define_y >= esi.pda[esi.observed_object["bjet_bx"]].define_y){esi.pda[esi.observed_object["bjet"]].define_y = esi.pda[esi.observed_object["bjet_bx"]].define_y;}
    else {esi.pda[esi.observed_object["bjet"]].define_y = esi.pda[esi.observed_object["bjet_b"]].define_y;}
  }
  /*
  // adapt define_y (jet-ljet/bjet) !
  if (esi.pda[esi.observed_object["ljet"]].define_y == 1.e99 && esi.pda[esi.observed_object["jet"]].define_y < 1.e99){esi.pda[esi.observed_object["ljet"]].define_y = esi.pda[esi.observed_object["jet"]].define_y;}
  if (esi.pda[esi.observed_object["bjet"]].define_y == 1.e99 && esi.pda[esi.observed_object["jet"]].define_y < 1.e99){esi.pda[esi.observed_object["bjet"]].define_y = esi.pda[esi.observed_object["jet"]].define_y;}

  if (esi.pda[esi.observed_object["jet"]].define_y == 0. && esi.pda[esi.observed_object["bjet"]].define_y == 0. && esi.pda[esi.observed_object["ljet"]].define_y > 0.){esi.pda[esi.observed_object["jet"]].define_y = esi.pda[esi.observed_object["ljet"]].define_y;}
  if (esi.pda[esi.observed_object["jet"]].define_y == 0. && esi.pda[esi.observed_object["ljet"]].define_y == 0. && esi.pda[esi.observed_object["bjet"]].define_y > 0.){esi.pda[esi.observed_object["jet"]].define_y = esi.pda[esi.observed_object["bjet"]].define_y;}

  if (esi.pda[esi.observed_object["jet"]].define_y < esi.pda[esi.observed_object["ljet"]].define_y &&
      esi.pda[esi.observed_object["jet"]].define_y < esi.pda[esi.observed_object["bjet"]].define_y){
    if (esi.pda[esi.observed_object["ljet"]].define_y >= esi.pda[esi.observed_object["bjet"]].define_y){esi.pda[esi.observed_object["jet"]].define_y = esi.pda[esi.observed_object["bjet"]].define_y;}
    else {esi.pda[esi.observed_object["jet"]].define_y = esi.pda[esi.observed_object["ljet"]].define_y;}
  }
  */
  // adapt define_y (lep-lm/lp) !
  if (esi.pda[esi.observed_object["lm"]].define_y == 1.e99 && esi.pda[esi.observed_object["lep"]].define_y < 1.e99){esi.pda[esi.observed_object["lm"]].define_y = esi.pda[esi.observed_object["lep"]].define_y;}
  if (esi.pda[esi.observed_object["lp"]].define_y == 1.e99 && esi.pda[esi.observed_object["lep"]].define_y < 1.e99){esi.pda[esi.observed_object["lp"]].define_y = esi.pda[esi.observed_object["lep"]].define_y;}
  if (esi.pda[esi.observed_object["lep"]].define_y < esi.pda[esi.observed_object["lm"]].define_y &&
      esi.pda[esi.observed_object["lep"]].define_y < esi.pda[esi.observed_object["lp"]].define_y){
    if (esi.pda[esi.observed_object["lm"]].define_y >= esi.pda[esi.observed_object["lp"]].define_y){esi.pda[esi.observed_object["lep"]].define_y = esi.pda[esi.observed_object["lp"]].define_y;}
    else {esi.pda[esi.observed_object["lep"]].define_y = esi.pda[esi.observed_object["lm"]].define_y;}
  }

  // adapt define_y (lep-e/mu/tau) !
  if (esi.pda[esi.observed_object["e"]].define_y == 1.e99 && esi.pda[esi.observed_object["lep"]].define_y < 1.e99){esi.pda[esi.observed_object["e"]].define_y = esi.pda[esi.observed_object["lep"]].define_y;}
  if (esi.pda[esi.observed_object["mu"]].define_y == 1.e99 && esi.pda[esi.observed_object["lep"]].define_y < 1.e99){esi.pda[esi.observed_object["mu"]].define_y = esi.pda[esi.observed_object["lep"]].define_y;}
  if (esi.pda[esi.observed_object["tau"]].define_y == 1.e99 && esi.pda[esi.observed_object["lep"]].define_y < 1.e99){esi.pda[esi.observed_object["tau"]].define_y = esi.pda[esi.observed_object["lep"]].define_y;}
  if (esi.pda[esi.observed_object["lep"]].define_y < esi.pda[esi.observed_object["e"]].define_y &&
      esi.pda[esi.observed_object["lep"]].define_y < esi.pda[esi.observed_object["mu"]].define_y &&
      esi.pda[esi.observed_object["lep"]].define_y < esi.pda[esi.observed_object["tau"]].define_y){
    if (esi.pda[esi.observed_object["tau"]].define_y >= esi.pda[esi.observed_object["e"]].define_y &&
	esi.pda[esi.observed_object["mu"]].define_y >= esi.pda[esi.observed_object["e"]].define_y){esi.pda[esi.observed_object["lep"]].define_y = esi.pda[esi.observed_object["e"]].define_y;}
    else if (esi.pda[esi.observed_object["tau"]].define_y >= esi.pda[esi.observed_object["mu"]].define_y){esi.pda[esi.observed_object["lep"]].define_y = esi.pda[esi.observed_object["mu"]].define_y;}
    else {esi.pda[esi.observed_object["lep"]].define_y = esi.pda[esi.observed_object["tau"]].define_y;}
  }

  // adapt define_y (e-em/ep) !
  if (esi.pda[esi.observed_object["em"]].define_y == 1.e99 && esi.pda[esi.observed_object["e"]].define_y < 1.e99){esi.pda[esi.observed_object["em"]].define_y = esi.pda[esi.observed_object["e"]].define_y;}
  if (esi.pda[esi.observed_object["ep"]].define_y == 1.e99 && esi.pda[esi.observed_object["e"]].define_y < 1.e99){esi.pda[esi.observed_object["ep"]].define_y = esi.pda[esi.observed_object["e"]].define_y;}
  if (esi.pda[esi.observed_object["e"]].define_y < esi.pda[esi.observed_object["em"]].define_y &&
      esi.pda[esi.observed_object["e"]].define_y < esi.pda[esi.observed_object["ep"]].define_y){
    if (esi.pda[esi.observed_object["em"]].define_y > esi.pda[esi.observed_object["ep"]].define_y){esi.pda[esi.observed_object["e"]].define_y = esi.pda[esi.observed_object["ep"]].define_y;}
    else {esi.pda[esi.observed_object["e"]].define_y = esi.pda[esi.observed_object["em"]].define_y;}
  }

  // adapt define_y (mu-mum/mup) !
  if (esi.pda[esi.observed_object["mum"]].define_y == 1.e99 && esi.pda[esi.observed_object["mu"]].define_y < 1.e99){esi.pda[esi.observed_object["mum"]].define_y = esi.pda[esi.observed_object["mu"]].define_y;}
  if (esi.pda[esi.observed_object["mup"]].define_y == 1.e99 && esi.pda[esi.observed_object["mu"]].define_y < 1.e99){esi.pda[esi.observed_object["mup"]].define_y = esi.pda[esi.observed_object["mu"]].define_y;}
  if (esi.pda[esi.observed_object["mu"]].define_y < esi.pda[esi.observed_object["mum"]].define_y &&
      esi.pda[esi.observed_object["mu"]].define_y < esi.pda[esi.observed_object["mup"]].define_y){
    if (esi.pda[esi.observed_object["mum"]].define_y > esi.pda[esi.observed_object["mup"]].define_y){esi.pda[esi.observed_object["mu"]].define_y = esi.pda[esi.observed_object["mup"]].define_y;}
    else {esi.pda[esi.observed_object["mu"]].define_y = esi.pda[esi.observed_object["mum"]].define_y;}
  }

  // adapt define_y (tau-taum/taup) !
  if (esi.pda[esi.observed_object["taum"]].define_y == 1.e99 && esi.pda[esi.observed_object["tau"]].define_y < 1.e99){esi.pda[esi.observed_object["taum"]].define_y = esi.pda[esi.observed_object["tau"]].define_y;}
  if (esi.pda[esi.observed_object["taup"]].define_y == 1.e99 && esi.pda[esi.observed_object["tau"]].define_y < 1.e99){esi.pda[esi.observed_object["taup"]].define_y = esi.pda[esi.observed_object["tau"]].define_y;}
  if (esi.pda[esi.observed_object["tau"]].define_y < esi.pda[esi.observed_object["taum"]].define_y &&
      esi.pda[esi.observed_object["tau"]].define_y < esi.pda[esi.observed_object["taup"]].define_y){
    if (esi.pda[esi.observed_object["taum"]].define_y > esi.pda[esi.observed_object["taup"]].define_y){esi.pda[esi.observed_object["tau"]].define_y = esi.pda[esi.observed_object["taup"]].define_y;}
    else {esi.pda[esi.observed_object["tau"]].define_y = esi.pda[esi.observed_object["taum"]].define_y;}
  }

  // adapt define_y (nua-nu/nux) !
  if (esi.pda[esi.observed_object["nu"]].define_y == 1.e99 && esi.pda[esi.observed_object["nua"]].define_y < 1.e99){esi.pda[esi.observed_object["nu"]].define_y = esi.pda[esi.observed_object["nua"]].define_y;}
  if (esi.pda[esi.observed_object["nux"]].define_y == 1.e99 && esi.pda[esi.observed_object["nua"]].define_y < 1.e99){esi.pda[esi.observed_object["nux"]].define_y = esi.pda[esi.observed_object["nua"]].define_y;}
  if (esi.pda[esi.observed_object["nua"]].define_y < esi.pda[esi.observed_object["nu"]].define_y &&
      esi.pda[esi.observed_object["nua"]].define_y < esi.pda[esi.observed_object["nux"]].define_y){
    if (esi.pda[esi.observed_object["nu"]].define_y >= esi.pda[esi.observed_object["nux"]].define_y){esi.pda[esi.observed_object["nua"]].define_y = esi.pda[esi.observed_object["nux"]].define_y;}
    else {esi.pda[esi.observed_object["nua"]].define_y = esi.pda[esi.observed_object["nu"]].define_y;}
  }

  // adapt define_y (nua-nea/nma/nta) !
  if (esi.pda[esi.observed_object["nea"]].define_y == 1.e99 && esi.pda[esi.observed_object["nua"]].define_y < 1.e99){esi.pda[esi.observed_object["nea"]].define_y = esi.pda[esi.observed_object["nua"]].define_y;}
  if (esi.pda[esi.observed_object["nma"]].define_y == 1.e99 && esi.pda[esi.observed_object["nua"]].define_y < 1.e99){esi.pda[esi.observed_object["nma"]].define_y = esi.pda[esi.observed_object["nua"]].define_y;}
  if (esi.pda[esi.observed_object["nta"]].define_y == 1.e99 && esi.pda[esi.observed_object["nua"]].define_y < 1.e99){esi.pda[esi.observed_object["nta"]].define_y = esi.pda[esi.observed_object["nua"]].define_y;}
  if (esi.pda[esi.observed_object["nua"]].define_y < esi.pda[esi.observed_object["nea"]].define_y &&
      esi.pda[esi.observed_object["nua"]].define_y < esi.pda[esi.observed_object["nma"]].define_y &&
      esi.pda[esi.observed_object["nua"]].define_y < esi.pda[esi.observed_object["nta"]].define_y){
    if (esi.pda[esi.observed_object["nta"]].define_y >= esi.pda[esi.observed_object["nea"]].define_y &&
	esi.pda[esi.observed_object["nma"]].define_y >= esi.pda[esi.observed_object["nea"]].define_y){esi.pda[esi.observed_object["nua"]].define_y = esi.pda[esi.observed_object["nea"]].define_y;}
    else if (esi.pda[esi.observed_object["nta"]].define_y >= esi.pda[esi.observed_object["nma"]].define_y){esi.pda[esi.observed_object["nua"]].define_y = esi.pda[esi.observed_object["nma"]].define_y;}
    else {esi.pda[esi.observed_object["nua"]].define_y = esi.pda[esi.observed_object["nta"]].define_y;}
  }

  // adapt define_y (nea-ne/nex) !
  if (esi.pda[esi.observed_object["ne"]].define_y == 1.e99 && esi.pda[esi.observed_object["nea"]].define_y < 1.e99){esi.pda[esi.observed_object["ne"]].define_y = esi.pda[esi.observed_object["nea"]].define_y;}
  if (esi.pda[esi.observed_object["nex"]].define_y == 1.e99 && esi.pda[esi.observed_object["nea"]].define_y < 1.e99){esi.pda[esi.observed_object["nex"]].define_y = esi.pda[esi.observed_object["nea"]].define_y;}
  if (esi.pda[esi.observed_object["nea"]].define_y < esi.pda[esi.observed_object["ne"]].define_y &&
      esi.pda[esi.observed_object["nea"]].define_y < esi.pda[esi.observed_object["nex"]].define_y){
    if (esi.pda[esi.observed_object["ne"]].define_y > esi.pda[esi.observed_object["nex"]].define_y){esi.pda[esi.observed_object["nea"]].define_y = esi.pda[esi.observed_object["nex"]].define_y;}
    else {esi.pda[esi.observed_object["nea"]].define_y = esi.pda[esi.observed_object["ne"]].define_y;}
  }

  // adapt define_y (nma-nm/nmx) !
  if (esi.pda[esi.observed_object["nm"]].define_y == 1.e99 && esi.pda[esi.observed_object["nma"]].define_y < 1.e99){esi.pda[esi.observed_object["nm"]].define_y = esi.pda[esi.observed_object["nma"]].define_y;}
  if (esi.pda[esi.observed_object["nmx"]].define_y == 1.e99 && esi.pda[esi.observed_object["nma"]].define_y < 1.e99){esi.pda[esi.observed_object["nmx"]].define_y = esi.pda[esi.observed_object["nma"]].define_y;}
  if (esi.pda[esi.observed_object["nma"]].define_y < esi.pda[esi.observed_object["nm"]].define_y &&
      esi.pda[esi.observed_object["nma"]].define_y < esi.pda[esi.observed_object["nmx"]].define_y){
    if (esi.pda[esi.observed_object["nm"]].define_y > esi.pda[esi.observed_object["nmx"]].define_y){esi.pda[esi.observed_object["nma"]].define_y = esi.pda[esi.observed_object["nmx"]].define_y;}
    else {esi.pda[esi.observed_object["nma"]].define_y = esi.pda[esi.observed_object["nm"]].define_y;}
  }

  // adapt define_y (nta-nt/ntx) !
  if (esi.pda[esi.observed_object["nt"]].define_y == 1.e99 && esi.pda[esi.observed_object["nta"]].define_y < 1.e99){esi.pda[esi.observed_object["nt"]].define_y = esi.pda[esi.observed_object["nta"]].define_y;}
  if (esi.pda[esi.observed_object["ntx"]].define_y == 1.e99 && esi.pda[esi.observed_object["nta"]].define_y < 1.e99){esi.pda[esi.observed_object["ntx"]].define_y = esi.pda[esi.observed_object["nta"]].define_y;}
  if (esi.pda[esi.observed_object["nta"]].define_y < esi.pda[esi.observed_object["nt"]].define_y &&
      esi.pda[esi.observed_object["nta"]].define_y < esi.pda[esi.observed_object["ntx"]].define_y){
    if (esi.pda[esi.observed_object["nt"]].define_y > esi.pda[esi.observed_object["ntx"]].define_y){esi.pda[esi.observed_object["nta"]].define_y = esi.pda[esi.observed_object["ntx"]].define_y;}
    else {esi.pda[esi.observed_object["nta"]].define_y = esi.pda[esi.observed_object["nt"]].define_y;}
  }


  // adapt n_observed_min (tjet-top/atop) !
  if (esi.pda[esi.observed_object["tjet"]].n_observed_min < esi.pda[esi.observed_object["top"]].n_observed_min + esi.pda[esi.observed_object["atop"]].n_observed_min){esi.pda[esi.observed_object["tjet"]].n_observed_min = esi.pda[esi.observed_object["top"]].n_observed_min + esi.pda[esi.observed_object["atop"]].n_observed_min;}

  // adapt n_observed_min (bjet-bjet_b/bjet_bx) !
  if (esi.pda[esi.observed_object["bjet"]].n_observed_min < esi.pda[esi.observed_object["bjet_b"]].n_observed_min + esi.pda[esi.observed_object["bjet_bx"]].n_observed_min){esi.pda[esi.observed_object["bjet"]].n_observed_min = esi.pda[esi.observed_object["bjet_b"]].n_observed_min + esi.pda[esi.observed_object["bjet_bx"]].n_observed_min;}

  // adapt n_observed_min (jets-ljets/bjets) !
  ///!X!  if (esi.pda[esi.observed_object["jet"]].n_observed_min < esi.pda[esi.observed_object["ljet"]].n_observed_min + esi.pda[esi.observed_object["bjet"]].n_observed_min){esi.pda[esi.observed_object["jet"]].n_observed_min = esi.pda[esi.observed_object["ljet"]].n_observed_min + esi.pda[esi.observed_object["bjet"]].n_observed_min;}

  // adapt n_observed_min (lm-em/mum/taum) !
  if (esi.pda[esi.observed_object["lm"]].n_observed_min < esi.pda[esi.observed_object["em"]].n_observed_min + esi.pda[esi.observed_object["mum"]].n_observed_min + esi.pda[esi.observed_object["taum"]].n_observed_min){
    esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["em"]].n_observed_min + esi.pda[esi.observed_object["mum"]].n_observed_min + esi.pda[esi.observed_object["taum"]].n_observed_min;
  }
  // adapt n_observed_min (lp-ep/mup/taup) !
  if (esi.pda[esi.observed_object["lp"]].n_observed_min < esi.pda[esi.observed_object["ep"]].n_observed_min + esi.pda[esi.observed_object["mup"]].n_observed_min + esi.pda[esi.observed_object["taup"]].n_observed_min){
    esi.pda[esi.observed_object["lp"]].n_observed_min = esi.pda[esi.observed_object["ep"]].n_observed_min + esi.pda[esi.observed_object["mup"]].n_observed_min + esi.pda[esi.observed_object["taup"]].n_observed_min;
  }
  // adapt n_observed_min (lep-lm/lp) !
  if (esi.pda[esi.observed_object["lep"]].n_observed_min < esi.pda[esi.observed_object["lm"]].n_observed_min + esi.pda[esi.observed_object["lp"]].n_observed_min){
    esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["lm"]].n_observed_min + esi.pda[esi.observed_object["lp"]].n_observed_min;
  }

  // adapt n_observed_min (e-ep/em) !
  if (esi.pda[esi.observed_object["e"]].n_observed_min < esi.pda[esi.observed_object["em"]].n_observed_min + esi.pda[esi.observed_object["ep"]].n_observed_min){
    esi.pda[esi.observed_object["e"]].n_observed_min = esi.pda[esi.observed_object["em"]].n_observed_min + esi.pda[esi.observed_object["ep"]].n_observed_min;
  }
  // adapt n_observed_min (mu-mup/mum) !
  if (esi.pda[esi.observed_object["mu"]].n_observed_min < esi.pda[esi.observed_object["mum"]].n_observed_min + esi.pda[esi.observed_object["mup"]].n_observed_min){
    esi.pda[esi.observed_object["mu"]].n_observed_min = esi.pda[esi.observed_object["mum"]].n_observed_min + esi.pda[esi.observed_object["mup"]].n_observed_min;
  }
  // adapt n_observed_min (tau-taup/taum) !
  if (esi.pda[esi.observed_object["tau"]].n_observed_min < esi.pda[esi.observed_object["taum"]].n_observed_min + esi.pda[esi.observed_object["taup"]].n_observed_min){
    esi.pda[esi.observed_object["tau"]].n_observed_min = esi.pda[esi.observed_object["taum"]].n_observed_min + esi.pda[esi.observed_object["taup"]].n_observed_min;
  }

  // adapt n_observed_min (lep-e/mu/tau) !
  if (esi.pda[esi.observed_object["lep"]].n_observed_min < esi.pda[esi.observed_object["e"]].n_observed_min + esi.pda[esi.observed_object["mu"]].n_observed_min + esi.pda[esi.observed_object["tau"]].n_observed_min){
    esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["e"]].n_observed_min + esi.pda[esi.observed_object["mu"]].n_observed_min + esi.pda[esi.observed_object["tau"]].n_observed_min;
  }
  // adapt n_observed_min (nu-ne/nm/nt) !
  if (esi.pda[esi.observed_object["nu"]].n_observed_min < esi.pda[esi.observed_object["ne"]].n_observed_min + esi.pda[esi.observed_object["nm"]].n_observed_min + esi.pda[esi.observed_object["nt"]].n_observed_min){
    esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["ne"]].n_observed_min + esi.pda[esi.observed_object["nm"]].n_observed_min + esi.pda[esi.observed_object["nt"]].n_observed_min;
  }
  // adapt n_observed_min (nux-nex/nmx/ntx) !
  if (esi.pda[esi.observed_object["nux"]].n_observed_min < esi.pda[esi.observed_object["nex"]].n_observed_min + esi.pda[esi.observed_object["nmx"]].n_observed_min + esi.pda[esi.observed_object["ntx"]].n_observed_min){
    esi.pda[esi.observed_object["nux"]].n_observed_min = esi.pda[esi.observed_object["nex"]].n_observed_min + esi.pda[esi.observed_object["nmx"]].n_observed_min + esi.pda[esi.observed_object["ntx"]].n_observed_min;
  }

  // adapt n_observed_min (nua-nu/nux) !
  if (esi.pda[esi.observed_object["nua"]].n_observed_min < esi.pda[esi.observed_object["nu"]].n_observed_min + esi.pda[esi.observed_object["nux"]].n_observed_min){
    esi.pda[esi.observed_object["nua"]].n_observed_min = esi.pda[esi.observed_object["nu"]].n_observed_min + esi.pda[esi.observed_object["nux"]].n_observed_min;
  }

  // adapt n_observed_min (nea-nex/ne) !
  if (esi.pda[esi.observed_object["nea"]].n_observed_min < esi.pda[esi.observed_object["ne"]].n_observed_min + esi.pda[esi.observed_object["nex"]].n_observed_min){
    esi.pda[esi.observed_object["nea"]].n_observed_min = esi.pda[esi.observed_object["ne"]].n_observed_min + esi.pda[esi.observed_object["nex"]].n_observed_min;
  }
  // adapt n_observed_min (nma-nmx/nm) !
  if (esi.pda[esi.observed_object["nma"]].n_observed_min < esi.pda[esi.observed_object["nm"]].n_observed_min + esi.pda[esi.observed_object["nmx"]].n_observed_min){
    esi.pda[esi.observed_object["nma"]].n_observed_min = esi.pda[esi.observed_object["nm"]].n_observed_min + esi.pda[esi.observed_object["nmx"]].n_observed_min;
  }
  // adapt n_observed_min (nta-ntx/nt) !
  if (esi.pda[esi.observed_object["nta"]].n_observed_min < esi.pda[esi.observed_object["nt"]].n_observed_min + esi.pda[esi.observed_object["ntx"]].n_observed_min){
    esi.pda[esi.observed_object["nta"]].n_observed_min = esi.pda[esi.observed_object["nt"]].n_observed_min + esi.pda[esi.observed_object["ntx"]].n_observed_min;
  }

  // adapt n_observed_min (nua-nea/nma/nta) !
  if (esi.pda[esi.observed_object["nua"]].n_observed_min < esi.pda[esi.observed_object["nea"]].n_observed_min + esi.pda[esi.observed_object["nma"]].n_observed_min + esi.pda[esi.observed_object["nta"]].n_observed_min){
    esi.pda[esi.observed_object["nua"]].n_observed_min = esi.pda[esi.observed_object["nea"]].n_observed_min + esi.pda[esi.observed_object["nma"]].n_observed_min + esi.pda[esi.observed_object["nta"]].n_observed_min;
  }


  // adapt n_observed_max (tjet-top/atop) !
  if (esi.pda[esi.observed_object["tjet"]].n_observed_max > esi.pda[esi.observed_object["top"]].n_observed_max + esi.pda[esi.observed_object["atop"]].n_observed_max){esi.pda[esi.observed_object["tjet"]].n_observed_max = esi.pda[esi.observed_object["top"]].n_observed_max + esi.pda[esi.observed_object["atop"]].n_observed_max;}

  // adapt n_observed_max (bjet-bjet_b/bjet_bx) !
  if (esi.pda[esi.observed_object["bjet"]].n_observed_max > esi.pda[esi.observed_object["bjet_b"]].n_observed_max + esi.pda[esi.observed_object["bjet_bx"]].n_observed_max){esi.pda[esi.observed_object["bjet"]].n_observed_max = esi.pda[esi.observed_object["bjet_b"]].n_observed_max + esi.pda[esi.observed_object["bjet_bx"]].n_observed_max;}
  // !!! xxxxx
  if (esi.pda[esi.observed_object["bjet"]].n_observed_max < esi.pda[esi.observed_object["bjet_b"]].n_observed_max){esi.pda[esi.observed_object["bjet_b"]].n_observed_max = esi.pda[esi.observed_object["bjet"]].n_observed_max;}
  if (esi.pda[esi.observed_object["bjet"]].n_observed_max < esi.pda[esi.observed_object["bjet_bx"]].n_observed_max){esi.pda[esi.observed_object["bjet_bx"]].n_observed_max = esi.pda[esi.observed_object["bjet"]].n_observed_max;}

  // adapt n_observed_max (jets-ljets/bjets) !
  ///!X!  if (esi.pda[esi.observed_object["jet"]].n_observed_max > esi.pda[esi.observed_object["ljet"]].n_observed_max + esi.pda[esi.observed_object["bjet"]].n_observed_max){esi.pda[esi.observed_object["jet"]].n_observed_max = esi.pda[esi.observed_object["ljet"]].n_observed_max + esi.pda[esi.observed_object["bjet"]].n_observed_max;}

  // adapt n_observed_max (lm-em/mum/taum) !
  if (esi.pda[esi.observed_object["lm"]].n_observed_max > esi.pda[esi.observed_object["em"]].n_observed_max + esi.pda[esi.observed_object["mum"]].n_observed_max + esi.pda[esi.observed_object["taum"]].n_observed_max){
    esi.pda[esi.observed_object["lm"]].n_observed_max = esi.pda[esi.observed_object["em"]].n_observed_max + esi.pda[esi.observed_object["mum"]].n_observed_max + esi.pda[esi.observed_object["taum"]].n_observed_max;
  }
  // adapt n_observed_max (lp-ep/mup/taup) !
  if (esi.pda[esi.observed_object["lp"]].n_observed_max > esi.pda[esi.observed_object["ep"]].n_observed_max + esi.pda[esi.observed_object["mup"]].n_observed_max + esi.pda[esi.observed_object["taup"]].n_observed_max){
    esi.pda[esi.observed_object["lp"]].n_observed_max = esi.pda[esi.observed_object["ep"]].n_observed_max + esi.pda[esi.observed_object["mup"]].n_observed_max + esi.pda[esi.observed_object["taup"]].n_observed_max;
  }

  // adapt n_observed_max (lep-lm/lp) !
  if (esi.pda[esi.observed_object["lep"]].n_observed_max > esi.pda[esi.observed_object["lm"]].n_observed_max + esi.pda[esi.observed_object["lp"]].n_observed_max){
    esi.pda[esi.observed_object["lep"]].n_observed_max = esi.pda[esi.observed_object["lm"]].n_observed_max + esi.pda[esi.observed_object["lp"]].n_observed_max;
  }

  // adapt n_observed_max (e-ep/em) !
  if (esi.pda[esi.observed_object["e"]].n_observed_max > esi.pda[esi.observed_object["em"]].n_observed_max + esi.pda[esi.observed_object["ep"]].n_observed_max){
    esi.pda[esi.observed_object["e"]].n_observed_max = esi.pda[esi.observed_object["em"]].n_observed_max + esi.pda[esi.observed_object["ep"]].n_observed_max;
  }
  // adapt n_observed_max (mu-mup/mum) !
  if (esi.pda[esi.observed_object["mu"]].n_observed_max > esi.pda[esi.observed_object["mum"]].n_observed_max + esi.pda[esi.observed_object["mup"]].n_observed_max){
    esi.pda[esi.observed_object["mu"]].n_observed_max = esi.pda[esi.observed_object["mum"]].n_observed_max + esi.pda[esi.observed_object["mup"]].n_observed_max;
  }
  // adapt n_observed_max (tau-taup/taum) !
  if (esi.pda[esi.observed_object["tau"]].n_observed_max > esi.pda[esi.observed_object["taum"]].n_observed_max + esi.pda[esi.observed_object["taup"]].n_observed_max){
    esi.pda[esi.observed_object["tau"]].n_observed_max = esi.pda[esi.observed_object["taum"]].n_observed_max + esi.pda[esi.observed_object["taup"]].n_observed_max;
  }

  // adapt n_observed_max (lep-e/mu/tau) !
  if (esi.pda[esi.observed_object["lep"]].n_observed_max > esi.pda[esi.observed_object["e"]].n_observed_max + esi.pda[esi.observed_object["mu"]].n_observed_max + esi.pda[esi.observed_object["tau"]].n_observed_max){
    esi.pda[esi.observed_object["lep"]].n_observed_max = esi.pda[esi.observed_object["e"]].n_observed_max + esi.pda[esi.observed_object["mu"]].n_observed_max + esi.pda[esi.observed_object["tau"]].n_observed_max;
  }
  // adapt n_observed_max (nu-ne/nm/nt) !
  if (esi.pda[esi.observed_object["nu"]].n_observed_max > esi.pda[esi.observed_object["ne"]].n_observed_max + esi.pda[esi.observed_object["nm"]].n_observed_max + esi.pda[esi.observed_object["nt"]].n_observed_max){
    esi.pda[esi.observed_object["nu"]].n_observed_max = esi.pda[esi.observed_object["ne"]].n_observed_max + esi.pda[esi.observed_object["nm"]].n_observed_max + esi.pda[esi.observed_object["nt"]].n_observed_max;
  }
  // adapt n_observed_max (nux-nex/nmx/ntx) !
  if (esi.pda[esi.observed_object["nux"]].n_observed_max > esi.pda[esi.observed_object["nex"]].n_observed_max + esi.pda[esi.observed_object["nmx"]].n_observed_max + esi.pda[esi.observed_object["ntx"]].n_observed_max){
    esi.pda[esi.observed_object["nux"]].n_observed_max = esi.pda[esi.observed_object["nex"]].n_observed_max + esi.pda[esi.observed_object["nmx"]].n_observed_max + esi.pda[esi.observed_object["ntx"]].n_observed_max;
  }

  // adapt n_observed_max (nua-nu/nux) !
  if (esi.pda[esi.observed_object["nua"]].n_observed_max > esi.pda[esi.observed_object["nu"]].n_observed_max + esi.pda[esi.observed_object["nux"]].n_observed_max){
    esi.pda[esi.observed_object["nua"]].n_observed_max = esi.pda[esi.observed_object["nu"]].n_observed_max + esi.pda[esi.observed_object["nux"]].n_observed_max;
  }
  // adapt n_observed_max (nea-nex/ne) !
  if (esi.pda[esi.observed_object["nea"]].n_observed_max > esi.pda[esi.observed_object["ne"]].n_observed_max + esi.pda[esi.observed_object["nex"]].n_observed_max){
    esi.pda[esi.observed_object["nea"]].n_observed_max = esi.pda[esi.observed_object["ne"]].n_observed_max + esi.pda[esi.observed_object["nex"]].n_observed_max;
  }
  // adapt n_observed_max (nma-nmx/nm) !
  if (esi.pda[esi.observed_object["nma"]].n_observed_max > esi.pda[esi.observed_object["nm"]].n_observed_max + esi.pda[esi.observed_object["nmx"]].n_observed_max){
    esi.pda[esi.observed_object["nma"]].n_observed_max = esi.pda[esi.observed_object["nm"]].n_observed_max + esi.pda[esi.observed_object["nmx"]].n_observed_max;
  }
  // adapt n_observed_max (nta-ntx/nt) !
  if (esi.pda[esi.observed_object["nta"]].n_observed_max > esi.pda[esi.observed_object["nt"]].n_observed_max + esi.pda[esi.observed_object["ntx"]].n_observed_max){
    esi.pda[esi.observed_object["nta"]].n_observed_max = esi.pda[esi.observed_object["nt"]].n_observed_max + esi.pda[esi.observed_object["ntx"]].n_observed_max;
  }

  // adapt n_observed_max (nua-nea/nma/nta) !
  if (esi.pda[esi.observed_object["nua"]].n_observed_max > esi.pda[esi.observed_object["nea"]].n_observed_max + esi.pda[esi.observed_object["nma"]].n_observed_max + esi.pda[esi.observed_object["nta"]].n_observed_max){
    esi.pda[esi.observed_object["nua"]].n_observed_max = esi.pda[esi.observed_object["nea"]].n_observed_max + esi.pda[esi.observed_object["nma"]].n_observed_max + esi.pda[esi.observed_object["nta"]].n_observed_max;
  }


  logger << LOG_DEBUG << endl << "Object definition after internal correlations: " << endl;
  logger << LOG_DEBUG << setw(20) << "object" << setw(20) << "n_observed_min" << setw(20) << "n_observed_max" << setw(20) << "define_pT" << setw(20) << "define_ET" << setw(20) << "define_eta"  << setw(20) << "define_y"<< endl;
  for (int i = 1; i < esi.object_list.size(); i++){
    //    logger << LOG_DEBUG << setw(20) << esi.object_list[i] << setw(20) << n_observed_min[i] << setw(20) << n_observed_max[i] << setw(20) << define_pT[i] << setw(20) << define_ET[i] << setw(20) << define_eta[i] << setw(20) << define_y[i] << endl;
    logger << LOG_DEBUG << setw(5) << i << setw(20) << esi.object_list[i] << setw(20) << esi.pda[i].n_observed_min << setw(20) << esi.pda[i].n_observed_max << setw(20) << esi.pda[i].define_pT << setw(20) << esi.pda[i].define_ET << setw(20) << esi.pda[i].define_eta << setw(20) << esi.pda[i].define_y << endl;
  }
  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::determine_object_partonlevel(){
  static  Logger logger("observable_set::determine_object_partonlevel");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // !!! Should be modified such that user input is never overwritten !!!

  // adapt n_observed_max trivially to parton-level particles
  // needs to be modified if partons can change their categories 
  // (like mis-identification of bjets as ljets etc.)
  for (int i_p = 1; i_p < esi.object_list.size(); i_p++){
    if (esi.pda[i_p].n_partonlevel < esi.pda[i_p].n_observed_max){
      if (esi.object_list[i_p] == "ljet"){
	if (esi.pda[esi.observed_object["jet"]].n_partonlevel < esi.pda[i_p].n_observed_max){esi.pda[i_p].n_observed_max = esi.pda[esi.observed_object["jet"]].n_partonlevel;}
      }
      else if (esi.object_category[i_p] == -1){logger << LOG_DEBUG << esi.object_list[i_p] << " is a user-defined particle -> no modification." << endl;}
      else {esi.pda[i_p].n_observed_max = esi.pda[i_p].n_partonlevel;}
    }
  }
  //  if (esi.pda[esi.observed_object["jet"]] < esi.pda[esi.observed_object["jet"]]){esi.pda[esi.observed_object["jet"]] = esi.pda[esi.observed_object["jet"]].n_observed_min;}


  // adapt n_observed_max to parton-level particles (tjet-top/atop)
  if (esi.pda[esi.observed_object["atop"]].n_partonlevel < esi.pda[esi.observed_object["atop"]].n_observed_max){esi.pda[esi.observed_object["atop"]].n_observed_max = esi.pda[esi.observed_object["atop"]].n_partonlevel;}
  if (esi.pda[esi.observed_object["top"]].n_partonlevel < esi.pda[esi.observed_object["top"]].n_observed_max){esi.pda[esi.observed_object["top"]].n_observed_max = esi.pda[esi.observed_object["top"]].n_partonlevel;}

  if (esi.pda[esi.observed_object["tjet"]].n_partonlevel == esi.pda[esi.observed_object["tjet"]].n_observed_max){
    esi.pda[esi.observed_object["atop"]].n_observed_max = esi.pda[esi.observed_object["atop"]].n_partonlevel;
    esi.pda[esi.observed_object["top"]].n_observed_max = esi.pda[esi.observed_object["top"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["top"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["atop"]].n_observed_max <= esi.pda[esi.observed_object["tjet"]].n_observed_max){esi.pda[esi.observed_object["tjet"]].n_observed_max = esi.pda[esi.observed_object["atop"]].n_observed_max;}
    else {esi.pda[esi.observed_object["atop"]].n_observed_max = esi.pda[esi.observed_object["tjet"]].n_observed_max;}
  }
  else if (esi.pda[esi.observed_object["atop"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["top"]].n_observed_max <= esi.pda[esi.observed_object["tjet"]].n_observed_max){esi.pda[esi.observed_object["tjet"]].n_observed_max = esi.pda[esi.observed_object["top"]].n_observed_max;}
    else {esi.pda[esi.observed_object["top"]].n_observed_max = esi.pda[esi.observed_object["tjet"]].n_observed_max;}
  }

  // adapt n_observed_min to parton-level particles (tjet-top/atop)
  if (esi.pda[esi.observed_object["tjet"]].n_partonlevel == esi.pda[esi.observed_object["tjet"]].n_observed_min){
    esi.pda[esi.observed_object["atop"]].n_observed_min = esi.pda[esi.observed_object["atop"]].n_partonlevel;
    esi.pda[esi.observed_object["top"]].n_observed_min = esi.pda[esi.observed_object["top"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["top"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["atop"]].n_observed_min >= esi.pda[esi.observed_object["tjet"]].n_observed_min){esi.pda[esi.observed_object["tjet"]].n_observed_min = esi.pda[esi.observed_object["atop"]].n_observed_min;}
    else {esi.pda[esi.observed_object["atop"]].n_observed_min = esi.pda[esi.observed_object["tjet"]].n_observed_min;}
  }
  else if (esi.pda[esi.observed_object["atop"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["top"]].n_observed_min >= esi.pda[esi.observed_object["tjet"]].n_observed_min){esi.pda[esi.observed_object["tjet"]].n_observed_min = esi.pda[esi.observed_object["top"]].n_observed_min;}
    else {esi.pda[esi.observed_object["top"]].n_observed_min = esi.pda[esi.observed_object["tjet"]].n_observed_min;}
  }


  // adapt n_observed_max to parton-level particles (bjet-bjet_b/bjet_bx)
  if (esi.pda[esi.observed_object["bjet_bx"]].n_partonlevel < esi.pda[esi.observed_object["bjet_bx"]].n_observed_max){esi.pda[esi.observed_object["bjet_bx"]].n_observed_max = esi.pda[esi.observed_object["bjet_bx"]].n_partonlevel;}
  if (esi.pda[esi.observed_object["bjet_b"]].n_partonlevel < esi.pda[esi.observed_object["bjet_b"]].n_observed_max){esi.pda[esi.observed_object["bjet_b"]].n_observed_max = esi.pda[esi.observed_object["bjet_b"]].n_partonlevel;}

  if (esi.pda[esi.observed_object["bjet"]].n_partonlevel == esi.pda[esi.observed_object["bjet"]].n_observed_max){
    esi.pda[esi.observed_object["bjet_bx"]].n_observed_max = esi.pda[esi.observed_object["bjet_bx"]].n_partonlevel;
    esi.pda[esi.observed_object["bjet_b"]].n_observed_max = esi.pda[esi.observed_object["bjet_b"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["bjet_b"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["bjet_bx"]].n_observed_max <= esi.pda[esi.observed_object["bjet"]].n_observed_max){esi.pda[esi.observed_object["bjet"]].n_observed_max = esi.pda[esi.observed_object["bjet_bx"]].n_observed_max;}
    else {esi.pda[esi.observed_object["bjet_bx"]].n_observed_max = esi.pda[esi.observed_object["bjet"]].n_observed_max;}
  }
  else if (esi.pda[esi.observed_object["bjet_bx"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["bjet_b"]].n_observed_max <= esi.pda[esi.observed_object["bjet"]].n_observed_max){esi.pda[esi.observed_object["bjet"]].n_observed_max = esi.pda[esi.observed_object["bjet_b"]].n_observed_max;}
    else {esi.pda[esi.observed_object["bjet_b"]].n_observed_max = esi.pda[esi.observed_object["bjet"]].n_observed_max;}
  }

  // adapt n_observed_min to parton-level particles (bjet-bjet_b/bjet_bx)
  if (esi.pda[esi.observed_object["bjet"]].n_partonlevel == esi.pda[esi.observed_object["bjet"]].n_observed_min){
    esi.pda[esi.observed_object["bjet_bx"]].n_observed_min = esi.pda[esi.observed_object["bjet_bx"]].n_partonlevel;
    esi.pda[esi.observed_object["bjet_b"]].n_observed_min = esi.pda[esi.observed_object["bjet_b"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["bjet_b"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["bjet_bx"]].n_observed_min >= esi.pda[esi.observed_object["bjet"]].n_observed_min){esi.pda[esi.observed_object["bjet"]].n_observed_min = esi.pda[esi.observed_object["bjet_bx"]].n_observed_min;}
    else {esi.pda[esi.observed_object["bjet_bx"]].n_observed_min = esi.pda[esi.observed_object["bjet"]].n_observed_min;}
  }
  else if (esi.pda[esi.observed_object["bjet_bx"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["bjet_b"]].n_observed_min >= esi.pda[esi.observed_object["bjet"]].n_observed_min){esi.pda[esi.observed_object["bjet"]].n_observed_min = esi.pda[esi.observed_object["bjet_b"]].n_observed_min;}
    else {esi.pda[esi.observed_object["bjet_b"]].n_observed_min = esi.pda[esi.observed_object["bjet"]].n_observed_min;}
  }


  // adapt n_observed_max to parton-level particles (lep-lp/lm) !
  if (esi.pda[esi.observed_object["lep"]].n_partonlevel < esi.pda[esi.observed_object["lep"]].n_observed_max){esi.pda[esi.observed_object["lep"]].n_observed_max = esi.pda[esi.observed_object["lep"]].n_partonlevel;}

  if (esi.pda[esi.observed_object["lep"]].n_partonlevel == esi.pda[esi.observed_object["lep"]].n_observed_max){
    esi.pda[esi.observed_object["lm"]].n_observed_max = esi.pda[esi.observed_object["lm"]].n_partonlevel;
    esi.pda[esi.observed_object["lp"]].n_observed_max = esi.pda[esi.observed_object["lp"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["lp"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["lm"]].n_observed_max <= esi.pda[esi.observed_object["lep"]].n_observed_max){esi.pda[esi.observed_object["lep"]].n_observed_max = esi.pda[esi.observed_object["lm"]].n_observed_max;}
    else {esi.pda[esi.observed_object["lm"]].n_observed_max = esi.pda[esi.observed_object["lep"]].n_observed_max;}
  }
  else if (esi.pda[esi.observed_object["lm"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["lp"]].n_observed_max <= esi.pda[esi.observed_object["lep"]].n_observed_max){esi.pda[esi.observed_object["lep"]].n_observed_max = esi.pda[esi.observed_object["lp"]].n_observed_max;}
    else {esi.pda[esi.observed_object["lp"]].n_observed_max = esi.pda[esi.observed_object["lep"]].n_observed_max;}
  }

  // adapt n_observed_min to parton-level particles (lep-lp/lm) !
  if (esi.pda[esi.observed_object["lep"]].n_partonlevel == esi.pda[esi.observed_object["lep"]].n_observed_min){
    if (esi.pda[esi.observed_object["lm"]].n_observed_min < esi.pda[esi.observed_object["lm"]].n_partonlevel){esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["lm"]].n_partonlevel;}
    if (esi.pda[esi.observed_object["lp"]].n_observed_min < esi.pda[esi.observed_object["lp"]].n_partonlevel){esi.pda[esi.observed_object["lp"]].n_observed_min = esi.pda[esi.observed_object["lp"]].n_partonlevel;}
  }
  /*
  else if (esi.pda[esi.observed_object["lp"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["lm"]].n_observed_min >= esi.pda[esi.observed_object["lep"]].n_observed_min){esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["lm"]].n_observed_min;}
    else {esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["lep"]].n_observed_min;}
  }
  else if (esi.pda[esi.observed_object["lm"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["lp"]].n_observed_min >= esi.pda[esi.observed_object["lep"]].n_observed_min){esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["lp"]].n_observed_min;}
    else {esi.pda[esi.observed_object["lp"]].n_observed_min = esi.pda[esi.observed_object["lep"]].n_observed_min;}
  }
  */
  /*
  // adapt n_observed_min to parton-level particles (lep-lp/lm) !
  if (esi.pda[esi.observed_object["lep"]].n_partonlevel == esi.pda[esi.observed_object["lep"]].n_observed_min){
    esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["lm"]].n_partonlevel;
    esi.pda[esi.observed_object["lp"]].n_observed_min = esi.pda[esi.observed_object["lp"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["lp"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["lm"]].n_observed_min >= esi.pda[esi.observed_object["lep"]].n_observed_min){esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["lm"]].n_observed_min;}
    else {esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["lep"]].n_observed_min;}
  }
  else if (esi.pda[esi.observed_object["lm"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["lp"]].n_observed_min >= esi.pda[esi.observed_object["lep"]].n_observed_min){esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["lp"]].n_observed_min;}
    else {esi.pda[esi.observed_object["lp"]].n_observed_min = esi.pda[esi.observed_object["lep"]].n_observed_min;}
  }
  */

  // adapt n_observed_min to parton-level particles (lep-e/mu/tau) !
  if (esi.pda[esi.observed_object["lep"]].n_partonlevel == esi.pda[esi.observed_object["lep"]].n_observed_min){
    if (esi.pda[esi.observed_object["e"]].n_observed_min < esi.pda[esi.observed_object["e"]].n_partonlevel){esi.pda[esi.observed_object["e"]].n_observed_min = esi.pda[esi.observed_object["e"]].n_partonlevel;}
    if (esi.pda[esi.observed_object["mu"]].n_observed_min < esi.pda[esi.observed_object["mu"]].n_partonlevel){esi.pda[esi.observed_object["mu"]].n_observed_min = esi.pda[esi.observed_object["mu"]].n_partonlevel;}
    if (esi.pda[esi.observed_object["tau"]].n_observed_min < esi.pda[esi.observed_object["tau"]].n_partonlevel){esi.pda[esi.observed_object["tau"]].n_observed_min = esi.pda[esi.observed_object["tau"]].n_partonlevel;}
    /*
    // modify: This could prevent one from switching off contributions that have less partons off that flavour than individually required by input
    esi.pda[esi.observed_object["e"]].n_observed_min = esi.pda[esi.observed_object["e"]].n_partonlevel;
    esi.pda[esi.observed_object["mu"]].n_observed_min = esi.pda[esi.observed_object["mu"]].n_partonlevel;
    esi.pda[esi.observed_object["tau"]].n_observed_min = esi.pda[esi.observed_object["tau"]].n_partonlevel;
    */
  }
  else if ((esi.pda[esi.observed_object["e"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mu"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["tau"]].n_observed_min >= esi.pda[esi.observed_object["lep"]].n_observed_min){esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["tau"]].n_observed_min;}
    else {esi.pda[esi.observed_object["tau"]].n_observed_min = esi.pda[esi.observed_object["lep"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["e"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["tau"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mu"]].n_observed_min > esi.pda[esi.observed_object["lep"]].n_observed_min)){
    if (esi.pda[esi.observed_object["mu"]].n_observed_min >= esi.pda[esi.observed_object["lep"]].n_observed_min){esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["mu"]].n_observed_min;}
    else {esi.pda[esi.observed_object["mu"]].n_observed_min = esi.pda[esi.observed_object["lep"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["mu"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["tau"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["e"]].n_observed_min > esi.pda[esi.observed_object["lep"]].n_observed_min)){
    if (esi.pda[esi.observed_object["e"]].n_observed_min >= esi.pda[esi.observed_object["lep"]].n_observed_min){esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["e"]].n_observed_min;}
    else {esi.pda[esi.observed_object["e"]].n_observed_min = esi.pda[esi.observed_object["lep"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["e"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mu"]].n_observed_min + esi.pda[esi.observed_object["tau"]].n_observed_min > esi.pda[esi.observed_object["lep"]].n_observed_min)){
    esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["mu"]].n_observed_min + esi.pda[esi.observed_object["tau"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["mu"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["e"]].n_observed_min + esi.pda[esi.observed_object["tau"]].n_observed_min > esi.pda[esi.observed_object["lep"]].n_observed_min)){
    esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["e"]].n_observed_min + esi.pda[esi.observed_object["tau"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["tau"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["e"]].n_observed_min + esi.pda[esi.observed_object["mu"]].n_observed_min > esi.pda[esi.observed_object["lep"]].n_observed_min)){
    esi.pda[esi.observed_object["lep"]].n_observed_min = esi.pda[esi.observed_object["e"]].n_observed_min + esi.pda[esi.observed_object["mu"]].n_observed_min;
  }

  // adapt n_observed_max to parton-level particles (lep-e/mu/tau) !
  if (esi.pda[esi.observed_object["lep"]].n_partonlevel == esi.pda[esi.observed_object["lep"]].n_observed_max){
    esi.pda[esi.observed_object["e"]].n_observed_max = esi.pda[esi.observed_object["e"]].n_partonlevel;
    esi.pda[esi.observed_object["mu"]].n_observed_max = esi.pda[esi.observed_object["mu"]].n_partonlevel;
    esi.pda[esi.observed_object["tau"]].n_observed_max = esi.pda[esi.observed_object["tau"]].n_partonlevel;
  }
  else if ((esi.pda[esi.observed_object["e"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mu"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["tau"]].n_observed_max <= esi.pda[esi.observed_object["lep"]].n_observed_max){esi.pda[esi.observed_object["lep"]].n_observed_max = esi.pda[esi.observed_object["tau"]].n_observed_max;}
    else {esi.pda[esi.observed_object["tau"]].n_observed_max = esi.pda[esi.observed_object["lep"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["e"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["tau"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mu"]].n_observed_max > esi.pda[esi.observed_object["lep"]].n_observed_max)){
    if (esi.pda[esi.observed_object["mu"]].n_observed_max <= esi.pda[esi.observed_object["lep"]].n_observed_max){esi.pda[esi.observed_object["lep"]].n_observed_max = esi.pda[esi.observed_object["mu"]].n_observed_max;}
    else {esi.pda[esi.observed_object["mu"]].n_observed_max = esi.pda[esi.observed_object["lep"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["mu"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["tau"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["e"]].n_observed_max > esi.pda[esi.observed_object["lep"]].n_observed_max)){
    if (esi.pda[esi.observed_object["e"]].n_observed_max <= esi.pda[esi.observed_object["lep"]].n_observed_max){esi.pda[esi.observed_object["lep"]].n_observed_max = esi.pda[esi.observed_object["e"]].n_observed_max;}
    else {esi.pda[esi.observed_object["e"]].n_observed_max = esi.pda[esi.observed_object["lep"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["e"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mu"]].n_observed_max + esi.pda[esi.observed_object["tau"]].n_observed_max < esi.pda[esi.observed_object["lep"]].n_observed_max)){
    esi.pda[esi.observed_object["lep"]].n_observed_max = esi.pda[esi.observed_object["mu"]].n_observed_max + esi.pda[esi.observed_object["tau"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["mu"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["e"]].n_observed_max + esi.pda[esi.observed_object["tau"]].n_observed_max < esi.pda[esi.observed_object["lep"]].n_observed_max)){
    esi.pda[esi.observed_object["lep"]].n_observed_max = esi.pda[esi.observed_object["e"]].n_observed_max + esi.pda[esi.observed_object["tau"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["tau"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["e"]].n_observed_max + esi.pda[esi.observed_object["mu"]].n_observed_max < esi.pda[esi.observed_object["lep"]].n_observed_max)){
    esi.pda[esi.observed_object["lep"]].n_observed_max = esi.pda[esi.observed_object["e"]].n_observed_max + esi.pda[esi.observed_object["mu"]].n_observed_max;
  }

  // adapt n_observed_min to parton-level particles (lm-em/mum/taum) !
  if (esi.pda[esi.observed_object["lm"]].n_partonlevel == esi.pda[esi.observed_object["lm"]].n_observed_min){
    esi.pda[esi.observed_object["em"]].n_observed_min = esi.pda[esi.observed_object["em"]].n_partonlevel;
    esi.pda[esi.observed_object["mum"]].n_observed_min = esi.pda[esi.observed_object["mum"]].n_partonlevel;
    esi.pda[esi.observed_object["taum"]].n_observed_min = esi.pda[esi.observed_object["taum"]].n_partonlevel;
  }
  else if ((esi.pda[esi.observed_object["em"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mum"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["taum"]].n_observed_min >= esi.pda[esi.observed_object["lm"]].n_observed_min){esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["taum"]].n_observed_min;}
    else {esi.pda[esi.observed_object["taum"]].n_observed_min = esi.pda[esi.observed_object["lm"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["em"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["taum"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mum"]].n_observed_min > esi.pda[esi.observed_object["lm"]].n_observed_min)){
    if (esi.pda[esi.observed_object["mum"]].n_observed_min >= esi.pda[esi.observed_object["lm"]].n_observed_min){esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["mum"]].n_observed_min;}
    else {esi.pda[esi.observed_object["mum"]].n_observed_min = esi.pda[esi.observed_object["lm"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["mum"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["taum"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["em"]].n_observed_min > esi.pda[esi.observed_object["lm"]].n_observed_min)){
    if (esi.pda[esi.observed_object["em"]].n_observed_min >= esi.pda[esi.observed_object["lm"]].n_observed_min){esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["em"]].n_observed_min;}
    else {esi.pda[esi.observed_object["em"]].n_observed_min = esi.pda[esi.observed_object["lm"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["em"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mum"]].n_observed_min + esi.pda[esi.observed_object["taum"]].n_observed_min > esi.pda[esi.observed_object["lm"]].n_observed_min)){
    esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["mum"]].n_observed_min + esi.pda[esi.observed_object["taum"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["mum"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["em"]].n_observed_min + esi.pda[esi.observed_object["taum"]].n_observed_min > esi.pda[esi.observed_object["lm"]].n_observed_min)){
    esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["em"]].n_observed_min + esi.pda[esi.observed_object["taum"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["taum"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["em"]].n_observed_min + esi.pda[esi.observed_object["mum"]].n_observed_min > esi.pda[esi.observed_object["lm"]].n_observed_min)){
    esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["em"]].n_observed_min + esi.pda[esi.observed_object["mum"]].n_observed_min;
  }

  // adapt n_observed_max to parton-level particles (lm-em/mum/taum) !
  if (esi.pda[esi.observed_object["lm"]].n_partonlevel == esi.pda[esi.observed_object["lm"]].n_observed_max){
    esi.pda[esi.observed_object["em"]].n_observed_max = esi.pda[esi.observed_object["em"]].n_partonlevel;
    esi.pda[esi.observed_object["mum"]].n_observed_max = esi.pda[esi.observed_object["mum"]].n_partonlevel;
    esi.pda[esi.observed_object["taum"]].n_observed_max = esi.pda[esi.observed_object["taum"]].n_partonlevel;
  }
  else if ((esi.pda[esi.observed_object["em"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mum"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["taum"]].n_observed_max <= esi.pda[esi.observed_object["lm"]].n_observed_max){esi.pda[esi.observed_object["lm"]].n_observed_max = esi.pda[esi.observed_object["taum"]].n_observed_max;}
    else {esi.pda[esi.observed_object["taum"]].n_observed_max = esi.pda[esi.observed_object["lm"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["em"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["taum"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mum"]].n_observed_max > esi.pda[esi.observed_object["lm"]].n_observed_max)){
    if (esi.pda[esi.observed_object["mum"]].n_observed_max <= esi.pda[esi.observed_object["lm"]].n_observed_max){esi.pda[esi.observed_object["lm"]].n_observed_max = esi.pda[esi.observed_object["mum"]].n_observed_max;}
    else {esi.pda[esi.observed_object["mum"]].n_observed_max = esi.pda[esi.observed_object["lm"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["mum"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["taum"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["em"]].n_observed_max > esi.pda[esi.observed_object["lm"]].n_observed_max)){
    if (esi.pda[esi.observed_object["em"]].n_observed_max <= esi.pda[esi.observed_object["lm"]].n_observed_max){esi.pda[esi.observed_object["lm"]].n_observed_max = esi.pda[esi.observed_object["em"]].n_observed_max;}
    else {esi.pda[esi.observed_object["em"]].n_observed_max = esi.pda[esi.observed_object["lm"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["em"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mum"]].n_observed_max + esi.pda[esi.observed_object["taum"]].n_observed_max < esi.pda[esi.observed_object["lm"]].n_observed_max)){
    esi.pda[esi.observed_object["lm"]].n_observed_max = esi.pda[esi.observed_object["mum"]].n_observed_max + esi.pda[esi.observed_object["taum"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["mum"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["em"]].n_observed_max + esi.pda[esi.observed_object["taum"]].n_observed_max < esi.pda[esi.observed_object["lm"]].n_observed_max)){
    esi.pda[esi.observed_object["lm"]].n_observed_max = esi.pda[esi.observed_object["em"]].n_observed_max + esi.pda[esi.observed_object["taum"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["taum"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["em"]].n_observed_max + esi.pda[esi.observed_object["mum"]].n_observed_max < esi.pda[esi.observed_object["lm"]].n_observed_max)){
    esi.pda[esi.observed_object["lm"]].n_observed_max = esi.pda[esi.observed_object["em"]].n_observed_max + esi.pda[esi.observed_object["mum"]].n_observed_max;
  }

  // adapt n_observed_min to parton-level particles (lp-ep/mup/taup) !
  if (esi.pda[esi.observed_object["lp"]].n_partonlevel == esi.pda[esi.observed_object["lp"]].n_observed_min){
    esi.pda[esi.observed_object["ep"]].n_observed_min = esi.pda[esi.observed_object["ep"]].n_partonlevel;
    esi.pda[esi.observed_object["mup"]].n_observed_min = esi.pda[esi.observed_object["mup"]].n_partonlevel;
    esi.pda[esi.observed_object["taup"]].n_observed_min = esi.pda[esi.observed_object["taup"]].n_partonlevel;
  }
  else if ((esi.pda[esi.observed_object["ep"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mup"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["taup"]].n_observed_min >= esi.pda[esi.observed_object["lm"]].n_observed_min){esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["taup"]].n_observed_min;}
    else {esi.pda[esi.observed_object["taup"]].n_observed_min = esi.pda[esi.observed_object["lm"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["ep"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["taup"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mup"]].n_observed_min > esi.pda[esi.observed_object["lm"]].n_observed_min)){
    if (esi.pda[esi.observed_object["mup"]].n_observed_min >= esi.pda[esi.observed_object["lm"]].n_observed_min){esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["mup"]].n_observed_min;}
    else {esi.pda[esi.observed_object["mup"]].n_observed_min = esi.pda[esi.observed_object["lm"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["mup"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["taup"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ep"]].n_observed_min > esi.pda[esi.observed_object["lm"]].n_observed_min)){
    if (esi.pda[esi.observed_object["ep"]].n_observed_min >= esi.pda[esi.observed_object["lm"]].n_observed_min){esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["ep"]].n_observed_min;}
    else {esi.pda[esi.observed_object["ep"]].n_observed_min = esi.pda[esi.observed_object["lm"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["ep"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mup"]].n_observed_min + esi.pda[esi.observed_object["taup"]].n_observed_min > esi.pda[esi.observed_object["lm"]].n_observed_min)){
    esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["mup"]].n_observed_min + esi.pda[esi.observed_object["taup"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["mup"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ep"]].n_observed_min + esi.pda[esi.observed_object["taup"]].n_observed_min > esi.pda[esi.observed_object["lm"]].n_observed_min)){
    esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["ep"]].n_observed_min + esi.pda[esi.observed_object["taup"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["taup"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ep"]].n_observed_min + esi.pda[esi.observed_object["mup"]].n_observed_min > esi.pda[esi.observed_object["lm"]].n_observed_min)){
    esi.pda[esi.observed_object["lm"]].n_observed_min = esi.pda[esi.observed_object["ep"]].n_observed_min + esi.pda[esi.observed_object["mup"]].n_observed_min;
  }

  // adapt n_observed_max to parton-level particles (lp-ep/mup/taup) !
  if (esi.pda[esi.observed_object["lp"]].n_partonlevel == esi.pda[esi.observed_object["lp"]].n_observed_max){
    esi.pda[esi.observed_object["ep"]].n_observed_max = esi.pda[esi.observed_object["ep"]].n_partonlevel;
    esi.pda[esi.observed_object["mup"]].n_observed_max = esi.pda[esi.observed_object["mup"]].n_partonlevel;
    esi.pda[esi.observed_object["taup"]].n_observed_max = esi.pda[esi.observed_object["taup"]].n_partonlevel;
  }
  else if ((esi.pda[esi.observed_object["ep"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mup"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["taup"]].n_observed_max <= esi.pda[esi.observed_object["lp"]].n_observed_max){esi.pda[esi.observed_object["lp"]].n_observed_max = esi.pda[esi.observed_object["taup"]].n_observed_max;}
    else {esi.pda[esi.observed_object["taup"]].n_observed_max = esi.pda[esi.observed_object["lp"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["ep"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["taup"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mup"]].n_observed_max > esi.pda[esi.observed_object["lp"]].n_observed_max)){
    if (esi.pda[esi.observed_object["mup"]].n_observed_max <= esi.pda[esi.observed_object["lp"]].n_observed_max){esi.pda[esi.observed_object["lp"]].n_observed_max = esi.pda[esi.observed_object["mup"]].n_observed_max;}
    else {esi.pda[esi.observed_object["mup"]].n_observed_max = esi.pda[esi.observed_object["lp"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["mup"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["taup"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ep"]].n_observed_max > esi.pda[esi.observed_object["lp"]].n_observed_max)){
    if (esi.pda[esi.observed_object["ep"]].n_observed_max <= esi.pda[esi.observed_object["lp"]].n_observed_max){esi.pda[esi.observed_object["lp"]].n_observed_max = esi.pda[esi.observed_object["ep"]].n_observed_max;}
    else {esi.pda[esi.observed_object["ep"]].n_observed_max = esi.pda[esi.observed_object["lp"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["ep"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["mup"]].n_observed_max + esi.pda[esi.observed_object["taup"]].n_observed_max < esi.pda[esi.observed_object["lp"]].n_observed_max)){
    esi.pda[esi.observed_object["lp"]].n_observed_max = esi.pda[esi.observed_object["mup"]].n_observed_max + esi.pda[esi.observed_object["taup"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["mup"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ep"]].n_observed_max + esi.pda[esi.observed_object["taup"]].n_observed_max < esi.pda[esi.observed_object["lp"]].n_observed_max)){
    esi.pda[esi.observed_object["lp"]].n_observed_max = esi.pda[esi.observed_object["ep"]].n_observed_max + esi.pda[esi.observed_object["taup"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["taup"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ep"]].n_observed_max + esi.pda[esi.observed_object["mup"]].n_observed_max < esi.pda[esi.observed_object["lp"]].n_observed_max)){
    esi.pda[esi.observed_object["lp"]].n_observed_max = esi.pda[esi.observed_object["ep"]].n_observed_max + esi.pda[esi.observed_object["mup"]].n_observed_max;
  }

  // adapt n_observed_max to parton-level particles (e-ep/em) !
  if (esi.pda[esi.observed_object["e"]].n_partonlevel == esi.pda[esi.observed_object["e"]].n_observed_max){
    esi.pda[esi.observed_object["em"]].n_observed_max = esi.pda[esi.observed_object["em"]].n_partonlevel;
    esi.pda[esi.observed_object["ep"]].n_observed_max = esi.pda[esi.observed_object["ep"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["ep"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["em"]].n_observed_max <= esi.pda[esi.observed_object["e"]].n_observed_max){esi.pda[esi.observed_object["e"]].n_observed_max = esi.pda[esi.observed_object["em"]].n_observed_max;}
    else {esi.pda[esi.observed_object["em"]].n_observed_max = esi.pda[esi.observed_object["e"]].n_observed_max;}
  }
  else if (esi.pda[esi.observed_object["em"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["ep"]].n_observed_max <= esi.pda[esi.observed_object["e"]].n_observed_max){esi.pda[esi.observed_object["e"]].n_observed_max = esi.pda[esi.observed_object["ep"]].n_observed_max;}
    else {esi.pda[esi.observed_object["ep"]].n_observed_max = esi.pda[esi.observed_object["e"]].n_observed_max;}
  }

  // adapt n_observed_min to parton-level particles (e-ep/em) !
  if (esi.pda[esi.observed_object["e"]].n_partonlevel == esi.pda[esi.observed_object["e"]].n_observed_min){
    esi.pda[esi.observed_object["em"]].n_observed_min = esi.pda[esi.observed_object["em"]].n_partonlevel;
    esi.pda[esi.observed_object["ep"]].n_observed_min = esi.pda[esi.observed_object["ep"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["ep"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["em"]].n_observed_min >= esi.pda[esi.observed_object["e"]].n_observed_min){esi.pda[esi.observed_object["e"]].n_observed_min = esi.pda[esi.observed_object["em"]].n_observed_min;}
    else {esi.pda[esi.observed_object["em"]].n_observed_min = esi.pda[esi.observed_object["e"]].n_observed_min;}
  }
  else if (esi.pda[esi.observed_object["em"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["ep"]].n_observed_min >= esi.pda[esi.observed_object["e"]].n_observed_min){esi.pda[esi.observed_object["e"]].n_observed_min = esi.pda[esi.observed_object["ep"]].n_observed_min;}
    else {esi.pda[esi.observed_object["ep"]].n_observed_min = esi.pda[esi.observed_object["e"]].n_observed_min;}
  }

  // adapt n_observed_max to parton-level particles (mu-mup/mum) !
  if (esi.pda[esi.observed_object["mu"]].n_partonlevel == esi.pda[esi.observed_object["mu"]].n_observed_max){
    esi.pda[esi.observed_object["mum"]].n_observed_max = esi.pda[esi.observed_object["mum"]].n_partonlevel;
    esi.pda[esi.observed_object["mup"]].n_observed_max = esi.pda[esi.observed_object["mup"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["mup"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["mum"]].n_observed_max <= esi.pda[esi.observed_object["mu"]].n_observed_max){esi.pda[esi.observed_object["mu"]].n_observed_max = esi.pda[esi.observed_object["mum"]].n_observed_max;}
    else {esi.pda[esi.observed_object["mum"]].n_observed_max = esi.pda[esi.observed_object["mu"]].n_observed_max;}
  }
  else if (esi.pda[esi.observed_object["mum"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["mup"]].n_observed_max <= esi.pda[esi.observed_object["mu"]].n_observed_max){esi.pda[esi.observed_object["mu"]].n_observed_max = esi.pda[esi.observed_object["mup"]].n_observed_max;}
    else {esi.pda[esi.observed_object["mup"]].n_observed_max = esi.pda[esi.observed_object["mu"]].n_observed_max;}
  }

  // adapt n_observed_min to parton-level particles (mu-mup/mum) !
  if (esi.pda[esi.observed_object["mu"]].n_partonlevel == esi.pda[esi.observed_object["mu"]].n_observed_min){
    esi.pda[esi.observed_object["mum"]].n_observed_min = esi.pda[esi.observed_object["mum"]].n_partonlevel;
    esi.pda[esi.observed_object["mup"]].n_observed_min = esi.pda[esi.observed_object["mup"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["mup"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["mum"]].n_observed_min >= esi.pda[esi.observed_object["mu"]].n_observed_min){esi.pda[esi.observed_object["mu"]].n_observed_min = esi.pda[esi.observed_object["mum"]].n_observed_min;}
    else {esi.pda[esi.observed_object["mum"]].n_observed_min = esi.pda[esi.observed_object["mu"]].n_observed_min;}
  }
  else if (esi.pda[esi.observed_object["mum"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["mup"]].n_observed_min >= esi.pda[esi.observed_object["mu"]].n_observed_min){esi.pda[esi.observed_object["mu"]].n_observed_min = esi.pda[esi.observed_object["mup"]].n_observed_min;}
    else {esi.pda[esi.observed_object["mup"]].n_observed_min = esi.pda[esi.observed_object["mu"]].n_observed_min;}
  }

  // adapt n_observed_max to parton-level particles (tau-taup/taum) !
  if (esi.pda[esi.observed_object["tau"]].n_partonlevel == esi.pda[esi.observed_object["tau"]].n_observed_max){
    esi.pda[esi.observed_object["taum"]].n_observed_max = esi.pda[esi.observed_object["taum"]].n_partonlevel;
    esi.pda[esi.observed_object["taup"]].n_observed_max = esi.pda[esi.observed_object["taup"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["taup"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["taum"]].n_observed_max <= esi.pda[esi.observed_object["tau"]].n_observed_max){esi.pda[esi.observed_object["tau"]].n_observed_max = esi.pda[esi.observed_object["taum"]].n_observed_max;}
    else {esi.pda[esi.observed_object["taum"]].n_observed_max = esi.pda[esi.observed_object["tau"]].n_observed_max;}
  }
  else if (esi.pda[esi.observed_object["taum"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["taup"]].n_observed_max <= esi.pda[esi.observed_object["tau"]].n_observed_max){esi.pda[esi.observed_object["tau"]].n_observed_max = esi.pda[esi.observed_object["taup"]].n_observed_max;}
    else {esi.pda[esi.observed_object["taup"]].n_observed_max = esi.pda[esi.observed_object["tau"]].n_observed_max;}
  }

  // adapt n_observed_min to parton-level particles (tau-taup/taum) !
  if (esi.pda[esi.observed_object["tau"]].n_partonlevel == esi.pda[esi.observed_object["tau"]].n_observed_min){
    esi.pda[esi.observed_object["taum"]].n_observed_min = esi.pda[esi.observed_object["taum"]].n_partonlevel;
    esi.pda[esi.observed_object["taup"]].n_observed_min = esi.pda[esi.observed_object["taup"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["taup"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["taum"]].n_observed_min >= esi.pda[esi.observed_object["tau"]].n_observed_min){esi.pda[esi.observed_object["tau"]].n_observed_min = esi.pda[esi.observed_object["taum"]].n_observed_min;}
    else {esi.pda[esi.observed_object["taum"]].n_observed_min = esi.pda[esi.observed_object["tau"]].n_observed_min;}
  }
  else if (esi.pda[esi.observed_object["taum"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["taup"]].n_observed_min >= esi.pda[esi.observed_object["tau"]].n_observed_min){esi.pda[esi.observed_object["tau"]].n_observed_min = esi.pda[esi.observed_object["taup"]].n_observed_min;}
    else {esi.pda[esi.observed_object["taup"]].n_observed_min = esi.pda[esi.observed_object["tau"]].n_observed_min;}
  }

  // adapt n_observed_max to parton-level particles (nua-nux/nu) !
  if (esi.pda[esi.observed_object["nua"]].n_partonlevel < esi.pda[esi.observed_object["nua"]].n_observed_max){esi.pda[esi.observed_object["nua"]].n_observed_max = esi.pda[esi.observed_object["nua"]].n_partonlevel;}

  if (esi.pda[esi.observed_object["nua"]].n_partonlevel == esi.pda[esi.observed_object["nua"]].n_observed_max){
    esi.pda[esi.observed_object["nu"]].n_observed_max = esi.pda[esi.observed_object["nu"]].n_partonlevel;
    esi.pda[esi.observed_object["nux"]].n_observed_max = esi.pda[esi.observed_object["nux"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["nux"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nu"]].n_observed_max <= esi.pda[esi.observed_object["nua"]].n_observed_max){esi.pda[esi.observed_object["nua"]].n_observed_max = esi.pda[esi.observed_object["nu"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nu"]].n_observed_max = esi.pda[esi.observed_object["nua"]].n_observed_max;}
  }
  else if (esi.pda[esi.observed_object["nu"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nux"]].n_observed_max <= esi.pda[esi.observed_object["nua"]].n_observed_max){esi.pda[esi.observed_object["nua"]].n_observed_max = esi.pda[esi.observed_object["nux"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nux"]].n_observed_max = esi.pda[esi.observed_object["nua"]].n_observed_max;}
  }

  // adapt n_observed_min to parton-level particles (nua-nux/nu) !
  if (esi.pda[esi.observed_object["nua"]].n_partonlevel == esi.pda[esi.observed_object["nua"]].n_observed_min){
    esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["nu"]].n_partonlevel;
    esi.pda[esi.observed_object["nux"]].n_observed_min = esi.pda[esi.observed_object["nux"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["nux"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nu"]].n_observed_min >= esi.pda[esi.observed_object["nua"]].n_observed_min){esi.pda[esi.observed_object["nua"]].n_observed_min = esi.pda[esi.observed_object["nu"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["nua"]].n_observed_min;}
  }
  else if (esi.pda[esi.observed_object["nu"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nux"]].n_observed_min >= esi.pda[esi.observed_object["nua"]].n_observed_min){esi.pda[esi.observed_object["nua"]].n_observed_min = esi.pda[esi.observed_object["nux"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nux"]].n_observed_min = esi.pda[esi.observed_object["nua"]].n_observed_min;}
  }

  // adapt n_observed_min to parton-level particles (nua-nea/nma/nta) !
  if (esi.pda[esi.observed_object["nua"]].n_partonlevel == esi.pda[esi.observed_object["nua"]].n_observed_min){
    esi.pda[esi.observed_object["nea"]].n_observed_min = esi.pda[esi.observed_object["nea"]].n_partonlevel;
    esi.pda[esi.observed_object["nma"]].n_observed_min = esi.pda[esi.observed_object["nma"]].n_partonlevel;
    esi.pda[esi.observed_object["nta"]].n_observed_min = esi.pda[esi.observed_object["nta"]].n_partonlevel;
  }
  else if ((esi.pda[esi.observed_object["nea"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nma"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["nta"]].n_observed_min >= esi.pda[esi.observed_object["nua"]].n_observed_min){esi.pda[esi.observed_object["nua"]].n_observed_min = esi.pda[esi.observed_object["nta"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nta"]].n_observed_min = esi.pda[esi.observed_object["nua"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["nea"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nta"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nma"]].n_observed_min > esi.pda[esi.observed_object["nua"]].n_observed_min)){
    if (esi.pda[esi.observed_object["nma"]].n_observed_min >= esi.pda[esi.observed_object["nua"]].n_observed_min){esi.pda[esi.observed_object["nua"]].n_observed_min = esi.pda[esi.observed_object["nma"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nma"]].n_observed_min = esi.pda[esi.observed_object["nua"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["nma"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nta"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nea"]].n_observed_min > esi.pda[esi.observed_object["nua"]].n_observed_min)){
    if (esi.pda[esi.observed_object["nea"]].n_observed_min >= esi.pda[esi.observed_object["nua"]].n_observed_min){esi.pda[esi.observed_object["nua"]].n_observed_min = esi.pda[esi.observed_object["nea"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nea"]].n_observed_min = esi.pda[esi.observed_object["nua"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["nea"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nma"]].n_observed_min + esi.pda[esi.observed_object["nta"]].n_observed_min > esi.pda[esi.observed_object["nua"]].n_observed_min)){
    esi.pda[esi.observed_object["nua"]].n_observed_min = esi.pda[esi.observed_object["nma"]].n_observed_min + esi.pda[esi.observed_object["nta"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["nma"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nea"]].n_observed_min + esi.pda[esi.observed_object["nta"]].n_observed_min > esi.pda[esi.observed_object["nua"]].n_observed_min)){
    esi.pda[esi.observed_object["nua"]].n_observed_min = esi.pda[esi.observed_object["nea"]].n_observed_min + esi.pda[esi.observed_object["nta"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["nta"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nea"]].n_observed_min + esi.pda[esi.observed_object["nma"]].n_observed_min > esi.pda[esi.observed_object["nua"]].n_observed_min)){
    esi.pda[esi.observed_object["nua"]].n_observed_min = esi.pda[esi.observed_object["nea"]].n_observed_min + esi.pda[esi.observed_object["nma"]].n_observed_min;
  }

  // adapt n_observed_max to parton-level particles (nua-nea/nma/nta) !
  if (esi.pda[esi.observed_object["nua"]].n_partonlevel == esi.pda[esi.observed_object["nua"]].n_observed_max){
    esi.pda[esi.observed_object["nea"]].n_observed_max = esi.pda[esi.observed_object["nea"]].n_partonlevel;
    esi.pda[esi.observed_object["nma"]].n_observed_max = esi.pda[esi.observed_object["nma"]].n_partonlevel;
    esi.pda[esi.observed_object["nta"]].n_observed_max = esi.pda[esi.observed_object["nta"]].n_partonlevel;
  }
  else if ((esi.pda[esi.observed_object["nea"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nma"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["nta"]].n_observed_max <= esi.pda[esi.observed_object["nua"]].n_observed_max){esi.pda[esi.observed_object["nua"]].n_observed_max = esi.pda[esi.observed_object["nta"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nta"]].n_observed_max = esi.pda[esi.observed_object["nua"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["nea"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nta"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nma"]].n_observed_max > esi.pda[esi.observed_object["nua"]].n_observed_max)){
    if (esi.pda[esi.observed_object["nma"]].n_observed_max <= esi.pda[esi.observed_object["nua"]].n_observed_max){esi.pda[esi.observed_object["nua"]].n_observed_max = esi.pda[esi.observed_object["nma"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nma"]].n_observed_max = esi.pda[esi.observed_object["nua"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["nma"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nta"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nea"]].n_observed_max > esi.pda[esi.observed_object["nua"]].n_observed_max)){
    if (esi.pda[esi.observed_object["nea"]].n_observed_max <= esi.pda[esi.observed_object["nua"]].n_observed_max){esi.pda[esi.observed_object["nua"]].n_observed_max = esi.pda[esi.observed_object["nea"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nea"]].n_observed_max = esi.pda[esi.observed_object["nua"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["nea"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nma"]].n_observed_max + esi.pda[esi.observed_object["nta"]].n_observed_max < esi.pda[esi.observed_object["nua"]].n_observed_max)){
    esi.pda[esi.observed_object["nua"]].n_observed_max = esi.pda[esi.observed_object["nma"]].n_observed_max + esi.pda[esi.observed_object["nta"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["nma"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nea"]].n_observed_max + esi.pda[esi.observed_object["nta"]].n_observed_max < esi.pda[esi.observed_object["nua"]].n_observed_max)){
    esi.pda[esi.observed_object["nua"]].n_observed_max = esi.pda[esi.observed_object["nea"]].n_observed_max + esi.pda[esi.observed_object["nta"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["nta"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nea"]].n_observed_max + esi.pda[esi.observed_object["nma"]].n_observed_max < esi.pda[esi.observed_object["nua"]].n_observed_max)){
    esi.pda[esi.observed_object["nua"]].n_observed_max = esi.pda[esi.observed_object["nea"]].n_observed_max + esi.pda[esi.observed_object["nma"]].n_observed_max;
  }

  // adapt n_observed_min to parton-level particles (nu-ne/nm/nt) !
  if (esi.pda[esi.observed_object["nu"]].n_partonlevel == esi.pda[esi.observed_object["nu"]].n_observed_min){
    esi.pda[esi.observed_object["ne"]].n_observed_min = esi.pda[esi.observed_object["ne"]].n_partonlevel;
    esi.pda[esi.observed_object["nm"]].n_observed_min = esi.pda[esi.observed_object["nm"]].n_partonlevel;
    esi.pda[esi.observed_object["nt"]].n_observed_min = esi.pda[esi.observed_object["nt"]].n_partonlevel;
  }
  else if ((esi.pda[esi.observed_object["ne"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nm"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["nt"]].n_observed_min >= esi.pda[esi.observed_object["nu"]].n_observed_min){esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["nt"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nt"]].n_observed_min = esi.pda[esi.observed_object["nu"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["ne"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nt"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nm"]].n_observed_min > esi.pda[esi.observed_object["nu"]].n_observed_min)){
    if (esi.pda[esi.observed_object["nm"]].n_observed_min >= esi.pda[esi.observed_object["nu"]].n_observed_min){esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["nm"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nm"]].n_observed_min = esi.pda[esi.observed_object["nu"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["nm"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nt"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ne"]].n_observed_min > esi.pda[esi.observed_object["nu"]].n_observed_min)){
    if (esi.pda[esi.observed_object["ne"]].n_observed_min >= esi.pda[esi.observed_object["nu"]].n_observed_min){esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["ne"]].n_observed_min;}
    else {esi.pda[esi.observed_object["ne"]].n_observed_min = esi.pda[esi.observed_object["nu"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["ne"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nm"]].n_observed_min + esi.pda[esi.observed_object["nt"]].n_observed_min > esi.pda[esi.observed_object["nu"]].n_observed_min)){
    esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["nm"]].n_observed_min + esi.pda[esi.observed_object["nt"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["nm"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ne"]].n_observed_min + esi.pda[esi.observed_object["nt"]].n_observed_min > esi.pda[esi.observed_object["nu"]].n_observed_min)){
    esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["ne"]].n_observed_min + esi.pda[esi.observed_object["nt"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["nt"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ne"]].n_observed_min + esi.pda[esi.observed_object["nm"]].n_observed_min > esi.pda[esi.observed_object["nu"]].n_observed_min)){
    esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["ne"]].n_observed_min + esi.pda[esi.observed_object["nm"]].n_observed_min;
  }

  // adapt n_observed_max to parton-level particles (nu-ne/nm/nt) !
  if (esi.pda[esi.observed_object["nu"]].n_partonlevel == esi.pda[esi.observed_object["nu"]].n_observed_max){
    esi.pda[esi.observed_object["ne"]].n_observed_max = esi.pda[esi.observed_object["ne"]].n_partonlevel;
    esi.pda[esi.observed_object["nm"]].n_observed_max = esi.pda[esi.observed_object["nm"]].n_partonlevel;
    esi.pda[esi.observed_object["nt"]].n_observed_max = esi.pda[esi.observed_object["nt"]].n_partonlevel;
  }
  else if ((esi.pda[esi.observed_object["ne"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nm"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["nt"]].n_observed_max <= esi.pda[esi.observed_object["nu"]].n_observed_max){esi.pda[esi.observed_object["nu"]].n_observed_max = esi.pda[esi.observed_object["nt"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nt"]].n_observed_max = esi.pda[esi.observed_object["nu"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["ne"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nt"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nm"]].n_observed_max > esi.pda[esi.observed_object["nu"]].n_observed_max)){
    if (esi.pda[esi.observed_object["nm"]].n_observed_max <= esi.pda[esi.observed_object["nu"]].n_observed_max){esi.pda[esi.observed_object["nu"]].n_observed_max = esi.pda[esi.observed_object["nm"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nm"]].n_observed_max = esi.pda[esi.observed_object["nu"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["nm"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nt"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ne"]].n_observed_max > esi.pda[esi.observed_object["nu"]].n_observed_max)){
    if (esi.pda[esi.observed_object["ne"]].n_observed_max <= esi.pda[esi.observed_object["nu"]].n_observed_max){esi.pda[esi.observed_object["nu"]].n_observed_max = esi.pda[esi.observed_object["ne"]].n_observed_max;}
    else {esi.pda[esi.observed_object["ne"]].n_observed_max = esi.pda[esi.observed_object["nu"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["ne"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nm"]].n_observed_max + esi.pda[esi.observed_object["nt"]].n_observed_max < esi.pda[esi.observed_object["nu"]].n_observed_max)){
    esi.pda[esi.observed_object["nu"]].n_observed_max = esi.pda[esi.observed_object["nm"]].n_observed_max + esi.pda[esi.observed_object["nt"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["nm"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ne"]].n_observed_max + esi.pda[esi.observed_object["nt"]].n_observed_max < esi.pda[esi.observed_object["nu"]].n_observed_max)){
    esi.pda[esi.observed_object["nu"]].n_observed_max = esi.pda[esi.observed_object["ne"]].n_observed_max + esi.pda[esi.observed_object["nt"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["nt"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ne"]].n_observed_max + esi.pda[esi.observed_object["nm"]].n_observed_max < esi.pda[esi.observed_object["nu"]].n_observed_max)){
    esi.pda[esi.observed_object["nu"]].n_observed_max = esi.pda[esi.observed_object["ne"]].n_observed_max + esi.pda[esi.observed_object["nm"]].n_observed_max;
  }

  // adapt n_observed_min to parton-level particles (nux-nex/nmx/ntx) !
  if (esi.pda[esi.observed_object["nux"]].n_partonlevel == esi.pda[esi.observed_object["nux"]].n_observed_min){
    esi.pda[esi.observed_object["nex"]].n_observed_min = esi.pda[esi.observed_object["nex"]].n_partonlevel;
    esi.pda[esi.observed_object["nmx"]].n_observed_min = esi.pda[esi.observed_object["nmx"]].n_partonlevel;
    esi.pda[esi.observed_object["ntx"]].n_observed_min = esi.pda[esi.observed_object["ntx"]].n_partonlevel;
  }
  else if ((esi.pda[esi.observed_object["nex"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nmx"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["ntx"]].n_observed_min >= esi.pda[esi.observed_object["nu"]].n_observed_min){esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["ntx"]].n_observed_min;}
    else {esi.pda[esi.observed_object["ntx"]].n_observed_min = esi.pda[esi.observed_object["nu"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["nex"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nmx"]].n_observed_min > esi.pda[esi.observed_object["nu"]].n_observed_min)){
    if (esi.pda[esi.observed_object["nmx"]].n_observed_min >= esi.pda[esi.observed_object["nu"]].n_observed_min){esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["nmx"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nmx"]].n_observed_min = esi.pda[esi.observed_object["nu"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["nmx"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nex"]].n_observed_min > esi.pda[esi.observed_object["nu"]].n_observed_min)){
    if (esi.pda[esi.observed_object["nex"]].n_observed_min >= esi.pda[esi.observed_object["nu"]].n_observed_min){esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["nex"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nex"]].n_observed_min = esi.pda[esi.observed_object["nu"]].n_observed_min;}
  }
  else if ((esi.pda[esi.observed_object["nex"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nmx"]].n_observed_min + esi.pda[esi.observed_object["ntx"]].n_observed_min > esi.pda[esi.observed_object["nu"]].n_observed_min)){
    esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["nmx"]].n_observed_min + esi.pda[esi.observed_object["ntx"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["nmx"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nex"]].n_observed_min + esi.pda[esi.observed_object["ntx"]].n_observed_min > esi.pda[esi.observed_object["nu"]].n_observed_min)){
    esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["nex"]].n_observed_min + esi.pda[esi.observed_object["ntx"]].n_observed_min;
  }
  else if ((esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nex"]].n_observed_min + esi.pda[esi.observed_object["nmx"]].n_observed_min > esi.pda[esi.observed_object["nu"]].n_observed_min)){
    esi.pda[esi.observed_object["nu"]].n_observed_min = esi.pda[esi.observed_object["nex"]].n_observed_min + esi.pda[esi.observed_object["nmx"]].n_observed_min;
  }

  // adapt n_observed_max to parton-level particles (nux-nex/nmx/ntx) !
  if (esi.pda[esi.observed_object["nux"]].n_partonlevel == esi.pda[esi.observed_object["nux"]].n_observed_max){
    esi.pda[esi.observed_object["nex"]].n_observed_max = esi.pda[esi.observed_object["nex"]].n_partonlevel;
    esi.pda[esi.observed_object["nmx"]].n_observed_max = esi.pda[esi.observed_object["nmx"]].n_partonlevel;
    esi.pda[esi.observed_object["ntx"]].n_observed_max = esi.pda[esi.observed_object["ntx"]].n_partonlevel;
  }
  else if ((esi.pda[esi.observed_object["nex"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nmx"]].n_partonlevel == 0)){
    if (esi.pda[esi.observed_object["ntx"]].n_observed_max <= esi.pda[esi.observed_object["nux"]].n_observed_max){esi.pda[esi.observed_object["nux"]].n_observed_max = esi.pda[esi.observed_object["ntx"]].n_observed_max;}
    else {esi.pda[esi.observed_object["ntx"]].n_observed_max = esi.pda[esi.observed_object["nux"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["nex"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nmx"]].n_observed_max > esi.pda[esi.observed_object["nux"]].n_observed_max)){
    if (esi.pda[esi.observed_object["nmx"]].n_observed_max <= esi.pda[esi.observed_object["nux"]].n_observed_max){esi.pda[esi.observed_object["nux"]].n_observed_max = esi.pda[esi.observed_object["nmx"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nmx"]].n_observed_max = esi.pda[esi.observed_object["nux"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["nmx"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nex"]].n_observed_max > esi.pda[esi.observed_object["nux"]].n_observed_max)){
    if (esi.pda[esi.observed_object["nex"]].n_observed_max <= esi.pda[esi.observed_object["nux"]].n_observed_max){esi.pda[esi.observed_object["nux"]].n_observed_max = esi.pda[esi.observed_object["nex"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nex"]].n_observed_max = esi.pda[esi.observed_object["nux"]].n_observed_max;}
  }
  else if ((esi.pda[esi.observed_object["nex"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nmx"]].n_observed_max + esi.pda[esi.observed_object["ntx"]].n_observed_max < esi.pda[esi.observed_object["nux"]].n_observed_max)){
    esi.pda[esi.observed_object["nux"]].n_observed_max = esi.pda[esi.observed_object["nmx"]].n_observed_max + esi.pda[esi.observed_object["ntx"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["nmx"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nex"]].n_observed_max + esi.pda[esi.observed_object["ntx"]].n_observed_max < esi.pda[esi.observed_object["nux"]].n_observed_max)){
    esi.pda[esi.observed_object["nux"]].n_observed_max = esi.pda[esi.observed_object["nex"]].n_observed_max + esi.pda[esi.observed_object["ntx"]].n_observed_max;
  }
  else if ((esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0) && (esi.pda[esi.observed_object["nex"]].n_observed_max + esi.pda[esi.observed_object["nmx"]].n_observed_max < esi.pda[esi.observed_object["nux"]].n_observed_max)){
    esi.pda[esi.observed_object["nux"]].n_observed_max = esi.pda[esi.observed_object["nex"]].n_observed_max + esi.pda[esi.observed_object["nmx"]].n_observed_max;
  }

  // adapt n_observed_max to parton-level particles (nea-nex/ne) !
  if (esi.pda[esi.observed_object["nea"]].n_partonlevel == esi.pda[esi.observed_object["nea"]].n_observed_max){
    esi.pda[esi.observed_object["ne"]].n_observed_max = esi.pda[esi.observed_object["ne"]].n_partonlevel;
    esi.pda[esi.observed_object["nex"]].n_observed_max = esi.pda[esi.observed_object["nex"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["nex"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["ne"]].n_observed_max <= esi.pda[esi.observed_object["nea"]].n_observed_max){esi.pda[esi.observed_object["nea"]].n_observed_max = esi.pda[esi.observed_object["ne"]].n_observed_max;}
    else {esi.pda[esi.observed_object["ne"]].n_observed_max = esi.pda[esi.observed_object["nea"]].n_observed_max;}
  }
  else if (esi.pda[esi.observed_object["ne"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nex"]].n_observed_max <= esi.pda[esi.observed_object["nea"]].n_observed_max){esi.pda[esi.observed_object["nea"]].n_observed_max = esi.pda[esi.observed_object["nex"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nex"]].n_observed_max = esi.pda[esi.observed_object["nea"]].n_observed_max;}
  }

  // adapt n_observed_min to parton-level particles (nea-nex/ne) !
  if (esi.pda[esi.observed_object["nea"]].n_partonlevel == esi.pda[esi.observed_object["nea"]].n_observed_min){
    esi.pda[esi.observed_object["ne"]].n_observed_min = esi.pda[esi.observed_object["ne"]].n_partonlevel;
    esi.pda[esi.observed_object["nex"]].n_observed_min = esi.pda[esi.observed_object["nex"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["nex"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["ne"]].n_observed_min >= esi.pda[esi.observed_object["nea"]].n_observed_min){esi.pda[esi.observed_object["nea"]].n_observed_min = esi.pda[esi.observed_object["ne"]].n_observed_min;}
    else {esi.pda[esi.observed_object["ne"]].n_observed_min = esi.pda[esi.observed_object["nea"]].n_observed_min;}
  }
  else if (esi.pda[esi.observed_object["ne"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nex"]].n_observed_min >= esi.pda[esi.observed_object["nea"]].n_observed_min){esi.pda[esi.observed_object["nea"]].n_observed_min = esi.pda[esi.observed_object["nex"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nex"]].n_observed_min = esi.pda[esi.observed_object["nea"]].n_observed_min;}
  }

  // adapt n_observed_max to parton-level particles (nma-nmx/nm) !
  if (esi.pda[esi.observed_object["nma"]].n_partonlevel == esi.pda[esi.observed_object["nma"]].n_observed_max){
    esi.pda[esi.observed_object["nm"]].n_observed_max = esi.pda[esi.observed_object["nm"]].n_partonlevel;
    esi.pda[esi.observed_object["nmx"]].n_observed_max = esi.pda[esi.observed_object["nmx"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["nmx"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nm"]].n_observed_max <= esi.pda[esi.observed_object["nma"]].n_observed_max){esi.pda[esi.observed_object["nma"]].n_observed_max = esi.pda[esi.observed_object["nm"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nm"]].n_observed_max = esi.pda[esi.observed_object["nma"]].n_observed_max;}
  }
  else if (esi.pda[esi.observed_object["nm"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nmx"]].n_observed_max <= esi.pda[esi.observed_object["nma"]].n_observed_max){esi.pda[esi.observed_object["nma"]].n_observed_max = esi.pda[esi.observed_object["nmx"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nmx"]].n_observed_max = esi.pda[esi.observed_object["nma"]].n_observed_max;}
  }

  // adapt n_observed_min to parton-level particles (nma-nmx/nm) !
  if (esi.pda[esi.observed_object["nma"]].n_partonlevel == esi.pda[esi.observed_object["nma"]].n_observed_min){
    esi.pda[esi.observed_object["nm"]].n_observed_min = esi.pda[esi.observed_object["nm"]].n_partonlevel;
    esi.pda[esi.observed_object["nmx"]].n_observed_min = esi.pda[esi.observed_object["nmx"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["nmx"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nm"]].n_observed_min >= esi.pda[esi.observed_object["nma"]].n_observed_min){esi.pda[esi.observed_object["nma"]].n_observed_min = esi.pda[esi.observed_object["nm"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nm"]].n_observed_min = esi.pda[esi.observed_object["nma"]].n_observed_min;}
  }
  else if (esi.pda[esi.observed_object["nm"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nmx"]].n_observed_min >= esi.pda[esi.observed_object["nma"]].n_observed_min){esi.pda[esi.observed_object["nma"]].n_observed_min = esi.pda[esi.observed_object["nmx"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nmx"]].n_observed_min = esi.pda[esi.observed_object["nma"]].n_observed_min;}
  }

  // adapt n_observed_max to parton-level particles (nta-ntx/nt) !
  if (esi.pda[esi.observed_object["nta"]].n_partonlevel == esi.pda[esi.observed_object["nta"]].n_observed_max){
    esi.pda[esi.observed_object["nt"]].n_observed_max = esi.pda[esi.observed_object["nt"]].n_partonlevel;
    esi.pda[esi.observed_object["ntx"]].n_observed_max = esi.pda[esi.observed_object["ntx"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nt"]].n_observed_max <= esi.pda[esi.observed_object["nta"]].n_observed_max){esi.pda[esi.observed_object["nta"]].n_observed_max = esi.pda[esi.observed_object["nt"]].n_observed_max;}
    else {esi.pda[esi.observed_object["nt"]].n_observed_max = esi.pda[esi.observed_object["nta"]].n_observed_max;}
  }
  else if (esi.pda[esi.observed_object["nt"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["ntx"]].n_observed_max <= esi.pda[esi.observed_object["nta"]].n_observed_max){esi.pda[esi.observed_object["nta"]].n_observed_max = esi.pda[esi.observed_object["ntx"]].n_observed_max;}
    else {esi.pda[esi.observed_object["ntx"]].n_observed_max = esi.pda[esi.observed_object["nta"]].n_observed_max;}
  }

  // adapt n_observed_min to parton-level particles (nta-ntx/nt) !
  if (esi.pda[esi.observed_object["nta"]].n_partonlevel == esi.pda[esi.observed_object["nta"]].n_observed_min){
    esi.pda[esi.observed_object["nt"]].n_observed_min = esi.pda[esi.observed_object["nt"]].n_partonlevel;
    esi.pda[esi.observed_object["ntx"]].n_observed_min = esi.pda[esi.observed_object["ntx"]].n_partonlevel;
  }
  else if (esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["nt"]].n_observed_min >= esi.pda[esi.observed_object["nta"]].n_observed_min){esi.pda[esi.observed_object["nta"]].n_observed_min = esi.pda[esi.observed_object["nt"]].n_observed_min;}
    else {esi.pda[esi.observed_object["nt"]].n_observed_min = esi.pda[esi.observed_object["nta"]].n_observed_min;}
  }
  else if (esi.pda[esi.observed_object["nt"]].n_partonlevel == 0){
    if (esi.pda[esi.observed_object["ntx"]].n_observed_min >= esi.pda[esi.observed_object["nta"]].n_observed_min){esi.pda[esi.observed_object["nta"]].n_observed_min = esi.pda[esi.observed_object["ntx"]].n_observed_min;}
    else {esi.pda[esi.observed_object["ntx"]].n_observed_min = esi.pda[esi.observed_object["nta"]].n_observed_min;}
  }


  logger << LOG_INFO << endl << "Object definition after taking into account parton content: " << endl;
  logger << LOG_INFO << setw(20) << "object" << setw(20) << "n_observed_min" << setw(20) << "n_observed_max" << setw(20) << "define_pT" << setw(20) << "define_ET" << setw(20) << "define_eta"  << setw(20) << "define_y"<< endl;
  for (int i = 1; i < esi.object_list.size(); i++){
    if (esi.pda[i].n_observed_max > 0){
      logger << LOG_INFO << setw(20) << esi.object_list[i] << setw(20) << esi.pda[i].n_observed_min << setw(20) << esi.pda[i].n_observed_max << endl;// << setw(20) << define_pT[i] << setw(20) << define_ET[i] << setw(20) << define_eta[i] << setw(20) << define_y[i] << endl;
      logger << LOG_INFO << setw(20) << esi.object_list[i] << setw(20) << esi.pda[i].n_observed_min << setw(20) << esi.pda[i].n_observed_max << endl;// << setw(20) << define_pT[i] << setw(20) << define_ET[i] << setw(20) << define_eta[i] << setw(20) << define_y[i] << endl;
    }
    else {
      logger << LOG_DEBUG << setw(20) << esi.object_list[i] << setw(20) << esi.pda[i].n_observed_min << setw(20) << esi.pda[i].n_observed_max << endl;
      logger << LOG_DEBUG << setw(20) << esi.object_list[i] << setw(20) << esi.pda[i].n_observed_min << setw(20) << esi.pda[i].n_observed_max << endl;
    }
  }
  logger.newLine(LOG_INFO);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

  

void observable_set::determine_equivalent_object(){
  static  Logger logger("observable_set::determine_equivalent_object");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  map<string, int> protect_observable_object;
  
  for (int i = 1; i < esi.object_list.size(); i++){
    if (esi.object_list[i] == "ljet" && esi.pda[esi.observed_object["jet"]].n_partonlevel != 0){equivalent_object["ljet"] = "ljet";}
    else if (esi.pda[esi.observed_object[esi.object_list[i]]].n_partonlevel == 0 && esi.object_category[i] != -1){equivalent_object[esi.object_list[i]] = "none";}
    else {equivalent_object[esi.object_list[i]] = esi.object_list[i];}
  }

  // check if "tjet" is identical to "top" or "atop" etc.
  if (esi.pda[esi.observed_object["tjet"]].n_partonlevel == esi.pda[esi.observed_object["atop"]].n_partonlevel && 
      esi.pda[esi.observed_object["atop"]].n_partonlevel > 0 && esi.pda[esi.observed_object["top"]].n_partonlevel == 0){equivalent_object["atop"] = equivalent_object["tjet"];}
  else if (esi.pda[esi.observed_object["tjet"]].n_partonlevel == esi.pda[esi.observed_object["top"]].n_partonlevel && 
	   esi.pda[esi.observed_object["top"]].n_partonlevel > 0 && esi.pda[esi.observed_object["atop"]].n_partonlevel == 0){equivalent_object["top"] = equivalent_object["tjet"];}

  /*
  // too involved - better never identify them !!!
  // only if less than one jet - otherwise jet algorithm might change the flavour...
  if (esi.pda[esi.observed_object["jet"]].n_partonlevel < 2){
    // check if "jet" is identical to "ljet" or "bjet" etc.
    if (esi.pda[esi.observed_object["jet"]].n_partonlevel == esi.pda[esi.observed_object["bjet"]].n_partonlevel && 
	esi.pda[esi.observed_object["bjet"]].n_partonlevel > 0 && esi.pda[esi.observed_object["ljet"]].n_partonlevel == 0){equivalent_object["bjet"] = equivalent_object["jet"];}
  }

  if (esi.pda[esi.observed_object["jet"]].n_partonlevel == esi.pda[esi.observed_object["ljet"]].n_partonlevel && 
      esi.pda[esi.observed_object["ljet"]].n_partonlevel > 0 && esi.pda[esi.observed_object["bjet"]].n_partonlevel == 0){equivalent_object["ljet"] = equivalent_object["jet"];}
  */

  // check if "bjet" is identical to "bjet_b" or "bjet_bx" etc.
  if (esi.pda[esi.observed_object["bjet"]].n_partonlevel == esi.pda[esi.observed_object["bjet_bx"]].n_partonlevel && 
      esi.pda[esi.observed_object["bjet_bx"]].n_partonlevel > 0 && esi.pda[esi.observed_object["bjet_b"]].n_partonlevel == 0){equivalent_object["bjet_bx"] = equivalent_object["bjet"];}
  else if (esi.pda[esi.observed_object["bjet"]].n_partonlevel == esi.pda[esi.observed_object["bjet_b"]].n_partonlevel && 
	   esi.pda[esi.observed_object["bjet_b"]].n_partonlevel > 0 && esi.pda[esi.observed_object["bjet_bx"]].n_partonlevel == 0){equivalent_object["bjet_b"] = equivalent_object["bjet"];}


  // check if "lep" is identical to "lm" or "lp" etc.
  if (esi.pda[esi.observed_object["lep"]].n_partonlevel == esi.pda[esi.observed_object["lp"]].n_partonlevel && 
      esi.pda[esi.observed_object["lp"]].n_partonlevel > 0 && esi.pda[esi.observed_object["lm"]].n_partonlevel == 0){equivalent_object["lp"] = equivalent_object["lep"];}
  else if (esi.pda[esi.observed_object["lep"]].n_partonlevel == esi.pda[esi.observed_object["lm"]].n_partonlevel && 
	   esi.pda[esi.observed_object["lm"]].n_partonlevel > 0 && esi.pda[esi.observed_object["lp"]].n_partonlevel == 0){equivalent_object["lm"] = equivalent_object["lep"];}

  // check if "lep" is identical to "e" or "mu" or "tau" etc.
  if (esi.pda[esi.observed_object["lep"]].n_partonlevel == esi.pda[esi.observed_object["e"]].n_partonlevel && 
      esi.pda[esi.observed_object["e"]].n_partonlevel > 0 && 
      esi.pda[esi.observed_object["mu"]].n_partonlevel == 0 && 
      esi.pda[esi.observed_object["tau"]].n_partonlevel == 0){equivalent_object["e"] = equivalent_object["lep"];}
  else if (esi.pda[esi.observed_object["lep"]].n_partonlevel == esi.pda[esi.observed_object["mu"]].n_partonlevel && 
	   esi.pda[esi.observed_object["mu"]].n_partonlevel > 0 && 
	   esi.pda[esi.observed_object["e"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["tau"]].n_partonlevel == 0){equivalent_object["mu"] = equivalent_object["lep"];}
  else if (esi.pda[esi.observed_object["lep"]].n_partonlevel == esi.pda[esi.observed_object["tau"]].n_partonlevel && 
	   esi.pda[esi.observed_object["e"]].n_partonlevel > 0 && 
	   esi.pda[esi.observed_object["mu"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["tau"]].n_partonlevel == 0){equivalent_object["tau"] = equivalent_object["lep"];}

  // check if "e" is identical to "ep" or "em" etc.
  if (esi.pda[esi.observed_object["e"]].n_partonlevel == esi.pda[esi.observed_object["em"]].n_partonlevel && 
      esi.pda[esi.observed_object["em"]].n_partonlevel > 0 && esi.pda[esi.observed_object["ep"]].n_partonlevel == 0){equivalent_object["em"] = equivalent_object["e"];}
  else if (esi.pda[esi.observed_object["e"]].n_partonlevel == esi.pda[esi.observed_object["ep"]].n_partonlevel && 
	   esi.pda[esi.observed_object["ep"]].n_partonlevel > 0 && esi.pda[esi.observed_object["em"]].n_partonlevel == 0){equivalent_object["ep"] = equivalent_object["e"];}

  // check if "mu" is identical to "mup" or "mum" etc.
  if (esi.pda[esi.observed_object["mu"]].n_partonlevel == esi.pda[esi.observed_object["mum"]].n_partonlevel && 
      esi.pda[esi.observed_object["mum"]].n_partonlevel > 0 && esi.pda[esi.observed_object["mup"]].n_partonlevel == 0){equivalent_object["mum"] = equivalent_object["mu"];}
  else if (esi.pda[esi.observed_object["mu"]].n_partonlevel == esi.pda[esi.observed_object["mup"]].n_partonlevel && 
	   esi.pda[esi.observed_object["mup"]].n_partonlevel > 0 && esi.pda[esi.observed_object["mum"]].n_partonlevel == 0){equivalent_object["mup"] = equivalent_object["mu"];}

  // check if "tau" is identical to "taup" or "taum" etc.
  if (esi.pda[esi.observed_object["tau"]].n_partonlevel == esi.pda[esi.observed_object["taum"]].n_partonlevel && 
      esi.pda[esi.observed_object["taum"]].n_partonlevel > 0 && esi.pda[esi.observed_object["taup"]].n_partonlevel == 0){equivalent_object["taum"] = equivalent_object["tau"];}
  else if (esi.pda[esi.observed_object["tau"]].n_partonlevel == esi.pda[esi.observed_object["taup"]].n_partonlevel && 
	   esi.pda[esi.observed_object["taup"]].n_partonlevel > 0 && esi.pda[esi.observed_object["taum"]].n_partonlevel == 0){equivalent_object["taup"] = equivalent_object["tau"];}

  // check if "lm" is identical to "em" or "mum" or "taum" etc.
  if (esi.pda[esi.observed_object["lm"]].n_partonlevel == esi.pda[esi.observed_object["em"]].n_partonlevel && 
      esi.pda[esi.observed_object["em"]].n_partonlevel > 0 && 
      esi.pda[esi.observed_object["mum"]].n_partonlevel == 0 && 
      esi.pda[esi.observed_object["taum"]].n_partonlevel == 0){equivalent_object["lm"] = equivalent_object["em"];}
  else if (esi.pda[esi.observed_object["lm"]].n_partonlevel == esi.pda[esi.observed_object["mum"]].n_partonlevel && 
	   esi.pda[esi.observed_object["mum"]].n_partonlevel > 0 && 
	   esi.pda[esi.observed_object["em"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["taum"]].n_partonlevel == 0){equivalent_object["lm"] = equivalent_object["mum"];}
  else if (esi.pda[esi.observed_object["lm"]].n_partonlevel == esi.pda[esi.observed_object["taum"]].n_partonlevel && 
	   esi.pda[esi.observed_object["em"]].n_partonlevel > 0 && 
	   esi.pda[esi.observed_object["mum"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["taum"]].n_partonlevel == 0){equivalent_object["lm"] = equivalent_object["taum"];}

  // check if "lp" is identical to "ep" or "mup" or "taup" etc.
  if (esi.pda[esi.observed_object["lp"]].n_partonlevel == esi.pda[esi.observed_object["ep"]].n_partonlevel && 
      esi.pda[esi.observed_object["ep"]].n_partonlevel > 0 && 
      esi.pda[esi.observed_object["mup"]].n_partonlevel == 0 && 
      esi.pda[esi.observed_object["taup"]].n_partonlevel == 0){equivalent_object["lp"] = equivalent_object["ep"];}
  else if (esi.pda[esi.observed_object["lp"]].n_partonlevel == esi.pda[esi.observed_object["mup"]].n_partonlevel && 
	   esi.pda[esi.observed_object["mup"]].n_partonlevel > 0 && 
	   esi.pda[esi.observed_object["ep"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["taup"]].n_partonlevel == 0){equivalent_object["lp"] = equivalent_object["mup"];}
  else if (esi.pda[esi.observed_object["lp"]].n_partonlevel == esi.pda[esi.observed_object["taup"]].n_partonlevel && 
	   esi.pda[esi.observed_object["ep"]].n_partonlevel > 0 && 
	   esi.pda[esi.observed_object["mup"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["taup"]].n_partonlevel == 0){equivalent_object["lp"] = equivalent_object["taup"];}


  // check if "nea" is identical to "nex" or "ne" etc.
  if (esi.pda[esi.observed_object["nea"]].n_partonlevel == esi.pda[esi.observed_object["ne"]].n_partonlevel && 
      esi.pda[esi.observed_object["ne"]].n_partonlevel > 0 && esi.pda[esi.observed_object["nex"]].n_partonlevel == 0){equivalent_object["nea"] = equivalent_object["ne"];}
  else if (esi.pda[esi.observed_object["nea"]].n_partonlevel == esi.pda[esi.observed_object["nex"]].n_partonlevel && 
	   esi.pda[esi.observed_object["nex"]].n_partonlevel > 0 && esi.pda[esi.observed_object["ne"]].n_partonlevel == 0){equivalent_object["nea"] = equivalent_object["nex"];}
  if (esi.pda[esi.observed_object["nma"]].n_partonlevel == esi.pda[esi.observed_object["nm"]].n_partonlevel && 
      esi.pda[esi.observed_object["nm"]].n_partonlevel > 0 && esi.pda[esi.observed_object["nmx"]].n_partonlevel == 0){equivalent_object["nma"] = equivalent_object["nm"];}
  else if (esi.pda[esi.observed_object["nma"]].n_partonlevel == esi.pda[esi.observed_object["nmx"]].n_partonlevel && 
	   esi.pda[esi.observed_object["nmx"]].n_partonlevel > 0 && esi.pda[esi.observed_object["nm"]].n_partonlevel == 0){equivalent_object["nma"] = equivalent_object["nmx"];}
  if (esi.pda[esi.observed_object["nta"]].n_partonlevel == esi.pda[esi.observed_object["nt"]].n_partonlevel && 
      esi.pda[esi.observed_object["nt"]].n_partonlevel > 0 && esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0){equivalent_object["nta"] = equivalent_object["nt"];}
  else if (esi.pda[esi.observed_object["nta"]].n_partonlevel == esi.pda[esi.observed_object["ntx"]].n_partonlevel && 
	   esi.pda[esi.observed_object["ntx"]].n_partonlevel > 0 && esi.pda[esi.observed_object["nt"]].n_partonlevel == 0){equivalent_object["nta"] = equivalent_object["ntx"];}

  if (esi.pda[esi.observed_object["nu"]].n_partonlevel == esi.pda[esi.observed_object["ne"]].n_partonlevel && 
      esi.pda[esi.observed_object["ne"]].n_partonlevel > 0 && esi.pda[esi.observed_object["nm"]].n_partonlevel == 0 && 
      esi.pda[esi.observed_object["nt"]].n_partonlevel == 0){equivalent_object["nu"] = equivalent_object["ne"];}
  else if (esi.pda[esi.observed_object["nu"]].n_partonlevel == esi.pda[esi.observed_object["nm"]].n_partonlevel && 
	   esi.pda[esi.observed_object["nm"]].n_partonlevel > 0 && esi.pda[esi.observed_object["ne"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["nt"]].n_partonlevel == 0){equivalent_object["nu"] = equivalent_object["nm"];}
  else if (esi.pda[esi.observed_object["nu"]].n_partonlevel == esi.pda[esi.observed_object["nt"]].n_partonlevel && 
	   esi.pda[esi.observed_object["ne"]].n_partonlevel > 0 && esi.pda[esi.observed_object["nm"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["nt"]].n_partonlevel == 0){equivalent_object["nu"] = equivalent_object["nt"];}

  if (esi.pda[esi.observed_object["nux"]].n_partonlevel == esi.pda[esi.observed_object["nex"]].n_partonlevel && 
      esi.pda[esi.observed_object["nex"]].n_partonlevel > 0 && 
      esi.pda[esi.observed_object["nmx"]].n_partonlevel == 0 && 
      esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0){equivalent_object["nux"] = equivalent_object["nex"];}
  else if (esi.pda[esi.observed_object["nux"]].n_partonlevel == esi.pda[esi.observed_object["nmx"]].n_partonlevel && 
	   esi.pda[esi.observed_object["nmx"]].n_partonlevel > 0 && 
	   esi.pda[esi.observed_object["nex"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0){equivalent_object["nux"] = equivalent_object["nmx"];}
  else if (esi.pda[esi.observed_object["nux"]].n_partonlevel == esi.pda[esi.observed_object["ntx"]].n_partonlevel && 
	   esi.pda[esi.observed_object["nex"]].n_partonlevel > 0 && 
	   esi.pda[esi.observed_object["nmx"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["ntx"]].n_partonlevel == 0){equivalent_object["nux"] = equivalent_object["ntx"];}

  if (esi.pda[esi.observed_object["nua"]].n_partonlevel == esi.pda[esi.observed_object["nea"]].n_partonlevel && 
      esi.pda[esi.observed_object["nea"]].n_partonlevel > 0 && 
      esi.pda[esi.observed_object["nma"]].n_partonlevel == 0 && 
      esi.pda[esi.observed_object["nta"]].n_partonlevel == 0){equivalent_object["nua"] = equivalent_object["nea"];}
  else if (esi.pda[esi.observed_object["nua"]].n_partonlevel == esi.pda[esi.observed_object["nma"]].n_partonlevel && 
	   esi.pda[esi.observed_object["nma"]].n_partonlevel > 0 && 
	   esi.pda[esi.observed_object["nea"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["nta"]].n_partonlevel == 0){equivalent_object["nua"] = equivalent_object["nma"];}
  else if (esi.pda[esi.observed_object["nua"]].n_partonlevel == esi.pda[esi.observed_object["nta"]].n_partonlevel && 
	   esi.pda[esi.observed_object["nea"]].n_partonlevel > 0 && 
	   esi.pda[esi.observed_object["nma"]].n_partonlevel == 0 && 
	   esi.pda[esi.observed_object["nta"]].n_partonlevel == 0){equivalent_object["nua"] = equivalent_object["nta"];}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_relevant_object(){
  static  Logger logger("observable_set::determine_relevant_object");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_INFO << "esi.pds.size() = " << esi.pds.size() << endl;

  esi.object_list_selection.push_back("all");
  esi.pds.push_back(define_particle_set ());

  esi.object_list_selection.push_back("none");
  esi.pds.push_back(define_particle_set ());

  logger << LOG_DEBUG << setw(5) << "no" << setw(10) << "object" << setw(10) << "n_min" << setw(10) << "n_max" << setw(10) << "pT" << setw(10) << "ET" << setw(10) << "|eta|"  << setw(10) << "|y|"<< endl;
  for (int i = 1; i < esi.object_list.size(); i++){
    if (equivalent_object[esi.object_list[i]] != "none"){
      //      logger << LOG_INFO << "esi.object_list[" << setw(3) << i << "] = " << setw(10) << esi.object_list[i] << "   ->   equivalent_object[" << setw(3) << setw(10) << esi.object_list[i] << "] = " << setw(10) << equivalent_object[esi.object_list[i]] << endl;
    }
    else {
      logger << LOG_DEBUG_VERBOSE << "esi.object_list[" << setw(3) << i << "] = " << setw(10) << esi.object_list[i] << "   ->   equivalent_object[" << setw(3) << setw(10) << esi.object_list[i] << "] = " << setw(10) << equivalent_object[esi.object_list[i]] << endl;
    }

    int flag = 0;
    if (equivalent_object[esi.object_list[i]] != ""){ 
      for (int j = 0; j < esi.object_list_selection.size(); j++){
	if (equivalent_object[esi.object_list[i]] == esi.object_list_selection[j]){flag = 1; break;}
      }
      if (flag == 0){
	logger << LOG_DEBUG_VERBOSE << "esi.object_list_selection.size() = " << esi.object_list_selection.size() << endl;

	int x_i = 0;
	for (int i_i = 0; i_i < esi.object_list.size(); i_i++){
	  if (esi.object_list[i_i] == equivalent_object[esi.object_list[i]]){
	    x_i = i_i; break;
	  }
	}

	logger << LOG_DEBUG << setw(5) << i << setw(10) << esi.object_list[i] << setw(10) << esi.pda[i].n_observed_min << setw(10) << esi.pda[i].n_observed_max << setw(10) << esi.pda[i].define_pT << setw(10) << esi.pda[i].define_ET << setw(10) << esi.pda[i].define_eta << setw(10) << esi.pda[i].define_y << endl;

	esi.object_list_selection.push_back(esi.object_list[x_i]);
	esi.pds.push_back(esi.pda[x_i]);
	esi.object_category_selection.push_back(x_i);  //  ???
      }
    }
  }

  // add objects without particle equivalent (in order to exclude contributions that do not fulfill the corresponding requirements)
  // Could be avoided by checking for full object_list to select vanishing contributions !!!
  for (int x_i = 1; x_i < esi.object_list.size(); x_i++){
    if ((esi.pda[x_i].n_observed_min > esi.pda[x_i].n_partonlevel) && equivalent_object[esi.object_list[x_i]] == "none"){
      logger << LOG_INFO << "    esi.object_list[" << setw(3) << x_i << "] = " << setw(10) << esi.object_list[x_i] << "   ->   equivalent_object[" << setw(3) << setw(10) << esi.object_list[x_i] << "] = " << setw(10) << equivalent_object[esi.object_list[x_i]] << "   n_observed_min[" << setw(3) << x_i << "] = " << esi.pda[x_i].n_observed_max  << "   n_observed_min[" << setw(3) << x_i << "] = " << esi.pda[x_i].n_observed_max << endl;
      logger << LOG_INFO << "No contribution can result from this subprocess - modification needed due to type conversion (e.g. jet/photon recombination alogorithm)" << endl;
      //      exit(1);

      esi.object_list_selection.push_back(esi.object_list[x_i]);
      esi.pds.push_back(esi.pda[x_i]);
      esi.object_category_selection.push_back(x_i);  //  ???
    }
  }

  logger.newLine(LOG_INFO);

  for (int i = 0; i < esi.object_list_selection.size(); i++){no_relevant_object[esi.object_list_selection[i]] = i;}

  //  equivalent_no_object.resize(esi.object_list.size());
  for (int i_o = 0; i_o < esi.object_list.size(); i_o++){
    for (int j_o = 0; j_o < esi.object_list_selection.size(); j_o++){
      if (esi.object_list[i_o] == esi.object_list_selection[j_o]){
	equivalent_no_object[j_o] = i_o;
      }
    }
  }

  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "esi.object_list_selection.size() = " << esi.object_list_selection.size() << endl;
  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG << "Relevant objects:" << endl;
  logger << LOG_DEBUG << setw(5) << "no" << setw(10) << "object" << setw(10) << "n_min" << setw(10) << "n_max" << setw(10) << "pT" << setw(10) << "ET" << setw(10) << "|eta|"  << setw(10) << "|y|"<< endl;
  for (int i = 0; i < esi.object_list_selection.size(); i++){
    logger << LOG_DEBUG << setw(5) << i << setw(10) << esi.object_list_selection[i] << setw(10) << esi.pds[i].n_observed_min << setw(10) << esi.pds[i].n_observed_max << setw(10) << esi.pds[i].define_pT << setw(10) << esi.pds[i].define_ET << setw(10) << esi.pds[i].define_eta << setw(10) << esi.pds[i].define_y << endl;
  }
  logger.newLine(LOG_DEBUG);

  /*
  // maybe later...
  for (int i = 0; i < esi.object_list.size(); i++){
    if (no_relevant_object[esi.object_list[i]] == 0){
      no_relevant_object[esi.object_list[i]] = no_relevant_object[equivalent_object[esi.object_list[i]]];
    }
  }
  */
  //  for (int i = 0; i < relevant_object_list.size(); i++){no_relevant_object[equivalent_object[relevant_object_list[i]]] = i;}

  for (int i_o = 0; i_o < esi.object_list_selection.size(); i_o++){
    esi.observed_object_selection[esi.object_list_selection[i_o]] = i_o;
  }

  for (int i_o = 0; i_o < esi.object_list_selection.size(); i_o++){
    logger << LOG_INFO << "esi.object_list_selection[" << i_o << "] = " << esi.object_list_selection[i_o] << "   ---   " << esi.observed_object_selection[esi.object_list_selection[i_o]] << endl;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_runtime_object(){
  static  Logger logger("observable_set::determine_runtime_object");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  // relevance for jet algorithm:
  // 1 = only (one type of) jets
  // 2, 3, ... = more than one type of jets
  // -1 = only (one type of) jets and photons 
  // -2, -3, ... =  more than one type of jets and photons

  //  vector<int> temp_original(relevant_object_list.size(), 0);
  vector<int> temp_original(esi.object_list_selection.size(), 0);

  //  if (jet_algorithm || frixione_isolation){
  {
    int temp_sign = 1;
    for (int i = 1; i < esi.object_list_selection.size(); i++){
      if (esi.object_list_selection[i] == "tjet" || 
	  esi.object_list_selection[i] == "top"){
	int flag = 0;
	for (int i_j = 0; i_j < jet_algorithm_list.size(); i_j++){
	  if (jet_algorithm_list[i_j] == 6){flag++; break;}
	}
	if (flag){
	  runtime_jet_recombination.push_back(i);
	  temp_original[i]++;
	  runtime_jet_algorithm++;
	}
      }
      else if (esi.object_list_selection[i] == "tjet" || 
	       esi.object_list_selection[i] == "atop"){
	int flag = 0;
	for (int i_j = 0; i_j < jet_algorithm_list.size(); i_j++){
	  if (jet_algorithm_list[i_j] == -6){flag++; break;}
	}
	if (flag){
	  runtime_jet_recombination.push_back(i);
	  temp_original[i]++;
	  runtime_jet_algorithm++;
	}
      }
      else if (esi.object_list_selection[i] == "bjet" || 
	  esi.object_list_selection[i] == "bjet_b"){
	int flag = 0;
	for (int i_j = 0; i_j < jet_algorithm_list.size(); i_j++){
	  if (jet_algorithm_list[i_j] == 5){flag++; break;}
	}
	if (flag){
	  runtime_jet_recombination.push_back(i);
	  temp_original[i]++;
	  runtime_jet_algorithm++;
	}
      }
      else if (esi.object_list_selection[i] == "bjet" || 
	       esi.object_list_selection[i] == "bjet_bx"){
	int flag = 0;
	for (int i_j = 0; i_j < jet_algorithm_list.size(); i_j++){
	  if (jet_algorithm_list[i_j] == -5){flag++; break;}
	}
	if (flag){
	  runtime_jet_recombination.push_back(i);
	  temp_original[i]++;
	  runtime_jet_algorithm++;
	}
      }
      else if (esi.object_list_selection[i] == "jet" || 
	       esi.object_list_selection[i] == "ljet" || 
	       (esi.object_list_selection[i] == "photon" && !frixione_isolation && (jet_algorithm || photon_jet_algorithm))){
	//	       (esi.object_list_selection[i] == "photon" && frixione_isolation == 0)){
	runtime_jet_recombination.push_back(i);
	temp_original[i]++;
	runtime_jet_algorithm++;
      }
      if (esi.object_list_selection[i] == "photon"){temp_sign = -1;}
    }
    runtime_jet_algorithm = temp_sign * runtime_jet_algorithm;
  }

  if (photon_recombination != 0){
    if (esi.pda[esi.observed_object["photon"]].n_partonlevel != 0){
      for (int i = 1; i < esi.object_list_selection.size(); i++){
	if (esi.object_list_selection[i] == "jet" || 
	    esi.object_list_selection[i] == "ljet" || 
	    esi.object_list_selection[i] == "bjet" || 
	    esi.object_list_selection[i] == "bjet_b" || 
	    esi.object_list_selection[i] == "bjet_bx" || 
	    esi.object_list_selection[i] == "lep" || 
	    esi.object_list_selection[i] == "lm" || 
	    esi.object_list_selection[i] == "lp" || 
	    esi.object_list_selection[i] == "e" || 
	    esi.object_list_selection[i] == "mu" || 
	    esi.object_list_selection[i] == "tau" || 
	    esi.object_list_selection[i] == "em" || 
	    esi.object_list_selection[i] == "mum" || 
	    esi.object_list_selection[i] == "taum" || 
	    esi.object_list_selection[i] == "ep" || 
	    esi.object_list_selection[i] == "mup" || 
	    esi.object_list_selection[i] == "taup" || 
	    esi.object_list_selection[i] == "photon"){
	  runtime_photon_recombination.push_back(i);

	  // leptons should be written to "particles" after the photon-recombination phase
	  if (esi.object_list_selection[i] == "lep" || 
	      esi.object_list_selection[i] == "lm" || 
	      esi.object_list_selection[i] == "lp" || 
	      esi.object_list_selection[i] == "e" || 
	      esi.object_list_selection[i] == "mu" || 
	      esi.object_list_selection[i] == "tau" || 
	      esi.object_list_selection[i] == "em" || 
	      esi.object_list_selection[i] == "mum" || 
	      esi.object_list_selection[i] == "taum" || 
	      esi.object_list_selection[i] == "ep" || 
	      esi.object_list_selection[i] == "mup" || 
	      esi.object_list_selection[i] == "taup"){
	    // Handled differently at the moment...
	  }
	  // temp_original statement: why should all particles but leptons be affected by this ???
	  else {temp_original[i]++;}
	  //      runtime_photon_algorithm++;
	}
	//    if (esi.object_list_selection[i] == "photon"){temp_sign = -1;}
      }
    }
  }
  //  runtime_jet_algorithm = temp_sign * runtime_jet_algorithm;

  if (frixione_isolation != 0){
    for (int i = 1; i < esi.object_list_selection.size(); i++){
      if (esi.object_list_selection[i] == "photon"){
	runtime_photon_isolation.push_back(i);
	temp_original[i]++;
	//      runtime_photon_algorithm++;
      }
      // new: partons could be changed by photon isolation (e.g. removed from list)
	if (esi.object_list_selection[i] == "jet" || 
	    esi.object_list_selection[i] == "ljet" || 
	    esi.object_list_selection[i] == "bjet" || 
	    esi.object_list_selection[i] == "bjet_b" || 
	    esi.object_list_selection[i] == "bjet_bx"){temp_original[i]++;}

      //    if (esi.object_list_selection[i] == "photon"){temp_sign = -1;}
    }
  }

  for (int i = 1; i < esi.object_list_selection.size(); i++){
    if (esi.object_list_selection[i] == "missing"){
      //      runtime_missing.push_back(i);
      temp_original[i]++;
      //      runtime_photon_algorithm++;
    }
    //    if (esi.object_list_selection[i] == "photon"){temp_sign = -1;}
  }


  vector<vector<fourvector> > momentum_object;
  momentum_object.resize(esi.object_list_selection.size());

  // lepton--photon recombination ignored at the moment -> all other partons enter event selection with their original momenta
  for (int i = 1; i < esi.object_list_selection.size(); i++){
    //    cout << "temp_original[" << setw(3) << i << "] = " << temp_original[i] << endl;
    if (temp_original[i] == 0 && esi.object_category[equivalent_no_object[i]] != -1){
      runtime_original.push_back(i);
    }
  }

  runtime_order.resize(esi.object_list_selection.size());
  runtime_order_inverse.resize(csi->type_parton[0].size());
  for (int i = 1; i < esi.object_list_selection.size(); i++){
    for (int i_p = 3; i_p < csi->type_parton[0].size(); i_p++){
      if (esi.object_list_selection[i] == "jet"){
	//      if (esi.object_list_selection[i] == "jet" && abs(csi->type_parton[0][i_p]) < 6){
	int i_a = 0;  //  temporary ???
	int flag = 0;
	for (int i_j = 0; i_j < jet_algorithm_list.size(); i_j++){
	  if (jet_algorithm_list[i_j] == csi->type_parton[i_a][i_p]){flag++; break;}
	}
	if (flag){
	  runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);
	}
      }
      if (esi.object_list_selection[i] == "ljet" && abs(csi->type_parton[0][i_p]) < 5){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "bjet" && abs(csi->type_parton[0][i_p]) == 5){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "bjet_b" && csi->type_parton[0][i_p] == 5){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "bjet_bx" && csi->type_parton[0][i_p] == -5){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "lep" && (abs(csi->type_parton[0][i_p]) == 11 ||
					       abs(csi->type_parton[0][i_p]) == 13 ||
					       abs(csi->type_parton[0][i_p]) == 15)){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "lm" && (csi->type_parton[0][i_p] == 11 ||
						 csi->type_parton[0][i_p] == 13 ||
						 csi->type_parton[0][i_p] == 15)){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "lp" && (csi->type_parton[0][i_p] == -11 ||
						 csi->type_parton[0][i_p] == -13 ||
						 csi->type_parton[0][i_p] == -15)){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "e" && abs(csi->type_parton[0][i_p]) == 11){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "mu" && abs(csi->type_parton[0][i_p]) == 13){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "tau" && abs(csi->type_parton[0][i_p]) == 15){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "em" && csi->type_parton[0][i_p] == 11){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "mum" && csi->type_parton[0][i_p] == 13){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "taum" && csi->type_parton[0][i_p] == 15){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "ep" && csi->type_parton[0][i_p] == -11){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "mup" && csi->type_parton[0][i_p] == -13){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "taup" && csi->type_parton[0][i_p] == -15){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}

      if (esi.object_list_selection[i] == "nua" && (abs(csi->type_parton[0][i_p]) == 12 ||
					       abs(csi->type_parton[0][i_p]) == 14 ||
					       abs(csi->type_parton[0][i_p]) == 16)){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "nu" && (csi->type_parton[0][i_p] == 12 ||
						 csi->type_parton[0][i_p] == 14 ||
						 csi->type_parton[0][i_p] == 16)){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "nux" && (csi->type_parton[0][i_p] == -12 ||
						 csi->type_parton[0][i_p] == -14 ||
						 csi->type_parton[0][i_p] == -16)){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "nea" && abs(csi->type_parton[0][i_p]) == 12){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "nma" && abs(csi->type_parton[0][i_p]) == 14){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "nta" && abs(csi->type_parton[0][i_p]) == 16){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "ne" && csi->type_parton[0][i_p] == 12){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "nm" && csi->type_parton[0][i_p] == 14){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "nt" && csi->type_parton[0][i_p] == 16){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "nex" && csi->type_parton[0][i_p] == -12){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "nmx" && csi->type_parton[0][i_p] == -14){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "ntx" && csi->type_parton[0][i_p] == -16){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "photon" && csi->type_parton[0][i_p] == 22){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "tjet" && abs(csi->type_parton[0][i_p]) == 6){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "top" && csi->type_parton[0][i_p] == 6){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "atop" && csi->type_parton[0][i_p] == -6){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "wm" && csi->type_parton[0][i_p] == 24){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "wp" && csi->type_parton[0][i_p] == -24){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "z" && csi->type_parton[0][i_p] == 23){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
      if (esi.object_list_selection[i] == "h" && csi->type_parton[0][i_p] == 25){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}

      if (esi.object_list_selection[i] == "missing" && (abs(csi->type_parton[0][i_p]) == 12 ||
						   abs(csi->type_parton[0][i_p]) == 14 ||
						   abs(csi->type_parton[0][i_p]) == 16)){runtime_order[i].push_back(i_p); runtime_order_inverse[i_p].push_back(i);}
    }
  }

  output_determine_runtime_object();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::output_determine_runtime_object(){
  static  Logger logger("observable_set::output_determine_runtime_object");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger.newLine(LOG_INFO);
  logger << LOG_INFO << "Relevant objects:" << endl;
  logger << LOG_INFO << setw(5) << "no" << setw(10) << "object" << setw(10) << "n_min" << setw(10) << "n_max" << setw(10) << "pT" << setw(10) << "ET" << setw(10) << "|eta|"  << setw(10) << "|y|"<< endl;
  for (int i = 0; i < esi.object_list_selection.size(); i++){
    logger << LOG_INFO << setw(5) << i << setw(10) << esi.object_list_selection[i] << setw(10) << esi.pds[i].n_observed_min << setw(10) << esi.pds[i].n_observed_max << setw(10) << esi.pds[i].define_pT << setw(10) << esi.pds[i].define_ET << setw(10) << esi.pds[i].define_eta << setw(10) << esi.pds[i].define_y << endl;
  }
  logger.newLine(LOG_INFO);

  logger << LOG_INFO << "Relevant objects:   runtime_jet_recombination   collects partons that enter jet recombination." << endl;
  for (int i = 0; i < runtime_jet_recombination.size(); i++){
    logger << LOG_INFO << "runtime_jet_recombination[" << setw(3) << i << "] = " << setw(3) << runtime_jet_recombination[i] << "[" << setw(10) <<  left << esi.object_list_selection[runtime_jet_recombination[i]] << "]" << endl;
  }
  logger.newLine(LOG_INFO);
  
  logger << LOG_INFO << "Relevant objects:   runtime_photon_isolation   collects partons that could become isolated photons a la Frixione." << endl;
  for (int i = 0; i < runtime_photon_isolation.size(); i++){
    logger << LOG_INFO << "runtime_photon_isolation[" << setw(3) << i << "] = " << setw(3) << runtime_photon_isolation[i] << "[" << setw(10) <<   left << esi.object_list_selection[runtime_photon_isolation[i]] << "]" << endl;
  }
  logger.newLine(LOG_INFO);
  
  /*
  //  runtime_photon_recombination is never used !!!
  logger << LOG_INFO << "Relevant objects:   runtime_photon_recombination   collects partons that enter photon recombination." << endl;
  for (int i = 0; i < runtime_photon_recombination.size(); i++){
    logger << LOG_INFO << "runtime_photon_recombination[" << setw(3) << i << "] = " << setw(3) << runtime_photon_recombination[i] << "[" << setw(10) << left << esi.object_list_selection[runtime_photon_recombination[i]] << "]" << endl;
  }
  logger.newLine(LOG_INFO);
  */
  /*
  logger << LOG_INFO << "Relevant objects:   runtime_missing   collects partons that enter event selection with their parton-level momenta." << endl;
  for (int i = 0; i < runtime_missing.size(); i++){
    logger << LOG_INFO << "runtime_missing[" << setw(3) << i << "] = " << setw(3) << runtime_missing[i] << "[" << setw(10) << left << esi.object_list_selection[runtime_missing[i]] << "]" << endl;
  }
  logger.newLine(LOG_INFO);
  */

  logger << LOG_INFO << "Relevant objects:   runtime_original   collects partons that enter event selection with their parton-level momenta." << endl;
  for (int i = 0; i < runtime_original.size(); i++){
    logger << LOG_INFO << "runtime_original[" << setw(3) << i << "] = " << setw(3) << runtime_original[i] << "[" << setw(10) << left << esi.object_list_selection[runtime_original[i]] << "]" << endl;
  }
  logger.newLine(LOG_INFO);

  logger << LOG_INFO << "Relevant objects:   runtime_order   collects partons that belong to the respective object class." << endl;

  for (int i = 0; i < esi.object_list_selection.size(); i++){
    stringstream sstemp;
    for (int j = 0; j < runtime_order[i].size(); j++){sstemp << setw(3) << runtime_order[i][j] << "[" << setw(6) << csi->type_parton[0][runtime_order[i][j]] << "]";}
    logger << LOG_INFO << setw(10) << esi.object_list_selection[i] << " -> " << sstemp.str() << endl;
  }
  logger.newLine(LOG_INFO);

  logger << LOG_INFO << "Relevant objects:   runtime_order_inverse   collects object classes a parton belongs to." << endl;

  for (int i = 3; i < csi->type_parton[0].size(); i++){
    stringstream sstemp;
    for (int j = 0; j < runtime_order_inverse[i].size(); j++){sstemp << setw(10) << esi.object_list_selection[runtime_order_inverse[i][j]] << " [" << setw(3) << runtime_order_inverse[i][j] << "]";}
    logger << LOG_INFO << "pa[" << setw(3) << i << "] = " << setw(5) << csi->type_parton[0][i] << " -> " << sstemp.str() << endl;
  }
  logger.newLine(LOG_INFO);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
