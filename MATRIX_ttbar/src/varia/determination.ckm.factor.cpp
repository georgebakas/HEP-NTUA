#include "../include/classes.cxx"
//#include "../include/definitions.cxx"
#include "../include/definitions.observable.set.cxx"
double_complex determine_CKM_quarkchain(vector<int> type_parton, observable_set & oset){
  // all (anti-)quarks are defined as incoming!
  // A useful result CKM factor is only expected if there is only one flavour-changing coupling to a W boson involved!
  Logger logger("determine_CKM_quarkchain");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  vector<int> list_q_up;
  vector<int> list_a_down;
  vector<int> list_a_up;
  vector<int> list_q_down;
  for (int i_p = 1; i_p < type_parton.size(); i_p++){
    if (abs(type_parton[i_p]) == 1 || abs(type_parton[i_p]) == 3 || abs(type_parton[i_p]) == 5){
      if      (i_p < 3 && type_parton[i_p] > 0){list_q_down.push_back(type_parton[i_p]);}
      else if (i_p > 3 && type_parton[i_p] > 0){list_a_down.push_back(-type_parton[i_p]);}
      else if (i_p > 3 && type_parton[i_p] < 0){list_q_down.push_back(-type_parton[i_p]);}
      else if (i_p < 3 && type_parton[i_p] < 0){list_a_down.push_back(type_parton[i_p]);}
    }
    else if (abs(type_parton[i_p]) == 2 || abs(type_parton[i_p]) == 4 || abs(type_parton[i_p]) == 6){
      if      (i_p < 3 && type_parton[i_p] > 0){list_q_up.push_back(type_parton[i_p]);}
      else if (i_p > 3 && type_parton[i_p] > 0){list_a_up.push_back(-type_parton[i_p]);}
      else if (i_p > 3 && type_parton[i_p] < 0){list_q_up.push_back(-type_parton[i_p]);}
      else if (i_p < 3 && type_parton[i_p] < 0){list_a_up.push_back(type_parton[i_p]);}
    }
  }
  // quark-number conservation
  logger << LOG_DEBUG_VERBOSE << "list_q_up.size() = " << list_q_up.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "list_q_down.size() = " << list_q_down.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "list_a_up.size() = " << list_a_up.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "list_a_down.size() = " << list_a_down.size() << endl;
  assert((list_q_up.size() + list_q_down.size()) - (list_a_up.size() + list_a_down.size()) == 0);
  for (int i_q = list_q_up.size() - 1; i_q >= 0; i_q--){
    for (int i_a = 0; i_a < list_a_up.size(); i_a++){
      if (list_q_up[i_q] == -list_a_up[i_a]){
	list_q_up.erase(list_q_up.begin() + i_q, list_q_up.begin() + i_q + 1);
	list_a_up.erase(list_a_up.begin() + i_a, list_a_up.begin() + i_a + 1);
	break;
      }
    }
  }
  for (int i_q = list_q_down.size() - 1; i_q >= 0; i_q--){
    for (int i_a = 0; i_a < list_a_down.size(); i_a++){
      if (list_q_down[i_q] == -list_a_down[i_a]){
	list_q_down.erase(list_q_down.begin() + i_q, list_q_down.begin() + i_q + 1);
	list_a_down.erase(list_a_down.begin() + i_a, list_a_down.begin() + i_a + 1);
	break;
      }
    }
  }
  logger << LOG_DEBUG_VERBOSE << "list_q_up.size() = " << list_q_up.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "list_q_down.size() = " << list_q_down.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "list_a_up.size() = " << list_a_up.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "list_a_down.size() = " << list_a_down.size() << endl;

  assert(list_q_up.size() == list_a_down.size());
  assert(list_q_down.size() == list_a_up.size());

  double_complex factor_CKM = 1.;
  if (list_q_up.size() == 0 && list_a_down.size() == 0){}
  else if (list_q_up.size() == 1 && list_a_down.size() == 1){
    factor_CKM *= osi_msi.V_ckm[list_q_up[0] / 2][(-list_a_down[0] + 1) / 2];
    logger << LOG_DEBUG_VERBOSE << "CKM factor: " << list_q_up[0] << ", " << -list_a_down[0] << " = " << osi_msi.V_ckm[list_q_up[0] / 2][(-list_a_down[0] + 1) / 2] << endl;
  }
  else if (list_q_up.size() == list_a_down.size()){logger << LOG_INFO << "Two many flavour-changing up-down chains!" << endl;}
  else {}

  if (list_q_down.size() == 0 && list_a_up.size() == 0){}
  else if (list_q_down.size() == 1 && list_a_up.size() == 1){
    logger << LOG_DEBUG_VERBOSE << "-list_a_up[0] / 2 = " << -list_a_up[0] / 2 << endl;
    logger << LOG_DEBUG_VERBOSE << "(list_q_down[0] + 1) / 2 = " << (list_q_down[0] + 1) / 2 << endl;
    factor_CKM *= osi_msi.V_ckm[-list_a_up[0] / 2][(list_q_down[0] + 1) / 2];
    logger << LOG_DEBUG_VERBOSE << "CKM factor: " << -list_a_up[0] << ", " << list_q_down[0] << " = " << osi_msi.V_ckm[-list_a_up[0] / 2][(list_q_down[0] + 1) / 2] << endl;
  }
  else if (list_q_down.size() == list_a_up.size()){logger << LOG_INFO << "Two many flavour-changing down-up chains!" << endl;}
  else {}

  logger << LOG_DEBUG_VERBOSE << "CKM factor = " << factor_CKM << endl;
  if (norm(factor_CKM) == 0){cout << "vanishing CKM factor, exiting" << endl; exit(0);}
  return factor_CKM;
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
