#include "../include/classes.cxx"

////////////////////
//  constructors  //
////////////////////

//event_set::event_set(){}
event_set::event_set(){
  Logger logger("event_set::event_set");
  logger << LOG_DEBUG << "started" << endl;

  define_basic_object_list();
  pda.resize(object_list.size());
  n_parton_nu = 0;

  fiducial_cut_list.clear();
  fiducial_cut_list.push_back("M");
  fiducial_cut_list.push_back("pT");
  fiducial_cut_list.push_back("|eta|");
  fiducial_cut_list.push_back("|y|");
  fiducial_cut_list.push_back("deta");
  fiducial_cut_list.push_back("dy");
  fiducial_cut_list.push_back("absdeta");
  fiducial_cut_list.push_back("absdy");
  fiducial_cut_list.push_back("dReta");
  fiducial_cut_list.push_back("dRy");
  fiducial_cut_list.push_back("dphi");

  logger << LOG_DEBUG << "finished" << endl;
}



// This function is probably not needed at all !!!
void event_set::determine_n_partonlevel(vector<vector<int> > & type_parton){
  Logger logger("event_set::determine_n_partonlevel");
  logger << LOG_DEBUG << "started" << endl;

  for (int i_p = 3; i_p < type_parton[0].size(); i_p++){
    if (abs(type_parton[0][i_p]) == 22){pda[observed_object["photon"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 11 || abs(type_parton[0][i_p]) == 13 || abs(type_parton[0][i_p]) == 15){pda[observed_object["lep"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 11 || type_parton[0][i_p] == 13 || type_parton[0][i_p] == 15){pda[observed_object["lm"]].n_partonlevel++;}
    if (type_parton[0][i_p] == -11 || type_parton[0][i_p] == -13 || type_parton[0][i_p] == -15){pda[observed_object["lp"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 11){pda[observed_object["e"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 13){pda[observed_object["mu"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 15){pda[observed_object["tau"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 11){pda[observed_object["em"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 13){pda[observed_object["mum"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 15){pda[observed_object["taum"]].n_partonlevel++;}
    if (type_parton[0][i_p] == -11){pda[observed_object["ep"]].n_partonlevel++;}
    if (type_parton[0][i_p] == -13){pda[observed_object["mup"]].n_partonlevel++;}
    if (type_parton[0][i_p] == -15){pda[observed_object["taup"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 12 || abs(type_parton[0][i_p]) == 14 || abs(type_parton[0][i_p]) == 16){pda[observed_object["nua"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 12 || type_parton[0][i_p] == 14 || type_parton[0][i_p] == 16){pda[observed_object["nu"]].n_partonlevel++;}
    if (type_parton[0][i_p] == -12 || type_parton[0][i_p] == -14 || type_parton[0][i_p] == -16){pda[observed_object["nux"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 12){pda[observed_object["nea"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 14){pda[observed_object["nma"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 16){pda[observed_object["nta"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 12){pda[observed_object["ne"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 14){pda[observed_object["nm"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 16){pda[observed_object["nt"]].n_partonlevel++;}
    if (type_parton[0][i_p] == -12){pda[observed_object["nex"]].n_partonlevel++;}
    if (type_parton[0][i_p] == -14){pda[observed_object["nmx"]].n_partonlevel++;}
    if (type_parton[0][i_p] == -16){pda[observed_object["ntx"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 12 || abs(type_parton[0][i_p]) == 14 || abs(type_parton[0][i_p]) == 16){pda[observed_object["missing"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 12 || abs(type_parton[0][i_p]) == 14 || abs(type_parton[0][i_p]) == 16){n_parton_nu++;}
    if (abs(type_parton[0][i_p]) < 5){pda[observed_object["ljet"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 5){pda[observed_object["bjet_b"]].n_partonlevel++;}
    if (type_parton[0][i_p] == -5){pda[observed_object["bjet_bx"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 5){pda[observed_object["bjet"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 6){pda[observed_object["top"]].n_partonlevel++;}
    if (type_parton[0][i_p] == -6){pda[observed_object["atop"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) == 6){pda[observed_object["tjet"]].n_partonlevel++;}
    if (abs(type_parton[0][i_p]) < 6){pda[observed_object["jet"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 24){pda[observed_object["wm"]].n_partonlevel++;}
    if (type_parton[0][i_p] == -24){pda[observed_object["wp"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 23){pda[observed_object["z"]].n_partonlevel++;}
    if (type_parton[0][i_p] == 25){pda[observed_object["h"]].n_partonlevel++;}
  }

  logger << LOG_DEBUG << "finished" << endl;
}



void event_set::define_basic_object_list(){
  //  0 -> jet
  //  1 -> charged particle (without jets; lepton, W)
  //  2 -> photon 
  //  3 -> neutral particle (neutrino, Higgs, Z)
  // -1 -> user-defined particle (no automatic modifications)

  object_list.push_back("");
  object_category.push_back(0);
  object_list.push_back("jet");
  object_category.push_back(0);
  object_list.push_back("ljet"); // order l and b changed
  object_category.push_back(0);
  object_list.push_back("bjet");
  object_category.push_back(0);
  
  object_list.push_back("bjet_b");
  object_category.push_back(0);
  object_list.push_back("bjet_bx");
  object_category.push_back(0);
  
  object_list.push_back("lep");
  object_category.push_back(1);
  object_list.push_back("lm");
  object_category.push_back(1);
  object_list.push_back("lp");
  object_category.push_back(1);
  object_list.push_back("e");
  object_category.push_back(1);
  object_list.push_back("mu");
  object_category.push_back(1);
  object_list.push_back("tau");
  object_category.push_back(1);
  object_list.push_back("em");
  object_category.push_back(1);
  object_list.push_back("mum");
  object_category.push_back(1);
  object_list.push_back("taum");
  object_category.push_back(1);
  object_list.push_back("ep");
  object_category.push_back(1);
  object_list.push_back("mup");
  object_category.push_back(1);
  object_list.push_back("taup");
  object_category.push_back(1);
  object_list.push_back("photon");
  object_category.push_back(2);
  object_list.push_back("tjet");
  object_category.push_back(0);
  object_list.push_back("top");
  object_category.push_back(0);
  object_list.push_back("atop");
  object_category.push_back(0);
  object_list.push_back("wp");
  object_category.push_back(1);
  object_list.push_back("wm");
  object_category.push_back(1);
  object_list.push_back("z");
  object_category.push_back(1);
  object_list.push_back("h");
  object_category.push_back(3);
  
  object_list.push_back("nua");
  object_category.push_back(4);
  object_list.push_back("nu");
  object_category.push_back(4);
  object_list.push_back("nux");
  object_category.push_back(4);
  object_list.push_back("nea");
  object_category.push_back(4);
  object_list.push_back("nma");
  object_category.push_back(4);
  object_list.push_back("nta");
  object_category.push_back(4);
  object_list.push_back("ne");
  object_category.push_back(4);
  object_list.push_back("nm");
  object_category.push_back(4);
  object_list.push_back("nt");
  object_category.push_back(4);
  object_list.push_back("nex");
  object_category.push_back(4);
  object_list.push_back("nmx");
  object_category.push_back(4);
  object_list.push_back("ntx");
  object_category.push_back(4);
 
  object_list.push_back("none");
  object_category.push_back(0);
  
  object_list.push_back("missing");
  object_category.push_back(0);
    
  for (int i = 0; i < object_list.size(); i++){observed_object[object_list[i]] = i;}
}


void event_set::define_specific_object_list(user_defined & user){
  Logger logger("event_set::define_specific_object_list");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  int old_n = object_list.size();

  for (int i_p = 0; i_p < user.particle_name.size(); i_p++){
    logger << LOG_DEBUG << "user.particle_name[" << i_p << "] = " << user.particle_name[i_p] << endl;
    object_list.push_back(user.particle_name[i_p]);
    object_category.push_back(-1);
  }

  pda.resize(object_list.size());

  for (int i = old_n; i < object_list.size(); i++){
    observed_object[object_list[i]] = i;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
