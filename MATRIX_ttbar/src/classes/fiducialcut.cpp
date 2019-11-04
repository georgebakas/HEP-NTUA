#include "../include/classes.cxx"

////////////////////
//  constructors  //
////////////////////

fiducialcut::fiducialcut(string & _name, event_set & _esi, observable_set & _osi){
  Logger logger("fiducialcut::fiducialcut");
  logger << LOG_DEBUG << "started" << endl;

  name = _name;
  esi = &_esi;
  osi = &_osi;

  cut_value_lower = -std::numeric_limits<double>::infinity();
  cut_value_upper = +std::numeric_limits<double>::infinity();

  vector<string> split_name(1, "");
  //  split_name.push_back("");
  logger << LOG_DEBUG << "split_name.size() = " << split_name.size() << endl;
  int counter = 0;
  int start = 0;
  for (int i_c = 0; i_c < name.size(); i_c++){
    if (name[i_c] == ' '){start = 0;}
    if (name[i_c] != ' ' && start == 0){
      counter++;
      split_name.push_back("");
      start = 1;
    }
    if (name[i_c] != ' '){
      split_name[counter].push_back(name[i_c]);
    }
  }
  for (int i_s = 0; i_s < split_name.size(); i_s++){
    logger << LOG_DEBUG << "split_name[" << i_s << "] = " << split_name[i_s] << endl;
  }
  
  vector<int> identification(split_name.size(), 0);
  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "name = " << name << endl;
  logger.newLine(LOG_DEBUG);

  for (int i_s = 0; i_s < split_name.size(); i_s++){
    if (split_name[i_s] == ""){identification[i_s] = -1; continue;}

    //    logger << LOG_DEBUG_VERBOSE << "split_name[" << i_s << "] = " << split_name[i_s] << endl;

    if (identification[i_s - 1] == -1 || identification[i_s - 1] == 4){
    for (int i_l = 0; i_l < esi->fiducial_cut_list.size(); i_l++){
      //      logger << LOG_DEBUG << "split_name[" << i_s << "] == esi->fiducial_cut_list[" << i_l << "]" << endl;
      logger << LOG_DEBUG << "split_name[" << i_s << "] = " << split_name[i_s] << "  ?=  esi->fiducial_cut_list[" << i_l << "] = " << esi->fiducial_cut_list[i_l] << endl;
      if (split_name[i_s] == esi->fiducial_cut_list[i_l]){
	logger << LOG_DEBUG_VERBOSE << "split_name[" << i_s << "] = " << split_name[i_s] << " = esi->fiducial_cut_list[" << i_l << "] = " << esi->fiducial_cut_list[i_l] << endl;
	identification[i_s] = 1;
	break;
      }
    }
    }

    if (identification[i_s]){continue;}

    if (identification[i_s - 1] == 1 || identification[i_s - 1] == 2 || identification[i_s - 1] == 3){
      for (int i_l = 0; i_l < esi->object_list.size(); i_l++){
	if (split_name[i_s] == esi->object_list[i_l]){
	  //	logger << LOG_DEBUG_VERBOSE << "split_name[" << i_s << "] = " << split_name[i_s] << " = esi->object_list[" << i_l << "] = " << esi->object_list[i_l] << endl;
	  identification[i_s] = 2;
	  break;
	}
      }
    }
    
    if (identification[i_s]){continue;}
    
    if (split_name[i_s] == ">" || split_name[i_s] == "<"){
      identification[i_s] = 4;
    }
    
    if (identification[i_s]){continue;}
    
    if (identification[i_s - 1] == 2){
      if ((split_name[i_s][0] == '(' && split_name[i_s][split_name[i_s].size() - 1] == ')') ||
	  (split_name[i_s][0] == '[' && split_name[i_s][split_name[i_s].size() - 1] == ']') ||
	  (split_name[i_s][0] == '{' && split_name[i_s][split_name[i_s].size() - 1] == '}')){
	split_name[i_s].erase(split_name[i_s].size() - 1, 1);
	split_name[i_s].erase(0, 1);
      }
      else {
      }
      identification[i_s] = 3;
    }

    if (identification[i_s]){continue;}
    
    int check_number = 0;
    for (int i_c = 0; i_c < split_name[i_s].size(); i_c++){
      if (split_name[i_s][i_c] == '+' || 
	  split_name[i_s][i_c] == '-' || 
	  split_name[i_s][i_c] == '.' || 
	  split_name[i_s][i_c] == 'e' || 
	  split_name[i_s][i_c] == '0' || 
	  split_name[i_s][i_c] == '1' || 
	  split_name[i_s][i_c] == '2' || 
	  split_name[i_s][i_c] == '3' || 
	  split_name[i_s][i_c] == '4' || 
	  split_name[i_s][i_c] == '5' || 
	  split_name[i_s][i_c] == '6' || 
	  split_name[i_s][i_c] == '7' || 
	  split_name[i_s][i_c] == '8' || 
	  split_name[i_s][i_c] == '9'){ 
	check_number++;
      }
    }
    if (split_name[i_s].size() == check_number){
      identification[i_s] = 5;
    }
  }


  for (int i_s = 0; i_s < identification.size(); i_s++){
    logger << LOG_DEBUG << "split_name[" << i_s << "] = " << setw(10) << split_name[i_s] << "   identification[" << i_s << "] = " << identification[i_s] << endl;
    if (identification[i_s] == 0){
      logger << LOG_FATAL << "Input cannot be interpreted: " << split_name[i_s] << endl;
      exit(1);
    }
  }

  
  int no_cut = -1;
  // sanity checks and assignment from input to actual cut to be applied:
  for (int i_s = 0; i_s < identification.size(); i_s++){
    if (identification[i_s] == 1 && no_cut == -1){
      no_cut = i_s;
      observable = split_name[i_s];
    }
    else if (identification[i_s] == 1 && no_cut > -1){
      logger << LOG_FATAL << "Two cuts defined -> input cannot be interpreted." << endl;
      logger << LOG_FATAL << "Will be used for alternative cuts (like rapidity gaps)." << endl;
      exit(1);
    }
  }

  //  logger << LOG_DEBUG_VERBOSE << "Defined cut type identified." << endl;

  int no_lower = -1;
  int no_upper = -1;

  for (int i_s = 0; i_s < identification.size(); i_s++){
    if (identification[i_s] == 4){
      if ((i_s < no_cut && split_name[i_s] == "<") ||
	  (i_s > no_cut && split_name[i_s] == ">")){
	if (no_lower == -1){
	  no_lower = i_s;
	  if (i_s < no_cut && i_s > 0){
	    if (identification[i_s - 1] == 5){
	      cut_value_lower = atof(split_name[i_s - 1].c_str());
	    }
	  }
	  if (i_s > no_cut && i_s < identification.size() - 1){
	    if (identification[i_s + 1] == 5){
	      cut_value_lower = atof(split_name[i_s + 1].c_str());
	    }
	  }
	}
	else {
	  logger << LOG_FATAL << "Two limits defined -> input cannot be interpreted." << endl;
	  exit(1);
	}
      }

      if ((i_s < no_cut && split_name[i_s] == ">") ||
	  (i_s > no_cut && split_name[i_s] == "<")){
	if (no_upper == -1){
	  no_upper = i_s;
	  if (i_s < no_cut && i_s > 0){
	    if (identification[i_s - 1] == 5){
	      cut_value_upper = atof(split_name[i_s - 1].c_str());
	    }
	  }
	  if (i_s > no_cut && i_s < identification.size() - 1){
	    if (identification[i_s + 1] == 5){
	      cut_value_upper = atof(split_name[i_s + 1].c_str());
	    }
	  }
	}
	else {
	  logger << LOG_FATAL << "Two limits defined -> input cannot be interpreted." << endl;
	  exit(1);
	}
      }
    }
  }

  if (cut_value_lower == -std::numeric_limits<double>::infinity() && cut_value_upper == +std::numeric_limits<double>::infinity()){
    logger << LOG_FATAL << "No cut specified." << endl;
    exit(1);
  }
  else if (cut_value_lower > -std::numeric_limits<double>::infinity() && cut_value_upper == +std::numeric_limits<double>::infinity()){
    type = LOWER;
  }
  else if (cut_value_lower == -std::numeric_limits<double>::infinity() && cut_value_upper < +std::numeric_limits<double>::infinity()){
    type = UPPER;
  }
  else if (cut_value_lower < cut_value_upper){
    type = WINDOW;
  }
  else if (cut_value_lower > cut_value_upper){
    type = GAP;
  }
  else if (cut_value_lower == cut_value_upper){
    logger << LOG_FATAL << "No cut interval specified (lower limit equals upper limit)." << endl;
    exit(1);
  }
  

  
  //  logger << LOG_DEBUG_VERBOSE << "Upper and lower limits identified." << endl;

  for (int i_s = 0; i_s < identification.size(); i_s++){
    if (identification[i_s] == 2){
      int temp = 0;
      type_particle.push_back(split_name[i_s]);
      if (i_s < identification.size() - 1){
	if (identification[i_s + 1] == 3){
	  order_min.push_back(atoi(split_name[i_s + 1].c_str()));
	  order_max.push_back(atoi(split_name[i_s + 1].c_str()));
	  temp = 1;
	}
      }
      if (!temp){
	order_min.push_back(1);
	///	order_min.push_back(0);
	//	cout << "i_s = " << i_s << "   temp = " << temp << endl;
	//	cout << "type_particle[" << type_particle.size() - 1 << "]  = " << type_particle[type_particle.size() - 1] << endl;
	//	cout << "esi->observed_object.size() = " << esi->observed_object.size() << endl;
	order_max.push_back(esi->pda[esi->observed_object[type_particle[type_particle.size() - 1]]].n_observed_max);
	///	order_max.push_back(esi->pda[esi->observed_object[type_particle[type_particle.size() - 1]]].n_observed_max - 1);
	//	order_max.push_back(2);
      }
    }
  }

  //  logger << LOG_DEBUG_VERBOSE << "Involved particle types and specification identified." << endl;


  vector<int> n_particles_in_type(type_particle.size());

  logger << LOG_DEBUG_VERBOSE << "type_particle.size() = " << type_particle.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "order_min.size() = " << order_min.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "order_max.size() = " << order_max.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "n_particles_in_type.size() = " << n_particles_in_type.size() << endl;
  int n_combination_max = 1.;
  for (int i_p = 0; i_p < type_particle.size(); i_p++){
    n_particles_in_type[i_p] = (order_max[i_p] - order_min[i_p] + 1);
    n_combination_max *= (order_max[i_p] - order_min[i_p] + 1);
  }

  logger << LOG_DEBUG << "n_combination_max = " << n_combination_max << endl;
  //  logger << LOG_DEBUG << "n_combination_max = " << n_combination_max << endl;
  // e.g.:
  // 1 -> 3 entries
  // 2 -> 2 entries
  // 3 -> 2 entries

  // 0 - 0 - 0
  // 0 - 0 - 1
  // 0 - 1 - 0
  // 0 - 1 - 1
  // 1 - 0 - 0
  // 1 - 0 - 1
  // 1 - 1 - 0
  // 1 - 1 - 1
  // 2 - 0 - 0
  // 2 - 0 - 1
  // 2 - 1 - 0
  // 2 - 1 - 1

  for (int i_c = 0; i_c < n_combination_max; i_c++){
    vector<int> temp_entry(type_particle.size());
    int rest = i_c;
    //    for (int i_p = type_particle.size() - 1; i_p >= 0; i_p--){
    for (int i_p = 0; i_p < type_particle.size(); i_p++){
      temp_entry[i_p] = rest % n_particles_in_type[i_p] + 1;
      ///      temp_entry[i_p] = rest % n_particles_in_type[i_p];
      logger << LOG_DEBUG << "i_c = " << i_c << "   rest = " << rest << "   rest - temp_entry[" << i_p << "] = " << rest - temp_entry[i_p] << "   new_rest = " << (rest - temp_entry[i_p]) / n_particles_in_type[i_p] << "   n_particles_in_type[" << i_p << "] = " << n_particles_in_type[i_p] << endl;
      rest = (rest - temp_entry[i_p] + 1) / n_particles_in_type[i_p];
      ///      rest = (rest - temp_entry[i_p]) / n_particles_in_type[i_p];
      //      logger << LOG_DEBUG << "temp_entry[" << i_p << "] = " <<  temp_entry[i_p] << endl;
      //      logger << LOG_DEBUG << "rest = " << rest << endl;
    }
    // n_particles_in_type[i_p] == 1 -> no container, i.e. '==' instead of '>=' in removal procedure:
    // maybe switch meaning of temp_entry here: selected number instead of ordering number (which would always be zero in that case):

    for (int i_p = 0; i_p < type_particle.size(); i_p++){
      if (n_particles_in_type[i_p] == 1){temp_entry[i_p] = order_min[i_p];}
    }

    int remove_combination = 0;
    for (int i_p = 0; i_p < type_particle.size(); i_p++){
      for (int j_p = i_p + 1; j_p < type_particle.size(); j_p++){
	//	logger << LOG_DEBUG << "type_particle[" << i_p << "] = " << type_particle[i_p] << "   type_particle[" << j_p << "] = " << type_particle[j_p] << "     " << (type_particle[i_p] == type_particle[j_p]) << endl; 
	logger << LOG_DEBUG << "i_c = " << i_c << "   temp_entry[" << i_p << "] = " << temp_entry[i_p] << "   temp_entry[" << j_p << "] = " << temp_entry[j_p] << "     " << (temp_entry[i_p] >= temp_entry[j_p]) << endl; 
	if (type_particle[i_p] == type_particle[j_p] && 
	    ((n_particles_in_type[i_p] > 1 && n_particles_in_type[j_p] > 1 && temp_entry[i_p] >= temp_entry[j_p]) ||
	     (n_particles_in_type[i_p] == 1 && n_particles_in_type[j_p] > 1 && temp_entry[i_p] == temp_entry[j_p]) ||
	     (n_particles_in_type[i_p] > 1 && n_particles_in_type[j_p] == 1 && temp_entry[i_p] == temp_entry[j_p]))){
	  remove_combination = 1;
	  //	  logger << LOG_DEBUG << "i_p = " << i_p << "   j_p = " << j_p << "   remove_combination = " << remove_combination << endl;
	  break;
	}
	//	logger << LOG_DEBUG << "i_p = " << i_p << "   j_p = " << j_p << endl; //"   remove_combination = " << remove_combination << endl;
      }
    }

    if (!(remove_combination)){
      for (int i_p = 0; i_p < type_particle.size(); i_p++){
	logger << LOG_DEBUG << "temp_entry[" << i_p << "] = " <<  temp_entry[i_p] << endl;
      }

      // doesn't work so far, if (at least some) particle numbers are given !!!
      // version of temp_entry needed that contains not 0, but the specified ordering number for each group !!!

      // Reduce particle numbers by one to simplify container access:
      for (int i_p = 0; i_p < type_particle.size(); i_p++){temp_entry[i_p]--;}
      
      no_particles.push_back(temp_entry);
    }
  }

  n_combination = no_particles.size();
  cut_combination.resize(n_combination, 0);

  logger.newLine(LOG_DEBUG);

  // Validity checks for cuts should be implemented !!!  

  for (int i_s = 0; i_s < split_name.size(); i_s++){
    logger << LOG_DEBUG << "split_name[" << i_s << "] = " << setw(10) << split_name[i_s] << "   identification = " << identification[i_s] << endl;
  }
  logger << LOG_DEBUG << "cut_type = " << split_name[no_cut] << " @ " << no_cut << endl;
  logger << LOG_DEBUG << "cut_value_lower = " << cut_value_lower << " @ " << no_lower << endl;
  logger << LOG_DEBUG << "cut_value_upper = " << cut_value_upper << " @ " << no_upper << endl;

  logger << LOG_DEBUG << "n_combination = " << no_particles.size() << ":" << endl;
  /*
  stringstream temp_particles_ss;
  for (int i_p = 0; i_p < type_particle.size(); i_p++){
    temp_particles_ss << type_particle[i_p] << "   ";
  }
  temp_particles_ss << "   ";
  */
  for (int i_c = 0; i_c < no_particles.size(); i_c++){
    stringstream temp_no_particles_ss;
    for (int i_p = 0; i_p < type_particle.size(); i_p++){
      temp_no_particles_ss << type_particle[i_p] << " " << no_particles[i_c][i_p] << "   ";
    }
    logger << LOG_DEBUG << "particles = " << temp_no_particles_ss.str() << endl;
  }

  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG << "finished" << endl;
}



void fiducialcut::apply_fiducialcut(int i_a){
  Logger logger("fiducialcut::apply_fiducialcut");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  cut = 0;
  cut_combination.resize(n_combination, 0);

  for (int i_c = 0; i_c < no_particles.size(); i_c++){
    logger << LOG_DEBUG_VERBOSE << "observable [" << i_c << "]" << endl;
    // calculate observable:
    
    double value_observable = 0.;
    if (observable == "M"){
      fourvector temp_fv;
      for (int i_p = 0; i_p < type_particle.size(); i_p++){
	temp_fv = temp_fv + osi->particle_event[osi->access_object[type_particle[i_p]]][i_a][no_particles[i_c][i_p]].momentum;
      }
      value_observable = temp_fv.m();
    }
    else if (observable == "pT"){
      if (type_particle.size() == 1){
	value_observable = osi->particle_event[osi->access_object[type_particle[0]]][i_a][no_particles[i_c][0]].pT;
      }
      else {
	fourvector temp_fv;
	for (int i_p = 0; i_p < type_particle.size(); i_p++){
	  temp_fv = temp_fv + osi->particle_event[osi->access_object[type_particle[i_p]]][i_a][no_particles[i_c][i_p]].momentum;
	}
	value_observable = temp_fv.pT();
      }
    }
    else if (observable == "|eta|"){
      if (type_particle.size() == 1){
	value_observable = osi->particle_event[osi->access_object[type_particle[0]]][i_a][no_particles[i_c][0]].eta;
      }
      else {
	fourvector temp_fv;
	for (int i_p = 0; i_p < type_particle.size(); i_p++){
	  temp_fv = temp_fv + osi->particle_event[osi->access_object[type_particle[i_p]]][i_a][no_particles[i_c][i_p]].momentum;
	}
	value_observable = abs(temp_fv.eta());
      }
    }
    else if (observable == "|y|"){
      if (type_particle.size() == 1){
	value_observable = osi->particle_event[osi->access_object[type_particle[0]]][i_a][no_particles[i_c][0]].rapidity;
      }
      else {
	fourvector temp_fv;
	for (int i_p = 0; i_p < type_particle.size(); i_p++){
	  temp_fv = temp_fv + osi->particle_event[osi->access_object[type_particle[i_p]]][i_a][no_particles[i_c][i_p]].momentum;
	}
	value_observable = abs(temp_fv.rapidity());
      }
    }
    else if (observable == "deta"){
      if (type_particle.size() == 2){
	value_observable = osi->particle_event[osi->access_object[type_particle[0]]][i_a][no_particles[i_c][0]].eta - osi->particle_event[osi->access_object[type_particle[1]]][i_a][no_particles[i_c][1]].eta;

      }
      else {
	logger << LOG_FATAL << "Fiducial cut for " << observable << " with " << type_particle.size() << " is not defined." << endl;
	exit(1);
      }
    }
    else {
      logger << LOG_FATAL << "Fiducial cut for " << observable << " is not defined." << endl;
      exit(1);
    }
    logger << LOG_DEBUG_VERBOSE << "value_observable = " << value_observable << endl;

    // apply cut on observable:
    
    if ((type == UPPER && (value_observable > cut_value_upper)) ||
	(type == LOWER && (value_observable < cut_value_lower)) ||
	(type == WINDOW && (value_observable < cut_value_lower || value_observable > cut_value_upper)) ||
	(type == GAP && (value_observable > cut_value_upper && value_observable > cut_value_lower))){

      // Only in some optimization phase ???
      //      cut_combination[i_c]++;
      logger << LOG_DEBUG << "i_a = " << setw(2) << i_a << "   Cut applied due to " << name << " [" << i_c << "]: value = " << value_observable << endl;
      osi->cut_ps[i_a] = -1;
    }
    else if (value_observable != value_observable){
      logger << LOG_DEBUG << "i_a = " << setw(2) << i_a << "   Cut applied due to " << name << " [" << i_c << "]: value = " << value_observable << endl;
      osi->cut_ps[i_a] = -1;
    }
    else {
      logger << LOG_DEBUG << "i_a = " << setw(2) << i_a << "   Cut passed due to " << name << " [" << i_c << "]: value = " << value_observable << endl;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




ostream & operator << (ostream &s, const fiducialcut & fc){
  //  if (fv.crossing == -1){s << "crossed ";}
  string unit = "";
  if (fc.observable == "M" || fc.observable == "pT"){unit = " GeV";}
        
  //  vector<string> full_observable(fc.no_particles.size());
  for (int i_c = 0; i_c < fc.no_particles.size(); i_c++){
    //    s << setw(30) << fc.name << "   " << setw(6) << fc.type << "   ";
    stringstream full_observable_ss;
    full_observable_ss << fc.observable << " [ ";
    //    stringstream temp_no_particles_ss;
    for (int i_p = 0; i_p < fc.type_particle.size(); i_p++){
      full_observable_ss << fc.type_particle[i_p] << " " << fc.no_particles[i_c][i_p] + 1;
      if (i_p < fc.type_particle.size() - 1){full_observable_ss << "   ";}
    }
    full_observable_ss << " ]";
    
    if (fc.type == UPPER){
      s << "UPPER    " << full_observable_ss.str() << " < " << fc.cut_value_upper << unit << endl;
    }
    if (fc.type == LOWER){
      s << "LOWER    " << full_observable_ss.str() << " > " << fc.cut_value_lower << unit << endl;
    }
    if (fc.type == WINDOW){
      s << "WINDOW   " << fc.cut_value_lower << unit << " < " << full_observable_ss.str() << " < " << fc.cut_value_upper << unit << endl;
    }
    if (fc.type == GAP){
      s << "GAP      " << full_observable_ss.str() << " < " << fc.cut_value_upper << unit << "   OR   " << full_observable_ss.str() << " > " << fc.cut_value_lower << unit << endl;
    }
  }
  /*
  logger << LOG_DEBUG << "cut_type = " << split_name[no_cut] << " @ " << no_cut << endl;

  s << "( " << setw(21) << setprecision(15) << fv.t << " ; " << setw(21) << setprecision(15) << fv.x << " , " << setw(21) << setprecision(15) << fv.y << " , " << setw(21) << setprecision(15) << fv.z <<" )";
  //  s << "( " << setw(14) << setprecision(8) << fv.t << " ; " << setw(14) << setprecision(8) << fv.x << " , " << setw(14) << setprecision(8) << fv.y << " , " << setw(14) << setprecision(8) << fv.z <<" )";
  */
  return s;
}
