#include "../include/classes.cxx"

void observable_set::initialization_distribution(){
  Logger logger("observable_set::initialization_distribution");
  logger << LOG_DEBUG << "started" << endl;

// **************************************************************************
// *                                                                        *
// *  reading parameters from 'file_distribution.dat'                       *
// *                                                                        *
// **************************************************************************
  logger << LOG_DEBUG << "switch_distribution = " << switch_distribution << endl;
  if (switch_distribution){readin_file_distribution();}
  logger << LOG_INFO << "List of differential distributions:" << endl << endl;
  for (int i_d = 0; i_d < dat.size(); i_d++){logger << LOG_INFO << "No. " << setw(3) << i_d << ":" << dat[i_d] << endl;}
  logger.newLine(LOG_INFO);


// **************************************************************************
// *                                                                        *
// *  reading parameters from 'file_dddistribution.dat'                     *
// *                                                                        *
// **************************************************************************
  if (switch_distribution){readin_file_dddistribution();}
  logger << LOG_INFO << "List of doubly-differential distributions:" << endl << endl;
  for (int i_ddd = 0; i_ddd < dddat.size(); i_ddd++){logger << LOG_INFO << "dddistribution " << setw(3) << i_ddd << ":" << dddat[i_ddd] << endl;}
  logger.newLine(LOG_INFO);



  int isystem = 0; 
  string xorder;
  string temp_directory;

  string dir_distribution = "distribution";
  if (switch_output_distribution){system_execute(logger, "mkdir " + dir_distribution, isystem);}
  filename_distribution.resize(dat.size());
  filename_dddistribution.resize(dddat.size());
  for (int i = 0; i < dat.size(); i++){
    filename_distribution[i] = dir_distribution + "/" + dat[i].xdistribution_name + "_" + name_process + ".dat";
  }
  for (int i = 0; i < dddat.size(); i++){
    filename_dddistribution[i] = dir_distribution + "/" + dddat[i].name + "_" + name_process + ".dat";
  }
  
  if (switch_CV){
    if (switch_output_distribution){system_execute(logger, "mkdir " + dir_distribution + "/CV", isystem);}
    filename_distribution_CV.resize(n_scales_CV, filename_distribution);
    filename_dddistribution_CV.resize(n_scales_CV, filename_dddistribution);

    logger << LOG_DEBUG_VERBOSE << "directory_name_scale_CV.size() = " << directory_name_scale_CV.size() << endl;
    for (int i_s = 0; i_s < n_scales_CV; i_s++){
      logger << LOG_DEBUG_VERBOSE << "directory_name_scale_CV[" << i_s << "] = " << directory_name_scale_CV[i_s] << endl;
      if (switch_output_distribution){system_execute(logger, "mkdir " + dir_distribution + "/CV/" + directory_name_scale_CV[i_s], isystem);}
      logger << LOG_DEBUG_VERBOSE << xorder << "   execution status: " << isystem << "." << endl;
    
      filename_distribution_all_CV = dir_distribution + "/CV" + "/distribution_" + name_process + ".txt";
      // Should replace the following files at some point...

      for (int i_d = 0; i_d < dat.size(); i_d++){
	filename_distribution_CV[i_s][i_d]  = dir_distribution + "/CV/" + directory_name_scale_CV[i_s] + "/" + dat[i_d].xdistribution_name + "_" + name_process + ".dat";
      }
      for (int i_d = 0; i_d < dddat.size(); i_d++){
	filename_dddistribution_CV[i_s][i_d]  = dir_distribution + "/CV/" + directory_name_scale_CV[i_s] + "/" + dddat[i_d].name + "_" + name_process + ".dat";
      }
    }
  }

  if (switch_TSV){
    logger << LOG_DEBUG_VERBOSE << "before bin_count_TSV" << endl;
    bin_count_TSV.resize(n_qTcut, vector<vector<double> > (dat.size() + dddat.size()));
    for (int i_q = 0; i_q < n_qTcut; i_q++){
      for (int i_d = 0; i_d < dat.size(); i_d++){
	bin_count_TSV[i_q][i_d].resize(dat[i_d].n_bins, 0);
      }
      for (int i_d = 0; i_d < dddat.size(); i_d++){
	int i_ddd = dat.size() + i_d;
	bin_count_TSV[i_q][i_ddd].resize(dddat[i_d].n_bins, 0);
      }
    }

    logger << LOG_DEBUG_VERBOSE << "before bin_weight_TSV" << endl;

    // n_qTcut will typically be by far too large
    // define a selection of a handful of qTcut values instead and use their number n_qTcut_distribution (or so...) !!!

    bin_weight_TSV.resize(n_extended_set_TSV, vector<vector<vector<vector<vector<double> > > > > (n_qTcut));
    bin_weight2_TSV.resize(n_extended_set_TSV, vector<vector<vector<vector<vector<double> > > > > (n_qTcut));
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      for (int i_q = 0; i_q < n_qTcut; i_q++){
	bin_weight_TSV[i_s][i_q].resize(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (dat.size() + dddat.size())));
	bin_weight2_TSV[i_s][i_q].resize(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (dat.size() + dddat.size())));
	for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	    for (int i_d = 0; i_d < dat.size(); i_d++){
	      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d].resize(dat[i_d].n_bins, 0.);
	      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d].resize(dat[i_d].n_bins, 0.);
	    }
	    for (int i_d = 0; i_d < dddat.size(); i_d++){
	      int i_ddd = dat.size() + i_d;
	      bin_weight_TSV[i_s][i_q][i_r][i_f][i_ddd].resize(dddat[i_d].n_bins, 0.);
	      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_ddd].resize(dddat[i_d].n_bins, 0.);
	    }
	  }
	}
      }
    }
    
    logger << LOG_DEBUG_VERBOSE << "after bin_weight_TSV" << endl;

    filename_distribution_TSV.resize(n_extended_set_TSV);
    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      string temp_directory = dir_distribution + "/" + name_extended_set_TSV[i_s];
      if (switch_output_distribution){system_execute(logger, "mkdir " + temp_directory, isystem);}
      filename_distribution_TSV[i_s] = temp_directory + "/distribution_" + name_process + ".txt";
    }

    if (switch_distribution_at_all_TSV){
      change_weight_TSV.resize(n_extended_set_TSV);
      change_weight2_TSV.resize(n_extended_set_TSV);
      change_weight_2nd_TSV.resize(n_extended_set_TSV);
      change_weight2_2nd_TSV.resize(n_extended_set_TSV);

      change_qTcut_bin_weight_TSV.resize(n_extended_set_TSV);
      change_qTcut_bin_weight2_TSV.resize(n_extended_set_TSV);
      change_qTcut_bin_weight_2nd_TSV.resize(n_extended_set_TSV);
      change_qTcut_bin_weight2_2nd_TSV.resize(n_extended_set_TSV);

      for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	if (switch_distribution_TSV[i_s]){
	  change_weight_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.));
	  change_weight2_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.));
	  change_weight_2nd_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.));
	  change_weight2_2nd_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<double> (n_scale_fact_TSV[i_s], 0.));

	  change_qTcut_bin_weight_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (0)));
	  change_qTcut_bin_weight2_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (0)));
	  change_qTcut_bin_weight_2nd_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (0)));
	  change_qTcut_bin_weight2_2nd_TSV[i_s].resize(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (0)));
	}
      }
    }

  }




  logger << LOG_DEBUG << "finished" << endl;
}

void observable_set::readin_file_distribution(){
  Logger logger("observable_set::readin_file_distribution");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

// **************************************************************************
// *                                                                        *
// *  reading parameters from 'file_distribution.dat'                       *
// *                                                                        *
// **************************************************************************

  char LineBuffer[512];
  string s0;
  //  int counter = 0;
  string readindata;
  vector<string> readin, vs0;
  vector<string> user_variable;
  vector<string> user_variable_counter;
  vector<string> user_value;
  string pathname = "" + path_to_main + "/file_distribution.dat";
  logger << LOG_DEBUG << "path_to_main = " << path_to_main << endl;
  logger << LOG_DEBUG << "pathname = " << pathname << endl;
  ifstream in_new_distribution_gen(pathname.c_str());
  while (in_new_distribution_gen.getline(LineBuffer, 512)){readin.push_back(LineBuffer);}
  in_new_distribution_gen.close();
  
  if (readin.size() == 0) {
    logger << LOG_INFO << "file_distributions.dat could not be read in" << endl;
    //    assert(readin.size() > 0);
  }

  vector<string> user_variable_additional;
  get_userinput_extravalue_from_readin(user_variable, user_variable_additional, user_value, readin);

  vector<string> distributionname;
  vector<string> distributiontype;
  
  vector<vector<vector<int> > > all_particle;
  vector<vector<vector<int> > > required_subparticle;
  vector<vector<vector<int> > > all_particle_group;
  vector<vector<vector<int> > > all_particle_order;

  vector<double> startpoint;
  vector<double> endpoint;
  vector<double> binwidth;
  vector<string> edges;
  vector<int> binnumber;
  vector<string> binningtype;
  int current_counter = 0;
  vector<int> vi0(0);
  vector<vector<int> > vvi0(0);
  vector<vector<int> > vvi3(3);
  current_counter = -1;

  for (int i = 0; i < user_variable.size(); i++){
    if (user_variable[i] == "distributionname"){
      current_counter++;
      distributionname.push_back(s0);
      distributiontype.push_back(s0);
      all_particle.push_back(vvi0);
      required_subparticle.push_back(vvi0);
      all_particle_group.push_back(vvi0);
      all_particle_order.push_back(vvi0);
      startpoint.push_back(0.);
      endpoint.push_back(0.);
      binwidth.push_back(0.);
      edges.push_back("");
      binnumber.push_back(0);
      binningtype.push_back("linear");

      distributionname[current_counter] = user_value[i];
    }
    else if (user_variable[i] == "distributiontype"){distributiontype[current_counter] = user_value[i];}
    else if (user_variable[i] == "startpoint"){startpoint[current_counter] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "endpoint"){endpoint[current_counter] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "binnumber"){binnumber[current_counter] = atoi(user_value[i].c_str());}
    else if (user_variable[i] == "binwidth"){binwidth[current_counter] = atof(user_value[i].c_str());}
    else if (user_variable[i] == "edges"){edges[current_counter] = user_value[i];}
    else if (user_variable[i] == "binningtype"){binningtype[current_counter] = user_value[i].c_str();}
    else if (user_variable[i] == "particle"){
      //      logger << LOG_DEBUG << "new: distribution " << current_counter << endl;
      //      logger << LOG_DEBUG << "user_variable[" << setw(3) << i << "] = " << user_variable[i] << "   user_variable_additional[" << setw(3) << i << "] = " << user_variable_additional[i] << "   user_value[" << setw(3) << i << "] = " << user_value[i] << endl;
      int temp_i = 1;
      int temp_r = 2; // former default: 1 !!!
      if (user_variable_additional[i] != ""){
	vector<string> temp_s(1);
	for (int i_s = 0; i_s < user_variable_additional[i].size(); i_s++){
	  if (user_variable_additional[i][i_s] == ' ' || user_variable_additional[i][i_s] == char(9)){
	    if (temp_s[temp_s.size() - 1].size() == 0){}
	    else {temp_s.push_back("");}
	  }
	  else {temp_s[temp_s.size() - 1].push_back(user_variable_additional[i][i_s]);}
	}	
	if (temp_s.size() > 0){
	  temp_i = atoi(temp_s[0].c_str());
	  if (temp_s.size() == 2){
	    if (temp_s[1] == "+"){temp_r = 1;}
	    else if (temp_s[1] == "x"){temp_r = 2;}
	    else if (temp_s[1] == "0"){temp_r = 0;}
	    else if (temp_s[1] == "-"){temp_r = -1;}
	    else {logger << LOG_ERROR<< "entry '" << temp_s[1] << "' not allowed!" << endl; assert(false); exit(1);}
	  }
	}
	// before introduction of '+' and '-'
	//	temp_i = atoi(user_variable_additional[i].c_str());
      }
      if (temp_i > all_particle[current_counter].size()){
	all_particle[current_counter].resize(temp_i);
	all_particle_group[current_counter].resize(temp_i);
	required_subparticle[current_counter].resize(temp_i);
      }
      required_subparticle[current_counter][temp_i - 1].push_back(temp_r);

      vector<string> temp_s(1);
      for (int i_s = 0; i_s < user_value[i].size(); i_s++){
	if (user_value[i][i_s] == ' ' || user_value[i][i_s] == char(9)){
	  if (temp_s[temp_s.size() - 1].size() == 0){}
	  else{temp_s.push_back("");}
	}
	else {temp_s[temp_s.size() - 1].push_back(user_value[i][i_s]);}
      }
      //      logger << LOG_DEBUG << "current_counter = " << current_counter << endl;
      //      logger << LOG_DEBUG << "temp_i = " << temp_i << endl;
      //      logger << LOG_DEBUG << "temp_s.size() = " << temp_s.size() << endl;
      //      for (int i_t = 0; i_t < temp_s.size(); i_t++){logger << LOG_DEBUG << "temp_s[" << i_t << "] = " << temp_s[i_t] << endl;}
      if (temp_s.size() == 1){
	all_particle_group[current_counter][temp_i - 1].push_back(0);
	all_particle[current_counter][temp_i - 1].push_back(atoi(temp_s[0].c_str()));
      }
      else if (temp_s.size() == 2){
	if (esi.observed_object[temp_s[0]] != 0 || temp_s[0] == "all"){
	  //	  all_particle_group[current_counter][temp_i - 1].push_back(esi.observed_object[temp_s[0]]);
	  all_particle_group[current_counter][temp_i - 1].push_back(esi.observed_object[equivalent_object[temp_s[0]]]);

	  //	  logger << LOG_DEBUG << "new: distribution no. " << current_counter << "   temp_s[0] = " << setw(7) << temp_s[0] << "   esi.observed_object[" << setw(7) << temp_s[0] << "] = " << setw(7) << esi.observed_object[temp_s[0]] << "   no_relevant_object[equivalent_object[" << setw(7) << temp_s[0] << "] = " << setw(7) << equivalent_object[temp_s[0]] << "] = " << no_relevant_object[equivalent_object[temp_s[0]]] << endl;

	  //	  logger << LOG_DEBUG << "new: distribution no. " << current_counter << "   temp_s[0] = " << setw(7) << temp_s[0] << "   esi.observed_object[" << setw(7) << temp_s[0] << "] = " << setw(7) << esi.observed_object[temp_s[0]] << "   esi.observed_object[equivalent_object[" << setw(7) << temp_s[0] << "] = " << setw(7) << equivalent_object[temp_s[0]] << "] = " << setw(7) << esi.observed_object[equivalent_object[temp_s[0]]] << endl;
	  //equivalent_object
	  //equivalent_object, map<int, int> & equivalent_no_object, map<string, int> & no_relevant_object
	  all_particle[current_counter][temp_i - 1].push_back(atoi(temp_s[1].c_str()));
	}
	else {
          logger << LOG_ERROR << "object not known: " << temp_s[0] << endl;
          assert(false);
	  exit (1);
	}
      }
      else {
        logger << LOG_ERROR << "Input not understood!" << endl;
        logger << LOG_ERROR << "user_variable[" << setw(3) << i << "] = " << user_variable[i] << "   user_variable_additional[" << setw(3) << i << "] = " << user_variable_additional[i] << "   user_value[" << setw(3) << i << "] = " << user_value[i] << endl;
        assert(false);
	exit(1);
      }
    }
  }
  
  for (int i_d = 0; i_d < current_counter + 1; i_d++){
    logger << LOG_DEBUG << "distributionname[" << i_d << "] = " << distributionname[i_d] << endl;
    xdistribution one_dist(distributionname[i_d], distributiontype[i_d], all_particle[i_d], all_particle_group[i_d], required_subparticle[i_d], startpoint[i_d], endpoint[i_d], binnumber[i_d], binwidth[i_d], edges[i_d], binningtype[i_d], this); 
    dat.push_back(one_dist);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


//vector<dddistribution> & all_dddistributions, vector<xdistribution> & all_distributions, observable_set & oset
void observable_set::readin_file_dddistribution(){
  Logger logger("observable_set::readin_file_dddistributions");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

// **************************************************************************
// *                                                                        *
// *  reading parameters from 'file_distribution.dat'                       *
// *                                                                        *
// **************************************************************************
  char LineBuffer[256];
  string s0;
  //  int counter = 0;
  string readindata;
  vector<string> readin, vs0;
  vector<string> user_variable;
  vector<string> user_variable_counter;
  vector<string> user_value;
  string pathname = "" + path_to_main + "/file_dddistribution.dat";
  logger << LOG_DEBUG << "pathname = " << pathname << endl;

  ifstream in_new_distribution_gen(pathname.c_str());
  while (in_new_distribution_gen.getline(LineBuffer, 256)) {
    readin.push_back(LineBuffer);
  }
  in_new_distribution_gen.close();

  for (int i = 0; i < readin.size(); i++){
    readindata = readin[i][0];
    if (readindata != "/" && readindata != "#" && readindata != "%"){
      int start = 0;
      user_variable.push_back(s0);
      user_variable_counter.push_back(s0);
      user_value.push_back(s0);
      
      for (int j = 0; j < readin[i].size(); j++){
	if (start == 0 || start == 1){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 0){}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    user_variable[user_variable.size() - 1].push_back(readin[i][j]);
	    if (start != 1){start = 1;}
	  }
	  else {start++;}
	}
	else if (start == 2){
	  if (readin[i][j] == '='){start++;}
	  else if ((readin[i][j] == ' ') || (readin[i][j] == char(9))){}
	  else {logger << LOG_DEBUG << "Incorrect input in line " << i << endl;
	    user_variable.erase(user_variable.end(), user_variable.end());
	    user_value.erase(user_value.end(), user_value.end());
	    break;
	  }
	}
	else if (start == 3 || start == 4){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 3){}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    if (readin[i][j] != char(13)){user_value[user_value.size() - 1].push_back(readin[i][j]);}
	    if (start != 4){start = 4;}
	  }
	  else {start++;}
	}
	else {break;}
      }
    }
  }


  vector<string> dddistributionname;
  vector<string> distributionname_1;
  vector<string> distributionname_2;
  vector<int> d_1;
  vector<int> d_2;
  vector<int> n_bins;
  int current_counter = 0;
  vector<int> vi0(0);
  current_counter = -1;
  for (int i = 0; i < user_variable.size(); i++){
      logger << LOG_DEBUG << "user_variable[" << i << "] = " << user_variable[i] << "   user_value[" << i << "] = " << user_value[i] << endl;
    if (user_variable[i] == "dddistributionname"){
      current_counter++;
      dddistributionname.push_back(s0);
      distributionname_1.push_back(s0);
      distributionname_2.push_back(s0);
      //      d_1.push_back(0);
      //      d_2.push_back(0);
      //      n_bins.push_back(0);
      dddistributionname[current_counter] = user_value[i];
    }
    else if (user_variable[i] == "distributionname_1"){distributionname_1[current_counter] = user_value[i];}
    else if (user_variable[i] == "distributionname_2"){distributionname_2[current_counter] = user_value[i];}
  }
    
  for (int d = 0; d < current_counter + 1; d++){
      logger << LOG_DEBUG << "dddistributionname[" << d << "] = " << dddistributionname[d] << endl;
    int d_1 = -1;
    int d_2 = -1;
    //    int n_bins;
    for (int id = 0; id < dat.size(); id++){
      if (dat[id].xdistribution_name == distributionname_1[d]){d_1 = id;}
      if (dat[id].xdistribution_name == distributionname_2[d]){d_2 = id;}
    }
    if (d_1 != -1 && d_2 != -1){
      //      n_bins = dat[d_1].n_bins * dat[d_2].n_bins;
      logger << LOG_DEBUG << "dddistributionname[" << d << "] = " << dddistributionname[d] << "   d_1 = " << d_1 << "   d_2 = " << d_2 << endl;
      dddistribution one_dddist(dddistributionname[d], d_1, d_2, this);
      //      dddistribution one_dddist(dddistributionname[d], d_1, d_2, dat, this);
      //      dddistribution one_dddist(dddistributionname[d], d_1, d_2, dat[d_1], dat[d_2]);
      dddat.push_back(one_dddist);
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_extended_distribution(){//, observable_set & oset
  Logger logger("determine_extended_distribution");
  logger << LOG_DEBUG << "called" << endl;
  //  static int initialization = 1;

  //  int plotmode = 1;
  //  int n_subgroup = subgroup.size();
  // !!! subgroup not needed !!!
  string vs;
  string s0;
  vector<string> vs0;
  string filename;
  vector<string> readin;
  //  char LineBuffer[128];

  //  int QCD_order = 0;

  //#include "../include/declarations.cxx"
  /*
  map<string, int> observed_object;
  vector<string> object_list;
  object_list.reserve(30);
  vector<int> object_category;
  object_category.reserve(30);
  //  vector<string> relevant_object_list;
  //  vector<int> relevant_object_category;

  temp_define_object_list(observed_object, object_list, object_category);

  assert(esi.observed_object == observed_object);
  assert(esi.observed_list == observed_list);
  assert(esi.observed_category == observed_category);


  //  map<string, string> equivalent_object;
  //  map<int, int> equivalent_no_object;
  //  map<string, int> no_relevant_object;
  */

// **************************************************************************
// *                                                                        *
// *  reading parameters from 'file_distribution.dat'                       *
// *                                                                        *
// **************************************************************************
//  vector<xdistribution> all_distribution;
//  simplified_readin_distributions(all_distribution, oset);
  readin_file_distribution();
  vector<xdistribution> all_distribution = dat;
  logger << LOG_DEBUG << "all_distribution.size() = " << all_distribution.size() << endl;
  for (int d = 0; d < all_distribution.size(); d++){logger << LOG_DEBUG << "distribution[" << d << "] = " << endl << all_distribution[d] << endl;}

// **************************************************************************
// *                                                                        *
// *  reading parameters from 'file_dddistribution.dat'                     *
// *                                                                        *
// **************************************************************************
//  vector<dddistribution> all_dddistribution;
//  simplified_readin_dddistributions(all_dddistribution, all_distribution, oset);
  readin_file_dddistribution();
  vector<dddistribution> all_dddistribution = dddat;

  logger << LOG_DEBUG << "all_dddistribution.size() = " << all_dddistribution.size() << endl;
  for (int d = 0; d < all_dddistribution.size(); d++){logger << LOG_DEBUG << "all_dddistribution[" << d << "] = " << endl << all_dddistribution[d] << endl;}
 
  vector<int> mapping_dddistribution(all_distribution.size(), 0);
  // should be available for all subroutines !!!

  //  vector<xdistribution> extended_distribution = all_distribution;
  extended_distribution = all_distribution;
  extended_distribution.reserve(all_distribution.size() + all_dddistribution.size() * 2);
  //  extended_distribution.reserve(all_distribution.size() + all_dddistribution.size());
    logger << LOG_DEBUG << "extended_distribution.size() = " << extended_distribution.size() << endl;

  for (int i_d = 0; i_d < all_dddistribution.size(); i_d++){
    vector<int> vi;
    string type;
    //    double asym_n_bins;
    if ((((all_dddistribution[i_d].distribution_1).xdistribution_type == "y" ||
	  (all_dddistribution[i_d].distribution_1).xdistribution_type == "absydiff" ||
	  (all_dddistribution[i_d].distribution_1).xdistribution_type == "eta" ||
	  (all_dddistribution[i_d].distribution_1).xdistribution_type == "absetadiff") &&
	 ((all_dddistribution[i_d].distribution_1).n_bins == 2)) &&
	(((all_dddistribution[i_d].distribution_2).xdistribution_type == "y" ||
	  (all_dddistribution[i_d].distribution_2).xdistribution_type == "absydiff" ||
	  (all_dddistribution[i_d].distribution_2).xdistribution_type == "eta" ||
	  (all_dddistribution[i_d].distribution_2).xdistribution_type == "absetadiff") &&
	 ((all_dddistribution[i_d].distribution_2).n_bins == 2))){
      type = "doubleasym";
    }
    else if ((((all_dddistribution[i_d].distribution_1).xdistribution_type == "y" ||
	       (all_dddistribution[i_d].distribution_1).xdistribution_type == "absydiff" ||
	       (all_dddistribution[i_d].distribution_1).xdistribution_type == "eta" ||
	       (all_dddistribution[i_d].distribution_1).xdistribution_type == "absetadiff")) &&
	     ((all_dddistribution[i_d].distribution_1).n_bins == 2)){
      xdistribution temp(all_dddistribution[i_d].name, "asym1", (all_dddistribution[i_d].distribution_2).all_particle, (all_dddistribution[i_d].distribution_2).all_particle_group, (all_dddistribution[i_d].distribution_2).required_subparticle, all_dddistribution[i_d].distribution_2.start, 0, all_dddistribution[i_d].distribution_2.n_bins, all_dddistribution[i_d].distribution_2.step, "", "linear");//, 0, 3);
      extended_distribution.push_back(temp);
      mapping_dddistribution.push_back(i_d);
      xdistribution tempsum("sum" + all_dddistribution[i_d].name, "sumasym1", (all_dddistribution[i_d].distribution_2).all_particle, (all_dddistribution[i_d].distribution_2).all_particle_group, (all_dddistribution[i_d].distribution_2).required_subparticle, all_dddistribution[i_d].distribution_2.start, 0, all_dddistribution[i_d].distribution_2.n_bins, all_dddistribution[i_d].distribution_2.step, "", "linear");//, 0, 3);
      extended_distribution.push_back(tempsum);
      mapping_dddistribution.push_back(i_d);
    }
    else if ((((all_dddistribution[i_d].distribution_2).xdistribution_type == "y" ||
	       (all_dddistribution[i_d].distribution_2).xdistribution_type == "absydiff" ||
	       (all_dddistribution[i_d].distribution_2).xdistribution_type == "eta" ||
	       (all_dddistribution[i_d].distribution_2).xdistribution_type == "absetadiff")) &&
	     ((all_dddistribution[i_d].distribution_2).n_bins == 2)){
      xdistribution temp(all_dddistribution[i_d].name, "asym2", (all_dddistribution[i_d].distribution_1).all_particle, (all_dddistribution[i_d].distribution_1).all_particle_group, (all_dddistribution[i_d].distribution_1).required_subparticle, all_dddistribution[i_d].distribution_1.start, 0, all_dddistribution[i_d].distribution_1.n_bins, all_dddistribution[i_d].distribution_1.step, "", "linear");//, 0, 3);
      extended_distribution.push_back(temp);
      mapping_dddistribution.push_back(i_d);
      xdistribution tempsum("sum" + all_dddistribution[i_d].name, "sumasym2", (all_dddistribution[i_d].distribution_1).all_particle, (all_dddistribution[i_d].distribution_1).all_particle_group, (all_dddistribution[i_d].distribution_1).required_subparticle, all_dddistribution[i_d].distribution_1.start, 0, all_dddistribution[i_d].distribution_1.n_bins, all_dddistribution[i_d].distribution_1.step, "", "linear");//, 0, 3);
      extended_distribution.push_back(tempsum);
      mapping_dddistribution.push_back(i_d);
    }
    else {
      type = "dddistribution";
      cout << "Seems like double-distributions have not been properly installed yet!" << endl;

      xdistribution temp = get_fake_distribution_from_dddistribution(all_dddistribution[i_d]);
      //(all_dddistribution[i_d].name, type, (all_dddistribution[i_d].distribution_1).all_particle, (all_dddistribution[i_d].distribution_1).all_particle_group, (all_dddistribution[i_d].distribution_1).required_subparticle, all_dddistribution[i_d].distribution_1.start, 0, all_dddistribution[i_d].distribution_1.n_bins * all_dddistribution[i_d].distribution_2.n_bins, all_dddistribution[i_d].distribution_1.step, "", all_dddistribution[i_d].distribution_1.type_binning, 0, 3);
      /*
      temp.xdistribution_name = all_dddistribution[i_d].name;
      temp.n_bins = all_dddistribution[i_d].n_bins;
      */
      extended_distribution.push_back(temp);
      mapping_dddistribution.push_back(i_d);

    }
  }

  cout << "extended_distribution.size() = " << extended_distribution.size() << endl;
  for (int d = 0; d < extended_distribution.size(); d++){cout << "extended_distribution[" << d << "] = " << endl << extended_distribution[d] << endl;}

  //  vector<double> fakeasymfactor(extended_distribution.size());
  fakeasymfactor.resize(extended_distribution.size(), 1.);
  for (int i_d = 0; i_d < extended_distribution.size(); i_d++){
    // everything below this line is repeated for all (extended) distributions !!!
    logger << LOG_DEBUG << "extended_distribution[" << i_d << "] = " << extended_distribution[i_d].xdistribution_name << endl;
    int ddd = mapping_dddistribution[i_d];
    //      if (i_d >= all_distribution.size()){ddd = i_d - all_distribution.size();}


    if (i_d >= all_distribution.size()){
      if (extended_distribution[i_d].xdistribution_type == "asym1"){fakeasymfactor[ddd] = (all_dddistribution[ddd].distribution_1).step;}
      else if (extended_distribution[i_d].xdistribution_type == "asym2"){fakeasymfactor[ddd] = (all_dddistribution[ddd].distribution_2).step;}
      if (extended_distribution[i_d].xdistribution_type == "sumasym1"){fakeasymfactor[ddd] = (all_dddistribution[ddd].distribution_1).step;}
      else if (extended_distribution[i_d].xdistribution_type == "sumasym2"){fakeasymfactor[ddd] = (all_dddistribution[ddd].distribution_2).step;}
    }
    else {fakeasymfactor[ddd] = extended_distribution[i_d].step;}
  }
}



void observable_set::get_userinput_extravalue_from_readin(vector<string> & user_variable, vector<string> & user_variable_additional, vector<string> & user_value, vector<string> & readin){
  Logger logger("get_userinput_extravalue_from_readin");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string readindata;
  for (int i = 0; i < readin.size(); i++){
    readindata = readin[i][0];
    if (readindata != "/" && readindata != "#" && readindata != "%"){
      int start = 0;
      user_variable.push_back("");
      user_variable_additional.push_back("");
      user_value.push_back("");
      for (int j = 0; j < readin[i].size(); j++){
	if (start == 0 || start == 1){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 0){}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    user_variable[user_variable.size() - 1].push_back(readin[i][j]);
	    if (start != 1){start = 1;}
	  }
	  else {start++;}
	}
	else if (start == 2){
	  if (readin[i][j] == '='){start = 5;}
	  else if ((readin[i][j] == ' ') || (readin[i][j] == char(9))){}
	  else if (start == 2){
	    start++;
	    j--;
	  }
	  else {
	    logger << LOG_ERROR << "Incorrect input in line " << i << endl;
	    user_variable.erase(user_variable.end(), user_variable.end());
	    user_variable_additional.erase(user_variable_additional.end(), user_variable_additional.end());
	    user_value.erase(user_value.end(), user_value.end());
	    break;
	  }
	}
	else if (start == 3){
	  if ((readin[i][j] != '=')){
	    user_variable_additional[user_variable_additional.size() - 1].push_back(readin[i][j]);
	  }
	  else {start++; j--;} // should be the same as shifting (start == 4) here !!!
	}
	else if (start == 4){
	  // additional: should be allowed to contain ' ' !!!
	  if (readin[i][j] == '='){
	    start = 5;
	    logger << LOG_DEBUG << "before: user_variable_aditional: ---" << user_variable_additional[user_variable_additional.size() - 1] << "---" << endl;
	    for (int i_s = user_variable_additional[user_variable_additional.size() - 1].size() - 1; i_s > 0; i_s--){
	      if (user_variable_additional[user_variable_additional.size() - 1][i_s] == ' ' ||
		  user_variable_additional[user_variable_additional.size() - 1][i_s] == char(9)){
		user_variable_additional[user_variable_additional.size() - 1].erase(user_variable_additional[user_variable_additional.size() - 1].end() - 1, user_variable_additional[user_variable_additional.size() - 1].end());
	      }
	      else {break;}
	    }
	    logger << LOG_DEBUG << "after:  user_variable_aditional: ---" << user_variable_additional[user_variable_additional.size() - 1] << "---" << endl;
	  }
	  //	  else if ((readin[i][j] == ' ') || (readin[i][j] == char(9))){}
	  else {
	    logger << LOG_ERROR << "Incorrect input in line " << i << endl;
	    user_variable.erase(user_variable.end(), user_variable.end());
	    user_variable_additional.erase(user_variable_additional.end(), user_variable_additional.end());
	    user_value.erase(user_value.end(), user_value.end());
	    break;
	  }
	}
	else if (start == 5 || start == 6){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 5){}
	  else if (readindata != "/" && readindata != "#" && readindata != "%"){
	    //	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    user_value[user_value.size() - 1].push_back(readin[i][j]);
	    if (start != 6){start = 6;}
	  }
	  else {start++;}
	}
	else {break;}
      }
      logger << LOG_DEBUG << "i = " << setw(3) << i << "   ---" << user_value[user_value.size() - 1] << "---" << endl;
      for (int i_s = user_value[user_value.size() - 1].size() - 1; i_s > 0; i_s--){
	if (user_value[user_value.size() - 1][i_s] == ' ' ||
	    user_value[user_value.size() - 1][i_s] == char(9)){
	  user_value[user_value.size() - 1].erase(user_value[user_value.size() - 1].end() - 1, user_value[user_value.size() - 1].end());
	}
	else {break;}
      }
      logger << LOG_DEBUG << "i = " << setw(3) << i << "   ---" << user_value[user_value.size() - 1] << "---" << endl;
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



// !!! not used: !!!
void observable_set::distribution_back_to_back_configuration(int i_a, int k_g){
  Logger logger("observable_set::distribution_back_to_back_configuration");
  logger << LOG_DEBUG << "started" << endl;

  if (abs(particle_event[k_g][i_a][0].pT - particle_event[k_g][i_a][1].pT) / (particle_event[k_g][i_a][0].pT + particle_event[k_g][i_a][1].pT) < 1.e-12){
    if (particle_event[k_g][0].size() > 1){
      if (particle_event[k_g][i_a][0].momentum.x1() > 0. && particle_event[k_g][0][0].momentum.x1() < 0.){
	swap(particle_event[k_g][i_a][0], particle_event[k_g][i_a][1]);
	logger << LOG_DEBUG_VERBOSE << "Momenta have been swapped in order to match real phase space." << endl;
      } 
    }
  }

  logger << LOG_DEBUG << "finished" << endl;
}



// !!! Was only partially used from here on: !!!  It actually is used - for Born and RA !!!
// !!! Functions from  routines.distribution.cpp  are used. !!!
// !!! -> need to be transferred here...

void observable_set::determine_distribution_complete(){
  Logger logger("observable_set::determine_distribution_complete");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  // done already when calling this function
  //  if (osi_switch_distribution == 0){return;}

  // replaces determine_distribution_bin(osi_dat, osi_bin, osi_bin_max, oset);
  for (int i_d = 0; i_d < dat.size(); i_d++){dat[i_d].determineBin();}
  logger << LOG_DEBUG_VERBOSE << "single-differential distributions:   determine_distribution_bin done !" << endl;

  // to be introduced (now part of determine_double_distribution):
  
  for (int i_ddd = 0; i_ddd < dddat.size(); i_ddd++){dddat[i_ddd].determineBin();}
  logger << LOG_DEBUG_VERBOSE << "double-differential distributions:   determine_distribution_bin done !" << endl;
  
  
  // Should be replaced since only dat[i_d].bin[i_a] should be used !!!
  /*
  for (int i_d = 0; i_d < osi_dat.size(); i_d++){
    for (int i_a = 0; i_a < osi_n_ps; i_a++){
      // use only new implementation: osi_bin[i_d][i_a] etc. could be generically replaced by osi_dat[i_d].bin[i_a] !!!
      osi_bin[i_d][i_a] = osi_dat[i_d].bin[i_a];
      if (osi_dat[i_d].typeCumulative != CUMULATIVE_NONE){osi_bin_max[i_d][i_a] = osi_dat[i_d].bin_max[i_a];}
    }
  }
  */
  
  determine_single_distribution();
  logger << LOG_DEBUG_VERBOSE << "determine_single_distribution  done !" << endl;

  determine_double_distribution();
  logger << LOG_DEBUG_VERBOSE << "determine_dddistribution  done !" << endl;

  logger << LOG_DEBUG_VERBOSE << "n_extended_set_TSV = " << n_extended_set_TSV << "   switch_distribution_at_all_TSV = " << switch_distribution_at_all_TSV << endl;
  
  if (n_extended_set_TSV && switch_distribution_at_all_TSV){
    preparation_distribution_TSV();

    determine_single_distribution_TSV();
    determine_double_distribution_TSV();
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::determine_single_distribution(){
  static Logger logger("observable_set::determine_single_distribution");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  vector<vector<int> > bin_phasespace(3);

  for (int i_d = 0; i_d < dat.size(); i_d++){
    if (dat[i_d].typeCumulative != CUMULATIVE_NONE){
      // cumulative distributions
      // only non-symmetric (symm == 0) distributions possible here !!!
      if (n_ps == 1){
	if (dat[i_d].bin[0] != -1 && dat[i_d].bin_max[0] != -1){
	  for (int presentbin = dat[i_d].bin[0]; presentbin < dat[i_d].bin_max[0]; presentbin++){
	    bin_weight[i_d][presentbin] += integrand_D[1][0];
	    if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][presentbin] += integrand_D_CV[i_s][1][0];}}
	    bin_counts[i_d][presentbin]++;

	    if (csi->type_parton[0][1] != csi->type_parton[0][2]){
	      // asymmetric distribution
	      // check if always correct, e.g. for uc->... !!!
	      bin_weight[i_d][presentbin] += integrand_D[2][0];
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][presentbin] += integrand_D_CV[i_s][2][0];}}
	      bin_counts[i_d][presentbin]++;
	    }
	  
	    if (integrand_D[2][0] == 0.){
	      // should be the same as integrand_D_CV[i_s][1][0] == 0.
	      bin_weight2[i_d][presentbin] += pow(integrand_D[1][0], 2);
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight2_CV[i_s][i_d][presentbin] += pow(integrand_D_CV[i_s][1][0], 2);}}
	    }
	    else {
	      bin_weight2[i_d][presentbin] += pow(integrand_D[0][0], 2);
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight2_CV[i_s][i_d][presentbin] += pow(integrand_D_CV[i_s][0][0], 2);}}
	    }
	  }
	}
      }
      else {
	vector<vector<double> > xbin_dev(dat[i_d].n_bins);
	vector<vector<vector<double> > > xbin_dev_CV(n_scales_CV, vector<vector<double> > (dat[i_d].n_bins));
	for (int i_a = 0; i_a < n_ps; i_a++){
	  //	  cout << "bin[" << setw(3) << d << "][" << setw(2) << i_a << "] = " << dat[i_d].bin[i_a] << "   bin_max[" << setw(3) << d << "][" << setw(2) << i_a << "] = " <<  dat[i_d].bin_max[i_a] << endl;
	  if (dat[i_d].bin[i_a] != -1 && dat[i_d].bin_max[i_a] != -1){
	    for (int presentbin = dat[i_d].bin[i_a]; presentbin < dat[i_d].bin_max[i_a]; presentbin++){
	      bin_weight[i_d][presentbin] += integrand_D[1][i_a];
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][presentbin] += integrand_D_CV[i_s][1][i_a];}}
	      bin_counts[i_d][presentbin]++;

	      if (csi->type_parton[0][1] != csi->type_parton[0][2]){  // asymmetric distribution
		bin_weight[i_d][presentbin] += integrand_D[2][i_a];
		if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][presentbin] += integrand_D_CV[i_s][2][i_a];}}
		bin_counts[i_d][presentbin]++;	      
	      }

	      if (integrand_D[2][i_a] == 0.){  // should be the same as integrand_D_CV[i_s][1][0] == 0.
		xbin_dev[presentbin].push_back(integrand_D[1][i_a]);
		if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_CV[i_s][presentbin].push_back(integrand_D_CV[i_s][1][i_a]);}}
	      }
	      else {
		xbin_dev[presentbin].push_back(integrand_D[0][i_a]);
		if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_CV[i_s][presentbin].push_back(integrand_D_CV[i_s][0][i_a]);}}
	      }
	    }
	  }
	}
	for (int x = 0; x < bin_weight2[i_d].size(); x++){
	  if (xbin_dev[x].size() != 0){bin_weight2[i_d][x] += pow(accumulate(xbin_dev[x].begin(), xbin_dev[x].end(), 0.), 2);}
	  if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){if (xbin_dev_CV[i_s][x].size() != 0){bin_weight2_CV[i_s][i_d][x] += pow(accumulate(xbin_dev_CV[i_s][x].begin(), xbin_dev_CV[i_s][x].end(), 0.), 2);}}}
	}
      }
    }


    else if (dat[i_d].typeCumulative == CUMULATIVE_NONE){  // non-cumulative distributions
      if (n_ps == 1){  // only one phase-space
	int presentbin = dat[i_d].bin[0];
	int mirrorbin = dat[i_d].n_bins - 1 - presentbin;
	if (presentbin != -1){
	  bin_weight[i_d][presentbin] += integrand_D[1][0];
	  if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][presentbin] += integrand_D_CV[i_s][1][0];}}
	  bin_counts[i_d][presentbin]++;

	  if (csi->type_parton[0][1] != csi->type_parton[0][2]){
	    if (dat[i_d].symm == 0){ // asymmetric distribution
	      bin_weight[i_d][presentbin] += integrand_D[2][0];
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][presentbin] += integrand_D_CV[i_s][2][0];}}
	      bin_counts[i_d][presentbin]++;
	    }
	    else if (dat[i_d].symm == 1) { // symmetric distribution
	      bin_weight[i_d][mirrorbin] += integrand_D[2][0];
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][mirrorbin] += integrand_D_CV[i_s][2][0];}}
	      bin_counts[i_d][mirrorbin]++;
	    }
	  }
	  
	  if (integrand_D[2][0] == 0.){  // should be the same as integrand_D_CV[i_s][1][0] == 0.
	    bin_weight2[i_d][presentbin] += pow(integrand_D[1][0], 2);
	    if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight2_CV[i_s][i_d][presentbin] += pow(integrand_D_CV[i_s][1][0], 2);}}
	  }
	  else {
	    if (dat[i_d].symm == 0){  // asymmetric distribution
	      bin_weight2[i_d][presentbin] += pow(integrand_D[0][0], 2);
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight2_CV[i_s][i_d][presentbin] += pow(integrand_D_CV[i_s][0][0], 2);}}
	    }
	    else if (dat[i_d].symm == 1){  // symmetric distribution
	      bin_weight2[i_d][presentbin] += pow(integrand_D[1][0], 2);
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight2_CV[i_s][i_d][presentbin] += pow(integrand_D_CV[i_s][1][0], 2);}}
	      bin_weight2[i_d][mirrorbin] += pow(integrand_D[2][0], 2);
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight2_CV[i_s][i_d][mirrorbin] += pow(integrand_D_CV[i_s][2][0], 2);}}
	    }
	  }
	}
      }
      else {  // more than one phase-space (n_ps > 1)
	vector<int> xbin_count(0);
	vector<double> xbin_dev(0);
	vector<vector<double> > xbin_dev_CV(n_scales_CV), xbin_dev_2_CV(n_scales_CV);
	for (int i_a = 0; i_a < n_ps; i_a++){
	  int presentbin = dat[i_d].bin[i_a];
	  int mirrorbin = dat[i_d].n_bins - 1 - presentbin;
	  if (presentbin != -1){
	    bin_weight[i_d][presentbin] += integrand_D[1][i_a];
	    if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][presentbin] += integrand_D_CV[i_s][1][i_a];}}
	    bin_counts[i_d][presentbin]++;
	    if (csi->type_parton[0][1] != csi->type_parton[0][2]){
	      if (dat[i_d].symm == 0){ // asymmetric distribution
		bin_weight[i_d][presentbin] += integrand_D[2][i_a];
		if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][presentbin] += integrand_D_CV[i_s][2][i_a];}}
		bin_counts[i_d][presentbin]++;
	      }
	      else if (dat[i_d].symm == 1){ // symmetric distribution
		bin_weight[i_d][mirrorbin] += (integrand_D[2][i_a]);
		if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][mirrorbin] += integrand_D_CV[i_s][2][i_a];}}
		bin_counts[i_d][mirrorbin]++;
	      }
	    }
	    //
	    // deviation
	    //
	    int exist = -1;
	    for (int x = 0; x < xbin_count.size(); x++){
	      if (presentbin == xbin_count[x]){exist = x; break;}
	    }
	    if (exist == -1){
	      xbin_count.push_back(presentbin);
	      xbin_dev.push_back(integrand_D[1][i_a]);
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_CV[i_s].push_back(integrand_D_CV[i_s][1][i_a]);}}
	    }
	    else {
	      xbin_dev[exist] += integrand_D[1][i_a];
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_CV[i_s][exist] += integrand_D_CV[i_s][1][i_a];}}
	    }
	    if (integrand_D[2][i_a] != 0.){  // should be the same as integrand_D_CV[i_s][1][0] != 0.
	      exist = -1;
	      if (dat[i_d].symm == 0){ // asymmetric distribution
		for (int x = 0; x < xbin_count.size(); x++){
		  if (presentbin == xbin_count[x]){exist = x; break;}
		}
		if (exist == -1){
		  xbin_count.push_back(presentbin);
		  xbin_dev.push_back(integrand_D[2][i_a]);
		  if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_CV[i_s].push_back(integrand_D_CV[i_s][2][i_a]);}}
		}
		else {
		  xbin_dev[exist] += integrand_D[2][i_a];
		  if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_CV[i_s][exist] += integrand_D_CV[i_s][2][i_a];}}
		}
	      }
	      else if (dat[i_d].symm == 1){ // symmetric distribution
		for (int x = 0; x < xbin_count.size(); x++){
		  if (mirrorbin == xbin_count[x]){exist = x; break;}
		}
		if (exist == -1){
		  xbin_count.push_back(mirrorbin);
		  xbin_dev.push_back(integrand_D[2][i_a]);
		  if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_CV[i_s].push_back(integrand_D_CV[i_s][2][i_a]);}}
		}
		else{
		  xbin_dev[exist] += integrand_D[2][i_a];
		  if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_CV[i_s][exist] += integrand_D_CV[i_s][2][i_a];}}
		}
	      }
	    }
	  }
	}
	for (int x = 0; x < xbin_count.size(); x++){
	  //	  bin_dev[i_d][xbin_count[x]] += pow(xbin_dev[x], 2);
	  bin_weight2[i_d][xbin_count[x]] += pow(xbin_dev[x], 2);
	  if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight2_CV[i_s][i_d][xbin_count[x]] += pow(xbin_dev_CV[i_s][x], 2);}}

	    // why different from double-differential distrbution ???
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


    // maybe introduce a   vector<vector<vector<vector<int> > > > bin_phasespace(dat.size(), vector<vector<vector<int> > > > (3)) , which contains only the i_a values that contribute with 0, 1 or 2 value at a particular point ???



// mirror-bin modifications not jet implemented for dddistributions !!!

void observable_set::determine_double_distribution(){
  static Logger logger("observable_set::determine_double_distribution");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_ddd = 0; i_ddd < dddat.size(); i_ddd++){
    int i_d = dat.size() + i_ddd;
    logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << "   i_ddd = " << i_ddd << "   n_ps = " << n_ps << endl;
    
    // temporary solution:
    // dddat[i_ddd].distribution_1).bin -> dat[dddat[i_ddd].d_1].bin
    // dddat[i_ddd].distribution_2).bin -> dat[dddat[i_ddd].d_2].bin
    // distribution_1/2 should be references to the respective distributions instead of copies !!!
    
    if (n_ps == 1){  // only one phase-space
      logger << LOG_DEBUG_VERBOSE << "CASE   n_ps == 1" << endl;
      if ((dat[dddat[i_ddd].d_1].bin[0] == -1) || (dat[dddat[i_ddd].d_2].bin[0] == -1)){dddat[i_ddd].bin[0] = -1;}
      else {dddat[i_ddd].bin[0] = dat[dddat[i_ddd].d_1].bin[0] * (dddat[i_ddd].distribution_2).n_bins + dat[dddat[i_ddd].d_2].bin[0];}
      //  to be generalized for multi-differential distributions

      int presentbin = dddat[i_ddd].bin[0];
      
      //  check if dddat[i_ddd].step factor is correct here !!! (shifted to summary for single-differential distributions)
      if (presentbin != -1){
	int mirrorbin = -1;
	logger << LOG_DEBUG_VERBOSE << "i_d = " << i_d << "   i_ddd = " << i_ddd << "   presentbin = " << presentbin << "   mirrorbin = " << mirrorbin << endl;
	logger << LOG_DEBUG_VERBOSE << "bin_weight.size() = " << bin_weight.size() << endl;
	logger << LOG_DEBUG_VERBOSE << "bin_weight[" << i_d << "].size() = " << bin_weight[i_d].size() << endl;
	bin_weight[i_d][presentbin] += integrand_D[1][0] / dddat[i_ddd].step;
	if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][presentbin] += integrand_D_CV[i_s][1][0] / dddat[i_ddd].step;}}
	bin_counts[i_d][presentbin]++;
	if (csi->type_parton[0][1] != csi->type_parton[0][2]){
	  //  to be generalized for multi-differential distributions !!!
	  if (((dddat[i_ddd].distribution_1).symm == 0) && ((dddat[i_ddd].distribution_2).symm == 0)){ // asymmetric distribution 1, asymmetric distribution 2
	    mirrorbin = presentbin;
	  }
	  else if (((dddat[i_ddd].distribution_1).symm == 0) && ((dddat[i_ddd].distribution_2).symm == 1)){ // asymmetric distribution 1, symmetric distribution 2
	    mirrorbin = dat[dddat[i_ddd].d_1].bin[0] * (dddat[i_ddd].distribution_2).n_bins + ((dddat[i_ddd].distribution_2).n_bins - 1 - dat[dddat[i_ddd].d_2].bin[0]);
	  }
	  else if (((dddat[i_ddd].distribution_1).symm == 1) && ((dddat[i_ddd].distribution_2).symm == 0)){ // symmetric distribution 1, asymmetric distribution 2
	    mirrorbin = ((dddat[i_ddd].distribution_1).n_bins - 1 - dat[dddat[i_ddd].d_1].bin[0]) * (dddat[i_ddd].distribution_2).n_bins + dat[dddat[i_ddd].d_2].bin[0];
	  }
	  else if (((dddat[i_ddd].distribution_1).symm == 1) && ((dddat[i_ddd].distribution_2).symm == 1)){ // symmetric distribution 1, symmetric distribution 2
	    mirrorbin = dat[i_ddd].n_bins - 1 - presentbin;
	  }
	  else {logger << LOG_FATAL << "No allowed case chosen." << endl; exit(1);}
	  
	  bin_weight[i_d][mirrorbin] += (integrand_D[2][0]) / dddat[i_ddd].step;
	  if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][presentbin] += integrand_D_CV[i_s][2][0] / dddat[i_ddd].step;}}
	  bin_counts[i_d][mirrorbin]++;
	}
	bin_weight2[i_d][presentbin] += pow(integrand_D[1][0] / dddat[i_ddd].step, 2);
	if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight2_CV[i_s][i_d][presentbin] += pow(integrand_D_CV[i_s][1][0] / dddat[i_ddd].step, 2);}}
	
	if (integrand_D[2][0] != 0.){
	  bin_weight2[i_d][presentbin] += pow(integrand_D[2][0] / dddat[i_ddd].step, 2);
	  if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight2_CV[i_s][i_d][presentbin] += pow(integrand_D_CV[i_s][2][0] / dddat[i_ddd].step, 2);}}
	}
      }
    }
    else if (n_ps > 1){  // more than one phase-space (n_ps > 1)
      logger << LOG_DEBUG_VERBOSE << "CASE   n_ps > 1" << endl;
      vector<int> xbin_count(0), xbin_count_2(0);
      vector<double> xbin_dev(0), xbin_dev_2(0);
      vector<vector<double> > xbin_dev_CV(n_scales_CV), xbin_dev_2_CV(n_scales_CV);
      
      for (int i_a = 0; i_a < n_ps; i_a++){
	//x//	if ((dat[dddat[i_ddd].d_1].bin[i_a] == -1) || (dat[dddat[i_ddd].d_2].bin[i_a] == -1)){dddat[i_ddd].bin[i_a] = -1;}
	//x//	else {dddat[i_ddd].bin[i_a] = dat[dddat[i_ddd].d_1].bin[i_a] * (dddat[i_ddd].distribution_2).n_bins + dat[dddat[i_ddd].d_2].bin[i_a];}
	logger << LOG_DEBUG_VERBOSE << "dddat.size() = " << dddat.size() << endl;
	logger << LOG_DEBUG_VERBOSE << "dddat[" << i_ddd << "].bin.size() = " << dddat[i_ddd].bin.size() << endl;
	logger << LOG_DEBUG_VERBOSE << "dddat[" << i_ddd << "].mirror_bin.size() = " << dddat[i_ddd].mirror_bin.size() << endl;
	logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << endl;
	int presentbin = dddat[i_ddd].bin[i_a];
	logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << endl;
	int mirrorbin = dddat[i_ddd].mirror_bin[i_a];
	logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << endl;
	logger << LOG_DEBUG_VERBOSE << "dddat[" << setw(2) << i_ddd << "].bin[" << setw(2) << i_a << "] = " << setw(3) << dddat[i_ddd].bin[i_a] << "   dddat[" << setw(2) << i_ddd << "].mirror_bin[" << setw(2) << i_a << "] = " << setw(3) << dddat[i_ddd].mirror_bin[i_a] << endl;
	//x//	int mirrorbin = -1;
	if (presentbin != -1){
	  bin_weight[i_d][presentbin] += (integrand_D[1][i_a]) / dddat[i_ddd].step;
	  if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][presentbin] += (integrand_D_CV[i_s][1][i_a]) / dddat[i_ddd].step;}}
	  bin_counts[i_d][presentbin]++;

	  if (csi->type_parton[0][1] != csi->type_parton[0][2]){
	    /*//x//
	    if (((dddat[i_ddd].distribution_1).symm == 0) && ((dddat[i_ddd].distribution_2).symm == 0)){ // asymmetric distribution 1, asymmetric distribution 2
	      mirrorbin = presentbin;
	    }
	    else if (((dddat[i_ddd].distribution_1).symm == 0) && ((dddat[i_ddd].distribution_2).symm == 1)){ // asymmetric distribution 1, symmetric distribution 2
	      mirrorbin = dat[dddat[i_ddd].d_1].bin[i_a] * (dddat[i_ddd].distribution_2).n_bins + ((dddat[i_ddd].distribution_2).n_bins - 1 - dat[dddat[i_ddd].d_2].bin[i_a]);
	    }
	    else if (((dddat[i_ddd].distribution_1).symm == 1) && ((dddat[i_ddd].distribution_2).symm == 0)){ // symmetric distribution 1, asymmetric distribution 2
	      mirrorbin = ((dddat[i_ddd].distribution_1).n_bins - 1 - dat[dddat[i_ddd].d_1].bin[i_a]) * (dddat[i_ddd].distribution_2).n_bins + dat[dddat[i_ddd].d_2].bin[i_a];
	    }
	    else if (((dddat[i_ddd].distribution_1).symm == 1) && ((dddat[i_ddd].distribution_2).symm == 1)){ // symmetric distribution 1, symmetric distribution 2
	      mirrorbin = dddat[i_ddd].n_bins - 1 - presentbin;
	    }
	    else {logger << LOG_FATAL << "No allowed case chosen." << endl; exit(1);}
	    *///x//
	    // may not be called if mirror_bin not set !!!
	    // introduce some variable in xdistribution !!!
	    bin_weight[i_d][mirrorbin] += (integrand_D[2][i_a]) / dddat[i_ddd].step;
	    if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight_CV[i_s][i_d][mirrorbin] += (integrand_D_CV[i_s][2][i_a]) / dddat[i_ddd].step;}}
	    bin_counts[i_d][mirrorbin]++;
	  }

	  int exist = -1;
	  for (int x = 0; x < xbin_count.size(); x++){
	    if (presentbin == xbin_count[x]){exist = x; break;}
	  }
	  if (exist == -1){
	    xbin_count.push_back(presentbin);
	    xbin_dev.push_back(integrand_D[1][i_a] / dddat[i_ddd].step);
	    if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_CV[i_s].push_back(integrand_D_CV[i_s][1][i_a] / dddat[i_ddd].step);}}
	  }
	  else {
	    xbin_dev[exist] += integrand_D[1][i_a] / dddat[i_ddd].step;
	    if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_CV[i_s][exist] += integrand_D_CV[i_s][1][i_a] / dddat[i_ddd].step;}}
	  }
	  if (integrand_D[2][i_a] != 0.){  // should be the same as integrand_D_CV[i_s][2][0] != 0.
	    exist = -1;
	    // Should actually work for all non-cumulative distributions
	    //	    if (dat[i_ddd].symm == 0){ // asymmetric distribution
	    for (int x = 0; x < xbin_count_2.size(); x++){
	      if (mirrorbin == xbin_count_2[x]){exist = x; break;}
	    }
	    if (exist == -1){
	      xbin_count_2.push_back(mirrorbin);
	      xbin_dev_2.push_back(integrand_D[2][i_a] / dddat[i_ddd].step);
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_2_CV[i_s].push_back(integrand_D_CV[i_s][2][i_a] / dddat[i_ddd].step);}}
	    }
	    else{
	      xbin_dev_2[exist] += integrand_D[2][i_a] / dddat[i_ddd].step;
	      if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){xbin_dev_2_CV[i_s][exist] += integrand_D_CV[i_s][2][i_a] / dddat[i_ddd].step;}}
	    }
	  }
	}
      }
      
      // Check what gives the better error estimate !!!
      for (int x = 0; x < xbin_count.size(); x++){
	bin_weight2[i_d][xbin_count[x]] += pow(xbin_dev[x], 2);
	if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight2_CV[i_s][i_d][xbin_count[x]] += pow(xbin_dev_CV[i_s][x], 2);}}
      }
      for (int x = 0; x < xbin_count_2.size(); x++){
	bin_weight2[i_d][xbin_count_2[x]] += pow(xbin_dev_2[x], 2);
	if (switch_CV){for (int i_s = 0; i_s < n_scales_CV; i_s++){bin_weight2_CV[i_s][i_d][xbin_count_2[x]] += pow(xbin_dev_2_CV[i_s][x], 2);}}
      }
    }
    else {
      logger << LOG_FATAL << "No allowed n_ps value chosen." << endl; exit(1);
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::preparation_distribution_TSV(){
  static Logger logger("observable_set::preparation_distribution_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  distribution_no_qTcut   (contains all relevant qTcut numbers)
  //  1st index -> counter of relevant no_qTcut 
  //  entry     -> no_qTcut (finally in ascending order)

  //  distribution_no_qTcut_phasespace:
  //  1st index -> number of qTcut from "distribution_no_qTcut" (no_qTcut available via distribution_no_qTcut[<1st index>])
  //  2nd index -> counter of relevant phasespaces (i_a) with the respective cut value
  //  entry     -> no_phasespace (i_a (in ascending order?))

  distribution_no_qTcut = cut_ps;

  if (n_ps == 1){
    distribution_no_qTcut_phasespace.resize(1, vector<int> (1, 0)); // -> (1, 0) only 1 phasespace value, namely 0; resize(1, ...) -> only 1 qTcut value
  }
  else if (n_ps > 1){
    logger << LOG_DEBUG_VERBOSE << "CASE   n_ps > 1" << endl;
    for (int i_q = 0; i_q < distribution_no_qTcut.size(); i_q++){logger << LOG_DEBUG_VERBOSE << "distribution_no_qTcut[" << i_q << "] = " << distribution_no_qTcut[i_q] << endl;}
    // i_q > 0, because the 0 entry always remains unchanged!
    logger << LOG_DEBUG_VERBOSE << "sort   distribution_no_qTcut" << endl;
    sort(distribution_no_qTcut.begin(), distribution_no_qTcut.end());
    for (int i_q = 0; i_q < distribution_no_qTcut.size(); i_q++){logger << LOG_DEBUG_VERBOSE << "distribution_no_qTcut[" << i_q << "] = " << distribution_no_qTcut[i_q] << endl;}
    for (int i_q = n_ps - 1; i_q > 0; i_q--){  // i_q > 0, because the 0 entry always remains unchanged!
      if (distribution_no_qTcut[i_q] == distribution_no_qTcut[i_q - 1]){distribution_no_qTcut.erase(distribution_no_qTcut.begin() + i_q);}
    }
    logger << LOG_DEBUG_VERBOSE << "remove identical entries from   distribution_no_qTcut" << endl;
    logger << LOG_DEBUG_VERBOSE << "distribution_no_qTcut.size() = " << distribution_no_qTcut.size() << endl;
    for (int i_q = 0; i_q < distribution_no_qTcut.size(); i_q++){logger << LOG_DEBUG_VERBOSE << "distribution_no_qTcut[" << i_q << "] = " << distribution_no_qTcut[i_q] << endl;}
    logger << LOG_DEBUG_VERBOSE << "remove -1 entries from   distribution_no_qTcut" << endl;
    if (distribution_no_qTcut[0] == -1){distribution_no_qTcut.erase(distribution_no_qTcut.begin());}
    
    // !!!
    logger << LOG_DEBUG_VERBOSE << "distribution_no_qTcut.size() = " << distribution_no_qTcut.size() << endl;
    for (int j_q = 0; j_q < distribution_no_qTcut.size(); j_q++){
      logger << LOG_DEBUG_VERBOSE << "distribution_no_qTcut[" << j_q << "] = " << distribution_no_qTcut[j_q] << endl;
    }
    
    distribution_no_qTcut_phasespace.clear();
    distribution_no_qTcut_phasespace.resize(distribution_no_qTcut.size());
    for (int i_a = 0; i_a < n_ps; i_a++){
      for (int i_q = 0; i_q < distribution_no_qTcut.size(); i_q++){
	if (distribution_no_qTcut[i_q] == cut_ps[i_a]){distribution_no_qTcut_phasespace[i_q].push_back(i_a);}
      }
    }
    
    for (int i_q = 0; i_q < distribution_no_qTcut.size(); i_q++){
      stringstream out_dnq;
      for (int i_qa = 0; i_qa < distribution_no_qTcut_phasespace[i_q].size(); i_qa++){
	out_dnq << distribution_no_qTcut_phasespace[i_q][i_qa] << " ";
      }
      logger << LOG_DEBUG << "distribution_no_qTcut_phasespace [" << i_q << "] @ " << setw(3) << distribution_no_qTcut[i_q] << " = " << out_dnq.str() << endl;
    }

    // distribution_no_qTcut
    // 1st entry: all relevant values of cut_ps, in ascending (?) order, appearing only once each (-1 entries removed)

    /*
      e.g.:
      distribution_no_qTcut[0] = 7
      distribution_no_qTcut[1] = 15
      -->
      i_q = 0 - 7: bin_weight_TSV[i_q] = ...[distribution_no_qTcut[0]]
      i_q = 8 - 15: bin_weight_TSV[i_q] = ...[distribution_no_qTcut[1]]
      i_q = 16 - output_n_qTcut: bin_weight_TSV[i_q] = 0
    */
    
    // distribution_no_qTcut_phasespace:
    // 1st entry: cut value from "distribution_no_qTcut"
    // 2nd entry: phasespace (i_a) with that cut value
  }


  // maybe clearing of other variables for TSV distributions here ???

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_single_distribution_TSV(){
  static Logger logger("observable_set::determine_single_distribution_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_d = 0; i_d < dat.size(); i_d++){dat[i_d].initialization_distribution_bin();}

  for (int i_d = 0; i_d < dat.size(); i_d++){
    if (n_ps == 1){
      logger << LOG_DEBUG_VERBOSE << "CASE   n_ps == 1" << endl;
      if (dat[i_d].bin[0] != -1){
	if (csi->class_contribution_IRcut_implicit){
	  //	if (type_contribution == "CT" || type_contribution == "CT2" || type_contribution == "L2CT"){
	  // special case wrt. qTcut, as all qTcut values are simultaneously filled with explicitly qTcut-dependent results
	  for (int i_q = 0; i_q <= cut_ps[0]; i_q++){
	    int min = 0, max = 0;
	    if (dat[i_d].typeCumulative != CUMULATIVE_NONE){
	      min = dat[i_d].distribution_no_bin[0]; 
	      max = dat[i_d].distribution_no_bin_max[0];
	    }
	    else {
	      min = dat[i_d].distribution_no_bin[0];
	      max = dat[i_d].distribution_no_bin[0] + 1;
	    }
	    
	    for (int i_b = min; i_b < max; i_b++){ // was '<=' before !!!

	      bin_count_TSV[i_q][i_d][i_b]++;

	      for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
		if (switch_distribution_TSV[i_s] > 0){
		  for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		    for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		      // maybe improve for cumulative distributions ???
		      // i_b -> only one bin: dat[i_d].bin[0]
		      // i_q -> n_qTcut different values: cut_ps[0] is always -1  or  n_qTcut (-1?)
		      if (dat[i_d].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2]){
			bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][0][0];
			bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += pow(ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][0][0], 2);
		      }
		      else {
			bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][1][0];
			bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += pow(ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][1][0], 2);
			bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].n_bins - 1 - i_b] += ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][2][0];
			bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].n_bins - 1 - i_b] += pow(ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][2][0], 2);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	else {
	  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	    if (switch_distribution_TSV[i_s] > 0){
	      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		  if (dat[i_d].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2]){
		    change_weight_TSV[i_s][i_r][i_f] = ps_integrand_TSV[i_s][i_r][i_f][0][0];
		    change_weight2_TSV[i_s][i_r][i_f] = pow(change_weight_TSV[i_s][i_r][i_f], 2);
		  }
		  else {
		    change_weight_TSV[i_s][i_r][i_f] = ps_integrand_TSV[i_s][i_r][i_f][1][0];
		    change_weight2_TSV[i_s][i_r][i_f] = pow(change_weight_TSV[i_s][i_r][i_f], 2);
		    change_weight_2nd_TSV[i_s][i_r][i_f] = ps_integrand_TSV[i_s][i_r][i_f][2][0];
		    change_weight2_2nd_TSV[i_s][i_r][i_f] = pow(change_weight_2nd_TSV[i_s][i_r][i_f], 2);
		  }
		}
	      }
	    }
	  }

	  for (int i_q = 0; i_q <= cut_ps[0]; i_q++){  //  !!! check '<' or '<=' !!! 2nd ?
	    int min = 0, max = 0;
	    if (dat[i_d].typeCumulative != CUMULATIVE_NONE){
	      min = dat[i_d].distribution_no_bin[0];
	      max = dat[i_d].distribution_no_bin_max[0];
	    }
	    else {
	      min = dat[i_d].distribution_no_bin[0];
	      max = dat[i_d].distribution_no_bin[0] + 1;
	    }
	    
	    for (int i_b = min; i_b < max; i_b++){ // was '<=' before !!!
	      bin_count_TSV[i_q][i_d][i_b]++;
      	      for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
		if (switch_distribution_TSV[i_s] > 0){
		  for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		    for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += change_weight_TSV[i_s][i_r][i_f];
		      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += change_weight2_TSV[i_s][i_r][i_f];
		      if (!(dat[i_d].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
			bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].n_bins - 1 - i_b] += change_weight_2nd_TSV[i_s][i_r][i_f];
			bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].n_bins - 1 - i_b] += change_weight2_2nd_TSV[i_s][i_r][i_f];
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else if (n_ps > 1){
      logger << LOG_DEBUG_VERBOSE << "CASE   n_ps > 1" << endl;

      if (dat[i_d].typeCumulative != CUMULATIVE_NONE){
	logger << LOG_DEBUG_VERBOSE << "CASE   n_ps > 1 && dat[" <<  i_d << "].typeCumulative != CUMULATIVE_NONE" << endl;
	logger << LOG_DEBUG_VERBOSE << "dat[" << i_d << "].xdistribution_name = " << dat[i_d].xdistribution_name << endl;
	logger << LOG_DEBUG_VERBOSE << "dat[i_d].distribution_no_bin_all.size() = " << dat[i_d].distribution_no_bin_all.size() << endl;
	
	if (dat[i_d].distribution_no_bin_all.size() != 0){
	  for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	    if (switch_distribution_TSV[i_s] > 0){
	      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		  change_qTcut_bin_weight_TSV[i_s][i_r][i_f].clear();
		  change_qTcut_bin_weight2_TSV[i_s][i_r][i_f].clear();
		  change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f].clear();
		  change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f].clear();
		  
		  change_qTcut_bin_weight_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dat[i_d].distribution_no_bin_all.size() - 1, 0.));
		  change_qTcut_bin_weight2_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dat[i_d].distribution_no_bin_all.size() - 1, 0.));
		  change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dat[i_d].distribution_no_bin_all.size() - 1, 0.));
		  change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dat[i_d].distribution_no_bin_all.size() - 1, 0.));
		}
	      }
	    }
	  }

	  for (int i_b = 0; i_b < dat[i_d].distribution_no_bin_all.size() - 1; i_b++){
	    for (int i_q = distribution_no_qTcut.size() - 1; i_q >= 0; i_q--){
	      for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
		if (switch_distribution_TSV[i_s] > 0){
		  logger << LOG_DEBUG_VERBOSE << "i_b = " << setw(4) << i_b << "   " << "i_q = " << setw(4) << i_q << "   " << "i_s = " << setw(4) << i_s << "   " << "before" << endl;
		  for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		    for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		      if (dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][i_q].size() > 0){
			for (int i_a = 0; i_a < dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][i_q].size(); i_a++){
			  if (dat[i_d].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2]){
			    change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][0][dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			  }
			  else {
			    change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][1][dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			    change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][2][dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			  }
			}
			change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][i_q][i_b] = pow(change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b], 2);
			if (!(dat[i_d].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
			  change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][i_q][i_b] = pow(change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q][i_b], 2);
			}
		      }
		      else if (change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] == 0.){
			// can happen in case of non-overlapping intervals, e.g. [0:3] v [7:21]
		      }
		      else {
			logger << LOG_WARN << "n_ps != 1   ---   dat[" << i_d << "].typeCumulative != CUMULATIVE_NONE   ---   Should not happen now!" << endl;
			change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q + 1][i_b];
			change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][i_q + 1][i_b];
			if (!(dat[i_d].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
			  change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q + 1][i_b];
			  change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][i_q + 1][i_b];
			}
		      }
		    }
		  }
		  logger << LOG_DEBUG_VERBOSE << "i_b = " << setw(4) << i_b << "   " << "i_q = " << setw(4) << i_q << "   " << "i_s = " << setw(4) << i_s << "   " << "after" << endl;
		}
	      }
	    }
	  }
	
		    
	  for (int i_b = 0; i_b < dat[i_d].distribution_no_bin_all.size() - 1; i_b++){
	    int i_b_min = dat[i_d].distribution_no_bin_all[i_b];
	    int i_b_max = dat[i_d].distribution_no_bin_all[i_b + 1];
	    logger << LOG_DEBUG_VERBOSE << "i_b = " << setw(4) << i_b << "   " << "i_b_min = " << setw(4) << i_b_min << "   " << "i_b_max = " << setw(4) << i_b_max << "   " << "before" << endl;
	    for (int j_b = i_b_min; j_b < i_b_max; j_b++){
	      int j_q = 0;
	      int start_i_q = distribution_no_qTcut[0];
	      for (int i_q = start_i_q; i_q < output_n_qTcut; i_q++){
		if (i_q > start_i_q && i_q == distribution_no_qTcut[j_q]){j_q++;}
		logger << LOG_DEBUG_VERBOSE << "i_q = " << setw(4) << i_q << "   " << "start_i_q = " << setw(4) << start_i_q << "   " << "j_q = " << setw(4) << j_q << "   " << "distribution_no_qTcut[" << j_q << "] = " << distribution_no_qTcut[j_q] << endl;
		logger << LOG_DEBUG_VERBOSE << "dat[" << i_d << "].distribution_no_bin_no_qTcut_phasespace.size() = " << dat[i_d].distribution_no_bin_no_qTcut_phasespace.size() << endl;
		logger << LOG_DEBUG_VERBOSE << "dat[" << i_d << "].distribution_no_bin_no_qTcut_phasespace[" << i_b << "].size() = " << dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b].size() << endl;
		logger << LOG_DEBUG_VERBOSE << "dat[" << i_d << "].distribution_no_bin_no_qTcut_phasespace[" << i_b << "][" << j_q << "].size() = " << dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][j_q].size() << endl;
		//		logger << LOG_DEBUG_VERBOSE << " = " <<  << endl;
		logger << LOG_DEBUG_VERBOSE << "bin_count_TSV.size() = " << bin_count_TSV.size() << endl;
		logger << LOG_DEBUG_VERBOSE << "bin_count_TSV[" << i_q << "].size() = " << bin_count_TSV[i_q].size() << endl;
		logger << LOG_DEBUG_VERBOSE << "bin_count_TSV[" << i_q << "][" << i_d << "].size() = " << bin_count_TSV[i_q][i_d].size() << endl;

		bin_count_TSV[i_q][i_d][j_b] += dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][j_q].size();
		
		logger << LOG_DEBUG_VERBOSE << "distribution_no_qTcut.size() = " << distribution_no_qTcut.size() << endl;
		logger << LOG_DEBUG_VERBOSE << "distribution_no_qTcut[" << j_q << "] = " << distribution_no_qTcut[j_q] << "   i_b = " << i_b << endl;
		for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
		  if (switch_distribution_TSV[i_s] > 0){
		    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
			logger << LOG_DEBUG_VERBOSE << "change_qTcut_bin_weight_TSV.size() = " << change_qTcut_bin_weight_TSV.size() << endl;
			
			logger << LOG_DEBUG_VERBOSE << "change_qTcut_bin_weight_TSV[" << i_s << "][" << i_r << "][" << i_f << "].size() = " << change_qTcut_bin_weight_TSV[i_s][i_r][i_f].size() << endl;
			logger << LOG_DEBUG_VERBOSE << "change_qTcut_bin_weight_TSV[" << i_s << "][" << i_r << "][" << i_f << "][" << j_q << "].size() = " << change_qTcut_bin_weight_TSV[i_s][i_r][i_f][j_q].size() << endl;
			//			logger << LOG_DEBUG_VERBOSE << "change_qTcut_bin_weight_TSV[" << i_s << "][" << i_r << "][" << i_f << "][" << distribution_no_qTcut[j_q]<< "].size() = " << change_qTcut_bin_weight_TSV[i_s][i_r][i_f][distribution_no_qTcut[j_q]].size() << endl;

			// 1st: asymmetric (no symmetry wrt. rotating the positive beam axis to the negative one and vice versa) distribution
			// 2nd: identical initial-state partons -> no contribution from rotated point (double-counting)
			bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][j_b] += change_qTcut_bin_weight_TSV[i_s][i_r][i_f][j_q][i_b];
			//			bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][j_b] += change_qTcut_bin_weight_TSV[i_s][i_r][i_f][distribution_no_qTcut[j_q]][i_b];
			bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][j_b] += change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][j_q][i_b];
			//			bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][j_b] += change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][distribution_no_qTcut[j_q]][i_b];
			if (!(dat[i_d].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
			  bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].n_bins - 1 - j_b] += change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][j_q][i_b];
			  //			  bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].n_bins - 1 - j_b] += change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][distribution_no_qTcut[j_q]][i_b];
			  bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].n_bins - 1 - j_b] += change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][j_q][i_b];
			  //			  bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].n_bins - 1 - j_b] += change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][distribution_no_qTcut[j_q]][i_b];
			}
		      }
		    }
		  }
		}
		// different from CUMULATIVE_NONE case !!!
		if (j_q == distribution_no_qTcut.size() - 1){break;}// ???
	      }
	    }
	    logger << LOG_DEBUG_VERBOSE << "i_b = " << setw(4) << i_b << "   " << "i_b_min = " << setw(4) << i_b_min << "   " << "i_b_max = " << setw(4) << i_b_max << "   " << "after" << endl;
	  }
	}
      }
      else if (dat[i_d].typeCumulative == CUMULATIVE_NONE){
	logger << LOG_DEBUG_VERBOSE << "CASE   n_ps > 1 && dat[" <<  i_d << "].typeCumulative == CUMULATIVE_NONE" << endl;

	for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	  if (switch_distribution_TSV[i_s] > 0){
	    for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
	      for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		change_qTcut_bin_weight_TSV[i_s][i_r][i_f].clear();
		change_qTcut_bin_weight2_TSV[i_s][i_r][i_f].clear();
		change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f].clear();
		change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f].clear();
		
		change_qTcut_bin_weight_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dat[i_d].distribution_no_bin.size(), 0.));
		change_qTcut_bin_weight2_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dat[i_d].distribution_no_bin.size(), 0.));
		change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dat[i_d].distribution_no_bin.size(), 0.));
		change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dat[i_d].distribution_no_bin.size(), 0.));
	      }
	    }
	  }
	}

	for (int i_b = 0; i_b < dat[i_d].distribution_no_bin.size(); i_b++){
	  for (int i_q = distribution_no_qTcut.size() - 1; i_q >= 0; i_q--){
	    for (int i_ba = 0; i_ba < dat[i_d].distribution_no_bin_phasespace[i_b].size(); i_ba++){
	      for (int i_qa = 0; i_qa < distribution_no_qTcut_phasespace[i_q].size(); i_qa++){
		if (dat[i_d].distribution_no_bin_phasespace[i_b][i_ba] == distribution_no_qTcut_phasespace[i_q][i_qa]){
		  for (int j_q = 0; j_q <= distribution_no_qTcut[i_q]; j_q++){
		    bin_count_TSV[j_q][i_d][dat[i_d].bin[dat[i_d].distribution_no_bin_phasespace[i_b][i_ba]]]++;
		  }
		  // check if correct !!!
		}
	      }
	    }




	    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	      if (switch_distribution_TSV[i_s] > 0){
		for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		    // !!! change to fix bug !!!
		    if (i_q < distribution_no_qTcut.size() - 1){
		      // start with weight from higher qTcut values:
		      change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q + 1][i_b];
		      change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q + 1][i_b];
		    }
		    // !!! change to fix bug end !!!
		    if (dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][i_q].size() > 0){
		      for (int i_a = 0; i_a < dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][i_q].size(); i_a++){
			if (dat[i_d].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2]){
			  change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][0][dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			}
			else {
			  change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][1][dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			  change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][2][dat[i_d].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			}
		      }
		      change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][i_q][i_b] = pow(change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b], 2);
		      if (!(dat[i_d].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
			change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][i_q][i_b] = pow(change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q][i_b], 2);
		      }
		    }
		    else if (change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] == 0.){
		      // can happen in case of non-overlapping intervals, e.g. [0:3] v [7:21], at least in CUMULATIVE case
		      logger << LOG_DEBUG << "Why would this happen here (i_d = " << i_d << ") ??? Maybe because of dipoles with precisely opposite signs?" << endl;
		      logger << LOG_DEBUG << "change_qTcut_bin_weight_TSV[i_s][" << i_r << "][" << i_f << "][" << i_q << "][" << i_b << "] = " << change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] << endl;
		    }
		    else {
		      double temp = change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b];
		      change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q + 1][i_b];
		      change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][i_q + 1][i_b];
		      if (!(dat[i_d].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
			change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q + 1][i_b];
			change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][i_q + 1][i_b];
		      }
		      if (temp != change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b]){
			logger << LOG_DEBUG << "Why would this happen (i_d = " << i_d << ") ???" << endl;
			logger << LOG_DEBUG << "change_qTcut_bin_weight_TSV[i_s][" << i_r << "][" << i_f << "][" << i_q << "][" << i_b << "] = " << temp << " --> " << change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] << endl;
		      }
		      logger << LOG_DEBUG << "change_qTcut_bin_weight_TSV[" << i_s << "][" << i_r << "][" << i_f << "][" << i_q << "][" << i_b << "] = " << change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] << endl;
		    }
		  }
		}
	      }
	    }
	  }
	}
	
	for (int i_b = 0; i_b < dat[i_d].distribution_no_bin.size(); i_b++){
	  int j_q = 0;  // j_q  counts the different values of  change_qTcut_bin_weight_TSV[i_s]  etc. that belong to the respective  i_q  values.
	  for (int i_q = 0; i_q < output_n_qTcut; i_q++){

	    ///	    bin_count_TSV[j_q][i_d][i_b]++;

	    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	      if (switch_distribution_TSV[i_s] > 0){
		for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		    // 1st: asymmetric (no symmetry wrt. rotating the positive beam axis to the negative one and vice versa) distribution
		    // 2nd: identical initial-state partons -> no contribution from rotated point (double-counting)
		    bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].distribution_no_bin[i_b]] += change_qTcut_bin_weight_TSV[i_s][i_r][i_f][j_q][i_b];
		    bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].distribution_no_bin[i_b]] += change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][j_q][i_b];
		    if (!(dat[i_d].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
		      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].n_bins - 1 - dat[i_d].distribution_no_bin[i_b]] += change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][j_q][i_b];
		      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][dat[i_d].n_bins - 1 - dat[i_d].distribution_no_bin[i_b]] += change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][j_q][i_b];
		    }
		  }
		}
	      }
	    }
	    if (i_q == distribution_no_qTcut[j_q]){j_q++;}
	    if (j_q == distribution_no_qTcut.size()){break;}
	  }
	}
      }
      else {
	logger << LOG_FATAL << "Wrong typeCumulative" << endl; exit(1);
      }
    }
    else {
      logger << LOG_FATAL << "Wrong n_ps value (< 1)" << endl; exit(1);
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_double_distribution_TSV(){
  static Logger logger("observable_set::determine_double_distribution_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_ddd = 0; i_ddd < dddat.size(); i_ddd++){dddat[i_ddd].initialization_distribution_bin();}
  
  logger << LOG_DEBUG_VERBOSE << "n_ps = " << n_ps << "   type_contribution = " << type_contribution << endl;
  
  if (n_ps == 1){
    /*///
    for (int i_ddd = 0; i_ddd < dddat.size(); i_ddd++){
      int i_d = dat.size() + i_ddd;
      if ((dddat[i_ddd].distribution_1).typeCumulative != CUMULATIVE_NONE && (dddat[i_ddd].distribution_2).typeCumulative != CUMULATIVE_NONE){
	// to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
      }
      else if ((dddat[i_ddd].distribution_1).typeCumulative == CUMULATIVE_NONE && (dddat[i_ddd].distribution_2).typeCumulative != CUMULATIVE_NONE){
	// to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
      }
      else if ((dddat[i_ddd].distribution_1).typeCumulative != CUMULATIVE_NONE && (dddat[i_ddd].distribution_2).typeCumulative == CUMULATIVE_NONE){
	// to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
      }
      else if ((dddat[i_ddd].distribution_1).typeCumulative == CUMULATIVE_NONE && (dddat[i_ddd].distribution_2).typeCumulative == CUMULATIVE_NONE){
	// 'bin_count_TSV' will need to remain part of this function - or could be made part of corresponding one in  dddistribution
	for (int i_q = 0; i_q <= cut_ps[0]; i_q++){bin_count_TSV[i_q][i_d][dddat[i_ddd].bin[0]]++;}
	// ??? < or <= ???  // mirrored contributions ???
      }
    }
    *////
  }
  
  else if (n_ps > 1){
    for (int i_ddd = 0; i_ddd < dddat.size(); i_ddd++){
      int i_d = dat.size() + i_ddd;
      
      if ((dddat[i_ddd].distribution_1).typeCumulative != CUMULATIVE_NONE && (dddat[i_ddd].distribution_2).typeCumulative != CUMULATIVE_NONE){
	// to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
      }
      else if ((dddat[i_ddd].distribution_1).typeCumulative == CUMULATIVE_NONE && (dddat[i_ddd].distribution_2).typeCumulative != CUMULATIVE_NONE){
	// to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
      }
      else if ((dddat[i_ddd].distribution_1).typeCumulative != CUMULATIVE_NONE && (dddat[i_ddd].distribution_2).typeCumulative == CUMULATIVE_NONE){
	// to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
      }
      else if ((dddat[i_ddd].distribution_1).typeCumulative == CUMULATIVE_NONE && (dddat[i_ddd].distribution_2).typeCumulative == CUMULATIVE_NONE){
	for (int i_b = 0; i_b < dddat[i_ddd].distribution_no_bin.size(); i_b++){
	  for (int i_q = 0; i_q < distribution_no_qTcut.size(); i_q++){
	    for (int i_ba = 0; i_ba < dddat[i_ddd].distribution_no_bin_phasespace[i_b].size(); i_ba++){
	      for (int i_qa = 0; i_qa < distribution_no_qTcut_phasespace[i_q].size(); i_qa++){
		if (dddat[i_ddd].distribution_no_bin_phasespace[i_b][i_ba] == distribution_no_qTcut_phasespace[i_q][i_qa]){
		  // 'bin_count_TSV' will need to remain part of this function - or could be made part of corresponding one in  dddistribution
		  //		  bin_count_TSV[i_q][i_d][dddat[i_ddd].bin[dddat[i_ddd].distribution_no_bin_phasespace[i_b][i_ba]]]++;
		  // ??? < or <= ???  // mirrored contributions ???
		}
	      }
	    }
	  }
	}
      }
    }
  }
  

  for (int i_ddd = 0; i_ddd < dddat.size(); i_ddd++){
    int i_d = dat.size() + i_ddd;
    if (n_ps == 1){
      logger << LOG_DEBUG_VERBOSE << "CASE   n_ps == 1" << endl;
      if (dddat[i_ddd].bin[0] != -1){
	int mirror_bin = dddat[i_ddd].mirror_bin[0];
	int min = 0;
	int max = 0;
	if ((dddat[i_ddd].distribution_1).typeCumulative != CUMULATIVE_NONE || (dddat[i_ddd].distribution_2).typeCumulative != CUMULATIVE_NONE){
	  // to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
	  ///		      if (dddat[i_ddd].xdistribution_number > 100){min = distribution_no_bin[i_d][0]; max = distribution_no_bin_max[i_d][0];}
	  ///		      else if (dddat[i_ddd].xdistribution_number > 0){min = distribution_no_bin[i_d][0]; max = distribution_no_bin[i_d][0] + 1;}
	}
	else if ((dddat[i_ddd].distribution_1).typeCumulative == CUMULATIVE_NONE && (dddat[i_ddd].distribution_2).typeCumulative == CUMULATIVE_NONE){
	  logger << LOG_DEBUG_VERBOSE << "((dddat[i_ddd].distribution_1).typeCumulative == CUMULATIVE_NONE && (dddat[i_ddd].distribution_2).typeCumulative == CUMULATIVE_NONE)" << endl;
	  logger << LOG_DEBUG_VERBOSE << "distribution_no_bin[" << i_ddd << "].size() = " << dddat[i_ddd].distribution_no_bin.size() << endl;
	  min = dddat[i_ddd].distribution_no_bin[0]; 
	  max = dddat[i_ddd].distribution_no_bin[0] + 1;
	}
	
	if (csi->class_contribution_IRcut_implicit){
	  for (int i_q = 0; i_q <= cut_ps[0]; i_q++){
	    for (int i_b = min; i_b < max; i_b++){ // was '<=' before !!!

	      bin_count_TSV[i_q][i_d][i_b]++;

	      for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
		if (!switch_distribution_TSV[i_s]){continue;}
		for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		    // i_b -> only one bin: dddat[i_ddd].bin[0]
		    // i_q -> all qTcut's:   result simulataneously calculated for all values
		    if (dddat[i_ddd].mirror_type < 1){
		      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][0][0];
		      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += pow(ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][0][0], 2);
		    }
		    else {
		      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][1][0];
		      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += pow(ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][1][0], 2);
		      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][mirror_bin] += ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][2][0];
		      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][mirror_bin] += pow(ps_integrand_qTcut_TSV[i_q][i_s][i_r][i_f][2][0], 2);
		    }
		  }
		}
	      }
	    }
	  }
	}
	else { //  !csi->class_contribution_IRcut_implicit
	  for (int i_q = 0; i_q <= cut_ps[0]; i_q++){
	    for (int i_b = min; i_b < max; i_b++){ // was '<=' before !!!

	      bin_count_TSV[i_q][i_d][i_b]++;

	      for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
		if (!switch_distribution_TSV[i_s]){continue;}
		for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		    if (dddat[i_ddd].mirror_type < 1){
		      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += ps_integrand_TSV[i_s][i_r][i_f][0][0];
		      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += pow(ps_integrand_TSV[i_s][i_r][i_f][0][0], 2);
		    }
		    else {
		      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += ps_integrand_TSV[i_s][i_r][i_f][1][0];
		      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][i_b] += pow(ps_integrand_TSV[i_s][i_r][i_f][1][0], 2);
		      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][mirror_bin] += ps_integrand_TSV[i_s][i_r][i_f][2][0];
		      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][mirror_bin] += pow(ps_integrand_TSV[i_s][i_r][i_f][2][0], 2);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    else if (n_ps > 1){
      /*      
      for (int i_ddd = 0; i_ddd < dddat.size(); i_ddd++){
	int i_d = dat.size() + i_ddd;
	logger << LOG_DEBUG_VERBOSE << "n_ps = " << n_ps << ":   dddat[" << i_ddd << " - " << i_ddd << "] : " << dddat[i_ddd].name << endl;

	if ((dddat[i_ddd].distribution_1).typeCumulative != CUMULATIVE_NONE || (dddat[i_ddd].distribution_2).typeCumulative != CUMULATIVE_NONE){
	  // to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
	  
	}
	else if ((dddat[i_ddd].distribution_1).typeCumulative == CUMULATIVE_NONE || (dddat[i_ddd].distribution_2).typeCumulative == CUMULATIVE_NONE){
	  logger << LOG_DEBUG_VERBOSE << "dddat[i_ddd].distribution_no_bin.size() = " << dddat[i_ddd].distribution_no_bin.size() << endl;
	  if (dddat[i_ddd].distribution_no_bin.size() != 0){
	    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
	      if (!switch_distribution_TSV[i_s]){continue;}
	      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		  logger << LOG_DEBUG_VERBOSE << "change_qTcut_bin_weight_TSV[" << i_s << "][" << i_r << "][" << i_f << "].size() = " << change_qTcut_bin_weight_TSV[i_s][i_r][i_f].size() << endl;
		  logger << LOG_DEBUG_VERBOSE << "change_qTcut_bin_weight2_TSV[" << i_s << "][" << i_r << "][" << i_f << "].size() = " << change_qTcut_bin_weight2_TSV[i_s][i_r][i_f].size() << endl;
		  logger << LOG_DEBUG_VERBOSE << "change_qTcut_bin_weight_2nd_TSV[" << i_s << "][" << i_r << "][" << i_f << "].size() = " << change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f].size() << endl;
		  logger << LOG_DEBUG_VERBOSE << "change_qTcut_bin_weight2_2nd_TSV[" << i_s << "][" << i_r << "][" << i_f << "].size() = " << change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f].size() << endl;

		  
		  change_qTcut_bin_weight_TSV[i_s][i_r][i_f].clear();
	      logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
		  change_qTcut_bin_weight2_TSV[i_s][i_r][i_f].clear();
	      logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
		  change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f].clear();
	      logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << endl;
		  change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f].clear();
		  logger << LOG_DEBUG_VERBOSE << "change_qTcut_bin_weight_TSV[" << i_s << "][" << i_r << "][" << i_f << "].size() = " << change_qTcut_bin_weight_TSV[i_s][i_r][i_f].size() << endl;
		  		  
		  change_qTcut_bin_weight_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dddat[i_ddd].distribution_no_bin.size() - 1, 0.));
		  change_qTcut_bin_weight2_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dddat[i_ddd].distribution_no_bin.size() - 1, 0.));
		  change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dddat[i_ddd].distribution_no_bin.size() - 1, 0.));
		  change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f].resize(distribution_no_qTcut.size(), vector<double> (dddat[i_ddd].distribution_no_bin.size() - 1, 0.));
		  logger << LOG_DEBUG_VERBOSE << "change_qTcut_bin_weight_TSV[" << i_s << "][" << i_r << "][" << i_f << "].size() = " << change_qTcut_bin_weight_TSV[i_s][i_r][i_f].size() << endl;
		}
	      }
	    }
	  }
	  logger << LOG_DEBUG_VERBOSE << "dddat[i_ddd].distribution_no_bin.size() = " << dddat[i_ddd].distribution_no_bin.size() << endl;
	  
	  // 1st: asymmetric (no symmetry wrt. rotating the positive beam axis to the negative one and vice versa) distribution
	  // 2nd: identical initial-state partons -> no contribution from rotated point (double-counting)
		    logger << LOG_DEBUG_VERBOSE << "i_ddd = " << i_ddd << endl;
	  for (int i_b = 0; i_b < dddat[i_ddd].distribution_no_bin.size(); i_b++){
	    // start with weight from higher qTcut values:
	    for (int i_q = distribution_no_qTcut.size() - 1; i_q >= 0; i_q--){
	      
	      logger << LOG_DEBUG_VERBOSE << "i_q = " << i_q << endl;
	      for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
		if (!switch_distribution_TSV[i_s]){continue;}
		for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		    if (i_q < distribution_no_qTcut.size() - 1){
		      logger << LOG_DEBUG_VERBOSE << "before 1" << endl;
		      change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q + 1][i_b];
		      change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][i_q + 1][i_b];
		      // new:
		      // mirror (2nd) contribution ??? !!!
		      logger << LOG_DEBUG_VERBOSE << "after 1" << endl;
		    }
		    
		    if (dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q].size() > 0){
		      logger << LOG_DEBUG_VERBOSE << "before 2" << endl;
		      for (int i_a = 0; i_a < dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q].size(); i_a++){
			if (dddat[i_ddd].mirror_type < 1){
			  // check if only [0] component of type_parton is relevant !!!
			  //			  if (csi->type_parton[0][1] == csi->type_parton[0][2] || ((dddat[i_ddd].distribution_1).symm == 0 && (dddat[i_ddd].distribution_2).symm == 0)){
			  //			  if (dat[i_ddd].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2]){
			  //			  logger << LOG_DEBUG_VERBOSE << "distribution_no_bin_no_qTcut_phasespace[" << i_ddd << "][" << i_b << "][" << i_q << "][" << i_q << "] = " << dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a] << endl;
			  change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][0][dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			}
			else {
			  change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][1][dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			  change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][2][dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			}
		      }
		      logger << LOG_DEBUG_VERBOSE << "after 2" << endl;

		      change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][i_q][i_b] = pow(change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b], 2);
		      logger << LOG_DEBUG_VERBOSE << "before 3" << endl;
		      if (dddat[i_ddd].mirror_type == 1){
			//if (!(csi->type_parton[0][1] == csi->type_parton[0][2] || ((dddat[i_ddd].distribution_1).symm == 0 && (dddat[i_ddd].distribution_2).symm == 0))){
			//			if (!(dat[i_ddd].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
			change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][i_q][i_b] = pow(change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q][i_b], 2);
		      }
		      logger << LOG_DEBUG_VERBOSE << "after 3" << endl;

		    }
//		    else if (change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] == 0.){} // What about "2nd" bin ???
		    else {
		      logger << LOG_DEBUG_VERBOSE << "before 4" << endl;
		      // Why actually not ??? !!!
		      // Should actually happen when only other bins change for this qTcut value !!!
		      // -> distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]  could be adapted correspondingly (
		      //			logger << LOG_WARN << "dddistribution   ---   n_ps != 1   ---   dddat[" << i_ddd << "].typeCumulative == CUMULATIVE_NONE   ---   Should not happen now!" << endl;
		      //			logger << LOG_WARN << "distribution_no_bin_no_qTcut_phasespace[" << i_ddd << "][" << i_b << "][" << i_q << "].size() = " << dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q].size() << endl;
		      if (i_q < distribution_no_qTcut.size() - 1){// new !!! check if anything changes ...
			change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight_TSV[i_s][i_r][i_f][i_q + 1][i_b];
			change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][i_q + 1][i_b];
			if (dddat[i_ddd].mirror_type == 1){
			  //if (!(csi->type_parton[0][1] == csi->type_parton[0][2] || ((dddat[i_ddd].distribution_1).symm == 0 && (dddat[i_ddd].distribution_2).symm == 0))){
			  //			if (!(dat[i_ddd].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
			  change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][i_q + 1][i_b];
			  change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][i_q][i_b] = change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][i_q + 1][i_b];
			}
		      }
		      logger << LOG_DEBUG_VERBOSE << "after 4" << endl;
		    }
		  }
		}
	      }
	      logger << LOG_DEBUG_VERBOSE << "i_q = " << i_q << endl;
	    }
	  }
	
	  logger << LOG_DEBUG_VERBOSE << "before 5" << endl;
		      
	  for (int i_b = 0; i_b < dddat[i_ddd].distribution_no_bin.size(); i_b++){
	    logger << LOG_DEBUG_VERBOSE << " dddat[" << i_ddd << "].distribution_no_bin[" << i_b << "] = " <<  dddat[i_ddd].distribution_no_bin[i_b] << endl;

	    int mirror_bin = -1;
	    if (dddat[i_ddd].mirror_type == 1){
	      //		      if (!(csi->type_parton[0][1] == csi->type_parton[0][2] || ((dddat[i_ddd].distribution_1).symm == 0 && (dddat[i_ddd].distribution_2).symm == 0))){
	      // introduce mirror bins (needed for symmetric disributions) !!!
	      if (((dddat[i_ddd].distribution_1).symm == 0) && ((dddat[i_ddd].distribution_2).symm == 0)){ // asymmetric distribution 1, asymmetric distribution 2
		mirror_bin = dddat[i_ddd].distribution_no_bin[i_b];
	      }
	      else if (((dddat[i_ddd].distribution_1).symm == 0) && ((dddat[i_ddd].distribution_2).symm == 1)){ // asymmetric distribution 1, symmetric distribution 2
		// check if correct !!! - doesn't seem to be ... 
		mirror_bin = bin[dddat[i_ddd].d_1][0] * (dddat[i_ddd].distribution_2).n_bins + ((dddat[i_ddd].distribution_2).n_bins - 1 - bin[dddat[i_ddd].d_2][0]);
	      }
	      else if (((dddat[i_ddd].distribution_1).symm == 1) && ((dddat[i_ddd].distribution_2).symm == 0)){ // symmetric distribution 1, asymmetric distribution 2
		// check if correct !!! - doesn't seem to be ... 
		mirror_bin = ((dddat[i_ddd].distribution_1).n_bins - 1 - bin[dddat[i_ddd].d_1][0]) * (dddat[i_ddd].distribution_2).n_bins + bin[dddat[i_ddd].d_2][0];
	      }
	      else if (((dddat[i_ddd].distribution_1).symm == 1) && ((dddat[i_ddd].distribution_2).symm == 1)){ // symmetric distribution 1, symmetric distribution 2
		mirror_bin = dddat[i_ddd].n_bins - 1 - dddat[i_ddd].distribution_no_bin[i_b];
		// mirror_bin ?????
	      }
	    }
	    int j_q = 0;
	    for (int i_q = 0; i_q < output_n_qTcut; i_q++){

	      bin_count_TSV[i_q][i_d][dddat[i_ddd].distribution_no_bin[i_b]]++;

	      for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
		if (!switch_distribution_TSV[i_s]){continue;}
		for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		  for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		    bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][dddat[i_ddd].distribution_no_bin[i_b]] += change_qTcut_bin_weight_TSV[i_s][i_r][i_f][j_q][i_b];
		    bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][dddat[i_ddd].distribution_no_bin[i_b]] += change_qTcut_bin_weight2_TSV[i_s][i_r][i_f][j_q][i_b];
		    if (dddat[i_ddd].mirror_type == 1){
		      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][mirror_bin] += change_qTcut_bin_weight_2nd_TSV[i_s][i_r][i_f][j_q][i_b];
		      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][mirror_bin] += change_qTcut_bin_weight2_2nd_TSV[i_s][i_r][i_f][j_q][i_b];
		    }
		  }
		}
	      }
	      if (i_q == distribution_no_qTcut[j_q]){j_q++;}
	      if (j_q == distribution_no_qTcut.size()){break;}
	    }
	    logger << LOG_DEBUG_VERBOSE << "break   j_q = " << j_q << endl;		    
	  }
	}
      }
*/
    }
  }


  if (switch_distribution_at_all_TSV == 1){
    // temporary !!!
    int x_s = no_reference_TSV;
    int x_r = (n_scale_ren_TSV[x_s] - 1) / 2;
    int x_f = (n_scale_fact_TSV[x_s] - 1) / 2;

    for (int i_s = 0; i_s < n_extended_set_TSV; i_s++){
      if (switch_distribution_TSV[i_s]){
	for (int i_ddd = 0; i_ddd < dddat.size(); i_ddd++){
	  int i_d = dat.size() + i_ddd;
	  logger << LOG_DEBUG_VERBOSE << "n_ps = " << n_ps << ":   dddat[" << i_ddd << " - " << i_ddd << "] : " << dddat[i_ddd].name << endl;
  	  if (n_ps == 1){
	  }
  
	  else if (n_ps > 1){
	    
	    if ((dddat[i_ddd].distribution_1).typeCumulative != CUMULATIVE_NONE || (dddat[i_ddd].distribution_2).typeCumulative != CUMULATIVE_NONE){
	      // to be filled !!! dddistributions don't seem to be implemented for cumulative distribution

	    }
	    else if ((dddat[i_ddd].distribution_1).typeCumulative == CUMULATIVE_NONE || (dddat[i_ddd].distribution_2).typeCumulative == CUMULATIVE_NONE){
	      vector<vector<vector<vector<double> > > > change_bin_weight(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (distribution_no_qTcut.size(), vector<double> (dddat[i_ddd].distribution_no_bin.size(), 0.)))); // could be done statically via oset (with maximum length n_ps)
	      vector<vector<vector<vector<double> > > > change_bin_weight2(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (distribution_no_qTcut.size(), vector<double> (dddat[i_ddd].distribution_no_bin.size(), 0.)))); // could be done statically via oset (with maximum length n_ps)
	      vector<vector<vector<vector<double> > > > change_bin_weight_2nd(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (distribution_no_qTcut.size(), vector<double> (dddat[i_ddd].distribution_no_bin.size(), 0.)))); // could be done statically via oset (with maximum length n_ps)
	      vector<vector<vector<vector<double> > > > change_bin_weight2_2nd(n_scale_ren_TSV[i_s], vector<vector<vector<double> > > (n_scale_fact_TSV[i_s], vector<vector<double> > (distribution_no_qTcut.size(), vector<double> (dddat[i_ddd].distribution_no_bin.size(), 0.)))); // could be done statically via oset (with maximum length n_ps)

	      for (int i_r = 0; i_r < n_scale_ren_TSV[i_s]; i_r++){
		for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
		  // 1st: asymmetric (no symmetry wrt. rotating the positive beam axis to the negative one and vice versa) distribution
		  // 2nd: identical initial-state partons -> no contribution from rotated point (double-counting)
		  for (int i_b = 0; i_b < dddat[i_ddd].distribution_no_bin.size(); i_b++){
		    // start with weight from higher qTcut values:
		    for (int i_q = distribution_no_qTcut.size() - 1; i_q >= 0; i_q--){
		      if (i_q < distribution_no_qTcut.size() - 1){
			change_bin_weight[i_r][i_f][i_q][i_b] = change_bin_weight[i_r][i_f][i_q + 1][i_b];
			change_bin_weight_2nd[i_r][i_f][i_q][i_b] = change_bin_weight_2nd[i_r][i_f][i_q + 1][i_b];
		      }

		      if (dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q].size() > 0){
			for (int i_a = 0; i_a < dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q].size(); i_a++){
			  if (dddat[i_ddd].mirror_type < 1){
			    // check if only [0] component of type_parton is relevant !!!
			    //			  if (csi->type_parton[0][1] == csi->type_parton[0][2] || ((dddat[i_ddd].distribution_1).symm == 0 && (dddat[i_ddd].distribution_2).symm == 0)){
			    //			  if (dat[i_ddd].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2]){
			    logger << LOG_DEBUG_VERBOSE << "distribution_no_bin_no_qTcut_phasespace[" << i_ddd << "][" << i_b << "][" << i_q << "][" << i_q << "] = " << dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a] << endl;
			    change_bin_weight[i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][0][dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			  }
			  else {
			    change_bin_weight[i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][1][dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			    change_bin_weight_2nd[i_r][i_f][i_q][i_b] += ps_integrand_TSV[i_s][i_r][i_f][2][dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]];
			  }
			}
			change_bin_weight2[i_r][i_f][i_q][i_b] = pow(change_bin_weight[i_r][i_f][i_q][i_b], 2);
			if (dddat[i_ddd].mirror_type == 1){
			  //if (!(csi->type_parton[0][1] == csi->type_parton[0][2] || ((dddat[i_ddd].distribution_1).symm == 0 && (dddat[i_ddd].distribution_2).symm == 0))){
			  //			if (!(dat[i_ddd].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
			  change_bin_weight2_2nd[i_r][i_f][i_q][i_b] = pow(change_bin_weight_2nd[i_r][i_f][i_q][i_b], 2);
			}
		      }
		      //		      else if (change_bin_weight[i_r][i_f][i_q][i_b] == 0.){}
		      else {
			// Why actually not ??? !!!
			// Should actually happen when only other bins change for this qTcut value !!!
			// -> distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_a]  could be adapted correspondingly (
			//			logger << LOG_WARN << "dddistribution   ---   n_ps != 1   ---   dddat[" << i_ddd << "].typeCumulative == CUMULATIVE_NONE   ---   Should not happen now!" << endl;
			//			logger << LOG_WARN << "distribution_no_bin_no_qTcut_phasespace[" << i_ddd << "][" << i_b << "][" << i_q << "].size() = " << dddat[i_ddd].distribution_no_bin_no_qTcut_phasespace[i_b][i_q].size() << endl;
			if (i_q < distribution_no_qTcut.size() - 1){
			  change_bin_weight[i_r][i_f][i_q][i_b] = change_bin_weight[i_r][i_f][i_q + 1][i_b];
			  change_bin_weight2[i_r][i_f][i_q][i_b] = change_bin_weight2[i_r][i_f][i_q + 1][i_b];
			  if (dddat[i_ddd].mirror_type == 1){
			    //if (!(csi->type_parton[0][1] == csi->type_parton[0][2] || ((dddat[i_ddd].distribution_1).symm == 0 && (dddat[i_ddd].distribution_2).symm == 0))){
			    //			if (!(dat[i_ddd].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
			    change_bin_weight_2nd[i_r][i_f][i_q][i_b] = change_bin_weight_2nd[i_r][i_f][i_q + 1][i_b];
			    change_bin_weight2_2nd[i_r][i_f][i_q][i_b] = change_bin_weight2_2nd[i_r][i_f][i_q + 1][i_b];
			  }
			}
		      }
		    }
		  }
		  for (int i_b = 0; i_b < dddat[i_ddd].distribution_no_bin.size(); i_b++){
		    int j_q = 0;

		    for (int i_q = 0; i_q < output_n_qTcut; i_q++){
		      if (i_d == 0 && i_r == 0 && i_f == 0){
			logger << LOG_DEBUG_VERBOSE << "i_q = " << setw(3) << i_q << "   change_bin_weight[" << i_r << "][" << i_f << "][j_q = " << j_q << "][" << i_b << "] = " << change_bin_weight[i_r][i_f][j_q][i_b] << endl;
		      }
		    
		      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][dddat[i_ddd].distribution_no_bin[i_b]] += change_bin_weight[i_r][i_f][j_q][i_b];
		      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][dddat[i_ddd].distribution_no_bin[i_b]] += change_bin_weight2[i_r][i_f][j_q][i_b];

		      if (i_s == x_s && i_r == x_r && i_f == x_f){
			bin_count_TSV[i_q][i_d][dddat[i_ddd].distribution_no_bin[i_b]]++;
		      }

		      //		    logger << LOG_DEBUG_VERBOSE << "dat[" << i_ddd << "].symm = " << dat[i_ddd].symm << endl;
		      //		    if (!(dat[i_ddd].symm == 0 || csi->type_parton[0][1] == csi->type_parton[0][2])){
		      if (!(csi->type_parton[0][1] == csi->type_parton[0][2] || ((dddat[i_ddd].distribution_1).symm == 0 && (dddat[i_ddd].distribution_2).symm == 0))){
			// introduce mirror bins (needed for symmetric disributions) !!!
			int mirror_bin = -1;
			if (((dddat[i_ddd].distribution_1).symm == 0) && ((dddat[i_ddd].distribution_2).symm == 0)){ // asymmetric distribution 1, asymmetric distribution 2
			  mirror_bin = dddat[i_ddd].distribution_no_bin[i_b];
			}
			else if (((dddat[i_ddd].distribution_1).symm == 0) && ((dddat[i_ddd].distribution_2).symm == 1)){ // asymmetric distribution 1, symmetric distribution 2
			  // check if correct !!! - doesn't seem to be ... 
			  mirror_bin = bin[dddat[i_ddd].d_1][0] * (dddat[i_ddd].distribution_2).n_bins + ((dddat[i_ddd].distribution_2).n_bins - 1 - bin[dddat[i_ddd].d_2][0]);
			}
			else if (((dddat[i_ddd].distribution_1).symm == 1) && ((dddat[i_ddd].distribution_2).symm == 0)){ // symmetric distribution 1, asymmetric distribution 2
			  // check if correct !!! - doesn't seem to be ... 
			  mirror_bin = ((dddat[i_ddd].distribution_1).n_bins - 1 - bin[dddat[i_ddd].d_1][0]) * (dddat[i_ddd].distribution_2).n_bins + bin[dddat[i_ddd].d_2][0];
			}
			else if (((dddat[i_ddd].distribution_1).symm == 1) && ((dddat[i_ddd].distribution_2).symm == 1)){ // symmetric distribution 1, symmetric distribution 2
			  mirror_bin = dddat[i_ddd].n_bins - 1 - dddat[i_ddd].distribution_no_bin[i_b];
			  // mirror_bin ?????
			}
			
			bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][mirror_bin] += change_bin_weight_2nd[i_r][i_f][j_q][i_b];
			bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][mirror_bin] += change_bin_weight2_2nd[i_r][i_f][j_q][i_b];
			//		      bin_weight_TSV[i_s][i_q][i_r][i_f][i_d][dddat[i_ddd].n_bins - 1 - dddat[i_ddd].distribution_no_bin[i_b]] += change_bin_weight_2nd[i_r][i_f][j_q][i_b];
			//		      bin_weight2_TSV[i_s][i_q][i_r][i_f][i_d][dddat[i_ddd].n_bins - 1 - dddat[i_ddd].distribution_no_bin[i_b]] += change_bin_weight2_2nd[i_r][i_f][j_q][i_b];
		      }
		      if (i_q == distribution_no_qTcut[j_q]){j_q++;}
		      if (j_q == distribution_no_qTcut.size()){break;}
		    }
		    logger << LOG_DEBUG_VERBOSE << "break   j_q = " << j_q << endl;		    
		  }
		}
	      }
	    }
	    else {
	      logger << LOG_FATAL << "n_ps != 1:   dddat[" << i_ddd << "]   Wrong typeCumulative" << endl; exit(1);
	    }
	    
	  }
	  else {
	    logger << LOG_FATAL << "Wrong n_ps value (< 1)" << endl; exit(1);
	  }
	  logger << LOG_DEBUG_VERBOSE << "2nd step   dddat[" << i_ddd << "] = " << dddat[i_ddd].name << "   done" << endl;
	}
	logger << LOG_DEBUG_VERBOSE << "i_s = " << i_s << " done" << endl;
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
