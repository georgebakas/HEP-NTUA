#include "../include/classes.cxx"
dddistribution::dddistribution(){
  name = "";
  d_1 = 0;
  d_2 = 0;
  n_bins = 0;
  step = 0.;
  xdistribution temp;
  distribution_1 = temp;
  distribution_2 = temp;
  mirror_type = -1;
}
dddistribution::dddistribution(string _name, int _d_1, int _d_2, observable_set *_oset){
  Logger logger("dddistribution::dddistribution");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  //, vector<xdistribution> _dat
  oset = _oset;

  name = _name;
  d_1 = _d_1;
  d_2 = _d_2;
  n_bins = oset->dat[d_1].n_bins * oset->dat[d_2].n_bins;
  step = oset->dat[d_1].step * oset->dat[d_2].step;
  distribution_1 = oset->dat[d_1];
  distribution_2 = oset->dat[d_2];
  // required to use bins etc. from distributions...
  //  distribution_1 = &oset->dat[d_1];
  //  distribution_2 = &oset->dat[d_2];

  logger << LOG_DEBUG_VERBOSE << "before mirror_type" << endl;
  logger << LOG_DEBUG_VERBOSE << "oset->csi->type_parton.size() = " << oset->csi->type_parton.size() << endl;
  logger << LOG_DEBUG_VERBOSE << "oset->csi->type_parton[0].size() = " << oset->csi->type_parton[0].size() << endl;
  if (oset->csi->type_parton[0].size() > 2){
      if (oset->csi->type_parton[0][1] != oset->csi->type_parton[0][2]){mirror_type = -1;}
      else if (distribution_1.symm == 0 && distribution_2.symm == 0){mirror_type = 0;}
      else if (distribution_1.symm != 0 || distribution_2.symm != 0){mirror_type = 1;}
      else {logger << LOG_ERROR << "dddistribution  " << name << " : Invalid symmtetry structure!" << endl; exit(1);}
  }
  logger << LOG_DEBUG_VERBOSE << "mirror_type = " << mirror_type << endl;

  bin.resize(oset->n_ps);
  bin_max.resize(oset->n_ps);
   
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void dddistribution::determineBin() {
  Logger logger("dddistribution::determineBin");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;
  
  /*
  observable.resize(oset->n_ps);
  observable_max.resize(oset->n_ps);
  */
  bin.resize(oset->n_ps);
  bin_max.resize(oset->n_ps);

  mirror_bin.resize(oset->n_ps);
  mirror_bin_max.resize(oset->n_ps);

  /*
  if (oset->n_ps == 1){
    if (oset->dat[d_1].bin[0] == -1 || oset->dat[d_2].bin[0] == -1){
      bin[0] = -1;
    }
    else {
      bin[0] = oset->dat[d_1].bin[0] * distribution_2.n_bins + oset->dat[d_2].bin[0];
    }

    if (oset->csi->type_parton[0][1] != oset->csi->type_parton[0][2]){
      if (distribution_1.symm == 0 && distribution_2.symm == 0){ //  asymmetric distribution 1, asymmetric distribution 2
	mirror_bin[0] = bin[0];
      }
      else if (distribution_1.symm == 0 && distribution_2.symm == 1){ //  asymmetric distribution 1, symmetric distribution 2
	mirror_bin[0] = oset->dat[d_1].bin[0] * distribution_2.n_bins + (distribution_2.n_bins - 1 - oset->dat[d_2].bin[0]);
      }
      else if (distribution_1.symm == 1 && distribution_2.symm == 0){ //  symmetric distribution 1, asymmetric distribution 2
	mirror_bin[0] = (distribution_1.n_bins - 1 - oset->dat[d_1].bin[0]) * distribution_2.n_bins + oset->dat[d_2].bin[0];
      }
      else if (distribution_1.symm == 1 && distribution_2.symm == 1){ //  symmetric distribution 1, symmetric distribution 2
	mirror_bin[0] = n_bins - 1 - bin[0];
      }
      else {logger << LOG_FATAL << "No allowed case chosen." << endl; exit(1);}
    }
  }
  else if (oset->n_ps > 1){  // more than one phase-space (n_ps > 1)
    for (int i_a = 0; i_a < oset->n_ps; i_a++){
      if (oset->dat[d_1].bin[i_a] == -1 || oset->dat[d_2].bin[i_a] == -1){
	bin[i_a] = -1;
      }
      else {
      bin[i_a] = oset->dat[d_1].bin[i_a] * distribution_2.n_bins + oset->dat[d_2].bin[i_a];
      }
      
      if (oset->csi->type_parton[i_a][1] != oset->csi->type_parton[i_a][2]){
	if (distribution_1.symm == 0 && distribution_2.symm == 0){ //  asymmetric distribution 1, asymmetric distribution 2
	  mirror_bin[i_a] = bin[i_a];
	}
	else if (distribution_1.symm == 0 && distribution_2.symm == 1){ //  asymmetric distribution 1, symmetric distribution 2
	  mirror_bin[i_a] = oset->dat[d_1].bin[i_a] * distribution_2.n_bins + (distribution_2.n_bins - 1 - oset->dat[d_2].bin[i_a]);
	}
	else if (distribution_1.symm == 1 && distribution_2.symm == 0){ //  symmetric distribution 1, asymmetric distribution 2
	  mirror_bin[i_a] = (distribution_1.n_bins - 1 - oset->dat[d_1].bin[i_a]) * distribution_2.n_bins + oset->dat[d_2].bin[i_a];
	}
	else if (distribution_1.symm == 1 && distribution_2.symm == 1){ //  symmetric distribution 1, symmetric distribution 2
	  mirror_bin[i_a] = n_bins - 1 - bin[i_a];
	}
	else {logger << LOG_FATAL << "No allowed case chosen." << endl; exit(1);}
      }
    }
  }
  */
  // Actually no need to split...
  logger << LOG_DEBUG_VERBOSE << "distribution_1.name = " << distribution_1.xdistribution_name << endl;
  logger << LOG_DEBUG_VERBOSE << "distribution_2.name = " << distribution_2.xdistribution_name << endl;
  //    logger << LOG_DEBUG_VERBOSE << "distribution_1.bin.size() = " << distribution_1.bin.size() << endl;
  for (int i_a = 0; i_a < oset->n_ps; i_a++){
    logger << LOG_DEBUG_VERBOSE << "distribution_1.bin.size() = " << distribution_1.bin.size() << endl;
    logger << LOG_DEBUG_VERBOSE << "distribution_2.bin.size() = " << distribution_2.bin.size() << endl;
    logger << LOG_DEBUG_VERBOSE << "oset->dat[" << d_1 << "].bin.size() = " << oset->dat[d_1].bin.size() << endl;
    logger << LOG_DEBUG_VERBOSE << "oset->dat[" << d_2 << "].bin.size() = " << oset->dat[d_2].bin.size() << endl;
    if (oset->dat[d_1].bin[i_a] == -1 || oset->dat[d_2].bin[i_a] == -1){
      bin[i_a] = -1;
    }
    else {
      bin[i_a] = oset->dat[d_1].bin[i_a] * distribution_2.n_bins + oset->dat[d_2].bin[i_a];
      logger << LOG_DEBUG_VERBOSE << "distribution_2.n_bins = " << distribution_2.n_bins << endl;
    }
    
    if (mirror_type == -1){ //  actually not mirror_bin needed !!!
      mirror_bin[i_a] = bin[i_a];
    }
    else if (mirror_type == 0){
      //    if (oset->csi->type_parton[i_a][1] != oset->csi->type_parton[i_a][2]){
      //      if (distribution_1.symm == 0 && distribution_2.symm == 0){ //  asymmetric distribution 1, asymmetric distribution 2
      mirror_bin[i_a] = bin[i_a];
    }
    else if (mirror_type == 1){
      //      else if (distribution_1.symm == 0 && distribution_2.symm == 1){ //  asymmetric distribution 1, symmetric distribution 2
      if (distribution_1.symm == 0 && distribution_2.symm == 1){ //  asymmetric distribution 1, symmetric distribution 2
	mirror_bin[i_a] = oset->dat[d_1].bin[i_a] * distribution_2.n_bins + (distribution_2.n_bins - 1 - oset->dat[d_2].bin[i_a]);
      }
      else if (distribution_1.symm == 1 && distribution_2.symm == 0){ //  symmetric distribution 1, asymmetric distribution 2
	mirror_bin[i_a] = (distribution_1.n_bins - 1 - oset->dat[d_1].bin[i_a]) * distribution_2.n_bins + oset->dat[d_2].bin[i_a];
      }
      else if (distribution_1.symm == 1 && distribution_2.symm == 1){ //  symmetric distribution 1, symmetric distribution 2
	mirror_bin[i_a] = n_bins - 1 - bin[i_a];
      }
      else {logger << LOG_FATAL << "No allowed case chosen." << endl; exit(1);}
    }
    else {logger << LOG_ERROR << "dddistribution  " << name << " : Invalid symmtetry structure!" << endl; exit(1);}
    
    logger << LOG_DEBUG_VERBOSE << setw(20) << name << "   bin[" << i_a << "] = " << bin[i_a] << "   mirror_bin[" << i_a << "] = " << mirror_bin[i_a] << endl;

  }
   
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void dddistribution::initialization_distribution_bin(){
  static Logger logger("dddistribution::initialization_distribution_bin");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int i_q = 0; i_q < oset->distribution_no_qTcut.size(); i_q++){
    stringstream out_dnq;
    for (int i_qa = 0; i_qa < oset->distribution_no_qTcut_phasespace[i_q].size(); i_qa++){
      out_dnq << oset->distribution_no_qTcut_phasespace[i_q][i_qa] << " ";
    }
    logger << LOG_DEBUG << "oset->distribution_no_qTcut_phasespace [" << i_q << "] @ " << setw(3) << oset->distribution_no_qTcut[i_q] << " = " << out_dnq.str() << endl;
  }
  
  distribution_no_bin.clear();
  distribution_no_bin_phasespace.clear();
  distribution_no_bin_no_qTcut_phasespace.clear();
  
  if (oset->n_ps == 1){
    logger << LOG_DEBUG_VERBOSE << name << "   oset->n_ps = 1" << endl;
    /////    int i_d = dat.size() + i_ddd;
    if (distribution_1.typeCumulative != CUMULATIVE_NONE && distribution_2.typeCumulative != CUMULATIVE_NONE){
      // to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
    }
    else if (distribution_1.typeCumulative == CUMULATIVE_NONE && distribution_2.typeCumulative != CUMULATIVE_NONE){
      // to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
    }
    else if (distribution_1.typeCumulative != CUMULATIVE_NONE && distribution_2.typeCumulative == CUMULATIVE_NONE){
      // to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
    }
    else if (distribution_1.typeCumulative == CUMULATIVE_NONE && distribution_2.typeCumulative == CUMULATIVE_NONE){
      logger << LOG_DEBUG_VERBOSE << name << "   bin[0] = " << bin[0] << endl;
      if (bin[0] != -1){
	distribution_no_bin.resize(1, bin[0]);
	//	  distribution_no_bin = bin; // could still be -1 --- does the same ???
	distribution_no_bin_phasespace.resize(1, vector<int> (1, 0)); // -> (1, 0) only 1 phasespace, namely 0; resize(1, ...) -> only 1 bin
	distribution_no_bin_no_qTcut_phasespace.resize(1, vector<vector<int> > (1, vector<int> (1, 0)));
	
	/////	for (int i_q = 0; i_q <= cut_ps[0]; i_q++){bin_count_TSV[i_q][i_d][dddat[i_ddd].bin[0]]++;}// ??? < or <= ???  // mirrored contributions ???
      }
    }
    else {
      logger << LOG_FATAL << "dddat(" << name << "):   oset->n_ps = 1   typeCumulative combination not implemented." << endl; exit(1);
    }
  }
  else if (oset->n_ps > 1){
    if (distribution_1.typeCumulative != CUMULATIVE_NONE || distribution_2.typeCumulative != CUMULATIVE_NONE){
      // to be filled !!! dddistributions don't seem to be implemented for cumulative distribution
    }
    else if (distribution_1.typeCumulative == CUMULATIVE_NONE && distribution_2.typeCumulative == CUMULATIVE_NONE){
      ///////////////////////////////////////////////////////////////////////////////////////
      //  oset->n_ps > 1  //  (dddat[i_ddd].distribution_1&2).typeCumulative == CUMULATIVE_NONE  //
      ///////////////////////////////////////////////////////////////////////////////////////
      distribution_no_bin = bin;

      sort(distribution_no_bin.begin(), distribution_no_bin.end());
      for (int i_b = oset->n_ps - 1; i_b > 0; i_b--){
	if (distribution_no_bin[i_b] == distribution_no_bin[i_b - 1]){distribution_no_bin.erase(distribution_no_bin.begin() + i_b);}
      }
      if (distribution_no_bin[0] == -1){distribution_no_bin.erase(distribution_no_bin.begin());}
      logger << LOG_DEBUG_VERBOSE << "distribution_no_bin(" << name << ").size() = " << distribution_no_bin.size() << endl;
      
      stringstream out_dnb;
      for (int i_b = 0; i_b < distribution_no_bin.size(); i_b++){out_dnb << distribution_no_bin[i_b] << "   ";}
      logger << LOG_DEBUG << "distribution_no_bin (" << name << ") = " << out_dnb.str() << endl;
      
      distribution_no_bin_phasespace.resize(distribution_no_bin.size());
      for (int i_a = 0; i_a < oset->n_ps; i_a++){
	for (int i_b = 0; i_b < distribution_no_bin.size(); i_b++){
	  if (distribution_no_bin[i_b] == bin[i_a]){distribution_no_bin_phasespace[i_b].push_back(i_a); break;}
	}
      }


      for (int x_b = 0; x_b < distribution_no_bin_phasespace.size(); x_b++){
	stringstream out_dnbp;
	for (int x_a = 0; x_a < distribution_no_bin_phasespace[x_b].size(); x_a++){
	  out_dnbp << distribution_no_bin_phasespace[x_b][x_a] << " ";
	}
	logger << LOG_DEBUG << "distribution_no_bin_phasespace (" << name << ")[" << x_b << "] @ " << setw(3) << distribution_no_bin[x_b] << " = " << out_dnbp.str() << endl;
      }

      distribution_no_bin_no_qTcut_phasespace.resize(distribution_no_bin.size(), vector<vector<int> > (oset->distribution_no_qTcut.size()));
	    
      for (int i_b = 0; i_b < distribution_no_bin.size(); i_b++){
	for (int i_q = 0; i_q < oset->distribution_no_qTcut.size(); i_q++){
	  for (int i_ba = 0; i_ba < distribution_no_bin_phasespace[i_b].size(); i_ba++){
	    for (int i_qa = 0; i_qa < oset->distribution_no_qTcut_phasespace[i_q].size(); i_qa++){
	      //	      logger << LOG_DEBUG_VERBOSE << "i_b = " << i_b << "   i_q = " << i_q << "   i_ba = " << i_ba << "   i_qa = " << i_qa << endl;
	      if (distribution_no_bin_phasespace[i_b][i_ba] == oset->distribution_no_qTcut_phasespace[i_q][i_qa]){
		distribution_no_bin_no_qTcut_phasespace[i_b][i_q].push_back(distribution_no_bin_phasespace[i_b][i_ba]);
	      }
	    }
	  }
	}
      }

      for (int i_b = 0; i_b < distribution_no_bin.size(); i_b++){
	for (int i_q = 0; i_q < oset->distribution_no_qTcut.size(); i_q++){
	  stringstream out_dnbnqp;
	  for (int i_x = 0; i_x < distribution_no_bin_no_qTcut_phasespace[i_b][i_q].size(); i_x++){
	    out_dnbnqp << distribution_no_bin_no_qTcut_phasespace[i_b][i_q][i_x] << " ";
	  }
	  logger << LOG_DEBUG << "distribution_no_bin_no_qTcut_phasespace [" << i_b << "][" << i_q << "] @ bin -> " << setw(3) << distribution_no_bin[i_b] << " @ qTcut -> " << setw(3) << oset->distribution_no_qTcut[i_q] << " = " << out_dnbnqp.str() << endl;
	}
      }
    }
    else {
      logger << LOG_FATAL << "dddat(" << name << "):   typeCumulative combination not implemented." << endl; exit(1);
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}






ostream & operator << (ostream &s, const dddistribution &sp){
  static Logger logger("dddistribution");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  stringstream temp_ss;
  temp_ss << setw(20) << sp.name;
  temp_ss << "   n_bins = " << setw(3) << sp.n_bins;
  temp_ss << "   mirror_type = " << setw(3) << sp.mirror_type;
  temp_ss << "   step = " << setw(3) << sp.step << " (unused)";
  temp_ss << "   d1 = " << setw(20) << (sp.distribution_1).xdistribution_name << " (" << setw(2) << sp.d_1 << ")";
  temp_ss << sp.distribution_1;
  temp_ss << "   d2 = " << setw(20) << (sp.distribution_2).xdistribution_name << " (" << setw(2) << sp.d_2 << ")";
  temp_ss << sp.distribution_2;
  temp_ss << endl;

  s << temp_ss.str() << endl;
  logger << LOG_DEBUG_VERBOSE << temp_ss.str() << endl;
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  return s;
}
