#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

void observable_set::initialization_LHAPDF(){
  static Logger logger("initialization_LHAPDF");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  const int SUBSET = LHAPDFsubset;

  cout << "LHAPDFname.length() = " << LHAPDFname.length() << endl;
  if (LHAPDFname.length()>=6) {
    if (LHAPDFname.substr(LHAPDFname.size() - 6, 6) == ".LHpdf"){
      const string NAME = LHAPDFname.substr(0, LHAPDFname.size() - 6);
      LHAPDF::initPDFSet(NAME, LHAPDF::LHPDF, SUBSET);
      //    LHAPDF::initPDFSet(NAME, LHAPDF::LHPDF, SUBSET);
      cout << "XXX initialization done here!!! LHAPDFname.length()>=6 .LHpdf" << endl;
      return;
    }
  }
  if (LHAPDFname.length()>=7) {
    if (LHAPDFname.substr(LHAPDFname.size() - 7, 7) == ".LHgrid"){
      const string NAME = LHAPDFname.substr(0, LHAPDFname.size() - 7);
      LHAPDF::initPDFSet(NAME, LHAPDF::LHGRID, SUBSET);
      cout << "XXX initialization done here!!! LHAPDFname.length()>=7 .LHgrid" << endl;
      return;
    }
  }

  const string NAME = LHAPDFname;
  LHAPDF::initPDFSet(NAME, LHAPDF::LHGRID, SUBSET);
  logger << LOG_DEBUG_VERBOSE << "LHAPDF::initPDFSet(NAME = " << NAME << ", LHAPDF::LHGRID, SUBSET = " << SUBSET << ");" << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::determine_selection_content_pdf(){
  static Logger logger("observable_set::determine_selection_content_pdf");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  logger << LOG_DEBUG << "N_f_active = " << N_f_active << endl;
  for (int i_p = 0; i_p < pdf_selection.size(); i_p++){
    logger << LOG_DEBUG << "pdf_selection[" << i_p << "] = " << pdf_selection[i_p] << endl;
  }
  for (int i_d = 0; i_d < pdf_disable.size(); i_d++){
    logger << LOG_DEBUG << "pdf_disable[" << i_d << "] = " << pdf_disable[i_d] << endl;
  }
  vector<int> content_p(14);
  content_p[0] = 1;
  content_p[1] = 2;
  content_p[2] = 3;
  content_p[3] = 4;
  content_p[4] = 5;
  content_p[5] = 6;
  content_p[6] = -1;
  content_p[7] = -2;
  content_p[8] = -3;
  content_p[9] = -4;
  content_p[10] = -5;
  content_p[11] = -6;
  content_p[12] = 0;
  content_p[13] = 7;
  vector<int> content_Q(12);
  content_Q[0] = 1;
  content_Q[1] = 2;
  content_Q[2] = 3;
  content_Q[3] = 4;
  content_Q[4] = 5;
  content_Q[5] = 6;
  content_Q[6] = -1;
  content_Q[7] = -2;
  content_Q[8] = -3;
  content_Q[9] = -4;
  content_Q[10] = -5;
  content_Q[11] = -6;
  vector<int> content_q(6);
  content_q[0] = 1;
  content_q[1] = 2;
  content_q[2] = 3;
  content_q[3] = 4;
  content_q[4] = 5;
  content_q[5] = 6;
  vector<int> content_qx(6);
  content_qx[0] = -1;
  content_qx[1] = -2;
  content_qx[2] = -3;
  content_qx[3] = -4;
  content_qx[4] = -5;
  content_qx[5] = -6;
  vector<int> content_g(1);
  content_g[0] = 0;
  vector<int> content_a(1);
  content_a[0] = 7;
  vector<int> content_D(2);
  content_D[0] = 1;
  content_D[1] = -1;
  vector<int> content_U(2);
  content_U[0] = 2;
  content_U[1] = -2;
  vector<int> content_S(2);
  content_S[0] = 3;
  content_S[1] = -3;
  vector<int> content_C(2);
  content_C[0] = 4;
  content_C[1] = -4;
  vector<int> content_B(2);
  content_B[0] = 5;
  content_B[1] = -5;
  vector<int> content_T(2);
  content_T[0] = 6;
  content_T[1] = -6;
  vector<int> content_d(1);
  content_d[0] = 1;
  vector<int> content_u(1);
  content_u[0] = 2;
  vector<int> content_s(1);
  content_s[0] = 3;
  vector<int> content_c(1);
  content_c[0] = 4;
  vector<int> content_b(1);
  content_b[0] = 5;
  vector<int> content_t(1);
  content_t[0] = 6;
  vector<int> content_dx(1);
  content_dx[0] = -1;
  vector<int> content_ux(1);
  content_ux[0] = -2;
  vector<int> content_sx(1);
  content_sx[0] = -3;
  vector<int> content_cx(1);
  content_cx[0] = -4;
  vector<int> content_bx(1);
  content_bx[0] = -5;
  vector<int> content_tx(1);
  content_tx[0] = -6;

  vector<vector<int> > content_parton(24);
  content_parton[0] = content_p;
  content_parton[1] = content_Q;
  content_parton[2] = content_q;
  content_parton[3] = content_qx;
  content_parton[4] = content_g;
  content_parton[5] = content_a;
  content_parton[6] = content_D;
  content_parton[7] = content_U;
  content_parton[8] = content_S;
  content_parton[9] = content_C;
  content_parton[10] = content_B;
  content_parton[11] = content_T;
  content_parton[12] = content_d;
  content_parton[13] = content_u;
  content_parton[14] = content_s;
  content_parton[15] = content_c;
  content_parton[16] = content_b;
  content_parton[17] = content_t;
  content_parton[18] = content_dx;
  content_parton[19] = content_ux;
  content_parton[20] = content_sx;
  content_parton[21] = content_cx;
  content_parton[22] = content_bx;
  content_parton[23] = content_tx;

  // pdf_disable is a vector of all disabled particles, which may contain "+"'s that are resolved here:

  //  vector<string> temp_pdf_disable = pdf_disable;
  //  vector<string> temp_pdf_disable_split;
  // pdf_disable is a vector of all disabled particles, whichout "+"'s
  //  int counter = 0;
  //    int start = 0;

  vector<string> temp_pdf_disable;
  for (int i_d = 0; i_d < pdf_disable.size(); i_d++){
    string temp_split = "";
    for (int i_s = 0; i_s < pdf_disable[i_d].size(); i_s++){
      if (pdf_disable[i_d][i_s] == ' '){continue;}
      else if (pdf_disable[i_d][i_s] == '+'){
	if (temp_split != ""){temp_pdf_disable.push_back(temp_split); temp_split = "";} // counter++; 
      }
      else {temp_split.push_back(pdf_disable[i_d][i_s]);}
    }
    if (temp_split != ""){temp_pdf_disable.push_back(temp_split);}
  }

  // according to N_f_active, further particles are added to the list of temp_pdf_disable:
  if (N_f_active < 6){temp_pdf_disable.push_back("top");}
  if (N_f_active < 5){temp_pdf_disable.push_back("bottom");}
  if (N_f_active < 4){temp_pdf_disable.push_back("charm");}
  if (N_f_active < 3){temp_pdf_disable.push_back("strange");}
  if (N_f_active < 2){temp_pdf_disable.push_back("up");}
  if (N_f_active < 1){temp_pdf_disable.push_back("down");}

  // according to chosen PDF set, photons are added to the list of temp_pdf_disable:
  if (!(LHAPDF::hasPhoton())){temp_pdf_disable.push_back("photon");}

  logger << LOG_DEBUG << "temp_pdf_disable.size() = " << temp_pdf_disable.size() << endl;
  for (int i_d = 0; i_d < temp_pdf_disable.size(); i_d++){
    logger << LOG_DEBUG << "temp_pdf_disable[" << i_d << "] = " << temp_pdf_disable[i_d] << endl;
  }


  for (int i_d = 0; i_d < temp_pdf_disable.size(); i_d++){
    if (temp_pdf_disable[i_d] == "Q" ||
	temp_pdf_disable[i_d] == "q" ||
	temp_pdf_disable[i_d] == "qx" ||
	temp_pdf_disable[i_d] == "g" ||
	temp_pdf_disable[i_d] == "a" ||
	temp_pdf_disable[i_d] == "D" ||
	temp_pdf_disable[i_d] == "U" ||
	temp_pdf_disable[i_d] == "S" ||
	temp_pdf_disable[i_d] == "C" ||
	temp_pdf_disable[i_d] == "B" ||
	temp_pdf_disable[i_d] == "T" ||
	temp_pdf_disable[i_d] == "d" ||
	temp_pdf_disable[i_d] == "u" ||
	temp_pdf_disable[i_d] == "s" ||
	temp_pdf_disable[i_d] == "c" ||
	temp_pdf_disable[i_d] == "b" ||
	temp_pdf_disable[i_d] == "t" ||
	temp_pdf_disable[i_d] == "dx" ||
	temp_pdf_disable[i_d] == "ux" ||
	temp_pdf_disable[i_d] == "sx" ||
	temp_pdf_disable[i_d] == "cx" ||
	temp_pdf_disable[i_d] == "bx" ||
	temp_pdf_disable[i_d] == "tx"){
      for (int i_c = 0; i_c < content_parton.size(); i_c++){
	for (int i_p = 0; i_p < content_parton[i_c].size(); i_p++){
	  if ((content_parton[i_c][i_p] == 0 && temp_pdf_disable[i_d] == "g") ||
	      (content_parton[i_c][i_p] == 7 && temp_pdf_disable[i_d] == "a") ||
	      (content_parton[i_c][i_p] == 1 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "q" || temp_pdf_disable[i_d] == "d" || temp_pdf_disable[i_d] == "D")) ||
	      (content_parton[i_c][i_p] == 2 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "q" || temp_pdf_disable[i_d] == "u" || temp_pdf_disable[i_d] == "U")) ||
	      (content_parton[i_c][i_p] == 3 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "q" || temp_pdf_disable[i_d] == "s" || temp_pdf_disable[i_d] == "S")) ||
	      (content_parton[i_c][i_p] == 4 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "q" || temp_pdf_disable[i_d] == "c" || temp_pdf_disable[i_d] == "C")) ||
	      (content_parton[i_c][i_p] == 5 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "q" || temp_pdf_disable[i_d] == "b" || temp_pdf_disable[i_d] == "B")) ||
	      (content_parton[i_c][i_p] == 6 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "q" || temp_pdf_disable[i_d] == "t" || temp_pdf_disable[i_d] == "T")) ||
	      (content_parton[i_c][i_p] == -1 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "qx" || temp_pdf_disable[i_d] == "dx" || temp_pdf_disable[i_d] == "D")) ||
	      (content_parton[i_c][i_p] == -2 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "qx" || temp_pdf_disable[i_d] == "ux" || temp_pdf_disable[i_d] == "U")) ||
	      (content_parton[i_c][i_p] == -3 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "qx" || temp_pdf_disable[i_d] == "sx" || temp_pdf_disable[i_d] == "S")) ||
	      (content_parton[i_c][i_p] == -4 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "qx" || temp_pdf_disable[i_d] == "cx" || temp_pdf_disable[i_d] == "C")) ||
	      (content_parton[i_c][i_p] == -5 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "qx" || temp_pdf_disable[i_d] == "bx" || temp_pdf_disable[i_d] == "B")) ||
	      (content_parton[i_c][i_p] == -6 && (temp_pdf_disable[i_d] == "Q" || temp_pdf_disable[i_d] == "qx" || temp_pdf_disable[i_d] == "tx" || temp_pdf_disable[i_d] == "T"))){
	    content_parton[i_c].erase(content_parton[i_c].begin() + i_p);
	    i_p--;
	  }
	}
      }
    }
    else if (temp_pdf_disable[i_d] == "quark" ||
	temp_pdf_disable[i_d] == "antiquark" ||
	temp_pdf_disable[i_d] == "gluon" ||
	temp_pdf_disable[i_d] == "photon" ||
	temp_pdf_disable[i_d] == "down" ||
	temp_pdf_disable[i_d] == "up" ||
	temp_pdf_disable[i_d] == "strange" ||
	temp_pdf_disable[i_d] == "charm" ||
	temp_pdf_disable[i_d] == "bottom" ||
	temp_pdf_disable[i_d] == "top" ||
	temp_pdf_disable[i_d] == "down-quark" ||
	temp_pdf_disable[i_d] == "up-quark" ||
	temp_pdf_disable[i_d] == "strange-quark" ||
	temp_pdf_disable[i_d] == "charm-quark" ||
	temp_pdf_disable[i_d] == "bottom-quark" ||
	temp_pdf_disable[i_d] == "top-quark" ||
	temp_pdf_disable[i_d] == "down-antiquark" ||
	temp_pdf_disable[i_d] == "up-antiquark" ||
	temp_pdf_disable[i_d] == "strange-antiquark" ||
	temp_pdf_disable[i_d] == "charm-antiquark" ||
	temp_pdf_disable[i_d] == "bottom-antiquark" ||
	temp_pdf_disable[i_d] == "top-antiquark"){
      for (int i_c = 0; i_c < content_parton.size(); i_c++){
	for (int i_p = 0; i_p < content_parton[i_c].size(); i_p++){
	  if ((content_parton[i_c][i_p] == 0 && temp_pdf_disable[i_d] == "gluon") ||
	      (content_parton[i_c][i_p] == 7 && temp_pdf_disable[i_d] == "photon") ||
	      (content_parton[i_c][i_p] == 1 && (temp_pdf_disable[i_d] == "quark" || temp_pdf_disable[i_d] == "down" || temp_pdf_disable[i_d] == "down-quark")) ||
	      (content_parton[i_c][i_p] == 2 && (temp_pdf_disable[i_d] == "quark" || temp_pdf_disable[i_d] == "up" || temp_pdf_disable[i_d] == "up-quark")) ||
	      (content_parton[i_c][i_p] == 3 && (temp_pdf_disable[i_d] == "quark" || temp_pdf_disable[i_d] == "strange" || temp_pdf_disable[i_d] == "strange-quark")) ||
	      (content_parton[i_c][i_p] == 4 && (temp_pdf_disable[i_d] == "quark" || temp_pdf_disable[i_d] == "charm" || temp_pdf_disable[i_d] == "charm-quark")) ||
	      (content_parton[i_c][i_p] == 5 && (temp_pdf_disable[i_d] == "quark" || temp_pdf_disable[i_d] == "bottom" || temp_pdf_disable[i_d] == "bottom-quark")) ||
	      (content_parton[i_c][i_p] == 6 && (temp_pdf_disable[i_d] == "quark" || temp_pdf_disable[i_d] == "top" || temp_pdf_disable[i_d] == "top-quark")) ||
	      (content_parton[i_c][i_p] == -1 && (temp_pdf_disable[i_d] == "antiquark" || temp_pdf_disable[i_d] == "down" || temp_pdf_disable[i_d] == "down-antiquark")) ||
	      (content_parton[i_c][i_p] == -2 && (temp_pdf_disable[i_d] == "antiquark" || temp_pdf_disable[i_d] == "up" || temp_pdf_disable[i_d] == "up-antiquark")) ||
	      (content_parton[i_c][i_p] == -3 && (temp_pdf_disable[i_d] == "antiquark" || temp_pdf_disable[i_d] == "strange" || temp_pdf_disable[i_d] == "strange-antiquark")) ||
	      (content_parton[i_c][i_p] == -4 && (temp_pdf_disable[i_d] == "antiquark" || temp_pdf_disable[i_d] == "charm" || temp_pdf_disable[i_d] == "charm-antiquark")) ||
	      (content_parton[i_c][i_p] == -5 && (temp_pdf_disable[i_d] == "antiquark" || temp_pdf_disable[i_d] == "bottom" || temp_pdf_disable[i_d] == "bottom-antiquark")) ||
	      (content_parton[i_c][i_p] == -6 && (temp_pdf_disable[i_d] == "antiquark" || temp_pdf_disable[i_d] == "top" || temp_pdf_disable[i_d] == "top-antiquark"))){
	    content_parton[i_c].erase(content_parton[i_c].begin() + i_p);
	    i_p--;
	  }
	}
      }
    }
    else {
      logger << LOG_FATAL << "Wrong input: pdf_disable cannot be " << temp_pdf_disable[i_d] << " !" << endl;
      exit(1);
    }
  }

  for (int i_c = 0; i_c < content_parton.size(); i_c++){
    stringstream temp;
    temp.str("");
    temp << setw(16) << "content_parton[" << setw(2) << i_c << "] = ";
    for (int i_p = 0; i_p < content_parton[i_c].size(); i_p++){
      temp << setw(4) << content_parton[i_c][i_p];
    }
    logger << LOG_DEBUG << temp.str() << endl;
  }

  map<string, vector<int> > relation_content_parton;
  relation_content_parton["p"] = content_parton[0];
  relation_content_parton["Q"] = content_parton[1];
  relation_content_parton["q"] = content_parton[2];
  relation_content_parton["qx"] = content_parton[3];
  relation_content_parton["g"] = content_parton[4];
  relation_content_parton["a"] = content_parton[5];
  relation_content_parton["D"] = content_parton[6];
  relation_content_parton["U"] = content_parton[7];
  relation_content_parton["S"] = content_parton[8];
  relation_content_parton["C"] = content_parton[9];
  relation_content_parton["B"] = content_parton[10];
  relation_content_parton["T"] = content_parton[11];
  relation_content_parton["d"] = content_parton[12];
  relation_content_parton["u"] = content_parton[13];
  relation_content_parton["s"] = content_parton[14];
  relation_content_parton["c"] = content_parton[15];
  relation_content_parton["b"] = content_parton[16];
  relation_content_parton["t"] = content_parton[17];
  relation_content_parton["dx"] = content_parton[18];
  relation_content_parton["ux"] = content_parton[19];
  relation_content_parton["sx"] = content_parton[20];
  relation_content_parton["cx"] = content_parton[21];
  relation_content_parton["bx"] = content_parton[22];
  relation_content_parton["tx"] = content_parton[23];


  //  pdf_selection
  /*
  vector<string> temp_pdf_selection(1);
  int start = 0;
  int counter = 0;
  for (int i_s = 0; i_s < pdf_selection.size(); i_s++){
    if (pdf_selection[i_s] == ' '){continue;}
    else if (pdf_selection[i_s] == '+'){counter++; temp_pdf_selection.push_back("");}
    else {temp_pdf_selection[counter].push_back(pdf_selection[i_s]);}
  }
  */
  vector<string> temp_pdf_selection;
  for (int i_d = 0; i_d < pdf_selection.size(); i_d++){
    string temp_split = "";
    for (int i_s = 0; i_s < pdf_selection[i_d].size(); i_s++){
      if (pdf_selection[i_d][i_s] == ' '){continue;}
      else if (pdf_selection[i_d][i_s] == '+'){
	if (temp_split != ""){temp_pdf_selection.push_back(temp_split); temp_split = "";} // counter++; 
      }
      else {temp_split.push_back(pdf_selection[i_d][i_s]);}
    }
    if (temp_split != ""){temp_pdf_selection.push_back(temp_split);}
  }

  if (temp_pdf_selection.size() == 0){temp_pdf_selection.push_back("pp");}

  for (int i_d = 0; i_d < temp_pdf_selection.size(); i_d++){
    logger << LOG_DEBUG << "temp_pdf_selection[" << i_d << "] = " << temp_pdf_selection[i_d] << endl;
  }


  vector<vector<string> > temp_pdf_selection_pair(temp_pdf_selection.size(), vector<string> (2));

  for (int i_p = 0; i_p < temp_pdf_selection.size(); i_p++){
    int start = 1;
    for (int i_s = temp_pdf_selection[i_p].size() - 1; i_s >= 0; i_s--){
      if (i_s > 0 && start == -1){
	logger << LOG_FATAL << "Wrong input: pdf_selection cannot be " << temp_pdf_selection[i_p] << " !" << endl;
	exit(1);
      }
      if (temp_pdf_selection[i_p][i_s] == 'x' && i_s > 0){
	temp_pdf_selection_pair[i_p][start] = temp_pdf_selection[i_p].substr(i_s - 1, 2);
	logger << LOG_DEBUG << "temp_pdf_selection_pair[i_p = " << i_p << "][start = " << start << "] = " << temp_pdf_selection_pair[i_p][start] << endl;
	start--;
	i_s--;
      }
      else {
	temp_pdf_selection_pair[i_p][start] = temp_pdf_selection[i_p].substr(i_s, 1);
	logger << LOG_DEBUG << "temp_pdf_selection_pair[i_p = " << i_p << "][start = " << start << "] = " << temp_pdf_selection_pair[i_p][start] << endl;
	start--;
      }
    }
  }

  for (int i_p = 0; i_p < temp_pdf_selection_pair.size(); i_p++){
    stringstream temp;
    temp.str("");
    temp << setw(26) << "temp_pdf_selection_pair[" << setw(2) << i_p << "] = ";
    for (int i_s = 0; i_s < temp_pdf_selection_pair[i_p].size(); i_s++){
      temp << setw(4) << temp_pdf_selection_pair[i_p][i_s];
    }
    logger << LOG_DEBUG << temp.str() << endl;
  }

  allowed_all_pdf.resize(2);
  vector<int> temp_all_pdf(3);
  int j_x;
  for (int i_x = 0; i_x < 2; i_x++){
    //  for (int i_x = 1; i_x >= -1; i_x -= 2){
    if (i_x == 0){j_x = 1;}
    else if (i_x == 1){j_x = -1;}
    temp_all_pdf[0] = j_x;
    for (int i_p = 0; i_p < temp_pdf_selection_pair.size(); i_p++){
      //      logger << LOG_DEBUG << "temp_pdf_selection_pair[" << i_p << "][0] = " << temp_pdf_selection_pair[i_p][0] << "  ->  " << relation_content_parton[temp_pdf_selection_pair[i_p][0]].size() << endl;
      //      logger << LOG_DEBUG << "temp_pdf_selection_pair[" << i_p << "][1] = " << temp_pdf_selection_pair[i_p][1] << "  ->  " << relation_content_parton[temp_pdf_selection_pair[i_p][1]].size() << endl;
      for (int i_i = 0; i_i < relation_content_parton[temp_pdf_selection_pair[i_p][0]].size(); i_i++){
	for (int j_i = 0; j_i < relation_content_parton[temp_pdf_selection_pair[i_p][1]].size(); j_i++){
	  if (j_x == 1){
	    temp_all_pdf[1] = relation_content_parton[temp_pdf_selection_pair[i_p][0]][i_i];
	    temp_all_pdf[2] = relation_content_parton[temp_pdf_selection_pair[i_p][1]][j_i];
	  }
	  else if (j_x == -1){
	    temp_all_pdf[1] = relation_content_parton[temp_pdf_selection_pair[i_p][1]][j_i];
	    temp_all_pdf[2] = relation_content_parton[temp_pdf_selection_pair[i_p][0]][i_i];
	  }
	  allowed_all_pdf[i_x].push_back(temp_all_pdf);
	}
      }
    }
  }

  for (int i_x = 0; i_x < 2; i_x++){
    for (int i_p = 0; i_p < allowed_all_pdf[i_x].size(); i_p++){
      stringstream temp;
      temp.str("");
      temp << setw(26) << "allowed_all_pdf[" << setw(2) << i_x << "][" << setw(2) << i_p << "] = ";
      for (int i_s = 0; i_s < allowed_all_pdf[i_x][i_p].size(); i_s++){
	temp << setw(4) << allowed_all_pdf[i_x][i_p][i_s];
      }
      logger << LOG_DEBUG << temp.str() << endl;
    }
  }
  //  pdf_selection
  //  pdf_disable

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}




void observable_set::perform_selection_content_pdf(){
  static Logger logger("observable_set::perform_selection_content_pdf");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (!coll_choice){return;}

  //  determine_selection_pdf_content(oset);
  determine_selection_content_pdf();
  logger << LOG_DEBUG << "Before perform_selection_content_pdf:" << endl;
  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){
    stringstream temp;
    temp.str("");
    temp << "combination_pdf[" << setw(2) << i_x << "] = "; 
    for (int i_y = 0; i_y < combination_pdf[i_x].size(); i_y++){temp << setw(2) << right << combination_pdf[i_x][i_y] << "  ";}
    logger << LOG_DEBUG << temp.str() << endl;
  }

  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){
    int match = 0;
    int x_d = 0;
    if (combination_pdf[i_x][0] == 1){x_d = 0;}
    else if (combination_pdf[i_x][0] == -1){x_d = 1;}
    for (int i_p = 0; i_p < allowed_all_pdf[x_d].size(); i_p++){
      if (combination_pdf[i_x] == allowed_all_pdf[x_d][i_p]){
	match = 1; 
	break;
      }
    }
    if (match == 0){
      combination_pdf.erase(combination_pdf.begin() + i_x); 
      i_x--;
    }
  }

  logger << LOG_DEBUG << "After perform_selection_content_pdf:" << endl;
  for (int i_x = 0; i_x < combination_pdf.size(); i_x++){
    stringstream temp;
    temp.str("");
    temp << "combination_pdf[" << setw(2) << i_x << "] = "; 
    for (int i_y = 0; i_y < combination_pdf[i_x].size(); i_y++){temp << setw(2) << right << combination_pdf[i_x][i_y] << "  ";}
    logger << LOG_DEBUG << temp.str() << endl;
  }
  if (combination_pdf.size() == 0){
    logger << LOG_FATAL << "According to pdf-selection criteria, no partonic channels remain." << endl;
    int_end = 1;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::perform_selection_content_pdf_collinear(){
  static Logger logger("observable_set::perform_selection_content_pdf_collinear");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (!coll_choice){int_end = 1; return;}

  //  determine_selection_pdf_content(oset);
  determine_selection_content_pdf();
  logger << LOG_INFO << "Before perform_selection_content_pdf_collinear:" << endl;
  logger << LOG_DEBUG_VERBOSE << "(*CA_collinear).size() = " << (*CA_collinear).size() << endl;
  logger << LOG_DEBUG_VERBOSE << "CA_combination_pdf.size() = " << CA_combination_pdf.size() << endl;
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      logger << LOG_DEBUG << "pdf contributions: " << setw(15) << left << (*CA_collinear)[i_c][0].all_name()[i_i] << ":   " << endl;
      for (int i_x = 0; i_x < CA_combination_pdf[i_c][i_i].size(); i_x++){
	stringstream temp;
	temp.str("");
	temp << setw(48) << "";
	temp << "CA_combination_pdf[" << setw(2) << i_c << "][" << setw(2) << i_i << "][" << setw(2) << i_x << "] = "; 
	for (int i_y = 0; i_y < CA_combination_pdf[i_c][i_i][i_x].size(); i_y++){temp << setw(2) << right << CA_combination_pdf[i_c][i_i][i_x][i_y] << "  ";}
	logger << LOG_DEBUG << temp.str() << endl;
      }
    }
  }
  logger.newLine(LOG_DEBUG);

  for (int i_c = 0; i_c < CA_combination_pdf.size(); i_c++){
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      for (int i_x = 0; i_x < CA_combination_pdf[i_c][i_i].size(); i_x++){
	int match = 0;
	int x_d = 0;
	if (CA_combination_pdf[i_c][i_i][i_x][0] == 1){x_d = 0;}
	else if (CA_combination_pdf[i_c][i_i][i_x][0] == -1){x_d = 1;}
	for (int i_p = 0; i_p < allowed_all_pdf[x_d].size(); i_p++){
	  if (CA_combination_pdf[i_c][i_i][i_x] == allowed_all_pdf[x_d][i_p]){
	    match = 1; 
	    break;
	  }
	}
	if (match == 0){
	  CA_combination_pdf[i_c][i_i].erase(CA_combination_pdf[i_c][i_i].begin() + i_x); 
	  CA_Q2f[i_c][i_i].erase(CA_Q2f[i_c][i_i].begin() + i_x); 
	  i_x--;
	}
      }
    }
  }

  logger << LOG_INFO << "After selection, before removing empty entries, in perform_selection_content_pdf_collinear:" << endl;
  logger.newLine(LOG_DEBUG);
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      logger << LOG_DEBUG << "pdf contributions: " << setw(15) << left << (*CA_collinear)[i_c][0].all_name()[i_i] << ":   " << endl;
      for (int i_x = 0; i_x < CA_combination_pdf[i_c][i_i].size(); i_x++){
	stringstream temp;
	temp.str("");
	temp << setw(48) << "";
	temp << "CA_combination_pdf[" << setw(2) << i_c << "][" << setw(2) << i_i << "][" << setw(2) << i_x << "] = "; 
	for (int i_y = 0; i_y < CA_combination_pdf[i_c][i_i][i_x].size(); i_y++){temp << setw(2) << right << CA_combination_pdf[i_c][i_i][i_x][i_y] << "  ";}
	logger << LOG_DEBUG << temp.str() << endl;
      }
    }
  }
  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG << "Remove now those collinear 'dipoles' that do not have any left-over contribution:" << endl;
  logger.newLine(LOG_DEBUG);

  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    int count_remaining = 0;
    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
      count_remaining += CA_combination_pdf[i_c][i_i].size();
    }
    if (count_remaining == 0){
      CA_combination_pdf.erase(CA_combination_pdf.begin() + i_c);
      (*CA_collinear).erase((*CA_collinear).begin() + i_c);
      CA_Q2f.erase(CA_Q2f.begin() + i_c);
      i_c--;
    }
  }
  logger.newLine(LOG_DEBUG);

  n_pc = (*CA_collinear).size();


  logger << LOG_INFO << "After perform_selection_content_pdf_collinear:" << endl;

  logger.newLine(LOG_INFO);
  output_collinear();
  logger.newLine(LOG_INFO);
  output_collinear_pdf();

  if ((*CA_collinear).size() == 0){
    logger << LOG_FATAL << "According to pdf-selection criteria, no partonic channels remain." << endl;
    int_end = 1;
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::perform_selection_content_pdf_list(){
  static Logger logger("observable_set::perform_selection_content_pdf_list");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  determine_selection_content_pdf();

  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "Before perform_selection_content_pdf_list:   list_combination_pdf.size() = " << list_combination_pdf.size() << endl;
  logger.newLine(LOG_DEBUG);
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    logger << LOG_DEBUG << "list_combination_pdf[" << setw(2) << i_l << "]:" << endl; 
    for (int i_i = 0; i_i < list_combination_pdf[i_l].size(); i_i++){
      for (int i_z = 0; i_z < list_combination_pdf[i_l][i_i].size(); i_z++){
	stringstream temp;
	temp.str("");
	temp << setw(32) << "";
	temp << "list_combination_pdf[" << setw(2) << i_l << "][" << setw(2) << i_i << "][" << setw(2) << i_z << "] = "; 
	for (int i_y = 0; i_y < list_combination_pdf[i_l][i_i][i_z].size(); i_y++){temp << setw(2) << right << list_combination_pdf[i_l][i_i][i_z][i_y] << "  ";}
	logger << LOG_DEBUG << temp.str() << endl;
      }
    }
  }
  logger.newLine(LOG_DEBUG);

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    for (int i_i = 0; i_i < list_combination_pdf[i_l].size(); i_i++){
      for (int i_x = 0; i_x < list_combination_pdf[i_l][i_i].size(); i_x++){
	int match = 0;
	int x_d = 0;
	if (list_combination_pdf[i_l][i_i][i_x][0] == 1){x_d = 0;}
	else if (list_combination_pdf[i_l][i_i][i_x][0] == -1){x_d = 1;}
	for (int i_p = 0; i_p < allowed_all_pdf[x_d].size(); i_p++){
	  if (list_combination_pdf[i_l][i_i][i_x] == allowed_all_pdf[x_d][i_p]){
	    match = 1; 
	    break;
	  }
	}
	if (match == 0){
	  list_combination_pdf[i_l][i_i].erase(list_combination_pdf[i_l][i_i].begin() + i_x); 
	  i_x--;
	}
      }
    }
  }

  logger << LOG_DEBUG << "Remove now those entries of 'list_combination_pdf[i_l]' that do not have any left-over contribution." << endl;
  logger.newLine(LOG_DEBUG);
  
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    for (int i_i = list_combination_pdf[i_l].size() - 1; i_i >= 0; i_i--){
      if (list_combination_pdf[i_l][i_i].size() == 0){
	list_combination_pdf[i_l].erase(list_combination_pdf[i_l].begin() + i_i);
      }
    }
  }
  
  logger << LOG_DEBUG << "After perform_selection_content_pdf_list:   list_combination_pdf.size() = " << list_combination_pdf.size() << endl;
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    logger << LOG_DEBUG << "list_combination_pdf[" << setw(2) << i_l << "]:" << endl; 
    for (int i_i = 0; i_i < list_combination_pdf[i_l].size(); i_i++){
      for (int i_z = 0; i_z < list_combination_pdf[i_l][i_i].size(); i_z++){
	stringstream temp;
	temp.str("");
	temp << setw(32) << "";
	temp << "list_combination_pdf[" << setw(2) << i_l << "][" << setw(2) << i_i << "][" << setw(2) << i_z << "] = "; 
	for (int i_y = 0; i_y < list_combination_pdf[i_l][i_i][i_z].size(); i_y++){temp << setw(2) << right << list_combination_pdf[i_l][i_i][i_z][i_y] << "  ";}
	logger << LOG_DEBUG << temp.str() << endl;
      }
    }
  }
  logger.newLine(LOG_DEBUG);


  logger.newLine(LOG_DEBUG);
  logger << LOG_DEBUG << "Remove now those entries of 'ncollinear' that do not have any left-over contribution:" << endl;
  logger.newLine(LOG_DEBUG);
  int counter_ncollinear = 0;
  for (int i_c = ncollinear.size() - 1; i_c >= 0; i_c--){
    logger << LOG_DEBUG_VERBOSE << "ncollinear[i_c].no_pdf = " << ncollinear[i_c].no_pdf << endl;
    logger << LOG_DEBUG_VERBOSE << "list_combination_pdf[ncollinear[i_c].no_pdf].size() = " << list_combination_pdf[ncollinear[i_c].no_pdf].size() << endl;
    for (int i_i = 0; i_i < list_combination_pdf[ncollinear[i_c].no_pdf].size(); i_i++){
      logger << LOG_DEBUG_VERBOSE << "list_combination_pdf[ncollinear[i_c].no_pdf][i_i].size() = " << list_combination_pdf[ncollinear[i_c].no_pdf][i_i].size() << endl;
    }

    if (list_combination_pdf[ncollinear[i_c].no_pdf].size() == 0){
      logger << LOG_DEBUG << ncollinear[i_c].name << " [" << setw(2) << i_c << "] is removed due to parton content." << endl;

      //      ncollinear.erase(ncollinear.begin() + i_c); // !!! check if this works -> doesn't work yet if i_c == 0 is directly accessed !!!

      counter_ncollinear++;
    }
  }

  if (ncollinear.size() == 0){
    logger << LOG_FATAL << "According to pdf-selection criteria, no partonic channels remain." << endl;
    int_end = 1;
  }
  // 
  if (counter_ncollinear == ncollinear.size()){
    logger << LOG_FATAL << "According to pdf-selection criteria, no partonic channels remain." << endl;
    int_end = 1;
  }


  
  // !!! Has nothing to do with pdfs: Only relevant on Born phase-space (e.g. CT[2], VT[2])
  
  initial_gg = 0;
  initial_qqx = 0;
  if (csi->type_parton[0][1] == 0 && 
      csi->type_parton[0][2] == 0){initial_gg = 1;}  // hard process: IS: gg
  else if (csi->type_parton[0][1] != 0 && 
	   csi->type_parton[0][1] != 22 && 
	   csi->type_parton[0][1] == -csi->type_parton[0][2]){initial_qqx = 1;}  // hard process: IS: qqx (SF, neutral case)
  else if (csi->type_parton[0][1] != 0 && 
	   csi->type_parton[0][1] != 22 && 
	   csi->type_parton[0][1] * csi->type_parton[0][2] < 0 &&
	   (abs(csi->type_parton[0][1]) + abs(csi->type_parton[0][2])) % 2 == 1){initial_qqx = 2;}  // hard process: IS: qqx (all flavour combinations, charged case) - might need refinement !!!)



  initial_pdf_gg = 0;
  initial_pdf_qqx = 0;
  initial_pdf_gq = 0;
  initial_pdf_rest = 0;
  initial_pdf_diag = 0;

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    for (int i_i = 0; i_i < list_combination_pdf[i_l].size(); i_i++){
      for (int i_z = 0; i_z < list_combination_pdf[i_l][i_i].size(); i_z++){
	if (list_combination_pdf[i_l][i_i][i_z][1] == 0 && 
	    list_combination_pdf[i_l][i_i][i_z][2] == 0){initial_pdf_gg = 1;}  // PDF initial state: gg
	else if ((list_combination_pdf[i_l][i_i][i_z][1] == 0 && 
		  list_combination_pdf[i_l][i_i][i_z][2] != 0 && 
		  list_combination_pdf[i_l][i_i][i_z][2] != 7) ||
		 (list_combination_pdf[i_l][i_i][i_z][2] == 0 && 
		  list_combination_pdf[i_l][i_i][i_z][1] != 0 && 
		  list_combination_pdf[i_l][i_i][i_z][1] != 7)){initial_pdf_gq = 1;}  // PDF initial state: gQ or Qg
	else if (list_combination_pdf[i_l][i_i][i_z][1] != 0 && 
		 list_combination_pdf[i_l][i_i][i_z][1] != 7 && 
		 list_combination_pdf[i_l][i_i][i_z][1] == -list_combination_pdf[i_l][i_i][i_z][2]){initial_pdf_qqx = 1;}  // PDF initial state: qqx (SF only: might need refinement for charged case !!!)
	else {initial_pdf_rest = 1;}  // PDF initial state: all other combinations, i.e. qq, qxqx, qqx (DF only: might need refinement for charged case !!!)

	if (list_combination_pdf[i_l][i_i][i_z][1] == csi->type_parton[0][1] && 
	    list_combination_pdf[i_l][i_i][i_z][2] == csi->type_parton[0][2]){initial_pdf_diag = 1;}  // PDF IS equal to hard-process IS: relevant for e.g. if H2 is needed, etc.
      }
    }
  }

  initial_diag = 0;
  initial_diag_gg = 0;
  initial_diag_qqx = 0;
  
  if (initial_gg && initial_pdf_gg){
    initial_diag_gg = 1;
    initial_diag = 1;
  }
  if (initial_qqx && initial_pdf_qqx){
    initial_diag_qqx = 1;
    initial_diag = 2;
  }

  initial_channel = 0;
  if (initial_gg){initial_channel = 1;}
  else if (initial_qqx){initial_channel = 2;}
  //  else {initial_channel = 0;}



  /*
  if (initial_gg && initial_pdf_gg){initial_diag = 1;}
  if (initial_qqx && initial_pdf_qqx){initial_diag = 2;}
  //  if (initial_channel){initial_diag = 1;}
  */

  logger << LOG_INFO << "initial_channel  = " << initial_channel << endl;

  logger << LOG_INFO << "initial_gg       = " << initial_gg << endl;
  logger << LOG_INFO << "initial_qqx      = " << initial_qqx << endl;

  logger << LOG_INFO << "initial_pdf_gg   = " << initial_pdf_gg << endl;
  logger << LOG_INFO << "initial_pdf_gq   = " << initial_pdf_gq << endl;
  logger << LOG_INFO << "initial_pdf_qqx  = " << initial_pdf_qqx << endl;
  logger << LOG_INFO << "initial_pdf_rest = " << initial_pdf_rest << endl;

  logger << LOG_INFO << "initial_pdf_diag = " << initial_pdf_diag << endl;

  logger << LOG_INFO << "initial_diag     = " << initial_diag << endl;
  logger << LOG_INFO << "initial_diag_gg  = " << initial_diag_gg << endl;
  logger << LOG_INFO << "initial_diag_qqx = " << initial_diag_qqx << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::calculate_pdf_LHAPDF_list_CV(){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_list_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    if (value_mu_fact[sd].size() == 0){continue;}
    for (int ss = 0; ss < value_mu_fact[sd].size(); ss++){
      calculate_pdf_LHAPDF_scale_list_CV(sd, ss);
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_pdf_LHAPDF_scale_list_CV(int sd, int ss){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_scale_list_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static vector<vector<vector<double> > > lhapdf_result;
  static int hasPhoton = LHAPDF::hasPhoton();
  static vector<int> max_emission(3, 0);
  if (initialization == 1){
    if (hasPhoton == 0){lhapdf_result.resize(3, vector<vector<double> > (3, vector<double> (13)));}
    else {lhapdf_result.resize(3, vector<vector<double> > (3, vector<double> (14)));}
    for (int i_x = 1; i_x < 3; i_x++){
      for (int i_c = 0; i_c < ncollinear.size(); i_c++){
	if (ncollinear[i_c].n_emission_leg_type[i_x][0] > 0){max_emission[i_x] = 1;}
      }
    }
    initialization = 0;
  }

  // new meaning for multicollinear !!!
  // xz_pdf: 0 - 0 : empty - could be filled with tau = x1 * x2
  //         1 - 0 : x1 -> no collinear emission from 1st parton
  //         2 - 0 : x2 -> no collinear emission from 2nd parton
  //         0 - 1 : empty
  //         1 - 1 : x1 / z1 -> collinear emission from 1st parton
  //         2 - 1 : x2 / z2 -> collinear emission from 2nd parton

  for (int i_x = 1; i_x < 3; i_x++){
    for (int i_z = 0; i_z <= max_emission[i_x]; i_z++){
      if (hasPhoton == 0){lhapdf_result[i_x][i_z] = LHAPDF::xfx(xz_pdf[i_x][i_z], value_mu_fact[sd][ss]);}
      else {lhapdf_result[i_x][i_z] = LHAPDF::xfxphoton(xz_pdf[i_x][i_z], value_mu_fact[sd][ss]);}
      //      modify_pdf_content(lhapdf_result[i_x][i_z], *this);
      for (int i_p = 0; i_p < lhapdf_result[i_x][i_z].size(); i_p++){logger << LOG_DEBUG_VERBOSE << "xz_pdf[" << i_x << "][" << i_z << "] = " << setw(23) << setprecision(15) << xz_pdf[i_x][i_z] << "   lhapdf_result[" << i_x << "][" << i_z << "][" << i_p << "] = " << lhapdf_result[i_x][i_z][i_p] << endl;}
    }
  }
  
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    value_list_pdf_factor[sd][ss][i_l] = vector<double> (3, 0.);
    if (list_combination_pdf[i_l].size() == 0){continue;}
    int i_z1 = list_combination_pdf_emission[i_l][1];
    int i_z2 = list_combination_pdf_emission[i_l][2];
    //    cout << "ncollinear[" << i_c << "].emission[1] = " << ncollinear[i_c].emission[1] << endl;
    //    cout << "ncollinear[" << i_c << "].emission[2] = " << ncollinear[i_c].emission[2] << endl;
    for (int i_i = 0; i_i < list_combination_pdf[i_l].size(); i_i++){
      for (int i_x = 0; i_x < list_combination_pdf[i_l][i_i].size(); i_x++){
	double factor1 = 0.;
	double factor2 = 0.;
	int direction = 1;
	if (list_combination_pdf[i_l][i_i][i_x][0] == -1){direction = 2;}
	
	if (list_combination_pdf[i_l][i_i][i_x][1] == 7){factor1 = lhapdf_result[1][i_z1][13];}
	else if (coll_choice == 1){factor1 = lhapdf_result[1][i_z1][6 + list_combination_pdf[i_l][i_i][i_x][1]];}
	else if (coll_choice == 2){
	  if      (list_combination_pdf[i_l][i_i][i_x][0] ==  1){factor1 = lhapdf_result[1][i_z1][6 + list_combination_pdf[i_l][i_i][i_x][1]];}
	  else if (list_combination_pdf[i_l][i_i][i_x][0] == -1){factor1 = lhapdf_result[1][i_z1][6 - list_combination_pdf[i_l][i_i][i_x][1]];}
	}
	if (list_combination_pdf[i_l][i_i][i_x][2] == 7){factor2 = lhapdf_result[2][i_z2][13];}
	else if (coll_choice == 1){factor2 = lhapdf_result[2][i_z2][6 + list_combination_pdf[i_l][i_i][i_x][2]];}
	else if (coll_choice == 2){
	  if      (list_combination_pdf[i_l][i_i][i_x][0] ==  1){factor2 = lhapdf_result[2][i_z2][6 - list_combination_pdf[i_l][i_i][i_x][2]];}
	  else if (list_combination_pdf[i_l][i_i][i_x][0] == -1){factor2 = lhapdf_result[2][i_z2][6 + list_combination_pdf[i_l][i_i][i_x][2]];}
	}
	value_list_pdf_factor[sd][ss][i_l][direction] += factor1 * factor2;
      }
    }
    
    for (int i_y = 1; i_y < 3; i_y++){value_list_pdf_factor[sd][ss][i_l][i_y] = value_list_pdf_factor[sd][ss][i_l][i_y] / (xz_pdf[1][i_z1] * xz_pdf[2][i_z2]);}
    value_list_pdf_factor[sd][ss][i_l][0] = value_list_pdf_factor[sd][ss][i_l][1] + value_list_pdf_factor[sd][ss][i_l][2];
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::calculate_pdf_LHAPDF_list_TSV(){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_list_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (!switch_TSV){return;}

  /*
  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    if (value_mu_fact[sd].size() == 0){continue;}
    for (int ss = 0; ss < value_mu_fact[sd].size(); ss++){
      calculate_pdf_LHAPDF_scale_list_TSV(sd, ss);
    }
  }
  */
  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i_v = 0; i_v < max_dyn_fact + 1; i_v++){
      for (int i_m = 0; i_m < n_scale_dyn_fact[i_v]; i_m++){
	calculate_pdf_LHAPDF_scale_list_TSV(i_a, i_v, i_m);
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_pdf_LHAPDF_scale_list_TSV(int i_a, int i_v, int i_m){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_scale_list_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static vector<vector<vector<double> > > lhapdf_result;
  static int hasPhoton = LHAPDF::hasPhoton();
  static vector<int> max_emission(3, 0);
  if (initialization == 1){
    if (hasPhoton == 0){lhapdf_result.resize(3, vector<vector<double> > (3, vector<double> (13)));}
    else {lhapdf_result.resize(3, vector<vector<double> > (3, vector<double> (14)));}
    for (int i_x = 1; i_x < 3; i_x++){
      for (int i_c = 0; i_c < ncollinear.size(); i_c++){
	if (ncollinear[i_c].n_emission_leg_type[i_x][0] > 0){max_emission[i_x] = 1;}
      }
    }
    initialization = 0;
  }

  // new meaning for multicollinear !!!
  // xz_pdf: 0 - 0 : empty - could be filled with tau = x1 * x2
  //         1 - 0 : x1 -> no collinear emission from 1st parton
  //         2 - 0 : x2 -> no collinear emission from 2nd parton
  //         0 - 1 : empty
  //         1 - 1 : x1 / z1 -> collinear emission from 1st parton
  //         2 - 1 : x2 / z2 -> collinear emission from 2nd parton

  for (int i_x = 1; i_x < 3; i_x++){
    for (int i_z = 0; i_z <= max_emission[i_x]; i_z++){
      if (hasPhoton == 0){lhapdf_result[i_x][i_z] = LHAPDF::xfx(xz_pdf[i_x][i_z], value_scale_fact[i_a][i_v][i_m]);}
      else {lhapdf_result[i_x][i_z] = LHAPDF::xfxphoton(xz_pdf[i_x][i_z], value_scale_fact[i_a][i_v][i_m]);}
      //      modify_pdf_content(lhapdf_result[i_x][i_z], *this);
      for (int i_p = 0; i_p < lhapdf_result[i_x][i_z].size(); i_p++){logger << LOG_DEBUG_VERBOSE << "xz_pdf[" << i_x << "][" << i_z << "] = " << setw(23) << setprecision(15) << xz_pdf[i_x][i_z] << "   lhapdf_result[" << i_x << "][" << i_z << "][" << i_p << "] = " << lhapdf_result[i_x][i_z][i_p] << endl;}
    }
  }
  
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    value_list_pdf_factor_TSV[i_a][i_v][i_m][i_l] = vector<double> (3, 0.);
    if (list_combination_pdf[i_l].size() == 0){continue;}
    int i_z1 = list_combination_pdf_emission[i_l][1];
    int i_z2 = list_combination_pdf_emission[i_l][2];
    //    cout << "ncollinear[" << i_c << "].emission[1] = " << ncollinear[i_c].emission[1] << endl;
    //    cout << "ncollinear[" << i_c << "].emission[2] = " << ncollinear[i_c].emission[2] << endl;
    for (int i_i = 0; i_i < list_combination_pdf[i_l].size(); i_i++){
      for (int i_x = 0; i_x < list_combination_pdf[i_l][i_i].size(); i_x++){
	double factor1 = 0.;
	double factor2 = 0.;
	int direction = 1;
	if (list_combination_pdf[i_l][i_i][i_x][0] == -1){direction = 2;}
	
	if (list_combination_pdf[i_l][i_i][i_x][1] == 7){factor1 = lhapdf_result[1][i_z1][13];}
	else if (coll_choice == 1){factor1 = lhapdf_result[1][i_z1][6 + list_combination_pdf[i_l][i_i][i_x][1]];}
	else if (coll_choice == 2){
	  if      (list_combination_pdf[i_l][i_i][i_x][0] ==  1){factor1 = lhapdf_result[1][i_z1][6 + list_combination_pdf[i_l][i_i][i_x][1]];}
	  else if (list_combination_pdf[i_l][i_i][i_x][0] == -1){factor1 = lhapdf_result[1][i_z1][6 - list_combination_pdf[i_l][i_i][i_x][1]];}
	}
	if (list_combination_pdf[i_l][i_i][i_x][2] == 7){factor2 = lhapdf_result[2][i_z2][13];}
	else if (coll_choice == 1){factor2 = lhapdf_result[2][i_z2][6 + list_combination_pdf[i_l][i_i][i_x][2]];}
	else if (coll_choice == 2){
	  if      (list_combination_pdf[i_l][i_i][i_x][0] ==  1){factor2 = lhapdf_result[2][i_z2][6 - list_combination_pdf[i_l][i_i][i_x][2]];}
	  else if (list_combination_pdf[i_l][i_i][i_x][0] == -1){factor2 = lhapdf_result[2][i_z2][6 + list_combination_pdf[i_l][i_i][i_x][2]];}
	}
	value_list_pdf_factor_TSV[i_a][i_v][i_m][i_l][direction] += factor1 * factor2;
      }
    }
    
    for (int i_y = 1; i_y < 3; i_y++){
      value_list_pdf_factor_TSV[i_a][i_v][i_m][i_l][i_y] = value_list_pdf_factor_TSV[i_a][i_v][i_m][i_l][i_y] / (xz_pdf[1][i_z1] * xz_pdf[2][i_z2]);
    }
    value_list_pdf_factor_TSV[i_a][i_v][i_m][i_l][0] = value_list_pdf_factor_TSV[i_a][i_v][i_m][i_l][1] + value_list_pdf_factor_TSV[i_a][i_v][i_m][i_l][2];
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


/*
void observable_set::xxx_calculate_pdf_LHAPDF_CX_ncollinear_TSV(vector<vector<double> > & all_xz_coll_pdf){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_CA_collinear_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (osi_coll_choice == 0){cout << "parton--parton scattering CA contribution doesn't make sense!" << endl; exit(1);}
  static int initialization = 1;
  static vector<vector<vector<vector<vector<vector<double> > > > > > lhapdf_result(osi_n_ps, vector<vector<vector<vector<vector<double> > > > > (osi_max_dyn_fact + 1));
  static vector<vector<int> > in_collinear(osi_CA_collinear.size());
  static vector<vector<vector<int> > > all_pdf(osi_CA_collinear.size());
  static int hasPhoton = LHAPDF::hasPhoton();
  if (initialization == 1){
    for (int i_a = 0; i_a < osi_n_ps; i_a++){
      for (int i_v = 0; i_v < osi_max_dyn_fact + 1; i_v++){
	logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << "   i_v = " << i_v << endl;
	if (hasPhoton == 0){lhapdf_result[i_a][i_v].resize(osi_n_scale_dyn_fact[i_v], vector<vector<vector<double> > > (3, vector<vector<double> > (3, vector<double> (13))));}
	else {lhapdf_result[i_a][i_v].resize(osi_n_scale_dyn_fact[i_v], vector<vector<vector<double> > > (3, vector<vector<double> > (3, vector<double> (14))));}
      }
    }
    for (int i_a = 0; i_a < osi_CA_collinear.size(); i_a++){
      logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << endl;
      in_collinear[i_a] = osi_CA_collinear[i_a][0].in_collinear();
      all_pdf[i_a] = osi_CA_collinear[i_a][0].all_pdf();
    }
    initialization = 0;
    logger << LOG_DEBUG_VERBOSE << "initialization finished!" << endl;
  }

  for (int i_a = 0; i_a < osi_n_ps; i_a++){
    for (int i_v = 0; i_v < osi_max_dyn_fact + 1; i_v++){
      logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << "   i_v = " << i_v << endl;
      for (int i_m = 0; i_m < osi_n_scale_dyn_fact[i_v]; i_m++){
	for (int i_x = 1; i_x < 3; i_x++){
	  for (int i_z = 0; i_z < i_x + 1; i_z += i_x){
	    if (hasPhoton == 0){lhapdf_result[i_a][i_v][i_m][i_x][i_z] = LHAPDF::xfx(all_xz_coll_pdf[i_x][i_z], osi_value_scale_fact[i_a][i_v][i_m]);}
	    else {lhapdf_result[i_a][i_v][i_m][i_x][i_z] = LHAPDF::xfxphoton(all_xz_coll_pdf[i_x][i_z], osi_value_scale_fact[i_a][i_v][i_m]);}
 	  }
	}
      }
    }
  }

  for (int i_c = 0; i_c < osi_n_pc; i_c++){
    if (osi_cut_ps[osi_relation_pc_ps[i_c]] == -1){continue;}
    for (int i_z = 0; i_z < osi_n_pz; i_z++){
      if (in_collinear[i_c][i_z] == 0){continue;}
      for (int i_v = 0; i_v < osi_max_dyn_fact + 1; i_v++){
	// !!! check that this is calculated only if osi_cut_ps[...] == 0 !!!
	if (i_v == 0 && osi_relation_pc_ps[i_c] > 0){osi_value_pdf_factor[i_c][i_z][i_v] = osi_value_pdf_factor[0][i_z][i_v];}
	else {
	  for (int i_m = 0; i_m < osi_n_scale_dyn_fact[i_v]; i_m++){
	    for (int i_h = 0; i_h < 3; i_h++){osi_value_pdf_factor[i_c][i_z][i_v][i_m][i_h] = 0.;}
	    //	    for (int i_i = 0; i_i < all_pdf[i_c].size(); i_i++){
	    for (int i_i = 0; i_i < osi_CA_combination_pdf[i_c].size(); i_i++){
	      int i_z1;
	      int i_z2;
	      if      (i_z == 0){i_z1 = 0; i_z2 = 0;}
	      else if (i_z == 1){i_z1 = 1; i_z2 = 0;}
	      else if (i_z == 2){i_z1 = 0; i_z2 = 2;}

	      double factor1 = 0.;
	      double factor2 = 0.;
	      //	      if (all_pdf[i_c][i_i][0] == -1){direction = 2;}
	      for (int i_x = 0; i_x < osi_CA_combination_pdf[i_c][i_i].size(); i_x++){
		int direction = 1;
		if (osi_CA_combination_pdf[i_c][i_i][i_x][0] == -1){direction = 2;}

		// ??? osi_relation_pc_ps[i_c]
		if (all_pdf[i_c][i_i][1] == 7){factor1 = lhapdf_result[osi_relation_pc_ps[i_c]][i_v][i_m][1][i_z1][13];}
		else if (osi_coll_choice == 1){factor1 = lhapdf_result[osi_relation_pc_ps[i_c]][i_v][i_m][1][i_z1][6 + osi_CA_combination_pdf[i_c][i_i][i_x][1]];}
		else if (osi_coll_choice == 2){
		  if      (osi_CA_combination_pdf[i_c][i_i][i_x][0] ==  1){factor1 = lhapdf_result[osi_relation_pc_ps[i_c]][i_v][i_m][1][i_z1][6 + osi_CA_combination_pdf[i_c][i_i][i_x][1]];}
		  else if (osi_CA_combination_pdf[i_c][i_i][i_x][0] == -1){factor1 = lhapdf_result[osi_relation_pc_ps[i_c]][i_v][i_m][1][i_z1][6 - osi_CA_combination_pdf[i_c][i_i][i_x][1]];}
		}
		if (osi_CA_combination_pdf[i_c][i_i][i_x][2] == 7){factor2 = lhapdf_result[osi_relation_pc_ps[i_c]][i_v][i_m][2][i_z2][13];}
		else if (osi_coll_choice == 1){factor2 = lhapdf_result[osi_relation_pc_ps[i_c]][i_v][i_m][2][i_z2][6 + osi_CA_combination_pdf[i_c][i_i][i_x][2]];}
		else if (osi_coll_choice == 2){
		  if      (osi_CA_combination_pdf[i_c][i_i][i_x][0] ==  1){factor2 = lhapdf_result[osi_relation_pc_ps[i_c]][i_v][i_m][2][i_z2][6 - osi_CA_combination_pdf[i_c][i_i][i_x][2]];}
		  else if (osi_CA_combination_pdf[i_c][i_i][i_x][0] == -1){factor2 = lhapdf_result[osi_relation_pc_ps[i_c]][i_v][i_m][2][i_z2][6 + osi_CA_combination_pdf[i_c][i_i][i_x][2]];}
		}

		osi_value_pdf_factor[i_c][i_z][i_v][i_m][direction] += factor1 * factor2;
	      }
	    }
	    for (int i_h = 1; i_h < 3; i_h++){
	      if (i_z == 0){osi_value_pdf_factor[i_c][i_z][i_v][i_m][i_h] = osi_value_pdf_factor[i_c][i_z][i_v][i_m][i_h] / all_xz_coll_pdf[0][0];}
	      else {osi_value_pdf_factor[i_c][i_z][i_v][i_m][i_h] = osi_value_pdf_factor[i_c][i_z][i_v][i_m][i_h] / all_xz_coll_pdf[i_z % 2 + 1][i_z];}
	    }
	    osi_value_pdf_factor[i_c][i_z][i_v][i_m][0] = osi_value_pdf_factor[i_c][i_z][i_v][i_m][1] + osi_value_pdf_factor[i_c][i_z][i_v][i_m][2];
	    //	    cout << "osi_value_pdf_factor[" << i_c << "][" << i_z << "][" << i_v << "][" << i_m << "][0]  = " << osi_value_pdf_factor[i_c][i_z][i_v][i_m][0]  << endl;
	  }
	}
      }
    }
  }

  // all_xz_coll_pdf: 0 - 0 : tau = x1 * x2
  //                  1 - 0 : x1 -> enters 1st pdf if ()_+ or delta distributions are present
  //                  2 - 0 : x2 -> enters 2nd pdf if ()_+ or delta distributions are present
  //                  0 - 1 : z1 -> collinear emission from a
  //                  1 - 1 : x1 / z1 -> enters 1st pdf for regular and ()_+ distributions
  //                  2 - 1 : empty - could be filled with (x1 / z1) * x2
  //                  0 - 2 : z2 -> collinear emission from b
  //                  1 - 2 : empty - could be filled with (x2 / z2) * x1
  //                  2 - 2 : x2 / z2 -> enters 1st pdf for regular and ()_+ distributions
  // should maybe replace x_pdf, z_coll, and xz_coll_pdf later...

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
*/













// temporarily !!! only until pdf evaluation for collinear contributions has been shifted to the 'list' formulation...





void observable_set::calculate_pdf_LHAPDF_CA_collinear_TSV(vector<vector<double> > & all_xz_coll_pdf){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_CA_collinear_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (!switch_TSV){return;}

  if (coll_choice == 0){cout << "parton--parton scattering CA contribution doesn't make sense!" << endl; exit(1);}
  static int initialization = 1;
  static vector<vector<vector<vector<vector<vector<double> > > > > > lhapdf_result(n_ps, vector<vector<vector<vector<vector<double> > > > > (max_dyn_fact + 1));
  static vector<vector<int> > in_collinear((*CA_collinear).size());
  static vector<vector<vector<int> > > all_pdf((*CA_collinear).size());
  static int hasPhoton = LHAPDF::hasPhoton();
  static map<int, double> charge_particle;
  static vector<double> charge2_particle(26);
  if (initialization == 1){
    fill_charge_particle(charge_particle);
    for (int i_p = 0; i_p < 26; i_p++){
      charge2_particle[i_p] = pow(charge_particle[i_p], 2);
    }
    for (int i_a = 0; i_a < n_ps; i_a++){
      for (int i_v = 0; i_v < max_dyn_fact + 1; i_v++){
	logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << "   i_v = " << i_v << endl;
	if (hasPhoton == 0){lhapdf_result[i_a][i_v].resize(n_scale_dyn_fact[i_v], vector<vector<vector<double> > > (3, vector<vector<double> > (3, vector<double> (13))));}
	else {lhapdf_result[i_a][i_v].resize(n_scale_dyn_fact[i_v], vector<vector<vector<double> > > (3, vector<vector<double> > (3, vector<double> (14))));}
      }
    }
    for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
      logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << endl;
      in_collinear[i_a] = (*CA_collinear)[i_a][0].in_collinear();
      all_pdf[i_a] = (*CA_collinear)[i_a][0].all_pdf();
    }
    initialization = 0;
    logger << LOG_DEBUG_VERBOSE << "initialization finished!" << endl;
  }

  for (int i_a = 0; i_a < n_ps; i_a++){
    for (int i_v = 0; i_v < max_dyn_fact + 1; i_v++){
      logger << LOG_DEBUG_VERBOSE << "i_a = " << i_a << "   i_v = " << i_v << endl;
      for (int i_m = 0; i_m < n_scale_dyn_fact[i_v]; i_m++){
	logger << LOG_DEBUG_VERBOSE << "value_scale_fact[" << i_a << "][" << i_v << "][" << i_m << "] = " << value_scale_fact[i_a][i_v][i_m] << endl;
	for (int i_x = 1; i_x < 3; i_x++){
	  for (int i_z = 0; i_z < i_x + 1; i_z += i_x){
	    if (hasPhoton == 0){lhapdf_result[i_a][i_v][i_m][i_x][i_z] = LHAPDF::xfx(all_xz_coll_pdf[i_x][i_z], value_scale_fact[i_a][i_v][i_m]);}
	    else {lhapdf_result[i_a][i_v][i_m][i_x][i_z] = LHAPDF::xfxphoton(all_xz_coll_pdf[i_x][i_z], value_scale_fact[i_a][i_v][i_m]);}
	    //      modify_pdf_content(lhapdf_result[i_a][i_v][i_m][i_x][i_z], *this);
	    logger << LOG_DEBUG_VERBOSE << "i_m = " << i_m << "   i_x = " << i_x << "   i_z = " << i_z << endl;
	    /*
	    for (int i_p = 0; i_p < lhapdf_result[i_a][i_v][i_m][i_x][i_z].size(); i_p++){
	      logger << LOG_DEBUG_VERBOSE << "lhapdf_result[" << i_a << "][" << i_v << "][" << i_m << "][" << i_x << "][" << i_z << "][" << i_p << "] = " << setprecision(15) << setw(23) << lhapdf_result[i_a][i_v][i_m][i_x][i_z][i_p] << endl;
	    }
	    */

	    for (int i_p = 0; i_p < lhapdf_result[i_a][i_v][i_m][i_x][i_z].size(); i_p++){
	      logger << LOG_DEBUG_VERBOSE << "pdf[" << i_a << "][" << i_v << "][" << i_m << "][" << i_x << "][" << i_z << "][" << i_p << "] = " << setprecision(15) << setw(23) << lhapdf_result[i_a][i_v][i_m][i_x][i_z][i_p] << " / " << setprecision(15) << setw(17) << all_xz_coll_pdf[i_x][i_z] << " = " << lhapdf_result[i_a][i_v][i_m][i_x][i_z][i_p] / all_xz_coll_pdf[i_x][i_z] << endl;
	      
	      if (i_z == 0){
		//		logger << LOG_DEBUG_VERBOSE << "lhapdf_result[" << i_a << "][" << i_v << "][" << i_m << "][" << i_x << "][" << i_z << "][" << i_p << "] / " << setprecision(15) << setw(23) << all_xz_coll_pdf[i_x][i_z] << " = " << lhapdf_result[i_a][i_v][i_m][i_x][i_z][i_p] / all_xz_coll_pdf[i_x][i_z] << endl;
	      }
	      else {
		//		logger << LOG_DEBUG_VERBOSE << "lhapdf_result[" << i_a << "][" << i_v << "][" << i_m << "][" << i_x << "][" << i_z << "][" << i_p << "] / " << setprecision(15) << setw(23) << all_xz_coll_pdf[i_x][i_z] << " divided by " << all_xz_coll_pdf[0][i_z] << " = " << lhapdf_result[i_a][i_v][i_m][i_x][i_z][i_p] / all_xz_coll_pdf[i_x][i_z] / all_xz_coll_pdf[0][i_z] << endl;
		//		logger << LOG_DEBUG_VERBOSE << "lhapdf_result[" << i_a << "][" << i_v << "][" << i_m << "][" << i_x << "][" << i_z << "][" << i_p << "] / " << setprecision(15) << setw(23) << all_xz_coll_pdf[i_x][i_z] << " = " << lhapdf_result[i_a][i_v][i_m][i_x][i_z][i_p] / all_xz_coll_pdf[i_x][i_z] << endl;
		logger << LOG_DEBUG_VERBOSE << "pdf[" << i_a << "][" << i_v << "][" << i_m << "][" << i_x << "][" << i_z << "][" << i_p << "] / " << setprecision(15) << setw(17) << all_xz_coll_pdf[i_x][i_z] << " divided by " << setprecision(15) << setw(17) << all_xz_coll_pdf[0][i_z] << " = " << lhapdf_result[i_a][i_v][i_m][i_x][i_z][i_p] / all_xz_coll_pdf[i_x][i_z] / all_xz_coll_pdf[0][i_z] << endl;
	      }

	    }
 	  }
	}
      }
    }
  }
  /*
  for (int i_c = 0; i_c < n_pc; i_c++){
    if (cut_ps[relation_pc_ps[i_c]] == -1){continue;}
    for (int i_z = 0; i_z < n_pz; i_z++){
      if (in_collinear[i_c][i_z] == 0){continue;}
      for (int i_v = 0; i_v < max_dyn_fact + 1; i_v++){
	if (i_v == 0 && relation_pc_ps[i_c] > 0){
	  value_pdf_factor_combination_1[i_c][i_z][i_v] = value_pdf_factor_combination_1[0][i_z][i_v];
	  value_pdf_factor_combination_2[i_c][i_z][i_v] = value_pdf_factor_combination_2[0][i_z][i_v];
	  // check
	  for (int i_m = 0; i_m < n_scale_dyn_fact[i_v]; i_m++){
	    for (int i_h = 0; i_h < 3; i_h++){
	      value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].resize(CA_combination_pdf[i_c].size(), 0.);
	      value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][i_h].resize(CA_combination_pdf[i_c].size(), 0.);
	      logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].size() = " << value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].size() << endl;
	    }
	  }
	  // !!!
	}
	else {
	  value_pdf_factor_combination_1[i_c][i_z][i_v].resize(n_scale_dyn_fact[i_v], vector<vector<double> > (3));
	  value_pdf_factor_combination_2[i_c][i_z][i_v].resize(n_scale_dyn_fact[i_v], vector<vector<double> > (3));
	  for (int i_m = 0; i_m < n_scale_dyn_fact[i_v]; i_m++){
	    for (int i_h = 0; i_h < 3; i_h++){
	      value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].resize(CA_combination_pdf[i_c].size(), 0.);
	      value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][i_h].resize(CA_combination_pdf[i_c].size(), 0.);
	      logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].size() = " << value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h].size() << endl;
	    }
	  }
	}
      }
    }
  }
*/
  /*
  logger << LOG_DEBUG_VERBOSE << "before pointer" << endl;
  for (int i_c = 0; i_c < n_pc; i_c++){
    ///    if (cut_ps[relation_pc_ps[i_c]] == -1){continue;}
    for (int i_z = 0; i_z < n_pz; i_z++){
      ///      if (in_collinear[i_c][i_z] == 0){continue;}
      for (int i_s = 0; i_s < n_set_TSV; i_s++){
	for (int i_f = 0; i_f < n_scale_fact_TSV[i_s]; i_f++){
	  for (int i_h = 0; i_h < 3; i_h++){
	    logger << LOG_DEBUG_VERBOSE << "CA_combination_pdf[i_c].size() = " << CA_combination_pdf[i_c].size() << endl;
	    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
	      /*
	      logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << "   i_z = " << i_z << "   i_s = " << i_s << "   i_f = " << i_f << "   i_h = " << i_h << "   i_i = " << i_i << endl;
	      logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_1.size() = " << value_pdf_factor_combination_1.size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_1[" << i_c << "].size() = " << value_pdf_factor_combination_1[i_c].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_1[" << i_c << "][" << i_z << "].size() = " << value_pdf_factor_combination_1[i_c][i_z].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_1[" << i_c << "][" << i_z << "][" << dynamic_scale_fact_TSV[i_s] << "].size() = " << value_pdf_factor_combination_1[i_c][i_z][dynamic_scale_fact_TSV[i_s]].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_1[" << i_c << "][" << i_z << "][" << dynamic_scale_fact_TSV[i_s] << "][" << no_value_fact_TSV[i_s][i_f] << "].size() = " << value_pdf_factor_combination_1[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "value_pdf_factor_combination_1[" << i_c << "][" << i_z << "][" << dynamic_scale_fact_TSV[i_s] << "][" << no_value_fact_TSV[i_s][i_f] << "][" << i_h << "].size() = " << value_pdf_factor_combination_1[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h].size() << endl;

	      logger << LOG_DEBUG_VERBOSE << "pointer_pdf_factor_combination_1.size() = " << pointer_pdf_factor_combination_1.size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "pointer_pdf_factor_combination_1[" << i_c << "].size() = " << pointer_pdf_factor_combination_1[i_c].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "pointer_pdf_factor_combination_1[" << i_c << "][" << i_z << "].size() = " << pointer_pdf_factor_combination_1[i_c][i_z].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "pointer_pdf_factor_combination_1[" << i_c << "][" << i_z << "][" << i_s << "].size() = " << pointer_pdf_factor_combination_1[i_c][i_z][i_s].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "pointer_pdf_factor_combination_1[" << i_c << "][" << i_z << "][" << i_s << "][" << i_f << "].size() = " << pointer_pdf_factor_combination_1[i_c][i_z][i_s][i_f].size() << endl;
	      logger << LOG_DEBUG_VERBOSE << "pointer_pdf_factor_combination_1[" << i_c << "][" << i_z << "][" << i_s << "][" << i_f << "][" << i_h << "].size() = " << pointer_pdf_factor_combination_1[i_c][i_z][i_s][i_f][i_h].size() << endl;
*/
	      /*

	      pointer_pdf_factor_combination_1[i_c][i_z][i_s][i_f][i_h][i_i] = &value_pdf_factor_combination_1[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h][i_i];
	      logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << "   i_z = " << i_z << "   i_s = " << i_s << "   i_f = " << i_f << "   i_h = " << i_h << "   i_i = " << i_i << endl;
	      pointer_pdf_factor_combination_2[i_c][i_z][i_s][i_f][i_h][i_i] = &value_pdf_factor_combination_2[i_c][i_z][dynamic_scale_fact_TSV[i_s]][no_value_fact_TSV[i_s][i_f]][i_h][i_i];
	      logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << "   i_z = " << i_z << "   i_s = " << i_s << "   i_f = " << i_f << "   i_h = " << i_h << "   i_i = " << i_i << endl;
	      *//*
	    }
	  }
	}
      }
    }
  }
*/
  logger << LOG_DEBUG << "after old pointer" << endl;

  // n_pc -> number of (dipole) contributions: 
  for (int i_c = 0; i_c < n_pc; i_c++){
    logger << LOG_DEBUG << "cut_ps[relation_pc_ps[i_c]] = " << cut_ps[relation_pc_ps[i_c]] << endl;
    if (cut_ps[relation_pc_ps[i_c]] == -1){continue;}
    // n_pz -> number of x1-x2 combinations for pdfs (like x1-x2, z1x1-x2, x1-z2x2, ...)
    for (int i_z = 0; i_z < n_pz; i_z++){
      if (in_collinear[i_c][i_z] == 0){continue;}
      logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << "   i_z = " << i_z << endl;
      for (int i_v = 0; i_v < max_dyn_fact + 1; i_v++){
	// !!! check that this is calculated only if cut_ps[...] == 0 !!!
	if (i_v == 0 && relation_pc_ps[i_c] > 0){
	  value_pdf_factor[i_c][i_z][i_v] = value_pdf_factor[0][i_z][i_v];
	  // for debugging against Sherpa
	  //	  value_pdf_factor_1[i_c][i_z][i_v] = value_pdf_factor_1[0][i_z][i_v];
	  //	  value_pdf_factor_2[i_c][i_z][i_v] = value_pdf_factor_2[0][i_z][i_v];
	  // end
	}
	else {
	  int i_z1 = 0;
	  int i_z2 = 0;
	  if      (i_z == 0){i_z1 = 0; i_z2 = 0;}
	  else if (i_z == 1){i_z1 = 1; i_z2 = 0;}
	  else if (i_z == 2){i_z1 = 0; i_z2 = 2;}

	  for (int i_m = 0; i_m < n_scale_dyn_fact[i_v]; i_m++){
	    for (int i_h = 0; i_h < 3; i_h++){value_pdf_factor[i_c][i_z][i_v][i_m][i_h] = 0.;}
	    // for debugging against Sherpa
	    for (int i_h = 0; i_h < 3; i_h++){
	      for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
		value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h][i_i] = 0.;
		value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][i_h][i_i] = 0.;
	      }
	    }
	    // end

	    // i_i counts different quarks/antiquarks in case of q->a r q->g splittings - otherwise only one component !!!
	    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){

	      double factor1 = 0.;
	      double factor2 = 0.;
	      //	      if (all_pdf[i_c][i_i][0] == -1){direction = 2;}
	      for (int i_x = 0; i_x < CA_combination_pdf[i_c][i_i].size(); i_x++){
		int direction = 1;
		if (CA_combination_pdf[i_c][i_i][i_x][0] == -1){direction = 2;}

		// ??? relation_pc_ps[i_c]
		// all_pdf[i_c][i_i][1] instead of CA_combination_pdf[i_c][i_i][i_x][1] ???
		if (all_pdf[i_c][i_i][1] == 7){factor1 = lhapdf_result[relation_pc_ps[i_c]][i_v][i_m][1][i_z1][13];}
		else if (coll_choice == 1){factor1 = lhapdf_result[relation_pc_ps[i_c]][i_v][i_m][1][i_z1][6 + CA_combination_pdf[i_c][i_i][i_x][1]];}
		else if (coll_choice == 2){
		  if      (CA_combination_pdf[i_c][i_i][i_x][0] ==  1){factor1 = lhapdf_result[relation_pc_ps[i_c]][i_v][i_m][1][i_z1][6 + CA_combination_pdf[i_c][i_i][i_x][1]];}
		  else if (CA_combination_pdf[i_c][i_i][i_x][0] == -1){factor1 = lhapdf_result[relation_pc_ps[i_c]][i_v][i_m][1][i_z1][6 - CA_combination_pdf[i_c][i_i][i_x][1]];}
		}
		if (CA_combination_pdf[i_c][i_i][i_x][2] == 7){factor2 = lhapdf_result[relation_pc_ps[i_c]][i_v][i_m][2][i_z2][13];}
		else if (coll_choice == 1){factor2 = lhapdf_result[relation_pc_ps[i_c]][i_v][i_m][2][i_z2][6 + CA_combination_pdf[i_c][i_i][i_x][2]];}
		else if (coll_choice == 2){
		  if      (CA_combination_pdf[i_c][i_i][i_x][0] ==  1){factor2 = lhapdf_result[relation_pc_ps[i_c]][i_v][i_m][2][i_z2][6 - CA_combination_pdf[i_c][i_i][i_x][2]];}
		  else if (CA_combination_pdf[i_c][i_i][i_x][0] == -1){factor2 = lhapdf_result[relation_pc_ps[i_c]][i_v][i_m][2][i_z2][6 + CA_combination_pdf[i_c][i_i][i_x][2]];}
		}
		double factor = factor1 * factor2;
		logger << LOG_DEBUG_VERBOSE << "factor1 = " << factor1 << endl;
		logger << LOG_DEBUG_VERBOSE << "factor2 = " << factor2 << endl;
		logger << LOG_DEBUG_VERBOSE << "CA_Q2f[" << i_c << "][" << i_x << "] = " << CA_Q2f[i_c][i_i][i_x] << endl;

		if ((*CA_collinear)[i_c][i_i].type_correction() == 2){factor *= CA_Q2f[i_c][i_i][i_x];}// charge factors !!!
		value_pdf_factor[i_c][i_z][i_v][i_m][direction] += factor;

		// for debugging against Sherpa
		double temp_Sherpa_factor1 = factor1;
		// factor1 always from beam 1 (not necessarily the one with collinear emission)
		double temp_Sherpa_factor2 = factor2;
		// factor2 always from beam 2 (not necessarily the one with collinear emission)
		if ((*CA_collinear)[i_c][i_i].type_correction() == 2){
		  if ((*CA_collinear)[i_c][i_i].no_emitter() == 1){temp_Sherpa_factor1 *= CA_Q2f[i_c][i_i][i_x];}
		  else if ((*CA_collinear)[i_c][i_i].no_emitter() == 2){temp_Sherpa_factor2 *= CA_Q2f[i_c][i_i][i_x];}
		}// charge factors !!!

		//		logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << "   direction = " << direction << "   temp_Sherpa_factor1 = " << temp_Sherpa_factor1 << "   temp_Sherpa_factor2 = " << temp_Sherpa_factor2 << endl;
		//		logger << LOG_DEBUG_VERBOSE << "(*CA_collinear)[" << i_c << "][" << 0 << "].no_emitter() = " << (*CA_collinear)[i_c][0].no_emitter() << endl;
		if ((*CA_collinear)[i_c][0].no_emitter() == 1){
		  value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][direction][i_i] += temp_Sherpa_factor1;
		  value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][direction][i_i] = temp_Sherpa_factor2;
		}
		else if ((*CA_collinear)[i_c][0].no_emitter() == 2){
		  value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][direction][i_i] = temp_Sherpa_factor1;
		  value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][direction][i_i] += temp_Sherpa_factor2;
		}
		// end


		//		if (){
		//		cout << "(*CA_collinear)[" << i_c << "][" << i_i << "].type_correction() = " << (*CA_collinear)[i_c][i_i].type_correction() << endl;
		//		cout << "(*CA_collinear)[" << i_c << "][" << i_i << "].no_emitter() = " << (*CA_collinear)[i_c][i_i].no_emitter() << endl;
		//		cout << "(*CA_collinear)[" << i_c << "][" << i_i << "].type_parton()[" << (*CA_collinear)[i_c][i_i].no_emitter() << "] = " << (*CA_collinear)[i_c][i_i].type_parton()[(*CA_collinear)[i_c][i_i].no_emitter()] << endl;
		

		/*
		if ((*CA_collinear)[i_c][i_i].type_correction() == 2 && (*CA_collinear)[i_c][i_i].type_parton()[(*CA_collinear)[i_c][i_i].no_emitter()] == 22){
		  value_pdf_factor[i_c][i_z][i_v][i_m][direction] *= charge2_particle[abs(CA_combination_pdf[i_c][i_i][i_x][(*CA_collinear)[i_c][i_i].no_emitter()])];// charge factors !!!
		}
		*/
		  /*
		  cout << "CA_combination_pdf[i_c = " << i_c << "][i_i = " << i_i << "][i_x = " << i_x << "] = ";
		  for (int i = 0; i < 3; i++){
		    cout << CA_combination_pdf[i_c][i_i][i_x][i] << "   ";
		  }
		  cout << endl;
		  cout << "charge2_particle[abs(CA_combination_pdf[i_c = " << i_c << "][i_i = " << i_i << "][i_x = " << i_x << "][(*CA_collinear)[i_a].no_emitter()] = " << (*CA_collinear)[i_c][i_i].no_emitter() << "]) = " << charge2_particle[abs(CA_combination_pdf[i_c][i_i][i_x][(*CA_collinear)[i_c][i_i].no_emitter()])] << endl;
		  */		//		}
	      }
	    }
	    
	    for (int i_h = 1; i_h < 3; i_h++){
	      // divide by the arguments of the pdfs, as LHAPDF delivers x*f(x)
	      if (i_z == 0){value_pdf_factor[i_c][i_z][i_v][i_m][i_h] = value_pdf_factor[i_c][i_z][i_v][i_m][i_h] / all_xz_coll_pdf[0][0];}
	      else {value_pdf_factor[i_c][i_z][i_v][i_m][i_h] = value_pdf_factor[i_c][i_z][i_v][i_m][i_h] / all_xz_coll_pdf[i_z % 2 + 1][i_z];}
	    }
	    value_pdf_factor[i_c][i_z][i_v][i_m][0] = value_pdf_factor[i_c][i_z][i_v][i_m][1] + value_pdf_factor[i_c][i_z][i_v][i_m][2];

	    // for debugging against Sherpa
	    for (int i_h = 1; i_h < 3; i_h++){
	      for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
		// divide only by x_pdf[1,2] (LHAPDF delivers x*fx); with collinear emission, the result is pdf / z_coll[1,2] (corresponds to 1 / z_coll[1,2] factor absorbsed into g_z_coll[1,2]
		value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h][i_i] = value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][i_h][i_i] / all_xz_coll_pdf[1][0];
		value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][i_h][i_i] = value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][i_h][i_i] / all_xz_coll_pdf[2][0];
	      }
	    }
	    for (int i_i = 0; i_i < CA_combination_pdf[i_c].size(); i_i++){
	      value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][0][i_i] = value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][1][i_i] + value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][2][i_i];
	      value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][0][i_i] = value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][1][i_i] + value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][2][i_i];

















	      logger << LOG_DEBUG_POINT << "value_pdf_factor_combination_1[" << i_c << "][" << i_z << "][" << i_v << "][" << i_m << "][0][" << i_i << "] = " << value_pdf_factor_combination_1[i_c][i_z][i_v][i_m][0][i_i] << endl;
	      logger << LOG_DEBUG_POINT << "value_pdf_factor_combination_2[" << i_c << "][" << i_z << "][" << i_v << "][" << i_m << "][0][" << i_i << "] = " << value_pdf_factor_combination_2[i_c][i_z][i_v][i_m][0][i_i] << endl;

	    }
	    // end

	    //	    cout << "value_pdf_factor[" << i_c << "][" << i_z << "][" << i_v << "][" << i_m << "][0]  = " << value_pdf_factor[i_c][i_z][i_v][i_m][0]  << endl;
	  }
	}
      }
    }
  }

  // all_xz_coll_pdf: 0 - 0 : tau = x1 * x2
  //                  1 - 0 : x1 -> enters 1st pdf if ()_+ or delta distributions are present
  //                  2 - 0 : x2 -> enters 2nd pdf if ()_+ or delta distributions are present
  //                  0 - 1 : z1 -> collinear emission from a
  //                  1 - 1 : x1 / z1 -> enters 1st pdf for regular and ()_+ distributions
  //                  2 - 1 : empty - could be filled with (x1 / z1) * x2
  //                  0 - 2 : z2 -> collinear emission from b
  //                  1 - 2 : empty - could be filled with (x2 / z2) * x1
  //                  2 - 2 : x2 / z2 -> enters 1st pdf for regular and ()_+ distributions
  // should maybe replace x_pdf, z_coll, and xz_coll_pdf later...

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void observable_set::calculate_pdf_LHAPDF_CA_collinear_CV(vector<vector<double> > & all_xz_coll_pdf){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_CA_collinear_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    if (value_mu_fact[sd].size() == 0){continue;}
    //    CA_value_pdf_factor[sd].resize(value_mu_fact[sd].size(), vector<vector<vector<double> > > ((*CA_collinear)[0][0].all_pdf().size())); // ???
    //    cout << "collinear[0][0].all_pdf().size() = " << (*CA_collinear)[0][0].all_pdf().size() << endl;
    for (int ss = 0; ss < value_mu_fact[sd].size(); ss++){
      logger << LOG_DEBUG_VERBOSE << "sd = " << sd << "   ss = " << ss << endl;
      logger << LOG_DEBUG_VERBOSE << "value_mu_fact.size() = " << value_mu_fact.size() << endl;
      logger << LOG_DEBUG_VERBOSE << "value_mu_fact[" << sd << "].size() = " << value_mu_fact[sd].size() << endl;
      logger << LOG_DEBUG_VERBOSE << "value_mu_fact[" << sd << "][" << ss << "] = " << value_mu_fact[sd][ss] << endl;
      logger << LOG_DEBUG_VERBOSE << "CA_value_pdf_factor.size() = " << CA_value_pdf_factor.size() << endl;
      logger << LOG_DEBUG_VERBOSE << "CA_value_pdf_factor[" << sd << "].size() = " << CA_value_pdf_factor[sd].size() << endl;
      logger << LOG_DEBUG_VERBOSE << "CA_value_pdf_factor[" << sd << "][" << ss << "].size() = " << CA_value_pdf_factor[sd][ss].size() << endl;

      if (coll_choice == 0){cout << "parton--parton scattering CA contribution doesn't make sense!" << endl; exit(1);}
      else {calculate_pdf_LHAPDF_scale_CA_collinear_CV(all_xz_coll_pdf, value_mu_fact[sd][ss], CA_value_pdf_factor[sd][ss]);}
      //      if (coll_choice == 1){calculate_pdf_LHAPDF_scale_collinear_CV((*CA_collinear), all_xz_coll_pdf, sd, ss, value_mu_fact[sd][ss], CA_value_pdf_factor[sd][ss], oset);}
      //      else if (coll_choice == 2){calculate_pdf_LHAPDF_scale_collinear_CV((*CA_collinear), all_xz_coll_pdf, sd, ss, value_mu_fact[sd][ss], CA_value_pdf_factor[sd][ss], oset);} // should work with the same routine !!!
      // maybe include "calculate_pdf_LHAPDF_CA_LHC_collinear" here:
      // mu_fact -> value_mu_fact[sd][ss]
      // all_pdf_factor -> CA_value_pdf_factor[sd][ss]
      // old //	  else if (coll_choice == 2){calculate_pdf_LHAPDF_CA_Tevatron(all_pdf, x_pdf, CA_value_mu_fact[ib][sd][ss], CA_value_pdf_factor[ib][sd][ss]);}
    }
  }
  // CA_pdf_factor: [ia][iz1][iz2][0-1-2]
  // CA_pdf_factor_CV: [ia][s][iz1][iz2][0-1-2]

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::calculate_pdf_LHAPDF_scale_CA_collinear_CV(vector<vector<double> > & all_xz_coll_pdf, double & mu_fact, vector<vector<vector<double> > > & all_pdf_factor){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_scale_CA_collinear_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static vector<vector<vector<double> > > lhapdf_result;
  static vector<vector<int> > in_collinear((*CA_collinear).size());
  static vector<vector<vector<int> > > all_pdf((*CA_collinear).size());
  static int hasPhoton = LHAPDF::hasPhoton();
  static map<int, double> charge_particle;
  static vector<double> charge2_particle(26);
  //  static vector<vector<vector<double> > > empty_all_pdf_factor((*CA_collinear)[0][0].all_pdf().size(), vector<vector<double> > (3, vector<double> (3, 0.))); // ???
  if (initialization == 1){
    fill_charge_particle(charge_particle);
    for (int i_p = 0; i_p < 26; i_p++){
      charge2_particle[i_p] = pow(charge_particle[i_p], 2);
    }
    //    all_pdf_factor = empty_all_pdf_factor;
    if (hasPhoton == 0){lhapdf_result.resize(3, vector<vector<double> > (3, vector<double> (13)));}
    else {lhapdf_result.resize(3, vector<vector<double> > (3, vector<double> (14)));}
    for (int i_a = 0; i_a < (*CA_collinear).size(); i_a++){
      in_collinear[i_a] = (*CA_collinear)[i_a][0].in_collinear();
      all_pdf[i_a] = (*CA_collinear)[i_a][0].all_pdf();
    }
    initialization = 0;
  }

  for (int i_x = 1; i_x < 3; i_x++){
    for (int i_z = 0; i_z < i_x + 1; i_z += i_x){
      if (hasPhoton == 0){lhapdf_result[i_x][i_z] = LHAPDF::xfx(all_xz_coll_pdf[i_x][i_z], mu_fact);}
      else {lhapdf_result[i_x][i_z] = LHAPDF::xfxphoton(all_xz_coll_pdf[i_x][i_z], mu_fact);}
      //      modify_pdf_content(lhapdf_result[i_x][i_z], *this);
      //      for (int i_p = 0; i_p < lhapdf_result[i_x][i_z].size(); i_p++){cout << "lhapdf_result[" << i_x << "][" << i_z << "][" << i_p << "] = " << lhapdf_result[i_x][i_z][i_p] << endl;}
    }
  }
  static vector<vector<vector<double> > > empty_all_pdf_factor(all_pdf_factor.size(), vector<vector<double> > (3, vector<double> (3, 0.))); // ???
  all_pdf_factor = empty_all_pdf_factor;
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_i = 0; i_i < all_pdf[i_c].size(); i_i++){
      for (int i_z = 0; i_z < 3; i_z++){
	if (in_collinear[i_c][i_z] == 1){
	  //	  cout << "i_i = " << i_i << "   all_pdf_factor[" << i_c << "][" << i_z << "][1] = " << setw(25) << setprecision(16) << all_pdf_factor[i_c][i_z][1] << "   all_pdf_factor[" << i_c << "][" << i_z << "][2] = " << setw(25) << setprecision(16) << all_pdf_factor[i_c][i_z][2] << endl;
	  int i_z1;
	  int i_z2;
	  if      (i_z == 0){i_z1 = 0; i_z2 = 0;}
	  else if (i_z == 1){i_z1 = 1; i_z2 = 0;}
	  else if (i_z == 2){i_z1 = 0; i_z2 = 2;}
	  // new version (all (q/q~ -> g splittings combined) start
	  double factor1 = 0.;
	  double factor2 = 0.;
	  //	  int direction = 1;
	  //	  if (all_pdf[i_c][i_i][0] == -1){direction = 2;}
	  for (int i_x = 0; i_x < CA_combination_pdf[i_c][i_i].size(); i_x++){
	    int direction = 1;
	    if (CA_combination_pdf[i_c][i_i][i_x][0] == -1){direction = 2;}

	    if (CA_combination_pdf[i_c][i_i][i_x][1] == 7){factor1 = lhapdf_result[1][i_z1][13];}
	    else if (coll_choice == 1){factor1 = lhapdf_result[1][i_z1][6 + CA_combination_pdf[i_c][i_i][i_x][1]];}
	    else if (coll_choice == 2){
	      if      (CA_combination_pdf[i_c][i_i][i_x][0] ==  1){factor1 = lhapdf_result[1][i_z1][6 + CA_combination_pdf[i_c][i_i][i_x][1]];}
	      else if (CA_combination_pdf[i_c][i_i][i_x][0] == -1){factor1 = lhapdf_result[1][i_z1][6 - CA_combination_pdf[i_c][i_i][i_x][1]];}
	    }
	    if (CA_combination_pdf[i_c][i_i][i_x][2] == 7){factor2 = lhapdf_result[2][i_z2][13];}
	    else if (coll_choice == 1){factor2 = lhapdf_result[2][i_z2][6 + CA_combination_pdf[i_c][i_i][i_x][2]];}
	    else if (coll_choice == 2){
	      if      (CA_combination_pdf[i_c][i_i][i_x][0] ==  1){factor2 = lhapdf_result[2][i_z2][6 - CA_combination_pdf[i_c][i_i][i_x][2]];}
	      else if (CA_combination_pdf[i_c][i_i][i_x][0] == -1){factor2 = lhapdf_result[2][i_z2][6 + CA_combination_pdf[i_c][i_i][i_x][2]];}
	    }

	    double factor = factor1 * factor2;
	    if ((*CA_collinear)[i_c][i_i].type_correction() == 2){factor *= CA_Q2f[i_c][i_i][i_x];}// charge factors !!!
	    all_pdf_factor[i_c][i_z][direction] += factor;


	    //	    cout << "(*CA_collinear)[" << i_c << "][" << i_i << "].type_correction() = " << (*CA_collinear)[i_c][i_i].type_correction() << endl;
	    //	    cout << "(*CA_collinear)[" << i_c << "][" << i_i << "].no_emitter() = " << (*CA_collinear)[i_c][i_i].no_emitter() << endl;
	    //	    cout << "(*CA_collinear)[" << i_c << "][" << i_i << "].type_parton()[" << (*CA_collinear)[i_c][i_i].no_emitter() << "] = " << (*CA_collinear)[i_c][i_i].type_parton()[(*CA_collinear)[i_c][i_i].no_emitter()] << endl;
	    	  
	    //	    if ((*CA_collinear)[i_c][i_i].type_correction() == 2){all_pdf_factor[i_c][i_z][direction] *= CA_Q2f[i_c][i_i][i_x];}// charge factors !!!
	    /*
	    if ((*CA_collinear)[i_c][i_i].type_correction() == 2 && (*CA_collinear)[i_c][i_i].type_parton()[(*CA_collinear)[i_c][i_i].no_emitter()] == 22){
	      // charge factors !!!
	      all_pdf_factor[i_c][i_z][direction] *= charge2_particle[abs(CA_combination_pdf[i_c][i_i][i_x][(*CA_collinear)[i_c][i_i].no_emitter()])];
	    }
	    */
	      /*
	      cout << "CA_combination_pdf[i_c = " << i_c << "][i_i = " << i_i << "][i_x = " << i_x << "] = ";
	      for (int i = 0; i < 3; i++){
		cout << CA_combination_pdf[i_c][i_i][i_x][i] << "   ";
	      }
	      cout << endl;
	      
	      cout << "charge2_particle[abs(CA_combination_pdf[i_c = " << i_c << "][i_i = " << i_i << "][i_x = " << i_x << "][(*CA_collinear)[i_c].no_emitter()] = " << (*CA_collinear)[i_c][i_i].no_emitter() << "]) = " << charge2_particle[abs(CA_combination_pdf[i_c][i_i][i_x][(*CA_collinear)[i_c][i_i].no_emitter()])] << endl;
	      */
	  }

	  
	  /*
	  if (all_pdf[i_c][i_i][1] == 10){
	    for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	      if (i_q != 0){factor1 += lhapdf_result[1][i_z1][6 + i_q];}
	    }
	  }
	  else {
	    if (all_pdf[i_c][i_i][1] == 7){factor1 = lhapdf_result[1][i_z1][13];}
	    else if (coll_choice == 1){factor1 = lhapdf_result[1][i_z1][6 + all_pdf[i_c][i_i][1]];}
	    else if (coll_choice == 2){
	      if      (all_pdf[i_c][i_i][0] ==  1){factor1 = lhapdf_result[1][i_z1][6 + all_pdf[i_c][i_i][1]];}
	      else if (all_pdf[i_c][i_i][0] == -1){factor1 = lhapdf_result[1][i_z1][6 - all_pdf[i_c][i_i][1]];}
	    }
	  }
	  if (all_pdf[i_c][i_i][2] == 10){
	    for (int i_q = -N_f_active; i_q < N_f_active + 1; i_q++){
	      if (i_q != 0){factor2 += lhapdf_result[2][i_z2][6 + i_q];}
	    }
	  }
	  else {
	    if (all_pdf[i_c][i_i][2] == 7){factor2 = lhapdf_result[2][i_z2][13];}
	    else if (coll_choice == 1){factor2 = lhapdf_result[2][i_z2][6 + all_pdf[i_c][i_i][2]];}
	    else if (coll_choice == 2){
	      if      (all_pdf[i_c][i_i][0] ==  1){factor2 = lhapdf_result[2][i_z2][6 - all_pdf[i_c][i_i][2]];}
	      else if (all_pdf[i_c][i_i][0] == -1){factor2 = lhapdf_result[2][i_z2][6 + all_pdf[i_c][i_i][2]];}
	    }
	  }
	  int direction = 1;
	  if (all_pdf[i_c][i_i][0] == -1){direction = 2;}
	  all_pdf_factor[i_c][i_z][direction] += factor1 * factor2;
	  */
	  
	  // new version (all (q/q~ -> g splittings combined) end
	}
      }
    }
  }
  for (int i_c = 0; i_c < (*CA_collinear).size(); i_c++){
    for (int i_z = 0; i_z < 3; i_z++){
      if (in_collinear[i_c][i_z] == 0){continue;}
      for (int i_y = 1; i_y < 3; i_y++){
	if (i_z == 0){all_pdf_factor[i_c][i_z][i_y] = all_pdf_factor[i_c][i_z][i_y] / all_xz_coll_pdf[0][0];}
	else {all_pdf_factor[i_c][i_z][i_y] = all_pdf_factor[i_c][i_z][i_y] / all_xz_coll_pdf[i_z % 2 + 1][i_z];}
      }
      all_pdf_factor[i_c][i_z][0] = all_pdf_factor[i_c][i_z][1] + all_pdf_factor[i_c][i_z][2];
    }
  }
  
  // all_xz_coll_pdf: 0 - 0 : tau = x1 * x2
  //                  1 - 0 : x1 -> enters 1st pdf if ()_+ or delta distributions are present
  //                  2 - 0 : x2 -> enters 2nd pdf if ()_+ or delta distributions are present
  //                  0 - 1 : z1 -> collinear emission from a
  //                  1 - 1 : x1 / z1 -> enters 1st pdf for regular and ()_+ distributions
  //                  2 - 1 : empty - could be filled with (x1 / z1) * x2
  //                  0 - 2 : z2 -> collinear emission from b
  //                  1 - 2 : empty - could be filled with (x2 / z2) * x1
  //                  2 - 2 : x2 / z2 -> enters 1st pdf for regular and ()_+ distributions
  // should maybe replace x_pdf, z_coll, and xz_coll_pdf later...

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
  }




void observable_set::calculate_pdf_LHAPDF_TSV(){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_TSV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  if (!switch_TSV){return;}

  static int initialization = 1;
  static int hasPhoton = LHAPDF::hasPhoton();
  static vector<vector<double> > lhapdf_result;
  if (initialization == 1){
    if (hasPhoton == 0){lhapdf_result.resize(3, vector<double> (13));}
    else {lhapdf_result.resize(3, vector<double> (14));}
    initialization = 0;
  }
  for (int i_c = 0; i_c < n_pc; i_c++){
    if (cut_ps[relation_pc_ps[i_c]] != -1){
      for (int i_v = 0; i_v < max_dyn_fact + 1; i_v++){
	if (i_v == 0 && i_c > first_non_cut_ps){value_pdf_factor[i_c][0][i_v] = value_pdf_factor[first_non_cut_ps][0][i_v];}
	else {
	  for (int i_m = 0; i_m < n_scale_dyn_fact[i_v]; i_m++){
	    for (int i_h = 1; i_h < 3; i_h++){
	      value_pdf_factor[i_c][0][i_v][i_m][i_h] = 0.;
	      if (hasPhoton == 0){lhapdf_result[i_h] = LHAPDF::xfx(x_pdf[i_h], value_scale_fact[i_c][i_v][i_m]);}
	      else {lhapdf_result[i_h] = LHAPDF::xfxphoton(x_pdf[i_h], value_scale_fact[i_c][i_v][i_m]);}
	      //	      if (pdf_content_modify != -1){modify_pdf_content(lhapdf_result[i_h], *this);}
	    }
	    for (int i_p = 0; i_p < combination_pdf.size(); i_p++){
	      if (coll_choice == 1){
		if      (combination_pdf[i_p][0] ==  1){value_pdf_factor[i_c][0][i_v][i_m][1] += lhapdf_result[1][6 + combination_pdf[i_p][1]] * lhapdf_result[2][6 + combination_pdf[i_p][2]];}
		else if (combination_pdf[i_p][0] == -1){value_pdf_factor[i_c][0][i_v][i_m][2] += lhapdf_result[1][6 + combination_pdf[i_p][1]] * lhapdf_result[2][6 + combination_pdf[i_p][2]];}
		else {cout << "No valid pdf selection!" << endl; exit(1);}
	      }
	      else if (coll_choice == 2){
		if      (combination_pdf[i_p][0] ==  1){value_pdf_factor[i_c][0][i_v][i_m][1] += lhapdf_result[1][6 + combination_pdf[i_p][1]] * lhapdf_result[2][6 - combination_pdf[i_p][2]];}
		else if (combination_pdf[i_p][0] == -1){value_pdf_factor[i_c][0][i_v][i_m][2] += lhapdf_result[1][6 - combination_pdf[i_p][1]] * lhapdf_result[2][6 + combination_pdf[i_p][2]];}
		else {cout << "No valid pdf selection!" << endl; exit(1);}
	      }
	    }
	    for (int i_h = 1; i_h < 3; i_h++){value_pdf_factor[i_c][0][i_v][i_m][i_h] = value_pdf_factor[i_c][0][i_v][i_m][i_h] / x_pdf[0];}
	    value_pdf_factor[i_c][0][i_v][i_m][0] = value_pdf_factor[i_c][0][i_v][i_m][1] + value_pdf_factor[i_c][0][i_v][i_m][2];
	  }
	}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::calculate_pdf_LHAPDF_CV(){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<vector<vector<double> > > value_pdf_factor(value_mu_fact.size());
  for (int sd = 0; sd < value_mu_fact.size(); sd++){
    value_pdf_factor[sd].resize(value_mu_fact[sd].size(), vector<double> (3));
    for (int ss = 0; ss < value_mu_fact[sd].size(); ss++){
      if (coll_choice == 0){value_pdf_factor[sd][ss][1] = 1.; value_pdf_factor[sd][ss][2] = 0.; value_pdf_factor[sd][ss][0] = value_pdf_factor[sd][ss][1];}
      else {calculate_pdf_LHAPDF_scale_CV(value_mu_fact[sd][ss], value_pdf_factor[sd][ss]);}
      //      else if (coll_choice == 1){calculate_pdf_LHAPDF_LHC_CV(combination_pdf, x_pdf, value_mu_fact[sd][ss], value_pdf_factor[sd][ss], oset);}
      //      else if (coll_choice == 2){calculate_pdf_LHAPDF_Tevatron_CV(combination_pdf, x_pdf, value_mu_fact[sd][ss], value_pdf_factor[sd][ss], oset);}
    }
  }
  pdf_factor = value_pdf_factor[dynamic_scale][map_value_scale_fact];
  if (switch_CV){for (int s = 0; s < n_scales_CV; s++){pdf_factor_CV[s] = value_pdf_factor[dynamic_scale_CV][map_value_scale_fact_CV[s]];}}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::calculate_pdf_LHAPDF_RA_CV(){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_RA_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<vector<vector<vector<double> > > > RA_value_pdf_factor(RA_value_mu_fact.size());
  for (int ia = 0; ia < RA_value_pdf_factor.size(); ia++){
    RA_value_pdf_factor[ia].resize(RA_value_mu_fact[ia].size());
    for (int sd = 0; sd < RA_value_mu_fact[ia].size(); sd++){
      if (sd == 0 && ia > 0){
	RA_value_pdf_factor[ia][0] = RA_value_pdf_factor[0][0];
      }
      else {
	RA_value_pdf_factor[ia][sd].resize(RA_value_mu_fact[ia][sd].size(), vector<double> (3));
	for (int ss = 0; ss < RA_value_mu_fact[ia][sd].size(); ss++){
	  if (coll_choice == 0){RA_value_pdf_factor[ia][sd][ss][1] = 1.; RA_value_pdf_factor[ia][sd][ss][2] = 0.; RA_value_pdf_factor[ia][sd][ss][0] = RA_value_pdf_factor[ia][sd][ss][1];}
	  else {calculate_pdf_LHAPDF_scale_CV(RA_value_mu_fact[ia][sd][ss], RA_value_pdf_factor[ia][sd][ss]);}
	  //	  else if (coll_choice == 1){calculate_pdf_LHAPDF_LHC_CV(combination_pdf, x_pdf, RA_value_mu_fact[ia][sd][ss], RA_value_pdf_factor[ia][sd][ss], oset);}
	  //	  else if (coll_choice == 2){calculate_pdf_LHAPDF_Tevatron_CV(combination_pdf, x_pdf, RA_value_mu_fact[ia][sd][ss], RA_value_pdf_factor[ia][sd][ss], oset);}
	}
      }
    }
    RA_pdf_factor[ia] = RA_value_pdf_factor[ia][dynamic_scale][map_value_scale_fact];
    if (switch_CV){for (int s = 0; s < n_scales_CV; s++){RA_pdf_factor_CV[ia][s] = RA_value_pdf_factor[ia][dynamic_scale_CV][map_value_scale_fact_CV[s]];}}
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::calculate_pdf_LHAPDF_scale_CV(double & mu_fact, vector<double> & pdf_factor){
  static Logger logger("observable_set::calculate_pdf_LHAPDF_scale_CV");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  static int initialization = 1;
  static int hasPhoton = LHAPDF::hasPhoton();
  static vector<vector<double> > lhapdf_result;
  if (initialization == 1){
    if (hasPhoton == 0){lhapdf_result.resize(3, vector<double> (13));}
    else {lhapdf_result.resize(3, vector<double> (14));}
    initialization = 0;
  }
  for (int i_h = 1; i_h < 3; i_h++){
    if (hasPhoton == 0){lhapdf_result[i_h] = LHAPDF::xfx(x_pdf[i_h], mu_fact);}
    else {lhapdf_result[i_h] = LHAPDF::xfxphoton(x_pdf[i_h], mu_fact);}
    //    if (pdf_content_modify != -1){modify_pdf_content(lhapdf_result[i_h], *this);}
  }
  pdf_factor.resize(3, 0.);
  for (int i = 0; i < combination_pdf.size(); i++){
    if (coll_choice == 1){
      if      (combination_pdf[i][0] ==  1){pdf_factor[1] += lhapdf_result[1][6 + combination_pdf[i][1]] * lhapdf_result[2][6 + combination_pdf[i][2]];}
      else if (combination_pdf[i][0] == -1){pdf_factor[2] += lhapdf_result[1][6 + combination_pdf[i][1]] * lhapdf_result[2][6 + combination_pdf[i][2]];}
      else {cout << "No valid pdf selection!" << endl; exit(1);}
    }
    else if (coll_choice == 2){
      if      (combination_pdf[i][0] ==  1){pdf_factor[1] += lhapdf_result[1][6 + combination_pdf[i][1]] * lhapdf_result[2][6 - combination_pdf[i][2]];}
      else if (combination_pdf[i][0] == -1){pdf_factor[2] += lhapdf_result[1][6 - combination_pdf[i][1]] * lhapdf_result[2][6 + combination_pdf[i][2]];}
      else {cout << "No valid pdf selection!" << endl; exit(1);}
    }
    else {cout << "No valid pdf selection!" << endl; exit(1);}
  }
  pdf_factor[0] = pdf_factor[1] + pdf_factor[2];
  for (int i_x = 0; i_x < 3; i_x++){
    pdf_factor[i_x] = pdf_factor[i_x] / x_pdf[0];
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
