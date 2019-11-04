#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
//observable_set::
void determine_selection_content_particle(vector<int> & content_list, vector<string> & content_selection, vector<string> & content_disable){
  //  static Logger logger("observable_set::determine_selection_content_particle");
  static Logger logger("determine_selection_content_particle");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  //  logger << LOG_DEBUG << "N_f_active = " << N_f_active << endl;
  for (int i_p = 0; i_p < content_selection.size(); i_p++){
    logger << LOG_DEBUG << "content_selection[" << i_p << "] = " << content_selection[i_p] << endl;
  }
  for (int i_d = 0; i_d < content_disable.size(); i_d++){
    logger << LOG_DEBUG << "content_disable[" << i_d << "] = " << content_disable[i_d] << endl;
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
  content_p[12] = 21;
  content_p[13] = 22;

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


  vector<int> content_L(6);
  content_L[0] = 11;
  content_L[1] = 13;
  content_L[2] = 15;
  content_L[3] = -11;
  content_L[4] = -13;
  content_L[5] = -15;
  vector<int> content_l(3);
  content_l[0] = 11;
  content_l[1] = 13;
  content_l[2] = 15;
  vector<int> content_lx(3);
  content_lx[0] = -11;
  content_lx[1] = -13;
  content_lx[2] = -15;
  vector<int> content_E(2);
  content_E[0] = 11;
  content_E[1] = -11;
  vector<int> content_M(2);
  content_M[0] = 13;
  content_M[1] = -13;
  vector<int> content_Y(2);
  content_Y[0] = 15;
  content_Y[1] = -15;
  vector<int> content_e(1);
  content_e[0] = 11;
  vector<int> content_m(1);
  content_m[0] = -11;
  vector<int> content_y(1);
  content_y[0] = 13;
  vector<int> content_ex(1);
  content_ex[0] = -13;
  vector<int> content_mx(1);
  content_mx[0] = -15;
  vector<int> content_yx(1);
  content_yx[0] = -15;


  vector<int> content_N(6);
  content_N[0] = 12;
  content_N[1] = 14;
  content_N[2] = 16;
  content_N[3] = -12;
  content_N[4] = -14;
  content_N[5] = -16;
  vector<int> content_n(3);
  content_n[0] = 12;
  content_n[1] = 14;
  content_n[2] = 16;
  vector<int> content_nx(3);
  content_nx[0] = -12;
  content_nx[1] = -14;
  content_nx[2] = -16;
  vector<int> content_NE(2);
  content_NE[0] = 12;
  content_NE[1] = -12;
  vector<int> content_NM(2);
  content_NM[0] = 14;
  content_NM[1] = -14;
  vector<int> content_NY(2);
  content_NY[0] = 16;
  content_NY[1] = -16;
  vector<int> content_ne(1);
  content_ne[0] = 12;
  vector<int> content_nm(1);
  content_nm[0] = 14;
  vector<int> content_ny(1);
  content_ny[0] = 16;
  vector<int> content_nex(1);
  content_nex[0] = -12;
  vector<int> content_nmx(1);
  content_nmx[0] = -14;
  vector<int> content_nyx(1);
  content_nyx[0] = -16;


  vector<int> content_g(1);
  content_g[0] = 21;
  vector<int> content_a(1);
  content_a[0] = 22;
  vector<int> content_z(1);
  content_z[0] = 23;
  vector<int> content_W(2);
  content_W[0] = 24;
  content_W[0] = -24;
  vector<int> content_w(1);
  content_w[0] = 24;
  vector<int> content_wx(1);
  content_wx[0] = -24;
  vector<int> content_h(1);
  content_h[0] = 25;


  vector<vector<int> > content_parton(53);
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

  content_parton[24] = content_L;
  content_parton[25] = content_l;
  content_parton[26] = content_lx;
  content_parton[27] = content_E;
  content_parton[28] = content_M;
  content_parton[29] = content_Y;
  content_parton[30] = content_e;
  content_parton[31] = content_m;
  content_parton[32] = content_y;
  content_parton[33] = content_ex;
  content_parton[34] = content_mx;
  content_parton[35] = content_yx;

  content_parton[36] = content_N;
  content_parton[37] = content_n;
  content_parton[38] = content_nx;
  content_parton[39] = content_NE;
  content_parton[40] = content_NM;
  content_parton[41] = content_NY;
  content_parton[42] = content_ne;
  content_parton[43] = content_nm;
  content_parton[44] = content_ny;
  content_parton[45] = content_nex;
  content_parton[46] = content_nmx;
  content_parton[47] = content_nyx;

  content_parton[48] = content_z;

  content_parton[49] = content_W;
  content_parton[50] = content_w;
  content_parton[51] = content_wx;

  content_parton[52] = content_h;

  // content_disable is a vector of all disabled particles, which may contain "+"'s that are resolved here:

  //  vector<string> temp_content_disable = content_disable;
  //  vector<string> temp_content_disable_split;
  // content_disable is a vector of all disabled particles, whichout "+"'s
  //  int counter = 0;
  //    int start = 0;

  vector<string> temp_content_disable;
  for (int i_d = 0; i_d < content_disable.size(); i_d++){
    string temp_split = "";
    for (int i_s = 0; i_s < content_disable[i_d].size(); i_s++){
      if (content_disable[i_d][i_s] == ' '){continue;}
      else if (content_disable[i_d][i_s] == '+'){
	if (temp_split != ""){temp_content_disable.push_back(temp_split); temp_split = "";} // counter++; 
      }
      else {temp_split.push_back(content_disable[i_d][i_s]);}
    }
    if (temp_split != ""){temp_content_disable.push_back(temp_split);}
  }



  logger << LOG_DEBUG << "temp_content_disable.size() = " << temp_content_disable.size() << endl;
  for (int i_d = 0; i_d < temp_content_disable.size(); i_d++){
    logger << LOG_DEBUG << "temp_content_disable[" << i_d << "] = " << temp_content_disable[i_d] << endl;
  }

  for (int i_d = 0; i_d < temp_content_disable.size(); i_d++){
    if (temp_content_disable[i_d] == "Q" ||
	temp_content_disable[i_d] == "q" ||
	temp_content_disable[i_d] == "qx" ||
	temp_content_disable[i_d] == "g" ||
	temp_content_disable[i_d] == "a" ||
	temp_content_disable[i_d] == "D" ||
	temp_content_disable[i_d] == "U" ||
	temp_content_disable[i_d] == "S" ||
	temp_content_disable[i_d] == "C" ||
	temp_content_disable[i_d] == "B" ||
	temp_content_disable[i_d] == "T" ||
	temp_content_disable[i_d] == "d" ||
	temp_content_disable[i_d] == "u" ||
	temp_content_disable[i_d] == "s" ||
	temp_content_disable[i_d] == "c" ||
	temp_content_disable[i_d] == "b" ||
	temp_content_disable[i_d] == "t" ||
	temp_content_disable[i_d] == "dx" ||
	temp_content_disable[i_d] == "ux" ||
	temp_content_disable[i_d] == "sx" ||
	temp_content_disable[i_d] == "cx" ||
	temp_content_disable[i_d] == "bx" ||
	temp_content_disable[i_d] == "tx"){
      for (int i_c = 0; i_c < content_parton.size(); i_c++){
	for (int i_p = 0; i_p < content_parton[i_c].size(); i_p++){
	  if ((content_parton[i_c][i_p] == 0 && temp_content_disable[i_d] == "g") ||
	      (content_parton[i_c][i_p] == 22 && temp_content_disable[i_d] == "a") ||
	      (content_parton[i_c][i_p] == 1 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "q" || temp_content_disable[i_d] == "d" || temp_content_disable[i_d] == "D")) ||
	      (content_parton[i_c][i_p] == 2 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "q" || temp_content_disable[i_d] == "u" || temp_content_disable[i_d] == "U")) ||
	      (content_parton[i_c][i_p] == 3 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "q" || temp_content_disable[i_d] == "s" || temp_content_disable[i_d] == "S")) ||
	      (content_parton[i_c][i_p] == 4 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "q" || temp_content_disable[i_d] == "c" || temp_content_disable[i_d] == "C")) ||
	      (content_parton[i_c][i_p] == 5 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "q" || temp_content_disable[i_d] == "b" || temp_content_disable[i_d] == "B")) ||
	      (content_parton[i_c][i_p] == 6 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "q" || temp_content_disable[i_d] == "t" || temp_content_disable[i_d] == "T")) ||
	      (content_parton[i_c][i_p] == -1 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "qx" || temp_content_disable[i_d] == "dx" || temp_content_disable[i_d] == "D")) ||
	      (content_parton[i_c][i_p] == -2 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "qx" || temp_content_disable[i_d] == "ux" || temp_content_disable[i_d] == "U")) ||
	      (content_parton[i_c][i_p] == -3 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "qx" || temp_content_disable[i_d] == "sx" || temp_content_disable[i_d] == "S")) ||
	      (content_parton[i_c][i_p] == -4 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "qx" || temp_content_disable[i_d] == "cx" || temp_content_disable[i_d] == "C")) ||
	      (content_parton[i_c][i_p] == -5 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "qx" || temp_content_disable[i_d] == "bx" || temp_content_disable[i_d] == "B")) ||
	      (content_parton[i_c][i_p] == -6 && (temp_content_disable[i_d] == "Q" || temp_content_disable[i_d] == "qx" || temp_content_disable[i_d] == "tx" || temp_content_disable[i_d] == "T"))){
	    content_parton[i_c].erase(content_parton[i_c].begin() + i_p);
	    i_p--;
	  }
	}
      }
    }
    else if (temp_content_disable[i_d] == "L" ||
	     temp_content_disable[i_d] == "l" ||
	     temp_content_disable[i_d] == "lx" ||
	     temp_content_disable[i_d] == "E" ||
	     temp_content_disable[i_d] == "M" ||
	     temp_content_disable[i_d] == "Y" ||
	     temp_content_disable[i_d] == "e" ||
	     temp_content_disable[i_d] == "m" ||
	     temp_content_disable[i_d] == "y" ||
	     temp_content_disable[i_d] == "ex" ||
	     temp_content_disable[i_d] == "mx" ||
	     temp_content_disable[i_d] == "yx"){
      for (int i_c = 0; i_c < content_parton.size(); i_c++){
	for (int i_p = 0; i_p < content_parton[i_c].size(); i_p++){
	  if ((content_parton[i_c][i_p] == 11 && (temp_content_disable[i_d] == "L" || temp_content_disable[i_d] == "l" || temp_content_disable[i_d] == "E" || temp_content_disable[i_d] == "e")) ||
	      (content_parton[i_c][i_p] == 13 && (temp_content_disable[i_d] == "L" || temp_content_disable[i_d] == "l" || temp_content_disable[i_d] == "M" || temp_content_disable[i_d] == "m")) ||
	      (content_parton[i_c][i_p] == 15 && (temp_content_disable[i_d] == "L" || temp_content_disable[i_d] == "l" || temp_content_disable[i_d] == "Y" || temp_content_disable[i_d] == "y")) ||
	      (content_parton[i_c][i_p] == -11 && (temp_content_disable[i_d] == "L" || temp_content_disable[i_d] == "lx" || temp_content_disable[i_d] == "E" || temp_content_disable[i_d] == "ex")) ||
	      (content_parton[i_c][i_p] == -13 && (temp_content_disable[i_d] == "L" || temp_content_disable[i_d] == "lx" || temp_content_disable[i_d] == "M" || temp_content_disable[i_d] == "mx")) ||
	      (content_parton[i_c][i_p] == -15 && (temp_content_disable[i_d] == "L" || temp_content_disable[i_d] == "lx" || temp_content_disable[i_d] == "Y" || temp_content_disable[i_d] == "yx"))){
	    content_parton[i_c].erase(content_parton[i_c].begin() + i_p);
	    i_p--;
	  }
	}
      }
    }
    else if (temp_content_disable[i_d] == "N" ||
	     temp_content_disable[i_d] == "n" ||
	     temp_content_disable[i_d] == "nx" ||
	     temp_content_disable[i_d] == "NE" ||
	     temp_content_disable[i_d] == "NM" ||
	     temp_content_disable[i_d] == "NY" ||
	     temp_content_disable[i_d] == "ne" ||
	     temp_content_disable[i_d] == "nm" ||
	     temp_content_disable[i_d] == "ny" ||
	     temp_content_disable[i_d] == "nex" ||
	     temp_content_disable[i_d] == "nmx" ||
	     temp_content_disable[i_d] == "nyx"){
      for (int i_c = 0; i_c < content_parton.size(); i_c++){
	for (int i_p = 0; i_p < content_parton[i_c].size(); i_p++){
	  if ((content_parton[i_c][i_p] == 11 && (temp_content_disable[i_d] == "N" || temp_content_disable[i_d] == "n" || temp_content_disable[i_d] == "NE" || temp_content_disable[i_d] == "ne")) ||
	      (content_parton[i_c][i_p] == 13 && (temp_content_disable[i_d] == "N" || temp_content_disable[i_d] == "n" || temp_content_disable[i_d] == "NM" || temp_content_disable[i_d] == "nm")) ||
	      (content_parton[i_c][i_p] == 15 && (temp_content_disable[i_d] == "N" || temp_content_disable[i_d] == "n" || temp_content_disable[i_d] == "NY" || temp_content_disable[i_d] == "ny")) ||
	      (content_parton[i_c][i_p] == -11 && (temp_content_disable[i_d] == "N" || temp_content_disable[i_d] == "nx" || temp_content_disable[i_d] == "NE" || temp_content_disable[i_d] == "nex")) ||
	      (content_parton[i_c][i_p] == -13 && (temp_content_disable[i_d] == "N" || temp_content_disable[i_d] == "nx" || temp_content_disable[i_d] == "NM" || temp_content_disable[i_d] == "nmx")) ||
	      (content_parton[i_c][i_p] == -15 && (temp_content_disable[i_d] == "N" || temp_content_disable[i_d] == "nx" || temp_content_disable[i_d] == "NY" || temp_content_disable[i_d] == "nyx"))){
	    content_parton[i_c].erase(content_parton[i_c].begin() + i_p);
	    i_p--;
	  }
	}
      }
    }
    else if (temp_content_disable[i_d] == "g" ||
	     temp_content_disable[i_d] == "a" ||
	     temp_content_disable[i_d] == "z" ||
	     temp_content_disable[i_d] == "W" ||
	     temp_content_disable[i_d] == "w" ||
	     temp_content_disable[i_d] == "wx" ||
	     temp_content_disable[i_d] == "h"){
      for (int i_c = 0; i_c < content_parton.size(); i_c++){
	for (int i_p = 0; i_p < content_parton[i_c].size(); i_p++){
	  if ((content_parton[i_c][i_p] == 21 && temp_content_disable[i_d] == "g") ||
	      (content_parton[i_c][i_p] == 22 && temp_content_disable[i_d] == "a") ||
	      (content_parton[i_c][i_p] == 23 && temp_content_disable[i_d] == "z") ||
	      (content_parton[i_c][i_p] == 24 && (temp_content_disable[i_d] == "W" || temp_content_disable[i_d] == "w")) ||
	      (content_parton[i_c][i_p] == -24 && (temp_content_disable[i_d] == "W" || temp_content_disable[i_d] == "wx")) ||
	      (content_parton[i_c][i_p] == 25 && temp_content_disable[i_d] == "h")){
	    content_parton[i_c].erase(content_parton[i_c].begin() + i_p);
	    i_p--;
	  }
	}
      }
    }
    else {
      logger << LOG_FATAL << "Wrong input: content_disable cannot be " << temp_content_disable[i_d] << " !" << endl;
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

  relation_content_parton["L"] = content_parton[24];
  relation_content_parton["l"] = content_parton[25];
  relation_content_parton["lx"] = content_parton[26];
  relation_content_parton["E"] = content_parton[27];
  relation_content_parton["M"] = content_parton[28];
  relation_content_parton["Y"] = content_parton[29];
  relation_content_parton["e"] = content_parton[30];
  relation_content_parton["m"] = content_parton[31];
  relation_content_parton["y"] = content_parton[32];
  relation_content_parton["ex"] = content_parton[33];
  relation_content_parton["mx"] = content_parton[34];
  relation_content_parton["yx"] = content_parton[35];

  relation_content_parton["N"] = content_parton[36];
  relation_content_parton["n"] = content_parton[37];
  relation_content_parton["nx"] = content_parton[38];
  relation_content_parton["NE"] = content_parton[39];
  relation_content_parton["NM"] = content_parton[40];
  relation_content_parton["NY"] = content_parton[41];
  relation_content_parton["ne"] = content_parton[42];
  relation_content_parton["nm"] = content_parton[43];
  relation_content_parton["ny"] = content_parton[44];
  relation_content_parton["nex"] = content_parton[45];
  relation_content_parton["nmx"] = content_parton[46];
  relation_content_parton["nyx"] = content_parton[47];

  relation_content_parton["z"] = content_parton[48];

  relation_content_parton["W"] = content_parton[49];
  relation_content_parton["w"] = content_parton[50];
  relation_content_parton["wx"] = content_parton[51];

  relation_content_parton["h"] = content_parton[52];




  //  content_selection
  /*
  vector<string> temp_content_selection(1);
  int start = 0;
  int counter = 0;
  for (int i_s = 0; i_s < content_selection.size(); i_s++){
    if (content_selection[i_s] == ' '){continue;}
    else if (content_selection[i_s] == '+'){counter++; temp_content_selection.push_back("");}
    else {temp_content_selection[counter].push_back(content_selection[i_s]);}
  }
  */
  vector<string> temp_content_selection;
  for (int i_d = 0; i_d < content_selection.size(); i_d++){
    string temp_split = "";
    for (int i_s = 0; i_s < content_selection[i_d].size(); i_s++){
      if (content_selection[i_d][i_s] == ' '){continue;}
      else if (content_selection[i_d][i_s] == '+'){
	if (temp_split != ""){temp_content_selection.push_back(temp_split); temp_split = "";} // counter++; 
      }
      else {temp_split.push_back(content_selection[i_d][i_s]);}
    }
    if (temp_split != ""){temp_content_selection.push_back(temp_split);}
  }

  ///  if (temp_content_selection.size() == 0){temp_content_selection.push_back("pp");}

  for (int i_d = 0; i_d < temp_content_selection.size(); i_d++){
    logger << LOG_DEBUG << "temp_content_selection[" << i_d << "] = " << temp_content_selection[i_d] << endl;
  }

  vector<string> temp_content_list;
  for (int i_s = 0; i_s < temp_content_selection.size(); i_s++){
    logger << LOG_DEBUG << "temp_content_selection[" << i_s << "] = " << temp_content_selection[i_s] << endl;
    int flag = -1;
    for (int i_d = 0; i_d < temp_content_disable.size(); i_d++){
      if (temp_content_selection[i_s] == temp_content_disable[i_d]){flag = i_d; break;}
    }
    if (flag == -1){
      for (int i_l = 0; i_l < temp_content_list.size(); i_l++){
	if (temp_content_selection[i_s] == temp_content_list[i_l]){flag = i_l; break;}
      }
      if (flag == -1){
	temp_content_list.push_back(temp_content_selection[i_s]);
      }
    }
  }

  for (int i_l = 0; i_l < temp_content_list.size(); i_l++){
    logger << LOG_DEBUG << "temp_content_list[" << i_l << "] = " << temp_content_list[i_l] << endl;
  }
  for (int i_l = 0; i_l < temp_content_list.size(); i_l++){
    for (int j_l = 0; j_l < relation_content_parton[temp_content_list[i_l]].size(); j_l++){
      int flag = -1;
      for (int i_n = 0; i_n < content_list.size(); i_n++){
	if (relation_content_parton[temp_content_list[i_l]][j_l] == content_list[i_n]){flag = i_n; break;}
      }
      if (flag == -1){content_list.push_back(relation_content_parton[temp_content_list[i_l]][j_l]);}
    }
  }
  
  for (int i_l = 0; i_l < content_list.size(); i_l++){
    if (content_list[i_l] == 21){content_list[i_l] = 0;}
    //    if (content_list[i_l] == 22){content_list[i_l] = 7;}
  }

  sort(content_list.begin(), content_list.end());
  for (int i_l = 0; i_l < content_list.size(); i_l++){
    logger << LOG_DEBUG << "content_list[" << i_l << "] = " << content_list[i_l] << endl;
  }

  //  content_selection
  //  content_disable

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

