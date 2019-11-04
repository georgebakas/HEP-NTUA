#include "ppttx20.header.cxx"
void ppttx20_determination_no_subprocess_born(int & no_map, vector<int> & o_map, int & no_prc, vector<int> & o_prc, double & symmetry_factor, vector<int> & tp, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference){
  static Logger logger("ppttx20_determination_no_subprocess_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  no_map = -1;
  vector<int> out(5);
  for (int i = 0; i < 2; i++){
    if (i == 0){out[1] = 1; out[2] = 2;}
    if (i == 1){out[1] = 2; out[2] = 1;}
    for (int o = 0; o < 1; o++){
      if (o == 0){}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6)){no_map = 1; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6)){no_map = 2; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6)){no_map = 2; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6)){no_map = 3; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6)){no_map = 3; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6)){no_map = 4; break;}
    }
    if (no_map != -1){break;}
  }
  for (int o = 0; o < out.size(); o++){
    if (out[o] == 0){o_map[o] = o;}
    else {o_map[o] = out[o];}
  }
  no_prc = no_map;
  o_prc = o_map;
  if (no_map == 0){}
  else if (no_map == 1){symmetry_factor = 1;}
  else if (no_map == 2){symmetry_factor = 1;}
  else if (no_map == 3){symmetry_factor = 1;}
  else if (no_map == 4){symmetry_factor = 1;}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
void ppttx20_determination_subprocess_born(int i_a, phasespace_set & psi){
  static Logger logger("ppttx20_determination_subprocess_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<int> tp = psi.csi->basic_type_parton[i_a];
  //  ppttx20_determination_no_subprocess_born(psi_no_map[i_a], psi_o_map[i_a], psi_no_prc[i_a], psi_o_prc[i_a], psi_symmetry_factor, tp, psi_phasespace_order_alpha_s[i_a], psi_phasespace_order_alpha_e[i_a], psi_phasespace_order_interference[i_a]);

  psi_no_map[i_a] = -1;
  vector<int> out(5);
  for (int i = 0; i < 2; i++){
    if (i == 0){out[1] = 1; out[2] = 2;}
    if (i == 1){out[1] = 2; out[2] = 1;}
    for (int o = 0; o < 1; o++){
      if (o == 0){}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6)){psi_no_map[i_a] = 1; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6)){psi_no_map[i_a] = 2; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6)){psi_no_map[i_a] = 2; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6)){psi_no_map[i_a] = 3; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6)){psi_no_map[i_a] = 3; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6)){psi_no_map[i_a] = 4; break;}
    }
    if (psi_no_map[i_a] != -1){break;}
  }
  for (int o = 0; o < out.size(); o++){
    if (out[o] == 0){psi_o_map[i_a][o] = o;}
    else {psi_o_map[i_a][o] = out[o];}
  }
  psi_no_prc[i_a] = psi_no_map[i_a];
  psi_o_prc[i_a] = psi_o_map[i_a];
  if (psi_no_map[i_a] == 0){}
  else if (psi_no_map[i_a] == 1){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 2){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 3){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 4){psi_symmetry_factor = 1;}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
void ppttx20_combination_subprocess_born(int i_a, phasespace_set & psi, observable_set & oset){
  static Logger logger("ppttx20_combination_subprocess_born");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (psi_no_map[i_a] == 0){}
  else if (psi_no_map[i_a] == 1){
    int pdf_1[] = { 1,   0,   0};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  g   -> t  tx 
  }
  else if (psi_no_map[i_a] == 2){
    int pdf_1[] = { 1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  dx  -> t  tx 
    int pdf_2[] = { 1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  sx  -> t  tx 
    int pdf_3[] = {-1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx d   -> t  tx 
    int pdf_4[] = {-1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx s   -> t  tx 
  }
  else if (psi_no_map[i_a] == 3){
    int pdf_1[] = { 1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  ux  -> t  tx 
    int pdf_2[] = { 1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  cx  -> t  tx 
    int pdf_3[] = {-1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux u   -> t  tx 
    int pdf_4[] = {-1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx c   -> t  tx 
  }
  else if (psi_no_map[i_a] == 4){
    int pdf_1[] = { 1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // b  bx  -> t  tx 
    int pdf_2[] = {-1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // bx b   -> t  tx 
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
