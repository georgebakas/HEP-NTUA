#include "ppttx20.header.cxx"
void ppttx20_determination_no_subprocess_real(int & no_map, vector<int> & o_map, int & no_prc, vector<int> & o_prc, double & symmetry_factor, vector<int> & tp, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference){
  static Logger logger("ppttx20_determination_no_subprocess_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  no_map = -1;
  vector<int> out(6);
  for (int i = 0; i < 2; i++){
    if (i == 0){out[1] = 1; out[2] = 2;}
    if (i == 1){out[1] = 2; out[2] = 1;}
    for (int o = 0; o < 1; o++){
      if (o == 0){out[5] = 5;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){no_map = 1; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1)){no_map = 2; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3)){no_map = 2; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2)){no_map = 3; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4)){no_map = 3; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5)){no_map = 4; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1)){no_map = 5; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3)){no_map = 5; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -2)){no_map = 6; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -4)){no_map = 6; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -5)){no_map = 7; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){no_map = 8; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){no_map = 8; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){no_map = 9; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){no_map = 9; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){no_map = 10; break;}
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
  else if (no_map == 5){symmetry_factor = 1;}
  else if (no_map == 6){symmetry_factor = 1;}
  else if (no_map == 7){symmetry_factor = 1;}
  else if (no_map == 8){symmetry_factor = 1;}
  else if (no_map == 9){symmetry_factor = 1;}
  else if (no_map == 10){symmetry_factor = 1;}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
void ppttx20_determination_subprocess_real(int i_a, phasespace_set & psi){
  static Logger logger("ppttx20_determination_subprocess_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<int> tp = psi.csi->basic_type_parton[i_a];
  //  ppttx20_determination_no_subprocess_real(psi_no_map[i_a], psi_o_map[i_a], psi_no_prc[i_a], psi_o_prc[i_a], psi_symmetry_factor, tp, psi_phasespace_order_alpha_s[i_a], psi_phasespace_order_alpha_e[i_a], psi_phasespace_order_interference[i_a]);

  psi_no_map[i_a] = -1;
  vector<int> out(6);
  for (int i = 0; i < 2; i++){
    if (i == 0){out[1] = 1; out[2] = 2;}
    if (i == 1){out[1] = 2; out[2] = 1;}
    for (int o = 0; o < 1; o++){
      if (o == 0){out[5] = 5;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){psi_no_map[i_a] = 1; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1)){psi_no_map[i_a] = 2; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3)){psi_no_map[i_a] = 2; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2)){psi_no_map[i_a] = 3; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4)){psi_no_map[i_a] = 3; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5)){psi_no_map[i_a] = 4; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1)){psi_no_map[i_a] = 5; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3)){psi_no_map[i_a] = 5; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -2)){psi_no_map[i_a] = 6; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -4)){psi_no_map[i_a] = 6; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -5)){psi_no_map[i_a] = 7; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){psi_no_map[i_a] = 8; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){psi_no_map[i_a] = 8; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){psi_no_map[i_a] = 9; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){psi_no_map[i_a] = 9; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0)){psi_no_map[i_a] = 10; break;}
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
  else if (psi_no_map[i_a] == 5){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 6){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 7){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 8){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 9){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 10){psi_symmetry_factor = 1;}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
void ppttx20_combination_subprocess_real(int i_a, phasespace_set & psi, observable_set & oset){
  static Logger logger("ppttx20_combination_subprocess_real");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (psi_no_map[i_a] == 0){}
  else if (psi_no_map[i_a] == 1){
    int pdf_1[] = { 1,   0,   0};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  g   -> t  tx g  
  }
  else if (psi_no_map[i_a] == 2){
    int pdf_1[] = { 1,   0,   1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  d   -> t  tx d  
    int pdf_2[] = { 1,   0,   3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // g  s   -> t  tx s  
    int pdf_3[] = {-1,   0,   1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // d  g   -> t  tx d  
    int pdf_4[] = {-1,   0,   3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // s  g   -> t  tx s  
  }
  else if (psi_no_map[i_a] == 3){
    int pdf_1[] = { 1,   0,   2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  u   -> t  tx u  
    int pdf_2[] = { 1,   0,   4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // g  c   -> t  tx c  
    int pdf_3[] = {-1,   0,   2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // u  g   -> t  tx u  
    int pdf_4[] = {-1,   0,   4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // c  g   -> t  tx c  
  }
  else if (psi_no_map[i_a] == 4){
    int pdf_1[] = { 1,   0,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  b   -> t  tx b  
    int pdf_2[] = {-1,   0,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // b  g   -> t  tx b  
  }
  else if (psi_no_map[i_a] == 5){
    int pdf_1[] = { 1,   0,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  dx  -> t  tx dx 
    int pdf_2[] = { 1,   0,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // g  sx  -> t  tx sx 
    int pdf_3[] = {-1,   0,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx g   -> t  tx dx 
    int pdf_4[] = {-1,   0,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx g   -> t  tx sx 
  }
  else if (psi_no_map[i_a] == 6){
    int pdf_1[] = { 1,   0,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  ux  -> t  tx ux 
    int pdf_2[] = { 1,   0,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // g  cx  -> t  tx cx 
    int pdf_3[] = {-1,   0,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux g   -> t  tx ux 
    int pdf_4[] = {-1,   0,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx g   -> t  tx cx 
  }
  else if (psi_no_map[i_a] == 7){
    int pdf_1[] = { 1,   0,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  bx  -> t  tx bx 
    int pdf_2[] = {-1,   0,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // bx g   -> t  tx bx 
  }
  else if (psi_no_map[i_a] == 8){
    int pdf_1[] = { 1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  dx  -> t  tx g  
    int pdf_2[] = { 1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  sx  -> t  tx g  
    int pdf_3[] = {-1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx d   -> t  tx g  
    int pdf_4[] = {-1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx s   -> t  tx g  
  }
  else if (psi_no_map[i_a] == 9){
    int pdf_1[] = { 1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  ux  -> t  tx g  
    int pdf_2[] = { 1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  cx  -> t  tx g  
    int pdf_3[] = {-1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux u   -> t  tx g  
    int pdf_4[] = {-1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx c   -> t  tx g  
  }
  else if (psi_no_map[i_a] == 10){
    int pdf_1[] = { 1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // b  bx  -> t  tx g  
    int pdf_2[] = {-1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // bx b   -> t  tx g  
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
