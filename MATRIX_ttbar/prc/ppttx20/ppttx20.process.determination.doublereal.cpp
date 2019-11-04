#include "ppttx20.header.cxx"
void ppttx20_determination_no_subprocess_doublereal(int & no_map, vector<int> & o_map, int & no_prc, vector<int> & o_prc, double & symmetry_factor, vector<int> & tp, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference){
  static Logger logger("ppttx20_determination_no_subprocess_doublereal");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  no_map = -1;
  vector<int> out(7);
  for (int i = 0; i < 2; i++){
    if (i == 0){out[1] = 1; out[2] = 2;}
    if (i == 1){out[1] = 2; out[2] = 1;}
    for (int o = 0; o < 2; o++){
      if (o == 0){out[5] = 5; out[6] = 6;}
      if (o == 1){out[5] = 6; out[6] = 5;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){no_map = 1; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){no_map = 2; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){no_map = 2; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){no_map = 3; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){no_map = 3; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){no_map = 4; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   1)){no_map = 5; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   3)){no_map = 5; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   2)){no_map = 6; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   4)){no_map = 6; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   5)){no_map = 7; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==  -1)){no_map = 8; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==  -3)){no_map = 8; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==  -2)){no_map = 9; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==  -4)){no_map = 9; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==  -5)){no_map = 10; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==   1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==   1)){no_map = 11; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==   3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==   3)){no_map = 11; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==   2)){no_map = 12; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==   4)){no_map = 12; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==   3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==   3)){no_map = 13; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==   1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==   1)){no_map = 13; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==   4)){no_map = 14; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==   2)){no_map = 14; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==   5)){no_map = 16; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==   5)){no_map = 16; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){no_map = 17; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){no_map = 17; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){no_map = 18; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){no_map = 18; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){no_map = 19; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){no_map = 19; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){no_map = 20; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){no_map = 20; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){no_map = 21; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){no_map = 21; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){no_map = 22; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){no_map = 22; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -2)){no_map = 23; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -4)){no_map = 23; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -3)){no_map = 25; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -1)){no_map = 25; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -4)){no_map = 27; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -2)){no_map = 27; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -5)){no_map = 28; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -5)){no_map = 28; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==   2)){no_map = 29; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==   4)){no_map = 29; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==   4)){no_map = 30; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==   2)){no_map = 30; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==   5)){no_map = 31; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==   5)){no_map = 31; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -1)){no_map = 32; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -3)){no_map = 32; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){no_map = 34; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){no_map = 34; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){no_map = 35; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){no_map = 35; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){no_map = 36; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){no_map = 36; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){no_map = 37; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){no_map = 37; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){no_map = 38; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){no_map = 38; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){no_map = 39; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){no_map = 39; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -3)){no_map = 40; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -1)){no_map = 40; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -4)){no_map = 42; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -2)){no_map = 42; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -5)){no_map = 43; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -5)){no_map = 43; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==   5)){no_map = 44; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -1)){no_map = 45; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -3)){no_map = 45; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -2)){no_map = 46; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -4)){no_map = 46; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){no_map = 47; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){no_map = 48; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){no_map = 48; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){no_map = 49; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){no_map = 49; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){no_map = 50; break;}
      if ((tp[out[1]] ==  -1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1) && (tp[out[6]] ==  -1)){no_map = 51; break;}
      if ((tp[out[1]] ==  -3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3) && (tp[out[6]] ==  -3)){no_map = 51; break;}
      if ((tp[out[1]] ==  -1) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1) && (tp[out[6]] ==  -2)){no_map = 52; break;}
      if ((tp[out[1]] ==  -3) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3) && (tp[out[6]] ==  -4)){no_map = 52; break;}
      if ((tp[out[1]] ==  -1) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1) && (tp[out[6]] ==  -3)){no_map = 53; break;}
      if ((tp[out[1]] ==  -3) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3) && (tp[out[6]] ==  -1)){no_map = 53; break;}
      if ((tp[out[1]] ==  -1) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1) && (tp[out[6]] ==  -4)){no_map = 54; break;}
      if ((tp[out[1]] ==  -3) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3) && (tp[out[6]] ==  -2)){no_map = 54; break;}
      if ((tp[out[1]] ==  -1) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1) && (tp[out[6]] ==  -5)){no_map = 56; break;}
      if ((tp[out[1]] ==  -3) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3) && (tp[out[6]] ==  -5)){no_map = 56; break;}
      if ((tp[out[1]] ==  -2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -2) && (tp[out[6]] ==  -2)){no_map = 57; break;}
      if ((tp[out[1]] ==  -4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -4) && (tp[out[6]] ==  -4)){no_map = 57; break;}
      if ((tp[out[1]] ==  -2) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -2) && (tp[out[6]] ==  -4)){no_map = 58; break;}
      if ((tp[out[1]] ==  -4) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -4) && (tp[out[6]] ==  -2)){no_map = 58; break;}
      if ((tp[out[1]] ==  -2) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -2) && (tp[out[6]] ==  -5)){no_map = 59; break;}
      if ((tp[out[1]] ==  -4) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -4) && (tp[out[6]] ==  -5)){no_map = 59; break;}
      if ((tp[out[1]] ==  -5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -5) && (tp[out[6]] ==  -5)){no_map = 60; break;}
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
  else if (no_map == 1){symmetry_factor = 2;}
  else if (no_map == 2){symmetry_factor = 1;}
  else if (no_map == 3){symmetry_factor = 1;}
  else if (no_map == 4){symmetry_factor = 1;}
  else if (no_map == 5){symmetry_factor = 1;}
  else if (no_map == 6){symmetry_factor = 1;}
  else if (no_map == 7){symmetry_factor = 1;}
  else if (no_map == 8){symmetry_factor = 1;}
  else if (no_map == 9){symmetry_factor = 1;}
  else if (no_map == 10){symmetry_factor = 1;}
  else if (no_map == 11){symmetry_factor = 2;}
  else if (no_map == 12){symmetry_factor = 1;}
  else if (no_map == 13){symmetry_factor = 1;}
  else if (no_map == 14){symmetry_factor = 1;}
  else if (no_map == 16){symmetry_factor = 1;}
  else if (no_map == 17){symmetry_factor = 2;}
  else if (no_map == 18){symmetry_factor = 1;}
  else if (no_map == 19){symmetry_factor = 1;}
  else if (no_map == 20){symmetry_factor = 1;}
  else if (no_map == 21){symmetry_factor = 1;}
  else if (no_map == 22){symmetry_factor = 1;}
  else if (no_map == 23){symmetry_factor = 1;}
  else if (no_map == 25){symmetry_factor = 1;}
  else if (no_map == 27){symmetry_factor = 1;}
  else if (no_map == 28){symmetry_factor = 1;}
  else if (no_map == 29){symmetry_factor = 2;}
  else if (no_map == 30){symmetry_factor = 1;}
  else if (no_map == 31){symmetry_factor = 1;}
  else if (no_map == 32){symmetry_factor = 1;}
  else if (no_map == 34){symmetry_factor = 2;}
  else if (no_map == 35){symmetry_factor = 1;}
  else if (no_map == 36){symmetry_factor = 1;}
  else if (no_map == 37){symmetry_factor = 1;}
  else if (no_map == 38){symmetry_factor = 1;}
  else if (no_map == 39){symmetry_factor = 1;}
  else if (no_map == 40){symmetry_factor = 1;}
  else if (no_map == 42){symmetry_factor = 1;}
  else if (no_map == 43){symmetry_factor = 1;}
  else if (no_map == 44){symmetry_factor = 2;}
  else if (no_map == 45){symmetry_factor = 1;}
  else if (no_map == 46){symmetry_factor = 1;}
  else if (no_map == 47){symmetry_factor = 2;}
  else if (no_map == 48){symmetry_factor = 1;}
  else if (no_map == 49){symmetry_factor = 1;}
  else if (no_map == 50){symmetry_factor = 1;}
  else if (no_map == 51){symmetry_factor = 2;}
  else if (no_map == 52){symmetry_factor = 1;}
  else if (no_map == 53){symmetry_factor = 1;}
  else if (no_map == 54){symmetry_factor = 1;}
  else if (no_map == 56){symmetry_factor = 1;}
  else if (no_map == 57){symmetry_factor = 2;}
  else if (no_map == 58){symmetry_factor = 1;}
  else if (no_map == 59){symmetry_factor = 1;}
  else if (no_map == 60){symmetry_factor = 2;}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
void ppttx20_determination_subprocess_doublereal(int i_a, phasespace_set & psi){
  static Logger logger("ppttx20_determination_subprocess_doublereal");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  vector<int> tp = psi.csi->basic_type_parton[i_a];
  //  ppttx20_determination_no_subprocess_doublereal(psi_no_map[i_a], psi_o_map[i_a], psi_no_prc[i_a], psi_o_prc[i_a], psi_symmetry_factor, tp, psi_phasespace_order_alpha_s[i_a], psi_phasespace_order_alpha_e[i_a], psi_phasespace_order_interference[i_a]);

  psi_no_map[i_a] = -1;
  vector<int> out(7);
  for (int i = 0; i < 2; i++){
    if (i == 0){out[1] = 1; out[2] = 2;}
    if (i == 1){out[1] = 2; out[2] = 1;}
    for (int o = 0; o < 2; o++){
      if (o == 0){out[5] = 5; out[6] = 6;}
      if (o == 1){out[5] = 6; out[6] = 5;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){psi_no_map[i_a] = 1; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 2; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 2; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 3; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 3; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   0) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 4; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   1)){psi_no_map[i_a] = 5; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   3)){psi_no_map[i_a] = 5; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   2)){psi_no_map[i_a] = 6; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   4)){psi_no_map[i_a] = 6; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   5)){psi_no_map[i_a] = 7; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 8; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 8; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 9; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 9; break;}
      if ((tp[out[1]] ==   0) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 10; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==   1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==   1)){psi_no_map[i_a] = 11; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==   3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==   3)){psi_no_map[i_a] = 11; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==   2)){psi_no_map[i_a] = 12; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==   4)){psi_no_map[i_a] = 12; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==   3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==   3)){psi_no_map[i_a] = 13; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==   1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==   1)){psi_no_map[i_a] = 13; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==   4)){psi_no_map[i_a] = 14; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==   2)){psi_no_map[i_a] = 14; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==   5)){psi_no_map[i_a] = 16; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==   5)){psi_no_map[i_a] = 16; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){psi_no_map[i_a] = 17; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){psi_no_map[i_a] = 17; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 18; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 18; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 19; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 19; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 20; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 20; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 21; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 21; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 22; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 22; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 23; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 23; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 25; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 25; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 27; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 27; break;}
      if ((tp[out[1]] ==   1) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 28; break;}
      if ((tp[out[1]] ==   3) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 28; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==   2)){psi_no_map[i_a] = 29; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==   4)){psi_no_map[i_a] = 29; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==   4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==   4)){psi_no_map[i_a] = 30; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==   2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==   2)){psi_no_map[i_a] = 30; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==   5)){psi_no_map[i_a] = 31; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==   5)){psi_no_map[i_a] = 31; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 32; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 32; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){psi_no_map[i_a] = 34; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){psi_no_map[i_a] = 34; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 35; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 35; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 36; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 36; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 37; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 37; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 38; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 38; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 39; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 39; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 40; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 40; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 42; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 42; break;}
      if ((tp[out[1]] ==   2) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 43; break;}
      if ((tp[out[1]] ==   4) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 43; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==   5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==   5)){psi_no_map[i_a] = 44; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 45; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 45; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 46; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 46; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   0) && (tp[out[6]] ==   0)){psi_no_map[i_a] = 47; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   1) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 48; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   3) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 48; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   2) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 49; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   4) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 49; break;}
      if ((tp[out[1]] ==   5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==   5) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 50; break;}
      if ((tp[out[1]] ==  -1) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 51; break;}
      if ((tp[out[1]] ==  -3) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 51; break;}
      if ((tp[out[1]] ==  -1) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 52; break;}
      if ((tp[out[1]] ==  -3) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 52; break;}
      if ((tp[out[1]] ==  -1) && (tp[out[2]] ==  -3) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1) && (tp[out[6]] ==  -3)){psi_no_map[i_a] = 53; break;}
      if ((tp[out[1]] ==  -3) && (tp[out[2]] ==  -1) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3) && (tp[out[6]] ==  -1)){psi_no_map[i_a] = 53; break;}
      if ((tp[out[1]] ==  -1) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 54; break;}
      if ((tp[out[1]] ==  -3) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 54; break;}
      if ((tp[out[1]] ==  -1) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -1) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 56; break;}
      if ((tp[out[1]] ==  -3) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -3) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 56; break;}
      if ((tp[out[1]] ==  -2) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -2) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 57; break;}
      if ((tp[out[1]] ==  -4) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -4) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 57; break;}
      if ((tp[out[1]] ==  -2) && (tp[out[2]] ==  -4) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -2) && (tp[out[6]] ==  -4)){psi_no_map[i_a] = 58; break;}
      if ((tp[out[1]] ==  -4) && (tp[out[2]] ==  -2) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -4) && (tp[out[6]] ==  -2)){psi_no_map[i_a] = 58; break;}
      if ((tp[out[1]] ==  -2) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -2) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 59; break;}
      if ((tp[out[1]] ==  -4) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -4) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 59; break;}
      if ((tp[out[1]] ==  -5) && (tp[out[2]] ==  -5) && (tp[3] ==   6) && (tp[4] ==  -6) && (tp[out[5]] ==  -5) && (tp[out[6]] ==  -5)){psi_no_map[i_a] = 60; break;}
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
  else if (psi_no_map[i_a] == 1){psi_symmetry_factor = 2;}
  else if (psi_no_map[i_a] == 2){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 3){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 4){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 5){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 6){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 7){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 8){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 9){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 10){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 11){psi_symmetry_factor = 2;}
  else if (psi_no_map[i_a] == 12){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 13){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 14){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 16){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 17){psi_symmetry_factor = 2;}
  else if (psi_no_map[i_a] == 18){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 19){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 20){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 21){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 22){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 23){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 25){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 27){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 28){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 29){psi_symmetry_factor = 2;}
  else if (psi_no_map[i_a] == 30){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 31){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 32){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 34){psi_symmetry_factor = 2;}
  else if (psi_no_map[i_a] == 35){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 36){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 37){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 38){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 39){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 40){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 42){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 43){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 44){psi_symmetry_factor = 2;}
  else if (psi_no_map[i_a] == 45){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 46){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 47){psi_symmetry_factor = 2;}
  else if (psi_no_map[i_a] == 48){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 49){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 50){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 51){psi_symmetry_factor = 2;}
  else if (psi_no_map[i_a] == 52){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 53){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 54){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 56){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 57){psi_symmetry_factor = 2;}
  else if (psi_no_map[i_a] == 58){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 59){psi_symmetry_factor = 1;}
  else if (psi_no_map[i_a] == 60){psi_symmetry_factor = 2;}

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
void ppttx20_combination_subprocess_doublereal(int i_a, phasespace_set & psi, observable_set & oset){
  static Logger logger("ppttx20_combination_subprocess_doublereal");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (psi_no_map[i_a] == 0){}
  else if (psi_no_map[i_a] == 1){
    int pdf_1[] = { 1,   0,   0};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  g   -> t  tx g  g  
  }
  else if (psi_no_map[i_a] == 2){
    int pdf_1[] = { 1,   0,   0};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  g   -> t  tx d  dx 
    int pdf_2[] = { 1,   0,   0};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // g  g   -> t  tx s  sx 
  }
  else if (psi_no_map[i_a] == 3){
    int pdf_1[] = { 1,   0,   0};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  g   -> t  tx u  ux 
    int pdf_2[] = { 1,   0,   0};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // g  g   -> t  tx c  cx 
  }
  else if (psi_no_map[i_a] == 4){
    int pdf_1[] = { 1,   0,   0};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  g   -> t  tx b  bx 
  }
  else if (psi_no_map[i_a] == 5){
    int pdf_1[] = { 1,   0,   1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  d   -> t  tx g  d  
    int pdf_2[] = { 1,   0,   3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // g  s   -> t  tx g  s  
    int pdf_3[] = {-1,   0,   1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // d  g   -> t  tx g  d  
    int pdf_4[] = {-1,   0,   3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // s  g   -> t  tx g  s  
  }
  else if (psi_no_map[i_a] == 6){
    int pdf_1[] = { 1,   0,   2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  u   -> t  tx g  u  
    int pdf_2[] = { 1,   0,   4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // g  c   -> t  tx g  c  
    int pdf_3[] = {-1,   0,   2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // u  g   -> t  tx g  u  
    int pdf_4[] = {-1,   0,   4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // c  g   -> t  tx g  c  
  }
  else if (psi_no_map[i_a] == 7){
    int pdf_1[] = { 1,   0,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  b   -> t  tx g  b  
    int pdf_2[] = {-1,   0,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // b  g   -> t  tx g  b  
  }
  else if (psi_no_map[i_a] == 8){
    int pdf_1[] = { 1,   0,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  dx  -> t  tx g  dx 
    int pdf_2[] = { 1,   0,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // g  sx  -> t  tx g  sx 
    int pdf_3[] = {-1,   0,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx g   -> t  tx g  dx 
    int pdf_4[] = {-1,   0,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx g   -> t  tx g  sx 
  }
  else if (psi_no_map[i_a] == 9){
    int pdf_1[] = { 1,   0,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  ux  -> t  tx g  ux 
    int pdf_2[] = { 1,   0,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // g  cx  -> t  tx g  cx 
    int pdf_3[] = {-1,   0,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux g   -> t  tx g  ux 
    int pdf_4[] = {-1,   0,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx g   -> t  tx g  cx 
  }
  else if (psi_no_map[i_a] == 10){
    int pdf_1[] = { 1,   0,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // g  bx  -> t  tx g  bx 
    int pdf_2[] = {-1,   0,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // bx g   -> t  tx g  bx 
  }
  else if (psi_no_map[i_a] == 11){
    int pdf_1[] = { 1,   1,   1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  d   -> t  tx d  d  
    int pdf_2[] = { 1,   3,   3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  s   -> t  tx s  s  
  }
  else if (psi_no_map[i_a] == 12){
    int pdf_1[] = { 1,   1,   2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  u   -> t  tx d  u  
    int pdf_2[] = { 1,   3,   4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  c   -> t  tx s  c  
    int pdf_3[] = {-1,   1,   2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // u  d   -> t  tx d  u  
    int pdf_4[] = {-1,   3,   4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // c  s   -> t  tx s  c  
  }
  else if (psi_no_map[i_a] == 13){
    int pdf_1[] = { 1,   1,   3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  s   -> t  tx d  s  
    int pdf_2[] = { 1,   3,   1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  d   -> t  tx d  s  
  }
  else if (psi_no_map[i_a] == 14){
    int pdf_1[] = { 1,   1,   4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  c   -> t  tx d  c  
    int pdf_2[] = { 1,   3,   2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  u   -> t  tx u  s  
    int pdf_3[] = {-1,   3,   2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // u  s   -> t  tx u  s  
    int pdf_4[] = {-1,   1,   4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // c  d   -> t  tx d  c  
  }
  else if (psi_no_map[i_a] == 16){
    int pdf_1[] = { 1,   1,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  b   -> t  tx d  b  
    int pdf_2[] = { 1,   3,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  b   -> t  tx s  b  
    int pdf_3[] = {-1,   1,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // b  d   -> t  tx d  b  
    int pdf_4[] = {-1,   3,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // b  s   -> t  tx s  b  
  }
  else if (psi_no_map[i_a] == 17){
    int pdf_1[] = { 1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  dx  -> t  tx g  g  
    int pdf_2[] = { 1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  sx  -> t  tx g  g  
    int pdf_3[] = {-1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx d   -> t  tx g  g  
    int pdf_4[] = {-1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx s   -> t  tx g  g  
  }
  else if (psi_no_map[i_a] == 18){
    int pdf_1[] = { 1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  dx  -> t  tx d  dx 
    int pdf_2[] = { 1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  sx  -> t  tx s  sx 
    int pdf_3[] = {-1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx d   -> t  tx d  dx 
    int pdf_4[] = {-1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx s   -> t  tx s  sx 
  }
  else if (psi_no_map[i_a] == 19){
    int pdf_1[] = { 1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  dx  -> t  tx u  ux 
    int pdf_2[] = { 1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  sx  -> t  tx c  cx 
    int pdf_3[] = {-1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx d   -> t  tx u  ux 
    int pdf_4[] = {-1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx s   -> t  tx c  cx 
  }
  else if (psi_no_map[i_a] == 20){
    int pdf_1[] = { 1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  dx  -> t  tx s  sx 
    int pdf_2[] = { 1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  sx  -> t  tx d  dx 
    int pdf_3[] = {-1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx d   -> t  tx s  sx 
    int pdf_4[] = {-1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx s   -> t  tx d  dx 
  }
  else if (psi_no_map[i_a] == 21){
    int pdf_1[] = { 1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  dx  -> t  tx c  cx 
    int pdf_2[] = { 1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  sx  -> t  tx u  ux 
    int pdf_3[] = {-1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx d   -> t  tx c  cx 
    int pdf_4[] = {-1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx s   -> t  tx u  ux 
  }
  else if (psi_no_map[i_a] == 22){
    int pdf_1[] = { 1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  dx  -> t  tx b  bx 
    int pdf_2[] = { 1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  sx  -> t  tx b  bx 
    int pdf_3[] = {-1,   1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx d   -> t  tx b  bx 
    int pdf_4[] = {-1,   3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx s   -> t  tx b  bx 
  }
  else if (psi_no_map[i_a] == 23){
    int pdf_1[] = { 1,   1,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  ux  -> t  tx d  ux 
    int pdf_2[] = { 1,   3,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  cx  -> t  tx s  cx 
    int pdf_3[] = {-1,   1,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux d   -> t  tx d  ux 
    int pdf_4[] = {-1,   3,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx s   -> t  tx s  cx 
  }
  else if (psi_no_map[i_a] == 25){
    int pdf_1[] = { 1,   1,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  sx  -> t  tx d  sx 
    int pdf_2[] = { 1,   3,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  dx  -> t  tx s  dx 
    int pdf_3[] = {-1,   3,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx s   -> t  tx s  dx 
    int pdf_4[] = {-1,   1,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx d   -> t  tx d  sx 
  }
  else if (psi_no_map[i_a] == 27){
    int pdf_1[] = { 1,   1,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  cx  -> t  tx d  cx 
    int pdf_2[] = { 1,   3,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  ux  -> t  tx s  ux 
    int pdf_3[] = {-1,   3,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux s   -> t  tx s  ux 
    int pdf_4[] = {-1,   1,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx d   -> t  tx d  cx 
  }
  else if (psi_no_map[i_a] == 28){
    int pdf_1[] = { 1,   1,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // d  bx  -> t  tx d  bx 
    int pdf_2[] = { 1,   3,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // s  bx  -> t  tx s  bx 
    int pdf_3[] = {-1,   1,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // bx d   -> t  tx d  bx 
    int pdf_4[] = {-1,   3,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // bx s   -> t  tx s  bx 
  }
  else if (psi_no_map[i_a] == 29){
    int pdf_1[] = { 1,   2,   2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  u   -> t  tx u  u  
    int pdf_2[] = { 1,   4,   4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  c   -> t  tx c  c  
  }
  else if (psi_no_map[i_a] == 30){
    int pdf_1[] = { 1,   2,   4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  c   -> t  tx u  c  
    int pdf_2[] = { 1,   4,   2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  u   -> t  tx u  c  
  }
  else if (psi_no_map[i_a] == 31){
    int pdf_1[] = { 1,   2,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  b   -> t  tx u  b  
    int pdf_2[] = { 1,   4,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  b   -> t  tx c  b  
    int pdf_3[] = {-1,   2,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // b  u   -> t  tx u  b  
    int pdf_4[] = {-1,   4,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // b  c   -> t  tx c  b  
  }
  else if (psi_no_map[i_a] == 32){
    int pdf_1[] = { 1,   2,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  dx  -> t  tx u  dx 
    int pdf_2[] = { 1,   4,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  sx  -> t  tx c  sx 
    int pdf_3[] = {-1,   2,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx u   -> t  tx u  dx 
    int pdf_4[] = {-1,   4,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx c   -> t  tx c  sx 
  }
  else if (psi_no_map[i_a] == 34){
    int pdf_1[] = { 1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  ux  -> t  tx g  g  
    int pdf_2[] = { 1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  cx  -> t  tx g  g  
    int pdf_3[] = {-1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux u   -> t  tx g  g  
    int pdf_4[] = {-1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx c   -> t  tx g  g  
  }
  else if (psi_no_map[i_a] == 35){
    int pdf_1[] = { 1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  ux  -> t  tx d  dx 
    int pdf_2[] = { 1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  cx  -> t  tx s  sx 
    int pdf_3[] = {-1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux u   -> t  tx d  dx 
    int pdf_4[] = {-1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx c   -> t  tx s  sx 
  }
  else if (psi_no_map[i_a] == 36){
    int pdf_1[] = { 1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  ux  -> t  tx u  ux 
    int pdf_2[] = { 1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  cx  -> t  tx c  cx 
    int pdf_3[] = {-1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux u   -> t  tx u  ux 
    int pdf_4[] = {-1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx c   -> t  tx c  cx 
  }
  else if (psi_no_map[i_a] == 37){
    int pdf_1[] = { 1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  ux  -> t  tx s  sx 
    int pdf_2[] = { 1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  cx  -> t  tx d  dx 
    int pdf_3[] = {-1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux u   -> t  tx s  sx 
    int pdf_4[] = {-1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx c   -> t  tx d  dx 
  }
  else if (psi_no_map[i_a] == 38){
    int pdf_1[] = { 1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  ux  -> t  tx c  cx 
    int pdf_2[] = { 1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  cx  -> t  tx u  ux 
    int pdf_3[] = {-1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux u   -> t  tx c  cx 
    int pdf_4[] = {-1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx c   -> t  tx u  ux 
  }
  else if (psi_no_map[i_a] == 39){
    int pdf_1[] = { 1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  ux  -> t  tx b  bx 
    int pdf_2[] = { 1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  cx  -> t  tx b  bx 
    int pdf_3[] = {-1,   2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux u   -> t  tx b  bx 
    int pdf_4[] = {-1,   4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx c   -> t  tx b  bx 
  }
  else if (psi_no_map[i_a] == 40){
    int pdf_1[] = { 1,   2,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  sx  -> t  tx u  sx 
    int pdf_2[] = { 1,   4,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  dx  -> t  tx c  dx 
    int pdf_3[] = {-1,   4,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx c   -> t  tx c  dx 
    int pdf_4[] = {-1,   2,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx u   -> t  tx u  sx 
  }
  else if (psi_no_map[i_a] == 42){
    int pdf_1[] = { 1,   2,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  cx  -> t  tx u  cx 
    int pdf_2[] = { 1,   4,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  ux  -> t  tx c  ux 
    int pdf_3[] = {-1,   4,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux c   -> t  tx c  ux 
    int pdf_4[] = {-1,   2,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx u   -> t  tx u  cx 
  }
  else if (psi_no_map[i_a] == 43){
    int pdf_1[] = { 1,   2,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // u  bx  -> t  tx u  bx 
    int pdf_2[] = { 1,   4,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // c  bx  -> t  tx c  bx 
    int pdf_3[] = {-1,   2,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // bx u   -> t  tx u  bx 
    int pdf_4[] = {-1,   4,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // bx c   -> t  tx c  bx 
  }
  else if (psi_no_map[i_a] == 44){
    int pdf_1[] = { 1,   5,   5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // b  b   -> t  tx b  b  
  }
  else if (psi_no_map[i_a] == 45){
    int pdf_1[] = { 1,   5,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // b  dx  -> t  tx b  dx 
    int pdf_2[] = { 1,   5,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // b  sx  -> t  tx b  sx 
    int pdf_3[] = {-1,   5,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // dx b   -> t  tx b  dx 
    int pdf_4[] = {-1,   5,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // sx b   -> t  tx b  sx 
  }
  else if (psi_no_map[i_a] == 46){
    int pdf_1[] = { 1,   5,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // b  ux  -> t  tx b  ux 
    int pdf_2[] = { 1,   5,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // b  cx  -> t  tx b  cx 
    int pdf_3[] = {-1,   5,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux b   -> t  tx b  ux 
    int pdf_4[] = {-1,   5,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx b   -> t  tx b  cx 
  }
  else if (psi_no_map[i_a] == 47){
    int pdf_1[] = { 1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // b  bx  -> t  tx g  g  
    int pdf_2[] = {-1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // bx b   -> t  tx g  g  
  }
  else if (psi_no_map[i_a] == 48){
    int pdf_1[] = { 1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // b  bx  -> t  tx d  dx 
    int pdf_2[] = { 1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // b  bx  -> t  tx s  sx 
    int pdf_3[] = {-1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // bx b   -> t  tx d  dx 
    int pdf_4[] = {-1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // bx b   -> t  tx s  sx 
  }
  else if (psi_no_map[i_a] == 49){
    int pdf_1[] = { 1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // b  bx  -> t  tx u  ux 
    int pdf_2[] = { 1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // b  bx  -> t  tx c  cx 
    int pdf_3[] = {-1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // bx b   -> t  tx u  ux 
    int pdf_4[] = {-1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // bx b   -> t  tx c  cx 
  }
  else if (psi_no_map[i_a] == 50){
    int pdf_1[] = { 1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // b  bx  -> t  tx b  bx 
    int pdf_2[] = {-1,   5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // bx b   -> t  tx b  bx 
  }
  else if (psi_no_map[i_a] == 51){
    int pdf_1[] = { 1,  -1,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // dx dx  -> t  tx dx dx 
    int pdf_2[] = { 1,  -3,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // sx sx  -> t  tx sx sx 
  }
  else if (psi_no_map[i_a] == 52){
    int pdf_1[] = { 1,  -1,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // dx ux  -> t  tx dx ux 
    int pdf_2[] = { 1,  -3,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // sx cx  -> t  tx sx cx 
    int pdf_3[] = {-1,  -1,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux dx  -> t  tx dx ux 
    int pdf_4[] = {-1,  -3,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx sx  -> t  tx sx cx 
  }
  else if (psi_no_map[i_a] == 53){
    int pdf_1[] = { 1,  -1,  -3};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // dx sx  -> t  tx dx sx 
    int pdf_2[] = { 1,  -3,  -1};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // sx dx  -> t  tx dx sx 
  }
  else if (psi_no_map[i_a] == 54){
    int pdf_1[] = { 1,  -1,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // dx cx  -> t  tx dx cx 
    int pdf_2[] = { 1,  -3,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // sx ux  -> t  tx ux sx 
    int pdf_3[] = {-1,  -3,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // ux sx  -> t  tx ux sx 
    int pdf_4[] = {-1,  -1,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // cx dx  -> t  tx dx cx 
  }
  else if (psi_no_map[i_a] == 56){
    int pdf_1[] = { 1,  -1,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // dx bx  -> t  tx dx bx 
    int pdf_2[] = { 1,  -3,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // sx bx  -> t  tx sx bx 
    int pdf_3[] = {-1,  -1,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // bx dx  -> t  tx dx bx 
    int pdf_4[] = {-1,  -3,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // bx sx  -> t  tx sx bx 
  }
  else if (psi_no_map[i_a] == 57){
    int pdf_1[] = { 1,  -2,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // ux ux  -> t  tx ux ux 
    int pdf_2[] = { 1,  -4,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // cx cx  -> t  tx cx cx 
  }
  else if (psi_no_map[i_a] == 58){
    int pdf_1[] = { 1,  -2,  -4};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // ux cx  -> t  tx ux cx 
    int pdf_2[] = { 1,  -4,  -2};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // cx ux  -> t  tx ux cx 
  }
  else if (psi_no_map[i_a] == 59){
    int pdf_1[] = { 1,  -2,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // ux bx  -> t  tx ux bx 
    int pdf_2[] = { 1,  -4,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_2, 3));   // cx bx  -> t  tx cx bx 
    int pdf_3[] = {-1,  -2,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_3, 3));   // bx ux  -> t  tx ux bx 
    int pdf_4[] = {-1,  -4,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_4, 3));   // bx cx  -> t  tx cx bx 
  }
  else if (psi_no_map[i_a] == 60){
    int pdf_1[] = { 1,  -5,  -5};    osi_combination_pdf.push_back(get_vector_from_array_int(pdf_1, 3));   // bx bx  -> t  tx bx bx 
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
