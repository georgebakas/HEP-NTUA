#ifndef CORRELATIONOPERATOR_H
#define CORRELATIONOPERATOR_H

class correlationoperator {
private:

public:
////////////////////
//  constructors  //
////////////////////
  correlationoperator();
  correlationoperator(string name, int type, vector<int> pair, int no_prc, vector<int> type_parton, double charge_factor, int no_emitter, int no_spectator, int type_correction, int type_combination, int massive);

///////////////////////
//  access elements  //
///////////////////////

  string name;
  int type;
  vector<int> pair;
  int no_prc;
  vector<int> type_parton;
  double charge_factor;
  int no_emitter;
  int no_spectator;
  int type_combination;
  int type_correction;
  int massive;
  int no_OL_entry;
  int no_BLHA_entry;
  double colour_factor;
  int process_id;
};
#endif
