#ifndef MULTICOLLINEAR_SET_H
#define MULTICOLLINEAR_SET_H

class multicollinear_set {
private:
  /*
  string _name;
  vector<string> _all_name;
  int _type;
  vector<int> _in_collinear;
  int _no_prc;
  vector<int> _type_parton;
  vector<vector<int> > _pdf;
  //  vector<int> _pdf;
  double _charge_factor;
  int _no_emitter;
  int _no_spectator;
  vector<int> _pair;
  int _no_OL_entry;
  int _type_correction;
  int _massive;
  */

public:
  string name;
  vector<string> all_name;
  int no_prc;
  vector<int> type_parton;
  vector<vector<int> > pdf;
  int process_id;

  //  formerly "type":
  //  int type_splitting;
  vector<int> type_splitting;
  //  three entries for each step 0 (irregular part) 1,2 (emission from parton 1,2)
  //  vector<int> in_collinear;
  vector<vector<int> > in_collinear;
  vector<int> no_emitter;
  //  could be different if QCD and QEW emissions are combined !!!
  vector<int> type_correction;
  vector<int> massive;
  vector<int> n_emission;
  vector<int> emission;

  // maybe modify this part: spectators are method-specific !!!
  vector<int> no_spectator;
  vector<vector<int> > pair;
  vector<int> no_OL_entry;
  vector<int> no_BLHA_entry;

  // needs modification anyway...
  vector<double> charge_factor;

  vector<int> previous_splitting;

  vector<int> no_endpoint;
  vector<int> type_splitting_full;
  vector<vector<int> > type_splitting_leg;
  vector<vector<int> > n_emission_leg_type;
  int leg_emission;
  int no_pdf;
  //  vector<int> order_emitter;
  int x_a;
  int x_b;

////////////////////
//  constructors  //
////////////////////
  multicollinear_set();
  //  multicollinear_set(multicollinear_set ms);
  multicollinear_set(int _no_prc, vector<int> _type_parton, vector<vector<int> > _pdf);
  //  multicollinear_set(string _name, vector<string> _all_name, int _type_splitting, vector<int> _in_collinear, int _no_prc, vector<int> _type_parton, vector<vector<int> > _pdf, double _charge_factor, int _no_emitter, int _no_spectator, vector<int> _pair, int _type_correction, int _massive);
  multicollinear_set(string _name, vector<string> _all_name, vector<int> _type_splitting, vector<vector<int> > _in_collinear, int _no_prc, vector<int> _type_parton, vector<vector<int> > _pdf, vector<double> _charge_factor, vector<int> _no_emitter, vector<int> _no_spectator, vector<vector<int> > _pair, vector<int> _type_correction, vector<int> _massive);
///////////////////////
//  access elements  //
///////////////////////
/*
  string name() const;
  const vector<string> all_name() const;
  int type() const;
  const vector<int> in_collinear() const;
  int no_prc() const;
  const vector<int> type_parton() const;
  //  const vector<int> pdf() const;
  const vector<vector<int> > pdf() const;
  double charge_factor() const;
  int no_emitter() const;
  int no_spectator() const;
  int no_OL_entry() const;
  const vector<int> pair() const;

  int type_correction() const;
  int massive() const;
*/
};
#endif
