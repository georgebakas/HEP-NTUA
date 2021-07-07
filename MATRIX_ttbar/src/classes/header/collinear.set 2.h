#ifndef COLLINEAR_SET_H
#define COLLINEAR_SET_H

class collinear_set {
private:
  string _name;
  vector<string> _all_name;
  int _type;
  vector<int> _in_collinear;
  int _no_prc;
  vector<int> _type_parton;
  vector<vector<int> > _all_pdf;
  //  vector<int> _all_pdf;
  double _charge_factor;
  double _colour_factor;
  int _no_emitter;
  int _no_spectator;
  vector<int> _pair;
  int _no_OL_entry;
  int _type_correction;
  int _massive;

  /*
  double _symmetry_factor;
  vector<vector<int> > _dx_pa;
  vector<int> _phasespace;
  vector<int> _o_prc;
  vector<int> _o_map;
  vector<vector<double> > _colourmatrix;
  vector<vector<int> > _spinorder;
  vector<int> _fckm;
  vector<int> _data;
  int _no_prc;
  int _no_map;
  int _n_channel;
  int _sum_channel;
  double _exp_pdf;
*/
  //  vector<int> _type_parton;
  //  int _type_splitting;
  //  int _binary_emitter;
  //  int _binary_spectator;
  
public:
////////////////////
//  constructors  //
////////////////////
  collinear_set();
  collinear_set(string name, vector<string> all_name, int type, vector<int> in_collinear, int no_prc, vector<int> type_parton, vector<vector<int> > all_pdf, double charge_factor, int no_emitter, int no_spectator, vector<int> pair, int type_correction, int massive);
  //vector<vector<int> > all_pdf, 
  /*
  collinear_set(const collinear_set & QEW_collinear_candidate, int no_map, vector<int> o_map, int no_prc, vector<int> o_prc, int n_channel, int sum_channel, double charge_factor, double symmetry_factor);

  collinear_set(string name, vector<int> type_parton, int type, int type_splitting, int no_R_emitter_1, int no_R_emitter_2, int no_R_spectator, int no_A_emitter, int no_A_spectator);

  //  collinear_set(string name, vector<int> type_parton, vector<vector<int> > dx_pa, double charge_factor, double symmetry_factor, int no_map, vector<int> o_map, int no_prc, vector<int> o_prc, int n_channel, int type, int type_splitting, int no_R_emitter_1, int no_R_emitter_2, int no_R_spectator, int no_A_emitter, int no_A_spectator);
  collinear_set(string name, vector<int> type_parton, double symmetry_factor, int no_map, vector<int> o_map, int no_prc, vector<int> o_prc, int n_channel);
  //  collinear_set(string name, vector<int> type_parton, vector<vector<int> > dx_pa, double charge_factor, double symmetry_factor, int no_map, vector<int> o_map, int no_prc, vector<int> o_prc, int n_channel);

  collinear_set(string name, double charge_factor, double symmetry_factor, vector<vector<int> > dx_pa, vector<int> phasespace, vector<int> o_prc, vector<int> o_map, vector<vector<double> > colourmatrix, vector<vector<int> > spinorder, vector<int> fckm, vector<int> data, int no_prc, int no_map, int n_channel, double exp_pdf, vector<int> type_parton, int type, int type_splitting, int no_R_emitter_1, int no_R_emitter_2, int no_R_spectator, int no_A_emitter, int no_A_spectator);
  */

///////////////////////
//  access elements  //
///////////////////////
  string name() const;
  const vector<string> all_name() const;
  int type() const;
  const vector<int> in_collinear() const;
  int no_prc() const;
  const vector<int> type_parton() const;
  //  const vector<int> all_pdf() const;
  const vector<vector<int> > all_pdf() const;
  double charge_factor() const;
  double colour_factor() const;
  int no_emitter() const;
  int no_spectator() const;
  int no_OL_entry() const;
  const vector<int> pair() const;

  int type_correction() const;
  int massive() const;
  /*
  double symmetry_factor() const;
  const vector<vector<int> > dx_pa() const;
  const vector<int> phasespace() const;
  const vector<int> o_prc() const;
  const vector<int> o_map() const;
  const vector<vector<double> > colourmatrix() const;
  const vector<vector<int> > spinorder() const;
  const vector<int> fckm() const;
  const vector<int> data() const;
  int no_prc() const;
  int no_map() const;
  int n_channel() const;
  int sum_channel() const;
  double exp_pdf() const;
  int type_splitting() const;
  */

  /*
  int binary_emitter() const;
  int binary_spectator() const;
  */
  int no_BLHA_entry;
  int process_id;
};
#endif
