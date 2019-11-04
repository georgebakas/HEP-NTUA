#ifndef DIPOLE_SET_H
#define DIPOLE_SET_H

class dipole_set {
private:
  string _name;
  double _charge_factor;
  double _colour_factor;
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
  vector<int> _type_parton;
  vector<int> _basic_type_parton;
  int _type_dipole;
  int _type_splitting;
  int _no_R_emitter_1;
  int _no_R_emitter_2;
  int _no_R_spectator;
  int _no_A_emitter;
  int _no_A_spectator;
  int _binary_R_emitter_1;
  int _binary_R_emitter_2;
  int _binary_R_spectator;
  int _binary_A_emitter;
  int _binary_A_spectator;
  int _type_correction;
  int _massive;

public:
////////////////////
//  constructors  //
////////////////////
  dipole_set();

  dipole_set(const dipole_set & QEW_dipole_candidate, int no_map, vector<int> o_map, int no_prc, vector<int> o_prc, int n_channel, int sum_channel, double charge_factor, double symmetry_factor, int massive);

  dipole_set(string name, vector<int> type_parton, vector<int> basic_type_parton, int type_dipole, int type_splitting, int no_R_emitter_1, int no_R_emitter_2, int no_R_spectator, int no_A_emitter, int no_A_spectator, int type_correction);

  dipole_set(string name, vector<int> type_parton, vector<int> basic_type_parton, double symmetry_factor, int no_map, vector<int> o_map, int no_prc, vector<int> o_prc, int n_channel);

  dipole_set(string name, double charge_factor, double symmetry_factor, vector<vector<int> > dx_pa, vector<int> phasespace, vector<int> o_prc, vector<int> o_map, vector<vector<double> > colourmatrix, vector<vector<int> > spinorder, vector<int> fckm, vector<int> data, int no_prc, int no_map, int n_channel, double exp_pdf, vector<int> type_parton, vector<int> basic_type_parton, int type_dipole, int type_splitting, int no_R_emitter_1, int no_R_emitter_2, int no_R_spectator, int no_A_emitter, int no_A_spectator, int type_correction, int massive);
  

///////////////////////
//  access elements  //
///////////////////////
  string name() const;
  double charge_factor() const;
  double colour_factor() const;
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
  const vector<int> type_parton() const;
  const vector<int> basic_type_parton() const;
  int type_dipole() const;
  int type_splitting() const;
  int no_R_emitter_1() const;
  int no_R_emitter_2() const;
  int no_R_spectator() const;
  int no_A_emitter() const;
  int no_A_spectator() const;
  int binary_R_emitter_1() const;
  int binary_R_emitter_2() const;
  int binary_R_spectator() const;
  int binary_A_emitter() const;
  int binary_A_spectator() const;

  int type_correction() const;
  int massive() const;

  int no_BLHA_entry;
  int process_id;

  int contribution_order_alpha_s;
  int contribution_order_alpha_e;
  int contribution_order_interference;

  double xy;
  double zuv;
};
#endif
