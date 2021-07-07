#ifndef RUNRESUMPTIONSET_H
#define RUNRESUMPTIONSET_H

#include <map>
#include <string>
#include <vector>

using namespace std;

class runresumption_set{
 private:

 public:
////////////////////
//  constructors  //
////////////////////
  runresumption_set();
  runresumption_set(observable_set & _osi, phasespace_set & _psi, observable_set & _save_osi, phasespace_set & _save_psi);

  void perform_proceeding_step();
  void perform_proceeding_in();
  void perform_proceeding_out();
  void perform_proceeding_check();

  void perform_iteration_step();

  int size_proc_generic;
  vector<int> size_proceeding;

  observable_set *resumption_osi;
  phasespace_set *resumption_psi;

  observable_set *osi;
  phasespace_set *psi;
};
#endif
