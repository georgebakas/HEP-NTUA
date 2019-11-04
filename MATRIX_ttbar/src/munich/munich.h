#ifndef MUNICH_H
#define MUNICH_H

//#include "../include/classes.cxx"

#include <string>
#include <vector>

//#include "randommanager.h"
#include "../classes/header/logger.h"
//#include "../../include/ftypes.h"
//#include "fourvector.h"
//#include "particle.h"
//#include "phasespace.set.h"
//#include "dddistribution.h"
//#include "xdistribution.h"

using namespace std;

//class summary_generic;

class munich{
 private:

 public:
////////////////////
//  constructors  //
////////////////////
  munich();
  munich(int argc, char *argv[], string basic_process_class);

  void run_initialization();
  void run_integration();

  void get_summary();

  void calculate_p_parton();

  void integration_born();

  void integration_VA_QCD();
  void integration_CA_QCD();
  void integration_RA_QCD();
  void integration_VA_QEW();
  void integration_CA_QEW();
  void integration_RA_QEW();
  void integration_VA_MIX();
  void integration_RA_MIX();

  void integration_VT_QCD();
  void integration_CT_QCD();
  void integration_VT2_QCD();
  void integration_CT2_QCD();

  void integration_CJ_QCD();
  void integration_VJ_QCD();
  void integration_VJ2_QCD();
  void integration_CJ2_QCD();

  //  void handling_cut_psp();
  void handling_vanishing_me2();
  void errorhandling_c_psp();
  void errorhandling_me2();
  void errorhandling_alpha_S();
  void errorhandling_pdf();
  void errorhandling_gtot();
  void errorhandling_OL();
  void errorhandling_cut_technical();
  void errorhandling_c_xbpsp();
  void errorhandling_c_psp_initial();
  void errorhandling_RA_me2(int xswitch);
  void errorhandling_RA_gtot();
  void errorhandling_collinear_me2();
  

  string subprocess;
  string order;
  string infilename;
  string infilename_scaleband;

  inputparameter_set isi;
  contribution_set csi;
  model_set msi;
  event_set esi;
  observable_set osi;
  phasespace_set psi;
  user_defined user;

  // Catani--Seymour subtraction
  vector<vector<ioperator_set> > ioperator;
  vector<vector<collinear_set> > collinear;
  vector<dipole_set> dipole;
  vector<dipole_set> dipole_candidate;

  call_generic generic;

  summary_generic ysi;

};

#endif
