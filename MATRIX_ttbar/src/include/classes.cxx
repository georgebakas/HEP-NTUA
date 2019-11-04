#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <map>
#include <complex>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <time.h>
#include <dirent.h>
#include <cstdlib>
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <limits>
using namespace std;
typedef complex<double> double_complex;
#include "LHAPDF/LHAPDF.h"
#include "classes.h"
#include "routines.h"
#include "phasespace.h"
#include "varia.h"
#include "dipolesubtraction.h"
#include "qTsubtraction.h"
#include "NJsubtraction.h"
#include "const.cxx"
#include "gsl/gsl_sf_dilog.h"
#include "gsl/gsl_minmax.h"
#include "openloops.h"
using namespace mathconst;
using namespace physconst;

#define setdr right << setw(25) << setprecision(16) << showpoint
#define setdl left << setw(25) << setprecision(16) << showpoint
#define setsr right << setw(15) << setprecision(8) << showpoint
#define setsl left << setw(15) << setprecision(8) << showpoint

/*
extern "C" {
  void maincall_();
  void initialization_(double *mwin, double *mzin, double *wwin, double *wzin, double *mhin, double *whin, double *gfin, double *alfasin, double *mudim, int *scheme);
  void ol_start();
  void ol_finish();
  int ol_register_process(char* name, int amptype);
  int ol_n_external(int id);
  void ol_phase_space_point(int id, double sqrt_s, double* pp);
  void ol_evaluate_tree(int id, double* P, double* res);
  void ol_evaluate_full(int id, double* P, double* m2l0, double* m2l1, double* ir1, double* m2l2, double* ir2, double* acc);
  void ol_evaluate_loop(int id, double* pp, double* m2_tree, double* m2_loop, double* acc);
  void ol_evaluate_loop2(int id, double* P, double* m2l2, double* acc);
  void ol_evaluate_ct(int id, double* P, double* m2l0, double* m2ct);
  //  void ol_evaluate_cc(int id, double* P, double* m2cc);
  void ol_evaluate_cc(int id, double* psp, double* tree, double* cc, double* ewcc);
  void ol_evaluate_sc(int id, double* P, int emitter, double* polvect, double* m2sc);
  void ol_evaluate_cc2(int id, double* psp, double* tree, double* cc, double* ewcc);
  void ol_evaluate_sc2(int id, double* P, int emitter, double* polvect, double* m2sc);
  void ol_evaluate_loopcc(int id, double* psp, double* tree, double* cc, double* ewcc);
  void ol_setparameter_int_c_(char* param, int* value, int* error);
  void ol_setparameter_double_c_(char* param, double* value, int* error);
  void ol_setparameter_string_c_(char* param, char* value, int* error);
  void ol_setparameter_double(char* param, double value);
  void ol_getparameter_double(char* param, double* value);
  void ol_setparameter_int(char* param, int value);
  void ol_getparameter_int(char* param, int* value);
  void ol_setparameter_string(char* param, char* value);

  void ol_tree_colbasis_dim(int id, int* ncolb, int* colelemsz, int* nheltot);
  void ol_tree_colbasis(int id, int* basis, int* needed);
  void ol_evaluate_tree_colvect(int id, double* pp, double* amp, int* nhel);

  void OLP_PrintParameter(char* name);

  void ol_parameters_flush_();
  void parameters_init_(double *Mass_E, double *Mass_M, double *Mass_L, double *Mass_U, double *Mass_D, double *Mass_S, double *Mass_C, double *Width_C, double *Mass_B, double *Width_B, double *Mass_T, double *Width_T, double *Mass_W, double *Width_W, double *Mass_Z, double *Width_Z, double *Mass_H, double *Width_H, double *Coupl_Alpha_QED, double *Coupl_Alpha_QCD, int *last_switch, int *amp_switch, int *amp_switch_rescue, int *use_coli_cache, int *check_Ward_tree, int *check_Ward_loop, int *out_symmetry, int *leading_colour);
  void olloopparametersinit_(double *pole1_UV, double *pole1_IR, double *pole2_IR, int *N_quarks, int *N_light, int *CT_on, int *R2_on, int *IR_on, double *renscale, int *polenorm_swi, int *polecheck);
  void olloopparametersdeltascale_(double *pole1_UV, double *pole1_IR, double *pole2_IR, double *renscale);
  void olloopparametersscales_(double *renscale, double *fact_UV, double *fact_IR);
  void olparametersleadingcolour_(int *leading_colour);


  /*
  void H1qqdeltaext(double* P, double* H1);
  void xtest_();
*//*
}
  */
