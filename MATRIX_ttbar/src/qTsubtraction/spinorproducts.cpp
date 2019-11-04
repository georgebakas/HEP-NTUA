#include "../include/classes.cxx"
/*
#include<assert.h>
#include "header/amplitudes.h"

#include <iostream>
#include<vector>

//#include "../classes/header/fourvector.h"

#include "header/spinorproducts.h"

using namespace std;

*/

// takes a list of momenta p[i] and calculates the spinor products <ij>=za[ij], [ij]=zb[ij] and the invariants s_ij
//const
void calcSpinorProducts(vector<fourvector> &p, vector<vector<double_complex > > &za, vector<vector<double_complex > > &zb, vector<vector <double> > &s) {
  static Logger logger("calcSpinorProducts");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  za.resize(p.size());
  zb.resize(p.size());
  s.resize(p.size());

  vector<double > sqrt_pplus(p.size());
  vector<double_complex > pplustrans(p.size());
  vector<double_complex > c_sign(p.size());

  for (int i=0; i<p.size(); i++) {
    logger << LOG_DEBUG_VERBOSE << "i = " << i << endl;
    za[i].resize(p.size());
    zb[i].resize(p.size());
    s[i].resize(p.size());

    if (p[i].x0()>=0) { // positive energy case
      sqrt_pplus[i]=sqrt(p[i].x0()+p[i].x1());
      pplustrans[i]=complex<double>(p[i].x3(),-p[i].x2());

//      sqrt_pplus[i]=sqrt(p[i].x0()+p[i].x2());
//      pplustrans[i]=complex<double>(p[i].x1(),-p[i].x3());
//
//      sqrt_pplus[i]=sqrt(p[i].x0()+p[i].x3());
//      pplustrans[i]=complex<double>(p[i].x1(),p[i].x2());

      c_sign[i]=1;
    } else { // negative energy case
      sqrt_pplus[i]=sqrt(-p[i].x0()-p[i].x1());
      pplustrans[i]=complex<double>(-p[i].x3(),p[i].x2());

//      sqrt_pplus[i]=sqrt(-p[i].x0()-p[i].x2());
//      pplustrans[i]=complex<double>(-p[i].x1(),p[i].x3());
//
//      sqrt_pplus[i]=sqrt(-p[i].x0()-p[i].x3());
//      pplustrans[i]=complex<double>(-p[i].x1(),-p[i].x2());

      c_sign[i]=double_complex(0,1);
    }
  }

  for (int i=1; i<p.size(); i++) {
    for (int j=1; j<i; j++) {
      s[i][j]=2.0*p[i]*p[j];

      if (sqrt_pplus[i]!=0 && sqrt_pplus[j]!=0) {
        za[i][j]=pplustrans[i]*sqrt_pplus[j]/sqrt_pplus[i]-pplustrans[j]*sqrt_pplus[i]/sqrt_pplus[j];
        za[i][j] = c_sign[i]*c_sign[j]*za[i][j];
      } else if (sqrt_pplus[j]!=0) {
        za[i][j]=sqrt(2*p[i].x0())*sqrt_pplus[j]-pplustrans[j]*sqrt_pplus[i]/sqrt_pplus[j];
        za[i][j] = c_sign[i]*c_sign[j]*za[i][j];
      } else if (sqrt_pplus[i]!=0) {
        za[i][j]=pplustrans[i]*sqrt_pplus[j]/sqrt_pplus[i]-sqrt(2*p[j].x0())*sqrt_pplus[i];
        za[i][j] = c_sign[i]*c_sign[j]*za[i][j];
      } else {
        cout << i << ", " << j << endl;
        assert(false);
      }

      if (fabs(s[i][j])>1e-9) {
        zb[i][j]=-s[i][j]/za[i][j];
      } else {
        zb[i][j] = -c_sign[i]*c_sign[j]*c_sign[i]*c_sign[j]*conj((za[i][j]));
      }

      za[j][i]=-za[i][j];
      zb[j][i]=-zb[i][j];
      s[j][i]=s[i][j];
    }
  }
//  for (int i=2; i<p.size(); i++) {
//    for (int j=1; j<i; j++) {
//      cout << i << "," << j << ": " << za[i][j] << ", " << zb[i][j] << endl;
//    }
//  }
//
//  // test Schouten identities
//  for (int i=1; i<p.size(); i++) {
//    for (int j=1; j<p.size(); j++) {
//      for (int k=1; k<p.size(); k++) {
//        for (int l=1; l<p.size(); l++) {
//          double schouten = abs(za[i][j]*za[k][l]-za[i][k]*za[j][l]-za[i][l]*za[k][j]);
//          assert(schouten<1e-6);
//        }
//      }
//    }
//  }
//
//  // test momentum conservation
//  cout << "input: " << (p[1]+p[2]+p[5]-p[3]-p[4]).r() << endl;;
//
//  for (int j=1; j<p.size(); j++) {
//    for (int k=1; k<p.size(); k++) {
//      double_complex momconv=0;
//      momconv=zb[j][1]*za[1][k]+zb[j][2]*za[2][k];
//      for (int i=3; i<p.size(); i++) {
//        if (p[i].x0()<0)
//          momconv+=zb[j][i]*za[i][k];
//        else
//          momconv-=zb[j][i]*za[i][k];
//      }
//      cout << "momentum conservation " << j << ", " << k << ": " << abs(momconv) << endl;
//      assert(abs(momconv)<1e-6);
//    }
//  }
}
