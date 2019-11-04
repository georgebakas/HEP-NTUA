#include "../include/classes.cxx"
/*
#include <assert.h>
#include <iostream>

#include "header/hpl.h"

using namespace std;
*/
struct fort_complex
{double dr; double di;};

extern"C" {
void hplog_(double *x, int *nw, fort_complex *Hc1, fort_complex *Hc2, fort_complex *Hc3, fort_complex *Hc4, double *Hr1, double *Hr2, double *Hr3, double *Hr4, double *Hi1, double *Hi2, double *Hi3, double *Hi4, int *n1, int *n2);
}

extern"C" {
void tdhpl_(double *y, double *z, int *nmax, double *GYZ1, double *GYZ2, double *GYZ3, double *GYZ4, double *HZ1, double *HZ2, double *HZ3, double *HZ4);
}

//extern"C" {
//void tdhpl_(double *y, double *z, int *nmax, double [], double [][4], double [][4][4], double [][4][4][4], double [], double [][2], double [][2][2], double [][2][2][2]);
//}

int test_hpl()
{
    double x=0.5;
    int nw=2;
    fort_complex Hc1[2],Hc2[2*2],Hc3[2*2*2],Hc4[2*2*2*2];
    double Hr1[2],Hr2[2*2],Hr3[2*2*2],Hr4[2*2*2*2];
    double Hi1[2],Hi2[2*2],Hi3[2*2*2],Hi4[2*2*2*2];
    int n1=0;
    int n2=1;


   hplog_(&x, &nw, Hc1, Hc2, Hc3, Hc4, Hr1, Hr2, Hr3, Hr4, Hi1, Hi2, Hi3, Hi4, &n1, &n2);

   cout << x << ": " << Hc1[0].dr << ", " << Hr1[0] << endl;

   return 0;
}

void compute_HPLs_weight2(double v, double *HZ1, double *HZ2, double *HZ3) {
  int nw=2;
  fort_complex Hc1[2],Hc2[2*2],Hc3[2*2*2],Hc4[2*2*2*2];
  double Hr1[2],Hr2[2*2],Hr3[2*2*2],Hr4[2*2*2*2];
  double Hi1[2],Hi2[2*2],Hi3[2*2*2],Hi4[2*2*2*2];
  int n1=0;
  int n2=1;

  hplog_(&v, &nw, Hc1, Hc2, Hc3, Hc4, Hr1, Hr2, Hr3, Hr4, Hi1, Hi2, Hi3, Hi4, &n1, &n2);

//  HPLs["H(0,v)"] = double_complex(Hc1[0].dr,Hc1[0].di);
//  HPLs["H(0,0,v)"] = double_complex(Hc2[0*2+0].dr,Hc2[0*2+0].di);
//  HPLs["H(1,0,v)"] = double_complex(Hc2[0*2+1].dr,Hc2[0*2+1].di);
}

void compute_GHPLs(int weight, double u, double v, double *GYZ1, double *GYZ2, double *GYZ3, double *GYZ4, double *HZ1, double *HZ2, double *HZ3, double *HZ4) {
  int nw=weight;

  assert(u>=0);
  assert(u<=1-v);
  assert(v>=0);
  assert(v<=1);

  tdhpl_(&u, &v, &nw, GYZ1, GYZ2, GYZ3, GYZ4, HZ1, HZ2, HZ3, HZ4);

//  if (weight>2) {
//    fort_complex Hc1[2],Hc2[2*2],Hc3[2*2*2],Hc4[2*2*2*2];
//    double Hr1[2],Hr2[2*2],Hr3[2*2*2],Hr4[2*2*2*2];
//    double Hi1[2],Hi2[2*2],Hi3[2*2*2],Hi4[2*2*2*2];
//    int n1=0;
//    int n2=1;

//    hplog_(&v, &nw, Hc1, Hc2, Hc3, Hc4, Hr1, Hr2, Hr3, Hr4, Hi1, Hi2, Hi3, Hi4, &n1, &n2);
//    for (int i=0; i<2*2*2; i++) {
//        cout << Hc3[i].dr << ", " << Hc3[i].di << ", " << Hr3[i] << ", " << Hi3[i] << endl;
//        cout << HZ3[i] << endl;
//    }
//    cout << endl;
//  }

//  GHPLs["G(0,u)"] = GYZ1[0];
//  GHPLs["G(1,u)"] = GYZ1[1];
//  GHPLs["G(1 - v,u)"] = GYZ1[2];
//  GHPLs["G( - v,u)"] = GYZ1[3];
////  cout << log(1+u/v) << ", " << GYZ1[3] << endl;
////  cout << GYZ1[0] << ", " << log(u) << endl;
//
//  GHPLs["G(1,1 - v,u)"] = GYZ2[2*4+1];
//  GHPLs["G( - v,0,u)"] = GYZ2[0*4+3];

//  // only needed at NNLO
//  GHPLs["G(0,0,u)"] = GYZ2[0*4+0];
//  GHPLs["G(1,0,u)"] = GYZ2[0*4+1];
//  GHPLs["G(1 - v,0,u)"] = GYZ2[0*4+2];
//
//  GHPLs["G(0,1,u)"] = GYZ2[1*4+0];
//  GHPLs["G(1,1,u)"] = GYZ2[1*4+1];
//  GHPLs["G(1 - v,1,u)"] = GYZ2[1*4+2];
//  GHPLs["G( - v,1,u)"] = GYZ2[1*4+3];
//
//  GHPLs["G(0,1 - v,u)"] = GYZ2[2*4+0];
//  GHPLs["G(1,1 - v,u)"] = GYZ2[2*4+1];
//  GHPLs["G(1 - v,1 - v,u)"] = GYZ2[2*4+2];
//  GHPLs["G( - v,1 - v,u)"] = GYZ2[2*4+3];
//
//  GHPLs["G(0, - v,u)"] = GYZ2[3*4+0];
//  GHPLs["G(1, - v,u)"] = GYZ2[3*4+1];
//  GHPLs["G(1 - v, - v,u)"] = GYZ2[3*4+2];
//  GHPLs["G( - v, - v,u)"] = GYZ2[3*4+3];

//  HPLs["H(0,v)"] = HZ1[0];
//  HPLs["H(1,v)"] = HZ1[1];
//  HPLs["H(0,0,v)"] = HZ2[0*2+0];
//  HPLs["H(1,0,v)"] = HZ2[0*2+1];

//  cout << GHPLs.size() << endl;
//  cout << HPLs.size() << endl;

//  cout << GHPLs["G( - v,u)"] << ", " << GHPLs["G(0,u)"] << ", " << GHPLs["G( - v,0,u)"] << endl;
//  cout << HPLs["H(0,v)"] << ", " << HPLs["H(0,0,v)"] << ", " << HPLs["H(1,0,v)"] << endl;

//  GHPLs["G(0,1,1)"] = double_complex(GYZ2[1*4+0],0);
////  cout << GHPLs["G(0,1,1)"] << endl;
//
//  GHPLs["G(-1,1"] = double_complex(GYZ1[1],0);
//  cout << GHPLs["G(-1,1"] << endl;
//  cout << double_complex(GYZ1[0].dr,GYZ1[0].di) << endl;
//  cout << GYZ1[0] << endl;
//  cout << double_complex(GYZ1[2].dr,GYZ1[2].di) << endl;
//  cout << double_complex(GYZ1[3].dr,GYZ1[3].di) << endl;

//cout << "u=" << u << ", v=" << v << endl;
//
// for (int i=0; i<2; i++) {
//    cout << "H(" << i << "): " << HZ1[i] << endl;
//  }
//
//   for (int i=0; i<4; i++) {
//    cout << "G(" << i << "): " << GYZ1[i] << endl;
//  }
//
//  for (int i=0; i<2; i++) {
//    for (int j=0; j<2; j++) {
//    cout << "H(" << i << "," << j << "): " << HZ2[i+2*j] << endl;
//    }
//  }
//    for (int i=0; i<4; i++) {
//    for (int j=0; j<4; j++) {
//    cout << "G(" << i << "," << j << "): " << GYZ2[i+4*j] << endl;
//    }
//  }
//
//  for (int i=0; i<2; i++) {
//    for (int j=0; j<2; j++) {
//      for (int k=0; k<2; k++) {
//    cout << "H(" << i << "," << j << "," << k << "): " << HZ3[i+2*j+2*2*k] << endl;
//      }
//    }
//  }
//    for (int i=0; i<4; i++) {
//    for (int j=0; j<4; j++) {
//      for (int k=0; k<4; k++) {
//    cout << "G(" << i << "," << j << "," << k << "): " << GYZ3[i+4*j+4*4*k] << endl;
//      }
//    }
//  }
//
//  for (int i=0; i<2; i++) {
//    for (int j=0; j<2; j++) {
//      for (int k=0; k<2; k++) {
//        for (int l=0; l<2; l++) {
//    cout << "H(" << i << "," << j << "," << k << "," << l << "): " << HZ4[i+2*j+2*2*k+2*2*2*l] << endl;
//        }
//      }
//    }
//  }
//    for (int i=0; i<4; i++) {
//    for (int j=0; j<4; j++) {
//      for (int k=0; k<4; k++) {
//        for (int l=0; l<4; l++) {
//    cout << "G(" << i << "," << j << "," << k << "," << l << "): " << GYZ4[i+4*j+4*4*k+4*4*4*l] << endl;
//        }
//      }
//    }
//  }
//  assert(0);
}
