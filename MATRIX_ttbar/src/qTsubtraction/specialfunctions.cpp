#include "../include/classes.cxx"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_integration.h"

// multiplied by z!
double zbesselk0(double z) {
    return z*gsl_sf_bessel_Knu(1.0,z);
}

double zbesselk1(double z) {
    return gsl_sf_bessel_Knu(0.0,z);
}

double zbesselk2(double z) {
  // TODO: implement more efficient version
  const double a[14] = { 1.15443132980306572,
      1.97811199065594511,
      0.154431329803065721,
      4.801792651508824500,
      0.806235643470665767,
      -0.672784335098467139,
      3.285072828402112960,
      -1.945338757678943440,
      -0.181575166960855634,
      0.694195147571435559,
      -0.607655744858515573,
      -0.019182189839330562,
      0.068894530444636532,
      -0.070514317816328185 };

  const double r[12] = { 3.19461756880e4,
      -5.82903207466e3,
      1.17069096329e3,
      -2.61456867712e2,
      6.57620334072e1,
      -18.9305966582,
      6.37010269165,
      -2.57905883789,
      1.30957031250,
      -0.888020833333,
      0.875,
      1 };
//double tmp;
  // TODO: check!
  if (z<2.0) {
      double zm2=z*z/4;
      double loz=log(z/2);

      return loz*loz+a[0]*loz+a[1]+zm2*(2.0*loz*loz*loz/3+a[2]*loz*loz+a[3]*loz+a[4])+pow(zm2,2)*(loz*loz*loz/3+a[5]*loz*loz+a[6]*loz+a[7])+pow(zm2,3)*(loz*loz*loz/18+a[8]*loz*loz+a[9]*loz+a[10])+pow(zm2,4)*(loz*loz*loz/216+a[11]*loz*loz+a[12]*loz+a[13]);

//  zbesselk2=loz**2+a(0)*loz+a(1)
//     &       +zm**2*(2*loz**3/3d0+a(2)*loz**2+a(3)*loz+a(4))
//     &       +zm**4*(loz**3/3d0+a(5)*loz**2+a(6)*loz+a(7))
//     &       +zm**6*(loz**3/18d0+a(8)*loz**2+a(9)*loz+a(10))
//     &       +zm**8*(loz**3/216d0+a(11)*loz**2+a(12)*loz+a(13))
  } else if (z>4.0) {
      return sqrt(M_PI/2)*pow(z,-11.5)*exp(-z)*(r[0]+r[1]*z+r[2]*z*z+r[3]*z*z*z+r[4]*pow(z,4)+r[5]*pow(z,5)+r[6]*pow(z,6)+r[7]*pow(z,7)+r[8]*pow(z,8)+r[9]*pow(z,9)+r[10]*pow(z,10)+r[11]*pow(z,11));
//      zbesselk2=dsqrt(Pi/2d0)*(z)**(-11.5d0)*dexp(-z)*
//     &    (r(0)+r(1)*z+r(2)*z**2+r(3)*z**3+r(4)*z**4+r(5)*z**5+
//     &     r(6)*z**6+r(7)*z**7+r(8)*z**8+r(9)*z**9+r(10)*z**10+
//     &     r(11)*z**11)
  }


  double h=1e-5;

  // error seems to be of the order of 10^-6
  return z*(gsl_sf_bessel_Knu(1.0-h,z)-2.0*gsl_sf_bessel_Knu(1.0,z)+gsl_sf_bessel_Knu(1.0+h,z))/(h*h);

//  if (z>4.0)
//    cout << tmp << ", " << tmp2 << endl;
//    return tmp2;
}

double zbesselk3(double z) {
  // FIXME: implement more efficient version

   const double b[15] = { 1.731646994704598580,
      5.934335971967835330,
      5.444874456485317730,
      -1.268353005295401420,
      8.471041982558638170,
      -3.026167526073320430,
      -0.692088251323850355,
      2.809848746963509900,
      -2.161466255000085060,
      -0.104676472369316706,
      0.381989731242156681,
      -0.367492827636283900,
      -0.007844362856415627,
      0.027796539630842606,
      -0.029917436634978395 };

  const double r[10] = {-3.19152148877e3,
      7.05641513542e2,
      -1.75295543138e2,
      4.96775524480e1,
      -1.63798988342e1,
      6.45276489258,
      -3.1533203125,
      2.0234375,
      -1.875,
      3 };

      if(z<2.0) {
        double zm2=z*z/4;
        double loz=log(z/2);

        return -(pow(loz,3)+b[0]*loz*loz+b[1]*loz+b[2]+zm2*(pow(loz,3)+b[3]*loz*loz+b[4]*loz+b[5])+pow(zm2,2)*(pow(loz,3)/4+b[6]*loz*loz+b[7]*loz+b[8])+pow(zm2,3)*(pow(loz,3)/36+b[9]*loz*loz+b[10]*loz+b[11])+pow(zm2,4)*(pow(loz,3)/576+b[12]*loz*loz+b[13]*loz+b[14]));
      } else if (z>4.0) {
        // tested!
        return sqrt(M_PI/2)*pow(z,-10.5)*exp(-z)*(r[0]+r[1]*z+r[2]*z*z+r[3]*z*z*z+r[4]*pow(z,4)+r[5]*pow(z,5)+r[6]*pow(z,6)+r[7]*pow(z,7)+r[8]*pow(z,8)+r[9]*pow(z,9));
      }


  double h=5e-4;

  // error seems to be of the order of 10^-6
  return z*(-0.5*gsl_sf_bessel_Knu(1.0-2.0*h,z)+gsl_sf_bessel_Knu(1.0-h,z)-gsl_sf_bessel_Knu(1.0+h,z)+0.5*gsl_sf_bessel_Knu(1.0+2.0*h,z))/(h*h*h);
}

void Itilde(const double xmio, const int order, double &LL1, double &LL2, double &LL3, double &LL4) {
  //    const double Eulergamma=0.577215664902;
  //    const double z2=1.64493406685;
  //    const double z3=1.20205690316;
  //    const double b0_Eulergamma=2*exp(-Eulergamma);

  const double z2 = zeta2;
  const double z3 = zeta3;
  //const double b0 = b0_Eulergamma;

    double argum=b0_Eulergamma*xmio;
    double logx=log(xmio);

    double zbk0,zbk1,zbk2,zbk3;
    zbk0=zbesselk0(argum);
    zbk1=zbesselk1(argum);
    zbk2=zbesselk2(argum);
    zbk3=zbesselk3(argum);

    LL1=-zbk0/pow(xmio,2);
    LL2= 2.0/pow(xmio,2)*(zbk0*logx-zbk1);
    if (order==2) {
      LL3=-3.0/pow(xmio,2)*(zbk0*(logx*logx-z2)-2*zbk1*logx+zbk2);
      LL4=4.0/pow(xmio,2)*(zbk0*(logx*logx*logx-3*z2*logx+2*z3)-3*zbk1*(logx*logx-z2)+3*zbk2*logx-zbk3);
    }
}

double myli3(double x) {
  //  const double Z3=1.20205690315959429;
  const double Z3 = zeta3;

  double myLI3;

  if (x<-20) { // asymptotic expansion, suggested by Mathematica
    myLI3=-1.0/6*(pow(log(-x),3)+M_PI*M_PI*log(-x))+1.0/x+1.0/8/x/x+1.0/27/x/x/x;
  } else if (x<0.5) {
    double xlog=log(1-x);
    myLI3 = -xlog-(3*xlog*xlog)/8-17*pow(xlog,3)/216-5*pow(xlog,4)/576-7*pow(xlog,5)/54000+7*pow(xlog,6)/86400+19*pow(xlog,7)/5556600-pow(xlog,8)/752640-11*pow(xlog,9)/127008000+11*pow(xlog,10)/435456000;
  } else if (x<1.0) {
    double xlog=log(x);
    myLI3 = Z3+M_PI*M_PI*xlog/6+(3.0/4-log(-xlog)/2)*xlog*xlog-pow(xlog,3)/12-pow(xlog,4)/288+pow(xlog,6)/86400-pow(xlog,8)/10160640;
  } else if (x==1.0) {
    myLI3 = Z3;
  } else {
    cout << "myli3: wrong argument x=" << x << endl;
    assert(false);
    return 0.0;
  }

  return myLI3;
}

double I_tilde1(double xmio) {
  //    const double Eulergamma=0.577215664902;
  //    const double b0_Eulergamma=2*exp(-Eulergamma);
  //const double b0 = b0_Eulergamma;

    xmio = 10.;

    double argum;
    argum=b0_Eulergamma*xmio;

    double zbk0;
    zbk0=zbesselk0(argum);

    //    cout << "xmio = " << xmio << " -> " << setprecision(15) << setw(23) << -zbk0/pow(xmio,2) << endl;
    //    cout << "xmio = " << xmio << " -> Itilde1 = " << setprecision(15) << setw(23) << -zbk0/pow(xmio,2)*2*xmio << endl;
 
    return -zbk0/pow(xmio,2)*2*xmio;
}

double I_tilde2(double xmio) {
  //    const double Eulergamma=0.577215664902;
  //    const double b0_Eulergamma=2*exp(-Eulergamma);
  //const double b0 = b0_Eulergamma;

    xmio = 10.;

    double argum;
    double logx;
    argum=b0_Eulergamma*xmio;
    logx=log(xmio);

    double zbk0,zbk1;
    zbk0=zbesselk0(argum);
    zbk1=zbesselk1(argum);

    //    cout << "xmio = " << xmio << " -> " << setprecision(15) << setw(23) << 2.0/pow(xmio,2)*(zbk0*logx-zbk1) << endl;
    //    cout << "xmio = " << xmio << " -> Itilde2 = " << setprecision(15) << setw(23) << 2.0/pow(xmio,2)*(zbk0*logx-zbk1)*2*xmio << endl;
 
    return 2.0/pow(xmio,2)*(zbk0*logx-zbk1)*2*xmio;
}

double Itilde1(double xmio, void *params) {
  double argum = b0_Eulergamma * xmio;
  double zbk0 = zbesselk0(argum);
  return -zbk0 / pow(xmio, 2) * 2 * xmio;
}

double Itilde2(double xmio, void *params) {
  double argum = b0_Eulergamma * xmio;
  double logx = log(xmio);
  double zbk0 = zbesselk0(argum);
  double zbk1 = zbesselk1(argum);
  return 2. / pow(xmio, 2) * (zbk0 * logx - zbk1) * 2 * xmio;
}

double Itilde3(double xmio, void *params) {
  double argum = b0_Eulergamma * xmio;
  double logx = log(xmio);
  double zbk0 = zbesselk0(argum);
  double zbk1 = zbesselk1(argum);
  double zbk2 = zbesselk2(argum);
  return -3. / pow(xmio, 2) * (zbk0 * (logx * logx - zeta2) - 2 * zbk1 * logx + zbk2) * 2 * xmio;
}

double Itilde4(double xmio, void *params) {
  double argum = b0_Eulergamma * xmio;
  double logx = log(xmio);

  //    const double Eulergamma=0.577215664902;
  //    const double z2=1.64493406685;
  //    const double z3=1.20205690316;
  //  const double b0_Eulergamma=2*exp(-Eulergamma);
  const double z2 = zeta2;
  const double z3 = zeta3;
  //const double b0 = b0_Eulergamma;


  /*
    double argum;
    double logx;
    argum=b0_Eulergamma*xmio;
    logx=log(xmio);
  */

  //    double zbk0,zbk1,zbk2,zbk3;
  double zbk0 = zbesselk0(argum);
  double zbk1 = zbesselk1(argum);
  double zbk2 = zbesselk2(argum);
  double zbk3 = zbesselk3(argum);
  return 4. / pow(xmio, 2) * (zbk0 * (logx * logx * logx - 3 * z2 * logx + 2 * z3) - 3 * zbk1 * (logx * logx - z2) + 3 * zbk2 * logx - zbk3) * 2 * xmio;
}

void computeItildaIntegrals(double xmin, double xstep, int steps, vector<double> &I1_int, vector<double> &I2_int, vector<double> &I3_int, vector<double> &I4_int) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);

  double result, error;

  gsl_function F;
  F.params = NULL;

  //  double xmax=60;
  //  double xmax=500; // unchanged I1, I2, ...
  //  double xmax=600; // unchanged I1, I2, ...
  double xmax=600.; // unchanged I1, I2, ...
  //  double xmax=20; // minimal changes I1, I2, ...

//  for (int i=0; i<steps; i++) {
//    double xcut=xmin+i*xstep;
//
//    cout << "xcut = " << xcut << endl;
//
//    F.function = &Itilde1;
//    gsl_integration_qags (&F, xcut, xmax, 0, 1e-7, 10000,
//                        w, &result, &error);
//    I1_int[i]=result;
//    cout << "I1: " << result << " +- " << error << endl;
//
//    F.function = &Itilde2;
//    gsl_integration_qags (&F, xcut, xmax, 0, 1e-7, 10000,
//                        w, &result, &error);
//    I2_int[i]=result;
//    cout << "I2: " << result << " +- " << error << endl;
//
//    F.function = &Itilde3;
//    gsl_integration_qags (&F, xcut, xmax, 0, 1e-7, 10000,
//                        w, &result, &error);
//    I3_int[i]=result;
//    cout << "I3: " << result << " +- " << error << endl;
//
//    F.function = &Itilde4;
//    gsl_integration_qag (&F, xcut, xmax, 0, 1e-8, 10000,6,
//                        w, &result, &error);
//    I4_int[i]=result;
//    cout << "I4: " << result << " +- " << error << endl;
//  }

  cout << "precomputing Itilde integrals" << endl;

  double error_I1=0,error_I2=0,error_I3=0,error_I4=0;

  for (int i=steps-1; i>=0; i--) {
    double xcut=xmin+i*xstep;

    cout << "xcut = " << xcut << endl;

    if (i<steps-1)
      xmax=xmin+(i+1)*xstep;

    F.function = &Itilde1;
    gsl_integration_qags (&F, xcut, xmax, 0, 1e-7, 10000,
                        w, &result, &error);
    I1_int[i]=result;
    if (i<steps-1) {
      I1_int[i] += I1_int[i+1];
    }
    error_I1=sqrt(error_I1*error_I1+error*error);

    //    cout << "I1 (rcut= " << xcut << ": " << setprecision(15) << setw(23) << I_tilde1(xcut) << endl;

    cout << "I1: " << setprecision(15) << setw(23) << I1_int[i] << " +- " << error_I1 << ", " << error << endl;

    F.function = &Itilde2;
    gsl_integration_qags (&F, xcut, xmax, 0, 1e-7, 10000,
                        w, &result, &error);
    I2_int[i]=result;
    if (i<steps-1) {
      I2_int[i] += I2_int[i+1];
    }
    error_I2=sqrt(error_I2*error_I2+error*error);

    //    cout << "I2 (rcut= " << xcut << ": " << setprecision(15) << setw(23) << I_tilde2(xcut) << endl;

    cout << "I2: " << setprecision(15) << setw(23) << I2_int[i] << " +- " << error_I2 << ", " << error << endl;

    F.function = &Itilde3;
    gsl_integration_qags (&F, xcut, xmax, 0, 1e-7, 10000,
                        w, &result, &error);
    I3_int[i]=result;
    if (i<steps-1) {
      I3_int[i] += I3_int[i+1];
    }
    error_I3=sqrt(error_I3*error_I3+error*error);
    cout << "I3: " << setprecision(15) << setw(23) << I3_int[i] << " +- " << error_I3 << ", " << error << endl;

    F.function = &Itilde4;
    gsl_integration_qag (&F, xcut, xmax, 0, 1e-8, 10000,6,
                        w, &result, &error);
    I4_int[i]=result;
    if (i<steps-1) {
      I4_int[i] += I4_int[i+1];
    }
    error_I4=sqrt(error_I4*error_I4+error*error);
    cout << "I4: " << setprecision(15) << setw(23) << I4_int[i] << " +- " << error_I4 << ", " << error << endl;
  }
}
