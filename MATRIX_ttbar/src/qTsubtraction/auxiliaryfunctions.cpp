#include "../include/classes.cxx"
/*
#include "../include/definitions.cxx"
#include <cmath>
#include "gsl/gsl_sf_dilog.h"

#include "header/auxiliaryfunctions.h"
#include "header/specialfunctions.h"

#include <iostream>
*/
using namespace std;

// takes the finite part of the one-loop matrix element and computes the delta coefficient of C1qq
double C1qqdelta(double A_F) {
  // see notes, Eq. (15)
  double CF=4.0/3;
  return 0.5*(CF*M_PI*M_PI/6 + 0.5*A_F);
}

double C2qqreg(double z, double nf) {
  const double CA=3.0;
  const double CF=4.0/3;
  const double Z3=1.20205690316;
  const double z2=M_PI*M_PI/6;

  return
     (CF*(-344+24*M_PI*M_PI+974*z-1600*CA*z+2052*CF*z+148*nf*z-60*M_PI*M_PI*z+
     108*CF*M_PI*M_PI*z-1188*pow(z,2)+1584*CA*pow(z,2)-4104*CF*pow(z,2)-72*nf*pow(z,2)+
     72*M_PI*M_PI*pow(z,2)-216*CF*M_PI*M_PI*pow(z,2)+830*pow(z,3)+16*CA*pow(z,3)+2052*CF*pow(z,3)-
     76*nf*pow(z,3)-60*M_PI*M_PI*pow(z,3)+108*CF*M_PI*M_PI*pow(z,3)-
     272*pow(z,4)+24*M_PI*M_PI*pow(z,4)+
     324*CA*z*z2-1728*CF*z*z2 - 648*CA*pow(z,2)*z2 + 3456*CF*pow(z,2)*z2+
     324*CA*pow(z,3)*z2-1728*CF*pow(z,3)*z2 + 1188*CA*z*Z3 + 864*CF*z*Z3-
     324*CA*pow(z,3)*Z3-864*CF*pow(z,3)*Z3 - 108*CA*pow(z,2)*log(1-z) +
     108*CF*pow(z,2)*log(1-z)+108*CA*pow(z,3)*log(1-z)-
     108*CF*pow(z,3)*log(1-z)-
     216*CA*z*z2*log(1-z) + 216*CF*z*z2*log(1-z) -
     216*CA*pow(z,3)*z2*log(1-z)+216*CF*pow(z,3)*z2*log(1-z)-252*z*log(z)+
     348*CA*z*log(z)-540*CF*z*log(z)-
     60*nf*z*log(z)+612*pow(z,2)*log(z)-
     432*CA*pow(z,2)*log(z)+1404*CF*pow(z,2)*log(z)-744*pow(z,3)*log(z)+
     996*CA*pow(z,3)*log(z)-1728*CF*pow(z,3)*log(z)-60*nf*pow(z,3)*log(z)+
     384*pow(z,4)*log(z)-144*log(1-z)*log(z)+360*z*log(1-z)*log(z)-
     216*CA*z*log(1-z)*log(z) + 648*CF*z*log(1-z)*log(z) -
     432*pow(z,2)*log(1-z)*log(z) + 432*CA*pow(z,2)*log(1-z)*log(z) -
     1296*CF*pow(z,2)*log(1-z)*log(z) + 360*pow(z,3)*log(1-z)*log(z) -
     216*CA*pow(z,3)*log(1-z)*log(z) + 648*CF*pow(z,3)*log(1-z)*log(z) -
     144*pow(z,4)*log(1-z)*log(z)+216*CA*z*pow(log(1-z),2)*log(z) -
     324*CF*z*pow(log(1-z),2)*log(z)+216*CA*pow(z,3)*pow(log(1-z),2)*log(z) -
     324*CF*pow(z,3)*pow(log(1-z),2)*log(z)+27*z*pow(log(z),2)+
     99*CA*z*pow(log(z),2)-
     162*CF*z*pow(log(z),2)-18*nf*z*pow(log(z),2)+108*CA*pow(z,2)*pow(log(z),2) -
     108*CF*pow(z,2)*pow(log(z),2)+45*pow(z,3)*pow(log(z),2)-9*CA*pow(z,3)*pow(log(z),2)+
     108*CF*pow(z,3)*pow(log(z),2)-18*nf*pow(z,3)*pow(log(z),2)-72*pow(z,4)*pow(log(z),2)-
     108*CF*z*log(1-z)*pow(log(z),2)-108*CF*pow(z,3)*log(1-z)*pow(log(z),2)-
     18*z*pow(log(z),3) + 18*CA*z*pow(log(z),3)-
     18*CF*z*pow(log(z),3)+18*pow(z,3)*pow(log(z),3) +
     18*CA*pow(z,3)*pow(log(z),3)+18*CF*pow(z,3)*pow(log(z),3) -
     72*(pow(-1 + z,2)*(2 + (-1+3*CA-6*CF)*z + 2*pow(z,2)) -
     3*(CA-CF)*z*(1 + pow(z,2))*log(1-z)-3*(CA-3*CF)*z*(1+pow(z,2))*log(z))*
     gsl_sf_dilog(z) + 216*CA*z*myli3(1-z)-216*CF*z*myli3(1-z) +
     216*CA*pow(z,3)*myli3(1-z)-216*CF*pow(z,3)*myli3(1-z) -
     432*CA*z*myli3(z) + 1080*CF*z*myli3(z) -
     432*CA*pow(z,3)*myli3(z)+1080*CF*pow(z,3)*myli3(z)-1944*CF*z*Z3-
     216*CF*pow(z,3)*Z3))/(864*(-1 + z)*z)
     -1.0/4*CF*(1-z)*CF*(M_PI*M_PI/2-4); // correction due to different resummation scheme
}

double C2qg(double z) {
  const double CA=3.0;
  const double CF=4.0/3;
  const double Z3=1.20205690316;

  return (-35*CA)/48. - (13*CF)/32. +
        (43*CA)/(108.*z) + (43*CA*z)/48. +
        (43*CF*z)/32. + (CA*M_PI*M_PI*z)/24. + (CF*M_PI*M_PI*z)/12. -
        (149*CA*z*z)/216. - (5*CF*z*z)/4. - (CF*M_PI*M_PI*z*z)/12. -
        (3*CA*Z3)/8. + CF*Z3 - 2*CF*z*Z3 - (3*CA*z*z*Z3)/4. +
        2*CF*z*z*Z3 - (CA*M_PI*M_PI*log(1-z))/24.
        + (3*CA*z*log(1-z))/16. -
        (3*CF*z*log(1-z))/16. - (CA*M_PI*M_PI*z*log(1-z))/12. -
        (CA*z*z*log(1-z))/4. + (CF*z*z*log(1-z))/4. -
        (CA*M_PI*M_PI*z*z*log(1-z))/12. - (CA*z*pow(log(1-z),2))/8. +
        (CF*z*pow(log(1-z),2))/8. + (CA*z*z*pow(log(1-z),2))/8. -
        (CF*z*z*pow(log(1-z),2))/8. - (CA*pow(log(1-z),3))/48. -
        (CF*pow(log(1-z),3))/48. - (CA*z*pow(log(1-z),3))/8. +
        (CF*z*pow(log(1-z),3))/24. - (CA*z*z*pow(log(1-z),3))/24. -
        (CF*z*z*pow(log(1-z),3))/24. +(7*CA*log(z))/24.+(CF*log(z))/4.+
        (CF*M_PI*M_PI*log(z))/48.-(5*CA*z*log(z))/12.+(15*CF*z*log(z))/32.-
        (CA*M_PI*M_PI*z*log(z))/12. - (CF*M_PI*M_PI*z*log(z))/24. +
        (17*CA*z*z*log(z))/18. - (CF*z*z*log(z))/4. +
        (CF*M_PI*M_PI*z*z*log(z))/24. - (CF*z*log(1-z)*log(z))/4. +
        (CF*z*z*log(1-z)*log(z))/4. + (CF*pow(log(1-z),2)*log(z))/16.-
        (CF*z*pow(log(1-z),2)*log(z))/8.
        + (CF*z*z*pow(log(1-z),2)*log(z))/8. -
        (CA*pow(log(z),2))/32. + (CF*pow(log(z),2))/64.+(CA*z*pow(log(z),2))/8.+
        (3*CF*z*pow(log(z),2))/16. - (11*CA*z*z*pow(log(z),2))/24. -
        (CF*z*z*pow(log(z),2))/8. - (CF*log(1-z)*pow(log(z),2))/16. +
        (CA*z*log(1-z)*pow(log(z),2))/2.+(CF*z*log(1-z)*pow(log(z),2))/8.-
        (CF*z*z*log(1-z)*pow(log(z),2))/8. + (CA*pow(log(z),3))/48. -
        (CF*pow(log(z),3))/96. + (CA*z*pow(log(z),3))/24.
        + (CF*z*pow(log(z),3))/48.-
        (CF*z*z*pow(log(z),3))/24. + (CA*M_PI*M_PI*log(1 + z))/16. +
        (CA*M_PI*M_PI*z*log(1 + z))/8. + (CA*M_PI*M_PI*z*z*log(1 + z))/8. +
        (CA*pow(log(1-z),2)*log(1 + z))/8. +
        (CA*z*pow(log(1-z),2)*log(1 + z))/4. +
        (CA*z*z*pow(log(1-z),2)*log(1 + z))/4. +
        (CA*z*log(z)*log(1 + z))/4. + (CA*z*z*log(z)*log(1 +z))/4. +
        (CA*pow(log(z),2)*log(1 + z))/16. + (CA*z*pow(log(z),2)*log(1+z))/8.+
        (CA*z*z*pow(log(z),2)*log(1 + z))/8. -
        (CA*log(1-z)*pow(log(1+z),2))/8. -
        (CA*z*log(1-z)*pow(log(1+z),2))/4. -
        (CA*z*z*log(1-z)*pow(log(1+z),2))/4. + (CA*gsl_sf_dilog(1-z))/4-
        (CA*gsl_sf_dilog(1-z))/(6.*z) - CA*z*gsl_sf_dilog(1-z) +
        (11*CA*z*z*gsl_sf_dilog(1-z))/12.
         - (CA*log(1-z)*gsl_sf_dilog(1-z))/8. +
        (CF*log(1-z)*gsl_sf_dilog(1-z))/8. +
        (CA*z*log(1-z)*gsl_sf_dilog(1-z))/4. -
        (CF*z*log(1-z)*gsl_sf_dilog(1-z))/4. -
        (CA*z*z*log(1-z)*gsl_sf_dilog(1-z))/4. +
        (CF*z*z*log(1-z)*gsl_sf_dilog(1-z))/4.
         - (CF*log(z)*gsl_sf_dilog(1-z))/8. +
        (CA*z*log(z)*gsl_sf_dilog(1-z))/2. + (CF*z*log(z)*gsl_sf_dilog(1-z))/4-
        (CF*z*z*log(z)*gsl_sf_dilog(1-z))/4. + (CA*z*gsl_sf_dilog(-z))/4. +
        (CA*z*z*gsl_sf_dilog(-z))/4. - (CA*log(z)*gsl_sf_dilog(-z))/8. -
        (CA*z*log(z)*gsl_sf_dilog(-z))/4. - (CA*z*z*log(z)*gsl_sf_dilog(-z))/4. +
        (CA*myli3(1-z))/8. - (CF*myli3(1-z))/8. -
        (CA*z*myli3(1-z))/4. + (CF*z*myli3(1-z))/4. +
        (CA*z*z*myli3(1-z))/4. - (CF*z*z*myli3(1-z))/4 +
        (3*CA*myli3(-z))/8. + (3*CA*z*myli3(-z))/4. +
        (3*CA*z*z*myli3(-z))/4. - (CF*myli3(z))/8. + CA*z*myli3(z)+
        (CF*z*myli3(z))/4. - (CF*z*z*myli3(z))/4. +
        (CA*myli3(1.0/(1 + z)))/4. + (CA*z*myli3(1.0/(1 + z)))/2 +
        (CA*z*z*myli3(1.0/(1 + z)))/2.
         - (CA*myli3((-1 + z)/(1 + z)))/4. -
        (CA*z*myli3((-1 + z)/(1 + z)))/2. -
        (CA*z*z*myli3((-1 + z)/(1 + z)))/2. +
        (CA*myli3((1 + z)/(-1 + z)))/4. +
        (CA*z*myli3((1 + z)/(-1 + z)))/2. +
        (CA*z*z*myli3((1 + z)/(-1 + z)))/2.
        -1.0/4*z*(1-z)*CF*(M_PI*M_PI/2-4); // correction due to different resummation scheme
}

double C2qqb(double z, double nf) {
  const double CF=4.0/3;
  const double CA=3;
  const double Z3=1.20205690316;

  return C2qqp(z,nf)+
       (CF*(-CA+2*CF)*(45-3*M_PI*M_PI-2*M_PI*M_PI*z-45*z*z+M_PI*M_PI*z*z+9*
       log(z)+42*z*log(z)+33*z*z*log(z)+12*log(1-z)*log(z)-
       12*z*z*log(1-z)*log(z)-pow(log(z),3)-z*z*
       pow(log(z),3)+2*M_PI*M_PI*log(1+z) +
       2*M_PI*M_PI*z*z*log(1+z)-12*log(z)*
       log(1+z)-24*z*log(z)*log(1+z)-
       12*z*z*log(z)*log(1+z)+6*pow(log(z),2)*log(1+z)+
       6*z*z*pow(log(z),2)*log(1+z)-4*pow(log(1+z),3)-4*z*z*pow(log(1+z),3)-
       12*(pow(1+z,2)+(1+z*z)*log(z))*gsl_sf_dilog(-z)-
       12*(-1+z*z+log(z)+z*z*log(z))*gsl_sf_dilog(z)+36*myli3(-z)+
       36*z*z*myli3(-z)+24*myli3(z)+24*z*z*myli3(z)+
       24*myli3(1.0/(1+z))+24*z*z*myli3(1.0/(1+z))-
       18*Z3-18*z*z*Z3))/(48*(1+z));
}

//      double precision function C2qqb(z)
//      implicit none
//      real *8 Pi,Z3,myli2,myli3,z,CF,CA,C2qqp
//
//
//      external myli2,myli3,C2qqp
//
//      Pi=3.14159265358979d0
//      Z3=1.20205690316d0
//
//      CF=4d0/3
//      CA=3d0
//
//
//       C2qqb=C2qqp(z)+
//     &  (CF*(-CA+2*CF)*(45-3*Pi**2-2*Pi**2*z-45*z**2+Pi**2*z**2+9*
//     &  dlog(z)+42*z*dlog(z)+33*z**2*dlog(z)+12*dlog(1-z)*dlog(z)-
//     &  12*z**2*dlog(1-z)*dlog(z)-dlog(z)**3-z**2*
//     &  dlog(z)**3+2*Pi**2*dlog(1+z) +
//     &  2*Pi**2*z**2*dlog(1+z)-12*dlog(z)*
//     &  dlog(1+z)-24*z*dlog(z)*dlog(1+z)-
//     &  12*z**2*dlog(z)*dlog(1+z)+6*dlog(z)**2*dlog(1+z)+
//     &  6*z**2*dlog(z)**2*dlog(1+z)-4*dlog(1+z)**3-4*z**2*dlog(1+z)**3-
//     &  12*((1+z)**2+(1+z**2)*dlog(z))*myli2(-z)-
//     &  12*(-1+z**2+dlog(z)+z**2*dlog(z))*myli2(z)+36*myli3(-z)+
//     &  36*z**2*myli3(-z)+24*myli3(z)+24*z**2*myli3(z)+
//     &  24*myli3(1d0/(1+z))+24*z**2*myli3(1d0/(1+z))-
//     &  18*Z3-18*z**2*Z3))/(48*(1+z))
//
//
//      return
//      end

double C2qqp(double z, double nf) {
  const double CF=4.0/3;

  return (CF*(2*(-1+z)*(-172+143*z-136*z*z+6*M_PI*M_PI*(2-z+2*z*z))-
          12*(z*(-21 + 30*z - 32*z*z)+
         6*(-2+3*z-3*z*z+2*z*z*z)*log(1-z))*
         log(z)-9*z*(3+3*z+8*z*z)*pow(log(z),2)+18*z*(1 + z)*pow(log(z),3)-
         72*(-2 + 3*z - 3*z*z + 2*z*z*z)*gsl_sf_dilog(z)))/(864.0*z);
}

//double precision function C2qqp(z)
//      implicit none
//      real *8 Pi,myli2,z,CF
//      integer nf
//
//      external myli2
//
//      Pi=3.14159265358979d0
//
//
//      CF=4d0/3
//
//
//      C2qqp=(CF*(2*(-1+z)*(-172+143*z-136*z**2+6*Pi**2*(2-z+2*z**2))-
//     &   12*(z*(-21 + 30*z - 32*z**2)+
//     &  6*(-2+3*z-3*z**2+2*z**3)*dlog(1-z))*
//     &  dlog(z)-9*z*(3+3*z+8*z**2)*dlog(z)**2+18*z*(1 + z)*dlog(z)**3-
//     &  72*(-2 + 3*z - 3*z**2 + 2*z**3)*myli2(z)))/(864d0*z)
//
//
//      return
//      end

//CC Integral of 1/(1-x) from 0 to z
double D0int(double z) {
    return -log(1-z);
}

//CC Integral of log(1-x)/(1-x) from 0 to z
double D1int(double z) {
  return -0.5*log(1-z)*log(1-z);
//  D1int=-0.5d0*dlog(1-z)**2
}

////CC gg splitting function: regular part (with asopi normalization)
double Pggreg(double z) {
    return 3*((1-2*z)/z+z*(1-z));
}

////CC gq splitting function (with asopi normalization)
double Pgq(double z) {
  return 2.0/3*(1+(1-z)*(1-z))/z;
}

//CC qq splitting function (with asopi normalization)
double Pqq(double z) {
  return 2.0/3*(1+z*z)/(1-z);
}


//CC qg splitting function (with asopi normalization)
double Pqg(double z) {
  return 0.25*(1-2*z*(1-z));
}

//CC Non delta term in Cqq coefficient (with asopi normalization)
double Cqq(double z) {
  return 2.0/3*(1-z);
}

//CC Cqg coefficient (with asopi normalization)
double Cqg(double z) {
  return 0.5*z*(1-z);
}

// Cgq
double Cgq(double z)
{
  return 2.0/3*z;
}

//CC Integral of Pqq=1/2 CF (1+x^2)/(1-x) from 0 to z
double Pqqint(double z) {
  return -2.0/3*(z+z*z/2+2*log(1-z));
}

//C                P*P convolutions

//CC Regular part of Pqq*Pqq (checked !)
double Pqqqq(double z) {
  return 4.0/9*(-4.0*log(z)/(1-z)-2.0*(1-z)+(1+z)*(3.0*log(z)-4.0*log(1-z)-3));
//  Pqqqq=4d0/9*(-4*dlog(z)/(1-z)-2*(1-z)
//     &  +(1+z)*(3*dlog(z)-4*dlog(1-z)-3))
}

//CC Pqq*Pqg (checked !)
double Pqqqg(double z) {
  return 1.0/3*((z*z+(1-z)*(1-z))*log((1-z)/z)-(z-0.5)*log(z)+z-0.25);
//        Pqqqg=1d0/3*((z**2+(1-z)**2)*dlog((1-z)/z)     &  -(z-0.5d0)*dlog(z)+z-0.25d0)
}

//CC Pqg*Pgq (checked !)
double Pqggq(double z) {
  return 1.0/3*(2.0/3/z+(1+z)*log(z)-2.0/3*z*z-0.5*(z-1));
//  Pqggq=1d0/3*(2d0/3/z+(1+z)*dlog(z)-2d0/3*z**2-0.5d0*(z-1))
}

//CC Full Pqg*Pgg (checked !)
double Pqggg(double z, double beta0) {
  return 1.5*(1.0/3/z+(z*z-z+0.5)*log(1-z)+(2*z+0.5)*log(z)+0.25+2*z-31.0/12*z*z)+beta0*Pqg(z);

//  beta0=(33-2*nf)/12d0
//      Pqggg=1.5d0*(1/3d0/z+(z**2-z+0.5d0)*dlog(1-z)
//     &     +(2*z+0.5d0)*dlog(z)+0.25d0+2*z-31d0/12*z**2)
//
//      Pqggg=Pqggg+beta0*Pqg(z)
}

//C           Two loop AP:  pqq of ESW is my 3/2 Pqq

//C     Pqq NS: Eq. (4.107) ESW (no 1/(1-x)_+ and delta term)



//----------------------------------------------------------------------------------------------
// new
double Pgqqg( double z, int nf )
{
  // Pgqqg = nf/3*( 1 + 4d0/3/x - x - 4*x**2 /3 + 2*(1+x)*dlog(x) )
  return   nf/3.0*( 1 + 4.0/3/z - z - 4.0*z*z/3 + 2*(1+z)* log(z) );
}

double Pggggreg( double z, double beta0)
{
  //  Pggggreg = -9*dlog(x)/(1-x) + 2*beta0*Pggreg(x) + 9* ( 3*(1-x) + 11d0/3/x*(x**3 -1d0) + 2d0/3*dlog(1-x)*Pggreg(x) + dlog(x)*(x**2-3*x-1d0/x) )
  return         -9* log(z)/(1-z) + 2*beta0*Pggreg(z) + 9* ( 3*(1-z) + 11.0/3/z*(z*z*z-1.0) + 2.0/3* log(1-z)*Pggreg(z) +  log(z)*(z*z -3*z-1.0/z) ) ;
}

double Pgqqq( double z )
{
  //  Pgqqq = 4d0/9*( (2-x)*dlog(x) + dlog(1-x)*(2*x+4/x  -4) + 2 - x/2.0 )
  return      4.0/9*( (2-z)* log(z) +  log(1-z)*(2*z+4.0/z-4) + 2 - z/2.0 );
}

double Pgggq( double z, double beta0 )
{
  //  Pgggq = 2 * ( (1+(1-x)**2)   /x*dlog(1-x) - 2*(1+x+1d0/x)*dlog(x) + 4 - 31d0/6/x + x/2   + 2d0/3*x**2 ) + beta0*Pgq(x)
  return      2 * ( (1+(1-z)*(1-z))/z* log(1-z) - 2*(1+z+1.0/z)* log(z) + 4 - 31.0/6/z + z/2.0 + 2.0/3*z*z  ) + beta0*Pgq(z);
}

double CgqPqg( double z, int nf )
{
  // CgqPqg = 2*nf*1d0/6*( 1 + x - 2*x**2 + 2*x*dlog(x) )
  return      nf/3.0    *( 1 + z - 2*z*z  + 2*z* log(z) );
}

double CgqPqq( double z )
{
  // CgqPqq = 2d0/9*( 2 + x + 4*x*dlog(1-x) - 2*x*dlog(x) )
  return      2.0/9*( 2 + z + 4*z* log(1-z) - 2*z* log(z) );
}

double Pgg( double z )
{
  // Pgg = 1d0/(1-x)+1d0/x-2+x*(1-x)
  return   1.0/(1-z)+1.0/z-2+z*(1-z);
}

double P2gg( double z, int nf )
{
  double P2gg;

//P2gg = 2d0/3*nf* [ -16 + 8*x + 20d0/3*x**2 + 4d0/3/x - (6+10*x)*dlog(x) - (2+2*x)*dlog(x)**2     ] + 1.5d0*nf* [ 2 - 2*x + 26d0/9*(x**2-1d0/x) - 4d0/3*(1+x)*dlog(x) - 20d0/9*(1/x-2+x*(1-x)) ]
  P2gg = 2.0/3*nf* ( -16 + 8*z + 20.0/3*z*z  + 4.0/3/z - (6+10*z)* log(z) - (2+2*z)* log(z)*log(z) ) + 1.5  *nf* ( 2 - 2*z + 26.0/9*(z*z -1.0/z) - 4.0/3*(1+z)* log(z) - 20.0/9*(1/z-2+z*(1-z)) );

//            + 9*[ 27d0/2*(1-x) + 67d0/9*(x**2-1d0/x) - (25d0/3-11d0/3*x+44d0/3*x**2)*dlog(x) + 4*(1+x)*dlog(x)**2     + 2*pgg(-x)*S2(x) + (67d0/9-pi**2    /3)*(1/x-2+x*(1-x)) + (-4*dlog(x)*dlog(1-x)+dlog(x)**2    )*pgg(x) ]    
  P2gg = P2gg + 9*( 27.0/2*(1-z) + 67.0/9*(z*z -1.0/z) - (25.0/3-11.0/3*z+44.0/3*z*z )* log(z) + 4*(1+z)* log(z)*log(z) + 2*Pgg(-z)*S2(z) + (67.0/9-M_PI*M_PI/3)*(1/z-2+z*(1-z)) + (-4* log(z)* log(1-z)+ log(z)*log(z))*Pgg(z) );

//P2gg = P2gg/4d0
  P2gg = P2gg/4.0;

  return P2gg;
}

double P2gq( double z, int nf )
{
  double P2gq;  // logx=dlog(x)  logomx=dlog(1-x)

//P2gq = 16d0/9*( -2.5d0 - 3.5d0*x + (2+3.5d0*x)*logx   - (1-0.5d0*x)*logx**2       - 2*x*logomx   - (3*logomx  +logomx**2        )*1.5d0*Pgq(x) )
  P2gq = 16.0/9*( -2.5   - 3.5  *z + (2+3.5  *z)*log(z) - (1-0.5  *z)*log(z)*log(z) - 2*z*log(1-z) - (3*log(1-z)+log(1-z)*log(1-z))*1.5  *Pgq(z) );

//         +4d0*(28d0/9+65d0/18*x+44d0/9*x**2-(12+5*x+8d0/3*x**2)*logx  +(4+x)*logx**2      +2*x*logomx  +S2(x)*1.5d0*Pgq(-x)+(0.5d0-2*logx  *logomx  +0.5d0*logx**2      +11d0/3*logomx  +logomx**2        -pi**2/6    )*1.5d0*Pgq(x))
  P2gq=P2gq+4.0*(28.0/9+65.0/18*z+44.0/9*z*z -(12+5*z+8.0/3*z*z )*log(z)+(4+z)*log(z)*log(z)+2*z*log(1-z)+S2(z)*1.5  *Pgq(-z)+(0.5  -2*log(z)*log(1-z)+0.5  *log(z)*log(z)+11.0/3*log(1-z)+log(1-z)*log(1-z)-M_PI*M_PI/6)*1.5  *Pgq(z));

//            + 2d0/3*nf*( -4d0*x/3 - (20d0/9+4d0/3*logomx  )*1.5d0*Pgq(x))
  P2gq = P2gq + 2.0/3*nf*( -4.0*z/3 - (20.0/9+4.0/3*log(1-z))*1.5  *Pgq(z));

//P2gq = P2gq/4
  P2gq = P2gq/4;
  
  return P2gq;
}

double C1ggdelta( double A_F )
{
  double CA=3.0;
  return 0.5*(CA*M_PI*M_PI/6 + 0.5*A_F);
}

 
double C2ggreg( double z, int nf )
{
  double CF = 4.0/3, CA = 3.0;
  
  return 
    ( CA*CA*( pow(1+z+z*z,2)/(z*(z+1))*(2*myli3(z/(1+z))-myli3(-z)) + (2-17*z-22*z*z-10*z*z*z-12*z*z*z*z)*zeta3/(2*z*(1+z)) - (5-z+5*pow(z,2)+pow(z,3)-5*pow(z,4)+pow(z,5))/(z*(1-z)*(1+z))*(myli3(z)-zeta3)
	      + gsl_sf_dilog(z)*log(z)/(1-z)*(3-z+3*pow(z,2)+pow(z,3)-3*pow(z,4)+pow(z,5))/(z*(1+z)) + pow(1+z+z*z,2)/(z*(1+z))*(log(z)*gsl_sf_dilog(-z)-pow(log(1+z),3)/3.0+zeta2*log(1+z)) 
	      + (1-z)/(3*z)*(11-z+11*z*z)*gsl_sf_dilog(1-z)+1.0/12*z*log(1-z) - 1.0/6*pow(log(z),3)/(1-z)*pow(1+z-z*z,2)/(1+z)+pow(log(z),2)*(pow(1-z+z*z,2)/(2*z*(1-z))*log(1-z)-pow(1+z+z*z,2)/(2*z*(1+z))*log(1+z)+(25-11*z+44*z*z)/24.0)
	      + log(z)*( pow(1+z+z*z,2)/(z*(1+z))*pow(log(1+z),2)+pow(1-z+z*z,2)/(2*z*(1-z))*pow(log(1-z),2)-(72+773*z+149*z*z+536*z*z*z)/(72.0*z)) + 517.0/27 - 449.0/(27*z) -380.0*z/27 + 835.0*z*z/54 ) 
      + CA*nf*( (1+z)/12.0*pow(log(z),2) + 1.0/36*(13+10*z)*log(z) - z/12.0*log(1-z) -83.0/54 + 121.0/(108*z) + 55.0/54*z - 139.0/108*z*z ) 
      + CF*nf*( (1+z)/12.0*pow(log(z),3) + 1.0/8*(3+z)*pow(log(z),2) +3.0/2*(1+z)*log(z) - (1-z)/(6*z)*(1-23*z+z*z) ) + CA*CA*((1+z)*log(z)/z+2*(1-z)/z) )/2.0 ;

  /*
  return 
    -(-10776 + 226*nf + 396*M_PI*M_PI + 12408*z + 52*nf*z - 432*M_PI*M_PI*z + 
      1656*pow(z,2) - 390*nf*pow(z,2) + 36*M_PI*M_PI*pow(z,2) - 2388*pow(z,3) - 
      314*nf*pow(z,3) + 36*M_PI*M_PI*pow(z,3) + 9120*pow(z,4) + 164*nf*pow(z,4) - 
      432*M_PI*M_PI*pow(z,4) - 10020*pow(z,5) + 262*nf*pow(z,5) + 396*M_PI*M_PI*pow(z,5) + 
      54*pow(z,2)*log(1-z) - 18*nf*pow(z,2)*log(1-z) - 54*pow(z,4)*log(1-z) + 
      18*nf*pow(z,4)*log(1-z) - 648*log(z) - 6957*z*log(z) + 
      222*nf*z*log(z) - 693*pow(z,2)*log(z) + 204*nf*pow(z,2)*log(z) + 
      2133*pow(z,3)*log(z) - 222*nf*pow(z,3)*log(z) + 1341*pow(z,4)*log(z) - 
      204*nf*pow(z,4)*log(z) + 4824*pow(z,5)*log(z) - 2376*log(1-z)*log(z) + 
      2592*z*log(1-z)*log(z) - 216*pow(z,2)*log(1-z)*log(z) - 
      216*pow(z,3)*log(1-z)*log(z) + 2592*pow(z,4)*log(1-z)*log(z) - 
      2376*pow(z,5)*log(1-z)*log(z) + 324*log(1-z)*log(1-z)*log(z) - 
      324*z*log(1-z)*log(1-z)*log(z) + 324*pow(z,2)*log(1-z)*log(1-z)*log(z) + 
      324*pow(z,3)*log(1-z)*log(1-z)*log(z) - 324*pow(z,4)*log(1-z)*log(1-z)*log(z) + 
      324*pow(z,5)*log(1-z)*log(1-z)*log(z) + 675*z*log(z)*log(z) + 
      54*nf*z*log(z)*log(z) - 297*pow(z,2)*log(z)*log(z) + 30*nf*z*z*log(z)*log(z) + 
      513*pow(z,3)*log(z)*log(z) - 54*nf*pow(z,3)*log(z)*log(z) + 297*pow(z,4)*log(z)*log(z) - 
      30*nf*pow(z,4)*log(z)*log(z) - 1188*pow(z,5)*log(z)*log(z) + 
      324*log(1-z)*log(z)*log(z) - 324*z*log(1-z)*log(z)*log(z) + 
      324*pow(z,2)*log(1-z)*log(z)*log(z) + 324*pow(z,3)*log(1-z)*log(z)*log(z) - 
      324*pow(z,4)*log(1-z)*log(z)*log(z) + 324*pow(z,5)*log(1-z)*log(z)*log(z) - 
      108*z*pow(log(z),3) + 8*nf*z*pow(log(z),3) - 216*pow(z,2)*pow(log(z),3) + 
      8*nf*pow(z,2)*pow(log(z),3) + 108*pow(z*log(z),3) - 8*nf*pow(z*log(z),3) + 
      216*pow(z,4)*pow(log(z),3) - 8*nf*pow(z,4)*pow(log(z),3) - 108*pow(z,5)*pow(log(z),3) + 
      108*M_PI*M_PI*log(1+z) + 108*M_PI*M_PI*z*log(1+z) + 
      108*M_PI*M_PI*pow(z,2)*log(1+z) - 108*M_PI*M_PI*pow(z,3)*log(1+z) - 
      108*M_PI*M_PI*pow(z,4)*log(1+z) - 108*M_PI*M_PI*pow(z,5)*log(1+z) - 
      324*pow(log(z),2)*log(1+z) - 324*z*pow(log(z),2)*log(1+z) - 
      324*pow(log(z)*z,2)*log(1+z) + 324*pow(z,3)*pow(log(z),2)*log(1+z) + 
      324*pow(z,4)*pow(log(z),2)*log(1+z) + 324*pow(z,5)*pow(log(z),2)*log(1+z) + 
      648*log(z)*pow(log(1+z),2) + 648*z*log(z)*pow(log(1+z),2) + 
      648*pow(z,2)*log(z)*pow(log(1+z),2) - 648*pow(z,3)*log(z)*pow(log(1+z),2) - 
      648*pow(z,4)*log(z)*pow(log(1+z),2) - 648*pow(z,5)*log(z)*pow(log(1+z),2) - 
      216*pow(log(1+z),3) - 216*z*pow(log(1+z),3) - 216*pow(z,2)*pow(log(1+z),3) + 
      216*pow(z,3)*pow(log(1+z),3) + 216*pow(z,4)*pow(log(1+z),3) + 
      216*pow(z,5)*pow(log(1+z),3) - 
      648*(-1 + z)*pow((1 + z + z*z),2)*log(z)*gsl_sf_dilog(-z) + 
      216*(-(pow((-1 + z),2)*(11 + 10*z + 10*z*z + 11*z*z*z)) + 
	   3*(3 - z + 3*pow(z,2) + pow(z,3) - 3*pow(z,4) + pow(z,5))*log(z))*gsl_sf_dilog(z) -
      648*myli3(-z) - 648*z*myli3(-z) - 
      648*pow(z,2)*myli3(-z) + 648*pow(z,3)*myli3(-z) + 
      648*pow(z,4)*myli3(-z) + 648*pow(z,5)*myli3(-z) - 
      3240*myli3(z) + 648*z*myli3(z) - 3240*pow(z,2)*myli3(z) - 
      648*pow(z,3)*myli3(z) + 3240*pow(z,4)*myli3(z) - 
      648*pow(z,5)*myli3(z) + 1296*myli3(z/(1 + z)) + 
      1296*z*myli3(z/(1 + z)) + 1296*pow(z,2)*myli3(z/(1 + z)) - 
      1296*pow(z,3)*myli3(z/(1 + z)) - 1296*pow(z,4)*myli3(z/(1 + z)) - 
      1296*pow(z,5)*myli3(z/(1 + z)) + 3888*zeta3 - 6804*z*zeta3 + 
      1620*pow(z,2)*zeta3 + 4536*pow(z,3)*zeta3 - 3888*pow(z,4)*zeta3 + 
      4536*pow(z,5)*zeta3)/(72.0*z*(-1 + pow(z,2))) -
      (-9)*(2*(1-z)+(1+z)*log(z))/z;*/
}


double C2gq( double z, int nf )
{
  double CF = 4.0/3, CA = 3.0;

  return 
    CF*CF*( (2-z)*pow(log(z),3)/48.0 - (3*z+4)*pow(log(z),2)/32.0 + 5*(z-3)*log(z)/16.0 + 1.0/12*(1.0/z+z/2.0-1)*pow(log(1-z),3) + 1.0/16*(z+6.0/z-6)*pow(log(1-z),2) + (5*z/8.0+2.0/z-2)*log(1-z) + 5.0/8 - 13.0/16*z )
    + CF*nf*( 1/(24.0*z)*(z*z-2*z+2)*pow(log(1-z),2) + 1/18.0*(z+5.0/z-5)*log(1-z) - 14.0/27 + 14.0/(27*z) + 13*z/108.0 )
    + CF*CA*( -(z*z+2*z+2)/(2*z)*myli3(1.0/(1+z)) + (1.0/2-5.0/(2*z)-5*z/4.0)*myli3(z) - 3.0/(4*z)*(z*z+2*z+2)*myli3(-z) + (2-11.0/(6*z)-z/2.0+z*z/3.0+(-1.0/2+3.0/(2*z)+3*z/4.0)*log(z))*gsl_sf_dilog(z)
	      + (z/4.0+(z*z+2*z+2)/(4.0*z)*log(z))*gsl_sf_dilog(-z) + (z*z+2*z+2)/(12.0*z)*pow(log(1+z),3) - 1.0/(24*z)*((z*z+2*z+2)*(3*pow(log(z),2)+M_PI*M_PI)-6*z*z*log(z))*log(1+z)
	      - (z*z-2*z+2)/(24.0*z)*pow(log(1-z),3) + 1.0/(48*z)*(6*(z*z-2*z+2)*log(z) - 5*z*z - 22*(1-z))*pow(log(1-z),2) + 1.0/(72*z)*(-152+152*z-43*z*z+6*(-22+24*z-9*z*z+4*z*z*z)*log(z)+9*(z*z-2*z+2)*pow(log(z),2))*log(1-z)
	      - 1.0/12*(1+z/2.0)*pow(log(z),3) + 1.0/48*(36+9*z+8*z*z)*pow(log(z),2) + (-107.0/24-1.0/z+z/12.0-11*z*z/9)*log(z) + 1.0/z*(4*zeta3-503.0/54+11.0/36*M_PI*M_PI)+1007.0/108 - M_PI*M_PI/3.0 - 5.0/2*zeta3 
	      + z*(M_PI*M_PI/3.0+2*zeta3-133.0/108) + z*z*(38.0/27-M_PI*M_PI/18.0) )
    + CF*CF*z*3.0/4 + CF*CA/z*((1+z)*log(z)+2*(1-z)-(5.0+M_PI*M_PI)/4*z*z) ;

  
  /*
    return
      -(12072 - 224*nf - 396*M_PI*M_PI - 12444*z + 224*nf*z + 432*M_PI*M_PI*z +
      2064*pow(z,2) - 52*nf*pow(z,2) - 432*M_PI*M_PI*pow(z,2) - 1824*pow(z,3) +
      72*M_PI*M_PI*pow(z,3) + 1584*log(1-z) - 120*nf*log(1-z) -
      162*M_PI*M_PI*log(1-z) - 1584*z*log(1-z) + 120*nf*z*log(1-z) -
      324*M_PI*M_PI*z*log(1-z) + 414*pow(z,2)*log(1-z) -
      24*nf*pow(z,2)*log(1-z) - 162*M_PI*M_PI*pow(z,2)*log(1-z) +
      378*pow(log(1-z),2) - 36*nf*pow(log(1-z),2) - 378*z*pow(log(1-z),2) +
      36*nf*z*pow(log(1-z),2) + 99*pow(z*log(1-z),2) -
      18*nf*pow(z*log(1-z),2) + 60*pow(log(1-z),3) -
      60*z*pow(log(1-z),3) + 30*pow(z,2)*pow(log(1-z),3) + 1296*log(z) +
      6318*z*log(z) - 288*pow(z,2)*log(z) + 1584*pow(z,3)*log(z) +
      2376*log(1-z)*log(z) - 2592*z*log(1-z)*log(z) +
      972*pow(z,2)*log(1-z)*log(z) - 432*pow(z,3)*log(1-z)*log(z) +
      324*pow(log(1-z),2)*log(z) + 1620*z*pow(log(1-z),2)*log(z) +
      486*pow(z*log(1-z),2)*log(z) - 900*z*pow(log(z),2) -
      189*pow(z*log(z),2) - 216*pow(z,3)*pow(log(z),2) -
      324*log(1-z)*pow(log(z),2) + 324*z*log(1-z)*pow(log(z),2) -
      162*pow(z,2)*log(1-z)*pow(log(z),2) + 84*z*pow(log(z),3) +
      66*pow(z,2)*pow(log(z),3) + 270*M_PI*M_PI*log(1+z) +
      432*M_PI*M_PI*z*log(1+z) + 216*M_PI*M_PI*pow(z,2)*log(1+z) -
      324*pow(z,2)*log(z)*log(1+z) - 1296*log(1-z)*log(z)*log(1+z) -
      2592*z*log(1-z)*log(z)*log(1+z) -
      1296*pow(z,2)*log(1-z)*log(z)*log(1+z) +
      324*pow(log(z),2)*log(1+z) + 324*z*pow(log(z),2)*log(1+z) +
      162*pow(z*log(z),2)*log(1+z) + 648*log(z)*pow(log(1+z),2)+
      1296*z*log(z)*pow(log(1+z),2) + 648*pow(z,2)*log(z)*pow(log(1+z),2) -
      216*pow(log(1+z),3) - 216*z*pow(log(1+z),3) -
      108*pow(z,2)*pow(log(1+z),3) -
      324*(z*z + 2*pow((1+z),2)*log(1-z) + (2+2*z+z*z)*log(z) -
	   2*log(1+z) - 4*z*log(1+z) - 2*pow(z,2)*log(1+z))*gsl_sf_dilog(-z)-108*
      (-22 + 24*z - 6*z*z + 4*pow(z,3) - 6*pow((1+z),2)*log(1-z) +
       3*(6 - 2*z + 3*z*z)*log(z) + 6*log(1+z) + 12*z*log(1+z) +
       6*z*z*log(1+z))*gsl_sf_dilog(z) +
      648*log(1-z)*gsl_sf_dilog((1-z)/(1+z)) +
      1296*z*log(1-z)*gsl_sf_dilog((1-z)/(1+z)) +
      648*z*z*log(1-z)*gsl_sf_dilog((1-z)/(1+z)) -
      648*log(1+z)*gsl_sf_dilog((1-z)/(1+z)) -
      1296*z*log(1+z)*gsl_sf_dilog((1-z)/(1+z)) -
      648*z*z*log(1+z)*gsl_sf_dilog((1-z)/(1+z)) -
      648*log(1-z)*gsl_sf_dilog((-1+z)/(1+z)) -
      1296*z*log(1-z)*gsl_sf_dilog((-1+z)/(1+z)) -
      648*z*z*log(1-z)*gsl_sf_dilog((-1+z)/(1+z)) +
      648*log(1+z)*gsl_sf_dilog((-1+z)/(1+z)) +
      1296*z*log(1+z)*gsl_sf_dilog((-1+z)/(1+z)) +
      648*z*z*log(1+z)*gsl_sf_dilog((-1+z)/(1+z)) +
      1944*myli3(-z)+1944*z*myli3(-z) +
      972*z*z*myli3(-z)+3240*myli3(z)-648*z*myli3(z) +
      1620*z*z*myli3(z)+1296*myli3(1/(1+z)) +
      1296*z*myli3(1/(1+z)) + 648*z*z*myli3(1/(1+z)) -
      5184*zeta3 + 3240*z*zeta3 - 2592*z*z*zeta3)/(324.0*z) -
      (-4)*(2*(1-z)+(1+z)*log(z))/z;*/
}

double Ggg( double z )
{
  return 3*(1-z)/z;
  ///  Only in order to validate without Gga terms !!!
  ///  return 0.;
}

double Ggq( double z )
{
  return 4.0/3*(1-z)/z;
  ///  Only in order to validate without Gga terms !!!
  ///  return 0.;
}


//------------------------------------------------------------------------------------------------

double P2qqV(double x, int nf) {
  return 0.25*(16.0/9*(-(2*log(x)*log(1-x)+1.5*log(x))*3.0/2*Pqq(x)-(1.5+3.5*x)*log(x)-0.5*(1+x)*log(x)*log(x)-5*(1-x))+4*((0.5*log(x)*log(x)+11.0/6*log(x))*3.0/2*Pqq(x)-(67.0/18-M_PI*M_PI/6)*(1+x)+(1+x)*log(x)+20.0/3*(1-x))+2.0/3*nf*(-log(x)*Pqq(x)+10.0/9*(1+x)-4.0/3*(1-x)));

//  P2qqV=16d0/9*(-(2*dlog(x)*dlog(1-x)+1.5d0*dlog(x))*3d0/2*Pqq(x)
//     &     -(1.5d0+3.5d0*x)*dlog(x)-0.5d0*(1+x)*dlog(x)**2-5*(1-x))
//     &     +4*((0.5d0*dlog(x)**2+11d0/6*dlog(x))*3d0/2*Pqq(x)
//     &     -(67d0/18-pi**2/6)*(1+x)
//     &     +(1+x)*dlog(x)+20d0/3*(1-x))
//     &     +2d0/3d0*nf*(-dlog(x)*Pqq(x)+10d0/9*(1+x)-4d0/3*(1-x))
//
//c     Change to as/pi normalization
//
//      P2qqV=P2qqV/4
}

//C    S2: Eq. (4.114) ESW
double S2(double x) {
  return -2*gsl_sf_dilog(-x)+0.5*log(x)*log(x)-2*log(x)*log(1+x)-M_PI*M_PI/6;
//  S2=-2*myli2(-x)+0.5d0*dlog(x)**2-2*dlog(x)*dlog(1+x)-pi**2/6
}

//C    Pqqb NS: Eq. (4.108) ESW
double P2qqbV(double x) {
  return 0.25*(-2.0/9*(3*Pqq(-x)*S2(x)+2*(1+x)*log(x)+4*(1-x)));

//  P2qqbV=-2d0/9*(3d0*Pqq(-x)*S2(x)+2*(1+x)*dlog(x)+4*(1-x))
//
//c     Change to as/pi normalization
//
//      P2qqbV=P2qqbV/4
}

//C    Pqg Singlet: Eq. (4.110) ESW (ESW Pqg is 4 times my Pqg)
double P2qg(double x) {
  double logomxsx=log((1-x)/x);

  return 0.5*0.25*(2.0/3*(4-9*x-(1-4*x)*log(x)-(1-2*x)*log(x)*log(x)+4*log(1-x)+(2*logomxsx*logomxsx-4*logomxsx-2.0/3*M_PI*M_PI+10)*4*Pqg(x))+1.5*(182.0/9+14.0/9*x+40.0/9/x+(136.0/3*x-38.0/3)*log(x)-4*log(1-x)-(2+8*x)*log(x)*log(x)+8*Pqg(-x)*S2(x)+(-log(x)*log(x)+44.0/3*log(x)-2*log(1-x)*log(1-x)+4*log(1-x)+M_PI*M_PI/3-218.0/9)*4*Pqg(x)));

//  logx=dlog(x)
//      logomxsx=dlog((1-x)/x)
//
//      P2qg=2d0/3*(4-9*x-(1-4*x)*logx-(1-2*x)*logx**2+4*dlog(1-x)
//     &    +(2*logomxsx**2-4*logomxsx-2d0/3*pi**2+10d0)*4*Pqg(x))
//     &    +1.5d0*(182d0/9+14d0/9*x+40d0/9/x+(136d0/3*x-38d0/3)*logx
//     &    -4*dlog(1-x)-(2+8*x)*logx**2+8*Pqg(-x)*S2(x)
//     &    +(-logx**2+44d0/3*logx-2*dlog(1-x)**2+4*dlog(1-x)+pi**2/3
//     &    -218d0/9)*4*Pqg(x))
//
//c     Change to as/pi normalization
//
//      P2qg=P2qg/4d0
//c     Divide by 2 to eliminate 2nf factor
//
//      P2qg=P2qg/2d0
}

//C     Pqq Pure Singlet appearing in ESW Eq. (4.95)
//C     PSqq=PSqqb
//C     Obtained through Eq.(4.101)
//C     PSqq=1/2/nf (P2qq-P2qqbV-P2qqV) (contains only CF TR=2/3)
// okay
double P2qqS(double x) {
  return 0.25*(2.0/3*(20-18*x+54*x*x-56*x*x*x+3*x*(3+15*x+8*x*x)*log(x)-9*x*(1+x)*log(x)*log(x))/(9*x));
//  P2qqS=2d0/3*(20 - 18*x + 54*x**2 - 56*x**3
//     &    +3*x*(3 + 15*x + 8*x**2)*dlog(x)
//     &    - 9*x*(1 + x)*dlog(x)**2)/(9*x)
//
//      P2qqS=P2qqS/4
}

//C                C*P convolutions
// all okay
//CC Cqq*Pqq (without delta term in Cqq) (checked !)
double CqqPqq(double z) {
  return 2.0/9*(1-z)*(4*log(1-z)-2*log(z)-1);
//  CqqPqq=2d0/9*(1-z)*(4*dlog(1-z)-2*dlog(z)-1)
}

//CC Cqq*Pqg (without delta term in Cqq) (checked !)
double CqqPqg(double z) {
  return (-2+z+z*z-(1+2*z)*log(z))/6;
//  CqqPqg=(-2+z+z**2-(1+2*z)*dlog(z))/6d0
}

//CC Cqg*Pgq (checked !)
double CqgPgq(double z) {
  return (1.0/3/z-1+2*z*z/3-z*log(z))/3;
//  CqgPgq=(1d0/3/z-1+2*z**2/3-z*dlog(z))/3d0
}

//CC Cqg*Pgg (checked !)
double CqgPgg(double z, double beta0) {
  return 3.0/4*(2*z*(1-z)*log(1-z)-4*z*log(z)+1.0/3/z-1-5*z+17*z*z/3)+beta0/2*z*(1-z);
//  CqgPgg=3d0/4*(2*z*(1-z)*dlog(1-z)-4*z*dlog(z)
//     &      +1d0/3/z-1-5*z+17d0*z**2/3)+beta0/2*z*(1-z)
}
