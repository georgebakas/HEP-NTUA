#include "../include/classes.cxx"
//
// constructors
//
fourvector::fourvector(){
  t = 0;
  x = 0;
  y = 0;
  z = 0;
  //  crossing = 1;
}
fourvector::fourvector(double r0, double r1, double r2, double r3){
  t = r0;
  x = r1;
  y = r2;
  z = r3;
  //  crossing = 1;
}
fourvector::fourvector(double r0, double r1, double r2, double r3, int cv){
  t = r0;
  x = r1;
  y = r2;
  z = r3;
  //  crossing = cv;
}
fourvector::fourvector(double phi, double costheta, fourvector pz){
  t = pz.t;
  x = cos(phi) * sqrt(1. - pow(costheta, 2.)) * pz.z;
  y = sin(phi) * sqrt(1. - pow(costheta, 2.)) * pz.z;
  z = costheta * pz.z;
  //  crossing = pz.crossing;
  //  crossing = 1;
}
/*
fourvector::fourvector(bispinor_dud bs){
  t = 0.5 * real(bs.x11() + bs.x22());
  x = 0.5 * real(bs.x12() + bs.x21());
  y = 0.5 * imag(bs.x12() - bs.x21());
  z = 0.5 * real(bs.x11() - bs.x22());
  crossing = 1;
}
*/
//fourvector::~fourvector();

//
// access elements
//
double fourvector::x0() const {return t;}
double fourvector::x1() const {return x;}
double fourvector::x2() const {return y;}
double fourvector::x3() const {return z;}
//int fourvector::cr() const {return crossing;}
//
// methods
//
//
// Betrag des raeumlichen Anteils
//
double fourvector::r() const {
  //  return (sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
  return (sqrt(x * x + y * y + z * z));
}
double fourvector::r2() const {
  //  return pow(x, 2) + pow(y, 2) + pow(z, 2);
  return x * x + y * y + z * z;
}

//
// P - Transformation
//
fourvector fourvector::Pinv() const{
  fourvector tv;
  tv.t = t;
  tv.x = -x;
  tv.y = -y;
  tv.z = -z;
  //  tv.crossing = crossing;
  return tv;
}

//
// Spiegelung
//
fourvector fourvector::Prot() const{
  fourvector tv;
  tv.t = t;
  tv.x = x;
  tv.y = -y;
  tv.z = -z;
  //  tv.crossing = crossing;
  return tv;
}

//
// Crossing
//
/*
fourvector fourvector::crossed() const{
  fourvector tv;
  tv.t = t;
  tv.x = x;
  tv.y = y;
  tv.z = z;
  //  tv.crossing = -crossing;
  return tv;
}
*/

double fourvector::theta() const {
  return acos(z / r());
}
double fourvector::costheta() const {
  return z / r();
}
/*
double fourvector::sintheta() const {
  double sintheta0;
  sintheta0 = sqrt(1. - pow(z / r(), 2));
  return sintheta0;
}
*/
double fourvector::eta() const {
  /*
  double eta = -log(tan((theta() / 2.)));
  return eta;
  */
  double eps=(x*x+y*y)/z/z;
  if (eps==0) {
    if (z>=0) {
      return +numeric_limits<double>::infinity();
    } else {
      return -numeric_limits<double>::infinity();
    }
  }
  double _eta=0.5*log((r()+z)/(r()-z));
 
  if (eps<1e-13) {
    if (z>0) {
      _eta=0.5*log((r()+z)/(0.5*z*eps));
    } else {
      _eta=0.5*log(-0.5*z*eps/(r()-z));
    }
  }
  if (munich_isnan(_eta) || abs(_eta)>1.e99) {
    cout << "fourvector::eta: computation unstable; eta=" << _eta << ", eps=" << eps << endl;
    cout << (*this) << endl;
    cout << "r=" << r() << ", z=" << z << ", old implementation: eta=" << -log(tan((theta() / 2.))) << endl;
  }
  return _eta;
}

double fourvector::rapidity() const {
  double qq = z / t;
  return .5 * log((qq + 1) / (1 - qq));
}


//
// Winkel phi des raeumlichen Anteils
//
double fourvector::phi() const {
  double phi0 = atan2(y,x);
  if (phi0<0)
    phi0+=2*M_PI;

  return phi0;
}
double fourvector::pT() const {
  return sqrt(x * x + y * y);
}

double fourvector::pT2() const {
  return x * x + y * y;
}

double fourvector::ET() const {
  return sqrt(t * t - z * z);
}

double fourvector::ET2() const {
  return t * t - z * z;
}

double fourvector::cosphi() const {
  if (y != 0. || x != 0.){
    return x / sqrt(x * x + y * y);
    //    return x / sqrt(pow(x, 2) + pow (y, 2));
  }
  else {
    return 1.;
  }
}

double fourvector::sinphi() const {
  if (y != 0){
    return y / sqrt(x * x + y * y);
    //    return y / sqrt(pow(x, 2) + pow (y, 2));
  }
  else {
    return 0.;
  }
}

double fourvector::m2() const {
  //  return t * t - (x * x + y * y + z * z);
  return t * t - x * x - y * y - z * z;
}

double fourvector::m() const {
  return sqrt(t * t - (x * x + y * y + z * z));
}

fourvector fourvector::operator - () const{
  fourvector tv;
  tv.t = -t;
  tv.x = -x;
  tv.y = -y;
  tv.z = -z;
  //  tv.crossing = crossing;
  return tv;
}

//
// Rotation eines x3-orientierten Vektors zu phi / theta
//
fourvector fourvector::rot3(double phi, double costheta) const{
  fourvector tv;
  double sintheta = sqrt(1. - pow(costheta, 2.));
  tv.t = t;
  tv.x = cos(phi) * sintheta * z;
  tv.y = sin(phi) * sintheta * z;
  tv.z = costheta * z;
  //  tv.crossing = crossing;
  return tv;
}

//
// Rotation eines Vektors um phi / theta
//
// Seems to be never used: 
fourvector fourvector::newrot(double phi, double costheta) const{
  fourvector tv;
  double sintheta = sqrt(1. - pow(costheta, 2.));
  tv.t = t;
  tv.x = cos(phi) * costheta * x + sin(phi) * costheta * y - sintheta * z;
  tv.y = -sin(phi) * x + cos(phi) * y;
  tv.z = cos(phi) * sintheta * x + sin(phi) * sintheta * y + costheta * z;
  //  tv.crossing = crossing;
  return tv;
}

fourvector fourvector::rot(double phi, double costheta) const{
  fourvector tv;
  double sintheta = sqrt(1. - pow(costheta, 2.));
  double sinphi = sin(phi);
  double cosphi = cos(phi);
  tv.t = t;
  tv.x = cosphi * costheta * x - sinphi * y + cosphi * sintheta * z;
  tv.y = sinphi * costheta * x + cosphi * y + sinphi * sintheta * z;
  tv.z = -sintheta * x + costheta * z;
  //  tv.crossing = crossing;
  return tv;
}

fourvector fourvector::rot_inverse(double phi, double costheta) const{
  fourvector tv;
  double sintheta = sqrt(1. - pow(costheta, 2.));
  double sinphi = sin(phi);
  double cosphi = cos(phi);
  tv.t = t;
  tv.x = cosphi * costheta * x + sinphi * costheta * y - sintheta * z;
  tv.y = -sinphi * x + cosphi * y;
  tv.z = sintheta * cosphi * x + sinphi * sintheta * y + costheta * z;
  //  tv.crossing = crossing;
  return tv;
}


//
//  method to rotate fourvector "frot" into fourvector in x3-direction
//
fourvector fourvector::rotate(fourvector frot) const{
  double phi;
  double costheta;
  if (frot.r() > 0.){
    costheta = frot.z / frot.r();
    if (frot.x > 0.){phi = atan(frot.y / frot.x);}
    else if (frot.x < 0.){phi = atan(frot.y / frot.x) + pi;}
    else {
      if (frot.y > 0.){phi = pi / 2.;}
      else if (frot.y < 0.){phi = 3. * pi / 2.;}
      else {phi = 0.;}
    }
  }
  else {costheta = 0; phi = 0;}
  double sintheta = sqrt(1. - costheta * costheta);
  double sinphi = sin(phi);
  double cosphi = cos(phi);
  fourvector tv;
  tv.t = t;
  tv.x = cosphi * costheta * x + sinphi * costheta * y - sintheta * z;
  tv.y = -sinphi * x + cosphi * y;
  tv.z = sintheta * cosphi * x + sintheta * sinphi * y + costheta * z;
  //  tv.crossing = crossing;
  return tv;
}

//
//  method to rotate fourvector in x3-direction to fourvector "frot"
//
fourvector fourvector::rotateback(fourvector frot) const{
  double phi;
  double costheta;
  if (frot.r() > 0.){
    costheta = frot.z / frot.r();
    if (frot.x > 0.){phi = atan(frot.y / frot.x);}
    else if (frot.x < 0.){phi = atan(frot.y / frot.x) + pi;}
    else {
      if (frot.y > 0.){phi = pi / 2.;}
      else if (frot.y < 0.){phi = 3. * pi / 2.;}
      else {phi = 0.;}
    }
  }
  else{costheta = 0; phi = 0;}
  double sintheta = sqrt(1. - costheta * costheta);
  double sinphi = sin(phi);
  double cosphi = cos(phi);
  fourvector tv;
  tv.t = t;
  tv.x = cosphi * costheta * x - sinphi * y + cosphi * sintheta * z;
  tv.y = sinphi * costheta * x + cosphi * y + sinphi * sintheta * z;
  tv.z = -sintheta * x + costheta * z;
  //  tv.crossing = crossing;
  return tv;
}

// TODO: check all the *rot* functions for redundancy
// fourvector::rotateback_inverse  should be identical to  fourvector::rotate !!!
fourvector fourvector::rotateback_inverse(fourvector frot) const{
  double phi;
  double costheta;
  if (frot.r() > 0.){
    costheta = frot.z / frot.r();
    if (frot.x > 0.){phi = atan(frot.y / frot.x);}
    else if (frot.x < 0.){phi = atan(frot.y / frot.x) + pi;}
    else {
      if (frot.y > 0.){phi = pi / 2.;}
      else if (frot.y < 0.){phi = 3. * pi / 2.;}
      else {phi = 0.;}
    }
  }
  else {costheta = 0; phi = 0;}
  double sintheta = sqrt(1. - costheta * costheta);
  double sinphi = sin(phi);
  double cosphi = cos(phi);
  fourvector tv;
  tv.t = t;
  tv.x = cosphi * costheta * x + sinphi * costheta * y - sintheta * z;
  tv.y = -sinphi * x + cosphi * y;
  tv.z = sintheta * cosphi * x + sinphi * sintheta * y + costheta * z;
  //  tv.crossing = crossing;
  return tv;
}

//
// Boost eines Vierervektors in das Ruhesystem des Vektors "boost"
//
fourvector fourvector::boost(fourvector boost) const{
  static Logger logger("fourvector::boost");
  
  double m2 = boost.m2();
  double m = sqrt(m2);
  //  double gamma = boost.t / m;
  double sp_iii = sp3(boost);
  //  double spg = sp_iii / (m2 * (gamma + 1.));
  //  double spg = sp_iii / (m * (m + boost.t));
  fourvector tv;
  /*
  tv.t = gamma * t - sp_iii / m;
  tv.x = x + boost.x * (spg - (t / m));
  tv.y = y + boost.y * (spg - (t / m));
  tv.z = z + boost.z * (spg - (t / m));
  */
  //  return fourvector(gamma * t - sp_iii / m, x + boost.x * (spg - (t / m)), y + boost.y * (spg - (t / m)), z + boost.z * (spg - (t / m)));
  //  tv.crossing = crossing;
  
  tv.t = (t * boost.t - sp_iii) / m;
  tv.x = x + boost.x / m * (sp_iii / (m + boost.t) - t);
  tv.y = y + boost.y / m * (sp_iii / (m + boost.t) - t);
  tv.z = z + boost.z / m * (sp_iii / (m + boost.t) - t);
  //  return fourvector((t * boost.t - sp_iii) / m, x + boost.x / m * (sp_iii / (m + boost.t) - t), y + boost.y / m * (sp_iii / (m + boost.t) - t), z + boost.z / m * (sp_iii / (m + boost.t) - t));
  
  /*
  //  return tv;
  
  double t_tmp = (t*boost.t-sp_iii)/m;
  if (fabs(tv.t/t_tmp-1)>1e-3) {
    logger << LOG_DEBUG << "boost(p): t = " << setprecision(15) << setw(23) << tv.t << ", " << setprecision(15) << setw(23) << t_tmp << "   ratio = " << setprecision(15) << setw(23) << fabs(tv.t/t_tmp) << ", m2 = " << setprecision(15) << setw(23) << m*m << endl;
  }

  double x_tmp = x + boost.x/m*(sp_iii/(m+boost.t)-t);
  if (fabs(tv.x/x_tmp-1)>1e-3) {
    //    cout << "x=" << tv.x<< ", " << x_tmp << endl;
    //    cout << fabs(tv.x/x_tmp) << ", " << m*m << endl;
    logger << LOG_DEBUG << "boost(p): x = " << setprecision(15) << setw(23) << tv.x << ", " << setprecision(15) << setw(23) << x_tmp << "   ratio = " << setprecision(15) << setw(23) << fabs(tv.x/x_tmp) << ", m2 = " << setprecision(15) << setw(23) << m*m << endl;
  }

  double y_tmp = y + boost.y/m*(sp_iii/(m+boost.t)-t);
  if (fabs(tv.y/y_tmp-1)>1e-3) {
    //    cout << "y=" << tv.y<< ", " << y_tmp << endl;
    //    cout << fabs(tv.y/y_tmp) << ", " << m*m << endl;
    logger << LOG_DEBUG << "boost(p): y = " << setprecision(15) << setw(23) << tv.y << ", " << setprecision(15) << setw(23) << y_tmp << "   ratio = " << setprecision(15) << setw(23) << fabs(tv.y/y_tmp) << ", m2 = " << setprecision(15) << setw(23) << m*m << endl;
  }

  double z_tmp = z + boost.z/m*(sp_iii/(m+boost.t)-t);
  if (fabs(tv.z/z_tmp-1)>1e-3) {
    //    cout << "z=" << tv.z<< ", " << z_tmp << endl;
    //    cout << fabs(tv.z/z_tmp) << ", " << m*m << endl;
    logger << LOG_DEBUG << "boost(p): z = " << setprecision(15) << setw(23) << tv.z << ", " << setprecision(15) << setw(23) << z_tmp << "   ratio = " << setprecision(15) << setw(23) << fabs(tv.z/z_tmp) << ", m2 = " << setprecision(15) << setw(23) << m*m << endl;
  }

  tv.t=t_tmp;
  tv.x=x_tmp;
  tv.y=y_tmp;
  tv.z=z_tmp;
  */
  return tv;
}
fourvector fourvector::boost(fourvector boost, double m2) const{
  static Logger logger("fourvector::boost");
  
  double m = sqrt(m2);
  //  double gamma = boost.t / m;
  double sp_iii = sp3(boost);
  //  double spg = sp_iii / (m2 * (gamma + 1.));
  //  double spg = sp_iii / (m * (m + boost.t));
  fourvector tv;
  /*
  tv.t = gamma * t - sp_iii / m;
  tv.x = x + boost.x * (spg - (t / m));
  tv.y = y + boost.y * (spg - (t / m));
  tv.z = z + boost.z * (spg - (t / m));
  */
  //  tv.crossing = crossing;
  
  tv.t = (t * boost.t - sp_iii) / m;
  tv.x = x + boost.x / m * (sp_iii / (m + boost.t) - t);
  tv.y = y + boost.y / m * (sp_iii / (m + boost.t) - t);
  tv.z = z + boost.z / m * (sp_iii / (m + boost.t) - t);
  
  
  //  return tv;
  /*
  // later version
  double t_tmp = (t*boost.t-sp_iii)/m;
  if (fabs(tv.t/t_tmp-1)>1e-3) {
    logger << LOG_DEBUG << "boost(p, s): t = " << setprecision(15) << setw(23) << tv.t << ", " << setprecision(15) << setw(23) << t_tmp << "   ratio = " << setprecision(15) << setw(23) << fabs(tv.t/t_tmp) << ", m2 = " << setprecision(15) << setw(23) << m*m << endl;
  }

  double x_tmp = x + boost.x/m*(sp_iii/(m+boost.t)-t);
  if (fabs(tv.x/x_tmp-1)>1e-3) {
    //    cout << "x=" << tv.x<< ", " << x_tmp << endl;
    //    cout << fabs(tv.x/x_tmp) << ", " << m*m << endl;
    logger << LOG_DEBUG << "boost(p, s): x = " << setprecision(15) << setw(23) << tv.x << ", " << setprecision(15) << setw(23) << x_tmp << "   ratio = " << setprecision(15) << setw(23) << fabs(tv.x/x_tmp) << ", m2 = " << setprecision(15) << setw(23) << m*m << endl;
  }

  double y_tmp = y + boost.y/m*(sp_iii/(m+boost.t)-t);
  if (fabs(tv.y/y_tmp-1)>1e-3) {
    //    cout << "y=" << tv.y<< ", " << y_tmp << endl;
    //    cout << fabs(tv.y/y_tmp) << ", " << m*m << endl;
    logger << LOG_DEBUG << "boost(p, s): y = " << setprecision(15) << setw(23) << tv.y << ", " << setprecision(15) << setw(23) << y_tmp << "   ratio = " << setprecision(15) << setw(23) << fabs(tv.y/y_tmp) << ", m2 = " << setprecision(15) << setw(23) << m*m << endl;
  }

  double z_tmp = z + boost.z/m*(sp_iii/(m+boost.t)-t);
  if (fabs(tv.z/z_tmp-1)>1e-3) {
    //    cout << "z=" << tv.z<< ", " << z_tmp << endl;
    //    cout << fabs(tv.z/z_tmp) << ", " << m*m << endl;
    logger << LOG_DEBUG << "boost(p, s): z = " << setprecision(15) << setw(23) << tv.z << ", " << setprecision(15) << setw(23) << z_tmp << "   ratio = " << setprecision(15) << setw(23) << fabs(tv.z/z_tmp) << ", m2 = " << setprecision(15) << setw(23) << m*m << endl;
  }
  */
  /*
  // older version
  double t_tmp = (t*boost.t-sp_iii)/m;
  if (fabs(tv.t/t_tmp-1)>1e-3) {
    cout << "t=" << tv.t<< ", " << t_tmp << endl;
    cout << fabs(tv.t/t_tmp) << ", " << m*m << endl;
  }

  double x_tmp = x + boost.x/m*(sp_iii/(m+boost.t)-t);
  if (fabs(tv.x/x_tmp-1)>1e-3) {
    cout << "x=" << tv.x<< ", " << x_tmp << endl;
    cout << fabs(tv.x/x_tmp) << ", " << m*m << endl;
  }

  double y_tmp = y + boost.y/m*(sp_iii/(m+boost.t)-t);
  if (fabs(tv.y/y_tmp-1)>1e-3) {
    cout << "y=" << tv.y<< ", " << y_tmp << endl;
    cout << fabs(tv.y/y_tmp) << ", " << m*m << endl;
  }

  double z_tmp = z + boost.z/m*(sp_iii/(m+boost.t)-t);
  if (fabs(tv.z/z_tmp-1)>1e-3) {
    cout << "z=" << tv.z<< ", " << z_tmp << endl;
    cout << fabs(tv.z/z_tmp) << ", " << m*m << endl;
  }
  */

  /*
  tv.t=t_tmp;
  tv.x=x_tmp;
  tv.y=y_tmp;
  tv.z=z_tmp;
  */
  
  return tv;
}

// Should be removed !!!
fourvector fourvector::spacelikeboost(fourvector boost) const{
  double m = sqrt(-boost.m2());
  double gamma = boost.t / m;
  double sp_iii = sp3(boost);
  double spg = sp_iii / (boost.m2() * (gamma + 1));
  fourvector tv;
  tv.t = gamma * t - sp_iii / m;
  tv.x = x + boost.x * (spg - (t / m));
  tv.y = y + boost.y * (spg - (t / m));
  tv.z = z + boost.z * (spg - (t / m));
  //  tv.crossing = crossing;
  return tv;
}
//
// z-Boost eines Vierervektors mit Parameter beta
//
fourvector fourvector::zboost(double beta) const{
  double gamma = 1. / sqrt(1. - pow(beta, 2.));
  fourvector tv;
  tv.t = gamma * (t - beta * z);
  tv.x = x;
  tv.y = y;
  tv.z = gamma * (-beta * t + z);
  //  tv.crossing = crossing;
  return tv;
  //  return fourvector(gamma * (t - beta * z), x, y, gamma * (-beta * t + z));
}

//
// x-Boost eines Vierervektors mit Parameter beta
//
fourvector fourvector::xboost(double beta) const{
  double gamma = 1. / sqrt(1. - pow(beta, 2.));
  fourvector tv;
  tv.t = gamma * (t - beta * x);
  tv.x = gamma * (-beta * t + x);
  tv.y = y;
  tv.z = z;
  //  tv.crossing = crossing;
  return tv;
  //  return fourvector(, , , );
}

//
// Rotation eines Vektors (x3-orientierten Vektor zu "frot")
//
fourvector fourvector::rot(fourvector frot) const{
  double phi, costheta;
  if (frot.r() > 0.){
    costheta = frot.z / frot.r();
    if (frot.x > 0.){phi = atan(frot.y / frot.x);}
    else if (frot.x < 0.){phi = atan(frot.y / frot.x) + pi;}
    else{
      if (frot.y > 0.){phi = pi / 2.;}
      else if(frot.y < 0.){phi = 3. * pi / 2.;}
      else{phi = 0.;}
    }
  }
  else{costheta = 0; phi = 0;}

  fourvector t = rot(phi, costheta);

  return t;
  //  return rot(phi, costheta);
}

// direct calculation of sinphi/cosphi from tanphi ???

fourvector fourvector::rot_inverse(fourvector frot) const{
  double phi, costheta;
  if (frot.r() > 0.){
    costheta = frot.z / frot.r();
    if (frot.x > 0.){phi = atan(frot.y / frot.x);}
    else if (frot.x < 0.){phi = atan(frot.y / frot.x) + pi;}
    else {
      if (frot.y > 0.){phi = pi / 2.;}
      else if(frot.y < 0.){phi = 3. * pi / 2.;}
      else{phi = 0.;}
    }
  }
  else {costheta = 0; phi = 0;}

  fourvector t = rot_inverse(phi, costheta);

  return t;
}

// is in principle working and gives a minor speed up
//fourvector fourvector::rot_opt(fourvector frot) const{
//  double costheta;
//  double sinphi,phi;
//
//  if (frot.r() > 0.){
//    costheta = frot.z / frot.r();
//    if (frot.x!=0) {
//      sinphi=frot.y/frot.x/sqrt(1+pow(frot.y/frot.x,2));
////      cout << sinphi << endl;
//      if (frot.x<0) sinphi*=-1;
////      cout << sinphi << ", " << sin(phi) << endl;
//    }
//     else {
//      if (frot.y > 0.){sinphi = 1;}
//      else if(frot.y < 0.){sinphi = -1;}
//      else{sinphi = 0.;}
//    }
//  }
//  else{costheta = 0; sinphi = 0;}
//
//  if (frot.r() > 0.){
//    //costheta = frot.z / frot.r();
//    if (frot.x > 0.){phi = atan(frot.y / frot.x);}
//    else if (frot.x < 0.){phi = atan(frot.y / frot.x) + pi;}
//    else{
//      if (frot.y > 0.){phi = pi / 2.;}
//      else if(frot.y < 0.){phi = 3. * pi / 2.;}
//      else{phi = 0.;}
//    }
//  }
//  else{/*costheta = 0;*/ phi = 0;}
//
////  if (abs(sin(phi)-sinphi)>1e-5) {
////    cout << "trouble: sinus=" <<sinphi << endl;
////  }
//
////  cout << sinphi << ", " << costheta << endl;
//
//  fourvector tv;
//  double sintheta = sqrt(1. - pow(costheta, 2.));
//  double cosphi=sqrt(1. - pow(sinphi, 2.));
//  if (frot.x<0.0)
//    cosphi*=-1;
//
////  if (abs(cos(phi)-cosphi)>1e-5) {
////    cout << "trouble" << endl;
////    cout << cos(phi) << ", " << cosphi << ", " << frot.x << endl;
////  }
//
//  tv.t = t;
//  tv.x = cosphi * costheta * x - sinphi * y + cosphi * sintheta * z;
//  tv.y = sinphi * costheta * x + cosphi * y + sinphi * sintheta * z;
//  tv.z = -sintheta * x + costheta * z;
//  tv.crossing = crossing;
//  return tv;
//}

//
// Skalarprodukt der Dreier-Vektor-Komponenten zweier Vierervektoren
//
double fourvector::sp3(fourvector fv) const{
  return (x * fv.x + y * fv.y + z * fv.z);
  /*
  double t;
  t = x * fv.x + y * fv.y + z * fv.z;
  return t;
  */
  //  return x * fv.x + y * fv.y + z * fv.z;
  ///  if (crossing != 1 || fv.crossing != 1){cout << "crossing != 1 || fv.crossing != 1" << endl;}

  /*
  cout << "sp3 x*x' = " << x * fv.x << endl;
  cout << "sp3 y*y' = " << y * fv.y << endl;
  cout << "sp3 z*z' = " << z * fv.z << endl;
  */
  ///  return crossing * fv.crossing * (x * fv.x + y * fv.y + z * fv.z);
}

//
// Summe zweier Vierervektoren
//
fourvector fourvector::operator + (const fourvector  & fv) const{
  fourvector tv;
  tv.t = t + fv.t;
  tv.x = x + fv.x;
  tv.y = y + fv.y;
  tv.z = z + fv.z;
  /*
  tv.t = crossing * t + fv.crossing * fv.t;
  tv.x = crossing * x + fv.crossing * fv.x;
  tv.y = crossing * y + fv.crossing * fv.y;
  tv.z = crossing * z + fv.crossing * fv.z;
  tv.crossing = 1;
  */
  return tv;
  //  return fourvector(, , , );
}
//
// Differenz zweier Vierervektoren
//
fourvector fourvector::operator - (const fourvector  & fv) const{
  fourvector tv;
  tv.t = t - fv.t;
  tv.x = x - fv.x;
  tv.y = y - fv.y;
  tv.z = z - fv.z;
  /*
  tv.t = crossing * t - fv.crossing * fv.t;
  tv.x = crossing * x - fv.crossing * fv.x;
  tv.y = crossing * y - fv.crossing * fv.y;
  tv.z = crossing * z - fv.crossing * fv.z;
  tv.crossing = 1;
  */
  return tv;
  //  return fourvector(, , , );
}
//
//Skalarprodukt zweier Vierervektoren
//
double fourvector::operator * (const fourvector  & fv) const{
  //  double sp;
  //  sp = t * fv.t - (x * fv.x + y * fv.y + z * fv.z);
  //  sp = crossing * t * fv.crossing * fv.t - (crossing * x * fv.crossing * fv.x + crossing * y * fv.crossing * fv.y + crossing * z * fv.crossing * fv.z);
  return t * fv.t - (x * fv.x + y * fv.y + z * fv.z);
  //  return sp;
  /*
  double sp;
  vector<double> dv;
  dv.push_back(crossing * t * fv.crossing * fv.t);
  dv.push_back(-crossing * x * fv.crossing * fv.x);
  dv.push_back(-crossing * y * fv.crossing * fv.y);
  dv.push_back(-crossing * z * fv.crossing * fv.z);
  sp = addreal(dv);
  return sp;
  */
}
//
// Gleichheit zweier Vierervektoren
//
int fourvector::operator == (const fourvector & fv ) const{
  if ((t == fv.t) && (x == fv.x) && (y == fv.y) && (z == fv.z)){return 1;}
  //  if ((t == fv.t) && (x == fv.x) && (y == fv.y) && (z == fv.z) && (crossing == fv.crossing)){return 1;}
  else {return 0;}
}
//
// Ungleichheit zweier Vierervektoren
//
int fourvector::operator != (const fourvector & fv ) const{
  if ((t != fv.t) || (x != fv.x) || (y != fv.y) || (z != fv.z)){return 1;}
  //  if ((t != fv.t) || (x != fv.x) || (y != fv.y) || (z != fv.z) || (crossing != fv.crossing)){return 1;}
  else {return 0;}
}
//
// Produkt von reeller Zahl und Vierervektor
//
fourvector operator * (const double & rn, const fourvector & fv ){
fourvector tv;
  tv.t = fv.t * rn;
  tv.x = fv.x * rn;
  tv.y = fv.y * rn;
  tv.z = fv.z * rn;
  //  tv.crossing = fv.crossing;
  return tv;
  //  return fourvector(, , , );
}
//
// Produkt von Vierervektor und reeller Zahl
//
fourvector operator * (const fourvector & fv, const double & rn ){
  fourvector tv;
  tv.t = fv.t * rn;
  tv.x = fv.x * rn;
  tv.y = fv.y * rn;
  tv.z = fv.z * rn;
  //  tv.crossing = fv.crossing;
  return tv;
  //  return fourvector(, , , );
}
//
// Quotient von Vierervektor und reeller Zahl
//
fourvector operator / (const fourvector & fv, const double & rn ){
  fourvector tv;
  tv.t = fv.t / rn;
  tv.x = fv.x / rn;
  tv.y = fv.y / rn;
  tv.z = fv.z / rn;
  //  tv.crossing = fv.crossing;
  return tv;
  //  return fourvector(, , , );
}

double angle(const fourvector & fv1, const fourvector & fv2){
  double angle, cosangle;
  //  cosangle = (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  //  if (fv1.crossing != 1 || fv2.crossing != 1){cout << "crossing != 1 || fv.crossing != 1" << endl;}
  cosangle = (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  //  cosangle = fv1.crossing * fv2.crossing * (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  angle = acos(cosangle);
  return angle;
  //  return acos((fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r()));
}

double R(const fourvector & fv1, const fourvector & fv2){
  double delta_phi, R;
  //  cosangle = (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  delta_phi = abs(fv1.phi() - fv2.phi());
  if (delta_phi > pi){delta_phi = 2. * pi - delta_phi;}
  R = sqrt(pow(fv1.eta() - fv2.eta(), 2.) + pow(delta_phi, 2.));
  return R;
  //  return sqrt(pow(fv1.eta() - fv2.eta(), 2.) + pow(delta_phi, 2.));
}

double R2_eta(const fourvector & fv1, const fourvector & fv2){
  double R2_eta;
  double delta_phi;
  //  cosangle = (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  delta_phi = abs(fv1.phi() - fv2.phi());
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  R2_eta = pow(fv1.eta() - fv2.eta(), 2) + pow(delta_phi, 2);
  // TODO: could be optimized further
  return R2_eta;
  //  return pow(fv1.eta() - fv2.eta(), 2) + pow(delta_phi, 2);
}

double R2_rapidity(const fourvector & fv1, const fourvector & fv2){
  double R2_rapidity;
  double delta_phi;
  //  cosangle = (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  delta_phi = abs(fv1.phi() - fv2.phi());
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  R2_rapidity = pow(fv1.rapidity() - fv2.rapidity(), 2) + pow(delta_phi, 2);
  return R2_rapidity;
  //  return pow(fv1.rapidity() - fv2.rapidity(), 2) + pow(delta_phi, 2);
}

double R2_cosheta_cosphi(const fourvector & fv1, const fourvector & fv2){
  double R2_cosheta_cosphi;
  double delta_phi;
  //  cosangle = (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  delta_phi = abs(fv1.phi() - fv2.phi());
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  //  R2_cosheta_cosphi = cosh(fv1.eta() - fv2.eta()) - cos(delta_phi);
  R2_cosheta_cosphi = 2 * (cosh(fv1.eta() - fv2.eta()) - cos(delta_phi));
  return R2_cosheta_cosphi;
  //  return 2 * (cosh(fv1.eta() - fv2.eta()) - cos(delta_phi));
}

double R2_coshrapidity_cosphi(const fourvector & fv1, const fourvector & fv2){
  double R2_coshrapidity_cosphi;
  double delta_phi;
  //  cosangle = (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  delta_phi = abs(fv1.phi() - fv2.phi());
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  //  R2_coshrapidity_cosphi = cosh(fv1.rapidity() - fv2.rapidity()) - cos(delta_phi);
  R2_coshrapidity_cosphi = 2 * (cosh(fv1.rapidity() - fv2.rapidity()) - cos(delta_phi));
  return R2_coshrapidity_cosphi;
  //  return 2 * (cosh(fv1.rapidity() - fv2.rapidity()) - cos(delta_phi));
}

double cosangle(const fourvector & fv1, const fourvector & fv2){
  double cosangle;
  //  cosangle = (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  //  if (fv1.crossing != 1 || fv2.crossing != 1){cout << "crossing != 1 || fv.crossing != 1" << endl;}
  cosangle = (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  //  cosangle = fv1.crossing * fv2.crossing * (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  return cosangle;
  //  return (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
}

double sinangle(const fourvector & fv1, const fourvector & fv2){
  double sinangle;
  //  cosangle = (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r());
  //  if (fv1.crossing != 1 || fv2.crossing != 1){cout << "crossing != 1 || fv.crossing != 1" << endl;}
  sinangle = sqrt(1. - pow((fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r()), 2));
  //  sinangle = sqrt(1. - pow(fv1.crossing * fv2.crossing * (fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r()), 2));
  return sinangle;
  //  return sqrt(1. - pow((fv1.x * fv2.x + fv1.y * fv2.y + fv1.z * fv2.z) / (fv1.r() * fv2.r()), 2));
}


double fourvector::angle(const fourvector & fv) const{
  //double angle(const fourvector & fv) const{
  /*
  double angle, cosangle;
  cosangle = (x * fv.x + y * fv.y + z * fv.z) / (r() * fv.r());
  angle = acos(cosangle);
  return angle;
*/
  //  if (crossing != 1 || fv.crossing != 1){cout << "crossing != 1 || fv.crossing != 1" << endl;}
  return acos((x * fv.x + y * fv.y + z * fv.z) / (r() * fv.r()));
  //  return acos(crossing * fv.crossing * (x * fv.x + y * fv.y + z * fv.z) / (r() * fv.r()));
}

//
// Standardausgabe eines Vierervektors
//
//std::ostream & operator << (std::ostream &s, const fourvector & fv){
ostream & operator << (ostream &s, const fourvector & fv){
  //  if (fv.crossing == -1){s << "crossed ";}
  s << "( " << setw(21) << setprecision(15) << fv.t << " ; " << setw(21) << setprecision(15) << fv.x << " , " << setw(21) << setprecision(15) << fv.y << " , " << setw(21) << setprecision(15) << fv.z <<" )";
  //  s << "( " << setw(14) << setprecision(8) << fv.t << " ; " << setw(14) << setprecision(8) << fv.x << " , " << setw(14) << setprecision(8) << fv.y << " , " << setw(14) << setprecision(8) << fv.z <<" )";
  return s;
}
