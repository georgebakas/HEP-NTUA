#include "../include/classes.cxx"
////////////////////
//  constructors  //
////////////////////
particle::particle(){
  momentum = fourvector(0., 0., 0., 0.);
  pT2 = 0.;
  pT = 0.;
  ET2 = 0.;
  ET = 0.;
  m2 = 0.;
  m = 0.;
  eta = 0.;
  rapidity = 0.;
  phi = 0.;
}

particle::~particle(){}

particle::particle(fourvector _momentum){
  momentum = _momentum;
  pT2 = momentum.pT2();
  pT = sqrt(pT2);
  ET2 = momentum.ET2();
  ET = sqrt(ET2);
  m2 = momentum.m2();
  m = sqrt(m2);
  if (pT == 0.){
    eta = 0.;
    rapidity = 0.;
    phi = 0.;
  }
  else {
    eta = momentum.eta();
    rapidity = momentum.rapidity();
    phi = momentum.phi();
  }
}
particle::particle(fourvector _momentum, double _M){
  momentum = _momentum;
  pT2 = momentum.pT2();
  pT = sqrt(pT2);
  ET2 = momentum.ET2();
  ET = sqrt(ET2);
  m = _M;
  m2 = pow(m, 2);
  if (pT == 0.){
    eta = 0.;
    rapidity = 0.;
    phi = 0.;
  }
  else {
    eta = momentum.eta();
    rapidity = momentum.rapidity();
    phi = momentum.phi();
  }
}
particle::particle(fourvector _momentum, double _M, double _M2){
  momentum = momentum;
  pT2 = momentum.pT2();
  pT = sqrt(pT2);
  ET2 = momentum.ET2();
  ET = sqrt(ET2);
  m = _M;
  m2 = _M2;
  if (pT == 0.){
    eta = 0.;
    rapidity = 0.;
    phi = 0.;
  }
  else {
    eta = momentum.eta();
    rapidity = momentum.rapidity();
    phi = momentum.phi();
  }
}
///////////////////////
//  access elements  //
///////////////////////

///////////////
//  methods  //
///////////////
double particle::cosphi() const {return cos(phi);}
double particle::sinphi() const {return sin(phi);}



particle particle::operator + (const particle & p2) const{
  particle psum;
  psum.momentum = momentum + p2.momentum;
  psum.pT2 = (psum.momentum).pT2();
  psum.pT = sqrt(psum.pT2);
  psum.ET2 = (psum.momentum).ET2();
  psum.ET = sqrt(psum.ET2);
  psum.m2 = (psum.momentum).m2();
  psum.m = sqrt(psum.m2);
  psum.eta = (psum.momentum).eta();
  psum.rapidity = (psum.momentum).rapidity();
  psum.phi = (psum.momentum).phi();
  return psum;
}

double R2_eta(const particle & p1, const particle & p2){
  double delta_phi = abs(p1.phi - p2.phi);
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return pow(p1.eta - p2.eta, 2) + pow(delta_phi, 2);
}
double R2_rapidity(const particle & p1, const particle & p2){
  double delta_phi = abs(p1.phi - p2.phi);
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return pow(p1.rapidity - p2.rapidity, 2) + pow(delta_phi, 2);
}
double R2_cosheta_cosphi(const particle & p1, const particle & p2){
  double delta_phi = abs(p1.phi - p2.phi);
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return 2 * (cosh(p1.eta - p2.eta) - cos(delta_phi));
}
double R2_coshrapidity_cosphi(const particle & p1, const particle & p2){
  double delta_phi = abs(p1.phi - p2.phi);
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return 2 * (cosh(p1.rapidity - p2.rapidity) - cos(delta_phi));
}


double R2_eta(const particle & p1, const fourvector & fv2){
  double delta_phi = abs(p1.phi - fv2.phi());
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return pow(p1.eta - fv2.eta(), 2) + pow(delta_phi, 2);
}
double R2_rapidity(const particle & p1, const fourvector & fv2){
  double delta_phi = abs(p1.phi - fv2.phi());
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return pow(p1.rapidity - fv2.rapidity(), 2) + pow(delta_phi, 2);
}
double R2_cosheta_cosphi(const particle & p1, const fourvector & fv2){
  double delta_phi = abs(p1.phi - fv2.phi());
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return 2 * (cosh(p1.eta - fv2.eta()) - cos(delta_phi));
}
double R2_coshrapidity_cosphi(const particle & p1, const fourvector & fv2){
  double delta_phi = abs(p1.phi - fv2.phi());
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return 2 * (cosh(p1.rapidity - fv2.rapidity()) - cos(delta_phi));
}



double R_eta(const particle & p1, const particle & p2){
  double delta_phi = abs(p1.phi - p2.phi);
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return sqrt(pow(p1.eta - p2.eta, 2) + pow(delta_phi, 2));
}
double R_rapidity(const particle & p1, const particle & p2){
  double delta_phi = abs(p1.phi - p2.phi);
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return sqrt(pow(p1.rapidity - p2.rapidity, 2) + pow(delta_phi, 2));
}


double R_eta(const particle & p1, const fourvector & fv2){
  double delta_phi = abs(p1.phi - fv2.phi());
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return sqrt(pow(p1.eta - fv2.eta(), 2) + pow(delta_phi, 2));
}
double R_rapidity(const particle & p1, const fourvector & fv2){
  double delta_phi = abs(p1.phi - fv2.phi());
  if (delta_phi > pi){delta_phi = 2 * pi - delta_phi;}
  return sqrt(pow(p1.rapidity - fv2.rapidity(), 2) + pow(delta_phi, 2));
}


ostream & operator << (ostream &s, const particle & p){
  s << p.momentum << endl << setw(10) << "" << "m = " << setw(21) << setprecision(15) << p.m << endl << "   pT = " << setw(21) << setprecision(15) << p.pT << "  ET = " << setw(21) << setprecision(15) << p.ET << "  eta = " << setw(21) << setprecision(15) << p.eta << "  y = " << setw(21) << setprecision(15) << p.rapidity << "  phi = " << setw(21) << setprecision(15) << p.phi;
  //  s << "( " << setw(14) << setprecision(8) << fv.t << " ; " << setw(14) << setprecision(8) << fv.x << " , " << setw(14) << setprecision(8) << fv.y << " , " << setw(14) << setprecision(8) << fv.z <<" )";
  return s;
}




bool lessBypT(const particle & p1, const particle & p2){return p1.pT < p2.pT;}
bool greaterBypT(const particle & p1, const particle & p2){return p1.pT > p2.pT;}
/*
bool greaterBypT(const particle & p1, const particle & p2){
  if (abs(p1.pT - p2.pT) / (p1.pT + p2.pT) < 1.e-12){
    //    cout << "Two back-to-back particles in greaterBypT:  p1.pT = " << p1.pT << "   p2.pT = " << p2.pT << endl;
    if (p1.momentum.x1() > p2.momentum.x2()){return true;}
    else {return false;}
  }
  else if (abs(p1.pT - p2.pT) / (p1.pT + p2.pT) < 1.e-09){cout << "WARNING: p1.pT and p2.pT are almost identical, but not identified !!!   " << abs(p1.pT - p2.pT) / (p1.pT + p2.pT) << endl;}
  else {return p1.pT > p2.pT;}
}
*/

bool idgreaterBypT(const particle_id & p1, const particle_id & p2){return (*p1.xparticle).pT > (*p2.xparticle).pT;}

bool idgreaterByET(const particle_id & p1, const particle_id & p2){return (*p1.xparticle).ET > (*p2.xparticle).ET;}


