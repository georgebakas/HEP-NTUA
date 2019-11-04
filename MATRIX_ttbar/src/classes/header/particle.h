#ifndef PARTICLE_H
#define PARTICLE_H

class particle {
private:

public:
////////////////////
//  constructors  //
////////////////////
  particle();
  ~particle();
  particle(fourvector momentum);
  particle(fourvector momentum, double M);
  particle(fourvector momentum, double M, double M2);
///////////////////////
//  access elements  //
///////////////////////
  particle operator + (const particle & p2) const;
  double cosphi() const;
  double sinphi() const;

  friend double R2_eta(const particle & fv1, const particle & fv2);
  friend double R2_rapidity(const particle & fv1, const particle & fv2);
  friend double R2_cosheta_cosphi(const particle & fv1, const particle & fv2);
  friend double R2_coshrapidity_cosphi(const particle & fv1, const particle & fv2);

  friend double R2_eta(const particle & p1, const fourvector & fv2);
  friend double R2_rapidity(const particle & p1, const fourvector & fv2);
  friend double R2_cosheta_cosphi(const particle & p1, const fourvector & fv2);
  friend double R2_coshrapidity_cosphi(const particle & p1, const fourvector & fv2);

  friend double R_eta(const particle & fv1, const particle & fv2);
  friend double R_rapidity(const particle & fv1, const particle & fv2);

  friend double R_eta(const particle & p1, const fourvector & fv2);
  friend double R_rapidity(const particle & p1, const fourvector & fv2);

  friend bool lessBypT(const particle & p1, const particle & p2);
  friend bool greaterBypT(const particle & p1, const particle & p2);

  friend std::ostream & operator << ( std::ostream &, const particle  &);

  
  fourvector momentum;
  double pT;
  double pT2;
  double ET;
  double ET2;
  double eta;
  double rapidity;
  double phi;
  double m;
  double m2;
  

};

bool lessBypT(const particle & p1, const particle & p2);
bool greaterBypT(const particle & p1, const particle &p2);

//bool lessByET(const particle & p1, const particle & p2);
//bool greaterByET(const particle & p1, const particle &p2);

struct particle_id{
  int identity;
  vector<int> content;
  particle *xparticle;
};
bool idgreaterBypT(const particle_id & p1, const particle_id & p2);
bool idgreaterByET(const particle_id & p1, const particle_id & p2);

#endif
