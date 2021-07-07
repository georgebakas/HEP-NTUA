#ifndef FOURVECTOR_H
#define FOURVECTOR_H

class fourvector {
private:
// coordinates, hidden
  double t;
  double x;
  double y;
  double z;
  //  int crossing;
  
 public:
// constructors
  fourvector();
  fourvector(double tv, double xv, double yv, double zv);
  fourvector(double tv, double xv, double yv, double zv, int cv);
  fourvector(double phi, double costheta, fourvector pz);
  // fourvector(bispinor_dud bs);
  // ~fourvector();
  
  // access elements
  double x0() const;
  double x1() const;
  double x2() const;
  double x3() const;
  //  int cr() const;
  // methods
  fourvector rot3(double phi, double costheta) const;
  fourvector newrot(double phi, double costheta) const;
  fourvector rot(double phi, double costheta) const;
  fourvector rot_inverse(double phi, double costheta) const;
  fourvector rot(fourvector frot) const;
  fourvector rot_inverse(fourvector frot) const;
  fourvector rotate(fourvector frot) const;
  fourvector rotateback(fourvector frot) const;
  fourvector rotateback_inverse(fourvector frot) const;
  fourvector zboost(double beta) const;
  fourvector xboost(double beta) const;
  fourvector boost(fourvector boost) const;
  fourvector boost(fourvector boost, double m2) const;
  fourvector spacelikeboost(fourvector boost) const;
  fourvector crossed() const;
  fourvector Pinv() const;
  fourvector Prot() const;
  double sp3(fourvector fv) const;
  double r() const;
  double r2() const;
  double eta() const;
  double rapidity() const;
  double theta() const;
  double costheta() const;
  // double sintheta() const;
  double pT() const;
  double pT2() const;
  double ET() const;
  double ET2() const;
  double phi() const;
  double cosphi() const;
  double sinphi() const;
  
  double m2() const;
  double m() const;
  fourvector operator - () const;
  
  fourvector operator + (const fourvector & fv ) const;
  fourvector operator - (const fourvector & fv ) const;
  double operator * (const fourvector & fv ) const;
  int operator == (const fourvector & fv ) const;
  int operator != (const fourvector & fv ) const;
  
  friend fourvector operator * (const double & rn, const fourvector & fv );
  friend fourvector operator * (const fourvector & fv, const double & rn );
  friend fourvector operator / (const fourvector & fv, const double & rn );
  friend double angle(const fourvector & fv1, const fourvector & fv2);
  friend double cosangle(const fourvector & fv1, const fourvector & fv2);
  friend double sinangle(const fourvector & fv1, const fourvector & fv2);
  friend double R(const fourvector & fv1, const fourvector & fv2);
  friend double R2_eta(const fourvector & fv1, const fourvector & fv2);
  friend double R2_rapidity(const fourvector & fv1, const fourvector & fv2);
  friend double R2_cosheta_cosphi(const fourvector & fv1, const fourvector & fv2);
  friend double R2_coshrapidity_cosphi(const fourvector & fv1, const fourvector & fv2);
  
  double angle(const fourvector & fv) const;
  
  friend std::ostream & operator << ( std::ostream &, const fourvector  &);
};

#endif
