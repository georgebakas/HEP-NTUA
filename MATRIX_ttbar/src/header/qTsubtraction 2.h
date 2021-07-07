// qTsubtraction/qTsubtraction.ME2.VT.cpp
void calculate_ME2_VT_QCD(      observable_set & oset, call_generic & generic);
void calculate_ME2check_VT_QCD( observable_set & oset, call_generic & generic);
void calculate_ME2_VT2_QCD(     observable_set & oset, call_generic & generic);
void calculate_ME2check_VT2_QCD(observable_set & oset, call_generic & generic);

// qTsubtraction/hpl.cpp
void compute_HPLs_weight2(double v, double *HZ1, double *HZ2, double *HZ3);
void compute_GHPLs(int weight, double u, double v, double *GYZ1, double *GYZ2, double *GYZ3, double *GYZ4, double *HZ1, double *HZ2, double *HZ3, double *HZ4);

// qTsubtraction/specialfunctions.cpp
double zbesselk0( double z);
double zbesselk1( double z);
double zbesselk2( double z);
double zbesselk3( double z);
void Itilde(const double xmio, const int order, double &LL1, double &LL2, double &LL3, double &LL4);
double myli3(     double x);
double Itilde1(   double xmio, void *params);
double Itilde2(   double xmio, void *params);
double Itilde3(   double xmio, void *params);
double Itilde4(   double xmio, void *params);
void computeItildaIntegrals(double xmin, double xstep, int steps, vector<double> &I1_int, vector<double> &I2_int, vector<double> &I3_int, vector<double> &I4_int);

// qTsubtraction/auxiliaryfunctions.cpp
double C1qqdelta(double A_F);
double C2qqreg(  double z, double nf);
double C2qg(     double z);
double C2qqb(    double z, double nf);
double C2qqp(    double z, double nf);
double D0int(    double z);
double D1int(    double z);
double Pggreg(   double z);
double Pgq(      double z);
double Pqq(      double z);
double Pqg(      double z);
double Cqq(      double z);
double Cqg(      double z);
double Cgq(      double z);
double Pqqint(   double z);
double Pqqqq(    double z);
double Pqqqg(    double z);
double Pqggq(    double z);
double Pqggg(    double z, double beta0);
double P2qqV(    double x, int nf);
double S2(       double x);
double P2qqbV(   double x);
double P2qg(     double x);
double P2qqS(    double x);
double CqqPqq(   double z);
double CqqPqg(   double z);
double CqgPgq(   double z);
double CqgPgg(   double z, double beta0);
// new          
double Pgqqg (   double z, int nf);
double Pggggreg( double z, double beta0);
double Pgqqq(    double z);
double Pgggq(    double z, double beta0);
double CgqPqg(   double z, int nf);
double CgqPqq(   double z);
double Pgg(      double z);
double P2gg(     double z, int nf);
double P2gq(     double z, int nf );
double C1ggdelta(double A_F);
double C2ggreg(  double z, int nf );
double C2gq(     double z, int nf );
double Ggg(      double z );
double Ggq(      double z );


// qTsubtraction/spinorproducts.cpp
void calcSpinorProducts(vector<fourvector> &p, vector<vector<double_complex > > &za, vector<vector<double_complex > > &zb, vector<vector <double> > &s);

// qTsubtraction/counterterm.cpp
double sigma12_qqbar(double pdf_factor_x1x2);
double sigma11_qqbar(double pdf_factor_x1x2, double tH1F, double LQ);
double sigma24_qqbar(double pdf_factor_x1x2);
double sigma23_qqbar(double pdf_factor_x1x2, double sig11, double beta0);
void calcIntermediateTerms(double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_z1z2, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb, int nf, double A_F, double beta0, double &tgaga, double &tcga, double &tgamma2);
double sigma22_qqbar(double pdf_factor_x1x2, double beta0, double sig11, double tH1, double tH1F, double H1full, double tgaga, double LR, double LF, double LQ, int nf);
double sigma21_qqbar(double pdf_factor_x1x2, double beta0, double sig11, double tH1, double tH1F, double H1full, double tgaga, double tcga, double tgamma2, double LR, double LF, double LQ, double A_F, int nf);

// new
double sigma12_gg(double pdf_factor_x1x2, double A1g);
double sigma11_gg(double pdf_factor_x1x2, double B1g, double tH1F);
double sigma24_gg(double pdf_factor_x1x2, double A1g);
double sigma23_gg(double pdf_factor_x1x2, double A1g, double sig11, double beta0);
double sigma22_gg(double pdf_factor_x1x2, double A1g, double B1g, double A2g, double beta0, double sig11, double tH1F, double H1full, double tgaga, double LR, int nf);
double sigma21_gg(double pdf_factor_x1x2, double B1g, double A_F, double beta0, double sig11, double tH1, double tH1F, double H1full, double tgaga, double tcga, double tgamma2, double LR, double LF,  int nf);
void calcIntermediateTerms_gg( double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_qz2, double pdf_factor_z1q, double pdf_factor_qq, int nf, double A_F, double beta0, double Kappa, double &tgaga, double &tcga, double &tgamma2 );
double H1_M_gg( double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_qx2, double pdf_factor_x1q, double H1_delta );
double H1_gg( double pdf_factor_x1x2, double tH1, double tH1F, double LR, double LF, double beta0 );
double H2_M_gg( double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_qz2, double pdf_factor_z1q,  double pdf_factor_qq, int nf, double H1_delta, double H2_delta );
double H2_gg(   double z1, double z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_qz2, double pdf_factor_z1q, double pdf_factor_qq, int nf, double H1_delta, double H2_delta, double LR, double LF, double tH1, double tH1F );


double H1_M_qqbar(  double z1, double z2, double g_z1, double g_z2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double H1_delta);
double H1F(   double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g);
double H1F_gg(double z1, double z2,                            double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_qx2, double pdf_factor_x1q, double beta0);
double H2_M_qqbar(  double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb, double nf, double H2_delta);
double H1_qqbar(double pdf_factor_x1x2, double tH1, double tH1F, double LR, double LF, double LQ);
double H2_qqbar(double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb,double tH1, double tH1F, double nf, double H1_delta, double H2_delta, double LR, double LF, double LQ);





/*
double H1_M(double z1, double z2, double g_z1, double g_z2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double H1_delta);
double H2_M(double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb, double nf, double H2_delta);
double H1(double pdf_factor_x1x2, double tH1, double tH1F, double LR, double LF, double LQ);
double H2(double z1, double z2, double g_z1, double g_z2, double x1, double x2, double pdf_factor_x1x2, double pdf_factor_z1x2, double pdf_factor_x1z2, double pdf_factor_gx2, double pdf_factor_x1g, double pdf_factor_gg, double pdf_factor_qx2, double pdf_factor_x1q, double pdf_factor_z1z2, double pdf_factor_gz2, double pdf_factor_z1g, double pdf_factor_qbx2, double pdf_factor_x1qb,double tH1, double tH1F, double nf, double H1_delta, double H2_delta, double LR, double LF, double LQ);
*/
