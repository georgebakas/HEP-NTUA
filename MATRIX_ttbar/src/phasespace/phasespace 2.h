//double max_subset_from_binary(vector<double> & xbsqrtsmin_opt, int b);
//vector<int> vectorint_from_binary(int b);
//int calc_strange_number(int b, int &size);
//vector<int> vectorbinary_from_binary(int b);
//double calc_strange_sum(int b);
//double multdouble(vector<double> mult);
//double lambda(double x, double y, double z);
//double ran(vector<double> & s);
//void randomvector(vector<double> & s, int n, vector<double> & r);

//double h_newmap(double r1, double r2, double smin, int & channel, vector<double> & tau_alpha, vector<double> & tau_beta);
//double g_newmap(double s, double smin, int channel, vector<double> & tau_alpha);

// phasespace/phasespace.initial.cpp
void initialization_fake_dipole_mapping_RT(phasespace_set & psi, call_generic & generic);
void initialization_fake_dipole_mapping_RRA(vector<dipole_set> & dipole, phasespace_set & psi, call_generic & generic);
//void ps_initial_collinear_z1_z2(vector<double> & s, vector<double> & x_pdf, vector<double> & z_coll, vector<double> & g_z_coll, vector<vector<double> > & all_xz_coll_pdf, vector<int> & z1z2_channel, vector<vector<double> > & z1z2_alpha, vector<vector<double> > & z1z2_beta, phasespace_set & psi);


// phasespace/phasespace.generator.cpp
double h_pot_min(double r, double smin, double smax, double n);
double g_pot_min(double s, double smin, double smax, double n);
double h_pot_max(double r, double smin, double smax, double n);
double g_pot_max(double s, double smin, double smax, double n);
double h_propto_pot(double r, double smin, double n);
double g_propto_pot(double s, double smin, double n);
double h_propto_pot(double r, double smin, double smax, double nu, double cut_technical);
double g_propto_pot(double s, double smin, double smax, double nu, double cut_technical);
double inv_propto_pot(double s, double smin, double smax, double nu, double cut_technical);

double h_propto_pot_mod(double r, double smin, double smax, double nu, double cut_technical);
double g_propto_pot_mod(double s, double smin, double smax, double nu, double cut_technical);
double h_propto_pot_symm(double r, double smin, double smax, double nu, double cut_technical);
double g_propto_pot_symm(double s, double smin, double smax, double nu, double cut_technical);

/*
double c_phi(double r);
double inv_phi(double phi);

double c_costheta(double r);
double inv_costheta(double costheta);
*/

/*
void c_propagator_narrow_width();
double g_propagator_narrow_width(int m, phasespace_set & psi);
void inv_propagator_narrow_width();

void c_propagator(double & r, int m, int out, vector<double> & s, vector<double> & sqrts, double smin, double smax, phasespace_set & psi);
double g_propagator(int m, int out, vector<double> & s, double smin, double smax, phasespace_set & psi);
void inv_propagator(double & r, int m, int out, vector<double> & s, vector<double> & sqrts, double smin, double smax, phasespace_set & psi);
double c_propagator_Breit_Wigner(double & r, int m, double & smin, double & smax, phasespace_set & psi);
double g_propagator_Breit_Wigner(double & s, int m, double & smin, double & smax, phasespace_set & psi);
double inv_propagator_Breit_Wigner(double s, int m, double & smin, double & smax, phasespace_set & psi);
double c_propagator_vanishing_width(double r, double mass2, double smin, double smax, double nu);
double g_propagator_vanishing_width(double s, double mass2, double smin, double smax, double nu);
double inv_propagator_vanishing_width(double s, double mass2, double smin, double smax, double nu);
void c_timelikeinvariant(double r, int out, vector<double> & s, vector<double> & sqrts, double smin, double smax);
void inv_timelikeinvariant(double & r, int out, vector<double> & s, vector<double> & sqrts, double smin, double smax);
double g_timelikeinvariant(double smin, double smax);
double c_t_propagator(double r, int m, double smin, double smax, phasespace_set & psi);
double g_t_propagator(double s, int m, double smin, double smax, phasespace_set & psi);
double inv_t_propagator(double s, int m, double smin, double smax, phasespace_set & psi);
void c_decay(double r1, double r2, int in, int out1, int out2, vector<fourvector> & p, vector<double> & s, vector<double> & sqrts);
double g_decay(int in, int out1, int out2, vector<double> & s);
void inv_decay(double & r1, double & r2, int in, int out1, int out2, vector<fourvector> & p, vector<double> & s, vector<double> & sqrts);
void c_tchannel(double r1, double r2, int m, int in, int in1, int in2, int out1, int out2, vector<fourvector> & p, vector<double> & s, phasespace_set & psi);
double g_tchannel(int m, int in, int in1, int in2, int out1, int out2, vector<fourvector> & p, vector<double> & s, phasespace_set & psi);
void inv_tchannel(double & r1, double & r2, int m, int in, int in1, int in2, int out1, int out2, vector<fourvector> & p, vector<double> & s, phasespace_set & psi);
void c_tchannel_opt(double r1, double r2, int m, int out, int in1, int in2, int out1, int out2, vector<fourvector> & p, vector<double> & s, vector<double> & sqrts, phasespace_set & psi);
double g_tchannel_opt(int m, int out, int in1, int in2, int out1, int out2, vector<fourvector> & p, vector<double> & s, vector<double> & sqrts, phasespace_set & psi);
void c_tchannel_opt_TO(double r1, double r2, int m, int out, int in1, int in2, int out1, int out2, vector<fourvector> & p, vector<double> & s, vector<double> & sqrts, phasespace_set & psi);
double g_tchannel_opt_TO(int m, int out, int in1, int in2, int out1, int out2, vector<fourvector> & p, vector<double> & s, vector<double> & sqrts, phasespace_set & psi);
void inv_tchannel_opt(double & r1, double & r2, int m, int out, int in1, int in2, int out1, int out2, vector<fourvector> & p, vector<double> & s, vector<double> & sqrts, phasespace_set & psi);
double c_subtraction_variable(double r, double mass2, double xmin, double xmax, double nu);
double g_subtraction_variable(double x, double mass2, double xmin, double xmax, double nu);
*/
