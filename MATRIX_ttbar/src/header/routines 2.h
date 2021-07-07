//  routines/routines.generic.cpp
string get_path();
void system_execute(Logger & logger, string xorder);
void system_execute(Logger & logger, string xorder, int isystem);
int intpow(int basis, int exponent);
double g_from_alpha(double alpha);
double lambda(double x, double y, double z);
string double2hexastr(double d);
double hexastr2double(const string & s);
char * stch(string temp_s);
vector<int> get_vector_from_array_int(int * temp_array, int size_array);
string time_hms_from_double(double & time);

double max_subset_from_binary(vector<double> & xbsqrtsmin_opt, int b);
vector<int> vectorint_from_binary(int b);
vector<int> vectorbinary_from_binary(int b);
double ran(vector<double> & s);
void randomvector(vector<double> & s, int n, vector<double> & r);
int sign(int x);

bool munich_isnan(double & temp);
bool munich_isinf(double & temp);
void get_userinput_from_readin(vector<string> & user_variable, vector<string> & user_variable_additional, vector<string> & user_value, vector<string> & readin);


//  routines/routines.output.cpp
string output_subprocess(vector<int> & subprocess);
void output_weight_vegas(string & filename_alpha, vector<double> & alpha);
void output_weight_optimization(phasespace_set & psi, observable_set & oset, int & size_proc_generic, vector<int> & size_proceeding);
void output_momenta(ofstream & out_comparison, observable_set & oset);
void output_result_VA(ofstream & out_comparison, double ME2, double V_ME2, double X_ME2, double I_ME2);
void output_gnuplotfile(string & filename, int icount, double reldeviation, phasespace_set & psi);
void get_latex_subprocess(vector<string> & latex_subprocess, vector<string> & subprocess);
void get_latex_subgroup(vector<string> & latex_subgroup, vector<string> & subgroup);
string output_result_deviation(double & av_result, double & av_deviation, int digits_deviation);
string output_percent(double & av_result, double & result, int digits_deviation);
string output_latex_percent(double & av_result, double & result, int digits_deviation);
string output_commadigits(double & result, int digits);
void table_4_start(ofstream & out_res, int is, int scale_number, int columns, vector<string> v_dir_directory);
void table_4_end(ofstream & out_res);
void write_infile_int(ofstream & out_file, string name, int value_i);
void write_infile_long_long(ofstream & out_file, string name, long long value_i);
void write_infile_double(ofstream & out_file, string name, double value_d);
void write_infile_string(ofstream & out_file, string name, string value_s);



//  routines/routines.particle.code.cpp
void fill_code_particle(map<string, int> & code_particle);
void fill_name_particle(map<int, string> & name_particle);
//void fill_pname(map<int, string> & pname);
void fill_datname(map<int, string> & datname);
//void fill_pmadgraph(map<int, string> & pmadgraph);
void fill_charge_particle(map<int, double> & charge_particle);
void fill_latexname(map<int, string> & latexname);
void fill_gnuplotname(map<int, string> & gnuplotname);
void determine_process(string & subprocess, vector<vector<string> > & particles_name, int & process_type, int & n_particle, vector<int> & pa, map<string, int> & code_particle);
void subprocess_readin(string & process_class, string & subprocess, vector<string> & decay, int & process_type, int & n_particle, vector<vector<int> > & type_parton);



//  routines/routines.distribution.cpp
xdistribution get_fake_distribution_from_dddistribution(dddistribution & dddist);
void determine_distribution_complete(observable_set & oset);
double perform_binning(double value, xdistribution & this_distribution);
void determine_distribution_bin(vector<xdistribution> & dat, vector<vector<int> > & bin, vector<vector<int> > & bin_max, observable_set & oset);
void determine_distribution(vector<xdistribution> & dat, vector<vector<double> > & bin_value, vector<vector<double> > & bin_dev, int switch_bin_count, vector<vector<long long> > & bin_counts, vector<vector<double> > & var_x_integrand, vector<vector<int> > & bin, vector<vector<int> > & bin_max, observable_set & oset);
void determine_dddistribution(vector<dddistribution> & dddat, vector<xdistribution> & dat, vector<vector<double> > & bin_value, vector<vector<double> > & bin_dev, int switch_bin_count, vector<vector<long long> > & bin_counts, vector<vector<double> > & integrand_D, vector<vector<int> > & bin, vector<vector<int> > & bin_max, observable_set & oset);
void determine_distribution_TSV(vector<vector<int> > & bin, vector<vector<int> > & bin_max, vector<xdistribution> & dat, observable_set & oset);


//  routines/oldversion.routines.pdf.cpp
void calculate_pdf_LHAPDF_QT(vector<vector<int> > & all_pdf, phasespace_set & psi, observable_set & oset, int order);
void calculate_pdf_LHAPDF_LHC_QT(vector<vector<int> > & all_pdf, int sd, int ss, phasespace_set & psi, observable_set & oset, int order);
void calculate_pdf_LHAPDF_Tevatron_QT(vector<vector<int> > & all_pdf, int sd, int ss, phasespace_set & psi, observable_set & oset, int order);
