// varia/readin.help.cpp
void system_execute(Logger & logger, string xorder);
void system_execute(Logger & logger, string xorder, int isystem);
void fill_charge_particle(map<int, double> & charge_particle);
void fill_code_particle(map<string, int> & code_particle);
void fill_name_particle(map<int, string> & name_particle);
void fill_latexname(map<int, string> & latexname);
void fill_gnuplotname(map<int, string> & gnuplotname);
void fill_pname(map<int, string> & pname);
void fill_datname(map<int, string> & datname);
void fill_pmadgraph(map<int, string> & pmadgraph);
void determine_process(string & subprocess, vector<vector<string> > & particles_name, int & process_type, int & n_particles, vector<int> & pa, map<string, int> & code_particle);
void get_simple_userinput_from_readin(vector<string> & user_variable, vector<string> & user_value, vector<string> & readin);
void get_userinput_from_readin(vector<string> & user_variable, vector<string> & user_variable_additional, vector<string> & user_value, vector<string> & readin);
void get_userinput_extravalue_from_readin(vector<string> & user_variable, vector<string> & user_variable_additional, vector<string> & user_value, vector<string> & readin);
void temp_define_object_list(map<string, int> & observed_object, vector<string> & object_list, vector<int> & object_category);



// varia/readin.model.cpp
//parameter_model readin_model(int N_f, int N_f_active, int perturbative_QCD_order, int max_perturbative_QCD_order, vector<string> & contribution_LHAPDFname, vector<int> & contribution_LHAPDFsubset, int switch_alpha_CMS, int switch_costheta_real);



// varia/readin.process.cpp
void parameter_readin_file(string filename, vector<string> & readin, bool essential, int retry=1);
void parameter_readin(int iswitch, string & subprocess, vector<string> & readin);
void parameter_readin(vector<string> & readin);
void parameter_readin(vector<string> & readin, string & contribution_directory);



// varia/readin.subprocess.cpp
void subprocess_readin(string & process_class, string & subprocess, vector<string> & decay, int & process_type, int & n_particle, vector<vector<int> > & type_parton);



// varia/readin.define.object.cpp
// void determine_object_definition(vector<int> & n_observed_min, vector<int> & n_observed_max, vector<double> & define_pT, vector<double> & define_ET, vector<double> & define_eta, vector<double> & define_y, map<string, int> & observed_object, vector<string> & object_list);
// void determine_object_partonlevel(vector<int> & n_observed_min, vector<int> & n_observed_max, vector<int> & n_partonlevel, int & n_parton_nu, vector<vector<int> > & type_parton, map<string, int> & observed_object, vector<string> & object_list);
// void determine_equivalent_object(vector<int> & n_partonlevel, map<string, int> & observed_object, vector<string> & object_list, map<string, string> & equivalent_object);
// void determine_relevant_object(vector<int> & n_observed_min, vector<int> & n_observed_max, vector<double> & define_pT, vector<double> & define_ET, vector<double> & define_eta, vector<double> & define_y, vector<string> & object_list, vector<int> & object_category, vector<int> & relevant_n_observed_min, vector<int> & relevant_n_observed_max, vector<double> & relevant_define_pT, vector<double> & relevant_define_ET, vector<double> & relevant_define_eta, vector<double> & relevant_define_y, vector<string> & relevant_object_list, vector<int> & relevant_object_category, map<string, string> & equivalent_object, map<int, int> & equivalent_no_object, map<string, int> & no_relevant_object);
// void determine_runtime_object(vector<int> & n_observed_min, vector<int> & n_observed_max, vector<double> & define_pT, vector<double> & define_ET, vector<double> & define_eta, vector<double> & define_y, vector<int> & n_partonlevel, vector<vector<int> > & type_parton, map<string, int> & observed_object, vector<string> & object_list, vector<int> & object_category, vector<int> & relevant_n_observed_min, vector<int> & relevant_n_observed_max, vector<double> & relevant_define_pT, vector<double> & relevant_define_ET, vector<double> & relevant_define_eta, vector<double> & relevant_define_y, vector<int> & relevant_n_partonlevel, vector<string> & relevant_object_list, vector<int> & relevant_object_category, int & n_parton_nu, int & runtime_jet_algorithm, vector<int> & runtime_jet_recombination, vector<int> & runtime_photon_recombination, vector<int> & runtime_photon_isolation, vector<int> & runtime_original, vector<int> & runtime_missing, vector<int> & runtime_object, vector<vector<int> > & runtime_order, vector<vector<int> > & runtime_order_inverse, map<string, string> & equivalent_object, map<int, int> & equivalent_no_object, map<string, int> & no_relevant_object, int frixione_isolation, int jet_algorithm, int photon_recombination);



// varia/photon.algorithm.cpp
void photon_recombination(vector<int> & no_unrecombined_photon, int i_a, observable_set & oset);
double frixione_discr(double delta, observable_set & oset);
void frixione_isolation(int & number_photon, vector<particle> & isolated_photon, particle & photon, vector<particle> & protojet, observable_set & oset);



// varia/jet.algorithm.cpp
void jet_algorithm_flavour(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin, observable_set & oset);
void sc_Ellis_Soper_flavour(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin, observable_set & oset);
void sc_Ellis_Soper_mod_flavour(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin, observable_set & oset);
void kTrun2_flavour(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin, observable_set & oset);
void antikT_flavour(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin, observable_set & oset);
void transfer_protojet_to_jet(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, int pj1, vector<particle> & jet, vector<vector<int> > & jet_flavour, vector<vector<int> > & jet_parton_origin);
void combine_protojet(vector<particle> & protojet, vector<vector<int> > & protojet_flavour, vector<vector<int> > & protojet_parton_origin, int pj1, int pj2);


// varia/event.selection.cpp
void determine_selection_content_particle(vector<int> & content_list, vector<string> & content_selection, vector<string> & content_disable);
void perform_event_selection(observable_set & oset, call_generic & generic);
void perform_event_selection_phasespace(int i_a, observable_set & oset);

// varia/determination.ckm.factor.cpp
double_complex determine_CKM_quarkchain(vector<int> type_parton, observable_set & oset);

// tools/extrapolation.cpp
void extrapolation_quadratic(vector<double> & xvalue, double & min_extrapolation, double & max_extrapolation, vector<double> & data_result, vector<double> & data_deviation2, vector<double> & X, vector<double> & dX);
//void extrapolation_quadratic(vector<double> & xvalue, vector<double> & data_result, vector<double> & data_deviation2, vector<double> & X, vector<double> & dX);
//void extrapolation_quadratic_TSV(vector<double> & xvalue, vector<vector<vector<vector<double> > > > & data_result_TSV, vector<vector<vector<vector<double> > > > & data_deviation_TSV, vector<vector<vector<double> > > & extrapolation_result_TSV, vector<vector<vector<double> > > & extrapolation_deviation_TSV);
