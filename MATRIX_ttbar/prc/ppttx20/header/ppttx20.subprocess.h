void ppttx20_list_subprocess_born(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
void ppttx20_list_subprocess_C_QCD(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
void ppttx20_list_subprocess_V_QCD(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
void ppttx20_list_subprocess_C2_QCD(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
void ppttx20_list_subprocess_V2_QCD(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
void ppttx20_list_subprocess_R_QCD(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
void ppttx20_list_subprocess_RC_QCD(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
void ppttx20_list_subprocess_RV_QCD(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
void ppttx20_list_subprocess_RR_QCD(vector<string> & subprocess, vector<vector<int> > & subgroup_no_member, int contribution_order_alpha_s, int contribution_order_alpha_e, int contribution_order_interference);
void ppttx20_determination_no_subprocess_born(int & no_map, vector<int> & o_map, int & no_prc, vector<int> & o_prc, double & symmetry_factor, vector<int> & tp, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference);
void ppttx20_determination_subprocess_born(int i_a, phasespace_set & psi);

void ppttx20_combination_subprocess_born(int i_a, phasespace_set & psi, observable_set & oset);
void ppttx20_determination_no_subprocess_real(int & no_map, vector<int> & o_map, int & no_prc, vector<int> & o_prc, double & symmetry_factor, vector<int> & tp, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference);
void ppttx20_determination_subprocess_real(int i_a, phasespace_set & psi);

void ppttx20_combination_subprocess_real(int i_a, phasespace_set & psi, observable_set & oset);
void ppttx20_determination_no_subprocess_doublereal(int & no_map, vector<int> & o_map, int & no_prc, vector<int> & o_prc, double & symmetry_factor, vector<int> & tp, int basic_order_alpha_s, int basic_order_alpha_e, int basic_order_interference);
void ppttx20_determination_subprocess_doublereal(int i_a, phasespace_set & psi);

void ppttx20_combination_subprocess_doublereal(int i_a, phasespace_set & psi, observable_set & oset);
