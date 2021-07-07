#include "../include/classes.cxx"
//#include "../include/definitions.cxx"
#include "../include/definitions.observable.set.cxx"
#include "../include/definitions.phasespace.set.cxx"

//#include "move.routines.output.cpp"

string output_subprocess(vector<int> & subprocess){
  stringstream temp_ss;
  for (int i_p = 1; i_p < subprocess.size(); i_p++){
    temp_ss << setw(4) << subprocess[i_p] << "";
    if (i_p == 2){temp_ss << " -> ";}
  }
  return temp_ss.str();
}

//////////////////////////////////
//  output weight optimization  //
//////////////////////////////////

void output_weight_vegas(string & filename_alpha, vector<double> & alpha){
  Logger logger("output_weight_vegas");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  ofstream out_alpha;
  out_alpha.open(filename_alpha.c_str(), ofstream::out | ofstream::trunc);
  for (int j = 0; j < alpha.size(); j++){
    out_alpha << setprecision(15) << alpha[j] << endl;
  }
  out_alpha.close();

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void output_weight_optimization(phasespace_set & psi, observable_set & oset, int & size_proc_generic, vector<int> & size_proceeding){
  Logger logger("output_weight_optimization");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  if (psi.MC_phasespace.active_optimization == 1){
    if (psi.MC_phasespace.end_optimization < 2 && *psi.i_step_mode % psi.MC_phasespace.n_alpha_events == 0){
      //  if (psi_MC_opt_end < 2 && *psi.i_step_mode % psi_n_alpha_events == 0){
      ofstream out_integration;
      out_integration.open(osi_filename_integration.c_str(), ofstream::out | ofstream::app);
      out_integration << "************************************************* weight optimization step " << setw(2) << psi.MC_phasespace.i_alpha_it << " *************************************************" << endl;
      /*
      if (psi.MC_phasespace.counter_minimum_weight != 0){
	out_integration << "   " << psi.MC_phasespace.counter_minimum_weight << " alpha's are set to " << psi.MC_phasespace.reserved_minimum_weight / psi.MC_phasespace.counter_minimum_weight << "." << endl;
	///	out_integration << "   " << psi.MC_phasespace.counter_minimum_weight << " alpha's are set to " << psi_a_reserved_min / psi.MC_phasespace.counter_minimum_weight << "." << endl;
      }
      */
      out_integration << "   D_min[alpha(i),alpha(j)](step " << psi.MC_phasespace.i_alpha_it - 1 << ") = " << setw(20) << setprecision(15) << psi.MC_phasespace.diff_w[psi.MC_phasespace.i_alpha_it - 1] << endl;
      if (!(*psi.i_step_mode == psi.MC_phasespace.n_alpha_events  * psi.MC_phasespace.n_alpha_steps)){
	//    if (!(*psi.i_step_mode == psi_n_alpha_steps * psi_n_alpha_events)){
      if (psi.MC_phasespace.counter_minimum_weight != 0){
	out_integration << "   switch_minimum_weight:   " << psi.MC_phasespace.counter_minimum_weight << " alpha's are set to " << psi.MC_phasespace.reserved_minimum_weight << " in the next step." << endl;
	  ///	  out_integration << "   Remaining " << psi.MC_phasespace.counter_minimum_weight << " alpha's are set to " << psi.MC_phasespace.reserved_minimum_weight / psi.MC_phasespace.counter_minimum_weight << "." << endl;
	  ///	  out_integration << "   Remaining " << psi.MC_phasespace.counter_minimum_weight << " alpha's are set to " << psi_a_reserved_min / psi.MC_phasespace.counter_minimum_weight << "." << endl;
      }
      //      if (!(*psi.i_step_mode == psi.MC_phasespace.n_alpha_events  * psi.MC_phasespace.n_alpha_steps)){
	out_integration << "*******************************************************************************************************************************" << endl;
      }
      out_integration.close();
    }
  
    if (psi.MC_phasespace.end_optimization < 2 && *psi.i_step_mode == psi.MC_phasespace.n_alpha_events  * psi.MC_phasespace.n_alpha_steps){
      ofstream out_integration;
      out_integration.open(osi_filename_integration.c_str(), ofstream::out | ofstream::app);
      out_integration << "************************************************* weight optimization result **************************************************" << endl;
      for (int j = 0; j < psi.MC_phasespace.diff_w.size(); j++){
	out_integration << "diff_w[" << j << "] = " << showpoint << setprecision(15) << setw(25) << psi.MC_phasespace.diff_w[j] << "   " << psi.MC_phasespace.alpha_it[j][0] << "   " << psi.MC_phasespace.alpha_it[j][psi.MC_phasespace.alpha_it[j].size() - 1] << endl;
      }
      out_integration << "************************************************* weight optimization result **************************************************" << endl;
      out_integration << "   min(D_min[alpha(i),alpha(j)]) = " << showpoint << setprecision(15) << setw(25) << psi.MC_phasespace.diff_w[psi.MC_phasespace.x_minimum_diff_w] << "   at optimization step "<< setw(2) << psi.MC_phasespace.x_minimum_diff_w << endl;
      out_integration << "************************************************* end of weight optimization **************************************************" << endl;
      out_integration.close();
      
      size_proceeding[2] = 0;
      size_proc_generic = accumulate(size_proceeding.begin(), size_proceeding.end(), 0);
      ///    psi_MC_opt_end = 2;
    }
  }

  if (psi.IS_tau.active_optimization == 1){
    if (psi.IS_tau.end_optimization < 2 && (*psi.i_step_mode % psi_n_tau_events == 0 || *psi.i_step_mode != 0)){
      if (*psi.i_step_mode == psi.IS_tau.n_event_per_step * psi.IS_tau.n_optimization_step){
	size_proceeding[3] = 0;
	size_proc_generic = accumulate(size_proceeding.begin(), size_proceeding.end(), 0);
      }
    }
  }
  
  if (psi.IS_x1x2.active_optimization == 1){
    if (psi.IS_x1x2.end_optimization < 2 && (*psi.i_step_mode % psi_n_x1x2_events == 0 || *psi.i_step_mode != 0)){
      if (*psi.i_step_mode == psi.IS_x1x2.n_event_per_step * psi.IS_x1x2.n_optimization_step){
	size_proceeding[4] = 0;
	size_proc_generic = accumulate(size_proceeding.begin(), size_proceeding.end(), 0);
      }
    }
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



//////////////////////
//  output momenta  //
//////////////////////

void output_momenta(ofstream & out_comparison, observable_set & oset){
  Logger logger("output_momenta");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  for (int ib = 1; ib < osi_p_parton[0].size(); ib++){
    out_comparison << "   p[" << setw(3) << osi_type_parton[0][ib] << "] = (" <<
      setprecision(15) << setw(19) << osi_p_parton[0][ib].x0() << "; " <<
      setprecision(15) << setw(19) << osi_p_parton[0][ib].x1() << ", " <<
      setprecision(15) << setw(19) << osi_p_parton[0][ib].x2() << ", " <<
      setprecision(15) << setw(19) << osi_p_parton[0][ib].x3() << ")" <<
      "   " << sqrt(abs(osi_p_parton[0][ib].m2())) << endl;
  }
  out_comparison << endl;
  fourvector check;
  for (int ib = 1; ib < osi_p_parton[0].size(); ib++){
    if (ib < 3){check = check + osi_p_parton[0][ib];}
    else{check = check - osi_p_parton[0][ib];}
  }
  out_comparison << "momentum conservation check = (" <<
    setprecision(15) << setw(19) << check.x0() << "; " <<
    setprecision(15) << setw(19) << check.x1() << ", " <<
    setprecision(15) << setw(19) << check.x2() << ", " <<
    setprecision(15) << setw(19) << check.x3() << ")" <<
    "   " << sqrt(abs(check.m2())) << endl;
  out_comparison << "sqrt(abs((p(in) - p(out))^2 / p(in)^2)) = " << sqrt(abs(check.m2() / osi_p_parton[0][0].m2())) << endl;
  out_comparison << endl;

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}





/////////////////////
//  output ME2 VA  //
/////////////////////

void output_result_VA(ofstream & out_comparison, double ME2, double V_ME2, double X_ME2, double I_ME2){
  Logger logger("output_result_VA");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  out_comparison << setw(12) << " ME2_born = " << setprecision(15) << setw(23) << ME2 << endl;
  out_comparison << setw(12) << "ME2_V+X+I = " << setprecision(15) << setw(23) << (V_ME2 + X_ME2 + I_ME2) << endl;
  out_comparison << endl;
  out_comparison << setw(12) << "    ME2_V = " << setprecision(15) << setw(23) << V_ME2 << endl;
  out_comparison << setw(12) << "    ME2_X = " << setprecision(15) << setw(23) << X_ME2 << endl;
  out_comparison << setw(12) << "    ME2_I = " << setprecision(15) << setw(23) << I_ME2 << endl;
  out_comparison << endl;
}





//////////////////////
//  gnuplot output  //
//////////////////////

void output_gnuplotfile(string & filename, int icount, double reldeviation, phasespace_set & psi){
  Logger logger("output_gnuplotfile");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  char LineBuffer[256];
  //  string filename = osi_filename_makegnuplot;

  //  cout << filename << endl;
  /*
  filename.erase(filename.end() - 4, filename.end());
  cout << filename << endl;
  */
  ifstream in_gnuplot(filename.c_str());
  vector<string> readin;
  while (in_gnuplot.getline(LineBuffer, 256)){
    readin.push_back(LineBuffer);
  }
  in_gnuplot.close();
  ofstream out_makegnuplot;

  out_makegnuplot.open(filename.c_str(), ofstream::out | ofstream::trunc);

  double xvalue = double(icount) / double(psi_n_events_max);
  double yvalue = reldeviation;
  double yvalue_min = yvalue * sqrt(xvalue);
  //  cout << icount << endl;
  //  cout << "xvalue = " << xvalue << "   " << "yvalue  = " << yvalue << "   " << "yvalue_min = " << yvalue_min << endl;
  int log_min = int(log10(yvalue_min)) - 1;
  //  cout << "log_min = " << log_min << endl;
  //  cout << "estimated final relative error: " << yvalue_min << endl;
  for (int i = 0; i < readin.size(); i++){
    if (i == 12){
      out_makegnuplot  << "set yrange [" << pow(10.,log_min) << ":" << 5.e03 * pow(10.,log_min) << "]" << endl;
    }
    else{out_makegnuplot << readin[i] << endl;}
  }
  out_makegnuplot.close();
  //  string xorder = "chmod 700 " + xfilename_sgnuplot;
  //  cout << xorder << endl;
  //  system(xorder.c_str());

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void get_latex_subprocess(vector<string> & latex_subprocess, vector<string> & subprocess){
  Logger logger("get_latex_subprocess");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  int bc;
  char bs = char(92);
  string bss;
  bss.push_back(bs);
  for (int j = 0; j < latex_subprocess.size(); j++){
    bc = 0;
    //    cout << "latex_subprocess["<< j << "] = " << latex_subprocess[j] << endl;
    for (int k = 0; k < latex_subprocess[j].size(); k++){
      if      (latex_subprocess[j][k] == 'd' && latex_subprocess[j][k+1] == '~'){latex_subprocess[j].replace(k, 2, bss + "mathrm{" + bss + "bar{d}}"); k+=14;}
      else if (latex_subprocess[j][k] == 'u' && latex_subprocess[j][k+1] == '~'){latex_subprocess[j].replace(k, 2, bss + "mathrm{" + bss + "bar{u}}"); k+=14;}
      else if (latex_subprocess[j][k] == 's' && latex_subprocess[j][k+1] == '~'){latex_subprocess[j].replace(k, 2, bss + "mathrm{" + bss + "bar{s}}"); k+=14;}
      else if (latex_subprocess[j][k] == 'c' && latex_subprocess[j][k+1] == '~'){latex_subprocess[j].replace(k, 2, bss + "mathrm{" + bss + "bar{c}}"); k+=14;}
      else if (latex_subprocess[j][k] == 'b' && latex_subprocess[j][k+1] == '~'){latex_subprocess[j].replace(k, 2, bss + "mathrm{" + bss + "bar{b}}"); k+=14;}
      else if (latex_subprocess[j][k] == 't' && latex_subprocess[j][k+1] == '~'){latex_subprocess[j].replace(k, 2, bss + "mathrm{" + bss + "bar{t}}"); k+=14;}
      else if (latex_subprocess[j][k] == 'd'){latex_subprocess[j].replace(k, 1, bss + "mathrm{d}"); k+=9;}
      else if (latex_subprocess[j][k] == 'u'){latex_subprocess[j].replace(k, 1, bss + "mathrm{u}"); k+=9;}
      else if (latex_subprocess[j][k] == 's'){latex_subprocess[j].replace(k, 1, bss + "mathrm{s}"); k+=9;}
      else if (latex_subprocess[j][k] == 'c'){latex_subprocess[j].replace(k, 1, bss + "mathrm{c}"); k+=9;}
      else if (latex_subprocess[j][k] == 'b'){latex_subprocess[j].replace(k, 1, bss + "mathrm{b}"); k+=9;}
      else if (latex_subprocess[j][k] == 't'){latex_subprocess[j].replace(k, 1, bss + "mathrm{t}"); k+=9;}
      else if (latex_subprocess[j][k] == 'e' && latex_subprocess[j][k+1] == 'p'){latex_subprocess[j].replace(k, 2, bss + "mathrm{e^+}"); k+=10;}
      else if (latex_subprocess[j][k] == 'e' && latex_subprocess[j][k+1] == 'm'){latex_subprocess[j].replace(k, 2, bss + "mathrm{e^-}"); k+=10;}
      else if (latex_subprocess[j][k] == 'v' && latex_subprocess[j][k+1] == 'e' && latex_subprocess[j][k+2] == 'b'){latex_subprocess[j].replace(k, 3, bss + "bar{" + bss + "nu}_" + bss + "mathrm{e}"); k+=17;}
      else if (latex_subprocess[j][k] == 'v' && latex_subprocess[j][k+1] == 'e'){latex_subprocess[j].replace(k, 2, bss + "nu_" + bss + "mathrm{e}"); k+=12;}
      else if (latex_subprocess[j][k] == 'm' && latex_subprocess[j][k+1] == 'u' && latex_subprocess[j][k+2] == 'p'){latex_subprocess[j].replace(k, 3, bss + "mathrm{" + bss + "mu^+}"); k+=11;}
      else if (latex_subprocess[j][k] == 'm' && latex_subprocess[j][k+1] == 'u' && latex_subprocess[j][k+2] == 'm'){latex_subprocess[j].replace(k, 3, bss + "mathrm{" + bss + "mu^-}"); k+=11;}
      else if (latex_subprocess[j][k] == 'v' && latex_subprocess[j][k+1] == 'm' && latex_subprocess[j][k+2] == 'b'){latex_subprocess[j].replace(k, 3, bss + "bar{" + bss + "nu}_" + bss + "mu "); k+=13;}
      else if (latex_subprocess[j][k] == 'v' && latex_subprocess[j][k+1] == 'm'){latex_subprocess[j].replace(k, 2, bss + "nu_" + bss + "mu"); k+=5;}
      else if (latex_subprocess[j][k] == 'w' && latex_subprocess[j][k+1] == 'p'){latex_subprocess[j].replace(k, 2, bss + "mathrm{W^+}"); k+=11;}
      else if (latex_subprocess[j][k] == 'w' && latex_subprocess[j][k+1] == 'm'){latex_subprocess[j].replace(k, 2, bss + "mathrm{W^-}"); k+=11;}
      else if (latex_subprocess[j][k] == 'z'){latex_subprocess[j].replace(k, 1, bss + "mathrm{Z}"); k+=9;}
      else if (latex_subprocess[j][k] == '_' && bc == 0){latex_subprocess[j].replace(k, 1, bss + "rightarrow "); bc = 1; k+=11;}
    }
    latex_subprocess[j] = '$' + latex_subprocess[j] + '$';
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}


void get_latex_subgroup(vector<string> & latex_subgroup, vector<string> & subgroup){
  Logger logger("get_latex_subgroup");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;
  char bs = char(92);
  string bss;
  bss.push_back(bs);
  for (int j = 0; j < latex_subgroup.size(); j++){
    for (int k = 0; k < latex_subgroup[j].size(); k++){
     if (latex_subgroup[j][k] == 'q' && latex_subgroup[j][k+1] == '~'){latex_subgroup[j].replace(k, 2, bss + "bar{q}");}
     if (latex_subgroup[j][k] == 'd' && latex_subgroup[j][k+1] == '~'){latex_subgroup[j].replace(k, 2, bss + "bar{d}");}
     if (latex_subgroup[j][k] == 'u' && latex_subgroup[j][k+1] == '~'){latex_subgroup[j].replace(k, 2, bss + "bar{u}");}
     if (latex_subgroup[j][k] == 's' && latex_subgroup[j][k+1] == '~'){latex_subgroup[j].replace(k, 2, bss + "bar{s}");}
     if (latex_subgroup[j][k] == 'c' && latex_subgroup[j][k+1] == '~'){latex_subgroup[j].replace(k, 2, bss + "bar{c}");}
     if (latex_subgroup[j][k] == 'b' && latex_subgroup[j][k+1] == '~'){latex_subgroup[j].replace(k, 2, bss + "bar{b}");}
     if (latex_subgroup[j][k] == 't' && latex_subgroup[j][k+1] == '~'){latex_subgroup[j].replace(k, 2, bss + "bar{t}");}
     if (latex_subgroup[j][k] == 'w' && latex_subgroup[j][k+1] == 'p'){latex_subgroup[j].replace(k, 2, "W^+");}
     if (latex_subgroup[j][k] == 'e' && latex_subgroup[j][k+1] == 'p'){latex_subgroup[j].replace(k, 2, bss + "e^+");}
     if (latex_subgroup[j][k] == 'w' && latex_subgroup[j][k+1] == 'm'){latex_subgroup[j].replace(k, 2, "W^-");}
     if (latex_subgroup[j][k] == '_'){latex_subgroup[j].replace(k, 1, bss + "rightarrow ");}
    }
    latex_subgroup[j] = '$' + latex_subgroup[j] + '$';
  }
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

string output_result_deviation(double & av_result, double & av_deviation, int digits_deviation){
  double log_result = log10(abs(av_result));//!!!
  double log_deviation = log10(av_deviation);
  int ilog_result = int(log_result);
  if (log_result < 0.){ilog_result--;}
  int ilog_deviation = int(log_deviation);
  if (log_deviation < 0.){ilog_deviation--;}
  int precision = ilog_result - ilog_deviation + digits_deviation;
  double deviation2digits = av_deviation / pow(10.,ilog_deviation + 1 - digits_deviation);
  if (deviation2digits >= pow(10.,digits_deviation) - 0.5){deviation2digits = pow(10.,digits_deviation - 1); precision--;}
  stringstream out;
  out << setw(10) << setprecision(precision) << showpoint << av_result << "(" << setw(digits_deviation) << setprecision(digits_deviation) << noshowpoint << deviation2digits << ")";
  return out.str();
}
string output_percent(double & av_result, double & result, int digits_deviation){
  double deltapercent = (result / av_result - 1.) * 100.;
  double log_dpc = log10(abs(deltapercent));
  int ilog_dpc = int(log_dpc);
  if (log_dpc < 0.){ilog_dpc--;}
  int precision = digits_deviation + ilog_dpc + 1;
  if (result == av_result){precision = digits_deviation;}
  stringstream out;
  //  int precision = 3;
  out << setw(11 - 2) << setprecision(precision) << showpoint << showpos << deltapercent << " %";
  return out.str();
}

string output_latex_percent(double & av_result, double & result, int digits_deviation){
  double deltapercent = (result / av_result - 1.) * 100.;
  double log_dpc = log10(abs(deltapercent));
  int ilog_dpc = int(log_dpc);
  if (log_dpc < 0.){ilog_dpc--;}
  int precision = digits_deviation + ilog_dpc + 1;
  if (result == av_result){precision = digits_deviation;}
  if (precision < 1){precision = 1;}
  stringstream out;
  //  int precision = 3;
  out << setw(11 - 2) << setprecision(precision) << showpoint << showpos << deltapercent << " " << char(92) << "%";
  return out.str();
}

string output_commadigits(double & result, int digits){
  double log_dpc = log10(abs(result));
  int ilog_dpc = int(log_dpc);
  if (log_dpc < 0.){ilog_dpc--;}
  int precision = digits + ilog_dpc + 1;
  stringstream out;
  //  int precision = 3;
  out << setprecision(precision) << showpoint << showpos << result;
  return out.str();
}





void table_4_start(ofstream & out_res, int is, int scale_number, int columns, vector<string> v_dir_directory){
  //  out_res << char(92) << "begin{table}[p]" << endl;
  out_res << char(92) << "centering" << endl;
  out_res << char(92) << "begin{tabularx}{" << char(92) << "linewidth}{|U}" << endl;
  out_res << char(92) << "hline" << endl;
  for (int j = 0; j < columns; j++){out_res << "&";}
  out_res << "&" << char(92) << char(92) << "[-2.ex]" << endl;
  /*
  for (int j = 0; j < columns; j++){
    for (int k = 0; k < v_dir_directory[j + is * columns].size(); k++){
      if (v_dir_directory[j][k] == '_'){v_dir_directory[j + is * columns].replace(k, 1, "-");}
    }
  }
  */
  out_res << "no. & subprocess ";
  for (int j = 0; j < columns; j++){
    if (j + is * columns < v_dir_directory.size()){out_res << " & " << v_dir_directory[j + is * columns];}
    else{out_res << " & " << endl;}
  }
  out_res << char(92) << char(92) << endl;
  out_res << char(92) << "multicolumn{" << columns + 1 << "}{c}{}" << char(92) << char(92) << "[-4.5ex]" << endl;
  for (int j = 0; j < columns; j++){out_res << "&";}
  out_res << "&" << char(92) << char(92) << endl;
  out_res  << char(92) << "hline" << endl;
  out_res << char(92) << "multicolumn{" << columns + 1 << "}{c}{}" << char(92) << char(92) << "[-4.5ex]" << endl;
  for (int j = 0; j < columns; j++){out_res << "&";}
  out_res << "&" << char(92) << char(92) << endl;
}

void table_4_end(ofstream & out_res){
  out_res << char(92) << "multicolumn{8}{c}{}" << char(92) << char(92) << "[-4.5ex]" << endl;
  for (int j = 0; j < 7; j++){out_res << "&";}
  out_res << "&" << char(92) << char(92) << endl;
  out_res << char(92) << "hline" << endl;
  out_res << char(92) << "end{tabularx}" << endl;
  //  out_res << char(92) << "end{table}" << endl;
}

void write_infile_int(ofstream & out_file, string name, int value_i){
  stringstream tab_name;
  if (name.size() < 8){tab_name << char(9) << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 16){tab_name << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 24){tab_name << char(9) << char(9) << char(9);}
  else if (name.size() < 32){tab_name << char(9) << char(9);}
  else if (name.size() < 40){tab_name << char(9);}
  /*
  if (name.size() < 8){tab_name << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 16){tab_name << char(9) << char(9) << char(9);}
  else if (name.size() < 24){tab_name << char(9) << char(9);}
  else if (name.size() < 32){tab_name << char(9);}
  */
  stringstream value_ss;
  value_ss << value_i;
  string value_s = value_ss.str();
  stringstream tab_value;
  if (value_s.size() < 8){tab_value << char(9) << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 16){tab_value << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 24){tab_value << char(9) << char(9) << char(9);}
  else if (value_s.size() < 32){tab_value << char(9) << char(9);}
  else if (value_s.size() < 40){tab_value << char(9);}
  /*
  if (value_s.size() < 8){tab_value << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 16){tab_value << char(9) << char(9) << char(9);}
  else if (value_s.size() < 24){tab_value << char(9) << char(9);}
  else if (value_s.size() < 32){tab_value << char(9);}
  */
  out_file << name << tab_name.str() << "=" << char(9) << value_s << tab_value.str() << "%" << endl;
}

void write_infile_long_long(ofstream & out_file, string name, long long value_i){
  stringstream tab_name;
  if (name.size() < 8){tab_name << char(9) << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 16){tab_name << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 24){tab_name << char(9) << char(9) << char(9);}
  else if (name.size() < 32){tab_name << char(9) << char(9);}
  else if (name.size() < 40){tab_name << char(9);}
  /*
  if (name.size() < 8){tab_name << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 16){tab_name << char(9) << char(9) << char(9);}
  else if (name.size() < 24){tab_name << char(9) << char(9);}
  else if (name.size() < 32){tab_name << char(9);}
  */
  stringstream value_ss;
  value_ss << value_i;
  string value_s = value_ss.str();
  stringstream tab_value;
  if (value_s.size() < 8){tab_value << char(9) << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 16){tab_value << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 24){tab_value << char(9) << char(9) << char(9);}
  else if (value_s.size() < 32){tab_value << char(9) << char(9);}
  else if (value_s.size() < 40){tab_value << char(9);}
  /*
  if (value_s.size() < 8){tab_value << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 16){tab_value << char(9) << char(9) << char(9);}
  else if (value_s.size() < 24){tab_value << char(9) << char(9);}
  else if (value_s.size() < 32){tab_value << char(9);}
  */
  out_file << name << tab_name.str() << "=" << char(9) << value_s << tab_value.str() << "%" << endl;
}

void write_infile_double(ofstream & out_file, string name, double value_d){
  stringstream tab_name;
  if (name.size() < 8){tab_name << char(9) << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 16){tab_name << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 24){tab_name << char(9) << char(9) << char(9);}
  else if (name.size() < 32){tab_name << char(9) << char(9);}
  else if (name.size() < 40){tab_name << char(9);}
  /*
  if (name.size() < 8){tab_name << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 16){tab_name << char(9) << char(9) << char(9);}
  else if (name.size() < 24){tab_name << char(9) << char(9);}
  else if (name.size() < 32){tab_name << char(9);}
  */
  stringstream value_ss;
  value_ss << setw(8) << value_d;
  string value_s = value_ss.str();
  stringstream tab_value;
  if (value_s.size() < 8){tab_value << char(9) << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 16){tab_value << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 24){tab_value << char(9) << char(9) << char(9);}
  else if (value_s.size() < 32){tab_value << char(9) << char(9);}
  else if (value_s.size() < 40){tab_value << char(9);}
  /*
  if (value_s.size() < 8){tab_value << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 16){tab_value << char(9) << char(9) << char(9);}
  else if (value_s.size() < 24){tab_value << char(9) << char(9);}
  else if (value_s.size() < 32){tab_value << char(9);}
  */
  out_file << name << tab_name.str() << "=" << char(9) << value_s << tab_value.str() << "%" << endl;
}

void write_infile_string(ofstream & out_file, string name, string value_s){
  stringstream tab_name;
  if (name.size() < 8){tab_name << char(9) << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 16){tab_name << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 24){tab_name << char(9) << char(9) << char(9);}
  else if (name.size() < 32){tab_name << char(9) << char(9);}
  else if (name.size() < 40){tab_name << char(9);}
  /*
  stringstream tab_value;
  if (value_s.size() < 8){tab_value << char(9) << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 16){tab_value << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 24){tab_value << char(9) << char(9) << char(9);}
  else if (value_s.size() < 32){tab_value << char(9) << char(9);}
  else if (value_s.size() < 40){tab_value << char(9);}
  */
  out_file << name << tab_name.str() << "=" << char(9) << value_s << endl;
  // << tab_value.str() << "%"
}


void old_write_infile_string(ofstream & out_file, string name, string value_s){
  stringstream tab_name;
  if (name.size() < 8){tab_name << char(9) << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 16){tab_name << char(9) << char(9) << char(9) << char(9);}
  else if (name.size() < 24){tab_name << char(9) << char(9) << char(9);}
  else if (name.size() < 32){tab_name << char(9) << char(9);}
  else if (name.size() < 40){tab_name << char(9);}
  stringstream tab_value;
  if (value_s.size() < 8){tab_value << char(9) << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 16){tab_value << char(9) << char(9) << char(9) << char(9);}
  else if (value_s.size() < 24){tab_value << char(9) << char(9) << char(9);}
  else if (value_s.size() < 32){tab_value << char(9) << char(9);}
  else if (value_s.size() < 40){tab_value << char(9);}
  out_file << name << tab_name.str() << "=" << char(9) << value_s << tab_value.str() << "%" << endl;
}




//#include "OV.routines.output.cpp"
