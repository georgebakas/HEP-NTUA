#include "../include/classes.cxx"

string get_path(){
  static Logger logger("get_path");
  logger << LOG_DEBUG << "called" << endl;

  //  char arg1[20];
  char exepath[PATH_MAX + 1] = {0};
  
  //        sprintf(arg1, "/proc/%d/exe", getpid());
  //	cout << "arg1 = " << arg1 << endl;
  //        readlink(arg1, exepath, 1024);
  //	cout << "exepath = " << exepath << endl;
  
  int check = readlink("/proc/self/exe", exepath, 1024);
  logger << LOG_DEBUG_VERBOSE << "check = " << check << endl;
  logger << LOG_DEBUG_VERBOSE << "exepath = " << exepath << endl;
  //	logger << LOG_DEBUG_VERBOSE << "check = " << check << endl;
  string s_exepath = string(exepath);
  string path_MUNICH;
  int counter = 0;
  for (int i_s = s_exepath.size() - 1; i_s >=0; i_s--){
    if (s_exepath[i_s] == '/'){counter++;}
    if (counter == 2){
      path_MUNICH = s_exepath.substr(0, i_s + 1);
      break;
    }
  }
  
  logger << LOG_DEBUG << "finished - " << path_MUNICH << endl;
  
  return path_MUNICH;
}

void system_execute(Logger & logger, string xorder){
  logger << LOG_DEBUG << "xorder = " << xorder << endl;

  int isystem = 0;
  system_execute(logger, xorder, isystem);

  logger << LOG_DEBUG_VERBOSE << "xorder = " << xorder << "   execution status: " << isystem << "." << endl;
}



void system_execute(Logger & logger, string xorder, int isystem){
  logger << LOG_DEBUG << "xorder = " << xorder << endl;

  string order = xorder.substr(0, 5); 
  if (order == "mkdir"){
    DIR *pDir;
    string directory = xorder.substr(6, xorder.size() - 6); 
    logger << LOG_DEBUG_VERBOSE << "xorder = " << xorder << endl;
    logger << LOG_DEBUG_VERBOSE << "directory = " << directory << endl;
    pDir = opendir (directory.c_str());
    if (pDir == NULL) {
      logger << LOG_DEBUG << "pDir is NULL." << endl;
      vector<string> subdirectory;
      for (int i_c = 0; i_c < directory.size(); i_c++){
	if (directory[i_c] == '/'){subdirectory.push_back(directory.substr(0, i_c));}
      }
      subdirectory.push_back(directory);
      for (int i_s = 0; i_s < subdirectory.size(); i_s++){
	logger << LOG_DEBUG << "subdirectory[" << i_s << "] = " << subdirectory[i_s] << endl;
	pDir = opendir (subdirectory[i_s].c_str());
	if (pDir == NULL) {
	  xorder = "mkdir " + subdirectory[i_s];
	  isystem = system(xorder.c_str());
	}
      }
    }
    else {
      logger << LOG_DEBUG << "main pDir is not NULL." << endl;
    }
    closedir (pDir);
  }
  else {
    isystem = system(xorder.c_str());
  }

  logger << LOG_DEBUG << "xorder = " << xorder << "   execution status: " << isystem << "." << endl;
}



int intpow(int basis, int exponent){
  if (exponent < 0){return 0;}
  else {
    int result = 1;
    for (int i = 0; i < exponent; i++){
      result = result * basis;
    }
    return result;
  }
}



double g_from_alpha(double alpha){
  return sqrt(4 * pi * alpha);
}



double lambda(double x, double y, double z){
  double t;
  if      (x == 0.){t = pow(y - z, 2);}
  else if (y == 0.){t = pow(x - z, 2);}
  else if (z == 0.){t = pow(y - x, 2);}
  else {t = pow(x, 2) + pow(y, 2) + pow(z, 2) - 2 * x * y - 2 * x * z - 2 * y * z;}
  return t;
}



string double2hexastr(double d){
  char buffer[25] = { 0 };
  ::snprintf(buffer, 25, "%A", d); // TODO Check for errors
  //  ::snprintf(buffer, 21, "%A", d); // TODO Check for errors
  return buffer;
}



double hexastr2double(const string & s){
  double d = 0.0;
  ::sscanf(s.c_str(), "%lA", &d); // TODO Check for errors
  return d;
}



char * stch(string temp_s){
  char * temp_ch = new char[temp_s.size() + 1];
  std::copy(temp_s.begin(), temp_s.end(), temp_ch);
  temp_ch[temp_s.size()] = '\0';
  //  cout << "stch: " << temp_s << "   " << temp_ch << endl;
  return temp_ch;
}



vector<int> get_vector_from_array_int(int * temp_array, int size_array){
  //  cout << "sizeof(temp_array) = " << sizeof(temp_array) << endl;
  //  cout << "sizeof(int) = " << sizeof(int) << endl;
  //  cout << "size_array = " << size_array << endl;
  vector<int> temp_vector(temp_array, temp_array + size_array);
  //  cout << "temp_vector.size() = " << temp_vector.size() << endl;
  return temp_vector;
}



string time_hms_from_double(double & time){
  stringstream time_ss;
  long long temp_d = (long long)(time / 24 / 3600);
  long long temp_h = (long long)(time / 3600) - temp_d * 24;
  long long temp_min = (long long)(time / 60) - temp_d * 24 * 60 - temp_h * 60;
  long long temp_sec = (long long)(time) - temp_d * 24 * 3600 - temp_h * 3600 - temp_min * 60;
  time_ss << setw(4) << temp_d << " d " << setw(2) << temp_h << " h " << setw(2) << temp_min << " m " << setw(2) << temp_sec << " s";
  return time_ss.str();
}



double max_subset_from_binary(vector<double> & xbsqrtsmin_opt, int b){
  vector<int> temp = vectorbinary_from_binary(b);
  vector<int> xbcomb(0);
  vector<vector<int> > xcomb(0);
  vector<vector<int> > ycomb(0);
  for (int i = 4; i < b / 2; i += 4){
    vector<int> temp_part = vectorbinary_from_binary(i);
    int temp_counter = 0;
    for (int j = 0; j < temp_part.size(); j++){
      for (int k = 0; k < temp.size(); k++){
	if (temp_part[j] == temp[k]){
	  temp_counter++;
	  break;
	}
      }
    }
    if (temp_counter < temp_part.size()){continue;}
    vector<int> temp_part2 = vectorbinary_from_binary(b - i);
    xcomb.push_back(temp_part);
    ycomb.push_back(temp_part2);
    xbcomb.push_back(i);
  }
  for (int i = 0; i < xbcomb.size(); i++){
    double temp = xbsqrtsmin_opt[accumulate(xcomb[i].begin(), xcomb[i].end(), 0)] + xbsqrtsmin_opt[accumulate(ycomb[i].begin(), ycomb[i].end(), 0)];
    //    cout << "temp = " << setw(25) << temp << "   xbsqrtsmin_opt[" << setw(3) << b << "] = " << setw(25) << xbsqrtsmin_opt[b] << endl;
    if (xbsqrtsmin_opt[b] < temp){xbsqrtsmin_opt[b] = temp;}
  }
  //  cout << "max_subset_from_binary end" << endl;
  return xbsqrtsmin_opt[b];
}

vector<int> vectorint_from_binary(int b){
  vector<int> temp;
  int counter = 1;
  while (b != 0){
    if ((b % 2) == 1){temp.push_back(counter);}
    counter++;
    b = b / 2;
  }
  return temp;
}

vector<int> vectorbinary_from_binary(int b){
  vector<int> temp;
  int counter = 0;
  while (b != 0){
    if ((b % 2) == 1){temp.push_back(intpow(2, counter));}
    counter++;
    b = b / 2;
  }
  return temp;
}

double ran(vector<double> & s){
  s[0] = fmod(s[0] + s[1] + s[2], 1.);
  s[1] = fmod(s[0] + s[1] + s[2], 1.);
  s[2] = fmod(s[0] + s[1] + s[2], 1.);
  return s[0];
}

void randomvector(vector<double> & s, int n, vector<double> & r){
  for (int i = 0; i < n + 1; i++){
    r[i] = ran(s);
  }
}

int sign(int x){
  if (x > 0){return 1;}
  else if (x < 0){return -1;}
  else {return 0;}
}

bool munich_isnan(double & temp){
  //  return isnan(temp);
  return std::isnan(temp);
}

bool munich_isinf(double & temp){
  //  return isinf(temp);
  return std::isinf(temp);
}



void get_userinput_from_readin(vector<string> & user_variable, vector<string> & user_variable_additional, vector<string> & user_value, vector<string> & readin){
  Logger logger("get_userinput_from_readin");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  string readindata;
  for (int i = 0; i < readin.size(); i++){
    logger << LOG_DEBUG << "readin[" << setw(4) << i << "].size() = " << readin[i].size() << endl;
    readindata = readin[i][0];
    if (readindata != "/" && readindata != "#" && readindata != "%"){
      int start = 0;
      user_variable.push_back("");
      user_variable_additional.push_back("");
      user_value.push_back("");
      for (int j = 0; j < readin[i].size(); j++){
	if (start == 0 || start == 1){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 0){}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    user_variable[user_variable.size() - 1].push_back(readin[i][j]);
	    if (start != 1){start = 1;}
	  }
	  else {start++;}
	}
	else if (start == 2){
	  if (readin[i][j] == '='){start = 5;}
	  else if ((readin[i][j] == ' ') || (readin[i][j] == char(9))){}
	  else if (start == 2){
	    start++;
	    j--;
	  }
	  else {
	    logger << LOG_ERROR << "Incorrect input in line " << i << endl;
	    user_variable.erase(user_variable.end(), user_variable.end());
	    user_variable_additional.erase(user_variable_additional.end(), user_variable_additional.end());
	    user_value.erase(user_value.end(), user_value.end());
	    break;
	  }
	}
	else if (start == 3){
	  if ((readin[i][j] != '=')){
	    user_variable_additional[user_variable_additional.size() - 1].push_back(readin[i][j]);
	  }
	  else {start++; j--;} // should be the same as shifting (start == 4) here !!!
	}
	else if (start == 4){
	  // additional: should be allowed to contain ' ' !!!
	  if (readin[i][j] == '='){
	    start = 5;
	    logger << LOG_DEBUG << "before: ---" << user_variable_additional[user_variable_additional.size() - 1] << "---" << endl;
	    for (int i_s = user_variable_additional[user_variable_additional.size() - 1].size() - 1; i_s > 0; i_s--){
	      if (user_variable_additional[user_variable_additional.size() - 1][i_s] == ' ' ||
		  user_variable_additional[user_variable_additional.size() - 1][i_s] == char(9)){
		user_variable_additional[user_variable_additional.size() - 1].erase(user_variable_additional[user_variable_additional.size() - 1].end() - 1, user_variable_additional[user_variable_additional.size() - 1].end());
	      }
	      else {break;}
	    }
	    logger << LOG_DEBUG << "after:  ---" << user_variable_additional[user_variable_additional.size() - 1] << "---" << endl;
	  }
	  //	  else if ((readin[i][j] == ' ') || (readin[i][j] == char(9))){}
	  else {
	    logger << LOG_ERROR << "Incorrect input in line " << i << endl;
	    user_variable.erase(user_variable.end(), user_variable.end());
	    user_variable_additional.erase(user_variable_additional.end(), user_variable_additional.end());
	    user_value.erase(user_value.end(), user_value.end());
	    break;
	  }
	}
	else if (start == 5 || start == 6){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 5){}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    user_value[user_value.size() - 1].push_back(readin[i][j]);
	    if (start != 6){start = 6;}
	  }
	  else {start++;}
	}

	else if (start == 5 || start == 6 || start == 7){
	  // start == 5: before beginning of user_value
	  // start == 6: user_value has started
	  // start == 7: 
	  logger << LOG_DEBUG_VERBOSE << "readin[" << i << "][" << j << "] = " << readin[i][j] << endl;
	  if (readin[i][j] == '%'){break;}
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 5){}
	  else if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start > 5){start = 7;}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    if (start == 7){user_value[user_value.size() - 1].push_back(' ');}
	    user_value[user_value.size() - 1].push_back(readin[i][j]);
	    start = 6;
	  }
	  else {start++;}
	}

	/*
	else if (start == 5 || start == 6){
	  if (((readin[i][j] == ' ') || (readin[i][j] == char(9))) && start == 5){}
	  else if ((readin[i][j] != ' ') && (readin[i][j] != char(9))){
	    user_value[user_value.size() - 1].push_back(readin[i][j]);
	    if (start != 6){start = 6;}
	  }
	  else {start++;}
	}
	*/

	else {break;}
      }
    }
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

