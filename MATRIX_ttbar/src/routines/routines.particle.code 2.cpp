#include "../include/classes.cxx"

void fill_code_particle(map<string, int> & code_particle){
  code_particle["d"]   = 1;
  code_particle["u"]   = 2;
  code_particle["s"]   = 3;
  code_particle["c"]   = 4;
  code_particle["b"]   = 5;
  code_particle["t"]   = 6;
  code_particle["d~"]  = -1;
  code_particle["u~"]  = -2;
  code_particle["s~"]  = -3;
  code_particle["c~"]  = -4;
  code_particle["b~"]  = -5;
  code_particle["t~"]  = -6;
  code_particle["em"]  = 11;
  code_particle["mum"] = 13;
  code_particle["tam"] = 15;
  code_particle["ep"]  = -11;
  code_particle["mup"] = -13;
  code_particle["tap"] = -15;
  code_particle["ve"]  = 12;
  code_particle["vm"]  = 14;
  code_particle["vt"]  = 16;
  code_particle["ve~"] = -12;
  code_particle["vm~"] = -14;
  code_particle["vt~"] = -16;
  code_particle["g"]   = 0;
  code_particle["a"]   = 22;
  code_particle["wp"]  = -24;
  code_particle["wm"]  = 24;
  code_particle["z"]   = 23;
  code_particle["h"]   = 25;

  code_particle["q"]  = 10;
  code_particle["p"]  = 101;
  code_particle["p~"]  = -101;
  code_particle["jet"]  = 100;
}



void fill_name_particle(map<int, string> & name_particle){
  name_particle[  0] = "g";
  name_particle[ 22] = "a";
  name_particle[-24] = "wp";
  name_particle[ 24] = "wm";
  name_particle[ 23] = "z";
  name_particle[ 25] = "h";
  name_particle[  1] = "d";
  name_particle[  2] = "u";
  name_particle[  3] = "s";
  name_particle[  4] = "c";
  name_particle[  5] = "b";
  name_particle[  6] = "t";
  name_particle[ -1] = "d~";
  name_particle[ -2] = "u~";
  name_particle[ -3] = "s~";
  name_particle[ -4] = "c~";
  name_particle[ -5] = "b~";
  name_particle[ -6] = "t~";
  name_particle[ 11] = "em";
  name_particle[ 12] = "ve";
  name_particle[ 13] = "mum";
  name_particle[ 14] = "vm";
  name_particle[ 15] = "tam";
  name_particle[ 16] = "vt";
  name_particle[-11] = "ep";
  name_particle[-12] = "ve~";
  name_particle[-13] = "mup";
  name_particle[-14] = "vm~";
  name_particle[-15] = "tap";
  name_particle[-16] = "vt~";

  name_particle[10] = "q";
  name_particle[101] = "p";
  name_particle[-101] = "p~";
  name_particle[100] = "jet";
}


/*
void fill_pname(map<int, string> & pname){
  pname[10] = "q";
  pname[0] = "g";
  pname[22] = "A";
  pname[-24] = "W+";
  pname[24] = "W-";
  pname[23] = "Z";
  pname[25] = "h";
  pname[1] = "d";
  pname[2] = "u";
  pname[3] = "s";
  pname[4] = "c";
  pname[5] = "b";
  pname[6] = "t";
  pname[-1] = "d~";
  pname[-2] = "u~";
  pname[-3] = "s~";
  pname[-4] = "c~";
  pname[-5] = "b~";
  pname[-6] = "t~";
  pname[11] = "e-";
  pname[12] = "ve";
  pname[13] = "mu-";
  pname[14] = "vm";
  pname[15] = "ta-";
  pname[16] = "vt";
  pname[-11] = "e+";
  pname[-12] = "ve~";
  pname[-13] = "mu+";
  pname[-14] = "vm~";
  pname[-15] = "ta+";
  pname[-16] = "vt~";
}
*/


void fill_datname(map<int, string> & datname){
  datname[0] = "g";
  datname[22] = "a";
  datname[-24] = "wp";
  datname[24] = "wm";
  datname[23] = "z";
  datname[25] = "h";
  datname[1] = "d";
  datname[2] = "u";
  datname[3] = "s";
  datname[4] = "c";
  datname[5] = "b";
  datname[6] = "t";
  datname[11] = "em";
  datname[12] = "ve";
  datname[13] = "mum";
  datname[14] = "vm";
  datname[15] = "tam";
  datname[16] = "vt";
  datname[-1] = "d~";
  datname[-2] = "u~";
  datname[-3] = "s~";
  datname[-4] = "c~";
  datname[-5] = "b~";
  datname[-6] = "t~";
  datname[-11] = "ep";
  datname[-12] = "ve~";
  datname[-13] = "mup";
  datname[-14] = "vm~";
  datname[-15] = "tap";
  datname[-16] = "vt~";
}


/*
void fill_pmadgraph(map<int, string> & pmadgraph){
  pmadgraph[0] = "G";
  pmadgraph[1] = "D";
  pmadgraph[2] = "U";
  pmadgraph[3] = "S";
  pmadgraph[4] = "C";
  pmadgraph[5] = "B";
  pmadgraph[6] = "T";
  pmadgraph[11] = "EM";
  pmadgraph[12] = "VE";
  pmadgraph[13] = "MUM";
  pmadgraph[14] = "VM";
  pmadgraph[15] = "TAUM";
  pmadgraph[16] = "VT";
  pmadgraph[-1] = "DB";
  pmadgraph[-2] = "UB";
  pmadgraph[-3] = "SB";
  pmadgraph[-4] = "CB";
  pmadgraph[-5] = "BB";
  pmadgraph[-6] = "TB";
  pmadgraph[-11] = "EP";
  pmadgraph[-12] = "VEB";
  pmadgraph[-13] = "MUP";
  pmadgraph[-14] = "VMB";
  pmadgraph[-15] = "TAUP";
  pmadgraph[-16] = "VTB";
  pmadgraph[22] = "A";
  pmadgraph[-24] = "WP";
  pmadgraph[24] = "WM";
  pmadgraph[23] = "Z";
  pmadgraph[25] = "H";
}
*/


void fill_charge_particle(map<int, double> & charge_particle){
  charge_particle[1] = -1. / 3.;
  charge_particle[2] = 2. / 3.;
  charge_particle[3] = -1. / 3.;
  charge_particle[4] = 2. / 3.;
  charge_particle[5] = -1. / 3.;
  charge_particle[6] = 2. / 3.;

  charge_particle[-1] = 1. / 3.;
  charge_particle[-2] = -2. / 3.;
  charge_particle[-3] = 1. / 3.;
  charge_particle[-4] = -2. / 3.;
  charge_particle[-5] = 1. / 3.;
  charge_particle[-6] = -2. / 3.;
  
  charge_particle[11] = -1.;
  charge_particle[12] = 0.;
  charge_particle[13] = -1.;
  charge_particle[14] = 0.;
  charge_particle[15] = -1.;
  charge_particle[16] = 0.;

  charge_particle[-11] = 1.;
  charge_particle[-12] = 0.;
  charge_particle[-13] = 1.;
  charge_particle[-14] = 0.;
  charge_particle[-15] = 1.;
  charge_particle[-16] = 0.;
  
  charge_particle[0] = 0.;
  charge_particle[21] = 0.;
  charge_particle[22] = 0.;
  charge_particle[23] = 0.;
  charge_particle[25] = 0.;

  charge_particle[24] = -1.;
  charge_particle[-24] = 1.;
}



void fill_latexname(map<int, string> & latexname){
  char bs = char(92);
  string bss;
  bss.push_back(bs);
  latexname[0] = "g";
  latexname[22] = bss + "gamma";
  latexname[-24] = "W^+";
  latexname[24] = "W^-";
  latexname[23] = "Z";
  latexname[25] = "H";
  latexname[1] = "d";
  latexname[2] = "u";
  latexname[3] = "s";
  latexname[4] = "c";
  latexname[5] = "b";
  latexname[6] = "t";
  latexname[11] = "e";
  latexname[12] = bss + "nu_e";
  latexname[13] = bss + "mu";
  latexname[14] = bss + "nu_" + bss + "mu";
  latexname[15] = bss + "tau";
  latexname[16] =  bss + "nu_" + bss + "tau";
  latexname[-1] = bss + "bar{d}";
  latexname[-2] = bss + "bar{u}";
  latexname[-3] = bss + "bar{s}";
  latexname[-4] = bss + "bar{c}";
  latexname[-5] = bss + "bar{b}";
  latexname[-6] = bss + "bar{t}";
  latexname[-11] = bss + "bar{e}";
  latexname[-12] = bss + "bar{" + bss + "nu_e}";
  latexname[-13] = bss + "bar{" + bss + "mu}";
  latexname[-14] = bss + "bar{" + bss + "nu_" + bss + "mu}";
  latexname[-15] = bss + "bar{" + bss + "tau}";
  latexname[-16] = bss + "bar{" + bss + "nu_" + bss + "tau}";
}



void fill_gnuplotname(map<int, string> & gnuplotname){
  gnuplotname[0] = "&mathrm{g}";
  gnuplotname[22] = "&gamma";
  gnuplotname[-24] = "&mathrm{W^+}";
  gnuplotname[24] = "&mathrm{W^-}";
  gnuplotname[23] = "&mathrm{Z}";
  gnuplotname[25] = "&mathrm{H}";
  gnuplotname[1] = "&mathrm{d}";
  gnuplotname[2] = "&mathrm{u}";
  gnuplotname[3] = "&mathrm{s}";
  gnuplotname[4] = "&mathrm{c}";
  gnuplotname[5] = "&mathrm{b}";
  gnuplotname[6] = "&mathrm{t}";
  gnuplotname[11] = "&mathrm{e^-}";
  gnuplotname[12] = "&mathrm{&nu_e}";
  gnuplotname[13] = "&mathrm{&mu^-}";
  gnuplotname[14] = "&mathrm{&nu_&mu}";
  gnuplotname[15] = "&mathrm{&tau^-}";
  gnuplotname[16] = "&mathrm{&nu_&tau}";
  gnuplotname[-1] = "&mathrm{&bar{d}}";
  gnuplotname[-2] = "&mathrm{&bar{u}}";
  gnuplotname[-3] = "&mathrm{&bar{s}}";
  gnuplotname[-4] = "&mathrm{&bar{c}}";
  gnuplotname[-5] = "&mathrm{&bar{b}}";
  gnuplotname[-6] = "&mathrm{&bar{t}}";
  gnuplotname[-11] = "&mathrm{e^+}";
  gnuplotname[-12] = "&mathrm{&bar{&nu}_e}";
  gnuplotname[-13] = "&mathrm{&mu^+}";
  gnuplotname[-14] = "&mathrm{&bar{&nu}_&mu}";
  gnuplotname[-15] = "&mathrm{&tau^+}";
  gnuplotname[-16] = "&mathrm{&bar{&nu}_&tau}";
}



void determine_process(string & subprocess, vector<vector<string> > & particles_name, int & process_type, int & n_particle, vector<int> & pa, map<string, int> & code_particle){
  Logger logger("process readin");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  int particle_counter = 0;
  int exit = 0;
  logger << LOG_DEBUG_VERBOSE << "begin 'determine_process' for " << subprocess << endl;
  subprocess = subprocess + "     ";
  for (int i = 0; i < subprocess.size(); i++){
    logger << LOG_DEBUG_VERBOSE << i << "   " << subprocess[i] << "   " << particle_counter << endl;
    if (subprocess[i] == '_' || subprocess[i] == '-'){particle_counter++;}
    else if (subprocess[i] == 'p'){
      if (subprocess[i + 1] == '~'){particles_name[particle_counter].push_back("p~"); i++;}
      else {particles_name[particle_counter].push_back("p");}
    }
    else if (subprocess[i] == 'j'){
      if (subprocess[i + 1] == 'e' && subprocess[i + 2] == 't'){particles_name[particle_counter].push_back("jet"); i += 2;}
      else {exit = 1; break;}
    }
    else if (subprocess[i] == '2'){
      if (subprocess[i + 1] == 'j' && subprocess[i + 2] == 'e' && subprocess[i + 3] == 't' && subprocess[i + 4] == 's'){
	for (int j = 0; j < 2; j++){particles_name[particle_counter].push_back("jet");}
	i += 4;
      }
      else {exit = 1; break;}
    }
    else if (subprocess[i] == '3'){
      if (subprocess[i + 1] == 'j' && subprocess[i + 2] == 'e' && subprocess[i + 3] == 't' && subprocess[i + 4] == 's'){
	for (int j = 0; j < 3; j++){particles_name[particle_counter].push_back("jet");}
	i += 4;
      }
      else {exit = 1; break;}
    }
    else if (subprocess[i] == '4'){
      if (subprocess[i + 1] == 'j' && subprocess[i + 2] == 'e' && subprocess[i + 3] == 't' && subprocess[i + 4] == 's'){
	for (int j = 0; j < 4; j++){particles_name[particle_counter].push_back("jet");}
	i += 4;
      }
      else {exit = 1; break;}
    }
    else if (subprocess[i] == 'f'){
      if (subprocess[i + 1] == '~'){particles_name[particle_counter].push_back("f~"); i++;}
      else {particles_name[particle_counter].push_back("f");}
    }

    else if (subprocess[i] == 'd'){
      if (subprocess[i + 1] == '~'){particles_name[particle_counter].push_back("d~"); i++;}
      else {particles_name[particle_counter].push_back("d");}
    }
    else if (subprocess[i] == 'u'){
      if (subprocess[i + 1] == '~'){particles_name[particle_counter].push_back("u~"); i++;}
      else {particles_name[particle_counter].push_back("u");}
    }
    else if (subprocess[i] == 's'){
      if (subprocess[i + 1] == '~'){particles_name[particle_counter].push_back("s~"); i++;}
      else {particles_name[particle_counter].push_back("s");}
    }
    else if (subprocess[i] == 'c'){
      if (subprocess[i + 1] == '~'){particles_name[particle_counter].push_back("c~"); i++;}
      else {particles_name[particle_counter].push_back("c");}
    }
    else if (subprocess[i] == 'b'){
      if (subprocess[i + 1] == '~'){particles_name[particle_counter].push_back("b~"); i++;}
      else {particles_name[particle_counter].push_back("b");}
    }
    else if (subprocess[i] == 't'){
      if (subprocess[i + 1] == '~'){particles_name[particle_counter].push_back("t~"); i++;}
      else if (subprocess[i + 1] == 'a')
	if      (subprocess[i + 2] == 'm'){particles_name[particle_counter].push_back("tam"); i += 2;}
	else if (subprocess[i + 2] == 'p'){particles_name[particle_counter].push_back("tap"); i += 2;}
	else {exit = 1; break;}
      else {particles_name[particle_counter].push_back("t");}
    }
    else if (subprocess[i] == 'e'){
      if      (subprocess[i + 1] == 'm'){particles_name[particle_counter].push_back("em"); i++;}
      else if (subprocess[i + 1] == 'p'){particles_name[particle_counter].push_back("ep"); i++;}
      else {exit = 1; break;}
    }
    else if (subprocess[i] == 'm'){
      if      (subprocess[i + 1] == 'u' && subprocess[i + 2] == 'm'){particles_name[particle_counter].push_back("mum"); i += 2;}
      else if (subprocess[i + 1] == 'u' && subprocess[i + 2] == 'p'){particles_name[particle_counter].push_back("mup"); i += 2;}
      else {exit = 1; break;}
    }
    else if (subprocess[i] == 'v'){
      if      (subprocess[i + 1] == 'e'){
	if (subprocess[i + 2] == '~'){particles_name[particle_counter].push_back("ve~"); i += 2;}
	else {particles_name[particle_counter].push_back("ve"); i++;}
      }
      else if (subprocess[i + 1] == 'm'){
	if (subprocess[i + 2] == '~'){particles_name[particle_counter].push_back("vm~"); i += 2;}
	else {particles_name[particle_counter].push_back("vm"); i++;}
      }
      else if (subprocess[i + 1] == 't'){
	if (subprocess[i + 2] == '~'){particles_name[particle_counter].push_back("vt~"); i += 2;}
	else {particles_name[particle_counter].push_back("vt"); i++;}
      }
      else {exit = 1; break;}
    }
    else if (subprocess[i] == 'g'){particles_name[particle_counter].push_back("g");}
    else if (subprocess[i] == 'w'){
      if      (subprocess[i + 1] == 'm'){particles_name[particle_counter].push_back("wm"); i++;}
      else if (subprocess[i + 1] == 'p'){particles_name[particle_counter].push_back("wp"); i++;}
      else {exit = 1; break;}
    }
    else if (subprocess[i] == 'z'){particles_name[particle_counter].push_back("z");}
    else if (subprocess[i] == 'h'){particles_name[particle_counter].push_back("h");}
    else if (subprocess[i] == 'a'){particles_name[particle_counter].push_back("a");}
    else if (subprocess[i] == ' '){}
    else if (subprocess[i] == '+'){}
    else if (subprocess[i] == 'X'){}
    else {exit = 1; break;}
  }
  if (exit){} // ??? use of exit ???

  logger << LOG_DEBUG_VERBOSE << "particles_name.size() = " << particles_name.size() << endl;
  n_particle = particles_name[1].size();
  logger << LOG_DEBUG_VERBOSE << "n_particle = " << n_particle << endl;
    for (int i = 0; i < particles_name.size(); i++){
      for (int j = 0; j < particles_name[i].size(); j++){
	logger << LOG_DEBUG_VERBOSE << "particles_name[" << i << "][" << j << "] = " << particles_name[i][j] << endl;
      }
    }
  logger << LOG_DEBUG_VERBOSE << "particles_name[0].size() = " << particles_name[0].size() << endl;
  if (particles_name[0].size() == 1){}
  else if (particles_name[0].size() == 2){pa.push_back(0);}
  logger << LOG_DEBUG_VERBOSE << "particles_name[0].size() = " << particles_name[0].size() << endl;
  //  else {exit = 1; break;}
  process_type = particles_name[0].size();
  logger << LOG_DEBUG_VERBOSE << "process_type = " << process_type << endl;
  for (int i = 0; i < 2; i++){
    logger << LOG_DEBUG_VERBOSE << "particles_name[" << i << "].size() = " << particles_name[i].size() << endl;
    for (int j = 0; j < particles_name[i].size(); j++){
    logger << LOG_DEBUG_VERBOSE << i << "   " << j << endl;
      pa.push_back(code_particle[particles_name[i][j]]);
      logger << LOG_DEBUG_VERBOSE << "particles_name[" << right << setw(2) << i << "][" << setw(2) << j << "] = " << left << setw(5) << particles_name[i][j] << "   code_particle = " << right << setw(4) << code_particle[particles_name[i][j]] << endl;
    }
  }
  if (particles_name[0].size() == 1){for (int p = 0; p < pa.size(); p++){logger << LOG_DEBUG_VERBOSE << "pa[" << right << setw(2) << p << "] = " << right << pa[p] << endl;}}
  else if (particles_name[0].size() == 2){for (int p = 1; p < pa.size(); p++){logger << LOG_DEBUG_VERBOSE << "pa[" << right << setw(2) << p << "] = " << right << pa[p] << endl;}}
  subprocess.erase(subprocess.end() - 5, subprocess.end());
  logger << LOG_DEBUG_VERBOSE << "end 'determine_process' for " << subprocess << endl;
}



void subprocess_readin(string & process_class, string & subprocess, vector<string> & decay, int & process_type, int & n_particle, vector<vector<int> > & type_parton){
  Logger logger("subprocess_readin");
  logger << LOG_DEBUG_VERBOSE << "called" << endl;

  map<string,int> code_particle;
  fill_code_particle(code_particle);
  map<int,string> name_particle;
  fill_name_particle(name_particle);

  logger << LOG_DEBUG << "Partonic process:" << endl;
  vector<vector<string> > particles_name(2);
  determine_process(subprocess, particles_name, process_type, n_particle, type_parton[0], code_particle);
  logger << LOG_DEBUG << "n_particle = " << n_particle << endl;
  for (int i_p = 1; i_p < type_parton[0].size(); i_p++){logger << LOG_DEBUG << "type_parton[0][" << setw(2) << i_p << "] = " << setw(5) << right << type_parton[0][i_p] << "   " << setw(5) << left << name_particle[type_parton[0][i_p]] << endl;}

  int dummy;
  logger << LOG_DEBUG << "Hadronic process:" << endl;
  int hadron_n_particle;
  vector<vector<int> > hadron_type_parton(1);
  vector<vector<string> > hadron_particles_name(2);
  logger << LOG_DEBUG << "process_class = " << process_class << endl;

  determine_process(process_class, hadron_particles_name, dummy, hadron_n_particle, hadron_type_parton[0], code_particle);
  logger << LOG_DEBUG << "hadron_n_particle = " << hadron_n_particle << endl;

  for (int i_p = 1; i_p < hadron_type_parton[0].size(); i_p++){logger << LOG_DEBUG << "hadron_type_parton[0][" << setw(2) << i_p << "] = " << right << setw(5) << hadron_type_parton[0][i_p] << "   " << setw(5) << left << name_particle[hadron_type_parton[0][i_p]] << endl;}

  logger << LOG_DEBUG << process_class << endl;
  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG << "process_type = " << process_type << endl;

  //  logger << LOG_DEBUG << "process file" << endl;
  logger << LOG_DEBUG << "process_class           = " << process_class << endl;
  logger << LOG_DEBUG << "n_particle              = " << n_particle << endl;
  //  cout << "process_number          = " << process_number << endl;
  if (process_type == 2){
    //  if (type_parton[0][0] == 0){
    for (int i2 = 1; i2 < n_particle + 3; i2++){
      logger << LOG_DEBUG << "type_parton[0][" << i2 << "]                  = " << type_parton[0][i2] << endl;
    }
  }
  else if (process_type == 1){
    for (int i2 = 0; i2 < n_particle + 1; i2++){
      logger << LOG_DEBUG << "type_parton[0][" << i2 << "]                  = " << type_parton[0][i2] << endl;
    }
  }
  logger.newLine(LOG_DEBUG);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



