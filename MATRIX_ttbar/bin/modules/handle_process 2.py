


import sys
import os
import shutil
import glob
import readline
import tarfile
import subprocess
import time
from os.path import join as pjoin
from sys import platform as _platform
from initialize_classes_matrix import out, inp, bcolors, tail

#{{{ def: make_tarfile(output_filename, source_dir)
def make_tarfile(output_filename, source_dir):
  try:
      tar = tarfile.open(output_filename, "w:gz")
      tar.add(source_dir, arcname=os.path.basename(source_dir))
  except:
    out.print_error("Could not create tarball (output_filename="+output_filename+", sourcedir="+source_dir+")")
  else:
    tar.close()
#}}}
#{{{ class: CustomCompleter(object)
class CustomCompleter(object):  # Creates custom completer
#{{{ def: __init__(self, options)
    def __init__(self, options):
        self.options = sorted(options)
#}}}
#{{{ def: complete(self, text, state)
    def complete(self, text, state):
        if state == 0:  # on first trigger, build possible matches
            if text:  # cache matches (entries that start with entered text)
                self.matches = [s for s in self.options 
                                    if s and s.startswith(text)]
            else:  # no text entered, all matches possible
                self.matches = self.options[:]

        # return match indexed by state
        try: 
            return self.matches[state]
        except IndexError:
            return None
#}}}
#}}}
#{{{ class: select_process()
class select_process():
#{{{ def: __init__(self)
    def __init__(self):
        # this class handles the naming of the run
        self.available_processes = self.get_available_processes()
        self.process_order = []
        self.process_description = {}
#}}}
#{{{ def: get_available_processes(self)
    def get_available_processes(self):
        # get list of processes that are available in ./prc/ folder
        return [x.split("/",1)[1] for x in glob.glob("prc/pp*") if x[-1].isdigit() or x.endswith("nockm") or x.endswith("EW") or x.endswith("NLOgg")]
#}}}
#{{{ def: sort_available_processes(self)
    def sort_available_processes(self):
        # get list of processes that are available in ./prc/ folder
        unsorted_processes = self.available_processes
        sorted_processes = []
        for process in self.process_order:
            if process in unsorted_processes:
                sorted_processes.append(process)
                unsorted_processes.remove(process)
        sorted_processes = sorted_processes + unsorted_processes 
        self.available_processes = sorted_processes
#}}}
#{{{ def: list_processes(self)
    def list_processes(self):
        # get longest process_id and process by looping through available processes
        self.sort_available_processes()
        length_id = 0
        length_process = 0
        length_description = 0
        for process_name in self.available_processes:
            if process_name == "": continue
            length_id = max(length_id,len(process_name)+7)
            length_process = max(length_process,len(self.convert_to_nice_process_name(process_name))+7)
            length_description = max(length_description,len(self.process_description.get(process_name,"no description available, PROCESS POSSIBLY NOT SUPPORTED")))
        print "-"*(length_id+length_process+length_description)
        print "process_id"+" "*(length_id-len("process_id")-5)+"||   process"+" "*(length_process-len("process")-5)+"||   description"
        print "-"*(length_id+length_process+length_description)
#        for process_name in sorted(self.available_processes, key=lambda x: int(x.replace("hms","").replace("wms","").replace("zms","")[-1])): # sorted by QCD order
        for process_name in self.available_processes:
            print process_name+" "*(length_id-len(process_name)-5)+">>   "+self.convert_to_nice_process_name(process_name)+" "*(length_process-len(self.convert_to_nice_process_name(process_name))-5)+">>   "+self.process_description.get(process_name,"no description available, PROCESS POSSIBLY NOT SUPPORTED")
#}}}
#{{{ def: readin_process(self)
    def readin_process(self): # reads in the process as user input
        # use available_processes list for completion in readin
        completer = CustomCompleter(self.available_processes)
        readline.set_completer(completer.complete)
        readline.parse_and_bind('tab: complete')

        # readin process to be compiled/created
        input = raw_input("|============>> ")
#                         "                "  
        return input
#}}}
#{{{ def: check_process(self,folder)
    def check_process(self,process): # checks user defined process
        check_successful = True
        # process not in available_process list 
        if process not in self.available_processes:
            return False
        return check_successful
#}}}
#{{{ def: convert_to_nice_process_name(self,string)
    def convert_to_nice_process_name(self,string):
        # pre-define replacements so that output string is in nice format
        replacements = []
        # first remove the numbers from the string
        string = ''.join([i for i in string if not i.isdigit()])
        # careful the order is very important!
        replacements.append(["pp","p p --> "]) # easy one
        replacements.append(["heft","(HEFT)"])
        replacements.append(["zms","ZMS"])     # on-shell Z
        replacements.append(["wms","WMS"])     # on-shell W
        replacements.append(["hms","HMS"])     # on-shell Higgs
        replacements.append(["z","Z "])        # easy one
        replacements.append(["h","H "])        # easy one
        replacements.append(["e","e^- "])      # tricky
        replacements.append(["e^- x","e^+ "])  # tricky
        replacements.append(["m","mu^- "])      # tricky
        replacements.append(["mu^- x","mu^+ "])  # tricky
        replacements.append(["w","W^- "])        
        replacements.append(["W^- x","W^+ "])  
        replacements.append(["ZMS"," (Z on-shell) "]) # on-shell Z
        replacements.append(["WMS"," (W on-shell) "]) # on-shell W
        replacements.append(["HMS"," (Higgs on-shell) "]) # on-shell Higgs
        replacements.append(["a","gamma "])    # easy one
        replacements.append(["t","top "])      # tricky
        replacements.append(["top x","anti-top"])  # tricky
        replacements.append(["y","tau^- "])      # tricky
        replacements.append(["tau^- x","tau^+ "])  # tricky
        replacements.append(["nmu^+","v_mu^+"])  # recover the neutrinos
        replacements.append(["nmu^-","v_mu^-"])  # recover the neutrinos
        replacements.append(["ne^+", "v_e^+"])  # recover the neutrinos
        replacements.append(["ne^-", "v_e^-"])  # recover the neutrinos
        replacements.append(["nockmu^-"," (no CKM)"]) # CKM
        replacements.append(["j","j "])

        for item in replacements:
            find    = item[0]
            replace = item[1]
            string = string.replace(find,replace)
        return string
#}}}    
#}}}
#{{{ class: licence_agreement()
class licence_agreement():
#{{{ def: __init__(self,licence_order,licence_dict,agreed)
    def __init__(self,licence_order,licence_dict,agreed):
        # this class handles the licence agreement before the compilation
        self.licence_order = licence_order
        self.licence_dict = licence_dict
        self.agreed = agreed
        self.start_agreement()
        out.print_make("You have agreed with all MATRIX usage terms.")
#}}}
#{{{ def: start_agreement(self)
    def start_agreement(self):
        for licence in self.licence_order:
            if licence == "standard":
                out.print_make("MATRIX is based on several computations, studies and tools from various people and groups. When using results obtained by MATRIX these efforts must be acknowledged by citing the list of references in the CITATION.bib file, which is created in the result folder with every run.")
            else:
                link    = self.licence_dict[licence][0]
                authors = self.licence_dict[licence][1]
                arxiv   = self.licence_dict[licence][2]
                out.print_make("This compilation of MATRIX uses directly the code %s from %s. You have to cite %s from %s, when using results obtained with this installation." %(licence, link, arxiv, authors))
            while True:
                if self.agreed:
                    break
                out.print_read("Do you agree with these terms? Type \"y\" to agree, or \"n\" to abort the code.")
                # read user input
                input = raw_input("|============>> ")
                # abort script if not pressed ENTER
                if input.strip() == "y":
                    break
                elif input.strip() == "n":
                    out.print_error("You have not agree with the MATRIX usage terms. Stopping the Code...")
                    return
#}}}
#}}}
#{{{ class: compile_process()
class compile_process():
#{{{ def: __init__(self,process_in,matrix_dir_in,nr_cores_in)
    def __init__(self,process_in,matrix_dir_in,nr_cores_in,openloops_amplitudes_in,do_resum_in):
        # this class handles the compilation of the process and all related dependencies
        if _platform == "darwin": # MAC OS X
            self.on_mac = True
        else:
            self.on_mac = False
        self.process = process_in
        if self.process.endswith("EW"):
            self.OL_extra_string = "_ew" 
        elif self.process.endswith("NLOgg"):
            self.OL_extra_string = "_beta"
        else:
            self.OL_extra_string = "_beta" # now because of loop-induced amplitude calls, must always install OL public beta
        self.matrix_dir = matrix_dir_in
        self.nr_cores = nr_cores_in
        self.path_to_lhapdf = ""
        self.path_to_gsl = ""
        self.path_to_openloops = ""
        self.ginac_dir = ""
        self.cln_dir = ""
	self.chaplin_dir=""
        self.qqvvamp_dir = ""
        self.ggvvamp_dir = ""
        self.ampzz_dir = ""
        self.ampww_dir = ""
        self.libgfortran_dir = ""
        self.clean_process = False # per default set to false, might be changed by flag
        self.openloops_amplitudes = openloops_amplitudes_in # openloops amplitudes to be downloaded+compiled for the given process
        self.do_resum = do_resum_in

    def create_makefile(self):
        if self.do_resum:
            cloud_makefile = pjoin(self.matrix_dir,"Makefile.clean_res")
        else:
            cloud_makefile = pjoin(self.matrix_dir,"Makefile.clean")
        self.makefile = pjoin(self.matrix_dir,"Makefile")
        shutil.copy(cloud_makefile,self.makefile)
        inp.input_change_entry(self.makefile,"HOMEPATH",self.matrix_dir)
        inp.input_change_entry(self.makefile,"LHAPDF_CONFIG",self.path_to_lhapdf)
        inp.input_change_entry(self.makefile,"GSL_CONFIG",self.path_to_gsl)
        inp.input_change_entry(self.makefile,"OpenLoops_CONFIG",self.path_to_openloops)
        inp.input_change_entry(self.makefile,"GINAC_DIR",self.ginac_dir)
        inp.input_change_entry(self.makefile,"CLN_DIR",self.cln_dir)
        inp.input_change_entry(self.makefile,"CHAPLIN_DIR",self.chaplin_dir)
        inp.input_change_entry(self.makefile,"FORTRAN_LIB_PATH",self.libgfortran_dir)
        if self.do_resum:
            inp.input_change_entry(self.makefile,"MOREDIR",self.more_dir)
        inp.input_change_entry(self.makefile,"HOMEPATH",self.matrix_dir)

    def download_and_compile_openloops(self,in_dir, openloops_tar):
        openloops_dir = pjoin(in_dir,"OpenLoops-install"+self.OL_extra_string)

        self.path_to_openloops = pjoin(openloops_dir,"openloops")

        if self.openloops_compiled():
            out.print_make("OpenLoops already downloaded and compiled. Remove folder %s if you want to re-download and re-compile..." % openloops_dir)
            return
#        os.makedirs(openloops_dir)
        os.chdir(in_dir)
        out.print_make("Extracting OpenLoops...")
        tar = tarfile.open(openloops_tar)
        tar.extractall()
        tar.close()
        os.chdir(openloops_dir)
        try:
            os.makedirs(self.openloops_dir)
        except:
            pass

        openloops_out_path = pjoin(self.matrix_dir,"openloops.log")

        # install OpenLoops
        self.create_openloops_cfg(openloops_dir) # create config file with compile_extra = 1, so that reals are compiled as well
        out.print_make("Compiling OpenLoops...")
        with open(openloops_out_path,'a') as make_out:
            make_out.write("OpenLoops compilation (scons):\n\n")
        with open(openloops_out_path,'a') as make_out:
            download = subprocess.Popen(["./scons"], stdout=make_out, stderr=make_out)
            while True:
                if download.poll() != None:
                    break # stops while loop because job finished
        with open(openloops_out_path,'a') as make_out:
            download = subprocess.Popen(["./scons","lt=pptt","generator=0"], stdout=make_out, stderr=make_out)
            while True:
                if download.poll() != None:
                    break # stops while loop because job finished
        with open(openloops_out_path,'a') as make_out:
            download = subprocess.Popen(["./scons","lt=ppttj","generator=0"], stdout=make_out, stderr=make_out)
            while True:
                if download.poll() != None:
                    break # stops while loop because job finished
        with open(openloops_out_path,'a') as make_out:
            download = subprocess.Popen(["./scons","lt=pptt_im","generator=0"], stdout=make_out, stderr=make_out)
            while True:
                if download.poll() != None:
                    break # stops while loop because job finished
        with open(openloops_out_path,'a') as make_out:
            download = subprocess.Popen(["./scons","ls=pptt2","generator=0"], stdout=make_out, stderr=make_out)
            while True:
                if download.poll() != None:
                    break # stops while loop because job finished
        if not self.openloops_compiled():
            out.print_error_no_stop("Compilation failed: OpenLoops library was not created.")
            out.print_error_no_stop("Check \"OpenLoops.log\" file for errors. Last 20 lines of log:")
            print "\n" + bcolors.FAIL + tail( openloops_out_path, 20 )+"\n"
            exit(0)
#}}}
#{{{ def: openloops_compiled(self, path_to_openloops = self.path_to_openloops)
    def openloops_compiled(self, path_to_openloops = None):
        if path_to_openloops is None:
            path_to_openloops = self.path_to_openloops
        openloops_dir = os.path.dirname(self.path_to_openloops)
        compiled = False
        if os.path.isfile(path_to_openloops) and glob.glob(pjoin(openloops_dir,"lib","libopenloops.*")) and glob.glob(pjoin(openloops_dir,"lib","librambo.*")) and glob.glob(pjoin(openloops_dir,"lib","libcollier.*")) and glob.glob(pjoin(openloops_dir,"lib","liboneloop.*")) and glob.glob(pjoin(openloops_dir,"lib","libcuttools.*")) and glob.glob(pjoin(openloops_dir,"lib","libolcommon.*")):
            compiled = True
        return compiled
#}}}
#{{{ def: download_openloops_amplitudes(self)
    def download_openloops_amplitudes(self):
        openloops_dir = os.path.dirname(self.path_to_openloops)
        os.chdir(openloops_dir)
        self.create_openloops_cfg(openloops_dir)
        # download amplitudes with OpenLoops
        openloops_out_path = pjoin(self.matrix_dir,"OpenLoops_amplitude.log")
        open(openloops_out_path,'w').close()
        firsttime = True
        for amplitude in self.openloops_amplitudes:
            if self.openloops_amplitude_compiled(amplitude):
                out.print_make("OpenLoops %s amplitude already downloaded and compiled. Checking wether up-to-date..." % amplitude)
                with open(openloops_out_path,'a') as make_out:
                    if not firsttime:
                        make_out.write("\n\n")
                    make_out.write("OpenLoops %s amplitude download:\n\n" % amplitude)
                with open(openloops_out_path,'a') as make_out:
                    download = subprocess.Popen(["./openloops","libinstall",amplitude], stdout=make_out, stderr=make_out)
                    while True:
                        if download.poll() != None:
                            break # stops while loop because job finished
                skipped = False
                for line in tail( openloops_out_path, 20 ).splitlines(): # loop through last 20 lines of OpenLoops.log file
                    if line.startswith("- process: %s ... skipped" % amplitude):
                        out.print_make("...%s amplitude already installed and up-to-date." % amplitude)
                        skipped = True
                if not skipped:
                    out.print_make("...%s amplitude updated." % amplitude)
            else:
                out.print_make("Downloading and compiling %s amplitude with OpenLoops..." % amplitude)
                with open(openloops_out_path,'a') as make_out:
                    if not firsttime:
                        make_out.write("\n\n")
                    make_out.write("OpenLoops %s amplitude download:\n\n" % amplitude)
                with open(openloops_out_path,'a') as make_out:
                    download = subprocess.Popen(["./openloops","libinstall",amplitude], stdout=make_out, stderr=make_out)
                    while True:
                        if download.poll() != None:
                            break # stops while loop because job finished
            if not self.openloops_amplitude_compiled(amplitude):
                out.print_error_no_stop("Amplitude download/compilation failed: library for OpenLoops %s amplitude was not created." % amplitude)
                out.print_error_no_stop("Check \"OpenLoops_amplitude.log\" file for errors. Last 20 lines of log:")
                print "\n" + bcolors.FAIL + tail( openloops_out_path, 20 )+"\n"
                exit(0)            
            firsttime = False
        if self.on_mac:
            OL_lib_folder = pjoin(openloops_dir,"lib")
            OL_proclib_folder = pjoin(openloops_dir,"proclib")
            out.print_info("Running on Mac. Trying to make relative paths of linked OpenLoops dylibs absolute. Please consider using")
            print "export DYLD_LIBRARY_PATH=DYLD_LIBRARY_PATH:%s:%s" % (OL_lib_folder,OL_proclib_folder)
            out.print_info("in your terminal and possibly adding it to your .bashrc/.bash_profile, in case you still experience linking errors when running the code.")
            self.mac_make_absolute_path_in_dylib_linking(pjoin(self.matrix_dir,"bin",self.process))
            self.mac_make_absolute_path_in_dylib_linking(pjoin(openloops_dir,"lib","libopenloops.dylib"))
            self.mac_make_absolute_path_in_dylib_linking(pjoin(openloops_dir,"lib","libolcommon.dylib"))
            self.mac_make_absolute_path_in_dylib_linking(pjoin(openloops_dir,"lib","liboneloop.dylib"))
            self.mac_make_absolute_path_in_dylib_linking(pjoin(openloops_dir,"lib","libcollier.dylib"))
            self.mac_make_absolute_path_in_dylib_linking(pjoin(openloops_dir,"lib","librambo.dylib"))
            self.mac_make_absolute_path_in_dylib_linking(pjoin(openloops_dir,"lib","libcuttools.dylib"))
            for amplitude in self.openloops_amplitudes:
                try:
                    amp_libfile = glob.glob(pjoin(openloops_dir,"proclib","libopenloops_%s_*.dylib" % amplitude))[0]
                except:
                    amp_libfile = ""
                    pass
                if amp_libfile:
                    self.mac_make_absolute_path_in_dylib_linking(amp_libfile)
#}}}

#{{{ def: mac_make_absolute_path_in_dylib_linking(self,file_path)
    def mac_make_absolute_path_in_dylib_linking(self,file_path):
        otool = subprocess.Popen(["otool","-L",file_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
        linked_dylibs = [x.split()[0] for x in otool.split("\n\t") if "compatibility" in x and not x.startswith("/")]
#        print otool
        for dylib in linked_dylibs:
            openloops_dir = os.path.dirname(self.path_to_openloops)
            abs_path = pjoin(openloops_dir,"lib",dylib.split("/")[1])
            install_name_tool = subprocess.Popen(["install_name_tool","-change",dylib,abs_path,file_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
#            print install_name_tool
#}}}
#{{{ def: create_openloops_cfg(self,path)
    def create_openloops_cfg(self,path):
        cfg_path = pjoin(path,"openloops.cfg")
        if not os.path.isfile(cfg_path):
            with open(cfg_path,'w') as openloops_cfg:
                openloops_cfg.write("[OpenLoops]\n")
                openloops_cfg.write("compile_extra = 1\n")
                if self.on_mac:
                    openloops_cfg.write("link_flags = %(common_flags)s -headerpad_max_install_names\n")
                if self.OL_extra_string == "_ew":
                    openloops_cfg.write("process_repositories = public, public_ew\n")
                elif self.OL_extra_string == "_beta":
                    openloops_cfg.write("process_repositories = powheg, public, public_beta\n")
        else:
            inp.input_set_entry(cfg_path,"compile_extra","1")
            if self.on_mac:
                inp.input_set_entry(cfg_path,"link_flags","%(common_flags)s -headerpad_max_install_names")
            if self.OL_extra_string == "_ew":
                inp.input_set_entry(cfg_path,"process_repositories","public, public_ew")
            elif self.OL_extra_string == "_beta":
                inp.input_set_entry(cfg_path,"process_repositories","powheg, public, public_beta")
#}}}
#{{{ def: openloops_amplitude_compiled(self, amplitude, path_to_openloops = self.path_to_openloops)
    def openloops_amplitude_compiled(self, amplitude, path_to_openloops = None):
        if path_to_openloops is None:
            path_to_openloops = self.path_to_openloops

        compiled = False
        if glob.glob(pjoin(os.path.dirname(path_to_openloops),"proclib","libopenloops_%s_*" % amplitude)) and os.path.exists(pjoin(os.path.dirname(path_to_openloops),"process_src",amplitude)) and os.path.exists(pjoin(os.path.dirname(path_to_openloops),"process_obj",amplitude)):
            compiled = True
        return compiled
#}}}
#{{{ def: compile_MoRe(self,more_tar)
    def compile_MoRe(self,more_tar):
        self.more_dir = pjoin(os.path.dirname(more_tar),"MoRe-v1.0.0")
        if self.MoRe_compiled():
            out.print_make("MoRe already compiled. Remove folder %s if you want to re-compile..." % self.more_dir)
            return
        os.chdir(os.path.dirname(more_tar))
        out.print_make("Extracting and Compiling MoRe from %s into %s..." % (more_tar,self.more_dir))
        tar = tarfile.open(more_tar)
        tar.extractall()
        tar.close()
        os.chdir(self.more_dir)

        more_out_path = pjoin(self.matrix_dir,"MoRe.log")
        with open(more_out_path,'w') as make_out:
            make_out.write("MoRe make:\n\n")
        with open(more_out_path,'a') as make_out:
            make = subprocess.Popen(["make"], stdout=make_out, stderr=make_out)
            while True:
                if make.poll() != None:
                    break # stops while loop because job finished
        if not self.MoRe_compiled():
            out.print_error_no_stop("Compilation failed: MoRe library was not created.")
            out.print_error_no_stop("Check \"more.log\" file for errors. Last 20 lines of log:")
            print "\n" + bcolors.FAIL + tail( more_out_path, 20 )+"\n"
            exit(0)
#}}}
#{{{ def: MoRe_compiled(self, more_dir = self.more_dir)
    def MoRe_compiled(self, more_dir = None):
        if more_dir is None:
            more_dir = self.more_dir

        compiled = False
        if glob.glob(pjoin(more_dir,"libmore.*")):
            compiled = True
        return compiled
#}}}
#{{{ def: compile_qqvvamp(self,qqvvamp_tar)
    def compile_qqvvamp(self,qqvvamp_tar):
    # does actually only the extraction of qqVVamp, compilation is done with main code
        self.qqvvamp_dir = qqvvamp_tar.rsplit(".",1)[0]
        if self.qqvvamp_compiled():
            out.print_make("qqVVamp already extracted. Remove folder %s if you want to extract it again..." % self.qqvvamp_dir)
            return
        os.chdir(os.path.dirname(qqvvamp_tar))
        out.print_make("Extracting qqVVamp from %s into %s..." % (qqvvamp_tar,self.qqvvamp_dir))
        tar = tarfile.open(qqvvamp_tar)
        tar.extractall()
        tar.close()
#}}}
#{{{ def: qqvvamp_compiled(self, qqvvamp_dir = self.qqvvamp_dir)
    def qqvvamp_compiled(self, qqvvamp_dir = None):
    # does actually only test weither Qqvvamp is extracted
        if qqvvamp_dir is None:
            qqvvamp_dir = self.qqvvamp_dir

        compiled = False
        if os.path.exists(qqvvamp_dir):
            compiled = True
        return compiled
#}}}
#{{{ def: compile_ggvvamp(self,ggvvamp_tar)
    def compile_ggvvamp(self,ggvvamp_tar):
    # does actually only the extraction of ggVVamp, compilation is done with main code
        self.ggvvamp_dir = ggvvamp_tar.rsplit(".",1)[0]
        if self.ggvvamp_compiled():
            out.print_make("ggVVamp already extracted. Remove folder %s if you want to extract it again..." % self.ggvvamp_dir)
            return
        os.chdir(os.path.dirname(ggvvamp_tar))
        out.print_make("Extracting ggVVamp from %s into %s..." % (ggvvamp_tar,self.ggvvamp_dir))
        tar = tarfile.open(ggvvamp_tar)
        tar.extractall()
        tar.close()
#}}}
#{{{ def: ggvvamp_compiled(self, ggvvamp_dir = self.ggvvamp_dir)
    def ggvvamp_compiled(self, ggvvamp_dir = None):
    # does actually only test weither Ggvvamp is extracted
        if ggvvamp_dir is None:
            ggvvamp_dir = self.ggvvamp_dir

        compiled = False
        if os.path.exists(ggvvamp_dir):
            compiled = True
        return compiled
#}}}
#{{{ def: compile_ampzz(self,ampzz_tar)
    def compile_ampzz(self,ampzz_tar):
    # does actually only the extraction of Ampzz, compilation is done with main code
        self.ampzz_dir = ampzz_tar.rsplit(".",1)[0]
        if self.ampzz_compiled():
            out.print_make("On-shell ZZ amplitudes already extracted. Remove folder %s if you want to extract it again..." % self.ampzz_dir)
            return
        os.chdir(os.path.dirname(ampzz_tar))
        out.print_make("Extracting on-shell ZZ amplitudes from %s into %s..." % (ampzz_tar,self.ampzz_dir))
        tar = tarfile.open(ampzz_tar)
        tar.extractall()
        tar.close()
#}}}
#{{{ def: ampzz_compiled(self, ampzz_dir = self.ampzz_dir)
    def ampzz_compiled(self, ampzz_dir = None):
    # does actually only test weither Ampzz is extracted
        if ampzz_dir is None:
            ampzz_dir = self.ampzz_dir

        compiled = False
        if os.path.exists(ampzz_dir):
            compiled = True
        return compiled
#}}}
#{{{ def: compile_ampww(self,ampww_tar)
    def compile_ampww(self,ampww_tar):
    # does actually only the extraction of Ampww, compilation is done with main code
        self.ampww_dir = ampww_tar.rsplit(".",1)[0]
        if self.ampww_compiled():
            out.print_make("On-shell WW amplitudes already extracted. Remove folder %s if you want to extract it again..." % self.ampww_dir)
            return
        os.chdir(os.path.dirname(ampww_tar))
        out.print_make("Extracting on-shell WW amplitudes from %s into %s..." % (ampww_tar,self.ampww_dir))
        tar = tarfile.open(ampww_tar)
        tar.extractall()
        tar.close()
#}}}
#{{{ def: ampww_compiled(self, ampww_dir = self.ampww_dir)
    def ampww_compiled(self, ampww_dir = None):
    # does actually only test weither Ampww is extracted
        if ampww_dir is None:
            ampww_dir = self.ampww_dir

        compiled = False
        if os.path.exists(ampww_dir):
            compiled = True
        return compiled
#}}}
#{{{ def: compile_cln(self,cln_tar)
    def compile_cln(self,cln_tar):
        self.cln_dir = pjoin(os.path.dirname(cln_tar),"cln-install")
        if self.cln_compiled():
            out.print_make("Cln already compiled. Remove folder %s if you want to re-compile..." % self.cln_dir)
            return
        os.chdir(os.path.dirname(cln_tar))
        out.print_make("Extracting and Compiling Cln from %s into %s..." % (cln_tar,self.cln_dir))
        tar = tarfile.open(cln_tar)
        tar.extractall()
        tar.close()
        os.chdir(cln_tar.rsplit(".",1)[0])
        try:
            os.makedirs(self.cln_dir)
        except:
            pass

        cln_out_path = pjoin(self.matrix_dir,"cln.log")
        with open(cln_out_path,'w') as make_out:
            make_out.write("Cln configure:\n\n")
        with open(cln_out_path,'a') as make_out:
            configure = subprocess.Popen(["./configure","--prefix=%s" % self.cln_dir], stdout=make_out, stderr=make_out)
            while True:
                if configure.poll() != None:
                    break # stops while loop because job finished
        with open(cln_out_path,'a') as make_out:
            make_out.write("\n\n\n\nCln make:\n\n")
        with open(cln_out_path,'a') as make_out:
            make = subprocess.Popen(["make","-j",str(self.nr_cores)], stdout=make_out, stderr=make_out)
            while True:
                if make.poll() != None:
                    break # stops while loop because job finished
        with open(cln_out_path,'a') as make_out:
            make_out.write("\n\n\n\nCln make install:\n\n")
        with open(cln_out_path,'a') as make_out:
            make = subprocess.Popen(["make","install"], stdout=make_out, stderr=make_out)
            while True:
                if make.poll() != None:
                    break # stops while loop because job finished
        if not self.cln_compiled():
            out.print_error_no_stop("Compilation failed: Cln library was not created.")
            out.print_error_no_stop("Check \"cln.log\" file for errors. Last 20 lines of log:")
            print "\n" + bcolors.FAIL + tail( cln_out_path, 20 )+"\n"
            exit(0)
#}}}
#{{{ def: cln_compiled(self, cln_dir = self.cln_dir)
    def cln_compiled(self, cln_dir = None):
        if cln_dir is None:
            cln_dir = self.cln_dir

        compiled = False
        if (glob.glob(pjoin(cln_dir,"lib","libcln*")) or glob.glob(pjoin(cln_dir,"lib64","libcln*"))) and glob.glob(pjoin(cln_dir,"include","cln","cln.h")):
            compiled = True
        return compiled
#}}}
#{{{ def: compile_ginac(self,ginac_tar)
    def compile_ginac(self,ginac_tar):
        self.ginac_dir = pjoin(os.path.dirname(ginac_tar),"ginac-install")
        if self.ginac_compiled():
            out.print_make("Ginac already compiled. Remove folder %s if you want to re-compile..." % self.ginac_dir)
            return
        os.chdir(os.path.dirname(ginac_tar))
        out.print_make("Extracting and Compiling Ginac from %s into %s..." % (ginac_tar,self.ginac_dir))
        tar = tarfile.open(ginac_tar)
        tar.extractall()
        tar.close()
        os.chdir(ginac_tar.rsplit(".",1)[0])
        try:
            os.makedirs(self.ginac_dir)
        except:
            pass

        ginac_out_path = pjoin(self.matrix_dir,"ginac.log")

        my_env = os.environ.copy()
        my_env["CLN_CFLAGS"] = "-I"+pjoin(self.cln_dir,"include")
        try:
            my_env["CLN_LIBS"] = glob.glob(pjoin(self.cln_dir,"lib*","libcln*"))[0]
        except:
            out.print_error("Could not find lbcln* library in lib folder of cln install folder: %s" % self.cln_dir)
#        os.system("export CLN_CFLAGS=%s" % self.cln_dir)
#        os.system("export CLN_LIBS=%s" % self.cln_dir)
        with open(ginac_out_path,'w') as make_out:
            make_out.write("Ginac configure:\n\n")
        with open(ginac_out_path,'a') as make_out:
            configure = subprocess.Popen(["./configure","--prefix=%s" % self.ginac_dir,"CPPFLAGS=-I%s" % pjoin(self.cln_dir,"include")], stdout=make_out, stderr=make_out, env=my_env)
            while True:
                if configure.poll() != None:
                    break # stops while loop because job finished
        with open(ginac_out_path,'a') as make_out:
            make_out.write("\n\n\n\nGinac make:\n\n")
        with open(ginac_out_path,'a') as make_out:
            make = subprocess.Popen(["make","-j",str(self.nr_cores)], stdout=make_out, stderr=make_out, env=my_env)
            while True:
                if make.poll() != None:
                    break # stops while loop because job finished
        with open(ginac_out_path,'a') as make_out:
            make_out.write("\n\n\n\nGinac make install:\n\n")
        with open(ginac_out_path,'a') as make_out:
            make = subprocess.Popen(["make","install"], stdout=make_out, stderr=make_out, env=my_env)
            while True:
                if make.poll() != None:
                    break # stops while loop because job finished
        # check again if process compiled now, if not print end of ginac.log on screen
        if not self.ginac_compiled():
            out.print_error_no_stop("Compilation failed: Ginac library was not created.")
            out.print_error_no_stop("Check \"ginac.log\" file for errors. Last 20 lines of log:")
            print "\n" + bcolors.FAIL + tail( ginac_out_path, 20 )+"\n"
            exit(0)
#}}}
#{{{ def: ginac_compiled(self, ginac_dir = self.ginac_dir)
    def ginac_compiled(self, ginac_dir = None):
        if ginac_dir is None:
            ginac_dir = self.ginac_dir
        compiled = False
        if (glob.glob(pjoin(ginac_dir,"lib","libginac*")) or glob.glob(pjoin(ginac_dir,"lib64","libginac*"))) and glob.glob(pjoin(ginac_dir,"include","ginac","ginac.h")):
            compiled = True
        return compiled
#}}}

    def compile_chaplin(self,chaplin_tar):
        self.chaplin_dir = pjoin(os.path.dirname(chaplin_tar),"chaplin-install")
        if self.chaplin_compiled():
            out.print_make("Chaplin already compiled. Remove folder %s if you want to re-compile..." % self.chaplin_dir)
            return
        os.chdir(os.path.dirname(chaplin_tar))
        out.print_make("Extracting and Compiling Chaplin from %s into %s..." % (chaplin_tar,self.chaplin_dir))
        tar = tarfile.open(chaplin_tar)
        tar.extractall()
        tar.close()
        os.chdir(chaplin_tar.rsplit(".",1)[0])
        try:
            os.makedirs(self.chaplin_dir)
        except:
            pass

        chaplin_out_path = pjoin(self.matrix_dir,"chaplin.log")

        with open(chaplin_out_path,'w') as make_out:
            make_out.write("Chaplin configure:\n\n")
        with open(chaplin_out_path,'a') as make_out:
            configure = subprocess.Popen(["./configure","--prefix=%s" % self.chaplin_dir], stdout=make_out, stderr=make_out)
            while True:
                if configure.poll() != None:
                    break # stops while loop because job finished

        with open(chaplin_out_path,'a') as make_out:
            make_out.write("\n\n\n\nChaplin make install:\n\n")
        with open(chaplin_out_path,'a') as make_out:
            make = subprocess.Popen(["make","install"], stdout=make_out, stderr=make_out)
            while True:
                if make.poll() != None:
                    break # stops while loop because job finished
        if not self.chaplin_compiled():
            out.print_error_no_stop("Compilation failed: Chaplin library was not created.")
            out.print_error_no_stop("Check \"chaplin.log\" file for errors. Last 20 lines of log:")
            print "\n" + bcolors.FAIL + tail( chaplin_out_path, 20 )+"\n"
            exit(0)

# moving the lib to lib64
        if os.path.exists(pjoin(self.chaplin_dir,"lib")):
            shutil.move(pjoin(self.chaplin_dir,"lib"),pjoin(self.chaplin_dir,"lib64"))

    def chaplin_compiled(self, chaplin_dir = None):
        if chaplin_dir is None:
            chaplin_dir = self.chaplin_dir
        compiled = False
        if (glob.glob(pjoin(chaplin_dir,"lib","libchaplin*")) or glob.glob(pjoin(chaplin_dir,"lib64","libchaplin*"))):
            compiled = True
        return compiled


#{{{ def: compile_ginac(self,ginac_tar)
    # def compile_ginac(self,ginac_tar):
    #     self.ginac_dir = pjoin(os.path.dirname(ginac_tar),"ginac-install")
    #     if self.ginac_compiled():
    #         out.print_make("Ginac already compiled. Remove folder %s if you want to re-compile..." % self.ginac_dir)
    #         return
    #     os.chdir(os.path.dirname(ginac_tar))
    #     out.print_make("Extracting and Compiling Ginac from %s into %s..." % (ginac_tar,pjoin(os.path.dirname(ginac_tar),"ginac-install")))
    #     tar = tarfile.open(ginac_tar)
    #     tar.extractall()
    #     tar.close()
    #     os.chdir(ginac_tar.rsplit(".",1)[0])
    #     try:
    #         os.makedirs(self.ginac_dir)
    #     except:
    #         pass

    #     ginac_out_path = pjoin(self.matrix_dir,"ginac.log")
    #     with open(ginac_out_path,'w') as make_out:
    #         make_out.write("Ginac configure:\n\n")
    #     with open(cln_out_path,'a') as make_out:
    #         configure = subprocess.Popen(["./configure","--prefix=%s" % self.ginac_dir], stdout=make_out, stderr=make_out)
    #         while True:
    #             if configure.poll() != None:
    #                 break # stops while loop because job finished
    #     with open(ginac_out_path,'a') as make_out:
    #         make_out.write("\n\n\n\nGinac make:\n\n")
    #     with open(cln_out_path,'a') as make_out:
    #         make = subprocess.Popen(["make","-j",str(self.nr_cores)], stdout=make_out, stderr=make_out)
    #         while True:
    #             if make.poll() != None:
    #                 break # stops while loop because job finished
    #     with open(ginac_out_path,'a') as make_out:
    #         make_out.write("\n\n\n\nGinac make install:\n\n")
    #     with open(cln_out_path,'a') as make_out:
    #         make = subprocess.Popen(["make","install"], stdout=make_out, stderr=make_out)
    #         while True:
    #             if make.poll() != None:
    #                 break # stops while loop because job finished
    #     exit(0)
#}}}
#{{{ def: compile_matrix_process(self)
    def compile_matrix_process(self):
        if self.clean_process:
            self.clean_matrix_process()
        elif self.process_compiled():
            out.print_make("Process <%s> already compiled." % self.process)
            while True:
                out.print_read("Do you want to re-compile? Type \"y\" to re-compile, or \"n\" to continue without.")
                # read user input
                input = raw_input("|============>> ")
                # abort script if not pressed ENTER
                if input.strip() == "y":
                    self.clean_matrix_process()
                    break
                elif input.strip() == "n":
                    out.print_make("Continuing without re-compililation...")
                    return
        out.print_make("Compiling process <%s>, this may take a while...            (see make.log file to monitor the progress)" % self.process)
        # first unpack relevant phase space files if not yet unpacked
        if not os.path.exists(pjoin(self.matrix_dir,"psg",self.process)):
            os.chdir(pjoin(self.matrix_dir,"psg"))
            tar = tarfile.open("psg."+self.process+".tar")
            tar.extractall()
            tar.close()
        os.chdir(self.matrix_dir)
        make_out_path = pjoin(self.matrix_dir,"make.log")
        with open(make_out_path,'w+') as make_out:
            make = subprocess.Popen(["make","-j",str(self.nr_cores),"-f",self.makefile,self.process], stdout=make_out, stderr=make_out)
            while True:
                if make.poll() != None:
                    break # stops while loop because job finished
        # check again if process compiled now, if not print end of make.log on screen
        if not self.process_compiled():
            out.print_error_no_stop("Compilation failed: executable in \"bin/%s\" does not exist." % self.process)
            out.print_error_no_stop("Check \"make.log\" file for errors. Last 20 lines of log:")
            print "\n" + bcolors.FAIL + tail( make_out_path, 20 ) + bcolors.ENDC + "\n"
            exit(0)
#}}}
#{{{ def: clean_matrix_process(self)
    def clean_matrix_process(self):
        out.print_make("Cleaning process <%s>, this may take a while..." % self.process)
        os.chdir(self.matrix_dir)
        make_out_path = pjoin(self.matrix_dir,"clean.log")
        exe_path = pjoin(self.matrix_dir,"bin",self.process)
        try:
            os.remove(exe_path)
        except:
            pass
        # can be removed ??? --->>
        # twoloopID = self.openloops_amplitudes[0]
        # # remove releveant parts of 2-loop amplitudes
        # if twoloopID in ["ppllll","ppzz","ppww"]: # in these cases do it by hand
        #     if twoloopID == "ppzz":
        #         twploopID = "ppzz22" # otherwisethis we remove more than just the 2-loop object files for ppzz
        #     # don't remove the library, because it may effect running processes; removing the object file should be 
        #     # sufficient to force the rebuild
        #     # try:
        #     #     remove_files = glob.glob("lib/%s*/lib%s*" % (twoloopID,twoloopID))
        #     #     for remove_file in remove_files:
        #     #         os.remove(remove_file)
        #     # except:
        #     #     pass
        #     try:
        #         remove_files = glob.glob("obj/%s*/%s*" % (twoloopID,twoloopID))
        #         for remove_file in remove_files:
        #             os.remove(remove_file)
        #     except:
        #         pass
        #     clean_by_makefile = False
        # else: # otherwise use cleaning routine of Makefile
        #     clean_by_makefile = True
        # <<--- can be removed ???
        with open(make_out_path,'w+') as make_out:
#            make_out.write("make clean log:\n")
            make = subprocess.Popen(["make","-j",str(self.nr_cores),"-f",self.makefile,"clean"], stdout=make_out, stderr=make_out)
            while True:
                if make.poll() != None:
                    break # stops while loop because job finished
#            make_out.write("\n\nmake %s log:\n" % (self.process+"-clean-prc"))
            make = subprocess.Popen(["make","-j",str(self.nr_cores),"-f",self.makefile,self.process+"-clean-prc"], stdout=make_out, stderr=make_out)
            while True:
                if make.poll() != None:
                    break # stops while loop because job finished
#            make_out.write("\n\nmake %s log:\n" % (self.process+"-clean-psg"))
            make = subprocess.Popen(["make","-j",str(self.nr_cores),"-f",self.makefile,self.process+"-clean-psg"], stdout=make_out, stderr=make_out)
            while True:
                if make.poll() != None:
                    break # stops while loop because job finished
#            if clean_by_makefile:
            make = subprocess.Popen(["make","-j",str(self.nr_cores),"-f",self.makefile,self.process+"-clean-amp"], stdout=make_out, stderr=make_out)
            while True:
                if make.poll() != None:
                    break # stops while loop because job finished
#}}}
#{{{ def: process_compiled(self)
    def process_compiled(self):
        compiled = False
        exe_path = pjoin(self.matrix_dir,"bin",self.process)
        if os.path.isfile(exe_path): compiled = True
        return compiled
#}}}
#}}}
#{{{ class: process_folder()
class process_folder():
#{{{ def: __init__(self,process_in)
    def __init__(self,process_in,matrix_dir_in,standalone_in,folder_name_extension):
        # this class handles the naming of the run
        self.process = process_in
        self.matrix_dir = matrix_dir_in
        self.matrix_configuration_file = pjoin(self.matrix_dir,"config","MATRIX_configuration")
        self.matrix_run_process_script = pjoin(self.matrix_dir,"bin","run_process")
        self.process_executable   = pjoin(self.matrix_dir,"bin",self.process)
        self.setup_file_parameter = pjoin(self.matrix_dir,"run","setup","file_parameter.dat")
        self.process_folder_path = pjoin(self.matrix_dir,"run",self.process+folder_name_extension)
        self.process_default_folder_path = pjoin(self.matrix_dir,"run",self.process)
        self.process_default_folder_tar  = pjoin(self.matrix_dir,"run","run."+self.process+".tar")
        self.input_folder = pjoin(self.matrix_dir,"run","input_files",self.process)
        self.standalone = standalone_in # in standalone mode copy everything instead of creating links
#}}}
#{{{ def: create_process_folder(self)
    def create_process_folder(self):
        if os.path.exists(self.process_folder_path):
            out.print_warning("Process folder \"%s\" already exists. If you want to create it new, (re)move this folder and try again. Skipping process folder creation..." % self.process_folder_path)
        elif os.path.exists(self.process_default_folder_path):
            out.print_error("The default process folder \"%s\" already exists. This folder should not exist before creating a new MATRIX process folder. You can (re)move this folder and try again. Exiting..." % self.process_default_folder_path)
        else:
            os.chdir(pjoin(self.matrix_dir,"run"))
            tar = tarfile.open(self.process_default_folder_tar)
            tar.extractall()
            tar.close()
            shutil.move(self.process_default_folder_path,self.process_folder_path)
            os.chdir(self.process_folder_path)
            os.makedirs("bin")
            os.makedirs("input")
            if self.standalone:
                os.symlink(self.matrix_run_process_script,"bin/run_process")
                shutil.copy(self.matrix_configuration_file,"input/MATRIX_configuration")
            else: # default
                os.symlink(self.matrix_run_process_script,"bin/run_process")
                os.symlink(self.matrix_configuration_file,"input/MATRIX_configuration")
            # the following can be removed once the tar file have the right content >>>
            try:
                shutil.rmtree(pjoin(self.process_folder_path,"default.grid.final"))
            except:
                pass
            try:
                shutil.rmtree(pjoin(self.process_folder_path,"default.MUNICH"))
            except:
                pass
            try:
                shutil.rmtree(pjoin(self.process_folder_path,"batch"))
            except:
                pass
            # shutil.copy(pjoin(self.matrix_dir,"run",self.process+"_script","default.grid.final.tar"),pjoin(self.process_folder_path,"default.grid.final.tar"))
            # tar = tarfile.open(pjoin(self.process_folder_path,"default.grid.final.tar"))
            # tar.extractall()
            # tar.close()
            
            # copy MATRIX inputs
            shutil.copytree(pjoin(self.input_folder,"default.input.MATRIX"),pjoin(self.process_folder_path,"input","default.input.MATRIX"),symlinks=True)
            # # create internal folder that should not be touched by the user
            # os.makedirs(pjoin(self.process_folder_path,"bin","internal"))
            # copy MUNICH inputs directly in default MATRIX folder
            shutil.copy(pjoin(self.input_folder,"file_parameter.dat"),pjoin("default.MATRIX","file_parameter.dat"))
            shutil.copy(pjoin(self.input_folder,"file_model.dat"),pjoin("default.MATRIX","file_model.dat"))
            shutil.copy(pjoin(self.input_folder,"file_distribution.dat"),pjoin("default.MATRIX","file_distribution.dat"))
            os.remove(pjoin("default.MATRIX",self.process))
            os.symlink(self.process_executable,pjoin("default.MATRIX",self.process)) # we don't have to copy the executable, the script takes care of that
                                                                                     # reply: for debugging it makes sense to have it there anyways...
            os.remove(pjoin("default.MATRIX","setup","file_parameter.dat"))
            if self.standalone:
                shutil.copy(self.setup_file_parameter,pjoin("default.MATRIX","setup","file_parameter.dat"))
            else: # default
                os.symlink(self.setup_file_parameter,pjoin("default.MATRIX","setup","file_parameter.dat"))
            make_tarfile("default.MATRIX.tar","default.MATRIX")
            shutil.rmtree("default.MATRIX")
            out.print_info("Process folder successfully created.")
#}}}
#}}}
