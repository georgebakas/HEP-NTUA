#{{{ imports
import os
import shutil
import glob
import tarfile
import time
from os.path import join as pjoin
# own modules
from initialize_classes import out, run_name
#}}}
#{{{ which(program)
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
#}}}
#{{{ class: folder_structure()
class folder_structure(): # class that takes care of the matrix folder structure and creates
#{{{ def: __init__(self,mode,grid_folder,main_run_folder,NLO_subtraction,order,set_parallel_runs,grid_assignment,include_loop_induced,config_list)
    def __init__(self):
        # initialize all the relevant folders for the matrix code
        self.run_folder_path    = ""  # this is the main folder with the MUNICH code, where the runs are done
        self.input_folder_path  = ""  # MATRIX inputs
        self.log_folder_path    = ""  # MATRIX logs
        self.result_folder_path = ""  # MATRIX results
        self.exe_path           = ""  # MATRIX process executable
        self.default_input_path = ""  # default folder of the MATRIX inputs
        self.input_file_dir     = ""  # default folder for the MUNICH inputs to be copied in the main folder (file_parameter/model/distribution.dat)
        # try to create the needed folder structures
        try:
            os.makedirs(result_folder_path)
        except:
            pass
        try:
            os.makedirs(log_folder_path)
        except:
            pass
        # hard-coded lists:
        # a clean run dir must contain the following files and folders
        self.contained_in_clean_run_dir = ["file_parameter.dat","log","clean.contribution"]
        # hack until EW executable corrected
        self.exe_pathQCD        = ""  # MATRIX process executable for EW corrections
        self.exe_pathEW         = ""  # MATRIX process executable for QCD corrections
#}}}
#{{{ def: check_default_input_path(self)
    def check_default_input_path(self):
        # this checks wether the self.default_input_path is suitable (folder exists and contains parameter/model/distribution.dat
        if not os.path.exists(self.default_input_path):
            out.print_error("Path %s for default inputs does not exist. Either the default folder structure is broken or you have chosen an \"--input dir\" inside the input folder that does not exist. Try to recover the default.input.MATRIX folder or choose a proper default input_dir and restart. Exiting..." % self.default_input_path)
        for input_file in ["parameter.dat","model.dat","distribution.dat"]:
            if not os.path.isfile(pjoin(self.default_input_path,input_file)):
                out.print_error("Cannot find default input file %s under path %s for default inputs. Either the default folder structure is broken or you have chosen an \"--input dir\" inside the input folder that does contain the %s file. Try to add the a proper file to the respective folder or choose a proper default input folder and restart. Exiting..." % (input_file,self.default_input_path,input_file))
#}}}
#{{{ def: remove_run(self)
    def remove_run(self):
        # this deletes the run folder including its logs, inputs and results
        out.print_info("Deleting run folder (%s) its inputs, logs and results." % self.run_folder_path)
        try:
            shutil.rmtree(self.run_folder_path)
        except:
            pass
        try:
            shutil.rmtree(self.input_folder_path)
        except:
            pass
        try:
            shutil.rmtree(self.log_folder_path)
        except:
            pass
        try:
            shutil.rmtree(self.result_folder_path)
        except:
            pass
        out.print_info("Exiting...")
#}}}
#{{{ def: tar_run(self)
    def tar_run(self):
        # this tar's the run folder including its log, input and result folder
        out.print_info("Creating .tar archive of run folder (%s) including its input, log and result folder." % self.run_folder_path)
        with tarfile.open(self.run_folder_path+".tar", "w:gz") as tar:
            tar.add(self.run_folder_path, arcname=os.path.basename(self.run_folder_path))    
            tar.add(self.log_folder_path, arcname=pjoin("log",os.path.basename(self.log_folder_path)))    
            tar.add(self.result_folder_path, arcname=pjoin("result",os.path.basename(self.result_folder_path)))
            tar.add(self.input_folder_path, arcname=pjoin("input",os.path.basename(self.input_folder_path)))
        out.print_info("Exiting...")
#}}}
#{{{ def: rename_run_to(self,new_name)
    def rename_run_to(self,new_name):
        # this renames the run folder to new_name
        try:
            shutil.move(self.run_folder_path,pjoin(os.path.dirname(self.run_folder_path),new_name))
        except:
            out.print_error("Could not move the run folder to a new name.")
        try:
            shutil.move(self.log_folder_path,pjoin(os.path.dirname(self.log_folder_path),new_name))
        except:
            out.print_error("Could not move the log of run to a new name.")
        try:
            shutil.move(self.result_folder_path,pjoin(os.path.dirname(self.result_folder_path),new_name))
        except:
            out.print_error("Could not move the result of run to a new name.")
        try:
            shutil.move(self.input_folder_path,pjoin(os.path.dirname(self.input_folder_path),new_name))
        except:
            out.print_error("Could not move the input of run to a new name.")
#}}}
#{{{ def: copy_run_from(self,existing_run)
    def copy_run_from(self,existing_run):
        # this copies a run folder from an existing run
        try:
            shutil.copytree(pjoin(os.path.dirname(self.run_folder_path),existing_run),self.run_folder_path,symlinks=True)
        except:
            out.print_error("Could not copy the run folder from the existing run.")
        try:
            shutil.copytree(pjoin(os.path.dirname(self.log_folder_path),existing_run),self.log_folder_path,symlinks=True)
        except:
            out.print_error("Could not copy the log of the existing run.")
        try:
            os.remove(pjoin(self.log_folder_path,"main.running"))
        except:
            pass
        try:
            shutil.copytree(pjoin(os.path.dirname(self.result_folder_path),existing_run),self.result_folder_path,symlinks=True)
        except:
            out.print_error("Could not copy the result of the existing run.")
        try:
            shutil.copytree(pjoin(os.path.dirname(self.input_folder_path),existing_run),self.input_folder_path,symlinks=True)
        except:
            out.print_error("Could not copy the input of the existing run.")
#}}}
#{{{ def: copy_all_files_in_folder_with_ending(self,src_folder,dest_folder,ending)
    def copy_all_files_in_folder_with_ending(self,src_dir,dst_dir,ending):
        # this deletes the run folder including its logs, inputs and results
        for file_with_ending in glob.iglob(pjoin(src_dir,ending)):
            shutil.copy(file_with_ending, dst_dir)
#}}}
#{{{ def: clean_run_dir(self,folder_path):
    def clean_run_dir(self,folder_path):
        # this cleans up the a run.0, grid.IS0, etc. directory
        # returns True if something was cleaned
        something_cleaned = False
        for files_and_folders in glob.iglob(pjoin(folder_path,"*")):
            if not files_and_folders.rsplit("/",1)[1] in self.contained_in_clean_run_dir:
                try: # remove directories
                    shutil.rmtree(files_and_folders)
                    something_cleaned = True
                except:
                    pass
                try: # remove files
                    os.remove(files_and_folders)
                    something_cleaned = True
                except:
                    pass
        # clean logs in log dir, and only keep .in files
        for files in glob.iglob(pjoin(folder_path,"log","*")):
            if not files.endswith(".in"):
                try: # remove files
                    os.remove(files)
                    something_cleaned = True
                except:
                    pass
        return something_cleaned
#}}}
#{{{ def: add_dir_identifier(self,path,identifier):
    def add_dir_identifier(self,path,identifier):
        # this routine adds an identifier to the folder under path by creating a file with the name identifier and content identifier
        # it also removes all older identifier
        old_ident = glob.iglob(pjoin(path,"*.ident"))
        for ident in old_ident:
            if os.path.isfile(ident):
                os.remove(ident)
        new_ident = pjoin(path,identifier+".ident")
        with open(new_ident, 'w') as new_file:
            new_file.write(identifier+"\n")
            new_file.write("do not remove, this identifies the folder, crucial for result combination and more\n")
#}}}
#{{{ def: get_identifier(self,path):
    def get_identifier(self,path):
        # this routine returns the identifier of a folder under path
        if not os.path.isdir(path):
            out.print_error("Asking for identifier of path %s in routine get_identifier which is no folder." % path)
        try:
            files = glob.glob(pjoin(path,"*.ident"))
            if len(files) > 1:
                out.print_error("Folder %s in get_identifier has has more than one identifier." % path)
            identifier = files[0].rsplit('/',1)[1].replace(".ident","")
        except:
#            out.print_error("Folder %s in routine get_identifier has no identifier." % path)
            identifier = "NOIDENTIFIER"
            pass
        return identifier
#}}}
#{{{ def: get_dirs_with_identifier(self,path,identifier):
    def get_dirs_with_identifier(self,path,identifier):
        # this routine returns a list of dirs (full path) which have the identifier "main"
        if not os.path.isdir(path):
            out.print_error("Asking for main dirs under path %s in routine get_main_dirs which is no folder." % path)
        all_main_dirs = [d for d in glob.iglob(pjoin(path,"run.*")) if os.path.isdir(d) and self.get_identifier(d) == identifier]
        return all_main_dirs
#}}}
#{{{ def: create_dir(self,directory,channel):
    def create_dir(self,directory,channel):
        # this routine creates a run/grid folder including a log folder with all channel.in files
        try:
            os.makedirs(directory)
        except:
            pass
        try:
            os.makedirs(pjoin(directory,"log"))
        except:
            pass
        for chan in channel:
            time_before = time.time()
            while abs(time_before-time.time()) < 3600:
                try:
                    with open(pjoin(directory,"log",chan+".in"),"w+") as f:
                        f.write(chan)
                    break
                except:
                    time.sleep(1)
#}}}
#}}}
