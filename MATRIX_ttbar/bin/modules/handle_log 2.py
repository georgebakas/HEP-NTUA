#{{{ imports
import os
import shutil
import glob
from os.path import join as pjoin
# own modules
from initialize_classes import out, fold
#}}}

#{{{ class: run_log()
class run_log(): # class that takes of the logs during the runs and of the main code
#{{{ def: __init__(self)
    def __init__(self):
#        self.log_folder        = None # initialize and set later in main code
        # the following variables will be returned in the summary of the run
        self.total_runtime     = None
        self.min_time_grid     = None
        self.max_time_grid     = None
        self.min_time_pre      = None
        self.max_time_pre      = None
        self.min_time_main     = None
        self.max_time_main     = None
        self.cumulated_time    = None
        self.events_list_LO    = {}
        self.min_events_LO     = float('inf')
        self.max_events_LO     = 0
        self.total_events_LO   = 0
        self.events_list_NLO   = {}
        self.min_events_NLO    = float('inf')
        self.max_events_NLO    = 0
        self.total_events_NLO  = 0
        self.events_list_NNLO  = {}
        self.min_events_NNLO   = float('inf')
        self.max_events_NNLO   = 0
        self.total_events_NNLO = 0
        self.parallel_list_LO    = {}
        self.min_parallel_LO     = float('inf')
        self.max_parallel_LO     = 0
        self.total_parallel_LO   = 0
        self.parallel_list_NLO   = {}
        self.min_parallel_NLO    = float('inf')
        self.max_parallel_NLO    = 0
        self.total_parallel_NLO  = 0
        self.parallel_list_NNLO  = {}
        self.min_parallel_NNLO   = float('inf')
        self.max_parallel_NNLO   = 0
        self.total_parallel_NNLO = 0
        self.total_parallel      = 0
#}}}
#{{{ def: move_to_folder(self,path)
    def move_to_folder(self,path):
        # this routines moves all files from the log folder into a new folder given by path
        try: # first try to create the new folder
            os.makedirs(path)
        except:
            pass
        # then move all the *.log files into the new folder
        for log_file in glob.iglob(pjoin(fold.log_folder_path,"*.log")):
            if os.path.isfile(log_file):
                shutil.move(log_file,path)
        # also move the folders with logs of the successfull and failed runs into the new folder
        shutil.move(pjoin(fold.log_folder_path,"successful"),path)
        shutil.move(pjoin(fold.log_folder_path,"failed"),path)
        # finallye create new (empty) successful/failed folders for the next run phase
        os.makedirs(pjoin(fold.log_folder_path,"successful"))
        os.makedirs(pjoin(fold.log_folder_path,"failed"))
#}}}
#{{{ def: save_previous(self)
    def save_previous(self,run_mode): # must be dependend on run_mode; only move files and folders which will be overwritten in this runmode
        # this routine checks wether ther are logs in the log folder of this run and moves them to saved_log_$i
        if not fold.log_folder_path:
            out.print_error("Log folder is not yet set, but already saving previous results.")
        log_files = glob.glob(pjoin(fold.log_folder_path,"*.log"))
        log_files += glob.glob(pjoin(fold.log_folder_path,"*.success"))
        if run_mode in ["run","run_grid"]:
            log_files += glob.glob(pjoin(fold.log_folder_path,"grid_run"))
        if run_mode in ["run","run_pre","run_pre_and_main"]:
            log_files += glob.glob(pjoin(fold.log_folder_path,"pre_run"))
        if run_mode in ["run","run_main","run_pre_and_main"]:
            log_files += glob.glob(pjoin(fold.log_folder_path,"main_run"))
        log_files += glob.glob(pjoin(fold.log_folder_path,"failed","*"))
        log_files += glob.glob(pjoin(fold.log_folder_path,"successful","*"))

        save_folder = "saved_log_1" # start with result 1
        save_folders = glob.glob(pjoin(fold.log_folder_path,"saved_log_*"))
        max_number = 0
        for folder in save_folders:
            try:
                if "_" in folder.rsplit('_', 1)[1]:
                    next
                current_number = int(folder.rsplit('_', 1)[1])
            except:
                current_number = 0
                pass
            if current_number > max_number:
                max_number = current_number
        out_folder = pjoin(fold.log_folder_path,"saved_log_%s" % (max_number+1))
        
        if log_files:
            out.print_info("Saving previous log...")
            if os.path.exists(out_folder):
                out.print_error("Folder %s that should be used for saving the previous log already exists." % out_folder)
            else:
                os.makedirs(out_folder)
            for files in log_files:
                shutil.move(files,out_folder)
#}}}    
#{{{ def: code_running_add(self)
    def code_running_add(self):
        # this signals that the main code is running so that no second run can be performed in the same folder
        if not fold.log_folder_path:
            out.print_error("Log folder is not yet set, but already adding main.running.")
        code_running = pjoin(fold.log_folder_path,"main.running")
        pid = os.getpid()
        with open(code_running, 'a') as lock_file:
            lock_file.write(str(pid))
#}}}
#{{{ def: code_running_remove(self)
    def code_running_remove(self):
        # this removes the signal that the main code is running (should be done when aborted or at the end of the code)
        if not fold.log_folder_path:
            exit() # don't think we need an error message in that case
#            out.print_error("Log folder is not yet set, but already trying to remove main.running.")
        code_running = pjoin(fold.log_folder_path,"main.running")
        pid = os.getpid()
        with open(code_running, 'r') as lock_file:
            pid_lock_file = lock_file.readline()
        if int(pid_lock_file) == pid:
            try:
                os.remove(code_running)
            except:
                pass
#}}}
#{{{ def: code_running(self)
    def code_running(self):
        # this checks wether another instance of the main code is already running in that folder
        if not fold.log_folder_path:
            out.print_error("Log folder is not yet set, but already testing wether code is running.")
        code_running = pjoin(fold.log_folder_path,"main.running")
        if os.path.exists(code_running):
            return True
        else:
            return False
#}}}
#{{{ def: clear_list(self,list_name)
    def clear_list(self,list_name):
    # this routine adds a line to the file list_name
        list_path = pjoin(fold.log_folder_path,list_name)
        try:
            os.remove(list_path)
        except:
            pass
#}}}
#{{{ def: add_to_list(self,list_name,line)
    def add_to_list(self,list_name,line):
    # this routine adds a line to the file list_name
        list_path = pjoin(fold.log_folder_path,list_name)
        if os.path.isfile(list_path):
            list_file = open(list_path, 'a') 
        else:
            list_file = open(list_path, 'w') 
        list_file.write(line+"\n")
        list_file.close()
#}}}
#{{{ def: list_exists(self,list_name)
    def list_exists(self,list_name):
    # this routine checks if the list list_name is in log folder
        exists = False
        list_path = pjoin(fold.log_folder_path,list_name)
        if os.path.isfile(list_path):
            exists = True
        return exists
#}}}
#{{{ def: summary_add_events(self,events,path)
    def summary_add_events(self,events,contribution):
    # this routine saves the number of events for each contribution of the main run for the summary (contribution is a path)
        # add to list of events
        if "/LO" in contribution:
            self.events_list_LO[contribution] = events
            self.min_events_LO = min(self.min_events_LO,events)
            self.max_events_LO = max(self.max_events_LO,events)
            self.total_events_LO += events
        if "/NLO" in contribution:
            self.events_list_NLO[contribution] = events
            self.min_events_NLO = min(self.min_events_NLO,events)
            self.max_events_NLO = max(self.max_events_NLO,events)
            self.total_events_NLO += events
        if "/NNLO" in contribution:
            self.events_list_NNLO[contribution] = events
            self.min_events_NNLO = min(self.min_events_NNLO,events)
            self.max_events_NNLO = max(self.max_events_NNLO,events)
            self.total_events_NNLO += events
#}}}
#{{{ def: summary_add_parallel(self,parallel,path)
    def summary_add_parallel(self,parallel,contribution):
    # this routine saves the parallelization for each contribution of the main run for the summary (contribution is a path)
        # add to list of parallel
        if "/LO" in contribution:
            self.parallel_list_LO[contribution] = parallel
            self.min_parallel_LO = min(self.min_parallel_LO,parallel)
            self.max_parallel_LO = max(self.max_parallel_LO,parallel)
            self.total_parallel_LO += parallel
        if "/NLO" in contribution:
            self.parallel_list_NLO[contribution] = parallel
            self.min_parallel_NLO = min(self.min_parallel_NLO,parallel)
            self.max_parallel_NLO = max(self.max_parallel_NLO,parallel)
            self.total_parallel_NLO += parallel
        if "/NNLO" in contribution:
            self.parallel_list_NNLO[contribution] = parallel
            self.min_parallel_NNLO = min(self.min_parallel_NNLO,parallel)
            self.max_parallel_NNLO = max(self.max_parallel_NNLO,parallel)
            self.total_parallel_NNLO += parallel
        self.total_parallel += parallel
#}}}
#{{{ def: summary_clear_parallel(self)
    def summary_clear_parallel(self):
    # this routine clears the previously saved parallelization
        self.parallel_list_LO    = {}
        self.min_parallel_LO     = float('inf')
        self.max_parallel_LO     = 0
        self.total_parallel_LO   = 0
        self.parallel_list_NLO   = {}
        self.min_parallel_NLO    = float('inf')
        self.max_parallel_NLO    = 0
        self.total_parallel_NLO  = 0
        self.parallel_list_NNLO  = {}
        self.min_parallel_NNLO   = float('inf')
        self.max_parallel_NNLO   = 0
        self.total_parallel_NNLO = 0
        self.total_parallel      = 0
#}}}
#}}}
#{{{ class: cluster_log()
class cluster_log(): # class that takes care of cluster logs of the jobs
#{{{ def: __init__(self)
    def __init__(self):
        self.folder = "bla"
#}}}
#}}}
