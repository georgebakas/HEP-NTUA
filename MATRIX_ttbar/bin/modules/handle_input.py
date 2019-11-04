#{{{ imports
import os
import readline
import cmd
from os.path import join as pjoin
# own modules
from initialize_classes import out
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
#{{{ class: run_folder_class()
class run_folder_class():
#{{{ def: __init__(self,process_dir_in,no_run_folder_in)
    def __init__(self,process_dir_in,no_run_folder_in):
        # this class handles the naming of the run
        self.process_dir = process_dir_in
        self.no_run_folder = no_run_folder_in
#}}}
#{{{ def: get_run_folder_list(self)
    def get_run_folder_list(self):
        folder_list = [x for x in os.listdir(self.process_dir) if x.startswith("run_") and not x.endswith(".tar")]
        for item in self.no_run_folder:
            if item in folder_list: folder_list.remove(item)
        return folder_list
#}}}
#{{{ def: list_run_folders(self)
    def list_run_folders(self):
        l = self.get_run_folder_list()
        cols = 4
        for i in range(cols - len(l) % 4):
            l.append(" ")
        
        split=[l[i:i+len(l)/cols] for i in range(0,len(l),len(l)/cols)]
        for row in zip(*split):
            print "".join(str.ljust(i,30) for i in row)
#}}}
#{{{ def: readin_run_folder(self)
    def readin_run_folder(self): # reads in run folder as user input
        # get list of folders which can be used for the run
        run_folders = self.get_run_folder_list()
        # use this list for completion in readin
        completer = CustomCompleter(run_folders)
        readline.set_completer(completer.complete)
        readline.parse_and_bind('tab: complete')
        
        # readin folder to be run in
        input = raw_input("|============>> ")
#                         "                "  
        return input
#}}}
#{{{ def: check_run_folder(self,folder)
    def check_run_folder(self,folder): # checks user defined run folder
        check_successful = True
        # if folder is an empty string next_run_folder will be created and used
        if folder == "":
            return check_successful
        # the following folders should not be used
        if not folder.startswith("run_"):
            out.print_info("The name of the run-folder has start with \"run_\", but you have chosen: \"%s\". Try again..." % folder)
            return False
        if folder.endswith(".tar"):
            out.print_info("The name of the run-folder cannot end with \".tar\", since it is no archive, but you have chosen: \"%s\". Try again..." % folder)
            return False
        # if folder does not exist it will be created and everything is fine
        if not os.path.exists(folder):
            return check_successful
        # the folder should not be a file which already exists
        if os.path.isfile(folder):
            out.print_info("\"%s\" is a file and not a folder, choose another one" % folder)
            return False
        # the following folders should not be used
        if folder in self.no_run_folder:
            out.print_info("Foldername \"%s\" can not be used, choose another one" % folder)
            return False
        # existing folders should contain these files
        for filename in ["file_parameter.dat","file_model.dat","file_distribution.dat"]:
            if not os.path.isfile(pjoin(self.process_dir,folder,filename)):
                out.print_info("Folder \"%s\" already exists, but does not appear to be a run-folder" % folder)
                return False
        # stronger checks could be applied wether folder is indeed a run-folder
        return check_successful
#}}}
#{{{ def: get_highest_existing_run(self)
    def get_highest_existing_run(self):
        xx = 00
        for item in self.get_run_folder_list():
            if "run_" in item: 
                run_xx=item
                xx = max(xx,[int(s) for s in run_xx.split("_") if s.isdigit()])
        if xx:
            xx = xx[0]
        else:
            xx = 00
        if(xx == 99):
            out.print_error("Folder \"run_99\" found, cannot create more run folders. Move/Rename/Remove run_99 and restart script.")
        return xx
#}}}
#}}}
#{{{ class: edit_input_cmd(cmd.Cmd)
class edit_input_cmd(cmd.Cmd):
    """Command processor to choose inputs and start the run afterwards."""    
#{{{ def: __init__(self,process_name_in)
    def __init__(self,process_name_in):
        # this class handles the interactive modification of the inputs
        cmd.Cmd.__init__(self)
        self.process_name = process_name_in # get process name like that, since prc instance is not accessible, coming from same modul
        self.input_folder = ""
        self.editor = ""
        self.general_commands = ["help","help <command>","list","exit","quit"]
        self.inputs = ["parameter","model","distribution", "dddistribution"]
        self.run_modes = []
        self.command_description = {}
        self.command_description["help"] = "Show help menu."
        self.command_description["help <command>"] = "Show help message for specific <command>."
        self.command_description["exit"] = "Stop the code."
        self.command_description["quit"] = "Stop the code."
        self.command_description["list"] = "List available commands again."
        self.command_description["parameter"] = "Modify \"parameter.dat\" input file in editor."
        self.command_description["model"] = "Modify \"model.dat\" input file in editor."
        self.command_description["distribution"] = "Modify \"distribution.dat\" input file in editor."
        self.command_description["dddistribution"] = "Modify \"dddistribution.dat\" input file in editor."
        self.command_description["run"] = "Start cross section computation in standard mode."
        self.command_description["run_grid"] = "Start only grid setup phase."
        self.command_description["run_pre"] = "Start only extrapolation (grid must be already done)."
        self.command_description["run_pre_and_main"] = "Start after grid setup (grid must be already done)."
        self.command_description["run_main"] = "Start only main run (other runs must be already done)."
        self.command_description["run_results"] = "Start only result combination."
        self.command_description["run_gnuplot"] = "Start only gnuplotting the results."
        self.command_description["setup_run"] = "Setup the run folder, but not start running."
        self.command_description["delete_run"] = "Remove run folder (including input/log/result)."
        self.command_description["tar_run"] = "Create <run_folder>.tar (including input/log/result)."
#}}}
#{{{ def: preloop(self)
    def preloop(self):
        out.print_read("Type one of the following commands: (\"TAB\" for auto-completion)")
        self.list_general_commands()
        self.list_inputs()
        self.list_run_modes()
#}}}
#{{{ def default(self,line):
    def default(self,line):
        out.print_warning("Command \"%s\" not recognized, type a command from below:" % line)
        self.list_general_commands()
        self.list_inputs()
        self.list_run_modes()        
#}}}
#{{{ def: emptyline(self)
    def emptyline(self):
        pass
#}}}
#{{{ def: help_parameter(self)
    def help_parameter(self):
        out.print_info("Type \"parameter\" to open input file \"parameter.dat\" in editor.")
#}}}
#{{{ def: help_model(self)
    def help_model(self):
        out.print_info("Type \"model\" to open input file \"model.dat\" in editor.")
#}}}
#{{{ def: help_distribution(self)
    def help_distribution(self):
        out.print_info("Type \"distribution\" to open input file \"distribution.dat\" in editor.")

    def help_dddistribution(self):
        out.print_info("Type \"dddistribution\" to open input file \"dddistribution.dat\" in editor.")
#}}}
#{{{ def: help_run(self)
    def help_run(self):
        out.print_info("Type \"run\" to close this command line and start the run which computes the cross section of the process %s according to the input parameters in the input files specified before." % self.process_name)
#}}}
#{{{ def: help_delete_run(self)
    def help_delete_run(self):
        out.print_info("Type \"delete_run\" to close this command line and remove the specified run folder including all inputs, logs and results.")
#}}}
#{{{ def: help_tar_run(self)
    def help_tar_run(self):
        out.print_info("Type \"tar_run\" to close this command line and create a .tar archive of the run folder including all inputs, logs and results.")
#}}}
#{{{ def: help_setup_run(self)
    def help_setup_run(self):
        out.print_info("Type \"setup_run\" to close this command line and start ONLY the set up of the run folder, no run will be done.")
#}}}
#{{{ def: help_run_pre(self)
    def help_run_pre(self):
        out.print_info("Type \"run_pre\" to close this command line and start ONLY the runtime extrapolation run (pre run) without grid setup (which is assumed to already exist).")
#}}}
#{{{ def: help_run_main(self)
    def help_run_main(self):
        out.print_info("Type \"run_main\" to close this command line and start ONLY the main run (and the result collection) without grid setup and runtime extrapolation (which are assumed to already exist).")
#}}}
#{{{ def: help_run_pre_and_main(self)
    def help_run_pre_and_main(self):
        out.print_info("Type \"run_pre_and_main\" to close this command line and start ONLY the runtime extrapolation and main run (and the result collection) without grid setup (which is assumed to already exist).")
#}}}
#{{{ def: help_run_grid(self)
    def help_run_grid(self):
        out.print_info("Type \"run_grid\" to close this command line and start ONLY the grid run which will only set up the grids, but no main run to get physical results.")
#}}}
#{{{ def: help_run_results(self)
    def help_run_results(self):
        out.print_info("Type \"run_results\" to close this command line and start ONLY results run which will try to combination the results of the main run (which is assumed to be already done).")
#}}}
#{{{ def: help_run_gnuplot(self)
    def help_run_gnuplot(self):
        out.print_info("Type \"run_gnuplot\" to close this command line and start ONLY gnuplotting the results (which are assemed to be already there).")
#}}}
#{{{ def: help_exit(self)
    def help_exit(self):
        out.print_info("Type \"exit\" to stop the code.")
#}}}
#{{{ def: help_quit(self)
    def help_quit(self):
        out.print_info("Type \"quit\" to stop the code.")
#}}}
#{{{ def: help_list(self)
    def help_list(self):
        out.print_info("Type \"list\" to list all options with explanations.")
#}}}
#{{{ def: edit_file(self,filename):
    def edit_file(self,filename):
        if self.editor:
            editor = self.editor
        else:
            editor = os.getenv('EDITOR')
        if not editor:
            if os.path.isfile('/usr/bin/editor'):
                editor = '/usr/bin/editor'
        if editor:
            os.system('%s %s' % (editor, filename))
        else:
            out.print_warning("EDITOR not set; set \"default_editor\" variable in MATRIX_configuration filevia or via \"export EDITOR=XXX\" and restart the code.")
#}}}
#{{{ def: do_parameter(self, line)
    def do_parameter(self, line):
        filename=pjoin(self.input_folder,"parameter.dat")
        self.edit_file(filename)
#}}}
#{{{ def: do_model(self, line)
    def do_model(self, line):
        filename=pjoin(self.input_folder,"model.dat")
        self.edit_file(filename)
#}}}
#{{{ def: do_distribution(self, line)
    def do_distribution(self, line):
        filename=pjoin(self.input_folder,"distribution.dat")
        self.edit_file(filename)

    def do_dddistribution(self, line):
        filename=pjoin(self.input_folder,"dddistribution.dat")
        self.edit_file(filename)
#}}}
#{{{ def: do_run(self, line)
    def do_run(self, line):
        self.run_mode = "run"
        return True
#}}}
#{{{ def: do_delete_run(self, line)
    def do_delete_run(self, line):
        self.run_mode = "delete_run"
        return True
#}}}
#{{{ def: do_tar_run(self, line)
    def do_tar_run(self, line):
        self.run_mode = "tar_run"
        return True
#}}}
#{{{ def: do_setup_run(self, line)
    def do_setup_run(self, line):
        self.run_mode = "setup_run"
        return True
#}}}
#{{{ def: do_run_pre(self, line)
    def do_run_pre(self, line):
        self.run_mode = "run_pre"
        return True
#}}}
#{{{ def: do_run_main(self, line)
    def do_run_main(self, line):
        self.run_mode = "run_main"
        return True
#}}}
#{{{ def: do_run_pre_and_main(self, line)
    def do_run_pre_and_main(self, line):
        self.run_mode = "run_pre_and_main"
        return True
#}}}
#{{{ def: do_run_grid(self, line)
    def do_run_grid(self, line):
        self.run_mode = "run_grid"
        return True
#}}}
#{{{ def: do_run_results(self, line)
    def do_run_results(self, line):
        self.run_mode = "run_results"
        return True
#}}}
#{{{ def: do_run_gnuplot(self, line)
    def do_run_gnuplot(self, line):
        self.run_mode = "run_gnuplot"
        return True
#}}}
#{{{ def: do_exit(self, line)
    def do_exit(self, line):
        out.print_info("Exiting...")
        exit(0)
#}}}
#{{{ def: do_quit(self, line)
    def do_quit(self, line):
        out.print_info("Exiting...")
        exit(0)
#}}}
#{{{ def: do_list(self, line)
    def do_list(self,line):
        out.print_read("Type one of the following commands: (\"TAB\" for auto-completion)")
        self.list_general_commands()
        self.list_inputs()
        self.list_run_modes()
#}}}
#{{{ def: list_commands(self)
    def list_general_commands(self):
        # get longest process_id and process by looping through available processes
        length_id = 25
        # for command in self.general_commands:
        #     length_id = max(length_id,len(command)+7)
        print "-"*80
        print "General commands"+" "*(length_id-len("General commands")-5)+"||   description"
        print "-"*80
        for command in self.general_commands:
             print command+" "*(length_id-len(command)-5)+">>   "+self.command_description.get(command,"no description available")
#}}}

#{{{ def: list_commands(self)
    def list_inputs(self):
        # get longest process_id and process by looping through available processes
        length_id = 25
        for command in self.inputs:
            length_id = max(length_id,len(command)+7)
        print "-"*80
        print "Input to modify"+" "*(length_id-len("Input to modify")-5)+"||   description"
        print "-"*80
        for command in self.inputs:
            print command+" "*(length_id-len(command)-5)+">>   "+self.command_description.get(command,"no description available")
#}}}
#{{{ def: list_commands(self)
    def list_run_modes(self):
        # get longest process_id and process by looping through available processes
        length_id = 25
        # for command in self.run_modes:
        #     length_id = max(length_id,len(command)+7)
        print "-"*80
        print "Run-mode to start"+" "*(length_id-len("Run-mode to start")-5)+"||   description"
        print "-"*80
        for command in self.run_modes:
            if not self.command_description.get(command,None) == None:
                print command+" "*(length_id-len(command)-5)+">>   "+self.command_description.get(command,"no description available")
#}}}

#}}}
#{{{ class: process_id()
class process_id(): # class that takes care of everything related to the process name
#{{{ def: __init__(self,script_link_dir,proper_process_names)
    def __init__(self,script_link_dir,proper_process_names):
        # this is how you get the name/id of the process
        try:
            self.process_name = os.path.split(os.path.dirname(script_link_dir))[1].split("_")[0]
        except:
            self.process_name = os.path.split(os.path.dirname(script_link_dir))[1]
            pass
        # check if that is a proper process name
        if not self.process_name in proper_process_names:
            out.print_error("Process id \"%s\" determined by name of current directory \"%s\" (string before first occurence of \"_\" if \"_\" in current directory name) is not in the list of proper_process_names." % (self.process_name,os.path.dirname(script_link_dir).rsplit("/",1)[1]))
#}}}
#{{{ def: get_nice_process_name(self)
    def get_nice_process_name(self):
        # for string in ["ppeeexex04","ppeex02zms","ppeexmxnm04","ppemexmx04","ppenex02","ppexne02","ppexnej","ppmmx22","ppnenex02zms","ppw01","ppyyx22","ppzz02","ppeeexnex04","ppeexa03","ppeexnenex04","ppemexnmx04","ppenex02wms","ppexne02wms","pph21","ppmmx22hms","ppnenexa03","ppwx01","ppyyx22hms","ppaa02","ppeex02","ppeexexne04","ppeexnmnmx04","ppemxnmnex04","ppenexa03","ppexnea03","pphh22","ppnenex02","ppttx20","ppbbx20","ppwxw02","ppz01"]:
        #     length = 12
        #     print string," " * (length-len(string)),"  ===>   ",self.convert_to_nice_process_name(string)
        # exit(0)
        return self.convert_to_nice_process_name(self.process_name)
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
