import os
import shutil
from initialize_classes_matrix import out


#{{{ class: inputs
class inputs():
    """Class to readin user inputs, wrap them and adjust MUNICH inputs"""
#{{{ def: __init__(self)
# in case any default inputs are required
#    def __init__(self):
#}}}
#{{{ def: input_change_entry(self,file_path,parameter,value)
    def input_change_entry(self,file_path,parameter,value):
# function to change a single parameter (all occurences) of MUNICH input files
# if paramater does not exist in this file it is not set
# prop 2do: depending on parameter one could allow only for certain value types
        with open(file_path+".replace",'w') as new_file:
            with open(file_path, 'r') as in_file:
                for in_line in in_file:
                    line=in_line.strip() # strip removes all spaces (including tabs and newlines)
                    # if any line starts with %, # or is an emtpy line (disregarding spaces) it is a comment line and should be skipped
                    if line=="" or line[0]=="%" or line[0]=="#" or not "=" in line: 
                        new_file.write(in_line)
                        continue
                    # split line by "=" and remove spaces again from string
                    # [0] gives you the parameter (before "="-sign)
                    try:
                        file_parameter=line.split('=')[0].strip()
                        # [1] gives you the value (after "="-sign); remove trailing % or # and everything after in this case
                        file_value=line.split('=')[1].split("%")[0].split("#")[0].strip()
                    except:
                        new_file.write(in_line)
                        continue
                    if file_parameter == parameter:
                        new_file.write("%s = %s\n" %(file_parameter,value))
                    else:
                        new_file.write(in_line)
        os.remove(file_path)
        shutil.move(file_path+".replace",file_path)
#}}}
#{{{ def: input_set_entry(self,file_path,parameter,value)
    def input_set_entry(self,file_path,parameter,value):
# function to set a single parameter (all occurences) of MUNICH input files
# if paramater does not exist in this file it added at the end of 
# prop 2do: depending on parameter one could allow only for certain value types
        if not os.path.isfile(file_path):
            open(file_path, 'a').close()
        parameter_set = False
        with open(file_path+".replace",'w') as new_file:
            with open(file_path, 'r') as in_file:
                for in_line in in_file:
                    line=in_line.strip() # strip removes all spaces (including tabs and newlines)
                    # if any line starts with %, # or is an emtpy line (disregarding spaces) it is a comment line and should be skipped
                    if line=="" or line[0]=="%" or line[0]=="#": 
                        new_file.write(in_line)
                        continue
                    # split line by "=" and remove spaces again from string
                    # [0] gives you the parameter (before "="-sign)
                    try:
                        file_parameter=line.split('=')[0].strip()
                        # [1] gives you the value (after "="-sign); remove trailing % or # and everything after in this case
                        file_value=line.split('=')[1].split("%")[0].split("#")[0].strip()
                    except:
                        new_file.write(in_line)
                        continue
                    if file_parameter == parameter:
                        new_file.write("%s = %s \n" %(file_parameter,value))
                        parameter_set = True
                    else:
                        new_file.write(in_line)
                if not parameter_set:
                    new_file.write("%s = %s \n" %(parameter,value))
        os.remove(file_path)
        shutil.move(file_path+".replace",file_path)
#}}}
#{{{ def: input_add_entry(self,file_path,parameter,value)
    def input_add_entry(self,file_path,parameter,value):
# function to add a single parameter at the end of a MUNICH input file
        with open(file_path,'w') as in_file:
            in_file.write("%s = %s \n" %(parameter,value))
#}}}
#{{{ def: input_remove_entry(self,file_path,parameter)
    def input_remove_entry(self,file_path,parameter):
# function to remove a single parameter (all occurences) from MUNICH input files
# if paramater does not exists nothing happens
        with open(file_path+".replace",'w') as new_file:
            with open(file_path, 'r') as in_file:
                for in_line in in_file:
                    line=in_line.strip() # strip removes all spaces (including tabs and newlines)
                    # if any line starts with %, # or is an emtpy line (disregarding spaces) it is a comment line and should be skipped
                    if line=="" or line[0]=="%" or line[0]=="#": 
                        new_file.write(in_line)
                        continue
                    # split line by "=" and remove spaces again from string
                    # [0] gives you the parameter (before "="-sign)
                    file_parameter=line.split('=')[0].strip()
                    # [1] gives you the value (after "="-sign); remove trailing % or # and everything after in this case
                    file_value=line.split('=')[1].split("%")[0].split("#")[0].strip()
                    if file_parameter == parameter:
                        continue
                    else:
                        new_file.write(in_line)
        os.remove(file_path)
        shutil.move(file_path+".replace",file_path)
#}}}
#{{{ def: input_search_keyword_set_entry_below(self,file_path,parameter,value)
    def input_search_keyword_set_entry_below(self,file_path,parameter,value):
# function to set single parameter (first occurence) of MUNICH input files
# prop 2do: depending on parameter one could allow only for certain value types
        # initial condition: keyword_needs to be found, and parameter needs to be set
        keyword_found = False
        parameter_set = False
        with open(file_path+".replace",'w') as new_file:
            with open(file_path, 'r') as in_file:
                for in_line in in_file:
                    line=in_line.strip() # strip removes all spaces (including tabs and newlines)
                    # if any line starts with %, # or is an emtpy line (disregarding spaces) it is a comment line and should be skipped
                    # or once the parameter was set, stop searching for more occurences and simply write out the rest of the file
                    if line=="" or line[0]=="%" or line[0]=="#" or parameter_set: 
                        new_file.write(in_line)
                        continue
                    # split line by "=" and remove spaces again from string
                    # [0] gives you the parameter (before "="-sign)
                    file_parameter=line.split('=')[0].strip()
                    # [1] gives you the value (after "="-sign); remove trailing % or # and everything after in this case
                    file_value=line.split('=')[1].split("%")[0].split("#")[0].strip()
                    # check parameter names until you find the one of the key with the right value
                    if file_parameter == ord_params_keyword[parameter][0] and file_value == ord_params_keyword[parameter][1]:
                        keyword_found = True
                    # when keyword found search below the parameter you want to set
                    if keyword_found == True and file_parameter == ord_params_keyword[parameter][2]:
                        new_file.write(file_parameter+" = "+value+" \n")
                        # once the parameter was set, stop searching for more occurences
                        parameter_set = True
                    else:
                        new_file.write(in_line)
        os.remove(file_path)
        shutil.move(file_path+".replace",file_path)
        if not keyword_found:
            out.print_warning("Keyword \"%s\" with value \"%s\" in file_parameter.dat was not found, could not set parameter \"%s\" below to value \"%s\", continuing..."%(ord_params_keyword[parameter][0],ord_params_keyword[parameter][1],parameter,value))
        if not parameter_set:
            out.print_warning("Keyword \"%s\" with value \"%s\" in file_parameter.dat was found, but could not find parameter \"%s\" below and set it to value \"%s\", continuing..."%(ord_params_keyword[parameter][0],ord_params_keyword[parameter][1],ord_params_keyword[parameter][2],value))

#}}}
#{{{ def: input_read_parameter_dat(self,file_path,parameter_list)
    def input_read_parameter_dat(self,file_path,parameter_list):
# function to read all parameters from the MATRIX input file (file_path, "parameter.dat")
# and write it into a dictionary (parameter_list)
# works also for MATRIX configuration file
        if not os.path.isfile(file_path):
            out.print_error("file "+file_path+" in function input_read_parameter_dat does not exist!")
        with open(file_path, 'r') as param_file:
            for in_line in param_file:
                line=in_line.strip() # strip removes all spaces (including tabs and newlines)
                # if any line starts with %, # or is an emtpy line (disregarding spaces) it is a comment line and should be skipped
                if line=="" or line[0]=="%" or line[0]=="#": 
                    continue
                # split line by "=" and remove spaces again from string
                # [0] gives you the parameter (before "="-sign)
                file_parameter=line.split('=')[0].strip()
                # [1] gives you the value (after "="-sign); remove trailing % or # and everything after in this case
                file_value=line.split('=')[1].split("%")[0].split("#")[0].strip()
                parameter_list[file_parameter]=file_value
#}}}
#{{{ def: input_set_file_parameter_from_list(self,file_path,parameter_list)
    def input_set_file_parameter_from_list(self,file_path,parameter_list):
# function to set all parameters in MUNICH input ("file_parameter.dat") from parameter_list 
# read from the MATRIX input file
        for parameter in parameter_list:
            value = parameter_list[parameter] 
            if parameter in unique_parameters:
                # all parameters uniquely defined can be directly set as MUNICH input
                self.input_change_entry(file_path,parameter,value)
            elif parameter in ordered_parameters:
                # these parameters require a certain order, the following function sets the dependend on some keyword
                # the keywords for specific parameters are defined as a dictionary in the very beginning
                self.input_search_keyword_set_entry_below(file_path,parameter,value)
            elif parameter in renamed_parameters:
                parameter = renamed_parameter_mappings[parameter]
                self.input_change_entry(file_path,parameter,value)
            elif parameter in MATRIX_parameters:
                # nothing to do in that case
                continue
            elif parameter in special_parameters:
                # these are special cases for the paramter inputs
                if parameter == "switch_distribution": 
                    self.input_change_entry(file_path,parameter,value)
                    parameter = "switch_output_distribution"
                    self.input_change_entry(file_path,parameter,value)
                elif parameter == "dynamic_scale":
                    self.input_change_entry(file_path,parameter,value)
                    parameter = "dynamic_scale_CV" 
                    self.input_change_entry(file_path,parameter,value)
                elif parameter == "scale_variation":
                    if value == "1": # 7-point variation
                        self.input_change_entry(file_path,"switch_CV","5")
                        self.input_change_entry(file_path,"n_scales_CV","7")
                    if value == "2": # 9-point variation
                        self.input_change_entry(file_path,"switch_CV","6")
                        self.input_change_entry(file_path,"n_scales_CV","9")
                elif parameter == "frixione_fixed_ET_max":
                    parameter = "frixione_epsilon" 
                    value_new = float(value)/10
                    self.input_change_entry(file_path,parameter,value_new)
                elif parameter == "flavor_scheme":
                    if value == "1": # 4FS
                        self.input_change_entry(file_path,"N_f","4")
                        self.input_change_entry(file_path,"N_f_active","4")
                    elif value == "2": # 5FS
                        self.input_change_entry(file_path,"N_f","5")
                        self.input_change_entry(file_path,"N_f_active","5")
                continue
            elif parameter.startswith("user_"): # use this to include all possibly newly user-defined parameters
                # user defined parameters are process dependent and case dependent
                if prc.process_name in ["ppeeexex04"] and parameter in ["user_switch M_leplep","user_cut min_M_leplep","user_cut max_M_leplep","user_cut min_M_leplep_IR"]: # very special case to choose either ATLAS or CMS pairing
                    if not parameter == "user_switch M_leplep":
                        self.input_set_entry(file_path,parameter,value)                        
                    elif value == "0":
                        self.input_set_entry(file_path,parameter,value)                        
                    elif value == "1": # this case set ATLAS paring
                        self.input_set_entry(file_path,parameter,value)                        
                        self.input_set_entry(file_path,"user_switch M_leplep_ATLAS","1")
                        self.input_set_entry(file_path,"user_switch M_leplep_CMS","0")
                    elif value == "2": # this case set CMS paring
                        value = "1"
                        self.input_set_entry(file_path,parameter,value)                        
                        self.input_set_entry(file_path,"user_switch M_leplep_ATLAS","0")
                        self.input_set_entry(file_path,"user_switch M_leplep_CMS","1")
                else: # standard behavior
                    self.input_set_entry(file_path,parameter,value)
                    # not needed anymore since switches are now explicit in input file
#                     # and its switch
#                     if parameter in user_parameter_switches: # pre-defined hard-coded switch mappings
#                         switch = user_parameter_switches[parameter]
#                     else: # default behaviour
#                         if parameter.split()[1].startswith("max_"):
#                             switch = "user_switch "+parameter.split()[1][parameter.split()[1].startswith("max_") and len("max_"):]
#                         elif parameter.split()[1].startswith("min_"):
#                             switch = "user_switch "+parameter.split()[1][parameter.split()[1].startswith("min_") and len("min_"):]
#                         else:
#                             switch = "user_switch "+parameter.split()[1]

#                     switch_value = "1"
# #                    print "switch   ",switch,switch_value
#                     self.input_set_entry(file_path,switch,switch_value)
                continue
            else:
                out.print_error("Parameter \"%s\" is not listed as proper input parameter" % parameter)
#}}}
#{{{ def: input_user_cuts_default_and_consistency(self,file_path,parameter_list)
# 2do: do we need this again, to check user inputs, or can we use the normal check_parameter_consistencies routine ?

#     def input_user_cuts_default_and_consistency(self,file_path,parameter_list): # incomplete (not crucial)
#         # this function sets the default parameters (all switches off) for the user cuts
#         # and checks wether there are inconsistencies user inputs
#         out.print_info("Checking parameter input in parameter.dat file...")
#         # first check if all mandatory parameter are set
#         # for parameter in mandatory_parameters:
#         #     if not parameter in parameter_list:
#         #         out.print_error("Parameter \"%s\" is mandatory, but not defined in parameter.dat file." % parameter)
#         # turn of switches for parameters not defined in parameter.dat
#         for parameter in user_parameters[prc.process_name]:
#             if not parameter in parameter_list:
#                 switch_need_for_other_parameter = False
#                 if parameter in user_parameter_switches: # pre-defined hard-coded mappings
#                     # for other_parameter in user_parameters[prc.process_name]:
#                     #     if other_parameter in parameter_list and user_parameter_switches[other_parameter] == user_parameter_switches[parameter]:
#                     #         switch_need_for_other_parameter = True
#                     switch = user_parameter_switches[parameter]
#                 else: # default behaviour
#                     if parameter.split()[1].startswith("max_"):
#                         switch = "user_switch "+parameter.split()[1][parameter.split()[1].startswith("max_") and len("max_"):]
#                     elif parameter.split()[1].startswith("min_"):
#                         switch = "user_switch "+parameter.split()[1][parameter.split()[1].startswith("min_") and len("min_"):]
#                     else:
#                         switch = "user_switch "+parameter.split()[1]
#                     # for other_parameter in user_parameters[prc.process_name]:
#                     #     if other_parameter in parameter_list and not other_parameter in user_parameter_switches:
#                     #         if parameter.split()[1].startswith("max_"):
#                     #             other_switch = "user_switch "+parameter.split()[1][parameter.split()[1].startswith("max_") and len("max_"):]
#                     #         elif parameter.split()[1].startswith("min_"):
#                     #             other_switch = "user_switch "+parameter.split()[1][parameter.split()[1].startswith("min_") and len("min_"):]
#                     #         else:
#                     #             other_switch = "user_switch "+parameter.split()[1]
#                     #         if other_switch == switch:
#                     #             switch_need_for_other_parameter = True
#                 switch_value = "0"
# #                print "switch   ",switch,switch_value,switch_need_for_other_parameter
# #                if not switch_need_for_other_parameter: # not needed because inputs are set and overwritten afterwards anyway
#                 self.input_set_entry(file_path,switch,switch_value)
#}}}
#{{{ def: input_check_parameter_consistencies_from_list(self,file_path,parameter_list)
    def input_check_parameter_consistencies_from_list(self,file_path,parameter_list): # incomplete (not crucial)
        # this function checks wether there are inconsistencies in the input of the parameter.dat file
        # this will be a long file with if clauses; it will also set defaults to parameters that are not mandatory
        # first check if all mandatory parameter are set
        for parameter in mandatory_parameters:
            if not parameter in parameter_list:
                out.print_error("Parameter \"%s\" is mandatory, but not set in parameter.dat file." % parameter)
        # handle all pre-defined user defined parameters as mandatory; otherwise unexpected behavior could occur when commenting 
        # a parameter, which is actually set in the default file_parameter.dat
        for parameter in user_parameters[prc.process_name]: # these are the !pre-defined! user parameters
            if not parameter in parameter_list:
                out.print_error("Parameter \"%s\" is a mandatory user parameter, but not set in parameter.dat file." % parameter)
        # then set default values for all parameters that are not given
        for parameter in default_parameters:
            if not parameter in parameter_list:
                parameter_list[parameter] = default_parameters[parameter]
                out.print_warning("Parameter \"%s\" not set in parameter.dat file. Setting it to default value: %s." % (parameter,default_parameters[parameter]))
        # then check wether the given values meet the requirements
        # 2do

        # catch specific cases
        if parameter_list.get("frixione_isolation") in ["1","2"]: # frixione isolation turned on 
            frixione_list = ["frixione_n","frixione_delta_0"] # required parameters for frixione isolation
            for parameter in frixione_list:
                if not parameter in parameter_list: # catch if these parameters are not given
                    out.print_error("\"frixione_isolation\" requires input of \"%s\". Please specify in parameter.dat file and restart." % parameter)            
        if parameter_list.get("frixione_isolation") == "1": # frixione isolation ATLAS setup
            if not "frixione_epsilon" in parameter_list:
                out.print_error("\"frixione_isolation\" set to \"1\" (ATLAS setup) requires input of \"frixione_epsilon\". Please specify in parameter.dat file and restart.")
            if "frixione_fixed_ET_max" in parameter_list:
                out.print_error("\"frixione_isolation\" set to \"1\" (ATLAS setup) does not allow for input of \"frixione_fixed_ET_max\". Please remove (comment) from parameter.dat file and restart.")
        elif parameter_list.get("frixione_isolation") == "2": # frixione isolation CMS setup
            if not "frixione_fixed_ET_max" in parameter_list:
                out.print_error("\"frixione_isolation\" set to \"2\" (CMS setup) requires input of \"frixione_fixed_ET_max\". Please specify in parameter.dat file and restart.")
            if "frixione_epsilon" in parameter_list:
                out.print_error("\"frixione_isolation\" set to \"2\" (CMS setup) does not allow for input of \"frixione_epsilon\". Please remove (comment) from parameter.dat file and restart.")
        if float(parameter_list.get("max_time_per_job")) < 1: # max_time_per_job become unreliable below 1h
            out.print_warning("Parameter max_time_per_job (set to %s hours) chosen below 1 hour. This is fine to tune the degree of parallelization, but does not constitute a realistic maximal run time of the jobs." % parameter_list["max_time_per_job"])
        if not parameter_list.get("flavor_scheme","not-set") in ["1","2","not-set"]: # max_time_per_job become unreliable below 1h
            out.print_error("Parameter flavor_scheme in parameter.dat (set to %s) chosen different from \"1\" (four-flavor scheme) and \"2\" (five-flavor scheme). Please choose either of these two values and restart..." % parameter_list["flavor_scheme"])
#}}}
#{{{ def: input_read_SLHA(self,file_path,SLHA_list)
    def input_read_SLHA(self,file_path,SLHA_list):
# function to read a file in the SLHA format, returning a dictionary with:
# BLOCK[number]: [value, commet]
        if not os.path.isfile(file_path):
            out.print_error("file "+file_path+" in function input_read_SLHAdoes not exist!")
        with open(file_path, 'r') as model_file:
            for in_line in model_file:
                line=in_line.strip() # strip removes all spaces (including tabs and newlines)
                # if any line starts with %, # or is an emtpy line (disregarding spaces) it is a comment line and should be skipped
                if line=="" or line[0]=="%" or line[0]=="#": 
                    continue
#                print line

                # split line by spaces and remove otherspaces again from string
                # [0] gives you the number (before the first space)
                file_number=line.split()[0].strip() # 2do multiple spaces
#                print file_number
                file_value=line.split()[1].strip()
#                print file_value
                try:
                    file_comment=line.split('#')[1].strip()
                except: 
                    file_comment=""
                    pass
#                print file_comment

                if file_number.upper() == "DECAY":
                    Block=file_number
                    file_number=line.split()[1].strip()
                    file_value=line.split()[2].strip()
                    SLHA_list[Block][int(file_number)]=file_value
                    SLHA_list[Block+" comment"][int(file_number)]=file_comment
                    continue
 
                if file_number.upper() == "BLOCK":
                    Block=file_value
                else:
                    SLHA_list[Block][int(file_number)]=file_value
                    SLHA_list[Block+" comment"][int(file_number)]=file_comment
#}}}
#{{{ def: input_set_file_model_from_SLHA(self,file_path,SLHA_list)
    def input_set_file_model_from_SLHA(self,file_path,SLHA_list):
# function to set all model paramaters in MUNICH input ("file_model.dat") from SLHA_list 
# read from the MATRIX input file
        for Block in SLHA_list:
            try:
                split=Block.split()[1].strip()
                if not split=="comment":
                    out.print_error("Block in input_set_file_model_from_SLHA contains spaces.")
            except:
                for number in SLHA_list[Block]:              
                    parameter = model_mappings_to_MUNICH[Block][number]
                    value = SLHA_list[Block][number]
                    # add some sepecial case consistency checks
                    if parameter_list.get("flavor_scheme") == "1" and parameter == "M_b" and float(value) == 0.: # 4FS
                        out.print_warning("Bottom mass (block %s entry %s) in model.dat set to M_b=0, but four-flavor scheme (flavor_scheme=1) with massive bottom quarks chosen. Using default value for bottom mass (M_b=4.75 GeV)." %(Block,number))
                        value = 4.75
                    elif parameter_list.get("flavor_scheme") == "2" and parameter == "M_b" and float(value) != 0.: # 5FS
                        out.print_warning("Bottom mass (block %s entry %s) in model.dat set to M_b=%s GeV, but five-flavor scheme (flavor_scheme=1) with massless bottom quarks chosen. Using M_b=0 for bottom mass." %(Block,number,str(float(value))))
                        value = 0.
                    self.input_set_entry(file_path,parameter,value)
#}}}
#{{{ def: add_directories_to_result_files(self,file_path)
    def add_directories_to_result_files(self,file_path,runtime_table,phase):
# function to set add the directories that have to be taken into account
        keyword_found = False
        with open(file_path+".replace",'w') as new_file:
            with open(file_path, 'r') as in_file:
                for in_line in in_file:
                    line=in_line.strip() # strip removes all spaces (including tabs and newlines)
                    # if any line starts with %, # or is an emtpy line (disregarding spaces) it is a comment line and should be skipped
                    # or once the parameter was set, stop searching for more occurences and simply write out the rest of the file
                    if line=="" or line[0]=="%" or line[0]=="#": 
                        new_file.write(in_line)
                        continue
                    # split line by "=" and remove spaces again from string
                    # [0] gives you the parameter (before "="-sign)
                    file_parameter=line.split('=')[0].strip()
                    # [1] gives you the value (after "="-sign); remove trailing % or # and everything after in this case
                    file_value=line.split('=')[1].split("%")[0].split("#")[0].strip()
                    # check parameter names until you find the one of the key with the right value
                    if file_parameter == "directory" and not keyword_found:
                        continue
                    # when keyword found search below the parameter you want to set
                    if (keyword_found and file_parameter == "directory") or (file_parameter == "type_contribution" and keyword_found):
                        # remove from right to "/" (remove "/run.0")
                        # alternative to determine directory directly:
                        directory = pjoin(file_path.rsplit('/',1)[1].rsplit('.',2)[0],file_path.rsplit('/',1)[1].rsplit('.',2)[1],type_contribution+"%s" % (".QCD" if not type_contribution=="born" and not type_contribution=="loop" else ""))
                        # get directory from the value in the file
                        # directory = file_value.rsplit('/', 1)[0]
                        # print directory
                        # exit(0)
                        path = pjoin(fold.run_folder_path,directory)
                        run_dirs_in_path = [ f for f in os.listdir(path) if f.startswith("run.") and os.path.isdir(pjoin(path,f)) ]
                        max_parallel_runs = len(run_dirs_in_path) # result combination now uses always all folders which are there
                        # determine wether to include the extrapolation runs or not
                        # default is to include them
                        start_folder_index = 0
                        # in case the number of events of all extrapolation events is 100 times smaller than the events in a single main run folder exclude them
                        run_dir_0 = pjoin(path,"run.0")
                        pre_events   = run.get_pre_run_min_events_for_contribution(run_dir_0)
                        parameter_files = []
                        main_dirs = fold.get_dirs_with_identifier(path,"main")
                        for main_dir in main_dirs:
                            parameter_files += glob.iglob(pjoin(main_dir,"log","file_parameter*"))
                        event_sum_main = 0
                        for parameter_file in list(parameter_files):
                            parameter = {}
                            self.input_read_parameter_dat(parameter_file,parameter)
                            event_sum_main += float(parameter.get("n_events_min",0))
                        # since we computed the sum of events for this contribution here, save it to the log for the final summary
                        if event_sum_main > 0: # make sure to only include contributions where events were run (otherwise we will end up with logs of contributions that were not run)
                            log.summary_add_events(event_sum_main,path)
                        if event_sum_main/pre_events > 500:
                            start_folder_index = max_parallel_runs - len(main_dirs)
                        if start_folder_index >= max_parallel_runs:
                            start_folder_index = 0 # to avoid that the directory gets removed from the file
                            max_parallel_runs = 1 # to avoid that the directory gets removed from the file
                        # use pre-defined switch if set
                        if "include_pre_in_results" in parameter_list and len(main_dirs)>0:
                            if parameter_list["include_pre_in_results"] == "0" and phase > 0: # take only main runs into account (not in extrapolation phase=-1)
                                start_folder_index = max_parallel_runs - len(main_dirs)
                            elif parameter_list["include_pre_in_results"] == "1" or phase == -1: # take also all pre runs into account
                                start_folder_index = 0
                            else:
                                out.print_error("Parameter \"include_pre_in_results\" in parameter.dat can only have values \"0\" and \"1\", give value: %s." % parameter_list["include_pre_in_results"])
                        for i in range(start_folder_index,max_parallel_runs): # add all folders of the parallel runs
                            new_file.write("directory"+" = "+pjoin(directory,"run.%s"%i)+" \n")
                        # once we know the directory add it as often as needed
                    if file_parameter == "type_contribution":
                        new_file.write(in_line)
                        type_contribution = file_value
                        keyword_found = True
                    elif file_parameter == "directory":
                        keyword_found = False
                    else:
                        new_file.write(in_line)
                    if file_parameter == "type_contribution" and not keyword_found:
                        keyword_found = True
                        type_contribution = file_value
                else:
                    if keyword_found:
                        # remove from right to "/" (remove "/run.0")
                        # alternative to determine directory directly:
                        directory = pjoin(file_path.rsplit('/',1)[1].rsplit('.',2)[0],file_path.rsplit('/',1)[1].rsplit('.',2)[1],type_contribution+"%s" % (".QCD" if not type_contribution=="born" and not type_contribution=="loop" else ""))
                        # get directory from the value in the file
                        # directory = file_value.rsplit('/', 1)[0]
                        # print directory
                        # exit(0)
                        path = pjoin(fold.run_folder_path,directory)
                        run_dirs_in_path = [ f for f in os.listdir(path) if f.startswith("run.") and os.path.isdir(pjoin(path,f)) ]
                        max_parallel_runs = len(run_dirs_in_path) # result combination now uses always all folders which are there
                        # determine wether to include the extrapolation runs or not
                        # default is to include them
                        start_folder_index = 0
                        # in case the number of events of all extrapolation events is 100 times smaller than the events in a single main run folder exclude them
                        run_dir_0 = pjoin(path,"run.0")
                        pre_events   = run.get_pre_run_min_events_for_contribution(run_dir_0)
                        parameter_files = []
                        main_dirs = fold.get_dirs_with_identifier(path,"main")
                        for main_dir in main_dirs:
                            parameter_files += glob.iglob(pjoin(main_dir,"log","file_parameter*"))
                        event_sum_main = 0
                        for parameter_file in list(parameter_files):
                            parameter = {}
                            self.input_read_parameter_dat(parameter_file,parameter)
                            event_sum_main += float(parameter.get("n_events_min",0))
                        # since we computed the sum of events for this contribution here, save it to the log for the final summary
                        if event_sum_main > 0: # make sure to only include contributions where events were run (otherwise we will end up with logs of contributions that were not run)
                            log.summary_add_events(event_sum_main,path)
                        if event_sum_main/pre_events > 500:
                            start_folder_index = int(sorted(main_dirs)[0].rsplit('.',1)[1])
                        if start_folder_index >= max_parallel_runs:
                            start_folder_index = 0 # to avoid that the directory gets removed from the file
                            max_parallel_runs = 1 # to avoid that the directory gets removed from the file
                        # use pre-defined switch if set
                        if "include_pre_in_results" in parameter_list:
                            if parameter_list["include_pre_in_results"] == 0 and phase > 0: # take only main runs into account (not in extrapolation phase=-1)
                                start_folder_index = int(sorted(main_dirs)[0].rsplit('.',1)[1])
                            elif parameter_list["include_pre_in_results"] == 1 or phase == -1: # take also all pre runs into account
                                start_folder_index = 0
                            else:
                                out.print_error("Parameter \"include_pre_in_results\" in parameter.dat can only have values \"0\" and \"1\", give value: %s." % parameter_list["include_pre_in_results"])
                        for i in range(start_folder_index,max_parallel_runs): # add all folders of the parallel runs
                            new_file.write("directory"+" = "+pjoin(directory,"run.%s"%i)+" \n")
                        # once we know the directory add it as often as needed
        os.remove(file_path)
        shutil.move(file_path+".replace",file_path)
#}}}
#{{{ def: get_cross_sections_from_file(self,file_path)
    def get_cross_sections_from_file(self,file_path):
# read total rates (possibly within cuts) from MUNICH output
        if not os.path.isfile(file_path):
            out.print_error("file "+file_path+" in function get_cross_sections_from_file does not exist!")
        variation = []
        with open(file_path, 'r') as param_file:
            for in_line in param_file:
                line=in_line.strip() # strip removes all spaces (including tabs and newlines)
                # if any line starts with %, # or is an emtpy line (disregarding spaces) it is a comment line and should be skipped
                if line=="" or line[0]=="%" or line[0]=="#": 
                    continue
                muRmuF=float(line.split()[0].strip())
                cross_section=float(line.split()[1].strip())
                err=float(line.split()[2].strip())
                if muRmuF==1.:
                    central     = cross_section
                    central_err = err
                variation.append(cross_section)
        up   = max(variation)
        down = min(variation)
        return central, central_err, up, down
#}}}
#{{{ def: input_read_distribution_dat(self,file_path,content)
    def input_read_distribution_dat(self,file_path):
# this function reads in the whole file and returns it
        with open(file_path, "r") as f:
            content = f.read()
        return content
#}}}
#{{{ def: input_set_file_distribution_dat(self,file_path,content)
    def input_set_file_distribution_dat(self,file_path,content):
# this function appends content to the file file_path
        with open(file_path, "a") as f:
            f.write(content)
#}}}
#}}}
