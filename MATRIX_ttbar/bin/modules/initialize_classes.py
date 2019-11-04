# import some basic crucial information from main (run_process)
from __main__ import proper_process_names, script_link_dir, no_run_folder, process_dir # define in main function so that directly accessible

# here we initialize instances of classes needed everywhere
# import order is important !!!

# first and most important instance
from handle_output import *
out = print_output(process_dir) # class for handling the on-screen output

from handle_input import *
prc = process_id(script_link_dir,proper_process_names) # class for handling the process id
run_name = run_folder_class(process_dir,no_run_folder) # class for handling the name of the run folder
edit_input = edit_input_cmd(prc.process_name) # class for handling the interactive user inputs

from handle_folder import *
fold = folder_structure() # class for handling the required structure of folder

from handle_log import *
log = run_log() # class for handling the log of the run

from handle_result import *
res  = result() # class for handling the result
cite = citations(prc.process_name) # class for handling citations
