#{{{ imports
import os
import sys
import glob
import textwrap
from os.path import join as pjoin
#}}}

#{{{ class: banner
class banner():
    def __init__(self,insymb,outsymb,size,intend):
        self.insymb = insymb
        self.outsymb = outsymb
        self.size = size
        self.intend = intend
    def initial_print(self):
        print " "*self.intend+"/"+"-"*self.size+"\\"
    def final_print(self):
        print " "*self.intend+"\\"+"-"*self.size+"/"
    def separator_print(self):
        print " "*self.intend+self.insymb+"-"*(self.size+2-len(self.outsymb)-len(self.insymb))+self.outsymb
    def banner_print(self,line):
        print " "*self.intend+self.insymb+" "+line+" "*(self.size-len(line)-len(self.outsymb)-len(self.insymb))+" "+self.outsymb
    def print_center(self,line):
        
        print " "*self.intend+self.insymb+" "*int((self.size-len(line))/2)+line+" "*(self.size-len(line)-int((self.size-len(line))/2))+self.outsymb
#{{{ def: print_matrix(self)
    def print_matrix(self):
        self.initial_print()
#        self.banner_print("MATRIX: A fully-differential NNLO(+NNLL) process library")
        self.banner_print("")
        self.banner_print("          __  __     ___     ____    ___         _     _")
        self.banner_print("         | _\/_ |   / _ \   |____|  / _ \   ||   \\\\   //")
        self.banner_print("         || \/ ||  | |_| |    ||    ||_||   ||    \\\\ //  ")
        self.banner_print("         ||    ||  | ___ |    ||    ||\\\    ||    // \\\\  ")
        self.banner_print("         ||    ||  ||   ||    ||    || \\\   ||   //   \\\\  ")
        self.banner_print("")
#        self.banner_print("         arXiv:1711.06631                               ")
#        self.banner_print("         Version: 1.0.0, Nov 2017       arXiv:1711.06631")
#        self.banner_print("         Version: 1.0.0                         Nov 2017")
        self.banner_print("         Version: Top pair production - Private Release")
        self.banner_print("")
#        self.banner_print("")
#        self.banner_print("         Reference: EPJCXX(2017),X,XXX [arXiv:1711.06631]")
#        self.banner_print("         Reference: JHEP XXXX(2017)XXX [arXiv:1711.06631]")
        self.banner_print("         References: arXiv:1711.06631           Nov 2019")
        self.banner_print("                     arXiv:1901.04005")
        self.banner_print("                     arXiv:1906.06535")
#        self.banner_print("         Reference: arXiv:1711.06631")

        self.banner_print("")
        self.banner_print("Munich -- the MUlti-chaNnel Integrator at swiss (CH) precision --")
        self.banner_print("Automates qT-subtraction and Resummation to Integrate X-sections")
        #self.separator_print()
        self.banner_print("")
        self.banner_print("\          \ ____      \          \ ____      \ ____      \  ")
        self.banner_print(" \          \          |\          \          |\          |\  ")
        self.banner_print("  )====  +   )====  +  | )====  +   )====  +  | )====  +  |-)====")
        self.banner_print(" /          /          |/          /____      |/          |/  ")
        self.banner_print("/          /           /          /           /           / ")
        # self.banner_print("  \             \ ____        \            \ ____       \  ")
        # self.banner_print("   \             \            |\            \           |\  ")
        # self.banner_print("    )====  -|-    )====   -|- | )====  -|-   )==== -|-  |-)====")
        # self.banner_print("   /             /            |/            /____       |/  ")
        # self.banner_print("  /             /             /            /            / ")
        self.banner_print("")
        #self.banner_print("Version:   alpha-0.0.1   Dec 2015")
#        self.banner_print("                        arXiv:1711.06631                        ")
#        self.banner_print("Reference: JHEP XX (2017) XXX [arXiv:1711.06631]")
        self.banner_print("M. Grazzini                              (grazzini@physik.uzh.ch)")
        self.banner_print("S. Kallweit                             (stefan.kallweit@cern.ch)")
        self.banner_print("M. Wiesemann                           (marius.wiesemann@cern.ch)")
#        self.banner_print("Reference: JHEP XX (2017) XXX [arXiv:1711.06631]")
        self.separator_print()
#        self.banner_print("Reference: JHEP XX (2017) XXX [arXiv:1711.06631]")
        self.banner_print("MATRIX is based on a number of different computations and tools")
        self.banner_print("from various people and groups. Please acknowledge their efforts")
        self.banner_print("by citing the list of references which is created with every run.")
        self.final_print()

        # self.initial_print()
        # self.banner_print("MATRIX: A fully-differential NNLO(+NNLL) process library")
        # #self.banner_print("")
        # self.banner_print("          __  __     ___     ____    ___         _     _")
        # self.banner_print("         | _\/_ |   / _ \   |____|  / _ \   ||   \\\\   //")
        # self.banner_print("         || \/ ||  | |_| |    ||    ||_||   ||    \\\\ //  ")
        # self.banner_print("         ||    ||  | ___ |    ||    ||\\\    ||    // \\\\  ")
        # self.banner_print("         ||    ||  ||   ||    ||    || \\\   ||   //   \\\\  ")
        # self.banner_print("")
        # self.banner_print("         Version: alpha-0.0.1                   Dec 2015")
        # self.banner_print("")
        # self.banner_print("Munich -- the MUlti-chaNnel Integrator at swiss (CH) precision --")
        # self.banner_print("Automates qT-subtraction and Resummation to Integrate X-sections")
        # #self.separator_print()
        # self.banner_print("")
        # #self.banner_print("Version:   alpha-0.0.1   Dec 2015")
        # self.banner_print("S. Kallweit                               (kallweit@uni-mainz.de)")
        # self.banner_print("M. Grazzini                              (grazzini@physik.uzh.ch)")
        # self.banner_print("D. Rathlev                                (rathlev@physik.uzh.ch)")
        # self.banner_print("M. Wiesemann                              (mariusw@physik.uzh.ch)")
        # self.banner_print("")
        # self.banner_print("  \             \ ____        \            \ ____       \  ")
        # self.banner_print("   \             \            |\            \           |\  ")
        # self.banner_print("    )====  -|-    )====   -|- | )====  -|-   )==== -|-  |-)====")
        # self.banner_print("   /             /            |/            /____       |/  ")
        # self.banner_print("  /             /             /            /            / ")
        # self.banner_print("")
        # self.separator_print()
        # self.banner_print("MATRIX is based on a number of different computations and tools")
        # self.banner_print("from various people and groups. Please acknowledge their efforts")
        # self.banner_print("by citing the list of references which is created with every run.")
        # self.final_print()
#}}}        
#}}}
#{{{ class: Tee(object)
class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush() # If you want the output to be visible immediately
    def flush(self) :
        for f in self.files:
            f.flush()
#}}}
#{{{ def: crop_file(in_file,out_file)
def crop_file(in_file,out_file):
# remove all empty lines from end of the in_file and write it to out_file
# in_file and out_file can be the same
    """ Remove empty lines at end of file
    """  
    with open(in_file,"r") as infile:
        lines = infile.read()

    while lines.endswith("\n\n"):
        lines = lines[:-1]

    with open(out_file, 'w') as outfile:
        for line in lines:
            outfile.write(line)
#}}}
#{{{ def: tail(f, lines)
def tail(file_name, lines):
    with open(file_name,"r") as f: 
        total_lines_wanted = lines

        BLOCK_SIZE = 1024
        f.seek(0, 2)
        block_end_byte = f.tell()
        lines_to_go = total_lines_wanted
        block_number = -1
        blocks = [] # blocks of size BLOCK_SIZE, in reverse order starting
                    # from the end of the file
        while lines_to_go > 0 and block_end_byte > 0:
            if (block_end_byte - BLOCK_SIZE > 0):
                # read the last block we haven't yet read
                f.seek(block_number*BLOCK_SIZE, 2)
                blocks.append(f.read(BLOCK_SIZE))
            else:
                # file too small, start from begining
                f.seek(0,0)
                # only read what was not read
                blocks.append(f.read(block_end_byte))
            lines_found = blocks[-1].count('\n')
            lines_to_go -= lines_found
            block_end_byte -= BLOCK_SIZE
            block_number -= 1
        all_read_text = ''.join(reversed(blocks))
        return '\n'.join(all_read_text.splitlines()[-total_lines_wanted:])
#}}}
#{{{ class: bcolors
class bcolors:
    if sys.stdout.isatty():
        HEADER = "\033[95m"
        OKBLUE = "\033[94m"
        OKGREEN = "\033[92m"
        WARNING = "\033[93m"
        FAIL = "\033[91m"
        ENDC = "\033[0m"
        BOLD = "\033[1m"
        UNDERLINE = "\033[4m"
    else:
        HEADER = ""
        OKBLUE = ""
        OKGREEN = ""
        WARNING = ""
        FAIL = ""
        ENDC = ""
        BOLD = ""
        UNDERLINE = ""
#}}}
#{{{ class: print_output()
class print_output():
#{{{ def: __init__(self,process_dir_in="")
    def __init__(self,process_dir_in=""):
        global wrapper
        wrapper = textwrap.TextWrapper()
        wrapper.width = 80
        self.process_dir = process_dir_in
        self.run_folder = ""
#}}}
#{{{ def: print_error_no_stop(self,string)
    def print_error_no_stop(self,string):
        wrapper.initial_indent    = "<<MATRIX-ERROR>> "
        wrapper.subsequent_indent = "                 "
        try:
            string = "%s" % "\n".join(wrapper.wrap(string))
            print bcolors.FAIL + string + bcolors.ENDC #if self.on_terminal else string
        except:
            print "%s" % "\n".join(wrapper.wrap(string))
            pass
#}}}
#{{{ def: print_error(self,string)
    def print_error(self,string):
        wrapper.initial_indent    = "<<MATRIX-ERROR>> "
        wrapper.subsequent_indent = "                 "
        try:
            string = "%s" % "\n".join(wrapper.wrap(string))
            print bcolors.FAIL + string + bcolors.ENDC #if self.on_terminal else string
        except:
            print "%s" % "\n".join(wrapper.wrap(string))
            pass
        sys.exit(1)
#}}}
#{{{ def: print_warning(self,string)
    def print_warning(self,string):
        wrapper.initial_indent    = "<<MATRIX-WARN>> "
        wrapper.subsequent_indent = "                "
        try:
            string = "%s" % "\n".join(wrapper.wrap(string))
            print bcolors.WARNING + string + bcolors.ENDC #if self.on_terminal else string
        except:
            print "%s" % "\n".join(wrapper.wrap(string))
            pass
#}}}
#{{{ def: print_info(self,string)
    def print_info(self,string):
        wrapper.initial_indent    = "<<MATRIX-INFO>> "
        wrapper.subsequent_indent = "                "
        try:
            string = "%s" % "\n".join(wrapper.wrap(string))
            print bcolors.OKGREEN + string + bcolors.ENDC #if self.on_terminal else string
        except:
            print "%s" % "\n".join(wrapper.wrap(string))
            pass
#}}}
#{{{ def: print_test(self,string)
    def print_test(self,string):
        wrapper.initial_indent    = "<<MATRIX-TEST>> "
        wrapper.subsequent_indent = "                "
        try:
            string = "%s" % "\n".join(wrapper.wrap(string))
            print bcolors.UNDERLINE + string + bcolors.ENDC #if self.on_terminal else string
        except:
            print "%s" % "\n".join(wrapper.wrap(string))
            pass
#}}}
#{{{ def: print_test_error(self,string)
    def print_test_error(self,string):
        wrapper.initial_indent    = "<<MATRIX-TEST>> "
        wrapper.subsequent_indent = "                "
        try:
            string = "%s" % "\n".join(wrapper.wrap(string))
            print bcolors.UNDERLINE+bcolors.FAIL + string + bcolors.ENDC #if self.on_terminal else string
        except:
            print "%s" % "\n".join(wrapper.wrap(string))
            pass
#}}}
#{{{ def: print_test_warning(self,string)
    def print_test_warning(self,string):
        wrapper.initial_indent    = "<<MATRIX-TEST>> "
        wrapper.subsequent_indent = "                "
        try:
            string = "%s" % "\n".join(wrapper.wrap(string))
            print bcolors.UNDERLINE+bcolors.WARNING + string + bcolors.ENDC #if self.on_terminal else string
        except:
            print "%s" % "\n".join(wrapper.wrap(string))
            pass
#}}}
#{{{ def: print_jobs(self,string)
    def print_jobs(self,string):
        wrapper.width = 200
        wrapper.initial_indent    = "<<MATRIX-JOBS>> "
        wrapper.subsequent_indent = "                "
        print "%s" % "\n".join(wrapper.wrap(string))
        wrapper.width = 80
#}}}
#{{{ def: print_result(self,string)
    def print_result(self,string):
        wrapper.width = 200
        wrapper.initial_indent    = "<MATRIX-RESULT> "
        wrapper.subsequent_indent = "                "
        print "%s" % "\n".join(wrapper.wrap(string))
        wrapper.width = 80
#}}}
#{{{ def: print_read(self,string)
    def print_read(self,string):
        wrapper.initial_indent    = "<<MATRIX-READ>> "
        wrapper.subsequent_indent = "                "
        print "%s" % "\n".join(wrapper.wrap(string))
#}}}
#{{{ def: print_make(self,string)
    def print_make(self,string):
        wrapper.initial_indent    = "<<MATRIX-MAKE>> "
        wrapper.subsequent_indent = "                "
        try:
            string = "%s" % "\n".join(wrapper.wrap(string))
            print bcolors.OKBLUE + string + bcolors.ENDC #if self.on_terminal else string
        except:
            print "%s" % "\n".join(wrapper.wrap(string))
            pass
#}}}
#{{{ def: print_failed_runs(self)
    def print_failed_runs(self):
        failed_folder = pjoin(self.process_dir,"log",self.run_folder,"failed")
        file_list = os.listdir(failed_folder)
        if file_list:
            for entry in file_list:
                print bcolors.FAIL + "|============>> "+entry.replace("-","/").replace("QT/CS","QT-CS").replace(".failed"," failed") + bcolors.ENDC
                print bcolors.FAIL + "|= log-file =>> "+ "last 5 lines:" + bcolors.ENDC
                logfiles = list(reversed(sorted(glob.iglob(pjoin(self.process_dir,self.run_folder,entry.replace("-","/").replace("QT/CS","QT-CS").replace(">>","/log/").replace(".failed","_try*.log"))))))
                logfile = logfiles[0]
                crop_file(logfile,logfile)
                print "\n" + bcolors.FAIL + tail( logfile, 5 ) + bcolors.ENDC + "\n"
        else:
            self.print_info("No failed runs there. What do you want?")
#}}}
#{{{ def: print_list(self,list_name,output_type="default")
    def print_list(self,list_name,output_type="default"):
        # this function does both print a list from the path under list_name, or the list when list_name is a list object
        # define color of output
        color_end = bcolors.ENDC
        if output_type == "error":
            color_start = bcolors.FAIL
        elif output_type == "warning":
            color_start = bcolors.WARNING
        elif output_type == "info":
            color_start = bcolors.OKGREEN
        elif output_type == "default":
            color_start = ""
            color_end = ""
        else:
            self.print_error("Given output_type %s in function print_list not known." % output_type)

            
        if isinstance(list_name, basestring):
            if not os.path.isfile(list_name):
                self.print_warning("Trying to print a list from a file which does not exist. Path: %s" % list_name)
            with open(list_name,'r') as list_file:
                for entry in list_file:
                    print color_start + "|============>> " + entry.strip() + color_end
        else:
            for entry in list_name:
                print color_start + "|============>> " + entry.strip() + color_end
#}}}
#{{{ def: print_list(self,list_name,output_type="default")
    def print_list_no_comments(self,list_name,output_type="default"):
        # this function does both print a list from the path under list_name, or the list when list_name is a list object
        # define color of output
        color_end = bcolors.ENDC
        if output_type == "error":
            color_start = bcolors.FAIL
        elif output_type == "warning":
            color_start = bcolors.WARNING
        elif output_type == "info":
            color_start = bcolors.OKGREEN
        elif output_type == "default":
            color_start = ""
            color_end = ""
        else:
            self.print_error("Given output_type %s in function print_list not known." % output_type)
        
        if isinstance(list_name, basestring):
            if not os.path.isfile(list_name):
                self.print_warning("Trying to print a list from a file which does not exist. Path: %s" % list_name)
            with open(list_name,'r') as list_file:
                for entry in list_file:
                    if "#" in entry:
                        entry = entry.split("#")[0]
                    entry = entry.strip()
                    if entry:
                        print color_start + "|============>> " + entry.strip() + color_end
        else:
            for entry in list_name:
                if "#" in entry:
                    entry = entry.split("#")[0].strip()
                if entry:
                    print color_start + "|============>> " + entry.strip() + color_end
#}}}
#{{{ def: print_last_five_lines_of_file(self,file_name)
    def print_last_five_lines_of_file(self,file_name):
        if not os.path.exists(file_name):
            self.print_error("File %s does not exist. Exiting..." % file_name)
        print bcolors.FAIL + "|============>> "+file_name.replace("-","/").replace("QT/CS","QT-CS").replace(".failed"," failed") + bcolors.ENDC
        if file_name.endswith(".err"):
            print bcolors.FAIL + "|= err-file =>> "+ "last 5 lines:" + bcolors.ENDC
        elif file_name.endswith(".log"):
            print bcolors.FAIL + "|= log-file =>> "+ "last 5 lines:" + bcolors.ENDC
        else:
            print bcolors.FAIL + "|=== file ===>> "+ "last 5 lines:" + bcolors.ENDC
        print "\n" + bcolors.FAIL + tail( file_name, 5 ) + bcolors.ENDC + "\n"
#}}}
#}}}
#{{{ class: output_saver
class output_saver():
    def __init__(self,outfile):
        self.terminal = sys.stdout
        self.log = open(outfile,"w")

    def write(self,message):
        self.terminal.write(message)
        self.log.write(self.remove_colors(message))

    def write_only_onscreen(self,message):
        self.terminal.write(message)

    def write_only_file(self,message):
        self.log.write(self.remove_colors(message))

    def remove_colors(self,message):
        message = message.replace("\033[95m","").replace("\033[94m","").replace("\033[93m","").replace("\033[92m","").replace("\033[91m","").replace("\033[0m","").replace("\033[1m","").replace("\033[4m","")
        return message

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass    
#}}}
