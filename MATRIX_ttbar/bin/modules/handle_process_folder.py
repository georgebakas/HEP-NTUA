#{{{ class: process_folder()
class process_folder(process_in,process_folder_path_in):
#{{{ def: __init__(self)
    def __init__(self):
        # this class handles the naming of the run
        self.process = process_in
        self.process_folder_path = process_folder_path_in
#}}}
    def create_process_folder():
        if os.path.exists(process_folder_path):
            out.print_warning("Process folder \"%s\" already exists. If you want to create it anew, (re)move this folder and try again. Skipping process folder creation..." % process_folder_path)
        elif os.path.exists(process_default_folder_path):
            out.print_error("The default process folder \"%s\" already exists. This folder should not exist before creating a new MATRIX process folder. You can (re)move this folder and try again. Exiting..." % process_default_folder_path)
        else:
            os.chdir(pjoin(matrix_dir,"run"))
        tar = tarfile.open(process_default_folder_tar)
        tar.extractall()
        tar.close()
        shutil.move(process_default_folder_path,process_folder_path)
        os.chdir(process_folder_path)
        os.makedirs("bin")
        os.makedirs("input")
        os.symlink(matrix_run_process_script,"bin/run_process")
        os.symlink(matrix_configuration_file,"input/MATRIX_configuration")
        # the following can be removed once the tar file have the right content >>>
        try:
            shutil.rmtree(pjoin(process_folder_path,"default.grid.final"))
        except:
            pass
        try:
            shutil.rmtree(pjoin(process_folder_path,"default.MUNICH"))
        except:
            pass
        try:
            shutil.rmtree(pjoin(process_folder_path,"default.MATRIX"))
        except:
            pass
        shutil.rmtree(pjoin(process_folder_path,"batch"))
        shutil.copy(pjoin(matrix_dir,"run",process+"_script","default.grid.final.tar"),pjoin(process_folder_path,"default.grid.final.tar"))
        tar = tarfile.open(pjoin(process_folder_path,"default.grid.final.tar"))
        tar.extractall()
        tar.close()
        shutil.copytree(pjoin(matrix_dir,"run",process+"_script","input","default.input.final"),pjoin(process_folder_path,"input","default.input.final"),symlinks=True)
        os.remove(pjoin("default.grid.final",process))
        os.symlink(process_executable,pjoin("default.grid.final",process))
        os.remove(pjoin("default.grid.final","setup","file_parameter.dat"))
        os.symlink(setup_file_parameter,pjoin("default.grid.final","setup","file_parameter.dat"))
        make_tarfile("default.grid.final.tar","default.grid.final")
        shutil.rmtree("default.grid.final")
        out.print_info("Process folder successfully created.")
#}}}
