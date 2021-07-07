#{{{ imports
import os
import shutil
import glob
import math
from os.path import join as pjoin 
# own modules
from initialize_classes import out, fold
#}}}
#{{{ class: result()
class result(): # class to deal with the results, adding files, copying them to the results folder etc.
#{{{ def: save_input_with_results(self)
    def save_input_with_result(self):
        # this routine saves the input (at beginning of run) in the result folder
        dest_folder = pjoin(fold.result_folder_path,"input_of_run")
        try:
            os.makedirs(dest_folder)
        except:
            pass
        file_to_copy = pjoin(fold.input_folder_path,"parameter.dat")
        shutil.copy(file_to_copy,dest_folder)
        file_to_copy = pjoin(fold.input_folder_path,"distribution.dat")
        shutil.copy(file_to_copy,dest_folder)
        file_to_copy = pjoin(fold.input_folder_path,"model.dat")
        shutil.copy(file_to_copy,dest_folder)
#}}}    
#{{{ def: save_previous(self)
    def save_previous(self):
        # this routine checks wether there are results in the result folder of this run and moves them to saved_result_$i
        result_sub_folder = glob.glob(pjoin(fold.result_folder_path,"*-run"))
        if os.path.isdir(pjoin(fold.result_folder_path,"summary")):
            result_sub_folder.append(pjoin(fold.result_folder_path,"summary"))
#        if os.path.isdir(pjoin(fold.result_folder_path,"input_of_run")):
#            result_sub_folder.append(pjoin(fold.result_folder_path,"input_of_run"))
        if os.path.isdir(pjoin(fold.result_folder_path,"gnuplot")):
            result_sub_folder.append(pjoin(fold.result_folder_path,"gnuplot"))
        save_folder = "saved_result_1" # start with result 1
        save_folders = glob.glob(pjoin(fold.result_folder_path,"saved_result_*"))
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
        out_folder = pjoin(fold.result_folder_path,"saved_result_%s" % (max_number+1))
        if result_sub_folder:
            out.print_info("Saving previous result...")
            if os.path.exists(out_folder):
                out.print_error("Folder %s that should be used for saving the previous result already exists." % out_folder)
            else:
                os.makedirs(out_folder)
            for folder in result_sub_folder:
                shutil.move(folder,out_folder)
            try:
                shutil.copytree(pjoin(fold.result_folder_path,"input_of_run"),pjoin(out_folder,"input_of_run"))
            except:
                pass
        try:
            shutil.move(pjoin(fold.result_folder_path,"CITATIONS.bib"),pjoin(out_folder,"CITATIONS.bib"))
        except:
            pass
#}}}    
#{{{ def: add_result_files(self,file_list,out_file)
    def add_result_files(self,file_list,out_file):
        # this routine adds a number of files (that MUST have the same structure), keeping the first column fixed,
        # adding the second column linear, and adding the third in quadature, and writes it to out_file
        firsttime = True
        for file_path in file_list:
            new_data = self.readin_file(file_path) # reads space-separated matrix from file
            if firsttime: # initialize data array
                data = [[0 for col in range(len(new_data[0]))] for row in range(len(new_data))]
                firsttime = False
            for row in range(len(new_data)):
                for column in range(len(new_data[row])):
                    if column == 0:
                        data[row][column] = float(new_data[row][column])
                    elif column == 1:
                        data[row][column] += float(new_data[row][column])
                    elif column == 2:
                        data[row][column] += float(new_data[row][column])**2
                    else:
                        out.print_error("The rate file added in add_rate_files contains more than 3 columns")
        with open(out_file,'w') as f:
            for row in range(len(data)):
                line = ""
                for column in range(len(data[row])):
                    length = 15
                    if column == 2:
                        line = line+" "*(length-len(str("%.8g" % math.sqrt(data[row][column]))))+"%.8g" % math.sqrt(data[row][column])
                    else:
                        line = line+" "*(length-len("%.8g" % data[row][column]))+"%.8g" % data[row][column]
                f.write(line+"\n")
            f.close()
#}}}    
#{{{ def: get_central_min_max_from_distribution_files(self,file_list,out_file)
    def get_central_min_max_from_distribution_files(self,file_list,out_file):
        # this routine combines a number of files (that MUST have the same structure), keeping the first column fixed,
        # having the central (first element of file list) in the second and third column, and computing the minumum 
        # and maximum for the columns 4-7
        firsttime = True
        for file_path in file_list:
            new_data = self.readin_file(file_path) # reads space-separated matrix from file
            if firsttime: # initialize data array
                data = [[0 for col in range(7)] for row in range(len(new_data))] # 7, because 7 columns are needed
            for row in range(len(new_data)):
                if firsttime: 
                    data[row][0] = float(new_data[row][0]) # first column stays always the same (take it from first file)
                    data[row][1] = float(new_data[row][1]) # second column is the central (first file)
                    data[row][2] = float(new_data[row][2]) # third column is the cerror of the entral (first file)
                    data[row][3] = float(new_data[row][1]) # this is the starting point to get the minimu
                    data[row][4] = float(new_data[row][2]) # the error of the minimum
                    data[row][5] = float(new_data[row][1]) # this is the starting point to get the maximum
                    data[row][6] = float(new_data[row][2]) # the error of the maximimum
                else:
                    if(float(new_data[row][1]) < data[row][3]): # if the cross section in the new file is smaller
                        data[row][3] = float(new_data[row][1]) # set the minimum to this value
                        data[row][4] = float(new_data[row][2]) # the error of the minimum
                    if(float(new_data[row][1]) > data[row][5]): # if the cross section in the new file is bigger
                        data[row][5] = float(new_data[row][1]) # set the maximum to this value
                        data[row][6] = float(new_data[row][2]) # the error of the minimum
            firsttime = False
        with open(out_file,'w') as f:
            for row in range(len(data)):
                line = ""
                for column in range(7): # 7, because 7 columns were created
                    length = 15
                    line = line+" "*(length-len("%.8g" % data[row][column]))+"%.8g" % data[row][column]
                f.write(line+"\n")
            f.close()
#}}}    
#{{{ def: readin_file(self,file_path)
    def readin_file(self,file_path):
        # this reads a file and writes it into a table
        with open(file_path,'r') as f:
            table = [row.strip().split() for row in f if not row.split()[0].startswith("#") and not row.split()[0].startswith("%")]
        return table
#}}}    
#{{{ def: seven_point_variation(self,path)
    def seven_point_variation(self,path):
        # this check wether the current folder contains a 7-point variation
        n_scales = 0
        try:
            file_path = glob.glob(pjoin(path,"*","*CV*.dat"))[0]
        except:
            out.print_error("There is no *.dat file in path %s that could tell you what kind of variation was used in routine seven_point_variation." % path)
        with open(file_path) as in_file:
            for line in in_file:
                if line.strip():
                    n_scales += 1
        if n_scales == 7:
            return True
        else:
            return False
#}}}    
#{{{ def: nine_point_variation(self,path)
    def nine_point_variation(self,path):
        # this check wether the current folder contains a 9-point variation
        n_scales = 0
        try:
            file_path = glob.glob(pjoin(path,"*","*CV*.dat"))[0]
        except:
            out.print_error("There is no *.dat file in path %s that could tell you what kind of variation was used in routine nine_point_variation.")
        with open(file_path) as in_file:
            for line in in_file:
                if line.strip():
                    n_scales += 1
        if n_scales == 9:
            return True
        else:
            return False
#}}}    
#{{{ def: seven_point_variation_distribution(self,path)
    def seven_point_variation_distribution(self,path):
        # this check wether the current folder contains a 7-point variation for distributions
        scale_files = glob.glob(pjoin(path,"scale.*"))
        if len(scale_files) == 0:
            out.print_error("There are no scale.* files in path %s that could tell you what kind of variation was used in routine seven_point_variation_distribution" % path)
        n_scales = len(scale_files)
        if n_scales == 7:
            return True
        else:
            return False
#}}}    
#{{{ def: nine_point_variation_distribution(self,path)
    def nine_point_variation_distribution(self,path):
        # this check wether the current folder contains a 9-point variation for distributions
        scale_files = glob.glob(pjoin(path,"scale.*"))
        if len(scale_files) == 0:
            out.print_error("There are no scale.* files in path %s that could tell you what kind of variation was used in routine nine_point_variation_distribution" % path)
        n_scales = len(scale_files)
        if n_scales == 9:
            return True
        else:
            return False
#}}}    
#{{{ def: seven_point_mapping(self,factor)
    def seven_point_mapping(self,factor):
        # this returns an array that contains up and down variation for muR and muF with a factor
        mapping = []
        mapping.append([1./factor,1./factor]) # first row
        mapping.append([1./factor,1.])        # second row
        mapping.append([factor,1.])           # ...
        mapping.append([1.,1.])               # ..
        mapping.append([1.,factor])           # .
        mapping.append([1.,1/factor])       
        mapping.append([factor,factor])       
        return mapping
#}}}    
#{{{ def: nine_point_mapping(self,factor)
    def nine_point_mapping(self,factor):
        # this returns an array that contains up and down variation for muR and muF with a factor
        mapping = []
        mapping.append([1./factor,factor])    # frist row
        mapping.append([1./factor,1./factor]) # second row
        mapping.append([1./factor,1.])        # ...
        mapping.append([factor,1.])           # ..
        mapping.append([1.,1.])               # .
        mapping.append([1.,factor])       
        mapping.append([1.,1./factor])       
        mapping.append([factor,factor])       
        mapping.append([factor,1./factor])       
        return mapping
#}}}    
#{{{ def: convert_to_independent_scales(self,file_path)
    def convert_to_independent_scales(self,file_path):
        # this routine changes the first column of a rate file from muR=muF to independent muR and muF column
        firsttime = True
        new_data = self.readin_file(file_path) # reads space-separated matrix from file
        factor = 1/float(new_data[0][0]) # this is the factor with respect to the central scale
        if firsttime: # initialize data array
            data = [[0 for col in range(len(new_data[0])+1)] for row in range(len(new_data))] # add one column to split up muR and muF
            firsttime = False
        for row in range(len(new_data)):
            for column in range(len(new_data[row])):
                if column == 0:
                    if len(new_data) == 7: # 7-point variation
                        data[row][column]   = self.seven_point_mapping(factor)[row][0] # muR
                        data[row][column+1] = self.seven_point_mapping(factor)[row][1] # muF
                    elif len(new_data) == 9: # 9-point variation1
                        data[row][column]   = self.nine_point_mapping(factor)[row][0] # muR
                        data[row][column+1] = self.nine_point_mapping(factor)[row][1] # muF
                elif column == 1:
                    data[row][column+1] = float(new_data[row][column])
                elif column == 2:
                    data[row][column+1] = float(new_data[row][column])
                else:
                    out.print_error("The rate file converted in convert_to_independent_scales contains more than 3 columns")
        with open(file_path+".replace",'w') as f:
            for row in range(len(data)):
                line = ""
                for column in range(len(data[row])):
                    if column in [0,1]:
                        length = 10
                        line = line+" "*(length-len("%.4g" % data[row][column]))+"%.4g" % data[row][column]
                    else:
                        length = 15
                        line = line+" "*(length-len("%.8g" % data[row][column]))+"%.8g" % data[row][column]
                f.write(line+"\n")
        os.remove(file_path)
        shutil.move(file_path+".replace",file_path)
#}}}    
#{{{ def: extrapolate_qT_dependence(self,file_path)
    def extrapolate_qT_dependence(self,file_path,max_qt_value = 99):
        # this routine reads in a qT_cut-dependence file, extrapolates the central result (either 7 or 9 point variation) to rcut=0 and returns an uncertainty (difference between extrapolated result and at lowest qT_cut value)
        relative_uncertainty = 0 

        firsttime = True
        new_data = self.readin_file(file_path) # reads space-separated matrix from file

        # prepare date to for 7- or 9-point variation
        data_set = [[float(row[0]),float(row[(len(new_data[0])-1)/2])] for row in new_data if abs(float(row[0])) <= abs(max_qt_value)]
        # calculate coefficients of linear fit
        linear_fit = data_point_fitter(data_set)
        b0, b1 = linear_fit.coefficients()
        # take as measure of uncertainty: scope of linear regression times missing interval from lowest rcut value (usually rcut=0.15) to rcut=0 
        # normalized to cross section at lowerst rcut value; in percent; extra factor of 2 to be conservative
        relative_uncertainty = abs(b1/data_set[0][1]*data_set[0][0] * 100) * 2
 
        return relative_uncertainty
#}}}
#}}}
#{{{ class: citations()
class citations(): # class to deal with the results, adding files, copying them to the results folder etc.
#{{{ def: __init__(self,process_in)
    def __init__(self,process_in):
        self.physics_first_list = ["ppemxnmnex04","ppemexnmx04","ppeexmxnm04","ppeeexnex04","ppeexexne04","ppeexnenex04","ppzz02","ppwxw02","ppeeexex04","ppemexmx04","ppeexa03","ppnenexa03","ppenexa03","ppexnea03"]
        self.process = process_in
        self.banner_start = """%          /-------------------------------------------------------------------\\
%          |           __  __     ___     ____    ___         _     _          |
%          |          | _\/_ |   / _ \   |____|  / _ \   ||   \\\\   //          |
%          |          || \/ ||  | |_| |    ||    ||_||   ||    \\\\ //           |
%          |          ||    ||  | ___ |    ||    ||\\\\    ||    // \\\\           |
%          |          ||    ||  ||   ||    ||    || \\\\   ||   //   \\\\          |
%          |                                                                   |
%          |-------------------------------------------------------------------|
%          | MATRIX is based on a number of different computations and tools   |
%          | from various people and groups. Please acknowledge their efforts  |
%          | by citing the references for this run, which are given below.     |
%          \-------------------------------------------------------------------/
"""
        self.banner_end = """
%          /-------------------------------------------------------------------\\
%          | MATRIX is based on a number of different computations and tools   |
%          | from various people and groups. Please acknowledge their efforts  |
%          | by citing the references for this run, which are given above.     |
%          \-------------------------------------------------------------------/
"""
        self.citation_list_standard = self.get_citation_list_standard()
        self.citation_list_process  = self.get_citation_list_process()
        self.citation_list_amplitudes  = self.get_citation_list_amplitudes()
#}}}
#{{{ def: write_citations(self,file_path,citation_list_run, citation_list_standard = self.citation_list_standard,citation_list_process = self.citation_list_process)
    def write_citations(self,file_path,citation_list_run, citation_list_standard = None, citation_list_process = None, citation_list_amplitudes = None):
        # this routine writes out the relevant citations for a process
        if citation_list_standard is None:
            citation_list_standard = self.citation_list_standard
        if citation_list_process is None:
            citation_list_process = self.citation_list_process
        if citation_list_amplitudes is None:
            citation_list_amplitudes = self.citation_list_amplitudes
        with open(file_path,'w') as cite_file:
                cite_file.write(self.banner_start)
                cite_file.write("\n")
                cite_file.write("%-------------\\\n")
                cite_file.write("% Spires keys |\n")
                cite_file.write("%-------------/\n")
                if self.process in self.physics_first_list:
                    if citation_list_process:
                        cite_file.write("% Physics publications:\n")
                    else:
                        out.print_warning("No physics publications available for process \"%s\". Please check yourself wether relevant publications exist." % self.process)
                    for spireskey, bibtex in citation_list_process.iteritems():
                        cite_file.write("%% - %s\n" % spireskey)
                    cite_file.write("%\n")
                    cite_file.write("% MATRIX release publication:\n")
                    for spireskey, bibtex in citation_list_standard.iteritems():
                        cite_file.write("%% - %s\n" % spireskey)                    
                else:
                    cite_file.write("% MATRIX release publication:\n")
                    for spireskey, bibtex in citation_list_standard.iteritems():
                        cite_file.write("%% - %s\n" % spireskey)                    
                    if citation_list_process:
                        cite_file.write("%\n")
                        cite_file.write("% Physics publications:\n")
                    else:
                        out.print_warning("No physics publications available for process \"%s\". Please check yourself wether relevant publications exist." % self.process)
                    for spireskey, bibtex in citation_list_process.iteritems():
                        cite_file.write("%% - %s\n" % spireskey)
                if citation_list_amplitudes:
                    cite_file.write("%\n")
                    cite_file.write("% Amplitude citations:\n")
                for spireskey, bibtex in citation_list_amplitudes.iteritems():
                    cite_file.write("%% - %s\n" % spireskey)                    
                if citation_list_run:
                    cite_file.write("%\n")
                    cite_file.write("% Subtraction-method citations:\n")
                for spireskey, bibtex in citation_list_run.iteritems():
                    cite_file.write("%% - %s\n" % spireskey)                    
                cite_file.write("\n")
                cite_file.write("\n")
                cite_file.write("%----------------\\\n")
                cite_file.write("% BibTeX entries |\n")
                cite_file.write("%----------------/\n")
                if self.process in self.physics_first_list:
                    if citation_list_process:
                        cite_file.write("% Physics publications:\n")
                    for spireskey, bibtex in citation_list_process.iteritems():
                        cite_file.write("%s\n" % bibtex)                    
                        cite_file.write("\n")
                    cite_file.write("% MATRIX release publication:\n")
                    for spireskey, bibtex in citation_list_standard.iteritems():
                        cite_file.write("%s\n" % bibtex)                    
                        cite_file.write("\n")
                else:
                    cite_file.write("% MATRIX release publication:\n")
                    for spireskey, bibtex in citation_list_standard.iteritems():
                        cite_file.write("%s\n" % bibtex)                    
                        cite_file.write("\n")
                    if citation_list_process:
                        cite_file.write("% Physics publications:\n")
                    for spireskey, bibtex in citation_list_process.iteritems():
                        cite_file.write("%s\n" % bibtex)                    
                        cite_file.write("\n")
                if citation_list_amplitudes:
                    cite_file.write("% Amplitude citations:\n")
                for spireskey, bibtex in citation_list_amplitudes.iteritems():
                    cite_file.write("%s\n" % bibtex)                    
                    cite_file.write("\n")
                if citation_list_run:
                    cite_file.write("% Subtraction-method citations:\n")
                for spireskey, bibtex in citation_list_run.iteritems():
                    cite_file.write("%s\n" % bibtex)                    
                    cite_file.write("\n")
                cite_file.write(self.banner_end)
#}}}        
#{{{ def: get_citation_list_standard(self)
    def get_citation_list_standard(self):
        standard_list = {}
        standard_list["Grazzini:2017mhc"] = """
@article{Grazzini:2017mhc,
      author         = "Grazzini, Massimiliano and Kallweit, Stefan and
                        Wiesemann, Marius",
      title          = "{Fully differential NNLO computations with MATRIX}",
      year           = "2017",
      eprint         = "1711.06631",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "CERN-TH-2017-232, ZU-TH-30-17",
      SLACcitation   = "%%CITATION = ARXIV:1711.06631;%%"
}"""
        return standard_list
#}}}        
#{{{ def: get_citation_list_process(self, process = self.process)
    def get_citation_list_process(self, process = None):
        if process is None:
            process = self.process

        ttx_list = {}
	ttx_list["Catani:2019iny"] = """
@article{Catani:2019iny,
      author         = "Catani, Stefano and Devoto, Simone and Grazzini,
                        Massimiliano and Kallweit, Stefan and Mazzitelli, Javier
                        and Sargsyan, Hayk",
      title          = "{Top-quark pair hadroproduction at
                        next-to-next-to-leading order in QCD}",
      journal        = "Phys. Rev.",
      volume         = "D99",
      year           = "2019",
      number         = "5",
      pages          = "051501",
      doi            = "10.1103/PhysRevD.99.051501",
      eprint         = "1901.04005",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH 02/19",
      SLACcitation   = "%%CITATION = ARXIV:1901.04005;%%"
}"""
	ttx_list["Catani:2019hip"] = """
@article{Catani:2019hip,
      author         = "Catani, Stefano and Devoto, Simone and Grazzini,
                        Massimiliano and Kallweit, Stefan and Mazzitelli, Javier",
      title          = "{Top-quark pair production at the LHC: Fully differential
                        QCD predictions at NNLO}",
      journal        = "JHEP",
      volume         = "07",
      year           = "2019",
      pages          = "100",
      doi            = "10.1007/JHEP07(2019)100",
      eprint         = "1906.06535",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH 31/19",
      SLACcitation   = "%%CITATION = ARXIV:1906.06535;%%"
}"""

        Z_list = {}
        Z_list["Catani:2009sm"] = """
@article{Catani:2009sm,
      author         = "Catani, Stefano and Cieri, Leandro and Ferrera, Giancarlo
                        and de Florian, Daniel and Grazzini, Massimiliano",
      title          = "{Vector boson production at hadron colliders: a fully
                        exclusive QCD calculation at NNLO}",
      journal        = "Phys. Rev. Lett.",
      volume         = "103",
      year           = "2009",
      pages          = "082001",
      doi            = "10.1103/PhysRevLett.103.082001",
      eprint         = "0903.2120",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      SLACcitation   = "%%CITATION = ARXIV:0903.2120;%%"
}"""
        W_list = Z_list

        HH_list = {}
        HH_list["deFlorian:2016uhr"] = """
@article{deFlorian:2016uhr,
      author         = "de Florian, Daniel and Grazzini, Massimiliano and Hanga,
                        Catalin and Kallweit, Stefan and Lindert, Jonas M. and
                        Maierh\"ofer, Philipp and Mazzitelli, Javier and Rathlev,
                        Dirk",
      title          = "{Differential Higgs Boson Pair Production at
                        Next-to-Next-to-Leading Order in QCD}",
      journal        = "JHEP",
      volume         = "09",
      year           = "2016",
      pages          = "151",
      doi            = "10.1007/JHEP09(2016)151",
      eprint         = "1606.09519",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "DESY-16-107, FR-PHENO-2016-007, ICAS-08-16, MITP-16-061,
                        ZU-TH-20-16",
      SLACcitation   = "%%CITATION = ARXIV:1606.09519;%%"
}
"""

        gammagamma_list = {}
        gammagamma_list["Catani:2011qz"] = """
@article{Catani:2011qz,
      author         = "Catani, Stefano and Cieri, Leandro and de Florian, Daniel
                        and Ferrera, Giancarlo and Grazzini, Massimiliano",
      title          = "{Diphoton production at hadron colliders: a
                        fully-differential QCD calculation at NNLO}",
      journal        = "Phys. Rev. Lett.",
      volume         = "108",
      year           = "2012",
      pages          = "072001",
      doi            = "10.1103/PhysRevLett.108.072001,
                        10.1103/PhysRevLett.117.089901",
      note           = "[Erratum: Phys. Rev. Lett.117,no.8,089901(2016)]",
      eprint         = "1110.2375",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-18-11, IFUM-984-FT",
      SLACcitation   = "%%CITATION = ARXIV:1110.2375;%%"
}"""

        Zgamma_list = {}
        Zgamma_list["Grazzini:2013bna"] = """
@article{Grazzini:2013bna,
      author         = "Grazzini, Massimiliano and Kallweit, Stefan and Rathlev,
                        Dirk and Torre, Alessandro",
      title          = "{$Z\gamma$ production at hadron colliders in NNLO QCD}",
      journal        = "Phys. Lett.",
      volume         = "B731",
      year           = "2014",
      pages          = "204-207",
      doi            = "10.1016/j.physletb.2014.02.037",
      eprint         = "1309.7000",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-21-13",
      SLACcitation   = "%%CITATION = ARXIV:1309.7000;%%"
}
"""
        Zgamma_list["Grazzini:2015nwa"] = """
@article{Grazzini:2015nwa,
      author         = "Grazzini, Massimiliano and Kallweit, Stefan and Rathlev,
                        Dirk",
      title          = "{$W\gamma$ and $Z\gamma$ production at the LHC in NNLO
                        QCD}",
      journal        = "JHEP",
      volume         = "07",
      year           = "2015",
      pages          = "085",
      doi            = "10.1007/JHEP07(2015)085",
      eprint         = "1504.01330",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-04-15, MITP-15-021",
      SLACcitation   = "%%CITATION = ARXIV:1504.01330;%%"
}"""

        Wgamma_list = {}
        Wgamma_list["Grazzini:2015nwa"] = """
@article{Grazzini:2015nwa,
      author         = "Grazzini, Massimiliano and Kallweit, Stefan and Rathlev,
                        Dirk",
      title          = "{$W\gamma$ and $Z\gamma$ production at the LHC in NNLO
                        QCD}",
      journal        = "JHEP",
      volume         = "07",
      year           = "2015",
      pages          = "085",
      doi            = "10.1007/JHEP07(2015)085",
      eprint         = "1504.01330",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-04-15, MITP-15-021",
      SLACcitation   = "%%CITATION = ARXIV:1504.01330;%%"
}"""

        ZZ_list = {}
        ZZ_list["Cascioli:2014yka"] = """
@article{Cascioli:2014yka,
      author         = "Cascioli, F. and Gehrmann, T. and Grazzini, M. and
                        Kallweit, S. and Maierh\"ofer, P. and von Manteuffel, A.
                        and Pozzorini, S. and Rathlev, D. and Tancredi, L. and
                        Weihs, E.",
      title          = "{ZZ production at hadron colliders in NNLO QCD}",
      journal        = "Phys. Lett.",
      volume         = "B735",
      year           = "2014",
      pages          = "311-313",
      doi            = "10.1016/j.physletb.2014.06.056",
      eprint         = "1405.2219",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-19-14, MITP-14-033",
      SLACcitation   = "%%CITATION = ARXIV:1405.2219;%%"
}"""
        ZZ_list["Grazzini:2015hta"] = """
@article{Grazzini:2015hta,
      author         = "Grazzini, Massimiliano and Kallweit, Stefan and Rathlev,
                        Dirk",
      title          = "{ZZ production at the LHC: fiducial cross sections and
                        distributions in NNLO QCD}",
      journal        = "Phys. Lett.",
      volume         = "B750",
      year           = "2015",
      pages          = "407-410",
      doi            = "10.1016/j.physletb.2015.09.055",
      eprint         = "1507.06257",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      SLACcitation   = "%%CITATION = ARXIV:1507.06257;%%"
}"""

        WW_list = {}
        WW_list["Gehrmann:2014fva"] = """
@article{Gehrmann:2014fva,
      author         = "Gehrmann, T. and Grazzini, M. and Kallweit, S. and
                        Maierh\"ofer, P. and von Manteuffel, A. and Pozzorini, S.
                        and Rathlev, D. and Tancredi, L.",
      title          = "{$W^+W^-$ Production at Hadron Colliders in Next to Next
                        to Leading Order QCD}",
      journal        = "Phys. Rev. Lett.",
      volume         = "113",
      year           = "2014",
      number         = "21",
      pages          = "212001",
      doi            = "10.1103/PhysRevLett.113.212001",
      eprint         = "1408.5243",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-29-14, MITP-14-053",
      SLACcitation   = "%%CITATION = ARXIV:1408.5243;%%"
}"""
        WW_list["Grazzini:2016ctr"] = """
@article{Grazzini:2016ctr,
      author         = "Grazzini, Massimiliano and Kallweit, Stefan and
                        Pozzorini, Stefano and Rathlev, Dirk and Wiesemann,
                        Marius",
      title          = "{$W^{+}W^{-}$ production at the LHC: fiducial cross
                        sections and distributions in NNLO QCD}",
      journal        = "JHEP",
      volume         = "08",
      year           = "2016",
      pages          = "140",
      doi            = "10.1007/JHEP08(2016)140",
      eprint         = "1605.02716",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-16-16, MITP-16-038, NSF-KITP-16-047, DESY-16-075",
      SLACcitation   = "%%CITATION = ARXIV:1605.02716;%%"
}"""
        WZ_list = {}
        WZ_list["Grazzini:2016swo"] = """
@article{Grazzini:2016swo,
      author         = "Grazzini, Massimiliano and Kallweit, Stefan and Rathlev,
                        Dirk and Wiesemann, Marius",
      title          = "{$W^{\pm}Z$ production at hadron colliders in NNLO QCD}",
      journal        = "Phys. Lett.",
      volume         = "B761",
      year           = "2016",
      pages          = "179-183",
      doi            = "10.1016/j.physletb.2016.08.017",
      eprint         = "1604.08576",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-15-16, MITP-16-037, NSF-KITP-16-046, DESY-16-074",
      SLACcitation   = "%%CITATION = ARXIV:1604.08576;%%"
}"""
        WZ_list["Grazzini:2017ckn"] = """
@article{Grazzini:2017ckn,
      author         = "Grazzini, Massimiliano and Kallweit, Stefan and Rathlev,
                        Dirk and Wiesemann, Marius",
      title          = "{$W^\pm Z$ production at the LHC: fiducial cross sections
                        and distributions in NNLO QCD}",
      year           = "2017",
      eprint         = "1703.09065",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-06-17, CERN-TH-2017-065",
      SLACcitation   = "%%CITATION = ARXIV:1703.09065;%%"
}"""
        ZZplusWW_list = ZZ_list.copy()
        ZZplusWW_list.update(WW_list)

        citation_list = {}
#        citation_list["pph21"] = H_list
	citation_list["ppttx20"] = ttx_list
        citation_list["ppz01"] = Z_list
        citation_list["ppw01"] = W_list
        citation_list["ppwx01nockm"] = W_list
        citation_list["ppw01nockm"] = W_list
        citation_list["ppwx01"] = W_list
        citation_list["ppeex02"] = Z_list
        citation_list["ppnenex02"] = Z_list
        citation_list["ppexne02nockm"] = W_list
        citation_list["ppenex02nockm"] = W_list
        citation_list["ppexne02"] = W_list
        citation_list["ppenex02"] = W_list
        citation_list["pphh22"] = HH_list
        citation_list["ppaa02"] = gammagamma_list
        citation_list["ppeexa03"] = Zgamma_list
        citation_list["ppnenexa03"] = Zgamma_list
        citation_list["ppexnea03"] = Wgamma_list
        citation_list["ppenexa03"] = Wgamma_list
        citation_list["ppzz02"] = ZZ_list
        citation_list["ppwxw02"] = WW_list
        citation_list["ppeeexex04"] = ZZ_list
        citation_list["ppemexmx04"] = ZZ_list
        citation_list["ppeexnenex04"] = ZZplusWW_list
        citation_list["ppeexnmnmx04"] = ZZ_list
        citation_list["ppemxnmnex04"] = WW_list
        citation_list["ppemexnmx04"] = WZ_list
        citation_list["ppeexmxnm04"] = WZ_list
        citation_list["ppeeexnex04"] = WZ_list
        citation_list["ppeexexne04"] = WZ_list
        return citation_list.get(process,{})
#}}}        
#{{{ def: get_citation_list_amplitudes(self, process = self.process)
    def get_citation_list_amplitudes(self, process = None):
        if process is None:
            process = self.process
        amplitude_list = {}
	
        amplitude_list["Cascioli:2011va"] = """
@article{Cascioli:2011va,
      author         = "Cascioli, Fabio and Maierh\"ofer, Philipp and Pozzorini,
                        Stefano",
      title          = "{Scattering Amplitudes with Open Loops}",
      journal        = "Phys. Rev. Lett.",
      volume         = "108",
      year           = "2012",
      pages          = "111601",
      doi            = "10.1103/PhysRevLett.108.111601",
      eprint         = "1111.5206",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-23-11, LPN11-66",
      SLACcitation   = "%%CITATION = ARXIV:1111.5206;%%"
}
"""
	amplitude_list["Buccioni:2019sur"] ="""
@article{Buccioni:2019sur,
      author         = "Buccioni, Federico and Lang, Jean-Nicolas and Lindert,
                        Jonas M. and Maierh\"ofer, Philipp and Pozzorini, Stefano
                        and Zhang, Hantian and Zoller, Max F.",
      title          = "{OpenLoops 2}",
      year           = "2019",
      eprint         = "1907.13071",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "IPPP/19/62, FR-PHENO-2019-12, PSI-PR-19-15, ZU-TH 37/19",
      SLACcitation   = "%%CITATION = ARXIV:1907.13071;%%"
}"""
        amplitude_list["Denner:2016kdg"] = """
@article{Denner:2016kdg,
      author         = "Denner, Ansgar and Dittmaier, Stefan and Hofer, Lars",
      title          = "{Collier: a fortran-based Complex One-Loop LIbrary in
                        Extended Regularizations}",
      journal        = "Comput. Phys. Commun.",
      volume         = "212",
      year           = "2017",
      pages          = "220-238",
      doi            = "10.1016/j.cpc.2016.10.013",
      eprint         = "1604.06792",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "FR-PHENO-2016-003, ICCUB-16-016",
      SLACcitation   = "%%CITATION = ARXIV:1604.06792;%%"
}
"""
        if process == "pph21":
            amplitude_list["Harlander:2000mg"] = """
@article{Harlander:2000mg,
      author         = "Harlander, Robert V.",
      title          = "{Virtual corrections to g g ---> H to two loops in the
                        heavy top limit}",
      journal        = "Phys. Lett.",
      volume         = "B492",
      year           = "2000",
      pages          = "74-80",
      doi            = "10.1016/S0370-2693(00)01042-X",
      eprint         = "hep-ph/0007289",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "BNL-HET-00-29",
      SLACcitation   = "%%CITATION = HEP-PH/0007289;%%"
}"""
        elif process in ["ppw01","ppwx01","ppnenex02","ppexne02nockm","ppenex02nockm","ppexne02","ppenex02"]:
            amplitude_list["Matsuura:1988sm"] = """
@article{Matsuura:1988sm,
      author         = "Matsuura, T. and van der Marck, S. C. and van Neerven, W.
                        L.",
      title          = "{The Calculation of the Second Order Soft and Virtual
                        Contributions to the Drell-Yan Cross-Section}",
      journal        = "Nucl. Phys.",
      volume         = "B319",
      year           = "1989",
      pages          = "570-622",
      doi            = "10.1016/0550-3213(89)90620-2",
      reportNumber   = "Print-88-0798 (LEIDEN)",
      SLACcitation   = "%%CITATION = NUPHA,B319,570;%%"
}"""
        elif process == "pphh22":
            amplitude_list["deFlorian:2013uza"] = """
@article{deFlorian:2013uza,
      author         = "de Florian, Daniel and Mazzitelli, Javier",
      title          = "{Two-loop virtual corrections to Higgs pair production}",
      journal        = "Phys. Lett.",
      volume         = "B724",
      year           = "2013",
      pages          = "306-309",
      doi            = "10.1016/j.physletb.2013.06.046",
      eprint         = "1305.5206",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      SLACcitation   = "%%CITATION = ARXIV:1305.5206;%%"
}"""
        elif process == "ppaa02":
            amplitude_list["Anastasiou:2002zn"] = """
@article{Anastasiou:2002zn,
      author         = "Anastasiou, C. and Glover, E. W. Nigel and
                        Tejeda-Yeomans, M. E.",
      title          = "{Two loop QED and QCD corrections to massless fermion
                        boson scattering}",
      journal        = "Nucl. Phys.",
      volume         = "B629",
      year           = "2002",
      pages          = "255-289",
      doi            = "10.1016/S0550-3213(02)00140-2",
      eprint         = "hep-ph/0201274",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "SLAC-PUB-9130, DCPT-02-10, IPPP-02-05, YITP-SB-02-04",
      SLACcitation   = "%%CITATION = HEP-PH/0201274;%%"
}
"""
        elif process in ["ppeexa03","ppnenexa03","ppexnea03","ppenexa03"]:
            amplitude_list["Gehrmann:2011ab"] = """
@article{Gehrmann:2011ab,
      author         = "Gehrmann, Thomas and Tancredi, Lorenzo",
      title          = "{Two-loop QCD helicity amplitudes for $q\bar q \to W^\pm
                        \gamma$ and $q\bar q \to Z^0 \gamma$}",
      journal        = "JHEP",
      volume         = "02",
      year           = "2012",
      pages          = "004",
      doi            = "10.1007/JHEP02(2012)004",
      eprint         = "1112.1531",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-26-11",
      SLACcitation   = "%%CITATION = ARXIV:1112.1531;%%"
}"""
        elif process in ["ppzz02","ppwxw02","ppeeexex04","ppemexmx04","ppeexnenex04","ppeexnmnmx04","ppemxnmnex04","ppemexnmx04","ppeexmxnm04","ppeeexnex04","ppeexexne04"]:
            amplitude_list["Gehrmann:2015ora"] = """
@article{Gehrmann:2015ora,
      author         = "Gehrmann, Thomas and von Manteuffel, Andreas and
                        Tancredi, Lorenzo",
      title          = "{The two-loop helicity amplitudes for $
                        q\overline{q}^{\prime}\to {V}_1{V}_2\to 4 $ leptons}",
      journal        = "JHEP",
      volume         = "09",
      year           = "2015",
      pages          = "128",
      doi            = "10.1007/JHEP09(2015)128",
      eprint         = "1503.04812",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-03-15, MITP-15-011, TTP15-011",
      SLACcitation   = "%%CITATION = ARXIV:1503.04812;%%"
}
"""
        return amplitude_list
#}}}        
#{{{ def: get_citation_list_run(self, order, NLO_subtraction_method)
    def get_citation_list_run(self, order, NLO_subtraction_method):
        run_list = {}
        if "NNLO" in order or ("NLO" in order and NLO_subtraction_method != "1"):            
            run_list["Catani:2007vq"] = """
@article{Catani:2007vq,
      author         = "Catani, Stefano and Grazzini, Massimiliano",
      title          = "{An NNLO subtraction formalism in hadron collisions and
                        its application to Higgs boson production at the LHC}",
      journal        = "Phys. Rev. Lett.",
      volume         = "98",
      year           = "2007",
      pages          = "222002",
      doi            = "10.1103/PhysRevLett.98.222002",
      eprint         = "hep-ph/0703012",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      SLACcitation   = "%%CITATION = HEP-PH/0703012;%%"
}
"""
        if "NNLO" in order:
            if self.process in ["pph21","pphh22"]:
                run_list["Catani:2011kr"] = """
@article{Catani:2011kr,
      author         = "Catani, S. and Grazzini, M.",
      title          = "{Higgs Boson Production at Hadron Colliders:
                        Hard-Collinear Coefficients at the NNLO}",
      journal        = "Eur. Phys. J.",
      volume         = "C72",
      year           = "2012",
      pages          = "2013",
      doi            = "10.1140/epjc/s10052-012-2013-2,
                        10.1140/epjc/s10052-012-2132-9",
      note           = "[Erratum: Eur. Phys. J.C72,2132(2012)]",
      eprint         = "1106.4652",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-12-11",
      SLACcitation   = "%%CITATION = ARXIV:1106.4652;%%"
}
"""
            else:
                run_list["Catani:2012qa"] = """
@article{Catani:2012qa,
      author         = "Catani, Stefano and Cieri, Leandro and de Florian, Daniel
                        and Ferrera, Giancarlo and Grazzini, Massimiliano",
      title          = "{Vector boson production at hadron colliders:
                        hard-collinear coefficients at the NNLO}",
      journal        = "Eur. Phys. J.",
      volume         = "C72",
      year           = "2012",
      pages          = "2195",
      doi            = "10.1140/epjc/s10052-012-2195-7",
      eprint         = "1209.0158",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "ZU-TH-16-12",
      SLACcitation   = "%%CITATION = ARXIV:1209.0158;%%"
}
"""


#         if "NNLO" in order:
#             run_list["Catani:2013tia"] = """
# @article{Catani:2013tia,
#       author         = "Catani, Stefano and Cieri, Leandro and de Florian, Daniel
#                         and Ferrera, Giancarlo and Grazzini, Massimiliano",
#       title          = "{Universality of transverse-momentum resummation and hard
#                         factors at the NNLO}",
#       journal        = "Nucl. Phys.",
#       volume         = "B881",
#       year           = "2014",
#       pages          = "414-443",
#       doi            = "10.1016/j.nuclphysb.2014.02.011",
#       eprint         = "1311.1654",
#       archivePrefix  = "arXiv",
#       primaryClass   = "hep-ph",
#       reportNumber   = "ZU-TH-25-13",
#       SLACcitation   = "%%CITATION = ARXIV:1311.1654;%%"
# }
# """
#         if "NNLO" in order or ("NLO" in order and NLO_subtraction_method == "1"):            
#             run_list["Catani:1996vz"] = """
# @article{Catani:1996vz,
#       author         = "Catani, S. and Seymour, M. H.",
#       title          = "{A General algorithm for calculating jet cross-sections
#                         in NLO QCD}",
#       journal        = "Nucl. Phys.",
#       volume         = "B485",
#       year           = "1997",
#       pages          = "291-419",
#       doi            = "10.1016/S0550-3213(96)00589-5,
#                         10.1016/S0550-3213(98)81022-5",
#       note           = "[Erratum: Nucl. Phys.B510,503(1998)]",
#       eprint         = "hep-ph/9605323",
#       archivePrefix  = "arXiv",
#       primaryClass   = "hep-ph",
#       reportNumber   = "CERN-TH-96-029, CERN-TH-96-29",
#       SLACcitation   = "%%CITATION = HEP-PH/9605323;%%"
# }
# """
        return run_list
#}}}
#}}}        
#{{{ class: data_point_fitter()
class data_point_fitter(): # class to deal with the fitting of data points, first implementation: linear regression
#{{{ def: __init__(self,data_points)
    def __init__(self,data_points):
        self.data_points = data_points
#}}}
#{{{ def: mean(values)
    def mean(self,values):
    # Calculate the mean value of a list of numbers
        return sum(values) / float(len(values))
#}}}
#{{{ def: covariance(x, mean_x, y, mean_y)
    def covariance(self,x, mean_x, y, mean_y):
    # Calculate covariance between x and y
        covar = 0.0
        for i in range(len(x)):
            covar += (x[i] - mean_x) * (y[i] - mean_y)
        return covar
#}}}
#{{{ def: variance(values, mean)
    def variance(self,values, mean):
    # Calculate the variance of a list of numbers
        return sum([(x-mean)**2 for x in values])
#}}}
#{{{ def: coefficients(dataset)
    def coefficients(self,mode="linear"):
    # Calculate coefficients
        if mode == "linear":
            x = [row[0] for row in self.data_points]
            y = [row[1] for row in self.data_points]
            x_mean, y_mean = self.mean(x), self.mean(y)
            b1 = self.covariance(x, x_mean, y, y_mean) / self.variance(x, x_mean)
            b0 = y_mean - b1 * x_mean
        else:
            out.print_error("Other fitting procedures than \"linear\" not implemented...")
        return [b0, b1]
#}}}
#}}}        
