#{{{ imports
import os
import subprocess
import glob
import tarfile
import urllib
from os.path import join as pjoin
# own modules
from initialize_classes import out
#}}}

#{{{ class: lhapdf()
class lhapdf(): # class that takes care of PDF sets
#{{{ def: __init__(self,path,lo_set_in,nlo_set_in,nnlo_set_in,,order_in)
    def __init__(self,config_list,lo_set_in,nlo_set_in,nnlo_set_in,order_in):
        try: # first use path from MATRIX_configuration
            self.lhapdf_config = config_list["path_to_lhapdf"]
        except:
            try: # otherwise try using which to identify path to lhapdf-config
                self.lhapdf_config = self.get_lhapdf_config_path()
            except:
                try: # sometime which does not work, test wether simply using lhapdf-config directly works
                    self.lhapdf_config = "lhapdf-config"
                    subprocess.Popen([self.lhapdf_config,"--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip()
                except:
                    out.print_error("Path to \"lhapdf-config\" is not set in MATRIX_configuration and does not appear to be configured as executable. Set \"path_to_lhapdf\" in MATRIX_configuration file or add path to \"lhapdf-config\" to your environmental $PATH variable.")
        self.lhapdf_version = self.get_lhapdf_version(self.lhapdf_config) # this is also a check wether lhapdf-config works fine
        self.pdf_sets_path    = self.get_pdf_sets_path(self.lhapdf_config)
        self.lhapdf_get_data = self.get_lhpdf_get_data_path(self.lhapdf_config) # in 6 this is take care of by the "lhapdf" executable, which is set in that case
        if self.lhapdf_version.startswith("6."):
            self.installed_sets = [os.path.basename(x) for x in glob.iglob(pjoin(self.pdf_sets_path,"*")) if os.path.isdir(x)]
        elif self.lhapdf_version.startswith("5."):
            self.installed_sets = [os.path.splitext(os.path.basename(x))[0] for x in glob.iglob(pjoin(self.pdf_sets_path,"*.LHgrid"))]
        else:
            out.print_error("LHAPDF version gives %s, but must be either 5 or 6. Check wether your \"lhapdf-config\" is the right one." % self.lhapdf_version)
        self.lhapdf_set = {}
        self.lhapdf_set["LO"]   = lo_set_in
        self.lhapdf_set["NLO"]  = nlo_set_in
        self.lhapdf_set["NNLO"] = nnlo_set_in
        self.order = order_in
#}}}
#{{{ def: get_lhapdf_config_path(self)
    def get_lhapdf_config_path(self):
        lhapdf_config_path = subprocess.Popen(["which","lhapdf-config"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip()
        if lhapdf_config_path == "": # if empty string cause exception
            lhapdf_config_path = 1/lhapdf_config_path
        return lhapdf_config_path
#}}}    
#{{{ def: get_lhapdf_version(self,lhapdf_config_path)
    def get_lhapdf_version(self,lhapdf_config_path):
        try:
            version = subprocess.Popen([lhapdf_config_path,"--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip()
        except:
            out.print_error("Cannot extract version from path to \"lhapdf-config\": %s" % lhapdf_config_path)
        return version
#}}}    
#{{{ def: print_lhapdf_version(self,lhapdf_config_path)
    def print_lhapdf_version(self):
        out.print_info("Using LHAPDF version %s..." % self.lhapdf_version )
#}}}    
#{{{ def: get_pdf_sets_path(self,lhapdf_config_path)
    def get_pdf_sets_path(self,lhapdf_config_path):
        if self.lhapdf_version.startswith("6."):
            return self.get_pdf_sets_path_6(lhapdf_config_path)
        elif self.lhapdf_version.startswith("5."):
            return self.get_pdf_sets_path_5(lhapdf_config_path)
        else:
            out.print_error("LHAPDF version gives %s, but must be either 5 or 6. Check wether your \"lhapdf-config\" is the right one." % self.lhapdf_version)
#}}}    
#{{{ def: get_pdf_sets_path_5(self,lhapdf_config_path)
    def get_pdf_sets_path_5(self,lhapdf_config_path):
        pdf_sets_path = subprocess.Popen([lhapdf_config_path,"--pdfsets-path"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip()
        return pdf_sets_path
#}}}    
#{{{ def: get_pdf_sets_path_6(self,lhapdf_config_path)
    def get_pdf_sets_path_6(self,lhapdf_config_path):
        pdf_sets_path = subprocess.Popen([lhapdf_config_path,"--datadir"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip()
        return pdf_sets_path
#}}}    
#{{{ def: get_lhpdf_get_data_path(self,lhapdf_config_path)
    def get_lhpdf_get_data_path(self,lhapdf_config_path):
        if self.lhapdf_version.startswith("6."):
            return self.get_lhpdf_get_data_path_6(lhapdf_config_path)
        elif self.lhapdf_version.startswith("5."):
            return self.get_lhpdf_get_data_path_5(lhapdf_config_path)
        else:
            out.print_error("LHAPDF version gives %s, but must be either 5 or 6. Check wether your \"lhapdf-config\" is the right one." % self.lhapdf_version)
#}}}    
#{{{ def: get_lhpdf_get_data_path_6(self,lhapdf_config_path)
    def get_lhpdf_get_data_path_5(self,lhapdf_config_path):
        try:
            lhapdf_get_data_path = subprocess.Popen(["which","lhapdf-getdata"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip()
        except:
            try: 
                lhapdf_get_data_path = lhapdf_config_path.rsplit("-",1)[0]+"-getdata"
            except:
                out.print_error("Could not determine path to \"lhapdf-getdata\" file. Does not seem to be in the same folder as \"lhapdf-config\" nor in the environmental $PATH variable.")
        return lhapdf_get_data_path
#}}}    
#{{{ def: get_lhpdf_get_data_path_6(self,lhapdf_config_path)
    def get_lhpdf_get_data_path_6(self,lhapdf_config_path):
        try:
            lhapdf_get_data_path = subprocess.Popen(["which","lhapdf"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip()
        except:
            try: 
                raise ValueError('A very specific bad thing happened')
                lhapdf_get_data_path = lhapdf_config_path.rsplit("-",1)[0]
            except:
                out.print_error("Could not determine path to \"lhapdf-getdata\" file. Does not seem to be in the same folder as \"lhapdf-config\" nor in the environmental $PATH variable.")
        return lhapdf_get_data_path
#}}}    
#{{{ def: get_missing_sets(self)
    def get_missing_sets(self):
        missing_sets = []
        for this_order in self.order:
            if not self.lhapdf_set[this_order] in self.installed_sets:
                missing_sets.append(self.lhapdf_set[this_order])
        return missing_sets
#}}}
#{{{ def: download_pdf_sets(self)
    def download_pdf_sets(self,missing_sets):
        if self.lhapdf_version.startswith("6."):
            return self.download_pdf_sets_6(missing_sets)
        elif self.lhapdf_version.startswith("5."):
            return self.download_pdf_sets_5(missing_sets)
        else:
            out.print_error("LHAPDF version gives %s, but must be either 5 or 6. Check wether your \"lhapdf-config\" is the right one." % self.lhapdf_version)
#}}}
#{{{ def: download_pdf_sets_5(self)
    def download_pdf_sets_5(self,missing_sets):
        out.print_info("Downloading missing LHAPDF sets...")
        for missing_set in missing_sets:
            if "NNPDF30" in missing_set:
#                try:
                    self.download_set_by_hand("http://pcteserver.mi.infn.it/~nnpdf/nnpdf30/",missing_set+".LHgrid.tgz")
                    out.print_info("LHAPDF set \"%s.LHgrid\" successfully downloaded." % missing_set)
#                except:
#                    out.print_error("PDF set \"%s\" is no standard set in LHAPDF 5. Downloading it directly from the developer's website failed. Try to download it by hand or change PDF set or upgrade to LHAPDF 6, and restart code!" % missing_set)
            elif "MMHT2014" in missing_set:
                try:
                    self.download_set_by_hand("https://www.hep.ucl.ac.uk/mmht/LHAPDF5/",missing_set+".LHgrid")
                    out.print_info("LHAPDF set \"%s.LHgrid\" successfully downloaded." % missing_set)
                except:
                    out.print_error("PDF set \"%s\" is no standard set in LHAPDF 5. Downloading it directly from the developer's website failed. Try to download it by hand or change PDF set or upgrade to LHAPDF 6, and restart code!" % missing_set)
            else:
                download_set = subprocess.Popen(["lhapdf-getdata","--dest="+self.pdf_sets_path,missing_set+".LHgrid"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[1].strip().rsplit("\n",1)[1]
                print "LHAPDF: ",download_set
                if "Getting PDF set from" in download_set:
                    out.print_info("LHAPDF set \"%s.LHgrid\" successfully downloaded." % missing_set)
                elif "No sets match the arguments given" in download_set:
                    out.print_error("Failed to download LHAPDF set \"%s\" with \"lhapdf-getdata\" script. Try installing it manually to folder %s and restart the code." % (missing_set+".LHgrid",self.pdf_sets_path))
                else:
                    out.print_warning("LHAPDF printout \"%s\" not recognized. Stop the code manually if set \"%s\" was not correctly installed in folder %s." % (download_set,missing_set+".LHgrid",self.pdf_sets_path))
#}}}
#{{{ def: download_set_by_hand(self,url_link,PDF_set)
    def download_set_by_hand(self,url_link,pdf_set):
        cwd = os.getcwd()
        pdfsets_path = subprocess.Popen(["lhapdf-config","--pdfsets-path"],stdout=subprocess.PIPE).communicate()[0].strip()
        os.chdir(pdfsets_path)
        pdf_file = urllib.URLopener()
        pdf_file_url = url_link+pdf_set
        pdf_file.retrieve(pdf_file_url,pdf_set)
        if pdf_set.endswith(".tgz") or pdf_set.endswith(".tar")or pdf_set.endswith(".tar.gz"):
            tar = tarfile.open(pdf_set)
            tar.extractall()
            tar.close()
            os.remove(pdf_set)
        os.chdir(cwd)
#}}}
#{{{ def: download_pdf_sets_6(self)
    def download_pdf_sets_6(self,missing_sets):
        out.print_info("Downloading missing LHAPDF sets...")
        for missing_set in missing_sets:
            download_set = filter(None, subprocess.Popen(["lhapdf","install",missing_set], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate())
            print download_set
            try:
                download_set = list(download_set)
            except:
                pass
            try:
                download_set = (download_set+["noreturn"])[0].strip()
            except:
                pass
            if self.lhapdf_version.startswith("6."):
                installed_sets = [os.path.basename(x) for x in glob.iglob(pjoin(self.pdf_sets_path,"*")) if os.path.isdir(x)]
            elif self.lhapdf_version.startswith("5."):
                installed_sets = [os.path.splitext(os.path.basename(x))[0] for x in glob.iglob(pjoin(self.pdf_sets_path,"*.LHgrid"))]
            else:
                out.print_error("LHAPDF version gives %s, but must be either 5 or 6. Check wether your \"lhapdf-config\" is the right one." % self.lhapdf_version)
            if download_set.startswith("%s.tar.gz" % missing_set) or missing_set in installed_sets:
                out.print_info("LHAPDF set \"%s.LHgrid\" successfully downloaded." % missing_set)
            elif "WARNING: No matching PDFs for pattern:" in download_set:
                print "LHAPDF",download_set
                out.print_error("Failed to download LHAPDF set \"%s\" with \"lhapdf get\" script. Try installing it manually to folder %s and restart the code." % (missing_set+".LHgrid",self.pdf_sets_path))
            else:
                print "LHAPDF",download_set
                out.print_warning("LHAPDF printout \"%s\" not recognized. Stop the code manually if set \"%s\" was not correctly installed in folder %s." % (download_set,missing_set+".LHgrid",self.pdf_sets_path))
#}}}
#}}}
