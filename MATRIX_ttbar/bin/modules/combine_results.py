#!/usr/bin/env python

import numpy as np
import os
import sys
import copy

print_debug=False
print_channel_info=True

def printDebug(line):
  if print_debug:
    print(line)

class ResultData:
  def __init__(self,n_scales,cutnumber):
    self.n_scales=n_scales
    self.cutnumber=cutnumber
    self.N_tot=0
    self.bin_N=np.zeros((n_scales,cutnumber))
    self.bin_values=np.zeros((n_scales,cutnumber))
    self.bin_deviations=np.zeros((n_scales,cutnumber))
    self.values=np.zeros((n_scales,cutnumber))
    self.sigmas=np.zeros((n_scales,cutnumber))
  
  def readinData(self,filename):
    try:
      f = open(filename,'r')
    except:
      printDebug("could not open "+filename)
      return
    printDebug("reading in "+filename)
    lines = f.readlines()
    f.close()
    
    cut_dep = False
    if len(lines)==self.cutnumber*self.n_scales*2+3+4:
      cut_dep = True
      cutnumber = (len(lines)-7)/2/self.n_scales
    elif len(lines)==self.n_scales*2+3+4:
      cutnumber = 1
    else:
      print(filename+" is misshaped (length is "+np.str(len(lines))+", expected "+np.str(self.cutnumber*self.n_scales*2+3+4)+" or "+np.str(self.n_scales*2+3+4)+"), setting to zero")
      return


    printDebug(np.str(cutnumber)+" bins in file "+filename)

    self.N_tot=int(lines[0])
    for j in range(self.n_scales):
      for i in range(self.cutnumber):
        if cut_dep:
          self.values[j,i]=np.float(lines[3+2*i*self.n_scales+2*j])
          self.sigmas[j,i]=np.float(lines[3+2*i*self.n_scales+2*j+1])
        else:
          self.values[j,i]=np.float(lines[3+2*0*self.n_scales+2*j])
          self.sigmas[j,i]=np.float(lines[3+2*0*self.n_scales+2*j+1])

        # we now need to reconstruct the sum of weights in each bin
        self.bin_values[j,i]=self.N_tot*self.values[j,i]
        self.bin_deviations[j,i]=np.square((self.N_tot-1)*self.sigmas[j,i])+np.square(self.bin_values[j,i])/self.N_tot

  def printData(self):
    for i in range(len(self.bin_N)):
      print(str(i)+": " + np.str(self.bin_N[i])+", "+np.str(self.bin_values[i])+", "+np.str(self.bin_deviations[i]))

class ResultCombination:
  def __init__(self, n_scales, cutnumber, mincut, cutstep):
    self.n_scales=n_scales
    self.cutnumber=cutnumber
    self.contributions={}
    self.orders={}
    self.mincut=mincut
    self.cutstep=cutstep
  
  def readin(self):
    for c in self.contributions:
      self.contributions[c]["results"] = {}
      for directory in self.contributions[c]["dirs"]:
        if not os.path.exists("../"+directory):
          print("../"+directory+" does no exist, skipping folder.")
          assert(False)
          continue

        channels = os.listdir("../"+directory+"/log/")
        for channel in channels:
          printDebug(channel)
          if len(channel)>3 and channel[len(channel)-3:len(channel)]==".in":
            channel_name=channel[0:-3]
            printDebug(channel_name)
            if channel_name not in self.contributions[c]["results"]:
              self.contributions[c]["results"][channel_name] = []
            
            data = ResultData(self.n_scales,self.cutnumber)
            
            filename_in="../"+directory+"/result/result_"+channel_name+".dat"
            #try:
            data.readinData(filename_in)
            self.contributions[c]["results"][channel_name].append(data)
            #except:
              #print("could not readin from "+filename_in)
             # continue
            
            #self.contributions[0][c]["results"][channel_name][0].printData()
          
  def printResults(self):
    print(self.name)
    for c in self.contributions:
      print (c)
      for channel in self.contributions[c]["results"]:
        print(channel)
        for data in self.contributions[c]["results"][channel]:
          data.printData()

  def average(self):
    printDebug("averaging..")
    for c in self.contributions:
      self.contributions[c]["average"]={}
      for channel in self.contributions[c]["results"]:
        averageData=ResultData(self.n_scales,self.cutnumber)
        for data in self.contributions[c]["results"][channel]:
          averageData.N_tot += data.N_tot
          for scale in range(self.n_scales):
            printDebug(c + ", " +channel+", scale="+np.str(scale))
            
            # first average, then compute
            averageData.bin_values[scale] += data.bin_values[scale]
            averageData.bin_deviations[scale] += data.bin_deviations[scale]
          #print("averaged:")
          #averageData.printData()
        self.contributions[c]["average"][channel] = averageData

  def compute(self):
    for c in self.contributions:
      self.contributions[c]["final"]={}
      result_all_channels=ResultData(self.n_scales,self.cutnumber)
      for channel in self.contributions[c]["average"]:
        result = ResultData(self.n_scales,self.cutnumber)

        # first average subsets of the runs conservatively to gain sufficient statistics. The perform a weighted average afterwards to smoothen the histograms
        data_preaveraged = []
        data_tmp = ResultData(self.n_scales,self.cutnumber)
        num_av_datas = 0
        num_to_av = len(self.contributions[c]["results"][channel])/self.contributions[c]["averaging_factor"]

        for data in self.contributions[c]["results"][channel]:
          data_tmp.N_tot += data.N_tot
          for scale in range(self.n_scales):
            data_tmp.bin_values[scale] += data.bin_values[scale]
            data_tmp.bin_deviations[scale] += data.bin_deviations[scale]
          num_av_datas += 1
          if num_av_datas >= num_to_av:
            data_preaveraged.append(data_tmp)
            num_av_datas = 0
            data_tmp = ResultData(self.n_scales,self.cutnumber)
        if num_av_datas > 0:
          data_preaveraged.append(data_tmp)

        # now compute cross section for pre-averaged distributions
        for scale in range(self.n_scales):
          for data in data_preaveraged:
            if data.N_tot==0:
              printDebug(c+", "+channel+" seems to be empty, setting to zero")
              data.values[scale]=np.zeros(self.cutnumber)
              data.sigmas[scale]=np.zeros(self.cutnumber)
            else:
              data.values[scale]=data.bin_values[scale]/data.N_tot
              data.sigmas[scale]=np.sqrt((data.bin_deviations[scale]-np.square(data.bin_values[scale])/data.N_tot))/(data.N_tot-1)

        # first compute, then average
        for scale in range(self.n_scales):
          for data in data_preaveraged:
          #print(data.values)
          #print(len(data_preaveraged))
            with np.errstate(divide='ignore', invalid='ignore'):
              result.values[scale] += np.where(data.sigmas[scale] != 0., data.values[scale] / np.square(data.sigmas[scale]), 0)

            with np.errstate(divide='ignore', invalid='ignore'):
              result.sigmas[scale] += np.where(data.sigmas[scale] != 0., 1.0/np.square(data.sigmas[scale]), 0)
            
          with np.errstate(divide='ignore', invalid='ignore'):
            result.values[scale] = np.where(result.sigmas[scale] != 0., result.values[scale] / result.sigmas[scale], 0)
          
          with np.errstate(divide='ignore', invalid='ignore'):
            result.sigmas[scale] = 1.0/np.sqrt(result.sigmas[scale])
          
          result.sigmas[scale][np.isinf(result.sigmas[scale])] = 0
        
        #print(np.sum(result.values))
        
        
        if (np.isinf(np.sum(result.values))):
          print(result.values)
          print(result.sigmas)
          print(2.0*(result.sigmas==0))
          sys.exit(0)
        if (np.isinf(np.sum(result.sigmas))):
          print("inf")
          print(result.values)
          print(result.sigmas)
          #sys.exit(0)

        self.contributions[c]["final"][channel] = result

        for scale in range(self.n_scales):
          result_all_channels.values[scale] += result.values[scale]
          result_all_channels.sigmas[scale] = np.sqrt(np.square(result_all_channels.sigmas[scale])+np.square(result.sigmas[scale]))

        self.contributions[c]["finalsum"] = result_all_channels

  def printResult(self):
    for scale in range(self.n_scales):
      for c in self.contributions[scale]:
        for channel in self.contributions[scale][c]["final"]:
          print(self.name+": "+c+", "+channel)
          self.contributions[scale][c]["final"][channel].printDist()
        print(self.name+": "+c)
        self.contributions[scale][c]["finalsum"].printDist()
       
      for order in self.orders[scale]:
        print(self.name+": "+order)
        self.orders[scale][order].printDist()
       
  def computeOrders(self):
    print("computing final result")
    for c in self.contributions:
      order=c.split("_")[0]
      type_correction=c.split("_")[1]
      printDebug(order+", "+c)
      if order not in self.orders:
        self.orders[order]=ResultData(self.n_scales,self.cutnumber)
      
      for scale in range(self.n_scales):
        self.orders[order].values[scale] += self.contributions[c]["finalsum"].values[scale]
        self.orders[order].sigmas[scale] = np.sqrt(np.square(self.orders[order].sigmas[scale])+np.square(self.contributions[c]["finalsum"].sigmas[scale]))

    for order in self.orders:
      print(order+" CV")
      for scale in range(self.n_scales):
        print(np.str(scale)+": "+np.str(self.orders[order].values[scale][0])+" +/- "+np.str(self.orders[order].sigmas[scale][0]))

    # compute envelope of scale variations
#    for order in self.orders:
#      self.envelopes_up[order] = copy.deepcopy(self.orders[0][order])
#      self.envelopes_down[order] =  copy.deepcopy(self.orders[0][order])
#      
#      for i in range(self.binnumber):
#        for scale in range(1,self.n_scales):
#          if self.orders[scale][order].values[i] > self.envelopes_up[order].values[i]:
#            self.envelopes_up[order].values[i] = self.orders[scale][order].values[i]
#            self.envelopes_up[order].sigmas[i] = self.orders[scale][order].sigmas[i]
#          if self.orders[scale][order].values[i] < self.envelopes_down[order].values[i]:
#            self.envelopes_down[order].values[i] = self.orders[scale][order].values[i]
#            self.envelopes_down[order].sigmas[i] = self.orders[scale][order].sigmas[i]

  def saveDist(self,filename,dist):
    printDebug("saving "+filename)
    f = open(filename,'w')
    for i in range(self.cutnumber):
#      if self.binningtype=="linear":
#        edge = i*self.binwidth+self.startpoint
#        binwidth = self.binwidth
#      elif self.binningtype=="logarithmic":
#        step = np.log10(self.endpoint/self.startpoint)/self.binnumber
#        edge = np.power(10,np.log10(self.startpoint)+i*step)
#        binwidth = np.power(10,np.log10(self.startpoint)+(i+1)*step)-np.power(10,np.log10(self.startpoint)+i*step)
#      else:
#        print("Binning type "+self.binningtype+" not supported!")
#        assert(False)
#
#      assert(binwidth == self.binwidths[i])
      line_to_write = np.str(i*self.cutstep+self.mincut)
      for scale in range(self.n_scales):
        line_to_write = line_to_write + '\t' + np.str(dist.values[scale][i]) + '\t' + np.str(dist.sigmas[scale][i])
      
      f.write(line_to_write+'\n')
      #f.write(str(i])+'\t'+str((dist.values/self.binwidths)[i])+'\t'+str((dist.sigmas/self.binwidths)[i])+'\n')

    f.close()

  def save(self,resultdirectory):
    for c in self.contributions:
      contribution=c.split("_")[2]
      type_correction=c.split("_")[1]
      if type_correction=="---":
        type_correction=""
      else:
        type_correction = '.'+type_correction
      printDebug("saving "+c)

      if print_channel_info:
        for channel in self.contributions[c]["final"]:
          output_dir=resultdirectory+"/"+self.contributions[c]["output_dir"]+"/"+contribution+type_correction+"/channels"
          if not os.path.exists(output_dir):
            os.makedirs(output_dir)

          self.saveDist(output_dir+"/plot."+contribution+"."+channel+".dat",self.contributions[c]["final"][channel])

      # output of contribution wise distributions
      output_dir=resultdirectory+"/"+self.contributions[c]["output_dir"]+"/"+contribution+type_correction
      if not os.path.exists(output_dir):
        os.makedirs(output_dir)

      filename=output_dir+"/plot."+self.contributions[c]["output_dir"]+"."+contribution+type_correction+".dat"
      self.saveDist(filename,self.contributions[c]["finalsum"])
      
      #filename=output_dir+"/plot."+self.name+"."+self.contributions[scale][c]["output_dir"]+"."+contribution+type_correction+".cons.dat"
      #self.saveDist(filename,self.contributions[scale][c]["finalsum_cons"])
#      printDebug("integral: "+np.str(np.sum(self.contributions[scale][c]["finalsum"].values))+" +/- "+np.str(np.sqrt(np.sum(np.square(self.contributions[scale][c]["finalsum"].sigmas)))))
      
    for order in self.orders:
      if order in ["NLO","NNLO"]:
        filename=resultdirectory+"/plot.qTcut."+order+".QCD.dat"
      else:
        filename=resultdirectory+"/plot.qTcut."+order+".dat"
      self.saveDist(filename,self.orders[order])

# save CV files
      if order in ["NLO","NNLO"]:
        filename=resultdirectory+"/plot.CV."+order+".QCD.dat"
      else:
        filename=resultdirectory+"/plot.CV."+order+".dat"
      f = open(filename,'w')
      for scale in range(self.n_scales):
        f.write(np.str(scale)+'\t'+np.str(self.orders[order].values[scale][0])+'\t'+np.str(self.orders[order].sigmas[scale][0])+'\n')
      f.close()

      #print("integral: "+np.str(np.sum(self.orders[scale][order].values)*self.binwidth)+" +/- "+np.str(np.sqrt(np.sum(np.square(self.orders[scale][order].sigmas)))*self.binwidth))
      #f = open(resultdirectory+"_py/scale."+str(scale)+"/plot."+self.name+"."+order+".dat",'w')
      #for i in range(self.binnumber):
        #f.write(str(i*self.binwidth+self.startpoint)+'\t'+str(self.orders[scale][order].values[i])+'\t'+str(self.orders[scale][order].sigmas[i])+'\n')
      #f.write(str(self.endpoint)+'\t'+str(self.orders[scale][order].values[-1])+'\t'+str(self.orders[scale][order].sigmas[-1])+'\n')
      #f.close()
    
  def printSpecs(self):
    print("name="+self.name)
    print("startpoint="+str(self.startpoint))
    print("endpoint="+str(self.endpoint))
    print("binwidth="+str(self.binwidth))
    print("binnumber="+str(self.binnumber))

  def cleanup(self):
    del(self.contributions)


def parseContributionFile(filename, contributions, averaging_factor):
  f = open(filename,'r')
  for line in f:
    variable = ""
    value = ""
    variable,value = parseLine(line)
    if value=="":
      continue
    printDebug("parsing: "+ variable+" = " + value)
    if variable=="type_correction":
      printDebug("new contribution: " +value)
      contribution_name=order+"_"+value+"_"+type_contribution
      contributions[contribution_name] = {}
      contributions[contribution_name]["input_dirs"] = []
      contributions[contribution_name]["output_dir"] = output_dir # should have been encountered before
      contributions[contribution_name]["averaging_mode"] = "hybrid" # use hybrid averaging by default
      contributions[contribution_name]["averaging_factor"] = averaging_factor # use default from main infile
    elif variable=="type_contribution":
      type_contribution = value
    elif variable=="directory":
      contributions[contribution_name]["input_dirs"].append(value)
    elif variable=="resultdirectory":
      output_dir = value
    elif variable=="type_perturbative_order":
      order = value
    elif variable=="type_correction":
      type_correction = value
#    elif variable=="averaging_mode":
#      contributions[contribution_name]["averaging_mode"] = value
    elif variable=="averaging_factor":
      contributions[contribution_name]["averaging_factor"] = value


  f.close()

  return

def parseInfileDistribution(filename, contributions):
  contribution_name = ""
  averaging_factor = 4 # average down each contribution to four runs by default
  f = open(filename,'r')
  for line in f:
    variable,value = parseLine(line)
    if value=="":
      continue
    printDebug("parsing: "+ variable+" = " + value)
    if variable=="final_resultdirectory":
      resultDirectory = value
    if variable=="averaging_factor" or variable=="average_factor":
      averaging_factor = np.int(value)
    elif variable=="contribution_file":
      printDebug("parsing contribution file "+value)
      parseContributionFile(value,contributions,averaging_factor)
    else:
      print("do not understand >"+variable+"<")

  f.close()

  return resultDirectory

def parseFileParameter(filename):
  n_scales_CV = 0
  n_qTcut = 0
  min_qTcut = 0
  step_qTcut = 0

  f = open(filename,'r')
    #line = f.read()
  for line in f:
    if len(line)==1:
      continue
    if line[0]=='#':
      continue;
    #print(line)
    if len(line.split())<2:
      continue
    variable = line.split()[0]
    value = line.split()[2]
    #print("var: "+variable)
    #print("val: "+value)
    #print("result: "+variable+" = " + value)
    if variable=="n_scales_CV":
      printDebug("n_scales_CV = "+value)
      n_scales_CV = int(value)

    if variable=="n_qTcut":
      printDebug("n_qTcut = "+value)
      n_qTcut = int(value)

    if variable=="min_qTcut":
      printDebug("min_qTcut = "+value)
      min_qTcut = float(value)

    if variable=="step_qTcut":
      printDebug("step_qTcut = "+value)
      step_qTcut = float(value)

  f.close()
  
  return n_scales_CV,n_qTcut,min_qTcut,step_qTcut

def parseLine(line):
  printDebug("parsing line "+line)
  if len(line)==1:
    return "",""
  if line[0]=='#':
    return "",""

  variable = line.split()[0]
  printDebug("variable = >"+variable+"<")
  value = line.split("=")[1]
  value = value.split()[0]
  printDebug("value = >"+value+"<")
  return variable,value

def combine_results(directory,infile):
  os.chdir(directory)
  n_scales_CV,n_qTcut,step_qTcut,min_qTcut = parseFileParameter('../file_parameter.dat')

  resultDirectory=""
  contributions = {}
  resultDirectory = parseInfileDistribution(infile, contributions)
  
  if resultDirectory=="":
    print("result directory is not set, something went wrong!")
    exit()
  
  if not os.path.exists(resultDirectory+"_py"):
    os.makedirs(resultDirectory+"_py")

  if not os.path.exists(resultDirectory+"_py/CV/"):
    os.makedirs(resultDirectory+"_py/CV")
  
  for c in contributions:
    printDebug(c+": ")
    printDebug(contributions[c])
    
  results = ResultCombination(n_scales_CV,n_qTcut,step_qTcut,min_qTcut)

  for c in contributions:
    for scale in range(n_scales_CV):
      new_c={}
      new_c["dirs"]=contributions[c]["input_dirs"]
      if not new_c["dirs"]:
        print("no input directories specified for contribution "+c+", aborting.")
        assert(False)

      new_c["output_dir"]=contributions[c]["output_dir"]
      new_c["averaging_mode"]=contributions[c]["averaging_mode"]
      new_c["averaging_factor"]=contributions[c]["averaging_factor"]
      results.contributions[c] = new_c
#      printDebug(distributions[d].contributions[scale][c])

  results.readin()
#    #distributions[d].printDist()
  results.average()
  results.compute()
  results.computeOrders()
#    #distributions[d].printResult()
#    print("writing results to "+resultDirectory+"_py/CV")
  results.save(resultDirectory+"_py/CV")
#    distributions[d].cleanup()

if __name__=="__main__":
    infile = sys.argv[1]
    current_dir = os.getcwd()
    combine_results(current_dir,infile)
