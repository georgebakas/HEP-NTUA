#!/usr/bin/env python

import numpy as np
import os
import sys
import copy

print_debug=False
print_channel_info=False
remove_execution_files=False

#{{{ def: printDebut(line)
def printDebug(line):
  if print_debug:
    print(line)
#}}}
#{{{ class: DistributionData
class DistributionData:
  def __init__(self,binnumber):
    self.N_tot=0
    self.bin_N=np.zeros(binnumber)
    self.bin_values=np.zeros(binnumber)
    self.bin_deviations=np.zeros(binnumber)
    self.values=np.zeros(binnumber)
    self.sigmas=np.zeros(binnumber)
  
  def readinData(self,filename):
    try:
      f = open(filename,'r')
    except:
      #print("could not open "+filename)
      return
    printDebug("reading in "+filename)
    lines = f.readlines()
    f.close()
    
    binnumber = (len(lines)-1)/3
    printDebug(np.str(binnumber)+" bins in file "+filename)
    if binnumber != len(self.bin_N):
      print(filename+" is misshaped (length is "+np.str(binnumber)+", expected "+np.str(len(self.bin_N))+"), setting to zero")
      if remove_execution_files:      
          try:
              os.remove(os.path.join(filename.rsplit("/",4)[0],"execution",filename.rsplit("/",1)[1]).replace("pTveto.1jet_log_","execution_"))
              print "Removing the corresponding exection file: "+os.path.join(filename.rsplit("/",4)[0],"execution",filename.rsplit("/",1)[1]).replace("pTveto.1jet_log_","execution_")
          except:
              pass
      return

    self.N_tot=int(lines[0])
    for i in range((len(lines)-1)/3):
      self.bin_N[i]=np.int(lines[3*i+1])
      self.bin_values[i]=np.float(lines[3*i+2])
      self.bin_deviations[i]=np.float(lines[3*i+3])

  def printData(self):
    for i in range(len(self.bin_N)):
      print(str(i)+": " + np.str(self.bin_N[i])+", "+np.str(self.bin_values[i])+", "+np.str(self.bin_deviations[i]))

#}}}
#{{{ class: DistributionResult
class DistributionResult:
  def __init__(self,binnumber=0):
    self.binnumber=binnumber
    self.values=np.zeros(binnumber)
    self.sigmas=np.zeros(binnumber)
    self.chi_squared=np.zeros(binnumber)
    
  def printDist(self):
    for i in range(self.binnumber):
      print(str(i)+": "+str(self.values[i])+" +/- "+str(self.sigmas[i]))
#}}}
#{{{ class: DistributionType
class DistributionType:
  def __init__(self, n_scales):
    self.name = ""
    self.binsize = 0.0
    self.startpoint = 0.0
    self.endpoint = 0.0
    self.binnumber = 0
    self.n_scales = n_scales
    self.binningtype = ""
    self.edges = np.empty(0)
    self.binwidths = np.empty(0)
    self.binwidth = 0.0
    
    self.contributions=[]
    for i in range(n_scales):
      self.contributions.append({})
    self.orders=[]
    for i in range(n_scales):
      self.orders.append({})

    self.envelopes_up={}
    self.envelopes_down={}
  
  def readin(self):
    for scale in range(self.n_scales):
      printDebug("scale "+str(scale))
      for c in self.contributions[scale]:
        self.contributions[scale][c]["results"] = {}
        for directory in self.contributions[scale][c]["dirs"]:
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
              if channel_name not in self.contributions[scale][c]["results"]:
                self.contributions[scale][c]["results"][channel_name] = []
              
              data = DistributionData(self.binnumber)
              
              filename_in="../"+directory+"/distribution/CV/scale."+str(scale)+"."+str(scale)+"/"+self.name+"_"+channel_name+".dat"
              #try:
              data.readinData(filename_in)
              self.contributions[scale][c]["results"][channel_name].append(data)
              #except:
                #print("could not readin from "+filename_in)
               # continue
              
              #self.contributions[0][c]["results"][channel_name][0].printData()
            
  def printDist(self):
    print(self.name)
    for scale in range(self.n_scales):
      print("scale "+str(scale))
      for c in self.contributions[scale]:
        print (c)
        for channel in self.contributions[scale][c]["results"]:
          print(channel)
          for data in self.contributions[scale][c]["results"][channel]:
            data.printData()

  def average(self):
    printDebug("averaging..")
    for scale in range(self.n_scales):
      for c in self.contributions[scale]:
        self.contributions[scale][c]["average"]={}
        for channel in self.contributions[scale][c]["results"]:
          averageData=DistributionData(self.binnumber)
          for data in self.contributions[scale][c]["results"][channel]:
# can be removed later
            # first compute, then average
            if data.N_tot==0:
              printDebug(self.name+", "+c+", "+channel+" seems to be empty, setting to zero")
              data.values=np.zeros(self.binnumber)
              data.sigmas=np.zeros(self.binnumber)
            else:
              data.values=data.bin_values/data.N_tot
              data.sigmas=np.sqrt(data.bin_deviations-np.square(data.bin_values)/data.N_tot)/(data.N_tot-1)
            
            printDebug(self.name + ": " + c + ", " +channel+", scale="+np.str(scale))
            
            # first average, then compute
            averageData.N_tot += data.N_tot
            averageData.bin_N += data.bin_N
            averageData.bin_values += data.bin_values
            averageData.bin_deviations += data.bin_deviations
          #print("averaged:")
          #averageData.printData()
          self.contributions[scale][c]["average"][channel] = averageData

  def compute(self):
    for scale in range(self.n_scales):
      for c in self.contributions[scale]:
        self.contributions[scale][c]["final"]={}
        self.contributions[scale][c]["final_cons"]={}
        result_all_channels=DistributionResult(self.binnumber)
        result_all_channels_cons=DistributionResult(self.binnumber)
        for channel in self.contributions[scale][c]["average"]:
          result = DistributionResult(self.binnumber)
          result_cons = DistributionResult(self.binnumber)
          

          # first average subsets of the runs conservatively to gain sufficient statistics. The perform a weighted average afterwards to smoothen the histograms
          data_preaveraged = []
          data_tmp = DistributionData(self.binnumber)
          num_av_datas = 0
          num_to_av = len(self.contributions[scale][c]["results"][channel])/self.contributions[scale][c]["averaging_factor"]

          for data in self.contributions[scale][c]["results"][channel]:
            data_tmp.N_tot += data.N_tot
            data_tmp.bin_N += data.bin_N
            data_tmp.bin_values += data.bin_values
            data_tmp.bin_deviations += data.bin_deviations
            num_av_datas += 1
            if num_av_datas >= num_to_av:
              data_preaveraged.append(data_tmp)
              num_av_datas = 0
              data_tmp = DistributionData(self.binnumber)
          if num_av_datas > 0:
            data_preaveraged.append(data_tmp)

          # now compute cross section for pre-averaged distributions
          for data in data_preaveraged:
            if data.N_tot==0:
              printDebug(self.name+", "+c+", "+channel+" seems to be empty, setting to zero")
              data.values=np.zeros(self.binnumber)
              data.sigmas=np.zeros(self.binnumber)
            else:
              data.values=data.bin_values/data.N_tot
              data.sigmas=np.sqrt((data.bin_deviations-np.square(data.bin_values)/data.N_tot))/(data.N_tot-1)

          # first compute, then average
          #for data in self.contributions[scale][c]["results"][channel]:
          for data in data_preaveraged:
            #print(data.values)
            #print(len(data_preaveraged))
            with np.errstate(divide='ignore', invalid='ignore'):
              result.values += np.where(data.sigmas != 0., data.values / np.square(data.sigmas), 0)

            with np.errstate(divide='ignore', invalid='ignore'):
              result.sigmas += np.where(data.sigmas != 0., 1.0/np.square(data.sigmas), 0)
            
          with np.errstate(divide='ignore', invalid='ignore'):
            result.values = np.where(result.sigmas != 0., result.values / result.sigmas, 0)
          
          with np.errstate(divide='ignore', invalid='ignore'):
            result.sigmas = 1.0/np.sqrt(result.sigmas)
          
          result.sigmas[np.isinf(result.sigmas)] = 0
          
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

          # compute reduced chi squared for each bin
          for data in self.contributions[scale][c]["results"][channel]:
            with np.errstate(divide='ignore', invalid='ignore'):
              result.chi_squared += np.where(data.sigmas != 0., np.square((data.values-result.values))/np.square(data.sigmas), 0)
            
          if len(self.contributions[scale][c]["results"][channel]) > 1:
            result.chi_squared /= (len(self.contributions[scale][c]["results"][channel])-1.0)
          else:
            result.chi_squared = np.zeros(self.binnumber)

          self.contributions[scale][c]["final"][channel] = result

          # first average, then compute
          if self.contributions[scale][c]["average"][channel].N_tot==0:
            #print (self.name+", "+c+", "+channel+" seems to be empty, setting to zero")
            result_cons.values=np.zeros(self.binnumber)
            result_cons.sigmas=np.zeros(self.binnumber)
          else:
            result_cons.values=self.contributions[scale][c]["average"][channel].bin_values/self.contributions[scale][c]["average"][channel].N_tot
            result_cons.sigmas=np.sqrt((self.contributions[scale][c]["average"][channel].bin_deviations-np.square(self.contributions[scale][c]["average"][channel].bin_values)/self.contributions[scale][c]["average"][channel].N_tot))/(self.contributions[scale][c]["average"][channel].N_tot-1)
            
          self.contributions[scale][c]["final_cons"][channel] = result_cons
          #print(self.contributions[scale][c]["average"][channel].N_tot)
          
          result_all_channels.values += result.values
          result_all_channels.sigmas = np.sqrt(np.square(result_all_channels.sigmas)+np.square(result.sigmas))

          result_all_channels_cons.values += result_cons.values
          result_all_channels_cons.sigmas = np.sqrt(np.square(result_all_channels_cons.sigmas)+np.square(result_cons.sigmas))
          self.contributions[scale][c]["finalsum"] = result_all_channels
          self.contributions[scale][c]["finalsum_cons"] = result_all_channels_cons
  
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
    print(self.name+": computing final result")
    for scale in range(self.n_scales):
      for c in self.contributions[scale]:
        order=c.split("_")[0]
        type_correction=c.split("_")[1]
        printDebug(self.name+", "+order)
        if order not in self.orders[scale]:
          self.orders[scale][order]=DistributionResult(self.binnumber)
        
        if self.contributions[scale][c]["averaging_mode"] != "conservative":
          self.orders[scale][order].values += self.contributions[scale][c]["finalsum"].values
          self.orders[scale][order].sigmas = np.sqrt(np.square(self.orders[scale][order].sigmas)+np.square(self.contributions[scale][c]["finalsum"].sigmas))
        else:
          self.orders[scale][order].values += self.contributions[scale][c]["finalsum_cons"].values
          self.orders[scale][order].sigmas = np.sqrt(np.square(self.orders[scale][order].sigmas)+np.square(self.contributions[scale][c]["finalsum_cons"].sigmas))

    # compute envelope of scale variations
    for order in self.orders[0]:
      self.envelopes_up[order] = copy.deepcopy(self.orders[0][order])
      self.envelopes_down[order] =  copy.deepcopy(self.orders[0][order])
      
      for i in range(self.binnumber):
        for scale in range(1,self.n_scales):
          if self.orders[scale][order].values[i] > self.envelopes_up[order].values[i]:
            self.envelopes_up[order].values[i] = self.orders[scale][order].values[i]
            self.envelopes_up[order].sigmas[i] = self.orders[scale][order].sigmas[i]
          if self.orders[scale][order].values[i] < self.envelopes_down[order].values[i]:
            self.envelopes_down[order].values[i] = self.orders[scale][order].values[i]
            self.envelopes_down[order].sigmas[i] = self.orders[scale][order].sigmas[i]

  def saveDist(self,filename,dist,with_chi_squared=False):
    printDebug("saving "+filename)
    f = open(filename,'w')
    for i in range(self.binnumber):
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
      if with_chi_squared:
        f.write(str(self.edges[i])+'\t'+str((dist.values/self.binwidths)[i])+'\t'+str((dist.sigmas/self.binwidths)[i])+'\t'+str(dist.chi_squared[i])+'\n')
      else:
        f.write(str(self.edges[i])+'\t'+str((dist.values/self.binwidths)[i])+'\t'+str((dist.sigmas/self.binwidths)[i])+'\n')


    if with_chi_squared:
      f.write(str(self.endpoint)+'\t'+str((dist.values/self.binwidths)[-1])+'\t'+str((dist.sigmas/self.binwidths)[-1])+'\t'+str(dist.chi_squared[-1])+'\n')
    else:
      f.write(str(self.endpoint)+'\t'+str((dist.values/self.binwidths)[-1])+'\t'+str((dist.sigmas/self.binwidths)[-1])+'\n')

    f.close()

  def save(self,resultdirectory):
    for scale in range(self.n_scales):
      printDebug("")
      printDebug("scale="+np.str(scale))
      for c in self.contributions[scale]:
        contribution=c.split("_")[2]
        type_correction=c.split("_")[1]
        if type_correction=="---":
          type_correction=""
        else:
          type_correction = '.'+type_correction
        printDebug(self.name+", "+c)

        if print_channel_info:
          # channel wise output of reduced chi squared
          for channel in self.contributions[scale][c]["final"]:
            output_dir=resultdirectory+"/scale."+str(scale)+"."+str(scale)+"/"+self.contributions[scale][c]["output_dir"]+"/"+contribution+type_correction+"/channels"
            if not os.path.exists(output_dir):
              os.makedirs(output_dir)

            self.saveDist(output_dir+"/plot."+self.name+"."+contribution+"."+channel+".dat",self.contributions[scale][c]["final"][channel],True)

        # output of contribution wise distributions
        output_dir=resultdirectory+"/scale."+str(scale)+"."+str(scale)+"/"+self.contributions[scale][c]["output_dir"]+"/"+contribution+type_correction
        if not os.path.exists(output_dir):
          os.makedirs(output_dir)

        filename=output_dir+"/plot."+self.name+"."+self.contributions[scale][c]["output_dir"]+"."+contribution+type_correction+".dat"
        self.saveDist(filename,self.contributions[scale][c]["finalsum"])
        
        #filename=output_dir+"/plot."+self.name+"."+self.contributions[scale][c]["output_dir"]+"."+contribution+type_correction+".cons.dat"
        #self.saveDist(filename,self.contributions[scale][c]["finalsum_cons"])
        printDebug("integral: "+np.str(np.sum(self.contributions[scale][c]["finalsum"].values))+" +/- "+np.str(np.sqrt(np.sum(np.square(self.contributions[scale][c]["finalsum"].sigmas)))))
      
      for order in self.orders[scale]:
        if order in ["NLO","NNLO"]:
          filename=resultdirectory+"/scale."+str(scale)+"."+str(scale)+"/plot."+self.name+".."+order+".QCD.dat"
        else:
          filename=resultdirectory+"/scale."+str(scale)+"."+str(scale)+"/plot."+self.name+".."+order+".dat"
        self.saveDist(filename,self.orders[scale][order])
        #print("integral: "+np.str(np.sum(self.orders[scale][order].values)*self.binwidth)+" +/- "+np.str(np.sqrt(np.sum(np.square(self.orders[scale][order].sigmas)))*self.binwidth))
        #f = open(resultdirectory+"_py/scale."+str(scale)+"/plot."+self.name+"."+order+".dat",'w')
        #for i in range(self.binnumber):
          #f.write(str(i*self.binwidth+self.startpoint)+'\t'+str(self.orders[scale][order].values[i])+'\t'+str(self.orders[scale][order].sigmas[i])+'\n')
        #f.write(str(self.endpoint)+'\t'+str(self.orders[scale][order].values[-1])+'\t'+str(self.orders[scale][order].sigmas[-1])+'\n')
        #f.close()
    for order in self.orders[0]:
      printDebug(self.name+", "+order+", envelope")
      self.saveDist(resultdirectory+"/plot."+self.name+"."+order+".up.dat",self.envelopes_up[order])
      self.saveDist(resultdirectory+"/plot."+self.name+"."+order+".down.dat",self.envelopes_down[order])
      
      print(order+" (central, up, down), (from integral):")
      print("integral: "+np.str(np.sum(self.orders[self.n_scales/2][order].values))+" +/- "+np.str(np.sqrt(np.sum(np.square(self.orders[self.n_scales/2][order].sigmas)))))
      print("integral: "+np.str(np.sum(self.envelopes_up[order].values))+" +/- "+np.str(np.sqrt(np.sum(np.square(self.envelopes_up[order].sigmas)))))
      print("integral: "+np.str(np.sum(self.envelopes_down[order].values))+" +/- "+np.str(np.sqrt(np.sum(np.square(self.envelopes_down[order].sigmas))))+"\n")
  
  def printSpecs(self):
    print("name="+self.name)
    print("startpoint="+str(self.startpoint))
    print("endpoint="+str(self.endpoint))
    print("binwidth="+str(self.binwidth))
    print("binnumber="+str(self.binnumber))

  def cleanup(self):
    del(self.contributions)

#}}}
def completeNewDistribution(new_distribution):
# compute the non-specified properties
  print(new_distribution.name)
  if new_distribution.binningtype=="linear":
    if new_distribution.binnumber == 0:
      new_distribution.binnumber = np.int((new_distribution.endpoint - new_distribution.startpoint)/new_distribution.binwidth)
    elif new_distribution.binwidth == 0:
      new_distribution.binwidth = (new_distribution.endpoint - new_distribution.startpoint)/new_distribution.binnumber
    elif new_distribution.endpoint == 0:
      new_distribution.endpoint = new_distribution.startpoint + new_distribution.binnumber*new_distribution.binwidth
    elif new_distribution.startpoint == 0:
      new_distribution.startpoint = new_distribution.endpoint - new_distribution.binnumber*new_distribution.binwidth

    new_distribution.binwidths = np.empty(new_distribution.binnumber)
    new_distribution.binwidths.fill(new_distribution.binwidth)
    new_distribution.edges=np.arange(new_distribution.startpoint,new_distribution.endpoint,new_distribution.binwidth)
  elif new_distribution.binningtype=="logarithmic":
    new_distribution.binwidths = np.empty(new_distribution.binnumber)
    new_distribution.edges=np.empty(new_distribution.binnumber+1)
    step = np.log10(new_distribution.endpoint/new_distribution.startpoint)/new_distribution.binnumber
    for i in range(new_distribution.binnumber):
      new_distribution.edges[i] = np.power(10,np.log10(new_distribution.startpoint)+i*step)
      new_distribution.binwidths[i] = np.power(10,np.log10(new_distribution.startpoint)+(i+1)*step)-np.power(10,np.log10(new_distribution.startpoint)+i*step)

    new_distribution.edges[-1]=new_distribution.endpoint
  elif new_distribution.binningtype=="irregular":
    new_distribution.binnumber=len(new_distribution.edges)-1
    new_distribution.startpoint = new_distribution.edges[0]
    new_distribution.endpoint = new_distribution.edges[-1]
    new_distribution.binwidths = np.empty(new_distribution.binnumber)
    for i in range(new_distribution.binnumber):
      new_distribution.binwidths[i] = new_distribution.edges[i+1]-new_distribution.edges[i]
  else:
    print("binngtype "+new_distribution.binningtype+" not supported, aborting")
    assert(False)

  return new_distribution

##{{{ def: parseFileDistributions(filename, n_scales, distributions)
def parseFileDistributions(filename, n_scales, distributions):
  new_distribution = DistributionType(n_scales)
  f = open(filename,'r')
    #line = f.read()
  for line in f:
    if len(line)==1:
      continue
    if line[0]=='#':
      continue;
    if len(line.split()) < 2 or line.split()[1] != '=':
      print("could not parse \""+line+"\" in "+filename+", ignoring line")
      continue

    print(line)

    variable = line.split()[0]
    value = line.split()[2]
    printDebug(variable+" = " + value)
    if variable=="distributionname":
      printDebug("new distribution: "+value)
      #append old distribution, if any
      if new_distribution.name!="":
        new_distribution = completeNewDistribution(new_distribution)
      	distributions[new_distribution.name] = new_distribution
#	distributions[new_distribution.name].printSpecs()
	new_distribution = DistributionType(n_scales)
      new_distribution.name=value
      new_distribution.binningtype="linear" # default binning type
      
    elif variable=="startpoint":
      new_distribution.startpoint=float(value)
    elif variable=="endpoint":
      new_distribution.endpoint=float(value)
    elif variable=="binwidth":
      new_distribution.binwidth=float(value)
    elif variable=="binnumber":
      new_distribution.binnumber=int(value)
    elif variable=="binningtype":
      new_distribution.binningtype=value
    elif variable=="edges":
      edges = value.split(':')
      for edge in edges:
        new_distribution.edges = np.append(new_distribution.edges,np.float(edge))
      
  if new_distribution.name!="":
    new_distribution = completeNewDistribution(new_distribution)
    distributions[new_distribution.name] = new_distribution

  f.close()
#}}}
#{{{ def: parseFileParameter(filename)
def parseFileParameter(filename):
  n_scales_CV = 0

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

  f.close()
  
  return n_scales_CV
#}}}
#{{{ def: parseLine(line)
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
#}}}
#{{{ def: parseContributionFile(filename, contributions)
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
#}}}
#{{{ def: parseInfileDistribution(filename, contributions)
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
#}}}
#{{{ def: combine_distributions(filename)
def combine_distributions(directory,infile):
  os.chdir(directory)
  n_scales_CV = parseFileParameter('../file_parameter.dat')

  distributions = {}
  parseFileDistributions('../file_distribution.dat', n_scales_CV, distributions)
  print("readin file_distributions.dat done")
  if print_debug:
    for d in distributions:
      distributions[d].printSpecs()
    
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

  for scale in range(n_scales_CV):
    if not os.path.exists(resultDirectory+"_py/CV/scale."+str(scale)+"."+str(scale)):
      os.makedirs(resultDirectory+"_py/CV/scale."+str(scale)+"."+str(scale))
  
  
  for c in contributions:
    printDebug(c+": ")
    printDebug(contributions[c])
    
  for d in distributions:
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
        distributions[d].contributions[scale][c] = new_c
        printDebug(distributions[d].contributions[scale][c])
    distributions[d].readin()
    #distributions[d].printDist()
    distributions[d].average()
    distributions[d].compute()
    distributions[d].computeOrders()
    #distributions[d].printResult()
    print("writing results to "+resultDirectory+"_py/CV")
    distributions[d].save(resultDirectory+"_py/CV")
    distributions[d].cleanup()
#}}}

if __name__=="__main__":
    infile = sys.argv[1]
    current_dir = os.getcwd()
    combine_distributions(current_dir,infile)
