#!/usr/bin/python

#import htcondor
#import classad
import os
import sys


def CreateCondorSubmitFile(executable, arguments, inputFiles, outputFile, flavour='testmatch'):
    jobName = outputFile.replace('.root', '')
    f = open(jobName + ".sub", "w")
    f.write("universe = vanilla\n")
    f.write("executable = " + executable + "\n")
    f.write("arguments = " + arguments + "\n")
    f.write("error = " + jobName + ".err\n")
    f.write("output = " + jobName + ".out\n")
    f.write("log = " + jobName + ".log\n")
    f.write("+JobFlavour = \""+flavour+"\"\n")
    f.write("+AccountingGroup = \"group_u_CMS.u_zh.users\"\n")
    f.write("\n")
    f.write("should_transfer_files = yes\n")
    f.write("transfer_input_files = " + inputFiles + "\n")
    f.write("when_to_transfer_output = ON_EXIT\n")
    f.write("transfer_output_files = " + outputFile + "\n")
    f.write("queue\n")

    f.close()


def submitCondorJobs(executable, arguments, inputFiles, outputFile):
    CreateCondorSubmitFile(executable, arguments, inputFiles, outputFile)
    os.system("condor_submit -name bigbird09.cern.ch " + outputFile.replace('.root', '') + ".sub")
