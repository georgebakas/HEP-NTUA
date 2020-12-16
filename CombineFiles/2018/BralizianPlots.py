import sys
import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend
import CMS_lumi, tdrstyle
import subprocess # to execute shell command
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# CMS style
CMS_lumi.cmsText = "CMS"
CMS_lumi.extraText = "Preliminary"
CMS_lumi.cmsTextSize = 0.65
CMS_lumi.outOfFrame = True
tdrstyle.setTDRStyle()


# CREATE datacards
def createDataCardsThetaB(masses):

    datacard_lines = [ "imax * number of bins",
                       "jmax * number of processes minus 1",
                       "kmax * number of nuisance parameters",
                       "----------------------------------------------------------------------------",
                       "shapes *         SR_C       ProcessesFile_1500.root h_chi_$PROCESS",
                       "shapes Zprime    SR_C       ",
                       "shapes data_obs  SR_C       DataFile_1500.root h_Data",
                       "----------------------------------------------------------------------------",
                       "bin          SR_C",
                       "observation  -1.0",
                       "----------------------------------------------------------------------------",
                       "bin                               SR_C        SR_C        SR_C         SR_C",
                       "process                           Zprime      qcd         Subdominant  ttbar",
                       "process                           0           1           2            3",
                       "rate                              -1.0        -1.0        -1.0         -1.0",
                       "----------------------------------------------------------------------------",
                       "yield_ttbar           lnN         -           -           -            1.5",
                       "yield_qcd              lnN        -           1.5         -            -",
                       "yield_Subdominant      lnN        -           -          1.5            -",
                       "* autoMCStats 10 0 1",
                      ]

    # make datacards for differents values of theta_B
    for mass in masses:
        datacard = open("datacard_chi_SR_"+mass+".txt", 'w')
        for iline, line in enumerate(datacard_lines):
            if(iline == 5):
                 rootfile = "ZprimeFile"+mass+"_massCut1500.root"
                 rootfile = rootfile.replace("ZprimeFilemZ", "ZprimeFile")
                 print rootfile
                 datacard.write(line+"\t"+rootfile+" h_chi_$PROCESS \n")
            else:
                datacard.write(line+"\n")
        datacard.close()
        print "datacard_chi_SR_"+mass+".txt"

    return 1


# EXECUTE datacards
def executeDataCards(masses):

    for mass in masses:
        file_name = "datacard_chi_SR_"+mass+".txt"
        combine_command = "combine -M AsymptoticLimits -d "+file_name+" -n ."+mass
        print ""
        print ">>> " + combine_command
        p = subprocess.Popen(combine_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            print line.rstrip("\n")
        retval = p.wait()


# GET limits from root file
def getLimits(file_name):
    print file_name
    file = TFile(file_name)
    tree = file.Get("limit")

    limits = [ ]
    for quantile in tree:
        limits.append(tree.limit)
        print ">>>   %.2f" % limits[-1]

    return limits


# PLOT upper limits
def plotUpperLimits(file_names, values, width):
    # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/

    if width == 1:
        #here you have 2000_20, 2500_25, 3000_30, 3500_35, 4000_40
        cross_section = [0.165, 0.04715, 0.01495, 0.005108, 0.001908]
    elif width == 10:
        #here we got 2000_200, 2500_250, 3000_300, 3500_350, 4000_400
        cross_section = [0.01825, 0.005708, 0.002056, 0.0008352, 0.0003779]
    else:
        #30% we got 2000_600, 2500_750, 3000_900, 3500_1050, 4000_1200
        cross_section = [0.006399, 0.002255, 0.0009167, 0.0004248, 0.0002186]

    N = len(file_names)
    yellow = TGraph(2*N)    # yellow band
    green = TGraph(2*N)     # green band
    median = TGraph(N)      # median line
    observed = TGraph(N)

    up2s = [ ]
    for i in range(len(file_names)):
        file_name = "higgsCombine."+file_names[i]+".AsymptoticLimits.mH120.root"
        limit = getLimits(file_name)
        limit = [x * cross_section[i] for x in limit]
        up2s.append(limit[4])
        observed.SetPoint(  i,    values[i], limit[5] ) # observed
        yellow.SetPoint(    i,    values[i], limit[4] ) # + 2 sigma
        green.SetPoint(     i,    values[i], limit[3] ) # + 1 sigma
        median.SetPoint(    i,    values[i], limit[2] ) # median
        green.SetPoint(  2*N-1-i, values[i], limit[1] ) # - 1 sigma
        yellow.SetPoint( 2*N-1-i, values[i], limit[0] ) # - 2 sigma

    W = 800
    H  = 600
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    R = 0.04*W
    c = TCanvas("c","c",100,100,W,H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetGrid()
    c.SetLogy(True)
    c.cd()
    frame = c.DrawFrame(1.4,0.0001, 4.1, 10)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("95% upper limit on #sigma #times BR")
    #frame.GetYaxis().SetTitle("95% upper limit on #sigma #times BR / (#sigma #times BR)_{SM}")
    #frame.GetXaxis().SetTitle("background systematic uncertainty [%]")
    frame.GetXaxis().SetTitle("M_{Z} (GeV)")
    frame.SetMinimum(0.01)
    frame.SetMaximum(max(up2s)*10.5)
    frame.GetXaxis().SetLimits(min(values),max(values))

    yellow.SetFillColor(ROOT.kOrange)
    yellow.SetLineColor(ROOT.kOrange)
    yellow.SetFillStyle(1001)
    yellow.Draw('F')

    green.SetFillColor(ROOT.kGreen+1)
    green.SetLineColor(ROOT.kGreen+1)
    green.SetFillStyle(1001)
    green.Draw('Fsame')

    median.SetLineColor(1)
    median.SetLineWidth(2)
    median.SetLineStyle(2)
    median.Draw('Lsame')

    observed.SetLineColor(1)
    observed.SetLineWidth(2)
    observed.SetLineStyle(1)
    observed.SetMarkerStyle(20)
    observed.SetMarkerColor(1)
    observed.Draw('Lsame')

    CMS_lumi.CMS_lumi(c,14,11)
    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')

    x1 = 0.15
    x2 = x1 + 0.24
    y2 = 0.76
    y1 = 0.60
    legend = TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)
    legend.AddEntry(median, "Asymptotic CL_{s} Expected",'L')
    legend.AddEntry(observed, "Asymptotic CL_{s} Observed",'L')
    legend.AddEntry(green, "#pm 1 std. deviation",'f')
    legend.AddEntry(yellow,"#pm 2 std. deviation",'f')
#    legend.AddEntry(green, "Asymptotic CL_{s} #pm 2 std. deviation",'f')
    legend.Draw()

    print 'UpperLimit_w%s.png' % (width)
    label = 'UpperLimit_w%s.png' % (width)
    c.SaveAs(label)
    c.Close()


# RANGE of floats
def frange(start, stop, step):
    i = start
    while i <= stop:
        yield i
        i += step


# MAIN
def main(width):

    #labels = [ ]
    #values = [ ]

    if width == 1:
        file_names = ["mZ_2000_20", "mZ_2500_25", "mZ_3000_30", "mZ_3500_35", "mZ_4000_40"]
    elif width ==10:
        file_names = ["mZ_2000_200", "mZ_2500_250", "mZ_3000_300", "mZ_3500_350", "mZ_4000_400"]
    else:
        file_names = ["mZ_2000_600", "mZ_2500_750", "mZ_3000_900", "mZ_3500_1050", "mZ_4000_1200"]

    values = [2000, 2500, 3000, 3500, 4000];

    print(file_names)
    #createDataCardsThetaB(file_names)
    executeDataCards(file_names)
    plotUpperLimits(file_names, values, width)



if __name__ == '__main__':
   width = int(sys.argv[1])
   print type(width)
   main(width)
