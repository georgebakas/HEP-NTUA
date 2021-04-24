import sys
import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend
import CMS_lumi, tdrstyle
import subprocess # to execute shell command
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# CMS style
CMS_lumi.cmsText = "CMS"
CMS_lumi.extraText = "Work in Progress"
CMS_lumi.cmsTextSize = 0.65
CMS_lumi.outOfFrame = True
tdrstyle.setTDRStyle()


# CREATE datacards
def createDataCardsThetaB(masses, mass_cut):

    datacard_lines = [ "imax * number of bins",
                       "jmax * number of processes minus 1",
                       "kmax * number of nuisance parameters",
                       "----------------------------------------------------------------------------",
                       "shapes *         SR_C       ProcessesFile_"+mass_cut+".root h_chi_$PROCESS h_chi_$PROCESS_$SYSTEMATIC",
                       "shapes Zprime    SR_C       ",
                       "shapes data_obs  SR_C       DataFile_"+mass_cut+".root h_Data",
                       "----------------------------------------------------------------------------",
                       "bin          SR_C",
                       "observation  -1.0",
                       "----------------------------------------------------------------------------",
                       "bin                               SR_C        SR_C        SR_C         SR_C",
                       "process                           Zprime      qcd         Subdominant  ttbar",
                       "process                           0           1           2            3",
                       "rate                              -1.0        -1.0        -1.0         -1.0",
                       "----------------------------------------------------------------------------",
                       "yield_ttbar            lnN         -           -           -            1.5",
                       "yield_qcd              lnN        -           1.5         -            -",
                       "yield_Subdominant      lnN        -           -          1.5            -",
                       "lumi_13TeV             lnN         -         1.025      1.025        1.025",
                       "scale                  shapeN       -            -           -           1.0",
                       "pdf                    shapeN       -            -           -           1.0",
                       "fsr                    shapeN       -            -           -           1.0",
                       "isr                    shapeN       -            -           -           1.0",
                       "btag                   shapeN       -            -           -           1.0",
                       "*  autoMCStats  10  1"
                      ]

    # make datacards for differents values of theta_B
    for mass in masses:
        datacard = open("datacard_chi_SR_"+mass+"_cut_"+mass_cut+".txt", 'w')
        for iline, line in enumerate(datacard_lines):
            if(iline == 5):
                 rootfile = "ZprimeFile"+mass+"_massCut"+mass_cut+".root"
                 rootfile = rootfile.replace("ZprimeFilemZ", "ZprimeFile")
                 print rootfile
                 datacard.write(line+"\t"+rootfile+" h_chi_$PROCESS \n")
            else:
                datacard.write(line+"\n")
        datacard.close()
        print "datacard_chi_SR_"+mass+"_cut_"+mass_cut+".txt"

    return 1

# EXECUT FitDiagnostics
def executFitDiagnostics(masses):

    for mass in masses:
        file_name = "datacard_chi_SR_"+mass+".txt"
        fit_command = 'combine -M FitDiagnostics -d '+file_name+' -n _fit_result --saveShapes --saveWithUncertainties'
        print ""
        print ">>> " + fit_command
        p = subprocess.Popen(fit_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            print line.rstrip("\n")
        retval = p.wait()

# EXECUTE datacards
def executeDataCards(masses, mass_cut):

    for mass in masses:
        file_name = "datacard_chi_SR_"+mass+"_cut_"+mass_cut+".txt"
        combine_command = "combine -M AsymptoticLimits -d "+file_name+" -n ."+mass+"Cut"+mass_cut
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
def plotUpperLimits(file_names, values, width, mass_cut):
    # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/

    if width == 1:
        #here you have 1200_12, 1400_14, 1600_16, 1800_18, 2000_20,
        #              2500_25, 3000_30, 3500_35, 4000_40, 4500_45
        cross_section = [1.736e+00, 9.096e-01, 5.742e-01, 2.830e-01, 1.662e-01,
                        4.749e-02, 1.494e-02, 5.105e-03, 1.900e-03, 7.613e-04]
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
    theoretical_ = TGraph(N)

    up2s = [ ]
    for i in range(len(file_names)):
        print "HERE:!!", file_names[i], mass_cut[file_names[i]]
        file_name = "higgsCombine."+file_names[i]+"Cut"+str(mass_cut[file_names[i]])+".AsymptoticLimits.mH120.root"
        limit = getLimits(file_name)
        limit = [x * cross_section[i] for x in limit]
        print 'cross_section:', cross_section[i]
        for ilim in range(6):
            print i, values[i], limit[ilim]

        up2s.append(limit[4])
        observed.SetPoint(  i,    values[i], limit[5] ) # observed
        yellow.SetPoint(    i,    values[i], limit[4] ) # + 2 sigma
        green.SetPoint(     i,    values[i], limit[3] ) # + 1 sigma
        median.SetPoint(    i,    values[i], limit[2] ) # median
        green.SetPoint(  2*N-1-i, values[i], limit[1] ) # - 1 sigma
        yellow.SetPoint( 2*N-1-i, values[i], limit[0] ) # - 2 sigma

        theoretical_.SetPoint(i, values[i], cross_section[i])


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
    #frame.SetMinimum(0.05)
    #frame.SetMaximum(max(up2s)*10.5)
    #frame.SetMaximum(5)
    #frame.SetMinimum(0.01)
    #frame.SetMaximum(max(up2s)*50.5)
    frame.SetMinimum(0.001)
    frame.SetMaximum(50)
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

    theoretical_.SetLineColor(6)
    theoretical_.SetLineWidth(3)
    theoretical_.SetLineStyle(5)
    theoretical_.SetMarkerStyle(20)
    theoretical_.SetMarkerColor(6)
    theoretical_.Draw('Lsame')

    CMS_lumi.CMS_lumi(c,14,11)
    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')

    x1 = 0.40
    x2 = x1 + 0.24
    y1 = 0.70
    y2 = y1 + 0.16
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

    label = 'UpperLimit_w%s_slidingMassCut_test.png' % (width)
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
        file_names = ["mZ_1200_12", "mZ_1400_14", "mZ_1600_16", "mZ_1800_18", "mZ_2000_20",
                    "mZ_2500_25", "mZ_3000_30", "mZ_3500_35", "mZ_4000_40", "mZ_4500_45"]
    elif width ==10:
        file_names = ["mZ_2000_200", "mZ_2500_250", "mZ_3000_300", "mZ_3500_350", "mZ_4000_400"]
    else:
        file_names = ["mZ_2000_600", "mZ_2500_750", "mZ_3000_900", "mZ_3500_1050", "mZ_4000_1200"]

    values = [1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500]

    # I need to map these with the correct mass Cuts
    mass_cut = {"mZ_1200_12":1000, "mZ_1400_14":1200, "mZ_1600_16":1400, "mZ_1800_18":1600, "mZ_2000_20":1600,
                    "mZ_2500_25":2000, "mZ_3000_30":2000, "mZ_3500_35":2000, "mZ_4000_40":2000, "mZ_4500_45":2000}

    print(file_names)
    #createDataCardsThetaB(file_names, mass_cut)
    #executFitDiagnostics(file_names)
    #executeDataCards(file_names, mass_cut)
    plotUpperLimits(file_names, values, width, mass_cut)



def execute_all_masses(width):

    if width == 1:
        file_names = ["mZ_1200_12", "mZ_1400_14", "mZ_1600_16", "mZ_1800_18", "mZ_2000_20",
                    "mZ_2500_25", "mZ_3000_30", "mZ_3500_35", "mZ_4000_40", "mZ_4500_45"]
    elif width == 10:
        file_names=["mZ_2000_200", "mZ_2500_250", "mZ_3000_300", "mZ_3500_350", "mZ_4000_400"]
    else:
        file_names = ["mZ_2000_600", "mZ_3000_900", "mZ_3500_1050", "mZ_4000_1200"]

    values = [2000, 2500, 3000, 3500, 4000]

    mJJ_cuts = ['1000', '1200', '1400', '1600', '1800']

    for mJJCut in mJJ_cuts:
        executeDataCards(file_names, mJJCut)


if __name__ == '__main__':
   width = int(sys.argv[1])
   print type(width)
   main(width)
   #execute_all_masses(width)
