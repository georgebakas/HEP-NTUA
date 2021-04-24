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


# PLOT upper limits
def plotUpperLimits(values, width):
    # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/


    N = len(values)
    print N
    yellow = TGraph(2*N)    # yellow band
    green = TGraph(2*N)     # green band
    median = TGraph(N)      # median line
    observed = TGraph(N)

    up2s = [ ]
    for i, ival in enumerate(sorted(values)):
        print i, ival, values[ival][5]
        up2s.append(values[ival][4])
        observed.SetPoint(  i,    ival,     values[ival][5] ) # observed
        yellow.SetPoint(    i,    ival,     values[ival][4] ) # + 2 sigma
        green.SetPoint(     i,    ival,     values[ival][3] ) # + 1 sigma
        median.SetPoint(    i,    ival,     values[ival][2] ) # median
        green.SetPoint(  2*N-1-i, ival,     values[ival][1] ) # - 1 sigma
        yellow.SetPoint( 2*N-1-i, ival,     values[ival][0] ) # - 2 sigma

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
    frame.SetMaximum(max(up2s)*50.5)
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

    values = {}
    #-2 sigma [0],  -1 sigma, median, +1 sigma, +2 sigma [4],  observed[5]
    if width == 1:
        values = {
            1000:[0.746, 1.04, 1.47, 2.15, 3.01, 1.8],
            1250:[0.264, 0.377, 0.534, 0.778, 1.16, 1.14],
            1500:[0.145, 0.202, 0.291, 0.425, 0.617, 0.239],
            2000:[0.0568, 0.08, 0.117, 0.17, 0.235, 0.104],
            2500:[0.0314, 0.0443, 0.0614, 0.09, 0.132, 0.0464],
            3000:[0.0244, 0.033, 0.0469, 0.0708, 0.0992, 0.0462],
            3500:[0.0192, 0.0257, 0.036, 0.0554, 0.0813, 0.0248],
            4000:[0.0163, 0.022, 0.0318, 0.0488, 0.0749, 0.0224]
        }
    elif width == 10:
        values = {
            2000:[0.378, 0.252, 0.176, 0.12, 0.0896, 0.149],
            2500:[0.2, 0.145, 0.0965, 0.0663, 0.0478, 0.0684],
            3000:[0.156, 0.109, 0.0739, 0.0527, 0.0394, 0.0667],
            3500:[0.156, 0.101, 0.0686, 0.0467, 0.0341, 0.0507],
            4000:[0.171, 0.108, 0.0718, 0.0476, 0.0351, 0.0495]
        }
    else:
        values = {
            2000:[0.611, 0.409, 0.28, 0.191, 0.133, 0.216],
            3000:[0.318, 0.22, 0.151, 0.106, 0.079, 0.116],
            4000:[0.335, 0.228, 0.149, 0.102, 0.0754, 0.104]
        }

    plotUpperLimits(values, width)



if __name__ == '__main__':
   width = int(sys.argv[1])
   print type(width)
   main(width)
