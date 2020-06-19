/*
    This is a piece of code dedicated for Signal Extraction:
    For every variable x_reco we need to extract the signal using the formula:

    S(x) = D(x) - (Ryield * Ryield_correction * Nbkg) * QCDShapeCorrection * Q(x) - B(x)
    where:
    1.S(x) is the extracted signal
    2.Ryield is the NSR / NSRA taken from the file--> year/TransferFactor.root
    3.Nbkg is the number of QCD events in the SR taken from the simultaneous fit leaving the eb free
    4.QCDShapeCorrection is taken from the corrected shape (only for 2017 and 2018 for specific variables) from files: FitOutput.root
    5.Q(x) is taken from the CR from data
    6.B(x) is the signal region of the Subdominant bkg (exp yield) taken from MC
*/

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TGraphErrors.h"

using namespace std;

using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstants.h"

TH1F *getRebinned(TH1F *h, float BND[], int N);
void SignalExtractionSpecific(TString year = "2016", TString variable = "jetPt0", TString fitRecoVar = "leadingJetPt");

TH1F *getRebinned(TH1F *h, float BND[], int N)
{
  TString name = TString(h->GetName())+"_Rebinned";
  TH1F *hNew = new TH1F(name, name, N, BND);

  for(int i=0;i<h->GetNbinsX();i++) {
    float x = h->GetBinCenter(i+1);
    float y = h->GetBinContent(i+1);
    float e = h->GetBinError(i+1);
    for(int j=0;j<hNew->GetNbinsX();j++) {
      float x1 = hNew->GetBinLowEdge(j+1);
      float x2 = x1+hNew->GetBinWidth(j+1);
      if ((x>x1) && (x<x2)) {
        float yNew = hNew->GetBinContent(j+1);
        float eNew = hNew->GetBinError(j+1);
        hNew->SetBinContent(j+1,yNew+y);
        hNew->SetBinError(j+1,sqrt(e*e+eNew*eNew));
        break;
      }
    }
  }
  return hNew;

}

void SignalExtraction_UnequalBins(TString year)
{
    TString vars[] = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1", "jetY0", "jetY1"};
    TString fitRecoVar[] = {"mJJ", "ptJJ", "yJJ", "leadingJetPt","subleadingJetPt", "leadingJetY", "subleadingJetY"};
    for(int i=0; i<sizeof(vars)/sizeof(vars[0]); i++)
    {
        SignalExtractionSpecific(year, vars[i], fitRecoVar[i]);
        //break;
    }
}


void SignalExtractionSpecific(TString year = "2016", TString variable = "jetPt0", TString fitRecoVar = "leadingJetPt")
{
    initFilesMapping();
    cout<<luminosity[year]<<endl;

    gStyle->SetOptStat(0);
    //open the signal file: get D(x) and Q(x) for every variable
    TFile *infDataMedium = TFile::Open(TString::Format("%s/Histo_Data_%s_100_reduced_UnequalBinning.root", year.Data(), year.Data()));
    //cout<<TString::Format("%s/Histo_Data_%s_100_reduced.root", year.Data(), year.Data())<<endl;
    TH1F *hD = (TH1F*)infDataMedium->Get(TString::Format("hWt_%s_2btag", variable.Data()));
    TH1F *hQ = (TH1F*)infDataMedium->Get(TString::Format("hWt_%s_0btag_expYield", variable.Data()));
    //hQ has to be scaled to integral, because we need the shape

    //open the file to get the Ryield
    TFile *infRyield = TFile::Open(TString::Format("%s/Ryield/TransferFactor_HT300toInf_100_%s.root",year.Data(),variable.Data()));
    //TH1F *hRyield = (TH1F*)infRyield->Get("ClosureTest_TransferFactor");
    TH1F *hRyield = (TH1F*)infRyield->Get("dataTransferFactor");

    float Ryield = hRyield->GetBinContent(1);
    float Ryield_error = hRyield->GetBinError(1);

    float r_yield_correction;
    TH1F *hRyieldMC = (TH1F*)infRyield->Get("ClosureTest_TransferFactor");
    r_yield_correction = (hRyieldMC->GetBinContent(2)/ hRyieldMC->GetBinContent(1));
    //float r_yield_correction_error = TMath::Sqrt(TMath::Power(hRyieldMC->GetBinError(2)/hRyieldMC->GetBinContent(1),2)
    //                                  + TMath::Power((hRyieldMC->GetBinContent(2)*hRyieldMC->GetBinError(2))/TMath::Power(hRyieldMC->GetBinContent(1),2),2));

    float r_yield_correction_error = 0.104911;
    cout<<"-------------------------"<<endl;
    cout<<"Ryield_data (0): "<<Ryield<<" ± "<<Ryield_error<<endl;
    cout<<"r_yield_correction: "<<r_yield_correction<<" ± "<<r_yield_correction_error<<endl;

    float corrected_rYield = r_yield_correction * Ryield;
    float corrected_error = TMath::Sqrt(TMath::Power(r_yield_correction*Ryield_error,2) + TMath::Power(r_yield_correction_error*Ryield,2));
    //cout<<"Ryield_data (2): "<<hRyield->GetBinContent(2)<<endl;
    cout<<"corrected Ryield: "<<corrected_rYield<<" ± "<<corrected_error<<endl;
    //return;
    //open the file to get the Nbkg
    float NQCD = Nbkg2Constants[TString::Format("Nbkg%s",year.Data())];
    float NQCD_error = Nbkg2ConstantsErrors[TString::Format("Nbkg%s_error",year.Data())];

    //Subdominant bkgs files
    TFile *infSub = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100_reduced_UnequalBinning.root", year.Data()));
    TH1F *hSub = (TH1F*)infSub->Get(TString::Format("hWt_%s_2btag_expYield", variable.Data()));
    TH1F *hSub_0 = (TH1F*)infSub->Get(TString::Format("hWt_%s_0btag_expYield", variable.Data()));
    //here I will import correction factors for QCD if needed...

    TFile *infSignalMC;
    infSignalMC = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_100_reduced_UnequalBinning.root", year.Data()));
    //else infSignalMC = TFile::Open(TString::Format("%s/Histo_TT_Mtt-700toInf_100_reduced_UnequalBinning.root", year.Data()));
    TH1F *hSMC = (TH1F*)infSignalMC->Get(TString::Format("hWt_%s_2btag_expYield", variable.Data()));
    TH1F *hSMC_0= (TH1F*)infSignalMC->Get(TString::Format("hWt_%s_0btag_expYield", variable.Data()));


    //now i need the rebin function for aaaaaall my hists so that I am conistent
    //I will include binning in the TemplateConstants.h

    std::vector< std::vector <Float_t> > const BND = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3700, 5000}, //mjj
  													                          {0,60,150,300,450,600,750,900,1300}, //ptjj
  													                          {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
  		   	                                            {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}, //jetPt0 21
  													                          {400,450,500,570,650,750,850,1000,1200,1500}, //jetPt1
  													                          {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
                                                      {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1
                            

    int nBins[BND.size()];
    for (int i = 0; i<BND.size(); i++)
    	nBins[i] = BND[i].size()-1;

    //get the respected integer from the mapping in the TemplateConstants.h
    int selVar = variableConstant[variable];

    //use this template for the initialization of the BND source given as input for the rebinned histos
    float tempBND[nBins[selVar]+1];
    std::copy(BND[selVar].begin(), BND[selVar].end(), tempBND);

    //rebin all histograms with the same method
    TH1F *hD_rebinned, *hQ_rebinned, *hSub_rebinned;
    TH1F *hQ_copy = (TH1F*)hQ->Clone(TString::Format("hD_copy_%s", variable.Data()));
    cout<<hQ_copy->GetNbinsX()<<endl;
    //before remove subdominat and ttbar contribution from Data 0btag (IS IT NEEDED???)
    //hQ->Add(hSub_0, -1);
    //hQ->Add(hSMC_0, -1);

    /*
    cout<<"hData 0:"<<hQ->Integral()<<endl;
    cout<<"hD_copy 0:"<<hQ_copy->Integral()<<endl;
    cout<<"hSub 0:"<<hSub_0->Integral()<<endl;
    cout<<"httbar 0:"<<hSMC_0->Integral()<<endl;

    TCanvas *c1 = new TCanvas(TString::Format("c_copy_%s", variable.Data()), TString::Format("c_copy_%s", variable.Data()), 800, 600);
    hQ->Scale(1./hQ->Integral());
    hQ_copy->Scale(1./hQ_copy->Integral());

    hQ->Draw("hist");
    hQ->SetLineColor(kRed);
    hQ_copy->SetLineColor(kBlue);
    hQ_copy->Draw("hist same");
    */

    /* DO NOT REBIN!!!!!
    hD_rebinned = getRebinned(hD, tempBND, nBins[selVar]);
    hQ_rebinned = getRebinned(hQ, tempBND, nBins[selVar]);
    hSub_rebinned = getRebinned(hSub, tempBND, nBins[selVar]);

    */
    hD_rebinned = (TH1F*)hD->Clone();
    hQ_rebinned = (TH1F*)hQ->Clone();
    hSub_rebinned = (TH1F*)hSub->Clone();

    TH1F *hSignal = (TH1F*)hD_rebinned->Clone(TString::Format("hSignal_%s",variable.Data()));
    //work on the elements for the QCD so that we have the right Q(x)
    //Ryield * Nbkg  * Q(x)
    cout<<variable<<endl;

    float SF[hQ_rebinned->GetNbinsX()];
    //QCD correction factor for shape
    cout<<variable.Data()<<endl;
    if(variable.EqualTo("jetPt0") || (variable.EqualTo("jetPt1") && !year.EqualTo("2016")))
    {
        TFile *fitFile =  TFile::Open(TString::Format("../QCD_ClosureTests_All/closureTest_fitResults_%s_reduced.root",year.Data()));
        //TF1 *fitResult = (TF1*)fitFile->Get(TString::Format("func_%s",variable.Data()));
        TF1 *fitResult = (TF1*)fitFile->Get(TString::Format("FitFunction_%s",fitRecoVar.Data()));
        for(int i=0; i<hQ_rebinned->GetNbinsX(); i++)
        {
            float chi = hQ_rebinned->GetBinCenter(i+1);
            SF[i] = fitResult->Eval(chi);
        }
    }
    else
     for(int i=0; i<hQ_rebinned->GetNbinsX(); i++) SF[i] = 1;

    hQ_rebinned->Scale(1./hQ_rebinned->Integral());  //this is how you get the shape
    cout<<"--------"<<endl;
    cout<<"NQCD: "<<NQCD<<endl;
    cout<<"---"<<endl;
    for(int i =0; i<hQ_rebinned->GetNbinsX(); i++)
    {
        float oldContent = hQ_rebinned->GetBinContent(i+1);
        float oldError = hQ_rebinned->GetBinError(i+1);
        float newContent;
        //cout<<"i: "<<i+1<<" old content* SF: "<<oldContent* SF[i]<<endl; //* SF[i]<<endl;
        //cout<<oldContent * SF[i] * corrected_rYield * NQCD<<" ± "<<endl;
        newContent = oldContent * corrected_rYield * NQCD * SF[i];
        //cout<<Ryield * r_yield_correction * NQCD * oldContent * SF[i]<<endl;
        //cout<<NQCD2_reduced[year.Data()] * oldContent *SF[i]<<endl;
        //cout<<"i: "<<i+1<<", with content: "<<newContent<<endl;
        float newError   = TMath::Sqrt(TMath::Power(oldError*NQCD*corrected_rYield,2)+
                                       TMath::Power(NQCD_error*oldContent*SF[i]*corrected_rYield,2) +
                                       TMath::Power(corrected_error*NQCD*oldContent*SF[i],2));

        /*cout<<"bin: "<<i+1<<endl;
        cout<<"oldContent: "<<oldContent *SF[i]<<" ± "<<oldError<<endl;
        cout<<"NQCD: "<<NQCD<<" ± "<<NQCD_error<<endl;
        cout<<"corrected_rYield: "<<corrected_rYield<<" ± "<<corrected_error<<endl;
        double x1 = TMath::Power(oldContent */

        //cout<<"-----"<<endl;
        hQ_rebinned->SetBinContent(i+1, newContent);
        hQ_rebinned->SetBinError(i+1, newError);
        cout<<"newContent: "<<newContent<<" ± "<<newError<<endl;
        //cout<<"NQCD_error: "<<NQCD_error<<endl;
        //now setThe content for the hSignal
    }
    cout<<"-----"<<endl;
    hSignal->Add(hQ_rebinned,-1);
    hSignal->Add(hSub_rebinned,-1);
    cout<<hD_rebinned->Integral()<<endl;
    cout<<hSub_rebinned->Integral()<<endl;
    cout<<hQ_rebinned->Integral()<<endl;


    cout<<"-------"<<endl;
    hSignal->SetLineColor(kBlue);
    hSMC->SetLineColor(kRed);
    hSignal->SetMarkerStyle(20);
    hSignal->SetMarkerColor(kBlue);

    TLegend *leg = new TLegend(0.65,0.7,0.9,0.9);
    leg->AddEntry(hSignal, "Data", "lpe");
    leg->AddEntry(hSMC, "MC", "lp");

    TCanvas *can = new TCanvas(TString::Format("can_%s",variable.Data()),TString::Format("can_%s",variable.Data()) , 800,600);
    auto *closure_padRatio = new TPad("closure_pad2","closure_pad2",0.,0.,1.,0.3);
    closure_padRatio->Draw();
    closure_padRatio->SetTopMargin(0.05);
    closure_padRatio->SetBottomMargin(0.3);
    closure_padRatio->SetGrid();

    auto *closure_pad1 = new TPad("closure_pad1","closure_pad1",0.,0.3,1.,1.);
    closure_pad1->Draw();
    closure_pad1->SetBottomMargin(0.005);
    closure_pad1->cd();


    TH1F *hSignal_noScale = (TH1F*)hSignal->Clone("hSignal_noScale");
    TH1F *hSMC_noScale = (TH1F*)hSMC->Clone("hSMC_noScale");

    hSignal->Scale(1/luminosity[year],"width");
    hSMC->Scale(1/luminosity[year], "width");

    hSMC->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi} [pb]");
    hSMC->SetTitle(TString::Format("Data vs MC %s for %s ",year.Data(), variable.Data()));
    if(!variable.EqualTo("yJJ") && !variable.EqualTo("jetY0") && !variable.EqualTo("jetY1") ) gPad->SetLogy();

    closure_pad1->cd();
    hSMC->Draw();
    hSignal->Draw("same");
    leg->Draw();

    //theory - data)/data.
    closure_padRatio->cd();
    TH1F *hMCClone[2];
    TH1F *hSig_temp[2];

    hMCClone[0] = (TH1F*)hSMC->Clone("hMCClone0");
    hSig_temp[0] = (TH1F*)hSignal->Clone("hSignal_Clone");

    hSig_temp[0]->SetTitle("");
    hSig_temp[0]->GetYaxis()->SetTitle("#frac{Data}{MC}");
    hSig_temp[0]->GetYaxis()->SetTitleSize(14);
    hSig_temp[0]->GetYaxis()->SetTitleFont(43);
    hSig_temp[0]->GetYaxis()->SetTitleOffset(1.55);
    hSig_temp[0]->GetYaxis()->SetLabelFont(43);
    hSig_temp[0]->GetYaxis()->SetLabelSize(15);
    hSig_temp[0]->GetXaxis()->SetTitleSize(0.09);

    hSig_temp[0]->Divide(hMCClone[0]);
    hSig_temp[0]->SetLineColor(kRed);
    hSig_temp[0]->SetMarkerStyle(20);
    hSig_temp[0]->SetMarkerColor(kRed);
    hSig_temp[0]->Draw();
    hSig_temp[0]->GetXaxis()->SetLabelSize(0.09);

    TString path;
    TString method = "simpleMassFit";
    path = TString::Format("%s/FiducialMeasurement/UnequalBinning/fiducial_%s.pdf",year.Data(),variable.Data());
    can->Print(path,"pdf");


    TFile *outf;
    outf = new TFile(TString::Format("%s/FiducialMeasurement/UnequalBinning/SignalHistograms_%s.root",year.Data(),variable.Data()), "RECREATE");
    hSignal_noScale->Write(TString::Format("hSignal_%s", variable.Data()));
    hSMC_noScale->Write(TString::Format("hSMC_%s", variable.Data()));
    outf->Close();

    cout<<variable.Data()<<endl;
    cout<<hSignal_noScale->Integral()<<endl;


}
