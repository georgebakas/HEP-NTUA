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
#include "../CMS_plots/tdrstyle.C"
#include "../CMS_plots/CMS_lumi.C"

using namespace std;

using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstants.h"
using namespace RooFit;

TH1F *getRebinned(TH1F *h, float BND[], int N);
void SignalExtractionSpecific(TString year = "2016", TString variable = "chi");
bool normalised;

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

void SignalExtraction_SR2TeV(TString year, bool isNormalised)
{
    setTDRStyle();
    normalised = isNormalised;
    TString vars[] = {"chi", "cosTheta_0", "cosTheta_1"};

    int selectedYear;
    if(year.EqualTo("2016")) selectedYear = 0;
    else if(year.EqualTo("2017")) selectedYear = 1;
    else selectedYear =2;

    for(int i=0; i<sizeof(vars)/sizeof(vars[0]); i++)
    {
        SignalExtractionSpecific(year, vars[i]);
        //break;
    }
}


void SignalExtractionSpecific(TString year = "2016", TString variable = "chi")
{
    initFilesMapping();
    cout<<luminosity[year]<<endl;
    RooMsgService::instance().setSilentMode(kTRUE);
    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);

    gStyle->SetOptStat(0);
    //open the signal file: get D(x) and Q(x) for every variable
    TFile *infDataMedium = TFile::Open(TString::Format("%s/Histo_Data_%s_reduced_1500.root", year.Data(), year.Data()));
    //cout<<TString::Format("%s/Histo_Data_%s_100_reduced.root", year.Data(), year.Data())<<endl;
    TH1F *hD = (TH1F*)infDataMedium->Get(TString::Format("hWt_%s_2btag", variable.Data()));
    TH1F *hQ = (TH1F*)infDataMedium->Get(TString::Format("hWt_%s_0btag_expYield", variable.Data()));
    //hQ has to be scaled to integral, because we need the shape

    //open the file to get the Ryield for the SRA--> SR regions
    //------------------------------------------------------------------------------------------------------------------
    TFile *infRyield = TFile::Open(TString::Format("%s/Ryield/TransferFactor_%s.root",year.Data(),variable.Data()));
    //TH1F *hRyield = (TH1F*)infRyield->Get("ClosureTest_TransferFactor");
    TH1F *hRyield = (TH1F*)infRyield->Get("dataTransferFactor");

    float Ryield = hRyield->GetBinContent(1);
    float Ryield_error = hRyield->GetBinError(1);

    float r_yield_correction;
    TH1F *hRyieldMC = (TH1F*)infRyield->Get("ClosureTest_TransferFactor");
    r_yield_correction = (hRyieldMC->GetBinContent(2)/ hRyieldMC->GetBinContent(1));

    float r_yield_correction_error = TMath::Sqrt(TMath::Power(hRyieldMC->GetBinError(2)/hRyieldMC->GetBinContent(1),2) +
                                          TMath::Power((hRyieldMC->GetBinContent(2)*hRyieldMC->GetBinError(2))/TMath::Power(hRyieldMC->GetBinContent(1),2),2));

    cout<<"-------------------------"<<endl;
    cout<<"Ryield_data (0): "<<Ryield<<" ± "<<Ryield_error<<endl;
    cout<<"r_yield_correction: "<<r_yield_correction<<" ± "<<r_yield_correction_error<<endl;

    float corrected_rYield = r_yield_correction * Ryield;
    float corrected_error(0);
    corrected_error = TMath::Sqrt(TMath::Power(r_yield_correction*Ryield_error,2) + TMath::Power(r_yield_correction_error*Ryield,2));
    //cout<<"Ryield_data (2): "<<hRyield->GetBinContent(2)<<endl;
    cout<<"corrected Ryield: "<<corrected_rYield<<" ± "<<corrected_error<<endl;
    //------------------------------------------------------------------------------------------------------------------

    //now calculate transfer factor from mJJ> 1TeV --> 2TeV which is now our SR!
    //------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------
    //open the file to get the Ryield for the SRA--> SR regions
    TFile *infRyield_2TeV = TFile::Open(TString::Format("%s/Ryield/TransferFactor_%s_1.5TeV.root",year.Data(),variable.Data()));
    TH1F *hRyield_2TeV = (TH1F*)infRyield->Get("dataTransferFactor");

    float Ryield_2TeV = hRyield->GetBinContent(1);
    float Ryield_error_2TeV = hRyield->GetBinError(1);

    TH1F *hRyieldMC_2TeV = (TH1F*)infRyield_2TeV->Get("ClosureTest_TransferFactor");
    float r_yield_correction_2TeV = (hRyieldMC_2TeV->GetBinContent(2)/ hRyieldMC_2TeV->GetBinContent(1));

    float r_yield_correction_error_2TeV = TMath::Sqrt(TMath::Power(hRyieldMC_2TeV->GetBinError(2)/hRyieldMC_2TeV->GetBinContent(1),2) +
                                          TMath::Power((hRyieldMC_2TeV->GetBinContent(2)*hRyieldMC_2TeV->GetBinError(2))/TMath::Power(hRyieldMC_2TeV->GetBinContent(1),2),2));

    cout<<"-------------------------"<<endl;
    cout<<"Ryield_data _2TeV(0): "<<Ryield_2TeV<<" ± "<<Ryield_error_2TeV<<endl;
    cout<<"r_yield_correction _2TeV: "<<r_yield_correction_2TeV<<" ± "<<r_yield_correction_error_2TeV<<endl;

    float corrected_rYield_2TeV = r_yield_correction_2TeV * Ryield_2TeV;
    float corrected_error_2TeV(0);
    corrected_error_2TeV = TMath::Sqrt(TMath::Power(r_yield_correction_2TeV*Ryield_error_2TeV,2) + TMath::Power(r_yield_correction_error_2TeV*Ryield_2TeV,2));
    //cout<<"Ryield_data (2): "<<hRyield->GetBinContent(2)<<endl;
    cout<<"corrected Ryield _2TeV: "<<corrected_rYield_2TeV<<" ± "<<corrected_error_2TeV<<endl;
    //------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------



    //open the file to get the Nbkg
    TFile *fitFile = TFile::Open(TString::Format("../MassFit/%s/fitResults_%s.root", year.Data(), year.Data()));
    RooFitResult  *fitResult = (RooFitResult*)fitFile->Get(TString::Format("fitResults_%s", year.Data()));
    //float NQCD = Nbkg2Constants[TString::Format("Nbkg%s",year.Data())];
    //float NQCD_error = Nbkg2ConstantsErrors[TString::Format("Nbkg%s_error",year.Data())];
    RooRealVar *value = (RooRealVar*)fitResult->floatParsFinal().find("nFitQCD_2b");
    float NQCD = value->getVal();
    float NQCD_error = value->getError();

    //Subdominant bkgs files
      TFile *infSub = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_reduced_1500.root", year.Data()));
    TH1F *hSub = (TH1F*)infSub->Get(TString::Format("hWt_%s_2btag_expYield", variable.Data()));
    TH1F *hSub_0 = (TH1F*)infSub->Get(TString::Format("hWt_%s_0btag_expYield", variable.Data()));
    //here I will import correction factors for QCD if needed...

    TFile *infSignalMC;
    infSignalMC = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_reduced_1500.root", year.Data()));
    TH1F *hSMC = (TH1F*)infSignalMC->Get(TString::Format("hWt_%s_2btag_expYield", variable.Data()));
    TH1F *hSMC_0= (TH1F*)infSignalMC->Get(TString::Format("hWt_%s_0btag_expYield", variable.Data()));

    std::vector< std::vector <Float_t> > const BND = {{1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,11.5,13,14.5,16}, //chi
                                                       {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1}, //|cosTheta*| leading
                                                       {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1}};  //|cosTheta*| leading

    int nBins[BND.size()];
    for (int i = 0; i<BND.size(); i++)
    	nBins[i] = BND[i].size()-1;

    //get the respected integer from the mapping in the TemplateConstants.h
    int selVar = variableConstant[variable];

    //use this template for the initialization of the BND source given as input for the rebinned histos
    float tempBND[nBins[selVar]+1];
    std::copy(BND[selVar].begin(), BND[selVar].end(), tempBND);

    //rebin all histograms with the same method
    TH1F *hQ_copy = (TH1F*)hQ->Clone(TString::Format("hD_copy_%s", variable.Data()));
    cout<<hQ_copy->GetNbinsX()<<endl;

    TH1F *hSignal = (TH1F*)hD->Clone(TString::Format("hSignal_%s",variable.Data()));
    //work on the elements for the QCD so that we have the right Q(x)
    //Ryield * Nbkg  * Q(x)
    cout<<variable<<endl;

    float SF[hQ->GetNbinsX()];
    //QCD correction factor for shape
    cout<<variable.Data()<<endl;

    hQ->Scale(1./hQ->Integral());  //this is how you get the shape
    cout<<"--------"<<endl;
    cout<<"NQCD: "<<NQCD<<" ± "<<NQCD_error<<endl;
    cout<<"---"<<endl;
    for(int i =0; i<hQ->GetNbinsX(); i++)
    {
        float oldContent = hQ->GetBinContent(i+1);
        float oldError = hQ->GetBinError(i+1);

        float newContent = oldContent * corrected_rYield * NQCD * corrected_rYield_2TeV;
        float newError   = TMath::Sqrt(TMath::Power(oldError*NQCD*corrected_rYield*corrected_rYield_2TeV/*SF[i]*/,2)+
                                       TMath::Power(NQCD_error*oldContent*corrected_rYield*corrected_rYield_2TeV/*SF[i]*/,2) +
                                       TMath::Power(corrected_error*NQCD*oldContent*corrected_rYield_2TeV/**SF[i]*/,2) +
                                       TMath::Power(corrected_error_2TeV*NQCD*oldContent*corrected_rYield/**SF[i]*/,2));
        //cout<<"-----"<<endl;
        hQ->SetBinContent(i+1, newContent);
        hQ->SetBinError(i+1, newError);
        cout<<"newContent: "<<newContent<<" ± "<<newError<<endl;
        //cout<<
        //now setThe content for the hSignal
    }
    cout<<"-----"<<endl;
    hSignal->Add(hQ,-1);
    hSignal->Add(hSub,-1);
    cout<<hD->Integral()<<endl;
    cout<<hSub->Integral()<<endl;
    cout<<hQ->Integral()<<endl;


    cout<<"-------"<<endl;
    hSignal->SetLineColor(kBlue);
    hSMC->SetLineColor(kRed);
    hSignal->SetMarkerStyle(20);
    hSignal->SetMarkerColor(kBlue);

    TLegend *leg = new TLegend(0.65,0.7,0.9,0.9);
    leg->AddEntry(hSignal, "Data", "lpe");
    leg->AddEntry(hSMC, "MC", "lp");

    hSignal->GetXaxis()->SetTitle(variable);
    hSMC->GetXaxis()->SetTitle(variable);


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

    float sigmaTotalMC = hSMC->Integral() / luminosity[year];
    float sigmaTotal = hSignal->Integral() / luminosity[year];

    hSMC->Scale(1/luminosity[year], "width");
    hSignal->Scale(1/luminosity[year],"width");

    //if you don't want normalised xsec comment these lines
    if(normalised)
    {
      hSignal->Scale(1/sigmaTotal);
      hSMC->Scale(1/sigmaTotalMC);
      hSMC->GetYaxis()->SetTitle("(#frac{1}{#sigma})#frac{d#sigma}{d#chi} [pb]");
    }
    else
    {
      hSMC->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi} [pb]");
    }

    hSMC->SetTitle(TString::Format("Data vs MC %s for %s ",year.Data(), variable.Data()));
    //if(!variable.EqualTo("yJJ") && !variable.EqualTo("jetY0") && !variable.EqualTo("jetY1") ) gPad->SetLogy();

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

    lumi_13TeV = "35.9 fb^{-1}";
    //lumi_sqrtS = "13 TeV";
    int iPeriod = 4;
    int iPos = 10;
    writeExtraText=true;
    CMS_lumi(closure_pad1, iPeriod, iPos);
    
    TString path;
    TString method = "simpleMassFit";
    TString strNorm = "";
    if(normalised) strNorm = "_Norm";
    path = TString::Format("%s/FiducialMeasurement_2TeV/EqualBinning/fiducial_%s%s.pdf",year.Data(),variable.Data(), strNorm.Data());
    can->Print(path,"pdf");


    if(!normalised)
    {
      TFile *outf;
      outf = new TFile(TString::Format("%s/FiducialMeasurement_2TeV/EqualBinning/SignalHistograms_%s.root",year.Data(),variable.Data()), "RECREATE");
      hSignal_noScale->Write(TString::Format("hSignal_%s", variable.Data()));
      hSMC_noScale->Write(TString::Format("hSMC_%s", variable.Data()));
      hQ->Write(TString::Format("hQCD_%s", variable.Data())); //this is expected yield for QCD so I can use it with my stack
      outf->Close();
    }
    cout<<variable.Data()<<endl;
    cout<<hSignal_noScale->Integral()<<endl;


}
