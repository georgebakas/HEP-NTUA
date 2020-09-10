#include "TemplateConstants.h"
using namespace RooFit;

void drawPlot(TString nuisance, float values[3], float errors[3], float values_noBtag[3], float errors_noBtag[3]);

void NuisancesComparison()
{
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  const int N =3;
  TString years[N] = {"2016", "2017", "2018"};
  TFile *file[N], *file_noBtagSF[N];
  RooWorkspace *w[N], *w_noBtagSF[N];
  float yieldTT[N], yieldQCD[N], yieldSub[N], kMassScale[N], kMassResol[N], kQCD[N];
  float yieldTT_e[N], yieldQCD_e[N], yieldSub_e[N], kMassScale_e[N], kMassResol_e[N], kQCD_e[N];

  const int N_nuisances = 6;
  TString nuisances[N_nuisances] = {"nFitSig2b", "nFitQCD_2b","nFitBkg_2b","kMassScale","kMassResol","kQCD_2b"};
  float nuisanceValues[N_nuisances][N], nuisanceErrors[N_nuisances][N];
  float nuisanceValues_noBtag[N_nuisances][N], nuisanceErrors_noBtag[N_nuisances][N];
  for(int iy=0; iy<N; iy++)
  {
    file[iy] = TFile::Open(TString::Format("%s/MassFitResults__.root", years[iy].Data()));
    file_noBtagSF[iy] = TFile::Open(TString::Format("%s_noBTagSF/MassFitResults__.root", years[iy].Data()));
    w[iy] = (RooWorkspace*)file[iy]->Get("w");
    w_noBtagSF[iy] = (RooWorkspace*)file_noBtagSF[iy]->Get("w");

    for(int in=0; in<N_nuisances; in++)
    {
      nuisanceValues[in][iy]= ((RooRealVar*)w[iy]->var(nuisances[in].Data()))->getVal();
      nuisanceErrors[in][iy]= ((RooRealVar*)w[iy]->var(nuisances[in].Data()))->getError();
      nuisanceValues_noBtag[in][iy]= ((RooRealVar*)w_noBtagSF[iy]->var(nuisances[in].Data()))->getVal();
      nuisanceErrors_noBtag[in][iy]= ((RooRealVar*)w_noBtagSF[iy]->var(nuisances[in].Data()))->getError();
    }

  }
  //create temp arrays that get each nuisance for each year
  float tempValue[N], tempError[N];
  float tempValue_noBtag[N], tempError_noBtag[N];
  for(int in=0; in<N_nuisances;in++)
  {
    //assign the value for each year
    for(int iy=0; iy<N; iy++)
    {
      tempValue[iy]=nuisanceValues[in][iy];
      tempError[iy]=nuisanceErrors[in][iy];
      tempValue_noBtag[iy]=nuisanceValues_noBtag[in][iy];
      tempError_noBtag[iy]=nuisanceErrors_noBtag[in][iy];
    }
    cout<<"Drawing "<<nuisances[in]<<endl;
    drawPlot(nuisances[in], tempValue, tempError, tempValue_noBtag, tempError_noBtag);
  }

  //return;

}

void SignalStrengthComparison()
{
  initFilesMapping();
  float sigStrength[3] = {
    ttbarSigStrength["2016"],
    ttbarSigStrength["2017"],
    ttbarSigStrength["2018"]
  };
  //get errors
  float sigStrengthError[3] = {
    ttbarSigStrengthError["2016"],
    ttbarSigStrengthError["2017"],
    ttbarSigStrengthError["2018"]
  };

  float sigStrength_noBTagSF[3] = {
    ttbarSigStrength_noBTagSF["2016"],
    ttbarSigStrength_noBTagSF["2017"],
    ttbarSigStrength_noBTagSF["2018"]
  };

  float sigStrengthError_noBTagSF[3] = {
    ttbarSigStrengthError_noBTagSF["2016"],
    ttbarSigStrengthError_noBTagSF["2017"],
    ttbarSigStrengthError_noBTagSF["2018"]
  };
  drawPlot("t#bar{t} Sig.Strength", sigStrength, sigStrengthError, sigStrength_noBTagSF, sigStrengthError_noBTagSF);
}


void drawPlot(TString nuisance, float values[3], float errors[3], float values_noBtag[3], float errors_noBtag[3])
{

  const int n = 3;
  float e_x[n] = {0,0,0};
  float x[n] = {0.5,1.5,2.5};
  float x_noBTag[n] = {0.7,1.7,2.7};

  TGraphErrors *gr = new TGraphErrors(n,x, values, e_x, errors);
  TGraphErrors *gr_noBTagSF = new TGraphErrors(n,x_noBTag, values_noBtag,e_x, errors_noBtag);

  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(22);
  gr_noBTagSF->SetMarkerColor(kBlue);
  gr_noBTagSF->SetMarkerStyle(20);

  //use TMultiGraph to draw the two graphs together
  TCanvas *can = new TCanvas("can"+nuisance, "can"+nuisance, 800, 600);
  TLegend *leg = new TLegend(0.4,0.88,0.6,0.98);
  leg->AddEntry(gr, "bTag SFs", "lep");
  leg->AddEntry(gr_noBTagSF, "NO bTag SFs", "lep");

  TString Xname[n] = {"2016","2017","2018"};
  //gr->GetHistogram()->GetXaxis()->Set(3, 0, 3);
  //gr_noBTagSF->GetHistogram()->GetXaxis()->Set(3, 0, 3);
  //for(int i=0; i<n; i++){
  //  gr->GetHistogram()->GetXaxis()->SetBinLabel(i+1, Xname[i]);
  //  gr_noBTagSF->GetHistogram()->GetXaxis()->SetBinLabel(i+1, Xname[i]);
  //}


  TMultiGraph  *mg  = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr_noBTagSF);
  mg->GetHistogram()->GetXaxis()->Set(3, 0, 3);
  for(int i=0; i<n; i++){
    mg->GetHistogram()->GetXaxis()->SetBinLabel(i+1, Xname[i]);
  }
  mg->GetYaxis()->SetTitle(nuisance);

  if(nuisance.EqualTo("t#bar{t} Sig.Strength"))mg->GetYaxis()->SetRangeUser(0.5, 0.8);

  mg->Draw("AP");
  leg->Draw();

  can->Print(TString::Format("NuisancePlots_BeforeAfterBtagSFs/Comparison_%s.pdf",nuisance.Data()), "pdf");

}
