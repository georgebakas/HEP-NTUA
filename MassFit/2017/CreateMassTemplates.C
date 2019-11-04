TH1F *GetTruncatedHisto(TH1F *hIN,float XMIN,float XMAX);
TH1F *GetTruncatedHisto(TH1F *hIN,float XMIN,float XMAX)
{
  int N = (XMAX-XMIN)/hIN->GetBinWidth(1);
  TH1F *hOUT = new TH1F(TString(hIN->GetName())+"_Truncated",TString(hIN->GetName())+"_Truncated",N,XMIN,XMAX);
  hOUT->Sumw2();
  for(int i=0;i<hIN->GetNbinsX();i++) {
    int bin = hOUT->FindBin(hIN->GetBinCenter(i+1));
    hOUT->SetBinContent(bin,hIN->GetBinContent(i+1));
    hOUT->SetBinError(bin,hIN->GetBinError(i+1));
  }
  return hOUT;
}

void CreateMassTemplates()
{
  gROOT->ForceStyle();
  
  float normMC;
  float XSEC(832.);
  float LUMI(37000.);
  const int N = 18; 
  float norm;

  TFile *inf_bkgSF = TFile::Open("../ScaleFactor_Cut5_mTop_boosted.root");
  TF1 *bkgSF = (TF1*)inf_bkgSF->Get("fit_Cut5_2btag");

  TString SAMPLE[N] = {
    "JetHT",//data
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",//nominal
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",//JESUp
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",//JESDown
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",//JERUp
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",//JERDown
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",//BtagUp
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",//BtagDown
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",//PileupUp
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",//PileupDown 
    "TT_TuneCUETP8M2T4_13TeV-powheg-fsrup-pythia8",//FSRUp
    "TT_TuneCUETP8M2T4_13TeV-powheg-fsrdown-pythia8",//FSRDown
    "TT_TuneCUETP8M2T4_13TeV-powheg-isrup-pythia8",//ISRUp
    "TT_TuneCUETP8M2T4_13TeV-powheg-isrdown-pythia8",//ISRDown
    "TT_TuneCUETP8M2T4up_13TeV-powheg-pythia8",//TuneUp
    "TT_TuneCUETP8M2T4down_13TeV-powheg-pythia8",//TuneDown
    "TT_hdampUP_TuneCUETP8M2T4_13TeV-powheg-pythia8",//hdampUp
    "TT_hdampDOWN_TuneCUETP8M2T4_13TeV-powheg-pythia8"//hdampDown
  };

  TString ALIAS[N] = {"data","signal","signal_JESUp","signal_JESDown","signal_JERUp","signal_JERDown",
                      "signal_BtagUp","signal_BtagDown","signal_PileupUp","signal_PileupDown",
                      "signal_FSRUp","signal_FSRDown","signal_ISRUp","signal_ISRDown",
                      "signal_TuneUp","signal_TuneDown","signal_hdampUp","signal_hdampDown"};

  TString WEIGHT[N] = {"","Wt","Wt","Wt","Wt","Wt","WtBtagUp","WtBtagDown","WtPileupUp","WtPileupDown",
  "Wt","Wt","Wt","Wt","Wt","Wt","Wt","Wt"};

  TString DIR[N] = {
     "boosted","boosted","boostedShiftedUp","boostedShiftedDown",
     "boostedSmearedUp","boostedSmearedDown","boosted","boosted","boosted","boosted","boosted",
     "boosted","boosted","boosted","boosted","boosted","boosted","boosted"};

  TFile *inf[N];
  TH1F *h[N][3],*hRange[N][3];
  TH1F *bkg,*bkg_ShapeUp,*bkg_ShapeDown;

  TFile *outf = new TFile("mTop_shapes.root","RECREATE");

  for(int i=0;i<N;i++) {
    cout<<ALIAS[i]<<" "<<SAMPLE[i]<<endl;
    inf[i] = TFile::Open("../Histo_"+SAMPLE[i]+".root");
    for(int icat=0;icat<3;icat++) {
      TString CAT = TString::Format("%dbtag",icat);
      h[i][icat] = (TH1F*)inf[i]->Get(DIR[i]+"/h"+WEIGHT[i]+"_mTop_Cut5_"+CAT);
      h[i][icat]->Rebin(5);
      hRange[i][icat] = (TH1F*)GetTruncatedHisto(h[i][icat],50,300);
      if (i==0) {
        if (icat==0) {
          outf->cd();
          bkg           = (TH1F*)hRange[i][icat]->Clone("background");
          bkg_ShapeUp   = (TH1F*)hRange[i][icat]->Clone("background_BkgShapeUp");
          bkg_ShapeDown = (TH1F*)hRange[i][icat]->Clone("background_BkgShapeDown"); 
          for(int ibin=0;ibin<bkg->GetNbinsX();ibin++) {
            float y = hRange[i][icat]->GetBinContent(ibin+1);
            float e = hRange[i][icat]->GetBinError(ibin+1);
            float c = bkgSF->Eval(hRange[i][icat]->GetBinCenter(ibin+1));
            float c_up = c+fabs(1-c);
            float c_down = c-fabs(1-c);
            bkg->SetBinContent(ibin+1,y*c);
            bkg->SetBinError(ibin+1,e*c);  
            //bkg_ShapeUp->SetBinContent(ibin+1,y*(2.5*c-1.5));
            //bkg_ShapeUp->SetBinError(ibin+1,e*(2.5*c-1.5));
            bkg_ShapeUp->SetBinContent(ibin+1,y*c_up);
            bkg_ShapeUp->SetBinError(ibin+1,e*c_up);
            bkg_ShapeDown->SetBinContent(ibin+1,y*c_down);
            bkg_ShapeDown->SetBinError(ibin+1,e*c_down);
          }
          bkg->Write();
          bkg_ShapeUp->Write();
          bkg_ShapeDown->Write();          
        }
        if (i==0 && icat==2) {
          outf->cd();
          hRange[i][icat]->Write("data_obs");
        }
      }
      else {
        norm = ((TH1F*)inf[i]->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
        if (icat==2) {
          outf->cd();
          cout<<"Nexp = "<<(LUMI*XSEC/norm)*hRange[i][icat]->Integral()<<endl;
          hRange[i][icat]->Scale(LUMI*XSEC/norm);
          hRange[i][icat]->Write(ALIAS[i]);
        }
      }
      
    }
    inf[i]->Close();  
  }
  outf->Close();
}                            

