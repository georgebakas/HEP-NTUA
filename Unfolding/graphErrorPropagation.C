
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TGraphErrors.h"

TH1F *getRebinned(TH1F *h, float BND[], int N);

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
void graphErrorPropagation(TString year="2016", bool isParton = true, bool isData=false)
{
   std::vector< std::vector <Float_t> > const BND_gen = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}, //mjj
                                                         {0,60,150,300,450,600,750,1000,1300}, //ptjj
                                                         {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
                                                         {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt0
                                                         {400,450,500,570,650,750,850,1000,1200,1500}, //jetPt1
                                                         {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
                                                         {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1
    int NBINS[BND_gen.size()];
  	for (int i = 0; i<BND_gen.size(); i++)
  		NBINS[i] = BND_gen[i].size()-1;

	TString phaseSpace = "Particle";
	if(isParton) phaseSpace = "Parton";

	TString dataMC = "MC";
	if(isData) dataMC = "Data";

	TFile *infStd = TFile::Open(TString::Format("%s/%sMeasurements/%s/OutputFile.root",year.Data(),phaseSpace.Data(),dataMC.Data()));
	TFile *infRho = TFile::Open(TString::Format("%s/%sMeasurements/%s/OutputFile_RhoMethod.root",year.Data(),phaseSpace.Data(),dataMC.Data()));
    TFile *infLCurve = TFile::Open(TString::Format("%s/%sMeasurements/%s/OutputFile_LCurveMethod.root",year.Data(),phaseSpace.Data(),dataMC.Data()));

	//hErrorAfter[ivar]->Write(TString::Format("hErrorAfter_%s", variable[ivar].Data()));
    //hErrorBefore[ivar]->Write(TString::Format("hErrorBefore_%s", variable[ivar].Data()));

    const int NVAR = 7;
    TString variable[NVAR] = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1","jetY0", "jetY1"};

    TH1F *hErrorBeforeStd[NVAR], *hErrorBeforeRho[NVAR], *hErrorBeforeLCurve[NVAR];
    TH1F *hErrorBeforeStdRebinned[NVAR], *hErrorBeforeRhoRebinned[NVAR], *hErrorBeforeLCurveRebinned[NVAR];
    TH1F *hErrorAfterStd[NVAR], *hErrorAfterRho[NVAR], *hErrorAfterLCurve[NVAR];

    TGraph *globalCorrGraph[NVAR];
    //TFile *globalCorrFile = TFile::Open(TString::Format("%s/%sMeasurements/%s/OutputFileRhoGraphs.root",year.Data(), phaseSpace.Data(), dataMC.Data()));

    TCanvas *canStd[NVAR], *canRho[NVAR], *canGr[NVAR], *canCombined[NVAR], *canLCurve[NVAR];
    TLegend *legComb[NVAR];

    for(int ivar=0; ivar<NVAR; ivar++)
    {
    	int sizeBins = NBINS[ivar];
    	float tempBND[NBINS[ivar]+1];
    	std::copy(BND_gen[ivar].begin(), BND_gen[ivar].end(), tempBND);

    	hErrorBeforeStd[ivar] = (TH1F*)infStd->Get(TString::Format("hErrorBefore_%s", variable[ivar].Data()));
    	hErrorBeforeStd[ivar]->SetName(TString::Format("hErrorBeforeStd_%s", variable[ivar].Data()));
    	hErrorBeforeRho[ivar] = (TH1F*)infRho->Get(TString::Format("hErrorBefore_%s", variable[ivar].Data()));
    	cout<<hErrorBeforeRho[ivar]->GetNbinsX()<<endl;
    	hErrorBeforeRho[ivar]->SetName(TString::Format("hErrorBeforeRho_%s", variable[ivar].Data()));
        hErrorBeforeLCurve[ivar] = (TH1F*)infLCurve->Get(TString::Format("hErrorBefore_%s", variable[ivar].Data()));
        hErrorBeforeLCurve[ivar]->SetName(TString::Format("hErrorBeforeLCurve_%s", variable[ivar].Data()));

    	hErrorAfterStd[ivar] = (TH1F*)infStd->Get(TString::Format("hErrorAfter_%s", variable[ivar].Data()));
    	hErrorAfterRho[ivar] = (TH1F*)infRho->Get(TString::Format("hErrorAfter_%s", variable[ivar].Data()));
        hErrorAfterLCurve[ivar] = (TH1F*)infLCurve->Get(TString::Format("hErrorAfter_%s", variable[ivar].Data()));
    	hErrorAfterStd[ivar]->GetYaxis()->SetTitle("#frac{Error After}{Error Before}");
    	hErrorAfterRho[ivar]->GetYaxis()->SetTitle("#frac{Error After}{Error Before}");
    	hErrorAfterLCurve[ivar]->GetYaxis()->SetTitle("#frac{Error After}{Error Before}");
    	//hErrorAfterStd[ivar]->GetYaxis()->SetTitle("Error After");
    	//hErrorAfterRho[ivar]->GetYaxis()->SetTitle("Error After");
        //hErrorAfterLCurve[ivar]->GetYaxis()->SetTitle("Error After");
    	hErrorAfterRho[ivar]->SetLineColor(kBlue);
    	hErrorAfterRho[ivar]->SetMarkerColor(kBlue);
    	hErrorAfterRho[ivar]->SetMarkerStyle(23);
    	hErrorAfterStd[ivar]->SetLineColor(kRed);
    	hErrorAfterStd[ivar]->SetMarkerColor(kRed);
    	hErrorAfterStd[ivar]->SetMarkerStyle(20);

        hErrorAfterLCurve[ivar]->SetLineColor(kGreen);
        hErrorAfterLCurve[ivar]->SetMarkerColor(kGreen);
        hErrorAfterLCurve[ivar]->SetMarkerStyle(21);

        
    	if(!variable[ivar].EqualTo("yJJ"))
    	{
    		hErrorBeforeStdRebinned[ivar]=getRebinned(hErrorBeforeStd[ivar], tempBND, NBINS[ivar]);
    		hErrorBeforeRhoRebinned[ivar]=getRebinned(hErrorBeforeRho[ivar], tempBND, NBINS[ivar]);
            hErrorBeforeLCurveRebinned[ivar]=getRebinned(hErrorBeforeLCurve[ivar], tempBND, NBINS[ivar]);

    		hErrorAfterStd[ivar]->Divide(hErrorBeforeStdRebinned[ivar]);
    		hErrorAfterRho[ivar]->Divide(hErrorBeforeRhoRebinned[ivar]);
            hErrorAfterLCurve[ivar]->Divide(hErrorBeforeLCurveRebinned[ivar]);
    	}
    	else
		{
			hErrorAfterStd[ivar]->Divide(hErrorBeforeStd[ivar]);
    		hErrorAfterRho[ivar]->Divide(hErrorBeforeRho[ivar]);
            hErrorAfterLCurve[ivar]->Divide(hErrorBeforeRho[ivar]);
		}
        
/*
    	canStd[ivar] = new TCanvas(TString::Format("canStd_%s", variable[ivar].Data()), TString::Format("canStd_%s", variable[ivar].Data()), 800, 600);
    	canStd[ivar]->cd();
    	hErrorAfterStd[ivar]->Draw();

    	canStd[ivar]->Print(TString::Format("%s/%sMeasurements/%s/Errors/ErrorProp_%sStd.pdf",year.Data(),phaseSpace.Data(),dataMC.Data(),
    										 variable[ivar].Data()), "pdf");

    	canRho[ivar] = new TCanvas(TString::Format("canRho_%s", variable[ivar].Data()), TString::Format("canRho_%s", variable[ivar].Data()), 800, 600);
    	canRho[ivar]->cd();
    	hErrorAfterRho[ivar]->Draw();

    	canRho[ivar]->Print(TString::Format("%s/%sMeasurements/%s/Errors/ErrorProp_%sRhoMethod.pdf",year.Data(),phaseSpace.Data(),dataMC.Data(),
    										 variable[ivar].Data()), "pdf");

        canLCurve[ivar] = new TCanvas(TString::Format("canLCurve_%s", variable[ivar].Data()), TString::Format("canLCurve_%s", variable[ivar].Data()), 800, 600);
        canLCurve[ivar]->cd();
        hErrorAfterLCurve[ivar]->Draw();

        canLCurve[ivar]->Print(TString::Format("%s/%sMeasurements/%s/Errors/ErrorProp_%sLCurveMethod.pdf",year.Data(),phaseSpace.Data(),dataMC.Data(),
                                             variable[ivar].Data()), "pdf");
*/
    canCombined[ivar] = new TCanvas(TString::Format("canRhoComb_%s", variable[ivar].Data()), TString::Format("canRhoComb_%s", variable[ivar].Data()), 800, 600);
    legComb[ivar] = new TLegend(0.65,0.7,0.9,0.9);
    canCombined[ivar]->cd();
    legComb[ivar]->AddEntry(hErrorAfterRho[ivar],"Global Corr.", "l");
    legComb[ivar]->AddEntry(hErrorAfterStd[ivar],"Simple Matrix Inv.", "l");
    legComb[ivar]->AddEntry(hErrorAfterLCurve[ivar],"Scan LCurve", "l");
    hErrorAfterStd[ivar]->GetYaxis()->SetRangeUser(0.5,1.5);
    hErrorAfterStd[ivar]->Draw("hist ");
    hErrorAfterRho[ivar]->Draw("hist same");
    hErrorAfterLCurve[ivar]->Draw("hist same");
    hErrorAfterStd[ivar]->GetYaxis()->SetRangeUser(0.5,1.5);
    
    //cout<<"-----"<<variable[ivar]<<"-----"<<endl;
    //cout<<hErrorAfterLCurve[ivar]->GetNbinsX()<<endl;
    //cout<<hErrorAfterRho[ivar]->GetNbinsX()<<endl;
    //cout<<hErrorAfterStd[ivar]->GetNbinsX()<<endl;
    legComb[ivar]->Draw();

    canCombined[ivar]->Print(TString::Format("%s/%sMeasurements/%s/Errors/ErrorPropCombinedRatio_%s.pdf",year.Data(),phaseSpace.Data(),dataMC.Data(),
    										 variable[ivar].Data()), "pdf");
 /*
    globalCorrGraph[ivar] = (TGraph*)globalCorrFile->Get(TString::Format("globalCorrGraph_%s",variable[ivar].Data()));
 		canGr[ivar] = new TCanvas (TString::Format("globalCorrGraph_%s",variable[ivar].Data()), TString::Format("globalCorrGraph_%s",variable[ivar].Data()),
                                 800,600);
 		canGr[ivar]->cd();
  		globalCorrGraph[ivar]->Draw();
  		canGr[ivar]->Print(TString::Format("%s/%sMeasurements/%s/GlobalCorrelationGraph_%s.pdf",year.Data(),phaseSpace.Data(),dataMC.Data(),
    										 variable[ivar].Data()), "pdf");
*/
    }

}
