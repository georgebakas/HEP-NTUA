/*
  This is a file that creates parton/particle templates so that we can compare the results with the unfolded
  Parton selection:
  pt > 400
  |eta| < 2.4
  mTTbarParton > 1000

  Particle:
  Njets > 1
  Pt > 400
  |eta| < 2.4
  jetMassSoftDrop > 120 && jetMassSoftDrop < 220
  mJJ > 1000
*/

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TemplateConstantsUnfold.h"
using namespace std;

using std::cin;
using std::cout;
using std::endl;

std::vector<float> XSEC;
std::vector<TString> histoNames;
std::vector<TString> fileNames;
TString globalYear;
bool globalIsNominalMC;

void initXsections()
{
  if(!globalIsNominalMC)
  {
  	XSEC.push_back(69.64);
  	XSEC.push_back(16.74);
  }
  else
  {
   	if(globalYear.EqualTo("2016"))
   		XSEC.push_back(832.);
   	else
   	{
   		XSEC.push_back(377.96);
   		XSEC.push_back(365.34);
      XSEC.push_back(88.29);
 	}
  }
}

void initHistoNames()
{
  if(!globalIsNominalMC)
  {
	  histoNames.push_back("Signal_histo_Mtt_700_1000"); 
	  histoNames.push_back("Signal_histo_Mtt_1000_Inf");

	  fileNames.push_back("700-1000");
	  fileNames.push_back("1000-Inf");
  }
  else
  {
	  if(globalYear.EqualTo("2016"))
	  {
	  	fileNames.push_back("TTNominal");
	  	histoNames.push_back("Signal_histo_Nominal");
	  }
	  else
	  {
	  	fileNames.push_back("TTHadronic_0");
	  	fileNames.push_back("TTSemiLeptonic_0");
      fileNames.push_back("TTTo2L2Nu_0");

	  	histoNames.push_back("Signal_histo_TTHadronic_0");
	  	histoNames.push_back("Signal_histo_TTSemiLeptonic_0");
      histoNames.push_back("Signal_histo_TTTo2L2Nu_0");

	  }
  }
}

void CreateTheoryTemplates(TString inYear = "2016", bool isNominalMC= true)
{
  globalIsNominalMC = isNominalMC;
  globalYear = inYear;
  initFilesMapping();
  initFilesMapping();
  initHistoNames();
  initXsections();
  TString year = inYear;

  std::vector< std::vector <Float_t> > const BND = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}, //mjj
                                                        {0,60,150,300,450,600,750,950,1100,1300}, //ptjj
                                                        {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
                                                        {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt0     
                                                        {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt1
                                                        {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
                                                        {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1
  float LUMI = luminosity[year];
  //define histograms:
  int NBINS[BND.size()];
  const int NVAR = 7;
  for (int i = 0; i<BND.size(); i++) NBINS[i] = BND[i].size()-1;
  TString var[NVAR]  = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1"}; 
  TH1F *hParton[fileNames.size()][NVAR], *hParticle[fileNames.size()][NVAR];
  std::vector<float> weights;
  //loop on all files:
  for(int f=0; f<fileNames.size(); f++)
  {

  for(int ivar =0; ivar<NVAR; ivar++)
  {	
  		 int sizeBins = NBINS[ivar];
         float tempBND[NBINS[ivar]+1];
         std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND);
  		 
  		 //denominators for parton efficiency (hParton vs parton), particle eff (hParticle vs particle) and acceptance (same for both hReco vs reco)
  		 hParton[f][ivar] = new TH1F(TString::Format("hParton_%s_%s",histoNames[f].Data(),var[ivar].Data()), TString::Format("hParton_%s_%s",histoNames[f].Data(),var[ivar].Data()), sizeBins, tempBND);
         hParticle[f][ivar] = new TH1F(TString::Format("hParticle_%s_%s",histoNames[f].Data() ,var[ivar].Data()), TString::Format("hParticle_%s_%s",histoNames[f].Data() ,var[ivar].Data()), sizeBins, tempBND);
         hParton[f][ivar]->Sumw2();
         hParticle[f][ivar]->Sumw2();
  }
  cout<<"ok with histograms!"<<endl;
  

  TFile *file = TFile::Open(eospath[year.Data()]+files[year.Data()][fileNames[f].Data()]);
  TTree *trIN = (TTree*)file->Get("boosted/events");
  cout<<"working in file:"<<eospath[year.Data()]+files[year.Data()][fileNames[f].Data()]<<endl;
  //parton
  //std::vector<float> *partonPt(0), *partonPhi(0), *partonY(0), *partonMass(0),*partonMatchDR(0),*partonEta(0);
  float genEvtWeight(0);
  
  float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  float weight = XSEC[f]/norm;
  weights.push_back(weight);
  
  //particle
  std::vector<float> *genjetPt(0), *genjetY(0), *genjetEta(0), *genjetMassSoftDrop(0);
  int nJetsGen(0);
  float mJJGen(0), ptJJGen(0), yJJGen(0);

  trIN->SetBranchAddress("nJetsGen" ,&nJetsGen);
  trIN->SetBranchAddress("mJJGen"   ,&mJJGen);
  trIN->SetBranchAddress("ptJJGen"  ,&ptJJGen);
  trIN->SetBranchAddress("yJJGen"   ,&yJJGen);
  trIN->SetBranchAddress("genjetPt" ,&genjetPt);
  trIN->SetBranchAddress("genjetEta",&genjetEta);
  trIN->SetBranchAddress("genjetY"  ,&genjetY);
  trIN->SetBranchAddress("genjetMassSoftDrop", &genjetMassSoftDrop);
  trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);   

  bool particleCuts;
  
  long NN = trIN->GetEntries();
  //NN = 1000;
  std::cout<<"Entries: "<<NN<<std::endl;
  std::vector<float> xParticleAll(0);
  int decade(0);

  for(int iev=0;iev<NN;iev++) 
  {
   double progress = 10.0*iev/(1.0*NN);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
        cout<<10*k<<" %"<<endl;
      decade = k;
   trIN->GetEntry(iev);  
   xParticleAll.clear();
   int genSubleadingPt = 1;
   int genLeadingPt = 0;
 // cout<<"ok"<<endl;  
  if(nJetsGen>1)
   {

    if((*genjetPt)[0] < (*genjetPt)[1])
    {
        genSubleadingPt =0;
        genLeadingPt = 1;
    }
	  xParticleAll.push_back(mJJGen);
	  xParticleAll.push_back(ptJJGen);
	  xParticleAll.push_back(yJJGen);
	  xParticleAll.push_back((*genjetPt)[genLeadingPt]);
	  xParticleAll.push_back((*genjetPt)[genSubleadingPt]);
      xParticleAll.push_back(fabs((*genjetY)[genLeadingPt]));
      xParticleAll.push_back(fabs((*genjetY)[genSubleadingPt]));
  
    particleCuts = fabs((*genjetEta)[0]) < 2.4 && fabs((*genjetEta)[1]) < 2.4 && (*genjetPt)[0] > 400 && (*genjetPt)[1] > 400 && mJJGen > 1000 && nJetsGen >1 &&
      (*genjetMassSoftDrop)[0] > 120 && (*genjetMassSoftDrop)[0] < 220 && (*genjetMassSoftDrop)[1] > 120 && (*genjetMassSoftDrop)[1] < 220;
    
    if(particleCuts)
    { 
		for(int ivar = 0; ivar < xParticleAll.size(); ivar++) 
		{
		  //cout<<xParticleAll[ivar]<<endl; 
	      hParticle[f][ivar]->Fill(xParticleAll[ivar], genEvtWeight);
	    }
    }
  }//nJetsGen
  } //end of tree entries for-loop for gen level

  
  //loop over other tree -> eventCounter
  TTree *trCnt = (TTree*)file->Get("eventCounter/events");
  float ptTTbarPartonCnt(0), mTTbarPartonCnt(0), yTTbarPartonCnt(0);
  float partonPtCnt[2], partonEtaCnt[2],partonYCnt[2];  
  float genEvtWeightCnt;
  //tree for eventCounter   
  trCnt->SetBranchAddress("ptTopParton"    ,&partonPtCnt);
  trCnt->SetBranchAddress("etaTopParton"   ,&partonEtaCnt);
  trCnt->SetBranchAddress("yTopParton"     ,&partonYCnt);
  trCnt->SetBranchAddress("mTTbarParton"   ,&mTTbarPartonCnt);
  trCnt->SetBranchAddress("yTTbarParton"   ,&yTTbarPartonCnt);
  trCnt->SetBranchAddress("ptTTbarParton"  ,&ptTTbarPartonCnt);
  trCnt->SetBranchAddress("genEvtWeight"   ,&genEvtWeightCnt);   
  int NNCnt = trCnt->GetEntries();
  cout<<"NNCnt: "<<NNCnt<<endl;
  for(int iev = 0; iev < NNCnt; iev++)
  {
    trCnt->GetEntry(iev);
    std::vector<float> xPartonAllCnt(0);
    xPartonAllCnt.clear();
    bool partonCuts = fabs(partonEtaCnt[0]) < 2.4 && fabs(partonEtaCnt[1]) <2.4 && partonPtCnt[0] > 400 && partonPtCnt[1] > 400 && mTTbarPartonCnt > 1000;
    
    xPartonAllCnt.push_back(mTTbarPartonCnt);
    xPartonAllCnt.push_back(ptTTbarPartonCnt);
    xPartonAllCnt.push_back(yTTbarPartonCnt);

    if(partonPtCnt[0] >= partonPtCnt[1])
    {
      xPartonAllCnt.push_back(partonPtCnt[0]);
      xPartonAllCnt.push_back(partonPtCnt[1]);

      xPartonAllCnt.push_back(fabs(partonYCnt[0]));
      xPartonAllCnt.push_back(fabs(partonYCnt[1]));
    }
    else
    {
      xPartonAllCnt.push_back(partonPtCnt[1]);
      xPartonAllCnt.push_back(partonPtCnt[0]);

      xPartonAllCnt.push_back(fabs(partonYCnt[1]));
      xPartonAllCnt.push_back(fabs(partonYCnt[0]));
    }


    for(int ivar = 0; ivar < xPartonAllCnt.size(); ivar++)
    {
      if(partonCuts) hParton[f][ivar]->Fill(xPartonAllCnt[ivar], genEvtWeightCnt);
    }
  }//end of for -eventCounter loop
  

  }// end of files loop
  
  //for all vars in all files, scale
  for(int ivar=0; ivar<NVAR; ivar++) 
  {
	  hParton[0][ivar]->Scale(weights[0]*LUMI);
	  hParticle[0][ivar]->Scale(weights[0]*LUMI);
	  cout<<hParticle[0][ivar]->GetEntries()<<endl;
	  for(int f=1; f<fileNames.size(); f++) 
	  {
	    hParton[f][ivar]->Scale(weights[f]*LUMI);
		  hParton[0][ivar]->Add(hParton[f][ivar]);

		  hParticle[f][ivar]->Scale(weights[f]*LUMI);
		  hParticle[0][ivar]->Add(hParticle[f][ivar]);
	  }
  }


  //scale to luminosity and bin width
  for(int ivar=0; ivar<NVAR; ivar++) 
  {
    hParton[0][ivar]->Scale(1/LUMI, "width");
    hParticle[0][ivar]->Scale(1/LUMI, "width");
  }

  TString nominal = "NominalMC";
  if(!globalIsNominalMC) nominal = "HighMtt";
  TFile *outf = new TFile(TString::Format("%s/TheoryTemplates%s.root", year.Data(), nominal.Data()), "RECREATE");
  outf->cd();
  //write them to file
  for(int ivar=0; ivar<NVAR; ivar++) 
  {
    hParton[0][ivar]->Write(TString::Format("hParton_%s", var[ivar].Data()));
    hParticle[0][ivar]->Write(TString::Format("hParticle_%s", var[ivar].Data()));
  }

  outf->Close();
  fileNames.clear();
  histoNames.clear();
  XSEC.clear();
}

