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


void CreateTheoryTemplates(TString inYear = "2016", bool isParton = true)
{
  initFilesMapping();
  TString year = inYear;

  std::vector< std::vector <Float_t> > const BND = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}, //mjj
                                                        {0,60,150,300,450,600,750,950,1100,1300}, //ptjj
                                                        {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
                                                        {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt0     
                                                        {400,450,500,570,650,750,850,950,1100,1300,1500}}; //jetPt1
                                                        //{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
                                                        //{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1
  float LUMI = luminosity[year];
  //define histograms:
  int NBINS[BND.size()];
  const int NVAR = 7;
  for (int i = 0; i<BND.size(); i++) NBINS[i] = BND[i].size()-1;
  TString var[NVAR]  = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1"}; 

  TH1F *hParton[NVAR], *hParticle[NVAR];
  for(int ivar =0; ivar<NVAR-2; ivar++)
  {	
  		 int sizeBins = NBINS[ivar];
         float tempBND[NBINS[ivar]+1];
         std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND);
  		 
  		 //denominators for parton efficiency (hParton vs parton), particle eff (hParticle vs particle) and acceptance (same for both hReco vs reco)
  		 hParton[ivar] = new TH1F(TString::Format("hParton_%s", var[ivar].Data()), TString::Format("hParton_%s",var[ivar].Data()), sizeBins, tempBND);
         hParticle[ivar] = new TH1F(TString::Format("hParticle_%s", var[ivar].Data()), TString::Format("hParticle_%s", var[ivar].Data()), sizeBins, tempBND);
         hParton[ivar]->Sumw2();
         hParticle[ivar]->Sumw2();
  }
 cout<<"ok with histograms!"<<endl;
  TFile *inf = TFile::Open("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Signal/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");


  TTree *trIN = (TTree*)inf->Get("boosted/events");
  //parton
  std::vector<float> *partonPt(0), *partonPhi(0), *partonMass(0),*partonMatchDR(0),*partonEta(0);
  float yTTbarParton(0), ptTTbarParton(0), mTTbarParton(0);
  vector<int> *partonMatchIdx(0);
  float genEvtWeight(0);
  int nJets;

  trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
  trIN->SetBranchAddress("mTTbarParton" ,&mTTbarParton);
  trIN->SetBranchAddress("yTTbarParton" ,&yTTbarParton);
  trIN->SetBranchAddress("ptTTbarParton",&ptTTbarParton);
  trIN->SetBranchAddress("partonPt"     ,&partonPt);
  trIN->SetBranchAddress("partonEta"    ,&partonEta);
  trIN->SetBranchAddress("partonPhi"    ,&partonPhi);
  trIN->SetBranchAddress("partonMass"    ,&partonMass);
  trIN->SetBranchAddress("partonMatchDR" ,&partonMatchDR);
  trIN->SetBranchAddress("partonMatchIdx",&partonMatchIdx);
  trIN->SetBranchAddress("nJets"          ,&nJets);
  
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

  bool partonCuts, particleCuts;
  
  long NN = trIN->GetEntries();
    //NN = 100000;
  std::cout<<"Entries: "<<NN<<std::endl;
  std::vector<float> xPartonAll(0);
  std::vector<float> xParticleAll(0);
  int decade(0);
  
  for(int iev=1;iev<NN;iev++) 
  {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;
   trIN->GetEntry(iev);
   xPartonAll.clear();
   xParticleAll.clear();

	 int leadingPt =0;
   int genLeadingPt = 0;
   int subleadingPt = 1;
   int genSubleadingPt = 1;

  if((*partonPt)[0] < (*partonPt)[1])
    {
        subleadingPt =0;
        leadingPt = 1;
    }
    xPartonAll.push_back(mTTbarParton);
    xPartonAll.push_back(ptTTbarParton);
    xPartonAll.push_back(yTTbarParton);
    xPartonAll.push_back((*partonPt)[leadingPt]);
    xPartonAll.push_back((*partonPt)[subleadingPt]);
  partonCuts = fabs((*partonEta)[0]) < 2.4 && fabs((*partonEta)[1]) <2.4 && (*partonPt)[0] > 400 && (*partonPt)[1] > 400 && mTTbarParton > 1000;
	
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
  
    particleCuts = fabs((*genjetEta)[0]) < 2.4 && fabs((*genjetEta)[1]) && (*genjetPt)[0] > 400 && (*genjetPt)[1] > 400 && mJJGen > 1000 && nJetsGen >1 &&
      (*genjetMassSoftDrop)[0] > 120 && (*genjetMassSoftDrop)[0] < 220 && (*genjetMassSoftDrop)[1] > 120 && (*genjetMassSoftDrop)[1] < 220;
    
  }

	for(int ivar = 0; ivar < xParticleAll.size(); ivar++)
	{
		if(partonCuts) hParton[ivar]->Fill(xPartonAll[ivar], genEvtWeight);
		if(particleCuts) hParticle[ivar]->Fill(xParticleAll[ivar], genEvtWeight);				
	}//end of ivar fill histograms for-loop


  } //end of tree entries for-loop


  //scale to luminosity and bin width
  for(int ivar=0; ivar<xParticleAll.size(); ivar++) 
  {
    hParton[ivar]->Scale(1/LUMI, "width");
    hParticle[ivar]->Scale(1/LUMI, "width");
  }

  TFile *outf = new TFile("TheoryTemplates.root", "RECREATE");
  outf->cd();
  //write them to file
  for(int ivar=0; ivar<xParticleAll.size(); ivar++) 
  {
    hParton[ivar]->Write(TString::Format("hParton_%s", var[ivar].Data()));
    hParticle[ivar]->Write(TString::Format("hParticle_%s", var[ivar].Data()));
  }

  outf->Close();



}

