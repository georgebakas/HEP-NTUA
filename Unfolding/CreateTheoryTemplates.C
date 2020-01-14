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
#include "TGraphErrors.h"
#include "TUnfold.h"

using namespace std;

using std::cin;
using std::cout;
using std::endl;


void CreateTheoryTemplates(TString inYear = "2016", bool isParton = true)
{
  year = inYear;
  std::vector< std::vector <Float_t> > const BND_gen = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}, //mjj
                                                        {0,60,150,300,450,600,750,950,1100,1300}, //ptjj
                                                        {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
                                                        {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt0     
                                                        {400,450,500,570,650,750,850,950,1100,1300,1500}}; //jetPt1
                                                        //{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
                                                        //{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1

  float LUMI = luminosity[year];

  TFile *inf = TFile::Open("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Signal/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");


  TTree *trIn = (TTree*)inf->Get("boosted/events");
  //parton
  std::vector<float> *partonPt(0), *partonPhi(0), *partonMass(0),*partonMatchDR(0),*partonEta(0);
  float yTTbarParton(0), ptTTbarParton(0), mTTbarParton(0);
  vector<int> *partonMatchIdx(0);

  trIN->SetBranchAddress("mTTbarParton" ,&mTTbarParton);
  trIN->SetBranchAddress("yTTbarParton" ,&yTTbarParton);
  trIN->SetBranchAddress("ptTTbarParton",&ptTTbarParton);
  trIN->SetBranchAddress("partonPt"     ,&partonPt);
  trIN->SetBranchAddress("partonEta"    ,&partonEta);
  trIN->SetBranchAddress("partonPhi"    ,&partonPhi);
  trIN->SetBranchAddress("partonMass"    ,&partonMass);
  trIN->SetBranchAddress("partonMatchDR" ,&partonMatchDR);
  trIN->SetBranchAddress("partonMatchIdx",&partonMatchIdx);
  
  //particle
  std::vector<float> *genjetPt(0), *genjetY(0), *genjetEta(0), *genjetMassSoftDrop(0);
  int nJetsGen(0);
  float `(0), ptJJGen(0), yJJGen(0);

  trIN->SetBranchAddress("nJetsGen" ,&nJetsGen);
  trIN->SetBranchAddress("mJJGen"   ,&mJJGen);
  trIN->SetBranchAddress("ptJJGen"  ,&ptJJGen);
  trIN->SetBranchAddress("yJJGen"   ,&yJJGen);
  trIN->SetBranchAddress("genjetPt" ,&genjetPt);
  trIN->SetBranchAddress("genjetEta",&genjetEta);
  trIN->SetBranchAddress("genjetY"  ,&genjetY);
  trIN->SetBranchAddress("genjetMassSoftDrop", &genjetMassSoftDrop);  

  bool partonCuts, particleCuts;


   partonCuts = fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1]) <2.4 && (*partonPt_)[0] > 400 && (*partonPt_)[1] > 400 && mTTbarParton > 1000;
   particleCuts = fabs((*genjetEta)[0]) < 2.4 && fabs((*genjetEta)[1]) && (*genjetPt)[0] > 400 && (*genjetPt)[1] > 400 && mJJGen > 1000 && nJetsGen >1 &&
	  				 (*genjetMassSoftDrop)[0] > 120 && (*genjetMassSoftDrop)[0] < 220 && (*genjetMassSoftDrop)[1] > 120 && (*genjetMassSoftDrop)[1] < 220;


 

}

