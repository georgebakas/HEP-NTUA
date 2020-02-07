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
   		XSEC.push_back(377.96);
   		XSEC.push_back(365.34);
      XSEC.push_back(365.34);
   		XSEC.push_back(88.29);
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
	  	fileNames.push_back("TTHadronic_1");
	  	fileNames.push_back("TTSemiLeptonic_0");
	  	fileNames.push_back("TTSemiLeptonic_1");
      fileNames.push_back("TTTo2L2Nu_0");
      fileNames.push_back("TTTo2L2Nu_1");

	  	histoNames.push_back("Signal_histo_TTHadronic_0");
	  	histoNames.push_back("Signal_histo_TTHadronic_1");
	  	histoNames.push_back("Signal_histo_TTSemiLeptonic_0");
	  	histoNames.push_back("Signal_histo_TTSemiLeptonic_1");
      histoNames.push_back("Signal_histo_TTTo2L2Nu_0");
      histoNames.push_back("Signal_histo_TTTo2L2Nu_1");

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
                                                        {400,450,500,570,650,750,850,950,1100,1300,1500}}; //jetPt1
                                                        //{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
                                                        //{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1
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

  for(int ivar =0; ivar<NVAR-2; ivar++)
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
  

  //cout<<"Working in file: "<<files[year.Data()][fileNames[f].Data()]<<endl;
  TFile *file = TFile::Open(eospath[year.Data()]+files[year.Data()][fileNames[f].Data()]);
  TTree *trIN = (TTree*)file->Get("boosted/events");
  cout<<"working in file:"<<eospath[year.Data()]+files[year.Data()][fileNames[f].Data()]<<endl;
  //parton
  std::vector<float> *partonPt(0), *partonPhi(0), *partonMass(0),*partonMatchDR(0),*partonEta(0);
  float yTTbarParton(0), ptTTbarParton(0), mTTbarParton(0);
  vector<int> *partonMatchIdx(0);
  float genEvtWeight(0);
  int nJets;
  
  float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  float weight = XSEC[f]/norm;
  weights.push_back(weight);

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
		if(partonCuts) hParton[f][ivar]->Fill(xPartonAll[ivar], genEvtWeight);
		if(particleCuts) hParticle[f][ivar]->Fill(xParticleAll[ivar], genEvtWeight);				
	}//end of ivar fill histograms for-loop


  } //end of tree entries for-loop

  }// end of files loop


  //for all vars in all files, scale
  for(int ivar=0; ivar<NVAR-2; ivar++) 
  {
	  hParton[0][ivar]->Scale(weights[0]*LUMI);
	  hParticle[0][ivar]->Scale(weights[0]*LUMI);
	  for(int f=1; f<fileNames.size(); f++) 
	  {
	    hParton[f][ivar]->Scale(weights[f]*LUMI);
		hParton[0][ivar]->Add(hParton[f][ivar]);

		hParticle[f][ivar]->Scale(weights[f]*LUMI);
		hParticle[0][ivar]->Add(hParticle[f][ivar]);
	  }
  }


  //scale to luminosity and bin width
  for(int ivar=0; ivar<NVAR-2; ivar++) 
  {
    hParton[0][ivar]->Scale(1/LUMI, "width");
    hParticle[0][ivar]->Scale(1/LUMI, "width");
  }

  TFile *outf = new TFile(TString::Format("%s/TheoryTemplates.root", year.Data()), "RECREATE");
  outf->cd();
  //write them to file
  for(int ivar=0; ivar<NVAR-2; ivar++) 
  {
    hParton[0][ivar]->Write(TString::Format("hParton_%s", var[ivar].Data()));
    hParticle[0][ivar]->Write(TString::Format("hParticle_%s", var[ivar].Data()));
  }

  outf->Close();
  fileNames.clear();
  histoNames.clear();
  XSEC.clear();
}

