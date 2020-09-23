#include "2018/BTagCalibrationStandalone.h"
#include "2018/BTagCalibrationStandalone.cpp"

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"

typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;


void CalculateBTagSF(TString year = "2016", TString infStr = "")
{

   //btaging scale factors are applied only on MC samples
   //TString *eospath = TString::Format("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-%s/Signal",year.Data());
   TFile *inf = TFile::Open("../../WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8.root", "UPDATE");
   //TFile *inf = TFile::Open(eospath+infStr, "UPDATE");

   TTree *trIN    = (TTree*)inf->Get("boosted/events");
   //cout<<"XSEC: "<<XSEC[f]<<endl;
   //cout<<"weight: "<<weights[f]<<endl;
   int decade(0);
   int NN = trIN->GetEntries();

   int nJets,nLeptons;
   float genEvtWeight;
   vector<float> *jetPt(0),*tau3(0),*tau2(0),*tau1(0);
   vector<float> *jetMassSub0(0), *jetMassSub1(0);
   vector<float> *jetFlavorSub0(0), *jetFlavorSub1(0);
   vector<float> *jetPtSub0(0), *jetPtSub1(0);
   vector<float> *jetEtaSub0(0), *jetEtaSub1(0);
   vector<float> *jetMassSoftDrop(0);

   float mva(0);
   vector<float> *jetTtag(0);
   vector<bool> *bit = new vector<bool>;
   float mTTbarParton(0),mJJ(0), yJJ(0), ptJJ(0), yTTbarParton(0), ptTTbarParton(0);
   int  category(0);
   //matching info
   vector<float> *jetPhi(0), *jetEta(0), *jetY(0);
   vector<int> *partonId(0), *partonMatchIdx(0);

   vector<float> *partonEta(0), *partonPhi(0), *partonMatchDR(0),  *partonPt(0), *partonE(0), *partonMass(0), *deepAK8(0);
   std::vector<int> *addedIndexes = new std::vector<int>(0);
   std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);
   std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
   std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);
   //------- input tree --------------
   trIN->SetBranchAddress("nJets"          ,&nJets);
   trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
   trIN->SetBranchAddress("jetPt"          ,&jetPt);
   trIN->SetBranchAddress("jetEta"         ,&jetEta);
   trIN->SetBranchAddress("jetY"           ,&jetY);
   trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
   trIN->SetBranchAddress("jetTau3"        ,&tau3);
   trIN->SetBranchAddress("jetTau2"        ,&tau2);
   trIN->SetBranchAddress("jetTau1"        ,&tau1);
   trIN->SetBranchAddress("triggerBit"     ,&bit);
   trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
   trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
   trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
   trIN->SetBranchAddress("jetFlavorSub0"  ,&jetFlavorSub0);
   trIN->SetBranchAddress("jetFlavorSub1"  ,&jetFlavorSub1);
   trIN->SetBranchAddress("jetPtSub0"      ,&jetPtSub0);
   trIN->SetBranchAddress("jetPtSub1"      ,&jetPtSub1);
   trIN->SetBranchAddress("jetEtaSub0"     ,&jetEtaSub0);
   trIN->SetBranchAddress("jetEtaSub1"     ,&jetEtaSub1);
   trIN->SetBranchAddress("mJJ"            ,&mJJ);
   trIN->SetBranchAddress("yJJ"            ,&yJJ);
   trIN->SetBranchAddress("ptJJ"           ,&ptJJ);
   trIN->SetBranchAddress("jetBtagSub0"    ,&jetBtagSub0);
   //trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
   trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
   trIN->SetBranchAddress("mva"          ,&mva);
   trIN->SetBranchAddress("category"       ,&category);
   trIN->SetBranchAddress("deepAK8Tagger"  ,&deepAK8);
   trIN->SetBranchAddress("jetTtagCategory",&jetTtag);

   //deepCSV
   trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
   trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
   trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
   trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);

   //create the new Branch
   double ievDeepCSVWgt;
   //trOut->Branch("bTagEvntWeight",&ievDeepCSVWgt, "bTagEvntWeight/D");
   TBranch *trIN_br = trIN->Branch("bTagEvntWeight",&ievDeepCSVWgt,"bTagEvntWeight/D");
   //=================================================================================================================================
   //Load the .csv file containing the SFs using BTagCalibration Standalone tool
   //=================================================================================================================================

   std::cout << "===> Loading the input .csv SF file..." << std::endl;

   std::string inputCSVfile;
   if(year.EqualTo("2016"))inputCSVfile = "2016/DeepCSV_2016LegacySF_V1.csv";
   else if(year.EqualTo("2017")) inputCSVfile= "2017/DeepCSV_94XSF_V5_B_F.csv";
   else inputCSVfile= "2018/DeepCSV_102XSF_V2.csv";
   //std::string inputCSVfile = "2017/subjet_DeepCSV_94XSF_V4_B_F.csv";
   std::string measType = "comb";
   std::string sysType = "central";

   BTagCalibration calib("deepCSV", inputCSVfile);
   BTagCalibrationReader reader(BTagEntry::OP_MEDIUM, sysType) ;

   reader.load(calib, BTagEntry::FLAV_B, measType);
   reader.load(calib, BTagEntry::FLAV_C, measType);
   reader.load(calib, BTagEntry::FLAV_UDSG, measType);

   std::cout << "\tInput CSV weight file = " << inputCSVfile << "; measurementType = " << measType << "; sysType = " << sysType << std::endl;
   cout<<"Reading "<<NN<<" entries"<<endl;

   for(int iev=0;iev<NN;iev++)
   {
     ievDeepCSVWgt = 1;
     double progress = 10.0*iev/(1.0*NN);
     int k = TMath::FloorNint(progress);
     if (k > decade)
       cout<<10*k<<" %"<<endl;
     decade = k;
     trIN->GetEntry(iev);

     //how to create a new Branch and fill it
     //I have to calculate the Btagging SF in some way!

     //loop on all the jets:
     //I need jetPtSub0,1 and jetEtaSub0,1
    for(int ijet=0; ijet<nJets; ijet++)
    {
      double ievDeepCSVWgtJet =1;
      //jetFlavorSub0,1 [ijet]
      float pT[2], eta[2];
      int flavor[2];
      float dCSVScore[2];
      //pT[0] is the pt of the 1st subjet
      pT[0] = (*jetPtSub0)[ijet];
      pT[1]= (*jetPtSub1)[ijet];
      eta[0] = fabs((*jetEtaSub0)[ijet]);
      eta[1] = fabs((*jetEtaSub1)[ijet]);


      dCSVScore[0] = (*jetBtagSub0DCSVbb)[ijet] + (*jetBtagSub0DCSVbbb)[ijet];
      dCSVScore[1] = (*jetBtagSub1DCSVbb)[ijet] + (*jetBtagSub1DCSVbbb)[ijet];


      flavor[0] = fabs((*jetFlavorSub0)[ijet]);
      flavor[1] = fabs((*jetFlavorSub1)[ijet]);

      for(int isub =0; isub<2; isub++)
      {
        if(flavor[isub] == 5)
        {
          double ievDeepCSVWgtTemp = reader.eval(BTagEntry::FLAV_B, eta[isub], pT[isub], dCSVScore[isub]);
          if(ievDeepCSVWgtTemp!=0) ievDeepCSVWgtJet *=ievDeepCSVWgtTemp;
        }
      }//end of subjet loop
      ievDeepCSVWgt *= ievDeepCSVWgtJet;

    }//end of jet loop
    //cout<<"event: "<<iev<<" has weight: "<<ievDeepCSVWgt<<endl;
    //now fill the branch with the ievDeepCSVWgt
    trIN_br->Fill();

  }//-----end of iev loop-----
  inf->cd("boosted");
  //trIN->Write("boosted");
  trIN->Write();

/*
   // Using BTagCalibrationReader
   double csvWgtHF = 1., csvWgtLF = 1., csvWgtC = 1.;
   double csvWgtTotal = 1.;

   if ( insample>=0 ) {		//Do it for MC only, not data
   for( int iJet=0; iJet<int(jet_vect_TLV.size()); iJet++ ){
     TLorentzVector myJet = jet_vect_TLV[iJet];

     float csv = jet_CSV[iJet];
     float jetPt = myJet.Pt();
     float jetAbsEta = fabs(myJet.Eta());
     int flavor = jet_flavour[iJet];

     if (abs(flavor) == 5 ){    //HF
       double iCSVWgtHF = reader.eval(BTagEntry::FLAV_B, jetAbsEta, jetPt, csv);
       if( iCSVWgtHF!=0 ) csvWgtHF *= iCSVWgtHF;
     }
     else if( abs(flavor) == 4 ){  //C
       double iCSVWgtC = reader.eval(BTagEntry::FLAV_C, jetAbsEta, jetPt, csv);
       if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
     }
     else { //LF
       double iCSVWgtLF = reader.eval(BTagEntry::FLAV_UDSG, jetAbsEta, jetPt, csv);
       if( iCSVWgtLF!=0 ) csvWgtLF *= iCSVWgtLF;
     }
    }
       csvWgtTotal = csvWgtHF * csvWgtC * csvWgtLF;
   }

   void upd()
   {
     TFile *f = new TFile("hs.root","update");
     TTree *T = (TTree*)f->Get("ntuple");
     float px,py;
     float pt;
     TBranch *bpt = T->Branch("pt",&pt,"pt/F");
     T->SetBranchAddress("px",&px);
     T->SetBranchAddress("py",&py);
      Long64_t nentries = T->GetEntries();
      for (Long64_t i=0;i<nentries;i++)
      {
        T->GetEntry(i);
        pt = TMath::Sqrt(px*px+py*py);
        bpt->Fill();
      }
      T->Print();
      T->Write();
      delete f;
    }
    */

}
