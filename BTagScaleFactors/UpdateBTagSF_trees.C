#include "2016/BTagCalibrationStandalone.h"
#include "2016/BTagCalibrationStandalone.cpp"

void UpdateBTagSF_trees(TString oldFileName, TString year = "2017")
{
  /*for(int i=0; i<files.size(); i++)
  {*/
  //TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8_v1.root
    TString eosPath = TString::Format("/eos/cms/store/user/ipapakri/ttbar/MC/Signal/%s/variations/",year.Data());
    std::cout<<"Openning: "<<eosPath+oldFileName<<std::endl;
    TFile *oldFile = TFile::Open(eosPath+oldFileName);


    TString newFileName = oldFileName;
    TString newEosPath = TString::Format("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-%s/variations/",year.Data());
    std::cout<<"Output file is: "<<newEosPath+newFileName<<std::endl;
    TFile *newFile = TFile::Open(newEosPath + newFileName, "RECREATE");

    //for( auto&& keyAsObj: *oldFile->GetListOfKeys()
    for(auto&& keyAsObj: *oldFile->GetListOfKeys())
    {

      auto key = (TKey*)keyAsObj;
      cout<<key->GetName()<<" "<< key->GetClassName() << endl;

      TH1F* triggerNames = (TH1F*) oldFile->Get(TString::Format("%s/TriggerNames",key->GetName()));
      TH1F* cutFlow = (TH1F*) oldFile->Get(TString::Format("%s/CutFlow",key->GetName()));
      TH1F* triggerPass = (TH1F*) oldFile->Get(TString::Format("%s/TriggerPass",key->GetName()));
      TTree* tr = (TTree*) oldFile->Get(TString::Format("%s/events",key->GetName()));
      TString evtCntStr = "eventCounter";
      TString source = "boostedShiftedSrc6Up";
      if(key->GetName() == evtCntStr)
      {
        TH1F* GenEventWeight = (TH1F*) oldFile->Get("eventCounter/GenEventWeight");
        TH1F* pileup = (TH1F*) oldFile->Get("eventCounter/pileup");
        TTree* events = (TTree*) oldFile->Get("eventCounter/events");

        newFile->mkdir("eventCounter");
        newFile->cd("eventCounter");

        GenEventWeight->Write("GenEventWeight");
        pileup->Write("pileup");

        TTree *newTree = events->CloneTree(-1, "fast");
        newTree->Write();
        newFile->cd();
      }
      else if(!(key->GetName() == source))
      {
        newFile->mkdir(key->GetName());
        newFile->cd(key->GetName());

        triggerNames->Write("TriggerNames");
        cutFlow->Write("CutFlow");
        triggerPass->Write("TriggerPass");

        //tr->SetBranchStatus("jetTtagCategory", 0);
        TTree *outputTree = tr->CloneTree(-1, "fast");

        float ht;
        int nJets;
        //TTree *tr = (TTree*) f->Get("boosted/events");
        vector<float> *jetPt(0);
        vector<float> *jetMassSub0(0), *jetMassSub1(0);
        vector<float> *jetFlavorSub0(0), *jetFlavorSub1(0);
        vector<float> *jetPtSub0(0), *jetPtSub1(0);
        vector<float> *jetEtaSub0(0), *jetEtaSub1(0);

        outputTree->SetBranchAddress("jetPt", &jetPt);
        outputTree->SetBranchAddress("jetMassSub0", &jetMassSub0);
        outputTree->SetBranchAddress("jetMassSub1", &jetMassSub1);
        outputTree->SetBranchAddress("jetFlavorSub0"  ,&jetFlavorSub0);
        outputTree->SetBranchAddress("jetFlavorSub1"  ,&jetFlavorSub1);
        outputTree->SetBranchAddress("jetPtSub0"      ,&jetPtSub0);
        outputTree->SetBranchAddress("jetPtSub1"      ,&jetPtSub1);
        outputTree->SetBranchAddress("jetEtaSub0"     ,&jetEtaSub0);
        outputTree->SetBranchAddress("jetEtaSub1"     ,&jetEtaSub1);

        std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
        std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);
        outputTree->SetBranchAddress("nJets", &nJets);
        outputTree->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
        outputTree->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
        outputTree->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
        outputTree->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);

        //TBranch *bpt = outputTree->Branch("jetTtagCategoryMLP", "vector<float>", &jetTtagCategory);
        //create the new Branch
        double ievDeepCSVWgt;
        //trOut->Branch("bTagEvntWeight",&ievDeepCSVWgt, "bTagEvntWeight/D");
        TBranch *bpt = outputTree->Branch("bTagEvntWeight",&ievDeepCSVWgt,"bTagEvntWeight/D");

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

        int decade(0);
        int NN = tr->GetEntries();

        for(int evt = 0; evt<NN; evt++)
        {
          ievDeepCSVWgt = 1;
          double progress = 10.0*evt/(1.0*NN);
          int k = TMath::FloorNint(progress);
          if (k > decade)
            cout<<10*k<<" %"<<endl;
          decade = k;
          outputTree->GetEntry(evt);

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

          bpt->Fill();
        }//end of evt loop

        newFile->cd(key->GetName());
      }//end of else



    }//end of for loop over keys
    newFile->Write();
    oldFile->Close();
    newFile->Close();
  /*}*/
}
