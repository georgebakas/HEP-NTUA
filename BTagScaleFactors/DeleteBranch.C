void DeleteBranch(TString oldFileName)
{
  	
    //TString weightsFile = "/afs/cern.ch/work/i/ipapakri/private/analysis/TopAnalysis/MVA_Trainings/Training_2016/BoostedMVA/weights/boosted_MVA_MLPCat.weights.xml";
    std::cout<<"Openning: "<<oldFileName<<std::endl;
    TFile *oldFile = TFile::Open(oldFileName+".root", "update");

    TTree* tr = (TTree*) oldFile->Get("boosted/events");

    TBranch *b_temp = tr->GetBranch("jetTtagCategory");
    tr->GetListOfBranches()->Remove(b_temp);	
    TLeaf *l = tr->GetLeaf("jetTtagCategory");
    tr->GetListOfLeaves()->Remove(l);
    cout<<"deleted jetTtagCategory branch!"<<endl;
    //tr->update();
    oldFile->Write();
    oldFile->Close();
}