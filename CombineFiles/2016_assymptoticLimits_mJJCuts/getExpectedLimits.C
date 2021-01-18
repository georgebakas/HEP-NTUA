

void getExpectedLimits(int Z_mass, int Z_mass_width)
{
  float massCut[] = {1000,1200,1400,1600,1800};
  float massCutErrors[] = {0,0,0,0,0};
  float expected_limits[5];
  float expected_limits_errors[5];


  for (int im = 0; im < 5; im++)
  {
    cout<<"Mass Cut "<<massCut[im]<<endl;
    TFile *inf = TFile::Open(TString::Format("higgsCombine.mZ_%d_%dCut%d.AsymptoticLimits.mH120.root",
                            Z_mass, Z_mass_width, (int)massCut[im]));

    TTree *tr = (TTree*)inf->Get("limit");
    double limit, limitErr;
    tr->SetBranchAddress("limit",&limit);
    tr->SetBranchAddress("limitErr",&limitErr);
    //expected is median--> ilim == 2,
    //5 is observed
    for (int ilim =0; ilim < tr->GetEntries(); ilim++ )
    {
      tr->GetEntry(ilim);
      //cout<<limit<<"±"<<limitErr<<endl;
      if (ilim==2)
      {
        cout<<"Expected limit is: "<<limit<<endl; //error only in observed limit
        expected_limits[im] = limit;
        expected_limits_errors[im] = limitErr;
        //cout<<limit<<"±"<<limitErr<<endl;
      }
    } // end of ilim for loop

  } //end of im for loop

  //TGraph
  int n = 5;
  TCanvas *can = new TCanvas("can_", "can_", 800, 600);
  TGraph* gr = new TGraphErrors(n,massCut,expected_limits, massCutErrors, expected_limits_errors);
  //gr->SetNameTitle(TString::Format("Expected Limits for Z' (M:%d, W:%d)", Z_mass, Z_mass_width));
  gr->SetName(TString::Format("Expected Limits for Z' (M:%d, W:%d)", Z_mass, Z_mass_width));
  gr->SetTitle(TString::Format("Expected Limits for Z' (M:%d, W:%d)", Z_mass, Z_mass_width));
  gr->GetXaxis()->SetTitle("mJJCut (GeV)");
  gr->GetYaxis()->SetTitle("Exp. Limit");
  gr->Draw("A*");

  can->Print(TString::Format("ExpLim_Z_(M:%dW:%d).pdf", Z_mass, Z_mass_width), "pdf");

}
