#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1.h"
#include "TCanvas.h"

#include "TemplateConstants.h"

using std::cin;
using std::cout;
using std::endl;
void Limit_ZprimeMass_Width(TString year="2016", int z_mass=2000, float z_width=0.01, float xsec=0.0);


void AsymptoticLimitsPlots_mJJCuts(TString year="2016")
{
  std::vector< std::vector <Float_t> > const cross_section = {
                      {1.730e+00, 9.095e-01, 5.002e-01, 2.833e-01, 1.662e-01, 4.749e-02, 1.494e-02, 5.105e-03, 1.900e-03, 7.613e-04}, // width 1%
                      {0.01842, 0.00563, 0.002006, 0.0008262, 0.0003791}, // width 10%
                      {0.006504, 0.002282, 0.0009163, 0.0004251, 0.0002185}}; // width 30%
  const int number_of_masses = 10;
  int masses[] = {1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500};
  float widths[] = {.01};//, .1, .3};
  TFile *infZprime;
  TFile *outf_Zprime;
  TH1F *hZprime;

  for (int imass = 0; imass<number_of_masses; imass++)
  {
    if (masses[imass] == 3500 && year.EqualTo("2016")) continue;
    for(int iw = 0; iw<1; iw++)
    {
      //if(imass==1 && iw ==2) continue;
      float xsec = cross_section[iw][imass];
      Limit_ZprimeMass_Width(year, masses[imass], widths[iw], xsec);
    }
  }
}

void Limit_ZprimeMass_Width(TString year="2016", int z_mass=2000, float z_width=0.01, float xsec=0.0)
{
  initFilesMapping(false);
  int n = 6;
  float limit_values[6], limit_values_err[6];
  float mJJCuts[] = {1000,1200,1400,1600,1800,2000};
  float mJJCutErrors[] = {0,0,0,0,0,0};
  float width = z_mass*z_width;

  for (int icut=0; icut<n; icut++)
  {
    int mJJCut = (int)mJJCuts[icut];
    //---------------------------------- START reading files --------------------------------------------
    // GET limits from root file 
    /*
    observed.SetPoint(  i,    values[i], limit[5] ) # observed
    yellow.SetPoint(    i,    values[i], limit[4] ) # + 2 sigma
    green.SetPoint(     i,    values[i], limit[3] ) # + 1 sigma
    median.SetPoint(    i,    values[i], limit[2] ) # median
    green.SetPoint(  2*N-1-i, values[i], limit[1] ) # - 1 sigma
    yellow.SetPoint( 2*N-1-i, values[i], limit[0] ) # - 2 sigma
    
   if width == 1:
        #here you have 1200_12, 1400_14, 1600_16, 1800_18, 2000_20,
        #              2500_25, 3000_30, 3500_35, 4000_40, 4500_45
        cross_section = [1.730e+00, 9.095e-01, 5.002e-01, 2.833e-01, 1.662e-01,
                        4.749e-02, 1.494e-02, 5.105e-03, 1.900e-03, 7.613e-04]
    elif width == 10:
        #here we got 2000_200, 2500_250, 3000_300, 3500_350, 4000_400
        cross_section = [0.01842, 0.00563, 0.002006, 0.0008262, 0.0003791]
    else:
        #30% we got 2000_600, 2500_750, 3000_900, 3500_1050, 4000_1200
        cross_section = [0.006504, 0.002282, 0.0009163, 0.0004251, 0.0002185]
    */
    
  
    TFile *inf = TFile::Open(TString::Format("%s_test/higgsCombine.mZ_%d_%dCut%d.AsymptoticLimits.mH120.root", year.Data(), 
                              z_mass, (int)width, mJJCut));

    TTree *tree = (TTree*)inf->Get("limit");
    Double_t limit, limit_err;
    tree->SetBranchAddress("limit", &limit);
    tree->SetBranchAddress("limitErr", &limit_err);

    tree->GetEntry(2);
    float expected_ = limit * xsec;
    float expected_error = limit_err * xsec;
    
    cout<<z_mass<<" "<<z_width<<endl;
    cout<<mJJCut<<endl;
    cout<<expected_<<" +/-" <<expected_error<<endl; 

    //---------------------------------- END of file reading --------------------------------------------;

    limit_values[icut] = expected_;
    limit_values_err[icut] = expected_error;
  }

  //now draw a tgraph

  TCanvas *c1 = new TCanvas(TString::Format("Asymptotic_%d_%d", z_mass, (int)width),TString::Format("Asymptotic_%d_%d", z_mass, (int)width),800,600);
  TGraphErrors* gr = new TGraphErrors(n,mJJCuts,limit_values, limit_values_err, mJJCutErrors);
  gr->SetTitle(TString::Format("Asymptotic M%d W%d", z_mass, (int)width));
  gr->GetYaxis()->SetTitleOffset(1.56);
  gr->GetXaxis()->SetTitle("mJJCut (GeV)");
  gr->GetYaxis()->SetTitle("Expected Limit");
  gr->GetYaxis()->SetRangeUser(0, 4);
  cout<<"Maximum at: "<<gr->GetMaximum()<<endl;
  gr->Draw("A*");
  c1->Print(TString::Format("%s_test/AsymptoticPlots/Asymptotic_%d_%d.pdf", year.Data(), z_mass, (int)width), "pdf");
  c1->Print(TString::Format("%s_test/AsymptoticPlots/Asymptotic_%d_%d.png", year.Data(), z_mass, (int)width), "png");

}
