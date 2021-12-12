#include "TemplateConstants.h"
#include "../CMS_plots/tdrstyle.C"
#include "../CMS_plots/CMS_lumi.C"

//We split the uncertaintied in 5 groups and we create a histogram for each group
//1. Statistical which is basically the error bars of the unfolded histograms
//2. JES+JER+Pileup which is the variation (difference from the nominal restuls)
//   we get from these source
//3. Flavor tagging
//4. Parton Shower
//5. Hard scattering

//For every group we create a histogram that will containt the combined uncertainty
//that resulsts from the addition of the uncertainties from this gorup.
//The combination of the uncertainties comes by adding the in qudrature. Although
//we split the uncertainties in up and down depending on wheter they result in a
//upper or lower value for the XS in that bin we in the end present one combined
//which is given by 0.5 * 100 * (eUp + eDown)
std::vector<TString> listFiles(const char *dirname="", const char *var="", const char *ext=".root")
{
  std::vector<TString> list_of_files;
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
     TSystemFile *file;
     TString fname;
     TIter next(files);
     while ((file=(TSystemFile*)next())) {
       fname = file->GetName();
       if (!file->IsDirectory() && fname.EndsWith(ext) && fname.Contains(var)) {
         //cout << fname.Data() << endl;
         list_of_files.push_back(fname.Data());
       }
     }
   }
   return list_of_files;
}

void SystematicsUnfolding_levels(TString isParton = "Particle")
{
  gStyle->SetOptStat(0);
  initFilesMapping();
  std::vector<TString> dirs = {"Nominal", "bTagVariation", "PDFWeights", "ScaleWeights"};
  std::vector<TString> groups = {"Stat. Uncertainty", "Flavor Tagging", "PDF", "Scale"};
  std::vector<int> groupColors = {kBlack, kRed, kGreen, kBlue, kOrange};

  const int NVAR = 10;
  TString variable[NVAR]   = {"jetPt0","jetPt1","yJJ", "ptJJ", "mJJ", "jetY0", "jetY1", "chi", "cosTheta_0", "cosTheta_1"};
  TString variableParton[NVAR] = {"partonPt0", "partonPt1", "yTTbarParton", "ptTTbarParton", "mTTbarParton", "partonY0", "partonY1","chiParton", "cosThetaParton_0", "cosThetaParton_1"};
  TString variableGen[NVAR] = {"genjetPt0", "genjetPt1", "yJJGen", "ptJJGen", "mJJGen", "genjetY0", "genjetY1", "chiParticle", "cosThetaParticle_0", "cosThetaParticle_1"};

  //TString baseInputDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/";
  TString baseInputDir = "/Users/georgebakas/Documents/HEP-NTUA_ul/VariationHandling_Theory_amc@NLO/2018";
  TString outputDirectory = TString::Format("/SystematicsBreakdown/");
  //CheckAndCreateDirectory(outputDirectory);
  TString histo_name = "hUnfoldFinal_";
  //TString histo_name = "hSigInit_";
  
  

  TString fileName = TString::Format("%s/Nominal/Histograms_TTJets.root", baseInputDir.Data());
  TFile *fNominal = TFile::Open(fileName);

  for (int i = 0; i<NVAR; i++)
  {
    TString histoName;
    if(isParton.EqualTo("Parton"))
    {   
        histoName = "hParton_";
        vars[i] = variableParton[i];
    }
    else
    {
        histoName = "hParticle_";
        vars[i] = variableGen[i];
    }

    
    std::cout << "Variable: " << vars[i]<< std::endl;
    std::vector<TH1F *> groupHistogramsUp;
    std::vector<TH1F *> groupHistogramsDown;
    std::vector<TH1F *> groupHistogramsSym;

    TH1F *hNominal = (TH1F *)fNominal->Get(TString::Format("h%s_%s", isParton.Data() ,vars[i].Data()));

    //initialize group histograms
    for (int group = 0; group < groups.size(); group++)
    {
      TH1F *hSystematicsUp = (TH1F *)hNominal->Clone(TString::Format("%s %s Up", groups[group].Data(),
                                                                            vars[i].Data()));
      TH1F *hSystematicsDown = (TH1F *)hNominal->Clone(TString::Format("%s %s Down", groups[group].Data(),
                                                                              vars[i].Data()));
      TH1F *hSystematicsSym = (TH1F *)hNominal->Clone(TString::Format("%s %s Sym", groups[group].Data(),
                                                                              vars[i].Data()));
      hSystematicsUp->Reset();
      hSystematicsDown->Reset();
      hSystematicsSym->Reset();
      hSystematicsUp->SetLineColor(groupColors[group-1]);
      hSystematicsDown->SetLineColor(groupColors[group-1]);
      hSystematicsSym->SetLineColor(groupColors[group-1]);

      hSystematicsUp->SetLineWidth(3);
      hSystematicsDown->SetLineWidth(3);
      hSystematicsSym->SetLineWidth(3);

      groupHistogramsUp.push_back(hSystematicsUp);
      groupHistogramsDown.push_back(hSystematicsDown);
      groupHistogramsSym.push_back(hSystematicsSym);
    }

    //Start of nominal handling
    for (int bin = 0; bin < groupHistogramsUp[0]->GetNbinsX(); bin++)
    {
      double nominalValue = hNominal->GetBinContent(bin + 1);
      double nominalError = hNominal->GetBinError(bin + 1);
      if (nominalValue != 0)
      {
        groupHistogramsUp[0]->SetBinContent(bin + 1, (nominalError * nominalError) / (nominalValue * nominalValue));
        groupHistogramsDown[0]->SetBinContent(bin + 1, (nominalError * nominalError) / (nominalValue * nominalValue));
      }
      //std::cout << "nominalValue: " << nominalValue << " nominalError: " << nominalError << " " << (nominalError * nominalError) / (nominalValue * nominalValue) << std::endl;
    }
    //end of nominal handling

    //start variation handling

    for (int j = 0; j < dirs.size(); j++)
    {
      //std::cout << "Variation: " << dirs[j] << std::endl;
      TString variation = dirs[j];
      std::vector<TString> variationFiles = listFiles(TString::Format("2018/%s/", dirs[j].Data()), "");

      int group = -1;
      if (variation.Contains("Nominal"))
        group = 0;
      if (variation.Contains("Scale"))
        group = 1;
      if (variation.Contains("bTagVariation"))
        group = 2;
      if (variation.Contains("PDF"))
        group = 3;


      cout<<variation<<" group: "<<group <<endl;

      if (group == -1)
        continue;

      for(int jfile=0; jfile<variationFiles.size(); jfile++)
      {

        cout<<variation<<endl;
        cout<<variationFiles[jfile].Data()<<endl;

        if (variationFiles[jfile].Contains("nom0") || variationFiles[jfile].Contains("nom1")) continue;
        fileName = TString::Format("2018/%s/%s",
                                  variation.Data(),
                                  variationFiles[jfile].Data());
        
        if (!fileName.Contains("Def") && variation.Contains("PS"))
          continue;
        

        TFile *f = TFile::Open(fileName);

         //do this to get coherent name
        TString hname = variationFiles[jfile];
        hname.ReplaceAll("Histograms_TTJets_", "");
        hname.ReplaceAll(".root", "");
        
        TH1F *hVariation;
        if (variation.Contains("scale_9")) continue;
        if (variation.Contains("bTag") || variation.EqualTo("Nominal"))
        {   
            hVariation = (TH1F *)f->Get(TString::Format("%s%s",
                                        histoName.Data(),
                                        vars[i].Data()));
            cout<<TString::Format("%s%s",histoName.Data(), 
                                        vars[i].Data())<<endl;
        }
        else
        {   
            hVariation = (TH1F *)f->Get(TString::Format("%s%s_%s",
                                        histoName.Data(),
                                        vars[i].Data(), 
                                        hname.Data()));
            cout<<TString::Format("%s%s_%s",histoName.Data(), 
                                        vars[i].Data(), 
                                        hname.Data())<<endl;
            hVariation->Scale(2);
        }

        for (int bin = 0; bin < groupHistogramsUp[group]->GetNbinsX(); bin++)
        {
          double nominalValue = hNominal->GetBinContent(bin + 1);
          double nominalError = hNominal->GetBinError(bin + 1);
          double variationValue = hVariation->GetBinContent(bin + 1);
          double valuePull = (variationValue - nominalValue) / nominalValue;
          double variationErrorUp = groupHistogramsUp[group]->GetBinContent(bin + 1);
          double variationErrorDown = groupHistogramsDown[group]->GetBinContent(bin + 1);

          //std::cout << nominalValue << " " << nominalError << endl;
          //std::cout << " " << variationValue << " " << valuePull << endl;
          std::cout << fabs(variationValue - nominalValue) /nominalValue << endl;
          if (fabs(variationValue - nominalValue) /nominalValue > 0.7 ) 
            std::cout << fileName <<endl;
          
          //std::cout << " " << variationErrorUp << " " << variationErrorDown;
          if (nominalValue != 0)
          {
            if (valuePull >= 0)
            {
              //std::cout << " Variation up:";
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1);
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
              groupHistogramsUp[group]->SetBinContent(bin + 1, variationErrorUp + (valuePull * valuePull));
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1) <<endl;
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
            }
            else
            {
              //std::cout << " Variation down:";
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1);
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
              groupHistogramsDown[group]->SetBinContent(bin + 1, variationErrorDown + (valuePull * valuePull));
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1) <<endl;
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
            }
          }
          //std::cout << endl;
        }
        //std::cout << endl;
        f->Close();
      }//end of for loop (over specific variation files)
    } //end of for loop (over all variation types)


    for (int group = 0; group < groups.size(); group++)
    {
      std::cout << "Group: " << groups[group] << std::endl;
      for (int bin = 0; bin < groupHistogramsUp[group]->GetNbinsX(); bin++)
      {
        double errorUp = TMath::Sqrt(groupHistogramsUp[group]->GetBinContent(bin + 1));
        double errorDown = TMath::Sqrt(groupHistogramsDown[group]->GetBinContent(bin + 1));
        std::cout <<"bin "<<bin<<" " << errorUp << " " << errorDown;
        groupHistogramsSym[group]->SetBinContent(bin + 1, 0.5 * 100 * (errorUp + errorDown));
        std::cout << " " << groupHistogramsSym[group]->GetBinContent(bin + 1) << std::endl;
      }
    }
    TCanvas *c1 = new TCanvas(vars[i], vars[i], 800, 600);

    groupHistogramsSym[0]->SetFillColor(kGray);
    //groupHistogramsSym[0]->GetYaxis()->SetRangeUser(0, 120);
    groupHistogramsSym[0]->SetTitle("");
    groupHistogramsSym[0]->SetLineWidth(0);
    groupHistogramsSym[0]->GetXaxis()->SetTitle(vars[i]);
    groupHistogramsSym[0]->GetXaxis()->SetLabelSize(0.035);
    std::cout << groupHistogramsSym[0]->GetYaxis()->GetLabelSize() << std::endl;
    groupHistogramsSym[0]->Draw("hist");
    groupHistogramsSym[0]->GetYaxis()->SetRangeUser(0, 120);
    groupHistogramsSym[0]->GetYaxis()->SetTitle("Relative Uncertainty [%]");
    //groupHistogramsSym[0]->Draw();
    if (groupHistogramsSym[0]->GetMaximum() < 120)
    {
      groupHistogramsSym[0]->GetYaxis()->SetRangeUser(0, 120);
    }
    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9);
    leg->AddEntry(groupHistogramsSym[0], groups[0], "f");

    for (int group = 1; group < groups.size(); group++)
    {
      groupHistogramsSym[group]->Draw("hist same");
      leg->AddEntry(groupHistogramsSym[group], groups[group], "l");
    }

    leg->Draw();

    lumi_13TeV = TString::Format("%0.1f fb^{-1}", luminosity["luminosityAll"]/1000);
    writeExtraText = true;
    CMS_lumi(c1, 4, 10);
    c1->SaveAs(TString::Format("%s%s/Systematics%s%s_%s.png", 
                        baseInputDir.Data(), outputDirectory.Data(), 
                        histo_name.Data(), isParton.Data(), vars[i].Data()));
  }
}
