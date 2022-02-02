#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include <TParameter.h>

using namespace std;

#include <map>
#include <string>
#include <iostream>


void PlotMassDiff(TString year)
{
    gStyle->SetOptStat(0);

    TFile *inf = TFile::Open(TString::Format("%s/MassDifferences.root", year.Data()));

    // Setup structure for data and desired labels
    struct data_struct {
        int xVal;
        float yVal;
        float yError;
        TString xLabel;
    };
    map<int,data_struct> data_map;
    
    // get all objects within the file 
    TObject *obj;
    TKey *key;
    TIter next(inf->GetListOfKeys());
    int i =0;


    while ((key = (TKey *) next())) {

        obj = inf->Get(key->GetName()); // copy object to memory
        // do something with obj
        printf(" found object:%s\n",key->GetName());
        TString object_name = key->GetName();
        if (!(object_name.Contains("error")))
        {
            TParameter<float> *p   = (TParameter<float>*)inf->Get(key->GetName()); 
            TParameter<float> *p_e = (TParameter<float>*)inf->Get(TString::Format("%s_error", key->GetName())); 
            cout<<p->GetVal()<<" ± "<<p_e->GetVal()<<endl;
            float val = p->GetVal();
            float error = p_e->GetVal();
            TString name = key->GetName();
            name.ReplaceAll("MassFitResults_SignalTemplates_", "");
            name.ReplaceAll(".root", "");
            if (name.Contains("nom") || name.Length() ==0)
                continue;
            // for pdf variations
            if (name.IsDigit())
                name = TString::Format("pdf_%s", name.Data());
            data_struct tmp;  
            tmp = {i, val, error, name};
            
            data_map.insert(pair<int,data_struct>(i, tmp));
            // fill the graph
            i++;
        }
    }
    
    // Print the data and labels
    for (map<int,data_struct>::iterator it = data_map.begin(); it != data_map.end(); ++it) {
        cout << it->first << " => " << it->second.xVal << " " << it->second.xVal << " " << it->second.yVal << " " << it->second.xLabel << endl;
    }

    // Instantiate graph and set data
    TGraphErrors *g_original = new TGraphErrors();
    g_original->SetMarkerStyle(20);
    for (map<int,data_struct>::iterator it = data_map.begin(); it != data_map.end(); ++it) {
        g_original->SetPoint(it->first, it->second.xVal, it->second.yVal);
        g_original->SetPointError(it->first, 0, it->second.yError);
    } 

    // Clone graph, get axis object from clone, and change axis bin labels
    TGraphErrors *g_altered = (TGraphErrors*)g_original->Clone();
    TAxis *a = g_altered->GetXaxis();
    int bin;
    for (map<int,data_struct>::iterator it = data_map.begin(); it != data_map.end(); ++it) {
        bin = a->FindBin(it->second.xVal);
        cout<<it->second.xVal<< " "<<bin<<endl;
        a->SetBinLabel(bin, it->second.xLabel);
    }

    // Draw altered plot
    TCanvas *c_altered = new TCanvas("c1", "c1", 1200, 600);
    c_altered->cd();
    g_altered->GetYaxis()->SetRangeUser(0.8, 1.2);
    g_altered->GetXaxis()->SetRangeUser(0, i);
    g_altered->GetXaxis()->SetLabelOffset(0.005);
    g_altered->GetXaxis()->SetLabelSize(0.03);
    g_altered->GetXaxis()->LabelsOption("v");
    g_altered->Draw("AP");

    c_altered -> Print(TString::Format("%s/MassDifferences.pdf", year.Data()), "pdf");
    
}
