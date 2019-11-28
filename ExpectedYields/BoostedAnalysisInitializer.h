#ifndef BoostedAnalysisInitializer_h
#define BoostedAnalysisInitializer_h
#include <iostream>
#include <TString.h>

class BoostedAnalysisInitializer{
  public:
    std::vector<TString> listOfFiles;
    std::vector<float> XSEC;
    std::vector<TString> histoNames;
    float LUMI;
    TString eosPath;
    BoostedAnalysisInitializer(TString year, bool isSignal);
    
  private:
    bool isSignal;
    TString year;
    
    void initialize();
    void initSignalFiles();
    void initBkgFiles();
    void initHistoNames();
    void initXsections();
};

#endif