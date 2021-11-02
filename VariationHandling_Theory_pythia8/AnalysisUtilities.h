#ifndef _ANALYSISUTILITIES_H_
#define _ANALYSISUTILITIES_H_
#include "BASE.h"
//#include <sys/stat.h>
#include <iostream>

#include <TFile.h>
#include <TObjString.h>
#include <TH1F.h>

/*
 * This function finds all the histograms that match a spesific pattern
 * scales with the cross secions and luminosity found in the analysis
 * constants file and returns the merged histogram.
 * 
 * @param f: the file where it will look for the histograms
 * @param year: the year for which we will use the luminosities and cross sections
 * @param pattern: the pattern that the histogram names must contain
 * @param scale: not used
 */
template <typename T>
T *MergeHistograms(TFile *f, TString year, TString pattern, float lumi = 1., bool scale = false)
{
  if (AnalysisConstants::debug)
  {
    std::cout << "Searching for pattern: " << pattern << " in file: " << f->GetName() << std::endl;
  }

  TList *keys = f->GetListOfKeys();

  T *histo = NULL;

  for (int i = 0; i < keys->GetSize(); i++)
  {
    TString name = keys->At(i)->GetName();
    if (name.Contains(pattern))
    {
      //for the first entry we get in here
      if (!histo)
      {
        histo = (T *)f->Get(name);
        //split the name into tokens baed on '_', the last one is the file identifier
        TObjArray *arr = name.Tokenize("_");
        //see if there is a cross section for this file and if yes, scale it
        std::map<TString, float>::iterator it = AnalysisConstants::crossSections[year].find(arr->Last()->GetName());
        if (it != AnalysisConstants::crossSections[year].end())
        {
          if (AnalysisConstants::debug)
          {
            std::cout << "Found: " << it->first << std::endl;
            std::cout << it->second << " " << lumi << " " << it->second * lumi << std::endl;
          }

          histo->Scale(it->second * lumi);
        }
      }
      else
      {
        T *tempHisto = (T *)f->Get(name);
        //split the name into tokens baed on '_', the last one is the file identifier
        TObjArray *arr = name.Tokenize("_");
        //see if there is a cross section for this file and if yes, scale it
        std::map<TString, float>::iterator it = AnalysisConstants::crossSections[year].find(arr->Last()->GetName());
        if (it != AnalysisConstants::crossSections[year].end())
        {
          if (AnalysisConstants::debug)
          {
            std::cout << "Found: " << it->first << std::endl;
            std::cout << it->second << " " << lumi << " " << it->second * lumi << std::endl;
          }

          tempHisto->Scale(it->second * lumi);
          histo->Add(tempHisto);
        }
      }
    }
  }

  if (scale)
    histo->Scale(1. / histo->Integral());

  return histo;
}

/*
* recursively searches the firectories below @dirname@ and returns the file that match
* the @fileToSearchFor@. The last argument specifies whether the matching should be exact
* meaning only return files that have name that is exactly the @fileToSearchFor@ parameter
* or if the file name should just contain the @fileToSearchFor@
*/

std::vector<TString> list_files(TString dirname = "C:/root/folder/", TString fileToSearchFor = ".root", bool recursive = true, bool exactMatch = false)
{
  if (AnalysisConstants::debug)
  {
    std::cout << "Scanning directory: " << dirname << std::endl;
  }
  std::vector<TString> listOfFiles;

  if (!dirname.EndsWith("/"))
  {
    dirname = TString::Format("%s/", dirname.Data());
  }

  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();

  if (files)
  {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file = (TSystemFile *)next()))
    {
      fname = file->GetName();

      if (!fname.EqualTo(".") && !fname.EqualTo(".."))
      {
        std::vector<TString> temp;
        if (recursive && file->IsDirectory())
        {
          temp = list_files(TString::Format("%s%s", dirname.Data(), fname.Data()), fileToSearchFor, exactMatch);
          listOfFiles.insert(std::end(listOfFiles), std::begin(temp), std::end(temp));
        }
        else
        {
          if (AnalysisConstants::debug)
          {
            std::cout << "-->" << fname << std::endl;
          }

          if (exactMatch)
          {
            if (fname.EqualTo(fileToSearchFor))
            {
              listOfFiles.push_back(TString::Format("%s%s", dirname.Data(), fname.Data()));
            }
          }
          else
          {
            if (fname.Contains(fileToSearchFor))
            {
              listOfFiles.push_back(TString::Format("%s%s", dirname.Data(), fname.Data()));
            }
          }
        }
      }
    }
  }
  return listOfFiles;
}

/*
* Starting from the current directory it will check of the specified directory exists and if not
* will create it. It works recursively meaning that it will create all intermediate directories
* that do not exist
*/
void CheckAndCreateDirectory(TString directoryName)
{
  if (AnalysisConstants::debug)
  {
    TObjArray *arr = directoryName.Tokenize("/");

    for (int i = 0; i < arr->GetEntries(); i++)
    {
      TString s = ((TObjString *)arr->At(i))->String();

      std::cout << s << std::endl;
    }
  }

  TString command = TString::Format("mkdir -p %s", directoryName.Data());
  gSystem->Exec(command);
}

void AnalyzeEventHistogram(TFile *file, TString year, TString variable)
{
  AnalysisConstants::clearConstants();
  AnalysisConstants::initConstants();

  if (file)
  {
    TH1F *signal = MergeHistograms<TH1F>(file, year,
                                         TString::Format("2btag_%s",
                                                         variable.Data()),
                                         AnalysisConstants::luminositiesSR[year]);
    std::cout << "2btag" << std::endl;
    std::cout << "Entries: " << signal->GetEntries();
    std::cout << " Integral: " << signal->Integral() << std::endl;
    TH1F *bkg = MergeHistograms<TH1F>(file, year,
                                      TString::Format("0btag_%s",
                                                      variable.Data()),
                                      AnalysisConstants::luminositiesCR[year]);
    std::cout << "0btag" << std::endl;
    std::cout << "Entries: " << bkg->GetEntries();
    std::cout << " Integral: " << bkg->Integral() << std::endl;
  }
  else
  {
    std::cout << "File does not exist" << std::endl;
  }
}

Float_t *GetHistogramBins(TH1F *histogram)
{
  Float_t *bins = new Float_t[histogram->GetNbinsX() + 1];

  for (int i = 1; i <= histogram->GetNbinsX() + 1; i++)
  {
    bins[i - 1] = histogram->GetBinLowEdge(i);
  }

  return bins;
}

TString GetInputFilesPath(TString year)
{
  TString path = TString::Format("%s/AnalysisFiles/%s%s%s",
                                 AnalysisConstants::baseDir.Data(),
                                 year.Data(),
                                 (AnalysisConstants::isUL ? "/UL" : ""),
                                 AnalysisConstants::currentlyWorkingDirectory[year].Data());

  return path;
}

/*
TString GetOutputDirectory(TString year, TString dir, bool )
{
  TString path = TString::Format("%s/%s",
                                 AnalysisConstants::baseDir.Data(),
                                 dir.Data());

  return path
}*/

int findIndex(TString variable)
{
  std::vector<TString>::iterator it = std::find(AnalysisConstants::variables.begin(), AnalysisConstants::variables.end(), variable);
  if (it == AnalysisConstants::variables.end())
  {
    return -1;
  }
  int index = std::distance(AnalysisConstants::variables.begin(), it);
  return index;
}

int findIndexParton(TString variable)
{
  std::vector<TString>::iterator it = std::find(AnalysisConstants::partonVariables.begin(), AnalysisConstants::partonVariables.end(), variable);
  if (it == AnalysisConstants::partonVariables.end())
  {
    return -1;
  }
  int index = std::distance(AnalysisConstants::partonVariables.begin(), it);
  return index;
}

#endif