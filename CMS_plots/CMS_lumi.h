#ifndef __CMS_LUMI_H__
#define __CMS_LUMI_H__

#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TString.h"
//#include "TASImage.h"

#include <iostream>
#include "TString.h"

//
// Global variables
//

TString cmsText = "CMS";
float cmsTextFont = 61; // default is helvetic-bold

bool writeExtraText = false;
TString extraText = "Work In Progress";
float extraTextFont = 52; // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize = 0.6;
float lumiTextOffset = 0.2;
float cmsTextSize = 0.75;   // 1.25;
float cmsTextOffset = 0.01; // only used in outOfFrame version
float extraTextFactor = 0.15;

float relPosX = 0.045;
float relPosY = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize = 0.76;

// TString lumi_13TeV = "20.1 fb^{-1}";
TString lumi_13TeV = "41.5 fb^{-1}";
TString lumi_8TeV = "19.7 fb^{-1}";
TString lumi_7TeV = "5.1 fb^{-1}";
TString lumi_sqrtS = "";
TString lumi_2016_preVFP = "19.5 fb^{-1}";
TString lumi_2016_postVFP = "16.5 fb^{-1}";
TString lumi_2017 = "41.5 fb^{-1}";
TString lumi_2018 = "59.7 fb^{-1}";
TString lumi_RunII = "137.1 fb^{-1}";

bool drawLogo = false;

void CMS_lumi(TPad *pad, int iPeriod = 3, int iPosX = 10);

#endif