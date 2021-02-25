#include "TString.h"
/*
This code provides the qcd bkg parameters for the template fit, taken from the template fit on the MC QCD
*/
map<TString, map<TString, Double_t>> qcdParams;


void initQCDParams()
{
/* 2016
   Floating Parameter    FinalValue +/-  Error
  --------------------  --------------------------
                qcd_b0    3.0981e-01 +/-  4.26e-01
                qcd_b1    9.0739e-01 +/-  1.43e+00
                qcd_b2    2.1589e-07 +/-  3.57e-01
                qcd_b3    2.2265e-03 +/-  6.51e-02
                qcd_b4    4.2621e-03 +/-  8.19e-03
                qcd_f1    7.6541e-01 +/-  4.49e-02
              qcd_mean    1.5378e+02 +/-  5.00e+00
             qcd_sigma    2.7498e+01 +/-  3.25e+00
*/
//----------------------------------------------------------------------------------------------------------------
	map<TString, Double_t> qcdParams16 = {{"qcd_b0", 3.0981e-01},
									  								 {"qcd_b1", 9.0739e-01},
									  							 	 {"qcd_b2", 2.1589e-07},
									  							   {"qcd_b3", 2.2265e-03},
	                                   {"qcd_b4", 4.2621e-03},
	                              	   {"qcd_f1", 7.6541e-01},
	                              	   {"qcd_mean", 1.5378e+02},
	                              	   {"qcd_sigma", 2.7498e+01}};

/* 2017
    Floating Parameter    FinalValue +/-  Error
  --------------------  --------------------------
                qcd_b0    9.7042e-01 +/-  1.48e-01
                qcd_b1    1.9814e+00 +/-  2.43e-02
                qcd_b2    6.5869e-06 +/-  3.77e-02
                qcd_b3    8.7545e-02 +/-  9.95e-03
                qcd_b4    3.0788e-02 +/-  5.94e-03
                qcd_f1    6.4990e-01 +/-  4.20e-02
              qcd_mean    1.4635e+02 +/-  2.21e+00
             qcd_sigma    3.6270e+01 +/-  1.51e+00
*/
//----------------------------------------------------------------------------------------------------------------
	map<TString, Double_t> qcdParams17 = {{"qcd_b0", 9.7042e-01},
									  {"qcd_b1", 1.9814e+00},
									  {"qcd_b2", 6.5869e-06},
									  {"qcd_b3", 8.7545e-02},
	                                  {"qcd_b4", 3.0788e-02},
	                              	  {"qcd_f1", 6.4990e-01},
	                              	  {"qcd_mean", 1.4635e+02},
	                              	  {"qcd_sigma", 3.6270e+01}};

/* 2018
    Floating Parameter    FinalValue +/-  Error
  --------------------  --------------------------
                qcd_b0    2.7840e-01 +/-  5.86e-02
                qcd_b1    5.6995e-01 +/-  1.20e-01
                qcd_b2    4.6676e-02 +/-  1.52e-02
                qcd_b3    1.5639e-03 +/-  4.30e-03
                qcd_b4    5.8243e-03 +/-  1.43e-03
                qcd_f1    7.4736e-01 +/-  1.33e-02
              qcd_mean    1.5079e+02 +/-  6.62e-01
             qcd_sigma    2.9477e+01 +/-  8.10e-01
*/
//----------------------------------------------------------------------------------------------------------------
	map<TString, Double_t> qcdParams18 = {{"qcd_b0", 2.7840e-01},
									  {"qcd_b1", 5.6995e-01},
									  {"qcd_b2", 4.6676e-02},
									  {"qcd_b3", 1.5639e-03},
	                                  {"qcd_b4", 5.8243e-03},
	                              	  {"qcd_f1", 7.4736e-01},
	                              	  {"qcd_mean", 1.5079e+02},
	                              	  {"qcd_sigma", 2.9477e+01}};

	qcdParams.insert(pair<TString, map<TString, Double_t>>("2016_preVFP",qcdParams16));
  qcdParams.insert(pair<TString, map<TString, Double_t>>("2016_postVFP",qcdParams16));
	qcdParams.insert(pair<TString, map<TString, Double_t>>("2017",qcdParams17));
	qcdParams.insert(pair<TString, map<TString, Double_t>>("2018",qcdParams18));

//----------------------------------------------------------------------------------------------------------------

}
