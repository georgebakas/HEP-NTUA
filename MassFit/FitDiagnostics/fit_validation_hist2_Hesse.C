#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include <iostream>
#include <fstream>
using namespace RooFit;

void fit_validation_hist2_Hesse()
{ 

	// the signal will have known yield 
	// Signal factor = how much is the pseudodata mean yield of signal

	Float_t sig_factor = 1.;
	
	// Name of variable
	TString VAR_NAME = "mTtop";
	Float_t lim_down = 80.;
	Float_t lim_up = 600.;
	//Rebin 
	Int_t REBIN = 10;
    //qcd yield initial
	Float_t N_qcd_in;
	// directory in which the results are saved
	TString VAR_DIR = "/home/theodore/Documents/Diplomatikh/VariablesDiscrimination/2J1T/"+VAR_NAME+"/NLO_new/validation/";

	// number of processes
	const int N = 9;
	TString PROCESS[N]{
	"diboson",//0
	"data_e_cr",//1
	"Wjets_new",//2
	"DY",//3
	"tt",//4
	"s_channel",//5
	"Wt",//6
	"t_channel",//7
	"data_e"//8
	};
	

	//input files
    TFile *inf[N];
	TH1F *hist[N];
	for (int i = 0;i<N;i++){ 
		inf[i]= TFile::Open(PROCESS[i]+".root", "READ");
		//hist[i] -> GetXaxis() -> SetRange(lim_down,lim_up);
		hist[i] = (TH1F*) inf[i] -> Get(VAR_NAME);
		hist[i] -> Rebin(REBIN);
    }
	
	// take the histograms of var from files 
	TH1F *h_diboson = hist[0];
	TH1F *h_qcd = hist[1];
	TH1F *h_wjets = hist[2];
	TH1F *h_dy = hist[3];
	TH1F *h_tt = hist[4];
	TH1F *h_schannel = hist[5];
	TH1F *h_wt = hist[6];
	TH1F *h_tchannel = hist[7];

	TH1F *h_data = hist[8];

	// ---------- initial yields ------------------------
    Float_t N_diboson_in = (h_diboson->Integral());
	Float_t N_wjets_in = (h_wjets->Integral());
	Float_t N_dy_in = (h_dy-> Integral());
	Float_t N_tt_in = (h_tt->Integral());
	Float_t N_schannel_in = (h_schannel->Integral());
	Float_t N_wt_in = (h_wt->Integral());
    Float_t N_tchannel_in = (h_tchannel ->Integral());

    Float_t N_data= (h_data->Integral()); 
	
	N_qcd_in = N_data - N_diboson_in - N_wjets_in - N_dy_in - N_tt_in - N_schannel_in - N_wt_in - N_tchannel_in;
	
    h_qcd -> Scale(N_qcd_in/(h_qcd -> Integral()));

	
	//------------ Ranges of variables ------------------
	// diboson
	Float_t N_diboson_low = 0.*N_diboson_in;
	Float_t N_diboson_high =10.*N_diboson_in;
	// qcd
	Float_t N_qcd_low = 0.*N_qcd_in;
	Float_t N_qcd_high =10.*N_qcd_in;
	// wjets
	Float_t N_wjets_low = 0.*N_wjets_in;
	Float_t N_wjets_high = 10.*N_wjets_in;
	// dy
	Float_t N_dy_low = 0.*N_dy_in;
	Float_t N_dy_high = 10.*N_dy_in;
	// tt 
	Float_t N_tt_low = 0.*N_tt_in;
	Float_t N_tt_high = 10.*N_tt_in;
	// schannel 
	Float_t N_schannel_low = 0.*N_schannel_in;
	Float_t N_schannel_high = 10.*N_schannel_in;
	// wt
	Float_t N_wt_low = 0.*N_wt_in;
	Float_t N_wt_high = 10.*N_wt_in;
	// tchannel 
	Float_t N_tchannel_low = 0.*N_tchannel_in;
	Float_t N_tchannel_high = 10.*N_tchannel_in;
	
	// -------- Define number of events as roofit variables -----------
	RooRealVar N_diboson("N_diboson", "N_diboson", N_diboson_in, N_diboson_low, N_diboson_high);
	RooRealVar N_qcd("N_qcd","N_qcd",N_qcd_in, N_qcd_low, N_qcd_high);
	RooRealVar N_wjets("N_wjets","N_wjets",N_wjets_in, N_wjets_low, N_wjets_high);
	RooRealVar N_dy("N_dy","N_dy",N_dy_in, N_dy_low, N_dy_high);
	
	RooRealVar N_tt("N_tt","N_tt",N_tt_in, N_tt_low, N_tt_high);
	RooRealVar N_schannel("N_schannel","N_schannel",N_schannel_in, N_schannel_low, N_schannel_high);
	RooRealVar N_wt("N_wt","N_wt",N_wt_in, N_wt_low, N_wt_high);
	RooRealVar N_tchannel("N_tchannel","N_tchannel",N_tchannel_in, N_tchannel_low,N_tchannel_high);

	
	// ----- Define external constraints on fit variables -----
	RooGaussian N_tt_constraint("N_tt_constraint","N_tt_constraint",N_tt, RooConst(N_tt_in),RooConst(0.1*N_tt_in));

	RooGaussian N_wjets_constraint("N_wjets_constraint","N_wjets_constraint",N_wjets, RooConst(N_wjets_in),RooConst(0.2*N_wjets_in));
	
	RooGaussian N_dy_constraint("N_dy_constraint","N_dy_constraint",N_dy, RooConst(N_dy_in),RooConst(0.2*N_dy_in));
	
	RooGaussian N_wt_constraint("N_wt_constraint","N_wt_constraint",N_wt, RooConst(N_wt_in),RooConst(0.2*N_wt_in));
	
	RooGaussian N_schannel_constraint("N_schannel_constraint","N_schannel_constraint",N_schannel, RooConst(N_schannel_in),RooConst(0.2*N_schannel_in));
	
	RooGaussian N_diboson_constraint("N_diboson_constraint","N_diboson_constraint",N_diboson, RooConst(N_diboson_in),RooConst(0.2*N_diboson_in));

	//------- define the variable and the range ------
	RooRealVar variable(VAR_NAME,VAR_NAME,lim_down,lim_up);
	

    //---------- make histograms ------------------ 
	RooDataHist diboson("diboson","diboson", variable, h_diboson);
	RooDataHist qcd("qcd", "qcd", variable, h_qcd);
	RooDataHist wjets("wjets", "wjets", variable, h_wjets);
	RooDataHist dy("dy", "dy", variable, h_dy);
	RooDataHist tt("tt", "tt", variable, h_tt);
	RooDataHist schannel("schannel", "schannel", variable, h_schannel);
	RooDataHist wt("wt", "wt",variable, h_wt);
	RooDataHist tchannel("tchannel", "tchannel", variable, h_tchannel);
	//RooDataHist data("data","data",variable, h_data);

    //make pdfs
	RooHistPdf diboson_pdf("diboson_pdf", "diboson_pdf", variable, diboson);
	RooHistPdf qcd_pdf("qcd_pdf", "qcd_pdf", variable, qcd);
	RooHistPdf wjets_pdf("wjets_pdf", "wjets_pdf", variable, wjets);
	RooHistPdf dy_pdf("dy_pdf", "dy_pdf", variable, dy);
	RooHistPdf tt_pdf("tt_pdf", "tt_pdf", variable, tt);
	RooHistPdf schannel_pdf("schannel_pdf", "schannel_pdf", variable, schannel);
	RooHistPdf wt_pdf("wt_pdf", "wt_pdf", variable, wt);
	RooHistPdf tchannel_pdf("tchannel_pdf", "tchannel_pdf", variable, tchannel);


    //Make the fit model := Sum_[i]{ N[i]*Process_pdf[i] }

	RooAddPdf model("model", "model", RooArgSet(tchannel_pdf, wt_pdf, schannel_pdf,tt_pdf,wjets_pdf,dy_pdf,diboson_pdf,qcd_pdf), RooArgSet(N_tchannel, N_wt, N_schannel,N_tt, N_wjets,N_dy,N_diboson, N_qcd ));

	variable.setBins(tchannel.numEntries());

	// generated yields for bkgs    
	Double_t N_gen_diboson,N_gen_qcd,N_gen_wjets,N_gen_tt,N_gen_dy,
	N_gen_schannel,N_gen_wt;
	Double_t N_gen_data;

	TRandom3 ran;
	// ----- Histograms to save the results ------------
	TH1F *h_tchannel_fit = new TH1F("tchannel yield","tchannel yield",500,-0.5*sig_factor*N_tchannel_in,0.5*sig_factor*N_tchannel_in); // histogram of fit results
	TH1F *h_tchannel_fit_pull = new TH1F("tchannel yield pull","tchannel yield pull",50,-10.,10.); 
	TH1F *h_tchannel_fit_error = new TH1F("tchannel yield error","error",1000,0.,0.1*sig_factor*N_tchannel_in);
	TH1F *h_cor_tchan_wjets = new TH1F("correlation t-channel/W+jets","correlation t-channel/W+jets",1000,-1.,1.);
	TH1F *h_cor_tchan_tt = new TH1F("correlation t-channel/tt","correlation t-channel/tt",1000,-1.,1.);

	// wjets
	TH1F *h_wjets_fit = new TH1F("wjets yield","wjets yield",500,-0.5*N_wjets_in,0.5*N_wjets_in); // histogram of fit results
	TH1F *h_wjets_fit_pull = new TH1F("wjets yield pull","wjets yield pull",50,-10.,10.); 
	TH1F *h_wjets_fit_error = new TH1F("wjets yield error","error",1000,0.,0.1*N_tt_in);

	//tt
	TH1F *h_tt_fit = new TH1F("tt yield","tt yield",500,-0.5*N_tt_in,0.5*N_tt_in); // histogram of fit results
	TH1F *h_tt_fit_pull = new TH1F("tt yield pull","tt yield pull",50,-10.,10.); 
	TH1F *h_tt_fit_error = new TH1F("tt yield error","error",1000,0.,0.1*N_tt_in);

	//qcd
	TH1F *h_qcd_fit = new TH1F("qcd yield","qcd yield",500,-0.5*N_qcd_in,0.5*N_qcd_in); // histogram of fit results
	TH1F *h_qcd_fit_pull = new TH1F("qcd yield pull","qcd yield pull",50,-10.,10.); 
	TH1F *h_qcd_fit_error = new TH1F("qcd yield error","error",1000,0.,0.1*N_qcd_in);

	//check if the fitted yield is on limits
	// on any variable
	bool lim_check;
	TH1F *h_lim_check = new TH1F("lim_check","lim_check",2,0,2);
	// --------------------------------------------------
	
	for (int i=0;i<400;i++){
		cout << "runNo:  " << i << endl; 
		ran.SetSeed(0); // change seed
		// generate random yields for bkgs
		
		/*
		N_gen_diboson = ran.Gaus(N_diboson_in,(N_diboson_high-N_diboson_in)/4.);
		N_gen_qcd = ran.Gaus(N_qcd_in,N_qcd_in*0.2);
		N_gen_wjets = ran.Gaus(N_wjets_in,0.05*N_wjets_in);
		N_gen_dy = ran.Gaus(N_dy_in,0.05*N_dy_in);
		N_gen_tt = ran.Gaus(N_tt_in,0.025*N_tt_in);
		N_gen_schannel = ran.Gaus(N_schannel_in,0.05*N_schannel_in);
		N_gen_wt = ran.Gaus(N_wt_in,0.05*N_wt_in);
		*/
		N_gen_diboson = ran.Poisson(N_diboson_in);
		N_gen_qcd = ran.Poisson(N_qcd_in);
		N_gen_wjets = ran.Poisson(N_wjets_in);
		N_gen_dy = ran.Poisson(N_dy_in);
		N_gen_tt = ran.Poisson(N_tt_in);
		N_gen_schannel = ran.Poisson(N_schannel_in);
		N_gen_wt = ran.Poisson(N_wt_in);
		/*
		N_gen_diboson = ran.Gaus(N_diboson_in,sqrt(N_diboson_in));
		N_gen_qcd = ran.Gaus(N_qcd_in,sqrt(N_qcd_in));
		N_gen_wjets = ran.Gaus(N_wjets_in,sqrt(N_wjets_in));
		N_gen_dy = ran.Gaus(N_dy_in,sqrt(N_dy_in));
		N_gen_tt = ran.Gaus(N_tt_in,sqrt(N_tt_in));
		N_gen_schannel = ran.Gaus(N_schannel_in,sqrt(N_schannel_in));
		N_gen_wt = ran.Gaus(N_wt_in,sqrt(N_wt_in));
		*/
		N_gen_data = N_gen_diboson + N_gen_qcd + N_gen_wjets + N_gen_dy + N_gen_tt + N_gen_schannel + N_gen_wt + sig_factor*N_tchannel_in;

		//-------  generate pseudodata -----------------------
    	RooAddPdf pseudodata("pseudodata","pseudodata", RooArgSet(tchannel_pdf, wt_pdf, schannel_pdf,tt_pdf,wjets_pdf,dy_pdf,diboson_pdf,qcd_pdf), RooArgSet( RooConst(sig_factor*N_tchannel_in), RooConst(N_gen_wt), RooConst(N_gen_schannel) , RooConst(N_gen_tt), RooConst(N_gen_wjets) ,RooConst(N_gen_dy),RooConst(N_gen_diboson), RooConst(N_gen_qcd) ) );
	
		//sample N_gen_data events
		RooDataSet *data_gen = pseudodata.generate(variable,N_gen_data)	;
		// create a binned dataset
		RooDataHist *data = data_gen -> binnedClone();
	
   	    //Fit the model to data
		RooFitResult *fit_result = model.fitTo(*data,ExternalConstraints( RooArgSet(N_tt_constraint,  N_wjets_constraint,  N_dy_constraint,  N_wt_constraint ,  N_schannel_constraint,  N_diboson_constraint) ),Save(),Minos(kFALSE),Hesse(kTRUE),Verbose(kFALSE) , Timer(kFALSE),  PrintLevel(-1), PrintEvalErrors(-1) );
		//RooFitResult *fit_result = model.fitTo(*data,Save());

		
		lim_check =( 
		(N_tt.getVal() <= 1.05*(1.-3*0.1)*N_tt_in) ||
		(N_tt.getVal() >= 0.95*(1.+3*0.1)*N_tt_in) ||
		(N_wjets.getVal() <=  1.05*(1.-3*0.2)*N_wjets_in) ||
		(N_wjets.getVal() >= 0.95*(1.+3*0.2)*N_wjets_in) 
		 );

		h_tchannel_fit -> Fill((N_tchannel.getVal()-sig_factor*N_tchannel_in));
		h_tchannel_fit_pull -> Fill((N_tchannel.getVal()-sig_factor*N_tchannel_in)/N_tchannel.getError());
		h_tchannel_fit_error -> Fill(N_tchannel.getError());
		h_cor_tchan_wjets -> Fill(fit_result-> correlation("N_tchannel","N_wjets"));
		h_cor_tchan_tt -> Fill(fit_result-> correlation("N_tchannel","N_tt"));	

		h_wjets_fit -> Fill((N_wjets.getVal()-N_wjets_in));
		h_wjets_fit_pull -> Fill((N_wjets.getVal()-N_wjets_in)/N_wjets.getError());
		h_wjets_fit_error -> Fill(N_wjets.getError());

		h_tt_fit -> Fill((N_tt.getVal()- N_tt_in));
		h_tt_fit_pull -> Fill((N_tt.getVal()- N_tt_in)/N_tt.getError());
		h_tt_fit_error -> Fill(N_tt.getError());

		h_qcd_fit -> Fill((N_qcd.getVal()- N_qcd_in));
		h_qcd_fit_pull -> Fill((N_qcd.getVal()- N_qcd_in)/N_qcd.getError());
		h_qcd_fit_error -> Fill(N_qcd.getError());
		
		h_lim_check -> Fill(lim_check);

		// Reset the RooRealVal's of the yields to the initial given yields
		// so that we don't have influence of one fit to another	
		N_diboson.setVal(N_diboson_in);
		N_qcd.setVal(N_qcd_in);
		N_wjets.setVal(N_wjets_in);
		N_dy.setVal(N_dy_in);
		N_tt.setVal(N_tt_in);
		N_schannel.setVal(N_schannel_in);
		N_wt.setVal(N_wt_in);
		N_tchannel.setVal(N_tchannel_in);
	}	
	
	// ----------- Draw histograms --------------
	 
	TF1 *f = new TF1("gaus_fit","ROOT::Math::gaussian_pdf(x,1.,0.)",-4.,4.); 

	TCanvas *can1 = new TCanvas("can1","can1");
	h_tchannel_fit->SetTitle("");
	h_tchannel_fit -> GetXaxis() -> SetTitle("N^{fit}_{t-chan}-N^{in}_{t-chan}");
	h_tchannel_fit -> Draw("hist e");

	can1 -> SaveAs(VAR_DIR+"fit_hist.root");
	
	h_tchannel_fit_pull -> Scale(1./(h_tchannel_fit_pull -> Integral("width")));
	
	TCanvas *can2 = new TCanvas("can2","can2");
	h_tchannel_fit_pull->SetTitle("");
	h_tchannel_fit_pull -> GetXaxis() -> SetTitle("pull N^{fit}_{t-chan}");
	h_tchannel_fit_pull -> Draw("hist e");
	f -> Draw("same");
	can2 -> SaveAs(VAR_DIR+"fit_pull.root");

	TCanvas *can3 = new TCanvas("can3","can3");
	h_tchannel_fit_error->SetTitle("");
	h_tchannel_fit_error -> GetXaxis() -> SetTitle("#deltaN^{fit}_{t-chan}");
	h_tchannel_fit_error -> Draw("hist e");
	
	can3 -> SaveAs(VAR_DIR+"fit_error.root");

	TCanvas *can4 = new TCanvas("can4","can4");
	h_cor_tchan_wjets->SetTitle("");
	h_cor_tchan_wjets -> GetXaxis() -> SetTitle("correlation t-chan./W+jets");
	h_cor_tchan_wjets -> Draw("hist e");
	can4 -> SaveAs(VAR_DIR+"cor_tchan_wjets.root");

	TCanvas *can5 = new TCanvas("can5","can5");
	h_cor_tchan_tt->SetTitle("");
	h_cor_tchan_tt -> GetXaxis() -> SetTitle("correlation t-chan./tt");
	h_cor_tchan_tt -> Draw("hist e");
	can5 -> SaveAs(VAR_DIR+"cor_tchan_tt.root");
	
	TCanvas *can6 = new TCanvas("can6","can6");
	h_lim_check ->SetTitle("");
	h_lim_check -> GetXaxis() -> SetTitle("lim_check");
	h_lim_check -> Draw("hist e");
	can6 -> SaveAs(VAR_DIR+"lim_check_test.root");

	TCanvas *can7 = new TCanvas("can7","can7");
	h_wjets_fit->SetTitle("");
	h_wjets_fit -> GetXaxis() -> SetTitle("N^{fit}_{wjets}-N^{in}_{wjets}");
	h_wjets_fit -> Draw("hist e");
	
	can7 -> SaveAs(VAR_DIR+"fit_hist_wjets.root");
	
	h_wjets_fit_pull -> Scale(1./(h_wjets_fit_pull -> Integral("width")));
	TCanvas *can8 = new TCanvas("can8","can8");
	h_wjets_fit_pull->SetTitle("");
	h_wjets_fit_pull -> GetXaxis() -> SetTitle("pull N^{fit}_{wjets}");
	h_wjets_fit_pull -> Draw("hist e");
	f -> Draw("same");
	can8 -> SaveAs(VAR_DIR+"fit_pull_wjets.root");

	TCanvas *can9 = new TCanvas("can9","can9");
	h_tt_fit->SetTitle("");
	h_tt_fit -> GetXaxis() -> SetTitle("N^{fit}_{tt}-N^{in}_{tt}");
	h_tt_fit -> Draw("hist e");
	
	can9 -> SaveAs(VAR_DIR+"fit_hist_tt.root");
	
	h_tt_fit_pull -> Scale(1./(h_tt_fit_pull -> Integral("width")));

	TCanvas *can10 = new TCanvas("can10","can10");
	h_tt_fit_pull->SetTitle("");
	h_tt_fit_pull -> GetXaxis() -> SetTitle("pull N^{fit}_{tt}");
	h_tt_fit_pull -> Draw("hist e");
	f -> Draw("same");
	can10 -> SaveAs(VAR_DIR+"fit_pull_tt.root");

	TCanvas *can11 = new TCanvas("can11","can11");
	h_qcd_fit->SetTitle("");
	h_qcd_fit -> GetXaxis() -> SetTitle("N^{fit}_{qcd}-N^{in}_{qcd}");
	h_qcd_fit -> Draw("hist e");
	
	can11 -> SaveAs(VAR_DIR+"fit_hist_qcd.root");
	
	h_qcd_fit_pull -> Scale(1./(h_qcd_fit_pull -> Integral("width")));

	TCanvas *can12 = new TCanvas("can12","can12");
	h_qcd_fit_pull->SetTitle("");
	h_qcd_fit_pull -> GetXaxis() -> SetTitle("pull N^{fit}_{qcd}");
	h_qcd_fit_pull -> Draw("hist e");
	f -> Draw("same");
	can12 -> SaveAs(VAR_DIR+"fit_pull_qcd.root");
	
	TCanvas *can13 = new TCanvas("can13","can13");
	h_wjets_fit_error->SetTitle("");
	h_wjets_fit_error -> GetXaxis() -> SetTitle("#deltaN^{fit}_{wjets}");
	h_wjets_fit_error -> Draw("hist e");

	can13 -> SaveAs(VAR_DIR+"fit_error_wjets.root");
		
	TCanvas *can14 = new TCanvas("can14","can14");
	h_tt_fit_error->SetTitle("");
	h_tt_fit_error -> GetXaxis() -> SetTitle("#deltaN^{fit}_{tt}");
	h_tt_fit_error -> Draw("hist e");

	can14 -> SaveAs(VAR_DIR+"fit_error_tt.root");

	TCanvas *can15 = new TCanvas("can15","can15");
	h_qcd_fit_error->SetTitle("");
	h_qcd_fit_error -> GetXaxis() -> SetTitle("#deltaN^{fit}_{qcd}");
	h_qcd_fit_error -> Draw("hist e");

	can15 -> SaveAs(VAR_DIR+"fit_error_qcd.root");
}
