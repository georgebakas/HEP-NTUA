#include "qcdClosure.C"

/*
This is a function used to compute the response matrices for the two files using the Mtt cuts samples
So I will use this 
*/


void executeQCDClosure()
{
	const int NVAR =5;
	const int N_MJJ = 10;	
	//float BND_MJJ[N_MJJ+1] = {900,970,1070,1200,1350,1500,1650,1850,2100,2400,2800,4000};
	//float BND_MJJ[N_MJJ+1] = {800,880,1000,1120,1250,1400,1600,1800,2000,2300,2800,4000};
	float BND_MJJ[N_MJJ+1] = {1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000};
	
	const int N_PTJJ = 9;
	float BND_PTJJ[N_PTJJ+1] = {0,60,150,300,450,600,750,950,1100,1300};	

	const int N_YJJ = 8;
	float BND_YJJ[N_YJJ+1] = {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4};

	//jetPt will not be square
	const int N_PT = 10;
	//float BND_PT[N_PT+1] = {450,500,600,700,825,950,1100,1500};
	float BND_PT[N_PT+1] = {400,450,500,570,650,750,850,950,1100,1300,1500};

	const int N_ETA = 8;
	float BND_ETA[N_ETA+1] = {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4};
	
	
	int NBINS[NVAR] = {N_MJJ, N_PTJJ, N_YJJ, N_PT, N_ETA};
	std::vector< std::vector <Float_t> > const BND = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}, //mjj
											 {0,60,150,300,450,600,750,950,1100,1300}, //ptjj
											 {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
											 {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt
											 {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}}; //jetEta
	
  TString varReco[5]   = {"mJJ", "ptJJ", "yJJ", "jetPt", "jetEta"}; 
  TString varParton[5] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton", "partonPt", "partonEta"}; 

  TString hist;
  bool isNominal = true;
  
  const int sizeDAK8 = 5;
  
  float selMVACut = 0.2;
  float btagWP = 0.8838;
  bool isSignal = false; //false if you want QCD, true if you want Signal
  float deepAK8Wp[sizeDAK8]  = {0.2,0.3,0.4,0.5,0.6};  
  
  for(int ivar=0; ivar<1; ivar++)
  { 
	  cout<<varReco[ivar]<<endl;
	  cout<<varParton[ivar]<<endl;
	  float tempBND[NBINS[ivar]+1];
	  std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND);	
	  
	  for(int idAK8=4; idAK8<sizeDAK8; idAK8++)
	  {
		  cout<<"ENTERING idAK8 "<<idAK8<<endl;
		if(idAK8 ==4) qcdClosure(varParton[ivar], varReco[ivar], (ivar+1), NBINS[ivar], tempBND, true,  selMVACut, btagWP,isSignal,deepAK8Wp[idAK8]); 
		else qcdClosure(varParton[ivar], varReco[ivar], (ivar+1), NBINS[ivar], tempBND,false,  selMVACut, btagWP,isSignal,deepAK8Wp[idAK8]); 
	  }
  }
  
}
