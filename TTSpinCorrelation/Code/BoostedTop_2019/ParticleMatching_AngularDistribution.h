#ifndef ParticleMatching_AngularDistribution_h
#define ParticleMatching_AngularDistribution_h

class ParticleMatching_AngularDistribution
{
  public:
    explicit ParticleMatching_AngularDistribution();
	virtual ~ParticleMatching_AngularDistribution();
	
  private:  
    std::vector<float> *pt_Jet,*eta_Jet,*phi_Jet,*mass_Jet, *partonMatchDR, *pt_SubJet0, *pt_SubJet1, *eta_SubJet0, *eta_SubJet1,
						*phi_SubJet0, *phi_SubJet1, *mass_Subjet0, *mass_Subjet1, *matched_jetBtagSub0, *matched_jetBtagSub1;
    std::vector<int> *partonMatchIdx;
    
};

#endif
