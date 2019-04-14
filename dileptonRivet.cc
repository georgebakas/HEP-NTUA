	// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/RivetMT2.hh"
#include <bitset>
#include <array>
#include <unordered_set>
//#include "mt2_bisect.h"

namespace  Rivet {
  

    // @brief 
    /// @author A. Grohsjean <alexander.grohsjean@desy.de>
    
    class TTbarSpinDensityMatrix : public Analysis {
    public:
	//constructor 
	TTbarSpinDensityMatrix() : Analysis("TTbarSpinDensityMatrix") {}  
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	void init() {
	    
	    _dilepton = 0;
	    _selected = 0 ;

	    //cut and cutting related values 
	    _njets = 2;
	    _nbtags = 2; 
	    _mllmin = 0.;
	    _metmin = 0.; 
	    _mllwindow = -1.;
	    
	    //add projections 
	    
	    // full phase space in eta and pt
	    Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT > 1.0*MeV;

            // all final state particles
	    FinalState fs(eta_full);

            // get photons to dress leptons
	    IdentifiedFinalState photons(fs);
            photons.acceptIdPair(PID::PHOTON);

            // dressed and vetoed electrons
	    IdentifiedFinalState el_id(fs);
            el_id.acceptIdPair(PID::ELECTRON);
            PromptFinalState electrons(el_id);
            electrons.acceptTauDecays(true);
            addProjection(electrons, "electrons");
            DressedLeptons dressedelectrons(photons, electrons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 20*GeV, true, true);
            addProjection(dressedelectrons, "dressedelectrons");
            DressedLeptons vetodressedelectrons(photons, electrons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 15*GeV, true, true);
            addProjection(vetodressedelectrons, "vetodressedelectrons");
            DressedLeptons ewdressedelectrons(photons, electrons, 0.1, eta_full, true, true);
            addProjection(ewdressedelectrons, "ewdressedelectrons");

            // dressed and vetoed muons
	    IdentifiedFinalState mu_id(fs);
            mu_id.acceptIdPair(PID::MUON);
            PromptFinalState muons(mu_id);
            muons.acceptTauDecays(true);
            addProjection(muons, "muons");
            vector<pair<double, double> > eta_muon;
            DressedLeptons dressedmuons(photons, muons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 20*GeV, true, true);
            addProjection(dressedmuons, "dressedmuons");
            DressedLeptons vetodressedmuons(photons, muons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 15*GeV, true, true);
            addProjection(vetodressedmuons, "vetodressedmuons");
            DressedLeptons ewdressedmuons(photons, muons, 0.1, eta_full, true, true);
            addProjection(ewdressedmuons, "ewdressedmuons");

            // neutrinos, DM particles and MET
	    // neutrinos 
	    IdentifiedFinalState nu_id;
            nu_id.acceptNeutrinos();
            PromptFinalState neutrinos(nu_id);
            neutrinos.acceptTauDecays(true);
            addProjection(neutrinos, "neutrinos");
	    // DM particles 
	    IdentifiedFinalState chi_id(fs);
	    //agrohsje : NLO and default cms model orthogonal in IDs -> no conflict when using both 
	    vector<PdgId> chi_pid;
            chi_pid.clear();
	    chi_pid.push_back(9100022);
	    chi_pid.push_back(52); // 2hdm:52; dmsimp:52; (pseudo)scalar 9100022 
	    chi_id.acceptIdPairs(chi_pid); 
            PromptFinalState chis(chi_id);
            addProjection(chis, "chis");

	    
	    // jet clustering
	    VetoedFinalState vfs;
            vfs.addVetoOnThisFinalState(ewdressedelectrons);
            vfs.addVetoOnThisFinalState(ewdressedmuons);
            vfs.addVetoOnThisFinalState(neutrinos);
            FastJets jets(vfs, FastJets::ANTIKT, 0.4);
            jets.useInvisibles();
            addProjection(jets, "jets");
	    
	    
	    //book histograms
	    _bookhistograms();
	    
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	void analyze(const Event& event) {
	    
	  //agrohsje std::cout<<"---------------NEW EVENT-------------------"<<std::endl;
	  
	  //clear all vectors 
	  _clearvectors();
	  
	  //get event weight 
	  const double weight = event.weight();
	  //std::cout<<" weight is " << weight << std::endl;
	  // get selected objects using projections
	  _dressedelectrons     = sortByPt(applyProjection<DressedLeptons>(event, "dressedelectrons").dressedLeptons());
	  _vetodressedelectrons = applyProjection<DressedLeptons>(event, "vetodressedelectrons").dressedLeptons();
	  _dressedmuons         = sortByPt(applyProjection<DressedLeptons>(event, "dressedmuons").dressedLeptons());
	  _vetodressedmuons     = applyProjection<DressedLeptons>(event, "vetodressedmuons").dressedLeptons();
	  _neutrinos            = applyProjection<PromptFinalState>(event, "neutrinos").particlesByPt();
	  _chis                 = applyProjection<PromptFinalState>(event, "chis").particlesByPt();
	  _alljets              = applyProjection<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 2.5);
	  
	  // split all jets in light and b-tagged jets 
	  foreach (Jet jet, _alljets)
	    if (!jet.bTags().empty()) _bjets.push_back(jet);
	    else _ljets.push_back(jet);
	  
	  // calculate missing et 
	  _pmet.setE(0.); _pmet.setPx(0.);  _pmet.setPy(0.);  _pmet.setPz(0.); 
	  foreach(const Particle& p, _neutrinos) _pmet+=p.momentum(); 
	  foreach(const Particle& p, _chis)      _pmet+=p.momentum();
	  
	  // fill DM particles 
	  foreach(const Particle& p, _chis) 
	    _fmchi.push_back(p.momentum());
	    
	  // search for last copies of tops and DM mediators 
	  foreach (const GenParticle* gp, particles(event.genEvent())) {
	    const PdgId pid = gp->pdg_id();
	    if( (abs(pid) == 9100000 || abs(pid) == 54 || abs(pid) == PID::TQUARK) && _islastcopy(*gp) ) 
	      _selectedgenparticles.push_back(const_cast<GenParticle*> (gp));
	  }
	  
	  // search for the final state decay particles of the selected intermediate particles 	    
	  foreach (const GenParticle* gp, _selectedgenparticles){
	    const PdgId pid = gp->pdg_id();
	    
	    //agrohsje NLO and default cms model orthogonal in IDs -> no conflict when using both 
	    if(abs(pid) == 9100000 || abs(pid) == 54 || abs(pid) == 55 || abs(pid) == 36 ){
	      _fmphi.push_back(gp->momentum());
	    }// if DM particle 

	    // top part starts here 
	    if(abs(pid) == PID::TQUARK){
	      _fmtop.push_back(gp->momentum());
	      //search for stable decay particles from top 
	      HepMC::GenVertex* gv = const_cast<HepMC::GenVertex*>(gp->end_vertex());
	      for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::descendants); it!=gv->particles_end(HepMC::descendants); ++it) {
		//select stable, final state electrons/muons 
		if ((*it)->status() == 1 && (*it)->end_vertex() == NULL && 
		    //agrohsje added taus for debugging 
		    (abs((*it)->pdg_id()) == PID::ELECTRON || abs((*it)->pdg_id()) == PID::MUON )) {
		  //std::cout<<" found lepton "<<(*it)->pdg_id()<<" w./ status "<< (*it)->status() << std::endl;
		  //check if prompt 
		  bool isprompt(true); 
		  GenVertex* prodVtx = (*it)->production_vertex();
		  // orphaned particle, has to be assume false
		  if (prodVtx == NULL) isprompt=false;
		  const pair<GenParticle*, GenParticle*> beams = prodVtx->parent_event()->beam_particles();
		  foreach (const GenParticle* ancestor, Rivet::particles(prodVtx, HepMC::ancestors)) {
		    const PdgId pid_anc = ancestor->pdg_id();
		    // no non-standard statuses or beams to be used in decision making
		    if (ancestor->status() != 2) continue; 
		    // PYTHIA6 uses status 2 for beams, I think... (sigh)
		    if (ancestor == beams.first || ancestor == beams.second) continue; 
		    // PYTHIA6 also uses status 2 for some partons, I think... (sigh)
		    if (PID::isParton(pid_anc)) continue; 
		    // prompt particles can't be from hadron decays
		    if (PID::isHadron(pid_anc)) isprompt=false;
		    //ban particles from muon decays (permitting muon copies)
		    if (abs(pid_anc) == PID::MUON && abs((*it)->pdg_id()) != PID::MUON ) isprompt=false; 
		    //disable electrons/muons from tau decay 
		    //agrohsje for debugging 
		    if (abs(pid_anc) == PID::TAU) isprompt=false;
		  }
		  //store top and iHepMC::descendantsts final decay lepton 
		  if(isprompt && (*it)->pdg_id()*pid < 0) _toplep.push_back(make_pair(Particle(*gp),Particle(*it)));
		}// if final state lepton 
	      }// loop over all final state particles
	    }// if top  
	  }// loop all gen particles 
	  if (_toplep.size()==2) _dilepton++;
	  if (_toplep.size()>2) std::cout<<"WARNING:::Leptons size is "<<_toplep.size()<<std::endl; 
	  // fill in _toplepmomenta to calculate spin observables 
	  if(_toplep.size()>=2){
	    if(_toplep[0].first.pid()>0 && _toplep[1].first.pid()<0){
	      _toplepmomenta.push_back(make_pair(_toplep[0].first.momentum(),_toplep[0].second.momentum()));  
	      _toplepmomenta.push_back(make_pair(_toplep[1].first.momentum(),_toplep[1].second.momentum()));  
	    }
	    else if(_toplep[0].first.pid()<0 && _toplep[1].first.pid()>0){
	      //iter_swap(_toplep.begin(), _toplep.begin() + 1);
	      _toplepmomenta.push_back(make_pair(_toplep[1].first.momentum(),_toplep[1].second.momentum()));  
	      _toplepmomenta.push_back(make_pair(_toplep[0].first.momentum(),_toplep[0].second.momentum()));  
	    }
	  }

	  // calculate event observables 
	  if ( _toplep.size()>1 ){
	    // mt2_bisect::mt2 mt2;
	    // double l1[3] = { 0.0, _toplep[0].second.momentum().px(), _toplep[0].second.momentum().py() };
	    // double l2[3] = { 0.0, _toplep[1].second.momentum().px(), _toplep[1].second.momentum().py() };
	    // double ptmiss[3]={ 0.0, _pmet.px(), _pmet.py()};  // [mass, px, py]
	    // mt2.set_momenta(l1,l2,ptmiss); // give the lepton momenta and missing energy                                                                                                                
	    // mt2.set_mn( 0.0 ); // invisible particle is neutrino, so massless
	    _mt2ll = mT2::mT2(_toplep[0].second.momentum(), _toplep[1].second.momentum(), _pmet, 0.0); // zero mass of invisibles 
	    //std::cout<<" desy routine " << mt2.get_mt2() <<" rivet " << _mt2ll << std::endl;
	    _cem = _mt2ll + 0.2 * ( 200. - _pmet.pT() );
	  }
	  else{
	    _mt2ll = -1.;
	    _cem = -1.; 
	  }
	  // calculate spin density observables
	  
	  _reconstructspinobservables();
	  
	  //fill all histos according to event selection 
	  _fillhistos(0, weight);
	  unsigned emu_bits = 0, ee_bits = 0, mumu_bits = 0;
	  if (_passemu(emu_bits)||_passee(ee_bits)||_passmumu(mumu_bits)) {
	    _fillhistos(1, weight);
	    _selected++;
		if(!_overlap()) _fillhistos(2, weight);
	  }
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	void finalize() {
	    
	    // agrohsje 
	  
	  std::cout<<" Number of reconstructable dilepton events " << _dilepton<<std::endl;
	    std::cout<<" Number of selected dilepton events " << _selected <<std::endl;
	    std::cout<<" cos cos before selection:: "
		     << " mean " << _h_coscos[0]->xMean() 
	      //<< " variance " << _h_coscos[0]->xVariance() 
	      //<< " std dev " << _h_coscos[0]->xStdDev() 
		     <<std::endl;
	    std::cout<<" cos phi before selection:: "
		     << " mean " << _h_cosphi[0]->xMean() 
	      //<< " variance " << _h_cosphi[0]->xVariance() 
	      //<< " std dev " << _h_cosphi[0]->xStdDev() 
		     <<std::endl;
	  	
	}
      
    private:
      
	/// @name Physics object helper functions
	//@{
	
	void _clearvectors(){
	    
	    _selectedgenparticles.clear();
	    _bjets.clear();
            _ljets.clear();
	    _toplep.clear();
	    _toplepmomenta.clear();
	    _fmphi.clear();
	    _fmchi.clear();
	    _fmtop.clear();
	    
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	
	void _bookhistograms(){
	    
	  
	  // book remaining histograms 
	  for(unsigned i=0; i<_h_mtop.size(); i++){
	    // agrohsje 
	    _h_met[i]      = bookHisto1D(("met_cut_"+std::to_string(i)).c_str(),      20,   0.,  600.);
	    // _h_mt2ll[i]    = bookHisto1D(("mt2ll_cut_"+std::to_string(i)).c_str(),    20,   0.,  200.);
	    // _h_cem[i]      = bookHisto1D(("cem_cut_"+std::to_string(i)).c_str(),      20,   0.,  400.);
	    
	    // // top properties 
	    std::vector< double > binedges;
	    //binning kidonakis
	    // binedges.push_back(0.);
	    // binedges.push_back(65.);
	    // binedges.push_back(125.); 
	    // binedges.push_back(200.);
	    // binedges.push_back(290.);
	    // binedges.push_back(400.);
	    // binedges.push_back(550.);
	    
	    //binning mitov 
	    // binedges.push_back(0.);
	    // binedges.push_back(50.);
	    // binedges.push_back(100.);
	    // binedges.push_back(150.);
	    // binedges.push_back(200.);
	    // binedges.push_back(250.);
	    // binedges.push_back(300.);
	    // binedges.push_back(350.);
	    // binedges.push_back(400.);
	    // binedges.push_back(450.);
	    // binedges.push_back(500.);
	    // binedges.push_back(550.);
	    // binedges.push_back(600.);
	    // binedges.push_back(650.);
	    // binedges.push_back(700.);
	    // binedges.push_back(800.);
	    // binedges.push_back(900.);
	    // binedges.push_back(1000.);
	    // binedges.push_back(1100.);
	    // binedges.push_back(1200.);
	    // binedges.push_back(1400.);
	    // binedges.push_back(1600.);
	    // binedges.push_back(1800.);
	    // binedges.push_back(2000.);
	    // binedges.push_back(2200.);
	    // binedges.push_back(2600.);
	    // binedges.push_back(3000.);
	    
	    // _h_pttop[i]    = bookHisto1D(("pttop_cut_"+std::to_string(i)).c_str(), binedges); 
	    
	    binedges.clear();
	    binedges.push_back(340.);
	    binedges.push_back(380.);
	    binedges.push_back(420.);
	    binedges.push_back(460.);
	    binedges.push_back(500.);
	    binedges.push_back(540.);
	    binedges.push_back(580.);
	    binedges.push_back(620.);
	    binedges.push_back(660.);
	    binedges.push_back(700.);
	    binedges.push_back(780.);
	    binedges.push_back(860.);
	    binedges.push_back(940.);
	    binedges.push_back(1020.);
	    binedges.push_back(1100.);
	    binedges.push_back(1220.);
	    binedges.push_back(1340.);
	    binedges.push_back(1460.);
	    binedges.push_back(1580.);
	    binedges.push_back(1700.);
	    binedges.push_back(1900.);
	    binedges.push_back(2100.);
	    binedges.push_back(2300.);
	    binedges.push_back(2500.);
	    binedges.push_back(3000.);
	    binedges.push_back(3500.);
	    binedges.push_back(4000.);
	    binedges.push_back(5000.);
	    binedges.push_back(6000.);
	    
	    _h_mttbar[i]   = bookHisto1D(("mttbar_cut_"+std::to_string(i)).c_str(),  binedges); 
	    _h_ptttbar[i]  = bookHisto1D(("ptttbar_cut_"+std::to_string(i)).c_str(),  52,  0.1,  3.25); 
	    _h_mtop[i]     = bookHisto1D(("mtop_cut_"+std::to_string(i)).c_str(),     50, 160.,  185.); 
	    _h_pttop[i]    = bookHisto1D(("pttop_cut_"+std::to_string(i)).c_str(),    70,   0.,  350.); 
	    _h_etatop[i]   = bookHisto1D(("etatop_cut_"+std::to_string(i)).c_str(),   50,  -6.,    6.); 
	    _h_phitop[i]   = bookHisto1D(("phitop_cut_"+std::to_string(i)).c_str(),   10,    0,  2*PI);  
	   
	    // agrohsje 
	    // // DM properties 
	    // _h_mchi[i]     = bookHisto1D(("mchi_cut_"+std::to_string(i)).c_str(),     20,   0.,   20.); 
	    // _h_ptchi[i]    = bookHisto1D(("ptchi_cut_"+std::to_string(i)).c_str(),    30,   0.,  150.); 
	    // _h_etachi[i]   = bookHisto1D(("etachi_cut_"+std::to_string(i)).c_str(),   40,  -5.,    5.); 
	    // _h_phichi[i]   = bookHisto1D(("phichi_cut_"+std::to_string(i)).c_str(),   10,    0,  2*PI);  
	    
	    // // mediator properties from dark matter or directly 
	     _h_mphi[i]     = bookHisto1D(("mphi_cut_"+std::to_string(i)).c_str(),    160,   0.,  800.); 
	     _h_ptphi[i]    = bookHisto1D(("ptphi_cut_"+std::to_string(i)).c_str(),    50,   0.,  250.); 
	     _h_etaphi[i]   = bookHisto1D(("etaphi_cut_"+std::to_string(i)).c_str(),   50,  -5.,    5.); 
	    // _h_phiphi[i]   = bookHisto1D(("phiphi_cut_"+std::to_string(i)).c_str(),   10,    0,  2*PI); 
	    
	    // spin density observables 
	    _h_coscos[i]   = bookHisto1D(("coscos_cut_"+std::to_string(i)).c_str(),   20,  -1.,    1.); 
	    _h_cosphi[i]   = bookHisto1D(("cosphi_cut_"+std::to_string(i)).c_str(),   20,  -1.,    1.);
	    _h_deltaphi[i] = bookHisto1D(("deltaphi_cut_"+std::to_string(i)).c_str(), 20,   0.,    1.);
	    //agrohsje check kelly's implementation 
	    //_h_deltaphi[i] = bookHisto1D(("deltaphi_cut_"+std::to_string(i)).c_str(),  5,  -1.,    1.);
	    _h_deltaeta[i] = bookHisto1D(("deltaeta_cut_"+std::to_string(i)).c_str(), 20,   0.,    1.);
	    _h_c1[i]       = bookHisto1D(("c1_cut_"+std::to_string(i)).c_str(),       20,  -1.,    1.);
	    _h_c2[i]       = bookHisto1D(("c2_cut_"+std::to_string(i)).c_str(),       20,  -1.,    1.);
	    _h_c3[i]       = bookHisto1D(("c3_cut_"+std::to_string(i)).c_str(),       20,  -1.,    1.);
	    _h_c4[i]       = bookHisto1D(("c4_cut_"+std::to_string(i)).c_str(),       40,  -2.,    2.);
	    _h_c5[i]       = bookHisto1D(("c5_cut_"+std::to_string(i)).c_str(),       20,  -1.,    1.);
	    _h_c6[i]       = bookHisto1D(("c6_cut_"+std::to_string(i)).c_str(),       20,  -1.,    1.);
	    _h_b1_cpc[i]   = bookHisto1D(("b1cpc_cut_"+std::to_string(i)).c_str(),    20,  -1.,    1.);
	    _h_b2_cpc[i]   = bookHisto1D(("b2cpc_cut_"+std::to_string(i)).c_str(),    20,  -1.,    1.);
	    _h_b3_cpc[i]   = bookHisto1D(("b3cpc_cut_"+std::to_string(i)).c_str(),    20,  -1.,    1.);
	    _h_b1_cpv[i]   = bookHisto1D(("b1cpv_cut_"+std::to_string(i)).c_str(),    20,  -1.,    1.);
	    _h_b2_cpv[i]   = bookHisto1D(("b2cpv_cut_"+std::to_string(i)).c_str(),    20,  -1.,    1.);
	    _h_b3_cpv[i]   = bookHisto1D(("b3cpv_cut_"+std::to_string(i)).c_str(),    20,  -1.,    1.);
	    _h_c_nn[i]     = bookHisto1D(("c_nn_cut"+std::to_string(i)).c_str(),      20,  -1.,    1.);
	    _h_c_rr[i]     = bookHisto1D(("c_rr_cut"+std::to_string(i)).c_str(),      20,  -1.,    1.);
	    _h_c_kk[i]     = bookHisto1D(("c_kk_cut"+std::to_string(i)).c_str(),      20,  -1.,    1.);
	    _h_c_rkpkr[i]  = bookHisto1D(("c_rkpkr_cut"+std::to_string(i)).c_str(),   20,  -1.,    1.);
	    _h_c_nrprn[i]  = bookHisto1D(("c_nrprn_cut"+std::to_string(i)).c_str(),   20,  -1.,    1.);
	    _h_c_nkpkn[i]  = bookHisto1D(("c_nkpkn_cut"+std::to_string(i)).c_str(),   20,  -1.,    1.);
	    _h_c_rkmkr[i]  = bookHisto1D(("c_rkmkr_cut"+std::to_string(i)).c_str(),   20,  -1.,    1.);
	    _h_c_nrmrn[i]  = bookHisto1D(("c_nrmrn_cut"+std::to_string(i)).c_str(),   20,  -1.,    1.);
	    _h_c_nkmkn[i]  = bookHisto1D(("c_nkmkn_cut"+std::to_string(i)).c_str(),   20,  -1.,    1.);
	    // agrohsje 
	    _h_b_ntop[i]      = bookHisto1D(("b_ntop_cut"+std::to_string(i)).c_str(),      20,  -1.,    1.);
	    _h_b_ntbar[i]     = bookHisto1D(("b_ntbar_cut"+std::to_string(i)).c_str(),     20,  -1.,    1.);
	    _h_b_rtop[i]      = bookHisto1D(("b_rtop_cut"+std::to_string(i)).c_str(),      20,  -1.,    1.);
	    _h_b_rtbar[i]     = bookHisto1D(("b_rtbar_cut"+std::to_string(i)).c_str(),     20,  -1.,    1.);
	    _h_b_ktop[i]      = bookHisto1D(("b_ktop_cut"+std::to_string(i)).c_str(),      20,  -1.,    1.);
	    _h_b_ktbar[i]     = bookHisto1D(("b_ktbar_cut"+std::to_string(i)).c_str(),     20,  -1.,    1.);
	    _h_b_rstartop[i]  = bookHisto1D(("b_rstartop_cut"+std::to_string(i)).c_str(),  20,  -1.,    1.);
	    _h_b_rstartbar[i] = bookHisto1D(("b_rstartbar_cut"+std::to_string(i)).c_str(), 20,  -1.,    1.);
	    _h_b_kstartop[i]  = bookHisto1D(("b_kstartop_cut"+std::to_string(i)).c_str(),  20,  -1.,    1.);
	    _h_b_kstartbar[i] = bookHisto1D(("b_kstartbar_cut"+std::to_string(i)).c_str(), 20,  -1.,    1.);
	  }
	}
		
	///////////////////////////////////////////////////////////////////////////////////////////////

	void _fillhistos(int icut, const double weight){
	    
	  //agrohsje 
	  _h_met[icut]->fill(_pmet.pT(), weight);
	  // _h_mt2ll[icut]->fill(_mt2ll, weight);
	  // _h_cem[icut]->fill(_cem, weight);
	  // // add some debugging
	  if(_fmtop.size()!=2){
	    std::cout<<" Found "<< _fmtop.size() << " tops with status 62!" <<std::endl;
	    //throw std::runtime_error("There should be 2 tops!");   
	  }
	  // if(_fmphi.size()>1){
	  //   std::cout<<" Found "<< _fmphi.size() << " DM mediators with status 62 " <<std::endl;
	  //   throw std::runtime_error("There shouldn't be more than 1 DM mediator!");   
	  // }
	  // agrohsje remove check later 
	  if(_fmtop.size()==2){
	    // not possible for now as we only deal with 4-momenta 
	    //if(_fmtop[0].charge()*_fmtop[1].charge()<0){
	    for(unsigned i=0;i<_fmtop.size();i++) {
	      _h_mtop[icut]->fill(_fmtop[i].mass(),weight);
	      _h_etatop[icut]->fill(_fmtop[i].eta(),weight);
	      _h_phitop[icut]->fill(_fmtop[i].phi(ZERO_2PI),weight);
	      _h_pttop[icut]->fill(_fmtop[i].perp(),weight);
	    }
	    _h_mttbar[icut]->fill((_fmtop[0]+_fmtop[1]).mass(), weight);	  
	    _h_ptttbar[icut]->fill(std::log10((_fmtop[0]+_fmtop[1]).perp()), weight);	  
	  }
	  //   else 
	  //     std::cout<<" Found 2 tops but with same charge!" <<std::endl;
	  // }
	  // agrohsje 
	  if(_fmchi.size()==2){
	    _h_mphi[icut]->fill((_fmchi[0]+_fmchi[1]).mass(),weight);
	    _h_etaphi[icut]->fill((_fmchi[0]+_fmchi[1]).eta(),weight);
	    _h_ptphi[icut]->fill((_fmchi[0]+_fmchi[1]).perp(),weight);
	  }
	  //else{
	  //std::cout<<" Found "<<_fmchi.size()<<" DM particles instead of 2!" <<std::endl;
	  //}

	  // for(unsigned i=0;i<_fmchi.size();i++){
	  //   _h_mchi[icut]->fill(_fmchi[i].mass(),weight);
	  //   _h_etachi[icut]->fill(_fmchi[i].eta(),weight);
	  //   _h_phichi[icut]->fill(_fmchi[i].phi(ZERO_2PI),weight);
	  //   _h_ptchi[icut]->fill(_fmchi[i].perp(),weight);
	  // }
	  // for(unsigned i=0;i<_fmphi.size();i++) {
	  //   _h_mphi[icut]->fill(_fmphi[i].mass(),weight);
	  //   _h_etaphi[icut]->fill(_fmphi[i].eta(),weight);
	  //   _h_phiphi[icut]->fill(_fmphi[i].phi(ZERO_2PI),weight);
	  //   _h_ptphi[icut]->fill(_fmphi[i].perp(),weight);
	  // }


	  if(_recospinobservables){
	    //std::cout<<" check my " << _cosphi << " and kellys " << _deltaphi   << " implementation?" << std::endl; 
	    _h_coscos[icut]->fill(_coscos, weight);  
	    _h_cosphi[icut]->fill(_cosphi, weight);  
	    _h_deltaphi[icut]->fill(_deltaphi, weight);
	    _h_deltaeta[icut]->fill(_deltaeta, weight);
	    _h_c1[icut]->fill(_c1, weight);      
	    _h_c2[icut]->fill(_c2, weight);      
	    _h_c3[icut]->fill(_c3, weight);      
	    _h_c4[icut]->fill(_c4, weight);      
	    _h_c5[icut]->fill(_c5, weight);      
	    _h_c6[icut]->fill(_c6, weight);      
	    
	    // assume no CP violation 
	    _h_b1_cpc[icut]->fill( _b1_lp, weight);  
	    _h_b1_cpc[icut]->fill( _b1_lm, weight);  
	    _h_b2_cpc[icut]->fill( _b2_lp, weight);  
	    _h_b2_cpc[icut]->fill( _b2_lm, weight);  
	    _h_b3_cpc[icut]->fill( _b3_lp, weight);  
	    _h_b3_cpc[icut]->fill(-_b3_lm, weight);  
	    
	    //assume CP violation 
	    _h_b1_cpv[icut]->fill( _b1_lp, weight);  
	    _h_b1_cpv[icut]->fill(-_b1_lm, weight);  
	    _h_b2_cpv[icut]->fill( _b2_lp, weight);  
	    _h_b2_cpv[icut]->fill(-_b2_lm, weight);  
	    _h_b3_cpv[icut]->fill( _b3_lp, weight);  
	    _h_b3_cpv[icut]->fill( _b3_lm, weight);  
	   
	    _h_c_nn[icut]->fill(_c_nn, weight);
	    _h_c_rr[icut]->fill(_c_rr, weight);
	    _h_c_kk[icut]->fill(_c_kk, weight);
	    _h_c_rkpkr[icut]->fill(_c_rkpkr, weight);
	    _h_c_nrprn[icut]->fill(_c_nrprn, weight);
	    _h_c_nkpkn[icut]->fill(_c_nkpkn, weight);
	    _h_c_rkmkr[icut]->fill(_c_rkmkr, weight);
	    _h_c_nrmrn[icut]->fill(_c_nrmrn, weight);
	    _h_c_nkmkn[icut]->fill(_c_nkmkn, weight);
	    _h_b_ntop[icut]->fill(_b_ntop, weight);
	    _h_b_ntbar[icut]->fill(_b_ntbar, weight);
	    _h_b_rtop[icut]->fill(_b_rtop, weight);
	    _h_b_rtbar[icut]->fill(_b_rtbar, weight);
	    _h_b_ktop[icut]->fill(_b_ktop, weight);
	    _h_b_ktbar[icut]->fill(_b_ktbar, weight);
	    _h_b_rstartop[icut]->fill(_b_rstartop, weight);
	    _h_b_rstartbar[icut]->fill(_b_rstartbar, weight);
	    _h_b_kstartop[icut]->fill(_b_kstartop, weight);
	    _h_b_kstartbar[icut]->fill(_b_kstartbar, weight);

	  }
	}
      ///////////////////////////////////////////////////////////////////////////////////////////////
      
      bool _overlap() {
	bool overlap = false;
	for (unsigned i = 0; i < _alljets.size(); i++) {
	  const Jet& jet = _alljets[i];
	  foreach (const DressedLepton& el, _dressedelectrons) 
	    if (deltaR(jet, el) < 0.4) overlap = true;
	  foreach (const DressedLepton& mu, _dressedmuons) 
	    if (deltaR(jet, mu) < 0.4) overlap = true;
	  // agrohsje 
	  /*for (unsigned j = i+1; j < _alljets.size(); j++) {
	    const Jet& jet2 = _alljets[j];
		    if (deltaR(jet, jet2) < 0.5) overlap = true;
		    }*/
	}    
	return overlap;
      }
      
      ///////////////////////////////////////////////////////////////////////////////////////////////
	
	bool _passee(unsigned int& cutBits) {
	    // 1. Exactly 2 good electrons
	    cutBits += 1; if (_dressedelectrons.size() != 2) return false;
            // 2. No additional electrons passing the veto selection
	    cutBits += 1 << 1; if (_vetodressedelectrons.size() > 2) return false;
            // 3. No muons passing the veto selection
	    cutBits += 1 << 2; if (_vetodressedmuons.size() > 0) return false;
            // 4. Selected Leptons need to have opposite sign
	    cutBits += 1 << 3; if (_dressedelectrons[0].charge() == _dressedelectrons[1].charge()) return false;
	    double m_inv = (_dressedelectrons[0].momentum()+_dressedelectrons[1].momentum()).mass();
            // 5. Invariant dilepton mass > 20 GeV
	    cutBits += 1 << 4; if (m_inv < _mllmin*GeV) return false;
            // 6. Invariant dilepton mass not in Z-mass window
	    cutBits += 1 << 5; if (fabs(m_inv - 91.2*GeV) < _mllwindow*GeV) return false;
	    // 7. At least two good jets
	    cutBits += 1 << 6; if (_alljets.size() < _njets) return false;
            // 8. total neutrino pT > 60 GeV
	    cutBits += 1 << 7; if (_pmet.pT() <= _metmin*GeV) return false;
            // 9. At least one b-tagged jet
	    cutBits += 1 << 8; if (_bjets.size() < _nbtags) return false;
            // after all cuts
	    cutBits += 1 << 9;

            return true;
        }
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	bool _passmumu(unsigned int& cutBits) {

            // 1. Exactly 2 good muons
	    cutBits += 1; if (_dressedmuons.size() != 2) return false;
            // 2. No additional electrons passing the veto selection
	    cutBits += 1 << 1; if (_vetodressedmuons.size() > 2) return false;
            // 3. No muons passing the veto selection
	    cutBits += 1 << 2; if (_vetodressedelectrons.size() > 0) return false;
            // 4. Selected Leptons need to have opposite sign
	    cutBits += 1 << 3; if (_dressedmuons[0].charge() == _dressedmuons[1].charge()) return false;
            double m_inv = (_dressedmuons[0].momentum()+_dressedmuons[1].momentum()).mass();
            // 5. Invariant dilepton mass > 20 GeV
	    cutBits += 1 << 4; if (m_inv < _mllmin*GeV) return false;
            // 6. Invariant dilepton mass not in Z-mass window
	    cutBits += 1 << 5; if (fabs(m_inv - 91.2*GeV) < _mllwindow*GeV) return false;
            // 7. At least two good jets
	    cutBits += 1 << 6; if (_alljets.size() < _njets) return false;
            // 8. total neutrino pT > 60 GeV
	    cutBits += 1 << 7; if (_pmet.pT() <= _metmin*GeV) return false;
	    // 9. At least one b-tagged jet
	    cutBits += 1 << 8; if (_bjets.size() < _nbtags) return false;
            // after all cuts
	    cutBits += 1 << 9;
	    
            return true;
        }
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	bool _passemu(unsigned int& cutBits){
            
            // 1. Exactly 1 good muon
	    cutBits += 1; if (_dressedmuons.size()          != 1) return false;
            // 2. Exactly 1 good electron
	    cutBits += 1 << 1; if (_dressedelectrons.size() != 1) return false;
            // 3. No additional electrons passing the veto selection
	    cutBits += 1 << 2; if (_vetodressedmuons.size() > 1) return false;
            // 4. No muons passing the veto selection
	    cutBits += 1 << 3; if (_vetodressedelectrons.size() > 1) return false;
	    // 5. Selected Leptons need to have opposite sign
	    cutBits += 1 << 4; if (_dressedmuons[0].charge() == _dressedelectrons[0].charge()) return false;
            // 6. mll > 20 GeV
	    double m_inv = (_dressedmuons[0].momentum()+_dressedelectrons[0].momentum()).mass();
	    cutBits += 1 << 5; if ( m_inv <= _mllmin*GeV) return false;
            // 7. At least two good jets
	    cutBits += 1 << 6; if (_alljets.size() < _njets) return false;
            // 8. At least one b-tagged jet
	    cutBits += 1 << 7; if (_bjets.size() < _nbtags) return false;
	    // events after all cuts
	    cutBits += 1 << 8;
	    
            return true;
        }
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	void _reconstructspinobservables() {
	    
	    _recospinobservables = true;
	    //_toplepmomenta already sorted, first is top, second tbar 
	    if(_toplepmomenta.size()!=2) { _recospinobservables = false; return; }
	    
	    const FourMomentum top_lab(_toplepmomenta[0].first);
	    const FourMomentum lbar_lab(_toplepmomenta[0].second);
	    const FourMomentum tbar_lab(_toplepmomenta[1].first);
	    const FourMomentum lep_lab(_toplepmomenta[1].second);	    
	    
	    // const FourMomentum top_lab( 436.663253 , -81.089152 , 232.906579 , -316.382081 );
	    // const FourMomentum tbar_lab( 456.329093 , 6.416976 , -148.730683 , -395.370709 );
	    // const FourMomentum lep_lab( 66.830857 , 1.020382 , -59.641327 , -30.136928 );
	    // const FourMomentum lbar_lab( 140.992880 , -36.308579 , 60.157785 , -122.236330 );
	    


	    
	    // agrohsje stuff from james 
	    /*
	    const FourMomentum lP4_1 = lbar_lab;
	    const FourMomentum lP4_2 = lep_lab;
	    const FourMomentum llP4 =  lP4_1  +  lP4_2 ;
	    const FourMomentum t1P4 = top_lab;
	    const FourMomentum t2P4 = tbar_lab;
	    const FourMomentum ttP4 = t1P4 + t2P4;
	    const FourMomentum t1P4AtCM = LorentzTransform::mkFrameTransformFromBeta(ttP4.betaVec()).transform(t1P4);
	    const FourMomentum t2P4AtCM = LorentzTransform::mkFrameTransformFromBeta(ttP4.betaVec()).transform(t2P4);
	    const FourMomentum l1P4Att1 = LorentzTransform::mkFrameTransformFromBeta(t1P4.betaVec()).transform(lP4_1);
	    const FourMomentum l2P4Att2 = LorentzTransform::mkFrameTransformFromBeta(t2P4.betaVec()).transform(lP4_2);
	    double Cos_theta_1 = cos (    t1P4AtCM.vector3().angle(l1P4Att1.vector3()) );
	    double Cos_theta_2 = cos (    t2P4AtCM.vector3().angle(l2P4Att2.vector3()) );
	    //double Cos_theta1_Cos_theta2 = Cos_theta_1 * Cos_theta_2;
	    double Cos_theta_star = cos ( l1P4Att1.vector3().angle(l2P4Att2.vector3()) );
	    //std::cout<<" cos theta 1 " << Cos_theta_1 << " Cos_theta_2 " << Cos_theta_2 << " Cos_theta_star " << Cos_theta_star << std::endl;
	    // end of added stuff from james 
	    */
	    //agrohsje
	    /*
	    std::cout<<" top pT "<< top_lab.perp() 
		     <<" top eta "<< top_lab.eta() 
		     << " tbar pT  " << tbar_lab.perp() 		     
		     << " tbar eta  " << tbar_lab.eta() 
		     << std::endl;	    

	    std::cout<<" lep pT "<< lep_lab.perp() 
		     <<" lep eta "<< lep_lab.eta() 		     
		     << " lbar pT  " << lbar_lab.perp() 
		     << " lbar eta  " << lbar_lab.eta() 
		     << std::endl;	    
	    */
	    const FourMomentum ttbar_lab(top_lab+tbar_lab);
	    // boost tops and leptons to ttbar ZMF
	    LorentzTransform boosttottbar;
	    boosttottbar.setBetaVec(-ttbar_lab.boostVector());
	    // agrohsje debug boosting in madana
	    const FourMomentum top_ZMFttbar = boosttottbar.transform(top_lab);
	    const FourMomentum tbar_ZMFttbar = boosttottbar.transform(tbar_lab);
	    const FourMomentum lep_ZMFttbar = boosttottbar.transform(lep_lab);
	    const FourMomentum lbar_ZMFttbar = boosttottbar.transform(lbar_lab);

	    // std::cout
	    //   <<" top in ttbar " << top_ZMFttbar.vector3().unit().x() << "  "  << top_ZMFttbar.vector3().unit().y() << "  "  << top_ZMFttbar.vector3().unit().z() << " \n "  
	    //   <<" lbar in ttbar " << lbar_ZMFttbar.vector3().unit().x() << "  "  << lbar_ZMFttbar.vector3().unit().y() << "  "  << lbar_ZMFttbar.vector3().unit().z() << " \n "  
	    //   <<endl;

	    // boost lepton/antilepton to corresponding antitop/top  
	    LorentzTransform boosttotop;
            boosttotop.setBetaVec(-top_ZMFttbar.boostVector());
	    const FourMomentum lbar_ZMFtop = boosttotop.transform(lbar_ZMFttbar);
	    LorentzTransform boosttotbar;
            boosttotbar.setBetaVec(-tbar_ZMFttbar.boostVector());
	    const FourMomentum lep_ZMFtbar = boosttotbar.transform(lep_ZMFttbar);
	    

	    std::cout<<" after boosting top_ZMFttbar = ("<< top_ZMFttbar.E() << " , "<< top_ZMFttbar.px() << " , "<< top_ZMFttbar.py() << " , "<< top_ZMFttbar.pz() << " ) " << std::endl;
	    std::cout<<" after boosting tbar_ZMFttbar = ("<< tbar_ZMFttbar.E() << " , "<< tbar_ZMFttbar.px() << " , "<< tbar_ZMFttbar.py() << " , "<< tbar_ZMFttbar.pz() << " ) " << std::endl;
	    std::cout<<" after boosting lep_ZMFtbar = ("<< lep_ZMFtbar.E() << " , "<< lep_ZMFtbar.px() << " , "<< lep_ZMFtbar.py() << " , "<< lep_ZMFtbar.pz() << " ) " << std::endl;
	    std::cout<<" after boosting lbar_ZMFtop = ("<< lbar_ZMFtop.E() << " , "<< lbar_ZMFtop.px() << " , "<< lbar_ZMFtop.py() << " , "<< lbar_ZMFtop.pz() << " ) " << std::endl;
	    
	    // linear combinations of spin density observables 
	    // agrohsje debug james,mykola 
	    //std::cout<<" my angle cos 1 "<< cos(lep_ZMFtbar.angle(tbar_ZMFttbar)) << " cos 2 " << cos(lbar_ZMFtop.angle(top_ZMFttbar)) 
	    //	     <<" cos theta star "<< cos(lep_ZMFtbar.angle(lbar_ZMFtop)) << std::endl;
	    // cos cos heli basis 
	    _coscos = cos(lep_ZMFtbar.angle(tbar_ZMFttbar))*cos(lbar_ZMFtop.angle(top_ZMFttbar));
	    //agrohsje std::cout << " m_coscos is " << _coscos << std::endl;


	    // cos between leptons 
	    _cosphi = cos(lep_ZMFtbar.angle(lbar_ZMFtop));
	    
	    //agrohsje std::cout<<" m_cosphi is " << _cosphi << std::endl;
	    
	    // azimuthal difference between leptons 
	    _deltaphi=deltaPhi(lep_lab, lbar_lab)/PI;
	    
	    // agrohsje check kellys implementation 
	    // LorentzTransform boosttotoplab;
            // boosttotoplab.setBetaVec(-top_lab.boostVector());
	    // LorentzTransform boosttotbarlab;
	    // boosttotbarlab.setBetaVec(-tbar_lab.boostVector());
	    // const FourMomentum lep_ZMFtbarlab = boosttotbarlab.transform(lep_lab);
	    // const FourMomentum lbar_ZMFtoplab = boosttotoplab.transform(lbar_lab);
	    // _deltaphi=cos(lep_ZMFtbarlab.angle(lbar_ZMFtoplab));
	    // agrohsje end 

	    // pseudorapidity difference in lab frame  
	    //_deltaeta = abs(lep_lab.eta()-lbar_lab.eta()); 
	    _deltaeta = abs(tanh((lep_lab.eta()-lbar_lab.eta())/2.)); 
	    //agrohsje std::cout<<" m_deltaeta is " << _deltaeta << " " << std::endl;
	    
	    // initializing P
	    Vector3 pbone(0.,0.,1.); 
	    Vector3 kbone = top_ZMFttbar.vector3().unit();
	    Vector3 nbone = pbone.cross(kbone);
	    // agrohsje bug fix : 
	    //Vector3 lp= lbar_ZMFttbar.vector3().unit();
	    //Vector3 lm= lep_ZMFttbar.vector3().unit();
	    Vector3 lp= lbar_ZMFtop.vector3().unit();
	    Vector3 lm= lep_ZMFtbar.vector3().unit();
	    
	    // fill spin correlation part of density matrix 
	    _c1 = lp.dot(lm);
	    _c2 = pbone.dot(lp)*pbone.dot(lm);
	    _c3 = kbone.dot(lp)*kbone.dot(lm);
	    _c4 = pbone.dot(lp)*kbone.dot(lm)+kbone.dot(lp)*pbone.dot(lm);
	    _c5 = (lp.cross(lm)).dot(pbone);
	    _c6 = (lp.cross(lm)).dot(kbone);

	    // fill polarization part of density matrix 
	    _b1_lp = pbone.dot(lp);
	    _b1_lm = pbone.dot(lm);
	    _b2_lp = kbone.dot(lp);
	    _b2_lm = kbone.dot(lm);
	    _b3_lp = nbone.dot(lp);
	    _b3_lm = nbone.dot(lm);
	    //agrohsje debug 
	    //std::cout<<" compare the two "<< cos(lep_ZMFtbar.angle(tbar_ZMFttbar)) << "  ich  " << k.dot(lm) << std::endl;
	    
	    //agrohsje std::cout << " " << _c1 << " " << _c2 << " " << _c3 << " " << _c4 << " " << _c5 << " " << _c6 << std::endl;
	    //agrohsje std::cout << " " << _b1_lp << " " << _b1_lm << " " << _b2_lp << " " << _b2_lm << " " << _b3_lp << " " << _b3_lm << std::endl;
	    
	    // implement observables based on https://arxiv.org/pdf/1508.05271v2.pdf
	    // second base choice for spin correlation and polarization 
	    double yp = pbone.dot(kbone);
	    double rp = sqrt(1-yp*yp);
	    if(yp*yp>1)	
	      std::cout << "TTbarSpinDensityMatru::ReconstructSpinObservables rp not well defined : " << rp << std::endl;
	    
	    double signyp = (yp>0) ? 1. : -1.;
	    Vector3 pbtwo(0.,0.,1.); 
	    Vector3 rbtwo = 1/rp*(pbone-yp*kbone); 
	    Vector3 nbtwo = 1/rp*(pbone.cross(kbone));
	    
	    _c_nn    = lbar_ZMFtop.vector3().unit().dot(signyp*nbtwo) * lep_ZMFtbar.vector3().unit().dot(-1*signyp*nbtwo);
	    _c_rr    = lbar_ZMFtop.vector3().unit().dot(signyp*rbtwo) * lep_ZMFtbar.vector3().unit().dot(-1*signyp*rbtwo);
	    _c_kk    = lbar_ZMFtop.vector3().unit().dot(kbone) * lep_ZMFtbar.vector3().unit().dot(-1*kbone);
	    
	    _c_rkpkr = lbar_ZMFtop.vector3().unit().dot(signyp*rbtwo) * lep_ZMFtbar.vector3().unit().dot(-1*kbone) + 
	      lbar_ZMFtop.vector3().unit().dot(kbone) * lep_ZMFtbar.vector3().unit().dot(-1*signyp*rbtwo);
	    _c_nrprn = lbar_ZMFtop.vector3().unit().dot(signyp*nbtwo) * lep_ZMFtbar.vector3().unit().dot(-1*signyp*rbtwo) +
	      lbar_ZMFtop.vector3().unit().dot(signyp*rbtwo) * lep_ZMFtbar.vector3().unit().dot(-1*signyp*nbtwo);
	    _c_nkpkn = lbar_ZMFtop.vector3().unit().dot(signyp*nbtwo) * lep_ZMFtbar.vector3().unit().dot(-1*kbone) +
	      lbar_ZMFtop.vector3().unit().dot(kbone) * lep_ZMFtbar.vector3().unit().dot(-1*signyp*nbtwo);
	    
	    _c_rkmkr = lbar_ZMFtop.vector3().unit().dot(signyp*rbtwo) * lep_ZMFtbar.vector3().unit().dot(-1*kbone) -
	      lbar_ZMFtop.vector3().unit().dot(kbone) * lep_ZMFtbar.vector3().unit().dot(-1*signyp*rbtwo);
	    _c_nrmrn = lbar_ZMFtop.vector3().unit().dot(signyp*nbtwo) * lep_ZMFtbar.vector3().unit().dot(-1*signyp*rbtwo) -
	      lbar_ZMFtop.vector3().unit().dot(signyp*rbtwo) * lep_ZMFtbar.vector3().unit().dot(-1*signyp*nbtwo);
	    _c_nkmkn = lbar_ZMFtop.vector3().unit().dot(signyp*nbtwo) * lep_ZMFtbar.vector3().unit().dot(-1*kbone) -
	      lbar_ZMFtop.vector3().unit().dot(kbone) * lep_ZMFtbar.vector3().unit().dot(-1*signyp*nbtwo);
	    
	    _b_ntop  = lbar_ZMFtop.vector3().unit().dot(signyp*nbtwo); 
	    _b_ntbar = lep_ZMFtbar.vector3().unit().dot(-1*signyp*nbtwo);
	    
	    _b_rtop  = lbar_ZMFtop.vector3().unit().dot(signyp*rbtwo); 
	    _b_rtbar = lep_ZMFtbar.vector3().unit().dot(-1*signyp*rbtwo);
	    
	    _b_ktop  = lbar_ZMFtop.vector3().unit().dot(kbone); 
	    _b_ktbar = lep_ZMFtbar.vector3().unit().dot(-1*kbone);
	    
	    
	    // third base for spin polarization 
	    double signdrap = ( abs(top_lab.rapidity())>abs(tbar_lab.rapidity()) ) ? 1. : -1; 
	    Vector3 kbthree = signdrap * kbone;
	    Vector3 rbthree = signdrap * signyp * rbtwo;
	    
	    _b_rstartop  = lbar_ZMFtop.vector3().unit().dot(rbthree);
	    _b_rstartbar = lep_ZMFtbar.vector3().unit().dot(-1*rbthree);
	    
	    _b_kstartop  = lbar_ZMFtop.vector3().unit().dot(kbthree); 
	    _b_kstartbar = lep_ZMFtbar.vector3().unit().dot(-1*kbthree);
	    
	    return;
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	bool _islastcopy(const GenParticle& gp){
	    
	    const GenParticle* pcopy = &gp;
	    std::unordered_set<const GenParticle*> dupcheck;	    
	    dupcheck.clear();
	    //std::cout<<" original particle is " << pcopy->pdg_id() << " " << pcopy->status() << " " << pcopy->momentum().perp() << std::endl;  
	    while (_nextcopy(*pcopy)) {
		dupcheck.insert(pcopy);
		pcopy = _nextcopy(*pcopy);
		//std::cout<<" pcopy " << pcopy->pdg_id() << " " << pcopy->status() << " " << pcopy->momentum().perp() << " xcheck " <<  dupcheck.count(pcopy)  << std::endl;  
		if (dupcheck.count(pcopy)) return false;
	    }
	    return &gp == pcopy;	
	}
	///////////////////////////////////////////////////////////////////////////////////////////////
	const GenParticle* _nextcopy(const GenParticle& gp){
	    
	    unsigned ndau = gp.end_vertex() ? gp.end_vertex()->particles_out_size() : 0 ;
	    for(unsigned idau = 0; idau<ndau; ++idau) {
		const GenParticle* dau(*(gp.end_vertex()->particles_out_const_begin() + idau));
		if(dau->pdg_id()==gp.pdg_id()) return dau;
	    }
	    return 0;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////
	
    private:
      
      /// @name Objects that are used by the event selection decisions
      //@{
      unsigned _dilepton;
      unsigned _selected;
      
      bool _recospinobservables;
      
      vector<DressedLepton> _dressedelectrons;
      vector<DressedLepton> _vetodressedelectrons;
      vector<DressedLepton> _dressedmuons;
      vector<DressedLepton> _vetodressedmuons;
      Particles _neutrinos;
      Particles _chis;
      Jets _bjets, _ljets, _alljets;
      FourMomentum _pmet;
      vector<GenParticle*> _selectedgenparticles;

      // cut and cutting related values
      unsigned _njets, _nbtags; 
      double _mllmin, _mllwindow, _metmin;
      
      // object four momenta and event variables
      vector<FourMomentum> _fmtop, _fmchi, _fmphi;
      double _mt2ll, _cem;      

      // spin observables  
      double _coscos, _cosphi, _deltaphi, _deltaeta;
      double _c1, _c2, _c3, _c4, _c5, _c6;
      double _b1_lp, _b1_lm, _b2_lp, _b2_lm, _b3_lp, _b3_lm; 
      double _c_nn, _c_rr, _c_kk;
      double _c_rkpkr, _c_nrprn , _c_nkpkn, _c_rkmkr, _c_nrmrn, _c_nkmkn; 
      double _b_ntop, _b_ntbar, _b_rtop, _b_rtbar, _b_ktop, _b_ktbar; 
      double _b_rstartop, _b_rstartbar, _b_kstartop,  _b_kstartbar;
      
      std::vector<std::pair<Particle,Particle> > _toplep; 
      std::vector<std::pair<FourMomentum,FourMomentum> > _toplepmomenta; 
      
      // all histo arrays should have the same length corresponding to the different cuts 
      // event properties 
      std::array<Histo1DPtr,  3>  _h_met;
      std::array<Histo1DPtr,  3>  _h_mt2ll;
      std::array<Histo1DPtr,  3>  _h_cem; 
      
      // top properties 
      std::array<Histo1DPtr, 3>  _h_mttbar;
      std::array<Histo1DPtr, 3>  _h_ptttbar;
      std::array<Histo1DPtr, 3>  _h_mtop;
      std::array<Histo1DPtr, 3>  _h_pttop;
      std::array<Histo1DPtr, 3>  _h_etatop;
      std::array<Histo1DPtr, 3>  _h_phitop;
      //agrohsje 
      // DM properties 
      std::array<Histo1DPtr, 3>  _h_mchi;
      std::array<Histo1DPtr, 3>  _h_ptchi;
      std::array<Histo1DPtr, 3>  _h_etachi;
      std::array<Histo1DPtr, 3>  _h_phichi;
      // mediator properties 
      std::array<Histo1DPtr, 3>  _h_mphi;
      std::array<Histo1DPtr, 3>  _h_ptphi;
      std::array<Histo1DPtr, 3>  _h_etaphi;
      std::array<Histo1DPtr, 3>  _h_phiphi;
      
      // spin density observables 
      std::array<Histo1DPtr, 3>  _h_coscos;
      std::array<Histo1DPtr, 3>  _h_cosphi;
      std::array<Histo1DPtr, 3>  _h_deltaphi;
      std::array<Histo1DPtr, 3>  _h_deltaeta;
      std::array<Histo1DPtr, 3>  _h_c1;
      std::array<Histo1DPtr, 3>  _h_c2;
      std::array<Histo1DPtr, 3>  _h_c3;
      std::array<Histo1DPtr, 3>  _h_c4;
      std::array<Histo1DPtr, 3>  _h_c5;
      std::array<Histo1DPtr, 3>  _h_c6;
      std::array<Histo1DPtr, 3>  _h_b1_cpc;
      std::array<Histo1DPtr, 3>  _h_b2_cpc;
      std::array<Histo1DPtr, 3>  _h_b3_cpc;
      std::array<Histo1DPtr, 3>  _h_b1_cpv;
      std::array<Histo1DPtr, 3>  _h_b2_cpv;
      std::array<Histo1DPtr, 3>  _h_b3_cpv;
      std::array<Histo1DPtr, 3>  _h_c_nn;
      std::array<Histo1DPtr, 3>  _h_c_rr;
      std::array<Histo1DPtr, 3>  _h_c_kk;
      std::array<Histo1DPtr, 3>  _h_c_rkpkr;
      std::array<Histo1DPtr, 3>  _h_c_nrprn;
      std::array<Histo1DPtr, 3>  _h_c_nkpkn;
      std::array<Histo1DPtr, 3>  _h_c_rkmkr;
      std::array<Histo1DPtr, 3>  _h_c_nrmrn;
      std::array<Histo1DPtr, 3>  _h_c_nkmkn;
      std::array<Histo1DPtr, 3>  _h_b_ntop;
      std::array<Histo1DPtr, 3>  _h_b_ntbar;
      std::array<Histo1DPtr, 3>  _h_b_rtop;
      std::array<Histo1DPtr, 3>  _h_b_rtbar;
      std::array<Histo1DPtr, 3>  _h_b_ktop;
      std::array<Histo1DPtr, 3>  _h_b_ktbar;
      std::array<Histo1DPtr, 3>  _h_b_rstartop;
      std::array<Histo1DPtr, 3>  _h_b_rstartbar;
      std::array<Histo1DPtr, 3>  _h_b_kstartop;
      std::array<Histo1DPtr, 3>  _h_b_kstartbar;

      //std::array<Histo1DPtr, 1>  _h_debug;
      
      //@}
      
    };
  
    

    // Declare the class as a hook for the plugin system
    DECLARE_RIVET_PLUGIN(TTbarSpinDensityMatrix);
    
}