#include <limits>

enum signal_type {
    s_none, s_4e, s_4mu, s_2e2mu  
};

#define fill_def_mu(var,def) (if(i<=RECOMU_ ## var->size()) { var=RECOMU_ ## var->at(i); } else { var=def; })

//Get all the information from recomu that is actually used in the analysis

//These structs should be defined in the header so they can access members of the analyzer directly
struct muon {
    bool RECOMU_isPFMu;
    bool RECOMU_isGlobalMu;
    bool RECOMU_isTrackerMu;
    bool RECOMU_isTrackerHighPtMu;
    float RECOMU_PT;
    float RECOMU_ETA;
    float RECOMU_THETA;
    float RECOMU_PHI;
    float RECOMU_MASS;
    float RECOMU_CHARGE;
    float RECOMU_TRACKISO;
    double RECOMU_PFchHad;
    double RECOMU_PFneuHad;
    double RECOMU_PFphoton;
    double RECOMU_PFPUchAllPart;
    float RECOMU_SIP;
    float RECOMU_IP;
    float RECOMU_IPERROR;
    int RECOMU_numberOfMatches;
    int RECOMU_trkmuArbitration;
    float RECOMU_mutrkNMuonHits;
    float RECOMU_muInnertrktrackerLayersWithMeasurement;
    float RECOMU_muInnertrkPTError;
    float RECOMU_muInnertrkNPixHits;
    int RECOMU_mubesttrkType;
    float RECOMU_mubesttrkDxy;
    float RECOMU_mubesttrkDz;
    float RECOMU_mubesttrkPTError;
    bool loose;
    bool tight;
    bool loose_plus_SIP;
    bool tight_plus_SIP;
    muon(int i) const {
        //These branches are always filled
        isPFMu=RECOMU_isPFMu->at(i);
        isGlobalMu=RECOMU_isGlobalMu->at(i);
        isTrackerMu=RECOMU_isTrackerMu->at(i);
        isTrackerHighPtMu=RECOMU_isTrackerHighPtMu->at(i);
        PT=RECOMU_PT->at(i);
        ETA=RECOMU_ETA->at(i);
        THETA=RECOMU_THETA->at(i);
        PHI=RECOMU_PHI->at(i);
        MASS=RECOMU_MASS->at(i);
        CHARGE=RECOMU_CHARGE->at(i);
        TRACKISO=RECOMU_TRACKISO->at(i);
        PFchHad=RECOMU_PFchHad->at(i);
        PFneuHad=RECOMU_PFneuHad->at(i);
        PFphoton=RECOMU_PFphoton->at(i);
        PFPUchAllPart=RECOMU_PFPUchAllPart->at(i);
        SIP=RECOMU_SIP->at(i);
        IP=RECOMU_IP->at(i);
        IPERROR=RECOMU_IPERROR->at(i);
        numberOfMatches=RECOMU_numberOfMatches->at(i);
        trkmuArbitration=RECOMU_trkmuArbitration->at(i);
        muInnertrkPTError=RECOMU_muInnertrkPTError->at(i);
        muInnertrkNPixHits=RECOMU_muInnertrkNPixHits->at(i);
        //These branches are always filled, or not, together.
        if(i<=mutrkNMuonHits->size()) {
            mutrkNMuonHits=RECOMU_mutrkNMuonHits->at(i);
            muInnertrktrackerLayersWithMeasurement=RECOMU_muInnertrktrackerLayersWithMeasurement->at(i);
        } else {
            mutrkNMuonHits=-999;
        }   
        if(i<=RECOMU_mubesttrkType->size()) {
            mubesttrkType=RECOMU_mubesttrkType->at(i);
            mubesttrkDxy=RECOMU_mubesttrkDxy->at(i);
            mubesttrkDz=RECOMU_mubesttrkDz->at(i);
            mubesttrkPTError=RECOMU_mubesttrkPTError->at(i);
        } else {
            mubesttrkType=-999;
            mubesttrkDxy=-999;
            mubesttrkDz=-999;
            mubesttrkPTError=-999;
        }
        loose=is_loose();
        tight=is_tight();
        loose_plus_SIP=loose && is_SIP();
        tight_plus_SIP=tight && is_SIP();
    }
private:
    bool is_loose() const {
        return (isGlobalMu || (isTrackerMu && numberOfMatches>0)) &&
        mubesttrkType!=2 &&
        PT > 5. &&
        ETA < 2.4 &&
        fabs(mubesttrkDz) < 1.;
    }
    bool is_tight() const {
        return is_loose() &&
               (isPFMu || (isTrackerHighPtMu && PT > 200.)) &&
               fabs(mubesttrkDxy) < .5;
    }
    bool is_SIP() const {
        return fabs(RECOMU_SIP) < 4.;
    }
};

struct electron {
    float E;
    float PT;
    float PTError;
    float P;
    float ETA;
    float PHI;
    float CHARGE;
    float ID;
    int isGap;
    double PFchHad;
    double PFneuHad;
    double PFphoton;
    double PFPUchAllPart;
    float SIP;
    float IP;
    float IPERROR;
    double mvaNonTrigV0;
    electron(int i) {
        E=RECOELE_E->at(i);
        PT=RECOELE_PT->at(i);
        PTError=RECOELE_PTError->at(i);
        P=RECOELE_P->at(i);
        ETA=RECOELE_ETA->at(i);
        PHI=RECOELE_PHI->at(i);
        CHARGE=RECOELE_CHARGE->at(i);
        ID=RECOELE_ID->at(i);
        //Only filled if true. Seems to be only branch that isn't always filled among used branches
        if(i<=isGap->size()) {
            isGap=RECOELE_isGap->at(i);
        } else {
            isGap=0;
        }
        PFchHad=RECOELE_PFchHad->at(i);
        PFneuHad=RECOELE_PFneuHad->at(i);
        PFphoton=RECOELE_PFphoton->at(i);
        PFPUchAllPart=RECOELE_PFPUchAllPart->at(i);
        SIP=RECOELE_SIP->at(i);
        IP=RECOELE_IP->at(i);
        IPERROR=RECOELE_IPERROR->at(i);
        mvaNonTrigV0=RECOELE_mvaNonTrigV0->at(i);
        loose=is_loose();
        tight=is_tight();
        loose_and_sip=loose && is_SIP();
    }
private:
    bool is_loose() const {
        return PT > 7. &&
               fabs(ETA) < 2.5 &&
               fabs(RECOELE_gsftrack_dxy) < .5 &&
               fabs(RECOELE_gsftrack_dz) < 1.;
    }
    bool is_tight() const {
        return is_loose() && ID==1; //BDT used for tight lepton (RECOELE_ID)
    }
    bool is_SIP() const {
        return SIP < 4;
    }
};

struct lepton {
    lepton() : electron_(nullptr), is_a_muon(-1) {}
    reset(muon* mu) {
        muon_=mu;
        is_a_muon=1;
    }
    reset(electron* ele) {
        electron_=ele;
        is_a_muon=0;
    }
private:
    union {
        muon* muon_;
        electron* electron_;
    };
    int is_a_muon; //-1 for empty, 0 for electron, 1 for muon
};

struct photon {
//Not filled in ntuple, so no sense reading them
    //double RECOPFPHOT_PFchHad;
    //double RECOPFPHOT_PFneuHad;
    //double RECOPFPHOT_PFphoton;
    //double RECOPFPHOT_PFPUchAllPart;
    double RECOPFPHOT_PFX_rho;
    float RECOPFPHOT_PT;
    //float RECOPFPHOT_PTError;
    float RECOPFPHOT_ETA;
    float RECOPFPHOT_PHI;
    //double RECOPFPHOT_PT_uncorr;
    photon(int i) {
        //RECOPFPHOT_PFchHad=RECOPFPHOT_PFchHad->at(i);
        //RECOPFPHOT_PFneuHad=RECOPFPHOT_PFneuHad->at(i);
        //RECOPFPHOT_PFphoton=RECOPFPHOT_PFphoton->at(i);
        //RECOPFPHOT_PFPUchAllPart=RECOPFPHOT_PFPUchAllPart->at(i);
        RECOPFPHOT_PFX_rho=RECOPFPHOT_PFX_rho->at(i);
        RECOPFPHOT_PT=RECOPFPHOT_PT->at(i);
        //RECOPFPHOT_PTError=RECOPFPHOT_PTError->at(i);
        RECOPFPHOT_ETA=RECOPFPHOT_ETA->at(i);
        RECOPFPHOT_PHI=RECOPFPHOT_PHI->at(i);
        //RECOPFPHOT_PT_uncorr=RECOPFPHOT_PT_uncorr->at(i);
    }
};

signal_type get_signal_type() {
      //MC_PDGID contains MC truth particle info on:
      //MC_PDGID[0]  the Higgs 25
      //MC_PDGID[1]  the Z1    23
      //MC_PDGID[2]  the Z2    23
      //MC_PDGID[3]  the l1+   mu = 13   e = 11
      //MC_PDGID[4]  the l2-   mu = -13  e = -11
      //MC_PDGID[5]  the l3+   mu = 13   e = 11
      //MC_PDGID[6]  the l4-   mu = -13  e = -11

    int lepton_1=safeAccess(MC_PDGID)[3], lepton_2=safeAccess(MC_PDGID)[4], lepton_3=safeAccess(MC_PDGID)[5], lepton_4=safeAccess(MC_PDGID)[4];

    if(lepton_1!=-lepton_2) return s_none;
    if(lepton_3!=-lepton_4) return s_none;
    if(lepton_1==11 && lepton_2==11) return s_4e;
    if((lepton_1==13 && lepton_2==11) || (lepton_1==11 && lepton_2==13)) return s_2e2mu;
    if(lepton_1==13 && lepton_2==13) return s_4mu;
}

//In event loop
//Step 0:
++N_0; 
N_0_w+=newweight;
hweight->Fill(newweight);

//Step 1:
event_signal_type=get_signal_type();
switch(event_signal_type) {
    case s_4e:
        ++N_1_4e;
        N_1_4e_w+=newweight;
        break;
    case s_4mu:
        ++N_1_4mu;
        N_1_4mu_w+=newweight;
        break;
    case s_2e2mu:
        ++N_1_2e2mu;
        N_1_2e2mu_w+=newweight;
        break;
    case s_none:
        ++N_1_none;
        N_1_none+=newweight;
}

bool isInAcceptance=true;
for(int i=3;i<7;++i) {
    if(
        //Muon eta < 2.5
        (std::abs(safeAccess(MC_PDGID)[i])==13 && std::abs(safeAccess(MC_ETA)[i]) >= 2.5) ||
        //Electron eta < 2.4
        (std::abs(safeAccess(MC_PDGID)[i])==11 && std::abs(safeAccess(MC_ETA)[i]) >= 2.4)
    ) {
        isInAcceptance=false;
        break;
    }
}

if(isInAcceptance) {
    ++N_2;
    N_2+=newweight;
    switch(event_signal_type) {
        case s_4e:
            ++N_2_4e;
            N_2_4e_w+=newweight;
            break;
        case s_4mu:
            ++N_2_4mu;
            N_2_4mu_w+=newweight;
            break;
        case s_2e2mu:
            ++N_2_2e2mu;
            N_2_2e2mu_w+=newweight;
            break;
        case s_none:
            ++N_2_none;
            N_2_none+=newweight;
            break;
    }
}

//lepton IDs. Also fill std::vector of leptons so we don't have to access through RECO... branches
std::vector<muon> muons;
for(std::size_t i=0;i!=RECOMU_PT->size();++i) {
    muons.emplace_back(i);
}

//Vector of pointers quality vectors to avoid copying
std::vector<muon*> loose_muons;
std::vector<muon*> tight_muons;
std::vector<muon*> loose_plus_SIP_muons;
std::vector<muon*> tight_plus_SIP_muons;
for(auto& mu : muons) {
    if(mu.loose) loose_muons.push_back(&mu);
    if(mu.tight) tight_muons.push_back(&mu);
    if(mu.loose_plus_SIP) loose_plus_SIP_muons.push_back(&mu);
    if(mu.tight_plus_SIP) tight_plus_SIP_muons.push_back(&mu);
}

std::vector<electron> electrons;
for(std::size_t i=0;i!=RECOELE_PT->size();++i) {
    electron ele(i);
    //Cross cleaning
    //Electrons are removed from consideration if they are too close to a tight+SIP muon
    bool cross_clean=true;
    for(muon* mu : tight_plus_SIP_muons) {
        double deltaR=sqrt(pow(DELTAPHI(mu->PHI,ele.PHI),2)+pow(mu->ETA-ele.ETA,2));
        if(deltaR <= .05) {
            cross_clean=false;
            break;
        }
    }
    if(cross_clean) electrons.push_back(ele);
}

std::vector<electron*> loose_electrons;
std::vector<electron*> loose_and_SIP_electrons;
std::vector<electron*> tight_electrons;
for(auto &ele : electrons) {
    if(ele.loose) loose_electrons.push_back(&ele);
    if(ele.loose_and_SIP) loose_and_SIP_electrons.push_back(&ele);
    if(ele.tight) loose_electrons.push_back(&ele);
}

std::vector<photon> photons;
std::vector<lepton> closest_lepton_to_photon;
std::vector<double> min_deltaR;
for(std::size_t i=0;i!=RECOPFPHOT_PT->size();++i) {
    photon phot(i);
    bool clean=1;
    for(electron* ele : loose_and_SIP_electrons) {
        double deltaPhi=fabs(DELTAPHI(phot.PHI,ele->scl_Phi));
        double deltaEta=fabs(phot.ETA-ele->scl_Eta);
        double deltaR=sqrt(pow(deltaPhi,2)+pow(deltaEta,2));
        if((deltaPhi < 2 && deltaEta < .05) || deltaR <= .15) {
            clean=0;
            break;
        }
    }
    if(clean) photons.push_back(phot);
    //Find closest lepton (loose), subject to deltaR/photonPT**2 < .12 and deltaR < .5
    bool min_is_muon=1;
    double minDeltaR=.5;
    lepton min;
    for(muon* mu : loose_muons) {
        double deltaR=sqrt(pow(DELTAPHI(phot.PHI,mu->PHI),2)+pow(phot.ETA-mu->ETA,2));
        if(deltaR/pow(phot.PT,2)<0.012 && deltaR<minDeltaR) {
            minDeltaR=deltaR;
            min.reset(mu);
        }
    }
    //Should this use supercluster eta, etc?
    for(electron* ele : loose_electrons) {
        double deltaR=sqrt(pow(DELTAPHI(phot.PHI,ele->PHI),2)+pow(phot.ETA-ele->ETA,2));
        if(deltaR/pow(phot.PT,2) <.012 && deltaR<minDeltaR) {
            minDeltaR=deltaR;
            min.reset(ele);
        }
    }
    min_deltaR.push_back(minDeltaR);
}
