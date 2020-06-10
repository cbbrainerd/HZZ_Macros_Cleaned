std::vector<muon> get_loose_muons(NewNtuple* nt) {
    std::vector<muon> loose_muons;
    int year=nt->year;
    float bdt_low,bdt_high;
    switch(year) {
        case 2016:
            bdt_low=0.8847169876098633;
            bdt_high=-0.19389629721641488;
            break;
        case 2017:
            bdt_low=0.883555161952972;
            bdt_high=-0.3830992293357821;
            break;
        case 2018:
            bdt_low=0.9506129026412962;
            bdt_high=-0.3629065185785282;
            break;
        default:
            std::cerr << "Invalid year " << year << std::endl;
            exit(1);
    }
    std::size_t number_of_muons=nt->RECOMU_PT->size();
    for(std::size_t i;i!=number_of_muons;++i) {
        double pt=nt->RECOMU_PT[i];
        //All muons must be at least loose to be considered
        if(!(
            pt > 5 && 
            fabs(nt->RECOMU_ETA[i]) < 2.4 && 
            nt->RECOMU_mubesttrkDxy[i] < 0.5 &&
            nt->RECOMU_mubesttrkDz[i] < 1.0 &&
            (nt->RECOMU_isGlobalMu[i] || (nt->RECOMU_isTrackerMu[i] && nt->RECOMU_numberOfMatches[i] > 0)) &&
            nt->RECOMU_mubesttrkType[i] != 2
        )) continue;
        float bdt=nt->RECOMU_BDT_Id[i];
        bool tight=
            (pt <= 10 && bdt > bdt_low) ||
            (pt >  10 && bdt > bdt_high) ||
            (pt > 200 && nt->RECOMU_isTrackerHighPtMu[i]);
        bool SIP=fabs(RECOMU_SIP[i] < 4
            
    }
}

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
        return fabs(SIP) < 4;
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

