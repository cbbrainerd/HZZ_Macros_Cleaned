#include <limits>
#include "particle.h"

enum signal_type {
    s_none, s_4e, s_4mu, s_2e2mu  
};

#define fill_def_mu(var,def) (if(i<=RECOMU_ ## var->size()) { var=RECOMU_ ## var->at(i); } else { var=def; })
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

void analyze() {

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
}
