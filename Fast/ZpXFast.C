#include <limits>
#include "ZpXanalyzer.h"
#include "Math/Vector4D.h"
#include <unordered_map>
#include <algorithm>

typedef ROOT::Math::PtEtaPhiMVector four_vector;
float Z_nominal_mass=91.188;

void ZpXanalyzer::Loop() {
   if (fChain == 0) return;

   for (Long64_t jentry=0; jentry < fChain->GetEntriesFast() ;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);
      analyze();
      // if (Cut(ientry) < 0) continue;
   }
}

//Lepton ID
//isMuon true: muon, false: lepton
//isTight true: tight ID (sip, iso, etc) false: loose + SIP
void ZpXanalyzer::analyze() {

    //lepton IDs. Also fill std::vector of leptons so we don't have to access through RECO... branches
    std::vector<muon> muons;
    for(std::size_t i=0;i!=RECOMU_PT->size();++i) {
        muons.emplace_back(this,i);
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
        electron ele(this,i);
        //Cross cleaning
        //Electrons are removed from consideration if they are too close to a tight+SIP muon
        bool cross_clean=true;
        for(muon* mu : tight_plus_SIP_muons) {
            double deltaRsq=(pow(DELTAPHI(mu->PHI,ele.PHI),2)+pow(mu->ETA-ele.ETA,2));
            if(deltaRsq <= 0.0025) {
                cross_clean=false;
                break;
            }
        }
        if(cross_clean) electrons.push_back(ele);
    }
    std::vector<electron*> loose_electrons;
    std::vector<electron*> loose_plus_SIP_electrons;
    std::vector<electron*> tight_electrons;
    std::vector<electron*> tight_plus_SIP_electrons;
    for(auto &ele : electrons) {
        if(ele.loose) loose_electrons.push_back(&ele);
        if(ele.loose_plus_SIP) loose_plus_SIP_electrons.push_back(&ele);
        if(ele.tight) loose_electrons.push_back(&ele);
        if(ele.loose_plus_SIP && ele.tight) tight_plus_SIP_electrons.push_back(&ele);
    }
    
    //If we don't have exactly 3 loose+SIP leptons, we can stop now
    if(loose_plus_SIP_electrons.size()+loose_plus_SIP_muons.size()!=3) return;
    //bkg_type is equivalent to number of electrons in event
    bkg_type type=(bkg_type)loose_plus_SIP_electrons.size();
    //We can also stop if we don't have an OSSF pair
    int number_tight=-1;
    int charge[3];
    switch(type) {
        case b_3mu:
        case b_1e2mu:
            number_tight=tight_plus_SIP_muons.size();
            charge[0]=tight_plus_SIP_muons[0]->CHARGE;
            charge[1]=tight_plus_SIP_muons[1]->CHARGE;
            if(number_tight>2) charge[2]=tight_plus_SIP_muons[2]->CHARGE;
            break;
        case b_2e1mu:
        case b_3e:
            number_tight=tight_plus_SIP_electrons.size();
            charge[0]=tight_plus_SIP_electrons[0]->CHARGE;
            charge[1]=tight_plus_SIP_electrons[1]->CHARGE;
            if(number_tight>2) charge[2]=tight_plus_SIP_electrons[2]->CHARGE;
    }
    //Need at least two leptons of opposite sign
    if(number_tight<2) return;
    if(charge[0]==charge[1] && (number_tight==2 || charge[0]==charge[2])) return;
    //To store fsr photons for each lepton
    std::unordered_map<void*,photon> fsr_photons;
    for(std::size_t i=0;i!=RECOPFPHOT_PT->size();++i) {
        photon phot(this,i);
        bool clean=1;
        for(electron* ele : loose_plus_SIP_electrons) {
            double deltaPhi=fabs(DELTAPHI(phot.PHI,ele->scl_Phi));
            double deltaEta=fabs(phot.ETA-ele->scl_Eta);
            double deltaR=sqrt(pow(deltaPhi,2)+pow(deltaEta,2));
            if((deltaPhi < 2 && deltaEta < .05) || deltaR <= .15) {
                clean=0;
                break;
            }
        }
        //Only photons that match a lepton are kept
//        if(clean) photons.push_back(phot);
        //Find closest lepton (loose+SIP), subject to deltaR/photonPT**2 < .12 and deltaR < .5
        double minDeltaR=.5;
        double minDeltaRoverPtSq;
        //lepton min{};
        void* nearest_lepton=nullptr;
        //Map FSR photon to each lepton
        for(muon* mu : loose_plus_SIP_muons) {
            double deltaR=sqrt(pow(DELTAPHI(phot.PHI,mu->PHI),2)+pow(phot.ETA-mu->ETA,2));
            double deltaRoverPtSq=deltaR/pow(phot.PT,2);
            if(deltaRoverPtSq<0.012 && deltaR<minDeltaR) {
                minDeltaR=deltaR;
                //Save the deltaRoverPtSq to the photon for later
                phot.deltaR=deltaR;
                phot.deltaRoverPtSq=deltaRoverPtSq;
                nearest_lepton=(void*)mu;
                //min.reset(mu);
            }
        }
        //Should this use supercluster eta, etc?
        for(electron* ele : loose_plus_SIP_electrons) {
            double deltaR=sqrt(pow(DELTAPHI(phot.PHI,ele->PHI),2)+pow(phot.ETA-ele->ETA,2));
            double deltaRoverPtSq=deltaR/pow(phot.PT,2);
            if(deltaRoverPtSq<.012 && deltaR<minDeltaR) {
                minDeltaR=deltaR;
                phot.deltaR=deltaR;
                phot.deltaRoverPtSq=deltaRoverPtSq;
                nearest_lepton=(void*)ele;
                //min.reset(ele);
            }
        }
        //Add photon to the matched lepton
        //First check if the given lepton is already matched to another photon
        auto prematched=fsr_photons.find(nearest_lepton);
        if(prematched!=fsr_photons.end()) {
            //One already matched, so pick the one with lowest deltaR/pt^2
            photon &old=prematched->second;
            if(phot.deltaRoverPtSq < old.deltaRoverPtSq) {
                //Replace the fsr photon
                prematched->second=phot;
            }
        } else { //First photon match, so keep it
            fsr_photons.insert(std::make_pair(nearest_lepton,phot));
        }
    }
    //Compute isolation for all muons of interest, subtracting fsr photon if it exists
    //Iso=(PFchHad+max(0.,PFneuHad+PFphoton-.5PFPUchAllPart))/PT
    //Where PFphoton subtracts off the fsr photon if it exists within the cone .01 < deltaR <= .3
    for(muon* mu : loose_plus_SIP_muons) {
        auto fsr_photon=fsr_photons.find((void*)mu);
        //This should be the same as the one stored in the tree, but we can double check
        mu->PFX_dB_prefsr=(mu->PFchHad+std::max(0.,mu->PFneuHad+mu->PFphoton-.5*mu->PFPUchAllPart))/mu->PT;
        //By default don't subtract fsr
        mu->PFX_dB_fsr=mu->PFX_dB_prefsr;
        if(fsr_photon!=fsr_photons.end()) {
            photon const& phot=fsr_photon->second;
            //Muon iso cone
            if (phot.deltaR < .3 && phot.deltaR > .01) {
                mu->PFX_dB_fsr=(mu->PFchHad+std::max(0.,mu->PFneuHad+(mu->PFphoton-phot.PT)-.5*mu->PFPUchAllPart))/mu->PT;
            }
        }
    }
    for(electron* ele : loose_plus_SIP_electrons) {
        auto fsr_photon=fsr_photons.find((void*)ele);
        //This should be the same as the one stored in the tree, but we can double check
        ele->PFX_dB_prefsr=(ele->PFchHad+std::max(0.,ele->PFneuHad+ele->PFphoton-.5*ele->PFPUchAllPart))/ele->PT;
        //By default don't subtract fsr
        ele->PFX_dB_fsr=ele->PFX_dB_prefsr;
        if(fsr_photon!=fsr_photons.end()) {
            photon const& phot=fsr_photon->second;
            //Muon iso cone
            if (phot.deltaR < .3 && phot.deltaR > .01) {
                ele->PFX_dB_fsr=(ele->PFchHad+std::max(0.,ele->PFneuHad+(ele->PFphoton-phot.PT)-.5*ele->PFPUchAllPart))/ele->PT;
            }
        }
    }
    std::vector<muon*> muon_tight_SIP_iso;
    for(muon* mu : tight_plus_SIP_muons) {
        if(mu->PFX_dB_fsr < .35) muon_tight_SIP_iso.push_back(mu);
    }
    //If we are using muons for the Z candidate, check that we still have a possible Z candidate after isolation is applied
    switch(type) {
        case b_3mu:
        case b_1e2mu:
            number_tight=muon_tight_SIP_iso.size();
            charge[0]=muon_tight_SIP_iso[0]->CHARGE;
            charge[1]=muon_tight_SIP_iso[1]->CHARGE;
            if(number_tight>2) charge[2]=muon_tight_SIP_iso[2]->CHARGE;
            if(number_tight<2) return;
            if(charge[0]==charge[1] && (number_tight==2 || charge[0]==charge[2])) return;
            break;
        case b_2e1mu:
        case b_3e:
            ;
    }
    std::array<four_vector,3> leps;
    std::array<four_vector,3> leps_fsr;
    std::array<void*,3> cands;
    std::array<int,3> charges;
    int OS_num;
    switch(number_tight) {
        default: //?!?!?!?
            std::cout << "?!?!?!?!?!?!\n";
            exit(1);
        case 2: //No ambiguities, just need to actually fill things in
            //Fill first two first (same for two different cases)
            switch(type) {
                case b_3mu:
                case b_1e2mu:
                    //Only two tight muons, so just use both
                    leps[0]=muon_tight_SIP_iso[0]->fv();
                    leps[1]=muon_tight_SIP_iso[1]->fv();
                    charges[0]=muon_tight_SIP_iso[0]->CHARGE;
                    charges[1]=muon_tight_SIP_iso[1]->CHARGE;
                    cands[0]=(void*)(muon_tight_SIP_iso[0]);
                    cands[1]=(void*)(muon_tight_SIP_iso[1]);
                    break;
                case b_3e:
                case b_2e1mu:
                    //Only two tight muons, so just use both
                    leps[0]=tight_plus_SIP_electrons[0]->fv();
                    leps[1]=tight_plus_SIP_electrons[1]->fv();
                    charges[0]=tight_plus_SIP_electrons[0]->CHARGE;
                    charges[1]=tight_plus_SIP_electrons[1]->CHARGE;
                    cands[0]=(void*)(tight_plus_SIP_electrons[0]);
                    cands[1]=(void*)(tight_plus_SIP_electrons[1]);
                    //Same as last but for electrons
            }
            //Now fill in the last one
            switch(type) {
                case b_3mu: //The loose one that is not a Z candidate is leps[2]
                    for(int i=0;i<3;++i) {
                        void* cand=(void*)loose_plus_SIP_muons[i];
                        if(cand==cands[0] || cand==cands[1]) continue;
                        leps[2]=((muon*)cand)->fv();
                        charges[2]=((muon*)cand)->CHARGE;
                        cands[2]=cand;
                        break;
                    }
                    break;
                case b_3e:
                    for(int i=0;i<3;++i) {
                        void* cand=(void*)loose_plus_SIP_electrons[i];
                        if(cand==cands[0] || cand==cands[1]) continue;
                        leps[2]=((electron*)cand)->fv();
                        charges[2]=((electron*)cand)->CHARGE;
                        cands[2]=cand;
                        break;
                    }
                    break;
                case b_1e2mu: //Only one loose electron
                    leps[2]=loose_plus_SIP_electrons[0]->fv();
                    charges[2]=loose_plus_SIP_electrons[0]->CHARGE;
                    cands[2]=(void*)(loose_plus_SIP_electrons[0]);
                    break;
                case b_2e1mu:
                    leps[2]=loose_plus_SIP_muons[0]->fv();
                    charges[2]=loose_plus_SIP_muons[0]->CHARGE;
                    cands[2]=(void*)(loose_plus_SIP_muons[0]);
                    break;
            } //End populating four vectors for number_tight==2
            //Now add fsr to lepton four vectors
            for(int i=0;i<3;++i) {
                auto fsr_photon=fsr_photons.find(cands[i]);
                leps_fsr[i]=leps[i];
                if(fsr_photon!=fsr_photons.end()) {
                    leps_fsr[i]+=fsr_photon->second.fv();
                }
            }
            break;
        case 3: //We need to find the pair closest to the nominal Z mass
            switch(type) {
                default:
                    std::cout << "?!?!?!?!?!\n";
                    exit(1);
                case b_3mu:
                    //Muon that is OS of other two must be included in Z candidate, so assign it to lep1. Assign the other two to lep2 and lep3 (we will swap later if the order is wrong)
                    if(muon_tight_SIP_iso[0]->CHARGE==muon_tight_SIP_iso[1]->CHARGE) {
                        OS_num=2;
                    } else if(muon_tight_SIP_iso[0]->CHARGE==muon_tight_SIP_iso[2]->CHARGE) {
                        OS_num=1;
                    } else {
                        OS_num=0;
                    }
                    for(int i=0;i<3;++i) {
                        leps[i]=muon_tight_SIP_iso[(OS_num+i)%3]->fv();
                        cands[i]=(void*)(muon_tight_SIP_iso[(OS_num+i)%3]);
                        charges[i]=(muon_tight_SIP_iso[(OS_num+i)%3])->CHARGE;
                        auto fsr_photon=fsr_photons.find(cands[i]);
                        //Add fsr in for Z reconstruction
                        leps_fsr[i]=leps[i];
                        if(fsr_photon!=fsr_photons.end()) {
                            leps_fsr[i]+=fsr_photon->second.fv();
                        }
                    }
                    break;
                case b_3e:
                    //Electron that is OS of other two must be included in Z candidate, so assign it to lep1. Assign the other two to lep2 and lep3 (we will swap later if the order is wrong)
                    if(tight_plus_SIP_electrons[0]->CHARGE==tight_plus_SIP_electrons[1]->CHARGE) {
                        OS_num=2;
                    } else if(tight_plus_SIP_electrons[0]->CHARGE==tight_plus_SIP_electrons[2]->CHARGE) {
                        OS_num=1;
                    } else {
                        OS_num=0;
                    }
                    for(int i=0;i<3;++i) {
                        leps[i]=tight_plus_SIP_electrons[(OS_num+i)%3]->fv();
                        cands[i]=(void*)(tight_plus_SIP_electrons[(OS_num+i)%3]);
                        charges[i]=(tight_plus_SIP_electrons[(OS_num+i)%3])->CHARGE;
                        auto fsr_photon=fsr_photons.find(cands[i]);
                        //Add fsr in for Z reconstruction
                        leps_fsr[i]=leps[i];
                        if(fsr_photon!=fsr_photons.end()) {
                            leps_fsr[i]+=fsr_photon->second.fv();
                        }
                    }
            } //End switch on type for case 3
            //Compute the two alternative z masses
            float Zmass=(leps_fsr[0]+leps_fsr[1]).M();
            float altZmass=(leps_fsr[0]+leps_fsr[2]).M();
            //If the alternative Z mass is closer to nominal, switch the candidates. After this leps[0] and leps[1] contain the OSSF pair that is closest to the nominal Z mass, including FSR
            if(fabs(Zmass-Z_nominal_mass) > fabs(altZmass-Z_nominal_mass)) {
                std::swap(leps[1],leps[2]);
                std::swap(leps_fsr[1],leps_fsr[2]);
                std::swap(charges[1],charges[2]); //This shouldn't do anything
                std::swap(cands[1],cands[2]);
            }
        //End case 3
    } //End switch(number_tight)
    //We are now left with the Lorentz vectors: 0 and 1 correspond to the Z candidate and 2 to the remaining lepton
    //We also have the charge and a void* to the original object

    //Now fill our tree
   
    //Store Z mass
    tree_Zmass_nofsr=(leps[0]+leps[1]).M();
    tree_Zmass_fsr=(leps_fsr[0]+leps_fsr[1]).M();
    //Store bkg type
    tree_bkg_type=type;
    tree_lep1_4v=&leps[0];
    tree_lep2_4v=&leps[1];
    tree_lep3_4v=&leps[2];
    tree_lep1_fsr_4v=&leps_fsr[0];
    tree_lep2_fsr_4v=&leps_fsr[1];
    tree_lep3_fsr_4v=&leps_fsr[2];
    tree_lep1_q=charges[0];
    tree_lep2_q=charges[0];
    tree_lep3_q=charges[0];

    tree_out->Fill();
}
