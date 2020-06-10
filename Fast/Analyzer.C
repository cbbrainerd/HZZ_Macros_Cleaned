#include <limits>
#include "Analyzer.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include <unordered_map>
#include <algorithm>
#include <iostream>

void break_here(volatile int x) { volatile int y; 
    std::cout << y << std::endl; 
if(y) exit(5); }
#define ERROR do {std::cout << "Error: line " << __LINE__ << "\n"; return cut_ERROR;} while(0)

//To do: forgot ghost removal deltaR > .02 between any leptons (it's fine we can do that in the ntuples)

typedef ROOT::Math::PtEtaPhiMVector four_vector;
using ROOT::Math::VectorUtil::DeltaR2;

float Z_nominal_mass=91.188;

void Analyzer::Loop() {
    if (fChain == 0) return;
    Long64_t jentry;
    double cuts[cut_NCUTS]={};
    for (jentry=0; jentry < fChain->GetEntriesFast(); jentry++) {
       if(jentry % 1000==0) {
         std::cout << "Processing entry " << jentry << std::endl;
         verbose=true;
       }
       Long64_t ientry=LoadTree(jentry);
       fChain->GetEntry(jentry);
       tree_cutflow=analyze();
       if(tree_cutflow==-1) {
           std::cout << "Error!\n";
       } else {
           ++cuts[tree_cutflow];
       }
       if(tree_cutflow==0 || tree_cutflow > cut_NOZCAND) {
           tree_out->Fill();
       } 
       // if (Cut(ientry) < 0) continue;
    }
    //cutflow histogram
    const char* cutflow_labels[cut_NCUTS]={"Less than 2 leptons","Less than 2 tight leptons","No Z candidate","No 2 OSSF pairs","No ZZ selected","PASS"};
    //If a label is missing (i.e. nullptr), fill in a placeholder
    for(int i=0;i<cut_NCUTS;++i) {
        if(!cutflow_labels[i]) {
            std::cout << "Warning: missing cutflow label.\n";
            cutflow_labels[i]="Unlabeled.";
        }
    }
    cutflow_hist=new TH1D("cutflow","cutflow",cut_NCUTS,0,1);
    for(int i=0;i<cut_NCUTS;++i) {
        //jentry is initially the full number of events. Subtract off the number lost at each point in the cutflow
        jentry-=cuts[i];
        //Note that the first non-underflow bin of a TH1 is 1, not 0
        cutflow_hist->SetBinContent(i+1,jentry);
        cutflow_hist->GetXaxis()->SetBinLabel(i+1,cutflow_labels[i]);
    }
}

bool Analyzer::Notify() {
    enable_branches();
}

//Lepton ID
//isMuon true: muon, false: lepton
//isTight true: tight ID (sip, iso, etc) false: loose + SIP
//Returns value of failed cut, or -1 for error condition
int Analyzer::analyze() {
    //lepton IDs. Also fill std::vector of leptons so we don't have to access through RECO... branches
    if(is_MC) {
        tree_pileup_weight=pileup_corr.get_pileup_weight(num_PU_vertices);
    } else {
        tree_pileup_weight=1;
    }
    //Set weight to zero for later steps as easy way to ignore these steps in plotter
    tree_Z_step_weight=0;
    tree_ZZ_step_weight=0;
    //This vector is of loose+SIP muons, the loosest that is useful in the full analysis
    std::vector<muon> muons;
    for(std::size_t i=0;i!=RECOMU_PT->size();++i) {
        muon tmp_muon(this,i);
        if(tmp_muon.loose_plus_SIP) {
            muons.push_back(tmp_muon);
        }
    }
    
    //Populate the electrons and do cross cleaning
    std::vector<electron> electrons;
    {
        //Vector of pointers quality vectors to avoid copying
        std::vector<muon*> tight_muons;
        for(auto& mu : muons) {
            if(mu.tight_plus_SIP) tight_muons.push_back(&mu);
        }
        
        for(std::size_t i=0;i!=RECOELE_PT->size();++i) {
            electron ele(this,i);
            //Only loose plus SIP electrons or tighter are useful in this analysis
            if(!ele.loose_plus_SIP) continue;
            //Cross cleaning
            //Electrons are removed from consideration if they are too close to a tight+SIP muon
            bool cross_clean=true;
            for(muon* mu : tight_muons) {
                double deltaRsq=(pow(DELTAPHI(mu->PHI,ele.PHI),2)+pow(mu->ETA-ele.ETA,2));
                if(deltaRsq <= 0.0025) {
                    cross_clean=false;
                    break;
                }
            }
            if(cross_clean) electrons.push_back(ele);
        }
    }

    if(std::max(electrons.size(),muons.size()) < 2) return cut_NLEPTONS;
    //To store fsr photons for each lepton
    std::unordered_map<void*,photon> fsr_photons;
    for(std::size_t i=0;i!=RECOPFPHOT_PT->size();++i) {
        photon phot(this,i);
        bool clean=1;
        for(electron& ele : electrons) {
            double deltaPhi=fabs(DELTAPHI(phot.PHI,ele.scl_Phi));
            double deltaEta=fabs(phot.ETA-ele.scl_Eta);
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
        for(muon& mu : muons) {
            double deltaR=sqrt(pow(DELTAPHI(phot.PHI,mu.PHI),2)+pow(phot.ETA-mu.ETA,2));
            double deltaRoverPtSq=deltaR/pow(phot.PT,2);
            if(deltaRoverPtSq<0.012 && deltaR<minDeltaR) {
                minDeltaR=deltaR;
                //Save the deltaRoverPtSq to the photon for later
                phot.deltaR=deltaR;
                phot.deltaRoverPtSq=deltaRoverPtSq;
                nearest_lepton=(void*)&mu;
            }
        }
        //Should this use supercluster eta, etc?
        for(electron& ele : electrons) {
            double deltaR=sqrt(pow(DELTAPHI(phot.PHI,ele.PHI),2)+pow(phot.ETA-ele.ETA,2));
            double deltaRoverPtSq=deltaR/pow(phot.PT,2);
            if(deltaRoverPtSq<.012 && deltaR<minDeltaR) {
                minDeltaR=deltaR;
                phot.deltaR=deltaR;
                phot.deltaRoverPtSq=deltaRoverPtSq;
                nearest_lepton=(void*)&ele;
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
    //Compute isolation for all muons of interest, subtracting *all* fsr photons within isolation cone
    //Iso=(PFchHad+max(0.,PFneuHad+PFphoton-.5PFPUchAllPart))/PT
    //Where PFphoton subtracts off the fsr photon if it exists within the cone .01 < deltaR <= .3
    for(muon& mu : muons) {
        //Compute isolation before FSR subtraction
        mu.PFX_dB_prefsr=(mu.PFchHad+std::max(0.,mu.PFneuHad+mu.PFphoton-.5*mu.PFPUchAllPart))/mu.PT;
        mu.PFphoton_fsr=mu.PFphoton;
        for(auto fsr_photon : fsr_photons) {
            photon const& phot=fsr_photon.second;
            //If this is the matched for this muon, save a pointer to it for later use
            if(&mu==fsr_photon.first) mu.matched_photon=&phot;
            //Subtract off the photon if it is within the isolation cone
            double deltaR2=pow(DELTAPHI(phot.PHI,mu.PHI),2)+pow(phot.ETA-mu.ETA,2);
            if (deltaR2 < .09 && deltaR2 > .0001) {
                mu.PFphoton_fsr-=phot.PT;
            }
        }
        mu.PFX_dB_fsr=(mu.PFchHad+std::max(0.,mu.PFneuHad+(mu.PFphoton_fsr)-.5*mu.PFPUchAllPart))/mu.PT;
    }
    for(electron& ele : electrons) {
        ele.PFX_dB_prefsr=(ele.PFchHad+std::max(0.,ele.PFneuHad+ele.PFphoton-.5*ele.PFPUchAllPart))/ele.PT;
        ele.PFphoton_fsr=ele.PFphoton;
        for(auto fsr_photon : fsr_photons) {
            photon const& phot=fsr_photon.second;
            //If this is the matched for this electron, save a pointer to it for later use
            if(&ele==fsr_photon.first) ele.matched_photon=&phot;
            //Electron isolation < .3 delta R && (supercluster eta < 1.479 OR deltaR > .08)
            double deltaR2=pow(DELTAPHI(phot.PHI,ele.PHI),2)+pow(phot.ETA-ele.ETA,2);
            if (deltaR2 < .09 && (ele.scl_Eta < 1.479 || deltaR2 > .0064)) {
                ele.PFphoton_fsr-=phot.PT;
            }
        }
        ele.PFX_dB_fsr=(ele.PFchHad+std::max(0.,ele.PFneuHad+(ele.PFphoton_fsr)-.5*ele.PFPUchAllPart))/ele.PT;
    }
    //From here on, only fully selected leptons are of any use, so discard the rest
    muons.erase(std::remove_if(muons.begin(),muons.end(),[](muon const& mu) { return !(mu.tight_plus_SIP && mu.PFX_dB_fsr < .35); }),muons.end());
    //Note that 
    electrons.erase(std::remove_if(electrons.begin(),electrons.end(),[](electron const& ele) { return !ele.tight; }),electrons.end());

    if(std::max(muons.size(),electrons.size()) < 2) return cut_NLEPTONS_TIGHT;

    std::vector<four_vector> positive_muons;
    std::vector<four_vector> negative_muons;
    std::vector<four_vector> positive_muons_fsr;
    std::vector<four_vector> negative_muons_fsr;
    for(muon& mu : muons) {
        if(mu.CHARGE > 0) {
            positive_muons.push_back(mu.fv());
            positive_muons_fsr.push_back(mu.fv_fsr());
        } else {
            negative_muons.push_back(mu.fv());
            negative_muons_fsr.push_back(mu.fv_fsr());
        }
    }
    std::vector<four_vector> positive_electrons;
    std::vector<four_vector> negative_electrons;
    std::vector<four_vector> positive_electrons_fsr;
    std::vector<four_vector> negative_electrons_fsr;
    for(electron& ele : electrons) {
        if(ele.CHARGE > 0) {
            positive_electrons.push_back(ele.fv());
            positive_electrons_fsr.push_back(ele.fv_fsr());
        } else {
            negative_electrons.push_back(ele.fv());
            negative_electrons_fsr.push_back(ele.fv_fsr());
        }
    }
    assert(positive_muons.size()==positive_muons_fsr.size());
    assert(negative_muons.size()==negative_muons_fsr.size());
    assert(positive_electrons.size()==positive_electrons_fsr.size());
    assert(negative_electrons.size()==negative_electrons_fsr.size());
    std::size_t zcands_n=(positive_muons.size()*negative_muons.size())+(positive_electrons.size()*negative_electrons.size());
    if(zcands_n==0) return cut_NOZCAND;
    std::size_t zcands_i;
    enum particle_type {
        pt_electron=-1,
        pt_unknown=0,
        pt_muon=1
    };
    particle_type best_cand=pt_unknown; //-1 electron, 0 error condition (somehow no cand?), 1 muon
    std::size_t best_cand_p_i=-1;
    std::size_t best_cand_n_i=-1; 
    double bestZmass=-1, bestZmass_nofsr=-1;
    double best_mass_diff=std::numeric_limits<double>::max();
    for(std::size_t p_mu_i=0;p_mu_i!=positive_muons.size();++p_mu_i) {
        for(std::size_t n_mu_i=0;n_mu_i!=negative_muons.size();++n_mu_i) {
            double Zcand_mass=(positive_muons_fsr.at(p_mu_i)+negative_muons_fsr.at(n_mu_i)).M();
            double Zcand_mass_nofsr=(positive_muons.at(p_mu_i)+negative_muons.at(n_mu_i)).M();
            double Zmass_diff=fabs(Zcand_mass-Z_nominal_mass);
            if(Zmass_diff < best_mass_diff) { 
                best_cand_p_i=p_mu_i;
                best_cand_n_i=n_mu_i;
                bestZmass=Zcand_mass;
                bestZmass_nofsr=Zcand_mass_nofsr;
                best_mass_diff=Zmass_diff;
                best_cand=pt_muon;
            }
        }
    }
    for(std::size_t p_ele_i=0;p_ele_i!=positive_electrons.size();++p_ele_i) {
        for(std::size_t n_ele_i=0;n_ele_i!=negative_electrons.size();++n_ele_i) {
            double Zcand_mass=(positive_electrons_fsr.at(p_ele_i)+negative_electrons_fsr.at(n_ele_i)).M();
            if(Zcand_mass < 0) {
                std::cout << Zcand_mass << std::endl;
                break_here(5);
            }
            double Zcand_mass_nofsr=(positive_electrons.at(p_ele_i)+negative_electrons.at(n_ele_i)).M();
            double Zmass_diff=fabs(Zcand_mass-Z_nominal_mass);
            if(Zmass_diff < best_mass_diff) { 
                best_cand_p_i=p_ele_i;
                best_cand_n_i=n_ele_i;
                bestZmass=Zcand_mass;
                bestZmass_nofsr=Zcand_mass_nofsr;
                best_mass_diff=Zmass_diff;
                best_cand=pt_electron;
            }
        }
    }
    tree_leading_Z_mass_fsr=bestZmass;
    if(tree_leading_Z_mass_fsr < 0) {
        //if(best_cand==pt_electron) std::cout << (positive_electrons[0]+negative_electrons[0]).M() << std::endl;
        break_here(5);
    }
    tree_leading_Z_mass_nofsr=bestZmass_nofsr;
    tree_leading_Z_type=best_cand;
    //Shouldn't be possible
    if(best_cand==pt_unknown) ERROR;
    //Compute scale factors for Z
    tree_Z_step_weight=tree_pileup_weight;
    if(is_MC) {
        if(best_cand==pt_electron) {
            std::size_t e_p=0,e_n=0;
            bool match=false;
            for(electron &ele : electrons) {
                if(ele.CHARGE > 0) {
                    if(best_cand_p_i==e_p) {
                        match=true;
                    }
                    ++e_p;
                } else {
                    if(best_cand_n_i==e_n) {
                        match=true;
                    }
                    ++e_n;
                }
                if(match) {
                    tree_Z_step_weight*=scale_factors_ele.get_scale_factor(ele.SCL_ETA,ele.PT,ele.isGap);
                    tree_Z_step_weight*=scale_factors_ele.get_efficiency(ele.SCL_ETA,ele.PT);
                    match=false;
                }
            }
        } else if (best_cand==pt_muon) {
            std::size_t m_p=0,m_n=0;
            bool match=false;
            for(muon &mu : muons) {
                if(mu.CHARGE > 0) {
                    if(best_cand_p_i==m_p) {
                        match=true;
                    }
                    ++m_p;
                } else {
                    if(best_cand_n_i==m_n) {
                        match=true;
                    }
                    ++m_n;
                }
                if(match) {
                    tree_Z_step_weight*=scale_factors_mu.get_scale_factor(mu.ETA,mu.PT);
                    match=false;
                }
            }
        }
    }
    //Now do ZZ reconstruction. Pick the ZZ candidate with the highest sum pt. (For 4e, 4mu extra step after)
    if(verbose) std::cout << "Beginning ZZ reconstruction\n";
    verbose=false;
    double sum_pt=-1;
    channel_type best_ZZ=c_unknown;
    particle_type z1_type=pt_unknown;
    particle_type z2_type=pt_unknown;
    //For storing indices of matched particles
    std::size_t z1_pi=-1,z1_ni=-1,z2_pi=-1,z2_ni=-1;
    //Loop over all 2+/2- combinations
    std::size_t num_pm=positive_muons.size();
    std::size_t num_nm=negative_muons.size();
    std::size_t num_pe=positive_electrons.size();
    std::size_t num_ne=negative_electrons.size();
    std::size_t num_p=num_pm+num_pe;
    std::size_t num_n=num_nm+num_ne;
    if(verbose) std::cout << "Positive leptons: " << num_p << "Negative leptons: " << num_n << std::endl;
    if(std::min(num_pm,num_nm)+std::min(num_pe,num_ne) < 2) return cut_NO_OSSF_PAIRS;
    for(std::size_t np1=0;np1!=num_p;++np1) {
        for(std::size_t nn1=0;nn1!=num_n;++nn1) {
            for(std::size_t np2=np1+1;np2<num_p;++np2) {
                for(std::size_t nn2=nn1+1;nn2<num_n;++nn2) {
                    particle_type ptp1=(np1 >= positive_muons.size() ? pt_electron : pt_muon);
                    particle_type ptp2=(np2 >= positive_muons.size() ? pt_electron : pt_muon);
                    particle_type ptn1=(nn1 >= negative_muons.size() ? pt_electron : pt_muon);
                    particle_type ptn2=(nn2 >= negative_muons.size() ? pt_electron : pt_muon);
                    //Require particle types to match (all same flavor, or one of each)
                    if(ptp1!=ptn1 || ptp2!=ptn2) continue;
                    channel_type this_channel;
                    if(ptp1==pt_muon && ptp2==pt_muon) this_channel=c_4mu;
                    else if(ptp1==pt_electron && ptp2==pt_electron) this_channel=c_4e;
                    else this_channel=c_2e2mu;
                    four_vector const &p1=(ptp1 == pt_muon ? positive_muons.at(np1) : positive_electrons.at(np1-num_pm));
                    four_vector const &p2=(ptp2 == pt_muon ? positive_muons.at(np2) : positive_electrons.at(np2-num_pm));
                    four_vector const &n1=(ptn1 == pt_muon ? negative_muons.at(nn1) : negative_electrons.at(nn1-num_nm));
                    four_vector const &n2=(ptn2 == pt_muon ? negative_muons.at(nn2) : negative_electrons.at(nn2-num_nm));
                    four_vector const &p1_fsr=(ptp1 == pt_muon ? positive_muons_fsr.at(np1) : positive_electrons_fsr.at(np1-num_pm));
                    four_vector const &p2_fsr=(ptp2 == pt_muon ? positive_muons_fsr.at(np2) : positive_electrons_fsr.at(np2-num_pm));
                    four_vector const &n1_fsr=(ptn1 == pt_muon ? negative_muons_fsr.at(nn1) : negative_electrons_fsr.at(nn1-num_nm));
                    four_vector const &n2_fsr=(ptn2 == pt_muon ? negative_muons_fsr.at(nn2) : negative_electrons_fsr.at(nn2-num_nm));
                    double my_sum_pt=(p1_fsr+p2_fsr+n1_fsr+n2_fsr).Pt();
                    if(my_sum_pt < sum_pt) continue; //Pick ZZ candidate with greatest sum pt. This will be zero if a ZZ candidate has been selected previously.
                    //Pt cuts: 20/10. Don't think fsr recovery makes sense here since those wouldn't trigger
                    if(!(p1.pt() > 20 || p2.pt() > 20 || n1.pt() > 20 || n2.pt() > 20)) continue;
                    if(((p1.pt() > 10)+(p2.pt() > 10)+(n1.pt() > 10)+(n2.pt() > 10)) < 2) continue;
                    //4l mass > 70
                    if(!((p1+p2+n1+n2).M() > 70)) continue;
                    //Ghost removal, deltaR > .02 between any pair of leptons, regardless of sign/flavor
                    if(!(DeltaR2(p1,p2) > .0004 && DeltaR2(n1,n2) > .0004 && DeltaR2(p1,n1) > .0004 && DeltaR2(p1,n2) > .0004 && DeltaR2(p2,n1) > .0004 && DeltaR2(p2,n2) > .0004)) continue;
                    //QCD suppression, all OS pairs > 4 mass
                    if(!((p1+n1).M() > 4 && (p1+n2).M() > 4 && (p2+n1).M() > 4 && (p2+n2).M() > 4)) continue;
                    //Select Z1. Both Z1 and Z2 must be between 12 and 120 (i.e. fabs(m-66) < 54)
                    
                    double candidate_mass=-1;
                    double new_candidate_mass=(p1_fsr+n1_fsr).M();
                    double new_alternate_mass=(p2_fsr+n2_fsr).M();
                    double mass_difference=std::numeric_limits<double>::max();
                    //Identify leptons to which Z they are connected to
                    four_vector const* z1p; four_vector const* z1n; four_vector const* z2p; four_vector const* z2n; four_vector const* z1p_fsr; four_vector const* z1n_fsr; four_vector const* z2p_fsr; four_vector const* z2n_fsr;
                    int candidate_ordering=-1; //Remember which ordering we pick
                    if(fabs(new_candidate_mass-66) > 54 && fabs(new_alternate_mass-66) > 54) {
                        mass_difference=fabs(new_candidate_mass-Z_nominal_mass);
                        double alt_mass_difference=fabs(new_alternate_mass-Z_nominal_mass);
                        if(mass_difference < alt_mass_difference) {
                            z1p=&p1; z1n=&n1; z2p=&p2; z2n=&n2;
                            z1p_fsr=&p1_fsr; z1n_fsr=&n1_fsr; z2p_fsr=&p2_fsr; z2n_fsr=&n2_fsr;
                            candidate_ordering=0;
                        } else {
                            z1p=&p2; z2p=&p1; z1n=&n1; z2n=&n2;
                            z1p_fsr=&p2_fsr; z2p_fsr=&p1_fsr; z1n_fsr=&n1_fsr; z2n_fsr=&n2_fsr;
                            candidate_ordering=1;
                            mass_difference=alt_mass_difference;
                        }
                    }
                    //Try the other two pairs if 4e or 4mu
                    if(this_channel != c_2e2mu) {
                        //p1+n2 gives z1
                        new_candidate_mass=fabs((p1_fsr+n2_fsr).M()-Z_nominal_mass);
                        new_alternate_mass=fabs((p2_fsr+n1_fsr).M()-Z_nominal_mass);
                        //Both Z1 and Z2 must be within 12 to 120 to be considered
                        if(fabs(new_candidate_mass-66) > 54 && fabs(new_alternate_mass-66) > 54) {
                            double new_mass_difference=fabs(new_candidate_mass-Z_nominal_mass);
                            if(new_mass_difference < mass_difference) {
                                z1p=&p1; z2p=&p2; z1n=&n2; z2n=&n1;
                                z1p_fsr=&p1_fsr; z2p_fsr=&p2_fsr; z1n_fsr=&n2_fsr; z2n_fsr=&n1_fsr;
                                candidate_mass=new_candidate_mass;
                                mass_difference=new_mass_difference;
                                candidate_ordering=2;
                            }
                            new_mass_difference=fabs(new_alternate_mass-Z_nominal_mass);
                            if(new_candidate_mass < candidate_mass) {
                                z1p=&p2; z2p=&p1; z1n=&n1; z2n=&n2;
                                z1p_fsr=&p2_fsr; z2p_fsr=&p1_fsr; z1n_fsr=&n1_fsr; z2n_fsr=&n2_fsr;
                                candidate_mass=new_candidate_mass;
                                mass_difference=new_mass_difference;
                                candidate_ordering=3;
                            }   
                        }
                    }
                    if(candidate_ordering==-1) continue; //Failed to find any matching candidate
                    four_vector z1=((*z1p)+(*z1n));
                    four_vector z1_fsr=((*z1p_fsr)+(*z1n_fsr));
                    //Require Z1 mass greater than 40
                    if(!(z1_fsr.M() > 40)) continue;
                    four_vector z2=((*z2p)+(*z2n));
                    four_vector z2_fsr=((*z2p_fsr)+(*z2n_fsr));
                    double mZ1_fsr=z1_fsr.M();
                    //Smart cut
                    if(this_channel != c_2e2mu) {
                        //Za is the Z under the alternate pairing closest to Z mass. Construct the pairing
                        double mZa=((*z1p_fsr)+(*z2n_fsr)).M();
                        double mZb=((*z2p_fsr)+(*z1n_fsr)).M();
                        //Swap Za and Zb if we picked incorrectly
                        if(fabs(mZb-Z_nominal_mass) < fabs(mZa-Z_nominal_mass)) std::swap(mZa,mZb);
                        if(fabs(mZa-mZ1_fsr) < fabs(mZ1_fsr-Z_nominal_mass) && mZb < 12) continue; //Failed smart cut (looks like on-shell Z plus low mass ll)
                    }
                    //If we reach here, we have accepted a ZZ candidate (we may accept more later, so be sure to overwrite everything if that's the case
                    if(verbose) std::cout << "ZZ candidate accepted\n";
                    tree_Z1_mass_fsr=z1_fsr.M();
                    if(tree_Z1_mass_fsr > 200) std::cout << tree_Z1_mass_fsr << std::endl;
                    tree_Z2_mass_fsr=z2_fsr.M();
                    tree_Z1_mass_nofsr=z1.M();
                    tree_Z2_mass_nofsr=z2.M();
                    if(sum_pt > my_sum_pt) ERROR;
                    sum_pt=my_sum_pt; //Record the sum pt as the highest so far
                    tree_channel_type=this_channel;
                    //Record indices of which leptons match
                    switch(candidate_ordering) {
                        case 0:
                            z1_pi=np1;
                            z1_ni=nn1;
                            z2_pi=np2;
                            z2_ni=nn2;
                            z1_type=ptp1;
                            z2_type=ptp2;
                            break;
                        case 1:
                            z1_pi=np2;
                            z1_ni=nn2;
                            z2_pi=np1;
                            z2_ni=nn1;
                            z1_type=ptp2;
                            z2_type=ptp1;
                            break;
                        case 2:
                            z1_pi=np1;
                            z1_ni=nn2;
                            z2_pi=np2;
                            z2_ni=nn1;
                            z1_type=ptp1;
                            z2_type=ptp2;
                            break;
                        case 3:
                            z1_pi=np1;
                            z1_ni=nn2;
                            z2_pi=np2;
                            z2_ni=nn1;
                            z1_type=ptp1;
                            z2_type=ptp2;
                            break;
                        default:
                            ERROR;
                    }
                }
            }
        }
    } //End ZZ candidate selection
    if(sum_pt < 0) { //No ZZ candidate selected
        return cut_NO_ZZ_SELECTED;
    }
    tree_ZZ_step_weight=tree_pileup_weight;
    particle_type types[2]={z1_type,z2_type};
    std::size_t index_p[2]={z1_pi,z2_pi};
    std::size_t index_n[2]={z1_ni,z2_ni};
    for(int i=0;i!=2;++i) {
        particle_type best_cand=types[i];
        std::size_t best_cand_p_i=index_p[i];
        std::size_t best_cand_n_i=index_n[i];
        if(best_cand==pt_electron) {
            best_cand_p_i-=positive_muons.size();
            best_cand_n_i-=negative_muons.size();
        }
        if(is_MC) {
            if(best_cand==pt_electron) {
                std::size_t e_p=0,e_n=0;
                bool match=false;
                for(electron &ele : electrons) {
                    if(ele.CHARGE > 0) {
                        if(best_cand_p_i==e_p) {
                            match=true;
                        }
                        ++e_p;
                    } else {
                        if(best_cand_n_i==e_n) {
                            match=true;
                        }
                        ++e_n;
                    }
                    if(match) {
                        tree_ZZ_step_weight*=scale_factors_ele.get_scale_factor(ele.SCL_ETA,ele.PT,ele.isGap);
                        tree_ZZ_step_weight*=scale_factors_ele.get_efficiency(ele.SCL_ETA,ele.PT);
                        match=false;
                    }
                }
            } else if (best_cand==pt_muon) {
                std::size_t m_p=0,m_n=0;
                bool match=false;
                for(muon &mu : muons) {
                    if(mu.CHARGE > 0) {
                        if(best_cand_p_i==m_p) {
                            match=true;
                        }
                        ++m_p;
                    } else {
                        if(best_cand_n_i==m_n) {
                            match=true;
                        }
                        ++m_n;
                    }
                    if(match) {
                        tree_ZZ_step_weight*=scale_factors_mu.get_scale_factor(mu.ETA,mu.PT);
                        match=false;
                    }
                }
            }
        }
    }
    if(verbose) std::cout << tree_Z1_mass_fsr << std::endl;
    return cut_PASS;
}
