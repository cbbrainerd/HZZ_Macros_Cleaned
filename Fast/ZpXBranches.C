#include "ZpXanalyzer.h"
#include <iostream>

//branch->SetStatus(1);
//fChain->SetBranchAddress(branch_name,&variable,&branch);

//Gets branch from name, sets status to 1, and 
template <class T>
void activate_branch_impl(TTree* fChain,T& variable,TBranch*& tb,const char* vn,const char* tbn=nullptr) {
    auto myBranch=fChain->GetBranch(vn);
    if(!myBranch) {
        std::cout << "Could not activate branch " << tbn << " with name " << vn << "\n";
        std::cout << "Tree info: ";
        fChain->Print();
        std::cout << "===Tree info ends===\n";
        exit(1);
    }
    myBranch->SetStatus(1);
    fChain->SetBranchAddress(vn,&variable,&tb);
//    std::cout << "Activating branch for var " << vn << ", branch " << tbn << ", of name " << vn << ".\n";
}


//Otherwise
#define activate_branch_name(vn,bn) do{activate_branch_impl(fChain,(vn),(bn),#vn,#bn);}while(0)
//Assumes if variable name is `x` then branch is called `b_x`. This is the case for every branch except run, lumi, event
#define activate_branch(vn) activate_branch_name(vn,b_ ## vn)
//One branch has an extra space in its name. This would screw things up but I'm not using that branch

void ZpXanalyzer::enable_branches() {
    activate_branch_name(Run, b_irun);
    activate_branch_name(Event, b_ievt);
    activate_branch_name(LumiSection, b_ils);
    activate_branch(num_PU_vertices);
    activate_branch(MC_weighting);
    activate_branch(RECOELE_CHARGE);
    activate_branch(RECOELE_E);
    activate_branch(RECOELE_ETA);
    activate_branch(RECOELE_ID);
    activate_branch(RECOELE_IP);
    activate_branch(RECOELE_IPERROR);
    activate_branch(RECOELE_P);
    activate_branch(RECOELE_PFPUchAllPart);
    activate_branch(RECOELE_PFX_dB);
    activate_branch(RECOELE_PFchHad);
    activate_branch(RECOELE_PFneuHad);
    activate_branch(RECOELE_PFphoton);
    activate_branch(RECOELE_PHI);
    activate_branch(RECOELE_PT);
    activate_branch(RECOELE_PTError);
    activate_branch(RECOELE_SIP);
    activate_branch(RECOELE_gsftrack_dxy);
    activate_branch(RECOELE_gsftrack_dz);
    activate_branch(RECOELE_isGap);
    activate_branch(RECOELE_mvaNonTrigV0);
    activate_branch(RECOELE_scl_Eta);
    activate_branch(RECOELE_scl_Phi);
    activate_branch(RECOMU_CHARGE);
    activate_branch(RECOMU_ETA);
    activate_branch(RECOMU_IP);
    activate_branch(RECOMU_IPERROR);
    activate_branch(RECOMU_MASS);
    activate_branch(RECOMU_PFPUchAllPart);
    activate_branch(RECOMU_PFX_dB);
    activate_branch(RECOMU_PFchHad);
    activate_branch(RECOMU_PFneuHad);
    activate_branch(RECOMU_PFphoton);
    activate_branch(RECOMU_PHI);
    activate_branch(RECOMU_PT);
    activate_branch(RECOMU_SIP);
    activate_branch(RECOMU_THETA);
    activate_branch(RECOMU_TRACKISO);
    activate_branch(RECOMU_isGlobalMu);
    activate_branch(RECOMU_isPFMu);
    activate_branch(RECOMU_isTrackerHighPtMu);
    activate_branch(RECOMU_isTrackerMu);
    activate_branch(RECOMU_muInnertrkNPixHits);
    activate_branch(RECOMU_muInnertrkPTError);
    activate_branch(RECOMU_muInnertrktrackerLayersWithMeasurement);
    activate_branch(RECOMU_mubesttrkDxy);
    activate_branch(RECOMU_mubesttrkDz);
    activate_branch(RECOMU_mubesttrkPTError);
    activate_branch(RECOMU_mubesttrkType);
    activate_branch(RECOMU_mutrkNMuonHits);
    activate_branch(RECOMU_numberOfMatches);
    activate_branch(RECOMU_trkmuArbitration);
    activate_branch(RECOPFPHOT_ETA);
    activate_branch(RECOPFPHOT_PFPUchAllPart);
    activate_branch(RECOPFPHOT_PFX_rho);
    activate_branch(RECOPFPHOT_PFchHad);
    activate_branch(RECOPFPHOT_PFneuHad);
    activate_branch(RECOPFPHOT_PFphoton);
    activate_branch(RECOPFPHOT_PHI);
    activate_branch(RECOPFPHOT_PT);
    activate_branch(RECOPFPHOT_PTError);
    activate_branch(RECOPFPHOT_PT_uncorr);
    activate_branch(RECO_PFMET);
    activate_branch(RECO_PFMET_PHI);
    activate_branch(RECO_PFMET_PHI_xycorr);
    activate_branch(RECO_PFMET_xycorr);
    activate_branch(RECO_PUPPIMET);
    activate_branch(RECO_PUPPIMET_PHI);
}
