#include "NewNtuple.h"
#include "Math/Vector4D.h"
#include "TH1I.h"

//3mu,1e2mu,2e1mu,3e
enum bkg_type : unsigned char {
    b_3mu=0, 
    b_1e2mu=1, 
    b_2e1mu=2, 
    b_3e=3
};

typedef ROOT::Math::PtEtaPhiMVector four_vector;

class ZpXanalyzer : public NewNtuple {
    void analyze();
public:
    void Loop() override;
    TTree* tree_out;
    float tree_Zmass_nofsr;
    float tree_Zmass_fsr;
    bkg_type tree_bkg_type;
    four_vector *tree_lep1_4v;
    four_vector *tree_lep2_4v;
    four_vector *tree_lep3_4v;
    four_vector *tree_lep1_fsr_4v;
    four_vector *tree_lep2_fsr_4v;
    four_vector *tree_lep3_fsr_4v;
    int tree_lep1_q;
    int tree_lep2_q;
    int tree_lep3_q;
    int tree_lep3_tight;
    TFile* outfile;
    TH1* type_hist;
    void enable_branches();
    ZpXanalyzer(TTree* intree,TFile* outfile_) : NewNtuple(intree) , outfile(outfile_) {
        //Disable all branches for performance reasons
        intree->SetBranchStatus("*",0);
        Init(intree);
        //Enable used branches. Be careful to add any additional branches to ZpXBranches.C!
        enable_branches();
        outfile->cd();
        type_hist=new TH1I("Bkg Types","Bkg Types",4,-.5,3.5);
        const char *types[4]={"Z(mu+mu)+mu","Z(mu+mu)+e","Z(e+e)+mu","Z(e+e)+e"};
        for(int i=1;i<=4;++i) type_hist->GetXaxis()->SetBinLabel(i,types[i-1]);
        tree_out=new TTree("ZpXTree","ZpXTree");
        tree_out->Branch("Run",&Run);
        tree_out->Branch("Event",&Event);
        tree_out->Branch("LumiSection",&LumiSection);
        tree_out->Branch("Avginstlumi",&Avginstlumi);
        tree_out->Branch("num_PU_vertices",&num_PU_vertices);
        tree_out->Branch("PU_BunchCrossing",&PU_BunchCrossing);
        tree_out->Branch("MC_weighting",&MC_weighting);
        tree_out->Branch("Zmass_nofsr",&tree_Zmass_nofsr);
        tree_out->Branch("Zmass_fsr",&tree_Zmass_fsr);
        tree_out->Branch("bkg_type",(unsigned char*)&tree_bkg_type);
        tree_out->Branch("lep1_4v",&tree_lep1_4v);
        tree_out->Branch("lep2_4v",&tree_lep2_4v);
        tree_out->Branch("lep3_4v",&tree_lep3_4v);
        tree_out->Branch("lep1_fsr_4v",&tree_lep1_fsr_4v);
        tree_out->Branch("lep2_fsr_4v",&tree_lep2_fsr_4v);
        tree_out->Branch("lep3_fsr_4v",&tree_lep3_fsr_4v);
        tree_out->Branch("lep1_q",&tree_lep1_q);
        tree_out->Branch("lep2_q",&tree_lep2_q);
        tree_out->Branch("lep3_q",&tree_lep3_q);
        tree_out->Branch("passes_tight",&tree_lep3_tight);
        tree_out->Branch("pfmet",&RECO_PFMET);
        tree_out->Branch("pfmet_phi",&RECO_PFMET_PHI);
        tree_out->Branch("puppimet",&RECO_PUPPIMET);
        tree_out->Branch("puppimet_phi",&RECO_PUPPIMET_PHI);
        tree_out->Branch("pfmet_xycorr",&RECO_PFMET_xycorr);
        tree_out->Branch("pfmet_xycorr_phi",&RECO_PFMET_PHI_xycorr);
    }  
    ~ZpXanalyzer() {
        tree_out->Write();
        type_hist->Write();
    }
};

//Get all the information from recomu that is actually used in the analysis
struct muon {
    bool isPFMu;
    bool isGlobalMu;
    bool isTrackerMu;
    bool isTrackerHighPtMu;
    float PT;
    float ETA;
    float THETA;
    float PHI;
    float MASS;
    float CHARGE;
    float TRACKISO;
    double PFchHad;
    double PFneuHad;
    double PFphoton;
    double PFPUchAllPart;
    double PFX_dB;
    double PFX_dB_prefsr;
    double PFX_dB_fsr;
    float SIP;
    float IP;
    float IPERROR;
    int numberOfMatches;
    int trkmuArbitration;
    float mutrkNMuonHits;
    float muInnertrktrackerLayersWithMeasurement;
    float muInnertrkPTError;
    float muInnertrkNPixHits;
    int mubesttrkType;
    float mubesttrkDxy;
    float mubesttrkDz;
    float mubesttrkPTError;
    bool loose;
    bool tight;
    bool loose_plus_SIP;
    bool tight_plus_SIP;
    muon(ZpXanalyzer* x,int i) {
        //These branches are always filled
        isPFMu=x->RECOMU_isPFMu->at(i);
        isGlobalMu=x->RECOMU_isGlobalMu->at(i);
        isTrackerMu=x->RECOMU_isTrackerMu->at(i);
        isTrackerHighPtMu=x->RECOMU_isTrackerHighPtMu->at(i);
        PT=x->RECOMU_PT->at(i);
        ETA=x->RECOMU_ETA->at(i);
        THETA=x->RECOMU_THETA->at(i);
        PHI=x->RECOMU_PHI->at(i);
        MASS=x->RECOMU_MASS->at(i);
        CHARGE=x->RECOMU_CHARGE->at(i);
        TRACKISO=x->RECOMU_TRACKISO->at(i);
        PFchHad=x->RECOMU_PFchHad->at(i);
        PFneuHad=x->RECOMU_PFneuHad->at(i);
        PFphoton=x->RECOMU_PFphoton->at(i);
        PFPUchAllPart=x->RECOMU_PFPUchAllPart->at(i);
        PFX_dB=x->RECOMU_PFX_dB->at(i);
        SIP=x->RECOMU_SIP->at(i);
        IP=x->RECOMU_IP->at(i);
        IPERROR=x->RECOMU_IPERROR->at(i);
        numberOfMatches=x->RECOMU_numberOfMatches->at(i);
        trkmuArbitration=x->RECOMU_trkmuArbitration->at(i);
        muInnertrkPTError=x->RECOMU_muInnertrkPTError->at(i);
        muInnertrkNPixHits=x->RECOMU_muInnertrkNPixHits->at(i);
        //These branches are always filled, or not, together.
        if(i<x->RECOMU_mutrkNMuonHits->size()) {
            mutrkNMuonHits=x->RECOMU_mutrkNMuonHits->at(i);
            muInnertrktrackerLayersWithMeasurement=x->RECOMU_muInnertrktrackerLayersWithMeasurement->at(i);
        } else {
            mutrkNMuonHits=-999;
        }   
        if(i<x->RECOMU_mubesttrkType->size()) {
            mubesttrkType=x->RECOMU_mubesttrkType->at(i);
            mubesttrkDxy=x->RECOMU_mubesttrkDxy->at(i);
            mubesttrkDz=x->RECOMU_mubesttrkDz->at(i);
            mubesttrkPTError=x->RECOMU_mubesttrkPTError->at(i);
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
    four_vector fv() const {
        return four_vector(PT,ETA,PHI,0.106);
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
        return fabs(SIP) < 4.;
    }
};

struct electron{
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
    double PFX_dB;
    double PFX_dB_prefsr;
    double PFX_dB_fsr;
    float SIP;
    float IP;
    float IPERROR;
    float scl_Phi;
    float scl_Eta;
    float gsftrack_dxy,gsftrack_dz;
    double mvaNonTrigV0;
    bool loose, tight, loose_plus_SIP;
    electron(ZpXanalyzer *x,int i) {
        scl_Phi=x->RECOELE_scl_Phi->at(i);
        scl_Eta=x->RECOELE_scl_Eta->at(i);
        E=x->RECOELE_E->at(i);
        PT=x->RECOELE_PT->at(i);
        PTError=x->RECOELE_PTError->at(i);
        P=x->RECOELE_P->at(i);
        ETA=x->RECOELE_ETA->at(i);
        PHI=x->RECOELE_PHI->at(i);
        CHARGE=x->RECOELE_CHARGE->at(i);
        ID=x->RECOELE_ID->at(i);
        PFX_dB=x->RECOELE_PFX_dB->at(i);
        gsftrack_dxy=x->RECOELE_gsftrack_dxy->at(i);
        gsftrack_dz=x->RECOELE_gsftrack_dz->at(i);
        //Only filled if true. Seems to be only branch that isn't always filled among used branches
        if(i<x->RECOELE_isGap->size()) {
            isGap=x->RECOELE_isGap->at(i);
        } else {
            isGap=0;
        }
        PFchHad=x->RECOELE_PFchHad->at(i);
        PFneuHad=x->RECOELE_PFneuHad->at(i);
        PFphoton=x->RECOELE_PFphoton->at(i);
        PFPUchAllPart=x->RECOELE_PFPUchAllPart->at(i);
        SIP=x->RECOELE_SIP->at(i);
        IP=x->RECOELE_IP->at(i);
        IPERROR=x->RECOELE_IPERROR->at(i);
        mvaNonTrigV0=x->RECOELE_mvaNonTrigV0->at(i);
        loose=is_loose();
        tight=is_tight();
        loose_plus_SIP=loose && is_SIP();
    }
    four_vector fv() const {
        return four_vector(PT,ETA,PHI,0.000511);
    }
private:
    bool is_loose() const {
        return PT > 7. &&
               fabs(ETA) < 2.5 &&
               fabs(gsftrack_dxy) < .5 &&
               fabs(gsftrack_dz) < 1.;
    }
    bool is_tight() const {
        return is_loose() && ID==1; //BDT used for tight lepton (x->RECOELE_ID)
    }
    bool is_SIP() const {
        return fabs(SIP) < 4;
    }
};

struct lepton {
    lepton() : electron_(nullptr), is_a_muon(-1) {}
    void reset(muon* mu) {
        muon_=mu;
        is_a_muon=1;
    }
    void reset(electron* ele) {
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

struct photon{
//Not filled in ntuple, so no sense reading them
    //double x->RECOPFPHOT_PFchHad;
    //double x->RECOPFPHOT_PFneuHad;
    //double x->RECOPFPHOT_PFphoton;
    //double x->RECOPFPHOT_PFPUchAllPart;
    double PFX_rho;
    float PT;
    //float x->RECOPFPHOT_PTError;
    float ETA;
    float PHI;
    double deltaRoverPtSq;
    double deltaR;
    //double x->RECOPFPHOT_PT_uncorr;
    photon(ZpXanalyzer *x,int i) {
        //PFchHad=x->RECOPFPHOT_PFchHad->at(i);
        //PFneuHad=x->RECOPFPHOT_PFneuHad->at(i);
        //PFphoton=x->RECOPFPHOT_PFphoton->at(i);
        //PFPUchAllPart=x->RECOPFPHOT_PFPUchAllPart->at(i);
        PFX_rho=x->RECOPFPHOT_PFX_rho->at(i);
        PT=x->RECOPFPHOT_PT->at(i);
        //PTError=x->RECOPFPHOT_PTError->at(i);
        ETA=x->RECOPFPHOT_ETA->at(i);
        PHI=x->RECOPFPHOT_PHI->at(i);
        //PT_uncorr=x->RECOPFPHOT_PT_uncorr->at(i);
    }
    four_vector fv() const {
        return four_vector(PT,ETA,PHI,0.);
    }
};

#include <cmath>

template <typename T>
constexpr T reduceRange(T x) {
  constexpr T o2pi = 1. / (2. * M_PI);
  if (std::abs(x) <= T(M_PI))
    return x;
  T n = std::round(x * o2pi);
  return x - n * T(2. * M_PI);
}

inline float DELTAPHI(float a,float b) {
    return reduceRange(a-b);
}
