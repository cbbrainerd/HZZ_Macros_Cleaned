#include "NewNtuple.h"
#include "Math/Vector4D.h"
#include "TH1I.h"

//MC weights
#include "pileup_corrector.h"
#include "scale_factors.h"

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
    bool verbose;
    Bool_t Notify() override;
public:
    void Loop() override;
    TTree* tree_out;
    float tree_Zmass_nofsr;
    float tree_Zmass_fsr;
    bool is_MC;
    bkg_type tree_bkg_type;
    four_vector *tree_lep1_4v;
    four_vector *tree_lep2_4v;
    four_vector *tree_lep3_4v;
    four_vector *tree_lep1_fsr_4v;
    four_vector *tree_lep2_fsr_4v;
    four_vector *tree_lep3_fsr_4v;
    int cutflow_debug;
    int tree_lep1_q;
    int tree_lep2_q;
    int tree_lep3_q;
    int tree_lep3_tight;
    TFile* outfile;
    TH1* type_hist;
    TH1* cutflow_hist;
    void enable_branches();
    //Pileup reweighting and scale factors (NB: kfactor weight is not applied to data or WZ so we can skip it)
    pileup_corrector pileup_corr;
    scale_factors scale_factors_mu;
    scale_factors_and_efficiencies scale_factors_ele;
    float pileup_weight;
    float scale_factor[3];
    float efficiency[3];
    float total_weight;
    float drs12,drs13,drs23;
    std::size_t f_nleptons, f_ntight, f_no_os, f_ntight_iso, f_no_os_iso, f_ghost_removal;

    ZpXanalyzer(TTree* intree,TFile* outfile_,bool isMC,std::string era) : 
        NewNtuple(intree), 
        outfile(outfile_),
        is_MC(isMC),
        pileup_corr(isMC,era), 
        scale_factors_mu(isMC,era,std::string("muon")),
        scale_factors_ele(isMC,era,std::string("electron")),
        f_nleptons(0),
        f_ntight(0),
        f_no_os(0),
        f_ntight_iso(0),
        f_no_os_iso(0),
        f_ghost_removal(0)
    {
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
        tree_out->Branch("cutflow_debug",&cutflow_debug);
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
        tree_out->Branch("deltaR2_12",&drs12);
        tree_out->Branch("deltaR2_13",&drs13);
        tree_out->Branch("deltaR2_23",&drs23);
        if(isMC) {
            tree_out->Branch("pileup_weight",&pileup_weight);
            tree_out->Branch("scale_factor_1",&scale_factor[0]);
            tree_out->Branch("scale_factor_2",&scale_factor[1]);
            tree_out->Branch("scale_factor_3",&scale_factor[2]);
            tree_out->Branch("efficiency_1",&efficiency[0]);
            tree_out->Branch("efficiency_2",&efficiency[1]);
            tree_out->Branch("efficiency_3",&efficiency[2]);
            tree_out->Branch("total_weight",&total_weight);
        }
    }  
    inline void Write() {
        tree_out->Write();
        type_hist->Write();
        cutflow_hist->Write();
    }
};

//Particle struct definition
#include "particle.h"

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
