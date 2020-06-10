#include "NewNtuple.h"
#include "Math/Vector4D.h"
#include "TH1I.h"

//MC weights
#include "pileup_corrector.h"
#include "scale_factors.h"

enum cuts {
    cut_ERROR=-1,
    cut_NLEPTONS=0,
    cut_NLEPTONS_TIGHT,
    cut_NOZCAND,
    cut_NO_OSSF_PAIRS,
    cut_NO_ZZ_SELECTED,   
    cut_PASS,
    cut_NCUTS
};

enum channel_type : unsigned char {
    c_unknown=0,
    c_4e=1, 
    c_2e2mu=2,
    c_4mu=3
};

typedef ROOT::Math::PtEtaPhiMVector four_vector;

class Analyzer : public NewNtuple {
    int analyze();
    bool verbose;
    Bool_t Notify() override;
public:
    void Loop() override;
    void enable_branches();
    //Class members
    //MC/data scale factor classes
    pileup_corrector pileup_corr;
    scale_factors scale_factors_mu;
    scale_factors_and_efficiencies scale_factors_ele;
    //Output hists
    TH1 *all_Z_masses;
    TH1 *cutflow_hist;
    //Output file and tree
    TFile* outfile;
    TTree* tree_out;
    //Branch members
    int tree_cutflow;
    float tree_leading_lep_pt;
    float tree_subleading_lep_pt;
    float tree_leading_Z_mass_nofsr;
    float tree_leading_Z_mass_fsr;
    int tree_leading_Z_type;
    float tree_Z1_mass_nofsr;
    float tree_Z1_mass_fsr;
    int tree_Z1_type;
    float tree_Z2_mass_nofsr;
    float tree_Z2_mass_fsr;
    int tree_Z2_type;
    channel_type tree_channel_type;
    float tree_pileup_weight;
    float tree_Z_step_weight;
    float tree_ZZ_step_weight;
    bool is_MC;

    Analyzer(TTree* intree,TFile* outfile_,bool isMC,std::string era) : 
        NewNtuple(intree), 
        outfile(outfile_),
        is_MC(isMC),
        pileup_corr(isMC,era), 
        scale_factors_mu(isMC,era,std::string("muon")),
        scale_factors_ele(isMC,era,std::string("electron"))
    {
        //Disable all branches for performance reasons
        intree->SetBranchStatus("*",0);
        Init(intree);
        //Enable used branches. Be careful to add any additional branches to ZpXBranches.C!
        enable_branches();
        outfile->cd();
        //Histograms
        //type_hist=new TH1I("Bkg Types","Bkg Types",4,-.5,3.5);
        //const char *types[4]={"Z(mu+mu)+mu","Z(mu+mu)+e","Z(e+e)+mu","Z(e+e)+e"};
        //for(int i=1;i<=4;++i) type_hist->GetXaxis()->SetBinLabel(i,types[i-1]);

        //Output tree
        tree_out=new TTree("HZZTree","HZZTree");
        //Branches
        tree_out->Branch("Run",&Run);
        tree_out->Branch("Event",&Event);
        tree_out->Branch("LumiSection",&LumiSection);
        tree_out->Branch("num_PU_vertices",&num_PU_vertices);
        tree_out->Branch("PU_BunchCrossing",&PU_BunchCrossing);
        tree_out->Branch("MC_weighting",&MC_weighting);
        tree_out->Branch("cutflow_number",&tree_cutflow);
        tree_out->Branch("leading_lep_pt",&tree_leading_lep_pt);
        tree_out->Branch("subleading_lep_pt",&tree_subleading_lep_pt);
        tree_out->Branch("leading_Z_mass_nofsr",&tree_leading_Z_mass_nofsr);
        tree_out->Branch("leading_Z_mass_fsr",&tree_leading_Z_mass_fsr);
        tree_out->Branch("leading_Z_type",&tree_leading_Z_type);
        tree_out->Branch("Z1_mass_nofsr",&tree_Z1_mass_nofsr);
        tree_out->Branch("Z1_mass_fsr",&tree_Z1_mass_fsr);
        tree_out->Branch("Z1_type",&tree_Z1_type);
        tree_out->Branch("Z2_mass_nofsr",&tree_Z2_mass_nofsr);
        tree_out->Branch("Z2_mass_fsr",&tree_Z2_mass_fsr);
        tree_out->Branch("Z2_type",&tree_Z2_type);
        tree_out->Branch("pfmet",&RECO_PFMET);
        tree_out->Branch("pfmet_phi",&RECO_PFMET_PHI);
        tree_out->Branch("puppimet",&RECO_PUPPIMET);
        tree_out->Branch("puppimet_phi",&RECO_PUPPIMET_PHI);
        tree_out->Branch("pfmet_xycorr",&RECO_PFMET_xycorr);
        tree_out->Branch("pfmet_xycorr_phi",&RECO_PFMET_PHI_xycorr);
        if(isMC) {
            tree_out->Branch("pileup_weight",&tree_pileup_weight);
            tree_out->Branch("Z_step_weight",&tree_Z_step_weight);
            tree_out->Branch("ZZ_step_weight",&tree_ZZ_step_weight);
        }
    }  
    inline void Write() {
        tree_out->Write();
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
