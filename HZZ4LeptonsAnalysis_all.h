//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun  7 12:14:04 2012 by ROOT version 5.32/00
// from TTree HZZ4LeptonsAnalysis/HZZ4Leptons Analysis Tree
// found on file: roottree_leptons_Fall11_0706.root
//////////////////////////////////////////////////////////

#ifndef HZZ4LeptonsAnalysis_h
#define HZZ4LeptonsAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>

#if product_4mu

#elif product_2e2mu

#elif product_4e

#else
#error "Product must be one of \"4mu\", \"2e2mu\", and \"4e\".\n" 
#endif

#include "good_lumis.h" //Forgot to run over golden json for ntupler: do it now

//using namespace std;

//Parse the command line options
//Takes various formats
/*
class options {
public:
    bool valid;
    std::string output_file_name;
    std::string MC_year;
    std::string Data_year;
//    std::string 
    options(int argc,char **argv) {
        std::vector<std::string> opts;
        for(int x=0;x<argc;++x) opts.push_back(argv[x]);
        bool next_arg=false;
        for(auto const &opt : opts) {
            size_t equals_position=is_equal_opt(opt);
            if is_short_opt(opt) {
                parse_short_opt(opt,equals_position);
            } else if is_long_opt(opt) {
                parse_long_opt(opt,equals_position)
            } else if(equals_position) {
                parse_equals_opt(opt,equals_position);
            } else {
                usage();
                valid=false;
                return;
            }
        }
    }
private:
    bool is_short_opt(std::string opt) const { return opt.size() > 1 && opt[0]=='-' && opt[1] != '-'; }
    bool is_long_opt(std::string opt) const { return opt.size() > 1 && opt[0]=='-' && opt[1]=='-' }
    size_t is_equal_opt(std::string opt) const { size_t pos=opt.find('=',1); if(pos==std::string::npos) { return 0; } else { return pos; } }
    bool parse_short_opt(std::string opt,size_t pos) {
        
    }
}*/
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
class HZZ4LeptonsAnalysis {
public:
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          Run;
   UInt_t          Event;
   UInt_t          LumiSection;
   Float_t         Avginstlumi;
   Int_t           num_PU_vertices;
   Int_t           PU_BunchCrossing;
   Float_t         MC_weighting;
   Int_t           RECO_nMuHLTMatch;
   std::vector<float>   *RECOMU_PT_MuHLTMatch;
   std::vector<float>   *RECOMU_ETA_MuHLTMatch;
   std::vector<float>   *RECOMU_PHI_MuHLTMatch;
   Int_t           RECO_nEleHLTMatch;
   std::vector<float>   *RECOELE_PT_EleHLTMatch;
   std::vector<float>   *RECOELE_ETA_EleHLTMatch;
   std::vector<float>   *RECOELE_PHI_EleHLTMatch;
   Char_t          HLTPathsFired;
   std::vector<float>   *MC_E;
   std::vector<float>   *MC_PT;
   std::vector<float>   *MC_ETA;
   std::vector<float>   *MC_THETA;
   std::vector<float>   *MC_PHI;
   std::vector<float>   *MC_MASS;
   std::vector<float>   *MC_PDGID;
   std::vector<float>   *MC_LEPT_PT;
   std::vector<float>   *MC_LEPT_ETA;
   std::vector<float>   *MC_LEPT_PHI;
   std::vector<float>   *MC_LEPT_THETA;
   std::vector<float>   *MC_LEPT_PDGID;
   std::vector<std::vector<float> > *MC_Z_PT;
   std::vector<std::vector<float> > *MC_Z_ETA;
   std::vector<std::vector<float> > *MC_Z_PHI;
   std::vector<std::vector<float> > *MC_Z_THETA;
   std::vector<std::vector<float> > *MC_Z_MASS;
   std::vector<std::vector<float> > *MC_Z_PDGID;
   std::vector<std::vector<float> > *MC_fourl_MASS;
   std::vector<std::vector<float> > *MC_fourl_PT;
   std::vector<std::vector<float> > *MC_fourl_PDGID;
   std::vector<float>   *MC_ZZ_MASS;
   std::vector<float>   *MC_ZZ_PT;
   std::vector<float>   *MC_ZZ_ETA;
   std::vector<float>   *MC_ZZ_PHI;
   std::vector<float>   *MC_ZZ_THETA;
   std::vector<float>   *MC_ZZ_PDGID;
   std::vector<std::vector<float> > *MC_ZZ_D_MASS;
   std::vector<std::vector<float> > *MC_ZZ_D_PT;
   std::vector<std::vector<float> > *MC_ZZ_D_ETA;
   std::vector<std::vector<float> > *MC_ZZ_D_PHI;
   std::vector<std::vector<float> > *MC_ZZ_D_THETA;
   std::vector<std::vector<float> > *MC_ZZ_D_PDGID;
   std::vector<std::vector<std::vector<float> > > *MC_ZZ_DD_MASS;
   std::vector<std::vector<std::vector<float> > > *MC_ZZ_DD_PT;
   std::vector<std::vector<std::vector<float> > > *MC_ZZ_DD_ETA;
   std::vector<std::vector<std::vector<float> > > *MC_ZZ_DD_PHI;
   std::vector<std::vector<std::vector<float> > > *MC_ZZ_DD_THETA;
   std::vector<std::vector<std::vector<float> > > *MC_ZZ_DD_PDGID;
   std::vector<float>   *MC_GENJET_PT;
   std::vector<float>   *MC_GENJET_ETA;
   std::vector<float>   *MC_GENJET_PHI;
   Float_t         MC_GENMET;
   std::vector<double>  *RECORF_2e2mu_cosTheta1_spin;
   std::vector<double>  *RECORF_2e2mu_cosTheta2_spin;
   std::vector<double>  *RECORF_2e2mu_cosThetaStar_spin;
   std::vector<double>  *RECORF_2e2mu_Phi_spin;
   std::vector<double>  *RECORF_2e2mu_Phi1_spin;
   std::vector<double>  *RECORF_2e2mu_Phi2_spin;
   std::vector<double>  *RECORF_2e2mu_phi1RF_spin;
   std::vector<double>  *RECORF_2e2mu_phi2RF_spin;
   std::vector<double>  *RECORF_2e2mu_MELA;
   std::vector<double>  *RECORF_4e_cosTheta1_spin;
   std::vector<double>  *RECORF_4e_cosTheta2_spin;
   std::vector<double>  *RECORF_4e_cosThetaStar_spin;
   std::vector<double>  *RECORF_4e_Phi_spin;
   std::vector<double>  *RECORF_4e_Phi1_spin;
   std::vector<double>  *RECORF_4e_Phi2_spin;
   std::vector<double>  *RECORF_4e_phi1RF_spin;
   std::vector<double>  *RECORF_4e_phi2RF_spin;
   std::vector<double>  *RECORF_4e_MELA;
   std::vector<double>  *RECORF_4mu_cosTheta1_spin;
   std::vector<double>  *RECORF_4mu_cosTheta2_spin;
   std::vector<double>  *RECORF_4mu_cosThetaStar_spin;
   std::vector<double>  *RECORF_4mu_Phi_spin;
   std::vector<double>  *RECORF_4mu_Phi1_spin;
   std::vector<double>  *RECORF_4mu_Phi2_spin;
   std::vector<double>  *RECORF_4mu_phi1RF_spin;
   std::vector<double>  *RECORF_4mu_phi2RF_spin;
   std::vector<double>  *RECORF_4mu_MELA;
   std::vector<float>   *RECO_ZMM_MASS;
   std::vector<float>   *RECO_ZEE_MASS;
   std::vector<float>   *RECO_DiLep_MASS;
   std::vector<std::vector<float> > *RECO_ZMM_PT;
   std::vector<std::vector<float> > *RECO_ZEE_PT;
   std::vector<std::vector<float> > *RECO_DiLep_PT;
   std::vector<std::vector<float> > *RECO_ZMM_ETA;
   std::vector<std::vector<float> > *RECO_ZEE_ETA;
   std::vector<std::vector<float> > *RECO_DiLep_ETA;
   std::vector<std::vector<float> > *RECO_ZMM_PHI;
   std::vector<std::vector<float> > *RECO_ZEE_PHI;
   std::vector<std::vector<float> > *RECO_DiLep_PHI;
   std::vector<float>   *RECO_ZMMss_MASS;
   std::vector<float>   *RECO_ZEEss_MASS;
   std::vector<float>   *RECO_ZEM_MASS;
   std::vector<std::vector<float> > *RECO_ZMMss_PT;
   std::vector<std::vector<float> > *RECO_ZEEss_PT;
   std::vector<std::vector<float> > *RECO_ZEM_PT;
   std::vector<std::vector<float> > *RECO_ZMMss_ETA;
   std::vector<std::vector<float> > *RECO_ZEEss_ETA;
   std::vector<std::vector<float> > *RECO_ZEM_ETA;
   std::vector<std::vector<float> > *RECO_ZMMss_PHI;
   std::vector<std::vector<float> > *RECO_ZEEss_PHI;
   std::vector<std::vector<float> > *RECO_ZEM_PHI;
   std::vector<std::vector<float> > *RECO_MMMM_MASS;
   std::vector<std::vector<float> > *RECO_MMMM_PT;
   std::vector<std::vector<float> > *RECO_MMMM_ETA;
   std::vector<std::vector<float> > *RECO_MMMM_PHI;
   std::vector<float>   *RECO_MMMM_MASS_REFIT;
   std::vector<std::vector<float> > *RECO_EEEE_MASS;
   std::vector<std::vector<float> > *RECO_EEEE_PT;
   std::vector<std::vector<float> > *RECO_EEEE_ETA;
   std::vector<std::vector<float> > *RECO_EEEE_PHI;
   std::vector<float>   *RECO_EEEE_MASS_REFIT;
   std::vector<std::vector<float> > *RECO_EEMM_MASS;
   std::vector<std::vector<float> > *RECO_EEMM_PT;
   std::vector<std::vector<float> > *RECO_EEMM_ETA;
   std::vector<std::vector<float> > *RECO_EEMM_PHI;
   std::vector<float>   *RECO_EEMM_MASS_REFIT;
   std::vector<float>   *RECO_LLL0_MASS;
   std::vector<float>   *RECO_LLL1_MASS;
   std::vector<float>   *RECO_LLL2_MASS;
   std::vector<float>   *RECO_LLL3_MASS;
   std::vector<std::vector<float> > *RECO_LLL0_PT;
   std::vector<std::vector<float> > *RECO_LLL1_PT;
   std::vector<std::vector<float> > *RECO_LLL2_PT;
   std::vector<std::vector<float> > *RECO_LLL3_PT;
   std::vector<float>   *RECO_LLLl0_MASS;
   std::vector<float>   *RECO_LLLl1_MASS;
   std::vector<std::vector<float> > *RECO_LLLl0_PT;
   std::vector<std::vector<float> > *RECO_LLLl1_PT;
   std::vector<float>   *RECO_LLLL0ss_MASS;
   std::vector<std::vector<float> > *RECO_LLLL0ss_PT;
   std::vector<float>   *RECO_LLLL1ss_MASS;
   std::vector<std::vector<float> > *RECO_LLLL1ss_PT;
   std::vector<float>   *RECO_LLLL2ss_MASS;
   std::vector<std::vector<float> > *RECO_LLLL2ss_PT;
   std::vector<std::vector<float> > *RECO_LLLL_MASS;
   std::vector<std::vector<float> > *RECO_LLLL_PT;
   std::vector<std::vector<float> > *RECO_LLLL_ETA;
   std::vector<std::vector<float> > *RECO_LLLL_PHI;
   std::vector<float>   *RECOELE_E;
   std::vector<float>   *RECOELE_PT;
   std::vector<float>   *RECOELE_PTError;
   std::vector<float>   *RECOELE_P;
   std::vector<float>   *RECOELE_ETA;
   std::vector<float>   *RECOELE_THETA;
   std::vector<float>   *RECOELE_PHI;
   std::vector<float>   *RECOELE_MASS;
   std::vector<float>   *RECOELE_CHARGE;
   std::vector<float>   *RECOELE_ID;
   std::vector<float>   *RECOELE_PT_uncorr;
   std::vector<int>     *RECOELE_isEcalDriven;
   std::vector<int>     *RECOELE_isTrackerDriven;
   std::vector<float>   *RECOELE_gsftrack_NPixHits;
   std::vector<float>   *RECOELE_gsftrack_NStripHits;
   std::vector<float>   *RECOELE_gsftrack_chi2;
   std::vector<float>   *RECOELE_gsftrack_dxyB;
   std::vector<float>   *RECOELE_gsftrack_dxy;
   std::vector<float>   *RECOELE_gsftrack_dxyError;
   std::vector<float>   *RECOELE_gsftrack_dzB;
   std::vector<float>   *RECOELE_gsftrack_dz;
   std::vector<float>   *RECOELE_gsftrack_dzError;
   std::vector<int>     *RECOELE_gsftrack_losthits;
   std::vector<int>     *RECOELE_gsftrack_validhits;
   std::vector<int>     *RECOELE_gsftrack_expected_inner_hits;
   std::vector<float>   *RECOELE_scl_E;
   std::vector<float>   *RECOELE_scl_Et;
   std::vector<float>   *RECOELE_scl_Eta;
   std::vector<float>   *RECOELE_scl_Phi;
   std::vector<float>   *RECOELE_ecalEnergy;
   std::vector<float>   *RECOELE_ep;
   std::vector<float>   *RECOELE_eSeedp;
   std::vector<float>   *RECOELE_eSeedpout;
   std::vector<float>   *RECOELE_eElepout;
   std::vector<float>   *RECOELE_deltaEtaIn;
   std::vector<float>   *RECOELE_deltaEtaSeed;
   std::vector<float>   *RECOELE_deltaEtaEle;
   std::vector<float>   *RECOELE_deltaPhiIn;
   std::vector<float>   *RECOELE_deltaPhiSeed;
   std::vector<float>   *RECOELE_deltaPhiEle;
   std::vector<int>     *RECOELE_isbarrel;
   std::vector<int>     *RECOELE_isendcap;
   std::vector<int>     *RECOELE_isGap;
   std::vector<int>     *RECOELE_isEBetaGap;
   std::vector<int>     *RECOELE_isEBphiGap;
   std::vector<int>     *RECOELE_isEEdeeGap;
   std::vector<int>     *RECOELE_isEEringGap;
   std::vector<float>   *RECOELE_sigmaIetaIeta;
   std::vector<float>   *RECOELE_sigmaEtaEta;
   std::vector<float>   *RECOELE_e15;
   std::vector<float>   *RECOELE_e25max;
   std::vector<float>   *RECOELE_e55;
   std::vector<float>   *RECOELE_he;
   std::vector<float>   *RECOELE_r9;
   std::vector<float>   *RECOELE_mva;
   std::vector<float>   *RECOELE_fbrem;
   std::vector<int>     *RECOELE_nbrems;
   std::vector<int>     *RECOELE_Class;
   std::vector<float>   *RECOELE_fbrem_mode;
   std::vector<float>   *RECOELE_fbrem_mean;
   std::vector<float>   *RECOELE_EGMTRACKISO;
   std::vector<float>   *RECOELE_EGMHCALISO;
   std::vector<float>   *RECOELE_EGMECALISO;
   std::vector<float>   *RECOELE_EGMX;
   std::vector<double>  *RECOELE_PFchAllPart;
   std::vector<double>  *RECOELE_PFchHad;
   std::vector<double>  *RECOELE_PFneuHad;
   std::vector<double>  *RECOELE_PFphoton;
   std::vector<double>  *RECOELE_PFPUchAllPart;
   std::vector<double>  *RECOELE_PFX_dB;
   std::vector<double>  *RECOELE_PFX_rho;
   std::vector<double>  *RECOELE_regEnergy;
   std::vector<double>  *RECOELE_regEnergyError;
   std::vector<float>   *RECOELE_SIP;
   std::vector<float>   *RECOELE_IP;
   std::vector<float>   *RECOELE_IPERROR;
   std::vector<float>   *RECOELE_SIP_KF;
   std::vector<float>   *RECOELE_IP_KF;
   std::vector<float>   *RECOELE_IPERROR_KF;
   std::vector<float>   *RECOELE_SIP_GD;
   std::vector<float>   *RECOELE_SIP_GDEEEE;
   std::vector<float>   *RECOELE_SIP_Std;
   std::vector<float>   *RECOELE_SIP_StdEEEE;
   std::vector<float>   *RECOELE_SIP_Kin;
   std::vector<float>   *RECOELE_SIP_KinEEEE;
   std::vector<float>   *RECOELE_STIP;
   std::vector<float>   *RECOELE_SLIP;
   std::vector<float>   *RECOELE_TIP;
   std::vector<float>   *RECOELE_LIP;
   std::vector<float>   *RECOELE_TIPERROR;
   std::vector<float>   *RECOELE_LIPERROR;
   std::vector<double>  *RECOELE_sclRawE;
   std::vector<double>  *RECOELE_sclX;
   std::vector<double>  *RECOELE_sclY;
   std::vector<double>  *RECOELE_sclZ;
   std::vector<int>     *RECOELE_seedSubdet1;
   std::vector<double>  *RECOELE_seedDphi1;
   std::vector<double>  *RECOELE_seedDrz1;
   std::vector<int>     *RECOELE_seedSubdet2;
   std::vector<double>  *RECOELE_seedDphi2;
   std::vector<double>  *RECOELE_seedDrz2;
   std::vector<double>  *RECOELE_eidVeryLoose;
   std::vector<double>  *RECOELE_eidLoose;
   std::vector<double>  *RECOELE_eidMedium;
   std::vector<double>  *RECOELE_eidTight;
   std::vector<double>  *RECOELE_eidHZZVeryLoose;
   std::vector<double>  *RECOELE_eidHZZLoose;
   std::vector<double>  *RECOELE_eidHZZMedium;
   std::vector<double>  *RECOELE_eidHZZTight;
   std::vector<double>  *RECOELE_mvaTrigV0;
   std::vector<double>  *RECOELE_mvaNonTrigV0;
   std::vector<std::vector<std::vector<double> > > *RECOELE_COV;
   std::vector<float>   *RECOELE_TLE_ParentSC_X;
   std::vector<float>   *RECOELE_TLE_ParentSC_Y;
   std::vector<float>   *RECOELE_TLE_ParentSC_Z;
   std::vector<float>   *RECOELE_ecalTrkEnergyPreCorr;
   std::vector<float>   *RECOELE_ecalTrkEnergyErrPreCorr;
   std::vector<float>   *RECOELE_ecalTrkEnergyErrPostCorr;
   std::vector<float>   *RECOELE_energyScaleValue;
   std::vector<float>   *RECOELE_energySigmaValue;
   std::vector<float>   *RECOELE_energyScaleUp;
   std::vector<float>   *RECOELE_energyScaleDown;
   std::vector<float>   *RECOELE_energyScaleStatUp;
   std::vector<float>   *RECOELE_energyScaleStatDown;
   std::vector<float>   *RECOELE_energyScaleSystUp;
   std::vector<float>   *RECOELE_energyScaleSystDown;
   std::vector<float>   *RECOELE_energyScaleGainUp;
   std::vector<float>   *RECOELE_energyScaleGainDown;
   std::vector<float>   *RECOELE_energyScaleEtUp;
   std::vector<float>   *RECOELE_energyScaleEtDown;
   std::vector<float>   *RECOELE_energySigmaUp;
   std::vector<float>   *RECOELE_energySigmaDown;
   std::vector<float>   *RECOELE_energySigmaPhiUp;
   std::vector<float>   *RECOELE_energySigmaPhiDown;
   std::vector<float>   *RECOELE_energySigmaRhoUp;
   std::vector<float>   *RECOELE_energySigmaRhoDown;
   std::vector<int>     *RECOMU_isPFMu;
   std::vector<int>     *RECOMU_isGlobalMu;
   std::vector<int>     *RECOMU_isStandAloneMu;
   std::vector<int>     *RECOMU_isTrackerMu;
   std::vector<int>     *RECOMU_isCaloMu;
   std::vector<int>     *RECOMU_isTrackerHighPtMu;
   std::vector<int>     *RECOMU_isME0Muon;
   std::vector<float>   *RECOMU_E;
   std::vector<float>   *RECOMU_PT;
   std::vector<float>   *RECOMU_P;
   std::vector<float>   *RECOMU_ETA;
   std::vector<float>   *RECOMU_THETA;
   std::vector<float>   *RECOMU_PHI;
   std::vector<float>   *RECOMU_MASS;
   std::vector<float>   *RECOMU_CHARGE;
   std::vector<double>  *RECOMU_PT_uncorr;
   std::vector<std::vector<std::vector<double> > > *RECOMU_COV;
   std::vector<float>   *RECOMU_TRACKISO;
   std::vector<float>   *RECOMU_TRACKISO_SUMPT;
   std::vector<float>   *RECOMU_HCALISO;
   std::vector<float>   *RECOMU_ECALISO;
   std::vector<float>   *RECOMU_X;
   std::vector<double>  *RECOMU_PFchHad;
   std::vector<double>  *RECOMU_PFneuHad;
   std::vector<double>  *RECOMU_PFphoton;
   std::vector<double>  *RECOMU_PFPUchAllPart;
   std::vector<double>  *RECOMU_PFX_dB;
   std::vector<double>  *RECOMU_PFX_rho;
   std::vector<double>  *RECOPFPHOT_PFchHad;
   std::vector<double>  *RECOPFPHOT_PFneuHad;
   std::vector<double>  *RECOPFPHOT_PFphoton;
   std::vector<double>  *RECOPFPHOT_PFPUchAllPart;
   std::vector<double>  *RECOPFPHOT_PFX_rho;
   std::vector<float>   *RECOMU_SIP;
   std::vector<float>   *RECOMU_IP;
   std::vector<float>   *RECOMU_IPERROR;
   std::vector<float>   *RECOMU_SIP_KF;
   std::vector<float>   *RECOMU_IP_KF;
   std::vector<float>   *RECOMU_IPERROR_KF;
   std::vector<float>   *RECOMU_SIP_GD;
   std::vector<float>   *RECOMU_SIP_GDMMMM;
   std::vector<float>   *RECOMU_SIP_Std;
   std::vector<float>   *RECOMU_SIP_StdMMMM;
   std::vector<float>   *RECOMU_SIP_Kin;
   std::vector<float>   *RECOMU_SIP_KinMMMM;
   std::vector<float>   *RECOMU_STIP;
   std::vector<float>   *RECOMU_SLIP;
   std::vector<float>   *RECOMU_TIP;
   std::vector<float>   *RECOMU_LIP;
   std::vector<float>   *RECOMU_TIPERROR;
   std::vector<float>   *RECOMU_LIPERROR;
   std::vector<float>   *RECOMU_caloCompatibility;
   std::vector<float>   *RECOMU_segmentCompatibility;
   std::vector<int>     *RECOMU_numberOfMatches;
   std::vector<int>     *RECOMU_numberOfMatchedStations;
   std::vector<int>     *RECOMU_glbmuPromptTight;
   std::vector<int>     *RECOMU_trkmuArbitration;
   std::vector<int>     *RECOMU_trkmu2DCompatibilityLoose;
   std::vector<int>     *RECOMU_trkmu2DCompatibilityTight;
   std::vector<int>     *RECOMU_trkmuOneStationLoose;
   std::vector<int>     *RECOMU_trkmuOneStationTight;
   std::vector<int>     *RECOMU_trkmuLastStationLoose;
   std::vector<int>     *RECOMU_trkmuLastStationTight;
   std::vector<int>     *RECOMU_trkmuOneStationAngLoose;
   std::vector<int>     *RECOMU_trkmuOneStationAngTight;
   std::vector<int>     *RECOMU_trkmuLastStationAngLoose;
   std::vector<int>     *RECOMU_trkmuLastStationAngTight;
   std::vector<int>     *RECOMU_trkmuLastStationOptimizedLowPtLoose;
   std::vector<int>     *RECOMU_trkmuLastStationOptimizedLowPtTight;
   std::vector<float>   *RECOMU_mutrkPT;
   std::vector<float>   *RECOMU_mutrkPTError;
   std::vector<float>   *RECOMU_mutrkDxy;
   std::vector<float>   *RECOMU_mutrkDxyError;
   std::vector<float>   *RECOMU_mutrkDxyB;
   std::vector<float>   *RECOMU_mutrkDz;
   std::vector<float>   *RECOMU_mutrkDzError;
   std::vector<float>   *RECOMU_mutrkDzB;
   std::vector<float>   *RECOMU_mutrkChi2PerNdof;
   std::vector<float>   *RECOMU_mutrkCharge;
   std::vector<float>   *RECOMU_mutrkNHits;
   std::vector<float>   *RECOMU_mutrkNStripHits;
   std::vector<float>   *RECOMU_mutrkNPixHits;
   std::vector<float>   *RECOMU_mutrkNMuonHits;
   std::vector<float>   *RECOMU_mutrktrackerLayersWithMeasurement;
   std::vector<float>   *RECOMU_muInnertrkDxy;
   std::vector<float>   *RECOMU_muInnertrkDxyError;
   std::vector<float>   *RECOMU_muInnertrkDxyB;
   std::vector<float>   *RECOMU_muInnertrkDz;
   std::vector<float>   *RECOMU_muInnertrkDzError;
   std::vector<float>   *RECOMU_muInnertrkDzB;
   std::vector<float>   *RECOMU_muInnertrkChi2PerNdof;
   std::vector<float>   *RECOMU_muInnertrktrackerLayersWithMeasurement;
   std::vector<float>   *RECOMU_muInnertrkPT;
   std::vector<float>   *RECOMU_muInnertrkPTError;
   std::vector<float>   *RECOMU_muInnertrkCharge;
   std::vector<float>   *RECOMU_muInnertrkNHits;
   std::vector<float>   *RECOMU_muInnertrkNStripHits;
   std::vector<float>   *RECOMU_muInnertrkNPixHits;
   std::vector<int>     *RECOMU_mubesttrkType;
   std::vector<float>   *RECOMU_mubesttrkDxy;
   std::vector<float>   *RECOMU_mubesttrkDxyError;
   std::vector<float>   *RECOMU_mubesttrkDz;
   std::vector<float>   *RECOMU_mubesttrkDzError;
   std::vector<float>   *RECOMU_mubesttrkPTError;
   std::vector<float>   *RECOMU_Rochester_Error;
   std::vector<double>  *ftsigma;
   std::vector<double>  *gdX;
   std::vector<double>  *gdY;
   std::vector<double>  *gdZ;
   std::vector<double>  *ftsigmalag;
   std::vector<double>  *gdlagX;
   std::vector<double>  *gdlagY;
   std::vector<double>  *gdlagZ;
   std::vector<double>  *gdlagProb;
   std::vector<double>  *gdlagNdof;
   std::vector<double>  *ftsigmaMMMM;
   std::vector<double>  *gdXMMMM;
   std::vector<double>  *gdYMMMM;
   std::vector<double>  *gdZMMMM;
   std::vector<double>  *ftsigmalagMMMM;
   std::vector<double>  *gdlagXMMMM;
   std::vector<double>  *gdlagYMMMM;
   std::vector<double>  *gdlagZMMMM;
   std::vector<double>  *gdlagProbMMMM;
   std::vector<double>  *gdlagNdofMMMM;
   std::vector<double>  *ftsigmaEEEE;
   std::vector<double>  *gdXEEEE;
   std::vector<double>  *gdYEEEE;
   std::vector<double>  *gdZEEEE;
   std::vector<double>  *ftsigmalagEEEE;
   std::vector<double>  *gdlagXEEEE;
   std::vector<double>  *gdlagYEEEE;
   std::vector<double>  *gdlagZEEEE;
   std::vector<double>  *gdlagProbEEEE;
   std::vector<double>  *gdlagNdofEEEE;
   std::vector<double>  *StdFitVertexX;
   std::vector<double>  *StdFitVertexY;
   std::vector<double>  *StdFitVertexZ;
   std::vector<double>  *StdFitVertexChi2r;
   std::vector<double>  *StdFitVertexProb;
   std::vector<std::vector<float> > *StdFitVertexTrack_PT;
   std::vector<std::vector<float> > *StdFitVertexTrack_ETA;
   std::vector<std::vector<float> > *StdFitVertexTrack_PHI;
   std::vector<double>  *KinFitVertexX;
   std::vector<double>  *KinFitVertexY;
   std::vector<double>  *KinFitVertexZ;
   std::vector<double>  *KinFitVertexChi2r;
   std::vector<double>  *KinFitVertexProb;
   std::vector<double>  *StdFitVertexXMMMM;
   std::vector<double>  *StdFitVertexYMMMM;
   std::vector<double>  *StdFitVertexZMMMM;
   std::vector<double>  *StdFitVertexChi2rMMMM;
   std::vector<double>  *StdFitVertexProbMMMM;
   std::vector<std::vector<float> > *StdFitVertexTrackMMMM_PT;
   std::vector<std::vector<float> > *StdFitVertexTrackMMMM_ETA;
   std::vector<std::vector<float> > *StdFitVertexTrackMMMM_PHI;
   std::vector<double>  *KinFitVertexXMMMM;
   std::vector<double>  *KinFitVertexYMMMM;
   std::vector<double>  *KinFitVertexZMMMM;
   std::vector<double>  *KinFitVertexChi2rMMMM;
   std::vector<double>  *KinFitVertexProbMMMM;
   std::vector<double>  *StdFitVertexXEEEE;
   std::vector<double>  *StdFitVertexYEEEE;
   std::vector<double>  *StdFitVertexZEEEE;
   std::vector<double>  *StdFitVertexChi2rEEEE;
   std::vector<double>  *StdFitVertexProbEEEE;
   std::vector<std::vector<float> > *StdFitVertexTrackEEEE_PT;
   std::vector<std::vector<float> > *StdFitVertexTrackEEEE_ETA;
   std::vector<std::vector<float> > *StdFitVertexTrackEEEE_PHI;
   std::vector<double>  *KinFitVertexXEEEE;
   std::vector<double>  *KinFitVertexYEEEE;
   std::vector<double>  *KinFitVertexZEEEE;
   std::vector<double>  *KinFitVertexChi2rEEEE;
   std::vector<double>  *KinFitVertexProbEEEE;
   std::vector<double>  *StdFitVertexChi2rMMM;
   std::vector<double>  *StdFitVertexProbMMM;
   std::vector<double>  *StdFitVertexChi2rMME;
   std::vector<double>  *StdFitVertexProbMME;
   std::vector<double>  *StdFitVertexChi2rEEE;
   std::vector<double>  *StdFitVertexProbEEE;
   std::vector<double>  *StdFitVertexChi2rMEE;
   std::vector<double>  *StdFitVertexProbMEE;
   std::vector<double>  *StdFitVertexChi2rDiLep;
   std::vector<double>  *StdFitVertexProbDiLep;
   std::vector<float>   *ConvMapDist;
   std::vector<float>   *ConvMapDcot;
   std::vector<int>     *RECOMU_MatchingMCTruth;
   std::vector<float>   *RECOMU_MatchingMCpT;
   std::vector<float>   *RECOMU_MatchingMCEta;
   std::vector<float>   *RECOMU_MatchingMCPhi;
   std::vector<int>     *RECOELE_MatchingMCTruth;
   std::vector<float>   *RECOELE_MatchingMCpT;
   std::vector<float>   *RECOELE_MatchingMCEta;
   std::vector<float>   *RECOELE_MatchingMCPhi;
   std::vector<int>     *RECOPHOT_MatchingMCTruth;
   std::vector<float>   *RECOPHOT_MatchingMCpT;
   std::vector<float>   *RECOPHOT_MatchingMCEta;
   std::vector<float>   *RECOPHOT_MatchingMCPhi;
   std::vector<int>     *RECOzMuMu_MatchingMCTruth;
   std::vector<float>   *RECOzMuMu_MatchingMCpT;
   std::vector<float>   *RECOzMuMu_MatchingMCmass;
   std::vector<float>   *RECOzMuMu_MatchingMCEta;
   std::vector<float>   *RECOzMuMu_MatchingMCPhi;
   std::vector<int>     *RECOzEE_MatchingMCTruth;
   std::vector<float>   *RECOzEE_MatchingMCpT;
   std::vector<float>   *RECOzEE_MatchingMCmass;
   std::vector<float>   *RECOzEE_MatchingMCEta;
   std::vector<float>   *RECOzEE_MatchingMCPhi;
   std::vector<int>     *RECOHzzEEEE_MatchingMCTruth;
   std::vector<float>   *RECOHzzEEEE_MatchingMCpT;
   std::vector<float>   *RECOHzzEEEE_MatchingMCmass;
   std::vector<float>   *RECOHzzEEEE_MatchingMCEta;
   std::vector<float>   *RECOHzzEEEE_MatchingMCPhi;
   std::vector<int>     *RECOHzzEEMM_MatchingMCTruth;
   std::vector<float>   *RECOHzzEEMM_MatchingMCpT;
   std::vector<float>   *RECOHzzEEMM_MatchingMCmass;
   std::vector<float>   *RECOHzzEEMM_MatchingMCEta;
   std::vector<float>   *RECOHzzEEMM_MatchingMCPhi;
   std::vector<int>     *RECOHzzMMMM_MatchingMCTruth;
   std::vector<float>   *RECOHzzMMMM_MatchingMCpT;
   std::vector<float>   *RECOHzzMMMM_MatchingMCmass;
   std::vector<float>   *RECOHzzMMMM_MatchingMCEta;
   std::vector<float>   *RECOHzzMMMM_MatchingMCPhi;
   Int_t           RECO_NMU;
   Int_t           RECO_NELE;
   Int_t           RECO_NTRACK;
   Int_t           RECO_NPHOT;
   std::vector<float>   *RECOPHOT_PT;
   std::vector<float>   *RECOPHOT_ETA;
   std::vector<float>   *RECOPHOT_PHI;
   std::vector<float>   *RECOPHOT_THETA;
   std::vector<float>   *RECOPHOT_TLE_ParentSC_X;
   std::vector<float>   *RECOPHOT_TLE_ParentSC_Y;
   std::vector<float>   *RECOPHOT_TLE_ParentSC_Z;
   Int_t           RECO_NPFPHOT;
   std::vector<float>   *RECOPFPHOT_PT;
   std::vector<float>   *RECOPFPHOT_PTError;
   std::vector<float>   *RECOPFPHOT_ETA;
   std::vector<float>   *RECOPFPHOT_PHI;
   std::vector<float>   *RECOPFPHOT_THETA;
   std::vector<double>  *RECOPFPHOT_PT_uncorr;
   std::vector<float>   *RECOPFPHOT_ecalEnergyPreCorr;
   std::vector<float>   *RECOPFPHOT_ecalEnergyErrPreCorr;
   std::vector<float>   *RECOPFPHOT_ecalEnergyErrPostCorr;
   std::vector<float>   *RECOPFPHOT_energyScaleValue;
   std::vector<float>   *RECOPFPHOT_energySigmaValue;
   std::vector<float>   *RECOPFPHOT_energyScaleUp;
   std::vector<float>   *RECOPFPHOT_energyScaleDown;
   std::vector<float>   *RECOPFPHOT_energyScaleStatUp;
   std::vector<float>   *RECOPFPHOT_energyScaleStatDown;
   std::vector<float>   *RECOPFPHOT_energyScaleSystUp;
   std::vector<float>   *RECOPFPHOT_energyScaleSystDown;
   std::vector<float>   *RECOPFPHOT_energyScaleGainUp;
   std::vector<float>   *RECOPFPHOT_energyScaleGainDown;
   std::vector<float>   *RECOPFPHOT_energyScaleEtUp;
   std::vector<float>   *RECOPFPHOT_energyScaleEtDown;
   std::vector<float>   *RECOPFPHOT_energySigmaUp;
   std::vector<float>   *RECOPFPHOT_energySigmaDown;
   std::vector<float>   *RECOPFPHOT_energySigmaPhiUp;
   std::vector<float>   *RECOPFPHOT_energySigmaPhiDown;
   std::vector<float>   *RECOPFPHOT_energySigmaRhoUp;
   std::vector<float>   *RECOPFPHOT_energySigmaRhoDown;
   Double_t        BeamSpot_X;
   Double_t        BeamSpot_Y;
   Double_t        BeamSpot_Z;
   Int_t           RECO_NVTX;
   Float_t         RECO_VERTEX_x;
   Float_t         RECO_VERTEX_y;
   Float_t         RECO_VERTEX_z;
   Float_t         RECO_VERTEX_ndof;
   Float_t         RECO_VERTEX_chi2;
   Int_t           RECO_VERTEX_ntracks;
   Float_t         RECO_VERTEXPROB;
   Int_t           RECO_VERTEX_isValid;
   Int_t           RECO_PFJET_N;
   std::vector<int>     *RECO_PFJET_CHARGE;
   std::vector<float>   *RECO_PFJET_ET;
   std::vector<float>   *RECO_PFJET_PT;
   std::vector<float>   *RECO_PFJET_ETA;
   std::vector<float>   *RECO_PFJET_PHI;
   std::vector<int>     *RECO_PFJET_PUID_loose;
   std::vector<int>     *RECO_PFJET_PUID_medium;
   std::vector<int>     *RECO_PFJET_PUID;
   std::vector<float>   *RECO_PFJET_PUID_MVA;
   std::vector<float>   *RECO_PFJET_QG_Likelihood;
   std::vector<float>   *RECO_PFJET_QG_axis2;
   std::vector<float>   *RECO_PFJET_QG_ptd;
   std::vector<float>   *RECO_PFJET_QG_mult;
   Double_t        RHO_ele;
   Double_t        RHO_mu;
   Int_t           LHE_PARTON_N;
   std::vector<int>     *LHE_PARTON_CLEAR;
   std::vector<int>     *LHE_PARTON_PDGID;
   std::vector<float>   *LHE_PARTON_PT;
   std::vector<float>   *LHE_PARTON_ETA;
   std::vector<float>   *LHE_PARTON_PHI;
   std::vector<float>   *LHE_PARTON_E;
   std::vector<float>   *RECO_PFJET_PT_UncUp;
   std::vector<float>   *RECO_PFJET_PT_UncDn;
   std::vector<float>   *RECO_PFJET_AREA;
   std::vector<float>   *RECO_PFJET_PTD;
   std::vector<float>   *RECO_PFJET_CHARGED_HADRON_ENERGY;
   std::vector<float>   *RECO_PFJET_NEUTRAL_HADRON_ENERGY;
   std::vector<float>   *RECO_PFJET_PHOTON_ENERGY;
   std::vector<float>   *RECO_PFJET_ELECTRON_ENERGY;
   std::vector<float>   *RECO_PFJET_MUON_ENERGY;
   std::vector<float>   *RECO_PFJET_HF_HADRON_ENERGY;
   std::vector<float>   *RECO_PFJET_HF_EM_ENERGY;
   std::vector<float>   *RECO_PFJET_CHARGED_EM_ENERGY;
   std::vector<float>   *RECO_PFJET_CHARGED_MU_ENERGY;
   std::vector<float>   *RECO_PFJET_NEUTRAL_EM_ENERGY;
   std::vector<int>     *RECO_PFJET_CHARGED_HADRON_MULTIPLICITY;
   std::vector<int>     *RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY;
   std::vector<int>     *RECO_PFJET_PHOTON_MULTIPLICITY;
   std::vector<int>     *RECO_PFJET_ELECTRON_MULTIPLICITY;
   std::vector<int>     *RECO_PFJET_MUON_MULTIPLICITY;
   std::vector<int>     *RECO_PFJET_HF_HADRON_MULTIPLICTY;
   std::vector<int>     *RECO_PFJET_HF_EM_MULTIPLICITY;
   std::vector<int>     *RECO_PFJET_CHARGED_MULTIPLICITY;
   std::vector<int>     *RECO_PFJET_NEUTRAL_MULTIPLICITY;
   std::vector<int>     *RECO_PFJET_NCOMPONENTS;
   std::vector<std::vector<int> > *RECO_PFJET_COMPONENT_PDGID;
   std::vector<std::vector<float> > *RECO_PFJET_COMPONENT_PT;
   std::vector<std::vector<float> > *RECO_PFJET_COMPONENT_ETA;
   std::vector<std::vector<float> > *RECO_PFJET_COMPONENT_PHI;
   std::vector<std::vector<float> > *RECO_PFJET_COMPONENT_E;
   std::vector<std::vector<float> > *RECO_PFJET_COMPONENT_CHARGE;
   std::vector<std::vector<float> > *RECO_PFJET_COMPONENT_TRANSVERSE_MASS;
   std::vector<std::vector<float> > *RECO_PFJET_COMPONENT_XVERTEX;
   std::vector<std::vector<float> > *RECO_PFJET_COMPONENT_YVERTEX;
   std::vector<std::vector<float> > *RECO_PFJET_COMPONENT_ZVERTEX;
   std::vector<std::vector<float> > *RECO_PFJET_COMPONENT_VERTEX_CHI2;
   Float_t         RECO_CALOMET;
   Float_t         RECO_PFMET;
   Float_t         RECO_PFMET_X;
   Float_t         RECO_PFMET_Y;
   Float_t         RECO_PFMET_PHI;
   Float_t         RECO_PFMET_THETA;
   Float_t         RECO_PFMET_uncorr;
   Float_t         RECO_PFMET_X_uncorr;
   Float_t         RECO_PFMET_Y_uncorr;
   Float_t         RECO_PFMET_PHI_uncorr;
   Float_t         RECO_PFMET_THETA_uncorr;
   Float_t         RECO_PFMET_JetEnUp;
   Float_t         RECO_PFMET_JetEnDn;
   Float_t         RECO_PFMET_ElectronEnUp;
   Float_t         RECO_PFMET_ElectronEnDn;
   Float_t         RECO_PFMET_MuonEnUp;
   Float_t         RECO_PFMET_MuonEnDn;
   Float_t         RECO_PFMET_JetResUp;
   Float_t         RECO_PFMET_JetResDn;
   Float_t         RECO_PFMET_UnclusteredEnUp;
   Float_t         RECO_PFMET_UnclusteredEnDn;
   Float_t         RECO_PFMET_PhotonEnUp;
   Float_t         RECO_PFMET_PhotonEnDn;
   Float_t         RECO_PFMET_TauEnUp ;
   Float_t         RECO_PFMET_TauEnDown;
   Int_t           RECO_PFMET_GoodVtxNoiseFilter;
   Int_t           RECO_PFMET_GlobalSuperTightHalo2016NoiseFilter;
   Int_t           RECO_PFMET_HBHENoiseFilter;
   Int_t           RECO_PFMET_HBHENoiseIsoFilter;
   Int_t           RECO_PFMET_EcalDeadCellTriggerPrimitiveNoiseFilter;
   Int_t           RECO_PFMET_BadPFMuonFilter;
   Int_t           RECO_PFMET_BadChargedCandidateFilter;
   Int_t           RECO_PFMET_EEBadScNoiseFilter;
   Int_t           RECO_PFMET_EcalBadCalibFilter;
   Float_t         RECO_TCMET;
   Float_t         RECO_CORMETMUONS;
   std::vector<float>   *tCHighEff_BTagJet_PT;
   std::vector<float>   *tCHighPur_BTagJet_PT;
   std::vector<float>   *cSV_BTagJet_PT;
   std::vector<float>   *tCHighEff_BTagJet_ETA;
   std::vector<float>   *tCHighPur_BTagJet_ETA;
   std::vector<float>   *cSV_BTagJet_ETA;
   std::vector<float>   *tCHighEff_BTagJet_PHI;
   std::vector<float>   *tCHighPur_BTagJet_PHI;
   std::vector<float>   *cSV_BTagJet_PHI;
   std::vector<float>   *cSV_BTagJet_ET;
   std::vector<float>   *tCHighEff_BTagJet_DISCR;
   std::vector<float>   *tCHighPur_BTagJet_DISCR;
   std::vector<float>   *cSV_BTagJet_DISCR;

   // List of branches
   TBranch        *b_irun;   //!
   TBranch        *b_ievt;   //!
   TBranch        *b_ils;   //!
   TBranch        *b_Avginstlumi;   //!
   TBranch        *b_num_PU_vertices;   //!
   TBranch        *b_PU_BunchCrossing;   //!
   TBranch        *b_MC_weighting;   //!
   TBranch        *b_RECO_nMuHLTMatch;   //!
   TBranch        *b_RECOMU_PT_MuHLTMatch;   //!
   TBranch        *b_RECOMU_ETA_MuHLTMatch;   //!
   TBranch        *b_RECOMU_PHI_MuHLTMatch;   //!
   TBranch        *b_RECO_nEleHLTMatch;   //!
   TBranch        *b_RECOELE_PT_EleHLTMatch;   //!
   TBranch        *b_RECOELE_ETA_EleHLTMatch;   //!
   TBranch        *b_RECOELE_PHI_EleHLTMatch;   //!
   TBranch        *b_HLTPathsFired;   //!
   TBranch        *b_MC_E;   //!
   TBranch        *b_MC_PT;   //!
   TBranch        *b_MC_ETA;   //!
   TBranch        *b_MC_THETA;   //!
   TBranch        *b_MC_PHI;   //!
   TBranch        *b_MC_MASS;   //!
   TBranch        *b_MC_PDGID;   //!
   TBranch        *b_MC_LEPT_PT;   //!
   TBranch        *b_MC_LEPT_ETA;   //!
   TBranch        *b_MC_LEPT_PHI;   //!
   TBranch        *b_MC_LEPT_THETA;   //!
   TBranch        *b_MC_LEPT_PDGID;   //!
   TBranch        *b_MC_Z_PT;   //!
   TBranch        *b_MC_Z_ETA;   //!
   TBranch        *b_MC_Z_PHI;   //!
   TBranch        *b_MC_Z_THETA;   //!
   TBranch        *b_MC_Z_MASS;   //!
   TBranch        *b_MC_Z_PDGID;   //!
   TBranch        *b_MC_fourl_MASS;   //!
   TBranch        *b_MC_fourl_PT;   //!
   TBranch        *b_MC_fourl_PDGID;   //!
   TBranch        *b_MC_ZZ_MASS;   //!
   TBranch        *b_MC_ZZ_PT;   //!
   TBranch        *b_MC_ZZ_ETA;   //!
   TBranch        *b_MC_ZZ_PHI;   //!
   TBranch        *b_MC_ZZ_THETA;   //!
   TBranch        *b_MC_ZZ_PDGID;   //!
   TBranch        *b_MC_ZZ_D_MASS;   //!
   TBranch        *b_MC_ZZ_D_PT;   //!
   TBranch        *b_MC_ZZ_D_ETA;   //!
   TBranch        *b_MC_ZZ_D_PHI;   //!
   TBranch        *b_MC_ZZ_D_THETA;   //!
   TBranch        *b_MC_ZZ_D_PDGID;   //!
   TBranch        *b_MC_ZZ_DD_MASS;   //!
   TBranch        *b_MC_ZZ_DD_PT;   //!
   TBranch        *b_MC_ZZ_DD_ETA;   //!
   TBranch        *b_MC_ZZ_DD_PHI;   //!
   TBranch        *b_MC_ZZ_DD_THETA;   //!
   TBranch        *b_MC_ZZ_DD_PDGID;   //!
   TBranch        *b_MC_GENJET_PT;   //!
   TBranch        *b_MC_GENJET_ETA;   //!
   TBranch        *b_MC_GENJET_PHI;   //!
   TBranch        *b_MC_GENMET;   //!
   TBranch        *b_RECORF_2e2mu_cosTheta1_spin;   //!
   TBranch        *b_RECORF_2e2mu_cosTheta2_spin;   //!
   TBranch        *b_RECORF_2e2mu_cosThetaStar_spin;   //!
   TBranch        *b_RECORF_2e2mu_Phi_spin;   //!
   TBranch        *b_RECORF_2e2mu_Phi1_spin;   //!
   TBranch        *b_RECORF_2e2mu_Phi2_spin;   //!
   TBranch        *b_RECORF_2e2mu_phi1RF_spin;   //!
   TBranch        *b_RECORF_2e2mu_phi2RF_spin;   //!
   TBranch        *b_RECORF_2e2mu_MELA;   //!
   TBranch        *b_RECORF_4e_cosTheta1_spin;   //!
   TBranch        *b_RECORF_4e_cosTheta2_spin;   //!
   TBranch        *b_RECORF_4e_cosThetaStar_spin;   //!
   TBranch        *b_RECORF_4e_Phi_spin;   //!
   TBranch        *b_RECORF_4e_Phi1_spin;   //!
   TBranch        *b_RECORF_4e_Phi2_spin;   //!
   TBranch        *b_RECORF_4e_phi1RF_spin;   //!
   TBranch        *b_RECORF_4e_phi2RF_spin;   //!
   TBranch        *b_RECORF_4e_MELA;   //!
   TBranch        *b_RECORF_4mu_cosTheta1_spin;   //!
   TBranch        *b_RECORF_4mu_cosTheta2_spin;   //!
   TBranch        *b_RECORF_4mu_cosThetaStar_spin;   //!
   TBranch        *b_RECORF_4mu_Phi_spin;   //!
   TBranch        *b_RECORF_4mu_Phi1_spin;   //!
   TBranch        *b_RECORF_4mu_Phi2_spin;   //!
   TBranch        *b_RECORF_4mu_phi1RF_spin;   //!
   TBranch        *b_RECORF_4mu_phi2RF_spin;   //!
   TBranch        *b_RECORF_4mu_MELA;   //!
   TBranch        *b_RECO_ZMM_MASS;   //!
   TBranch        *b_RECO_ZEE_MASS;   //!
   TBranch        *b_RECO_DiLep_MASS;   //!
   TBranch        *b_RECO_ZMM_PT;   //!
   TBranch        *b_RECO_ZEE_PT;   //!
   TBranch        *b_RECO_DiLep_PT;   //!
   TBranch        *b_RECO_ZMM_ETA;   //!
   TBranch        *b_RECO_ZEE_ETA;   //!
   TBranch        *b_RECO_DiLep_ETA;   //!
   TBranch        *b_RECO_ZMM_PHI;   //!
   TBranch        *b_RECO_ZEE_PHI;   //!
   TBranch        *b_RECO_DiLep_PHI;   //!
   TBranch        *b_RECO_ZMMss_MASS;   //!
   TBranch        *b_RECO_ZEEss_MASS;   //!
   TBranch        *b_RECO_ZEM_MASS;   //!
   TBranch        *b_RECO_ZMMss_PT;   //!
   TBranch        *b_RECO_ZEEss_PT;   //!
   TBranch        *b_RECO_ZEM_PT;   //!
   TBranch        *b_RECO_ZMMss_ETA;   //!
   TBranch        *b_RECO_ZEEss_ETA;   //!
   TBranch        *b_RECO_ZEM_ETA;   //!
   TBranch        *b_RECO_ZMMss_PHI;   //!
   TBranch        *b_RECO_ZEEss_PHI;   //!
   TBranch        *b_RECO_ZEM_PHI;   //!
   TBranch        *b_RECO_MMMM_MASS;   //!
   TBranch        *b_RECO_MMMM_PT;   //!
   TBranch        *b_RECO_MMMM_ETA;   //!
   TBranch        *b_RECO_MMMM_PHI;   //!
   TBranch        *b_RECO_MMMM_MASS_REFIT;   //!
   TBranch        *b_RECO_EEEE_MASS;   //!
   TBranch        *b_RECO_EEEE_PT;   //!
   TBranch        *b_RECO_EEEE_ETA;   //!
   TBranch        *b_RECO_EEEE_PHI;   //!
   TBranch        *b_RECO_EEEE_MASS_REFIT;   //!
   TBranch        *b_RECO_EEMM_MASS;   //!
   TBranch        *b_RECO_EEMM_PT;   //!
   TBranch        *b_RECO_EEMM_ETA;   //!
   TBranch        *b_RECO_EEMM_PHI;   //!
   TBranch        *b_RECO_EEMM_MASS_REFIT;   //!
   TBranch        *b_RECO_LLL0_MASS;   //!
   TBranch        *b_RECO_LLL1_MASS;   //!
   TBranch        *b_RECO_LLL2_MASS;   //!
   TBranch        *b_RECO_LLL3_MASS;   //!
   TBranch        *b_RECO_LLL0_PT;   //!
   TBranch        *b_RECO_LLL1_PT;   //!
   TBranch        *b_RECO_LLL2_PT;   //!
   TBranch        *b_RECO_LLL3_PT;   //!
   TBranch        *b_RECO_LLLl0_MASS;   //!
   TBranch        *b_RECO_LLLl1_MASS;   //!
   TBranch        *b_RECO_LLLl0_PT;   //!
   TBranch        *b_RECO_LLLl1_PT;   //!
   TBranch        *b_RECO_LLLL0ss_MASS;   //!
   TBranch        *b_RECO_LLLL0ss_PT;   //!
   TBranch        *b_RECO_LLLL1ss_MASS;   //!
   TBranch        *b_RECO_LLLL1ss_PT;   //!
   TBranch        *b_RECO_LLLL2ss_MASS;   //!
   TBranch        *b_RECO_LLLL2ss_PT;   //!
   TBranch        *b_RECO_LLLL_MASS;   //!
   TBranch        *b_RECO_LLLL_PT;   //!
   TBranch        *b_RECO_LLLL_ETA;   //!
   TBranch        *b_RECO_LLLL_PHI;   //!
   TBranch        *b_RECOELE_E;   //!
   TBranch        *b_RECOELE_PT;   //!
   TBranch        *b_RECOELE_PTError;   //!
   TBranch        *b_RECOELE_P;   //!
   TBranch        *b_RECOELE_ETA;   //!
   TBranch        *b_RECOELE_THETA;   //!
   TBranch        *b_RECOELE_PHI;   //!
   TBranch        *b_RECOELE_MASS;   //!
   TBranch        *b_RECOELE_CHARGE;   //!
   TBranch        *b_RECOELE_ID;   //!
   TBranch        *b_RECOELE_PT_uncorr;   //!
   TBranch        *b_RECOELE_isEcalDriven;   //!
   TBranch        *b_RECOELE_isTrackerDriven;   //!
   TBranch        *b_RECOELE_gsftrack_NPixHits;   //!
   TBranch        *b_RECOELE_gsftrack_NStripHits;   //!
   TBranch        *b_RECOELE_gsftrack_chi2;   //!
   TBranch        *b_RECOELE_gsftrack_dxyB;   //!
   TBranch        *b_RECOELE_gsftrack_dxy;   //!
   TBranch        *b_RECOELE_gsftrack_dxyError;   //!
   TBranch        *b_RECOELE_gsftrack_dzB;   //!
   TBranch        *b_RECOELE_gsftrack_dz;   //!
   TBranch        *b_RECOELE_gsftrack_dzError;   //!
   TBranch        *b_RECOELE_gsftrack_losthits;   //!
   TBranch        *b_RECOELE_gsftrack_validhits;   //!
   TBranch        *b_RECOELE_gsftrack_expected_inner_hits;   //!
   TBranch        *b_RECOELE_scl_E;   //!
   TBranch        *b_RECOELE_scl_Et;   //!
   TBranch        *b_RECOELE_scl_Eta;   //!
   TBranch        *b_RECOELE_scl_Phi;   //!
   TBranch        *b_RECOELE_ecalEnergy;   //!
   TBranch        *b_RECOELE_ep;   //!
   TBranch        *b_RECOELE_eSeedp;   //!
   TBranch        *b_RECOELE_eSeedpout;   //!
   TBranch        *b_RECOELE_eElepout;   //!
   TBranch        *b_RECOELE_deltaEtaIn;   //!
   TBranch        *b_RECOELE_deltaEtaSeed;   //!
   TBranch        *b_RECOELE_deltaEtaEle;   //!
   TBranch        *b_RECOELE_deltaPhiIn;   //!
   TBranch        *b_RECOELE_deltaPhiSeed;   //!
   TBranch        *b_RECOELE_deltaPhiEle;   //!
   TBranch        *b_RECOELE_isbarrel;   //!
   TBranch        *b_RECOELE_isendcap;   //!
   TBranch        *b_RECOELE_isGap;   //!
   TBranch        *b_RECOELE_isEBetaGap;   //!
   TBranch        *b_RECOELE_isEBphiGap;   //!
   TBranch        *b_RECOELE_isEEdeeGap;   //!
   TBranch        *b_RECOELE_isEEringGap;   //!
   TBranch        *b_RECOELE_sigmaIetaIeta;   //!
   TBranch        *b_RECOELE_sigmaEtaEta;   //!
   TBranch        *b_RECOELE_e15;   //!
   TBranch        *b_RECOELE_e25max;   //!
   TBranch        *b_RECOELE_e55;   //!
   TBranch        *b_RECOELE_he;   //!
   TBranch        *b_RECOELE_r9;   //!
   TBranch        *b_RECOELE_mva;   //!
   TBranch        *b_RECOELE_fbrem;   //!
   TBranch        *b_RECOELE_nbrems;   //!
   TBranch        *b_RECOELE_Class;   //!
   TBranch        *b_RECOELE_fbrem_mode;   //!
   TBranch        *b_RECOELE_fbrem_mean;   //!
   TBranch        *b_RECOELE_EGMTRACKISO;   //!
   TBranch        *b_RECOELE_EGMHCALISO;   //!
   TBranch        *b_RECOELE_EGMECALISO;   //!
   TBranch        *b_RECOELE_EGMX;   //!
   TBranch        *b_RECOELE_PFchAllPart;   //!
   TBranch        *b_RECOELE_PFchHad;   //!
   TBranch        *b_RECOELE_PFneuHad;   //!
   TBranch        *b_RECOELE_PFphoton;   //!
   TBranch        *b_RECOELE_PFPUchAllPart;   //!
   TBranch        *b_RECOELE_PFX_dB;   //!
   TBranch        *b_RECOELE_PFX_rho;   //!
   TBranch        *b_RECOELE_regEnergy;   //!
   TBranch        *b_RECOELE_regEnergyError;   //!
   TBranch        *b_RECOELE_SIP;   //!
   TBranch        *b_RECOELE_IP;   //!
   TBranch        *b_RECOELE_IPERROR;   //!
   TBranch        *b_RECOELE_SIP_KF;   //!
   TBranch        *b_RECOELE_IP_KF;   //!
   TBranch        *b_RECOELE_IPERROR_KF;   //!
   TBranch        *b_RECOELE_SIP_GD;   //!
   TBranch        *b_RECOELE_SIP_GDEEEE;   //!
   TBranch        *b_RECOELE_SIP_Std;   //!
   TBranch        *b_RECOELE_SIP_StdEEEE;   //!
   TBranch        *b_RECOELE_SIP_Kin;   //!
   TBranch        *b_RECOELE_SIP_KinEEEE;   //!
   TBranch        *b_RECOELE_STIP;   //!
   TBranch        *b_RECOELE_SLIP;   //!
   TBranch        *b_RECOELE_TIP;   //!
   TBranch        *b_RECOELE_LIP;   //!
   TBranch        *b_RECOELE_TIPERROR;   //!
   TBranch        *b_RECOELE_LIPERROR;   //!
   TBranch        *b_RECOELE_sclRawE;   //!
   TBranch        *b_RECOELE_sclX;   //!
   TBranch        *b_RECOELE_sclY;   //!
   TBranch        *b_RECOELE_sclZ;   //!
   TBranch        *b_RECOELE_seedSubdet1;   //!
   TBranch        *b_RECOELE_seedDphi1;   //!
   TBranch        *b_RECOELE_seedDrz1;   //!
   TBranch        *b_RECOELE_seedSubdet2;   //!
   TBranch        *b_RECOELE_seedDphi2;   //!
   TBranch        *b_RECOELE_seedDrz2;   //!
   TBranch        *b_RECOELE_eidVeryLoose;   //!
   TBranch        *b_RECOELE_eidLoose;   //!
   TBranch        *b_RECOELE_eidMedium;   //!
   TBranch        *b_RECOELE_eidTight;   //!
   TBranch        *b_RECOELE_eidHZZVeryLoose;   //!
   TBranch        *b_RECOELE_eidHZZLoose;   //!
   TBranch        *b_RECOELE_eidHZZMedium;   //!
   TBranch        *b_RECOELE_eidHZZTight;   //!
   TBranch        *b_RECOELE_mvaTrigV0;   //!
   TBranch        *b_RECOELE_mvaNonTrigV0;   //!
   TBranch        *b_RECOELE_COV;   //!
   TBranch        *b_RECOELE_TLE_ParentSC_X;   //!
   TBranch        *b_RECOELE_TLE_ParentSC_Y;   //!
   TBranch        *b_RECOELE_TLE_ParentSC_Z;   //!
   TBranch        *b_RECOELE_ecalTrkEnergyPreCorr;   //!
   TBranch        *b_RECOELE_ecalTrkEnergyErrPreCorr;   //!
   TBranch        *b_RECOELE_ecalTrkEnergyErrPostCorr;   //!
   TBranch        *b_RECOELE_energyScaleValue;   //!
   TBranch        *b_RECOELE_energySigmaValue;   //!
   TBranch        *b_RECOELE_energyScaleUp;   //!
   TBranch        *b_RECOELE_energyScaleDown;   //!
   TBranch        *b_RECOELE_energyScaleStatUp;   //!
   TBranch        *b_RECOELE_energyScaleStatDown;   //!
   TBranch        *b_RECOELE_energyScaleSystUp;   //!
   TBranch        *b_RECOELE_energyScaleSystDown;   //!
   TBranch        *b_RECOELE_energyScaleGainUp;   //!
   TBranch        *b_RECOELE_energyScaleGainDown;   //!
   TBranch        *b_RECOELE_energyScaleEtUp;   //!
   TBranch        *b_RECOELE_energyScaleEtDown;   //!
   TBranch        *b_RECOELE_energySigmaUp;   //!
   TBranch        *b_RECOELE_energySigmaDown;   //!
   TBranch        *b_RECOELE_energySigmaPhiUp;   //!
   TBranch        *b_RECOELE_energySigmaPhiDown;   //!
   TBranch        *b_RECOELE_energySigmaRhoUp;   //!
   TBranch        *b_RECOELE_energySigmaRhoDown;   //!
   TBranch        *b_RECOMU_isPFMu;   //!
   TBranch        *b_RECOMU_isGlobalMu;   //!
   TBranch        *b_RECOMU_isStandAloneMu;   //!
   TBranch        *b_RECOMU_isTrackerMu;   //!
   TBranch        *b_RECOMU_isCaloMu;   //!
   TBranch        *b_RECOMU_isTrackerHighPtMu;   //!
   TBranch        *b_RECOMU_isME0Muon;   //!
   TBranch        *b_RECOMU_E;   //!
   TBranch        *b_RECOMU_PT;   //!
   TBranch        *b_RECOMU_P;   //!
   TBranch        *b_RECOMU_ETA;   //!
   TBranch        *b_RECOMU_THETA;   //!
   TBranch        *b_RECOMU_PHI;   //!
   TBranch        *b_RECOMU_MASS;   //!
   TBranch        *b_RECOMU_CHARGE;   //!
   TBranch        *b_RECOMU_PT_uncorr;   //!
   TBranch        *b_RECOMU_COV;   //!
   TBranch        *b_RECOMU_TRACKISO;   //!
   TBranch        *b_RECOMU_TRACKISO_SUMPT;   //!
   TBranch        *b_RECOMU_HCALISO;   //!
   TBranch        *b_RECOMU_ECALISO;   //!
   TBranch        *b_RECOMU_X;   //!
   TBranch        *b_RECOMU_PFchHad;   //!
   TBranch        *b_RECOMU_PFneuHad;   //!
   TBranch        *b_RECOMU_PFphoton;   //!
   TBranch        *b_RECOMU_PFPUchAllPart;   //!
   TBranch        *b_RECOMU_PFX_dB;   //!
   TBranch        *b_RECOMU_PFX_rho;   //!
   TBranch        *b_RECOPFPHOT_PFchHad;   //!
   TBranch        *b_RECOPFPHOT_PFneuHad;   //!
   TBranch        *b_RECOPFPHOT_PFphoton;   //!
   TBranch        *b_RECOPFPHOT_PFPUchAllPart;   //!
   TBranch        *b_RECOPFPHOT_PFX_rho;   //!
   TBranch        *b_RECOMU_SIP;   //!
   TBranch        *b_RECOMU_IP;   //!
   TBranch        *b_RECOMU_IPERROR;   //!
   TBranch        *b_RECOMU_SIP_KF;   //!
   TBranch        *b_RECOMU_IP_KF;   //!
   TBranch        *b_RECOMU_IPERROR_KF;   //!
   TBranch        *b_RECOMU_SIP_GD;   //!
   TBranch        *b_RECOMU_SIP_GDMMMM;   //!
   TBranch        *b_RECOMU_SIP_Std;   //!
   TBranch        *b_RECOMU_SIP_StdMMMM;   //!
   TBranch        *b_RECOMU_SIP_Kin;   //!
   TBranch        *b_RECOMU_SIP_KinMMMM;   //!
   TBranch        *b_RECOMU_STIP;   //!
   TBranch        *b_RECOMU_SLIP;   //!
   TBranch        *b_RECOMU_TIP;   //!
   TBranch        *b_RECOMU_LIP;   //!
   TBranch        *b_RECOMU_TIPERROR;   //!
   TBranch        *b_RECOMU_LIPERROR;   //!
   TBranch        *b_RECOMU_caloCompatibility;   //!
   TBranch        *b_RECOMU_segmentCompatibility;   //!
   TBranch        *b_RECOMU_numberOfMatches;   //!
   TBranch        *b_RECOMU_numberOfMatchedStations;   //!
   TBranch        *b_RECOMU_glbmuPromptTight;   //!
   TBranch        *b_RECOMU_trkmuArbitration;   //!
   TBranch        *b_RECOMU_trkmu2DCompatibilityLoose;   //!
   TBranch        *b_RECOMU_trkmu2DCompatibilityTight;   //!
   TBranch        *b_RECOMU_trkmuOneStationLoose;   //!
   TBranch        *b_RECOMU_trkmuOneStationTight;   //!
   TBranch        *b_RECOMU_trkmuLastStationLoose;   //!
   TBranch        *b_RECOMU_trkmuLastStationTight;   //!
   TBranch        *b_RECOMU_trkmuOneStationAngLoose;   //!
   TBranch        *b_RECOMU_trkmuOneStationAngTight;   //!
   TBranch        *b_RECOMU_trkmuLastStationAngLoose;   //!
   TBranch        *b_RECOMU_trkmuLastStationAngTight;   //!
   TBranch        *b_RECOMU_trkmuLastStationOptimizedLowPtLoose;   //!
   TBranch        *b_RECOMU_trkmuLastStationOptimizedLowPtTight;   //!
   TBranch        *b_RECOMU_mutrkPT;   //!
   TBranch        *b_RECOMU_mutrkPTError;   //!
   TBranch        *b_RECOMU_mutrkDxy;   //!
   TBranch        *b_RECOMU_mutrkDxyError;   //!
   TBranch        *b_RECOMU_mutrkDxyB;   //!
   TBranch        *b_RECOMU_mutrkDz;   //!
   TBranch        *b_RECOMU_mutrkDzError;   //!
   TBranch        *b_RECOMU_mutrkDzB;   //!
   TBranch        *b_RECOMU_mutrkChi2PerNdof;   //!
   TBranch        *b_RECOMU_mutrkCharge;   //!
   TBranch        *b_RECOMU_mutrkNHits;   //!
   TBranch        *b_RECOMU_mutrkNStripHits;   //!
   TBranch        *b_RECOMU_mutrkNPixHits;   //!
   TBranch        *b_RECOMU_mutrkNMuonHits;   //!
   TBranch        *b_RECOMU_mutrktrackerLayersWithMeasurement;   //!
   TBranch        *b_RECOMU_muInnertrkDxy;   //!
   TBranch        *b_RECOMU_muInnertrkDxyError;   //!
   TBranch        *b_RECOMU_muInnertrkDxyB;   //!
   TBranch        *b_RECOMU_muInnertrkDz;   //!
   TBranch        *b_RECOMU_muInnertrkDzError;   //!
   TBranch        *b_RECOMU_muInnertrkDzB;   //!
   TBranch        *b_RECOMU_muInnertrkChi2PerNdof;   //!
   TBranch        *b_RECOMU_muInnertrktrackerLayersWithMeasurement;   //!
   TBranch        *b_RECOMU_muInnertrkPT;   //!
   TBranch        *b_RECOMU_muInnertrkPTError;   //!
   TBranch        *b_RECOMU_muInnertrkCharge;   //!
   TBranch        *b_RECOMU_muInnertrkNHits;   //!
   TBranch        *b_RECOMU_muInnertrkNStripHits;   //!
   TBranch        *b_RECOMU_muInnertrkNPixHits;   //!
   TBranch        *b_RECOMU_mubesttrkType;   //!
   TBranch        *b_RECOMU_mubesttrkDxy;   //!
   TBranch        *b_RECOMU_mubesttrkDxyError;   //!
   TBranch        *b_RECOMU_mubesttrkDz;   //!
   TBranch        *b_RECOMU_mubesttrkDzError;   //!
   TBranch        *b_RECOMU_mubesttrkPTError;   //!
   TBranch        *b_RECOMU_Rochester_Error;   //!
   TBranch        *b_ftsigma;   //!
   TBranch        *b_gdX;   //!
   TBranch        *b_gdY;   //!
   TBranch        *b_gdZ;   //!
   TBranch        *b_ftsigmalag;   //!
   TBranch        *b_gdlagX;   //!
   TBranch        *b_gdlagY;   //!
   TBranch        *b_gdlagZ;   //!
   TBranch        *b_gdlagProb;   //!
   TBranch        *b_gdlagNdof;   //!
   TBranch        *b_ftsigmaMMMM;   //!
   TBranch        *b_gdXMMMM;   //!
   TBranch        *b_gdYMMMM;   //!
   TBranch        *b_gdZMMMM;   //!
   TBranch        *b_ftsigmalagMMMM;   //!
   TBranch        *b_gdlagXMMMM;   //!
   TBranch        *b_gdlagYMMMM;   //!
   TBranch        *b_gdlagZMMMM;   //!
   TBranch        *b_gdlagProbMMMM;   //!
   TBranch        *b_gdlagNdofMMMM;   //!
   TBranch        *b_ftsigmaEEEE;   //!
   TBranch        *b_gdXEEEE;   //!
   TBranch        *b_gdYEEEE;   //!
   TBranch        *b_gdZEEEE;   //!
   TBranch        *b_ftsigmalagEEEE;   //!
   TBranch        *b_gdlagXEEEE;   //!
   TBranch        *b_gdlagYEEEE;   //!
   TBranch        *b_gdlagZEEEE;   //!
   TBranch        *b_gdlagProbEEEE;   //!
   TBranch        *b_gdlagNdofEEEE;   //!
   TBranch        *b_StdFitVertexX;   //!
   TBranch        *b_StdFitVertexY;   //!
   TBranch        *b_StdFitVertexZ;   //!
   TBranch        *b_StdFitVertexChi2r;   //!
   TBranch        *b_StdFitVertexProb;   //!
   TBranch        *b_StdFitVertexTrack_PT;   //!
   TBranch        *b_StdFitVertexTrack_ETA;   //!
   TBranch        *b_StdFitVertexTrack_PHI;   //!
   TBranch        *b_KinFitVertexX;   //!
   TBranch        *b_KinFitVertexY;   //!
   TBranch        *b_KinFitVertexZ;   //!
   TBranch        *b_KinFitVertexChi2r;   //!
   TBranch        *b_KinFitVertexProb;   //!
   TBranch        *b_StdFitVertexXMMMM;   //!
   TBranch        *b_StdFitVertexYMMMM;   //!
   TBranch        *b_StdFitVertexZMMMM;   //!
   TBranch        *b_StdFitVertexChi2rMMMM;   //!
   TBranch        *b_StdFitVertexProbMMMM;   //!
   TBranch        *b_StdFitVertexTrackMMMM_PT;   //!
   TBranch        *b_StdFitVertexTrackMMMM_ETA;   //!
   TBranch        *b_StdFitVertexTrackMMMM_PHI;   //!
   TBranch        *b_KinFitVertexXMMMM;   //!
   TBranch        *b_KinFitVertexYMMMM;   //!
   TBranch        *b_KinFitVertexZMMMM;   //!
   TBranch        *b_KinFitVertexChi2rMMMM;   //!
   TBranch        *b_KinFitVertexProbMMMM;   //!
   TBranch        *b_StdFitVertexXEEEE;   //!
   TBranch        *b_StdFitVertexYEEEE;   //!
   TBranch        *b_StdFitVertexZEEEE;   //!
   TBranch        *b_StdFitVertexChi2rEEEE;   //!
   TBranch        *b_StdFitVertexProbEEEE;   //!
   TBranch        *b_StdFitVertexTrackEEEE_PT;   //!
   TBranch        *b_StdFitVertexTrackEEEE_ETA;   //!
   TBranch        *b_StdFitVertexTrackEEEE_PHI;   //!
   TBranch        *b_KinFitVertexXEEEE;   //!
   TBranch        *b_KinFitVertexYEEEE;   //!
   TBranch        *b_KinFitVertexZEEEE;   //!
   TBranch        *b_KinFitVertexChi2rEEEE;   //!
   TBranch        *b_KinFitVertexProbEEEE;   //!
   TBranch        *b_StdFitVertexChi2rMMM;   //!
   TBranch        *b_StdFitVertexProbMMM;   //!
   TBranch        *b_StdFitVertexChi2rMME;   //!
   TBranch        *b_StdFitVertexProbMME;   //!
   TBranch        *b_StdFitVertexChi2rEEE;   //!
   TBranch        *b_StdFitVertexProbEEE;   //!
   TBranch        *b_StdFitVertexChi2rMEE;   //!
   TBranch        *b_StdFitVertexProbMEE;   //!
   TBranch        *b_StdFitVertexChi2rDiLep;   //!
   TBranch        *b_StdFitVertexProbDiLep;   //!
   TBranch        *b_ConvMapDist;   //!
   TBranch        *b_ConvMapDcot;   //!
   TBranch        *b_RECOMU_MatchingMCTruth;   //!
   TBranch        *b_RECOMU_MatchingMCpT;   //!
   TBranch        *b_RECOMU_MatchingMCEta;   //!
   TBranch        *b_RECOMU_MatchingMCPhi;   //!
   TBranch        *b_RECOELE_MatchingMCTruth;   //!
   TBranch        *b_RECOELE_MatchingMCpT;   //!
   TBranch        *b_RECOELE_MatchingMCEta;   //!
   TBranch        *b_RECOELE_MatchingMCPhi;   //!
   TBranch        *b_RECOPHOT_MatchingMCTruth;   //!
   TBranch        *b_RECOPHOT_MatchingMCpT;   //!
   TBranch        *b_RECOPHOT_MatchingMCEta;   //!
   TBranch        *b_RECOPHOT_MatchingMCPhi;   //!
   TBranch        *b_RECOzMuMu_MatchingMCTruth;   //!
   TBranch        *b_RECOzMuMu_MatchingMCpT;   //!
   TBranch        *b_RECOzMuMu_MatchingMCmass;   //!
   TBranch        *b_RECOzMuMu_MatchingMCEta;   //!
   TBranch        *b_RECOzMuMu_MatchingMCPhi;   //!
   TBranch        *b_RECOzEE_MatchingMCTruth;   //!
   TBranch        *b_RECOzEE_MatchingMCpT;   //!
   TBranch        *b_RECOzEE_MatchingMCmass;   //!
   TBranch        *b_RECOzEE_MatchingMCEta;   //!
   TBranch        *b_RECOzEE_MatchingMCPhi;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCTruth;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCpT;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCmass;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCEta;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCPhi;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCTruth;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCpT;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCmass;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCEta;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCPhi;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCTruth;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCpT;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCmass;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCEta;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCPhi;   //!
   TBranch        *b_RECO_NMU;   //!
   TBranch        *b_RECO_NELE;   //!
   TBranch        *b_RECO_NTRACK;   //!
   TBranch        *b_RECO_NPHOT;   //!
   TBranch        *b_RECOPHOT_PT;   //!
   TBranch        *b_RECOPHOT_ETA;   //!
   TBranch        *b_RECOPHOT_PHI;   //!
   TBranch        *b_RECOPHOT_THETA;   //!
   TBranch        *b_RECOPHOT_TLE_ParentSC_X;   //!
   TBranch        *b_RECOPHOT_TLE_ParentSC_Y;   //!
   TBranch        *b_RECOPHOT_TLE_ParentSC_Z;   //!
   TBranch        *b_RECO_NPFPHOT;   //!
   TBranch        *b_RECOPFPHOT_PT;   //!
   TBranch        *b_RECOPFPHOT_PTError;   //!
   TBranch        *b_RECOPFPHOT_ETA;   //!
   TBranch        *b_RECOPFPHOT_PHI;   //!
   TBranch        *b_RECOPFPHOT_THETA;   //!
   TBranch        *b_RECOPFPHOT_PT_uncorr;   //!
   TBranch        *b_RECOPFPHOT_ecalEnergyPreCorr;   //!
   TBranch        *b_RECOPFPHOT_ecalEnergyErrPreCorr;   //!
   TBranch        *b_RECOPFPHOT_ecalEnergyErrPostCorr;   //!
   TBranch        *b_RECOPFPHOT_energyScaleValue;   //!
   TBranch        *b_RECOPFPHOT_energySigmaValue;   //!
   TBranch        *b_RECOPFPHOT_energyScaleUp;   //!
   TBranch        *b_RECOPFPHOT_energyScaleDown;   //!
   TBranch        *b_RECOPFPHOT_energyScaleStatUp;   //!
   TBranch        *b_RECOPFPHOT_energyScaleStatDown;   //!
   TBranch        *b_RECOPFPHOT_energyScaleSystUp;   //!
   TBranch        *b_RECOPFPHOT_energyScaleSystDown;   //!
   TBranch        *b_RECOPFPHOT_energyScaleGainUp;   //!
   TBranch        *b_RECOPFPHOT_energyScaleGainDown;   //!
   TBranch        *b_RECOPFPHOT_energyScaleEtUp;   //!
   TBranch        *b_RECOPFPHOT_energyScaleEtDown;   //!
   TBranch        *b_RECOPFPHOT_energySigmaUp;   //!
   TBranch        *b_RECOPFPHOT_energySigmaDown;   //!
   TBranch        *b_RECOPFPHOT_energySigmaPhiUp;   //!
   TBranch        *b_RECOPFPHOT_energySigmaPhiDown;   //!
   TBranch        *b_RECOPFPHOT_energySigmaRhoUp;   //!
   TBranch        *b_RECOPFPHOT_energySigmaRhoDown;   //!
   TBranch        *b_BeamSpot_X;   //!
   TBranch        *b_BeamSpot_Y;   //!
   TBranch        *b_BeamSpot_Z;   //!
   TBranch        *b_RECO_NVTX;   //!
   TBranch        *b_RECO_VERTEX_x;   //!
   TBranch        *b_RECO_VERTEX_y;   //!
   TBranch        *b_RECO_VERTEX_z;   //!
   TBranch        *b_RECO_VERTEX_ndof;   //!
   TBranch        *b_RECO_VERTEX_chi2;   //!
   TBranch        *b_RECO_VERTEX_ntracks;   //!
   TBranch        *b_RECO_VERTEXPROB;   //!
   TBranch        *b_RECO_VERTEX_isValid;   //!
   TBranch        *b_RECO_PFJET_N;   //!
   TBranch        *b_RECO_PFJET_CHARGE;   //!
   TBranch        *b_RECO_PFJET_ET;   //!
   TBranch        *b_RECO_PFJET_PT;   //!
   TBranch        *b_RECO_PFJET_ETA;   //!
   TBranch        *b_RECO_PFJET_PHI;   //!
   TBranch        *b_RECO_PFJET_PUID_loose;   //!
   TBranch        *b_RECO_PFJET_PUID_medium;   //!
   TBranch        *b_RECO_PFJET_PUID;   //!
   TBranch        *b_RECO_PFJET_PUID_MVA;   //!
   TBranch        *b_RECO_PFJET_QG_Likelihood;   //!
   TBranch        *b_RECO_PFJET_QG_axis2;   //!
   TBranch        *b_RECO_PFJET_QG_ptd;   //!
   TBranch        *b_RECO_PFJET_QG_mult;   //!
   TBranch        *b_RHO_ele;   //!
   TBranch        *b_RHO_mu;   //!
   TBranch        *b_LHE_PARTON_N;   //!
   TBranch        *b_LHE_PARTON_CLEAR;   //!
   TBranch        *b_LHE_PARTON_PDGID;   //!
   TBranch        *b_LHE_PARTON_PT;   //!
   TBranch        *b_LHE_PARTON_ETA;   //!
   TBranch        *b_LHE_PARTON_PHI;   //!
   TBranch        *b_LHE_PARTON_E;   //!
   TBranch        *b_RECO_PFJET_PT_UncUp;   //!
   TBranch        *b_RECO_PFJET_PT_UncDn;   //!
   TBranch        *b_RECO_PFJET_AREA;   //!
   TBranch        *b_RECO_PFJET_PTD;   //!
   TBranch        *b_RECO_PFJET_CHARGED_HADRON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_NEUTRAL_HADRON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_PHOTON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_ELECTRON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_MUON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_HF_HADRON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_HF_EM_ENERGY;   //!
   TBranch        *b_RECO_PFJET_CHARGED_EM_ENERGY;   //!
   TBranch        *b_RECO_PFJET_CHARGED_MU_ENERGY;   //!
   TBranch        *b_RECO_PFJET_NEUTRAL_EM_ENERGY;   //!
   TBranch        *b_RECO_PFJET_CHARGED_HADRON_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_PHOTON_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_ELECTRON_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_MUON_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_HF_HADRON_MULTIPLICTY;   //!
   TBranch        *b_RECO_PFJET_HF_EM_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_CHARGED_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_NEUTRAL_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_NCOMPONENTS;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_PDGID;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_PT;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_ETA;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_PHI;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_E;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_CHARGE;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_TRANSVERSE_MASS;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_XVERTEX;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_YVERTEX;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_ZVERTEX;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_VERTEX_CHI2;   //!
   TBranch        *b_RECO_CALOMET;   //!
   TBranch        *b_RECO_PFMET;   //!
   TBranch        *b_RECO_PFMET_X;   //!
   TBranch        *b_RECO_PFMET_Y;   //!
   TBranch        *b_RECO_PFMET_PHI;   //!
   TBranch        *b_RECO_PFMET_THETA;   //!
   TBranch        *b_RECO_PFMET_uncorr;   //!
   TBranch        *b_RECO_PFMET_X_uncorr;   //!
   TBranch        *b_RECO_PFMET_Y_uncorr;   //!
   TBranch        *b_RECO_PFMET_PHI_uncorr;   //!
   TBranch        *b_RECO_PFMET_THETA_uncorr;   //!
   TBranch        *b_RECO_PFMET_JetEnUp;   //!
   TBranch        *b_RECO_PFMET_JetEnDn;   //!
   TBranch        *b_RECO_PFMET_ElectronEnUp;   //!
   TBranch        *b_RECO_PFMET_ElectronEnDn;   //!
   TBranch        *b_RECO_PFMET_MuonEnUp;   //!
   TBranch        *b_RECO_PFMET_MuonEnDn;   //!
   TBranch        *b_RECO_PFMET_JetResUp;   //!
   TBranch        *b_RECO_PFMET_JetResDn;   //!
   TBranch        *b_RECO_PFMET_UnclusteredEnUp;   //!
   TBranch        *b_RECO_PFMET_UnclusteredEnDn;   //!
   TBranch        *b_RECO_PFMET_PhotonEnUp;   //!
   TBranch        *b_RECO_PFMET_PhotonEnDn;   //!
   TBranch        *b_RECO_PFMET_TauEnUp ;   //!
   TBranch        *b_RECO_PFMET_TauEnDown;   //!
   TBranch        *b_RECO_PFMET_GoodVtxNoiseFilter;   //!
   TBranch        *b_RECO_PFMET_GlobalSuperTightHalo2016NoiseFilter;   //!
   TBranch        *b_RECO_PFMET_HBHENoiseFilter;   //!
   TBranch        *b_RECO_PFMET_HBHENoiseIsoFilter;   //!
   TBranch        *b_RECO_PFMET_EcalDeadCellTriggerPrimitiveNoiseFilter;   //!
   TBranch        *b_RECO_PFMET_BadPFMuonFilter;   //!
   TBranch        *b_RECO_PFMET_BadChargedCandidateFilter;   //!
   TBranch        *b_RECO_PFMET_EEBadScNoiseFilter;   //!
   TBranch        *b_RECO_PFMET_EcalBadCalibFilter;   //!
   TBranch        *b_RECO_TCMET;   //!
   TBranch        *b_RECO_CORMETMUONS;   //!
   TBranch        *b_tCHighEff_BTagJet_PT;   //!
   TBranch        *b_tCHighPur_BTagJet_PT;   //!
   TBranch        *b_cSV_BTagJet_PT;   //!
   TBranch        *b_tCHighEff_BTagJet_ETA;   //!
   TBranch        *b_tCHighPur_BTagJet_ETA;   //!
   TBranch        *b_cSV_BTagJet_ETA;   //!
   TBranch        *b_tCHighEff_BTagJet_PHI;   //!
   TBranch        *b_tCHighPur_BTagJet_PHI;   //!
   TBranch        *b_cSV_BTagJet_PHI;   //!
   TBranch        *b_cSV_BTagJet_ET;   //!
   TBranch        *b_tCHighEff_BTagJet_DISCR;   //!
   TBranch        *b_tCHighPur_BTagJet_DISCR;   //!
   TBranch        *b_cSV_BTagJet_DISCR;   //!

   HZZ4LeptonsAnalysis(TTree *tree=0,Double_t weight_=1.,std::string DATA_type_="DATA",std::string MC_type_="MC");
   virtual ~HZZ4LeptonsAnalysis();
   Double_t weight;
   std::string DATA_type,MC_type;
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const Char_t *name);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   std::ofstream bnn_file;
   double EAele(int ,bool );
   double masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr );
   #if (product_4e || product_2e2mu)
   void printelebnn(int i);
   #endif
   #if (product_4mu || product_2e2mu)
   void printmubnn(int i);
   #endif
   float kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ, int finalState);
   float kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState);
   float kfactor_qqZZ_qcd_Pt(float GENpTZZ, int finalState);
   float kfactor_ggZZ(float GENmassZZ, int finalState);
};

#endif

#ifdef HZZ4LeptonsAnalysis_cxx
HZZ4LeptonsAnalysis::HZZ4LeptonsAnalysis(TTree *tree,Double_t weight_, std::string DATA_type_, std::string MC_type_) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   weight = weight_;
   DATA_type = DATA_type_;
   MC_type = MC_type_;

   if (tree == 0) { 
      exit(9);
      TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
      chain->Add("/lustre/cms/store/user/ndefilip/Summer12_52X_merged/roottree_leptons_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.root_a");
      chain->Add("/lustre/cms/store/user/ndefilip/Summer12_52X_merged/roottree_leptons_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_1.root_b");

      tree = chain;

   }
   Init(tree);
}

HZZ4LeptonsAnalysis::~HZZ4LeptonsAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HZZ4LeptonsAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HZZ4LeptonsAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HZZ4LeptonsAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   RECOMU_PT_MuHLTMatch = 0;
   RECOMU_ETA_MuHLTMatch = 0;
   RECOMU_PHI_MuHLTMatch = 0;
   RECOELE_PT_EleHLTMatch = 0;
   RECOELE_ETA_EleHLTMatch = 0;
   RECOELE_PHI_EleHLTMatch = 0;
   MC_E = 0;
   MC_PT = 0;
   MC_ETA = 0;
   MC_THETA = 0;
   MC_PHI = 0;
   MC_MASS = 0;
   MC_PDGID = 0;
   MC_LEPT_PT = 0;
   MC_LEPT_ETA = 0;
   MC_LEPT_PHI = 0;
   MC_LEPT_THETA = 0;
   MC_LEPT_PDGID = 0;
   MC_Z_PT = 0;
   MC_Z_ETA = 0;
   MC_Z_PHI = 0;
   MC_Z_THETA = 0;
   MC_Z_MASS = 0;
   MC_Z_PDGID = 0;
   MC_fourl_MASS = 0;
   MC_fourl_PT = 0;
   MC_fourl_PDGID = 0;
   MC_ZZ_MASS = 0;
   MC_ZZ_PT = 0;
   MC_ZZ_ETA = 0;
   MC_ZZ_PHI = 0;
   MC_ZZ_THETA = 0;
   MC_ZZ_PDGID = 0;
   MC_ZZ_D_MASS = 0;
   MC_ZZ_D_PT = 0;
   MC_ZZ_D_ETA = 0;
   MC_ZZ_D_PHI = 0;
   MC_ZZ_D_THETA = 0;
   MC_ZZ_D_PDGID = 0;
   MC_ZZ_DD_MASS = 0;
   MC_ZZ_DD_PT = 0;
   MC_ZZ_DD_ETA = 0;
   MC_ZZ_DD_PHI = 0;
   MC_ZZ_DD_THETA = 0;
   MC_ZZ_DD_PDGID = 0;
   MC_GENJET_PT = 0;
   MC_GENJET_ETA = 0;
   MC_GENJET_PHI = 0;
   RECORF_2e2mu_cosTheta1_spin = 0;
   RECORF_2e2mu_cosTheta2_spin = 0;
   RECORF_2e2mu_cosThetaStar_spin = 0;
   RECORF_2e2mu_Phi_spin = 0;
   RECORF_2e2mu_Phi1_spin = 0;
   RECORF_2e2mu_Phi2_spin = 0;
   RECORF_2e2mu_phi1RF_spin = 0;
   RECORF_2e2mu_phi2RF_spin = 0;
   RECORF_2e2mu_MELA = 0;
   RECORF_4e_cosTheta1_spin = 0;
   RECORF_4e_cosTheta2_spin = 0;
   RECORF_4e_cosThetaStar_spin = 0;
   RECORF_4e_Phi_spin = 0;
   RECORF_4e_Phi1_spin = 0;
   RECORF_4e_Phi2_spin = 0;
   RECORF_4e_phi1RF_spin = 0;
   RECORF_4e_phi2RF_spin = 0;
   RECORF_4e_MELA = 0;
   RECORF_4mu_cosTheta1_spin = 0;
   RECORF_4mu_cosTheta2_spin = 0;
   RECORF_4mu_cosThetaStar_spin = 0;
   RECORF_4mu_Phi_spin = 0;
   RECORF_4mu_Phi1_spin = 0;
   RECORF_4mu_Phi2_spin = 0;
   RECORF_4mu_phi1RF_spin = 0;
   RECORF_4mu_phi2RF_spin = 0;
   RECORF_4mu_MELA = 0;
   RECO_ZMM_MASS = 0;
   RECO_ZEE_MASS = 0;
   RECO_DiLep_MASS = 0;
   RECO_ZMM_PT = 0;
   RECO_ZEE_PT = 0;
   RECO_DiLep_PT = 0;
   RECO_ZMM_ETA = 0;
   RECO_ZEE_ETA = 0;
   RECO_DiLep_ETA = 0;
   RECO_ZMM_PHI = 0;
   RECO_ZEE_PHI = 0;
   RECO_DiLep_PHI = 0;
   RECO_ZMMss_MASS = 0;
   RECO_ZEEss_MASS = 0;
   RECO_ZEM_MASS = 0;
   RECO_ZMMss_PT = 0;
   RECO_ZEEss_PT = 0;
   RECO_ZEM_PT = 0;
   RECO_ZMMss_ETA = 0;
   RECO_ZEEss_ETA = 0;
   RECO_ZEM_ETA = 0;
   RECO_ZMMss_PHI = 0;
   RECO_ZEEss_PHI = 0;
   RECO_ZEM_PHI = 0;
   RECO_MMMM_MASS = 0;
   RECO_MMMM_PT = 0;
   RECO_MMMM_ETA = 0;
   RECO_MMMM_PHI = 0;
   RECO_MMMM_MASS_REFIT = 0;
   RECO_EEEE_MASS = 0;
   RECO_EEEE_PT = 0;
   RECO_EEEE_ETA = 0;
   RECO_EEEE_PHI = 0;
   RECO_EEEE_MASS_REFIT = 0;
   RECO_EEMM_MASS = 0;
   RECO_EEMM_PT = 0;
   RECO_EEMM_ETA = 0;
   RECO_EEMM_PHI = 0;
   RECO_EEMM_MASS_REFIT = 0;
   RECO_LLL0_MASS = 0;
   RECO_LLL1_MASS = 0;
   RECO_LLL2_MASS = 0;
   RECO_LLL3_MASS = 0;
   RECO_LLL0_PT = 0;
   RECO_LLL1_PT = 0;
   RECO_LLL2_PT = 0;
   RECO_LLL3_PT = 0;
   RECO_LLLl0_MASS = 0;
   RECO_LLLl1_MASS = 0;
   RECO_LLLl0_PT = 0;
   RECO_LLLl1_PT = 0;
   RECO_LLLL0ss_MASS = 0;
   RECO_LLLL0ss_PT = 0;
   RECO_LLLL1ss_MASS = 0;
   RECO_LLLL1ss_PT = 0;
   RECO_LLLL2ss_MASS = 0;
   RECO_LLLL2ss_PT = 0;
   RECO_LLLL_MASS = 0;
   RECO_LLLL_PT = 0;
   RECO_LLLL_ETA = 0;
   RECO_LLLL_PHI = 0;
   RECOELE_E = 0;
   RECOELE_PT = 0;
   RECOELE_PTError = 0;
   RECOELE_P = 0;
   RECOELE_ETA = 0;
   RECOELE_THETA = 0;
   RECOELE_PHI = 0;
   RECOELE_MASS = 0;
   RECOELE_CHARGE = 0;
   RECOELE_ID = 0;
   RECOELE_PT_uncorr = 0;
   RECOELE_isEcalDriven = 0;
   RECOELE_isTrackerDriven = 0;
   RECOELE_gsftrack_NPixHits = 0;
   RECOELE_gsftrack_NStripHits = 0;
   RECOELE_gsftrack_chi2 = 0;
   RECOELE_gsftrack_dxyB = 0;
   RECOELE_gsftrack_dxy = 0;
   RECOELE_gsftrack_dxyError = 0;
   RECOELE_gsftrack_dzB = 0;
   RECOELE_gsftrack_dz = 0;
   RECOELE_gsftrack_dzError = 0;
   RECOELE_gsftrack_losthits = 0;
   RECOELE_gsftrack_validhits = 0;
   RECOELE_gsftrack_expected_inner_hits = 0;
   RECOELE_scl_E = 0;
   RECOELE_scl_Et = 0;
   RECOELE_scl_Eta = 0;
   RECOELE_scl_Phi = 0;
   RECOELE_ecalEnergy = 0;
   RECOELE_ep = 0;
   RECOELE_eSeedp = 0;
   RECOELE_eSeedpout = 0;
   RECOELE_eElepout = 0;
   RECOELE_deltaEtaIn = 0;
   RECOELE_deltaEtaSeed = 0;
   RECOELE_deltaEtaEle = 0;
   RECOELE_deltaPhiIn = 0;
   RECOELE_deltaPhiSeed = 0;
   RECOELE_deltaPhiEle = 0;
   RECOELE_isbarrel = 0;
   RECOELE_isendcap = 0;
   RECOELE_isGap = 0;
   RECOELE_isEBetaGap = 0;
   RECOELE_isEBphiGap = 0;
   RECOELE_isEEdeeGap = 0;
   RECOELE_isEEringGap = 0;
   RECOELE_sigmaIetaIeta = 0;
   RECOELE_sigmaEtaEta = 0;
   RECOELE_e15 = 0;
   RECOELE_e25max = 0;
   RECOELE_e55 = 0;
   RECOELE_he = 0;
   RECOELE_r9 = 0;
   RECOELE_mva = 0;
   RECOELE_fbrem = 0;
   RECOELE_nbrems = 0;
   RECOELE_Class = 0;
   RECOELE_fbrem_mode = 0;
   RECOELE_fbrem_mean = 0;
   RECOELE_EGMTRACKISO = 0;
   RECOELE_EGMHCALISO = 0;
   RECOELE_EGMECALISO = 0;
   RECOELE_EGMX = 0;
   RECOELE_PFchAllPart = 0;
   RECOELE_PFchHad = 0;
   RECOELE_PFneuHad = 0;
   RECOELE_PFphoton = 0;
   RECOELE_PFPUchAllPart = 0;
   RECOELE_PFX_dB = 0;
   RECOELE_PFX_rho = 0;
   RECOELE_regEnergy = 0;
   RECOELE_regEnergyError = 0;
   RECOELE_SIP = 0;
   RECOELE_IP = 0;
   RECOELE_IPERROR = 0;
   RECOELE_SIP_KF = 0;
   RECOELE_IP_KF = 0;
   RECOELE_IPERROR_KF = 0;
   RECOELE_SIP_GD = 0;
   RECOELE_SIP_GDEEEE = 0;
   RECOELE_SIP_Std = 0;
   RECOELE_SIP_StdEEEE = 0;
   RECOELE_SIP_Kin = 0;
   RECOELE_SIP_KinEEEE = 0;
   RECOELE_STIP = 0;
   RECOELE_SLIP = 0;
   RECOELE_TIP = 0;
   RECOELE_LIP = 0;
   RECOELE_TIPERROR = 0;
   RECOELE_LIPERROR = 0;
   RECOELE_sclRawE = 0;
   RECOELE_sclX = 0;
   RECOELE_sclY = 0;
   RECOELE_sclZ = 0;
   RECOELE_seedSubdet1 = 0;
   RECOELE_seedDphi1 = 0;
   RECOELE_seedDrz1 = 0;
   RECOELE_seedSubdet2 = 0;
   RECOELE_seedDphi2 = 0;
   RECOELE_seedDrz2 = 0;
   RECOELE_eidVeryLoose = 0;
   RECOELE_eidLoose = 0;
   RECOELE_eidMedium = 0;
   RECOELE_eidTight = 0;
   RECOELE_eidHZZVeryLoose = 0;
   RECOELE_eidHZZLoose = 0;
   RECOELE_eidHZZMedium = 0;
   RECOELE_eidHZZTight = 0;
   RECOELE_mvaTrigV0 = 0;
   RECOELE_mvaNonTrigV0 = 0;
   RECOELE_COV = 0;
   RECOELE_TLE_ParentSC_X = 0;
   RECOELE_TLE_ParentSC_Y = 0;
   RECOELE_TLE_ParentSC_Z = 0;
   RECOELE_ecalTrkEnergyPreCorr = 0;
   RECOELE_ecalTrkEnergyErrPreCorr = 0;
   RECOELE_ecalTrkEnergyErrPostCorr = 0;
   RECOELE_energyScaleValue = 0;
   RECOELE_energySigmaValue = 0;
   RECOELE_energyScaleUp = 0;
   RECOELE_energyScaleDown = 0;
   RECOELE_energyScaleStatUp = 0;
   RECOELE_energyScaleStatDown = 0;
   RECOELE_energyScaleSystUp = 0;
   RECOELE_energyScaleSystDown = 0;
   RECOELE_energyScaleGainUp = 0;
   RECOELE_energyScaleGainDown = 0;
   RECOELE_energyScaleEtUp = 0;
   RECOELE_energyScaleEtDown = 0;
   RECOELE_energySigmaUp = 0;
   RECOELE_energySigmaDown = 0;
   RECOELE_energySigmaPhiUp = 0;
   RECOELE_energySigmaPhiDown = 0;
   RECOELE_energySigmaRhoUp = 0;
   RECOELE_energySigmaRhoDown = 0;
   RECOMU_isPFMu = 0;
   RECOMU_isGlobalMu = 0;
   RECOMU_isStandAloneMu = 0;
   RECOMU_isTrackerMu = 0;
   RECOMU_isCaloMu = 0;
   RECOMU_isTrackerHighPtMu = 0;
   RECOMU_isME0Muon = 0;
   RECOMU_E = 0;
   RECOMU_PT = 0;
   RECOMU_P = 0;
   RECOMU_ETA = 0;
   RECOMU_THETA = 0;
   RECOMU_PHI = 0;
   RECOMU_MASS = 0;
   RECOMU_CHARGE = 0;
   RECOMU_PT_uncorr = 0;
   RECOMU_COV = 0;
   RECOMU_TRACKISO = 0;
   RECOMU_TRACKISO_SUMPT = 0;
   RECOMU_HCALISO = 0;
   RECOMU_ECALISO = 0;
   RECOMU_X = 0;
   RECOMU_PFchHad = 0;
   RECOMU_PFneuHad = 0;
   RECOMU_PFphoton = 0;
   RECOMU_PFPUchAllPart = 0;
   RECOMU_PFX_dB = 0;
   RECOMU_PFX_rho = 0;
   RECOPFPHOT_PFchHad = 0;
   RECOPFPHOT_PFneuHad = 0;
   RECOPFPHOT_PFphoton = 0;
   RECOPFPHOT_PFPUchAllPart = 0;
   RECOPFPHOT_PFX_rho = 0;
   RECOMU_SIP = 0;
   RECOMU_IP = 0;
   RECOMU_IPERROR = 0;
   RECOMU_SIP_KF = 0;
   RECOMU_IP_KF = 0;
   RECOMU_IPERROR_KF = 0;
   RECOMU_SIP_GD = 0;
   RECOMU_SIP_GDMMMM = 0;
   RECOMU_SIP_Std = 0;
   RECOMU_SIP_StdMMMM = 0;
   RECOMU_SIP_Kin = 0;
   RECOMU_SIP_KinMMMM = 0;
   RECOMU_STIP = 0;
   RECOMU_SLIP = 0;
   RECOMU_TIP = 0;
   RECOMU_LIP = 0;
   RECOMU_TIPERROR = 0;
   RECOMU_LIPERROR = 0;
   RECOMU_caloCompatibility = 0;
   RECOMU_segmentCompatibility = 0;
   RECOMU_numberOfMatches = 0;
   RECOMU_numberOfMatchedStations = 0;
   RECOMU_glbmuPromptTight = 0;
   RECOMU_trkmuArbitration = 0;
   RECOMU_trkmu2DCompatibilityLoose = 0;
   RECOMU_trkmu2DCompatibilityTight = 0;
   RECOMU_trkmuOneStationLoose = 0;
   RECOMU_trkmuOneStationTight = 0;
   RECOMU_trkmuLastStationLoose = 0;
   RECOMU_trkmuLastStationTight = 0;
   RECOMU_trkmuOneStationAngLoose = 0;
   RECOMU_trkmuOneStationAngTight = 0;
   RECOMU_trkmuLastStationAngLoose = 0;
   RECOMU_trkmuLastStationAngTight = 0;
   RECOMU_trkmuLastStationOptimizedLowPtLoose = 0;
   RECOMU_trkmuLastStationOptimizedLowPtTight = 0;
   RECOMU_mutrkPT = 0;
   RECOMU_mutrkPTError = 0;
   RECOMU_mutrkDxy = 0;
   RECOMU_mutrkDxyError = 0;
   RECOMU_mutrkDxyB = 0;
   RECOMU_mutrkDz = 0;
   RECOMU_mutrkDzError = 0;
   RECOMU_mutrkDzB = 0;
   RECOMU_mutrkChi2PerNdof = 0;
   RECOMU_mutrkCharge = 0;
   RECOMU_mutrkNHits = 0;
   RECOMU_mutrkNStripHits = 0;
   RECOMU_mutrkNPixHits = 0;
   RECOMU_mutrkNMuonHits = 0;
   RECOMU_mutrktrackerLayersWithMeasurement = 0;
   RECOMU_muInnertrkDxy = 0;
   RECOMU_muInnertrkDxyError = 0;
   RECOMU_muInnertrkDxyB = 0;
   RECOMU_muInnertrkDz = 0;
   RECOMU_muInnertrkDzError = 0;
   RECOMU_muInnertrkDzB = 0;
   RECOMU_muInnertrkChi2PerNdof = 0;
   RECOMU_muInnertrktrackerLayersWithMeasurement = 0;
   RECOMU_muInnertrkPT = 0;
   RECOMU_muInnertrkPTError = 0;
   RECOMU_muInnertrkCharge = 0;
   RECOMU_muInnertrkNHits = 0;
   RECOMU_muInnertrkNStripHits = 0;
   RECOMU_muInnertrkNPixHits = 0;
   RECOMU_mubesttrkType = 0;
   RECOMU_mubesttrkDxy = 0;
   RECOMU_mubesttrkDxyError = 0;
   RECOMU_mubesttrkDz = 0;
   RECOMU_mubesttrkDzError = 0;
   RECOMU_mubesttrkPTError = 0;
   RECOMU_Rochester_Error = 0;
   ftsigma = 0;
   gdX = 0;
   gdY = 0;
   gdZ = 0;
   ftsigmalag = 0;
   gdlagX = 0;
   gdlagY = 0;
   gdlagZ = 0;
   gdlagProb = 0;
   gdlagNdof = 0;
   ftsigmaMMMM = 0;
   gdXMMMM = 0;
   gdYMMMM = 0;
   gdZMMMM = 0;
   ftsigmalagMMMM = 0;
   gdlagXMMMM = 0;
   gdlagYMMMM = 0;
   gdlagZMMMM = 0;
   gdlagProbMMMM = 0;
   gdlagNdofMMMM = 0;
   ftsigmaEEEE = 0;
   gdXEEEE = 0;
   gdYEEEE = 0;
   gdZEEEE = 0;
   ftsigmalagEEEE = 0;
   gdlagXEEEE = 0;
   gdlagYEEEE = 0;
   gdlagZEEEE = 0;
   gdlagProbEEEE = 0;
   gdlagNdofEEEE = 0;
   StdFitVertexX = 0;
   StdFitVertexY = 0;
   StdFitVertexZ = 0;
   StdFitVertexChi2r = 0;
   StdFitVertexProb = 0;
   StdFitVertexTrack_PT = 0;
   StdFitVertexTrack_ETA = 0;
   StdFitVertexTrack_PHI = 0;
   KinFitVertexX = 0;
   KinFitVertexY = 0;
   KinFitVertexZ = 0;
   KinFitVertexChi2r = 0;
   KinFitVertexProb = 0;
   StdFitVertexXMMMM = 0;
   StdFitVertexYMMMM = 0;
   StdFitVertexZMMMM = 0;
   StdFitVertexChi2rMMMM = 0;
   StdFitVertexProbMMMM = 0;
   StdFitVertexTrackMMMM_PT = 0;
   StdFitVertexTrackMMMM_ETA = 0;
   StdFitVertexTrackMMMM_PHI = 0;
   KinFitVertexXMMMM = 0;
   KinFitVertexYMMMM = 0;
   KinFitVertexZMMMM = 0;
   KinFitVertexChi2rMMMM = 0;
   KinFitVertexProbMMMM = 0;
   StdFitVertexXEEEE = 0;
   StdFitVertexYEEEE = 0;
   StdFitVertexZEEEE = 0;
   StdFitVertexChi2rEEEE = 0;
   StdFitVertexProbEEEE = 0;
   StdFitVertexTrackEEEE_PT = 0;
   StdFitVertexTrackEEEE_ETA = 0;
   StdFitVertexTrackEEEE_PHI = 0;
   KinFitVertexXEEEE = 0;
   KinFitVertexYEEEE = 0;
   KinFitVertexZEEEE = 0;
   KinFitVertexChi2rEEEE = 0;
   KinFitVertexProbEEEE = 0;
   StdFitVertexChi2rMMM = 0;
   StdFitVertexProbMMM = 0;
   StdFitVertexChi2rMME = 0;
   StdFitVertexProbMME = 0;
   StdFitVertexChi2rEEE = 0;
   StdFitVertexProbEEE = 0;
   StdFitVertexChi2rMEE = 0;
   StdFitVertexProbMEE = 0;
   StdFitVertexChi2rDiLep = 0;
   StdFitVertexProbDiLep = 0;
   ConvMapDist = 0;
   ConvMapDcot = 0;
   RECOMU_MatchingMCTruth = 0;
   RECOMU_MatchingMCpT = 0;
   RECOMU_MatchingMCEta = 0;
   RECOMU_MatchingMCPhi = 0;
   RECOELE_MatchingMCTruth = 0;
   RECOELE_MatchingMCpT = 0;
   RECOELE_MatchingMCEta = 0;
   RECOELE_MatchingMCPhi = 0;
   RECOPHOT_MatchingMCTruth = 0;
   RECOPHOT_MatchingMCpT = 0;
   RECOPHOT_MatchingMCEta = 0;
   RECOPHOT_MatchingMCPhi = 0;
   RECOzMuMu_MatchingMCTruth = 0;
   RECOzMuMu_MatchingMCpT = 0;
   RECOzMuMu_MatchingMCmass = 0;
   RECOzMuMu_MatchingMCEta = 0;
   RECOzMuMu_MatchingMCPhi = 0;
   RECOzEE_MatchingMCTruth = 0;
   RECOzEE_MatchingMCpT = 0;
   RECOzEE_MatchingMCmass = 0;
   RECOzEE_MatchingMCEta = 0;
   RECOzEE_MatchingMCPhi = 0;
   RECOHzzEEEE_MatchingMCTruth = 0;
   RECOHzzEEEE_MatchingMCpT = 0;
   RECOHzzEEEE_MatchingMCmass = 0;
   RECOHzzEEEE_MatchingMCEta = 0;
   RECOHzzEEEE_MatchingMCPhi = 0;
   RECOHzzEEMM_MatchingMCTruth = 0;
   RECOHzzEEMM_MatchingMCpT = 0;
   RECOHzzEEMM_MatchingMCmass = 0;
   RECOHzzEEMM_MatchingMCEta = 0;
   RECOHzzEEMM_MatchingMCPhi = 0;
   RECOHzzMMMM_MatchingMCTruth = 0;
   RECOHzzMMMM_MatchingMCpT = 0;
   RECOHzzMMMM_MatchingMCmass = 0;
   RECOHzzMMMM_MatchingMCEta = 0;
   RECOHzzMMMM_MatchingMCPhi = 0;
   RECOPHOT_PT = 0;
   RECOPHOT_ETA = 0;
   RECOPHOT_PHI = 0;
   RECOPHOT_THETA = 0;
   RECOPHOT_TLE_ParentSC_X = 0;
   RECOPHOT_TLE_ParentSC_Y = 0;
   RECOPHOT_TLE_ParentSC_Z = 0;
   RECOPFPHOT_PT = 0;
   RECOPFPHOT_PTError = 0;
   RECOPFPHOT_ETA = 0;
   RECOPFPHOT_PHI = 0;
   RECOPFPHOT_THETA = 0;
   RECOPFPHOT_PT_uncorr = 0;
   RECOPFPHOT_ecalEnergyPreCorr = 0;
   RECOPFPHOT_ecalEnergyErrPreCorr = 0;
   RECOPFPHOT_ecalEnergyErrPostCorr = 0;
   RECOPFPHOT_energyScaleValue = 0;
   RECOPFPHOT_energySigmaValue = 0;
   RECOPFPHOT_energyScaleUp = 0;
   RECOPFPHOT_energyScaleDown = 0;
   RECOPFPHOT_energyScaleStatUp = 0;
   RECOPFPHOT_energyScaleStatDown = 0;
   RECOPFPHOT_energyScaleSystUp = 0;
   RECOPFPHOT_energyScaleSystDown = 0;
   RECOPFPHOT_energyScaleGainUp = 0;
   RECOPFPHOT_energyScaleGainDown = 0;
   RECOPFPHOT_energyScaleEtUp = 0;
   RECOPFPHOT_energyScaleEtDown = 0;
   RECOPFPHOT_energySigmaUp = 0;
   RECOPFPHOT_energySigmaDown = 0;
   RECOPFPHOT_energySigmaPhiUp = 0;
   RECOPFPHOT_energySigmaPhiDown = 0;
   RECOPFPHOT_energySigmaRhoUp = 0;
   RECOPFPHOT_energySigmaRhoDown = 0;
   RECO_PFJET_CHARGE = 0;
   RECO_PFJET_ET = 0;
   RECO_PFJET_PT = 0;
   RECO_PFJET_ETA = 0;
   RECO_PFJET_PHI = 0;
   RECO_PFJET_PUID_loose = 0;
   RECO_PFJET_PUID_medium = 0;
   RECO_PFJET_PUID = 0;
   RECO_PFJET_PUID_MVA = 0;
   RECO_PFJET_QG_Likelihood = 0;
   RECO_PFJET_QG_axis2 = 0;
   RECO_PFJET_QG_ptd = 0;
   RECO_PFJET_QG_mult = 0;
   LHE_PARTON_CLEAR = 0;
   LHE_PARTON_PDGID = 0;
   LHE_PARTON_PT = 0;
   LHE_PARTON_ETA = 0;
   LHE_PARTON_PHI = 0;
   LHE_PARTON_E = 0;
   RECO_PFJET_PT_UncUp = 0;
   RECO_PFJET_PT_UncDn = 0;
   RECO_PFJET_AREA = 0;
   RECO_PFJET_PTD = 0;
   RECO_PFJET_CHARGED_HADRON_ENERGY = 0;
   RECO_PFJET_NEUTRAL_HADRON_ENERGY = 0;
   RECO_PFJET_PHOTON_ENERGY = 0;
   RECO_PFJET_ELECTRON_ENERGY = 0;
   RECO_PFJET_MUON_ENERGY = 0;
   RECO_PFJET_HF_HADRON_ENERGY = 0;
   RECO_PFJET_HF_EM_ENERGY = 0;
   RECO_PFJET_CHARGED_EM_ENERGY = 0;
   RECO_PFJET_CHARGED_MU_ENERGY = 0;
   RECO_PFJET_NEUTRAL_EM_ENERGY = 0;
   RECO_PFJET_CHARGED_HADRON_MULTIPLICITY = 0;
   RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY = 0;
   RECO_PFJET_PHOTON_MULTIPLICITY = 0;
   RECO_PFJET_ELECTRON_MULTIPLICITY = 0;
   RECO_PFJET_MUON_MULTIPLICITY = 0;
   RECO_PFJET_HF_HADRON_MULTIPLICTY = 0;
   RECO_PFJET_HF_EM_MULTIPLICITY = 0;
   RECO_PFJET_CHARGED_MULTIPLICITY = 0;
   RECO_PFJET_NEUTRAL_MULTIPLICITY = 0;
   RECO_PFJET_NCOMPONENTS = 0;
   RECO_PFJET_COMPONENT_PDGID = 0;
   RECO_PFJET_COMPONENT_PT = 0;
   RECO_PFJET_COMPONENT_ETA = 0;
   RECO_PFJET_COMPONENT_PHI = 0;
   RECO_PFJET_COMPONENT_E = 0;
   RECO_PFJET_COMPONENT_CHARGE = 0;
   RECO_PFJET_COMPONENT_TRANSVERSE_MASS = 0;
   RECO_PFJET_COMPONENT_XVERTEX = 0;
   RECO_PFJET_COMPONENT_YVERTEX = 0;
   RECO_PFJET_COMPONENT_ZVERTEX = 0;
   RECO_PFJET_COMPONENT_VERTEX_CHI2 = 0;
   tCHighEff_BTagJet_PT = 0;
   tCHighPur_BTagJet_PT = 0;
   cSV_BTagJet_PT = 0;
   tCHighEff_BTagJet_ETA = 0;
   tCHighPur_BTagJet_ETA = 0;
   cSV_BTagJet_ETA = 0;
   tCHighEff_BTagJet_PHI = 0;
   tCHighPur_BTagJet_PHI = 0;
   cSV_BTagJet_PHI = 0;
   cSV_BTagJet_ET = 0;
   tCHighEff_BTagJet_DISCR = 0;
   tCHighPur_BTagJet_DISCR = 0;
   cSV_BTagJet_DISCR = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_irun);
   fChain->SetBranchAddress("Event", &Event, &b_ievt);
   fChain->SetBranchAddress("LumiSection", &LumiSection, &b_ils);
   fChain->SetBranchAddress("Avginstlumi", &Avginstlumi, &b_Avginstlumi);
   fChain->SetBranchAddress("num_PU_vertices", &num_PU_vertices, &b_num_PU_vertices);
   fChain->SetBranchAddress("PU_BunchCrossing", &PU_BunchCrossing, &b_PU_BunchCrossing);
   fChain->SetBranchAddress("MC_weighting", &MC_weighting, &b_MC_weighting);
   fChain->SetBranchAddress("RECO_nMuHLTMatch", &RECO_nMuHLTMatch, &b_RECO_nMuHLTMatch);
   fChain->SetBranchAddress("RECOMU_PT_MuHLTMatch", &RECOMU_PT_MuHLTMatch, &b_RECOMU_PT_MuHLTMatch);
   fChain->SetBranchAddress("RECOMU_ETA_MuHLTMatch", &RECOMU_ETA_MuHLTMatch, &b_RECOMU_ETA_MuHLTMatch);
   fChain->SetBranchAddress("RECOMU_PHI_MuHLTMatch", &RECOMU_PHI_MuHLTMatch, &b_RECOMU_PHI_MuHLTMatch);
   fChain->SetBranchAddress("RECO_nEleHLTMatch", &RECO_nEleHLTMatch, &b_RECO_nEleHLTMatch);
   fChain->SetBranchAddress("RECOELE_PT_EleHLTMatch", &RECOELE_PT_EleHLTMatch, &b_RECOELE_PT_EleHLTMatch);
   fChain->SetBranchAddress("RECOELE_ETA_EleHLTMatch", &RECOELE_ETA_EleHLTMatch, &b_RECOELE_ETA_EleHLTMatch);
   fChain->SetBranchAddress("RECOELE_PHI_EleHLTMatch", &RECOELE_PHI_EleHLTMatch, &b_RECOELE_PHI_EleHLTMatch);
   fChain->SetBranchAddress("HLTPathsFired", &HLTPathsFired, &b_HLTPathsFired);
   fChain->SetBranchAddress("MC_E", &MC_E, &b_MC_E);
   fChain->SetBranchAddress("MC_PT", &MC_PT, &b_MC_PT);
   fChain->SetBranchAddress("MC_ETA", &MC_ETA, &b_MC_ETA);
   fChain->SetBranchAddress("MC_THETA", &MC_THETA, &b_MC_THETA);
   fChain->SetBranchAddress("MC_PHI", &MC_PHI, &b_MC_PHI);
   fChain->SetBranchAddress("MC_MASS", &MC_MASS, &b_MC_MASS);
   fChain->SetBranchAddress("MC_PDGID", &MC_PDGID, &b_MC_PDGID);
   fChain->SetBranchAddress("MC_LEPT_PT", &MC_LEPT_PT, &b_MC_LEPT_PT);
   fChain->SetBranchAddress("MC_LEPT_ETA", &MC_LEPT_ETA, &b_MC_LEPT_ETA);
   fChain->SetBranchAddress("MC_LEPT_PHI", &MC_LEPT_PHI, &b_MC_LEPT_PHI);
   fChain->SetBranchAddress("MC_LEPT_THETA", &MC_LEPT_THETA, &b_MC_LEPT_THETA);
   fChain->SetBranchAddress("MC_LEPT_PDGID", &MC_LEPT_PDGID, &b_MC_LEPT_PDGID);
   fChain->SetBranchAddress("MC_Z_PT", &MC_Z_PT, &b_MC_Z_PT);
   fChain->SetBranchAddress("MC_Z_ETA", &MC_Z_ETA, &b_MC_Z_ETA);
   fChain->SetBranchAddress("MC_Z_PHI", &MC_Z_PHI, &b_MC_Z_PHI);
   fChain->SetBranchAddress("MC_Z_THETA", &MC_Z_THETA, &b_MC_Z_THETA);
   fChain->SetBranchAddress("MC_Z_MASS", &MC_Z_MASS, &b_MC_Z_MASS);
   fChain->SetBranchAddress("MC_Z_PDGID", &MC_Z_PDGID, &b_MC_Z_PDGID);
   fChain->SetBranchAddress("MC_fourl_MASS", &MC_fourl_MASS, &b_MC_fourl_MASS);
   fChain->SetBranchAddress("MC_fourl_PT", &MC_fourl_PT, &b_MC_fourl_PT);
   fChain->SetBranchAddress("MC_fourl_PDGID", &MC_fourl_PDGID, &b_MC_fourl_PDGID);
   fChain->SetBranchAddress("MC_ZZ_MASS", &MC_ZZ_MASS, &b_MC_ZZ_MASS);
   fChain->SetBranchAddress("MC_ZZ_PT", &MC_ZZ_PT, &b_MC_ZZ_PT);
   fChain->SetBranchAddress("MC_ZZ_ETA", &MC_ZZ_ETA, &b_MC_ZZ_ETA);
   fChain->SetBranchAddress("MC_ZZ_PHI", &MC_ZZ_PHI, &b_MC_ZZ_PHI);
   fChain->SetBranchAddress("MC_ZZ_THETA", &MC_ZZ_THETA, &b_MC_ZZ_THETA);
   fChain->SetBranchAddress("MC_ZZ_PDGID", &MC_ZZ_PDGID, &b_MC_ZZ_PDGID);
   fChain->SetBranchAddress("MC_ZZ_D_MASS", &MC_ZZ_D_MASS, &b_MC_ZZ_D_MASS);
   fChain->SetBranchAddress("MC_ZZ_D_PT", &MC_ZZ_D_PT, &b_MC_ZZ_D_PT);
   fChain->SetBranchAddress("MC_ZZ_D_ETA", &MC_ZZ_D_ETA, &b_MC_ZZ_D_ETA);
   fChain->SetBranchAddress("MC_ZZ_D_PHI", &MC_ZZ_D_PHI, &b_MC_ZZ_D_PHI);
   fChain->SetBranchAddress("MC_ZZ_D_THETA", &MC_ZZ_D_THETA, &b_MC_ZZ_D_THETA);
   fChain->SetBranchAddress("MC_ZZ_D_PDGID", &MC_ZZ_D_PDGID, &b_MC_ZZ_D_PDGID);
   fChain->SetBranchAddress("MC_ZZ_DD_MASS", &MC_ZZ_DD_MASS, &b_MC_ZZ_DD_MASS);
   fChain->SetBranchAddress("MC_ZZ_DD_PT", &MC_ZZ_DD_PT, &b_MC_ZZ_DD_PT);
   fChain->SetBranchAddress("MC_ZZ_DD_ETA", &MC_ZZ_DD_ETA, &b_MC_ZZ_DD_ETA);
   fChain->SetBranchAddress("MC_ZZ_DD_PHI", &MC_ZZ_DD_PHI, &b_MC_ZZ_DD_PHI);
   fChain->SetBranchAddress("MC_ZZ_DD_THETA", &MC_ZZ_DD_THETA, &b_MC_ZZ_DD_THETA);
   fChain->SetBranchAddress("MC_ZZ_DD_PDGID", &MC_ZZ_DD_PDGID, &b_MC_ZZ_DD_PDGID);
   fChain->SetBranchAddress("MC_GENJET_PT", &MC_GENJET_PT, &b_MC_GENJET_PT);
   fChain->SetBranchAddress("MC_GENJET_ETA", &MC_GENJET_ETA, &b_MC_GENJET_ETA);
   fChain->SetBranchAddress("MC_GENJET_PHI", &MC_GENJET_PHI, &b_MC_GENJET_PHI);
   fChain->SetBranchAddress("MC_GENMET", &MC_GENMET, &b_MC_GENMET);
   fChain->SetBranchAddress("RECORF_2e2mu_cosTheta1_spin", &RECORF_2e2mu_cosTheta1_spin, &b_RECORF_2e2mu_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_cosTheta2_spin", &RECORF_2e2mu_cosTheta2_spin, &b_RECORF_2e2mu_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_cosThetaStar_spin", &RECORF_2e2mu_cosThetaStar_spin, &b_RECORF_2e2mu_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi_spin", &RECORF_2e2mu_Phi_spin, &b_RECORF_2e2mu_Phi_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi1_spin", &RECORF_2e2mu_Phi1_spin, &b_RECORF_2e2mu_Phi1_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi2_spin", &RECORF_2e2mu_Phi2_spin, &b_RECORF_2e2mu_Phi2_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_phi1RF_spin", &RECORF_2e2mu_phi1RF_spin, &b_RECORF_2e2mu_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_phi2RF_spin", &RECORF_2e2mu_phi2RF_spin, &b_RECORF_2e2mu_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_MELA", &RECORF_2e2mu_MELA, &b_RECORF_2e2mu_MELA);
   fChain->SetBranchAddress("RECORF_4e_cosTheta1_spin", &RECORF_4e_cosTheta1_spin, &b_RECORF_4e_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_4e_cosTheta2_spin", &RECORF_4e_cosTheta2_spin, &b_RECORF_4e_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_4e_cosThetaStar_spin", &RECORF_4e_cosThetaStar_spin, &b_RECORF_4e_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi_spin", &RECORF_4e_Phi_spin, &b_RECORF_4e_Phi_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi1_spin", &RECORF_4e_Phi1_spin, &b_RECORF_4e_Phi1_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi2_spin", &RECORF_4e_Phi2_spin, &b_RECORF_4e_Phi2_spin);
   fChain->SetBranchAddress("RECORF_4e_phi1RF_spin", &RECORF_4e_phi1RF_spin, &b_RECORF_4e_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_4e_phi2RF_spin", &RECORF_4e_phi2RF_spin, &b_RECORF_4e_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_4e_MELA", &RECORF_4e_MELA, &b_RECORF_4e_MELA);
   fChain->SetBranchAddress("RECORF_4mu_cosTheta1_spin", &RECORF_4mu_cosTheta1_spin, &b_RECORF_4mu_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_4mu_cosTheta2_spin", &RECORF_4mu_cosTheta2_spin, &b_RECORF_4mu_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_4mu_cosThetaStar_spin", &RECORF_4mu_cosThetaStar_spin, &b_RECORF_4mu_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi_spin", &RECORF_4mu_Phi_spin, &b_RECORF_4mu_Phi_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi1_spin", &RECORF_4mu_Phi1_spin, &b_RECORF_4mu_Phi1_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi2_spin", &RECORF_4mu_Phi2_spin, &b_RECORF_4mu_Phi2_spin);
   fChain->SetBranchAddress("RECORF_4mu_phi1RF_spin", &RECORF_4mu_phi1RF_spin, &b_RECORF_4mu_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_4mu_phi2RF_spin", &RECORF_4mu_phi2RF_spin, &b_RECORF_4mu_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_4mu_MELA", &RECORF_4mu_MELA, &b_RECORF_4mu_MELA);
   fChain->SetBranchAddress("RECO_ZMM_MASS", &RECO_ZMM_MASS, &b_RECO_ZMM_MASS);
   fChain->SetBranchAddress("RECO_ZEE_MASS", &RECO_ZEE_MASS, &b_RECO_ZEE_MASS);
   fChain->SetBranchAddress("RECO_DiLep_MASS", &RECO_DiLep_MASS, &b_RECO_DiLep_MASS);
   fChain->SetBranchAddress("RECO_ZMM_PT", &RECO_ZMM_PT, &b_RECO_ZMM_PT);
   fChain->SetBranchAddress("RECO_ZEE_PT", &RECO_ZEE_PT, &b_RECO_ZEE_PT);
   fChain->SetBranchAddress("RECO_DiLep_PT", &RECO_DiLep_PT, &b_RECO_DiLep_PT);
   fChain->SetBranchAddress("RECO_ZMM_ETA", &RECO_ZMM_ETA, &b_RECO_ZMM_ETA);
   fChain->SetBranchAddress("RECO_ZEE_ETA", &RECO_ZEE_ETA, &b_RECO_ZEE_ETA);
   fChain->SetBranchAddress("RECO_DiLep_ETA", &RECO_DiLep_ETA, &b_RECO_DiLep_ETA);
   fChain->SetBranchAddress("RECO_ZMM_PHI", &RECO_ZMM_PHI, &b_RECO_ZMM_PHI);
   fChain->SetBranchAddress("RECO_ZEE_PHI", &RECO_ZEE_PHI, &b_RECO_ZEE_PHI);
   fChain->SetBranchAddress("RECO_DiLep_PHI", &RECO_DiLep_PHI, &b_RECO_DiLep_PHI);
   fChain->SetBranchAddress("RECO_ZMMss_MASS", &RECO_ZMMss_MASS, &b_RECO_ZMMss_MASS);
   fChain->SetBranchAddress("RECO_ZEEss_MASS", &RECO_ZEEss_MASS, &b_RECO_ZEEss_MASS);
   fChain->SetBranchAddress("RECO_ZEM_MASS", &RECO_ZEM_MASS, &b_RECO_ZEM_MASS);
   fChain->SetBranchAddress("RECO_ZMMss_PT", &RECO_ZMMss_PT, &b_RECO_ZMMss_PT);
   fChain->SetBranchAddress("RECO_ZEEss_PT", &RECO_ZEEss_PT, &b_RECO_ZEEss_PT);
   fChain->SetBranchAddress("RECO_ZEM_PT", &RECO_ZEM_PT, &b_RECO_ZEM_PT);
   fChain->SetBranchAddress("RECO_ZMMss_ETA", &RECO_ZMMss_ETA, &b_RECO_ZMMss_ETA);
   fChain->SetBranchAddress("RECO_ZEEss_ETA", &RECO_ZEEss_ETA, &b_RECO_ZEEss_ETA);
   fChain->SetBranchAddress("RECO_ZEM_ETA", &RECO_ZEM_ETA, &b_RECO_ZEM_ETA);
   fChain->SetBranchAddress("RECO_ZMMss_PHI", &RECO_ZMMss_PHI, &b_RECO_ZMMss_PHI);
   fChain->SetBranchAddress("RECO_ZEEss_PHI", &RECO_ZEEss_PHI, &b_RECO_ZEEss_PHI);
   fChain->SetBranchAddress("RECO_ZEM_PHI", &RECO_ZEM_PHI, &b_RECO_ZEM_PHI);
   fChain->SetBranchAddress("RECO_MMMM_MASS", &RECO_MMMM_MASS, &b_RECO_MMMM_MASS);
   fChain->SetBranchAddress("RECO_MMMM_PT", &RECO_MMMM_PT, &b_RECO_MMMM_PT);
   fChain->SetBranchAddress("RECO_MMMM_ETA", &RECO_MMMM_ETA, &b_RECO_MMMM_ETA);
   fChain->SetBranchAddress("RECO_MMMM_PHI", &RECO_MMMM_PHI, &b_RECO_MMMM_PHI);
   fChain->SetBranchAddress("RECO_MMMM_MASS_REFIT", &RECO_MMMM_MASS_REFIT, &b_RECO_MMMM_MASS_REFIT);
   fChain->SetBranchAddress("RECO_EEEE_MASS", &RECO_EEEE_MASS, &b_RECO_EEEE_MASS);
   fChain->SetBranchAddress("RECO_EEEE_PT", &RECO_EEEE_PT, &b_RECO_EEEE_PT);
   fChain->SetBranchAddress("RECO_EEEE_ETA", &RECO_EEEE_ETA, &b_RECO_EEEE_ETA);
   fChain->SetBranchAddress("RECO_EEEE_PHI", &RECO_EEEE_PHI, &b_RECO_EEEE_PHI);
   fChain->SetBranchAddress("RECO_EEEE_MASS_REFIT", &RECO_EEEE_MASS_REFIT, &b_RECO_EEEE_MASS_REFIT);
   fChain->SetBranchAddress("RECO_EEMM_MASS", &RECO_EEMM_MASS, &b_RECO_EEMM_MASS);
   fChain->SetBranchAddress("RECO_EEMM_PT", &RECO_EEMM_PT, &b_RECO_EEMM_PT);
   fChain->SetBranchAddress("RECO_EEMM_ETA", &RECO_EEMM_ETA, &b_RECO_EEMM_ETA);
   fChain->SetBranchAddress("RECO_EEMM_PHI", &RECO_EEMM_PHI, &b_RECO_EEMM_PHI);
   fChain->SetBranchAddress("RECO_EEMM_MASS_REFIT", &RECO_EEMM_MASS_REFIT, &b_RECO_EEMM_MASS_REFIT);
   fChain->SetBranchAddress("RECO_LLL0_MASS", &RECO_LLL0_MASS, &b_RECO_LLL0_MASS);
   fChain->SetBranchAddress("RECO_LLL1_MASS", &RECO_LLL1_MASS, &b_RECO_LLL1_MASS);
   fChain->SetBranchAddress("RECO_LLL2_MASS", &RECO_LLL2_MASS, &b_RECO_LLL2_MASS);
   fChain->SetBranchAddress("RECO_LLL3_MASS", &RECO_LLL3_MASS, &b_RECO_LLL3_MASS);
   fChain->SetBranchAddress("RECO_LLL0_PT", &RECO_LLL0_PT, &b_RECO_LLL0_PT);
   fChain->SetBranchAddress("RECO_LLL1_PT", &RECO_LLL1_PT, &b_RECO_LLL1_PT);
   fChain->SetBranchAddress("RECO_LLL2_PT", &RECO_LLL2_PT, &b_RECO_LLL2_PT);
   fChain->SetBranchAddress("RECO_LLL3_PT", &RECO_LLL3_PT, &b_RECO_LLL3_PT);
   fChain->SetBranchAddress("RECO_LLLl0_MASS", &RECO_LLLl0_MASS, &b_RECO_LLLl0_MASS);
   fChain->SetBranchAddress("RECO_LLLl1_MASS", &RECO_LLLl1_MASS, &b_RECO_LLLl1_MASS);
   fChain->SetBranchAddress("RECO_LLLl0_PT", &RECO_LLLl0_PT, &b_RECO_LLLl0_PT);
   fChain->SetBranchAddress("RECO_LLLl1_PT", &RECO_LLLl1_PT, &b_RECO_LLLl1_PT);
   fChain->SetBranchAddress("RECO_LLLL0ss_MASS", &RECO_LLLL0ss_MASS, &b_RECO_LLLL0ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL0ss_PT", &RECO_LLLL0ss_PT, &b_RECO_LLLL0ss_PT);
   fChain->SetBranchAddress("RECO_LLLL1ss_MASS", &RECO_LLLL1ss_MASS, &b_RECO_LLLL1ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL1ss_PT", &RECO_LLLL1ss_PT, &b_RECO_LLLL1ss_PT);
   fChain->SetBranchAddress("RECO_LLLL2ss_MASS", &RECO_LLLL2ss_MASS, &b_RECO_LLLL2ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL2ss_PT", &RECO_LLLL2ss_PT, &b_RECO_LLLL2ss_PT);
   fChain->SetBranchAddress("RECO_LLLL_MASS", &RECO_LLLL_MASS, &b_RECO_LLLL_MASS);
   fChain->SetBranchAddress("RECO_LLLL_PT", &RECO_LLLL_PT, &b_RECO_LLLL_PT);
   fChain->SetBranchAddress("RECO_LLLL_ETA", &RECO_LLLL_ETA, &b_RECO_LLLL_ETA);
   fChain->SetBranchAddress("RECO_LLLL_PHI", &RECO_LLLL_PHI, &b_RECO_LLLL_PHI);
   fChain->SetBranchAddress("RECOELE_E", &RECOELE_E, &b_RECOELE_E);
   fChain->SetBranchAddress("RECOELE_PT", &RECOELE_PT, &b_RECOELE_PT);
   fChain->SetBranchAddress("RECOELE_PTError", &RECOELE_PTError, &b_RECOELE_PTError);
   fChain->SetBranchAddress("RECOELE_P", &RECOELE_P, &b_RECOELE_P);
   fChain->SetBranchAddress("RECOELE_ETA", &RECOELE_ETA, &b_RECOELE_ETA);
   fChain->SetBranchAddress("RECOELE_THETA", &RECOELE_THETA, &b_RECOELE_THETA);
   fChain->SetBranchAddress("RECOELE_PHI", &RECOELE_PHI, &b_RECOELE_PHI);
   fChain->SetBranchAddress("RECOELE_MASS", &RECOELE_MASS, &b_RECOELE_MASS);
   fChain->SetBranchAddress("RECOELE_CHARGE", &RECOELE_CHARGE, &b_RECOELE_CHARGE);
   fChain->SetBranchAddress("RECOELE_ID", &RECOELE_ID, &b_RECOELE_ID);
   fChain->SetBranchAddress("RECOELE_PT_uncorr", &RECOELE_PT_uncorr, &b_RECOELE_PT_uncorr);
   fChain->SetBranchAddress("RECOELE_isEcalDriven", &RECOELE_isEcalDriven, &b_RECOELE_isEcalDriven);
   fChain->SetBranchAddress("RECOELE_isTrackerDriven", &RECOELE_isTrackerDriven, &b_RECOELE_isTrackerDriven);
   fChain->SetBranchAddress("RECOELE_gsftrack_NPixHits", &RECOELE_gsftrack_NPixHits, &b_RECOELE_gsftrack_NPixHits);
   fChain->SetBranchAddress("RECOELE_gsftrack_NStripHits", &RECOELE_gsftrack_NStripHits, &b_RECOELE_gsftrack_NStripHits);
   fChain->SetBranchAddress("RECOELE_gsftrack_chi2", &RECOELE_gsftrack_chi2, &b_RECOELE_gsftrack_chi2);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxyB", &RECOELE_gsftrack_dxyB, &b_RECOELE_gsftrack_dxyB);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxy", &RECOELE_gsftrack_dxy, &b_RECOELE_gsftrack_dxy);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxyError", &RECOELE_gsftrack_dxyError, &b_RECOELE_gsftrack_dxyError);
   fChain->SetBranchAddress("RECOELE_gsftrack_dzB", &RECOELE_gsftrack_dzB, &b_RECOELE_gsftrack_dzB);
   fChain->SetBranchAddress("RECOELE_gsftrack_dz", &RECOELE_gsftrack_dz, &b_RECOELE_gsftrack_dz);
   fChain->SetBranchAddress("RECOELE_gsftrack_dzError", &RECOELE_gsftrack_dzError, &b_RECOELE_gsftrack_dzError);
   fChain->SetBranchAddress("RECOELE_gsftrack_losthits", &RECOELE_gsftrack_losthits, &b_RECOELE_gsftrack_losthits);
   fChain->SetBranchAddress("RECOELE_gsftrack_validhits", &RECOELE_gsftrack_validhits, &b_RECOELE_gsftrack_validhits);
   fChain->SetBranchAddress("RECOELE_gsftrack_expected_inner_hits", &RECOELE_gsftrack_expected_inner_hits, &b_RECOELE_gsftrack_expected_inner_hits);
   fChain->SetBranchAddress("RECOELE_scl_E", &RECOELE_scl_E, &b_RECOELE_scl_E);
   fChain->SetBranchAddress("RECOELE_scl_Et", &RECOELE_scl_Et, &b_RECOELE_scl_Et);
   fChain->SetBranchAddress("RECOELE_scl_Eta", &RECOELE_scl_Eta, &b_RECOELE_scl_Eta);
   fChain->SetBranchAddress("RECOELE_scl_Phi", &RECOELE_scl_Phi, &b_RECOELE_scl_Phi);
   fChain->SetBranchAddress("RECOELE_ecalEnergy", &RECOELE_ecalEnergy, &b_RECOELE_ecalEnergy);
   fChain->SetBranchAddress("RECOELE_ep", &RECOELE_ep, &b_RECOELE_ep);
   fChain->SetBranchAddress("RECOELE_eSeedp", &RECOELE_eSeedp, &b_RECOELE_eSeedp);
   fChain->SetBranchAddress("RECOELE_eSeedpout", &RECOELE_eSeedpout, &b_RECOELE_eSeedpout);
   fChain->SetBranchAddress("RECOELE_eElepout", &RECOELE_eElepout, &b_RECOELE_eElepout);
   fChain->SetBranchAddress("RECOELE_deltaEtaIn", &RECOELE_deltaEtaIn, &b_RECOELE_deltaEtaIn);
   fChain->SetBranchAddress("RECOELE_deltaEtaSeed", &RECOELE_deltaEtaSeed, &b_RECOELE_deltaEtaSeed);
   fChain->SetBranchAddress("RECOELE_deltaEtaEle", &RECOELE_deltaEtaEle, &b_RECOELE_deltaEtaEle);
   fChain->SetBranchAddress("RECOELE_deltaPhiIn", &RECOELE_deltaPhiIn, &b_RECOELE_deltaPhiIn);
   fChain->SetBranchAddress("RECOELE_deltaPhiSeed", &RECOELE_deltaPhiSeed, &b_RECOELE_deltaPhiSeed);
   fChain->SetBranchAddress("RECOELE_deltaPhiEle", &RECOELE_deltaPhiEle, &b_RECOELE_deltaPhiEle);
   fChain->SetBranchAddress("RECOELE_isbarrel", &RECOELE_isbarrel, &b_RECOELE_isbarrel);
   fChain->SetBranchAddress("RECOELE_isendcap", &RECOELE_isendcap, &b_RECOELE_isendcap);
   fChain->SetBranchAddress("RECOELE_isGap", &RECOELE_isGap, &b_RECOELE_isGap);
   fChain->SetBranchAddress("RECOELE_isEBetaGap", &RECOELE_isEBetaGap, &b_RECOELE_isEBetaGap);
   fChain->SetBranchAddress("RECOELE_isEBphiGap", &RECOELE_isEBphiGap, &b_RECOELE_isEBphiGap);
   fChain->SetBranchAddress("RECOELE_isEEdeeGap", &RECOELE_isEEdeeGap, &b_RECOELE_isEEdeeGap);
   fChain->SetBranchAddress("RECOELE_isEEringGap", &RECOELE_isEEringGap, &b_RECOELE_isEEringGap);
   fChain->SetBranchAddress("RECOELE_sigmaIetaIeta", &RECOELE_sigmaIetaIeta, &b_RECOELE_sigmaIetaIeta);
   fChain->SetBranchAddress("RECOELE_sigmaEtaEta", &RECOELE_sigmaEtaEta, &b_RECOELE_sigmaEtaEta);
   fChain->SetBranchAddress("RECOELE_e15", &RECOELE_e15, &b_RECOELE_e15);
   fChain->SetBranchAddress("RECOELE_e25max", &RECOELE_e25max, &b_RECOELE_e25max);
   fChain->SetBranchAddress("RECOELE_e55", &RECOELE_e55, &b_RECOELE_e55);
   fChain->SetBranchAddress("RECOELE_he", &RECOELE_he, &b_RECOELE_he);
   fChain->SetBranchAddress("RECOELE_r9", &RECOELE_r9, &b_RECOELE_r9);
   fChain->SetBranchAddress("RECOELE_mva", &RECOELE_mva, &b_RECOELE_mva);
   fChain->SetBranchAddress("RECOELE_fbrem", &RECOELE_fbrem, &b_RECOELE_fbrem);
   fChain->SetBranchAddress("RECOELE_nbrems", &RECOELE_nbrems, &b_RECOELE_nbrems);
   fChain->SetBranchAddress("RECOELE_Class", &RECOELE_Class, &b_RECOELE_Class);
   fChain->SetBranchAddress("RECOELE_fbrem_mode", &RECOELE_fbrem_mode, &b_RECOELE_fbrem_mode);
   fChain->SetBranchAddress("RECOELE_fbrem_mean", &RECOELE_fbrem_mean, &b_RECOELE_fbrem_mean);
   fChain->SetBranchAddress("RECOELE_EGMTRACKISO", &RECOELE_EGMTRACKISO, &b_RECOELE_EGMTRACKISO);
   fChain->SetBranchAddress("RECOELE_EGMHCALISO", &RECOELE_EGMHCALISO, &b_RECOELE_EGMHCALISO);
   fChain->SetBranchAddress("RECOELE_EGMECALISO", &RECOELE_EGMECALISO, &b_RECOELE_EGMECALISO);
   fChain->SetBranchAddress("RECOELE_EGMX", &RECOELE_EGMX, &b_RECOELE_EGMX);
   fChain->SetBranchAddress("RECOELE_PFchAllPart", &RECOELE_PFchAllPart, &b_RECOELE_PFchAllPart);
   fChain->SetBranchAddress("RECOELE_PFchHad", &RECOELE_PFchHad, &b_RECOELE_PFchHad);
   fChain->SetBranchAddress("RECOELE_PFneuHad", &RECOELE_PFneuHad, &b_RECOELE_PFneuHad);
   fChain->SetBranchAddress("RECOELE_PFphoton", &RECOELE_PFphoton, &b_RECOELE_PFphoton);
   fChain->SetBranchAddress("RECOELE_PFPUchAllPart", &RECOELE_PFPUchAllPart, &b_RECOELE_PFPUchAllPart);
   fChain->SetBranchAddress("RECOELE_PFX_dB", &RECOELE_PFX_dB, &b_RECOELE_PFX_dB);
   fChain->SetBranchAddress("RECOELE_PFX_rho", &RECOELE_PFX_rho, &b_RECOELE_PFX_rho);
   fChain->SetBranchAddress("RECOELE_regEnergy", &RECOELE_regEnergy, &b_RECOELE_regEnergy);
   fChain->SetBranchAddress("RECOELE_regEnergyError", &RECOELE_regEnergyError, &b_RECOELE_regEnergyError);
   fChain->SetBranchAddress("RECOELE_SIP", &RECOELE_SIP, &b_RECOELE_SIP);
   fChain->SetBranchAddress("RECOELE_IP", &RECOELE_IP, &b_RECOELE_IP);
   fChain->SetBranchAddress("RECOELE_IPERROR", &RECOELE_IPERROR, &b_RECOELE_IPERROR);
   fChain->SetBranchAddress("RECOELE_SIP_KF", &RECOELE_SIP_KF, &b_RECOELE_SIP_KF);
   fChain->SetBranchAddress("RECOELE_IP_KF", &RECOELE_IP_KF, &b_RECOELE_IP_KF);
   fChain->SetBranchAddress("RECOELE_IPERROR_KF", &RECOELE_IPERROR_KF, &b_RECOELE_IPERROR_KF);
   fChain->SetBranchAddress("RECOELE_SIP_GD", &RECOELE_SIP_GD, &b_RECOELE_SIP_GD);
   fChain->SetBranchAddress("RECOELE_SIP_GDEEEE", &RECOELE_SIP_GDEEEE, &b_RECOELE_SIP_GDEEEE);
   fChain->SetBranchAddress("RECOELE_SIP_Std", &RECOELE_SIP_Std, &b_RECOELE_SIP_Std);
   fChain->SetBranchAddress("RECOELE_SIP_StdEEEE", &RECOELE_SIP_StdEEEE, &b_RECOELE_SIP_StdEEEE);
   fChain->SetBranchAddress("RECOELE_SIP_Kin", &RECOELE_SIP_Kin, &b_RECOELE_SIP_Kin);
   fChain->SetBranchAddress("RECOELE_SIP_KinEEEE", &RECOELE_SIP_KinEEEE, &b_RECOELE_SIP_KinEEEE);
   fChain->SetBranchAddress("RECOELE_STIP", &RECOELE_STIP, &b_RECOELE_STIP);
   fChain->SetBranchAddress("RECOELE_SLIP", &RECOELE_SLIP, &b_RECOELE_SLIP);
   fChain->SetBranchAddress("RECOELE_TIP", &RECOELE_TIP, &b_RECOELE_TIP);
   fChain->SetBranchAddress("RECOELE_LIP", &RECOELE_LIP, &b_RECOELE_LIP);
   fChain->SetBranchAddress("RECOELE_TIPERROR", &RECOELE_TIPERROR, &b_RECOELE_TIPERROR);
   fChain->SetBranchAddress("RECOELE_LIPERROR", &RECOELE_LIPERROR, &b_RECOELE_LIPERROR);
   fChain->SetBranchAddress("RECOELE_sclRawE", &RECOELE_sclRawE, &b_RECOELE_sclRawE);
   fChain->SetBranchAddress("RECOELE_sclX", &RECOELE_sclX, &b_RECOELE_sclX);
   fChain->SetBranchAddress("RECOELE_sclY", &RECOELE_sclY, &b_RECOELE_sclY);
   fChain->SetBranchAddress("RECOELE_sclZ", &RECOELE_sclZ, &b_RECOELE_sclZ);
   fChain->SetBranchAddress("RECOELE_seedSubdet1", &RECOELE_seedSubdet1, &b_RECOELE_seedSubdet1);
   fChain->SetBranchAddress("RECOELE_seedDphi1", &RECOELE_seedDphi1, &b_RECOELE_seedDphi1);
   fChain->SetBranchAddress("RECOELE_seedDrz1", &RECOELE_seedDrz1, &b_RECOELE_seedDrz1);
   fChain->SetBranchAddress("RECOELE_seedSubdet2", &RECOELE_seedSubdet2, &b_RECOELE_seedSubdet2);
   fChain->SetBranchAddress("RECOELE_seedDphi2", &RECOELE_seedDphi2, &b_RECOELE_seedDphi2);
   fChain->SetBranchAddress("RECOELE_seedDrz2", &RECOELE_seedDrz2, &b_RECOELE_seedDrz2);
   fChain->SetBranchAddress("RECOELE_eidVeryLoose", &RECOELE_eidVeryLoose, &b_RECOELE_eidVeryLoose);
   fChain->SetBranchAddress("RECOELE_eidLoose", &RECOELE_eidLoose, &b_RECOELE_eidLoose);
   fChain->SetBranchAddress("RECOELE_eidMedium", &RECOELE_eidMedium, &b_RECOELE_eidMedium);
   fChain->SetBranchAddress("RECOELE_eidTight", &RECOELE_eidTight, &b_RECOELE_eidTight);
   fChain->SetBranchAddress("RECOELE_eidHZZVeryLoose", &RECOELE_eidHZZVeryLoose, &b_RECOELE_eidHZZVeryLoose);
   fChain->SetBranchAddress("RECOELE_eidHZZLoose", &RECOELE_eidHZZLoose, &b_RECOELE_eidHZZLoose);
   fChain->SetBranchAddress("RECOELE_eidHZZMedium", &RECOELE_eidHZZMedium, &b_RECOELE_eidHZZMedium);
   fChain->SetBranchAddress("RECOELE_eidHZZTight", &RECOELE_eidHZZTight, &b_RECOELE_eidHZZTight);
   fChain->SetBranchAddress("RECOELE_mvaTrigV0", &RECOELE_mvaTrigV0, &b_RECOELE_mvaTrigV0);
   fChain->SetBranchAddress("RECOELE_mvaNonTrigV0", &RECOELE_mvaNonTrigV0, &b_RECOELE_mvaNonTrigV0);
   fChain->SetBranchAddress("RECOELE_COV", &RECOELE_COV, &b_RECOELE_COV);
   fChain->SetBranchAddress("RECOELE_TLE_ParentSC_X", &RECOELE_TLE_ParentSC_X, &b_RECOELE_TLE_ParentSC_X);
   fChain->SetBranchAddress("RECOELE_TLE_ParentSC_Y", &RECOELE_TLE_ParentSC_Y, &b_RECOELE_TLE_ParentSC_Y);
   fChain->SetBranchAddress("RECOELE_TLE_ParentSC_Z", &RECOELE_TLE_ParentSC_Z, &b_RECOELE_TLE_ParentSC_Z);
   fChain->SetBranchAddress("RECOELE_ecalTrkEnergyPreCorr", &RECOELE_ecalTrkEnergyPreCorr, &b_RECOELE_ecalTrkEnergyPreCorr);
   fChain->SetBranchAddress("RECOELE_ecalTrkEnergyErrPreCorr", &RECOELE_ecalTrkEnergyErrPreCorr, &b_RECOELE_ecalTrkEnergyErrPreCorr);
   fChain->SetBranchAddress("RECOELE_ecalTrkEnergyErrPostCorr", &RECOELE_ecalTrkEnergyErrPostCorr, &b_RECOELE_ecalTrkEnergyErrPostCorr);
   fChain->SetBranchAddress("RECOELE_energyScaleValue", &RECOELE_energyScaleValue, &b_RECOELE_energyScaleValue);
   fChain->SetBranchAddress("RECOELE_energySigmaValue", &RECOELE_energySigmaValue, &b_RECOELE_energySigmaValue);
   fChain->SetBranchAddress("RECOELE_energyScaleUp", &RECOELE_energyScaleUp, &b_RECOELE_energyScaleUp);
   fChain->SetBranchAddress("RECOELE_energyScaleDown", &RECOELE_energyScaleDown, &b_RECOELE_energyScaleDown);
   fChain->SetBranchAddress("RECOELE_energyScaleStatUp", &RECOELE_energyScaleStatUp, &b_RECOELE_energyScaleStatUp);
   fChain->SetBranchAddress("RECOELE_energyScaleStatDown", &RECOELE_energyScaleStatDown, &b_RECOELE_energyScaleStatDown);
   fChain->SetBranchAddress("RECOELE_energyScaleSystUp", &RECOELE_energyScaleSystUp, &b_RECOELE_energyScaleSystUp);
   fChain->SetBranchAddress("RECOELE_energyScaleSystDown", &RECOELE_energyScaleSystDown, &b_RECOELE_energyScaleSystDown);
   fChain->SetBranchAddress("RECOELE_energyScaleGainUp", &RECOELE_energyScaleGainUp, &b_RECOELE_energyScaleGainUp);
   fChain->SetBranchAddress("RECOELE_energyScaleGainDown", &RECOELE_energyScaleGainDown, &b_RECOELE_energyScaleGainDown);
   fChain->SetBranchAddress("RECOELE_energyScaleEtUp", &RECOELE_energyScaleEtUp, &b_RECOELE_energyScaleEtUp);
   fChain->SetBranchAddress("RECOELE_energyScaleEtDown", &RECOELE_energyScaleEtDown, &b_RECOELE_energyScaleEtDown);
   fChain->SetBranchAddress("RECOELE_energySigmaUp", &RECOELE_energySigmaUp, &b_RECOELE_energySigmaUp);
   fChain->SetBranchAddress("RECOELE_energySigmaDown", &RECOELE_energySigmaDown, &b_RECOELE_energySigmaDown);
   fChain->SetBranchAddress("RECOELE_energySigmaPhiUp", &RECOELE_energySigmaPhiUp, &b_RECOELE_energySigmaPhiUp);
   fChain->SetBranchAddress("RECOELE_energySigmaPhiDown", &RECOELE_energySigmaPhiDown, &b_RECOELE_energySigmaPhiDown);
   fChain->SetBranchAddress("RECOELE_energySigmaRhoUp", &RECOELE_energySigmaRhoUp, &b_RECOELE_energySigmaRhoUp);
   fChain->SetBranchAddress("RECOELE_energySigmaRhoDown", &RECOELE_energySigmaRhoDown, &b_RECOELE_energySigmaRhoDown);
   fChain->SetBranchAddress("RECOMU_isPFMu", &RECOMU_isPFMu, &b_RECOMU_isPFMu);
   fChain->SetBranchAddress("RECOMU_isGlobalMu", &RECOMU_isGlobalMu, &b_RECOMU_isGlobalMu);
   fChain->SetBranchAddress("RECOMU_isStandAloneMu", &RECOMU_isStandAloneMu, &b_RECOMU_isStandAloneMu);
   fChain->SetBranchAddress("RECOMU_isTrackerMu", &RECOMU_isTrackerMu, &b_RECOMU_isTrackerMu);
   fChain->SetBranchAddress("RECOMU_isCaloMu", &RECOMU_isCaloMu, &b_RECOMU_isCaloMu);
   fChain->SetBranchAddress("RECOMU_isTrackerHighPtMu", &RECOMU_isTrackerHighPtMu, &b_RECOMU_isTrackerHighPtMu);
   fChain->SetBranchAddress("RECOMU_isME0Muon", &RECOMU_isME0Muon, &b_RECOMU_isME0Muon);
   fChain->SetBranchAddress("RECOMU_E", &RECOMU_E, &b_RECOMU_E);
   fChain->SetBranchAddress("RECOMU_PT", &RECOMU_PT, &b_RECOMU_PT);
   fChain->SetBranchAddress("RECOMU_P", &RECOMU_P, &b_RECOMU_P);
   fChain->SetBranchAddress("RECOMU_ETA", &RECOMU_ETA, &b_RECOMU_ETA);
   fChain->SetBranchAddress("RECOMU_THETA", &RECOMU_THETA, &b_RECOMU_THETA);
   fChain->SetBranchAddress("RECOMU_PHI", &RECOMU_PHI, &b_RECOMU_PHI);
   fChain->SetBranchAddress("RECOMU_MASS", &RECOMU_MASS, &b_RECOMU_MASS);
   fChain->SetBranchAddress("RECOMU_CHARGE", &RECOMU_CHARGE, &b_RECOMU_CHARGE);
   fChain->SetBranchAddress("RECOMU_PT_uncorr", &RECOMU_PT_uncorr, &b_RECOMU_PT_uncorr);
   fChain->SetBranchAddress("RECOMU_COV", &RECOMU_COV, &b_RECOMU_COV);
   fChain->SetBranchAddress("RECOMU_TRACKISO", &RECOMU_TRACKISO, &b_RECOMU_TRACKISO);
   fChain->SetBranchAddress("RECOMU_TRACKISO_SUMPT", &RECOMU_TRACKISO_SUMPT, &b_RECOMU_TRACKISO_SUMPT);
   fChain->SetBranchAddress("RECOMU_HCALISO", &RECOMU_HCALISO, &b_RECOMU_HCALISO);
   fChain->SetBranchAddress("RECOMU_ECALISO", &RECOMU_ECALISO, &b_RECOMU_ECALISO);
   fChain->SetBranchAddress("RECOMU_X", &RECOMU_X, &b_RECOMU_X);
   fChain->SetBranchAddress("RECOMU_PFchHad", &RECOMU_PFchHad, &b_RECOMU_PFchHad);
   fChain->SetBranchAddress("RECOMU_PFneuHad", &RECOMU_PFneuHad, &b_RECOMU_PFneuHad);
   fChain->SetBranchAddress("RECOMU_PFphoton", &RECOMU_PFphoton, &b_RECOMU_PFphoton);
   fChain->SetBranchAddress("RECOMU_PFPUchAllPart", &RECOMU_PFPUchAllPart, &b_RECOMU_PFPUchAllPart);
   fChain->SetBranchAddress("RECOMU_PFX_dB", &RECOMU_PFX_dB, &b_RECOMU_PFX_dB);
   fChain->SetBranchAddress("RECOMU_PFX_rho", &RECOMU_PFX_rho, &b_RECOMU_PFX_rho);
   fChain->SetBranchAddress("RECOPFPHOT_PFchHad", &RECOPFPHOT_PFchHad, &b_RECOPFPHOT_PFchHad);
   fChain->SetBranchAddress("RECOPFPHOT_PFneuHad", &RECOPFPHOT_PFneuHad, &b_RECOPFPHOT_PFneuHad);
   fChain->SetBranchAddress("RECOPFPHOT_PFphoton", &RECOPFPHOT_PFphoton, &b_RECOPFPHOT_PFphoton);
   fChain->SetBranchAddress("RECOPFPHOT_PFPUchAllPart", &RECOPFPHOT_PFPUchAllPart, &b_RECOPFPHOT_PFPUchAllPart);
   fChain->SetBranchAddress("RECOPFPHOT_PFX_rho", &RECOPFPHOT_PFX_rho, &b_RECOPFPHOT_PFX_rho);
   fChain->SetBranchAddress("RECOMU_SIP", &RECOMU_SIP, &b_RECOMU_SIP);
   fChain->SetBranchAddress("RECOMU_IP", &RECOMU_IP, &b_RECOMU_IP);
   fChain->SetBranchAddress("RECOMU_IPERROR", &RECOMU_IPERROR, &b_RECOMU_IPERROR);
   fChain->SetBranchAddress("RECOMU_SIP_KF", &RECOMU_SIP_KF, &b_RECOMU_SIP_KF);
   fChain->SetBranchAddress("RECOMU_IP_KF", &RECOMU_IP_KF, &b_RECOMU_IP_KF);
   fChain->SetBranchAddress("RECOMU_IPERROR_KF", &RECOMU_IPERROR_KF, &b_RECOMU_IPERROR_KF);
   fChain->SetBranchAddress("RECOMU_SIP_GD", &RECOMU_SIP_GD, &b_RECOMU_SIP_GD);
   fChain->SetBranchAddress("RECOMU_SIP_GDMMMM", &RECOMU_SIP_GDMMMM, &b_RECOMU_SIP_GDMMMM);
   fChain->SetBranchAddress("RECOMU_SIP_Std", &RECOMU_SIP_Std, &b_RECOMU_SIP_Std);
   fChain->SetBranchAddress("RECOMU_SIP_StdMMMM", &RECOMU_SIP_StdMMMM, &b_RECOMU_SIP_StdMMMM);
   fChain->SetBranchAddress("RECOMU_SIP_Kin", &RECOMU_SIP_Kin, &b_RECOMU_SIP_Kin);
   fChain->SetBranchAddress("RECOMU_SIP_KinMMMM", &RECOMU_SIP_KinMMMM, &b_RECOMU_SIP_KinMMMM);
   fChain->SetBranchAddress("RECOMU_STIP", &RECOMU_STIP, &b_RECOMU_STIP);
   fChain->SetBranchAddress("RECOMU_SLIP", &RECOMU_SLIP, &b_RECOMU_SLIP);
   fChain->SetBranchAddress("RECOMU_TIP", &RECOMU_TIP, &b_RECOMU_TIP);
   fChain->SetBranchAddress("RECOMU_LIP", &RECOMU_LIP, &b_RECOMU_LIP);
   fChain->SetBranchAddress("RECOMU_TIPERROR", &RECOMU_TIPERROR, &b_RECOMU_TIPERROR);
   fChain->SetBranchAddress("RECOMU_LIPERROR", &RECOMU_LIPERROR, &b_RECOMU_LIPERROR);
   fChain->SetBranchAddress("RECOMU_caloCompatibility", &RECOMU_caloCompatibility, &b_RECOMU_caloCompatibility);
   fChain->SetBranchAddress("RECOMU_segmentCompatibility", &RECOMU_segmentCompatibility, &b_RECOMU_segmentCompatibility);
   fChain->SetBranchAddress("RECOMU_numberOfMatches", &RECOMU_numberOfMatches, &b_RECOMU_numberOfMatches);
   fChain->SetBranchAddress("RECOMU_numberOfMatchedStations", &RECOMU_numberOfMatchedStations, &b_RECOMU_numberOfMatchedStations);
   fChain->SetBranchAddress("RECOMU_glbmuPromptTight", &RECOMU_glbmuPromptTight, &b_RECOMU_glbmuPromptTight);
   fChain->SetBranchAddress("RECOMU_trkmuArbitration", &RECOMU_trkmuArbitration, &b_RECOMU_trkmuArbitration);
   fChain->SetBranchAddress("RECOMU_trkmu2DCompatibilityLoose", &RECOMU_trkmu2DCompatibilityLoose, &b_RECOMU_trkmu2DCompatibilityLoose);
   fChain->SetBranchAddress("RECOMU_trkmu2DCompatibilityTight", &RECOMU_trkmu2DCompatibilityTight, &b_RECOMU_trkmu2DCompatibilityTight);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationLoose", &RECOMU_trkmuOneStationLoose, &b_RECOMU_trkmuOneStationLoose);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationTight", &RECOMU_trkmuOneStationTight, &b_RECOMU_trkmuOneStationTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationLoose", &RECOMU_trkmuLastStationLoose, &b_RECOMU_trkmuLastStationLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationTight", &RECOMU_trkmuLastStationTight, &b_RECOMU_trkmuLastStationTight);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationAngLoose", &RECOMU_trkmuOneStationAngLoose, &b_RECOMU_trkmuOneStationAngLoose);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationAngTight", &RECOMU_trkmuOneStationAngTight, &b_RECOMU_trkmuOneStationAngTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationAngLoose", &RECOMU_trkmuLastStationAngLoose, &b_RECOMU_trkmuLastStationAngLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationAngTight", &RECOMU_trkmuLastStationAngTight, &b_RECOMU_trkmuLastStationAngTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationOptimizedLowPtLoose", &RECOMU_trkmuLastStationOptimizedLowPtLoose, &b_RECOMU_trkmuLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationOptimizedLowPtTight", &RECOMU_trkmuLastStationOptimizedLowPtTight, &b_RECOMU_trkmuLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("RECOMU_mutrkPT", &RECOMU_mutrkPT, &b_RECOMU_mutrkPT);
   fChain->SetBranchAddress("RECOMU_mutrkPTError", &RECOMU_mutrkPTError, &b_RECOMU_mutrkPTError);
   fChain->SetBranchAddress("RECOMU_mutrkDxy", &RECOMU_mutrkDxy, &b_RECOMU_mutrkDxy);
   fChain->SetBranchAddress("RECOMU_mutrkDxyError", &RECOMU_mutrkDxyError, &b_RECOMU_mutrkDxyError);
   fChain->SetBranchAddress("RECOMU_mutrkDxyB", &RECOMU_mutrkDxyB, &b_RECOMU_mutrkDxyB);
   fChain->SetBranchAddress("RECOMU_mutrkDz", &RECOMU_mutrkDz, &b_RECOMU_mutrkDz);
   fChain->SetBranchAddress("RECOMU_mutrkDzError", &RECOMU_mutrkDzError, &b_RECOMU_mutrkDzError);
   fChain->SetBranchAddress("RECOMU_mutrkDzB", &RECOMU_mutrkDzB, &b_RECOMU_mutrkDzB);
   fChain->SetBranchAddress("RECOMU_mutrkChi2PerNdof", &RECOMU_mutrkChi2PerNdof, &b_RECOMU_mutrkChi2PerNdof);
   fChain->SetBranchAddress("RECOMU_mutrkCharge", &RECOMU_mutrkCharge, &b_RECOMU_mutrkCharge);
   fChain->SetBranchAddress("RECOMU_mutrkNHits", &RECOMU_mutrkNHits, &b_RECOMU_mutrkNHits);
   fChain->SetBranchAddress("RECOMU_mutrkNStripHits", &RECOMU_mutrkNStripHits, &b_RECOMU_mutrkNStripHits);
   fChain->SetBranchAddress("RECOMU_mutrkNPixHits", &RECOMU_mutrkNPixHits, &b_RECOMU_mutrkNPixHits);
   fChain->SetBranchAddress("RECOMU_mutrkNMuonHits", &RECOMU_mutrkNMuonHits, &b_RECOMU_mutrkNMuonHits);
   fChain->SetBranchAddress("RECOMU_mutrktrackerLayersWithMeasurement", &RECOMU_mutrktrackerLayersWithMeasurement, &b_RECOMU_mutrktrackerLayersWithMeasurement);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxy", &RECOMU_muInnertrkDxy, &b_RECOMU_muInnertrkDxy);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxyError", &RECOMU_muInnertrkDxyError, &b_RECOMU_muInnertrkDxyError);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxyB", &RECOMU_muInnertrkDxyB, &b_RECOMU_muInnertrkDxyB);
   fChain->SetBranchAddress("RECOMU_muInnertrkDz", &RECOMU_muInnertrkDz, &b_RECOMU_muInnertrkDz);
   fChain->SetBranchAddress("RECOMU_muInnertrkDzError", &RECOMU_muInnertrkDzError, &b_RECOMU_muInnertrkDzError);
   fChain->SetBranchAddress("RECOMU_muInnertrkDzB", &RECOMU_muInnertrkDzB, &b_RECOMU_muInnertrkDzB);
   fChain->SetBranchAddress("RECOMU_muInnertrkChi2PerNdof", &RECOMU_muInnertrkChi2PerNdof, &b_RECOMU_muInnertrkChi2PerNdof);
   fChain->SetBranchAddress("RECOMU_muInnertrktrackerLayersWithMeasurement", &RECOMU_muInnertrktrackerLayersWithMeasurement, &b_RECOMU_muInnertrktrackerLayersWithMeasurement);
   fChain->SetBranchAddress("RECOMU_muInnertrkPT", &RECOMU_muInnertrkPT, &b_RECOMU_muInnertrkPT);
   fChain->SetBranchAddress("RECOMU_muInnertrkPTError", &RECOMU_muInnertrkPTError, &b_RECOMU_muInnertrkPTError);
   fChain->SetBranchAddress("RECOMU_muInnertrkCharge", &RECOMU_muInnertrkCharge, &b_RECOMU_muInnertrkCharge);
   fChain->SetBranchAddress("RECOMU_muInnertrkNHits", &RECOMU_muInnertrkNHits, &b_RECOMU_muInnertrkNHits);
   fChain->SetBranchAddress("RECOMU_muInnertrkNStripHits", &RECOMU_muInnertrkNStripHits, &b_RECOMU_muInnertrkNStripHits);
   fChain->SetBranchAddress("RECOMU_muInnertrkNPixHits", &RECOMU_muInnertrkNPixHits, &b_RECOMU_muInnertrkNPixHits);
   fChain->SetBranchAddress("RECOMU_mubesttrkType", &RECOMU_mubesttrkType, &b_RECOMU_mubesttrkType);
   fChain->SetBranchAddress("RECOMU_mubesttrkDxy", &RECOMU_mubesttrkDxy, &b_RECOMU_mubesttrkDxy);
   fChain->SetBranchAddress("RECOMU_mubesttrkDxyError", &RECOMU_mubesttrkDxyError, &b_RECOMU_mubesttrkDxyError);
   fChain->SetBranchAddress("RECOMU_mubesttrkDz", &RECOMU_mubesttrkDz, &b_RECOMU_mubesttrkDz);
   fChain->SetBranchAddress("RECOMU_mubesttrkDzError", &RECOMU_mubesttrkDzError, &b_RECOMU_mubesttrkDzError);
   fChain->SetBranchAddress("RECOMU_mubesttrkPTError", &RECOMU_mubesttrkPTError, &b_RECOMU_mubesttrkPTError);
   fChain->SetBranchAddress("RECOMU_Rochester_Error", &RECOMU_Rochester_Error, &b_RECOMU_Rochester_Error);
   fChain->SetBranchAddress("ftsigma", &ftsigma, &b_ftsigma);
   fChain->SetBranchAddress("gdX", &gdX, &b_gdX);
   fChain->SetBranchAddress("gdY", &gdY, &b_gdY);
   fChain->SetBranchAddress("gdZ", &gdZ, &b_gdZ);
   fChain->SetBranchAddress("ftsigmalag", &ftsigmalag, &b_ftsigmalag);
   fChain->SetBranchAddress("gdlagX", &gdlagX, &b_gdlagX);
   fChain->SetBranchAddress("gdlagY", &gdlagY, &b_gdlagY);
   fChain->SetBranchAddress("gdlagZ", &gdlagZ, &b_gdlagZ);
   fChain->SetBranchAddress("gdlagProb", &gdlagProb, &b_gdlagProb);
   fChain->SetBranchAddress("gdlagNdof", &gdlagNdof, &b_gdlagNdof);
   fChain->SetBranchAddress("ftsigmaMMMM", &ftsigmaMMMM, &b_ftsigmaMMMM);
   fChain->SetBranchAddress("gdXMMMM", &gdXMMMM, &b_gdXMMMM);
   fChain->SetBranchAddress("gdYMMMM", &gdYMMMM, &b_gdYMMMM);
   fChain->SetBranchAddress("gdZMMMM", &gdZMMMM, &b_gdZMMMM);
   fChain->SetBranchAddress("ftsigmalagMMMM", &ftsigmalagMMMM, &b_ftsigmalagMMMM);
   fChain->SetBranchAddress("gdlagXMMMM", &gdlagXMMMM, &b_gdlagXMMMM);
   fChain->SetBranchAddress("gdlagYMMMM", &gdlagYMMMM, &b_gdlagYMMMM);
   fChain->SetBranchAddress("gdlagZMMMM", &gdlagZMMMM, &b_gdlagZMMMM);
   fChain->SetBranchAddress("gdlagProbMMMM", &gdlagProbMMMM, &b_gdlagProbMMMM);
   fChain->SetBranchAddress("gdlagNdofMMMM", &gdlagNdofMMMM, &b_gdlagNdofMMMM);
   fChain->SetBranchAddress("ftsigmaEEEE", &ftsigmaEEEE, &b_ftsigmaEEEE);
   fChain->SetBranchAddress("gdXEEEE", &gdXEEEE, &b_gdXEEEE);
   fChain->SetBranchAddress("gdYEEEE", &gdYEEEE, &b_gdYEEEE);
   fChain->SetBranchAddress("gdZEEEE", &gdZEEEE, &b_gdZEEEE);
   fChain->SetBranchAddress("ftsigmalagEEEE", &ftsigmalagEEEE, &b_ftsigmalagEEEE);
   fChain->SetBranchAddress("gdlagXEEEE", &gdlagXEEEE, &b_gdlagXEEEE);
   fChain->SetBranchAddress("gdlagYEEEE", &gdlagYEEEE, &b_gdlagYEEEE);
   fChain->SetBranchAddress("gdlagZEEEE", &gdlagZEEEE, &b_gdlagZEEEE);
   fChain->SetBranchAddress("gdlagProbEEEE", &gdlagProbEEEE, &b_gdlagProbEEEE);
   fChain->SetBranchAddress("gdlagNdofEEEE", &gdlagNdofEEEE, &b_gdlagNdofEEEE);
   fChain->SetBranchAddress("StdFitVertexX", &StdFitVertexX, &b_StdFitVertexX);
   fChain->SetBranchAddress("StdFitVertexY", &StdFitVertexY, &b_StdFitVertexY);
   fChain->SetBranchAddress("StdFitVertexZ", &StdFitVertexZ, &b_StdFitVertexZ);
   fChain->SetBranchAddress("StdFitVertexChi2r", &StdFitVertexChi2r, &b_StdFitVertexChi2r);
   fChain->SetBranchAddress("StdFitVertexProb", &StdFitVertexProb, &b_StdFitVertexProb);
   fChain->SetBranchAddress("StdFitVertexTrack_PT", &StdFitVertexTrack_PT, &b_StdFitVertexTrack_PT);
   fChain->SetBranchAddress("StdFitVertexTrack_ETA", &StdFitVertexTrack_ETA, &b_StdFitVertexTrack_ETA);
   fChain->SetBranchAddress("StdFitVertexTrack_PHI", &StdFitVertexTrack_PHI, &b_StdFitVertexTrack_PHI);
   fChain->SetBranchAddress("KinFitVertexX", &KinFitVertexX, &b_KinFitVertexX);
   fChain->SetBranchAddress("KinFitVertexY", &KinFitVertexY, &b_KinFitVertexY);
   fChain->SetBranchAddress("KinFitVertexZ", &KinFitVertexZ, &b_KinFitVertexZ);
   fChain->SetBranchAddress("KinFitVertexChi2r", &KinFitVertexChi2r, &b_KinFitVertexChi2r);
   fChain->SetBranchAddress("KinFitVertexProb", &KinFitVertexProb, &b_KinFitVertexProb);
   fChain->SetBranchAddress("StdFitVertexXMMMM", &StdFitVertexXMMMM, &b_StdFitVertexXMMMM);
   fChain->SetBranchAddress("StdFitVertexYMMMM", &StdFitVertexYMMMM, &b_StdFitVertexYMMMM);
   fChain->SetBranchAddress("StdFitVertexZMMMM", &StdFitVertexZMMMM, &b_StdFitVertexZMMMM);
   fChain->SetBranchAddress("StdFitVertexChi2rMMMM", &StdFitVertexChi2rMMMM, &b_StdFitVertexChi2rMMMM);
   fChain->SetBranchAddress("StdFitVertexProbMMMM", &StdFitVertexProbMMMM, &b_StdFitVertexProbMMMM);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_PT", &StdFitVertexTrackMMMM_PT, &b_StdFitVertexTrackMMMM_PT);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_ETA", &StdFitVertexTrackMMMM_ETA, &b_StdFitVertexTrackMMMM_ETA);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_PHI", &StdFitVertexTrackMMMM_PHI, &b_StdFitVertexTrackMMMM_PHI);
   fChain->SetBranchAddress("KinFitVertexXMMMM", &KinFitVertexXMMMM, &b_KinFitVertexXMMMM);
   fChain->SetBranchAddress("KinFitVertexYMMMM", &KinFitVertexYMMMM, &b_KinFitVertexYMMMM);
   fChain->SetBranchAddress("KinFitVertexZMMMM", &KinFitVertexZMMMM, &b_KinFitVertexZMMMM);
   fChain->SetBranchAddress("KinFitVertexChi2rMMMM", &KinFitVertexChi2rMMMM, &b_KinFitVertexChi2rMMMM);
   fChain->SetBranchAddress("KinFitVertexProbMMMM", &KinFitVertexProbMMMM, &b_KinFitVertexProbMMMM);
   fChain->SetBranchAddress("StdFitVertexXEEEE", &StdFitVertexXEEEE, &b_StdFitVertexXEEEE);
   fChain->SetBranchAddress("StdFitVertexYEEEE", &StdFitVertexYEEEE, &b_StdFitVertexYEEEE);
   fChain->SetBranchAddress("StdFitVertexZEEEE", &StdFitVertexZEEEE, &b_StdFitVertexZEEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rEEEE", &StdFitVertexChi2rEEEE, &b_StdFitVertexChi2rEEEE);
   fChain->SetBranchAddress("StdFitVertexProbEEEE", &StdFitVertexProbEEEE, &b_StdFitVertexProbEEEE);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_PT", &StdFitVertexTrackEEEE_PT, &b_StdFitVertexTrackEEEE_PT);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_ETA", &StdFitVertexTrackEEEE_ETA, &b_StdFitVertexTrackEEEE_ETA);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_PHI", &StdFitVertexTrackEEEE_PHI, &b_StdFitVertexTrackEEEE_PHI);
   fChain->SetBranchAddress("KinFitVertexXEEEE", &KinFitVertexXEEEE, &b_KinFitVertexXEEEE);
   fChain->SetBranchAddress("KinFitVertexYEEEE", &KinFitVertexYEEEE, &b_KinFitVertexYEEEE);
   fChain->SetBranchAddress("KinFitVertexZEEEE", &KinFitVertexZEEEE, &b_KinFitVertexZEEEE);
   fChain->SetBranchAddress("KinFitVertexChi2rEEEE", &KinFitVertexChi2rEEEE, &b_KinFitVertexChi2rEEEE);
   fChain->SetBranchAddress("KinFitVertexProbEEEE", &KinFitVertexProbEEEE, &b_KinFitVertexProbEEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rMMM", &StdFitVertexChi2rMMM, &b_StdFitVertexChi2rMMM);
   fChain->SetBranchAddress("StdFitVertexProbMMM", &StdFitVertexProbMMM, &b_StdFitVertexProbMMM);
   fChain->SetBranchAddress("StdFitVertexChi2rMME", &StdFitVertexChi2rMME, &b_StdFitVertexChi2rMME);
   fChain->SetBranchAddress("StdFitVertexProbMME", &StdFitVertexProbMME, &b_StdFitVertexProbMME);
   fChain->SetBranchAddress("StdFitVertexChi2rEEE", &StdFitVertexChi2rEEE, &b_StdFitVertexChi2rEEE);
   fChain->SetBranchAddress("StdFitVertexProbEEE", &StdFitVertexProbEEE, &b_StdFitVertexProbEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rMEE", &StdFitVertexChi2rMEE, &b_StdFitVertexChi2rMEE);
   fChain->SetBranchAddress("StdFitVertexProbMEE", &StdFitVertexProbMEE, &b_StdFitVertexProbMEE);
   fChain->SetBranchAddress("StdFitVertexChi2rDiLep", &StdFitVertexChi2rDiLep, &b_StdFitVertexChi2rDiLep);
   fChain->SetBranchAddress("StdFitVertexProbDiLep", &StdFitVertexProbDiLep, &b_StdFitVertexProbDiLep);
   fChain->SetBranchAddress("ConvMapDist", &ConvMapDist, &b_ConvMapDist);
   fChain->SetBranchAddress("ConvMapDcot", &ConvMapDcot, &b_ConvMapDcot);
   fChain->SetBranchAddress("RECOMU_MatchingMCTruth", &RECOMU_MatchingMCTruth, &b_RECOMU_MatchingMCTruth);
   fChain->SetBranchAddress("RECOMU_MatchingMCpT", &RECOMU_MatchingMCpT, &b_RECOMU_MatchingMCpT);
   fChain->SetBranchAddress("RECOMU_MatchingMCEta", &RECOMU_MatchingMCEta, &b_RECOMU_MatchingMCEta);
   fChain->SetBranchAddress("RECOMU_MatchingMCPhi", &RECOMU_MatchingMCPhi, &b_RECOMU_MatchingMCPhi);
   fChain->SetBranchAddress("RECOELE_MatchingMCTruth", &RECOELE_MatchingMCTruth, &b_RECOELE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOELE_MatchingMCpT", &RECOELE_MatchingMCpT, &b_RECOELE_MatchingMCpT);
   fChain->SetBranchAddress("RECOELE_MatchingMCEta", &RECOELE_MatchingMCEta, &b_RECOELE_MatchingMCEta);
   fChain->SetBranchAddress("RECOELE_MatchingMCPhi", &RECOELE_MatchingMCPhi, &b_RECOELE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCTruth", &RECOPHOT_MatchingMCTruth, &b_RECOPHOT_MatchingMCTruth);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCpT", &RECOPHOT_MatchingMCpT, &b_RECOPHOT_MatchingMCpT);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCEta", &RECOPHOT_MatchingMCEta, &b_RECOPHOT_MatchingMCEta);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCPhi", &RECOPHOT_MatchingMCPhi, &b_RECOPHOT_MatchingMCPhi);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCTruth", &RECOzMuMu_MatchingMCTruth, &b_RECOzMuMu_MatchingMCTruth);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCpT", &RECOzMuMu_MatchingMCpT, &b_RECOzMuMu_MatchingMCpT);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCmass", &RECOzMuMu_MatchingMCmass, &b_RECOzMuMu_MatchingMCmass);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCEta", &RECOzMuMu_MatchingMCEta, &b_RECOzMuMu_MatchingMCEta);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCPhi", &RECOzMuMu_MatchingMCPhi, &b_RECOzMuMu_MatchingMCPhi);
   fChain->SetBranchAddress("RECOzEE_MatchingMCTruth", &RECOzEE_MatchingMCTruth, &b_RECOzEE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOzEE_MatchingMCpT", &RECOzEE_MatchingMCpT, &b_RECOzEE_MatchingMCpT);
   fChain->SetBranchAddress("RECOzEE_MatchingMCmass", &RECOzEE_MatchingMCmass, &b_RECOzEE_MatchingMCmass);
   fChain->SetBranchAddress("RECOzEE_MatchingMCEta", &RECOzEE_MatchingMCEta, &b_RECOzEE_MatchingMCEta);
   fChain->SetBranchAddress("RECOzEE_MatchingMCPhi", &RECOzEE_MatchingMCPhi, &b_RECOzEE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCTruth", &RECOHzzEEEE_MatchingMCTruth, &b_RECOHzzEEEE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCpT", &RECOHzzEEEE_MatchingMCpT, &b_RECOHzzEEEE_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCmass", &RECOHzzEEEE_MatchingMCmass, &b_RECOHzzEEEE_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCEta", &RECOHzzEEEE_MatchingMCEta, &b_RECOHzzEEEE_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCPhi", &RECOHzzEEEE_MatchingMCPhi, &b_RECOHzzEEEE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCTruth", &RECOHzzEEMM_MatchingMCTruth, &b_RECOHzzEEMM_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCpT", &RECOHzzEEMM_MatchingMCpT, &b_RECOHzzEEMM_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCmass", &RECOHzzEEMM_MatchingMCmass, &b_RECOHzzEEMM_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCEta", &RECOHzzEEMM_MatchingMCEta, &b_RECOHzzEEMM_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCPhi", &RECOHzzEEMM_MatchingMCPhi, &b_RECOHzzEEMM_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCTruth", &RECOHzzMMMM_MatchingMCTruth, &b_RECOHzzMMMM_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCpT", &RECOHzzMMMM_MatchingMCpT, &b_RECOHzzMMMM_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCmass", &RECOHzzMMMM_MatchingMCmass, &b_RECOHzzMMMM_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCEta", &RECOHzzMMMM_MatchingMCEta, &b_RECOHzzMMMM_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCPhi", &RECOHzzMMMM_MatchingMCPhi, &b_RECOHzzMMMM_MatchingMCPhi);
   fChain->SetBranchAddress("RECO_NMU", &RECO_NMU, &b_RECO_NMU);
   fChain->SetBranchAddress("RECO_NELE", &RECO_NELE, &b_RECO_NELE);
   fChain->SetBranchAddress("RECO_NTRACK", &RECO_NTRACK, &b_RECO_NTRACK);
   fChain->SetBranchAddress("RECO_NPHOT", &RECO_NPHOT, &b_RECO_NPHOT);
   fChain->SetBranchAddress("RECOPHOT_PT", &RECOPHOT_PT, &b_RECOPHOT_PT);
   fChain->SetBranchAddress("RECOPHOT_ETA", &RECOPHOT_ETA, &b_RECOPHOT_ETA);
   fChain->SetBranchAddress("RECOPHOT_PHI", &RECOPHOT_PHI, &b_RECOPHOT_PHI);
   fChain->SetBranchAddress("RECOPHOT_THETA", &RECOPHOT_THETA, &b_RECOPHOT_THETA);
   fChain->SetBranchAddress("RECOPHOT_TLE_ParentSC_X", &RECOPHOT_TLE_ParentSC_X, &b_RECOPHOT_TLE_ParentSC_X);
   fChain->SetBranchAddress("RECOPHOT_TLE_ParentSC_Y", &RECOPHOT_TLE_ParentSC_Y, &b_RECOPHOT_TLE_ParentSC_Y);
   fChain->SetBranchAddress("RECOPHOT_TLE_ParentSC_Z", &RECOPHOT_TLE_ParentSC_Z, &b_RECOPHOT_TLE_ParentSC_Z);
   fChain->SetBranchAddress("RECO_NPFPHOT", &RECO_NPFPHOT, &b_RECO_NPFPHOT);
   fChain->SetBranchAddress("RECOPFPHOT_PT", &RECOPFPHOT_PT, &b_RECOPFPHOT_PT);
   fChain->SetBranchAddress("RECOPFPHOT_PTError", &RECOPFPHOT_PTError, &b_RECOPFPHOT_PTError);
   fChain->SetBranchAddress("RECOPFPHOT_ETA", &RECOPFPHOT_ETA, &b_RECOPFPHOT_ETA);
   fChain->SetBranchAddress("RECOPFPHOT_PHI", &RECOPFPHOT_PHI, &b_RECOPFPHOT_PHI);
   fChain->SetBranchAddress("RECOPFPHOT_THETA", &RECOPFPHOT_THETA, &b_RECOPFPHOT_THETA);
   fChain->SetBranchAddress("RECOPFPHOT_PT_uncorr", &RECOPFPHOT_PT_uncorr, &b_RECOPFPHOT_PT_uncorr);
   fChain->SetBranchAddress("RECOPFPHOT_ecalEnergyPreCorr", &RECOPFPHOT_ecalEnergyPreCorr, &b_RECOPFPHOT_ecalEnergyPreCorr);
   fChain->SetBranchAddress("RECOPFPHOT_ecalEnergyErrPreCorr", &RECOPFPHOT_ecalEnergyErrPreCorr, &b_RECOPFPHOT_ecalEnergyErrPreCorr);
   fChain->SetBranchAddress("RECOPFPHOT_ecalEnergyErrPostCorr", &RECOPFPHOT_ecalEnergyErrPostCorr, &b_RECOPFPHOT_ecalEnergyErrPostCorr);
   fChain->SetBranchAddress("RECOPFPHOT_energyScaleValue", &RECOPFPHOT_energyScaleValue, &b_RECOPFPHOT_energyScaleValue);
   fChain->SetBranchAddress("RECOPFPHOT_energySigmaValue", &RECOPFPHOT_energySigmaValue, &b_RECOPFPHOT_energySigmaValue);
   fChain->SetBranchAddress("RECOPFPHOT_energyScaleUp", &RECOPFPHOT_energyScaleUp, &b_RECOPFPHOT_energyScaleUp);
   fChain->SetBranchAddress("RECOPFPHOT_energyScaleDown", &RECOPFPHOT_energyScaleDown, &b_RECOPFPHOT_energyScaleDown);
   fChain->SetBranchAddress("RECOPFPHOT_energyScaleStatUp", &RECOPFPHOT_energyScaleStatUp, &b_RECOPFPHOT_energyScaleStatUp);
   fChain->SetBranchAddress("RECOPFPHOT_energyScaleStatDown", &RECOPFPHOT_energyScaleStatDown, &b_RECOPFPHOT_energyScaleStatDown);
   fChain->SetBranchAddress("RECOPFPHOT_energyScaleSystUp", &RECOPFPHOT_energyScaleSystUp, &b_RECOPFPHOT_energyScaleSystUp);
   fChain->SetBranchAddress("RECOPFPHOT_energyScaleSystDown", &RECOPFPHOT_energyScaleSystDown, &b_RECOPFPHOT_energyScaleSystDown);
   fChain->SetBranchAddress("RECOPFPHOT_energyScaleGainUp", &RECOPFPHOT_energyScaleGainUp, &b_RECOPFPHOT_energyScaleGainUp);
   fChain->SetBranchAddress("RECOPFPHOT_energyScaleGainDown", &RECOPFPHOT_energyScaleGainDown, &b_RECOPFPHOT_energyScaleGainDown);
   fChain->SetBranchAddress("RECOPFPHOT_energyScaleEtUp", &RECOPFPHOT_energyScaleEtUp, &b_RECOPFPHOT_energyScaleEtUp);
   fChain->SetBranchAddress("RECOPFPHOT_energyScaleEtDown", &RECOPFPHOT_energyScaleEtDown, &b_RECOPFPHOT_energyScaleEtDown);
   fChain->SetBranchAddress("RECOPFPHOT_energySigmaUp", &RECOPFPHOT_energySigmaUp, &b_RECOPFPHOT_energySigmaUp);
   fChain->SetBranchAddress("RECOPFPHOT_energySigmaDown", &RECOPFPHOT_energySigmaDown, &b_RECOPFPHOT_energySigmaDown);
   fChain->SetBranchAddress("RECOPFPHOT_energySigmaPhiUp", &RECOPFPHOT_energySigmaPhiUp, &b_RECOPFPHOT_energySigmaPhiUp);
   fChain->SetBranchAddress("RECOPFPHOT_energySigmaPhiDown", &RECOPFPHOT_energySigmaPhiDown, &b_RECOPFPHOT_energySigmaPhiDown);
   fChain->SetBranchAddress("RECOPFPHOT_energySigmaRhoUp", &RECOPFPHOT_energySigmaRhoUp, &b_RECOPFPHOT_energySigmaRhoUp);
   fChain->SetBranchAddress("RECOPFPHOT_energySigmaRhoDown", &RECOPFPHOT_energySigmaRhoDown, &b_RECOPFPHOT_energySigmaRhoDown);
   fChain->SetBranchAddress("BeamSpot_X", &BeamSpot_X, &b_BeamSpot_X);
   fChain->SetBranchAddress("BeamSpot_Y", &BeamSpot_Y, &b_BeamSpot_Y);
   fChain->SetBranchAddress("BeamSpot_Z", &BeamSpot_Z, &b_BeamSpot_Z);
   fChain->SetBranchAddress("RECO_NVTX", &RECO_NVTX, &b_RECO_NVTX);
   fChain->SetBranchAddress("RECO_VERTEX_x", &RECO_VERTEX_x, &b_RECO_VERTEX_x);
   fChain->SetBranchAddress("RECO_VERTEX_y", &RECO_VERTEX_y, &b_RECO_VERTEX_y);
   fChain->SetBranchAddress("RECO_VERTEX_z", &RECO_VERTEX_z, &b_RECO_VERTEX_z);
   fChain->SetBranchAddress("RECO_VERTEX_ndof", &RECO_VERTEX_ndof, &b_RECO_VERTEX_ndof);
   fChain->SetBranchAddress("RECO_VERTEX_chi2", &RECO_VERTEX_chi2, &b_RECO_VERTEX_chi2);
   fChain->SetBranchAddress("RECO_VERTEX_ntracks", &RECO_VERTEX_ntracks, &b_RECO_VERTEX_ntracks);
   fChain->SetBranchAddress("RECO_VERTEXPROB", &RECO_VERTEXPROB, &b_RECO_VERTEXPROB);
   fChain->SetBranchAddress("RECO_VERTEX_isValid", &RECO_VERTEX_isValid, &b_RECO_VERTEX_isValid);
   fChain->SetBranchAddress("RECO_PFJET_N", &RECO_PFJET_N, &b_RECO_PFJET_N);
   fChain->SetBranchAddress("RECO_PFJET_CHARGE", &RECO_PFJET_CHARGE, &b_RECO_PFJET_CHARGE);
   fChain->SetBranchAddress("RECO_PFJET_ET", &RECO_PFJET_ET, &b_RECO_PFJET_ET);
   fChain->SetBranchAddress("RECO_PFJET_PT", &RECO_PFJET_PT, &b_RECO_PFJET_PT);
   fChain->SetBranchAddress("RECO_PFJET_ETA", &RECO_PFJET_ETA, &b_RECO_PFJET_ETA);
   fChain->SetBranchAddress("RECO_PFJET_PHI", &RECO_PFJET_PHI, &b_RECO_PFJET_PHI);
   fChain->SetBranchAddress("RECO_PFJET_PUID_loose", &RECO_PFJET_PUID_loose, &b_RECO_PFJET_PUID_loose);
   fChain->SetBranchAddress("RECO_PFJET_PUID_medium", &RECO_PFJET_PUID_medium, &b_RECO_PFJET_PUID_medium);
   fChain->SetBranchAddress("RECO_PFJET_PUID", &RECO_PFJET_PUID, &b_RECO_PFJET_PUID);
   fChain->SetBranchAddress("RECO_PFJET_PUID_MVA", &RECO_PFJET_PUID_MVA, &b_RECO_PFJET_PUID_MVA);
   fChain->SetBranchAddress("RECO_PFJET_QG_Likelihood", &RECO_PFJET_QG_Likelihood, &b_RECO_PFJET_QG_Likelihood);
   fChain->SetBranchAddress("RECO_PFJET_QG_axis2", &RECO_PFJET_QG_axis2, &b_RECO_PFJET_QG_axis2);
   fChain->SetBranchAddress("RECO_PFJET_QG_ptd", &RECO_PFJET_QG_ptd, &b_RECO_PFJET_QG_ptd);
   fChain->SetBranchAddress("RECO_PFJET_QG_mult", &RECO_PFJET_QG_mult, &b_RECO_PFJET_QG_mult);
   fChain->SetBranchAddress("RHO_ele", &RHO_ele, &b_RHO_ele);
   fChain->SetBranchAddress("RHO_mu", &RHO_mu, &b_RHO_mu);
   fChain->SetBranchAddress("LHE_PARTON_N", &LHE_PARTON_N, &b_LHE_PARTON_N);
   fChain->SetBranchAddress("LHE_PARTON_CLEAR", &LHE_PARTON_CLEAR, &b_LHE_PARTON_CLEAR);
   fChain->SetBranchAddress("LHE_PARTON_PDGID", &LHE_PARTON_PDGID, &b_LHE_PARTON_PDGID);
   fChain->SetBranchAddress("LHE_PARTON_PT", &LHE_PARTON_PT, &b_LHE_PARTON_PT);
   fChain->SetBranchAddress("LHE_PARTON_ETA", &LHE_PARTON_ETA, &b_LHE_PARTON_ETA);
   fChain->SetBranchAddress("LHE_PARTON_PHI", &LHE_PARTON_PHI, &b_LHE_PARTON_PHI);
   fChain->SetBranchAddress("LHE_PARTON_E", &LHE_PARTON_E, &b_LHE_PARTON_E);
   fChain->SetBranchAddress("RECO_PFJET_PT_UncUp", &RECO_PFJET_PT_UncUp, &b_RECO_PFJET_PT_UncUp);
   fChain->SetBranchAddress("RECO_PFJET_PT_UncDn", &RECO_PFJET_PT_UncDn, &b_RECO_PFJET_PT_UncDn);
   fChain->SetBranchAddress("RECO_PFJET_AREA", &RECO_PFJET_AREA, &b_RECO_PFJET_AREA);
   fChain->SetBranchAddress("RECO_PFJET_PTD", &RECO_PFJET_PTD, &b_RECO_PFJET_PTD);
   fChain->SetBranchAddress("RECO_PFJET_CHARGED_HADRON_ENERGY", &RECO_PFJET_CHARGED_HADRON_ENERGY, &b_RECO_PFJET_CHARGED_HADRON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_NEUTRAL_HADRON_ENERGY", &RECO_PFJET_NEUTRAL_HADRON_ENERGY, &b_RECO_PFJET_NEUTRAL_HADRON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_PHOTON_ENERGY", &RECO_PFJET_PHOTON_ENERGY, &b_RECO_PFJET_PHOTON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_ELECTRON_ENERGY", &RECO_PFJET_ELECTRON_ENERGY, &b_RECO_PFJET_ELECTRON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_MUON_ENERGY", &RECO_PFJET_MUON_ENERGY, &b_RECO_PFJET_MUON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_HF_HADRON_ENERGY", &RECO_PFJET_HF_HADRON_ENERGY, &b_RECO_PFJET_HF_HADRON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_HF_EM_ENERGY", &RECO_PFJET_HF_EM_ENERGY, &b_RECO_PFJET_HF_EM_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_CHARGED_EM_ENERGY", &RECO_PFJET_CHARGED_EM_ENERGY, &b_RECO_PFJET_CHARGED_EM_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_CHARGED_MU_ENERGY", &RECO_PFJET_CHARGED_MU_ENERGY, &b_RECO_PFJET_CHARGED_MU_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_NEUTRAL_EM_ENERGY", &RECO_PFJET_NEUTRAL_EM_ENERGY, &b_RECO_PFJET_NEUTRAL_EM_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_CHARGED_HADRON_MULTIPLICITY", &RECO_PFJET_CHARGED_HADRON_MULTIPLICITY, &b_RECO_PFJET_CHARGED_HADRON_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY", &RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY, &b_RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_PHOTON_MULTIPLICITY", &RECO_PFJET_PHOTON_MULTIPLICITY, &b_RECO_PFJET_PHOTON_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_ELECTRON_MULTIPLICITY", &RECO_PFJET_ELECTRON_MULTIPLICITY, &b_RECO_PFJET_ELECTRON_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_MUON_MULTIPLICITY", &RECO_PFJET_MUON_MULTIPLICITY, &b_RECO_PFJET_MUON_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_HF_HADRON_MULTIPLICTY", &RECO_PFJET_HF_HADRON_MULTIPLICTY, &b_RECO_PFJET_HF_HADRON_MULTIPLICTY);
   fChain->SetBranchAddress("RECO_PFJET_HF_EM_MULTIPLICITY", &RECO_PFJET_HF_EM_MULTIPLICITY, &b_RECO_PFJET_HF_EM_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_CHARGED_MULTIPLICITY", &RECO_PFJET_CHARGED_MULTIPLICITY, &b_RECO_PFJET_CHARGED_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_NEUTRAL_MULTIPLICITY", &RECO_PFJET_NEUTRAL_MULTIPLICITY, &b_RECO_PFJET_NEUTRAL_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_NCOMPONENTS", &RECO_PFJET_NCOMPONENTS, &b_RECO_PFJET_NCOMPONENTS);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_PDGID", &RECO_PFJET_COMPONENT_PDGID, &b_RECO_PFJET_COMPONENT_PDGID);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_PT", &RECO_PFJET_COMPONENT_PT, &b_RECO_PFJET_COMPONENT_PT);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_ETA", &RECO_PFJET_COMPONENT_ETA, &b_RECO_PFJET_COMPONENT_ETA);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_PHI", &RECO_PFJET_COMPONENT_PHI, &b_RECO_PFJET_COMPONENT_PHI);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_E", &RECO_PFJET_COMPONENT_E, &b_RECO_PFJET_COMPONENT_E);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_CHARGE", &RECO_PFJET_COMPONENT_CHARGE, &b_RECO_PFJET_COMPONENT_CHARGE);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_TRANSVERSE_MASS", &RECO_PFJET_COMPONENT_TRANSVERSE_MASS, &b_RECO_PFJET_COMPONENT_TRANSVERSE_MASS);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_XVERTEX", &RECO_PFJET_COMPONENT_XVERTEX, &b_RECO_PFJET_COMPONENT_XVERTEX);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_YVERTEX", &RECO_PFJET_COMPONENT_YVERTEX, &b_RECO_PFJET_COMPONENT_YVERTEX);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_ZVERTEX", &RECO_PFJET_COMPONENT_ZVERTEX, &b_RECO_PFJET_COMPONENT_ZVERTEX);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_VERTEX_CHI2", &RECO_PFJET_COMPONENT_VERTEX_CHI2, &b_RECO_PFJET_COMPONENT_VERTEX_CHI2);
   fChain->SetBranchAddress("RECO_CALOMET", &RECO_CALOMET, &b_RECO_CALOMET);
   fChain->SetBranchAddress("RECO_PFMET", &RECO_PFMET, &b_RECO_PFMET);
   fChain->SetBranchAddress("RECO_PFMET_X", &RECO_PFMET_X, &b_RECO_PFMET_X);
   fChain->SetBranchAddress("RECO_PFMET_Y", &RECO_PFMET_Y, &b_RECO_PFMET_Y);
   fChain->SetBranchAddress("RECO_PFMET_PHI", &RECO_PFMET_PHI, &b_RECO_PFMET_PHI);
   fChain->SetBranchAddress("RECO_PFMET_THETA", &RECO_PFMET_THETA, &b_RECO_PFMET_THETA);
   fChain->SetBranchAddress("RECO_PFMET_uncorr", &RECO_PFMET_uncorr, &b_RECO_PFMET_uncorr);
   fChain->SetBranchAddress("RECO_PFMET_X_uncorr", &RECO_PFMET_X_uncorr, &b_RECO_PFMET_X_uncorr);
   fChain->SetBranchAddress("RECO_PFMET_Y_uncorr", &RECO_PFMET_Y_uncorr, &b_RECO_PFMET_Y_uncorr);
   fChain->SetBranchAddress("RECO_PFMET_PHI_uncorr", &RECO_PFMET_PHI_uncorr, &b_RECO_PFMET_PHI_uncorr);
   fChain->SetBranchAddress("RECO_PFMET_THETA_uncorr", &RECO_PFMET_THETA_uncorr, &b_RECO_PFMET_THETA_uncorr);
   fChain->SetBranchAddress("RECO_PFMET_JetEnUp", &RECO_PFMET_JetEnUp, &b_RECO_PFMET_JetEnUp);
   fChain->SetBranchAddress("RECO_PFMET_JetEnDn", &RECO_PFMET_JetEnDn, &b_RECO_PFMET_JetEnDn);
   fChain->SetBranchAddress("RECO_PFMET_ElectronEnUp", &RECO_PFMET_ElectronEnUp, &b_RECO_PFMET_ElectronEnUp);
   fChain->SetBranchAddress("RECO_PFMET_ElectronEnDn", &RECO_PFMET_ElectronEnDn, &b_RECO_PFMET_ElectronEnDn);
   fChain->SetBranchAddress("RECO_PFMET_MuonEnUp", &RECO_PFMET_MuonEnUp, &b_RECO_PFMET_MuonEnUp);
   fChain->SetBranchAddress("RECO_PFMET_MuonEnDn", &RECO_PFMET_MuonEnDn, &b_RECO_PFMET_MuonEnDn);
   fChain->SetBranchAddress("RECO_PFMET_JetResUp", &RECO_PFMET_JetResUp, &b_RECO_PFMET_JetResUp);
   fChain->SetBranchAddress("RECO_PFMET_JetResDn", &RECO_PFMET_JetResDn, &b_RECO_PFMET_JetResDn);
   fChain->SetBranchAddress("RECO_PFMET_UnclusteredEnUp", &RECO_PFMET_UnclusteredEnUp, &b_RECO_PFMET_UnclusteredEnUp);
   fChain->SetBranchAddress("RECO_PFMET_UnclusteredEnDn", &RECO_PFMET_UnclusteredEnDn, &b_RECO_PFMET_UnclusteredEnDn);
   fChain->SetBranchAddress("RECO_PFMET_PhotonEnUp", &RECO_PFMET_PhotonEnUp, &b_RECO_PFMET_PhotonEnUp);
   fChain->SetBranchAddress("RECO_PFMET_PhotonEnDn", &RECO_PFMET_PhotonEnDn, &b_RECO_PFMET_PhotonEnDn);
   fChain->SetBranchAddress("RECO_PFMET_TauEnUp ", &RECO_PFMET_TauEnUp , &b_RECO_PFMET_TauEnUp );
   fChain->SetBranchAddress("RECO_PFMET_TauEnDown", &RECO_PFMET_TauEnDown, &b_RECO_PFMET_TauEnDown);
   fChain->SetBranchAddress("RECO_PFMET_GoodVtxNoiseFilter", &RECO_PFMET_GoodVtxNoiseFilter, &b_RECO_PFMET_GoodVtxNoiseFilter);
   fChain->SetBranchAddress("RECO_PFMET_GlobalSuperTightHalo2016NoiseFilter", &RECO_PFMET_GlobalSuperTightHalo2016NoiseFilter, &b_RECO_PFMET_GlobalSuperTightHalo2016NoiseFilter);
   fChain->SetBranchAddress("RECO_PFMET_HBHENoiseFilter", &RECO_PFMET_HBHENoiseFilter, &b_RECO_PFMET_HBHENoiseFilter);
   fChain->SetBranchAddress("RECO_PFMET_HBHENoiseIsoFilter", &RECO_PFMET_HBHENoiseIsoFilter, &b_RECO_PFMET_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("RECO_PFMET_EcalDeadCellTriggerPrimitiveNoiseFilter", &RECO_PFMET_EcalDeadCellTriggerPrimitiveNoiseFilter, &b_RECO_PFMET_EcalDeadCellTriggerPrimitiveNoiseFilter);
   fChain->SetBranchAddress("RECO_PFMET_BadPFMuonFilter", &RECO_PFMET_BadPFMuonFilter, &b_RECO_PFMET_BadPFMuonFilter);
   fChain->SetBranchAddress("RECO_PFMET_BadChargedCandidateFilter", &RECO_PFMET_BadChargedCandidateFilter, &b_RECO_PFMET_BadChargedCandidateFilter);
   fChain->SetBranchAddress("RECO_PFMET_EEBadScNoiseFilter", &RECO_PFMET_EEBadScNoiseFilter, &b_RECO_PFMET_EEBadScNoiseFilter);
   fChain->SetBranchAddress("RECO_PFMET_EcalBadCalibFilter", &RECO_PFMET_EcalBadCalibFilter, &b_RECO_PFMET_EcalBadCalibFilter);
   fChain->SetBranchAddress("RECO_TCMET", &RECO_TCMET, &b_RECO_TCMET);
   fChain->SetBranchAddress("RECO_CORMETMUONS", &RECO_CORMETMUONS, &b_RECO_CORMETMUONS);
   fChain->SetBranchAddress("tCHighEff_BTagJet_PT", &tCHighEff_BTagJet_PT, &b_tCHighEff_BTagJet_PT);
   fChain->SetBranchAddress("tCHighPur_BTagJet_PT", &tCHighPur_BTagJet_PT, &b_tCHighPur_BTagJet_PT);
   fChain->SetBranchAddress("cSV_BTagJet_PT", &cSV_BTagJet_PT, &b_cSV_BTagJet_PT);
   fChain->SetBranchAddress("tCHighEff_BTagJet_ETA", &tCHighEff_BTagJet_ETA, &b_tCHighEff_BTagJet_ETA);
   fChain->SetBranchAddress("tCHighPur_BTagJet_ETA", &tCHighPur_BTagJet_ETA, &b_tCHighPur_BTagJet_ETA);
   fChain->SetBranchAddress("cSV_BTagJet_ETA", &cSV_BTagJet_ETA, &b_cSV_BTagJet_ETA);
   fChain->SetBranchAddress("tCHighEff_BTagJet_PHI", &tCHighEff_BTagJet_PHI, &b_tCHighEff_BTagJet_PHI);
   fChain->SetBranchAddress("tCHighPur_BTagJet_PHI", &tCHighPur_BTagJet_PHI, &b_tCHighPur_BTagJet_PHI);
   fChain->SetBranchAddress("cSV_BTagJet_PHI", &cSV_BTagJet_PHI, &b_cSV_BTagJet_PHI);
   fChain->SetBranchAddress("cSV_BTagJet_ET", &cSV_BTagJet_ET, &b_cSV_BTagJet_ET);
   fChain->SetBranchAddress("tCHighEff_BTagJet_DISCR", &tCHighEff_BTagJet_DISCR, &b_tCHighEff_BTagJet_DISCR);
   fChain->SetBranchAddress("tCHighPur_BTagJet_DISCR", &tCHighPur_BTagJet_DISCR, &b_tCHighPur_BTagJet_DISCR);
   fChain->SetBranchAddress("cSV_BTagJet_DISCR", &cSV_BTagJet_DISCR, &b_cSV_BTagJet_DISCR);
   std::cout << "Init(TTree*) called.\n";
   Notify();
}

Bool_t HZZ4LeptonsAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HZZ4LeptonsAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HZZ4LeptonsAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#if (product_4e || product_2e2mu)
void HZZ4LeptonsAnalysis::printelebnn(int i){    
  bnn_file 
                << (*RECOELE_PT)[i]  << " "  
                << (*RECOELE_ETA)[i] << " "  
                << (*RECOELE_PHI)[i] << " "  
                << (*RECOELE_CHARGE)[i] << " "
                << (*RECOELE_PFX_rho)[i] << " "
                << (*RECOELE_SIP)[i] << " ";
}
#endif

#if (product_4mu || product_2e2mu)
void HZZ4LeptonsAnalysis::printmubnn(int i){    
  bnn_file 
                << (*RECOMU_PT)[i]  << " "  
                << (*RECOMU_ETA)[i] << " "  
                << (*RECOMU_PHI)[i] << " "  
                << (*RECOMU_CHARGE)[i] << " "
                << (*RECOMU_PFX_dB)[i] << " "
                << (*RECOMU_SIP)[i] << " ";
}
#endif
//https://stackoverflow.com/q/17074324

#include <numeric>
#include <algorithm>

template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
    const std::vector<T>& vec,
    const Compare& compare)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
        [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
    return p;
}

template <typename T>
void apply_permutation(
    std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
        if (done[i])
        {
            continue;
        }
        done[i] = true;
        std::size_t prev_j = i;
        std::size_t j = p[i];
        while (i != j)
        {
            std::swap(vec[prev_j], vec[j]);
            done[j] = true;
            prev_j = j;
            j = p[j];
        }
    }
}

#include <cassert>

template <class T>
class defaults {
    static T def;
public:
    static T& get_default() { 
        if(def!=-999) { std::cout << "Warning! Value of vector changed out of bounds. Probably a bug.\n"; def=-999; }
        return def;
    }
};

template <class T>
T defaults<T>::def = -999;

template <class T,class=decltype(std::declval<T&>() = std::declval<int>())>
T& get_default_helper(T*) { return defaults<T>::get_default(); }
template <class T>
T& get_default_helper(...) { std::cout << "Called on vector. Exiting.\n" ; assert(false) ; }
template <class T>
T& get_default() { return get_default_helper<T>((T*)0); }
//Proxy class that returns -999 if a vector of int, double, etc., is accessed out of bounds. Hack to get macro code working quickly

struct pair_hash
{
    template<class T,class U>
    std::size_t operator()(const std::pair<T,U>& pair) const
    {
        return std::hash<T>()(pair.first) ^ std::hash<std::size_t>()(std::hash<U>()(pair.second));
    }
};

template <
    class T,
    class Hash = std::hash<T>,
    class KeyEqual = std::equal_to<T>,
    class Allocator = std::allocator<T>
>
class doOnce {
private:
    std::unordered_set<T,Hash,KeyEqual,Allocator> done;
    bool (*equal_to)(const T& x,const T& y);
    void (*to_do)(const T& x);
public:
    doOnce(void (*to_do_f)(const T&)) : to_do(to_do_f) , done() {}
    bool operator()(const T& x) {
        //Insert the element into set of done elements, return if insertion fails.
        if(!(done.insert(x).second)) return false;
        to_do(x);
        return true;
    }
};

doOnce<std::pair<const char*,std::size_t>,pair_hash> warn([](const std::pair<const char*,std::size_t>& x)->void{ std::cout << "Warning: vector \"" << x.first << "\" accessed out of bounds, line number " << x.second << "." << std::endl; });

template <class T>
class safeAccessHelper {
    std::vector<T>* myVector;
    const char* var_name;
    std::size_t line_number;
public:
    safeAccessHelper(std::vector<T>* x,const char* name,std::size_t line) : myVector(x),var_name(name),line_number(line) {}
    typename std::vector<T>::reference operator[](std::size_t n) {
        warn(std::make_pair(var_name,line_number));
        if(n < myVector->size()) {
            return (*myVector)[n];
        } else {
            return get_default<T>();
        }
    }
    typename std::vector<T>::const_reference operator[](std::size_t n) const {
        warn(std::make_pair(var_name,line_number));
        if(n < myVector->size()) {
            return (*myVector)[n];
        } else {
            return get_default<T>();
        }
    }
};

//Wrapper for class, so it can be called without cumbersome boilerplate like safeAccess<decltype(vector_name)::value_type>(vector_name)

template <class T>
safeAccessHelper<T> safeAccessF(std::vector<T>* vec,const char* var_name,std::size_t line) { return safeAccessHelper<T>(vec,var_name,line); }

#define safeAccess(var) (safeAccessF((var),#var,__LINE__))

//Dummy cout to suppress useless output
#include <iostream>
class coutImpl : public std::basic_ostream<char> {};
template<typename T>
coutImpl& operator<<(coutImpl& cout,T& t) { return cout; }

#endif // #ifdef HZZ4LeptonsAnalysis_cxx
