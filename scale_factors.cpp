#include "scale_factors.h"
#include <iostream>
#include "TFile.h"

#define FAIL() fail(__LINE__)

static void fail(int line) {
    std::cout << "Scale factors module failed on line " << line << ".\n";
    exit(1);
}

static inline float get_from_th2f(TH2F* hist,double x,double y) {
    Int_t bin_x=hist->GetXaxis()->FindFixBin(x);
    Int_t bin_y=hist->GetYaxis()->FindFixBin(y);
    //Overflow/underflow checks
    Int_t nbins_x=hist->GetXaxis()->GetNbins();
    Int_t nbins_y=hist->GetYaxis()->GetNbins();
    if(bin_x == 0) bin_x=1;
    if(bin_x > nbins_x) bin_x=nbins_x;
    if(bin_y == 0) bin_y=1;
    if(bin_y > nbins_y) bin_y=nbins_y;
    return hist->GetBinContent(bin_x,bin_y);
}

/*
template <class T>
T* get_from_root_file_and_put_in_pointer(const char* filename,const char* get_name,T*& put_here) {
    TFile* tfile=TFile::Open(
}
*/

void scale_factors::select_scale_factor_files(const char* gap,const char* no_gap) {
    TFile* sf=TFile::Open(gap); if(!sf) FAIL();
    scale_factor_is_gap=(TH2F*)(sf->Get(SFHistName)); if(!scale_factor_is_gap) FAIL();
    scale_factor_is_gap->SetDirectory(0);
    sf->Close();
    //For muons, just one hist rather than is_gap vs no_gap
    if(no_gap==nullptr) {
        scale_factor_no_gap=scale_factor_is_gap;
    } else {
        sf=TFile::Open(no_gap); if(!sf) FAIL();
        scale_factor_no_gap=(TH2F*)(sf->Get(SFHistName)); if(!scale_factor_no_gap) FAIL();
        scale_factor_no_gap->SetDirectory(0);
        sf->Close();
    }
}

void scale_factors_and_efficiencies::select_efficiency_files(const char* low_pt,const char* high_pt) {
    TFile* efficiency=TFile::Open(low_pt); if(!efficiency) FAIL();
    efficiency_low_pt=(TH2F*)(efficiency->Get(EffHistName)); if(!efficiency_low_pt) FAIL();
    efficiency_low_pt->SetDirectory(0);
    efficiency->Close();
    efficiency=TFile::Open(high_pt); if(!efficiency) FAIL();
    efficiency_high_pt=(TH2F*)(efficiency->Get(EffHistName)); if(!efficiency_high_pt) FAIL();
    efficiency_high_pt->SetDirectory(0);
    efficiency->Close();
}

scale_factors_and_efficiencies::scale_factors_and_efficiencies(bool isMC,std::string era,std::string particle_type) : scale_factors(isMC,era,particle_type) , efficiency_low_pt(nullptr) , efficiency_high_pt(nullptr) {
    std::cout << "Setting up efficiencies for " << particle_type << "s on " << (isMC?"MC":"data") << " for era " << era << ".\n";
    if(particle_type=="electron" or particle_type=="e") {
        EffHistName="EGamma_SF2D";
        if(isMC&&era=="Fall17") {
            select_efficiency_files("Ele_Reco_LowEt_2017.root","Ele_Reco_2017.root");
            return;
        } else if(isMC&&era=="Autumn18") {
            select_efficiency_files("Ele_Reco_LowEt_2018.root","Ele_Reco_2018.root");
            return;
        } else if((!isMC)&&(era=="2017" || era=="2018")) {
            return;
        }
    }
    std::cerr << "Efficiencies not defined for particle type \"" << particle_type << "\" on " << (isMC?"MC":"data") << " for era \"" << era << "\".\n";
    exit(1);
}

scale_factors::scale_factors(bool isMC,std::string era,std::string particle_type) : isMC_(isMC) , scale_factor_is_gap(nullptr) , scale_factor_no_gap(nullptr) {
    std::cout << "Setting up scale factors for " << particle_type << "s on " << (isMC?"MC":"data") << " for era " << era << ".\n";
    if(particle_type=="electron" || particle_type=="e") {
        SFHistName="EGamma_SF2D";
        if(isMC&&era=="Fall17") {
            select_scale_factor_files("ElectronSF_Legacy_2017_Gap.root","ElectronSF_Legacy_2017_NoGap.root");
            return;
        } else if(isMC&&era=="Autumn18") {
            select_scale_factor_files("ElectronSF_Legacy_2018_Gap.root","ElectronSF_Legacy_2018_NoGap.root");
            return;
        } else if((!isMC)&&(era=="2017" || era=="2018")) {
            return;
        }
    } else if(particle_type=="muon" || particle_type=="m") {
        SFHistName="FINAL";
        if(isMC&&era=="Fall17") {
            select_scale_factor_files("final_HZZ_SF_2017_rereco_mupogsysts_3010.root");
            return;
        } else if(isMC&&era=="Autumn18") {
            select_scale_factor_files("final_HZZ_SF_2018_rereco_mupogsysts_3010.root");
            return;
        } else if((!isMC)&&(era=="2017" || era=="2018")) {
            return;
        }
    }
    std::cerr << "Scale factors not defined for particle type \"" << particle_type << "\" on " << (isMC?"MC":"data") << " for era \"" << era << "\".\n";
    exit(1);
}

float scale_factors::get_scale_factor(double eta,double pt,bool isGap) {
    //if(pt>500) pt=500; //Overflow taken care of in get_from_th2f
    if(isGap && scale_factor_is_gap) {
        return get_from_th2f(scale_factor_is_gap,eta,pt);
    } else if (!isGap) {
        return get_from_th2f(scale_factor_no_gap,eta,pt);
    } else {
        std::cout << "Error: scale factors called with isGap==True but without a gap histogram.\n";
        exit(1);
    }
}

float scale_factors_and_efficiencies::get_efficiency(double eta,double pt) {
    //if(pt>500)pt=500; // Overflow taken care of in get_from_th2f
    if(pt>20) {
        return get_from_th2f(efficiency_high_pt,eta,pt);
    } else {
        return get_from_th2f(efficiency_low_pt,eta,pt);
    }
}

scale_factors::~scale_factors() {
    delete scale_factor_no_gap;
    delete scale_factor_is_gap;
}

scale_factors_and_efficiencies::~scale_factors_and_efficiencies() {
    delete efficiency_low_pt;
    delete efficiency_high_pt;
}
