#include "pileup_corrector.h"
#include <iostream>
#include <TFile.h>

#define FAIL() fail(__LINE__)

void pileup_corrector::fail(int line) {
    std::cout << "Pileup corrector failed on line " << line << ".\n";
    exit(1);
}

double pileup_corrector::get_pileup_weight(int num_pu_vertices) {
    if(isMC_) {
        int bin=pileup_ratio_->GetXaxis()->FindBin(num_pu_vertices);
        return pileup_ratio_->GetBinContent(bin);
    } else {
        return 1.;
    }
}
pileup_corrector::pileup_corrector(bool isMC,std::string era) : isMC_(isMC), pileup_ratio_(nullptr) {
    if(isMC_) {
        std::cout << "Getting pileup corrections for " << era << " MC.\n";
        if(era=="Fall17") { //Fall17 corrections
            TFile* puf=TFile::Open("PU_Reweight_2017.root");
            if(!puf) FAIL();
            pileup_ratio_=(TH1D*)(puf->Get("PU_Ratio"));
            if(!pileup_ratio_) FAIL();
            pileup_ratio_->SetDirectory(0);
            puf->Close();
        } else if(era=="Autumn18") { //Autumn18 corrections
            TFile* data_pu_file=TFile::Open("DataPileupHistogram2018_69200_100bins.root");
            TFile* mc_pu_file=TFile::Open("pu_weights_2018.root");
            if(!(data_pu_file && mc_pu_file)) FAIL();
            //Get data histogram
            pileup_ratio_=(TH1D*)(data_pu_file->Get("pileup"));
            if(!pileup_ratio_) FAIL();
            pileup_ratio_->SetDirectory(0);
            //Normalize data histogram to unity
            pileup_ratio_->Scale(1./pileup_ratio_->GetSum());
            //Get MC histogram
            TH1F* mcpu=(TH1F*)(mc_pu_file->Get("MC_out_of_the_box"));
            if(!mcpu) FAIL();
            pileup_ratio_->Divide(mcpu);
            pileup_ratio_->SetDirectory(0);
            data_pu_file->Close();
            mc_pu_file->Close();
        } else {
            std::cout << "Pileup corrections for MC era " << era << " not available\n.";
            exit(1);
        }
    } else {
        std::cout << "Running on data: no pileup reweighting required.\n";
    }
}
pileup_corrector::~pileup_corrector() {
    delete(pileup_ratio_);
}
