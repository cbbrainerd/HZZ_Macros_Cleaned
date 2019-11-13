void HZZ4LeptonsAnalysis::fill_xycorrected_met() {
    //Currently applied on Type1 MET
    auto xy_corrected_met=METXYCorr_Met_MetPhi(RECO_PFMET,RECO_PFMET_PHI,Run,year,isMC,RECO_NVTX);
    RECO_PFMET_xcorr=xy_corrected_met.first;
    RECO_PFMET_PHI_xcorr=xy_corrected_met.second;
}
