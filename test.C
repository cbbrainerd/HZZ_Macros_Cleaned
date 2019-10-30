TTree* test() {
    TChain *tc=new TChain("HZZ4LeptonsAnalysis","");
    tc->Add("root://cms-xrd-global.cern.ch//store/user/cbrainer/DoubleEG/crab_HZZ4LeptonsAnalysis_DoubleEG_Run2017B-31Mar2018-v1/190815_145730/0000/roottree_leptons_1.root");
    TTree *tt=tc;
    TTree *ttc=tt->CloneTree(0);
    return ttc;
}

