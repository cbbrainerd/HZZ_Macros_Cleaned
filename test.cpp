#include "TChain.h"
#include <iostream>

int main() {
    TChain *tc=new TChain("HZZ4LeptonsAnalysis","");
    tc->Add("roottree_leptons_754.root");
    TTree *tt=tc;
    TTree *ttc=tt->CloneTree(0);
}

