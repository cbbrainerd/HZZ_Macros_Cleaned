#include "TTree.h"
#include <array>

template <class T>
using array=std::array<T,3>;

enum bkg_type : unsigned char {
    b_3mu=0,
    b_1e2mu=1,
    b_2e1mu=2,
    b_3e=3
};

//Fake rate, OS method
//Twiki versus analysis note unclear: is 20/10 any two leptons or the Z candidate leptons? 3 leptons in event or 2 tight+exactly 1 loose? etc.
class ZpXOutputTree {
    TTree*          ttree;

    UInt_t          run;
    Uint_t          lumi;
    UInt_t          event;

    float           pt[3];
    float           eta[3];
    float           phi[3];
    int             charge1;
    int             is_muon1;

    float           pt2;
    float           eta2;
    float           phi2;
    int             charge2;
    int             is_muon2;

    float           pt2;
    float           eta2;
    float           phi2;
    int             charge2;
    int             is_muon2;

    float Z_mass;
    float Z_pt;
    float Z_eta;

    float MET_T1;
    float MET_T1xy;
    float MET_puppi;
    float mT;

    bool smart_cut; //Smart cut passed or not

    bool initialized;

    TTree* initialize_out(TFile* tfile);
    TTree* initialize_out(const char* filename);
}
