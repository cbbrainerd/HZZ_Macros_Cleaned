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
#include <utility>
#if product_4mu

#elif product_2e2mu

#elif product_4e

#else
#error "Product must be one of \"4mu\", \"2e2mu\", and \"4e\".\n" 
#endif

bool good_lumi(int run,int lumi);

#include "pileup_corrector.h"
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
// Header file for the classes stored in the TTree if any. Separate this from the rest of the logic so that further changes to the Ntuple structure can be handled more easily
#include "NewNtuple.h"
// Fixed size dimensions of array or collections stored in the TTree if any.
class HZZ4LeptonsAnalysis : public NewNtuple {
public:
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
   bool isMC;
   int year;
   void common_loop();
   pileup_corrector pileup_corr;
   Float_t         RECO_PFMET_xycorr;
   Float_t         RECO_PFMET_PHI_xycorr;
};

#endif

#ifdef HZZ4LeptonsAnalysis_cxx

HZZ4LeptonsAnalysis::common_loop() {
    //Fix
    while(RECOELE_isGap->size() < RECOELE_PT->size()) RECOELE_isGap->push_back(false);
}

HZZ4LeptonsAnalysis::HZZ4LeptonsAnalysis(TTree *tree,Double_t weight_, std::string DATA_type_, std::string MC_type_) : NewNtuple(tree), pileup_corr(MC_type_!="NO",(MC_type_=="NO")?MC_type:DATA_type_)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   weight = weight_;
   DATA_type = DATA_type_;
   MC_type = MC_type_;
   if(!((DATA_type=="NO")^((MC_type=="NO")))) {
     std::cout << "Invalid types: DATA_type: " << DATA_type << " MC_type: " << MC_type << ". One of these should be \"NO\".\n";
     exit(1);
   }
   isMC=(MC_type!="NO");
   if(isMC) {
    auto len=MC_type.length();
    year=std::stoi(MC_type.substr(len-2))+2000;
   } else {
    year=std::stoi(DATA_type);
   }
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
   std::cout << "Init(TTree*) called.\n";
//Call the (automatically generated) NewNtuple::Init function to setup branches, etc.
   NewNtuple::Init(tree);
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

struct placeholder {
    placeholder()=delete;
    template<typename T>
    operator T&() {};
};

template <class T,class=decltype(std::declval<T&>()=std::declval<int>())>
typename std::vector<T>::reference get_default_helper(T*) {
    static T def=-999;
    if(def!=-999) { std::cout << "Warning! Value of vector changed out of bounds. Probably a bug.\n"; def=-999; }
    return def;
}

template <>
typename std::vector<bool>::reference get_default_helper<bool>(bool*) {
    static std::vector<bool> def={false};
    if(def[0]) { std::cout << "Warning! Value of vector changed out of bounds. Probably a bug.\n"; def[0]=false; }
    return def[0]; 
}

placeholder& get_default_helper(...) {
    std::cerr << "default called on vector type. Nothing is implemented for this, so it's a fatal error.\n";
    exit(1);
}

template <class T>
auto get_default() -> decltype(get_default_helper((T*)nullptr)) { return get_default_helper((T*)nullptr); }

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
