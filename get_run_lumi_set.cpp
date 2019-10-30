#include <iostream>
#include <sstream>
#include <string>


template<class T>
std::string join(std::set<T> x) {
    std::ostringstream out;
    T first;
    T last;
    bool init=true;
    for(const auto& element : x) {
        if(init) {
            first=element;
            last=element;
            init=false;
            continue;
        }
        if(element==last+1) {
            last=element;
        } else {
            if(first==last) out << first << ",";
            else out << first << "-" << last << ",";
            first=element;
            last=element;
        }
    }
    if(last==first) out << first;
    else out << first << "-" << last;
    return out.str();
}

int get_lumis_for_file(const char* fn,const char* treename) {
    TFile* tf=TFile::Open(fn);
    std::unordered_map<UInt_t,std::set<UInt_t> > run_lumis;
    auto tree=(TTree*)(tf->Get(treename));
    if(!tree) return 1;
    auto trun=tree->GetBranch("Run");
    UInt_t nrun;
    auto tlumi=tree->GetBranch("LumiSection");
    UInt_t nlumi;
    for(std::size_t i=0;i<tree->GetEntriesFast();++i) {
        tree->GetEvent(i);
        std::set<UInt_t>& lumis=run_lumis.insert(std::make_pair(nrun,std::set<UInt_t>())).first->second;
        lumis.insert(nlumi);
    }
    std::cout << fn << ":\n";
    for(auto& element : run_lumis) {
            std::cout << element.first << ":" <<  join(element.second) << ";";
    }
    std::cout << "\n";
    return 0;
}
