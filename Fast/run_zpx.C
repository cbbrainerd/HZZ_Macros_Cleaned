#if !defined(__CINT__) || defined(__MAKECINT__)

#include "ZpXanalyzer.h"
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include "TChain.h"
#include <fstream>

#endif

//Arguments, in order: mc_or_data (just the string "mc" or "data"), dataset (e.g. DoubleEG or ZpBaryonic_MZp_100_MChi_1_13TeV_amcatnlo), year/configuration (e.g. 2017 for data, Fall17 for mc), and list of files to run over
bool good_lumi(int,int) { return true; }

int main (int argc, char ** argv){
    if(argc != 5) {
      std::cout << "Usage: ./RunReference mc_or_data dataset year/config list_of_files" << std::endl;
      exit(2);
    }
    std::cout << "Ciao" << std::endl;
    std::string mc_or_data=argv[1];
    bool is_mc;
    if(mc_or_data=="mc") {
      is_mc=true;
    } else if (mc_or_data=="data") {
      is_mc=false;
    } else {
      std::cout << "The first argument should be exactly \"mc\" or \"data\"\n";
    }
    std::string dataset=argv[2];
    std::string configuration=argv[3];
    //Takes in list of file names to run over
    std::ifstream files;
    files.open(argv[4]);
    if(!files.good()) {
      std::cout << "Failed to open file \"" << argv[4] << "\".\n";
      exit(3);
    }
    std::vector<std::string> filenames;
    std::string line;
    while(std::getline(files,line)) { 
      filenames.push_back(line);
      std::cout << line;
    }
    std::cout << "Running over " << (is_mc ? "MC" : "DATA") << " dataset \"" << dataset << "\" with configuration \"" << configuration << "\"." << std::endl;
   
    TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
    for(auto& fn : filenames) {
      //Remove carriage returns that may have snuck in from windows
      while(!fn.empty() && fn.back()=='\r') fn.pop_back();
      if(fn.empty()) continue;
      chain->Add(fn.c_str());
      std::cout << "Adding " << fn << ".\n";
    }
    //TFile *tf=TFile::Open(filenames[0].c_str());
    //TTree *chain=(TTree*)(tf->Get("HZZ4LeptonsAnalysis"));
    //if(chain==0) exit(5);
    //if (site.find("FNAL")<5){
    //  file3 = new  TDCacheFile (nome,"READ","ROOT file with trees",0);
    //}
    //else {
    //  file3 = TFile::Open(nome);
    //}
  
    //cout << "Read file with name: " << nome << std::endl;
    //TTree *tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");
    
    //HZZ4LeptonsAnalysis make3(tree3,1.,dataconf,mcconf); replace with next line
    //Does this work right?
    //Need to add weights
    //HZZ4LeptonsAnalysis make3(chain,1.,(is_mc?"NO":configuration),(is_mc?configuration:"NO"));
  
    char *path=NULL;
    size_t size=300;
    path=getcwd(path,size);
  
    std::string nome;
    nome+=path;nome+="/output_";
    nome+=dataset;
    nome+=".root";
    TFile *outfile=TFile::Open(nome.c_str(),"RECREATE");
    ZpXanalyzer(chain,outfile).Loop();
    outfile->Close();
    //delete tree3;
    //file3 -> Close();
    
    return 0; 
}
