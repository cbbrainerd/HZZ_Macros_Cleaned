#if !defined(__CINT__) || defined(__MAKECINT__)

#include "HZZ4LeptonsAnalysis_all.h"
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
#include <iostream>
#include <TSystem.h>
#include <TH2.h>
#include "TChain.h"
#include <stdlib.h>
#include <TDCacheFile.h>

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
  /*if (filenames.size() > 1) {
    std::cout << "Temporarily removed TChain support due to problem with TTree::CloneTree when run on TChains. Runs over one file only.\n";
    exit(4);
  }*/
  std::cout << "Running over " << (is_mc ? "MC" : "DATA") << " dataset \"" << dataset << "\" with configuration \"" << configuration << "\"." << std::endl;
 
  //const int ndata = atoi(argv[6]);
  //string datasamples[ndata];
  
  // 
  
  float lumifb=4710./1000.;
  
  Double_t mH=150.;
  std::cout << "mH= " << mH << std::endl;
  
  
  // data
  bool runondata=true;
  
  if (runondata==true){
    
    bool itera=false;
    
    //if (mH<=160.){
    //  if (i==0) itera=true;  
    //}
    //else {
    itera=true;
    //}
    
    //for(int i=0;i<ndata && itera==true;i++){
      
      //string name= "roottree_leptons_crab_";
      //name.append(datasamples[i]);
      //name.append(".root");
      
//      cout << "Input directory is:" << dirInput<< std::endl;

//      string nome_tmp=dirInput.Data();
//      nome_tmp.append("/");
//      nome_tmp.append(name.c_str());
      //sprintf(nome,"%s/%s",dirInput.Data(),name.c_str());
//      Char_t *nome = new Char_t[nome_tmp.length()+1];
//     strcpy(nome,nome_tmp.c_str());    

//      TFile *file3;
//      file3 = TFile::Open(nome);
//      if(!file3) { 
//        std::cout << "Failed to open file \"" << nome << "\". Exiting." << std::endl; 
//        exit(1);
//      }
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
      HZZ4LeptonsAnalysis make3(chain,1.,(is_mc?"NO":configuration),(is_mc?configuration:"NO"));

      char *path=NULL;
      size_t size=300;
      path=getcwd(path,size);

      //char nome[300];
      std::string nome;
      nome+=path;nome+="/output_";
      nome+=dataset;
      //nome+="DoubleMuon_Run2017B-31Mar2018-v1";
      nome+=".root";
      //sprintf(nome,"%s/output_%s.root",path,dataset.c_str())//datasamples[i].c_str());
      make3.Loop(nome.c_str());
      
      std::cout << "Create file with name: " << nome << std::endl;
      //delete tree3;
      //file3 -> Close();
      
    //}
    
  }
  
  return 0; 

}

