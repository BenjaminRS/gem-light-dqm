#include "TApplication.h"
#include "TProof.h"
#include "TList.h"
#include "TChain.h"
#include "TError.h"
#include "TFile.h"

#include"db_interface.cxx"

#include <iostream>

#define DEBUG 1

using namespace std;

int main(int argc, char** argv)
{
  //gErrorIgnoreLevel = kError;
  gProofDebugMask = TProofDebug::kAll;
  gProofDebugLevel = 5;
//  gProofDebugLevel = 1;
  std::cout << "--==DQM Main==--" << endl;
  if (argc<2)
    {
      cout << "Please provide input filenames" << endl;
      cout << "Usage: <path>/rundqm inputFile.root" << endl;
      return 0;
    }
  string ifilename = argv[1];
  string ofilename;
  string RunName;
  bool print_hist = false;
  if (argc > 2){
    string option = argv[2];
    if ((option == "--with-print") || (option == "-p")) {
      print_hist = true;
      std::cout << "[MAIN]: will print histograms" << std::endl;
    }
  }

  //if(ifilename.find("chunk") == string::npos) {
  //  RunName = ifilename.substr(ifilename.find("run"),ifilename.find(".raw.root")-ifilename.find("run"));
  //}
  //else {
  //  RunName = ifilename.substr(ifilename.find("run"),ifilename.find("_chunk")-ifilename.find("run"));
  //}
  std::string tmp = ifilename.substr(ifilename.size()-9, ifilename.size());
  if (tmp != ".raw.root") throw std::runtime_error("Wrong input filename (should end with '.raw.root'): "+ifilename);
  ofilename = ifilename.substr(0,ifilename.size()-9);
  ofilename += ".analyzed.root";

  gSystem->Load("libEvent.so");
  TChain *ch = new TChain ("GEMtree","chain");
  ch->Add(ifilename.c_str());

  ofilename = "BRSTest.root";
  TProof::Open("workers=4",0);
//  TProof::Open("workers=4","","",1); //For extra debugging
  gProof->AddInput(new TNamed("PROOF_OUTPUTFILE", ofilename.c_str()));
  gProof->GetInputList()->Print();
  ch->SetProof();;
  string sel;
  string path;
  sel = getenv("BUILD_HOME");
  path = sel+"/gem-light-dqm/proof_packages/gemevent.par";
  gProof->UploadPackage(path.c_str());
  gProof->ShowPackages(kTRUE);
  gProof->EnablePackage("gemevent");
  gProof->ShowPackages(kTRUE);

  path = sel+"/gem-light-dqm/dqm-root/include";
  gProof->AddIncludePath(path.c_str());
  path = sel+"/gem-light-dqm/dqm-root/src/common";
  gProof->AddIncludePath(path.c_str());
  path = sel+"/gem-light-dqm/gemtreewriter/include";
  gProof->AddIncludePath(path.c_str());

  gSystem->SetIncludePath("-I$BUILD_HOME/gem-light-dqm/dqm-root/include -I$BUILD_HOME/gem-light-dqm/dqm-root/src/common -I$BUILD_HOME/gem-light-dqm/gemtreewriter/include");
//  std::cout << "##Running gemTreeReader##" << endl;
//  sel += "/gem-light-dqm/dqm-root/src/common/gemTreeReader.cxx++g";
//  ch->Process(sel.c_str());
//  std::cout << "##Finished gemTreeReader##" << endl;
//  gSystem->SysLog(Error);
  std::cout << "##Running gemTreeTranslator##" << endl;
  sel += "/gem-light-dqm/dqm-root/src/common/gemTreeTranslator.cxx++g";
  //ch->Process(sel.c_str());   
  ch->Process(sel.c_str(),"",10000);   
  std::cout << "##Finished gemTreeTranslator##" << endl;

//  std::cout << "##Running gemTreeTranslator Lite##" << endl;
//  ch->Process(getenv("BUILD_HOME")+"/gem-light-dqm/dqm-root/src/common/gemTreeTranslator.cxx++g","",100);
//  std::cout << "##Finished gemTreeTranslator Lite##" << endl;


  //===================================================
  //    SPLIT OUTPUT INTO MANY TFILES BASED ON NUMBER OF TTREES
  //
  //===================================================
    TFile *inputfile = new TFile(ofilename.c_str(),"READ");
    for (auto keys : *inputfile->GetListOfKeys()) {
        //cout << keys->GetName() << endl;
        TFile *temp_file = new TFile(TString(keys->GetName())+".root","RECREATE");
        TTree *temp_tree = (TTree*)inputfile->Get(keys->GetName());
        temp_file->cd();
        temp_tree->CloneTree()->Write();
        temp_file->Close();
        delete temp_tree;
        temp_tree = 0;
        delete temp_file;
    }
    inputfile->Close();
  //===================================================


  //TProofLite::Mgr("__lite__")->GetSessionLogs()->Display("*");

  if (DEBUG) std::cout << "DQM Complete." << endl;
  return 0;
}
