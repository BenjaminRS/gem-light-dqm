#ifndef DEBUG
#define DEBUG 1
#endif
#include <TChain.h>
#include <TProofOutputFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <TSystem.h>
#include <bitset>
#include <utility> 

#include "AMC13_histogram.cxx"

unsigned int countBits(uint64_t n){
    unsigned int count = 0;
        while (n){
            count += n & 1;
            n >>= 1;
        }
    return count;
}

class gemTreeTranslator: public TSelector {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain
        TTreeReader     fReader;  //!the tree reader

        // Declaration of leaf types
        Event           *GEMEvents;
        // List of branches
        TBranch        *b_GEMEvents;   //!
        gemTreeTranslator(TTree * /*tree*/ =0) : fChain(0) { }
        //gemTreeTranslator(TTree * /*tree*/ =0) : { }
        virtual ~gemTreeTranslator() { }
        virtual Int_t   Version() const { return 2; }
        virtual void    Begin(TTree *tree);
        virtual void    SlaveBegin(TTree *tree);
        virtual void    Init(TTree *tree);
        virtual Bool_t  Notify();
        virtual Bool_t  Process(Long64_t entry);
        virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
        virtual void    SetOption(const char *option) { fOption = option; }
        virtual void    SetObject(TObject *obj) { fObject = obj; }
        virtual void    SetInputList(TList *input) { fInput = input; }
        virtual TList  *GetOutputList() const { return fOutput; }
        virtual void    SlaveTerminate();
        virtual void    Terminate();
    private :
        vector<AMC13Event> v_amc13;    ///<Vector of AMC13Event
        vector<AMCdata> v_amc;         ///<Vector of AMCdata
        vector<GEBdata> v_geb;         ///<Vector of GEBdata
        vector<VFATdata> v_vfat;       ///Vector of VFATdata

        TTree* treeTranslator; //Should make this private
        std::vector<TTree*> tree_list;
        //int amcID;
        //int amcBID;
        //int gebID;
        //int infoAMCVEC;
        //int infoGEBVEC;
        //int infoVFATVEC;

//        int m_calPhase; //- Not requested
        //int Dly;
//        int m_l1aTime; //- Not requested
        int m_Latency;
        int m_link;
        //int pDel;
        //int mode;
        //int mspl;
        int m_Nev;
        int m_Nhits;
        //int trimDAC;
        //int utime;
        //int vcal;
        uint64_t m_vfatCH_L;
        uint64_t m_vfatCH_M;
        std::pair<uint64_t,uint64_t> m_vfatCH;
//        std::map
        int m_vfatID;
        int m_vfatN;
        int m_vth;
        int m_vth1;
        //int vth2;
        //float ztrim;
        int m_RunType;
        int m_shelf;
        int m_slot;

        //int VFATMap[12][12][24];
        //int m_deltaV;
        //long long int m_OrbitNumber;
        //long long int m_RelOrbitNumber;
        //TDirectory* onlineHistsDir;

        TFile *fFile;
        TProofOutputFile *fProofFile;
        ClassDef(gemTreeTranslator,2);
};

void gemTreeTranslator::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    if (DEBUG) std::cout << "MASTER BEGIN"<< std::endl;
    TString option = GetOption();
}

void gemTreeTranslator::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    if (DEBUG) std::cout << "SLAVE BEGIN"<< std::endl;
    TString option = GetOption();
    gSystem->Load("libEvent.so");

    // The file for merging
    fProofFile = new TProofOutputFile("SimpleFile.root", TProofOutputFile::kMerge);
    TNamed *out = (TNamed *) fInput->FindObject("PROOF_OUTPUTFILE"); //analyzed.root
    if (out) fProofFile->SetOutputFileName(out->GetTitle());
    TDirectory *savedir = gDirectory;
    fFile = fProofFile->OpenFile("RECREATE");
    if (fFile && fFile->IsZombie()) SafeDelete(fFile);
    savedir->cd();

    // Cannot continue
    if (!fFile) {
        TString amsg = TString::Format("ProofSimpleFile::SlaveBegin: could not create '%s':"
                " instance is invalid!", fProofFile->GetName());
        Abort(amsg, kAbortProcess);
        return;
    }
    fFile->cd();

    // BRS: Lets make the Tfile here:
    int numOfAMC=3;
    int numOfGEB=5;
    int numOfVFAT=24;

    tree_list.clear(); // clear the vector having list of trees
    for (auto amc_count=0; amc_count<numOfAMC; ++amc_count){
        for (auto geb_count=0; geb_count<numOfGEB; ++geb_count){
            TString treeName="gemTree_"+std::to_string(amc_count)+"_"+std::to_string(geb_count);
            treeTranslator    = new TTree(treeName,"Tree holding gem info");
//            auto branchcalPhase  = treeTranslator->Branch("calPhase",&m_calPhase,"calPhase/I"); //calPhase - Not requested
            //auto branchDly  = treeTranslator->Branch("Dly",&Dly,"Dly/I"); //amcID
//            auto branchl1aTime  = treeTranslator->Branch("l1aTime",&m_l1aTime,"l1aTime/I"); //l1aTime - Not requested
            auto branchlatency  = treeTranslator->Branch("latency",&m_Latency,"latency/I"); //amcID
            auto branchlink  = treeTranslator->Branch("link",&m_link,"link/I"); //amcID
            //auto branchpDel  = treeTranslator->Branch("pDel",&pDel,"pDel/I"); //amcID
            //auto branchmode  = treeTranslator->Branch("mode",&mode,"mode/I"); //amcID
            //auto branchmspl  = treeTranslator->Branch("mspl",&mspl,"mspl/I"); //amcID
            auto branchNev  = treeTranslator->Branch("Nev",&m_Nev,"Nev/I"); //amcID
            auto branchNhits  = treeTranslator->Branch("Nhits",&m_Nhits,"Nhits/I"); //amcID
            //auto branchtrimDAC  = treeTranslator->Branch("trimDAC",&trimDAC,"trimDAC/I"); //amcID
            //auto branchutime  = treeTranslator->Branch("utime",&utime,"utime/I"); //amcID
            //auto branchvcal  = treeTranslator->Branch("vcal",&vcal,"vcal/I"); //amcID
            auto branchvfatCH  = treeTranslator->Branch("vfatCH",&m_vfatCH,"vfatCH/I"); //amcID
            auto branchvfatID  = treeTranslator->Branch("vfatID",&m_vfatID,"vfatID/I"); //amcID
            auto branchvfatN  = treeTranslator->Branch("vfatN",&m_vfatN,"vfatN/I"); //amcID
            auto branchvth  = treeTranslator->Branch("vth",&m_vth,"vth/I"); //amcID
            auto branchvth1  = treeTranslator->Branch("vth1",&m_vth1,"vth1/I"); //amcID
            //auto branchvth2  = treeTranslator->Branch("vth2",&vth2,"vth2/I"); //amcID
            //auto branchztrim  = treeTranslator->Branch("ztrim",&ztrim,"ztrim/F"); //amcID
            auto branchRType = treeTranslator->Branch("m_RunType",&m_RunType,"m_RunType/I"); //Run Type
            auto branchshelf  = treeTranslator->Branch("shelf",&m_shelf,"shelf/I"); //amcID
            auto branchslot  = treeTranslator->Branch("slot",&m_slot,"slot/I"); //amcID
            tree_list.push_back(treeTranslator);
        } // End of GEB loop
    }  // End of AMC loop

    //fFile->Add(treeTranslator);
    for (auto obj : tree_list) {
        fFile->Add(obj);
    }

    if (fChain){
        fChain->GetEntry(0);
        v_amc13 = GEMEvents->amc13s();
        if (DEBUG) std::cout << "SlaveBegin: v_amc13 = "<< v_amc13.size() << std::endl;
    }
    else{
        Error("Begin", "no fChain!");
    }
    ///// AMCdata::BID() or AMCdata::AMCNum() //AMC ID..
    ///// VFATdata::ChipID() vfatID
    ///// shelf and slot?

    if (DEBUG) std::cout << "SLAVE END"<< std::endl;
}

//!Fills the histograms that were book from bookAllHistograms
Bool_t gemTreeTranslator::Process(Long64_t entry)
{
    //fReader.SetEntry(entry);
    fChain->GetEntry(entry);
    if (DEBUG) std::cout<<"########################################"<<std::endl;
    if (DEBUG) std::cout<<"\n\tEvent loop\n"<<std::endl;
    if (DEBUG) std::cout<<"########################################"<<std::endl;
    int amc13_count=0;    //counter through AMC13s
    int amc_count=0;      //counter through AMCs
    int geb_count=0;      //counter through GEBs
    int vfat_count=0;      //counter through VFATs

    /* local variables that will be filled to the tree  */
//    int local_calPhase = -1; //Not requested
//    int local_l1aTime = -1; //- Not requested
    int local_link = -1;
    int local_Nev = 0;
    int local_Nhits = 0;
    uint64_t local_vfatCH_L = -1;
    uint64_t local_vfatCH_M = -1;
    std::pair<uint64_t,uint64_t> local_vfatCH(local_vfatCH_M,local_vfatCH_L);
    int local_vfatID = -1;
    int local_vfatN = -1;
    int local_vth = -1;
    int local_vth1 = -1;
    int local_shelf = -1;
    int local_slot = -1;
    /* END local variable*/

    v_amc13 = GEMEvents->amc13s();
    if (DEBUG) cout << "Found "<< v_amc13.size() <<" AMC13s" << endl;
    /* LOOP THROUGH AMC13s */
    for(auto amc13 = v_amc13.begin(); amc13!=v_amc13.end(); amc13++) {
        v_amc = amc13->amcs();
        if (DEBUG) cout << "Found "<<v_amc.size() << " AMCs" << endl;
        local_shelf=amc13->Source_id();
        /* LOOP THROUGH AMCs */
        for(auto amc=v_amc.begin(); amc!=v_amc.end(); amc++){
            v_geb = amc->gebs();
            if (DEBUG) cout << "Found "<< v_geb.size()<<" GEBs "<< endl;
            local_slot =amc->AMCnum();
//            local_l1aTime = amc->L1AT();
            //-- EXTRACT RUN_PARAMS --
            m_RunType=amc->Rtype();
            int runParams = ((amc->Param1() & 0xff) << 16) | ((amc->Param2() & 0xff) << 8) | amc->Param3();
            if (m_RunType == 0x2) m_Latency = runParams & 0x3ff;
            else if(m_RunType == 0x3) m_Latency = (runParams >> 13) & 0xff;
            else if(m_RunType == 0x4) m_Latency = (runParams >> 13) & 0xff;
    
            if (local_slot == 2){
                int l1ANum = amc->L1A();
                int pulseStretch = (runParams >> 10) & 0x7;
                int isCurrent = (runParams >> 21) & 0x1;
                if (m_RunType == 0x2){
                    int isInt = (runParams >> 22) & 0x1;
                    int calDac = (runParams >> 13) & 0xff;
                    if (DEBUG) std::cout<<"AMC" << local_slot << " L1Anum " << l1ANum <<" RType 0x"<< std::hex <<  m_RunType 
                                        <<": RunParams = 0x" << std::hex << runParams
                                        << ": Decoded (isInt | isC | CalDac | P.S. | Lat) = (" 
                                        << std::dec << isInt << " | " << isCurrent << " | " << calDac << " | " << pulseStretch << " | " << m_Latency << " )" << std::endl;
                }
                else if(m_RunType == 0x3){
                    int selCompMode = (runParams >> 21) & 0x3;
                    int thrZccDac = runParams & 0xff;
                    if (DEBUG) std::cout<<"AMC" << local_slot << " L1Anum " << l1ANum <<" RType 0x"<< std::hex <<  m_RunType
                                        <<": RunParams = 0x" << std::hex << runParams
                                        << ": Decoded ( selComp | armDac | P.S. | zccDac) = ("
                                        << std::dec << selCompMode << " | " << m_Latency << " | " << pulseStretch << " | " << thrZccDac << " )" << std::endl;
                }
                else if(m_RunType == 0x4){
                    int thrArmDac = runParams & 0xff;
                    if (DEBUG) std::cout<<"AMC" << local_slot << " L1Anum " << l1ANum <<" RType 0x"<< std::hex <<  m_RunType
                                        <<": RunParams = 0x" << std::hex << runParams
                                        << ": Decoded ( isC | CalDac | P.S. | armDac) = ("
                                        << std::dec << isCurrent << " | " << m_Latency << " | " << pulseStretch << " | " << thrArmDac << " )" << std::endl;
                }
            }
           //-- END EXTRACT RUN_PARAMS --

            geb_count=0;
            //if (DEBUG) std::cout<<"amc:"<< amc_count << ": local_slot = "<<local_slot<<std::endl;
//            if (DEBUG) std::cout<<"amc:"<< amc_count << ": local_l1aTime = "<<local_l1aTime<<std::endl;
            if (DEBUG) std::cout<<"amc:"<< amc_count << ": m_RunType= "<<m_RunType<<std::endl;
            if (DEBUG) std::cout<<"amc:"<< amc_count << ": m_Latency = "<<m_Latency<<std::endl;
            //if (DEBUG) std::cout<<"amc:"<< amc_count << ": local_vth = "<<local_vth<<std::endl;
            //if (DEBUG) std::cout<<"amc:"<< amc_count << ": local_vth1 = "<<local_vth1<<std::endl;
            //if (DEBUG) std::cout<<"amc:  = "<<<<std::endl;
            /* LOOP THROUGH GEBs */
            for(auto g=v_geb.begin(); g!=v_geb.end();g++) {
                v_vfat = g->vfats();
                local_link = g->InputID();
                vfat_count = 0;
                if (DEBUG) std::cout<<"-----"<<endl;
                if (DEBUG) std::cout<<"GEB:" << geb_count << ": local_link = "<<local_link<<std::endl;
                if (DEBUG) std::cout << "\tFound "<< v_vfat.size()<< " VFATs" << std::endl;
                /* LOOP THROUGH VFATs */
                for(auto v=v_vfat.begin(); v!=v_vfat.end();v++) {
                    //if (DEBUG) std::cout << "Found VFATs " << std::endl;
                    TString treeName="gemTree_"+std::to_string(amc_count)+"_"+std::to_string(geb_count);
                    //local_amcID = amc->AMCnum();
//                    local_calPhase = v->isBlockGood();
                    //local_vfatID = v->chipID();
                    local_vfatN = v->Pos();
//                    local_Nev  = v->EC(); //Not correct
                    local_vfatCH_L = v->lsData();
                    local_vfatCH_M = v->msData();
                    local_vfatCH.first = v->msData();
                    local_vfatCH.second = v->lsData();
                    local_Nhits += countBits(local_vfatCH_L) + countBits(local_vfatCH_M);
//                    if (DEBUG) std::cout<< "\t==========================" << std::endl;
//                    if (DEBUG) std::cout<<"\tVFAT:" << vfat_count << ": treeName = "<<treeName<<std::endl;
////                    if (DEBUG) std::cout<<"\tVFAT:" << vfat_count << ": local_amcID = "<<local_amcID<<std::endl;
////                    if (DEBUG) std::cout<<"\tVFAT:" << vfat_count << ": local_calPhase = "<<local_calPhase<<std::endl;
//                    if (DEBUG) std::cout<<"\tVFAT:" << vfat_count << ": local_Pos = "<<local_vfatN<<std::endl;
//                    if (DEBUG) std::cout<<"\tVFAT:" << vfat_count << ": local_SlotNumber = "<<v->SlotNumber()<<std::endl;
////                    if (DEBUG) std::cout<<"\tVFAT:" << vfat_count << ": countBits(binChannel_0_63) = "<<countBits(binChannel_0_63)<<std::endl;
////                    if (DEBUG) std::cout<<"\tVFAT:" << vfat_count << ": countBits(binChannel_64_127) = "<<countBits(binChannel_64_127)<<std::endl;
//                    //if (binChannel_0_63 != 0 || binChannel_64_127 != 0) {
//                    if (DEBUG) std::cout<<"\tVFAT:" << vfat_count << ": local_Pos = "<<local_vfatN<<std::endl;
//                    if (DEBUG) std::cout<<"\tVFAT:" << vfat_count <<  ": lsData = " << local_vfatCH_L << std::endl;
//                    if (DEBUG) std::cout<<"\tVFAT:" << vfat_count <<  ": msData = " << local_vfatCH_M << std::endl;
//                    if (DEBUG) std::cout<<"\tVFAT:" << vfat_count <<  ": local_Nhits = " << local_Nhits << std::endl;
//                    //}
//                    if (DEBUG) std::cout<< "\t==========================" << std::endl;

                    for (auto obj : tree_list) {
                        if (obj->GetName() == treeName) {
                            static_cast<TTree*>(obj)->Fill();
                            //if (DEBUG) cout << "Object name : " << obj->GetName() << endl;
                        }
                    }
                    //
                    vfat_count++;
                } /* END VFAT LOOP */
                ++geb_count;
            } /* END GEB LOOP */
            amc_count++;
        } /* END AMC LOOP */
        amc13_count++;
    } /* END AMC13 LOOP */
    /*  Save Branches */
//    m_calPhase = local_calPhase;
//    m_l1aTime = local_l1aTime;
    m_link = local_link;
    m_Nev = local_Nev;
    m_Nhits = local_Nhits;
    m_vfatCH_L=local_vfatCH_L;
    m_vfatCH_M=local_vfatCH_M;
    //vfatID = local_vfatID;
    m_vfatN = local_vfatN;
    m_vth = local_vth;
    m_vth1 = local_vth1;
    m_shelf = local_shelf;
    m_slot = local_slot;

   if (DEBUG) std::cout<< "==========================" << std::endl;
   if (DEBUG) std::cout<<"m_Latency:" << m_Latency <<std::endl;
   if (DEBUG) std::cout<<"m_link:" << m_link <<std::endl;
   if (DEBUG) std::cout<<"m_Nev:" << m_Nev <<std::endl;
   if (DEBUG) std::cout<<"m_Nhits:" << m_Nhits <<std::endl;
   if (DEBUG) std::cout<<"m_vfatCH_L:" << m_vfatCH_L <<std::endl;
   if (DEBUG) std::cout<<"m_vfatCH_M:" << m_vfatCH_M <<std::endl;
   if (DEBUG) std::cout<<"m_vfatN:" << m_vfatN <<std::endl;
   if (DEBUG) std::cout<<"m_vth:" << m_vth <<std::endl;
   if (DEBUG) std::cout<<"m_vth1:" << m_vth1 <<std::endl;
   if (DEBUG) std::cout<<"m_shelf:" << m_shelf <<std::endl;
   if (DEBUG) std::cout<<"m_slot:" << m_slot <<std::endl;
   if (DEBUG) std::cout<< "==========================" << std::endl;
//   if (DEBUG) std::cout<<":" <<  <<std::endl;

    //gebID = local_gebID;
    //amcID = local_amcID;

    treeTranslator->Fill();
    return kTRUE;
}

void gemTreeTranslator::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.
    //
    // Should we fill summary canvases here?
    //
    //for (unsigned int i = 0; i < 12; i++){
    //    if (auto a = m_amc13H->amcsH(i)){
    //        for (unsigned int j = 0; j < 12; j++){
    //            if (auto g = a->gebsH(j)){
    //                g->fillSummaryCanvases();
    //                for (unsigned int k = 0; k < 24; k++){
    //                    if (auto v = g->vfatsH(k)){
    //                        v->fillWarnings();
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}
    TDirectory *savedir = gDirectory;
    fFile->cd();
    fFile->Write();
    fOutput->Add(fProofFile);
    fFile->Close();
}
void gemTreeTranslator::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.
    Printf("Testing if gemTree was created successfully ....");
    TTree* gemTree = dynamic_cast<TTree*>(fOutput->FindObject("gemTree"));
    if (gemTree){
        Printf("gemTree exists!");
    }

    //  TTree* gemTreeB = dynamic_cast<TTree*>(fProofFile->FindObject("gemTree"));
    //if (treeTranslator){
    //Printf("gemTree exists! in memory at Terminate");
    //   treeTranslator->Print();
    //   treeTranslator->Show(0);
    //}
    //for (auto obj : tree_list)
    for(unsigned int i=0; i<tree_list.size(); i++) {
        treeTranslator = (TTree*)tree_list.at(i);
        treeTranslator->Write("", TObject::kOverwrite);
        delete treeTranslator;
        treeTranslator = 0;
    }

    if ((fProofFile = dynamic_cast<TProofOutputFile*>(fOutput->FindObject("SimpleFile.root")))) {
        TString outputFile(fProofFile->GetOutputFileName());
        TString outputName(fProofFile->GetName());
        Printf("outputFile: %s", outputFile.Data());
    }
}
//int gemTreeTranslator::slotFromMap(int a, int g, int cid)
//{
//  int res = -1;
//  for (int i=0; i<24; i++){
//    if (VFATMap[a][g][i] == cid) {res = i;}
//  }
//  return res;
//}
void gemTreeTranslator::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fChain->SetMakeClass(0);

    GEMEvents = new Event();

    fChain->SetBranchAddress("GEMEvents", &GEMEvents, &b_GEMEvents);
    v_amc13 = GEMEvents->amc13s();
    if (DEBUG) std::cout << "v_amc13 size = " << v_amc13.size() << std::endl;
    m_RunType=-1;
    m_Latency=-1;
//    m_calPhase=-1; //- Not requested
//    m_l1aTime=-1; //- Not requested
    m_Latency=-1;
    m_link=-1;
    m_Nev=-1;
    m_Nhits=-1;
    m_vfatCH_L=-1;
    m_vfatCH_M=-1;
    m_vfatID=-1;
    m_vfatN=-1;
    m_vth=-1;
    m_vth1=-1;
    m_RunType=-1;
    m_shelf=-1;
    m_slot=-1;
}

Bool_t gemTreeTranslator::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}
