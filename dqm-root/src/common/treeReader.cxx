#ifndef DEBUG
#define DEBUG 1
#endif
#define NVFAT 24
#include <iomanip> 
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TNtuple.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TBenchmark.h>
#include <TInterpreter.h>
#include <TApplication.h>
#include <TString.h>
#include <Event.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TRefArray.h>
#include <TRef.h>
#include <TH1.h>
#include <TBits.h>
#include <TMath.h>
#include <TFile.h>
#include <TClassTable.h>
#include <TTree.h>
#include <TBranch.h>

#include "gem/datachecker/GEMDataChecker.h"
#include "gem/readout/GEMslotContents.h"
#include "plotter.cxx"

/**
* GEM Tree Reader example (reader) application 
*/

/*! 
  \brief GEM Tree reader reads file with the GEM Event Tree and fills the histogram with number of vfats per event
*/

using namespace std;

uint16_t gem::readout::GEMslotContents::slot[NVFAT] = {
      0xfff,0xfff,0xfff,0xfff,0xfff,0xfff,0xfff,0xfff,
      0xfff,0xfff,0xfff,0xfff,0xfff,0xfff,0xfff,0xfff,
      0xfff,0xfff,0xfff,0xfff,0xfff,0xfff,0xfff,0xfff,
    };
bool gem::readout::GEMslotContents::isFileRead = false;

class gemTreeReader {
  public:
    gemTreeReader(const std::string &ifilename)
    {
      std::string tmp = ifilename.substr(ifilename.size()-9, ifilename.size());
      if (tmp != ".raw.root") throw std::runtime_error("Wrong input filename (should end with '.raw.root'): "+ifilename);
      ifile = new TFile(ifilename.c_str(), "READ");
      std::string ofilename = ifilename.substr(0,ifilename.size()-9);
      ofilename += ".analyzed.root";
      ofile = new TFile(ofilename.c_str(), "RECREATE");
      if (DEBUG) std::cout << std::dec << "[gemTreeReader]: File for histograms created, start to book histograms " << std::endl;   
      this->bookHistograms();
    }
    ~gemTreeReader(){}
    void createHistograms()
    {
      this->fillHistograms();
    }
  private:
    TFile *ifile;
    TFile *ofile;

    TH1F* hiVFAT                     [3]; 
    TH1F* hiVFATsn                   [3]; 
    TH1C* hiChip                     [3]; 
    TH1C* hi1010                     [3]; 
    TH1C* hi1100                     [3]; 
    TH1C* hi1110                     [3]; 
    TH1C* hiFlag                     [3]; 
    TH1C* hiCRC                      [3]; 
    TH1C* hiDiffCRC                  [3]; 
    TH1C* hiFake                     [3]; 
    TH1F* hiCh128                    [3];     
    TH2C* hi2DCRC                    [3];     
    TH1C* hiDiffBXandBC              [3];
    TH1C* hiRatioBXandBC             [3];
    TH1C* hiSignal                   [3];
    TH1C* hichfired                  [3];
    TH1C* hichnotfired               [3];
    TH1F* hiCh_notfired              [3];
    TH1F* hiVFATfired_perevent[3][NVFAT];
    TH1F* hiCh128chipFired    [3][NVFAT];

    TDirectory *dir[3];

    int counters_[4]; // [0] - total events
                      // [1] - good events
                      // [2] - bad events 
                      // [3] - good events failed CRC check

    void bookHistograms()
    {
      std::string dirname[3] = {"AllEvents", "GoodEvents", "BadEvents"};
      char name[128], title[500],name1[128], title1[500];
      std::string type[NVFAT] = {"Slot0" , "Slot1" , "Slot2" , "Slot3" , "Slot4" , "Slot5" , "Slot6" , "Slot7", 
                              "Slot8" , "Slot9" , "Slot10", "Slot11", "Slot12", "Slot13", "Slot14", "Slot15", 
                              "Slot16", "Slot17", "Slot18", "Slot19", "Slot20", "Slot21", "Slot22", "Slot23"};
      if (DEBUG) std::cout << std::dec << "[gemTreeReader]: Make directories and create histograms inside" << std::endl;   
      for (int i = 0; i < 3; i++) {
        dir[i] = ofile->mkdir(dirname[i].c_str());
        if (DEBUG) std::cout << std::dec << "[gemTreeReader]: Directory " << i+1 << " created" << std::endl;   
        dir[i]->cd();
        hiVFAT         [i] = new TH1F((dirname[i]+"_VFAT").c_str(), "Number VFAT blocks per event", 24,  0., 24. );
        hiVFATsn       [i] = new TH1F((dirname[i]+"_VFATsn").c_str(), "VFAT slot number", 24,  0., 24. );
        hiDiffBXandBC  [i] = new TH1C((dirname[i]+"_DiffBXandBC").c_str(), "Difference of BX and BC", 100000, 0x0, 0x1869F );
        hiRatioBXandBC [i] = new TH1C((dirname[i]+"_RatioBXandBC").c_str(), "Ratio of BX and BC", 1000, 0x0, 0xa );
        hiChip         [i] = new TH1C((dirname[i]+"_ChipID").c_str(), "ChipID",         4096, 0x0, 0xfff );
        hi1010         [i] = new TH1C((dirname[i]+"_1010").c_str(), "Control Bits 1010", 16, 0x0, 0xf );
        hi1100         [i] = new TH1C((dirname[i]+"_1100").c_str(), "Control Bits 1100", 16, 0x0, 0xf );
        hi1110         [i] = new TH1C((dirname[i]+"_1110").c_str(), "Control Bits 1110", 16, 0x0, 0xf );
        hiFlag         [i] = new TH1C((dirname[i]+"_Flag").c_str()  , "Flag",            16, 0x0, 0xf );
        hiCRC          [i] = new TH1C((dirname[i]+"_CRC").c_str(),     "CRC",             100, 0x0, 0xffff );
        hiDiffCRC      [i] = new TH1C((dirname[i]+"_DiffCRC").c_str(), "CRC_Diff",    100, 0xffff, 0xffff );
        hichfired      [i] = new TH1C((dirname[i]+"_chfired").c_str(), "Channels fired per event",      500, 0., 500. );
        hichnotfired   [i] = new TH1C((dirname[i]+"_chnotfired").c_str(), "Channels not fired per event",      500, 0., 500. );
        hiCh_notfired  [i] = new TH1F((dirname[i]+"_Ch_notfired").c_str(), "Strips",          128, 0., 128. );
        hiCh128        [i] = new TH1F((dirname[i]+"_Ch128").c_str(), "Strips",          128, 0., 128. );
        hi2DCRC        [i] = new TH2C((dirname[i]+"_CRC1_vs_CRC2").c_str(), "CRC1 vs_ RC2", 100, 0x0000, 0xffff, 100, 0x0000, 0xffff);
        if (DEBUG) std::cout << std::dec << "[gemTreeReader]: Main histograms ["<<i<<"] created" << std::endl;   
        for(int j=0; j < NVFAT; j++){
          if (DEBUG) std::cout << std::dec << "[gemTreeReader]: Start 2d array of histograms ["<<i<<"]["<<j<<"] creation" << std::endl;   
          sprintf (name  , (dirname[i]+"_hiVFATfired_perevent_%s").c_str(), type[j].c_str());
          sprintf (title , "VFAT chip %s fired per event", type[j].c_str());
          if (DEBUG) std::cout << std::dec << "[gemTreeReader]: Get name and title for hiVFATfired_perevent["<<i<<"]["<<j<<"]: name : " << name << " title : " << title << std::endl;   
          sprintf (name1 , (dirname[i]+"_hiCh128chipFired_%s").c_str(), type[j].c_str());
          sprintf (title1, "Strips fired for VFAT chip %s", type[j].c_str());
          if (DEBUG) std::cout << std::dec << "[gemTreeReader]: Get name and title for hiCh128chipFired["<<i<<"]["<<j<<"]" << std::endl;   
          hiVFATfired_perevent[i][j] = new TH1F(name, title, 20, 0., 20.);
          hiCh128chipFired    [i][j] = new TH1F(name1, title1, 128, 0., 128.);
          if (DEBUG) std::cout << std::dec << "[gemTreeReader]: 2d array of histograms ["<<i<<"]["<<j<<"] created" << std::endl;   
        }
      }
    }
    void fillHistograms()
    {
      if (DEBUG) std::cout << "[gemTreeReader]: Starting filling the histograms" << std::endl;
      TTree *tree = (TTree*)ifile->Get("GEMtree");
      Event *event = new Event();
      TBranch *branch = tree->GetBranch("GEMEvents");
      branch->SetAddress(&event);
      Int_t nentries = tree->GetEntries();
      if (DEBUG) std::cout << "[gemTreeReader]: Number of entries in the TTree : " << nentries << std::endl;
      // init counters
      for (int c_=0; c_<4; c_++) counters_[c_] = 0;
      // loop over tree entries
      for (Int_t i = 0; i < nentries; i++)
      {
        if (DEBUG) std::cout << "[gemTreeReader]: Processing event " << i << std::endl;
        // clear number of VFATs
        int nVFAT[3] = {0,0,0};
        int firedchannels[3] = {0,0,0};
        int notfiredchannels[3] = {0,0,0};
        int vfatId[3][NVFAT];
        for (int g=0; g < 3; g++){
          for (int l = 0; l<NVFAT; l++){
            vfatId[g][l] = 0;
            if (DEBUG) std::cout << std::dec << "[gemTreeReader]: Fired chip counter when initializing event"  <<  vfatId[g][l]  << std::endl;   
          }
        }
        // Retrieve next entry
        branch->GetEntry(i);
        bool eventIsOK = event->isEventGood();
        // retrieve bunch crossing from evet
        uint32_t BX = event->BXID();
        uint16_t BC;
        // create vector of GEBdata. For data format details look at Event.h
        vector<GEBdata> v_geb;
        v_geb = event->gebs();
        // loop over gebs
        for (Int_t j = 0; j < v_geb.size(); j++)
        {
          // create vector of VFATdata. For data format details look at Event.h
          vector<VFATdata> v_vfat;
          v_vfat = v_geb.at(j).vfats();
          // Increment number of VFATs in the given event
          nVFAT[0] += v_vfat.size();
          if (eventIsOK) { nVFAT[1] += v_vfat.size();} else { nVFAT[2] += v_vfat.size();}
          // loop over vfats
          for (Int_t k = 0; k < v_vfat.size(); k++)
          {
            if (DEBUG) std::cout << std::dec << "[gemTreeReader]: VFAT # "  <<  k << std::endl;   
            if (DEBUG) std::cout << std::dec << "[gemTreeReader]: EC of the vfat inside loop===> "  <<  static_cast<uint32_t>(v_vfat.at(k).EC()) << std::hex << std::endl;   
            if (DEBUG) std::cout << std::dec << "[gemTreeReader]: BC of the vfat inside loop===> "  <<  v_vfat.at(k).BC() << std::hex << std::endl;   
            if (DEBUG) std::cout << std::dec << "[gemTreeReader]: Slot number of responded chip "  <<  v_vfat.at(k).SlotNumber()  << std::endl;   
            if (DEBUG) std::cout << std::dec << "[gemTreeReader]: Fired chip counter before incrementing"  <<  vfatId[v_vfat.at(k).SlotNumber()]  << std::endl;   
            if (v_vfat.at(k).SlotNumber()>(-1)){ 
              vfatId[0][v_vfat.at(k).SlotNumber()]++;
              if (eventIsOK) {vfatId[1][v_vfat.at(k).SlotNumber()]++;} else {vfatId[2][v_vfat.at(k).SlotNumber()]++;}
            }
            if (DEBUG) std::cout << std::dec << "[gemTreeReader]: Fired chip counter "  <<  vfatId[v_vfat.at(k).SlotNumber()]  << std::endl;   
            // fill histograms for all events
            dir[0]->cd();
            this->fillVFATHistograms(&v_vfat.at(k), hiVFATsn[0], hiCh128[0], hiCh_notfired[0], hiChip[0], hi1010[0], hi1100[0], hi1110[0], hiFlag[0], hiCRC[0], hiDiffCRC[0], hi2DCRC[0], hiCh128chipFired[0], firedchannels[0], notfiredchannels[0]);
            // fill histograms for good and bad events
            if (eventIsOK){
              dir[1]->cd();
              this->fillVFATHistograms(&v_vfat.at(k), hiVFATsn[1], hiCh128[1], hiCh_notfired[1], hiChip[1], hi1010[1], hi1100[1], hi1110[1], hiFlag[1], hiCRC[1], hiDiffCRC[1], hi2DCRC[1], hiCh128chipFired[1], firedchannels[1], notfiredchannels[1]);
            } else {
              dir[2]->cd();
              this->fillVFATHistograms(&v_vfat.at(k), hiVFATsn[2], hiCh128[2], hiCh_notfired[2], hiChip[2], hi1010[2], hi1100[2], hi1110[2], hiFlag[2], hiCRC[2], hiDiffCRC[2], hi2DCRC[2], hiCh128chipFired[2], firedchannels[2], notfiredchannels[2]);
            }
            
            BC = v_vfat.at(k).BC();

          }// end of loop over VFATs
        }// end of loop over GEBs
        dir[0]->cd();
        this->fillEventHistograms(BX, BC, nVFAT[0], firedchannels[0], notfiredchannels[0], hiDiffBXandBC[0], hiRatioBXandBC[0], hiVFAT[0], hichfired[0], hichnotfired[0], hiVFATfired_perevent[0], vfatId[0]);
        if (eventIsOK){
          dir[1]->cd();
          this->fillEventHistograms(BX, BC, nVFAT[1], firedchannels[1], notfiredchannels[1], hiDiffBXandBC[1], hiRatioBXandBC[1], hiVFAT[1], hichfired[1], hichnotfired[1], hiVFATfired_perevent[1], vfatId[1]);
        } else {
          dir[2]->cd();
          this->fillEventHistograms(BX, BC, nVFAT[2], firedchannels[2], notfiredchannels[2], hiDiffBXandBC[2], hiRatioBXandBC[2], hiVFAT[2], hichfired[2], hichnotfired[2], hiVFATfired_perevent[2], vfatId[2]);
        }
      }// end of loop over events
      
      //setTitles(t_hiVFAT, "Number VFAT blocks per Event", "Number of Events");   
      //setTitles(t_hiChip, "ChipID value, max 0xfff", "Number of VFAT blocks");
      //setTitles(t_hi1010, "1010 marker, max 0xf", "Number of VFAT blocks");   
      //setTitles(t_hi1100, "1100 marker, max 0xf", "Number of VFAT blocks");   
      //setTitles(t_hi1110, "1110 marker, max 0xf", "Number of VFAT blocks");   
      //setTitles(t_hiFlag, "Flag marker value, max 0xf", "Number of VFAT blocks");   
      //setTitles(t_hiCRC, "CRC value, max 0xffff", "Number of VFAT blocks");
      //setTitles(t_hiDiffCRC, "CRC difference", "Number of VFAT blocks");
      ////setTitles(hiFake, "Fake events", "Number of Events");
      //setTitles(t_hiCh128, "Strips, max 128", "Number of VFAT blocks"); 
      //setTitles(t_hi2DCRC, "CRC VFAT", "CRC calc");  
      
      ofile->Write();
    }

    void fillVFATHistograms(VFATdata *m_vfat, TH1F* m_hiVFATsn, TH1F* m_hiCh128, TH1F* m_hiCh_notfired, TH1C* m_hiChip, TH1C* m_hi1010, TH1C* m_hi1100, TH1C* m_hi1110, TH1C* m_hiFlag, TH1C* m_hiCRC, TH1C* m_hiDiffCRC, TH2C* m_hi2DCRC, TH1F* m_hiCh128chipFired[], int & firedchannels, int & notfiredchannels)
    {
      // fill the control bits histograms
      m_hi1010->Fill(m_vfat->b1010());
      m_hi1100->Fill(m_vfat->b1100());
      m_hi1110->Fill(m_vfat->b1110());
      // fill Flag and chip id histograms
      m_hiFlag->Fill(m_vfat->Flag());
      m_hiChip->Fill(m_vfat->ChipID());
      // calculate and fill VFAT slot number
      int sn_ = m_vfat->SlotNumber();
      m_hiVFATsn->Fill(sn_);
      // calculate and fill the crc and crc_diff
      m_hiCRC->Fill(m_vfat->crc());
      // CRC check
      uint16_t b1010 = (0x000f & m_vfat->b1010());
      uint16_t b1100 = (0x000f & m_vfat->b1100());
      uint16_t b1110 = (0x000f & m_vfat->b1110());
      uint16_t flag  = (0x000f & m_vfat->Flag());
      uint16_t ec    = (0x00ff & m_vfat->EC());
      //BC             = m_vfat->BC();
      if (DEBUG) std::cout << "[gemTreeReader]: CRC read from vfat : " << std::hex << m_vfat->crc() << std::endl;
      if (DEBUG) std::cout << "[gemTreeReader]: CRC recalculated   : " << std::hex << m_vfat->crc_calc() << std::endl;
      m_hiDiffCRC->Fill(m_vfat->crc()-m_vfat->crc_calc());
      m_hi2DCRC->Fill(m_vfat->crc(), m_vfat->crc_calc());
      //I think it would be nice to time this...
      uint16_t chan0xf = 0;
      for (int chan = 0; chan < 128; ++chan) {
        if (chan < 64){
          chan0xf = ((m_vfat->lsData() >> chan) & 0x1);
          if(chan0xf) {
             m_hiCh128->Fill(chan);
             firedchannels++;
          } else {
            m_hiCh_notfired->Fill(chan);
            notfiredchannels++;
          }
        } else {
          chan0xf = ((m_vfat->msData() >> (chan-64)) & 0x1);
          if(chan0xf) {
             m_hiCh128->Fill(chan);
             firedchannels++;
          } else {
            m_hiCh_notfired->Fill(chan);
            notfiredchannels++;
          }
        }
      }// end of loop over channels 

      for(int m=0; m < NVFAT; m++){
        if(sn_ == m){
          if (DEBUG) std::cout << "[gemTreeReader]: Starting to fill hiCh128chipFired for slot : " << m << std::dec << std::endl;
          if (DEBUG) std::cout << "[gemTreeReader]: LS data           " << std::bitset<64>(m_vfat->lsData()) <<  std::endl;
          if (DEBUG) std::cout << "[gemTreeReader]: MS data           " << std::bitset<64>(m_vfat->msData()) <<  std::endl;
          uint16_t chan0xfFiredchip = 0;
          for (int chan = 0; chan < 128; ++chan) {
            if (chan < 64){
              chan0xfFiredchip = ((m_vfat->lsData() >> chan) & 0x1);
              if(chan0xfFiredchip) m_hiCh128chipFired[m]->Fill(chan);
            } else {
              chan0xfFiredchip = ((m_vfat->msData() >> (chan-64)) & 0x1);
              if(chan0xfFiredchip) m_hiCh128chipFired[m]->Fill(chan);
            }
          }
        } 
      }
    }
    void fillEventHistograms(const int& m_BX, const int& m_BC, const int & nVFAT, const int & m_firedchannels, const int& m_notfiredchannels, TH1C* m_hiDiffBXandBC, TH1C* m_hiRatioBXandBC, TH1F* m_hiVFAT, TH1C* m_hichfired, TH1C* m_hichnotfired, TH1F* m_hiVFATfired_perevent[], int vfatId[])
    {
      int diffBXandBC =  fabs(m_BX - m_BC);  
      double ratioBXandBC = (double) m_BX / m_BC;
      m_hiDiffBXandBC->Fill(diffBXandBC); 
      m_hiRatioBXandBC->Fill(ratioBXandBC);
      m_hiVFAT->Fill(nVFAT);
      for(Int_t x=0; x<NVFAT; x++) {
        m_hiVFATfired_perevent[x]->Fill(vfatId[x]);
        if (DEBUG) std::cout << std::dec << "[gemTreeReader]: Fired chip counter when filling event"  <<  vfatId[x]  << std::endl;   
      }
      m_hichfired->Fill(m_firedchannels);
      m_hichnotfired->Fill(m_notfiredchannels);
    }
};