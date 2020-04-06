#include <iostream>
#define ZBlowpu_cxx
#include "ZBlowpu.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include <set>
#include "TRandom3.h"
#include <vector>
#include "basic.h"
#include <string>
#include <map>
#include <utility>
#include "TString.h"
#include "TApplication.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace std;

void ZBlowpu::Loop()
{
    if (fChain == 0) return;
    TChain *chain=new TChain("ak4/ProcessedTree","analizpu"); // ak4'un icindeki Processed Tree deki treeleri chaine ekliyoruz
    
    chain->Add("HEG_*.root");      // -> 130GB
    // chain->Add("LEG_*.root");    // -> 180GB,263 milyon olay
    //chain->Add("ZB_*.root");      // -> 1GB
    
    ZBlowpu *t=new ZBlowpu(chain);
    
    //histogramımızı chainden cekelim
    TH1F *triggers = dynamic_cast<TH1F*>(fChain->GetCurrentFile()->Get("ak4/TriggerNames"));
    TAxis *xax = triggers->GetXaxis(); //histogramdan x axisindeki degerleri cektik
    vector <string> vectrigname; //trigger isimlerini koyacagımız bir vector olusturduk
    for (int trgidx = xax->GetFirst(); trgidx <= xax->GetLast(); trgidx++)
    {
        
        string trgName = xax->GetBinLabel(trgidx); //triggerları bin bin alıyor
        if (trgName.compare("")==0) continue; // bos binleri atla
        vectrigname.push_back(trgName);
        cout<<"Trigger "<<trgidx<<" :"<<trgName<<endl;
        
    }
    
    Long64_t nentries = t->fChain->GetEntries();
    cout<<"nentries"<<nentries<<endl;
    //nentries = 200000;
    Long64_t nbytes = 0, nb = 0;
    // TFile myFile("lowpuHEGtrigturnon.root", "RECREATE");
    TFile myFile("6Nisan_l2l3alleta_withlumi_Pt_HEG.root", "RECREATE");
    
    //CondFormat'ın icinden cektigimiz text dosyaları,degistirdiklerim var!!///
    JetCorrectorParameters *pfchs_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V1_SimpleL1_DATA_L1FastJet_AK4PFchs.txt"); //burayı degistirdim//
    JetCorrectorParameters *pfchs_l2 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V1_SimpleL1_DATA_L2Relative_AK4PFchs.txt");//burayı degistirdim//
    JetCorrectorParameters *pfchs_l3 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V1_SimpleL1_DATA_L3Absolute_AK4PFchs.txt"); //burayı degistirdim//
    JetCorrectorParameters *pfchs_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V1_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt"); //burayı degistirdim//
    
    //cout<<"Kontrol 1.Cdata txt dosyalarından sonra "<<endl;
    vector<JetCorrectorParameters> vParam_pfchs;
    vParam_pfchs.push_back(*pfchs_l1);
    vParam_pfchs.push_back(*pfchs_l2);
    vParam_pfchs.push_back(*pfchs_l3);
    vParam_pfchs.push_back(*pfchs_l2l3res);
    FactorizedJetCorrector *pfchs_jec = new FactorizedJetCorrector(vParam_pfchs);
    
    //cout<<"Kontrol 2: buraya geliyor mu ki "<<endl;
    
    bool isdijet_p;
    Float_t pf_jtpt[100],pf_jt_eta[100],pf_jt_phi[100],ptdijet_pf;
    static const int netabins = 10;
    static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,4.7};
    //static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,3.7,4.2,4.7};
    //static const double etabins[netabins+1] = {2.5,3.0,3.2,4.7};
    const double x[10][65]=
    { 
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389}, // Eta_0.0-0.5
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3637, 5220, 5492}, // Eta_0.5-1.0
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2941, 3832}, // Eta_1.0-1.5
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2500, 2640}, // Eta_1.5-2.0
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684}, // Eta_2.0-2.5
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_2.5-3.0
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_3.0-3.5
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032},// Eta_3.5-4.0
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_4.0-4.5
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032} // Eta_4.5-5.0
  
    };
    
   //const int nx[3] = {34,34,34};
    const int nx[10] = {64,63,58,54,48,40,40,40,40,40};
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    std::vector<TH1F*> vectptt[6];
    
    //butun eventlerın histogramlarını buraya tanımlıyoruz//
    TH1F *hpt[netabins];
    
    for ( int k = 0; k<netabins; ++k)
    {
        //TH1F *hlt_0   = new TH1F(Form("hlt_0_%d",k),vectrigname[0].c_str(),nx[k],&x[k][0]);
        /* TH1F *hlt_0   = new TH1F(Form("hlt_0_%d",k),vectrigname[0].c_str(),nx[k],&x[k][0]);
         vectptt[0].push_back (hlt_0); /// trigger histolarının oldugu vektor
         hlt_0->Sumw2();
         
         TH1F *hlt_1   = new TH1F(Form("hlt_1_%d",k),vectrigname[1].c_str(),nx[k],&x[k][0]);
         vectptt[1].push_back (hlt_1);
         hlt_1->Sumw2();
         
         TH1F *hlt_2   = new TH1F(Form("hlt_2_%d",k),vectrigname[2].c_str(),nx[k],&x[k][0]);
         vectptt[2].push_back (hlt_2);
         hlt_2->Sumw2();
         
         TH1F *hlt_3   = new TH1F(Form("hlt_3_%d",k),vectrigname[3].c_str(),nx[k],&x[k][0]);
         vectptt[3].push_back (hlt_3);
         hlt_3->Sumw2();*/
        
        TH1F *hlt_4   = new TH1F(Form("hlt_4_%d",k),vectrigname[4].c_str(),nx[k],&x[k][0]);
        vectptt[4].push_back (hlt_4);
        hlt_4->Sumw2();
        
        TH1F *hlt_5   = new TH1F(Form("hlt_5_%d",k),vectrigname[5].c_str(),nx[k],&x[k][0]);
        vectptt[5].push_back (hlt_5);
        hlt_5->Sumw2();
        
    }
    
    
    TLorentzVector p4pf;
    TLorentzVector p4pf2;
    std::vector<double>pfchs_ptcorr;
    //    cout << "Number of jets : " << PFJets__ << endl;
    int pf_njt = t->PFJetsCHS__;
    
    //---------------------------events dongusu icine girdik---------
    for (int i=0; i<nentries; i++)
    {
        
        cout<<"\r"<<"event number: "<<i<<"/ "<<nentries<<flush;
        chain->GetEntry(i);
        
        if (t->TriggerDecision_.empty()) continue;
        
        int pf_njt =t->PFJetsCHS__;	//corr
        pfchs_ptcorr.resize(pf_njt);	//corr
        float pf_jtpt[100],pf_jt_eta[100];	//corr
        int itrig=0;
        
        
        // cout<<"Kontrol 3: buraya geliyor mu ki "<<endl;
        
        
        for(int j=0; j<t->PFJetsCHS__; j++)
        {
            
            p4pf.SetPxPyPzE(t->PFJetsCHS__P4__fCoordinates_fX[j],t->PFJetsCHS__P4__fCoordinates_fY[j],
                            t->PFJetsCHS__P4__fCoordinates_fZ[j],t->PFJetsCHS__P4__fCoordinates_fT[j]);
            
            
            ////////////correction icin eklediklerimmm///
            pf_jtpt[j]   = p4pf.Pt();
            pf_jt_eta[j] = p4pf.Eta();
            
            float rho=t->EvtHdr__mPFRho;
            double area=t->PFJetsCHS__area_[j];
            //double ptraw = pf_jtpt[j]/PFJetsCHS__cor_[j]+rho*area;
            double ptraw = pf_jtpt[j]+rho*area;
            //cout<< "   rho: "<<rho<<"   area:  "<<area<<"   raw:  "<<ptraw<<endl;
            pfchs_jec->setJetEta(pf_jt_eta[j]);
            pfchs_jec->setJetPt(ptraw);
            pfchs_jec->setRho(rho);
            pfchs_jec->setJetA(area);
            
            pfchs_ptcorr[j]= ptraw*pfchs_jec->getCorrection();
            //cout<<"pfchs_ptcorr :"<< pfchs_ptcorr[j]<<endl;
            
            double pf_pt=pfchs_ptcorr[j];
            // double pf_pt=pf_jtpt[j];
            double pf_jteta=pf_jt_eta[j];
            //double pf_jtphi=pf_jt_phi[j];
            //cout<<"pfpt: "<< pf_pt<<"   pf_jteta: "<<pf_jteta<<endl;
            ///////correction bolgesi bitti////////
            
            
            //trigger PT aralıklarını buraya yerlestırdım
            vector<double>vecttrigcuts={49,49,64,97,133,196,300};   //yeni trigger rangeler
            //vector<double>vecttrigcuts={28,43,64,84,114,196,272};   //yeni trigger rangeler
            //vector<double>vectlumi={11702.608,70534.718,449652.826,2016829.117,6734615.464,199413819.206}; //Hannunun hesapladıkları
            vector<double>vectlumi={12410.999433624,74815.505012496,476831.591488795,2138811.439403298,7142065.860181809,211511974.809555233}; //Subat hesaplari
            // artık cutları buraya yerlestırebilirim ;)
            if( t->PFJetsCHS__tightID_[j] && (t->PFMet__et_<0.3*t->PFMet__sumEt_) && (t->PFJetsCHS__nhf_[j]<0.9) && (t->PFJetsCHS__nemf_[j]<0.9) && (t->PFJetsCHS__muf_[j]<0.9))
            {
                for ( int k = 0; k<netabins; ++k)
                {
                    
                    if (fabs(pf_jteta)>= etabins[k] && fabs(pf_jteta)< etabins[k+1])
                    {
                    
                        for(int trnameindex=4; trnameindex<6 ; trnameindex++)
                            
                        {
                           
                            for(int trdecindex=0; trdecindex<t->TriggerDecision_.size(); trdecindex++)
                            {
                               
                                if(t->TriggerDecision_[trdecindex]==trnameindex && ((pf_pt)>=vecttrigcuts[trnameindex]&& (pf_pt)<vecttrigcuts[trnameindex+1]))
                                {
                                     //vectptt[trnameindex][k]->Fill(pf_pt);
                                    double w=vectlumi[trnameindex];
                                    //vectptt[trnameindex][k]->Fill(pf_pt);
                                    vectptt[trnameindex][k]->Fill(pf_pt, 1./w);
                                }
                                
                            } //trdecindex
                            
                        } //trnameindex
                        
                    } /// fabs ifi bitiyor
                } ///etabin bitiyor
            } ///tightid li olan if cutı burda bitiyor
        } //jetler bitiyor
        
    } ///event bitiyor
    
    myFile.cd();
    myFile.Write();
    myFile.Close();
}


