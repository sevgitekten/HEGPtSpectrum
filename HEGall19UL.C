#include <iostream>
#define HEGall19UL_cxx
#include "HEGall19UL.h"
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

//#include <time.h>

using namespace std;

void HEGall19UL::Loop()
{
    //time_t start,end; double dif;  time (&start);
    
    if (fChain == 0) return; //MakeClass
    Long64_t nentries = fChain->GetEntriesFast();//MakeClass
    Long64_t nbytes = 0, nb = 0; //MakeClass
    
    cout<<"nentries"<<nentries<<endl;
    //nentries = 1000000; //72372097 olay var!
    //TH2D *hotzonemap = (TH2D*)hotzone->Get("h2hotfilter");
    TH2D *hotzonemap = (TH2D*)hotzone->Get("h2hot_ul17_plus_hep17");
    TFile myFile("8aralik_rowPt_HEG_turonstudy.root", "RECREATE");
    
    
    //CondFormat'ın icinden cektigimiz text dosyaları,degistirdiklerim var!!///
    JetCorrectorParameters *pfchs_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V4_DATA_L1FastJet_AK4PFchs.txt");
    JetCorrectorParameters *pfchs_l2 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V4_DATA_L2Relative_AK4PFchs.txt");
    JetCorrectorParameters *pfchs_l3 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V4_DATA_L3Absolute_AK4PFchs.txt");
    JetCorrectorParameters *pfchs_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Summer19UL17_RunF_V2M5_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt");
    
    vector<JetCorrectorParameters> vParam_pfchs;
    vParam_pfchs.push_back(*pfchs_l1);
    vParam_pfchs.push_back(*pfchs_l2);
    vParam_pfchs.push_back(*pfchs_l3);
    vParam_pfchs.push_back(*pfchs_l2l3res);
    FactorizedJetCorrector *pfchs_jec = new FactorizedJetCorrector(vParam_pfchs);
    
    static const int netabins = 3;
    //static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,4.7};
    //static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
    static const double etabins[netabins+1] = {0,1.3,3.2,4.7};
    const double x[3][65]=
    {
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389}, // Eta_0.0-0.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3637, 5220, 5492}, // Eta_0.5-1.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2941, 3832}, // Eta_1.0-1.5
        /*{10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2500, 2640}, // Eta_1.5-2.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684}, // Eta_2.0-2.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_2.5-3.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_3.0-3.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032},// Eta_3.5-4.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_4.0-4.5
         {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032} // Eta_4.5-5.0
         */
    };
    
    const int nx[3] = {64,63,58};
    //const int nx[8] = {64,63,58,54,48,40,40,40};
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    std::vector<TH1F*> vectptt[7];
    
    for ( int k = 0; k<netabins; ++k)
    {
        
        TH1F *hlt_4   = new TH1F(Form("hlt_4_%d",k),"cross section",nx[k],&x[k][0]);
        vectptt[4].push_back (hlt_4);
        hlt_4->Sumw2();
        
        TH1F *hlt_5   = new TH1F(Form("hlt_5_%d",k),"cross section",nx[k],&x[k][0]);
        vectptt[5].push_back (hlt_5);
        hlt_5->Sumw2();
        
        TH1F *hlt_6   = new TH1F(Form("hlt_6_%d",k),"cross section",nx[k],&x[k][0]);
        vectptt[6].push_back (hlt_6);
        hlt_6->Sumw2();
    }
    
    TLorentzVector p4pf;
    
    //trigger PT aralıklarını buraya yerlestırdım
    vector<double>vecttrigcuts={56,84,84,97,133,272,6389};   //22 Kasim yeni trigger rangeler
    //vector<double>vecttrigcuts={49,49,64,97,133,196,6389};   //yeni trigger rangeler
    
    //========luminosity valuess
    vector<double>vectlumi={12410.999433624,74815.505012496,476831.591488795,2138811.439403298,7142065.860181809,211511974.809555233}; //Subat hesaplari
    
    //=================Define the variables for Correction========================
    
    float pf_jtpt[100],pf_jt_eta[100],pf_jt_phi[100];    //corr definition (disarda kalabilir)
    vector<float>pfchs_ptcorr;
    
    //---------------------------events dongusu icine girdik---------
    
    for (Long64_t jentry=0; jentry<nentries;jentry++)   //MakeClass
    {
        Long64_t ientry = LoadTree(jentry);             //MakeClass
        if (ientry < 0) break;                          //MakeClass
        nb = fChain->GetEntry(jentry);   nbytes += nb;  //MakeClass
        
        cout<<"\r"<<"event number: "<<jentry<<"/ "<<nentries<<flush;
        
        if (TriggerDecision_.empty()) continue;
        //if (FilterDecision_.size()!=0)  continue;
        
        int pf_njt =PFJetsCHS__;        //==========corr definition
        pfchs_ptcorr.resize(pf_njt);    //=======corr definition
        
        for(int j=0; j<PFJetsCHS__; j++)
        {
            
            p4pf.SetPxPyPzE(PFJetsCHS__P4__fCoordinates_fX[j],PFJetsCHS__P4__fCoordinates_fY[j],
                            PFJetsCHS__P4__fCoordinates_fZ[j],PFJetsCHS__P4__fCoordinates_fT[j]);
            
            
            ////////////correction icin eklediklerimmm///
            pf_jtpt[j]   = p4pf.Pt();
            pf_jt_eta[j] = p4pf.Eta();
            pf_jt_phi[j] = p4pf.Phi();
            
            float rho=EvtHdr__mPFRho;
            double area=PFJetsCHS__area_[j];
            double ptraw = pf_jtpt[j]/PFJetsCHS__cor_[j];
            //double ptraw = pf_jtpt[j]+rho*area; //undo yapmadan ptraw hesabi
            pfchs_jec->setJetEta(pf_jt_eta[j]);
            pfchs_jec->setJetPt(ptraw);
            pfchs_jec->setRho(rho);
            pfchs_jec->setJetA(area);
            
            pfchs_ptcorr[j]= ptraw*pfchs_jec->getCorrection();
            double pf_pt=pfchs_ptcorr[j];
            double pf_jteta=pf_jt_eta[j];
            double pf_phi=pf_jt_phi[j];
            ///////correction bolgesi bitti////////
            
            
            if( PFJetsCHS__tightID_[j] && (PFMet__et_<0.3*PFMet__sumEt_) && (PFJetsCHS__cemf_[j]<0.9) && (PFJetsCHS__muf_[j]<0.9))
            {
                for ( int k = 0; k<netabins; ++k)
                {
                    
                    if (fabs(pf_jteta)>= etabins[k] && fabs(pf_jteta)< etabins[k+1])
                    {
                            bool trg4=false;
                        for(int trnameindex=4; trnameindex<6 ; trnameindex++)
    
                        {
                           
                            //cout<<endl<<"1 Size oncesi : "<<endl;
                            for(int trdecindex=0; trdecindex<TriggerDecision_.size(); trdecindex++)
                            {
                                 //cout<<endl<<"2 Size sonrasi : "<<endl;
                                if(TriggerDecision_[trdecindex]==trnameindex)// && ((pf_pt)>=vecttrigcuts[trnameindex]&& (pf_pt)<vecttrigcuts[trnameindex+1]))
                                {
                                    //cout<<"3 decision sonrasi : "<<TriggerDecision_[trdecindex]<<endl;
                                    double w=vectlumi[trnameindex];
                                    //vectptt[trnameindex][k]->Fill(pf_pt);
                                    int hotzonebin=hotzonemap->FindBin(pf_jteta,pf_phi);
                                    if(hotzonemap->GetBinContent(hotzonebin)==0.0)
                                    {
                                        if (TriggerDecision_[trdecindex]==4) trg4=true;
                                        if (trg4 && TriggerDecision_[trdecindex]==5) {
                                            //cout<<endl<<"4 ve 5 i gecenler : "<<endl;
                                            vectptt[6][k]->Fill(pf_pt);
                                            //cout<<"Pt45 ler : "<<pf_pt<<endl;
                                        }
                                        //cout<<" out of hotzone map"<<endl;
                                    vectptt[trnameindex][k]->Fill(pf_pt);
                                     //cout<<"Pt ler : "<<pf_pt<<endl;
                                    } //hotzone
                                }
                                
                            } //trdecindex
                        } //trnameindex
                    } // fabs
                } //etabin
            } //tightid
        } //jets
    } //events
    
    myFile.cd();
    myFile.Write();
    myFile.Close();
    hotzone->Close();
    
    //time (&end); dif = difftime (end,start); cout<< endl<< "zaman:"<< dif<<endl;
}


