
#include <cstdlib>
#include <iostream>
#include "ThinTargetMC.h"
#include "MakeReweight.h"     
#include <boost/property_tree/ptree.hpp> 
#include <boost/property_tree/xml_parser.hpp>  
#include "MakeReweight.h"

namespace NeutrinoFluxReweight{ 

  ThinTargetMC* ThinTargetMC::instance = 0;
  
  ThinTargetMC::ThinTargetMC(){
   

    const char* ppfxDir = getenv("PPFX_DIR");
    char dirData[400]; 
    sprintf(dirData,"%s/data",ppfxDir);
    spart_prod.push_back("pip");
    spart_prod.push_back("pim");
    spart_prod.push_back("kap");
    spart_prod.push_back("kam");
    spart_prod.push_back("prt");
    spart_prod.push_back("klong");
    spart_prod.push_back("kshort");
    spart_prod.push_back("neu"); const int idx_neu = 7;
    mom_inc.push_back(12);mom_inc.push_back(20);
    mom_inc.push_back(31);mom_inc.push_back(40);
    mom_inc.push_back(50);mom_inc.push_back(60);
    mom_inc.push_back(70);mom_inc.push_back(80);
    mom_inc.push_back(90);mom_inc.push_back(100);
    mom_inc.push_back(110);mom_inc.push_back(120);

    int part_size = spart_prod.size();
    int mom_size  = mom_inc.size();
    
    int expValue = MakeReweight::getInstance()->getExp();
    std::cout << "Exp Value in ThinTargetMC: " << expValue << std::endl;
   
   std::string model;
   if (expValue == 1) {
       model = "FTFP";
   } else if (expValue == 2) {
       model = "QGSP";
   } else {
       throw std::runtime_error("Invalid expValue: Expected 1 or 2");
}
   
    for(int i=0;i<part_size;i++){
      if(i!=idx_neu){
	fpC_x[i] = new TFile(Form("%s/MC/FTFP/%d/invxs_%s_%s_BERT.root",dirData,expValue,spart_prod[i].c_str(), model.c_str()),"read");
       //std::string fileName = Form("%s/MC/FTFP/%d/invxs_%s_%s_BERT.root", dirData, expValue, spart_prod[i].c_str(), model.c_str());

	//std::cout<<"The histogram being used is:"<< fileName<< std::endl;
	for(int j=0;j<mom_size;j++)vpC_x[i].push_back((TH2D*)fpC_x[i]->Get(Form("xF_pT_%dGeV",mom_inc[j])));
      }
      else if(i==idx_neu){
	fpC_x[i] = new TFile(Form("%s/MC/FTFP/%d/yield_%s_%s_BERT.root",dirData,expValue,spart_prod[i].c_str(), model.c_str()),"read");
	for(int j=0;j<mom_size;j++)vpC_n.push_back((TH1D*)fpC_x[i]->Get(Form("dndxf_%dGeV",mom_inc[j])));
      }
    }
    //QEL fraction:
    spart_qe_corr.push_back("prt");
    spart_qe_corr.push_back("neu");
    int qe_size = spart_qe_corr.size();
    for(int i=0;i<qe_size;i++){
      if(i<qe_size-1){
	fqe_corr[i] = new TFile(Form("%s/MC/FTFP/%d/invxs_qe_corr_%s.root",dirData,expValue,spart_qe_corr[i].c_str()),"read");
	for(size_t j=0;j<mom_inc.size();j++)vqe_corr_p.push_back((TH2D*)fqe_corr[i]->Get(Form("frac_prod_xF_pT_%dGeV",mom_inc[j])));
      }
      else if(i==(qe_size-1)){
	fqe_corr[i] = new TFile(Form("%s/MC/FTFP/%d/yield_qe_corr_%s.root",dirData,expValue,spart_qe_corr[i].c_str()),"read");
	for(int j=0;j<mom_size;j++)vqe_corr_n.push_back((TH1D*)fqe_corr[i]->Get(Form("frac_prod_xf_%dGeV",mom_inc[j])));
      }      
    }

    //Scaling:
    const int Nscl = 5;
    const char* cscl_parts[Nscl] = {"pip","pim","kap","kam","prt"};
    const int Nmomscl = 11;
    const int mom_scale[Nmomscl] = {12,20,31,40,50,60,70,80,100,120,158};
    for(Int_t ipart=0;ipart<Nscl;ipart++){
      ThinTargetMC::fTTscale.push_back(new TFile(Form("%s/SCALING/%s_scaling.root",dirData,cscl_parts[ipart])));
    }
    ThinTargetMC::fTTscale.push_back(new TFile(Form("%s/SCALING/%s_scaling.root",dirData,"neu")));
    //Loading scale histograms:
    for(Int_t ipart=0;ipart<Nscl;ipart++){
      std::vector<TH2F*> tmp_scl;
      for(int imom=0;imom<Nmomscl;imom++){
	tmp_scl.push_back((TH2F*)ThinTargetMC::fTTscale[ipart]->Get(Form("xF_pT_%dGeV",mom_scale[imom])));
      }
      hTTScl.push_back(tmp_scl);
      tmp_scl.clear(); 
    }
    //neutron scaling
    for(int imom=0;imom<Nmomscl;imom++){
      hTTScl_n.push_back((TH1F*)ThinTargetMC::fTTscale[Nscl]->Get(Form("xF_%dGeV",mom_scale[imom])));
    }

  }
  
  double ThinTargetMC::getMCval_pC_X(double incP, double xf,double pt, int pdgcode){

    //check:
    if(incP<12)return -1;    
    if(pdgcode!=211 && pdgcode!=-211 && pdgcode!=321 && pdgcode!=-321 && pdgcode!=2212 && pdgcode!=2112 && pdgcode!=130 && pdgcode!=310)return -1;    
    //idx:
    int idx_part = -1;
    //int idx_qe_corr = -1; // unused
    int idx_lowp = -1;
    int idx_hip  = -1;
    for(size_t i=0;i<mom_inc.size()-1;i++){
      if(incP>=double(mom_inc[i]) && incP<=double(mom_inc[i+1])){
	idx_lowp=i;
	idx_hip =i+1;
      }
    }
    if(idx_lowp==-1){
      idx_lowp = mom_inc.size()-2;
      idx_hip  = mom_inc.size()-1;
    }
    if(pdgcode ==  211)idx_part=0;
    if(pdgcode == -211)idx_part=1;
    if(pdgcode ==  321)idx_part=2;
    if(pdgcode == -321)idx_part=3;
    if(pdgcode == 2212){
      idx_part    = 4;
      //idx_qe_corr = 0; // unused
    }
    if(pdgcode == 2112){
       idx_part    = 7;
       //idx_qe_corr = 1; // unused
    }
    if(pdgcode ==  130)idx_part=5;
    if(pdgcode ==  310)idx_part=6;
    
    double mcval    = 0.0;
    // some of this code exists to set qe_corr, which is not used
    //double qe_corr = 1.0;
    if(idx_part<7){
      int binp     = vpC_x[idx_part][idx_lowp]->FindBin(xf,pt);
      double mclow = vpC_x[idx_part][idx_lowp]->GetBinContent(binp);
      double mchi  = vpC_x[idx_part][idx_hip]->GetBinContent(binp);
      mcval = mclow + (incP-double(mom_inc[idx_lowp]))*(mchi-mclow)/(double(mom_inc[idx_hip])-double(mom_inc[idx_lowp]));
      //if(idx_qe_corr==0){
	//int binqe       = vqe_corr_p[idx_lowp]->FindBin(xf,pt);
	//double qelow    = vqe_corr_p[idx_lowp]->GetBinContent(binqe);
	//double qehi     = vqe_corr_p[idx_hip]->GetBinContent(binqe);
	 //qe_corr  = qelow + (incP-double(mom_inc[idx_lowp]))*(qehi-qelow)/(double(mom_inc[idx_hip])-double(mom_inc[idx_lowp]));
      //}
    }
    else if(idx_part==7){
      int binp     = vpC_n[idx_lowp]->FindBin(xf);
      double mclow = vpC_n[idx_lowp]->GetBinContent(binp);
      double mchi  = vpC_n[idx_hip]->GetBinContent(binp);
      mcval = mclow + (incP-double(mom_inc[idx_lowp]))*(mchi-mclow)/(double(mom_inc[idx_hip])-double(mom_inc[idx_lowp]));
      //int binqe    = vqe_corr_n[idx_lowp]->FindBin(xf);
      //double qelow = vqe_corr_n[idx_lowp]->GetBinContent(binqe);
      //double qehi  = vqe_corr_n[idx_hip]->GetBinContent(binqe);
      //qe_corr  = qelow + (incP-double(mom_inc[idx_lowp]))*(qehi-qelow)/(double(mom_inc[idx_hip])-double(mom_inc[idx_lowp]));
    }
    
//    mcval /=qe_corr;
    
    //check:
    if(mcval<1.e-12 || mcval!=mcval)return -1;
    
    return mcval;
    
  }
  
  ////////////////////
double ThinTargetMC::getMCxs_pC_piK(int genid, double inc_mom){
    double xx[13] ={12,20,31,40,50,60,70,80,90,100,110,120,158};
    double yy[13] ={153386793./197812683.,160302538./197811564.,164508480./197831250.,166391359./197784915.,
		    167860919./197822312.,168882647./197807739.,169681805./197803099.,170311264./197811098.,
		    170860912./197822002.,171309291./197834756.,171651963./197811822.,171991260./197823012.,
		  172902228./197804669.};
    
    int idx_lowp = -1;
    int idx_hip  = -1;
    for(int i=0;i<12;i++){
      if(inc_mom>=double(xx[i]) && inc_mom<double(xx[i+1])){
	idx_lowp=i;
	idx_hip =i+1;
      }
    }
    double frac_low = yy[idx_lowp];
    double frac_hi  = yy[idx_hip];
    double frac_m   =  frac_low + (inc_mom-double(xx[idx_lowp]))*(frac_hi-frac_low)/(double(xx[idx_hip])-double(xx[idx_lowp]));
    
    if(genid==0)return frac_m*243.2435;
    else if(genid>0)return frac_m;
    else{
      std::cout<<"Something is wrong with gen "<<std::endl;
      return 1.0;
    }
    
  }
  /////
double ThinTargetMC::getMCxs_pC_nucleon(int genid, int pdg, double inc_mom){
    double xx[13] ={12,20,31,40,50,60,70,80,90,100,110,120,158};
    double yy[13] ={153386793./197812683.,160302538./197811564.,164508480./197831250.,166391359./197784915.,
		    167860919./197822312.,168882647./197807739.,169681805./197803099.,170311264./197811098.,
		    170860912./197822002.,171309291./197834756.,171651963./197811822.,171991260./197823012.,
		  172902228./197804669.};
    
    int idx_lowp = -1;
    int idx_hip  = -1;
    for(int i=0;i<12;i++){
      if(inc_mom>=double(xx[i]) && inc_mom<double(xx[i+1])){
	idx_lowp=i;
	idx_hip =i+1;
      }
    }
    double frac_low = yy[idx_lowp];
    double frac_hi  = yy[idx_hip];
    double frac_m   =  frac_low + (inc_mom-double(xx[idx_lowp]))*(frac_hi-frac_low)/(double(xx[idx_hip])-double(xx[idx_lowp]));
    
    if(genid==0 && pdg==2212)return frac_m*243.2435;
    if(genid>0  && pdg==2212)return frac_m;
    if(genid==0 && pdg==2112)return frac_m;
    if(genid>0  && pdg==2112)return frac_m/243.2435;
    return 1.0;

  }

  ThinTargetMC* ThinTargetMC::getInstance(){
    if (instance == 0) instance = new ThinTargetMC;
    return instance;
  }
  
}
