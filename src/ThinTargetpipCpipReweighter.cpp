
#include "ThinTargetpipCpipReweighter.h"
#include <iostream>
#include <fstream>
#include "CentralValuesAndUncertainties.h"
#include "ThinTargetpipCpipMC.h"
#include "ThinTargetpipCpipBins.h"
#include <cstdlib>

namespace NeutrinoFluxReweight{
  ThinTargetpipCpipReweighter::ThinTargetpipCpipReweighter(int iuniv, const ParameterTable& cv_pars, const ParameterTable& univ_pars):iUniv(iuniv),cvPars(cv_pars),univPars(univ_pars){
    
    ThinTargetpipCpipBins* Thinbins =  ThinTargetpipCpipBins::getInstance();
    
     vbin_data_pipC_pip.reserve(Thinbins->GetNbins_pipC_pip());
     vbin_data_pipC_pim.reserve(Thinbins->GetNbins_pipC_pim());
     vbin_data_pipC_kp.reserve(Thinbins->GetNbins_pipC_kp());
     vbin_data_pipC_km.reserve(Thinbins->GetNbins_pipC_km());
     vbin_data_pipC_p.reserve(Thinbins->GetNbins_pipC_p()); 
     vbin_data_pipC_k0s.reserve(Thinbins->GetNbins_pipC_k0s()); 
     vbin_data_pipC_lam.reserve(Thinbins->GetNbins_pipC_lam()); 
     vbin_data_pipC_alam.reserve(Thinbins->GetNbins_pipC_alam()); 
  
  
  
    data_prod_xs = univPars.getParameterValue("inel_piC_xsec_scale"); 
 
    char namepar[100];
    for(int ii=0;ii<Thinbins->GetNbins_pipC_pip();ii++){
      
      sprintf(namepar,"ThinTarget_pipC_%s_stat_%d","pip",ii);
      double data_cv  = cvPars.getParameterValue(std::string(namepar));
      double data_sta = univPars.getParameterValue(std::string(namepar));

    sprintf(namepar,"ThinTarget_pipC_%s_sys_%d","pip",ii);
    double data_sys = univPars.getParameterValue(std::string(namepar));

    vbin_data_pipC_pip.push_back(data_sys + data_sta - data_cv);
     
    }
   

  for(int ii=0;ii<Thinbins->GetNbins_pipC_pim();ii++){

      sprintf(namepar,"ThinTarget_pipC_%s_stat_%d","pim",ii);
      double data_cv  = cvPars.getParameterValue(std::string(namepar));
      double data_sta = univPars.getParameterValue(std::string(namepar));
    sprintf(namepar,"ThinTarget_pipC_%s_sys_%d","pim",ii);
    double data_sys = univPars.getParameterValue(std::string(namepar));

    vbin_data_pipC_pim.push_back(data_sys + data_sta - data_cv);
 }

         
 for(int ii=0;ii<Thinbins->GetNbins_pipC_kp();ii++){

      sprintf(namepar,"ThinTarget_pipC_%s_stat_%d","kp",ii);
      double data_cv  = cvPars.getParameterValue(std::string(namepar));
      double data_sta = univPars.getParameterValue(std::string(namepar));
    sprintf(namepar,"ThinTarget_pipC_%s_sys_%d","kp",ii);
    double data_sys = univPars.getParameterValue(std::string(namepar));

    vbin_data_pipC_kp.push_back(data_sys + data_sta - data_cv);
 }

 for(int ii=0;ii<Thinbins->GetNbins_pipC_km();ii++){

      sprintf(namepar,"ThinTarget_pipC_%s_stat_%d","km",ii);
      double data_cv  = cvPars.getParameterValue(std::string(namepar));
      double data_sta = univPars.getParameterValue(std::string(namepar));
    sprintf(namepar,"ThinTarget_pipC_%s_sys_%d","km",ii);
    double data_sys = univPars.getParameterValue(std::string(namepar));

    vbin_data_pipC_km.push_back(data_sys + data_sta - data_cv);
 }

 for(int ii=0;ii<Thinbins->GetNbins_pipC_p();ii++){

      sprintf(namepar,"ThinTarget_pipC_%s_stat_%d","p",ii);
      double data_cv  = cvPars.getParameterValue(std::string(namepar));
      double data_sta = univPars.getParameterValue(std::string(namepar));
    sprintf(namepar,"ThinTarget_pipC_%s_sys_%d","p",ii);
    double data_sys = univPars.getParameterValue(std::string(namepar));

    vbin_data_pipC_p.push_back(data_sys + data_sta - data_cv);
 }

 for(int ii=0;ii<Thinbins->GetNbins_pipC_k0s();ii++){

      sprintf(namepar,"ThinTarget_pipC_%s_stat_%d","k0s",ii);
      double data_cv  = cvPars.getParameterValue(std::string(namepar));
      double data_sta = univPars.getParameterValue(std::string(namepar));
    sprintf(namepar,"ThinTarget_pipC_%s_sys_%d","k0s",ii);
    double data_sys = univPars.getParameterValue(std::string(namepar));

    vbin_data_pipC_k0s.push_back(data_sys + data_sta - data_cv);
 }
 for(int ii=0;ii<Thinbins->GetNbins_pipC_lam();ii++){

      sprintf(namepar,"ThinTarget_pipC_%s_stat_%d","lam",ii);
      double data_cv  = cvPars.getParameterValue(std::string(namepar));
      double data_sta = univPars.getParameterValue(std::string(namepar));
    sprintf(namepar,"ThinTarget_pipC_%s_sys_%d","lam",ii);
    double data_sys = univPars.getParameterValue(std::string(namepar));

    vbin_data_pipC_lam.push_back(data_sys + data_sta - data_cv);
 }
 for(int ii=0;ii<Thinbins->GetNbins_pipC_alam();ii++){

      sprintf(namepar,"ThinTarget_pipC_%s_stat_%d","alam",ii);
      double data_cv  = cvPars.getParameterValue(std::string(namepar));
      double data_sta = univPars.getParameterValue(std::string(namepar));
    sprintf(namepar,"ThinTarget_pipC_%s_sys_%d","alam",ii);
    double data_sys = univPars.getParameterValue(std::string(namepar));

    vbin_data_pipC_alam.push_back(data_sys + data_sta - data_cv);
 }

}
       
  ThinTargetpipCpipReweighter::~ThinTargetpipCpipReweighter(){
    
  }
  bool ThinTargetpipCpipReweighter::canReweight(const InteractionData& aa){
    std::string mode(getenv("MODE"));
   if(aa.Proc.find("Inelastic")>100){       
     
    return false; 
     } 
    
   
     bool corr_inc = (aa.Inc_pdg ==211 && aa.Inc_P >= 12);  
     bool corr_pro = (aa.Prod_pdg == 211 ||aa.Prod_pdg == -211 || aa.Prod_pdg == 321|| aa.Prod_pdg == -321 || aa.Prod_pdg == 310 || aa.Prod_pdg == 2212 || aa.Prod_pdg == 3122 || aa.Prod_pdg == -3122 );  
    if(corr_inc && corr_pro)return true;
    else return false;
   bool is_wrong_volume = aa.Vol != "TGT1" && aa.Vol != "BudalMonitor" && aa.Vol != "Budal_HFVS" && aa.Vol != "Budal_VFHS";
    if( (mode=="REF") || (mode=="OPT") ){
      is_wrong_volume = aa.Vol != "TargetFinHorizontal" && aa.Vol != "TargetNoSplitSegment" && aa.Vol!="tCoreLog";
    }
    if(is_wrong_volume)return false;
  bool can_reweight = false;
 //int bin = -1; // this "bin" is not used
 ThinTargetpipCpipBins*  Thinbins =  ThinTargetpipCpipBins::getInstance();
if(aa.Prod_pdg == 211) {
      int bin = Thinbins->pipC_pip_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
      if(bin >=0 && bin< 200) can_reweight = true;
    }
if(aa.Prod_pdg == -211) {
      int bin = Thinbins->pipC_pim_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
      if(bin >=0 && bin<200) can_reweight = true;
    }
if(aa.Prod_pdg == 321) {
      int bin = Thinbins->pipC_kp_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
      if(bin >=0 && bin <35)can_reweight = true;
    }
if(aa.Prod_pdg == -321) {
      int bin = Thinbins->pipC_km_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
      if(bin >=0 && bin< 35) can_reweight = true;
    }
if(aa.Prod_pdg == 2212) {
      int bin = Thinbins->pipC_p_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
      if(bin >=0 && bin< 72) can_reweight = true;
    }
if(aa.Prod_pdg == 310) {
      int bin = Thinbins->pipC_k0s_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
      if(bin >=0 && bin< 35) can_reweight = true;
    }
if(aa.Prod_pdg == 3122) {
      int bin = Thinbins->pipC_lam_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
      if(bin >=0 && bin<21) can_reweight = true;
    }
if(aa.Prod_pdg == -3122) {
      int bin = Thinbins->pipC_alam_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
      if(bin >=0 && bin< 11) can_reweight = true;
    }

return can_reweight;
  }
  
  double ThinTargetpipCpipReweighter::calculateWeight(const InteractionData& aa){

	  
    double wgt = 1.0;
    double low_value = 1.e-18;     
    int bin = -1;
    //bool right_inc = aa.Inc_pdg == 211; // right_inc not used
   if(aa.Inc_pdg != 211){std::cout<<"Can reweight is not working properly, setting wgt 1"<<std::endl;
   return wgt; } 

      ThinTargetpipCpipBins*  Thinbins =  ThinTargetpipCpipBins::getInstance();
        ThinTargetpipCpipMC*  mc =  ThinTargetpipCpipMC::getInstance();
    
      if (aa.Prod_pdg == 211) {
        bin = Thinbins->pipC_pip_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
    } else if (aa.Prod_pdg == -211) {
        bin = Thinbins->pipC_pim_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
    } else if (aa.Prod_pdg == 321) {
        bin = Thinbins->pipC_kp_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
    } else if (aa.Prod_pdg == -321) {
        bin = Thinbins->pipC_km_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
    } else if (aa.Prod_pdg == 2212) {
        bin = Thinbins->pipC_p_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
    } else if (aa.Prod_pdg == 310) {
        bin = Thinbins->pipC_k0s_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
    } else if (aa.Prod_pdg == 3122) {
        bin = Thinbins->pipC_lam_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
    } else if (aa.Prod_pdg == -3122) {
        bin = Thinbins->pipC_alam_BinID(aa.Prod_P,aa.Theta,aa.Inc_pdg,aa.Prod_pdg);
    } else {
        std::cerr << "Invalid Prod_pdg: " << aa.Prod_pdg << std::endl;
        return wgt;
    }
    
    if (bin < 0) {
        std::cerr << "Invalid bin: " << bin <<"Momentum,Theta of Produced particle are"<<" "<< aa.Prod_P<<" "<< aa.Theta<<std::endl;
        return wgt;
    } 
    
    double data_wgt = 1.0;
    if (aa.Prod_pdg == 211) data_wgt = vbin_data_pipC_pip[bin];
    else if (aa.Prod_pdg == -211) data_wgt = vbin_data_pipC_pim[bin];
    else if (aa.Prod_pdg == 321) data_wgt = vbin_data_pipC_kp[bin];
    else if (aa.Prod_pdg == -321) data_wgt = vbin_data_pipC_km[bin];
    else if (aa.Prod_pdg == 2212) data_wgt = vbin_data_pipC_p[bin];
    else if (aa.Prod_pdg == 310) data_wgt = vbin_data_pipC_k0s[bin];
    else if (aa.Prod_pdg == 3122) data_wgt =vbin_data_pipC_lam[bin];
    else if (aa.Prod_pdg == -3122) data_wgt = vbin_data_pipC_alam[bin];
    
    double mc_wgt = 1.0;
  if (aa.Prod_pdg == 211) {
        mc_wgt = mc->getMCval(aa.Prod_P, aa.Theta, 211, 211);
    } else if (aa.Prod_pdg == -211) {
        mc_wgt = mc->getMCval(aa.Prod_P, aa.Theta, 211, -211);
    } else if (aa.Prod_pdg == 321) {
        mc_wgt = mc->getMCval(aa.Prod_P, aa.Theta, 211, 321);
    } else if (aa.Prod_pdg == -321) {
        mc_wgt = mc->getMCval(aa.Prod_P, aa.Theta, 211, -321);
    } else if (aa.Prod_pdg == 2212) {
        mc_wgt = mc->getMCval(aa.Prod_P, aa.Theta, 211, 2212);
    } else if (aa.Prod_pdg == 310) {
        mc_wgt = mc->getMCval(aa.Prod_P, aa.Theta, 211, 310);
    } else if (aa.Prod_pdg == 3122) {
        mc_wgt = mc->getMCval(aa.Prod_P, aa.Theta, 211, 3122);
    } else if (aa.Prod_pdg == -3122) {
        mc_wgt = mc->getMCval(aa.Prod_P, aa.Theta, 211, -3122);
    }
        
   if (mc_wgt < low_value){
           std::cerr << "Low MC Value "<< std::endl;
   }
   
       wgt = data_wgt / mc_wgt;



  if (wgt < 0) {
        return 1.0;
    }
    if (wgt > 10) {
        std::cerr << "BIG WGT IN TTMESONINC " << iUniv << " " << wgt << " " << aa.Inc_P << " " << aa.xF << " " << aa.Pt << " " << aa.Prod_pdg << std::endl;
        return 1.0;
        }


//std::cout<<"Weight"<<wgt<<std::endl;

    return wgt;
    }
    
    }
    
  


