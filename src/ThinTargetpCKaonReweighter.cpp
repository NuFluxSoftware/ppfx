
#include "ThinTargetpCKaonReweighter.h"
#include <iostream>
#include "MakeReweight.h"
#include "ThinTargetBins.h"
#include "ThinTargetMC.h"

namespace NeutrinoFluxReweight{
  
  ThinTargetpCKaonReweighter::ThinTargetpCKaonReweighter(int iuniv, const ParameterTable& cv_pars, const ParameterTable& univ_pars):iUniv(iuniv),cvPars(cv_pars),univPars(univ_pars){
  
    ThinTargetBins* Thinbins =  ThinTargetBins::getInstance();
    vbin_data_kap.reserve(Thinbins->GetNbins_pC_KX_NA49());
    vbin_data_kam.reserve(Thinbins->GetNbins_pC_KX_NA49());
    mipp_vbin_data_kap_pip.reserve(Thinbins->GetNbins_pC_KX_MIPP());
    mipp_vbin_data_kam_pim.reserve(Thinbins->GetNbins_pC_KX_MIPP());

    // const boost::interprocess::flat_map<std::string, double>& cv_table   = cvPars.getMap();
    // const boost::interprocess::flat_map<std::string, double>& univ_table = univPars.getMap();
       
    data_prod_xs = univPars.getParameterValue("prod_prtC_xsec");
    char namepar[100];
    for(int ii=0;ii<Thinbins->GetNbins_pC_KX_NA49();ii++){
      
      sprintf(namepar,"ThinTargetLowxF_pC_%s_sys_%d","kap",ii);
      double data_cv  = cvPars.getParameterValue(std::string(namepar));
      double data_sys = univPars.getParameterValue(std::string(namepar));
      sprintf(namepar,"ThinTargetLowxF_pC_%s_stats_%d","kap",ii);
      double data_sta = univPars.getParameterValue(std::string(namepar));
      vbin_data_kap.push_back(data_sta + data_sys - data_cv);
      
      sprintf(namepar,"ThinTargetLowxF_pC_%s_sys_%d","kam",ii);
      data_cv  = cvPars.getParameterValue(std::string(namepar));
      data_sys = univPars.getParameterValue(std::string(namepar));
      sprintf(namepar,"ThinTargetLowxF_pC_%s_stats_%d","kam",ii);
      data_sta = univPars.getParameterValue(std::string(namepar));
      vbin_data_kam.push_back(data_sta + data_sys - data_cv);
      
    }
    
    for(int ii=0;ii<Thinbins->GetNbins_pC_KX_MIPP();ii++){
      
      sprintf(namepar,"ThinTarget_kap_pip_sys_%d",ii);
      double data_cv  = cvPars.getParameterValue(std::string(namepar));
      double data_sys = univPars.getParameterValue(std::string(namepar));
      sprintf(namepar,"ThinTarget_kap_pip_stats_%d",ii);
      double data_sta = univPars.getParameterValue(std::string(namepar));
      mipp_vbin_data_kap_pip.push_back(data_sta + data_sys - data_cv);

      sprintf(namepar,"ThinTarget_kam_pim_sys_%d",ii);
      data_cv  = cvPars.getParameterValue(std::string(namepar));
      data_sys = univPars.getParameterValue(std::string(namepar));
      sprintf(namepar,"ThinTarget_kam_pim_stats_%d",ii);
      data_sta = univPars.getParameterValue(std::string(namepar));
      mipp_vbin_data_kam_pim.push_back(data_sta + data_sys - data_cv);

    } 
    
  }
  
  ThinTargetpCKaonReweighter::~ThinTargetpCKaonReweighter(){
    
  }
  bool ThinTargetpCKaonReweighter::canReweight(const InteractionData& aa){
    //checking:
    std::string mode(getenv("MODE"));
    if(aa.Inc_pdg != 2212)return false;
    if(aa.Inc_P < 12.0)return false;
    //volume check: 
    bool is_wrong_volume = aa.Vol != "TGT1" && aa.Vol != "BudalMonitor" && aa.Vol != "Budal_HFVS" && aa.Vol != "Budal_VFHS";
    if( (mode=="REF") || (mode=="OPT") ){
      is_wrong_volume = aa.Vol != "TargetFinHorizontal" && aa.Vol != "TargetNoSplitSegment" && aa.Vol!="tCoreLog";
    }
    if(is_wrong_volume)return false;
    //
    if(aa.Prod_pdg != 321 && aa.Prod_pdg != -321 && aa.Prod_pdg != 310 && aa.Prod_pdg != 130)return false;
    
    //Looking for low pz kaon:
    ThinTargetBins*  Thinbins =  ThinTargetBins::getInstance();
    int bin      = Thinbins->BinID_pC_k(aa.xF,aa.Pt,aa.Prod_pdg);
    if(bin>=0)return true; //found NA49 kaon 
    
    //Looking for high pz kaon:
    int mipp_bin = Thinbins->mipp_BinID_pC_k(aa.Pz,aa.Pt,aa.Prod_pdg);
    if(mipp_bin<0)return false; //no mipp thin target kaon  

    double inc_mom[3]  = {aa.Inc_P4[0], aa.Inc_P4[1], aa.Inc_P4[2]};
    double prod_mom[3] = {aa.Prod_P4[0],aa.Prod_P4[1],aa.Prod_P4[2]};
    double vtx_int[3]  = {aa.Vtx[0],aa.Vtx[1],aa.Vtx[2]};
    
    MakeReweight*  makerew =  MakeReweight::getInstance();
    if(iUniv==-1)tt_pCPionRew = (makerew->cv_rw)->THINTARGET_PC_PION_Universe;
    else tt_pCPionRew = (makerew->vec_rws[iUniv])->THINTARGET_PC_PION_Universe;   

    InteractionData aux_intData(aa.gen, inc_mom,2212,prod_mom,211,aa.Vol,aa.Proc,vtx_int);   
    bool there_is_pip = tt_pCPionRew->canReweight(aux_intData);    
    
    InteractionData aux_intData2(aa.gen, inc_mom,2212,prod_mom,-211,aa.Vol,aa.Proc,vtx_int);   
    bool there_is_pim = tt_pCPionRew->canReweight(aux_intData2); 

    if(aa.Prod_pdg == 321)return there_is_pip;
    else if(aa.Prod_pdg ==-321)return there_is_pim;
    else if(aa.Prod_pdg == 130 || aa.Prod_pdg == 310)return there_is_pim && there_is_pip;
    else{
      return false;
    }
  }
  
  double ThinTargetpCKaonReweighter::calculateWeight(const InteractionData& aa){
    
    double wgt = 1.0;

    ThinTargetBins*  Thinbins =  ThinTargetBins::getInstance();
    int bin      = Thinbins->BinID_pC_k(aa.xF,aa.Pt,aa.Prod_pdg);
    int mipp_bin = Thinbins->mipp_BinID_pC_k(aa.Pz,aa.Pt,aa.Prod_pdg);
    if(bin<0 && mipp_bin<0)return wgt;  //double check
    
    double dataval = -1;
    
    if(bin>=0){
      if(aa.Prod_pdg == 321)dataval = vbin_data_kap[bin];
      else if(aa.Prod_pdg == -321) dataval = vbin_data_kam[bin];
      else if(aa.Prod_pdg == 130 || aa.Prod_pdg == 310){
        dataval = 0.25*(vbin_data_kap[bin]+3.*vbin_data_kam[bin]);
      }
      else{
        std::cout<<"Something is wrong, pdg_prod: "<< aa.Prod_pdg  <<std::endl;
        return wgt;
      }
    }
    else if(mipp_bin>=0){
   
      MakeReweight*  makerew =  MakeReweight::getInstance();
      if(iUniv==-1)tt_pCPionRew = (makerew->cv_rw)->THINTARGET_PC_PION_Universe;
      else tt_pCPionRew = (makerew->vec_rws[iUniv])->THINTARGET_PC_PION_Universe; 
      
      double inc_mom[3]  = {aa.Inc_P4[0], aa.Inc_P4[1], aa.Inc_P4[2]};
      double prod_mom[3] = {aa.Prod_P4[0],aa.Prod_P4[1],aa.Prod_P4[2]};
      double vtx_int[3]  = {aa.Vtx[0],aa.Vtx[1],aa.Vtx[2]};
      double pip_data = -1.0; 
      double pim_data = -1.0; 
      int bin_pi = -1;
      InteractionData aux_aa(aa.gen, inc_mom,2212,prod_mom,211,aa.Vol,aa.Proc,vtx_int);    
      if(aux_aa.xF<=0.5){
        bin_pi = Thinbins->BinID_pC_pi(aux_aa.xF,aux_aa.Pt,aux_aa.Prod_pdg);
        if(bin_pi<0)return 1.0;
        pip_data = tt_pCPionRew->vbin_data_pip[bin_pi];
        pim_data = tt_pCPionRew->vbin_data_pim[bin_pi];	
      }
      else if(aux_aa.xF>0.5){
        bin_pi = Thinbins->barton_BinID_pC_pi(aux_aa.xF,aux_aa.Pt,aux_aa.Prod_pdg);
        if(bin_pi<0)return 1.0;
        pip_data = tt_pCPionRew->bart_vbin_data_pip[bin_pi];
        pim_data = tt_pCPionRew->bart_vbin_data_pim[bin_pi];	
      }

      if(aa.Prod_pdg == 321){
        dataval = mipp_vbin_data_kap_pip[mipp_bin]*pip_data;
      }
      else if(aa.Prod_pdg ==-321){
        dataval = mipp_vbin_data_kam_pim[mipp_bin]*pim_data;
      }
      else if(aa.Prod_pdg == 130 || aa.Prod_pdg == 310){
        dataval = 0.25*(mipp_vbin_data_kap_pip[mipp_bin]*pip_data +3.0*(mipp_vbin_data_kam_pim[mipp_bin]*pim_data));
      }
    }
         
    //checking if this is the first interaction:
    if(aa.gen==0)dataval /= data_prod_xs;
    else if(aa.gen>0)dataval /= 1.0;
    else{
      std::cout<<"Something is wrong with gen "<<std::endl;
      return wgt;
    }

    double data_scale = -1.0;
    if(aa.Prod_pdg == 321 || aa.Prod_pdg == -321){
      data_scale = calculateDataScale(aa.Inc_pdg,aa.Inc_P,aa.Prod_pdg,aa.xF,aa.Pt);   
    }
    else if(aa.Prod_pdg == 130 || aa.Prod_pdg == 310){
      data_scale = 0.25*(calculateDataScale(aa.Inc_pdg,aa.Inc_P,321,aa.xF,aa.Pt)+3.*calculateDataScale(aa.Inc_pdg,aa.Inc_P,-321,aa.xF,aa.Pt));  
      
    }
    dataval *= data_scale;

    ThinTargetMC*  mc =  ThinTargetMC::getInstance();
    double mc_cv   = mc->getMCval_pC_X(aa.Inc_P,aa.xF,aa.Pt,aa.Prod_pdg);
    double mc_prod = mc->getMCxs_pC_piK(aa.gen,aa.Inc_P);

    mc_cv /= mc_prod;

    if(mc_cv<1.e-12)return wgt;
    wgt = dataval/mc_cv;

    if(wgt<0){
      //std::cout<<"TTPCK check wgt(<0) "<<iUniv<<" "<<aa.Inc_P<<" "<<aa.xF<<" "<<aa.Pt<<" "<<aa.Prod_pdg<<std::endl;
      return 1.0;
    }
     if(wgt>10){
      //std::cout<<"TTPCK check wgt(<0) "<<iUniv<<" "<<aa.Inc_P<<" "<<aa.xF<<" "<<aa.Pt<<" "<<aa.Prod_pdg<<std::endl;
      return 1.0;
    }
    return wgt;
    
  }

  double ThinTargetpCKaonReweighter::calculateDataScale(int inc_pdg, double inc_mom, int prod_pdg,double xf, double pt){
    double scaling_violation = 1.0;
    ThinTargetMC*  dtH =  ThinTargetMC::getInstance();
    //temporary:
    const int Nscl = 11;
    const int moms[Nscl] = {12,20,31,40,50,60,70,80,100,120,158};
    
    int idx_part = -1;
    if(prod_pdg == 321)idx_part = 2;
    if(prod_pdg ==-321)idx_part = 3;
    if(idx_part<0){
      std::cout<<"Error in the prod particle"<<std::endl;
      return 1.0;
    }
    
    int binid = dtH->hTTScl[idx_part][Nscl-1]->FindBin(xf,pt);
    double scl_ref158 = dtH->hTTScl[idx_part][Nscl-1]->GetBinContent(binid);    
    
    int idx_lowp = -1;
    int idx_hip  = -1;
    for(int i=0;i<Nscl-1;i++){
      if(inc_mom>=double(moms[i]) && inc_mom<double(moms[i+1])){
        idx_lowp=i;
        idx_hip =i+1;
      }
    }
    if(idx_lowp<0 || idx_hip<0){
      std::cout<<"Error calculating the scaling"<<std::endl;
      return 1.0;
    }
    double scl_low = dtH->hTTScl[idx_part][idx_lowp]->GetBinContent(binid);
    double scl_hi  = dtH->hTTScl[idx_part][idx_hip]->GetBinContent(binid);
    double scl_m   =  scl_low + (inc_mom-double(moms[idx_lowp]))*(scl_hi-scl_low)/(double(moms[idx_hip])-double(moms[idx_lowp]));
    if(scl_ref158<1.e-10){
      // std::cout<<"ref158 zero!!! "<<scl_ref158<<std::endl;
      return 1.0;
    }
    scaling_violation = scl_m/scl_ref158;
    return scaling_violation;
  }

}
