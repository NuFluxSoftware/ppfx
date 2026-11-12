
#include "TargetAttenuationReweighter.h"
#include <iostream>
#include <locale>
#include <cstdlib>
#include <stdexcept>
#include "AttenuationMC.h"
#include "MIPPNumiYieldsBins.h"
#include "MakeReweight.h"
#include "ReweightDriver.h"


// FIXME: hzpos may be uninitialized in TargetAttenuationReweighter.cpp
#if defined __clang__
#elif defined __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
namespace NeutrinoFluxReweight{
  
  TargetAttenuationReweighter::TargetAttenuationReweighter(int iuniv, const ParameterTable& cv_pars, const ParameterTable& univ_pars)
    :cvPars(cv_pars),univPars(univ_pars),iUniv(iuniv){
    
    // const boost::interprocess::flat_map<std::string, double>& this_table = univPars.getMap();
     prod_prtC_xsec = univPars.getParameterValue("prod_prtC_xsec");
     qe_prtC_xsec   = univPars.getParameterValue("qe_prtC_xsec");
     delta_sigma_piC_xsec = univPars.getParameterValue("inel_piC_xsec");
     delta_sigma_kapC_xsec_lowP  = univPars.getParameterValue("inel_kapC_xsec_lowP");
     delta_sigma_kapC_xsec_highP = univPars.getParameterValue("inel_kapC_xsec_highP");
     delta_sigma_kamC_xsec_lowP  = univPars.getParameterValue("inel_kamC_xsec_lowP");
     delta_sigma_kamC_xsec_highP = univPars.getParameterValue("inel_kamC_xsec_highP");

  }
  TargetAttenuationReweighter::~TargetAttenuationReweighter(){
    
  }
  std::vector<bool> TargetAttenuationReweighter::canReweight(const InteractionChainData& aa){
    std::vector<bool> can_rws;
    std::vector<InteractionData> vec_inter = aa.interaction_chain;
    int ninter = vec_inter.size();
    std::string mode(getenv("MODE"));
    //Looking if there is a proton-Target interaction:
    for(int ii=0;ii<ninter;ii++){          
        if(vec_inter[ii].Inc_P==0.0){
	  can_rws.push_back(false);
	  return can_rws;
	}
    if(ii==0){
      //first interaction in the target or a primary proton passing through the target 

	bool is_tgt_int = vec_inter[0].Vol == "BudalMonitor" || vec_inter[0].Vol == "TGT1" || vec_inter[0].Vol == "Budal_HFVS"  || vec_inter[0].Vol == "Budal_VFHS";
	if((mode=="REF")||(mode=="OPT")){
	  is_tgt_int = vec_inter[0].Vol == "TargetFinHorizontal" || vec_inter[0].Vol == "TargetNoSplitSegment" || vec_inter[0].Vol=="tCoreLog";
	}
	if(is_tgt_int)can_rws.push_back(true);
	else if(vec_inter[0].Inc_pdg == 2212)can_rws.push_back(true);
	else can_rws.push_back(false);
      }
      //Absorption in the target of the secondaries:
      else{
        bool starts_tgt = vec_inter[ii-1].Vol == "BudalMonitor" || vec_inter[ii-1].Vol == "TGT1" || vec_inter[ii-1].Vol == "Budal_HFVS"  || vec_inter[ii-1].Vol == "Budal_VFHS";
	bool ends_tgt   = vec_inter[ii].Vol   == "BudalMonitor" || vec_inter[ii].Vol   == "TGT1" || vec_inter[ii-1].Vol == "Budal_HFVS"  || vec_inter[ii-1].Vol == "Budal_VFHS";
	if((mode=="REF")||(mode=="OPT")){
	  starts_tgt = vec_inter[ii-1].Vol == "TargetFinHorizontal" || vec_inter[ii-1].Vol == "TargetNoSplitSegment"|| vec_inter[0].Vol=="tCoreLog";
	  ends_tgt = vec_inter[ii-1].Vol == "TargetFinHorizontal" || vec_inter[ii-1].Vol == "TargetNoSplitSegment"||vec_inter[0].Vol=="tCoreLog";
	}
	if(starts_tgt && ends_tgt){
	  can_rws.push_back(true);
	}
	else if(starts_tgt && !ends_tgt){
	  can_rws.push_back(true);
	}
	else{
	  can_rws.push_back(false);
	}
      }
    }  

    return can_rws;
  }
  double TargetAttenuationReweighter::calculateWeight(const InteractionChainData& aa){
    std::string mode(getenv("MODE"));
    double wgt =1.0;
    
    MakeReweight*  makerew =  MakeReweight::getInstance();
    bool domipp = (makerew->cv_rw)->doMIPPNumi;

    //Some constants:
    const double graphite_density=1.78;// g/cc
    const double avog_x_mb2cm2 = 6.02214e-4; //useful constant
    const double graphite_A    = 12.01; // g/mole    

    std::vector<InteractionData> vec_inter = aa.interaction_chain;

    //Finding the Data and MC total cross sections:    
    AttenuationMC* dtH = AttenuationMC::getInstance();
    MIPPNumiYieldsBins*  MIPPbins =  MIPPNumiYieldsBins::getInstance();
   
    bool there_is_MIPP  = false;
    bool it_is_survival = false;
    
    //MIPP:
    TargetData tar = aa.tar_info;
    int binID = MIPPbins->BinID(tar.Pz,tar.Pt,tar.Tar_pdg);   
    bool is_le = isLE(aa.target_config);
    bool is_me = isME(aa.target_config);
    if( !is_le &&  !is_me ){
      throw std::runtime_error("cannot determine if it's LE or ME beam");
    }

    TH1D* hzpos;
    if(binID>=0){
       if(tar.Tar_pdg == 211 || tar.Tar_pdg ==- 211){
	there_is_MIPP = true;
	if(tar.Tar_pdg == 211 && is_le)hzpos = dtH->hzpostgt_pip_le[binID];
	if(tar.Tar_pdg == 211 && is_me)hzpos = dtH->hzpostgt_pip_me[binID];
	if(tar.Tar_pdg ==-211 && is_le)hzpos = dtH->hzpostgt_pim_le[binID];
	if(tar.Tar_pdg ==-211 && is_me)hzpos = dtH->hzpostgt_pim_me[binID];
       }
       else if(tar.Tar_pdg == 321){
	 int aux_binID = MIPPbins->BinID(tar.Pz,tar.Pt,211);
	 if(aux_binID>=0){
	   there_is_MIPP = true;
	   //I substracted 78 because we store all histos in AttenuationMC but the kaon histograms make 
	   // sense after P>20 GeV/c and we are using bin ID convention here.
	   if(is_le)hzpos = dtH->hzpostgt_kap_le[aux_binID-78];
	   if(is_me)hzpos = dtH->hzpostgt_kap_me[aux_binID-78];
	 }
       }
       else if(tar.Tar_pdg == -321){
	 int aux_binID = MIPPbins->BinID(tar.Pz,tar.Pt,-211);
	 if(aux_binID>=0){
	   there_is_MIPP = true;
	   if(is_le)hzpos = dtH->hzpostgt_kam_le[aux_binID-78];
	   if(is_me)hzpos = dtH->hzpostgt_kam_me[aux_binID-78];
	 }
       }
       else if(tar.Tar_pdg == 130 || tar.Tar_pdg == 310){
	 int aux_binID = MIPPbins->BinID(tar.Pz,tar.Pt,211);
	 if(aux_binID>=0){
	   there_is_MIPP = true;
	   if(is_le)hzpos = dtH->hzpostgt_kap_le[aux_binID-78];
	   if(is_me)hzpos = dtH->hzpostgt_kap_me[aux_binID-78];
	 }
       }
       else{
	 std::cout<<"=> There is an in MIPPNumiYieldsBins"<<std::endl;
       }
    }
    
    //Survival:
    if(!there_is_MIPP && vec_inter[0].Vol!="BudalMonitor" && vec_inter[0].Vol!="TGT1" && vec_inter[0].Vol!="Budal_HFVS" && vec_inter[0].Vol!="Budal_VFHS"){
      it_is_survival = true;
    }    
    if((mode=="REF")||(mode=="OPT")){
      if(!there_is_MIPP && vec_inter[0].Vol!="TargetFinHorizontal" && vec_inter[0].Vol!="TargetNoSplitSegment" && vec_inter[0].Vol!="tCoreLog"){

	it_is_survival = true;
      }
    }
    
    double delta_sigma = 0.0; // sigma_data - sigma_mc
    double ratio_sigma = 0.0; // sigma_data / sigma_mc

    delta_sigma  = prod_prtC_xsec;
    delta_sigma  += qe_prtC_xsec; 
    
    int binpartC = (dtH->hXS_prtC)->FindBin(vec_inter[0].Inc_P);
    double mcval = (dtH->hXS_prtC)->GetBinContent(binpartC);
    if(mcval<1.e-12){
      std::cout<<"This seems especial cases in which the proton interact before the target (Air?). we need to investigate it."<<" Inc_P "<<vec_inter[0].Inc_P<<std::endl;
      return 1.;
    }
    ratio_sigma = delta_sigma / mcval;
    delta_sigma -= mcval;
    
    //Finding the target penetration for the primary proton beam:
   
    double startZ        = targetStartZ(aa.target_config) + shiftPlaylist(aa.playlist);
    double endZ = (vec_inter[0].Vtx)[2];
    double totmatZ = 0.;  //this will be the amount of material passed. 

    if( is_le){
      //check of initial z position:
      if(vec_inter[0].Vol=="BudalMonitor"){
	if(endZ<startZ || endZ>(startZ+2.0))std::cout<<"Potential error of BudalMonitor=> startZ, endZ: "<<startZ<<" "<<endZ<<std::endl;   
      }      
      if((mode=="OPT")||(mode=="REF")){
	if(vec_inter[0].Vol=="TargetFinHorizontal"){
	  if(endZ<startZ || endZ>(startZ+2.0))std::cout<<"Potential error of BudalMonitor=> startZ, endZ: "<<startZ<<" "<<endZ<<std::endl;   
	}
      }
      totmatZ = getTargetPenetrationLE(startZ,endZ,startZ);
     // std::cout<<"TargetAttenuation:: Start Z "<<startZ<<" End Z "<<endZ<<" totMatZ "<<totmatZ<<" "<<vec_inter[0].Vol<<std::endl;
    }
    else if( is_me ){ 
      //check of initial z position:
      if(vec_inter[0].Vol=="Budal_HFVS"){
	if(endZ<startZ || endZ>(startZ+2.4))std::cout<<"Potential error of Budal_HFVS => startZ, endZ: "<<startZ<<" "<<endZ<<std::endl;   
      }
      totmatZ = getTargetPenetrationME(startZ,endZ,startZ);
    }   

    if(totmatZ<0)totmatZ=0;
    
    totmatZ *= avog_x_mb2cm2;
    totmatZ /= graphite_A;
    totmatZ *= graphite_density;

    //Calculating the weight:
    double norm = 1.0;
    if(there_is_MIPP && domipp){
      int nbins = hzpos->GetXaxis()->GetNbins();
      double integral_mc = 0;
      double integral_wt = 0;
      for(int ib=1;ib<=nbins;ib++){
	double cont = hzpos->GetBinContent(ib);
	double crr  = hzpos->GetXaxis()->GetBinCenter(ib);
	integral_mc += cont;	
	integral_wt += cont*ratio_sigma*exp(-1.0*crr*avog_x_mb2cm2*graphite_density*delta_sigma/graphite_A);
      }
      if(integral_wt<=0 || integral_mc<=0){
	throw std::runtime_error("yield is zero or negative... check!");
      }     
      norm = integral_mc/integral_wt;
    }
    
    //Calculating the weight for the primary:
    double wgt_pri = 1.0;
    
    if(it_is_survival){
      wgt_pri = exp(-1.0*totmatZ*delta_sigma);
    }
    else{
      wgt_pri = norm*ratio_sigma*exp(-1.0*totmatZ*delta_sigma);
    }
    
    //Calculating the weight for secondaries, tertiaries, etc:
    double wgt_sec = 1.0;
    if(!domipp){
      for(size_t ii=1;ii<vec_inter.size();ii++){
	bool starts_tgt = vec_inter[ii-1].Vol == "BudalMonitor" || vec_inter[ii-1].Vol == "TGT1" || vec_inter[ii-1].Vol == "Budal_HFVS"  || vec_inter[ii-1].Vol == "Budal_VFHS";
	if((mode=="REF")||(mode=="OPT")){
	  starts_tgt = vec_inter[ii-1].Vol == "TargetFinHorizontal" || vec_inter[ii-1].Vol == "TargetNoSplitSegment"|| vec_inter[ii-1].Vol == "tCoreLog";
	}
	if(!starts_tgt)continue;
	bool ends_tgt = vec_inter[ii].Vol   == "BudalMonitor" || vec_inter[ii].Vol   == "TGT1" || vec_inter[ii].Vol   == "Budal_HFVS"  || vec_inter[ii].Vol   == "Budal_VFHS";
	if((mode=="REF")||(mode=="OPT")){
	  ends_tgt = vec_inter[ii].Vol   == "TargetFinHorizontal" || vec_inter[ii].Vol   == "TargetNoSplitSegment" || vec_inter[ii].Vol == "tCoreLog";
	}
	//double totmatR  = 0.0;
	double dsigma   = 0.0;
	double fact_int = 1.0;	
	double fact     = 1.0;	
	//dsigma and fact_int:
	if(vec_inter[ii].Inc_pdg==2212){
	  dsigma   = delta_sigma;
	  binpartC = (dtH->hXS_prtC)->FindBin(vec_inter[ii].Inc_P);
	  mcval    = (dtH->hXS_prtC)->GetBinContent(binpartC);
	  fact_int = (dsigma+mcval)/mcval;
	}
	else if(vec_inter[ii].Inc_pdg==211 || vec_inter[ii].Inc_pdg==-211){
	  dsigma   = delta_sigma_piC_xsec;
	  binpartC = (dtH->hXS_piC)->FindBin(vec_inter[ii].Inc_P);
	  mcval    = (dtH->hXS_piC)->GetBinContent(binpartC);
	  fact_int = (dsigma+mcval)/mcval;
	}
	else if(vec_inter[ii].Inc_pdg == 321){
	  if(vec_inter[ii].Inc_P<2.0) dsigma = delta_sigma_kapC_xsec_lowP;
	  if(vec_inter[ii].Inc_P>=2.0)dsigma = delta_sigma_kapC_xsec_highP;
	  binpartC = (dtH->hXS_kapC)->FindBin(vec_inter[ii].Inc_P);
	  mcval    = (dtH->hXS_kapC)->GetBinContent(binpartC);
	  fact_int = (dsigma+mcval)/mcval;
	}
	else if(vec_inter[ii].Inc_pdg ==-321){
	  if(vec_inter[ii].Inc_P<2.0) dsigma = delta_sigma_kamC_xsec_lowP;
	  if(vec_inter[ii].Inc_P>=2.0)dsigma = delta_sigma_kamC_xsec_highP;
	  binpartC = (dtH->hXS_kamC)->FindBin(vec_inter[ii].Inc_P);
	  mcval    = (dtH->hXS_kamC)->GetBinContent(binpartC);
	  fact_int = (dsigma+mcval)/mcval;
	}
	else if(totmatZ<1.e-12){
	  dsigma = 0.0;
	  fact_int = 1.0;
	}
	else{
	  if(vec_inter[ii].Inc_P<2.0) dsigma = delta_sigma_kapC_xsec_lowP;
	  if(vec_inter[ii].Inc_P>=2.0)dsigma = delta_sigma_kapC_xsec_highP;
	  binpartC = (dtH->hXS_kapC)->FindBin(vec_inter[ii].Inc_P);
	  mcval    = (dtH->hXS_kapC)->GetBinContent(binpartC);
	  fact_int = (dsigma+mcval)/mcval;
	}
       if(mcval<1.e-12){
       
      //std::cout<<"This seems especial cases in which the proton interact before the target (Air?). we need to investigate it."<<" Inc_P "<<vec_inter[0].Inc_P<<std::endl;
      
      return 1.;
      
    }
	//Two cases: 1) ending in th target and leaving the target.
	double zi = (vec_inter[ii-1].Vtx)[2];
	double zf = 0.0;
	if(ends_tgt){	  
	  zf = (vec_inter[ii].Vtx)[2];	  
	  fact = fact_int;
	}
	else if(!ends_tgt){
	  double startpart[3] = {(vec_inter[ii-1].Vtx)[0],(vec_inter[ii-1].Vtx)[1],(vec_inter[ii-1].Vtx)[2]};
	  double momtar[3]    = {tar.Px,tar.Py,tar.Pz};
	  zf  = getZTgtExit(startpart,momtar,is_le,is_me);
	}
	if(is_le) totmatZ = getTargetPenetrationLE(zi,zf,startZ);
	if(is_me) totmatZ = getTargetPenetrationME(zi,zf,startZ);
      
        if(totmatZ<1.e-12){
	  dsigma = 0.0;
	  fact = 1.0;
	}
	
	totmatZ *= avog_x_mb2cm2;
	totmatZ /= graphite_A;
	totmatZ *= graphite_density;	
	wgt_sec *= fact*exp(-1.0*totmatZ*dsigma);
	if(std::isinf(wgt_sec))std::cout<<"BAD "<<zi<<" "<<zf<<" "<<totmatZ<<" "<<startZ<<
	" "<<dsigma<<" "<<fact<<std::endl;
	
      }
    }
   // if(isinf(wgt_pri*wgt_sec))std::cout<<wgt_pri<<" "<<wgt_sec<<" "<<zi<<" "<<zf<<std::endl;
    return wgt = wgt_pri*wgt_sec;
    
  }

  double TargetAttenuationReweighter::targetStartZ(const std::string& tgtcfg){
  
     std::string mode(getenv("MODE"));
    //check this (Leo)
    double z0=-51.1; // position of the upstream edge of the budal monitor in 000z config
    // check to see if we are in LE config and adjust accordingly
    if( isLE(tgtcfg) ) {
      z0=-51.72; // determined by ntuple tomography
      if(mode=="OPT")z0=0.0;//-27.3347;  
      if(mode=="REF")z0=-64.7002;
      // check to see if we are in ME config and adjust accordingly
    }else if( isME(tgtcfg) ){ 
      z0=-143.3; // determined by ntuple tomography
      // compared to LE target: recall, LE target is stuffed into the horn.
    }
    else{
      throw std::runtime_error("cannot determine if it's LE or ME beam");
    }
    // pull out all the numbers, recording them into a string
    std::string just_nums;
    for (size_t i=0; i<tgtcfg.length(); ++i) {
      if(isdigit(tgtcfg[i])) just_nums.push_back(tgtcfg[i]);
      if(tgtcfg[i]=='z') break;
    }
    // now convert the string to a number
    double dz=atof(just_nums.c_str());
    //
    return z0-dz;
  }
  
  bool TargetAttenuationReweighter::isLE(const std::string& tgtcfg){
    return (tgtcfg.find("LE")==0)|| (tgtcfg.find("le")==0) || (tgtcfg.find("Le")==0);
  }


  bool TargetAttenuationReweighter::isME(const std::string& tgtcfg){
    return (tgtcfg.find("ME")==0) || (tgtcfg.find("me")==0) || (tgtcfg.find("Me")==0);
  }

  double TargetAttenuationReweighter::shiftPlaylist(const int ipl){
    //Minerva definies dk2nu.vint[1] as the playlist ID.
    //If this variable is not filled, then the value is filled as -1 in InteractionChainData
    //If the experiment definied dk2nu.vint[1], this function needs to be remade.
    double shift_for_playlist = 0;    
    if( ipl == -1 ){
      shift_for_playlist = 0.;
    }
    else if( ipl == 0 || ipl == 1 ){
      shift_for_playlist = 0.50;
    }
    else if( ipl == 7 || ipl == 10 || ipl == 6 ){
      shift_for_playlist = 0.82;
    }
    else if( ipl == 9 ){
      shift_for_playlist = -0.40;
    }
    else if( ipl == 13){
      shift_for_playlist =  0.83;
    }
    else if( ipl == 5 ){
      shift_for_playlist =  1.15;
    }
    else if( ipl == 2 || ipl == 3 ){
      shift_for_playlist = 0.43;
    }
    else if( ipl == 11 || ipl == 12 ){
      shift_for_playlist = 0.83;
    }
    else if( ipl == 4 ){
      shift_for_playlist = 0.43;
    }
    else if( ipl == 8 ){
      shift_for_playlist = - 0.09;
    }
    
    return shift_for_playlist;
    
  }
  
  double TargetAttenuationReweighter::getTargetPenetrationLE(double z_start, double z_end, double z0_budal){
    /*!
     * Returns the amount of target material penetrated by a particle starting at
     * z_start and interacting or exiting the target at z_end.
     * the upstream edge of the budal monitor must be supplied.
     * All inputs must be in units of cm.
     */

    // fin edges determined by printing their position from inside g4numi. 
    // coordinate system will need to be translated such that the upstream
    // edge of the budal monitor is at the value returned by target
    // start z.
    // coordinate here in cm
    std::string mode(getenv("MODE"));
    const int nfins=48;
    const double us_edges[nfins]={25.4484,   42.1683,   44.1983,   46.2283, 48.2583,
				  50.2883,   52.3183,   54.3483,   56.3783,
				  58.4083,   60.4383,   62.4683,   64.4983,
				  66.5283,   68.5583,   70.5883,   72.6183,
				  74.6483,   76.6783,   78.7083,   80.7383,
				  82.7683,   84.7983,   86.8283,   88.8583,
				  90.8883,   92.9183,   94.9483,   96.9783,
				  99.0083,  101.0383,  103.0683,  105.0983,
				  107.1283, 109.1583,  111.1883,  113.2183,
				  115.2483, 117.2783,  119.3083,  121.3383,
				  123.3683, 125.3983,  127.4283,  129.4583,
				  131.4883, 133.5183,  135.5483};
    const int nfins_ref=49;
    //This is for the DUNE reference design target. A. Bashyal				  
    const double us_edges_ref[nfins_ref]={-64.7002,-46.9176,-44.8909, -42.8609, -40.8309, -38.8009, -36.7709, 
					  -34.7409, -32.7109, -30.6809, -28.6509, -26.6209, -24.5909, 
					  -22.5609, -20.5309, -18.5009, -16.4709, -14.4409, -12.4109, -10.3809, -8.3509, -6.3209, 
					  -4.2909, -2.2609, -0.2309, 1.7991, 3.8291, 5.8591, 7.8891, 9.9191, 11.9491, 13.9791, 16.0091, 
					  18.0391, 20.0691, 22.0991, 24.1291, 26.1591, 28.1891, 30.2191,
					  32.2491, 34.2791, 36.3091, 38.3391, 40.3691, 42.3991, 44.4291, 46.4591, 48.4891};
    
    const int nfins_opt = 119;
    const double us_edges_opt[nfins_opt] = {-27.3347,-0.0054, 4.0917, 6.0876, 8.1176, 10.1476, 12.1776, 14.2076, 16.2376, 18.2676, 20.2976, 22.3276, 24.3576, 26.3876, 28.4176, 30.4476, 32.4776, 34.5076, 36.5376, 38.5676, 40.5976,
					    42.6276, 44.6576, 46.6876, 48.7176, 50.7476, 52.7776, 54.8076, 56.8376, 58.8676, 60.8976, 62.9276, 64.9576, 66.9876, 69.0176, 71.0476, 73.0776, 75.1076, 77.1376, 79.1676, 81.1976,
					    83.2276, 85.2576, 87.2876, 89.3176, 91.3476, 93.3776, 95.4076, 97.4376, 99.4676, 101.498, 103.528, 105.558, 107.588, 109.618, 111.648, 113.678, 115.708, 117.738, 119.768, 121.798,
					    123.828, 125.858, 127.888, 129.918, 131.948, 133.978, 136.008, 138.038, 140.068, 142.098, 144.128, 146.158, 148.188, 150.218, 152.248, 154.278, 156.308, 158.338, 160.368, 162.398,
					    164.428, 166.458, 168.488, 170.518, 172.548, 174.578, 176.608, 178.638, 180.668, 182.698, 184.728, 186.758, 188.788, 190.818, 192.848, 194.878, 196.908, 198.938, 200.968, 202.998,
					    205.028, 207.058, 209.088, 211.118, 213.148, 215.178, 217.208, 219.238, 221.268, 223.298, 225.328, 227.358, 229.388,
					    231.418, 233.448, 235.478, 237.508};
    

   // const int nfins_opt = 1;
   // const double us_edges_opt[nfins_opt] = {0.0};
	
    std::vector<double> vus_edges;
    for(int i = 0;i<nfins;i++){
      vus_edges.push_back(us_edges[i]);
    }
    if(mode=="REF"){
      vus_edges.clear();
      for(int i = 0;i<nfins_ref-1;i++){
  	vus_edges.push_back(us_edges_ref[i]);
      }
    }
    if(mode=="OPT"){
      vus_edges.clear();
      for(int i = 0;i<nfins_opt;i++){ //this was still nfins_opt
  	vus_edges.push_back(us_edges_opt[i]);
      }      
    }
    
    const double budal_us_edge=vus_edges.at(0); // fin[0] is the budal
    double fin_width=2.0;
    //just to be not confused...this is fin thickness....the number taken from GEANT4. //Amit Bashyal
    if(mode=="OPT")fin_width = 1.95;
    if(mode=="REF")fin_width = 2.02;
 
    // budal_us_edge + z_trans = z0_budal
    // z_trans = z0_budal - budal_us_edge
    const double z_trans=z0_budal - budal_us_edge;
    
    // translated upstream edge of budal
    //    const double z_up=budal_us_edge+z_trans;
    
    // now, count the amount of material between z_up and z_start
    // for primary protons this will be zero, but the routine
    // should be written generally enough to handle the case
    // in which a trajectory starts in the target
    // we will end up subtracting this from the material between
    // z_up and z_end to get the material traversed.

    double mat_start=0.0;
    double mat_end = 0.0;

    for(unsigned int ifin=0; ifin<vus_edges.size(); ifin++){
      const double fin_us_edge=vus_edges.at(ifin)+z_trans;
      const double fin_ds_edge=fin_us_edge+fin_width;
      if(z_start<=fin_us_edge) break; // no more material to add up
      else if(z_start<fin_ds_edge){ // z_start in this fin
	mat_start+=(z_start-fin_us_edge);
	break;
      }
      else{ // z_start after this fin (z_start>fin_ds_edge)
	mat_start+=fin_width;
      }
    }
    
    // now do the same thing for z_end
 
    for(unsigned int ifin=0; ifin<vus_edges.size(); ifin++){
      const double fin_us_edge=vus_edges.at(ifin)+z_trans;
      const double fin_ds_edge=fin_us_edge+fin_width;
      if(z_end<=fin_us_edge) break; // no more material to add up
      else if(z_end<fin_ds_edge){ // z_end in this fin
	mat_end+=(z_end-fin_us_edge);
	break;
      }
      else{ // z_end after this fin (z_end>fin_ds_edge)
	mat_end+=fin_width;
      }
    }
    

    double mat_traversed=mat_end-mat_start;
    
   // std::cout<<"material traversed is "<<mat_traversed<<std::endl;
    
    return mat_traversed;
    
  }

  double TargetAttenuationReweighter::getTargetPenetrationME(double z_start, double z_end, double z0_budal){  

        /*!
     * Returns the amount of target material penetrated by a particle starting at
     * z_start and interacting or exiting the target at z_end.
     * the upstream edge of the budal monitor must be supplied.
     * All inputs must be in units of cm.
     */

    // fin edges determined by printing their position from inside g4numi. 
    // coordinate system will need to be translated such that the upstream
    // edge of the budal monitor is at the value returned by target
    // start z.
    // coordinate here in cm
    const int nfins=50;
    const double us_edges[nfins]={ 19.285, 22.185, 25.035, 27.485, 
				   29.935, 32.385, 34.835, 37.285, 
				   39.735, 42.185, 44.635, 47.085,
				   49.535, 51.985, 54.435, 56.885,
				   59.335, 61.785, 64.235, 66.685,
				   69.135, 71.585, 74.035, 76.485,
				   78.935, 81.385, 83.835, 86.285,
				   88.735, 91.185, 93.635, 96.085,
				   98.535, 100.985, 103.435, 105.885,
				   108.335, 110.785, 113.235, 115.685,
				   118.135, 120.585, 123.035, 125.485,
				   127.935, 130.385, 132.835, 135.285,
				   137.735, 140.185};
    
    const double ds_edges[nfins]={21.685, 24.585, 27.435, 29.885,
				  32.335, 34.785, 37.235, 39.685,
				  42.135, 44.585, 47.035, 49.485,
				  51.935, 54.385, 56.835, 59.285,
				  61.735, 64.185, 66.635, 69.085,
				  71.535, 73.985, 76.435, 78.885,
				  81.335, 83.785, 86.235, 88.685,
				  91.135, 93.585, 96.035, 98.485,
				  100.935, 103.385, 105.835, 108.285,
				  110.735, 113.185, 115.635, 118.085,
				  120.535, 122.985, 125.435, 127.885,
				  130.335, 132.785, 135.235, 137.685,
				  140.135, 142.585};
    
    // fin[0] is the first budal (HFVS) (horizontal fin, vertical scan)
    // fin[1] is the second budal (VFHS)
    const double budal_us_edge=us_edges[0]; 

    // budal_us_edge + z_trans = z0_budal
    // z_trans = z0_budal - budal_us_edge
    const double z_trans=z0_budal - budal_us_edge;
    
    // translated upstream edge of budal
    //    const double z_up=budal_us_edge+z_trans;
    
    // now, count the amount of material between z_up and z_start
    // for primary protons this will be zero, but the routine
    // should be written generally enough to handle the case
    // in which a trajectory starts in the target
    // we will end up subtracting this from the material between
    // z_up and z_end to get the material traversed.
    double mat_start=0.0;
    for(int ifin=0; ifin<nfins; ifin++){
      const double fin_us_edge=us_edges[ifin]+z_trans;
      const double fin_ds_edge=ds_edges[ifin]+z_trans;
      if(z_start<=fin_us_edge) break; // no more material to add up
      else if(z_start<fin_ds_edge){ // z_start in this fin
	mat_start+=(z_start-fin_us_edge);
	break;
      }
      else{ // z_start after this fin (z_start>fin_ds_edge)
	mat_start+=(fin_ds_edge-fin_us_edge);
      }
    }
    
    // now do the same thing for z_end
    double mat_end=0.0;
    for(int ifin=0; ifin<nfins; ifin++){
      const double fin_us_edge=us_edges[ifin]+z_trans;
      const double fin_ds_edge=ds_edges[ifin]+z_trans;
      if(z_end<=fin_us_edge) break; // no more material to add up
      else if(z_end<fin_ds_edge){ // z_end in this fin
	mat_end+=(z_end-fin_us_edge);
	break;
      }
      else{ // z_end after this fin (z_end>fin_ds_edge)
	mat_end+=(fin_ds_edge-fin_us_edge);
      }
    }
    double mat_traversed=mat_end-mat_start;
    

    return mat_traversed;

  }
  
  double TargetAttenuationReweighter::getZTgtExit(double pos_start[], double mom_start[], bool leflag, bool meflag){
    //check:
    //    if(mom_start[2]<0)std::cout<<"alert momz<0"<<std::endl;
    //this function approximates the z position of the particle that exit the target
    //for LE assuming 1.5 cm x 0.64 cm. xy view.
    
    //First, find the xy point where the partice leaves the target:
    std::string mode(getenv("MODE"));
    double a = -1;
    double b = -1;
    if(leflag){
      a = 0.32;
      b = 0.75;
    }
    if(meflag){
      a = 0.37;
      b = 5.93;
    }    
    if(mode=="REF"){
      a = 0.52;
      b = 1.34;
    }
    if(mode=="OPT"){
      a = 0.49;  //0.37; These numbers changed after discussion with Laura Fields. 06/30/2016
      b = 1.34; //1.05;      
    }
    
    double shift = -1;
    
    //X:
    double mx  = mom_start[0]/mom_start[2];
    if(mx>0)shift = (a-pos_start[0])/mx;
    else shift = (-1.*a-pos_start[0])/mx;    
    double z_usingx = pos_start[2] + shift;
    
    //Y:
    double my  = mom_start[1]/mom_start[2];
    if(my>0){
      if(leflag)shift = (b-pos_start[1])/my;
      if(meflag)shift = (a-pos_start[1])/my;
    }
    else{
      shift = (-1.*b-pos_start[1])/my;        
    }
    double z_usingy = pos_start[2] + shift;  
    
       
    
    double zout = z_usingx;
    if(z_usingy<z_usingx)zout = z_usingy;
    
    return zout;
    
  }

}
#if defined __clang__
#elif defined __GNUC__
  #pragma GCC diagnostic pop
#endif
