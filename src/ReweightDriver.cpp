#include "ReweightDriver.h"
#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace NeutrinoFluxReweight{
  
  ReweightDriver::ReweightDriver(int iuniv, const ParameterTable& cv_pars, const ParameterTable& univ_pars, std::string fileIn)
    : iUniv(iuniv), cvPars(cv_pars), univPars(univ_pars), fileOptions(fileIn)
  {
    ParseOptions();
    Configure();
    
  }
  
  void ReweightDriver::Configure(){
    
    //Creating the vector of reweighters:
   //This script will select which reweighters will be used. Reads input variables provided and proceeds accordingly.
    
    if(doMIPPNumi){
      MIPP_NUMI_PION_Universe = new MIPPNumiPionYieldsReweighter(iUniv,cvPars,univPars);
      MIPP_NUMI_KAON_Universe = new MIPPNumiKaonYieldsReweighter(iUniv,cvPars,univPars);
    }
    if(doNA61){
      THINTARGET_pipC_pip_Universe = new ThinTargetpipCpipReweighter(iUniv,cvPars,univPars);                       
    }
    TARG_ATT_Universe = new TargetAttenuationReweighter(iUniv,cvPars,univPars);    
    VOL_ABS_IC_Universe = new AbsorptionICReweighter(iUniv,cvPars,univPars);
    VOL_ABS_DPIP_Universe = new AbsorptionDPIPReweighter(iUniv,cvPars,univPars);
    VOL_ABS_DVOL_Universe = new AbsorptionDVOLReweighter(iUniv,cvPars,univPars);
    VOL_ABS_NUCLEON_Universe = new NucleonAbsorptionOutOfTargetReweighter(iUniv,cvPars,univPars);
    VOL_ABS_OTHER_Universe = new OtherAbsorptionOutOfTargetReweighter(iUniv,cvPars,univPars);

    THINTARGET_PC_PION_Universe = new ThinTargetpCPionReweighter(iUniv,cvPars,univPars);
    THINTARGET_PC_KAON_Universe = new ThinTargetpCKaonReweighter(iUniv,cvPars,univPars);
    THINTARGET_NC_PION_Universe = new ThinTargetnCPionReweighter(iUniv,cvPars,univPars);
    THINTARGET_PC_NUCLEON_Universe = new ThinTargetpCNucleonReweighter(iUniv,cvPars,univPars);

    THINTARGET_MESON_INCIDENT_Universe = new ThinTargetMesonIncidentReweighter(iUniv,cvPars,univPars);

    THINTARGET_NUCLEON_A_Universe = new ThinTargetnucleonAReweighter(iUniv,cvPars,univPars);
    OTHER_Universe = new OtherReweighter(iUniv,cvPars,univPars);
    
  }
  
  void ReweightDriver::ParseOptions(){
    //Parsing the file input:
    using boost::property_tree::ptree;
    ptree top;
    std::string val = "";
    read_xml(fileOptions.c_str(),top,2); // option 2 removes comment strings
    ptree& options = top.get_child("inputs.Settings");
    
    val = options.get<std::string>("Reweighters");
    if(val=="MIPPNuMIOn"){
                doMIPPNumi = true; 
                doNA61 = false;
                  }
    else if(val=="NA61On"){doMIPPNumi = false; doNA61=true; 
         // std::cout<<"NA61 is on and the other (MIPP) is OFF"<<std::endl;    
}
    else if(val== "BothOff"){doMIPPNumi = false; doNA61=false;
         std::cout<<"NA61 and MIPP are OFF, only the NA49 weights would act now"<<std::endl;

}
   
    else {doMIPPNumi = false;doNA61=false;}
    
  }
  double ReweightDriver::calculateWeight(const InteractionChainData& icd){

    double tot_wgt = 1.0;
    
    //Boolean flags: 
    const int nnodes=icd.interaction_chain.size();
    std::vector<bool> interaction_nodes(nnodes,false);
    std::vector<bool> attenuation_nodes(nnodes,false);
    std::vector<bool> absorption_nodes(nnodes,false);

    /// ----- PROCESS INTERACTION NODES ----- ///
    
    //MIPP NuMI Pions:
    bool has_mipp = false;
    mipp_pion_wgt = 1.0;
    if(doMIPPNumi){
      interaction_nodes = MIPP_NUMI_PION_Universe->canReweight(icd);
      for(size_t ii=0;ii<interaction_nodes.size();ii++){
	if(interaction_nodes[ii]==true){
	  has_mipp = true;
	  mipp_pion_wgt = MIPP_NUMI_PION_Universe->calculateWeight(icd);
	  break; 
	}
      }
      tot_wgt *= mipp_pion_wgt;
    }
    
    //MIPP NuMI Kaons:
    mipp_kaon_wgt = 1.0;
    if(!has_mipp && doMIPPNumi){
      interaction_nodes = MIPP_NUMI_KAON_Universe->canReweight(icd);

      for(size_t ii=0;ii<interaction_nodes.size();ii++){
	if(interaction_nodes[ii]==true){
	  has_mipp = true;
	  mipp_kaon_wgt = MIPP_NUMI_KAON_Universe->calculateWeight(icd);
	  break; 
	}
      }
      tot_wgt *= mipp_kaon_wgt;
    }

//Thin Target pipC-->pip Interaction:
  pipC_pip_wgt = 1.0;
  if(doNA61)
   {
    for(int ii=(interaction_nodes.size()-1);ii>=0;ii--){
      if(interaction_nodes[ii]==false){
        bool is_rew = THINTARGET_pipC_pip_Universe->canReweight((icd.interaction_chain)[ii]);
   
        if(is_rew){
          double rewval = THINTARGET_pipC_pip_Universe->calculateWeight((icd.interaction_chain)[ii]);
          pipC_pip_wgt *= rewval;
          interaction_nodes[ii]=true;
        }
      }
    }

   tot_wgt *= pipC_pip_wgt;
 
}   
 //Thin Target pC->piX:
    pC_pi_wgt = 1.0;
    for(int ii=(interaction_nodes.size()-1);ii>=0;ii--){	
      if(interaction_nodes[ii]==false){
	bool is_rew = THINTARGET_PC_PION_Universe->canReweight((icd.interaction_chain)[ii]);
	if(is_rew){
	  double rewval = THINTARGET_PC_PION_Universe->calculateWeight((icd.interaction_chain)[ii]);
	  pC_pi_wgt *= rewval;
	  interaction_nodes[ii]=true;
	}
      }
    }
    tot_wgt *= pC_pi_wgt;
       
    //Thin Target pC->KX:
    pC_k_wgt = 1.0;    
    for(int ii=(interaction_nodes.size()-1);ii>=0;ii--){	
      if(interaction_nodes[ii]==false){
	bool is_rew = THINTARGET_PC_KAON_Universe->canReweight((icd.interaction_chain)[ii]);
	if(is_rew){
	  double rewval = THINTARGET_PC_KAON_Universe->calculateWeight((icd.interaction_chain)[ii]);
	  pC_k_wgt *= rewval;
	  interaction_nodes[ii]=true;
	}
      }
    }
    tot_wgt *= pC_k_wgt;
    
    //Thin Target nC->piX:
    nC_pi_wgt = 1.0;
    for(int ii=(interaction_nodes.size()-1);ii>=0;ii--){	
      if(interaction_nodes[ii]==false){
	bool is_rew = THINTARGET_NC_PION_Universe->canReweight((icd.interaction_chain)[ii]);
	if(is_rew){
	  double rewval = THINTARGET_NC_PION_Universe->calculateWeight((icd.interaction_chain)[ii]);
	  nC_pi_wgt *= rewval;
	  interaction_nodes[ii]=true;
	}
      }
    }
    tot_wgt *= nC_pi_wgt;
    
    //Thin Target pC->nucleonX:
    pC_nu_wgt = 1.0;
    for(int ii=(interaction_nodes.size()-1);ii>=0;ii--){	
      if(interaction_nodes[ii]==false){
	bool is_rew = THINTARGET_PC_NUCLEON_Universe->canReweight((icd.interaction_chain)[ii]);
	if(is_rew){
	  double rewval = THINTARGET_PC_NUCLEON_Universe->calculateWeight((icd.interaction_chain)[ii]);
	  pC_nu_wgt *= rewval;
	  interaction_nodes[ii]=true;
	}
      }
    }
    tot_wgt *= pC_nu_wgt;

  

    
    //Thin Target Meson Incident:
    meson_inc_wgt = 1.0; 
    for(int ii=(interaction_nodes.size()-1);ii>=0;ii--){	
      if(interaction_nodes[ii]==false){
	bool is_rew = THINTARGET_MESON_INCIDENT_Universe->canReweight((icd.interaction_chain)[ii]);
	if(is_rew){
	  double rewval = THINTARGET_MESON_INCIDENT_Universe->calculateWeight((icd.interaction_chain)[ii]);
	  meson_inc_wgt *= rewval;
	  interaction_nodes[ii]=true;
	}
      }
    }
    tot_wgt *= meson_inc_wgt;    
 
 

  
   
    //Thin Target Nucleon Incident not hanldle NA49 or Barton:
    nuA_wgt = 1.0; 
    for(int ii=(interaction_nodes.size()-1);ii>=0;ii--){	
      if(interaction_nodes[ii]==false){
	bool is_rew = THINTARGET_NUCLEON_A_Universe->canReweight((icd.interaction_chain)[ii]);
	if(is_rew){
	  double rewval = THINTARGET_NUCLEON_A_Universe->calculateWeight((icd.interaction_chain)[ii]);
	  nuA_wgt *= rewval;
	  interaction_nodes[ii]=true;
	}
      }
    }
    tot_wgt *= nuA_wgt;
    
    //Any other interaction not handled yet:
    other_wgt = 1.0;
    for(int ii=(interaction_nodes.size()-1);ii>=0;ii--){	
      if(interaction_nodes[ii]==false){
	bool is_rew = OTHER_Universe->canReweight((icd.interaction_chain)[ii]);
	if(is_rew){
	  double rewval = OTHER_Universe->calculateWeight((icd.interaction_chain)[ii]);
	  other_wgt *= rewval;
	  interaction_nodes[ii]=true;
	}
      }
    }
    tot_wgt *= other_wgt;
    
    //Target attenuation correction:
    att_wgt = 1.0;
    attenuation_nodes = TARG_ATT_Universe->canReweight(icd);   
    //we just see for the first position (prmary proton)
    if(attenuation_nodes.size()>0 && attenuation_nodes[0]==true){
      att_wgt *= TARG_ATT_Universe->calculateWeight(icd);
    }
    tot_wgt *= att_wgt;
        
    //Absorption correction:
    tot_abs_wgt = 1.0;
    
    // Correction of the pi & K absorption in volumes (Al)    
    abs_ic_wgt = 1.0;    
    absorption_nodes = VOL_ABS_IC_Universe->canReweight(icd);
    //std::cout<<"size of absorption_nodes is "<<absorption_nodes.size()<<std::endl;

    if(absorption_nodes.size()>0 && absorption_nodes[0]==true){
      abs_ic_wgt *= VOL_ABS_IC_Universe->calculateWeight(icd);
    }
    tot_wgt     *= abs_ic_wgt;
    tot_abs_wgt *= abs_ic_wgt;
 
    //Correction of the pi & K absorption in volumes (Fe)
    abs_dpip_wgt = 1.0;    
    absorption_nodes = VOL_ABS_DPIP_Universe->canReweight(icd);
    if(absorption_nodes.size()>0 && absorption_nodes[0]==true){
      abs_dpip_wgt *= VOL_ABS_DPIP_Universe->calculateWeight(icd);
    }
    tot_wgt     *= abs_dpip_wgt;
    tot_abs_wgt *= abs_dpip_wgt;
    
    
    //Correction of the pi & K absorption in volumes (He)
    abs_dvol_wgt = 1.0;    
    absorption_nodes = VOL_ABS_DVOL_Universe->canReweight(icd);
    if(absorption_nodes.size()>0 && absorption_nodes[0]==true){
      abs_dvol_wgt *= VOL_ABS_DVOL_Universe->calculateWeight(icd);
    }
    tot_wgt     *= abs_dvol_wgt;
    tot_abs_wgt *= abs_dvol_wgt;
    
    //Correction of nucleons on Al, Fe and He.
    abs_nucleon_wgt = 1.0;    
    absorption_nodes = VOL_ABS_NUCLEON_Universe->canReweight(icd);
    if(absorption_nodes.size()>0 && absorption_nodes[0]==true){
      abs_nucleon_wgt *= VOL_ABS_NUCLEON_Universe->calculateWeight(icd);
    }
    tot_wgt     *= abs_nucleon_wgt;
    tot_abs_wgt *= abs_nucleon_wgt;
    
    //Correction of any other particle on Al, Fe and He.
    abs_other_wgt = 1.0;    
    absorption_nodes = VOL_ABS_OTHER_Universe->canReweight(icd);
    if(absorption_nodes.size()>0 && absorption_nodes[0]==true){
      abs_other_wgt *= VOL_ABS_OTHER_Universe->calculateWeight(icd);
    }
    tot_wgt     *= abs_other_wgt;
    tot_abs_wgt *= abs_other_wgt;
    
    if(tot_wgt!=tot_wgt){
      std::cout<<"Alert nan total wgt... check!!!"<<std::endl;
      return 1.0;
    }
    return tot_wgt;
  }
  
  ReweightDriver::~ReweightDriver(){
 
    if(doMIPPNumi){
      delete MIPP_NUMI_PION_Universe;
      delete MIPP_NUMI_KAON_Universe;
    }   
   if (doNA61){
      delete THINTARGET_pipC_pip_Universe;

} 
    delete TARG_ATT_Universe;
    delete VOL_ABS_IC_Universe;
    delete VOL_ABS_DPIP_Universe;
    delete VOL_ABS_DVOL_Universe;
    delete VOL_ABS_NUCLEON_Universe;
    delete VOL_ABS_OTHER_Universe;
    delete THINTARGET_PC_PION_Universe;
    delete THINTARGET_PC_KAON_Universe;
    delete THINTARGET_NC_PION_Universe;
    delete THINTARGET_PC_NUCLEON_Universe;
    delete THINTARGET_MESON_INCIDENT_Universe;
    delete THINTARGET_NUCLEON_A_Universe;
    delete OTHER_Universe;
  }

}
