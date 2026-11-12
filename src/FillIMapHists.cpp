#include "FillIMapHists.h"
#include "InteractionChainData.h"
#include "InteractionData.h"
#include "CentralValuesAndUncertainties.h"
#include "MIPPNumiYieldsBins.h"
#include "MIPPNumiMC.h"
#include "ThinTargetBins.h"
#include "ThinTargetpipCpipBins.h"
#include "ThinTargetpipCpipMC.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include <iostream>
#include <math.h>
#include "MakeReweight.h"

double FillIMapHists(TChain* tdk2nu, TChain* tdkmeta, HistList* hists, const FillIMapHistsOpts* opts){
  using namespace NeutrinoFluxReweight;
  
  // setup the event loop, filling Dk2Nu and DkMeta objects
  bsim::Dk2Nu*  dk2nu  = new bsim::Dk2Nu;  
  bsim::DkMeta* dkmeta = new bsim::DkMeta;
  tdk2nu->SetBranchAddress("dk2nu",&dk2nu);
  Long64_t nentries  = tdk2nu->GetEntries();
  Long64_t ntrees    = tdk2nu->GetNtrees();
  tdkmeta->SetBranchAddress("dkmeta",&dkmeta);

  ////////////////////// initializing reweighters ///////////////////////////
  //
  // here, we need to initialize the reweighers
  // and pass them into FillOneEntry
  // so it can use them to determine if a particular interaction
  // or chain of interactions can be reweighted

  ///// Some inputs file are needed to ensure correct operation
  ///// of the reweight drivers and also, some reweigthwers 
  /////(for instance: ThinTargetnCPionReweighter) use other reweighters. 
  ///// The singleton MakeReweight makes this initialization and pass 
  ///// th information between reweighters.

  const char* ppfxDir = getenv("PPFX_DIR");
  MakeReweight* makerew = MakeReweight::getInstance();
  makerew->SetOptions(Form("%s/scripts/inputs_imap.xml",ppfxDir));    // restored for MINERvA/MIPP comparison, see scripts/README.md

  FillIMapHistsReweighters reweighters;
  reweighters.NumiPions               = (makerew->cv_rw)->MIPP_NUMI_PION_Universe;
  reweighters.NumiKaons               = (makerew->cv_rw)->MIPP_NUMI_KAON_Universe;
  reweighters.ThinTargetpCPion        = (makerew->cv_rw)->THINTARGET_PC_PION_Universe;
  reweighters.ThinTargetpCKaon        = (makerew->cv_rw)->THINTARGET_PC_KAON_Universe;
  reweighters.ThinTargetnCPion        = (makerew->cv_rw)->THINTARGET_NC_PION_Universe;
  reweighters.ThinTargetpCNucleon     = (makerew->cv_rw)->THINTARGET_PC_NUCLEON_Universe;
  reweighters.ThinTargetMesonIncident = (makerew->cv_rw)->THINTARGET_MESON_INCIDENT_Universe;
  reweighters.ThinTargetpipCpip       = (makerew->cv_rw)->THINTARGET_pipC_pip_Universe;
  reweighters.ThinTargetnucleonA      = (makerew->cv_rw)->THINTARGET_NUCLEON_A_Universe;

  std::cout<<"FillIMapHists looping over "<<ntrees<<" trees with a total of "<<nentries<<" entries."<<std::endl;
  double total_weight=0.0;
  for(Long64_t ientry=0;ientry<nentries;ientry++){
    if(ientry%100000==0)std::cout<<"ientry "<<ientry/1000<<" k evts"<<std::endl;
    tdk2nu->GetEntry(ientry);
    tdkmeta->GetEntry(ientry);    
    
    total_weight+=FillOneEntry(dk2nu,dkmeta,hists,opts,&reweighters); 
    // std::cout<<"tot wgt: "<<total_weight<<" "<<dk2nu->decay.ntype<<std::endl;
  }
  //Releasing memory:
  makerew->resetInstance();
  
  return total_weight;
}

//#define DEBUG

double FillOneEntry(bsim::Dk2Nu* dk2nu, bsim::DkMeta* dkmeta, HistList* hists, const FillIMapHistsOpts* opts, FillIMapHistsReweighters* reweighters){
  double weight=0.0;
  TDatabasePDG* pdg = TDatabasePDG::Instance();
  // check that the neutrino is of the requested type and that
  // the energy is within range
  const int nu_type=dk2nu->decay.ntype;
  const int nuray_idx=1; // this corresponds to the location of minerva
  const double enu=dk2nu->nuray[nuray_idx].E; // energy at MINOS ND

#ifdef DEBUG
  std::cout<<"FillOneEntry() for nu_type= "<<nu_type
	   <<" and energy "<<enu<<std::endl;
#endif
  // setting opts.nuid=0 will result in all neutrino species being plotted
  if( (opts->nuid!=0 && nu_type!= opts->nuid) 
      || (enu<opts->elow || enu>opts->ehigh) ){
#ifdef DEBUG
    std::cout<<"Fails cut on nu_type or energy"<<std::endl;
#endif 
    return 0;
  }
  const double nwtnear=dk2nu->nuray[nuray_idx].wgt;
  const double nimpwt=dk2nu->decay.nimpwt;
  weight=nwtnear*nimpwt;
  
  hists->_h_nuflux->Fill(enu,weight/pival);
  
  /*
  if(isnan(weight)){
    std::cout<<"Encountered a NaN weight, dk2nu follows"<<std::endl;
    std::cout<<(*dk2nu)<<std::endl;    
  }
  */
  NeutrinoFluxReweight::InteractionChainData icd(dk2nu,dkmeta);
  const int ninter=icd.interaction_chain.size();
  std::vector<bool> numi_pion_nodes=reweighters->NumiPions->canReweight(icd);
  std::vector<bool> numi_kaon_nodes=reweighters->NumiKaons->canReweight(icd);

#ifdef DEBUG
  std::cout<<"Passes energy cut and has a weight of "<<weight
	   <<" with "<<ninter<<" entries in ancestry chain"<<std::endl;
#endif 

  int ninter_all=0; // a variable to count all non-Decay interactions
  int ninter_cuts=0;// ... and only those passing the MIPP/NA49/etc cuts
  
  for(int iinter=0; iinter<ninter; iinter++){
    
    const NeutrinoFluxReweight::InteractionData& interdata
      =icd.interaction_chain[iinter];
#ifdef DEBUG
    std::cout<<"Processing interaction "<<iinter<<endl;
    interdata.print(std::cout);
#endif 
    // check to see if this entry is a decay
    if(interdata.Proc=="Decay"){
#ifdef DEBUG
      std::cout<<"   This is a decay, skip it"<<std::endl;
#endif
      continue; // if so, don't histogram
    }
    ninter_all++;
    /////////////////////////////////////////////////////////////////////
    // check here if this interaction is covered by NA49, MIPP, etc
    /////////////////////////////////////////////////////////////////////
    
    if(opts->cut_mipp && numi_pion_nodes[iinter]) continue;
    if(opts->cut_mipp && numi_kaon_nodes[iinter]) continue;
    // Thin target reweighters are based on data and theoretical motivated data extensions.
    // bool covered_by_na61 = false; // unused
    // ThinTargetpipCpip is only constructed when the options file selects NA61On
    // (ReweightDriver::Configure, doNA61 gate) -- null in MIPPNuMIOn mode, restored here
    if(reweighters->ThinTargetpipCpip && reweighters->ThinTargetpipCpip->canReweight(interdata)){
       // covered_by_na61 = true; // unused
       if(! opts->cut_na61) hists->_h_aveint_vs_enu_pipCpip->Fill(enu,weight);
                     }
     
    bool covered_by_thintarget = false;
    if(reweighters->ThinTargetpCPion->canReweight(interdata)){
      covered_by_thintarget = true;
      if(! opts->cut_thintarget) hists->_h_aveint_vs_enu_thin_pCpion->Fill(enu,weight);
      if(interdata.Prod_pdg==211){
	hists->_h_occ_xfpt_pc_pip->Fill(interdata.xF,interdata.Pt,weight);
	double hpweight=reweighters->ThinTargetpCPion->calculateWeight(interdata);
	hists->_h_hpwgt_xfpt_pc_pip->Fill(interdata.xF,interdata.Pt,weight*hpweight);	
      }
    }
    else if(reweighters->ThinTargetpCKaon->canReweight(interdata)){
      covered_by_thintarget = true;
      if(! opts->cut_thintarget) hists->_h_aveint_vs_enu_thin_pCkaon->Fill(enu,weight);
      if(interdata.Prod_pdg==321){
	hists->_h_occ_xfpt_pc_kp->Fill(interdata.xF,interdata.Pt,weight);
	double hpweight=reweighters->ThinTargetpCKaon->calculateWeight(interdata);
	hists->_h_hpwgt_xfpt_pc_kp->Fill(interdata.xF,interdata.Pt,weight*hpweight);	
      }

    }
    else if(reweighters->ThinTargetnCPion->canReweight(interdata)){
      covered_by_thintarget = true;
      if(! opts->cut_thintarget) hists->_h_aveint_vs_enu_thin_nCpion->Fill(enu,weight);
    }
    else if(reweighters->ThinTargetpCNucleon->canReweight(interdata)){
      covered_by_thintarget = true;
      if(! opts->cut_thintarget) hists->_h_aveint_vs_enu_thin_pCnucleon->Fill(enu,weight);
    }
    else if(reweighters->ThinTargetMesonIncident->canReweight(interdata)){
      covered_by_thintarget = true;
      if(! opts->cut_thintarget) hists->_h_aveint_vs_enu_thin_mesoninc->Fill(enu,weight);
    }
     //else if(reweighters->ThinTargetpipCpip->canReweight(interdata)){
      //covered_by_thintarget = true;
      //if(! opts->cut_thintarget) hists->_h_aveint_vs_enu_pipCpip->Fill(enu,weight);
    //}
    else if(reweighters->ThinTargetnucleonA->canReweight(interdata)){
      covered_by_thintarget = true; //Amit Changed this...default was true.
      if(! opts->cut_thintarget) hists->_h_aveint_vs_enu_thin_nucleona->Fill(enu,weight);
    }
    else{
      covered_by_thintarget = false; 
      hists->_h_aveint_vs_enu_others->Fill(enu,weight);
    }
    
    if(! opts->cut_thintarget)hists->_h_aveint_vs_enu_tot->Fill(enu,weight);
    else{
	if(!covered_by_thintarget)hists->_h_aveint_vs_enu_tot->Fill(enu,weight);
    }
    
    if(opts->cut_na61 && covered_by_thintarget) continue;                                 //BHUMIKA, replaced mipp by na61
    
 
    ninter_cuts++;
    // get an index into the large arrays listing the volume names
    // and the material of each volume.
    int mv_idx=FindIndexFromVolume(interdata.Vol);
    std::string mode = getenv("MODE");
    std::string volume_name = IMap::volume[mv_idx];
    std::string material_name = IMap::materials[mv_idx];
    
    if(mode=="REF"||mode=="OPT"){
      volume_name = IMap::volumedune[mv_idx];
      material_name = IMap::materialsdune[mv_idx];
    }
    
    if(mv_idx==-1){
      std::cout<<"Skipping unknown volume "<< interdata.Vol
	       <<" for interaction "<<iinter<<std::endl;
    }
    
    // fill a 2D histogram of projectile vs. material
    const string proj_name=pdg->GetParticle(interdata.Inc_pdg)->GetName();
    const string prod_name=pdg->GetParticle(interdata.Prod_pdg)->GetName();
    //if(covered_by_thintarget) //Uncomment this to get the thin target coverage after turning off the nucleonA reweighter off above.
    hists->_h_in_vs_mat->Fill(material_name.c_str(),proj_name.c_str(),weight);
    // figure out if the produced particle is one that we want
    // to record in histograms
    // The list of such particles is in IMap::popparticle
    // the strange name is apparently a contraction: "popular particles"
    const int prod_pop_idx=FindIndexFromParticleName(prod_name);
    const int proj_pop_idx=FindIndexFromParticleName(proj_name);
#ifdef DEBUG
    std::cout<<"   Projectile: "<<proj_name<<" with popidx "<<proj_pop_idx<<std::endl;
    std::cout<<"   Produced  : "<<prod_name<<" with popidx "<<prod_pop_idx<<std::endl;
#endif
    

    // look at things from the produced particles standpoint
    if(prod_pop_idx!=-1){ // for each of the commonly produced particles.

      // histogram kinetic energy, 3-momentum, and xF,pT
      const double produced_KE=interdata.Prod_P4[3]-interdata.Prod_Mass;
      hists->_hkepop_tot[prod_pop_idx]->Fill(produced_KE,weight);
      hists->_htmpop_tot[prod_pop_idx]->Fill(interdata.Prod_P,weight);
      hists->_hxfpt_tot[prod_pop_idx]->Fill(interdata.xF,interdata.Pt,weight);

      // histogram the material that the interaction occured in
      // along with the projectile that made the particle in question
      hists->_hmatbkw[prod_pop_idx]->Fill(material_name.c_str(),proj_name.c_str(),weight);
      
      // now, dig deeper
      if(proj_pop_idx!=-1){ // for each of the common *projectiles*
	// histogram kinetic energy, 3-momentum, and xF,pT
	// of the produced particle
	hists->_hkepop[prod_pop_idx][proj_pop_idx]->Fill(produced_KE,weight);
	hists->_htmpop[prod_pop_idx][proj_pop_idx]->Fill(interdata.Prod_P,weight);
	hists->_hxfpt[prod_pop_idx][proj_pop_idx]->Fill(interdata.xF,interdata.Pt,weight);
      }
    }
    
    // now look at things from the projectile's standpoint
    if(proj_pop_idx!=-1){ // for each of the common projectiles
      // histogram the kinetic energy of the projectile
      const double projectile_KE=interdata.Inc_P4[3]-interdata.Inc_Mass;
      hists->_henergytotal[proj_pop_idx]->Fill(projectile_KE,weight);
      // histogram the volume/material and the produced particle
      hists->_hmat[proj_pop_idx]->Fill(material_name.c_str(),prod_name.c_str(),weight);
      hists->_hvol[proj_pop_idx]->Fill(volume_name.c_str(),prod_name.c_str(),weight);
      // histogram the energy of the projectile for each volume
      // This may be overkill!
      if(mv_idx!=-1) hists->_henergyvolume[mv_idx][proj_pop_idx]->Fill(projectile_KE,weight);
      if(projectile_KE>118 and proj_pop_idx==1){
	std::cout<<"Oh noes!"<<std::endl;
      }
    }

  }

  // now fill the # of interactions vs enu
  hists->_h_nint_vs_enu->Fill(enu,ninter_all,weight);
  hists->_h_nint_vs_enu_cuts->Fill(enu,ninter_cuts,weight);
  
  // the end
  return weight;
}


int FindIndexFromVolume(const std::string& wanted){
  std::string mode(getenv("MODE"));
  if((mode=="REF"||mode=="OPT")){
  for(int i=0; i<IMap::nvoldune; i++){
    if(wanted== std::string(IMap::volumedune[i])) return i;
  }
  }
  else{
  for(int i=0; i<IMap::nvol; i++){
    if(wanted== std::string(IMap::volume[i])) return i;
  }  
  }
  return -1;
}

int FindIndexFromParticleName(const std::string& wanted){
  for(int i=0; i<IMap::npop; i++){
    if(wanted== std::string(IMap::popparticle[i])) return i;
  }
  return -1;
}

