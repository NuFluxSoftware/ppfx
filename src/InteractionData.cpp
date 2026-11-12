
#include "InteractionData.h"
#include <iostream>
#include <iomanip>

namespace NeutrinoFluxReweight{ 
  
  InteractionData::InteractionData(){
    
    particle = TDatabasePDG::Instance();
    
    gen     = 0;
    Inc_pdg = 0;
    Prod_pdg= 0;
    
    Inc_P     = -1000.;
    Prod_P    = -1000.;
    Inc_Mass  = -1000.;
    Prod_Mass = -1000.;

    for(int i=0; i<4; i++){
      Inc_P4[i]=0;
      Prod_P4[i]=0;
      if(i<3) Vtx[i]=0;
    }
	
    xF      = -1000.;
    Pz      = -1000.;
    Theta   = -1000.;
    Pt      = -1000.;
    
    Ecm     = -1000.;
    Betacm  = -1000.;
    Gammacm = -1000.;
    
    Vol = "NoDefinied";
    
  }
  
  InteractionData::InteractionData(int genid, double incMom[], int incPdg, double prodMom[], int prodPdg, std::string volname, std::string procname, double vtx[]){

    particle = TDatabasePDG::Instance();
    // Z direction along the direction of the incident particle
    // Cos between incMom and prodMom:
    // The units are in GeV
    
    InteractionData::gen      = genid;
    
    InteractionData::Inc_pdg  = incPdg;
    InteractionData::Prod_pdg = prodPdg;
    
    InteractionData::Inc_P  = std::sqrt(incMom[0]*incMom[0] + incMom[1]*incMom[1] +incMom[2]*incMom[2]);
    InteractionData::Prod_P = std::sqrt(prodMom[0]*prodMom[0] + prodMom[1]*prodMom[1] +prodMom[2]*prodMom[2]);
    
    double cos_theta = (incMom[0]*prodMom[0]+incMom[1]*prodMom[1]+incMom[2]*prodMom[2])/(Inc_P*Prod_P);
    double sin_theta = std::sqrt(1.-pow(cos_theta,2.0));
    
    //Theta in rads:  
    InteractionData::Theta = std::acos(cos_theta);
    InteractionData::Pt = Prod_P*sin_theta;
    InteractionData::Pz = Prod_P*cos_theta; 
    
    if(Inc_pdg != 1000010020)Inc_Mass = particle->GetParticle(Inc_pdg)->Mass();
    else {Inc_Mass = 1.875;}
    
    if(Prod_pdg != 1000010020)Prod_Mass = particle->GetParticle(Prod_pdg)->Mass();
    else{Prod_Mass = 1.875;}
     
    InteractionData::Inc_Mass  = particle->GetParticle(Inc_pdg)->Mass();
    InteractionData::Prod_Mass = particle->GetParticle(Prod_pdg)->Mass();

    //Ecm, gamma:
    /**
     * Center of mass of the system of projectile
     * and nucleon (not the nucleus!)
     */
    static const double NUCLEON_MASS = 
      (particle->GetParticle("proton")->Mass()
      + particle->GetParticle("neutron")->Mass())/2;
    static const double NUCLEON_MASS2 = NUCLEON_MASS * NUCLEON_MASS;

    double inc_E_lab = std::sqrt(Inc_P*Inc_P + pow(Inc_Mass,2));
    InteractionData::Ecm       = std::sqrt(pow(Inc_Mass,2) + NUCLEON_MASS2 + 2.*inc_E_lab*NUCLEON_MASS); 
    InteractionData::Betacm    = std::sqrt(pow(inc_E_lab,2)-pow(Inc_Mass,2.0))/(inc_E_lab + NUCLEON_MASS); 
    InteractionData::Gammacm   = 1./std::sqrt(1.-pow(Betacm,2.0));
    
    //xF:
    double prod_E_lab  = std::sqrt(Prod_P*Prod_P + pow(Prod_Mass,2));
    double PL          = Gammacm*(Pz-Betacm*prod_E_lab);  // PL is measured in CM frame
    InteractionData::xF  = PL*2./Ecm;

    //4 momenta:
    Inc_P4[3]=inc_E_lab;
    Prod_P4[3]=prod_E_lab;
    for(int i=0; i<3; i++) {Inc_P4[i]=incMom[i]; Prod_P4[i]=prodMom[i];}

    
    //Volume:
    InteractionData::Vol = volname;

    //Process:
    InteractionData::Proc = procname;

    //Vertex:
    for(int i=0; i<3; i++) Vtx[i]=vtx[i];

  }

  //----------------------------------------------------------------------
  InteractionData::~InteractionData(){
    
  }
  
  std::ostream& InteractionData::print(std::ostream& os) const {
    using namespace std;
    os<<"in:"<<setw(5)<<Inc_pdg
      <<"|p3:";
    for(int i=0; i<3; i++) {
      os<<setiosflags(ios::fixed) << setprecision(2)<<setw(6)<<Inc_P4[i]<<" ";
    }
    os<<"||out:"<<setw(5)<<Prod_pdg
      <<"|p3:"<<setiosflags(ios::fixed) << setprecision(2);
    for(int i=0; i<3; i++) {
      os<<setiosflags(ios::fixed) << setprecision(2)<<setw(6)<<Prod_P4[i]<<" ";
    }
    os <<"|v3:";
    for(int i=0; i<3; i++) {
      os<<setiosflags(ios::fixed) << setprecision(2)<<setw(5)<<Vtx[i]<<" ";
    }
    os<<"xF,pT:"<<xF<<","<<Pt;
    os<<endl;
    return os;
  }

}
