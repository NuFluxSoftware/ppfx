#ifndef THINTARGETPIPCPIPMC_H
#define THINTARGETPIPCPIPMC_H

#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <iostream>
#include <sstream>

namespace NeutrinoFluxReweight{
  
  /* classThinTargetpipCpipMC brief A class to manage the MC value for  NA61 60GeV pip incident NuMI. 
   */
  class ThinTargetpipCpipMC{

  private:
   ThinTargetpipCpipMC();

  public:
    static ThinTargetpipCpipMC* getInstance();
    
    //! Read a xml file name to get the mc value for pip
    void pipC_pip_mc_from_xml(const char* filename);
    void pipC_pim_mc_from_xml(const char* filename);
    void pipC_kp_mc_from_xml(const char* filename); 
    void pipC_km_mc_from_xml(const char* filename);
    void pipC_p_mc_from_xml(const char* filename);    
    void pipC_k0s_mc_from_xml(const char* filename);
    void pipC_lam_mc_from_xml(const char* filename);
    void pipC_alam_mc_from_xml(const char* filename);
       
       //! MC value for this HP production
    double getMCval(double Prod_P,double Theta, int Inc_pdg, int Prod_pdg);

  private:    
    std::vector<double> pipC_pip_cv;
    std::vector<double> pipC_pim_cv;
    std::vector<double> pipC_kp_cv;
    std::vector<double> pipC_km_cv;
    std::vector<double> pipC_p_cv; 
    std::vector<double> pipC_k0s_cv; 
    std::vector<double> pipC_lam_cv;
    std::vector<double> pipC_alam_cv;

    
    std::vector<double> v_pip_pmin,v_pip_pmax,v_pip_thetamin,v_pip_thetamax;
    std::vector<double> v_pim_pmin,v_pim_pmax,v_pim_thetamin,v_pim_thetamax;
    std::vector<double> v_kp_pmin,v_kp_pmax,v_kp_thetamin,v_kp_thetamax;
    std::vector<double> v_p_pmin,v_p_pmax,v_p_thetamin,v_p_thetamax;
    std::vector<double> v_km_pmin,v_km_pmax,v_km_thetamin,v_km_thetamax;
    std::vector<double> v_k0s_pmin,v_k0s_pmax,v_k0s_thetamin,v_k0s_thetamax;
    std::vector<double> v_lam_pmin,v_lam_pmax,v_lam_thetamin,v_lam_thetamax;
    std::vector<double> v_alam_pmin,v_alam_pmax,v_alam_thetamin,v_alam_thetamax;
    bool ranges_pip_already_filled;
    bool ranges_pim_already_filled;
    bool ranges_kp_already_filled;
    bool ranges_km_already_filled;
    bool ranges_k0s_already_filled;
    bool ranges_p_already_filled;
    bool ranges_lam_already_filled;
    bool ranges_alam_already_filled;
    bool ranges_already_filled;
    //double proton_no_interacting; // unused

    static ThinTargetpipCpipMC* instance;
    
  };

  
  
}
#endif
