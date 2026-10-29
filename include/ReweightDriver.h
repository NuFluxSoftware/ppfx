#ifndef REWEIGHTDRIVER_H
#define REWEIGHTDRIVER_H

#include "MIPPNumiPionYieldsReweighter.h"
#include "MIPPNumiKaonYieldsReweighter.h"
#include "TargetAttenuationReweighter.h"

#include "AbsorptionICReweighter.h"
#include "AbsorptionDPIPReweighter.h"
#include "AbsorptionDVOLReweighter.h"
#include "NucleonAbsorptionOutOfTargetReweighter.h"
#include "OtherAbsorptionOutOfTargetReweighter.h"

#include "OtherReweighter.h"

#include "ThinTargetpCPionReweighter.h"
#include "ThinTargetpCKaonReweighter.h"
#include "ThinTargetnCPionReweighter.h"
#include "ThinTargetpCNucleonReweighter.h"
#include "ThinTargetpipCpipReweighter.h"                                  
#include "ThinTargetMesonIncidentReweighter.h"

#include "ThinTargetnucleonAReweighter.h"

#include "InteractionChainData.h"

namespace NeutrinoFluxReweight{
    
  /*! \class ReweightDriver
   *  \brief A class to manage and drive the weight calculation procedure. 
   *
   * There is a one-to-one correspondence between ReweightDrivers and universes.
   */

  class ReweightDriver{
  public:  
    /*!
     * the constructor
    */
    ReweightDriver(int iuniv, const ParameterTable& cv_pars, const ParameterTable& univ_pars,std::string fileIn);
    ~ReweightDriver();

    /*!
     * Calculate weight for this event, given a set of central value parameters
     * and a set of randomly varied parameters corresponding to the universe
     * that this ReweightDriver is responsible for.
     *
     * - Call mipp_yields on input chain. record interactions covered with a weight in interaction_covered.
     * - Call attenuation on input chain
     * - Call thin target reweighters. 
     *    -# Call na49 for interactions not yet covered. record in interaction_covered
     *    -# ditto. call MIPP, then theory
     */
    double calculateWeight(const InteractionChainData& icd);
    
    /*!  MIPP NuMI yield weight   */
    double mipp_pion_wgt;
 
     /*!   MIPP NuMI kaons yield weight   */
    double mipp_kaon_wgt;
    
    /*!  Target attenuation weight   */
    double att_wgt;

    /*!  Absorption weight   */
    double abs_ic_wgt;
    double abs_dpip_wgt;
    double abs_dvol_wgt;
    double abs_nucleon_wgt;
    double abs_other_wgt;
    double tot_abs_wgt;
       
    /*! Thin target proton on carbon producing pions weights   */
    double pC_pi_wgt;

    /*! Thin target proton on carbon producing kaons weights   */
    double pC_k_wgt;

    /*! Thin target neutron on carbon producing pions weights   */
    double nC_pi_wgt;
    
    /*! Thin target proton on carbon producing nucleons weights   */
    double pC_nu_wgt;

    /*! nuA    */
    double nuA_wgt;

   /*! NA61 pip incident weights   */
    double pipC_pip_wgt;                                                                
    
    /*! Meson incident weights   */
    double meson_inc_wgt;

 


    
     /*!   Any other hadronic interaction not corrected yet   */
    double other_wgt;
    
    MIPPNumiPionYieldsReweighter* MIPP_NUMI_PION_Universe = nullptr;
    MIPPNumiKaonYieldsReweighter* MIPP_NUMI_KAON_Universe = nullptr;
    TargetAttenuationReweighter* TARG_ATT_Universe;
    AbsorptionICReweighter*   VOL_ABS_IC_Universe;
    AbsorptionDPIPReweighter* VOL_ABS_DPIP_Universe;
    AbsorptionDVOLReweighter* VOL_ABS_DVOL_Universe;
    NucleonAbsorptionOutOfTargetReweighter* VOL_ABS_NUCLEON_Universe;
    OtherAbsorptionOutOfTargetReweighter* VOL_ABS_OTHER_Universe;
    ThinTargetpCPionReweighter* THINTARGET_PC_PION_Universe;
    ThinTargetpCKaonReweighter* THINTARGET_PC_KAON_Universe;
    ThinTargetnCPionReweighter* THINTARGET_NC_PION_Universe;
    ThinTargetpCNucleonReweighter* THINTARGET_PC_NUCLEON_Universe;
    ThinTargetpipCpipReweighter* THINTARGET_pipC_pip_Universe = nullptr;
    ThinTargetMesonIncidentReweighter* THINTARGET_MESON_INCIDENT_Universe;

    ThinTargetnucleonAReweighter* THINTARGET_NUCLEON_A_Universe;

    OtherReweighter* OTHER_Universe;

     //Flag to select the reweighters:
    bool doMIPPNumi;
    bool doNA61;
    bool doBoth;
  private:    
    /*!
     * Configures each of the reweighing tools.
     */
    void Configure();
    /*!
     * Parse the option to select which reweighter use
     */
    void ParseOptions();
    int iUniv;
    const ParameterTable& cvPars;
    const ParameterTable& univPars;

    std::string fileOptions;
    
  };

  
}
#endif
