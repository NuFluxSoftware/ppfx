
#ifndef THINTARGETPIPCPIPREWEIGHTER_H
#define THINTARGETPIPCPIPREWEIGHTER_H

#include "IInteractionReweighting.h"
#include "ParameterTable.h"

#if defined __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-private-field"
#endif
namespace NeutrinoFluxReweight{
  
  /*! \class ThinTargetpipCpiponReweighter
   *  \brief Reweighter of thin target pion production
   */
  
  class ThinTargetpipCpipReweighter : public IInteractionReweighting{
  public:  
    ThinTargetpipCpipReweighter(int iuniv, const ParameterTable& cv_pars, const ParameterTable& univ_pars);
    virtual ~ThinTargetpipCpipReweighter();
    virtual bool canReweight(const InteractionData& aa);
    virtual double calculateWeight(const InteractionData& aa);    
    //double calculateDataScale(int inc_pdg, double inc_mom, int prod_pdg,double xf, double pt);

    double data_prod_xs;
    std::vector<float> vbin_data_pipC_pip;
    std::vector<float> vbin_data_pipC_pim;   
    std::vector<float> vbin_data_pipC_kp;
    std::vector<float> vbin_data_pipC_km;
    std::vector<float> vbin_data_pipC_p;
    std::vector<float> vbin_data_pipC_k0s;
    std::vector<float> vbin_data_pipC_lam;     
    std::vector<float> vbin_data_pipC_alam;
    //std::vector<float> vbin_datasys_pipC_pip;
    
  private:
    int iUniv;
    const ParameterTable& cvPars;
    const ParameterTable& univPars;
  //  float aux_par;
  };

  
}
#if defined __clang__
  #pragma clang diagnostic pop
#endif

#endif
