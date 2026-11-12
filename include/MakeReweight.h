#ifndef MAKEREWEIGHT_H
#define MAKEREWEIGHT_H

#include "ReweightDriver.h"
#include "InteractionData.h"
#include "TargetData.h"
#include "ParticlesThroughVolumesData.h"
#include "InteractionChainData.h"
#include "ParameterTable.h"
#include "CentralValuesAndUncertainties.h"
#include "MIPPNumiYieldsBins.h"
#include "ThinTargetBins.h"

//NEW ADDED BY BHUMIKA
#include "ThinTargetpipCpipBins.h"
//#include "ThinTargetpipCpimBins.h"
//#include "ThinTargetpipCkpBins.h"
//#include "ThinTargetpipCkmpBins.h"
//#include "ThinTargetpipCpBins.h"
//#include "ThinTargetpipCk0sBins.h"
//#include "ThinTargetpipClamBins.h"
//#include "ThinTargetpipCalamBins.h"

#include "MIPPNumiMC.h"
#include "ThinTargetpipCpipMC.h"


#include "nu_g4numi.h"
#include "Numi2Pdg.h"

#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/dk2nu.h"

#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

namespace NeutrinoFluxReweight{
    
  /*! \class MakeReweight
   *  \brief A class to make the reweight event by event.
   */
  class MakeReweight{
 
  private:   
    MakeReweight();
    ~MakeReweight();
    
  public:
    void SetOptions(std::string fileIn);    
    static MakeReweight* getInstance();
    static void resetInstance();

    int getExp() const;     // New function to return Exp
    
    /*! calculate the weights for this event using the old flux ntuple format
     * (filled in nu_g4numi object),  the tgtcfg (the target configuration, example: le010z) 
     * and the horncfg (horn configuration, example: 185i)
     */
    void calculateWeights(nu_g4numi* nu, const char* tgtcfg, const char* horncfg);

    //! create an interaction chain from the new dk2nu(dkmeta) format
    void calculateWeights(bsim::Dk2Nu* nu, bsim::DkMeta* meta);

    //! total weights
    std::vector<double> GetTotalWeights();

    //! number of universes used in this run
    int GetNumberOfUniversesUsed();

    //! get the vector of the weights for a given reweighter
    std::vector<double> GetWeights(std::string nameReweighter);

    //! get the cv weights
    double GetCVWeight();
    
    //! vector of Reweighter Drivers, one per universe
    std::vector<ReweightDriver*> vec_rws;
    
    //! Reweighter Drivers for the central value
    ReweightDriver* cv_rw;

  private:    
    /*!
     * Inititalize the job and configurethe ReweighterDrivier
     */
    void Initialize();
    /*! Parse seeting options */
    void ParseOptions();
    /*! Parse seeting options */
    void doTheJob(InteractionChainData* icd);
    
    std::vector<ParameterTable> cvPars,univPars;
    int Nuniverses;
    int Exp;
    std::string mippCorrOption;
    std::string fileOptions;    
    std::vector<double> vec_wgts;
    std::map<std::string,std::vector<double> > map_rew_wgts;
    double cv_wgt;

    static MakeReweight* instance;
    
  };

  
}
#endif
