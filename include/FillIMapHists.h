#ifndef FILLIMAPHISTS_H
#define FILLIMAPHISTS_H

#include "CommonIMapIncludes.h"
#include "TChain.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include <string>

#include "MIPPNumiKaonYieldsReweighter.h"
#include "MIPPNumiPionYieldsReweighter.h"
#include "ThinTargetpCPionReweighter.h"
#include "ThinTargetpCKaonReweighter.h"
#include "ThinTargetnCPionReweighter.h"
#include "ThinTargetpCNucleonReweighter.h"
#include "ThinTargetMesonIncidentReweighter.h"
#include "ThinTargetnucleonAReweighter.h"
#include "ThinTargetpipCpipReweighter.h"   //BHUMIKA


struct FillIMapHistsOpts {
  float elow, ehigh;
  Int_t nuid;
  bool cut_thintarget, cut_mipp, cut_na61;                       //BHUMIKA, introduced na61 cut
};

struct FillIMapHistsReweighters {

  NeutrinoFluxReweight::MIPPNumiPionYieldsReweighter* NumiPions;     // restored for MINERvA/MIPP comparison, see scripts/README.md
  NeutrinoFluxReweight::MIPPNumiKaonYieldsReweighter* NumiKaons;

  NeutrinoFluxReweight::ThinTargetpCPionReweighter*         ThinTargetpCPion;
  NeutrinoFluxReweight::ThinTargetpCKaonReweighter*         ThinTargetpCKaon;
  NeutrinoFluxReweight::ThinTargetnCPionReweighter*         ThinTargetnCPion;
  NeutrinoFluxReweight::ThinTargetpCNucleonReweighter*      ThinTargetpCNucleon;
  NeutrinoFluxReweight::ThinTargetMesonIncidentReweighter*  ThinTargetMesonIncident;
  NeutrinoFluxReweight::ThinTargetpipCpipReweighter*        ThinTargetpipCpip;      //Bhumika

  NeutrinoFluxReweight::ThinTargetnucleonAReweighter*       ThinTargetnucleonA;

};


/*! \fn FillIMapHists
 * \brief Fills the interaction map hists, looping over dk2nu input
 * Loop over the events in the Dk2Nu chain passed in and  call FillOneEntry()
 * returns the total weight = total number of neutrino interactions which 
 * were used to fill the plots (non integer since g4numi is a weighted MC)
 */
double FillIMapHists(TChain* tdk2nu, TChain* dkmeta, HistList* hists, const FillIMapHistsOpts* opts);

/*! \fn FillOneEntry
 * \brief Fill hists with the data inside of the dk2nu and dkmeta ntuples
 * This is a reimplementation of the routine originally written by 
 * Maria (Majo) Bustamante and Leo Aliaga, and improved by others.
 *
 * this function returns the weight used in filling the histograms,
 * which is zero if the event was skipped, but otherwise corresponds
 * to the product of the decay and importance weights.
 */
double FillOneEntry(bsim::Dk2Nu* dk2nu, bsim::DkMeta* dkmeta, HistList* hists, const FillIMapHistsOpts* opts, FillIMapHistsReweighters* reweighters);


int FindIndexFromVolume(const std::string& volname);

int FindIndexFromParticleName(const std::string& wanted);

#endif // FILLIMAPHISTS_H
