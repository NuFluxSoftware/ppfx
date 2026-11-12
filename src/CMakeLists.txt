
cet_make_library(LIBRARY_NAME ppfx
                 SOURCE AbsorptionDPIPReweighter.cpp
                        AbsorptionDVOLReweighter.cpp
                        AbsorptionICReweighter.cpp
                        AttenuationMC.cpp
                        CentralValuesAndUncertainties.cpp
                        FillIMapHists.cpp
                        InteractionChainData.cpp
                        InteractionData.cpp
                        MakeReweight.cpp
                        MIPPNumiKaonYieldsReweighter.cpp
                        MIPPNumiMC.cpp
                        MIPPNumiPionYieldsReweighter.cpp
                        MIPPNumiYieldsBins.cpp
                        NucleonAbsorptionOutOfTargetReweighter.cpp
                        nu_g4numi.cpp
                        Numi2Pdg.cpp
                        NuWeight.cpp
                        OtherAbsorptionOutOfTargetReweighter.cpp
                        OtherReweighter.cpp
                        ParameterTable.cpp
                        ParticlesThroughVolumesData.cpp
                        ReweightDriver.cpp
                        TargetAttenuationReweighter.cpp
                        TargetData.cpp
                        ThinTargetBins.cpp
                        ThinTargetMC.cpp
                        ThinTargetMesonIncidentReweighter.cpp
                        ThinTargetnCPionReweighter.cpp
                        ThinTargetnucleonAReweighter.cpp
                        ThinTargetpCKaonReweighter.cpp
                        ThinTargetpCNucleonReweighter.cpp
                        ThinTargetpCPionReweighter.cpp
                        ThinTargetpipCpipBins.cpp
                        ThinTargetpipCpipMC.cpp
                        ThinTargetpipCpipReweighter.cpp
                 LIBRARIES PUBLIC
                        dk2nu::Tree
                        ROOT::EG
                        ROOT::Tree
                        Boost::filesystem
)

# executables
set( ppfx_execs
CreateIMapHists
CreateIMapPlots
doReweight_dk2nu_originalpip
make_basic_fluxhistograms
minos2nova
PrintMaterials
trial_ivol
)
 
foreach( ppfx_exec_name ${ppfx_execs} )
cet_make_exec( NAME ${ppfx_exec_name}
               SOURCE ${ppfx_exec_name}.C
               LIBRARIES PRIVATE ppfx::ppfx
               )
endforeach( ppfx_exec_name )

install_source()

