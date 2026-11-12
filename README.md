# PPFX - a Package to Predict the FluX

This is new repository for official PPFX development intended to support experiments over the long term.

It is maintained by the original authors (Leo Aliaga, Mike Kordosky) with help from others.

It was cloned on 1/13/23 from the Fermilab repository at

ssh://p-ppfx@cdcvs.fnal.gov/cvs/projects/ppfx 

as originally referenced in the Fermilab redmine link

https://cdcvs.fnal.gov/redmine/projects/ppfx/wiki/PPFX



This branch is specifically for NOvA and DUNE with NA61 60 GeV Reweighter. We have removed MIPP data covariances and introduced SVD decomposition for NA61 data.
Using the file scripts/inputa_na61.xml will switch on the NA61 reweighter. There is a new mode/variable in scripts/inputa_na61.xml, this will allow to switch on NA61 reweighter for the different experiments.
You need to carefullt specify which detector you are working with. If NOvA, use 1 else for DUNE use 2 in <Exp> variable.
The NA49 simulations are also for the Geant4 v_10.4, earlier v9.2 was used.


