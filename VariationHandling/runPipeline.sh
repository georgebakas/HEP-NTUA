#!/bin/sh

# Runs the entire pipeline if we want to re-run out analysis!

# Create all files for extended and reduced
#./runAll.sh

# Create all response matrices and efficiency/acceptance
#./runAllResponses.sh

# Create all combined files for reduced and extended 
#./runCombineFiles.sh

# Create the mass fit signal templates from ttbar
./runCreateSignalTemplates.sh

# Create final Mass Fit for the 2 btag reduced region
./runMassFit.sh

# Create the signal extraction files (fiducial cross sections) 
./runSignalExtraction.sh

# Pending

# Unfolding 
#Parton

#Particle