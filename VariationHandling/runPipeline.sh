#!/bin/sh

# Runs the entire pipeline if we want to re-run out analysis!

# Create all files for extended and reduced
#./runAll.sh

# Create all response matrices and efficiency/acceptance
#./runAllResponses.sh

# Create all combined files for reduced and extended 
# This step is needed only if you combine files BEFORE signal extraction!!!
# ./runCombineFiles.sh

# Create the mass fit signal templates from ttbar
./runCreateSignalTemplates.sh

# Create final Mass Fit for the 2 btag reduced region
./runMassFit.sh

# Create the signal extraction files (fiducial cross sections) for every year and variation
./runSignalExtraction.sh

# Combine Fiducial Measurements into 1 files per variation
# I need to manually create this step so this is tbd 
root -l CombineFiducialMeasurements.cxx

# Combine Acceptance and Efficiency 
python3 runEfficiency_Acceptance_Comb.py

# Unfolding Parton and Particle 
python3 Unfold_Combined.py

