#******************************************************************************
#	Filename:	RunPlotter.sh
#	Authors:	Andres Florez 
#	Last update:	Aug. 23, 2012
#
#******************************************************************************

#--------------------------------------------------------------------------------------------------------------
# User can modify path to directory containing input ROOT files
inputPath=$PWD"/data"
#--------------------------------------------------------------------------------------------------------------

# Store command-line arguments and paths
configFileList=$PWD"/config/LogFiles.in"

root -l -b -q 'src/makeEfficiencyTable.cpp("'$inputPath'", "'$configFileList'")'
