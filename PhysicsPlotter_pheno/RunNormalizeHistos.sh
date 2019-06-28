#******************************************************************************
#	Filename:	RunPlotter.sh
#	Authors:	Andres Florez 
#	Last update:	Aug. 23, 2012
#
#	This shell script is the driver for NormalizeHistos.cpp. This script makes an output
#	e.g.		sh RunPlotter.sh General Test
#******************************************************************************

#--------------------------------------------------------------------------------------------------------------
# User can modify path to directory containing input ROOT files
inputPath=$PWD"/data"
OutPath="/home/ma.segura10/PhysicsPlotter_pheno"
#--------------------------------------------------------------------------------------------------------------

# Store command-line arguments and paths
outputDir=$1						# e.g. Test
outputPath=$PWD"/"$outputDir				# e.g. ~/RootFiles/Test
configFileList=$PWD"/config/NormalizeHistos.in"
configFilePath=$PWD"/config"

tempFilePath=$PWD"/config/NormalizeHistos.in"	# Path to temporary file

# Create directory to store output files
if [ ! -d $outputPath ]; then
    mkdir $outputPath
fi

thefile=$PWD"/data"

# Open the file "regions.txt" and cycle through each region
for region in `cat $PWD"/config/directories.in"`; do
    
    regionPath=$outputPath"/"$region # e.g. ~/RootFiles/Test/sel_leptons
    
    # Create a subdirectory with the same name as the region
    if [ ! -d $regionPath ]; then
	mkdir $regionPath
    fi
    
    for conf in `cat $configFileList`; do
        
        cp $configFilePath"/"$conf $regionPath
        
	# Substitute all "REGION" words in each ROOT file's cutflow path with the region
	#sed -i s/"REGION"/$region/g $tempFilePath
	#-------------------------------------------------------------------------------------------------------
	# May need to create a backup file for MAC operating systems
        sed -i s:"FILE":$thefile:g $regionPath"/"$conf
        sed -i s:"PATH":$OutPath:g $regionPath"/"$conf
        sed -i s:"DIRECTORY":$region:g $regionPath"/"$conf
	#-------------------------------------------------------------------------------------------------------
    done 	
	# Run NormalizeHistos.cpp (l: load fast, b: batch mode, q: quit afterward)
	root -l -b -q 'src/NormalizeHistos.cpp("'$inputPath'", "'$configFileList'", "'$regionPath'", "'$region'")'
	
        mv *.root $regionPath
done	
