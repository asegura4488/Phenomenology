#******************************************************************************
#	Filename:	RunPlotter.sh
#	Author:		Andres Florez 
#	Last update:	Aug. 23, 2012
#
#	This shell script is the driver for StackPlotter.cpp. This script makes an output
#	directory, calls the file StackPlotter.cpp, and renames and moves the two
#	output ROOT files and log file to the output directory. The output directory is a
#	subdirectory of the input directory.
#******************************************************************************

#--------------------------------------------------------------------------------------------------------------
# User can modify path to directory containing input ROOT files
inputPath=$PWD
inputDataPath=$PWD"/data" 
#--------------------------------------------------------------------------------------------------------------

# Store command-line arguments and paths
inputDir=$1				        # e.g. General
outputDir=$2					# e.g. Test
outputPath=$inputPath"/"$outputDir		# e.g. ~/RootFiles/Test
configFileList=$PWD"/config/StackPlotter.in"	# List of processes to run over 
configFilePath=$PWD"/config"

saveHistos=$PWD"/SaveHistos/"$outputDir
# Create directory to store output files
if [ ! -d $outputPath ]; then
	mkdir $outputPath
fi

# Open the file "directories.in" and cycle through each region
for region in `cat $inputPath"/config/directories.in"`; do
	
	regionPath=$outputPath"/"$region	# e.g. ~/RootFiles/Test/sel_leptons
	
	# Create a subdirectory with the same name as the region
	if [ ! -d $regionPath ]; then
		mkdir $regionPath
	fi

        inputfile=$inputDir"\/"$region

        for conf in `cat $configFileList`; do

            cp $configFilePath"/"$conf $regionPath
            cp $configFilePath"/style_unit.in" $regionPath
            # Substitute "INPUTDIREC" words in all the config files for the input root file
            #-------------------------------------------------------------------------------------------------------
            # May need to create a backup file for MAC operating systems
            sed -i s/"INPUTDIREC"/$inputfile/g $regionPath"/"$conf
            #------------------------------------------------------------------------------------------------------- 
        done 

        root -l -b -q 'src/UnitPlotter.cpp("'$inputPath'", "'$configFileList'", "'$regionPath'", "'$saveHistos'")'
        mv UnitPlots.root $regionPath  	
	# Rename and move the output ROOT files and log file
	
done

