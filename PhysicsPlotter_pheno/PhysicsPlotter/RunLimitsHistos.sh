#******************************************************************************
#	Filename:	RunPlotter.sh
#	Author:		Andres Florez 
#	Last update:	Aug. 23, 2012
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
configFileList=$PWD"/config/StackPlotter.in"	# Path to temporary file
configFilePath=$PWD"/config"

# Create directory to store output files
if [ ! -d $outputPath ]; then
	mkdir $outputPath
fi

# Open the file "directories.in" and cycle through each region
for region in `cat $inputPath"/config/directories_limits.in"`; do
	
	regionPath=$outputPath"/"$region	# e.g. ~/RootFiles/Test/sel_leptons
	
	# Create a subdirectory with the same name as the region
	if [ ! -d $regionPath ]; then
		mkdir $regionPath
	fi

        inputfile=$inputDir"\/"$region

        for conf in `cat $configFileList`; do

            cp $configFilePath"/"$conf $regionPath
            cp $configFilePath"/style_limits.in" $regionPath
            # Substitute "INPUTDIREC" words in all the config files for the input root file
            #-------------------------------------------------------------------------------------------------------
            # May need to create a backup file for MAC operating systems
            sed -i .bak s/"INPUTDIREC"/$inputfile/g $regionPath"/"$conf
            #------------------------------------------------------------------------------------------------------- 
        done 

	root -l -b -q 'src/LimitsRootFileMaker.cpp("'$inputPath'", "'$configFileList'", "'$regionPath'")'
        mv InputHistosForLimitsVBF.root $regionPath  	
	# Rename and move the output ROOT files and log file
	
done

