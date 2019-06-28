#!/bin/bash
# -------------------------------------------------
# -------     Universidad de los Andes      -------
# -------      Departamento de Física       -------
# -------        Joven Investigador         -------
# -------  Andrés Felipe García Albarracín  -------
# -------    Juan Carlos Sanabria Arenas    -------
# -------------------------------------------------
# This file executes parallel simulations with the programs: MadGraph 5.2 + Pythia 8.2 + Delphes 3.2
# Based on Diego Sanz's execution file: scriptMGParallelV2.sh

# Load the parameter file
source config_Integration.ini
## make the RunCards Folder in the EVENTSFOLDER
mkdir ${EVENTSFOLDER}/RunCards
## make the ParamCard Folder in the EVENTSFOLDER
mkdir ${EVENTSFOLDER}/ParamCard
## copy the param card supplied to the EVENTSFOLDER/ParamCard and name it param_card.dat
cp ${PARAMCARDFOLDER}/${PARAMCARDFILE} ${EVENTSFOLDER}/ParamCard/param_card.dat

## first sequence for each run, where the madgraph files and the run cards are created
sequ () {
	## copy the run card frame to the RunCards directory and append the seed (counter $i)
	cp ${RUNCARDFOLDER}/${RUNCARDFILE} ${EVENTSFOLDER}/RunCards/run_card_$i.dat
	## copy the MadGraph file to the RunCards directory as mgParallelFile_$i
	cp ${MADGRAPHFILEFOLDER}/${MADGRAPHFILE} ${EVENTSFOLDER}/RunCards/mgFile_$i.mg5
	## copy the parameter pythia file to the RunCards directory
	cp ${PYTHIAPARAMFOLDER}/${PYTHIAPARAM} ${EVENTSFOLDER}/RunCards/input_pythia_$i.cmnd
	## copy the delphes card to the RunCards directory *** Delphes card is the same for all runs
	cp ${DELPHESCARDFOLDER}/${DELPHESCARD} ${EVENTSFOLDER}/RunCards/${DELPHESCARD}
	## change all the instances of SEED to the counter $i on the file run_card_$i.dat
	sed -i "s/SEED/$i/g" ${EVENTSFOLDER}/RunCards/run_card_$i.dat
	## change all the instances of SEED to the counter $i on the file mgParallelFile_$i.mg5
	sed -i "s/SEED/$i/g" ${EVENTSFOLDER}/RunCards/mgFile_$i.mg5
	## change all the instances of SEED to the counter $i on the file input_pythia_$i.cmnd
	sed -i "s/SEED/$i/g" ${EVENTSFOLDER}/RunCards/input_pythia_$i.cmnd
	## change all the instances of RUNEVENTSNUM to $NUMEVENTSRUN on the file run_card_$i.dat
	sed -i "s/RUNEVENTSNUM/$NUMEVENTSRUN/g" ${EVENTSFOLDER}/RunCards/run_card_$i.dat
	## change all the instances of FOLDEREVENTS to $EVENTSFOLDER on the file mgParallelFile.mg5
	sed -i "s|FOLDEREVENTS|$EVENTSFOLDER|g" ${EVENTSFOLDER}/RunCards/mgFile_$i.mg5
	## change all the instances of NUMBERCORES to $CORESNUMBER on the file mgParallelFile.mg5
	sed -i "s|NUMBERCORES|$CORESNUMBER|g" ${EVENTSFOLDER}/RunCards/mgFile_$i.mg5
	## change all the instances of SUBFOLDERNAME to $NAMESUBFOLDER on the file mgParallelFile_$i.mg5
	sed -i "s|SUBFOLDERNAME|$NAMESUBFOLDER|g" ${EVENTSFOLDER}/RunCards/mgFile_$i.mg5
	## change all the instances of RESULTSFOLDER to the name of the folder where the results are located
	sed -i "s|RESULTSFOLDER|${EVENTSFOLDER}/${NAMESUBFOLDER}_$i/Events/run_01|g" ${EVENTSFOLDER}/RunCards/input_pythia_$i.cmnd
	## change all the instances of RUNEVENTSNUM to $NUMEVENTSRUN on the file parameter pythia file 
	sed -i "s/RUNEVENTSNUM/$NUMEVENTSRUN/g" ${EVENTSFOLDER}/RunCards/input_pythia_$i.cmnd
}

## second sequence for each run, where the madgraph is called for each of the madgraph files (mgParallelFile_i.mg5). Pythia8 and Delphes are also executed
sequ2 () {
	source config_Integration.ini
        ## run madgraph with the corresponding madgraph file .mg5. all the messages are thrown to /dev/null
	## Madgraph execution
	$1/bin/mg5_aMC -f $2/RunCards/mgFile_$4.mg5 # &> /dev/null
        ## sleep for 1s. Important, for the wait order to work
        sleep 1s
        ## wait for previous subprocesses to finish
        wait
	# Uncompress .lhe.gz file
	gzip -d $2/$3_$4/Events/run_01/unweighted_events.lhe.gz

        ${PYTHIAPARAMFOLDER}/monotop_fixcolor.exe $2/$3_$4/Events/run_01/unweighted_events.lhe $2/$3_$4/Events/run_01/unweighted_events_1.lhe 
	## Pythia 8 execution
	${PYTHIA8FOLDER}/${PYTHIA8EXE} $2/RunCards/input_pythia_$4.cmnd $2/$3_$4/Events/run_01/output_pythia8.hep # &> /dev/null

	## Delphes execution
	${DELPHESFOLDER}/${DELPHESEXE} $2/RunCards/${DELPHESCARD} $2/$3_$4/Events/run_01/output_delphes.root $2/$3_$4/Events/run_01/output_pythia8.hep

	## ExRootAnalysis execution
	${EXROOTFOLDER}/${EXROOTEXE} $2/$3_$4/Events/run_01/output_pythia8.hep $2/$3_$4/Events/run_01/output_pythia8.root

	## The following lines convert the "Delphes" tree .root file created by Delphes to a "LHCO" tree .root file
	#${DELPHESFOLDER}/${DELPHESROOT2LCHO} $2/$3_$4/Events/run_01/output_delphes.root $2/$3_$4/Events/run_01/output_delphes.lhco
	#${EXROOTFOLDER}/${EXROOTLHCO} $2/$3_$4/Events/run_01/output_delphes.lhco $2/$3_$4/Events/run_01/output_delphes_lhco.root

	## Remove unnecessary files
#	rm $2/$3_$4/Events/run_01/output_delphes.lhco
#	rm $2/$3_$4/Events/run_01/output_pythia8.hep

}

## Uncompress .lhe.gz file
#gzip -d ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/unweighted_events.lhe.gz

## Pythia 8 execution
#${PYTHIA8FOLDER}/${PYTHIA8EXE} ${EVENTSFOLDER}/RunCards/${PYTHIAPARAM} ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/output_pythia8.hep # &> /dev/null

## Delphes execution
#${DELPHESFOLDER}/${DELPHESEXE} ${EVENTSFOLDER}/RunCards/${DELPHESCARD} ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/output_delphes.root ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/output_pythia8.hep

## ExRootAnalysis execution
#${EXROOTFOLDER}/${EXROOTEXE} ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/output_pythia8.hep ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/output_pythia8.root

## The following lines convert the "Delphes" tree .root file created by Delphes to a "LHCO" tree .root file
#${DELPHESFOLDER}/${DELPHESROOT2LCHO} ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/output_delphes.root ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/output_delphes.lhco
#${EXROOTFOLDER}/${EXROOTLHCO} ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/output_delphes.lhco ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/output_delphes_lhco.root

## Remove unnecessary files
#rm ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/output_delphes.lhco
#rm ${EVENTSFOLDER}/${NAMESUBFOLDER}_${SEED}/Events/run_01/output_pythia8.hep

export -f sequ
export -f sequ2
## start PARAMETERS variable
PARAMETERS=""
## loop to execute sequence "sequ" for all the values from $INIRUN to $ENDRUN
for i in `seq ${INIRUN} ${ENDRUN}`; do  # {21,28}; do ## `seq ${INIRUN} ${ENDRUN}`; do
        ## execute sequ
        sequ
        ## concatenate the variable PARAMETERS with the current value of $i
        PARAMETERS="$PARAMETERS ${i}"
done

## execute gnuparallel. Use %% as the replacement string instead of {}.
parallel -0 -I %% --gnu "sequ2 ${MADGRAPHFOLDER} ${EVENTSFOLDER} ${NAMESUBFOLDER} %%" ::: $PARAMETERS

