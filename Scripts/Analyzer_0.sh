# -------------------------------------------------
# -------     Universidad de los Andes      -------
# -------      Departamento de Física       -------
# -------     Manuel Alejandro Segura D     -------
# -------          Mayo - 2016              -------
# -------------------------------------------------

#--------------------------------------------------------------
#Directories 
#--------------------------------------------------------------
ANALYZERFOLDER="/home/ma.segura10/Analyzer"
PHENOANALYZERFOLDER="/home/ma.segura10/code_analysis/monotop_chirality"
PHENOANALYZEREXE="PhenoAnalyzer"
INROOTFILE="Events/run_01/output_delphes.root"
TEMPORALFOLDER="/home/ma.segura10/Analyzer/temporal"
OUTPUTFOLDER="/home/ma.segura10/Analyzer/output_folder"
#--------------------------------------------------------------
# Processes
#--------------------------------------------------------------- 
PROCESSFOLDER[2]="/Disco1/Pheno/monotop/monotop1jet_full"
PROCESSSSUBFOLDER[2]="monotop1jet_full"
RUNS[2]=11
TIMES[2]=1


PROCESSFOLDER[3]="/Disco1/Pheno/monotop/monotop1jet_left"
PROCESSSSUBFOLDER[3]="monotop1jet_left"
RUNS[3]=30
TIMES[3]=1


PROCESSFOLDER[4]="/Disco1/Pheno/monotop/monotop1jet_right"
PROCESSSSUBFOLDER[4]="monotop1jet_right"
RUNS[4]=30
TIMES[4]=1

PROCESSFOLDER[5]="/Disco2/Pheno/Backgrounds/DY+jets"
PROCESSSSUBFOLDER[5]="DY+jets"
RUNS[5]=3
TIMES[5]=1


PROCESSFOLDER[6]="/Disco1/Pheno/BackgroudSamples/singletop"
PROCESSSSUBFOLDER[6]="singletop"
RUNS[6]=50
TIMES[6]=6

PROCESSFOLDER[7]="/Disco2/Pheno/Backgrounds/W+jets"
PROCESSSSUBFOLDER[7]="W+jets"
RUNS[7]=50
TIMES[7]=14


PROCESSFOLDER[9]="/Disco2/Pheno/Backgrounds/W+jets-hadronic"
PROCESSSSUBFOLDER[9]="W+jets-hadronic"
RUNS[9]=50
TIMES[9]=6


#---------------------------------------------------------------

# Index process to run
#-----------------------
INDEX[1]=5
INDEX[2]=5
#----------------------

# Cut value
#--------------
VARIABLE="none"
#--------------

typeset -i start_cut=0
typeset -i end_cut=0
typeset -i delta=1
#-------------------------
# loop over the cut values
#------------------------
for ((i=$start_cut;i<=$end_cut;i=i+$delta))
do
parameters="sed s/mm/$i/ $PHENOANALYZERFOLDER/initial_parameters.in"
$parameters > $PHENOANALYZERFOLDER/config.in
cd $PHENOANALYZERFOLDER
make compile_ROOT_Delphes
cp *.h  $ANALYZERFOLDER
cp *.in $ANALYZERFOLDER 
cd -
	#-----------------------
	#loop over the processes
	#----------------------
        for j in `seq ${INDEX[1]} ${INDEX[2]}` ; 
        do
        	typeset -i start_times=1
        	typeset -i end_times=${TIMES[$j]} #Number of repetitions
                typeset -i start_runs=1
                typeset -i end_runs=${RUNS[$j]}   #Number of runs
		#loop over repetitions 
        	for k in `seq $start_times $end_times` ;
        	do
			#loop over runs
		    	for l in `seq $start_runs $end_runs` ;
                	do      
			$PHENOANALYZERFOLDER/$PHENOANALYZEREXE ${PROCESSFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$l/$INROOTFILE $TEMPORALFOLDER/${PROCESSSSUBFOLDER[$j]}_$l.root &   
   	   		echo  Running ${PROCESSSSUBFOLDER[$j]}_$l.root file
                	done
                	wait
        	start_runs=start_runs+${RUNS[$j]}
        	end_runs=end_runs+${RUNS[$j]}
                done
	hadd -f $OUTPUTFOLDER/${PROCESSSSUBFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$VARIABLE-$i.root $TEMPORALFOLDER/${PROCESSSSUBFOLDER[$j]}_*.root
	cd $TEMPORALFOLDER
        rm *
	cd - 
        echo FINISH ${PROCESSSSUBFOLDER[$j]} 
	done
done
