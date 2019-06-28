# -------------------------------------------------
# -------     Universidad de los Andes      -------
# -------      Departamento de Física       -------
# -------     Manuel Alejandro Segura D     -------
# -------          Mayo - 2016              -------
# -------------------------------------------------

#-------------------------------------------------------------
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

OCESSFOLDER[1]="/Disco2/Pheno/monotop_left"
PROCESSSSUBFOLDER[1]="monotop_left"
RUNS[1]=20
TIMES[1]=1
PROCESSFOLDER[2]="/Disco2/Pheno/monotop_right"
PROCESSSSUBFOLDER[2]="monotop_right"
RUNS[2]=20
TIMES[2]=1
PROCESSFOLDER[3]="/Disco2/Pheno/monotop_full"
PROCESSSSUBFOLDER[3]="monotop_full"
RUNS[3]=20
TIMES[3]=1
PROCESSFOLDER[4]="/Disco2/Pheno/monotop1jet_left"
PROCESSSSUBFOLDER[4]="monotop1jet_left"
RUNS[4]=30
TIMES[4]=1
PROCESSFOLDER[5]="/Disco2/Pheno/monotop1jet_right"
PROCESSSSUBFOLDER[5]="monotop1jet_right"
RUNS[5]=30
TIMES[5]=1
PROCESSFOLDER[6]="/Disco2/Pheno/monotop1jet_full"
PROCESSSSUBFOLDER[6]="monotop1jet_full"
RUNS[6]=30
TIMES[6]=1

PROCESSFOLDER[7]="/Disco2/Pheno/Backgrounds/DY+jets"
PROCESSSSUBFOLDER[7]="DY+jets"
RUNS[7]=50
TIMES[7]=10
PROCESSFOLDER[8]="/Disco2/Pheno/Backgrounds/W+jets"
PROCESSSSUBFOLDER[8]="W+jets"
RUNS[8]=50
TIMES[8]=14
PROCESSFOLDER[9]="/Disco1/Pheno/BackgroudSamples/singletop"
PROCESSSSUBFOLDER[9]="singletop"
RUNS[9]=50
TIMES[9]=1
#---------------------------------------------------------------

# Index process to run
#-----------------------
INDEX[1]=1
INDEX[2]=1
#----------------------

# Cut value
#--------------
VARIABLE="none"
#--------------

typeset -i start_cut=3
typeset -i end_cut=4
typeset -i delta=1
#-------------------------
# loop over the cut values
#------------------------
for ((i=$start_cut;i<=$end_cut;i=i+$delta))
do
parameters="sed s/mm/$i/ $PHENOANALYZERFOLDER/initial_parameters.in"
$parameters > $PHENOANALYZERFOLDER/config_$i.in
	typeset -i start_decimal=0
        typeset -i end_decimal=9     ##set decimal points
	typeset -i delta_decimal=1
	#loop over decimal points
	for ((d=$start_decimal;d<=$end_decimal;d=d+$delta_decimal))
	do
	decimal="sed s/nn/$d/ $PHENOANALYZERFOLDER/config_$i.in"
        $decimal > $PHENOANALYZERFOLDER/config.in
	cd $PHENOANALYZERFOLDER
	#make compile_ROOT_Delphes
	cp *.h  $ANALYZERFOLDER
	cp config.in $ANALYZERFOLDER 
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
##				$PHENOANALYZERFOLDER/$PHENOANALYZEREXE ${PROCESSFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$l/$INROOTFILE $TEMPORALFOLDER/${PROCESSSSUBFOLDER[$j]}_$l.root &   
   	   			echo  Running ${PROCESSSSUBFOLDER[$j]}_$l.root file
                		done
                		wait
        		start_runs=start_runs+${RUNS[$j]}
        		end_runs=end_runs+${RUNS[$j]}
                	done
		hadd -f $OUTPUTFOLDER/${PROCESSSSUBFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$VARIABLE-$i.$d.root $TEMPORALFOLDER/${PROCESSSSUBFOLDER[$j]}_*.root
		cd $TEMPORALFOLDER
       		rm *
		cd - 
        	echo FINISH ${PROCESSSSUBFOLDER[$j]} 
		done
	done
	rm $PHENOANALYZERFOLDER/config_$i.in
done
