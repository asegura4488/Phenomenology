###### eMuPlotter running script
# Pass one of the following arguments to run the eMuPlotter over the sample 
# General (ttbar, Ztautau, Zprime_tautau_SSM*, Diboson), QCD_J*, WenuNp*, WmunuNp*

SAMPLE=$1
Directory=$2

for region in `cat config/directories.log`; do
    cd ${Directory}/${region}
	rm normalizedHistos_${SAMPLE}.root 
	cd ../../
done    



