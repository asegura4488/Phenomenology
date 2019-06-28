typeset -i n=1
typeset -i k=1
typeset -i events=1
INICIO=2000
FINAL=2000
for ((k=$INICIO;k<=$FINAL;k=k+1))
#for k in {1..3}
do 
#misedt="sed s/Pt/$k/ /home/ma.segura10/background/RunCard_Template_BG/pythia_cardinitial.dat"
#$misedt > /home/ma.segura10/background/RunCard_Template_BG/pythia_card.dat
mised="sed s/INIRUN=1/INIRUN=$k/ config_Integration.ini"
$mised > config_Integration$k.ini
mised2="sed s/ENDRUN=1/ENDRUN=$k/ config_Integration$k.ini"
$mised2 > config_Integration100.$k.ini
rm config_Integration$k.ini
mised3="sed s/config_Integration.ini/config_Integration100.$k.ini/ script_Integration.sh"
$mised3 > script_Integration$k.sh
chmod +x config_Integration100.$k.ini
chmod +x script_Integration$k.sh
time ./script_Integration$k.sh
rm  config_Integration100.$k.ini script_Integration$k.sh
done
