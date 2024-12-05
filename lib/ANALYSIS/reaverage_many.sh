Nmols=108
starting_step=18000000

function select_profiles {
    tail -n +4 singleconfs_monomers_profile_mgel${1}.dat > tmpall.$1
    while (( `wc -l tmpall.$1 | awk '{print $1}'` > 0 ))
    do
	head -n $2 tmpall.$1 > tmp.$1
	tail -n +$(( $2+1 )) tmpall.$1 > tmpall2.$1 &&
	    mv tmpall2.$1 tmpall.$1
	(( `head -n 1 tmp.$1 | awk '{print $1}'` >= ${starting_step} )) && cat tmp.$1 >> single-prof_${1}.dat
    done
}

[[ -f average_hydrodynamic-r_mgels.dat ]]  &&  rm average_hydrodynamic-r_mgels.dat
[[ -f average_asphericity_mgels.dat ]]  &&  rm average_asphericity_mgels.dat
for i in `echo $Nmols | awk '{for(i=1; i<=$1; i++) printf "%d ", i;}'`
do
    echo $i
    awk -v start=${starting_step} '(NR==1) {print $0} (NR>1 && $1>=start) {print $0}' convex_hulls_${i}.dat > tmp
    mv tmp convex_hulls_${i}.dat
    head -n 1 convex_hulls_${i}.dat > .tmp2.dat
    tail -n +2 convex_hulls_${i}.dat > .tmp.dat
    Rh_avg=`avg .tmp.dat 7`
    asph1_avg=`avg .tmp.dat 9`
    cat .tmp.dat | sort -n >> .tmp2.dat
    mv .tmp2.dat convex_hulls_${i}.dat
    rm .tmp.dat
    echo $i $Rh_avg >> average_hydrodynamic-r_mgels.dat
    echo $i $asph1_avg >> average_asphericity_mgels.dat
done
tf=`tail -n 1 convex_hulls_1.dat | awk '{print $1}'`
awk -v tf=$tf '{mol[NR]=$1; rh[NR]=$2; rh_avg+=$2} END {printf "# final_t mgels: "; for(i=1; i<=NR; i++) printf " %d", mol[i]; printf " average_rh\n"; printf "%d", tf; for(i=1; i<=NR; i++) printf " %f", rh[i]; printf " %f\n", rh_avg/NR; }' average_hydrodynamic-r_mgels.dat > .tmp.dat
mv .tmp.dat average_hydrodynamic-r_mgels.dat
awk -v tf=$tf '{mol[NR]=$1; asph[NR]=$2; asph_avg+=$2} END {printf "# final_t mgels: "; for(i=1; i<=NR; i++) printf " %d", mol[i]; printf " average_asphericity\n"; printf "%d", tf; for(i=1; i<=NR; i++) printf " %f", asph[i]; printf " %f\n", asph_avg/NR; }' average_asphericity_mgels.dat > .tmp.dat
mv .tmp.dat average_asphericity_mgels.dat

Nconfs2el=`grep -A 1 TIMESTEP com_mgels.dat | grep [0-9] | awk -v start=${starting_step} 'BEGIN {n=0} ($1<start) {n+=1} END {print n}'`
if (( $Nconfs2el > 0 ))
then
    tail -n +$(( $Nconfs2el*($Nmols+9)+1 )) com_mgels.dat > tmp
    mv tmp com_mgels.dat
fi

awk -v start=${starting_step} '($1>=start) {print $0}' rg_mgels.dat > tmp
mv tmp rg_mgels.dat
head -n 1 average_rg_mgels.dat > tmp
mv tmp average_rg_mgels.dat
awk -v Nmols=$Nmols '{for(i=2; i<=NF; i++) {rg[i]+=$i; rg_avg+=$i}} END {printf "%d", NR ; for(i=2; i<=NF; i++) printf " %f", rg[i]/NR; printf " %f\n", rg_avg/NR/Nmols;}' rg_mgels.dat >> average_rg_mgels.dat

for i in `echo $Nmols | awk '{for(i=1; i<=$1; i++) printf "%d ", i;}'`
do
    Nlines=`head monomers_profile_mgel${i}.dat | awk '(NR==4) {print $2+1}'`
    select_profiles $i $Nlines &
done
wait

for i in `echo $Nmols | awk '{for(i=1; i<=$1; i++) printf "%d ", i;}'`
do
    echo ${i}
    head -n 3 singleconfs_monomers_profile_mgel${i}.dat > prova.dat
    cat single-prof_${i}.dat >> prova.dat
    average_radial_distributions.py -lmp_out prova.dat
    mv prova.dat monomers_profile_mgel${i}.dat
    mv single-prof_${i}.dat singleconfs_monomers_profile_mgel${i}.dat
done
rm prova_windows.dat tmp.* tmpall.* 
