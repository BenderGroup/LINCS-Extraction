declare -a cell_lines=("--cell_line=A375" "--cell_line=A549" "--cell_line=ASC" "--cell_line=HA1E" "--cell_line=HCC515" "--cell_line=HEPG2" "--cell_line=HT29" "--cell_line=MCF7" "--cell_line=NEU" "--cell_line=NPC" "--cell_line=PC3" "--cell_line=SKB" "--cell_line=VCAP")
declare -a pert_times=("--pert_time=24" "--pert_time=6")
declare -a pert_doses=("--pert_dose=10" "--pert_dose=5")

for i in "${cell_lines[@]}"; do
	for j in "${pert_times[@]}"; do
		for k in "${pert_doses[@]}"; do
			print $i
			print $j
			print $k
		       	python Extract_And_TAS.py --ncores=20 "$i" "$j" "$k" 
			wait
		done
	done
done

