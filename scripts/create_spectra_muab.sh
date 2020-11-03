ntrajs=1028
touch spectra_muab.in
> spectra_muab.in
for ((i=1;i<=$ntrajs;++i));do
	directory=pulse_pump/traj_$i/restart_0
	amberout=$directory/nasqm_pulse_pump_t${i}_r0.out
	s1_energy=$(grep -A2 "Frequencies (eV)" $amberout |
		awk '{if (NR == 3){print $2}}')
	awk -v s1=$s1_energy '
	BEGIN {}
	{ if ($1 == 1) {printf "%8.4f%10.5f", $3+s1, $7}}
	END {printf "\n"} ' $directory/muab.out >> spectra_muab.in
done

