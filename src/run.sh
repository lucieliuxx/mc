#!/bin/bash

if [ $# -ne 2 ]
then
	echo "run like this: src/run.sh RUN_ID simulate/madelung/converge/analysis/profile/movie"
	exit 1
fi 

run_id=$1
if [[ $2 == "bleach" ]]
then
	rm $run_id*.out
	rm $run_id*.trj
	rm $run_id*.config
	rm $run_id*.xyz
	rm $run_id*.sc
	rm $run_id*.uw
fi

if [[ $2 == "prepare" ]]
then
    #export OMP_NUM_THREADS=4
	gcc-9 -o main.o src/main.c src/utilities.c src/initialize.c src/fileio.c src/energycalc.c src/montecarlo.c src/cluster.c src/dipole.c -lm -ffast-math -fopenmp -O3
	./main.o $run_id.prep
fi

if [[ $2 == "simulate" ]]
then
    #export OMP_NUM_THREADS=4
	gcc-9 -o main.o src/main.c src/utilities.c src/initialize.c src/fileio.c src/energycalc.c src/montecarlo.c src/cluster.c src/dipole.c -lm -ffast-math -fopenmp -O3
	./main.o $run_id.in
fi


if [[ $2 == "madelung" ]]
then
	gcc-9 -o madelung.o src/madelung.c src/utilities.c src/initialize.c src/fileio.c src/energycalc.c src/montecarlo.c src/cluster.c src/dipole.c -lm -ffast-math
	./madelung.o
fi


if [[ $2 == "converge" ]]
then
	gnuplot -persist <<-EOFMarker
	    set key font ",25"
	    set term x11 font "Times, 25"
	    plot "converge.dat" index 0:5 using 1:4 with linespoints
	EOFMarker
fi


if [[ $2 == "analysis" ]]
then
	gcc-9 -o analysis.o src/analysis.c src/utilities.c src/initialize.c src/fileio.c src/energycalc.c src/montecarlo.c src/cluster.c src/dipole.c -lm -ffast-math
	./analysis.o $run_id.in 
fi

if [[ $2 == "plot" ]]
then
	python3 src/plot_data.py $run_id.in
fi

if [[ $2 == "profile" ]]
then
	gcc-9 -o main.o src/main.c src/utilities.c src/initialize.c src/fileio.c src/energycalc.c src/montecarlo.c src/cluster.c src/dipole.c -lm -ffast-math -fopenmp -g
	rm callgrind.out*
	valgrind --tool=Callgrind ./main.o $run_id.in
	qcachegrind callgrind.out*
fi

if [[ $2 == "movie" ]]
then
	read -p "which trajectory segment?: " ans
	
	if [[ "$ans" = "full" ]]
	then
		echo "Generating concatenated movie for all segments"
		xyz_file_name="$run_id""full.xyz"
		colorfile="$run_id""full.sc"
		rm $xyz_file_name
		for f in $run_id"seq"*.trj; do
            gcc-9 -o convert_trj_to_xyz.o src/convert_trj_to_xyz.c src/utilities.c src/initialize.c src/fileio.c src/energycalc.c src/montecarlo.c src/cluster.c -lm -ffast-math
            ./convert_trj_to_xyz.o $run_id.in $f $xyz_file_name 0 10000 1
		done
	else	
		echo "Generating movie for selected segment $ans"
		xyz_file_name="$run_id""seq$ans.xyz"
		colorfile="$run_id""seq$ans.sc"
		rm $xyz_file_name
        gcc-9 -o convert_trj_to_xyz.o src/convert_trj_to_xyz.c src/utilities.c src/initialize.c src/fileio.c src/energycalc.c src/montecarlo.c src/cluster.c -lm -ffast-math
        ./convert_trj_to_xyz.o $run_id.in $run_id"seq"$ans.trj $xyz_file_name 0 10000 10

	fi	

	# make rendering parameter file for vmd
	echo color Display Background white > $colorfile
	echo display depthcue off >> $colorfile
	echo display projection Orthographic >> $colorfile
	echo axes location off >> $colorfile
	echo topo readvarxyz $xyz_file_name >> $colorfile
	echo mol delete 0 >> $colorfile

	echo  >> $colorfile
	echo mol addrep 1 >> $colorfile
	echo mol modselect 1 1 name 1 >> $colorfile # name 1 = water
	echo mol modcolor 1 1 ColorID 23 >> $colorfile # ColorID 0 = blue2
	echo mol modmaterial 1 1 Transparent >> $colorfile
	echo mol modstyle 1 1 VDW 0.4 12.0 >> $colorfile

	echo  >> $colorfile
	echo mol addrep 1 >> $colorfile
	echo mol modselect 2 1 name 2 >> $colorfile # name 2 = oil
	echo mol modcolor 2 1 ColorID 1 >> $colorfile # ColorID 1 = red
	echo mol modmaterial 2 1 AOChalky >> $colorfile
	echo mol modstyle 2 1 VDW 0.4 12.0 >> $colorfile

	echo  >> $colorfile
	echo mol addrep 1 >> $colorfile
	echo mol modselect 3 1 name 3 >> $colorfile # name 3 = vac
	echo mol modcolor 3 1 ColorID 8 >> $colorfile # ColorID 8 = white
	echo mol modmaterial 3 1 Transparent >> $colorfile
	echo mol modstyle 3 1 VDW 0.4 12.0 >> $colorfile

	echo  >> $colorfile
	echo mol addrep 1 >> $colorfile
	echo mol modselect 4 1 name 4 >> $colorfile # name 4 = head
	echo mol modcolor 4 1 ColorID 15 >> $colorfile # ColorID 15 = ice blue
	echo mol modmaterial 4 1 AOChalky >> $colorfile
	echo mol modstyle 4 1 VDW 0.4 12.0 >> $colorfile

	echo  >> $colorfile
	echo mol addrep 1 >> $colorfile
	echo mol modselect 5 1 name 5 >> $colorfile # name 5 = tail
	echo mol modcolor 5 1 ColorID 32 >> $colorfile # ColorID 32 = orange3
	echo mol modmaterial 5 1 AOChalky >> $colorfile
	echo mol modstyle 5 1 VDW 0.4 12.0 >> $colorfile

	echo mol delrep 0 1 >> $colorfile

	vmd $xyz_file_name -e $colorfile
fi
