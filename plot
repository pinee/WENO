#!/bmein/bash
 
for i in `seq 1001 1695`;
do
name="datarho${i}"
gnuplot << plotting_many_files
set xlabel "x"
set ylabel "rho"
xy="./$name"
set terminal pngcairo
set output "$name.png"
plot xy using 1:2 with lines
plotting_many_files
done
