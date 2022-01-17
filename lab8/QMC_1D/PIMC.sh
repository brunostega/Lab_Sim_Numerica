cp input.pimc input.dat
make 
./qmc1d
mv probability.dat PIMC/
mv potential.dat PIMC/
mv kinetic.dat PIMC/