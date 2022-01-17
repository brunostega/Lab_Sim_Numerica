cp PIMC/temperatures/input.pimc2.dat input.dat
make 
./qmc1d
mv probability.dat probability_T2.dat
mv probability_T2.dat PIMC/temperatures
mv potential.dat PIMC/
mv kinetic.dat PIMC/


cp PIMC/temperatures/input.pimc3.dat input.dat
make 
./qmc1d
mv probability.dat probability_T3.dat
mv probability_T3.dat PIMC/temperatures
mv potential.dat PIMC/
mv kinetic.dat PIMC/

cp PIMC/temperatures/input.pimc5.dat input.dat
make 
./qmc1d
mv probability.dat probability_T5.dat
mv probability_T5.dat PIMC/temperatures
mv potential.dat PIMC/
mv kinetic.dat PIMC/