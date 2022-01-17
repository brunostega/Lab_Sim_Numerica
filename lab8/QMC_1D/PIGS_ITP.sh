cp PIGS/ITP/input.pigsITP2.dat input.dat
make 
./qmc1d
mv probability.dat probability_ITP_2.dat
mv probability_ITP_2.dat PIGS/ITP
mv potential.dat PIMC/
mv kinetic.dat PIMC/


cp PIGS/ITP/input.pigsITP4.dat input.dat
make 
./qmc1d
mv probability.dat probability_ITP_4.dat
mv probability_ITP_4.dat PIGS/ITP
mv potential.dat PIMC/
mv kinetic.dat PIMC/

cp PIGS/ITP/input.pigsITP1_2.dat input.dat
make 
./qmc1d
mv probability.dat probability_ITP_1_2.dat
mv probability_ITP_1_2.dat PIGS/ITP
mv potential.dat PIMC/
mv kinetic.dat PIMC/